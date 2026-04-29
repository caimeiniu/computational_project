#!/usr/bin/env python3
"""Post-process one HMC run from `data/decks/hmc_AlMg.lammps`.

Inputs (paired by --stub):
  <stub>.log   thermo block including f_hmc[1]=n_attempts, f_hmc[2]=n_accepts
  <stub>.dump  per-frame `id type` (sorted by id)
  --gb-mask    boolean array of length N_atoms; True = GB site

Outputs:
  <stub>_xgb.json  thermo summary, swap stats, X_GB(t) series, block-bootstrap CI
  <stub>_xgb.png   4-panel diagnostic (T, PE, instantaneous accept, X_GB)

Burn-in: drops the first `--burnin-frac` of the production frames before
computing the stationary mean and bootstrap CI. Default 0.2.
"""
from __future__ import annotations
import argparse, json, re, sys
from pathlib import Path
import numpy as np


def parse_thermo_blocks(log_path: Path):
    """Return list of (header_cols, data_array). LAMMPS prints one header
    line per `run`; we keep all blocks and tag the last one (production).
    """
    text = log_path.read_text().splitlines()
    blocks = []
    i = 0
    while i < len(text):
        line = text[i].strip()
        if line.startswith("Step ") and len(line.split()) >= 4:
            cols = line.split()
            data = []
            j = i + 1
            while j < len(text):
                row = text[j].split()
                if not row:
                    break
                if len(row) != len(cols):
                    break
                try:
                    data.append([float(x) for x in row])
                except ValueError:
                    break
                j += 1
            if data:
                blocks.append((cols, np.array(data)))
            i = j
        else:
            i += 1
    return blocks


def parse_dump_types(dump_path: Path, n_atoms_expected: int):
    """Yield (timestep:int, types:np.ndarray[uint8] indexed by id-1) per frame.

    `dump_modify sort id` is set in the deck so `id` increases monotonically
    within each frame; we still build types[id-1] explicitly to be safe.
    """
    with dump_path.open() as f:
        while True:
            line = f.readline()
            if not line:
                return
            if line.strip() != "ITEM: TIMESTEP":
                continue
            ts = int(f.readline())
            assert f.readline().strip() == "ITEM: NUMBER OF ATOMS"
            n = int(f.readline())
            if n != n_atoms_expected:
                raise ValueError(f"frame at t={ts}: N={n} vs mask N={n_atoms_expected}")
            # box bounds
            assert f.readline().startswith("ITEM: BOX BOUNDS")
            for _ in range(3):
                f.readline()
            atoms_hdr = f.readline().split()
            assert atoms_hdr[:2] == ["ITEM:", "ATOMS"]
            cols = atoms_hdr[2:]
            id_idx = cols.index("id")
            type_idx = cols.index("type")
            # Bulk read of n atom lines via np.loadtxt — orders of magnitude
            # faster than per-line Python int conversion.
            block = np.loadtxt(f, dtype=np.int64, max_rows=n,
                               usecols=(id_idx, type_idx))
            types = np.empty(n, dtype=np.uint8)
            types[block[:, 0] - 1] = block[:, 1].astype(np.uint8)
            yield ts, types


def block_bootstrap_ci(x: np.ndarray, block: int, n_boot: int = 2000,
                        ci: float = 0.95, seed: int = 0):
    """Stationary-block bootstrap mean CI. `block` ~ 2 * integrated autocorrelation.
    With our short series (~50 frames) this is a coarse estimate; report it as such.
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    if n < 2 * block:
        # fallback: naive sample SE
        m = float(x.mean())
        se = float(x.std(ddof=1) / np.sqrt(n))
        return m, m - 1.96 * se, m + 1.96 * se
    n_blocks = n // block
    means = np.empty(n_boot)
    for b in range(n_boot):
        starts = rng.integers(0, n - block + 1, size=n_blocks)
        sample = np.concatenate([x[s:s + block] for s in starts])
        means[b] = sample.mean()
    lo = float(np.quantile(means, (1 - ci) / 2))
    hi = float(np.quantile(means, 1 - (1 - ci) / 2))
    return float(x.mean()), lo, hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stub", required=True,
                    help="Path stem; expects <stub>.log and <stub>.dump")
    ap.add_argument("--gb-mask", required=True, type=Path)
    ap.add_argument("--xc", type=float, required=True, help="bulk Mg fraction")
    ap.add_argument("--temp", type=float, required=True, help="target T (K)")
    ap.add_argument("--fd-pred", type=float, default=None,
                    help="X_GB^FD predicted at this (T, Xc); for gap report")
    ap.add_argument("--burnin-frac", type=float, default=0.2)
    ap.add_argument("--block", type=int, default=5,
                    help="block-bootstrap block size in frames")
    ap.add_argument("--out-prefix", type=Path, default=None)
    args = ap.parse_args()

    stub = Path(args.stub)
    # NB: stub names contain dots (e.g. "hmc_dry_T500_Xc0.005"), so append
    # extensions textually rather than using Path.with_suffix.
    log_path = Path(str(stub) + ".log")
    dump_path = Path(str(stub) + ".dump")
    out_prefix = args.out_prefix or Path(str(stub) + "_xgb")

    if not log_path.exists() or not dump_path.exists():
        sys.exit(f"missing {log_path} or {dump_path}")

    gb_mask = np.load(args.gb_mask).astype(bool)
    n_gb = int(gb_mask.sum())
    n_atoms = int(gb_mask.size)

    # ---- thermo: take the LAST block (production with f_hmc columns) ----
    blocks = parse_thermo_blocks(log_path)
    if not blocks:
        sys.exit("no thermo blocks parsed")
    cols, prod = blocks[-1]
    if "f_hmc[1]" not in cols:
        sys.exit(f"production block missing f_hmc cols: {cols}")
    cstep = cols.index("Step")
    ctemp = cols.index("Temp")
    cpe = cols.index("PotEng")
    catt = cols.index("f_hmc[1]")
    cacc = cols.index("f_hmc[2]")

    step = prod[:, cstep]
    temp = prod[:, ctemp]
    pe = prod[:, cpe]
    n_att = prod[:, catt]
    n_acc = prod[:, cacc]

    # instantaneous acceptance from differences
    d_att = np.diff(n_att, prepend=n_att[0])
    d_acc = np.diff(n_acc, prepend=n_acc[0])
    inst_accept = np.where(d_att > 0, d_acc / np.maximum(d_att, 1), np.nan)
    total_attempts = int(n_att[-1])
    total_accepts = int(n_acc[-1])
    accept_overall = total_accepts / total_attempts if total_attempts else 0.0

    # ---- X_GB(t) + fwd/rev swap decomposition by streaming through the dump ----
    # Forward (Mg INTO GB)  per inter-frame interval = #(GB sites flipping 1→2).
    # Reverse (Mg OUT of GB) per interval            = #(GB sites flipping 2→1).
    # Net (fwd − rev) equals the exact change in n_Mg_GB and is mass-conserved.
    # GB-GB swaps add 1 to BOTH fwd and rev (cancel in net) — bias the ratio
    # toward 1, so a deviation from 1 in the observed ratio is a robust
    # one-sided indicator of "still drifting" vs "near detailed balance".
    # bulk-bulk swaps and intra-bulk flips are invisible to these counters.
    timesteps = []
    x_gb = []
    x_total = []
    fwd_gb_per = []  # inter-frame: 1→2 flips inside GB
    rev_gb_per = []  # inter-frame: 2→1 flips inside GB
    prev_types = None
    for ts, types in parse_dump_types(dump_path, n_atoms):
        is_mg = (types == 2)
        x_gb.append(int(is_mg[gb_mask].sum()) / n_gb)
        x_total.append(int(is_mg.sum()) / n_atoms)
        timesteps.append(ts)
        if prev_types is not None:
            # restrict to GB sites once, then count flips
            gb_was = prev_types[gb_mask]
            gb_now = types[gb_mask]
            fwd_gb_per.append(int(((gb_was == 1) & (gb_now == 2)).sum()))
            rev_gb_per.append(int(((gb_was == 2) & (gb_now == 1)).sum()))
        prev_types = types.copy()
    timesteps = np.array(timesteps)
    x_gb = np.array(x_gb)
    x_total = np.array(x_total)
    fwd_gb_per = np.array(fwd_gb_per, dtype=np.int64)
    rev_gb_per = np.array(rev_gb_per, dtype=np.int64)

    # ---- burn-in + bootstrap ----
    n_frames = len(x_gb)
    n_burn = int(np.ceil(n_frames * args.burnin_frac))
    x_gb_stat = x_gb[n_burn:]
    mean, lo, hi = block_bootstrap_ci(x_gb_stat, block=args.block)

    enrichment = mean / args.xc if args.xc > 0 else None
    gap_vs_fd = (mean - args.fd_pred) if args.fd_pred is not None else None

    # ---- fwd/rev swap decomposition aggregate ----
    def _ratio(f, r):
        return float(f) / float(r) if r > 0 else None
    def _imbalance(f, r):
        s = f + r
        return float(f - r) / float(s) if s > 0 else None

    fwd_total, rev_total = int(fwd_gb_per.sum()), int(rev_gb_per.sum())
    # Per-frame indices i = 0..n_frames-2 correspond to transitions frame i→i+1.
    # Drop transitions originating before n_burn (use n_burn..end).
    if n_burn < len(fwd_gb_per):
        f_post = int(fwd_gb_per[n_burn:].sum())
        r_post = int(rev_gb_per[n_burn:].sum())
    else:
        f_post = r_post = 0

    swap_decomp = {
        "note": ("inter-frame type-flip counts in GB sites: forward = 1→2 "
                 "(Mg into GB), reverse = 2→1. Each accepted GB-bulk swap "
                 "contributes exactly +1 to one of {fwd, rev}; GB-GB swaps "
                 "add +1 to BOTH so they bias the ratio toward 1. Net "
                 "(fwd-rev) equals the exact change in n_Mg_GB."),
        "all_frames": {
            "fwd_gb_total": fwd_total,
            "rev_gb_total": rev_total,
            "net_delta_n_mg_gb": fwd_total - rev_total,
            "fwd_over_rev": _ratio(fwd_total, rev_total),
            "imbalance_signed": _imbalance(fwd_total, rev_total),
        },
        "post_burnin": {
            "n_burn_frames_dropped": n_burn,
            "fwd_gb": f_post,
            "rev_gb": r_post,
            "net_delta_n_mg_gb": f_post - r_post,
            "fwd_over_rev": _ratio(f_post, r_post),
            "imbalance_signed": _imbalance(f_post, r_post),
        },
    }

    summary = {
        "stub": str(stub),
        "T": args.temp,
        "X_c": args.xc,
        "n_atoms": n_atoms,
        "n_gb_atoms": n_gb,
        "gb_fraction": n_gb / n_atoms,
        "swap": {
            "total_attempts": total_attempts,
            "total_accepts": total_accepts,
            "accept_rate_overall": accept_overall,
        },
        "thermo_prod": {
            "T_mean": float(np.mean(temp)),
            "T_std": float(np.std(temp)),
            "PE_initial": float(pe[0]),
            "PE_final": float(pe[-1]),
            "PE_drift_eV": float(pe[-1] - pe[0]),
        },
        "x_gb": {
            "n_frames_total": n_frames,
            "n_frames_burnin_dropped": n_burn,
            "block_bootstrap_block": args.block,
            "mean": mean,
            "ci95_lo": lo,
            "ci95_hi": hi,
            "x_total_mean": float(x_total[n_burn:].mean()),
            "enrichment_vs_bulk": enrichment,
        },
        "fd_predicted": args.fd_pred,
        "gap_hmc_minus_fd": gap_vs_fd,
        "swap_decomposition": swap_decomp,
        "series": {
            "timestep": timesteps.tolist(),
            "x_gb": x_gb.tolist(),
            "x_total": x_total.tolist(),
            # length n_frames-1; entry i = transition timesteps[i] → timesteps[i+1]
            "fwd_gb_per_frame": fwd_gb_per.tolist(),
            "rev_gb_per_frame": rev_gb_per.tolist(),
        },
    }

    out_json = Path(str(out_prefix) + ".json")
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"wrote {out_json}")
    print(f"  X_GB^HMC = {mean:.4f}  CI95 = [{lo:.4f}, {hi:.4f}]"
          f"  (n_prod_frames = {len(x_gb_stat)})")
    if args.fd_pred is not None:
        print(f"  X_GB^FD  = {args.fd_pred:.4f}   gap = {gap_vs_fd:+.4f}")
    print(f"  swap accept = {accept_overall:.3%}  ({total_accepts}/{total_attempts})")
    pb = swap_decomp["post_burnin"]
    if pb["fwd_over_rev"] is not None:
        print(f"  fwd/rev (post-burnin) = {pb['fwd_over_rev']:.3f}  "
              f"(fwd={pb['fwd_gb']}, rev={pb['rev_gb']}, "
              f"net={pb['net_delta_n_mg_gb']:+d}, "
              f"imbalance={pb['imbalance_signed']:+.3f})")
    else:
        print("  fwd/rev: insufficient post-burnin transitions")

    # ---- diagnostic plot ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; skipping plot")
        return

    fig, axes = plt.subplots(5, 1, figsize=(8, 11), sharex=False)

    axes[0].plot(step, temp, lw=0.8)
    axes[0].axhline(args.temp, ls="--", color="k", lw=0.6)
    axes[0].set_ylabel("T (K)")
    axes[0].set_title(f"HMC diagnostics — T={args.temp:g} K, "
                      f"X_c={args.xc:g}  ({stub.name})")

    axes[1].plot(step, pe, lw=0.8, color="C1")
    axes[1].set_ylabel("PE (eV)")

    axes[2].plot(step, inst_accept, lw=0.6, color="C2", marker=".", ms=2)
    axes[2].axhline(accept_overall, ls="--", color="k", lw=0.6,
                    label=f"overall {accept_overall:.2%}")
    axes[2].set_ylabel("swap accept (per call)")
    axes[2].set_xlabel("step")
    axes[2].set_ylim(0, max(0.5, np.nanmax(inst_accept) * 1.1) if np.any(np.isfinite(inst_accept)) else 1)
    axes[2].legend(fontsize=8)

    axes[3].plot(timesteps, x_gb, marker="o", ms=3, lw=0.8,
                 label=f"X_GB^HMC = {mean:.4f} [{lo:.4f},{hi:.4f}]")
    if n_burn > 0 and n_burn < n_frames:
        axes[3].axvspan(timesteps[0], timesteps[n_burn - 1], alpha=0.15,
                        color="gray", label=f"burn-in ({n_burn} frames)")
    if args.fd_pred is not None:
        axes[3].axhline(args.fd_pred, ls="--", color="r",
                        label=f"X_GB^FD = {args.fd_pred:.4f}")
    axes[3].axhline(args.xc, ls=":", color="b",
                    label=f"X_c (bulk) = {args.xc:g}")
    axes[3].set_ylabel("X_GB")
    axes[3].legend(fontsize=8, loc="best")

    # fwd/rev per inter-frame interval: stem-style on a shared time axis
    if len(fwd_gb_per) > 0:
        # transitions are between consecutive frames; place dots at the
        # later timestep so the x-axis aligns with the X_GB(t) panel.
        t_trans = timesteps[1:]
        axes[4].plot(t_trans, fwd_gb_per, "o-", ms=3, lw=0.6, color="C0",
                     label=f"fwd (1→2 in GB), Σ={fwd_total}")
        axes[4].plot(t_trans, rev_gb_per, "s-", ms=3, lw=0.6, color="C3",
                     label=f"rev (2→1 in GB), Σ={rev_total}")
        if pb["fwd_over_rev"] is not None:
            axes[4].set_title(
                f"swap fwd/rev decomposition  "
                f"(post-burnin fwd/rev = {pb['fwd_over_rev']:.3f}, "
                f"imbalance = {pb['imbalance_signed']:+.3f})",
                fontsize=10)
        if n_burn > 0 and n_burn < n_frames:
            axes[4].axvspan(timesteps[0], timesteps[n_burn - 1], alpha=0.15,
                            color="gray")
        axes[4].set_ylabel("flips per interval")
        axes[4].legend(fontsize=8, loc="best")
    axes[4].set_xlabel("timestep")

    fig.tight_layout()
    out_png = Path(str(out_prefix) + ".png")
    fig.savefig(out_png, dpi=130)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
