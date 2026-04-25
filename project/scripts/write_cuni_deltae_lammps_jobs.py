#!/usr/bin/env python3
"""Write LAMMPS minimization jobs for sampled Cu->Ni segregation energies."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


LAMMPS_TEMPLATE = """log {log_name}
units metal
dimension 3
boundary p p p
atom_style atomic

read_data "{data_file}"
mass 1 63.546
mass 2 58.693

pair_style eam/alloy
pair_coeff * * "{potential_file}" Cu Ni
neighbor 0.3 bin
neigh_modify every 1 delay 0 check yes

set atom {atom_id} type 2

thermo 50
thermo_style custom step pe press
min_style cg
minimize 1.0e-12 1.0e-12 10000 100000

variable final_pe equal pe
print "FINAL_PE_EV {role} {atom_id} ${{final_pe}}"
write_data {relaxed_name}
"""


def read_sites(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, default=Path("data/cuni_3d/cuni_3d_deltae_sites.csv"))
    parser.add_argument("--data", type=Path, default=Path("data/cuni_3d/cuni_3d_polycrystal.lammps"))
    parser.add_argument("--potential", type=Path, default=Path("data/examples/nc_swap_CuNi/Cu_Ni_Fischer_2018.eam.alloy"))
    parser.add_argument("--jobs-dir", type=Path, default=Path("data/cuni_3d/lammps/deltae_jobs"))
    parser.add_argument("--lammps-command", default="lmp")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_sites(args.sites)
    args.jobs_dir.mkdir(parents=True, exist_ok=True)
    data_file = os.path.relpath(args.data.resolve(), args.jobs_dir.resolve())
    potential_file = os.path.relpath(args.potential.resolve(), args.jobs_dir.resolve())

    commands = []
    for idx, row in enumerate(rows):
        role = row["role"]
        atom_id = int(row["atom_id"])
        prefix = "bulk_reference" if role == "bulk_reference" else f"gb_{idx:04d}_atom_{atom_id}"
        input_path = args.jobs_dir / f"in_{prefix}.lammps"
        log_name = f"log_{prefix}.lammps"
        relaxed_name = f"relaxed_{prefix}.lammps"
        input_path.write_text(
            LAMMPS_TEMPLATE.format(
                log_name=log_name,
                data_file=data_file,
                potential_file=potential_file,
                atom_id=atom_id,
                role=role,
                relaxed_name=relaxed_name,
            ),
            encoding="utf-8",
        )
        commands.append(f"{args.lammps_command} -in {input_path.name}")

    run_all = args.jobs_dir / "run_all_deltae.sh"
    run_all.write_text("#!/usr/bin/env bash\nset -euo pipefail\n\n" + "\n".join(commands) + "\n", encoding="utf-8")
    run_all.chmod(0o755)
    print(f"Wrote {len(rows)} LAMMPS inputs to {args.jobs_dir}")
    print(f"Run from that directory with: ./run_all_deltae.sh")


if __name__ == "__main__":
    main()
