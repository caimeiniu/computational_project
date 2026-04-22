Included is a commented lammps input script "in.lammps" as well as the file describing a W potential. These should be in the same directory
of the simulations. If one calls the command to run lammps as "lmp" then one would run a job as:

./lmp -var fn 0 -var nfi 10 -var seed 22031968 -echo screen < in.lammps

-var allows you to give variable values to the in.lammps script. In the above, the two relevant variables are:

fn: if zero, this starts a new CRA run in which the file relax.0.lammps is the perfect bcc lattice
nfi: is the number of Frenkel insertions, which is this case is 10

the supplied number 22031968 set the random seed, you should vary this everytime you want to run a new simulation

The output from lammps will be in this case stored in log.0.10.lammps

In addition, the configuration files

relax.0.lammps
relax.2.lammps
relax.4.lammps
....
relax.10.lammps

are outputted by the given script.

To run another 100 (say) FIs starting from relax.10.lammps, one would
run a job as:

./lmp -var fn 10 -var nfi 100 -var seed 04041973 -echo screen < in.lammps

Here the non-zero fn means that one would start from relax.10.lammps.

The remaining variable that is given to the lammps script is the random seed obtained from the system via the command line "echo $RANDOM"

Two further options for which in.lammps has to be modified directly concerns the system size and the frequency of configuration output.

When fn=0, the size of the bcc cell created is via the lammps script line:

region box block 0 nx 0 ny 0 nz

which would create a nx*ny*nz unit cell system.

The frequency of configuration output is via the lammps script line:

variable a equal ${i}%x

where every x'th configuration is outputed to relax.*.lammps

This functionality should make CRA simulations using lammps rather straight forward. The current script fixes the sigma_xx to be zero and
the in-plane strain to be zero.

Note, you will see that I do not worry about how close the displaced atom is to another atom (which is difficult to do entirely within a
lammps script). Most potentials can handle this extreme initial condition. If not, it will crash. 

