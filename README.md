# XMD
Molecular Dynamics Package of staging/NM PIMD, PILD, (CMD, RPMD, ELD, ECD, EHD etc.) simulation, mainly one dimensional.

# IO files
* \_put.in gives the parameter information, such as __temp sfreq	coeff nstep	nsmp methd scheme thermo mirror opt bead dtime__
* \_cofg.in give the initial atom masses and positions
* \_pot.in control the parameters of plot (by matplotlib)

# usage
* \*.exe can be pimd/pild/md3h2o/etc., as the xtag you set in make.inc when you compile it, we note in __exe__ here
* -x e, run for a configuration equilibrium
* -x r, run a trajectory from restart file __end.rst__
* -f n, run a trajectory from n-th sampling in the equilibrium file __e.trj__, note this argument is bound to cancel.

# function
* classMD, \[need add, but with setting bead=1 of pimd, we can get classMD\], except this one, else all to use the Path Ihtegral technologies, known as PIMD etc.
* spimd/spild, simple one that can compiles just as one file, only do 1-Dimension with Oscillator/Quardric potential with staging transform only, but need \_put.in also
* MD1/MD3, universal MD for 1 or 3 dimension simulation, with more settings of __methd__\[staging-pimd, normal-mode-pimd\], __thermo__\[langevin, Andersen, Nose-Hoover Chain\]
* MD1Atom/MD3Atom, for simple H/He atom calculation
* MD3H2O, for water calculation, use h2opes force field, now just for one molecule, and we now treat the IR information.

# bugs
* the trajectory of PILD is not correct.
* the elimination of transiton and rotation of water is not total correct.

# futures
* optimize the IO files
* fix the bugs
