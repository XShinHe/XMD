# <-    PIMD    verson 1.43     ->
# <-    2018-04-08              ->
# <-    manual  verson 1.43     ->
# <-    Copyright by Shin He, <1500011805@pku.edu.cn>   ->


# to compile the procedure, do:
>>> gfortran MyDef.f90 AM_script.f90 md_info.f90 myFF.f90 methd_plus.f90 thermo_plus.f90 h2opes.f90 md_rout.f90 main.f90 -o pimd
# or to debug by
>>> gfortran -fcheck=all -Wall MyDef.f90 AM_script.f90 liblapack.a md_info.f90 myFF.f90 methd_plus.f90 thermo_plus.f90 elim_plus.f90 h2opes.f90 md_rout.f90 main.f90 -o tmppimd

# features:

# Mydef.f90
#   1)  precision and length conrtol
#   2)  mathematics constant
#   3)  physics constant and Atomic Unit
#   4)  procedure control data
#   5)  shared variables
#   a)  mathematics random routine
#   b)  error control routine

# AM_script.f90
#   * this module is only for 1 dimension system!
#   1)  element table
#   2)  atom # here extend to high diemnsion by allocate 3N beads to an atom, each dimension has N beads!
#   3)  mole
#   a)  build mole routine

# md_info.f90
#   1)  ensemble control data
#   2)  dynamics control data
#   3)  adjust parameter
#   4)  step and time control data
#   5)  bead set control
#   6)  IO optimal optional choice
#   7)  system character frequence, (here is bound to abenden, move to AM_scrpit.f90)
#   8)  temporary values, (bound to abenden)
#   a)  init_md routine
#   b)  init_traj routine, (just init momenta, without calculate the position, force and energy)
#       or init from end.rst (md_mod=2) or from e.xpe (md_mod=3)

# thermo_plus.f90
#   1)  NHC related data
#   a)  generate_NHC, allocate related data
#   b)  NHC routines, for NHChain (using exponential function), NHChainP (using polynomial approxiamtion),
#       NHChain_partG, update G (virtual force) delayed!

# md_rout.f90
#   a)  md_ctrl, total arrangement
#   b)  md_smp, for sampling
#   c)  md_run_x, md_run_p, md_run_t
#   d)  staging transformation routines, calc_ks, calc_x, calc_fx, calc_fks, calc_pe
#       * calc_pe, here calculate of pe, ke, kev, te, tev
#       * (is bound to change this feature to methd_plus.f90)

# md_main.f90
#   *   this is the main routine that the hole program beign with.

# IO-files discription
#   1)  .xpe .xpf:  nstep, atomindex, nbead, ks, p, x, fks, fx
#   2)  .ena .ana:  nstep, atomindex, aver_x, pe, ke, kev , tote, totev

# update for version 1.43
#   1) add 3-dimension calculation
#   2) add methd_plus.f90, move staging, normal-mode modules to methd_plus, add parameter [md_methd] to md_info
#   



