include make.inc
FC = gfortran
EXE = $(xflag)

xbs_file = $(xbs)/MyDef.f90 $(xbs)/AM_script.f90 $(xbs)/md_info.f90 $(xbs)/methd_plus.f90 $(xbs)/thermo_plus.f90 $(xbs)/elim_plus.f90
xff_file = $(xff)/h2opes.f90 $(xff)/myFF.f90
xflag_file = ./src/$(xflag)/md_rout.f90 ./src/$(xflag)/main.f90

default:
	$(FC) -fcheck=all -Wall  $(xbs_file) $(xff_file) $(xflag_file) -L./$(xliblapcak) -llapack -o $(EXE)
clean:
	rm -rf *.o *.mod $(EXE)
