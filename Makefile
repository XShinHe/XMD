include make.inc
FC = gfortran
EXE = $(xtag)

xlib_file = $(xlib)/MyDef.f90 $(xlib)/AM_script.f90 $(xlib)/md_info.f90 $(xlib)/methd_plus.f90 $(xlib)/thermo_plus.f90 $(xlib)/elim_plus.f90
xff_file = $(xff)/h2opes.f90 $(xff)/myFF.f90
xtag_file = ./src/$(xtag)/md_rout.f90 ./src/$(xtag)/main.f90

default:
	$(FC) -fcheck=all -Wall  $(xlib_file) $(xff_file) $(xtag_file) -L./$(xlapcak) -llapack -o $(EXE)
clean:
	rm -rf *.o *.mod $(EXE)
