# if you see this line, it say that this makefile can not work!

include make.inc
FC = gfortran
EXE = pimd3
srccodes = MyDef.f90 AM_script.f90 md_info.f90 myFF.f90 methd_plus.f90 thermo_plus.f90 elim_plus.f90 h2opes.f90 md_rout.f90 main.f90

default:
	$(FC) -fcheck=all -Wall  $(srccodes) -L./Lapack -llapack -lrefblas -o $(EXE)
clean:
	rm -rf *.o *.mod $(EXE)
