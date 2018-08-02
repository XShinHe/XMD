FF = gfortran
EXE = pimd
default:
	$(FF) MyDef.f90 AM_script.f90 md_info.f90 myFF.f90 pimd_plus.f90 thermo_plus.f90 md_rout.f90 main.f90 -o $(EXE)
clean:
	rm *.mod
sweep all:
	rm *.png *.out *.ana *.mod xmd.rst xmd.smp xmd.lst $(EXE) 
