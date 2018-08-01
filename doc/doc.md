
<h1 style="text-align:center"><font face="Times Roman" size=8>XMD ver-0.1</font></h1> 
***
<h4 style="text-align:right"><font face="Times Roman">Author &nbsp;XShinHe<br><1500011805@pku.edu.cn></font> </h2>
<h2 style="text-align:left"><font face="Times Roman"> It is a molecular dynamics procedure of MD/PIMD, PILD simulation. </font></h2>
Latest update on 2018.7
***

## Argument files
 * put.rc  
 The file contains arguments handling molecule dynamics simulation, all should in Atom Unit

| name  | implement                       |  notes                     |
|-------|---------------------------------|----------------------------|
|temp	| temperature                     | \[a.u.\]                   |
|sfreq	| system characteristic frequency | related to spectrum data   |
|coeff	| coeffienct of friction/collide  | (OPT)                      |
|gammaAD| adiabatic coefficienet of PILD  | (OPT)                      |
|nstep	| steps of a trajectory           | (OPT)                      |
|nsmp	| intervals for sampling          | avoid correlation          |
|methd	| transform method                | primitive, stag, NM        |
|scheme	| thermostat scheme               | middle, end, side          |
|thermo	| thermostat method               | Langevin, Andersen, NHC    |
|mirror	| virtual/real dynamics           | (default 0)                |
|opt	| optimization for arguments      | (default 0)                |
|bead	| the bead of PIMD/PILD           | bead=1; bead=2N, PIMD/PILD |
|dtime	| the time interval for each step | (OPT)                      |

* cnfg.rc  
This is the initial input information of mass and positions of given atoms to build the molecule object.  
An example file is as following, N is the number of atoms, and M is the dimension of space (=1,2,3)
```
	!system-name    N       M
	A1     mass1    A1x1    ...     A1xM
	A2     mass2    A2x1    ...     A2xM
	...    ...      ...     ...     ...
	AN     massN    ANx1    ...     ANxM
```

* info.now  
This is the file record some information of simulation or analysis of the current directory.  
(mainly used by `matplotlib`)  
If you carry out a simulation under a new directory, you'd better carefully edit information of info.now.  

## Compile and usage
* compile: by `make` command in each subdirectory, default using `gfortran` compiler, it will generate `pimd` executable file (for convenience, we alias it as `main`)
* usage:
`main [-arg1 arg2]`, note arg1 and arg2 should be pairs. The detailed choise see as following:  
```shell
		main -h  h    # show the help information 
		main -x  n    # run n-th trajectory                    need xmd.rst or not
		main -x  e    # run equilibrium                        need nothing
		main -x  r    # run a restart trajectory               need xmd.rst
		main -x  s    # run from a sampling file               need xmd.smp(default)
		main -x  l    # run from a filename list               need xmd.lst(default)
		main -y  i    # the way for smapling                   need xmd.rst
		main -m  i    # from new/old steps (i=0/i=1)           need xmd.rst
		main -e  i    # set mode (i=0/i=1)                     need -x e/n
		main -r  *    # alternative for xmd.rst                need -x r
		main -s  *    # alternative for xmd.smp                need -x s
		main -l  *    # alternative for xmd.lst                need -x f
		main -sn n    # sample from a n-th of xmd.smp/xmd.lst  need -x s
		main -p  *    # alternative for putrc_file
		main -c  *    # alternative for cnfgrc_file
		main -o  *    # namefor out_file
		main -a  *    # namefor ana_file
```

## Program and file structure
### Subprojects  
the XMD project contains several subproject, for different functions.  
> * MD1  
This is the one-dimensional molecular dynamics simulation program.  
> * MDX  
This is the multi-dimensional molecular dynamics simulation program. (whose dimension is decided from `cnfg.rc` file)  
> * MDI  
The vector version of MD1, removing OOP (Object Orient Programming) feature!  
> * MDV  
The vector version of MDX, removing OOP (Object Orient Programming) feature!  
> * MD3H2O  
The version do with the spectroscropy of H2O.  
> * sPIMD/sPILD  
The simple version of PIMD/PILD, with only 1-dimensional staging transform method and Langevin thermostat. All codes are written in a single file. (but they are fast than MD1)

### Files  
Here contains fortran files and python files. Fortran files mainly run molecular dynamics and output the result, meanwhile the python files mainly process the analyze the result, and plot figures.  
> * MyDef.f90  
provide global definitions, contains:  
>> 1. defines precision: __sp, dp, len0, len1, len2, len3__ 
>> 2. defines mathematics constant: __pi, twopi, pi2. sqrtpi, e__
>> 3. defines physics constant: __hb, kb, au\_m, au\_e, au\_hb, au\_ke, au\_c, au\_a0, au\_eh, au\_t, au\_temp, au\_kb, au\_beta__  
>> 4. defines IO setting: __out\_unit, ana\_unit, out\_file, ana\_file, putrc\_file, cnfgrc\_file, rst\_file, smp\-file, lst\_file__  
>> 5. defines sharing valuables: __my\_iostat, my\_exist, rand\_throw__  
>> 6 defines public routines; __init\_seed, random\_norm, send\_err__

> * AM\_script.f90  
One should note, the OOP may slow down the procedure, which may be not good.  
Note, the program doesn't use OOP properly. This can be improved better (e.g. inheit)
>> 1. defines elemental list  
>> 2. defines atom  
>> 3. defines mole (molecular object)

> * md\_info.f90  
> * md\_rout.f90  
> * main.f90  
> * dynamic\_stat.f90  
> * methd\_plus.f90  
> * thermo\_plus.f90  
> * restrc\_plus.f90  
> * elim\_plus.f90  

Python parts (All under the `./py` directory) 
> * setinfo.py  (need info.now for the current directory)
> * rvs.py
> * stat.py  
For plotting histogram of n-th column of xxx file. You can use as:  
`python3 $XMD/py/stat.py xxx n [N/Y]`, when using `N`, the figure doesn't show in terminal, while shows for using `Y` 
> * cct.py
For calculating correlation time function of n-th column of xxx file. You can use as:  
`python3 $XMD/py/cct.py xxx n [N/Y]`, when using `N`, the figure doesn't show in terminal, while shows for using `Y` 
> * cctpild.py
> * fft.py
> * IR.py

## Nomination principle


## Bugs
* PILD under testing  
* MD3H2O, the trajectory be confused, need elimination of transition and rotation properly.  

## future
* improve the OPP version, mainly `MD1` and `MDX`  
* complete the vector version, mainly `MDI` and `MDV`  
* think about the multi-PILD whether can be realized?  
* fix the bug of MD3H2O, and get the IR of H2O  
* add interface to VMD
