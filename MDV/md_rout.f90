! --------------------------- mdbsc module -------------------------
module md_rout
use MyDef
use AM_script
use md_info
use thermo_plus
use methd_plus
implicit none

contains

	!-- md_ctrl   :  control each step of propagation
	!-- argument  :  mobj # a molecule object (see AM_script.f90 for defination)
	!-- algorithm :  
	!------------ 1) thermostat method: we provide [Langevin thermostat], [Andersen thermostat], [NHC thermostat], 
	!------------------- (see thermo_plus.f90 for more details)
	!------------ 2) transform  method: we provide [normal-mode transform], [staging transform], [primitive transform]
	!------------------- (see methd_plus.f90 for more details)
	!------------ 3)* difference method: we provide [vecolity-Verlet integration], and there might be more methods to try
	!-------------------- 1)) 4-th order spiltting
	!-------------------- 2)) Leapfrog method
	!-------------------- 3)) RK
	!-------------------- 4)) predictor and corrector, or Beemann splitting
	!-------------------- 5))* other numerical method for propagation
	!------------ 4)* restrict method: needed! (see restrc_plus.f90 for more details)
	!------------ 5)* eliminaion method: provide [classical picture], [SVD method], (see elim_plus.f90 for more details)
	!-------------------- this version (MD3) doesn't consider the elimination, it add properly in MD3H2O version
	subroutine md_ctrl(mobj)
	implicit none
		type(mole), intent(inout) :: mobj
		integer :: i,d,j

		!-- if therno method is NHC, set Nose-Hoover chain parameters
		!-- for other thermo method, we neddn't add any extra setttings
		if (md_thermo.eq.2) then
			!-- default set the chain length = 3
			!---------- set the RESPA algorithm factor = 5 (use for multi time scale)
		    call generate_NHC(mobj%nb, 3, 5)
		end if

		!-- for using transform method, such as staging tranform/normal-mode transformset transform, it relates a transform matrix
		!------  for md_methd = 1 : calling staging transform settings
		!------  for md_methd = 2 : calling normal-mode transform settings
		call set_trsfrm(md_bead)

		!-- for start for an equilibrium/(and some general), the initial generated coordinates is x (set from cnfg.rc/ set to zero), 
		!------ we should synchronize to ks and fx from x
		if(md_x_mode .eq. 1 .or. md_x_mode .eq. 0) then
		    do i=1,mobj%nb
		    	!-- for differnt dimensional direction, we seperate the transformation
		    	do d=1,md_dim
		        	call x2ks_stag(mobj%a(i)%ks(d::md_dim), mobj%a(i)%x(d::md_dim) )
		        end do
		    end do
		    call calc_realfx(mobj)
		end if

		!-- opening the IO reacording files
		!----- we use ana_unit to record analysis information, except the index, we require it contains 8 therms
		!----- we use out_unit to record trajectory information, but if for saving memory, we donn' output such information
		open(unit=ana_unit, file=ana_file, status='replace')
		write(ana_unit,*) 'nstep  ', 'atom  ','x    ','Ep   ','Ek    ', 'Ekv   ','Etot    '
		if(md_y_mode .eq. 1) then
			open(unit=out_unit, file=out_file, status='replace')
			write(out_unit,*) 'nstep    ', 'atom  ','nbead    ','ks   ','p    ','x    ','fks  ', 'fx   '
		end if

		!-- for choosing MD scheme for thermostat
		!------ 1) md_scheme = 1: middle schem
		!------ 2) md_scheme = 2: end schene
		!------ 3) md_scheme = 3: side scheme
		select case (md_scheme)
			!-- middle scheme
		    case (0)
		        call md_sampling(mobj)
		        do i = 1,md_nstep
		            call md_propagate_p(2,mobj)
		            call md_propagate_x(2,mobj)
		            call md_propagate_T(1,mobj)
		            call md_propagate_x(2,mobj)
		            call md_propagate_p(2,mobj)
		            md_npass = md_npass + 1
		            if(mod(i,md_nsmp) .eq. 0) then
		                call md_sampling(mobj)
		            end if
		        end do
		    !-- end scheme
		    case (1)
		        call md_sampling(mobj)
		        do i = 1,md_nstep
		            call md_propagate_p(2,mobj)
		            call md_propagate_x(1,mobj)
		            call md_propagate_p(2,mobj)
		            call md_propagate_T(1,mobj)
		            md_npass = md_npass + 1
		            if(mod(i,md_nsmp) .eq. 0) then
		                call md_sampling(mobj)
		            end if
		        end do
		    !-- side scheme
		    case (2)
		        call md_sampling(mobj)
		        do i = 1,md_nstep
		            call md_propagate_T(2,mobj)
		            call md_propagate_p(2,mobj)
		            call md_propagate_x(1,mobj)
		            call md_propagate_p(2,mobj)
		            call md_propagate_T(2,mobj)
		            md_npass = md_npass + 1
		            if(mod(i,md_nsmp) .eq. 0) then
		                call md_sampling(mobj)
		            end if
		        end do
		end select
		close(unit=ana_unit)
		if(md_y_mode .eq. 1) then
			close(unit=out_unit)
		end if
		
		!-- save a replica of restart file named 'xmd.rst'
		!------ * so here we define: [restart format: each bead each line]
		!------ * note reuse the out_unit for output
		open(unit=out_unit, file='xmd.rst', status='replace')
		do i=1,mobj%nb
		    do j=1,mdn_bead
		    	!-- we only output index-info (three terms) and tranformed trajectory information (five terms)
		        write (out_unit,*) md_npass, i, j, mobj%a(i)%ks(j), mobj%a(i)%p(j), &
		        				   mobj%a(i)%x(j), mobj%a(i)%fks(j), mobj%a(i)%fx(j)
		    end do
		end do
		close(unit=out_unit)
		
	end subroutine md_ctrl

	!-- md_sampling :  control each step of propagation
	!-- argument    :  mobj # a molecule object (see AM_script.f90 for defination)
	subroutine md_sampling(mobj)
	implicit none
		type(mole), intent(in) :: mobj
		integer :: i,d,j ! number of number, dimension, bead(in 1-diemsnion)
		do i=1,mobj%nb
			do j=1,md_bead
			    write (ana_unit,*) md_npass, i, j, mobj%a(i)%x(j*md_dim-2),  mobj%a(i)%x(j*md_dim-1),  mobj%a(i)%x(j*md_dim)
			end do
		end do
	end subroutine md_sampling

	!-- md_propagate_T :  propagation desribe the thermostat
	!-- argument : ndvd # splitting number of the original time increment
	!------------- mobj # a molecule object (see AM_script.f90 for defination)
	!-- algorithms:
	!------ 1) md_thermo = 0 : Langevin thermostat
	!------ 2) md_thermo = 1 : Andersen thermostat
	!------ 3) md_thermo = 2 : Nose-Hoover chain thermostat, (see thermo_plus.f90 for more details)
	subroutine md_propagate_T(n_dvd, mobj)
	implicit none
		integer, intent(in) :: n_dvd
		type(mole), intent(inout) :: mobj
		integer :: i,j
		real(dp) :: et ! local variable for white noise
		
		!-- choosing thermostat methods
		!-- note each dimension of each bead is indepenment
		select case(md_thermo)
			!-- Langevin themostat
		    case (0)
		        do i=1,mobj%nb
		            do j=1,mdn_bead
		                call random_norm(et, rand_throw, 1.0_dp)
		                mobj%a(i)%p(j) = langc(n_dvd,1) * mobj%a(i)%p(j)  &
		                    + sqrt( mobj%a(i)%m * masscoeff( int((j+md_dim-1)/md_dim) ) * md_temp ) * langc(n_dvd,2) * et
		            end do
		        end do
		    !-- Andersen thermostat
		    case (1)
		        do i=1,mobj%nb
		            do j=1,mdn_bead
		                call random_number(rand_throw)
		                if ( rand_throw < 1.- langc(n_dvd,1) ) then
		                    call random_norm(et, rand_throw, 1.0_dp)
		                    mobj%a(i)%p(j) = et * sqrt( mobj%a(i)%m * masscoeff( int((j+md_dim-1)/md_dim) ) * md_temp )
		                end if
		            end do
		        end do
		    !-- Nose-Hoover chain thermostat, originally used for only 1-dimension case, need debug for higher dimension
		    case (2)
		        call NHChain(n_dvd, mobj)
		end select
	end subroutine md_propagate_T

	!-- md_propagate_x :  propagation of the position (actually propagate on the tansformed position ks)
	!-- argument : ndvd # splitting number of the original time increment
	!------------- mobj # a molecule object (see AM_script.f90 for defination)
	subroutine md_propagate_x(n_dvd, mobj)
	implicit none
		integer, intent(in) :: n_dvd
		type(mole), intent(inout) :: mobj
		integer :: i, j
		!-- note each dimension of each bead is indepenment
		do i=1,mobj%nb
		    do j=1,mdn_bead
		        mobj%a(i)%ks(j) = mobj%a(i)%ks(j) + mobj%a(i)%p(j) * md_dtime / & 
		        				( real(n_dvd) * mobj%a(i)%m * masscoeff( int((j+md_dim-1)/md_dim) ) )
		    end do
		end do
		!-- synchronize the original and transformed position, and with the force information
		call md_update(mobj)
		
		!-- analysis term sampling in the position-space (we sampling here for default, note it is configuration space)
		call calc_ana(mobj)
		
	end subroutine md_propagate_x

	!-- md_propagate_p :  propagation of the momenta (actually the conrresponding tansformed momenta)
	!-- argument : ndvd # splitting number of the original time increment
	!------------- mobj # a molecule object (see AM_script.f90 for defination)
	subroutine md_propagate_p(n_dvd, mobj)
	implicit none
		integer, intent(in) :: n_dvd
		type(mole), intent(inout) :: mobj
		integer :: i, j
		
		!-- note each dimension of each bead is indepenment
		do i=1,mobj%nb
		    do j=1,mdn_bead
		        mobj%a(i)%p(j) = mobj%a(i)%p(j) - mobj%a(i)%fks(j) * md_dtime / real(n_dvd)
		    end do
		end do
	end subroutine md_propagate_p

	!-- md_update : synchronize the original position, transformed positon and force(fx)
	!-- argument : mobj # a molecule object (see AM_script.f90 for defination)
	!------ * force calcaulton notations
	!------------ 1) two-body's interaction
	!------------ 2) full parameter PES(potential energy surface)
	!------ * we pass all molecule object to force calculation, so that the transform can be seperated in each dimensions
	!------------ where we use integer d to indicate the dimension choise.
	subroutine md_update(mobj)
	implicit none
		type(mole), intent(inout) :: mobj
		integer :: i, d
		
		select case (md_methd)
		    case(1)
		        do i=1,mobj%nb
		        	do d=1,md_dim
		            	call ks2x_stag(mobj%a(i)%x(d::md_dim), mobj%a(i)%ks(d::md_dim) )
		            end do
		        end do
		        call calc_realfx(mobj)
		        do i=1,mobj%nb
		        	do d=1,md_dim
			            call fx2fks_stag(mobj%a(i)%fks(d::md_dim), mobj%a(i)%fx(d::md_dim), mobj%a(i)%ks(d::md_dim), mobj%a(i)%m )
		            end do
		        end do
		    case(2)
		        do i=1,mobj%nb
		        	do d=1,md_dim
		            	call ks2x_norm(mobj%a(i)%x(d::md_dim), mobj%a(i)%ks(d::md_dim) )
		            end do
		        end do
		        call calc_realfx(mobj)
		        do i=1,mobj%nb
		        	do d=1,md_dim
		            	call fx2fks_norm(mobj%a(i)%fks(d::md_dim), mobj%a(i)%fx(d::md_dim), mobj%a(i)%ks(d::md_dim), mobj%a(i)%m )
		            end do
		        end do
		end select
	end subroutine md_update

	!-- calc_ana : calculate relative quantities for analysis output, in configuration space
	!-- argument : mobj # a molecule object (see AM_script.f90 for defination)
	!------ * note it only work for dim = 3
	!------ terms:
	!------------ 1) pe : 
	!------------ 2) ke : 
	!------------ 3) kev: 
	!------------ 4) te : 
	!------------ 5) tev: 
	subroutine calc_ana(mobj)
	implicit none
		type(mole), intent(inout) :: mobj
		integer :: i,j,k
		real(dp) :: scalar
		do i=1,mobj%nb
		    mobj%a(i)%pe = 0.0_dp
		    mobj%a(i)%ke = 0.5_dp * md_dim * md_bead * md_temp
		    mobj%a(i)%kev = 0.0_dp
		end do
		do i=1,mobj%nb
			do k=i+1,mobj%nb
				do j=1,md_bead
				    call V3PN2_01(scalar, mobj%a(i)%x(md_dim*j-2:md_dim*j), mobj%a(k)%x(md_dim*j-2:md_dim*j) )
				    mobj%a(i)%pe = mobj%a(i)%pe + scalar
				    mobj%a(k)%pe = mobj%a(k)%pe + scalar
				    mobj%a(i)%kev = mobj%a(i)%kev + 0.5_dp * ( mobj%a(i)%x(md_dim*j-2) * mobj%a(i)%fx(md_dim*j-2) + &
				            mobj%a(i)%x(md_dim*j-1) * mobj%a(i)%fx(md_dim*j-1) + mobj%a(i)%x(md_dim*j) * mobj%a(i)%fx(md_dim*j) )
				    if(j.eq.1) cycle
				    mobj%a(i)%ke = mobj%a(i)%ke - 0.5_dp * mobj%a(i)%m * masscoeff(j)* md_bfreq2 * ( mobj%a(i)%ks(md_dim*j-2)**2 + &
				            mobj%a(i)%ks(md_dim*j-1)**2 + mobj%a(i)%ks(md_dim*j)**2 )
				end do
			end do
		end do
		do i=1,mobj%nb
			mobj%a(i)%pe = mobj%a(i)%pe / real(md_bead)
			mobj%a(i)%kev = mobj%a(i)%kev / real(md_bead)
			mobj%a(i)%te = mobj%a(i)%pe + mobj%a(i)%ke
			mobj%a(i)%tev = mobj%a(i)%pe + mobj%a(i)%kev
		end do
	end subroutine calc_ana

	!-- here we provide two force field for test
	!------ 1) two-body forcefield: use L-J potential for example
	!------ 2) H2O forcefield: (please see project [MD3H2O] for more details)
	
	!-- two-body forcefield, use L-J potential test, (see myFF.f90 module)
	subroutine calc_realfx(mobj)
	use myFF
	implicit none
		type(mole), intent(inout) :: mobj
		integer :: i, j, k
		real(dp), dimension(3) :: vec3
		do i=1,mobj%nb
		    mobj%a(i)%fx = 0.0_dp
		end do
		!!
		do i=1,mobj%nb
		    do k=i+1,mobj%nb
		        do j=1,md_bead
		            call F3PN2_01( vec3, mobj%a(i)%x(md_dim*j-2:md_dim*j), mobj%a(k)%x(md_dim*j-2:md_dim*j) )
		            mobj%a(i)%fx(md_dim*j-2:md_dim*j) = mobj%a(i)%fx(md_dim*j-2:md_dim*j) + vec3
		            mobj%a(k)%fx(md_dim*j-2:md_dim*j) = mobj%a(k)%fx(md_dim*j-2:md_dim*j) - vec3
		        end do
		    end do
		end do
	end subroutine calc_realfx
	
	!-- H2O PES forcefield, (see h2opes.f90 for more details)
	subroutine calc_realfxh2o(mobj)
	use h2opes
	implicit none
		type(mole), intent(inout) :: mobj
		integer :: i, j
		real(dp) :: v
		real(dp), dimension(3,3) :: dv
		real(dp), dimension(3,3,3,3) :: ddv
		real(dp), dimension(3,3) :: cart
		!-- here v, dv, ddv represent potential, gradient, hessian repectively
		!-- here cart(1,:) record cartesian coordinates of H1 
		!-- here cart(2,:) record cartesian coordinates of O
		!-- here cart(3,:) record cartesian coordinates of H2
		
		!-- note the index of atoms: 1:H1, 2:O, 3:H2,
		!----- if the atom number is not 3, this routine will stop the procedure
		if(mobj%nb.ne.3) stop
		do i=1,mobj%nb
		    mobj%a(i)%fx = 0.0_dp
		end do
		
		!-- diffenert index of bead is seperated well
		do j=1,md_bead
		    cart(1,:)=mobj%a(1)%x(md_dim*j-2:md_dim*j)
		    cart(2,:)=mobj%a(2)%x(md_dim*j-2:md_dim*j)
		    cart(3,:)=mobj%a(3)%x(md_dim*j-2:md_dim*j)
		    ! falg=1, meaning only calc v and dv from cartesian
		    call h2opot(v,dv,ddv,cart,1)
		    mobj%a(1)%fx(md_dim*j-2:md_dim*j) = dv(1,:)
		    mobj%a(2)%fx(md_dim*j-2:md_dim*j) = dv(2,:)
		    mobj%a(3)%fx(md_dim*j-2:md_dim*j) = dv(3,:)
		end do
	end subroutine calc_realfxh2o
	
end module md_rout

