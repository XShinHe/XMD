!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 08
!--- Acknowledgement to Liu group, of PKU
!*


module md_rout
use MyDef
use AM_script
use md_info
use myFF
use thermo_plus
use pimd_plus
implicit none

contains

	subroutine md_ctrl(mobj)
	implicit none
		type(mole), intent(inout) :: mobj
		integer :: i,j
		
		!-- if needed, set Nose-Hoover chain settings
		if (md_ithermo.eq.2) then
		    call generate_NHC(md_nsum, 3, 5)
		end if
		
		!-- set transform settings (staging transform / normal mode transform)
		call set_trsfrm(md_bead)
		
		!-- special for equilibrium, where from x to update other values
		call calc_fx(mobj)
		select case (md_ipimd)
		    case(0)
		    case(1)
		        call calc_ks_stag(mobj)
		        call calc_fks_stag(mobj)
		    case(2)
		        call calc_ks_norm(mobj)
		        call calc_fks_norm(mobj)
		end select
		call calc_ana(mobj)
		
		!-- open relative record files
		if(md_y_mode .ge. 1) then
			open(unit=out_unit, file=out_file, status='replace')
			write(out_unit,*) 'nstep    ', 'atom   ','nbead    ','ks   ','x    ','p    ','fx  ', 'fks   '
		end if
		open(unit=ana_unit, file=ana_file, status='replace')
		write(ana_unit,*) 'nstep    ', 'atom   ','<x>    ','Pe   ','Ke    ', 'Te   ','ana1    '
		
		!-- select md-scheme and do simulation
		select case (md_ischeme)
		    !-- middle scheme
		    case (0)
		        call md_sampling(mobj)
		        do i = 1,md_nstep
		            !-- PXTXP
		            call md_propagate_p(2,mobj)
		            call md_propagate_x(2,mobj)
		            call md_propagate_T(1,mobj)
		            call md_propagate_x(2,mobj)
		            call md_propagate_p(2,mobj)
					!!
		            md_npass = md_npass + 1
		            if(mod(i,md_nsamp) .eq. 0) then
		                call md_sampling(mobj)
		            end if
		        end do
		    !-- end scheme
		    case (1)
		        call md_sampling(mobj)
		        do i = 1,md_nstep
		            !-- PXPT
		            call md_propagate_p(2,mobj)
		            call md_propagate_x(1,mobj)
		            call md_propagate_p(2,mobj)
		            call md_propagate_T(1,mobj)
		            !!
		            md_npass = md_npass + 1
		            if(mod(i,md_nsamp) .eq. 0) then
		                call md_sampling(mobj)
		            end if
		        end do
		    !-- side scheme
		    case (2)
		        call md_sampling(mobj)
		        do i = 1,md_nstep
		            !-- TPXPT
		            call md_propagate_T(2,mobj)
		            call md_propagate_p(2,mobj)
		            call md_propagate_x(1,mobj)
		            call md_propagate_p(2,mobj)
		            call md_propagate_T(2,mobj)
		            !!
		            md_npass = md_npass + 1
		            if(mod(i,md_nsamp) .eq. 0) then
		                call md_sampling(mobj)
		            end if
		        end do
		end select
		if(md_y_mode .ge. 1) then
			close(unit=out_unit)
		end if
		close(unit=ana_unit)
		
		!-- save finial configuration for restartion, total 8 terms
		open(unit=out_unit, file='xmd.rst', status='replace')
		do i=1,md_nsum
		    do j=1,md_bead
		        write (out_unit,*) md_npass, i, j, mobj%a(i)%ks(j), mobj%a(i)%x(j),&
		            mobj%a(i)%p(j),mobj%a(i)%fx(j),mobj%a(i)%fks(j)
		    end do
		end do
		close(unit=out_unit)
	end subroutine md_ctrl

	!-- sampling and wwritting the result
	subroutine md_sampling(mobj)
	implicit none
		type(mole), intent(in) :: mobj
		integer :: i,j
		!!
		do i=1,md_nsum
			!-- write out_file
			if(md_y_mode.ge. 1) then
				do j=1,md_bead
				    write (out_unit,*) md_npass, i, j, mobj%a(i)%ks(j), mobj%a(i)%x(j), &
				        mobj%a(i)%p(j),mobj%a(i)%fx(j),mobj%a(i)%fks(j)                   ! 8 term
				end do
		    end if
		    !-- write ana_file
		    if (md_x_mode .ne. 1 .or. md_npass > md_nstep/4) then
		        if(md_ipimd.ne.2) then
		        	!-- for staging-pimd, averging the position
		            write (ana_unit,*) md_npass, i, sum( mobj%a(i)%x )/real(md_bead), &
		                mobj%a(i)%pe, mobj%a(i)%ke, mobj%a(i)%te, 0  ! 8 term
		        else
		        	!-- for normal-mode-pimd(CMD)
		            write (ana_unit,*) md_npass, i, mobj%a(i)%ks(1), &
		                mobj%a(i)%pe, mobj%a(i)%ke, mobj%a(i)%te, 0 ! 8 term
		        end if
		    end if
		end do
	end subroutine md_sampling


	!-- the thermostat
	!------ md_ithermo = 0 'lang' : Langevin thermostat
	!------ md_ithermo = 1 'ads'  : Anderson thermostat
	!------ md_ithermo = 2 'nhc'  : Nose-Hoover chain
	subroutine md_propagate_T(n_dvd, mobj)
	implicit none
		integer, intent(in) :: n_dvd
		type(mole), intent(inout) :: mobj
		integer :: i,j
		real(dp) :: et
		!!
		select case(md_ithermo)
		    case (0) ! langevin
		        do i=1,md_nsum
		            do j=1,md_bead
		                call random_norm(et, rand_throw, 1.0_dp)
		                mobj%a(i)%p(j) = langc(n_dvd,1) * mobj%a(i)%p(j)  &
		                    + sqrt( mobj%a(i)%m * masscoeff(j) * md_temp ) * langc(n_dvd,2) * et
		            end do
		        end do
		    case (1) ! Anderson
		        do i=1,md_nsum
		            do j=1,md_bead
		                call random_number(rand_throw)
		                if ( rand_throw < 1.- langc(n_dvd,1) ) then
		                    call random_norm(et, rand_throw, 1.0_dp)
		                    mobj%a(i)%p(j) = et * sqrt( mobj%a(i)%m * masscoeff(j) * md_temp )
		                end if
		            end do
		        end do
		    case (2) ! NHC, see more in <thermo_plus.f90>
		        call NHChain(n_dvd, mobj)
		end select
	end subroutine md_propagate_T

	!-- propagation of position
	subroutine md_propagate_x(n_dvd, mobj)
	implicit none
		integer, intent(in) :: n_dvd
		type(mole), intent(inout) :: mobj
		integer :: i, j
		!!
		do i=1,md_nsum
		    do j=1,md_bead
		        mobj%a(i)%ks(j) = mobj%a(i)%ks(j) + mobj%a(i)%p(j) * md_dtime / ( real(n_dvd) * mobj%a(i)%m * masscoeff(j) )
		    end do
		end do
		!!
		select case (md_ipimd)
		    case(0) ! primary-pimd
		    case(1) ! staging-pimd
		        call calc_x_stag(mobj)
		        call calc_fx(mobj)
		        call calc_fks_stag(mobj)
		    case(2) ! normal mode pimd
		        call calc_x_norm(mobj)
		        call calc_fx(mobj)
		        call calc_fks_norm(mobj)
		end select
		call calc_ana(mobj)
	end subroutine md_propagate_x


	!-- propagation of momenta
	subroutine md_propagate_p(n_dvd, mobj)
	implicit none
		integer, intent(in) :: n_dvd
		type(mole), intent(inout) :: mobj
		integer :: i, j
		!!
		do i=1,md_nsum
		    do j=1,md_bead
		        mobj%a(i)%p(j) = mobj%a(i)%p(j) - mobj%a(i)%fks(j) * md_dtime / real(n_dvd)
		    end do
		end do
	end subroutine md_propagate_p


	!-- calculation of force (v.s. x)
	subroutine calc_fx(mobj)
	implicit none
		type(mole), intent(inout) :: mobj
		integer :: i, j
		do i=1,md_nsum
		    do j=1, md_bead
		        mobj%a(i)%fx(j) = Fpn2( 0.5_dp*mobj%a(i)%m,  mobj%a(i)%x(j) )
		    end do
		end do
	end subroutine calc_fx


	!-- calculation of analyzer (mainly analysis fo energy)
	subroutine calc_ana(mobj)
	implicit none
		type(mole), intent(inout) :: mobj 
		integer :: i,j
		do i=1,md_nsum
		    mobj%a(i)%pe = 0.0_dp
		    mobj%a(i)%ke = 0.0_dp
		    !-- original estimitor of kinetc energy, IF NEED, TO REMOVE NEXT LINE
		    ! mobj%a(i)%ke = 0.5_dp*md_bead*md_temp + 0.5_dp * mobj%a(i)%m * masscoeff(1)* md_bfreq2 * mobj%a(i)%ks(1)**2
		    do j=1,md_bead
		        mobj%a(i)%pe = mobj%a(i)%pe + Vpn2( 0.5_dp*mobj%a(i)%m, mobj%a(i)%x(j))
		        mobj%a(i)%ke = mobj%a(i)%ke + 0.5_dp * mobj%a(i)%x(j) * mobj%a(i)%fx(j)
		        !-- original estimitor of kinetc energy, , IF NEED, TO REMOVE NEXT LINE
		        ! mobj%a(i)%ke = mobj%a(i)%ke - 0.5_dp * mobj%a(i)%m * masscoeff(j)* md_bfreq2 * mobj%a(i)%ks(j)**2
		    end do
		    mobj%a(i)%pe = mobj%a(i)%pe / real(md_bead)
		    mobj%a(i)%te = mobj%a(i)%pe + mobj%a(i)%ke
		end do
	end subroutine calc_ana

end module md_rout



