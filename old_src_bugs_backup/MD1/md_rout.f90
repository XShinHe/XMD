! --------------------------- mdbsc module -------------------------
module md_rout
use MyDef
use AM_script
use md_info
use myFF
use thermo_plus
use methd_plus
implicit none

contains

subroutine md_ctrl(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i,j
    !!
    if (md_thermo.eq.2) then
        call generate_NHC(mobj%nb, 3, 5)
    end if
    call set_trsfrm(md_bead)
    !!
    call calc_fx(mobj)
    select case (md_methd)
        case(0)
        case(1)
            call calc_ks_stag(mobj)
            call calc_fks_stag(mobj)
        case(2)
            call calc_ks_norm(mobj)
            call calc_fks_norm(mobj)
    end select
    call calc_pe(mobj)
    !!
    open(unit=out_unit, file=traj_outfile, status='replace')
    write(out_unit,*) 'nstep    ', 'atomidx ','nbead    ','ks   ','p    ','x    ','fks  ', 'fx   '
    open(unit=ana_unit, file=anal_file, status='replace')
    write(ana_unit,*) 'nstep    ', 'atomidx ','aver_x    ','pe   ','ke    ', 'kev   ','tote    ','is_eqb  '
    !!
    select case (md_scheme)
        case (0) ! middle scheme
            call md_smp(mobj)
            do i = 1,md_nstep
                !!
                call md_run_p(2,mobj)
                call md_run_x(2,mobj)
                call md_run_t(1,mobj)
                call md_run_x(2,mobj)
                call md_run_p(2,mobj)
                !!
                md_npass = md_npass + 1
                if(mod(i,md_nsmp) .eq. 0) then
                    call md_smp(mobj)
                end if
            end do
        case (1) ! end scheme
            call md_smp(mobj)
            !!
            do i = 1,md_nstep
                !!
                call md_run_p(2,mobj)
                call md_run_x(1,mobj)
                call md_run_p(2,mobj)
                call md_run_t(1,mobj)
                !!
                md_npass = md_npass + 1
                if(mod(i,md_nsmp) .eq. 0) then
                    call md_smp(mobj)
                end if
            end do
        case (2) ! side scheme
            call md_smp(mobj)
            !!
            do i = 1,md_nstep
                !!
                call md_run_t(2,mobj)
                call md_run_p(2,mobj)
                call md_run_x(1,mobj)
                call md_run_p(2,mobj)
                call md_run_t(2,mobj)
                !!
                md_npass = md_npass + 1
                if(mod(i,md_nsmp) .eq. 0) then
                    call md_smp(mobj)
                end if
            end do
    end select
    close(unit=out_unit)
    close(unit=ana_unit)
    !!
    open(unit=out_unit, file='end.rst', status='replace')
    do i=1,mobj%nb
        do j=1,md_bead
            write (out_unit,*) md_npass, i, j, mobj%a(i)%ks(j), mobj%a(i)%p(j), &
                mobj%a(i)%x(j),mobj%a(i)%fks(j),mobj%a(i)%fx(j)                     ! 8 term
        end do
    end do
    close(unit=out_unit)
end subroutine md_ctrl

subroutine md_smp(mobj)
implicit none
    type(mole), intent(in) :: mobj
    integer :: i,j
    !!
    do i=1,mobj%nb
        do j=1,md_bead
            write (out_unit,*) md_npass, i, j, mobj%a(i)%ks(j), mobj%a(i)%p(j), &
                mobj%a(i)%x(j),mobj%a(i)%fks(j),mobj%a(i)%fx(j)                     ! 8 term
        end do
        if (md_mod .ne. 1 .or. md_is_eqb .eq. 1 .or. md_npass > md_nstep/4) then
            if(md_methd.ne.2) then
                write (ana_unit,*) md_npass, i, sum( mobj%a(i)%x )/real(md_bead), &
                    mobj%a(i)%pe, mobj%a(i)%ke, mobj%a(i)%kev, mobj%a(i)%te, mobj%a(i)%tev  ! 8 term
            else
                write (ana_unit,*) md_npass, i, mobj%a(i)%ks(1), &
                    mobj%a(i)%pe, mobj%a(i)%ke, mobj%a(i)%kev, mobj%a(i)%te, mobj%a(i)%tev  ! 8 term
            end if
        end if
    end do
    !!
end subroutine md_smp
!------------------------------------------------------------------------
!subroutine check_equilibrium(mobj,pyd)
!implicit none
!    type(mole), intent(inout) :: mobj
!    type(pyrmd), intent(inout) :: pyd
!    integer :: oi, on, oj
    
!    oi = pyd%l(pyd%n-1)%c-1
!    on=size(pyd%l(pyd%n-1)%r)
!    oj = mod(oi-2+on, on ) + 1
    !!
!    if ( md_is_eqb .eq. 0) then
!        if ( comp( (/pyd%l(pyd%n-1)%r(oi),pyd%s(pyd%n-1)%r(oi)/),(/pyd%l(pyd%n-1)%r(oj),pyd%s(pyd%n-1)%r(oj)/), on) ) then
!            md_is_eqb = 1
!        end if
!    end if
!end check_equilibrium
!--------------------------------------------------------------------------

! an update of thermostat
! flag = 0 'lang' : langevin thermostat
! flag = 1 'ads' : andersen thermostat
! flag = 2 'nhc' : nose-hoover chain
subroutine md_run_t(n_dvd, mobj)
implicit none
    integer, intent(in) :: n_dvd
    type(mole), intent(inout) :: mobj
    integer :: i,j
    real(dp) :: et
    !!
    select case(md_thermo)
        case (0) ! langevin
            do i=1,mobj%nb
                do j=1,md_bead
                    call random_norm(et, rand_throw, 1.0_dp)
                    mobj%a(i)%p(j) = tmplst( 2*n_dvd - 1 ) * mobj%a(i)%p(j)  &
                        + sqrt( mobj%a(i)%m * masscoeff(j) * md_temp ) * tmplst( 2*n_dvd ) * et
                end do
            end do
        case (1) ! Anderson
            do i=1,mobj%nb
                do j=1,md_bead
                    call random_number(rand_throw)
                    if ( rand_throw < 1.- tmplst( 2*n_dvd - 1) ) then
                        call random_norm(et, rand_throw, 1.0_dp)
                        mobj%a(i)%p(j) = et * sqrt( mobj%a(i)%m * masscoeff(j) * md_temp )
                    end if
                end do
            end do
        case (2) ! NHC
            call NHChain(n_dvd, mobj)
    end select
end subroutine md_run_t

!! an update of position
subroutine md_run_x(n_dvd, mobj)
implicit none
    integer, intent(in) :: n_dvd
    type(mole), intent(inout) :: mobj
    integer :: i, j
    !!
    do i=1,mobj%nb
        do j=1,md_bead
            mobj%a(i)%ks(j) = mobj%a(i)%ks(j) + mobj%a(i)%p(j) * md_dtime / ( real(n_dvd) * mobj%a(i)%m * masscoeff(j) )
        end do
    end do
    !!
    select case (md_methd)
        case(0) ! primary-pimd
        case(1) ! staging-pimd
            call calc_x_stag(mobj)
            call calc_fx(mobj)
            call calc_fks_stag(mobj)
        case(2)! normal mode pimd
            call calc_x_norm(mobj)
            call calc_fx(mobj)
            call calc_fks_norm(mobj)
    end select
    call calc_pe(mobj)
end subroutine md_run_x

!! an update of momentum
subroutine md_run_p(n_dvd, mobj)
implicit none
    integer, intent(in) :: n_dvd
    type(mole), intent(inout) :: mobj
    integer :: i, j
    !!
    do i=1,mobj%nb
        do j=1,md_bead
            mobj%a(i)%p(j) = mobj%a(i)%p(j) - mobj%a(i)%fks(j) * md_dtime / real(n_dvd)
        end do
    end do
end subroutine md_run_p
!!
subroutine calc_fx(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i, j
    do i=1,mobj%nb
        do j=1, md_bead
            mobj%a(i)%fx(j) = Fpn2( 0.5_dp*mobj%a(i)%m,  mobj%a(i)%x(j) )
        end do
    end do
end subroutine calc_fx
!!
subroutine calc_pe(mobj)
implicit none
    type(mole), intent(inout) :: mobj 
    integer :: i,j
    do i=1,mobj%nb
        mobj%a(i)%pe = 0.0_dp
        mobj%a(i)%ke = 0.5_dp*md_bead*md_temp + 0.5_dp * mobj%a(i)%m * masscoeff(1)* md_bfreq2 * mobj%a(i)%ks(1)**2
        mobj%a(i)%kev = 0.0_dp
        do j=1,md_bead
            mobj%a(i)%pe = mobj%a(i)%pe + Vpn2( 0.5_dp*mobj%a(i)%m, mobj%a(i)%x(j))
            mobj%a(i)%ke = mobj%a(i)%ke - 0.5_dp * mobj%a(i)%m * masscoeff(j)* md_bfreq2 * mobj%a(i)%ks(j)**2
            mobj%a(i)%kev = mobj%a(i)%kev + 0.5_dp * mobj%a(i)%x(j) * mobj%a(i)%fx(j)
        end do
        mobj%a(i)%pe = mobj%a(i)%pe / real(md_bead)
        mobj%a(i)%kev = mobj%a(i)%kev / real(md_bead)
        mobj%a(i)%te = mobj%a(i)%pe + mobj%a(i)%ke
        mobj%a(i)%tev = mobj%a(i)%pe + mobj%a(i)%kev
    end do
end subroutine calc_pe
end module md_rout
