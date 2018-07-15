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
    
    !-- if needed, set Nose-Hoover chain settings
    if (md_thermo.eq.2) then
        call generate_NHC(mobj%nb, 3, 5)
    end if
    
    !-- set transform settings (staging transform / normal mode transform)
    call set_trsfrm(md_bead)
    
    !-- special for equilibrium, where from x to update other values
    if(md_mod.eq.1) then
        do i=1,mobj%nb
            call x2ks_stag(mobj%a(i)%ks(1::md_dim), mobj%a(i)%x(1::md_dim) )
            call x2ks_stag(mobj%a(i)%ks(2::md_dim), mobj%a(i)%x(2::md_dim) )
            call x2ks_stag(mobj%a(i)%ks(3::md_dim), mobj%a(i)%x(3::md_dim) )
        end do
        call calc_realfx(mobj)
    end if
    
    !-- open relative record files
    !open(unit=out_unit, file=traj_outfile, status='replace')
    !write(out_unit,*) 'nstep    ', 'atomidx ','nbead    ','ks   ','p    ','x    ','fks  ', 'fx   '
    open(unit=ana_unit, file=anal_file, status='replace')
    write(ana_unit,*) 'nstep    ', 'atomidx ','bar{x}    ','pe   ','ke    ', 'kev   ','tote    ','is_eqb  '
    
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
    close(unit=ana_unit)
    
    ! generate restart file
    open(unit=out_unit, file='end.rst', status='replace')
    do i=1,mobj%nb
        do j=1,mdn_bead
        ! 9 term, output trajectory information
            write (out_unit,*) md_npass, i, j, mobj%a(i)%ks(j), mobj%a(i)%p(j), &
                mobj%a(i)%x(j),mobj%a(i)%fks(j),mobj%a(i)%fx(j)
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
    	do d=1,md_dim
    	    write (ana_unit,*) md_npass, sum( mobj%a(i)%x(d::md_dim) )/real(md_bead)
    	end do
    end do
    !do i=1,mobj%nb
        !do j=1,md_bead
            ! 9 term, output trajectory information
            !write (trj_unit,*) md_npass, i, j, mobj%a(i)%ks(md_dim*j-2), mobj%a(i)%ks(md_dim*j-1), mobj%a(i)%ks(md_dim*j), &
            !    mobj%a(i)%p(md_dim*j-2), mobj%a(i)%p(md_dim*j-1), mobj%a(i)%p(md_dim*j)
            ! 9 term, output transformation information
            !write (trs_unit,*) md_npass, i, j, mobj%a(i)%x(md_dim*j-2), mobj%a(i)%x(md_dim*j-1), mobj%a(i)%x(md_dim*j), &
            !    mobj%a(i)%fx(md_dim*j-2), mobj%a(i)%fx(md_dim*j-1), mobj%a(i)%fx(md_dim*j)
        !end do
        ! 8 term, output energy information
        !write (erg_unit,*) md_npass, i, mobj%a(i)%pe, mobj%a(i)%ke, mobj%a(i)%kev, mobj%a(i)%te, mobj%a(i)%tev, md_iseqb
        !if(md_methd.ne.2) then
            ! 6 term, output else information, default
            !write (els_unit,*) md_npass, i, sum( mobj%a(i)%x(1::md_dim) )/real(md_bead), &
            !    sum( mobj%a(i)%x(2::md_dim) )/real(md_bead), sum( mobj%a(i)%x(3::md_dim) )/real(md_bead), md_iseqb
        !else
            ! 6 term, output else information, for normal mode
            !write (els_unit,*) md_npass, i, mobj%a(i)%ks(1), mobj%a(i)%ks(2), mobj%a(i)%ks(3), md_iseqb
        !end if
    !end do
end subroutine md_smp


! md_thermo = 0 'lang' : langevin thermostat
! md_thermo = 1 'ads' : andersen thermostat
! md_thermo = 2 'nhc' : nose-hoover chain
subroutine md_run_t(n_dvd, mobj)
implicit none
    integer, intent(in) :: n_dvd
    type(mole), intent(inout) :: mobj
    integer :: i,j
    real(dp) :: et
    !!
    select case(md_thermo)
        case (0)
            do i=1,mobj%nb
                do j=1,mdn_bead
                    call random_norm(et, rand_throw, 1.0_dp)
                    mobj%a(i)%p(j) = langc(n_dvd,1) * mobj%a(i)%p(j)  &
                        + sqrt( mobj%a(i)%m * masscoeff( int((j+md_dim-1)/md_dim) ) * md_temp ) * langc(n_dvd,2) * et
                end do
            end do
        case (1) ! Anderson
            do i=1,mobj%nb
                do j=1,mdn_bead
                    call random_number(rand_throw)
                    if ( rand_throw < 1.- langc(n_dvd,1) ) then
                        call random_norm(et, rand_throw, 1.0_dp)
                        mobj%a(i)%p(j) = et * sqrt( mobj%a(i)%m * masscoeff( int((j+md_dim-1)/md_dim) ) * md_temp )
                    end if
                end do
            end do
        case (2) ! NHC, not suit for high-dimension ???
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
        do j=1,mdn_bead
            mobj%a(i)%ks(j) = mobj%a(i)%ks(j) + mobj%a(i)%p(j) * &
                md_dtime / ( real(n_dvd) * mobj%a(i)%m * masscoeff( int((j+md_dim-1)/md_dim) ) )
        end do
    end do
    call md_update(mobj)
end subroutine md_run_x

!! an update of momentum
subroutine md_run_p(n_dvd, mobj)
implicit none
    integer, intent(in) :: n_dvd
    type(mole), intent(inout) :: mobj
    integer :: i, j
    !!
    do i=1,mobj%nb
        do j=1,mdn_bead
            mobj%a(i)%p(j) = mobj%a(i)%p(j) - mobj%a(i)%fks(j) * md_dtime / real(n_dvd)
        end do
    end do
end subroutine md_run_p
!!

subroutine md_update(mobj)
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i
    
    select case (md_methd)
        case(1)
            do i=1,mobj%nb
                call ks2x_stag(mobj%a(i)%x(1::md_dim), mobj%a(i)%ks(1::md_dim) )
                call ks2x_stag(mobj%a(i)%x(2::md_dim), mobj%a(i)%ks(2::md_dim) )
                call ks2x_stag(mobj%a(i)%x(3::md_dim), mobj%a(i)%ks(3::md_dim) )
            end do
            call calc_realfxh2o(mobj)
            do i=1,mobj%nb
                call fx2fks_stag(mobj%a(i)%fks(1::md_dim), mobj%a(i)%fx(1::md_dim), mobj%a(i)%ks(1::md_dim), mobj%a(i)%m )
                call fx2fks_stag(mobj%a(i)%fks(2::md_dim), mobj%a(i)%fx(2::md_dim), mobj%a(i)%ks(2::md_dim), mobj%a(i)%m )
                call fx2fks_stag(mobj%a(i)%fks(3::md_dim), mobj%a(i)%fx(3::md_dim), mobj%a(i)%ks(3::md_dim), mobj%a(i)%m )
            end do
        case(2)
            do i=1,mobj%nb
                call ks2x_norm(mobj%a(i)%x(1::md_dim), mobj%a(i)%ks(1::md_dim) )
                call ks2x_norm(mobj%a(i)%x(2::md_dim), mobj%a(i)%ks(2::md_dim) )
                call ks2x_norm(mobj%a(i)%x(3::md_dim), mobj%a(i)%ks(3::md_dim) )
            end do
            call calc_realfxh2o(mobj)
            do i=1,mobj%nb
                call fx2fks_norm(mobj%a(i)%fks(1::md_dim), mobj%a(i)%fx(1::md_dim), mobj%a(i)%ks(1::md_dim), mobj%a(i)%m )
                call fx2fks_norm(mobj%a(i)%fks(2::md_dim), mobj%a(i)%fx(2::md_dim), mobj%a(i)%ks(2::md_dim), mobj%a(i)%m )
                call fx2fks_norm(mobj%a(i)%fks(3::md_dim), mobj%a(i)%fx(3::md_dim), mobj%a(i)%ks(3::md_dim), mobj%a(i)%m )
            end do
    end select
    call calc_eg(mobj)
end subroutine md_update


subroutine calc_fx(mobj)    ! for 3-Dimension, it is the most tough to handle
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i, j
    do i=1,mobj%nb
        do j=1, mdn_bead
            mobj%a(i)%fx(j) = Fpn2( 0.5_dp*mobj%a(i)%m,  mobj%a(i)%x(j) )
        end do
    end do
end subroutine calc_fx


subroutine calc_eg(mobj)    ! 3-Dimension energy, how to solve
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i,j
    real(dp) :: scalar
    do i=1,mobj%nb
        mobj%a(i)%pe = 0.0_dp
        mobj%a(i)%ke = 0.5_dp*md_bead*md_temp
        mobj%a(i)%kev = 0.0_dp
        do j=1,md_bead
            call V3PN002(scalar, mobj%a(i)%x(md_dim*j-2:md_dim*j), 0.5_dp*mobj%a(i)%m )
            mobj%a(i)%pe = mobj%a(i)%pe + scalar
            mobj%a(i)%kev = mobj%a(i)%kev + 0.5_dp * (mobj%a(i)%x(md_dim*j-2) * mobj%a(i)%fx(md_dim*j-2) + &
                    mobj%a(i)%x(md_dim*j-1) * mobj%a(i)%fx(md_dim*j-1) + mobj%a(i)%x(md_dim*j) * mobj%a(i)%fx(md_dim*j) )
            if(j.eq.1) cycle
            mobj%a(i)%ke = mobj%a(i)%ke - 0.5_dp * mobj%a(i)%m * masscoeff(j)* md_bfreq2 * ( mobj%a(i)%ks(md_dim*j-2)**2 + &
                    mobj%a(i)%ks(md_dim*j-1)**2 + mobj%a(i)%ks(md_dim*j)**2 )
        end do
        mobj%a(i)%pe = mobj%a(i)%pe / real(md_bead)
        mobj%a(i)%kev = mobj%a(i)%kev / real(md_bead)
        mobj%a(i)%te = mobj%a(i)%pe + mobj%a(i)%ke
        mobj%a(i)%tev = mobj%a(i)%pe + mobj%a(i)%kev
    end do
end subroutine calc_eg

!!---------------------------------------------- for real system -----------------------------------------
subroutine calc_realfxh2o(mobj)    ! for 3-Dimension, it is the most tough to handle
use h2opes
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i, j
    real(dp) :: v,erg
    real(dp), dimension(3,3) :: dv
    real(dp), dimension(3,3,3,3) :: ddv
    real(dp), dimension(3,3) :: cart
    !!
    ! note: 1:H1, 2:O, 3:H2
    if(mobj%nb.ne.3) stop
    do i=1,mobj%nb
        mobj%a(i)%fx = 0.0_dp
    end do
    !!
    erg=0.0_dp
    do j=1,md_bead
        cart(1,:)=mobj%a(1)%x(md_dim*j-2:md_dim*j)
        cart(2,:)=mobj%a(2)%x(md_dim*j-2:md_dim*j)
        cart(3,:)=mobj%a(3)%x(md_dim*j-2:md_dim*j)
        !!
        call h2opot(v,dv,ddv,cart,1)
        !！
        mobj%a(1)%fx(md_dim*j-2:md_dim*j) = dv(1,:)
        mobj%a(2)%fx(md_dim*j-2:md_dim*j) = dv(2,:)
        mobj%a(3)%fx(md_dim*j-2:md_dim*j) = dv(3,:)
        erg=erg+v
    end do
    mobj%pes = erg / real(md_bead)
end subroutine calc_realfxh2o



subroutine calc_realfx(mobj)    ! for 3-Dimension, it is the most tough to handle
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
                call F3LJ001( vec3, mobj%a(i)%x(md_dim*j-2:md_dim*j), mobj%a(k)%x(md_dim*j-2:md_dim*j) )
                mobj%a(i)%fx(md_dim*j-2:md_dim*j) = mobj%a(i)%fx(md_dim*j-2:md_dim*j) + vec3
                mobj%a(k)%fx(md_dim*j-2:md_dim*j) = mobj%a(k)%fx(md_dim*j-2:md_dim*j) - vec3
            end do
        end do
    end do
end subroutine calc_realfx


!-- calculation values for analysis
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
