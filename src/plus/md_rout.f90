!*
! ---Copyright by--- He Xin <1500011805@pku.edu.cn>
!*

! --------------------------- mdbsc module -------------------------
module md_rout
use MyDef
use AM_script
use md_info
use myFF
use thermo_plus
use methd_plus
use elim_plus
implicit none

contains
subroutine md_ctrl(mobj)
use h2opes
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i,j
    !!
    real(dp), dimension(3) :: geom=(/1.81_dp,1.81_dp,1.8239_dp/)
    real(dp), dimension(3,3) :: cart
    !!
    !! if use HNC, initialize NHC component
    if (md_thermo.eq.2) then
        call generate_NHC(mobj%nb, 3, 5) 
    end if
    !! use transform method, initialize related component
    call set_trsfrm(md_bead)
    
    ! open trajectory record file
    !open(unit=trj_unit, file=traj_file, status='replace')
    !write(trj_unit,*) 'nstep    ', 'atomidx ','nbead    ','ks1   ','ks2    ','ks3    ','p1  ', 'p2   ', 'p3     ' 
    ! open transform record file
    !open(unit=trs_unit, file=trsf_file, status='replace')
    !write(trs_unit,*) 'nstep    ', 'atomidx ','nbead    ','x1   ','x2    ','x3    ','fx1  ', 'fx2   ', 'fx3     ' 
    ! open energy record file
    !open(unit=erg_unit, file=enrg_file, status='replace')
    !write(erg_unit,*) 'nstep    ', 'atomidx     ','pe   ','ke    ', 'kev   ','te    ','tev     ' 
    ! open elsesomthing record file, cedroid information
    !open(unit=els_unit, file=else_file, status='replace')
    !write(els_unit,*) 'nstep    ', 'atomidx ','aver_x1    ','aver_x2    ','aver_x3    ' 
    !！
    open(unit=ana_unit, file=anal_file, status='replace')
    write(ana_unit,*) 'nstep    ','(H-O)1x    ','(H-O)1y    ','(H-O)1z     ','(H-O)2x    ','(H-O)2y    ','(H-O)2z     ' 
    !!
    ! correct of h2o molecule, if use mod=1(e)
    if(md_mod.eq.1) then
        call convert_geom_to_cart(cart,geom)
        do j=1,md_bead
            mobj%a(1)%x(3*j-2:3*j)=cart(1,:)
            mobj%a(2)%x(3*j-2:3*j)=cart(2,:)
            mobj%a(3)%x(3*j-2:3*j)=cart(3,:)
        end do
        do i=1,mobj%nb
            call x2ks_stag(mobj%a(i)%ks(1::3), mobj%a(i)%x(1::3) )
            call x2ks_stag(mobj%a(i)%ks(2::3), mobj%a(i)%x(2::3) )
            call x2ks_stag(mobj%a(i)%ks(3::3), mobj%a(i)%x(3::3) )
        end do
        call calc_realfxh2o(mobj)
    end if
    
    !! calculate for initial configurations
    !! call md_update(mobj) ! we needn't
    call set_elim(mobj)
    !!
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
                call elim_trs(mobj)
                call md_run_p(2,mobj)
                call elim_rot(mobj)
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
    close(unit=trj_unit)
    close(unit=erg_unit)
    close(unit=trs_unit)
    close(unit=els_unit)
    
    ! generate restart file
    open(unit=trj_unit, file='end.rst', status='replace')
    do i=1,mobj%nb
        do j=1,md_bead
        ! 9 term, output trajectory information
            write (trj_unit,*) md_npass, i, j, mobj%a(i)%ks(3*j-2), mobj%a(i)%ks(3*j-1), mobj%a(i)%ks(3*j), &
                mobj%a(i)%p(3*j-2), mobj%a(i)%p(3*j-1), mobj%a(i)%p(3*j)
        end do
    end do
    close(unit=trj_unit)
end subroutine md_ctrl

subroutine md_smp(mobj)
implicit none
    type(mole), intent(in) :: mobj
    integer :: i,j
    !!
    write (ana_unit,*) md_npass, (sum( mobj%a(1)%x(1::3) )/real(md_bead)-sum( mobj%a(2)%x(1::3) )/real(md_bead)),&
                 (sum( mobj%a(1)%x(2::3) )/real(md_bead)-sum( mobj%a(2)%x(2::3) )/real(md_bead)),&
                  (sum( mobj%a(1)%x(3::3) )/real(md_bead)-sum( mobj%a(2)%x(3::3) )/real(md_bead)),&
                   (sum( mobj%a(3)%x(1::3) )/real(md_bead)-sum( mobj%a(2)%x(1::3) )/real(md_bead)),&
                    (sum( mobj%a(3)%x(2::3) )/real(md_bead)-sum( mobj%a(2)%x(2::3) )/real(md_bead)),&
                     (sum( mobj%a(3)%x(3::3) )/real(md_bead)-sum( mobj%a(2)%x(3::3) )/real(md_bead))
    !do i=1,mobj%nb
        !do j=1,md_bead
            ! 9 term, output trajectory information
            !write (trj_unit,*) md_npass, i, j, mobj%a(i)%ks(3*j-2), mobj%a(i)%ks(3*j-1), mobj%a(i)%ks(3*j), &
            !    mobj%a(i)%p(3*j-2), mobj%a(i)%p(3*j-1), mobj%a(i)%p(3*j)
            ! 9 term, output transformation information
            !write (trs_unit,*) md_npass, i, j, mobj%a(i)%x(3*j-2), mobj%a(i)%x(3*j-1), mobj%a(i)%x(3*j), &
            !    mobj%a(i)%fx(3*j-2), mobj%a(i)%fx(3*j-1), mobj%a(i)%fx(3*j)
        !end do
        ! 8 term, output energy information
        !write (erg_unit,*) md_npass, i, mobj%a(i)%pe, mobj%a(i)%ke, mobj%a(i)%kev, mobj%a(i)%te, mobj%a(i)%tev, md_iseqb
        !if(md_methd.ne.2) then
            ! 6 term, output else information, default
            !write (els_unit,*) md_npass, i, sum( mobj%a(i)%x(1::3) )/real(md_bead), &
            !    sum( mobj%a(i)%x(2::3) )/real(md_bead), sum( mobj%a(i)%x(3::3) )/real(md_bead), md_iseqb
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
                do j=1,md3_bead
                    call random_norm(et, rand_throw, 1.0_dp)
                    mobj%a(i)%p(j) = langc(n_dvd,1) * mobj%a(i)%p(j)  &
                        + sqrt( mobj%a(i)%m * masscoeff( int((j+2)/3) ) * md_temp ) * langc(n_dvd,2) * et
                end do
            end do
        case (1) ! Anderson
            do i=1,mobj%nb
                do j=1,md3_bead
                    call random_number(rand_throw)
                    if ( rand_throw < 1.- langc(n_dvd,1) ) then
                        call random_norm(et, rand_throw, 1.0_dp)
                        mobj%a(i)%p(j) = et * sqrt( mobj%a(i)%m * masscoeff( int((j+2)/3) ) * md_temp )
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
        do j=1,md3_bead
            mobj%a(i)%ks(j) = mobj%a(i)%ks(j) + mobj%a(i)%p(j) * &
                md_dtime / ( real(n_dvd) * mobj%a(i)%m * masscoeff( int((j+2)/3) ) )
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
        do j=1,md3_bead
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
                call ks2x_stag(mobj%a(i)%x(1::3), mobj%a(i)%ks(1::3) )
                call ks2x_stag(mobj%a(i)%x(2::3), mobj%a(i)%ks(2::3) )
                call ks2x_stag(mobj%a(i)%x(3::3), mobj%a(i)%ks(3::3) )
            end do
            call calc_realfxh2o(mobj)
            do i=1,mobj%nb
                call fx2fks_stag(mobj%a(i)%fks(1::3), mobj%a(i)%fx(1::3), mobj%a(i)%ks(1::3), mobj%a(i)%m )
                call fx2fks_stag(mobj%a(i)%fks(2::3), mobj%a(i)%fx(2::3), mobj%a(i)%ks(2::3), mobj%a(i)%m )
                call fx2fks_stag(mobj%a(i)%fks(3::3), mobj%a(i)%fx(3::3), mobj%a(i)%ks(3::3), mobj%a(i)%m )
            end do
        case(2)
            do i=1,mobj%nb
                call ks2x_norm(mobj%a(i)%x(1::3), mobj%a(i)%ks(1::3) )
                call ks2x_norm(mobj%a(i)%x(2::3), mobj%a(i)%ks(2::3) )
                call ks2x_norm(mobj%a(i)%x(3::3), mobj%a(i)%ks(3::3) )
            end do
            call calc_realfxh2o(mobj)
            do i=1,mobj%nb
                call fx2fks_norm(mobj%a(i)%fks(1::3), mobj%a(i)%fx(1::3), mobj%a(i)%ks(1::3), mobj%a(i)%m )
                call fx2fks_norm(mobj%a(i)%fks(2::3), mobj%a(i)%fx(2::3), mobj%a(i)%ks(2::3), mobj%a(i)%m )
                call fx2fks_norm(mobj%a(i)%fks(3::3), mobj%a(i)%fx(3::3), mobj%a(i)%ks(3::3), mobj%a(i)%m )
            end do
    end select
    call calc_eg(mobj)
end subroutine md_update


subroutine calc_fx(mobj)    ! for 3-Dimension, it is the most tough to handle
implicit none
    type(mole), intent(inout) :: mobj
    integer :: i, j
    do i=1,mobj%nb
        do j=1, md3_bead
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
            call V3PN002(scalar, mobj%a(i)%x(3*j-2:3*j), 0.5_dp*mobj%a(i)%m )
            mobj%a(i)%pe = mobj%a(i)%pe + scalar
            mobj%a(i)%kev = mobj%a(i)%kev + 0.5_dp * (mobj%a(i)%x(3*j-2) * mobj%a(i)%fx(3*j-2) + &
                    mobj%a(i)%x(3*j-1) * mobj%a(i)%fx(3*j-1) + mobj%a(i)%x(3*j) * mobj%a(i)%fx(3*j) )
            if(j.eq.1) cycle
            mobj%a(i)%ke = mobj%a(i)%ke - 0.5_dp * mobj%a(i)%m * masscoeff(j)* md_bfreq2 * ( mobj%a(i)%ks(3*j-2)**2 + &
                    mobj%a(i)%ks(3*j-1)**2 + mobj%a(i)%ks(3*j)**2 )
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
        cart(1,:)=mobj%a(1)%x(3*j-2:3*j)
        cart(2,:)=mobj%a(2)%x(3*j-2:3*j)
        cart(3,:)=mobj%a(3)%x(3*j-2:3*j)
        !!
        call h2opot(v,dv,ddv,cart,1)
        !！
        mobj%a(1)%fx(3*j-2:3*j) = dv(1,:)
        mobj%a(2)%fx(3*j-2:3*j) = dv(2,:)
        mobj%a(3)%fx(3*j-2:3*j) = dv(3,:)
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
                call F3LJ001( vec3, mobj%a(i)%x(3*j-2:3*j), mobj%a(k)%x(3*j-2:3*j) )
                mobj%a(i)%fx(3*j-2:3*j) = mobj%a(i)%fx(3*j-2:3*j) + vec3
                mobj%a(k)%fx(3*j-2:3*j) = mobj%a(k)%fx(3*j-2:3*j) - vec3
            end do
        end do
    end do
end subroutine calc_realfx


subroutine calc_realeg(mobj)    ! 3-Dimension energy, how to solve
implicit none
    type(mole), intent(inout) :: mobj
    real(dp) :: scalar
    integer :: i, j, k
    do i=1,mobj%nb
        mobj%a(i)%pe = 0.0_dp
        mobj%a(i)%ke = 0.5_dp*md3_bead*md_temp
        mobj%a(i)%kev = 0.0_dp
    end do
    !!
    do i=1,mobj%nb
        do j=1,md_bead
            do k=i+1, mobj%nb
                call V3LJ001(scalar, mobj%a(i)%x(3*j-2:3*j), mobj%a(k)%x(3*j-2:3*j) )
                mobj%a(i)%pe = mobj%a(i)%pe + scalar
                mobj%a(k)%pe = mobj%a(k)%pe + scalar
            end do
            mobj%a(i)%kev = mobj%a(i)%kev + 0.5_dp * (mobj%a(i)%x(3*j-2) * mobj%a(i)%fx(3*j-2) + &
                    mobj%a(i)%x(3*j-1) * mobj%a(i)%fx(3*j-1) + mobj%a(i)%x(3*j) * mobj%a(i)%fx(3*j) )
            if(j.eq.1) cycle
            mobj%a(i)%ke = mobj%a(i)%ke - 0.5_dp * mobj%a(i)%m * masscoeff(j)* md_bfreq2 *( mobj%a(i)%ks(3*j-2)**2 + &
                    mobj%a(i)%ks(3*j-1)**2 + mobj%a(i)%ks(3*j)**2 )
        end do
        mobj%a(i)%pe = mobj%a(i)%pe / real(md_bead)
        mobj%a(i)%kev = mobj%a(i)%kev / real(md_bead)
        mobj%a(i)%te = mobj%a(i)%pe + mobj%a(i)%ke
        mobj%a(i)%tev = mobj%a(i)%pe + mobj%a(i)%kev
    end do
!!-------------------------------------- molecule properties --------------------------------
!    mobj%pe=0
!    mobj%ke=0
!    mobj%kev=0
!    do i=1,mobj%nb
!        mobj%pe=mobj%pe + mobj%a(i)%pe
!        mobj%ke=mobj%ke + mobj%a(i)%ke
!        mobj%kev=mobj%kev + mobj%a(i)%kev
!    end do
!    mobj%pe = mobj%pe /2.0_dp       ! for interaction is mutual
!    mobj%ke = mobj%ke /2.0_dp
!    mobj%kev = mobj%kev /2.0_dp
!    mobj%te = mobj%pe + mobj%ke
!    mobj%tev = mobj%pe + mobj%kev     
end subroutine calc_realeg

end module md_rout



