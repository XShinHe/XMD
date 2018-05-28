! --- Copyright by He Xin <hx0824916@pku.edu.cn> ---

module methd_plus
use MyDef
use AM_script
use md_info
implicit none
    real(dp), dimension(:,:), allocatable, private :: norm_OM
    real(dp), dimension(:), allocatable, public :: masscoeff
    
contains
    subroutine set_trsfrm(bead)
        integer, intent(in) :: bead
        integer :: i,j
        allocate( masscoeff( bead ) )
        
        if (md_methd.eq.1) then
            masscoeff(1) = 1.0_dp
            do i=2,bead
                masscoeff(i) = real(i)/(real(i) - 1.0_dp)
            end do
        else
            if( mod(bead,2) .ne. 0) then
                call send_err('norm pimd, beads should be even!')
            end if
            !!
            allocate( norm_OM(bead, bead) )
            masscoeff(1) = real(bead,kind=dp)
            norm_OM(1,:) = 1.0_dp/ DSQRT( real(bead,kind=dp) )
            do i=2,bead/2
                masscoeff(2*i-2) = 2.0_dp*( 1 - DCOS(twopi * real(i-1)/real(bead,kind=dp) ) )*real(bead,kind=dp)
                masscoeff(2*i-1) = 2.0_dp*( 1 - DCOS(twopi * real(i-1)/real(bead,kind=dp) ) )*real(bead,kind=dp)
                do j=1,bead
                    norm_OM(2*i-2,j) = DCOS( twopi* real(i-1)*real(j-1)/real(bead,kind=dp) ) *&
                         DSQRT( 2.0_dp/real(bead,kind=dp) )
                    norm_OM(2*i-1,j) = -DSIN( twopi* real(i-1)*real(j-1)/real(bead,kind=dp) ) *&
                         DSQRT( 2.0_dp/real(bead,kind=dp) )
                end do
            end do
            masscoeff(bead) = 4.0_dp*real(bead,kind=dp)
            norm_OM(bead,1::2) = 1.0_dp/ DSQRT( real(bead,kind=dp) )
            norm_OM(bead,2::2) = -1.0_dp/ DSQRT( real(bead,kind=dp) )
        end if
    end subroutine set_trsfrm

!! ----------------------------------------- normal mode transformation ---------------------------
    subroutine x2ks_norm(ks, x)
    implicit none
        real(dp), dimension(:), intent(inout) :: ks,x
        integer :: j,k
        do j=1,md_bead
            ks(j) = 0.0_dp
            do k=1,md_bead
                ks(j) = ks(j) + x(k) * norm_OM(j,k)
            end do
            ks(j) = ks(j) / DSQRT(real(md_bead,kind=dp))
        end do
    end subroutine x2ks_norm
!!
    subroutine ks2x_norm(x,ks)
    implicit none
        real(dp), dimension(:), intent(inout) :: ks,x
        integer :: j, k
        !!
        do j=1,md_bead
            x(j) = 0.0_dp
            do k=1, md_bead
                x(j) = x(j) + ks(k) * norm_OM(k,j)
            end do
            x(j) = x(j) * DSQRT(real(md_bead,kind=dp))
        end do
    end subroutine ks2x_norm

    subroutine fx2fks_norm(fks, fx, ks, mi)
    implicit none
        real(dp), dimension(:), intent(inout) :: fks,fx,ks
        real(dp) :: mi
        integer :: j, k
        !!
        do j=1,md_bead,1
            fks(j) = 0.0_dp
            do k=1,md_bead
                fks(j) = fks(j) + fx(k) * norm_OM(j,k)
            end do
            fks(j) = fks(j) * DSQRT( real(md_bead,kind=dp) )
        end do
        fks(1) = fks(1)/DSQRT(real(md_bead,kind=dp))
        do j=2,md_bead
            fks(j) = fks(j)/DSQRT(real(md_bead,kind=dp)) + masscoeff(j) &
                * mi * md_bfreq2 * ks(j)
        end do
    end subroutine fx2fks_norm
!!------------------------------------------------ staging transformation -----------------------------------
    subroutine x2ks_stag(ks,x)
    implicit none
        real(dp), dimension(:), intent(inout) :: ks,x
        integer :: j
        ks(1) = x(1)
        if(md_bead.eq.1) return
        ks(md_bead) = x(md_bead) - x(1)
        do j=2,md_bead-1,1
            ks(j) = x(j) - ( x(j+1) * real(j - 1) + x(1) ) / real(j)
        end do
    end subroutine x2ks_stag
!!
    subroutine ks2x_stag(x,ks)
    implicit none
        real(dp), dimension(:), intent(inout) :: ks,x
        integer :: j
        !!
        x(1) = ks(1)
        if(md_bead.eq.1) return
        x( md_bead ) = ks( md_bead ) + x(1)
        do j=md_bead-1,2,-1
            x(j) = ks(j) + ( real( j - 1 ) * x(j+1) + ks(1) ) / real(j)
        end do
    end subroutine ks2x_stag
!!
    subroutine fx2fks_stag(fks, fx, ks, mi)
    implicit none
        real(dp), dimension(:), intent(inout) :: fks,fx,ks
        real(dp) :: mi
        integer :: j
        !!
        fks(1) = 0.0_dp
        do j=1,md_bead
            fks(1) = fks(1) + fx(j)
        end do
        do j=2,md_bead,1
            fks(j) = fx(j) + fks(j-1) * real(j-2)/real(j-1) 
        end do
        fks(1) = fks(1)/real(md_bead,kind=dp)
        do j=2,md_bead
            fks(j) = fks(j)/real(md_bead,kind=dp) + masscoeff(j) * mi * md_bfreq2 * ks(j)
        end do
    end subroutine fx2fks_stag

end module methd_plus


