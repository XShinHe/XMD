!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 08
!--- Acknowledgement to Liu group, of PKU
!*


module pimd_plus
use MyDef
use AM_script
use myFF
use md_info
implicit none
    real(dp), dimension(:,:), allocatable, private :: norm_OM
    real(dp), dimension(:), allocatable, public :: masscoeff
    
contains
    subroutine set_trsfrm(bead)
        integer, intent(in) :: bead
        integer :: i,j
        allocate( masscoeff( bead ) )
        
        if (md_ipimd.eq.1) then
            masscoeff(1) = 1.0_dp
            do i=2,bead
                masscoeff(i) = real(i)/(real(i) - 1.0_dp)
            end do
        else
            if( mod(bead,2) .ne. 0) then
                call send_err('norm mode pimd, beads should be even!')
            end if
            !!
            allocate( norm_OM(bead, bead) )
            masscoeff(1) = real(bead,kind=dp)
            norm_OM(1,:) = 1.0_dp/ DSQRT( real(bead,kind=dp) )
            do i=2,bead/2
                masscoeff(2*i-2) = 2.0_dp*( 1 - DCOS(twopi * real(i-1)/real(bead,kind=dp) ) )*real(bead,kind=dp)
                masscoeff(2*i-1) = 2.0_dp*( 1 - DCOS(twopi * real(i-1)/real(bead,kind=dp) ) )*real(bead,kind=dp)
                do j=1,bead
                    norm_OM(2*i-2,j) = DCOS( twopi* real(i-1)*real(j-1)/real(bead,kind=dp) ) * DSQRT(2.0_dp/real(bead,kind=dp) )
                    norm_OM(2*i-1,j) = -DSIN( twopi* real(i-1)*real(j-1)/real(bead,kind=dp) ) * DSQRT(2.0_dp/ real(bead,kind=dp) )
                end do
            end do
            masscoeff(bead) = 4.0_dp*real(bead,kind=dp)
            norm_OM(bead,1::2) = 1.0_dp/ DSQRT( real(bead,kind=dp) )
            norm_OM(bead,2::2) = -1.0_dp/ DSQRT( real(bead,kind=dp) )
        end if
        !!
    end subroutine set_trsfrm

	!-- normal mode transformation
    subroutine calc_ks_norm(mobj)
    implicit none
        type(mole), intent(inout) :: mobj
        integer :: i,j,k
        do i=1,md_nsum
            do j=1,md_bead
                mobj%v(i,1,j) = 0.0_dp
                do k=1,md_bead
                    mobj%v(i,1,j) = mobj%v(i,1,j) + mobj%v(i,2,k)*norm_OM(j,k)
                end do
                mobj%v(i,1,j) = mobj%v(i,1,j) / DSQRT(real(md_bead,kind=dp))
            end do
        end do
    end subroutine calc_ks_norm

    subroutine calc_x_norm(mobj)
    implicit none
        type(mole), intent(inout) :: mobj
        integer :: i, j, k
        !!
        do i=1,md_nsum
            do j=1,md_bead
                mobj%v(i,2,j) = 0.0_dp
                do k=1, md_bead
                    mobj%v(i,2,j) = mobj%v(i,2,j) + mobj%v(i,1,k) * norm_OM(k,j)
                end do
                mobj%v(i,2,j) = mobj%v(i,2,j) * DSQRT(real(md_bead,kind=dp))
            end do
        end do
    end subroutine calc_x_norm

    subroutine calc_fks_norm(mobj)
    implicit none
        type(mole), intent(inout) :: mobj
        integer :: i, j, k
        !!
        do i=1,md_nsum
            do j=1,md_bead,1
                mobj%v(i,5,j) = 0.0_dp
                do k=1,md_bead
                    mobj%v(i,5,j) = mobj%v(i,5,j) + mobj%v(i,4,k) * norm_OM(j,k)
                end do
                mobj%v(i,5,j) = mobj%v(i,5,j)
            end do
            mobj%v(i,5,1) = mobj%v(i,5,1) / DSQRT( real(md_bead,kind=dp) )
            do j=2,md_bead
                mobj%v(i,5,j) = mobj%v(i,5,j)/ DSQRT(real(md_bead,kind=dp)) + masscoeff(j) &
                    * mobj%a(i)%m * md_bfreq2 * mobj%v(i,1,j)
            end do
        end do
    end subroutine calc_fks_norm
	
	!-- staging transformation
    subroutine calc_ks_stag(mobj)
    implicit none
        type(mole), intent(inout) :: mobj
        integer :: i,j
        do i=1,md_nsum
            mobj%v(i,1,1) = mobj%v(i,2,1)
            if(md_bead.eq.1) cycle
            mobj%v(i,1,md_bead) = mobj%v(i,2,md_bead) - mobj%v(i,2,1)
            do j=2,md_bead-1,1
                mobj%v(i,1,j) = mobj%v(i,2,j) - ( mobj%v(i,2,j+1) * real(j - 1) + mobj%v(i,2,1) ) / real(j)
            end do
        end do
    end subroutine calc_ks_stag

    subroutine calc_x_stag(mobj)
    implicit none
        type(mole), intent(inout) :: mobj
        integer :: i, j
        !!
        do i=1,md_nsum
            mobj%v(i,2,1) = mobj%v(i,1,1)
            if(md_bead.eq.1) cycle
            mobj%v(i,2,md_bead) = mobj%v(i,1,md_bead) + mobj%v(i,1,1)
            do j=md_bead-1,2,-1
                mobj%v(i,2,j) = mobj%v(i,1,j) + ( real( j - 1 ) * mobj%v(i,2,j+1) + mobj%v(i,1,1) ) / real(j)
            end do
        end do
    end subroutine calc_x_stag

    subroutine calc_fks_stag(mobj)
    implicit none
        type(mole), intent(inout) :: mobj
        integer :: i, j
        !!
        do i=1,md_nsum
            mobj%v(i,5,1) = 0.0_dp
            do j=1,md_bead
                mobj%v(i,5,1) = mobj%v(i,5,1) + mobj%v(i,4,j)
            end do
            do j=2,md_bead,1
                mobj%v(i,5,j) = mobj%v(i,4,j) + mobj%v(i,5,j-1) * real(j-2)/real(j-1) 
            end do
            mobj%v(i,5,1) = mobj%v(i,5,1)/real(md_bead,kind=dp)
            do j=2,md_bead
                mobj%v(i,5,j) = mobj%v(i,5,j)/real(md_bead,kind=dp) + masscoeff(j) * mobj%a(i)%m * md_bfreq2 * mobj%v(i,1,j)
            end do
        end do
    end subroutine calc_fks_stag

	!-- public function
    subroutine calc_fx_ofmethd(mobj)
    implicit none
        type(mole), intent(inout) :: mobj
        integer :: i, j
        do i=1,md_nsum
            do j=1, md_bead
                mobj%v(i,4,j) = Fpn2( 0.5_dp,  mobj%v(i,2,j) )
            end do
        end do
    end subroutine calc_fx_ofmethd
        
end module pimd_plus



