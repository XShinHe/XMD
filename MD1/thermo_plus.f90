!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 08
!--- Acknowledgement to Liu group, of PKU
!*


module thermo_plus
use MyDef
use AM_script
use pimd_plus
use md_info
implicit none
    integer, private :: N_nhc=2 ! for default
    integer, private :: N_respa=1 ! for default
    integer, private :: N_sy=7
    real(dp), dimension(7), private :: Wsy=(/0.784513610477560_dp,0.235573213359357_dp,-1.17767998417887_dp,&
                                        1.3151863206839063_dp,-1.17767998417887_dp,0.235573213359357_dp,0.784513610477560_dp/)
    real(dp), dimension(7) :: delta
    real(dp), dimension(:,:,:), allocatable :: nx, np, nG
    real(dp), dimension(:), allocatable :: nQ
contains
    subroutine generate_NHC(natom, nnhc, nrespa)
        integer, intent(in) :: natom, nnhc, nrespa
        integer :: alpha
        if(nnhc>2) then
            N_nhc=nnhc
        end if
        N_respa = nrespa
        allocate( nx(N_nhc, natom, md_bead), np(N_nhc, natom, md_bead), nG(N_nhc, natom, md_bead), nQ(N_nhc) )
        nx=0.0_dp
        np=0.0_dp
        nG=0.0_dp
        if(md_coeff < 20*md_dtime) then
            md_coeff = 20*md_dtime
        end if
        ! optimal choice
        nQ = md_temp*md_coeff**2
        nQ(1) = natom * md_temp*md_coeff**2
        do alpha=1,N_sy
            delta(alpha) = Wsy(alpha) * md_dtime / N_respa
        end do
    end subroutine generate_NHC
    
    subroutine NHChainP(n_dvd, mobj)
        type(mole), intent(inout) :: mobj
        integer, intent(in) :: n_dvd
        integer :: h, i, j, k, alpha
        do i=1,md_nsum
            do j=1,md_bead
                do k=1,N_respa
                    do alpha=1,N_sy
                        np(N_nhc,i,j) = np(N_nhc,i,j) + ( np(N_nhc-1,i,j)*np(N_nhc-1,i,j)/nQ(N_nhc-1) - md_temp ) &
                                        *(delta(alpha) / real(2*n_dvd))
                        do h=N_nhc-1,2,-1
                            np(h,i,j) = np(h,i,j) - np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd)
                            np(h,i,j) = np(h,i,j) + ( np(h-1,i,j)*np(h-1,i,j)/nQ(h-1) - md_temp ) * (delta(alpha) / real(2*n_dvd))
                            np(h,i,j) = np(h,i,j) - np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd)
                        end do
                        np(1,i,j) = np(1,i,j) - np(2,i,j)/nQ(2) * delta(alpha) / real(4*n_dvd) 
                        np(1,i,j) = np(1,i,j) + ( mobj%a(i)%p(j)**2/(mobj%a(i)%m * masscoeff(j)) - md_temp ) &
                                     * (delta(alpha) / real(2*n_dvd))
                        np(1,i,j) = np(1,i,j) - np(2,i,j)/nQ(2) * delta(alpha) / real(4*n_dvd) 
                        !!
                        do h=1,N_nhc
                            nx(h,i,j) = nx(h,i,j) + np(h,i,j) / nQ(h)* delta(alpha)/real(n_dvd)
                        end do
                        mobj%a(i)%p(j) = mobj%a(i)%p(j) - np(1,i,j)/nQ(1)*(delta(alpha)/real(n_dvd))
                        print *,'i j k a p p1', i, j, k, alpha, mobj%a(i)%p(j), np(1,i,j)
                        !!
                        np(1,i,j) = np(1,i,j) - np(2,i,j)/nQ(2) * delta(alpha) / real(4*n_dvd)
                        np(1,i,j) = np(1,i,j) + ( mobj%a(i)%p(j)**2/(mobj%a(i)%m* masscoeff(j)) - md_temp ) &
                                    * (delta(alpha) / real(2*n_dvd))
                        np(1,i,j) = np(1,i,j) - np(2,i,j)/nQ(2) * delta(alpha) / real(4*n_dvd)
                        do h=2,N_nhc-1,1
                            np(h,i,j) = np(h,i,j) - np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd)
                            np(h,i,j) = np(h,i,j) + ( np(h-1,i,j)*np(h-1,i,j)/nQ(h-1) - md_temp ) * (delta(alpha) / real(2*n_dvd))
                            np(h,i,j) = np(h,i,j) - np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd)
                        end do
                        np(N_nhc,i,j) = np(N_nhc,i,j) + ( np(N_nhc-1,i,j)*np(N_nhc-1,i,j)/nQ(N_nhc-1) - md_temp ) &
                                        *(delta(alpha) / real(2*n_dvd))
                    end do
                end do
            end do
        end do
    end subroutine NHChainP
    
    subroutine NHChain(n_dvd, mobj)
        type(mole), intent(inout) :: mobj
        integer, intent(in) :: n_dvd
        integer :: h, i, j, k, alpha
        do i=1,md_nsum
            do j=1,md_bead
                do k=1,N_respa
                    do alpha=1,N_sy
                        np(N_nhc,i,j) = np(N_nhc,i,j) + ( np(N_nhc-1,i,j)*np(N_nhc-1,i,j)/nQ(N_nhc-1) - md_temp ) &
                                        *(delta(alpha) / real(2*n_dvd))
                        do h=N_nhc-1,2,-1
                            np(h,i,j) = np(h,i,j) * EXP( -np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd) )
                            np(h,i,j) = np(h,i,j) + ( np(h-1,i,j)*np(h-1,i,j)/nQ(h-1) - md_temp ) * (delta(alpha) / real(2*n_dvd))
                            np(h,i,j) = np(h,i,j) * EXP( -np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd) )
                        end do
                        np(1,i,j) = np(1,i,j) * EXP( -np(2,i,j)/nQ(2) * delta(alpha) / real(4*n_dvd) )
                        np(1,i,j) = np(1,i,j) + ( mobj%a(i)%p(j)**2/(mobj%a(i)%m* masscoeff(j)) - md_temp ) &
                                    * (delta(alpha) / real(2*n_dvd))
                        np(1,i,j) = np(1,i,j) * EXP( -np(2,i,j)/nQ(2) * delta(alpha) / real(4*n_dvd) )
                        !!
                        do h=1,N_nhc
                            nx(h,i,j) = nx(h,i,j) + np(h,i,j) / nQ(h)* delta(alpha)/real(n_dvd)
                        end do
                        mobj%a(i)%p(j) = mobj%a(i)%p(j) * EXP( -np(1,i,j)/nQ(1)*delta(alpha)/real(n_dvd) )
                        !!
                        np(1,i,j) = np(1,i,j) * EXP( -np(2,i,j)/nQ(2) * delta(alpha) / real(4*n_dvd) )
                        np(1,i,j) = np(1,i,j) + ( mobj%a(i)%p(j)**2/(mobj%a(i)%m* masscoeff(j)) - md_temp ) & 
                                    * (delta(alpha) / real(2*n_dvd))
                        np(1,i,j) = np(1,i,j) * EXP( -np(2,i,j)/nQ(2) * delta(alpha) / real(4*n_dvd) )
                        do h=2,N_nhc-1,1
                            np(h,i,j) = np(h,i,j) * EXP( -np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd) )
                            np(h,i,j) = np(h,i,j) + ( np(h-1,i,j)*np(h-1,i,j)/nQ(h-1) - md_temp ) * (delta(alpha) / real(2*n_dvd))
                            np(h,i,j) = np(h,i,j) * EXP( -np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd) )
                        end do
                        np(N_nhc,i,j) = np(N_nhc,i,j) + ( np(N_nhc-1,i,j)*np(N_nhc-1,i,j)/nQ(N_nhc-1) - md_temp ) &
                                        *(delta(alpha) / real(2*n_dvd))
                    end do
                end do
            end do
        end do
    end subroutine NHChain
    
    subroutine NHChain_partG(n_dvd, mobj)
        type(mole), intent(inout) :: mobj
        integer, intent(in) :: n_dvd
        integer :: h, i, j, k, alpha
        do i=1,md_nsum
            do j=1,md_bead
                do k=1,N_respa
                    do alpha=1,N_sy
                        nG(1,i,j) = mobj%a(i)%p(j)**2/(mobj%a(i)%m* masscoeff(j)) - md_temp
                        do h=2,N_nhc-1,1
                            nG(h,i,j) = np(h-1,i,j)*np(h-1,i,j)/nQ(h-1) - md_temp
                        end do
                        !!
                        np(N_nhc,i,j) = np(N_nhc,i,j) + nG(N_nhc,i,j) * delta(alpha) / real(2*n_dvd)
                        do h=N_nhc-1,1,-1
                            np(h,i,j) = np(h,i,j) * EXP( -np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd) )
                            np(h,i,j) = np(h,i,j) + nG(h,i,j) * delta(alpha) / real(2*n_dvd)
                            np(h,i,j) = np(h,i,j) * EXP( -np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd) )
                        end do
                        !!
                        do h=1,N_nhc
                            nx(h,i,j) = nx(h,i,j) + np(h,i,j) / nQ(h)* delta(alpha)/real(n_dvd)
                        end do
                        mobj%a(i)%p(j) = mobj%a(i)%p(j) * EXP( -np(1,i,j)/nQ(1)*delta(alpha)/real(n_dvd) )
                        print *,'i j k a p p1', i, j, k, alpha, mobj%a(i)%p(j), np(1,i,j)
                        !!
                        do h=1,N_nhc-1,1
                            np(h,i,j) = np(h,i,j) * EXP( -np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd) )
                            np(h,i,j) = np(h,i,j) + nG(h,i,j) * delta(alpha) / real(2*n_dvd)
                            np(h,i,j) = np(h,i,j) * EXP( -np(h+1,i,j)/nQ(h+1) * delta(alpha) / real(4*n_dvd) )
                        end do
                        np(N_nhc,i,j) = np(N_nhc,i,j) + nG(N_nhc,i,j) *delta(alpha) / real(2*n_dvd)
                        !!
                    end do
                end do
            end do
        end do
    end subroutine NHChain_partG
    
end module thermo_plus



