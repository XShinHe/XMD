! the my defined force field file
module myFF
use MyDef
implicit none
    type ffbind
        integer :: ffnumber
    end type ffbind
    
contains
    subroutine show_myFF()
        print *,'the my defined force fields are all list here:&
            &1) Harmonic Osscillitor:    Vpn2,   Fpn2 (stand for polynomial)&
            &2) Quardic Potenial:        Vpn4,   Fpn4 &
            &3) General Polynomial:      Vpng,   Fpng &
            &4) Double Well Potential:   Vdbw,   Fdbw &
            &5) Morse Potential:         Vmor,   Fmor &
            &6) Van de Waals (Lennard-Jones) Potential:&
            &                            Vvdw,   Fvdw &
            &7) Others'
    end subroutine show_myFF
    
    ! Harmonic Oscillitor potential and force
    real(dp) function Vpn2(k,x)
        real(dp), intent(in) :: k, x
        Vpn2 = k*x*x
    end function Vpn2
    
    real(dp) function Fpn2(k,x)
        real(dp), intent(in) :: k, x
        Fpn2 = 2_dp*k*x
    end function Fpn2
    
    ! Quardic potential and force 
    real(dp) function Vpn4(k,x)
        real(dp), intent(in) :: k, x
        Vpn4 = k*x*x*x*x
    end function Vpn4
    
    real(dp) function Fpn4(k,x)
        real(dp), intent(in) :: k, x
        Fpn4 = 4_dp*k*x*x*x
    end function Fpn4
    
    ! General polynomial potential and force, neglect the constant term for such term can be always set to zero
    real(dp) function Vpng(k,x)
        real(dp), intent(in), dimension(:) :: k
        real(dp), intent(in) :: x
        real(dp) :: term
        integer :: i
        Vpng=0
        term=1
        do i=1,size(k)
            term=term*x
            Vpng=Vpng+term*k(i)
        end do
    end function Vpng
    
    real(dp) function Fpng(k,x)
        real(dp), intent(in), dimension(:) :: k
        real(dp), intent(in) :: x
        real(dp) :: term
        integer :: i
        Fpng=0
        term=1
        do i=1,size(k)
            Fpng=Fpng + real(i,kind=dp)*term*k(i)
            term=term*x
        end do
    end function Fpng
    
    ! the doubel well potential and force
    real(dp) function Vdbw(k, a, x)
        real(dp), intent(in) :: k, a, x
        Vdbw = k*(1-a*x*x)*(1-a*x*x)
    end function Vdbw
    
    real(dp) function Fdbw(k, a, x)
        real(dp), intent(in) :: k, a, x
        Fdbw = -4_dp*k*(1-a*x*x)*a*x
    end function Fdbw
    
    ! Poschl-Teller potential, without force
    real(dp) function Vpt0(k, a, x)
        real(dp), intent(in) :: k, a, x
        Vpt0 = -k*(k+1_dp)/(2_dp*DCOSH(x/a)*DCOSH(x/a))
    end function Vpt0
    
    ! Morse Potential
    real(dp) function Vmor(D, a, r0, x)
        real(dp), intent(in) :: D, a, r0, x
        Vmor = D*(1-DEXP(-a*(x-r0)))*(1-DEXP(-a*(x-r0)))
    end function Vmor
    
    real(dp) function Fmor(D, a, r0, x)
        real(dp), intent(in) :: D, a, r0, x
        Fmor = 2_dp*D*(1-DEXP(-a*(x-r0)))*DEXP(-a*(x-r0))*a
    end function Fmor
    
    ! Lennard Jones potential and force
    real(dp) function Vvdw(e, s, x)
        real(dp), intent(in) :: e, s, x
        Vvdw = 4_dp*e*((s/x)**12-(s/x)**6)
    end function Vvdw
    
    real(dp) function Fvdw(e, s, x)
        real(dp), intent(in) :: e, s, x
        Fvdw = -24_dp*e*(2*(s/x)**13-(s/x)**7)
    end function Fvdw
    
    subroutine F3LJ001( f3, xa3, xb3 )
        real(dp), dimension(3), intent(in) :: xa3, xb3
        real(dp), dimension(3), intent(inout) :: f3
        real(dp), dimension(3) :: r3
        real(dp) :: radius
        r3 = xa3 - xb3
        radius = DSQRT( dot_product( r3, r3 ) )
        f3 = r3 / radius * Fvdw(1.0_dp,1.0_dp, radius )
    end subroutine F3LJ001
    
    subroutine V3LJ001( v, xa3, xb3 )
        real(dp), dimension(3), intent(in) :: xa3, xb3
        real(dp), intent(inout) :: v
        real(dp), dimension(3) :: r3
        real(dp) :: radius
        r3 = xa3 - xb3
        radius = DSQRT( dot_product( r3, r3 ) )
        v = Vvdw(1.0_dp,1.0_dp, radius )
    end subroutine V3LJ001
    
    subroutine V3PN002( v, xa3, a_par )
        real(dp), dimension(3), intent(in) :: xa3
        real(dp), intent(inout) :: v
        real(dp), intent(in) :: a_par
        real(dp) :: radius
        radius = DSQRT( dot_product( xa3, xa3 ) )
        v = Vpn2(a_par, radius )
    end subroutine V3PN002
    
end module myFF




