!!------------------------------module statbsc---------------------------------------------
module dynamic_stat
use MyDef
implicit none
type lyr
    integer :: c = 0
    integer :: lth
    real(dp), dimension(:), allocatable :: r
end type lyr
!! use pyrmd-like structure for huge data statbscs
type pyrmd
    type(lyr), dimension(:), allocatable :: l
    type(lyr), dimension(:), allocatable :: s
    integer :: c = 0
    integer :: n
    real(dp) :: aver = 0
    real(dp) :: stdv = 0
end type pyrmd
!!
contains
    subroutine creat_pyd(pyd,lst)                       ! lst: a list to discribe the structure of the pyrmd
    implicit none
        type(pyrmd), intent(inout) :: pyd
        integer, dimension(:), intent(in) :: lst
        integer :: n_lyr, i
        n_lyr = size(lst)
        pyd%n = n_lyr
        allocate( pyd%l(n_lyr) )
        allocate( pyd%s(n_lyr - 1) )
        do i=1,n_lyr-1
            allocate( pyd%l(i)%r( lst(i) ) )
            pyd%l(i)%r = 0
            allocate( pyd%s(i)%r( lst(i) ) )
            pyd%s(i)%r = 0
            pyd%l(i)%lth = lst(i)
            pyd%s(i)%lth = lst(i)
        end do
        allocate( pyd%l(n_lyr)%r( lst(n_lyr) ) )
        pyd%l(n_lyr)%lth = lst(n_lyr)
        !!
    end subroutine creat_pyd
    !!
    subroutine dat_pyrmd( pyd , dat )
    implicit none
        type(pyrmd), intent(inout) :: pyd
        real(dp), intent(in) :: dat
        integer :: i
        !!
        i=pyd%n
        do while (.true.)                               ! find the highest lyr had to be change
            if (pyd%l(i)%c < pyd%l(i)%lth ) then
                pyd%l(i)%c = pyd%l(i)%c + 1
                if (i .eq. pyd%n) then
                    pyd%l(i)%r( pyd%l(i)%c ) = dat
                else
                    pyd%s(i)%c = pyd%s(i)%c + 1
                    if (i .eq. pyd%n - 1) then
                        pyd%s(i)%r( pyd%s(i)%c ) = variance( pyd%l(i+1)%r )
                    else
                        pyd%s(i)%r( pyd%s(i)%c ) = average( pyd%s(i+1)%r )
                    end if
                    pyd%l(i)%r( pyd%l(i)%c ) = average( pyd%l(i+1)%r )                    
                end if
                exit
            end if
            i = i - 1
            if (i .eq. 0) then
                pyd%aver = ( pyd%aver * pyd%c + average( pyd%l(1)%r ) ) / ( pyd%c + 1 )
                pyd%stdv = ( pyd%stdv * pyd%c + average( pyd%s(1)%r ) ) / ( pyd%c + 1 )
                pyd%c = pyd%c + 1
                exit
            end if
        end do
        !!
        i = i + 1
        do while( i < pyd%n + 1 )
            if (i .eq. pyd%n) then
                pyd%l(i)%c = 1
                pyd%l(i)%r( pyd%l(i)%c ) = dat
                exit
            end if
            pyd%l(i)%c = 1
            pyd%l(i)%r( pyd%l(i)%c ) = average( pyd%l(i+1)%r )
            pyd%s(i)%c = 1
            if (i .eq. pyd%n - 1) then
                pyd%s(i)%r( pyd%s(i)%c ) = variance( pyd%l(i+1)%r )
            else
                pyd%s(i)%r( pyd%l(i)%c ) = average( pyd%s(i+1)%r )
            end if
            i = i + 1
        end do
        !!
        if (pyd%c .ge. 100) then
            print *, "stack warning!"
        end if
    end subroutine dat_pyrmd
    !!
    subroutine show_pyrmd(pyd)
    implicit none
        type(pyrmd), intent(in) :: pyd
        integer :: i
        print *,pyd%c
        do i=1,pyd%n
            print *,pyd%l(i)%c
        end do
    end subroutine show_pyrmd
    !!
    real(dp) function stack_aver( pyd )
    implicit none
        type(pyrmd), intent(in) :: pyd
        integer :: i
        real(dp) :: s
        real(dp) :: power
        real(dp) :: sum_power
        !!
        s = 0
        power = 1
        sum_power = 0
        !!
        ! find the highest lyr had to be change
        !!
        if (pyd%c > 0) then
            s = pyd%aver * pyd%c * power
            sum_power = sum_power + pyd%c * power
        else
            i = 1
            do while ( pyd%l(i)%c .eq. 0  .and. (i +1 < pyd%n) )
                i = i + 1
            end do
            s = s + sum( pyd%l(i)%r(1:pyd%l(i)%c) ) * power
            sum_power = sum_power + pyd%l(i)%c * power
        end if
        i = i + 1
        do while (i < pyd%n)
            power = power / size( pyd%l(i)%r )
            s = s + sum( pyd%l(i)%r(1:pyd%l(i)%c) ) * power
            sum_power = sum_power + pyd%l(i)%c * power
            i = i + 1
        end do
        !!
        stack_aver = s / sum_power        
    end function stack_aver
    !!
    real(dp) function stack_stdv( pyd )
    implicit none
        type(pyrmd), intent(in) :: pyd
        integer :: i
        real(dp) :: s
        real(dp) :: power
        real(dp) :: sum_power
        !!
        s = 0
        power = 1
        sum_power = 0
        !!
        if (pyd%c > 0) then
            s = pyd%stdv * pyd%c * power
            sum_power = sum_power + pyd%c * power
        else
            i = 1
            do while ( pyd%l(i)%c .eq. 0  .and. (i +1 < pyd%n) )
                i = i + 1
            end do
            s = s + sum( pyd%s(i)%r(1:pyd%s(i)%c) ) * power
            sum_power = sum_power + pyd%s(i)%c * power
        end if
        i = i + 1
        do while (i < pyd%n)
            power = power / size( pyd%s(i)%r )
            s = s + sum( pyd%s(i)%r(1:pyd%s(i)%c) ) * power
            sum_power = sum_power + pyd%s(i)%c * power
            i = i + 1
        end do
        !!
        stack_stdv = s / sum_power        
    end function stack_stdv
    !!
    real(dp) function average(xn)
        real(dp),dimension(:) :: xn
        real(dp) :: sumup
        integer :: i
        ! 
        sumup = 0.
        do i=1,size(xn)
            sumup = sumup + xn(i)
        end do
        average = sumup / size(xn)
    end function average
    ! 
     real(dp) function mxx(x1,x2)
        real(dp), dimension(:) :: x1, x2
        real(dp), dimension(size(x1)) :: xx
        integer :: i

        do i=1,size(x1)
            xx(i) = x1(i) * x2(i)
        end do
        mxx = average(xx)
    end function mxx
    !!
    real(dp) function variance(x1)
        real(dp), dimension(:), intent(in) :: x1
        real(dp), dimension(size(x1)) :: xx
        real(dp) :: tmpav
        integer :: i
        !!
        variance = 0.
        tmpav = average(x1)
        do i=1,size(x1)
            xx(i) = (x1(i) - tmpav) * (x1(i) - tmpav)
        end do
        variance = average(xx)
    end function variance
    !!
    real(dp) function covar(x1,x2)
        real(dp), dimension(:), intent(in) :: x1, x2
        real(dp), dimension(size(x1)) :: xx
        real(dp) :: tmp1, tmp2
        integer :: i
        !!
        do i=1,size(x1)
            xx(i) = x1(i) * x2(i)
        end do
        covar = average(xx) - average(x1) * average(x2)
    end function covar
    !!
    real(dp) function stdev(x1)
        real(dp), dimension(:), intent(in) :: x1
        !
        stdev = sqrt( variance(x1) )
    end function stdev
    !!
    real(dp) function stdevs(x1)
        real(dp), dimension(:), intent(in) :: x1
        !
        stdevs = sqrt( covar(x1, x1) * size(x1) / (size(x1)-1) )
    end function stdevs
    !
    logical function comp( ms1, ms2 ,lth )
        real(dp), dimension(2), intent(in) :: ms1, ms2
        integer :: lth
        real(dp) :: ttest, ftest
        !!
        ! size is 1000
        !!
        ttest = abs( (ms1(1) - ms2(1)) / sqrt( ( ms1(2) + ms2(2) )/ (2*lth) ) )
        if (ms1(2) > ms2(2)) then
            ftest = ms1(2) / ms2(2)
        else
            ftest = ms2(2) / ms1(2)
        end if
        !!
        print *, ms1, ms2
        print *, 'ttest', ttest, '    ', 'ftest',ftest
        comp = .false.
        if (ttest < 1.2) then
            if (ftest < 1.2) then
                comp = .true.
            end if
        end if
        !!
    end function comp
end module dynamic_stat
