!*
! ---Copyright by--- He Xin <1500011805@pku.edu.cn>
!*

module MyDef ! mathematics, physics, chemistry basic
implicit none
    ! control the precision
    integer, parameter :: sp=kind(0.0)
    integer, parameter :: dp=kind(0.d0)
    integer, parameter :: len0=2
    integer, parameter :: len1=10
    integer, parameter :: len2=50
    integer, parameter :: len3=100
    ! mathematics constant
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: twopi = 2_dp*pi
    real(dp), parameter :: pi2 = pi**2
    real(dp), parameter :: sqrtpi = sqrt(pi)
    real(dp), parameter :: e = 2.718281828459045_dp
    ! physics constant
    real(dp), parameter :: hb = 1.054571726E-19_dp
    real(dp), parameter :: kb = 1.38064852E-23_dp
    real(dp), parameter :: au_m = 9.10938291E-31
    real(dp), parameter :: au_e = 1.602176565E-19
    real(dp), parameter :: au_hb = 1.054571726E-34
    real(dp), parameter :: au_ke = 8.9875517873681E+9
    real(dp), parameter :: au_c = 137.0359991
    real(dp), parameter :: au_a0 = 5.2917721092E-11
    real(dp), parameter :: au_eh = 4.35974417E-18
    real(dp), parameter :: au_t = 2.418884326505E-17
    real(dp), parameter :: au_temp = 3.1577464E+5
    real(dp), parameter :: au_kb = 1.00000000000
    real(dp), parameter :: au_beta = 3.166815367E-6
    !!
    ! procedure control
    
    integer, dimension(10) :: exe_vi
    real(dp), dimension(10) :: exe_vr
    character(len=len1), dimension(10,2) :: exe_args
    character(len=len1), dimension(2) :: exe_arg2
    !!
    integer, parameter :: trj_unit =13, erg_unit =14, trs_unit=23, els_unit=24, ana_unit=25
    character(len=len2) :: traj_file
    character(len=len2) :: enrg_file
    character(len=len2) :: trsf_file
    character(len=len2) :: else_file
    character(len=len2) :: anal_file
    character(len=len2) :: mdin_file='_put.in'
    character(len=len2) :: cnfg_file='_cnfg.in'
    !!
    ! shared varibles
    integer :: my_iostat
    logical :: my_exist
    logical :: my_tf_1, my_tf_2, my_tf_3
    real(dp) :: rand_throw
!!
contains
!! 1) especially for some mathematics procedure!

subroutine init_seed()
    integer :: n, ival(8), v(3), i
    integer, allocatable :: seed(:)
    call date_and_time(values=ival)
    v(1) = ival(8) + 2048*ival(7)
    v(2) = ival(6) + 64*ival(5)     ! value(4) isn't real(dp)ly 'random'
    v(3) = ival(3) + 32*ival(2) + 32*8*ival(1)
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed()   ! Give the seed an implementation-dependent kick
    call random_seed(get=seed)
    do i=1, n
        seed(i) = seed(i) + v(mod(i-1, 3) + 1)
    enddo
    call random_seed(put=seed)
    deallocate(seed)
end subroutine init_seed

subroutine random_norm(n_norm, n_u, n_sigma)
implicit none
    real(dp), intent(inout) :: n_norm, n_u
    real(dp), intent(in) :: n_sigma
    real(dp) :: tmp_c

    call random_number(n_u)
    tmp_c = sqrt(-2*log(n_u))*n_sigma
    call random_number(n_u)
    n_norm = tmp_c*dcos(2*pi*n_u)       ! the defualt occassion is use cos() function, alternatively can be substitute by sin()
    !!
end subroutine random_norm

!! 2) the manager procedure general
subroutine send_err(msg)
implicit none
    character(*), intent(in) :: msg
    print *, msg
    stop
end subroutine send_err

end module MyDef



