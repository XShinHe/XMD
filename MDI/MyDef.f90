!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 08
!--- Acknowledgement to Liu group, of PKU
!*


!-- Gobal settings for numerical, mathematics, physics, and chemisty constants
module MyDef
implicit none

    !-- precision of float number and length of strings
    integer, parameter :: sp=kind(0.0)
    integer, parameter :: dp=kind(0.d0)
    integer, parameter :: len0=4
    integer, parameter :: len1=20
    integer, parameter :: len2=100
    integer, parameter :: len3=500
    
    !-- mathematics constant
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: twopi = 2_dp*pi
    real(dp), parameter :: pi2 = pi**2
    real(dp), parameter :: sqrtpi = sqrt(pi)
    real(dp), parameter :: e = 2.718281828459045_dp
    
    !-- physics constant
    real(dp), parameter :: hb    = 1.054571726E-19_dp
    real(dp), parameter :: kb    = 1.38064852E-23_dp
    real(dp), parameter :: au_m  = 9.10938291E-31
    real(dp), parameter :: au_e  = 1.602176565E-19
    real(dp), parameter :: au_hb = 1.054571726E-34
    real(dp), parameter :: au_ke = 8.9875517873681E+9
    real(dp), parameter :: au_c  = 137.0359991
    real(dp), parameter :: au_a0 = 5.2917721092E-11
    real(dp), parameter :: au_eh = 4.35974417E-18
    real(dp), parameter :: au_t  = 2.418884326505E-17
    real(dp), parameter :: au_temp = 3.1577464E+5
    real(dp), parameter :: au_kb = 1.00000000000
    real(dp), parameter :: au_beta = 3.166815367E-6
    
    !-- IO name settings
    integer, parameter  :: out_unit =10, ana_unit =11
    character(len=len2) :: out_file, ana_file
    character(len=len2) :: putrc_file  = 'put.rc'
    character(len=len2) :: cnfgrc_file = 'cnfg.rc'
    character(len=len2) :: rst_file = 'xmd.rst'
    character(len=len2) :: smp_file = 'xmd.smp'
    character(len=len2) :: lst_file = 'xmd.lst'

    !-- sharing varibles
    integer :: my_iostat ! indicates status of reading from a file
    logical :: my_exist  ! indicates existence of a given file
    real(dp) :: rand_throw ! share a memory to generate a random number

contains
	!-- Mathematics procedures
	!-- for initial random seed from corrent date and time
	subroutine init_seed()
		integer :: n, ival(8), v(3), i
		integer, allocatable :: seed(:)
		call date_and_time(values=ival)
		v(1) = ival(8) + 2048*ival(7)
		v(2) = ival(6) + 64*ival(5)     ! value(4) isn't really 'random'
		v(3) = ival(3) + 32*ival(2) + 32*8*ival(1)
	 	call random_seed(size=n)
		allocate(seed(n))
		call random_seed()              ! give the seed an implementation-dependent kick
		call random_seed(get=seed)
		do i=1, n
			seed(i) = seed(i) + v(mod(i-1, 3) + 1)
	  	enddo
	  	call random_seed(put=seed)
	  	deallocate(seed)
	end subroutine init_seed

	!-- random_norm: generate a normal distribution random number
	!------ we can use many methods, such as
	!------ 1) Box-Muller algorithm
	!------ 2) Marsaglia algorithm
	!------ 3) Zigguart algorithm (a kind of rejection method) 
	!------ * note that here we just use the most easy one, Box-Muller algorithm
	subroutine random_norm(n_norm, n_u, n_sigma)
	implicit none
		real(dp), intent(inout) :: n_norm, n_u
		real(dp), intent(in) :: n_sigma
		real(dp) :: radius

		call random_number(n_u)
		radius = sqrt(-2*log(n_u))*n_sigma
		call random_number(n_u)
		!-- the defualt occassion is use cos() function, alternatively can be substitute by sin()
		n_norm = radius*dcos(2*pi*n_u)
	end subroutine random_norm

	!-- messager procedures
	subroutine send_err(msg)
	implicit none
		character(*), intent(in) :: msg
		print *, msg
		stop
	end subroutine send_err

end module MyDef



