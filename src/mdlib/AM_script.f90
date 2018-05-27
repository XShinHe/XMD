!*
! ---Copyright by--- XShinHe <hx0824916@163.com/1500011805@pku.edu.cn>
!*

! use atomic and molecular script
module AM_script
use MyDef
!use fsmiles
implicit none
    integer,parameter,private :: nelement = 54
    character(len=len0),parameter,private :: element_list(nelement) =                       &
    (/   ' H',                                                                                'He', &  !  2
         'Li','Be',                                                  ' B',' C',' N',' O',' F','Ne', &  ! 10
         'Na','Mg',                                                  'Al','Si',' P',' S','Cl','Ar', &  ! 18
         ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &  ! 36
         'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe'  /) ! 54
         
type atom
    character(len=len0) :: elem
    integer :: bd ! bead for atom, for different atom can use different scale, Advanced feature
    real(dp) :: pe, ke, kev, te, tev
    real(dp) :: m
    real(dp), dimension(:), allocatable :: x
    real(dp), dimension(:), allocatable :: p
    real(dp), dimension(:), allocatable :: ks
    real(dp), dimension(:), allocatable :: fx
    real(dp), dimension(:), allocatable :: fks
end type atom
!!
type mole
    character(len=len3) :: nm ! necessary use fsmiles.f90
    integer :: nb ! number of atom
    integer :: dm ! dimension, 1, 2, 3
    real(dp) :: pes
    type(atom), dimension(:) , allocatable :: a
end type mole

contains
    subroutine build_mole(mobj, bdnumber)
    implicit none
        type(mole), intent(inout) :: mobj
        integer, intent(in) :: bdnumber
        character(len=len3) :: input_string
        character(len=len0) :: input_chars
        real(dp) :: input_real
        integer :: input_int
        integer :: i
        
        inquire(file=cnfg_file, exist=my_exist)
        if (my_exist .eqv. .false.) then
            call send_err("error, configuration file (_cnfg.in) opening wrong")
        end if
        open(unit=10,file=cnfg_file)
        read(10,*) input_string, input_int
        !!
        read(input_string(2:2),*) mobj%dm
        mobj%nm=input_string(4:len(input_string))
        mobj%nb=input_int
        allocate(mobj%a( mobj%nb ))
        do i=1,mobj%nb
            allocate( mobj%a(i)%x( 3*bdnumber ) )
            allocate( mobj%a(i)%ks( 3*bdnumber ) )
            allocate( mobj%a(i)%p( 3*bdnumber ) )
            allocate( mobj%a(i)%fx( 3*bdnumber ) )
            allocate( mobj%a(i)%fks( 3*bdnumber ) )
        end do
        do i=1,mobj%nb
            read(10,*,iostat=my_iostat) input_chars, input_real, input_int
            if(my_iostat<0) exit
            mobj%a(i)%elem = input_chars
            mobj%a(i)%m = input_real
            if(input_int > 0) then
                mobj%a(i)%bd = input_int
            else
                mobj%a(i)%bd = bdnumber
            endif
        end do
        close(unit=10)
    end subroutine build_mole
end module AM_script
