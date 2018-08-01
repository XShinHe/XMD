!*
! ---Copyright by--- He Xin <1500011805@pku.edu.cn>
!*

module AM_script
use MyDef
!use fsmiles
implicit none
    integer,parameter,private		:: nelement_max = 54
	character(len=len0),parameter,private :: element_list(nelement_max) =                       &
  	(/   ' H',                                                                                'He', &  !  2
         'Li','Be',                                                  ' B',' C',' N',' O',' F','Ne', &  ! 10
         'Na','Mg',                                                  'Al','Si',' P',' S','Cl','Ar', &  ! 18
         ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &  ! 36
         'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe'  /) ! 54

type atom
    character(len=len0) :: elem
    !!
    real(dp) :: pe, ke, kev, te, tev
    real(dp) :: m
    real(dp), dimension(:), allocatable :: x
    real(dp), dimension(:), allocatable :: p
    real(dp), dimension(:), allocatable :: ks
    !real(dp), dimension(:), allocatable :: ps
    real(dp), dimension(:), allocatable :: fx
    real(dp), dimension(:), allocatable :: fks
end type atom
!!
type mole
    character(len=len3) :: nm ! necessary use smiles
    integer :: nb ! number of atom
    type(atom), dimension(:) , allocatable :: a 
end type mole

contains
    subroutine build_mole(mobj, anumber)
    implicit none
        type(mole), intent(inout) :: mobj
        integer, intent(in) :: anumber
        character(len=len3) :: ins_s
        character(len=len0) :: ins_c
        real(dp) :: ins_r
        integer :: ins_i
        integer :: i,j
        
        cnfg_file = 'cnfg.rc'
        inquire(file=cnfg_file, exist=my_exist)
        if (my_exist .eqv. .false.) then
            call send_err("error, configuration file (cnfg.rc) opening wrong")
        end if
        open(unit=10,file=cnfg_file)
        read(10,*) ins_s, ins_i
        !!
        !if (.not. ins_s(1,1) .eq.'!') then
        !    call send_err('error: format error of _cnfg.in')
        !end if
        mobj%nm=ins_s(2:len(ins_s))
        mobj%nb=ins_i
        allocate(mobj%a( mobj%nb ))
        do i=1,mobj%nb
            allocate( mobj%a(i)%x( anumber ) )
            allocate( mobj%a(i)%ks( anumber ) )
            allocate( mobj%a(i)%p( anumber ) )
            allocate( mobj%a(i)%fx( anumber ) )
            allocate( mobj%a(i)%fks( anumber ) )
        end do
        do i=1,anumber
            read(10,*,iostat=my_iostat) ins_c, ins_r
            if(my_iostat<0) exit
            mobj%a(i)%elem = ins_c
            mobj%a(i)%m = ins_r
        end do
        close(unit=10)
    end subroutine build_mole
end module AM_script






