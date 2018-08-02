!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 08
!--- Acknowledgement to Liu group, of PKU
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

type, public :: atom
    character(len=len0) :: iname
    real(dp) :: pe, ke, te
    real(dp) :: m
    real(dp), dimension(:), allocatable :: ks
    real(dp), dimension(:), allocatable :: x
    real(dp), dimension(:), allocatable :: p
    real(dp), dimension(:), allocatable :: fx
    real(dp), dimension(:), allocatable :: fks
end type atom

type, public, extends(atom) :: mole
    character(len=len3) :: formula      ! specified for smiles-formula
    type(atom), dimension(:) , allocatable :: a
end type mole

end module AM_script




