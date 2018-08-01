!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 07
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

type atom
    character(len=len0) :: elem
    real(dp) :: pe, ke, kev, te, tev
    real(dp) :: m
    real(dp), dimension(:), allocatable :: x
    real(dp), dimension(:), allocatable :: p
    real(dp), dimension(:), allocatable :: ks
    !!real(dp), dimension(:), allocatable :: ps
    real(dp), dimension(:), allocatable :: fx
    real(dp), dimension(:), allocatable :: fks
end type atom

type mole
    character(len=len3) :: nm   ! necessary if use smiles formula
    integer :: nb               ! number of atom for a molecule
    real(dp) :: pe, ke, kev, te, tev
    real(dp) :: m
    real(dp) :: x, p, ks, fx, fks
    type(atom), dimension(:) , allocatable :: a
end type mole
    
end module AM_script






