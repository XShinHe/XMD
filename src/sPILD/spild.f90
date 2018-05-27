!! an almost simplest version of pimd
!!
! ------------------------- basic module -------------------------
!!
module md_spild
implicit none
    real(kind=8), parameter :: pi = 3.14159265358979323846
    real(kind=8), parameter :: e = 2.718281828459045
    real(kind=8), parameter :: hb = 1.054571726E-19
    real(kind=8), parameter :: kb = 1.38064852E-23
    ! atom unit information
    real(kind=8), parameter :: au_m = 9.10938291E-31
    real(kind=8), parameter :: au_e = 1.602176565E-19
    real(kind=8), parameter :: au_hb = 1.054571726E-34
    real(kind=8), parameter :: au_ke = 8.9875517873681E+9
    real(kind=8), parameter :: au_c = 137.0359991
    real(kind=8), parameter :: au_a0 = 5.2917721092E-11
    real(kind=8), parameter :: au_eh = 4.35974417E-18
    real(kind=8), parameter :: au_t = 2.418884326505E-17
    real(kind=8), parameter :: au_temp = 3.1577464E+5
    real(kind=8), parameter :: au_kb = 1.00000000000
    real(kind=8), parameter :: au_beta = 3.166815367E-6
    !!
    integer :: md_nstep
    integer :: md_npass
    integer :: md_is_eqb = 0
    integer :: md_scheme
    integer :: md_thermo
    integer :: md_virtual
    integer :: md_nsmp = 1000
    !!
    integer :: md_bead
    real :: md_bfreq
    real :: md_bf2
    !!
    integer :: md_opt = 0
    integer :: md_distr_x = 0
    integer :: md_distr_p = 1
    !!
    real(kind=8) :: md_temp
    real(kind=8) :: md_beta
    real(kind=8) :: md_mass
    real(kind=8) :: md_dtime
    real(kind=8) :: md_cfreq
    real(kind=8) :: md_coeff
    !!
    real(kind=8) :: md_rand
    !!
    real(kind=8) :: md_pe
    real(kind=8) :: md_ke
    real(kind=8) :: md_kev
    real(kind=8) :: md_tote
    real(kind=8), dimension(:,:), allocatable :: tmplst
    !!
    integer, parameter :: out_unit = 10, ana_unit=11
    integer :: md_erst = 0
    integer :: md_irst = 0
! a struct for beadectory
    !!
    real(kind=8), dimension(:), allocatable :: bd_x
    real(kind=8), dimension(:), allocatable :: bd_ks
    real(kind=8), dimension(:), allocatable :: bd_p
    real(kind=8), dimension(:), allocatable :: bd_fx
    real(kind=8), dimension(:), allocatable :: bd_fks
    real(kind=8), dimension(:), allocatable :: bd_m
! for pild process
    real(kind=8) :: md_mth
    real(kind=8) :: md_gammaAD=0.0001_8 ! here gammaAD is bound to be small
                                        ! to seperate the motion of the first bead with others
    real(kind=8) :: mda_x2x2
    real(kind=8) :: mda_pp
!!
contains
!!
! %fun% pential energy surface
!!
real(kind=8) function fun_pes(x)
implicit none
    real(kind=8), intent(in) :: x
    fun_pes = 0.500000 * ( x * x )
end function fun_pes
!!
subroutine calc_ana()
implicit none
    integer :: i
    mda_x2x2 = bd_x(1)**2 + 0.25*md_beta/md_mth - 0.25*(bd_p(1)*md_beta/md_mth)**2
    mda_pp = bd_m(1)/md_mth*bd_p(1)
end subroutine calc_ana
!!
subroutine calc_ks()
implicit none
    integer :: i
    bd_ks(1) = bd_x(1)
    !if(md_bead.eq.1) return
    bd_ks(md_bead) = bd_x(md_bead) - bd_x(1)
    do i=2,md_bead-1
        bd_ks(i) = bd_x(i) - ( bd_x(i+1) * real(i - 1) + bd_x(1) ) / real(i)
    end do
end subroutine calc_ks
!!
subroutine calc_x() ! attention bug that md_bead should be 1, for repeated!
implicit none
    integer :: i
    bd_x(1) = bd_ks(1)
    !if(md_bead.eq.1) return
    bd_x(md_bead) = bd_ks(md_bead) + bd_ks(1)
    do i=md_bead-1,2,-1
        bd_x(i) = bd_ks(i) + ( real( i - 1 ) * bd_x(i+1) + bd_ks(1) ) / real(i)
    end do
end subroutine calc_x
!!
subroutine calc_fx()
implicit none
    integer :: i
    do i=1,md_bead
        bd_fx(i) = 2*bd_x(i) - 0.3*bd_x(i) + 0.4*bd_x(i)**3
    end do
end subroutine calc_fx
!!
subroutine calc_fks()
implicit none
    integer :: i
    !!
    call calc_fx()
    bd_fks(1) = 0.0
    do i=1,md_bead
        bd_fks(1) = bd_fks(1) + bd_fx(i)
    end do
    bd_fks(1) = bd_fks(1)
    !!
    do i=2,md_bead,1
        bd_fks(i) = bd_fx(i) + bd_fks(i-1) * real(i-2)/real(i-1)
    end do
    !!
    bd_fks(1) = bd_fks(1)/real(md_bead)
    do i=2,md_bead,1
        bd_fks(i) = bd_fks(i)/real(md_bead) + bd_m(i) * md_bf2 * bd_ks(i)
    end do
    !!
end subroutine calc_fks
!
! %sub% a substitute initial random_seed
!
subroutine init_seed()
	integer :: n, ival(8), v(3), i
	integer, allocatable :: seed(:)
	call date_and_time(values=ival)
	v(1) = ival(8) + 2048*ival(7)
	v(2) = ival(6) + 64*ival(5)     ! value(4) isn't real(kind=8)ly 'random'
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
!
! %sub% a generator of a gaussian random number
! n_u1, n_u2 need to be (different) unifrom distributions
!
subroutine random_norm(n_norm, n_u, n_sigma, flag)
implicit none
    real(kind=8), intent(inout) :: n_norm, n_u
    real(kind=8), intent(in) :: n_sigma
    character, optional :: flag
    real(kind=8) :: tmp_c

    call random_number(n_u)
    tmp_c = sqrt(-2*log(n_u))*n_sigma
    call random_number(n_u)
    !!
    if (present(flag) .and. flag == 'c') then
        n_norm = tmp_c*cos(2*pi*n_u)
    else
        n_norm = tmp_c*sin(2*pi*n_u)
    end if
    !!
end subroutine random_norm
!!
!!
integer function init_md()
implicit none
    real(kind=8) :: tmp1, tmp2, tmp_ppT
    character(len=8) :: tmpin1, tmpin2
    real(kind=8) :: tmpin3
    integer :: istat, n_init
    logical :: alive
    !!
    n_init = 0
    inquire(file='_put.in', exist=alive)
    if (alive .eqv. .true.) then
        open(unit=10, file='./_put.in')
        do while (.true.)
            read(10,*,iostat=istat) tmpin1, tmpin2
            if (istat < 0) exit
            select case (tmpin1)
                case ('temp')
                    read(tmpin2,*) md_temp
                    md_beta = 1.0_8 / md_temp
                case ('mass')
                    read(tmpin2,*) md_mass
                case ('sfreq')
                    read(tmpin2,*) md_cfreq
                case ('dtime')
                    read(tmpin2,*) md_dtime
                case ('coeff')
                    read(tmpin2,*) md_coeff
                case ('nstep')
                    read(tmpin2,*) md_nstep
                case ('scheme')
                    read(tmpin2,*) md_scheme
                case ('thermo')
                    read(tmpin2,*) md_thermo
                case ('mirror')
                    read(tmpin2,*) md_virtual
                case ('nsmp')
                    read(tmpin2,*) md_nsmp
                    n_init = n_init - 1
                case ('gammaAD')
                    read(tmpin2,*) md_gammaAD
                    n_init = n_init - 1
                case ('bead')
                    read(tmpin2,*) md_bead
                    n_init = n_init - 1
                case ('opt')
                    read(tmpin2,*) md_opt
                    n_init = n_init - 1
                case ('distrx')
                    read(tmpin2,*) md_distr_x
                    n_init = n_init - 1
                case ('distrp')
                    read(tmpin2,*) md_distr_p
                    n_init = n_init - 1
                case default
                    n_init = n_init - 1
            end select
            n_init = n_init + 1
        end do
    else
        print *, 'file ./_put.in is not existing, please check it'
    end if
    !!
    ! input the evaluated the momenta info
    open(unit=10,file='.mth.tmp',status='old')
        read(10,*) tmp_ppT
    close(unit=10)
    md_mth=tmp_ppT*md_beta
    print *, 'mthermo', md_mth
    !!
    if (md_opt .eq. 1) then
        md_dtime = 0.2 / md_cfreq
        md_nstep = 10000 / md_cfreq
        select case (md_thermo)
            case (0)
                md_coeff = md_cfreq
            case (1)
                md_coeff = sqrt(2.) * md_cfreq
            case (2)
                md_coeff = 20 * md_dtime
        end select
    end if
    !!
    md_npass = 0
    md_is_eqb = 0
    md_bfreq = sqrt( real(md_bead) ) * md_temp
    md_bf2 = real(md_bead) * md_temp * md_temp
    !!
    if (n_init < 9) then
        init_md = 0
    else
        init_md = 1
    end if
end function init_md

end module md_spild

! ------------------------- main -------------------------
!
program main
use md_spild
implicit none
    interface
        subroutine init_tj()
        use md_spild
        implicit none
            integer :: i
        end subroutine init_tj
        !!
        subroutine md_control( savetofile, ana_file )
        use md_spild
        implicit none
            character(len=20), intent(in) :: savetofile, ana_file
            integer :: i
        end subroutine md_control
    end interface

    character(len=20) :: traj_outfile, ana_file                       ! filename to discribe the trajectory
    character(len=20) :: para_infile                        ! filename input the parameters in as 2nd args
    character(len=8) :: flg,arg                                 ! get n-th trajectory as 1st args
    integer :: n, i, j, init_md_completed
    !!
    !!
    n = command_argument_count()
    if ((n .ne. 2)) then
        print *,"error, arguments mismatch"
        stop
    end if
    !! if nessessary check for some arguments
    call get_command_argument(1,flg)
    call get_command_argument(2,arg)
    select case (flg)
        case ("-x")
            select case (arg)
                case ('e')
                    md_erst=1
                    traj_outfile = trim(arg)//'.xpe'
                    ana_file = trim(arg)//'.ena'
                case ('r')
                    md_erst=2
                    traj_outfile = trim(arg)//'.xpr'
                    ana_file = trim(arg)//'.ana'
                case default
                    md_erst=0
                    traj_outfile = trim(arg)//'.xpf'
                    ana_file = trim(arg)//'.ana'
            end select
        case ("-f")
            md_erst=3
            read(arg,*) md_irst
            traj_outfile = trim(arg)//'.xpf'
            ana_file = trim(arg)//'.ana'
        case default
            print *,"error, arguments mismarch"
    end select
    
    !!
    init_md_completed = init_md()
    if ( init_md_completed .eq. 0 ) then
        stop
    end if
    !!
    ! initial and then run
    call init_seed()
    call init_tj()
    call md_control(traj_outfile, ana_file)
!!
end program main


! --------------------------(init) subr -------------------------
! (init_tj )
!!
! initialization and optimization control of (x,p)
!!
subroutine init_tj()
use md_spild
implicit none
    ! do some parameters optimization
    ! 
    integer :: i
    real(kind=8), dimension(7) :: ins
    allocate(bd_x(md_bead))
    allocate(bd_ks(md_bead))
    allocate(bd_p(md_bead))
    allocate(bd_fx(md_bead))
    allocate(bd_fks(md_bead))
    allocate(bd_m(md_bead))
    allocate(tmplst(10,md_bead))
    !!
    do i=1,md_bead
        !!
        if (i .eq. 1) then
            bd_m(i) = md_mass
        else
            bd_m(i) = md_mass * real(i) / real( i - 1 )
        end if
        !!
        tmplst(1,i) = sqrt( bd_m(i) * md_temp )                             ! variance of p
        tmplst(2,i) = tmplst(1,i) * md_mass * md_cfreq                      ! variance of x
        tmplst(3,i) = exp( - md_coeff * md_dtime )                          ! e^(-g*dt)
        tmplst(4,i) = tmplst(1,i) * sqrt( 1 - tmplst(3,i) * tmplst(3,i) )
        tmplst(5,i) = exp( - 0.5 * md_coeff * md_dtime )                    ! e^(-0.5*g*dt)
        tmplst(6,i) = tmplst(1,i) * sqrt( 1 - tmplst(5,i) * tmplst(5,i) )
        !!
    end do
    ! init x(or xi) and p
    select case (md_erst)
        case(2)
            open(unit=10,file='end.rst')
            do i=1,md_bead
                read(10,*) ins
                bd_ks(i) = ins(3)
                bd_p(i) = ins(4)
                bd_x(i) = ins(5)
            end do
            close(unit=10)
        case(3)
            open(unit=10,file='e.xpe')
            read(10,*)
            do i=1,md_irst
                read(10,*)
            end do
            do i=1,md_bead
                read(10,*) ins
                bd_ks(i) = ins(3)
                bd_p(i) = ins(4)
                bd_x(i) = ins(5)
            end do
            close(unit=10)
        case default
            do i=1,md_bead
                if (md_distr_x .eq. 1) then
                    call random_norm(bd_ks(i), md_rand, tmplst(2,i))
                    if (md_distr_p .eq. 0) then
                        bd_p(i) = 0.0
                    else
                        call random_norm(bd_p(i), md_rand, tmplst(1,i))
                    end if
                else
                    bd_ks(i) = 0.0
                    call random_norm(bd_p(i), md_rand, tmplst(1,i))
                end if
            end do
    end select
    call calc_x()
    call calc_fks()
    call calc_ana()
    !!
end subroutine init_tj


! ------------------------- (run) subr -------------------------
!
subroutine md_control(savetofile, ana_file)
use md_spild
implicit none
    interface
        subroutine md_run_t(n_dvd)
        use md_spild
        implicit none
            integer, intent(in) :: n_dvd
            integer :: i
            real(kind=8) :: et
        end subroutine md_run_t
        !!
        subroutine md_run_x(n_dvd)
        use md_spild
        implicit none
            integer, intent(in) :: n_dvd
            integer :: i
        end subroutine md_run_x
        !!
        subroutine md_run_p(n_dvd)
        use md_spild
        implicit none
            integer, intent(in) :: n_dvd
            integer :: i
        end subroutine md_run_p
        !!
        subroutine md_sample()
        use md_spild
        implicit none
            integer :: i
        end subroutine
        !!
    end interface

    character(len=20), intent(in) :: savetofile, ana_file
    integer :: i
    !!
    open(unit=out_unit, file=savetofile, status='replace')
    write(out_unit,*) 'nstep    ','nbead    ','ks   ','p    ','x    ','fks  ', 'fx   '
    open(unit=ana_unit, file=ana_file, status='replace')
    write(ana_unit,*) 'md_npass     ', 'bd_x1   ', 'bd_p1   ', 'mda_x2x2    ', 'mda_pp  ', 'md_pe   '
    !!
    select case (md_scheme)
        case (0) ! middle scheme
            call md_sample()
            !!
            do i  = 1,md_nstep
                !!
                call md_run_p(2)
                call md_run_x(2)
                call md_run_t(1)
                call md_run_x(2)
                call md_run_p(2)
                !!
                md_npass = md_npass + 1
                if(mod(i,md_nsmp) .eq. 0) then
                    call md_sample()
                end if
                !!
            end do
        case (1) ! end scheme
            call md_sample()
            !!
            do i  = 1,md_nstep
                !!
                call md_run_p(2)
                call md_run_x(1)
                call md_run_p(2)
                call md_run_t(1)
                !!
                md_npass = md_npass + 1
                call md_sample()
                !!
            end do
        case (2) ! side scheme
            call md_sample()
            !!
            do i  = 1,md_nstep
                !!
                call md_run_t(2)
                call md_run_p(2)
                call md_run_x(1)
                call md_run_p(2)
                call md_run_t(2)
                !!
                md_npass = md_npass + 1
                call md_sample()
                !!
            end do
    end select
    close(unit=out_unit)
    close(unit=ana_unit)
    !!
    open(unit=10,file='end.rst',status='replace')
    do i=1,md_bead
        write(10,*) md_npass, i, bd_ks(i), bd_p(i), bd_x(i), bd_fks(i), bd_fx(i)
    end do
    close(unit=10)
end subroutine md_control

! shouldn't always open and close
subroutine md_sample()
use md_spild
implicit none
    integer :: i,j,n
    !!
    if (md_erst .ne. 1 .or. md_npass > md_nstep/2 ) then
        write (ana_unit,*) md_npass, bd_x(1), bd_p(1), mda_x2x2, mda_pp, md_pe
    end if
end subroutine md_sample


! an update of thermostat
! flag = 0 'lang' : langevin thermostat
! flag = 1 'ads' : andersen thermostat
! flag = 2 'nhc' : nose-hoover chain
!
subroutine md_run_t(n_dvd)
use md_spild
implicit none
    integer, intent(in) :: n_dvd
    integer :: i
    real(kind=8) :: et
    !!
    !!
    do i=2,md_bead ! for PILD, the first bead is free of thermostat
        select case (md_thermo)
            case (0)
                call random_norm(et, md_rand, 1.0_8)
                bd_p(i) = tmplst( 2*n_dvd + 1 , i ) * bd_p(i) + sqrt(md_gammaAD) * tmplst( 2*n_dvd + 2, i ) * et
                ! note md_gammaAD here!
            case (1)
                call random_number(md_rand)
                if ( md_rand < 1.- tmplst( 2*n_dvd + 1 , i ) ) then
                    call random_norm(et, md_rand, 1.0_8)
                    bd_p(i) = et * tmplst(1,i) * sqrt(md_gammaAD)
                end if
            case (2)
                print *,'is need to promove'
        end select
    end do
end subroutine md_run_t


! an update of position
!
subroutine md_run_x(n_dvd)
use md_spild
implicit none
    integer, intent(in) :: n_dvd
    integer :: i
    !!
    bd_ks(1) = bd_ks(1) + bd_p(1) * md_dtime / ( real(n_dvd) * bd_m(1) )
    do i=2,md_bead
        bd_ks(i) = bd_ks(i) + bd_p(i) * md_dtime / ( real(n_dvd) * bd_m(i) * md_gammaAD)
    end do
    call calc_x()
    call calc_fks()
    call calc_ana()
!!
end subroutine md_run_x

! an update of momentum
subroutine md_run_p(n_dvd)
use md_spild
implicit none
    integer, intent(in) :: n_dvd
    integer :: i
    !!
    bd_p(1) = bd_p(1) - md_mth/bd_m(1)*bd_fks(1) * md_dtime / real(n_dvd)
    do i=2,md_bead
        bd_p(i) = bd_p(i) - bd_fks(i) * md_dtime / real(n_dvd)
    end do
!!
end subroutine md_run_p
