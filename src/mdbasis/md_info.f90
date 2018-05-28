module md_info
use MyDef
use AM_script
implicit none
    !! esemble control
    real(dp) :: md_temp
    real(dp) :: md_beta
    real(dp) :: md_pres
    real(dp) :: md_volm
    !!
    !! dynamic choose
    integer :: md_scheme
    integer :: md_thermo
    integer :: md_virtual
    integer :: md_methd ! primitive, staging, normal-mode, 012
    !!
    !! outer changalbe for dynamics control
    real(dp) :: md_coeff
    !!
    !! time/step control
    integer :: md_nstep
    integer :: md_npass
    real(dp) :: md_dtime
    !!
    !! bead set control
    integer :: md_bead = 1              ! default vaule
    integer :: md3_bead
    real(dp) :: md_bfreq                 ! default
    real(dp) :: md_bfreq2
    !!
    !! IO related control
    integer :: md_opt = 0
    integer :: md_mod=0 ! 0,run; 1, eqb; 2, restart; 3, from file to import
    integer :: md_mod3n=0
    integer :: md_iseqb = 0
    integer :: md_nsmp=1000
    !!
    !! system instinct parameter  
    real(dp) :: md_sysfreq                    ! mass and atomfreq are should be arranged to _cnfg.in file, not here
    !!
    !! use for temporary values
    real(dp), dimension(4,2) :: langc
    !!    
contains
subroutine init_md()
implicit none
    integer :: n_init
    !!
    n_init = 0
    inquire(file=mdin_file, exist=my_exist)
    if (my_exist .eqv. .true.) then
        open(unit=10, file=mdin_file)
        do while (.true.)
            read(10,*,iostat=my_iostat) exe_arg2
            if (my_iostat < 0) exit
            select case (exe_arg2(1))
                case ('temp')
                    read(exe_arg2(2),*) md_temp
                    md_beta = 1_dp / md_temp
                case ('sfreq')
                    read(exe_arg2(2),*) md_sysfreq
                case ('dtime')
                    read(exe_arg2(2),*) md_dtime
                case ('coeff')
                    read(exe_arg2(2),*) md_coeff
                case ('nstep')
                    read(exe_arg2(2),*) md_nstep
                case ('methd')
                    read(exe_arg2(2),*) md_methd
                case ('scheme')
                    read(exe_arg2(2),*) md_scheme
                case ('thermo')
                    read(exe_arg2(2),*) md_thermo
                case ('mirror')
                    read(exe_arg2(2),*) md_virtual
                case ('bead')
                    read(exe_arg2(2),*) md_bead
                    n_init = n_init - 1
                case ('nsmp')
                    read(exe_arg2(2),*) md_nsmp
                    n_init = n_init - 1
                case ('opt')
                    read(exe_arg2(2),*) md_opt
                    n_init = n_init - 1
                case default
                    n_init = n_init - 1
            end select
            n_init = n_init + 1
        end do
        close(unit=10)
    else
        call send_err('err:init_md, the _put.in is not exist in current directory!')
    end if
    !!
    if (md_opt .eq. 1) then
        md_dtime = 0.2 / md_sysfreq
        md_nstep = int(10000 / md_sysfreq)
        select case (md_thermo)
            case (0)
                md_coeff = md_sysfreq
            case (1)
                md_coeff = sqrt(2.0_dp) * md_sysfreq
            case (2)
                md_coeff = 20 * md_dtime
        end select
    end if
    !!
    langc(1,1) = exp( - md_coeff * md_dtime )                           ! c1=e^(-g*dt), c2=sqrt(1-c1**2), for dt
    langc(1,2) = sqrt( 1 - langc(1,1) * langc(1,1) )                    ! ;
    langc(2,1) = exp( - 0.5_dp * md_coeff * md_dtime )                     ! c1=e^(-g*dt), c2=sqrt(1-c1**2), for dt/2
    langc(2,2) = sqrt( 1 - langc(2,1) * langc(2,1) )                    ! ;
    langc(3,1) = exp( - 0.25_dp * md_coeff * md_dtime )                     ! c1=e^(-g*dt), c2=sqrt(1-c1**2), for dt/4
    langc(3,2) = sqrt( 1 - langc(3,1) * langc(3,1) )                    ! ;
    langc(4,1) = exp( - 0.125_dp * md_coeff * md_dtime )                     ! c1=e^(-g*dt), c2=sqrt(1-c1**2), for dt/8
    langc(4,2) = sqrt( 1 - langc(4,1) * langc(4,1) )                    ! ;
    !!
    md_npass = 0
    md_iseqb = 0
    md3_bead = 3*md_bead
    md_bfreq = sqrt( real(md_bead) ) * md_temp
    md_bfreq2 = md_bfreq**2
    !!
    if (n_init < 8) then
        call send_err('err:init_md, the number of file-parameters is mismatch or less!')
    end if
end subroutine init_md

subroutine init_traj(mobj)
implicit none
    type(mole), intent(inout) :: mobj 
    real(dp), dimension(9) :: input_reals !----------------------------------, move 9?
    integer :: i, j, k
    !!
    call build_mole(mobj, md_bead)
    !!
    select case (md_mod)
        case(2) ! restart from end.rst, recommend!
            open(unit=10,file='end.rst')
            do i=1,mobj%nb
                do j=1,md_bead
                    read(10,*) input_reals
                    mobj%a(i)%ks(3*j-2) = input_reals(4)
                    mobj%a(i)%ks(3*j-1) = input_reals(5)
                    mobj%a(i)%ks(3 * j) = input_reals(6)
                    mobj%a(i)%p(3*j-2) = input_reals(7)
                    mobj%a(i)%p(3*j-1) = input_reals(8)
                    mobj%a(i)%p(3 * j) = input_reals(9)
                end do
            end do
            close(unit=10)
        case(3) ! restart from e.etj
            if(md_mod3n .le. 0) then
                call send_err("error: init traj from e.xpe failed!")
            end if
            !skip the used lines
            open(unit=10,file='e.etj')
            do k=2,md_mod3n
                do i=1,mobj%nb
                    do j=1,md_bead
                        read(10,*)
                    end do
                end do
            end do
            ! init from md_mod3n-th lines of e.etj
            do i=1,mobj%nb ! here a line an atom
                do j=1,md_bead
                    read(10,*) input_reals
                    mobj%a(i)%ks(3*j-2) = input_reals(4)
                    mobj%a(i)%ks(3*j-1) = input_reals(5)
                    mobj%a(i)%ks(3 * j) = input_reals(6)
                    mobj%a(i)%p(3*j-2) = input_reals(7)
                    mobj%a(i)%p(3*j-1) = input_reals(8)
                    mobj%a(i)%p(3 * j) = input_reals(9)
                end do
            end do
            close(unit=10)
        case default ! for md_mod = 0, or md_mod = 1
            do i=1,mobj%nb
                do j=1,md3_bead
                    mobj%a(i)%ks(j) = 0
                    call random_norm( mobj%a(i)%p(j), rand_throw, sqrt( mobj%a(i)%m * md_temp) )
                end do
            end do
    end select
end subroutine init_traj

subroutine show_traj(mobj)
implicit none
    type(mole), intent(inout) :: mobj 
    integer :: i
    do i=1,mobj%nb
        print *,'-----'
        print *, mobj%a(i)%ks
        print *, mobj%a(i)%p
        print *, mobj%a(i)%fks
        print *,'-----'
    end do
end subroutine show_traj

end module md_info
