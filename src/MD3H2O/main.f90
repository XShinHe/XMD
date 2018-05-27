! --- Copyright by Shin He <hx0824916@pku.edu.cn> ---

program main
use MyDef
use AM_script
use md_info
use md_rout
implicit none
    type(mole) :: mobj
    integer :: n
    n = command_argument_count()
    if ((n .ne. 2)) then
        call send_err( "error, arguments mismatch" )    ! filename input the parameters in as 2nd args
    end if

    call get_command_argument(1,exe_arg2(1))
    call get_command_argument(2,exe_arg2(2))
    enrg_file = trim(exe_arg2(2))//'.erg'
    trsf_file = trim(exe_arg2(2))//'.trs'
    else_file = trim(exe_arg2(2))//'.els'
    anal_file = trim(exe_arg2(2))//'.ana'
    select case (exe_arg2(1))
        case ("-x")
            select case (exe_arg2(2))
                case ('e')
                    md_mod=1
                    traj_file = trim(exe_arg2(2))//'.etj'
                case ('r')
                    md_mod=2
                    traj_file = trim(exe_arg2(2))//'.rtj'
                case default
                    md_mod=0
                    traj_file = trim(exe_arg2(2))//'.xtj'
            end select
        case ("-f")
            md_mod = 3
            read(exe_arg2(2),*) md_mod3n
            traj_file = trim(exe_arg2(2))//'.ftj'
        case default
            call send_err("error, arguments mismarch")
    end select
!---------------------------------------- next arrangement -------------------------------

    call init_md()
    call init_seed()
    call init_traj(mobj)
    call md_ctrl(mobj)
    
end program main


