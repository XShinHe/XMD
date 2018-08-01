!*
! ---Copyright by--- He Xin <1500011805@pku.edu.cn>
!*

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
        call send_err( "error, arguments mismatch" )
    end if

	! handling the arguments, now support following usage
	! main -x e,   equilibrium mode (fisrt run)
	! main -x r,   restart mode
	! main -f n,   from sampling file to start 
    call get_command_argument(1,exe_arg2(1))
    call get_command_argument(2,exe_arg2(2))
    anal_file = trim(exe_arg2(2))//'.ana'
    select case (exe_arg2(1))
        case ("-x")
            select case (exe_arg2(2))
                case ('e')
                    md_mod=1
                    traj_outfile = trim(exe_arg2(2))//'.xpe'
                    anal_file = trim(exe_arg2(2))//'.ena'
                case ('r')
                    md_mod=2
                    traj_outfile = trim(exe_arg2(2))//'.xst'
                case default
                    traj_outfile = trim(exe_arg2(2))//'.xpf'
                    md_mod=0
            end select
        case ("-f")
            md_mod = 3
            read(exe_arg2(2),*) md_mod_fn
            traj_outfile = trim(exe_arg2(2))//'.xpf'
        case default
            call send_err("error, arguments mismarch")
    end select
	
	!------------- do some initialization --------------------
	!-- initialize control parameters of molecular dynamics
    call init_md()
    !-- initialize rondom seed by date-and-time
    call init_seed()
    !-- build the molecular configuration, and generate/read initial positions and momenta
    call init_traj(mobj)
    
    !------------- do MD task --------------------------------
    call md_ctrl(mobj)
    
end program main

