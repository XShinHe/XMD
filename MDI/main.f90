!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 07
!--- Acknowledgement to Liu group, of PKU
!*

program xmd
use MyDef
use AM_script
use md_info
use md_rout
implicit none
    type(mole) :: mobj
    character(len=len1), dimension(2) :: pairs
    integer :: n,i,j
    integer :: t1, t2, time_cal
    
    n = command_argument_count()
    if (mod(n,2) .ne. 0) then
        call send_err( "error, arguments mismatch, they should be pairs! try \' [main] -h h\' " )
    else
	    n = n/2
	end if
	
	!-- set default mooes
	md_x_mode = 0          ! general simulation
	md_y_mode = 0          ! output ana_file only
	md_m_mode = 0          ! refresh the start-point
	md_e_mode = 0          ! needn't strike on balancd
	md_sn_mode = 0         ! set n-th line zero
	md_b_mode = 0          ! build (from cnfg.rc) without coordinates
	
	!-- read from terminal
	do i=1,n
		call get_command_argument(2*i-1, pairs(1))
		call get_command_argument(2*i, pairs(2))
		select case (pairs(1))
			case ("-h")
				print *, "handling the arguments, now support following usage : "
				print *, "    main -h  h,   show the help information "
				print *, "    main -x  n,   run n-th trajectory,       need xmd.rst or not"
				print *, "    main -x  e,   run equilibrium,           need nothing"
				print *, "    main -x  r,   run a restart trajectory   need xmd.rst"
				print *, "    main -x  s,   run from a sampling file   need xmd.smp(default)"
				print *, "    main -x  l,   run from a filename list   need xmd.lst(default)"
				print *, "    main -y  i,   the way for smapling                  need xmd.rst"
				print *, "    main -m  i,   from new/old steps (i=0/i=1)          need xmd.rst"
				print *, "    main -e  i,   set mode (i=0/i=1)                    [need -x e/n]"
				print *, "    main -r  *,   alternative for xmd.rst               [need -x r]"
				print *, "    main -s  *,   alternative for xmd.smp               [need -x s]"
				print *, "    main -l  *,   alternative for xmd.lst               [need -x f]"
				print *, "    main -sn n,   sample from a n-th of xmd.smp/xmd.lst [need -x s]"   !-- you'd better not use such property
				print *, "    main -p  *,   alternative for putrc_file"
				print *, "    main -c  *,   alternative for cnfgrc_file"
				print *, "    main -o  *,   namefor out_file"
				print *, "    main -a  *,   namefor ana_file"
				stop
		    case ("-x")
		    	ana_file = trim(pairs(2))//'.ana'
		    	out_file = trim(pairs(2))//'.out'
		        select case (pairs(2))
		            case ('e')
		                md_x_mode = 1
		            case ('r')
		                md_x_mode = 2
		            case ('s')
		                md_x_mode = 3
		            case ('l')
		                md_x_mode = 4
		            case default
		            	md_x_mode = 0
		        end select
		    case ("-y")
				select case (pairs(2))
					case ('0')
						md_y_mode = 0
					case ('1')
						md_y_mode = 1
					case ('2')
						md_y_mode = 2
					case default
						call send_err("error, -y * arguments misch")
				end select
		    case ("-m")
		    	if(md_x_mode .ne. 0 .and. md_x_mode .ne. 2) then
		    		call send_err('error, -m * arguments need restart mode')
		    	end if
				select case (pairs(2))
					!-- refresh the step counting
					case ('0')
						md_m_mode = 0
					!-- use the old step counting
					case ('1')
						md_m_mode = 1
					case default
						call send_err('error, -m * arguments misch')
				end select
		    case ("-e")
		    	if(md_x_mode .ne. 1) then
		    		call send_err('error, -e * arguments need restart mode')
		    	end if
				select case (pairs(2))
					!-- don't rarely care about equilibrium
					case ('0')
						md_e_mode = 0
					!-- must srike on the equilibrium, use stattistic tools
					case ('1')
						md_e_mode = 1
					!-- 
					case ('2')
						md_e_mode = 2
					case default
						call send_err("error, -e * arguments misch")
				end select
		    case ("-r")
		    	if(md_x_mode .ne. 2) then
		    		call send_err('error, -r * arguments need restart mode')
		    	end if
		        rst_file = trim(pairs(2))
		    case ("-s")
		    	if(md_x_mode .ne. 3) then
		    		call send_err('error, -s * arguments need restart mode')
		    	end if
				smp_file = trim(pairs(2))
			case ("-l")
				if(md_x_mode .ne. 4) then
		    		call send_err('error, -l * arguments need restart mode')
		    	end if
				lst_file = trim(pairs(2))
		    case ("-sn")
		    	if(md_x_mode .ne. 3 .and. md_x_mode .ne. 4) then
		    		call send_err('error, -n * arguments need restart mode')
		    	end if
		        read(pairs(2),*) md_sn_mode
		    case ("-p")
		        putrc_file = trim(pairs(2))
		    case ("-c")
		        cnfgrc_file = trim(pairs(2))
		    case ("-o")
		       	out_file = trim(pairs(2))
		    case ("-a")
		        ana_file = trim(pairs(2))
		    case default
		        call send_err("error, arguments mismarch")
		end select
	end do
	
	!-- initializzation rondom seed by date-and-time
	call init_seed()
	
	!-- timer 1
	call system_clock(t1)
	!-- for initialization control parameters of molecular dynamics
    call init_md()
    !-- build the molecular configuration, and generate/read initial positions and momenta
    call init_traj(mobj)
    !-- timer 2
    call system_clock(t2)
    time_cal = t2 - t1
    print *, "initialization time: ", time_cal
    
    !-- timer 1
	call system_clock(t1)
    !-- Molecular dynamics simulations
    call md_ctrl(mobj)
    !-- timer 2
    call system_clock(t2)
    time_cal = t2 - t1
    print *, "running time: ", time_cal
    
end program xmd

