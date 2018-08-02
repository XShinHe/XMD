!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 08
!--- Acknowledgement to Liu group, of PKU
!*

module md_info
use MyDef
use AM_script
implicit none
    !-- Recording number(real/integer) information, using preffix 'md'
    !-- 1) Esemble parameters
    integer  :: md_nsum    ! N
    real(dp) :: md_esum    ! E
    real(dp) :: md_temp    ! T
    real(dp) :: md_beta    ! 1/T
    real(dp) :: md_pres    ! P
    real(dp) :: md_volm    ! V
    
    !-- 2) Simulation and sampling parameters
    real(dp) :: md_dtime
    integer  :: md_nstep
    integer  :: md_npass
    integer  :: md_nsamp
    
    !-- 3) Path integral technique parameters
    integer  :: md_ipimd                 ! path integral mode [ 0: primitive / 1: staging / 2: normal-mode ]
    integer  :: md_bead                  ! default vaule
    real(dp) :: md_bfreq                 ! default
    real(dp) :: md_bfreq2
    real(dp) :: md_gamma                 ! for PILD, adiabatic parameter
    
    !-- 4) Thermostat and spliting parameters
    integer  :: md_ithermo               ! thermostat parameter [ Langevin, Anderson, Nose-Hoover Chain]
    integer  :: md_ischeme               ! spilting parameter [ VV,LF,PC,RK / middle, end, side]
    real(dp) :: md_coeff                 ! (parameter describing the friction/collision process)
	real(dp), dimension(4,2) :: langc    ! parameters of Langevin/Anderson

	!-- 5) Other parameters
	!----- for example, the vir/real process (NOT the IMAGINARY/REAL-TIME process)
    integer  :: md_ivir
	!----- and for example, the optimization to initialization
    integer  :: md_iopt
    real(dp) :: md_cfreq                 ! giving an approximation of system's character frequency for optimizattion
    

    !-- 6) execute mode parameters
    integer  :: md_x_mode                ! 0, general; 1, equilibrium; 2, restart; 3, fr. sampling; 4, fr. list
    integer  :: md_y_mode                ! 0, output ana_file; 1, output ana and out file
    integer  :: md_m_mode                ! 0, set new npass to restart; 1, from last npass to restart
    integer  :: md_e_mode                ! 0 needn't strike equilibrium; 1, need that.
    integer  :: md_sn_mode               ! n-th example of sampling file
    integer  :: md_b_mode                ! build mode from the cnfg.rc
   
contains
	
	!-- initialize parameters from putrc_file
	subroutine init_md()
	implicit none
		integer :: n_init
		character(len=len1), dimension(2) :: pairs

		!-- set default parameters
		md_ipimd   = 1    ! staging path integral mode
		md_ischeme = 0    ! VV-middle
		md_ithermo = 0    ! Langevin thermostat
		md_ivir    = 0    ! don't use vir-dynamics
		md_iopt    = 0    ! don't optimize parameters
		md_bead    = 1    ! classical molecular dynamics
		md_nsamp   = 100  ! default sampling step
		
		!-- read and count the parameters
		n_init = 0
		inquire( file=trim(putrc_file), exist=my_exist )
		if (my_exist .eqv. .true.) then
		    open( unit=10, file=trim(putrc_file) )
		    do while (.true.)
		        read(10,*,iostat=my_iostat) pairs
		        if (my_iostat < 0) exit
		        select case (pairs(1))
		            case ('temp')
		                read(pairs(2),*) md_temp
		                md_beta = 1_dp / md_temp
		            case ('dtime')
		                read(pairs(2),*) md_dtime
		            case ('nstep')
		                read(pairs(2),*) md_nstep
		            case ('nsamp')
		                read(pairs(2),*) md_nsamp
		            case ('thermo')
		                read(pairs(2),*) md_ithermo
		            case ('scheme')
		                read(pairs(2),*) md_ischeme  
		            case ('coeff')
		                read(pairs(2),*) md_coeff
		            case ('pimode')
		                read(pairs(2),*) md_ipimd
		            case ('bead')
		                read(pairs(2),*) md_bead
		                n_init = n_init - 1
		            case ('gamma')
		                read(pairs(2),*) md_gamma
		                n_init = n_init - 1
		            case ('vir')
		                read(pairs(2),*) md_ivir
		            case ('opt')
		                read(pairs(2),*) md_iopt
		                n_init = n_init - 1
		            case ('cfreq')
		                read(pairs(2),*) md_cfreq
		                n_init = n_init - 1
		            case default
		                n_init = n_init - 1
		        end select
		        n_init = n_init + 1
		    end do
		    close(unit=10)
		else
		    call send_err('err:init_md, the put.rc is not exist in current directory!')
		end if
		
		!-- if allow to optimize
		if (md_iopt .eq. 1) then
		    md_dtime = 0.2 / md_cfreq
		    md_nstep = 10000 / md_cfreq
		    select case (md_ithermo)
		        case (0)
		            md_coeff = md_cfreq
		        case (1)
		            md_coeff = sqrt(2.) * md_cfreq
		        case (2)
		            md_coeff = 20 * md_dtime
		    end select
		end if

		!-- calculate temporary valuables for thermostat
		langc(1,1) = exp( - md_coeff * md_dtime )                           ! c1=e^(-g*dt), c2=sqrt(1-c1**2), for dt
		langc(1,2) = sqrt( 1 - langc(1,1) * langc(1,1) )                    ! ;
		langc(2,1) = exp( - 0.5_dp * md_coeff * md_dtime )                  ! c1=e^(-g*dt), c2=sqrt(1-c1**2), for dt/2
		langc(2,2) = sqrt( 1 - langc(2,1) * langc(2,1) )                    ! ;
		langc(3,1) = exp( - 0.25_dp * md_coeff * md_dtime )                 ! c1=e^(-g*dt), c2=sqrt(1-c1**2), for dt/4
		langc(3,2) = sqrt( 1 - langc(3,1) * langc(3,1) )                    ! ;
		langc(4,1) = exp( - 0.125_dp * md_coeff * md_dtime )                ! c1=e^(-g*dt), c2=sqrt(1-c1**2), for dt/8
		langc(4,2) = sqrt( 1 - langc(4,1) * langc(4,1) )                    ! ;
		
		!-- calculte some simplified valuables
		md_npass = 0
		md_bfreq = sqrt( real(md_bead) ) * md_temp
		md_bfreq2 = md_bfreq**2
		
		!-- check the initialization of md is whether complete
		if (n_init < 8) then
		    call send_err('err:init_md, the number of file-parameters is mismatch or less!')
		end if
	end subroutine init_md

	!-- creat the moleculer object
	!------ now the moleculer object is like a type, though they may not perform well in numerical calculation,
	!------ so, an alternative choise is using all-array to construct a molecular object, this might be a better choise
	subroutine build_mole(mobj)
    implicit none
        type(mole), intent(inout) :: mobj
        character(len=len3) :: ins_s
        character(len=len0) :: ins_c
        real(dp) :: ins_r, ins_rs
        integer :: ins_i, ins_j
        integer :: i,j
        
        inquire(file=cnfgrc_file, exist=my_exist)
        if (my_exist .eqv. .false.) then
            call send_err("error, configuration file (cnfg.rc) opening wrong")
        end if
        open(unit=10,file=cnfgrc_file)
        read(10,*) ins_s, ins_i, ins_j
        if (ins_s(1:1) .ne.'!') then
            call send_err('error: format error of cnfg.rc')
        end if
        mobj%iname=ins_s(2:len(ins_s))
        md_nsum=ins_i
        
        if(ins_j > 0) then
        	md_b_mode = 1
        else
        	md_b_mode = 0
        end if
        
        !-- allocations
        allocate(mobj%a( md_nsum ))
        do i=1,md_nsum
            allocate( mobj%a(i)%ks( md_bead ) )
            allocate( mobj%a(i)%x( md_bead) )
            allocate( mobj%a(i)%p( md_bead ) )
            allocate( mobj%a(i)%fx( md_bead) )
            allocate( mobj%a(i)%fks( md_bead) )
        end do
        
        !-- first initialization of position
        do i=1,md_nsum
        	if(ins_j .eq. 0) then
	            read(10,*,iostat=my_iostat) ins_c, ins_r
	        else if(ins_j > 0) then
	        	read(10,*,iostat=my_iostat) ins_c, ins_r, ins_rs
	        end if
            if(my_iostat<0) exit
            mobj%a(i)%iname = ins_c
    	    mobj%a(i)%m = ins_r
    	    if(ins_j>0) then
    	    	mobj%a(i)%x = ins_rs
	            !do j=1,md_bead
    	        !	mobj%a(i)%x(j) = ins_rs  !-------------- ugly
    	        !end do
    		end if
        end do
        close(unit=10)    
    end subroutine build_mole

	!-- initialization of a trajectory (with position and momenta)
	subroutine init_traj(mobj)
	implicit none
		type(mole), intent(inout) :: mobj 
		real(dp), dimension(8) :: ins_r !----------------------------------, move 8?
		character(len=10) :: ins_c
		integer :: i, j, k
		
		!-- building the molecule object
		call build_mole(mobj)

		select case (md_x_mode)
			!-- restart from nothing/cnfg.rc
			case(1)
		        do i=1,md_nsum
		            do j=1,md_bead
		                if(md_b_mode .eq. 0) then 
		                	mobj%a(i)%x(j) = 0
		                end if
		                call random_norm(mobj%a(i)%p(j), rand_throw, sqrt( mobj%a(i)%m * md_temp) )
		            end do
		        end do
			!-- restart from xmd.rst
		    case(2) 
		    	inquire(file=trim(rst_file), exist=my_exist)
        		if (my_exist .eqv. .false.) then
            		call send_err("error, rst_file doesn't exist")
        		end if
        		
		        open(unit=10,file=trim(rst_file))
		        do i=1,md_nsum
		            do j=1,md_bead
		                read(10,*) ins_r
		                mobj%a(i)%ks(j)=ins_r(4)
		                mobj%a(i)%x(j)=ins_r(5)
		                mobj%a(i)%p(j)=ins_r(6)
		                mobj%a(i)%fx(j)=ins_r(7)
		                mobj%a(i)%fks(j)=ins_r(8)
		            end do
		        end do
		        if(md_m_mode .eq. 1) then
		             md_npass = md_npass + int(ins_r(1))
		        end if
		        close(unit=10)
		    !-- restart from xmd.smp
		    case(3) 
		        if(md_sn_mode .le. 0) then
		            call send_err("error: init traj from xmd.smp/xmd.lst failed!")
		        end if
		        !-- skip the used lines
		        open(unit=10,file=trim(smp_file))
		        do k=2,md_sn_mode
		            do i=1,md_nsum
		                do j=1,md_bead
		                    read(10,*)
		                end do
		            end do
		        end do
		        !-- initial from n-th example of smp_file
		        do i=1,md_nsum
		            do j=1,md_bead
		                read(10,*) ins_r
		               	mobj%a(i)%ks(j)=ins_r(4)
		                mobj%a(i)%x(j)=ins_r(5)
		                mobj%a(i)%p(j)=ins_r(6)
		                mobj%a(i)%fx(j)=ins_r(7)
		                mobj%a(i)%fks(j)=ins_r(8)
		            end do
		        end do
		        close(unit=10)
		    !-- restart from xmd.lst
		    case(4)
		    	if(md_sn_mode .le. 0) then
		            call send_err("error: init traj from xmd.smp/xmd.lst failed!")
		        end if
		        
		        !-- read restart filename from lst_file
		    	inquire(file=lst_file, exist=my_exist)
        		if (my_exist .eqv. .false.) then
            		call send_err("error, lst_file opening wrong")
        		end if
        		open(unit=10,file=lst_file)
        		do k=2,md_sn_mode
        			read(10,*)
        		end do
        		read(10,*) rst_file
        		close(unit=10)
        		
        		!-- run from rst_file
        		inquire(file=trim(rst_file), exist=my_exist)
        		if (my_exist .eqv. .false.) then
            		call send_err("error, lst_file record wrong information")
        		end if
        		open(unit=10,file=trim(rst_file))
		        do i=1,md_nsum
		            do j=1,md_bead
		                read(10,*) ins_r
		                mobj%a(i)%ks(j)=ins_r(4)
		                mobj%a(i)%x(j)=ins_r(5)
		                mobj%a(i)%p(j)=ins_r(6)
		                mobj%a(i)%fx(j)=ins_r(7)
		                mobj%a(i)%fks(j)=ins_r(8)
		                if(md_m_mode .eq. 1) then
		                	md_npass = md_npass + int(ins_r(1))
		                end if
		            end do
		        end do
		        close(unit=10)
		    case default
		    	!-- first try to initial from cnfg.rc
		        do i=1,md_nsum
		            do j=1,md_bead
		                if(md_b_mode .eq. 0) then 
		                	mobj%a(i)%x(j) = 0
		                end if
		                call random_norm(mobj%a(i)%p(j), rand_throw, sqrt( mobj%a(i)%m * md_temp) )
		            end do
		        end do
		        
		        !-- if there is restart we can use, we then try to initial from rst_file
		        inquire(file=trim(rst_file), exist=my_exist)
        		if (my_exist .eqv. .true.) then
            		open(unit=10,file=trim(rst_file))
				    do i=1,md_nsum
				        do j=1,md_bead
				            read(10,*) ins_r
				            mobj%a(i)%ks(j)=ins_r(4)
				            mobj%a(i)%x(j)=ins_r(5)
				            mobj%a(i)%p(j)=ins_r(6)
				            mobj%a(i)%fx(j)=ins_r(7)
				            mobj%a(i)%fks(j)=ins_r(8)
				            if(md_m_mode .eq. 1) then
				            	md_npass = md_npass + int(ins_r(1))
				            end if
				        end do
				    end do
		        	close(unit=10)
        		end if	
		end select
		
	end subroutine init_traj

end module md_info
