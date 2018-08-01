!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 07
!--- Acknowledgement to Liu group, of PKU
!*

module md_info
use MyDef
use AM_script
implicit none
    !-- esemble settings
    real(dp) :: md_nsum    ! N
    real(dp) :: md_esum    ! E
    real(dp) :: md_temp    ! T
    real(dp) :: md_beta    ! 1/T
    real(dp) :: md_pres    ! P
    real(dp) :: md_volm    ! V
    integer  :: md_dim = 3
    
    !-- thermostat parameters (friction/collide argumnet)
    integer  :: md_thermo
    integer  :: md_scheme
    real(dp) :: md_coeff

	!-- dynamic parameters
    integer  :: md_methd ! primitive, staging, normal-mode
    integer  :: md_virtual
    integer  :: md_bead = 1              ! default vaule
    integer  :: mdn_bead                 ! consider the dimensions
    real(dp) :: md_bfreq                 ! default
    real(dp) :: md_bfreq2
    real(dp) :: md_gamma                 ! for PILD

    !-- for numerical simulation parameters
    integer  :: md_nstep
    integer  :: md_npass
    real(dp) :: md_dtime
    integer  :: md_nsmp = 1000  

    !-- execute mode parameters
    integer  :: md_x_mode = -1 ! 0, general; 1, eqb; 2, restart; 3, from sampling; 4, from list
    integer  :: md_y_mode = -1 ! 0, only output ana_file; 1, output ana and out file
    integer  :: md_m_mode = 0  ! clear to a new start step
    integer  :: md_e_mode = 0  ! 0 needn't strike equilibrium; 1, needs.
    integer  :: md_sn_mode = 0 ! n-th example of sampling file
    integer  :: md_b_mode = 0  ! build from the cnfg.rc

    !-- optimization controlling parameter
    integer  :: md_opt = 0
    real(dp) :: md_sysfreq                    ! mass and atomfreq are should be arranged to cnfg.in file, not here

    !-- for temporary values
    real(dp), dimension(4,2) :: langc
   
contains

	subroutine init_md()
	implicit none
		integer :: n_init
		character(len=len1), dimension(2) :: pairs

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
		            case ('sfreq')
		                read(pairs(2),*) md_sysfreq
		            case ('dtime')
		                read(pairs(2),*) md_dtime
		            case ('coeff')
		                read(pairs(2),*) md_coeff
		            case ('nstep')
		                read(pairs(2),*) md_nstep
		            case ('methd')
		                read(pairs(2),*) md_methd
		            case ('scheme')
		                read(pairs(2),*) md_scheme
		            case ('thermo')
		                read(pairs(2),*) md_thermo
		            case ('mirror')
		                read(pairs(2),*) md_virtual
		            case ('bead')
		                read(pairs(2),*) md_bead
		                n_init = n_init - 1
		            case ('nsmp')
		                read(pairs(2),*) md_nsmp
		                n_init = n_init - 1
		            case ('opt')
		                read(pairs(2),*) md_opt
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
		if (md_opt .eq. 1) then
		    md_dtime = 0.2 / md_sysfreq
		    md_nstep = 10000 / md_sysfreq
		    select case (md_thermo)
		        case (0)
		            md_coeff = md_sysfreq
		        case (1)
		            md_coeff = sqrt(2.) * md_sysfreq
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
		mdn_bead = md_bead * md_dim
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
        real(dp) :: ins_r
        real(dp), dimension(:), allocatable :: ins_rs
        integer :: ins_i, ins_j
        integer :: i,j,d
        
        inquire(file=cnfgrc_file, exist=my_exist)
        if (my_exist .eqv. .false.) then
            call send_err("error, configuration file (cnfg.rc) opening wrong")
        end if
        open(unit=10,file=cnfgrc_file)
        read(10,*) ins_s, ins_i, ins_j
        if (ins_s(1:1) .ne.'!') then
            call send_err('error: format error of cnfg.rc')
        end if
        mobj%nm=ins_s(2:len(ins_s))
        mobj%nb=ins_i
        
        if(ins_j > 0) then 
        	allocate(ins_rs(ins_j))
        	md_b_mode = 1
        	md_dim = ins_j
        	mdn_bead = md_dim * md_bead
        else
        	md_b_mode = 0
        	!-- so here position should next to initialize
        	!-- so here md_dim using the default settings
        end if
        
        !-- allocations
        allocate(mobj%a( mobj%nb ))
        do i=1,mobj%nb
            allocate( mobj%a(i)%x( mdn_bead ) )
            allocate( mobj%a(i)%ks( mdn_bead) )
            allocate( mobj%a(i)%p( mdn_bead ) )
            allocate( mobj%a(i)%fx( mdn_bead) )
            allocate( mobj%a(i)%fks( mdn_bead) )
        end do
        
        !-- pre-initialization of position
        do i=1,mobj%nb
        	if(ins_j .eq. 0) then
	            read(10,*,iostat=my_iostat) ins_c, ins_r
	        else if(ins_j > 0) then
	        	read(10,*,iostat=my_iostat) ins_c, ins_r, ins_rs
	        end if
            if(my_iostat<0) exit
            mobj%a(i)%elem = ins_c
    	    mobj%a(i)%m = ins_r
    	    if(ins_j>0) then
	            do d=1,md_dim
	            	do j=1,md_bead
    	        		mobj%a(i)%x(md_dim*(j-1)+d) = ins_rs(d)
    	        	end do
    		    end do
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
		        do i=1,mobj%nb
		            do j=1,mdn_bead
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
		        do i=1,mobj%nb
		            do j=1,mdn_bead
		                read(10,*) ins_r
		                mobj%a(i)%ks(j)=ins_r(4)
		                mobj%a(i)%p(j)=ins_r(5)
		                mobj%a(i)%x(j)=ins_r(6)
		                mobj%a(i)%fks(j)=ins_r(7)
		                mobj%a(i)%fx(j)=ins_r(8)
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
		            do i=1,mobj%nb
		                do j=1,mdn_bead
		                    read(10,*)
		                end do
		            end do
		        end do
		        !-- initial from n-th example of smp_file
		        do i=1,mobj%nb
		            do j=1,mdn_bead
		                read(10,*) ins_r
		                mobj%a(i)%ks(j)=ins_r(4)
		                mobj%a(i)%p(j)=ins_r(5)
		                mobj%a(i)%x(j)=ins_r(6)
		                mobj%a(i)%fks(j)=ins_r(7)
		                mobj%a(i)%fx(j)=ins_r(8)
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
		        do i=1,mobj%nb
		            do j=1,mdn_bead
		                read(10,*) ins_r
		                mobj%a(i)%ks(j)=ins_r(4)
		                mobj%a(i)%p(j)=ins_r(5)
		                mobj%a(i)%x(j)=ins_r(6)
		                mobj%a(i)%fks(j)=ins_r(7)
		                mobj%a(i)%fx(j)=ins_r(8)
		                if(md_m_mode .eq. 1) then
		                	md_npass = md_npass + int(ins_r(1))
		                end if
		            end do
		        end do
		        close(unit=10)
		    case default
		    	!-- first try to initial from cnfg.rc
		        do i=1,mobj%nb
		            do j=1,mdn_bead
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
				    do i=1,mobj%nb
				        do j=1,mdn_bead
				            read(10,*) ins_r
				            mobj%a(i)%ks(j)=ins_r(4)
				            mobj%a(i)%p(j)=ins_r(5)
				            mobj%a(i)%x(j)=ins_r(6)
				            mobj%a(i)%fks(j)=ins_r(7)
				            mobj%a(i)%fx(j)=ins_r(8)
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
