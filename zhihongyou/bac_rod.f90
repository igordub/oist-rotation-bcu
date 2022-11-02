module bac_rod_data
    implicit none
    
    ! Simulation control:
    integer,parameter:: dimen=2           ! Dimension of data, must be no less than 2.
    integer,parameter:: dimen_effe=2      ! Real dimension in the model, could be different from dimen.
    integer,parameter:: n_max=50000          ! Maximum number of cells, sets the size of x et al.
    ! Global switches:
    integer,parameter:: init_cond=1       ! 0-import, 1-one cell, 2-chain, 3-box-seed, 4-box-dumb, 5-battle
    integer,parameter:: boundary_cond=5   ! 0-free, 1-fixed, 2-periodic, 3-absorb+periodic, 4-periodic+dynamic, 5-absorb+fixed
    integer,parameter:: timing_scheme=1   ! Timing scheme during simulation: 1-time, 2-cell number, 3-cell division, 4/5-x/y extension
    integer,parameter:: stop_nmax=1       ! Stop simulation when cell number is n_max.
    integer,parameter:: stop_lmax=1       ! Stop simulation when colony gets the critical size.
    integer,parameter:: stop_buckle=0     ! Stop simulation when buckles.
    integer,parameter:: adap_dt=0         ! Time step control: 1-adaptive time step, 2-dt=dt_init.
    integer,parameter:: i_check_buckle=0     ! Stop simulation when buckles.
    integer,parameter:: print_runtime_status=1
    integer,parameter:: meas_stress=1     ! Measure stress or not.
    integer,parameter:: meas_force=0     ! Measure forces or not.
    integer,parameter:: exp_conf=1        ! Export configuration or not.
    integer,parameter:: exp_meas=1        ! Export measurement.
    integer,parameter:: exp_resu=0        ! Export other results at the final time step.
    integer,parameter:: dele_data=1       ! Delete data files at the beginning of simulation.
    integer,parameter:: rod_or_bead=1     ! Type of interaction: 1-rod, 2-bead, 3-combined.
    integer,parameter:: grow_type=1       ! Type of growth: 1-normal, 2-force dependent, 3-times g_time
    integer,parameter:: substrate_on=0    ! Turn on substrate.
    integer,parameter:: compress_on=0     ! Turn on compress force.
    integer,parameter:: propel_on=0       ! Turn on propel force.
    integer,parameter:: dumb_fixed=0      ! To fix the positions of dumb cells or not.
    integer,parameter:: random_seed_time=1  ! 0-fixed random seed, 1-initial random seed from system time.
    ! Timing
    real(8),parameter:: t_max=40.0         ! Time to stop simulations, valid when timing_scheme=1.
    real(8),parameter:: t0=0              ! Time to start export results, valid when timing_scheme=1.
    real(8),parameter:: dt_exp=0.02       ! Export when mod(t,dt_exp)=0, valid when timing_scheme=1.
    integer,parameter:: n0=0              ! Start export results at n=n0, valid when timing_scheme=2.
    integer,parameter:: dn_exp=100           ! Export when mod(n,dn_exp)=0, valid when timing_scheme=2.
    integer,parameter:: n0_divi=0         ! Number of new cells to start export, valid when timing_scheme=3.
    integer,parameter:: n_max_divi=100000 ! Number of new cells to stop simulation, valid when timing_scheme=3.
    real(8),parameter:: y_c0=20.0            ! Minimal colony size to export, valid when timing_scheme=4/5.          
    real(8),parameter:: dl_exp=1.0        ! Step of colony size between export, valid when timing_scheme=4/5.          
    integer,parameter:: n_max_clear=50000 ! Stop simulation when n_max_clear cells are cleared.
    real(8),parameter:: rho_max=10.5      ! Stop simulation when average cell density exceeds rho_max.
    ! Time step.
    real(8),parameter:: dt_init=2*0.1**6    ! Initial time step.
    real(8),parameter:: dt_max=0.1_8**5   ! Maximum time step.
    real(8),parameter:: dt_min=0.1_8**7   ! Minimum time step.
    real(8),parameter:: dx_maxc=0.0005_8   ! Maximum distance between iterations, used to control dt.
    ! Space control
    real(8),parameter:: l_box_x=70.0        ! Size of box in x.
    real(8),parameter:: l_box_y=50.0        ! Size of box in y.
    real(8),parameter:: dl_box=5.0          ! Length of mini-box in box decomposition.
    integer,parameter:: nx_box=floor(l_box_x/dl_box+0.001)+4  ! Number of mini-box in x.
    integer,parameter:: ny_box=floor(l_box_y/dl_box+0.001)+4  ! Number of mini-box in y.
    real(8),parameter:: v_box_x=0         ! Change of system size in unit time.
    real(8),parameter:: v_box_y=0        ! Change of system size in unit time.
    ! Initial condition
    real(8),parameter:: t_init=0          ! Time at the beginning of simulation.
    integer,parameter:: n_init=n_max          ! Cell number at the beginning of simulation, valid when init_cond>=2.
    character(len=80):: file_impo='../run2/conf_408.dat' ! file to import.
    integer,parameter:: impo_2to3=0       ! Import 2D data to 3D system.
    integer,parameter:: n_dumb=0          ! Number non-growing cells, fixed if dumb_fixed=1, valid when init_cond=4.
    integer,parameter:: n_seed=n_max          ! Number of seed cells, valid when init_cond=3.
    integer,parameter:: n_left=n_init/2          ! Number of seed cells, valid when init_cond=3.
    integer,parameter:: arra_leng=2       ! 1-random length, 2-maximum length, 3-maximum length but central cell.
    real(8),parameter:: l_dev=2.0         ! Deviation of cell length, valid when init_cond=3 and arra_leng=1.
    real(8),parameter:: l_cen=4.0         ! Length of central cell, valid when init_cond=3 and arra_leng=3.
    ! Export data and measurement
    character(len=80):: dire_expo='./'    ! directory to export
    integer,parameter:: t_measure=10      ! Number of time steps to get average measured data.
    integer,parameter:: stress_xy_n=1     ! Express stress in 1-xy coordinate, 2-director coordinate.
    integer,parameter:: n_force_max=100        ! Maximum number of forces when measure force of each cell
    
    ! Model parameters:
    ! Cell size
    real(8),parameter:: d_cell=1.0        ! Cell diameter.
    real(8),parameter:: l_max=2.0           ! Maximum cell length.
    integer,parameter:: nb_max=10         ! Maximum number of beads, valid when rod_or_bead=2&3.
    real(8),parameter:: dl_bead=l_max/(nb_max-1_8)  ! Distance between beads, valid when rod_or_bead=2&3.
    ! Growth
    real(8),parameter:: g_ave=2.0           ! Average growth rate.
    real(8),parameter:: g_time=5          ! Multiple of growth rate when grow_type=3.
    ! Cell-cell interaction.
    real(8),parameter:: e0=20000          ! Rigidity of cells.
    real(8),parameter:: power_fcc=1.5_8
    real(8),parameter:: k_cohe=0          ! Strength of cohesive force among cells.
    real(8),parameter:: r_cohe=0.1*d_cell ! Maximum range of cohesive force among cells.
    real(8),parameter:: r_cutoff=d_cell+r_cohe  ! Cut-off distance when calculating cell-cell interaction.
    ! Cell-substrate interaction.
    real(8),parameter:: e_subs=100000      ! Rigidity of the substrate.
    real(8),parameter:: power_fcs=1
    real(8),parameter:: e_agro=1000       ! Rigidity of the top agarose.
    real(8),parameter:: k_adhe=500.0         ! Strength of adhesive force per unit length from substrate.
    real(8),parameter:: power_adhe=1
    real(8),parameter:: k_grav_cent=10     ! Strength of gravation at cell center.
    real(8),parameter:: k_grav_head=0     ! Strength of gravation on the two ends.
    real(8),parameter:: r_adhe=0.01*d_cell ! Maximum range of adhesive force from substrate.
    ! Other forces.
    real(8),parameter:: f_compress=10000.0     ! Magnitude of the compress force.
    real(8),parameter:: vari_fc=0         
    real(8),parameter:: prop_ave=50       ! Average value of the propel force.
    real(8),parameter:: prop_dev=prop_ave*0.0        ! Deviation of the propel force.
    real(8),parameter:: kapp=5*d_cell     ! Not used at this moment.
    
    ! Other parameters
    character(len=80):: format_real_screen="(12f10.2)"   ! Format of number to print on screen.
    real(8),parameter:: pi=3.1415926535897932384626433_8
    real(8),parameter:: zeros=0.1_8**9
  
    ! List of global variables during run time:
    ! Configurations
    real(8),dimension(dimen,n_max):: x, q, dx, dq   ! Position and orientation
    real(8),dimension(n_max):: l, g       ! Cell length and growth rate.
    integer,dimension(n_max):: n_bead      ! Number of beads in each cell, valid when rod_or_bead>=2.
    real(8),dimension(dimen,nb_max,n_max):: xb  ! Positions of beads.
    real(8),dimension(n_max):: prop       ! Propel force on each cell.
    real(8),dimension(nx_box,ny_box)::p_first  ! First cell in each mini-box.
    real(8),dimension(n_max):: p_next     ! Next cell in mini-box.
    ! System status.
    integer:: n                           ! Number of cells.
    real(8):: t                           ! Time.
    real(8):: dt                          ! Time step.
    integer:: n_divi                      ! Number of new born cells.
    integer:: bx_min,bx_max,by_min,by_max ! Maximum(minimum) indices of non-empty mini-box.
    integer:: bb_l,bb_r,bb_b,bb_t         ! Indices of boxes where the boundary is
    real(8):: l_box_xt                    ! Real box size at run time
    real(8):: l_box_yt
    ! Simulation status.
    integer:: i_stop                      ! Stop simulation when i_stop=1.
    integer:: i_error                     ! Found error when i_error=1.
    integer:: i_buckle                    ! System buckles when i_buckle=1.
    integer:: p_buckle                    ! Which cell buckles first.
    integer:: i_measure_stress                   ! Measure stress when i_measure_stress=1.
    integer:: i_measure_force                   ! Measure force when i_measure_force=1.
    ! Measurement
    real(8),dimension(n_max):: f_axis     ! Axial forces.
    real(8),dimension(4,n_max):: stress_n_t  ! Stress at a time step.
    real(8),dimension(4,n_max):: stress_n    ! Average stress over t_measure time steps.
    real(8),dimension(dimen,n_max):: x0      ! Positions at previous time stamp
    real(8),dimension(dimen,n_max):: q0      ! Orientations at previous time stamp
    real(8),dimension(n_max):: l0      ! Cell lengths at previous time stamp
    real(8),dimension(2,n_max):: v           ! Velocity
    real(8),dimension(n_max):: omega         ! Angular velocity
    real(8),dimension(3*n_force_max,n_max):: force_n   ! Forces
    integer,dimension(n_max):: n_force       ! Number of forces on each cell
  
  end module bac_rod_data
  
  
  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------
  module bac_rod_func
    use bac_rod_data
  contains
  
    !------------------------------------------------------------------------
    subroutine init_parameter()
      ! Initialize run-time parameters.
      implicit none
      dt=dt_init
      i_error=0
      i_stop=0
      i_buckle=0
      i_measure_stress=0
      i_measure_force=0
      f_axis=0
      l_box_xt=l_box_x
      l_box_yt=l_box_y
    end subroutine init_parameter
  
    !------------------------------------------------------------------------
    subroutine initial()
      ! Initialize cell configuration.
      implicit none
      t=0
      if (init_cond==0) then
         call init_impo()                 ! load data from file
      elseif (init_cond==1) then
         call init_one()                  ! start with one cell
      elseif (init_cond==2) then
         call init_chain()                ! 1d array of cells
      elseif (init_cond==3) then
         call init_box_seed()                  ! arrange in a square box
      elseif (init_cond==4) then
         call init_box_dumb()                 ! let one cell grow to n_init cells
      elseif (init_cond==5) then
         call init_battle()
      endif
      if (t_init>0) then
         t=t_init                         ! Use t_init as the current time.
      end if
      if (i_stop==0) then
         write(*,*)"Initialization successfully!"
      end if
    end subroutine initial
  
    !------------------------------------------------------------------------
    subroutine init_impo()
      ! Import initial configuration from file.
      implicit none
      integer:: i, stat, impo_type
      impo_type=0   ! 0-import data as it is; 1-1d configuration; 2-import 2d to 3d.
      if (dimen_effe==1) then  ! Import 1D configuration
         impo_type=1
      elseif (dimen_effe==3.and.impo_2to3==1) then
         ! Import 2D configuration to 3D system.
         impo_type=2
      end if
      open(unit=1,file=trim(adjustl(file_impo)),status="old",action="read")
      i=1
      do
         if (impo_type==0) then
            read(1,*,iostat=stat)x(:,i),q(:,i),l(i),g(i),prop(i)
            prop(i)=prop_ave
         elseif (impo_type==1) then
            read(1,*,iostat=stat)x(1,i),l(i),g(i),prop(i)
         elseif (impo_type==2) then
            read(1,*,iostat=stat)x(1:2,i),q(1:2,i),l(i),g(i),prop(i)
            x(dimen_effe,i)=0.5_8*d_cell
            q(dimen_effe,i)=0
         end if
         if (stat/=0) exit          ! End of the file.
         i=i+1
      enddo
      n=i-1
      close(1)
  
    end subroutine init_impo
  
    !------------------------------------------------------------------------
    subroutine init_one()
      ! Start from one cell.
      implicit none
      integer:: i
      real(8):: ran1
      n=1
      x(1:dimen-1,1)=0
      x(dimen,1)=d_cell*0.5_8
      q(1,1)=1.0_8
      q(2:dimen,1)=0
      l(1)=l_max-0.02
      call random_number(ran1)
      g(1)=g_ave*(0.5_8+ran1)
      if (propel_on==1) then
         call random_number(ran1)
         prop(1)=prop_ave+prop_dev*2*(ran1-0.5)
      end if
    end subroutine init_one
  
    !------------------------------------------------------------------------
    subroutine init_chain()
      ! Arrange cells as a 1D chain along x.
      ! Overlap between neighboring cells is zero.
      implicit none
      integer:: i
      real(8):: ran1, lx
      ! Get distance between centers of cells at the two ends of the chain.
      n=n_init
      if (arra_leng==1) then       ! random cell length
         lx=(n-1)*(l_max+d_cell-0.5_8*l_dev)
      elseif (arra_leng==2.or.arra_leng==3) then   ! equal length
         lx=(n-1)*(l_max+d_cell)
      endif
      
      do i=1,n
         if (arra_leng==1) then         
            call random_number(ran1)
            l(i)=l_max-l_dev*ran1
         elseif (arra_leng==2) then
            l(i)=l_max
         elseif (arra_leng==3) then
            l(i)=l_max
            if (i==n_init/2) l(i)=l_cen
         endif
         x(2:dimen,i)=0
         if (i==1) then
            x(1,i)=-lx*0.5_8
         else
            x(1,i)=x(1,i-1)+0.5_8*(l(i-1)+l(i))+d_cell
         endif
         x(dimen,i)=d_cell*0.5_8
         q(2:dimen,i)=0
         q(1,i)=1.0_8
         call random_number(ran1)
         g(i)=g_ave*(0.5_8+ran1)
         if (propel_on==1) then
            call random_number(ran1)
            prop(i)=prop_ave+prop_dev*2*(ran1-0.5)
         end if
      enddo
    end subroutine init_chain
  
    !------------------------------------------------------------------------
    subroutine init_box_seed()
      ! Manually assign random positions (no overlaps) in the box to seed cells,
      !   and let the seed cells evolve (grow and move) to the expected configuration.
      ! If arra_leng=1, i.e., random cell length, number of seed cells is n_seed,
      !   and they have random lengths and grow rates.
      !   They will then grow (with g_time) until cell number is n_init.
      ! If arra_leng=2, i.e., fixed cell length, number of seed cells is n_init.
      !   They have minimum cell length and average growth rate g.
      !   They will elongate to the maximum cell length.
      implicit none
      real(8),dimension(3):: x_min1, x_min2
      real(8):: ran1, ran2, dr_ij, l_min, dr_min
      integer:: i, j, no_overlap, i_try, i_try_max
  
      ! Assign rods with minimum length and average growth rate,
      ! at random positions in the box, without overlapping with each other.
      i_try_max=10000
      l_min=0.5_8*(l_max-d_cell)
      if (arra_leng==1) then   ! Random length
         n=n_seed
      elseif (arra_leng==2) then   ! Equal length
         n=n_init
      end if
  
      ! Generating configurations of seed cells.
      ! Assign quantities except positions.
      q(dimen,1:n)=0
      do i=1,n
         call random_number(ran1)
         ran1=2*pi*ran1
         q(1,i)=cos(ran1)
         q(2,i)=sin(ran1)
         if (propel_on==1) then
            call random_number(ran1)
            prop(i)=prop_ave+prop_dev*2*(ran1-0.5)
         end if
         if (arra_leng==1) then
            call random_number(ran1)
            l(i)=l_min+(l_max-l_min)*ran1
            g(i)=g_ave+0.5_8*g_ave*(ran1-0.5_8)
         elseif (arra_leng==2) then
            l(i)=l_min
            g(i)=g_ave
         end if
      enddo
      ! Now assign positions.
      x(dimen,1:n)=0.5_8*d_cell
      call random_number(ran1)
      call random_number(ran2)
      x(1,1)=(l_box_xt-0.5*d_cell)*(ran1-0.5)
      x(2,1)=(l_box_yt-0.5*d_cell)*(ran2-0.5)
      do i=2,n
         no_overlap=0
         i_try=0
         do while (no_overlap==0)
            no_overlap=1
            i_try=i_try+1
            call random_number(ran1)
            call random_number(ran2)
            x(1,i)=(l_box_xt-0.5*d_cell)*(ran1-0.5)
            x(2,i)=(l_box_yt-0.5*d_cell)*(ran2-0.5)
            do j=1,i-1
               dr_ij=norm2(x(:,i)-x(:,j))
               if (dr_ij<d_cell+0.5*(l(i)+l(j))) then
                  call get_dr_min(x(:,i),x(:,j),q(:,i),q(:,j),l(i),l(j),x_min1,x_min2,dr_min)
                  if (dr_min<d_cell) then
                     no_overlap=0
                  end if
               end if
            end do
            if (i_try>=i_try_max) then
               no_overlap=1
               write(*,*)"Initialization failed! Too much trial!"
            end if
         end do
      end do
  
      ! Let the seed cells evolve to the expected initial condition.
      if (arra_leng==1) then
         do while (n<n_init)
            call run(3,n_init,propel_on,compress_on)
         end do
      elseif (arra_leng==2) then
         do while (l(1)<l_max)
            l(1:n)=l(1:n)+g(1:n)*dt
            call run(0,n_init,propel_on,compress_on)
         end do
      end if
  
    end subroutine init_box_seed
  
    !---------------------------------------------------------------------------
    subroutine init_box_dumb()
      ! Generate initial configuration with dumb cells from one cell.
      ! Step one:
      ! If n_dumb>0, it will let the single cell grow (with g_time) and divide to n_dumb.
      !    During this process, cells cannot propel, so the n_dumb cells finally forms an island at the box center.
      !    We then turn off the propulsion and growth of the first n_dumb cells, and create another cell,
      !    i.e., cell n_dumb+1 at the top-right corner of the box, which can grow and propel (if propel_on=1).
      ! If n_dumb=0, the system will generate one growing and propelling cell (if propel_on=1) at the box center.
      ! Step two:
      ! The system evolves until n=n_init.
      ! Step three:
      ! If cell lengths are expected to be equal, we let the system evolves,
      !    but cells with maximum length don't grow and divide.
      implicit none
      integer::i,n_short
      call init_one()
      i_stop=0   !Stop the simulation or not
      ! In case there are dummy cells.
      if (n_dumb>0) then
         do while (i_stop==0)
            t=t+dt
            call run(3,n_dumb,0,compress_on)
            if (n>=n_dumb.or.n>=n_max) then
               i_stop=1
            end if
            if (i_error==1) then
               write(*,*)'Some ERROR encounter!'
            endif
         enddo
         g(1:n)=0
         prop(1:n)=0
         i_stop=0
         ! Create another growing seed cell.
         n=n+1
         g(n)=g_ave
         if (propel_on==1) then
            prop(n)=prop_ave
         end if
         l(n)=l_max-0.1
         x(1,n)=l_box_xt/2.0-1.0
         x(2,n)=l_box_yt/2.0-1.0
         q(1,n)=1.0
         q(2,n)=0
      end if
      
      ! Continue growing until n=n_init.
      do while (i_stop==0)
         t=t+dt
         call run(3,n_init,propel_on,compress_on)
         if (n>=n_init.or.n>=n_max) then
            i_stop=1
         end if
         if (i_error==1) then
            write(*,*)'Some ERROR encounter!'
         endif
      enddo
  
      ! If cells have equal length
      i_stop=0
      if (arra_leng==2) then
         do while (i_stop==0)
            t=t+dt
            n_short=0
            do i=n_dumb+1,n
               if (l(i)<l_max) then
                  l(i)=l(i)+g(i)*dt*g_time
                  n_short=1
               end if
            end do
            call run(0,n_init,propel_on,compress_on)
            if (n_short==0) then
               i_stop=1
            end if
            if (i_error==1) then
               write(*,*)'Some ERROR encounter!'
            endif
         enddo
      end if
      t=0
    end subroutine init_box_dumb
  
    !------------------------------------------------------------------------
    subroutine init_battle()
      ! Battle formation with two groups of cells marching towards each other.
      ! Each group contains n_init/2 cells of identical length, and anti-parallel
      ! orientations.
      integer:: i,i1,i2
      real(8):: ran1, d_gap,xc
  
      n=n_init
      d_gap=1.0
      do i=1,n_left
         i1=n_left+1-i
         if (arra_leng==1) then         
            call random_number(ran1)
            l(i1)=l_max-l_dev*ran1
         elseif (arra_leng==2) then
            l(i1)=l_max
         endif
         x(2:dimen,i1)=0
         if (i==1) then
            x(1,i1)=-0.5_8*(d_gap+l(i1)+d_cell)
         else
            x(1,i1)=x(1,i1+1)-0.5_8*(l(i1+1)+l(i1))-d_cell
         endif
         x(dimen,i1)=d_cell*0.5_8
         q(2:dimen,i1)=0
         q(1,i1)=1.0_8
         call random_number(ran1)
         g(i1)=g_ave*(0.5_8+ran1)
         if (propel_on==1) then
            call random_number(ran1)
            prop(i1)=prop_ave+prop_dev*2*(ran1-0.5)
         end if
      enddo
  
      do i=1,n-n_left
         i2=n_left+i
         if (arra_leng==1) then         
            call random_number(ran1)
            l(i2)=l_max-l_dev*ran1
         elseif (arra_leng==2) then
            l(i2)=l_max
         endif
         x(2:dimen,i2)=0
         if (i==1) then
            x(1,i2)=0.5_8*(d_gap+l(i2)+d_cell)
         else
            x(1,i2)=x(1,i2-1)+0.5_8*(l(i2-1)+l(i2))+d_cell
         endif
         x(dimen,i2)=d_cell*0.5_8
         q(2:dimen,i2)=0
         q(1,i2)=-1.0_8
         call random_number(ran1)
         g(i2)=g_ave*(0.5_8+ran1)
         if (propel_on==1) then
            call random_number(ran1)
            prop(i2)=prop_ave+prop_dev*2*(ran1-0.5)
         end if
      enddo
  
      xc=sum(x(1,1:n))
      x(1,1:n)=x(1,1:n)-xc/n
      
    end subroutine init_battle
  
    !------------------------------------------------------------------------
    subroutine assign_box()
      ! Assign cells to mini-box, and identify the maximum(minimum) indices of non-empty boxes.
      implicit none
      integer::i,bx,by
      integer::p1
      ! Assign cells to box
      p_first=0
      p_next(1:n)=0
      bx_min=nx_box
      by_min=ny_box
      bx_max=0
      by_max=0
      do i=1,n
         if (mod(nx_box,2)==0) then
            bx=floor(x(1,i)/dl_box)+nx_box/2+1
         else
            bx=floor(x(1,i)/dl_box+0.5_8)+nx_box/2+1
         end if
         if (mod(ny_box,2)==0) then
            by=floor(x(2,i)/dl_box)+ny_box/2+1
         else
            by=floor(x(2,i)/dl_box+0.5_8)+ny_box/2+1
         end if
         if (bx<bx_min) then
            bx_min=bx
         endif
         if (bx>bx_max) then
            bx_max=bx
         endif
         if (by<by_min) then
            by_min=by
         endif
         if (by>by_max) then
            by_max=by
         endif
         p1=p_first(bx,by)
         p_first(bx,by)=i
         p_next(i)=p1
      enddo
  
      if (boundary_cond>=2.and.boundary_cond<=4) then    ! period boundary condition
         ! Assign boundaries to box
         if (mod(nx_box,2)==0) then
            bb_l=floor(-0.5_8*l_box_xt/dl_box)+nx_box/2+1
            bb_r=floor(0.5_8*l_box_xt/dl_box)+nx_box/2+1
         else
            bb_l=floor(-0.5_8*l_box_xt/dl_box+0.5_8)+nx_box/2+1
            bb_r=floor(0.5_8*l_box_xt/dl_box+0.5_8)+nx_box/2+1
         end if
         if (mod(ny_box,2)==0) then
            bb_b=floor(-0.5_8*l_box_yt/dl_box)+ny_box/2+1
            bb_t=floor(0.5_8*l_box_yt/dl_box)+ny_box/2+1
         else
            bb_b=floor(-0.5_8*l_box_yt/dl_box+0.5_8)+ny_box/2+1
            bb_t=floor(0.5_8*l_box_yt/dl_box+0.5_8)+ny_box/2+1
         end if
      
         if (bx_min<=bb_l+1.and.bx_max>=bb_r-1) then
            p_first(bx_min-2:bx_min-1,:)=p_first(bx_max-1:bx_max,:)
            p_first(bx_max+1:bx_max+2,:)=p_first(bx_min:bx_min+1,:)
         end if
         if (by_min<=bb_b+1.and.by_max>=bb_t-1) then
            p_first(:,by_min-2:by_min-1)=p_first(:,by_max-1:by_max)
            p_first(:,by_max+1:by_max+2)=p_first(:,by_min:by_min+1)
         end if
      end if
      ! write(*,*)bx_min,bx_max,by_min,by_max
      ! write(*,*)bb_l,bb_r,bb_b,bb_t
      
    end subroutine assign_box
  
    !------------------------------------------------------------------------
    subroutine assign_bead()
      ! Determine number and positions of beads in each cell.
      implicit none
      integer::p1,i,nb_i
      do p1=1,n
         nb_i=floor(l(p1)/dl_bead)+1
         if (nb_i==1) then
            xb(:,1,p1)=x(:,p1)
         else
            if (nb_i>nb_max) then
               nb_i=nb_max
            endif
            do i=1,nb_i/2
               xb(:,i,p1)=x(:,p1)+((i-1)*dl_bead-l(p1)*0.5)*q(:,p1)
               xb(:,nb_i-i+1,p1)=x(:,p1)-((i-1)*dl_bead-l(p1)*0.5)*q(:,p1)
            enddo
            if (mod(nb_i,2)==1) then
               xb(:,nb_i/2+1,p1)=xb(:,nb_i/2,p1)+dl_bead*q(:,p1)
            endif
         end if
         n_bead(p1)=nb_i
      end do
    end subroutine assign_bead
  
    !------------------------------------------------------------------------
     subroutine force_cell_cell()
      implicit none
      real(8),dimension(2,4)::i2
      integer::i,j,i1,j1,j1s,j2s,p1,p2,k
      integer::i10,i11,cont,p2_imag
      real(8),dimension(dimen)::dx_cc,x_p2_temp
      real(8),dimension(dimen,nb_max)::xb_p2_temp
      real(8)::l1h,l2h,dr_cc,nbc1,nbc2
      real(8)::rtemp1,rtemp2
  
      i2(:,1)=(/ 0, 1 /)
      i2(:,2)=(/ -1, 1 /)
      i2(:,3)=(/ -1, 1 /)
      i2(:,4)=(/ -1, 1 /)
      ! write(*,*)"start!"
      do j=by_min,by_max
         j1s=j+1
         if (j==by_max-1) then
            j1s=j+3
         elseif (j==by_max) then
            j1s=j+2
         end if
         do i=bx_min,bx_max
            p1=p_first(i,j)
            do while(p1>0)
               do j1=j,j1s
                  i10=i2(1,j1-j+1)
                  i11=i2(2,j1-j+1)
                  ! At right boundary
                  if (i==bx_max-1) then
                     i11=3
                  elseif (i==bx_max) then
                     i11=2
                  end if
                  ! At left boundary
                  if (j1>j) then  ! Only at higher rows 
                     if (i==bx_min+1) then
                        i10=-3
                     elseif (i==bx_min) then
                        i10=-2
                     end if
                  end if
                  do i1=i+i10,i+i11
                     
                     if (t>=0.5.and.i==6) then
                        !write(*,*)p1,j,i1,j1
                     end if
                     
                     if (j1==j.and.i1==i) then
                        p2=p_next(p1)
                     else
                        p2=p_first(i1,j1)
                     endif
                     do while(p2>0)
                        if (p1/=p2) then
                           ! p2 at the boundary of the box
                           if (i1<bx_min.or.i1>bx_max.or.j1>by_max) then
                              x_p2_temp=x(:,p2)
                              xb_p2_temp=xb(:,:,p2)
                              if (i1<bx_min) then
                                 x(1,p2)=x(1,p2)-l_box_xt
                                 xb(1,:,p2)=xb(1,:,p2)-l_box_xt
                              elseif (i1>bx_max) then
                                 x(1,p2)=x(1,p2)+l_box_xt
                                 xb(1,:,p2)=xb(1,:,p2)+l_box_xt
                              end if
                              if (j1>by_max) then
                                 x(2,p2)=x(2,p2)+l_box_yt
                                 xb(2,:,p2)=xb(2,:,p2)+l_box_yt
                              end if
                           end if
                           
                           l1h=l(p1)*0.5_8
                           l2h=l(p2)*0.5_8
                           dx_cc=x(:,p1)-x(:,p2)
                           dr_cc=norm2(dx_cc)
                           if (dr_cc<(l1h+l2h+r_cutoff)) then
                              rtemp1=abs(dot_product(dx_cc,q(:,p1)))
                              rtemp2=abs(dot_product(dx_cc,q(:,p2)))
                              dr_cc=dr_cc-(l1h*rtemp1+l2h*rtemp2)/dr_cc
                              if (dr_cc<r_cutoff) then
                                 if (rod_or_bead==1) then
                                    call force_rod_rod(p1,p2)
                                 elseif (rod_or_bead==2) then
                                    call force_bead_bead(p1,p2)
                                 elseif (rod_or_bead==3) then
                                    call force_rod_bead(p1,p2)
                                 end if
                              endif
                           endif
  
                           if (i1<bx_min.or.i1>bx_max.or.j1>by_max) then
                              x(:,p2)=x_p2_temp
                              xb(:,:,p2)=xb_p2_temp
                           end if
                        endif
                        p2=p_next(p2)
                     enddo
                  enddo
               enddo
               p1=p_next(p1)
            enddo
         enddo
      enddo
      ! write(*,*)"stop!"
      !write(*,*)k
    end subroutine force_cell_cell
    
    !------------------------------------------------------------------------
    subroutine force_rod_rod(p1,p2)
      implicit none
      integer,intent(in)::p1,p2
      real(8),dimension(dimen)::x_min1,x_min2,dx_min,f_2to1
      real(8),dimension(dimen)::dx_mc1,dx_mc2,vec_norm,vec_ave
      real(8)::dr_min,f_magnitude,norm_dx_min
  
      call get_dr_min(x(:,p1),x(:,p2),q(:,p1),q(:,p2),l(p1),l(p2),x_min1,x_min2,dr_min)
      if (dr_min<d_cell) then
         f_magnitude=e0*(d_cell-dr_min)**power_fcc
         !write(*,*)"contact",f_magnitude,dr_min
         if (dr_min>0.000000001_8) then
            dx_min=x_min1-x_min2
            norm_dx_min=norm2(dx_min)
            f_2to1=f_magnitude*dx_min/norm_dx_min
         else
            dx_mc1=x_min1-x(:,p1)
            dx_mc2=x_min2-x(:,p2)
            if (dimen==2) then
               vec_norm=cross_product_2d(dx_mc1,dx_mc2,1)
            elseif (dimen==3) then
               vec_norm=cross_product_3d(dx_mc1,dx_mc2)
            end if
            vec_ave=dx_mc1+dx_mc2
            if (dimen==2) then
               f_2to1=f_magnitude*cross_product_2d(vec_norm,vec_ave,2)
            elseif (dimen==3) then
               f_2to1=f_magnitude*cross_product_3d(vec_norm,vec_ave)
            end if
         endif
         call update_force(p1,x_min1,f_2to1,p2,x_min2)
      endif
    end subroutine force_rod_rod
  
    !------------------------------------------------------------------------
    subroutine force_bead_bead(p1,p2)
      implicit none
      integer::p1,p2,i,j
      real(8),dimension(dimen)::dx_2to1,f_2to1
      real(8)::dr_bb,f_magnitude,i_interact,h_overlap
  
      do i=1,n_bead(p1)
         do j=1,n_bead(p2)
            dx_2to1=xb(:,i,p1)-xb(:,j,p2)
            dr_bb=norm2(dx_2to1)
            h_overlap=d_cell-dr_bb
            if (h_overlap>0) then
               f_magnitude=e0*h_overlap**power_fcc
               f_2to1=f_magnitude*dx_2to1/dr_bb
               call update_force(p1,xb(:,i,p1),f_2to1,p2,xb(:,j,p2))
            elseif (h_overlap>-r_cohe) then
               f_magnitude=k_cohe*(4/r_cohe**2*(-h_overlap-r_cohe*0.5)**2-1)
               f_2to1=f_magnitude*dx_2to1/dr_bb
               call update_force(p1,xb(:,i,p1),f_2to1,p2,xb(:,j,p2))
            endif
         end do
      enddo
  
    end subroutine force_bead_bead
  
    !------------------------------------------------------------------------
    subroutine force_rod_bead(p1,p2)
      implicit none
      integer::p1,p2,i,j
      real(8),dimension(dimen)::dx_2to1
      real(8)::dr_bb,i_interact,h_overlap
      real(8),dimension(dimen)::x_min1,x_min2,dx_min,f_2to1
      real(8),dimension(dimen)::dx_mc1,dx_mc2,vec_norm,vec_ave
      real(8)::dr_min,f_magnitude,norm_dx_min
  
      call get_dr_min(x(:,p1),x(:,p2),q(:,p1),q(:,p2),l(p1),l(p2),x_min1,x_min2,dr_min)
      if (dr_min<d_cell) then
         f_magnitude=e0*(d_cell-dr_min)**power_fcc
         !write(*,*)"contact",f_magnitude,dr_min
         if (dr_min>0.000000001_8) then
            dx_min=x_min1-x_min2
            norm_dx_min=norm2(dx_min)
            f_2to1=f_magnitude*dx_min/norm_dx_min
         else
            dx_mc1=x_min1-x(:,p1)
            dx_mc2=x_min2-x(:,p2)
            if (dimen==2) then
               vec_norm=cross_product_2d(dx_mc1,dx_mc2,1)
            elseif (dimen==3) then
               vec_norm=cross_product_3d(dx_mc1,dx_mc2)
            end if
            vec_ave=dx_mc1+dx_mc2
            if (dimen==2) then
               f_2to1=f_magnitude*cross_product_2d(vec_norm,vec_ave,2)
            elseif (dimen==3) then
               f_2to1=f_magnitude*cross_product_3d(vec_norm,vec_ave)
            end if
         endif
         call update_force(p1,x_min1,f_2to1,p2,x_min2)
      endif
      
      do i=1,n_bead(p1)
         do j=1,n_bead(p2)
            dx_2to1=xb(:,i,p1)-xb(:,j,p2)
            dr_bb=norm2(dx_2to1)
            h_overlap=d_cell-dr_bb
            if (h_overlap<0.and.h_overlap>-r_cohe) then
               f_magnitude=k_cohe*(4/r_cohe**2*(-h_overlap-r_cohe*0.5)**2-1)
               f_2to1=f_magnitude*dx_2to1/dr_bb
               call update_force(p1,xb(:,i,p1),f_2to1,p2,xb(:,j,p2))
            endif
         end do
      enddo
  
    end subroutine force_rod_bead
  
    !------------------------------------------------------------------------
    subroutine force_random()
      implicit none
      real(8),dimension(dimen)::f_random,x_force
      real(8)::ran1,ran2,eta,f_magnitude
      integer::i
      eta=zeros
      do i=1,n
         call random_number(ran1)
         f_magnitude=eta*ran1
         call random_number(ran1)
         x_force=x(:,i)+(ran1-0.5_8)*q(:,i)
         call random_number(ran1)
         ran1=ran1*2*pi
         if (dimen==2) then
            f_random(1)=f_magnitude*cos(ran1)
            f_random(2)=f_magnitude*sin(ran1)
         elseif(dimen==3) then
            call random_number(ran2)
            ran2=ran2*pi
            f_random(1)=f_magnitude*sin(ran2)*cos(ran1)
            f_random(2)=f_magnitude*sin(ran2)*sin(ran1)
            f_random(dimen)=f_magnitude*cos(ran2)
         end if
         call update_force(i,x_force,f_random)
      end do
    end subroutine force_random
  
    !------------------------------------------------------------------------
    subroutine force_subs_rod1()
      implicit none
      real(8),dimension(dimen,3)::x_htc       ! coordinates of head, tail, and cell center
      real(8),dimension(dimen)::f_subs
      real(8)::h_overlap
      integer::i,i_low,j,i_high
      
      do i=1,n
         x_htc(:,1)=x(:,i)+0.5_8*l(i)*q(:,i)
         x_htc(:,2)=x(:,i)-0.5_8*l(i)*q(:,i)
         x_htc(:,3)=x(:,i)
         f_subs=0
         if (q(dimen,i)>0) then
            i_high=1
            i_low=2
         else
            i_high=2
            i_low=1
         endif
         
         ! Penetrate into the substrate
         h_overlap=0.5_8*d_cell-x_htc(dimen,i_low)
         if (h_overlap>0) then
            f_subs(dimen)=e_subs*h_overlap**power_fcs
            call update_force(i,x_htc(:,i_low),f_subs)
         end if
         
         ! Penetrate into the 'top agrose'
         f_subs(dimen)=-k_grav_head
         h_overlap=x_htc(dimen,i_high)-1.5_8*d_cell
         if (h_overlap>0) then
            f_subs(dimen)=f_subs(dimen)-e_agro*h_overlap**1.5_8
         end if
         call update_force(i,x_htc(:,i_high),f_subs)
         
         ! Adhesive force from substrate
         do j=1,3
            h_overlap=0.5_8*d_cell-x_htc(dimen,j)
            if (h_overlap<0.and.h_overlap>-r_adhe) then
               f_subs(dimen)=k_adhe*l(i)/3.0*(4/r_adhe**2*(-h_overlap-r_adhe*0.5)**2-1)
               call update_force(i,x_htc(:,j),f_subs)
            end if
         enddo
  
         ! constant gravity
         f_subs(dimen)=-k_grav_cent
         call update_force(i,x(:,i),f_subs)
         
      enddo
  
    end subroutine force_subs_rod1
  
    !------------------------------------------------------------------------
    subroutine force_subs_rod()
      implicit none
      real(8),dimension(dimen,3)::x_htc       ! coordinates of head, tail, and cell center
      real(8),dimension(dimen)::f_subs
      real(8)::h_overlap
      integer::i,i_low,j,i_high
      
      do i=1,n
         x_htc(:,1)=x(:,i)+0.5_8*l(i)*q(:,i)
         x_htc(:,2)=x(:,i)-0.5_8*l(i)*q(:,i)
         x_htc(:,3)=x(:,i)
         f_subs=0
         if (q(dimen,i)>0) then
            i_high=1
            i_low=2
         else
            i_high=2
            i_low=1
         endif
         
         !Penetrate into the substrate
         h_overlap=0.5_8*d_cell-x_htc(dimen,i_low)
         if (h_overlap>0) then
            f_subs(dimen)=e_subs*h_overlap**power_fcs
            call update_force(i,x_htc(:,i_low),f_subs)
         end if
  
         ! Adhesive force from substrate
         do j=1,2
            h_overlap=0.5_8*d_cell-x_htc(dimen,j)
            if (h_overlap<0.and.h_overlap>-r_adhe) then
               f_subs(dimen)=0.5_8*l(i)*k_adhe*h_overlap**power_adhe
               call update_force(i,x_htc(:,j),f_subs)
            elseif (h_overlap<-d_cell) then
               f_subs(dimen)=-e_agro*(-d_cell-h_overlap)
               call update_force(i,x_htc(:,j),f_subs)
            end if
         enddo
  
         ! Apply a torque proportional to 
         if (h_overlap<-d_cell*0.9) then
            f_subs(dimen)=-0.5*l(i)*prop(i)
            call update_force(i,x(:,i)+0.5_8*l(i)*q(:,i),f_subs)
         end if
         
      enddo
      
    end subroutine force_subs_rod
    
    !------------------------------------------------------------------------
    subroutine force_subs_bead()
      implicit none
      integer::p1,i
      real(8),dimension(dimen)::f_subs
      real(8)::h_overlap
  
      do p1=1,n
         do i=1,n_bead(p1)
            f_subs=0
            h_overlap=0.5_8*d_cell-xb(dimen,i,p1)
            if (h_overlap>0) then
               f_subs(dimen)=e_subs*h_overlap
            elseif (h_overlap<0.and.h_overlap>-r_adhe) then
               if (abs(q(dimen,p1))<0.3) then
                  f_subs(dimen)=k_adhe*l(i)/n_bead(p1)*(4/r_adhe**2*(-h_overlap-r_adhe*0.5)**2-1)
               endif
            endif
            call update_force(p1,xb(:,i,p1),f_subs)
         enddo
      enddo
  
    end subroutine force_subs_bead
  
    !------------------------------------------------------------------------
    subroutine force_boun_fixed()
      implicit none
      integer::p1,i
      real(8),dimension(dimen)::f_subs
      real(8)::h_overlap
  
  
    end subroutine force_boun_fixed
  
    !------------------------------------------------------------------------
    subroutine force_compress()
      implicit none
      real(8),dimension(dimen)::x_end,f_comp_vec
      real(8)::lih
      integer::i,i_minx,i_maxx
  
      f_comp_vec=0
      if (vari_fc==0) then
         f_comp_vec(1)=f_compress
      else
         f_comp_vec(1)=f_compress*t/t_max
      end if
         
      i_minx=1
      i_maxx=1
      do i=2,n
         if (x(1,i)<x(1,i_minx)) then
            i_minx=i
         endif
         if (x(1,i)>x(1,i_maxx)) then
            i_maxx=i
         endif
      enddo
  
      i=i_minx
      lih=l(i)*0.5_8
      if (q(1,i)>=0) then
         x_end=x(:,i)-lih*q(:,i)
      else
         x_end=x(:,i)+lih*q(:,i)
      endif
      call update_force(i,x_end,f_comp_vec)
  
      i=i_maxx
      lih=l(i)*0.5_8
      if (q(1,i)>=0) then
         x_end=x(:,i)+lih*q(:,i)
      else
         x_end=x(:,i)-lih*q(:,i)
      endif
      call update_force(i,x_end,-f_comp_vec)
  
    end subroutine force_compress
  
    !------------------------------------------------------------------------
    subroutine force_propel()
      implicit none
      integer::i
      do i=1,n
         call update_force(i,x(:,i),prop(i)*l(i)*q(:,i))
      end do
    end subroutine force_propel
  
    !------------------------------------------------------------------------
    subroutine update_force(p1,x_min1,f_vec,p2,x_min2)
      ! still need f_axis
      implicit none
      integer:: i_temp
      integer,intent(in)::p1
      real(8),dimension(dimen)::x_min1,f_vec,x_contact
      real(8),dimension(dimen)::f_n,x_cc
      real(8):: dr_cc
      integer,intent(in),optional::p2
      real(8),dimension(dimen),optional::x_min2
      dx(:,p1)=dx(:,p1)+f_vec
      if (dimen==2) then
         dq(:,p1)=dq(:,p1)+cross_product_2d(x_min1-x(:,p1),f_vec,1)
      elseif (dimen==3) then
         dq(:,p1)=dq(:,p1)+cross_product_3d(x_min1-x(:,p1),f_vec)
      end if
      if (grow_type==2) then
         x_cc(1:dimen_effe)=x_min1(1:dimen_effe)-x(1:dimen_effe,p1)
         dr_cc=norm2(x_cc)
         if (dr_cc>zeros) then
            f_axis(p1)=f_axis(p1)-0.5_8*dot_product(x_cc,f_vec)/dr_cc
         end if
      end if
      if (present(p2)) then
         dx(:,p2)=dx(:,p2)-f_vec
         if (dimen==2) then
            dq(:,p2)=dq(:,p2)+cross_product_2d(x_min2-x(:,p2),-f_vec,1)
         elseif (dimen==3) then
            dq(:,p2)=dq(:,p2)+cross_product_3d(x_min2-x(:,p2),-f_vec)
         end if
         if (grow_type==2) then
            x_cc(1:dimen_effe)=x_min2(1:dimen_effe)-x(1:dimen_effe,p2)
            dr_cc=norm2(x_cc)
            if (dr_cc>zeros) then
               f_axis(p2)=f_axis(p2)+0.5_8*dot_product(x_cc,f_vec)/dr_cc
            end if
         end if
         ! Measure stress
         if (i_measure_stress==1) then
            x_contact=(x_min1+x_min2)/2
            ! For cell p1
            if (stress_xy_n==1) then
               ! Measure stress in x-y coordinate system
               x_cc(1:dimen_effe)=x_contact(1:dimen_effe)-x(1:dimen_effe,p1)
               f_n=f_vec
            elseif (stress_xy_n==2) then
               ! Get stress in n-np coordinate.
               x_cc(1)=q(1,p1)*(x_contact(1)-x(1,p1))+q(2,p1)*(x_contact(2)-x(2,p1))
               x_cc(2)=-q(2,p1)*(x_contact(1)-x(1,p1))+q(1,p1)*(x_contact(2)-x(2,p1))
               f_n(1)=q(1,p1)*f_vec(1)+q(2,p1)*f_vec(2)
               f_n(2)=-q(2,p1)*f_vec(1)+q(1,p1)*f_vec(2)
            end if
            stress_n_t(1:2,p1)=stress_n_t(1:2,p1)+x_cc(1)*f_n(1:2)
            stress_n_t(3:4,p1)=stress_n_t(3:4,p1)+x_cc(2)*f_n(1:2)
            ! Cell p2
            if (stress_xy_n==1) then
               ! Measure stress in x-y coordinate system
               x_cc(1:dimen_effe)=x_contact(1:dimen_effe)-x(1:dimen_effe,p2)
               f_n=-f_vec
            elseif (stress_xy_n==2) then
               ! Get stress in n-np coordinate.
               x_cc(1)=q(1,p2)*(x_contact(1)-x(1,p2))+q(2,p2)*(x_contact(2)-x(2,p2))
               x_cc(2)=-q(2,p2)*(x_contact(1)-x(1,p2))+q(1,p2)*(x_contact(2)-x(2,p2))
               f_n(1)=-q(1,p2)*f_vec(1)-q(2,p2)*f_vec(2)
               f_n(2)=q(2,p2)*f_vec(1)-q(1,p2)*f_vec(2)
            end if
            stress_n_t(1:2,p2)=stress_n_t(1:2,p2)+x_cc(1)*f_n(1:2)
            stress_n_t(3:4,p2)=stress_n_t(3:4,p2)+x_cc(2)*f_n(1:2)
         end if
         ! Measure force
         if (i_measure_force==1) then
            ! Cell p1
            i_temp=n_force(p1)*3
            if (abs(q(1,p1)) > abs(q(2,p1))) then
               force_n(i_temp+1,p1)=(x_min1(1)-x(1,p1))/(l(p1)*q(1,p1))
            else
               force_n(i_temp+1,p1)=(x_min1(2)-x(2,p1))/(l(p1)*q(2,p1))
            end if
            force_n(i_temp+2,p1)=f_vec(1)
            force_n(i_temp+3,p1)=f_vec(2)
            ! Cell p2
            i_temp=n_force(p2)*3
            if (abs(q(1,p2)) > abs(q(2,p2))) then
               force_n(i_temp+1,p2)=(x_min2(1)-x(1,p2))/(l(p2)*q(1,p2))
            else
               force_n(i_temp+1,p2)=(x_min2(2)-x(2,p2))/(l(p2)*q(2,p2))
            end if
            force_n(i_temp+2,p2)=-f_vec(1)
            force_n(i_temp+3,p2)=-f_vec(2)
            if (n_force(p1)<n_force_max) then
               n_force(p1)=n_force(p1)+1
            end if
            if (n_force(p2)<n_force_max) then
               n_force(p2)=n_force(p2)+1
            end if
         end if
      end if
    end subroutine update_force
    
    !------------------------------------------------------------------------
    subroutine get_dr_min(xc1,xc2,q1,q2,l1,l2,x_min1,x_min2,dr_min)
      implicit none
      real(8),intent(in),dimension(dimen)::xc1,xc2,q1,q2
      real(8),intent(out),dimension(dimen)::x_min1,x_min2
      real(8),dimension(dimen)::x_min10,x_min20
      real(8),intent(in)::l1,l2
      real(8),intent(out)::dr_min
      integer::i
      real(8)::a,b1,b2,lambda_l1,lambda_l2,l1h,l2h,dr_min0
      real(8),dimension(4)::dr_mine,lambda_e1,lambda_e2
  
      l1h=0.5_8*l1
      l2h=0.5_8*l2
      a=dot_product(q1,q2)
      b1=dot_product(q1,xc1-xc2)
      b2=dot_product(q2,xc1-xc2)
      lambda_l1=(-a*b2+b1)/(a**2-1)
      lambda_l2=(a*b1-b2)/(a**2-1)
      if (abs(lambda_l1)<l1h.and.abs(lambda_l2)<l2h) then
         ! Two minimum points within the two segments.
         x_min1=xc1+lambda_l1*q1
         x_min2=xc2+lambda_l2*q2
         dr_min=norm2(x_min1-x_min2)
      else
         lambda_e1=[l1h,-l1h,-b1+a*l2h,-b1-a*l2h]
         lambda_e2=[b2+a*l1h,b2-a*l1h,l2h,-l2h]
         where (lambda_e1<-l1h)
            lambda_e1=-l1h
         elsewhere (lambda_e1>l1h)
            lambda_e1=l1h
         end where
         where (lambda_e2<-l2h)
            lambda_e2=-l2h
         elsewhere (lambda_e2>l2h)
            lambda_e2=l2h
         end where
         dr_min=100*d_cell
         do i=1,4
            x_min10=xc1+lambda_e1(i)*q1
            x_min20=xc2+lambda_e2(i)*q2
            dr_min0=norm2(x_min10-x_min20)
            if (dr_min0<dr_min) then
               dr_min=dr_min0
               x_min1=x_min10
               x_min2=x_min20
            endif
         end do
      end if
  
    end subroutine get_dr_min
  
    !------------------------------------------------------------------------
    subroutine grow(grow_type_t,n_max_t)
      ! grow_type_t controls the how bacteria are growing and dividing:
      ! 1 - standard growth and division
      ! 2 - cell division depends on local force to eliminate sudden change of local
      !     force due to cell division. This is done by add an overlap between the
      !     two daughter cells, creating an axial force to balance the external force.
      ! 3 - cell growth is changed by a factor of g_time
      ! n_max_t controls the maximum number of cells. Cell growth and division will be
      !     inhibited if n>=n_max_t
      implicit none
      integer::grow_type_t,n_max_t,n_temp,i,k,km,i_divide
      real(8)::ran1,ran2,dx_overlap,dx_overlap0,kx,dk,l_temp,l_temp0
      integer,dimension(8)::kmin
      
      if (n>=n_max_t) return
      n_temp=n
      dx_overlap=0
      kmin=(/ 0, 1, 4, 10, 20, 35, 56, 84 /)
      
      do i=1,n
  
         ! Cell growth.
         if (grow_type_t==3) then
            l(i)=l(i)+g(i)*dt*g_time
         else
            l(i)=l(i)+g(i)*dt
         end if
  
         ! Check cell division.
         i_divide=0
         if (n_temp<n_max_t) then
            if (grow_type_t==2) then
               if (rod_or_bead==1) then
                  dx_overlap=f_axis(i)/e0
                  ! f_axis times 0.5 because we count for both ends
               elseif (rod_or_bead==2) then
                  dx_overlap0=f_axis(i)/e0
                  kx=dx_overlap0/dl_bead
                  do km=1,7
                     if (kx>=kmin(km).and.kx<kmin(km+1)) then
                        dk=kx-kmin(km)
                        dx_overlap=(km-1)*dl_bead+2*dl_bead*dk/(km*(km+1))
                     endif
                  enddo
               endif
               if (l(i)+dx_overlap>l_max) then
                  ! Cell divides when the daughter cells have length 0.5*(l_max-d_cell).
                  i_divide=1
               endif
            else
               if (l(i)>=l_max) then
                  i_divide=1
               endif
            endif
         end if
  
         ! Cell division.
         if (i_divide==1) then
            n_temp=n_temp+1
            q(:,n_temp)=q(:,i)
            l(i)=0.5_8*(l(i)-d_cell+dx_overlap)
            l(n_temp)=l(i)
            l_temp=0.5_8*(l(i)+d_cell-dx_overlap)
            x(:,n_temp)=x(:,i)+l_temp*q(:,i)
            x(:,i)=x(:,i)-l_temp*q(:,i)
            ! Generate new growth rates and propulsions
            call random_number(ran1)
            call random_number(ran2)
            g(i)=g_ave*(0.5_8+ran1)
            g(n_temp)=g_ave*(0.5_8+ran2)
            call random_number(ran1)
            call random_number(ran2)
            prop(i)=prop_ave+prop_dev*2*(ran1-0.5)
            prop(n_temp)=prop_ave+prop_dev*2*(ran2-0.5)
            ! Estimate the configurations for the previous time stamp.
            q0(:,n_temp)=q0(:,i)
            l0(i)=0.5_8*(l0(i)-d_cell+dx_overlap)
            l0(n_temp)=l0(i)
            l_temp0=0.5_8*(l0(i)+d_cell-dx_overlap)
            x0(:,n_temp)=x0(:,i)+l_temp0*q0(:,i)
            x0(:,i)=x0(:,i)-l_temp0*q0(:,i)
            v(:,n_temp)=v(:,i)
            omega(n_temp)=omega(i)
            stress_n(:,n_temp)=stress_n(:,i)
         endif
      enddo
      n=n_temp
    end subroutine grow
  
    !------------------------------------------------------------------------
    subroutine rescale()
      implicit none
      integer::i
      real(8)::k_lx,k_ly
      k_lx=1+v_box_x*dt/l_box_xt
      K_ly=1+v_box_y*dt/l_box_yt
      x(1,1:n)=k_lx*x(1,1:n)
      x(2,1:n)=k_ly*x(2,1:n)
      l_box_xt=k_lx*l_box_xt
      l_box_yt=k_ly*l_box_yt
    end subroutine rescale
  
    !------------------------------------------------------------------------
    subroutine boun_fixed(dire)
      implicit none
      integer::i,dire,go_x,go_y,i_high,i_low,i_left,i_right
      real(8)::lbhx,lbhy,h_overlap
      real(8),dimension(dimen,2)::x_ht
      real(8),dimension(dimen)::f_subs
      lbhx=l_box_xt*0.5_8
      lbhy=l_box_yt*0.5_8
      go_x=0
      go_y=0
      if (dire==1) then
         go_x=1
      elseif (dire==2) then
         go_y=1
      elseif (dire==3) then
         go_x=1
         go_y=1
      end if
      do i=1,n
         x_ht(:,1)=x(:,i)+0.5_8*l(i)*q(:,i)
         x_ht(:,2)=x(:,i)-0.5_8*l(i)*q(:,i)
         if (q(dimen,i)>0) then
            i_high=1
            i_low=2
         else
            i_high=2
            i_low=1
         endif
         if (q(1,i)>0) then
            i_right=1
            i_left=2
         else
            i_right=2
            i_left=1
         endif
  
         if (go_x==1) then   ! Check x boundaries
            f_subs=0    ! Left boundary
            h_overlap=-lbhx+0.5_8*d_cell-x_ht(1,i_left)
            if (h_overlap>0) then
               f_subs(1)=e0*h_overlap**power_fcc
               call update_force(i,x_ht(:,i_left),f_subs)
            end if
            f_subs=0    ! Right boundary
            h_overlap=x_ht(1,i_right)+0.5_8*d_cell-lbhx
            if (h_overlap>0) then
               f_subs(1)=-e0*h_overlap**power_fcc
               call update_force(i,x_ht(:,i_right),f_subs)
            end if
         end if
         if (go_y==1) then   ! Check y boundaries
            f_subs=0    ! Low boundary
            h_overlap=-lbhy+0.5_8*d_cell-x_ht(dimen,i_low)
            if (h_overlap>0) then
               f_subs(dimen)=e0*h_overlap**power_fcc
               call update_force(i,x_ht(:,i_low),f_subs)
            end if
            f_subs=0    ! High boundary
            h_overlap=x_ht(dimen,i_high)+0.5_8*d_cell-lbhy
            if (h_overlap>0) then
               f_subs(dimen)=-e0*h_overlap**power_fcc
               call update_force(i,x_ht(:,i_high),f_subs)
            end if
         end if
            
      enddo
    end subroutine boun_fixed
  
    !------------------------------------------------------------------------
    subroutine boun_periodic(dire)
      implicit none
      integer::i,dire,go_x,go_y
      real(8)::lbhx,lbhy
      lbhx=l_box_xt*0.5_8
      lbhy=l_box_yt*0.5_8
      go_x=0
      go_y=0
      if (dire==1) then
         go_x=1
      elseif (dire==2) then
         go_y=1
      elseif (dire==3) then
         go_x=1
         go_y=1
      end if
      do i=1,n
         if (go_x==1) then
            if (x(1,i)<-lbhx) then
               x(1,i)=x(1,i)+l_box_xt
            elseif (x(1,i)>lbhx) then
               x(1,i)=x(1,i)-l_box_xt
            end if
         end if
         if (go_y==1) then
            if (x(2,i)<-lbhy) then
               x(2,i)=x(2,i)+l_box_yt
            elseif (x(2,i)>lbhy) then
               x(2,i)=x(2,i)-l_box_yt
            end if
         end if
      enddo
    end subroutine boun_periodic
  
    !------------------------------------------------------------------------
    subroutine boun_absorb(dire)
      implicit none
      integer::i,n_temp,dire,go_x,go_y
      real(8)::lbhx,lbhy
      lbhx=l_box_xt*0.5
      lbhy=l_box_yt*0.5
      n_temp=n
      go_x=0
      go_y=0
      if (dire==1) then
         go_x=1
      elseif (dire==2) then
         go_y=1
      elseif (dire==3) then
         go_x=1
         go_y=1
      end if
      do i=n_temp,1,-1
         if (go_x==1) then
            ! Absorbed boundary at lx/2, -lx/2
            if (abs(x(1,i))>lbhx-0.5*l_max) then
               call remove_cell_i(i)
            end if
         end if
         if (go_y==1) then
            if (abs(x(2,i))>lbhy-0.5*l_max) then
               call remove_cell_i(i)
            end if
         end if
      enddo
    end subroutine boun_absorb
  
    !------------------------------------------------------------------------
    subroutine remove_cell_i(i)
      ! Remove i-th cell from the system
      implicit none
      integer,intent(in)::i
      if (i<n) then
         x(:,i:n-1)=x(:,i+1:n)
         q(:,i:n-1)=q(:,i+1:n)
         l(i:n-1)=l(i+1:n)
         g(i:n-1)=g(i+1:n)
         prop(i:n-1)=prop(i+1:n)
         x0(:,i:n-1)=x0(:,i+1:n)
         q0(:,i:n-1)=q0(:,i+1:n)
         l0(i:n-1)=l0(i+1:n)
         if (rod_or_bead>=2) then
            n_bead(i:n-1)=n_bead(i+1:n)
            xb(:,:,i:n-1)=xb(:,:,i+1:n)
         end if
      end if
      n=n-1
    end subroutine remove_cell_i
  
    !------------------------------------------------------------------------
    subroutine get_dt()
      ! Get optimal dt when adaptive dt is on.
      implicit none
      integer::i
      real(8)::dx_max,dri,mv,ma,dx_h,dx_t,dx_maxi
      real(8),dimension(dimen)::dx_rot
  
      dx_max=0.1_8**7
      do i=1,n
         mv=1/l(i)
         ma=12.0_8/l(i)**3
         ! displacement due to rotation
         if (dimen==2) then
            dx_rot=0.5_8*l(i)*ma*cross_product_2d(dq(:,i),q(:,i),2) 
         elseif (dimen==3) then
            dx_rot=0.5_8*l(i)*ma*cross_product_3d(dq(:,i),q(:,i))
         end if
         ! total displacement
         dx_h=norm2(mv*dx(:,i)+dx_rot)
         dx_t=norm2(mv*dx(:,i)-dx_rot)
         dx_maxi=maxval([dx_h,dx_t])
         if (dx_maxi>dx_max) then
            dx_max=dx_maxi
         endif
      enddo
      dt=dx_maxc/dx_max
      if (dt>dt_max) then
         dt=dt_max
      elseif (dt<dt_min.and.exp_conf==1) then
         !write(*,*)'Warning! dt is smaller than dt_min!'
      endif
    end subroutine get_dt
  
    !------------------------------------------------------------------------
    subroutine measure_core(t_measure_t)
      ! Measure stress without considering local density
      implicit none
      integer::i,ti,t_measure_t
      real(8),dimension(dimen):: dqi
      real(8),dimension(dimen,n_max)::x_temp,q_temp
      real(8):: v_corr_x,v_corr_y
  
      ! Measure velocity and angular velocity
      v(1:2,1:n)=(x(1:2,1:n)-x0(1:2,1:n))/dt_exp
      v_corr_y=l_box_yt/dt_exp     ! Correct the velocity in case of periodic hopping to the other side of the system.
      do i=1,n
         if (v(2,i)<-0.8*v_corr_y) then
            v(2,i)=v(2,i)+v_corr_y
         elseif (v(2,i)>0.8*v_corr_y) then
            v(2,i)=v(2,i)-v_corr_y
         end if
         dqi=q(:,i)-q0(:,i)
         omega(i)=asin(q0(1,i)*dqi(2)-q0(2,i)*dqi(1))/dt_exp
      end do
      x0=x
      q0=q
      l0=l
      ! Measure stress
      i_measure_stress=1
      x_temp(:,1:n)=x(:,1:n)
      q_temp(:,1:n)=q(:,1:n)
      stress_n(:,1:n)=0
      do ti=1,t_measure_t
         call run(0,n_max,propel_on,compress_on)
         stress_n(:,1:n)=stress_n(:,1:n)+stress_n_t(:,1:n)
      end do
      do i=1,4
         stress_n(i,1:n)=stress_n(i,1:n)/ti/(2*(l+d_cell)*d_cell)
      end do
      x(:,1:n)=x_temp(:,1:n)
      q(:,1:n)=q_temp(:,1:n)
      i_measure_stress=0
      
    end subroutine measure_core
  
    !------------------------------------------------------------------------
    subroutine get_force(t_measure_t)
      ! Measure force without considering local density
      implicit none
      integer::i,ti,t_measure_t
      real(8),dimension(dimen,n_max)::x_temp,q_temp
      i_measure_force=1
      x_temp(:,1:n)=x(:,1:n)
      q_temp(:,1:n)=q(:,1:n)
      n_force(1:n)=0
      force_n(:,1:n)=0
      do ti=1,t_measure_t
         call run(0,n_max,propel_on,compress_on)
      end do
      x(:,1:n)=x_temp(:,1:n)
      q(:,1:n)=q_temp(:,1:n)
      i_measure_force=0
    end subroutine get_force
  
    !------------------------------------------------------------------------
    subroutine check_buckle()
      implicit none
      integer:: i
      real(8):: max_h
      character(len=80)::s_label
      i_buckle=0
      max_h=0
      do i=1,n
         if (x(dimen,i)>0.51_8*d_cell) then
            i_buckle=1
            if (x(dimen,i)>max_h) then
               max_h=x(dimen,i)
               p_buckle=i
            end if
         endif
      enddo
      if (i_buckle==1) then
         call export_results()
         write(s_label,*)"final"
         call export_conf(s_label)
         call export_meas(s_label)
      end if
    end subroutine check_buckle
  
    !------------------------------------------------------------------------
    subroutine export(s_label)
      implicit none
      character(len=80)::s_label
      if (exp_conf==1) then
         call export_conf(s_label)
      endif
      if (exp_meas==1) then
         call export_meas(s_label)
      endif
    end subroutine export
    
    !------------------------------------------------------------------------
    subroutine export_conf(s_label)
      implicit none
      integer::i,j
      character(len=80),intent(in)::s_label
      character(len=80)::file_conf
      write(file_conf,*)trim(adjustl(dire_expo))//"conf_"//trim(adjustl(s_label))//".dat"
      open(unit=10,file=adjustl(file_conf),status="unknown",action="write")
      do i=1,n
         write(10,*)x(1:dimen_effe,i),q(1:dimen_effe,i),l(i),g(i),prop(i)
      enddo
      close(10)
    end subroutine export_conf
  
    !------------------------------------------------------------------------
    subroutine export_meas(s_label)
      implicit none
      integer::i,j
      character(len=80),intent(in)::s_label
      character(len=80)::file_meas
      write(file_meas,*)trim(adjustl(dire_expo))//"meas_"//trim(adjustl(s_label))//".dat"
      open(unit=11,file=adjustl(file_meas),status="unknown",action="write")
      do i=1,n
         if (dimen_effe==1) then
            write(11,*)x(1,i),stress_n(1,i),f_axis(i)
         else
            if (force_n(1,i)<1) then
               !write(*,*) 'cell', i, force_n(1:4,i)
            end if
            if (meas_stress==1.and.meas_force==1) then
               write(11,*)stress_n(1:4,i),v(1:2,i),omega(i),n_force(i),force_n(1:n_force_max*3,i)
            elseif (meas_stress==1.and.meas_force==0) then
               write(11,*)stress_n(1:4,i),v(1:2,i),omega(i)
            elseif (meas_stress==0.and.meas_force==1) then
               write(11,*)force_n(:,i)
            end if
         end if
      enddo
      close(11)
    end subroutine export_meas
  
    !------------------------------------------------------------------------
    subroutine export_results()
      implicit none
      integer::i
      real(8)::r_colony,xc_colony,yc_colony
      character(len=80)::file_results
      write(file_results,*)trim(adjustl(dire_expo))//"results.dat"
      if (dimen_effe==2) then
         r_colony=0.5_8*( maxval(x(1,1:n)) - minval(x(1,1:n)) )
         xc_colony=0.5_8*( maxval(x(1,1:n)) + minval(x(1,1:n)) )
      elseif (dimen_effe==3) then
         r_colony=0.25_8*( maxval(x(1,1:n)) - minval(x(1,1:n)) )
         r_colony=r_colony+0.25_8*( maxval(x(2,1:n)) - minval(x(2,1:n)) )
         xc_colony=0.5_8*( maxval(x(1,1:n)) + minval(x(1,1:n)) )
         yc_colony=0.5_8*( maxval(x(2,1:n)) + minval(x(2,1:n)) )
      end if
      open(unit=12,file=adjustl(file_results),status="unknown",action="write")
      i=p_buckle
      write(12,*)t,n,r_colony,x(1,i)-xc_colony,l(i),g(i),stress_n(1,i)
      close(12)
    end subroutine export_results
    
    !------------------------------------------------------------------------
    subroutine init_random_seed()
      integer:: i, n_random_seed, clock
      integer,dimension(:),allocatable:: seed
      call random_seed(size = n_random_seed)
      allocate(seed(n_random_seed))
      call system_clock(count=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n_random_seed) /)
      call random_seed(put = seed)
      deallocate(seed)
    end subroutine init_random_seed
  
    !------------------------------------------------------------------------
    function cross_product_2d(a,b,i_mode)
      real(8),dimension(dimen)::cross_product_2d
      real(8),dimension(dimen),intent(in)::a,b
      integer,intent(in)::i_mode
      if (i_mode==1) then  ! x*y->z
         cross_product_2d(1)=a(1)*b(2)-a(2)*b(1)
      elseif (i_mode==2) then  ! z*x->y
         cross_product_2d(1)=-a(1)*b(2)
         cross_product_2d(2)=a(1)*b(1)
      end if
    end function cross_product_2d
    
    !------------------------------------------------------------------------
    function cross_product_3d(a,b)
      real(8),dimension(dimen):: cross_product_3d
      real(8),dimension(dimen),intent(in):: a,b
      cross_product_3d(1)= a(2)*b(dimen)-a(dimen)*b(2)
      cross_product_3d(2)= a(dimen)*b(1)-a(1)*b(dimen)
      cross_product_3d(dimen)= a(1)*b(2)-a(2)*b(1)
    end function cross_product_3d
  
    !-------------------------------------------------------------------------
    function vec_uni(a)
      real(8),dimension(3),intent(in):: a
      real(8),dimension(3)::vec_uni
      real(8)::dr
      dr=sqrt(a(1)**2+a(2)**2+a(3)**2)
      if (dr/=0) then
         vec_uni=a/dr
      else
         vec_uni=0
      endif
    end function vec_uni
  
    !-------------------------------------------------------------------------
    
  end module bac_rod_func
  
  
  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------
  program bac_rod
    use bac_rod_func
    implicit none
    real::cput0,cput1
    character(len=80)::s_vari
    real(8),dimension(dimen)::xt1,xt2,qt1,qt2,xm1,xm2
    real(8),dimension(dimen,10,100)::xbead
    real(8),dimension(dimen,100)::xcell
    real(8)::l1,l2,b,dr
    integer::i,j,k
  
    call cpu_time(cput0)
    if (random_seed_time==0) then
       call random_seed()
    else
       call init_random_seed()
    end if
    call init_parameter()
    call evolve()
    call cpu_time(cput1)
    write(s_vari,format_real_screen)cput1-cput0
    write(*,*)"Program finished, the time spent is:"//trim(adjustl(s_vari))//" s."
    
  end program bac_rod
  
  !---------------------------------------------------------------------------
  subroutine evolve()
    use bac_rod_func
    implicit none
    integer::i,i_expo,n_step,n_temp,n_divi_temp,n_while
    integer::t_exp_0,t_exp
    real::dt_ave,t_sum,dx_col,dy_col,dr_col,dr_max
    character(len=80)::s_t,s_n,s_label,command
    
    ! Generate initial configurations.
    call initial()
    ! Delete old data if demanded.
    if (dele_data==1) then
       if (exp_conf==1) then
          write(command,*)"rm "//trim(adjustl(dire_expo))//"conf_*"
          call system(command)
       endif
       if (exp_meas==1) then
          write(command,*)"rm "//trim(adjustl(dire_expo))//"meas_*"
          call system(command)
       endif
    endif
  
    n_while=0  !Number of while iterations.
    i_stop=0   !Stop the simulation or not
    n_step=0   !Time steps between two exports.
    t_sum=0
    n_divi=0
    t_exp_0=-1
    
    do while (i_stop==0)
       t=t+dt
       i_expo=0
       t_sum=t_sum+dt
       n_step=n_step+1
  
       ! Check to export or not.
       if (timing_scheme==1) then
          if (t>=t0) then
             t_exp=floor(t/dt_exp)
             if (t_exp>t_exp_0) then
                i_expo=1
                write(s_label,*)t_exp
                t_exp_0=t_exp
                write(s_t,*)floor(100*t/(t_max+0.0_8))
             end if
          endif
       elseif (timing_scheme==2) then
          if (n>=n0) then
             t_exp=n/dn_exp
             if (t_exp>t_exp_0) then
                i_expo=1
                t_exp_0=t_exp
                write(s_label,*)n
                write(s_t,*)floor(100*n/(n_max+0.0_8))
             end if
          endif
       elseif (timing_scheme==3) then
          if (n_divi>=n0_divi) then
             if (n_divi/dn_exp>n_divi_temp/dn_exp) then
                i_expo=1
                write(s_label,*)n_divi/dn_exp
                write(s_t,*)floor(100*n_divi/(n_max_divi+0.0_8))
             end if
          endif
       elseif (timing_scheme>=4.and.timing_scheme<=5) then
          dx_col=maxval(x(1,1:n))-minval(x(1,1:n))
          dy_col=maxval(x(2,1:n))-minval(x(2,1:n))
          if (timing_scheme==4) then  ! Expansion in x.
             dr_col=dx_col
             dr_max=l_box_x
          elseif (timing_scheme==5) then  ! Expansion in y.
             dr_col=dy_col
             dr_max=l_box_y
          elseif (timing_scheme==6) then  ! Expansion in both direction.
             dr_col=0.5*(dx_col+dy_col)
             dr_max=min(l_box_x,l_box_y)
          end if
          if (dr_col>=y_c0) then
             t_exp=floor(dr_col/dl_exp)
             if (dr_col>=dr_max-d_cell) then  ! If gets the boundary.
                if (t_exp==t_exp_0) then
                   t_exp=t_exp+1
                end if
                if (stop_lmax==1) then
                   i_stop=1
                   write(*,*)'Maximum colony size reached.'
                end if
             end if
             if (t_exp>t_exp_0) then
                i_expo=1
                t_exp_0=t_exp
                write(s_label,*)t_exp
                write(s_t,*)floor(100*t/(t_max+0.0_8))
             end if
          endif
       
       endif
       
       if (i_check_buckle==1) then
          call check_buckle()
          if (i_buckle==1) then
             write(*,*)'Buckling detected.'
             if (stop_buckle==1) then
                i_stop=1
             endif
          end if
       end if
       
       ! To export.
       if (i_expo==1) then
          dt_ave=t_sum/(n_step+0.0_8)
          t_sum=0
          n_step=0
          write(s_n,*)n
          if (print_runtime_status==1) then
             write(*,*)trim(adjustl(s_t))//" percent finished, "//trim(adjustl(s_n))//" cells, dt=",dt_ave
          end if
          if (meas_stress==1) then
             call measure_core(t_measure)
          end if
          if (meas_force==1) then
             call get_force(t_measure)
          end if
          call export(s_label)
       endif
  
       n_temp=n
       n_divi_temp=n_divi
       call run(grow_type,n_max,propel_on,compress_on)
       if (n>n_max_clear) then
          n_divi=n_divi+(n-n_max_clear)
          n=n_max_clear
       end if
      
       if (timing_scheme==1.and.t>=t_max) then
          i_stop=1
          write(*,*)'Time limit reached.'
       endif
       if (n>=n_max.and.stop_nmax==1) then
          i_stop=1
          write(*,*)'Maximum cell number reached.'
       end if
       if (timing_scheme==3.and.n_divi>=n_max_divi) then
          i_stop=1
          write(*,*)'Maximum division reached.'
       end if
       if (i_stop==1) then
          write(s_label,*)"final"
          call export(s_label)
       end if
       if (i_error==1) then
          write(*,*)'Some ERROR encounter!'
       endif
       
       n_while=n_while+1
    enddo
  end subroutine evolve
  
  !---------------------------------------------------------------------------
  subroutine run(grow_type_t,n_max_t,propel_on_t,compress_on_t)
    use bac_rod_func
    implicit none
    integer,intent(in)::grow_type_t,n_max_t,propel_on_t,compress_on_t
    real(8),dimension(dimen)::ran
    integer::i,j,bx,by,i0
    real(8)::ma,mv,eta,rho_t
  
    f_axis(1:n)=0
    dq(1:dimen_effe,1:n)=0
    dx(1:dimen_effe,1:n)=0
    if (i_measure_stress==1) then
       stress_n_t(:,1:n)=0
    end if
    call assign_box()
    if (boundary_cond==1) then
       !Fixed boundary condition
       call boun_fixed(3)
    elseif (boundary_cond==5) then
       call boun_fixed(2)
    end if
    if (rod_or_bead==1) then
       if (substrate_on==1) then
          call force_subs_rod()
       end if
    elseif (rod_or_bead==2) then
       call assign_bead()
       if (substrate_on==1) then
          call force_subs_bead()
       end if
    elseif (rod_or_bead==3) then
       call assign_bead()
       if (substrate_on==1) then
          call force_subs_bead()
       end if
    end if
    call force_cell_cell()
    call force_random()
  
    if (compress_on_t==1) then
       call force_compress()
    endif
    if (propel_on_t==1) then
       call force_propel()
    end if
  
    if (adap_dt==1) then
       ! this step has to be after the force calculation
       call get_dt()
    endif
  
    i0=1
    ! When dumb_fixed is on, and when cell number larger than n_dumb, fix
    ! configurations of dumb cells.
    if (dumb_fixed==1.and.n>n_dumb) then
       i0=n_dumb+1
    end if
    do i=i0,n
       mv=dt/l(i)
       ma=12.0_8/l(i)**3*dt
       if (dimen_effe==2) then
          q(:,i)=q(:,i)+ma*cross_product_2d(dq(:,i),q(:,i),2)
          q(:,i)=q(:,i)/norm2(q(:,i))
       elseif (dimen_effe==3) then
          q(:,i)=q(:,i)+ma*cross_product_3d(dq(:,i),q(:,i))
          q(:,i)=q(:,i)/norm2(q(:,i))
       end if
       x(1:dimen_effe,i)=x(1:dimen_effe,i)+dx(1:dimen_effe,i)*mv
    enddo
  
    rho_t=d_cell*(sum(l(1:n))+n*d_cell)/(l_box_xt*l_box_yt)
      
    if (grow_type_t>0.and.rho_t<rho_max) then
       call grow(grow_type_t,n_max_t)
    endif
    if (boundary_cond==1) then
       call boun_fixed(3)
    elseif (boundary_cond==2) then
       call boun_periodic(3)
    elseif (boundary_cond==3) then
       if (i_measure_stress==0.and.i_measure_force==0) then  ! Only apply absorb boundary when not measuring.
          call boun_absorb(1)
       end if
       call boun_periodic(2)
    elseif (boundary_cond==4) then
       if (i_measure_stress==0.and.i_measure_force==0) then  ! Only apply absorb boundary when not measuring.
          call rescale()
       end if
       call boun_periodic(3)
    elseif (boundary_cond==5) then
       if (i_measure_stress==0.and.i_measure_force==0) then  ! Only apply absorb boundary when not measuring.
          call boun_absorb(1)
       end if
    end if
  end subroutine run
  
  !---------------------------------------------------------------------------
  