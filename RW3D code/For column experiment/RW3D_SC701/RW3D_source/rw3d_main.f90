!******************************************************************************************************************
!                                        _/_/_/_/  _/      _/  _/_/_/_/  _/_/_/
!                                       _/    _/  _/      _/        _/  _/    _/
!                                      _/_/_/    _/  _/  _/    _/_/_/  _/    _/
!                                     _/    _/  _/  _/  _/        _/  _/    _/
!                                    _/    _/  _/_/_/_/_/  _/_/_/_/  _/_/_/
!
!    RW3D_Bimol_v5:
!    - 3D random walk particle tracking
!    - Transient flow (Modflow 2000)
!    - Reaction:   . Linear Sorption
!                   . Multirate mass transfer
!                   . First-order decay network
!                   . Kinetic reaction network
!
!    MAIN PROGRAM:
!******************************************************************************************************************
	 
	 program rw3d
     
     !--- Modules
     use constants
     use global_variables
	 use code_options
	 use loops_particles
     use cal_elapsed_time
	 use gslib
	 use particle_movements
	 use read_input_rw3d
	 use particle_injection
	 use checks_to_particle
	 use gslib
	 use bimolecular_reaction
	 use linear_reactions

	 !--- Classes
	 use source_vect_class
     use plane_vect_class
	 use well_vect_class
	 use plume_class
	 use geometry_class
	 use histogram_class
     use advection_class
	 use dispersion_class
	 use particle_class
	 use breakthru_class
     use list_class
     use plume_class
     use ctimes_class
	 use mass_trans_class
	 use reaction_class
	 use recirculation_class
!    _____________________________________________________________________________________________________________________

     implicit none

	 type(source_vect_cl)       :: source
	 type(plane_vect_cl)        :: plane
	 type(well_vect_cl)         :: well
	 type(plume_cl)             :: plume
	 type(geometry_cl)          :: geo, geoblock
	 type(histo_cl)             :: histo(2)
     type(advection_cl)         :: advection
	 type(dispersion_cl)        :: dispersion
	 type(mass_trans_cl)        :: mass_trans
	 type(reaction_cl)          :: reaction       
     type(particle_cl)          :: particle
	 type(ctimes_cl)            :: tc  
	 type(recirculation_cl)     :: recirculation      

     type(partID_cl),  pointer  :: plumepart

     real*8                     :: dxp_AD(3),dxp_GD(3),dxp_BR(3) !particle displacements due to ADVECTION, GRADIENT, DISPERSION
	 real*8                     :: dxp(3)                        !particle displacement in one time step
	 real*8                     :: dt                            !time step
	 integer                    :: ibdt(8)                       !clock variables to calculate elapsed run time
	 integer                    :: iseed, iseed2
	 integer                    :: ispe,izone,i,jp
	 logical                    :: remove, exist0
     
     integer                    :: seed(2), clock, seedbimo, nsim,osp,temp1,temp2, ii, ioerr
!    _____________________________________________________________________________________________________________________

     ! Read save summary of last step nitrate and CT numbers
     nsim = 1
     inquire(file="Nitrate_CT_Numbers.dat", exist=exist0) 
     if (exist0) then
         ! Read in the "Nitrate_CT_Numbers.dat" file and find the number of last file
         open(10004, file="Nitrate_CT_Numbers.dat")
         read(10004,*)
         ii = 1
         do
             read(10004,*,iostat=ioerr) osp,temp1,temp2
             IF (ioerr < 0) THEN        ! Find the last record to assign new number and exit
                 nsim = osp + 1
                 exit
             END IF
         end do
         close(10004)                   ! close the file after reading and re-open to write
         
         open(10004, file="Nitrate_CT_Numbers.dat", status="old", position="append", action="write")
     else
         open(10004, file="Nitrate_CT_Numbers.dat", status="new", action="write")
         write(10004,'(3A9)') 'Simulation', 'Nitrate', 'CT'
     end if
     
!....read parameters and initialize:

	 call print_version   ( version, date )
     call read_parameters ( geo, advection, dispersion, mass_trans, reaction, plume, well, plane, histo, source, recirculation, nsim ) 
	 call get_start_time  ( ibdt, fdbg )                                                                    !start clock and print 

!.....Start loops:
     ! Generate random seed for the random numbers, so that no need to re-compile every run
     CALL SYSTEM_CLOCK(COUNT=clock)
     seed(1) = clock/10+100000000
     seed(2) = 999999999 - clock/10
     seedbimo = clock+111000000                         ! Seed for bimolecular reaction
     iseed = mod(clock/10,38930000)
     iseed2 = 99999-mod(clock/10,38930000)             ! For re-circulation subroutine
     
          write(*,*) clock, seed,iseed

     call random_seed ( put = seed )                                                                        !set seed in random number generator

     call allocate_particle_ ( particle, mass_trans )                                                       !set default values to one dummy particle  

     moves: do nmove=1,nmaxmove     

        if ( plume%time > tsim ) exit moves
        if ( plume%time > EndTimeInjection .and. npMobile <=0 ) exit moves

        !...read flow and inject particles:

        call read_fluxes_               ( advection, geo )                                                  !read fluxes if flow changes 
        call inject_particles           ( plume, geo, advection, mass_trans, reaction, source )             !inject particles                           
        call print_number_of_particles_ ( plume )                                                           !print to screen
        call check_plume_inside_domain  ( plume , geo )                                                     !check if all particles are inside the domain (only once)       

        if (nmove   == 1) call print_number_of_particles_ ( plume, fdbg)        
        if (iwcshot == 1) call print_plume_               ( plume, files_nam(4) )                           !print snapshot if needed
        if (ixmom   == 1) call print_moments_plume_       ( plume, files_nam(6) )                           !print spatial moments if needed
 
        dt = calculate_plume_time_step ( tc, plume, geo, advection, dispersion, reaction, mass_trans )      !calculate time step

        np = get_plume_np_ (plume)

        !...loop over plume particles (organized by species and by zones): 
        
        !call reinitialize_characteristic_times ( tc )

        tc = UNEST   !initialize characteristic times

        npMobile = 0 !initialize counting mobile particles

        species: do ispe  = 1, plume%nspecie
          zones: do izone = 0, plume%nzone-1

            if (plume%np <= 0) exit species

            plumepart => plume%species(ispe)%zone(izone)%head

        !...move all particles one time step:
        
        particles: do jp=1,np(ispe,izone)

            !...exit if no more particles in plume:
            
            if (.not.associated(plumepart)) exit particles 

            remove =.FALSE.

            if (phasespecie(ispe)==0) then !skip moving particles if phase not zero

            !...create dummy particle from plume particle:
           
            call from_plumeparticle_to_particle_    ( particle, plume%time, plumepart, izone, ispe )

            !...update particle properties: 
            
            call update_cell_location_particle_     ( particle, geo )                                       !update cell location and local coordinates, save switchcell-flag
	        call update_properties_particle_        ( particle, geo, advection, dispersion, reaction )      !update properties particle poro,rpt,aL,aTH,aTV
            call update_velocity_particle_          ( particle, geo, advection, dispersion )                !update velocities qL,qT,qnode,qfaces
            call update_dispersion_nodes_particle_  ( particle, geo, dispersion )                           !update dispersion at nodes of cell where particle moves
            call update_characteristic_times_       ( particle, geo, advection, dispersion, reaction, tc )  !update characteristic times
	        call update_sorption_                   ( particle, geo, reaction )                             !update sorption parameters
            call update_mass_trans_                 ( particle, geo, mass_trans )                           !update mass transfer parameters
            call update_decay_                      ( particle, geo, reaction )                             !update decay parameters

	        call check_mobility_particle_           ( particle, advection, dispersion )                     !check if particle can move    
  	        call print_position_particle_           ( particle, files_nam(5), nmove )                       !print pathline position if needed

            call run_mass_trans_and_decay_network   ( particle, dt, iseed )                                 !run Linear Reaction: Mass transfer, sorption and first order decay network
            
            dxp = move_part_by_advection_dispersion ( advection, dispersion, particle, geo, dt )            !particle displacement in one time step
	        call add_move_to_particle_              ( particle,  dxp, dt )                                  !add displacement and incremental time to dummy particle

            call update_cell_location_particle_     ( particle, geo )                                       !update cell location and local coordinates, save switchcell-flag
            call check_boundary                     ( particle, geo )                                       !check if particle is out of system or bouncing at the boundaries
            call check_well_arrival                 ( particle, well,  reaction, recirculation )            !check if particle arrived at wells
            call check_plane_arrival                ( particle, plane, reaction )                           !check if particle arrived at planes
            call check_inactive_cell                ( particle, geo )                                       !check if inactive cell => reflect particle
            
            call from_particle_to_plumeparticle_    ( particle, plumepart )                                 !add displacement to plume particle
            call recirculate_particle_              ( particle, advection, geo, plume, well, recirculation, iseed2 )!recirculate particles when needed
            call update_particle_state_in_plume_    ( particle, plumepart, plume, remove )                  !relocate particles that have changed species or zone in the plume and delete particles when needed

            end if

            if (.not.remove) plumepart => plumepart%next                                                    !go to the next plume particle

        end do particles

        if (phasespecie(ispe)==0) npMobile = npMobile + plume%species(ispe)%zone(izone)%np

        end do zones  
        end do species

        plume%time = plume%time + dt

        !...kinetic reactions within the same zone region (mobile/immobile regions):
        
        call bimolecular_plume_reaction (plume,geo,advection,dispersion,reaction,dt,seedbimo)

        !...update velocity timeshot index:
        
        call update_advection_timeshot_ (advection,plume%time) 
                
        call print_number_of_particles_ ( plume, fdbg)
     
     end do moves
     
     ! Write the summary
     write(10004,'(3I9)') nsim, plume%species(1)%zone(0)%np, plume%species(4)%zone(0)%np

	 call analyze_BTCs_and_print_to_files   ( well, plane, histo )

!........delete memory in objects

!    call delete_plane_vect_  (plane)
!	 call delete_well_vect_   (well)
!	 call delete_plume_vect_  (plume)
!	 call delete_source_vect_ (source)
!	 call delete_geometry_    (geo)
!	 call delete_advection_   (advection)
!	 call delete_dispersion_  (dispersion)
	 call delete_mass_trans_  (mass_trans)
	 call delete_reaction_    (reaction)
	 call delete_particle_    (particle)

!........calculates elapsed run time:

	 call get_end_time (ibdt,fdbg)
     
     close(unit=10004)      ! Close the last step summary file 'Nitrate_CT_Numbers.dat'
	 
     end program rw3d





!********************************************************************************************************************************************
!
!                     PRINT RESULTS
!
!********************************************************************************************************************************************
	 subroutine analyze_BTCs_and_print_to_files (well,plane,histo)
	 
	 use gslib, only: open_fname,Delete_File
	 use global_variables
	 use code_options
	 use well_vect_class
	 use plane_vect_class
	 use histogram_class
	 use breakthru_class
	 implicit none

	 type(well_vect_cl),   intent(inout) :: well
	 type(plane_vect_cl),  intent(inout) :: plane
	 type(histo_cl),       intent(inout) :: histo(2)
	 type(breakthru_cl)                  :: btc
	 integer                             :: i,n,unit, ispe


!.....delete repeated particles in btc
!
!	      n= nwell_ (well) 
!	      do i=1,n
!	      do ispe=1,nspecie
!	      	 call delete_repeated_particles_in_btc_ (well%num(i)%btc(ispe))
!		  end do
!		  end do
!
!          n = nplane_ (plane)	  
!		  do i=1,n
!		  do ispe=1,nspecie
!             call delete_repeated_particles_in_btc_ (plane%num(i)%btc(ispe))
!		  end do
!		  end do
!
!.....analyze breakthru curves

      if (iwbtc == 1) then
	      n= nwell_ (well) 
	      do i=1,n
	      do ispe=1,nspecie
	      	 call print_well_ (well%num(i),files_nam(2),ispe)
			 !call print_well_ (well%num(i),files_nam(2))
             call print_pdf_breakthru_ (well%num(i)%btc(ispe),histo(1),0,files_nam(2))
			 call print_histo_ (histo(1),fdbg) 
		  end do
		  end do

          n = nplane_ (plane)	  
		  do i=1,n
		  do ispe=1,nspecie
		     call print_plane_ (plane%num(i),files_nam(2),ispe)
             call print_pdf_breakthru_ (plane%num(i)%btc(ispe),histo(1),0,files_nam(2))
             !call print_breakthru_ (plane%num(i)%btc(ispe),699)
			 call print_histo_ (histo(1),fdbg) 
		  end do
		  end do
	  end if

      if (iwDbtc == 1) then
	      n= nwell_ (well) 
	      do i=1,n
	      do ispe=1,nspecie
			 call print_well_ (well%num(i),files_nam(15))
			 call print_pdf_breakthru_ (well%num(i)%btc(ispe),histo(2),1,files_nam(15))
			 call print_histo_ (histo(2),fdbg) 
		  end do
		  end do

          n = nplane_ (plane)
		  do i=1,n
		  do ispe=1,nspecie
             call print_plane_ (plane%num(i),files_nam(15))
			 call print_pdf_breakthru_ (plane%num(i)%btc(ispe),histo(2),1,files_nam(15))
			 call print_histo_ (histo(2),fdbg) 
		  end do
		  end do
	  end if

      if (iwcbtc == 1) then
	      n= nwell_ (well) 
	      do i=1,n
	      do ispe=1,nspecie
	         if (well%num(i)%btc(ispe)%NP == 0) cycle 
			 call print_well_ (well%num(i),files_nam(3))
             call print_cdf_breakthru_ (well%num(i)%btc(ispe),histo(1),files_nam(3))
		  end do
		  end do

          n = nplane_ (plane)
		  do i=1,n
		  do ispe=1,nspecie
		     if (plane%num(i)%btc(ispe)%NP == 0) cycle 
			 call print_plane_ (plane%num(i),files_nam(3),ispe)
             call print_cdf_breakthru_ (plane%num(i)%btc(ispe),histo(1),files_nam(3))
		  end do
		  end do
	  end if

      if (itmom == 1) then
             call print_well_moments_vect_  (well,files_nam(11))
             call print_plane_moments_vect_ (plane,files_nam(11))
	  end if

      if (ipldisp == 1)   call print_plane_dispersivity_vect_ (plane,files_nam(12))
	  if (iwcshotpl == 1) call print_plane_breakthru_vect_    (plane,files_nam(8))

	 end subroutine


!***********************************************************************************************************
!  PRINT RESULTS ARRIVALS
!***********************************************************************************************************
	 subroutine print_summary_arrivals (plane,well,fname)
	 
	 use global_variables, only: calcul_time_method
	 use gslib, only: open_fname
	 use plane_vect_class
	 use well_vect_class
     implicit none
     type(plane_vect_cl), intent(in) :: plane
	 type(well_vect_cl),  intent(in) :: well
	 character(len=*),    intent(in) :: fname
	 integer                         :: iunit,nwell,nplane,i

      call open_fname(fname,iunit)

      write(iunit,*)
      write(iunit,'(5x,a30,a40)') 'Method Calcul Time Step.....: ',calcul_time_method
      write(iunit,*)

      nwell  = 0
	  nplane = 0

      if (associated(well%num))  nwell  = size(well%num)
      if (associated(plane%num)) nplane = size(plane%num)

	  if (nwell > 0 .or. nplane > 0 ) write(iunit,'(/,a/)') ' Particles Arriving at Control Surfaces:'

	  do i=1,nwell
	     write(iunit,'(10x,a8,i5,a12,i5)') 'well no.',i,' .........: ',well%num(i)%btc%np
	  end do

	  do i=1,nplane
 	  write(iunit,'(10x,a9,i5,a11,i5)') 'plane no.',i,' ........: ',plane%num(i)%btc%np
	  end do

	  close(iunit)

	 end subroutine


!*************************************************************************************
!   print version program
!*************************************************************************************
	 subroutine print_version(version1,date1)
	 
	 implicit none

	 character(len=*),intent(in) :: version1,date1
	 character(len=len(trim(version1))+len(trim(date1))+1) :: version

	 version = trim(adjustl(version1))//' '//trim(adjustl(date1))

   	  write(*,*)
	  write(*,*)  '                  RW3D-MRMT                    '
      write(*,*)  '      Reactive Solute Transport Program       '
	  write(*,10) '           Version ',version
	  write(*,*)  '          Daniel Fernandez-Garcia             '
	  write(*,*)

 10   format (a,a)

	 end subroutine


!****************************************************************************************************
!  PRINT PARTICLE NUMBER TO SCREEN
!****************************************************************************************************
	 subroutine print_to_screen (inj,ip)
	 
	 use global_variables, only: iwscreen
	 implicit none

	 integer, intent(in) :: ip,inj
	 integer, parameter  :: iunit = 6

	  if (iwscreen /= 0) then
	     write(iunit,*) 'injection: ',inj,' particle: ',ip
	  end if		   

	 end subroutine
!****************************************************************************************************


