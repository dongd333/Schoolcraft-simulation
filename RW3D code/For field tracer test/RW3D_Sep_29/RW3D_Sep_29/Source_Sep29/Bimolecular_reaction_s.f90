 module bimolecular_reaction

 public

 contains


!************************************************************************************************************************************************************
    subroutine bimolecular_plume_reaction (plume,geo,advection,dispersion,reaction,dt,seedbimo)
        use plume_class
        use geometry_class
        use advection_class
        use dispersion_class
        use reaction_class
        implicit none

        type (plume_cl),             intent (inout) :: plume
        type (geometry_cl),          intent (in)    :: geo
        type (advection_cl),         intent (in)    :: advection
        type (dispersion_cl),        intent (in)    :: dispersion
        type (reaction_cl),          intent (in)    :: reaction
        real*8,                      intent (in)    :: dt
        integer                                     :: izone, seedbimo

          do izone = 0, plume%nzone-1
           call bimolecular_plume_reaction_one_zone (plume,izone,geo,advection,dispersion,reaction,dt,seedbimo)
          end do
          
    end subroutine


!************************************************************************************************************************************************************
!    DAVE APPROACH
!************************************************************************************************************************************************************
    subroutine bimolecular_plume_reaction_one_zone (plume,izone,geo,advection,dispersion,reaction,dt,seedbimo)

        !USE IFPORT

        use array_class
        use plume_class
        use geometry_class
        use advection_class
        use dispersion_class
        use reaction_class
	    use particle_class
        use kdtree2_precision_module
        use kdtree2_module
        use list_class
        
		use loops_particles,    only: nmove,ip
        use constants,          only: pi        
        use global_variables,   only: UNEST,EPS,ndim,ActiveDim,phasespecie     ! Add phasespecie for reactions with immobile phase - Dong
        use gslib,              only: rand3
        use to_solve,           only: kineticACTION
        
        implicit none

        !------------------- variables-------------------------------------------------------------------------------------------------------------------------
        type (plume_cl),             intent (inout) :: plume
        type (geometry_cl),          intent (in)    :: geo
        type (advection_cl),         intent (in)    :: advection
        type (dispersion_cl),        intent (in)    :: dispersion
        type (reaction_cl),          intent (in)    :: reaction   
        type (particle_cl)                          :: particleA,particleB
        real*8,                      intent (in)    :: dt
        integer,                     intent (in)    :: izone  
        real(kdkind), allocatable                   :: Bloc(:,:)
        real(kdkind), allocatable                   :: Alocvec(:),Brad2(:)      ! These are types for the kdtree searches
        type(kdtree2),pointer                       :: Btree                    ! These are types for the kdtree searches       
        type(kdtree2_result), allocatable           :: nearest(:)               ! This (nearest) will have closest points by kdtree. %idx is indices and %dis is distances^2 
        integer                                     :: nf,idim
        real*8                                      :: maxsearch,gaussfact,vofs,sep(3,1),prob,R
        real*8                                      :: r3, yield, gy            ! Add for growth reaction - Dong
        integer                                     :: npA,npB,npAstart,npAend,npBstart,npBend
        integer                                     :: i,ii,numAnow,nreact,ncount,ipA,ipB,igr0
        real*8                                      :: h2A(3,3),h2B(3,3),h2Binv(3,3),det,kf,s2,s22, temp(3,1), temp2(3,1), temp3(1,1)   !particle support
        real*8                                      :: w(3,3),v(3)
        real*8, pointer                             :: xp(:,:)
        logical                                     :: testDIM(3) = .TRUE.
        integer                                     :: irxn,sp1,sp2,spA,spB,spC,spD, spE
        character(len=1)                            :: spcat
        logical                                     :: catalytic =.FALSE.
        real*8, pointer, dimension(:)               :: xpA
        type(partID_cl),  pointer                   :: partA,partB
        type(list_cl), pointer                      :: list
        logical                                     :: react
        logical, allocatable                        :: Breact(:)
        real*8                                      :: xpC,ypC,zpC,mpC
        real*8                                      :: xpD,ypD,zpD,mpD
        integer                                     :: ipC,ipD
        real*8                                      :: xpE,ypE,zpE,mpE          ! Add for growth reaction - Dong
        real*8                                      :: z1,z2,z3,zz1,zz2,zz3,dif_x,dif_y,dif_z,dif2_x,dif2_y,dif2_z               ! Add to place product particle with diffusion     - Dong
        real*8                                      :: dm,dmTH,dmTV
        
        integer                                     :: seedbimo
        !--------------------------------------------------------------------------------------------------------------------------------------------------------
        ! DAVE: gaussfact=(8.*pi*D*dt)**(real(ndim)/2.)
        ! DAVE: maxsearch=sqrt(24.0*D*dt)
        ! DAVE: vofs=exp(sep2/(-8.*D*dt))/gaussfact
        ! DANI: h2 = 2*pi*dt
        !--------------------------------------------------------------------------------------------------------------------------------------------------------
        !
        !NOTE: we constuct the kdtree of B particles once. This means that particle A can be deleted on the fly but particle B not. 
        !We keep track of how many B particles react and delete them in the end.

        if (.not.kineticACTION) return

        if (plume%np <= 0) return

        ! for each chemical reaction:
        igr0 = 0            ! Add for growth reactions      - Dong

        do irxn = 1,reaction%kinetic%nreact

          !-----------------------------------------------------------------------------------------------------------------------------------------------------
          ! ONE REACTANT AND ZERO PRODUCTS: A --> 0
          !-----------------------------------------------------------------------------------------------------------------------------------------------------

          if (reaction%kinetic%rxn(irxn)%nreactants == 1 .and. reaction%kinetic%rxn(irxn)%nproducts == 0) then

              !initialize
              
              spA = reaction%kinetic%rxn(irxn)%reactants(1)                                             !species number 
              npA = plume%species(spA)%zone(izone)%np                                                   !number of particles A in list
              
			  if( npA <=0 ) return
              
              !loop over particles in plume list 
              
              partA => plume%species(spA)%zone(izone)%head 
              
              do i=1,npA           
                 
                 call from_plumeparticle_to_particle_      ( particleA, plume%time, partA, izone, spA)  !create particle
                 call update_cell_location_particle_       ( particleA, geo )                           !update cell location and local coordinates, save switchcell-flag	                            
                 kf   = value_array_ (reaction%kinetic%kf(irxn), particleA%cell%num)                    !get reaction coefficient                 
				 prob = 1.d0-dexp(- dt*kf)                                                              !probability of reaction to occur          
                 R = rand3(seedbimo)                                                                    !random number between 1 and 0
                 select case (Prob > R) 
                    case (.TRUE.);  call delete_plumeparticle_  (plume,izone,spA,partA)                 !reaction occurrs, so kill particle and go to next one                                                                             
                    case (.FALSE.); partA => partA%next                                                 !reaction does not occur, go to next particle in plume list
                 end select                           
              
              end do
          !-----------------------------------------------------------------------------------------------------------------------------------------------------
          ! ONE REACTANT AND ONE PRODUCT: A --> C
          !-----------------------------------------------------------------------------------------------------------------------------------------------------

          elseif (reaction%kinetic%rxn(irxn)%nreactants == 1 .and. reaction%kinetic%rxn(irxn)%nproducts == 1) then
          
              !initialize
              
              spA = reaction%kinetic%rxn(irxn)%reactants(1)                                             !species number 
              npA = plume%species(spA)%zone(izone)%np                                                   !number of particles A in list
              
			  if( npA <=0 ) return
              
              !loop over particles in plume list 
              
              partA => plume%species(spA)%zone(izone)%head 
              
              do i=1,npA           
                 
                 call from_plumeparticle_to_particle_      ( particleA, plume%time, partA, izone, spA)  !create particle
                 call update_cell_location_particle_       ( particleA, geo )                           !update cell location and local coordinates, save switchcell-flag	                            
                 kf   = value_array_ (reaction%kinetic%kf(irxn), particleA%cell%num)                    !get reaction coefficient                 
				 prob = 1.d0-dexp(- dt*kf)                                                              !probability of reaction to occur          
                 R = rand3(seedbimo)                                                                    !random number between 1 and 0
                 select case (Prob > R) 
                 case (.TRUE.)    !reaction occurrs, so kill particle and go to next one
                           spC = reaction%kinetic%rxn(irxn)%products(1)
                           
                           ! Generate a random number at each dimension     - Dong
                           z1 = 2*rand3(seedbimo)-1
                           z2 = 2*rand3(seedbimo)-1
                           z3 = 2*rand3(seedbimo)-1
                           ! Calculate the one diffusion distance at each dimension     - Dong
                           dif_x = dsqrt( 2.d0 * dt * dispersion%dm%values(1,1,1) * dispersion%MultD(spC) ) 
                           dif_y = dsqrt( 2.d0 * dt * dispersion%dmTH%values(1,1,1) * dispersion%MultD(spC) ) 
                           dif_z = dsqrt( 2.d0 * dt * dispersion%dmTV%values(1,1,1) * dispersion%MultD(spC) ) 
                           
                           xpC = particleA%position%xp(1)+z1*dif_x
                           ypC = particleA%position%xp(2)+z2*dif_y
                           zpC = particleA%position%xp(3)+z3*dif_z
                           mpC = partA%mp                                                                       
                           call generate_plumeparticle_ID_ (ipC)
                           call add_particle_to_plume_ (plume,ipC,xpC,ypC,zpC,mpC,1.d0,izone,spC)
                           call delete_plumeparticle_  (plume,izone,spA,partA)  !kill particle A in list and go to next A particle                                                                                               
                    case (.FALSE.)
                        partA => partA%next                                     !reaction does not occur, go to next particle in plume list
                 end select                           
              
              end do

          !-----------------------------------------------------------------------------------------------------------------------------------------------------
          ! ONE REACTANT AND TWO PRODUCTs: A --> C + D
          !-----------------------------------------------------------------------------------------------------------------------------------------------------

          elseif (reaction%kinetic%rxn(irxn)%nreactants == 1 .and. reaction%kinetic%rxn(irxn)%nproducts == 2) then

              !initialize

              spA = reaction%kinetic%rxn(irxn)%reactants(1)                                             !species number 
              npA = plume%species(spA)%zone(izone)%np                                                   !number of particles A in list

			  if( npA <=0 ) return

              !loop over particles in plume list 

              partA => plume%species(spA)%zone(izone)%head 

              do i=1,npA           

                 call from_plumeparticle_to_particle_      ( particleA, plume%time, partA, izone, spA)  !create particle
                 call update_cell_location_particle_       ( particleA, geo )                           !update cell location and local coordinates, save switchcell-flag	                            
                 kf   = value_array_ (reaction%kinetic%kf(irxn), particleA%cell%num)                    !get reaction coefficient                 
				 prob = 1.d0-dexp(- dt*kf)                                                              !probability of reaction to occur          
                 R = rand3(seedbimo)                                                                    !random number between 1 and 0
                 select case (Prob > R) 
                    case (.TRUE.)    !reaction occurrs, so kill particle and go to next one                                                                            
                           spC = reaction%kinetic%rxn(irxn)%products(1)
                           ! Generate a random number at each dimension
                           z1 = 2*rand3(seedbimo)-1
                           z2 = 2*rand3(seedbimo)-1
                           z3 = 2*rand3(seedbimo)-1
                           ! Calculate the one diffusion distance at each dimension
                           dif_x = dsqrt( 2.d0 * dt * dispersion%dm%values(1,1,1) * dispersion%MultD(spC) ) 
                           dif_y = dsqrt( 2.d0 * dt * dispersion%dmTH%values(1,1,1) * dispersion%MultD(spC) ) 
                           dif_z = dsqrt( 2.d0 * dt * dispersion%dmTV%values(1,1,1) * dispersion%MultD(spC) ) 
                           xpC = particleA%position%xp(1)+z1*dif_x
                           ypC = particleA%position%xp(2)+z2*dif_y
                           zpC = particleA%position%xp(3)+z3*dif_z
                           mpC = partA%mp                                                                       
                           
                           spD = reaction%kinetic%rxn(irxn)%products(2)
                           ! Generate a random number at each dimension
                           zz1 = 2*rand3(seedbimo)-1
                           zz2 = 2*rand3(seedbimo)-1
                           zz3 = 2*rand3(seedbimo)-1
                           ! Calculate the one diffusion distance at each dimension
                           dif2_x = dsqrt( 2.d0 * dt * dispersion%dm%values(1,1,1) * dispersion%MultD(spD) ) 
                           dif2_y = dsqrt( 2.d0 * dt * dispersion%dmTH%values(1,1,1) * dispersion%MultD(spD) ) 
                           dif2_z = dsqrt( 2.d0 * dt * dispersion%dmTV%values(1,1,1) * dispersion%MultD(spD) )
                           xpD = particleA%position%xp(1)+zz1*dif2_x  
                           ypD = particleA%position%xp(2)+zz2*dif2_y
                           zpD = particleA%position%xp(3)+zz3*dif2_z
                           mpD = mpC                                                                       
                           
                           call generate_plumeparticle_ID_ (ipC)                         
                           call add_particle_to_plume_ (plume,ipC,xpC,ypC,zpC,mpC,1.d0,izone,spC)
                           call generate_plumeparticle_ID_ (ipD)                         
                           call add_particle_to_plume_ (plume,ipD,xpD,ypD,zpD,mpD,1.d0,izone,spD)
                           call delete_plumeparticle_  (plume,izone,spA,partA)  !kill particle A in list and go to next A particle                                                                                               
                    case (.FALSE.)
                        partA => partA%next                                     !reaction does not occur, go to next particle in plume list
                 end select                           
              
              end do
        
              
          !-----------------------------------------------------------------------------------------------------------------------------------------------------
              !! Test adding growth _DD
          ! ONE REACTANT AND Three PRODUCTs: A --> C + D + Z
              ! Z is a hypothetical species, not used. It is for simulating the growth of C. The reaction resembles A -> (1+yield)*C + D
          !-----------------------------------------------------------------------------------------------------------------------------------------------------
              
          elseif (reaction%kinetic%rxn(irxn)%nreactants == 1 .and. reaction%kinetic%rxn(irxn)%nproducts == 3) then
              !write(*,*) 'Three product!'
              igr0 = igr0 + 1
              !initialize

              spA = reaction%kinetic%rxn(irxn)%reactants(1)                                             !species number 
              npA = plume%species(spA)%zone(izone)%np                                                   !number of particles A in list

			  if( npA <=0 ) return

              !loop over particles in plume list 

              partA => plume%species(spA)%zone(izone)%head 

              do i=1,npA           

                 call from_plumeparticle_to_particle_      ( particleA, plume%time, partA, izone, spA)  !create particle
                 call update_cell_location_particle_       ( particleA, geo )                           !update cell location and local coordinates, save switchcell-flag	                            
                 kf   = value_array_ (reaction%kinetic%kf(irxn), particleA%cell%num)                    !get reaction coefficient                 
				 prob = 1.d0-dexp(- dt*kf)                                                              !probability of reaction to occur          
                 R = rand3(seedbimo)                                                                    !random number between 1 and 0
                 
                 gy   = value_array_ (reaction%kinetic%gy(igr0), particleA%cell%num)                    !get growth yield coefficient
                 !yield = 1.d0-dexp(- dt*gy)                                           ! To simulate growth yield
                 r3 = rand3(seedbimo)
                 select case (Prob > R) 
                    case (.TRUE.)    !reaction occurrs, so kill particle and go to next one                                                                            
                           spC = reaction%kinetic%rxn(irxn)%products(1)
                           ! Generate a random number at each dimension
                           z1 = 2*rand3(seedbimo)-1
                           z2 = 2*rand3(seedbimo)-1
                           z3 = 2*rand3(seedbimo)-1
                           ! Calculate the one diffusion distance at each dimension
                           dif_x = dsqrt( 2.d0 * dt * dispersion%dm%values(1,1,1) * dispersion%MultD(spC) ) 
                           dif_y = dsqrt( 2.d0 * dt * dispersion%dmTH%values(1,1,1) * dispersion%MultD(spC) ) 
                           dif_z = dsqrt( 2.d0 * dt * dispersion%dmTV%values(1,1,1) * dispersion%MultD(spC) ) 
                           xpC = particleA%position%xp(1)+z1*dif_x
                           ypC = particleA%position%xp(2)+z2*dif_y
                           zpC = particleA%position%xp(3)+z3*dif_z
                           mpC = partA%mp                                                                       
                           
                           spD = reaction%kinetic%rxn(irxn)%products(2)
                           ! Generate a random number at each dimension
                           zz1 = 2*rand3(seedbimo)-1
                           zz2 = 2*rand3(seedbimo)-1
                           zz3 = 2*rand3(seedbimo)-1
                           ! Calculate the one diffusion distance at each dimension
                           dif2_x = dsqrt( 2.d0 * dt * dispersion%dm%values(1,1,1) * dispersion%MultD(spD) ) 
                           dif2_y = dsqrt( 2.d0 * dt * dispersion%dmTH%values(1,1,1) * dispersion%MultD(spD) ) 
                           dif2_z = dsqrt( 2.d0 * dt * dispersion%dmTV%values(1,1,1) * dispersion%MultD(spD) )
                           xpD = particleA%position%xp(1)+zz1*dif2_x  
                           ypD = particleA%position%xp(2)+zz2*dif2_y
                           zpD = particleA%position%xp(3)+zz3*dif2_z
                           mpD = mpC
                           
                            ! Test calling diffusion coefficient
                           !write(*,*) dispersion%dm%values(1,1,1), dispersion%dm%values(1,1,1) * dispersion%MultD(spC), dispersion%dm%values(1,1,1) * dispersion%MultD(spD)
                           
                           call generate_plumeparticle_ID_ (ipC)                         
                           call add_particle_to_plume_ (plume,ipC,xpC,ypC,zpC,mpC,1.d0,izone,spC)
                           call generate_plumeparticle_ID_ (ipD)                         
                           call add_particle_to_plume_ (plume,ipD,xpD,ypD,zpD,mpD,1.d0,izone,spD)
                           If (gy > r3) then
                               xpE = particleA%position%xp(1)
                               ypE = particleA%position%xp(2)
                               zpE = particleA%position%xp(3)
                               mpE = mpC
                               spE = reaction%kinetic%rxn(irxn)%products(1)
                               call generate_plumeparticle_ID_ (ipC)
                               call add_particle_to_plume_ (plume,ipC,xpE,ypE,zpE,mpE,1.d0,izone,spE)
                           End IF
                           call delete_plumeparticle_  (plume,izone,spA,partA)  !kill particle A in list and go to next A particle
                           
                    case (.FALSE.)
                        partA => partA%next                                     !reaction does not occur, go to next particle in plume list
                 end select                           
                            
              end do
              !write(*,*) gy, igr0            

          !-----------------------------------------------------------------------------------------------------------------------------------------------------
          ! TWO REACTANTS AND SEVERAL PRODUCTS: A + B --> C + D
          !-----------------------------------------------------------------------------------------------------------------------------------------------------
          
          elseif (reaction%kinetic%rxn(irxn)%nreactants == 2) then

             spA = reaction%kinetic%rxn(irxn)%reactants(1)
             spB = reaction%kinetic%rxn(irxn)%reactants(2)
             
             If (phasespecie(spA) >0) stop '>> The first species of the bimolecular reaction is immobile'       ! Add a simple check    - Dong
             
             npA = plume%species(spA)%zone(izone)%np 
             npB = plume%species(spB)%zone(izone)%np 
             
             if( npA <=0 .or. npB<=0 ) return                   

             allocate(Alocvec(ndim))
             allocate(nearest(npB))
             
             allocate(Breact(npB))
             
             Breact = .FALSE.
             
             nearest%idx = int(UNEST)
             nearest%dis = UNEST      

             ! Assign position of B particles:
             
             allocate(Bloc(ndim,npB))

             partB => plume%species(spB)%zone(izone)%head
             			                 
             do i=1,npB
                !if (i==2639) then
                !  continue
                !end if
             
                    testDIM = .TRUE.
                    do idim=1,ndim
                           if (ActiveDim(1).and.testDIM(1)) then
                                    Bloc(idim,i) = partB%xp
                                    testDIM(1) = .FALSE.
                           else if (ActiveDim(2).and.testDIM(2)) then
                                    Bloc(idim,i) = partB%yp
                                    testDIM(2) = .FALSE.
                           else if (ActiveDim(3).and.testDIM(3)) then
                                    Bloc(idim,i) = partB%zp
                                    testDIM(3) = .FALSE.
                           end if
                     end do
                     partB => partB%next                   
             end do

             nullify(partB)

             ! Create Kdtree associated to B particles:
        
             Btree => kdtree2_create(Bloc,sort=.false.,rearrange=.false.)  

             ! For each A particle see whether particle A reacts with B:
             
             partA => plume%species(spA)%zone(izone)%head 

      loopA: do i=1,npA !loop particles A

                 if( plume%species(spA)%zone(izone)%np <=0 .or. plume%species(spB)%zone(izone)%np <=0 ) exit 
                 
                 call from_plumeparticle_to_particle_ ( particleA, plume%time, partA, izone, spA) 

                 call update_cell_location_particle_  ( particleA, geo )                                    !update cell location and local coordinates, save switchcell-flag	        
	             call update_properties_particle_     ( particleA, geo, advection, dispersion, reaction   ) !update properties particle poro,rpt,aL,aTH,aTV
                 call update_velocity_particle_       ( particleA, geo, advection, dispersion )             !update velocities qL,qT,qnode,qfaces            

                 testDIM = .TRUE.
                 do idim=1,ndim
                        if (ActiveDim(1).and.testDIM(1)) then
                                    alocvec(idim) = partA%xp 
                                    testDIM(1) = .FALSE.
                        else if (ActiveDim(2).and.testDIM(2)) then
                                    alocvec(idim) = partA%yp 
                                    testDIM(2) = .FALSE.
                        else if (ActiveDim(3).and.testDIM(3)) then
                                    alocvec(idim) = partA%zp
                                    testDIM(3) = .FALSE.
                        end if
                 end do
                 
                 h2A = calculate_particle_support_square ( particleA, dispersion, dt )
                 !write(*,*) spA, h2A(1,1), dispersion%dm%values(1,1,1) * dispersion%MultD(spA)

                 maxsearch = dsqrt(12.0d0 * maxval(h2A))
                 
                 ! search nearest particles
                 
                 call kdtree2_r_nearest(tp=Btree,qv=Alocvec,r2=maxsearch*maxsearch,nfound=nf,nalloc=npB,results=nearest)
                
                 ! check if any B particle reacts       
          
                 react = .FALSE.
          
          loopB: do ii=1,nf 
                    
                    ipB = nearest(ii)%idx
                    
                    call allocate_particle_ (particleB)
                    
 	                particleB % position % xp(1) = partA%xp    
		            particleB % position % xp(2) = partA%yp     
		            particleB % position % xp(3) = partA%zp 
                    
                    particleB % prop % aTH  = dispersion%aTH%values(1,1,1) * dispersion%MultA(spB)
                    particleB % prop % aTV  = dispersion%aTV%values(1,1,1) * dispersion%MultA(spB)
                    particleB % prop % dm   = dispersion%dm%values(1,1,1) * dispersion%MultD(spB)
                    particleB % specie % num = spB
		                                
                    testDIM = .TRUE.
                    do idim=1,ndim
                        if (ActiveDim(1).and.testDIM(1)) then
                                    particleB % position % xp(1) = Btree%the_data(idim,ipB) 
                                    testDIM(1) = .FALSE.
                        else if (ActiveDim(2).and.testDIM(2)) then
                                    particleB % position % xp(2) = Btree%the_data(idim,ipB)
                                    testDIM(2) = .FALSE.
                        else if (ActiveDim(3).and.testDIM(3)) then
                                    particleB % position % xp(3) = Btree%the_data(idim,ipB)
                                    testDIM(3) = .FALSE.
                        end if
                    end do
                    
                    if (Breact(ipB)) cycle
            
                    call update_cell_location_particle_         ( particleB, geo )                                    !update cell location and local coordinates, save switchcell-flag	        
	                call update_properties_particle_            ( particleB, geo, advection, dispersion, reaction   ) !update properties particle poro,rpt,aL,aTH,aTV
                    call update_velocity_particle_              ( particleB, geo, advection, dispersion )             !update velocities qL,qT,qnode,qfaces            

                    sep(:,1) = particleA%position%xp - particleB%position%xp  !distance vector

                    ! If the reactant is immoible, the array of h2B is all 0        - Dong
                    ! This loop is not needed because the species-dependent dispersion is updated.
                    !If (phasespecie(spB) >0) then
                    !    h2B      = 0
                    !Else
                        h2B      = calculate_particle_support_square (particleB,dispersion,dt)
                    !End If
                    !write(*,*) spB, h2B(1,1), dispersion%dm%values(1,1,1) * dispersion%MultD(spB)

                    h2B      = h2A + h2B
                  
                    if (ndim<3) then				  
				         if (h2B(1,1)<EPS) h2B(1,1) = 1.d0 
                         if (h2B(2,2)<EPS) h2B(2,2) = 1.d0                       
                         if (h2B(3,3)<EPS) h2B(3,3) = 1.d0				                       
				    end if
				                   
				    det      = dabs(M33DET(h2B))
                  
				    h2Binv   = M33INV(h2B)

				    gaussfact = dsqrt( det * ( (2.d0 * pi )**(dble(ndim)) ) ) 

				    s2 = maxval(matmul(transpose(sep),matmul(h2Binv,sep))) 
				    !temp2 = matmul(h2Binv,sep)
				    !temp3 = matmul(transpose(sep),temp2)
				    !call DGEMM('N','N',3,1,3,1.d0,h2Binv,3,sep,3,0.d0,temp,3)
				    !call DGEMM('T','N',1,1,3,1.d0,sep,3,temp,3,0.d0,temp2,1)
				    !s22 = maxval(temp2)
                  
				    vofs = dexp( - 0.5d0 * s2 ) / gaussfact

                    kf   = value_array_ (reaction%kinetic%kf(irxn), particleB%cell%num)
                  
				    prob = partA%mp * vofs * dt * kf

                    R = rand3(seedbimo)

                               if( prob > R )then
                               
                                  react = .TRUE.  !reaction occurrs, so transform A into products
                                                                                                    
                                    select case (reaction%kinetic%rxn(irxn)%nproducts)

                                        !--------------------------------
                                        !  ZERO PRODUCTS: A + B --> 0
                                        !--------------------------------
                                    
                                        case (0) !no products
                                             
                                              call delete_plumeparticle_  (plume,izone,spA,partA) !kill particle A                                                                            
                                              Breact(ipB) = .TRUE.                                !mark particle B to be killed after loop                                                                                                                           
                                              exit loopB                                                                              

                                        !--------------------------------
                                        !  ONE PRODUCT: A + B --> C
                                        !--------------------------------
                                        
                                        case (1) !we have only one product
                                              ! Generate a random number at each dimension
                                              z1 = 2*rand3(seedbimo)-1
                                              z2 = 2*rand3(seedbimo)-1
                                              z3 = 2*rand3(seedbimo)-1
                                              ! Randomly distributed in between
                                              xpC = abs(partA%xp+(particleB%position%xp(1)-partA%xp)*z1)
                                              ypC = abs(partA%yp+(particleB%position%xp(2)-partA%yp)*z2)
                                              zpC = abs(partA%zp+(particleB%position%xp(3)-partA%zp)*z3)
                                              mpC = partA%mp
                                              !zone = kinetic%zone(kinetic%rxn(irxn)%products(1))                                              
                                              spC = reaction%kinetic%rxn(irxn)%products(1)
                                              call generate_plumeparticle_ID_ (ipC)
                                              call add_particle_to_plume_ (plume,ipC,xpC,ypC,zpC,mpC,1.d0,izone,spC)
                                              call delete_plumeparticle_  (plume,izone,spA,partA)  !kill particle A in list and go to next A particle                                                                            
                                              Breact(ipB)  = .TRUE.                                   !mark particle B to be killed after loop
                                              exit loopB
                                    
                                        !--------------------------------
                                        !  TWO PRODUCTS: A + B --> C + D
                                        !--------------------------------
                                        
                                        case (2) !we have two products

                                              ! Generate a random number at each dimension
                                              z1 = 2*rand3(seedbimo)-1
                                              z2 = 2*rand3(seedbimo)-1
                                              z3 = 2*rand3(seedbimo)-1
                                              ! Randomly distributed in between
                                              xpC = abs(partA%xp+(particleB%position%xp(1)-partA%xp)*z1)
                                              ypC = abs(partA%yp+(particleB%position%xp(2)-partA%yp)*z2)
                                              zpC = abs(partA%zp+(particleB%position%xp(3)-partA%zp)*z3)
                                              mpC = partA%mp
                                              !zone = kinetic%zone(kinetic%rxn(irxn)%products(2))
                                              spC = reaction%kinetic%rxn(irxn)%products(1)
                                              call generate_plumeparticle_ID_ (ipC)                  
                                              call add_particle_to_plume_ (plume,ipC,xpC,ypC,zpC,mpC,1.d0,izone,spC)
                                              
                                              ! Generate a random number at each dimension
                                              zz1 = 2*rand3(seedbimo)-1
                                              zz2 = 2*rand3(seedbimo)-1
                                              zz3 = 2*rand3(seedbimo)-1
                                              ! Randomly distributed in between
                                              xpD = abs(partA%xp+(particleB%position%xp(1)-partA%xp)*zz1)
                                              ypD = abs(partA%yp+(particleB%position%xp(2)-partA%yp)*zz2)
                                              zpD = abs(partA%zp+(particleB%position%xp(3)-partA%zp)*zz3)
                                              mpD = mpC
                                              !zone = kinetic%zone(kinetic%rxn(irxn)%products(2))                                              
                                              spD = reaction%kinetic%rxn(irxn)%products(2)                                                
                                              call generate_plumeparticle_ID_ (ipD)
                                              call add_particle_to_plume_ (plume,ipD,xpD,ypD,zpD,mpD,1.d0,izone,spD)
                                              
                                              call delete_plumeparticle_  (plume,izone,spA,partA)     !kill particle A in list and go to next A particle                                                                            
                                              Breact(ipB)  = .TRUE.                                   !mark particle B to be killed after loop
                                             
                                              exit loopB                                  
                                    
                                    end select                                                                                                                   
                               
                               end if
                
                end do loopB
                
                if (.not.react) partA => partA%next ! go to next particle

             end do loopA

             call delete_plumeparticle_  (plume,Breact,spB,izone)  !kill all particles B that were removed in loopA
             
             deallocate (Breact)

             call kdtree2_destroy(Btree)
             deallocate(Alocvec,Bloc,nearest)
        
          end if
        
        end do !loop reactions
 
    end subroutine


!*************************************************************************************
!      Estimate h with dispersive motion
!*************************************************************************************
!    A = P*D*P^t
!-------------------------------------------------------------------------------------
   function calculate_particle_support_square (particle,dispersion,dt) result (A)
       use dispersion_class
       use particle_class
	   implicit none
       type(dispersion_cl), intent(in) :: dispersion 
       type(particle_cl),   intent(in) :: particle
	   real*8,              intent(in) :: dt
       real*8 :: poro,qx,qy,qz,al,ath,atv,dm,rpt,dmTH,dmTV
	   real*8 :: vx,vy,vz,v,vh,z1x,z1y,z1z,z2x,z2y,z3x,z3y,z3z,z1,z2,z3
	   real*8 :: v2,vh2
	   real*8 :: P(3,3),D(3,3),A(3,3), temp(3,3)

       D = 0.d0
       P = 0.d0
       A = 0.d0

       if ( .not.dispersion%action ) return

       if ( dt <= 0.d0 ) return 

       poro  = particle % prop % poro 
	   qx   = particle % vel  % qpT(1)
	   qy   = particle % vel  % qpT(2)
	   qz   = particle % vel  % qpT(3)

	   dm     = particle % prop % dm
	   dmTH   = particle % prop % dmTH	   
	   dmTV   = particle % prop % dmTV	   
	   
	   rpt  = particle % prop % rp

       vx = qx / poro
	   vy = qy / poro
	   vz = qz / poro

       v2 = vx * vx + vy * vy + vz * vz

       if (v2 == 0.d0) then; v = 0.d0
	                   else; v = dsqrt(vx * vx + vy * vy + vz * vz); endif

       vh2 = vx * vx + vy * vy

       if (vh2 == 0.d0) then; vh = 1.d-21
	                    else; vh = dsqrt( vx * vx + vy * vy); endif 


!------- special cases:
       
	   if ( v == 0.d0 .and. dm /=0.d0 ) then   !if no velocity allow molecular diffusion

             A(1,1) = (2.d0*(dm)/rpt*dt)     
             A(2,2) = (2.d0*(dmTH)/rpt*dt)
             A(3,3) = (2.d0*(dmTV)/rpt*dt)       
	         return
	   
	   else if ( particle%zone%num > 0 .and. dm /=0.d0 ) then !if immobile region allow molecular diffusion

             A(1,1) = (2.d0*(dm)/rpt*dt)     
             A(2,2) = (2.d0*(dmTH)/rpt*dt)
             A(3,3) = (2.d0*(dmTV)/rpt*dt)       
	         return

	   else if ( v  == 0.d0 .and. dm ==0.d0 ) then
	        
			 return
	   
	   else if ( particle%zone%num > 0 .and. dm ==0.d0 ) then
	   
	         return

	   end if

!--------- General Case:
!
	   al   = particle % prop % aL
	   ath  = particle % prop % aTH
	   atv  = particle % prop % aTV

!--------- eigen values:

       D(1,1) = (2.d0*(al *v+dm)/rpt*dt)
       D(2,2) = (2.d0*(ath*v+dmTH)/rpt*dt)
       D(3,3) = (2.d0*(atv*v+dmTV)/rpt*dt)        

!------- calculates transformation matrix:

        P(1,1) = vx/v
        P(2,1) = vy/v
        P(3,1) = vz/v
        
        P(1,2) =  vy*vh/v/v +  vz*vz*vy/v/v/vh
        P(2,2) = -vx*vh/v/v - vz*vz*vx/v/v/vh
!       P(3,2) = 0.

        P(1,3) = -vz*vx/v/vh
        P(2,3) = -vz*vy/v/vh
        P(3,3) =  vh/v         

!------- calculates Dispersion matrix in fixed (Eulerian) coordinates:
        
        
        !A = matmul(matmul(P,D),transpose(P))
        call DGEMM('N','N',3,3,3,1.d0,P,3,D,3,0.d0,temp,3)
        call DGEMM('N','t',3,3,3,1.d0,temp,3,P,3,0.d0,A,3)


   end function

!***********************************************************************************************************************************
!
!                                                       M 3 3 D E T _ M A I N
!
!  Program:      M33DET_MAIN
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         July 22, 2005
!
!  Language:     Fortran-90
!
!  Version:      1.00a
!
!  Description:  This program is a short "driver" to call function M33DET, which computes the determinant a 3x3 matrix.
!
!  Files:        Source files:
!
!                   m33det.f90                   Main program
!
!********************************************************************************************************************************
!***********************************************************************************************************************************
!  M33DET  -  Compute the determinant of a 3x3 matrix.
!***********************************************************************************************************************************

      FUNCTION M33DET (A) RESULT (DET)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A

      DOUBLE PRECISION :: DET


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      RETURN

      END FUNCTION M33DET
      
!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

      FUNCTION M33INV (A) RESULT (AINV)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3) :: AINV
      LOGICAL :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         STOP 'COULD NOT MAKE THE INVERSION WITH M33INV'
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END FUNCTION M33INV
      
!-------------------------------------------------------------------------------------------------------
!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!-------------------------------------------------------------------------------------------------------
REAL FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
	REAL, DIMENSION(n,n) :: matrix
    REAL :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
    
END FUNCTION FindDet


 end module