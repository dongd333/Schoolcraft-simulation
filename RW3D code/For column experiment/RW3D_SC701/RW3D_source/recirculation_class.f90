  module recirculation_class   ! class defining equilibrium reaction (it is only implemented linear sorption)
  use timefunction_class
  
  implicit none

  private
  public :: recirculation_cl                    !classes
  public ::                                 &       !methods
            read_recirculation_          ,  &
            alloc_recirculation_         ,  &
            list_recirculation_          ,  &
            read_one_connection_         ,  &
            recirculate_particle_           

  type connection_cl
      integer                        :: nExt                    !total number of extraction wells
      integer                        :: nInj                    !total number of injection wells
      integer, pointer, dimension(:) :: ExtIndex => null()      !Index of extraction wells 
      integer, pointer, dimension(:) :: InjIndex  => null()     !Index of Injection wells
      real, pointer, dimension(:)    :: Inj_Rate  => null()     !Injection rate of Injection wells      - DD
      real, pointer, dimension(:)    :: Inj_Ratio  => null()    !Cumulative ratio of of Injection rates for wells       -DD
      type(timefunction_cl)          :: timefunct
  end type

  type recirculation_cl 
      logical                                    :: action 
      integer,             pointer               :: nconnect => null()   !number of connections
      type(connection_cl), pointer, dimension(:) :: connection  => null()   !connection
  end type

  contains
  
  
    !************************************
    subroutine alloc_recirculation_ (this,nconnect)
	  use global_variables, only: nspecie
	  implicit none

	  type(recirculation_cl), intent(inout) :: this
	  integer,                intent(in)    :: nconnect	 
         if (.not. associated (this % nconnect)) then
               allocate (this % nconnect)
               this%nconnect = nconnect
         end if
         if (.not. associated (this % connection))        allocate (this % connection(nconnect))
    end subroutine


    !************************************
    function read_recirculation_ (unit,well) result(this)
	  use global_variables
      use to_solve
	  use well_vect_class
      use timefunction_class
      
	  implicit none

      integer,            intent(in)    :: unit
	  type(well_vect_cl), intent(inout) :: well
	  integer                           :: nconnect,iconnect 
	  type(recirculation_cl)            :: this
      character(len=100)                :: file
      real*8                            :: const

      read(unit,*) this%action
      recirculationACTION = this%action
	  read(unit,*) nconnect

      select case (this%action)
	  case (.TRUE.)
	     call alloc_recirculation_ (this,nconnect)
         do iconnect=1,nconnect
            call read_one_connection_ (this,iconnect,well,unit)
            read(unit,*) file,const
            call read_timefunction_ (file,const,1,this%connection(iconnect)%timefunct)
         end do
      case (.FALSE.)
         do iconnect=1,nconnect
            read(unit,*) !connections
            read(unit,*) !Injection rates
            read(unit,*) !time function
         end do
	  end select
	  
    end function
  
  
    !************************************
    subroutine read_one_connection_ (this,iconnect,well,iunit)
     use gslib, only: generate_unit,count_words,locate_blanks,upper_case
	 use well_vect_class

     implicit none

	 type(recirculation_cl),     intent(inout) :: this
	 type(well_vect_cl),         intent(inout) :: well
	 integer,                    intent(in)    :: iconnect
     integer,                    intent(in)    :: iunit
     character (len=20), allocatable           :: words(:)
     logical                                   :: exists,match
     character (len=150)                       :: string
     integer                                   :: iword,nwords,nExt,nInj,iExt,iInj
     integer                                   :: locarrow,locend,locstart,iwell
     integer, allocatable                      :: locblanks(:)
     integer                                   :: kk
     real*8                                    :: Sum_Rate          ! Sum of injection rates

     inquire (unit=iunit,exist=exists)
	 if(.not.exists) then
	     print *, 'file does not exist during read_assign_plane'
		 stop 'NOT a normal termination' 
	 end if
     read(iunit,'(a150)',err=10,end=11) string
     nwords = count_words (string)
     allocate(words(nwords))
     allocate(locblanks(1+nwords))
	 locblanks = locate_blanks(string)	  
	 do iword=1,nwords
		  words(iword)=string(locblanks(iword):locblanks(iword+1))
		  words(iword)=upper_case(words(iword))
		  if (trim(adjustl(words(iword)))== '-->') locarrow = iword
		  if (trim(adjustl(words(iword)))== '[')   locstart = iword
		  if (trim(adjustl(words(iword)))== ']')   locend   = iword
	 end do
     nExt = (locarrow - locstart)/2
     nInj  = (locend -locarrow)/2
     this%connection(iconnect)%nExt = nExt
     this%connection(iconnect)%nInj  = nInj
     if (nExt > 0) allocate (this%connection(iconnect)%ExtIndex(nExt))
     if (nInj  > 0) allocate (this%connection(iconnect)%InjIndex (nInj ))
     if (nInj  > 0) allocate (this%connection(iconnect)%Inj_Rate (nInj ))               ! Allocate Injection rate array - DD
     if (nInj  > 0) allocate (this%connection(iconnect)%Inj_Ratio (nInj ))               ! Allocate cumulative ratio array - DD
     if (nExt > 0) then
        iExt = 0
        do iword=locstart+1,locarrow,2
           iExt = iExt + 1
           match = .FALSE.
           do iwell=1,nwell_(well)
              if (trim(adjustl(words(iword)))== trim(adjustl(well%num(iwell)%name))) then
                  this%connection(iconnect)%ExtIndex(iExt) = iwell
                  well%num(iwell)%recirculation =.TRUE.
                  !well%num(iwell)%out = 1
                  match = .TRUE.
              end if
           end do
           if (.not.match) stop '>> COULD NOT MATCH RECIRCULATION WITH EXTRACTION'
        end do
     end if
     if (nInj>0) then
        iInj  = 0
        do iword=locarrow+1,locend,2
          iInj = iInj + 1
          match = .FALSE.
          do iwell=1,nwell_(well)
              if (trim(adjustl(words(iword)))== trim(adjustl(well%num(iwell)%name))) then
                 this%connection(iconnect)%InjIndex(iInj) = iwell
                 !well%num(iwell)%out = 0
                 match = .TRUE.
              end if       
          end do
          if (.not.match) write(*,*) '>> COULD NOT MATCH RECIRCULATION WITH INJECTION'
          write(*,*) this%connection(iconnect)%InjIndex(iInj), iword        ! Check which one has probem.   - Dong
        end do
        
        ! Read the injection rate for every injection well
        read(iunit,*,err=10,end=12) (this%connection(iconnect)%Inj_Rate (kk),kk=1,nInj)
        Sum_Rate = sum(this%connection(iconnect)%Inj_Rate)    ! Sum of injection rates
        this%connection(iconnect)%Inj_Ratio(1) = this%connection(iconnect)%Inj_Rate(1)/Sum_Rate     ! The ratio of injection rate for the first well
        ! Calculate cumulative ratio of injection rate, so that the injection wells are distributed from 0 to 1 based on their injection rate
        do kk=2,nInj
            this%connection(iconnect)%Inj_Ratio(kk) = this%connection(iconnect)%Inj_Ratio(kk-1) + this%connection(iconnect)%Inj_Rate(kk)/Sum_Rate
        end do
        
     end if
     return
     10   write(*,*) '>> error reading recirculation input file'
     11   write(*,*) '>> end-of-file reading recirculation input file'
     12     write(*,*) '>> Missing injection rates in the input file'           ! - DD
    stop 'NOT a normal termination'   
        
    end subroutine
  
  
    !************************************
    subroutine list_recirculation_ (this,fname)   !printer
  	  use gslib, only: open_fname
	  use global_variables, only: fdbg
	  implicit none

	  character(len=*), intent(in)            :: fname
	  type (recirculation_cl), intent (inout) :: this
	  integer                                 :: unit,iconnect,i,nconnect,nExt,nInj

      call open_fname (fdbg,unit)
      write(unit,*)
      write(unit,*) ' RECIRCULATION OF PARTICLES'
      write(unit,*)

      nconnect = this%nconnect

      do iconnect=1,nconnect
         nExt = this%connection(iconnect)%nExt
         nInj = this%connection(iconnect)%nInj
         write(unit,11) (this%connection(iconnect)%ExtIndex(i),i=1,nExt),' --> ',(this%connection(iconnect)%InjIndex(i),i=1,nInj)
   11    format( <nExt>(i2,x),a5,<nInj>(i2,x)  )   
      end do


    end subroutine
  
  
    !************************************
    subroutine recirculate_particle_ (particle,advection,geo,plume,well,recirculation,iseed)
     use particle_class
     use advection_class
     use geometry_class
     use well_vect_class
     use array_class
     use plume_class
     use gslib, only: rand2,sortem
     implicit none
     type(particle_cl),      intent(inout) :: particle
     type(advection_cl),     intent(in)    :: advection
     type(geometry_cl),      intent(in)    :: geo
     type(plume_cl),         intent(inout) :: plume
     type(well_vect_cl),     intent(in)    :: well
     type(recirculation_cl), intent(in)    :: recirculation
     integer                               :: k2,k3,k,iwell,wcol,wrow,iconnect,iExt,ninj,inj
     integer                               :: i,nd,kk,nz,zone,spe,ip
     real*8                                :: xp,yp,zp,mp,rp
     real*8                                :: sumq,qx,qy,qz,aa(1),rand,dz,cdf
     real*8, allocatable                   :: qm(:),num(:)
     integer                               :: iseed
     real*8                                :: active
     integer, dimension(1)                 :: ii
     real*8                                :: randw, check0
     integer                               :: ii0
     
         
         if (.not.recirculation%action) return
         !!if (.not.particle%control%remove) return
         !!if (particle%control%passwell == 0) return 
         !!if (.not.well%num(particle%control%passwell)%recirculation) return 
                                       
         do iconnect=1,recirculation%nconnect
         
            !check if connection is active at this time
            
            active = evaluate_timefunction_upwind_ (recirculation%connection(iconnect)%timefunct,plume%time)
            
            if (active <= 0.d0) cycle
            
         
         do iExt=1,recirculation%connection(iconnect)%nExt
             
                if (particle%control%passwell/=recirculation%connection(iconnect)%ExtIndex(iExt)) cycle
             
                ninj = recirculation%connection(iconnect)%nInj

                !estimate the mass of the particle
                mp = particle%prop%mp
                !mp = particle%prop%mp/float(ninj)
                
                ! Evoke a random number to find which well the particle will go
                randw = rand2(iseed)
                ! Find in the ratio array where the random number fall in
                ! MINLOC find the index for the minimum value
                ii= MINLOC(abs(randw-recirculation%connection(iconnect)%Inj_Ratio))
                
                ! Check if the minimum value is less than the random number, if yes, index-1
                check0 = randw-recirculation%connection(iconnect)%Inj_Ratio(ii(1))
                !write(*,*) ii, check0, randw, recirculation%connection(iconnect)%Inj_Ratio(ii(1))
                if (check0 < 0) then
                    ii0 = ii(1)
                else
                    ii0 = ii(1)+1
                end if
                
                !do inj=1,ninj
                
                   iwell = recirculation%connection(iconnect)%InjIndex(ii0)
                
                   k2 = well%num(iwell)%botlay
                   k3 = well%num(iwell)%toplay
             
                   wcol = well%num(iwell)%col
                   wrow = well%num(iwell)%row
                        
                   ! get velocities
  
                   sumq = 0.d0
       
                   allocate (qm(k2:k3),num(k2:k3))
       
                   dz = 1.d0
       
                   do k=k2,k3
                      num(k)=float(k)
                      
                      qx = value_array_ (advection%qx,wcol+1,wrow,k) - value_array_ (advection%qx,wcol,wrow,k)
                      qy = value_array_ (advection%qy,wcol,wrow+1,k) - value_array_ (advection%qy,wcol,wrow,k)
                      qz = value_array_ (advection%qz,wcol,wrow,k+1) - value_array_ (advection%qz,wcol,wrow,k)
                      qm(k) = qx+qy+qz
                      sumq = sumq + qm(k)
                  
                      !qx = 0.5d0*( value_array_ (advection%qx,wcol,wrow,k) + value_array_ (advection%qx,wcol+1,wrow,k) )
                      !qy = 0.5d0*( value_array_ (advection%qy,wcol,wrow,k) + value_array_ (advection%qy,wcol,wrow+1,k) )
                      !if (nz /= 1) dz = value_array_ (geo%dz,1,1,k)
                      !qm(k) = dsqrt(qx*qx+qy*qy)*dz
                      !sumq = sumq + qm(k)
                   end do
     
                   !construct cdf of q
             
                   rand = rand2(iseed)

                   nd = k3-k2+1

		           call sortem(1,nd,qm,1,num,aa,aa,aa,aa,aa,aa)

                   !choose layer

	               cdf   = 0.0d0
	         
	               !do i=1,nd
	               do i=k2,nd
	                   cdf = float(i)/float(nd)
	                   if (rand < cdf) then
	                        kk = dint(num(i))
	                        exit
	                   end if
	               end do

                   deallocate (qm,num)
         
                   !write(*,*) ii0, iwell, wcol, wrow, kk
                   !estimate new particle location
             
                   xp   = well%num(iwell)%xwell 
                   yp   = well%num(iwell)%ywell
                   zp   = 0.5d0*(geo%zmesh%values(1,1,kk)+geo%zmesh%values(1,1,kk+1))            
                   rp   = particle%prop%rp 
                   zone = particle%zone%num
                   spe  = particle%specie%num
                   
                   !create particle in plume
                   
                   call generate_plumeparticle_ID_ (ip)
                   call add_particle_to_plume_ (plume,ip,xp,yp,zp,mp,rp,zone,spe)
                  
                !end do

                !---DANI

                particle%control%remove = .TRUE.
                
                return               
                
         end do
         end do

    end subroutine
  
  
  !************************************************************************
  !************************************************************************

  end module recirculation_class

