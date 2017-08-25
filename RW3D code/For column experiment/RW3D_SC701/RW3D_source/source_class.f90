

!******************************************************************************
!  SOURCE CLASS 
!******************************************************************************
 module source_class
   use timefunction_class
   implicit none

   private
   public :: source_cl             !class
   public ::                        & ! methods
                print_source_     , &
				delete_source_    , &
				read_source_      
				  
   type parameters_cl
	   character(len=200), pointer :: file => null() 
	   real*8,  pointer :: xinj   => null(), yinj   => null(), zinj   => null()
	   real*8,  pointer :: zbot   => null(), ztop   => null(), rcyr   => null()
	   real*8,  pointer :: rcp    => null(), xdist  => null(), width  => null(), height => null()
       real*8,  pointer :: xinj_1 => null(), yinj_1 => null(), zinj_1 => null(), xinj_2 => null()
	   real*8,  pointer :: yinj_2 => null(), zinj_2 => null()
	   integer, pointer :: idwn   => null(), jdwn   => null(), kdwn   => null()
	   integer, pointer :: iup    => null(), jup    => null(), kup    => null()
	   real*8,  pointer :: np11x  => null(), np11y  => null(), np11z  => null()
	   real*8,  pointer :: const
	   integer, pointer :: ix(:)  => null(), iy(:)  => null(), iz(:)  => null()
   end type 
   
   
   type source_cl
	   character(len=200)     :: name 
       character(len=10)      :: TypeInj
       real*8                 :: TimeStartInj
       real*8                 :: TimeStopInj
       integer                :: np 
	   real*8                 :: pmass 
	   integer                :: zone
	   integer                :: specie
       type(parameters_cl)    :: par
       type(timefunction_cl)  :: timefunct
       integer                :: freq
   end type
   

   interface delete_source_
       module procedure source_null; end interface

   
   contains

!   --------------------------------------------------------------------------------
!   destructor
!	--------------------------------------------------------------------------------
	subroutine source_null (this)
    implicit none
	
	type(source_cl), intent(inout) :: this
	
	   this % name   = ' '
	   this % np     = 0
	   this % pmass  = 0.d0
	   this % zone   = 0
	   this % specie = 0
	   this % freq   = 1
	   this % TimeStartInj = 0.d0
	   this % TimeStopInj  = 0.d0
	   
	   nullify(this % par % xinj)
	   nullify(this % par % yinj)
	   nullify(this % par % zinj)
	   nullify(this % par % zbot)
	   nullify(this % par % ztop)
	   nullify(this % par % rcyr)
	   nullify(this % par % rcp )
	   nullify(this % par % xdist)
	   nullify(this % par % width )
	   nullify(this % par % height)
	   nullify(this % par % xinj_1)
	   nullify(this % par % yinj_1)
	   nullify(this % par % zinj_1)
	   nullify(this % par % xinj_2)
	   nullify(this % par % yinj_2)
	   nullify(this % par % zinj_2)
	   nullify(this % par % idwn  )
	   nullify(this % par % jdwn  )
	   nullify(this % par % kdwn  )
	   nullify(this % par % iup   )
	   nullify(this % par % jup   )
	   nullify(this % par % kup   )
	   nullify(this % par % np11x )
	   nullify(this % par % np11y )
	   nullify(this % par % np11z )
	   nullify(this % par % file  )
	   nullify(this % par % ix    )
	   nullify(this % par % iy    )
	   nullify(this % par % iz    )	
       
    end subroutine

!   --------------------------------------------------------------------------------------
!   reader
!	-------------------------------------------------------------------------------------- 			     

    subroutine read_source_ (this,namein,TypeInjIn,fname)
	   use gslib, only: upper_case, open_fname, open_fname_normal
	   use global_variables, only: UNEST
	   implicit none
	   type(source_cl),          intent(inout) :: this
	   character(len=*),         intent(in) :: namein,TypeInjIn
	   character(len=*),         intent(in) :: fname
	   character(len=len(trim(adjustl(namein))))    :: name
	   character(len=len(trim(adjustl(TypeInjIn)))) :: name2
	   character(len=100)                   :: string
	   real*8                               :: xinj,yinj,zinj,zbot,ztop,rcyr,rcp,xdist,width,height
	   real*8                               :: xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2,pmass,totmass
	   integer                              :: zone,specie
       integer                              :: idwn,jdwn,kdwn,iup,jup,kup,np
	   real*8                               :: np11x,np11y,np11z
       integer                              :: unit,unit2
	   logical                              :: existeix
       character(len=100)                   :: file
       real*8                               :: const
       integer                              :: flag,freq
       logical                              :: ReadFromFile
       
       integer,allocatable                  :: ix(:),iy(:),iz(:)

           call source_null(this)

           call open_fname (fname,unit)

           name  = trim(adjustl(namein))
           name2 = trim(adjustl(TypeInjIn))

		   name  = upper_case (name)
		   name2 = upper_case (name2)

		   this % name = name
		   this % TypeInj = name2	   
		   
		   np    = UNEST
		   
           ReadFromFile = .FALSE.
		   
		   select case (trim(adjustl(name)))

		    case ('BLOCK' )
			  read(unit,*) pmass,zone,specie
			  read(unit,*) idwn,jdwn,kdwn,iup,jup,kup !,np11x,np11y,np11z 
	               allocate(this % par % idwn); this % par % idwn  = idwn
	               allocate(this % par % jdwn); this % par % jdwn  = jdwn
	               allocate(this % par % kdwn); this % par % kdwn  = kdwn
                   allocate(this % par % iup);  this % par % iup   = iup
                   allocate(this % par % jup);  this % par % jup   = jup
                   allocate(this % par % kup);  this % par % kup   = kup
                   !allocate(this % np11x);  this % np11x = np11x
                   !allocate(this % np11y);  this % np11y = np11y
                   !allocate(this % np11z);  this % np11z = np11z
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie
 
		   case ( 'POINT' )        
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj,yinj,zinj 
			  	   allocate(this % par % xinj); this % par % xinj = xinj
	               allocate(this % par % yinj); this % par % yinj = yinj
	               allocate(this % par % zinj); this % par % zinj = zinj
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie

		   case ( 'LINE' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj,yinj,zbot,ztop 
	               allocate(this % par % xinj); this % par % xinj = xinj  
	               allocate(this % par % yinj); this % par % yinj = yinj
	               allocate(this % par % zbot); this % par % zbot = zinj 
                   allocate(this % par % ztop); this % par % ztop = ztop  
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie
		             
		   case ( 'CIRCLE' )        
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj,yinj,zbot,ztop,rcyr
	               allocate(this % par % xinj); this % par % xinj = xinj
	               allocate(this % par % yinj); this % par % yinj = yinj
	               allocate(this % par % zbot); this % par % zbot = zbot
                   allocate(this % par % ztop); this % par % ztop = ztop
                   allocate(this % par % rcyr); this % par % rcyr = rcyr
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie
	       
		   case ( 'RADIAL' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj,yinj,zbot,ztop,rcp 
  	               allocate (this % par % xinj);this % par % xinj = xinj  
	               allocate (this % par % yinj);this % par % yinj = yinj 
	               allocate (this % par % zbot);this % par % zbot = zbot
                   allocate (this % par % ztop);this % par % ztop = ztop
                   allocate (this % par % rcp); this % par % rcp  = rcp
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie
	       
		   case ( 'PLANE_RANDOM' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xdist,width,height
	               allocate (this % par % xdist) ; this % par % xdist  = xdist 
	               allocate (this % par % width) ; this % par % width  = width
	               allocate (this % par % height); this % par % height = height
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie
	       
		   case ( 'PLANE'  )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xdist,width,height
	               allocate (this % par % xdist) ; this % par % xdist  = xdist 
	               allocate (this % par % width) ; this % par % width  = width
	               allocate (this % par % height); this % par % height = height
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie
		   
		   case ( 'LINE_BY_POINTS' )			  
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2
	  	           allocate (this % par % xinj_1); this % par % xinj_1 = xinj_1   
	               allocate (this % par % yinj_1); this % par % yinj_1 = yinj_1
		           allocate (this % par % zinj_1); this % par % zinj_1 = zinj_1  
                   allocate (this % par % xinj_2); this % par % xinj_2 = xinj_2 
                   allocate (this % par % yinj_2); this % par % yinj_2 = yinj_2  
                   allocate (this % par % zinj_2); this % par % zinj_2 = zinj_2
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie

		   case ( 'LINE_BY_POINTS_RANDOM' )			  
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2
	  	           allocate (this % par % xinj_1); this % par % xinj_1 = xinj_1   
	               allocate (this % par % yinj_1); this % par % yinj_1 = yinj_1
		           allocate (this % par % zinj_1); this % par % zinj_1 = zinj_1  
                   allocate (this % par % xinj_2); this % par % xinj_2 = xinj_2 
                   allocate (this % par % yinj_2); this % par % yinj_2 = yinj_2  
                   allocate (this % par % zinj_2); this % par % zinj_2 = zinj_2
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie

		   case ( 'READ_PARTICLE_FILE' )
		       if (this%TypeInj /= 'DIRAC') then
		           stop '...reading from concentration file only accessible for GENERAL injections'
		       end if
		       ReadFromFile = .TRUE.
		       read(unit,*) file
			   inquire (file=file,exist=existeix)
			   if (.not. existeix) then
			            print *, '**** file does not exist: ',file
                        stop; end if
	           call open_fname_normal ( file, unit2 )
	           read(unit2,*)
			   read(unit2,*) np
		       allocate (this % par % file)
		       this % par % file = file
			   this % np     = np
			   this % pmass  =  UNEST
	           this % zone   = UNEST
	           this % specie = UNEST			   
               close(unit2)

		   case ( 'READ_CONCENTRATION_FILE' )
		       if (this%TypeInj /= 'DIRAC') then
		           stop '...reading from concentration file only accessible for DIRAC injections'
		       end if
		       ReadFromFile = .TRUE.
		       read(unit,*) file,const
			   inquire (file=file,exist=existeix)
			   if (.not. existeix) then
			            print *, '**** file does not exist: ',file
                        stop; end if
   			   read(unit,*) pmass,zone,specie
		       allocate (this % par % file)
		       allocate (this % par % const)
		       this % par % file = file
		       this % par % const = const
			   this % np     = UNEST
			   this % pmass  = pmass
	           this % zone   = zone
	           this % specie = specie

               
		   case ( 'LINE_FLUX_WEIGHTED' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2
	  	           allocate (this % par % xinj_1); this % par % xinj_1 = xinj_1   
	               allocate (this % par % yinj_1); this % par % yinj_1 = yinj_1
		           allocate (this % par % zinj_1); this % par % zinj_1 = zinj_1  
                   allocate (this % par % xinj_2); this % par % xinj_2 = xinj_2 
                   allocate (this % par % yinj_2); this % par % yinj_2 = yinj_2  
                   allocate (this % par % zinj_2); this % par % zinj_2 = zinj_2
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie

		   case ( 'VERTICAL_LINE_FLUX_WEIGHTED' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) idwn,jdwn,kdwn,kup  
	               allocate(this % par % idwn);  this % par % idwn  = idwn
	               allocate(this % par % jdwn);  this % par % jdwn  = jdwn
	               allocate(this % par % kdwn);  this % par % kdwn  = kdwn
                   allocate(this % par % kup);   this % par % kup   = kup
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie

		   case ( 'VERTICAL_BLOCK_FLUX_WEIGHTED' )
			  read(unit,*) pmass,zone,specie
			  read(unit,*) idwn,jdwn,kdwn,kup  
	               allocate(this % par % idwn);  this % par % idwn  = idwn
	               allocate(this % par % jdwn);  this % par % jdwn  = jdwn
	               allocate(this % par % kdwn);  this % par % kdwn  = kdwn
                   allocate(this % par % kup);   this % par % kup   = kup
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie

		   case ( 'CELLS_FILE_FLUX_WEIGHTED' )
			  read(unit,*) pmass,zone,specie
			  read(unit,*) file
			  call read_bocks_ (file,ix,iy,iz)
			       allocate(this % par % ix(size(ix))); this % par % ix = ix
			       allocate(this % par % iy(size(iy))); this % par % iy = iy
			       allocate(this % par % iz(size(iz))); this % par % iz = iz
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie

	       case default
	       
		      stop '**** could not match name of injection ****'

           end select
           
           !....read time function
           
           if (this%TypeInj == 'DIRAC') then             
              read(unit,*) this%TimeStartInj
              if(.not.ReadFromFile) read(unit,*) this%np
              this%TimeStopInj = this%TimeStartInj 
           end if
           if (this%typeInj == 'GENERAL') then
              read(unit,*) file,const
              call read_timefunction_ (file,const,1,this%timefunct)
              read(unit,*) freq
              this % freq = freq
              this % TimeStartInj = this%timefunct%time(1)
              this % TimeStopInj  = this%timefunct%time(this%timefunct%nt)
           end if
           

	end subroutine
!   --------------------------------------------------------------------------------------
!   printer
!	-------------------------------------------------------------------------------------- 			     

    subroutine read_bocks_ (fname,ix,iy,iz)
        use gslib, only: open_fname_normal
        implicit none
        character(len=*),      intent(in)    :: fname
        integer, allocatable,  intent(inout) :: ix(:),iy(:),iz(:)
        integer                              :: nblock,ivar,i,iunit

        call open_fname_normal(fname,iunit)
        read(iunit,*) nblock
        allocate(ix(nblock),iy(nblock),iz(nblock))
        read(iunit,*) ivar
        do i=1,ivar
            read(iunit,*)
        end do
        do i=1,nblock
            read(iunit,*) ix(i),iy(i),iz(i)
        end do
        close(iunit)
    end subroutine


    subroutine print_source_ (this,fname)
	   use gslib, only: upper_case, open_fname
	   implicit none
	   type(source_cl),       intent(inout) :: this
	   character(len=*),      intent(in)    :: fname
	   character(len=len(trim(adjustl(this%name)))) :: name
	   real*8                               :: xinj,yinj,zinj,zbot,ztop,rcyr,rcp,xdist,width,height
	   real*8                               :: xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2
	   integer                              :: zone,specie
       integer                              :: idwn,jdwn,kdwn,iup,jup,kup
       integer                              :: unit


           call open_fname (trim(adjustl(fname)),unit)

		   write(unit,'(a26,x,i5)')    'Mobile/Immobile Zone....: ',this%zone
		   write(unit,'(a26,x,i5,/)')  'Specie..................: ',this%specie


           name = trim(adjustl(this % name))
		   
		   if ( trim(adjustl(name)) == 'BLOCK' ) then
                   write(unit,'(/a26,x,i10)')  'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,11)  'type..: ',  this%name,  &
				                   'i-down: ',  this % par%idwn,  &
								   'j-down: ',  this % par%jdwn,  &
								   'k-down: ',  this % par%kdwn,  &
								   'i-up..: ',  this % par%iup,   &
								   'j-up..: ',  this % par%jup,   &
								   'k-up..: ',  this % par%kup    
               11  format (a8,a15,x,6(a8,i5,x)/)

		   elseif ( trim(adjustl(name)) == 'POINT' ) then
                   write(unit,'(/a26,x,i10)')   'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)')  'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,10) 'type: ',this%name,  &
				                  'xinj: ',this % par%xinj,  &
								  'yinj: ',this % par%yinj,  &
								  'zinj: ',this % par%zinj
               10  format (a6,a15,x,3(a6,g15.6,x)/)                	
 
		   elseif (trim(adjustl(name)) == 'LINE' ) then	          
                   write(unit,'(/a26,x,i10)')   'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)')  'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,12) 'type: ',this%name,  &
				                  'xinj: ',this % par%xinj,  &
								  'yinj: ',this % par%yinj,  &
								  'zbot: ',this % par%zbot,  &
								  'ztop: ',this % par%ztop 
               12  format (a6,a15,x,4(a6,g15.6,x)/)		             
		   elseif (trim(adjustl(name)) == 'CIRCLE' ) then        
                   write(unit,'(/a26,x,i10)')   'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,13) 'type: ',this%name,  &
				                  'xinj: ',this % par%xinj,  &
								  'yinj: ',this % par%yinj,  &
								  'zbot: ',this % par%zbot,  &
								  'ztop: ',this % par%ztop,  &
								  'rcyr: ',this % par%rcyr
               13  format (a6,a15,x,5(a6,g15.6,x)/)
	       
		   elseif ( trim(adjustl(name)) == 'RADIAL' ) then	          
                   write(unit,'(/a26,x,i10)')   'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,14)	'type: ',this%name,  &
				                    'xinj: ',this % par%xinj,  &
									'yinj: ',this % par%yinj,  &
									'zbot: ',this % par%zbot,  &
									'ztop: ',this % par%ztop,  &
									'rcp:  ',this % par%rcp
               14  format(a6,a15,x,5(a6,g15.6,x)/)
	       
		   elseif ( trim(adjustl(name)) == 'PLANE_RANDOM' ) then	          
                   write(unit,'(/a26,x,i10)')   'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,15)	'type: ',   this % name,  &
				                    'xdist: ',  this % par%xdist, &
									'width: ',  this % par%width, &
									'height: ', this % par%height
               15  format(a6,a,x,3(a7,g15.6,x)/)

		   elseif ( trim(adjustl(name)) == 'PLANE' ) then	          
                   write(unit,'(/a26,x,i10)')   'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,16)	'type: ',   this % name,  &
				                    'xdist: ',  this % par%xdist, &
									'width: ',  this % par%width, &
									'height: ', this % par%height
               16  format(a6,a15,x,3(a7,g15.6,x)/)
		   
		   elseif ( trim(adjustl(name)) == 'LINE_BY_POINTS' ) then			  
                   write(unit,'(/a26,x,i10)')   'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
                   write(unit,17) 'type: ',this % name,  &
				                  'x1: ',this % par%xinj_1,  &
								  'y1: ',this % par%yinj_1,  &
								  'z1: ',this % par%zinj_1,  &
                                  'x2: ',this % par%xinj_2,  &
								  'y2: ',this % par%yinj_2,  &
								  'z2: ',this % par%zinj_2
               17   format(a6,a15,x,3(a7,g15.6,x),3(a7,g15.6,x)/)

	       elseif ( trim(adjustl(name)) == 'READ_PARTICLE_FILE' ) then
		            write(unit,*) 
					write(unit,'(a40,a)') ' Reading from Particle File: ', trim(adjustl(this%par%file)) 
				    write(unit,*)
	       
	       elseif ( trim(adjustl(name)) == 'READ_CONCENTRATION_FILE' ) then
		            write(unit,*) 
					write(unit,'(a40,a)') ' Reading from Concentration File: ', trim(adjustl(this%par%file)) 
				    write(unit,*)

		   elseif ( trim(adjustl(name)) == 'LINE_FLUX_WEIGHTED' ) then	          
                   write(unit,'(/a26,x,i10)')   'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
                   write(unit,17) 'type: ',this%name,  &
				                  'x1: ',this % par%xinj_1,  &
								  'y1: ',this % par%yinj_1,  &
								  'z1: ',this % par%zinj_1,  &
                                  'x2: ',this % par%xinj_2,  &
								  'y2: ',this % par%yinj_2,  &
								  'z2: ',this % par%zinj_2
								  
		   elseif ( trim(adjustl(name)) == 'VERTICAL_LINE_FLUX_WEIGHTED' ) then
                   write(unit,'(/a26,x,i10)')  'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,11)  'type..: ',  this%name,  &
				                   'i-down: ',  this % par%idwn,  &
								   'j-down: ',  this % par%jdwn,  &
								   'k-down: ',  this % par%kdwn,  &
								   'k-up..: ',  this % par%kup

		   elseif ( trim(adjustl(name)) == 'VERTICAL_BLOCK_FLUX_WEIGHTED' ) then
                   write(unit,'(/a26,x,i10)')  'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,11)  'type..: ',  this%name,  &
				                   'i-down: ',  this % par%idwn,  &
								   'j-down: ',  this % par%jdwn,  &
								   'k-down: ',  this % par%kdwn,  &
								   'k-up..: ',  this % par%kup

		   elseif ( trim(adjustl(name)) == 'CELLS_FILE_FLUX_WEIGHTED' ) then
                   write(unit,'(/a26,x,i10)')  'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
	               write(unit,11)  'type..: ',  this%name
				                   !'i-down: ',  this % par%idwn,  &
								   !'j-down: ',  this % par%jdwn,  &
								   !'k-down: ',  this % par%kdwn,  &
								   !'k-up..: ',  this % par%kup


		   elseif ( trim(adjustl(name)) == 'LINE_BY_POINTS_RANDOM' ) then			  
                   write(unit,'(/a26,x,i10)')   'Number Particles........: ',this%np
		           write(unit,'(a26,x,g16.6)') 'Particle Mass...........: ',this%pmass
				   write(unit,*)
                   write(unit,17) 'type: ',this%name,  &
				                  'x1: ',this % par%xinj_1,  &
								  'y1: ',this % par%yinj_1,  &
								  'z1: ',this % par%zinj_1,  &
                                  'x2: ',this % par%xinj_2,  &
								  'y2: ',this % par%yinj_2,  &
								  'z2: ',this % par%zinj_2

		   else
		      stop '**** could not match name of injection ****'

           end if
 
	end subroutine
    



   end module

!********************************************************************************************
!  VECTOR OF SOURCES CLASS 
!********************************************************************************************
 module source_vect_class
 use source_class

 implicit none

 public !everything is public (inheritance of source_class)
 public :: source_vect_cl           !class
 public ::                                    & !methods
           alloc_source_vect_               , &
           delete_source_vect_            !  , &
          ! organize_source_by_species_
           
 type source_vect_cl
      type(source_cl),pointer :: num(:) => null()
 end type

 contains

! subroutine organize_source_by_species_ (this)
!    use source_class
!    use global_variables, only: nspecie
!    use gslib,            only: sortem
!    implicit none
!    type(source_vect_cl), intent(inout) :: this
!    type(source_vect_cl)                :: dummy
!    integer                             :: i,nsource
!    real*8                              :: sp(size(this%num)),b(size(this%num))
!    
!    sp=-999.d0
!    
!    nsource=size(this%num)
!    
!    call alloc_source_vect_ (dummy,nsource)
!   
!    do i=1,nsource
!       dummy%num(i) = copy_source_ (this%num(i))
!       sp(i) = this%num(i)%specie
!       b(i)  = i
!    end do
!    
!    call sortem(1,nsource,sp,1,b,b,b,b,b,b,b)  
!
!    do i=1,nsource
!       this%num(i) = copy_source_ (dummy%num(int(b(i))))
!    end do
!
!    call delete_source_vect_ (dummy)
! 
! end subroutine

 subroutine alloc_source_vect_ (this,n)
    implicit none
    type(source_vect_cl), intent(inout) :: this
    integer,              intent(in)    :: n
	allocate(this%num(n))
 end subroutine

 subroutine delete_source_vect_ (this)
    use source_class
	implicit none
    type(source_vect_cl), intent(inout) :: this
	integer                             :: i,n
	if (associated(this%num)) then
	n = size(this%num)
	do i=1,n
	   call delete_source_ (this%num(i)); end do
	end if
 end subroutine

 end module