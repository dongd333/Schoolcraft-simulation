
!****************************************************************************
!   WELL CLASS
!****************************************************************************
    module well_class
	use breakthru_class
	implicit none

	private
	public :: well_cl     ! class
	public ::                          & ! methods
              print_well_            , &
              distance_to_well_      , &
			  print_well_moments_    , &
			  read_assign_well_      , &
			  get_cell_well_                        

	type well_cl
		 character(len=20)              :: name
		 real*8                         :: xwell,ywell,rwell
		 real*8, pointer                :: zwbot => null()     ! bottom of a well
		 real*8, pointer                :: zwtop => null()     ! top of well
         integer                        :: col,row,botlay,toplay
		 integer                        :: out                 !1 = particle is removed when crossing
		 type(breakthru_cl), pointer    :: btc(:)              !arrival time
		 logical                        :: SaveBTC =.TRUE.
		 logical                        :: recirculation =.FALSE.
	end type well_cl


    interface print_well_moments_
	    module procedure print_well_moments_,print_well_moments_1_; end interface

    interface print_well_
	     module procedure print_well_1_,print_well_2_; end interface


	contains
   

    subroutine read_assign_well_ (this,unit)  !read well parameters from file
         use gslib, only: count_words,separate_words,char_to_dble,char_to_integer,upper_case 
		 implicit none
		 type(well_cl), intent(inout)   :: this
		 integer,       intent(in)      :: unit
         character(len=20)              :: name
         integer                        :: flag,botlay,toplay
		 real*8                         :: xw,yw,rw,zbot,ztop
		 logical                        :: SaveBTC
		 zbot = 1.d0
		 ztop = 0.d0
		 read(unit,*) name,xw,yw,rw,zbot,ztop,flag,SaveBTC 
	        name = upper_case (name)
	        this%name  = name
	        this%xwell = xw
   		    this%ywell = yw
		    this%rwell = rw
		    this%out   = flag
		    this%SaveBTC = SaveBTC
			if (ztop > zbot) then			   
               if (.not.associated(this%zwbot)) allocate (this%zwbot); this%zwbot = zbot 
		       if (.not.associated(this%zwtop)) allocate (this%zwtop); this%zwtop = ztop 
            end if
            
            !write(10002, '(A15, 4F9.4)') this%name, ztop, zbot, this%xwell, this%ywell
	end subroutine


    subroutine get_cell_well_ (well,geo)
         use gslib, only: locate
         use array_class
         use geometry_class
         implicit none
		 type(well_cl), intent(inout)   :: well
         type(geometry_cl),  intent(in) :: geo
         integer                        :: i,j,k2,k3,nx,ny,nz
   
            nx = geo%nx
            ny = geo%ny
            nz = geo%nz
            
           if (associated(well%zwtop)) then
			    if (nz == 1) then
				   k2  = 1
                   k3  = 1
			    else
				   call locate(geo%zmesh%values,nz+1,1,nz+1,well%zwbot,k2)
			       call locate(geo%zmesh%values,nz+1,1,nz+1,well%zwtop,k3)
                end if     
            else
                k2 = 1
                k3 = nz
            end if
            
            well%toplay = k3
            well%botlay = k2

			if (nx == 1) then
				i = 1
			else
				call locate(geo%xmesh%values,nx+1,1,nx+1,well%xwell,i)
            end if     

			if (ny == 1) then
				j = 1
			else
				call locate(geo%ymesh%values,ny+1,1,ny+1,well%ywell,j)
            end if     

            well%col = i
            well%row = j
            !write(10002, '(15x, 4I9)') well%toplay, well%botlay, well%col, well%row
            
    end subroutine

    subroutine print_well_moments_ (name,fname)
         use breakthru_class
         use gslib,             only: open_fname
         use global_variables,  only: nspecie
		 implicit none
		 type(well_cl),             intent(inout) :: name
         character(len=*),optional, intent(in)    :: fname
		 real*8                                   :: mean(3),var(3)
		 real*8                                   :: alfa(4),mom(4)
	     real*8                                   :: skew,kurt,xwell,ywell,rwell
         integer                                  :: iunit,np,ispe
         do ispe=1,nspecie
		    if (name%btc(ispe)%np >5) then
            call open_fname(fname,iunit)
		    call moments_breakthru_ (name%btc(ispe),alfa,mom,skew,kurt)
		    np = name%btc(ispe)%np
		      xwell = name%xwell
              ywell = name%ywell
			  rwell = name%rwell
		      write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np, &
			                'well ','xwell= ',xwell,'ywell= ',ywell,'rwell= ',rwell      
			end if
	     end do
	 1     format (9(g21.14,1x),i20,x,a6,3(a8,g15.6,x))
     end subroutine

    subroutine print_well_moments_1_ (name,btc,fname)
         use breakthru_class
         use gslib,             only: open_fname
         use global_variables,  only: nspecie
		 implicit none
		 type(well_cl),             intent(inout) :: name
		 type(breakthru_cl),        intent(in)    :: btc(nspecie)
         character(len=*),optional, intent(in)    :: fname
		 real*8                                   :: mean(3),var(3)
		 real*8                                   :: alfa(4),mom(4)
	     real*8                                   :: skew,kurt,xwell,ywell,rwell
         integer                                  :: iunit,np,ispe
         if (btc(1)%np>5) then
            call open_fname(fname,iunit)
            do ispe=1,nspecie 
		        call moments_breakthru_ (btc(ispe),alfa,mom,skew,kurt)
		        np = btc(ispe)%np
		        xwell = name%xwell
                ywell = name%ywell
			    rwell = name%rwell
		        write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np, &
			                'well ','xwell= ',xwell,'ywell= ',ywell,'rwell= ',rwell
            end do     
         end if
	 1   format (9(g21.14,1x),i20,x,a6,3(a8,g15.6,x))
     end subroutine



    function distance_to_well_ (name,xp,yp,zp) result(value)
         implicit none
		 type(well_cl),    intent(in) :: name
		 real*8,           intent(in) :: xp,yp
		 real*8, optional, intent(in) :: zp
		 real*8                       :: value	 
		    value=dsqrt((name%xwell-xp)**2+(name%ywell-yp)**2) 
	end function

    subroutine print_well_1_ (name,fname)
         use gslib, only: open_fname
		 implicit none
		 type(well_cl), intent(in)             :: name
		 character(len=*),optional, intent(in) :: fname
		 integer                               :: iunit
 		   call open_fname(fname,iunit)
		   if (associated(name%zwbot)) then
                write(iunit,2) 'ZONE T= "Xw: ',name%xwell,'Yw: ',name%ywell,'Rw: ',name%rwell,   &
		                       'Zwbot: ',name%zwbot,'Zwtop: ',name%zwtop,'OUT: ',name%out,'"'
           elseif (.not.associated(name%zwbot)) then
                write(iunit,3) 'ZONE T= "Xw: ',name%xwell,'Yw: ',name%ywell,'Rw: ',name%rwell,   &
		                       'OUT: ',name%out,'"'
		   end if
         2 format(a13,g15.6,2x,2(a4,g15.6,2x),2(a7,g15.6,2x),a5,i2,a2)
		 3 format(a13,g15.6,2x,2(a4,g15.6,2x),a5,i2,a2)                                                       
		   close(iunit)
	end subroutine

    subroutine print_well_2_ (name,fname,ispe)
         use gslib, only: open_fname
		 implicit none
		 type(well_cl), intent(in)             :: name
		 integer,       intent(in)             :: ispe
		 character(len=*),optional, intent(in) :: fname
		 integer                               :: iunit
 		   call open_fname(fname,iunit)
		   if (associated(name%zwbot)) then
                write(iunit,2) 'ZONE T= "Xw: ',name%xwell,'Yw: ',name%ywell,'Rw: ',name%rwell,   &
		                       'Zwbot: ',name%zwbot,'Zwtop: ',name%zwtop,' Specie: ',ispe,'"'
           elseif (.not.associated(name%zwbot)) then
                write(iunit,3) 'ZONE T= "Xw: ',name%xwell,'Yw: ',name%ywell,'Rw: ',name%rwell,   &
		                       ' Specie: ',ispe,'"'
		   end if
         2 format(a13,g15.6,2x,2(a4,g15.6,2x),2(a7,g15.6,2x),a9,i2,a2)
		 3 format(a13,g15.6,2x,2(a4,g15.6,2x),a9,i2,a2)                                                       
		   close(iunit)
	end subroutine

    end module

!******************************************************************************
!  VECTOR OF WELLS  
!******************************************************************************
 module well_vect_class
 use well_class

 implicit none

 public    !everytging is public 
 public :: well_vect_cl                      !class
 public ::                                 & !methods
              alloc_well_vect_           , &
			  nwell_                     , &
			  print_well_moments_vect_         

 type well_vect_cl
      type(well_cl), pointer :: num(:) => null()
 end type



 interface print_well_moments_vect_
      module procedure print_well_moments_vect_,print_well_moments_vect_1_; end interface

 contains



 subroutine alloc_well_vect_ (name,n,nspecie)
    use well_class
    use breakthru_class
	implicit none
    type(well_vect_cl),   intent(inout) :: name
    integer,              intent(in)    :: n
    integer,              intent(in)    :: nspecie
	integer                             :: i,ispe
	   allocate(name%num(n))
	   do i=1,n
	      allocate (name%num(i)%btc(nspecie))
	      do ispe=1,nspecie
	       call initialize_btc_ (name%num(i)%btc(ispe))
	     end do
	   end do
 end subroutine


 subroutine print_well_moments_vect_ (name,fname)
    use well_class
	implicit none
	type(well_vect_cl),         intent(inout) :: name
	character(len=*), optional, intent(in)    :: fname
	integer                                   :: i,n
	  if (associated(name%num)) then
	  n = size(name%num)
      do i=1,n
         call print_well_moments_ (name%num(i),fname);  end do
	  end if
 end subroutine

 
 subroutine print_well_moments_vect_1_ (name,i,fname)
    use well_class
	implicit none
	type(well_vect_cl),         intent(inout) :: name
	integer,                    intent(in)    :: i
	character(len=*), optional, intent(in)    :: fname
	integer                                   :: n
	  if (associated(name%num)) then
	  n = size(name%num)
         call print_well_moments_ (name%num(i),fname)
	  end if
 end subroutine



 function nwell_ (name) result (n)
         implicit none
		 type(well_vect_cl), intent(in) :: name
         integer                      :: n
		   if (associated(name % num)) then
		       n = size(name % num)
		   else
		       n = 0
		   end if
 end function

 end module

