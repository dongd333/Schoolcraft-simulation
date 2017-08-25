
!*****************************************************************************************************
!   PLANE CLASS
!*****************************************************************************************************   
    module Inner_bdy_class
	implicit none

	private
	public :: Inner_bdy_cl     ! class
	public ::                               & ! methods
              assign_Inner_bdy_xy               , & 
			  read_assign_Inner_bdy_          
	type Inner_bdy_cl
		 real*8                         :: xinner_bdy=0.d0,yinner_bdy=0.d0            !Position of an Inner-boundary
		 character(len=1)               :: side ='N'                         !'x'=parallel x direction,'y'=parallel to y direction
		 integer                        :: ls = 1                            !1 = particle larger than the value is removed when crossing, 0 = less than -> remove
	end type Inner_bdy_cl


	contains
   
	

	subroutine assign_Inner_bdy_xy (name,a,typ,flag)
		 use gslib, only: upper_case
		 implicit none
		 type(Inner_bdy_cl),                intent(inout) :: name
		 real*8,                        intent(in)    :: a
		 character(len=*),    optional, intent(inout) :: typ
	     integer,             optional, intent(in)    :: flag
         
         name%side = typ 
			if (name%side == 'X') then 
			       name%xinner_bdy = a
				   name%yinner_bdy = 0.d0
			else if (name%side == 'Y') then
			       name%xinner_bdy = 0.d0
				   name%yinner_bdy = a
            end if
            name%ls = flag
            !write(*,*) a, t
	end subroutine  

    subroutine read_assign_Inner_bdy_ (name,iunit)
	     use gslib, only: ifcharacter
		 implicit none
		 type(Inner_bdy_cl),  intent(inout) :: name
		 integer,         intent(in)    :: iunit
		 logical                        :: exists,oldinput
		 real*8                         :: ibdy
		 character(len=1)               :: side
		 character(len=3)               :: flagchar
		 integer                        :: flag
              inquire (unit=iunit,exist=exists)
			  if(.not.exists) then
			     print *, 'file does not exist during read_assign_Inner_bdy'
				 stop
			  end if    
	              read(10,*) ibdy,side,flag;  call assign_Inner_bdy_xy (name,ibdy,side,flag) 
                  !write(*,*) ibdy,side,flag
                  
	end subroutine


    end module
    
    !******************************************************************************
!  VECTOR OF PLANES  
!******************************************************************************
 module Inner_bdy_vect_class
 use Inner_bdy_class

 implicit none

 public    !everything is public (inheritance of plane_class)
 public :: Inner_bdy_vect_cl                        !class
 public ::                                    & !methods
           alloc_Inner_bdy_vect_                , &
		   n_Inner_bdy_

 type Inner_bdy_vect_cl
      type(Inner_bdy_cl),pointer :: num(:) => null()
 end type

 contains

 subroutine alloc_Inner_bdy_vect_ (name,n)
    use Inner_bdy_class
	implicit none
    type(Inner_bdy_vect_cl), intent(inout) :: name
    integer,             intent(in)    :: n
	integer                            :: i
	      allocate(name%num(n))
 end subroutine

 function n_Inner_bdy_ (name) result (n)
         implicit none
		 type(Inner_bdy_vect_cl), intent(in) :: name
         integer                         :: n
		   if (associated(name%num)) then
		       n = size(name%num)
		   else
		       n = 0
		   end if
 end function



 end module