  module massProfile_class
   use geometry_class
   implicit none

   private
   public :: massProfile_cl            !class
   public :: print_massProfile_  ! , &       !methods
                
    type massProfile_cl
	   type(geometry_cl), pointer :: geo => null()
	end type 

    contains

    subroutine print_massProfile_ (this,fname)
	     use gslib, only: open_fname
		 use array_class
		 use global_variables, only:fdbg
		 implicit none
		 type(massProfile_cl),       intent(in) :: this
		 character(len=*), optional, intent(in) :: fname
		 integer                                :: unit

         if (.not.present(fname)) then; unit = 6
		                          else; call open_fname (fname,unit); endif

		 write(unit,*)
		 if (associated(this % geo)) then
              write(unit,1) 'Mass Profile Nx...: ',this%geo%nx
	          write(unit,1) 'Mass Profile Ny...: ',this%geo%ny
	          write(unit,1) 'Mass Profile Nz...: ',this%geo%nz
		    1 format (x,a20,x,i7)
			  write(unit,*)
	          call list_array_ (this%geo%dx,fdbg,'Mass Profile Dx')
		      call list_array_ (this%geo%dy,fdbg,'Mass Profile Dy')
		      call list_array_ (this%geo%dz,fdbg,'Mass Profile Dz')
		 end if
		 write(unit,*)		  
	
	end subroutine    
    

  end module