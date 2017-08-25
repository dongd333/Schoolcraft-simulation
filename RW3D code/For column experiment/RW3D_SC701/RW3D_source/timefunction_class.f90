

module timefunction_class
   
	implicit none

	private
	public :: timefunction_cl     ! class
	public ::                                             & ! methods
              read_timefunction_                          , &
			  evaluate_timefunction_interpolate_          , &
			  evaluate_timefunction_upwind_                    
			                
    type timefunction_cl    
         logical           :: action
         character(len=50) :: name
         integer           :: nt
         real*8            :: scale 
         real*8,  pointer  :: time(:)  => null()  !initial time
         real*8,  pointer  :: val(:)   => null()  !final time      
    end type      
    
    contains
    
    subroutine read_timefunction_ (fname,scale,flag,func) !initialize the empty list
       use gslib, only: open_fname_normal
       implicit none
       character(len=*),      intent(in)    :: fname
       real*8,                intent(in)    :: scale
       integer,               intent(in)    :: flag
       type(timefunction_cl), intent(inout) :: func
       integer                              :: iunit,nt,i
       real*8                               :: tt
       
       if ( flag == 0 ) then
          func%action = .FALSE.
          func%nt = 1
          func%scale = scale
          return
       end if
       
       if (flag == 1) then
          func%action =.TRUE.
          call open_fname_normal(fname,iunit)
          read(iunit,*) !heading
          read(iunit,*) func%name          
          read(iunit,*) nt
          func%nt = nt
          func%scale = scale
          allocate(func%time(nt),func%val(nt))
          do i=1,nt
            read(iunit,*) func%time(i),func%val(i)
          end do
          close(iunit)
          return
      end if
       
    end subroutine

    function evaluate_timefunction_interpolate_ (func,time) result(val)
       use gslib, only: locate
       implicit none
       type(timefunction_cl), intent(in) :: func
       real*8,                intent(in) :: time
       integer                           :: it,nt
       real*8                            :: val
       nt = func%nt
       if (.not.associated(func%time)) then
          val = func%scale
          return
       end if
       if (time < func%time(1).or. time > func%time(nt)) then
          val = 0.d0
          return
       end if
       call locate(func%time,nt,1,nt,time,it)
       val = func%val(it) + (func%val(it+1)-func%val(it))/(func%time(it+1)-func%time(it))*(time-func%time(it)) 
       val = val * func%scale 
    end function

   function evaluate_timefunction_upwind_ (func,time) result(val)
       use gslib, only: locate
       implicit none
       type(timefunction_cl), intent(in) :: func
       real*8,                intent(in) :: time
       integer                           :: it,nt, test
       real*8                            :: val
       nt = func%nt
       if (.not.associated(func%time)) then
          val = func%scale
          return
       end if
       if (time < func%time(1).or. time > func%time(nt)) then
          val = 0.d0
          return
       else if (time == func%time(1)) then
          val = func%val(1)
          return
       end if
       call locate(func%time,nt,1,nt,time,it)
       val = func%val(it) 
       val = val * func%scale 
    end function

end module timefunction_class