program f_call_c
  implicit none


  interface
     subroutine cfun(x,len) bind( c )
       use,intrinsic :: iso_c_binding
       implicit none
       integer( c_int) :: len
       real(c_double) :: x(0:3,0:len)
     end subroutine   cfun

     subroutine vec(r,len) bind(c)
       use,intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: len
       real(c_double) :: r(0:len)
     end subroutine vec


  end interface


  double precision, allocatable :: x(:,:),r(:)
  integer :: len,row,i,j

  len = 7
  allocate(x(0:3,0:len))
  allocate(r(0:len))

  do i =0,len
     r(i) = i
  enddo


  do i = 0,len
     do j = 0,3
        x(j,i) = j+i
     enddo
  enddo

  call vec(r,%val(len) )

  row = 3
  call cfun(x,%val(len))
end program f_call_c

