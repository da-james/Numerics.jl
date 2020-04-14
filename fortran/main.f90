program main
  use module_values, only : PI
  use module_ode
  implicit none

  real, dimension(2) :: v,t,y0
  real, dimension(:), allocatable :: x,y
  real :: h
  integer :: n, i

  interface
     subroutine nhde(x,f)
     end subroutine nhde
  end interface

  v(1) = -2
  v(2) = 1

  t(1) = 0
  t(2) = 2

  y0(1) = 0
  y0(2) = -4

  h = 0.2
  n = floor((t(2)-t(1))/h)

  allocate(x(1:n),y(1:n))
  call linspace(t(1),t(2),x)

  call fdm(nhde,v,t,y0,n,h,x,y)
  do i=1,n
     print *, x(i), y(i), 1/6*x(i)**3*exp(x(i)) - 5/3*x(i)*exp(x(i)) + 2*exp(x(i)) - x(i) - 2
  enddo

end program

subroutine nhde(x,f)
  implicit none
  real, intent(in) :: x
  real, intent(out) :: f

  ! f = -1*sin(3*PI*x/2)
  f = x*exp(x) - x
end subroutine nhde
