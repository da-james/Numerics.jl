!==============================================================================
! Author: David James, 20200414
! sources:
! 1. Numerical Recipes; 2 ed.; Press, W.
! 2. Numerical Methods for Fluid Dynamics; 2 ed.; Durran, D.
! 3. Numerical Analysis; 10 ed.; Burden, R.
!------------------------------------------------------------------------------
! module: module_ode
! - subroutine rkdriver: driver for rk4 (1)
! - subroutine rk4: 4th-Order Runge-Kutta Method (1)
! - subroutine fdm: Finite Difference Method (3)
!==============================================================================
module module_ode
  use module_values, only : dp, i16, NMAX, NSTPMX
  implicit none

contains

  subroutine rkdriver(vstart,nvar,x1,x2,nstep,derivs, xx, y)
    implicit none
    integer(i16), intent(in) :: nvar, nstep
    real(dp), intent(in) :: x1, x2
    real(dp), dimension(nvar), intent(in) :: vstart
    real(dp), dimension(NSTPMX), intent(out) :: xx
    real(dp), dimension(NMAX,NSTPMX), intent(out) :: y
    external :: derivs
    !--------------------------------------------------------------------------
    !
    ! subroutine rkdriver: driver for rk4 (1)
    ! A driver routine to run the rk4 sequentionally.
    ! Refer to Numerical Recipes(pg. 707-708) for
    ! information about the subroutine.
    ! @param:
    !   vstart - real, array, size:nvar, initial condition
    !   nvar - int, number of dependent variables
    !   x1 - real, intial step
    !   x2 - real, final step
    !   nstep - int, number of steps rk4 will run
    !   derivs - subroutine, coupled equations, derivs(
    !          x - real, step
    !          v - real, array, size of NMAX, input conditions
    !          dv - real, array, size of NMAX, output values)
    !   xx - real, array(size:NSTPMX), output steps
    !   y - real, array(r:NMAX,c:NSTPMX), output values
    !
    !--------------------------------------------------------------------------
    integer(i16) :: i,k
    real(dp) :: h, x
    real(dp), dimension(NMAX) :: dv, v

    do i=1,nvar
       v(i) = vstart(i)
       y(i,1) = v(i)
    end do

    xx(1) = x1
    x = x1
    h = (x2-x1)/nstep

    do k=1,nstep
       call derivs(x,v,dv)
       call rk4(v,dv,nvar,x,h,v,derivs)

       x = x+h
       xx(k+1) = x
       do i=1,nvar
          y(i,k+1) = v(i)
       end do
    end do
  end subroutine rkdriver

  subroutine rk4(y, dydx, n, x, h, yout, derivs)
    implicit none
    integer(i16), intent(in) :: n
    real(dp), intent(in) :: h, x
    real(dp), dimension(n), intent(in) :: dydx, y
    real(dp), dimension(n), intent(out) :: yout
    external :: derivs
    !-------------------------------------------------------------------------
    !
    ! subroutine rk4: 4th-Order Runge-Kutta Method (1)
    ! This will solve out an ode numerically. Refer to
    ! Numerical Recipes (pg. 706-707) for information
    ! on the subroutine
    ! @param:
    !   y - real, array, initial start condition
    !   dydx - real, array, next step condition
    !   n - int, size of array
    !   x - real, independent variable
    !   h - real, step
    !   yout - real, array, next advanced solution
    !   derivs - subroutine, coupled equations, derivs(
    !          x - real, step
    !          v - real, array, size of NMAX, input conditions
    !          dv - real, array, size of NMAX, output values)
    !
    !-------------------------------------------------------------------------
    integer(i16) :: i
    real(dp), dimension(NMAX) :: dym, dyt, yt
    real(dp) :: h6, hh, xh

    hh = h*0.5
    h6 = h/6.
    xh = x+hh

    do i=1,n
       yt(i) = y(i) + hh*dydx(i)
    end do

    call derivs(xh,yt,dyt)

    do i=1,n
       yt(i) = y(i) + hh*dyt(i)
    end do

    call derivs(xh,yt,dym)

    do i=1,n
       yt(i) = y(i) + h*dym(i)
       dym(i) = dyt(i) + dym(i)
    end do

    call derivs(x+h,yt,dyt)

    do i=1,n
       yout(i) = y(i) + h6*(dydx(i) + dyt(i) + 2.*dym(i))
    end do
  end subroutine rk4

  subroutine fdm(f,v,t,y0,n,h,x,y)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: h
    real, dimension(2), intent(in) :: v, t ,y0
    real, dimension(n), intent(out) :: x, y
    external :: f
    !-------------------------------------------------------------------------
    !
    ! subroutine fdm: Finite Difference Method (3)
    ! This will solve out a 2nd order DE using the Finite Difference
    ! approach.
    ! Let D = d/dx then the equation should have the form
    ! (D^2)y + p*D(y) + q*y = f(x) where p,q are constants and f(x) is a
    ! function, and y(t0) = y0 and y(t1) = y1 such that t0 <= t <= t1.
    ! @param:
    !  f -  sub, external, (x,y): f(x) evaluates values and returns y
    !  v - real, in, dim(2)     : p, q as v(1) and v(2), respectively
    !  t - real, in, dim(2)     : t0, t1 as t(1) and t(2), respectively
    ! y0 - real, in, dim(2)     : y0, y1 as y0(1) and y0(2), respectively
    !  n -  int, in             : number of steps in t0 <= t <= t1
    !  h - real, in             : step size
    !  x - real, out, dim(n)    : x elements along t0 <= t <= t1
    !  y - real, out, dim(n)    : y elements evaluated on t0 <= t <= t1
    !
    !-------------------------------------------------------------------------
    real :: s1, s2, s3
    real, dimension(:), allocatable :: L, D, U, f_b
    integer :: i, m

    s1 = 1/h**2 - v(1)/(2*h)
    s2 = v(2) - 2/h**2
    s3 = 1/h**2 + v(1)/(2*h)

    m = n - 2
    allocate(L(m),D(m),U(m),f_b(m))
    L(:) = s1
    D(:) = s2
    U(:) = s3

    do i=1,m
       call f(x(i+1),f_b(i))
    enddo
    f_b(1) = f_b(1) - s1*y0(1)
    f_b(m) = f_b(n) - s3*y0(2)

    call tridiag(m,L,D,U,f_b,y)
    y(1) = y0(1)
    y(2:) = f_b
    y(n) = y0(2)

  end subroutine fdm

  subroutine tridiag(jmx,a,b,c,f,q)
    integer, intent(in) :: jmx
    real, dimension(jmx), intent(inout) :: a, b, c, f, q
    !-------------------------------------------------------------------------
    !
    ! subroutine tridiag: Thomas Algorithm
    ! This subroutine will solve out a tridiagonal matrix system of equations
    ! using the Thomas Algorithm. Where the system is T*x = b, where T is a
    ! matrix made up of U, D, and L, the upper diagonal, diagonal, and lower
    ! diagonal, respectively.
    ! @param:
    ! jmx -  int,    in          : size of the diagonal of the matrix
    !   a - real, inout, dim(jmx): lower diagonal
    !   b - real, inout, dim(jmx): diagonal
    !   c - real, inout, dim(jmx): upper diagonal
    !   f - real, inout, dim(jmx): constant terms
    !   q - real, inout, dim(jmx): working array
    !
    !-------------------------------------------------------------------------
    real :: p       ! working variable for performing algorithm
    integer :: j    ! index variable for loop

    c(jmx) = 0

    ! forward substitution
    q(1) = -c(1) / b(1)
    f(1) = f(1) / b(1)

    do j=2,jmx
       p = 1.0 / ( b(j) + a(j)*q(j-1) )
       q(j) = -c(j)*p
       f(j) = ( f(j) - a(j)*f(j-1) )*p
    enddo

    ! back substitution
    do j=jmx-1,1,-1
       f(j) = f(j) + q(j)*f(j+1)
    enddo
  end subroutine tridiag

  subroutine linspace(a,b,array)
    real, intent(in) :: a, b
    real, dimension(:), intent(inout) :: array

    integer :: n, i
    real :: range

    n = size(array)
    range = b - a

    if(n == 0) then
       return
    elseif(n == 1) then
       array(1) = a
       return
    else
       do i=1, n
          array(i) = a + range*(i-1)/(n-1)
       enddo
    endif
  end subroutine linspace

end module module_ode
