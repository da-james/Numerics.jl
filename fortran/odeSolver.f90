!==================================================
! Author: David James, 20190905
! background: Numerical Recipes, 2 ed., 92
! module: odeSolver
! - subroutine rkdriver: driver for rk4
! - subroutine rk4: 4th-Order Runge-Kutta Method
! - subroutine polint: polynomial interpolation
! - subroutine hunt: index finder
! - function smallerIndex: Array Shortener
!==================================================
module odeSolver
  implicit none
  integer, parameter, private :: ikind = selected_real_kind(p=15) ! precision for values
  integer, parameter, public :: NMAX = 10
  integer, parameter, public :: NSTPMX = 500

contains
  !==================================================
  ! subroutine rkdriver: driver for rk4
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
  !==================================================
  subroutine rkdriver(vstart,nvar,x1,x2,nstep,derivs, xx, y)
    integer(kind=4) :: nstep, nvar
    real(kind=ikind) :: x1, x2, vstart(nvar), xx(NSTPMX), y(NMAX,NSTPMX)
    external :: derivs

    integer(kind=4) :: i,k
    real(kind=ikind) :: h, x, dv(NMAX), v(NMAX)

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

  !==================================================
  ! subroutine rk4: 4th-Order Runge-Kutta Method
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
  !==================================================
  subroutine rk4(y, dydx, n, x, h, yout, derivs)
    implicit none
    integer(kind=4) :: i, n
    real(kind=ikind) :: h, x, dydx(n), y(n), yout(n)
    external :: derivs

    real(kind=ikind) :: h6, hh, xh, dym(NMAX), dyt(NMAX), yt(NMAX)

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

  !==================================================
  ! subroutine polint: polynomial interpolation
  ! Will interpolate values to a nth order polynomial.
  ! Refer to Numerical Recipes(pg. 103-104) for
  ! information about the subroutine.
  ! @param:
  !   xa - real, array, x values of function
  !   ya - real, array, y values of function
  !   n - int, size of array
  !   x - real, initial x value
  !   y - real, correlated y
  !   dy - real, error estimate
  !==================================================
  subroutine polint(xa,ya,n,x,y,dy)
    integer(kind=4), intent(in) ::  n
    real(kind=ikind), intent(in) :: xa(n), ya(n), x
    real(kind=ikind), intent(out) :: y, dy
    integer(kind=4) i, m, ns
    real(kind=ikind) den, dif, dift, ho, hp, w, c(NMAX), d(NMAX)

    ns = 1
    dif = abs(x-xa(1))

    ! we find the index ns of the closest table entry
    do i=1,n
       dift = abs(x-xa(i))
       if(dift < dif) then
          ns = i
          dif = dift
       end if

       ! initialize values for c(NMAX) and d(NMAX)
       c(i) = ya(i)
       d(i) = ya(i)
    end do

    ! initial approximation to y
    y = ya(ns)
    ns = ns - 1

    do m=1,n-1
       ! for each column of the table
       ! we loop over the current c's and d's and update
       do i=1,n-m
          ho = xa(i) - x
          hp = xa(i+m) - x
          w = c(i+1) - d(i)
          den = ho - hp

          den = w/den

          ! c's and d's are updated
          d(i) = hp*den
          c(i) = ho*den
       end do
       if(2*ns < n-m) then
          dy = c(ns+1)
       else
          dy = d(ns)
          ns = ns - 1
       end if
       y = y + dy
    end do
  end subroutine polint

  !==================================================
  ! subroutine hunt: index finder
  ! Will find the index of the desired x value given
  ! in an array
  ! Refer to Numerical Recipes(pg. 112-113) for
  ! information about the subroutine.
  ! @param:
  !   xx - real, array, x values of function
  !   n - int, size of array
  !   x - real, initial x value
  !   jlo - int, index of x such that x is between
  !              xx(jlo) and xx(jlo+1)
  !==================================================
  subroutine hunt(xx,n,x,jlo)
    integer(kind=4), intent(in) :: n
    real(kind=ikind), intent(in) :: x, xx(n)
    integer(kind=4), intent(out) :: jlo
    integer(kind=4) :: inc, jhi, jm
    logical :: ascnd

    ! True if ascendinf order of table,
    ! False otherwise
    ascnd = (xx(n) > xx(1))

    ! Input guess is not useful
    ! Go immediately to bisection
    if((jlo <= 0) .OR. (jlo > n)) then
       jlo = 0
       jhi = n + 1
       goto 3
    end if

    ! Set the hunting increment
    inc = 1

    ! Hunt Up:
    if((x >= xx(jlo)) .EQV. ascnd) then
1      jhi = jlo + inc
       ! Done hunting, since off end of table
       if(jhi > n) then
          jhi = n + 1
          ! Not done hunting
       else if((x >= xx(jhi)) .EQV. ascnd) then
          jlo = jhi
          ! double the increment and try again
          inc = inc + inc
          goto 1
       end if
       ! Done hunting, value bracketed
       ! Hunt Down:
    else
       jhi = jlo
2      jlo = jhi - inc
       ! Done hunting, since off end of table
       if(jlo < 1) then
          jlo = 0
          ! Not done hunting
       else if((x < xx(jlo)) .EQV. ascnd) then
          jhi = jlo
          ! double the increment and try again
          inc = inc + inc
          goto 2
       end if
       ! Done hunting, value bracketed
    end if
    ! Hunt is done, so begin the final bisection phase

3   if(jhi-jlo == 1) return
    jm = (jhi + jlo)/2

    if((x > xx(jm)) .EQV. ascnd) then
       jlo = jm
    else
       jhi = jm
    end if
    goto 3
  end subroutine hunt

  !==================================================
  ! function smallerIndex: Array Shortener
  ! Will give a value to offset an array
  ! Refer to Numerical Recipes(pg. 113) for
  ! information about the function.
  ! @param:
  !   j - int, index of desired x value
  !   m - int, shorter array length
  !   n - int, full array length
  ! @return:
  !   k - int, offset for array
  !==================================================
  function smallerIndex(j, m, n) result(k)
    integer(kind=4) :: j, m, n
    integer(kind=4) :: k

    k = min(max(j-(m-1)/2, 1), n+1-m)
  end function smallerIndex

end module odeSolver
