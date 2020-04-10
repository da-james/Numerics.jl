module module_ode
  implicit none

contains

  subroutine fdm(f,v,t,y0,n,h,x,y)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: h
    real, dimension(2), intent(in) :: v, t ,y0
    real, dimension(n), intent(out) :: x, y
    external :: f
    !======================================================================
    !
    ! function fdm: Finite Difference Method
    ! This will solve out a 2nd order DE using the
    ! Finite Difference approach
    !
    !======================================================================
    integer, parameter :: m = n - 2
    real :: s1, s2, s3
    real, dimension(m) :: L, D, U, f_b
    integer :: i


    s1 = 1/h**2 - v(1)/(2*h)
    s2 = v(2) - 2/h**2
    s3 = 1/h**2 + v(1)/(2*h)

    L(:) = s1
    D(:) = s2
    U(:) = s3

    do i=1,n-2
       call f(x(i),f_b(i))
    enddo
    f_b(1) = f_b(1) - s1*y0(1)
    f_b(n) = f_b(n) - s3*y0(2)

    call tridiag(n,L,D,U,f_b,y)
    y(2:) = f_b

  end subroutine fdm

  subroutine tridiag(jmx,a,b,c,f,q)
    integer, intent(in) :: jmx
    real, dimension(jmx), intent(inout) :: a, b, c, f, q
    !======================================================================
    !
    ! subroutine tridiag: Thomas Algorithm
    ! This subroutine will solve out a tridiagonal matrix system of equations
    ! using the Thomas Algorithm. Where the system is T*x = b, where T is a
    ! matrix made up of U, D, and L, the upper diagonal, diagonal, and lower
    ! diagonal, respectively.
    ! @param:
    ! n - int, intent(in); size of the diagonal of the matrix
    ! L - real, intent(inout), dimension(n-1): lower diagonal of the matrix
    ! D - real, intent(inout), dimension(n): diagonal of the matrix
    ! U - real, inout, dimension(n-1): upper diagonal of the matrix
    ! x - real, inout, dim(n): variables for system
    ! b - real, inout, dim(n): constant terms
    !
    !======================================================================
    real :: p
    integer :: j

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
