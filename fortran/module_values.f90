module module_values
  implicit none

  ! real precision values
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(30, 291)

  ! integer precision values
  ! 8  bit -128 to 127d2
  integer, parameter :: i8  = selected_int_kind(2)
  ! 16 bit -32768 to 32767d4
  integer, parameter :: i16 = selected_int_kind(4)
  ! 32 bit -2147483648 to 2147483647d9
  integer, parameter :: i32 = selected_int_kind(9)
  ! 64 bit -9223372036854775808 to 9223372036854775807
  integer, parameter :: i64 = selected_int_kind(15)

  ! constants used in math
  real(dp), parameter :: PI = 4.0*atan(1.0)      ! pi calculated
  real(dp), parameter :: GRAVG = 6.67430d-11     ! gravitational constant [N*m^2/kg^2]


end module module_values
