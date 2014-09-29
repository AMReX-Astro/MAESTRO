program compare_new_rates

  use rates_module
  use bl_types
  use bl_constants_module
  use network

  implicit none

  real(kind=dp_t), parameter :: temp = 1.1d9

  real(kind=dp_t) :: t9i23,t9i32, wwrate, rate


  ! build the common t factors
  call calc_tfactors(temp)
  t9i32 = t9i**(THREE/TWO)
  t9i23 = t9i**(TWO*THIRD)

  ! 15o(a,g)19ne
  call rate_he4_o15_to_ne19(rate)
  wwrate = t9i32*(exp(-5.849d0*t9i) + &
                  4.1d2*exp(-9.864d0*t9i) + &
                  2.3d3*exp(-11.837d0*t9i) + &
                  6.6d3*exp(-12.487d0*t9i))
  print *, '15o(a,g)19ne', rate/wwrate

  ! 19ne(p,g)20na
  call rate_p_ne19_to_na20(rate)
  wwrate = t9i32*(2.83d3*exp(-8.007d0*t9i) + &
                  3.16d3*exp(-8.936d0*t9i) + &
                  3.66d3*exp(-8.936d0*t9i))
  print *, '19ne(p,g)20na', rate/wwrate

  ! 14o(a,p)17f
  call rate_he4_o14_to_p_f17(rate)
  wwrate = 3.92d10*t9i23*(ONE + 0.0505d0*t9)**(-FIVE*SIXTH) &
       * exp(16.44d0-39.38d0*(ONE+0.0505d0*t9)**THIRD*t9i13)
  print *, '14o(a,p)17f', rate/wwrate

end program compare_new_rates
