module rkf45
    use interfaces
    implicit none
    contains
    
subroutine rk4vec ( t0, m, u0, dt, f, u )

!*****************************************************************************80
!
!! RK4VEC takes one Runge-Kutta step for a vector ODE.
!
!  Discussion:
!
!    Thanks  to Dante Bolatti for correcting the final function call to:
!      call f ( t3, m, u3, f3 )
!    18 August 2016.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 August 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T0, the current time.
!
!    Input, integer ( kind = 4 ) M, the dimension of the system.
!
!    Input, real ( kind = 8 ) U0(M), the solution estimate at the current time.
!
!    Input, real ( kind = 8 ) DT, the time step.
!    Input, a function of the form 
!      function f ( u, t ) result(uprime)
!    which evaluates the derivative UPRIME(1:M) given the time T and
!    solution vector U(1:M).
!
!    Output, real ( kind = 8 ) U(M), the fourth-order Runge-Kutta solution 
!    estimate at time T0+DT.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) dt
  real ( kind = 8 ) f0(m)
  real ( kind = 8 ) f1(m)
  real ( kind = 8 ) f2(m)
  real ( kind = 8 ) f3(m)
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) u0(m)
  real ( kind = 8 ) u1(m)
  real ( kind = 8 ) u2(m)
  real ( kind = 8 ) u3(m)
  procedure(odes) :: f
!
!  Get four sample values of the derivative.
!
  f0 = f ( u0, t0 )

  t1 = t0 + dt / 2.0D+00
  u1(1:m) = u0(1:m) + dt * f0(1:m) / 2.0D+00
  f1 = f ( u1, t1 )

  t2 = t0 + dt / 2.0D+00
  u2(1:m) = u0(1:m) + dt * f1(1:m) / 2.0D+00
  f2 = f ( u2, t2 )

  t3 = t0 + dt
  u3(1:m) = u0(1:m) + dt * f2(1:m)
  f3 = f ( u3, t3 )
!
!  Combine them to estimate the solution U at time T1.
!
  u(1:m) = u0(1:m) + ( dt / 6.0D+00 ) * ( &
                 f0(1:m) &
     + 2.0D+00 * f1(1:m) &
     + 2.0D+00 * f2(1:m) &
     +           f3(1:m) )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

end module
