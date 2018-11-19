!  orbit_integrator.f90 
!
!  FUNCTIONS:
!  orbit_integrator - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: orbit_integrator
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program orbit_integrator

    use cauchy
    use functions
    use numerical_methods
    
    implicit none
    
    integer, parameter :: n = 100
    real (kind = 8) :: tf = 2.0*acos(-1.0)*365.0*10
    real (kind = 8) :: u(0:n, 0:3), time(0:n)
    integer :: i

    u(0,:) = [1.0, 0.0, 0.0, 1.0]
    time = [(tf/n*i,i=0,n)]

    call cauchy_problem(time, kepler, euler_explicit, u)
    open(unit=13,file="out_ee.txt")
    do i = 0, n
        write(13,*) u(i, 0:1)
    end do
    close(13)

    call cauchy_problem(time, kepler, runge_kutta, u)
    open(unit=13,file="out_rk4.txt")
    do i = 0, n
        write(13,*) u(i, 0:1)
    end do
    close(13)

end program orbit_integrator

