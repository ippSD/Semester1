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
    
    implicit none
    
    integer, parameter :: n = 10000
    real (kind = 8) :: tf = 2.0*acos(-1.0)*365.0*10
    real (kind = 8) :: u(0:n, 0:3)
    integer :: i

    u(0,:) = [1.0, 0.0, 0.0, 1.0]

    write(*,*) u(0,0:1)

    call cauchy_problem(u, tf, kepler)
    open(unit=13,file="out.txt")
    do i = 0, n
        write(13,*) u(i, 0:1)
    end do
    write(*,*) u(0, 0:1)
    write(*,*) u(n, 0:1)
    close(13)

    end program orbit_integrator

