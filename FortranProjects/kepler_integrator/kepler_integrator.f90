program orbit_integrator
    use cauchy
    use orbit_functions
    use numerical_methods
    implicit none
    
    integer, parameter :: N = 100
    real, parameter :: PI = acos(-1d0)
    real :: tf = 2.0 * PI * 365d0 * 10d0
    real :: u(0:n, 0:3), time(0:n)
    integer :: i

    u(0,:) = [1d0, 0d0, 0d0, 1d0]
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

