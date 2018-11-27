program $3bodies
    use functions
    use numerical_methods
    use cauchy
    use objects
    implicit none
    
    integer, parameter :: m = 400, s = 2
    integer, parameter :: n = s * 7
    double precision, parameter :: PI = acos(-1d0)
    double precision, parameter :: PERIOD_MOON = 28d0 * 24d0 * 3600d0
    double precision, parameter :: MU_EARTH = 3.986d5
    double precision, parameter :: MU_MOON = 0.049d5
    double precision, parameter :: A_MOON = (MU_EARTH * (PERIOD_MOON / (2d0 * PI)) ** 2d0) ** (1d0/3d0)
    double precision, parameter :: V_MOON = 2d0*PI/PERIOD_MOON*A_MOON
    double precision, parameter :: R0_EARTH = (MU_MOON) / (MU_MOON + MU_EARTH) * A_MOON
    double precision, parameter :: V0_EARTH = (MU_MOON) / (MU_MOON + MU_EARTH) * V_MOON
    double precision, parameter :: R0_MOON = A_MOON - R0_EARTH
    double precision, parameter :: V0_MOON = V_MOON - V0_EARTH
    
    double precision, dimension(s) :: mu =  0d0

    double precision, dimension(3) :: x0e = [ -R0_EARTH, 0d0, 0d0]
    double precision, dimension(3) :: v0e = [ 0d0, -V0_EARTH, 0d0]
    double precision, dimension(3) :: x0m = [ R0_MOON, 0d0, 0d0]
    double precision, dimension(3) :: v0m = [ 0d0, V0_MOON ,0d0]
    
    double precision :: tf, time(0:m)
    double precision, target :: u(0:m,n)
    integer :: i
    type(GalaxyPointer) :: mygalaxy(0:m)

    tf = PERIOD_MOON * 20
    
    time = [(tf/m*i,i=0,m)]
    mu(1:2) = [ MU_EARTH, MU_MOON]
    
    u = 0d0
    u(0,:) = [mu, x0e, x0m, v0e, v0m]
    !mygalaxy(0) = Galaxy([mu, x0e, x0m, v0e, v0m])
    
    do i = 0, m
        call mygalaxy(i)%set_pointer(u(i,:))
    end do
    !mygalaxy(0).set_pointer(
    
    call cauchy_problem(time, n_bodies, runge_kutta, u)
    
    open(unit=13, file="out.txt", status="REPLACE")
    do i = 0, m
        write(13,"(6F)") mygalaxy(i)%bodies(1)%r(1), mygalaxy(i)%bodies(1)%r(2), mygalaxy(i)%bodies(1)%r(3), &
                         mygalaxy(i)%bodies(2)%r(1), mygalaxy(i)%bodies(2)%r(2), mygalaxy(i)%bodies(2)%r(3)
    end do
    close(13)

    contains
    
end program $3bodies