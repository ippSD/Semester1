program richardson
    use dislin
    use functions
    use numerical_methods
    use cauchy
    implicit none
    
    logical, parameter :: DO_GRAPHS = .true.
    integer, parameter :: m = 400, s = 4
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
    double precision, parameter :: W_EM = V0_EARTH / R0_EARTH
    
    double precision, parameter :: R_L1 =  326400.67d0, V_L1 = W_EM*R_L1
    double precision, parameter :: ANG_L4_DEG = atand(sqrt(3d0)*(R0_MOON+R0_EARTH)/(R0_MOON-R0_EARTH))
    double precision, parameter :: R_L4 = (R0_MOON+R0_EARTH)*sqrt(3d0)/2d0/sind(ANG_L4_DEG)
    double precision, parameter :: V_L4 = W_EM*R_L4
    
    double precision, dimension(s) :: mu =  0d0
    double precision :: r0e_s = R0_EARTH, r0m_s = R0_MOON, r0s_s = 3d0
    double precision :: v0e_s = V0_EARTH, v0m_s = V0_MOON, v0s_s = 3d0

    double precision, dimension(3) :: x0e = [ -R0_EARTH, 0d0, 0d0]
    double precision, dimension(3) :: v0e = [ 0d0, -V0_EARTH, 0d0]
    double precision, dimension(3) :: x0m = [ R0_MOON, 0d0, 0d0]
    double precision, dimension(3) :: v0m = [ 0d0, V0_MOON ,0d0]
    double precision, dimension(3) :: x0s_l1 = [R_L1, 0d0, 0d0]
    double precision, dimension(3) :: v0s_l1 = [0d0, V_L1, 0d0]
    double precision, dimension(3) :: x0s_l4 = R_L4 * [  cosd(ANG_L4_DEG), sind(ANG_L4_DEG), 0d0]
    double precision, dimension(3) :: v0s_l4 = V_L4 * [- sind(ANG_L4_DEG), cosd(ANG_L4_DEG), 0d0]
    
    double precision :: tf, time1(0:m), time2(0:m*2), u1(0:m,1:n), u2(0:m*2,1:n), er(1:m,1:n)
    integer :: i, j, k, io

    tf = PERIOD_MOON * 20
    
    time1 = [(tf/m*i,i=0,m)]
    time2 = [(tf/m/2d0*i,i=0,m*2)]
    
    mu(1:2) = [ MU_EARTH, MU_MOON]
    
    u1(0,1:s) = mu
    u1(0,1*s+1:4*s) = [x0e, x0m, x0s_l1, x0s_l4]
    u1(0,4*s+1:7*s) = [v0e, v0m, v0s_l1, v0s_l4]
    
    u2(0,1:s) = mu
    u2(0,1*s+1:4*s) = [x0e, x0m, x0s_l1, x0s_l4]
    u2(0,4*s+1:7*s) = [v0e, v0m, v0s_l1, v0s_l4]
    
    call cauchy_problem(time1, n_bodies, runge_kutta, u1)
    
    call cauchy_problem(time2, n_bodies, runge_kutta, u2)
    
    er = (u1(1:m,:)-u2(1:2*m:2,:))/(1d0-0.5d0**4d0)
    
    call qplot(time1,er,m)
    
end program