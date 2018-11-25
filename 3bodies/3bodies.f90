program $3bodies
    use functions
    use numerical_methods
    use cauchy
    use graphs
    implicit none
    
    logical, parameter :: DO_GRAPHS = .true.
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
    double precision, parameter :: R_L2 =  473631.53d0, V_L2 = W_EM*R_L2
    double precision, parameter :: R_L3 = -495443.51d0, V_L3 = W_EM*R_L3
    double precision, parameter :: ANG_L4_DEG = atand(sqrt(3d0)*(R0_MOON+R0_EARTH)/(R0_MOON-R0_EARTH))
    double precision, parameter :: R_L4 = (R0_MOON+R0_EARTH)*sqrt(3d0)/2d0/sind(ANG_L4_DEG)
    double precision, parameter :: V_L4 = W_EM*R_L4
    double precision, parameter :: ANG_L5_DEG = -ANG_L4_DEG, R_L5 = R_L4, V_L5 = W_EM*R_L5
    
    integer, parameter :: m = 400, s = 7
    integer, parameter :: n = s * 7
    
    double precision, dimension(s) :: mu =  0d0
    double precision :: r0e_s = R0_EARTH, r0m_s = R0_MOON, r0s_s = 3d0
    double precision :: v0e_s = V0_EARTH, v0m_s = V0_MOON, v0s_s = 3d0

    double precision, dimension(3) :: x0e = [ -R0_EARTH, 0d0, 0d0]
    double precision, dimension(3) :: v0e = [ 0d0, -V0_EARTH, 0d0]
    double precision, dimension(3) :: x0m = [ R0_MOON, 0d0, 0d0]
    double precision, dimension(3) :: v0m = [ 0d0, V0_MOON ,0d0]
    double precision, dimension(3) :: x0s_l1 = [R_L1, 0d0, 0d0]
    double precision, dimension(3) :: v0s_l1 = [0d0, V_L1, 0d0]
    double precision, dimension(3) :: x0s_l2 = [R_L2, 0d0, 0d0]
    double precision, dimension(3) :: v0s_l2 = [0d0, V_L2, 0d0]
    double precision, dimension(3) :: x0s_l3 = [R_L3, 0d0, 0d0]
    double precision, dimension(3) :: v0s_l3 = [0d0, V_L3, 0d0]
    double precision, dimension(3) :: xx_l4 = R_L4 * [  cosd(ANG_L4_DEG), sind(ANG_L4_DEG), 0d0]
    double precision, dimension(3) :: vv_l4 = V_L4 * [- sind(ANG_L4_DEG), cosd(ANG_L4_DEG), 0d0]
    double precision, dimension(3) :: xx_l5 = R_L5 * [  cosd(ANG_L5_DEG), sind(ANG_L5_DEG), 0d0]
    double precision, dimension(3) :: vv_l5 = V_L5 * [- sind(ANG_L5_DEG), cosd(ANG_L5_DEG), 0d0]
    
    double precision :: tf, time(0:m), u(0:m,1:n), x_cg(0:m,3), v_cg(0:m,3)
    integer :: i, k

    tf = PERIOD_MOON * 10
    time = [(tf/m*i,i=0,m)]
    mu(1:2) = [ MU_EARTH, MU_MOON ]
    u(0,1:s) = mu
    u(0,1*s+1:4*s) = [x0e, x0m, x0s_l1, x0s_l2, x0s_l3, xx_l4, xx_l5]
    u(0,4*s+1:7*s) = [v0e, v0m, v0s_l1, v0s_l2, v0s_l3, vv_l4, vv_l5]
        
    call cauchy_problem(time, n_bodies, runge_kutta, u)
    
    
    open(unit=13, file="out.txt", status="REPLACE")
    do i = 0, m
        write(13,"(9F)") u(i, s+1:4*s)
    end do
    close(13)
    
    
    if (DO_GRAPHS) then
        call plot_orbit_xy(u, 1)  ! Earth
        call plot_orbit_xy(u, 2)  ! Moon
        call plot_orbit_xy(u / R_L5, 3)  ! S/C L1
        call plot_orbit_xy(u / R_L2, 4)  ! S/C L2
        call plot_orbit_xy(u / R_L3, 5)  ! S/C L3
        call plot_orbit_xy(u / R_L4, 6)  ! S/C L4
        call plot_orbit_xy(u / R_L5, 7)  ! S/C L5
    end if

    contains
    
    subroutine plane_intersection_points(u, nor, inter_points)
        double precision, intent(in) :: u(:,:), nor(3)  ! u(:,3)
        double precision, allocatable, intent(out) :: inter_points(:,:)
        
        integer :: m
        integer :: i, j, n
        double precision :: r1, r2, p1(3), p2(3), lambda, p_corte(3)
        
        m = size(u(:,1))
        n = 0
        do i = 1, m - 1
            r1 = dot_product(nor, u(i  ,:))
            r2 = dot_product(nor, u(i+1,:))
            if (r1*r2 < 0d0) then
                n = n + 1
            end if
        end do
        
        allocate(inter_points(n, 3))
        
        j = 1
        do i = 1, m - 1
            r1 = dot_product(nor, u(i  ,:))
            r2 = dot_product(nor, u(i+1,:))
            if (r1*r2 < 0d0) then
                ! Cambio de signo, se ha atravesado el plano.
                ! r:  p = p1 + lambda * (p2 - p1)
                ! pi: dot(p, normal) = 0
                ! p_corte pertenece a ambos.
                
                p1 = u(i  ,:)
                p2 = u(i+1,:)
                lambda = -dot_product(nor, p1) / dot_product(nor, p2-p1)
                p_corte = p1 + lambda * (p2 - p1)
                
                inter_points(j,:) = p_corte
                j = j + 1
            end if
        end do
    end subroutine plane_intersection_points

end program $3bodies