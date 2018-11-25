program $3bodies
    use functions
    use numerical_methods
    use cauchy
    use graphs
    implicit none
    
    logical, parameter :: DO_GRAPHS = .true.
    integer, parameter :: m = 400, s = 7
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
    double precision, parameter :: R_L2 =  473631.53d0, V_L2 = W_EM*R_L2
    double precision, parameter :: R_L3 = -495443.51d0, V_L3 = W_EM*R_L3
    double precision, parameter :: ANG_L4_DEG = atand(sqrt(3d0)*(R0_MOON+R0_EARTH)/(R0_MOON-R0_EARTH))
    double precision, parameter :: R_L4 = (R0_MOON+R0_EARTH)*sqrt(3d0)/2d0/sind(ANG_L4_DEG)
    double precision, parameter :: V_L4 = W_EM*R_L4
    double precision, parameter :: ANG_L5_DEG = -ANG_L4_DEG, R_L5 = R_L4, V_L5 = W_EM*R_L5
    
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
    double precision, dimension(3) :: x0s_l4 = R_L4 * [  cosd(ANG_L4_DEG), sind(ANG_L4_DEG), 0d0]
    double precision, dimension(3) :: v0s_l4 = V_L4 * [- sind(ANG_L4_DEG), cosd(ANG_L4_DEG), 0d0]
    double precision, dimension(3) :: x0s_l5 = R_L5 * [  cosd(ANG_L5_DEG), sind(ANG_L5_DEG), 0d0]
    double precision, dimension(3) :: v0s_l5 = V_L5 * [- sind(ANG_L5_DEG), cosd(ANG_L5_DEG), 0d0]
    
    double precision :: tf, time(0:m), u(0:m,1:n)
    integer :: i, j, k, io

    tf = PERIOD_MOON * 20
    
    time = [(tf/m*i,i=0,m)]
    mu(1:2) = [ MU_EARTH, MU_MOON]
    
    u(0,1:s) = mu
    u(0,1*s+1:4*s) = [x0e, x0m, x0s_l1, x0s_l2, x0s_l3, x0s_l4, x0s_l5]
    u(0,4*s+1:7*s) = [v0e, v0m, v0s_l1, v0s_l2, v0s_l3, v0s_l4, v0s_l5]
        
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
        
        call generate_poincare(x0s_l1, [(V_L1*(1+(i-5d0)/5d0*4d-1),i=0,10)], "poincare_L1.txt")  ! S/C L1
        call generate_poincare(x0s_l2, [(V_L2*(1+(i-5d0)/5d0*4d-1),i=0,10)], "poincare_L2.txt")  ! S/C L2
        call generate_poincare(x0s_l3, [(V_L3*(1+(i-5d0)/5d0*4d-1),i=0,10)], "poincare_L3.txt")  ! S/C L3
        call generate_poincare(x0s_l4, [(V_L4*(1+(i-5d0)/5d0*4d-1),i=0,10)], "poincare_L4.txt")  ! S/C L4
        call generate_poincare(x0s_l5, [(V_L5*(1+(i-5d0)/5d0*4d-1),i=0,10)], "poincare_L5.txt")  ! S/C L5
        
    end if

    contains
    
    subroutine generate_poincare(x0s, velocities, name)
        character(len=*), intent(in) :: name
        double precision, intent(in) :: velocities(:), x0s(3)
        
        integer, parameter :: mm = 1000
        double precision :: tf, time(0:mm), u(0:mm,1:21), v0s(3), param
        double precision, allocatable :: u_poincare(:,:), params(:)
        
        integer :: i, j, k, l
        
        l = size(velocities(:))
        tf = PERIOD_MOON * 100
        time = [(tf/mm*i,i=0,mm)]
        
        open(unit=13, file=name, status="REPLACE")
        close(13)
        open(unit=13, file=name, status="OLD")
        do i = 1, l
            param = velocities(i)
            v0s = [0d0, param, 0d0] 
            u(0,:) = [MU_EARTH, MU_MOON, 0d0, x0e, x0m, x0s, v0e, v0m, v0s]
            
            call cauchy_problem(time, n_bodies, runge_kutta, u)
            call plane_intersection_points(u(:,10:12), [0d0,1d0,0d0], u_poincare)
            k = size(u_poincare(:,1))
            do j = 1, k
                write(13,"(4F)") param, u_poincare(j,:)
            end do
            deallocate(u_poincare)
        end do
        close(13)
        
        open(unit=13, file=name, status="OLD")
        
        k = 0
        do
            read(13, *, iostat=io)
            if (io/=0) then
                exit
            end if
            k = k + 1
        end do
        
        allocate(params(k), u_poincare(k,4))
        rewind(13)
        
        do i = 1, k
            read(13,"(4F)") params(i), u_poincare(i,1:3)
        end do
        close(13)
        
        call poincare_map(params, u_poincare(:,1))
        deallocate(u_poincare)
    end subroutine generate_poincare
    
    
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
            if (r1*r2 < 0d0 .and. r2 >= 0d0) then
                n = n + 1
            end if
        end do
        
        allocate(inter_points(n, 3))
        
        j = 1
        do i = 1, m - 1
            r1 = dot_product(nor, u(i  ,:))
            r2 = dot_product(nor, u(i+1,:))
            if (r1*r2 < 0d0 .and. r2 > 0d0) then
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