program poincare
    use n_bodies
!    use graphs
    implicit none
    
    logical, parameter :: DO_PLOT = .true.
    character(len=*), parameter :: INIT_FILE = &
        "namelist_lagpoints_earthmoon.dat"
    integer, parameter :: M = 400, S = 7
    integer, parameter :: N = 7 * S
    real, allocatable :: u0(:,:)

    call calc_lagrange_points(INIT_FILE, [0d0], u0)    
    
    if (DO_PLOT) then
        call plots(u0(0,:))
    end if

    contains
    
    subroutine plots(u0)
        real, intent(in) :: u0(:)
        integer, parameter :: NN = 501
        real, parameter :: EPS_X_MAX = 1d4, EPS_VY_MAX = 2d0
        integer :: i, l
        real :: tf, period, time(0:m), eps_x(NN), eps_vy(NN)
        
        period = 28d0 * 24d0 * 36d2
        tf = period * 6d0/4d0
        time = [(tf/M*i,i=0,M)]
        eps_x  = [( EPS_X_MAX*(i-(NN-1d0))/(NN-1d0),i=0,NN-1)]
        eps_vy = [(EPS_VY_MAX*(i-(NN-1d0))/(NN-1d0),i=0,NN-1)]
        
        call plot_poincare(time, eps_x,  u0,   7+(3-1)*3+1, 3)  ! Poincare on L1, X0_L1
        call plot_poincare(time, eps_vy, u0, 4*7+(3-1)*3+2, 3)  ! Poincare on L1, VY0_L1
        call plot_poincare(time, eps_x,  u0,   7+(6-1)*3+1, 6)  ! Poincare on L4, X0_L4
        call plot_poincare(time, eps_vy, u0, 4*7+(6-1)*3+2, 6)  ! Poincare on L4, VY0_L4 
    end subroutine
    
    subroutine plot_poincare(time_domain, eps, design_initial_conditions, var_idx, body_idx)
        integer, intent(in) :: var_idx, body_idx
        real, intent(in) :: time_domain(0:)
        real, intent(in) :: eps(:), design_initial_conditions(:)
        
        integer :: i, m, n, s, l, q_i, q_last
        real, allocatable :: u0(:), x_poincare(:,:), temp(:,:), temp2(:,:)
        real, allocatable, target :: u(:,:)
        real, pointer :: x_body(:,:)
        
        m = size(time_domain) - 1
        n = size(design_initial_conditions)
        s = n / 7
        l = size(eps)
        
        allocate(u0(n))
        
        do i = 1, l
            u0 = design_initial_conditions
            u0(var_idx) = u0(var_idx) + eps(i)
            call calc_n_bodies(u0, time_domain, u)
            
            x_body => u(:,s+(body_idx-1)*3+1:s+(body_idx-1)*3+3)
            call plane_intersection_points(x_body, [0d0,1d0,0d0], temp)
            
            if(allocated(x_poincare)) then
                q_last = size(x_poincare(:,1))
                q_i = size(temp(:,1))
                allocate(temp2(q_last+q_i,2))
            
                temp2(1:q_last,:) = x_poincare(:,:)
                temp2(q_last+1:,1) = u0(var_idx)
                temp2(q_last+1:,2) = temp(:,1)
                deallocate(x_poincare)
            else
                q_last = 0
                q_i = size(temp(:,1))
                allocate(temp2(q_last+q_i,2))
                temp2(q_last+1:,1) = u0(var_idx)
                temp2(q_last+1:,2) = temp(:,1)
            end if
            
            allocate(x_poincare(q_last+q_i,2))
            x_poincare = temp2
            nullify(x_body)
            deallocate(temp, temp2, u)
        end do
        
        call qplsca(x_poincare(:,1), x_poincare(:,2), q_last+q_i)
    end subroutine
    
    subroutine plane_intersection_points(u, nor, inter_points)
        real, intent(in) :: u(:,:), nor(3)  ! u(:,3)
        real, allocatable, intent(out) :: inter_points(:,:)
        
        integer :: m, i, j, n
        real :: r1, r2, p1(3), p2(3), lambda, p_corte(3)
        
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
    
    
end program poincare