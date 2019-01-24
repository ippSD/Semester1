!---------------------------------------------------------------------------!
!   poincare_map                                                            !
!---------------------------------------------------------------------------!
!   Draws a Celestial body's Poincaré map by disturbing its                 !
!   initial condition.                                                      !
!---------------------------------------------------------------------------!
    
program poincare_map
    use orbit_lagrange_points
    use orbit_pointers
    use n_bodies
    use dislin_mod
    implicit none
    
    logical, parameter :: DO_PLOT = .true.
    integer, parameter :: M = 400, S = 7
    integer, parameter :: N = 7 * S
    real, allocatable :: u0(:)
        
    call lagrange_points2initial_conditions(LP_EARTH_MOON, u0)
    if (DO_PLOT) then
        call plots(u0)
    end if

    contains
    
    subroutine plots(u0)
        real, intent(in) :: u0(:)
        integer, parameter :: NN = 501
        real, parameter :: EPS_X_MAX = 1d4, EPS_VY_MAX = 2d0
        integer :: i, l
        real :: tf, period, time(0:m), eps_x(NN), eps_vy(NN)
        
        integer :: plot_bodies(2) = [1, 4] + 2
        
        period = T_MOON
        tf = period * 6d0/4d0
        time = [(i*tf/M,i=0,M)]
        
        ! Set disturvance on initial values X and VY
        eps_x  = [( EPS_X_MAX*(i-(NN-1d0))/(NN-1d0),i=0,NN-1)]
        eps_vy = [(EPS_VY_MAX*(i-(NN-1d0))/(NN-1d0),i=0,NN-1)]
        
        call plot_poincare(                 &
            time                          , &
            eps_x                         , &
            u0                            , &
            7+(3-1)*3+1                   , &
            3                             , &
            "$X\ [km]$"                   , &
            "L1 Lag. Point Poincare over X" &
        )
            
        
        call plot_poincare(time, eps_vy, u0, 4*7+(3-1)*3+2, 3)  ! Poincare on L1, VY0_L1
        call plot_poincare(time, eps_x,  u0,   7+(6-1)*3+1, 6)  ! Poincare on L4, X0_L4
        call plot_poincare(time, eps_vy, u0, 4*7+(6-1)*3+2, 6)  ! Poincare on L4, VY0_L4 
    end subroutine
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) plot_poincare                                        !
    !-----------------------------------------------------------------------!
    !   Plots the Poincaré Map.                                             !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(0:M)) time_domain                               !
    !               ! Time domain in which the two-bodies problem evaluates.!
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(NN)) eps                                        !
    !               ! Disturbance vector. The disturbed initial condition   !
    !               ! is added the values: CI_i = CI_0 + EPS_i.             !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) design_initial_conditions                   !
    !               ! Undisturbed initial conditions.                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (integer) var_idx                                     !
    !               ! Index of the state vector whose initial condition     !
    !               ! will be disturbed.                                    !
    !---------------!-------------------------------------------------------!
    !   IN          ! (integer) body_idx                                    !
    !               ! Body whose Poincare Map will be draw.                 !
    !---------------!-------------------------------------------------------!
    subroutine plot_poincare(          &
            time_domain              , &
            eps                      , &
            design_initial_conditions, &
            var_idx                  , &
            body_idx                 , &
            title                    , &
            xlabel                     &
        )
        integer, intent(in) :: var_idx, body_idx
        real, intent(in) :: time_domain(0:)
        real, intent(in) :: eps(:), design_initial_conditions(:)
        character( len = * ), intent(in) :: xlabel, title
        
        integer :: i, m, n, s, l, q_i, q_last
        real, allocatable :: u0(:), x_poincare(:,:), temp(:,:), temp2(:,:)
        real, allocatable, target :: u(:,:)
        real, pointer :: r_body(:,:)
        
        m = size(time_domain) - 1
        n = size(design_initial_conditions)
        s = n / 7
        l = size(eps)
        
        allocate(u0(n))
        
        do i = 1, l
            u0 = design_initial_conditions
            u0(var_idx) = u0(var_idx) + eps(i)
            call calc_n_bodies(u0, time_domain, u)
            
            call pointer_to_body_r_vs_time( u, body_idx, r_body )
            ! x_body => u(:,s+(body_idx-1)*3+1:s+(body_idx-1)*3+3)
            call plane_intersection_points(r_body, [0d0,1d0,0d0], temp)
            
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
            nullify(r_body)
            deallocate(temp, temp2, u)
        end do
        
        call scatter (           &
                x_poincare(:,1), &
                x_poincare(:,2), &
                "RED"          , &
                xlabel = xlabel, &
                hold_on = .true. &
        )
        call plot_title( title )
        call plot_end()
    end subroutine
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) plane_intersection_points                            !
    !-----------------------------------------------------------------------!
    !   Computes the intersection points of a trayectory with a plane.      !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(M,N)) r                                         !
    !               ! Trayectory matrix. Rows contain the XYZ points of the !
    !               ! trayectory over the time and columns X,Y,Z components.!
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(3)) nor                                         !
    !               ! Normal vector of the plane.                           !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(?,3)) inter_points                              !
    !               ! Intersection points' matrix. Rows contain the         !
    !               ! intersection points, unknown at first.                !
    !---------------!-------------------------------------------------------!
    subroutine plane_intersection_points(r, nor, inter_points)
        real, intent(in) :: r(:,:), nor(3)
        real, allocatable, intent(out) :: inter_points(:,:)
        
        integer :: m, i, j, n
        real :: r1, r2, p1(3), p2(3), lambda, p_corte(3)
        
        m = size(r(:,1))
        n = 0
        do i = 1, m - 1
            r1 = dot_product(nor, r(i  ,:))
            r2 = dot_product(nor, r(i+1,:))
            if (r1*r2 < 0d0 .and. r2 >= 0d0) then
                n = n + 1
            end if
        end do
        
        allocate(inter_points(n, 3))
        
        j = 1
        do i = 1, m - 1
            r1 = dot_product(nor, r(i  ,:))
            r2 = dot_product(nor, r(i+1,:))
            if (r1*r2 < 0d0 .and. r2 > 0d0) then
                ! Cambio de signo, se ha atravesado el plano.
                ! r:  p = p1 + lambda * (p2 - p1)
                ! pi: dot(p, normal) = 0
                ! p_corte pertenece a ambos.
                
                p1 = r(i  ,:)
                p2 = r(i+1,:)
                lambda = -dot_product(nor, p1) / dot_product(nor, p2-p1)
                p_corte = p1 + lambda * (p2 - p1)
                
                inter_points(j,:) = p_corte
                j = j + 1
            end if
        end do
    end subroutine plane_intersection_points
    
    
end program poincare_map