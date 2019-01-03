module lagrange_points
    use math
    implicit none
    
    type LagrangePoints
        real :: mu_1
        real :: mu_2
        real :: sma
        real :: r_li_sma(1:5), theta_li(1:5)
    end type LagrangePoints
    
    real, parameter :: MU_EARTH = 3.986d5, MU_MOON = 0.049d5
    real, parameter :: T_MOON = 28d0 * 24d0 * 36d2
    real, private, parameter :: PI = acos(-1d0)
    
    type(LagrangePoints), parameter :: LP_EARTH_MOON = LagrangePoints( &
        mu_1 = MU_EARTH, &
        mu_2 = MU_MOON, &
        sma = 3.844d5, &
        r_li_sma = [84.9100d-2, 116.7800d-2, 99.2916d-2, 1d0, 1d0], &
        theta_li = [0d0, 0d0, 180d0, 60d0, -60d0] &
    )
    
    contains
    
    subroutine pointer_to_body_position(state_vector, l, r_pointer)
        real, target, intent(in) :: state_vector(:)
        integer, intent(in) :: l
        real, pointer, intent(out) :: r_pointer(:)
        
        integer :: n
        
        n = size(state_vector) / 7
        
        r_pointer => state_vector(n+3*(l-1)+1:n+3*(l-1)+3)
    end subroutine
    
    subroutine pointer_to_body_velocity(state_vector, l, v_pointer)
        real, target, intent(in) :: state_vector(:)
        integer, intent(in) :: l
        real, pointer, intent(out) :: v_pointer(:)
        
        integer :: n
        
        n = size(state_vector) / 7
        
        v_pointer => state_vector(4*n+3*(l-1)+1:4*n+3*(l-1)+3)
    end subroutine
    
    subroutine lagrange_points2initial_conditions(lagrange_points, u0, opt_bodies)
        type(LagrangePoints), intent(in) :: lagrange_points
        real, allocatable, target, intent(inout) :: u0(:)
        integer, optional, intent(in) :: opt_bodies(:)
        
        integer, allocatable :: bodies(:)
        integer :: n, l, i
        real :: mu, freq, omega(3), xcg(3), r_i
        real, pointer :: r_i_vec(:), v_i_vec(:)
        
        mu = lagrange_points%mu_1 + lagrange_points%mu_2
        freq = sqrt(mu / lagrange_points%sma ** 3d0)
        omega = freq * [0d0, 0d0, 1d0]
        xcg = lagrange_points%sma * lagrange_points%mu_2 / mu * [1d0, 0d0, 0d0]
        
        if(allocated(u0)) deallocate(u0)
        if(present(opt_bodies)) then
            if(all(opt_bodies <= 5)) then
                bodies = opt_bodies
            else
                bodies = [1,2,3,4,5]
            end if
        else
            bodies = [1,2,3,4,5]
        end if
        
        l = size(bodies)
        n = l + 2
        
        allocate(u0(7*n))        
        
        ! Allocate masses
        u0(1:2) = [lagrange_points%mu_1, lagrange_points%mu_2]
        u0(3:n) = 0d0
        
        ! Allocate positions
        call pointer_to_body_position(u0, 1, r_i_vec)
        r_i_vec = 0d0 - xcg
        nullify(r_i_vec)
        call pointer_to_body_position(u0, 2, r_i_vec)
        r_i_vec = lagrange_points%sma * [1d0, 0d0, 0d0] - xcg
        nullify(r_i_vec)
        do i = 1, l
            call pointer_to_body_position(u0, i+2, r_i_vec)
            r_i = lagrange_points%r_li_sma(bodies(i)) * lagrange_points%sma
            r_i_vec = [ &
                cosd(lagrange_points%theta_li(bodies(i))), &
                sind(lagrange_points%theta_li(bodies(i))), &
                0d0 ]
            r_i_vec = r_i_vec * r_i
            r_i_vec = r_i_vec - xcg
            
            nullify(r_i_vec)
        end do
        
        !! Allocate velocities
        !  Main bodies' velocities
        u0(4*n+1  :4*n+3  ) = 0d0
        u0(4*n+1+3:4*n+3+3) = 0d0
        
        !  Lagrange Points' velocities
        u0(4*n+7:7*n) = 0d0
        
        !  Apply rotation of the fixed frame around inertial axes
        do i = 1, n
            call pointer_to_body_position(u0, i, r_i_vec)
            call pointer_to_body_velocity(u0, i, v_i_vec)
            
            v_i_vec = v_i_vec + omega .times. r_i_vec
            nullify(r_i_vec, v_i_vec)
        end do
        
    end subroutine
    
end module lagrange_points
