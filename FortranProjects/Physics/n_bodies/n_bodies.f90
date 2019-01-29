module n_bodies
    use orbit_check, only: check_state_vector
    use orbit_functions, only: f_n_bodies
    use temporal_schemes, only: runge_kutta_4
    use cauchy_problem_solver, only: cauchy_problem
    implicit none
    contains
    
    subroutine calc_n_bodies(initial_conditions, time_domain, u)
        real, intent(in) :: time_domain(0:), initial_conditions(:)
        real, allocatable, intent(out) :: u(:,:)
        integer :: n, m
        real :: tf, period
        
        if ( .not. check_state_vector(initial_conditions) ) return
        n = size(initial_conditions)
        m = size(time_domain) - 1  ! size = m + 1, indexes in [0,m]
        allocate(u(0:m,n))
        
        u(0,:) = initial_conditions
        call cauchy_problem( &
            time_domain           = time_domain  , &
            differential_operator = f_n_bodies   , &
            temporal_scheme       = runge_kutta_4, &
            solution              = u              &
        )
    end subroutine calc_n_bodies
    
end module n_bodies