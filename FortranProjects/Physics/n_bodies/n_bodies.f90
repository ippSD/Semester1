module n_bodies
    use orbit_functions!, only : f_n_bodies
    use numerical_methods!, only: runge_kutta
    use cauchy!, only: cauchy_problem
    use math
    implicit none
    real, parameter :: PI = acos(-1d0)
    contains
    
    subroutine calc_n_bodies(initial_conditions, time_domain, u)
        real, intent(in) :: time_domain(0:), initial_conditions(:)
        real, allocatable, intent(out) :: u(:,:)
        integer :: n, m
        real :: tf, period
        n = size(initial_conditions)
        if (mod(n,7) /= 0) then
            write(*,*) "Error: initial conditions are not multiple of 7."
            return
        end if
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