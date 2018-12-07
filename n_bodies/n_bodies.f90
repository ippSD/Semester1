module n_bodies
    use orbit_functions, only : f_n_bodies => n_bodies
    use numerical_methods
    use cauchy
    use auxiliary_functions
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
        call cauchy_problem(time_domain, f_n_bodies, runge_kutta, u)
    end subroutine calc_n_bodies
    
    subroutine calc_lagrange_points(init_file, time_domain, u)
        character(len=*), intent(in) :: init_file
        real, intent(in) :: time_domain(0:)
        real, allocatable, intent(out) :: u(:,:)
        integer :: n, m
        
        n = 7 * 7  ! 7 bodies (2 bodies + 5 lagranges) x 7 state vars.
        m = size(time_domain) - 1  ! size = m + 1, indexes in [0,m]
        allocate(u(0:m,n))
    
        call init_lagrange_statevector(init_file, u(0,:))
        call cauchy_problem(time_domain, f_n_bodies, runge_kutta, u)
    end subroutine calc_lagrange_points
    
    subroutine save_orbits(unit, u)
        integer, intent(in) :: unit
        real, intent(in) :: u(0:,:)
        
        integer :: i, m, n, s, idx_i, idx_f
        
        n = size(u(0,:))
        m = size(u(:,1)) - 1
        if(mod(n,7) /= 0) then
            write(*,*) "Error: State vector is not multiple of 7."
            return
        end if
        s = n / 7
        
        idx_i = s + 1
        idx_f = 4 * s
        do i = 0, m
            write(unit,*) u(i, idx_i:idx_f)
        end do
        close(unit)
    
    end subroutine
    
    subroutine save_body_orbit(unit, u, l)
        integer, intent(in) :: unit, l
        real, intent(in) :: u(0:,:)
        
        integer :: i, m, n, s, idx_i, idx_f
        
        n = size(u(0,:))
        m = size(u(:,1)) - 1
        if(mod(n,7) /= 0) then
            write(*,*) "Error: State vector is not multiple of 7."
            return
        end if
        s = n / 7
        if(l > s) then
            write(*,*) "Error: Body does not exist."
            return
        end if
        
        idx_i = s + (l-1) * 3 + 1
        idx_f = s + (l-1) * 3 + 3
        do i = 0, m
            write(unit,"(3F)") u(i, idx_i:idx_f)
        end do
        close(unit)
    end subroutine
    
    subroutine init_lagrange_statevector(init_file, u0, period)
        character(len=*) :: init_file
        real, intent(out) :: u0(7*7)
        real, optional, intent(out) :: period
        
        integer :: i
        real :: mu(7)
        real :: dcg, x_cg(3), omega_sys, p
        real :: x0_1(3), x0_2(3), x0_l(5, 3)
        real :: v0_1(3), v0_2(3), v0_l(5, 3)
        real :: mu_1, mu_2, sma, l_sma(5), l_ang(5)
        namelist /lagrpoints/ mu_1, mu_2, sma, l_sma, l_ang
    
        open(13, file=trim(init_file))
        read(13, nml=lagrpoints)
        close(13)
        
        x0_1 = - sma / (1d0+mu_1/mu_2) * [1d0, 0d0, 0d0]
        x0_2 =   sma / (1d0+mu_2/mu_1) * [1d0, 0d0, 0d0]
                
        omega_sys = sqrt(mu_1/sma**3d0)        
        p = 2*PI/omega_sys
                
        v0_1 = [0d0, 0d0, omega_sys] .vec. x0_1
        v0_2 = [0d0, 0d0, omega_sys] .vec. x0_2
                
        do i = 1, 5
            x0_l(i,:) = sma*l_sma(i)*[cosd(l_ang(i)),sind(l_ang(i)),0d0]+x0_1
            v0_l(i,:) = [0d0, 0d0, omega_sys] .vec. x0_l(i,:)
        end do
        mu = 0d0
        mu(1:2) = [mu_1, mu_2]
        
        u0 = [mu, x0_1, x0_2, transpose(x0_l), v0_1, v0_2, transpose(v0_l)]
        
        if(present(period)) then
            period = p
        end if

    end subroutine
    
end module n_bodies