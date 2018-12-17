program richardson
    !use n_bodies
    use orbit_functions
    use numerical_methods
    use dislin

    !use orbit_functions only: f_n_bodies => n_bodies
    implicit none
    
    logical, parameter :: DO_PLOT = .false., DO_EXAMPLE = .false.
    character(len=10), parameter :: FILENAMES(7) = & 
        ["earth.dat", "moon.dat", &
         "L1.dat", "L2.dat", "L3.dat", "L4.dat", "L5.dat"]
    character(len=*), parameter :: INIT_FILE = &
        "namelist_lagpoints_earthmoon.dat"
    integer, parameter :: M = 1
    real, parameter :: DT = 1d-1;
    
    integer :: i, l
    real :: tf, period, u(7*7)
    real, allocatable :: k_rchs(:), u0(:), k_rchs2(:)
    
    if (DO_EXAMPLE) then
        call example()
    else
        
        tf = 1d-5
    
        !call init_lagrange_statevector(INIT_FILE, u, period)
    
        !allocate(u0(3*7))
        !u0( 1: 3) = u(1:3)
        !u0( 4:12) = u(1*7+1:1*7+9)
        !u0(13:21) = u(4*7+1:4*7+9)
        
        !deallocate(u0)
        allocate(u0(4))
        u0 = [1.0d0, 0.0d0, 0.0d0, 1.0d0]
        
        call richardson_extrapolation(differential_operator = kepler, &
            temporal_scheme = runge_kutta, &
            t0 = 0d0, &
            tf = 1000d0, &
            dt= DT, &
            u0 = u0, &
            epsilon0 = 1d3, &
            k = k_rchs)
        call richardson_extrapolation(differential_operator = kepler, &
            temporal_scheme = runge_kutta, &
            t0 = 0d0, &
            tf = 1000d0, &
            dt= DT/2.0d0, &
            u0 = u0, &
            epsilon0 = 1d3, &
            k = k_rchs2)
    
        call qplot([((i*DT),i=1,size(k_rchs))], k_rchs, size(k_rchs))
        deallocate(k_rchs)
        
        call qplot([((i*DT/2.0d0),i=1,size(k_rchs2))], k_rchs2, size(k_rchs2))
        deallocate(k_rchs2)
        
    end if

    contains
    
    subroutine richardson_extrapolation(differential_operator, temporal_scheme, t0, tf, dt, u0, epsilon0, k)
        procedure(odes) :: differential_operator
        procedure(scheme) :: temporal_scheme
        real, intent(in) :: u0(:), t0, tf, dt, epsilon0 ! N(variables)
        real, allocatable, intent(out) :: k(:)
        integer :: m, n, i
        real, allocatable :: time_domain(:), s1(:), s2(:)
        
        m = int((tf - t0) / dt)
        n = size(u0)
        write(*,*) m, n
        !allocate(time_domain(0:m)!, time_domain_2(0:2*m))
        !allocate(sol_1(0:m, n), sol_2(0:2*m, n))
        allocate(s1(n), s2(n))
        allocate(k(m))
        
        !time_domain_1 = [(t0 + i * dt    , i = 0, m    )]
        !time_domain_2 = [(t0 + i * dt/2d0, i = 0, m * 2)]
        
        !sol_1(0,:) = u0
        !sol_2(0,:) = u0
        s1 = u0
        s2 = u0
        do i = 0, m - 1
            call temporal_scheme(differential_operator, t0+i*dt, t0+(i+1)*dt, s1, s1)
            !write(*,*) "Sol_1 en i=", i+1, ": ", s1, "\n"
            call temporal_scheme(differential_operator, t0+(2*i+0)*dt/2d0,t0+(2*i+1)*dt/2d0, s2, s2)
            !write(*,*) "Sol_2 en i=", i+1, ": ", s2, "\n"
            call temporal_scheme(differential_operator, t0+(2*i+1)*dt/2d0,t0+(2*i+2)*dt/2d0, s2, s2)
            !write(*,*) "Sol_2 en i=", i+2, ": ", s2, "\n"
            
            k(i+1) = norm2(s1 - s2) / (1d0 - 2d0 ** (-4d0))! * dt ** (-4d0)

            !write(*,*) k(i+1)
            !read(*,*)
        end do
            
    end subroutine
    
    function example_function(u, t) result(f)
        real, intent(in) :: u(:)
        real, intent(in) :: t
        real :: f(size(u))
        
        f = -u
    end
    
    subroutine example()
        integer, parameter :: M = 400
        real, parameter :: DT = 1d-1
        integer :: i
        real :: u_scheme(0:M,1), time(0:M), err(1:M)
        
        time = [(DT*i,i=0,M)]
        u_scheme(0,1) = 1d0
    
        do i = 0, M - 1
            call runge_kutta(example_function, i*dt, (i+1)*dt, u_scheme(i,:), u_scheme(i+1,:))
            err(i+1) = exp(-1d0*(i+1)*DT) - u_scheme(i+1,1)

            !write(*,*) k(i+1)
            !read(*,*)
        end do
        call qplot(time(1:), err, M)
    end subroutine example
end program