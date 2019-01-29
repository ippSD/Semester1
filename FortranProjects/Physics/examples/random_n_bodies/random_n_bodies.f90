!---------------------------------------------------------------------------!
!   random_n_bodies                                                         !
!---------------------------------------------------------------------------!
!   Simulation of N-bodies problem on outter space with random              !
!   initial conditions.                                                     !
!---------------------------------------------------------------------------!
    
program random_n_bodies
    use n_bodies
    use orbit_plots
    use orbit_pointers
    use dislin_mod
    implicit none
    
    integer, parameter :: BODIES = 800
    integer, parameter :: M = 100
    integer, parameter :: N = BODIES * 7
    
    integer :: i
    real :: time(0:M)
    real, target :: u0(N)
    real, pointer :: r(:), v(:), x(:), y(:)
    real, allocatable :: u(:,:)
    
    time = [(i * 1000d0 / M, i = 0, M)]
    
    ! Get random initial values in [0, 1)
    call random_seed()
    call random_number(u0)
    
    ! Get negative values
    u0 = (u0 - 5d-1) * 2d0
    ! Positive mass with a classic factor of 1E6 km^3 s^-2
    u0(1:BODIES) = abs(u0(1:BODIES)) * 1d6
    ! Positions of order 1E4 km
    u0(BODIES+1:4*BODIES) = u0(BODIES+1:4*BODIES) * 1d4
    ! Velocities of order 1E0 km/s
    u0(4*BODIES+1:7*BODIES) = u0(4*BODIES+1:7*BODIES) * 5d-1
    
    nullify_z_coordinates: do i = 1, BODIES
        call pointer_to_body_position(u0, i, r)
        r(3) = 0d0
        nullify(r, v)
        call pointer_to_body_velocity(u0, i, v)
        v(3) = 0d0
        nullify(r, v)
    end do nullify_z_coordinates
    
    call calc_n_bodies(u0, time, u)
    
    plot_some_orbits: do i = 1, BODIES, BODIES / 10
        call pointer_to_body_x_vs_time(u, i, x)
        call pointer_to_body_y_vs_time(u, i, y)
        call plot(x, y, "RED")
        nullify(x, y)
    end do plot_some_orbits
        
end program random_n_bodies    