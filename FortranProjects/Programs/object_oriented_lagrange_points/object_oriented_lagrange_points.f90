!---------------------------------------------------------------------------!
!   object_oriented_lagrange_points                                         !
!---------------------------------------------------------------------------!
!   Propagates the Earth-Moon Lagrange Points using Object Oriented         !
!   Programming on Fortran.                                                 !
!---------------------------------------------------------------------------!

program object_oriented_lagrange_points
    use cauchy_problem_solver, only: cauchy_problem
    use n_bodies, only: calc_n_bodies!, only: f_n_bodies, cauchy_problem, runge_kutta
    use orbit_plots!, only: plot_orbit_xy
    use orbit_exports
    use orbit_lagrange_points
    use orbit_objects!, only: GalaxyPointer
    implicit none
    
    character(len=10), parameter :: FILENAMES(7) = & 
        ["earth.dat", &
         "moon.dat",  &
         "L1.dat",    &
         "L2.dat",    &
         "L3.dat",    &
         "L4.dat",    &
         "L5.dat"]
    integer, parameter :: BODIES(5) = [1,2,3,4,5]
    integer, parameter :: M = 8000, P = size(BODIES)
    integer, parameter :: N = 7 * (P + 2)
    
    integer :: i, l
    real :: tf, period, time(0:M)
    real, target, allocatable :: u(:,:)
    real, allocatable :: u0(:)
    type(GalaxyPointer) :: mygalaxy(0:M)
    
    ! Define time domain based on Moon period
    period = T_MOON
    tf = period * 6d0/4d0 * 10d0
    time = [(tf/M*i,i=0,M)]
    
    ! Get Earth, Moon and some of its Lagrange Points as initial conditions
    call lagrange_points2initial_conditions(LP_EARTH_MOON, u0, BODIES)
    ! Solve N-bodies cauchy problem
    call calc_n_bodies(u0, time, u)
    
    ! For each state of the galaxy over time, each body
    ! on the galaxy points its state (position, speed, mass)
    ! to its equivalent on the state vector.
    do i = 0, M
        call mygalaxy(i)%set_pointer(u(i,:))
    end do
    
    ! Export Main Bodies' position
    do l = 1, 2
        open(unit=13, file=trim(FILENAMES(l)), status="REPLACE")
        call save_body_orbit(13, u, l)
        close(13)
    end do
    
    ! Export Lagrange Points' position
    do l = 1, p
        open(unit=13, file=trim(FILENAMES(BODIES(l)+2)), status="REPLACE")
        call save_body_orbit(13, u, l)
        close(13)
    end do
    
    call plots(u)

    contains
    
    subroutine plots(u)
        real, intent(in) :: u(0:,:)
        integer :: p
        
        ! Plot main bodies
        call plot_orbit_xy(u, 1, "EARTH")
        call plot_orbit_xy(u, 2, "MOON")
        
        ! Plot Lagrange Points
        p = size(u(0,:))/7 - 2
        do i = 1, p
            call plot_orbit_xy( u, i+2, FILENAMES(i+2)(1:2) )
        end do
    end subroutine plots
    
end program object_oriented_lagrange_points