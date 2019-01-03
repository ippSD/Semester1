program $3bodies_oop
    use n_bodies, only: f_n_bodies, cauchy_problem, runge_kutta
    use graphs_n_bodies, only: plot_orbit_xy
    use lagrange_points
    use objects, only: GalaxyPointer
    implicit none
    
    logical, parameter :: DO_PLOT = .false.
    character(len=10), parameter :: FILENAMES(7) = & 
        ["earth.dat", &
         "moon.dat",  &
         "L1.dat",    &
         "L2.dat",    &
         "L3.dat",    &
         "L4.dat",    &
         "L5.dat"]
    integer, parameter :: M = 400, N = 7*7
    
    integer :: i, l
    real :: tf, period, time(0:M)
    real, target, allocatable :: u(:,:)
    real, allocatable :: u0(:)
    type(GalaxyPointer) :: mygalaxy(0:M)
    
    period = 28d0 * 24d0 * 36d2
    tf = period * 6d0/4d0
    time = [(tf/M*i,i=0,M)]
    
    allocate(u(0:M,1:N))
    call lagrange_points2initial_conditions(LP_EARTH_MOON, u0, [1,2,3,4,5])
    u(0,:) = u0
    deallocate(u0)
    
    call cauchy_problem(time, f_n_bodies, runge_kutta, u)
    
    ! For each state of the galaxy over time, each body
    ! on the galaxy points its state (position, speed, mass)
    ! to its equivalent on the state vector.
    do i = 0, M
        call mygalaxy(i)%set_pointer(u(i,:))
    end do
    
    do l = 1, 7
        open(unit=13, file=trim(FILENAMES(l)), status="REPLACE")
        call save_body_orbit(13, u, l)
    end do
    
    if (DO_PLOT) then
        call plots(u)
    end if

    contains
    
    subroutine plots(u)
        real, intent(in) :: u(0:,:)
        
        call plot_orbit_xy(u, 1)
        call plot_orbit_xy(u, 2)
        call plot_orbit_xy(u, 3)
        call plot_orbit_xy(u, 4)
        call plot_orbit_xy(u, 5)
        call plot_orbit_xy(u, 6)
        call plot_orbit_xy(u, 7)
    end subroutine plots
    
end program $3bodies_oop