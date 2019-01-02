program $3bodies_oop
    use n_bodies
    use objects
    !use graphs
    implicit none
    
    logical, parameter :: DO_PLOT = .false.
    character(len=10), parameter :: FILENAMES(7) = & 
        ["earth.dat", "moon.dat", &
         "L1.dat", "L2.dat", "L3.dat", "L4.dat", "L5.dat"]
    character(len=*), parameter :: INIT_FILE = &
        "namelist_lagpoints_earthmoon.dat"
    integer, parameter :: M = 400
    
    integer :: i, l
    real :: tf, period, time(0:m)
    real, target, allocatable :: u(:,:)
    type(GalaxyPointer) :: mygalaxy(0:M)
    
    period = 28d0 * 24d0 * 36d2
    tf = period * 6d0/4d0
    time = [(tf/M*i,i=0,M)]
    
    call calc_lagrange_points(INIT_FILE, time, u)
    
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