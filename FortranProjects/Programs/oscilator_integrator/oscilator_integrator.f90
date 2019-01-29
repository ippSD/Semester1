!---------------------------------------------------------------------------!
!   oscilator_integration                                                   !
!---------------------------------------------------------------------------!
!   Integrates the harmonic oscilator with different temporal schemes.      !
!---------------------------------------------------------------------------!
    
program oscilator_integrator
    use cauchy
    use ode_interfaces, only: propagator => temporal_scheme
    use hamiltonian_functions, only: f => oscilator
    use numerical_methods
    use dislin_mod
    implicit none
    
    logical, parameter :: DO_TIMEPLOT = .true.
    
    integer, parameter :: N = 100
    real, parameter :: PI = acos(-1d0)
    real :: period = 2.0 * PI, tf, dt
    real, target :: u(0:n, 2)
    real, target :: time(0:n)
    real, pointer :: x(:), p(:), x_plot(:), y_plot(:)
    integer :: i, j, k
    character(len=100) :: tit, subtit, subsubtit
    character(len=20) :: datafiles(3), legends(3), colors(3)
    character(len=20) :: plotfiles(2), xlabels(2), ylabels(2)
    
    procedure(propagator), pointer :: selected_propagator    
    
    tit       = "Integration of Oscilator equation."
    subtit    = "$x'' + x = 0$"
    subsubtit = "$x_0 = 1 \quad ,\quad x'_0 = 0$"
    
    legends   = [  "Explicit Euler",     "Runge Kutta 4",   "Implicit Euler"]
    colors    = [             "RED",             "GREEN",             "BLUE"]
    datafiles = ["oscilator_ee.dat", "oscilator_rk4.dat", "oscilator_ei.dat"]
    
    plotfiles = ["x_vs_time.png", "p_vs_x.png"]
    xlabels = ["$Time\ [s]$", "$X\ [m]$"  ]
    ylabels = ["$X\ [m]$"   , "$P\ [m/s]$"]
    
    x(0:) => u(0:n, 1)
    p(0:) => u(0:n, 2)
    
    u(0,:) = [1d0, 0d0]
    tf = period * 2d0
    time = [(tf/n*i,i=0,n)]
    
    do j = 1, 2
        if( j == 1 ) then
            x_plot(0:) => time(0:n)
            y_plot(0:) => x(0:)
        else
            x_plot(0:) => x(0:)
            y_plot(0:) => p(0:)
        end if
    
        do i = 1, 3
            if ( i == 1 ) selected_propagator => euler_explicit
            if ( i == 2 ) selected_propagator => runge_kutta_4
            if ( i == 3 ) selected_propagator => euler_implicit

            call cauchy_problem( time, f, selected_propagator, u )
            open( unit = 13, file = datafiles(i) )
            do k = 0, n
                write(13,*) time(k), x(k), p(k)
            end do
            close(13)
    
            call plot(                     &
                x_plot                   , &
                y_plot                   , &
                colors(i)                , &
                xlabel = xlabels(j)      , &
                ylabel = ylabels(j)      , &
                hold_on = .true.         , &
                filename = plotfiles(j)    &
            )
        end do
        
        call plot_legend(legends)
        call plot_title(tit, 1)
        call plot_subtitle(subtit)
        call plot_title(subsubtit, 3)
        call plot_end()
    end do
    
end program oscilator_integrator

