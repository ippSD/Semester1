program oscilator_integrator
    use cauchy
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
    real :: implicit_euler_matrix(2, 2)
    real, pointer :: x(:), p(:), x_plot(:), y_plot(:)
    integer :: i
    character(len=100) :: tit, subtit, subsubtit, filename
    character(len=20) :: legend_names(3)
    
    
    tit = "Integration of Oscilator equation."
    subtit = "$x'' + x = 0$"
    subsubtit = "$x_0 = 1 \quad ,\quad x'_0 = 0$"
    legend_names = ["Explicit Euler", "Runge Kutta", "Implicit Euler"]
    filename = "OscilatorVSNumericalIntegrators.png"
    
    x(0:) => u(0:n, 1)
    p(0:) => u(0:n, 2)
    
    if(DO_TIMEPLOT) then
        x_plot(0:) => time(0:n)
        y_plot(0:) => x(0:)
    else
        x_plot(0:) => x(0:)
        y_plot(0:) => p(0:)
    end if
    
    u(0,:) = [1d0, 0d0]
    tf = period * 2d0
    time = [(tf/n*i,i=0,n)]

    call cauchy_problem(time, f, euler_explicit, u)
    open(unit=13,file="oscilator_ee.txt")
    do i = 0, n
        write(13,*) x(i), p(i)
    end do
    close(13)
    
    call plot(x_plot, y_plot, plotcolor = "RED", hold_on = .true.)!, file_name = filename)
    
    call cauchy_problem(time, f, runge_kutta, u)
    open(unit=13,file="oscilator_rk4.txt")
    do i = 0, n
        write(13,*) x(i), p(i)
    end do
    close(13)
    
    call plot(x_plot, y_plot, plotcolor = "GREEN", hold_on = .true.)!, file_name = filename)

    
    !do i = 1, n
   !     dt = time(i) - time(i-1)
    !    implicit_euler_matrix(1,:) = [1d0, dt]
     !   implicit_euler_matrix(2,:) = [-dt, 1d0]
      !  implicit_euler_matrix = implicit_euler_matrix / (1d0 +dt**2d0)
        
    !    u(i,:) = matmul(implicit_euler_matrix, u(i-1,:))
    !end do
    call cauchy_problem(time, f, euler_implicit, u)
    open(unit=13,file="oscilator_ei.txt")
    do i = 0, n
        write(13,*) x(i), p(i)
    end do
    close(13)
    
    call plot(x_plot, y_plot, plotcolor = "BLUE", hold_on = .true.)
    
    call plot_legend(legend_names)
    call plot_title(tit)
    call plot_subtitle(subtit)
    call plot_title(subsubtit, 3)
    call plot_end()
    
end program oscilator_integrator

