!  plots_with_dislin.f90 
!
!  FUNCTIONS:
!  plots_with_dislin - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: plots_with_dislin
!
!  PURPOSE:  Example of how to plot with dislin_mod library.
!
!****************************************************************************

    program plots_with_dislin
    
    use dislin_mod

    implicit none

    ! Variables
    integer, parameter :: M = 100
    real, parameter :: xmax = 30d0
    integer :: i
    real :: x(M), y(M), z(M)
    
    ! Assign Plot Data
    x = [(i * xmax / M, i = 1, M)]
    y = [(bessel_jn(0, x(i)), i = 1, M)]
    z = [(bessel_jn(1, x(i)), i = 1, M)]

    ! Normal plot:
    call plot ( x, y )
    
    ! Plot with custom color:
    call plot ( x, y, "RED" )
    
    ! Plot with labels:
    call plot ( x, y, "RED", xlabel = "x", ylabel = "J0(x)" )
    
    ! Multiple plots:
    call plot ( x, y, "RED", xlabel = "x", ylabel = "Jn(x)", hold_on = .true. )
    call plot ( x, z, "BLUE" )
    
    ! Plots with titles:
    call plot ( x, y, "RED", xlabel = "x", ylabel = "Jn(x)", hold_on = .true. )
    call plot ( x, z, "BLUE", hold_on = .true. )
    call plot_title( "1st kind Bessel Functions" )
    call plot_subtitle( "x**2 * y'' + x * y' + (x**2 - a**2) * y = 0" )
    call plot_title ( "(Cont)", 3 )
    call plot_end()
    
    ! Plots with legend
    call plot ( x, y, "RED", xlabel = "x", ylabel = "Jn(x)", hold_on = .true. )
    call plot ( x, z, "BLUE", hold_on = .true. )
    call plot_title( "1st kind Bessel Functions" )
    call plot_subtitle( "x**2 * y'' + x * y' + (x**2 - a**2) * y = 0" )
    call plot_title ( "(Cont)", 3 )
    call plot_legend( [ "J0(x)", "J1(x)" ] )
    call plot_end()
    
    ! Latex support is enabled by default.
    call plot ( x, y, "RED", xlabel = "x", ylabel = "$J_n(x)$", hold_on = .true. )
    call plot ( x, z, "BLUE", hold_on = .true. )
    call plot_title( "1st kind Bessel Functions" )
    call plot_title( "$x^2 \cdot \frac{d^2y}{dx^2} + x \cdot \frac{dy}{dx} + (x^2 - \alpha^2) \cdot y = 0$", 3 )
    call plot_legend( [ "$J_0(x)$", "$J_1(x)$" ] )
    call plot_end()
    
    ! Scattered points.
    call scatter ( x, y, "RED", xlabel = "x", ylabel = "$J_n(x)$", hold_on = .true. )
    call scatter ( x, z, "BLUE", hold_on = .true. )
    call plot_title( "1st kind Bessel Functions" )
    call plot_title( "$x^2 \cdot \frac{d^2y}{dx^2} + x \cdot \frac{dy}{dx} + (x^2 - \alpha^2) \cdot y = 0$", 3 )
    call plot_legend( [ "$J_0(x)$", "$J_1(x)$" ] )
    call plot_end()

    end program plots_with_dislin

