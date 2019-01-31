!---------------------------------------------------------------------------!
!   richardson                                                              !
!---------------------------------------------------------------------------!
!   Use Richardson extrapolation on a sample function.                      !
!---------------------------------------------------------------------------!
    
program richardson_extrapolation_and_lyapunov_exp
    use richardson_extrapolation
    use cauchy_problem_solver
    use temporal_schemes
    use dislin_mod, only: plot, semilogy, plot_title, plot_legend, plot_end
    implicit none
    
    logical, parameter :: PLOT_TO_FILE = .true.
    integer, parameter :: M = 1000, M0 = 1, N = 1
    real, parameter :: TF = 1d2
    
    character( len = 50 ), parameter :: FILENAMES(5) = [  &
        "PlotRichardsonExtError.png"                    , &
        "PlotRichardsonAndCauchyAnalyticError.png"      , &
        "PlotRichardsonExtErrorVsRealError.png"         , &
        "PlotRichardsonExtAndRealLyapunovExp.png"       , &
        "PlotRichardsonExtVsRealLyapunovExpDeviation.png" &
    ]
    character( len = 50 ), parameter :: TITLES(5) = [  &
        "Richardson Error."                          , &
        "Committed real error:"                      , &
        "Expected Error:"                            , &
        "Evolution of Lyapunov Exponent"             , &
        "Lyapunov Exponent Deviation"                  &
    ]
    character( len = 50 ), parameter :: SUBTITLES(5) = [  &
        ""                                              , &
        "$Cauchy \quad VS \quad Richardson$"            , &
        "$Richardson \quad VS \quad Real$"              , &
        "$\epsilon_N = \epsilon_0^{\lambda_N\cdot t_N}$", &
        ""                                                &
    ]
    logical, parameter :: PLOT_SUBTITLES(5) = [ &
        .false., &
        .true. , &
        .true. , &
        .true. , &
        .false.  &
    ]
    character( len = 50 ), parameter :: LEGENDS(5, 2) = reshape( &
        [                                                       &
            "", ""                                            , &
            "Cauchy", "Richardson"                            , &
            "Richardson", "Real"                              , &
            "Real Lyapunov Exp.", "Richardson's Lyapunov Exp.", &
            "", ""                                              &
        ]     , &
        [5, 2], &
        order = [2, 1] &
    )
    logical, parameter :: PLOT_LEGENDS(5) = [ &
        .false., &
        .true. , &
        .true. , &
        .true. , &
        .false.  &
    ]
    character( len = 50 ), parameter :: LABELS(5, 2) = reshape( &
        [                                     &
            "Time [s]", "$Log10(ERR)$"      , &
            "Time [s]", "$Log10(ERR)$"      , &
            "Time [s]", "$Log10(ERR)$"      , &
            "Time [s]", "$\lambda_N$"       , &
            "Time [s]", "$\Delta \lambda_N$"  &
        ]     , &
        [5, 2], &
        order = [2, 1] &
    )
    
    
    integer :: i, l
    real :: u_cauchy(0:M, N), u_richard(0:M, N)
    real :: time(0:M)
    real, target :: richard_error(M0:M)
    real, target :: real_error_cauchy(M0:M), real_error_richard(M0:M)
    real, target :: richard_lyapunov(M0:M), real_lyapunov(M0:M)
    real, target :: deviation_lyapunov(M0:M)
    real, pointer :: y1(:), y2(:)
    
    ! Set time domain.
    time = [(i*TF/M, i = 0, M)]
    
    ! Set initial condition.
    u_cauchy(0,1)  = 1d0
    u_richard(0,1) = 1d0
    
    ! Solve Cauchy Problem
    call cauchy_problem(        &
        time                  , &
        cauchy_neg_exponential, &
        runge_kutta_4         , &
        u_cauchy                &
    )
        
    ! Solve Cauchy Problem with Richardson Extrapolation.
    call richardson_extrapolator (                      &
        time_domain           = time                  , &
        differential_operator = cauchy_neg_exponential, &
        temporal_scheme       = runge_kutta_4         , &
        order                 = 4                     , &
        solution              = u_richard             , &
        error                 = richard_error           &
    )
    
    ! Calculate Analytic Error
    real_error_cauchy  = abs( exp(-time(M0:M)) -  u_cauchy(M0:M, 1) )
    real_error_richard = abs( exp(-time(M0:M)) - u_richard(M0:M, 1) )
    
    ! Get Lyapunov Exp.
    richard_lyapunov = (richard_error/richard_error(1)) ** (log(10d0))
    real_lyapunov = (real_error_richard/real_error_richard(1)) ** (log(10d0))
    
    ! Get Lyapunov Exp. deviation Richardson Vs Real errors.
    deviation_lyapunov = richard_lyapunov / real_lyapunov
    plots_loop: do i = 1, 5
        if ( i == 1 ) then
            y1 => richard_error
            y2 => richard_error
        else if ( i == 2 ) then
            y1 => real_error_cauchy
            y2 => real_error_richard
        else if ( i == 3 ) then
            y1 => richard_error
            y2 => real_error_richard
        else if ( i == 4 ) then
            y1 => real_lyapunov
            y2 => richard_lyapunov
        else if ( i == 5 ) then
            y1 => deviation_lyapunov
            y2 => deviation_lyapunov
        end if
        save_plot: if ( PLOT_TO_FILE ) then
            call plot(                  &
                x = time(M0:M)        , &
                y = log10(y1)         , &
                color = "RED"         , &
                xlabel = LABELS(i, 1) , &
                ylabel = LABELS(i, 2) , &
                hold_on = .true.      , &
                filename = FILENAMES(i) &
            )
        else
            call plot(                 &
                x = time(M0:M)       , &
                y = log10(y1)        , &
                color = "RED"        , &
                xlabel = LABELS(i, 1), &
                ylabel = LABELS(i, 2), &
                hold_on = .true.       &
            )
        end if save_plot
        
        multiplot: if ( PLOT_LEGENDS(i) ) then
            call plot(           &
                x = time(M0:M) , &
                y = log10(y2)  , &
                color = "BLUE" , &
                hold_on = .true. &
            )
            call plot_legend(LEGENDS(i, 1:2))
        end if multiplot
        
        call plot_title(TITLES(i))
        if( PLOT_SUBTITLES(i) ) call plot_title( SUBTITLES(i), 3 )
        call plot_end()
    end do plots_loop
    
    contains
    
    !-----------------------------------------------------------------------!
    !   ( FUNCTION ) cauchy_neg_exponential                                 !
    !-----------------------------------------------------------------------!
    !   Derivative of the negative exponential for Cauchy Problem.          !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(1)) u                                           !
    !               ! State vector of the cauchy problem at time t.         !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t                                              !
    !               ! Time at which the state vector are given.             !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(1)) f                                           !
    !               ! Derivative of the state vector at time t.             !
    !---------------!-------------------------------------------------------!
    function cauchy_neg_exponential(u,t) result(f)
        real, intent(in) :: u(:), t
        real :: f(size(u))
        
        f(1) = - u(1)
    end function
    
    !-----------------------------------------------------------------------!
    !   ( FUNCTION ) neg_exponential                                        !
    !-----------------------------------------------------------------------!
    !   Derivative of the negative exponential for Cauchy Problem.          !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t                                              !
    !               ! Time argument.                                        !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(1)) f                                           !
    !               ! Negative exponential of time argument (exp(-t)).      !
    !---------------!-------------------------------------------------------!
    function neg_exponential(t) result(f)
        real, intent(in) :: t
        real :: f(1)
        
        f(1) = exp(-t)
    end function
    
end program richardson_extrapolation_and_lyapunov_exp