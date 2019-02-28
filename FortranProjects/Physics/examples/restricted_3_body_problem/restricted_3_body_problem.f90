program restricted_3_body_problem
    use orbit_functions, only: cr3bp, cr3bp_u, cr3bp_uu
    use cauchy_problem_solver, only: cp => cauchy_problem
    use Cauchy_Problem, only: System_matrix
    use temporal_schemes_my, only: runge_kutta_4
    use Non_Linear_Systems
    use Numerical_Recipes
    use dislin_mod
    implicit none
    
    integer, parameter :: N = 6, M = 10000
    real, parameter :: TF = 2d0 * acos(-1d0) * 300, MU = (0.049d5)/( 0.049d5+3.986d5)
    real :: u(0:M,N), time(0:M), u0(5,N), a(N,N)
    complex :: lambda(N)
    integer :: i
    
    u0(1,:) = 1d0*[5d-1, 0d0, 0d0, 0d0, 0d0, 0d0]
    u0(2,:) = 1d0*[1d0, 0d0, 0d0, 0d0, 0d0, 0d0]
    u0(3,:) = 1d0*[-1d0, 0d0, 0d0, 0d0, 0d0, 0d0]
    u0(4,:) = 1d0*[cosd(60d0), sind(60d0), 0d0, 0d0, 0d0, 0d0]
    u0(5,:) = 1d0*[cosd(-60d0), sind(-60d0), 0d0, 0d0, 0d0, 0d0]
    
    do i = 1, 5
        write(*,*) "U inicial", u0(i,1:3)
        call Newton(cr3bp_u, u0(i,1:3))
        write(*,*) "U critico", u0(i,1:3)
        
        
        call System_matrix(cr3bp_uu, u0(i,:), 0d0, a)
        call Eigenvalues_QR(a,lambda)
        write(*,*) "Lambda", lambda(:)
        write(*,*) ""
        
        call plot_cp(u0(i,:))
    end do
    write(*,*) "Press ENTER"
    read(*,*)
    contains
    
    subroutine plot_cp(u0)
    real, intent(in) :: u0(:)
        real :: u(0:M,N), time(0:M)
        
        time(0:M) = [(i*TF/M, i=0, M)]
        u(0,:) = u0 + 1e-3 * [1d0, 1d0, 0d0, 0d0, 0d0, 0d0]
    
        call cp( &
            time_domain = time, &
            temporal_scheme = runge_kutta_4, &
            differential_operator = cr3bp, &
            solution = u &
        )
        
        call plot(u(:,1) - u0(1), u(:,2) - u0(2), color = "RED", hold_on = .true.)
        !call plot(time, u(:,1), color = "RED", hold_on = .true.)
        !call plot(time, u(:,2), color = "GREEN", hold_on = .true.)
        !call plot(time, u(:,3), color = "BLUE", hold_on = .true.)
        call plot_legend(["PLOT"])
        call plot_end()
    end subroutine
    
    !write(*,*) u0
    !call System_matrix(cr3bp, u0, 1d-3, a)
    !call Eigenvalues_QR(a,lambda)
    
    !time(0:M) = [(i*TF/M, i=0, M)]
    !u(0,:) = u0
    
    !call cauchy_problem( &
    !    time_domain = time, &
    !    temporal_scheme = runge_kutta_4, &
    !    differential_operator = cr3bp, &
    !    solution = u &
    !)
    
    !call plot(time, u(:,1), color = "RED", hold_on = .true.)
    !call plot(time, u(:,2), color = "GREEN", hold_on = .true.)
    !call plot(time, u(:,3), color = "BLUE", hold_on = .true.)
    !call plot_legend(["NONE", "X", "Y", "Z"])
    !call plot_end()
    
end program restricted_3_body_problem