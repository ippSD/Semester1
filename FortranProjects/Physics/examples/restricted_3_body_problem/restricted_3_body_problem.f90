program restricted_3_body_problem
    use orbit_functions, only: cr3bp, cr3bp_u, cr3bp_uu
    use cauchy_problem_solver, only: cp => cauchy_problem
    use Cauchy_Problem, only: System_matrix
    use temporal_schemes_my, only: runge_kutta_4, euler_explicit
    use Non_Linear_Systems
    use Numerical_Recipes
    use dislin_mod
    implicit none
    
    integer, parameter :: N = 6, M = 10000
    real, parameter :: TF = 2d0 * acos(-1d0) * 4, MU = (0.049d5)/( 0.049d5+3.986d5), EPSIL = 1e-6
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
        
        call plot_cp(u0(i,:), a)
    end do
    write(*,*) "Press ENTER"
    read(*,*)
    contains
    
    subroutine plot_cp(u0,a)
        real, intent(in) :: u0(:), a(:,:)
        real :: u(0:M,N), time(0:M)
        real :: w(0:M,N), w0(N), xmin, xmax, ymin, ymax
        
        time(0:M) = [(i*TF/M, i=0, M)]
        call random_seed()
        call random_number(w0)
        !w0 = 0d0;
        u(0,:) = u0 + EPSIL * w0
        w(0,:) = w0
    
        call cp( &
            time_domain = time, &
            temporal_scheme = runge_kutta_4, &
            differential_operator = cr3bp, &
            solution = u &
        )
            
        call cp( &
            time_domain = time, &
            temporal_scheme = euler_explicit, &
            differential_operator = cr3bp_lineal, &
            solution = w &
        )
            
        xmin = minval([minval(u(:,1) - u0(1)), minval(EPSIL*w(:,1))]);
        xmax = maxval([maxval(u(:,1) - u0(1)), maxval(EPSIL*w(:,1))]);
        ymin = minval([minval(u(:,2) - u0(2)), minval(EPSIL*w(:,2))]);
        ymax = maxval([maxval(u(:,2) - u0(2)), maxval(EPSIL*w(:,2))]);
        
        call plot([xmin,xmax], [ymin,ymax], color = "BLACK", hold_on = .true.)
        
        call plot(u(:,1) - u0(1), u(:,2) - u0(2), color = "RED", hold_on = .true.)
        call plot(EPSIL*w(:,1), EPSIL*w(:,2), color = "BLUE", hold_on = .true.)
        call plot_legend(["PLOT", "LINEAL"])
        call plot_end()
    end subroutine
    
    function cr3bp_lineal(w,t) result(f)
        real, intent(in) :: w(:), t
        real :: f(size(w))
        
        f = matmul(a,w);
    end function
    
end program restricted_3_body_problem