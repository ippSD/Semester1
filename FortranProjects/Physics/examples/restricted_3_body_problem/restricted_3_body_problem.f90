program restricted_3_body_problem
    use math, only: jacobian_matrix
    use orbit_functions, only: cr3bp
    use ode_interfaces, only: propagator => temporal_scheme
    use cauchy_problem_solver, only: cp => cauchy_problem
    use Cauchy_Problem, only: System_matrix
    use temporal_schemes, only: dopri853, ode, odex
    use Non_Linear_Systems, only: Newton
    use Numerical_Recipes, only: Eigenvalues_QR
    use dislin_mod
    implicit none
    
    integer, parameter :: N = 6, M = 100
    real, parameter :: PI = acos(-1d0)
    real, parameter :: TFS(5) = 2d0 * PI * [5e-1, 5e-1, 1e0, 5e0, 5e0]
    real, parameter :: MU_E = 3.986d5, MU_M = 0.049d5
    real, parameter :: MU = 1d0 / ( 1d0 + MU_E / MU_M )
    real, parameter :: EPSIL = 1e-4
    
    integer :: i, j
    
    real, target :: u0_lagrange_points(5,N)
    real, pointer :: u0(:)
    real :: u(0:M,N), time(0:M), a(N,N)
    real :: w(0:M,N), w0(N)
    
    procedure(propagator), pointer :: temporal_scheme
    
    complex :: lambda(N)
    
    character(len = 100) :: filename
    
    temporal_scheme => dopri853
        
    ! Seed for Newton
    u0_lagrange_points(1,:) = 1d0*[5d-1, 0d0, 0d0, 0d0, 0d0, 0d0]
    u0_lagrange_points(2,:) = 1d0*[1d0, 0d0, 0d0, 0d0, 0d0, 0d0]
    u0_lagrange_points(3,:) = 1d0*[-1d0, 0d0, 0d0, 0d0, 0d0, 0d0]
    u0_lagrange_points(4,:) = 1d0*[cosd(+60d0), sind(+60d0), 0d0, 0d0, 0d0, 0d0]
    u0_lagrange_points(5,:) = 1d0*[cosd(-60d0), sind(-60d0), 0d0, 0d0, 0d0, 0d0]
    
    do i = 1, 5
        ! Init variables
        u0 => u0_lagrange_points(i,:)
        write(filename,"(1I,4A)") i, ".dat"
        write(*,"(A,I1)") "Punto Lagrange L", i
        
        ! Find Zeros
        call get_zeros(u0(1:3))
        ! System matrix around a zero
        call get_sysmatrix(u0, a)
        ! Eigenvalues
        call get_eigs(a, lambda)
        ! Propagate body
        call get_transients(u0, EPSIL, TFS(i), u, w)
        ! Export data to file
        call export_data(u, u0, w, EPSIL, filename)
        ! Plot results
        call plot_disturbance( &
            u(:,1) - u0(1), &
            EPSIL * w(:,1), &
            u(:,2) - u0(2), &
            EPSIL * w(:,2) &
        )
        nullify(u0)
    end do
    ! End
    write(*,*) "Press ENTER"
    read(*,*)
    
    contains
    
    subroutine get_zeros(u0)
        real, intent(inout) :: u0(:)
        write(*,"(4A,3E12.4)") char(9),  "U inicial:", char(10), char(9), u0
        call Newton(cr3bp_zeros_calculation, u0)
        write(*,"(4A,3E12.4)") char(9), "U solucion:", char(10), char(9), u0
    end subroutine
    
    subroutine get_sysmatrix(u0, a)
        real, intent(in) :: u0(:)
        real, intent(out) :: a(size(u0), size(u0))
        
        ! 2 Ways for calculation
        !call System_matrix(cr3bp_systemmatrix_compatible, u0, 0d0, a)
        call jacobian_matrix(cr3bp_notime, u0, a, 1d-5)
        
        ! Write output to term
        write(*,*) char(9), "Matriz del sistema"
        write(*,"(2A,6D12.4,A)") char(9), "X0 = [", u0, " ]"
        call write_matrix(a,1)
    end subroutine get_sysmatrix
    
    subroutine get_eigs(a, lambdas)
        real, intent(in) :: a(:,:)
        complex, intent(out) :: lambdas(size(a(1,:)))
        
        integer :: i
        
        call Eigenvalues_QR(a,lambdas)
        write(*,"(2A)") char(9), "Lambdas:"
        print_eigs: do i = 1, 6
            write(*,"(A,2E12.4)") char(9), real(lambdas(i)), imag(lambdas(i));
        end do print_eigs
    end subroutine get_eigs
    
    subroutine get_transients(u0, eps, tf, u, w)
        real, intent(in) :: tf, u0(N), eps
        real, intent(out) :: u(0:M,N), w(0:M,N)
        
        integer :: i
        real :: time(0:M)
        
        time = [(i * tf / M, i = 0, M)]
        call random_seed()
        call random_number(w0)
        w(0,:) = w0
        u(0,:) = u0 + eps * w(0,:)
        
        ! Propagate
        call cp( &
            time_domain = time, &
            temporal_scheme = temporal_scheme, &
            differential_operator = cr3bp, &
            solution = u &
        )
        
        call cp( &
            time_domain = time, &
            temporal_scheme = temporal_scheme, &
            differential_operator = cr3bp_perturbacion_lineal, &
            solution = w &
        )
    end subroutine
    
    subroutine export_data(u, u0, w, eps, filename)
        character(len=*), intent(in) :: filename
        real, intent(in) :: u(:,:), u0(:), w(:,:), eps
        
        character(len=30) :: fmt
        real :: linear_disturbance(2), nonlinear_disturbance(2)
        integer :: i, m
        
        m = size(u0)
        write(fmt, "(A)") "(3(D12.4,A1),D12.4)"
        write(*,*) "Export"
        open(unit = 13, file = trim(filename))
        export: do i = 1, M
            nonlinear_disturbance = u(i,1:2) - u0(1:2)
            linear_disturbance = eps * w(i,1:2)
            !write(*,*) u(i,1), u0(1)
            write(13, fmt) &
                nonlinear_disturbance(1), char(9), &
                nonlinear_disturbance(2), char(9), &
                linear_disturbance(1), char(9), &
                linear_disturbance(2)
        end do export
        close(13);
    end subroutine export_data
    
    subroutine write_matrix(a, s)
        real, intent(in) :: a(:,:)
        integer, intent(in), optional :: s
        
        character(len=100) :: fmt
        integer :: i, j, m, n, tabs = 0
        if (present(s)) tabs = s
        m = size(a(:,1))
        n = size(a(1,:))
        write(fmt,"(A,I2,A)") "(A,", n, "D12.4,A)"
        
        do j = 1, tabs; write(*,"(A1)", advance = "no") char(9); enddo
        write(*,*) "A = ["
        do i = 1, m
            do j = 1, tabs + 1; write(*,"(A1)", advance = "no") char(9); enddo
            write(*,fmt) "[ ", a(i,:), "]"
        end do
        do j = 1, tabs; write(*,"(A1)", advance = "no") char(9); enddo
        write(*,*) "    ]"
    end subroutine write_matrix
    
    subroutine plot_disturbance(x_nl, x_l, y_nl, y_l)
        real, intent(in) :: x_nl(:), x_l(:), y_nl(:), y_l(:)
        real :: xmin, xmax, ymin, ymax
        
        xmin = minval([x_nl, x_l])
        xmax = maxval([x_nl, x_l])
        ymin = minval([y_nl, y_l])
        ymax = maxval([y_nl, y_l])

        !xmin = minval([minval(u(:,1) - u0(1)), minval(EPSIL*w(:,1))]);
        !xmax = maxval([maxval(u(:,1) - u0(1)), maxval(EPSIL*w(:,1))]);
        !ymin = minval([minval(u(:,2) - u0(2)), minval(EPSIL*w(:,2))]);
        !ymax = maxval([maxval(u(:,2) - u0(2)), maxval(EPSIL*w(:,2))]);
        
        !call plot([xmin,xmax], [ymin,ymax], color = "BLACK", hold_on = .true.)
        
        !call plot(time, u(:,1), color = "RED")
        !call plot(time, u(:,2), color = "RED")
        
        !call plot(u(:,1) - u0(1), u(:,2) - u0(2), color =  "RED", hold_on = .true.)
        call plot([xmin, xmax], [ymin, ymax], color = "BLACK", hold_on = .true.)
        call plot(x_nl, y_nl, color =  "RED", hold_on = .true.)
        call plot( x_l,  y_l, color = "BLUE", hold_on = .true.)
        call plot_legend(["NONE", "NON-LINEAL", "LINEAL"])
        call plot_end()
                
    end subroutine
    
    function cr3bp_notime(u) result(f)
        real, intent(in) :: u(:)
        real :: f(size(u))
        
        f = cr3bp(u, 0d0)
    end function cr3bp_notime
    
    function cr3bp_zeros_calculation(p) result(f)
        real, intent(in) :: p(:)
        real :: f(size(p))
        real :: f2(6)
        
        f2 = cr3bp([p, 0d0, 0d0, 0d0], 0d0)
        f = f2(4:6)
    end function cr3bp_zeros_calculation
    
    function cr3bp_perturbacion_lineal(w,t) result(f)
        real, intent(in) :: w(:), t
        real :: f(size(w))
        
        f = matmul(a,w);
    end function
    
    function cr3bp_systemmatrix_compatible(u, t) result(f)
        real :: u(:), t
        real :: f(size(u))
        
        f = cr3bp(u, t)
    end function
    
end program restricted_3_body_problem