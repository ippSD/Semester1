module solver
    implicit none
    
    interface
    
        subroutine g_sol(alpha, m, x, y)
            integer, intent(in) :: m
            double precision, intent(in) :: alpha(:), x(m)
            double precision, intent(out) :: y(m)
        end procedure g_sol
    
    contains
    
    subroutine jac(alpha, eps, g, u_iter, jac)
        
    end subroutine
    
    subroutine sis_solver(alpha, m, g, u0, usol)
        double precision, intent(in) :: alpha(:)
        integer, intent(in) :: m
        procedure(g_sol) :: g
        double precision, intent(in) :: u0(m)
        double precision, intent(out) :: usol(m)
        
        integer :: i
        double precision :: u_iter(m),  jac(m,m), ident(m,m) = 0e0
        
        do i = 1, m
            ident(i,i) = 1e0
        end do
        
        do i = 1, 10
            call jacobian(alpha, eps, g, u_iter, jac)
            usol = (ident - jac) * u_iter
            u_iter = usol
        end do
        
end module solver