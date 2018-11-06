program prog_n_bodies
    use functions
    use numerical_methods
    implicit none
    
    real( kind = 8 ), dimension(2) :: mu = [8e1,1e0]
    real( kind = 8 ), dimension(3) :: x0e = [-1e0,0e0,0e0], &
        v0e = [-1e0,0e0,0e0], &
        x0m = [80e0,0e0,0e0], &
        v0m = [80e0,0e0,0e0]
    real( kind = 8 ),dimension(14) :: u0, u
    real( kind = 8 ) :: dt = 1e-3
    integer :: i, iters, n
    
    n = 2
    u0 = [8e1,1e0,-1e0,0e0,0e0,80e0,0e0,0e0,-1e0,0e0,0e0,80e0,0e0,0e0]
    u = u0
    iters = 500
    
    open(unit=13, file="out.txt", status="REPLACE")
    
    do i = 0, iters
        call runge_kutta ( i*dt, 14, u0, dt, sub_n_bodies, u )
        write(13,"(6F)") u(n+1:4*n)
        u0 = u
    end do
    
    close(13)
        
end program
