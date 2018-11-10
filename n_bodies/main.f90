program prog_n_bodies
    use functions
    use numerical_methods
    use dislin
    implicit none
    
    integer, parameter :: iters = 10000, n = 3
    integer, parameter :: s = n * 7
    double precision, dimension(n) :: mu =  [ 8d1, 1d0,0d0]

    double precision, dimension(3) :: x0e = [-1d0, 0d0 ,0d0]
    double precision, dimension(3) :: v0e = [ 0d0,-1d0, 0d0] / 8d1
    double precision, dimension(3) :: x0m = [ 8d1, 0d0 ,0d0]
    double precision, dimension(3) :: v0m = [ 0d0, 1d0 ,0d0]
    double precision, dimension(3) :: x0s = [ 0d0, 9d0 ,0d0]
    double precision, dimension(3) :: v0s = [ 3d0, 0d0 ,0d0]
    
    real(kind = 8) :: u(s, iters), x_cg(3, iters), v_cg(3, iters), time(iters)
    real(kind = 8) :: dt = 1d-1
    integer :: i
    
    u(:,1) = [mu, x0e, x0m, x0s, v0e, v0m, v0s]
    
    open(unit=13, file="out.txt", status="REPLACE")
    
    do i = 1, iters - 1
        call runge_kutta ( i*dt, s, u(:,i), dt, sub_n_bodies, u(:,i+1) )
        write(13,"(9F)") u(n+1:4*n, i)
    end do
    
    time = [(i*dt,i=1,iters)]
    x_cg = (mu(1)*u(  n+1:  n+3,:) + mu(2)*u(  n+4:  n+6,:)) / sum(mu)
    v_cg = (mu(1)*u(4*n+1:4*n+3,:) + mu(2)*u(4*n+4:4*n+6,:)) / sum(mu)
        
    call qplot(u(n+1,:), u(n+2,:), iters)  ! Earth
    call qplot(u(n+4,:), u(n+5,:), iters)  ! Moon
    call qplot(u(n+7,:), u(n+8,:), iters)  ! S/C
    !call qplot(time, x_cg(3,:), iters)  ! Zcg
    !call qplot(time, v_cg(3,:), iters)  ! VZcg
    
    close(13)
        
end program
