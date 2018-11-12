program prog_n_bodies
    use functions
    use numerical_methods
    use dislin
    implicit none
    
    integer, parameter :: iters = 6000, n = 3
    integer, parameter :: s = n * 7
    double precision, dimension(n) :: mu =  [ 8d1, 1d0,0d0]
    double precision, parameter :: v0m_s = 8d1/81d0, v0e_s = 1d0/81d0, v0s_s = 3d0

    double precision, dimension(3) :: x0e = [-1d0, 0d0 ,0d0]
    double precision, dimension(3) :: v0e = [ 0d0,-v0e_s, 0d0]
    double precision, dimension(3) :: x0m = [ 8d1, 0d0 ,0d0]
    double precision, dimension(3) :: v0m = [ 0d0, v0m_s ,0d0]
    double precision, dimension(3) :: x0s = [ 0d0, 9d0 ,0d0]
    double precision, dimension(3) :: v0s = [ v0s_s, 0d0 ,0d0]
    
    double precision :: u(s, iters), x_cg(3, iters), v_cg(3, iters), time(iters)
    double precision :: dt = 1d-1
    integer :: i, k

    u(:,1) = [mu, x0e, x0m, x0s, v0e, v0m, v0s]
    
    
    open(unit=13, file="out.txt", status="REPLACE")
    do i = 1, iters - 1
        call runge_kutta ( i*dt, s, u(:,i), dt, sub_n_bodies, u(:,i+1) )
        write(13,"(9F)") u(n+1:4*n, i)
    end do
    close(13)
    
    time = [(i*dt,i=1,iters)]
    x_cg = (mu(1)*u(  n+1:  n+3,:) + mu(2)*u(  n+4:  n+6,:)) / sum(mu)
    v_cg = (mu(1)*u(4*n+1:4*n+3,:) + mu(2)*u(4*n+4:4*n+6,:)) / sum(mu)
        
    call qplot(u(n+1,:), u(n+2,:), iters)  ! Earth
    call qplot(u(n+4,:), u(n+5,:), iters)  ! Moon
    call qplot(u(n+7,:), u(n+8,:), iters)  ! S/C
    call qplot(time, x_cg(3,:), iters)  ! Zcg
    call qplot(time, v_cg(1,:), iters)  ! VXcg
    
    
    call poincare_map(u(n+1  :n+3  ,:), v0e_s, iters)
    call poincare_map(u(n+1+3:n+3+3,:), v0m_s, iters)
    call poincare_map(u(n+1+6:n+3+6,:), v0s_s, iters)
    
    contains
    
    subroutine poincare_map(u, param, m)
        integer, intent(in) :: m
        double precision, intent(in) :: u(3,m)
        double precision, intent(in) :: param
        
        double precision, dimension(3) :: nor = [0d0,1d0,0d0]
        double precision, dimension(3) :: p1, p2, r_c
        double precision, allocatable :: r_pi(:,:), r_pi_2(:,:)
        double precision :: lambda
        integer :: i, k
        
        k = 1
        allocate(r_pi(3,k), r_pi_2(3,k))
        r_pi = 0
        r_pi_2 = 0
        do i = 1, m - 2
            if (dot_product(nor, u(:,i))*dot_product(nor, u(:,i+1)) < 0d0 .and. u(1,i)>=0d0) then
                p1 = u(:,i)
                p2 = u(:,i+1)
                lambda = -dot_product(nor,p1) / dot_product(nor, p2-p1)
                r_c = p1 + lambda * (p2 - p1)
            
                deallocate(r_pi)
                allocate(r_pi(3,k+1))
            
                r_pi(:,1:k) = r_pi_2
                r_pi(:,k+1) = r_c
            
                deallocate(r_pi_2)
                allocate(r_pi_2(3,k+1))
            
                r_pi_2 = r_pi
                k = k + 1
            end if
        end do
    
        call QPLSCA([(param,i=1,k-1)], r_pi(1,2:k), k-1)
        deallocate(r_pi, r_pi_2)
    end subroutine poincare_map
        
        
end program
