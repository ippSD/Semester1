program $3bodies
    use functions
    use numerical_methods
    use cauchy
    use dislin
    implicit none
    
    integer, parameter :: m = 2000, s = 3
    integer, parameter :: n = s * 7
    double precision, dimension(s) :: mu =  [ 8d1,1d0,0d0]
    double precision, parameter :: v0m_s = 8d1/81d0, v0e_s = 1d0/81d0, v0s_s = 3d0

    double precision, dimension(3) :: x0e = [-1d0, 0d0 ,0d0]
    double precision, dimension(3) :: v0e = [ 0d0,-v0e_s, 0d0]
    double precision, dimension(3) :: x0m = [ 8d1, 0d0 ,0d0]
    double precision, dimension(3) :: v0m = [ 0d0, v0m_s ,0d0]
    double precision, dimension(3) :: x0s = [ 0d0, 9d0 ,0d0]
    double precision, dimension(3) :: v0s = [ v0s_s, 0d0 ,0d0]
    
    double precision :: tf, time(0:m), u(0:m,1:n), x_cg(0:m,3), v_cg(0:m,3)
    integer :: i, k

    tf = 1d3
    time = [(tf/m*i,i=0,m)]
    u(0,:) = [mu, x0e, x0m, x0s, v0e, v0m, v0s]
    
    
    call cauchy_problem(time, n_bodies, runge_kutta, u)
    
    
    open(unit=13, file="out.txt", status="REPLACE")
    do i = 0, m
        write(13,"(9F)") u(i, s+1:4*s)
    end do
    close(13)
        
    x_cg = (mu(1)*u(:,  s+1:  s+3) + mu(2)*u(:,  s+4:  s+6)) / sum(mu)
    v_cg = (mu(1)*u(:,4*s+1:4*s+3) + mu(2)*u(:,4*s+4:4*s+6)) / sum(mu)
    call qplot(u(:,s+1), u(:,s+2), m)  ! Earth
    call qplot(u(:,s+4), u(:,s+5), m)  ! Moon
    call qplot(u(:,s+7), u(:,s+8), m)  ! S/C
    call qplot(time, x_cg(:,3), m)  ! Zcg
    call qplot(time, v_cg(:,1), m)  ! VXcg
    
end program $3bodies
    !time = [(i*dt,i=1,iters)]
    !x_cg = (mu(1)*u(  n+1:  n+3,:) + mu(2)*u(  n+4:  n+6,:)) / sum(mu)
    !v_cg = (mu(1)*u(4*n+1:4*n+3,:) + mu(2)*u(4*n+4:4*n+6,:)) / sum(mu)
        
    !call qplot(u(n+1,:), u(n+2,:), iters)  ! Earth
    !call qplot(u(n+4,:), u(n+5,:), iters)  ! Moon
    !call qplot(u(n+7,:), u(n+8,:), iters)  ! S/C
    !call qplot(time, x_cg(3,:), iters)  ! Zcg
    !call qplot(time, v_cg(1,:), iters)  ! VXcg
    
    
    !call poincare_map(u(n+1  :n+3  ,:), v0e_s, iters)
    !call poincare_map(u(n+1+3:n+3+3,:), v0m_s, iters)
    !call poincare_map(u(n+1+6:n+3+6,:), v0s_s, iters)
    
    !call lagrange(mu(1:2), u(  n+1:  n+3,:), u(  n+4:  n+6,:), u(4*n+1:4*n+3,:), u(4*n+4:4*n+6,:), iters)
    
    !contains
    
    !subroutine poincare_map(u, param, m)
    !    integer, intent(in) :: m
    !    double precision, intent(in) :: u(3,m)
    !    double precision, intent(in) :: param
        
    !    double precision, dimension(3) :: nor = [0d0,1d0,0d0]
    !    double precision, dimension(3) :: p1, p2, r_c
    !    double precision, allocatable :: r_pi(:,:), r_pi_2(:,:)
    !    double precision :: lambda
    !    integer :: i, k
        
    !    k = 1
    !    allocate(r_pi(3,k), r_pi_2(3,k))
    !    r_pi = 0
    !    r_pi_2 = 0
    !    do i = 1, m - 2
    !        if (dot_product(nor, u(:,i))*dot_product(nor, u(:,i+1)) < 0d0 .and. u(1,i)>=0d0) then
    !            p1 = u(:,i)
    !            p2 = u(:,i+1)
    !            lambda = -dot_product(nor,p1) / dot_product(nor, p2-p1)
    !            r_c = p1 + lambda * (p2 - p1)
            
    !            deallocate(r_pi)
    !            allocate(r_pi(3,k+1))
            
    !            r_pi(:,1:k) = r_pi_2
    !            r_pi(:,k+1) = r_c
            
    !            deallocate(r_pi_2)
    !            allocate(r_pi_2(3,k+1))
            
    !            r_pi_2 = r_pi
    !            k = k + 1
    !        end if
    !    end do
    
    !    call QPLSCA([(param,i=1,k-1)], r_pi(1,2:k), k-1)
    !    deallocate(r_pi, r_pi_2)
    !end subroutine poincare_map
    
    !subroutine lagrange(mu, x1, x2, v1, v2, m)
    !    integer, intent(in) :: m
    !    double precision, intent(in) :: x1(3,m), x2(3,m), v1(3,m), v2(3,m), mu(2)
        
    !    double precision, dimension(3,m) :: xcg, x1_rel, x2_rel
    !    double precision :: r(3,m), a(3,m), de(3), dm(3), alpha, w, y
    !    integer :: i, k
        
    !    xcg = (mu(1) * x1 + mu(2) * x2) / sum(mu)
    !    x1_rel = x1 - xcg
    !    x2_rel = x2 - xcg
    !    w = norm2(v2(:,1))/norm2(x2_rel(:,1))
    !    alpha = 0d0
    !    do i = 1, m
    !        y = 2d0 * i / 3d0
    !        if(i >= m / 2d0) then
    !            y = y + m/3d0
    !        end if
    !        r(:,i) = (y - m/2) * 2d0 / m * 100d0 * [cosd(alpha), sind(alpha), 0d0]
    !        dm = r(:,i) - 80d0 * [1d0,0d0,0d0]
    !        de = r(:,i) - 1d0 * [-1d0,0d0,0d0]
    !        a(:,i) = -mu(2)*dm/norm2(dm)**3d0 - mu(1)*de/norm2(de)**3d0 - w**2d0*r(:,i)
    !    end do
        
    !    call QPLSCA(r(1,:), a(1,:), m)
        
    !end subroutine lagrange
        

