module functions
  implicit none
contains

  function oscilador(u,t) result(f)
    real( kind = 8 ), intent(in) :: u(0:), t
    real( kind = 8 ) :: f(0:size(u) - 1)

    f = [u(1), -u(0)]
  end function

  function kepler(u,t) result(f)
    real( kind = 8 ), intent(in) :: u(0:), t
    real( kind = 8 ) :: f(0:size(u))

    f = [u(2), u(3), -u(0)/norm2(u(0:1)), -u(1)/norm2(u(0:1))]
  end function
  
  function n_bodies(u,t) result(f)
    real( kind = 8 ), intent(in) :: u(:)
    real( kind = 8 ), intent(in) :: t
    real( kind = 8 ) :: f(size(u))
    real( kind = 8 ) :: r_i(3), r_j(3), a_i(3)
    double precision, allocatable :: mu(:)
    integer :: i, j, n
    
    n = size(u) / 7
    allocate(mu(n))
    mu = u(1:n)
    
    ! u = [(mu) * N, (x,y,z) * N, (vx,vy,vz) * N] = 7 * N
    ! u = [1:N     , N+1:4N     ,4N+1:7N        ] = 7 * N
    
    
    f(1:n) = [ (0e0,i=1,n)]  ! Mass is constant
    f(n+1:4*n) = u(4*n+1:7*n)  ! dx/dt = v
    
    do i = 0, n-1
        r_i = u(n+1+3*i:n+4*i)
        a_i = [0e0,0e0,0e0]
        do j = 0, n-1
            r_j = u(n+j+3*j:n+4*j)
            if (i/=j) then
                a_i = a_i - mu(j+1) * (r_i - r_j) / norm2(r_i-r_j)**3e0
            end if
        end do
        f(4*n+1+3*i:4*n+4*i) = a_i
        
    end do
  end function
  
  subroutine sub_n_bodies( t, m, u, u_next)
    integer, intent(in) :: m
    real( kind = 8 ), intent(in) :: t, u(1:m)
    real( kind = 8 ), intent(out) :: u_next(size(u))
    
    u_next = n_bodies(u,t)
  end subroutine
  
end module
