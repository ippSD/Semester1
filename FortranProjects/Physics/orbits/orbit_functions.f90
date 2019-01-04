module orbit_functions
    implicit none
    contains
    
    function f_n_bodies(u, t) result(f)
        real, intent(in) :: u(:)
        real, intent(in) :: t
        real :: f(size(u))
    
        real, pointer :: r_i(:), r_j(:), mu_j, v_i(:)
        real, pointer :: dmdt_i, drdt_i(:), dvdt_i(:)
        integer :: i, j
        integer :: n
        
        n = size(u) / 7
        
        ! Mass is constant: dmdt = 0 forall body
        do i = 1, n
            call p2m(f, i, dmdt_i)
            dmdt_i = 0d0
        end do
        
        ! drdt = v
        do i = 1, n
            call p2v(u, i, v_i)
            call p2r(f, i, drdt_i)
            drdt_i = v_i
        end do
        
        ! dvdt = Fg
        do i = 1, n
            call p2v(f, i, dvdt_i)
            call p2r(u, i, r_i)
            
            dvdt_i = 0d0
            do j = 1, n
                call p2m(u, j, mu_j)
                call p2r(u, j,  r_j)
                if (i/=j) then
                    dvdt_i = dvdt_i - mu_j * (r_i - r_j) / norm2(r_i-r_j)**3d0
                end if
            end do
        end do
        
    contains
    
        ! u = [(mu) * N, (x,y,z) * N, (vx,vy,vz) * N] = 7 * N
        ! u = [1:N     , N+1:4N     ,4N+1:7N        ] = 7 * N
    
        subroutine p2m(state_vector, l, m_pointer)
            real, target, intent(in) :: state_vector(:)
            integer, intent(in) :: l
            real, pointer, intent(out) :: m_pointer
            m_pointer => state_vector(l)
        end subroutine
    
        subroutine p2r(state_vector, l, r_pointer)
            real, target, intent(in) :: state_vector(:)
            integer, intent(in) :: l
            real, pointer, intent(out) :: r_pointer(:)
        
            integer :: n
        
            n = size(state_vector) / 7
            r_pointer => state_vector(n+3*(l-1)+1:n+3*(l-1)+3)
        end subroutine
    
        subroutine p2v(state_vector, l, v_pointer)
            real, target, intent(in) :: state_vector(:)
            integer, intent(in) :: l
            real, pointer, intent(out) :: v_pointer(:)
        
            integer :: n
        
            n = size(state_vector) / 7
            v_pointer => state_vector(4*n+3*(l-1)+1:4*n+3*(l-1)+3)
        end subroutine
    
    end function
    
    subroutine sub_n_bodies( t, m, u, u_next)
        integer, intent(in) :: m
        real, intent(in) :: t, u(1:m)
        real, intent(out) :: u_next(size(u))
        
        u_next = f_n_bodies(u,t)
    end subroutine

end module orbit_functions
