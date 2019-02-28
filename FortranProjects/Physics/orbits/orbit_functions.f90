module orbit_functions
    implicit none
    contains
    
    function kepler_2d(u, t) result(f)
        real, intent(in) :: u(:)
        real, intent(in) :: t
        real :: f(size(u))
        
        real :: d
        real, pointer :: r(:), v(:)
        real, pointer :: dr_dt(:), dv_dt(:)
        
        call p2r(u, r); call p2r(f, dr_dt)
        call p2v(u, v); call p2v(f, dv_dt)
        
        ! d = sqrt(x^2+y^2)
        d = norm2(r)
        
        ! dr_dt = v
        dr_dt = v
        
        ! dv_dt = -r/d^2
        dv_dt = -r / d ** 2d0;
        
    contains
    
        ! u = [r = (x, y), v = (vx, vy)] = 4
        ! u = [dr_dt = (dx_dt, dy_dt), dv_dt = (dvx_dt, dvy_dt)] = 4
    
        subroutine p2r(state_vector, r_pointer)
            real, target, intent(in) :: state_vector(:)
            real, pointer, intent(out) :: r_pointer(:)
                
            r_pointer => state_vector(1:2)
        end subroutine
    
        subroutine p2v(state_vector, v_pointer)
            real, target, intent(in) :: state_vector(:)
            real, pointer, intent(out) :: v_pointer(:)
                
            v_pointer => state_vector(3:4)
        end subroutine
    
    end function
    
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
    
    function cr3bp(u, t) result(f)
        real, intent(in) :: u(:)
        real, intent(in) :: t
        real :: f(size(u))
        
        real, pointer :: r(:), v(:), dr_dt(:), dv_dt(:)
        
        real :: mu, d, r_mag
        
        ! Colocar los punteros
        ! Estado
        call p2r(u, r)
        call p2v(u, v)
        
        ! Derivadas
        call p2r(f, dr_dt)
        call p2v(f, dv_dt)
        
        ! Variables
        mu = (0.049d5)/( 0.049d5+3.986d5)
        d = sqrt((r(1)+mu)**2d0 + r(2)**2d0 + r(3)**2d0)
        r_mag = ((r(1)-1d0+mu)**2d0+r(2)**2d0+r(3)**2d0) ** (1d0/3d0)
        
        ! Ecuaciones
        dr_dt = v
        dv_dt(1) = r(1) + 2*v(2) - (1d0 - mu) * (r(1) + mu) / d**3d0 - mu * (r(1) - 1d0 + mu) / r_mag ** 3d0
        dv_dt(2) = r(2) - 2*v(1) - (1d0 - mu) * r(2) / d**3d0 - mu * r(2) / r_mag ** 3d0
        dv_dt(3) = -(1d0-mu)*r(3)/d**3d0 - mu * r(3) / r_mag**3d0
        
    contains
    
        subroutine p2r(state_vector, r_pointer)
            real, target, intent(in) :: state_vector(6)
            real, pointer, intent(out) :: r_pointer(:)
        
            r_pointer(1:3) => state_vector(1:3)
        end subroutine
    
        subroutine p2v(state_vector, v_pointer)
            real, target, intent(in) :: state_vector(6)
            real, pointer, intent(out) :: v_pointer(:)
        
            v_pointer(1:3) => state_vector(4:6)
        end subroutine
    
    end function cr3bp
    
    function cr3bp_u(u) result(f)
        real, intent(in) :: u(:)
        real :: f(size(u))
        real :: f2(6)
        
        f2 = cr3bp([u, 0d0, 0d0, 0d0], 0d0)
        f = f2(4:6)
    end function
    
    function cr3bp_uu(u, t) result(f)
        real :: u(:), t
        real :: f(size(u))
        
        f = cr3bp(u, t)
    end function

end module orbit_functions
