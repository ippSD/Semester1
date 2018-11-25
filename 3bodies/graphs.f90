module graphs
    use dislin
    implicit none
    
    contains
    
!!! CHECK_STATE_VECTOR
!!! Checks if acquired state vector is compatible with
!!! required input -> u = u(7N) = N mases + (x,y,z) x N + (vx,vy,vz) x N
!!! 
!!! INPUT:
!!!     u(m,n) state vector over 'm' time steps with 'n' elements.
!!! OUTPUT:
!!!     logical _ : true if n == 7N, else false.
    logical function check_state_vector(u) result(f)
        double precision, intent(in) :: u(:,:)  ! u(m,n)
        integer :: n
        
        n = size(u(1,:))
        
        if (mod(n, 7) /= 0) then
            write(*,*) "Error: invalid state vector."
            write(*,*) "State vector should be size '7N' with"
            write(*,*) "'N' being the number of bodies"
            f = .false.
        else
            f = .true.
        end if
    end function check_state_vector

!!! PLOT_ORBIT_XY
!!! Plot XY orbit of the 'l'-th body from the state vector.
!!! 
!!! INPUT:
!!!     u(m,n) : state vector over 'm' time steps with 'n' elements.
!!!     l : body to be plotted. 0 < l <= n/7
!!! OUTPUT:
!!!     None
    subroutine plot_orbit_xy(u, l)
        integer, intent(in) :: l
        double precision, intent(in) :: u(:,:)  ! u(m,n)
        
        integer :: m, n, nbodies, x_idx
        
        m = size(u(:,1))
        n = size(u(1,:))
        
        if (.not. check_state_vector(u)) then
            return
        else if(n == 0 .or. l > n/7) then
            write(*,*) "Error: invalid body number."
            write(*,*) "Number must be grater than 0 and lower or equal to ", n/7
            return
        end if
        
        nbodies = n/7
        x_idx = nbodies + (l-1) * 3       + 1
        !       mases   + previous bodies + X coor.
        
        call qplot(u(:, x_idx), u(:, x_idx + 1), m)
        !call qplot(time, v_cg(:,1), m)  ! VXcg
    end subroutine plot_orbit_xy

!!! PLOT_CDG_XY
!!! Plot XY gravity center position.
!!! 
!!! INPUT:
!!!     u(m,n) : state vector over 'm' time steps with 'n' elements.
!!! OUTPUT:
!!!     None
    subroutine plot_cdg_xy(u)
        double precision, intent(in) :: u(:,:)  ! u(m,n)
        
        integer :: l, m, n, nbodies, x_idx
        double precision, allocatable :: cdg(:,:)
        
        m = size(u(:,1))
        n = size(u(1,:))
        
        if (.not. check_state_vector(u)) then
            return
        end if
        
        nbodies = n / 7
        allocate(cdg(m,3))
        cdg = 0d0
        
        do l = 1, nbodies
            x_idx = nbodies + (l-1)*3 + 1
            cdg = cdg + u(1,l) * u(:,x_idx:x_idx+2) / sum(u(1,1:nbodies))
        end do
        
        call qplot(cdg(:,1), cdg(:,2), m)
    end subroutine plot_cdg_xy

!!! POINCARE_MAP
!!! Plot the Poincaré Map from data.
!!! 
!!! INPUT:
!!!     dat(m,2) : inicial condition parameter (dat(:,1)) and
!!!                 intersection points of orbit with Poincaré plane (dat(:,2))
!!! OUTPUT:
!!!     None        
    subroutine poincare_map(dat)
        double precision, intent(in) :: dat(:,:)
        
        integer :: m
        
        m = size(dat(:,1))
        call qplsca(dat(:,1), dat(:,2), m)
    end subroutine poincare_map
    
end module
    