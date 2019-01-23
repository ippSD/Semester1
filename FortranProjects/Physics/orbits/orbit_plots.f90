module orbit_plots
    use dislin_mod
    implicit none
    contains
    
    logical function check_state_vector(u) result(f)
        real, intent(in) :: u(:,:)  ! u(m,n)
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
        real, intent(in) :: u(:,:)  ! u(m,n)

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
        
        call plot(u(:, x_idx), u(:, x_idx + 1))
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
        real, intent(in) :: u(:,:)  ! u(m,n)
        
        integer :: l, m, n, nbodies, x_idx
        real, allocatable :: cdg(:,:)
        
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
        
        call plot(cdg(:,1), cdg(:,2))
    end subroutine plot_cdg_xy


end module orbit_plots