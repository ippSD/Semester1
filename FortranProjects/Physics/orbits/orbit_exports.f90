module orbit_exports
    use orbit_check
    implicit none
    contains
    
    subroutine save_orbits(unit, u)
        integer, intent(in) :: unit
        real, intent(in) :: u(0:,:)
        
        integer :: i, m, n, s, idx_i, idx_f
        
        if (.not. check_state_matrix(u)) return
        n = size(u(0,:))
        m = size(u(:,1)) - 1
        s = n / 7
        
        idx_i = s + 1
        idx_f = 4 * s
        do i = 0, m
            write(unit,*) u(i, idx_i:idx_f)
        end do
    end subroutine
    
    subroutine save_body_orbit(unit, u, l)
        integer, intent(in) :: unit, l
        real, intent(in) :: u(0:,:)
        
        integer :: i, m, n, s, idx_i, idx_f
        
        if (.not. check_state_matrix(u)) return
        if (.not. check_existing_body(u(0,:), l)) return
        n = size(u(0,:))
        m = size(u(:,1)) - 1
        s = n / 7
        
        idx_i = s + (l-1) * 3 + 1
        idx_f = s + (l-1) * 3 + 3
        do i = 0, m
            write(unit,"(3F)") u(i, idx_i:idx_f)
        end do
    end subroutine
    
end module orbit_exports