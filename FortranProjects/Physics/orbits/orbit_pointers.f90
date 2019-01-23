module orbit_pointers
    implicit none
    
    contains
    
    subroutine pointer_to_body_position(state_vector, l, r_pointer)
        real, target, intent(in) :: state_vector(:)
        integer, intent(in) :: l
        real, pointer, intent(out) :: r_pointer(:)
        
        integer :: n, lower_idx, upper_idx
        
        n = size(state_vector) / 7
        
        lower_idx = n+3*(l-1)+1
        upper_idx = lower_idx + 2
        r_pointer => state_vector(lower_idx:upper_idx)
    end subroutine
    
    subroutine pointer_to_body_x_vs_time(state_matrix, l, x_pointer)
        real, target, intent(in) :: state_matrix(0:,:)
        integer, intent(in) :: l
        real, pointer, intent(out) :: x_pointer(:)
        
        integer :: n, x_idx
        
        n = size(state_matrix(0,:)) / 7
        
        x_idx = n+3*(l-1)+1
        x_pointer(0:) => state_matrix(0:,x_idx)
    end subroutine
    
    subroutine pointer_to_body_y_vs_time(state_matrix, l, y_pointer)
        real, target, intent(in) :: state_matrix(0:,:)
        integer, intent(in) :: l
        real, pointer, intent(out) :: y_pointer(:)
        
        integer :: n, y_idx
        
        n = size(state_matrix(0,:)) / 7
        
        y_idx = n+3*(l-1)+2
        y_pointer(0:) => state_matrix(0:,y_idx)
    end subroutine
    
    subroutine pointer_to_body_velocity(state_vector, l, v_pointer)
        real, target, intent(in) :: state_vector(:)
        integer, intent(in) :: l
        real, pointer, intent(out) :: v_pointer(:)
        
        integer :: n, lower_idx, upper_idx
        
        n = size(state_vector) / 7
        lower_idx = 4*n+3*(l-1)+1
        upper_idx = lower_idx + 2
        
        v_pointer => state_vector(lower_idx:upper_idx)
    end subroutine
    
    
    
end module orbit_pointers
