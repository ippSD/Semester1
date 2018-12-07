module auxiliary_functions
    implicit none
    interface operator (.vec.)
        module procedure times_r3
    end interface
    
    contains
    
    function times_r3(v1, v2) result(v3)
        real, intent(in) :: v1(3), v2(3)
        real :: v3(3)
        
        v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
        v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
        v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
    end function times_r3
end module auxiliary_functions