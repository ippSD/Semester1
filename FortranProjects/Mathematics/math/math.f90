module math
    !use Jacobian_module
    !use Linear_systems
    !use Non_Linear_Systems
    implicit none
    
    abstract interface
        function f_rn(x) result(f)
            real, intent(in) :: x(:)
            real :: f(size(x))
        end function f_rn
    end interface
    
    
    !-----------------------------------------------------------------------!
    !   ( INTERFACE ) ( OPERATOR ) .times.                                  !
    !-----------------------------------------------------------------------!
    !   Operates the cross product of a pair of R3 vectors.                 !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(3)) v1                                          !
    !               ! First R3 vector.                                      !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(3)) v2                                          !
    !               ! Second R3 vector.                                     !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(3)) v3                                          !
    !               ! Cross product of v1 and v2: v3 = v1 x v2.             !
    !---------------!-------------------------------------------------------!
    interface operator (.times.)
        module procedure times_r3
    end interface
    
    private
    public :: operator (.times.)
    public :: jacobian_matrix
    
    contains
    
    !-----------------------------------------------------------------------!
    !   ( FUNCTION ) times_r3                                               !
    !-----------------------------------------------------------------------!
    !   Implementation of the cross product of a pair of R3 vectors.        !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(3)) v1                                          !
    !               ! First R3 vector.                                      !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(3)) v2                                          !
    !               ! Second R3 vector.                                     !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(3)) v3                                          !
    !               ! Cross product of v1 and v2: v3 = v1 x v2.             !
    !---------------!-------------------------------------------------------!
    function times_r3(v1, v2) result(v3)
        real, intent(in) :: v1(3), v2(3)
        real :: v3(3)
        
        v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
        v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
        v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
    end function times_r3
    
    subroutine jacobian_matrix(f, x0, a, dx)
        procedure(f_rn) :: f
        real, intent(in) :: x0(:)
        real, intent(out) :: a(size(x0), size(x0))
        real, intent(in), optional :: dx
        
        integer :: i, j, n
        real :: h = 1d-3
        
        n = size(x0)
        if ( present(dx) ) h = dx
        
        loop_column: do i = 1, n
            a(:,i) = (f(x0 + h * [(kronecker_delta(j, i), j = 1, n)]) - f(x0)) / h
        end do loop_column
    end subroutine jacobian_matrix
    
    function kronecker_delta(i,j) result(k)
        integer, intent(in) :: i, j
        integer :: k
        k = 0
        if ( i == j) k = 1
    end function kronecker_delta
    
end module math