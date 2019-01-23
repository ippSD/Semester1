module ode_interfaces
    implicit none
    
    abstract interface
    
        !-------------------------------------------------------------------!
        !   ( INTERFACE ) ( FUNCTION ) ode_function                         !
        !-------------------------------------------------------------------!
        !   Defines an standard interface for Cauchy Problem's derivative   !
        !   function compatible with the temporal schemes and cauchy        !
        !   problem api of this repository.                                 !
        !-------------------------------------------------------------------!
        !   Parameters: !                                                   !
        !---------------!---------------------------------------------------!
        !   IN          ! (real(N)) u                                       !
        !               ! State vector of the cauchy problem at time t.     !
        !---------------!---------------------------------------------------!
        !   IN          ! (real) t                                          !
        !               ! Time at which the state vector are given.         !
        !---------------!---------------------------------------------------!
        !   OUT         ! (real(N)) f                                       !
        !               ! Derivative of the state vector at time t.         !
        !---------------!---------------------------------------------------!
        function ode_function ( u, t ) result( f )
          real, intent(in) :: u(:),  t
          real :: f(size(u))
        end function ode_function
        
        !-------------------------------------------------------------------!
        !   ( INTERFACE ) ( SUBROUTINE ) ode_scalar_subroutine              !
        !-------------------------------------------------------------------!
        !   Defines an standard interface for Cauchy (scalar) Problem's     !
        !   derivative subroutine.                                          !
        !-------------------------------------------------------------------!
        !   Parameters: !                                                   !
        !---------------!---------------------------------------------------!
        !   IN          ! (real) u                                          !
        !               ! State of the cauchy problem at time t.            !
        !---------------!---------------------------------------------------!
        !   IN          ! (real) t                                          !
        !               ! Time at which the state vector are given.         !
        !---------------!---------------------------------------------------!
        !   OUT         ! (real) uprime                                     !
        !               ! Derivative of the state at time t.                !
        !---------------!---------------------------------------------------!
        subroutine ode_scalar_subroutine ( t, u, uprime )
            real, intent(in) :: t
            real, intent(in) :: u
            real, intent(out) :: uprime
        end subroutine ode_scalar_subroutine
        
        !-------------------------------------------------------------------!
        !   ( INTERFACE ) ( SUBROUTINE ) ode_scalar_subroutine              !
        !-------------------------------------------------------------------!
        !   Defines an standard interface for Cauchy Problem's              !
        !   derivative subroutine.                                          !
        !-------------------------------------------------------------------!
        !   Parameters: !                                                   !
        !---------------!---------------------------------------------------!
        !   IN          ! (real(N)) u                                       !
        !               ! State vector of the cauchy problem at time t.     !
        !---------------!---------------------------------------------------!
        !   IN          ! (real) t                                          !
        !               ! Time at which the state vector are given.         !
        !---------------!---------------------------------------------------!
        !   OUT         ! (real(N)) uprime                                  !
        !               ! Derivative of the state vector at time t.         !
        !---------------!---------------------------------------------------!
        subroutine ode_vectorial_subroutine ( t, m, u, uprime )
            real, intent(in) :: t
            integer, intent(in) :: m
            real, intent(in) :: u(m)
            real, intent(out) :: uprime(m)
        end subroutine ode_vectorial_subroutine

        !-------------------------------------------------------------------!
        !   ( INTERFACE ) ( SUBROUTINE ) temporal_scheme                    !
        !-------------------------------------------------------------------!
        !   Defines an standard interface for Cauchy Problem's              !
        !   temporal schemes.                                               !
        !-------------------------------------------------------------------!
        !   Parameters: !                                                   !
        !---------------!---------------------------------------------------!
        !   IN          ! (ode_function) f(u,t)                             !
        !               ! f: RN x R => RN                                   !
        !               !     u , t => du_dt                                !
        !               ! Derivative function of the Cauchy Problem where   !
        !               ! u is the state vector ant t the evaluation time.  !
        !---------------!---------------------------------------------------!
        !   IN          ! (real) t1                                         !
        !               ! Time at which the initial conditions are given.   !
        !---------------!---------------------------------------------------!
        !   IN          ! (real) t2                                         !
        !               ! Time at which the solution is evaluated.          !
        !---------------!---------------------------------------------------!
        !   IN          ! (real(N)) u1                                      !
        !               ! Initial conditions' state vector at time 't1'.    !
        !---------------!---------------------------------------------------!
        !   OUT         ! (real(N)) u2                                      !
        !               ! Solution's state vector at time 't2'.             !
        !---------------!---------------------------------------------------!
        subroutine temporal_scheme ( f, t1, t2, u1, u2 )
            import :: ode_function
            procedure(ode_function) :: f
            real, intent(in) :: t1, t2, u1(:)
            real, intent(out) :: u2(size(u1))
        end subroutine temporal_scheme
        
    end interface
    
end module ode_interfaces