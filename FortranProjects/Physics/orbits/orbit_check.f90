module orbit_check
    implicit none
    contains
    
    !-----------------------------------------------------------------------!
    !   ( FUNCTION ) check_state_matrix                                     !
    !-----------------------------------------------------------------------!
    !   Checks if the state matrix (state vector across the time domain)    !
    !   is compatible with the orbit functions of this repo.                !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(0:M,N)) u                                       !
    !               ! State matrix to be checked.                           !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (logical) f                                           !
    !               ! True if matrix is valid, else false.                  !
    !---------------!-------------------------------------------------------!
    logical function check_state_matrix(u) result(f)
        real, intent(in) :: u(0:,:)  ! u(0:m,n)
        
        f = check_state_vector( u(1,:) )
    end function check_state_matrix
    
    !-----------------------------------------------------------------------!
    !   ( FUNCTION ) check_state_vector                                     !
    !-----------------------------------------------------------------------!
    !   Checks if the state vector is compatible with the orbit functions   !
    !   of this repo.                                                       !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) u                                           !
    !               ! State vector to be checked.                           !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (logical) f                                           !
    !               ! True if vector is valid, else false.                  !
    !---------------!-------------------------------------------------------!
    logical function check_state_vector(u) result(f)
        real, intent(in) :: u(:)  ! u(m,n)
        integer :: n
        
        n = size(u(:))
        f = .true.
        
        if (mod(n, 7) /= 0) then
            write(*,*) "Error: invalid state vector."
            write(*,*) "State vector should be size '7N' with"
            write(*,*) "'N' being the number of bodies"
            f = .false.
        end if
    end function check_state_vector
    
    !-----------------------------------------------------------------------!
    !   ( FUNCTION ) check_existing_body                                    !
    !-----------------------------------------------------------------------!
    !   Checks if the selected bodystate vector is compatible with          !
    !   the orbit functions of this repo.                                   !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) u                                           !
    !               ! State vector to be checked.                           !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (logical) f                                           !
    !               ! True if vector is valid, else false.                  !
    !---------------!-------------------------------------------------------!
    logical function check_existing_body(u, l) result(f)
        real, intent(in) :: u(:)
        integer, intent(in) :: l
        integer :: n
        
        f = .false.
        
        if ( check_state_vector(u) ) then
            n = size(u) / 7
            if ( 0 < l .and. l <= n ) then
                f = .true.
            else
                write(*,*) "Error: invalid body number."
                write(*,*) "Maximum body index: ", n
                write(*,*) "Current: ", l
            end if
        end if
    end function check_existing_body


end module orbit_check