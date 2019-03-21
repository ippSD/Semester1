
!***********************************************************************
!   It integrates the Cauchy problem.   
!***********************************************************************    
module Cauchy_Problem



use Numerical_recipes
use External_Temporal_Schemes
implicit none   

 private 
 public :: Cauchy_ProblemS, Floquet_multipliers, System_matrix
 
 public :: Temporal_scheme, ODES
 public :: Euler, Runge_Kutta2, Runge_Kutta4, Leap_Frog, Adams_Bashforth, Adams_Bashforth3, Predictor_Corrector1, Inverse_Euler, Crank_Nicolson, Cash_Karp
 public :: RK2, RK4 


abstract interface 

 function F_R_RN(N, t) 
         integer, intent(in) :: N 
         real, intent(in) ::  t 
         real :: F_R_RN(N) 
      
 end function 

end interface 
    

contains   

!************************************************************************************************
! It integrates the following Cauchy problem 
!
!       dU/dt = F(U, t),      U(0) = U^0 (CI)   
!
!     Inputs: 
!            Time_Domain(:)        : time discretization 
!            Differential_operator : vector values function F(U, t)
!            Scheme                : Optional Temporal numerical scheme ( default: Runge_Kutta4)
!            Solution(0,:)         : Initial condition 
!
!     Outputs: 
!            Solution(:,:)         : first index represents time 
!                                    second index  represents the i-component of solution
!     
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!***************************************************************************************************
subroutine Cauchy_ProblemS( Time_Domain, Differential_operator, Scheme, & 
                            Solution ) 
     real, intent(in) :: Time_Domain(:) 
     procedure (ODES) :: Differential_operator
     procedure (Temporal_Scheme), optional :: Scheme
     real, intent(out) :: Solution(:,:) 
     
!  *** Initial and final time
       real :: start, finish         ! CPU time 
       real ::  t1, t2               ! simulation time step  
       integer ::  i, N_steps, ierr
       
       call cpu_time(start)
       N_steps = size(Time_Domain) - 1; 
     
!  *** loop for temporal integration 
       do i=1, N_steps 
           t1 = Time_Domain(i) 
           t2 = Time_Domain(i+1)
         
           if (present(Scheme))  then 
               
                  call Scheme(        Differential_operator, t1, t2,        & 
                                      Solution(i,:), Solution(i+1,:), ierr )  
           else 
                  call Runge_Kutta4(  Differential_operator, t1, t2,        & 
                                      Solution(i,:), Solution(i+1,:), ierr )  
           endif 
           if (ierr>0) exit 
       enddo 
       
 call cpu_time(finish)
 write(*, '("Cauchy_Problem, CPU Time = ",f6.3," seconds.")') finish - start
 
end subroutine 

 
  
!***********************************************************************
!  If determines the Floquet multipliers of the system dU/dt = F(U,t)  
!                            
!     Inputs: 
!             Differential_operator(U,t) :  function of Cauchy problem 
!             U0(t)                      :  function or periodic solution 
!             Period 
!                            
!     Outputs: 
!             rho(:)   complex Floquet multipliers 
!             Phi(:,:) fundamental matrix solution 
!                            
!***********************************************************************
subroutine Floquet_multipliers( Differential_operator, U0, Period, rho, Phi) 
            procedure (ODES) :: Differential_operator 
            procedure (F_R_RN) :: U0 
            real, intent(in) ::  Period
            real, intent(out) :: Phi(:,:) 
            complex, intent(out) :: rho(:) 
          
          
          integer, parameter :: M = 1000   
          real :: Time(0:M)
          real, allocatable :: U(:, :), P(:,:)   
          integer :: i, N  
          
          N = size( rho )
          allocate( U(0:M, N), P(N,N) ) 
          
          Time = [ ( Period*i/M, i=0, M) ] 
          
          do i=1, N 
              
              U(0, :) = 0 
              U(0, i) = 1
              call Cauchy_ProblemS( Time, Linear_system, Cash_Karp, U ) 
              
              Phi(:,i) = U(M, :) 
              
              !call qplot( U(:,1), U(:, 2), M+1) 
              
          end do 
          
          P = Phi
          call Eigenvalues_QR( P, rho )
          
          
contains 
!-----------------------------------------------------------------
! Linearised system around the periodic solution of function U0(t) 
!-----------------------------------------------------------------
function Linear_system( U, t ) result(F) 
                real :: U(:), t 
                real :: F( size(U) ) 
                
              real :: A(N, N), Up(N) 
              
      Up = U0(N, t)   
      
      call System_matrix( Differential_operator, Up, t, A) 
      
      F = matmul( A, U ) 
      
     

end function 
             
end subroutine 


!***********************************************************************
!*  If gives the linearized operator F(U,t) in some point U0 
!***********************************************************************
subroutine System_matrix( F, U0, t, A) 
             procedure (ODES) :: F 
             real, intent(in) :: U0(:), t 
             real, intent(out) :: A( size(U0), size(U0) ) 
          
          
          real :: F0( size(U0) ), delta( size(U0) )
          real ::  eps = 1d-6  
          integer :: i, N 
           
          
          F0 = F( U0, t ) 
          
          N = size(U0) 
          
          do i=1, N 
              
              delta = 0 
              delta(i) = eps 
              A(:, i) = ( F( U0 + delta, t ) - F0 )/eps
              
          end do 
          
             
end subroutine 
                          
                            
                            
                            
                            

end module 
