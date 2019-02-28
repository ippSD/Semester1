

    
!********************************************************************************
!*  Temporal scheme for the solution of the Cauchy problem 
!*
!       U^{n+1} = G( U^n... U^{n-1+p}, F^n... F^{n-1+p}, dt ) 
!*
!*        Inputs: 
!*                F(U) vector valued function of the system of ordinary differential equations 
!*                t1 : initil time 
!*                t2 : final time  
!*                U1 :  vector for the initial condition 
!*
!*        Outputs:
!*                U2   : vector solution for the final state 
!*                ierr : integer variable to inform of internal errors 
!* 
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 2016
!********************************************************************************
module Temporal_Schemes

use Non_Linear_Systems 
use ODE_Interface
implicit none


  
abstract interface 
  
  subroutine Temporal_Scheme(F, t1, t2,  U1, U2, ierr )
                use ODE_Interface
                procedure (ODES) :: F 
                
                
                real, intent(in)    :: t1, t2
                real, intent(in) ::  U1(:)
                
                real, intent(out) ::  U2(:)
                integer, intent(out) :: ierr 
                
 end subroutine Temporal_Scheme
       
end interface  
 
private 

public :: Euler, Runge_Kutta2, Runge_Kutta4, Leap_Frog, Adams_Bashforth, Adams_Bashforth3, Predictor_Corrector1, Inverse_Euler, Crank_Nicolson, Cash_Karp
public :: RK2, RK4 
public :: Temporal_Scheme, ODES

real  :: Tolerance = 1d-10

! ** Runge Kutta Butcher array 
real, save, allocatable :: k(:,:), Up(:) 
real, save, allocatable :: a(:,:), b(:), bs(:), c(:)  



contains 
  

!*******************************************************************************
 ! Explicit Euler  
!*******************************************************************************  
 subroutine Euler(F, t1, t2, U1, U2, ierr )      
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
  
          real :: t, dt 
               
       dt = t2 - t1
       t = t1 
       U2 = U1 + dt * F(U1, t) 
       ierr = 0    
 
 end subroutine 
 
   
!*******************************************************************************
! Runge Kutta 2 stages 
!*******************************************************************************  
subroutine Runge_Kutta2(F, t1, t2, U1, U2, ierr )
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
   
       real :: t, dt  
       real,save :: t_old
       integer :: N 
       real, save, allocatable :: k1(:), k2(:) 

       if (.not.allocated(k1)) then 
                         N = size(U1) 
                         allocate( k1(N), k2(N) )
       elseif (t1 < t_old) then
                         N = size(U1) 
                         deallocate( k1, k2 )
                         allocate( k1(N), k2(N) )
       endif 


       dt = t2-t1;   t = t1 

       k1 = F( U1, t)
       k2 = F( U1 + dt * k1, t + dt )
   
 
       U2 = U1 + dt * ( k1 + k2 )/2


      ierr = 0 
      t_old = t2

    end subroutine 

!*******************************************************************************
! Runge Kutta 4 stages 
!*******************************************************************************  
subroutine Runge_Kutta4(F, t1, t2, U1, U2, ierr )
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
       
       real :: t, dt  
       real,save :: t_old
       integer :: N 
       real, save, allocatable :: k1(:), k2(:), k3(:), k4(:) 

       if (.not.allocated(k1)) then 
                         N = size(U1) 
                         allocate( k1(N), k2(N), k3(N), k4(N) )
       elseif (t1 < t_old) then
                         N = size(U1) 
                         deallocate( k1, k2, k3, k4 )
                         allocate( k1(N), k2(N), k3(N), k4(N) )
       endif 


       dt = t2-t1;   t = t1 

       k1 = F( U1, t)
       k2 = F( U1 + dt * k1/2, t + dt/2 )
       k3 = F( U1 + dt * k2/2, t + dt/2 )
       k4 = F( U1 + dt * k3,   t + dt   )
 
       U2 = U1 + dt * ( k1 + 2*k2 + 2*k3 + k4 )/6


      ierr = 0 
      t_old = t2

    end subroutine 



!*******************************************************************************
!   Leap Frog 
!*******************************************************************************  
subroutine Leap_Frog(F, t1, t2, U1, U2, ierr )   
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
   

       real :: t, dt
       real, save :: t_old = 0;  
       integer ::  N 
       real, save, allocatable ::  U0(:)
       

       dt = t2 - t1;   t = t1
       
       N = size(U1)  
       
       if ( t1 < t_old .or. t_old == 0 ) then 
           
              if ( allocated(U0) ) then 
                               deallocate ( U0 )
              end if 
              allocate(  U0(N)  )
              U0 = U1 
              U2 = U1 + dt * F( U1, t) 
                   
       else 
          
            U2 = U0 + 2 * dt * F( U1, t) 
             
            U0 = U1      
                  
                               
       endif 

       ierr = 0
       t_old = t2
 
    end subroutine


!*******************************************************************************
! Second order Adams Bashforth    
!*******************************************************************************  
subroutine Adams_Bashforth(F, t1, t2, U1, U2, ierr )   
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
   


       real :: t, dt
       real, save :: t_old = 0; 
       integer :: N
       real, save, allocatable ::   F0(:)


       dt = t2-t1;   t = t1 
       N = size(U1) 

       if (.not.allocated(F0)) then 
                         allocate( F0(N)  ) 
                         F0 = F( U1, t)  
                         U2 = U1 + dt * F0
       
       elseif (t1 < t_old ) then
                         deallocate(F0)
                         allocate( F0(N) ) 
                         F0 = F( U1, t)
                         U2 = U1 + dt * F0
                        
       else                                                      
                  U2 = U1 + dt/2 * ( 3*F( U1,t ) - F0 ) 
                  F0 = F( U1, t); 

                               
       endif 

       ierr = 0 
       t_old = t2
 
    end subroutine


!*******************************************************************************
! Third order Adams-Bashforth    
!*******************************************************************************  
subroutine Adams_Bashforth3(F, t1, t2, U1, U2, ierr )   
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 

       real :: t, dt
       real, save :: t_old = 0; 
       integer :: N
       real, save, allocatable ::   F1(:), F0(:), Fm1(:) 
       integer, save :: ipass; 


       dt = t2-t1;   t = t1 
       N = size(U1) 
              
       if (.not.allocated(F1)) then 
                         N = size(U1) 
                         allocate( F1(N), F0(N),Fm1(N) )
                         ipass = - 1 
                         
       elseif (t1 < t_old) then
                         N = size(U1) 
                         deallocate( F1, F0, Fm1 )
                         allocate( F1(N), F0(N),Fm1(N) )
                         ipass = - 1; 
       endif 
       
       ipass = ipass + 1 
       
       
      if (ipass ==0) then
                         F1 = F( U1, t)
                         U2 = U1 + dt * F1
                         F0 = F1 
                         
       elseif (ipass==1) then 
       
                        F1 = F( U1, t)
                        U2 = U1 + dt/2 * ( 3*F1 - F0) 
                        Fm1 = F0
                        F0 = F1 
                    
       else 
                        F1 = F( U1, t)
                                                                         
                        U2 = U1 + dt/12 * ( 23*F1 - 16*F0 + 5*Fm1) 
                        Fm1 = F0
                        F0 = F1 
                               
       endif 

       ierr = 0 
       t_old = t2
 
    end subroutine




!*******************************************************************************
! Predictor- corrector with variable time step 
!*******************************************************************************  
subroutine Predictor_Corrector1(F, t1, t2, U1, U2, ierr )

       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 

       real :: t
       real, save :: t_old
      
       real, save, allocatable ::   Up(:), Uc(:), Kp(:), Kc(:) 
       real :: Residual

       integer, save :: n 
       real, save :: dt=0 


       t = t1 
       N = size(U1) 
       
       if (.not.allocated(Up)) then 
                         
                         allocate(  Up(n), Uc(N), Kp(N), Kc(N)   ) 
       elseif (t1<t_old) then
       
                         deallocate( Up, Uc, Kp, Kc)
                         allocate( Up(N), Uc(N), Kp(N), Kc(N) )
       endif
          
       dt = t2 - t1
       
      
       Up = U1  
       Kp = F( U1, t) 
       Kc = F( U1 + dt*U1, t)
       
       do while( t < t2 ) 
           
              Residual  =  dt * norm2(Kp-Kc)/N 

              if (Residual > Tolerance) then ! this time step is rejected  

                  !  write(*,*) " Rejected dt = ", dt 
                    dt = dt / 2 
                    
              else  
             
                    Up = Up + dt * Kp 
                    Uc = Up + dt * Kc 
                    t = t + dt 
                    
               ! ** new time step 
                    dt = dt * sqrt( Tolerance / Residual ) 
                    
                    Kp  =  F( Up, t ) 
                    Kc  =  F( Uc, t ) 
                  !  write(*,*) " dt = ", dt 
                  
              endif 

              if (t+dt>t2) then 
                                if (t<t2) then 
                                                dt = t2 - t
                                endif 
              endif 

            
       enddo 
    
       U2 = Uc 
       ierr = 0 
       t_old = t2
 
    end subroutine



!*******************************************************************************
! Implicit Inverse Euler    
!*******************************************************************************  
subroutine Inverse_Euler(F, t1, t2, U1, U2, ierr )   
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 

      real, save :: dt
      real, save, allocatable ::  a(:)
     

      ierr = 0 
      if (.not.allocated(a) )  allocate( a(size(U1)) )

      dt = t2-t1
    
      a = U1   
      U2 = U1 
       
    ! Try to find a zero of the residual of the inverse Euler  
      call Newton( F = Residual_IE, x0 = U2 )
              
      
     deallocate( a ) 
  

  contains 
!---------------------------------------------------------------------------
!  Residual of inverse Euler 
!---------------------------------------------------------------------------
function Residual_IE(Z) result(G) 
         real, intent(in) :: Z(:)
         real :: G(size(Z)) 
    
        G = Z - dt *  F( Z, t2) - a 

 end function 
 
end  subroutine 

!*******************************************************************************
! Crank Nicolson scheme    
!*******************************************************************************  
subroutine Crank_Nicolson(F, t1, t2, U1, U2, ierr )   
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 

      real, save :: dt
      real, save, allocatable ::  a(:)
     

      ierr = 0 
      if (.not.allocated(a) )  allocate( a(size(U1)) )

      dt = t2-t1
    
      a = U1  +  dt/2 * F( U1, t1)   
      U2 = U1 

      call Newton( F = Residual_CN, x0 = U2 )
      
      deallocate( a ) 
  

  contains 
!---------------------------------------------------------------------------
!  Residual of the Crank Nicolson scheme 
!---------------------------------------------------------------------------
function Residual_CN(Z) result(G)
         real, intent(in) :: Z(:)
         real :: G(size(Z)) 
       
         G = Z - dt/2 *  F( Z, t2) - a 
 
end function 

end subroutine 
  
  
!*******************************************************************************
!  Runge kutta. Six stages Cash- Karp  
!*******************************************************************************  
   subroutine Cash_Karp(F, t1, t2, U1, U2, ierr )
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
       
       call RKs("Cash_Karp", F, t1, t2, U1, U2, ierr )
	       

   end subroutine
   
!*******************************************************************************
!  Runge kutta 2 stages 
!*******************************************************************************  
   subroutine RK2(F, t1, t2, U1, U2, ierr )
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
       
       call RKs("RK2", F, t1, t2, U1, U2, ierr )
	       

   end subroutine 
   
!*******************************************************************************
!  Runge kutta 4 stages 
!*******************************************************************************  
   subroutine RK4(F, t1, t2, U1, U2, ierr )
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
       
       call RKs("RK4", F, t1, t2, U1, U2, ierr )
	       

   end subroutine 
   
   
  
   
!*******************************************************************************
!  Runge kutta. Explicit of s stages 
!*******************************************************************************  
   subroutine Rks( name, F, t1, t2, U1, U2, ierr )
   
       character(len=*), intent(in) :: name 
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
       
	   real :: t, dt  
       real,save :: t_old
      
       integer, save :: Ns
       integer :: N, i, j, n_subpasos = 1  
     
      
       if (.not.allocated(k)) then 
                         call Butcher_array(name, Ns) 
                         N = size(U1) 
						 allocate( k(Ns, N), Up(N) )
	   elseif (t1 < t_old) then
						 N = size(U1) 
                         deallocate(k, Up)
						 allocate( k(Ns, N), Up(N) )
       endif 
       
       
       dt = (t2-t1)/n_subpasos;   t = t1 
       
       do i=1, Ns 
           
           Up = U1 
           do j=1, i-1
             Up = Up + dt * a(i,j) * k(j, :)
           end do  
           
           k(i, :) =  F( Up, t + c(i) * dt ) 
            
       end do 
       
       U2 = U1
       do i=1, Ns 
           
            U2 = U2 + dt *  b(i) * k(i, :)
            
       end do 
       
	  ierr = 0 
	  t_old = t2
      

   end subroutine
   
!*******************************************************************************
!  Butcher array 
!*******************************************************************************   
subroutine Butcher_array( name, Ns ) 

   character(len=*), intent(in) :: name 
   integer, intent(out) :: Ns 
   
   
   if (name=="Cash_Karp") then 
   
       Ns = 6 
       allocate (a(Ns,Ns-1), b(Ns), bs(Ns), c(Ns)  ) 
       
   
       a(1,:) = [ 0.,          0.,       0.,         0.,            0. ] 
       a(2,:) = [ 1./5,        0.,       0.,         0.,            0. ] 
       a(3,:) = [ 3./40,       9./40,    0.,         0.,            0. ] 
       a(4,:) = [ 3./10,      -9./10,    6./5,       0.,            0. ] 
       a(5,:) = [ -11./54,     5./2,    -70./27,     35./27,        0. ] 
       a(6,:) = [ 1631./55296, 175./512, 575./13824, 44275./110592, 253./4096 ] 
       
       c(:) = [ 0., 1./5, 3./10, 3./5, 1., 7./8 ]
       b(:) = [ 37./378, 0., 250./621, 125./594, 0., 512./1771]
       bs(:) = [2825./27648, 0., 18575./48384, 13525./55296, 277./14336, 1./4 ]   
       
     else if (name=="RK2") then 
   
       Ns = 2 
       allocate (a(Ns,Ns-1), b(Ns), bs(Ns), c(Ns)  ) 
   
       a(1,:) = [ 0.,  0. ] 
       a(2,:) = [ 1.,  0. ] 
       
       c(:) = [ 0., 1. ]
       b(:) = [ 1./2,  1./2 ]
       bs(:) = [1./2, 1./2 ]  
       
     else if (name=="RK4") then 
   
       Ns = 4 
       allocate (a(Ns,Ns-1), b(Ns), bs(Ns), c(Ns)  ) 
       
   
       a(1,:) = [ 0.,          0.,       0. ] 
       a(2,:) = [ 1./2,        0.,       0. ] 
       a(3,:) = [ 0.,        1./2,       0. ] 
       a(4,:) = [ 0.,          0.,       1. ] 
      
       
       c(:) =  [ 0., 1./2, 1./2, 1. ]
       b(:)  = [ 1./6, 1./3, 1./3, 1./6 ]
       bs(:) = [ 1./6, 1./3, 1./3, 1./6 ] 
       
       
   else 
       
       write(*,*) " Error Butcher array " 
       stop 
       
   end if 
   
 end subroutine   
   

end module 





 
 
  
