module Non_Linear_Systems 

 use Jacobian_module 
 use Linear_systems
  
  implicit none 


contains 

!***************************************************************************
!*  Newton solver 
!*                x0   : initial guess and output value 
!*                F(X) : vector function 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!***************************************************************************
subroutine Newton(F, x0) 
 
 procedure (FunctionRN_RN) :: F 
 real, intent(inout) :: x0(:) 


   real ::   Dx( size(x0) ), b(size(x0)), eps
   real :: J( size(x0), size(x0) )
   integer :: iteration 

   
   integer :: N 

   N = size(x0) 

   Dx = 2 * x0 
   iteration = 0 
   eps = 1 
 
   do while ( eps > 1d-8 .and. iteration < 1000 )
    
      iteration = iteration + 1 
      J = Jacobian( F, x0 ) 
     
      call LU_factorization( J ) 
      b = F(x0);
      Dx = Solve_LU( J,  b ) 
      
      x0 = x0 - Dx;  

      eps = norm2( DX )  
    
   end do 
   
   if (iteration>900) then 
      write(*,*) " morm2(J) =", maxval(J),  minval(J) 
      write(*,*) " Norm2(Dx) = ", eps, iteration  
   endif 
 
   
end subroutine
!***************************************************************************
!*  Newton solver 
!*                x0   : initial guess and output value 
!*                F(X) : vector function 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!***************************************************************************
subroutine Newton_Debug(F, x0) 

 procedure (FunctionRN_RN) :: F 
 real, intent(inout) :: x0(:) 


   real ::   Dx( size(x0) ), b(size(x0)), eps
   real :: J( size(x0), size(x0) )
   integer :: iteration 

   integer :: N 
   

   N = size(x0); 
  
   

   Dx = 2 * x0; 
   iteration = 0; 
   eps = 1; 
  
   write(*,*) " Inside Newton " 
   do while ( eps > 1d-8 .and. iteration < 1000 )
    
      write(*,*) " iteration =",  iteration 


      iteration = iteration + 1 
      J = Jacobian( F, x0 )
      
   !  write(*,*) " morm2(Jacobian) =", norm2(J)
   !   do i=1, size(x0)
   !       write(*,*) " Jacobian(i,:)=", i, maxloc(J(i,:)) 
   !   end do 
          
      
      call LU_factorization( J ) 
      b = F(x0);
      
  !    write(*,*) " morm2(b) =", norm2(b) 
  !    write(*,*) " morm2(LU) =", norm2(J)
      
      Dx = Solve_LU( J,  b ) 
    
      x0 = x0 - Dx;  

      eps = norm2( DX );  
    
     write(*,*) " morm2(J) =", maxval(J),  minval(J) 
     write(*,*) " Norm2(Dx) = ", eps, iteration  
  !   read(*,*) 
    
     

   end do 
   
   if (iteration>900) then 
      write(*,*) " morm2(J) =", maxval(J),  minval(J) 
      write(*,*) " Norm2(Dx) = ", eps, iteration  
   endif 
  
 
 
end subroutine



end module 
