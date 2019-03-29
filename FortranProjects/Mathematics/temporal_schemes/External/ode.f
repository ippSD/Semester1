include "d1mach.f90"
      subroutine ode(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
      implicit real(a-h,o-z)

!   double precision subroutine ode integrates a system of neqn
!   first order ordinary differential equations of the form:
!             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
!             y(i) given at  t .
!   the subroutine integrates from  t  to  tout .  on return the
!   parameters in the call list are set for continuing the integration.
!   the user has only to define a new value  tout  and call  ode  again.
!
!   the differential equations are actually solved by a suite of codes
!   de ,  step , and  intrp .  ode  allocates virtual storage in the
!   arrays  work  and  iwork  and calls  de .  de  is a supervisor which
!   directs the solution.  it calls on the routines  step  and  intrp
!   to advance the integration and to interpolate at output points.
!   step  uses a modified divided difference form of the adams pece
!   formulas and local extrapolation.  it adjusts the order and step
!   size to control the local error per unit step in a generalized
!   sense.  normally each call to  step  advances the solution one step
!   in the direction of  tout .  for reasons of efficiency  de
!   integrates beyond  tout  internally, though never beyond
!   t+10*(tout-t), and calls  intrp  to interpolate the solution at
!   tout .  an option is provided to stop the integration at  tout  but
!   it should be used only if it is impossible to continue the
!   integration beyond  tout .
!
!   this code is completely explained and documented in the text,
!   computer solution of ordinary differential equations:  the initial
!   value problem  by l. f. shampine and m. k. gordon.
!
!   the parameters represent:
!      f -- double precision subroutine f(t,y,yp) to evaluate
!                derivatives yp(i)=dy(i)/dt
!      neqn -- number of equations to be integrated (integer*4)
!      y(*) -- solution vector at t                 (real*8)
!      t -- independent variable                    (real*8)
!      tout -- point at which solution is desired   (real*8)
!      relerr,abserr -- relative and absolute error tolerances for local
!           error test (real*8).  at each step the code requires
!             dabs(local error) .le. dabs(y)*relerr + abserr
!           for each component of the local error and solution vectors
!      iflag -- indicates status of integration     (integer*4)
!      work(*)  (real*8)  -- arrays to hold information internal to
!      iwork(*) (integer*4)    which is necessary for subsequent calls
!
!   first call to ode --
!
!   the user must provide storage in his calling program for the arrays
!   in the call list,
!      y(neqn), work(100+21*neqn), iwork(5),
!   declare  f  in an external statement, supply the double precision
!   subroutine f(t,y,yp)  to evaluate
!      dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn))
!   and initialize the parameters:
!      neqn -- number of equations to be integrated
!      y(*) -- vector of initial conditions
!      t -- starting point of integration
!      tout -- point at which solution is desired
!      relerr,abserr -- relative and absolute local error tolerances
!      iflag -- +1,-1.  indicator to initialize the code.  normal input
!           is +1.  the user should set iflag=-1 only if it is
!           impossible to continue the integration beyond  tout .
!   all parameters except  f ,  neqn  and  tout  may be altered by the
!   code on output so must be variables in the calling program.
!
!   output from  ode  --
!
!      neqn -- unchanged
!      y(*) -- solution at  t
!      t -- last point reached in integration.  normal return has
!           t = tout .
!      tout -- unchanged
!      relerr,abserr -- normal return has tolerances unchanged.  iflag=3
!           signals tolerances increased
!      iflag = 2 -- normal return.  integration reached  tout
!            = 3 -- integration did not reach  tout  because error
!                   tolerances too small.  relerr ,  abserr  increased
!                   appropriately for continuing
!            = 4 -- integration did not reach  tout  because more than
!                   500 steps needed
!            = 5 -- integration did not reach  tout  because equations
!                   appear to be stiff
!            = 6 -- invalid input parameters (fatal error)
!           the value of  iflag  is returned negative when the input
!           value is negative and the integration does not reach  tout ,
!           i.e., -3, -4, -5.
!      work(*),iwork(*) -- information generally of no interest to the
!           user but necessary for subsequent calls.
!
!   subsequent calls to  ode --
!
!   subroutine  ode  returns with all information needed to continue
!   the integration.  if the integration reached  tout , the user need
!   only define a new  tout  and call again.  if the integration did not
!   reach  tout  and the user wants to continue, he just calls again.
!   the output value of  iflag  is the appropriate input value for
!   subsequent calls.  the only situation in which it should be altered
!   is to stop the integration internally at the new  tout , i.e.,
!   change output  iflag=2  to input  iflag=-2 .  error tolerances may
!   be changed by the user before continuing.  all other parameters must
!   remain unchanged.
!
!***********************************************************************
!*  subroutines  de  and  step  contain machine dependent constants. *
!*  be sure they are set before using  ode .                          *
!***********************************************************************
!
      logical start,phase1,nornd
      dimension y(neqn),work(1),iwork(5)
      external f
      data ialpha,ibeta,isig,iv,iw,ig,iphase,ipsi,ix,ih,ihold,istart, &
          itold,idelsn/1,13,25,38,50,62,75,76,88,89,90,91,92,93/
      iyy = 100
      iwt = iyy + neqn
      ip = iwt + neqn
      iyp = ip + neqn
      iypout = iyp + neqn
      iphi = iypout + neqn
      if(iabs(iflag) .eq. 1) go to 1
      start = work(istart) .gt. 0.0d0
      phase1 = work(iphase) .gt. 0.0d0
      nornd = iwork(2) .ne. -1
    1 call de(f,neqn,y,t,tout,relerr,abserr,iflag,work(iyy), &
      work(iwt),work(ip),work(iyp),work(iypout),work(iphi),&
      work(ialpha),work(ibeta),work(isig),work(iv),work(iw),work(ig),&
      phase1,work(ipsi),work(ix),work(ih),work(ihold),start,&
      work(itold),work(idelsn),iwork(1),nornd,iwork(3),iwork(4),iwork(5))
      work(istart) = -1.0d0
      if(start) work(istart) = 1.0d0
      work(iphase) = -1.0d0
      if(phase1) work(iphase) = 1.0d0
      iwork(2) = -1
      if(nornd) iwork(2) = 1
      return
      end
      subroutine de(f,neqn,y,t,tout,relerr,abserr,iflag, &
          yy,wt,p,yp,ypout,phi,alpha,beta,sig,v,w,g,phase1,psi,x,h,hold, &
              start,told,delsgn,ns,nornd,k,kold,isnold)
      implicit real(a-h,o-z)
!
!   ode  merely allocates storage for  de  to relieve the user of the
!   inconvenience of a long call list.  consequently  de  is used as
!   described in the comments for  ode .
!
!   this code is completely explained and documented in the text,
!   computer solution of ordinary differential equations:  the initial
!   value problem  by l. f. shampine and m. k. gordon.
!
      logical stiff,crash,start,phase1,nornd
        dimension y(neqn),yy(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn), &
              ypout(neqn),psi(12),alpha(12),beta(12),sig(13),v(12),w(12), g(13)
      external f
!
!***********************************************************************
!*  the only machine dependent constant is based on the machine unit   *
!*  roundoff error  u  which is the smallest positive number such that *
!*  1.0+u .gt. 1.0 .  u  must be calculated and  fouru=4.0*u  inserted *
!*  in the following data statement before using  de .  the routine    *
!*  machin  calculates  u .  fouru  and  twou=2.0*u  must also be      *
!*  inserted in subroutine  step  before calling  de .                 *
!     data fouru/.888d-15/                                              ***
!***********************************************************************
!
!   the constant  maxnum  is the maximum number of steps allowed in one
!   call to  de .  the user may change this limit by altering the
!   following statement
      data maxnum/500/
!
!            ***            ***            ***
!   test for improper parameters
!
      fouru = 4.0 * d1mach(4)
      if(neqn .lt. 1) go to 10
      if(t .eq. tout) go to 10
      if(relerr .lt. 0.0d0  .or.  abserr .lt. 0.0d0) go to 10
      eps = dmax1(relerr,abserr)
      if(eps .le. 0.0d0) go to 10
      if(iflag .eq. 0) go to 10
      isn = isign(1,iflag)
      iflag = iabs(iflag)
      if(iflag .eq. 1) go to 20
      if(t .ne. told) go to 10
      if(iflag .ge. 2  .and.  iflag .le. 5) go to 20
   10 iflag = 6
      return
!
!   on each call set interval of integration and counter for number of
!   steps.  adjust input error tolerances to define weight vector for
!   subroutine  step
!
   20 del = tout - t
      absdel = dabs(del)
      tend = t + 10.0d0*del
      if(isn .lt. 0) tend = tout
      nostep = 0
      kle4 = 0
      stiff = .false.
      releps = relerr/eps
      abseps = abserr/eps
      if(iflag .eq. 1) go to 30
      if(isnold .lt. 0) go to 30
      if(delsgn*del .gt. 0.0d0) go to 50
!
!   on start and restart also set work variables x and yy(*), store the
!   direction of integration and initialize the step size
!
   30 start = .true.
      x = t
      do 40 l = 1,neqn
   40   yy(l) = y(l)
      delsgn = dsign(1.0d0,del)
      h = dsign(dmax1(dabs(tout-x),fouru*dabs(x)),tout-x)
!
!   if already past output point, interpolate and return
!
   50 if(dabs(x-t) .lt. absdel) go to 60
      call intrp(x,yy,tout,y,ypout,neqn,kold,phi,psi)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
!
!   if cannot go past output point and sufficiently close,
!   extrapolate and return
!
   60 if(isn .gt. 0  .or.  dabs(tout-x) .ge. fouru*dabs(x)) go to 80
      h = tout - x
      call f(x,yy,yp)
      do 70 l = 1,neqn
   70   y(l) = yy(l) + h*yp(l)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
!
!   test for too many steps
!
   80 if(nostep .lt. maxnum) go to 100
      iflag = isn*4
      if(stiff) iflag = isn*5
      do 90 l = 1,neqn
   90   y(l) = yy(l)
      t = x
      told = t
      isnold = 1
      return
!
!   limit step size, set weight vector and take a step
!
  100 h = dsign(dmin1(dabs(h),dabs(tend-x)),h)
      do 110 l = 1,neqn
  110 wt(l) = releps*dabs(yy(l)) + abseps
      call step(x,yy,f,neqn,h,eps,wt,start,&
          hold,k,kold,crash,phi,p,yp,psi,&
              alpha,beta,sig,v,w,g,phase1,ns,nornd)
!
!   test for tolerances too small
!
      if(.not.crash) go to 130
      iflag = isn*3
      relerr = eps*releps
      abserr = eps*abseps
      do 120 l = 1,neqn
  120   y(l) = yy(l)
      t = x
      told = t
      isnold = 1
      return
!
!   augment counter on number of steps and test for stiffness
!
  130 nostep = nostep + 1
      kle4 = kle4 + 1
      if(kold .gt. 4) kle4 = 0
      if(kle4 .ge. 50) stiff = .true.
      go to 50
      end
      subroutine step(x,y,f,neqn,h,eps,wt,start,&
          hold,k,kold,crash,phi,p,yp,psi,&
          alpha,beta,sig,v,w,g,phase1,ns,nornd)
      implicit real(a-h,o-z)
!
!   double precision subroutine  step
!   integrates a system of first order ordinary
!   differential equations one step, normally from x to x+h, using a
!   modified divided difference form of the adams pece formulas.  local
!   extrapolation is used to improve absolute stability and accuracy.
!   the code adjusts its order and step size to control the local error
!   per unit step in a generalized sense.  special devices are included
!   to control roundoff error and to detect when the user is requesting
!   too much accuracy.
!
!   this code is completely explained and documented in the text,
!   computer solution of ordinary differential equations:  the initial
!   value problem  by l. f. shampine and m. k. gordon.
!
!
!   the parameters represent:
!      x -- independent variable             (real*8)
!      y(*) -- solution vector at x          (real*8)
!      yp(*) -- derivative of solution vector at  x  after successful
!           step                             (real*8)
!      neqn -- number of equations to be integrated (integer*4)
!      h -- appropriate step size for next step.  normally determined by
!           code                             (real*8)
!      eps -- local error tolerance.  must be variable  (real*8)
!      wt(*) -- vector of weights for error criterion   (real*8)
!      start -- logical variable set .true. for first step,  .false.
!           otherwise                        (logical*4)
!      hold -- step size used for last successful step  (real*8)
!      k -- appropriate order for next step (determined by code)
!      kold -- order used for last successful step
!      crash -- logical variable set .true. when no step can be taken,
!           .false. otherwise.
!   the arrays  phi, psi  are required for the interpolation subroutine
!   intrp.  the array p is internal to the code.  all are real*8
!
!   input to  step
!
!      first call --
!
!   the user must provide storage in his driver program for all arrays
!   in the call list, namely
!
!     dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
!
!   the user must also declare  start  and  crash  logical variables
!   and  f  an external subroutine, supply the subroutine  f(x,y,yp)
!   to evaluate
!      dy(i)/dx = yp(i) = f(x,y(1),y(2),...,y(neqn))
!   and initialize only the following parameters:
!      x -- initial value of the independent variable
!      y(*) -- vector of initial values of dependent variables
!      neqn -- number of equations to be integrated
!      h -- nominal step size indicating direction of integration
!           and maximum size of step.  must be variable
!      eps -- local error tolerance per step.  must be variable
!      wt(*) -- vector of non-zero weights for error criterion
!      start -- .true.
!
!   step  requires the l2 norm of the vector with components
!   local error(l)/wt(l)  be less than  eps  for a successful step.  the
!   array  wt  allows the user to specify an error test appropriate
!   for his problem.  for example,
!      wt(l) = 1.0  specifies absolute error,
!            = dabs(y(l))  error relative to the most recent value of
!                 the l-th component of the solution,
!            = dabs(yp(l))  error relative to the most recent value of
!                 the l-th component of the derivative,
!            = dmax1(wt(l),dabs(y(l)))  error relative to the largest
!                 magnitude of l-th component obtained so far,
!            = dabs(y(l))*relerr/eps + abserr/eps  specifies a mixed
!                 relative-absolute test where  relerr  is relative
!                 error,  abserr  is absolute error and  eps =
!                 dmax1(relerr,abserr) .
!
!      subsequent calls --
!
!   subroutine  step  is designed so that all information needed to
!   continue the integration, including the step size  h  and the order
!   k , is returned with each step.  with the exception of the step
!   size, the error tolerance, and the weights, none of the parameters
!   should be altered.  the array  wt  must be updated after each step
!   to maintain relative error tests like those above.  normally the
!   integration is continued just beyond the desired endpoint and the
!   solution interpolated there with subroutine  intrp .  if it is
!   impossible to integrate beyond the endpoint, the step size may be
!   reduced to hit the endpoint since the code will not take a step
!   larger than the  h  input.  changing the direction of integration,
!   i.e., the sign of  h , requires the user set  start = .true. before
!   calling  step  again.  this is the only situation in which  start
!   should be altered.
!
!   output from  step
!
!      successful step --
!
!   the subroutine returns after each successful step with  start  and
!   crash  set .false. .  x  represents the independent variable
!   advanced one step of length  hold  from its value on input and  y
!   the solution vector at the new value of  x .  all other parameters
!   represent information corresponding to the new  x  needed to
!   continue the integration.
!
!      unsuccessful step --
!
!   when the error tolerance is too small for the machine precision,
!   the subroutine returns without taking a step and  crash = .true. .
!   an appropriate step size and error tolerance for continuing are
!   estimated and all other information is restored as upon input
!   before returning.  to continue with the larger tolerance, the user
!   just calls the code again.  a restart is neither required nor
!   desirable.
!
      logical start,crash,phase1,nornd
      dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
      dimension alpha(12),beta(12),sig(13),w(12),v(12),g(13),&
      gstr(13),two(13)
      external f
!***********************************************************************
!*  the only machine dependent constants are based on the machine unit *
!*  roundoff error  u  which is the smallest positive number such that *
!*  1.0+u .gt. 1.0  .  the user must calculate  u  and insert          *
!*  twou=2.0*u  and  fouru=4.0*u  in the data statement before calling *
!*  the code.  the routine  machin  calculates  u .                    *
!     data twou,fouru/.444d-15,.888d-15/                                ***
!***********************************************************************
      data two/2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,256.0d0,&
512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/
      data gstr/0.500d0,0.0833d0,0.0417d0,0.0264d0,0.0188d0,0.0143d0,&
 0.0114d0,0.00936d0,0.00789d0,0.00679d0,0.00592d0,0.00524d0,&
0.00468d0/
!     data g(1),g(2)/1.0,0.5/,sig(1)/1.0/
!
!
      twou = 2.0 * d1mach(4)
      fouru = 2.0 * twou
!       ***     begin block 0     ***
!   check if step size or error tolerance is too small for machine
!   precision.  if first step, initialize phi array and estimate a
!   starting step size.
!                   ***
!
!   if step size is too small, determine an acceptable one
!
      crash = .true.
      if(dabs(h) .ge. fouru*dabs(x)) go to 5
      h = dsign(fouru*dabs(x),h)
      return
    5 p5eps = 0.5d0*eps
!
!   if error tolerance is too small, increase it to an acceptable value
!
      round = 0.0d0
      do 10 l = 1,neqn
   10   round = round + (y(l)/wt(l))**2
      round = twou*dsqrt(round)
      if(p5eps .ge. round) go to 15
      eps = 2.0*round*(1.0d0 + fouru)
      return
   15 crash = .false.
      g(1)=1.0d0
      g(2)=0.5d0
      sig(1)=1.0d0
      if(.not.start) go to 99
!
!   initialize.  compute appropriate step size for first step
!
      call f(x,y,yp)
      sum = 0.0d0
      do 20 l = 1,neqn
        phi(l,1) = yp(l)
        phi(l,2) = 0.0d0
   20   sum = sum + (yp(l)/wt(l))**2
      sum = dsqrt(sum)
      absh = dabs(h)
      if(eps .lt. 16.0d0*sum*h*h) absh = 0.25d0*dsqrt(eps/sum)
      h = dsign(dmax1(absh,fouru*dabs(x)),h)
      hold = 0.0d0
      k = 1
      kold = 0
      start = .false.
      phase1 = .true.
      nornd = .true.
      if(p5eps .gt. 100.0d0*round) go to 99
      nornd = .false.
      do 25 l = 1,neqn
   25   phi(l,15) = 0.0d0
   99 ifail = 0
!       ***     end block 0     ***
!
!       ***     begin block 1     ***
!   compute coefficients of formulas for this step.  avoid computing
!   those quantities not changed when step size is not changed.
!                   ***
!
  100 kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
!
!   ns is the number of steps taken with size h, including the current
!   one.  when k.lt.ns, no coefficients change
!
      if(h .ne. hold) ns = 0
      if(ns.le.kold)   ns=ns+1
      nsp1 = ns+1
      if (k .lt. ns) go to 199
!
!   compute those components of alpha(*),beta(*),psi(*),sig(*) which
!   are changed
!
      beta(ns) = 1.0d0
      realns = ns
      alpha(ns) = 1.0d0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0d0
      if(k .lt. nsp1) go to 110
      do 105 i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
  105   sig(i+1) = reali*alpha(i)*sig(i)
  110 psi(k) = temp1
!
!   compute coefficients g(*)
!
!   initialize v(*) and set w(*).  g(2) is set in data statement
!
      if(ns .gt. 1) go to 120
      do 115 iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0d0/temp3
  115   w(iq) = v(iq)
      go to 140
!
!   if order was raised, update diagonal part of v(*)
!
  120 if(k .le. kold) go to 130
      temp4 = k*kp1
      v(k) = 1.0d0/temp4
      nsm2 = ns-2
      if(nsm2 .lt. 1) go to 130
      do 125 j = 1,nsm2
        i = k-j
  125   v(i) = v(i) - alpha(j+1)*v(i+1)
!
!   update v(*) and set w(*)
!
  130 limit1 = kp1 - ns
      temp5 = alpha(ns)
      do 135 iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
  135   w(iq) = v(iq)
      g(nsp1) = w(1)
!
!   compute the g(*) in the work vector w(*)
!
  140 nsp2 = ns + 2
      if(kp1 .lt. nsp2) go to 199
      do 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        do 145 iq = 1,limit2
  145     w(iq) = w(iq) - temp6*w(iq+1)
  150   g(i) = w(1)
  199   continue
!       ***     end block 1     ***
!
!       ***     begin block 2     ***
!   predict a solution p(*), evaluate derivatives using predicted
!   solution, estimate local error at order k and errors at orders k,
!   k-1, k-2 as if constant step size were used.
!                   ***
!
!   change phi to phi star
!
      if(k .lt. nsp1) go to 215
      do 210 i = nsp1,k
        temp1 = beta(i)
        do 205 l = 1,neqn
  205     phi(l,i) = temp1*phi(l,i)
  210   continue
!
!   predict solution and differences
!
  215 do 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0d0
  220   p(l) = 0.0d0
      do 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        do 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
  225     phi(l,i) = phi(l,i) + phi(l,ip1)
  230   continue
      if(nornd) go to 240
      do 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
  235   phi(l,16) = (p(l) - y(l)) - tau
      go to 250
  240 do 245 l = 1,neqn
  245   p(l) = y(l) + h*p(l)
  250 xold = x
      x = x + h
      absh = dabs(h)
      call f(x,p,yp)
!
!   estimate errors at orders k,k-1,k-2
!
      erkm2 = 0.0d0
      erkm1 = 0.0d0
      erk = 0.0d0
      do 265 l = 1,neqn
        temp3 = 1.0d0/wt(l)
        temp4 = yp(l) - phi(l,1)
        if(km2)265,260,255
  255   erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
  260   erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
  265   erk = erk + (temp4*temp3)**2
      if(km2)280,275,270
  270 erkm2 = absh*sig(km1)*gstr(km2)*dsqrt(erkm2)
  275 erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
  280 temp5 = absh*dsqrt(erk)
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      knew = k
!
!   test if order should be lowered
!
      if(km2)299,290,285
  285 if(dmax1(erkm1,erkm2) .le. erk) knew = km1
      go to 299
  290 if(erkm1 .le. 0.5d0*erk) knew = km1
!
!   test if step successful
!
  299 if(err .le. eps) go to 400
!       ***     end block 2     ***
!
!       ***     begin block 3     ***
!   the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
!   if third consecutive failure, set order to one.  if step fails more
!   than three times, consider an optimal step size.  double error
!   tolerance and return if estimated step size is too small for machine
!   precision.
!                   ***
!
!   restore x, phi(*,*) and psi(*)
!
      phase1 = .false.
      x = xold
      do 310 i = 1,k
        temp1 = 1.0d0/beta(i)
        ip1 = i+1
        do 305 l = 1,neqn
  305     phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
  310   continue
      if(k .lt. 2) go to 320
      do 315 i = 2,k
  315   psi(i-1) = psi(i) - h
!
!   on third failure, set order to one.  thereafter, use optimal step
!   size
!
  320 ifail = ifail + 1
      temp2 = 0.5d0
      if(ifail - 3) 335,330,325
  325 if(p5eps .lt. 0.25d0*erk) temp2 = dsqrt(p5eps/erk)
  330 knew = 1
  335 h = temp2*h
      k = knew
      if(dabs(h) .ge. fouru*dabs(x)) go to 340
      crash = .true.
      h = dsign(fouru*dabs(x),h)
      eps = eps + eps
      return
  340 go to 100
!       ***     end block 3     ***
!
!       ***     begin block 4     ***
!   the step is successful.  correct the predicted solution, evaluate
!   the derivatives using the corrected solution and update the
!   differences.  determine best order and step size for next step.
!                   ***
  400 kold = k
      hold = h
!
!   correct and evaluate
!
      temp1 = h*g(kp1)
      if(nornd) go to 410
      do 405 l = 1,neqn
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
  405   phi(l,15) = (y(l) - p(l)) - rho
      go to 420
  410 do 415 l = 1,neqn
  415   y(l) = p(l) + temp1*(yp(l) - phi(l,1))
  420 call f(x,y,yp)
!
!   update differences for next step
!
      do 425 l = 1,neqn
        phi(l,kp1) = yp(l) - phi(l,1)
  425   phi(l,kp2) = phi(l,kp1) - phi(l,kp2)
      do 435 i = 1,k
        do 430 l = 1,neqn
  430     phi(l,i) = phi(l,i) + phi(l,kp1)
  435   continue
!
!   estimate error at order k+1 unless:
!     in first phase when always raise order,
!     already decided to lower order,
!     step size not constant so estimate unreliable
!
      erkp1 = 0.0d0
      if(knew .eq. km1  .or.  k .eq. 12) phase1 = .false.
      if(phase1) go to 450
      if(knew .eq. km1) go to 455
      if(kp1 .gt. ns) go to 460
      do 440 l = 1,neqn
  440   erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
      erkp1 = absh*gstr(kp1)*dsqrt(erkp1)
!
!   using estimated error at order k+1, determine appropriate order
!   for next step
!
      if(k .gt. 1) go to 445
      if(erkp1 .ge. 0.5d0*erk) go to 460
      go to 450
  445 if(erkm1 .le. dmin1(erk,erkp1)) go to 455
      if(erkp1 .ge. erk  .or.  k .eq. 12) go to 460
!
!   here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would have
!   been lowered in block 2.  thus order is to be raised
!
!   raise order
!
  450 k = kp1
      erk = erkp1
      go to 460
!
!   lower order
!
  455 k = km1
      erk = erkm1
!
!   with new order determine appropriate step size for next step
!
  460 hnew = h + h
      if(phase1) go to 465
      if(p5eps .ge. erk*two(k+1)) go to 465
      hnew = h
      if(p5eps .ge. erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0d0/temp2)
      hnew = absh*dmax1(0.5d0,dmin1(0.9d0,r))
      hnew = dsign(dmax1(hnew,fouru*dabs(x)),h)
  465 h = hnew
      return
!       ***     end block 4     ***
      end
      subroutine intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)
      implicit real(a-h,o-z)
!
!   the methods in subroutine  step  approximate the solution near  x
!   by a polynomial.  subroutine  intrp  approximates the solution at
!   xout  by evaluating the polynomial there.  information defining this
!   polynomial is passed from  step  so  intrp  cannot be used alone.
!
!   this code is completely explained and documented in the text,
!   computer solution of ordinary differential equations:  the initial
!   value problem  by l. f. shampine and m. k. gordon.
!
!   input to intrp --
!
!   all floating point variables are double precision
!   the user provides storage in the calling program for the arrays in
!   the call list
       dimension y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)
!   and defines
!      xout -- point at which solution is desired.
!   the remaining parameters are defined in  step  and passed to  intrp
!   from that subroutine
!
!   output from  intrp --
!
!      yout(*) -- solution at  xout
!      ypout(*) -- derivative of solution at  xout
!   the remaining parameters are returned unaltered from their input
!   values.  integration with  step  may be continued.
!
      dimension g(13),w(13),rho(13)
      data g(1)/1.0d0/,rho(1)/1.0d0/
!
      hi = xout - x
      ki = kold + 1
      kip1 = ki + 1
!
!   initialize w(*) for computing g(*)
!
      do 5 i = 1,ki
        temp1 = i
    5   w(i) = 1.0d0/temp1
      term = 0.0d0
!
!   compute g(*)
!
      do 15 j = 2,ki
        jm1 = j - 1
        psijm1 = psi(jm1)
        gamma = (hi + term)/psijm1
        eta = hi/psijm1
        limit1 = kip1 - j
        do 10 i = 1,limit1
   10     w(i) = gamma*w(i) - eta*w(i+1)
        g(j) = w(1)
        rho(j) = gamma*rho(jm1)
   15   term = psijm1
!
!   interpolate
!
      do 20 l = 1,neqn
        ypout(l) = 0.0d0
   20   yout(l) = 0.0d0
      do 30 j = 1,ki
        i = kip1 - j
        temp2 = g(i)
        temp3 = rho(i)
        do 25 l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
   25     ypout(l) = ypout(l) + temp3*phi(l,i)
   30   continue
      do 35 l = 1,neqn
   35   yout(l) = y(l) + hi*yout(l)
      return
      end
