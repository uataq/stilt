!#######################################################################
FUNCTION rtsafe(x1_in,x2_in,PP,TTEST,xacc)
!  rtsafe() is used to find TP above LCL
!  Iteration to find Temperature of an airparcel at
!  a given sat. equ. pot. temp (TTEST) and a given
!  pressure (PP); returns "in-cloud" Temperature within
!  accuracy xacc
!  THIS FUNCTION IS TAKEN FROM W.H. PRESS', "NUMERICAL RECIPES"
!
!  $Id: rtsafe.f,v 1.1 2009/10/26 15:36:55 jel Exp $
!
      REAL             :: rtsafe
      REAL, INTENT(IN) :: x1_in,x2_in,xacc,PP,TTEST

! Maximum allowed number of iterations.
      INTEGER, PARAMETER :: MAXIT=100
! Using a combination of Newton-Raphson and bisection,
! find the root of a function bracketed between x1 and x2.
! The root, returned as the function value rtsafe,
! will be refined until its accuracy is known within +/- xacc.
!  funcd is a user-supplied subroutine which returns
! both the function value and the first derivative of the function.
! CHG: here funcd is THES_T
      INTEGER j
      REAL :: x1,x2,df,dx,dxold,f,fh,fl,temp,xh,xl
      EXTERNAL THES_T

!-------------------------------------------------------------------------------
      x1 = x1_in
      x2 = x2_in
      call THES_T(PP,x1,TTEST,fl,df)
      call THES_T(PP,x2,TTEST,fh,df)
!      if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.))
!    *    pause 'root must be bracketed in rtsafe'
!     Increase max. temperature, if root not bracketed
      if(fh.lt.0.)then
         x1=x1+10
      call THES_T(PP,x1,TTEST,fl,df)
      end if
      if(fh.lt.0.)then
         x1=x1+10
      call THES_T(PP,x1,TTEST,fl,df)
      end if
      if(fh.lt.0.)then
         x1=x1+10
      call THES_T(PP,x1,TTEST,fl,df)
      end if

      if(fl.eq.0.)then
         rtsafe=x1
         return
      else if(fh.eq.0.)then
         rtsafe=x2
         return
!     Orient the search so that f(xl) < 0.
      else if(fl.lt.0.)then
         xl=x1
         xh=x2
      else
         xh=x1
         xl=x2
      end if
!     Initialize the guess for root,
      rtsafe=.5*(x1+x2)
!     the "stepsize before last,"
!dwen(20090825)      dxold=dabs(x2-x1)
      dxold=abs(x2-x1)
!     and the last step.
      dx=dxold
      call THES_T(PP,rtsafe,TTEST,f,df)
!     Loop over allowed iterations.
      do j=1,MAXIT

! TEST ONLY
!         IF(j.EQ.1)THEN
!            WRITE(45,*) 'j  dx rtsafe temp xacc'
!         END IF
!         WRITE(45,*)j,dx,rtsafe,temp,xacc
! TEST ONLY

         if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).gt.0.                 &
!dwen(20090825)     &       .or. dabs(2.*f).gt.dabs(dxold*df) ) then
            .or. abs(2.*f).gt.abs(dxold*df) ) then
!     Bisect if Newton out of range,
!     or not decreasing fast enough.

! TEST ONLY
!         WRITE(45,*)'Bisect if Newton out of range'

            dxold=dx
            dx=0.5*(xh-xl)
            rtsafe=xl+dx
            if(xl.eq.rtsafe)return
!     Change in root is negligible.
            else
!     Newton step acceptable. Take it.
               dxold=dx
               dx=f/df
               temp=rtsafe
               rtsafe=rtsafe-dx
               if(temp.eq.rtsafe)return
         end if
!dwen(20090825)         if(dabs(dx).lt.xacc) return
         if(abs(dx).lt.xacc) return
!     Convergence criterion.
         call THES_T(PP,rtsafe,TTEST,f,df)
!     The one new function evaluation per iteration.
         if(f.lt.0.)then
!     Maintain the bracket on the root.
            xl=rtsafe
         else
            xh=rtsafe
         end if
      end do
!      pause 'rtsafe exceeding maximum iterations'
      RETURN
END FUNCTION rtsafe
