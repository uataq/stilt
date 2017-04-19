      subroutine proj_3d(stcprm, point, vect, enx,eny,enz)
!/*
! *  At a given point, resolves the components of a vector vect in the
! *  local coordinate system of the map projection.  It is assumed that
! *  vect and point is given in the 3-dimensional geocentric _map_
! *  coordinate system, rather than the North centered _geo_ system.
! *  returns these components as enx, eny, and enz.  Replaces vect with
! *  its projection on the tangent plane at point.
! */
! CHG(03/12/03)
!
! $Id: proj_3d.f,v 1.1 2010/10/26 19:03:41 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
      real*8 point(3),vect(3)
      dot_pr = 0.
      do k=1,3
        dot_pr = dot_pr + point(k)*vect(k)
      enddo
!/*
! *  dot_prod is the local vertical component of vect.  Note, point is
! *  assumed to be a unit vector.
! */
      enz = dot_pr
      do k=1,3
        vect(k) = vect(k) - dot_pr * point(k)
      enddo
!/*  vect has now been projected to a plane tangent to point.
! */
      fact = 1. + point(3)
      xi = vect(1) - vect(3) * point(1) / fact
      eta = vect(2) - vect(3) * point(2) / fact
!/*
! *  xi, eta represent components of vector on the projection plane,
! *  i.e. the 2-dimensional map system.
! */
      fact = stcprm(1) -1.
      if (fact .lt. 0.) then
!/*
! *  For Lambert Conformal and Mercator projections (gamma < 1.0) ,
! *  a rotation of the vector components is needed.
! */
        if ( dabs(point(3)) .lt. 1.) then
          theta = fact * datan2(point(2),point(1))
          cgthta = dcos(theta)
          sgthta = dsin(theta)
          fact = xi * cgthta - eta * sgthta
          eta = eta * cgthta + xi * sgthta
          xi = fact
        endif
      endif
!/*
! *  Now rotate xi, eta to the final map projection direction.
! */
      enx = eta * stcprm(13) - xi * stcprm(14)
      eny = - xi * stcprm(13) - eta * stcprm(14)
      return
      END
