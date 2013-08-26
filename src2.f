c
c      =======================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,t,dt)
c      =======================================================
c
      implicit double precision (a-h,o-z)
      dimension    q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension  yt(1:4), yout(1:4)
c
c     # dummy subroutine for use when equation has no source term.
c     # If method(5)=0 then this routine is never called, but its
c     # existence may be required by some compilers.
      do i=1-mbc,mx + mbc
        xcell = xlower + (i-0.5d0)*dx
        do j=1-mbc,my + mbc
            rcell = ylower + (j-0.5d0)*dy
            ! Equation of state variables             
            gamma = aux(1,i,j)
            pinf = aux(2,i,j)
            omega = aux(3,i,j)
            
            yt(1) = q(1,i,j)
            yt(2) = q(2,i,j)
            yt(3) = q(3,i,j)
            yt(4) = q(4,i,j)
            
            ! Do runge kutta fourth order
            CALL rk4(yt,yout,rcell,t,dt,gamma,pinf)
            yt = 1.0*yout
        end do
      end do
      return
      end

      ! RK4 subroutine, requires a functin FUNC that returns a vector
c      =======================================================
      subroutine rk4(yt,yout,r,t,dt,gamma,pinf)
      implicit double precision (a-h,o-z)
      dimension yt(1:4), yout(1:4)
      double precision, dimension(4) ::  k1,k2,k3,k4
      
      ! Returns k1,k2,k3 and k4 as k1 = f(y(t),t)...etc
      CALL RKFUNC(yt,t,r,k1,gamma,pinf)
      CALL RKFUNC(yt+k1/2.0, t + dt/2.0,r,k2,gamma,pinf)
      CALL RKFUNC(yt+k2/2.0, t + dt/2.0,r,k3,gamma,pinf)  
      CALL RKFUNC(yt+k3, t + dt,r,k4,gamma,pinf)
      ! Compute RK time step with the k's
      yout = yt + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0
      return
      end
c      =======================================================

c      =======================================================
      subroutine RKFUNC(yt,t,r,ki,gamma,pinf)
      implicit double precision (a-h,o-z)
      dimension yt(1:4)
      double precision, dimension(4) :: ki(1:4)
      
      ur = yt(2)/yt(1)
      uz = yt(3)/yt(1)
      p = (gamma - 1.0)*(yt(4) - 0.5*(ur**2 + uz**2)) - gamma*pinf
      
      ki(1) = -yt(2)/r
      ki(2) = -yt(2)**2/(yt(1)*r)
      ki(3) = -yt(2)*yt(3)/yt(1)
      ki(4) = ur*(yt(4) + p)/r
      
      return
      end
c      =======================================================

