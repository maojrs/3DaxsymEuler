c     ============================================
      subroutine b4step2(mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,t,dt,maux,aux)
c     ============================================
c
c     # called from claw2 before each call to step2.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.

c
c     # dummy routine 
c
c     
      ! for aux averaged time dependant EOS parameters in aux arrays
!       USE auxmodule
      
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      
      common /param2/ dt2, dx2
      
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ rhog,rhop,rhow
      
      dt2 = 1.d0*dt
      dx2 = 1.d0*dx
c
      ! GAUGES for pressure (choose point to see pressure as a function of time)
      ! Choose x value to obtain pressure as a function of time
      
      ! GAUGE 1
      xcell = -0.015
      ycell = 0.005
      i = floor((xcell - xlower)/dx + 0.5)
      j = floor((ycell - ylower)/dy + 0.5)

      ! Calculate pressure at point xcell
      gamma = aux(1,i,j)
      gamma1 = aux(1,i,j) - 1.0
      pinf = aux(2,i,j)
      rho = q(1,i,j)           ! density
      momx = q(2,i,j)           ! momentum
      momy = q(3,i,j)
      ene = q(4,i,j)           ! energy
      p = gamma1*(ene - 0.5*(momx*momx + momy*momy)/rho)
      p = p - gamma*pinf 
     
      ! Write Gauge data to file
      write (21,*) t, p
      
      ! GAUGE 2
      xcell = 0.0
      ycell = 0.005
      i = floor((xcell - xlower)/dx + 0.5)
      j = floor((ycell - ylower)/dy + 0.5)

      ! Calculate pressure at point xcell
      gamma = aux(1,i,j)
      gamma1 = aux(1,i,j) - 1.0
      pinf = aux(2,i,j)
      rho = q(1,i,j)           ! density
      momx = q(2,i,j)           ! momentum
      momy = q(3,i,j)
      ene = q(4,i,j)           ! energy
      p = gamma1*(ene - 0.5*(momx*momx + momy*momy)/rho)
      p = p - gamma*pinf 
     
      ! Write Gauge data to file
      write (22,*) t, p
      
!       ! DO average of moved aux arrays
!       do i=1,mx
!         do j=1,my
!             if (ustar_array(i,j,1) .ne. 0) then
!                 frac = ustar_array(i,j,1)*dt/dx
!                 if (ustar_array(i,j,1) > 0) then
!           aux(1,i,j) = frac*aux(1,i-1,j) + (1.0 - frac)*aux(1,i,j)
!           aux(2,i,j) = frac*aux(2,i-1,j) + (1.0 - frac)*aux(2,i,j)
!                 else
!           aux(1,i-1,j) = frac*aux(1,i,j) + (1.0 - frac)*aux(1,i-1,j)
!           aux(2,i-1,j) = frac*aux(2,i,j) + (1.0 - frac)*aux(2,i-1,j)
!                 end if
!             end if
!             if (ustar_array(i,j,2) .ne. 0) then
!                 frac = ustar_array(i,j,2)*dt/dy
!                 if (ustar_array(i,j,2) > 0) then
!           aux(1,i,j) = frac*aux(1,i,j-1) + (1.0 - frac)*aux(1,i,j)
!           aux(2,i,j) = frac*aux(2,i,j-1) + (1.0 - frac)*aux(2,i,j)
!                 else
!           aux(1,i,j-1) = frac*aux(1,i,j) + (1.0 - frac)*aux(1,i,j-1)
!           aux(2,i,j-1) = frac*aux(2,i,j) + (1.0 - frac)*aux(2,i,j-1)
!                 end if
!             end if
!         end do
!       end do


c
      return
      end

