c
c
c     =====================================================
      subroutine bc2(meqn,mbc,mx,my,xlower,ylower,
     &               dx,dy,q,maux,aux,t,dt,mthbc)
c     =====================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my)
c     # to the ghost cells outside the region:
c     #   (i, 1-jbc)   for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (i, my+jbc)  for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (1-ibc, j)   for ibc = 1,mbc,  j = 1-mbc, my+mbc
c     #   (mx+ibc, j)  for ibc = 1,mbc,  j = 1-mbc, my+mbc
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension mthbc(4)
      
      dimension timeIC(500), pressure(500), p(0:2)

      
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      p0 = 101325.0
      c0 = sqrt(gammagas*p0/rhog) 
      p = p0
      ! Read pressure and time data from file
      open (25,file="a-pICtime.dat",action="read",status="old") 
      do i=1,500
          read(25,*) timeIC(i), pressure(i)
      enddo
      close(25)
      
      ! Look for correct value for pressure in data file
      ddt = timeIC(100) - timeIC(101)
      delt = max(ddt,dt)
      do k=1,500
        time2 = (timeIC(k) - 0.018758) 
        if (abs(time2 - t) <  delt) then
          ! Convert PSI to Pascals
          p(0) = p0 + 6894.75729*pressure(k)!*exp(-0.1*ycell**2)
        end if
        ! EXTRA RESOLUTION IS UNNECESSARY
!         if (abs(time2 - (t + dt)) <  delt/2) then
!           ! Convert PSI to Pascals
!           p(1) = p0 + 6894.75729*pressure(k)!*exp(-0.1*ycell**2)
!         end if
!         if (abs(time2 - (t + 2*dt)) <  delt/2) then
!           ! Convert PSI to Pascals
!           p(2) = p0 + 6894.75729*pressure(k)!*exp(-0.1*ycell**2)
!         end if
      end do
!       if (p(0) == p0 .and. t .ge. .0001) then
!         p(0)=p0 + 6894.75729*12
!       end if
      p(2) = p(0)
      p(1) = p(0)
      ! Assign corresponding pressure values to left boudary ghost cells
      do j = 1-mbc, my+mbc
        do ibc = 0, mbc  
          q(1,1-ibc,j) = rhog*(p(ibc)/p0)**(1/gammagas)
          q(2,1-ibc,j) = (2.0/(gammagas - 1.0))*(-c0 + 
     & sqrt(gammagas*p(ibc)/q(1,1-ibc,j)))
          q(3,1-ibc,j) = 0.d0
          q(4,1-ibc,j) = (p(ibc) + gammagas*pinfgas)/(gammagas - 1.0) +
     & (q(2,1-ibc,j)**2 + q(3,1-ibc,j)**2)/(2.0*q(1,1-ibc,j))
        end do
      end do
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do 115 j = 1-mbc, my+mbc
         do 115 ibc=1,mbc
            do 115 m=1,meqn
               q(m,1-ibc,j) = q(m,1,j)
  115       continue
      go to 199

  120 continue
c     # periodic:  
      do 125 j = 1-mbc, my+mbc
         do 125 ibc=1,mbc
            do 125 m=1,meqn
               q(m,1-ibc,j) = q(m,mx+1-ibc,j)
  125       continue
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 j = 1-mbc, my+mbc
         do 135 ibc=1,mbc
            do 135 m=1,meqn
               q(m,1-ibc,j) = q(m,ibc,j)
  135       continue
c     # negate the normal velocity:
      do 136 j = 1-mbc, my+mbc
         do 136 ibc=1,mbc
            q(2,1-ibc,j) = -q(2,ibc,j)
  136    continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do 215 j = 1-mbc, my+mbc
         do 215 ibc=1,mbc
            do 215 m=1,meqn
               q(m,mx+ibc,j) = q(m,mx,j)
  215       continue
      go to 299

  220 continue
c     # periodic:  
      do 225 j = 1-mbc, my+mbc
         do 225 ibc=1,mbc
            do 225 m=1,meqn
               q(m,mx+ibc,j) = q(m,ibc,j)
  225       continue
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 j = 1-mbc, my+mbc
         do 235 ibc=1,mbc
            do 235 m=1,meqn
               q(m,mx+ibc,j) = q(m,mx+1-ibc,j)
  235       continue
c     # negate the normal velocity:
      do 236 j = 1-mbc, my+mbc
         do 236 ibc=1,mbc
            q(2,mx+ibc,j) = -q(2,mx+1-ibc,j)
  236    continue
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
      stop
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 jbc=1,mbc
         do 315 i = 1-mbc, mx+mbc
            do 315 m=1,meqn
               q(m,i,1-jbc) = q(m,i,1)
  315       continue
      go to 399

  320 continue
c     # periodic:  
      do 325 jbc=1,mbc
         do 325 i = 1-mbc, mx+mbc
            do 325 m=1,meqn
               q(m,i,1-jbc) = q(m,i,my+1-jbc)
  325       continue
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 jbc=1,mbc
         do 335 i = 1-mbc, mx+mbc
            do 335 m=1,meqn
               q(m,i,1-jbc) = q(m,i,jbc)
  335       continue
c     # negate the normal velocity:
      do 336 jbc=1,mbc
         do 336 i = 1-mbc, mx+mbc
            q(3,i,1-jbc) = -q(3,i,jbc)
  336    continue
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do 415 jbc=1,mbc
         do 415 i = 1-mbc, mx+mbc
            do 415 m=1,meqn
               q(m,i,my+jbc) = q(m,i,my)
  415       continue
      go to 499

  420 continue
c     # periodic:  
      do 425 jbc=1,mbc
         do 425 i = 1-mbc, mx+mbc
            do 425 m=1,meqn
               q(m,i,my+jbc) = q(m,i,jbc)
  425       continue
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 jbc=1,mbc
         do 435 i = 1-mbc, mx+mbc
            do 435 m=1,meqn
               q(m,i,my+jbc) = q(m,i,my+1-jbc)
  435       continue
c     # negate the normal velocity:
      do 436 jbc=1,mbc
         do 436 i = 1-mbc, mx+mbc
            q(3,i,my+jbc) = -q(3,i,my+1-jbc)
  436    continue
      go to 499

  499 continue

      return
      end

