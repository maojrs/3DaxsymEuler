
c
c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Acoustics with smooth radially symmetric profile to test accuracy
c
       ! for aux averaged time dependant EOS parameters in aux arrays
!        USE auxmodule

       implicit double precision (a-h,o-z)
       dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
       dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
       dimension p(1-mbc:mx+mbc)
       dimension dist(500), pressure(500)
       

       common /param/ gammagas, gammaplas, gammawat
       common /param/ pinfgas,pinfplas,pinfwat
       common /param/ omegas,omeplas,omewat
       common /param/ rhog,rhop,rhow
       common /param/ dx2,dy2

       dx2 = 1.0*dx
       dy2 = 1.0*dy
       pi = 4.d0*datan(1.d0)
       
       ! Important for averaged time dependat array
!        CALL init_auxmodule(mx,my,mbc,maux,aux)

!        do 20 i=1,mx
!           xcell = xlower + (i-0.5d0)*dx
!           do 20 j=1,my
!              ycell = ylower + (j-0.5d0)*dy
! !              r = dsqrt(xcell**2 + ycell**2)
! !              if (dabs(r-0.5d0) .le. 0.2) then
! !                  density = 2.d0 + dcos(pi*(r - 0.5d0)/0.2)
! !              else
! !                  density = 1.d0
! !              endif
! 
! 	     q(1,i,j) = 1.d0 + 0.5d0*dexp(-8.d0*((xcell + 4.0)**2))
! 	     q(2,i,j) = 0.d0
!              q(3,i,j) = 0.d0
! 	     q(4,i,j) = q(1,i,j) / (gammagas - 1)
!   20         continue

      open (23,file="a-pIC.dat",action="read",status="old")
      open (11,file="a-pgauge1.dat",action="write",
     & status="replace")
      open (12,file="a-pgauge2.dat",action="write",
     & status="replace")
      open (13,file="a-pgauge3.dat",action="write",
     & status="replace")
      open (14,file="a-pgauge4.dat",action="write",
     & status="replace")
      open (25,file="a-pICtime.dat",action="read",status="old") 
!       
      ! Read pressure data from file
      do i=1,500
          read(23,*) dist(i), pressure(i)
!         print*, dist(i),pressure(i)
      enddo


      ! Write initial conditions
!       rhog = 1.0 !kg/m^3
!       rhop = 1050.0 !kg/m^3
!       rhow = 1000.0 !kg/m^3
      p0 = 101325.d0 !
      p = p0      ! Energy in water an plastic in air-water-plastic arrangement
      EgwpWat = (Egas0*(gammagas - 1.d0) - gammagas*pinfgas
     & + gammawat*pinfwat)/(gammawat - 1.d0)         
      EgwpPlas = (EgwpWat*(gammawat - 1.d0) - gammawat*pinfwat
     & + gammaplas*pinfplas)/(gammaplas - 1.d0)   
      c0 = sqrt(gammagas*p0/rhog) 
      
      ! Defines jump in energies to make pressure equal accross interfaces using SGEOS
      ! Setady state gas Energy
      Egas0 = (p0 + gammagas*pinfgas)/(gammagas - 1.d0)  
      ! Energy in water an plastic in air-water-plastic arrangement
      EgwpWat = (Egas0*(gammagas - 1.d0) - gammagas*pinfgas
     & + gammawat*pinfwat)/(gammawat - 1.d0)         
      EgwpPlas = (EgwpWat*(gammawat - 1.d0) - gammawat*pinfwat
     & + gammaplas*pinfplas)/(gammaplas - 1.d0)   
      ! Energy in water an plastic in air-plastic-water arrangement
      EgpwPlas = (Egas0*(gammagas - 1.d0) - gammagas*pinfgas
     & + gammaplas*pinfplas)/(gammaplas - 1.d0)         
      EgpwWat = (EgpwPlas*(gammaplas - 1.d0) - gammaplas*pinfplas
     & + gammaat*pinfwat)/(gammawat - 1.d0)   
      
      do 150 i=1-mbc,mx+mbc
        xcell = xlower + (i-0.5d0)*dx
        do 151 j=1-mbc,my+mbc
          ycell = ylower + (j-0.5d0)*dy
          
!           if (xcell .le. -0.025) then
!             p(i) = 184062.08
!           end if
          
!           ! Look for correct value for pressure in data file
!           ddx = dist(101) - dist(100)
!           do k=1,500
!             dist2 = (dist(k) - 10.0) 
!             if (abs(dist2 - xcell) <  ddx/2) then
!               ! Convert PSI to Pascals
!               p(i) = p0 + 6894.75729*pressure(k)!*exp(-0.1*ycell**2)
!             end if
!           end do
          
          ! Adjust initial conditions depending if it's gas,PS or water
          if (aux(1,i,j) == gammagas) then
            q(1,i,j) = rhog*(p(i)/p0)**(1/gammagas) !+ 5.0d0*dexp(-200.d0*(xcell+1.0)**2)
            q(2,i,j) = (2/(gammagas - 1.0))*(-c0 + 
     & sqrt(gammagas*p(i)/q(1,i,j)))
            q(3,i,j) = 0.d0
            !q(2,i,j) = 0.d0
            q(4,i,j) = (p(i) + gammagas*pinfgas)/(gammagas - 1.0) +
     & (q(2,i,j)**2 + q(3,i,j)**2)/(2.0*q(1,i,j))
     
         else if (aux(1,i,j) == gammaplas) then
            q(1,i,j) = rhop
            q(2,i,j) = 0.d0
            q(3,i,j) = 0.d0
            ! make sure pressure jump is zero across interface using SGEOS (check correct order of interfaces)
            q(4,i,j) = 1.d0*EgwpPlas
     
         else if (aux(1,i,j) == gammawat) then
            q(1,i,j) = rhow
            q(2,i,j) = 0.d0
            q(3,i,j) = 0.d0
            ! make sure pressure jump is zero across interface using SGEOS (check correct order of interfaces)
            q(4,i,j) = 1.d0*EgwpWat

          end if
  151    continue	  
  150  continue
       return
       end

