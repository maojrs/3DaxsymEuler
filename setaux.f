c     ============================================
      subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # dummy routine when no auxiliary arrays
c
c     
      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
      
! !       Aux(i,1) plays the role of gamma, aux(i,2) could be another parameter
      ! Aux arrays establish values for the three free parameters in polystyrene EOS
      ! Polystyrene EOS based on Van der Waals approximation in Spender and Gilmore paper
      
      pwidth = 2.0
      do i=1-mbc,mx + mbc
        xcell = xlower + (i-0.5d0)*dx
        do j=1-mbc,my + mbc
  !         indexi = i + mbc
            ycell = ylower + (j-0.5d0)*dy
            
!             ! Plastic wall pwidth width air on left, water on right
!             if (abs(xcell) .le. pwidth) then
!                 aux(1,i,j) = gammaplas
!                 aux(2,i,j) = pinfplas
!                 aux(3,i,j) = omeplas
!             else if (xcell < -pwidth) then
!                 aux(1,i,j) = gammagas
!                 aux(2,i,j) = pinfgas
!                 aux(3,i,j) = omegas
!             else
!                 aux(1,i,j) = gammawat
!                 aux(2,i,j) = pinfwat
!                 aux(3,i,j) = omewat
!             end if
             
            ! Water box immersed in air
            if ((abs(xcell+0.0).le.0.0085).and. (ycell.le.0.0085)) then
                aux(1,i,j) = gammawat
                aux(2,i,j) = pinfwat
                aux(3,i,j) = omewat
            else 
                aux(1,i,j) = gammagas
                aux(2,i,j) = pinfgas
                aux(3,i,j) = omegas
            end if
            
!             ! Add hydrophone
!             if ((abs(xcell) .le. 0.005) .and. (j .eq. 1)) then
!                 aux(1,i,j) = gammaplas
!                 aux(2,i,j) = pinfplas
!                 aux(3,i,j) = omeplas
!             end if
            
            aux(4,i,j) = 1.0*i
            aux(5,i,j) = 1.0*j
            
            
!             if (xcell < -0.025) then
!                 aux(4,i,j) = 1.0
!             else
!                 aux(4,i,j) = 0.0
!             end if
          
        end do
      end do
c
      return
      end


