c
c
c
c =========================================================
      subroutine rpn2(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &           wave,s,amdq,apdq,num_aux)
c =========================================================
c
c     # solve Riemann problems for the 2D Euler equations (normal solver) 
c     # using HLLC approximate Riemann solver
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves, 
c     #            s the speeds, 
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
c
      implicit double precision (a-h,o-z)
      double precision pstar, Sr, Sl
      dimension   ql(meqn, 1-mbc:maxmx+mbc)
      dimension   qr(meqn, 1-mbc:maxmx+mbc)
      dimension   qml(meqn, 1-mbc:maxmx+mbc)
      dimension   qmr(meqn, 1-mbc:maxmx+mbc)
      dimension    s(mwaves, 1-mbc:maxmx+mbc)
      dimension wave(meqn, mwaves, 1-mbc:maxmx+mbc)
      dimension amdq(meqn, 1-mbc:maxmx+mbc)
      dimension apdq(meqn, 1-mbc:maxmx+mbc)
      dimension auxl(num_aux,1-mbc:maxmx+mbc)
      dimension auxr(num_aux,1-mbc:maxmx+mbc)

c     # local storage
c     ---------------
      parameter (max2 = 20002)  !# assumes at most 2000 grid points with mbc=2
      dimension u(-1:max2),v(-1:max2),u2v2(-1:max2),enth(-1:max2),
     &       al(-1:max2),ar(-1:max2),g1a2(-1:max2),euv(-1:max2),
     &       cROE(-1:max2)
     
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
      
      common /param2/ dt2, dx2

      common /comroe/ u,v,u2v2,enth,al,ar,g1a2,euv
c
c     # Dimensional splitting
      if (ixy.eq.1) then
          mu = 2
          mv = 3
      else
          mu = 3
          mv = 2
      endif
      
c     # Compute Roe-averaged quantities:
      do 20 i=2-mbc,mx+mbc
!          xcell = -10.0 + (i-0.5d0)*.052
!          if (ixy .eq. 2) then
!             xcell = (i-0.5d0)*.05
!          end if
         gammal = auxr(1,i-1)
         gammar = auxl(1,i)
         gamma1l = gammal - 1.0
         gamma1r = gammar - 1.0
         pinfl = auxr(2,i-1)
         pinfr = auxl(2,i)
         omel = auxr(3,i-1)
         omer = auxl(3,i)
         ! Densities
         rho_l = qr(1,i-1)
         rho_r = ql(1,i)
         ! Velocities
         ul = qr(mu,i-1)/rho_l
         ur = ql(mu,i)/rho_r
         vl = qr(mv,i-1)/rho_l
         vr = ql(mv,i)/rho_r
         ! Kinetic Energy
         ek_l = 0.5*rho_l*(ul**2 + vl**2)
         ek_r = 0.5*rho_r*(ur**2 + vr**2)
         ! Pressures (Use Tait EOS on water and/or plastic, SGEOS on air or SGEOS on both)
         pl = gamma1l*(qr(4,i-1) - ek_l) - pinfl*gammal
         pr = gamma1r*(ql(4,i) - ek_r) - pinfr*gammar
    
         cl = dsqrt(gammal*(pl + pinfl)/rho_l)
         cr = dsqrt(gammar*(pr + pinfr)/rho_r)
                         
!         ! Force zero speed at contact discontinuity
!         ! makes wave two not move without affecting it's influence
!         ! to wave 1 and 3 given by Sm
!          if (gammal .ne. gammar) then
!             !Sm = 0.0
!             s(2,i) = 0.0
!           end if

        ! Newton's ,method
          pstar = 0.5*(pl + pr)
          pold = pstar + 10
          do while (abs(pstar - pold) > 0.0001)
            pold = pstar
            CALL phi_exact(gammal,gammar,pr,pl,rho_r,rho_l,ul,ur,
     & pinfl,pinfr,pstar,phi,phi_prime,rhos_l,rhos_r,ustar)
            pstar = pstar - phi/phi_prime
          end do
              
!             ! Bisection method to find pressure in exact solution
!             pold = -10.0
!             pleft = pl
!             pright = pr
!             pstar = 0.5*(pleft + pright)
!             do while (abs(pstar - pold) > 0.0001)
!               pold = pstar
!               CALL phi_exact(gammal,gammar,pr,pl,rho_r,rho_l,pinfl,
!      & pinfr,pleft,phil,phi_prime,rhos_l,rhos_r,ustar)
!               CALL phi_exact(gammal,gammar,pr,pl,rho_r,rho_l,pinfl,
!      & pinfr,pright,phir,phi_prime,rhos_l,rhos_r,ustar)
!               CALL phi_exact(gammal,gammar,pr,pl,rho_r,rho_l,pinfl,
!      & pinfr,pstar,phi,phi_prime,rhos_l,rhos_r,ustar)
!               if (phi*phil > 0) then
!                 pleft = pstar
!               else
!                 pright = pstar
!               end if
!               pstar = 0.5*(pleft + pright)
!             end do
                       
            ! Compute the speed of left and right wave (See MJ IVINGS paper)
            betal = (pl + pinfl)*(gammal - 1.0)/(gammal + 1.0)
            betar = (pr + pinfr)*(gammar - 1.0)/(gammar + 1.0)
            alphal = 2.0/(rho_l*(gammal + 1.0))
            alphar = 2.0/(rho_r*(gammar + 1.0))
            Sl = ul - dsqrt((pstar + pinfl + betal)/alphal)/rho_l
            Sr = ur + dsqrt((pstar + pinfr + betar)/alphar)/rho_r
!             Sl2 = min(ul - cl,ur - cr) ! u(i) - a(i)
!             Sr2 = max(ul + cl,ur + cr) ! u(i) + a(i),
!             Sl = ul - cl
!             Sr = ur + cr

            
            s(1,i) = 1.d0*Sl
            s(2,i) = 0.0*ustar
            s(3,i) = 1.d0*Sr
            
!             if (gammal .ne. gammar) then
!                 s(2,i) = 0.0
!             end if
            
            
            bl = (gammal + 1.0)/(gammal - 1.0)
            br = (gammar + 1.0)/(gammar - 1.0)
            ! Calculate densities, momentums and energys
            qml(1,i) = rhos_l !rho_l*(1 + bl*pstar/pl)/(pstar/pl + bl)
            qmr(1,i) = rhos_r !rho_r*(1 + br*pstar/pr)/(pstar/pr + br)
            qml(mu,i) = qml(1,i)*ustar
            qmr(mu,i) = qmr(1,i)*ustar
            qml(mv,i) = qml(1,i)*vl !*rho_l*(Sl - ul)/(Sl - ustar)
            qmr(mv,i) = qmr(1,i)*vr!*rho_r*(Sr - ur)/(Sr - ustar)
            qml(4,i) = (pstar + gammal*pinfl)/(gammal - 1.0) + 
     & 0.5*(qml(mu,i)**2 + qml(mv,i)**2)/qml(1,i)
            qmr(4,i) = (pstar + gammar*pinfr)/(gammar - 1.0) + 
     & 0.5*(qmr(mu,i)**2 + qmr(mv,i)**2)/qmr(1,i)
            
!             if (qmr(mu,i) .ne. 0.0) then
!               print*, rhos_r, qmr(mu,i), qmr(4,i)
!             end if
           
            
!               if (isnan(pstar)) then
!                 print*, Sl,Sr
!                 write(*,'(" 1 ")') !BP
!                 read(*,*)
!               end if
!               
                         
        
c        # Compute the 3 waves.
c        j index over q variables
        do j=1,meqn
            q_l = qr(j,i-1)
            q_r = ql(j,i)
            wave(j,1,i) = qml(j,i) - q_l
            wave(j,2,i) = qmr(j,i) - qml(j,i)
            wave(j,3,i) = q_r - qmr(j,i) 
        end do
         
   20    continue
 

c
c
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves
c
      do 100 m=1,meqn
         do 100 i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mw=1,mwaves
               if (s(mw,i) .lt. 0.d0) then
                   amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
   90          continue
  100       continue
      go to 900
c
c-----------------------------------------------------

c

c
  900 continue
      return
      end

