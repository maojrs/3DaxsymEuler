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
         xcell = -10.0 + (i-0.5d0)*.052
         if (ixy .eq. 2) then
            xcell = (i-0.5d0)*.05
         end if
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
         ! Pressures
         pl = gamma1l*(qr(4,i-1) - ek_l) 
         pl = pl/(1.0 - omel*rho_l) - pinfl*gammal
         pr = gamma1r*(ql(4,i) - ek_r) 
         pr = pr/(1.0 - omer*rho_r) - pinfr*gammar
         

         ! Additional qunatites to pass to transverse solver (ROE averages)
         rhsqrtl = dsqrt(qr(1,i-1))
         rhsqrtr = dsqrt(ql(1,i))
         rhsq2 = rhsqrtl + rhsqrtr
         u(i) = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
         v(i) = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
         enth(i) = (((qr(4,i-1)+pl)/rhsqrtl
     &             + (ql(4,i)+pr)/rhsqrtr)) / rhsq2
         u2v2(i) = u(i)**2 + v(i)**2
         a2l = gamma1l*(enth(i) - .5d0*u2v2(i))
         a2r = gamma1r*(enth(i) - .5d0*u2v2(i))
         al(i) = dsqrt(a2l)
         ar(i) = dsqrt(a2r)
         g1a2(i) = 1.d0/(enth(i) - .5d0*u2v2(i))
         euv(i) = enth(i) - u2v2(i) 
         cROE(i) = (pl/rhsqrtl + pr/rhsqrtr) / rhsq2 + 
     &  0.5*((ur - ul)/rhsq2)**2
         gamma1ROE = (gamma1l*rhsqrtl + gamma1r*rhsqrtr) / rhsq2
         psiROE = (gamma1l*(qr(4,i-1) - ek_l)/rhsqrtl +
     & gamma1r*(ql(4,i) - ek_r)/rhsqrtr) / rhsq2
         cROE(i) = dsqrt(psiROE + gamma1ROE*cROE(i))
    

         ! Compute left and right speeds
!          enth_l = (qr(4,i-1) + pl)/rho_l
!          enth_r = (ql(4,i) + pr)/rho_r
!          cl = dsqrt(gamma1l*(enth_l - .5d0*ul**2))
!          cr = dsqrt(gamma1r*(enth_r - .5d0*ur**2))
!          cl = dsqrt(gammal*(pl + pinfl)/rho_l)
!          cr = dsqrt(gammar*(pr + pinfr)/rho_r)
         
         cl = dsqrt(gamma1l*(pl + qr(4,i-1) - ek_l)/rho_l)
         cr = dsqrt(gamma1r*(pr + ql(4,i) - ek_r)/rho_r)

         ! Compute the speed of left and right HLLC wave
!          Sl = min(ul - cl,ur - cr) ! u(i) - a(i)
!          Sr = max(ul + cl,ur + cr) ! u(i) + a(i),
         Sl = min(ul - cl,u(i) - cROE(i))
         Sr = max(u(i) + cROE(i), ur + cr) 
         s(1,i) = 1.d0*Sl
         s(3,i) = 1.d0*Sr

         ! Compute HLLC middle speed state (see research notebook)
         Sm = pr - pl + rho_r*ur*(ur - Sr) - rho_l*ul*(ul - Sl) 
         Sm = Sm/(rho_r*(ur - Sr) - rho_l*(ul - Sl))
         s(2,i) = 1.d0*Sm
         
         ! Average EOS parameter for two wave cut on interaface 
         ! and recalculate speeds
!          Smold = 1000.0
!          if (gammal .ne. gammar) then
!             ! Use iterative method modifying gamma averages and 
!             ! CD speed, Sm
!             do while (abs(Sm - Smold) < .0001)
!               Smold = 1.d0*Sm
!               frac = abs(Sm*dt2/dx2)
!               if (Sm > 0.0) then
!                 gammarn = gammal*frac + gammar*(1.0 - frac)
!                 gamma1rn = gammarn - 1.0
!                 pinfrn = pinfl*frac + pinfr*(1.0 - frac)
!               end if
!               
!               if (Sm < 0.0) then
!                 gammaln = gammar*frac + gammal*(1.0 - frac)
!                 gamma1ln = gammaln - 1.0
!                 pinfln = pinfr*frac + pinfl*(1.0 - frac)
!               end if
!               
!               ! Recompute pressures
!               pl = gamma1ln*(qr(4,i-1) - ek_l) 
!               pl = pl/(1.0 - omel*rho_l) - pinfln*gammaln
!               pr = gamma1rn*(ql(4,i) - ek_r) 
!               pr = pr/(1.0 - omer*rho_r) - pinfrn*gammarn
!               
!               ! Recompute speeds
!               cROE(i) = (pl/rhsqrtl + pr/rhsqrtr) / rhsq2 + 
!      &  0.5*((ur - ul)/rhsq2)**2
!               gamma1ROE = (gamma1ln*rhsqrtl + gamma1rn*rhsqrtr) / rhsq2
!               psiROE = (gamma1ln*(qr(4,i-1) - ek_l)/rhsqrtl +
!      & gamma1rn*(ql(4,i) - ek_r)/rhsqrtr) / rhsq2
!               cROE(i) = dsqrt(psiROE + gamma1ROE*cROE(i))
!               
!               cl = dsqrt(gamma1ln*(pl + qr(4,i-1) - ek_l)/rho_l)
!               cr = dsqrt(gamma1rn*(pr + ql(4,i) - ek_r)/rho_r)
!               
!               Sl = min(ul - cl,u(i) - cROE(i))
!               Sr = max(u(i) + cROE(i), ur + cr) 
!               s(1,i) = 1.d0*Sl
!               s(3,i) = 1.d0*Sr
!               
!               Sm = pr - pl + rho_r*ur*(ur - Sr) - rho_l*ul*(ul - Sl) 
!               Sm = Sm/(rho_r*(ur - Sr) - rho_l*(ul - Sl))
!               s(2,i) = 1.d0*Sm
!             end do
!             if (Sm > 0.0) then
!                 gammar = gammarn
!                 gamma1r = gamma1rn
!                 pinfr = pinfrn
!               end if
!               
!               if (Sm < 0.0) then
!                 gammal = gammaln
!                 gamma1l = gamma1ln
!                 pinfl = pinfln
!               end if
!           end if
          
!         Force zero speed at contact discontinuity
!         makes wave two not move without affecting it's influence
!         to wave 1 and 3 given by Sm
         if (gammal .ne. gammar) then
            !Sm = 0.0
            s(2,i) = 0.0
          end if
         
         ! Compute middle state pressure (both formulas yield the same value )
         pstar_l = pl + rho_l*(ul - Sm)*(ul - Sl)
         pstar_r = pr + rho_r*(ur - Sm)*(ur - Sr)

         ! Calculate ql* and qr* HLLC middle states (without pressure term)
         do j=1,meqn
             qml(j,i) = rho_l*(Sl - ul)/(Sl - Sm)
             qmr(j,i) = rho_r*(Sr - ur)/(Sr - Sm)
         end do
         
         ! Add multiplicative terms to momentum ones
         qml(mu,i) = Sm*qml(mu,i)
         qmr(mu,i) = Sm*qmr(mu,i)
         qml(mv,i) = vl*qml(mv,i)
         qmr(mv,i) = vr*qmr(mv,i)
         ! Add second terms to energy one (see Toro pg. 325)
         qml(4,i) = qml(4,i)*(qr(4,i-1)/rho_l + 
     & (Sm - ul)*(Sm + pl/(rho_l*(Sl - ul))))
         qmr(4,i) = qmr(4,i)*(ql(4,i)/rho_r + 
     & (Sm - ur)*(Sm + pr/(rho_r*(Sr - ur))))

c        # Compute the 3 waves.
c        j index over q variables
         do j=1,meqn
             q_l = qr(j,i-1)
             q_r = ql(j,i)
             wave(j,1,i) = qml(j,i) - q_l
             wave(j,2,i) = qmr(j,i) - qml(j,i)
             wave(j,3,i) = q_r - qmr(j,i) 
         end do
!         wave(1,1,i) = qml(1,i) - qr(1,i-1)
!         wave(mu,1,i) = qml(mu,i) - qr(mu,i-1)
!         wave(mv,1,i) = qml(mv,i) - qr(mv,i-1) 
!         wave(4,1,i) = qml(4,i) - qr(4,i-1)
!         
!         wave(1,2,i) = qmr(1,i) - qml(1,i)
!         wave(mu,2,i) = qmr(mu,i) - qml(mu,i)
!         wave(mv,2,i) = qmr(mv,i) - qml(mv,i) 
!         wave(4,2,i) = qmr(4,i) - qml(4,i)
!         
!         wave(1,3,i) = ql(1,i) - qmr(1,i)
!         wave(mu,3,i) = ql(mu,i) - qmr(mu,i)
!         wave(mv,3,i) = ql(mv,i) - qmr(mv,i) 
!         wave(4,3,i) = ql(4,i) - qmr(4,i)

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

