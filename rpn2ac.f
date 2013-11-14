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
      
      common /param2/ dt2, dx2,dy2

      common /comroe/ u,v,u2v2,enth,al,ar,g1a2,euv
c
c     # Dimensional splitting
      if (ixy.eq.1) then
          mu = 2
          mv = 3
          dxx = dx2
      else
          mu = 3
          mv = 2
          dxx = dy2
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
         ! Pressures (Use Tait EOS on water and/or plastic, SGEOS on air or SGEOS on both)
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
         cl = dsqrt(gammal*(pl + pinfl)/rho_l)
         cr = dsqrt(gammar*(pr + pinfr)/rho_r)
         
!          cl = dsqrt(gamma1l*(pl + qr(4,i-1) - ek_l)/rho_l)
!          cr = dsqrt(gamma1r*(pr + ql(4,i) - ek_r)/rho_r)

!          ! Compute the speed of left and right HLLC wave
         Sl = min(ul - cl,ur - cr) ! u(i) - a(i)
         Sr = max(ul + cl,ur + cr) ! u(i) + a(i),
! !          Sl = min(ul - cl,u(i) - cROE(i))
! !          Sr = max(u(i) + cROE(i), ur + cr) 
!          s(1,i) = 1.d0*Sl
!          s(3,i) = 1.d0*Sr

         ! Recompute Sm, Sl, and Sr using pressure wave speed estimates
         rhoavg = 0.5*(rho_l + rho_r)
         cavg = 0.5*(cl + cr)
         pstar = 0.5*(pl + pr) - 0.5*(ur - ul)*rhoavg*cavg
         
!          Sm = 0.5*(ul + ur) - 0.5*(pr - pl)/(rhoavg*cavg)
!          s(2,i) = 1.d0*Sm
         
!          !Pressure wave-speed estimates 
!          if (pl .ge. pstar) then
!               Sl = ul - cl
!          else
!               Sl = ul - cl*dsqrt(1.0 + (gammal + 1.0)/(2*gammal)*
!      & ((pstar + pinfl)/(pl + pinfl) - 1.0))
!          end if
!          if (pr .ge. pstar) then
!               Sr = ur + cr
!          else
!               Sr = ur + cr*dsqrt(1.0 + (gammar + 1.0)/(2*gammar)*
!      & ((pstar + pinfr)/(pr + pinfr) - 1.0))
!          end if
         s(1,i) = 1.d0*Sl
         s(3,i) = 1.d0*Sr

         ! Compute HLLC middle speed state (see research notebook)
         Sm = pr - pl + rho_r*ur*(ur - Sr) - rho_l*ul*(ul - Sl) 
         Sm = Sm/(rho_r*(ur - Sr) - rho_l*(ul - Sl))
         s(2,i) = 1.d0*Sm
!          if (ixy .eq. 1 .and. i .eq. floor(real(mx)/2.0)) then
!             print*, pl, pr, Sm, cl
!          end if
         
!          ! Average EOS parameter for two wave cut on interaface 
!          ! and recalculate speeds
!          Smold = 100000.0
!          if (gammal .ne. gammar) then
!             ! Use iterative method modifying gamma averages and 
!             ! CD speed, Sm
!             do while (abs(Sm - Smold) > .001)
!               Smold = 1.d0*Sm
!               frac = abs(Sm*dt2/dx2) !We dont know dt
!               if (Sm .ge. 0.0) then
!                 gammarn = gammal*frac + gammar*(1.0 - frac)
!                 gamma1rn = gammarn - 1.0
!                 pinfrn = pinfl*frac + pinfr*(1.0 - frac)
!                 gammaln = gammal
!                 gamma1ln = gamma1l
!                 pinfln = pinfl
!               end if
!               
!               if (Sm < 0.0) then
!                 gammaln = gammar*frac + gammal*(1.0 - frac)
!                 gamma1ln = gammaln - 1.0
!                 pinfln = pinfr*frac + pinfl*(1.0 - frac)
!                 gammarn = gammar
!                 gamma1rn = gamma1r
!                 pinfrn = pinfr
!               end if
!               
!               ! Recompute pressures
!               pl = gamma1ln*(qr(4,i-1) - ek_l) 
!               pl = pl - pinfln*gammaln
!               pr = gamma1rn*(ql(4,i) - ek_r) 
!               pr = pr - pinfrn*gammarn
!               
! ! !               ! Recompute speeds
! !               cROE(i) = (pl/rhsqrtl + pr/rhsqrtr) / rhsq2 + 
! !      &  0.5*((ur - ul)/rhsq2)**2
! !               gamma1ROE = (gamma1ln*rhsqrtl + gamma1rn*rhsqrtr) / rhsq2
! !               psiROE = (gamma1ln*(qr(4,i-1) - ek_l)/rhsqrtl +
! !      & gamma1rn*(ql(4,i) - ek_r)/rhsqrtr) / rhsq2
! !               cROE(i) = dsqrt(psiROE + gamma1ROE*cROE(i))
! !               
! !               cl = dsqrt(gamma1ln*(pl + qr(4,i-1) - ek_l)/rho_l)
! !               cr = dsqrt(gamma1rn*(pr + ql(4,i) - ek_r)/rho_r)
!               cl = dsqrt(gammal*(pl + pinfl)/rho_l)
!               cr = dsqrt(gammar*(pr + pinfr)/rho_r)
!               
! !               Sl = min(ul - cl,u(i) - cROE(i))
! !               Sr = max(u(i) + cROE(i), ur + cr)
!               Sl = min(ul - cl,ur - cr) ! u(i) - a(i)
!               Sr = max(ul + cl,ur + cr)
!               s(1,i) = 1.d0*Sl
!               s(3,i) = 1.d0*Sr
!               
!               Sm = pr - pl + rho_r*ur*(ur - Sr) - rho_l*ul*(ul - Sl) 
!               Sm = Sm/(rho_r*(ur - Sr) - rho_l*(ul - Sl))
!               s(2,i) = 1.d0*Sm
!             end do
! !             print*, gammar, gammarn
! !             gammar = gammarn
! !             gamma1r = gamma1rn
! !             pinfr = pinfrn
! !             gammal = gammaln
! !             gamma1l = gamma1ln
! !             pinfl = pinfln
!           end if
         
         
         ! Compute middle state pressure (both formulas yield the same value )
         pstar_l = pl + rho_l*(ul - Sm)*(ul - Sl)
         pstar_r = pr + rho_r*(ur - Sm)*(ur - Sr)
         
        ! Force zero speed at contact discontinuity
        ! makes wave two not move without affecting it's influence
        ! to wave 1 and 3 given by Sm
         if (gammal .ne. 1.4 .or. gammar .ne. 1.4) then
            !Sm = 0.0
!             go to 10
            s(2,i) = 0.0
            Sl = Sl - ustar
            Sr = Sr - ustar
            s(1,i) = 1.d0*Sl
            s(3,i) = 1.d0*Sr
          end if
!           s(2,i) = 0.0
          
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
 
          ! Do almost EXACT RIEMANN SOLVER for interface
   10     continue
!           if (gammal .ne. gammar) then
!   !           ! Newton's ,method
!               pstar = min(pl,pr) + 0.9*(max(pl,pr) - min(pl,pr)) !0.5*(pl + pr)
!               pold = pstar + 10
!               do while (abs(pstar - pold) > 0.0001)
!                 pold = pstar
!                 CALL phi_exact(gammal,gammar,pr,pl,rho_r,rho_l,
!      & pinfl,pinfr,pstar,phi,phi_prime,rhos_l,rhos_r,ustar)
!                 pstar = pstar - phi/phi_prime
!               end do
! 
! !             ! Bisection method to find pressure in exact solution
! !             pold = -10.0
! !             pleft = pl
! !             pright = pr
! !             pstar = 0.5*(pleft + pright)
! !             do while (abs(pstar - pold) > 0.01)
! !               pold = pstar
! !               CALL phi_exact(gammal,gammar,pr,pl,rho_r,rho_l,pinfl,
! !      & pinfr,pleft,phil,phi_prime,rhos_l,rhos_r,ustar)
! !               CALL phi_exact(gammal,gammar,pr,pl,rho_r,rho_l,pinfl,
! !      & pinfr,pright,phir,phi_prime,rhos_l,rhos_r,ustar)
! !               CALL phi_exact(gammal,gammar,pr,pl,rho_r,rho_l,pinfl,
! !      & pinfr,pstar,phi,phi_prime,rhos_l,rhos_r,ustar)
! !               if (phi*phil > 0) then
! !                 pleft = pstar
! !               else
! !                 pright = pstar
! !               end if
! !               pstar = 0.5*(pleft + pright)
! !             end do
!               
!               ! Compute the speed of left and right HLLC wave
!             betal = (pl + pinfl)*(gammal - 1.0)/(gammal + 1.0)
!             betar = (pr + pinfr)*(gammar - 1.0)/(gammar + 1.0)
!             alphal = 2.0/(rho_l*(gammal + 1.0))
!             alphar = 2.0/(rho_r*(gammar + 1.0))
!             Sl = ul - dsqrt((pstar + pinfl + betal)/alphal)/rho_l
!             Sr = ur + dsqrt((pstar + pinfr + betar)/alphar)/rho_r
! 
! !             ustar = 0.0
! 
!             s(1,i) = 1.d0*Sl
!             s(2,i) = 1.d0*ustar
!             s(3,i) = 1.d0*Sr              
!             
!             bl = (gammal + 1.0)/(gammal - 1.0)
!             br = (gammar + 1.0)/(gammar - 1.0)
! !             
! !             if (isnan(pstar)) then
! !                 print*, pl,pr
! !                 print*,1
! !                 read(*,*)
! !             end if
!             
!               ! Calculate densities, momentums and energys 
!             qml(1,i) = rhos_l !rho_l*(1 + bl*pstar/pl)/(pstar/pl + bl)
!             qmr(1,i) = rhos_r !rho_r*(1 + br*pstar/pr)/(pstar/pr + br)
!             qml(mu,i) = qml(1,i)*ustar
!             qmr(mu,i) = qmr(1,i)*ustar
!             qml(mv,i) = qml(1,i)*vl!*rho_l*(Sl - ul)/(Sl - ustar)
!             qmr(mv,i) = qmr(1,i)*vr!*rho_r*(Sr - ur)/(Sr - ustar)
!             qml(4,i) = (pstar + gammal*pinfl)/(gammal - 1.0) + 
!      & 0.5*(qml(mu,i)**2 + qml(mv,i)**2)/qml(1,i)
!             qmr(4,i) = (pstar + gammar*pinfr)/(gammar - 1.0) + 
!      & 0.5*(qmr(mu,i)**2 + qmr(mv,i)**2)/qmr(1,i)
! 
! !               qml(4,i) = rhos_l*(qr(4,i-1)/rho_l + 
! !      & (ustar - ul)*(ustar + pl/(rho_l*(Sl - ul))))
! !          qmr(4,i) = rhos_r*(ql(4,i)/rho_r + 
! !      & (ustar - ur)*(ustar + pr/(rho_r*(Sr - ur))))
!               
!           
! c        # Compute the 3 waves.
! c        j index over q variables
!           do j=1,meqn
!               q_l = qr(j,i-1)
!               q_r = ql(j,i)
!               wave(j,1,i) = qml(j,i) - q_l
!               wave(j,2,i) = qmr(j,i) - qml(j,i)
!               wave(j,3,i) = q_r - qmr(j,i) 
!           end do
!          end if
         
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

