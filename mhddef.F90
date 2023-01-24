        SUBROUTINE prim_to_cons(gamma,prim,cons,bx)
        implicit none        
        real, dimension(8):: prim
        real :: bx, gamma
        real, dimension(8):: cons

        cons(1)  = prim(1)
        cons(2) = prim(2)*prim(1)
        cons(3) = prim(3)*prim(1)
        cons(4) = prim(4)*prim(1)
        cons(5) = prim(5)/(gamma-1.) + prim(1)*0.5* &
             (prim(2)**2+prim(3)**2+prim(4)**2)  &
             + 0.5*(prim(6)**2+prim(7)**2+bx**2)
        cons(6) = prim(6)
        cons(7) = prim(7)
        cons(8) = bx ! Note that the code is written in such a way such that the
        ! magnetic field normal to the flow is the 8th variable.

        END SUBROUTINE prim_to_cons


        SUBROUTINE eigen_cons(gamma,d,vx,vy,vz,h,Bx,by,bz,Xfac,Yfac, &
                       lambda,rem,lem)
        IMPLICIT NONE
        real :: d, vx, vy, vz, h
        real :: Bx, by, bz, Xfac, Yfac
        real, dimension(7) :: lambda
        real, dimension(7,7) :: rem,lem
        ! LOCALS
        real :: btsq, bt_starsq, bt, bt_star, gamma
        real :: vsq, vax,  vaxsq, hp, twid_asq, ct2, tsum, tdif,cf2_cs2
        real :: cfsq, cfast, cssq, cslow
        real :: beta_y,beta_z,beta_ystar,beta_zstar,beta_starsq,vbeta
        real :: alpha_f, alpha_s, droot, s, twid_a
        real :: Qfast, Qslow, af_prime, as_prime, Afpbb, Aspbb, na
        real :: cff, css, af, as, Afpb, Aspb, vqstr, norm
        real :: Q_ystar, Q_zstar
        real,parameter :: eps=2.2204E-16 

        !-----------------------------------------------------------------------
        vsq = vx**2+vy**2+vz**2
        btsq = by**2+bz**2
        bt_starsq = (gamma-1. - (gamma-2.)*Yfac)*btsq
        bt=sqrt(btsq)
        bt_star = sqrt(bt_starsq)

        vaxsq = Bx**2/d
        vax = sqrt(vaxsq)

        hp = h - (vaxsq + btsq/d)
        twid_asq = ((gamma-1.)*(hp-.5*vsq)-(gamma-2.)*Xfac)
        twid_asq = MAX(twid_asq,eps)
        ct2 = bt_starsq/d
        tsum = vaxsq + ct2 + twid_asq
        tdif = vaxsq + ct2 - twid_asq
        cf2_cs2 = sqrt(tdif*tdif + 4.*twid_asq*ct2)
        cfsq = 0.5*(tsum+cf2_cs2)
        cfast = sqrt(cfsq)
        cssq = twid_asq*vaxsq/cfsq
        cslow = sqrt(cssq)

        if (bt .eq. 0) then
         beta_y = 1.
         beta_z = 0.
        else
         beta_y = by/bt
         beta_z = bz/bt
        endif
        beta_ystar = beta_y/sqrt((gamma-1.) - (gamma-2.)*Yfac)
        beta_zstar = beta_z/sqrt((gamma-1.) - (gamma-2.)*Yfac)

        beta_starsq = beta_ystar**2 + beta_zstar**2
        vbeta = vy*beta_ystar + vz*beta_zstar

        if ( (cfsq - cssq) .eq. 0.) then
         alpha_f = 1.0
         alpha_s = 0.0
        else if ( (twid_asq - cssq) .le. 0.) then
          alpha_f = 0.0
          alpha_s = 1.0
         else if ( (cfsq - twid_asq) .le. 0.) then
          alpha_f = 1.0
          alpha_s = 0.0
         else
          alpha_f = sqrt((twid_asq-cssq)/(cfsq-cssq))
          alpha_s = sqrt((cfsq-twid_asq)/(cfsq-cssq))
         endif

        ! compute Qs and As for eigenmatrices
         droot = sqrt(d)
         s  = SIGN(1.,bx)
         twid_a = sqrt(twid_asq)
         Qfast = s*cfast*alpha_f
         Qslow = s*cslow*alpha_s
         af_prime = twid_a*alpha_f/droot
         as_prime = twid_a*alpha_s/droot
         Afpbb = af_prime*bt_star*beta_starsq
         Aspbb = as_prime*bt_star*beta_starsq

        ! eigenvalues
        lambda(1) = vx - cfast
        lambda(2) = vx - vax
        lambda(3) = vx - cslow
        lambda(4) = vx
        lambda(5) = vx + cslow
        lambda(6) = vx + vax
        lambda(7) = vx + cfast        
       
        ! eigenmatrix
        rem(1,1) = alpha_f
        rem(1,2) = alpha_f*(vx-cfast)
        rem(1,3) = alpha_f*vy + Qslow*beta_ystar
        rem(1,4) = alpha_f*vz + Qslow*beta_zstar
        rem(1,5) = alpha_f*(hp-vx*cfast) + Qslow*vbeta+Aspbb
        rem(1,6) = as_prime*beta_ystar
        rem(1,7) = as_prime*beta_zstar

        rem(2,3) = -beta_z
        rem(2,4) =  beta_y
        rem(2,5) = -(vy*beta_z - vz*beta_y)
        rem(2,6) = -s*beta_z/droot
        rem(2,7) =  s*beta_y/droot

        rem(3,1) = alpha_s
        rem(3,2) = alpha_s*(vx-cslow)
        rem(3,3) = alpha_s*vy - Qfast*beta_ystar
        rem(3,4) = alpha_s*vz - Qfast*beta_zstar
        rem(3,5) = alpha_s*(hp - vx*cslow) - Qfast*vbeta - Afpbb
        rem(3,6) = -af_prime*beta_ystar
        rem(3,7) = -af_prime*beta_zstar

        rem(4,1) = 1.0
        rem(4,2) = vx
        rem(4,3) = vy
        rem(4,4) = vz
        rem(4,5) = 0.5*vsq + (gamma-2.)*Xfac/(gamma-1.)

        rem(5,1) =  alpha_s
        rem(5,2) =  alpha_s*(vx+cslow)
        rem(5,3) =  alpha_s*vy + Qfast*beta_ystar
        rem(5,4) =  alpha_s*vz + Qfast*beta_zstar
        rem(5,5) =  alpha_s*(hp+vx*cslow) + Qfast*vbeta - Afpbb
        rem(5,6) =  rem(3,6)
        rem(5,7) =  rem(3,7)

        rem(6,3) =  beta_z
        rem(6,4) = -beta_y
        rem(6,5) = -rem(2,5)
        rem(6,6) =  rem(2,6)
        rem(6,7) =  rem(2,7)

        rem(7,1) =  alpha_f
        rem(7,2) = alpha_f*(vx+cfast)
        rem(7,3) = alpha_f*vy - Qslow*beta_ystar
        rem(7,4) = alpha_f*vz - Qslow*beta_zstar
        rem(7,5) = alpha_f*(hp + vx*cfast) - Qslow*vbeta + Aspbb
        rem(7,6) =  rem(1,6)
        rem(7,7) =  rem(1,7)

        ! Left eignematrix
        ! normalize some of the quantities by 1/(2a^2)
        ! some by (gamma-1)/2a^2
        na  = 0.5/twid_asq
        cff = na*alpha_f*cfast
        css = na*alpha_s*cslow
        Qfast = Qfast*na
        Qslow = Qslow*na
        af = na*af_prime*d
        as = na*as_prime*d
        Afpb = na*af_prime*bt_star
        Aspb = na*as_prime*bt_star       

        alpha_f = (gamma-1.)*na*alpha_f
        alpha_s = (gamma-1.)*na*alpha_s
        Q_ystar = beta_ystar/beta_starsq
        Q_zstar = beta_zstar/beta_starsq
        vqstr   = (vy*Q_ystar+vz*Q_zstar)
        norm = (gamma-1.)*2.*na

        lem(1,1) = alpha_f*(vsq-hp) + cff*(cfast+vx)-Qslow*vqstr-Aspb
        lem(2,1) = -alpha_f*vx - cff
        lem(3,1) = -alpha_f*vy + Qslow*Q_ystar
        lem(4,1) = -alpha_f*vz + Qslow*Q_zstar
        lem(5,1) =  alpha_f
        lem(6,1) = as*Q_ystar - alpha_f*by
        lem(7,1) = as*Q_zstar - alpha_f*bz       

        lem(1,2) =  0.5*(vy*beta_z - vz*beta_y)
        lem(3,2) = -0.5*beta_z
        lem(4,2) =  0.5*beta_y
        lem(6,2) = -0.5*droot*beta_z*s
        lem(7,2) =  0.5*droot*beta_y*s

        lem(1,3) =  alpha_s*(vsq-hp) + css*(cslow+vx)+Qfast*vqstr+Afpb
        lem(2,3) = -alpha_s*vx - css
        lem(3,3) = -alpha_s*vy - Qfast*Q_ystar
        lem(4,3) = -alpha_s*vz - Qfast*Q_zstar
        lem(5,3) =  alpha_s
        lem(6,3) = -af*Q_ystar - alpha_s*by
        lem(7,3) = -af*Q_zstar - alpha_s*bz

        lem(1,4) =  1. - norm*(.5*vsq - (gamma-2.)*Xfac/(gamma-1.))
        lem(2,4) =  norm*vx
        lem(3,4) =  norm*vy
        lem(4,4) =  norm*vz
        lem(5,4) = -norm
        lem(6,4) =  norm*by
        lem(7,4) =  norm*bz

        lem(1,5) =  alpha_s*(vsq-hp) + css*(cslow-vx)-Qfast*vqstr+Afpb
        lem(2,5) = -alpha_s*vx + css
        lem(3,5) = -alpha_s*vy + Qfast*Q_ystar
        lem(4,5) = -alpha_s*vz + Qfast*Q_zstar
        lem(5,5) =  alpha_s
        lem(6,5) =  lem(6,3)
        lem(7,5) =  lem(7,3)

        lem(1,6) = -lem(1,2)
        lem(3,6) = -lem(3,2)
        lem(4,6) = -lem(4,2)
        lem(6,6) =  lem(6,2)
        lem(7,6) =  lem(7,2)

        lem(1,7) =  alpha_f*(vsq-hp) + cff*(cfast-vx)+Qslow*vqstr-Aspb
        lem(2,7) = -alpha_f*vx + cff
        lem(3,7) = -alpha_f*vy - Qslow*Q_ystar
        lem(4,7) = -alpha_f*vz - Qslow*Q_zstar
        lem(5,7) =  alpha_f
        lem(6,7) =  lem(6,1)
        lem(7,7) =  lem(7,1)

        END SUBROUTINE eigen_cons
