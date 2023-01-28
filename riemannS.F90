        subroutine riemannS(fluxMth,w,wL,wR,gammaa,normal,flux)
        implicit none
        character(len = 15):: fluxMth
        real:: gammaa
        real:: w(8),wL(8),wR(8),flux(8)
        real:: wL_int(8),wR_int(8),flux_int(8)
        integer:: normal(2)

        wL_int=wL
        wR_int=wR
        if (normal(1).eq.1) then        
        else if (normal(2).eq.1) then
        wL_int(2)=wL(3)
        wL_int(3)=wL(2)
        wL_int(5)=wL(6)
        wL_int(6)=wL(5)
        wR_int(2)=wR(3)
        wR_int(3)=wR(2)
        wR_int(5)=wR(6)
        wR_int(6)=wR(5)
        end if
        
        
        if (fluxMth.eq.'RUSA') then
         call RUSA(wL_int,wR_int,gammaa,[1,0],flux_int)
        else if (fluxMth.eq.'HLLE') then
         call HLLE(wL_int,wR_int,gammaa,[1,0],flux_int)
        else if (fluxMth.eq.'HLLC_Linde') then
         call HLLC_Linde(wL_int,wR_int,gammaa,[1,0],flux_int)
        else if (fluxMth.eq.'HLLC') then
         call HLLC(wL_int,wR_int,gammaa,[1,0],flux_int)
        else if (fluxMth.eq.'HLLD') then
         call HLLD(wL_int,wR_int,gammaa,[1,0],flux_int)
        else if (fluxMth.eq.'Roe') then
         write(6,*) 'Still not ready!'
         call RoeS(wL_int,wR_int,gammaa,[1,0],flux_int)
        else
         write(6,*) 'No selection chosen! Stop process...'
         stop
        end if        

        flux=flux_int
        if (normal(1).eq.1) then
        else if (normal(2).eq.1) then
        flux(2)=flux_int(3)
        flux(3)=flux_int(2)
        flux(5)=flux_int(6)
        flux(6)=flux_int(5)
        end if

        end subroutine riemannS


!-------------------------------------------------------------------
! Here place all the subroutines related to Riemann solver
!-------------------------------------------------------------------

! Rusanov flux--------------------------------------------------------
        subroutine RUSA(wL,wR,gammaa,normal,flux)
        implicit none
        real:: gammaa
        real:: wL(8),wR(8),flux(8)
        real:: vnL,BnL,pL,caL,aL,canL,cfL,eL,HL
        real:: vnR,BnR,pR,caR,aR,canR,cfR,eR,HR
        real:: psL,psR
        real:: BxLs,BxRs,ByLs,ByRs,BzLs,BzRs,BnLs,BnRs
        real:: ptL,ptR
        real:: qL(8),FL(8),qR(8),FR(8)
        integer:: nx,ny,nz,normal(2)

        ! Allocate zero
        flux=0

        ! normal vectors        
        nx=normal(1)
        ny=normal(2)
        nz=0

        ! Left state
        vnL=wL(2)*nx+wL(3)*ny+wL(4)*nz
        BnL=wL(5)*nx+wL(6)*ny+wL(7)*nz
        caL=((wL(5)**2+wL(6)**2+wL(7)**2)/wL(1))**0.5
        aL=(gammaa*wL(8)/wL(1))**0.5
        canL=((BnL**2)/wL(1))**0.5
        cfL=(0.5*(aL**2+caL**2)+0.5*((aL**2+caL**2)**2 &
                -4*(aL**2)*canL**2)**0.5)**0.5
        eL=(wL(8)/(wL(1)*(gammaa-1))) &
         +0.5*(wL(2)**2+wL(3)**2+wL(4)**2)+ &
         +0.5*(wL(5)**2+wL(6)**2+wL(7)**2)/wL(1)
        eL=eL*wL(1)
        ! Specific enthalpy H 
        HL=(eL+wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2))/wL(1)
        qL=[wL(1),wL(1)*wL(2),wL(1)*wL(3),wL(1)*wL(4), &
            wL(5),wL(6),wL(7),eL]

        ! Right state
        vnR=wR(2)*nx+wR(3)*ny+wR(4)*nz
        BnR=wR(5)*nx+wR(6)*ny+wR(7)*nz
        caR=((wR(5)**2+wR(6)**2+wR(7)**2)/wR(1))**0.5
        aR=(gammaa*wR(8)/wR(1))**0.5
        canR=((BnR**2)/wR(1))**0.5
        cfR=(0.5*(aR**2+caR**2)+0.5*((aR**2+caR**2)**2 &
               -4*(aR**2)*canR**2)**0.5)**0.5
        eR=(wR(8)/(wR(1)*(gammaa-1))) &
         +0.5*(wR(2)**2+wR(3)**2+wR(4)**2)+ &
         +0.5*(wR(5)**2+wR(6)**2+wR(7)**2)/wR(1)
        eR=eR*wR(1)
        ! Specific enthalpy H 
        HR=(eR+wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2))/wR(1)
        qR=[wR(1),wR(1)*wR(2),wR(1)*wR(3),wR(1)*wR(4), &
            wR(5),wR(6),wR(7),eR]
 
        ! Left fluxes 
        ptL=wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2)
        FL=[wL(1)*vnL, wL(1)*vnL*wL(2)+ptL*nx-BnL*wL(5), &
           wL(1)*vnL*wL(3)+ptL*ny-BnL*wL(6), &
           wL(1)*vnL*wL(4)+ptL*nz-BnL*wL(7), &
           vnL*wL(5)-BnL*wL(2),vnL*wL(6)-BnL*wL(3),vnL*wL(7)-BnL*wL(4), &
           wL(1)*vnL*HL-BnL*(wL(2)*wL(5)+wL(3)*wL(6)+wL(4)*wL(7))]
        ! Right fluxes
        ptR=wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2)
        FR=[wR(1)*vnR, wR(1)*vnR*wR(2)+ptR*nx-BnR*wR(5), &
           wR(1)*vnR*wR(3)+ptR*ny-BnR*wR(6), &
           wR(1)*vnR*wR(4)+ptR*nz-BnR*wR(7), &
           vnR*wR(5)-BnR*wR(2),vnR*wR(6)-BnR*wR(3),vnR*wR(7)-BnR*wR(4), &
           wR(1)*vnR*HR-BnR*(wR(2)*wR(5)+wR(3)*wR(6)+wR(4)*wR(7))]

        ! Compute the HLL flux
         flux=0.5*(FL+FR)-0.5*maxval([cfL+abs(wL(2)),cfR+abs(wR(2))]) &
                       *(qR-qL)

        if (nx.eq.1) then
         flux(5)=0.
        else 
         write(6,*) 'Directional error'
         stop
        end if

        end subroutine RUSA

! HLLE flux-----------------------------------------------------------
        subroutine HLLE(wL,wR,gammaa,normal,flux)
        implicit none
        real:: gammaa
        real:: wL(8),wR(8),flux(8)
        real:: RT,r,u,v,w,Bx,By,Bz,H,vn
        real:: vnL,BnL,pL,caL,aL,canL,cfL,eL,HL
        real:: vnR,BnR,pR,caR,aR,canR,cfR,eR,HR
        real:: psL,psR
        real:: BxLs,BxRs,ByLs,ByRs,BzLs,BzRs,BnLs,BnRs
        real:: SL,SR,ptL,ptR
        real:: qL(8),FL(8),qR(8),FR(8)
        integer:: nx,ny,nz,normal(2)
        real:: Xfactor,Yfactor,lambda(7)

        ! Allocate zero
        flux=0

        ! normal vectors        
        nx=normal(1)
        ny=normal(2)
        nz=0

        ! Left state
        vnL=wL(2)*nx+wL(3)*ny+wL(4)*nz
        BnL=wL(5)*nx+wL(6)*ny+wL(7)*nz
        caL=((wL(5)**2+wL(6)**2+wL(7)**2)/wL(1))**0.5
        aL=(gammaa*wL(8)/wL(1))**0.5
        canL=((BnL**2)/wL(1))**0.5
        cfL=(0.5*(aL**2+caL**2)+0.5*((aL**2+caL**2)**2 &
                -4*(aL**2)*canL**2)**0.5)**0.5
        eL=(wL(8)/(wL(1)*(gammaa-1))) &
         +0.5*(wL(2)**2+wL(3)**2+wL(4)**2)+ &
         +0.5*(wL(5)**2+wL(6)**2+wL(7)**2)/wL(1)
        eL=eL*wL(1)
        ! Specific enthalpy H 
        HL=(eL+wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2))/wL(1)
        qL=[wL(1),wL(1)*wL(2),wL(1)*wL(3),wL(1)*wL(4), &
            wL(5),wL(6),wL(7),eL]

        ! Right state
        vnR=wR(2)*nx+wR(3)*ny+wR(4)*nz
        BnR=wR(5)*nx+wR(6)*ny+wR(7)*nz
        caR=((wR(5)**2+wR(6)**2+wR(7)**2)/wR(1))**0.5
        aR=(gammaa*wR(8)/wR(1))**0.5
        canR=((BnR**2)/wR(1))**0.5
        cfR=(0.5*(aR**2+caR**2)+0.5*((aR**2+caR**2)**2 &
                -4*(aR**2)*canR**2)**0.5)**0.5
        eR=(wR(8)/(wR(1)*(gammaa-1))) &
         +0.5*(wR(2)**2+wR(3)**2+wR(4)**2)+ &
         +0.5*(wR(5)**2+wR(6)**2+wR(7)**2)/wR(1)
        eR=eR*wR(1)
        ! Specific enthalpy H 
        HR=(eR+wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2))/wR(1)
        qR=[wR(1),wR(1)*wR(2),wR(1)*wR(3),wR(1)*wR(4), &
            wR(5),wR(6),wR(7),eR]

        ! First compute the Roe Averages
        RT=sqrt(wR(1)/wL(1))
        r=sqrt(wR(1)*wL(1))
        u=(wL(2)+RT*wR(2))/(1+RT)
        v=(wL(3)+RT*wR(3))/(1+RT)
        w=(wL(4)+RT*wR(4))/(1+RT)
        Bx=(sqrt(wR(1))*wL(5)+sqrt(wL(1))*wR(5)) &  
           /(sqrt(wL(1))+sqrt(wR(1)))                    
        By=(sqrt(wR(1))*wL(6)+sqrt(wL(1))*wR(6)) &
           /(sqrt(wL(1))+sqrt(wR(1)))
        Bz=(sqrt(wR(1))*wL(7)+sqrt(wL(1))*wR(7)) &
           /(sqrt(wL(1))+sqrt(wR(1)))
        H=(HL+RT*HR)/(1+RT) ! correct same as Athena
        Xfactor=0.5*((wL(6)-wR(6))**2+(wL(7)-wR(7))**2)/ &
               (sqrt(wL(1))+sqrt(wR(1)))
        Yfactor=0.5*(wL(1)+wR(1))/r
        call RoeEigen(gammaa,r,u,v,w,H,Bx,By,Bz,Xfactor,Yfactor,lambda)

        vn=u*nx+v*ny+w*nz

        ! Wave speed estimates
        SL=minval([vnL-cfL,lambda(1),0.])
        SR=maxval([vnR+cfR,lambda(7),0.])

        ! Left fluxes 
        ptL=wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2)
        FL=[wL(1)*vnL, wL(1)*vnL*wL(2)+ptL*nx-BnL*wL(5), &
           wL(1)*vnL*wL(3)+ptL*ny-BnL*wL(6), &
           wL(1)*vnL*wL(4)+ptL*nz-BnL*wL(7), &
           vnL*wL(5)-BnL*wL(2), & 
           vnL*wL(6)-BnL*wL(3),vnL*wL(7)-BnL*wL(4), &
           wL(1)*vnL*HL-BnL*(wL(2)*wL(5)+wL(3)*wL(6)+wL(4)*wL(7))]
        ! Right fluxes
        ptR=wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2)
        FR=[wR(1)*vnR, wR(1)*vnR*wR(2)+ptR*nx-BnR*wR(5), &
           wR(1)*vnR*wR(3)+ptR*ny-BnR*wR(6), &
           wR(1)*vnR*wR(4)+ptR*nz-BnR*wR(7), &
           vnR*wR(5)-BnR*wR(2),vnR*wR(6)-BnR*wR(3),vnR*wR(7)-BnR*wR(4), &
           wR(1)*vnR*HR-BnR*(wR(2)*wR(5)+wR(3)*wR(6)+wR(4)*wR(7))]

        ! Compute the HLL flux 
        flux=(SR*FL-SL*FR+SL*SR*(qR-qL))/(SR-SL)
        if (nx.eq.1) then
         flux(5)=0.
        else
         write(6,*) 'Directional error'
         stop
        end if

        end subroutine HLLE


! HLLC flux from T.J.Linde--------------------------------------------
        subroutine HLLC_Linde(wL,wR,gammaa,normal,flux)
        implicit none
        integer:: version
        real:: gammaa
        real:: wL(8),wR(8),flux(8)
        real:: SM
        real:: interm
        real:: pL,caL,aL,canL,cfL,eL,HL
        real:: pR,caR,aR,canR,cfR,eR,HR
        real:: psL,psR
        real:: BxLs,BxRs,ByLs,ByRs,BzLs,BzRs,BnLs,BnRs
        real:: SL,SR,ptL,ptR
        real:: qL(8),FL(8),qR(8),FR(8)
        real:: qsL(8),qsR(8)
        real:: FLs(8),FRs(8),alpha
! Declaration of normal variables
        integer:: nx,ny,nz,normal(2)
        real:: vnL,BnL,vnR,BnR
! Declaration of Roe-averaged variables
        real:: u_roe,v_roe,w_roe,vn_roe,r,Bx,By,Bz,RT,H
        real:: Xfactor,Yfactor,lambda(7),cf_roe

        ! Allocate zero
        flux=0

        ! normal vectors        
        nx=normal(1)
        ny=normal(2)
        nz=0

        ! Left state
        vnL=wL(2)*nx+wL(3)*ny+wL(4)*nz
        BnL=wL(5)*nx+wL(6)*ny+wL(7)*nz
        caL=((wL(5)**2+wL(6)**2+wL(7)**2)/wL(1))**0.5
        aL=(gammaa*wL(8)/wL(1))**0.5
        canL=((BnL**2)/wL(1))**0.5
        cfL=(0.5*(aL**2+caL**2)+0.5*((aL**2+caL**2)**2 &
                -4*(aL**2)*canL**2)**0.5)**0.5
        eL=(wL(8)/(wL(1)*(gammaa-1))) &
         +0.5*(wL(2)**2+wL(3)**2+wL(4)**2)+ &
         +0.5*(wL(5)**2+wL(6)**2+wL(7)**2)/wL(1)
        eL=eL*wL(1)
        ! Specific enthalpy H 
        HL=(eL+wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2))/wL(1)
        qL=[wL(1),wL(1)*wL(2),wL(1)*wL(3),wL(1)*wL(4), &
            wL(5),wL(6),wL(7),eL]

        ! Right state
        vnR=wR(2)*nx+wR(3)*ny+wR(4)*nz
        BnR=wR(5)*nx+wR(6)*ny+wR(7)*nz
        caR=((wR(5)**2+wR(6)**2+wR(7)**2)/wR(1))**0.5
        aR=(gammaa*wR(8)/wR(1))**0.5
        canR=((BnR**2)/wR(1))**0.5
        cfR=(0.5*(aR**2+caR**2)+0.5*((aR**2+caR**2)**2 &
                -4*(aR**2)*canR**2)**0.5)**0.5
        eR=(wR(8)/(wR(1)*(gammaa-1))) &
         +0.5*(wR(2)**2+wR(3)**2+wR(4)**2)+ &
         +0.5*(wR(5)**2+wR(6)**2+wR(7)**2)/wR(1)
        eR=eR*wR(1)
        ! Specific enthalpy H 
        HR=(eR+wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2))/wR(1)
        qR=[wR(1),wR(1)*wR(2),wR(1)*wR(3),wR(1)*wR(4), &
            wR(5),wR(6),wR(7),eR]

        ! First compute the Roe Averages
        u_roe=(sqrt(wL(1))*wL(2)+sqrt(wR(1))*wR(2)) &
           /(sqrt(wL(1))+sqrt(wR(1)))
        v_roe=(sqrt(wL(1))*wL(3)+sqrt(wR(1))*wR(3)) &
           /(sqrt(wL(1))+sqrt(wR(1)))
        w_roe=(sqrt(wL(1))*wL(4)+sqrt(wR(1))*wR(4)) &
           /(sqrt(wL(1))+sqrt(wR(1)))
        vn_roe=u_roe*nx+v_roe*ny+w_roe*nz

        RT=sqrt(wR(1)/wL(1))
        r=sqrt(wR(1)*wL(1))
        Bx=(sqrt(wR(1))*wL(5)+sqrt(wL(1))*wR(5)) &
            /(sqrt(wL(1))+sqrt(wR(1)))                             
        By=(sqrt(wR(1))*wL(6)+sqrt(wL(1))*wR(6)) &
            /(sqrt(wL(1))+sqrt(wR(1)))
        Bz=(sqrt(wR(1))*wL(7)+sqrt(wL(1))*wR(7)) &
            /(sqrt(wL(1))+sqrt(wR(1)))
        H=(HL+RT*HR)/(1+RT)       
        Xfactor=0.5*((wL(6)-wR(6))**2+(wL(7)-wR(7))**2)/ &
                (sqrt(wL(1))+sqrt(wR(1)))
        Yfactor=0.5*(wL(1)+wR(1))/r
        call RoeEigen(gammaa,r,u_roe,v_roe,w_roe,H,Bx,By,Bz &
                ,Xfactor,Yfactor,lambda)

        ! Wave speed estimates
        SL=minval([vnL-cfL,lambda(1)])
        SR=maxval([vnR+cfR,lambda(7)])
        cf_roe=lambda(7)-vn_roe

        ! Left fluxes 
        ptL=wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2)
        FL=[wL(1)*vnL, wL(1)*vnL*wL(2)+ptL*nx-BnL*wL(5), &
           wL(1)*vnL*wL(3)+ptL*ny-BnL*wL(6), &
           wL(1)*vnL*wL(4)+ptL*nz-BnL*wL(7), &
           vnL*wL(5)-BnL*wL(2) & 
           ,vnL*wL(6)-BnL*wL(3),vnL*wL(7)-BnL*wL(4), &
           wL(1)*vnL*HL-BnL*(wL(2)*wL(5)+wL(3)*wL(6)+wL(4)*wL(7))]
        ! Right fluxes
        ptR=wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2)
        FR=[wR(1)*vnR, wR(1)*vnR*wR(2)+ptR*nx-BnR*wR(5), &
           wR(1)*vnR*wR(3)+ptR*ny-BnR*wR(6), &
           wR(1)*vnR*wR(4)+ptR*nz-BnR*wR(7), &
           vnR*wR(5)-BnR*wR(2),vnR*wR(6)-BnR*wR(3),vnR*wR(7)-BnR*wR(4), &
           wR(1)*vnR*HR-BnR*(wR(2)*wR(5)+wR(3)*wR(6)+wR(4)*wR(7))]

        ! As SM & alpha are loosely defined, one can construct various
        ! versions of the Linde solver by introducing adequate SM and
        ! alpha(=strength of middle wave)
        version=1 ! Choose type of Linde solver
!-------------------------------------------------------------------------------
        if (version.eq.1) then
         ! From Linde 1998 paper (use Roe average and heuristic linear
         ! function w.r.t. gap of jump)
         ! Code I wrote it based on Gurski's paper (2004)
         SM=vn_roe
         interm=1-sum(abs((FR-FL)- &
                vn_roe*(qR-qL)))/(cf_roe*sum(abs(qR-qL)))
         alpha=maxval([0.,interm])
         if (alpha.gt.1.or.alpha.lt.0.or.alpha.ne.alpha) &
                write(6,*)'Error as alpha:',alpha
         ! Get the star variable (U* & F*)        
         if ((SR-SL).eq.0) then
          write(6,*) 'Error as SR-SL=0 which will give NaN value'
          stop
         end if
         FLs=((SL*((1-alpha)*SR+alpha*vn_roe))/(SR-SL))*(qR-qL) &
            -(SL/(SR-SL))*FR+((SR)/(SR-SL))*FL
         FRs=((SR*((1-alpha)*SL+alpha*vn_roe))/(SR-SL))*(qR-qL) &
            -(SL/(SR-SL))*FR+((SR)/(SR-SL))*FL

!-------------------------------------------------------------------------------
        else if (version.eq.2) then
         ! From Linde 2002 paper(used least-square solution of jump with
         ! appropriate rescaling)
        
        else
         write(6,*) 'Incorrect version selected. Ending process...'
         stop
        end if

        if (SL.gt.0) then
         flux=FL
        else if (SL.le.0.and.SM.ge.0) then
         flux=FLs
        else if (SM.le.0.and.SR.ge.0) then
         flux=FRs
        else if (SR.lt.0) then
         flux=FR
        end if

        if (nx.eq.1) then
         flux(5)=0
        else
         write(6,*) 'Directional error'
         stop
        end if

        end subroutine HLLC_Linde


! HLLC flux from Shengtai Li--------------------------------------------
        subroutine HLLC(wL,wR,gammaa,normal,flux)
        implicit none 
        real:: gammaa
        real:: wL(8),wR(8),flux(8)
        real:: vnL,BnL,pL,caL,aL,canL,cfL,eL,HL
        real:: vnR,BnR,pR,caR,aR,canR,cfR,eR,HR
        real:: psL,psR
        real:: BxLs,BxRs,ByLs,ByRs,BzLs,BzRs,BnLs,BnRs
        real:: SL,SR,SM,ptL,ptR
        real:: HLL_B(3),HLL_v(3)
        real:: qL(8),qsL(8),FL(8),qR(8),qsR(8),FR(8)
        integer:: nx,ny,nz,normal(2)
! Declaration for Roe-averaged quantities
        real:: RT,r,u,v,w,Bx,By,Bz,H,vn
        real:: Xfactor,Yfactor,lambda(7)

        ! Allocate zero
        flux=0

        ! normal vectors        
        nx=normal(1)
        ny=normal(2)
        nz=0

        ! Left state
        vnL=wL(2)*nx+wL(3)*ny+wL(4)*nz
        BnL=wL(5)*nx+wL(6)*ny+wL(7)*nz
        caL=((wL(5)**2+wL(6)**2+wL(7)**2)/wL(1))**0.5
        aL=(gammaa*wL(8)/wL(1))**0.5
        canL=((BnL**2)/wL(1))**0.5
        cfL=(0.5*(aL**2+caL**2)+0.5*((aL**2+caL**2)**2 &
                -4*(aL**2)*canL**2)**0.5)**0.5
        eL=(wL(8)/(wL(1)*(gammaa-1))) &
         +0.5*(wL(2)**2+wL(3)**2+wL(4)**2)+ &
         +0.5*(wL(5)**2+wL(6)**2+wL(7)**2)/wL(1)
        eL=eL*wL(1)
        ! Specific enthalpy H 
        HL=(eL+wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2))/wL(1)
        qL=[wL(1),wL(1)*wL(2),wL(1)*wL(3),wL(1)*wL(4), &
            wL(5),wL(6),wL(7),eL]

        ! Right state
        vnR=wR(2)*nx+wR(3)*ny+wR(4)*nz
        BnR=wR(5)*nx+wR(6)*ny+wR(7)*nz
        caR=((wR(5)**2+wR(6)**2+wR(7)**2)/wR(1))**0.5
        aR=(gammaa*wR(8)/wR(1))**0.5
        canR=((BnR**2)/wR(1))**0.5
        cfR=(0.5*(aR**2+caR**2)+0.5*((aR**2+caR**2)**2 &
                -4*(aR**2)*canR**2)**0.5)**0.5
        eR=(wR(8)/(wR(1)*(gammaa-1))) &
         +0.5*(wR(2)**2+wR(3)**2+wR(4)**2)+ &
         +0.5*(wR(5)**2+wR(6)**2+wR(7)**2)/wR(1)
        eR=eR*wR(1)
        ! Specific enthalpy H 
        HR=(eR+wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2))/wR(1)
        qR=[wR(1),wR(1)*wR(2),wR(1)*wR(3),wR(1)*wR(4), &
            wR(5),wR(6),wR(7),eR]

        ! Try Davis speed - not recommended according to Toro!
        !SL=minval([vnL-cfL,vnR-cfR])
        !SR=maxval([vnL+cfL,vnR+cfR])
        
        ! First compute the Roe Averages
        RT=sqrt(wR(1)/wL(1))
        r=sqrt(wR(1)*wL(1))
        u=(wL(2)+RT*wR(2))/(1+RT)
        v=(wL(3)+RT*wR(3))/(1+RT)
        w=(wL(4)+RT*wR(4))/(1+RT)
        vn=u*nx+v*ny
        ! Roe average of magnetic field defined differently
        Bx=(sqrt(wR(1))*wL(5)+sqrt(wL(1))*wR(5)) &
            /(sqrt(wL(1))+sqrt(wR(1)))
        By=(sqrt(wR(1))*wL(6)+sqrt(wL(1))*wR(6)) &
            /(sqrt(wL(1))+sqrt(wR(1)))
        Bz=(sqrt(wR(1))*wL(7)+sqrt(wL(1))*wR(7)) &
            /(sqrt(wL(1))+sqrt(wR(1)))
        H=(HL+RT*HR)/(1+RT)       
        Xfactor=0.5*((wL(6)-wR(6))**2+(wL(7)-wR(7))**2)/ &
                (sqrt(wL(1))+sqrt(wR(1)))
        Yfactor=0.5*(wL(1)+wR(1))/r
        call RoeEigen(gammaa,r,u,v,w,H,Bx,By,Bz,Xfactor,Yfactor,lambda) 

        ! Wave speed estimates
        SL=minval([vnL-cfL,lambda(1)])
        SR=maxval([vnR+cfR,lambda(7)])

        SM=(wL(8)-wR(8)+wR(1)*vnR*(SR-vnR)-wL(1)*vnL*(SL-vnL)-BnL**2 &
           +BnR**2)/(wR(1)*(SR-vnR)-wL(1)*(SL-vnL))
        if (abs(SM).le.2.2204E-16) SM=0
        if (abs(SL).le.2.2204E-16) SL=0
        if (abs(SR).le.2.2204E-16) SR=0

        ! Left fluxes 
        ptL=wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2)
        FL=[wL(1)*vnL, wL(1)*vnL*wL(2)+ptL*nx-BnL*wL(5), &
           wL(1)*vnL*wL(3)+ptL*ny-BnL*wL(6), &
           wL(1)*vnL*wL(4)+ptL*nz-BnL*wL(7), &
           vnL*wL(5)-BnL*wL(2),vnL*wL(6)-BnL*wL(3),vnL*wL(7)-BnL*wL(4), &
           wL(1)*vnL*HL-BnL*(wL(2)*wL(5)+wL(3)*wL(6)+wL(4)*wL(7))]
        ! Right fluxes
        ptR=wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2)
        FR=[wR(1)*vnR, wR(1)*vnR*wR(2)+ptR*nx-BnR*wR(5), &
           wR(1)*vnR*wR(3)+ptR*ny-BnR*wR(6), &
           wR(1)*vnR*wR(4)+ptR*nz-BnR*wR(7), &
           vnR*wR(5)-BnR*wR(2),vnR*wR(6)-BnR*wR(3),vnR*wR(7)-BnR*wR(4), &
           wR(1)*vnR*HR-BnR*(wR(2)*wR(5)+wR(3)*wR(6)+wR(4)*wR(7))]

        ! Compute the HLL flux for energy quantity for HLLC energy flux
        HLL_B=(SR*wR(5:7)-SL*wL(5:7)+FR(5:7)-FL(5:7))/(SR-SL)
        HLL_v=(SR*wR(2:4)-SL*wL(2:4)+FR(2:4)-FL(2:4))/(SR-SL)

        ! Compute the HLLC flux
        if ((SL.lt.0.and.SM.ge.0).or.SL.ge.0) then
        ! This is the crucial part as this is the U*1 
          BxLs=(SR*wR(5)-SL*wL(5))/(SR-SL) !Note is same as right one
          ByLs=(SR*wR(6)-SL*wL(6)-(FR(6)-FL(6)))/(SR-SL)
          BzLs=(SR*wR(7)-SL*wL(7)-(FR(7)-FL(7)))/(SR-SL)
          psL=wL(1)*(SL-vnL)*(SM-vnL)+ptL-wL(5)**2+BxLs**2

         BnLs=BxLs*nx+ByLs*ny+BzLs*nz
         qsL=wL(1)*((SL-vnL)/(SL-SM))*[1.0, &
         SM*(nx)+wL(2)*abs(nx-1)-abs(ny)*((ByLs*BxLs-wL(6)*wL(5)) &
         /(SL-vnL))/wL(1)-abs(nz)*((BzLs*BxLs-wL(7)*wL(5)) &
         /(SL-vnL))/wL(1), &
         SM*(ny)+wL(3)*abs(ny-1)-abs(nx)*((BxLs*ByLs-wL(5)*wL(6)) &
         /(SL-vnL))/wL(1)-abs(nz)*((BzLs*ByLs-wL(7)*wL(6)) &
         /(SL-vnL))/wL(1), &
         SM*(nz)+wL(4)*abs(nz-1)-abs(nx)*((BxLs*BzLs-wL(5)*wL(7)) &
         /(SL-vnL))/wL(1)-abs(ny)*((ByLs*BzLs-wL(6)*wL(7)) &
         /(SL-vnL))/wL(1), &
         BxLs/((wL(1)*(SL-vnL))/(SL-SM)), &
         ByLs/((wL(1)*(SL-vnL))/(SL-SM)), &
         BzLs/((wL(1)*(SL-vnL))/(SL-SM)), &
         eL/wL(1) + (1/wL(1))*(((psL*SM-wL(8)*vnL)- &
         (BnLs*dot_product(HLL_B,HLL_v) &
         -BnL*dot_product(wL(5:7),wL(2:4))))/(SL-vnL))]
         if (SL.ge.0) then ! Right-going supersonic flow
          flux=FL
         else if (SL.lt.0.and.SM.ge.0) then ! Subsonic flow to the right
          flux=FL+SL*(qsL-qL)
         end if
        else if ((SM.lt.0.and.SR.ge.0).or.SR.lt.0) then
          BxRs=(SR*wR(5)-SL*wL(5))/(SR-SL)
          ByRs=(SR*wR(6)-SL*wL(6)-(FR(6)-FL(6)))/(SR-SL)
          BzRs=(SR*wR(7)-SL*wL(7)-(FR(7)-FL(7)))/(SR-SL)
          psR=wR(1)*(SR-vnR)*(SM-vnR)+ptR-wR(5)**2+BxRs**2
         
         BnRs=BxRs*nx+ByRs*ny+BzRs*nz
         qsR=wR(1)*((SR-vnR)/(SR-SM))*[1.0, &
         SM*(nx)+wR(2)*abs(nx-1)-abs(ny)*((ByRs*BxRs-wR(6)*wR(5)) &
         /(SR-vnR))/wR(1)-abs(nz)*((BzRs*BxRs-wR(7)*wR(5)) &
         /(SR-vnR))/wR(1), &
         SM*(ny)+wR(3)*abs(ny-1)-abs(nx)*((BxRs*ByRs-wR(5)*wR(6)) &
         /(SR-vnR))/wR(1)-abs(nz)*((BzRs*ByRs-wR(7)*wR(6)) &
         /(SR-vnR))/wR(1), &
         SM*(nz)+wR(4)*abs(nz-1)-abs(nx)*((BxRs*BzRs-wR(5)*wR(7)) &
         /(SR-vnR))/wR(1)-abs(ny)*((ByRs*BzRs-wR(6)*wR(7)) &
         /(SR-vnR))/wR(1), &
         BxRs/((wR(1)*(SR-vnR))/(SR-SM)), &
         ByRs/((wR(1)*(SR-vnR))/(SR-SM)), &
         BzRs/((wR(1)*(SR-vnR))/(SR-SM)), &
         eR/wR(1) + (1/wR(1))*(((psR*SM-wR(8)*vnR)- &
         (BnRs*dot_product(HLL_B,HLL_v) &
         -BnR*dot_product(wR(5:7),wR(2:4))))/(SR-vnR))]
         if (SM.lt.0.and.SR.ge.0) then ! Subsonic flow to the left
          flux=FR+SR*(qsR-qR)
         else if (SR.le.0) then ! Left-going supersonic flow
          flux=FR
         end if
        end if

        if (nx.eq.1) then
         flux(5)=0
        else
         write(6,*) 'Directional error'
         stop
        end if

        end subroutine HLLC


! HLLD flux from T.Miyoshi----------------------------------------
! Note here I used 1D fluxes
        subroutine HLLD(wL,wR,gammaa,normal,flux)
        implicit none
        real:: gammaa
        real:: wL(8),wR(8),flux(8)
        real:: vnL,BnL,pL,caL,aL,canL,cfL,eL,HL
        real:: vnR,BnR,pR,caR,aR,canR,cfR,eR,HR
        real:: psL,psR
        real:: BxLs,BxRs,ByLs,ByRs,BzLs,BzRs,BnLs,BnRs
        real:: SL,SR,SM,ptL,ptR
        real:: HLL_B(3),HLL_v(3)
        real:: qL(8),FL(8),qR(8),FR(8)
        integer:: nx,ny,nz,normal(2)
! Declare for * variables
        real:: SLs,SRs,wLs(8),wRs(8),eLs,eRs,qsL(8),qsR(8)
! Declare for ** variables 
        real:: pssL,pssR
        real:: BxLss,BxRss,ByLss,ByRss,BzLss,BzRss,BnLss,BnRss
        real:: ptLs,ptRs,ptLss,ptRss,pts
        real:: wLss(8),wRss(8),qssL(8),qssR(8),eLss,eRss
! Declare for Roe variables
        real:: RT,r,u,v,w_roe,Bx,By,Bz,H,vn,a

        ! Allocate zero
        flux=0

        ! normal vectors        
        nx=normal(1)
        ny=normal(2)
        nz=0

        ! Left state
        vnL=wL(2)*nx+wL(3)*ny+wL(4)*nz
        BnL=wL(5)*nx+wL(6)*ny+wL(7)*nz
        caL=((wL(5)**2+wL(6)**2+wL(7)**2)/wL(1))**0.5
        aL=(gammaa*wL(8)/wL(1))**0.5
        canL=((BnL**2)/wL(1))**0.5
        cfL=(0.5*(aL**2+caL**2)+0.5*((aL**2+caL**2)**2 &
               -4*(aL**2)*canL**2)**0.5)**0.5
        eL=(wL(8)/(wL(1)*(gammaa-1))) &
          +0.5*(wL(2)**2+wL(3)**2+wL(4)**2)+ &
          +0.5*(wL(5)**2+wL(6)**2+wL(7)**2)/wL(1)
        eL=eL*wL(1)
        ! Specific enthalpy H 
        HL=(eL+wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2))/wL(1)
        qL=[wL(1),wL(1)*wL(2),wL(1)*wL(3),wL(1)*wL(4), &
            wL(5),wL(6),wL(7),eL]
        ! Other required left variable
        ptL=wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2)

        ! Right state
        vnR=wR(2)*nx+wR(3)*ny+wR(4)*nz
        BnR=wR(5)*nx+wR(6)*ny+wR(7)*nz
        caR=((wR(5)**2+wR(6)**2+wR(7)**2)/wR(1))**0.5
        aR=(gammaa*wR(8)/wR(1))**0.5
        canR=((BnR**2)/wR(1))**0.5
        cfR=(0.5*(aR**2+caR**2)+0.5*((aR**2+caR**2)**2 &
                -4*(aR**2)*canR**2)**0.5)**0.5
        eR=(wR(8)/(wR(1)*(gammaa-1))) &
         +0.5*(wR(2)**2+wR(3)**2+wR(4)**2)+ &
         +0.5*(wR(5)**2+wR(6)**2+wR(7)**2)/wR(1)
        eR=eR*wR(1)
        ! Specific enthalpy H 
        HR=(eR+wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2))/wR(1)
        qR=[wR(1),wR(1)*wR(2),wR(1)*wR(3),wR(1)*wR(4), &
            wR(5),wR(6),wR(7),eR]
        ! Other required right variables
        ptR=wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2)

        ! Obtain waves speed 
        SL=minval([vnL,vnR])-maxval([cfL,cfR])
        SR=maxval([vnL,vnR])+maxval([cfL,cfR])
        SM=((SR-vnR)*wR(1)*vnR-(SL-vnL)*wL(1)*vnL-ptR+ptL)/ &
           ((SR-vnR)*wR(1)-(SL-vnL)*wL(1))
        if (abs(SM).le.2.2204E-16) SM=0
        if (abs(SL).le.2.2204E-16) SL=0
        if (abs(SR).le.2.2204E-16) SR=0
        
        ! Average total pressure in Riemann fan (from eqn (41))
        pts=((SR-vnR)*wR(1)*ptL-(SL-vnL)*wL(1)*ptR &
         +wL(1)*wR(1)*(SR-vnR)*(SL-vnL)*(vnR-vnL)) &
         /((SR-vnR)*wR(1)-(SL-vnL)*wL(1))
        ptLs=pts
        ptRs=pts
        ptLss=pts
        ptRss=pts

        ! Left fluxes 
        ptL=wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2)
        FL=[wL(1)*vnL, wL(1)*vnL*wL(2)+ptL*nx-BnL*wL(5), &
           wL(1)*vnL*wL(3)+ptL*ny-BnL*wL(6), &
           wL(1)*vnL*wL(4)+ptL*nz-BnL*wL(7), &
           vnL*wL(5)-BnL*wL(2), &
           vnL*wL(6)-BnL*wL(3),vnL*wL(7)-BnL*wL(4), &
           wL(1)*vnL*HL-BnL*(wL(2)*wL(5)+wL(3)*wL(6)+wL(4)*wL(7))]
        ! Right fluxes
        ptR=wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2)
        FR=[wR(1)*vnR, wR(1)*vnR*wR(2)+ptR*nx-BnR*wR(5), &
           wR(1)*vnR*wR(3)+ptR*ny-BnR*wR(6), &
           wR(1)*vnR*wR(4)+ptR*nz-BnR*wR(7), &
           vnR*wR(5)-BnR*wR(2),vnR*wR(6)-BnR*wR(3),vnR*wR(7)-BnR*wR(4), &
           wR(1)*vnR*HR-BnR*(wR(2)*wR(5)+wR(3)*wR(6)+wR(4)*wR(7))]

        ! Star quantities (*)
        wLs(1)=wL(1)*(SL-vnL)/(SL-SM)
         wLs(2)=SM
         wLs(5)=wL(5)
         if ((wL(1)*(SL-wL(2))*(SL-SM)-wL(5)**2).ne.0) then
          wLs(3)=wL(3)-wL(5)*wL(6)*((SM-wL(2))/(wL(1)*(SL-wL(2))*(SL-SM) &
                -wL(5)**2))
          wLs(4)=wL(4)-wL(5)*wL(7)*((SM-wL(2))/(wL(1)*(SL-wL(2))*(SL-SM) &
                -wL(5)**2))
          wLs(6)=wL(6)*((wL(1)*(SL-wL(2))**2-wL(5)**2)/(wL(1)*(SL-wL(2)) &
                *(SL-SM)-wL(5)**2))
          wLs(7)=wL(7)*((wL(1)*(SL-wL(2))**2-wL(5)**2)/(wL(1)*(SL-wL(2)) &
                *(SL-SM)-wL(5)**2))
         else
          wLs(3:4)=wL(3:4)
          wLs(6:7)=0
         end if
         eLs=(((SL-wL(2))*eL-ptL*wL(2)+pts*SM+wL(5) &
            *(dot_product(wL(2:4),wL(5:7)) &
            -dot_product(wLs(2:4),wLs(5:7))))/(SL-SM))
        

        wRs(1)=wR(1)*(SR-vnR)/(SR-SM)
         wRs(2)=SM
         wRs(5)=wR(5)
         if ((wR(1)*(SR-wR(2))*(SR-SM)-wR(5)**2).ne.0) then
          wRs(3)=wR(3)-wR(5)*wR(6)*((SM-wR(2))/(wR(1)*(SR-wR(2))*(SR-SM) &
                -wR(5)**2))
          wRs(4)=wR(4)-wR(5)*wR(7)*((SM-wR(2))/(wR(1)*(SR-wR(2))*(SR-SM) &
                -wR(5)**2))
          wRs(6)=wR(6)*((wR(1)*(SR-wR(2))**2-wR(5)**2)/(wR(1)*(SR-wR(2)) &
                *(SR-SM)-wR(5)**2))
          wRs(7)=wR(7)*((wR(1)*(SR-wR(2))**2-wR(5)**2)/(wR(1)*(SR-wR(2)) &
                *(SR-SM)-wR(5)**2))
         else
          wRs(3:4)=wR(3:4)
          wRs(6:7)=0
         end if
         eRs=(((SR-wR(2))*eR-ptR*wR(2)+pts*SM+wR(5) &
            *(dot_product(wR(2:4),wR(5:7)) &
            -dot_product(wRs(2:4),wRs(5:7))))/(SR-SM))

        ! Star star quantities (**)
        wLss(1)=wLs(1)
         wLss(2)=SM
         wLss(3)=(sqrt(wLs(1))*wLs(3)+sqrt(wRs(1))*wRs(3)+ &
          (wRs(6)-wLs(6))*sign(1.,wL(5)))/(sqrt(wLs(1))+sqrt(wRs(1)))
         wLss(4)=(sqrt(wLs(1))*wLs(4)+sqrt(wRs(1))*wRs(4)+ &
          (wRs(7)-wLs(7))*sign(1.,wL(5)))/(sqrt(wLs(1))+sqrt(wRs(1)))
         wLss(5)=wLs(5) 
         wLss(6)=(sqrt(wLs(1))*wRs(6)+sqrt(wRs(1))*wLs(6)+ &
          sqrt(wLs(1)*wRs(1))*(wRs(3)-wLs(3))*sign(1.,wL(5))) &
          /(sqrt(wLs(1))+sqrt(wRs(1)))
         wLss(7)=(sqrt(wLs(1))*wRs(7)+sqrt(wRs(1))*wLs(7)+ &
          sqrt(wLs(1)*wRs(1))*(wRs(4)-wLs(4))*sign(1.,wL(5))) &
          /(sqrt(wLs(1))+sqrt(wRs(1)))
         eLss=eLs-sqrt(wLs(1))*(dot_product(wLs(2:4),wLs(5:7))- &
         dot_product(wLss(2:4),wLss(5:7)))*sign(1.,wL(5))
         wRss=wLss
         eRss=eRs+sqrt(wRs(1))*(dot_product(wRs(2:4),wRs(5:7))- &
         dot_product(wRss(2:4),wRss(5:7)))*sign(1.,wR(5))
        
        ! Convert from primitive to conservative
        qsL=wLs
        qsL(2:4)=wLs(2:4)*wLs(1)
        qsL(8)=eLs
        qsR=wRs
        qsR(2:4)=wRs(2:4)*wRs(1)
        qsR(8)=eRs
        qssL=wLss
        qssL(2:4)=wLss(2:4)*wLss(1)
        qssL(8)=eLss
        qssR=wRss
        qssR(2:4)=wRss(2:4)*wRss(1)
        qssR(8)=eRss
        
        !Propagation of Alfven waves in intermediate states (eqn (51))
         SLs=SM-abs(wL(5))/sqrt(wLs(1))
         SRs=SM+abs(wR(5))/sqrt(wRs(1))
  
        ! Write up the fluxes
        if (SL.gt.0) then
         flux=FL 
        else if (SL.le.0.and.SLs.ge.0) then
         flux=FL+SL*qsL-SL*qL
        else if (SLs.le.0.and.SM.ge.0) then
         flux=FL+SLs*qssL-(SLs-SL)*qsL-SL*qL
        else if (SM.le.0.and.SRs.ge.0) then
         flux=FR+SRs*qssR-(SRs-SR)*qsR-SR*qR
        else if (SRs.le.0.and.SR.ge.0) then
         flux=FR+SR*qsR-SR*qR
        else if (SR.lt.0) then
         flux=FR
        end if

        if(sum(flux).ne.sum(flux)) then 
         write(6,*) 'FL',FL,'FR',FR
         write(6,*) 'qL',qL,'qR',qR
         write(6,*) 'qsL',qsL,'qsR',qsR
         write(6,*) 'qssL',qssL,'qssR',qssR
         write(6,*) SL,SR,SLs,SRs,wLs(1),wRs(1)
        end if

        if (nx.eq.1) then
         flux(5)=0
        else
         write(6,*) 'Directional error'
         stop
        end if

        end subroutine HLLD



! ROE flux-----------------------------------------------------------
        subroutine RoeS(wL,wR,gammaa,normal,flux)
        implicit none
        real:: gammaa
        real:: wL(8),wR(8),flux(8)
        real:: RT,r,u,v,w,Bx,By,Bz,H,vn
        real:: vnL,BnL,pL,caL,aL,canL,cfL,eL,HL
        real:: vnR,BnR,pR,caR,aR,canR,cfR,eR,HR
        real:: psL,psR
        real:: BxLs,BxRs,ByLs,ByRs,BzLs,BzRs,BnLs,BnRs
        real:: SL,SR,ptL,ptR
        real:: qL(8),FL(8),qR(8),FR(8)
        integer:: nx,ny,nz,normal(2),i
        real:: Xfactor,Yfactor,lambda(7)
        real :: lem(7,7), rem(7,7),a(7),coeff(7)
        ! Allocate zero
        flux=0

        ! normal vectors        
        nx=normal(1)
        ny=normal(2)
        nz=0

        ! Left state
        vnL=wL(2)*nx+wL(3)*ny+wL(4)*nz
        BnL=wL(5)*nx+wL(6)*ny+wL(7)*nz
        caL=((wL(5)**2+wL(6)**2+wL(7)**2)/wL(1))**0.5
        aL=(gammaa*wL(8)/wL(1))**0.5
        canL=((BnL**2)/wL(1))**0.5
        cfL=(0.5*(aL**2+caL**2)+0.5*((aL**2+caL**2)**2 &
                -4*(aL**2)*canL**2)**0.5)**0.5
        eL=(wL(8)/(wL(1)*(gammaa-1))) &
         +0.5*(wL(2)**2+wL(3)**2+wL(4)**2)+ &
         +0.5*(wL(5)**2+wL(6)**2+wL(7)**2)/wL(1)
        eL=eL*wL(1)
        ! Specific enthalpy H defined from google books
        HL=(eL+wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2))/wL(1)
        qL=[wL(1),wL(1)*wL(2),wL(1)*wL(3),wL(1)*wL(4), &
            wL(5),wL(6),wL(7),eL]

        ! Right state
        vnR=wR(2)*nx+wR(3)*ny+wR(4)*nz
        BnR=wR(5)*nx+wR(6)*ny+wR(7)*nz
        caR=((wR(5)**2+wR(6)**2+wR(7)**2)/wR(1))**0.5
        aR=(gammaa*wR(8)/wR(1))**0.5
        canR=((BnR**2)/wR(1))**0.5
        cfR=(0.5*(aR**2+caR**2)+0.5*((aR**2+caR**2)**2 &
                -4*(aR**2)*canR**2)**0.5)**0.5
        eR=(wR(8)/(wR(1)*(gammaa-1))) &
         +0.5*(wR(2)**2+wR(3)**2+wR(4)**2)+ &
         +0.5*(wR(5)**2+wR(6)**2+wR(7)**2)/wR(1)
        eR=eR*wR(1)
        ! Specific enthalpy H defined from google books
        HR=(eR+wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2))/wR(1)
        qR=[wR(1),wR(1)*wR(2),wR(1)*wR(3),wR(1)*wR(4), &
            wR(5),wR(6),wR(7),eR]

        ! First compute the Roe Averages
        RT=sqrt(wR(1)/wL(1))
        r=sqrt(wR(1)*wL(1))
        u=(wL(2)+RT*wR(2))/(1+RT)
        v=(wL(3)+RT*wR(3))/(1+RT)
        w=(wL(4)+RT*wR(4))/(1+RT)
        Bx=(sqrt(wR(1))*wL(5)+sqrt(wL(1))*wR(5)) &
           /(sqrt(wL(1))+sqrt(wR(1)))   
        By=(sqrt(wR(1))*wL(6)+sqrt(wL(1))*wR(6)) &
           /(sqrt(wL(1))+sqrt(wR(1)))
        Bz=(sqrt(wR(1))*wL(7)+sqrt(wL(1))*wR(7)) &
           /(sqrt(wL(1))+sqrt(wR(1)))
        H=(HL+RT*HR)/(1+RT) ! correct same as Athena
        Xfactor=0.5*((wL(6)-wR(6))**2+(wL(7)-wR(7))**2)/ &
               (sqrt(wL(1))+sqrt(wR(1)))
        Yfactor=0.5*(wL(1)+wR(1))/r
        ! Compute eigenvalues and eigematrices from Roe-averaged values
        call Eigen_cons(gammaa,r,u,v,w,H,Bx,By,Bz,Xfactor,Yfactor &
                ,lambda,rem,lem)
        ! Create intermediate states from eigenmtrices
        a=(qR(1)-qL(1))*lem(1,:)
        a=a+(qR(2)-qL(2))*lem(2,:)
        a=a+(qR(3)-qL(3))*lem(3,:)
        a=a+(qR(4)-qL(4))*lem(4,:)
        a=a+(qR(8)-qL(8))*lem(5,:)
        a=a+(qR(6)-qL(6))*lem(6,:)
        a=a+(qR(7)-qL(7))*lem(7,:)

        vn=u*nx+v*ny+w*nz

        ! Wave speed estimates
        SL=minval([vnL-cfL,lambda(1),0.])
        SR=maxval([vnR+cfR,lambda(7),0.])

        ! Left fluxes 
        ptL=wL(8)+0.5*(wL(5)**2+wL(6)**2+wL(7)**2)
        FL=[wL(1)*vnL, wL(1)*vnL*wL(2)+ptL*nx-BnL*wL(5), &
           wL(1)*vnL*wL(3)+ptL*ny-BnL*wL(6), &
           wL(1)*vnL*wL(4)+ptL*nz-BnL*wL(7), &
           vnL*wL(5)-BnL*wL(2), & 
           vnL*wL(6)-BnL*wL(3),vnL*wL(7)-BnL*wL(4), &
           wL(1)*vnL*HL-BnL*(wL(2)*wL(5)+wL(3)*wL(6)+wL(4)*wL(7))]
        ! Right fluxes
        ptR=wR(8)+0.5*(wR(5)**2+wR(6)**2+wR(7)**2)
        FR=[wR(1)*vnR, wR(1)*vnR*wR(2)+ptR*nx-BnR*wR(5), &
           wR(1)*vnR*wR(3)+ptR*ny-BnR*wR(6), &
           wR(1)*vnR*wR(4)+ptR*nz-BnR*wR(7), &
           vnR*wR(5)-BnR*wR(2),vnR*wR(6)-BnR*wR(3),vnR*wR(7)-BnR*wR(4), &
           wR(1)*vnR*HR-BnR*(wR(2)*wR(5)+wR(3)*wR(6)+wR(4)*wR(7))]

        ! Compute the ROE flux 
        flux=0.5*(FL+FR)
        do i=1,7
          coeff(i)=0.5*abs(lambda(i))*a(i)
        end do
        do i=1,7
          flux(1:4)=flux(1:4)-coeff(i)*rem(i,1:4)
          flux(6:7)=flux(6:7)-0.5*abs(lambda(i))*a(i)*rem(i,6:7)
          flux(8)=flux(8)-coeff(i)*rem(i,5)
        end do
        

        if (nx.eq.1) then
         flux(5)=0.
        else
         write(6,*) 'Directional error'
         stop
        end if

        end subroutine RoeS



        subroutine RoeEigen(gamma,d,vx,vy,vz,h,Bx,by,bz &
                ,Xfac,Yfac,lambda)
        IMPLICIT NONE
        real :: d, vx, vy, vz, h
        real :: Bx, by, bz, Xfac, Yfac
        real, dimension(7) :: lambda
        ! LOCALS
        real :: btsq, bt_starsq, bt, bt_star, gamma
        real :: vsq, vax,  vaxsq, hp, twid_asq, ct2, tsum, tdif,cf2_cs2
        real :: cfsq, cfast, cssq, cslow
        real,parameter :: eps=2.2204E-16

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

        ! eigenvalues
        lambda(1) = vx - cfast
        lambda(2) = vx - vax
        lambda(3) = vx - cslow
        lambda(4) = vx
        lambda(5) = vx + cslow
        lambda(6) = vx + vax
        lambda(7) = vx + cfast
        end subroutine
