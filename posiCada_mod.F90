        subroutine posiCada_mod(gammaa,dx,dy,warray &
                       ,wxl,wxr,wyl,wyr,Bi,lambdaPS)
        use param
        integer:: i,j,m 
        real:: warray(nx,ny,8)
        real:: dwdx(nx,ny,8),dwdy(nx,ny,8)
        real:: dwdxl(nx,ny,8),dwdxr(nx,ny,8)
        real:: dwdyl(nx,ny,8),dwdyr(nx,ny,8)
        real:: r,u,v,w,Bx,By,Bz,p
        real:: dx,dy,gammaa
        real:: glimit,alpha
        real:: drx,dux,dvx,dwx,dBx_x,dBy_x,dBz_x,dpx
        real:: dry,duy,dvy,dwy,dBx_y,dBy_y,dBz_y,dpy
        real:: cDx(8),cDy(8)
        real:: int_rhocx,int_rhocy,int_pcx,int_pcy
        real:: cDux,uDrx,cDvy,vDry
        real:: ruDux,rvDvy,rduDux,rdvDvy
        real:: wxl(nx,ny,8),wxr(nx,ny,8)
        real:: wyl(nx,ny,8),wyr(nx,ny,8)
        real:: rcx,rcy
        real:: rccx,rccy
        real:: rsx(8),rsy(8),rssx,rssy
        real:: Lx,Ly,rEx,rEy
        real:: intx,inty
        real:: du,udpx,dv,vdpy
        real:: fx,fy
        real:: lambdaPS
! Parameters - from PS.
        parameter(glimit=0.9,alpha=1.0/3.0) 
        real:: checkxl,checkxr
! From CT
        real:: Bi(nx,ny,2)

! Allocate zeros to arrays 
        cDx=0
        cDy=0
        
! Allocate zeros to variables
        cDux=0
        uDrx=0
        cDvy=0
        vDry=0
        ruDux=0
        rvDvy=0
        rduDux=0
        rdvDvy=0        
        du=0
        udpx=0
        dv=0
        vdpy=0
        intx=0
        inty=0
        
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nx,ny,warray,dwdx,dwdxl,dwdxr,gammaa, &
!$OMP dwdy,dwdyl,dwdyr,wxl,wxr,wyl,wyr,dx,dy,dt) &
!$OMP PRIVATE (i,j,m,cDx,cDy,checkxl,checkxr, &
!$OMP r,u,v,w,Bx,By,Bz,p,drx,dux,dvx,dwx,dBx_x,dBy_x,dBz_x, &
!$OMP dpx,dry,duy,dvy,dwy,dBx_y,dBy_y,dBz_y,dpy,rcx,rcy, &
!$OMP cDux,uDrx,cDvy,vDry,ruDux,rvDvy,rduDux,rdvDvy,rccx, &
!$OMP rccy,rsx,rsy,rssx,rssy,Lx,Ly,du,udpx,dv,vdpy,rEx,rEy, &
!$OMP intx,inty,fx,fy)

! Obtain dwdx & dwdy
        !$OMP DO COLLAPSE (3)
        do m=1,8
         do j=1,ny
          do i=1,nx
           dwdx(i,j,m)=wxl(i,j,m)-wxr(i,j,m)
           dwdy(i,j,m)=wyl(i,j,m)-wyr(i,j,m)
           dwdxl(i,j,m)=wxl(i,j,m)-warray(i,j,m)
           dwdxr(i,j,m)=warray(i,j,m)-wxr(i,j,m)
           dwdyl(i,j,m)=wyl(i,j,m)-warray(i,j,m)
           dwdyr(i,j,m)=warray(i,j,m)-wyr(i,j,m)        
          end do
         end do
        end do
        !$OMP END DO

!############## Start of positive scheme ########################
        !$OMP DO COLLAPSE (2)
        do i=2,nx-1
         do j=2,ny-1
! Start of nx,ny loop --------------------------------------------
          r=warray(i,j,1)
          u=warray(i,j,2)
          v=warray(i,j,3)
          w=warray(i,j,4)
          Bx=warray(i,j,5)
          By=warray(i,j,6)
          Bz=warray(i,j,7)
          p=warray(i,j,8)

! Step 2: (b) Find the flux version of DW
        ! i) For F() flux
          cDx(:)=warray(i,j,:)-0.5*(wxl(i,j,:)+wxr(i,j,:))
        ! ii) For G() flux
          cDy(:)=warray(i,j,:)-0.5*(wyl(i,j,:)+wyr(i,j,:))

        checkxl=maxval(abs(dwdxl(i,j,:)-(0.5*dwdx(i,j,:)-cDx(:))))
        checkxr=maxval(abs(dwdxr(i,j,:)-(0.5*dwdx(i,j,:)+cDx(:))))
! Step 2: (a) Limit dW        
         ! Density inequality
         if (abs(dwdx(i,j,1)).gt.(glimit*r)) then
          dwdx(i,j,1)=sign((glimit*r),dwdx(i,j,1))
          dwdxl(i,j,1)=0.5*dwdx(i,j,1)-cDx(1)
          dwdxr(i,j,1)=0.5*dwdx(i,j,1)+cDx(1)
         end if
         if (abs(dwdy(i,j,1)).gt.(glimit*r)) then
          dwdy(i,j,1)=sign((glimit*r),dwdy(i,j,1))
          dwdyl(i,j,1)=0.5*dwdy(i,j,1)-cDy(1)
          dwdyr(i,j,1)=0.5*dwdy(i,j,1)+cDy(1)
         end if
         
        ! Velocity inequality        
         if (abs(dwdx(i,j,2)).gt.glimit*dx/(dt)) then
          dwdx(i,j,2)=sign((glimit*dx/(dt)),dwdx(i,j,2))
          dwdxl(i,j,2)=0.5*dwdx(i,j,2)-cDx(2)
          dwdxr(i,j,2)=0.5*dwdx(i,j,2)+cDx(2)
         end if
         if (abs(dwdy(i,j,3)).gt.glimit*dy/(dt)) then
          dwdy(i,j,3)=sign((glimit*dy/(dt)),dwdy(i,j,3))
          dwdyl(i,j,3)=0.5*dwdy(i,j,3)-cDy(3)
          dwdyr(i,j,3)=0.5*dwdy(i,j,3)+cDy(3)
         end if        

        ! Pressure inequality
         if (abs(dwdx(i,j,8)).gt.glimit*p) then
           dwdx(i,j,8)=sign((glimit*p),dwdx(i,j,8))
           dwdxl(i,j,8)=0.5*dwdx(i,j,8)-cDx(8)
           dwdxr(i,j,8)=0.5*dwdx(i,j,8)+cDx(8)  
          end if
         if (abs(dwdy(i,j,8)).gt.glimit*p) then
           dwdy(i,j,8)=sign((glimit*p),dwdy(i,j,8))
           dwdyl(i,j,8)=0.5*dwdy(i,j,8)-cDy(8)
           dwdyr(i,j,8)=0.5*dwdy(i,j,8)+cDy(8)
          end if

        ! 1.(b) for differentiated primitive variables-rewritten for ease
        ! Note this must be after the inequality eqn 3.14
! For F() flux
        drx=dwdx(i,j,1)
        dux=dwdx(i,j,2) 
        dvx=dwdx(i,j,3)
        dwx=dwdx(i,j,4)
        dBx_x=dwdx(i,j,5) 
        dBy_x=dwdx(i,j,6)
        dBz_x=dwdx(i,j,7)
        dpx=dwdx(i,j,8)
! For G() flux
        dry=dwdy(i,j,1)
        duy=dwdy(i,j,2)       
        dvy=dwdy(i,j,3)
        dwy=dwdy(i,j,4)
        dBx_y=dwdy(i,j,5)       
        dBy_y=dwdy(i,j,6)
        dBz_y=dwdy(i,j,7)
        dpy=dwdy(i,j,8)
                
! 3 (i) For F() flux
        rcx=r-cDx(1)
! 3 (ii) For G() flux
        rcy=r-cDy(1)
! To obtain positive/negative values for eqn 3.30 & 3.31
        cDux=-0.5*(abs(dwdx(i,j,2))-dwdx(i,j,2))
        uDrx=-0.5*(abs(dwdx(i,j,1)*u)-(dwdx(i,j,1)*u))
        cDvy=-0.5*(abs(dwdy(i,j,3))-dwdy(i,j,3))
        vDry=-0.5*(abs(dwdy(i,j,1)*v)-(dwdy(i,j,1)*v))
        ruDux=0.5*(abs(dwdx(i,j,1)*dot_product([u,v,w],dwdx(i,j,2:4))) &
               +dwdx(i,j,1)*dot_product([u,v,w],dwdx(i,j,2:4)))
        rvDvy=0.5*(abs(dwdy(i,j,1)*dot_product([u,v,w],dwdy(i,j,2:4))) &
               +dwdy(i,j,1)*dot_product([u,v,w],dwdy(i,j,2:4)))
        rduDux=-0.5*(abs(dwdx(i,j,1)*dot_product(dwdx(i,j,2:4) &
               ,cDx(2:4))) &
               -dwdx(i,j,1)*dot_product(dwdx(i,j,2:4),cDx(2:4)))
        rdvDvy=-0.5*(abs(dwdy(i,j,1)*dot_product(dwdy(i,j,2:4) &
               ,cDy(2:4))) &
               -dwdy(i,j,1)*dot_product(dwdy(i,j,2:4),cDy(2:4)))
        

! Only internal cells - for eqn 3.31
        rccx=r-(-0.5*(abs(cDx(1))-cDx(1))) ! Mod
        if (rccx.lt.rcx-1e-4) then
         write(6,*) 'Error as inequality 3.31 for rccx not fulfilled!'
         stop
        end if
        rccy=r-(-0.5*(abs(cDy(1))-cDy(1)))
        if (rccy.lt.rcy-1e-4) then
         write(6,*) 'Error as inquality 3.31 for rccy not fulfilled!'
         stop
        end if

        rsx=r
        rsy=r
        
        rssx=3*r-2*rccx
        rssy=3*r-2*rccy
        if (rssx.gt.rsx(1)+1e-4) then
         write(6,*) 'Error as inequality 3.31 for rssx not fulfilled!'
         stop
        end if
        if (rssy.gt.rsy(1)+1e-4) then
         write(6,*) 'Error as inequality 3.31 for rssy not fulfilled!'
         write(6,*) rssy,rsy(1),r,rccy,3*r-2*rccy
         stop
        end if

        Lx=3*((rccx*r)/(rssx))*dot_product(cDx(2:4),cDx(2:4)) &
          +3*dot_product(cDx(5:7),cDx(5:7)) &
          +0.25*((rccx)                     &             ! Mod 4
          *dot_product(dwdx(i,j,2:4),dwdx(i,j,2:4)) &
          +dot_product(dwdx(i,j,5:7),dwdx(i,j,5:7))) &
          + ruDux-rduDux
        Ly=3*((rccy*r)/(rssy))*dot_product(cDy(2:4),cDy(2:4)) &
          +3*dot_product(cDy(5:7),cDy(5:7)) &
          +0.25*((rccy)                     &            ! Mod 4
          *dot_product(dwdy(i,j,2:4),dwdy(i,j,2:4)) &
          +dot_product(dwdy(i,j,5:7),dwdy(i,j,5:7))) &
          + rvDvy-rdvDvy

! rhoe(dW) from eqn 3.32
        du=-0.5*(abs(dwdx(i,j,2))-dwdx(i,j,2))
        udpx=-0.5*(abs(u*dwdx(i,j,8))-(u*dwdx(i,j,8)))
        dv=-0.5*(abs(dwdy(i,j,3))-dwdy(i,j,3))
        vdpy=-0.5*(abs(v*dwdy(i,j,8))-(v*dwdy(i,j,8))) 
        rEx=(p+2*(-0.5*(abs(cDx(8))-cDx(8))))/(gammaa-1)
        if (rEx.lt.0) rEx=(p+(dt/dx)*(udpx+gammaa*p*du))/(gammaa-1)
        rEy=(p+2*(-0.5*(abs(cDy(8))-cDy(8))))/(gammaa-1) 
        if (rEy.lt.0) rEy=(p+(dt/dy)*(vdpy+gammaa*p*dv))/(gammaa-1)
        
! To obtain the max value for eqn 3.33 
        intx=maxval([Lx,rEx])
        inty=maxval([Ly,rEy])
        ! For checking as reciprocal will go to zero
        if (intx.eq.0.or.inty.eq.0) then
         write(6,*) 'Potential error as intx/y is 0'
         stop
        end if

        ! for checking
        if (minval([rEx,rEy]).lt.0) then
         write(6,*) 'Negative rE will cause imaginary number.'
         write(6,*) 'rEx=',rEx,' rEy=',rEy,' p=',p,' at',i,j
        if (intx.ge.0) then
         write(6,*) 'intx not negative to cancel out -ve rEx'
         write(6,*) 'intx=',intx,' udpx=',udpx,' du=',du
         stop
        end if
        end if

! f(dW) from eqn 3.33
        fx=(rEx/intx)**0.5
        fy=(rEy/inty)**0.5
       
! Step 4: Set DW=f(dW)dW and DeltaW = f(dW)A(W)dW      
           ! 4.(i) For F() flux
           dwdx(i,j,:)=fx*dwdx(i,j,:)
           dwdxl(i,j,:)=fx*dwdxl(i,j,:)
           dwdxr(i,j,:)=fx*dwdxr(i,j,:)
           ! 4.(ii) For G() flux
           dwdy(i,j,:)=fy*dwdy(i,j,:)
           dwdyl(i,j,:)=fy*dwdyl(i,j,:)
           dwdyr(i,j,:)=fy*dwdyr(i,j,:)           
         end do
        end do
!################### EOL for positive scheme ########################
        
! Step 5: Compute Wc and Wr/l       
        !$OMP DO COLLAPSE (3)
        do m=1,8
         do j=1,ny
          do i=1,nx 
           wxl(i,j,m)=warray(i,j,m)+dwdxl(i,j,m)
           wxr(i,j,m)=warray(i,j,m)-dwdxr(i,j,m)
           wyl(i,j,m)=warray(i,j,m)+dwdyl(i,j,m)
           wyr(i,j,m)=warray(i,j,m)-dwdyr(i,j,m)
          end do
         end do
        end do
        !$OMP END DO
!$OMP END PARALLEL
! For getting dt during positive scheme implementation
        lambdaPS=maxval([maxval(wxl(:,:,2:4)),maxval(wxr(:,:,2:4)) &
                ,maxval(wyl(:,:,2:4)),maxval(wyr(:,:,2:4))])
        
        end subroutine 
