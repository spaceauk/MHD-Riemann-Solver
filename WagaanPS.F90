! Need this subroutine for K.Waagan's positive scheme
! Note that dwdx/y/z are inputs (after higher order) and outputs
        subroutine WagaanPS(int_posi,gammaa,dx,dy,dt,warray &
                   ,dwdx,dwdy,wxl,wxr,wyl,wyr,lambdaPS)
        use param
        integer:: int_posi ! 1=W, 2=p, 3=U -reconstruction
        integer:: i,j,m 
        real:: warray(nx,ny,8)
        real:: dwdx(nx,ny,8),dwdy(nx,ny,8)
        real:: r,u,v,w,Bx,By,Bz,p
        real:: dx,dy,dt,gammaa
        real:: glimit,alpha
        real:: drx,dux,dvx,dwx,dBx_x,dBy_x,dBz_x,dpx
        real:: dry,duy,dvy,dwy,dBx_y,dBy_y,dBz_y,dpy
        real:: cDx(nx,ny,8),cDy(nx,ny,8)
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

! Allocate zeros to arrays 
        cDx=0
        cDy=0
        wxl=0
        wxr=0
        wyl=0
        wyr=0
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
        
!############## Start of positive scheme ########################
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (dx,dy,dt,nx,ny,gammaa,int_posi,warray,dwdx,dwdy, &
!$OMP cDx,cDy) &
!$OMP PRIVATE (i,j,r,u,v,w,Bx,By,Bz,p,drx,dux,dvx,dwx,dpx,dBx_x, &
!$OMP dBy_x,dBz_x,dry,duy,dvy,dwy,dpy,dBx_y,dBy_y,dBz_y, &
!$OMP int_rhocx,int_rhocy,int_pcx,int_pcy,rcx,rcy,cDux,uDrx,cDvy, &
!$OMP cDry,ruDux,rvDvy,rduDux,rdvDvy,rccx,rccy,rsx,rsy,rssx,rssy, &
!$OMP Lx,Ly,du,udpx,dv,vdpy,rEx,rEy,intx,inty,fx,fy,vDry)        
        !$OMP DO         
        do i=2,nx-1
         do j=2,ny-1
! Start of nx,ny loop --------------------------------------------
          ! Define quantities for wase of coding
          r=warray(i,j,1)
          u=warray(i,j,2)
          v=warray(i,j,3)
          w=warray(i,j,4)
          Bx=warray(i,j,5)
          By=warray(i,j,6)
          Bz=warray(i,j,7)
          p=warray(i,j,8)
          
! Step 2: (a) Limit dW        
         ! Density inequality
         if (abs(dwdx(i,j,1)).gt.(glimit*0.5*r)) then
          dwdx(i,j,1)=sign((glimit*0.5*r),dwdx(i,j,1))
         end if
         if (abs(dwdy(i,j,1)).gt.(glimit*0.5*r)) then
          dwdy(i,j,1)=sign((glimit*0.5*r),dwdy(i,j,1))
         end if
        
        ! Velocity inequality
         if (int_posi.eq.1.or.int_posi.eq.2) then ! W or p reconst.
          if (abs(dwdx(i,j,2)).gt.glimit*dx/(dt*(gammaa+1))) then
           dwdx(i,j,2)=sign((glimit*dx/(dt*(gammaa+1))),dwdx(i,j,2))
          end if
         if (abs(dwdy(i,j,3)).gt.glimit*dy/(dt*(gammaa+1))) then
           dwdy(i,j,3)=sign((glimit*dy/(dt*(gammaa+1))),dwdy(i,j,3))
          end if
        
         else if (int_posi.eq.3) then ! U-reconst.
          if (abs(dwdx(i,j,2)).gt.glimit*dx/(dt*2)) then        
           dwdx(i,j,2)=sign((glimit*dx/(dt*2)),dwdx(i,j,2))
          end if
          if (abs(dwdy(i,j,3)).gt.glimit*dy/(dt*2)) then
           dwdy(i,j,3)=sign((glimit*dy/(dt*2)),dwdy(i,j,3))
          end if
         end if         
        
        ! Pressure inequality
        if (int_posi.eq.1.or.int_posi.eq.2) then ! W or p reconst.
         if (abs(dwdx(i,j,8)).gt.glimit*p/(gammaa+1)) then
           dwdx(i,j,8)=sign((glimit*p/(gammaa+1)),dwdx(i,j,8))
          end if
         if (abs(dwdy(i,j,8)).gt.glimit*p/(gammaa+1)) then
           dwdy(i,j,8)=sign((glimit*p/(gammaa+1)),dwdy(i,j,8))
          end if

        else if (int_posi.eq.3) then ! U reconst.
         ! Need pc & rhoc for this
         drx=dwdx(i,j,1)
         dux=dwdx(i,j,2)
         dpx=dwdx(i,j,8)
         dry=dwdy(i,j,1)
         dvy=dwdy(i,j,3)
         dpy=dwdy(i,j,8)
         cDx(i,j,1)=(dt/(2*dx))*(u*drx + r*dux)
         cDx(i,j,8)=(dt/(2*dx))*(u*dpx+gammaa*p*dux)
         cDy(i,j,1)=(dt/(2*dy))*(v*dry + r*dvy)
         cDy(i,j,8)=(dt/(2*dy))*(v*dpy+gammaa*p*dvy)
         int_rhocx = r-cDx(i,j,1)
         int_rhocy = r-cDy(i,j,1)
         int_pcx = p-cDx(i,j,8)
         int_pcy = p-cDy(i,j,8)
         if (abs(dwdx(i,j,8)).gt. &
            ( glimit*(2*int_pcx-(0.25*(gammaa-1))* &
             (int_rhocx*norm2(dwdx(i,j,2:4))**2 &
             +norm2(dwdx(i,j,5:7))**2)) ) ) then 
          dwdx(i,j,8)= sign(( glimit*(2*int_pcx-(0.25*(gammaa-1))* &
             (int_rhocx*norm2(dwdx(i,j,2:4))**2 &
             +norm2(dwdx(i,j,5:7))**2)) ), dwdx(i,j,8))
         end if
         if (abs(dwdy(i,j,8)).gt. &
            ( glimit*(2*int_pcy-(0.25*(gammaa-1))* &
             (int_rhocy*norm2(dwdy(i,j,2:4))**2 &
             +norm2(dwdy(i,j,5:7))**2)) ) ) then
          dwdy(i,j,8)= sign(( glimit*(2*int_pcy-(0.25*(gammaa-1))* &
             (int_rhocy*norm2(dwdy(i,j,2:4))**2 &
             +norm2(dwdy(i,j,5:7))**2)) ), dwdy(i,j,8))
         end if
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
        
! Step 2: (b) Find the flux version of DW
        ! i) For F() flux
        cDx(i,j,1)=(dt/(2*dx))*(u*drx + r*dux)
        cDx(i,j,2)=(dt/(2*dx))*(u*dux+(dpx &
                +dot_product([By,Bz],[dBy_x,dBz_x])+Bx*dBx_x)/r)
        cDx(i,j,3)=(dt/(2*dx))*(u*dvx-((Bx*dBy_x+By*dBx_x)/r))
        cDx(i,j,4)=(dt/(2*dx))*(u*dwx-((Bx*dBz_x+Bz*dBx_x)/r))
        cDx(i,j,5)=(dt/(2*dx))*u*dBx_x
        cDx(i,j,6)=(dt/(2*dx))*(u*dBy_x+dux*By-Bx*dvx)
        cDx(i,j,7)=(dt/(2*dx))*(u*dBz_x+dux*Bz-Bx*dwx)
        cDx(i,j,8)=(dt/(2*dx))*(u*dpx+gammaa*p*dux)
        ! ii) For G() flux
        cDy(i,j,1)=(dt/(2*dy))*(v*dry + r*dvy)
        cDy(i,j,2)=(dt/(2*dy))*(v*duy-((By*dBx_y+Bx*dBy_y)/r))
        cDy(i,j,3)=(dt/(2*dy))*(v*dvy+(dpy &
                +dot_product([Bx,Bz],[dBx_y,dBz_y])+By*dBy_y)/r)
        cDy(i,j,4)=(dt/(2*dy))*(v*dwy-((By*dBz_y+Bz*dBy_y)/r))
        cDy(i,j,5)=(dt/(2*dy))*(v*dBx_y+dvy*Bx-By*duy)
        cDy(i,j,6)=(dt/(2*dy))*v*dBy_y
        cDy(i,j,7)=(dt/(2*dy))*(v*dBz_y+dvy*Bz-By*dwy)
        cDy(i,j,8)=(dt/(2*dy))*(v*dpy+gammaa*p*dvy)
        
! 3 (i) For F() flux
        rcx=r-cDx(i,j,1)
! 3 (ii) For G() flux
        rcy=r-cDy(i,j,1)
! To obtain positive/negative values for eqn 3.30 & 3.31
        cDux=-0.5*(abs(dwdx(i,j,2))-dwdx(i,j,2))
        uDrx=-0.5*(abs(dwdx(i,j,1)*u)-(dwdx(i,j,1)*u))
        cDvy=-0.5*(abs(dwdy(i,j,3))-dwdy(i,j,3))
        vDry=-0.5*(abs(dwdy(i,j,1)*v)-(dwdy(i,j,1)*v))
        if (int_posi.eq.1) then
         ruDux=0.5*(abs(dwdx(i,j,1)*dot_product([u,v,w],dwdx(i,j,2:4))) &
               +dwdx(i,j,1)*dot_product([u,v,w],dwdx(i,j,2:4)))
         rvDvy=0.5*(abs(dwdy(i,j,1)*dot_product([u,v,w],dwdy(i,j,2:4))) &
               +dwdy(i,j,1)*dot_product([u,v,w],dwdy(i,j,2:4)))
         rduDux=-0.5*(abs(dwdx(i,j,1)*dot_product(dwdx(i,j,2:4) &
               ,cDx(i,j,2:4))) &
               -dwdx(i,j,1)*dot_product(dwdx(i,j,2:4),cDx(i,j,2:4)))
         rdvDvy=-0.5*(abs(dwdy(i,j,1)*dot_product(dwdy(i,j,2:4) &
               ,cDy(i,j,2:4))) &
               -dwdy(i,j,1)*dot_product(dwdy(i,j,2:4),cDy(i,j,2:4)))
        end if
        

! Only internal cells - for eqn 3.31
        rccx=r-(0.5*dt/dx)*(r*cDux+uDrx)
        if (rccx.lt.rcx-1e-4) then
         write(6,*) 'Error as inequality 3.31 for rccx not fulfilled!'
         stop
        end if
        rccy=r-(0.5*dt/dy)*(r*cDvy+vDry)
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

        if (int_posi.eq.1) then ! for W-reconstruction
         Lx=3*((rccx*r)/(rssx))*dot_product(cDx(i,j,2:4),cDx(i,j,2:4)) &
          +3*dot_product(cDx(i,j,5:7),cDx(i,j,5:7)) &
          +0.25*((rccx+dwdx(i,j,1)**2/(2*rssx)) &
          *dot_product(dwdx(i,j,2:4),dwdx(i,j,2:4)) &
          +dot_product(dwdx(i,j,5:7),dwdx(i,j,5:7))) &
          + 0.5*ruDux-(rccx/rssx+1)*rduDux
         Ly=3*((rccy*r)/(rssy))*dot_product(cDy(i,j,2:4),cDy(i,j,2:4)) &
          +3*dot_product(cDy(i,j,5:7),cDy(i,j,5:7)) &
          +0.25*((rccy+dwdy(i,j,1)**2/(2*rssy)) &
          *dot_product(dwdy(i,j,2:4),dwdy(i,j,2:4)) &
          +dot_product(dwdy(i,j,5:7),dwdy(i,j,5:7))) &
          + 0.5*rvDvy-(rccy/rssy+1)*rdvDvy
        else if (int_posi.eq.2) then ! for p-reconstruction
         Lx=3*((rccx*r)/(rssx))*dot_product(cDx(i,j,2:4),cDx(i,j,2:4)) &
          +3*dot_product(cDx(i,j,5:7),cDx(i,j,5:7)) &
          +0.25*(rccx*dot_product(dwdx(i,j,2:4),dwdx(i,j,2:4))  &
          +dot_product(dwdx(i,j,5:7),dwdx(i,j,5:7)))
         Ly=3*((rccy*r)/(rssy))*dot_product(cDy(i,j,2:4),cDy(i,j,2:4)) &
          +3*dot_product(cDy(i,j,5:7),cDy(i,j,5:7))  &
          +0.25*(rccy*dot_product(dwdy(i,j,2:4),dwdy(i,j,2:4))  &
          +dot_product(dwdy(i,j,5:7),dwdy(i,j,5:7)))
        else if (int_posi.eq.3) then ! for U-reconstruction
         Lx=3*((rccx*r)/(rssx))*dot_product(cDx(i,j,2:4),cDx(i,j,2:4)) &
          +3*dot_product(cDx(i,j,5:7),cDx(i,j,5:7))
         Ly=3*((rccy*r)/(rssy))*dot_product(cDy(i,j,2:4),cDy(i,j,2:4)) &
          +3*dot_product(cDy(i,j,5:7),cDy(i,j,5:7))
        end if

! rhoe(dW) from eqn 3.32
        du=-0.5*(abs(dwdx(i,j,2))-dwdx(i,j,2))
        udpx=-0.5*(abs(u*dwdx(i,j,8))-(u*dwdx(i,j,8)))
        dv=-0.5*(abs(dwdy(i,j,3))-dwdy(i,j,3))
        vdpy=-0.5*(abs(v*dwdy(i,j,8))-(v*dwdy(i,j,8)))
        rEx=(p+(dt/dx)*(udpx+gammaa*p*du))/(gammaa-1)
        rEy=(p+(dt/dy)*(vdpy+gammaa*p*dv))/(gammaa-1)
        
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
         write(6,*) 'rEx=',rEx,' rEy=',rEy,' p=',p
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
           cDx(i,j,:)=fx*cDx(i,j,:) 
           ! 4.(ii) For G() flux
           dwdy(i,j,:)=fy*dwdy(i,j,:)
           cDy(i,j,:)=fy*cDy(i,j,:)           
         end do
        end do
!$OMP END PARALLEL
!################### EOL for positive scheme ########################

! Step 5: Compute Wc and Wr/l
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nx,ny,wxl,wxr,warray,dwdx,cDx,wyl,wyr,dwdy,cDy) &        
!$OMP PRIVATE (i,j,m)
        !$OMP DO COLLAPSE (3)
        do m=1,8
         do j=2,ny-1
          do i=3,nx-1
           wxl(i-1,j,m)=warray(i-1,j,m)+0.5*dwdx(i-1,j,m)-cDx(i-1,j,m)
           wxr(i,j,m)=warray(i,j,m)-0.5*dwdx(i,j,m)-cDx(i,j,m)
          end do
         end do
        end do
        !$OMP END DO
        !$OMP DO COLLAPSE (3)
        do m=1,8
         do j=3,ny-1
          do i=2,nx-1
           wyl(i,j-1,m)=warray(i,j-1,m)+0.5*dwdy(i,j-1,m)-cDy(i,j-1,m)
           wyr(i,j,m)=warray(i,j,m)-0.5*dwdy(i,j,m)-cDy(i,j,m)
          end do
         end do
        end do
        !$OMP END DO
!$OMP END PARALLEL
        ! Account for B.C.
        wxr(nx-1,2:ny-1,:)=warray(nx-1,2:ny-1,:) &
               +0.5*dwdx(nx-1,2:ny-1,:)-cDx(nx-1,2:ny-1,:)
        wxl(nx-1,2:ny-1,:)=wxr(nx-1,2:ny-1,:)
        wyr(2:nx-1,ny-1,:)=warray(2:nx-1,ny-1,:) &
               +0.5*dwdy(2:nx-1,ny-1,:)-cDy(2:nx-1,ny-1,:)
        wyl(2:nx-1,ny-1,:)=wyr(2:nx-1,ny-1,:)
        wxr(2,2:ny-1,:)=warray(2,2:ny-1,:) &
               -0.5*dwdx(2,2:ny-1,:)-cDx(2,2:ny-1,:)
        wxl(2,2:ny-1,:)=wxr(2,2:ny-1,:)
        wyr(2:nx-1,2,:)=warray(2:nx-1,2,:) &
               -0.5*dwdy(2:nx-1,2,:)-cDy(2:nx-1,2,:)
        wyl(2:nx-1,2,:)=wyr(2:nx-1,2,:)
! For getting dt during positive scheme implementation
        lambdaPS=maxval([maxval(wxl(:,:,2:4)),maxval(wxr(:,:,2:4)) &
         ,maxval(wyl(:,:,2:4)),maxval(wyr(:,:,2:4))])
        
        end subroutine 
