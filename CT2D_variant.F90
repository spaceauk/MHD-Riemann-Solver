!###################################################################
!#              Constrained transport variant function
!##################################################################
! Where F is for flux along the x-direction, G is for flux along the
! y-direction 
        subroutine CT2D_variant(kstep,fluxX,fluxY,Bi,dt,dx,dy, &
                   Binew,q,energyc) ! outputs
        use omp_lib
        use param
        implicit none
        integer:: i,j,kstep,energyc
        real:: dx,dy,dt
        real:: q(nx,ny,8),Bi(nx,ny,2)
        real:: Binew(nx,ny,2)
        real:: EMFz(nx,ny)
        real:: fluxX(nx,ny,8)
        real:: fluxY(nx,ny,8)
        real:: divB,diffBx,diffBy,divBold
        real:: Bfcx(nx,ny),Bfcy(nx,ny)
        real:: p(nx,ny),C(nx,ny)
        real:: gammaa,psi,dpx,dpy,divV,beta,delta
        integer::check(2)
        logical:: sw1,sw2
        real:: tempC,tempP
! Allocate zeros first
        divB=0
        divBold=0
        EMFz=0
        diffBx=0
        diffBy=0
        C=0
        p=0
       
! Apply BC to flux (at Ghost cells) as it will affect the E since it
! depends on +1 too
        fluxX(1,:,:)=fluxX(2,:,:)
        fluxX(nx,:,:)=fluxX(nx-1,:,:)
        fluxY(1,:,:)=fluxY(2,:,:)
        fluxY(nx,:,:)=fluxY(nx-1,:,:)
        fluxX(:,1,:)=fluxX(:,2,:)
        fluxX(:,ny,:)=fluxX(:,ny-1,:)
        fluxY(:,1,:)=fluxY(:,2,:)
        fluxY(:,ny,:)=fluxY(:,ny-1,:)

!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nx,ny,gammaa,p,q,C) &
!$OMP PRIVATE (i,j)
        !$OMP DO COLLAPSE (2)
        do j=1,ny
         do i=1,nx        
          p(i,j)=(gammaa-1)*(q(i,j,8)-(0.5*(q(i,j,2)**2+q(i,j,3)**2+ &
            q(i,j,4)**2)/q(i,j,1))-0.5*(q(i,j,5)**2+ &
            q(i,j,6)**2+q(i,j,7)**2))
          C(i,j)=sqrt((gammaa*p(i,j)+(q(i,j,5)**2+q(i,j,6)**2) &
                +q(i,j,7)**2)/q(i,j,1))
         end do
        end do
!$OMP END PARALLEL

! Ensure BC
        Bi(1,:,:)=Bi(2,:,:)
        Bi(:,1,:)=Bi(:,2,:)
        Bi(nx,:,:)=Bi(nx-1,:,:)
        Bi(:,ny,:)=Bi(:,ny-1,:)
 
! main loop 
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP PRIVATE (i,j,beta,delta,dpx,dpy,divV,sw1,sw2,psi,tempP,tempC) &
!$OMP SHARED (nx,ny,p,q,EMFz,fluxX,fluxY,C) 
        !$OMP DO COLLAPSE (2)
        do j=2,ny 
         do i=2,nx         
           beta=0.5
           delta=0.1

           dpx=0.5*abs(p(i,j)+p(i,j-1)-p(i-1,j)-p(i-1,j-1))
           dpy=0.5*abs(p(i,j)+p(i-1,j)-p(i,j-1)-p(i-1,j-1))
           divV=0.5*(q(i,j,2)/q(i,j,1)+q(i,j-1,2)/q(i,j-1,1) &
               -q(i-1,j,2)/q(i-1,j,1)-q(i-1,j-1,2)/q(i-1,j-1,1)) &
               +0.5*(q(i,j,3)/q(i,j,1)+q(i-1,j,3)/q(i-1,j,1) &
               -q(i,j-1,3)/q(i,j-1,1)-q(i-1,j-1,3)/q(i-1,j-1,1))

           ! Switch 1: pick out situations in vicinity of strong
           ! magnetosonic shock 
           tempP=minval([p(i-1,j-1),p(i,j-1),p(i-1,j),p(i,j)])
           if ((dPx+dPy).gt.beta*tempP) then
            sw1=.true.
           else
            sw1=.false.
           end if
           ! Switch 2: pick out strongly compressive motions
           tempC=minval([C(i-1,j-1),C(i,j-1),C(i-1,j),C(i,j)])
           if ((-delta*tempC).gt.divV) then
            sw2=.true.
           else
            sw2=.false.
           end if
           ! Switched on sw1 & sw2 means that flow in local region has a
           ! shock 
           if (sw1.and.sw2) then
            psi=dPx/(dPx+dPy)
           else
            psi=0.5
           end if

           EMFz(i,j)=0.5*(1-psi)*(fluxY(i,j,5)+fluxY(i-1,j,5)) &
                    +0.5*psi*(-fluxX(i,j,6)-fluxX(i,j-1,6))
         end do
        end do
!$OMP END PARALLEL
        EMFz(1,:)=EMFz(2,:)
        EMFz(:,1)=EMFz(:,2)
        EMFz(nx,:)=EMFz(nx-1,:)
        EMFz(:,ny)=EMFz(:,ny-1)
         

        if (kstep.eq.1) then
        Binew(:,:,1)=Bi(:,:,1)-0.5*(dt/dy)* &
              (EOSHIFT(EMFz(:,:),DIM=2,SHIFT=+1) - EMFz(:,:))
        Binew(:,:,2)=Bi(:,:,2)+0.5*(dt/dx)* &
              (EOSHIFT(EMFz(:,:),DIM=1,SHIFT=+1) - EMFz(:,:))
        else if (kstep.eq.2) then
        Binew(:,:,1)=Bi(:,:,1)-(dt/dy)* &
               (EOSHIFT(EMFz(:,:),DIM=2,SHIFT=+1) - EMFz(:,:))
        Binew(:,:,2)=Bi(:,:,2)+(dt/dx)* &
               (EOSHIFT(EMFz(:,:),DIM=1,SHIFT=+1) - EMFz(:,:))
        end if
        Binew(1,:,:)=Binew(2,:,:)
        Binew(:,1,:)=Binew(:,2,:)
        Binew(nx,:,:)=Binew(nx-1,:,:)
        Binew(:,ny,:)=Binew(:,ny-1,:)

! Check whether magnetic field is solenoidal or not
        check=1
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP PRIVATE (i,j,divBold,check) &
!$OMP SHARED (nx,ny,dx,dy,Binew,Bfcx,Bfcy) &
!$OMP REDUCTION (max: divB,diffBx,diffBy)      
        !$OMP DO COLLAPSE (2)
        do j=2,ny-1
         do i=2,nx-1
           divBold=divB
           divB=maxval([divB,abs((Binew(i+1,j,1)-Binew(i,j,1))/dx &
               +(Binew(i,j+1,2)-Binew(i,j,2))/dy)]) 
           if (divB.gt.divBold.and.divB.gt.1e-6) then
                write(1000,*) 'CT error:',divB, &
                (Binew(i+1,j,1)-Binew(i,j,1)), &
                (Binew(i,j+1,2)-Binew(i,j,2)),i,j
                divBold=divB
           end if 
           if (abs(Binew(i+1,j,1)-Binew(i,j,1)).gt.diffBx) check=[i,j]
           diffBx=maxval([diffBx,abs(Binew(i+1,j,1)-Binew(i,j,1))])
           diffBy=maxval([diffBy,abs(Binew(i,j+1,2)-Binew(i,j,2))])
           Bfcx(i,j)=0.5*(Binew(i+1,j,1)+Binew(i,j,1))
           Bfcy(i,j)=0.5*(Binew(i,j+1,2)+Binew(i,j,2))
         end do
        end do
!$OMP END PARALLEL
        Bfcx(1,:)=Bfcx(2,:)
        Bfcx(:,1)=Bfcx(:,2)
        Bfcx(nx,:)=Bfcx(nx-1,:)
        Bfcx(:,ny)=Bfcx(:,ny-1)
        Bfcy(1,:)=Bfcy(2,:)
        Bfcy(:,1)=Bfcy(:,2)
        Bfcy(nx,:)=Bfcy(nx-1,:)
        Bfcy(:,ny)=Bfcy(:,ny-1)

        write(1000,*) divB,diffBx,diffBy,check
        if (divB.gt.0.01) then
         write(6,*) 'CT not working properly as divB=',divB
         stop
        end if
       
! Correct magnetic field & energy
        if (energyc.eq.1) q(:,:,8)=q(:,:,8) & ! energy correction
            +0.5*((Bfcx**2+Bfcy**2)-(q(:,:,5)**2+q(:,:,6)**2))
        write(666,*) 'Max diff for Bx=',maxval(abs(q(:,:,5)-Bfcx)), &
                ' By=',maxval(abs(q(:,:,6)-Bfcy))
        q(:,:,5)=Bfcx
        q(:,:,6)=Bfcy

        end subroutine       
