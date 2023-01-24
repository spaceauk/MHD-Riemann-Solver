!###################################################################
!#                 Constrained transport function
!##################################################################
! Where F is for flux along the x-direction, G is for flux along the
! y-direction 
        subroutine CT2D(kstep,fluxX,fluxY,Bi,dt,dx,dy, &
                  Binew,q,energyc) ! outputs
        use omp_lib
        use param
        implicit none
        integer:: i,j,kstep,energyc
        real:: dx,dy,dt
        real:: q(nx,ny,8),Bi(nx,ny,2)
        real:: Binew(nx,ny,2)
        real:: EMFz(nx,ny),EMFcc(nx,ny)
        real:: fluxX(nx,ny,8)
        real:: fluxY(nx,ny,8)
        real:: divB,diffBx,diffBy,divBold
        real:: Bfcx(nx,ny),Bfcy(nx,ny)
        integer::check(2)
! Allocate zeros first
        divB=0
        divBold=0
        EMFz=0
        EMFcc=0
        diffBx=0
        diffBy=0
       
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

! Ensure BC
        Bi(1,:,:)=Bi(2,:,:)
        Bi(:,1,:)=Bi(:,2,:)
        Bi(nx,:,:)=Bi(nx-1,:,:)
        Bi(:,ny,:)=Bi(:,ny-1,:)

! Compute EMF at zone center
!$OMP PARALLEL DEFAULT(NONE) SHARED(nx,ny,EMFcc,q) PRIVATE(i,j)
        !$OMP DO
        do j=1,ny
         do i=1,nx       
        EMFcc(i,j)=(q(i,j,5)*q(i,j,3)-q(i,j,6)*q(i,j,2))/q(i,j,1)
         end do
        end do
!$OMP END PARALLEL
 
! main loop
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP PRIVATE (i,j) &
!$OMP SHARED (nx,ny,EMFz,fluxX,fluxY)
        !$OMP DO COLLAPSE (2)
        do j=2,ny   
         do i=2,nx
          EMFz(i,j)=0.25*(fluxY(i,j,5)+fluxY(i-1,j,5)-fluxX(i,j,6) &
                  -fluxX(i,j-1,6))
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

        if (divB.gt.0.01) then
         write(6,*) 'CT not working properly as divB=',divB
         stop
        end if
       
! Correct magnetic field & energy
        !call check_2Ddata(nx,ny,Bfcy(:,:))
        !stop
        if (energyc.eq.1) q(:,:,8)=q(:,:,8) & ! energy correction
           +0.5*((Bfcx**2+Bfcy**2)-(q(:,:,5)**2+q(:,:,6)**2))
        q(:,:,5)=Bfcx
        q(:,:,6)=Bfcy

        end subroutine       
