! Compact third-order limiter functions for finite volume methods
! From Miroslav Cada *, Manuel Torrilhon
        subroutine Cada3rd(dx,dy,w,wxL,wxR,wyL,wyR)
        use param
        integer:: i,j,m
        real:: dx,dy
        real:: w(nx,ny,8)
        real:: wxL(nx,ny,8),wxR(nx,ny,8)
        real:: wyL(nx,ny,8),wyR(nx,ny,8)
        real:: delta(nx,ny),phiL(nx,ny),phiR(nx,ny)
        real, parameter:: r=1.  ! Range of values possible
        real, parameter:: eps=2.2204E-16
        real:: theta,eta,phi_hat
        
        do m=1,8
! For x-------------------------------------------------------------------
         delta=EOSHIFT(w(:,:,m),DIM=1,shift=+1)-w(:,:,m)
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nx,ny,dx,dy,m,w,wxL,wxR,wyL,wyR,delta,phiL,phiR) &
!$OMP PRIVATE (i,j,eta,theta)
        !$OMP DO COLLAPSE (2)
         do j=2,ny-1
          do i=2,nx-1
           theta=1. 
           if (delta(i,j).ne.0) theta=delta(i-1,j)/delta(i,j)           
           eta=(delta(i-1,j)**2+delta(i,j)**2)/((r*dx)**2)
           phiL(i,j)=phifunc(theta,eta,eps)
           theta=1.
           if (delta(i-1,j).ne.0) theta=delta(i,j)/delta(i-1,j)
           phiR(i,j)=phifunc(theta,eta,eps)
          end do
         end do
        !$OMP END DO
        !$OMP DO COLLAPSE (2)
         do j=3,ny-1
          do i=3,nx-1
           wxL(i-1,j,m)=w(i-1,j,m)+0.5*phiL(i-1,j)*delta(i-1,j)
           wxR(i,j,m)=w(i,j,m)-0.5*phiR(i,j)*delta(i-1,j) ! Note is -1 here!
          end do
         end do
        !$OMP END DO
         wxL(nx-1,2:ny-1,m)=wxR(nx-1,2:ny-1,m)
         wxR(2,2:ny-1,m)=wxL(2,2:ny-1,m)

! For y------------------------------------------------------------------------------
        delta=0.
        delta=EOSHIFT(w(:,:,m),DIM=2,shift=+1)-w(:,:,m)
        !$OMP DO COLLAPSE (2)
         do j=2,ny-1
          do i=2,nx-1
           theta=1.
           if (delta(i,j).ne.0) theta=delta(i,j-1)/delta(i,j)
           eta=(delta(i,j-1)**2+delta(i,j)**2)/((r*dy)**2)
           phiL(i,j)=phifunc(theta,eta,eps)
           theta=1.
           if (delta(i,j-1).ne.0) theta=delta(i,j)/delta(i,j-1)
           phiR(i,j)=phifunc(theta,eta,eps)
          end do
         end do
        !$OMP END DO
        !$OMP DO COLLAPSE (2)
         do j=3,ny-1
          do i=2,nx-1
           wyL(i,j-1,m)=w(i,j-1,m)+0.5*phiL(i,j-1)*delta(i,j-1)
           wyR(i,j,m)=w(i,j,m)-0.5*phiR(i,j)*delta(i,j-1)
          end do
         end do
        !$OMP END DO
!$OMP END PARALLEL
         wyL(2:nx-1,ny-1,m)=wyR(2:nx-1,ny-1,m)
         wyR(2:nx-1,2,m)=wyL(2:nx-1,2,m)

        end do
        
        end subroutine

! Obtain phi-------------------------------
        function phifunc(theta,eta,eps)       
        real:: theta
        real:: eta,phi_hat,phifunc
        
        phi_hat=max(-0.5*theta,minval([2.*theta,(2.+theta)/3.,1.6]))
        phi_hat=max(0.,min((2+theta)/3.,phi_hat))
        if (eta.le.1-eps) then
         phifunc=(2.+theta)/3.
        else if (eta.ge.1+eps) then
         phifunc=phi_hat
        else
         phifunc=0.5*((1.-(eta-1.)/eps)*(2.+theta)/3. + &
                (1+(eta-1.)/eps)*phi_hat)
        end if
        
        end 
