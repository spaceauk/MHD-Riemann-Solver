!########################################################################
!#   MUSCL Monotonic Upstreat Centered Scheme for Conservation Laws
!#   Van Leer's MUSCL reconstruction scheme using piece wise linear
!#   reconstruction
!#
!#   Flux at j+1/2
!# 
!#     j+1/2         Cell's grid:
!#   | wL|   |
!#   |  /|wR |           1   2   3   4        N-2 N-1  N
!#   | / |\  |   {x=0} |-o-|-o-|-o-|-o-| ... |-o-|-o-|-o-| {x=L}
!#   |/  | \ |         1   2   3   4   5        N-1  N  N+1
!#   |   |  \|
!#   |   |   |       NC: Here cells 1 and N are ghost cells
!#     j  j+1            faces 2 and N, are the real boundary faces.
!#######################################################################
        subroutine MUSCL2D(dx,dy,dt,gammaa,q,limiter,fluxMth &
         ,res & ! The output for this is residual (res).
         ,init_ct,Bi,fluxX,fluxY & ! for CT method
         ,diff,lambdaPS)  ! For yhPS/CT stuffs

        use omp_lib
        use param
        use timers
        implicit none
        integer:: i,j,m
        character(len = 15):: limiter,fluxMth
        real:: dx,dy,dt,gammaa
        real:: dww,dwe,dws,dwn,dwf,dwb,dwc
        real:: q(nx,ny,8) 
        real:: wxL(nx,ny,8),wxR(nx,ny,8)  ! for cell edge values
        real:: wyL(nx,ny,8),wyR(nx,ny,8)
        real:: w(nx,ny,8),flux(8),res(nx,ny,8)
        real:: dwdx(nx,ny,8),dwdy(nx,ny,8)
        real:: lambdaPS ! from positive scheme (can choose to use or not)
! Arrays for constrained transport
        integer:: init_ct
        real:: fluxX(nx,ny,8)
        real:: fluxY(nx,ny,8)
        real:: Bi(nx,ny,2) ! Both inputs and outputs
! For yhPS/CT stuffs
        real:: diff(nx,ny,2)
! For OpenMP implementation
        real:: cput0,cput1,ompt0,ompt1 

        ompt0=omp_get_wtime()        
! q is conserved quantity while w is primitive quantity
        w=0
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nx,ny,gammaa,w,q) &
!$OMP PRIVATE (i,j)
        !$OMP DO COLLAPSE (2)                        
        do j=1,ny
         do i=1,nx
          w(i,j,1)=q(i,j,1)
          w(i,j,2)=q(i,j,2)/q(i,j,1)
          w(i,j,3)=q(i,j,3)/q(i,j,1)
          w(i,j,4)=q(i,j,4)/q(i,j,1)
          w(i,j,5:7)=q(i,j,5:7)
          w(i,j,8)=(gammaa-1)*(q(i,j,8)-0.5*w(i,j,1)*(w(i,j,2)**2 &
           +w(i,j,3)**2+w(i,j,4)**2)-0.5*(w(i,j,5)**2+w(i,j,6)**2 &
           +w(i,j,7)**2)) ! gives you the pressure
         end do
        end do
!$OMP END PARALLEL

! NEED to allocate zero otherwise will slow down process
        dwdx=0
        dwdy=0
        flux=0
        res=0
        wxL=0
        wxR=0
        wyL=0
        wyR=0
        res=0
        fluxX=0
        fluxY=0

! Slope Limiter -
        ompt1=omp_get_wtime() 
        call slopelimiter(limiter,8,w,dwdx,dwdy)
        t_MUSCL(1)=t_MUSCL(1)+(omp_get_wtime()-ompt1)
        
! Positive scheme: placed in separate file cause need a lot inputs
        ompt1=omp_get_wtime()
        if (limiter.eq.'p-scheme') then
         call WagaanPS(1,gammaa,dx,dy,dt,w &
               ,dwdx,dwdy,wxl,wxr,wyl,wyr,lambdaPS)        
        end if 
        if (limiter.eq.'Cada3rd'.or.limiter.eq.'positiveCada') then
         call Cada3rd(dx,dy,w,wxl,wxr,wyl,wyr)
        if (limiter.eq.'positiveCada') then
         call posiCada_mod(gammaa,dx,dy,w &
                ,wxl,wxr,wyl,wyr,Bi,lambdaPS)
        end if
        end if
        t_MUSCL(2)=t_MUSCL(2)+(omp_get_wtime()-ompt1) 
! Prediction step: Evaluate cell edge values
        if (limiter.ne."p-scheme" &
           .and.limiter.ne.'PPM'.and.limiter.ne.'Cada3rd'.and. &
            limiter.ne.'positiveCada') then
        ompt1 = omp_get_wtime()
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nx,ny,w,wxL,wxR,wyL,wyR,dwdx,dwdy) &
!$OMP PRIVATE (i,j,m)
        !$OMP DO COLLAPSE (3) 
        do m=1,8
         do j=2,ny-1
          do i=3,nx-1 
           wxL(i-1,j,m)=w(i-1,j,m)+0.5*dwdx(i-1,j,m)
           wxR(i,j,m)=w(i,j,m)-0.5*dwdx(i,j,m)
          end do
         end do
        end do
        !$OMP END DO
        !$OMP DO COLLAPSE (3)
        do m=1,8
         do j=3,ny-1
          do i=2,nx-1
           wyL(i,j-1,m)=w(i,j-1,m)+0.5*dwdy(i,j-1,m)
           wyR(i,j,m)=w(i,j,m)-0.5*dwdy(i,j,m)
          end do
         end do
        end do
        !$OMP END DO
!$OMP END PARALLEL        

        ! Account for B.C. - MUST not involve ghost cells 
        wxL(nx-1,2:ny-1,:)=wxR(nx-1,2:ny-1,:)        
        wyL(2:nx-1,ny-1,:)=wyR(2:nx-1,ny-1,:)        
        wxR(2,2:ny-1,:)=wxL(2,2:ny-1,:)       
        wyR(2:nx-1,2,:)=wyL(2:nx-1,2,:)
        t_MUSCL(3)=t_MUSCL(3)+(omp_get_wtime()-ompt1)
        end if
        
! Constrained method transport
        ompt1 = omp_get_wtime()
        if (init_ct.eq.0) then
         Bi(:,:,1)=w(:,:,5)
         Bi(:,:,2)=w(:,:,6)
        else 
          wxL(:,:,5)=EOSHIFT(Bi(:,:,1),DIM=1,SHIFT=+1)
          wxR(:,:,5)=Bi(:,:,1)
          wxL(nx,:,5)=wxL(nx-1,:,5)
          wyL(:,:,6)=EOSHIFT(Bi(:,:,2),DIM=2,SHIFT=+1)
          wyR(:,:,6)=Bi(:,:,2)
          wyL(:,ny,6)=wyL(:,ny-1,6) 
        end if
        t_MUSCL(4)=t_MUSCL(4)+(omp_get_wtime()-ompt1)
         
! Residuals - 
        ompt1 = omp_get_wtime()
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (w,wxL,wxR,gammaa,res,fluxX,init_ct, &
!$OMP wyL,wyR,fluxY,fluxMth,dx,dy,nx,ny) &
!$OMP PRIVATE (i,j,flux)
        !$OMP DO COLLAPSE (2)
        ! Compute residuals in x-direction
        do j=2,ny-1
         do i=3,nx-1
           ! compute flux at i+1/2 using
           call riemannS(fluxMth,w(i-1,j,:),wxL(i-1,j,:),wxR(i,j,:) &
                         ,gammaa,[1,0],flux)
           ! Contribution to the residual of cell (i,j) 
           res(i-1,j,:)=res(i-1,j,:)+flux/dx
           res(i,j,:)=res(i,j,:)-flux/dx
           fluxX(i,j,:)=flux
         end do
        end do
        !$OMP END DO
        
        !$OMP DO COLLAPSE (2)
        ! Compute residuals in y-direction
        do j=3,ny-1
         do i=2,nx-1
           ! compute flux at j+1/2 using
           call riemannS(fluxMth,w(i,j-1,:),wyL(i,j-1,:),wyR(i,j,:) &
                         ,gammaa,[0,1],flux)
           if (init_ct.ne.0) then
           if (wyL(i,j-1,6).ne.wyR(i,j,6)) then
            write(6,*) 'Error as must be same for CT'
            stop
           end if
           end if
           ! Contribution to the residual of cell (i,j) and cells
           ! around it         
           res(i,j-1,:)=res(i,j-1,:)+flux/dy
           res(i,j,:)=res(i,j,:)-flux/dy
           fluxY(i,j,:)=flux 
         end do
        end do
        !$OMP END DO
!$OMP END PARALLEL
        t_MUSCL(5)=t_MUSCL(5)+(omp_get_wtime()-ompt1)
        
! Set BCs for residual
        ompt1 = omp_get_wtime()
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (w,wxL,wxR,gammaa,res,fluxX,init_ct, &
!$OMP wyL,wyR,fluxY,fluxMth,dx,dy,nx,ny) &
!$OMP PRIVATE (i,j,flux)
        
        ! Flux contribution of the most NORTH face: north face of cells
        ! j=M-1
        !$OMP DO
        do i=2,nx-1
          ! Obtain flux based on flux method selected          
          call riemannS(fluxMth,w(i,ny-1,:),wyL(i,ny-1,:),wyR(i,ny-1,:) &
                        ,gammaa,[0,1],flux)
          
          res(i,ny-1,:)=res(i,ny-1,:)+flux/dy
          fluxY(i,ny-1,:)=flux
        end do
        !$OMP END DO

        ! Flux contribution of the most SOUTH face: south face of cells
        ! j=2
        !$OMP DO
        do i=2,nx-1
          ! Obtain flux based on flux method selected         
          call riemannS(fluxMth,w(i,2,:),wyL(i,2,:),wyR(i,2,:) &
                        ,gammaa,[0,1],flux)
          
          res(i,2,:)=res(i,2,:)-flux/dy
          fluxY(i,2,:)=flux
        end do
        !$OMP END DO
        
        ! Flux contribution of the most EAST face: east face of cells
        ! i=N-1
        !$OMP DO
        do j=2,ny-1
          ! Obtain flux based on flux method selected          
          call riemannS(fluxMth,w(nx-1,j,:),wxL(nx-1,j,:),wxR(nx-1,j,:) &
                       ,gammaa,[1,0],flux)
          
          res(nx-1,j,:)=res(nx-1,j,:)+flux/dx
          fluxX(nx-1,j,:)=flux          
        end do
        !$OMP END DO

        ! Flux contribution of the most WEST face: west face of cells
        ! i=2
        !$OMP DO
        do j=2,ny-1
          ! Obtain flux based on flux method selected          
          call riemannS(fluxMth,w(2,j,:),wxL(2,j,:),wxR(2,j,:) &
                       ,gammaa,[1,0],flux)
       
          res(2,j,:)=res(2,j,:)-flux/dx
          fluxX(2,j,:)=flux
        end do
        !$OMP END DO
!$OMP END PARALLEL
        t_MUSCL(6)=t_MUSCL(6)+(omp_get_wtime()-ompt1)
        t_MUSCL(7)=t_MUSCL(7)+(omp_get_wtime()-ompt0)
        end subroutine 

