! Subroutine for slope limiter: outputs are dwdx,dwdy
        subroutine slopelimiter(limiter,nv,w,dwdx,dwdy)
        
        use omp_lib
        use param
        implicit none
        character(len = 15):: limiter
        integer:: nv ! where nv is no. of variables
        integer:: i,j,m,num_thr,ithr
        real:: dwdx(nx,ny,nv),dwdy(nx,ny,nv)
        real:: w(nx,ny,nv)
        real:: dww,dwe,dws,dwn,dwf,dwb,dwc,mm        
        real:: eps
        
        eps=2.2204E-16        

!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nv,nx,ny,w,dwdx,dwdy,limiter,eps) &
!$OMP PRIVATE (i,j,m,dww,dwe,dwc,dws,dwn)
! 1st order Godunov
        if (limiter.eq.'1st') then
        !$OMP DO COLLAPSE (3)
        do m=1,nv
         do j=2,ny-1
          do i=2,nx-1
           dwdx(i,j,m)=0
           dwdy(i,j,m)=0
          end do
         end do
        end do
        !$OMP END DO

! Minmod limiter
        else if (limiter.eq.'MM') then
        !$OMP DO COLLAPSE (3)
         do m=1,nv
          do j=2,ny-1
           do i=2,nx-1
             dww=2*(w(i,j,m)-w(i-1,j,m))
             dwe=2*(w(i+1,j,m)-w(i,j,m))
             call minmod([dww,dwe],2,dwdx(i,j,m))
             dws=2*(w(i,j,m)-w(i,j-1,m))
             dwn=2*(w(i,j+1,m)-w(i,j,m))
             call minmod([dws,dwn],2,dwdy(i,j,m))
           end do
          end do
         end do
        !$OMP END DO
        
! Monotonized Central (MC) limiter
        else if (limiter.eq.'MC') then
        !$OMP DO COLLAPSE (3)
         do m=1,nv
          do j=2,ny-1
           do i=2,nx-1
             dww=2*(w(i,j,m)-w(i-1,j,m))
             dwe=2*(w(i+1,j,m)-w(i,j,m))
             dwc=0.5*(w(i+1,j,m)-w(i-1,j,m))
             call minmod([dww,dwe,dwc],3,dwdx(i,j,m))
             dws=2*(w(i,j,m)-w(i,j-1,m))
             dwn=2*(w(i,j+1,m)-w(i,j,m))
             dwc=0.5*(w(i,j+1,m)-w(i,j-1,m))
             call minmod([dws,dwn,dwc],3,dwdy(i,j,m))
           end do
          end do
         end do
        !$OMP END DO
 
! Positive scheme (K. Waagan):
        else if (limiter.eq.'p-scheme') then
        ! Uses the normal MC limiter as highest order
        !$OMP DO COLLAPSE (3)
        do m=1,nv
         do j=2,ny-1
          do i=2,nx-1
             dww=2*(w(i,j,m)-w(i-1,j,m))
             dwe=2*(w(i+1,j,m)-w(i,j,m))
             dwc=0.5*(w(i+1,j,m)-w(i-1,j,m))
             dwdx(i,j,m)=minval(abs([dwe,dwc,dww]))
             if (dwe/2.gt.0.and.dww/2.gt.0) then
              dwdx(i,j,m)=dwdx(i,j,m)
             else if (dwe/2.lt.0.and.dww/2.lt.0) then
              dwdx(i,j,m)=-dwdx(i,j,m)
             else
              dwdx(i,j,m)=0
             end if
             dws=2*(w(i,j,m)-w(i,j-1,m))
             dwn=2*(w(i,j+1,m)-w(i,j,m))
             dwc=0.5*(w(i,j+1,m)-w(i,j-1,m))
             dwdy(i,j,m)=minval(abs([dwn,dwc,dws]))
             if (dwn/2.gt.0.and.dws/2.gt.0) then
              dwdy(i,j,m)=dwdy(i,j,m)
             else if (dwn/2.lt.0.and.dws/2.lt.0) then
              dwdy(i,j,m)=-dwdy(i,j,m)
             else
              dwdy(i,j,m)=0
             end if
           end do
          end do
         end do
        !$OMP END DO

! Albada-type limiter: a smooth function of w
        else if (limiter.eq.'Albada') then
        !$OMP DO COLLAPSE (3)
         do m=1,nv
          do j=2,ny-1
           do i=2,nx-1
             dwdx(i,j,m)= &
              ((2*(w(i+1,j,m)-w(i,j,m))*(w(i,j,m)-w(i-1,j,m))+eps)/ &
              ((w(i+1,j,m)-w(i,j,m))**2+(w(i,j,m)-w(i-1,j,m))**2+eps)) &
              *0.5*(w(i+1,j,m)-w(i-1,j,m))
             dwdy(i,j,m)= &
              ((2*(w(i,j+1,m)-w(i,j,m))*(w(i,j,m)-w(i,j-1,m))+eps)/ &
              ((w(i,j+1,m)-w(i,j,m))**2+(w(i,j,m)-w(i,j-1,m))**2+eps)) &
              *0.5*(w(i,j+1,m)-w(i,j-1,m))
           end do
          end do
         end do
        !$OMP END DO

! Van-Leer limiter:
        else if (limiter.eq.'VanLeer'.or.limiter.eq.'PPM') then       
        !$OMP DO COLLAPSE (3) 
         do m=1,nv
          do j=2,ny-1
           do i=2,nx-1
             dww=(w(i,j,m)-w(i-1,j,m))
             dwe=(w(i+1,j,m)-w(i,j,m))
             dwc=0.5*(w(i+1,j,m)-w(i-1,j,m))
             if (dww*dwe.gt.0) then
              dwdx(i,j,m)=sign((min(abs(dwc),2*min(abs(dww),abs(dwe))))&
                         ,dwc)
             else if (dww*dwe.le.0) then
              dwdx(i,j,m)=0
             end if
             dws=(w(i,j,m)-w(i,j-1,m))
             dwn=(w(i,j+1,m)-w(i,j,m))
             dwc=0.5*(w(i,j+1,m)-w(i,j-1,m))
             if (dws*dwn.gt.0) then
              dwdy(i,j,m)=sign((min(abs(dwc),2*min(abs(dws),abs(dwn))))&
                         ,dwc)
             else if (dws*dwn.le.0) then
              dwdy(i,j,m)=0
             end if
           end do
          end do
         end do
        !$OMP END DO

        else if (limiter.ne.'Cada3rd'.and.limiter.ne.'positiveCada')then
         write(6,*) 'No slope limiter selected. Stopping process...'
         stop
        end if
!$OMP END PARALLEL
        end subroutine slopelimiter


! Need this function for MC limiter
        subroutine minmod(v,dims,mm) 
        implicit none
        real:: mm
        integer:: i,s,dims 
        integer:: signed(dims)
        real:: v(dims)
        
        signed=0
        do i=1,dims
         if (v(i).lt.0) then
          signed(i)=-1
         else if (v(i).gt.0) then
          signed(i)=1
         else if (v(i).eq.0) then
          signed(i)=0
         end if
        end do
        s=sum(signed)/dims
        if (abs(s).eq.1) then
         mm=s*minval(abs(v))
        else 
         mm=0
        end if
        
        end subroutine minmod

