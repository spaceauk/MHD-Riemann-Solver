!#######################################################
!#
!#    Main file for MHD
!#
!########################################################

        program main

        use omp_lib
        use param
        use timers
        implicit none
        integer:: loopcount,backup
        integer:: i,j,n
        character(len = 15):: IC,fluxMth,limiter
        real:: Lx,Ly,dx,dy,dt,tEnd
        real(8):: t
        real:: gammaa
        real:: beta ! parameter for measuring plasma stability
        real:: lambda1,lambda2,lambda3,lambda4,lambda
! Measuring time
        integer:: rate,wstart,wfinish
        real:: start,finish,cput0,cput1

! Allocate memory spaces (size depends on test case used)
        real,allocatable:: x(:,:),y(:,:)
        real,allocatable:: u(:,:),v(:,:),w(:,:)
        real,allocatable:: speed(:,:)
        real,allocatable:: r(:,:),p(:,:) ! r is for density
        real,allocatable:: Bx(:,:),By(:,:),Bz(:,:)
        real,allocatable:: c(:,:) ! speed of sound
        real,allocatable:: E(:,:) ! energy
        real,allocatable:: ca(:,:) !Alfven waves
        real,allocatable:: cax(:,:),cay(:,:)!Alfven waves
        real,allocatable:: cfx(:,:),cfy(:,:)!Fast M waves
        real,allocatable:: q(:,:,:) ! MHD conserved quantities
        real,allocatable:: qs(:,:,:),qss(:,:,:),res1(:,:,:)
! For CT method
        integer:: init_ct,CT_energyc
        real,allocatable:: Bi(:,:,:),Bis(:,:,:),Biss(:,:,:)
        real,allocatable:: fluxX(:,:,:),fluxY(:,:,:)
! For yhPS/CT stuffs
        real,allocatable:: diff(:,:,:)
        real:: lambdaPS
! For OpenMP
        integer:: ithr,num_thr
        real:: ompt0,ompt1

#ifdef OPENMP
!$OMP PARALLEL private (ithr)
        num_thr=OMP_GET_NUM_THREADS()
        ithr=OMP_GET_THREAD_NUM()
        if (ithr.lt.5.or.ithr.gt.(num_thr-6)) then
         write(6,*) 'Check OpenMP: Thr ID, No. of thr =',ithr,num_thr
        end if
!$OMP END PARALLEL        
#endif

        call input
! Define terms
        tEnd=0.27     ! Final time; (Temp fix: 0.25 switch to 0.15) 
        n=3        ! Degrees of freedom: ideal air=5, monoatomic gas=3;
        gammaa=1.4 ! Ratio of specific heats for ideal di-atomic gas
        fluxMth='HLLC_Linde'    ! RUSA, HLLE, HLLC, HLLC_Linde, HLLD
        limiter='1st'    ! 1st, MC, MM, Albada, VanLeer, p-scheme, Cada3rd, positiveCada
        init_ct=2           ! 0(no CT),1(default),2(variant) 
        CT_energyc=1        ! energy correction for CT method
        loopcount=0         ! Do NOT change. Req to start CT!
        backup=0
!  BWx,BWy,rotor2D,Wagan_rotor2D,Sedov2D,extremeBlast,Iso2D_p100,CAFEQ2D
        IC = 'Wagan_rotor2D' ! Choose test case option

! Change domain size due to certain test cases
        Lx=1.
        Ly=Lx
        if (IC.eq.'Sedov2D') then
         Ly=Lx*1.5
         ny=nx*15/10 ! Cause is integer
        else if (IC.eq.'rotateBW') then
         if (mod(nx,25).ne.0) nx=nx-mod(nx,25)
         Lx=1.
         Ly=0.16
         ny=nx*4/25
         !write(6,*) nx,ny
        end if
! Account for 'ghost' cells
        nx=nx+2
        ny=ny+2
! Allocate memory based on test cases
        allocate(x(nx,ny),y(nx,ny))
        allocate(r(nx,ny),p(nx,ny),E(nx,ny))
        allocate(u(nx,ny),v(nx,ny),w(nx,ny))
        allocate(Bx(nx,ny),By(nx,ny),Bz(nx,ny))
        allocate(speed(nx,ny),c(nx,ny),ca(nx,ny))
        allocate(cax(nx,ny),cay(nx,ny))
        allocate(cfx(nx,ny),cfy(nx,ny))
        allocate(q(nx,ny,8),qs(nx,ny,8),qss(nx,ny,8))
        allocate(res1(nx,ny,8))
        allocate(Bi(nx,ny,2),Bis(nx,ny,2),Biss(nx,ny,2))
        allocate(fluxX(nx,ny,8),fluxY(nx,ny,8))
! Allocate memory for yhPS/CT stuffs
        allocate(diff(nx,ny,2))
       
! Types of initialization condition 
        call IC2Dtype(IC,Lx,Ly,gammaa,dx,dy,x,y &
                     ,tEnd,r,u,v,w,Bx,By,Bz,p,E)
        write(6,*) 'No of grid points used=',nx,ny
        write(6,*) 'Flux method =',fluxMth
        write(6,*) 'Slope limiter =',limiter
        if (init_ct.ne.0) write(6,*) 'Constrained transport used.'
! Calculate the various waves speed        
        c=sqrt(gammaa*p/r)
        ca=sqrt((Bx**2+By**2+Bz**2)/r)
        cax=sqrt((Bx**2)/r)
        cay=sqrt((By**2)/r)
        cfx=sqrt(0.5*(c**2+ca**2)+0.5*sqrt((c**2+ca**2)**2- &
                 4*(c**2*cax**2)))
        cfy=sqrt(0.5*(c**2+ca**2)+0.5*sqrt((c**2+ca**2)**2- &
                 4*(c**2*cay**2)))
        beta=minval((2*p/(Bx**2+By**2+Bz**2)))
! Concatenate above arrays to one array (conserved quantities q)
        q(:,:,1)=r
        q(:,:,2)=r*u
        q(:,:,3)=r*v
        q(:,:,4)=r*w
        q(:,:,5)=Bx
        q(:,:,6)=By
        q(:,:,7)=Bz
        q(:,:,8)=r*E
! For B.C.(Soft)
        q(1,:,:)=q(2,:,:)
        q(:,1,:)=q(:,2,:)
        q(nx,:,:)=q(nx-1,:,:)
        q(:,ny,:)=q(:,ny-1,:)

! For discretizing time 
        speed=sqrt(u**2+v**2+w**2)
        lambda1=maxval(abs(speed+cfx))
        lambda2=maxval(abs(speed-cfx))
        lambda3=maxval(abs(speed+cfy))
        lambda4=maxval(abs(speed-cfy))
        lambda=maxval([lambda1,lambda2,lambda3,lambda4])
        dt=CFL*minval([dx,dy])/lambda
        t=0 ! initialize time
        write(6,*)'Begin loop and expected to end at tEnd=' &
                ,tEnd,' where 1st dt=',dt
        write(6,*) 'Beta=',beta

        call system_clock(count_rate=rate)
        call system_clock(wstart)
        call cpu_time(start)
!###########################################
!#      Main part 
!###########################################
        do while (((t.lt.real(tEnd,8)).or.(dt.lt.0.001*tEnd)) )!&
                !.and.(loopcount.le.4))
        if (dt.lt.0.0001*tEnd) then 
         write(6,*)'Small dt encountered as dt=',dt
         stop
        end if
        if (t+dt.gt.real(tEnd,8)) then
         dt=tEnd-t
         t=t+real(dt,8)
        else
         t=t+real(dt,8)
        end if
! Main code with subroutines in it 
        ! RK3 1st step / update qs
        ompt1=omp_get_wtime()
        call MUSCL2D(dx,dy,dt,gammaa,q,limiter,fluxMth &
              ,res1 & ! Default output
              ,init_ct*loopcount,Bi,fluxX,fluxY  &! Output for CT
              ,diff,lambdaPS)  ! required inputs/outputs for yhPS/CT stuffs              
        tcpu(1)=tcpu(1)+(omp_get_wtime()-ompt1)
        
        qss=q-dt*res1      

        ompt1=omp_get_wtime()
        if (init_ct.eq.1) then
         call CT2D(1,fluxX,fluxY,Bi,dt,dx,dy, &
                   Biss,qss,CT_energyc)
        else if (init_ct.eq.2) then
         call CT2D_variant(1,fluxX,fluxY,Bi,dt,dx,dy, &
                   Biss,qss,CT_energyc)
        end if        
        tcpu(2)=tcpu(2)+(omp_get_wtime()-ompt1)  
        
        ! For natural B.C.
        qss(1,:,:)=qss(2,:,:)
        qss(:,1,:)=qss(:,2,:)
        qss(nx,:,:)=qss(nx-1,:,:)
        qss(:,ny,:)=qss(:,ny-1,:)
         
        ! To help with CT initialization (placed after 1st RK step)
        loopcount=loopcount+1

        ! RK3 2nd step / update q
        ompt1=omp_get_wtime()
        call MUSCL2D(dx,dy,dt,gammaa,qss,limiter,fluxMth &
               ,res1 & ! Default output
               ,init_ct,Biss,fluxX,fluxY & ! Output for CT
               ,diff,lambdaPS) ! requried input/output for yhPS/CT stuffs        
        tcpu(1)=tcpu(1)+(omp_get_wtime()-ompt1)  

        qs=0.75*q+0.25*qss-0.25*dt*res1        
        Biss=Bi 

        ompt1=omp_get_wtime()
        if (init_ct.eq.1) then
         call CT2D(2,fluxX,fluxY,Biss,dt,dx,dy, &
                   Bis,qs,CT_energyc)
        else if (init_ct.eq.2) then
         call CT2D_variant(2,fluxX,fluxY,Biss,dt,dx,dy, &
                   Bis,qs,CT_energyc)
        end if
        tcpu(2)=tcpu(2)+(omp_get_wtime()-ompt1)  
        ! For natural B.C.
        qs(1,:,:)=qs(2,:,:)
        qs(:,1,:)=qs(:,2,:)
        qs(nx,:,:)=qs(nx-1,:,:)
        qs(:,ny,:)=qs(:,ny-1,:)


        ! RK3 3rd step / update q
        ompt1=omp_get_wtime()
        call MUSCL2D(dx,dy,dt,gammaa,qs,limiter,fluxMth &
               ,res1 & ! Default output
               ,init_ct,Bis,fluxX,fluxY & ! Output for CT
               ,diff,lambdaPS) ! requried input/output for yhPS/CT stuffs        
        tcpu(1)=tcpu(1)+(omp_get_wtime()-ompt1)  

        q=(1./3.)*q+(2./3.)*qs-(2./3.)*dt*res1       
        Bis=Biss 

        ompt1=omp_get_wtime()
        if (init_ct.eq.1) then
         call CT2D(2,fluxX,fluxY,Bis,dt,dx,dy, &
                   Bi,q,CT_energyc)
        else if (init_ct.eq.2) then
         call CT2D_variant(2,fluxX,fluxY,Bis,dt,dx,dy, &
                   Bi,q,CT_energyc)
        end if
        tcpu(2)=tcpu(2)+(omp_get_wtime()-ompt1)  
        ! For natural B.C.
        q(1,:,:)=q(2,:,:)
        q(:,1,:)=q(:,2,:)
        q(nx,:,:)=q(nx-1,:,:)
        q(:,ny,:)=q(:,ny-1,:)


        
! Extract impt variables from new q: e.g. obtain pressure
        ompt1=omp_get_wtime()
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nx,ny,gammaa,q,p) PRIVATE (i,j) 
        !$OMP DO COLLAPSE (2)               
        do j=1,ny
         do i=1,nx
           p(i,j)=(gammaa-1)*(q(i,j,8)-(0.5*(q(i,j,2)**2+q(i,j,3)**2+ &
            q(i,j,4)**2)/q(i,j,1))-0.5*(q(i,j,5)**2+ &
            q(i,j,6)**2+q(i,j,7)**2))
         end do
        end do
!$OMP END PARALLEL
        tcpu(3)=tcpu(3)+(omp_get_wtime()-ompt1)  
        
! Backup data 
        if (t.ge.0.27.and.backup.eq.0) then
        write(6,*) 'Backup data at t=0.27 ...'
        call savedata(x,y,p,q,IC)       
        backup=1
        stop  
        end if

! Check for any errors
        ompt1=omp_get_wtime()
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nx,ny,p,q) &
!$OMP PRIVATE (i,j)        
        !$OMP DO
        do j=2,ny-1
         do i=2,nx-1
          if ((p(i,j).lt.0).or.(q(i,j,1).lt.0)) then
           write(6,*) '-ve p',p(i,j), 'at',i,j
           write(6,*) '-ve density',q(i,j,1),' at',i,j    
           stop
          end if        
          if (p(i,j).ne.p(i,j)) then
           write(6,*) 'NaN p values at',i,j
          end if
         end do
        end do
!$OMP END PARALLEL
        ! Check for any NaN values
        if (sum(p(2:nx-1,2:ny-1)).ne.sum(p(2:nx-1,2:ny-1))) then
         write(6,*) 'Error as NaN values obtained! Sum of p:' &
                     ,sum(p(2:nx-1,2:ny-1))
         write(6,*) sum(q(:,:,1)),sum(q(:,:,2)),sum(q(:,:,3)), &
              sum(q(:,:,5)),sum(q(:,:,6)),sum(q(:,:,8))
         call check_2Ddata(p)
         exit
        end if
        tcpu(4)=tcpu(4)+(omp_get_wtime()-ompt1)  

        ompt1=omp_get_wtime()
        if (limiter.eq.'p-scheme' & !) then !&
                .or.limiter.eq.'positiveCada') then
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nx,ny,q) PRIVATE (i,j) &
!$OMP REDUCTION (max: lambda1,lambda2,lambda3)
        !$OMP DO COLLAPSE (2)
        do j=1,ny
         do i=1,nx
          lambda1=max(lambda1, &
                q(i,j,2)/q(i,j,1)+abs(q(i,j,5))/sqrt(q(i,j,1)))
          lambda2=max(lambda2, &
                q(i,j,3)/q(i,j,1)+abs(q(i,j,6))/sqrt(q(i,j,1)))
          lambda3=max(lambda3, &
                q(i,j,4)/q(i,j,1)+abs(q(i,j,7))/sqrt(q(i,j,1)))
         end do
        end do
!$OMP END PARALLEL
        dt=CFL*min(dx,dy)/maxval([lambda1,lambda2,lambda3,lambdaPS])
        
        else
!
! Collect the different waves speed
        lambda1=0; lambda2=0; lambda3=0; lambda4=0;   
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (c,gammaa,q,p,ca,cax,cay,cfx,cfy,speed,nx,ny) &
!$OMP PRIVATE (i,j) &
!$OMP REDUCTION (max: lambda1,lambda2,lambda3,lambda4)
        !$OMP DO
        do j=1,ny
         do i=1,nx
          c(i,j)=sqrt(gammaa*p(i,j)/q(i,j,1))
          ca(i,j)=sqrt((q(i,j,5)**2+q(i,j,6)**2+q(i,j,7)**2)/ &
                  q(i,j,1))
          cax(i,j)=sqrt((q(i,j,5)**2)/q(i,j,1))
          cay(i,j)=sqrt((q(i,j,6)**2)/q(i,j,1))
          cfx(i,j)=sqrt(0.5*(c(i,j)**2+ca(i,j)**2)+0.5*sqrt((c(i,j)**2 &
                 +ca(i,j)**2)**2- &
                 4*(c(i,j)**2*cax(i,j)**2)))
          cfy(i,j)=sqrt(0.5*(c(i,j)**2+ca(i,j)**2)+0.5*sqrt((c(i,j)**2 &
                 +ca(i,j)**2)**2- &
                 4*(c(i,j)**2*cay(i,j)**2)))
          speed(i,j)=sqrt((q(i,j,2)**2+q(i,j,3)**2+q(i,j,4)**2) &
               /q(i,j,1))
          lambda1=max(lambda1,abs(speed(i,j)+cfx(i,j)))
          lambda2=max(lambda2,abs(speed(i,j)-cfx(i,j)))
          lambda3=max(lambda3,abs(speed(i,j)+cfy(i,j)))
          lambda4=max(lambda4,abs(speed(i,j)-cfy(i,j)))
         end do
        end do      
!$OMP END PARALLEL      

! Update time
        lambda=maxval([lambda1,lambda2,lambda3,lambda4])
        dt=CFL*minval([dx,dy])/lambda        
        end if
        tcpu(5)=tcpu(5)+(omp_get_wtime()-ompt1)  

        write(6,101)'Loop counts=',loopcount,'|  dt=',dt,'s|  t=',t,'s'
101     format(A,I4,A,E12.5,A,F18.16,A)
        
        end do
!#######################################################################
!#                     End main loop
!########################################################################
        call cpu_time(finish)
        call system_clock(wfinish)
        write(6,'(A,F11.6,A)') 'CPU time=',finish-start,'s'
        write(6,'(A,F11.6,A)') 'Wall time=', &
              real(wfinish-wstart)/real(rate),'s'
        tcpu(7)=real(wfinish-wstart)/real(rate)
        tcpu(6)=tcpu(7)-sum(tcpu(1:5))

        call write_timings
! For post-processing
        write(6,*) 'Loop completed. Saving data...'
        call savedata(x,y,p,q,IC)
      
        end program main

