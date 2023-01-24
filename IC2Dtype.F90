        subroutine IC2Dtype(IC,Lx,Ly &
           ,gammaa,dx,dy,x,y,tEnd,r,u,v,w,Bx,By,Bz,p,E)
!#########################################################################
!#
!#   1.0 +-----------+-----------+
!#       |           |           |       
!#       |   reg 2   |   reg 1   |
!#       |           |           |
!#   0.5 +-----------+-----------+
!#       |           |           |
!#       |   reg 3   |   reg 4   |
!#       |           |           |
!#   0.0 +-----------+-----------+
!#      0.0         0.5         1.0
!#
!#########################################################################
        use param
        implicit none
        character(len = 15):: IC
        integer:: i,j
        real:: Lx,Ly,dx,dy
        real:: pi
        real:: tEnd,gammaa
        real:: x(nx,ny),y(nx,ny)
        real:: radius(nx,ny),fr(nx,ny)
        real:: r(nx,ny),p(nx,ny),E(nx,ny)
        real:: u(nx,ny),v(nx,ny),w(nx,ny)
        real:: Bx(nx,ny),By(nx,ny),Bz(nx,ny)
        logical:: reg1(nx,ny),reg2(nx,ny),reg3(nx,ny),reg4(nx,ny)
! Allocate memory spaces for I.C.
        real:: rinit(4),pinit(4),Einit(4)
        real:: uinit(4),vinit(4),winit(4)
        real:: Bxinit(4),Byinit(4),Bzinit(4)
        real:: sound
! Define common variables
        pi=4.D0*DATAN(1.D0)
        
! Initialize all to zero first
        rinit=0
        pinit=0
        Einit=0
        uinit=0
        vinit=0
        winit=0
        Bxinit=0
        Byinit=0
        Bzinit=0
        
! Obtain dx & dy
        dx=Lx/(nx-2)
        dy=Ly/(ny-2)
        if (dx.ne.dy) write(6,*) 'dx \neq dy:',dx,dy

! Discretize spatial grid
        do i=1,nx
        x(i,:)=(i-1)*dx ! -1 to account for BC
        enddo
        do j=1,ny
        y(:,j)=(j-1)*dy
        enddo

! Set conditions for specific I.C.s selected
        if (IC.eq.'BWx') then
        write(6,*) '2D Brio & Wu shocktube configuration along x'
        tEnd=0.2
        gammaa=2.
        pinit=[0.1,1.0,1.0,0.1]
        rinit=[0.125,1.0,1.0,0.125]
        Bxinit=[0.75,0.75,0.75,0.75]
        Byinit=[-1.0,1.0,1.0,-1.0]

        else if (IC.eq.'BWy') then
        write(6,*) '2D Brio & Wu shocktube configuration along y'
        tEnd=0.15
        pinit=[0.1,0.1,1.0,1.0]
        rinit=[0.125,0.125,1.0,1.0]
        Byinit=[0.75,0.75,0.75,0.75]
        Bxinit=[-1.0,-1.0,1.0,1.0]
        
        else if (IC.eq.'rotateBW') then
        write(6,*) 'Rotated Brio & Wu Shocktube configuration'
        write(6,*) 'Ignore this as you need to use diagonal periodic BC' &
                   ,' in order to get this test case to work. But too' &
                   ,' troublesome for me as I need changed a lot.'
        gammaa=2
        tEnd=0.27
        pinit=[0.1,1.,0.,0.]
        rinit=[0.125,1.,0.,0.]
        Bxinit=[7.,-1.,0.,0.]/(4.*sqrt(2.))
        Byinit=[-1.,7.,0.,0.]/(4.*sqrt(2.))
      
        else if (IC.eq.'rotor2D'.or.IC.eq.'Wagan_rotor2D') then
        write(6,*) 'Rotor problem in 2D'
        gammaa=1.4
        tEnd=0.295
        if (IC.eq.'rotor2D') then
         pinit=0.5
        else if (IC.eq.'Wagan_rotor2D') then
         pinit=1E-8
        end if
        Bxinit=2.5/(sqrt(4*pi))
        
        else if (IC.eq.'extremeBlast') then
        write(6,*) 'Extreme Blast wave (p=1000) in 2D'
        gammaa=1.4
        tEnd=0.01
        rinit = 1
        Bxinit = 100/sqrt(4*pi)
        
        else if (IC.eq.'Sedov2D') then
        write(6,*) '2D MHD Sedov problem from Princeton'
        gammaa=5./3.
        tEnd=0.2
        rinit=1
        Bxinit=1./sqrt(2.)
        Byinit=1./sqrt(2.)
      
        else if (IC.eq.'Iso2D_p100') then
        write(6,*)'Isothermal blast wave in 2d (r&p=100) - sensitive to' &
                  ,'errors related to treatment of div.B'  
        gammaa=2
        tEnd=0.09
        Bxinit=5/sqrt(pi)        

        else if (IC.eq.'CAFEQ2D') then
        write(6,*) '2D MHD blast wave in 2d (p=100) - from CAFE-Q'
        gammaa=5./3.
        tEnd=0.02
        rinit=1
        Bxinit=10/sqrt(2.)
        Byinit=10/sqrt(2.)
       
        else
        write(6,*)'Incorrect initial condition selected! End process...'
        stop
        end if

        if (IC.eq.'rotor2D'.or.IC.eq.'Wagan_rotor2D') then
        radius=sqrt((x-0.5)**2+(y-0.5)**2)
        fr=(23.-200.*radius)/3.
        reg1 = (radius.le.0.1)
        reg2 = (radius.gt.0.1.and.radius.lt.0.115)
        reg3 = (radius.ge.0.115)
        ! IC for our 2D domain
        where(reg1)
                    r=10
                    p=pinit(1)
                    u=(-(y-0.5)/0.1)
                    v=((x-0.5)/0.1)
                    w=winit(1)
                    Bx=Bxinit(1)
                    By=Byinit(1)
                    Bz=Bzinit(1)
        endwhere
        where(reg2)
                    r=(1+9*fr)
                    p=pinit(2)
                    u=(-(y-0.5)*fr/0.1)
                    v=((x-0.5)*fr/0.1)
                    w=winit(2)
                    Bx=Bxinit(2)
                    By=Byinit(2)
                    Bz=Bzinit(2)
        endwhere
        where(reg3)
                    r=1
                    p=pinit(3)
                    u=0
                    v=0
                    w=winit(3)
                    Bx=Bxinit(3)
                    By=Byinit(3)
                    Bz=Bzinit(3)
        endwhere
        
        else if (IC.eq.'extremeBlast') then
        radius=sqrt((x-0.5)**2+(y-0.5)**2)
        reg1 = (radius.le.0.1)
        reg2 = (radius.gt.0.1)
        where(reg1)
                    r=rinit(1)
                    p=1000
                    u=uinit(1)
                    v=vinit(1)
                    w=winit(1)
                    Bx=Bxinit(1)
                    By=Byinit(1)
                    Bz=Bzinit(1)
        endwhere
        where(reg2)
                    r=r+rinit(2)
                    p=p+0.1
                    u=u+uinit(2)
                    v=v+vinit(2)
                    w=w+winit(2)
                    Bx=Bx+Bxinit(2)
                    By=By+Byinit(2)
                    Bz=Bz+Bzinit(2)
        endwhere
        
        else if (IC.eq.'Sedov2D') then
        radius =sqrt((x-0.5)**2+(y-0.75)**2)
        reg1 = (radius.le.0.1)
        reg2 = (radius.gt.0.1)
        where(reg1)
                    r=rinit(1)
                    p=10
                    u=uinit(1)
                    v=vinit(1)
                    w=winit(1)
                    Bx=Bxinit(1)
                    By=Byinit(1)
                    Bz=Bzinit(1)
        endwhere
        where(reg2)
                    r=r+rinit(2)
                    p=p+0.1
                    u=u+uinit(2)
                    v=v+vinit(2)
                    w=w+winit(2)
                    Bx=Bx+Bxinit(2)
                    By=By+Byinit(2)
                    Bz=Bz+Bzinit(2)
        endwhere
        
        else if (IC.eq.'Iso2D_p100') then
        radius =sqrt((x-0.5)**2+(y-0.5)**2)
        reg1 = (radius.le.0.05)
        reg2 = (radius.gt.0.05)
        sound=1. ! Sound speed set to 1
        where(reg1)
                    r=100
                    p=100*sound
                    u=uinit(1)
                    v=vinit(1)
                    w=winit(1)
                    Bx=Bxinit(1)
                    By=Byinit(1)
                    Bz=Bzinit(1)
        endwhere
        where(reg2)
                    r=r+1
                    p=p+1*sound**2
                    u=u+uinit(2)
                    v=v+vinit(2)
                    w=w+winit(2)
                    Bx=Bx+Bxinit(2)
                    By=By+Byinit(2)
                    Bz=Bz+Bzinit(2)
        endwhere

        else if (IC.eq.'CAFEQ2D') then
        radius =sqrt((x-0.5)**2+(y-0.5)**2)
        reg1 = (radius.le.0.125)
        reg2 = (radius.gt.0.125)
        where(reg1)
                    r=1.
                    p=100.
                    u=0
                    v=0
                    w=0
                    Bx=10./sqrt(2.)
                    By=10./sqrt(2.)
                    Bz=0
        endwhere
        where(reg2)
                    r=r+1
                    p=p+1
                    u=u+0
                    v=v+0
                    w=w+0
                    Bx=Bx+10./sqrt(2.)
                    By=By+10./sqrt(2.)
                    Bz=Bz+0
        endwhere

        else if (IC.eq.'rotateBW') then
        reg1 = ((y.gt.x-0.42).and.(x.le.0.58))
        reg2 = ((y.le.x-0.42).and.(x.ge.0.42))
        where(reg1)
                    r=rinit(1)
                    p=pinit(1)
                    u=uinit(1)
                    v=vinit(1)
                    w=winit(1)
                    Bx=Bxinit(1)
                    By=Byinit(1)
                    Bz=Bzinit(1)
        endwhere
        where(reg2)
                    r=r+rinit(2)
                    p=p+pinit(2)
                    u=u+uinit(2)
                    v=v+vinit(2)
                    w=w+winit(2)
                    Bx=Bx+Bxinit(2)
                    By=By+Byinit(2)
                    Bz=Bz+Bzinit(2)
        endwhere

        else
! Masking of regions for initialization
        reg1=(x.ge.0.5.and.y.ge.0.5)
        reg2=(x.lt.0.5.and.y.ge.0.5)
        reg3=(x.lt.0.5.and.y.lt.0.5)
        reg4=(x.ge.0.5.and.y.lt.0.5)

! Initialization
        where(reg1)
                    r=rinit(1)
                    p=pinit(1)
                    u=uinit(1)
                    v=vinit(1)
                    w=winit(1)
                    Bx=Bxinit(1)
                    By=Byinit(1)
                    Bz=Bzinit(1)
        endwhere
        where(reg2)
                    r=r+rinit(2)
                    p=p+pinit(2)
                    u=u+uinit(2)
                    v=v+vinit(2)
                    w=w+winit(2)
                    Bx=Bx+Bxinit(2)
                    By=By+Byinit(2)
                    Bz=Bz+Bzinit(2)
        endwhere
        where(reg3)
                    r=r+rinit(3)
                    p=p+pinit(3)
                    u=u+uinit(3)
                    v=v+vinit(3)
                    w=w+winit(3)
                    Bx=Bx+Bxinit(3)
                    By=By+Byinit(3)
                    Bz=Bz+Bzinit(3)
        endwhere
        where(reg4)
                    r=r+rinit(4)
                    p=p+pinit(4)
                    u=u+uinit(4)
                    v=v+vinit(4)
                    w=w+winit(4)
                    Bx=Bx+Bxinit(4)
                    By=By+Byinit(4)
                    Bz=Bz+Bzinit(4)
        endwhere
        end if
        E=(p/(r*(gammaa-1)))+0.5*(u**2+v**2+ &
             w**2)+0.5*(Bx**2+By**2+Bz**2)/r

        end subroutine
