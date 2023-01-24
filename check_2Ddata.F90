        subroutine check_2Ddata(u)
        use param
        integer :: i, j
        real :: u(nx,ny)
        
        open(50,file='data_check')
        do j=2,ny-1 ! Takes adv of Fortran's column major order
         do i=2,nx-1
          write(50,'(E15.8,A)',advance='no') u(i,j),' '
         end do
         write(50,*) ''
        end do
        close(50)
        end subroutine


        subroutine savedata(x,y,p,q,ictype)
        use param
        character (len=15)::ictype
        integer:: i,j
        real :: x(nx,ny),y(nx,ny)
        real :: q(nx,ny,8),p(nx,ny)
        open(30,file=trim(ictype)//'/data_x')
        do j=2,ny-1 ! Takes adv of Fortran's column major order
         do i=2,nx-1
          write(30,'(E15.8,A)',advance='no') x(i,j),' '
         end do
         write(30,*) ''
        end do
        close(30)
        open(31,file=trim(ictype)//'/data_y')
        do j=2,ny-1
         do i=2,nx-1
          write(31,'(E15.8,A)',advance='no') y(i,j),' '
         end do
         write(31,*) ''
        end do
        close(31)

        open(41,file=trim(ictype)//'/data_density')
        do j=2,ny-1
         do i=2,nx-1
          write(41,'(E15.8,A)',advance='no') q(i,j,1),' '
         end do
         write(41,*) ''
        end do
        close(41)
        open(42,file=trim(ictype)//'/data_u')
        do j=2,ny-1
         do i=2,nx-1
          write(42,'(E15.8,A)',advance='no') q(i,j,2)/q(i,j,1),' '
         end do
         write(42,*) ''
        end do
        close(42)
        open(44,file=trim(ictype)//'/data_v')
        do j=2,ny-1
         do i=2,nx-1
          write(44,'(E15.8,A)',advance='no') q(i,j,3)/q(i,j,1),' '
         end do
         write(44,*) ''
        end do
        close(44)
        open(43,file=trim(ictype)//'/data_pressure')
        do j=2,ny-1
         do i=2,nx-1
         write(43,'(E15.8,A)',advance='no') p(i,j),' '
         end do
         write(43,*) ''
        end do
        close(43)
        open(45,file=trim(ictype)//'/data_Bx')
        do j=2,ny-1
         do i=2,nx-1
          write(45,'(E15.8,A)',advance='no') q(i,j,5),' '
         end do
         write(45,*) ''
        end do
        close(45)
        open(46,file=trim(ictype)//'/data_By')
        do j=2,ny-1
         do i=2,nx-1
          write(46,'(E15.8,A)',advance='no') q(i,j,6),' '
         end do
         write(46,*) ''
        end do
        close(46)
        end subroutine
