        subroutine input
        
        use param
        implicit none
        integer:: lu

        lu=1111
        open(lu,file='input')
        read(lu,fmt='(1x)')
        read(lu,*) nx,ny
        read(lu,fmt='(1x)')
        read(lu,*) CFL
        write(6,*) 'CFL used for this simulation is',CFL

        end subroutine input
