        subroutine write_timings
        
        use param
        use timers
        implicit none
        character*6 filepos


        filepos='rewind'

        open(77,file='cpu_main',position=filepos)
        write(77,711) nx,ny,CFL
 711    format('MHD code:, nx,ny,CFL',2i6,1e10.3)
        write(77,*) 'In main loop:-------------------------------------'
        write(77,"(' MUSCL     CTmethod  Pressure  Error', &
                '     Timestep  Others    Total')")
        write(77,630) tcpu(1:7)
        write(77,*) 'In MUSCL:-----------------------------------------'
        write(77,"(' Slopelim  P-scheme  CellEdge  CT   ', & 
                '     RSmain    RS_BCs    Total')")
        write(77,630) t_MUSCL(1:7)
 630    format(7e10.3)
        close(77)       

        end subroutine
