        module param
        integer:: nx,ny
        real:: CFL

        end module param


        module timers
        real:: tcpu(10),t_MUSCL(10),t_CT(10)
        end module timers
