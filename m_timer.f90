module m_timer
implicit none

! saved vars needed between both subroutines to complete calculation.
integer(kind(1.0d0)),save :: values_init(8),values_final(8)
integer(kind(1.0d0)) :: time_elapsed

  contains
    subroutine set_timer
    character(8) :: date
    character(10) :: time
    character(5) :: zone

! get values_init for later use.
    call DATE_AND_TIME(date,time,zone,values_init)
        
    end subroutine

    subroutine elapsed_time(seconds)
    character(8) :: date
    character(10) :: time
    character(5) :: zone
    integer(kind(1.0d0)) :: values(8)
    integer(kind(1.0d0)) :: values_change(8)
    real(kind(1.0d0)),intent(out) :: seconds
! get values_final.
    call DATE_AND_TIME(date,time,zone,values_final)
! calculate change in each value entry.    
    values_change = values_final - values_init

! convert days, hours, minutes, and ms into seconds and accum.
    time_elapsed = 86400000 * values_change(3)
    time_elapsed = time_elapsed + 3600000 * values_change(5)
    time_elapsed = time_elapsed + 60000 * values_change(6)
    time_elapsed = time_elapsed + 1000*values_change(7)
    time_elapsed = time_elapsed + values_change(8)
! finally calculate real-valued number of seconds.
    seconds = real(time_elapsed) / 1000d0

    write(*,*) "elapsed time = ", int(real(seconds)/60), " minutes", mod(seconds,60d0), " seconds"
    !return seconds
    end subroutine 

End module
