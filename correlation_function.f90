module correlation_function

contains 
!Note: these subroutines calculate the correlation function <phi(t)phi*(t)> 

!-------------------------------------------------------------------------------
!------------- Simple direct computation of the correlation function -----------
!-------------------------------------------------------------------------------
subroutine simple_complex_corr_function(input, ACF, nsteps, tcorr)
 Implicit None 
 integer, intent(in) :: nsteps, tcorr
 double complex, dimension(nsteps), intent(in) :: input
 real(8), dimension(tcorr), intent(out) :: ACF
 integer :: t, t0, tmax

 ACF = 0 

 do t = 1, tcorr 
	!we want t = 1 to be t = 0
	tmax = nsteps - t - 1 
	do t0 = 1, tmax	
		ACF(t) = ACF(t) + dble( conjg(input(t0))*input(t0 + t - 1) )
	enddo

	!normalization  1
	if (tmax .gt. 0) ACF(t) = ACF(t)/tmax
 enddo
 
 !normalization  2
! ACF = ACF/ACF(1)

end subroutine simple_complex_corr_function



!--------------------------------------------------------------------------
!----------------  calcuates correlation function using FFTW3 ------------
!---------------- (this is much much faster, scales as n*log(n) ----------
!---------------- calls to fftw3 have been commented out (library issues)-
!--------------------------------------------------------------------------

subroutine calc_corr_function(input,output,N)
 Implicit none
 include 'fftw3.f90'
 double complex, dimension(N), intent(in) :: input
 real(8), dimension(N), intent(out) :: output
 integer, intent(in) :: N
 double complex, dimension(:), allocatable :: input_padded, transformed, transformed2, output_padded
 integer*8 :: plan=0, i, trun
  
 if (.not.(allocated(input_padded))) then 
 	!find closest power of 2 greater than number of steps
 	trun = 2**(    floor( dlog(  dble(N) )/dlog(2d0)  )  + 1  )
	allocate(input_padded(2*trun))
	allocate(output_padded(2*trun))
	allocate(transformed(2*trun))
	!allocate(transformed2(2*trun))
 endif 

 input_padded = 0 
 input_padded(1:N) = input 

! call dfftw_plan_dft_1d_( plan, 2*N, input_padded, transformed, FFTW_FORWARD, FFTW_ESTIMATE )
! call dfftw_execute_dft( plan, input_padded, transformed )

 transformed = transformed*conjg(transformed)

! call dfftw_plan_dft_1d_( plan, 2*N, transformed, output_padded, FFTW_BACKWARD, FFTW_ESTIMATE )
! call dfftw_execute_dft( plan, transformed, output_padded)	
	
 output = dble(output_padded(1:N))/N 

 !normalization 1
 do i = 1, N-1
	output(i) = output(i)/(N-i)
 enddo

 !normalization 2
! output = output/output(1)

end subroutine calc_corr_function



!------------------------------------------------------------------------
!----------------  Compute autocorrelation function using four1.f-------
!----------------  (proprietary code from Numerical Recipes) -----------
!---------------- (has not been fully tested) --------------------------
!------------------------------------------------------------------------ 
subroutine calc_corr_function2(input,output,N)
 Implicit none
 double complex, dimension(N), intent(in) :: input
 real(8), dimension(N), intent(out) :: output
 integer, intent(in) :: N
 double complex, dimension(:), allocatable :: input_padded,  output_padded
 complex, dimension(:), allocatable :: transformed
 integer*8 :: plan=0, i, trun
  
 if (.not.(allocated(input_padded))) then 
 	!find closest power of 2 greater than number of steps
 	trun = 2**(    floor( dlog(  dble(N) )/dlog(2d0)  )  + 1  )
	allocate(input_padded(2*trun))
	allocate(output_padded(2*trun))
	allocate(transformed(2*trun))
 endif 

 input_padded = 0 
 input_padded(1:N) = input 

 transformed = real(input_padded)

 call four1(transformed,2*trun,1)

 transformed = transformed*conjg(transformed)

 call four1(transformed,2*trun,-1)
	
 output = dble(transformed(1:N))/N 

 !normalization 1
 do i = 1, N-1
	output(i) = output(i)/(N-i)
 enddo

end subroutine calc_corr_function2






end module correlation_function

