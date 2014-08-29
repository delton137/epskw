program test_correlation_functions 
 use correlation_function
 Implicit none
 integer, parameter :: nsteps = 10000
 double complex, dimension(nsteps) :: input 		
 real(8), dimension(nsteps) :: output, output2
 real :: rand1, rand2
 INTEGER, DIMENSION (1) :: seed = (/3/) 
 integer :: i 

 CALL RANDOM_SEED(PUT=seed)

 do i = 1, nsteps
        !CALL RANDOM_NUMBER(rand1)
 	!CALL RANDOM_NUMBER(rand2)
	!input(i) = cmplx(rand1-.5,rand1-.5)
	input(i) = sin(2*3.14159*real(i)/1000.0)
 enddo  


 call simple_complex_corr_function(input, output, nsteps, nsteps)
 call calc_corr_function(input,output2,nsteps)


 do i = 1, nsteps
	write(*,'(1i,3f12.6)') i, real(input(i)), real(output(i)), real(output2(i))		
 enddo



end program test_correlation_functions 		



