!-------------------------------------------------------------------------------
!truncate k-vectors (average over identical k's with the same magnitude)  
!
!Copyright 2014 Daniel C. Elton
!-------------------------------------------------------------------------------

module truncate_datasets

contains

subroutine truncate 
 use main_stuff
 Implicit none 

 write(*,*) "Total number of k vectors (including diagonals) = ", Nk
!figure out how many of each magnitude there are
 n = 1
 num_this_mag = 1

 do i = 2, Nk
	if ( (magk(i) - magk(i-1)) .lt. .001 ) then 
		num_this_mag(n) = num_this_mag(n) + 1
	else
		n = n + 1
	endif 
 enddo

 num_ind_mags = n

 write(*,*) "Number of k vector magnitudes  = ", num_ind_mags

 allocate(magk_tr(num_ind_mags))
 allocate(chik0_tr(num_ind_mags))
 allocate(chik0_self_tr(num_ind_mags))

 allocate(eps0T_tr(num_ind_mags))
 allocate(chik0T_tr(num_ind_mags))
 allocate(str_fac_tr(num_ind_mags))

 allocate(chik0_err_tr(num_ind_mags))
 allocate(str_fac_err(num_ind_mags))

  allocate(str_fackt_tr(num_ind_mags,maxsteps))
 ! allocate(rhokt_tr(num_ind_mags,maxsteps))
 !allocate(polTkt_tr(num_ind_mags,maxsteps,3))
 allocate(phiL_tr(num_ind_mags,nsteps_out))
 allocate(phiT_tr(num_ind_mags,nsteps_out))


 write(*,'(a)',advance='no') "truncating datasets .."
!truncate all the results and figure out error bars for certain things
ix = 1
do n = 1, num_ind_mags

	magk_tr(n) =   sum( magk(ix:ix+num_this_mag(n)-1))/real(num_this_mag(n)) 

	chik0_tr(n) =  sum(chik0(ix:ix+num_this_mag(n)-1))/real(num_this_mag(n)) 

	!error in chi(k,0)
	if (num_this_mag(n) .gt. 3) then
		chik0_err_tr(n) = sqrt(sum( ( chik0(ix:ix+num_this_mag(n)-1)- chik0_tr(n) )**2))/real(num_this_mag(n))
	else 
		chik0_err_tr(n) =  0 
	endif

	chik0_self_tr(n) =  sum(chik0_self(ix:ix+num_this_mag(n)-1))/real(num_this_mag(n))  

 	str_fac_tr(n) = sum(str_fac(ix:ix+num_this_mag(n)-1))/real(num_this_mag(n))

	if (num_this_mag(n) .gt. 3) then
		str_fac_err(n) = real( sqrt(sum( ( str_fac(ix:ix+num_this_mag(n)-1) - str_fac_tr(n) )**2) )/real(num_this_mag(n)) )
	else 
		str_fac_err(n) =  0 
	endif

	str_fackt_tr(n,:) =  sum(str_fackt(ix:ix+num_this_mag(n)-1,:),1)/real(num_this_mag(n))
  
	!rhokt_tr(n,:)     =  sum(rhokt(ix:ix+num_this_mag(n)-1,:),1)/real(num_this_mag(n))  
	!polTkt_tr(n,:,:)  =  sum(polTkt(ix:ix+num_this_mag(n)-1,:,:),1)/real(num_this_mag(n))  

	phiL_tr(n,:) =  sum(phiL(ix:ix+num_this_mag(n)-1,1:nsteps_out),1)/real(num_this_mag(n))  
	phiT_tr(n,:) =  sum(phiT(ix:ix+num_this_mag(n)-1,1:nsteps_out),1)/real(num_this_mag(n))  


	ix = ix + num_this_mag(n)
enddo
 write(*,*) "... done"



end subroutine truncate




end module
