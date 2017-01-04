!--------------------------------------------------------------------------
! routines for managing calculation of correlation functions 
! and routine for distance decomposition calculation 
!
! Copyright 2014-2015, Daniel C. Elton <delton 17 at gmail .com> 
! 
! License: The MIT License
!--------------------------------------------------------------------------

module calc_phi

contains

!---------------------------------------------------------------------------------
!-------- Phi(k,t) calculations, truncation of data and normalization -----------
!---------------------------------------------------------------------------------
subroutine calc_phikt
 use main_stuff
 use correlation_function
 use truncate_datasets
Implicit none 

!----------------  Compute autocorrelation functions  

 write(*,'(a)',advance='no') "calculating correlation functions .."
 allocate(phicomponentT(nsteps))
 allocate(phiL(Nk,nsteps))
 allocate(phiT(Nk,nsteps))

 phiL = 0 
 phiT = 0 
 phicomponentT = 0 

 do i = 1, Nk
	!call simple_complex_corr_function(rhokt(i,1:nsteps), phiL(i,1:nsteps_out), nsteps, nsteps_out)
	!call calc_corr_function(rhokt(i,1:nsteps), phiL(i,1:nsteps), nsteps) 
	call calc_corr_function2(rhokt(i,1:nsteps), phiL(i,1:nsteps), nsteps) 

	do ix = 1,3
		!call simple_complex_corr_function(PolTkt(i,1:nsteps,ix), phiTcomponent, nsteps, nsteps_out)
		!call calc_corr_function(PolTkt(i,1:nsteps,ix), phiTcomponent, nsteps) 
		call calc_corr_function2(PolTkt(i,1:nsteps,ix), phicomponentT, nsteps) 
		phiT(i,1:nsteps_out) = phiT(i,1:nsteps_out) + phicomponentT(1:nsteps_out)
	enddo
 enddo
 write(*,*) "... done"

!----------------  truncate results to k with different magnitudes  
 if (.not. ALT_CALC) chik0 = phiL(:,1)

 call truncate

 Nk = num_ind_mags !Nk changes here!!

 !save static transverse part 
 eps0T_tr  = phiT_tr(:,1) 

do n = 1, num_ind_mags
	!2nd normalization of correlation fun
	phiT_tr(n,:) = phiT_tr(n,:)/phiT_tr(n,1)
	phiL_tr(n,:) = phiL_tr(n,:)/phiL_tr(n,1)
enddo
 
!----------------  Normalization & prefactors  
 vol = box(1)*box(2)*box(3)
 prefac = (e2C**2)/(eps_0*kb*temp*vol*a2m) 

 if (DYNAMIC_STR_FAC) prefac = 1/(dble(Nmol)*3d0) 
 
 !prefactors
 if (.not. ALT_CALC) then 
 chik0_tr      = prefac*chik0_tr 
 chik0_err_tr  = prefac*chik0_err_tr 
 endif
 if (ALT_CALC) then 
 chik0_tr      = prefac*chik0_tr/(dble(nsteps))  
 chik0_err_tr  = prefac*chik0_err_tr/(dble(nsteps))
 chik0_self_tr = prefac*chik0_self_tr/(dble(nsteps)) 
 str_fackt_tr  = str_fackt_tr/(dble(Nmol)*dble(nsteps))  
 endif


do n = 1, num_ind_mags
 	chik0_tr(n) = chik0_tr(n)/(magk_tr(n)**2)

	chik0_self_tr(n) = chik0_self_tr(n)/(magk_tr(n)**2)
	
 	chik0T_tr(n)  = prefac*eps0T_tr(n)/(magk_tr(n)**2)  

	eps0T_tr(n) = chik0T_tr(n)  + 1

	str_fac_tr(n) = sum(str_fackt_tr(n,1:nsteps))/3d0 !fudge factor of 1/3 needed here
enddo

end subroutine calc_phikt




!---------------------------------------------------------------------------------
!-------- Phi(k,t) calculations and normalization - DISTANCE-DEPENDENT CASE -----
!---------------------------------------------------------------------------------
subroutine calc_phikRt
 use main_stuff
 Implicit none 

 allocate(phicomponent(nsteps))
 allocate(phicomponentL(nsteps))
 allocate(phicomponentT(nsteps))
 allocate(phiL(Nr,nsteps))
 allocate(phiT(Nr,nsteps))
 allocate(numR(Nr))

 phicomponent = 0 
 phiL = 0 
 numR = 0  

 write(*,*) "calculating phi(k,R,t)"
 call calc_distdep
 !if (IRCALC) call calc_IR_R

 chik0  = phiL(:,1) 
 chik0T = phiT(:,1) 
	
!-------- Normalization & prefactors  
do i = 1, Nr
	phiL(i,:) = phiL(i,:)/phiL(i,1)
	phiT(i,:) = phiT(i,:)/phiT(i,1)
enddo

!----------------  Normalization & prefactors  
 vol = box(1)*box(2)*box(3)
 prefac = (e2C**2)/(eps_0*kb*temp*vol*a2m) 

 if (IRCALC) prefac = (e2C**2)/(eps_0*kb*temp*vol*a2m)/6/(2.99e8)/100 !IR prefactor (1/100 converts from 1/m to 1/cm)

 !prefactors and /k^2
 chik0           = prefac*chik0/(magk(1)**2) 
 chik0T          = prefac*chik0T/(magk(1)**2) 

end subroutine calc_phikRt




!--------------------------------------------------------------------------
!---- calculation with distance-dependent Kirkwood *spheres* -------------
!--------------------------------------------------------------------------
subroutine calc_distdep
 use correlation_function
 use main_stuff
 Implicit none 

!!$OMP PARALLEL 
!!$OMP DO  reduction(+:phiL,phiT)
do i = 1, Nmol 
	spheresT = 0 
	spheresL = 0 	

	do t = 1, nsteps
		!Kirkwood spheres calculation 
		NumR = 0 
		do j = 1, Ncalc
			Rij = rCMs(i,:,t) - rCMs(j,:,t)
			Rij = Rij - box*anint(Rij/box)!PBC
			dist = sqrt(dot_product(Rij,Rij))
			countR = floor(dist/delta_R) + 1

			do  k = countR, Nr
				!cutoff_fn = 1d0/(1d0 +  dexp( ( dist - (k-1)*delta_R )/.01d0   ) ) 
				spheresT(k, :, t) = spheresT(k, :, t) +  mPolsT(j, :, t)!*cutoff_fn
				spheresL(k, t)    = spheresL(k, t)    +  mPolsL(j, t)!*cutoff_fn
				numR(k) = numR(k) + 1
			enddo
		enddo!j = 1, Nmol, Nskip

		!!if sphere-sphere normalize (as recommended by Heyden)
		if (SPHERESPHERE) then
			do k = 1, Nr
				if (numR(k) .eq. 0) then 
					write(*,*) "ERROR, numR = 0 for shell", k
				else
					spheresT(k, :, t) = spheresT(k, :, t)/sqrt(numR(k))	
					spheresL(k, t)    = spheresL(k, t)/sqrt(numR(k))	
				endif	
			enddo
		endif

	enddo !do t = 1, nsteps
 

	!--------------------------------------------------------------------------------- 
	!-------- calculate cross-correlation of central molecule w/ sphere -----------
	!-------- (or calculate auto-correlation of sphere) -----------------------------
	!--------------------------------------------------------------------------------- 
	if (mod(i,100) .eq. 1) write(*,*) "calculating corr fun. ", (i-1*Nr )*4, " of ", Ncalc*Nr*4

	do j = 1, Nr	
		phicomponentL = 0 
		phicomponent  = 0 !shouldn't be necessary 
		if (.not. K_EQ_0_DIST_DEP) then 
			if (DIPSPHERE)    call  calc_cross_corr_function2(mPolsL(i,:)   ,spheresL(j, :), phicomponentL, nsteps) 
			if (SPHERESPHERE) call  calc_cross_corr_function2(spheresL(j, :),spheresL(j, :), phicomponentL, nsteps) 
		endif

		phicomponentT = 0
		do k = 1, 3
			if (DIPSPHERE)    call  calc_cross_corr_function2(mPolsT(i,k,:)    ,spheresT(j, k, :), phicomponent, nsteps) 
			if (SPHERESPHERE) call  calc_cross_corr_function2(spheresT(j, k, :),spheresT(j, k, :), phicomponent, nsteps) 

			!call  calc_cross_corr_function2(spheresT(i, j, k, :),spheresT(i, j, k, :), phicomponent, nsteps) 
			!call simple_complex_cross_corr_function(mPolsT(i,k,:),spheresT(i, j, k, :), phicomponent, nsteps, nsteps)
			 phicomponentT = phicomponentT + phicomponent
		enddo

		phiL(j,:) = phiL(j,:) + phicomponentL
		phiT(j,:) = phiT(j,:) + phicomponentT
 	enddo!j = 1, Nr	

enddo !do i = 1, Nmol, Nskip
!!$OMP END DO
!!$OMP END PARALLEL

end subroutine calc_distdep



end module calc_phi
