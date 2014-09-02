!----------------------------------------------------------------------------- 
! eps(k,w) program
!-----------------------------------------------------------------------------
! Computes correlation functions of the charge-charge structure factor
!        phi(k,t) = <rho(k,t)*rho(-k,t)>/<rho(k,0)>^2/k^2 
! It then attempts to Fourier Transform these to get Im{chi(k,w)}
!
! It outputs the longitudinal susceptibility defined by P(k,w) = \chi(k,w) D(k,w) and eps(k,w) = 1/(1-chi(k,w))
!
! It also outputs the static nonlocal susceptibility chi(k,0) and the static structure factor S(k,0)
!  
! REFERENCE: JCP 109 5 pg 1939
!
! This code is based off the eps_k code
!-----------------------------------------------------------------------------
! 2014 Dan Elton 
!----------------------------------------------------------------------------- 
Program epskw
	use m_timer 
	use main_stuff
	use correlation_function
	use chi_k
	use truncate_datasets
Implicit None

call set_timer
call read_input_file
call set_up_model
call open_trajectory_files
call setup_k_vectors

 chik0_self = 0 
 chik0 = 0 
 str_fac = 0 
 str_fackt = 0 
 rhokt = 0 

!------------------------------------------------------------------------------
!---------------------------------- main loop --------------------------------
!------------------------------------------------------------------------------
nsteps = 0
Do t = 1, maxsteps
  
	call read_trajectory_frame

	if ((filetype .eq. "xtc") .and. (RET.EQ.0)) then 
		write(*,*) "reached end of file at t = ", t 
		write(*,*) "length of traj is ", TIME, "ps "
		goto 1000 
	endif 

	if (TTM3F) then 
		call calc_chik_TTM3F
	else	
	 	call calc_chik 
	endif 

	call calc_chik_transverse

	if (mod(t,10) .eq. 0) write(*,*) t  

 	nsteps = nsteps + 1 

enddo 

1000 continue 

 write(*,*) "number of steps used: ", nsteps

do i = 1, Nk
	 write(*,*) 1d0 +  sum(PolTkt(:,:,1),2) 
enddo

!--------------------------------------------------------------------------------- 
!----------------  truncate results to k with different magnitudes -------------- 
!--------------------------------------------------------------------------------- 
 write(*,'(a)',advance='no') "truncating stuff .."
 call truncate
 write(*,*) "... done"

 Nk = num_ind_mags !Nk changes here!!


!--------------------------------------------------------------------------------- 
!----------------  Compute autocorrelation functions ---------------------------- 
!--------------------------------------------------------------------------------- 
 allocate(phiTcomponent(nsteps))
 allocate(phiL(Nk,nsteps))
 allocate(phiT(Nk,nsteps))

 phiL = 0 
 phiT = 0 
 phiTcomponent = 0 

 do i = 1, Nk
	call simple_complex_corr_function(rhokt_tr(i,1:nsteps), phiL(i,1:nsteps), nsteps, nsteps)
	!call calc_corr_function(rhokt(i,1:nsteps), phiL(i,1:nsteps), nsteps) 

	do ix = 1,3
		call simple_complex_corr_function(PolTkt_tr(i,1:nsteps,ix), phiTcomponent, nsteps, nsteps)
		!qcall calc_corr_function(PolTkt(i,:,ix), phiTcomponent, nsteps) 
		phiT(i,:) = phiT(i,:) + phiTcomponent
	enddo
 enddo

 !save static transverse part 
 eps0T_tr  = phiT(:,1) 


do n = 1, num_ind_mags
	!2nd normalization of correlation fun
	phiT(n,:) = phiT(n,:)/phiT(n,1)
	phiL(n,:) = phiL(n,:)/phiL(n,1)
enddo


!--------------------------------------------------------------------------------- 
!----------------  Normalization & prefactors ----------------------------------- 
!---------------------------------------------------------------------------------
 vol = box(1)*box(2)*box(3)
 prefac = (e2C**2)/(eps_0*kb*temp*vol*a2m) 

 !prefactors
 !PolTkt_tr     = prefac*PolTkt/(dble(nsteps)) 
 chik0_tr      = prefac*chik0_tr/(dble(nsteps))  
 chik0_self_tr = prefac*chik0_self_tr/(dble(nsteps)) 
 str_fackt_tr  = str_fackt_tr/(dble(Nmol)*dble(nsteps))  



do n = 1, num_ind_mags
 	chik0_tr(n) = chik0_tr(n)/(magk_tr(n)**2)

	chik0_self_tr(n) = chik0_self_tr(n)/(magk_tr(n)**2)
	
	eps0T_tr(n) = prefac*eps0T_tr(n)/(magk_tr(n)**2) + 1d0	

 	chik0T_tr(n)  = 1d0 - 1d0/eps0T_tr(n)

	str_fac_tr(n) = sum(str_fackt_tr(n,1:nsteps)) 
enddo

 

!-------------------------------------------------------------------------------
!--------------------------  calculate Im{chi(k,w)}----------------------------
!-------------------------------------------------------------------------------
!set up frequencies 
write(*,*) "calculating Im{chi(k,w)}"

Nw = 200
allocate(chikw(Nk,Nw))
allocate(omegas(Nw))

max_freq = 1d0/timestep 
write(*,*) "max_frequency (ps^-1) = ", max_freq
write(*,*) "max_frequency (cm^-1) = ", max_freq*(10.0/3.0)

delta = (Log10(real(max_freq))-4)/Nw !span 4 orders of magnitude

do i = 2, Nw
	omegas(i) = 2d0*pi*(10d0**(dble(Nw - i)*delta))
enddo	
	omegas(1) = 0.0


 chikw = 0d0

do i = 1, Nk
	do w = 1, Nw
		!calculate integral using Trapezoid rule (may be slightly more accurate)
		chikw(i,w) = chikw(i,w) + phiL(i,1)*dcos(omegas(w)*(0)*timestep)/2
		do t = 2, nsteps-1
			chikw(i,w) = chikw(i,w) + phiL(i,t)*dcos(omegas(w)*(t-1)*timestep)
		enddo
		chikw(i,w) = chikw(i,w) + phiL(i,nsteps)*dcos(omegas(w)*(nsteps)*timestep)/2
		chikw(i,w) = chikw(i,w)*timestep
	enddo
	chikw(i,:) = 2d0*chik0(i)*chikw(i,:)
enddo

 
call write_out

call elapsed_time(seconds) 


End Program epskw
