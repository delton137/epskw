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
	use m_timer !for timing	
	use main_stuff
	use correlation_function
	use chi_k
Implicit None

call set_timer
call read_input_file
call set_up_model
call open_trajectory_files

!------------------------------------------------------------------------------
!--------------------------- set up k vectors --------------------------------
!--------- because of PBCs, k values must be multiples of mink --------------
!------------------------------------------------------------------------------
ibox = 1d0/box
mink = minval(2d0*pi*ibox)
max_num = floor(maxk/mink)

write(*,*) "minimum k vector = ", mink
write(*,*) "number of k's = ", max_num

Nk = 6
!Nk = max_num
allocate(chik0(Nk))
allocate(chik0_self(Nk))

allocate(chik0T(Nk))
allocate(eps0T(Nk))

allocate(str_fac(Nk))

allocate(rhokt(Nk,maxsteps))
allocate(chikT(Nk,maxsteps))
allocate(str_fackt(Nk,maxsteps))
allocate(polTkt(Nk,maxsteps,3))


allocate(k(Nk))

!do i = 1, max_num
!	k(i) =  i*mink
!enddo

k(1) = mink
k(2) = 2*mink
k(3) = 3*mink
k(4) = 4*mink
k(5) = 5*mink
k(6) = 6*mink

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
	 	call calc_chik ! Compute chi(k)  
	endif 

	call calc_chik_transverse

	if (mod(t,10) .eq. 0) write(*,*) t  

 	nsteps = nsteps + 1 

enddo 

1000 continue 


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
	call simple_complex_corr_function(rhokt(i,1:nsteps), phiL(i,1:nsteps), nsteps, nsteps)
	!call calc_corr_function(rhokt(i,1:nsteps), phiL(i,1:nsteps), nsteps) 

	do ix = 1,3
		call simple_complex_corr_function(PolTkt(i,1:nsteps,ix), phiTcomponent, nsteps, nsteps)
		!qcall calc_corr_function(PolTkt(i,:,ix), phiTcomponent, nsteps) 
		phiT(i,:) = phiT(i,:) + phiTcomponent
	enddo
enddo



!--------------------------------------------------------------------------------- 
!----------------  Normalization & prefactors ----------------------------------- 
!---------------------------------------------------------------------------------
 write(*,*) "number of steps used: ", nsteps
 vol = box(1)*box(2)*box(3)
 prefac = (e2C**2)/(eps_0*kb*temp*vol*a2m) 

 !prefactors
 PolTkt = prefac*PolTkt/(3d0*dble(nsteps)) 
 chik0 = prefac*chik0/(3d0*dble(nsteps))  
 chik0_self = prefac*chik0_self/(3d0*dble(nsteps)) 
 str_fackt = str_fackt/(3d0*dble(Nmol)*dble(nsteps))  

 !static transverse  
 eps0T  = 1d0 + phiT(:,1) 
 chik0T = 1d0 - 1d0/eps0T

do n = 1, Nk
 	chik0(n) = chik0(n)/(k(n)**2)

	chik0_self(n) = chik0_self(n)/(k(n)**2)

	str_fac(n) = sum(str_fackt(n,1:nsteps)) 

	!2nd normalization of correlation fun
	phiT(n,:) = phiT(n,:)/phiT(n,1)
	phiL(n,:) = phiL(n,:)/phiL(n,1)
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
		chikw(i,w) = chikw(i,w) + phiL(i,0)*dcos(omegas(w)*(0)*timestep)/2
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
