!----------------------------------------------------------------------------- 
! epsk(k,w) program
! Computes correlation functions of the charge-charge structure factor
!        phi(k,t) = <rho(k,t)*rho(-k,t)>/<rho(k,0)>^2/k^2 
! It works for an arbitrary monomolecular system described by a point charge model
!
! It outputs the longitudinal susceptibility defined by P(k,w) = \chi(k,w) D(k,w) and eps(k,w) = 1/(1-chi(k,w))
!
! It also outputs the static nonlocal susceptibility chi(k,0) and the static structure factor S(k,0)
!  
! REFERENCE: JCP 109 5 pg 1939
!
! in 2015 the epskwR code was merged into this code 
!
! 2014-2015 Daniel C. Elton 
!----------------------------------------------------------------------------- 
Program epskw
	use m_timer 
	use main_stuff
	use InputOutput
	use chi_k
	use calc_phi
Implicit None

call set_timer
call read_input_file
call set_up_model
call open_trajectory_files

if (.not. DISTDEP)       call setup_k_vectors_and_allocate
if (DISTDEP) call setup_dist_dep

!------------------------------------------------------------------------------
!---------------------------------- main loop --------------------------------
!------------------------------------------------------------------------------
nsteps = 0
 do t = 1, maxsteps
  
	call read_trajectory_frame

	if ((filetype .eq. "xtc") .and. (RET.EQ.0)) then 
		write(*,*) "reached end of file at t = ", t, "length of traj is ", TIME, "ps "
		goto 1000 
	endif 

	if (K_EQ_0_DIST_DEP) then
		call calc_dip_vectors
	else if (ALT_CALC) then
		call calc_chikL_alternate 
	else 
	 	call calc_pol_vectors
	endif

	if (mod(t,10) .eq. 0) write(*,*) t  

 	nsteps = nsteps + 1 
 enddo 

 1000 continue 

 write(*,*) "number of steps used: ", nsteps
 if (nsteps_out .gt. nsteps) nsteps_out = nsteps
 call elapsed_time(seconds) 
 
 deallocate(Oxy)
 deallocate(Hydro) 
 
 if (DISTDEP) then 
	call calc_phikRt
	call write_out_dist_dependent
 else 
	call calc_phikt
	call write_out_phi_chik_raw
	call write_out_phi_chik_xmgrace
 endif
 
 call elapsed_time(seconds) 


End Program epskw
