!--------------------------------------------------------------------------
! Module for all of the file reading/writing subroutines
!
! General routines for reading/writing .xtc and .xyz can be found here 
!
! Copyright 2014-2015 Daniel C. Elton <delton 17 at gmail .com> 
! 
! License: The MIT License
!--------------------------------------------------------------------------
module InputOutput
use main_stuff

contains

!--------------------------------------------------------------------------
!----------------  Read input file ---------------------------------------
!--------------------------------------------------------------------------
subroutine read_input_file
Implicit none
 !read from input stream 
 read(5,*) fileinp
 read(5,*) filetype 
 read(5,*) fileheader
 read(5,*) box
 read(5,*) temp
 read(5,*) CHECK_TRAJECTORY_LENGTH
 read(5,*) maxk
 read(5,*) num_face_diagonals
 read(5,*) num_body_diagonals
 read(5,*) max_diag
 read(5,*) model
 read(5,*) maxsteps
 read(5,*) Na
 read(5,*) timestep
 read(5,*) SMALLKSET
 read(5,*) nsteps_out
 read(5,*) ALT_CALC
 read(5,*) DYNAMIC_STR_FAC
 read(5,*) DISTDEP
 read(5,*) DIPSPHERE
 read(5,*) SPHERESPHERE
 read(5,*) K_EQ_0_DIST_DEP

 if ((SPHERESPHERE .eqv. .true.) .and. (DIPSPHERE .eqv. .true.)) then 
	write(*,*) "ERROR IN INPUT FILE - please select either dip-sphere or sphere-sphere"
 endif
 if ((DISTDEP .eqv. .true.) .and. ((SPHERESPHERE .eqv. .false.) .and. (DIPSPHERE .eqv. .false.))) then 
	write(*,*) "ERROR IN INPUT FILE - please select either dip-sphere or sphere-sphere"
 endif
 if ((K_EQ_0_DIST_DEP .eqv. .true.) .and. (DISTDEP .eqv. .false.)) then
	write(*,*) "Since k=0 is only configured for dist dependent, setting it equal to true"
 endif

 LIMIT_CALC = .false.

end subroutine read_input_file

!------------------------------------------------------------------------------
!--------------------------- Open Files --------------------------------------
!------------------------------------------------------------------------------
subroutine open_trajectory_files
Implicit none
if (filetype .eq. "xyz") then
	open(12,file=fileinp,status="old",action="read",iostat=ierror)
 	read(12,*,iostat=ierror) NAT
	rewind(12)

        if ( .not. (NAT .eq. Na)) then
            write (*,*) 'ERROR: Number of Atoms in file is different than specified'
	    write(*,*)  NAT, " neq " , Na
            write(*,*) "using", Na
	     Na = NAT
        endif
	if (CHECK_TRAJECTORY_LENGTH) then
		write(*,*) "Checking length of trajectory.. this can take some time)..."
		Do
			Read(12,*,iostat=ierror) 
			if (ierror /= 0) then
   		    		exit
   			End if
  			npts = npts + 1 !Count number of lines
		enddo
		rewind(12)
		write(*,*) "...done"
		nsteps = floor(real(npts) / real((Na + 2)))! Number of timesteps

		write(*,*) "there are ", nsteps, " steps in the file"
		if (maxsteps .gt. nsteps) then
			maxsteps = nsteps
		endif 
	endif
	
	Nmol = Na/3
	write(*,*) "Number of molecules: ", Nmol
	write(*,*) "Box size (user specified) is :", box
endif
!-----------------------------------------------------------------------
!--------------------  Initial check of xtc  -------------------------- 
!-----------------------------------------------------------------------
if (filetype .eq. "xtc") then 
        allocate(X(3*Na))
     	NAT = Na
        CALL XTCOPEN(ID1,fileinp,"R",MAXAT)
        CALL READXTC(ID1,NAT,STEP,TIME,XTCBOX,X,PREC,RET)

	if ( .not. (NAT .eq. Na)) then
	        write (*,*) 'ERROR: Number of Atoms in file is different than specified'
		write(*,*)  NAT, " neq " , Na
		write(*,*) "please fix and rerun" 
		stop
        endif
	
 	write(*,*) "Box size from file is (nm)", XTCBOX(1),XTCBOX(5),XTCBOX(9)
	write(*,*) "Assuming rectangular box and using these dimensions."
	box(1) = XTCBOX(1)*10
	box(2) = XTCBOX(5)*10
	box(3) = XTCBOX(9)*10

	Nmol = Na/AtomsPerMol
	if (TIP4P) Nmol = Na/4
	if (model == 'methanolH1') Nmol = Na/6

    	write(*,*) "Reading XTC. If this is a  4 site model we assume all 4 sites are in the XTC."
     	write(*,*) "We assume the units are nm in the .xtc. They will be converted to Ang."
endif 


 allocate(atoms(3,AtomsPerMol,Nmol))

 write(*,*) "number of atoms = ", Na
 write(*,*) "atoms per molecule = ", AtomsPerMol
 write(*,*) "number of molecules =", Nmol

if (TTM3F) then 
	allocate(Pdip(3,NmoL))

	open(50,file=fileinp(1:LEN_TRIM(fileinp)-9)//"chgs.dat",status="old",action="read",iostat=ierror)
	if (ierror /= 0) then
		write(*,*) "ERROR opening TTM3F charge file"
	endif
	open(51,file=fileinp(1:LEN_TRIM(fileinp)-9)//"Edip.dat",status="old",action="read",iostat=ierror)
	if (ierror /= 0) then
		write(*,*) "ERROR opening TTM3F dipoles file"
	endif
endif

end subroutine open_trajectory_files

!--------------------------------------------------------------------------
!----------------  Read Trajectory Frame ---------------------------------
!--------------------------------------------------------------------------
subroutine read_trajectory_frame
Implicit none
 !----------- Reading xtc 
  if (filetype .eq. "xtc") then

	if (t .gt. 1) then
        	CALL READXTC(ID1,NAT,STEP,TIME,XTCBOX,X,PREC,RET)
		OLDTIME = time
	endif 
	if (t .eq. 2) then
		write(*,*) "measured timestep is ", OLDTIME - TIME , " ps" 
	endif 

	indx = 1 
	!four site model special case - the m sites are in the .xtc, these are considered the "oxy"
	!the first coordinates in each molecule (for Oxygen) are skipped 
	if (TIP4P .eqv. .true.) then
		do j = 1, Nmol
			atoms(:,2,j)  = X(indx+3:indx+5)
			atoms(:,3,j)  = X(indx+6:indx+8) 
			atoms(:,1,j)  = X(indx+9:indx+11)
			indx = indx + 12
		enddo
	!methanol H1 model special case
	!the 3 dummy hydrogens are skipped
	else if (model == 'methanolH1') then
		do j = 1, Nmol
			atoms(:,1,j)  = X(indx+0:indx+2)  !C 
			atoms(:,2,j)  = X(indx+12:indx+14)!O
			atoms(:,3,j)  = X(indx+15:indx+17)!H
			indx = indx + 18	
		enddo
	else
		!general case of an arbitrary molecule w charge on each site
		do j = 1, Nmol
 			do i = 1, AtomsPerMol
				atoms(:, i, j) = X(indx+0:indx+2)
				indx = indx + 3
			enddo 
		enddo
	endif

	!Convert to Angstroms for consistancy
	atoms = atoms*10
  endif

  !----------- Reading xyz------------------------
  !--- the xyz support for TIP4P/TTM3F & msites has not been tested
  
  if (filetype .eq. "xyz") then
        read(12,*)
	read(12,*)
	!If TI4P-like model than then negative charge ("Oxy") is m-site positions
	!These are typically not included in .xyz so we compute them here.
	!The Oxygen position is replaced with the Msite position
  	if (TIP4P .or. TTM3F) then
		do i = 1, Nmol
			do j = 1, 3
				read (12,*) sym, (atoms(ix, j, i), ix=1,3)
			enddo
			v1 = atoms(:,2,i)-atoms(:,1,i)
   			v1 = v1 - box*anint(v1*ibox)!PBC
			v2 = atoms(:,3,i)-atoms(:,1,i)
       			v2 = v2 - box*anint(v2*ibox)!PBC
			summ = v1 + v2  
			v3=(summ/sqrt(dot_product(summ,summ)))*rOM!vector from O to m-site
			atoms(:, 1, i) = atoms(:, 1, i) + v3 !shift the O position to the m-site
		enddo
	else
		do i = 1, Nmol
			do j = 1, AtomsPerMol
				read (12,*) sym, (atoms(ix, i, j), ix=1,3)
			enddo
		enddo
	endif
  endif

 if (TTM3F) then 
	do i = 1, Nmol
		read(50,*) qOs(i)
		read(50,*) qHs(2*i-0)
		read(50,*) qHs(2*i-1)
		read(51,*) Pdip(1:3,i), junkmag
	enddo
 endif

end subroutine read_trajectory_frame




!-------------------------------------------------------------------------------
!------------ write out raw phik, chik data -----------------------------------
!-------------------------------------------------------------------------------
subroutine write_out_phi_chik_raw
Implicit none

!-------Phi_L 
if (DYNAMIC_STR_FAC) then
 open(20,file=trim(fileheader)//"_F_L_raw.dat",status="unknown")
else
 open(20,file=trim(fileheader)//"_phikL_raw.dat",status="unknown")
endif

 do t = 1, nsteps_out
	write(20,'(1ES12.3)',advance='no') real(t-1)*timestep 
 	do n = 1, num_ind_mags-1
 		write(20,'(1f15.6)',advance='no') phiL_tr(n,t) 
	enddo
	 	write(20,'(1f15.6)',advance='yes') phiL_tr(num_ind_mags,t) 
 enddo
 close(20)

!---------Phi_T
if (DYNAMIC_STR_FAC) then
 open(20,file=trim(fileheader)//"_F_T_raw.dat",status="unknown")
else
 open(20,file=trim(fileheader)//"_phikT_raw.dat",status="unknown")
endif 

 do t = 1, nsteps_out
	write(20,'(1ES12.3)',advance='no') real(t-1)*timestep 
 	do n = 1, num_ind_mags-1
 		write(20,'(1f15.6)',advance='no') phiT_tr(n,t) 
	enddo
	 	write(20,'(1f15.6)',advance='yes') phiT_tr(num_ind_mags,t) 
 enddo
 close(20)

!--------Chik_L
if (DYNAMIC_STR_FAC) then
 open(21,file=trim(fileheader)//"_atom_str_fac_raw.dat",status="unknown")
else
 open(21,file=trim(fileheader)//"_chik_raw.dat",status="unknown")
endif

 do i = 1, num_ind_mags
 	write(21,'(1f10.4,1f16.4)')  magk_tr(i), chik0_tr(i)
 enddo
 close(21)


!---------------chi_T
if (.not. DYNAMIC_STR_FAC) then
 open(18,file=trim(fileheader)//"_chikT.dat",status="unknown")
 do i = 1, num_ind_mags
 	write(18,'(1f10.4,2f16.4)')  magk_tr(i), chik0T_tr(i),  eps0T_tr(i)
 enddo
 close(18)
endif

!----------Chik_L_err
if (DYNAMIC_STR_FAC) then
 open(21,file=trim(fileheader)//"_atom_str_fac_err.dat",status="unknown")
else
 open(21,file=trim(fileheader)//"_chik_err.dat",status="unknown")
endif

 do i = 1, num_ind_mags
        write(21,'(1f10.4,2f16.4)')  magk_tr(i), chik0_tr(i), chik0_err_tr(i)
 enddo
 close(21)

end subroutine write_out_phi_chik_raw


!-------------------------------------------------------------------------------
!------------ write out phik, chik data formatted for xmgrace -----------------
!-------------------------------------------------------------------------------
subroutine write_out_phi_chik_xmgrace
Implicit None

!------- Phi_L
if (DYNAMIC_STR_FAC) then
 open(20,file=trim(fileheader)//"_F_L.dat",status="unknown")
else
 open(20,file=trim(fileheader)//"_phikL.dat",status="unknown")
endif

 write(20,'(a)') '# This .xvg is formated for xmgrace "'
 write(20,'(a)') '@map color 0 to (255, 255, 255), "white"  '
 write(20,'(a)') '@map color 1 to (0, 0, 0), "black"  '
 write(20,'(a)') '@map color 2 to (255, 0, 0), "red"  '
 write(20,'(a)') '@map color 3 to (255, 165, 0), "orange" '   
 write(20,'(a)') '@map color 4 to (255, 255, 0), "yellow" '
 write(20,'(a)') '@map color 5 to (0, 139, 0), "green4"   '
 write(20,'(a)') '@map color 6 to (0, 255, 0), "green"    '
 write(20,'(a)') '@map color 7 to (0,   0, 255), "blue"   '
 write(20,'(a)') '@map color 8 to (0, 255, 255), "cyan"    '
 write(20,'(a)') '@map color 9 to (255, 0, 255), "magenta" '
 write(20,'(a)') '@map color 10 to (148, 0, 211), "violet" '
 write(20,'(a)') '@map color 11 to (188, 143, 143), "brown" '
 write(20,'(a)') '@map color 12 to (220, 220, 220), "grey" '
 write(20,'(a)') '@map color 13 to (103, 7, 72), "maroon"  '
 write(20,'(a)') '@map color 14 to (64, 224, 208), "turquoise" '
 write(20,'(a)') '@ xaxis label "t (ps\S-1\N)" '
 write(20,'(a)') '@ yaxis label "\f{Symbol}F\f{Times-Roman}\sL\N(k,t)" '
 write(20,'(a)') '@ xaxes scale Normal '
 write(20,'(a)') '@ yaxes scale Normal '
 write(20,'(a)') '@TYPE xy '
 write(20,'(a)') '@ view 0.100000, 0.150000, 0.900000, 0.850000 '
 write(20,'(a)') '@ legend on '
 write(20,'(a)') '@ legend box off '
 write(20,'(a)') '@ legend loctype view '
 write(20,'(a)') '@ legend 0.93, 0.85'
 write(20,'(a)') '@ legend char size 0.890000 '
 write(20,'(a)') '@ legend color 1  '
 write(20,'(a)') '@ legend length 4 '
 write(20,'(a)') '@ legend vgap 1  '
 write(20,'(a)') '@ legend hgap 0 '
 write(20,'(a,1I2)') '@ legend length ', num_ind_mags
do i = 1, num_ind_mags
	if (i .lt. 11) then 
 		write(20,'(a,I1,a,1f6.2,a,1f6.2,a)') '@ s', i-1, ' legend "k = ', magk_tr(i) , & 
				'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/magk_tr(i),'\cE\C" '
	else
 		write(20,'(a,I2,a,1f6.2,a,1f6.2,a)') '@ s', i-1, ' legend "k =', magk_tr(i)  ,  & 
				'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/(i),'\cE\C" '
	endif
enddo

 do t = 1, nsteps_out
	write(20,'(1ES15.6)',advance='no') real(t-1)*timestep 
 	do n = 1, num_ind_mags-1
 		write(20,'(1f15.6)',advance='no') phiL_tr(n,t) 
	enddo
	 	write(20,'(1f15.6)',advance='yes') phiL_tr(num_ind_mags,t) 
 enddo
 close(20)

!---------Phi_T
if (DYNAMIC_STR_FAC) then
 open(20,file=trim(fileheader)//"_F_T.dat",status="unknown")
else
 open(20,file=trim(fileheader)//"_phikT.dat",status="unknown")
endif

 write(20,'(a)') '# This .xvg is formated for xmgrace "'
 write(20,'(a)') '@map color 0 to (255, 255, 255), "white"  '
 write(20,'(a)') '@map color 1 to (0, 0, 0), "black"  '
 write(20,'(a)') '@map color 2 to (255, 0, 0), "red"  '
 write(20,'(a)') '@map color 3 to (255, 165, 0), "orange" '   
 write(20,'(a)') '@map color 4 to (255, 255, 0), "yellow" '
 write(20,'(a)') '@map color 5 to (0, 139, 0), "green4"   '
 write(20,'(a)') '@map color 6 to (0, 255, 0), "green"    '
 write(20,'(a)') '@map color 7 to (0,   0, 255), "blue"   '
 write(20,'(a)') '@map color 8 to (0, 255, 255), "cyan"    '
 write(20,'(a)') '@map color 9 to (255, 0, 255), "magenta" '
 write(20,'(a)') '@map color 10 to (148, 0, 211), "violet" '
 write(20,'(a)') '@map color 11 to (188, 143, 143), "brown" '
 write(20,'(a)') '@map color 12 to (220, 220, 220), "grey" '
 write(20,'(a)') '@map color 13 to (103, 7, 72), "maroon"  '
 write(20,'(a)') '@map color 14 to (64, 224, 208), "turquoise" '
 write(20,'(a)') '@ xaxis label "t (ps\S-1\N)" '
 write(20,'(a)') '@ yaxis label "\f{Symbol}F\f{Times-Roman}\sT\N(k,t)" '
 write(20,'(a)') '@ xaxes scale Normal '
 write(20,'(a)') '@ yaxes scale Normal '
 write(20,'(a)') '@TYPE xy '
 write(20,'(a)') '@ view 0.100000, 0.150000, 0.900000, 0.850000 '
 write(20,'(a)') '@ legend on '
 write(20,'(a)') '@ legend box off '
 write(20,'(a)') '@ legend loctype view '
 write(20,'(a)') '@ legend 0.93, 0.85'
 write(20,'(a)') '@ legend char size 0.890000 '
 write(20,'(a)') '@ legend color 1  '
 write(20,'(a)') '@ legend length 4 '
 write(20,'(a)') '@ legend vgap 1  '
 write(20,'(a)') '@ legend hgap 0 '
 write(20,'(a,1I2)') '@ legend length ', num_ind_mags
do i = 1, num_ind_mags
	if (i .lt. 11) then 
 		write(20,'(a,I1,a,1f6.2,a,1f6.2,a)') '@ s', i-1, ' legend "k = ', magk_tr(i) , & 
			'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/magk_tr(i),'\cE\C" '
	else
 		write(20,'(a,I2,a,1f6.2,a,1f6.2,a)') '@ s', i-1, ' legend "k =', magk_tr(i) ,  & 
			'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/(i),'\cE\C" '
	endif
enddo

 do t = 1, nsteps_out
	write(20,'(1ES15.6)',advance='no') real(t-1)*timestep 
 	do n = 1, num_ind_mags-1
 		write(20,'(1f15.6)',advance='no') phiT_tr(n,t) 
	enddo
	 	write(20,'(1f15.6)',advance='yes') phiT_tr(num_ind_mags,t) 
 enddo
 close(20)


!---------chik_L
if (DYNAMIC_STR_FAC) then
 open(21,file=trim(fileheader)//"_atom_str_fac.dat",status="unknown")
else
 open(21,file=trim(fileheader)//"_chik.dat",status="unknown")
endif 

write(21,'(a)') '# This .xvg is formated for xmgrace '
write(21,'(a)') '@ xaxis  label "k (\cE\C\S-1\N)" '
write(21,'(a)') '@ yaxis  label "\f{Symbol}c\f{Times-Roman}(k,0)" '
write(21,'(a)') '@ TYPE xy '
write(21,'(a)') '@ legend on '
write(21,'(a)') '@ legend box off '
write(21,'(a)') '@ legend loctype view '
write(21,'(a)') '@ legend 0.78, 0.8'
write(21,'(a)') '@ legend length 2'
write(21,'(a)') '@ s0 legend \" ", "\"" '
 if (ALT_CALC) then 
 	do i = 1, num_ind_mags
 		write(21,'(1f10.4,3f16.4)')  magk_tr(i), chik0_tr(i), chik0_self_tr(i), chik0_tr(i) - chik0_self_tr(i)
 	enddo
 else 
 	do i = 1, num_ind_mags
 		write(21,'(1f10.4,1f16.4)')  magk_tr(i), chik0_tr(i)
 	enddo

 endif
 close(21)


!------------epsk
if (.not. DYNAMIC_STR_FAC) then
 open(17,file=trim(fileheader)//"_epsk.dat",status="unknown")
 write(17,'(a)') '# This .xvg is formated for xmgrace "'
 write(17,'(a)') '@ xaxis  label "k (\cE\C\S-1\N)" '
 write(17,'(a)') '@ yaxis  label "\f{Symbol}e\f{Times-Roman}(k,\f{Symbol}w\f{Times-Roman})" '
 write(17,'(a)') '@ TYPE xy '
 write(17,'(a)') '@ legend on '
 write(17,'(a)') '@ legend box off '
 write(17,'(a)') '@ legend loctype view '
 write(17,'(a)') '@ legend 0.78, 0.8'
 write(17,'(a)') '@ legend length 2'
 write(17,'(a)') '@ s0 legend \" ", "\"" '
 do i = 1, num_ind_mags
 	write(17,'(1f10.4,2f16.4)')  magk_tr(i), 1d0/(1d0 - chik0_tr(i)), chik0_err_tr(i)/(1-chik0_err_tr(i))**2

 enddo
 close(17)
endif
end subroutine write_out_phi_chik_xmgrace

!-------------------------------------------------------------------------------
!------------ write out phik, chik data formatted for xmgrace -----------------
!-------------------------------------------------------------------------------
subroutine write_out_str_fac

 if (DYNAMIC_STR_FAC) then
   open(21,file=trim(fileheader)//"_mol_str_fac.dat",status="unknown")
 else
   open(18,file=trim(fileheader)//"_str_fac.dat",status="unknown")
 endif
 write(18,'(a)') '# This .xvg is formated for xmgrace "'
 write(18,'(a)') '@ xaxis  label "k (\cE\C\S-1\N)" '
 write(18,'(a)') '@ yaxis  label "S\smol\N(k,0)" '
 write(18,'(a)') '@ TYPE xy '
 write(18,'(a)') '@ legend on '
 write(18,'(a)') '@ legend box off '
 write(18,'(a)') '@ legend loctype view '
 write(18,'(a)') '@ legend 0.78, 0.8'
 write(18,'(a)') '@ legend length 2'
 write(18,'(a)') '@ s0 legend \" ", "\"" '
 do i = 1, num_ind_mags
 	write(18,'(1f10.4,2f16.4)')  magk_tr(i), str_fac_tr(i)
 enddo
 close(18)

end subroutine write_out_str_fac




!-------------------------------------------------------------------------------
!------------ write out distance-dependent stuff ------------------------------
!-------------------------------------------------------------------------------
subroutine write_out_dist_dependent
Implicit None

if (.not. (K_EQ_0_DIST_DEP)) then
	!----------chik0L(R)
	 open(21,file=trim(fileheader)//"_chik0R_L.dat",status="unknown")
 	 do i = 1, Nr
 		write(21,'(1f10.4,2f16.9)')  dble(i)*delta_R, chik0(i)
	 enddo
	 close(21)

	!---------PhiL(R)
	 open(20,file=trim(fileheader)//"_phikL_R_raw.dat",status="unknown")

	 do t = 1, nsteps_out
		write(20,'(1ES12.3)',advance='no') real(t-1)*timestep 
 		do n = 1, Nr - 1
 			write(20,'(1f15.9)',advance='no') phiL(n,t) 
		enddo	
	 	write(20,'(1f15.9)',advance='yes') phiL(Nr,t) 
	 enddo
	 close(20)
endif  

!----------chik0T(R) 
 open(21,file=trim(fileheader)//"_chik0R_T.dat",status="unknown")
 do i = 1, Nr
 	write(21,'(1f10.4,2f16.9)')  dble(i)*delta_R, chik0T(i)
 enddo
 close(21)

!-----------PhiT(R)
 open(20,file=trim(fileheader)//"_phikT_R_raw.dat",status="unknown")

 do t = 1, nsteps_out
	write(20,'(1ES12.3)',advance='no') real(t-1)*timestep 
 	do n = 1, Nr - 1
 		write(20,'(1f15.9)',advance='no') phiT(n,t) 
	enddo
	 	write(20,'(1f15.9)',advance='yes') phiT(Nr,t) 
 enddo
 close(20)

end subroutine write_out_dist_dependent





!-------------------------------------------------------------------------------
!------------------------ Extra code to out put Imag(chi(k,w)) ----------------
!-------------------------------------------------------------------------------
!subroutine write_Imagkw
! open(20,file=trim(fileheader)//"_chikwL.dat",status="unknown")
! write(20,'(a)') '# This .xvg is formated for xmgrace "'
!! write(20,'(a)') '@map color 0 to (255, 255, 255), "white"  '
! write(20,'(a)') '@map color 1 to (0, 0, 0), "black"  '
! write(20,'(a)') '@map color 2 to (255, 0, 0), "red"  '
! write(20,'(a)') '@map color 3 to (255, 165, 0), "orange" '   
! write(20,'(a)') '@map color 4 to (255, 255, 0), "yellow" '
! write(20,'(a)') '@map color 5 to (0, 139, 0), "green4"   '
! write(20,'(a)') '@map color 6 to (0, 255, 0), "green"    '
! write(20,'(a)') '@map color 7 to (0,   0, 255), "blue"   '
! write(20,'(a)') '@map color 8 to (0, 255, 255), "cyan"    '
! write(20,'(a)') '@map color 9 to (255, 0, 255), "magenta" '
! write(20,'(a)') '@map color 10 to (148, 0, 211), "violet" '
 !write(20,'(a)') '@map color 11 to (188, 143, 143), "brown" '
! write(20,'(a)') '@map color 12 to (220, 220, 220), "grey" '
! write(20,'(a)') '@map color 13 to (103, 7, 72), "maroon"  '
! write(20,'(a)') '@map color 14 to (64, 224, 208), "turquoise" '
! write(20,'(a)') '@ xaxis label "\f{Symbol} w \f{Times-Roman}(cm\S-1\N)" '
! write(20,'(a)') '@ yaxis label "Im{ \f{Symbol}c\f{Times-Roman}\sL\N(k,\f{Symbol}w\f{Times-Roman})}" '
! write(20,'(a)') '@ xaxes scale Logarithmic '
! write(20,'(a)') '@ yaxes scale Normal '
! write(20,'(a)') '@TYPE xy '
! write(20,'(a)') '@ view 0.100000, 0.150000, 0.900000, 0.850000 '
! write(20,'(a)') '@ legend on '
 !write(20,'(a)') '@ legend box off '
! write(20,'(a)') '@ legend loctype view '
! write(20,'(a)') '@ legend 0.93, 0.85'
! write(20,'(a)') '@ legend char size 0.890000 '
! write(20,'(a)') '@ legend color 1  '
! write(20,'(a)') '@ legend length 4 '
 !write(20,'(a)') '@ legend vgap 1  '
! write(20,'(a)') '@ legend hgap 0 '
! write(20,'(a,1I2)') '@ legend length ', num_ind_mags
!do i = 1, num_ind_mags
!	if (i .lt. 11) then 
!		write(20,'(a,I1,a,1f6.2,a,1f6.2,a)') '@ s', i, ' legend "k = ', magk_tr(i) ,'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/magk_tr(i),'\cE\C" '!
!	else!
! 		write(20,'(a,I2,a,1f6.2,a,1f6.2,a)') '@ s', i, ' legend "k =', magk_tr(i) ,'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/magk_tr(i),'\cE\C" '
!	endif
!enddo
!
!
! do w = 1, Nw
!	write(20,'(1ES12.3)',advance='no') omegas(w)*100.0/2.99792458d0 !convert to cm^-1 
! 	do n = 1, num_ind_mags-1
! 		write(20,'(1f)',advance='no') chikw(n,w) 
!	enddo
!	 	write(20,'(1f)',advance='yes') chikw(num_ind_mags,w) 
! enddo
! close(20)
!end subroutine write_Imagkw

end module InputOutput

