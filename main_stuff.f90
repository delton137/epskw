!Module to store the Input-Output part of the code and global variables
Module main_stuff
 Implicit none 
real(8),dimension(:,:),allocatable :: Oxy, Hydro, Msites
real(8),dimension(:)  ,allocatable :: k, chik0, omegas, chik0T, eps0T, chik0_self, str_fac, phiTcomponent
real(8),dimension(:,:),allocatable ::  phiL, phiT, chikw, chikT, str_fackt
double complex, dimension(:,:), allocatable :: rhokt 
double complex, dimension(:,:,:), allocatable :: polTkt
real(8), dimension(:),allocatable :: qHs, qOs
real(8),dimension(3) :: v1, v2, v3, summ, box, ibox
real(8) ::  tmpOr, tmpOc, tmpHr, tmpHc, Orp, Ocp, Hrp, Hcp, mink
 character(len=3) :: sym
 character(120)   :: fileinp
 character(120) :: TTM3F_dip_input,fileheader,model
real(8) :: vol,  maxk, qO, qH, qO2, qH2, temp,  qOqH, delta
real(8) :: prefac,  r, seconds, rOM, mag1, timestep, max_freq 
integer :: Na, Nmol, i, j, ia, ix, nsteps, w
integer :: npts, t, n, Nk, ierror, max_num, Nw
logical :: TIP4P, GRIDSAMPLE, TTM3F
real(8), parameter :: pi = 3.14159265d0 
real(8), parameter :: kb = 1.3806488d-23 ! Jul/Kelvin
real(8), parameter :: eps_0 = 8.85418781d-12 ! F/m = C^2/Jul/m
real(8), parameter :: a2m = 1e-10 ! ang to meters
real(8), parameter :: e2C = 1.60217657d-19 !electron charge in C
real(8), parameter :: mole = 6.0221413d+23 !a mole
real(8) :: rHH, rMH
!Variables for reading XTC
 character(len=3) :: filetype    
Integer               :: ID1,STEP,RET, NAT, indx, maxsteps
Integer, parameter   :: MAGIC=1995, MAXAT=1000000001
real(4)                :: XTCBOX(9),DT,TIME,PREC, OLDTIME
real(4), dimension(:),allocatable, save :: X
!Fourier transform stuff
complex, dimension(:), allocatable :: aux1,vcross 
integer :: trun, tread


 contains 

!--------------------------------------------------------------------------
!----------------  Read input file ---------------------------------------
!--------------------------------------------------------------------------
subroutine read_input_file

!----------------------- open input file and read --------------------------
!Open(5,file="eps_k.inp",status="old",action="read",iostat=ierror)
read(5,*) fileinp
read(5,*) filetype 
read(5,*) fileheader
read(5,*) box
read(5,*) maxk
read(5,*) Nk
read(5,*) temp
read(5,*) model
read(5,*) qO
read(5,*) qH
read(5,*) maxsteps
read(5,*) Na
read(5,*) timestep
!close(5)
end subroutine read_input_file

!------------------------------------------------------------------------------
!--------------------------- set up model ------------------------------------
!------------------------------------------------------------------------------
subroutine set_up_model
if (model == 'spce') then
	qO = -.8476
	qH = .4238
	TIP4P = .false. 
	write(*,*) "Model is SPC/E"
else if (model == 'tip3p') then
	qO = -0.834
	qH =  0.417 	
	rMH = 0.9572
	rHH = 1.5139
	TIP4P = .false. 
	write(*,*) "Model is TIP3P"
else if ((model == 'tip4eps') .or. (model == 'TIP4eps')) then
	qO = -1.054
	qH =  0.527
	rOM = .105
	rHH = 1.513900
	rMH = 0.896784
	TIP4P = .true. 
	write(*,*) "Model is TIP4eps"
else if (model == 'tip4p') then
	qO = -1.04
	qH = .52
	rOM = .15
	TIP4P = .true.
	write(*,*) "Model is TIP4P"
else if (model == 'tip4p2005') then
	qO = -1.1128
	qH = .5564
	rOM = .15555
	rHH = 1.513900
	rMH = 0.896784
	TIP4P = .true.
	write(*,*) "Model is TIP4P/2005"
else if (model == 'tip4p2005f') then
	qO = -1.1128	
	qH = .5564
	rOM = .1546
	rHH = 1.513900 !these are TIP4P/2005
	rMH = 0.896784
	TIP4P = .true.
	write(*,*) "Model is TIP4P/2005f"
else 
	write(*,*) "Model is generic 3 site"
	qO = -1	
	qH = .5
	rOM = .1546
	rHH = 1.513900  
	rMH = 0.896784
	TIP4P = .false.
endif 
if ((model == 'ttm3') .or. (model == 'ttm3f')) then
	write(*,*) "Model is TTM3F"
	qO = -1	
	qH = .5
	rOM = 0.4646
	TIP4P = .true.
	TTM3F = .true.
endif
if ( .not. (qO + 2*qH .eq. 0) ) then
	write(*,*) "ERROR in charge values!!"
	stop
endif
end subroutine set_up_model
   
!------------------------------------------------------------------------------
!--------------------------- Open Files --------------------------------------
!------------------------------------------------------------------------------
subroutine open_trajectory_files
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
	Do
		Read(12,*,iostat=ierror) 
		if (ierror /= 0) then
   	    		exit
   		End if
  		npts = npts + 1 !Count number of lines
	enddo
	rewind(12)
	close(12)
	open(12,file=fileinp,status="old",action="read",iostat=ierror)

	nsteps = floor(real(npts) / real((Na + 2)))! Number of timesteps
	write(*,*) "there are ", nsteps, " steps in the file"
	if (maxsteps .gt. nsteps) then
		maxsteps = nsteps
	endif 
		
	Nmol = Na/3
	write(*,*) "Number of molecules: ", Nmol
endif
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

	if (TIP4P) then
		Nmol = Na/4
	else
		Nmol = Na/3
	endif
    	write(*,*) "Reading XTC. If this is a  4 site model we assume all 4 sites are in the XTC."
     	write(*,*) "We assume the units are nm in the .xtc. They will be converted to Ang."
endif 

write(*,*) "number of atoms = ", Na
write(*,*) "number of molecules =", Nmol

if (TIP4P) allocate(Msites(3,Nmol))
allocate(Oxy(3,Nmol))
allocate(Hydro(3,Nmol*2))
if (TTM3F) allocate(qOs(Nmol))
if (TTM3F) allocate(qHs(2*Nmol))


if (TTM3F) then 
	open(50,file=fileinp(1:LEN_TRIM(fileinp)-9)//"chgs.dat",status="old",action="read",iostat=ierror)
	if (ierror /= 0) then
		write(*,*) "ERROR opening TTM3F charge file"
	endif
	!open(51,file=fileinp(1:LEN_TRIM(fileinp)-9)//"dip.dat",status="old",action="read",iostat=ierror)
endif


end subroutine open_trajectory_files



!--------------------------------------------------------------------------
!----------------  Read Trajectory Frame ---------------------------------
!--------------------------------------------------------------------------
subroutine read_trajectory_frame
 !----------- Reading xtc 
  if (filetype .eq. "xtc") then

	if (t .gt. 1) then
        	CALL READXTC(ID1,NAT,STEP,TIME,XTCBOX,X,PREC,RET)
	endif 

	!move coords from the X(:) array to the Oxy & Hydro arrays
	indx = 1
	!four site model case - the m sites are in the .xtc, these are considered the "oxy"
	if (TIP4P .eq. .true.) then
		do ia = 1, Nmol
			Oxy(:,ia)   	= X(indx+0:indx+2) 
			Hydro(:,2*ia-0) = X(indx+3:indx+5)   
			Hydro(:,2*ia-1) = X(indx+6:indx+8) 
			Oxy(:,ia)    = X(indx+9:indx+11)  
			indx = indx+12	
		enddo
	else
	!three site model case
		do ia = 1, Nmol
			Hydro(:,2*ia-0) = X(indx+0:indx+2)   
			Hydro(:,2*ia-1) = X(indx+3:indx+5)   
			Oxy(:,ia)       = X(indx+6:indx+8)  
			indx = ia*9  	
		enddo
	endif

	!Convert to Angstroms for consistancy
	Hydro = 10*Hydro
	Oxy = 10*Oxy
	if (TIP4P) Msites = 10*Msites
  endif
  !----------- Reading xyz------------------------
  if (filetype .eq. "xyz") then
        read(12,*)
	read(12,*)
	Do ia = 1, Nmol
		Read (12,*) sym, (Oxy(ix,ia),ix=1,3)
		Read (12,*) sym, (Hydro(ix,2*ia-0), ix=1,3)
		Read (12,*) sym, (Hydro(ix,2*ia-1), ix=1,3)
	enddo
	!If TI4P-like model than then negative charge ("Oxy") is m-site positions
	!These are typically not included in .xyz so we compute them here.
  	if (TIP4P .or. TTM3F) then
  		Do ia = 1, Nmol
			v1=Hydro(:,2*ia-0)-Oxy(:,ia)
   			v1=v1 - box*anint(v1*ibox)!PBC
			v2=Hydro(:,2*ia-1)-Oxy(:,ia)
       			v2=v2 - box*anint(v2*ibox)!PBC
			summ = v1 + v2  
			v3=(summ/sqrt(dot_product(summ,summ)))*rOM!vector from O to m-site
			Msites(:,ia) = Oxy(:,ia) + v3 !shift the O position to they m-site
		enddo
	endif
  endif
 if (TTM3F) then 
	do i = 1, Nmol
		read(50,*) qOs(i)
		read(50,*) qHs(2*i-0)
		read(50,*) qHs(2*i-1)
	enddo
 endif

end subroutine read_trajectory_frame

!-------------------------------------------------------------------------------
!--------------------------  Write out all the things! -------------------------
!-------------------------------------------------------------------------------
subroutine write_out


!-------------------------------------------------------------------------------
 open(20,file=trim(fileheader)//"_phikL.dat",status="unknown")

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
 write(20,'(a)') '@ yaxis label "\f{Symbol}f\f{Times-Roman}\sL\N(k,t)" '
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
 write(20,'(a,1I2)') '@ legend length ', Nk
do i = 1, Nk
	if (i .lt. 10) then 
 		write(20,'(a,I1,a,1f6.2,a,1f6.2,a)') '@ s', i-1, ' legend "k = ', k(i) ,'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/k(i),'\cE\C" '
	else
 		write(20,'(a,I2,a,1f6.2,a,1f6.2,a)') '@ s', i-1, ' legend "k =', k(i)  ,'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/(i),'\cE\C" '
	endif
enddo

 do t = 1, nsteps/2
	write(20,'(1ES12.3)',advance='no') real(t-1)*timestep 
 	do n = 1, Nk-1
 		write(20,'(1f12.4)',advance='no') phiL(n,t) 
	enddo
	 	write(20,'(1f12.4)',advance='yes') phiL(Nk,t) 
 enddo
 close(20)




!-------------------------------------------------------------------------------
 open(20,file=trim(fileheader)//"_phikT.dat",status="unknown")

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
 write(20,'(a)') '@ yaxis label "\f{Symbol}f\f{Times-Roman}\sT\N(k,t)" '
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
 write(20,'(a,1I2)') '@ legend length ', Nk
do i = 1, Nk
	if (i .lt. 10) then 
 		write(20,'(a,I1,a,1f6.2,a,1f6.2,a)') '@ s', i-1, ' legend "k = ', k(i) ,'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/k(i),'\cE\C" '
	else
 		write(20,'(a,I2,a,1f6.2,a,1f6.2,a)') '@ s', i-1, ' legend "k =', k(i) ,'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/(i),'\cE\C" '
	endif
enddo

 do t = 1, nsteps/2
	write(20,'(1ES12.3)',advance='no') real(t-1)*timestep 
 	do n = 1, Nk-1
 		write(20,'(1f12.4)',advance='no') phiT(n,t) 
	enddo
	 	write(20,'(1f12.4)',advance='yes') phiT(Nk,t) 
 enddo
 close(20)



!-------------------------------------------------------------------------------
 open(21,file=trim(fileheader)//"_chik.dat",status="unknown")
write(21,'(a)') '# This .xvg is formated for xmgrace '
write(21,'(a)') '@ xaxis  label "k(\cE\C\S-1\N)" '
write(21,'(a)') '@ yaxis  label "\f{Symbol}c\f{Times-Roman}(k,0)" '
write(21,'(a)') '@ TYPE xy '
write(21,'(a)') '@ legend on '
write(21,'(a)') '@ legend box off '
write(21,'(a)') '@ legend loctype view '
write(21,'(a)') '@ legend 0.78, 0.8'
write(21,'(a)') '@ legend length 2'
write(21,'(a)') '@ s0 legend \" ", "\"" '
 do i = 1, Nk
 	write(21,'(1f10.4,2f16.4)')  k(i), chik0(i), chik0_self(n)
 enddo
 close(21)

!-------------------------------------------------------------------------------
 open(17,file=trim(fileheader)//"_epsk.dat",status="unknown")
 write(17,'(a)') '# This .xvg is formated for xmgrace "'
 write(17,'(a)') '@ xaxis  label "k(\cE\C\S-1\N)" '
 write(17,'(a)') '@ yaxis  label "\f{Symbol}e\f{Times-Roman}(k,\f{Symbol}w\f{Times-Roman})" '
 write(17,'(a)') '@ TYPE xydy '
 write(17,'(a)') '@ legend on '
 write(17,'(a)') '@ legend box off '
 write(17,'(a)') '@ legend loctype view '
 write(17,'(a)') '@ legend 0.78, 0.8'
 write(17,'(a)') '@ legend length 2'
 write(17,'(a)') '@ s0 legend \" ", "\"" '
 do i = 1, Nk
 	write(17,'(1f10.4,2f16.4)')  k(i), 1/(1-chik0(i))
 enddo
 close(17)

!-------------------------------------------------------------------------------
 open(18,file=trim(fileheader)//"_epskT.dat",status="unknown")
 write(18,'(a)') '# This .xvg is formated for xmgrace "'
 write(18,'(a)') '@ xaxis  label "k(\cE\C\S-1\N)" '
 write(18,'(a)') '@ yaxis  label "\f{Symbol}e\f{Times-Roman}(k,\f{Symbol}w\f{Times-Roman})" '
 write(18,'(a)') '@ TYPE nxy '
 write(18,'(a)') '@ legend on '
 write(18,'(a)') '@ legend box off '
 write(18,'(a)') '@ legend loctype view '
 write(18,'(a)') '@ legend 0.78, 0.8'
 write(18,'(a)') '@ legend length 2'
 write(18,'(a)') '@ s0 legend \" ", "\"" '
 do i = 1, Nk
 	write(18,'(1f10.4,2f16.4)')  k(i), eps0T(i), chik0T(i)
 enddo
 close(18)



!-------------------------------------------------------------------------------
 open(18,file=trim(fileheader)//"_str_fac.dat",status="unknown")
 write(18,'(a)') '# This .xvg is formated for xmgrace "'
 write(18,'(a)') '@ xaxis  label "k(\cE\C\S-1\N)" '
 write(18,'(a)') '@ yaxis  label "\f{Symbol}S\f{Times-Roman}(k,0)" '
 write(18,'(a)') '@ TYPE nxy '
 write(18,'(a)') '@ legend on '
 write(18,'(a)') '@ legend box off '
 write(18,'(a)') '@ legend loctype view '
 write(18,'(a)') '@ legend 0.78, 0.8'
 write(18,'(a)') '@ legend length 2'
 write(18,'(a)') '@ s0 legend \" ", "\"" '
 do i = 1, Nk
 	write(18,'(1f10.4,2f16.4)')  k(i), str_fac(i)
 enddo
 close(18)



!-------------------------------------------------------------------------------
 open(20,file=trim(fileheader)//"_chikw.dat",status="unknown")

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
 write(20,'(a)') '@ xaxis label "\f{Symbol} w \f{Times-Roman}(cm\S-1\N)" '
 write(20,'(a)') '@ yaxis label "Im{ \f{Symbol}c\f{Times-Roman}(k,\f{Symbol}w\f{Times-Roman})}" '
 write(20,'(a)') '@ xaxes scale Logarithmic '
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
 write(20,'(a,1I2)') '@ legend length ', Nk
do i = 1, Nk
	if (i .lt. 10) then 
 		write(20,'(a,I1,a,1f6.2,a,1f6.2,a)') '@ s', i, ' legend "k = ', k(i) ,'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/k(i),'\cE\C" '
	else
 		write(20,'(a,I2,a,1f6.2,a,1f6.2,a)') '@ s', i, ' legend "k =', k(i) ,'\cE\C\S-1\N \f{Symbol}l\f{Times-Roman} =', (2.0*pi)/k(i),'\cE\C" '
	endif
enddo


 do w = 1, Nw
	write(20,'(1ES12.3)',advance='no') omegas(w)*100.0/2.99792458d0 !convert to cm^-1 
 	do n = 1, Nk-1
 		write(20,'(1f15.6)',advance='no') chikw(n,w) 
	enddo
	 	write(20,'(1f15.6)',advance='yes') chikw(Nk,w) 
 enddo
 close(20)

end subroutine write_out


end module main_stuff
