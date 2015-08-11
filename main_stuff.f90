!---------------------------------------------------------------------
!Module to store the global variables and various initialization subroutines
!
!Copyright 2014-2015, Daniel C. Elton
!---------------------------------------------------------------------
Module main_stuff
 Implicit none 
!-- kvec stuff --------------
integer, parameter :: max_num_kvecs=10000
real(8),dimension(3,max_num_kvecs) :: kvec
real(8),dimension(max_num_kvecs) :: mags
real(8),dimension(:), allocatable :: magk
integer, dimension(max_num_kvecs) :: num_this_mag=0
integer :: num_ind_mags
real(8), dimension(3) :: mink
real(8) :: mag1
integer, dimension(3) :: max_num
!-----------------------------
real(8),dimension(:,:),allocatable :: Oxy, Hydro, Msites
real(8),dimension(:)  ,allocatable :: omegas
real(8), dimension(:), allocatable :: phicomponent, phicomponentL, phicomponentT
real(8),dimension(:)  ,allocatable ::  chik0,   chik0_self, chik0_distinct
real(8),dimension(:)  ,allocatable ::  chik0T, chik0T_tr, eps0T_tr, str_fac
real(8),dimension(:)  ,allocatable :: magk_tr, chik0_tr, chik0_self_tr 
real(8),dimension(:)  ,allocatable :: str_fac_tr, chik0_err_tr, str_fac_err
real(8),dimension(:,:),allocatable ::  phiL, phiT, phiL_tr, phiT_tr, chikw, chikwT, chikT, str_fackt, str_fackt_tr
double complex, dimension(:,:), allocatable :: rhokt, rhokt_tr
double complex, dimension(:,:,:), allocatable :: polTkt, polTkt_tr
real(8), dimension(:),allocatable :: qHs, qOs
real(8), dimension(:,:),allocatable :: Pdip
real(8), dimension(3) :: v1, v2, v3, summ, box, ibox
real(8) ::    junkmag, muL
 character(len=3) :: sym
 character(120)   :: fileinp
 character(120) :: TTM3F_dip_input,TTM3F_input,fileheader,model
real(8) :: vol,  maxk,     temp,  qOqH, delta
real(8) :: prefac,  r, seconds, rOM, timestep, max_freq,  max_diag
integer :: Na, Nmol, i, j, k, l, ia, ix, nsteps, nsteps_out, w
integer :: npts, t, n, Nk, ierror,  Nw, num_face_diagonals, num_body_diagonals
logical :: TIP4P=.false., TTM3F, SMALLKSET
logical :: CHECK_TRAJECTORY_LENGTH, DYNAMIC_STR_FAC, ALT_CALC
real(8), parameter :: pi = 3.14159265d0 
real(8), parameter :: kb = 1.3806488d-23 ! Jul/Kelvin
real(8), parameter :: eps_0 = 8.85418781d-12 ! F/m = C^2/Jul/m
real(8), parameter :: a2m = 1e-10 ! ang to meters
real(8), parameter :: e2C = 1.60217657d-19 !electron charge in C
real(8), parameter :: mole = 6.0221413d+23 !a mole
real(8), parameter :: Cspeed=3.00d10 ! cm/s
real(8), parameter :: hbar=6.626d-34 ! J*s
real(8), parameter :: ps2s=1d-12! 1fs in s
real(8) :: rHH, rMH
!Variables for reading XTC
 character(len=3) :: filetype    
Integer               :: ID1,STEP,RET, NAT, indx, maxsteps
Integer, parameter   :: MAGIC=1995, MAXAT=1000000001
real(4)                :: XTCBOX(9),DT,TIME,PREC,OLDTIME
real(4), dimension(:),allocatable, save :: X
!Fourier transform stuff
complex, dimension(:), allocatable :: aux1,vcross 
integer :: trun, tread
!Setting up model / reading .xtc stuff
integer :: AtomsPerMol
real(8), dimension(:), allocatable :: qs
real(4), dimension(:,:,:), allocatable :: atoms
!distance-dependent stuff 
complex, dimension(:,:), allocatable :: spheresL 
complex, dimension(:,:,:), allocatable :: spheresT
complex, dimension(:,:), allocatable :: mPolsL
complex, dimension(:,:,:), allocatable :: mPolsT
logical :: DISTDEP, IRCALC, DIPSPHERE, SPHERESPHERE, LIMIT_CALC, K_EQ_0_DIST_DEP
integer :: Ncalc, Nskip
integer :: Nr, countR
real, dimension(:), allocatable :: numR
real(8) :: distance, delta_R, dist
double precision, dimension(3) :: tmp3, halfbox
real, dimension(:,:,:), allocatable :: rCMs
real, dimension(3) :: Rij

 contains 

!------------------------------------------------------------------------------
!--------------------------- set up model ------------------------------------
!------------------------------------------------------------------------------
subroutine set_up_model

if (model == 'spce') then
	write(*,*) "Model is SPC/E"
	AtomsPerMol = 3
	allocate(qs(AtomsPerMol))
 	qs(1)     = -.8476d0
	qs([2,3]) = .4238d0
 
else if (model == 'tip3p') then
	write(*,*) "Model is TIP3P"
	AtomsPerMol = 3
	allocate(qs(AtomsPerMol))
 	qs(1)     = -0.834d0
	qs([2,3]) =  0.417d0 	

else if ((model == 'tip4eps') .or. (model == 'TIP4eps')) then
	write(*,*) "Model is TIP4eps"
	AtomsPerMol = 3
	allocate(qs(AtomsPerMol))
	qs(1)     = -1.054d0
	qs([2,3]) =  0.527d0
	TIP4P = .true. 

else if (model == 'tip4p') then
	write(*,*) "Model is TIP4P"
	AtomsPerMol = 3
	allocate(qs(AtomsPerMol))
	qs(1) = -1.04d0
	qs([2,3]) = .52d0
	TIP4P = .true.

else if (model == 'tip4p2005') then
	write(*,*) "Model is TIP4P/2005"
	AtomsPerMol = 3
	allocate(qs(AtomsPerMol))
	qs(1)     = -1.1128d0
	qs([2,3]) = .5564d0
	TIP4P = .true.

else if (model == 'tip4p2005f') then
	write(*,*) "Model is TIP4P/2005f"
	AtomsPerMol = 3
	allocate(qs(AtomsPerMol))
	qs(1)     = -1.1128d0	
	qs([2,3]) = .5564d0
	TIP4P = .true.

else if ((model == 'ttm3') .or. (model == 'ttm3f')) then
	write(*,*) "Model is TTM3F"
	AtomsPerMol = 3
	allocate(qs(AtomsPerMol))
	qs(1)     = -1	
	qs([2,3]) = .5
	TIP4P = .true.
	TTM3F = .true.

else if (model == 'methanol') then
	write(*,*) "Model is methanol with GAFF parameters"
	AtomsPerMol = 6
	allocate(qs(AtomsPerMol))
	qs(1)  = .23940 !C atom
	qs(2)  = .0035  !CH1 atom
	qs(3)  = .0035  !CH2 atom
	qs(4)  = .0035  !CH3 atom
	qs(5)  =-.6676  !O atom
	qs(6)  =0.4177  !OH atom

else if (model == 'methanolH1') then
	write(*,*) "Model is methanol with H1+3 parameters"
	AtomsPerMol = 3
	allocate(qs(AtomsPerMol))
	qs(1)  = .297 !C atom
	qs(2)  = -.728!O atom
	qs(3)  = .431 !OH hydrogen

else if (model == 'acetonitrile') then
	write(*,*) "Model is acetonitrile with GAFF forcefield parameters"
	AtomsPerMol = 6
	allocate(qs(AtomsPerMol))
	qs(1)  = -.4008 !C3 
	qs(2)  = .15640 !H1 
	qs(3)  = .15640 !H2 	
	qs(4)  = .15640 !H3 	
	qs(5)  = .44840 !C2 	
	qs(6)  = -.5168 !N1 

else 
	write(*,*) "Model is generic 3 site water model"
	AtomsPerMol = 3
	allocate(qs(AtomsPerMol))
	qs(1)     = -1d0	
	qs([2,3]) = .5d0
	TIP4P = .false.
endif 
if ( .not. (sum(qs) .eq. 0.0) ) then
	write(*,*) "ERROR in charge values!! - sum of charges not equal to zero = ", sum(qs)
	stop
endif

if (DYNAMIC_STR_FAC) then
	qs = 1
endif

end subroutine set_up_model
   

!------------------------------------------------------------------------------
!--------------------------- set up k vectors --------------------------------
!--------- because of PBCs, k values must be multiples of mink --------------
!------------------------------------------------------------------------------
subroutine setup_k_vectors_and_allocate

 ibox = 1d0/box
 mink = 2d0*pi/box
 max_num = floor(maxk/mink)


 write(*,*) "minimum k's in each direction:"
 write(*,*) mink
 write(*,*) "maximum number of k's along each edge:"
 write(*,*) max_num 
 write(*,*) "maximum number of k's along edges", sum(max_num)

if (SMALLKSET) then

	Nk = 21*3

	n = 1
	do ix = 1,3
		kvec(ix,n+0)  =  1d0*mink (ix)
		kvec(ix,n+1)  =  2d0*mink (ix)
		kvec(ix,n+2)  =  5d0*mink (ix)
		kvec(ix,n+3)  =  nint(0.75d0/mink (ix))*mink (ix)  
		kvec(ix,n+4)  =  nint(1.05d0/mink (ix))*mink (ix) 
		kvec(ix,n+5)  =  nint(1.35d0/mink (ix))*mink (ix) 
		kvec(ix,n+6)  =  nint(1.65d0/mink (ix))*mink (ix) 
		kvec(ix,n+7)  =  nint(1.8d0/mink (ix))*mink (ix) 
		kvec(ix,n+8)  =  nint(1.95d0/mink (ix))*mink (ix) 
		kvec(ix,n+9)  =  nint(2.1d0/mink (ix))*mink (ix) 
		kvec(ix,n+10)  =  nint(2.25d0/mink (ix))*mink (ix) 
		kvec(ix,n+11)  =  nint(2.55d0/mink (ix))*mink (ix) 
		kvec(ix,n+12)  =  nint(3.0d0/mink (ix))*mink (ix) 
		kvec(ix,n+13)  =  nint(4.05d0/mink (ix))*mink (ix) 
		kvec(ix,n+14) =  nint(5.1d0/mink (ix))*mink (ix) 
		kvec(ix,n+15) =  nint(6.15d0/mink (ix))*mink (ix) 
		kvec(ix,n+16) =  nint(7.2d0/mink (ix))*mink (ix) 
		kvec(ix,n+17) =  nint(8.25d0/mink (ix))*mink (ix) 
		kvec(ix,n+18) =  nint(9.3d0/mink (ix))*mink (ix) 
		kvec(ix,n+19) =  nint(10.35d0/mink (ix))*mink (ix) 
		kvec(ix,n+20) =  nint(11.4d0/mink (ix))*mink (ix) 

		n = n + 21
	enddo

     mags(1:Nk) = sum(kvec(:,:),1)

      !add some of the small diagonals to the small k - set ( the face diagonals are especially great, you get 12 different directions!) 
     ! do i = 0,1
!	do j = 0,1
!		do k = 0,1
!			if ( i+j+k .gt. 1) then
!				mag1 = dsqrt( (i*mink (1))**2 + (j*mink (2))**2 + (k*mink (3))**2 )
!				kvec(:,n) = (/  i*mink (1), j*mink (2), k*mink (3) /)
!				mags(n) = mag1
!				Nk = Nk + 1
!			endif
!		enddo
!	enddo
 !    enddo
	
else
 !k vectors parallel to box edges 
 n = 1
 do ix = 1,3
 	do i = 1, max_num(ix)
                if ((i*mink (ix) .gt. 10) .and. (i*mink (ix) .lt. 15)) then
                        !between  k = 4 and k=5 , k points become more sparse
                        if (mod(i,1) .eq. 0) then       
                                kvec(ix,n) =  i*mink (ix)
                                mags(n) = i*mink (ix)
                                n = n + 1
                        endif
                endif
                if (i*mink (ix) .gt. 15) then
                        !at greater than k = 5 , k points become more sparse
                        if (mod(i,2) .eq. 0) then
                                kvec(ix,n) =  i*mink (ix)
                                mags(n) = i*mink (ix)
                !               write(*,*) "adding", mags(n)
                                n = n + 1
                        endif
                else
			kvec(ix,n) =  i*mink (ix)
			mags(n) = i*mink (ix)
			n = n + 1
		endif
	enddo 
 enddo
 write(*,*) "Using ", n-1, "k vectors parallel to the box edges"


if (max_diag .gt. 0) then
 write(*,*) "Building diagonals.."
 do i = 0,nint(5d0/mink (1))
	do j = 0,nint(5d0/mink (2))
		do k = 0,nint(5d0/mink (3))
			if ( i+j+k .gt. 1) then
				if ( (mag1 .lt. maxk) .and. (i+j+k .ne. 0 ) ) then
				mag1 = dsqrt( (i*mink (1))**2 + (j*mink (2))**2 + (k*mink (3))**2 )
					kvec(:,n) = (/  i*mink (1), j*mink (2), k*mink (3) /)
					mags(n) = mag1
					!!     write(*,*) "adding diag ", mags(n)
					n = n + 1

					if (n .gt. max_num_kvecs-1) then
						write(*,*) "limit of 10,000 k vecs reached"
						write(*,*) "code will be too slow/run out of memory if more are used"
						stop
					endif
				endif
			endif
		enddo
	enddo
 enddo
else
!!Construct only k vectors that lie along the faces that are of the form (l,l,0), (l,0,l) or (0,l,l) , where l is an integer
 write(*,*) "Building face diagonals.."
 do l = 1, floor(real(num_face_diagonals/3))
    do i = 0,l,l
       do j = 0,l,l
          do k = 0,l,l
              if  (i+j+k .gt. l ) then
		if ((i .eq. 0) .or. (j .eq. 0) .or. (k .eq. 0)) then 
                    mag1 = dsqrt( (i*mink(1))**2 +(j*mink(2))**2 + (k*mink(3))**2 )
		    if  (mag1 .lt. maxk) then 
                 	    kvec(:,n) = (/  i*mink(1), j*mink(2),k*mink(3) /)
                	    mags(n) = mag1
                             !  write(*,*) "adding face diag ", kvec(:,n)
                	    n = n + 1

		 	    if (n .gt. max_num_kvecs-1) then
				write(*,*) "limit of 10,000 k vecs reached"
				write(*,*) "code will be too slow/run out of memory if more are used"
				stop
			    endif
	 	    endif
	         endif
              endif
          enddo
       enddo
    enddo
 enddo
 do l = 1, floor(real(num_body_diagonals))
                   mag1 = dsqrt( (l*mink(1))**2 +(l*mink(2))**2 + (l*mink(3))**2 )
		    if  (mag1 .lt. maxk) then 
                 	    kvec(:,n) = (/  l*mink(1), l*mink(2),l*mink(3) /)
                	    mags(n) = mag1
                             !  write(*,*) "adding face diag ", kvec(:,n)
                	    n = n + 1
		 	    if (n .gt. max_num_kvecs-1) then
				write(*,*) "limit of 10,000 k vecs reached"
				write(*,*) "code will be too slow/run out of memory if more are used"
				stop
			    endif
	 	    endif
 enddo


endif!if (max_diag .gt. 0) 

 Nk = n - 2

endif! (SMALLKSET)

 write(*,*) "Total number of k vectors (including diagonals) = ", Nk

 allocate(magk (Nk))
 allocate(chik0(Nk))
 allocate(chik0_self(Nk))
 allocate(str_fac(Nk))
 allocate(rhokt(Nk,maxsteps))
 allocate(chikT(Nk,maxsteps))
 allocate(str_fackt(Nk,maxsteps))

 allocate(polTkt(Nk,maxsteps,3))

 chik0_self = 0 
 chik0 = 0 
 str_fac = 0 
 str_fackt = 0 
 rhokt = 0 

 magk = mags(1:Nk)

 call Bubble_Sort(magk, kvec, Nk, max_num_kvecs)



end subroutine setup_k_vectors_and_allocate


!------------------------------------------------------------------------------
!------------- setup all variables relevant to distance-dependence -----------
!------------- by default we use only one k vector in this case --------------
!------------------------------------------------------------------------------
subroutine setup_dist_dep

 if (.not.(LIMIT_CALC)) Ncalc = Nmol
 Nskip = floor(real(Nmol)/real(Ncalc))
 write(*,*) "Number of molecules used in calculation = ", Ncalc
 allocate(mPolsL(Ncalc, maxsteps))
 allocate(mPolsT(Ncalc, 3, maxsteps))
 allocate(rCMs(Ncalc,3,maxsteps))

 ibox = 1d0/box
 mink = 2d0*pi/box
 max_num = floor(maxk/mink)

 write(*,*) "minimum k in each direction:"
 write(*,*) mink
 
Nk = 1

kvec(1,1)  =  1d0*mink(1)
kvec(1,2)  =  0
kvec(1,3)  =  0

mags(1) = mink(1) 

write(*,*) "Total number of k vectors (including diagonals) = ", Nk
 
halfbox = box / 2.0d0

delta_R = 2d0 

Nr = floor(.5d0*sqrt(box(1)**2 + box(2)**2 + box(3)**2)/delta_R) 

write(*,*) "number of different R =", Nr
write(*,*) "spacing bewtween different R = ", delta_R, " Ang"

 allocate(magk(Nk))

 allocate(chik0(Nr))
 allocate(chik0T(Nr))
 allocate(spheresL(Nr, maxsteps))
 allocate(spheresT(Nr, 3, maxsteps))

 write(*,*) "Estimated memory usage for storing dipoles : ", (4*2*Ncalc*maxsteps + 4*2*Ncalc*maxsteps*3)/1000000.0, "Mb"
 write(*,*) "Estimated memory usage for storing spheres : ", (4*2*maxsteps*Nr + 4*2*maxsteps*3*Nr)/1000000.0, " Mb"
 write(*,*) "Estimated memory usage for storing COM     : ", (4*2*Ncalc*maxsteps*3)/1000000.0, " Mb"
 
 magk = mags(1:Nk)

end subroutine setup_dist_dep



!-------------------------------------------------------------------------------
!-------  bubble sort routine for sorting k-vectors ---------------------------
!-------------------------------------------------------------------------------
SUBROUTINE Bubble_Sort(magk, a, Nk, max_num_kvecs)
 Implicit None 
  Integer, intent(in) :: Nk, max_num_kvecs
  REAL(8), INTENT(inout), DIMENSION(Nk) :: magk
  real(8), intent(inout), Dimension(3,max_num_kvecs) :: a
  REAL(8) :: temp
  REAL(8), dimension(3) :: temp2
  INTEGER :: i, j
  LOGICAL :: swapped = .TRUE.
 
  DO j = SIZE(magk)-1, 1, -1
    swapped = .FALSE.
    DO i = 1, j
      IF (magk (i) > magk (i+1)) THEN

        temp = magk (i)
        magk (i) = magk (i+1)
        magk (i+1) = temp

        temp2 = a(:,i)
        a(:,i) = a(:,i+1)
        a(:,i+1) = temp2

        swapped = .TRUE.
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END SUBROUTINE Bubble_Sort

!-------------------------------------------------------------------------------
!-------------  calculate Im{chi_L(k,w)} and Im{chi_T(k,w)}  -------------------
!-------------  via a direct integration of the Fourier integral ---------------
!------------- THIS PART DOES NOT WORK  ----------------------------------------
!------------- (I have a MATLAB routine which performs this function) ----------
!-------------------------------------------------------------------------------

!Subroutine calc_Imagkw
 !Implicit none 
 !real(8) :: avgMag, omega
 !Integer :: NumAvgOver, PointsAvailable, MaxFreqOut
 !complex , dimension(nsteps) :: aux1

 !set up frequencies 
 !write(*,*) "calculating Im{chi(k,w)}"

 !Nw = 2000
 !allocate(chikw(Nk,Nw))
 !allocate(omegas(Nw))
 !allocate(chikwT(Nk,Nw))

 !max_freq = 1d0/timestep 
 !write(*,*) "max_frequency (ps^-1) = ", max_freq
 !write(*,*) "max_frequency (cm^-1) = ", max_freq*(10.0/3.0)

 !PointsAvailable = nsteps	
 !write(*,*) "points available = ", PointsAvailable

 !if (PointsAvailable .lt. Nw) Nw = PointsAvailable
		
 !NumAvgOver = PointsAvailable/Nw 

 !write(*,*) "Averaging over", NumAvgOver

 !do n = 1, Nk 
	!Fourier transform the ACF
!	aux1 = cmplx(phiL_tr(n,:))
!	call four1(aux1,nsteps,-1)

	!smoothing (simple window average) to remove noise
!	Do t = 0, Nw-1
!
!		avgMag = 0 
 ! 		do i = 1, NumAvgOver
!			omega=2d0*3.1415926d0*( t*NumAvgOver+i )/(timestep*nsteps*ps2s) !get freq in 1/s (Hz)
!			avgMag = avgMag + omega*real(aux1(t*NumAvgOver+i)) 
 !		enddo
!	  	chikw(n, t) = avgMag/real(NumAvgOver)
!		omegas(t) = ( floor((t+.5)*numAvgOver) )/(timestep*nsteps*ps2s)    !get central freq in 1/s (Hz
!		
!	enddo
!	chikw(n,:) = chik0(n)*(chikw(n,:) - 1d0)
! enddo

! omegas = omegas/Cspeed  	! convert frequency to cm-1


! Old method (non-FFT)

!delta = (Log10(real(max_freq))-4)/Nw !span 4 orders of magnitude

!do i = 2, Nw
!	omegas(i) = 2d0*pi*(10d0**(dble(Nw - i)*delta))
!enddo	
!	omegas(1) = 0.0


! chikw = 0d0

!do i = 1, Nk
!	do w = 1, Nw
!		!calculate integral using Trapezoid rule (may be slightly more accurate)
!		chikw(i,w)  = chikw(i,w)  + phiL_tr(i,1)*dcos(omegas(w)*(0)*timestep)/2d0
!		chikwT(i,w) = chikwT(i,w) + phiT_tr(i,1)*dcos(omegas(w)*(0)*timestep)/2d0
		
!		do t = 2, nsteps
!			chikw(i,w)  = chikw(i,w)  + phiL_tr(i,t)*dcos(omegas(w)*(t-1)*timestep)
!			chikwT(i,w) = chikwT(i,w) + phiT_tr(i,t)*dcos(omegas(w)*(t-1)*timestep)
!		enddo

!		chikw(i,w)  = chikw(i,w) + phiL_tr(i,nsteps)*dcos(omegas(w)*(nsteps)*timestep)/2d0
!		chikwT(i,w) = chikw(i,w) + phiL_tr(i,nsteps)*dcos(omegas(w)*(nsteps)*timestep)/2d0

!		chikw(i,w)  = omegas(w)*chikw(i,w)*timestep
!		chikwT(i,w) = omegas(w)*chikwT(i,w)*timestep
!	enddo
!	chikw(i,:) = chik0(i)*chikw(i,:)
!	chikwT(i,:) = chik0(i)*chikwT(i,:)
!enddo

!end subroutine calc_Imagkw





end module main_stuff





