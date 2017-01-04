!--------------------------------------------------------------------------
! subroutines to find k-dependent polarization vectors  
!
! Copyright 2014-2015, Daniel C. Elton <delton 17 at gmail .com> 
! 
! License: The MIT License
!--------------------------------------------------------------------------
module chi_k
Implicit none 

contains 

!--------------------------------------------------------------------------
!----------------Find polarization vectors -------------------------------
!--------------------------------------------------------------------------
subroutine calc_pol_vectors 
 use main_stuff
 Implicit none 
 real(8), dimension(3) :: rCM, raj 
 real(8) :: magraj, q, kdotr
 double complex, dimension(3) :: mPol !polarization vector for molecule
 double complex, dimension(3) :: Pol !total polarization vector at that k


do n = 1, Nk 
	Pol = 0
	do i = 1, Nmol
		!find geometrical charge center of each molecule
		!(not really center of mass)
		!rCM = 0 
		!do j = 1, AtomsPerMol		
		!	rCM = rCM + atoms(:,j,i)
		!enddo
		!rCM = rCM/AtomsPerMol
		!use first atom in system as reference
		rCM = atoms(:,1,i)
		if (DISTDEP) rCMs(i,:,t) = rCM 

		mPol = 0

		!if TTM3F load charges for each atom
		if (TTM3F) then 
			if (j .eq. 1) qs(1) = qOs(i)
			if (j .eq. 2) qs(2) = qHs(2*i-0)
			if (j .eq. 3) qs(3) = qHs(2*i-1)
		endif

		do j = 1, AtomsPerMol
			raj = atoms(:,j,i) - rCM
			raj = raj - box*anint(raj/box)!PBC

			kdotr = dot_product(kvec(:,n),raj)

			if (kdotr /= 0.0) then
				mPol = mPol - dcmplx(0, 1)*( qs(j)*raj/kdotr )*( exp( dcmplx(0, 1)*kdotr ) - 1d0 ) 
			endif

		enddo!j = 1, AtomsPermol

		!add k-dependence 
		kdotr = dot_product(kvec(:,n), rCM)
		mPol = mPol*exp( cmplx(0, -1)*kdotr )

		!if TTM3F add polarization vector (cf Bertolini Tani Mol Phys 75 1065)
		if (TTM3F) mPol = mPol + Pdip(:,i)*0.20819434d0*exp( dcmplx(0, -1)*dot_product(kvec(:,n),atoms(:,1,i)) )

		if (DISTDEP) then 
			mPolsL(i,t) = dot_product(kvec(:,n), mPol) 
			mPolsT(i,:,t) = cross_product(kvec(:,n), mPol) 
		endif

		Pol = Pol + mPol 
	
	enddo! i = 1, Nmol
 		
	if (.not. DISTDEP) then 
		PolTkt(n,t,:) = cross_product(kvec(:,n) , Pol)
		rhokt(n,t)    = dot_product(kvec(:,n) , Pol)
	endif

enddo! n = 1, Nk 

end subroutine calc_pol_vectors


!--------------------------------------------------------------------------
!----------------Find dipole vectors (for dist-dep, k=0 case) ------------
!----------------This only works for water!! -----------------------------
!--------------------------------------------------------------------------
subroutine calc_dip_vectors 
 use main_stuff
 Implicit none 
 double precision, dimension(3) :: d1 !polarization vector for molecule
 
 
do i = 1, Nmol
	rCMs(i,:,t) = atoms(:,1,i)
	!find dipole
	if (TTM3F) then 
		v1 = atoms(:,2,i) - atoms(:,1,i)
   		v1 = v1 - box*anint(v1/box)!PBC
		v2 = atoms(:,3,i) - atoms(:,1,i)
		v2 = v2 - box*anint(v2/box)!PBC
	      	d1 = (v1*qHs(2*i-0)+v2*qHs(2*i-1))
	      	d1 = d1 + Pdip(:,i)*0.20819434d0!convert Debye to eAng
	else 
		!(it is assumed that atoms(:,1,i) is actually the Msite positions here if TIP4P)
		v1 = atoms(:,2,i) - atoms(:,1,i)
   		v1 = v1 - box*anint(v1/box)!PBC
		v2 = atoms(:,3,i) - atoms(:,1,i)
		v2 = v2 - box*anint(v2/box)!PBC
		d1 = (v1*qH+v2*qH)
	endif 

	mPolsT(i, :, t) = d1 !store dipoles here
enddo

end subroutine calc_dip_vectors


!---------------------------------------------------------------------------
!- alternative longitudinal chi(k) & structure factor calculation for H2O -
!---------------------------------------------------------------------------
subroutine calc_chikL_alternate 
 use main_stuff
 Implicit none 
 real(8) :: qRP, qCP, rRP, rCP, molRp, molCP

do n = 1, Nk 
	qRP = 0 
	qCP = 0 
 	!longitudinal part 
	do j = 1, Nmol
		molRP = 0		
		molCP = 0 

		!if TTM3F load charges for each atom and add point dipole
		!(cf Bertolini Tani Mol Phys 75 1065)
		if (TTM3F) then 
			if (j .eq. 1) qs(1) = qOs(j)
			if (j .eq. 2) qs(2) = qHs(2*j-0)
			if (j .eq. 3) qs(3) = qHs(2*j-1)

			muL = dot_product(kvec(:,n),Pdip(:,i))*0.20819434d0  !convert Debye to eAng
			molRP = molRP + muL*dcos( dot_product(kvec(:,n),atoms(:,1,i)) )
			molCP = molcP + muL*dsin( dot_product(kvec(:,n),atoms(:,1,i)) )
		endif

		!molecular longitudinal polarization
		do i = 1, AtomsPerMol
			molRP = molRP + qs(i)*dcos( dot_product(kvec(:,n),atoms(:,i,j)) )
			molCP = molCP + qs(i)*dsin( dot_product(kvec(:,n),atoms(:,i,j)) ) 
		enddo
 
		!self part contribution for this moleucle
		chik0_self(n) = chik0_self(n) +  molRP**2 + molCP**2

		qRP = qRP + molRP
		qCP = qCP + molCP

 	enddo!do j = 1, Nmol 

	rhokt(n,t) = rhokt(n,t) + dcmplx(qRP, qCP)

	chik0(n)   = chik0(n) + qRP**2 + qCP**2

	str_fackt(n,t) = 0 !defunct! str_fackt(n,t) +  (tmpOr    +    tmpHr)**2 +  (tmpOc    +    tmpHc)**2
enddo
end subroutine calc_chikL_alternate 


!----------------------------------------------------------------------------------
!------------------- function to compute COMPLEX cross product ------------------- 
!----------------------------------------------------------------------------------
function cross_product(x,z)
 Implicit None
 double precision, dimension(3), intent(in) :: x
 double complex, dimension(3), intent(in) :: z 
 double complex, dimension(3) :: cross_product

  cross_product(1) = dcmplx(x(2))*z(3) - z(2)*dcmplx(x(3)) 
  cross_product(2) = z(1)*dcmplx(x(3)) - dcmplx(x(1))*z(3)
  cross_product(3) = dcmplx(x(1))*z(2) - z(1)*dcmplx(x(2))

end function cross_product




end module chi_k
