!--------------------------------------------------------------------------
!----- In the interest of maintaining the speed for non-TTM3F  -----------
!----- calculations we have a seperate subroutine for TTM3F --------------
!--------------------------------------------------------------------------

module chi_k
Implicit none 

contains 

!--------------------------------------------------------------------------
!----------------  Longitudinal chi(k) & structure factor calculation ----
!--------------------------------------------------------------------------
subroutine calc_chik
 use main_stuff
 Implicit none 

do n = 1, Nk 
	do ix = 1,3
		tmpOr = 0 
		tmpOc = 0
		tmpHr = 0
		tmpHc = 0 
		!longitudinal part 
		do i = 1, Nmol
			!Oxygens 
			Orp = dcos( k(n)*Oxy(ix,i) )
			Ocp = dsin( k(n)*Oxy(ix,i) ) 

			tmpOr = tmpOr + Orp 
			tmpOc = tmpOc + Ocp
 
			!Hydrogens
			Hrp = dcos( k(n)*Hydro(ix,2*i-0) ) + dcos( k(n)*Hydro(ix,2*i-1) )
			Hcp = dsin( k(n)*Hydro(ix,2*i-0) ) + dsin( k(n)*Hydro(ix,2*i-1) )	

			tmpHr = tmpHr + Hrp
			tmpHc = tmpHc + Hcp	
	
			!self part contribution for this moleucle
			chik0_self(n) = chik0_self(n) +  (qO*Orp + qH*Hrp)**2 + (qO*Ocp + qH*Hcp)**2 
		enddo

		rhokt(n,t) = rhokt(n,t) + dcmplx(qO*tmpOr + qH*tmpHr, qO*tmpOc + qH*tmpHc)

		chik0(n)   = chik0(n)   + (qO*tmpOr + qH*tmpHr)**2 +  (qO*tmpOc + qH*tmpHc)**2

		str_fackt(n,t) = str_fackt(n,t) +  (tmpOr    +    tmpHr)**2 +  (tmpOc    +    tmpHc)**2
	enddo
enddo


end subroutine calc_chik



!--------------------------------------------------------------------------
!----------------  Transverse chi(k) calculation ---------------------------
!--------------------------------------------------------------------------
subroutine calc_chik_transverse
 use main_stuff
 Implicit none 
 real(8), dimension(3) :: rCM, raj
 real(8) :: magraj, q
 double complex, dimension(3) :: mPol !polarization vector for molecule
 double complex, dimension(3) :: Pol !total polarization vector at that k

do n = 1, Nk 

	!transverse part
	do ix = 1,3

		Pol = 0
		do i = 1, Nmol
			rCM = (16d0*Oxy(:,i) +  Hydro(:,2*i)  + Hydro(:,2*i-1))/18d0

			mPol = 0
			do j = 1,3
				if (j .eq. 1) raj = Oxy(:,i) - rCM
				if (j .eq. 2) raj = Hydro(:,2*i-0) - rCM
				if (j .eq. 3) raj = Hydro(:,2*i-1) - rCM
				if (j .eq. 1) q = qO
				if (j .eq. 2) q = qH
				if (j .eq. 3) q = qH

				if (raj(ix) /= 0.0) then
					mPol = mPol - dcmplx(0, 1)*( q*raj/(k(n)*raj(ix)) )*( exp( dcmplx(0, 1)*k(n)*raj(ix) ) - (1d0,1d0) ) 
				endif

			enddo 

			if (Pol(1) /= Pol(1) ) then
				write(*,*) "ERROR in transverse polarization!! NaN"
				write(*,*) "i  ", i
				write(*,*) "j  ", j
				write(*,*) "ix ", ix
				write(*,*) "n ", n
				write(*,*) "Pol = ", Pol
				write(*,*) "mPol = ", mPol
				write(*,*) "product ",  mPol*cmplx( dcos(k(n)*rCM(ix)),  -dsin(k(n)*rCM(ix)) )
			endif
			Pol = Pol + mPol*cmplx( dcos(k(n)*rCM(ix)),  -dsin(k(n)*rCM(ix)) )
		

		enddo
 		
		if (ix .eq. 1) PolTkt(n,t,:) = PolTkt(n,t,:) + cross_product((/ k(n), 0d0,  0d0 /) , Pol)
		if (ix .eq. 2) PolTkt(n,t,:) = PolTkt(n,t,:) + cross_product((/ 0d0 , k(n), 0d0 /) , Pol)
		if (ix .eq. 3) PolTkt(n,t,:) = PolTkt(n,t,:) + cross_product((/ 0d0 , 0d0, k(n) /) , Pol)
	enddo

enddo


end subroutine calc_chik_transverse
 




!--------------------------------------------------------------------------
!----------------  TTM3F chi(k) calculation (variable charges) -----------
!--------------------------------------------------------------------------
subroutine calc_chik_TTM3F


end subroutine calc_chik_TTM3F



!----------------------------------------------------------------------------------
!------------------- function to compute COMPLEX cross product ---------------------------
!----------------------------------------------------------------------------------
function cross_product(x,z)
 Implicit None
 real(8), dimension(3), intent(in) :: x
 double complex, dimension(3), intent(in) :: z 
 double complex, dimension(3) :: cross_product

  cross_product(1) = cmplx(x(2))*z(3) - z(2)*cmplx(x(3)) 
  cross_product(2) = z(1)*cmplx(x(3)) - cmplx(x(1))*z(3)
  cross_product(3) = cmplx(x(1))*z(2) - z(1)*cmplx(x(2))

end function cross_product




end module chi_k
