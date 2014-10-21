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
		tmpOr = 0 
		tmpOc = 0
		tmpHr = 0
		tmpHc = 0 
		!longitudinal part 
		do i = 1, Nmol
			!Oxygens 
			Orp = dcos( dot_product(kvec(:,n),Oxy(:,i)) )
			Ocp = dsin( dot_product(kvec(:,n),Oxy(:,i)) ) 

			tmpOr = tmpOr + Orp 
			tmpOc = tmpOc + Ocp
 
			!Hydrogens
			Hrp = dcos( dot_product(kvec(:,n),Hydro(:,2*i-0)) )  & 
			    + dcos( dot_product(kvec(:,n),Hydro(:,2*i-1)) )
			Hcp = dsin( dot_product(kvec(:,n),Hydro(:,2*i-0)) )  &
			    + dsin( dot_product(kvec(:,n),Hydro(:,2*i-1)) )	

			tmpHr = tmpHr + Hrp
			tmpHc = tmpHc + Hcp	
	
			!self part contribution for this moleucle
			chik0_self(n) = chik0_self(n) +  (qO*Orp + qH*Hrp)**2 + (qO*Ocp + qH*Hcp)**2 
		enddo

		rhokt(n,t) = rhokt(n,t) + dcmplx(qO*tmpOr + qH*tmpHr, qO*tmpOc + qH*tmpHc)

		chik0(n)   = chik0(n)   + (qO*tmpOr + qH*tmpHr)**2 +  (qO*tmpOc + qH*tmpHc)**2

		str_fackt(n,t) = str_fackt(n,t) +  (tmpOr    +    tmpHr)**2 +  (tmpOc    +    tmpHc)**2
enddo


end subroutine calc_chik


!--------------------------------------------------------------------------
!----------------  Longitudinal chi(k) & structure factor calculation - TTM3F
!--------------------------------------------------------------------------
subroutine calc_chik_TTM3F
 use main_stuff
 Implicit none 

do n = 1, Nk 
		tmpOr = 0 
		tmpOc = 0
		tmpHr = 0
		tmpHc = 0 
		tmpDr = 0 
		tmpDc = 0 
		tmpHr_nocharge = 0 
		tmpHr_nocharge = 0 
		!longitudinal part 
		do i = 1, Nmol
			!Oxygens 
			Orp = dcos( dot_product(kvec(:,n),Oxy(:,i)) )
			Ocp = dsin( dot_product(kvec(:,n),Oxy(:,i)) ) 

			tmpOr = tmpOr + qOs(i)*Orp 
			tmpOc = tmpOc + qOs(i)*Ocp

			tmpOr_nocharge  = tmpOr + Orp 
			tmpOc_nocharge  = tmpOc + Ocp
 
			!Hydrogens
			!first hydro
			Hrp = dcos( dot_product(kvec(:,n),Hydro(:,2*i-0)) )   
			Hcp = dsin( dot_product(kvec(:,n),Hydro(:,2*i-0)) )  

			!second hydro
			Hrp2 = dcos( dot_product(kvec(:,n),Hydro(:,2*i-1)) )
			Hcp2 = dsin( dot_product(kvec(:,n),Hydro(:,2*i-1)) )	

			tmpHr = tmpHr + qHs(2*i-0)*Hrp + qHs(2*i-1)*Hrp2
			tmpHc = tmpHc + qHs(2*i-0)*Hcp + qHs(2*i-1)*Hcp2

			tmpHr_nocharge  = tmpHr_nocharge  + Hrp + Hrp2
			tmpHc_nocharge  = tmpHc_nocharge  + Hcp + Hcp2

			!contribution of the point dipole (cf Bertolini Tani Mol Phys 75 1065)
			muL = dot_product(kvec(:,n),Pdip(:,i))*0.20819434d0  !convert Debye to eAng

			Drp = muL*dcos( dot_product(kvec(:,n),Msites(:,i)) )
			Dcp = muL*dsin( dot_product(kvec(:,n),Msites(:,i)) )
			
			tmpDr = tmpDr + Drp 
			tmpDc = tmpDc + Dcp

			!self part contribution for this molecule
			chik0_self(n) = chik0_self(n) +  (qOs(i)*Orp  + qHs(2*i-0)*Hrp + qHs(2*i-1)*Hrp2 + Drp )**2 + (qOs(i)*Ocp + qHs(2*i-0)*Hcp + qHs(2*i-1)*Hcp2 + Dcp)**2 
		enddo

		rhokt(n,t) = rhokt(n,t) + dcmplx(tmpOr + tmpHr + tmpDr, tmpOc + tmpHc + tmpDc)

		chik0(n)   = chik0(n)   + (tmpOr + tmpHr + tmpDr)**2 +  (tmpOc + tmpHc + tmpDc)**2

		str_fackt(n,t) = str_fackt(n,t) +  (tmpOr_nocharge + tmpHr_nocharge )**2  +  (tmpOc_nocharge + tmpHc_nocharge)**2
enddo


end subroutine calc_chik_TTM3F



!--------------------------------------------------------------------------
!----------------  Transverse chi(k) calculation ---------------------------
!--------------------------------------------------------------------------
subroutine calc_chik_transverse
 use main_stuff
 Implicit none 
 real(8), dimension(3) :: rCM, raj
 real(8) :: magraj, q, kdotr
 double complex, dimension(3) :: mPol !polarization vector for molecule
 double complex, dimension(3) :: Pol !total polarization vector at that k


do n = 1, Nk 
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

			kdotr = dot_product(kvec(:,n),raj)
			if (kdotr /= 0.0) then
				mPol = mPol - dcmplx(0, 1)*( q*raj/kdotr )*( exp( dcmplx(0, 1)*kdotr ) - 1d0 ) 
			endif

		enddo !j = 1,3

		if (mPol(1) /= mPol(1) ) then
			write(*,*) "ERROR in transverse polarization!! NaN Debug Info:"
			write(*,*) "molecule ", i
			write(*,*) "atom ", j
			write(*,*) "n ", n
			write(*,*) "Pol = ", Pol
			write(*,*) "mPol = ", mPol
		endif

		kdotr = dot_product(kvec(:,n),rCM)

		Pol = Pol + mPol*exp( dcmplx(0, -1)*kdotr )

	enddo! i = 1, Nmol
 		
	PolTkt(n,t,:) = cross_product(kvec(:,n) , Pol)

enddo! n = 1, Nk 


end subroutine calc_chik_transverse






!----------------------------------------------------------------------------------
!------------------- function to compute COMPLEX cross product ------------------- 
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
