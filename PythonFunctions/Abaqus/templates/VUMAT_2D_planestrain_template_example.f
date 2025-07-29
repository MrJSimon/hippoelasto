C *******************************************************************
C  VUMAT (2D plane strain)
C     
C  ---------------------------------------------------
C  Authors:  Jamie Simon
C  Date:     2024-05-10
C  E-mail:   ...
C  Source:   
C *******************************************************************
C
C  Strain-energy function:
C
C  *** INPUT FROM PYTHON PROGRAM *** STRAIN ENERGY DEFINITION
C
C  Description: 
C
C *******************************************************************
C
      subroutine vumat(
     & nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     & stepTime, totalTime, dt, cmname, coordMp, charLength,
     & props, density, strainInc, relSpinInc,
     & tempOld, stretchOld, defgradOld, fieldOld,
     & stressOld, stateOld, enerInternOld, enerInelasOld,
     & tempNew, stretchNew, defgradNew, fieldNew,
     & stressNew, stateNew, enerInternNew, enerInelasNew )
C
      INCLUDE 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     & charLength(nblock), strainInc(nblock,ndir+nshr),
     & relSpinInc(nblock,nshr), tempOld(nblock),
     & stretchOld(nblock,ndir+nshr),
     & defgradOld(nblock,ndir+nshr+nshr),
     & fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     & stateOld(nblock,nstatev), enerInternOld(nblock),
     & enerInelasOld(nblock), tempNew(nblock),
     & stretchNew(nblock,ndir+nshr),
     & defgradNew(nblock,ndir+nshr+nshr),
     & fieldNew(nblock,nfieldv),
     & stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     & enerInternNew(nblock), enerInelasNew(nblock)
C
      CHARACTER*80 cmname
C
C     LOCAL VARIABLES
C     ---------------
      REAL*8 C10, C01, C20, D1, E, nu
      REAL*8 G1, k1
      REAL*8 B11, B22, B33, B12, Bbar11, Bbar22, Bbar33, Bbar12
      REAL*8 I1b, I2b, detJ, dWdI1, dWdI2, dWdJ1
      REAL*8 p1, p2, p3, u1
	  REAL*8 trace
	  REAL :: U(3,3)
	  
C
C     MATERIAL PARAMETERS
C     ----------------------------------------------------------------
      C10 = props(1)
      C01 = props(2)
      C20 = props(3)
      D   = props(4)
      E   = props(5)
      nu  = props(6) 
C
C	  COMPUTE SHEAR AND KAPPA MODULUS
C	  ----------------------------------------------------------------
      G1 = E/(2.d0*(1.d0 + nu))
	  k1 = E/(3.d0*(1.d0 - 2.d0*nu))
C
C     ***************************************************************
C     ------------ INITIALIZE MATERIAL AS LINEARLY ELASTIC ----------
C	   sigma_ij = 2*G1*epsilon_ij+(k1 - 2/3*G1) delta_ij*epsilon_kk
C     ***************************************************************
C  
      IF (totalTime.EQ.0.0) THEN
         DO k = 1,nblock
			trace =  strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
            stressNew(k,1) = stressOld(k,1) + 2.d0*G1*strainInc(k,1) + (k1-2.d0/3.d0 * G1) * trace
            stressNew(k,2) = stressOld(k,2) + 2.d0*G1*strainInc(k,2) + (k1-2.d0/3.d0 * G1) * trace
			stressNew(k,3) = stressOld(k,3) + 2.d0*G1*strainInc(k,3) + (k1-2.d0/3.d0 * G1) * trace
            stressNew(k,4) = stressOld(k,4) + 2.d0*G1*strainInc(k,4)
         END DO
C
      ELSE
C
C     ***************************************************************
C     ----------- START LOOP FOR MATERIAL POINT CALCULATIONS --------
C     ***************************************************************
C
      DO k = 1,nblock
C		  
C        CALCULATE LEFT CAUCHY-GREEN STRAIN TENSOR, B^star_ij = U_ij^2 = sum_k=1 U_ik*U_jk
C        ---------------------------------------------------------------------------------------------------------------
         B11 = stretchNew(k,1) * stretchNew(k,1) + stretchNew(k,4) * stretchNew(k,4)
         B22 = stretchNew(k,2) * stretchNew(k,2) + stretchNew(k,4) * stretchNew(k,4)
         B12 = stretchNew(k,1) * stretchNew(k,4) + stretchNew(k,4) * stretchNew(k,2)
		 B33 = stretchNew(k,3) * stretchNew(k,3)
C		  
C        CALCULATE THE RIGHT STRETCH TENSOR U (RECAL THE POLAR DECOMPOSITION F = RU)
C        ---------------------------------------------------------------------------------------------------------------
		 U(1,1) = stretchNew(k,1)
		 U(2,2) = stretchNew(k,2)
		 U(1,2) = stretchNew(k,4)
		 U(2,1) = stretchNew(k,4)
		 U(3,3) = stretchNew(k,3)
		 U(1,3) = 0.d0
		 U(3,1) = 0.d0		 
C		 
C        CALCULATE J = |F| = |U| = det(U) = gamma_ijk*U_1i*U_2j*U_3k, where gamma = Levi-Civita symbol
C        ---------------------------------------------------------------------------------------------------------------
		 detJ = U(3,3)*(U(1,1)*U(2,2) - U(1,2)*U(2,1))
C
C        CALCULATE MODIFIED STRAIN TENSOR, B^starbar_ij = J^(-2/3)*B^star_{ij}
C        ---------------------------------------------------------------------------------------------------------------
         Bbar11 = detJ**(-2.d0/3.d0) * B11
         Bbar22 = detJ**(-2.d0/3.d0) * B22
         Bbar12 = detJ**(-2.d0/3.d0) * B12
		 Bbar33 = detJ**(-2.d0/3.d0) * B33
C
C        CALCULATE FIRST AND SECOND INVARIANT of B^starbar. Please note these are the modified invariants !!!
C        ---------------------------------------------------------------------------------------------------------------
		 I1b = Bbar11 + Bbar22 + Bbar33
		 I2b = (1.d0/2.d0) * ((Bbar11 + Bbar22 + Bbar33)**(2.0d0) - (Bbar11**(2.d0) + Bbar22**(2.d0) + Bbar33**(2.d0)))
C
C        CALCULATE DERIVATIVES OF STRAIN-ENERGY FUNCTION
C        ---------------------------------------------------------------------------------------------------------------
		 dWdI1 = C10 + C20*(2.0d0*I1b - 6.0d0)
		 dWdI2 = C01
		 dWdJ  = 1d0/D*(2.0d0*detJ - 2.0d0)
C
C        CALCULATE THE COROTATIONAL STRESS 
C        ---------------------------------------------------------------------------------------------------------------
C        sigma_co_ji = 2/J *(dWdI1+dWdI2*I1b)*B^star_ij - 
C 					   2/J * dWdI2*B^star_ik*B^star_kj + 
C 					  (dWdJ1 - 2*I1b/(3*J) * dWdI1 - (4*I2b)/(3*J) * dWdI2)*delta_ij
C  		 ---------------------------------------------------------------------------------------------------------------
		 p1 = (2.d0/detJ)*(dWdI1+dWdI2*I1b)
		 p2 = (2.d0/detJ)*dWdI2
		 p3 = (dWdJ1 - (2.d0*I1b)/(3.d0*detJ) * dWdI1 - (4.d0*I2b)/(3.d0*detJ) * dWdI2)
		 stressNew(k,1) = p1 * Bbar11 - p2*(Bbar11*Bbar11+Bbar12*Bbar12) + p3
		 stressNew(k,2) = p1 * Bbar22 - p2*(Bbar12*Bbar12+Bbar22*Bbar22) + p3
		 stressNew(k,3) = p1 * Bbar33 - p2*(Bbar33*Bbar33) + p3
		 stressNew(k,4) = p1 * Bbar12 - p2*(Bbar11*Bbar12+Bbar12*Bbar22)		 
C
C        UPDATE SPECIFIC INTERNAL ENERGY
C        ---------------------------------------------------------------------------------------------------------------
         u1 = half * ( (stressOld(k,1)+stressNew(k,1))*strainInc(k,1) +
     $                 (stressOld(k,2)+stressNew(k,2))*strainInc(k,2) +
     $                 (stressOld(k,3)+stressNew(k,3))*strainInc(k,3) +
     $                  two * ( (stressOld(k,4) + stressNew(k,4))*
     $                           strainInc(k,4) ) )
C
         enerInternNew(k) = enerInternOld(k) + u1 / density(k) 
C
      END DO
C
	  END IF
C
      RETURN
C
      END
  