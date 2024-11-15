      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, 
     &                STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, 
     &                CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, 
     &                DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, 
     &                KSPT, KSTEP, KINC)


      INCLUDE 'ABA_PARAM.INC'

      ! Declare variables
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS), 
     & STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), 
     & PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), 
     & DDSDDT(NTENS), DRPLDE(NTENS), SSE(1), SPD(1), SCD(1), 
     & DRPLDT(1), DPRED(1), DFGRD1(3,3), JSTEP(4), 
     & eelas(NTENS), eplas(NTENS), flow(NTENS), olds(NTENS), 
     & oldpl(NTENS) 

      PARAMETER(toler=1.d-6, newton=20)

      ! Material properties
      REAL*8 E, NU, SIGY, G, K, STR, DEQPL, PLASTIC_MULTIPLIER, YIELD_FUNC

      INTEGER U, J, P, O, I, KEWTON

      ! Input variables
      REAL*8 :: strain_rate, coeff0
      REAL*8 :: coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7
      REAL*8 :: coeff8, coeff9

      ! Coefficients 
      coeff0 = 1042.805131
      coeff1 = -2730.61563  ! Constant term
      coeff2 = 373117.74     ! Coefficient for strain^3 373117.74
      coeff3 = -3.06        ! Coefficient for strain^2 * strain_rate
      coeff4 = 3019.59631   ! Coefficient for strain * strain_rate^2
      coeff5 = 84164.4140   ! Coefficient for strain_rate^3 84164.4140
      coeff6 = -1.242442
      coeff7 = 0 !-38352281
      coeff8 = -11.3477767
      coeff9 = 0.0279732101
      
      E = PROPS(1)
      NU = PROPS(2)
      SIGY = PROPS(3)
      edot = PROPS(4)
      g = PROPS(5)

      ! Initialize DDSDDE to zero
      DDSDDE = 0.D0
      
      ! Call rotsig for elastic and plastic components
      CALL rotsig(STATEV(1), DROT, eelas, 2, NDI, NSHR)
      CALL rotsig(STATEV(NTENS + 1), DROT, eplas, 2, NDI, NSHR)
      
      DEQPL = STATEV(1 + 2 * NTENS)

      olds = STRESS
      oldpl = eplas

      ! Elastic modulus calculations
      EG = E / (1.D0 + NU) / 2.D0
      ELAM = (E / (1.D0 - 2.D0 * NU) - 2.D0 * EG) / 3.D0

      ! Fill DDSDDE matrix
      DO I = 1, 3
          DO J = 1, 3
              DDSDDE(J, I) = ELAM
          END DO
          DDSDDE(I, I) = 2.D0 * EG + ELAM
      END DO
      
      DO I = 4, NTENS
          DDSDDE(I, I) = EG
      END DO
      
      ! Update stress and elastic strain
      STRESS = STRESS + MATMUL(DDSDDE, DSTRAN)
      eelas = eelas + DSTRAN
      
      ! Calculate effective stress and check for yielding
      SMises = (STRESS(1) - STRESS(2)) ** 2 
      SMises = SMises + (STRESS(2) - STRESS(3)) ** 2
      SMises = SMises + (STRESS(3) - STRESS(1)) ** 2
      DO I = 4, NTENS
          SMises = SMises + 6.D0 * STRESS(I) ** 2
      END DO
      SMises = SQRT(SMises / 2.D0)

      ! Calculate flow stress
      !Sf = coeff1 * eqplas + coeff2 * strain_rate
      !Sf = Sf + coeff3 * eqplas ** 2
      !Sf = Sf + coeff4 * eqplas * strain_rate
      !Sf = Sf + coeff5 * strain_rate ** 2                                                                                           
      !Sf = Sf + coeff6 * eqplas ** 3
      !Sf = Sf + coeff7 * eqplas ** 2 * strain_rate
      !Sf = Sf + coeff8 * eqplas * strain_rate ** 2
      !Sf = Sf + coeff9 * strain_rate ** 3
      Sf = 115.87000 + 581.17094 * eqplas + 24694.12860 * edot - 1168663.06000 * g
      Sf = Sf - 1368.07412 * eqplas**2 + 12543.13140 * eqplas * edot
      Sf = Sf - 29095.13080 * eqplas * g - 2102507.11000 * edot**2
      Sf = Sf - 9820971.89000 * edot * g + 8222685370.00000 * g**2

      
      
	  !Sf=1042.805131-2.73061563e+03*eqplas+3.73117740e+05*edot
	  !Sf=Sf-3.06874783e+00*g+3.01959631e+03*eqplas*eqplas
	  !Sf=Sf+8.41644140e+04*eqplas*edot-1.24244211e+00*eqplas*g
	  !Sf=Sf-3.83522810e+07*edot*edot-1.13477767e+01*edot*g
	  !Sf=Sf+2.79732101e-02*g**2 
	  !Ezil

	  !eqplas - strain
	  !edot - strain rate
	  !g - grain size

      IF (SMises > (1.D0 + toler) * Sf) THEN
          Sh = (STRESS(1) + STRESS(2) + STRESS(3)) / 3.D0
          flow(1:3) = (STRESS(1:3) - Sh) / SMises
          flow(4:NTENS) = STRESS(4:NTENS) / SMises
          DEQPL = 0.D0
          eqplas = DEQPL + DEQPL

          ! Compute tangent modulus Et
          !Et = coeff1 + 2.0 * coeff3 * eqplas
          !Et = Et + coeff4 * strain_rate
          !Et = Et + 3.0 * coeff6 * eqplas ** 2 + 2.0 * coeff7 * eqplas * strain_rate
          !Et = Et + coeff8 * strain_rate ** 2
          !Et = (Et/ SIGY)
		  
		  
		  Et = 12543.1314 * edot - 2736.14824 * eqplas - 29095.1308 * g + 581.17094
		  !Et = (-2730.61563 ) + 3019.59631*eqplas**2 + 84164.4140*edot 
		  !Et= Et+ (-1.242442)*g    !ours


          ! Newton-Raphson iteration
          DO KEWTON = 1, NEWTON
              rhs = SMises - (3.D0 * EG) * DEQPL - Sf
              DEQPL = DEQPL + rhs / ((3.D0 * EG) + Et)
              
              !Sf = coeff1 * eqplas + coeff2 * strain_rate
              !Sf = Sf + coeff3 * eqplas ** 2 + coeff4 * eqplas * strain_rate
              !Sf = Sf + coeff5 * strain_rate ** 2
              !Sf = Sf + coeff6 * eqplas ** 3
              !Sf = Sf + coeff7 * eqplas ** 2 * strain_rate
              !Sf = Sf + coeff8 * eqplas * strain_rate ** 2
              !Sf = Sf + coeff9 * strain_rate ** 3 + coeff0
            
              !Et = coeff1 + 2.0 * coeff3 * eqplas
              !Et = Et + coeff4 * strain_rate
              !Et = Et + 3.0 * coeff6 * eqplas ** 2 + 2.0 * coeff7 * eqplas * strain_rate
              
              Stress (SCAU11) = 115.87000 + (581.17094) * E11 + (24694.12860) * strain_rate + (-1168663.06000) * grain_size + (-1368    .07412) * E11^2 + (12543.13140) * E11 * strain_rate + (-29095.13080) * E11 * grain_size + (-2102507.11000) * strain_rate^2 + (-9820971.89000) * strain_rate * grain_size + (8222685370.00000) * grain_size^2

              !Et = Et + coeff8 * strain_rate ** 2
			 Sf = 115.87000 + 581.17094 * (eqplas + 1.d-12 + deqpl)
             Sf = Sf + 24694.12860 * edot - 1168663.06000 * g
             Sf = Sf - 1368.07412 * (eqplas + 1.d-12 + deqpl)**2
             Sf = Sf + 12543.13140 * (eqplas + 1.d-12 + deqpl) * edot
             Sf = Sf - 29095.13080 * (eqplas + 1.d-12 + deqpl) * g
             Sf = Sf - 2102507.11000 * edot**2 - 9820971.89000 * edot * g
             Sf = Sf + 8222685370.00000 * g**2

			  
			  
			  
			  !Sf=1042.805131-2.73061563e+03*(eqplas + 1.d-12 + deqpl)
			  !Sf=Sf+3.73117740e+05*edot-3.06874783e+00*g
			  !Sf=Sf+3.01959631e+03*(eqplas + 1.d-12 + deqpl)
			  !Sf=Sf+8.41644140e+04*(eqplas + 1.d-12 + deqpl)*edot
			  !Sf=Sf-1.24244211e+00*(eqplas + 1.d-12 + deqpl)*g
			  !Sf=Sf-3.83522810e+07*edot**2-1.13477767e+01*edot*g
			  !Sf=Sf+2.79732101e-02*g**2 !Ezil
				!Sf=1042.805131-2.73061563e+03*eqplas+3.73117740e+05*edot-3.06874783e+00*g+3.01959631e+03*eqplas*eqplas+8.41644140e+04*eqplas*edot-1.24244211e+00*eqplas*g-3.83522810e+07*edot*edot-1.13477767e+01*edot*g+2.79732101e-02*g*g
			  !1, 4*2e, 5*edot, 6*g
			  
			  ! (eqplas + 1.d-12 + deqpl) - for changing the values in loop
			  
			  
			  Et = 12543.1314 * edot - 2736.14824 * (eqplas+1.d-12+deqpl)**2 - 29095.1308 * g + 581.17094 
			  !Et = (-2730.61563 ) + 3019.59631*(eqplas + 1.d-12 + deqpl)**2 + 84164.4140*edot + (-1.242442)g!ezil
          
              
              IF (ABS(rhs) < toler * SIGY) EXIT
          END DO
          
          IF (KEWTON .EQ. NEWTON) WRITE(7,*) !'WARNING: plasticity loop failed'
          
          ! Update stresses and strains
          STRESS(1:3) = flow(1:3) * Sf + Sh
          eplas(1:3) = eplas(1:3) + 3.D0 / 2.D0 * flow(1:3) * DEQPL
          eelas(1:3) = eelas(1:3) - 3.D0 / 2.D0 * flow(1:3) * DEQPL
          STRESS(4:NTENS) = flow(4:NTENS) * Sf
          eplas(4:NTENS) = eplas(4:NTENS) + 3.D0 * flow(4:NTENS) * DEQPL
          eelas(4:NTENS) = eelas(4:NTENS) - 3.D0 * flow(4:NTENS) * DEQPL
          eqplas = eqplas + DEQPL
          
          ! Calculate plastic strain energy density
          DO I = 1, NTENS
              SPD = SPD + (STRESS(I) + olds(I)) * (eplas(I) - oldpl(I)) / 2.D0
          END DO

          ! Formulate the Jacobian (material tangent)
          EFFG = EG * Sf / SMises
          EFFLAM = (E / (1.D0 - 2.D0 * NU) - 2.D0 * EFFG) / 3.D0
          EFFHRD = 3.D0 * EG * Et / (3.D0 * EG + Et) - 3.D0 * EFFG
          
          DO I = 1, 3
              DO J = 1, 3
                  DDSDDE(J, I) = EFFLAM
              END DO
              DDSDDE(I, I) = 2.D0 * EFFG + EFFLAM
          END DO

          DO I = 4, NTENS
              DDSDDE(I, I) = EFFG
          END DO
          
          DO I = 1, NTENS
              DO J = 1, NTENS
                DDSDDE(J, I) = DDSDDE(J, I) + EFFHRD * flow(J) * flow(I)
              END DO
          END DO
      END IF
      
      ! Update state variables
      STATEV(1:NTENS) = eelas
      STATEV((NTENS + 1):2 * NTENS) = eplas
      STATEV(1 + 2 * NTENS) = eqplas

      RETURN
      END SUBROUTINE UMAT