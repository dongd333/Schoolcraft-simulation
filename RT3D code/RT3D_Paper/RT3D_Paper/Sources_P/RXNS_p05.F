c
c Reaction package discussed in Example-2 (method-1)
c Please refer to the RT3D user manual for further details
c
	SUBROUTINE rxns(ncomp,nvrxndata,j,i,k,y,dydt,
     &         poros,rhob,reta,rc,nlay,nrow,ncol,vrc,kper)
C*Block 1:**************************************************************
c List of calling arguments
c ncomp - Total number of components
c nvrxndata - Total number of variable reaction parameters to be input via RCT file
c J, I, K - node location (used if reaction parameters are spatially variable)
c y - Concentration value of all component at the node [array variable y(ncomp)]
c dydt - Computed RHS of your differential equation [array variable dydt(ncomp)]
c poros - porosity of the node
c reta -  Retardation factor [ignore dummy reta values of immobile species]
c rhob -  bulk density of the node
c rc - Stores spatially constant reaction parameters (can dimension upto 100 values)
c nlay, nrow, ncol - Grid size (used only for dimensioning purposes)
c vrc - Array variable that stores spatially variable reaction parameters
C*End of Block 1********************************************************

C*Block 2:**************************************************************
c*    *Please do not modify this standard interface block*
c      !DEC$ ATTRIBUTES DLLEXPORT :: rxns
      IMPLICIT NONE
      INTEGER ncol,nrow,nlay
      INTEGER ncomp,nvrxndata,j,i,k,ii,kk,kper
      INTEGER, SAVE :: First_time=1
      DOUBLE PRECISION y,dydt,poros,rhob,reta
      DOUBLE PRECISION rc,vrc
      Double Precision KCT,GAMMA,KATT,KDET,MUMAX,BKC,KD(44),RCT(44)
      DIMENSION y(ncomp),dydt(ncomp),rc(100)
      DIMENSION vrc(ncol,nrow,nlay,nvrxndata),reta(ncomp)
C*End of block 2********************************************************

C*Block 3:************************************************************** 
c     *Declare your problem-specific new variables here*
c     INTEGER 
      DOUBLE PRECISION RA,RN,MU,YA,YN,KSA,KSN,BETA
      DOUBLE PRECISION RHO,FRAC,POR,YNB,DELTA
      DOUBLE PRECISION CT,KCM,KCIMM,ACE,NITR,TRACER,SCT,AA,BB
C*End of block 3********************************************************

c     Write(*,*) 'KPER = ', KPER
c     Pause   'Check'
c     Units:  Time (Days), Length (cm), Weight (mg)

C*Block 4:**************************************************************
c      *Initilize reaction parameters here, if required*
      RA = 1.0
      RN = 1.0
      MUMAX = 0.31         ! per day (regular = 3.1 per day)
      YA = 0.4            ! mg cells / mg Acetate
      YN = 0.25           ! mg cells / mg Nitrate
      YNB = 0.456
      KATT = 9.0          ! per day (regular = 0.9 per day)
      !KDET = 0.04
      KSA = 1.0           ! mg/L
      KSN = 12.0          ! mg/L
      
      BETA = 0.36
      FRAC = 0.437
      RHO = 1.63D+06
      POR = 0.35          ! This leads v=15 cm/day with no fluid removal

C     New parameters from compact scheme
      KCT = 0.05         ! 0.05 0.12
      GAMMA = 3.84       ! 3.84 18.89
      KDET = 0.03        ! 0.03 0.048
      BKC = 0.185         ! 0.185 0.136
      
      If ((KPER .EQ. 45) .or. (KPER .EQ. 47)) Then
          KATT = KATT * 10.0
      End If
      
      IF (First_time .EQ. 1) THEN
          Do KK = 1,44
              if (KK .LE. 15) then
                  KD(KK) = 1.45D-07
                  RCT(KK) = 1.0 + RHO*KD(KK)*FRAC/POR
              elseif ((KK .GT. 15) .AND. (KK .LE. 28)) then
                  KD(KK) = 1.65D-07
                  RCT(KK) = 1.0 + RHO*KD(KK)*FRAC/POR
              else
                  KD(KK) = 3.53D-07
                  RCT(KK) = 1.0 + RHO*KD(KK)*FRAC/POR
              end if
          End Do
          First_time = 0 !reset First_time to skip this block later
      END IF
      
C*End of block 4********************************************************

C*Block 5:**************************************************************
c      *Assign or compute values for new variables, if required*
       CT = y(1)
       KCM = y(2)
       ACE = y(3)
       NITR = y(4)
       TRACER = y(5)
       SCT = y(6)
       KCIMM = y(7)
       
       AA = ACE/(KSA+ACE)
       BB = NITR/(KSN+NITR)
       MU = MUMAX*AA*BB
C*End of block 5********************************************************

C*Block 6:**************************************************************
c      *Differential Reaction Equations*
       dydt(1) = -KCT*CT*KCM/RCT(K) - 
     & RHO*BETA/(POR*RCT(K))*((1.0-FRAC)*KD(K)*CT - SCT)
     & -(RHO*FRAC*KD(K)/POR)*KCT*CT*KCIMM  
       dydt(2) = (MU-BKC*(1.0-AA)-KATT)*KCM+KDET*KCIMM*(1.0-AA)
       dydt(3) = -MU/(YA*RA)*(KCM+KCIMM)              ! Acetate
       dydt(4) = -(MU/YN)*(KCM+KCIMM)-BKC*(KCM+KCIMM)*(1.0-AA)/YNB
     1    - GAMMA*(KCM+KCIMM)*BB                      ! Nitrate
       dydt(5) = 0.0                                  ! Tracer
       dydt(6) = BETA*((1.0-FRAC)*KD(K)*CT - SCT)-KCT*KCIMM*SCT   ! SCT
	 dydt(7) = (MU-BKC*(1.0-AA)-KDET*(1.0-AA))*KCIMM + KATT*KCM ! Immob KC

C*End of block 6********************************************************

      RETURN
      END

