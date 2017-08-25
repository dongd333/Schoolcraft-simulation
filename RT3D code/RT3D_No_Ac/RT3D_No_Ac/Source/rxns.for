C***********************Reaction Model****************************C


C*********************************************************************C	

 	SUBROUTINE Rxns(ncomp,nvrxndata,j,i,k,y,dydt,
     &             poros,rhob,reta,rc,nlay,nrow,ncol,vrc,kper)
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
      !MS$ATTRIBUTES DLLEXPORT :: rxns
      IMPLICIT NONE
      INTEGER ncol,nrow,nlay
      INTEGER ncomp,nvrxndata,j,i,k, ii,KK, kper

      INTEGER, SAVE :: First_time=1
      DOUBLE PRECISION y,dydt,poros,rhob,reta
      DOUBLE PRECISION rc,vrc
	DOUBLE PRECISION, SAVE :: KCT,GAMMA,KATT,KDET,MUMAX,BKC, KD(44),
     & RCT(44)

      DIMENSION y(ncomp),dydt(ncomp),rc(50)
      DIMENSION vrc(ncol,nrow,nlay,nvrxndata),reta(ncomp)
C*End of Block 2********************************************************

C*Block 3:************************************************************** 
c      *Declare your problem-specific new variables here*
c	INTEGER 
       DOUBLE PRECISION RA,RN,MU,YA,
	1	YN,KSA,KSN,BETA,RHO,FRAC,POR, YNB, DELTA
       DOUBLE PRECISION  CT,KCM,KCIMM,ACE,NITR,TRACER,SCT,AA,BB


c	WRITE(*,*) 'KPER = ', KPER
c	PAUSE 'CHECK'
c
c	Unites are in time (days), length (cm), weight (milligram)
c	Area of Cross-section = 21.2371 cm2


	  RA=1.0
	  RN=1.0
	  MUMAX=3.1   ! per day (regular value = 3.1 per day)
c	  BKC=0.1
	  YA=0.40     !mg cells/mg acetate
	  YN=0.25	  !mg cells/mg nitrate
	  YNB=0.456
	  KATT=0.9    ! (regular value of 0.9  per day)
c	  KDET=0.18   !per day
	  KSA=1.0	  != 1.0 mg/L	
	  KSN=12.0	  != 12.0 mg/L

	  BETA=0.36
C	  KD=3.9D-07 
	  FRAC=0.437
	  RHO = 1.63D+06
	  POR = 0.35  !This gives a velocity of 15cm/day with no fluid removal

c
	 CT     = y(1)
	 KCM     = y(2)
	 ACE   = y(3)
	 NITR    = y(4)
	 TRACER   = y(5)
	 SCT    =   y(6)
	 KCIMM     = y(7)
	 
	AA=ACE/(KSA+ACE)
	BB=NITR/(KSN+NITR)
	MU=MUMAX*AA*BB*2

C	The following are the new parameters from compact scheme:
	 KCT =   0.05 
	 GAMMA=  18.89  
	 !KDET=   0.04  
	 !BKC=    0.136
       
       KDET=   0.04  
	 BKC=    0.068

	IF((KPER.EQ.45).OR.(KPER.EQ.47)) THEN
	KATT = KATT*10.0
	ENDIF



	IF(First_time.EQ.1) THEN	 
	DO KK=1,44
	IF (KK.LE.15) THEN
	KD(KK)=1.45D-07
	RCT(KK) = 1.0 + (RHO*KD(KK)*FRAC/POR)
	ELSEIF ((KK.GT.15).AND.(KK.LE.28)) THEN
	KD(KK)=1.65D-07
	RCT(KK) = 1.0 + (RHO*KD(KK)*FRAC/POR)
	ELSEIF ((KK.GT.28).AND.(KK.LE.44)) THEN
	KD(KK)=3.53D-07
	RCT(KK) = 1.0 + (RHO*KD(KK)*FRAC/POR)
	ENDIF
	ENDDO
	First_time = 0    !First_time + 1
	ENDIF

C	Monod Nitrate for degradation terms:


C	WRITE(*,*) K, KD(K), RCT(K)

       dydt(1) = -KCT*CT*KCM/RCT(K) - 
     & (RHO*BETA/(POR*RCT(K)))*((1.0-FRAC)*KD(K)*CT - SCT)
     & -(RHO*FRAC*KD(K)/POR)*KCT*CT*KCIMM

C      dydt(1) = -KCT*CT*(KCM+KCIMM)/RCT(K) - (RHO*BETA/(POR*RCT(K)))
C	1 *((1.0-FRAC)*KD(K)*CT - SCT)

 	 !dydt(2) = (MU-BKC*(1.0-AA)-KATT)*KCM+KDET*KCIMM*(1.0-AA)
       dydt(2) = (MU-BKC-KATT)*KCM+KDET*KCIMM

       dydt(3) = 	 -(MU/(YA*RA))*(KCM+KCIMM)

	 !dydt(4) = -(MU/YN)*(KCM+KCIMM) -(BKC*(1.0-AA)/YNB)*(KCM+KCIMM) 
	!1   - gamma*(KCM+KCIMM)*BB
       
       dydt(4)=-(MU/YN+gamma*BB)*(KCM+KCIMM)-BKC/YNB*(KCM+KCIMM)

	 dydt(5)=0.0

	 dydt(6) = BETA*((1.0-FRAC)*KD(K)*CT - SCT)-KCT*KCIMM*SCT

C	 dydt(6) = BETA*((1.0-FRAC)*KD(K)*CT - SCT)

	 !dydt(7) =  (MU-BKC*(1.0-AA)-KDET*(1.0-AA))*KCIMM + KATT*KCM
       dydt(7) =  (MU-BKC-KDET)*KCIMM + KATT*KCM



C*End of Block 6********************************************************

       RETURN
       END
