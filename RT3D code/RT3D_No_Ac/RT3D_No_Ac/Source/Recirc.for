      SUBROUTINE RECIRC(IN,IOUT,KPER,NCOL,NROW,NLAY,NCOMP,ICBUND,CNEW
     &,CRCH,CEVT,MXSS,NSS,SS,SSMC,SSMC0,ICOMP,MCOMP,TIME2)
C ********************************************************************
C THIS SUBROUTINE READS CONCENTRATIONS OF SOURCES OR SINKS NEEDED BY
C THE SINK AND SOURCE MIXING (SSM) PACKAGE.
C ********************************************************************
C last modified: 06-23-98

      USE SSMDECAY_ARRAY
      IMPLICIT	NONE

	INTEGER, SAVE :: FIRST_TIME = 1

      INTEGER   IN,IOUT,KPER,NCOL,NROW,NLAY,NCOMP,ICBUND,
     &          MXSS,NSS,JJ,II,KK,NUM,IQ,INCRCH,INCEVT,NTMP,INDEX,
     &          ierr, I,J,K, ICOMP, YESNO, KPER_OLD, MCOMP
      REAL     CRCH,CEVT,SS,SSMC,CSS,CNEW,CAVG, SUMCK,TIME2,SSMC0 

	REAL, SAVE ::  COND, CONDSUM

  
      LOGICAL	FWEL,FDRN,FRIV,FGHB,FRCH,FEVT
      CHARACTER ANAME*24,TYPESS(-1:5)*15
      DIMENSION SS(6,MXSS),SSMC(NCOMP,MXSS),CRCH(NCOL,NROW,NCOMP),
     &          CEVT(NCOL,NROW,NCOMP), SSMC0(NCOMP,MXSS),
     &          ICBUND(NCOL,NROW,NLAY,NCOMP),CNEW(NCOL,NROW,NLAY,NCOMP),
     &		  COND(150,150,50),CAVG(NCOMP), SUMCK(NCOMP)

      COMMON   /FC/FWEL,FDRN,FRCH,FEVT,FRIV,FGHB


	INTEGER*4 RESULT, NUMKOUNT

	IF(FIRST_TIME.EQ.1) THEN
        OPEN(UNIT=999, FILE='cond.dat')
 	DO K=1,NLAY
	DO I=1,NROW
	DO J=1,NCOL
	READ(999,*) COND(J,I,K) 
	ENDDO
	ENDDO
	ENDDO
	CLOSE(999)
	KPER_OLD=0
	FIRST_TIME=0
	ELSE
	ENDIF

C-----Calculate CONDSUM------------------------
	YESNO=0

C	IF(KPER.NE.KPER_OLD) THEN
	CONDSUM = 0.0
	do NUM = 1,NSS
	K=SS(1,NUM)
	I=SS(2,NUM)
	J=SS(3,NUM)
	IQ = SS(6,NUM)
	
	IF(IQ.EQ.-3) THEN

	CONDSUM=cond(j,i,k) + CONDSUM
	YESNO=1
	ELSE
	ENDIF
	enddo
	KPER_OLD=KPER
C	ENDIF


C----------------RECIRCULATION LOGIC-----------------------------

	IF(YESNO.EQ.0) GOTO 1996

C	WRITE(*,*) 'ICOMP, MCOMP=', ICOMP, MCOMP


	SUMCK(ICOMP)=0.0

	DO 123 NUM=1,NSS 
	K=SS(1,NUM)
	I=SS(2,NUM)
	J=SS(3,NUM)
	IQ=SS(6,NUM)


	IF (IQ.EQ.-3) THEN  !for all extraction wells


C	Calculate the average concentrations:

	SUMCK(ICOMP) =(cond(j,i,k)*cnew(j,i,k,ICOMP)
     &			+ SUMCK(ICOMP))

	ENDIF
123	CONTINUE


	cavg(ICOMP)=SUMCK(ICOMP)/CONDSUM


	IF (CAVG(1).LT.(0.032)) CAVG(1) = 0.032


444	FORMAT(5(I4,1X),5(G12.6,1X))

C	INJECTION LOGIC:
	open(unit=6912, file='cavg.txt',access='APPEND')

	NUMKOUNT=0
	DO 124 NUM=1,NSS 
	K=SS(1,NUM)
	I=SS(2,NUM)
	J=SS(3,NUM)
	IQ=SS(6,NUM)
	
	
	IF (IQ.EQ. 2) THEN   
	   
	NUMKOUNT=NUMKOUNT+1
	ssmc(ICOMP,num)=ssmc0(ICOMP,num)+cavg(ICOMP)    !+cold(j,i,k,ICOMP)

  	CSS = SSMC(1,num)
	SS(4,NUM)=CSS

	GOTO 9200

c	IF((KPER.NE.1).OR.(KPER.NE.3)) THEN
c	Force Bromide:
c	IF(ICOMP.EQ.5) THEN
c	ssmc(5,num)=25.0
cC	if ((kper.eq.45).or.(kper.eq.47)) ssmc(1,num) = 180.0
cC	if ((kper.eq.1).or.(kper.eq.3)) ssmc(5,num) = 180.0
c	ENDIF
c	ENDIF
9200	CONTINUE

	IF((NUMKOUNT.EQ.1).AND.(ICOMP.EQ.MCOMP)) THEN
	WRITE(6912, 6914) KPER,TIME2,I,J,K,(SSMC(ICOMP,NUM),ICOMP=1,MCOMP)
	ENDIF

	ENDIF                 


6914	FORMAT(I4,2X,G12.6,1X,3(I4),2X,<NCOMP>(E12.6,2X))
6911	FORMAT(5(I4,1X),2(G12.6,1X))
255	FORMAT(2(I3,1X),E12.6,1X,3(I3,1X),7(F12.6,1X))
6913	FORMAT(2(I4,1X),G12.6,1X,3(I4,1X),1(E12.6,2X),<3*NCOMP>(E12.6,2X))

124	CONTINUE

	CLOSE(6912)


1996	CONTINUE
      RETURN
      END