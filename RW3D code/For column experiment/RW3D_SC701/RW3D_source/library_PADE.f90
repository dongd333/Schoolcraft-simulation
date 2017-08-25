
  module library_PADE_expokit

  public

  contains

!----------------------------------------------------------------------|                                                                                                                                
!     EXPOKIT PADE APPROX. FOR A GENERAL MATRIX
!     AND ASSOCIATED SUBROUTINES AND FUNCTIONS FROM LAPACK AND BLAS
!     LIBRARIES
!----------------------------------------------------------------------|
  function expokit_Pade (A,m,t) result (p)
  implicit none
	integer, intent(in) :: m
	real*8,  intent(in) :: t
	real*8,  intent(in) :: A(m,m)
	real*8              :: P(m,m)
    integer, parameter  :: ideg=6
	integer             :: i,j,k,lda,ldh,ns,iflag,iexp,lwsp,liwsp
    integer             :: iwsp(m)
    real*8              :: wsp(4*m*m+ideg+1)
	real*8              :: Pcol
	lda   = m
	ldh   = m
	lwsp  = 4*ldh*ldh+ideg+1
	liwsp = ldh 
!
!....  Pade ...........................................................
      call DGPADM( ideg, m, t, A,lda, wsp,lwsp, iwsp,iexp,ns, iflag )
!......................................................................
      do j=1,m
	    do i=1,m
		P(i,j)=wsp(iexp+(j-1)*m+i-1)
		end do
	  end do
!......check we obtain a stochastic matrix: sum of a column equal to 1
      do j=1,m
	     Pcol = sum(P(1:m,j:j))
         P(1:m,j:j)=P(1:m,j:j)/pcol	  
	  end do      
  end function expokit_Pade
!----------------------------------------------------------------------|
      subroutine DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )                                                                                                                                  
                                                                                                                                                                                                        
      implicit none                                                                                                                                                                                     
      integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)                                                                                                                                             
      double precision t, H(ldh,m), wsp(lwsp)                                                                                                                                                           
                                                                                                                                                                                                        
!-----Purpose----------------------------------------------------------|                                                                                                                                
!                                                                                                                                                                                                       
!     Computes exp(t*H), the matrix exponential of a general matrix in                                                                                                                                  
!     full, using the irreducible rational Pade approximation to the                                                                                                                                    
!     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),                                                                                                                                  
!     combined with scaling-and-squaring.                                                                                                                                                               
!                                                                                                                                                                                                       
!-----Arguments--------------------------------------------------------|                                                                                                                                
!                                                                                                                                                                                                       
!     ideg      : (input) the degre of the diagonal Pade to be used.                                                                                                                                    
!                 a value of 6 is generally satisfactory.                                                                                                                                               
!                                                                                                                                                                                                       
!     m         : (input) order of H.                                                                                                                                                                   
!                                                                                                                                                                                                       
!     H(ldh,m)  : (input) argument matrix.                                                                                                                                                              
!                                                                                                                                                                                                       
!     t         : (input) time-scale (can be < 0).                                                                                                                                                      
!                                                                                                                                                                                                       
!     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.                                                                                                                                            
!                                                                                                                                                                                                       
!     ipiv(m)   : (workspace)                                                                                                                                                                           
!                                                                                                                                                                                                       
!>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)                                                                                                                                
!                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)                                                                                                                                
!                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                                                                                                                
!                 NOTE: if the routine was called with wsp(iptr),                                                                                                                                       
!                       then exp(tH) will start at wsp(iptr+iexph-1).                                                                                                                                   
!                                                                                                                                                                                                       
!     ns        : (output) number of scaling-squaring used.                                                                                                                                             
!                                                                                                                                                                                                       
!     iflag     : (output) exit flag.                                                                                                                                                                   
!                      0 - no problem                                                                                                                                                                   
!                     <0 - problem                                                                                                                                                                      
!                                                                                                                                                                                                       
!----------------------------------------------------------------------|                                                                                                                                
!     Roger B. Sidje (rbs@maths.uq.edu.au)                                                                                                                                                              
!     EXPOKIT: Software Package for Computing Matrix Exponentials.                                                                                                                                      
!     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998                                                                                                                                  
!----------------------------------------------------------------------|                                                                                                                                
!                                                                                                                                                                                                       
      integer mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget                                                                                                                                       
      double precision hnorm,scale,scale2,cp,cq                                                                                                                                                         
                                                                                                                                                                                                        
      intrinsic INT,ABS,DBLE,LOG,MAX                                                                                                                                                                    
                                                                                                                                                                                                        
!---  check restrictions on input parameters ...                                                                                                                                                        
      mm = m*m                                                                                                                                                                                          
      iflag = 0                                                                                                                                                                                         
      if ( ldh.lt.m ) iflag = -1                                                                                                                                                                        
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2                                                                                                                                                             
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DGPADM)'                                                                                                                                           
!                                                                                                                                                                                                       
!---  initialise pointers ...                                                                                                                                                                           
!                                                                                                                                                                                                       
      icoef = 1                                                                                                                                                                                         
      ih2 = icoef + (ideg+1)                                                                                                                                                                            
      ip  = ih2 + mm                                                                                                                                                                                    
      iq  = ip + mm                                                                                                                                                                                     
      ifree = iq + mm                                                                                                                                                                                   
!                                                                                                                                                                                                       
!---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;                                                                                                                                                    
!     and set scale = t/2^ns ...                                                                                                                                                                        
!                                                                                                                                                                                                       
      do i = 1,m                                                                                                                                                                                        
         wsp(i) = 0.0d0                                                                                                                                                                                 
      enddo                                                                                                                                                                                             
      do j = 1,m                                                                                                                                                                                        
         do i = 1,m                                                                                                                                                                                     
            wsp(i) = wsp(i) + ABS( H(i,j) )                                                                                                                                                             
         enddo                                                                                                                                                                                          
      enddo                                                                                                                                                                                             
      hnorm = 0.0d0                                                                                                                                                                                     
      do i = 1,m                                                                                                                                                                                        
         hnorm = MAX( hnorm,wsp(i) )                                                                                                                                                                    
      enddo                                                                                                                                                                                             
      hnorm = ABS( t*hnorm )                                                                                                                                                                            
      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of DGPADM.'                                                                                                                                   
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )                                                                                                                                                        
      scale = t / DBLE(2**ns)                                                                                                                                                                           
      scale2 = scale*scale                                                                                                                                                                              
!                                                                                                                                                                                                       
!---  compute Pade coefficients ...                                                                                                                                                                     
!                                                                                                                                                                                                       
      i = ideg+1                                                                                                                                                                                        
      j = 2*ideg+1                                                                                                                                                                                      
      wsp(icoef) = 1.0d0                                                                                                                                                                                
      do k = 1,ideg                                                                                                                                                                                     
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )                                                                                                                                    
      enddo                                                                                                                                                                                             
!                                                                                                                                                                                                       
!---  H2 = scale2*H*H ...                                                                                                                                                                               
!                                                                                                                                                                                                       
      call DGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )                                                                                                                                   
!                                                                                                                                                                                                       
!---  initialize p (numerator) and q (denominator) ...                                                                                                                                                  
!                                                                                                                                                                                                       
      cp = wsp(icoef+ideg-1)                                                                                                                                                                            
      cq = wsp(icoef+ideg)                                                                                                                                                                              
      do j = 1,m                                                                                                                                                                                        
         do i = 1,m                                                                                                                                                                                     
            wsp(ip + (j-1)*m + i-1) = 0.0d0                                                                                                                                                             
            wsp(iq + (j-1)*m + i-1) = 0.0d0                                                                                                                                                             
         enddo                                                                                                                                                                                          
         wsp(ip + (j-1)*(m+1)) = cp                                                                                                                                                                     
         wsp(iq + (j-1)*(m+1)) = cq                                                                                                                                                                     
      enddo                                                                                                                                                                                             
!                                                                                                                                                                                                       
!---  Apply Horner rule ...                                                                                                                                                                             
!                                                                                                                                                                                                       
      iodd = 1                                                                                                                                                                                          
      k = ideg - 1                                                                                                                                                                                      
 100  continue                                                                                                                                                                                          
      iused = iodd*iq + (1-iodd)*ip                                                                                                                                                                     
      call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iused),m,   &                                                                                                                                                
                   wsp(ih2),m, 0.0d0,wsp(ifree),m )                                                                                                                                                     
      do j = 1,m                                                                                                                                                                                        
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)                                                                                                                                 
      enddo                                                                                                                                                                                             
      ip = (1-iodd)*ifree + iodd*ip                                                                                                                                                                     
      iq = iodd*ifree + (1-iodd)*iq                                                                                                                                                                     
      ifree = iused                                                                                                                                                                                     
      iodd = 1-iodd                                                                                                                                                                                     
      k = k-1                                                                                                                                                                                           
      if ( k.gt.0 )  goto 100                                                                                                                                                                           
!                                                                                                                                                                                                       
!---  Obtain (+/-)(I + 2*(p\q)) ...                                                                                                                                                                     
!                                                                                                                                                                                                       
      if ( iodd .eq. 1 ) then                                                                                                                                                                           
         call DGEMM( 'n','n',m,m,m, scale,wsp(iq),m,   &                                                                                                                                                
                      H,ldh, 0.0d0,wsp(ifree),m )                                                                                                                                                       
         iq = ifree                                                                                                                                                                                     
      else                                                                                                                                                                                              
         call DGEMM( 'n','n',m,m,m, scale,wsp(ip),m,   &                                                                                                                                                
                      H,ldh, 0.0d0,wsp(ifree),m )                                                                                                                                                       
         ip = ifree                                                                                                                                                                                     
      endif                                                                                                                                                                                             
      call DAXPY( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )                                                                                                                                                     
      call DGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )                                                                                                                                              
      if ( iflag.ne.0 ) stop 'Problem in DGESV (within DGPADM)'                                                                                                                                         
      call DSCAL( mm, 2.0d0, wsp(ip), 1 )                                                                                                                                                               
      do j = 1,m                                                                                                                                                                                        
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0                                                                                                                                              
      enddo                                                                                                                                                                                             
      iput = ip                                                                                                                                                                                         
      if ( ns.eq.0 .and. iodd.eq.1 ) then                                                                                                                                                               
         call DSCAL( mm, -1.0d0, wsp(ip), 1 )                                                                                                                                                           
         goto 200                                                                                                                                                                                       
      endif                                                                                                                                                                                             
!                                                                                                                                                                                                       
!--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...                                                                                                                                                       
!                                                                                                                                                                                                       
      iodd = 1                                                                                                                                                                                          
      do k = 1,ns                                                                                                                                                                                       
         iget = iodd*ip + (1-iodd)*iq                                                                                                                                                                   
         iput = (1-iodd)*ip + iodd*iq                                                                                                                                                                   
         call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m,   &                                                                                                                                 
                      0.0d0,wsp(iput),m )                                                                                                                                                               
         iodd = 1-iodd                                                                                                                                                                                  
      enddo                                                                                                                                                                                             
 200  continue                                                                                                                                                                                          
      iexph = iput                                                                                                                                                                                      
      END SUBROUTINE                                                                                                                                                                                    
!----------------------------------------------------------------------|  




!----------------------------------------------------------------------|                                                                                                                                                                                               
!     Subset of BLAS and LAPACK routines used by EXPOKIT in the 
!                    General PADE approx                                                                                                                                                                                                                                                                                                                                                
!----------------------------------------------------------------------|                                                                                                                                

!----------------------------------------------------------------------|                                                                                                                                
      double precision function dcabs1(z)                                                                                                                                                                
      complex*16 z,zz                                                                                                                                                                                   
      double precision t(2)                                                                                                                                                                             
      equivalence (zz,t(1))                                                                                                                                                                             
      zz = z                                                                                                                                                                                            
      dcabs1 = dabs(t(1)) + dabs(t(2))                                                                                                                                                                                                                                                                                                                                                              
      end function                                                                                                                                                                                      
!----------------------------------------------------------------------|                                                                                                                                
      LOGICAL  FUNCTION LSAME( CA, CB )                                                                                                                                                         
!                                                                                                                                                                                                       
!  -- LAPACK auxiliary routine (version 1.1) --                                                                                                                                                         
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,                                                                                                                                       
!     Courant Institute, Argonne National Lab, and Rice University                                                                                                                                      
!     February 29, 1992                                                                                                                                                                                 
!                                                                                                                                                                                                       
!     .. Scalar Arguments ..                                                                                                                                                                            
      CHARACTER          CA, CB                                                                                                                                                                         
!     ..                                                                                                                                                                                                
!                                                                                                                                                                                                       
!  Purpose                                                                                                                                                                                              
!  =======                                                                                                                                                                                              
!                                                                                                                                                                                                       
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of                                                                                                                                    
!  case.                                                                                                                                                                                                
!                                                                                                                                                                                                       
!  Arguments                                                                                                                                                                                            
!  =========                                                                                                                                                                                            
!                                                                                                                                                                                                       
!  CA      (input) CHARACTER*1                                                                                                                                                                          
!  CB      (input) CHARACTER*1                                                                                                                                                                          
!          CA and CB specify the single characters to be compared.                                                                                                                                      
!                                                                                                                                                                                                       
!     .. Intrinsic Functions ..                                                                                                                                                                         
      INTRINSIC          ICHAR                                                                                                                                                                          
!     ..                                                                                                                                                                                                
!     .. Local Scalars ..                                                                                                                                                                               
      INTEGER            INTA, INTB, ZCODE                                                                                                                                                              
!     ..                                                                                                                                                                                                
!     .. Executable Statements ..                                                                                                                                                                       
!                                                                                                                                                                                                       
!     Test if the characters are equal                                                                                                                                                                  
!                                                                                                                                                                                                       
      LSAME = CA.EQ.CB                                                                                                                                                                                  
      IF( LSAME )   &                                                                                                                                                                                   
         RETURN                                                                                                                                                                                         
!                                                                                                                                                                                                       
!     Now test for equivalence if both characters are alphabetic.                                                                                                                                       
!                                                                                                                                                                                                       
      ZCODE = ICHAR( 'Z' )                                                                                                                                                                              
!                                                                                                                                                                                                       
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime                                                                                                                                    
!     machines, on which ICHAR returns a value with bit 8 set.                                                                                                                                          
!     ICHAR('A') on Prime machines returns 193 which is the same as                                                                                                                                     
!     ICHAR('A') on an EBCDIC machine.                                                                                                                                                                  
!                                                                                                                                                                                                       
      INTA = ICHAR( CA )                                                                                                                                                                                
      INTB = ICHAR( CB )                                                                                                                                                                                
!                                                                                                                                                                                                       
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN                                                                                                                                                          
!                                                                                                                                                                                                       
!        ASCII is assumed - ZCODE is the ASCII code of either lower or                                                                                                                                  
!        upper case 'Z'.                                                                                                                                                                                
!                                                                                                                                                                                                       
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32                                                                                                                                            
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32                                                                                                                                            
!                                                                                                                                                                                                       
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN                                                                                                                                                    
!                                                                                                                                                                                                       
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or                                                                                                                                
!        upper case 'Z'.                                                                                                                                                                                
!                                                                                                                                                                                                       
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.   &                                                                                                                                                     
             INTA.GE.145 .AND. INTA.LE.153 .OR.   &                                                                                                                                                     
             INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64                                                                                                                                           
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.   &                                                                                                                                                     
             INTB.GE.145 .AND. INTB.LE.153 .OR.   &                                                                                                                                                     
             INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64                                                                                                                                           
!                                                                                                                                                                                                       
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN                                                                                                                                                    
!                                                                                                                                                                                                       
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code                                                                                                                                  
!        plus 128 of either lower or upper case 'Z'.                                                                                                                                                    
!                                                                                                                                                                                                       
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32                                                                                                                                           
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32                                                                                                                                           
      END IF                                                                                                                                                                                            
      LSAME = INTA.EQ.INTB                                                                                                                                                                              
!                                                                                                                                                                                                       
!     RETURN                                                                                                                                                                                            
!                                                                                                                                                                                                       
!     End of LSAME                                                                                                                                                                                      
!                                                                                                                                                                                                       
      END FUNCTION        
!----------------------------------------------------------------------|                                                                                                                                
      SUBROUTINE XERBLA ( SRNAME, INFO )                                                                                                                                                                
!     ..    Scalar Arguments ..                                                                                                                                                                         
      INTEGER            INFO                                                                                                                                                                           
      CHARACTER*6        SRNAME                                                                                                                                                                         
!     ..                                                                                                                                                                                                
!                                                                                                                                                                                                       
!  Purpose                                                                                                                                                                                              
!  =======                                                                                                                                                                                              
!                                                                                                                                                                                                       
!  XERBLA  is an error handler for the Level 2 BLAS routines.                                                                                                                                           
!                                                                                                                                                                                                       
!  It is called by the Level 2 BLAS routines if an input parameter is                                                                                                                                   
!  invalid.                                                                                                                                                                                             
!                                                                                                                                                                                                       
!  Installers should consider modifying the STOP statement in order to                                                                                                                                  
!  call system-specific exception-handling facilities.                                                                                                                                                  
!                                                                                                                                                                                                       
!  Parameters                                                                                                                                                                                           
!  ==========                                                                                                                                                                                           
!                                                                                                                                                                                                       
!  SRNAME - CHARACTER*6.                                                                                                                                                                                
!           On entry, SRNAME specifies the name of the routine which                                                                                                                                    
!           called XERBLA.                                                                                                                                                                              
!                                                                                                                                                                                                       
!  INFO   - INTEGER.                                                                                                                                                                                    
!           On entry, INFO specifies the position of the invalid                                                                                                                                        
!           parameter in the parameter-list of the calling routine.                                                                                                                                     
!                                                                                                                                                                                                       
!                                                                                                                                                                                                       
!  Auxiliary routine for Level 2 Blas.                                                                                                                                                                  
!                                                                                                                                                                                                       
!  Written on 20-July-1986.                                                                                                                                                                             
!                                                                                                                                                                                                       
!     .. Executable Statements ..                                                                                                                                                                       
!                                                                                                                                                                                                       
      WRITE (*,99999) SRNAME, INFO                                                                                                                                                                      
!                                                                                                                                                                                                       
      STOP                                                                                                                                                                                              
!                                                                                                                                                                                                       
99999 FORMAT ( ' ** On entry to ', A6, ' parameter number ', I2,   &                                                                                                                                    
               ' had an illegal value' )                                                                                                                                                                
!                                                                                                                                                                                                       
!     End of XERBLA.                                                                                                                                                                                    
!                                                                                                                                                                                                       
      END SUBROUTINE  xerbla 

!----------------------------------------------------------------------|                                                                                                                                
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,   &                                                                                                                            
                         BETA, C, LDC )                                                                                                                                                                 
!     .. Scalar Arguments ..                                                                                                                                                                            
      CHARACTER*1        TRANSA, TRANSB                                                                                                                                                                 
      INTEGER            M, N, K, LDA, LDB, LDC                                                                                                                                                         
      DOUBLE PRECISION   ALPHA, BETA                                                                                                                                                                    
!     .. Array Arguments ..                                                                                                                                                                             
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )                                                                                                                                          
!     ..                                                                                                                                                                                                
!                                                                                                                                                                                                       
!  Purpose                                                                                                                                                                                              
!  =======                                                                                                                                                                                              
!                                                                                                                                                                                                       
!  DGEMM  performs one of the matrix-matrix operations                                                                                                                                                  
!                                                                                                                                                                                                       
!     C := alpha*op( A )*op( B ) + beta*C,                                                                                                                                                              
!                                                                                                                                                                                                       
!  where  op( X ) is one of                                                                                                                                                                             
!                                                                                                                                                                                                       
!     op( X ) = X   or   op( X ) = X',                                                                                                                                                                  
!                                                                                                                                                                                                       
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )                                                                                                                                
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.                                                                                                                                
!                                                                                                                                                                                                       
!  Parameters                                                                                                                                                                                           
!  ==========                                                                                                                                                                                           
!                                                                                                                                                                                                       
!  TRANSA - CHARACTER*1.                                                                                                                                                                                
!           On entry, TRANSA specifies the form of op( A ) to be used in                                                                                                                                
!           the matrix multiplication as follows:                                                                                                                                                       
!                                                                                                                                                                                                       
!              TRANSA = 'N' or 'n',  op( A ) = A.                                                                                                                                                       
!                                                                                                                                                                                                       
!              TRANSA = 'T' or 't',  op( A ) = A'.                                                                                                                                                      
!                                                                                                                                                                                                       
!              TRANSA = 'C' or 'c',  op( A ) = A'.                                                                                                                                                      
!                                                                                                                                                                                                       
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  TRANSB - CHARACTER*1.                                                                                                                                                                                
!           On entry, TRANSB specifies the form of op( B ) to be used in                                                                                                                                
!           the matrix multiplication as follows:                                                                                                                                                       
!                                                                                                                                                                                                       
!              TRANSB = 'N' or 'n',  op( B ) = B.                                                                                                                                                       
!                                                                                                                                                                                                       
!              TRANSB = 'T' or 't',  op( B ) = B'.                                                                                                                                                      
!                                                                                                                                                                                                       
!              TRANSB = 'C' or 'c',  op( B ) = B'.                                                                                                                                                      
!                                                                                                                                                                                                       
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  M      - INTEGER.                                                                                                                                                                                    
!           On entry,  M  specifies  the number  of rows  of the  matrix                                                                                                                                
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.                                                                                                                                
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  N      - INTEGER.                                                                                                                                                                                    
!           On entry,  N  specifies the number  of columns of the matrix                                                                                                                                
!           op( B ) and the number of columns of the matrix C. N must be                                                                                                                                
!           at least zero.                                                                                                                                                                              
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  K      - INTEGER.                                                                                                                                                                                    
!           On entry,  K  specifies  the number of columns of the matrix                                                                                                                                
!           op( A ) and the number of rows of the matrix op( B ). K must                                                                                                                                
!           be at least  zero.                                                                                                                                                                          
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  ALPHA  - DOUBLE PRECISION.                                                                                                                                                                           
!           On entry, ALPHA specifies the scalar alpha.                                                                                                                                                 
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is                                                                                                                                
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.                                                                                                                                        
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k                                                                                                                                
!           part of the array  A  must contain the matrix  A,  otherwise                                                                                                                                
!           the leading  k by m  part of the array  A  must contain  the                                                                                                                                
!           matrix A.                                                                                                                                                                                   
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  LDA    - INTEGER.                                                                                                                                                                                    
!           On entry, LDA specifies the first dimension of A as declared                                                                                                                                
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then                                                                                                                                
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at                                                                                                                                
!           least  max( 1, k ).                                                                                                                                                                         
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is                                                                                                                                
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.                                                                                                                                        
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n                                                                                                                                
!           part of the array  B  must contain the matrix  B,  otherwise                                                                                                                                
!           the leading  n by k  part of the array  B  must contain  the                                                                                                                                
!           matrix B.                                                                                                                                                                                   
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  LDB    - INTEGER.                                                                                                                                                                                    
!           On entry, LDB specifies the first dimension of B as declared                                                                                                                                
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then                                                                                                                                
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at                                                                                                                                
!           least  max( 1, n ).                                                                                                                                                                         
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  BETA   - DOUBLE PRECISION.                                                                                                                                                                           
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is                                                                                                                                
!           supplied as zero then C need not be set on input.                                                                                                                                           
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).                                                                                                                                             
!           Before entry, the leading  m by n  part of the array  C must                                                                                                                                
!           contain the matrix  C,  except when  beta  is zero, in which                                                                                                                                
!           case C need not be set on entry.                                                                                                                                                            
!           On exit, the array  C  is overwritten by the  m by n  matrix                                                                                                                                
!           ( alpha*op( A )*op( B ) + beta*C ).                                                                                                                                                         
!                                                                                                                                                                                                       
!  LDC    - INTEGER.                                                                                                                                                                                    
!           On entry, LDC specifies the first dimension of C as declared                                                                                                                                
!           in  the  calling  (sub)  program.   LDC  must  be  at  least                                                                                                                                
!           max( 1, m ).                                                                                                                                                                                
!           Unchanged on exit.                                                                                                                                                                          
!                                                                                                                                                                                                       
!                                                                                                                                                                                                       
!  Level 3 Blas routine.                                                                                                                                                                                
!                                                                                                                                                                                                       
!  -- Written on 8-February-1989.                                                                                                                                                                       
!     Jack Dongarra, Argonne National Laboratory.                                                                                                                                                       
!     Iain Duff, AERE Harwell.                                                                                                                                                                          
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.                                                                                                                                                   
!     Sven Hammarling, Numerical Algorithms Group Ltd.                                                                                                                                                  
!                                                                                                                                                                                                       
!                                                                                                                                                                                                       
!     .. External Functions ..                                                                                                                                                                          
!      LOGICAL            LSAME                                                                                                                                                                          
!      EXTERNAL           LSAME                                                                                                                                                                          
!     .. External Subroutines ..                                                                                                                                                                        
!      EXTERNAL           XERBLA                                                                                                                                                                         
!     .. Intrinsic Functions ..                                                                                                                                                                         
      INTRINSIC          MAX                                                                                                                                                                            
!     .. Local Scalars ..                                                                                                                                                                               
      LOGICAL            NOTA, NOTB                                                                                                                                                                     
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB                                                                                                                                             
      DOUBLE PRECISION   TEMP                                                                                                                                                                           
!     .. Parameters ..                                                                                                                                                                                  
      DOUBLE PRECISION   ONE         , ZERO                                                                                                                                                             
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )                                                                                                                                                  
!     ..                                                                                                                                                                                                
!     .. Executable Statements ..                                                                                                                                                                       
!                                                                                                                                                                                                       
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not                                                                                                                                
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows                                                                                                                                
!     and  columns of  A  and the  number of  rows  of  B  respectively.                                                                                                                                
!                                                                                                                                                                                                       
      NOTA  = LSAME( TRANSA, 'N' )                                                                                                                                                                      
      NOTB  = LSAME( TRANSB, 'N' )                                                                                                                                                                      
      IF( NOTA )THEN                                                                                                                                                                                    
         NROWA = M                                                                                                                                                                                      
         NCOLA = K                                                                                                                                                                                      
      ELSE                                                                                                                                                                                              
         NROWA = K                                                                                                                                                                                      
         NCOLA = M                                                                                                                                                                                      
      END IF                                                                                                                                                                                            
      IF( NOTB )THEN                                                                                                                                                                                    
         NROWB = K                                                                                                                                                                                      
      ELSE                                                                                                                                                                                              
         NROWB = N                                                                                                                                                                                      
      END IF                                                                                                                                                                                            
!                                                                                                                                                                                                       
!     Test the input parameters.                                                                                                                                                                        
!                                                                                                                                                                                                       
      INFO = 0                                                                                                                                                                                          
      IF(      ( .NOT.NOTA                 ).AND.   &                                                                                                                                                   
               ( .NOT.LSAME( TRANSA, 'C' ) ).AND.   &                                                                                                                                                   
               ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN                                                                                                                                                 
         INFO = 1                                                                                                                                                                                       
      ELSE IF( ( .NOT.NOTB                 ).AND.   &                                                                                                                                                   
               ( .NOT.LSAME( TRANSB, 'C' ) ).AND.   &                                                                                                                                                   
               ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN                                                                                                                                                 
         INFO = 2                                                                                                                                                                                       
      ELSE IF( M  .LT.0               )THEN                                                                                                                                                             
         INFO = 3                                                                                                                                                                                       
      ELSE IF( N  .LT.0               )THEN                                                                                                                                                             
         INFO = 4                                                                                                                                                                                       
      ELSE IF( K  .LT.0               )THEN                                                                                                                                                             
         INFO = 5                                                                                                                                                                                       
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN                                                                                                                                                             
         INFO = 8                                                                                                                                                                                       
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN                                                                                                                                                             
         INFO = 10                                                                                                                                                                                      
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN                                                                                                                                                             
         INFO = 13                                                                                                                                                                                      
      END IF                                                                                                                                                                                            
      IF( INFO.NE.0 )THEN                                                                                                                                                                               
         CALL XERBLA( 'DGEMM ', INFO )                                                                                                                                                                  
         RETURN                                                                                                                                                                                         
      END IF                                                                                                                                                                                            
!                                                                                                                                                                                                       
!     Quick return if possible.                                                                                                                                                                         
!                                                                                                                                                                                                       
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.   &                                                                                                                                                              
          ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )   &                                                                                                                             
         RETURN                                                                                                                                                                                         
!                                                                                                                                                                                                       
!     And if  alpha.eq.zero.                                                                                                                                                                            
!                                                                                                                                                                                                       
      IF( ALPHA.EQ.ZERO )THEN                                                                                                                                                                           
         IF( BETA.EQ.ZERO )THEN                                                                                                                                                                         
            DO 20, J = 1, N                                                                                                                                                                             
               DO 10, I = 1, M                                                                                                                                                                          
                  C( I, J ) = ZERO                                                                                                                                                                      
   10          CONTINUE                                                                                                                                                                                 
   20       CONTINUE                                                                                                                                                                                    
         ELSE                                                                                                                                                                                           
            DO 40, J = 1, N                                                                                                                                                                             
               DO 30, I = 1, M                                                                                                                                                                          
                  C( I, J ) = BETA*C( I, J )                                                                                                                                                            
   30          CONTINUE                                                                                                                                                                                 
   40       CONTINUE                                                                                                                                                                                    
         END IF                                                                                                                                                                                         
         RETURN                                                                                                                                                                                         
      END IF                                                                                                                                                                                            
!                                                                                                                                                                                                       
!     Start the operations.                                                                                                                                                                             
!                                                                                                                                                                                                       
      IF( NOTB )THEN                                                                                                                                                                                    
         IF( NOTA )THEN                                                                                                                                                                                 
!                                                                                                                                                                                                       
!           Form  C := alpha*A*B + beta*C.                                                                                                                                                              
!                                                                                                                                                                                                       
            DO 90, J = 1, N                                                                                                                                                                             
               IF( BETA.EQ.ZERO )THEN                                                                                                                                                                   
                  DO 50, I = 1, M                                                                                                                                                                       
                     C( I, J ) = ZERO                                                                                                                                                                   
   50             CONTINUE                                                                                                                                                                              
               ELSE IF( BETA.NE.ONE )THEN                                                                                                                                                               
                  DO 60, I = 1, M                                                                                                                                                                       
                     C( I, J ) = BETA*C( I, J )                                                                                                                                                         
   60             CONTINUE                                                                                                                                                                              
               END IF                                                                                                                                                                                   
               DO 80, L = 1, K                                                                                                                                                                          
                  IF( B( L, J ).NE.ZERO )THEN                                                                                                                                                           
                     TEMP = ALPHA*B( L, J )                                                                                                                                                             
                     DO 70, I = 1, M                                                                                                                                                                    
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )                                                                                                                                          
   70                CONTINUE                                                                                                                                                                           
                  END IF                                                                                                                                                                                
   80          CONTINUE                                                                                                                                                                                 
   90       CONTINUE                                                                                                                                                                                    
         ELSE                                                                                                                                                                                           
!                                                                                                                                                                                                       
!           Form  C := alpha*A'*B + beta*C                                                                                                                                                              
!                                                                                                                                                                                                       
            DO 120, J = 1, N                                                                                                                                                                            
               DO 110, I = 1, M                                                                                                                                                                         
                  TEMP = ZERO                                                                                                                                                                           
                  DO 100, L = 1, K                                                                                                                                                                      
                     TEMP = TEMP + A( L, I )*B( L, J )                                                                                                                                                  
  100             CONTINUE                                                                                                                                                                              
                  IF( BETA.EQ.ZERO )THEN                                                                                                                                                                
                     C( I, J ) = ALPHA*TEMP                                                                                                                                                             
                  ELSE                                                                                                                                                                                  
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )                                                                                                                                            
                  END IF                                                                                                                                                                                
  110          CONTINUE                                                                                                                                                                                 
  120       CONTINUE                                                                                                                                                                                    
         END IF                                                                                                                                                                                         
      ELSE                                                                                                                                                                                              
         IF( NOTA )THEN                                                                                                                                                                                 
!                                                                                                                                                                                                       
!           Form  C := alpha*A*B' + beta*C                                                                                                                                                              
!                                                                                                                                                                                                       
            DO 170, J = 1, N                                                                                                                                                                            
               IF( BETA.EQ.ZERO )THEN                                                                                                                                                                   
                  DO 130, I = 1, M                                                                                                                                                                      
                     C( I, J ) = ZERO                                                                                                                                                                   
  130             CONTINUE                                                                                                                                                                              
               ELSE IF( BETA.NE.ONE )THEN                                                                                                                                                               
                  DO 140, I = 1, M                                                                                                                                                                      
                     C( I, J ) = BETA*C( I, J )                                                                                                                                                         
  140             CONTINUE                                                                                                                                                                              
               END IF                                                                                                                                                                                   
               DO 160, L = 1, K                                                                                                                                                                         
                  IF( B( J, L ).NE.ZERO )THEN                                                                                                                                                           
                     TEMP = ALPHA*B( J, L )                                                                                                                                                             
                     DO 150, I = 1, M                                                                                                                                                                   
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )                                                                                                                                          
  150                CONTINUE                                                                                                                                                                           
                  END IF                                                                                                                                                                                
  160          CONTINUE                                                                                                                                                                                 
  170       CONTINUE                                                                                                                                                                                    
         ELSE                                                                                                                                                                                           
!                                                                                                                                                                                                       
!           Form  C := alpha*A'*B' + beta*C                                                                                                                                                             
!                                                                                                                                                                                                       
            DO 200, J = 1, N                                                                                                                                                                            
               DO 190, I = 1, M                                                                                                                                                                         
                  TEMP = ZERO                                                                                                                                                                           
                  DO 180, L = 1, K                                                                                                                                                                      
                     TEMP = TEMP + A( L, I )*B( J, L )                                                                                                                                                  
  180             CONTINUE                                                                                                                                                                              
                  IF( BETA.EQ.ZERO )THEN                                                                                                                                                                
                     C( I, J ) = ALPHA*TEMP                                                                                                                                                             
                  ELSE                                                                                                                                                                                  
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )                                                                                                                                            
                  END IF                                                                                                                                                                                
  190          CONTINUE                                                                                                                                                                                 
  200       CONTINUE                                                                                                                                                                                    
         END IF                                                                                                                                                                                         
      END IF                                                                                                                                                                                            
!                                                                                                                                                                                                       
      RETURN                                                                                                                                                                                            
!                                                                                                                                                                                                       
!     End of DGEMM .                                                                                                                                                                                    
!                                                                                                                                                                                                       
      END SUBROUTINE                                                                                                                                                                                    
!----------------------------------------------------------------------|       
      subroutine daxpy(n,da,dx,incx,dy,incy)                                                                                                                                                            
!                                                                                                                                                                                                       
!     constant times a vector plus a vector.                                                                                                                                                            
!     uses unrolled loops for increments equal to one.                                                                                                                                                  
!     jack dongarra, linpack, 3/11/78.                                                                                                                                                                  
!                                                                                                                                                                                                       
      double precision dx(1),dy(1),da                                                                                                                                                                   
      integer i,incx,incy,ix,iy,m,mp1,n                                                                                                                                                                 
!                                                                                                                                                                                                       
      if(n.le.0)return                                                                                                                                                                                  
      if (da .eq. 0.0d0) return                                                                                                                                                                         
      if(incx.eq.1.and.incy.eq.1)go to 20                                                                                                                                                               
!                                                                                                                                                                                                       
!        code for unequal increments or equal increments                                                                                                                                                
!          not equal to 1                                                                                                                                                                               
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      iy = 1                                                                                                                                                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                                                                                                                                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                                                                                                                                                 
      do 10 i = 1,n                                                                                                                                                                                     
        dy(iy) = dy(iy) + da*dx(ix)                                                                                                                                                                     
        ix = ix + incx                                                                                                                                                                                  
        iy = iy + incy                                                                                                                                                                                  
   10 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for both increments equal to 1                                                                                                                                                            
!                                                                                                                                                                                                       
!                                                                                                                                                                                                       
!        clean-up loop                                                                                                                                                                                  
!                                                                                                                                                                                                       
   20 m = mod(n,4)                                                                                                                                                                                      
      if( m .eq. 0 ) go to 40                                                                                                                                                                           
      do 30 i = 1,m                                                                                                                                                                                     
        dy(i) = dy(i) + da*dx(i)                                                                                                                                                                        
   30 continue                                                                                                                                                                                          
      if( n .lt. 4 ) return                                                                                                                                                                             
   40 mp1 = m + 1                                                                                                                                                                                       
      do 50 i = mp1,n,4                                                                                                                                                                                 
        dy(i) = dy(i) + da*dx(i)                                                                                                                                                                        
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)                                                                                                                                                            
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)                                                                                                                                                            
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)                                                                                                                                                            
   50 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end subroutine  daxpy
!----------------------------------------------------------------------|                                                                                                                                                                                               
      subroutine DGESV( N, M, A,LDA, IPIV, B,LDB, IFLAG )                                                                                                                                               
      integer N, M, LDA, LDB, IPIV(N), IFLAG                                                                                                                                                            
      double precision A(LDA,N), B(LDB,M)                                                                                                                                                               
      call DGEFA( A,LDA, N, IPIV, IFLAG )                                                                                                                                                               
      if ( IFLAG.ne.0 ) stop "Error in DGESV (LU factorisation)"                                                                                                                                        
      do j = 1,M                                                                                                                                                                                        
         call DGESL( A,LDA, N, IPIV,B(1,j), 0 )                                                                                                                                                         
      enddo                                                                                                                                                                                             
      end subroutine dgesv
!----------------------------------------------------------------------|                                                                                                                                                                                   
      subroutine  dscal(n,da,dx,incx)                                                                                                                                                                   
!                                                                                                                                                                                                       
!     scales a vector by a constant.                                                                                                                                                                    
!     uses unrolled loops for increment equal to one.                                                                                                                                                   
!     jack dongarra, linpack, 3/11/78.                                                                                                                                                                  
!     modified 3/93 to return if incx .le. 0.                                                                                                                                                           
!                                                                                                                                                                                                       
      double precision da,dx(1)                                                                                                                                                                         
      integer i,incx,m,mp1,n,nincx                                                                                                                                                                      
!                                                                                                                                                                                                       
      if( n.le.0 .or. incx.le.0 )return                                                                                                                                                                 
      if(incx.eq.1)go to 20                                                                                                                                                                             
!                                                                                                                                                                                                       
!        code for increment not equal to 1                                                                                                                                                              
!                                                                                                                                                                                                       
      nincx = n*incx                                                                                                                                                                                    
      do 10 i = 1,nincx,incx                                                                                                                                                                            
        dx(i) = da*dx(i)                                                                                                                                                                                
   10 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for increment equal to 1                                                                                                                                                                  
!                                                                                                                                                                                                       
!                                                                                                                                                                                                       
!        clean-up loop                                                                                                                                                                                  
!                                                                                                                                                                                                       
   20 m = mod(n,5)                                                                                                                                                                                      
      if( m .eq. 0 ) go to 40                                                                                                                                                                           
      do 30 i = 1,m                                                                                                                                                                                     
        dx(i) = da*dx(i)                                                                                                                                                                                
   30 continue                                                                                                                                                                                          
      if( n .lt. 5 ) return                                                                                                                                                                             
   40 mp1 = m + 1                                                                                                                                                                                       
      do 50 i = mp1,n,5                                                                                                                                                                                 
        dx(i) = da*dx(i)                                                                                                                                                                                
        dx(i + 1) = da*dx(i + 1)                                                                                                                                                                        
        dx(i + 2) = da*dx(i + 2)                                                                                                                                                                        
        dx(i + 3) = da*dx(i + 3)                                                                                                                                                                        
        dx(i + 4) = da*dx(i + 4)                                                                                                                                                                        
   50 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end subroutine dscal

!----------------------------------------------------------------------|                                                                                                                                
      subroutine dgefa(a,lda,n,ipvt,info)                                                                                                                                                               
      integer lda,n,ipvt(n),info                                                                                                                                                                        
      double precision a(lda,n)                                                                                                                                                                         
!                                                                                                                                                                                                       
!     dgefa factors a double precision matrix by gaussian elimination.                                                                                                                                  
!                                                                                                                                                                                                       
!     dgefa is usually called by dgeco, but it can be called                                                                                                                                            
!     directly with a saving in time if  rcond  is not needed.                                                                                                                                          
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .                                                                                                                                                   
!                                                                                                                                                                                                       
!     on entry                                                                                                                                                                                          
!                                                                                                                                                                                                       
!        a       double precision(lda, n)                                                                                                                                                               
!                the matrix to be factored.                                                                                                                                                             
!                                                                                                                                                                                                       
!        lda     integer                                                                                                                                                                                
!                the leading dimension of the array  a .                                                                                                                                                
!                                                                                                                                                                                                       
!        n       integer                                                                                                                                                                                
!                the order of the matrix  a .                                                                                                                                                           
!                                                                                                                                                                                                       
!     on return                                                                                                                                                                                         
!                                                                                                                                                                                                       
!        a       an upper triangular matrix and the multipliers                                                                                                                                         
!                which were used to obtain it.                                                                                                                                                          
!                the factorization can be written  a = l*u  where                                                                                                                                       
!                l  is a product of permutation and unit lower                                                                                                                                          
!                triangular matrices and  u  is upper triangular.                                                                                                                                       
!                                                                                                                                                                                                       
!        ipvt    integer(n)                                                                                                                                                                             
!                an integer vector of pivot indices.                                                                                                                                                    
!                                                                                                                                                                                                       
!        info    integer                                                                                                                                                                                
!                = 0  normal value.                                                                                                                                                                     
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error                                                                                                                                       
!                     condition for this subroutine, but it does                                                                                                                                        
!                     indicate that dgesl or dgedi will divide by zero                                                                                                                                  
!                     if called.  use  rcond  in dgeco for a reliable                                                                                                                                   
!                     indication of singularity.                                                                                                                                                        
!                                                                                                                                                                                                       
!     linpack. this version dated 08/14/78 .                                                                                                                                                            
!     cleve moler, university of new mexico, argonne national lab.                                                                                                                                      
!                                                                                                                                                                                                       
!     subroutines and functions                                                                                                                                                                         
!                                                                                                                                                                                                       
!     blas daxpy,dscal,idamax                                                                                                                                                                           
!                                                                                                                                                                                                       
!     internal variables                                                                                                                                                                                
!                                                                                                                                                                                                       
      double precision t                                                                                                                                                                                
      integer j,k,kp1,l,nm1                                                                                                                                                                      
!                                                                                                                                                                                                       
!                                                                                                                                                                                                       
!     gaussian elimination with partial pivoting                                                                                                                                                        
!                                                                                                                                                                                                       
      info = 0                                                                                                                                                                                          
      nm1 = n - 1                                                                                                                                                                                       
      if (nm1 .lt. 1) go to 70                                                                                                                                                                          
      do 60 k = 1, nm1                                                                                                                                                                                  
         kp1 = k + 1                                                                                                                                                                                    
!                                                                                                                                                                                                       
!        find l = pivot index                                                                                                                                                                           
!                                                                                                                                                                                                       
         l = idamax(n-k+1,a(k,k),1) + k - 1                                                                                                                                                             
         ipvt(k) = l                                                                                                                                                                                    
!                                                                                                                                                                                                       
!        zero pivot implies this column already triangularized                                                                                                                                          
!                                                                                                                                                                                                       
         if (a(l,k) .eq. 0.0d0) go to 40                                                                                                                                                                
!                                                                                                                                                                                                       
!           interchange if necessary                                                                                                                                                                    
!                                                                                                                                                                                                       
            if (l .eq. k) go to 10                                                                                                                                                                      
               t = a(l,k)                                                                                                                                                                               
               a(l,k) = a(k,k)                                                                                                                                                                          
               a(k,k) = t                                                                                                                                                                               
   10       continue                                                                                                                                                                                    
!                                                                                                                                                                                                       
!           compute multipliers                                                                                                                                                                         
!                                                                                                                                                                                                       
            t = -1.0d0/a(k,k)                                                                                                                                                                           
            call dscal(n-k,t,a(k+1,k),1)                                                                                                                                                                
!                                                                                                                                                                                                       
!           row elimination with column indexing                                                                                                                                                        
!                                                                                                                                                                                                       
            do 30 j = kp1, n                                                                                                                                                                            
               t = a(l,j)                                                                                                                                                                               
               if (l .eq. k) go to 20                                                                                                                                                                   
                  a(l,j) = a(k,j)                                                                                                                                                                       
                  a(k,j) = t                                                                                                                                                                            
   20          continue                                                                                                                                                                                 
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)                                                                                                                                                  
   30       continue                                                                                                                                                                                    
         go to 50                                                                                                                                                                                       
   40    continue                                                                                                                                                                                       
            info = k                                                                                                                                                                                    
   50    continue                                                                                                                                                                                       
   60 continue                                                                                                                                                                                          
   70 continue                                                                                                                                                                                          
      ipvt(n) = n                                                                                                                                                                                       
      if (a(n,n) .eq. 0.0d0) info = n                                                                                                                                                                   
      return                                                                                                                                                                                            
      end subroutine dgefa    
!----------------------------------------------------------------------|                                                                                                                                
      subroutine dgesl(a,lda,n,ipvt,b,job)                                                                                                                                                              
      integer lda,n,ipvt(n),job                                                                                                                                                                         
      double precision a(lda,n),b(n)                                                                                                                                                                    
!                                                                                                                                                                                                       
!     dgesl solves the double precision system                                                                                                                                                          
!     a * x = b  or  trans(a) * x = b                                                                                                                                                                   
!     using the factors computed by dgeco or dgefa.                                                                                                                                                     
!                                                                                                                                                                                                       
!     on entry                                                                                                                                                                                          
!                                                                                                                                                                                                       
!        a       double precision(lda, n)                                                                                                                                                               
!                the output from dgeco or dgefa.                                                                                                                                                        
!                                                                                                                                                                                                       
!        lda     integer                                                                                                                                                                                
!                the leading dimension of the array  a .                                                                                                                                                
!                                                                                                                                                                                                       
!        n       integer                                                                                                                                                                                
!                the order of the matrix  a .                                                                                                                                                           
!                                                                                                                                                                                                       
!        ipvt    integer(n)                                                                                                                                                                             
!                the pivot vector from dgeco or dgefa.                                                                                                                                                  
!                                                                                                                                                                                                       
!        b       double precision(n)                                                                                                                                                                    
!                the right hand side vector.                                                                                                                                                            
!                                                                                                                                                                                                       
!        job     integer                                                                                                                                                                                
!                = 0         to solve  a*x = b ,                                                                                                                                                        
!                = nonzero   to solve  trans(a)*x = b  where                                                                                                                                            
!                            trans(a)  is the transpose.                                                                                                                                                
!                                                                                                                                                                                                       
!     on return                                                                                                                                                                                         
!                                                                                                                                                                                                       
!        b       the solution vector  x .                                                                                                                                                               
!                                                                                                                                                                                                       
!     error condition                                                                                                                                                                                   
!                                                                                                                                                                                                       
!        a division by zero will occur if the input factor contains a                                                                                                                                   
!        zero on the diagonal.  technically this indicates singularity                                                                                                                                  
!        but it is often caused by improper arguments or improper                                                                                                                                       
!        setting of lda .  it will not occur if the subroutines are                                                                                                                                     
!        called correctly and if dgeco has set rcond .gt. 0.0                                                                                                                                           
!        or dgefa has set info .eq. 0 .                                                                                                                                                                 
!                                                                                                                                                                                                       
!     to compute  inverse(a) * c  where  c  is a matrix                                                                                                                                                 
!     with  p  columns                                                                                                                                                                                  
!           call dgeco(a,lda,n,ipvt,rcond,z)                                                                                                                                                            
!           if (rcond is too small) go to ...                                                                                                                                                           
!           do 10 j = 1, p                                                                                                                                                                              
!              call dgesl(a,lda,n,ipvt,c(1,j),0)                                                                                                                                                        
!        10 continue                                                                                                                                                                                    
!                                                                                                                                                                                                       
!     linpack. this version dated 08/14/78 .                                                                                                                                                            
!     cleve moler, university of new mexico, argonne national lab.                                                                                                                                      
!                                                                                                                                                                                                       
!     subroutines and functions                                                                                                                                                                         
!                                                                                                                                                                                                       
!     blas daxpy,ddot                                                                                                                                                                                   
!                                                                                                                                                                                                       
!     internal variables                                                                                                                                                                                
!                                                                                                                                                                                                       
      double precision t                                                                                                                                                                           
      integer k,kb,l,nm1                                                                                                                                                                                
!                                                                                                                                                                                                       
      nm1 = n - 1                                                                                                                                                                                       
      if (job .ne. 0) go to 50                                                                                                                                                                          
!                                                                                                                                                                                                       
!        job = 0 , solve  a * x = b                                                                                                                                                                     
!        first solve  l*y = b                                                                                                                                                                           
!                                                                                                                                                                                                       
         if (nm1 .lt. 1) go to 30                                                                                                                                                                       
         do 20 k = 1, nm1                                                                                                                                                                               
            l = ipvt(k)                                                                                                                                                                                 
            t = b(l)                                                                                                                                                                                    
            if (l .eq. k) go to 10                                                                                                                                                                      
               b(l) = b(k)                                                                                                                                                                              
               b(k) = t                                                                                                                                                                                 
   10       continue                                                                                                                                                                                    
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)                                                                                                                                                       
   20    continue                                                                                                                                                                                       
   30    continue                                                                                                                                                                                       
!                                                                                                                                                                                                       
!        now solve  u*x = y                                                                                                                                                                             
!                                                                                                                                                                                                       
         do 40 kb = 1, n                                                                                                                                                                                
            k = n + 1 - kb                                                                                                                                                                              
            b(k) = b(k)/a(k,k)                                                                                                                                                                          
            t = -b(k)                                                                                                                                                                                   
            call daxpy(k-1,t,a(1,k),1,b(1),1)                                                                                                                                                           
   40    continue                                                                                                                                                                                       
      go to 100                                                                                                                                                                                         
   50 continue                                                                                                                                                                                          
!                                                                                                                                                                                                       
!        job = nonzero, solve  trans(a) * x = b                                                                                                                                                         
!        first solve  trans(u)*y = b                                                                                                                                                                    
!                                                                                                                                                                                                       
         do 60 k = 1, n                                                                                                                                                                                 
            t = ddot(k-1,a(1,k),1,b(1),1)                                                                                                                                                               
            b(k) = (b(k) - t)/a(k,k)                                                                                                                                                                    
   60    continue                                                                                                                                                                                       
!                                                                                                                                                                                                       
!        now solve trans(l)*x = y                                                                                                                                                                       
!                                                                                                                                                                                                       
         if (nm1 .lt. 1) go to 90                                                                                                                                                                       
         do 80 kb = 1, nm1                                                                                                                                                                              
            k = n - kb                                                                                                                                                                                  
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)                                                                                                                                                 
            l = ipvt(k)                                                                                                                                                                                 
            if (l .eq. k) go to 70                                                                                                                                                                      
               t = b(l)                                                                                                                                                                                 
               b(l) = b(k)                                                                                                                                                                              
               b(k) = t                                                                                                                                                                                 
   70       continue                                                                                                                                                                                    
   80    continue                                                                                                                                                                                       
   90    continue                                                                                                                                                                                       
  100 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end subroutine dgesl
!----------------------------------------------------------------------|                                                                                                                                
      integer function idamax(n,dx,incx)                                                                                                                                                                
!                                                                                                                                                                                                       
!     finds the index of element having max. absolute value.                                                                                                                                            
!     jack dongarra, linpack, 3/11/78.                                                                                                                                                                  
!     modified 3/93 to return if incx .le. 0.                                                                                                                                                           
!                                                                                                                                                                                                       
      double precision dx(1),dmax                                                                                                                                                                       
      integer i,incx,ix,n                                                                                                                                                                               
!                                                                                                                                                                                                       
      idamax = 0                                                                                                                                                                                        
      if( n.lt.1 .or. incx.le.0 ) return                                                                                                                                                                
      idamax = 1                                                                                                                                                                                        
      if(n.eq.1)return                                                                                                                                                                                  
      if(incx.eq.1)go to 20                                                                                                                                                                             
!                                                                                                                                                                                                       
!        code for increment not equal to 1                                                                                                                                                              
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      dmax = dabs(dx(1))                                                                                                                                                                                
      ix = ix + incx                                                                                                                                                                                    
      do 10 i = 2,n                                                                                                                                                                                     
         if(dabs(dx(ix)).le.dmax) go to 5                                                                                                                                                               
         idamax = i                                                                                                                                                                                     
         dmax = dabs(dx(ix))                                                                                                                                                                            
    5    ix = ix + incx                                                                                                                                                                                 
   10 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for increment equal to 1                                                                                                                                                                  
!                                                                                                                                                                                                       
   20 dmax = dabs(dx(1))                                                                                                                                                                                
      do 30 i = 2,n                                                                                                                                                                                     
         if(dabs(dx(i)).le.dmax) go to 30                                                                                                                                                               
         idamax = i                                                                                                                                                                                     
         dmax = dabs(dx(i))                                                                                                                                                                             
   30 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end function idamax                                                                                                                                                                                                                                        
!----------------------------------------------------------------------|                                                                                                                                
      double precision function ddot(n,dx,incx,dy,incy)                                                                                                                                                 
!                                                                                                                                                                                                       
!     forms the dot product of two vectors.                                                                                                                                                             
!     uses unrolled loops for increments equal to one.                                                                                                                                                  
!     jack dongarra, linpack, 3/11/78.                                                                                                                                                                  
!                                                                                                                                                                                                       
      double precision dx(1),dy(1),dtemp                                                                                                                                                                
      integer i,incx,incy,ix,iy,m,mp1,n                                                                                                                                                                 
!                                                                                                                                                                                                       
      ddot = 0.0d0                                                                                                                                                                                      
      dtemp = 0.0d0                                                                                                                                                                                     
      if(n.le.0)return                                                                                                                                                                                  
      if(incx.eq.1.and.incy.eq.1)go to 20                                                                                                                                                               
!                                                                                                                                                                                                       
!        code for unequal increments or equal increments                                                                                                                                                
!          not equal to 1                                                                                                                                                                               
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      iy = 1                                                                                                                                                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                                                                                                                                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                                                                                                                                                 
      do 10 i = 1,n                                                                                                                                                                                     
        dtemp = dtemp + dx(ix)*dy(iy)                                                                                                                                                                   
        ix = ix + incx                                                                                                                                                                                  
        iy = iy + incy                                                                                                                                                                                  
   10 continue                                                                                                                                                                                          
      ddot = dtemp                                                                                                                                                                                      
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for both increments equal to 1                                                                                                                                                            
!                                                                                                                                                                                                       
!                                                                                                                                                                                                       
!        clean-up loop                                                                                                                                                                                  
!                                                                                                                                                                                                       
   20 m = mod(n,5)                                                                                                                                                                                      
      if( m .eq. 0 ) go to 40                                                                                                                                                                           
      do 30 i = 1,m                                                                                                                                                                                     
        dtemp = dtemp + dx(i)*dy(i)                                                                                                                                                                     
   30 continue                                                                                                                                                                                          
      if( n .lt. 5 ) go to 60                                                                                                                                                                           
   40 mp1 = m + 1                                                                                                                                                                                       
      do 50 i = mp1,n,5                                                                                                                                                                                 
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +   &                                                                                                                                         
         dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)                                                                                                                                
   50 continue                                                                                                                                                                                          
   60 ddot = dtemp                                                                                                                                                                                      
      return                                                                                                                                                                                            
      end function ddot 

!-----------------------------------------------------------------------------
      subroutine  dswap (n,dx,incx,dy,incy)                                                                                                                                                             
!                                                                                                                                                                                                       
!     interchanges two vectors.                                                                                                                                                                         
!     uses unrolled loops for increments equal one.                                                                                                                                                     
!     jack dongarra, linpack, 3/11/78.                                                                                                                                                                  
!                                                                                                                                                                                                       
      double precision dx(1),dy(1),dtemp                                                                                                                                                                
      integer i,incx,incy,ix,iy,m,mp1,n                                                                                                                                                                 
!                                                                                                                                                                                                       
      if(n.le.0)return                                                                                                                                                                                  
      if(incx.eq.1.and.incy.eq.1)go to 20                                                                                                                                                               
!                                                                                                                                                                                                       
!       code for unequal increments or equal increments not equal                                                                                                                                       
!         to 1                                                                                                                                                                                          
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      iy = 1                                                                                                                                                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                                                                                                                                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                                                                                                                                                 
      do 10 i = 1,n                                                                                                                                                                                     
        dtemp = dx(ix)                                                                                                                                                                                  
        dx(ix) = dy(iy)                                                                                                                                                                                 
        dy(iy) = dtemp                                                                                                                                                                                  
        ix = ix + incx                                                                                                                                                                                  
        iy = iy + incy                                                                                                                                                                                  
   10 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!       code for both increments equal to 1                                                                                                                                                             
!                                                                                                                                                                                                       
!                                                                                                                                                                                                       
!       clean-up loop                                                                                                                                                                                   
!                                                                                                                                                                                                       
   20 m = mod(n,3)                                                                                                                                                                                      
      if( m .eq. 0 ) go to 40                                                                                                                                                                           
      do 30 i = 1,m                                                                                                                                                                                     
        dtemp = dx(i)                                                                                                                                                                                   
        dx(i) = dy(i)                                                                                                                                                                                   
        dy(i) = dtemp                                                                                                                                                                                   
   30 continue                                                                                                                                                                                          
      if( n .lt. 3 ) return                                                                                                                                                                             
   40 mp1 = m + 1                                                                                                                                                                                       
      do 50 i = mp1,n,3                                                                                                                                                                                 
        dtemp = dx(i)                                                                                                                                                                                   
        dx(i) = dy(i)                                                                                                                                                                                   
        dy(i) = dtemp                                                                                                                                                                                   
        dtemp = dx(i + 1)                                                                                                                                                                               
        dx(i + 1) = dy(i + 1)                                                                                                                                                                           
        dy(i + 1) = dtemp                                                                                                                                                                               
        dtemp = dx(i + 2)                                                                                                                                                                               
        dx(i + 2) = dy(i + 2)                                                                                                                                                                           
        dy(i + 2) = dtemp                                                                                                                                                                               
   50 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end subroutine  dswap
!-------------------------------------------------------------------------
      integer function izamax(n,zx,incx)                                                                                                                                                                
!                                                                                                                                                                                                       
!     finds the index of element having max. absolute value.                                                                                                                                            
!     jack dongarra, 1/15/85.                                                                                                                                                                           
!     modified 3/93 to return if incx .le. 0.                                                                                                                                                           
!                                                                                                                                                                                                       
      double complex zx(1)                                                                                                                                                                              
      double precision smax                                                                                                                                                                             
      integer i,incx,ix,n                                                                                                                                                                               
!      double precision dcabs1                                                                                                                                                                           
!                                                                                                                                                                                                       
      izamax = 0                                                                                                                                                                                        
      if( n.lt.1 .or. incx.le.0 )return                                                                                                                                                                 
      izamax = 1                                                                                                                                                                                        
      if(n.eq.1)return                                                                                                                                                                                  
      if(incx.eq.1)go to 20                                                                                                                                                                             
!                                                                                                                                                                                                       
!        code for increment not equal to 1                                                                                                                                                              
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      smax = dcabs1(zx(1))                                                                                                                                                                              
      ix = ix + incx                                                                                                                                                                                    
      do 10 i = 2,n                                                                                                                                                                                     
         if(dcabs1(zx(ix)).le.smax) go to 5                                                                                                                                                             
         izamax = i                                                                                                                                                                                     
         smax = dcabs1(zx(ix))                                                                                                                                                                          
    5    ix = ix + incx                                                                                                                                                                                 
   10 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for increment equal to 1                                                                                                                                                                  
!                                                                                                                                                                                                       
   20 smax = dcabs1(zx(1))                                                                                                                                                                              
      do 30 i = 2,n                                                                                                                                                                                     
         if(dcabs1(zx(i)).le.smax) go to 30                                                                                                                                                             
         izamax = i                                                                                                                                                                                     
         smax = dcabs1(zx(i))                                                                                                                                                                           
   30 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end function izamax
!-------------------------------------------------------------------------------------                                                                                                                                                                                     
      subroutine  zscal(n,za,zx,incx)                                                                                                                                                                   
!                                                                                                                                                                                                       
!     scales a vector by a constant.                                                                                                                                                                    
!     jack dongarra, 3/11/78.                                                                                                                                                                           
!     modified 3/93 to return if incx .le. 0.                                                                                                                                                           
!                                                                                                                                                                                                       
      double complex za,zx(1)                                                                                                                                                                           
      integer i,incx,ix,n                                                                                                                                                                               
!                                                                                                                                                                                                       
      if( n.le.0 .or. incx.le.0 )return                                                                                                                                                                 
      if(incx.eq.1)go to 20                                                                                                                                                                             
!                                                                                                                                                                                                       
!        code for increment not equal to 1                                                                                                                                                              
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      do 10 i = 1,n                                                                                                                                                                                     
        zx(ix) = za*zx(ix)                                                                                                                                                                              
        ix = ix + incx                                                                                                                                                                                  
   10 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for increment equal to 1                                                                                                                                                                  
!                                                                                                                                                                                                       
   20 do 30 i = 1,n                                                                                                                                                                                     
        zx(i) = za*zx(i)                                                                                                                                                                                
   30 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end subroutine zscal
!----------------------------------------------------------------------|                                                                                                                                
      subroutine zaxpy(n,za,zx,incx,zy,incy)                                                                                                                                                            
!                                                                                                                                                                                                       
!     constant times a vector plus a vector.                                                                                                                                                            
!     jack dongarra, 3/11/78.                                                                                                                                                                           
!                                                                                                                                                                                                       
      complex*16 zx(1),zy(1),za                                                                                                                                                                         
!      double precision dcabs1                                                                                                                                                                           
      if(n.le.0)return                                                                                                                                                                                  
      if (dcabs1(za) .eq. 0.0d0) return                                                                                                                                                                 
      if (incx.eq.1.and.incy.eq.1)go to 20                                                                                                                                                              
!                                                                                                                                                                                                       
!        code for unequal increments or equal increments                                                                                                                                                
!          not equal to 1                                                                                                                                                                               
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      iy = 1                                                                                                                                                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                                                                                                                                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                                                                                                                                                 
      do 10 i = 1,n                                                                                                                                                                                     
        zy(iy) = zy(iy) + za*zx(ix)                                                                                                                                                                     
        ix = ix + incx                                                                                                                                                                                  
        iy = iy + incy                                                                                                                                                                                  
   10 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for both increments equal to 1                                                                                                                                                            
!                                                                                                                                                                                                       
   20 do 30 i = 1,n                                                                                                                                                                                     
        zy(i) = zy(i) + za*zx(i)                                                                                                                                                                        
   30 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end subroutine zaxpy                                                                                                                                                                                   
!---------------------------------------------------------------------------                                                                                                                                                                                   
      double complex function zdotc(n,zx,incx,zy,incy)                                                                                                                                                  
!                                                                                                                                                                                                       
!     forms the dot product of a vector.                                                                                                                                                                
!     jack dongarra, 3/11/78.                                                                                                                                                                           
!                                                                                                                                                                                                       
      double complex zx(1),zy(1),ztemp                                                                                                                                                                  
      ztemp = (0.0d0,0.0d0)                                                                                                                                                                             
      zdotc = (0.0d0,0.0d0)                                                                                                                                                                             
      if(n.le.0)return                                                                                                                                                                                  
      if(incx.eq.1.and.incy.eq.1)go to 20                                                                                                                                                               
!                                                                                                                                                                                                       
!        code for unequal increments or equal increments                                                                                                                                                
!          not equal to 1                                                                                                                                                                               
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      iy = 1                                                                                                                                                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                                                                                                                                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                                                                                                                                                 
      do 10 i = 1,n                                                                                                                                                                                     
        ztemp = ztemp + dconjg(zx(ix))*zy(iy)                                                                                                                                                           
        ix = ix + incx                                                                                                                                                                                  
        iy = iy + incy                                                                                                                                                                                  
   10 continue                                                                                                                                                                                          
      zdotc = ztemp                                                                                                                                                                                     
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for both increments equal to 1                                                                                                                                                            
!                                                                                                                                                                                                       
   20 do 30 i = 1,n                                                                                                                                                                                     
        ztemp = ztemp + dconjg(zx(i))*zy(i)                                                                                                                                                             
   30 continue                                                                                                                                                                                          
      zdotc = ztemp                                                                                                                                                                                     
      return                                                                                                                                                                                            
      end function zdotc  
!------------------------------------------------------------------------------                                                                                                                                                                                   
      subroutine  zdscal(n,da,zx,incx)                                                                                                                                                                  
!                                                                                                                                                                                                       
!     scales a vector by a constant.                                                                                                                                                                    
!     jack dongarra, 3/11/78.                                                                                                                                                                           
!     modified 3/93 to return if incx .le. 0.                                                                                                                                                           
!                                                                                                                                                                                                       
      double complex zx(1)                                                                                                                                                                              
      double precision da                                                                                                                                                                               
      integer i,incx,ix,n                                                                                                                                                                               
!                                                                                                                                                                                                       
      if( n.le.0 .or. incx.le.0 )return                                                                                                                                                                 
      if(incx.eq.1)go to 20                                                                                                                                                                             
!                                                                                                                                                                                                       
!        code for increment not equal to 1                                                                                                                                                              
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      do 10 i = 1,n                                                                                                                                                                                     
        zx(ix) = dcmplx(da,0.0d0)*zx(ix)                                                                                                                                                                
        ix = ix + incx                                                                                                                                                                                  
   10 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for increment equal to 1                                                                                                                                                                  
!                                                                                                                                                                                                       
   20 do 30 i = 1,n                                                                                                                                                                                     
        zx(i) = dcmplx(da,0.0d0)*zx(i)                                                                                                                                                                  
   30 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end subroutine zdscal                                                                                                                                                                                  
!-----------------------------------------------------------------------------
      subroutine  zswap (n,zx,incx,zy,incy)                                                                                                                                                             
!                                                                                                                                                                                                       
!     interchanges two vectors.                                                                                                                                                                         
!     jack dongarra, 3/11/78.                                                                                                                                                                           
!                                                                                                                                                                                                       
      complex*16 zx(1),zy(1),ztemp                                                                                                                                                                      
!                                                                                                                                                                                                       
      if(n.le.0)return                                                                                                                                                                                  
      if(incx.eq.1.and.incy.eq.1)go to 20                                                                                                                                                               
!                                                                                                                                                                                                       
!       code for unequal increments or equal increments not equal                                                                                                                                       
!         to 1                                                                                                                                                                                          
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      iy = 1                                                                                                                                                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                                                                                                                                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                                                                                                                                                 
      do 10 i = 1,n                                                                                                                                                                                     
        ztemp = zx(ix)                                                                                                                                                                                  
        zx(ix) = zy(iy)                                                                                                                                                                                 
        zy(iy) = ztemp                                                                                                                                                                                  
        ix = ix + incx                                                                                                                                                                                  
        iy = iy + incy                                                                                                                                                                                  
   10 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!       code for both increments equal to 1                                                                                                                                                             
   20 do 30 i = 1,n                                                                                                                                                                                     
        ztemp = zx(i)                                                                                                                                                                                   
        zx(i) = zy(i)                                                                                                                                                                                   
        zy(i) = ztemp                                                                                                                                                                                   
   30 continue                                                                                                                                                                                          
      return                                                                                                                                                                                            
      end subroutine zswap                                                     
!--------------------------------------------------------------------------------                                 
      double complex function zdotu(n,zx,incx,zy,incy)                                                                                                                                                  
!                                                                                                                                                                                                       
!     forms the dot product of a vector.                                                                                                                                                                
!     jack dongarra, 3/11/78.                                                                                                                                                                           
!                                                                                                                                                                                                       
      double complex zx(1),zy(1),ztemp                                                                                                                                                                  
      ztemp = (0.0d0,0.0d0)                                                                                                                                                                             
      zdotu = (0.0d0,0.0d0)                                                                                                                                                                             
      if(n.le.0)return                                                                                                                                                                                  
      if(incx.eq.1.and.incy.eq.1)go to 20                                                                                                                                                               
!                                                                                                                                                                                                       
!        code for unequal increments or equal increments                                                                                                                                                
!          not equal to 1                                                                                                                                                                               
!                                                                                                                                                                                                       
      ix = 1                                                                                                                                                                                            
      iy = 1                                                                                                                                                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                                                                                                                                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                                                                                                                                                 
      do 10 i = 1,n                                                                                                                                                                                     
        ztemp = ztemp + zx(ix)*zy(iy)                                                                                                                                                                   
        ix = ix + incx                                                                                                                                                                                  
        iy = iy + incy                                                                                                                                                                                  
   10 continue                                                                                                                                                                                          
      zdotu = ztemp                                                                                                                                                                                     
      return                                                                                                                                                                                            
!                                                                                                                                                                                                       
!        code for both increments equal to 1                                                                                                                                                            
!                                                                                                                                                                                                       
   20 do 30 i = 1,n                                                                                                                                                                                     
        ztemp = ztemp + zx(i)*zy(i)                                                                                                                                                                     
   30 continue                                                                                                                                                                                          
      zdotu = ztemp                                                                                                                                                                                     
      return                                                                                                                                                                                            
      end function zdotu
!---------------------------------------------------------------------------------


      end module