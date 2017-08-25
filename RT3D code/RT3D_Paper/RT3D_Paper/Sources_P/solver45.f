C THIS IS A DUMMY PACKAGE FOR ISOLVER 4 AND 5.  To use solvers 4 or 5, this 
C   package should be replaced by the real "solver45.f" code.  Unfortunately, 
C   the real code cannot be distributed via Internet due to potential copyright 
C   issues.  Please see the RT3Dv25_Readme.txt file.  Requests for the real 
C   code may be sent to rt3d@pnl.gov.

      SUBROUTINE dodeint(ystart,nvar,isolver,x1,x2,eps,h1,hmin,nok,
     &                   nbad,fex,jac)

c*    This is a modified double precision version of an adaptive stepsize
c*    controlled Runge-Kutta and a backward stiff solver driver

      IMPLICIT NONE
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX,isolver
      INTEGER i,kmax,kount,nstp
      DOUBLE PRECISION eps,h1,hmin,x1,x2,TINY
      DOUBLE PRECISION ystart(nvar)
      EXTERNAL fex,jac
      PARAMETER (MAXSTP=2000,NMAX=100,KMAXX=1,TINY=1.e-30)
      DOUBLE PRECISION dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),
     *                 y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp

        write (*,*)
        write (*,*) "Solvers 4 and 5 are not available..."
        write (*,*) "  Replace solver45.f with a real file."
        write (*,*)
        write (*,*) "Read 'RT3Dv25_Readme.txt' or send a request"
        write (*,*) "  to rt3d@pnl.gov for more details."
        write (*,*)

C-------EMRL JIG
        call stopfile
C-------EMRL JIG
        stop

      RETURN
      END
