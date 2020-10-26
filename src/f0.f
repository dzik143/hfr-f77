C*****************************************************************
      FUNCTION F0(ARG)
C                                                                *
C     CALCULATES THE ERROR FUNCTION ACCORDING TO A RATIONAL      *
C     APPROXIMATION FROM M. ABRAMOWITZ AND I.A. STEGUN,          *
C     ABSOLUTE ERROR IS LESS THAN 1.5*10**(-7)                   *
C     CAN BE REPLACED BY A BUILT-IN FUNCTION ON SOME MACHINES    *
C                                                                *
C*****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(5)
      DATA P/0.3275911D0/
      DATA A/0.254829592D0,-0.284496736D0,1.421413741D0,
     $ -1.453152027D0,1.061405429D0/
      real*8 SARG
      real*8 ERF

      if (ARG.lt.1.0D-6) then
        F0=1.0D0-ARG/3.0D0
        return
      endif

      SARG = sqrt(ARG)
      T=1.0D0/(1.0D0+P*SARG)
      TN=T
      POLY=A(1)*TN
      DO 10 I=2,5
      TN=TN*T
      POLY=POLY+A(I)*TN
   10 CONTINUE
      ERF=1.0D0-POLY*DEXP(-ARG)

      F0 = 0.886226925452758014 / SARG * ERF
      RETURN
      END
C***************************************************************
