      SUBROUTINE LINEARELASTICITY(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      implicit none

      CHARACTER*80 CMNAME
	  
      integer :: NTENS,NSTATEV,NPROPS,NDI,NSHR,NOEL,
     &           NPT,LAYER,KSPT,KSTEP,KINC
	 
      real(8) :: SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP,
     &           PNEWDT,CELENT
	 
      real(8) :: STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
	real(8), PARAMETER :: ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, NINE=9.0D0
	real(8) :: E,ANU,ALAMBDA,AMU
	integer :: I,J
	

      ! add your code here
      E = props(1)
      ANU = props(2)
      
      ALAMBDA=E*ANU/(ONE+ANU)/(ONE-TWO*ANU)
	AMU=E/(ONE+ANU)/TWO
C
	DO I=1,NTENS
	 DO J=1,NTENS
	  DDSDDE(I,J)=ZERO
	 ENDDO
      ENDDO
      
      DDSDDE(1,1) = ALAMBDA+TWO*AMU
	DDSDDE(2,2)=ALAMBDA+TWO*AMU
	DDSDDE(3,3)=ALAMBDA+TWO*AMU
	DDSDDE(4,4)=AMU
	DDSDDE(5,5)=AMU
	DDSDDE(6,6)=AMU
	DDSDDE(1,2)=ALAMBDA
	DDSDDE(1,3)=ALAMBDA
	DDSDDE(2,3)=ALAMBDA
	DDSDDE(2,1)=ALAMBDA
	DDSDDE(3,1)=ALAMBDA
	DDSDDE(3,2)=ALAMBDA
C
      ! add your code here
C
	DO I=1,NTENS
	 DO J=1,NTENS
	  STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
	 ENDDO
	ENDDO
	
      RETURN
      END