! ====================================================================
! 
!  MATERIAL MODELS
!
!  - Linear-Elasticity
!  - Hypo-Elasticity_c
!  - Hypo-Elasticity_s
!  - Mod-Cam-Clay
!
! ====================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
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


      CHARACTER*50 :: model

      ! ===================================================
      
      model = trim(CMNAME)
      model= TRIM(ADJUSTL(model))
      
	IF (index(model,'Linear-Elasticity') .GT. 0) THEN
	  
	  IF (NPROPS .NE. 2) THEN 
	    WRITE(*,*) 'ERROR: NPROPS not correct for Linear-Elastic model.'
	    STOP
	  ENDIF 
	
	  call LINEARELASTICITY(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      ELSE IF (index(model,'Hypo-Elasticity_c') .GT. 0) THEN
	  
	  IF (NPROPS .NE. 2) THEN 
	    WRITE(*,*) 'ERROR: NPROPS not correct for Hypo-Elastic clay model.'
	    STOP
	  ENDIF 
	
	  call HYPOELASTICITY_clay(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)  
	  
      ELSE IF (index(model,'Hypo-Elasticity_s') .GT. 0) THEN
	  
	  IF (NPROPS .NE. 3) THEN 
	    WRITE(*,*) 'ERROR: NPROPS not correct for Hypo-Elastic sand model.'
	    STOP
	  ENDIF 
	
	  call HYPOELASTICITY_sand(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)  
	  
      ELSE IF (index(model,'MCC') .GT. 0) THEN
	  
	  IF (NPROPS .NE. 4) THEN 
	    WRITE(*,*) 'ERROR: NPROPS not correct for Mod-Cam-Clay model.'
	    STOP
	  ENDIF 
	
	  call MCC(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)  

	ELSE
	  
	  WRITE(*,*) ' '
	  WRITE(*,*) 'ERROR: Material Model Identifier not correct.'
	  WRITE(*,*) ' '
	  WRITE(*,*) 'Choose one of the following models:'
	  WRITE(*,*) 'Linear-Elasticity'
	  STOP
	  
	ENDIF
      
      return
      end subroutine       
