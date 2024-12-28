      SUBROUTINE HYPOELASTICITY_clay(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      use tensor_tools
      
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
      
      real(8) :: Tb(3,3), epor, epsb(3,3), p
      real(8) :: kappa, nu, K, G
      real(8) :: Idev(3,3,3,3), EE(3,3,3,3), Jac(3,3,3,3)
	
      ! Defining material parameters   
      kappa   = props(1)
      nu = props(2)
      
      ! Initialization of state variables    
      epor = statev(1)
      Tb   = map2T(stress,ntens)
      epsb = map2D(dstran,ntens) 
      
      p = -tr(Tb)/3.0d0
      K = p*(1.0d0+(epor))/kappa           
      G = K*3.0d0*(1.0d0-2.0d0*(nu))/(2.0d0*(1.0d0+(nu)))      
      
      Idev=Idelta-1.0d0/3.0d0*((delta).out.(delta))
      
      EE = 3.0d0*K*(unit(delta).out.unit(delta)) + 2.0d0*G*(Idev)   
      
      Tb = Tb + (EE.xx.epsb)
      Jac = EE
      
      ! Update void ratio
      epor = epor +(1.0d0+epor)*tr(epsb) 
      
      ! Export state variables  
      stress = map2stress(Tb,ntens)
      ddsdde = map2ddsdde(Jac,ntens)
      statev(1) = epor  
	
      RETURN
      
      end subroutine HYPOELASTICITY_clay 