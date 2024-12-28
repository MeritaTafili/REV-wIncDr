C     Umat coded by William Fuentes and Merita Tafili 
C     Modified Cam Clay, Implicit integration
C     Inifinitesimal strains
      SUBROUTINE MCC(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
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

C ---------------------------------------------------------------     
        real(8) :: T(3,3), DEPS(3,3), JAC(3,3,3,3),
     1 evoid, epsP(3,3), pc, lambda, kappa, M
     2 ,nu, Kbulk, E, G, ELMOD(3,3,3,3), Idev(3,3,3,3), TrialSTRESS(3,3)
     3 , TrialNormS, TrialDevSTRESS(3,3), p , q
     4 , TrialNormSt,Trialp,Trialq, TrialF, DeltaPhi, theta
     5 ,TolG, TolF, Gpc, denom, Gpcderiv, Fyield, dfdp
     6 ,dfdq, dfdpc ,dpdDPhi, dqdDPhi, dpcdDPhi ,derivF
     7 ,normaldev(3,3) ,Depsp(3,3), NormDepsp
     8 , caa,ca1,ca2,ca3,ca4,ca5,ca6,cbb,cb1,cb2,tsi, NormSt, term1,pcn      
      integer maxiter, kcounter, j
C     
      CALL Initial(STRESS,T, DSTRAN, DEPS, NTENS,NDI, NSHR)       
C ---------------------------------------------------------------       
C     State variables 
C     Void ratio as state variable (scalar)
      evoid=STATEV(1)
C     Plastic strains as state variable (tensor R3xR3)
      CALL Vectortomatrix(STATEV(2:7),epsP) 
C     Overconsolidated mean stress state variable
      pc=abs(STATEV(8))
C ---------------------------------------------------------------      
C     Defining material parameters 
      lambda=PROPS(1) !compresion index
      kappa=PROPS(2) !swelling index
      M=PROPS(3) !critic state line slope
      nu=PROPS(4) !poisson modulus    
C ---------------------------------------------------------------
C     Mean stress
      p=-Tr(T)/3.0d0   
      if (abs(p)<1.0d-3) then 
      p=1.0d-3
      endif
      if (abs(p)>abs(pc)) then 
      pc=p
      endif
C --------------------------------------------------------------- 
C     ELASTIC MODULUS 
C     Material constants         
      Kbulk=abs((1.0d0+evoid)*p)/(kappa) !bulk modulus
      if (abs(Kbulk)<1.0d0) then 
      Kbulk=1.0d0
      endif 
      if (nu==0.5d0) then
      nu=0.4999999999999999d0
      endif
      E=Kbulk*3.0d0*(1.0d0-2.0d0*nu)!young modulus
      G=E/(2.0d0*(1.0d0+nu))!shear modulus
      
      CALL ELMOD1(E, nu, ELMOD) 
C ---------------------------------------------------------------  
C ---------------------------------------------------------------  
C     Trial Elastic trial step
C     Deviator fourth unit tensor
      CALL Idev1(Idev) 
C     Trial stress 
      TrialSTRESS=(ELMOD.xx.DEPS)+T
C     Trial deviator stress  
      TrialDevSTRESS=Idev.xx.TrialSTRESS   
C     Norm of trial stress 
      TrialNormS=norm(TrialSTRESS)  
      TrialNormSt= norm(TrialDevSTRESS)    
      call pq(TrialSTRESS,Trialp,Trialq)    
C     Trial yield function F
      TrialF=(Trialq**2.0d0/M**2.0d0)+(Trialp)*((Trialp)-(pc))
C ---------------------------------------------------------------       
      if (TrialF<0.0d0) then
C     Elastic Steps, Trial step is ok
      T=TrialSTRESS
      JAC=ELMOD
C     Void ratio 
      evoid=evoid+(1.0d0+evoid)*(Tr(Deps))
C --------------------------------------------------------------- 
C     Saving State variables 
C     Alpha isotropic as state variable (scalar)
      STATEV(1)=evoid
C     Plastic strains as state variable (tensor R3xR3)
      CALL Matrixtovector(STATEV(2:7),epsp)
      STATEV(8)=pc
C ---------------------------------------------------------------       
      else ! Plastic corrector step

      DeltaPhi=0.0d0     !initializing consistency parameter
      evoid=evoid+(1.0d0+evoid)*(tr(Deps))
      theta=(1.0d0+evoid)/(lambda-kappa) !state variable (for softening hardening law)
      pcn=pc !preconsolidation pressure at the beginning of the iterations 

C ---------------------------------------------------------------
C     Defining tolerances for G and F and maximum iterations
      TolG=1.0E-7*pcn           ! Tolerance for preconsolidation
      TolF=1.0E-7*TrialNormS    ! Tolerance for yield surface
      maxiter=50                ! maximum iterations
C ---------------------------------------------------------------
C     Cycle for iterating with DeltaPhi 
      do kcounter=1,maxiter
C ---------------------------------------------------------------
C     Cycle for iterating with pc   
      do j=0,maxiter  
C     Target function
      Gpc=pcn*exp(theta*DeltaPhi*(2.0d0*(Trialp)-(pc))/
     1 (1.0d0+2.0d0*DeltaPhi*Kbulk))-(pc)
      denom=2.0d0*Kbulk*DeltaPhi !constant 
C     Derivative of G(abs(pc)) with respect to time 
      Gpcderiv=-exp((-((pc)-2.0d0*(Trialp))*DeltaPhi*theta)/
     1 (denom+1.0d0))*pcn*DeltaPhi*Theta/(denom+1.0d0)-1.0d0 
C     Next value for abs(pc)     
      pc=pc-Gpc/Gpcderiv
      if (abs(Gpc)<TolG) then
      exit
      endif
      enddo
C --------------------------------------------------------------- 
C --------------------------------------------------------------- 
C     Current q 
      q=Trialq/(1.0d0+6.0d0*G*DeltaPhi/M**2.0d0)
C     Current p      
      p=((Trialp)+DeltaPhi*Kbulk*(pc))/(1.0d0+2.0d0*DeltaPhi*Kbulk)
C     Current yield function F     
      Fyield=q**2.0d0/M**2.0d0+(p)*((p)-(pc))      
C     Calculating derivatives of F with respect to p,q,pc
      dfdp=2.0d0*(p)-(pc)
      dfdq=2.0d0*q/M**2.0d0
      dfdpc=-(p)
C     Calculating derivatives of p ,q,and pc with respect to DeltaPhi (consistency parameter)
      dpdDPhi=-Kbulk*dfdp/(1.0d0+(2.0d0*Kbulk+theta*(pc))*DeltaPhi)
      dqdDPhi=-q/(DeltaPhi+M**2.0d0/(6.0d0*G))
      dpcdDPhi=theta*(pc)*dfdp/(1.0d0+(2.0d0*Kbulk+
     1 theta*(pc))*DeltaPhi)    
C     Calculating derivative of F with respect to DeltaPhi
      derivF=dfdp*dpdDPhi+dfdq*dqdDPhi+dfdpc*dpcdDPhi  
C     Calculating new DeltaPhi for next iteration
      DeltaPhi=DeltaPhi-Fyield/derivF  
      if (abs(Fyield)<TolF) then
      exit         
      endif  
      enddo 
	  
C     End of iterations 
C ---------------------------------------------------------------
C     Increment of plastic strains DepsP
C ---------------------------------------------------------------   
C     Delta plastic strains
      if  (TrialNormSt/=0.0d0) then
C     Flow rule direction
      normaldev=TrialDevSTRESS/TrialNormSt    
      Depsp=(DeltaPhi*(-(2.0d0*p-pc)*1.0d0/3.0d0*delta+sqrt(3.0d0/2.0d0)
     1 *(2.0d0*q/M**2.)*normaldev)) 
      NormDepsp=norm(Depsp)
      else
      Depsp=0.0d0
      endif
C ---------------------------------------------------------------
C     Stress
      T=TrialSTRESS-(ELMOD.xx.Depsp)
C     Plastic strain 
      epsP=epsP+Depsp
C     Void ratio 
c      evoid=evoid+(1.0d0+evoid)*Tr(Deps)
      call pq(T,p,q) 
C --------------------------------------------------------------- 
C     Saving State variables 
C     Alpha isotropic as state variable (scalar)
      STATEV(1)=evoid
C     Plastic strains as state variable (tensor R3xR3)
      CALL Matrixtovector(STATEV(2:7),epsp)
      if (p.gt.0.0d0) then
	  term1=(q/M)**2.0d0/p+p
	  endif
	  pc=term1
	  STATEV(8)=pc
C ---------------------------------------------------------------
      normaldev=-normaldev
C ---------------------------------------------------------------  
C     Consistent elastoplastic modulus
      if  (NormDepsp/=0) then       
      caa=1.0d0+2.0d0*Kbulk*DeltaPhi+abs(pc)*theta*DeltaPhi !constant
      ca1=(1.0d0+abs(pc)*theta*DeltaPhi)/caa !constant
      ca2=-(2.0d0*abs(p)-abs(pc))/caa !constant
      ca3=2.0d0*abs(pc)*theta*DeltaPhi/caa !constant
      ca4=theta*abs(pc)*(2.0d0*abs(p)-abs(pc))/(caa*Kbulk) !constant
      ca5=sqrt(3.0d0/2.0d0)*(1.0d0+6.0d0*G*DeltaPhi/M**2.0d0)**(-1.0d0) !constant
      ca6=-3.0d0*q/M**2.0d0*(1.0d0+6.0d0*G*DeltaPhi/M**2.0d0)**(-1.0d0) !constant
      cbb=-2.0d0*G*2.0d0*q*ca6/M**2.0d0-Kbulk*((2.0d0*ca2-ca4)*abs(p)
     1 -ca2*abs(pc))!const
      cb1=-Kbulk*((ca3-2.0d0*ca1)*abs(p)+ca1*abs(pc))/cbb !constant
      cb2=2.0d0*G*2.0d0*q/M**2.0d0*ca5/cbb !constant
      NormSt=q*sqrt(2.0d0/3.0d0)  
      tsi=NormSt/TrialNormSt !flow rule
C     Consistent elastoplastic modulus    
      JAC=2.0d0*G*tsi*(Idelta)+(Kbulk*(ca1+ca2*cb1)-1.0d0
     1 /3.0d0*2.0d0*G*tsi)*(delta.out.delta)+Kbulk*ca2*cb2*
     2 (delta.out.normaldev)+2.0d0*G*sqrt(2.0d0/3.0d0)
     3 *(ca6*cb1)*(normaldev.out.delta)+2.0d0*G
     4 *(sqrt(2.0d0/3.0d0)*(ca5+ca6*cb2)-tsi)*
     5 (normaldev.out.normaldev)   
       else
      JAC=ELMOD
       endif  
       endif
C ---------------------------------------------------------------
      Call Solution(NTENS, NDI, NSHR, T, STRESS, JAC
     1 , DDSDDE) 
	   if (isnan(stress(1))) then
	   term1=1.0d0
      endif
      
      END SUBROUTINE MCC
      
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Initial(STRESS,T, DSTRAN, DEPS, NTENS,NDI, NSHR)      
      double precision STRESS(ntens), T(3,3)
     1 ,DSTRAN(ntens), DEPS(3,3) 
      Integer ntens, nshr, ndi   
      DEPS=0.0D0
      T=0.0D0
C     
      do i=1,ndi
      T(i,i)=stress(i)
      DEPS(i,i)=DSTRAN(i)
      enddo 
C     
      if (nshr.ge.1) then
      T(1,2)=stress(4)
      T(2,1)=stress(4)    
      DEPS(1,2)=0.5d0*DSTRAN(4)
      DEPS(2,1)=0.5d0*DSTRAN(4)           
      endif
      if (nshr.ge.2) then
      T(1,3)=stress(5)
      T(3,1)=stress(5)   
      DEPS(1,3)=0.5d0*DSTRAN(5)
      DEPS(3,1)=0.5d0*DSTRAN(5)          
      endif
      if (nshr.ge.3) then
      T(2,3)=stress(6)
      T(3,2)=stress(6)    
      DEPS(2,3)=0.5d0*DSTRAN(6)
      DEPS(3,2)=0.5d0*DSTRAN(6)         
      endif   
      return          
      END SUBROUTINE Initial  
      
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Matrixtovector(vec1,mat1)      
      double precision vec1(6),  mat1(3,3)   
      vec1(1)=mat1(1,1)  
      vec1(2)=mat1(2,2) 
      vec1(3)=mat1(3,3) 
      vec1(4)=mat1(1,2)  
      vec1(5)=mat1(1,3) 
      vec1(6)=mat1(2,3)  
      return          
      END SUBROUTINE Matrixtovector
c------------------------------------------------------ 
c------------------------------------------------------ 
       SUBROUTINE Solution(NTENS, NDI, NSHR, T, STRESS, JAC, DDSDDE) 
       integer NTENS, NDI, NSHR, i, j, k, l
C     Subroutine for filling the stress and Jacobian matrix        
       double precision T(3,3), JAC(3,3,3,3), STRESS(NTENS),
     1 DDSDDE(NTENS,NTENS), JAC66(6,6)
      k=1
      l=1
c------------------------------------------------------
      do i=1,ndi
      stress(i)=T(i,i)
      enddo 
C     
      if (nshr.ge.1) then
      stress(ndi+1)=T(1,2)         
      endif
      if (nshr.ge.2) then
      stress(ndi+2)=T(1,3)         
      endif
      if (nshr.ge.3) then
      stress(ndi+3)=T(2,3)         
      endif   
      call tensortomatrix(jac,  jac66)   
        do i=1,ndi
        do j=1,ndi
          ddsdde(i,j)=jac66(i,j)
        enddo
      enddo 
      do i=ndi+1,ndi+nshr
        do j=1,ndi
          ddsdde(i,j)=jac66(3+k,j)
        enddo
        k=k+1
      enddo  
      do i=1,ndi
      l=1
        do j=ndi+1,ndi+nshr
          ddsdde(i,j)=jac66(i,3+l)
          l=l+1
        enddo
      enddo   
      k=1
      do i=ndi+1,ndi+nshr
        l=1
        do j=ndi+1,ndi+nshr
          ddsdde(i,j)=jac66(3+k,3+l)
          l=l+1
        enddo
        k=k+1
      enddo   
       Return       
      END SUBROUTINE Solution 
c------------------------------------------------------ 
c------------------------------------------------------       
      SUBROUTINE Vectortomatrix(vec1,mat1)      
      double precision vec1(6),  mat1(3,3)   
      mat1(1,1)=vec1(1)  
      mat1(2,2)=vec1(2)  
      mat1(3,3)=vec1(3)  
      mat1(1,2)=vec1(4)  
      mat1(2,1)=vec1(4)  
      mat1(1,3)=vec1(5)  
      mat1(3,1)=vec1(5)  
      mat1(2,3)=vec1(6)  
      mat1(3,2)=vec1(6)  
      return          
      END SUBROUTINE Vectortomatrix             
c------------------------------------------------------ 
c------------------------------------------------------ 
      subroutine tensortomatrix(a3333,  b66)   ! returns b(6,6)
      double precision a3333(3,3,3,3),b66(6,6)
      integer i,j,i9(6),j9(6)
      data i9/1,2,3,1,1,2/
     .     j9/1,2,3,2,3,3/
      do  i=1,6   !  switch to matrix notation
      do  j=1,6
      b66(i,j)=a3333(i9(i),j9(i),i9(j),j9(j))  
      enddo
      enddo   
      return
      end subroutine tensortomatrix
c------------------------------------------------------ 
c------------------------------------------------------ 
      subroutine pq(T,p,q)   ! p and q
      use tensor_tools
      double precision T(3,3),p,q, DevSTRESS(3,3),
     1 NormSt, Idev(3,3,3,3)
C      mean stress 
      p=-1.0d0/3.0d0*Tr(T)
C     deviator stress
      Call idev1(Idev)
      DevSTRESS=dev(T)
C     Norm of the deviator stress 
      NormSt=norm(DevSTRESS)
C     Trial q 
      q=sqrt(3.0d0/2.0d0)*NormSt    
      return
      end subroutine pq
      
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Idev1(Idev) 
      use tensor_tools
      DOUBLE PRECISION Idev(3,3,3,3), 
     1  Ivol(3,3,3,3), Isym(3,3,3,3)
      INTEGER I,J,K,L
      Idev=0.0d0
      Ivol=0.0d0
      Isym=0.0D0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      Isym(i,j,k,l)=1.0d0/2.0d0*(delta(i,k)*delta(j,l)+
     1 delta(i,l)*delta(j,k))              
      Enddo
      Enddo
      Enddo
      Enddo
      Ivol=0.0D0
      Ivol=1.0d0/3.0d0*delta.out.delta    
      Idev=Isym-Ivol         
      END SUBROUTINE Idev1
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE ELMOD1(E, nu, ELMOD) 
      use tensor_tools
      DOUBLE PRECISION ELMOD(3,3,3,3), E, nu,
     1 Idev(3,3,3,3), Ivol(3,3,3,3), Isym(3,3,3,3)
      INTEGER I,J,K,L
      Idev=0.0d0
      Ivol=0.0d0
      ELMOD=0.0D0
      Isym=0.0D0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      Isym(i,j,k,l)=1.0d0/2.0d0*(delta(i,k)*delta(j,l)+
     1 delta(i,l)*delta(j,k))              
      Enddo
      Enddo
      Enddo
      Enddo
      Ivol=0.0D0
      Ivol=1.0d0/3.0d0*delta.out.delta    
      Idev=Isym-Ivol         
      ELMOD=E/(1.0D0-2.0D0*nu)*Ivol+E/(1.0d0+nu)*Idev
      END SUBROUTINE ELMOD1  