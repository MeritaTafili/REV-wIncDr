      module tensor_tools
      ! taken from niemunis_tools.for
      implicit none
      
      real(8), parameter,dimension(3,3) ::
     &                  delta =RESHAPE([1,0,0,0,1,0,0,0,1],[3,3])  ! \com Kronecker $\delta_{ij}$
      real(8), parameter,dimension(3,3,3,3) ::
     &         Idelta= RESHAPE([ 2,0,0,0,0,0,0,0,0,0,               !  \com $ I_{ijkl} = \frac12(\delta_{ik} *\delta_{jl} + \delta_{il} *\delta_{jk}) $
     &                  1,0,1,0,0,0,0,0,0,0,     1,0,0,0,1,0,0,0,1,0,
     &                  1,0,0,0,0,0,0,0,0,0,     2,0,0,0,0,0,0,0,0,0,
     &                  1,0,1,0,0,0,1,0,0,0,     1,0,0,0,0,0,0,0,1,0,
     &                  1,0,0,0,0,0,0,0,0,0,     2],[3,3,3,3]) /2.0d0
      integer, parameter,dimension(1:6) :: i6=[1,2,3,1,1,2],         !  \com aux. tables for transformations tensor <--> matrix
     &                                     j6=[1,2,3,2,3,3]
      interface norm                                          ! \com double contraction, e.g. $\cL:\Db$, use (), interfaces have lowest priority
         module procedure norm33
      end interface
      
      interface dev                                          ! \com double contraction, e.g. $\cL:\Db$, use (), interfaces have lowest priority
         module procedure dev33
      end interface
      
      interface unit                                        ! \com double contraction, e.g. $\cL:\Db$, use (), interfaces have lowest priority
         module procedure unit33
      end interface
      
      interface operator(.xx.)                                          ! \com double contraction, e.g. $\cL:\Db$, use (), interfaces have lowest priority
         module procedure mal, mal2,mal3,mal4
       end interface
       
      interface operator(.out.)                                          ! \com double contraction, e.g. $\cL:\Db$, use (), interfaces have lowest priority
         module procedure outmal3333, outmal33 
         end interface 

      INTERFACE Inverse    
      MODULE PROCEDURE Inverse3333, inverse33
      END INTERFACE
      
      INTERFACE tpose23    
      MODULE PROCEDURE tpose23
      END INTERFACE
      
      INTERFACE tpose24 
      MODULE PROCEDURE tpose24
      END INTERFACE
       
       
      
      private
      public map2stran,map2D,map2stress,map2T,map2ddsdde,norm,dev,tr,
     &       Idelta,delta, operator(.xx.) , operator(.out.), unit, 
     &       inverse,tpose23,tpose24,Iunit1, mb
      
      contains
      !------------------------------------------------------------------------
      ! taken from niemunis: unsymmetric_module.for
      function map2stran(a,ntens)                                       ! \com ..........{\large MAP2STRAN}
      
        implicit none                                                   ! \com converts D(3,3)  to stran(6) with $\gamma_{12} = 2 \epsilon_{12}$ etc.
        real(8), intent(in), dimension(1:3,1:3) :: a
        integer, intent(in) :: ntens
        real(8),  dimension(1:ntens) :: map2stran
        integer :: i
        map2stran(1)=a(1,1)
        map2stran(2)=a(2,2)
        map2stran(3)=a(3,3)
        do i=4,ntens
        map2stran(i) = a(i6(i),j6(i)) * 2.0d0                           ! \com  abaqus needs gammas  as strains
        enddo                                                           ! \com   i6=(/ 1,2,3, 1,1,2/), j6=(/ 1,2,3, 2,3,3/)
      end function map2stran
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! taken from niemunis: unsymmetric_module.for
      function map2D(a,ntens)                                           ! \com ..........{\large MAP2D}
        implicit none                                                   ! \com convert strain rate from vector dstran(1:ntens) to  D(3,3)
        real(8),  dimension(1:3,1:3) :: map2D
        integer, intent(in) :: ntens
        real(8), intent(in), dimension(:) :: a
        integer :: i
        map2D=0
        map2D(1,1) = a(1)
        map2D(2,2) = a(2)
        map2D(3,3) = a(3)
        do i=4,ntens
         map2D(i6(i),j6(i))=a(i)/2.0d0                                   ! \com  abaqus uses  gammas  for shear strains
         map2D(j6(i),i6(i))=a(i)/2.0d0                                   ! \com   i6=(/ 1,2,3, 1,1,2/), j6=(/ 1,2,3, 2,3,3/)
        enddo
      end function map2D
      !------------------------------------------------------------------------
      ! 
      function mb(a)                            
        implicit none   
        real(8) :: mb                                               
        real(8), intent(in) :: a

        if (a.ge.0.0d0) then
        mb=a
        else
        mb=0.0d0
        endif                                                     
      end function mb
      
      
            
      !------------------------------------------------------------------------
      ! taken from niemunis: unsymmetric_module.for
      function map2stress(a,ntens)                                      ! \com ..........{\large MAP2STRESS}
        implicit none                                                   ! \com convert tensor T(3,3)  to matrix stress(ntens)
        real(8), intent(in), dimension(1:3,1:3) :: a
        integer, intent(in) :: ntens
        real(8),  dimension(1:ntens) :: map2stress
        integer :: i
        do i=1,ntens
          map2stress(i) = a(i6(i),j6(i))                                ! \com  abaqus stress
        enddo                                                           ! \com   i6=(/ 1,2,3, 1,1,2/), j6=(/ 1,2,3, 2,3,3/)
      end function map2stress
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! taken from niemunis: unsymmetric_module.for
      function map2T(a,ntens)                                           ! \com ..........{\large MAP2T}
        implicit none                                                   ! \com  convert  matrix stress(1:ntens)  to tensor T(3,3)
        real(8),  dimension(1:3,1:3) :: map2T
        integer, intent(in) :: ntens
        real(8), intent(in), dimension(:) :: a
        integer :: i
        map2T=0
        map2T(1,1) = a(1)
        map2T(2,2) = a(2)
        map2T(3,3) = a(3)
        do i=4,ntens
        map2T(i6(i),j6(i))=a(i)
        map2T(j6(i),i6(i))=a(i)
        enddo
      end function map2T
      !------------------------------------------------------------------------
        
c------------------------------------------------------ 
      function map2ddsdde(LL,ntens)                                     ! \com ..........{\large MAP2DDSDDE}
        implicit none                                                   ! \com  convert stiffness  3333 tensor to  abaqus's matrix ddsdde(ntens,ntens)
        real(8), intent(in), dimension(1:3,1:3,1:3,1:3) :: LL
        integer, intent(in) :: ntens
        real(8),  dimension(1:ntens,1:ntens) :: map2ddsdde
        integer :: i,j
        do i=1,ntens
        do j=1,ntens
          if (j <= 3) map2ddsdde(i,j) = LL(i6(i),j6(i),i6(j),j6(j))
          if (j >  3) map2ddsdde(i,j) =  0.5d0*
     &     (LL(i6(i),j6(i),i6(j),j6(j))+LL(i6(i),j6(i),j6(j),i6(j)) )
        enddo
        enddo
        end function map2ddsdde
        
c------------------------------------------------------ 
      function mal(a,b)                                               ! \com  usage:  (a.xx.b),  mind brackets
        implicit none
        real(8), intent(in), dimension(1:3,1:3)  :: a,b                 !\com            33.xx.33
        real(8) :: mal
         mal          =  a(1,1)*b(1,1)+
     &                    a(1,2)*b(1,2)+
     &                    a(1,3)*b(1,3)+
     &                    a(2,1)*b(2,1)+
     &                    a(2,2)*b(2,2)+
     &                    a(2,3)*b(2,3)+
     &                    a(3,1)*b(3,1)+
     &                    a(3,2)*b(3,2)+
     &                    a(3,3)*b(3,3)
         end function mal
        
c------------------------------------------------------ 
       function mal2(a,b)                                               ! \com  usage:  (a.xx.b),  mind brackets
       implicit none
           real(8), intent(in), dimension(1:3,1:3,1:3,1:3) :: a         !\com     3333.xx.33
           real(8), intent(in),  dimension(1:3,1:3) :: b
           real(8), dimension(1:3,1:3):: mal2
           integer :: i,j
          do  i=1,3
          do  j=1,3
          mal2(i,j)    =  a(i,j,1,1)*b(1,1)+
     &                    a(i,j,1,2)*b(1,2)+
     &                    a(i,j,1,3)*b(1,3)+
     &                    a(i,j,2,1)*b(2,1)+
     &                    a(i,j,2,2)*b(2,2)+
     &                    a(i,j,2,3)*b(2,3)+
     &                    a(i,j,3,1)*b(3,1)+
     &                    a(i,j,3,2)*b(3,2)+
     &                    a(i,j,3,3)*b(3,3)
          enddo
          enddo
        end function mal2
        
c------------------------------------------------------ 
       function mal3(a,b)                                               ! \com  usage:  (a.xx.b),  mind brackets
         implicit none
         real(8), intent(in),  dimension(1:3,1:3) :: a
         real(8), intent(in), dimension(1:3,1:3,1:3,1:3) :: b           !\com    33.xx.3333
         real(8), dimension(1:3,1:3):: mal3
         integer :: k,l
          do  k=1,3
          do  l=1,3
          mal3(k,l) =   a(1,1)*b(1,1,k,l)+
     &                  a(1,2)*b(1,2,k,l)+
     &                  a(1,3)*b(1,3,k,l)+
     &                  a(2,1)*b(2,1,k,l)+
     &                  a(2,2)*b(2,2,k,l)+
     &                  a(2,3)*b(2,3,k,l)+
     &                  a(3,1)*b(3,1,k,l)+
     &                  a(3,2)*b(3,2,k,l)+
     &                  a(3,3)*b(3,3,k,l)
          enddo
          enddo
       end function mal3
        
c------------------------------------------------------ 
       function mal4(a,b)                                               ! \com  usage:  (a.xx.b),  mind brackets
         implicit none
         real(8), intent(in), dimension(1:3,1:3,1:3,1:3):: a,b          !\com    3333.xx.3333
         real(8), dimension(1:3,1:3,1:3,1:3):: mal4
         integer :: i,j,k,l
          do  i=1,3
          do  j=1,3
          do  k=1,3
          do  l=1,3
          mal4(i,j,k,l)= a(i,j,1,1)*b(1,1,k,l)+
     &                  a(i,j,1,2)*b(1,2,k,l)+
     &                  a(i,j,1,3)*b(1,3,k,l)+
     &                  a(i,j,2,1)*b(2,1,k,l)+
     &                  a(i,j,2,2)*b(2,2,k,l)+
     &                  a(i,j,2,3)*b(2,3,k,l)+
     &                  a(i,j,3,1)*b(3,1,k,l)+
     &                  a(i,j,3,2)*b(3,2,k,l)+
     &                  a(i,j,3,3)*b(3,3,k,l)
          enddo
          enddo
          enddo
          enddo
        end function mal4
c------------------------------------------------------ 
        function outmal3333(a,b)                                            ! \com  usage:  (a.out.b) , mind brackets
        implicit none
        
        real(8), intent(in), dimension(1:3,1:3)  :: a,b
        real(8), dimension(1:3,1:3,1:3,1:3) :: outmal3333
        integer :: i,j,k,l
        do i=1,3
        do j=1,3
        do k=1,3
        do l=1,3
            outmal3333(i,j,k,l) =  a(i,j)*b(k,l)
        enddo
        enddo
        enddo
        enddo
        end function outmal3333
c------------------------------------------------------ 
        function outmal33(a,b)                                            ! \com  usage:  (a.out.b) , mind brackets
        implicit none
        real(8), intent(in), dimension(1:3)  :: a,b
        real(8), dimension(1:3,1:3) :: outmal33
        integer :: i,j
        do i=1,3
        do j=1,3
            outmal33(i,j) =  a(i)*b(j)
        enddo
        enddo
        end function outmal33
                
c------------------------------------------------------ 
       function norm33(a)                                             ! \com   usage:  norm(a)
         implicit none
          real(8), intent(in), dimension(3,3)  :: a
          real(8)  :: norm33
          norm33 =a(1,1)*a(1,1)+a(1,2)*a(1,2)+a(1,3)*a(1,3)+
     &          a(2,1)*a(2,1)+a(2,2)*a(2,2)+a(2,3)*a(2,3)+
     &          a(3,1)*a(3,1)+a(3,2)*a(3,2)+a(3,3)*a(3,3)
          if(norm33 > tiny(norm33)) norm33 = dsqrt(norm33)
       end function norm33 
        
c------------------------------------------------------ 
       function dev33(a)                                             ! \com   usage:  dev(a)
         implicit none
          real(8), intent(in), dimension(3,3)  :: a
          real(8)  :: dev33(3,3)
          dev33 = a-((1.0d0/3.0d0*(tr(a)))*delta)
       end function dev33  
             
c------------------------------------------------------ 
      function unit33(a)                                             ! \com   usage:  unit(a)
         implicit none
          real(8), intent(in), dimension(3,3)  :: a
          real(8)  :: unit33(3,3)
      if (norm(a).ne.0.0d0) then
      unit33=a/norm(a)
      else
      unit33=0.0d0    
      endif
      end function unit33
        
      
         function tr(a)                                                   ! \com ..........{\large TR}
          implicit none                                                  ! \com    returns  trace of a tensor
          real(8), intent(in), dimension(1:3,1:3)  :: a
          real(8) :: tr
          tr=a(1,1)+a(2,2)+a(3,3)
        end function tr
        
c------------------------------------------------------ 
c------------------------------------------------------
      FUNCTION inverse3333(a3333) RESULT (b3333) 
      DOUBLE PRECISION, INTENT(IN):: a3333(3,3,3,3)   
      DOUBLE PRECISION b3333(3,3,3,3)          
      DOUBLE PRECISION a(9,9),b(9,9),
     1 c(9,18),cc,dd
      integer i,j,k,inv,i9(9),j9(9)
      data i9/1,2,3,1,2,1,3,2,3/,
     .     j9/1,2,3,2,1,3,1,3,2/
      do 10 i=1,9   !  switch to matrix notation
      do 10 j=1,9
  10  a(i,j)=a3333(i9(i),j9(i),i9(j),j9(j))
c (1) Preparam. c(9,18)
      do i=1,9
      do j=1,9
      c(i,j+9)=0.0d0
      c(i,j)=a(i,j)
      enddo
      c(i,9+i)=1.0d0
      enddo
c (2) invert c with result  going to the right half of it
      do i=1,9
      cc = c(i,i)
      if (abs(cc).lt.1d-4) then
      b3333=0.0d0
      inv =0              ! inversion failed
      Print*, "UMAT: Inversion failed"
      pause
      CALL XIT
      return
      endif
      c(i,i)=cc-1.0d0
      do k=i+1,18
      dd=c(i,k)/cc
      do j=1,9
      c(j,k)=c(j,k)-dd*c(j,i)
      enddo
      enddo
      enddo
c (3) copy result to b()
      do i=1,9
      do j=1,9
      b(i,j)= c(i,j+9)
      enddo
      enddo
      inv=1                    ! inversion successfull
      do 20 i=1,9     ! switch to tensorial notation
      do 20 j=1,9
  20  b3333(i9(i),j9(i),i9(j),j9(j))=b(i,j)
      return
      END FUNCTION inverse3333
c------------------------------------------------------ 
c------------------------------------------------------ 
c     
      FUNCTION inverse33(a33) RESULT (b33)  
C
      integer m, n, i, j, k
      PARAMETER (M=3,N=3)
      Double precision a33(3,3),b33(3,3), p
C
      DO 5 I=1,M
      DO 5 J=1,N
      b33(I,J)=a33(I,J)
5     CONTINUE
      DO 10 K=1,M
      P=b33(K,K)
      b33(K,K)=1.0d0
      DO 20 J=1,N
      b33(K,J)=b33(K,J)/P
20    CONTINUE
      DO 10 I=1,M
      IF(I .EQ. K) GO TO 10
      P=b33(I,K)
      b33(I,K)=0.0d0
      DO 30 J=1,N
      b33(I,J)=b33(I,J)-b33(K,J)*P
30    CONTINUE
10    CONTINUE
      return
      END FUNCTION inverse33    
c-------------------------------      
c------------------------------------------------------ 
c------------------------------------------------------ 
      FUNCTION tpose23(a3333) RESULT (b3333)  
C
      integer m, n, p, q, i, j, k, l
      PARAMETER (m=3,n=3,p=3,q=3)
      Double precision a3333(3,3,3,3),b3333(3,3,3,3)
C
      do i=1,m
      do j=1,n
      do k=1,p
      do l=1,q
      b3333(i,j,k,l)=a3333(i,k,j,l)
      enddo
      enddo
      enddo
      enddo
      END FUNCTION tpose23    
c------------------------------------------------------ 
c------------------------------------------------------ 
      FUNCTION tpose24(a3333) RESULT (b3333)  
C
      integer m, n, p, q, i, j, k, l
      PARAMETER (m=3,n=3,p=3,q=3)
      Double precision a3333(3,3,3,3),b3333(3,3,3,3)
C
      do i=1,m
      do j=1,n
      do k=1,p
      do l=1,q
      b3333(i,j,k,l)=a3333(i,l,k,j)
      enddo
      enddo
      enddo
      enddo
      END FUNCTION tpose24   
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Iunit1(Iunit) 
      DOUBLE PRECISION Iunit(3,3,3,3)
      INTEGER I,J,K,L
      Iunit=0.0D0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      Iunit(i,j,k,l)=delta(i,k)*delta(j,l)           
      Enddo
      Enddo
      Enddo
      Enddo     
      END SUBROUTINE Iunit1     
c------------------------------------------------------ 
      end module tensor_tools