! Jordan to Einstein Frame
!
!  The quantities transformed, along with their transformation factors are:
!
!	g_ij^E   =  A^-2(phi)  g_ij^J
!
!	alp^E    =  A^-1(phi) alp^J  
!       
!       rho^E    =  A^4(phi) rho^J
!
!       beta^iE  =  beta^iJ
!
!	eps^E    =  A^4(phi) eps^J
!
!	p^E      =  A^4(phi) p^J
!
!	u_mu^E   =  A^-1(phi) u_mu^J 
!
!       T_munu^E =  A^2(phi) T_munu^J

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine  JordanToEinstein( CCTK_ARGUMENTS )
   implicit none
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
      
   !Variables
   CCTK_REAL                A_conf, phiJ
   CCTK_REAL                metricJtoE, alpJtoE, hydroJtoE, velJtoE, tmunuJtoE
   CCTK_INT                 i, j, k
   CCTK_INT                 a, b, c
   character(len=200) :: message


! Allocate memory for the character array
! Format and store values in the character array
   if (EJverbose.eq.1) then 
        write(message, *) rho(0,0,0) 

         call  CCTK_INFO("rho_center in J frame")
         call  CCTK_INFO(message)
         call  CCTK_INFO("going to Einstein frame")
   end if 
   !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(k,j,i, a,b,c, A_conf, metricJtoE, alpJtoE,hydroJtoE, velJtoE, tmunuJtoE)
       
   do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
     do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
       do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

        !---Get local variables--------
        phiJ = phi1(i,j,k)  
       
        !---Calculate conformal factor A(phi) 
        
        if (decoupling_limit.eq.1) then 
                A_conf = 1.0
        else        
                A_conf = EXP(0.5 * betaDEF * phiJ * phiJ   +  k0BD * phiJ)
        end if 

        AJtoE(i,j,k) = A_conf
                
        !---Conversion factors
        metricJtoE = 1/(A_conf*A_conf)
        alpJtoE      = 1/A_conf
        tmunuJtoE   = A_conf*A_conf

        !----Assign values to grid function
        gxx(i,j,k)   =  metricJtoE * gxx(i,j,k)   
        gxy(i,j,k)   =  metricJtoE * gxy(i,j,k) 
        gxz(i,j,k)   =  metricJtoE * gxz(i,j,k) 
        gyy(i,j,k)   =  metricJtoE * gyy(i,j,k) 
        gyz(i,j,k)   =  metricJtoE * gyz(i,j,k) 
        gzz(i,j,k)   =  metricJtoE * gzz(i,j,k) 
        alp(i,j,k)   =  alpJtoE * alp(i,j,k)


        if (decoupling_limit.eq.1) then


                eTtt(i,j,k)  =  tmunuJtoE * eTtt(i,j,k)   
                eTtx(i,j,k)  =  tmunuJtoE * eTtx(i,j,k)   
                eTty(i,j,k)  =  tmunuJtoE * eTty(i,j,k)   
                eTtz(i,j,k)  =  tmunuJtoE * eTtz(i,j,k)   
                eTxx(i,j,k)  =  tmunuJtoE * eTxx(i,j,k)   
                eTxy(i,j,k)  =  tmunuJtoE * eTxy(i,j,k)   
                eTxz(i,j,k)  =  tmunuJtoE * eTxz(i,j,k)   
                eTyy(i,j,k)  =  tmunuJtoE * eTyy(i,j,k)   
                eTyz(i,j,k)  =  tmunuJtoE * eTyz(i,j,k)   
                eTzz(i,j,k)  =  tmunuJtoE * eTzz(i,j,k)   
        
        else                   
        
                mTtt(i,j,k)  =  tmunuJtoE * mTtt(i,j,k)   
                mTtx(i,j,k)  =  tmunuJtoE * mTtx(i,j,k)   
                mTty(i,j,k)  =  tmunuJtoE * mTty(i,j,k)   
                mTtz(i,j,k)  =  tmunuJtoE * mTtz(i,j,k)   
                mTxx(i,j,k)  =  tmunuJtoE * mTxx(i,j,k)   
                mTxy(i,j,k)  =  tmunuJtoE * mTxy(i,j,k)   
                mTxz(i,j,k)  =  tmunuJtoE * mTxz(i,j,k)   
                mTyy(i,j,k)  =  tmunuJtoE * mTyy(i,j,k)   
                mTyz(i,j,k)  =  tmunuJtoE * mTyz(i,j,k)   
                mTzz(i,j,k)  =  tmunuJtoE * mTzz(i,j,k)   
        end if         

       end do
     end do
   end do
    


     !$OMP END PARALLEL DO

   if (EJverbose.eq.1) then 
        call   CCTK_INFO("rho_center in J frame")
        write(message, *) rho(0,0,0), AJtoE(0,0,0) 
        call   CCTK_INFO(message)
    end if

end subroutine JordanToEinstein


! Einstein to Jordan Frame
!
!  The quantities transformed, along with their transformation factors are:
!
!	g_ij^J   =  A^2(phi)  g_ij^E
!
!	alp^J    =  A^1(phi) alp^E  
!       
!       rho^J    =  A^-4(phi) rho^E
!
!       beta^iJ  =  beta^iE
!
!	eps^J    =  A^-4(phi) eps^E
!
!	p^J      =  A^-4(phi) p^E
!
!	u_mu^J   =  A^1(phi) u_mu^E
!
!       T_munu^J =  A^-2(phi) T_munu^E



subroutine  EinsteinToJordan( CCTK_ARGUMENTS )
   implicit none
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
      
   !Variables
   CCTK_REAL                A_conf, phiJ
   CCTK_REAL                metricEtoJ, alpEtoJ, hydroEtoJ, velEtoJ, tmunuEtoJ
   CCTK_INT                 i, j, k
   CCTK_INT                 a, b, c

   character(len=200) :: message
! Allocate memory for the character array
! Format and store values in the character array
   if (EJverbose.eq.1) then 
           write(message, *) rho(0,0,0)
           call CCTK_INFO("rho_center in E frame")
           call CCTK_INFO(message)
           call CCTK_INFO("Going now to Jordan frame")
   end if 

   !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(k,j,i, a,b,c, A_conf,metricEtoJ,alpEtoJ,hydroEtoJ,velEtoJ,tmunuEtoJ)
       
   do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
     do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
       do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

        !---Get local variables--------
        phiJ = phi1(i,j,k)  
       
        !---Calculate conformal factor A(phi) 

        if (decoupling_limit.eq.1) then 
                A_conf = 1.0
        else        
                A_conf = EXP(0.5 * betaDEF * phiJ * phiJ   +  k0BD * phiJ)
        end if 

        AEtoJ(i,j,k) = A_conf

        !---Conversion factors
        metricEtoJ = (A_conf*A_conf)
        alpEtoJ      = A_conf
        tmunuEtoJ   = 1/(A_conf*A_conf)

        !----Assign values to grid function
        gxx(i,j,k)   =  metricEtoJ * gxx(i,j,k)   
        gxy(i,j,k)   =  metricEtoJ * gxy(i,j,k) 
        gxz(i,j,k)   =  metricEtoJ * gxz(i,j,k) 
        gyy(i,j,k)   =  metricEtoJ * gyy(i,j,k) 
        gyz(i,j,k)   =  metricEtoJ * gyz(i,j,k) 
        gzz(i,j,k)   =  metricEtoJ * gzz(i,j,k) 
        alp(i,j,k)   =  alpEtoJ * alp(i,j,k)

        
        if (decoupling_limit.eq.1) then

                eTtt(i,j,k)  =  tmunuEtoJ * eTtt(i,j,k)   
                eTtx(i,j,k)  =  tmunuEtoJ * eTtx(i,j,k)   
                eTty(i,j,k)  =  tmunuEtoJ * eTty(i,j,k)   
                eTtz(i,j,k)  =  tmunuEtoJ * eTtz(i,j,k)   
                eTxx(i,j,k)  =  tmunuEtoJ * eTxx(i,j,k)   
                eTxy(i,j,k)  =  tmunuEtoJ * eTxy(i,j,k)   
                eTxz(i,j,k)  =  tmunuEtoJ * eTxz(i,j,k)   
                eTyy(i,j,k)  =  tmunuEtoJ * eTyy(i,j,k)   
                eTyz(i,j,k)  =  tmunuEtoJ * eTyz(i,j,k)   
                eTzz(i,j,k)  =  tmunuEtoJ * eTzz(i,j,k)   
        
        else                   
        
                mTtt(i,j,k)  =  tmunuEtoJ * mTtt(i,j,k)   
                mTtx(i,j,k)  =  tmunuEtoJ * mTtx(i,j,k)   
                mTty(i,j,k)  =  tmunuEtoJ * mTty(i,j,k)   
                mTtz(i,j,k)  =  tmunuEtoJ * mTtz(i,j,k)   
                mTxx(i,j,k)  =  tmunuEtoJ * mTxx(i,j,k)   
                mTxy(i,j,k)  =  tmunuEtoJ * mTxy(i,j,k)   
                mTxz(i,j,k)  =  tmunuEtoJ * mTxz(i,j,k)   
                mTyy(i,j,k)  =  tmunuEtoJ * mTyy(i,j,k)   
                mTyz(i,j,k)  =  tmunuEtoJ * mTyz(i,j,k)   
                mTzz(i,j,k)  =  tmunuEtoJ * mTzz(i,j,k)   
        end if         

       end do
     end do
   end do
     !$OMP END PARALLEL DO

   if (EJverbose.eq.1) then 
        call CCTK_INFO("rho_center in J frame")   
        write(message, *) rho(0,0,0), AEtoJ(0,0,0) 
        call CCTK_INFO(message)
    end if 
end subroutine EinsteinToJordan 
