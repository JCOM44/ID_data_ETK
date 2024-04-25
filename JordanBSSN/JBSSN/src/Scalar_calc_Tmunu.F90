#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine JordanFBSSN_calc_scalarTmunu( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3), betad(3)
  CCTK_REAL                gg(4,4), gu(4,4), deth
  CCTK_REAL                lphi, lKphi
  CCTK_REAL                Tab(4,4)


  ! First derivatives
  CCTK_REAL                d1_lphi(4)


  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12
  CCTK_REAL                aux
  CCTK_REAL, parameter ::  pi   = acos(-1.0)
  CCTK_INT                 i, j, k
  CCTK_INT                 a, b, c

  ! jacobian
  integer                  istat
  logical                  use_jacobian
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ11, lJ12, lJ13
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ21, lJ22, lJ23
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ31, lJ32, lJ33
  CCTK_POINTER             lJ11_ptr, lJ12_ptr, lJ13_ptr
  CCTK_POINTER             lJ21_ptr, lJ22_ptr, lJ23_ptr
  CCTK_POINTER             lJ31_ptr, lJ32_ptr, lJ33_ptr
  CCTK_REAL                jac(3,3)

  pointer (lJ11_ptr, lJ11), (lJ12_ptr, lJ12), (lJ13_ptr, lJ13)
  pointer (lJ21_ptr, lJ21), (lJ22_ptr, lJ22), (lJ23_ptr, lJ23)
  pointer (lJ31_ptr, lJ31), (lJ32_ptr, lJ32), (lJ33_ptr, lJ33)

  call CCTK_IsFunctionAliased(istat, "MultiPatch_GetDomainSpecification")
  if (istat == 0) then
     use_jacobian = .false.
  else
     use_jacobian = .true.
  end if

  if (use_jacobian) then
     call CCTK_VarDataPtr(lJ11_ptr, cctkGH, 0, "Coordinates::J11")
     call CCTK_VarDataPtr(lJ12_ptr, cctkGH, 0, "Coordinates::J12")
     call CCTK_VarDataPtr(lJ13_ptr, cctkGH, 0, "Coordinates::J13")
     call CCTK_VarDataPtr(lJ21_ptr, cctkGH, 0, "Coordinates::J21")
     call CCTK_VarDataPtr(lJ22_ptr, cctkGH, 0, "Coordinates::J22")
     call CCTK_VarDataPtr(lJ23_ptr, cctkGH, 0, "Coordinates::J23")
     call CCTK_VarDataPtr(lJ31_ptr, cctkGH, 0, "Coordinates::J31")
     call CCTK_VarDataPtr(lJ32_ptr, cctkGH, 0, "Coordinates::J32")
     call CCTK_VarDataPtr(lJ33_ptr, cctkGH, 0, "Coordinates::J33")
  end if

  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)


  !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(k,j,i, a,b, aux,&
  !$OMP                                 alph,beta,betad,&
  !$OMP                                 gg,gu,deth,&
  !$OMP                                 lphi,lKphi,&
  !$OMP                                 Tab, d1_lphi)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
     do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
        do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

           !------------ Get local variables ----------
           lphi     = phi1(i,j,k)

           lKphi    = Kphi1(i,j,k)


           alph      = alp(i,j,k)

           beta(1)   = betax(i,j,k)
           beta(2)   = betay(i,j,k)
           beta(3)   = betaz(i,j,k)

           gg(1,1)   = gxx(i,j,k)
           gg(1,2)   = gxy(i,j,k)
           gg(1,3)   = gxz(i,j,k)
           gg(2,2)   = gyy(i,j,k)
           gg(2,3)   = gyz(i,j,k)
           gg(3,3)   = gzz(i,j,k)
           gg(2,1)   = gg(1,2)
           gg(3,1)   = gg(1,3)
           gg(3,2)   = gg(2,3)

           if (use_jacobian) then
              jac(1,1) = lJ11(i,j,k)
              jac(1,2) = lJ12(i,j,k)
              jac(1,3) = lJ13(i,j,k)
              jac(2,1) = lJ21(i,j,k)
              jac(2,2) = lJ22(i,j,k)
              jac(2,3) = lJ23(i,j,k)
              jac(3,1) = lJ31(i,j,k)
              jac(3,2) = lJ32(i,j,k)
              jac(3,3) = lJ33(i,j,k)
           else
              jac      = 0.0d0
              jac(1,1) = 1.0d0
              jac(2,2) = 1.0d0
              jac(3,3) = 1.0d0
           end if

           ! now we compute beta_i (betad)
           betad = 0.0d0
           do a = 1, 3
              do b = 1, 3
                 betad(a) = betad(a) + gg(a,b) * beta(b)
              end do
           end do

           ! and finish with the rest of the 4-metric gg.
           ! we will use the 4th slot for time throughout
           gg(1,4)   = betad(1)
           gg(2,4)   = betad(2)
           gg(3,4)   = betad(3)
           gg(4,1)   = gg(1,4)
           gg(4,2)   = gg(2,4)
           gg(4,3)   = gg(3,4)

           gg(4,4)   = -alph * alph
           do a = 1, 3
              gg(4,4) = gg(4,4) + beta(a) * betad(a)
           end do

           !------------ Invert metric ----------------

           ! determinant of the 3-metric
           deth    =     gg(1,1) * gg(2,2) * gg(3,3)                              &
                   + 2 * gg(1,2) * gg(1,3) * gg(2,3)                              &
                   -     gg(1,1) * gg(2,3) ** 2                                   &
                   -     gg(2,2) * gg(1,3) ** 2                                   &
                   -     gg(3,3) * gg(1,2) ** 2

           gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / deth               &
                   - beta(1) * beta(1) / (alph * alph)
           gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / deth               &
                   - beta(2) * beta(2) / (alph * alph)
           gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / deth               &
                   - beta(3) * beta(3) / (alph * alph)
           gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / deth               &
                   - beta(1) * beta(2) / (alph * alph)
           gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / deth               &
                   - beta(1) * beta(3) / (alph * alph)
           gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / deth               &
                   - beta(2) * beta(3) / (alph * alph)

           gu(1,4) = beta(1) / (alph * alph)
           gu(2,4) = beta(2) / (alph * alph)
           gu(3,4) = beta(3) / (alph * alph)
           gu(4,4) = -1.0d0 / (alph * alph)

           gu(2,1) = gu(1,2)
           gu(3,1) = gu(1,3)
           gu(3,2) = gu(2,3)

           gu(4,1) = gu(1,4)
           gu(4,2) = gu(2,4)
           gu(4,3) = gu(3,4)


           !------------- Centered 1st derivatives -----------

           ! d1_lphi1(3)
           d1_lphi(1)  = (   -phi1(i+2,j,k) + 8*phi1(i+1,j,k)                        &
                           - 8*phi1(i-1,j,k) +   phi1(i-2,j,k) ) / dx12

           d1_lphi(2)  = (   -phi1(i,j+2,k) + 8*phi1(i,j+1,k)                        &
                           - 8*phi1(i,j-1,k) +   phi1(i,j-2,k) ) / dy12

           d1_lphi(3)  = (   -phi1(i,j,k+2) + 8*phi1(i,j,k+1)                        &
                           - 8*phi1(i,j,k-1) +   phi1(i,j,k-2) ) / dz12

           !-------------------------------------------
           if (use_jacobian) then
              call JordanFBSSN_apply_jacobian(d1_lphi, jac)
           end if
           !-------------------------------------------

           ! time derivatives
           d1_lphi(4)  = -2 * alph * lKphi

           do a = 1, 3
              d1_lphi(4) = d1_lphi(4) + beta(a) * d1_lphi(a)
           end do

           !-------------------------------------------

           aux = 0.0                                            

           do a = 1, 4
              do b = 1, 4
                 aux = aux + gu(a,b) *  d1_lphi(a) * d1_lphi(b)
              end do
           end do

           ! compute the stress-energy tensor
           do a = 1, 4
              do b = 1, 4
                 Tab(a,b) =    ( 2 * d1_lphi(a) * d1_lphi(b)   - gg(a,b) * aux )/ (8*pi)
              end do
           end do

           ! and finally store it in the Tmunu variables
           sfTtt(i,j,k) =  Tab(4,4) 
           sfTtx(i,j,k) =  Tab(4,1) 
           sfTty(i,j,k) =  Tab(4,2) 
           sfTtz(i,j,k) =  Tab(4,3) 
           sfTxx(i,j,k) =  Tab(1,1) 
           sfTxy(i,j,k) =  Tab(1,2) 
           sfTxz(i,j,k) =  Tab(1,3) 
           sfTyy(i,j,k) =  Tab(2,2) 
           sfTyz(i,j,k) =  Tab(2,3) 
           sfTzz(i,j,k) =  Tab(3,3) 



        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine JordanFBSSN_calc_scalarTmunu
