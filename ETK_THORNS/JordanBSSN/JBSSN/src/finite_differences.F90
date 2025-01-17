! fin_dif_deriv.F90
!
! Contains  calc_d1(array, i,j,k, order, dx,dy,dz)  returns d1_array(3)
!           calc_d2(array, i,j,k, order, dx, dy, dz, k_sum) returns d2_array(3,3)
!           kahan_sum(summ, c, input) returns summ
!           calc_ad1(array, i,j,k, order, dx,dy,dz, di,dj,dk)


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   1st derivative of  grid  function f(x,y,z)       
!     inputs:
!       array:    grid function  f(x,y,z)
!       nx,ny,nz: array size (can vary depending on MPI processes)
!       i,j,k:    grid indexing
!       order:    order of the scheme
!       dx,dy,dz: grid spacing
!      output:
!         d1_array(3): derivative of array 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_d1(array,nx,ny,nz, i, j, k, order, dx,dy,dz, d1_array)
    implicit none
    
    CCTK_INT, intent(in) :: nx, ny, nz
    CCTK_REAL, dimension(nx,ny,ny), intent(in) :: array
    CCTK_INT, intent(in) :: i, j, k, order
    CCTK_REAL  dx, dy, dz
    CCTK_REAL, intent(inout) ::  d1_array(3)
          
    select case(order)
    case(2)
      d1_array(1) = (array(i+1,j,k) - array(i-1,j,k)) / (2.0*dx)
      d1_array(2) = (array(i,j+1,k) - array(i,j-1,k)) / (2.0*dy)
      d1_array(3) = (array(i,j,k+1) - array(i,j,k-1)) / (2.0*dz)
    case(4)
      d1_array(1) = (-array(i+2,j,k) + 8.0*array(i+1,j,k) - 8.0*array(i-1,j,k) + array(i-2,j,k)) / (12.0*dx)
      d1_array(2) = (-array(i,j+2,k) + 8.0*array(i,j+1,k) - 8.0*array(i,j-1,k) + array(i,j-2,k)) / (12.0*dy)
      d1_array(3) = (-array(i,j,k+2) + 8.0*array(i,j,k+1) - 8.0*array(i,j,k-1) + array(i,j,k-2)) / (12.0*dz)
    case(6)
      d1_array(1) = (-array(i+3,j,k) + 9.0*array(i+2,j,k) - 45.0*array(i+1,j,k) + 45.0*array(i-1,j,k) - 9.0*array(i-2,j,k) + array(i-3,j,k)) / (60.0*dx)
      d1_array(2) = (-array(i,j+3,k) + 9.0*array(i,j+2,k) - 45.0*array(i,j+1,k) + 45.0*array(i,j-1,k) - 9.0*array(i,j-2,k) + array(i,j-3,k)) / (60.0*dy)
      d1_array(3) = (-array(i,j,k+3) + 9.0*array(i,j,k+2) - 45.0*array(i,j,k+1) + 45.0*array(i,j,k-1) - 9.0*array(i,j,k-2) + array(i,j,k-3)) / (60.0*dz)
    case default
      write(*,*) "Error: Unsupported derivative order:", order    
    end select
  end subroutine calc_d1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   2nd derivative of grid function f(x,y,z)       
!     inputs:
!       array:    grid function  f(x,y,z)
!       nx,ny,nz: array size (can vary depending on MPI processes)
!       i,j,k:    grid indexing
!       order:    order of the scheme
!       dx,dy,dz: grid spacing
!       k_sum:    Add Kahan sum correction
!      output:
!         d2_array(3,3): derivative of array 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_d2(array,nx,ny,nz, i, j, k, order, dx, dy, dz, k_sum, d2_array)
implicit none

    CCTK_INT, intent(in) :: nx, ny, nz
    CCTK_REAL, dimension(nx,ny,nz), intent(in) :: array
    CCTK_INT, intent(in) :: i, j, k, order
    logical, intent(in) :: k_sum
    CCTK_REAL, intent(in) :: dx, dy, dz
    CCTK_REAL, intent(inout) :: d2_array(3,3)    
    CCTK_REAL  dx2, dy2, dz2, dxdy, dxdz, dydz
    CCTK_REAL  temp1, temp2, temp3, temp4
    CCTK_REAL  summ, c, y, t ! Variables for Kahan summation


    d2_array = 0.0d0
    dx2 = dx * dx
    dy2 = dy * dy
    dz2 = dz * dz
    dxdy = dx * dy
    dxdz = dx * dz
    dydz = dy * dz

    select case(order)
    case(2)
      ! Second-order accurate 2nd derivatives using combination of first derivatives
      temp1 = (array(i+1,j,k) - array(i,j,k)) / dx
      temp2 = (array(i,j,k) - array(i-1,j,k)) / dx
      d2_array(1,1) = (temp1 - temp2) / dx

      temp1 = (array(i,j+1,k) - array(i,j,k)) / dy
      temp2 = (array(i,j,k) - array(i,j-1,k)) / dy
      d2_array(2,2) = (temp1 - temp2) / dy

      temp1 = (array(i,j,k+1) - array(i,j,k)) / dz
      temp2 = (array(i,j,k) - array(i,j,k-1)) / dz
      d2_array(3,3) = (temp1 - temp2) / dz

      ! Mixed derivatives (2nd order) 
      temp1 = (array(i+1,j+1,k) - array(i+1,j-1,k)) / (2.0d0*dy)
      temp2 = (array(i-1,j+1,k) - array(i-1,j-1,k)) / (2.0d0*dy)
      d2_array(1,2) = (temp1 - temp2) / (2.0d0*dx)

      temp1 = (array(i+1,j,k+1) - array(i+1,j,k-1)) / (2.0d0*dz)
      temp2 = (array(i-1,j,k+1) - array(i-1,j,k-1)) / (2.0d0*dz)
      d2_array(1,3) = (temp1 - temp2) / (2.0d0*dx)

      temp1 = (array(i,j+1,k+1) - array(i,j+1,k-1)) / (2.0d0*dz)
      temp2 = (array(i,j-1,k+1) - array(i,j-1,k-1)) / (2.0d0*dz)
      d2_array(2,3) = (temp1 - temp2) / (2.0d0*dy)

    case(4) ! Fourth order second derivatives
            if (k_sum) then  ! with summation compensation

                summ = 0.0d0
                c = 0.0d0
                call kahan_sum(summ, c, -array(i+2,j,k))
                call kahan_sum(summ, c, 16.0d0*array(i+1,j,k))
                call kahan_sum(summ, c, -30.0d0*array(i,j,k))
                call kahan_sum(summ, c, 16.0d0*array(i-1,j,k))
                call kahan_sum(summ, c, -array(i-2,j,k))
                d2_array(1,1) = summ / (12.0d0*dx2)

                summ = 0.0d0
                c = 0.0d0
                call kahan_sum(summ, c, -array(i,j+2,k))
                call kahan_sum(summ, c, 16.0d0*array(i,j+1,k))
                call kahan_sum(summ, c, -30.0d0*array(i,j,k))
                call kahan_sum(summ, c, 16.0d0*array(i,j-1,k))
                call kahan_sum(summ, c, -array(i,j-2,k))
                d2_array(2,2) = summ / (12.0d0*dy2)

                summ = 0.0d0
                c = 0.0d0
                call kahan_sum(summ, c, -array(i,j,k+2))
                call kahan_sum(summ, c, 16.0d0*array(i,j,k+1))
                call kahan_sum(summ, c, -30.0d0*array(i,j,k))
                call kahan_sum(summ, c, 16.0d0*array(i,j,k-1))
                call kahan_sum(summ, c, -array(i,j,k-2))
                d2_array(3,3) = summ / (12.0d0*dz2)

                        ! Mixed derivatives with Kahan summation compensation
                summ = 0.0d0
                c = 0.0d0
                call kahan_sum(summ, c, array(i+2,j+2,k))
                call kahan_sum(summ, c, -8.0d0*array(i+2,j+1,k))
                call kahan_sum(summ, c, 8.0d0*array(i+2,j-1,k))
                call kahan_sum(summ, c, -array(i+2,j-2,k))
                call kahan_sum(summ, c, -8.0d0*array(i+1,j+2,k))
                call kahan_sum(summ, c, 64.0d0*array(i+1,j+1,k))
                call kahan_sum(summ, c, -64.0d0*array(i+1,j-1,k))
                call kahan_sum(summ, c, 8.0d0*array(i+1,j-2,k))
                call kahan_sum(summ, c, 8.0d0*array(i-1,j+2,k))
                call kahan_sum(summ, c, -64.0d0*array(i-1,j+1,k))
                call kahan_sum(summ, c, 64.0d0*array(i-1,j-1,k))
                call kahan_sum(summ, c, -8.0d0*array(i-1,j-2,k))
                call kahan_sum(summ, c, -array(i-2,j+2,k))
                call kahan_sum(summ, c, 8.0d0*array(i-2,j+1,k))
                call kahan_sum(summ, c, -8.0d0*array(i-2,j-1,k))
                call kahan_sum(summ, c, array(i-2,j-2,k))
                d2_array(1,2) = summ / (144.0d0*dxdy)

                summ = 0.0d0
                c = 0.0d0
                call kahan_sum(summ, c, array(i+2,j,k+2))
                call kahan_sum(summ, c, -8.0d0*array(i+2,j,k+1))
                call kahan_sum(summ, c, 8.0d0*array(i+2,j,k-1))
                call kahan_sum(summ, c, -array(i+2,j,k-2))
                call kahan_sum(summ, c, -8.0d0*array(i+1,j,k+2))
                call kahan_sum(summ, c, 64.0d0*array(i+1,j,k+1))
                call kahan_sum(summ, c, -64.0d0*array(i+1,j,k-1))
                call kahan_sum(summ, c, 8.0d0*array(i+1,j,k-2))
                call kahan_sum(summ, c, 8.0d0*array(i-1,j,k+2))
                call kahan_sum(summ, c, -64.0d0*array(i-1,j,k+1))
                call kahan_sum(summ, c, 64.0d0*array(i-1,j,k-1))
                call kahan_sum(summ, c, -8.0d0*array(i-1,j,k-2))
                call kahan_sum(summ, c, -array(i-2,j,k+2))
                call kahan_sum(summ, c, 8.0d0*array(i-2,j,k+1))
                call kahan_sum(summ, c, -8.0d0*array(i-2,j,k-1))
                call kahan_sum(summ, c, array(i-2,j,k-2))
                d2_array(1,3) = summ / (144.0d0*dxdz)

                summ = 0.0d0
                c = 0.0d0
                call kahan_sum(summ, c, array(i,j+2,k+2))
                call kahan_sum(summ, c, -8.0d0*array(i,j+2,k+1))
                call kahan_sum(summ, c, 8.0d0*array(i,j+2,k-1))
                call kahan_sum(summ, c, -array(i,j+2,k-2))
                call kahan_sum(summ, c, -8.0d0*array(i,j+1,k+2))
                call kahan_sum(summ, c, 64.0d0*array(i,j+1,k+1))
                call kahan_sum(summ, c, -64.0d0*array(i,j+1,k-1))
                call kahan_sum(summ, c, 8.0d0*array(i,j+1,k-2))
                call kahan_sum(summ, c, 8.0d0*array(i,j-1,k+2))
                call kahan_sum(summ, c, -64.0d0*array(i,j-1,k+1))
                call kahan_sum(summ, c, 64.0d0*array(i,j-1,k-1))
                call kahan_sum(summ, c, -8.0d0*array(i,j-1,k-2))
                call kahan_sum(summ, c, -array(i,j-2,k+2))
                call kahan_sum(summ, c, 8.0d0*array(i,j-2,k+1))
                call kahan_sum(summ, c, -8.0d0*array(i,j-2,k-1))
                call kahan_sum(summ, c, array(i,j-2,k-2))
                d2_array(2,3) = summ / (144.0d0*dydz)

             else
               ! Fourth-order accurate second derivatives without summation compensation
                temp1 = (-array(i+2,j,k) + 16.0d0*array(i+1,j,k) - 30.0d0*array(i,j,k) + 16.0d0*array(i-1,j,k) - array(i-2,j,k))
                d2_array(1,1) = temp1 / (12.0d0*dx2)
             
                temp1 = (-array(i,j+2,k) + 16.0d0*array(i,j+1,k) - 30.0d0*array(i,j,k) + 16.0d0*array(i,j-1,k) - array(i,j-2,k))
                d2_array(2,2) = temp1 / (12.0d0*dy2)
                
                temp1 = (-array(i,j,k+2) + 16.0d0*array(i,j,k+1) - 30.0d0*array(i,j,k) + 16.0d0*array(i,j,k-1) - array(i,j,k-2))
                d2_array(3,3) = temp1 / (12.0d0*dz2)
                
                   ! Mixed derivatives (4th order accurate)
                temp1 = array(i+2,j+2,k) - 8.0d0*array(i+2,j+1,k) + 8.0d0*array(i+2,j-1,k) - array(i+2,j-2,k)
                temp2 = -8.0d0*array(i+1,j+2,k) + 64.0d0*array(i+1,j+1,k) - 64.0d0*array(i+1,j-1,k) + 8.0d0*array(i+1,j-2,k)
                temp3 = 8.0d0*array(i-1,j+2,k) - 64.0d0*array(i-1,j+1,k) + 64.0d0*array(i-1,j-1,k) - 8.0d0*array(i-1,j-2,k)
                temp4 = -array(i-2,j+2,k) + 8.0d0*array(i-2,j+1,k) - 8.0d0*array(i-2,j-1,k) + array(i-2,j-2,k)
                d2_array(1,2) = (temp1 + temp2 + temp3 + temp4) / (144.0d0*dxdy)
                
                temp1 = array(i+2,j,k+2) - 8.0d0*array(i+2,j,k+1) + 8.0d0*array(i+2,j,k-1) - array(i+2,j,k-2)
                temp2 = -8.0d0*array(i+1,j,k+2) + 64.0d0*array(i+1,j,k+1) - 64.0d0*array(i+1,j,k-1) + 8.0d0*array(i+1,j,k-2)
                temp3 = 8.0d0*array(i-1,j,k+2) - 64.0d0*array(i-1,j,k+1) + 64.0d0*array(i-1,j,k-1) - 8.0d0*array(i-1,j,k-2)
                temp4 = -array(i-2,j,k+2) + 8.0d0*array(i-2,j,k+1) - 8.0d0*array(i-2,j,k-1) + array(i-2,j,k-2)
                d2_array(1,3) = (temp1 + temp2 + temp3 + temp4) / (144.0d0*dxdz)
                
                temp1 = array(i,j+2,k+2) - 8.0d0*array(i,j+2,k+1) + 8.0d0*array(i,j+2,k-1) - array(i,j+2,k-2)
                temp2 = -8.0d0*array(i,j+1,k+2) + 64.0d0*array(i,j+1,k+1) - 64.0d0*array(i,j+1,k-1) + 8.0d0*array(i,j+1,k-2)
                temp3 = 8.0d0*array(i,j-1,k+2) - 64.0d0*array(i,j-1,k+1) + 64.0d0*array(i,j-1,k-1) - 8.0d0*array(i,j-1,k-2)
                temp4 = -array(i,j-2,k+2) + 8.0d0*array(i,j-2,k+1) - 8.0d0*array(i,j-2,k-1) + array(i,j-2,k-2)
                d2_array(2,3) = (temp1 + temp2 + temp3 + temp4) / (144.0d0*dydz)

         end if 
 
    case default
      write(*,*) "Error: Unsupported derivative order:", order      
    end select

    ! Fill in the symmetric part of the tensor
    d2_array(2,1) = d2_array(1,2)
    d2_array(3,1) = d2_array(1,3)
    d2_array(3,2) = d2_array(2,3)
    
  end subroutine calc_d2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Kahan summation subroutine
!       applies Kahan summation compensation to improve accuracy
!       useful, but more expensive 
!   input: 
!       sum: Accumulated value of sum
!       c: compensation
!       input: term to add        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine kahan_sum(sum, c, input)
    CCTK_REAL, intent(inout) :: sum, c
    CCTK_REAL, intent(in) :: input
    CCTK_REAL  y, t
    y = input - c
        !write(*,*) "y =", y
    t = sum + y
        !write(*,*) "t =", t
    c = (t - sum) - y
        !write(*,*) "c =", c
    sum = t
        !write(*,*) "sum =", sum
  end subroutine kahan_sum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Advective derivatives
!        5th order method to calculate advective 1st order derivatives
!       input: 
!         array:     grid function 
!         beta:      shift in local coordinates
!         nx,ny,nz:  array size (can vary depending on MPI processes)
!         i,j,k:     grid index 
!         dx,dy,dz:  grid spacing
!         di,dj,dk:  direction of advection
!       output: 
!         ad1_f: advective derivative of grid function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_ad1(array,beta_s,nx,ny,nz, i,j,k, dx,dy,dz, di,dj,dk,ad1_f)

    CCTK_INT, intent(in) :: nx, ny, nz
    CCTK_REAL, dimension(nx,ny,nz), intent(in) :: array    
    CCTK_REAL, intent(in) :: beta_s(3)
    CCTK_INT, intent(in) :: i, j, k, di, dj, dk
    CCTK_REAL, intent(in) :: dx, dy, dz
    CCTK_REAL   d1_f(3)
    CCTK_REAL, intent(inout) :: ad1_f

           d1_f(1) = di * ( -3.0d0*array(i-di,j,k) - 10.0d0*array(i,j,k) + 18.0d0*array(i+di,j,k)   &
                - 6.0d0*array(i+2*di,j,k) + array(i+3*di,j,k)) / (12.0d0*dx)
           d1_f(2) = dj * ( -3.0d0*array(i,j-dj,k) - 10.0d0*array(i,j,k) + 18.0d0*array(i,j+dj,k)   &
                - 6.0d0*array(i,j+2*dj,k) + array(i,j+3*dj,k)) / (12.0d0*dy)
           d1_f(3) = dk * ( -3.0d0*array(i,j,k-dk) - 10.0d0*array(i,j,k) + 18.0d0*array(i,j,k+dk)   &
                - 6.0d0*array(i,j,k+2*dk) + array(i,j,k+3*dk)) / (12.0d0*dz)
           ad1_f = beta_s(1)*d1_f(1) + beta_s(2)*d1_f(2) + beta_s(3)*d1_f(3)
  end subroutine calc_ad1



