! fin_dif_deriv.F90
!
! Contains  d1_scalar(array, i,j,k, order, dx)  returns d1_array(3)
!           d1_vector(array, i,j,k, order, dx)  returns d1_array(3)
!           d1_tensor(array, i,j,k, order, dx)      
!           d2_scalar(array, i,j,k, order, dx)
!           d2_vector(array, i,j,k, order, dx)
!           d2_tensor(array, i,j,k, order, dx)
!          ad1_scalar(array, i,j,k, order, dx)
!          ad1_vector(array, i,j,k, order, dx)
!          ad1_tensor(array, i,j,k, order, dx)      

      

module finite_difference_mod
  implicit none
contains
  function d1_scalar(array, i, j, k, order, dx) result(d1_array)
    implicit none
    real, dimension(:,:,:), intent(in) :: array
    integer, intent(in) :: i, j, k, order
    real, intent(in) :: dx
    real :: d1_array(3)
    
    if (order==6) then
    d1_array(1) = ( - array(i+2, j, k) + 16*array(i+1, j, k) - 30*array(i, j, k) &
             + 16*array(i-1, j, k) - array(i-2, j, k) ) / dxsq12
     

     else if (order==6) 
        
     d1_array(1) = ( array(i+3,j,k) - 9*array(i+2,j,k) + 45*array(i+1,j,k)  &
                    -array(i-3,j,k) + 9*array(i-2,j,k) - 45*array(i-1,j,k) )

     d1_array(2) = ( array(i,j+3,k) - 9*array(i,j+2,k) + 45*array(i,j,k+2)  &
                    -array(i,j-3,k) + 9*array(i,j-2,k) - 45*array(i,j,k-2) )
     
     d1_array(3) = ( array(i,j,k+3) - 9*array(i,j,k+2) + 45*array(i,j,k+1)  &
                    -array(i,j,k-3) + 9*array(i,j,k-2) - 45*array(i,j,k-1) )
     else 
        write(*,*) "---------ERROR in derivaritve order----"
     end if 

  end function compute_d2_hh

  function compute_d1_hh(array, i, j, k, dx) result(d1_hh)
    implicit none
    real, dimension(:,:,:), intent(in) :: array
    integer, intent(in) :: i, j, k
    real, intent(in) :: dx
    real :: d1_hh

    d1_hh = (array(i+1, j, k) - array(i-1, j, k)) / (2.0 * dx)
  end function compute_d1_hh
end module finite_difference_mod

