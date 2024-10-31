module bspline
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer,parameter::bspline_order=3 !3rd order

    contains

    subroutine make_knotsvector(knots,d,npoints,rmax)
        real(dp), intent(out) :: knots(:)
        integer,intent(in)::d,npoints
        integer::n,i
        real(dp),intent(in)::rmax
        real(dp)::dr
        n = npoints + 2*d

        
        dr = rmax/dble(npoints-1)
        do i=1,3
            knots(i) = 0d0
            knots(n-i+1) = rmax
        end do
    
        do i=1,npoints
            knots(i+d) =  dble(i-1)*dr
        end do
    end subroutine

    subroutine bspline_basis_functions(values,x,knots,d)
        implicit none
        real(8),intent(in)::x
        integer,intent(in)::d
        real(8), intent(in) :: knots(:)
        real(8),intent(out) :: values(:)
        integer::i,i0
        values = 0d0

        i0 = find_knot_position(knots,x,size(knots),d)
        do i=i0-3,i0+3
            if (i >= 1 .and. i <= size(values)) then
                values(i) = bspline_basis(d, knots, i, x)
            end if
        end do

        !do i=1,size(values)
        !    values(i) = bspline_basis(d, knots, i, x)
        !end do
    end subroutine

    subroutine bspline_basis_functions_deriv(values,dvalues,x,knots,d)
        implicit none
        real(8),intent(in)::x
        integer,intent(in)::d
        real(8), intent(in) :: knots(:)
        real(8),intent(out) :: values(:),dvalues(:)
        integer::i,i0,n
        values = 0d0
        dvalues = 0d0

        n = size(knots)

        i0 = find_knot_position(knots,x,n,d)

        do i=i0-3,i0+1
            if (i >= 1 .and. i <= size(values)) then
                values(i) = bspline_basis(d, knots, i, x)
                dvalues(i) = bspline_derivative(d, knots, i, x)
            end if
        end do

        !do i=1,size(values)
        !    values(i) = bspline_basis(d, knots, i, x)
        !    dvalues(i) = bspline_derivative(d, knots, i, x)
        !end do
    end subroutine

    ! B-spline
    recursive function bspline_basis(d, knots, i, x) result(basis_value)
      integer, intent(in) :: d, i
      real(dp), intent(in) :: knots(:), x
      real(dp) :: basis_value
      real(dp) :: left, right
  
      if (d == 0) then
          ! 0-th term
          if (knots(i) <= x .and. x < knots(i+1)) then
              basis_value = 1.0_dp
          else
              basis_value = 0.0_dp
          end if
      else
          ! recursive
          if (knots(i+d) - knots(i) /= 0.0_dp) then
              left = (x - knots(i)) / (knots(i+d) - knots(i)) * bspline_basis(d-1, knots, i, x)
          else
              left = 0.0_dp
          end if
          if (knots(i+d+1) - knots(i+1) /= 0.0_dp) then
              right = (knots(i+d+1) - x) / (knots(i+d+1) - knots(i+1)) * bspline_basis(d-1, knots, i+1, x)
          else
              right = 0.0_dp
          end if
          basis_value = left + right
      end if
    end function bspline_basis

    integer function find_knot_position(x,x00,nn,dd) result(low)
        integer,intent(in)::nn,dd
        real(8),intent(in)::x(nn)
        real(8),intent(in)::x00
        integer :: high, mid

        low = 1+dd
        high = nn-dd
        if(x(high) -x00 .eq. 0d0) then
            low = high
            return
        end if


        do while (high - low > 1)
            mid = (low + high) / 2
            if (x(mid) > x00) then
                high = mid
            else
                low = mid
            end if
        end do
    end function


  
    ! B-spline
    recursive function bspline_derivative(d, knots, i, x) result(deriv_value)
      integer, intent(in) :: d, i
      real(dp), intent(in) :: knots(:), x
      real(dp) :: deriv_value
      real(dp) :: left, right
  
      if (d == 0) then
          ! zero for zero-th term
          deriv_value = 0.0_dp
      else
          ! left hand side
          if (knots(i+d) - knots(i) /= 0.0_dp) then
              left = d / (knots(i+d) - knots(i)) * bspline_basis(d-1, knots, i, x)
          else
              left = 0.0_dp
          end if
  
          ! right hand side
          if (knots(i+d+1) - knots(i+1) /= 0.0_dp) then
              right = d / (knots(i+d+1) - knots(i+1)) * bspline_basis(d-1, knots, i+1, x)
          else
              right = 0.0_dp
          end if
  
          deriv_value = left - right
      end if
    end function bspline_derivative
end module