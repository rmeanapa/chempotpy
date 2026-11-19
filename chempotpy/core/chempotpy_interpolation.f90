!> @brief Provide reusable interpolation algorithms used by PES kernels.
module chempotpy_interpolation
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: cubic_spline_coefficients
    public :: evaluate_cubic_spline
    public :: evaluate_tabulated_cubic_spline
    public :: prepare_tabulated_cubic_spline

    contains

    !> @brief Build cubic-spline second derivatives with natural or clamped endpoints.
    subroutine cubic_spline_coefficients( coordinates, values, first_derivative, last_derivative, second_derivatives )
        real(wp), intent(in)  :: coordinates(:), values(:), first_derivative, last_derivative
        real(wp), intent(out) :: second_derivatives(size(coordinates))
        integer               :: index
        real(wp)              :: factor, last_factor, last_value, ratio, temporary(size(coordinates))
        if( first_derivative > 0.99e30_wp )then
            second_derivatives(1) = 0.0_wp
            temporary(1)          = 0.0_wp
        else
            second_derivatives(1) = -0.5_wp
            temporary(1)          = (3.0_wp/(coordinates(2)-coordinates(1)))* ((values(2)-values(1))/(coordinates(2)- &
                & coordinates(1))-first_derivative)
        endif
        do index = 2, size(coordinates)-1
            ratio                     = (coordinates(index)-coordinates(index-1))/ (coordinates(index+1)- &
                & coordinates(index-1))
            factor                    = ratio*second_derivatives(index-1)+2.0_wp
            second_derivatives(index) = (ratio-1.0_wp)/factor
            temporary(index)          = (6.0_wp*((values(index+1)-values(index))/ (coordinates(index+1)- &
                & coordinates(index))-(values(index)-values(index-1))/ (coordinates(index)-coordinates(index-1)))/ &
                & (coordinates(index+1)-coordinates(index-1)) - ratio*temporary(index-1))/factor
        enddo
        if( last_derivative > 0.99e30_wp )then
            last_factor = 0.0_wp
            last_value  = 0.0_wp
        else
            last_factor = 0.5_wp
            last_value  = (3.0_wp/(coordinates(size(coordinates))-coordinates(size(coordinates)-1)))* &
                & (last_derivative-(values(size(values))-values(size(values)-1))/ (coordinates(size(coordinates))- &
                & coordinates(size(coordinates)-1)))
        endif
        second_derivatives(size(coordinates)) = (last_value-last_factor*temporary(size(coordinates)-1))/ (last_factor* &
            & second_derivatives(size(coordinates)-1)+1.0_wp)
        do index = size(coordinates)-1, 1, -1
            second_derivatives(index) = second_derivatives(index)*second_derivatives(index+1)+temporary(index)
        enddo
    end subroutine cubic_spline_coefficients

    !> @brief Evaluate a cubic spline from its tabulated values and second derivatives.
    subroutine evaluate_cubic_spline( coordinates, values, second_derivatives, argument, value )
        real(wp), intent(in)  :: coordinates(:), values(:), second_derivatives(:), argument
        real(wp), intent(out) :: value
        integer               :: index, lower_index, upper_index
        real(wp)              :: interval, lower_weight, upper_weight
        lower_index = 1
        upper_index = size(coordinates)
        do while( upper_index-lower_index > 1 )
            index = (upper_index+lower_index)/2
            if( coordinates(index) > argument )then
                upper_index = index
            else
                lower_index = index
            endif
        enddo
        interval     = coordinates(upper_index)-coordinates(lower_index)
        lower_weight = (coordinates(upper_index)-argument)/interval
        upper_weight = (argument-coordinates(lower_index))/interval
        value        = lower_weight*values(lower_index)+upper_weight*values(upper_index) + ((lower_weight**3- &
            & lower_weight)*second_derivatives(lower_index) + (upper_weight**3-upper_weight)* &
            & second_derivatives(upper_index))*interval**2/6.0_wp
    end subroutine evaluate_cubic_spline

    !> @brief Prepare the legacy tabulated cubic-spline coefficients used by ClH2 surfaces.
    subroutine prepare_tabulated_cubic_spline( coordinates, values, initial_tolerance, second_derivatives, &
            & third_derivatives, first_derivatives )
        real(wp), intent(in)    :: coordinates(:), values(:), initial_tolerance
        real(wp), intent(inout) :: second_derivatives(:)
        real(wp), intent(out)   :: third_derivatives(:), first_derivatives(:)
        integer                 :: index, iteration, interval_count
        real(wp)                :: correction, curvature_change, expansion, half, omega, tolerance
        real(wp)                :: left_weight(size(coordinates)), interval(size(coordinates))
        real(wp)                :: interval_sum(size(coordinates)), target_curvature(size(coordinates))
        expansion      = 1.1_wp
        omega          = 1.0717968_wp
        half           = 0.5_wp
        tolerance      = initial_tolerance
        interval_count = size(coordinates)-1
        iteration      = 0
        do index = 1, interval_count
            interval(index)          = coordinates(index+1)-coordinates(index)
            first_derivatives(index) = (values(index+1)-values(index))/interval(index)
        enddo
        do index = 2, interval_count
            interval_sum(index)       = interval(index-1)+interval(index)
            left_weight(index)        = half*interval(index-1)/interval_sum(index)
            curvature_change          = (first_derivatives(index)-first_derivatives(index-1))/interval_sum(index)
            second_derivatives(index) = curvature_change+curvature_change
            target_curvature(index)   = curvature_change+second_derivatives(index)
        enddo
        do
            correction = 0.0_wp
            do index = 2, interval_count
                curvature_change = (target_curvature(index)-left_weight(index)*second_derivatives(index-1) - (half- &
                    & left_weight(index))*second_derivatives(index+1)-second_derivatives(index))*omega
                correction                = max(correction, abs(curvature_change))
                second_derivatives(index) = second_derivatives(index)+curvature_change
            enddo
            iteration = iteration+1
            if( iteration > 10 ) tolerance=tolerance*expansion
            if( correction < tolerance ) exit
        enddo
        do index = 1, interval_count
            third_derivatives(index) = (second_derivatives(index+1)-second_derivatives(index))/interval(index)
        enddo
    end subroutine prepare_tabulated_cubic_spline

    !> @brief Evaluate a prepared tabulated cubic spline and its first two derivatives.
    subroutine evaluate_tabulated_cubic_spline( coordinates, values, argument, second_derivatives, third_derivatives, &
            & first_derivatives, value, first_derivative, second_derivative )
        real(wp), intent(in)  :: coordinates(:), values(:), argument
        real(wp), intent(in)  :: second_derivatives(:), third_derivatives(:), first_derivatives(:)
        real(wp), intent(out) :: value, first_derivative, second_derivative
        integer               :: index
        real(wp)              :: curvature_sum, left_offset, offset_product, right_offset
        if( argument <= coordinates(1) )then
            value             = values(1)
            first_derivative  = first_derivatives(1)
            second_derivative = second_derivatives(1)
            return
        elseif( argument >= coordinates(size(coordinates)) )then
            value             = values(size(values))
            first_derivative  = first_derivatives(size(coordinates))
            second_derivative = second_derivatives(size(coordinates))
            return
        endif
        index = 1
        do while( argument > coordinates(index) )
            index = index+1
        enddo
        if( argument < coordinates(index) ) index=index-1
        left_offset       = argument-coordinates(index)
        right_offset      = argument-coordinates(index+1)
        offset_product    = left_offset*right_offset
        second_derivative = second_derivatives(index)+left_offset*third_derivatives(index)
        curvature_sum     = (second_derivatives(index)+second_derivatives(index+1)+second_derivative)/6.0_wp
        value             = values(index)+left_offset*first_derivatives(index)+offset_product*curvature_sum
        first_derivative   = first_derivatives(index)+(left_offset+right_offset)*curvature_sum + offset_product* &
            & third_derivatives(index)/6.0_wp
    end subroutine evaluate_tabulated_cubic_spline
end module chempotpy_interpolation
