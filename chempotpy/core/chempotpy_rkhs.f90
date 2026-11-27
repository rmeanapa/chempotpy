!> @brief Evaluate immutable reproducing-kernel Hilbert-space fits.
module chempotpy_rkhs
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_ts2_rp25_2d
    public :: evaluate_ts2_rp26_3d

    contains

    !> @brief Evaluate a 2D order-two Taylor-spline/reciprocal-power RKHS.
    pure subroutine evaluate_ts2_rp25_2d( point, angular_grid, radial_grid, coefficients, value )
        real(wp), intent(in)  :: point(2), angular_grid(:), radial_grid(:), coefficients(:)
        real(wp), intent(out) :: value
        real(wp)              :: angular_kernel(size(angular_grid)), correction(1), partial(1)
        real(wp)              :: radial_kernel(size(radial_grid)), term(1)
        integer               :: angular_index, coefficient_index, radial_index
        do angular_index=1, size(angular_grid)
            angular_kernel(angular_index) = taylor_spline_n2(point(1), angular_grid(angular_index))
        enddo
        do radial_index=1, size(radial_grid)
            radial_kernel(radial_index) = reciprocal_power_n2_m5(point(2), radial_grid(radial_index))
        enddo
        partial = 0.0_wp; correction = 0.0_wp; coefficient_index = 0
        do angular_index=1, size(angular_grid)
            do radial_index=1, size(radial_grid)
                coefficient_index = coefficient_index+1
                term(1)           = coefficients(coefficient_index)*angular_kernel(angular_index)* &
                    & radial_kernel(radial_index)
                call add_compensated(term, partial, correction)
            enddo
        enddo
        value = partial(1)+correction(1)
    end subroutine evaluate_ts2_rp25_2d

    !> @brief Evaluate a 3D tensor-product RKHS and its first derivatives.
    pure subroutine evaluate_ts2_rp26_3d( point, angular_grid, first_radial_grid, second_radial_grid, &
            coefficients, value, derivative )
        real(wp), intent(in)  :: point(3), angular_grid(:), first_radial_grid(:), second_radial_grid(:)
        real(wp), intent(in)  :: coefficients(:)
        real(wp), intent(out) :: value, derivative(3)
        real(wp)              :: angular_kernel(size(angular_grid)), angular_derivative(size(angular_grid))
        real(wp)              :: first_kernel(size(first_radial_grid)), first_derivative(size(first_radial_grid))
        real(wp)              :: second_kernel(size(second_radial_grid)), second_derivative(size(second_radial_grid))
        real(wp)              :: correction(4), partial(4), term(4)
        integer               :: coefficient_index, first_index, second_index, third_index
        do first_index = 1, size(angular_grid)
            angular_kernel(first_index)     = taylor_spline_n2(point(1), angular_grid(first_index))
            angular_derivative(first_index) = taylor_spline_n2_derivative(point(1), angular_grid(first_index))
        enddo
        do second_index = 1, size(first_radial_grid)
            first_kernel(second_index) = reciprocal_power_n2_m6(point(2), first_radial_grid(second_index))
            first_derivative(second_index) = reciprocal_power_n2_m6_derivative(point(2), &
                & first_radial_grid(second_index))
        enddo
        do third_index = 1, size(second_radial_grid)
            second_kernel(third_index) = reciprocal_power_n2_m6(point(3), second_radial_grid(third_index))
            second_derivative(third_index) = reciprocal_power_n2_m6_derivative(point(3), &
                & second_radial_grid(third_index))
        enddo
        partial = 0.0_wp; correction = 0.0_wp; coefficient_index = 0
        do first_index = 1, size(angular_grid)
            do second_index = 1, size(first_radial_grid)
                do third_index = 1, size(second_radial_grid)
                    coefficient_index = coefficient_index+1
                    term(1) = coefficients(coefficient_index)*angular_kernel(first_index)* first_kernel(second_index)* &
                        & second_kernel(third_index)
                    term(2) = coefficients(coefficient_index)*angular_derivative(first_index)* &
                        & first_kernel(second_index)*second_kernel(third_index)
                    term(3) = coefficients(coefficient_index)*angular_kernel(first_index)* &
                        & first_derivative(second_index)*second_kernel(third_index)
                    term(4) = coefficients(coefficient_index)*angular_kernel(first_index)* first_kernel(second_index)* &
                        & second_derivative(third_index)
                    call add_compensated(term, partial, correction)
                enddo
            enddo
        enddo
        value      = partial(1)+correction(1)
        derivative = partial(2:4)+correction(2:4)
    end subroutine evaluate_ts2_rp26_3d

    !> @brief Accumulate values with Neumaier compensation for cancellation-sensitive fits.
    pure subroutine add_compensated( values, partial, correction )
        real(wp), intent(in)    :: values(:)
        real(wp), intent(inout) :: partial(size(values)), correction(size(values))
        real(wp)                :: updated
        integer                 :: index
        do index = 1, size(values)
            updated = partial(index)+values(index)
            if( abs(partial(index)) >= abs(values(index)) )then
                correction(index) = correction(index)+(partial(index)-updated)+values(index)
            else
                correction(index) = correction(index)+(values(index)-updated)+partial(index)
            endif
            partial(index) = updated
        enddo
    end subroutine add_compensated

    !> @brief Evaluate the order-two Taylor-spline kernel on a bounded coordinate.
    pure function taylor_spline_n2( first, second ) result(value)
        real(wp), intent(in) :: first, second
        real(wp)             :: value, larger, smaller
        smaller = min(first, second); larger = max(first, second)
        value   = 1.0_wp+smaller*larger+2.0_wp*smaller**2*larger-2.0_wp*smaller**3/3.0_wp
    end function taylor_spline_n2

    !> @brief Differentiate the order-two Taylor-spline kernel by its first coordinate.
    pure function taylor_spline_n2_derivative( first, second ) result(value)
        real(wp), intent(in) :: first, second
        real(wp)             :: value
        if( first <= second )then
            value = second+4.0_wp*second*first-2.0_wp*first**2
        else
            value = second*(2.0_wp*second+1.0_wp)
        endif
    end function taylor_spline_n2_derivative

    !> @brief Evaluate the n=2, m=6 reciprocal-power radial kernel.
    pure function reciprocal_power_n2_m6( first, second ) result(value)
        real(wp), intent(in) :: first, second
        real(wp)             :: value, larger, smaller
        smaller = min(first, second); larger = max(first, second)
        value   = 1.0_wp/(14.0_wp*larger**7)-smaller/(18.0_wp*larger**8)
    end function reciprocal_power_n2_m6

    !> @brief Evaluate the n=2, m=5 reciprocal-power radial kernel.
    pure function reciprocal_power_n2_m5( first, second ) result(value)
        real(wp), intent(in) :: first, second
        real(wp)             :: value, larger, smaller
        smaller = min(first, second); larger = max(first, second)
        value   = 4.0_wp*larger**(-6)/(6.0_wp*7.0_wp)*(1.0_wp-smaller/larger*6.0_wp/8.0_wp)
    end function reciprocal_power_n2_m5

    !> @brief Differentiate the n=2, m=6 reciprocal-power kernel by its first coordinate.
    pure function reciprocal_power_n2_m6_derivative( first, second ) result(value)
        real(wp), intent(in) :: first, second
        real(wp)             :: value
        if( first <= second )then
            value = -1.0_wp/(18.0_wp*second**8)
        else
            value = 4.0_wp*second/(9.0_wp*first**9)-1.0_wp/(2.0_wp*first**8)
        endif
    end function reciprocal_power_n2_m6_derivative
end module chempotpy_rkhs
