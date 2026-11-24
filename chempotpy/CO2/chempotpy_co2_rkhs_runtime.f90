! System:                        CO2
! Functional form:               reproducing-kernel Hilbert-space interpolation
! Common name:                   CO2 five-state RKHS runtime
! Reference:                     J. C. San Vicente Veliz et al., PCCP 23, 11251 (2021)
! Number of derivatives:         1
! Number of bodies:              3
! Number of electronic surfaces: 1 per owning module
! Interface:                     shared implementation
! Data file:                     none; every fitted coefficient is embedded
! Notes:                         Atom order at the Cartesian boundary is O, C, O.

!> @brief Evaluate the shared global-coordinate construction for the five CO2 RKHS surfaces.
module chempotpy_co2_rkhs_runtime
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_co2_rkhs_surface
    real(wp), parameter :: carbon_mass = 12.00000_wp
    real(wp), parameter :: oxygen_mass = 15.99491462_wp

    contains

    !> @brief Adapt one immutable CO2 RKHS fit to the CHEMPOTPY Cartesian contract.
    subroutine evaluate_co2_rkhs_surface( x, igrad, p, g, d, angular_grid, oo_first_grid, oo_second_grid, &
            oo_coefficients, co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, status )
        use chempotpy_internal_coordinates, only: finalize_legacy_triatomic_surface, &
            & prepare_legacy_triatomic_coordinates
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(3,3), angular_grid(:), oo_first_grid(:), oo_second_grid(:)
        real(wp), intent(in)  :: oo_coefficients(:), co_first_grid(:), co_second_grid(:), co_coefficients(:)
        real(wp), intent(in)  :: switching_scales(3), channel_offset
        real(wp), intent(out) :: p(1), g(1,3,3), d(1,1,3,3)
        real(wp)              :: co2_distances(3), co2_jacobian(3,9), distance_jacobian(3,9), distances(3)
        real(wp)              :: energy, internal_gradient(3)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad < 0 .or. igrad > 1 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'CO2 RKHS surfaces provide energy and first derivatives')
            return
        endif
        call prepare_legacy_triatomic_coordinates(x, distances, distance_jacobian)
        co2_distances     = distances([3,1,2])
        co2_jacobian(1,:) = distance_jacobian(3,:)
        co2_jacobian(2,:) = distance_jacobian(1,:)
        co2_jacobian(3,:) = distance_jacobian(2,:)
        call evaluate_internal(co2_distances, angular_grid, oo_first_grid, oo_second_grid, oo_coefficients, &
            & co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, energy, &
            & internal_gradient)
        call finalize_legacy_triatomic_surface(igrad, energy, internal_gradient, co2_jacobian, p, g, d)
    end subroutine evaluate_co2_rkhs_surface

    !> @brief Evaluate one global CO2 surface in distance coordinates and atomic units.
    subroutine evaluate_internal( distances, angular_grid, oo_first_grid, oo_second_grid, oo_coefficients, &
            co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, energy, derivative )
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(wp), intent(in)  :: distances(3), angular_grid(:), oo_first_grid(:), oo_second_grid(:)
        real(wp), intent(in)  :: oo_coefficients(:), co_first_grid(:), co_second_grid(:), co_coefficients(:)
        real(wp), intent(in)  :: switching_scales(3), channel_offset
        real(wp), intent(out) :: energy, derivative(3)
        call blend_channels(distances, angular_grid, oo_first_grid, oo_second_grid, oo_coefficients, co_first_grid, &
            & co_second_grid, co_coefficients, switching_scales, channel_offset, energy, derivative)
        if( .not. all(ieee_is_finite(derivative)) )then
            call finite_difference_derivative(distances, angular_grid, oo_first_grid, oo_second_grid, oo_coefficients, &
                & co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, derivative)
        endif
    end subroutine evaluate_internal

    !> @brief Blend the three arrangement-channel RKHS evaluations and their distance derivatives.
    subroutine blend_channels( distances, angular_grid, oo_first_grid, oo_second_grid, oo_coefficients, &
            co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, energy, derivative )
        real(wp), intent(in)  :: distances(3), angular_grid(:), oo_first_grid(:), oo_second_grid(:)
        real(wp), intent(in)  :: oo_coefficients(:), co_first_grid(:), co_second_grid(:), co_coefficients(:)
        real(wp), intent(in)  :: switching_scales(3), channel_offset
        real(wp), intent(out) :: energy, derivative(3)
        real(wp)              :: arrangement_derivative(3,3), arrangement_energy(3), local_derivative(3)
        real(wp)              :: local_distances(3), masses(3), transform(3,3), weight(3), weight_derivative(3,3)
        local_distances = [distances(1), distances(3), distances(2)]
        masses          = [oxygen_mass, oxygen_mass, carbon_mass]
        call evaluate_arrangement(local_distances, masses, 1, angular_grid, oo_first_grid, oo_second_grid, &
            & oo_coefficients, local_derivative, arrangement_energy(1), transform)
        arrangement_derivative(1,1) = dot_product(local_derivative, transform(:,1))
        arrangement_derivative(1,2) = dot_product(local_derivative, transform(:,3))
        arrangement_derivative(1,3) = dot_product(local_derivative, transform(:,2))
        local_distances             = [distances(2), distances(3), distances(1)]
        masses                      = [oxygen_mass, carbon_mass, oxygen_mass]
        call evaluate_arrangement(local_distances, masses, 2, angular_grid, co_first_grid, co_second_grid, &
            & co_coefficients, local_derivative, arrangement_energy(2), transform)
        arrangement_energy(2)       = arrangement_energy(2)+channel_offset
        arrangement_derivative(2,1) = dot_product(local_derivative, transform(:,3))
        arrangement_derivative(2,2) = dot_product(local_derivative, transform(:,1))
        arrangement_derivative(2,3) = dot_product(local_derivative, transform(:,2))
        local_distances             = [distances(3), distances(2), distances(1)]
        call evaluate_arrangement(local_distances, masses, 2, angular_grid, co_first_grid, co_second_grid, &
            & co_coefficients, local_derivative, arrangement_energy(3), transform)
        arrangement_energy(3)       = arrangement_energy(3)+channel_offset
        arrangement_derivative(3,1) = dot_product(local_derivative, transform(:,3))
        arrangement_derivative(3,2) = dot_product(local_derivative, transform(:,2))
        arrangement_derivative(3,3) = dot_product(local_derivative, transform(:,1))
        call switching_weights(distances, switching_scales, weight, weight_derivative)
        energy = dot_product(arrangement_energy, weight)
        derivative(1) = dot_product(arrangement_derivative(:,1), weight) + dot_product(weight_derivative(:,1), &
            & arrangement_energy)
        derivative(2) = dot_product(arrangement_derivative(:,2), weight) + dot_product(weight_derivative(:,2), &
            & arrangement_energy)
        derivative(3) = dot_product(arrangement_derivative(:,3), weight) + dot_product(weight_derivative(:,3), &
            & arrangement_energy)
    end subroutine blend_channels

    !> @brief Evaluate one Jacobi-coordinate arrangement and transform its derivative to distances.
    subroutine evaluate_arrangement( distances, masses, asymptote, angular_grid, first_grid, second_grid, &
            coefficients, distance_derivative, energy, transform )
        use chempotpy_rkhs, only: evaluate_ts2_rp26_3d
        real(wp), intent(in)  :: distances(3), masses(3), angular_grid(:), first_grid(:), second_grid(:)
        real(wp), intent(in)  :: coefficients(:)
        integer, intent(in)   :: asymptote
        real(wp), intent(out) :: distance_derivative(3), energy, transform(3,3)
        real(wp)              :: asymptotic_derivative, asymptotic_value, jacobi_derivative(3)
        real(wp)              :: point(3), radial_distance, theta
        call jacobi_coordinates(distances, masses, radial_distance, theta, transform)
        point = [(1.0_wp-cos(theta))/2.0_wp, radial_distance, distances(1)]
        call evaluate_ts2_rp26_3d(point, angular_grid, first_grid, second_grid, coefficients, energy, jacobi_derivative)
        call asymptotic_correction(distances(1), asymptote, asymptotic_value, asymptotic_derivative)
        energy                 = energy+asymptotic_value
        distance_derivative(1) = jacobi_derivative(2)
        distance_derivative(2) = asymptotic_derivative+jacobi_derivative(3)
        distance_derivative(3) = jacobi_derivative(1)*sin(theta)/2.0_wp
    end subroutine evaluate_arrangement

    !> @brief Convert three distances to published Jacobi coordinates and their distance Jacobian.
    pure subroutine jacobi_coordinates( distances, masses, radial_distance, theta, derivative )
        real(wp), intent(in)  :: distances(3), masses(3)
        real(wp), intent(out) :: radial_distance, theta, derivative(3,3)
        real(wp)              :: center_first, center_second, cosine_theta, first_fraction, second_fraction
        first_fraction  = masses(2)/(masses(1)+masses(2))
        second_fraction = masses(1)/(masses(1)+masses(2))
        center_first    = distances(1)*first_fraction
        center_second   = distances(1)*second_fraction
        radial_distance = sqrt(distances(2)**2/distances(1)*center_first + distances(3)**2/distances(1)*center_second- &
            & center_first*center_second)
        radial_distance = max(radial_distance, sqrt(epsilon(1.0_wp)))
        cosine_theta    = (center_second**2+radial_distance**2-distances(2)**2)/ (2.0_wp*center_second*radial_distance)
        cosine_theta    = max(-1.0_wp, min(1.0_wp, cosine_theta)); theta = acos(cosine_theta)
        derivative(1,1) = -first_fraction*second_fraction*distances(1)/radial_distance
        derivative(1,2) = distances(2)*first_fraction/radial_distance
        derivative(1,3) = distances(3)*second_fraction/radial_distance
        derivative(2,:) = [1.0_wp, 0.0_wp, 0.0_wp]
        derivative(3,1) = (derivative(1,1)/radial_distance*cosine_theta+cosine_theta/distances(1) - (radial_distance* &
            & derivative(1,1)+center_second*second_fraction)/ (center_second*radial_distance))/sqrt(1.0_wp- &
            & cosine_theta**2)
        derivative(3,2) = (distances(2)/(center_second*radial_distance)-derivative(1,2)/center_second + cosine_theta/ &
            & radial_distance*derivative(1,2))/sqrt(1.0_wp-cosine_theta**2)
        derivative(3,3) = (cosine_theta/radial_distance*derivative(1,3)-derivative(1,3)/center_second)/ sqrt(1.0_wp- &
            & cosine_theta**2)
    end subroutine jacobi_coordinates

    !> @brief Evaluate the released arrangement switching weights and their derivatives.
    pure subroutine switching_weights( distances, scales, weight, derivative )
        real(wp), intent(in)  :: distances(3), scales(3)
        real(wp), intent(out) :: weight(3), derivative(3,3)
        real(wp)              :: adjusted(3), exponential(3), exponential_sum
        integer               :: first, second
        adjusted        = distances
        exponential_sum = 0.0_wp
        do while( exponential_sum < epsilon(1.0_wp) )
            exponential     = exp(-(adjusted/scales)**2)
            exponential_sum = sum(exponential)
            if( exponential_sum < epsilon(1.0_wp) ) adjusted = adjusted-0.1_wp
        enddo
        weight = exponential/exponential_sum
        do first = 1, 3
            do second = 1, 3
                if( first == second )then
                    derivative(first,second) = -2.0_wp*(adjusted(second)/scales(second))*exponential(first)* &
                        & (exponential_sum-exponential(first))/(scales(second)*exponential_sum**2)
                else
                    derivative(first,second) = 2.0_wp*(adjusted(second)/scales(second))*exponential(first)* &
                        & exponential(second)/(scales(second)*exponential_sum**2)
                endif
            enddo
        enddo
    end subroutine switching_weights

    !> @brief Evaluate the shared long-range correction for one asymptotic channel.
    pure subroutine asymptotic_correction( distance, asymptote, energy, derivative )
        use chempotpy_co2_asymptotic_parameters, only: co_asymptotic_coefficients, co_asymptotic_grid, &
            & oo_asymptotic_coefficients, oo_asymptotic_grid
        integer, intent(in)   :: asymptote
        real(wp), intent(in)  :: distance
        real(wp), intent(out) :: energy, derivative
        integer               :: index
        energy = 0.0_wp; derivative = 0.0_wp
        if( asymptote == 1 )then
            do index = 1, size(oo_asymptotic_grid)
                energy     = energy+radial_correction(distance, oo_asymptotic_grid(index))* &
                    & oo_asymptotic_coefficients(index)
                derivative = derivative+radial_correction_derivative(distance, oo_asymptotic_grid(index))* &
                    & oo_asymptotic_coefficients(index)
            enddo
        else
            do index = 1, size(co_asymptotic_grid)
                energy     = energy+radial_correction(distance, co_asymptotic_grid(index))* &
                    & co_asymptotic_coefficients(index)
                derivative = derivative+radial_correction_derivative(distance, co_asymptotic_grid(index))* &
                    & co_asymptotic_coefficients(index)
            enddo
        endif
    end subroutine asymptotic_correction

    !> @brief Evaluate the n=2, m=6 radial correction kernel.
    pure function radial_correction( first, second ) result(value)
        real(wp), intent(in) :: first, second
        real(wp)             :: value, larger, smaller
        larger = max(first, second); smaller = min(first, second)
        value  = 1.0_wp/(14.0_wp*larger**7)-smaller/(18.0_wp*larger**8)
    end function radial_correction

    !> @brief Differentiate the radial correction kernel by its first coordinate.
    pure function radial_correction_derivative( first, second ) result(value)
        real(wp), intent(in) :: first, second
        real(wp)             :: value
        if( first < second )then
            value = -1.0_wp/(18.0_wp*second**8)
        else
            value = -1.0_wp/(2.0_wp*first**8)+4.0_wp*second/(9.0_wp*first**9)
        endif
    end function radial_correction_derivative

    !> @brief Recover distance derivatives at coordinate singularities with Richardson differences.
    subroutine finite_difference_derivative( distances, angular_grid, oo_first_grid, oo_second_grid, &
            oo_coefficients, co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, &
            derivative )
        real(wp), intent(in)  :: distances(3), angular_grid(:), oo_first_grid(:), oo_second_grid(:)
        real(wp), intent(in)  :: oo_coefficients(:), co_first_grid(:), co_second_grid(:), co_coefficients(:)
        real(wp), intent(in)  :: switching_scales(3), channel_offset
        real(wp), intent(out) :: derivative(3)
        real(wp), parameter   :: step = 0.005_wp
        real(wp)              :: displaced(3), dummy(3), energy(4), first_difference, second_difference
        integer               :: coordinate
        do coordinate = 1, 3
            displaced = distances; displaced(coordinate) = distances(coordinate)-2.0_wp*step
            call blend_channels(displaced, angular_grid, oo_first_grid, oo_second_grid, oo_coefficients, &
                & co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, energy(1), dummy)
            displaced = distances; displaced(coordinate) = distances(coordinate)-step
            call blend_channels(displaced, angular_grid, oo_first_grid, oo_second_grid, oo_coefficients, &
                & co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, energy(2), dummy)
            displaced = distances; displaced(coordinate) = distances(coordinate)+step
            call blend_channels(displaced, angular_grid, oo_first_grid, oo_second_grid, oo_coefficients, &
                & co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, energy(3), dummy)
            displaced = distances; displaced(coordinate) = distances(coordinate)+2.0_wp*step
            call blend_channels(displaced, angular_grid, oo_first_grid, oo_second_grid, oo_coefficients, &
                & co_first_grid, co_second_grid, co_coefficients, switching_scales, channel_offset, energy(4), dummy)
            first_difference       = (energy(3)-energy(2))/(2.0_wp*step)
            second_difference      = (energy(4)-energy(1))/(4.0_wp*step)
            derivative(coordinate) = (4.0_wp*first_difference-second_difference)/3.0_wp
        enddo
    end subroutine finite_difference_derivative
end module chempotpy_co2_rkhs_runtime
