!> @brief Evaluate the shared three-dimensional Ne--NO interaction fit and vibrational averages.
module chempotpy_neno_runtime
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_neno_vibrational_components

    contains

    !> @brief Evaluate Vsum and Vdif in inverse centimeters for one NO vibrational state.
    subroutine evaluate_neno_vibrational_components( distance_bohr, angle_degree, vibrational_index, &
            v_sum, v_difference )
        use chempotpy_interpolation, only: cubic_spline_polynomial_coefficients, &
            evaluate_cubic_spline_polynomial
        use chempotpy_neno_parameters, only: no_vibrational_normalization, wavefunction_distance_bohr, &
            wavefunction_value_v0, wavefunction_value_v2
        use chempotpy_quadrature, only: gauss_legendre_rule
        real(wp), intent(in)  :: distance_bohr, angle_degree
        integer, intent(in)   :: vibrational_index
        real(wp), intent(out) :: v_sum, v_difference
        integer, parameter    :: quadrature_order = 50
        real(wp)              :: first(56), second(56), third(56), wavefunction_values(56)
        real(wp)              :: node_sum, node_difference, nodes(quadrature_order), wavefunction
        real(wp)              :: weights(quadrature_order)
        integer               :: index
        if( vibrational_index == 1 )then
            wavefunction_values = wavefunction_value_v0
        elseif( vibrational_index == 2 )then
            wavefunction_values = wavefunction_value_v2
        else
            v_sum = 0.0_wp; v_difference = 0.0_wp
            return
        endif
        call cubic_spline_polynomial_coefficients(wavefunction_distance_bohr, wavefunction_values, &
            first, second, third)
        call gauss_legendre_rule(1.8_wp, 2.7_wp, nodes, weights)
        v_sum = 0.0_wp; v_difference = 0.0_wp
        do index=1, quadrature_order
            wavefunction = evaluate_cubic_spline_polynomial(wavefunction_distance_bohr, wavefunction_values, first, &
                & second, third, nodes(index))
            call evaluate_three_dimensional_components(distance_bohr, nodes(index), angle_degree, &
                node_sum, node_difference)
            v_sum        = v_sum+weights(index)*node_sum*wavefunction**2
            v_difference = v_difference+weights(index)*node_difference*wavefunction**2
        enddo
        v_sum        = v_sum/no_vibrational_normalization(vibrational_index)
        v_difference = v_difference/no_vibrational_normalization(vibrational_index)
    end subroutine evaluate_neno_vibrational_components

    !> @brief Interpolate the Ne--NO fit between its three published NO bond distances.
    subroutine evaluate_three_dimensional_components( distance_bohr, no_distance_bohr, angle_degree, &
            v_sum, v_difference )
        use chempotpy_neno_parameters, only: no_equilibrium_distance_bohr, no_reference_distances_bohr
        real(wp), intent(in)  :: distance_bohr, no_distance_bohr, angle_degree
        real(wp), intent(out) :: v_sum, v_difference
        real(wp)              :: coefficient_1, coefficient_2, determinant, difference_1, difference_3
        real(wp)              :: fixed_difference(3), fixed_sum(3), z, z_1, z_3
        integer               :: bond_index
        do bond_index=1, 3
            call evaluate_fixed_bond_components(distance_bohr, angle_degree, bond_index, &
                fixed_sum(bond_index), fixed_difference(bond_index))
        enddo
        z_1           = (no_reference_distances_bohr(1)-no_equilibrium_distance_bohr)/no_equilibrium_distance_bohr
        z_3           = (no_reference_distances_bohr(3)-no_equilibrium_distance_bohr)/no_equilibrium_distance_bohr
        determinant   = z_1*z_3**2-z_3*z_1**2
        z             = (no_distance_bohr-no_equilibrium_distance_bohr)/no_equilibrium_distance_bohr
        difference_1  = fixed_sum(1)-fixed_sum(2)
        difference_3  = fixed_sum(3)-fixed_sum(2)
        coefficient_1 = (difference_1*z_3**2-difference_3*z_1**2)/determinant
        coefficient_2 = (z_1*difference_3-z_3*difference_1)/determinant
        v_sum         = fixed_sum(2)+coefficient_1*z+coefficient_2*z**2
        difference_1  = fixed_difference(1)-fixed_difference(2)
        difference_3  = fixed_difference(3)-fixed_difference(2)
        coefficient_1 = (difference_1*z_3**2-difference_3*z_1**2)/determinant
        coefficient_2 = (z_1*difference_3-z_3*difference_1)/determinant
        v_difference  = fixed_difference(2)+coefficient_1*z+coefficient_2*z**2
    end subroutine evaluate_three_dimensional_components

    !> @brief Evaluate angular expansions at one of the three fitted NO bond distances.
    subroutine evaluate_fixed_bond_components( distance_bohr, angle_degree, bond_index, v_sum, v_difference )
        use chempotpy_constants,       only: PI
        use chempotpy_neno_parameters, only: difference_angular_inverse, radial_fit_a_double_prime_r2, &
            radial_fit_a_prime_r2, radial_fit_a_double_prime_re, radial_fit_a_prime_re, &
            radial_fit_a_double_prime_r4, radial_fit_a_prime_r4, sum_angular_inverse
        use chempotpy_special_functions, only: associated_legendre_polynomial
        real(wp), intent(in)  :: distance_bohr, angle_degree
        integer, intent(in)   :: bond_index
        real(wp), intent(out) :: v_sum, v_difference
        ! Preserve the release's rounded microhartree-to-wavenumber factor for numerical compatibility.
        real(wp), parameter :: released_microhartree_to_inverse_centimeter = 1.0_wp/4.556335_wp
        real(wp)            :: basis_difference(11), basis_sum(13), coefficients_difference(11)
        real(wp)            :: coefficients_sum(13), cosine_theta, difference_values(11), state_values(2,13)
        real(wp)            :: sum_values(13)
        integer             :: angle_index, degree
        do angle_index=1, 13
            select case( bond_index )
            case( 1 )
                state_values(1,angle_index) = evaluate_radial_fit( radial_fit_a_double_prime_r2(:, angle_index), &
                    & distance_bohr)
                state_values(2,angle_index) = evaluate_radial_fit( radial_fit_a_prime_r2(:, angle_index), distance_bohr)
            case( 2 )
                state_values(1,angle_index) = evaluate_radial_fit( radial_fit_a_double_prime_re(:, angle_index), &
                    & distance_bohr)
                state_values(2,angle_index) = evaluate_radial_fit( radial_fit_a_prime_re(:, angle_index), distance_bohr)
            case( 3 )
                state_values(1,angle_index) = evaluate_radial_fit( radial_fit_a_double_prime_r4(:, angle_index), &
                    & distance_bohr)
                state_values(2,angle_index) = evaluate_radial_fit( radial_fit_a_prime_r4(:, angle_index), distance_bohr)
            end select
            sum_values(angle_index) = 0.5_wp*sum(state_values(:,angle_index))
        enddo
        difference_values       = 0.5_wp*(state_values(1,2:12)-state_values(2,2:12))
        coefficients_sum        = matmul(sum_angular_inverse, sum_values)
        coefficients_difference = matmul(difference_angular_inverse, difference_values)
        cosine_theta            = cos(angle_degree*PI/180.0_wp)
        do degree=0, 12
            basis_sum(degree+1) = associated_legendre_polynomial(degree, 0, cosine_theta)
        enddo
        do degree=2, 12
            basis_difference(degree-1) = associated_legendre_polynomial(degree, 2, cosine_theta)/ &
                & sqrt(real((degree-1)*degree*(degree+1)*(degree+2), wp))
        enddo
        v_sum        = dot_product(coefficients_sum, basis_sum)
        v_difference = dot_product(coefficients_difference, basis_difference)
        v_sum        = v_sum*released_microhartree_to_inverse_centimeter
        v_difference = v_difference*released_microhartree_to_inverse_centimeter
    end subroutine evaluate_fixed_bond_components

    !> @brief Evaluate one published angle-specific radial fit.
    pure function evaluate_radial_fit( parameters, input_distance_bohr ) result(energy)
        real(wp), intent(in) :: parameters(16), input_distance_bohr
        real(wp)             :: energy
        real(wp)             :: distance_bohr, dispersion, polynomial, switching
        integer              :: power
        distance_bohr = max(3.0_wp, input_distance_bohr)
        polynomial    = 0.0_wp
        do power=0, 8
            polynomial = polynomial+parameters(power+3)*distance_bohr**power
        enddo
        switching = 0.5_wp*(1.0_wp+tanh(1.0_wp+parameters(12)*distance_bohr))
        dispersion = parameters(13)/distance_bohr**6+parameters(14)/distance_bohr**7 + parameters(15)/distance_bohr** &
            & 8+parameters(16)/distance_bohr**9
        energy = exp(-parameters(1)*distance_bohr-parameters(2))*polynomial-switching*dispersion
    end function evaluate_radial_fit
end module chempotpy_neno_runtime
