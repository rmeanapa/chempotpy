!> @brief Evaluate the common fixed-Ni(111) EAM6 and EAM7 H-Ni algebra.
module chempotpy_hni_eam67_runtime
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_fixed_hydrogen_eam
    integer, parameter  :: nickel_type                  = 1, hydrogen_type = 2
    real(wp), parameter :: inner_cutoff_angstrom        = 5.0_wp, outer_cutoff_angstrom = 10.0_wp
    real(wp), parameter :: electrostatic_conversion     = 14.3888_wp
    integer, parameter  :: nickel_principal_numbers(10) = [1, 1, 2, 2, 3, 3, 4, 4, 3, 3]
    real(wp), parameter :: nickel_factorials(10)        = [ &
        & 2.0_wp, 2.0_wp, 24.0_wp, 24.0_wp, 720.0_wp, 720.0_wp, 40320.0_wp, 40320.0_wp, 720.0_wp, 720.0_wp ]
    real(wp), parameter :: nickel_coefficients(10) = [ &
        & -0.00389_wp, -0.02991_wp, -0.03189_wp, 0.15289_wp, -0.20048_wp, -0.05423_wp, 0.49292_wp, &
        & 0.61875_wp, 0.4212_wp, 0.70658_wp ]
    real(wp), parameter :: nickel_exponents(10) = [ &
        & 54.88885_wp, 38.48431_wp, 27.42703_wp, 20.88204_wp, 10.95707_wp, 7.31958_wp, 3.9265_wp, &
        & 2.15289_wp, 12.67582_wp, 5.43253_wp ]
    real(wp), parameter :: nickel_radial_coefficients(10) = nickel_coefficients* &
        & (2.0_wp*nickel_exponents)**(real(nickel_principal_numbers, wp)+0.5_wp)/sqrt(nickel_factorials)
    real(wp), parameter :: hydrogen_radial_coefficient = &
        & (2.0_wp*1.88972652_wp)**1.5_wp/sqrt(2.0_wp)

    contains

    !> @brief Evaluate one EAM6 or EAM7 energy and analytic active-H gradient in eV and eV/angstrom.
    pure subroutine evaluate_fixed_hydrogen_eam( coordinates, model_version, energy, hydrogen_gradient )
        real(wp), intent(in)  :: coordinates(:,:)
        integer, intent(in)   :: model_version
        real(wp), intent(out) :: energy, hydrogen_gradient(3)
        real(wp)              :: auxiliary_cubic, auxiliary_cubic_derivative
        real(wp)              :: auxiliary_root, auxiliary_root_derivative
        real(wp)              :: density_derivative, density_value, displacement(3), distance
        real(wp)              :: embedding_energy, h_aux_cubic_derivative
        real(wp)              :: h_aux_root_derivative, h_embedding_derivative, pair_derivative, pair_energy
        real(wp)              :: radial_derivative, raw_density, raw_density_derivative
        real(wp)              :: densities(size(coordinates, 1)), embedding_derivatives(size(coordinates, 1))
        integer               :: atom, first_atom, first_type, second_atom
        densities       = 0.0_wp
        auxiliary_cubic = 0.0_wp
        auxiliary_root  = 0.0_wp
        energy          = 0.0_wp
        do first_atom=1, size(coordinates, 1)-1
            first_type = merge(hydrogen_type, nickel_type, first_atom == 1)
            do second_atom=first_atom+1, size(coordinates, 1)
                displacement = coordinates(first_atom,:)-coordinates(second_atom,:)
                distance     = sqrt(sum(displacement*displacement))
                if( distance >= outer_cutoff_angstrom ) cycle
                call atomic_density(nickel_type, distance, coordinates(second_atom,3), density_value, &
                    & density_derivative, raw_density, raw_density_derivative)
                densities(first_atom) = densities(first_atom)+density_value
                if( first_atom == 1 .and. raw_density > 0.0_wp )then
                    auxiliary_cubic = auxiliary_cubic+(52.0834_wp*raw_density)**2
                    auxiliary_root  = auxiliary_root+sqrt(30.0_wp*raw_density)
                endif
                call atomic_density(first_type, distance, coordinates(first_atom,3), density_value, &
                    & density_derivative, raw_density, raw_density_derivative)
                densities(second_atom) = densities(second_atom)+density_value
                if( first_atom == 1 )then
                    call pair_potential(model_version, hydrogen_type, nickel_type, distance, pair_energy, pair_derivative)
                    energy = energy+pair_energy
                endif
            enddo
        enddo
        do atom=2, size(coordinates, 1)
            call nickel_embedding(densities(atom), embedding_energy, embedding_derivatives(atom))
            energy = energy+embedding_energy
        enddo
        call hydrogen_embedding(model_version, densities(1), auxiliary_cubic, auxiliary_root, embedding_energy, &
            & h_embedding_derivative, h_aux_cubic_derivative, h_aux_root_derivative)
        embedding_derivatives(1) = h_embedding_derivative
        energy                   = energy+embedding_energy
        hydrogen_gradient        = 0.0_wp
        do second_atom=2, size(coordinates, 1)
            displacement = coordinates(1,:)-coordinates(second_atom,:)
            distance     = sqrt(sum(displacement*displacement))
            if( distance >= outer_cutoff_angstrom ) cycle
            call atomic_density(nickel_type, distance, coordinates(second_atom,3), density_value, density_derivative, &
                & raw_density, raw_density_derivative)
            radial_derivative = h_embedding_derivative*density_derivative
            if( raw_density > 0.0_wp )then
                auxiliary_cubic_derivative = 2.0_wp*(52.0834_wp*raw_density)*raw_density_derivative*52.0834_wp
                auxiliary_root_derivative  = 0.5_wp*(30.0_wp*raw_density)**(-0.5_wp)*raw_density_derivative*30.0_wp
                radial_derivative = radial_derivative+h_aux_cubic_derivative*auxiliary_cubic_derivative+ &
                    & h_aux_root_derivative*auxiliary_root_derivative
            endif
            call atomic_density(hydrogen_type, distance, coordinates(1,3), density_value, density_derivative, &
                & raw_density, raw_density_derivative)
            radial_derivative = radial_derivative+embedding_derivatives(second_atom)*density_derivative
            call pair_potential(model_version, hydrogen_type, nickel_type, distance, pair_energy, pair_derivative)
            radial_derivative = radial_derivative+pair_derivative
            hydrogen_gradient = hydrogen_gradient+radial_derivative*displacement/distance
        enddo
    end subroutine evaluate_fixed_hydrogen_eam

    !> @brief Evaluate an EAM atomic density, its switched value, and radial derivatives.
    pure subroutine atomic_density( atom_type, distance, source_z, density, derivative, raw_density, raw_derivative )
        use chempotpy_constants, only: PI
        integer, intent(in)   :: atom_type
        real(wp), intent(in)  :: distance, source_z
        real(wp), intent(out) :: density, derivative, raw_density, raw_derivative
        real(wp)              :: cutoff_derivative, cutoff_value, d_sum_d, d_sum_s
        real(wp)              :: electron_count_d, electron_count_s, sum_d, sum_s, term
        integer               :: orbital
        call cutoff(distance, cutoff_value, cutoff_derivative)
        sum_s = 0.0_wp; d_sum_s = 0.0_wp; sum_d = 0.0_wp; d_sum_d = 0.0_wp
        if( atom_type == nickel_type )then
            do orbital=1, 8
                term = nickel_radial_coefficients(orbital)*distance**(nickel_principal_numbers(orbital)-1)* exp(- &
                    & nickel_exponents(orbital)*distance)
                sum_s = sum_s+term
                d_sum_s = d_sum_s+term*(real(nickel_principal_numbers(orbital)-1, wp)/distance- &
                    & nickel_exponents(orbital))
            enddo
            do orbital=9, 10
                term = nickel_radial_coefficients(orbital)*distance**(nickel_principal_numbers(orbital)-1)* exp(- &
                    & nickel_exponents(orbital)*distance)
                sum_d = sum_d+term
                d_sum_d = d_sum_d+term*(real(nickel_principal_numbers(orbital)-1, wp)/distance- &
                    & nickel_exponents(orbital))
            enddo
            electron_count_s = 2.0_wp
            if( abs(source_z) < 0.5_wp ) electron_count_s = 2.0_wp
            electron_count_d = 10.0_wp-electron_count_s
        else
            term             = hydrogen_radial_coefficient*exp(-1.88972652_wp*distance)
            sum_s            = term
            d_sum_s          = -1.88972652_wp*term
            electron_count_s = 1.0_wp
            electron_count_d = 0.0_wp
        endif
        raw_derivative = (electron_count_s*sum_s*d_sum_s+electron_count_d*sum_d*d_sum_d)/(2.0_wp*PI)
        raw_density    = (electron_count_s*sum_s**2+electron_count_d*sum_d**2)/(4.0_wp*PI)
        derivative     = raw_derivative*cutoff_value+cutoff_derivative*raw_density
        density        = raw_density*cutoff_value
    end subroutine atomic_density

    !> @brief Evaluate the EAM6/EAM7 nickel embedding function and density derivative.
    pure subroutine nickel_embedding( density, energy, derivative )
        use chempotpy_hni_eam67_parameters, only: parameters => nickel_embedding_parameters
        real(wp), intent(in)  :: density
        real(wp), intent(out) :: energy, derivative
        real(wp)              :: delta, delta_squared, scaled_density, scaled_power, term
        if( density > parameters(8) )then
            energy     = parameters(12)
            derivative = 0.0_wp
        elseif( density > parameters(7) )then
            delta         = density-parameters(8)
            delta_squared = delta*delta
            energy        = parameters(12)+delta_squared*delta*(parameters(11)+delta*(parameters(10)+delta* &
                & parameters(9)))
            derivative    = delta_squared*(3.0_wp*parameters(11)+delta*(4.0_wp*parameters(10)+ delta*5.0_wp* &
                & parameters(9)))
        else
            scaled_density = -parameters(2)*density
            term           = parameters(1)*exp(scaled_density)
            energy         = term*density
            derivative     = term*(1.0_wp+scaled_density)
            scaled_density = -parameters(4)*density
            term           = parameters(3)*density
            energy         = energy+term**3*exp(scaled_density)
            derivative     = derivative+term**2*exp(scaled_density)*(3.0_wp*parameters(3)-parameters(4)*term)
            scaled_density = -parameters(6)*density
            term           = parameters(5)*exp(scaled_density)
            energy         = energy+term*density
            derivative     = derivative+term*(1.0_wp+scaled_density)
            scaled_power   = parameters(13)*density
            term           = exp(-parameters(14)*density)
            energy         = energy+scaled_power**5*term
            derivative     = derivative+term*(5.0_wp*parameters(13)*scaled_power**4- parameters(14)*scaled_power**5)
        endif
    end subroutine nickel_embedding

    !> @brief Evaluate the version-specific H embedding and auxiliary-density corrections.
    pure subroutine hydrogen_embedding( model_version, density, cubic_density, root_density, energy, derivative, &
            cubic_derivative, root_derivative )
        integer, intent(in)   :: model_version
        real(wp), intent(in)  :: density, cubic_density, root_density
        real(wp), intent(out) :: energy, derivative, cubic_derivative, root_derivative
        real(wp)              :: exponential, term
        if( model_version == 6 )then
            exponential = exp(-5.0715_wp*density)
            energy      = 476.1208_wp*exponential*density
            derivative  = 476.1208_wp*exponential*(1.0_wp-5.0715_wp*density)
            exponential = exp(-5.2849_wp*density)
            energy      = energy-545.3935_wp*exponential*density
            derivative  = derivative-545.3935_wp*exponential*(1.0_wp-5.2849_wp*density)
        else
            call eam7_hydrogen_embedding(density, energy, derivative)
        endif
        cubic_derivative = 0.0_wp
        root_derivative  = 0.0_wp
        if( (model_version == 6 .and. density > 0.0_wp) .or. &
            & (model_version == 7 .and. density > 0.5_wp) )then
        call symmetric_descriptor(-0.081682_wp, 1.0_wp, 5.0_wp, cubic_density, term, cubic_derivative)
        energy = energy+term
        call symmetric_descriptor(-0.006505_wp, 10.0_wp, 4.4_wp, root_density, term, root_derivative)
        energy = energy+term
        endif
    end subroutine hydrogen_embedding

    !> @brief Evaluate the EAM7 piecewise H embedding function and density derivative.
    pure subroutine eam7_hydrogen_embedding( density, energy, derivative )
        real(wp), intent(in)  :: density
        real(wp), intent(out) :: energy, derivative
        real(wp), parameter   :: boundary = 0.009_wp, width = 0.007_wp
        real(wp)              :: exponential, lower_offset, upper_offset
        if( density <= boundary-width )then
            if( density > 0.0_wp )then
                exponential = exp(-34.462_wp*density)
                energy      = -4.8029_wp*density**0.38553_wp*exponential
                derivative  = -4.8029_wp*density**(0.38553_wp-1.0_wp)*exponential* (0.38553_wp-34.462_wp*density)
            else
                energy = 0.0_wp; derivative = 0.0_wp
            endif
        elseif( density <= boundary )then
            upper_offset = density-boundary
            lower_offset = density-(boundary-width)
            energy = 0.1840945234e12_wp*upper_offset**5+0.3276017408e10_wp*upper_offset**4+ 15102027.0_wp* &
                & upper_offset**3-0.1303885824e12_wp*lower_offset**5+ 0.2375889152e10_wp*lower_offset**4- &
                & 11946214.0_wp*lower_offset**3
            derivative = 5.0_wp*0.1840945234e12_wp*upper_offset**4+4.0_wp*0.3276017408e10_wp* upper_offset**3+3.0_wp* &
                & 15102027.0_wp*upper_offset**2-5.0_wp*0.1303885824e12_wp* lower_offset**4+4.0_wp*0.2375889152e10_wp* &
                & lower_offset**3-3.0_wp*11946214.0_wp* lower_offset**2
        else
            exponential = exp(-7.2784_wp*density)
            energy      = -74.556_wp*density**1.0154_wp*exponential
            derivative  = -74.556_wp*density**(1.0154_wp-1.0_wp)*exponential*(1.0154_wp-7.2784_wp*density)
        endif
    end subroutine eam7_hydrogen_embedding

    !> @brief Evaluate one symmetric reciprocal-exponential descriptor correction.
    pure subroutine symmetric_descriptor( amplitude, slope, center, value, energy, derivative )
        real(wp), intent(in)  :: amplitude, slope, center, value
        real(wp), intent(out) :: energy, derivative
        real(wp)              :: denominator, forward, reverse
        forward     = exp(slope*(value-center))
        reverse     = exp(slope*(center-value))
        denominator = forward+reverse
        energy      = amplitude/denominator
        derivative  = -amplitude*slope*(forward-reverse)/denominator**2
    end subroutine symmetric_descriptor

    !> @brief Evaluate the EAM6/EAM7 screened pair repulsion and radial derivative.
    pure subroutine pair_potential( model_version, first_type, second_type, distance, energy, derivative )
        use chempotpy_hni_eam67_parameters, only: nickel_pair_parameters
        integer, intent(in)   :: model_version, first_type, second_type
        real(wp), intent(in)  :: distance
        real(wp), intent(out) :: energy, derivative
        real(wp)              :: first_derivative, first_value, second_derivative, second_value
        call effective_charge(model_version, first_type, distance, first_value, first_derivative, nickel_pair_parameters)
        call effective_charge(model_version, second_type, distance, second_value, second_derivative, nickel_pair_parameters)
        energy = electrostatic_conversion*first_value*second_value/distance
        derivative = electrostatic_conversion*(first_value*second_derivative+ second_value*first_derivative- &
            & first_value*second_value/distance)/distance
    end subroutine pair_potential

    !> @brief Evaluate the Ni or version-specific H effective charge and radial derivative.
    pure subroutine effective_charge( model_version, atom_type, distance, value, derivative, nickel_parameters )
        integer, intent(in)   :: model_version, atom_type
        real(wp), intent(in)  :: distance, nickel_parameters(3)
        real(wp), intent(out) :: value, derivative
        real(wp)              :: cutoff_derivative, cutoff_value, exponential, forward, reverse, denominator
        if( atom_type == nickel_type )then
            call cutoff(distance, cutoff_value, cutoff_derivative)
            exponential = nickel_parameters(1)*exp(-nickel_parameters(3)*distance)
            value       = exponential*(1.0_wp+nickel_parameters(2)*distance)
            derivative  = (-nickel_parameters(3)*value+nickel_parameters(2)*exponential)*cutoff_value+ &
                & cutoff_derivative*value
            value = value*cutoff_value
        else
            value = (distance/0.93232_wp)**2.77511_wp*exp(-distance/0.37473_wp)+ (distance/0.54879_wp)**142.67761_wp* &
                & exp(-distance/0.01012_wp)+ (distance/0.37697_wp)**95.4552_wp*exp(-distance/0.01073_wp)
            derivative = exp(-distance/0.37473_wp)*(distance/0.93232_wp)**(2.77511_wp-1.0_wp)* (2.77511_wp/0.93232_wp- &
                & distance/(0.37473_wp*0.93232_wp))+ exp(-distance/0.01012_wp)*(distance/0.54879_wp)**(142.67761_wp- &
                & 1.0_wp)* (142.67761_wp/0.54879_wp-distance/(0.01012_wp*0.54879_wp))+ exp(-distance/0.01073_wp)* &
                & (distance/0.37697_wp)**(95.4552_wp-1.0_wp)* (95.4552_wp/0.37697_wp-distance/(0.01073_wp*0.37697_wp))
            forward     = exp(8.0_wp*(distance-2.4_wp))
            reverse     = exp(8.0_wp*(2.4_wp-distance))
            denominator = forward+reverse
            value       = value+0.015_wp/denominator
            derivative  = derivative-0.015_wp*8.0_wp*(forward-reverse)/denominator**2
            if( model_version == 7 .and. distance <= 1.09_wp )then
                value      = 1.0_wp-5.957_wp*distance**5+9.804_wp*distance**4-3.705_wp*distance**3
                derivative = -5.0_wp*5.957_wp*distance**4+4.0_wp*9.804_wp*distance**3- 3.0_wp*3.705_wp*distance**2
            endif
        endif
    end subroutine effective_charge

    !> @brief Evaluate the released smooth radial cutoff and derivative.
    pure subroutine cutoff( distance, value, derivative )
        use chempotpy_constants, only: PI
        real(wp), intent(in)  :: distance
        real(wp), intent(out) :: value, derivative
        real(wp)              :: cutoff_width, inner_offset, outer_offset
        if( distance <= inner_cutoff_angstrom )then
            value = 1.0_wp; derivative = 0.0_wp
        elseif( distance < outer_cutoff_angstrom )then
            cutoff_width = inner_cutoff_angstrom-outer_cutoff_angstrom
            inner_offset = distance-inner_cutoff_angstrom
            outer_offset = distance-outer_cutoff_angstrom
            value        = outer_offset/cutoff_width+sin(PI*(inner_offset+outer_offset)/cutoff_width)/(2.0_wp*PI)
            derivative   = (1.0_wp+cos(PI*(inner_offset+outer_offset)/cutoff_width))/cutoff_width
        else
            value = 0.0_wp; derivative = 0.0_wp
        endif
    end subroutine cutoff
end module chempotpy_hni_eam67_runtime
