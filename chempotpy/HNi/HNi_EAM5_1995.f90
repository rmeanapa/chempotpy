! System:                        HNi
! Functional form:               embedded-atom method
! Common name:                   EAM5 H on Ni(111)
! Reference:                     S. E. Wonchoba, W.-P. Hu, and D. G. Truhlar, Phys. Rev. B 51, 9985 (1995)
! Number of derivatives:         1
! Number of bodies:              1 active H plus a fixed 698-atom Ni(111) environment
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable POTLIB potfile1 parameters are embedded
! Notes:                         The public coordinate is the active H atom; the canonical Ni(111) slab is fixed.

!> @brief Evaluate the EAM5 interaction of one active H atom with the canonical Ni(111) environment.
module chempotpy_hni_eam5_1995_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
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

    !> @brief Return EAM5 energy and analytic H-atom gradient in the CHEMPOTPY eV/angstrom contract.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_hni_eam5_parameters,  only: canonical_asymptote_correction_electronvolt, energy_zero_hartree
        use chempotpy_hni_ni111_parameters, only: fixed_ni_coordinates_bohr
        use chempotpy_status,               only: chempotpy_status_t, status_derivative_unavailable
        use chempotpy_units,                only: bohr_to_angstrom, hartree_to_electronvolt
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms = 1, nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp), allocatable :: coordinates(:,:), embedding_derivatives(:), densities(:)
        real(wp)              :: energy_electronvolt, gradient(3)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad == 0 .or. igrad == 1 )then
            allocate(coordinates(1+size(fixed_ni_coordinates_bohr, 1),3), &
                & embedding_derivatives(1+size(fixed_ni_coordinates_bohr, 1)), &
                & densities(1+size(fixed_ni_coordinates_bohr, 1)))
            coordinates(1,:)  = x(1,:)
            coordinates(2:,:) = fixed_ni_coordinates_bohr*bohr_to_angstrom
            call evaluate_eam5(coordinates, energy_electronvolt, gradient, densities, embedding_derivatives)
            p(1) = energy_electronvolt+energy_zero_hartree*hartree_to_electronvolt+ &
                & canonical_asymptote_correction_electronvolt
            if( igrad == 1 ) g(1,1,:) = gradient
        elseif( present(status) )then
            call status%set(status_derivative_unavailable, 'HNi_EAM5_1995 provides energy and first derivatives')
        endif
    end subroutine evaluate

    !> @brief Contract the EAM5 density, embedding, and screened pair terms for H plus fixed Ni atoms.
    pure subroutine evaluate_eam5( coordinates, energy, hydrogen_gradient, densities, embedding_derivatives )
        real(wp), intent(in)  :: coordinates(:,:)
        real(wp), intent(out) :: energy, hydrogen_gradient(3), densities(size(coordinates, 1))
        real(wp), intent(out) :: embedding_derivatives(size(coordinates, 1))
        real(wp)              :: density_derivative, density_value, displacement(3), distance
        real(wp)              :: pair_derivative, pair_energy, radial_derivative
        integer               :: atom, first_atom, first_type, second_atom, second_type
        densities = 0.0_wp
        energy    = 0.0_wp
        do first_atom=1, size(coordinates, 1)-1
            first_type = merge(hydrogen_type, nickel_type, first_atom == 1)
            do second_atom=first_atom+1, size(coordinates, 1)
                second_type  = nickel_type
                displacement = coordinates(first_atom,:)-coordinates(second_atom,:)
                distance     = sqrt(sum(displacement*displacement))
                if( distance >= outer_cutoff_angstrom ) cycle
                call atomic_density(second_type, distance, coordinates(second_atom,3), density_value, density_derivative)
                densities(first_atom) = densities(first_atom)+density_value
                call atomic_density(first_type, distance, coordinates(first_atom,3), density_value, density_derivative)
                densities(second_atom) = densities(second_atom)+density_value
                if( first_atom == 1 )then
                    call pair_potential(first_type, second_type, distance, pair_energy, pair_derivative)
                    energy = energy+pair_energy
                endif
            enddo
        enddo
        do atom=1, size(coordinates, 1)
            if( atom == 1 )then
                call hydrogen_embedding(densities(atom), pair_energy, embedding_derivatives(atom))
            else
                call nickel_embedding(densities(atom), pair_energy, embedding_derivatives(atom))
            endif
            energy = energy+pair_energy
        enddo
        hydrogen_gradient = 0.0_wp
        do second_atom=2, size(coordinates, 1)
            displacement = coordinates(1,:)-coordinates(second_atom,:)
            distance     = sqrt(sum(displacement*displacement))
            if( distance >= outer_cutoff_angstrom ) cycle
            call atomic_density(nickel_type, distance, coordinates(second_atom,3), density_value, density_derivative)
            radial_derivative = embedding_derivatives(1)*density_derivative
            call atomic_density(hydrogen_type, distance, coordinates(1,3), density_value, density_derivative)
            radial_derivative = radial_derivative+embedding_derivatives(second_atom)*density_derivative
            call pair_potential(hydrogen_type, nickel_type, distance, pair_energy, pair_derivative)
            radial_derivative = radial_derivative+pair_derivative
            hydrogen_gradient = hydrogen_gradient+radial_derivative*displacement/distance
        enddo
    end subroutine evaluate_eam5

    !> @brief Evaluate the released atomic-density contribution and radial derivative.
    pure subroutine atomic_density( atom_type, distance, source_z, density, derivative )
        use chempotpy_constants, only: PI
        integer, intent(in)   :: atom_type
        real(wp), intent(in)  :: distance, source_z
        real(wp), intent(out) :: density, derivative
        real(wp)              :: cutoff_derivative, cutoff_value, d_sum_d, d_sum_s
        real(wp)              :: d_component, electron_count_d, electron_count_s, sum_d, sum_s, term
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
        d_component = (electron_count_s*sum_s*d_sum_s+electron_count_d*sum_d*d_sum_d)/(2.0_wp*PI)
        density     = (electron_count_s*sum_s**2+electron_count_d*sum_d**2)/(4.0_wp*PI)
        derivative  = d_component*cutoff_value+cutoff_derivative*density
        density     = density*cutoff_value
    end subroutine atomic_density

    !> @brief Evaluate the EAM5 nickel embedding function and its density derivative.
    pure subroutine nickel_embedding( density, energy, derivative )
        use chempotpy_hni_eam5_parameters, only: ni_embedding_parameters
        real(wp), intent(in)  :: density
        real(wp), intent(out) :: energy, derivative
        real(wp)              :: delta, delta_squared, scaled_density, term
        if( density > ni_embedding_parameters(8) )then
            energy     = ni_embedding_parameters(12)
            derivative = 0.0_wp
        elseif( density > ni_embedding_parameters(7) )then
            delta         = density-ni_embedding_parameters(8)
            delta_squared = delta*delta
            energy = ni_embedding_parameters(12)+delta_squared*delta*(ni_embedding_parameters(11)+ delta* &
                & (ni_embedding_parameters(10)+delta*ni_embedding_parameters(9)))
            derivative = delta_squared*(3.0_wp*ni_embedding_parameters(11)+delta*(4.0_wp* ni_embedding_parameters(10)+ &
                & delta*5.0_wp*ni_embedding_parameters(9)))
        else
            scaled_density = -ni_embedding_parameters(2)*density
            term           = ni_embedding_parameters(1)*exp(scaled_density)
            energy         = term*density
            derivative     = term*(1.0_wp+scaled_density)
            scaled_density = -ni_embedding_parameters(4)*density
            term           = ni_embedding_parameters(3)*density
            energy         = energy+term**3*exp(scaled_density)
            derivative     = derivative+term**2*exp(scaled_density)*(3.0_wp*ni_embedding_parameters(3)- &
                & ni_embedding_parameters(4)*term)
            scaled_density = -ni_embedding_parameters(6)*density
            term           = ni_embedding_parameters(5)*exp(scaled_density)
            energy         = energy+term*density
            derivative     = derivative+term*(1.0_wp+scaled_density)
        endif
    end subroutine nickel_embedding

    !> @brief Evaluate the EAM5 hydrogen embedding function and its density derivative.
    pure subroutine hydrogen_embedding( density, energy, derivative )
        use chempotpy_hni_eam5_parameters, only: h_embedding_parameters
        real(wp), intent(in)  :: density
        real(wp), intent(out) :: energy, derivative
        real(wp)              :: scaled_density, term
        scaled_density = -h_embedding_parameters(2)*density
        term           = h_embedding_parameters(1)*exp(scaled_density)
        energy         = term*density
        derivative     = term*(1.0_wp+scaled_density)
    end subroutine hydrogen_embedding

    !> @brief Evaluate the screened EAM5 pair repulsion and its radial derivative.
    pure subroutine pair_potential( first_type, second_type, distance, energy, derivative )
        use chempotpy_hni_eam5_parameters, only: h_pair_parameters, ni_pair_parameters
        integer, intent(in)   :: first_type, second_type
        real(wp), intent(in)  :: distance
        real(wp), intent(out) :: energy, derivative
        real(wp)              :: first_derivative, first_value, second_derivative, second_value
        if( first_type == nickel_type )then
            call screened_charge(ni_pair_parameters, distance, first_value, first_derivative)
        else
            call screened_charge(h_pair_parameters, distance, first_value, first_derivative)
        endif
        if( second_type == nickel_type )then
            call screened_charge(ni_pair_parameters, distance, second_value, second_derivative)
        else
            call screened_charge(h_pair_parameters, distance, second_value, second_derivative)
        endif
        energy = electrostatic_conversion*first_value*second_value/distance
        derivative = electrostatic_conversion*(first_value*second_derivative+ second_value*first_derivative- &
            & first_value*second_value/distance)/distance
    end subroutine pair_potential

    !> @brief Evaluate one atom type's switched effective charge and radial derivative.
    pure subroutine screened_charge( parameters, distance, value, derivative )
        real(wp), intent(in)  :: parameters(3), distance
        real(wp), intent(out) :: value, derivative
        real(wp)              :: cutoff_derivative, cutoff_value, exponential
        call cutoff(distance, cutoff_value, cutoff_derivative)
        exponential = parameters(1)*exp(-parameters(3)*distance)
        value       = exponential*(1.0_wp+parameters(2)*distance)
        derivative  = (-parameters(3)*value+parameters(2)*exponential)*cutoff_value+cutoff_derivative*value
        value       = value*cutoff_value
    end subroutine screened_charge

    !> @brief Evaluate the released smooth radial cutoff and derivative.
    pure subroutine cutoff( distance, value, derivative )
        use chempotpy_constants, only: PI
        real(wp), intent(in)  :: distance
        real(wp), intent(out) :: value, derivative
        real(wp)              :: cutoff_width, inner_offset, outer_offset
        if( distance <= inner_cutoff_angstrom )then
            value      = 1.0_wp
            derivative = 0.0_wp
        elseif( distance < outer_cutoff_angstrom )then
            cutoff_width = inner_cutoff_angstrom-outer_cutoff_angstrom
            inner_offset = distance-inner_cutoff_angstrom
            outer_offset = distance-outer_cutoff_angstrom
            value        = outer_offset/cutoff_width+sin(PI*(inner_offset+outer_offset)/cutoff_width)/(2.0_wp*PI)
            derivative   = (1.0_wp+cos(PI*(inner_offset+outer_offset)/cutoff_width))/cutoff_width
        else
            value      = 0.0_wp
            derivative = 0.0_wp
        endif
    end subroutine cutoff
end module chempotpy_hni_eam5_1995_pes

!> @brief Expose HNi_EAM5_1995 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_hni_eam5_1995_pes, only: evaluate
    use chempotpy_kinds,             only: wp
    implicit none
    integer, parameter    :: natoms = 1, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
