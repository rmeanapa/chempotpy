! System:                        H3ClO
! Functional form:               permutationally invariant polynomial
! Common name:                   HCl-H2O Bowman-Mancini PIP
! Reference:                     J. S. Mancini and J. M. Bowman, J. Chem. Phys. 138, 121102 (2013)
! Number of derivatives:         0
! Number of bodies:              5
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; released PCF coefficients are embedded
! Notes:                         Canonical atom order is H, H, H, O, Cl.

!> @brief Evaluate the full-dimensional Bowman-Mancini HCl-H2O PIP surface.
module chempotpy_h3clo_pip_2013_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the HCl-H2O energy in electronvolts.
    subroutine evaluate( x, derivative_order, energy, gradient, coupling, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable
        type(chempotpy_status_t), intent(out), optional :: status
        integer, intent(in)   :: derivative_order
        real(wp), intent(in)  :: x(5,3)
        real(wp), intent(out) :: energy(1), gradient(1,5,3), coupling(1,1,5,3)
        if( present(status) ) call status%clear()
        gradient = 0.0_wp
        coupling = 0.0_wp
        call evaluate_energy(x, energy(1))
        if( derivative_order > 0 .and. present(status) ) call status%set(status_derivative_unavailable, &
            & 'the released HCl-H2O surface provides energy only')
    end subroutine evaluate

    !> @brief Evaluate the released many-body expansion and convert hartree to electronvolts.
    pure subroutine evaluate_energy( x, energy )
        use chempotpy_geometry,               only: cartesian_distance
        use chempotpy_hcl_h2o_pip_basis,      only: evaluate_hcl_h2o_basis
        use chempotpy_hcl_h2o_pip_parameters, only: five_body_coefficients, five_body_cutoff_radius, &
            & hcl_coefficients, hh_coefficients, ho_coefficients, ocl_coefficients, pair_morse_length
        use chempotpy_units, only: angstrom_to_bohr, hartree_to_electronvolt
        real(wp), intent(in)  :: x(5,3)
        real(wp), intent(out) :: energy
        real(wp)              :: basis(0:1602), cutoff, distance_matrix(0:4,0:4), rms_distance
        real(wp)              :: transformed_distances(0:4,0:4), value_hartree
        integer               :: first_atom, second_atom
        distance_matrix = 0.0_wp
        do first_atom=0, 3
            do second_atom=first_atom+1, 4
                distance_matrix(first_atom,second_atom) = cartesian_distance(x(first_atom+1,:), x(second_atom+1,:))* &
                    & angstrom_to_bohr
                distance_matrix(second_atom,first_atom) = distance_matrix(first_atom,second_atom)
            enddo
        enddo
        value_hartree = 0.0_wp
        do first_atom=0, 1
            do second_atom=first_atom+1, 2
                value_hartree = value_hartree+pair_energy(distance_matrix(first_atom,second_atom), hh_coefficients)
            enddo
        enddo
        do first_atom=0, 2
            value_hartree = value_hartree+pair_energy(distance_matrix(first_atom,3), ho_coefficients) + &
                & pair_energy(distance_matrix(first_atom,4), hcl_coefficients)
        enddo
        value_hartree         = value_hartree+pair_energy(distance_matrix(3,4), ocl_coefficients)
        transformed_distances = exp(-distance_matrix/pair_morse_length)
        do first_atom=0, 4
            transformed_distances(first_atom,first_atom) = 0.0_wp
        enddo
        call evaluate_hcl_h2o_basis(transformed_distances, basis)
        rms_distance = 0.0_wp
        do first_atom=0, 3
            rms_distance = rms_distance+sum(distance_matrix(first_atom,first_atom+1:4)**2)
        enddo
        rms_distance = sqrt(rms_distance/10.0_wp)
        if( rms_distance < five_body_cutoff_radius )then
            cutoff        = (1.0_wp-rms_distance/five_body_cutoff_radius)**3
            value_hartree = value_hartree+cutoff*dot_product(five_body_coefficients, basis)
        endif
        energy = value_hartree*hartree_to_electronvolt
    end subroutine evaluate_energy

    !> @brief Evaluate one released pair contribution in hartree.
    pure function pair_energy( distance, coefficients ) result(value)
        use chempotpy_hcl_h2o_pip_parameters, only: pair_cutoff_radius, pair_morse_length
        real(wp), intent(in) :: distance, coefficients(3)
        real(wp)             :: transformed_distance, value
        if( distance >= pair_cutoff_radius )then
            value = 0.0_wp
        else
            transformed_distance = exp(-distance/pair_morse_length)
            value = (1.0_wp-distance/pair_cutoff_radius)**3/distance * (coefficients(1)+transformed_distance* &
                & (coefficients(2)+transformed_distance*coefficients(3)))
        endif
    end function pair_energy
end module chempotpy_h3clo_pip_2013_pes

!> @brief Expose H3ClO_PIP_2013 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3clo_pip_2013_pes, only: evaluate
    use chempotpy_kinds,              only: wp
    implicit none
    integer, parameter    :: natoms = 5, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
