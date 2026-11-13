!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_rehalp_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the reparameterized Halicioglu-Pound potential through Python.
    subroutine evaluate( coordinates_angstrom, atom_count, derivative_order, &
            energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
        use chempotpy_units, only: coordinates_to_bohr, &
            energy_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter    :: state_count = 1
        real(wp), intent(in)  :: coordinates_angstrom(10000,3)
        integer, intent(in)   :: atom_count, derivative_order
        real(wp), intent(out) :: energy_ev(state_count)
        real(wp), intent(out) :: gradient_ev_per_angstrom( &
            state_count,atom_count,3 )
        real(wp), intent(out) :: coupling_per_angstrom( &
            state_count,state_count,atom_count,3 )
        real(wp) :: coordinates_bohr(atom_count,3), energy_hartree
        energy_ev                = 0.0e0_wp
        gradient_ev_per_angstrom = 0.0e0_wp
        coupling_per_angstrom    = 0.0e0_wp
        call coordinates_to_bohr( &
            coordinates_angstrom(1:atom_count,:), coordinates_bohr)
        call evaluate_reparameterized_halicioglu_pound( &
            coordinates_bohr, atom_count, energy_hartree)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate the generalized Lennard-Jones aluminum cluster energy.
    subroutine evaluate_reparameterized_halicioglu_pound( &
            coordinates_bohr, atom_count, energy_hartree )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), parameter :: hartree_to_ev    = 27.2113961e0_wp
        real(wp), parameter :: bohr_to_angstrom = 0.529177249e0_wp
        ! Preserve the default-real parameter rounding in the released kernel.
        real(wp), parameter :: sigma_bohr = &
            2.5383489131927490234e0_wp/bohr_to_angstrom
        real(wp), parameter :: well_scale_hartree = &
            0.23957009613513946533e0_wp/hartree_to_ev
        real(wp), parameter :: attractive_exponent = &
            5.8041033744812011719e0_wp
        real(wp), parameter :: repulsive_exponent = &
            8.1436252593994140625e0_wp
        real(wp), intent(in)  :: coordinates_bohr(atom_count,3)
        integer, intent(in)   :: atom_count
        real(wp), intent(out) :: energy_hartree
        real(wp)              :: distance_bohr
        real(wp)              :: repulsive_term, attractive_term, pair_sum
        integer               :: first_atom, second_atom
        pair_sum = 0.0e0_wp
        do first_atom = 1, atom_count-1
            do second_atom = first_atom+1, atom_count
                distance_bohr   = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
                repulsive_term  = (sigma_bohr/distance_bohr) **repulsive_exponent
                attractive_term = (sigma_bohr/distance_bohr) **attractive_exponent
                pair_sum        = pair_sum+repulsive_term-attractive_term
            enddo
        enddo
        energy_hartree = well_scale_hartree*pair_sum
    end subroutine evaluate_reparameterized_halicioglu_pound
end module chempotpy_aln_rehalp_pes

!> @brief Expose Aln_reHalP through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_rehalp_pes, only: evaluate
    use chempotpy_kinds,          only: wp
    implicit none
    integer, parameter    :: state_count = 1
    real(wp), intent(in)  :: coordinates_angstrom(10000,3)
    integer, intent(in)   :: atom_count, derivative_order
    real(wp), intent(out) :: energy_ev(state_count)
    real(wp), intent(out) :: gradient_ev_per_angstrom(state_count,atom_count,3)
    real(wp), intent(out) :: coupling_per_angstrom(state_count,state_count,atom_count,3)
    call evaluate(coordinates_angstrom, atom_count, derivative_order, energy_ev, gradient_ev_per_angstrom, &
        & coupling_per_angstrom)
end subroutine pes
