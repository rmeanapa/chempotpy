!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_halp_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the Halicioglu-Pound potential through the Python contract.
    subroutine evaluate( coordinates_angstrom, atom_count, derivative_order, &
            energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
        use chempotpy_units, only: coordinates_to_bohr, &
            energy_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter    :: state_count = 1
        real(wp), intent(in)  :: coordinates_angstrom(10000,3)
        integer, intent(in)   :: atom_count
        integer, intent(in)   :: derivative_order
        real(wp), intent(out) :: energy_ev(state_count)
        real(wp), intent(out) :: gradient_ev_per_angstrom( &
            state_count,atom_count,3 )
        real(wp), intent(out) :: coupling_per_angstrom( &
            state_count,state_count,atom_count,3 )
        real(wp) :: coordinates_bohr(atom_count,3)
        real(wp) :: energy_hartree
        energy_ev                = 0.0e0_wp
        gradient_ev_per_angstrom = 0.0e0_wp
        coupling_per_angstrom    = 0.0e0_wp
        call coordinates_to_bohr( &
            coordinates_angstrom(1:atom_count,:), coordinates_bohr)
        call evaluate_halicioglu_pound( &
            coordinates_bohr, atom_count, energy_hartree)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate the pairwise Halicioglu-Pound aluminum cluster energy.
    subroutine evaluate_halicioglu_pound( &
            coordinates_bohr, atom_count, energy_hartree )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), parameter   :: hartree_to_ev       = 27.2113961e0_wp
        real(wp), parameter   :: bohr_to_angstrom    = 0.529177249e0_wp
        real(wp), parameter   :: sigma_angstrom      = 2.620e0_wp
        real(wp), parameter   :: well_depth_ev       = 0.3922e0_wp
        real(wp), parameter   :: repulsive_exponent  = 12.0e0_wp
        real(wp), parameter   :: attractive_exponent = 6.0e0_wp
        real(wp), intent(in)  :: coordinates_bohr(atom_count,3)
        integer, intent(in)   :: atom_count
        real(wp), intent(out) :: energy_hartree
        real(wp)              :: distance_bohr
        real(wp)              :: sigma_bohr, four_well_depth_hartree
        real(wp)              :: pair_sum
        integer               :: first_atom, second_atom
        sigma_bohr              = sigma_angstrom/bohr_to_angstrom
        four_well_depth_hartree = 4.0e0_wp*well_depth_ev/hartree_to_ev
        pair_sum                = 0.0e0_wp
        do first_atom = 1, atom_count-1
            do second_atom = first_atom+1, atom_count
                distance_bohr = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
                pair_sum = pair_sum +(sigma_bohr/distance_bohr)**repulsive_exponent -(sigma_bohr/distance_bohr)** &
                    & attractive_exponent
            enddo
        enddo
        energy_hartree = four_well_depth_hartree*pair_sum
    end subroutine evaluate_halicioglu_pound
end module chempotpy_aln_halp_pes

!> @brief Expose Aln_HalP through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_halp_pes, only: evaluate
    use chempotpy_kinds,        only: wp
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
