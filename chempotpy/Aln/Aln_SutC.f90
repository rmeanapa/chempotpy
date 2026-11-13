!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_sutc_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the Sutton-Chen aluminum potential through Python.
    subroutine evaluate( coordinates_angstrom, atom_count, derivative_order, &
            energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
        use chempotpy_kinds, only: wp
        use chempotpy_units, only: coordinates_to_bohr, &
            energy_to_electronvolt
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
        call evaluate_sutton_chen( &
            coordinates_bohr, atom_count, energy_hartree)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate the Sutton-Chen embedded-atom aluminum cluster energy.
    subroutine evaluate_sutton_chen( &
            coordinates_bohr, atom_count, energy_hartree )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), parameter :: hartree_to_ev     = 27.2113961e0_wp
        real(wp), parameter :: bohr_to_angstrom  = 0.529177249e0_wp
        real(wp), parameter :: pair_scale        = 0.0331467e0_wp/hartree_to_ev
        real(wp), parameter :: embedding_scale   = 0.54357e0_wp/hartree_to_ev
        real(wp), parameter :: length_scale_bohr = &
            4.05e0_wp/bohr_to_angstrom
        real(wp), parameter   :: pair_exponent    = 7.0e0_wp
        real(wp), parameter   :: density_exponent = 6.0e0_wp
        real(wp), intent(in)  :: coordinates_bohr(atom_count,3)
        integer, intent(in)   :: atom_count
        real(wp), intent(out) :: energy_hartree
        real(wp)              :: distance_bohr, pair_sum, density_sum, embedding_sum
        integer               :: center_atom, neighbor_atom, first_atom, second_atom
        pair_sum = 0.0e0_wp
        do first_atom = 1, atom_count-1
            do second_atom = first_atom+1, atom_count
                distance_bohr = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
                pair_sum      = pair_sum +(length_scale_bohr/distance_bohr)**pair_exponent
            enddo
        enddo
        embedding_sum = 0.0e0_wp
        do center_atom = 1, atom_count
            density_sum = 0.0e0_wp
            do neighbor_atom = 1, atom_count
                if( neighbor_atom == center_atom ) cycle
                distance_bohr = cartesian_distance( coordinates_bohr(center_atom,:), coordinates_bohr(neighbor_atom,:) )
                density_sum   = density_sum +(length_scale_bohr/distance_bohr)**density_exponent
            enddo
            embedding_sum = embedding_sum+sqrt(density_sum)
        enddo
        energy_hartree = pair_scale*pair_sum-embedding_scale*embedding_sum
    end subroutine evaluate_sutton_chen
end module chempotpy_aln_sutc_pes

!> @brief Expose Aln_SutC through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_sutc_pes, only: evaluate
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
