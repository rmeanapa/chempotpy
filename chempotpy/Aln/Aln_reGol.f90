!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_regol_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the reparameterized Gollisch potential through Python.
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
        call evaluate_reparameterized_gollisch( &
            coordinates_bohr, atom_count, energy_hartree)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate the reparameterized Gollisch embedded-atom energy.
    subroutine evaluate_reparameterized_gollisch( &
            coordinates_bohr, atom_count, energy_hartree )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), parameter   :: hartree_to_ev      = 27.2113961e0_wp
        real(wp), parameter   :: bohr_to_angstrom   = 0.529177249e0_wp
        real(wp), parameter   :: repulsive_scale    = 641.335e0_wp/hartree_to_ev
        real(wp), parameter   :: embedding_scale    = 12.2192e0_wp/hartree_to_ev
        real(wp), parameter   :: repulsive_decay    = 2.9913e0_wp*bohr_to_angstrom
        real(wp), parameter   :: density_decay      = 1.364e0_wp*bohr_to_angstrom
        real(wp), parameter   :: embedding_exponent = 0.66258e0_wp
        real(wp), intent(in)  :: coordinates_bohr(atom_count,3)
        integer, intent(in)   :: atom_count
        real(wp), intent(out) :: energy_hartree
        real(wp)              :: distance_bohr
        real(wp)              :: repulsive_sum, density_sum
        integer               :: center_atom, neighbor_atom
        energy_hartree = 0.0e0_wp
        do center_atom = 1, atom_count
            repulsive_sum = 0.0e0_wp
            density_sum   = 0.0e0_wp
            do neighbor_atom = 1, atom_count
                if( center_atom /= neighbor_atom )then
                    distance_bohr = cartesian_distance( coordinates_bohr(center_atom,:), &
                        & coordinates_bohr(neighbor_atom,:) )
                    repulsive_sum = repulsive_sum +exp(-distance_bohr*repulsive_decay)
                    density_sum   = density_sum +exp(-distance_bohr*density_decay)
                endif
            enddo
            energy_hartree = energy_hartree +repulsive_scale*repulsive_sum -embedding_scale*density_sum** &
                & embedding_exponent
        enddo
    end subroutine evaluate_reparameterized_gollisch
end module chempotpy_aln_regol_pes

!> @brief Expose Aln_reGol through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_regol_pes, only: evaluate
    use chempotpy_kinds,         only: wp
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
