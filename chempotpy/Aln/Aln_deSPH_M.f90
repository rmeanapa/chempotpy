!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_desph_m_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the deSPH/M aluminum potential through Python.
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
        call evaluate_desph_morse( &
            coordinates_bohr, atom_count, energy_hartree)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate the distance-dependent de Sainte Claire Morse energy.
    subroutine evaluate_desph_morse( &
            coordinates_bohr, atom_count, energy_hartree )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), parameter :: bohr_to_angstrom          = 0.529177249e0_wp
        real(wp), parameter :: hartree_to_kcal_per_mol   = 627.5095e0_wp
        real(wp), parameter :: equilibrium_distance_bohr = &
            2.834e0_wp/bohr_to_angstrom
        real(wp), parameter :: well_depth_hartree = &
            12.933e0_wp/hartree_to_kcal_per_mol
        real(wp), parameter :: zeroth_decay    = 0.5282e0_wp*bohr_to_angstrom
        real(wp), parameter :: quadratic_decay = &
            4.390e-3_wp*bohr_to_angstrom**3
        real(wp), parameter   :: cubic_decay = 5.109e0_wp*bohr_to_angstrom**4
        real(wp), intent(in)  :: coordinates_bohr(atom_count,3)
        integer, intent(in)   :: atom_count
        real(wp), intent(out) :: energy_hartree
        real(wp)              :: distance_bohr
        real(wp)              :: radial_offset, decay, pair_energy
        integer               :: first_atom, second_atom
        energy_hartree = 0.0e0_wp
        do first_atom = 1, atom_count-1
            do second_atom = first_atom+1, atom_count
                distance_bohr  = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
                radial_offset  = distance_bohr-equilibrium_distance_bohr
                decay          = zeroth_decay+quadratic_decay*radial_offset**2 +cubic_decay*radial_offset**3
                pair_energy    = well_depth_hartree *(1.0e0_wp-exp(-decay*radial_offset))**2-well_depth_hartree
                energy_hartree = energy_hartree+pair_energy
            enddo
        enddo
    end subroutine evaluate_desph_morse
end module chempotpy_aln_desph_m_pes

!> @brief Expose Aln_deSPH_M through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_desph_m_pes, only: evaluate
    use chempotpy_kinds,           only: wp
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
