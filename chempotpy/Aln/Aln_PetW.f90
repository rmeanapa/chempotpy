!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_petw_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the Pettifor-Ward aluminum potential through Python.
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
        call evaluate_pettifor_ward( &
            coordinates_bohr, atom_count, energy_hartree)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate the oscillatory Pettifor-Ward pair-potential energy.
    subroutine evaluate_pettifor_ward( &
            coordinates_bohr, atom_count, energy_hartree )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: real32, wp
        implicit none
        real(wp), parameter :: bohr_to_angstrom = 0.529177249e0_wp
        real(wp), parameter :: effective_charge = 3.0e0_wp
        real(wp), parameter :: amplitudes(3)    = [7.964e0_wp, 1.275e0_wp, 0.030e0_wp]
        ! Preserve default-real rounding of the released phase multipliers.
        real(wp), parameter   :: phase_factors(3) = real( [-0.441_real32, 0.832_real32, 0.431_real32], wp )
        real(wp), parameter   :: wave_factors(3)  = [0.156e0_wp, 0.644e0_wp, 0.958e0_wp]
        real(wp), parameter   :: decay_factors(3) = [0.793e0_wp, 0.698e0_wp, 0.279e0_wp]
        real(wp), intent(in)  :: coordinates_bohr(atom_count,3)
        integer, intent(in)   :: atom_count
        real(wp), intent(out) :: energy_hartree
        real(wp)              :: distance_bohr, pair_energy, pi, twice_fermi_wavevector
        real(wp)              :: phases(3), wavevectors(3), decay_rates(3)
        integer               :: first_atom, second_atom
        pi                     = acos(-1.0e0_wp)
        phases                 = phase_factors*pi
        twice_fermi_wavevector = 3.504e0_wp*bohr_to_angstrom
        wavevectors            = wave_factors*twice_fermi_wavevector
        decay_rates            = decay_factors*twice_fermi_wavevector
        energy_hartree         = 0.0e0_wp
        do first_atom = 1, atom_count-1
            do second_atom = first_atom+1, atom_count
                distance_bohr = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
                pair_energy = ((2.0e0_wp*effective_charge**2)/distance_bohr) *(amplitudes(1)*cos(wavevectors(1)* &
                    & distance_bohr+ phases(1)) *exp(-decay_rates(1)*distance_bohr) +amplitudes(2)*cos(wavevectors(2)* &
                    & distance_bohr+ phases(2)) *exp(-decay_rates(2)*distance_bohr) +amplitudes(3)*cos(wavevectors(3)* &
                    & distance_bohr+ phases(3)) *exp(-decay_rates(3)*distance_bohr))
                energy_hartree = energy_hartree+pair_energy
            enddo
        enddo
    end subroutine evaluate_pettifor_ward
end module chempotpy_aln_petw_pes

!> @brief Expose Aln_PetW through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_petw_pes, only: evaluate
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
