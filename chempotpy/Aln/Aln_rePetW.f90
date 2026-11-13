!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_repetw_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the reparameterized Pettifor-Ward potential through Python.
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
        call evaluate_reparameterized_pettifor_ward( &
            coordinates_bohr, atom_count, energy_hartree)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate the reparameterized oscillatory Pettifor-Ward energy.
    subroutine evaluate_reparameterized_pettifor_ward( &
            coordinates_bohr, atom_count, energy_hartree )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: real32, wp
        implicit none
        ! Preserve default-real rounding of every released fit parameter.
        real(wp), parameter :: amplitudes(3) = real([ 67.1568148510014709_real32, 0.164142647777234973_real32, &
            & 57.3815339521250607_real32], wp)
        real(wp), parameter :: phases(3) = real([ 47.5280898876404478_real32, 96.8734733756717077_real32, &
            & 75.0854909623839717_real32], wp)
        real(wp), parameter :: wavevectors(3) = real([ 1.32467024914509035_real32, 1.23111871030776743_real32, &
            & 0.274157303370786531_real32], wp)
        real(wp), parameter :: decay_rates(3) = real([ 3.50244259892525633_real32, 0.513678553981436248_real32, &
            & 3.57093307278944794_real32], wp)
        real(wp), intent(in)  :: coordinates_bohr(atom_count,3)
        integer, intent(in)   :: atom_count
        real(wp), intent(out) :: energy_hartree
        real(wp)              :: distance_bohr, pair_energy
        integer               :: first_atom, second_atom
        energy_hartree = 0.0e0_wp
        do first_atom = 1, atom_count-1
            do second_atom = first_atom+1, atom_count
                distance_bohr = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
                pair_energy = (18.0e0_wp/distance_bohr) *(amplitudes(1)*cos(wavevectors(1)*distance_bohr+phases(1)) * &
                    & exp(- decay_rates(1)*distance_bohr) +amplitudes(2)*cos(wavevectors(2)*distance_bohr+phases(2)) * &
                    & exp(- decay_rates(2)*distance_bohr) +amplitudes(3)*cos(wavevectors(3)*distance_bohr+phases(3)) * &
                    & exp(- decay_rates(3)*distance_bohr))
                energy_hartree = energy_hartree+pair_energy
            enddo
        enddo
    end subroutine evaluate_reparameterized_pettifor_ward
end module chempotpy_aln_repetw_pes

!> @brief Expose Aln_rePetW through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_repetw_pes, only: evaluate
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
