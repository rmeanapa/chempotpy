!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the Pettifor-Ward aluminum potential through Python.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_kinds, only: real64
    use chempotpy_legacy_units, only: legacy_coordinates_to_bohr, &
        legacy_energy_to_electronvolt
    implicit none



    integer, parameter        :: state_count = 1
    real(real64), intent(in)  :: coordinates_angstrom(10000,3)
    integer, intent(in)       :: atom_count, derivative_order
    real(real64), intent(out) :: energy_ev(state_count)
    real(real64), intent(out) :: gradient_ev_per_angstrom( &
        state_count,atom_count,3 )
    real(real64), intent(out) :: coupling_per_angstrom( &
        state_count,state_count,atom_count,3 )
    real(real64) :: coordinates_bohr(atom_count,3), energy_hartree

    energy_ev                = 0.0d0
    gradient_ev_per_angstrom = 0.0d0
    coupling_per_angstrom    = 0.0d0
    call legacy_coordinates_to_bohr( &
        coordinates_angstrom(1:atom_count,:), coordinates_bohr )
    call evaluate_pettifor_ward( &
        coordinates_bohr, atom_count, energy_hartree )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*,*) 'Only energy is available'
    endif
end subroutine pes


!> @brief Evaluate the oscillatory Pettifor-Ward pair-potential energy.
subroutine evaluate_pettifor_ward( &
        coordinates_bohr, atom_count, energy_hartree )
    use chempotpy_geometry, only: cartesian_distance
    use chempotpy_kinds, only: real32, real64
    implicit none



    real(real64), parameter :: bohr_to_angstrom = 0.529177249d0
    real(real64), parameter :: effective_charge = 3.0d0
    real(real64), parameter :: amplitudes(3)    = [7.964d0, 1.275d0, 0.030d0]
    ! Preserve default-real rounding of the released phase multipliers.
    real(real64), parameter   :: phase_factors(3) = real( [-0.441_real32, 0.832_real32, 0.431_real32], real64 )
    real(real64), parameter   :: wave_factors(3)  = [0.156d0, 0.644d0, 0.958d0]
    real(real64), parameter   :: decay_factors(3) = [0.793d0, 0.698d0, 0.279d0]
    real(real64), intent(in)  :: coordinates_bohr(atom_count,3)
    integer, intent(in)       :: atom_count
    real(real64), intent(out) :: energy_hartree
    real(real64)              :: distance_bohr, pair_energy, pi, twice_fermi_wavevector
    real(real64)              :: phases(3), wavevectors(3), decay_rates(3)
    integer                   :: first_atom, second_atom

    pi                     = acos(-1.0d0)
    phases                 = phase_factors*pi
    twice_fermi_wavevector = 3.504d0*bohr_to_angstrom
    wavevectors            = wave_factors*twice_fermi_wavevector
    decay_rates            = decay_factors*twice_fermi_wavevector
    energy_hartree         = 0.0d0
    do first_atom = 1, atom_count-1
        do second_atom = first_atom+1, atom_count
            distance_bohr = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
            pair_energy = ((2.0d0*effective_charge**2)/distance_bohr) *(amplitudes(1)*cos(wavevectors(1)* &
                & distance_bohr+ phases(1)) *exp(-decay_rates(1)*distance_bohr) +amplitudes(2)*cos(wavevectors(2)* &
                & distance_bohr+ phases(2)) *exp(-decay_rates(2)*distance_bohr) +amplitudes(3)*cos(wavevectors(3)* &
                & distance_bohr+ phases(3)) *exp(-decay_rates(3)*distance_bohr))
            energy_hartree = energy_hartree+pair_energy
        enddo
    enddo
end subroutine evaluate_pettifor_ward
