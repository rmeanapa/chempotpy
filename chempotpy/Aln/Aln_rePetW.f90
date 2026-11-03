!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the reparameterized Pettifor-Ward potential through Python.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_kinds,        only: real64
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
    call evaluate_reparameterized_pettifor_ward( &
        coordinates_bohr, atom_count, energy_hartree )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*, *) 'Only energy is available'
    end if
end subroutine pes

!> @brief Evaluate the reparameterized oscillatory Pettifor-Ward energy.
subroutine evaluate_reparameterized_pettifor_ward( &
        coordinates_bohr, atom_count, energy_hartree )
    use chempotpy_geometry, only: cartesian_distance
    use chempotpy_kinds,    only: real32, real64
    implicit none
    ! Preserve default-real rounding of every released fit parameter.
    real(real64), parameter :: amplitudes(3) = real([ 67.1568148510014709_real32, 0.164142647777234973_real32, &
        & 57.3815339521250607_real32], real64)
    real(real64), parameter :: phases(3) = real([ 47.5280898876404478_real32, 96.8734733756717077_real32, &
        & 75.0854909623839717_real32], real64)
    real(real64), parameter :: wavevectors(3) = real([ 1.32467024914509035_real32, 1.23111871030776743_real32, &
        & 0.274157303370786531_real32], real64)
    real(real64), parameter :: decay_rates(3) = real([ 3.50244259892525633_real32, 0.513678553981436248_real32, &
        & 3.57093307278944794_real32], real64)
    real(real64), intent(in)  :: coordinates_bohr(atom_count,3)
    integer, intent(in)       :: atom_count
    real(real64), intent(out) :: energy_hartree
    real(real64)              :: distance_bohr, pair_energy
    integer                   :: first_atom, second_atom

    energy_hartree = 0.0d0
    do first_atom = 1, atom_count-1
        do second_atom = first_atom+1, atom_count
            distance_bohr = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
            pair_energy = (18.0d0/distance_bohr) *(amplitudes(1)*cos(wavevectors(1)*distance_bohr+phases(1)) *exp(- &
                & decay_rates(1)*distance_bohr) +amplitudes(2)*cos(wavevectors(2)*distance_bohr+phases(2)) *exp(- &
                & decay_rates(2)*distance_bohr) +amplitudes(3)*cos(wavevectors(3)*distance_bohr+phases(3)) *exp(- &
                & decay_rates(3)*distance_bohr))
            energy_hartree = energy_hartree+pair_energy
        end do
    end do
end subroutine evaluate_reparameterized_pettifor_ward
