!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the reparameterized Sutton-Chen potential through Python.
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
    call evaluate_reparameterized_sutton_chen( &
        coordinates_bohr, atom_count, energy_hartree )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*, *) 'Only energy is available'
    end if
end subroutine pes

!> @brief Evaluate the reparameterized Sutton-Chen cluster energy.
subroutine evaluate_reparameterized_sutton_chen( &
        coordinates_bohr, atom_count, energy_hartree )
    use chempotpy_geometry, only: cartesian_distance
    use chempotpy_kinds,    only: real32, real64
    implicit none
    real(real64), parameter :: hartree_to_ev    = 27.2113961d0
    real(real64), parameter :: bohr_to_angstrom = 0.529177249d0
    ! Preserve default-real rounding of the released fit parameters.
    real(real64), parameter :: pair_scale = &
        real(0.0230405_real32, real64)/hartree_to_ev
    real(real64), parameter :: embedding_scale = &
        real(0.56485_real32, real64)/hartree_to_ev
    real(real64), parameter :: length_scale_bohr = &
        4.05d0/bohr_to_angstrom
    real(real64), parameter :: pair_exponent = &
        real(6.8529_real32, real64)
    real(real64), parameter :: density_exponent = &
        real(4.2056_real32, real64)
    real(real64), intent(in)  :: coordinates_bohr(atom_count,3)
    integer, intent(in)       :: atom_count
    real(real64), intent(out) :: energy_hartree
    real(real64)              :: distance_bohr, pair_sum, density_sum, embedding_sum
    integer                   :: center_atom, neighbor_atom, first_atom, second_atom
    pair_sum = 0.0d0
    do first_atom = 1, atom_count-1
        do second_atom = first_atom+1, atom_count
            distance_bohr = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
            pair_sum      = pair_sum +(length_scale_bohr/distance_bohr)**pair_exponent
        end do
    end do

    embedding_sum = 0.0d0
    do center_atom = 1, atom_count
        density_sum = 0.0d0
        do neighbor_atom = 1, atom_count
            if( neighbor_atom == center_atom ) cycle
            distance_bohr = cartesian_distance( coordinates_bohr(center_atom,:), coordinates_bohr(neighbor_atom,:) )
            density_sum   = density_sum +(length_scale_bohr/distance_bohr)**density_exponent
        end do
        embedding_sum = embedding_sum+sqrt(density_sum)
    end do
    energy_hartree = pair_scale*pair_sum-embedding_scale*embedding_sum
end subroutine evaluate_reparameterized_sutton_chen
