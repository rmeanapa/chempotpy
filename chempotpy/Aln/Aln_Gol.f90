!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the Gollisch aluminum potential through the Python contract.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_legacy_units, only: legacy_coordinates_to_bohr, &
        legacy_energy_to_electronvolt
    use chempotpy_kinds, only: real64
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
        coordinates_angstrom(1:atom_count,:), coordinates_bohr)
    call evaluate_gollisch(coordinates_bohr, atom_count, energy_hartree)
    if( derivative_order == 0 )then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*, *) 'Only energy is available'
    endif
end subroutine pes

!> @brief Evaluate the Gollisch embedded-atom aluminum cluster energy.
subroutine evaluate_gollisch( coordinates_bohr, atom_count, energy_hartree )
    use chempotpy_geometry, only: cartesian_distance
    use chempotpy_kinds,    only: real64
    implicit none
    real(real64), parameter :: hartree_to_ev    = 27.2113961d0
    real(real64), parameter :: bohr_to_angstrom = 0.529177249d0
    real(real64), parameter :: repulsive_scale  = 350.415d0/hartree_to_ev
    real(real64), parameter :: embedding_scale  = 10.3812d0/hartree_to_ev
    real(real64), parameter :: repulsive_decay  = 2.7338d0*bohr_to_angstrom
    ! The released source used a default-real literal here; retain its rounding.
    real(real64), parameter   :: density_decay      = 1.3681*bohr_to_angstrom
    real(real64), parameter   :: embedding_exponent = 0.6d0
    real(real64), intent(in)  :: coordinates_bohr(atom_count,3)
    integer, intent(in)       :: atom_count
    real(real64), intent(out) :: energy_hartree
    real(real64)              :: distance_bohr
    real(real64)              :: repulsive_sum, density_sum
    integer                   :: center_atom, neighbor_atom
    energy_hartree = 0.0d0
    do center_atom = 1, atom_count
        repulsive_sum = 0.0d0
        density_sum   = 0.0d0
        do neighbor_atom = 1, atom_count
            if( center_atom /= neighbor_atom )then
                distance_bohr = cartesian_distance( coordinates_bohr(center_atom,:), coordinates_bohr(neighbor_atom, &
                    & :) )
                repulsive_sum = repulsive_sum +exp(-distance_bohr*repulsive_decay)
                density_sum   = density_sum +exp(-distance_bohr*density_decay)
            endif
        enddo
        energy_hartree = energy_hartree +repulsive_scale*repulsive_sum -embedding_scale*density_sum**embedding_exponent
    enddo
end subroutine evaluate_gollisch
