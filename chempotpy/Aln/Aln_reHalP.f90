!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the reparameterized Halicioglu-Pound potential through Python.
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
        coordinates_angstrom(1:atom_count,:), coordinates_bohr )
    call evaluate_reparameterized_halicioglu_pound( &
        coordinates_bohr, atom_count, energy_hartree )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*,*) 'Only energy is available'
    endif
end subroutine pes


!> @brief Evaluate the generalized Lennard-Jones aluminum cluster energy.
subroutine evaluate_reparameterized_halicioglu_pound( &
        coordinates_bohr, atom_count, energy_hartree )
    use chempotpy_geometry, only: cartesian_distance
    use chempotpy_kinds, only: real64
    implicit none



    real(real64), parameter :: hartree_to_ev    = 27.2113961d0
    real(real64), parameter :: bohr_to_angstrom = 0.529177249d0
    ! Preserve the default-real parameter rounding in the released kernel.
    real(real64), parameter :: sigma_bohr = &
        2.5383489131927490234d0/bohr_to_angstrom
    real(real64), parameter :: well_scale_hartree = &
        0.23957009613513946533d0/hartree_to_ev
    real(real64), parameter :: attractive_exponent = &
        5.8041033744812011719d0
    real(real64), parameter :: repulsive_exponent = &
        8.1436252593994140625d0
    real(real64), intent(in)  :: coordinates_bohr(atom_count,3)
    integer, intent(in)       :: atom_count
    real(real64), intent(out) :: energy_hartree
    real(real64)              :: distance_bohr
    real(real64)              :: repulsive_term, attractive_term, pair_sum
    integer                   :: first_atom, second_atom

    pair_sum = 0.0d0
    do first_atom = 1, atom_count-1
        do second_atom = first_atom+1, atom_count
            distance_bohr   = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
            repulsive_term  = (sigma_bohr/distance_bohr) **repulsive_exponent
            attractive_term = (sigma_bohr/distance_bohr) **attractive_exponent
            pair_sum        = pair_sum+repulsive_term-attractive_term
        enddo
    enddo
    energy_hartree = well_scale_hartree*pair_sum
end subroutine evaluate_reparameterized_halicioglu_pound
