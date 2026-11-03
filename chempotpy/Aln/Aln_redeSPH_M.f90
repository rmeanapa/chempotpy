!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the reparameterized deSPH/M potential through Python.
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
    call evaluate_reparameterized_desph_morse( &
        coordinates_bohr, atom_count, energy_hartree )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*,*) 'Only energy is available'
    endif
end subroutine pes

!> @brief Evaluate the reparameterized distance-dependent Morse energy.
subroutine evaluate_reparameterized_desph_morse( &
        coordinates_bohr, atom_count, energy_hartree )
    use chempotpy_geometry, only: cartesian_distance
    use chempotpy_kinds, only: real64
    implicit none

    real(real64), parameter :: hartree_to_ev    = 27.2113961d0
    real(real64), parameter :: bohr_to_angstrom = 0.529177249d0
    ! Preserve the default-real parameter rounding in the released kernel.
    real(real64), parameter :: equilibrium_distance_bohr = &
        3.2823643684387207031d0/bohr_to_angstrom
    real(real64), parameter :: well_depth_hartree = &
        0.92818760871887207031d0/hartree_to_ev
    real(real64), parameter :: zeroth_decay = &
        0.85588663816452026367d0*bohr_to_angstrom
    real(real64), parameter :: quadratic_decay = &
        0.7552515864372253418d0*bohr_to_angstrom**3
    real(real64), parameter :: cubic_decay = &
        1.470933079719543457d0*bohr_to_angstrom**4
    real(real64), intent(in)  :: coordinates_bohr(atom_count,3)
    integer, intent(in)       :: atom_count
    real(real64), intent(out) :: energy_hartree
    real(real64)              :: distance_bohr
    real(real64)              :: radial_offset, decay, pair_energy
    integer                   :: first_atom, second_atom

    energy_hartree = 0.0d0
    do first_atom = 1, atom_count-1
        do second_atom = first_atom+1, atom_count
            distance_bohr  = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
            radial_offset  = distance_bohr-equilibrium_distance_bohr
            decay          = zeroth_decay+quadratic_decay*radial_offset**2 +cubic_decay*radial_offset**3
            pair_energy    = well_depth_hartree *(1.0d0-exp(-decay*radial_offset))**2-well_depth_hartree
            energy_hartree = energy_hartree+pair_energy
        enddo
    enddo
end subroutine evaluate_reparameterized_desph_morse
