!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the reparameterized Streitz-Mintmire potential through Python.
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
    call evaluate_reparameterized_streitz_mintmire( &
        coordinates_bohr, atom_count, energy_hartree )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*,*) 'Only energy is available'
    endif
end subroutine pes

!> @brief Evaluate the reparameterized Streitz-Mintmire cluster energy.
subroutine evaluate_reparameterized_streitz_mintmire( &
        coordinates_bohr, atom_count, energy_hartree )
    use chempotpy_geometry, only: cartesian_distance
    use chempotpy_kinds, only: real32, real64
    implicit none

    real(real64), parameter :: hartree_to_ev    = 27.2113961d0
    real(real64), parameter :: bohr_to_angstrom = 0.529177249d0
    ! Preserve default-real rounding of the released fit parameters.
    real(real64), parameter :: reference_distance_bohr = &
        real(2.85197850512945772_real32,real64)/bohr_to_angstrom
    real(real64), parameter :: pair_decay = &
        real(2.02589154860771847_real32,real64)*bohr_to_angstrom
    real(real64), parameter :: density_decay = &
        real(1.24035173424523704_real32,real64)*bohr_to_angstrom
    real(real64), parameter :: first_pair_scale = &
        real(0.850952613580850115e-1_real32,real64)/hartree_to_ev
    real(real64), parameter :: second_pair_scale = &
        real(0.410058622374206150_real32,real64)/hartree_to_ev
    real(real64), parameter :: embedding_scale = &
        real(0.687103077674645846_real32,real64)/hartree_to_ev
    real(real64), intent(in)  :: coordinates_bohr(atom_count,3)
    integer, intent(in)       :: atom_count
    real(real64), intent(out) :: energy_hartree
    real(real64)              :: distance_bohr, radial_offset, density
    real(real64)              :: pair_energy, embedding_energy, first_term, second_term
    integer                   :: center_atom, neighbor_atom

    pair_energy      = 0.0d0
    embedding_energy = 0.0d0
    do center_atom = 1, atom_count
        density = 0.0d0
        do neighbor_atom = 1, atom_count
            if( neighbor_atom == center_atom ) cycle
            distance_bohr = cartesian_distance( coordinates_bohr(center_atom,:), coordinates_bohr(neighbor_atom,:) )
            radial_offset = distance_bohr-reference_distance_bohr
            density       = density+exp(-density_decay*radial_offset)
            if( center_atom > neighbor_atom ) then
                first_term  = 2.0d0*first_pair_scale *exp(-0.5d0*density_decay*radial_offset)
                second_term = -second_pair_scale *(1.0d0+pair_decay*radial_offset) *exp(-pair_decay*radial_offset)
                pair_energy = pair_energy+first_term+second_term
            endif
        enddo
        embedding_energy = embedding_energy-embedding_scale*sqrt(density)
    enddo
    energy_hartree = pair_energy+embedding_energy
end subroutine evaluate_reparameterized_streitz_mintmire
