!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the ER aluminum potential through the Python contract.
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
    call evaluate_er(coordinates_bohr, atom_count, energy_hartree)
    if( derivative_order == 0 )then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*, *) 'Only energy is available'
    endif
end subroutine pes

!> @brief Evaluate the reparameterized cutoff extended-Rydberg pair potential.
subroutine evaluate_er( coordinates_bohr, atom_count, energy_hartree )
    use chempotpy_geometry, only: cartesian_distance
    use chempotpy_kinds,    only: real64
    implicit none
    real(real64), parameter :: hartree_to_ev    = 27.2113961d0
    real(real64), parameter :: bohr_to_angstrom = 0.529177249d0
    ! Preserve the default-real parameter rounding in the released kernel.
    real(real64), parameter :: well_depth_hartree = &
        1.1249634027481079102d0/hartree_to_ev
    real(real64), parameter :: equilibrium_distance_bohr = &
        3.0000977516174316406d0/bohr_to_angstrom
    real(real64), parameter :: linear_coefficient    = 1.6509038209915161133d0
    real(real64), parameter :: quadratic_coefficient = &
        0.86853933334350585938d0
    real(real64), parameter :: cubic_coefficient = &
        0.15108940005302429199d0
    real(real64), parameter   :: cutoff_shape         = 0.39692232012748718262d0
    real(real64), parameter   :: cutoff_distance_bohr = 6.5d0/bohr_to_angstrom
    real(real64), intent(in)  :: coordinates_bohr(atom_count,3)
    integer, intent(in)       :: atom_count
    real(real64), intent(out) :: energy_hartree
    real(real64)              :: distance_bohr
    real(real64)              :: radial_offset, polynomial, exponential, cutoff
    integer                   :: first_atom, second_atom
    energy_hartree = 0.0d0
    do first_atom = 1, atom_count-1
        do second_atom = first_atom+1, atom_count
            distance_bohr = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
            if( distance_bohr <= cutoff_distance_bohr )then
                cutoff        = exp(cutoff_shape*(1.0d0 -1.0d0/(1.0d0-distance_bohr/cutoff_distance_bohr)))
                radial_offset = distance_bohr-equilibrium_distance_bohr
                polynomial = 1.0d0+linear_coefficient*radial_offset +quadratic_coefficient*radial_offset**2 + &
                    & cubic_coefficient*radial_offset**3
                exponential    = exp(-linear_coefficient*radial_offset)
                energy_hartree = energy_hartree -well_depth_hartree*exponential*polynomial*cutoff
            endif
        enddo
    enddo
end subroutine evaluate_er
