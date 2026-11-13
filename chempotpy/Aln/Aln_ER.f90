!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_er_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the ER aluminum potential through the Python contract.
    subroutine evaluate( coordinates_angstrom, atom_count, derivative_order, &
            energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
        use chempotpy_units, only: coordinates_to_bohr, &
            energy_to_electronvolt
        use chempotpy_kinds, only: wp
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
        call evaluate_er(coordinates_bohr, atom_count, energy_hartree)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate the reparameterized cutoff extended-Rydberg pair potential.
    subroutine evaluate_er( coordinates_bohr, atom_count, energy_hartree )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), parameter :: hartree_to_ev    = 27.2113961e0_wp
        real(wp), parameter :: bohr_to_angstrom = 0.529177249e0_wp
        ! Preserve the default-real parameter rounding in the released kernel.
        real(wp), parameter :: well_depth_hartree = &
            1.1249634027481079102e0_wp/hartree_to_ev
        real(wp), parameter :: equilibrium_distance_bohr = &
            3.0000977516174316406e0_wp/bohr_to_angstrom
        real(wp), parameter :: linear_coefficient    = 1.6509038209915161133e0_wp
        real(wp), parameter :: quadratic_coefficient = &
            0.86853933334350585938e0_wp
        real(wp), parameter :: cubic_coefficient = &
            0.15108940005302429199e0_wp
        real(wp), parameter   :: cutoff_shape         = 0.39692232012748718262e0_wp
        real(wp), parameter   :: cutoff_distance_bohr = 6.5e0_wp/bohr_to_angstrom
        real(wp), intent(in)  :: coordinates_bohr(atom_count,3)
        integer, intent(in)   :: atom_count
        real(wp), intent(out) :: energy_hartree
        real(wp)              :: distance_bohr
        real(wp)              :: radial_offset, polynomial, exponential, cutoff
        integer               :: first_atom, second_atom
        energy_hartree = 0.0e0_wp
        do first_atom = 1, atom_count-1
            do second_atom = first_atom+1, atom_count
                distance_bohr = cartesian_distance( coordinates_bohr(first_atom,:), coordinates_bohr(second_atom,:) )
                if( distance_bohr <= cutoff_distance_bohr )then
                    cutoff        = exp(cutoff_shape*(1.0e0_wp -1.0e0_wp/(1.0e0_wp-distance_bohr/cutoff_distance_bohr)))
                    radial_offset = distance_bohr-equilibrium_distance_bohr
                    polynomial = 1.0e0_wp+linear_coefficient*radial_offset +quadratic_coefficient*radial_offset**2 + &
                        & cubic_coefficient*radial_offset**3
                    exponential    = exp(-linear_coefficient*radial_offset)
                    energy_hartree = energy_hartree -well_depth_hartree*exponential*polynomial*cutoff
                endif
            enddo
        enddo
    end subroutine evaluate_er
end module chempotpy_aln_er_pes

!> @brief Expose Aln_ER through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_er_pes, only: evaluate
    use chempotpy_kinds,      only: wp
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
