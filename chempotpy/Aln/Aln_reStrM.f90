!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_restrm_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the reparameterized Streitz-Mintmire potential through Python.
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
        call evaluate_reparameterized_streitz_mintmire( &
            coordinates_bohr, atom_count, energy_hartree)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate the reparameterized Streitz-Mintmire cluster energy.
    subroutine evaluate_reparameterized_streitz_mintmire( &
            coordinates_bohr, atom_count, energy_hartree )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: real32, wp
        implicit none
        real(wp), parameter :: hartree_to_ev    = 27.2113961e0_wp
        real(wp), parameter :: bohr_to_angstrom = 0.529177249e0_wp
        ! Preserve default-real rounding of the released fit parameters.
        real(wp), parameter :: reference_distance_bohr = &
            real(2.85197850512945772_real32, wp)/bohr_to_angstrom
        real(wp), parameter :: pair_decay = &
            real(2.02589154860771847_real32, wp)*bohr_to_angstrom
        real(wp), parameter :: density_decay = &
            real(1.24035173424523704_real32, wp)*bohr_to_angstrom
        real(wp), parameter :: first_pair_scale = &
            real(0.850952613580850115e-1_real32, wp)/hartree_to_ev
        real(wp), parameter :: second_pair_scale = &
            real(0.410058622374206150_real32, wp)/hartree_to_ev
        real(wp), parameter :: embedding_scale = &
            real(0.687103077674645846_real32, wp)/hartree_to_ev
        real(wp), intent(in)  :: coordinates_bohr(atom_count,3)
        integer, intent(in)   :: atom_count
        real(wp), intent(out) :: energy_hartree
        real(wp)              :: distance_bohr, radial_offset, density
        real(wp)              :: pair_energy, embedding_energy, first_term, second_term
        integer               :: center_atom, neighbor_atom
        pair_energy      = 0.0e0_wp
        embedding_energy = 0.0e0_wp
        do center_atom = 1, atom_count
            density = 0.0e0_wp
            do neighbor_atom = 1, atom_count
                if( neighbor_atom == center_atom ) cycle
                distance_bohr = cartesian_distance( coordinates_bohr(center_atom,:), coordinates_bohr(neighbor_atom,:) )
                radial_offset = distance_bohr-reference_distance_bohr
                density       = density+exp(-density_decay*radial_offset)
                if( center_atom > neighbor_atom )then
                    first_term = 2.0e0_wp*first_pair_scale *exp(-0.5e0_wp*density_decay*radial_offset)
                    second_term = -second_pair_scale *(1.0e0_wp+pair_decay*radial_offset) *exp(-pair_decay* &
                        & radial_offset)
                    pair_energy = pair_energy+first_term+second_term
                endif
            enddo
            embedding_energy = embedding_energy-embedding_scale*sqrt(density)
        enddo
        energy_hartree = pair_energy+embedding_energy
    end subroutine evaluate_reparameterized_streitz_mintmire
end module chempotpy_aln_restrm_pes

!> @brief Expose Aln_reStrM through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_restrm_pes, only: evaluate
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
