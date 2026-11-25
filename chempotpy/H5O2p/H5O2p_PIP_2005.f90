! System:                        H5O2p
! Functional form:               permutationally invariant polynomial
! Common name:                   H5O2+ PIP
! Reference:                     X. Huang, B. J. Braams, and J. M. Bowman, J. Chem. Phys. 122, 044308 (2005), doi:10.1063/1.1834500
! Number of derivatives:         1
! Number of bodies:              7
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable coefficients are embedded
! Notes:                         The public order is O, O, H, H, H, H, H.

!> @brief Implement the 2005 H5O2+ permutation-invariant polynomial surface.
module chempotpy_h5o2p_pip_2005_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H5O2p_PIP_2005 energy and optional analytic gradient.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_geometry, only: all_pair_distances_and_jacobian, reshape_cartesian_coordinates, &
            transform_internal_gradient_to_cartesian
        use chempotpy_h5o2p_pip_basis,      only: evaluate_pip
        use chempotpy_h5o2p_pip_parameters, only: pip_coefficients
        use chempotpy_status,               only: chempotpy_status_t, status_derivative_unavailable
        use chempotpy_units,                only: bohr_to_angstrom, hartree_per_bohr_to_ev_per_angstrom, &
            hartree_to_electronvolt
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms  = 7
        integer, parameter    :: nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: coordinates_bohr(natoms,3), distance_gradient(21), distances(21)
        real(wp)              :: distance_jacobian(21,21), energy_hartree, gradient_flat(21)
        real(wp)              :: gradient_internal(natoms,3)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad < 0.or.igrad > 1 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'H5O2p_PIP_2005 provides energies and first derivatives')
        else
            coordinates_bohr(1:5,:) = x(3:7,:)/bohr_to_angstrom
            coordinates_bohr(6:7,:) = x(1:2,:)/bohr_to_angstrom
            call all_pair_distances_and_jacobian(coordinates_bohr, distances, distance_jacobian)
            call evaluate_pip(distances, pip_coefficients, igrad == 1, energy_hartree, distance_gradient)
            p(1) = energy_hartree*hartree_to_electronvolt
            if( igrad == 1 )then
                call transform_internal_gradient_to_cartesian(distance_gradient, distance_jacobian, gradient_flat)
                call reshape_cartesian_coordinates(gradient_flat, gradient_internal)
                g(1,1:2,:) = gradient_internal(6:7,:)*hartree_per_bohr_to_ev_per_angstrom
                g(1,3:7,:) = gradient_internal(1:5,:)*hartree_per_bohr_to_ev_per_angstrom
            endif
        endif
    end subroutine evaluate
end module chempotpy_h5o2p_pip_2005_pes

!> @brief Expose H5O2p_PIP_2005 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h5o2p_pip_2005_pes, only: evaluate
    use chempotpy_kinds,              only: wp
    implicit none
    integer, parameter    :: natoms  = 7
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
