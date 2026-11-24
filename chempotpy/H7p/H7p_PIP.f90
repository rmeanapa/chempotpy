! System:                        H7p
! Functional form:               permutationally invariant polynomial
! Common name:                   H7+
! Reference:                     C. Qu, R. Prosmiti, and J. M. Bowman, Theor. Chem. Acc. 132, 1413 (2013)
! Number of derivatives:         0
! Number of bodies:              7
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable parameters are embedded
! Notes:                         Released coefficients are embedded in the owning parameter module.

!> @brief Implement H7p_PIP, a potential energy surface for H7p.
module chempotpy_h7p_pip_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H7p_PIP potential energy surface for H7p.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_geometry,           only: fill_distance_matrix
        use chempotpy_h7p_pip_basis,      only: h7p_invariant_basis
        use chempotpy_h7p_pip_parameters, only: pip_coefficients
        use chempotpy_status,             only: chempotpy_status_t, status_derivative_unavailable
        use chempotpy_units,              only: bohr_to_angstrom, hartree_to_electronvolt
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms  = 7
        integer, parameter    :: nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: coordinates_bohr(natoms,3), distances(natoms,natoms)
        real(wp)              :: basis(size(pip_coefficients))
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad /= 0 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'H7p_PIP provides energies but not analytic derivatives')
            return
        endif
        coordinates_bohr = x/bohr_to_angstrom
        call fill_distance_matrix(coordinates_bohr, distances)
        call h7p_invariant_basis(distances, basis)
        p(1) = dot_product(pip_coefficients, basis)*hartree_to_electronvolt
    end subroutine evaluate
end module chempotpy_h7p_pip_pes

!> @brief Expose H7p_PIP through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h7p_pip_pes, only: evaluate
    use chempotpy_kinds,       only: wp
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
