! System:                        C2H4O4
! Functional form:               permutationally invariant polynomial
! Common name:                   formic acid dimer
! Reference:                     C. Qu and J. M. Bowman, PCCP 18, 24835 (2016)
! Number of derivatives:         0
! Number of bodies:              10
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable fit coefficients are embedded
! Notes:                         Atom order is H4-O4-C2; total energy is returned in eV.

!> @brief Evaluate the Qu-Bowman permutationally invariant formic-acid dimer surface.
module chempotpy_c2h4o4_fad_pip_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the formic-acid dimer energy for Cartesian coordinates in angstrom.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_c2h4o4_invariants, only: evaluate_fad_fit
        use chempotpy_c2h4o4_parameters, only: fad_coefficients
        use chempotpy_geometry,          only: fill_distance_matrix
        use chempotpy_status,            only: chempotpy_status_t, status_derivative_unavailable
        use chempotpy_units,             only: bohr_to_angstrom, hartree_to_electronvolt
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms = 10, nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: coordinates_bohr(natoms,3), distances(natoms,natoms)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad == 0 )then
            coordinates_bohr = x/bohr_to_angstrom
            call fill_distance_matrix(coordinates_bohr, distances)
            p(1) = evaluate_fad_fit(distances, fad_coefficients)*hartree_to_electronvolt
        elseif( present(status) )then
            call status%set(status_derivative_unavailable, &
                'C2H4O4_FAD_PIP provides energies but not analytic derivatives')
        endif
    end subroutine evaluate
end module chempotpy_c2h4o4_fad_pip_pes

!> @brief Expose C2H4O4_FAD_PIP through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_c2h4o4_fad_pip_pes, only: evaluate
    use chempotpy_kinds,              only: wp
    implicit none
    integer, parameter    :: natoms = 10, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
