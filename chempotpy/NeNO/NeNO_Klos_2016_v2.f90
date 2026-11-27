! System:                        NeNO
! Functional form:               Vibrationally averaged analytic interaction fit
! Common name:                   Ne-NO(X) v=2
! Reference:                     J. Klos et al., J. Chem. Phys. 144, 114307 (2016); DOI 10.1063/1.4943797
! Number of derivatives:         0
! Number of bodies:              3
! Number of electronic surfaces: 2
! Interface:                     PES
! Data file:                     none; immutable parameters are embedded
! Notes:                         Released at DOI 10.5281/zenodo.46119; atom order Ne, N, O

!> @brief Implement the two-component Klos Ne--NO(X 2Pi), v=2 interaction potential.
module chempotpy_neno_klos_2016_v2_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    integer, parameter :: natoms = 3, nstates = 2

    contains

    !> @brief Evaluate the vibrationally averaged v=2 A-double-prime and A-prime surfaces.
    subroutine evaluate( coordinates, derivative_order, energy, gradient, coupling )
        use chempotpy_atomic_elements,      only: atomic_mass_u
        use chempotpy_constants,            only: PI
        use chempotpy_internal_coordinates, only: rigid_diatom_jacobi_coordinates
        use chempotpy_neno_runtime,         only: evaluate_neno_vibrational_components
        use chempotpy_units,                only: electronvolt_to_inverse_centimeter
        integer, intent(in)   :: derivative_order
        real(wp), intent(in)  :: coordinates(natoms,3)
        real(wp), intent(out) :: energy(nstates), gradient(nstates,natoms,3)
        real(wp), intent(out) :: coupling(nstates,nstates,natoms,3)
        logical               :: valid
        real(wp)              :: angle_degree, cosine_theta, distance_bohr, v_difference, v_sum
        gradient = 0.0_wp; coupling = 0.0_wp
        call rigid_diatom_jacobi_coordinates(coordinates, [atomic_mass_u('N'), atomic_mass_u('O')], &
            distance_bohr, cosine_theta, valid)
        if( .not. valid )then
            energy = huge(1.0_wp)
            return
        endif
        angle_degree = acos(cosine_theta)*180.0_wp/PI
        call evaluate_neno_vibrational_components(distance_bohr, angle_degree, 2, v_sum, v_difference)
        energy = [v_sum+v_difference, v_sum-v_difference]/electronvolt_to_inverse_centimeter
    end subroutine evaluate
end module chempotpy_neno_klos_2016_v2_pes

!> @brief Expose NeNO_Klos_2016_v2 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_kinds,                 only: wp
    use chempotpy_neno_klos_2016_v2_pes, only: evaluate
    implicit none
    integer, parameter    :: natoms = 3, nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
