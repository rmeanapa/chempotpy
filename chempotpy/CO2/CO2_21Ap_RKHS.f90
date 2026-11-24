! System:                        CO2
! Functional form:               reproducing-kernel Hilbert-space interpolation
! Common name:                   CO2 (2)1A'
! Reference:                     J. C. San Vicente Veliz et al., PCCP 23, 11251 (2021)
! Number of derivatives:         1
! Number of bodies:              3
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable grids and coefficients are embedded
! Notes:                         Atom order is O, C, O; energy is returned in eV.

!> @brief Evaluate the (2)1A' reactive RKHS surface for C + O2 and O + CO.
module chempotpy_co2_21ap_rkhs_pes
    use chempotpy_kinds,               only: wp
    use chempotpy_co2_21ap_parameters, only: angular_grid, co_channel_alpha, co_first_radial_grid, &
        & co_second_radial_grid, oo_channel_alpha, oo_first_radial_grid, oo_second_radial_grid
    implicit none
    private
    public :: evaluate
    real(wp), parameter :: switching_scales(3) = [1.05_wp, 1.2_wp, 1.2_wp]
    real(wp), parameter :: channel_offset      = 7.73416186112261816e-02_wp

    contains

    !> @brief Evaluate energy and optional Cartesian gradient for coordinates in angstrom.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_co2_rkhs_runtime, only: evaluate_co2_rkhs_surface
        use chempotpy_status,           only: chempotpy_status_t
        type(chempotpy_status_t), intent(out), optional :: status
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(3,3)
        real(wp), intent(out) :: p(1), g(1,3,3), d(1,1,3,3)
        call evaluate_co2_rkhs_surface(x, igrad, p, g, d, angular_grid, oo_first_radial_grid, &
            & oo_second_radial_grid, oo_channel_alpha, co_first_radial_grid, co_second_radial_grid, &
            & co_channel_alpha, switching_scales, channel_offset, status)
    end subroutine evaluate
end module chempotpy_co2_21ap_rkhs_pes

!> @brief Expose CO2_21Ap_RKHS through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_kinds,             only: wp
    use chempotpy_co2_21ap_rkhs_pes, only: evaluate
    implicit none
    integer, parameter    :: natoms = 3, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
