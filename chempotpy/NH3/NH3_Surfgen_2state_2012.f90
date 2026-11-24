! System:                        NH3
! Functional form:               fourth-order symmetry-adapted quasi-diabatic expansion
! Common name:                   NH3 two-state Surfgen model
! Number of derivatives:         1
! Number of bodies:              4
! Number of electronic surfaces: 2
! Interface:                     PES
! Data file:                     none; immutable hd.data coefficients are embedded
! Notes:                         Canonical atom order is N, H, H, H.
! Reference: X. Zhu and D. R. Yarkony, J. Chem. Phys. 137, 22A511 (2012).

!> @brief Implement the two-state NH3 Surfgen adiabatic potential released in 2012.
module chempotpy_nh3_surfgen_2state_2012_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate two NH3 adiabatic energies and optional first derivatives.
    subroutine evaluate( x, derivative_order, energy, gradient, coupling, status )
        use chempotpy_nh3_surfgen_parameters, only: coefficient_column, coefficient_order, coefficient_row, &
            & coefficient_term, coefficient_value
        use chempotpy_nh3_surfgen_runtime, only: evaluate_nh3_adiabatic
        use chempotpy_status,              only: chempotpy_status_t
        type(chempotpy_status_t), intent(out), optional :: status
        integer, intent(in)   :: derivative_order
        real(wp), intent(in)  :: x(4,3)
        real(wp), intent(out) :: energy(2), gradient(2,4,3), coupling(2,2,4,3)
        call evaluate_nh3_adiabatic(x, derivative_order, coefficient_order, coefficient_row, coefficient_column, &
            & coefficient_term, coefficient_value, energy, gradient, coupling, status)
    end subroutine evaluate
end module chempotpy_nh3_surfgen_2state_2012_pes

!> @brief Expose the two-state NH3 Surfgen model through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_nh3_surfgen_2state_2012_pes, only: evaluate
    use chempotpy_kinds,                       only: wp
    implicit none
    integer, parameter    :: natoms = 4, nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
