! System:                        CH2OH
! Functional form:               third-order symmetry-adapted quasi-diabatic expansion
! Common name:                   CH2OH three-state Surfgen model
! Number of derivatives:         1
! Number of bodies:              5
! Number of electronic surfaces: 3
! Interface:                     PES
! Data file:                     none; immutable hd.data coefficients are embedded
! Notes:                         Canonical atom order is O, C, H, H, hydroxyl H.
! Reference: C. L. Malbon and D. R. Yarkony, Journal of Chemical Physics
!            146, 134302 (2017).

!> @brief Implement the three-state CH2OH Surfgen adiabatic potential released in 2017.
module chempotpy_ch2oh_surfgen_3state_2017_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate three CH2OH adiabatic energies and optional first derivatives.
    subroutine evaluate( x, derivative_order, energy, gradient, coupling, status )
        use chempotpy_ch2oh_3state_parameters, only: coefficient_column, coefficient_order, coefficient_row, &
            & coefficient_term, coefficient_value
        use chempotpy_ch2oh_surfgen_runtime, only: evaluate_ch2oh_adiabatic
        use chempotpy_status,                only: chempotpy_status_t
        type(chempotpy_status_t), intent(out), optional :: status
        integer, intent(in)   :: derivative_order
        real(wp), intent(in)  :: x(5,3)
        real(wp), intent(out) :: energy(3), gradient(3,5,3), coupling(3,3,5,3)
        call evaluate_ch2oh_adiabatic(x, derivative_order, coefficient_order, coefficient_row, coefficient_column, &
            & coefficient_term, coefficient_value, energy, gradient, coupling, status)
    end subroutine evaluate
end module chempotpy_ch2oh_surfgen_3state_2017_pes

!> @brief Expose the three-state CH2OH Surfgen model through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_ch2oh_surfgen_3state_2017_pes, only: evaluate
    use chempotpy_kinds,                         only: wp
    implicit none
    integer, parameter    :: natoms = 5, nstates = 3
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
