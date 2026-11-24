! System:                        CH2OH
! Functional form:               third-order symmetry-adapted quasi-diabatic expansion
! Common name:                   CH2OH two-state Surfgen model
! Number of derivatives:         1
! Number of bodies:              5
! Number of electronic surfaces: 2
! Interface:                     PES
! Data file:                     none; immutable hd.data coefficients are embedded
! Notes:                         Canonical atom order is O, C, H, H, hydroxyl H.
! Reference: C. Xie, C. L. Malbon, D. R. Yarkony, D. Xie, and H. Guo,
!            Journal of Chemical Physics (2018), released through POTLIB.

!> @brief Implement the two-state CH2OH Surfgen adiabatic potential released in 2018.
module chempotpy_ch2oh_surfgen_2state_2018_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate two CH2OH adiabatic energies and optional first derivatives.
    subroutine evaluate( x, derivative_order, energy, gradient, coupling, status )
        use chempotpy_ch2oh_2state_parameters, only: coefficient_column, coefficient_order, coefficient_row, &
            & coefficient_term, coefficient_value
        use chempotpy_ch2oh_surfgen_runtime, only: evaluate_ch2oh_adiabatic
        use chempotpy_status,                only: chempotpy_status_t
        type(chempotpy_status_t), intent(out), optional :: status
        integer, intent(in)   :: derivative_order
        real(wp), intent(in)  :: x(5,3)
        real(wp), intent(out) :: energy(2), gradient(2,5,3), coupling(2,2,5,3)
        call evaluate_ch2oh_adiabatic(x, derivative_order, coefficient_order, coefficient_row, coefficient_column, &
            & coefficient_term, coefficient_value, energy, gradient, coupling, status)
    end subroutine evaluate
end module chempotpy_ch2oh_surfgen_2state_2018_pes

!> @brief Expose the two-state CH2OH Surfgen model through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_ch2oh_surfgen_2state_2018_pes, only: evaluate
    use chempotpy_kinds,                         only: wp
    implicit none
    integer, parameter    :: natoms = 5, nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
