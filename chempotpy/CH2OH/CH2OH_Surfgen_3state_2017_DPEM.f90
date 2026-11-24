! System:                        CH2OH
! Functional form:               third-order symmetry-adapted quasi-diabatic expansion
! Common name:                   CH2OH three-state Surfgen DPEM
! Number of derivatives:         1
! Number of bodies:              5
! Number of electronic surfaces: 3
! Interface:                     DPEM
! Data file:                     none; immutable hd.data coefficients are embedded
! Notes:                         Canonical atom order is O, C, H, H, hydroxyl H.
! Reference: C. L. Malbon and D. R. Yarkony, Journal of Chemical Physics
!            146, 134302 (2017).

!> @brief Implement the three-state CH2OH Surfgen diabatic potential matrix released in 2017.
module chempotpy_ch2oh_surfgen_3state_2017_dpem_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the CH2OH three-state diabatic Hamiltonian and optional first derivative.
    subroutine evaluate( x, derivative_order, hamiltonian, gradient, status )
        use chempotpy_ch2oh_3state_parameters, only: coefficient_column, coefficient_order, coefficient_row, &
            & coefficient_term, coefficient_value
        use chempotpy_ch2oh_surfgen_runtime, only: evaluate_ch2oh_diabatic
        use chempotpy_status,                only: chempotpy_status_t
        type(chempotpy_status_t), intent(out), optional :: status
        integer, intent(in)   :: derivative_order
        real(wp), intent(in)  :: x(5,3)
        real(wp), intent(out) :: hamiltonian(3,3), gradient(3,3,5,3)
        call evaluate_ch2oh_diabatic(x, derivative_order, coefficient_order, coefficient_row, coefficient_column, &
            & coefficient_term, coefficient_value, hamiltonian, gradient, status)
    end subroutine evaluate
end module chempotpy_ch2oh_surfgen_3state_2017_dpem_pes

!> @brief Expose the three-state CH2OH Surfgen model through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_ch2oh_surfgen_3state_2017_dpem_pes, only: evaluate
    use chempotpy_kinds,                              only: wp
    implicit none
    integer, parameter    :: natoms = 5, nstates = 3
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: u(nstates,nstates), ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem
