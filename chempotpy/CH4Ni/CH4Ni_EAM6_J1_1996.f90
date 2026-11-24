! System:                        CH4Ni
! Functional form:               EAM/DIM plus modified J1
! Common name:                   methane on Ni(111)
! Reference:                     S. E. Wonchoba and D. G. Truhlar, POTLIB nich4 release (1996)
! Number of derivatives:         1 (centered finite difference)
! Number of bodies:              5 active CH4 atoms plus a fixed 397-atom Ni(111) environment
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable POTLIB parameters and Ni coordinates are embedded
! Notes:                         Public atom order is H,H,H,H,C and Cartesian coordinates are in angstrom.

!> @brief Evaluate the released EAM/DIM plus modified-J1 CH4/Ni(111) potential.
module chempotpy_ch4ni_eam6_j1_1996_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Return energy and an optional numerical gradient in the CHEMPOTPY eV/angstrom contract.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_finite_difference, only: central_finite_difference_gradient
        use chempotpy_status,            only: chempotpy_status_t, status_derivative_unavailable
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms = 5, nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp), parameter   :: finite_difference_step_angstrom = 1.0e-5_wp
        type(chempotpy_status_t) :: evaluation_status
        real(wp) :: flat_gradient(nstates,3*natoms)
        integer  :: atom, direction
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        call evaluation_status%clear()
        if( igrad == 0 )then
            call evaluate_energy_ev(x, p, evaluation_status)
        elseif( igrad == 1 )then
            call central_finite_difference_gradient(evaluate_energy_ev, x, finite_difference_step_angstrom, &
                & p, flat_gradient, evaluation_status)
            if( .not. evaluation_status%failed() )then
                do atom=1, natoms
                    do direction=1, 3
                        g(1,atom,direction) = flat_gradient(1,3*(atom-1)+direction)
                    enddo
                enddo
            endif
        else
            call evaluation_status%set(status_derivative_unavailable, &
                & 'CH4Ni_EAM6_J1_1996 provides energy and a numerical first derivative')
        endif
        if( present(status) ) status = evaluation_status
    end subroutine evaluate

    !> @brief Evaluate the one-state CH4/Ni(111) energy in electronvolt for finite-difference reuse.
    pure subroutine evaluate_energy_ev( coordinates, values, status )
        use chempotpy_ch4ni_runtime, only: evaluate_ch4ni_energy
        use chempotpy_status,        only: chempotpy_status_t
        use chempotpy_units,         only: hartree_to_electronvolt
        real(wp), intent(in)  :: coordinates(:,:)
        real(wp), intent(out) :: values(:)
        type(chempotpy_status_t), intent(inout) :: status
        real(wp) :: energy_hartree
        if( status%failed() ) return
        call evaluate_ch4ni_energy(coordinates, energy_hartree)
        values(1) = energy_hartree*hartree_to_electronvolt
    end subroutine evaluate_energy_ev
end module chempotpy_ch4ni_eam6_j1_1996_pes

!> @brief Provide the Python-facing CH4/Ni(111) PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_ch4ni_eam6_j1_1996_pes, only: evaluate
    use chempotpy_kinds,                  only: wp
    implicit none
    integer, parameter    :: natoms = 5, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
