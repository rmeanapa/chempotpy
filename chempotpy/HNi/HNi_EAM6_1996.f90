! System:                        HNi
! Functional form:               embedded-atom method
! Common name:                   EAM6 H on Ni(111)
! Reference:                     S. E. Wonchoba and D. G. Truhlar, Phys. Rev. B 53, 11222 (1996)
! Number of derivatives:         1
! Number of bodies:              1 active H plus a fixed 698-atom Ni(111) environment
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable POTLIB potfile1 parameters are embedded
! Notes:                         The public coordinate is the active H atom; the canonical Ni(111) slab is fixed.

!> @brief Evaluate the EAM6 interaction of one active H atom with the canonical Ni(111) environment.
module chempotpy_hni_eam6_1996_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Return EAM6 energy and analytic H-atom gradient in the CHEMPOTPY eV/angstrom contract.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_hni_eam67_parameters, only: eam6_asymptote_correction_electronvolt, energy_zero_hartree
        use chempotpy_hni_eam67_runtime,    only: evaluate_fixed_hydrogen_eam
        use chempotpy_hni_ni111_parameters, only: fixed_ni_coordinates_bohr
        use chempotpy_status,               only: chempotpy_status_t, status_derivative_unavailable
        use chempotpy_units,                only: bohr_to_angstrom, hartree_to_electronvolt
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms = 1, nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp), allocatable :: coordinates(:,:)
        real(wp)              :: energy_electronvolt, gradient(3)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad == 0 .or. igrad == 1 )then
            allocate(coordinates(1+size(fixed_ni_coordinates_bohr, 1),3))
            coordinates(1,:)  = x(1,:)
            coordinates(2:,:) = fixed_ni_coordinates_bohr*bohr_to_angstrom
            call evaluate_fixed_hydrogen_eam(coordinates, 6, energy_electronvolt, gradient)
            p(1) = energy_electronvolt+energy_zero_hartree*hartree_to_electronvolt+ &
                & eam6_asymptote_correction_electronvolt
            if( igrad == 1 ) g(1,1,:) = gradient
        elseif( present(status) )then
            call status%set(status_derivative_unavailable, 'HNi_EAM6_1996 provides energy and first derivatives')
        endif
    end subroutine evaluate
end module chempotpy_hni_eam6_1996_pes

!> @brief Provide the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_hni_eam6_1996_pes, only: evaluate
    use chempotpy_kinds,             only: wp
    implicit none
    integer, parameter    :: natoms = 1, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
