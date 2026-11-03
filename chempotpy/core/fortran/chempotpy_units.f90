!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Provide canonical unit conversions for new atomic-unit interfaces.
module chempotpy_units
    use chempotpy_kinds, only: real64
    implicit none
    private
    public :: angstrom_to_bohr
    public :: bohr_to_angstrom
    public :: electronvolt_to_hartree
    public :: hartree_per_bohr_to_ev_per_angstrom
    public :: hartree_to_electronvolt
    ! These values match chempotpy.core.units and are used only by new APIs.
    real(real64), parameter :: bohr_to_angstrom = 0.529177210903_real64
    real(real64), parameter :: angstrom_to_bohr = 1.0_real64 &
        /bohr_to_angstrom
    real(real64), parameter :: hartree_to_electronvolt = &
        27.211386245988_real64
    real(real64), parameter :: electronvolt_to_hartree = 1.0_real64 &
        /hartree_to_electronvolt
    real(real64), parameter :: hartree_per_bohr_to_ev_per_angstrom = &
        hartree_to_electronvolt/bohr_to_angstrom
end module chempotpy_units
