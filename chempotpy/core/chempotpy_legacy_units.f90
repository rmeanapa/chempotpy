!> @brief Preserve conversion behavior of the released CHEMPOTPY wrappers.
module chempotpy_legacy_units
    use chempotpy_kinds, only: real32, wp
    implicit none
    private
    public :: legacy_bohr_to_angstrom
    public :: legacy_coordinates_to_bohr
    public :: legacy_energy_to_electronvolt
    public :: legacy_gradient_to_electronvolt_per_angstrom
    public :: legacy_hartree_to_electronvolt
    ! The released literals were default real and then promoted to wp.
    real(wp), parameter :: legacy_bohr_to_angstrom = &
        real(0.529177211_real32, wp)
    real(wp), parameter :: legacy_hartree_to_electronvolt = &
        real(27.211386_real32, wp)
    ! Preserve the independently rounded gradient factor used by released PES wrappers.
    real(wp), parameter :: legacy_gradient_to_electronvolt_per_angstrom = 51.422067_wp

    contains

    !> @brief Convert wrapper coordinates using the released rounded factor.
    pure subroutine legacy_coordinates_to_bohr( &
            coordinates_angstrom, coordinates_bohr )
        real(wp), intent(in)  :: coordinates_angstrom(:,:)
        real(wp), intent(out) :: coordinates_bohr( &
            size(coordinates_angstrom, 1),size(coordinates_angstrom, 2))
        coordinates_bohr = coordinates_angstrom/legacy_bohr_to_angstrom
    end subroutine legacy_coordinates_to_bohr

    !> @brief Convert wrapper energy using the released rounded factor.
    pure function legacy_energy_to_electronvolt( &
            energy_hartree ) result(energy_electronvolt)
        real(wp), intent(in) :: energy_hartree
        real(wp)             :: energy_electronvolt
        energy_electronvolt = energy_hartree *legacy_hartree_to_electronvolt
    end function legacy_energy_to_electronvolt
end module chempotpy_legacy_units
