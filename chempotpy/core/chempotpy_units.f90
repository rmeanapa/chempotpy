!> @brief Derive CHEMPOTPY unit conversions from NIST CODATA 2022 constants.
module chempotpy_units
    use chempotpy_kinds,              only: wp
    use chempotpy_physical_constants, only: avogadro_constant_per_mole, bohr_radius_angstrom, &
        elementary_charge_coulomb, hartree_energy_electronvolt, hartree_energy_joule, &
        planck_constant_joule_second, rydberg_constant_per_meter, speed_of_light_meter_per_second
    implicit none
    private
    public :: angstrom_to_bohr
    public :: bohr_to_angstrom
    public :: coordinates_to_bohr
    public :: electronvolt_to_inverse_centimeter
    public :: electronvolt_to_hartree
    public :: electronvolt_to_kilocalorie_per_mole
    public :: electronvolt_to_kilojoule_per_mole
    public :: energy_to_electronvolt
    public :: hartree_per_bohr_to_ev_per_angstrom
    public :: ev_per_angstrom_to_hartree_per_bohr
    public :: hartree_to_electronvolt
    public :: hartree_to_inverse_centimeter
    public :: hartree_to_joule
    public :: hartree_to_kilocalorie_per_mole
    public :: hartree_to_kilojoule_per_mole
    public :: joule_to_hartree
    real(wp), parameter :: joule_per_thermochemical_kilocalorie = 4184.0_wp
    real(wp), parameter :: bohr_to_angstrom                     = bohr_radius_angstrom
    real(wp), parameter :: angstrom_to_bohr                     = 1.0_wp/bohr_to_angstrom
    real(wp), parameter :: hartree_to_electronvolt              = hartree_energy_electronvolt
    real(wp), parameter :: electronvolt_to_hartree              = 1.0_wp/hartree_to_electronvolt
    real(wp), parameter :: hartree_per_bohr_to_ev_per_angstrom  = &
        hartree_to_electronvolt/bohr_to_angstrom
    real(wp), parameter :: ev_per_angstrom_to_hartree_per_bohr = &
        1.0_wp/hartree_per_bohr_to_ev_per_angstrom
    real(wp), parameter :: hartree_to_joule              = hartree_energy_joule
    real(wp), parameter :: joule_to_hartree              = 1.0_wp/hartree_to_joule
    real(wp), parameter :: hartree_to_kilojoule_per_mole = &
        hartree_energy_joule*avogadro_constant_per_mole/1000.0_wp
    real(wp), parameter :: hartree_to_kilocalorie_per_mole = &
        hartree_energy_joule*avogadro_constant_per_mole/joule_per_thermochemical_kilocalorie
    real(wp), parameter :: electronvolt_to_kilojoule_per_mole = &
        elementary_charge_coulomb*avogadro_constant_per_mole/1000.0_wp
    real(wp), parameter :: electronvolt_to_kilocalorie_per_mole = &
        elementary_charge_coulomb*avogadro_constant_per_mole/joule_per_thermochemical_kilocalorie
    real(wp), parameter :: electronvolt_to_inverse_centimeter = elementary_charge_coulomb/ &
        (planck_constant_joule_second*speed_of_light_meter_per_second*100.0_wp)
    real(wp), parameter :: hartree_to_inverse_centimeter = &
        2.0_wp*rydberg_constant_per_meter/100.0_wp

    contains

    !> @brief Convert Cartesian coordinates from angstrom to bohr.
    pure subroutine coordinates_to_bohr( coordinates_angstrom, coordinates_bohr )
        real(wp), intent(in)  :: coordinates_angstrom(:,:)
        real(wp), intent(out) :: coordinates_bohr( &
            size(coordinates_angstrom, 1),size(coordinates_angstrom, 2))
        coordinates_bohr = coordinates_angstrom*angstrom_to_bohr
    end subroutine coordinates_to_bohr

    !> @brief Convert one energy from hartree to electronvolt.
    pure function energy_to_electronvolt( energy_hartree ) result(energy_electronvolt)
        real(wp), intent(in) :: energy_hartree
        real(wp)             :: energy_electronvolt
        energy_electronvolt = energy_hartree*hartree_to_electronvolt
    end function energy_to_electronvolt
end module chempotpy_units
