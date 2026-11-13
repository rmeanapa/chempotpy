!> @brief Provide the NIST CODATA 2022 constants used by CHEMPOTPY.
module chempotpy_physical_constants
    use chempotpy_kinds, only: wp
    implicit none
    private
    integer, parameter, public          :: codata_adjustment_year = 2022
    character(len=*), parameter, public :: codata_source_name     = &
        "NIST CODATA 2022"
    character(len=*), parameter, public :: codata_source_url = &
        "https://physics.nist.gov/cuu/pdf/all.pdf"
    real(wp), parameter, public :: avogadro_constant_per_mole      = 6.02214076e23_wp
    real(wp), parameter, public :: bohr_radius_angstrom            = 0.529177210544_wp
    real(wp), parameter, public :: bohr_radius_meter               = 5.29177210544e-11_wp
    real(wp), parameter, public :: elementary_charge_coulomb       = 1.602176634e-19_wp
    real(wp), parameter, public :: hartree_energy_electronvolt     = 27.211386245981_wp
    real(wp), parameter, public :: hartree_energy_joule            = 4.3597447222060e-18_wp
    real(wp), parameter, public :: planck_constant_joule_second    = 6.62607015e-34_wp
    real(wp), parameter, public :: rydberg_constant_per_meter      = 10973731.568157_wp
    real(wp), parameter, public :: speed_of_light_meter_per_second = 299792458.0_wp
end module chempotpy_physical_constants
