!> @brief Store immutable EAM5 fit and energy-zero parameters for H on Ni(111).
module chempotpy_hni_eam5_parameters
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: canonical_asymptote_correction_electronvolt, energy_zero_hartree
    public :: h_embedding_parameters, h_pair_parameters, ni_embedding_parameters, ni_pair_parameters
    real(wp), parameter :: energy_zero_hartree                         = 178.49075_wp
    real(wp), parameter :: canonical_asymptote_correction_electronvolt = 0.0354010647943141521_wp
    real(wp), parameter :: ni_embedding_parameters(12)                 = [ &
        -126.5009308_wp, 0.3362141252_wp, 18.2047475_wp, 11.02211666_wp, &
        & -209.7682800_wp, 51.76818085_wp, 0.2_wp, 0.21_wp, &
        & -1.657208422e10_wp, -5.13226816e8_wp, -4.377938e6_wp, -19.23463058_wp ]
    real(wp), parameter :: h_embedding_parameters(2) = [-70.5460962_wp, 6.9507109_wp]
    real(wp), parameter :: ni_pair_parameters(3)     = [10.0_wp, 0.8957_wp, 1.8633_wp]
    real(wp), parameter :: h_pair_parameters(3)      = [0.1959_wp, 3.2108_wp, 1.7957_wp]
end module chempotpy_hni_eam5_parameters
