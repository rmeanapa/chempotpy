!> @brief Store immutable parameters shared by the EAM6 and EAM7 H-Ni surfaces.
module chempotpy_hni_eam67_parameters
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: eam6_asymptote_correction_electronvolt, eam7_asymptote_correction_electronvolt
    public :: energy_zero_hartree, nickel_embedding_parameters, nickel_pair_parameters
    real(wp), parameter :: energy_zero_hartree                    = 178.4662455_wp
    real(wp), parameter :: eam6_asymptote_correction_electronvolt = 0.034873145016717899_wp
    real(wp), parameter :: eam7_asymptote_correction_electronvolt = 0.034873145016717899_wp
    real(wp), parameter :: nickel_embedding_parameters(14)        = [ &
        & -303.28911_wp, 7.64744_wp, 87.98723_wp, 75.07472_wp, -522.77364_wp, 373.37826_wp, 0.1_wp, &
        & 0.11_wp, -2.014563533e10_wp, -538419456.0_wp, -4060125.0_wp, -11.03098054_wp, 39.42118_wp, &
        & 56.34512_wp ]
    real(wp), parameter :: nickel_pair_parameters(3) = [10.0_wp, 0.8957_wp, 1.8633_wp]
end module chempotpy_hni_eam67_parameters
