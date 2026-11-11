!> @brief Store immutable coefficients for F2H2_GEN_JBKKL_1990.
! Generated from the released bjkkl.par model data.
! Source SHA-256: 10e2c911d6975a0b8caae9aa55db283aec5990f287b306d57391d1a58856cb78
module chempotpy_f2h2_jbkkl_parameters
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: jbkkl_coefficients
    real(real64), parameter :: jbkkl_coefficients(359) = [ 1.7372668_real64, 4.9765668_real64, 1.174_real64, - &
        & 2994.4122_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 9.9365320_real64, - &
        & 27.070633_real64, 0.0_real64, 0.29240271_real64, 0.0_real64, -0.021857381_real64, 0.0_real64, - &
        & 0.94214268_real64, 0.92209662_real64, 0.0_real64, -0.094797709_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0022455114_real64, -0.065804962_real64, 0.0_real64, -0.023426045_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, -0.0091445247_real64, 0.00_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 1.7771732_real64, 0.0_real64, 0.0_real64, 0.0031791885_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, -0.0031266377_real64, &
        & 0.0_real64, 0.0_real64, 0.00_real64, 0.0_real64, -0.024630599_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, -0.0057064317_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0047412395_real64, 0.0_real64, 0.0_real64, 0.0_real64, - &
        & 0.027199911_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.029551319_real64, 0.0_real64, 0.014124155_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0028003588_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 2.2180410_real64, 0.0_real64, -0.0057110422_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, -0.0015976308_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, - &
        & 0.0011884005_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 1.7468293_real64, 0.0_real64, 0.0049164363_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0011594877_real64, 0.0_real64, &
        & 0.0_real64, 0.033895262_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 1.1714284_real64, &
        & 0.0_real64, 0.0048108200_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 2.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.28484311_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 2.0_real64, -0.0013708544_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, -0.0015690213_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0017600947_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        & 0.0_real64, 0.0_real64, 0.18437241_real64, 0.0_real64, 0.0_real64 ]
end module chempotpy_f2h2_jbkkl_parameters
