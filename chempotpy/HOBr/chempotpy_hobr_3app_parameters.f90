!> @brief Store immutable RKHS parameters for HOBr_3App_PIP.
! Generated from the released o_hbr_a_pp.par model data.
! Source SHA-256: ee5d8ff9a341a7c4cee6e2b725559711ba44ec72f81785f34b49f200ec476cf6
module chempotpy_hobr_3app_parameters
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: phi_oh, phi_hbr, rkhs_x_1, rkhs_x_2, rkhs_y_1, &
        & rkhs_coefficients
    real(real64), parameter :: phi_oh(14) = [ 1.83161656_real64, 0.1699472_real64, 0.340069645_real64, 0.77_real64, &
        & 7.00_real64, -7.666934234281819_real64, 4.001607924467653_real64, 0.8169114080858936_real64, &
        & 15.31611785643029_real64, 88.9759472634144_real64, -51.55643183366661_real64, -262.3584874455482_real64, - &
        & 48.44154116389853_real64, 320.5902107457166_real64 ]
    real(real64), parameter :: phi_hbr(14) = [ 2.68060879_real64, 0.1424201_real64, 1.16863709097667_real64, &
        & 0.77_real64, 7.35_real64, -7.024310697817183_real64, 4.812387376658926_real64, 5.627595600751125_real64, &
        & 69.80522802206754_real64, -155.2259704326473_real64, -983.5618355500748_real64, 7949.431079920515_real64, - &
        & 13614.0066694763_real64, 4515.364024267248_real64 ]
    real(real64), parameter :: rkhs_x_1(1110) = [ 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.35_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 1.8_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.5_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, &
        & 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.25_real64, &
        & 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.5_real64, 3.5_real64, &
        & 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, &
        & 4.0_real64, 4.0_real64, 4.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, &
        & 5.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, &
        & 2.75_real64, 2.75_real64, 2.75_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, &
        & 3.0_real64, 3.0_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, &
        & 3.25_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 4.0_real64, &
        & 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, &
        & 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, &
        & 6.0_real64, 6.0_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 3.0_real64, &
        & 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.25_real64, 3.25_real64, &
        & 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.5_real64, 3.5_real64, 3.5_real64, &
        & 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, &
        & 4.0_real64, 4.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, &
        & 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, &
        & 2.75_real64, 2.75_real64, 2.75_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, &
        & 3.0_real64, 3.0_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, &
        & 3.25_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 4.0_real64, &
        & 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, &
        & 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, &
        & 6.0_real64, 6.0_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 3.0_real64, &
        & 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.25_real64, 3.25_real64, &
        & 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.5_real64, 3.5_real64, 3.5_real64, &
        & 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, &
        & 4.0_real64, 4.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, &
        & 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, &
        & 2.75_real64, 2.75_real64, 2.75_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, &
        & 3.0_real64, 3.0_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, &
        & 3.25_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 4.0_real64, &
        & 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, &
        & 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, &
        & 6.0_real64, 6.0_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, 1.35_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, &
        & 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.5_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, 1.65_real64, &
        & 1.65_real64, 1.65_real64, 1.65_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, 1.8_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, &
        & 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, 2.25_real64, &
        & 2.25_real64, 2.25_real64, 2.25_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, 2.5_real64, &
        & 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 2.75_real64, 3.0_real64, &
        & 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 3.25_real64, 3.25_real64, &
        & 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.25_real64, 3.5_real64, 3.5_real64, 3.5_real64, &
        & 3.5_real64, 3.5_real64, 3.5_real64, 3.5_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, &
        & 4.0_real64, 4.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, &
        & 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 6.0_real64, 4.369_real64, &
        & 4.369_real64, 4.369_real64, 4.369_real64, 4.369_real64, 4.369_real64, 4.369_real64, 4.369_real64, &
        & 4.369_real64, 4.469_real64, 4.469_real64, 4.469_real64, 4.469_real64, 4.469_real64, 4.469_real64, &
        & 4.469_real64, 4.469_real64, 4.469_real64, 4.569_real64, 4.569_real64, 4.569_real64, 4.569_real64, &
        & 4.569_real64, 4.569_real64, 4.569_real64, 4.569_real64, 4.569_real64, 2.493_real64, 2.493_real64, &
        & 2.493_real64, 2.493_real64, 2.493_real64, 2.493_real64, 2.493_real64, 2.493_real64, 2.493_real64, &
        & 2.593_real64, 2.593_real64, 2.593_real64, 2.593_real64, 2.593_real64, 2.593_real64, 2.593_real64, &
        & 2.593_real64, 2.593_real64, 2.693_real64, 2.693_real64, 2.693_real64, 2.693_real64, 2.693_real64, &
        & 2.693_real64, 2.693_real64, 2.693_real64, 2.693_real64, 1.739_real64, 1.739_real64, 1.739_real64, &
        & 1.739_real64, 1.739_real64, 1.739_real64, 1.739_real64, 1.739_real64, 1.739_real64, 1.839_real64, &
        & 1.839_real64, 1.839_real64, 1.839_real64, 1.839_real64, 1.839_real64, 1.839_real64, 1.839_real64, &
        & 1.839_real64, 1.939_real64, 1.939_real64, 1.939_real64, 1.939_real64, 1.939_real64, 1.939_real64, &
        & 1.939_real64, 1.939_real64, 1.939_real64 ]
    real(real64), parameter :: rkhs_x_2(1110) = [ 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, &
        & 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, &
        & 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, &
        & 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, &
        & 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, &
        & 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, &
        & 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, &
        & 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, &
        & 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, &
        & 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, &
        & 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, &
        & 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, &
        & 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, &
        & 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, &
        & 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, &
        & 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, &
        & 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, &
        & 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, &
        & 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, &
        & 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, &
        & 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, &
        & 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, &
        & 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, &
        & 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, &
        & 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, &
        & 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, &
        & 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, &
        & 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, &
        & 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, &
        & 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, &
        & 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, &
        & 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, &
        & 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, &
        & 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, &
        & 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, &
        & 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, &
        & 6.85_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 3.6_real64, 3.85_real64, 4.1_real64, 4.35_real64, 4.85_real64, 5.85_real64, 6.85_real64, 2.2_real64, &
        & 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, &
        & 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, &
        & 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, &
        & 2.85_real64, 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, &
        & 3.1_real64, 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, &
        & 3.35_real64, 2.2_real64, 2.35_real64, 2.5_real64, 2.65_real64, 2.85_real64, 3.1_real64, 3.35_real64, &
        & 2.584_real64, 2.584_real64, 2.584_real64, 2.684_real64, 2.684_real64, 2.684_real64, 2.784_real64, &
        & 2.784_real64, 2.784_real64, 2.584_real64, 2.584_real64, 2.584_real64, 2.684_real64, 2.684_real64, &
        & 2.684_real64, 2.784_real64, 2.784_real64, 2.784_real64, 2.584_real64, 2.584_real64, 2.584_real64, &
        & 2.684_real64, 2.684_real64, 2.684_real64, 2.784_real64, 2.784_real64, 2.784_real64, 2.771_real64, &
        & 2.771_real64, 2.771_real64, 2.871_real64, 2.871_real64, 2.871_real64, 2.971_real64, 2.971_real64, &
        & 2.971_real64, 2.771_real64, 2.771_real64, 2.771_real64, 2.871_real64, 2.871_real64, 2.871_real64, &
        & 2.971_real64, 2.971_real64, 2.971_real64, 2.771_real64, 2.771_real64, 2.771_real64, 2.871_real64, &
        & 2.871_real64, 2.871_real64, 2.971_real64, 2.971_real64, 2.971_real64, 4.751_real64, 4.751_real64, &
        & 4.751_real64, 4.851_real64, 4.851_real64, 4.851_real64, 4.951_real64, 4.951_real64, 4.951_real64, &
        & 4.751_real64, 4.751_real64, 4.751_real64, 4.851_real64, 4.851_real64, 4.851_real64, 4.951_real64, &
        & 4.951_real64, 4.951_real64, 4.751_real64, 4.751_real64, 4.751_real64, 4.851_real64, 4.851_real64, &
        & 4.851_real64, 4.951_real64, 4.951_real64, 4.951_real64 ]
    real(real64), parameter :: rkhs_y_1(1110) = [ 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, 0.9999999899999999_real64, &
        & 0.9999999899999999_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, &
        & 0.96984631_real64, 0.96984631_real64, 0.96984631_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, &
        & 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.88302222_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, 0.75_real64, &
        & 0.75_real64, 0.75_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, &
        & 0.58682409_real64, 0.58682409_real64, 0.58682409_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, &
        & 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.41317591_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64, &
        & 0.25_real64, 0.25_real64, 0.96981646_real64, 0.99238872_real64, 0.9999999899999999_real64, &
        & 0.96981646_real64, 0.99238872_real64, 0.9999999899999999_real64, 0.96981646_real64, 0.99238872_real64, &
        & 0.9999999899999999_real64, 0.96981646_real64, 0.99238872_real64, 0.9999999899999999_real64, &
        & 0.96981646_real64, 0.99238872_real64, 0.9999999899999999_real64, 0.96981646_real64, 0.99238872_real64, &
        & 0.9999999899999999_real64, 0.96981646_real64, 0.99238872_real64, 0.9999999899999999_real64, &
        & 0.96981646_real64, 0.99238872_real64, 0.9999999899999999_real64, 0.96981646_real64, 0.99238872_real64, &
        & 0.9999999899999999_real64, 0.78964059_real64, 0.85601302_real64, 0.91156818_real64, 0.78964059_real64, &
        & 0.85601302_real64, 0.91156818_real64, 0.78964059_real64, 0.85601302_real64, 0.91156818_real64, &
        & 0.78964059_real64, 0.85601302_real64, 0.91156818_real64, 0.78964059_real64, 0.85601302_real64, &
        & 0.91156818_real64, 0.78964059_real64, 0.85601302_real64, 0.91156818_real64, 0.78964059_real64, &
        & 0.85601302_real64, 0.91156818_real64, 0.78964059_real64, 0.85601302_real64, 0.91156818_real64, &
        & 0.78964059_real64, 0.85601302_real64, 0.91156818_real64, 0.21895831_real64, 0.29503848_real64, &
        & 0.37734631_real64, 0.21895831_real64, 0.29503848_real64, 0.37734631_real64, 0.21895831_real64, &
        & 0.29503848_real64, 0.37734631_real64, 0.21895831_real64, 0.29503848_real64, 0.37734631_real64, &
        & 0.21895831_real64, 0.29503848_real64, 0.37734631_real64, 0.21895831_real64, 0.29503848_real64, &
        & 0.37734631_real64, 0.21895831_real64, 0.29503848_real64, 0.37734631_real64, 0.21895831_real64, &
        & 0.29503848_real64, 0.37734631_real64, 0.21895831_real64, 0.29503848_real64, 0.37734631_real64 ]
    real(real64), parameter :: rkhs_coefficients(1110) = [ -19458707678.36746_real64, 202683261721.6015_real64, - &
        & 554003880508.4917_real64, 578726187023.6848_real64, -259683752412.1888_real64, 109969130245.8585_real64, - &
        & 94905417712.93542_real64, 290017633467.4608_real64, -1470540908618.166_real64, 2945012862655.228_real64, - &
        & 2256436043370.615_real64, 584266219594.5244_real64, 170011572844.3399_real64, -437337891249.6375_real64, &
        & 345152661952.2921_real64, -1495433519479.391_real64, 2563162845575.633_real64, -1785774720563.554_real64, &
        & 194750054156.2031_real64, 7406165466.207526_real64, 640721915184.7566_real64, -2228851583133.63_real64, &
        & 7264379783108.392_real64, -11822484392594.95_real64, 7981326103245.849_real64, -1879269293795.278_real64, - &
        & 483305004286.3889_real64, 1393111997357.818_real64, -1119569714646.828_real64, 3693939649759.094_real64, - &
        & 4477601406732.44_real64, 1275933110750.474_real64, 1933881521857.858_real64, -698059934981.1385_real64, - &
        & 3821425423296.492_real64, 10649888488286.78_real64, -21112282499733.94_real64, 26965099809007.74_real64, - &
        & 17231787871608.31_real64, 5033124678991.855_real64, -849477997831.8839_real64, -1094724023178.205_real64, &
        & 1159157378119.521_real64, -3154418157925.162_real64, 2593102130193.639_real64, 1135365535424.367_real64, - &
        & 3732689032025.747_real64, 2698645551134.636_real64, -1404976706900.518_real64, -2897280525457.777_real64, &
        & 20663450842135.07_real64, -43791054965722.65_real64, 38924046472767.49_real64, -16843906501937.8_real64, &
        & 7226482323699.544_real64, -1940926456599.577_real64, -83789491301.16164_real64, -267392110535.4463_real64, &
        & 367743823628.8182_real64, 1204759883959.205_real64, -2998915138228.544_real64, -6372481298813.103_real64, &
        & 50493595052247.88_real64, -108228932673125.0_real64, 89388782547038.34_real64, 30366125389171.09_real64, - &
        & 96262006652250.61_real64, 60723522945210.59_real64, -31581666463110.38_real64, 13876553341246.64_real64, - &
        & 754226689788.6443_real64, 949741953376.8622_real64, 11337330414650.64_real64, -39776837858943.06_real64, &
        & 64775780211122.7_real64, -42850536470987.58_real64, -104202689973966.2_real64, 322033681747059.9_real64, - &
        & 344236849618248.8_real64, 51279528894508.66_real64, 185276965015900.4_real64, -150145515532535.0_real64, &
        & 85644857627701.72_real64, -42406918772176.33_real64, 1130651882219.643_real64, 10813383823176.39_real64, - &
        & 108205736992048.1_real64, 303046391981143.9_real64, -478591859657832.9_real64, 479373417622528.8_real64, - &
        & 136383029198662.8_real64, -271627003087678.3_real64, 319444435861356.8_real64, -80338843592569.62_real64, - &
        & 145363335580115.7_real64, 133802495172189.4_real64, -79659722553670.28_real64, 40669371609262.92_real64, - &
        & 10033878201706.29_real64, -31571109097855.2_real64, 431580364960078.5_real64, -1231549569165791.0_real64, &
        & 2010560829475508.0_real64, -2221779607724804.0_real64, 1148927910768610.0_real64, 57021369420617.12_real64, &
        & -35166621799184.29_real64, -803791945967081.4_real64, 2541035950521897.0_real64, -4282577940296464.0_real64, &
        & 4875880343330692.0_real64, -2570480518664088.0_real64, -140182003518927.5_real64, 262241225124471.1_real64, &
        & 599753172211870.0_real64, -2455424159391079.0_real64, 3947709020732038.0_real64, -3779893366818031.0_real64, &
        & 1596967149999965.0_real64, 140534712209753.0_real64, -331186567074309.0_real64, -27612684368133.71_real64, &
        & 807568352466261.1_real64, -752305891655654.9_real64, -824779226018026.2_real64, 1322622390031352.0_real64, - &
        & 62221628627373.88_real64, 163342325511559.5_real64, -146166347275316.4_real64, 145722353097885.0_real64, - &
        & 843411712178159.9_real64, 2381515043198356.0_real64, -1985385350950843.0_real64, -8239190144325.646_real64, &
        & 18571407945859.68_real64, -13633616774721.47_real64, -78193880315138.17_real64, 561763070396395.4_real64, - &
        & 1455728564659004.0_real64, 1172638012727957.0_real64, 39306169641629.77_real64, -108395142649045.8_real64, &
        & 120883936580926.5_real64, -19906462880417.72_real64, -295862236528195.4_real64, 783633317573486.4_real64, - &
        & 622394593208921.4_real64, 40641175405.9417_real64, -333239172710.0114_real64, 853591187969.2491_real64, - &
        & 864139015888.7233_real64, 369270705143.4178_real64, -147411102910.5788_real64, 124868390588.3616_real64, - &
        & 435010985232.0858_real64, 2369768652127.241_real64, -4807399750518.559_real64, 3698723444265.055_real64, - &
        & 975783068463.5348_real64, -201745106387.1486_real64, 608043070616.6604_real64, -550547295501.8901_real64, &
        & 2330285827638.399_real64, -3902321704307.315_real64, 2604112172867.84_real64, -118682929051.0298_real64, - &
        & 162470061203.4325_real64, -797558114573.457_real64, 3206177758677.063_real64, -11393904552657.67_real64, &
        & 18970171497357.2_real64, -12859719516855.5_real64, 3059406231271.781_real64, 552605740163.1603_real64, - &
        & 1912604670660.796_real64, 1693346240301.635_real64, -5561443185904.96_real64, 6594218891329.778_real64, - &
        & 1426227787014.306_real64, -3664856295914.929_real64, 1740918091272.632_real64, 4972631724287.952_real64, - &
        & 15022999512588.69_real64, 31477963532445.92_real64, -41137890949677.22_real64, 26260982846796.06_real64, - &
        & 7634014390245.524_real64, 1533962894857.0_real64, 1340585134470.04_real64, -1659523925998.06_real64, &
        & 4422465402955.525_real64, -3154401335021.637_real64, -3038572214452.183_real64, 7251101142516.306_real64, - &
        & 5307376301794.53_real64, 2749732294467.796_real64, 3760044757115.298_real64, -29577363798441.86_real64, &
        & 62879512241031.23_real64, -55981631479971.62_real64, 24683374254346.61_real64, -11199205762089.33_real64, &
        & 3392838141554.636_real64, -71757641086.8367_real64, 1185420078575.452_real64, -2327249313483.625_real64, &
        & 968494539263.6621_real64, 1548558152806.467_real64, 10865387377199.43_real64, -72302574192262.62_real64, &
        & 153818607562334.1_real64, -129293569561196.5_real64, -37983035619330.31_real64, 133886534980956.7_real64, - &
        & 87164411325191.95_real64, 46757676731979.64_real64, -21015901503534.53_real64, 1381783705138.662_real64, - &
        & 2658644844384.678_real64, -13104199991331.93_real64, 52097828131419.09_real64, -87787209616216.33_real64, &
        & 58334392772758.72_real64, 148144439470368.9_real64, -456104846907016.8_real64, 493067130580092.5_real64, - &
        & 86973856652767.78_real64, -252032252183376.9_real64, 211925235360703.4_real64, -123373205180961.9_real64, &
        & 61439867433397.6_real64, -1821343420783.782_real64, -14349293976629.24_real64, 151073396181412.6_real64, - &
        & 426757626369570.9_real64, 677390961174447.0_real64, -681013193067234.9_real64, 195126248895902.8_real64, &
        & 384158925942268.9_real64, -456187940091768.1_real64, 126273288915474.7_real64, 195161101775285.4_real64, - &
        & 187186694408324.6_real64, 113486948467279.1_real64, -57978614154087.9_real64, 14221520233847.88_real64, &
        & 44811546342892.51_real64, -612741991329674.6_real64, 1750010061563704.0_real64, -2860252478759343.0_real64, &
        & 3164444584872327.0_real64, -1637936350423231.0_real64, -80781527124914.83_real64, 49370460821256.74_real64, &
        & 1142946175151155.0_real64, -3611790661646813.0_real64, 6087466297418169.0_real64, - &
        & 6934147146898280.0_real64, 3658233214317036.0_real64, 199605985713150.2_real64, -374288558737004.7_real64, - &
        & 850134109541829.1_real64, 3485775995951914.0_real64, -5598956688087645.0_real64, 5353786195702064.0_real64, &
        & -2259549463532828.0_real64, -202002156400783.2_real64, 477714030123126.8_real64, 31761281387157.02_real64, - &
        & 1140986742309198.0_real64, 1051824895506076.0_real64, 1204577054698802.0_real64, -1901542451904349.0_real64, &
        & 92936226233952.56_real64, -243723556777136.1_real64, 215296131491116.7_real64, -208207232855920.8_real64, &
        & 1209550678462573.0_real64, -3413004586215728.0_real64, 2842937252864057.0_real64, 6416518868367.973_real64, &
        & -24313078557387.27_real64, 45949304234258.89_real64, 79990308867721.44_real64, -814967440979905.4_real64, &
        & 2142976557139290.0_real64, -1725917031284626.0_real64, -57331282518844.76_real64, 183058213982160.8_real64, &
        & -250517009589033.4_real64, 90327921150752.78_real64, 455077194439446.6_real64, -1249031166813277.0_real64, &
        & 992050788490171.0_real64, -48090874884.87218_real64, 233858611846.3445_real64, -468689125273.7759_real64, &
        & 416747379957.1271_real64, -153921727881.4044_real64, 52212603839.72861_real64, -32467420137.28289_real64, &
        & 172048248184.1093_real64, -1226320468269.89_real64, 2605596430117.238_real64, -2034741962123.177_real64, &
        & 563144733596.5966_real64, 4207556088.823776_real64, -184106361970.8206_real64, 350898111516.5049_real64, - &
        & 1359282817004.277_real64, 2091300454296.59_real64, -1253722296135.805_real64, -82915536103.66129_real64, &
        & 218465720084.7429_real64, 121526067084.1721_real64, -1098776503569.968_real64, 5477088439076.15_real64, - &
        & 9823910660907.643_real64, 6758248508666.674_real64, -1660857127361.974_real64, 37294827786.95081_real64, &
        & 534390099379.5754_real64, -895471992418.203_real64, 2901469819554.759_real64, -3284484448875.297_real64, &
        & 317178718106.1541_real64, 2503632572599.244_real64, -1670388586126.175_real64, -956362837976.2288_real64, &
        & 4771412526459.885_real64, -12868840909332.78_real64, 18392866456749.04_real64, -11766495277896.82_real64, &
        & 3424030260935.904_real64, -1127531098873.0_real64, -55152588372.5625_real64, 754967772857.9806_real64, - &
        & 1900856812415.522_real64, 723978207723.9702_real64, 3121653233167.161_real64, -5556953102058.064_real64, &
        & 4182579246941.428_real64, -2345357537017.144_real64, -340616408375.4686_real64, 10191458010678.38_real64, - &
        & 22563466243912.86_real64, 20408733331203.4_real64, -9899781032930.523_real64, 5660019119143.145_real64, - &
        & 2424452953174.829_real64, 257787310802.0989_real64, -1575379464427.345_real64, 3875207004431.055_real64, - &
        & 5164116954107.122_real64, 4265710308583.458_real64, -6865990390380.666_real64, 25910423837834.38_real64, - &
        & 52929416259170.5_real64, 47576295592948.08_real64, 5705480958082.125_real64, -42031555849933.1_real64, &
        & 31623606042757.4_real64, -19231687611276.18_real64, 9390427534723.066_real64, -995725691772.4391_real64, &
        & 3231021975950.436_real64, -1236732873921.426_real64, -9315853881184.518_real64, 22185504639604.15_real64, - &
        & 15420649472134.23_real64, -50925469133365.87_real64, 153307611034700.9_real64, -173343833885846.3_real64, &
        & 49267726345354.81_real64, 70205032838989.57_real64, -69859638323604.89_real64, 43425460890416.86_real64, - &
        & 21893468750700.88_real64, 971026922359.6981_real64, 3092672108804.019_real64, -46104301322264.75_real64, &
        & 137071341025042.3_real64, -224229555487697.6_real64, 230058743397561.0_real64, -67834286156745.91_real64, - &
        & 128052812447042.3_real64, 158179514277751.0_real64, -59759767157347.25_real64, -50192959879342.73_real64, &
        & 58489878128799.11_real64, -37365652539010.6_real64, 18984302270208.09_real64, -4554741373320.545_real64, - &
        & 15788562576372.66_real64, 207582222421536.1_real64, -593029164486064.9_real64, 973326107456852.1_real64, - &
        & 1082367003969684.0_real64, 562615348697410.1_real64, 26458287552227.02_real64, -13583994639239.95_real64, - &
        & 392783298396344.4_real64, 1229601140956742.0_real64, -2067219533976532.0_real64, 2356855992627805.0_real64, &
        & -1246737057738835.0_real64, -68296541061671.73_real64, 128288499411321.2_real64, 289808443642052.9_real64, - &
        & 1183849512769141.0_real64, 1886550276050333.0_real64, -1788856965086742.0_real64, 750223740744196.8_real64, &
        & 73810969394052.8_real64, -177570076346560.2_real64, 2691625710632.994_real64, 381733861859218.0_real64, - &
        & 337040679079796.8_real64, -454514335566799.8_real64, 677319846538472.8_real64, -41563484188548.29_real64, &
        & 109261813626975.3_real64, -91264431397218.97_real64, 70019807698957.54_real64, -422305510260621.1_real64, &
        & 1201116601964427.0_real64, -998734904931376.2_real64, 6997153813544.085_real64, 7719885008754.813_real64, - &
        & 68889951937553.84_real64, 29510799198680.64_real64, 312693842722917.4_real64, -871687328070400.5_real64, &
        & 701036805957422.1_real64, 26790852901811.92_real64, -132704991301392.4_real64, 245407257700899.6_real64, - &
        & 134008336678529.4_real64, -247169099895962.2_real64, 703154494823089.5_real64, -552290648294985.4_real64, &
        & 57175999193.93517_real64, -217041704165.8141_real64, 351833929688.4115_real64, -283938217300.7131_real64, &
        & 128183343340.5703_real64, -71127907710.21759_real64, 23370341633.47879_real64, 15364808746.3056_real64, &
        & 290600970767.678_real64, -802466372151.1106_real64, 674732779512.1943_real64, -200832496932.1319_real64, &
        & 33390464521.60381_real64, 14119153398.94608_real64, -305873473110.3616_real64, 1104450618662.138_real64, - &
        & 1626906997911.81_real64, 1112583461902.734_real64, -344631785858.34_real64, 191565759438.5793_real64, - &
        & 70596621601.22795_real64, -47382282878.43769_real64, -1216306339687.263_real64, 2851854605294.669_real64, - &
        & 2088575877805.589_real64, 552683566033.1357_real64, -139592290375.9377_real64, -4916745273.507141_real64, &
        & 690249891008.9135_real64, -2248632451698.138_real64, 2783052411608.549_real64, -1232290482799.419_real64, - &
        & 322111044882.2086_real64, 325449902340.8409_real64, -200105226225.3467_real64, 45103837654.46812_real64, &
        & 2325737956409.067_real64, -4626800089844.668_real64, 3157955010334.468_real64, -1058647942564.135_real64, &
        & 658018158271.2852_real64, -309693138067.6807_real64, -608371015776.8824_real64, 1637625013204.523_real64, - &
        & 1159187962858.363_real64, -937769060969.854_real64, 2200931782687.973_real64, -1843542662764.745_real64, &
        & 1875745960635.015_real64, -1985230439532.362_real64, -562624466431.7668_real64, 4175799740028.683_real64, - &
        & 4518626496396.447_real64, 3013615771067.893_real64, -2592732254162.141_real64, 1538137723687.404_real64, - &
        & 80444716312.28485_real64, 799982275625.8677_real64, -2734974264689.439_real64, 4389536260691.129_real64, - &
        & 3984154969534.838_real64, 3376148592189.282_real64, -6570781832058.247_real64, 11693213254677.64_real64, - &
        & 10924136756001.26_real64, 1267416455919.102_real64, 7206599850869.638_real64, -7341315468904.064_real64, &
        & 5436453649841.863_real64, -2953644640925.872_real64, 649410078741.6183_real64, -2527659225727.932_real64, &
        & 4977795574849.692_real64, -5238675930088.641_real64, 1862012385633.843_real64, -186317142144.7943_real64, &
        & 9813127007417.85_real64, -27138975170185.2_real64, 32006809257626.92_real64, -13465389897611.32_real64, - &
        & 7793352836669.79_real64, 10037422806247.28_real64, -5178004616304.312_real64, 2063334948429.758_real64, - &
        & 450273075288.9221_real64, 619406753638.7299_real64, 3996253939028.438_real64, -16370002361155.56_real64, &
        & 31656768735320.77_real64, -36150707936908.06_real64, 11613915324187.33_real64, 21116286017153.73_real64, - &
        & 27563782911565.13_real64, 14047863430594.85_real64, 3127454385056.92_real64, -5338743057881.025_real64, &
        & 1868058723655.095_real64, -461070854710.0752_real64, -114668533063.3199_real64, 5990017067420.994_real64, - &
        & 38676945417161.76_real64, 98427629554136.33_real64, -156460933384874.5_real64, 174369788046456.5_real64, - &
        & 91611763999118.19_real64, -1475392128261.351_real64, -9767920920907.08_real64, 85723984395409.31_real64, - &
        & 217712512681944.4_real64, 337633613123190.6_real64, -373663008173175.4_real64, 196482251992050.5_real64, &
        & 9142767378760.65_real64, -9054881219099.158_real64, -76735619656431.84_real64, 226297819245529.7_real64, - &
        & 321404697362039.4_real64, 284616507075024.4_real64, -114724930174195.7_real64, -16653804709682.27_real64, &
        & 39470183010894.73_real64, 1275533756514.251_real64, -85254341320945.73_real64, 78491376640581.12_real64, &
        & 64481671862338.91_real64, -108103066748943.9_real64, 17849431100688.22_real64, -49760739431335.57_real64, &
        & 41652845934623.25_real64, -8126797264388.672_real64, 54797262106044.71_real64, -194399199573094.1_real64, &
        & 167256142536158.5_real64, -9241479350831.498_real64, 6882725690738.27_real64, 22875473887787.39_real64, - &
        & 2828764653546.728_real64, -109498039283190.8_real64, 260457417484636.6_real64, -202879311310509.8_real64, - &
        & 14633618021045.86_real64, 82624184757069.11_real64, -122714455031347.0_real64, -6145593087540.097_real64, &
        & 215150502674674.2_real64, -384410727382357.5_real64, 274034442850728.5_real64, -52609803496.65534_real64, &
        & 200916101615.5557_real64, -327477267759.7859_real64, 274692596849.5021_real64, -148178923753.5764_real64, &
        & 83452708108.71599_real64, 31645656784.26268_real64, -287810805460.3381_real64, 547094592527.1306_real64, - &
        & 545229456945.5374_real64, 266711547165.34_real64, -53690135875.41321_real64, 14809919736.99799_real64, - &
        & 4128060461.162109_real64, 286190382259.998_real64, -1045935846171.252_real64, 1586149246759.851_real64, - &
        & 1227937058680.117_real64, 631505033763.3636_real64, -390057856928.2278_real64, -121395923926.1331_real64, &
        & 1222421910875.021_real64, -2172489801435.53_real64, 2139390353502.051_real64, -1101505444804.173_real64, &
        & 221642597389.4062_real64, -33390007572.6416_real64, 3694512671.084229_real64, -677588577709.3079_real64, &
        & 2247279140926.674_real64, -2966020601774.421_real64, 1896060544221.933_real64, -744530501344.9713_real64, &
        & 261483759884.8765_real64, 1028154610187.762_real64, -2932014653362.273_real64, 3549608020752.532_real64, - &
        & 2657612845636.726_real64, 1086515301048.285_real64, 26014117942.67188_real64, -244698553664.3711_real64, &
        & 139820501131.9922_real64, 695493405168.0909_real64, -2023404277416.681_real64, 1987667826650.335_real64, - &
        & 399089505519.939_real64, -589690272189.1674_real64, 1319138883223.862_real64, -4180285400009.534_real64, &
        & 6825260037571.881_real64, -4838763450695.762_real64, 466904644221.3926_real64, 1561574567015.445_real64, - &
        & 1498815269501.688_real64, 1310349806150.625_real64, -793981665413.7031_real64, -152105598566.8847_real64, &
        & 252109054666.2753_real64, 761554344766.5514_real64, -2306094461457.561_real64, 2411458804296.284_real64, - &
        & 2916498691481.502_real64, 8641434306356.987_real64, -15075360338256.05_real64, 9885224606450.998_real64, &
        & 3364893069377.702_real64, -7664100456899.543_real64, 3998857859950.102_real64, -1928086964137.109_real64, &
        & 907316382169.0469_real64, -299646744732.7592_real64, 713451237927.2633_real64, -1370398995574.423_real64, &
        & 2163850800393.248_real64, -1192427868873.928_real64, 373807438604.2028_real64, -8230938215115.179_real64, &
        & 22716567677830.46_real64, -17970698297311.37_real64, -8858176656441.328_real64, 16308162045210.91_real64, - &
        & 4007742132208.854_real64, -1273014099457.242_real64, 1218183237802.648_real64, -114660546044.2902_real64, &
        & 2196513151556.065_real64, -7116631215709.257_real64, 11177042350300.56_real64, -13662738533127.57_real64, &
        & 14226055226525.23_real64, -3804284686063.335_real64, -14434393729529.85_real64, 13953311758047.04_real64, &
        & 7790540597853.452_real64, -12934928429613.03_real64, 395023662692.7401_real64, 3092479093104.284_real64, - &
        & 1350614134747.351_real64, 2350607156364.812_real64, -14517754863526.07_real64, 40221296075326.15_real64, - &
        & 59617675267552.41_real64, 63391880763625.62_real64, -58286142717249.55_real64, 28999362417675.74_real64, - &
        & 7497830156366.76_real64, 42936420741679.73_real64, -110841876362891.1_real64, 153504704671037.7_real64, - &
        & 151203495260445.9_real64, 128056685729966.9_real64, -59137391631851.63_real64, 9939346965606.754_real64, - &
        & 61595682464957.52_real64, 160887181901289.2_real64, -221676656809169.1_real64, 209718093082172.6_real64, - &
        & 155987439942106.9_real64, 60788349203900.27_real64, -977181328413.0098_real64, 22542844206721.41_real64, - &
        & 84255529200696.64_real64, 141848738928730.7_real64, -144671345783652.9_real64, 92189654092136.02_real64, - &
        & 26473687470219.96_real64, -9563519969712.414_real64, 29927512917633.12_real64, -29559010995164.79_real64, - &
        & 9871669556571.482_real64, 38758877972583.14_real64, -16064666592189.65_real64, -6424923487108.289_real64, &
        & 8682578177030.484_real64, -37256064188540.66_real64, 95726500459560.08_real64, -133996856335626.5_real64, &
        & 113049802312584.2_real64, -99268864201393.61_real64, 66534240039666.85_real64, 3993845722108.063_real64, &
        & 4183874124025.487_real64, -122251578672148.7_real64, 283838729639258.1_real64, -298775633069599.6_real64, &
        & 233879969400935.1_real64, -123916925391030.2_real64, 30314570629.11138_real64, -119295179405.7641_real64, &
        & 197536354233.7276_real64, -158362577626.4174_real64, 65626234117.36263_real64, 3826918165.30455_real64, - &
        & 153652751125.8799_real64, 525239070104.2095_real64, -1071568710132.397_real64, 1286898331748.261_real64, - &
        & 753368513583.3779_real64, 181952807564.5647_real64, -46869742951.09985_real64, 8995354806.075897_real64, - &
        & 174604349566.2355_real64, 649256151964.894_real64, -969789661316.2656_real64, 681422015716.0975_real64, - &
        & 238081483875.9602_real64, -32095856007.9952_real64, 706129754268.1013_real64, -2281018983798.294_real64, &
        & 4300725275106.133_real64, -4924486653698.362_real64, 2814539705574.197_real64, -647318785738.978_real64, &
        & 161652654695.3506_real64, -40302405092.4314_real64, 443577990001.1881_real64, -1481897273386.696_real64, &
        & 1824751597413.643_real64, -823885770071.7155_real64, -109619777288.6033_real64, 720692615453.7644_real64, - &
        & 2561876853020.526_real64, 5337952351428.719_real64, -7131527825300.779_real64, 6586693153110.304_real64, - &
        & 3342906118595.006_real64, 586778173577.8672_real64, -56997769451.95312_real64, 15253974475.01953_real64, - &
        & 536817043926.046_real64, 1613328453472.174_real64, -1380551948073.51_real64, -450120364499.0804_real64, &
        & 1560769445574.167_real64, -2748396857218.895_real64, 7152831790964.069_real64, -11145161661117.42_real64, &
        & 8596541422809.369_real64, -2693037632789.879_real64, -396397601241.0859_real64, 709631181717.7734_real64, - &
        & 462027529122.6719_real64, 259711943933.0938_real64, 319886802447.3705_real64, -1058053945486.439_real64, &
        & 495555087709.2324_real64, 1723923246700.757_real64, -2751126771396.97_real64, 4206700874179.441_real64, - &
        & 12685535042893.48_real64, 21699247868439.83_real64, -13309754408688.94_real64, -5504051325422.426_real64, &
        & 9457180240209.793_real64, -3231561023635.625_real64, 851758207211.75_real64, -348767993977.8086_real64, - &
        & 239948856960.2482_real64, 1483561517720.125_real64, -2428924624925.948_real64, 917745267299.5765_real64, &
        & 332288542893.8776_real64, -1030037301586.244_real64, 11093379021716.85_real64, -29400518172782.14_real64, &
        & 20501543908836.0_real64, 18323522773558.84_real64, -24411370687759.07_real64, 4250428333811.362_real64, &
        & 414112527444.333_real64, 226820172435.5056_real64, 849161363974.5648_real64, -5588302059806.737_real64, &
        & 13716150060875.45_real64, -16472078626731.11_real64, 13480047939847.37_real64, -10619000885659.76_real64, &
        & 1259478548615.425_real64, 17254467872463.24_real64, -15205377921885.13_real64, -16858841627838.2_real64, &
        & 21318945972006.64_real64, -815094327941.095_real64, -1518100952289.731_real64, -944421809101.1404_real64, - &
        & 3863872137506.954_real64, 22588878375784.41_real64, -55711277381727.01_real64, 68410789868593.34_real64, - &
        & 54367000902171.73_real64, 39139901935640.09_real64, -17685779379782.15_real64, 12996698120866.0_real64, - &
        & 69239096957575.61_real64, 157562134921832.8_real64, -181930750135613.9_real64, 136285843402314.3_real64, - &
        & 88785584395423.03_real64, 34804406428070.3_real64, -22398042952218.18_real64, 115660400902257.4_real64, - &
        & 253422882512886.8_real64, 292295646157892.1_real64, -230249239545776.2_real64, 156077910765844.8_real64, - &
        & 60774229527053.9_real64, 13883995480320.56_real64, -73132089879508.52_real64, 165055478244254.6_real64, - &
        & 211765953121747.7_real64, 200962393751350.0_real64, -161934125042153.7_real64, 74840393272622.09_real64, &
        & 1841341793601.803_real64, -7223532450069.731_real64, 7070891649828.395_real64, 26791810495060.59_real64, - &
        & 76927608430300.45_real64, 88973899802431.84_real64, -46311234911146.3_real64, -7358888090003.234_real64, &
        & 50015292614374.16_real64, -148267091963641.9_real64, 186422193717209.5_real64, -100483285399963.3_real64, &
        & 24391169137213.53_real64, -11095685653543.25_real64, 5244523447365.818_real64, -59775392601112.02_real64, &
        & 238160267885126.5_real64, -371422809864660.4_real64, 279774540232977.3_real64, -129267746106498.7_real64, &
        & 45671295089443.82_real64, -7976939400.566559_real64, 32136929665.2681_real64, -52688918668.49712_real64, &
        & 35672542416.32325_real64, 9823115.813047694_real64, -32974305663.58474_real64, 104263959933.3175_real64, - &
        & 285724098268.1221_real64, 571854645220.8142_real64, -695644674652.8407_real64, 412283798096.7196_real64, - &
        & 101051717488.2069_real64, 26803132734.47705_real64, -5805481324.323227_real64, 48927867377.61115_real64, - &
        & 183800323533.2747_real64, 257713722113.248_real64, -126317539503.0742_real64, -50130684059.91687_real64, &
        & 178311555758.5972_real64, -494325576029.993_real64, 1253888606431.941_real64, -2303794374988.417_real64, &
        & 2660272592095.93_real64, -1533946910101.71_real64, 361194702982.5962_real64, -97764422411.146_real64, &
        & 27738617173.0791_real64, -135555345975.261_real64, 452389857160.4962_real64, -472319553993.6606_real64, - &
        & 22197009799.43947_real64, 428993007460.2854_real64, -710449002406.5356_real64, 1581016027313.134_real64, - &
        & 2909660019094.34_real64, 3833728441804.47_real64, -3592320544900.736_real64, 1874892894986.929_real64, - &
        & 388825345045.5742_real64, 91333629788.04492_real64, -38560309452.31982_real64, 198630321718.0901_real64, - &
        & 605597870661.3627_real64, 391759392243.2587_real64, 594293384928.2036_real64, -1175628358923.416_real64, &
        & 1749615354317.604_real64, -3932588480638.138_real64, 5898086533672.911_real64, -4560366076139.509_real64, &
        & 1558242571872.873_real64, -5156617381.714844_real64, -159135911767.25_real64, 52716972752.60547_real64, - &
        & 28918121854.59961_real64, -197156842768.2396_real64, 688326557893.5311_real64, -456880513067.037_real64, - &
        & 835890984054.347_real64, 1550395047587.787_real64, -2349103209641.531_real64, 6643457271867.084_real64, - &
        & 11186681738936.05_real64, 6803945407092.148_real64, 2833905664456.155_real64, -4674352211756.564_real64, &
        & 1413669219213.473_real64, -311001978549.3398_real64, 147657961867.0605_real64, 271608780005.5567_real64, - &
        & 1237472847819.18_real64, 1881176234042.927_real64, -867795909727.9719_real64, -188024111393.1676_real64, &
        & 807797678203.4397_real64, -5807129753277.667_real64, 14858786287696.7_real64, -10210842744531.41_real64, - &
        & 9755752597791.68_real64, 12684254929910.12_real64, -2252524429531.041_real64, 376549444966.0415_real64, - &
        & 678782281016.0558_real64, -583340503463.1307_real64, 3295440286703.083_real64, -7507142478090.427_real64, &
        & 8445971432802.623_real64, -6124089355255.35_real64, 4158809104875.749_real64, 47395429010.48811_real64, - &
        & 8569703985865.056_real64, 7516611201944.66_real64, 9014453982603.84_real64, -11331827310160.56_real64, &
        & 679297278450.7517_real64, 60811602156.72156_real64, 1121287966432.533_real64, 2035814783828.713_real64, - &
        & 11732524211600.85_real64, 28241421711841.04_real64, -33183537751915.34_real64, 24126120717338.44_real64, - &
        & 15682011333804.08_real64, 6751045686938.137_real64, -6851897175683.5_real64, 36120452486285.05_real64, - &
        & 80241312908620.8_real64, 88749171537490.5_real64, -61229847511202.94_real64, 36202053092558.47_real64, - &
        & 13295127433995.6_real64, 12436260257500.07_real64, -62453301340026.09_real64, 132279805527287.1_real64, - &
        & 145908344139976.2_real64, 108233117260355.3_real64, -70663966031871.03_real64, 27418692142600.68_real64, - &
        & 8790642241483.315_real64, 43062465553366.98_real64, -90704942853499.7_real64, 108814552881446.6_real64, - &
        & 98836069467505.14_real64, 81213651858690.8_real64, -39275044955224.55_real64, 770667921133.3744_real64, - &
        & 2000098700372.355_real64, 3274339728395.479_real64, -16828889856514.36_real64, 40318214839641.7_real64, - &
        & 49021922458276.9_real64, 27167578256914.8_real64, 2742859356661.146_real64, -21719352563561.19_real64, &
        & 66637393713235.23_real64, -81315363098858.94_real64, 37404405860103.91_real64, -902044536910.0802_real64, - &
        & 355527859063.2392_real64, -3372814512784.159_real64, 31126028713408.41_real64, -109514259690364.8_real64, &
        & 158135018578786.2_real64, -108704891613649.7_real64, 42320296644470.91_real64, -13422759149157.41_real64, &
        & 101430531996079.0_real64, -326790538616101.0_real64, 225877491819537.3_real64, -242800321323902.9_real64, &
        & 792601009427518.6_real64, -550970830764485.6_real64, 149263265844044.3_real64, -494574917310713.9_real64, &
        & 345981988297228.8_real64, -256241577614737.2_real64, 829140072278390.8_real64, -574177528521221.2_real64, &
        & 553657992888829.4_real64, -1772614979566383.0_real64, 1221853663534118.0_real64, -308552103434606.6_real64, &
        & 980695832664743.5_real64, -673809253909547.0_real64, 162865375176750.8_real64, -530698355910966.1_real64, &
        & 368618278954474.5_real64, -324238768315837.1_real64, 1025131503115705.0_real64, -702672154434622.2_real64, &
        & 164221277784023.2_real64, -501315886027908.1_real64, 338118793290529.9_real64, 30898647469079.16_real64, - &
        & 66637610843496.03_real64, 35791220316520.74_real64, -74184870259699.64_real64, 160186918117288.5_real64, - &
        & 86117602123316.02_real64, 44739627068780.4_real64, -96679398517191.83_real64, 52005372453284.39_real64, - &
        & 60315170075073.81_real64, 129928566576669.4_real64, -69722691913374.2_real64, 144359261121873.3_real64, - &
        & 311194269610522.1_real64, 167086109234861.0_real64, -86936418530190.59_real64, 187561711307933.8_real64, - &
        & 100769231925143.4_real64, 28736029651098.9_real64, -61792795777504.94_real64, 33114088188433.97_real64, - &
        & 68437735879849.92_real64, 147140928358626.7_real64, -78842175056850.45_real64, 41093874577429.59_real64, - &
        & 88457217540699.14_real64, 47442444688919.08_real64, 18494480720370.5_real64, -35766062712334.0_real64, &
        & 17259757479457.0_real64, -37455977557880.0_real64, 72476861759284.0_real64, -34995868607951.0_real64, &
        & 18943527130447.0_real64, -36682145629326.0_real64, 17725468854635.0_real64, -44095296001046.0_real64, &
        & 85409033100108.0_real64, -41280255458022.0_real64, 84544166950460.0_real64, -164548224407368.0_real64, &
        & 79912509870388.0_real64, -39876213126728.0_real64, 78122406976276.0_real64, -38185127967802.0_real64, &
        & 26897867090768.0_real64, -52187880034600.0_real64, 25265871321295.0_real64, -48207673797510.0_real64, &
        & 94525324481944.0_real64, -46238872178816.0_real64, 20584227074990.0_real64, -41057670841148.0_real64, &
        & 20414797991306.0_real64 ]
end module chempotpy_hobr_3app_parameters
