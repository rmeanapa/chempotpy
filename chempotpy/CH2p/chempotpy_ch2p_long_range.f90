!> @brief Provide the asymptotic CH+ + H and C+ + H2 terms used by the CH2+ DPEM.
module chempotpy_ch2p_long_range
    use chempotpy_geometry, only: all_pair_distances
    use chempotpy_kinds,    only: wp
    implicit none
    private
    public :: ch2p_long_range

    contains

    !> @brief Select the smooth long-range channel and return its Cartesian derivatives.
    subroutine ch2p_long_range( pos, lr_shift, der )
        ! expected order of pos:
        ! H H C
        implicit none
        real(wp), intent(in)  :: pos(3*3)
        real(wp), intent(out) :: lr_shift(3), der(3,9)
        real(wp)              :: lr, chp(6), h(3), cp(3), h2(6), v_ff(3), dis(3)
        real(wp)              :: der_chp(6), der_h(3), dener_dpos(9)
        integer               :: imin(1)
        lr_shift = 0e0_wp
        v_ff     = 0e0_wp
        der      = 0e0_wp
        call all_pair_distances(pos, dis)
        ! Evaluate the CHp + H long range for each permutation
        ! H(1) C -- H(2)
        chp(1:3) = pos(7:9)
        chp(4:6) = pos(1:3)
        h(1:3)   = pos(4:6)
        v_ff(1)  = sigmoid_fn(dis(2)) ! dCH1
        ! H(2) C -- H(1)
        chp(1:3) = pos(7:9)
        chp(4:6) = pos(4:6)
        h(1:3)   = pos(1:3)
        v_ff(2)  = sigmoid_fn(dis(3)) ! dCH2
        ! H(1) H(2) -- C
        cp      = pos(7:9)
        h2      = pos(1:6)
        v_ff(3) = sigmoid_fn(dis(1)) ! dHH
        imin    = minloc(v_ff)
        if( imin(1) == 1 )then
            chp(1:3) = pos(7:9)
            chp(4:6) = pos(1:3)
            h(1:3)   = pos(4:6)
            call evaluate_chp_h_long_range(chp, h, lr, der_chp, der_h)
            lr_shift(1) = lr ! (sigma) CH+ + H
            der(1,7:9)  = der_chp(1:3)
            der(1,1:3)  = der_chp(4:6)
            der(1,4:6)  = der_h(1:3)
        elseif( imin(1) == 2 )then
            chp(1:3) = pos(7:9)
            chp(4:6) = pos(4:6)
            h(1:3)   = pos(1:3)
            call evaluate_chp_h_long_range(chp, h, lr, der_chp, der_h)
            lr_shift(1) = lr ! (sigma) CH+ + H
            der(1,7:9)  = der_chp(1:3)
            der(1,4:6)  = der_chp(4:6)
            der(1,1:3)  = der_h(1:3)
        elseif( imin(1) == 3 )then
            cp = pos(7:9)
            h2 = pos(1:6)
            call evaluate_cp_h2_long_range(cp, h2, lr, dener_dpos)
            lr_shift(1) = lr
            lr_shift(2) = lr
            der(1,7:9)  = dener_dpos(7:9)
            der(1,1:6)  = dener_dpos(1:6)
            der(2,7:9)  = dener_dpos(7:9)
            der(2,1:6)  = dener_dpos(1:6)
        endif
    end subroutine ch2p_long_range

    !> @brief Return the logistic switching value used to select an asymptotic channel.
    pure function sigmoid_fn( x ) result(y)
        implicit none
        real(wp), intent(in) :: x
        real(wp)             :: y
        y = 1e0_wp / (1e0_wp + exp(-x))
    end function sigmoid_fn

    !> @brief Evaluate the charge-induced-dipole CH+ + H long-range interaction.
    subroutine evaluate_chp_h_long_range( chp, h, ener, dener_dchp, dener_dh )
        implicit none
        real(wp), parameter :: alpha = 9e0_wp / 2e0_wp, &
            mass_chp(2) = [12.0_wp, 1.0_wp], &
            re = 2.5e0_wp, r0 = 100e0_wp ! Bohr
        real(wp), intent(in)  :: chp(6), h(3)
        real(wp), intent(out) :: ener, dener_dchp(6), dener_dh(3)
        integer               :: iat
        real(wp)              :: com(3), r_vec(3), rnorm, drnorm_drvec(3), drnorm_dh(3), drnorm_dchp(6)
        real(wp)              :: dedrtilde, drvec_dh(3), drvec_dchp(6), rtilde, drtilde_drnorm
        ! Initialize outputs
        ener       = 0e0_wp
        dener_dchp = 0e0_wp
        dener_dh   = 0e0_wp
        ! Compute center of mass (COM)
        com = 0e0_wp
        do iat = 1, 2
            com = com + chp(3 * (iat - 1) + 1 : 3 * iat) * mass_chp(iat)
        enddo
        com = com / sum(mass_chp(1:2))
        ! Compute R_vec and Rnorm
        r_vec  = h - com
        rnorm  = sqrt(sum(r_vec**2))
        rtilde = rnorm + r0 * exp(-(rnorm - re))
        ! Energy expression
        ener = -0.5e0_wp * alpha / rtilde**4
        ! Derivatives of Rnorm
        drnorm_drvec   = r_vec / rnorm
        drtilde_drnorm = 1e0_wp - r0 * exp(-(rnorm - re))
        ! Derivatives of R_vec with respect to H and CHp
        drvec_dh        = 1e0_wp
        drvec_dchp(1:3) = -mass_chp(1) / sum(mass_chp(1:2))
        drvec_dchp(4:6) = -mass_chp(2) / sum(mass_chp(1:2))
        ! Derivatives of Rnorm with respect to H and CHp
        drnorm_dh        = drnorm_drvec * drvec_dh
        drnorm_dchp(1:3) = drnorm_drvec * drvec_dchp(1:3)
        drnorm_dchp(4:6) = drnorm_drvec * drvec_dchp(4:6)
        ! Derivative of energy with respect to Rnorm
        dedrtilde = 2e0_wp * alpha / rtilde**5
        ! Derivatives of energy with respect to H and CHp
        dener_dh   = dedrtilde * drtilde_drnorm * drnorm_dh
        dener_dchp = dedrtilde * drtilde_drnorm * drnorm_dchp
    end subroutine evaluate_chp_h_long_range

    !> @brief Evaluate the multipolar C+ + H2 long-range interaction.
    subroutine evaluate_cp_h2_long_range( cp, h2, ener, dener_dpos )
        implicit none
        real(wp), parameter   :: re = 2.5e0_wp, r0 = 100.0e0_wp, charge = +1e0_wp
        real(wp), intent(in)  :: cp(3), h2(6)
        real(wp), intent(out) :: ener, dener_dpos(9)
        real(wp)              :: rcm(3), r(3), rhh(3), normrhh, normr, ctheta, h2quad, axx, azz, a0, a2, pol2, normrtilde
        real(wp)              :: dh2quad_drhh, daxx_drhh, dazx_drhh, dpol2_dctheta
        real(wp)              :: dnormrtilde_dr(3), dedrhh(3), dedr(3), dedctheta
        real(wp)              :: drhh_dh2(6), dr_dcp(3), dr_dh2(6), dnormr_dr(3)
        real(wp)              :: dctheta_dr(3), dctheta_drhh(3), dnormrhh_drhh(3)
        ! Initialize outputs
        ener       = 0e0_wp
        dener_dpos = 0e0_wp
        ! Compute intermediate quantities
        rhh        = h2(1:3) - h2(4:6)
        rcm        = (h2(1:3) + h2(4:6)) / 2e0_wp
        r          = cp - rcm
        normrhh    = sqrt(sum(rhh**2))
        normr      = sqrt(sum(r**2))
        normrtilde = normr + r0 * exp(-(normr - re))
        ctheta     = sum(r * rhh) / (normr * normrhh)
        if( abs(ctheta) > 1e0_wp ) ctheta = sign(1e0_wp, ctheta)
        ! Derivatives of intermediates w.r.t. inputs
        ! rHH derivatives w.r.t. H2 (drHH/dH2 = [I, -I])
        drhh_dh2(1:3) = 1e0_wp
        drhh_dh2(4:6) = -1e0_wp
        ! R derivatives: dR/dCp = I, dR/dH2 = -drcm/dH2 = [-0.5I, -0.5I]
        dr_dcp = 1e0_wp
        dr_dh2 = -0.5e0_wp
        ! normR derivatives: d(normR)/dR = R / normR
        dnormr_dr     = r / normr
        dnormrhh_drhh = rhh / normrhh
        ! normRtilde derivative: d(normRtilde)/dR = (1 - R0 * exp(-(normR - Re))) * dnormR_dR
        dnormrtilde_dr = (1e0_wp - r0 * exp(-(normr - re))) * dnormr_dr
        ! ctheta derivatives: d(ctheta)/dR = (rHH / (normR * normRHH)) - (ctheta * R / normR^2)
        dctheta_dr = (rhh / (normr * normrhh)) - (ctheta * r / normr**2)
        ! d(ctheta)/drHH = (R / (normR * normRHH)) - (ctheta * rHH / normRHH**2)
        dctheta_drhh = (r / (normr * normrhh)) - (ctheta * rhh / normrhh**2)
        ! Call subroutines with derivatives
        call evaluate_h2_quadrupole(normrhh, h2quad, dh2quad_drhh)
        call evaluate_h2_transverse_polarizability(normrhh, axx, daxx_drhh)
        call evaluate_h2_axial_polarizability(normrhh, azz, dazx_drhh)
        ! Compute a0, a2, pol2 and their derivatives
        a0            = (azz + 2e0_wp * axx) / 3e0_wp
        a2            = azz - axx
        pol2          = 0.5e0_wp * (3e0_wp * ctheta**2 - 1e0_wp)
        dpol2_dctheta = 3e0_wp * ctheta
        ! Energy contributions and their derivatives
        ! Charge-quadrupole term
        ener = ener + charge * h2quad * pol2 / normrtilde**3
        ! Charge-polarizability term
        ener = ener - (0.5e0_wp * a0 + a2 / 3e0_wp * pol2) / normrtilde**4
        ! Derivatives of energy components
        ! Charge-quadrupole term derivatives
        dedrhh    = charge * dh2quad_drhh * pol2 / normrtilde**3 * dnormrhh_drhh
        dedr      = charge * h2quad * pol2 * (-3e0_wp / normrtilde**4) * dnormrtilde_dr
        dedctheta = charge * h2quad * dpol2_dctheta / normrtilde**3
        ! Charge-polarizability term derivatives
        dedrhh = dedrhh - ( (dazx_drhh + 2e0_wp*daxx_drhh)/3e0_wp * 0.5e0_wp + (dazx_drhh - daxx_drhh)/3e0_wp * &
            & pol2 ) / normrtilde**4 * dnormrhh_drhh
        dedr      = dedr + (0.5e0_wp * a0 + a2 / 3e0_wp * pol2) * 4e0_wp / normrtilde**5 * dnormrtilde_dr
        dedctheta = dedctheta - (a2 / 3e0_wp * dpol2_dctheta) / normrtilde**4
        ! Combine derivatives using chain rule
        ! dE/dCp = dE/dR * dR/dCp + dE/dctheta * dctheta/dR * dR/dCp
        dener_dpos(7:9) = dedr * dr_dcp + dedctheta * dctheta_dr * dr_dcp
        ! dE/dH2 = dE/drHH * drHH/dH2 + dE/dR * dR/dH2 + dE/dctheta * (dctheta/drHH * drHH/dH2 + dctheta/dR * dR/dH2)
        dener_dpos(1:3) = dedrhh * drhh_dh2(1:3) + dedr * dr_dh2(1:3) + dedctheta * (dctheta_drhh * drhh_dh2(1:3) + &
            & dctheta_dr * dr_dh2(1:3))
        dener_dpos(4:6) = dedrhh * drhh_dh2(4:6) + dedr * dr_dh2(4:6) + dedctheta * (dctheta_drhh * drhh_dh2(4:6) + &
            & dctheta_dr * dr_dh2(4:6))
    end subroutine evaluate_cp_h2_long_range

    !> @brief Evaluate the released H2 quadrupole fit and its bond-length derivative.
    pure subroutine evaluate_h2_quadrupole( x, y, dydx )
        implicit none
        real(wp), parameter :: cf(15) = [ &
            0.8166648734890371e-04_wp, -0.1881855960082434e+01_wp, 0.2421277875344904e+00_wp, &
            0.6085863990932506e+02_wp, -0.2757874922103535e+03_wp, 0.1512576449669221e+04_wp, &
            -0.1090442521871086e+05_wp, 0.1905266714338982e+05_wp, 0.1730115053828064e+06_wp, &
            -0.1000012086016814e+07_wp, 0.1762147199545549e+07_wp, 0.9342476802465859e+06_wp, &
            -0.8105382035645692e+07_wp, 0.1151121162307790e+08_wp, -0.5602134540904779e+07_wp ]
        real(wp), parameter   :: vex(2) = [0.1340130881463356e+01_wp, 0.8322962713014708e+00_wp]
        real(wp), intent(in)  :: x
        real(wp), intent(out) :: y, dydx
        real(wp)              :: expo, vari, d_vari_dx, d_expo, new_expo, new_d_expo, aux
        integer               :: i
        ! Initialize y and dydx with the first two terms
        aux  = cf(2) * exp(-vex(1)*x)
        y    = cf(1) + aux
        dydx = -aux * vex(1)
        ! Variables for the recursive terms
        expo      = 1.e0_wp
        aux       = exp(-vex(2)*x)
        vari      = x * aux
        d_vari_dx = aux * (1.e0_wp - vex(2)*x)
        d_expo    = 0.e0_wp
        ! Loop through coefficients 3 to 15
        do i = 3, 15
            new_expo   = expo * vari
            new_d_expo = d_expo * vari + expo * d_vari_dx
            y          = y + cf(i) * new_expo
            dydx       = dydx + cf(i) * new_d_expo
            ! Update expo and d_expo for next iteration
            expo   = new_expo
            d_expo = new_d_expo
        enddo
    end subroutine evaluate_h2_quadrupole

    !> @brief Evaluate the released transverse H2 polarizability fit and derivative.
    pure subroutine evaluate_h2_transverse_polarizability( x, y, dydx )
        implicit none
        real(wp), parameter :: cf(15) = [ &
            0.8974583249592978e+01_wp, -0.1290835863481878e+02_wp, -0.5379421065843753e+01_wp, &
            0.1765706276501767e+03_wp, -0.3385803556921954e+04_wp, 0.4493699029122313e+05_wp, &
            -0.3652369817756390e+06_wp, 0.1852753697133512e+07_wp, -0.5809312606038393e+07_wp, &
            0.9721353024158247e+07_wp, -0.7092806937216849e+06_wp, -0.3327729880967365e+08_wp, &
            0.6983684025829618e+08_wp, -0.6381343977493817e+08_wp, 0.2300922761580030e+08_wp ]
        real(wp), parameter   :: vex(2) = [0.1187879729111910e+01_wp, 0.7554684825669591e+00_wp]
        real(wp), intent(in)  :: x
        real(wp), intent(out) :: y, dydx
        real(wp)              :: expo, vari, d_vari_dx, d_expo, new_expo, new_d_expo, aux
        integer               :: i
        ! Initialize y and dydx with the first two terms
        aux  = cf(2) * exp(-vex(1)*x)
        y    = cf(1) + aux
        dydx = -aux * vex(1)
        ! Variables for the recursive terms
        expo      = 1.e0_wp
        aux       = exp(-vex(2)*x)
        vari      = x * aux
        d_vari_dx = aux * (1.e0_wp - vex(2)*x)
        d_expo    = 0.e0_wp
        ! Loop through coefficients 3 to 15
        do i = 3, 15
            new_expo   = expo * vari
            new_d_expo = d_expo * vari + expo * d_vari_dx
            y          = y + cf(i) * new_expo
            dydx       = dydx + cf(i) * new_d_expo
            ! Update expo and d_expo for next iteration
            expo   = new_expo
            d_expo = new_d_expo
        enddo
    end subroutine evaluate_h2_transverse_polarizability

    !> @brief Evaluate the released axial H2 polarizability fit and derivative.
    pure subroutine evaluate_h2_axial_polarizability( x, y, dydx )
        implicit none
        real(wp), parameter :: cf(15) = [ &
            0.9005313464224852e+01_wp, -0.4087039493711530e+02_wp, 0.1169970996784607e+02_wp, &
            -0.2570966837001311e+03_wp, 0.1115822309638538e+05_wp, -0.1611895201468293e+06_wp, &
            0.1476263886597057e+07_wp, -0.8725387879754866e+07_wp, 0.3255107028148777e+08_wp, &
            -0.7292647566134660e+08_wp, 0.7756856658321184e+08_wp, 0.3625591185348335e+08_wp, &
            -0.2179510419880131e+09_wp, 0.2531262189573988e+09_wp, -0.1047041830161959e+09_wp ]
        real(wp), parameter   :: vex(2) = [0.1826054733720348e+01_wp, 0.7429159069146456e+00_wp]
        real(wp), intent(in)  :: x
        real(wp), intent(out) :: y, dydx
        real(wp)              :: expo, vari, d_vari_dx, d_expo, new_expo, new_d_expo, aux
        integer               :: i
        ! Initialize y and dydx with the first two terms
        aux  = cf(2) * exp(-vex(1)*x)
        y    = cf(1) + aux
        dydx = -aux * vex(1)
        ! Variables for the recursive terms
        expo      = 1.e0_wp
        aux       = exp(-vex(2)*x)
        vari      = x * aux
        d_vari_dx = aux * (1.e0_wp - vex(2)*x)
        d_expo    = 0.e0_wp
        ! Loop through coefficients 3 to 15
        do i = 3, 15
            new_expo   = expo * vari
            new_d_expo = d_expo * vari + expo * d_vari_dx
            y          = y + cf(i) * new_expo
            dydx       = dydx + cf(i) * new_d_expo
            ! Update expo and d_expo for next iteration
            expo   = new_expo
            d_expo = new_d_expo
        enddo
    end subroutine evaluate_h2_axial_polarizability
end module chempotpy_ch2p_long_range
