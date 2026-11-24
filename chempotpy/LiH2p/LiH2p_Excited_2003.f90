! System:                        LiH2p
! Functional form:               many-body expansion with switched long range
! Common name:                   LiH2+ first excited state
! Reference:                     R. Martinazzo et al., J. Chem. Phys. 119, 11241 (2003)
! Number of derivatives:         0
! Number of bodies:              3
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable parameters are embedded
! Notes:                         The atom order is H, H, Li. Energy is relative to three-body breakup.

!> @brief Evaluate the Martinazzo first-excited-state surface for LiH2+.
module chempotpy_lih2p_excited_2003_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    real(wp), parameter :: three_body_coefficients(230) = [ &
        7.0139563146_wp, 3.0003108967_wp, -76.4867842937_wp, &
        581.0183707108_wp, -447.5847281554_wp, -99.4864290612_wp, &
        982.7398204780_wp, 362.5648698799_wp, -3393.2414511771_wp, &
        5178.1173784156_wp, -7863.7073601031_wp, 6560.3607376437_wp, &
        1113.8316789772_wp, -5017.3030972528_wp, -440.8591101941_wp, &
        25489.4179184216_wp, 11173.7792402553_wp, -35732.6711370254_wp, &
        30244.5599087393_wp, -48056.3359768260_wp, 77877.3998465644_wp, &
        -71808.6974729424_wp, 5764.8295676300_wp, 1006.7053378493_wp, &
        28332.5137592468_wp, -8092.2103234528_wp, -79498.1299204307_wp, &
        -31293.4097147527_wp, 164370.3230362090_wp, -191847.1898188904_wp, &
        -85926.7361277200_wp, 260275.4135666166_wp, -133673.0490531631_wp, &
        263387.6756272943_wp, -632639.9064772174_wp, 603041.5077466281_wp, &
        -236922.6937759927_wp, 6700.8245146343_wp, -113425.9637954208_wp, &
        60455.4406802392_wp, 378520.5104870343_wp, 86984.0196323107_wp, &
        55343.8072617274_wp, -503235.4170471951_wp, 285008.8945165113_wp, &
        516410.1142055998_wp, -1024788.4006496501_wp, 1354991.4595010185_wp, &
        -235933.4600355194_wp, -868190.4855862955_wp, 603123.2438464863_wp, &
        -1252224.6800955555_wp, 3314794.4823525730_wp, -3330758.8231016411_wp, &
        2355595.5164276436_wp, -207196.8778448466_wp, 101681.1108576819_wp, &
        257460.6576209367_wp, -201212.0670082429_wp, -580128.3865388532_wp, &
        -126806.8597396127_wp, 54164.8740416092_wp, 1010250.0667563326_wp, &
        -2188435.0117557826_wp, 797679.9871993139_wp, -1793697.7770774949_wp, &
        2611739.3747180039_wp, -883251.7685729633_wp, -1349145.2436505784_wp, &
        3009615.9249107814_wp, -7052322.0718244938_wp, 3874587.9290700238_wp, &
        1225199.8923960077_wp, -3507124.1024165265_wp, 5796231.3624979686_wp, &
        -8891846.5984294061_wp, 10832850.4398638811_wp, -12425853.1764493100_wp, &
        301736.6636589296_wp, -371127.7362507661_wp, -321645.7644776623_wp, &
        374599.8369838075_wp, 1973036.8288541059_wp, -225429.4068309222_wp, &
        581127.4474386526_wp, -459020.9465607482_wp, -1293190.5723435762_wp, &
        1565452.0461066423_wp, -1254345.0027045240_wp, 2769830.0552521981_wp, &
        -4183276.3368961909_wp, 9865817.5040407330_wp, -5687710.6781170880_wp, &
        6331461.7343169162_wp, -6199060.3498003324_wp, 2366805.4481070125_wp, &
        -295291.9032663772_wp, -4626557.7481302628_wp, 21978300.2186226137_wp, &
        -14309379.2452538181_wp, 96083.8516814366_wp, 14169688.6588801872_wp, &
        -20227745.2251465656_wp, 4905125.0786394291_wp, -17851102.8897421658_wp, &
        38924717.4709330797_wp, -576834.7877907237_wp, -46335.9803832678_wp, &
        504037.2236824127_wp, 202416.2680335442_wp, -401805.4777432226_wp, &
        -1465800.8280503147_wp, 738207.4260055233_wp, -986097.7081596345_wp, &
        840428.3375779393_wp, 991264.5669337681_wp, -5995665.7003214313_wp, &
        3238025.1562320245_wp, -198540.6010630105_wp, -2406737.5763311395_wp, &
        4214995.0096759945_wp, -2530526.3486393159_wp, 276406.6675839963_wp, &
        -4399887.9279605504_wp, 6257808.5376846725_wp, -38630207.8119526505_wp, &
        29587282.2978549488_wp, -17989117.3686880358_wp, 10206444.9309776053_wp, &
        -7956727.4088536510_wp, 11935639.2603963744_wp, 272587.5750498392_wp, &
        -33205003.3066417538_wp, 17814420.1847252585_wp, 1717949.1537269219_wp, &
        -29524642.9344762526_wp, 41664281.4065521657_wp, 30649831.9779503755_wp, &
        6411040.1861357680_wp, -72601162.1412622780_wp, 301573.8474337890_wp, &
        -151845.9726148361_wp, -312104.3259274737_wp, -41842.7503088380_wp, &
        232165.1746365801_wp, 2329864.8484232700_wp, -431076.9872573870_wp, &
        -81290.3065144457_wp, 668919.9198792472_wp, -695227.8774063820_wp, &
        -391985.4736881022_wp, -617857.7470953518_wp, 868434.0566439229_wp, &
        -1795080.1182449430_wp, 2099682.7952767937_wp, -2570819.2477941490_wp, &
        39547802.8107521236_wp, -26430084.7782766521_wp, 14040722.8520985711_wp, &
        -2630006.2582660136_wp, -2890426.7765964856_wp, -4285618.4932269054_wp, &
        -1483471.0828852160_wp, 9488176.4738190025_wp, -7094951.3043547636_wp, &
        89999375.1135660559_wp, -53602764.0800236166_wp, 16091618.2000474222_wp, &
        -6629563.5433165785_wp, -1253521.0871043403_wp, -11254990.7965670228_wp, &
        5098017.7156343600_wp, 34072707.4960278347_wp, -2872237.4304915103_wp, &
        -11445974.2646609601_wp, 25413934.4136321582_wp, -42223432.0668093711_wp, &
        -68299818.6303414106_wp, 19227087.4643655121_wp, 74483220.7264803648_wp, &
        -179309.5475248325_wp, 17617.4766226398_wp, 53230.0044518733_wp, &
        78078.6615530894_wp, -7811.8766089654_wp, -55943.5555708998_wp, &
        -267671.7329121784_wp, -117273.5346655442_wp, 160907.6128860498_wp, &
        -260296.6368125519_wp, 243850.4724256268_wp, 51209.2539915440_wp, &
        -21857101.1706669591_wp, 18320151.9379566535_wp, -9529829.7924672943_wp, &
        3690041.8095493489_wp, -1307964.2923558105_wp, 749871.0760208279_wp, &
        1844403.7695033511_wp, -6647893.7860741401_wp, 2730696.9130423600_wp, &
        -112201.7731048101_wp, 689349.4466145103_wp, -124924487.2901448756_wp, &
        99219617.4259525836_wp, -41801056.6200024039_wp, 6218141.7051764801_wp, &
        1148692.0085307148_wp, 8811748.7360878177_wp, -19204502.2225882821_wp, &
        7627010.5606003199_wp, 1193908.8264275370_wp, -96724132.1153392941_wp, &
        64398757.3842907548_wp, -15602214.6525493693_wp, 859044.0599557476_wp, &
        8829806.1955996510_wp, -4666527.9020765796_wp, -379501.8213393019_wp, &
        -24132122.5725693628_wp, -2759943.9894857286_wp, 10162478.5972349979_wp, &
        -4359214.8902296340_wp, 15410672.6411291305_wp, 43276836.9157600775_wp, &
        -19009630.1279357336_wp, -32327109.3021873273_wp, 0.8545339729_wp, &
        0.5385533051_wp, 8.6596046408_wp]
    integer, parameter :: first_exponent(227) = [ &
        0, 1, 0, 1, 1, 2, 0, 0, 1, 1, 2, 2, 3, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 10, 10, 11]
    integer, parameter :: second_exponent(227) = [ &
        1, 1, 2, 1, 2, 1, 2, 3, 2, 3, 1, 2, 1, 3, 4, 2, 3, 4, 2, 3, 1, 2, 1, 3, 4, 5, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 4, 5, 6, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 4, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 5, 6, 7, 8, 4, 5, 6, 7, 8, 4, 5, 6, 7, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 5, 6, 7, 8, 9, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 4, 5, 6, 7, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 6, 7, 8, 9, 10, 5, 6, 7, 8, 9, 10, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 4, 5, 6, 7, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 6, 7, 8, 9, 10, 11, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 4, 5, 6, 7, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1]
    integer, parameter :: third_exponent(227) = [ &
        1, 0, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 2, 1, 0, 1, 0, 1, 0, 0, 3, 2, 1, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 3, 2, 1, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 4, 3, 2, 1, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 4, 3, 2, 1, 4, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 5, 4, 3, 2, 1, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 5, 4, 3, 2, 1, 5, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 6, 5, 4, 3, 2, 1, 5, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0]
    real(wp), parameter :: dipole_coefficients(5) = [ &
        3.81107729_wp, -3.07621218_wp, 1.20552107_wp, &
        -0.213831318_wp, 0.0142348131_wp]
    real(wp), parameter :: quadrupole_coefficients(6) = [ &
        -4.05169074_wp, 7.010268_wp, -4.53005171_wp, &
        1.53830649_wp, -0.26118265_wp, 0.0176316218_wp]
    real(wp), parameter :: polarizability_coefficients(5) = [ &
        -137.958769_wp, 68.7249887_wp, -14.5416019_wp, &
        1.24211167_wp, -0.00608506371_wp]
    real(wp), parameter :: anisotropy_coefficients(5) = [ &
        -3.10727044_wp, -4.16298439_wp, 3.20124828_wp, &
        -0.724437507_wp, 0.0614349836_wp]
    real(wp), parameter :: h2p_coefficients(8) = [ &
        -0.433932717_wp, -1.46625756_wp, 26.1872609_wp, &
        -227.903359_wp, 1135.07033_wp, -3281.68774_wp, &
        5086.2199_wp, -3273.43527_wp]
    real(wp), parameter :: lih_coefficients(8) = [ &
        -10.7239499_wp, -117.330715_wp, -1700.19362_wp, &
        26575.4401_wp, -273771.612_wp, 1500872.94_wp, &
        -4367963.76_wp, 5024043.4_wp]
    real(wp), parameter :: lihp_coefficients(6) = [ &
        -0.161505456_wp, 0.151370422_wp, 1.09938176_wp, &
        -3.30073329_wp, 5.16501195_wp, -2.98494189_wp]

    contains

    !> @brief Evaluate the energy in eV for Cartesian coordinates in angstrom.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_status,   only: chempotpy_status_t, status_derivative_unavailable
        use chempotpy_units,    only: angstrom_to_bohr, hartree_to_electronvolt
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms = 3, nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
        real(wp)              :: distances(3)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad /= 0 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'LiH2p_Excited_2003 provides energies but not analytic derivatives')
            return
        endif
        distances(1) = cartesian_distance(x(1,:), x(2,:))*angstrom_to_bohr
        distances(2) = cartesian_distance(x(2,:), x(3,:))*angstrom_to_bohr
        distances(3) = cartesian_distance(x(1,:), x(3,:))*angstrom_to_bohr
        p(1)         = excited_energy(distances)*hartree_to_electronvolt
    end subroutine evaluate

    !> @brief Evaluate the complete excited-state fit in atomic units.
    pure function excited_energy( distances ) result(energy)
        real(wp), intent(in) :: distances(3)
        real(wp)             :: energy, shorter_lih, longer_lih
        shorter_lih = min(distances(2), distances(3))
        longer_lih  = max(distances(2), distances(3))
        energy = three_body_energy(distances(1), distances(2), distances(3)) + h2_ion_potential(distances(1)) + &
            & lih_potential(shorter_lih) + lih_ion_excited_potential(longer_lih)
    end function excited_energy

    !> @brief Evaluate the permutationally adapted three-body and long-range terms.
    pure function three_body_energy( hh_distance, first_lih_distance, second_lih_distance ) result(energy)
        real(wp), intent(in) :: hh_distance, first_lih_distance, second_lih_distance
        real(wp)             :: energy
        real(wp)             :: angular_cosine, denominator, factor, long_range, radial(3)
        real(wp)             :: scattering_distance, shorter_lih, longer_lih, switched_radius, symmetric_term
        integer              :: term
        radial(1) = hh_distance*exp(-three_body_coefficients(228)*hh_distance)
        radial(2) = first_lih_distance*exp(-three_body_coefficients(229)*first_lih_distance)
        radial(3) = second_lih_distance*exp(-three_body_coefficients(229)*second_lih_distance)
        energy    = 0.0_wp
        do term = 1, 227
            if( second_exponent(term) == third_exponent(term) )then
                symmetric_term = radial(1)**first_exponent(term)*radial(2)**second_exponent(term)* radial(3)** &
                    & third_exponent(term)
            else
                symmetric_term = radial(1)**first_exponent(term)*( radial(2)**second_exponent(term)*radial(3)** &
                    & third_exponent(term) + radial(2)**third_exponent(term)*radial(3)**second_exponent(term))
            endif
            energy = energy+three_body_coefficients(term)*symmetric_term
        enddo
        shorter_lih         = min(first_lih_distance, second_lih_distance)
        longer_lih          = max(first_lih_distance, second_lih_distance)
        scattering_distance = sqrt(0.5_wp*longer_lih**2 + 0.5_wp*hh_distance**2 - 0.25_wp*shorter_lih**2)
        if( scattering_distance <= abs(three_body_coefficients(230)) ) return
        angular_cosine  = (longer_lih**2-hh_distance**2)/(2.0_wp*shorter_lih*scattering_distance)
        switched_radius = abs(three_body_coefficients(230))/scattering_distance
        denominator     = 1.0_wp/(1.0_wp-switched_radius**4)
        factor          = exp(1.0_wp-denominator)
        long_range = -lih_dipole(shorter_lih)*angular_cosine/scattering_distance**2 + lih_quadrupole(shorter_lih)* &
            & legendre_two(angular_cosine)/scattering_distance**3 - 0.5_wp*lih_polarizability(shorter_lih)/ &
            & scattering_distance**4 - 0.5_wp*lih_anisotropic_polarizability(shorter_lih)* &
            & legendre_two(angular_cosine)/ scattering_distance**4
        energy = energy*(1.0_wp-factor)+long_range*factor
    end function three_body_energy

    !> @brief Return the second Legendre polynomial.
    pure function legendre_two( cosine ) result(value)
        real(wp), intent(in) :: cosine
        real(wp)             :: value
        value = 0.5_wp*(3.0_wp*cosine**2-1.0_wp)
    end function legendre_two

    !> @brief Evaluate a polynomial in the damped radial coordinate.
    pure function damped_polynomial( distance, beta, coefficients ) result(value)
        real(wp), intent(in) :: distance, beta, coefficients(:)
        real(wp)             :: value, radial
        integer              :: order
        radial = distance*exp(-beta*distance); value = 0.0_wp
        do order = 1, size(coefficients)
            value = value+coefficients(order)*radial**order
        enddo
    end function damped_polynomial

    !> @brief Evaluate a Fermi-switched damped polynomial.
    pure function switched_polynomial( distance, beta, center, steepness, coefficients ) result(value)
        real(wp), intent(in) :: distance, beta, center, steepness, coefficients(:)
        real(wp)             :: value
        value = damped_polynomial(distance, beta, coefficients)/ (1.0_wp+exp(steepness*(distance-center)))
    end function switched_polynomial

    !> @brief Evaluate the LiH dipole fit.
    pure function lih_dipole( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = switched_polynomial(distance, 0.0307255411_wp, 6.68701348_wp, 1.664664_wp, dipole_coefficients)
    end function lih_dipole

    !> @brief Evaluate the LiH quadrupole fit.
    pure function lih_quadrupole( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = switched_polynomial(distance, 0.0381326063_wp, 6.54214083_wp, 1.72807586_wp, quadrupole_coefficients)
    end function lih_quadrupole

    !> @brief Evaluate the spherical LiH three-body polarizability.
    pure function lih_polarizability( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = (damped_polynomial(distance, 0.00608598433_wp, polarizability_coefficients) - 52.7688091_wp* &
            & exp(0.0040014362_wp*distance))/ (1.0_wp+exp(1.64572959_wp*(distance-6.74611063_wp)))
    end function lih_polarizability

    !> @brief Evaluate the anisotropic LiH polarizability.
    pure function lih_anisotropic_polarizability( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = switched_polynomial(distance, -0.0274272462_wp, 6.69233524_wp, 1.77848435_wp, anisotropy_coefficients)
    end function lih_anisotropic_polarizability

    !> @brief Evaluate a smooth switch from a pair fit to an inverse-fourth-power tail.
    pure function ion_pair_potential( distance, beta, exponential_beta, exponential_scale, switch_radius, &
            switch_power, long_range_scale, coefficients ) result(value)
        real(wp), intent(in) :: distance, beta, exponential_beta, exponential_scale, switch_radius
        integer, intent(in)  :: switch_power
        real(wp), intent(in) :: long_range_scale, coefficients(:)
        real(wp)             :: value, factor, reduced_radius, short_range
        short_range = damped_polynomial(distance, beta, coefficients) + exponential_scale*exp(-exponential_beta* &
            & distance)/distance
        reduced_radius = switch_radius/distance; factor = 0.0_wp
        if( reduced_radius < 1.0_wp ) factor = exp(1.0_wp-1.0_wp/(1.0_wp-reduced_radius**switch_power))
        value = short_range*(1.0_wp-factor)-long_range_scale*factor/distance**4
    end function ion_pair_potential

    !> @brief Evaluate the H2+ ground-state pair potential.
    pure function h2_ion_potential( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = ion_pair_potential(distance, 0.900001396_wp, 0.899999999_wp, 0.616429616_wp, 10.9526825_wp, 8, &
            & 2.25_wp, h2p_coefficients)
    end function h2_ion_potential

    !> @brief Evaluate the neutral LiH ground-state pair potential.
    pure function lih_potential( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = damped_polynomial(distance, 1.2_wp, lih_coefficients) + 216.039947_wp*exp(-1.2_wp*distance)/distance
    end function lih_potential

    !> @brief Evaluate the excited LiH+ pair potential.
    pure function lih_ion_excited_potential( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = ion_pair_potential(distance, 0.479976824_wp, 0.480098406_wp, 0.976997582_wp, 8.81582567_wp, 2, &
            & 82.5_wp, lihp_coefficients)
    end function lih_ion_excited_potential
end module chempotpy_lih2p_excited_2003_pes

!> @brief Expose LiH2p_Excited_2003 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_kinds,                  only: wp
    use chempotpy_lih2p_excited_2003_pes, only: evaluate
    implicit none
    integer, parameter    :: natoms = 3, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
