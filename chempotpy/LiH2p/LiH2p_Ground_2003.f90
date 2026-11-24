! System:                        LiH2p
! Functional form:               many-body expansion with switched long range
! Common name:                   LiH2+ ground state
! Reference:                     R. Martinazzo et al., J. Chem. Phys. 119, 11241 (2003)
! Number of derivatives:         0
! Number of bodies:              3
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable parameters are embedded
! Notes:                         The atom order is H, H, Li. Energy is relative to three-body breakup.

!> @brief Evaluate the Martinazzo ground-state surface for LiH2+.
module chempotpy_lih2p_ground_2003_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    real(wp), parameter :: three_body_coefficients(143) = [ &
        59.9041369487_wp, -4.4190710855_wp, -1793.3690652080_wp, &
        -92.7339536285_wp, 231.0047518887_wp, 23.6312294844_wp, &
        42386.6962109884_wp, 20887.4910336156_wp, -101.8488429800_wp, &
        -4783.2820589398_wp, 5601.4299484618_wp, 647.9336802252_wp, &
        -2188.0964589478_wp, -333391.4339867483_wp, -50650.0836783211_wp, &
        -216249.6666120038_wp, -10968.2992516742_wp, 14746.3468957760_wp, &
        64623.8750484771_wp, -6585.4074148974_wp, -121946.1513253672_wp, &
        2470.5168826501_wp, 39388.8512400949_wp, 674663.0178284654_wp, &
        -623209.8463494780_wp, -625800.6376789842_wp, 5061315.2544497522_wp, &
        -499972.9861830364_wp, 738165.5689500275_wp, -3062084.0673502581_wp, &
        -890975.0918916499_wp, 37100.2541750794_wp, 481176.7795094966_wp, &
        68335.2125605610_wp, 860916.1659883368_wp, -182446.6634467566_wp, &
        -313469.8588537204_wp, 16296025.2245538831_wp, 24283257.8151341975_wp, &
        2226423.3841185221_wp, -55670137.2963673249_wp, -41861251.3156600073_wp, &
        12522078.8982328828_wp, -10780744.6234461721_wp, 19920454.5483050309_wp, &
        3725852.5719358511_wp, -888096.4455630139_wp, 5805044.5654820763_wp, &
        -1174676.0747827506_wp, -152672.8060593364_wp, -3582082.6546327011_wp, &
        252873.3936963925_wp, -3291637.8143871599_wp, 1409200.5075722407_wp, &
        1356914.6186123572_wp, -79156096.6236694306_wp, -134289948.7898918688_wp, &
        -144168137.8314403296_wp, 18010795.4982191287_wp, 184494908.8457103074_wp, &
        201704765.4667312205_wp, -99315107.5417414308_wp, 64030848.3463394046_wp, &
        62186512.0140229985_wp, -69323914.7188201100_wp, -18225109.7742339559_wp, &
        8767846.0950237475_wp, -71368402.0258482397_wp, 20651551.6130857319_wp, &
        1535688.4530311411_wp, 12342388.8021743763_wp, 5825471.3843415780_wp, &
        -3734791.5687767509_wp, 10453889.2599439640_wp, -704923.1248095465_wp, &
        7092463.0621291548_wp, -5229921.5270726196_wp, -3297509.8045103801_wp, &
        170414954.6685615480_wp, 442318060.4571669102_wp, 338310296.6504344344_wp, &
        -106295716.0848532170_wp, -243177154.1358072460_wp, -272160559.5405060053_wp, &
        -560001873.6517249346_wp, 345105141.7074456215_wp, -179453069.6148003042_wp, &
        -205181566.6644286513_wp, 178625484.6556275785_wp, 80765193.6117105931_wp, &
        -34816331.7087118998_wp, -22162585.9540813640_wp, 117201851.0785028338_wp, &
        -94646151.0332854241_wp, -10872024.4385482874_wp, 134794077.4815447628_wp, &
        -23336856.3831176758_wp, 15174814.7355045378_wp, -78927252.4897296727_wp, &
        -12871086.9000668079_wp, 4843936.1891660253_wp, -11668605.0056117158_wp, &
        426732.0261371761_wp, -8474699.3137837537_wp, 9738486.1953393687_wp, &
        4215023.2278999090_wp, 83709195.3309966177_wp, -81955642.6349996924_wp, &
        -509178802.8625016212_wp, -253164045.3674562573_wp, 145445180.5049311221_wp, &
        -261403144.5543840528_wp, -214500767.8142659366_wp, 597722827.9000399113_wp, &
        -439796210.5928990245_wp, 195491463.8727613091_wp, 412858125.1965236664_wp, &
        989693212.4453485012_wp, 179744133.8334875405_wp, -158573195.2249470055_wp, &
        47721187.7360739559_wp, -309115588.9086966515_wp, -865951061.6187733412_wp, &
        103449948.7291473746_wp, 25790340.4448547661_wp, 33253194.3017021790_wp, &
        415688777.1302601099_wp, 169521356.8767541051_wp, -23916729.3615326211_wp, &
        -214129802.4606905282_wp, -89171572.5754874349_wp, -11048394.8646354862_wp, &
        122655134.3360205442_wp, 37156079.7304429114_wp, 1487230.4260974966_wp, &
        -1885306.1160342137_wp, -502645.5527109082_wp, 5039172.7627611337_wp, &
        -7203742.8395571373_wp, -2200156.4276824687_wp, 1.0277126653_wp, &
        1.4179477712_wp, 2.5479621480_wp]
    integer, parameter :: first_exponent(140) = [ &
        0, 1, 0, 1, 1, 2, 0, 0, 1, 1, 2, 2, 3, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9]
    integer, parameter :: second_exponent(140) = [ &
        1, 1, 2, 1, 2, 1, 2, 3, 2, 3, 1, 2, 1, 3, 4, 2, 3, 4, 2, 3, 1, 2, 1, 3, 4, 5, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 4, 5, 6, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 4, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 5, 6, 7, 8, 4, 5, 6, 7, 8, 4, 5, 6, 7, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1, 5, 6, 7, 8, 9, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 4, 5, 6, 7, 3, 4, 5, 6, 3, 4, 5, 2, 3, 4, 2, 3, 1, 2, 1]
    integer, parameter :: third_exponent(140) = [ &
        1, 0, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 2, 1, 0, 1, 0, 1, 0, 0, 3, 2, 1, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 3, 2, 1, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 4, 3, 2, 1, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 4, 3, 2, 1, 4, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 5, 4, 3, 2, 1, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0]
    real(wp), parameter :: h2_quadrupole_coefficients(4) = [ &
        -0.0126271402_wp, +0.341987071_wp, -0.0249680149_wp, &
        +0.000251740527_wp]
    real(wp), parameter :: h2_polarizability_coefficients(4) = [ &
        +0.0682048466_wp, +1.69575813_wp, -0.0379060701_wp, &
        +0.0149673329_wp]
    real(wp), parameter :: h2_anisotropy_coefficients(5) = [ &
        0.0266251827_wp, 0.139081715_wp, 0.324064812_wp, &
        -0.0511107403_wp, 0.00279102422_wp]
    real(wp), parameter :: h2_potential_coefficients(8) = [ &
        -0.142978719_wp, -17.4496708_wp, 215.8658_wp, &
        -1763.38375_wp, 9338.55361_wp, -30639.411_wp, &
        56147.2305_wp, -43929.2685_wp]
    real(wp), parameter :: lih_potential_coefficients(6) = [ &
        -0.436998324_wp, -10.9156475_wp, 55.0379023_wp, &
        -437.23156_wp, 1443.14556_wp, -2340.82341_wp]

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
                'LiH2p_Ground_2003 provides energies but not analytic derivatives')
            return
        endif
        distances(1) = cartesian_distance(x(1,:), x(2,:))*angstrom_to_bohr
        distances(2) = cartesian_distance(x(2,:), x(3,:))*angstrom_to_bohr
        distances(3) = cartesian_distance(x(1,:), x(3,:))*angstrom_to_bohr
        p(1)         = ground_energy(distances)*hartree_to_electronvolt
    end subroutine evaluate

    !> @brief Evaluate the complete ground-state fit in atomic units.
    pure function ground_energy( distances ) result(energy)
        real(wp), intent(in) :: distances(3)
        real(wp)             :: energy
        energy = three_body_energy(distances(1), distances(2), distances(3)) + h2_potential(distances(1)) + &
            & lih_potential(distances(2)) + lih_potential(distances(3))
    end function ground_energy

    !> @brief Evaluate the permutationally adapted three-body and long-range terms.
    pure function three_body_energy( hh_distance, first_lih_distance, second_lih_distance ) result(energy)
        real(wp), intent(in) :: hh_distance, first_lih_distance, second_lih_distance
        real(wp)             :: energy
        real(wp)             :: angular_cosine, anisotropy, denominator, factor, long_range
        real(wp)             :: radial(3), scattering_distance, switched_radius, symmetric_term
        integer              :: term
        radial(1) = hh_distance*exp(-three_body_coefficients(141)*hh_distance)
        radial(2) = first_lih_distance*exp(-three_body_coefficients(142)*first_lih_distance)
        radial(3) = second_lih_distance*exp(-three_body_coefficients(142)*second_lih_distance)
        energy    = 0.0_wp
        do term = 1, 140
            if( second_exponent(term) == third_exponent(term) )then
                symmetric_term = radial(1)**first_exponent(term)*radial(2)**second_exponent(term)* radial(3)** &
                    & third_exponent(term)
            else
                symmetric_term = radial(1)**first_exponent(term)*( radial(2)**second_exponent(term)*radial(3)** &
                    & third_exponent(term) + radial(2)**third_exponent(term)*radial(3)**second_exponent(term))
            endif
            energy = energy+three_body_coefficients(term)*symmetric_term
        enddo
        scattering_distance = sqrt(0.5_wp*first_lih_distance**2 + 0.5_wp*second_lih_distance**2 - 0.25_wp* &
            & hh_distance**2)
        if( scattering_distance < abs(three_body_coefficients(143)) ) return
        angular_cosine  = (first_lih_distance**2-second_lih_distance**2)/ (2.0_wp*hh_distance*scattering_distance)
        switched_radius = abs(three_body_coefficients(143))/scattering_distance
        denominator     = 1.0_wp/(1.0_wp-switched_radius**2)
        factor          = exp(1.0_wp-denominator)
        anisotropy      = 0.5_wp*(3.0_wp*angular_cosine**2-1.0_wp)
        long_range = h2_quadrupole(hh_distance)*anisotropy/scattering_distance**3 - 0.5_wp* &
            & h2_polarizability(hh_distance)/scattering_distance**4 - 0.5_wp* &
            & h2_anisotropic_polarizability(hh_distance)*anisotropy/scattering_distance**4
        energy = energy*(1.0_wp-factor)+long_range*factor
    end function three_body_energy

    !> @brief Evaluate a polynomial in the damped radial coordinate.
    pure function damped_polynomial( distance, beta, coefficients ) result(value)
        real(wp), intent(in) :: distance, beta, coefficients(:)
        real(wp)             :: value, radial
        integer              :: order
        radial = distance*exp(-beta*distance)
        value  = 0.0_wp
        do order = 1, size(coefficients)
            value = value+coefficients(order)*radial**order
        enddo
    end function damped_polynomial

    !> @brief Evaluate the H2 quadrupole function in atomic units.
    pure function h2_quadrupole( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = damped_polynomial(distance, 0.0673276303_wp, h2_quadrupole_coefficients)/ (1.0_wp+exp(1.3899166_wp* &
            & (distance-3.30599025_wp)))
    end function h2_quadrupole

    !> @brief Evaluate the spherical three-body H2 polarizability.
    pure function h2_polarizability( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = (damped_polynomial(distance, 0.0567659808_wp, h2_polarizability_coefficients) - 7.4989518_wp*exp(- &
            & 0.0724081327_wp*distance))/ (1.0_wp+exp(1.55657199_wp*(distance-3.19559801_wp)))
    end function h2_polarizability

    !> @brief Evaluate the anisotropic H2 polarizability.
    pure function h2_anisotropic_polarizability( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = damped_polynomial(distance, -0.0729582362_wp, h2_anisotropy_coefficients)/ (1.0_wp+exp(1.59140975_wp* &
            & (distance-3.1416315_wp)))
    end function h2_anisotropic_polarizability

    !> @brief Evaluate the fitted H2 ground-state pair potential.
    pure function h2_potential( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        value = damped_polynomial(distance, 1.19999843_wp, h2_potential_coefficients) + 0.637669117_wp*exp(-1.2_wp* &
            & distance)/distance
    end function h2_potential

    !> @brief Evaluate the fitted LiH+ ground-state pair potential.
    pure function lih_potential( distance ) result(value)
        real(wp), intent(in) :: distance
        real(wp)             :: value
        real(wp)             :: factor, reduced_radius, short_range
        short_range = damped_polynomial(distance, 1.00000006_wp, lih_potential_coefficients) + 16.1279348_wp*exp(- &
            & distance)/distance
        reduced_radius = 5.85874466_wp/distance
        factor         = 0.0_wp
        if( reduced_radius < 1.0_wp ) factor = exp(1.0_wp-1.0_wp/(1.0_wp-reduced_radius**2))
        value = short_range*(1.0_wp-factor)-2.25_wp*factor/distance**4
    end function lih_potential
end module chempotpy_lih2p_ground_2003_pes

!> @brief Expose LiH2p_Ground_2003 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_kinds,                 only: wp
    use chempotpy_lih2p_ground_2003_pes, only: evaluate
    implicit none
    integer, parameter    :: natoms = 3, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
