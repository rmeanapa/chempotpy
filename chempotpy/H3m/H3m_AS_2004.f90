! System:                        H3m
! Functional form:               many-body expansion with asymptotic long-range correction
! Common name:                   AS H3- surface
! Reference:                     A. N. Panda and N. Sathyamurthy, J. Chem. Phys. 121, 9343 (2004)
! Number of derivatives:         0
! Number of bodies:              3
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable fit coefficients are embedded
! Notes:                         Repairs two undefined polarizability names in the released source.

!> @brief Evaluate the Panda-Sathyamurthy ground-state H3- potential energy surface.
module chempotpy_h3m_as_2004_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    real(wp), parameter :: diatomic_coefficients(10) = [ &
        & 0.613866431409e+00_wp, -0.190101798306e+01_wp, -0.125228951225e+02_wp, 0.462035129873e+03_wp, &
        & -0.719421334160e+04_wp, 0.675855488723e+05_wp, -0.402534224315e+06_wp, 0.148601194147e+07_wp, &
        & -0.310617990248e+07_wp, 0.281672666657e+07_wp ]
    real(wp), parameter :: three_body_coefficients(43) = [ &
        & -0.1112208493660213e+02_wp, 0.8397428790889362e+03_wp, 0.1094798235180020e+03_wp, -0.5709190211594033e+04_wp, &
        & -0.1298244769404030e+04_wp, -0.5972616424494424e+03_wp, 0.3340252211548539e+05_wp, 0.1928328064734260e+05_wp, &
        & 0.6458200480851449e+04_wp, 0.2435060288684594e+04_wp, -0.1434637341827002e+06_wp, -0.8873294052910528e+05_wp, &
        & -0.2560409973248071e+05_wp, -0.4155091370361293e+05_wp, -0.1908724741559524e+05_wp, -0.6916964417547845e+04_wp, &
        & 0.2700864073926629e+06_wp, 0.1761420315708862e+06_wp, 0.1285760798314841e+06_wp, 0.5241468442401338e+05_wp, &
        & 0.6179297668523817e+05_wp, 0.3613316158567466e+05_wp, 0.1167502561479660e+05_wp, -0.3483578214904613e+06_wp, &
        & -0.1987372527942004e+06_wp, -0.1655884747032285e+06_wp, -0.7617059629974826e+05_wp, -0.1202657194761244e+06_wp, &
        & -0.4853231398644997e+05_wp, -0.5201212017261994e+05_wp, -0.4169368997779521e+05_wp, -0.9856384447393651e+04_wp, &
        & 0.4031304431979369e+06_wp, 0.9745826736095028e+05_wp, 0.9724244007877808e+05_wp, 0.2355094235748628e+05_wp, &
        & 0.7589289221958217e+05_wp, 0.3505234464517287e+05_wp, 0.6811389394620169e+05_wp, 0.1432035524851620e+05_wp, &
        & 0.1201363089340074e+05_wp, 0.2214599648186910e+05_wp, 0.2937776762252248e+04_wp ]
    integer, parameter :: first_exponents(43) = [ &
        & 0, 1, 2, 2, 0, 3, 1, 3, &
        & 3, 4, 2, 3, 0, 4, 4, 5, &
        & 3, 1, 4, 4, 5, 5, 6, 2, &
        & 4, 4, 0, 5, 5, 6, 6, 7, &
        & 3, 4, 1, 5, 5, 5, 6, 6, &
        & 7, 7, 8 ]
    integer, parameter :: second_exponents(43) = [ &
        & 1, 1, 1, 1, 2, 1, 2, 1, &
        & 2, 1, 2, 2, 3, 1, 2, 1, &
        & 2, 3, 2, 3, 1, 2, 1, 3, &
        & 2, 3, 4, 2, 3, 1, 2, 1, &
        & 3, 3, 4, 2, 3, 4, 2, 3, &
        & 1, 2, 1 ]
    integer, parameter :: third_exponents(43) = [ &
        & 1, 1, 0, 1, 2, 0, 2, 1, &
        & 0, 0, 2, 1, 3, 1, 0, 0, &
        & 2, 3, 1, 0, 1, 0, 0, 3, &
        & 2, 1, 4, 1, 0, 1, 0, 0, &
        & 3, 2, 4, 2, 1, 0, 1, 0, &
        & 1, 0, 0 ]
    integer, parameter :: permutation_counts(43) = [ &
        & 3, 1, 6, 3, 3, 6, 3, 3, &
        & 6, 6, 1, 6, 3, 3, 6, 6, &
        & 3, 3, 6, 6, 3, 6, 6, 3, &
        & 3, 6, 3, 6, 6, 3, 6, 6, &
        & 1, 6, 3, 3, 6, 6, 6, 6, &
        & 3, 6, 6 ]

    contains

    !> @brief Evaluate the corrected, permutation-invariant H3m_AS_2004 energy in eV.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_geometry, only: fill_distance_matrix
        use chempotpy_status,   only: chempotpy_status_t, status_derivative_unavailable
        use chempotpy_units,    only: bohr_to_angstrom, hartree_to_electronvolt
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms = 3, nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: coordinates_bohr(natoms,3), distance_matrix(natoms,natoms)
        real(wp)              :: pair_distances(3), ordered_distances(3)
        integer               :: pair_index, remaining_index, shortest_index
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad == 0 )then
            coordinates_bohr = x/bohr_to_angstrom
            call fill_distance_matrix(coordinates_bohr, distance_matrix)
            pair_distances       = [distance_matrix(1,2), distance_matrix(2,3), distance_matrix(3,1)]
            shortest_index       = minloc(pair_distances, dim=1)
            ordered_distances(2) = pair_distances(shortest_index)
            remaining_index      = 0
            do pair_index=1, 3
                if( pair_index == shortest_index ) cycle
                remaining_index                        = remaining_index+1
                ordered_distances(2*remaining_index-1) = pair_distances(pair_index)
            enddo
            p(1) = corrected_as_energy(ordered_distances)*hartree_to_electronvolt
        elseif( present(status) )then
            call status%set(status_derivative_unavailable, &
                'H3m_AS_2004 provides energies but not analytic derivatives')
        endif
    end subroutine evaluate

    !> @brief Apply the repaired isotropic-polarizability and asymptotic switching expressions.
    function corrected_as_energy( distances ) result(energy_hartree)
        real(wp), intent(in) :: distances(3)
        real(wp)             :: energy_hartree
        real(wp)             :: alpha_isotropic, alpha_anisotropic, alpha_short, alpha_long
        real(wp)             :: base_three_body, center_distance, cosine_angle, legendre_p2
        real(wp)             :: quadrupole, short_distance, switch_argument, switch_weight
        real(wp)             :: long_range, long_range_induction, long_range_quadrupole
        short_distance  = distances(2)
        base_three_body = three_body_energy(distances)
        quadrupole = (-52.0714_wp*short_distance**2+21.6823_wp*short_distance-3.35413_wp+ 60.5523_wp*short_distance** &
            & 3-32.8347_wp*short_distance**4+8.11269_wp*short_distance**5)* exp(-2.09243_wp*short_distance)
        alpha_short = (11.0013_wp*short_distance**4-25.5122_wp*short_distance**3+ 37.9853_wp*short_distance**2- &
            & 15.3491_wp*short_distance+4.52484_wp)* exp(-1.25065_wp*short_distance)
        alpha_long = 9.14_wp+2.53_wp*exp(-2.6_wp*(short_distance-4.0_wp))
        if( short_distance <= 4.0_wp )then
            alpha_isotropic = alpha_short
        else
            alpha_isotropic = alpha_long
        endif
        alpha_anisotropic     = 0.67_wp*9.22521_wp*exp(-0.5022_wp*(short_distance-3.2_wp)**2)
        center_distance       = sqrt(0.5_wp*distances(1)**2+0.5_wp*distances(3)**2- 0.25_wp*short_distance**2)
        cosine_angle          = (distances(3)**2-distances(1)**2)/(2.0_wp*short_distance*center_distance)
        legendre_p2           = (3.0_wp*cosine_angle**2-1.0_wp)/2.0_wp
        long_range_induction  = -0.5_wp*(alpha_isotropic+alpha_anisotropic*legendre_p2)/center_distance**4
        long_range_quadrupole = -quadrupole*legendre_p2/center_distance**3
        long_range            = long_range_induction+long_range_quadrupole
        switch_argument       = 6.6_wp/center_distance
        if( switch_argument < 1.0_wp )then
            switch_weight   = exp(1.0_wp-1.0_wp/(1.0_wp-switch_argument))
            base_three_body = base_three_body*(1.0_wp-switch_weight)+long_range*switch_weight
        endif
        energy_hartree = base_three_body+diatomic_energy(distances(1))+ diatomic_energy(distances(2))+ &
            & diatomic_energy(distances(3))+0.17300816_wp
    end function corrected_as_energy

    !> @brief Evaluate the released ten-term H2 pair potential in hartree.
    function diatomic_energy( distance ) result(energy_hartree)
        real(wp), intent(in) :: distance
        real(wp)             :: energy_hartree
        real(wp)             :: basis_factor, basis_power
        integer              :: term
        energy_hartree = diatomic_coefficients(1)*exp(-distance)/distance
        basis_factor   = distance*exp(-1.499133_wp*distance)
        basis_power    = 1.0_wp
        do term=2, size(diatomic_coefficients)
            basis_power    = basis_power*basis_factor
            energy_hartree = energy_hartree+diatomic_coefficients(term)*basis_power
        enddo
    end function diatomic_energy

    !> @brief Evaluate the released 43-term symmetric three-body correction in hartree.
    function three_body_energy( distances ) result(energy_hartree)
        real(wp), intent(in) :: distances(3)
        real(wp)             :: energy_hartree
        real(wp)             :: powers(3,0:8), term_value
        integer, parameter   :: permutation_order(6,3) = reshape([ &
            1, 1, 2, 2, 3, 3, &
            2, 3, 1, 3, 1, 2, &
            3, 2, 3, 1, 2, 1], [6, 3])
        integer :: exponent_1, exponent_2, exponent_3, permutation_index, term
        powers(:,0) = 1.0_wp
        do exponent_1=1, 8
            powers(:,exponent_1) = powers(:,exponent_1-1)*distances*exp(-0.9100010000000002_wp*distances)
        enddo
        energy_hartree = 0.0_wp
        do term=1, size(three_body_coefficients)
            exponent_1 = first_exponents(term)
            exponent_2 = second_exponents(term)
            exponent_3 = third_exponents(term)
            term_value = 0.0_wp
            select case( permutation_counts(term) )
            case(1)
                term_value = powers(1,exponent_1)*powers(2,exponent_2)*powers(3,exponent_3)
            case(3)
                term_value = powers(1,exponent_1)*powers(2,exponent_2)*powers(3,exponent_3)+ powers(1,exponent_3)* &
                    & powers(2,exponent_1)*powers(3,exponent_2)+ powers(1,exponent_2)*powers(2,exponent_3)* powers(3, &
                    & exponent_1)
            case(6)
                term_value = sum([ (powers(permutation_order(permutation_index,1),exponent_1)* &
                    & powers(permutation_order(permutation_index,2),exponent_2)* &
                    & powers(permutation_order(permutation_index,3),exponent_3), permutation_index=1, 6) ])
            end select
            energy_hartree = energy_hartree+three_body_coefficients(term)*term_value
        enddo
    end function three_body_energy
end module chempotpy_h3m_as_2004_pes

!> @brief Expose H3m_AS_2004 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3m_as_2004_pes, only: evaluate
    use chempotpy_kinds,           only: wp
    implicit none
    integer, parameter    :: natoms = 3, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
