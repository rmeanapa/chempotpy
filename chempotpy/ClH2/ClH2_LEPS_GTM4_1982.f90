!   System:          ClH2
!   Functional form: Empirical triatomic potential surface defined
!                    over orthogonal bond order coordinates
!   Common name: AL/AB
!   References:      B. C. Garrett, D. G. Truhlar, and A. W. Magnuson
!                    J. Chem. Phys. 76, 2321 (1982)
!                    N. Agmon and R. D. Levine
!                    J. Chem. Phys. 71, 3034 (1979)
!   Calling Sequence:
!               must be called once before any calls to POT
!      POT    - driver for the evaluation of the energy and the derivati
!               of the energy with respect to the coordinates for a give
!               geometric configuration
!   Units:
!      energies    - hartrees
!      coordinates - bohrs
!      derivatives - hartrees/bohr
!   Surface(s):
!      ground electronic state
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Cl
!      infinitely far from the H2 diatomic and R(H2) set equal to the
!      H2 equilibrium diatomic value.
!   Parameters:
!      Set in the BLOCK DATA subprogram PTPARM
!   Coordinates:
!      Internal, Definition: R(1) = R(Cl-first H)
!                            R(2) = R(first H-second H)
!                            R(3) = R(Cl-second H)
!   Common Blocks (used between the calling program and this potential):
!      /PT1CM/ R(3), ENERGY, DEDR(3)
!        passes the coordinates, ground state electronic energy, and
!        derivatives of the ground electronic state energy with respect
!        to the coordinates.
!      /PT2CM/ NSURF, NDER, NDUM(21)
!        passes the control flags where
!        NDER  = 0 => no derivatives are computed
!              = 1 => derivatives of the energy for the ground electroni
!                     state with respect to the coordinates are computed
!        NDUM  - not used
!      /PT4CM/ IPTPRT, IDUM(19)
!        IPTPRT passes the FORTRAN unit number used for potential output
!        IDUM   not used
!      /PT5CM/ EASYAB, EASYBC, EASYAC
!        passes the energy in the three asymptotic valleys for an A + BC
!        The energy in the AB valley, EASYAB, is equal to the energy of
!        C atom "infinitely" far from the AB diatomic and R(AB) set equa
!        Re(AB), the equilibrium bond length for the AB diatomic.
!        Similarly, the terms EASYBC and EASYAC represent the energies
!        in the BC and the AC valleys, respectively.
!   Default Parameter Values:
!      Variable      Default value
!      NSURF            0
!      NDER             1
!      IPTPRT           6

!> @brief Implement ClH2_LEPS_GTM4_1982, a potential energy surface for ClH2.
module chempotpy_clh2_leps_gtm4_1982_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the ClH2_LEPS_GTM4_1982 potential energy surface for ClH2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,              only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_units,              only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: triatomic_distances_to_bohr
        use chempotpy_kinds,              only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: der(3), energy, r(3), r2(3), tx(9), dedx(9), drdx(3,9), v
        integer               :: iatom, idir, istate, j
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is ClHH
        call triatomic_distances_to_bohr(x, r)
        der = 0.e0_wp
        call pot(r, igrad, energy, der)
        energy = energy*hartree_to_electronvolt
        der    = der*hartree_per_bohr_to_ev_per_angstrom
        r2     = r*bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(der, drdx, dedx)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = energy
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = energy
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    j               = (iatom-1)*3+idir
                    g(1,iatom,idir) = dedx(j)
                enddo
            enddo
        elseif( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
    end subroutine evaluate
end module chempotpy_clh2_leps_gtm4_1982_pes

!> @brief Expose ClH2_LEPS_GTM4_1982 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_clh2_leps_gtm4_1982_pes, only: evaluate
    use chempotpy_kinds,                   only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX for the ClH2_LEPS_GTM4_1982 potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    use chempotpy_kinds,    only: wp
    implicit none
    integer, parameter    :: first_atoms(3)  = [1, 2, 1]
    integer, parameter    :: second_atoms(3) = [2, 3, 3]
    real(wp), intent(in)  :: x(9), r(3)
    real(wp), intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

module chempotpy_clh2_leps_gtm4_1982_pot_legacy_state_69fc7d3b
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_internal

    contains

    !> @brief Evaluate dispatch.
    subroutine evaluate_internal( rab, rbc, rac, nder, energy, der )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), parameter   :: bohr   = .52917706e0_wp
        real(wp), parameter   :: hrtree = 627.5095e0_wp
        real(wp), intent(in)  :: rab, rbc, rac
        integer, intent(in)   :: nder
        real(wp), intent(out) :: energy, der(3)
        real(wp)              :: a, bcoor, beta(3), c, coef(2), dbcoor(2), de(3), denom, dew, dvstar(2), &
            & easyab, easyac, easybc, gamma, junk(4), lamda, m, morse(2), nab, nab1, nbc, nbc1, num, omega(2), &
            & rbb, req(3), rhh, rob, sigma, stuff, v1, v2, v3, v4, vabbc, vac, vstar
        real(wp), external :: bigm
        integer            :: it
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !      Echo the potential parameters
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !       CONVERT TO ATOMIC UNITS
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        a      = 0.e0_wp
        bcoor  = 0.e0_wp
        beta   = 0.e0_wp
        c      = 0.e0_wp
        coef   = 0.e0_wp
        dbcoor = 0.e0_wp
        de     = 0.e0_wp
        denom  = 0.e0_wp
        dew    = 0.e0_wp
        dvstar = 0.e0_wp
        easyab = 0.e0_wp
        easyac = 0.e0_wp
        easybc = 0.e0_wp
        gamma  = 0.e0_wp
        junk   = 0.e0_wp
        lamda  = 0.e0_wp
        m      = 0.e0_wp
        morse  = 0.e0_wp
        nab    = 0.e0_wp
        nab1   = 0.e0_wp
        nbc    = 0.e0_wp
        nbc1   = 0.e0_wp
        num    = 0.e0_wp
        omega  = 0.e0_wp
        rbb    = 0.e0_wp
        req    = 0.e0_wp
        rhh    = 0.e0_wp
        rob    = 0.e0_wp
        sigma  = 0.e0_wp
        stuff  = 0.e0_wp
        v1     = 0.e0_wp
        v2     = 0.e0_wp
        v3     = 0.e0_wp
        v4     = 0.e0_wp
        vabbc  = 0.e0_wp
        vac    = 0.e0_wp
        vstar  = 0.e0_wp
        de     = [106.447e0_wp, 109.458e0_wp, 106.447e0_wp]
        req    = [1.2732e0_wp, 0.74127e0_wp, 1.2732e0_wp]
        beta   = [1.8674e0_wp, 1.9413e0_wp, 1.8674e0_wp]
        lamda  = 8.2444e0_wp
        a      = 0.272e0_wp
        gamma  = 0.5e0_wp
        rbb    = 0.74127e0_wp
        do it = 1, 3
            req(it)  = req(it)/bohr
            beta(it) = beta(it) * bohr
            de(it)   = de(it)/hrtree
        enddo
        a     = a/bohr
        lamda = lamda/hrtree
        rhh   = 1.4007977e0_wp
        rbb   = rbb /bohr
        !   Set the values for the potenetial in the three asymptotic valleys
        easyab = de(1)
        easybc = de(2)
        easyac = de(3)
        !CCCCCCCCCCC
        !CCCCCCCCCCC
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        MOVE ONTO THE MAIN PART OF THE ROUTINE
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !   Check the value of NDER
        !        FIRST CALCULATE SOME USEFUL NUMBERS
        nab = exp((req(1) - rab)/a)
        if( nab .lt. 1.e-15_wp ) nab = 1.e-15_wp
        nbc = exp((req(2) - rbc)/a)
        if( nbc .lt. 1.e-15_wp ) nbc = 1.e-15_wp
        !        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
        !        AB AND BC RESPECTIVELY. (EQUATION 2.2)
        c = 1.0e0_wp-nab/nbc
        if( abs(c) <= 1.0e-14_wp )then
            nab1 = 0.5e0_wp
            nbc1 = 0.5e0_wp
        else
            c    = c/nab
            nab1 = (2.0e0_wp+c-sqrt(4.0e0_wp+c*c))/(2.0e0_wp*c)
            nbc1 = 1.0e0_wp-nab1
        endif
        !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
        !       SEE EQUATION 2.7A, 2.7B
        m = nab + nbc
        !        M IS THE SUM OF THE BOND ORDERS (EQUATION 2.3)
        omega(1) = nab / m
        omega(2) = nbc / m
        !       THE OMEGAS ARE USED IN EQUATION 3.15A
        sigma = -a * log(m)
        !        SIGMA IS THE MYSTERY NUMBER REFERENCED IN EQUATION 3.15B
        !        SE EQUATION 2.13
        do it = 1, 2
            junk(it)  = exp(beta(it)*(-sigma))
            morse(it) = de(it)*(junk(it)*junk(it)-2.e0_wp*junk(it))
        enddo
        !        MORSE FUNCTION--EQUATION 4-16, GAS PHASE REACTION
        !        RATE THEORY, H. S. JOHNSTON.
        !        ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
        !          FIRST TRAP OUT ANY ZERO ARGUEMENTS
        if( nab1*nbc1 > 0.0e0_wp )then
            stuff = a * log(nab1*nbc1)
            rob   = (rbb - stuff)/(rhh - stuff)
        else
            rob = 1.0e0_wp
        endif
        !CCCCCCCCCCCCCCCCCCCCCCCCCC
        !       CALCULATE DE(OMEGA)
        !CCCCCCCCCCCCCCCCCCCCCCCCCC
        !        THIS IS FORMULA 3.17
        dew = 0.0e0_wp
        do it = 1, 2
            dew = omega(it)*de(it) + dew
        enddo
        dew = dew - lamda * bigm(omega(1))
        !CCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        CALCULATE V*(SIGMA)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCC
        vstar = 0.e0_wp
        num   = 0.e0_wp
        denom = 0.e0_wp
        do it = 1, 2
            num   = omega(it) * morse(it) + num
            denom = omega(it)*de(it) + denom
        enddo
        vstar = num/denom
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        CALCULATE BENDING CORRECTION
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        CALCULATE V(RAC)
        junk(3) = exp(-beta(3)*(rac-req(3))/rob)
        vac     = gamma * de(3) * junk(3)* (1.e0_wp + 0.5e0_wp*junk(3))
        !        CALCULATE V(RAB+RBC)
        junk(4) = exp(beta(3)*(req(3)-rbc-rab)/rob)
        vabbc   = gamma * de(3) * junk(4) * (1.e0_wp + 0.5e0_wp*junk(4))
        !        HERE IS THE BENDING CORRECTION
        bcoor = vac - vabbc
        !          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
        !          OF BCOOR WITH RESPECT TO RAB  RBC AND RAC.
        !        CALCULATE DER(3) (THE DERIVATIVE WITH RESPECT TO RAC)
        if( nder .eq. 1 )then
            der(3)    = -gamma*de(3)*beta(3)*junk(3)*(1.e0_wp+junk(3))/rob
            dbcoor(1) = gamma * de(3) * beta(3) * junk(4) * (1.e0_wp+ junk(4))/rob
            dbcoor(2) = dbcoor(1)
        endif
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        NOW CALCULATE THE ENERGY
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        energy = dew * vstar + bcoor + de(2)
        !           DE(2) ADJUSTS THE POTENTIAL SO THAT IT IS ZERO
        !           IN THE ASYMPTOTIC REACTANTS REGION.
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        NOW CALCULATE THE DERIVATIVES, IF NDER = 1
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        goto_block_920_428: block
            if( nder .ne. 1 ) exit goto_block_920_428
            !        NEXT CALCULATE DER(1) AND DER(2)
            do it = 1, 2
                !        FIRST CALCULATE DVSTAR(1) AND DVSTAR(2) SEE EQUATION B.5
                v1 = omega(it)/(de(1)*omega(1) + de(2) * omega(2))
                if( it == 1 )then
                    v2 = omega(2)*(morse(2) - morse(1) - vstar*(de(2)-de(1)))/a
                else
                    v2 = omega(1)*(morse(1) - morse(2) - vstar * (de(1)-de(2)))/a
                endif
                v3 = -omega(it)*2.e0_wp*de(it)*beta(it)*(junk(it)*junk(it)-junk(it))
                if( it == 1 )then
                    v4 = -omega(2)*2*de(2)*beta(2)*(junk(2)*junk(2) - junk(2))
                else
                    v4 = -omega(1)*2.e0_wp*de(1)*beta(1)*(junk(1)*junk(1) - junk(1))
                endif
                dvstar(it) = v1*(v2+v3+v4)
                !       NOW CALCULATE THE COEFFICIENT OF VSTAR IN EQUATION B.1.
                coef(it) = (nab*nbc/(m*m*a)) * ((de(1)-de(2)+lamda * log(nab/nbc)) * ((-1.e0_wp)**it))
                !       FINALLY PUT IT ALL TOGETHER TO CALCULATE THE DERIVITIVE
                !       SEE EQUATION B.1
                der(it) = coef(it) * vstar + dew * dvstar(it) + dbcoor(it)
            enddo
        end block goto_block_920_428
    end subroutine evaluate_internal
end module chempotpy_clh2_leps_gtm4_1982_pot_legacy_state_69fc7d3b

!> @brief Evaluate GTM4 for one set of internal distances.
subroutine pot( r, nder, energy, der )
    use chempotpy_clh2_leps_gtm4_1982_pot_legacy_state_69fc7d3b, only: evaluate_internal
    use chempotpy_kinds,                                         only: wp
    implicit none
    real(wp), intent(in)  :: r(3)
    integer, intent(in)   :: nder
    real(wp), intent(out) :: energy, der(3)
    call evaluate_internal(r(1), r(2), r(3), nder, energy, der)
end subroutine pot

!> @brief Evaluate BIGM for the ClH2_LEPS_GTM4_1982 potential.
function bigm( dummy )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), intent(in) :: dummy
    real(wp)             :: bigm
    if( dummy >= 1.e0_wp .or. dummy <= 0.e0_wp )then
        bigm = 0.e0_wp
    else
        bigm = -dummy*log(dummy) - (1.e0_wp-dummy)*log(1.e0_wp-dummy)
    endif
end function bigm

!> @brief Evaluate POTINFO for the ClH2_LEPS_GTM4_1982 potential.
