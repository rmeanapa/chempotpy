!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_clh2_leps_gtm4_1982_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_legacy_units,       only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: legacy_triatomic_distances_to_bohr
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: der(3), energy, r(3), r2(3), tx(9), dedx(9), drdx(3,9), v
        integer                       :: iatom, idir, istate, j
        !initialize
        v = 0.d0; g = 0.d0; d = 0.d0
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is ClHH
        call legacy_triatomic_distances_to_bohr(x, r)
        der = 0.d0
        call pot(r, igrad, energy, der)
        energy = energy*legacy_hartree_to_electronvolt
        der    = der*51.422067
        r2     = r*legacy_bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(der, drdx, dedx)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = energy
            enddo
        else if( igrad == 1 )then
            do istate=1, nstates
                p(istate) = energy
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    j               = (iatom-1)*3+idir
                    g(1,iatom,idir) = dedx(j)
                enddo
            enddo
        else if( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
    end subroutine evaluate
end module chempotpy_clh2_leps_gtm4_1982_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_clh2_leps_gtm4_1982_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX for the ClH2_LEPS_GTM4_1982 potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    implicit none
    integer, parameter            :: first_atoms(3)  = [1, 2, 1]
    integer, parameter            :: second_atoms(3) = [2, 3, 3]
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!
module chempotpy_clh2_leps_gtm4_1982_pot_legacy_state_69fc7d3b
    implicit none
    private
    public :: evaluate_internal
    contains
    !> @brief Evaluate dispatch.
    subroutine evaluate_internal( rab, rbc, rac, nder, energy, der )
        !
        !   System:          ClH2
        !   Functional form: Empirical triatomic potential surface defined
        !                    over orthogonal bond order coordinates
        !   Common name: AL/AB
        !   References:      B. C. Garrett, D. G. Truhlar, and A. W. Magnuson
        !                    J. Chem. Phys. 76, 2321 (1982)
        !                    N. Agmon and R. D. Levine
        !                    J. Chem. Phys. 71, 3034 (1979)
        !
        !   Calling Sequence:
        !               must be called once before any calls to POT
        !      POT    - driver for the evaluation of the energy and the derivati
        !               of the energy with respect to the coordinates for a give
        !               geometric configuration
        !
        !   Units:
        !      energies    - hartrees
        !      coordinates - bohrs
        !      derivatives - hartrees/bohr
        !
        !   Surface(s):
        !      ground electronic state
        !
        !   Zero of energy:
        !      The classical potential energy is set equal to zero for the Cl
        !      infinitely far from the H2 diatomic and R(H2) set equal to the
        !      H2 equilibrium diatomic value.
        !
        !   Parameters:
        !      Set in the BLOCK DATA subprogram PTPARM
        !
        !   Coordinates:
        !      Internal, Definition: R(1) = R(Cl-first H)
        !                            R(2) = R(first H-second H)
        !                            R(3) = R(Cl-second H)
        !
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
        !
        !   Default Parameter Values:
        !      Variable      Default value
        !      NSURF            0
        !      NDER             1
        !      IPTPRT           6
        !
        implicit none
        double precision, parameter   :: bohr   = .52917706d0
        double precision, parameter   :: hrtree = 627.5095d0
        double precision, intent(in)  :: rab, rbc, rac
        integer, intent(in)           :: nder
        double precision, intent(out) :: energy, der(3)
        double precision              :: a, bcoor, beta(3), c, coef(2), dbcoor(2), de(3), denom, dew, dvstar(2), &
            & easyab, easyac, easybc, gamma, junk(4), lamda, m, morse(2), nab, nab1, nbc, nbc1, num, omega(2), &
            & rbb, req(3), rhh, rob, sigma, stuff, v1, v2, v3, v4, vabbc, vac, vstar
        double precision, external :: bigm
        integer                    :: it
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !      Echo the potential parameters
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !WRITE (IPTPRT,600) DE, REQ, BETA
        !WRITE (IPTPRT,610) A, GAMMA, RBB, LAMDA
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !       CONVERT TO ATOMIC UNITS
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        a      = 0.d0
        bcoor  = 0.d0
        beta   = 0.d0
        c      = 0.d0
        coef   = 0.d0
        dbcoor = 0.d0
        de     = 0.d0
        denom  = 0.d0
        dew    = 0.d0
        dvstar = 0.d0
        easyab = 0.d0
        easyac = 0.d0
        easybc = 0.d0
        gamma  = 0.d0
        junk   = 0.d0
        lamda  = 0.d0
        m      = 0.d0
        morse  = 0.d0
        nab    = 0.d0
        nab1   = 0.d0
        nbc    = 0.d0
        nbc1   = 0.d0
        num    = 0.d0
        omega  = 0.d0
        rbb    = 0.d0
        req    = 0.d0
        rhh    = 0.d0
        rob    = 0.d0
        sigma  = 0.d0
        stuff  = 0.d0
        v1     = 0.d0
        v2     = 0.d0
        v3     = 0.d0
        v4     = 0.d0
        vabbc  = 0.d0
        vac    = 0.d0
        vstar  = 0.d0
        de     = [106.447d0, 109.458d0, 106.447d0]
        req    = [1.2732d0, 0.74127d0, 1.2732d0]
        beta   = [1.8674d0, 1.9413d0, 1.8674d0]
        lamda  = 8.2444d0
        a      = 0.272d0
        gamma  = 0.5d0
        rbb    = 0.74127d0
        do it = 1, 3
            req(it)  = req(it)/bohr
            beta(it) = beta(it) * bohr
            de(it)   = de(it)/hrtree
        enddo
        a     = a/bohr
        lamda = lamda/hrtree
        rhh   = 1.4007977d0
        rbb   = rbb /bohr
        !
        !   Set the values for the potenetial in the three asymptotic valleys
        easyab = de(1)
        easybc = de(2)
        easyac = de(3)
        !
        !CCCCCCCCCCC
        !CCCCCCCCCCC
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        MOVE ONTO THE MAIN PART OF THE ROUTINE
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !
        !   Check the value of NDER
        !
        !        FIRST CALCULATE SOME USEFUL NUMBERS
        !
        nab = exp((req(1) - rab)/a)
        !
        if( nab .lt. 1.d-15 ) nab = 1.d-15
        !
        nbc = exp((req(2) - rbc)/a)
        !
        if( nbc .lt. 1.d-15 ) nbc = 1.d-15
        !
        !        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
        !        AB AND BC RESPECTIVELY. (EQUATION 2.2)
        !
        c = 1.0d0-nab/nbc
        if( abs(c) <= 1.0d-14 )then
            nab1 = 0.5d0
            nbc1 = 0.5d0
        else
            c    = c/nab
            nab1 = (2.0d0+c-sqrt(4.0d0+c*c))/(2.0d0*c)
            nbc1 = 1.0d0-nab1
        endif
        !
        !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
        !       SEE EQUATION 2.7A, 2.7B
        !
        m = nab + nbc
        !
        !        M IS THE SUM OF THE BOND ORDERS (EQUATION 2.3)
        !
        omega(1) = nab / m
        omega(2) = nbc / m
        !
        !       THE OMEGAS ARE USED IN EQUATION 3.15A
        !
        sigma = -a * log(m)
        !
        !        SIGMA IS THE MYSTERY NUMBER REFERENCED IN EQUATION 3.15B
        !        SE EQUATION 2.13
        !
        do it = 1, 2
            junk(it)  = exp(beta(it)*(-sigma))
            morse(it) = de(it)*(junk(it)*junk(it)-2.d0*junk(it))
        enddo
        !
        !        MORSE FUNCTION--EQUATION 4-16, GAS PHASE REACTION
        !        RATE THEORY, H. S. JOHNSTON.
        !
        !        ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
        !
        !          FIRST TRAP OUT ANY ZERO ARGUEMENTS
        !
        if( nab1*nbc1 > 0.0d0 )then
            stuff = a * log(nab1*nbc1)
            rob   = (rbb - stuff)/(rhh - stuff)
        else
            rob = 1.0d0
        endif
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCC
        !       CALCULATE DE(OMEGA)
        !CCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !        THIS IS FORMULA 3.17
        !
        dew = 0.0d0
        do it = 1, 2
            dew = omega(it)*de(it) + dew
        enddo
        dew = dew - lamda * bigm(omega(1))
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        CALCULATE V*(SIGMA)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        vstar = 0.d0
        num   = 0.d0
        denom = 0.d0
        do it = 1, 2
            num   = omega(it) * morse(it) + num
            denom = omega(it)*de(it) + denom
        enddo
        vstar = num/denom
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        CALCULATE BENDING CORRECTION
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !        CALCULATE V(RAC)
        !
        junk(3) = exp(-beta(3)*(rac-req(3))/rob)
        vac     = gamma * de(3) * junk(3)* (1.d0 + 0.5d0*junk(3))
        !
        !        CALCULATE V(RAB+RBC)
        !
        junk(4) = exp(beta(3)*(req(3)-rbc-rab)/rob)
        vabbc   = gamma * de(3) * junk(4) * (1.d0 + 0.5d0*junk(4))
        !
        !        HERE IS THE BENDING CORRECTION
        !
        bcoor = vac - vabbc
        !
        !          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
        !          OF BCOOR WITH RESPECT TO RAB  RBC AND RAC.
        !
        !        CALCULATE DER(3) (THE DERIVATIVE WITH RESPECT TO RAC)
        !
        if( nder .eq. 1 )then
            der(3)    = -gamma*de(3)*beta(3)*junk(3)*(1.d0+junk(3))/rob
            dbcoor(1) = gamma * de(3) * beta(3) * junk(4) * (1.d0+ junk(4))/rob
            dbcoor(2) = dbcoor(1)
        endif
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        NOW CALCULATE THE ENERGY
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        energy = dew * vstar + bcoor + de(2)
        !
        !           DE(2) ADJUSTS THE POTENTIAL SO THAT IT IS ZERO
        !           IN THE ASYMPTOTIC REACTANTS REGION.
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        NOW CALCULATE THE DERIVATIVES, IF NDER = 1
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        goto_block_920_428: block
            if( nder .ne. 1 ) exit goto_block_920_428
            !
            !        NEXT CALCULATE DER(1) AND DER(2)
            !
            do it = 1, 2
                !
                !        FIRST CALCULATE DVSTAR(1) AND DVSTAR(2) SEE EQUATION B.5
                !
                v1 = omega(it)/(de(1)*omega(1) + de(2) * omega(2))
                !
                if( it == 1 )then
                    v2 = omega(2)*(morse(2) - morse(1) - vstar*(de(2)-de(1)))/a
                else
                    v2 = omega(1)*(morse(1) - morse(2) - vstar * (de(1)-de(2)))/a
                endif
                !
                v3 = -omega(it)*2.d0*de(it)*beta(it)*(junk(it)*junk(it)-junk(it))
                !
                if( it == 1 )then
                    v4 = -omega(2)*2*de(2)*beta(2)*(junk(2)*junk(2) - junk(2))
                else
                    v4 = -omega(1)*2.d0*de(1)*beta(1)*(junk(1)*junk(1) - junk(1))
                endif
                !
                dvstar(it) = v1*(v2+v3+v4)
                !
                !       NOW CALCULATE THE COEFFICIENT OF VSTAR IN EQUATION B.1.
                !
                coef(it) = (nab*nbc/(m*m*a)) * ((de(1)-de(2)+lamda * log(nab/nbc)) * ((-1.d0)**it))
                !
                !       FINALLY PUT IT ALL TOGETHER TO CALCULATE THE DERIVITIVE
                !       SEE EQUATION B.1
                !
                der(it) = coef(it) * vstar + dew * dvstar(it) + dbcoor(it)
                !
            enddo
        end block goto_block_920_428
        !
        !
    end subroutine evaluate_internal
end module chempotpy_clh2_leps_gtm4_1982_pot_legacy_state_69fc7d3b

!> @brief Evaluate GTM4 for one set of internal distances.
subroutine pot( r, nder, energy, der )
    use chempotpy_clh2_leps_gtm4_1982_pot_legacy_state_69fc7d3b, only: evaluate_internal
    implicit none
    double precision, intent(in)  :: r(3)
    integer, intent(in)           :: nder
    double precision, intent(out) :: energy, der(3)
    call evaluate_internal(r(1), r(2), r(3), nder, energy, der)
end subroutine pot

!
!> @brief Evaluate BIGM for the ClH2_LEPS_GTM4_1982 potential.
function bigm( dummy )
    implicit none
    double precision, intent(in) :: dummy
    double precision             :: bigm
    if( dummy >= 1.d0 .or. dummy <= 0.d0 )then
        bigm = 0.d0
    else
        bigm = -dummy*log(dummy) - (1.d0-dummy)*log(1.d0-dummy)
    endif
end function bigm

!
!> @brief Evaluate POTINFO for the ClH2_LEPS_GTM4_1982 potential.
