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
    implicit none
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    integer                       :: i, j
    ! Initialize dRdX(3,9)
    do i=1, 3
        do j=1, 9
            drdx(i,j) = 0.0d0
        enddo
    enddo
    drdx(1,1) = (x(1)-x(4))/r(1)
    drdx(1,2) = (x(2)-x(5))/r(1)
    drdx(1,3) = (x(3)-x(6))/r(1)
    drdx(1,4) = -drdx(1,1)
    drdx(1,5) = -drdx(1,2)
    drdx(1,6) = -drdx(1,3)
    drdx(2,4) = (x(4)-x(7))/r(2)
    drdx(2,5) = (x(5)-x(8))/r(2)
    drdx(2,6) = (x(6)-x(9))/r(2)
    drdx(2,7) = -drdx(2,4)
    drdx(2,8) = -drdx(2,5)
    drdx(2,9) = -drdx(2,6)
    drdx(3,1) = (x(1)-x(7))/r(3)
    drdx(3,2) = (x(2)-x(8))/r(3)
    drdx(3,3) = (x(3)-x(9))/r(3)
    drdx(3,7) = -drdx(3,1)
    drdx(3,8) = -drdx(3,2)
    drdx(3,9) = -drdx(3,3)
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
        do 50 it = 1, 3
            req(it)  = req(it)/bohr
            beta(it) = beta(it) * bohr
            de(it)   = de(it)/hrtree
        50 continue
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
        c = 1.d0 - nab/nbc
        if( abs(c) .gt. 1.d-14 ) goto 90
        nab1 = .5d0
        nbc1 = .5d0
        go to 95
        90 c = c/nab
        nab1 = (2.d0 + c - sqrt(4.d0+c*c))/(2.d0*c)
        nbc1 = 1.d0 - nab1
        !
        !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
        !       SEE EQUATION 2.7A, 2.7B
        !
        95 m = nab + nbc
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
        do 100 it = 1, 2
            junk(it)  = exp(beta(it)*(-sigma))
            morse(it) = de(it)*(junk(it)*junk(it)-2.d0*junk(it))
        100 continue
        !
        !        MORSE FUNCTION--EQUATION 4-16, GAS PHASE REACTION
        !        RATE THEORY, H. S. JOHNSTON.
        !
        !        ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
        !
        !          FIRST TRAP OUT ANY ZERO ARGUEMENTS
        !
        if( nab1*nbc1 .gt. 0.d0 ) go to 103
        rob = 1.d0
        go to 107
        !
        103 stuff = a * log(nab1*nbc1)
        rob = (rbb - stuff)/(rhh - stuff)
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCC
        !       CALCULATE DE(OMEGA)
        !CCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !        THIS IS FORMULA 3.17
        !
        107 dew = 0.d0
        do 110 it = 1, 2
            dew = omega(it)*de(it) + dew
        110 continue
        dew = dew - lamda * bigm(omega(1))
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        CALCULATE V*(SIGMA)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        vstar = 0.d0
        num   = 0.d0
        denom = 0.d0
        do 210 it = 1, 2
            num   = omega(it) * morse(it) + num
            denom = omega(it)*de(it) + denom
        210 continue
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
        if( nder .ne. 1 ) go to 920
        !
        !        NEXT CALCULATE DER(1) AND DER(2)
        !
        do 300 it = 1, 2
            !
            !        FIRST CALCULATE DVSTAR(1) AND DVSTAR(2) SEE EQUATION B.5
            !
            v1 = omega(it)/(de(1)*omega(1) + de(2) * omega(2))
            !
            if( it .eq. 1 ) go to 280
            v2 = omega(1)*(morse(1) - morse(2) - vstar * (de(1)-de(2)))/a
            go to 285
            280 v2 = omega(2)*(morse(2) - morse(1) - vstar*(de(2)-de(1)))/a
            285 continue
            !
            v3 = -omega(it)*2.d0*de(it)*beta(it)*(junk(it)*junk(it)-junk(it))
            !
            if( it .eq. 1 ) go to 290
            v4 = -omega(1)*2.d0*de(1)*beta(1)*(junk(1)*junk(1) - junk(1))
            go to 295
            290 v4 = -omega(2)*2*de(2)*beta(2)*(junk(2)*junk(2) - junk(2))
            295 continue
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
        300 continue
        310 continue
        920 continue
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
