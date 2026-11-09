!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_o2h_gen_pes
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
        double precision              :: der(3), energy, r(3), v
        integer                       :: iatom, idir, istate, j
        double precision              :: tx(9), r2(3), dedx(9), drdx(3,9)
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
end module chempotpy_o2h_gen_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_o2h_gen_pes, only: evaluate
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

!> @brief Evaluate EvdRdX for the O2H_GEN potential.
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
module chempotpy_o2h_gen_pot_legacy_state_fffab58a
    implicit none
    private
    public :: evaluate_internal
    contains
    !> @brief Evaluate dispatch.
    subroutine evaluate_internal( r, nder, energy, der )
        !
        !   System: HO2
        !   Reference: C. F. Melius and R. J. Blint
        !              Chem. Phys. Lett. 64, 183 (1979).
        !
        !   The potential parameters are included in DATA statements.
        !   and in the block data subprogram PPARMS.
        !
        !   Common Blocks (used between the calling program and this potential):
        !      /PT1CM/ R(3), ENERGY, DEDR(3)
        !        passes the coordinates, ground state electronic energy, and
        !        derivatives of the ground electronic state energy with respect
        !        to the coordinates.
        !      /PT2CM/ NSURF, NDER, NDUM(21)
        !        passes the control flags where
        !        NSURF = 0    ground electronic surface
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
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !         GENERAL INFORMATION
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !          THE POTENTIAL CONSISTS OF THREE MORSE TERMS FOR
        !          THE THREE DIATOMICS, HO,OO, AND HO RESPECTIVELY
        !          PLUS TWO INTERACTION TERMS.
        !
        !          EACH POTENTIAL TERM WILL BE CALCULATED INDIVIDUALLY
        !          TO SIMPLIFY READING THE CODE.
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          DEFINE VARIABLES, ARRAYS, AND COMMON BLOCK
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        implicit none
        intent(in)  :: r, nder
        intent(out) :: energy, der
        dimension r(3), der(3)
        double precision :: a, alpha, b, beta, c, csth, da, db, dcsth, de, der, dmew, dprod1, dprod2, &
            & dprod3, dsec, dthird, dv, easyab, easyac, easybc, energy, prod1, prod2, prod3, r, rab, rac, &
            & rbc, req, second, third, vab, vac, vbc, zeroad, zev, zkcal
        integer :: i, idum, iprt, it, it1, j, nder, ndum, nsurf
        !
        double precision :: num, denom, gamma, junk, expo, dexpo, mew
        !
        !
        !
        dimension dcsth(3), dv(3), gamma(3), junk(4), c(5,3,3), a(5), da(5,3), &
            & req(3), de(3), expo(3), dexpo(3), dmew(3), b(5, 3), db(5, 3, 3), dprod1(3 &
            & ), dprod2(3), dsec(3), dthird(3), dprod3(3)
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          Initialize the constants
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !   The array DE contains the dissociation energies in hartree atomic un
        !   The array GAMMA contains the Morse Betas in reciprocal bohrs.
        !   The array REQ contains the equilibrium bond lengths in bohr.
        !   The variable ZEROAD contains the constant added to the energy in har
        !   atomic units.
        !
        !          SET ALPHA, UNITS: RECIPROCAL BOHR
        a      = 0.d0
        alpha  = 0.d0
        b      = 0.d0
        beta   = 0.d0
        c      = 0.d0
        csth   = 0.d0
        da     = 0.d0
        db     = 0.d0
        dcsth  = 0.d0
        de     = 0.d0
        denom  = 0.d0
        dexpo  = 0.d0
        dmew   = 0.d0
        dprod1 = 0.d0
        dprod2 = 0.d0
        dprod3 = 0.d0
        dsec   = 0.d0
        dthird = 0.d0
        dv     = 0.d0
        easyab = 0.d0
        easyac = 0.d0
        easybc = 0.d0
        expo   = 0.d0
        gamma  = 0.d0
        junk   = 0.d0
        mew    = 0.d0
        num    = 0.d0
        prod1  = 0.d0
        prod2  = 0.d0
        prod3  = 0.d0
        rab    = 0.d0
        rac    = 0.d0
        rbc    = 0.d0
        req    = 0.d0
        second = 0.d0
        third  = 0.d0
        vab    = 0.d0
        vac    = 0.d0
        vbc    = 0.d0
        zeroad = 0.d0
        zev    = 0.d0
        zkcal  = 0.d0
        c = reshape([ 77.45d0, -0.4071d0, -0.508d0, 0.1489d0, 1.013d0, -10.495d0, 9.1484d0, 10.273d0, 0.0d0, 0.0d0, - &
            & 3.050d0, -33.78d0, -22.951d0, 0.0d0, 0.0d0, -1.3699d0, 0.4101d0, 0.3411d0, 16.906d0, 0.0005637d0, - &
            & 0.6359d0, -0.06253d0, -0.1225d0, 0.0d0, 0.0d0, -0.00906d0, 0.4435d0, 0.7439d0, 0.0d0, 0.0d0, -0.3498d0, &
            & 0.6617d0, 0.8677d0, 1.3858d0, 233.36d0, -4.756d0, -1.626d0, 0.6337d0, 0.0d0, 0.0d0, -476.32d0, - &
            & 1.2104d0, -1.2762d0, 0.0d0, 0.0d0], shape(c))
        de     = [ 0.1559d0, 0.1779d0, 0.1559d0]
        gamma  = [ 1.2670d0, 1.4694d0, 1.2670d0]
        req    = [ 1.8460d0, 2.3158d0, 1.8460d0]
        zeroad = -0.022000295721d0
        alpha  = 0.9172d0
        !
        !          SET BETA, UNITS: RECIPROCAL BOHR
        beta = 1.4694d0
        !
        !   Initialize the potential values in the three asymptotic valleys.
        easyab = de(1)
        easybc = de(2)
        easyac = de(3)
        !
        !   Echo potential paramters to unit IPRT.
        !
        !   Convert ZEROAD to eV and kcal/mol and echo to unit IPRT.
        zev   = zeroad*27.21106d0
        zkcal = zeroad*627.5095d0
        !
        !CCCCCCCCCCC
        !
        !
        !CCCCCCCCCCC
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        MOVE ONTO THE MAIN PART OF THE ROUTINE
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !
        !   Check the values of NSURF and NDER
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !         SET RAB,RBC,AND RAC
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        rbc = r(2)
        !
        !           RAB IS DEFINED AS THE SMALLER OF R(1) AND R(3)
        !           RAC IS DEFINED AS THE LARGER OF R(1) AND R(3)
        !
        if( r(1) <= r(3) )then
            rab = r(1)
            rac = r(3)
        else
            rab = r(3)
            rac = r(1)
        endif
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE USEFUL EXPONENTIALS
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        expo(1) = exp(-alpha*rab)
        expo(2) = exp(-beta*rbc)
        expo(3) = exp(-alpha*rac)
        !
        !          AND THEIR DERIVATIVES
        !
        if( nder .eq. 1 )then
            dexpo(1) = -alpha*expo(1)
            dexpo(2) = -beta*expo(2)
            dexpo(3) = -alpha*expo(3)
        endif
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE VAB AND DV(1)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !          THE POTENTIAL CALCULATED HERE AND THE POTENTIALS VAC AND VBC
        !          CORRESPOND TO THE VOH POTENTIAL DESCRIBED IN EQUATION 3.
        !          DV(1) IS THE DERIVATIVE OF THIS TERM WITH RESPECT TO RAB.
        !
        junk(1) = exp(gamma(1)*(req(1)-rab))
        vab     = de(1)*(junk(1)*junk(1)-2.0d0*junk(1))
        if( nder .eq. 1 ) &
            & dv(1) = 2.0d0*de(1)*(1.0d0-junk(1))*gamma(1)*junk(1)
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE VAC AND DV(3)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        junk(1) = exp(gamma(3)*(req(3)-rac))
        vac     = de(3)*(junk(1)*junk(1)-2.0d0*junk(1))
        if( nder .eq. 1 ) &
            & dv(3) = 2.0d0*de(3)*(1.0d0-junk(1))*gamma(3)*junk(1)
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE VBC AND DV(2)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        junk(1) = exp(gamma(2)*(req(2)-rbc))
        vbc     = de(2)*(junk(1)*junk(1)-2.0d0*junk(1))
        if( nder .eq. 1 ) &
            & dv(2) = -2.0d0*de(2)*(junk(1)-1.0d0)*gamma(2)*junk(1)
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE CSTH AND DERIVATIVES
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        !          TH IS THE BEND ANGLE OF THE THREE ATOM SYSTEM.
        !          IT IS NEVER USED USED EXPLICITELY WITHIN THE PROGRAM
        !          HOWEVER COSINE(TH) IS.  COSINE(TH) IS CALCULATED
        !          FROM THE LAW OF COSINES)
        !
        num   = rbc*rbc+rab*rab-rac*rac
        denom = 2.0d0*rab*rbc
        !
        csth = num/denom
        !
        if( nder .eq. 1 )then
            dcsth(1) = (2.0d0*rab/denom)-(2.0d0*rbc)*csth/denom
            dcsth(2) = (2.0d0*rbc/denom)-(2.0d0*rab)*csth/denom
            dcsth(3) = -(2.0d0*rac/denom)
        endif
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE THE B IJ 'S--SEE EQUATION 5
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        do i = 1, 5
            do j = 1, 3
                b(i, j) = c(i,j,1)*(1.0d0+c(i,j,2)*csth*(1.0d0+c(i,j,3)* csth))
                !
                !               CALCULATE THE DERIVATIVES TOO.
                !
                if( nder .eq. 1 )then
                    do it = 1, 3
                        db(i, j, it) = c(i,j,1)*c(i,j,2)*(csth*c(i,j,3)* dcsth(it)+dcsth(it)* (1.0d0+c(i,j,3)*csth))
                    enddo
                endif
            enddo
        enddo
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          NEXT CALCULATE MEW--SEE EQUATION FOUR
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        mew = expo(2)-0.02538745816335d0
        if( nder .eq. 1 )then
            dmew(2) = -beta*expo(2)
            dmew(1) = 0.0d0
            dmew(3) = 0.0d0
        endif
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          NOW CALCULATE A I 'S AND DERIVATIVES--SEE EQUATION 4
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        do it = 1, 5
            junk(1) = 1.0d0+b(it, 3)*mew
            junk(2) = mew*b(it, 2)
            junk(3) = 1.0d0+junk(2)*junk(1)
            !
            a(it) = b(it, 1)*junk(3)
            if( nder .eq. 1 )then
                do it1 = 1, 3
                    da(it,it1) = db(it, 1, it1)*junk(3)+b(it, 1)* ((db(it, 2, it1)*mew+b(it, 2)*dmew(it1))* junk(1)+ &
                        & junk(2)* (db(it, 3, it1)*mew+ b(it, 3)*dmew(it1)))
                enddo
            endif
        enddo
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE THE SECOND TERM IN EQUATION THREE
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        prod1 = expo(1)*expo(3)*a(1)
        !
        prod2 = 1.0d0+(a(3)*rab*rac)/(rab+rac)
        !
        if( nder .eq. 1 )then
            dprod1(1) = dexpo(1)*expo(3)*a(1)+expo(1)*expo(3)*da(1,1)
            dprod1(2) = expo(1)*expo(3)*da(1,2)
            dprod1(3) = dexpo(3)*expo(1)*a(1)+expo(1)*expo(3)*da(1,3)
            !
            dprod2(2) = (da(3,2)*rab*rac)/(rab+rac)
            dprod2(1) = ((da(3,1)*rab*rac+a(3)*rac)*(rab+rac)- (a(3)*rab*rac))/((rab+rac)*(rab+rac))
            dprod2(3) = ((da(3,3)*rab*rac+a(3)*rab)*(rab+rac)- (a(3)*rab*rac))/((rab+rac)*(rab+rac))
        endif
        !
        !          CALCULATE PROD3
        !
        prod3 = 1.0d0+a(2)*(rab+rac)*prod2
        !
        if( nder .eq. 1 )then
            do it = 1, 3
                dprod3(it) = da(2,it)*(rab+rac)*prod2+a(2)* (rab+rac)*dprod2(it)
                if( it .eq. 2 ) cycle
                dprod3(it) = dprod3(it)+a(2)*prod2
            enddo
        endif
        !
        second = prod1*prod3
        !
        if( nder .eq. 1 )then
            do it = 1, 3
                dsec(it) = prod1*dprod3(it)+dprod1(it)*prod3
            enddo
        endif
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE THE THIRD TERM
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        junk(1) = expo(1)+expo(3)
        junk(2) = (1.0d0+a(5)*(rab+rac))
        !
        third = junk(1)*expo(2)*a(4)*junk(2)
        !
        if( nder .eq. 1 )then
            do it1 = 1, 2
                it = 1
                if( it1 .eq. 2 ) it = 3
                dthird(it) = dexpo(it)*expo(2)*a(4)*junk(2)+ junk(1)*expo(2)*da(4,it)*junk(2)+ junk(1)*expo(2)* a(4)* &
                    & (da(5,it)* (rab+rac)+a(5))
            enddo
            dthird(2) = junk(1)*(dexpo(2)*a(4)*junk(2)+expo(2)*(da(4,2)* junk(2)+a(4)*(da(5,2)*(rab+rac))))
        endif
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          FINALLY SUM IT ALL UP
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        energy = vab+vac+vbc+second+third+de(2)
        energy = energy+zeroad
        !
        if( nder .eq. 1 )then
            do it = 1, 3
                der(it) = dv(it)+dsec(it)+dthird(it)
            enddo
            !
            goto_block_160_476: block
                if( r(1) .lt. r(3) ) exit goto_block_160_476
                junk(1) = der(3)
                der(3)  = der(1)
                der(1)  = junk(1)
            end block goto_block_160_476
        endif
        !
        !CCCCCCCCCCCCCCC
        !
        !
        !CCCCCCCCCCCCCCC
        !
        !
    end subroutine evaluate_internal
end module chempotpy_o2h_gen_pot_legacy_state_fffab58a

!> @brief Evaluate O2H for one set of internal distances.
subroutine pot( r, nder, energy, der )
    use chempotpy_o2h_gen_pot_legacy_state_fffab58a, only: evaluate_internal
    implicit none
    double precision, intent(in)  :: r(3)
    integer, intent(in)           :: nder
    double precision, intent(out) :: energy, der(3)
    call evaluate_internal(r, nder, energy, der)
end subroutine pot

!
!
