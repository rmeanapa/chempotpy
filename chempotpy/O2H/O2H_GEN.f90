!   System: HO2
!   Reference: C. F. Melius and R. J. Blint
!              Chem. Phys. Lett. 64, 183 (1979).
!   The potential parameters are included in DATA statements.
!   and in the block data subprogram PPARMS.
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
!   Default Parameter Values:
!      Variable      Default value
!      NSURF            0
!      NDER             1
!      IPTPRT           6
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         GENERAL INFORMATION
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          THE POTENTIAL CONSISTS OF THREE MORSE TERMS FOR
!          THE THREE DIATOMICS, HO,OO, AND HO RESPECTIVELY
!          PLUS TWO INTERACTION TERMS.
!          EACH POTENTIAL TERM WILL BE CALCULATED INDIVIDUALLY
!          TO SIMPLIFY READING THE CODE.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          DEFINE VARIABLES, ARRAYS, AND COMMON BLOCK
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!> @brief Implement O2H_GEN, a potential energy surface for O2H.
module chempotpy_o2h_gen_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the O2H_GEN potential energy surface for O2H.
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
        real(wp)              :: der(3), energy, r(3), v
        integer               :: iatom, idir, istate, j
        real(wp)              :: tx(9), r2(3), dedx(9), drdx(3,9)
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
end module chempotpy_o2h_gen_pes

!> @brief Expose O2H_GEN through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_o2h_gen_pes, only: evaluate
    use chempotpy_kinds,       only: wp
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

!> @brief Evaluate EvdRdX for the O2H_GEN potential.
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

module chempotpy_o2h_gen_pot_legacy_state_fffab58a
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_internal

    contains

    !> @brief Evaluate dispatch.
    subroutine evaluate_internal( r, nder, energy, der )
        use chempotpy_kinds, only: wp
        implicit none
        intent(in)  :: r, nder
        intent(out) :: energy, der
        real(wp) :: a(5), alpha, b(5,3), beta, c(5,3,3), csth, da(5,3), db(5,3,3), dcsth(3), de(3), der(3), dmew(3), &
            & dprod1(3), dprod2(3), dprod3(3), dsec(3), dthird(3), dv(3), easyab, easyac, easybc, energy, prod1, prod2, prod3, &
            & r(3), rab, rac, rbc, req(3), second, third, vab, vac, vbc, zeroad, zev, zkcal
        integer  :: i, idum, iprt, it, it1, j, nder, ndum, nsurf
        real(wp) :: num, denom, gamma(3), junk(4), expo(3), dexpo(3), mew
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          Initialize the constants
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !   The array DE contains the dissociation energies in hartree atomic un
        !   The array GAMMA contains the Morse Betas in reciprocal bohrs.
        !   The array REQ contains the equilibrium bond lengths in bohr.
        !   The variable ZEROAD contains the constant added to the energy in har
        !   atomic units.
        !          SET ALPHA, UNITS: RECIPROCAL BOHR
        a      = 0.e0_wp
        alpha  = 0.e0_wp
        b      = 0.e0_wp
        beta   = 0.e0_wp
        c      = 0.e0_wp
        csth   = 0.e0_wp
        da     = 0.e0_wp
        db     = 0.e0_wp
        dcsth  = 0.e0_wp
        de     = 0.e0_wp
        denom  = 0.e0_wp
        dexpo  = 0.e0_wp
        dmew   = 0.e0_wp
        dprod1 = 0.e0_wp
        dprod2 = 0.e0_wp
        dprod3 = 0.e0_wp
        dsec   = 0.e0_wp
        dthird = 0.e0_wp
        dv     = 0.e0_wp
        easyab = 0.e0_wp
        easyac = 0.e0_wp
        easybc = 0.e0_wp
        expo   = 0.e0_wp
        gamma  = 0.e0_wp
        junk   = 0.e0_wp
        mew    = 0.e0_wp
        num    = 0.e0_wp
        prod1  = 0.e0_wp
        prod2  = 0.e0_wp
        prod3  = 0.e0_wp
        rab    = 0.e0_wp
        rac    = 0.e0_wp
        rbc    = 0.e0_wp
        req    = 0.e0_wp
        second = 0.e0_wp
        third  = 0.e0_wp
        vab    = 0.e0_wp
        vac    = 0.e0_wp
        vbc    = 0.e0_wp
        zeroad = 0.e0_wp
        zev    = 0.e0_wp
        zkcal  = 0.e0_wp
        c = reshape([ 77.45e0_wp, -0.4071e0_wp, -0.508e0_wp, 0.1489e0_wp, 1.013e0_wp, -10.495e0_wp, 9.1484e0_wp, &
            & 10.273e0_wp, 0.0e0_wp, 0.0e0_wp, - 3.050e0_wp, -33.78e0_wp, -22.951e0_wp, 0.0e0_wp, 0.0e0_wp, - &
            & 1.3699e0_wp, 0.4101e0_wp, 0.3411e0_wp, 16.906e0_wp, 0.0005637e0_wp, - 0.6359e0_wp, -0.06253e0_wp, - &
            & 0.1225e0_wp, 0.0e0_wp, 0.0e0_wp, -0.00906e0_wp, 0.4435e0_wp, 0.7439e0_wp, 0.0e0_wp, 0.0e0_wp, - &
            & 0.3498e0_wp, 0.6617e0_wp, 0.8677e0_wp, 1.3858e0_wp, 233.36e0_wp, -4.756e0_wp, -1.626e0_wp, 0.6337e0_wp, &
            & 0.0e0_wp, 0.0e0_wp, -476.32e0_wp, - 1.2104e0_wp, -1.2762e0_wp, 0.0e0_wp, 0.0e0_wp], shape(c))
        de     = [ 0.1559e0_wp, 0.1779e0_wp, 0.1559e0_wp]
        gamma  = [ 1.2670e0_wp, 1.4694e0_wp, 1.2670e0_wp]
        req    = [ 1.8460e0_wp, 2.3158e0_wp, 1.8460e0_wp]
        zeroad = -0.022000295721e0_wp
        alpha  = 0.9172e0_wp
        !          SET BETA, UNITS: RECIPROCAL BOHR
        beta = 1.4694e0_wp
        !   Initialize the potential values in the three asymptotic valleys.
        easyab = de(1)
        easybc = de(2)
        easyac = de(3)
        !   Echo potential paramters to unit IPRT.
        !   Convert ZEROAD to eV and kcal/mol and echo to unit IPRT.
        zev   = zeroad*27.21106e0_wp
        zkcal = zeroad*627.5095e0_wp
        !CCCCCCCCCCC
        !CCCCCCCCCCC
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !        MOVE ONTO THE MAIN PART OF THE ROUTINE
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !   Check the values of NSURF and NDER
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !         SET RAB,RBC,AND RAC
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCC
        rbc = r(2)
        !           RAB IS DEFINED AS THE SMALLER OF R(1) AND R(3)
        !           RAC IS DEFINED AS THE LARGER OF R(1) AND R(3)
        if( r(1) <= r(3) )then
            rab = r(1)
            rac = r(3)
        else
            rab = r(3)
            rac = r(1)
        endif
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE USEFUL EXPONENTIALS
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        expo(1) = exp(-alpha*rab)
        expo(2) = exp(-beta*rbc)
        expo(3) = exp(-alpha*rac)
        !          AND THEIR DERIVATIVES
        if( nder .eq. 1 )then
            dexpo(1) = -alpha*expo(1)
            dexpo(2) = -beta*expo(2)
            dexpo(3) = -alpha*expo(3)
        endif
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE VAB AND DV(1)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          THE POTENTIAL CALCULATED HERE AND THE POTENTIALS VAC AND VBC
        !          CORRESPOND TO THE VOH POTENTIAL DESCRIBED IN EQUATION 3.
        !          DV(1) IS THE DERIVATIVE OF THIS TERM WITH RESPECT TO RAB.
        junk(1) = exp(gamma(1)*(req(1)-rab))
        vab     = de(1)*(junk(1)*junk(1)-2.0e0_wp*junk(1))
        if( nder .eq. 1 ) &
            & dv(1) = 2.0e0_wp*de(1)*(1.0e0_wp-junk(1))*gamma(1)*junk(1)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE VAC AND DV(3)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        junk(1) = exp(gamma(3)*(req(3)-rac))
        vac     = de(3)*(junk(1)*junk(1)-2.0e0_wp*junk(1))
        if( nder .eq. 1 ) &
            & dv(3) = 2.0e0_wp*de(3)*(1.0e0_wp-junk(1))*gamma(3)*junk(1)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE VBC AND DV(2)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        junk(1) = exp(gamma(2)*(req(2)-rbc))
        vbc     = de(2)*(junk(1)*junk(1)-2.0e0_wp*junk(1))
        if( nder .eq. 1 ) &
            & dv(2) = -2.0e0_wp*de(2)*(junk(1)-1.0e0_wp)*gamma(2)*junk(1)
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE CSTH AND DERIVATIVES
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          TH IS THE BEND ANGLE OF THE THREE ATOM SYSTEM.
        !          IT IS NEVER USED USED EXPLICITELY WITHIN THE PROGRAM
        !          HOWEVER COSINE(TH) IS.  COSINE(TH) IS CALCULATED
        !          FROM THE LAW OF COSINES)
        num   = rbc*rbc+rab*rab-rac*rac
        denom = 2.0e0_wp*rab*rbc
        csth  = num/denom
        if( nder .eq. 1 )then
            dcsth(1) = (2.0e0_wp*rab/denom)-(2.0e0_wp*rbc)*csth/denom
            dcsth(2) = (2.0e0_wp*rbc/denom)-(2.0e0_wp*rab)*csth/denom
            dcsth(3) = -(2.0e0_wp*rac/denom)
        endif
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE THE B IJ 'S--SEE EQUATION 5
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        do i = 1, 5
            do j = 1, 3
                b(i,j) = c(i,j,1)*(1.0e0_wp+c(i,j,2)*csth*(1.0e0_wp+c(i,j,3)* csth))
                !               CALCULATE THE DERIVATIVES TOO.
                if( nder .eq. 1 )then
                    do it = 1, 3
                        db(i,j,it) = c(i,j,1)*c(i,j,2)*(csth*c(i,j,3)* dcsth(it)+dcsth(it)* (1.0e0_wp+c(i,j,3)*csth))
                    enddo
                endif
            enddo
        enddo
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          NEXT CALCULATE MEW--SEE EQUATION FOUR
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        mew = expo(2)-0.02538745816335e0_wp
        if( nder .eq. 1 )then
            dmew(2) = -beta*expo(2)
            dmew(1) = 0.0e0_wp
            dmew(3) = 0.0e0_wp
        endif
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          NOW CALCULATE A I 'S AND DERIVATIVES--SEE EQUATION 4
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        do it = 1, 5
            junk(1) = 1.0e0_wp+b(it,3)*mew
            junk(2) = mew*b(it,2)
            junk(3) = 1.0e0_wp+junk(2)*junk(1)
            a(it)   = b(it,1)*junk(3)
            if( nder .eq. 1 )then
                do it1 = 1, 3
                    da(it,it1) = db(it,1,it1)*junk(3)+b(it,1)* ((db(it,2,it1)*mew+b(it,2)*dmew(it1))* junk(1)+ &
                        & junk(2)* (db(it,3,it1)*mew+ b(it,3)*dmew(it1)))
                enddo
            endif
        enddo
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE THE SECOND TERM IN EQUATION THREE
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        prod1 = expo(1)*expo(3)*a(1)
        prod2 = 1.0e0_wp+(a(3)*rab*rac)/(rab+rac)
        if( nder .eq. 1 )then
            dprod1(1) = dexpo(1)*expo(3)*a(1)+expo(1)*expo(3)*da(1,1)
            dprod1(2) = expo(1)*expo(3)*da(1,2)
            dprod1(3) = dexpo(3)*expo(1)*a(1)+expo(1)*expo(3)*da(1,3)
            dprod2(2) = (da(3,2)*rab*rac)/(rab+rac)
            dprod2(1) = ((da(3,1)*rab*rac+a(3)*rac)*(rab+rac)- (a(3)*rab*rac))/((rab+rac)*(rab+rac))
            dprod2(3) = ((da(3,3)*rab*rac+a(3)*rab)*(rab+rac)- (a(3)*rab*rac))/((rab+rac)*(rab+rac))
        endif
        !          CALCULATE PROD3
        prod3 = 1.0e0_wp+a(2)*(rab+rac)*prod2
        if( nder .eq. 1 )then
            do it = 1, 3
                dprod3(it) = da(2,it)*(rab+rac)*prod2+a(2)* (rab+rac)*dprod2(it)
                if( it .eq. 2 ) cycle
                dprod3(it) = dprod3(it)+a(2)*prod2
            enddo
        endif
        second = prod1*prod3
        if( nder .eq. 1 )then
            do it = 1, 3
                dsec(it) = prod1*dprod3(it)+dprod1(it)*prod3
            enddo
        endif
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          CALCULATE THE THIRD TERM
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        junk(1) = expo(1)+expo(3)
        junk(2) = (1.0e0_wp+a(5)*(rab+rac))
        third   = junk(1)*expo(2)*a(4)*junk(2)
        if( nder .eq. 1 )then
            do it1 = 1, 2
                it = 1
                if( it1 .eq. 2 ) it = 3
                dthird(it) = dexpo(it)*expo(2)*a(4)*junk(2)+ junk(1)*expo(2)*da(4,it)*junk(2)+ junk(1)*expo(2)* a(4)* &
                    & (da(5,it)* (rab+rac)+a(5))
            enddo
            dthird(2) = junk(1)*(dexpo(2)*a(4)*junk(2)+expo(2)*(da(4,2)* junk(2)+a(4)*(da(5,2)*(rab+rac))))
        endif
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !          FINALLY SUM IT ALL UP
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        energy = vab+vac+vbc+second+third+de(2)
        energy = energy+zeroad
        if( nder .eq. 1 )then
            do it = 1, 3
                der(it) = dv(it)+dsec(it)+dthird(it)
            enddo
            goto_block_160_476: block
                if( r(1) .lt. r(3) ) exit goto_block_160_476
                junk(1) = der(3)
                der(3)  = der(1)
                der(1)  = junk(1)
            end block goto_block_160_476
        endif
        !CCCCCCCCCCCCCCC
        !CCCCCCCCCCCCCCC
    end subroutine evaluate_internal
end module chempotpy_o2h_gen_pot_legacy_state_fffab58a

!> @brief Evaluate O2H for one set of internal distances.
subroutine pot( r, nder, energy, der )
    use chempotpy_o2h_gen_pot_legacy_state_fffab58a, only: evaluate_internal
    use chempotpy_kinds,                             only: wp
    implicit none
    real(wp), intent(in)  :: r(3)
    integer, intent(in)   :: nder
    real(wp), intent(out) :: energy, der(3)
    call evaluate_internal(r, nder, energy, der)
end subroutine pot
