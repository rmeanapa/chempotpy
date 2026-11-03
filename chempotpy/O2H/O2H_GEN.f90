!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the O2H_GEN potential.
subroutine pes( x, igrad, p, g, d )
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
    double precision              :: der, energy, r, v
    integer                       :: i, iatom, idir, istate, j, nder, ndum, nsurf
    common /pt1cm/ r(3), energy, der(3)
    common /pt2cm/ nsurf, nder, ndum(21)
    double precision :: tx(9), r2(3), dedx(9), drdx(3,9)
    logical, save    :: first_time_data=.true.
    !initialize
    v = 0.d0; g = 0.d0; d = 0.d0

    do iatom=1, natoms
        do idir=1, 3
            j     = (iatom-1)*3+idir
            tx(j) = x(iatom,idir)
        end do
    end do
    ! input cartesian is ClHH
    call legacy_triatomic_distances_to_bohr( x, r )

    nsurf = 0
    nder  = igrad
    der   = 0.d0

    if(first_time_data) then
        call prepot
        first_time_data = .false.
    end if

    call pot

    energy = energy*legacy_hartree_to_electronvolt
    der    = der*51.422067

    r2 = r*legacy_bohr_to_angstrom
    call evdrdx(tx, r2, drdx)
    dedx = 0.d0
    do i=1, 9
        do j=1, 3
            dedx(i) = dedx(i)+der(j)*drdx(j,i)
        end do
    end do

    if (igrad==0) then
        do istate=1, nstates
            p(istate) = energy
        end do
    else if (igrad==1) then
        do istate=1, nstates
            p(istate) = energy
        end do
        do iatom=1, natoms
            do idir=1, 3
                j               = (iatom-1)*3+idir
                g(1,iatom,idir) = dedx(j)
            end do
        end do
    else if (igrad==2) then
        write (*, *) 'Only energy and gradient are available'
    end if

end subroutine pes

!> @brief Evaluate EvdRdX for the O2H_GEN potential.
subroutine evdrdx( x, r, drdx )
    implicit none
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    integer                       :: i, j
    ! Initialize dRdX(3,9)
    do i=1, 3
        do j=1, 9
            drdx(i,j) = 0.0d0
        end do
    end do

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
!> @brief Evaluate PREPOT for the O2H_GEN potential.
subroutine prepot
    !
    !   System: HO2
    !   Reference: C. F. Melius and R. J. Blint
    !              Chem. Phys. Lett. 64, 183 (1979).
    !
    !   PREPOT must be called once before any calls to POT.
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
    double precision :: a, alpha, b, beta, c, csth, da, db, dcsth, de, der, dmew, dprod1, dprod2, &
        & dprod3, dsec, dthird, dv, easyab, easyac, easybc, energy, prod1, prod2, prod3, r, rab, rac, &
        & rbc, req, second, third, vab, vac, vbc, zeroad, zev, zkcal
    integer :: i, idum, iprt, it, it1, j, nder, ndum, nsurf
    !
    double precision :: num, denom, gamma, junk, expo, dexpo, mew
    !
    !
    common /pt1cm/ r(3), energy, der(3)
    common /pt2cm/ nsurf, nder, ndum(20)
    common /pt4cm/ iprt, idum(19)
    common /pt5cm/ easyab, easybc, easyac
    !
    dimension dcsth(3), dv(3), gamma(3), junk(4), c(5,3,3), a(5), da(5,3), &
        & req(3), de(3), expo(3), dexpo(3), dmew(3), b(5, 3), db(5, 3, 3), dprod1(3 &
        & ), dprod2(3), dsec(3), dthird(3), dprod3(3)
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !          Initialize the constants
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    data c / 77.45d0, -0.4071d0, -0.508d0, &
        & 0.1489d0, 1.013d0, -10.495d0, &
        & 9.1484d0, 10.273d0, 0.0d0, 0.0d0, &
        & -3.050d0, -33.78d0, -22.951d0, 0.0d0, &
        & 0.0d0, -1.3699d0, 0.4101d0, 0.3411d0, &
        & 16.906d0, 0.0005637d0, -0.6359d0, &
        & -0.06253d0, -0.1225d0, 0.0d0, 0.0d0, &
        & -0.00906d0, 0.4435d0, 0.7439d0, 0.0d0, 0.0d0, &
        & -0.3498d0, 0.6617d0, 0.8677d0, 1.3858d0, &
        & 233.36d0, -4.756d0, -1.626d0, &
        & 0.6337d0, 0.0d0, 0.0d0, -476.32d0, &
        & -1.2104d0, -1.2762d0, 0.0d0, 0.0d0 /
    !   The array DE contains the dissociation energies in hartree atomic un
    data de / 0.1559d0, 0.1779d0, 0.1559d0/
    !   The array GAMMA contains the Morse Betas in reciprocal bohrs.
    data gamma / 1.2670d0, 1.4694d0, 1.2670d0/
    !   The array REQ contains the equilibrium bond lengths in bohr.
    data req / 1.8460d0, 2.3158d0, 1.8460d0/
    !   The variable ZEROAD contains the constant added to the energy in har
    !   atomic units.
    data zeroad /-0.022000295721d0/
    save
    !
    !          SET ALPHA, UNITS: RECIPROCAL BOHR
    alpha = 0.9172d0
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
    return
    !
    !CCCCCCCCCCC
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !        MOVE ONTO THE MAIN PART OF THE ROUTINE
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    entry pot
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
    if (r(1).gt.r(3)) go to 40
    !
    rab = r(1)
    rac = r(3)
    go to 50
    !
    40 rab = r(3)
    rac = r(1)
    50 continue
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
    if (nder .eq. 1) then
        dexpo(1) = -alpha*expo(1)
        dexpo(2) = -beta*expo(2)
        dexpo(3) = -alpha*expo(3)
    end if
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
    if (nder .eq. 1) &
        & dv(1) = 2.0d0*de(1)*(1.0d0-junk(1))*gamma(1)*junk(1)
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !          CALCULATE VAC AND DV(3)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    junk(1) = exp(gamma(3)*(req(3)-rac))
    vac     = de(3)*(junk(1)*junk(1)-2.0d0*junk(1))
    if (nder .eq. 1) &
        & dv(3) = 2.0d0*de(3)*(1.0d0-junk(1))*gamma(3)*junk(1)
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !          CALCULATE VBC AND DV(2)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    junk(1) = exp(gamma(2)*(req(2)-rbc))
    vbc     = de(2)*(junk(1)*junk(1)-2.0d0*junk(1))
    if (nder .eq. 1) &
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
    if (nder .eq. 1) then
        dcsth(1) = (2.0d0*rab/denom)-(2.0d0*rbc)*csth/denom
        dcsth(2) = (2.0d0*rbc/denom)-(2.0d0*rab)*csth/denom
        dcsth(3) = -(2.0d0*rac/denom)
    end if
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !          CALCULATE THE B IJ 'S--SEE EQUATION 5
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    do 80 i = 1, 5
        do 70 j = 1, 3
            b(i, j) = c(i,j,1)*(1.0d0+c(i,j,2)*csth*(1.0d0+c(i,j,3)* csth))
            !
            !               CALCULATE THE DERIVATIVES TOO.
            !
            if (nder .eq. 1) then
                do 60 it = 1, 3
                    db(i, j, it) = c(i,j,1)*c(i,j,2)*(csth*c(i,j,3)* dcsth(it)+dcsth(it)* (1.0d0+c(i,j,3)*csth))
                60 continue
            end if
        70 continue
    80 continue
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !          NEXT CALCULATE MEW--SEE EQUATION FOUR
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    mew = expo(2)-0.02538745816335d0
    if (nder .eq. 1) then
        dmew(2) = -beta*expo(2)
        dmew(1) = 0.0d0
        dmew(3) = 0.0d0
    end if
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !          NOW CALCULATE A I 'S AND DERIVATIVES--SEE EQUATION 4
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    do 100 it = 1, 5
        junk(1) = 1.0d0+b(it, 3)*mew
        junk(2) = mew*b(it, 2)
        junk(3) = 1.0d0+junk(2)*junk(1)
        !
        a(it) = b(it, 1)*junk(3)
        if (nder .eq. 1) then
            do 90 it1 = 1, 3
                da(it,it1) = db(it, 1, it1)*junk(3)+b(it, 1)* ((db(it, 2, it1)*mew+b(it, 2)*dmew(it1))* junk(1)+ &
                    & junk(2)* (db(it, 3, it1)*mew+ b(it, 3)*dmew(it1)))
            90 continue
        end if
    100 continue
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !          CALCULATE THE SECOND TERM IN EQUATION THREE
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    prod1 = expo(1)*expo(3)*a(1)
    !
    prod2 = 1.0d0+(a(3)*rab*rac)/(rab+rac)
    !
    if (nder .eq. 1) then
        dprod1(1) = dexpo(1)*expo(3)*a(1)+expo(1)*expo(3)*da(1,1)
        dprod1(2) = expo(1)*expo(3)*da(1,2)
        dprod1(3) = dexpo(3)*expo(1)*a(1)+expo(1)*expo(3)*da(1,3)
        !
        dprod2(2) = (da(3,2)*rab*rac)/(rab+rac)
        dprod2(1) = ((da(3,1)*rab*rac+a(3)*rac)*(rab+rac)- (a(3)*rab*rac))/((rab+rac)*(rab+rac))
        dprod2(3) = ((da(3,3)*rab*rac+a(3)*rab)*(rab+rac)- (a(3)*rab*rac))/((rab+rac)*(rab+rac))
    end if
    !
    !          CALCULATE PROD3
    !
    prod3 = 1.0d0+a(2)*(rab+rac)*prod2
    !
    if (nder .eq. 1) then
        do 120 it = 1, 3
            dprod3(it) = da(2,it)*(rab+rac)*prod2+a(2)* (rab+rac)*dprod2(it)
            if (it.eq.2) go to 110
            dprod3(it) = dprod3(it)+a(2)*prod2
            110 continue
        120 continue
    end if
    !
    second = prod1*prod3
    !
    if (nder .eq. 1) then
        do 130 it = 1, 3
            dsec(it) = prod1*dprod3(it)+dprod1(it)*prod3
        130 continue
    end if
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
    if (nder .eq. 1) then
        do 140 it1 = 1, 2
            it = 1
            if (it1.eq.2) it = 3
            dthird(it) = dexpo(it)*expo(2)*a(4)*junk(2)+ junk(1)*expo(2)*da(4,it)*junk(2)+ junk(1)*expo(2)*a(4)* &
                & (da(5,it)* (rab+rac)+a(5))
        140 continue
        dthird(2) = junk(1)*(dexpo(2)*a(4)*junk(2)+expo(2)*(da(4,2)* junk(2)+a(4)*(da(5,2)*(rab+rac))))
    end if
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !          FINALLY SUM IT ALL UP
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    energy = vab+vac+vbc+second+third+de(2)
    energy = energy+zeroad
    !
    if (nder .eq. 1) then
        do 150 it = 1, 3
            der(it) = dv(it)+dsec(it)+dthird(it)
        150 continue
        !
        if (r(1).lt.r(3)) go to 160
        junk(1) = der(3)
        der(3)  = der(1)
        der(1)  = junk(1)
        160 continue
    end if
    !
    !CCCCCCCCCCCCCCC
    !
    !
    !CCCCCCCCCCCCCCC
    !
    !
end subroutine prepot

!
!*****
!
block data ptparm
    implicit none
    integer :: idum, iprt, nder, ndum, nsurf
    common /pt2cm/ nsurf, nder, ndum(21)
    common /pt4cm/ iprt, idum(19)
    !
    !   Initialize the flags and the I/O unit numbers for the potential
    data iprt /6/
    data nsurf /0/
    data nder /1/
    data ndum /21*0/
    !
end

!*****
