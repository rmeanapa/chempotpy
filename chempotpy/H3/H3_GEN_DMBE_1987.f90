!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the H3_GEN_DMBE_1987 potential.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_geometry,     only: cartesian_distance
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
    implicit none
    ! number of electronic state
    integer, parameter            :: nstates = 2
    integer, parameter            :: natoms  = 3
    integer, parameter            :: n3atom  = 75
    integer, parameter            :: natom   = 25
    integer, parameter            :: isurf   = 5
    integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    double precision              :: anuzero, cart, deesdr, degsdr, deijdr, descart, dgscart, dijcart, engyes, &
        & engygs, engyij, pengyes, pengygs, pengyij, r, v
    integer :: iatom, nasurf, nder, nflag, nulbl
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    logical, save    :: first_time_data=.true.
    double precision :: tx(9), r1(3), drdx(3,9)
    double precision :: hr(nstates,nstates,3), gr(nstates,3)
    double precision :: hx(nstates,nstates,9), gx(nstates,9)
    integer          :: istate, idir, i, j, k, l
    !initialize
    v = 0.d0; g = 0.d0; d = 0.d0

    cart = 0.d0
    do iatom=1, natoms
        do idir=1, 3
            cart(iatom, idir) = x(iatom,idir)/legacy_bohr_to_angstrom
        end do
    end do

    if (igrad.eq.0) then
        nder = 0
    else if (igrad.eq.1 .or. igrad.eq.2) then
        nder = 1
    end if

    nasurf = 1
    if(first_time_data) then
        call prepot
        first_time_data = .false.
    end if

    call pot

    p(1) = engygs*legacy_hartree_to_electronvolt
    do istate=2, nstates
        i         = istate-1
        p(istate) = engyes(i)*legacy_hartree_to_electronvolt
    end do

    ! now compute gradients
    do iatom=1, natoms
        do idir=1, 3
            j     = (iatom-1)*3+idir
            tx(j) = x(iatom,idir)
        end do
    end do
    ! input cartesian is HHH
    r1(1) = cartesian_distance(x(1,:), x(2,:))
    r1(2) = cartesian_distance(x(2,:), x(3,:))
    r1(3) = cartesian_distance(x(1,:), x(3,:))

    call evdrdx(tx, r1, drdx)
    hr = 0.d0
    gr = 0.d0
    do i=1, 3
        hr(1,2,i) = pengyij(i)/legacy_bohr_to_angstrom
        hr(2,1,i) = -pengyij(i)/legacy_bohr_to_angstrom
    end do
    do i=1, 3
        gr(1,i) = degsdr(i)*51.422067
        gr(2,i) = deesdr(i, 1)*51.422067
    end do

    hx = 0.d0
    gx = 0.d0
    do i=1, 9
        do j=1, 3
            hx(:,:,i) = hx(:,:,i)+hr(:,:,j)*drdx(j,i)
            gx(:,i)   = gx(:,i)+gr(:,j)*drdx(j,i)
        end do
    end do

    do iatom=1, natoms
        do idir=1, 3
            j                 = (iatom-1)*3+idir
            g(:,iatom,idir)   = gx(:,j)
            d(:,:,iatom,idir) = hx(:,:,j)
        end do
    end do

end subroutine pes

!> @brief Evaluate EvdRdX for the H3_GEN_DMBE_1987 potential.
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
!> @brief Evaluate PREPOT for the H3_GEN_DMBE_1987 potential.
subroutine prepot
    !
    !   System:          H3
    !   Functional Form: Double many-body expansion
    !   Common Name:     DMBE
    !   Reference:       A. J. C. Varandas, F. B. Brown, C. A. Mead,
    !                    D. G. Truhlar, and N. C. Blais
    !                    J. Chem. Phys. 86, 6258 (1987).
    !   Interface: potlib2001
    ! Number of bodies: 3
    ! Number of electronic states: 1
    ! Number of derivatives: 1

    !
    !   Calling Sequence:
    !      PREPOT - initializes the potential's variables and
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
    !   Surfaces:
    !      ground electronic state
    !
    !   Zero of energy:
    !      The classical potential energy is set equal to zero for the H
    !      infinitely far from the H2 diatomic and R(H2) set equal to the
    !      H2 equilibrium diatomic value.
    !
    !   Parameters:
    !      Set in the BLOCK DATA subprogram PTPACM
    !
    !   Coordinates:
    !      Internal, Definition: R(1) = R( first H -- second H)
    !                            R(2) = R(second H -- third H )
    !                            R(3) = R( third H -- first H )
    !
    !   Common Blocks (used between the calling program and this potential):
    !        passes the coordinates, ground state electronic energy, and
    !        derivatives of the ground electronic state energy with respect
    !        to the coordinates.
    !        passes the control flags where
    !              = 1, first excited electronic state energy and derivative
    !              = 2, ground and first excited electronic state energies a
    !                   derivatives
    !        NDER  - not used
    !        NFLAG - Control flags
    !
    !              derivatives for the excited electronic states
    !
    !
    !              derivatives for the nonadiabatic coupling between ground
    !              excited states.
    !
    !   Default Parameter Values:
    !      Variable      Default value
    !      NFLAG(18)        6
    !
    !***********************************************************************
    !     Calculates double many-body expansion of the H3 potential.
    !        ENGYGS  - energy of surface 1 (Ground State).        [ARRAY]
    !        DEGSDR  - derivatives for surface 1 (Ground State).  [ARRAY]
    !        ENGYES  - energy of surface 2 (Excited State).       [ARRAY]
    !        DEESDR  - derivatives for surface 2 (Excited State). [ARRAY]
    !        ENGYIJ  - nonadiabatic coupling energy.              [ARRAY]
    !        CAPF   - nonadiabatic coupling.  [SUBROUTINE TO CALCULATE ENGYI
    !***********************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, az2, bet0, bet1, beta1, beta2, beta3, &
        & cart, cartnu, cd0, cd1, chh, ck0, ck1, cnvrtchi, cnvrtd, cnvrte, cnvrtq, cnvrts, corrs, corrt, cphi3, damp, &
        & dampa, dampb, dcorrs, dcorrt, dcphi1, dcphi2, dcphi3, ddamp, dq1, dq2, dq3, drho1, drho2, drho3, ds1, ds2, ds3, &
        & ezero, h2r0, h2rm, h2rmt, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, per, per2, qcoord, &
        & r12, r22, r32, rho, rhol, s, s2, sqrt3, xpar
    integer :: i, icartr, icouple, id, indexes, irctnt, mder, msurf, nasurf, natoms, nder, nexp, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt3cm/  ezero(isurf+1)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    common /h3crcm/ corrs(3), dcorrs(3), &
        & corrt(3), dcorrt(3), &
        & damp(3, 3), ddamp(3, 3)
    !
    !**************************************************
    !                                                 *
    !  units conversion for anomalous adiabatic       *
    !  coupling pengyij                               *
    !                                                 *
    !**************************************************
    !
    common /ijconvert/cnvrtq, cnvrts, cnvrtchi, icouple
    !
    !     DIMENSION DLEP(3), DLEP2(3), DV2(3), DV3(3), DVA(3), DC3(3)
    !     DIMENSION R(3)
    !
    if(natoms.gt.25) then
        write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    end if
    !
    ! NON STANDARD COUPLING BETWEEN SURFACES
    !
    icouple      = 1
    nasurf(1, 2) = 0
    nasurf(2, 1) = 0
    cnvrtd       = 1.d0
    if(nflag(1).eq.2) cnvrtd = .529177249d0
    cnvrte = 1.d0
    if(nflag(2).eq.2) cnvrte = 1000.d0
    if(nflag(2).eq.3) cnvrte = 27.211395d0
    if(nflag(2).eq.4) cnvrte = 627.5095d0
    if(nflag(2).eq.5) cnvrte = 219474.6d0
    if(nflag(2).eq.6) cnvrte = 2625.500d0
    cnvrtq   = cnvrte/cnvrtd
    cnvrts   = cnvrte/cnvrtd
    cnvrtchi = cnvrte
    !
    nexp = 4
    do 10 id = 1, 3
        nexp      = nexp + 2
        dampa(id) = alph0/dble(nexp)**alph1
        dampb(id) = bet0*exp(-bet1*dble(nexp))
    10 continue
    !
    !
    do i=1, 5
        ref(i) = ' '
    end do
    !
    ref(1) = 'A. J. C. Varandas, F. B. Brown, C. A. Mead,'
    ref(2) = 'D. G. Truhlar, N. C. Blais,'
    ref(3) = 'J. Chem. Phys. 86, 6258(1987)'
    !
    indexes(1) = 1
    indexes(2) = 1
    indexes(3) = 1
    !
    irctnt = 2
    !
    !      CALL POTINFO
    !
    call ancvrt
    !
end subroutine prepot

!
!> @brief Evaluate POT for the H3_GEN_DMBE_1987 potential.
subroutine pot
    !
    !   Surfaces:
    !      ground electronic state
    !
    !   Zero of energy:
    !      The classical potential energy is set equal to zero for the H
    !      infinitely far from the H2 diatomic and R(H2) set equal to the
    !      H2 equilibrium diatomic value.
    !
    !   Coordinates:
    !      Internal, Definition: R(1) = R( first H -- second H)
    !                            R(2) = R(second H -- third H )
    !                            R(3) = R( third H -- first H )
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, az2, bet0, bet1, beta1, beta2, beta3, &
        & cart, cartnu, cd0, cd1, ce2, ce3, chh, ck0, ck1, cnvrtchi, cnvrtq, cnvrts, corrs, corrt, cphi3, damp, dampa, &
        & dampb, dc3, dcorrs, dcorrt, dcphi1, dcphi2, dcphi3, ddamp, deesdr, degsdr, deijdr, descart, dgscart, dijcart, &
        & dlep, dlep2, dq1, dq2, dq3, drho1, drho2, drho3, ds1, ds2, ds3, dv2, dv3, dva, engyes, engygs, engyij, ezero, &
        & h2r0, h2rm, h2rmt, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, pengyes, pengygs, pengyij, &
        & per, per2, qcoord, r, r12, r22, r32, rho, rhol, s, s2, sqrt3, t, va, vc, vii, viii, vleps, vleps2, xpar
    integer :: i, icartr, icouple, indes, indexes, irctnt, mder, msurf, nasurf, natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: natom  = 25
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    common /h3crcm/ corrs(3), dcorrs(3), &
        & corrt(3), dcorrt(3), &
        & damp(3, 3), ddamp(3, 3)
    !
    !**************************************************
    !                                                 *
    !  units conversion for anomalous adiabatic       *
    !  coupling pengyij                               *
    !                                                 *
    !**************************************************
    !
    common /ijconvert/cnvrtq, cnvrts, cnvrtchi, icouple
    !
    dimension dlep(3), dlep2(3), &
        & dv2(3), dv3(3), dva(3), dc3(3)
    !     DIMENSION R(3)
    !
    call cartou
    call carttor
    !
    ! Set up symmetry coordinates and their derivatives
    !
    call h3coor (r(1), r(2), r(3))
    !
    ! Set up 2-body correlation energies, dispersion damping terms, and
    !    their derivatives
    !
    call h3cor2
    !
    ! Get Leps potentials and derivatives
    !
    call h3leps (vleps, vleps2, dlep, dlep2)
    !
    ! Get VA term and its derivatives
    !
    call h3va (r(1), r(2), r(3), va, dva)
    !
    ! Get VII term and its derivatives
    !
    call h3vii (r(1), r(2), r(3), vii, dv2)
    !
    ! Get VIII term and its derivatives
    !
    call h3viii (viii, dv3)
    !
    ! Get 3-body correlation and its derivative
    !
    call h3cor3 (ce3, dc3)
    !
    ! 2-body correlation term
    !
    ce2 = corrs(1) + corrs(2) + corrs(3)
    !
    ! ground state
    !
    if (nasurf(1, 1) .gt. 0) then
        vc     = va + vii + viii + ce2 + ce3
        engygs = vleps + vc + ezero(1)
        if(nder.eq.1) then
            do 20 i = 1, 3
                t         = dva(i) + dv2(i) + dv3(i) + dcorrs(i) + dc3(i)
                degsdr(i) = ( dlep(i) + t)
            20 continue
        end if
    end if

    !
    ! first electronic state
    !
    if (nasurf(2, 2) .gt. 0) then
        indes         = 1
        vc            = va + vii + viii + ce2 + ce3
        engyes(indes) = vleps2 + vc + ezero(2)
        if(nder.eq.1) then
            do 30 i = 1, 3
                t                = dva(i) + dv2(i) + dv3(i) + dcorrs(i) + dc3(i)
                deesdr(i, indes) = ( dlep2(i) + t)
            30 continue
        end if
    end if

    !
    ! nonadiabatic coupling in an anomalous form where pengyij is treated as
    ! -  a vector with <psi|dpsi/dQ>,<psi|dpsi/ds>,<psi|dpsi/dchi>
    ! -  already transformed to the user's units
    !    (thus no RTOCART, DEDCOU calls via setting is nasurf(i,j))
    ! -  having no coded derivatives (no matter what nder is set to )
    !
    if (icouple .gt. 0) then
        call capf(r(1), r(2), r(3))
        pengyij(1) = engyij(1) * cnvrtq
        pengyij(2) = engyij(2) * cnvrts
        pengyij(3) = engyij(3) * cnvrtchi
    end if
    !
    call eunitzero
    if(nder.ne.0) then
        call rtocart
        call dedcou
    end if
    !
end subroutine pot

!
block data ptpacm
    !
    !***********************************************************************
    !     Data for the DMBE H3 surface.
    !***********************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, az2, bet0, bet1, beta1, beta2, beta3, &
        & cart, cartnu, cd0, cd1, chh, ck0, ck1, corrs, corrt, cphi3, damp, dampa, dampb, dcorrs, dcorrt, dcphi1, dcphi2, &
        & dcphi3, ddamp, dq1, dq2, dq3, drho1, drho2, drho3, ds1, ds2, ds3, ezero, h2r0, h2rm, h2rmt, h2ta, h2tb1, h2tb2, &
        & h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, per, per2, qcoord, r12, r22, r32, rho, rhol, s, s2, sqrt3, xpar
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt3cm/  ezero(isurf+1)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    common /h3crcm/ corrs(3), dcorrs(3), &
        & corrt(3), dcorrt(3), &
        & damp(3, 3), ddamp(3, 3)
    !
    data ezero(1) /0.174474112d0/
    data ezero(2) /0.174474112d0/
    data nasurf /1, 35*0/
    data nder /0/
    data nflag /1, 1, 15*0, 6, 0, 0/
    !
    data anuzero /0.0d0/
    data icartr, msurf, mder/3, 1, 1/
    data nulbl /25*0/
    data natoms /3/
    !
    data rhol /2.4848d0/
    data alph0, alph1, bet0, bet1 &
        & /25.9528d0, 1.1868d0, 15.7381d0, 0.09729d0/
    data alpha5, beta1, beta2, beta3 &
        & /8.2433033d-03, 0.53302897d0, 0.39156612d-01, 0.69996945d0/
    data alph2 /4.735364d-1/
    data al0, al1, al2, al3 &
        & /0.45024472d+01, -0.62467617d+01, 0.40966542d+01, &
        & 0.21813012d+01/
    data az2 /4.453649d-4/
    data cd0, cd1 /6.333404d-03, -1.726839d-03/
    data chh /  6.499027d0, 1.243991d+02, 3285.8d0/
    data ck0 / -1.372843d-01, -1.638459d-01, -1.973814d-01/
    data ck1 /  1.011204d0, 9.988099d-01, 9.399411d-01/
    data hfd, hfa1, hfa2, hfa3, hfgam &
        & /0.15796326d0, 2.1977034d0, 1.2932502d0, 0.64375666d0, &
        & 2.1835071d0/
    data h2rm, h2r0, h2rmt /1.401d0, 6.928203d0, 7.82d0/
    data h2ta, h2tb1, h2tb2, h2tb3, h2tb4 &
        & /0.448467d0, -0.056687d0, 0.246831d0, -0.018419d0, &
        & 0.000598d0/
    !      DATA PI / 3.14159265358979D0/
    data sqrt3 /1.73205080756887d0/
    data xpar /-0.9286579d-02, 0.2811592d-03, -0.4665659d-05, &
        & 0.2069800d-07, 0.2903613d+02, -0.2934824d+01, &
        & 0.7181886d0, -0.3753218d0, -0.1114538d+01, &
        & 0.2134221d+01, -0.4164343d0, 0.2022584d0, &
        & -0.4662687d-01, -0.4818623d+02, 0.2988468d0/
end

!
!> @brief Evaluate H3COOR for the H3_GEN_DMBE_1987 potential.
subroutine h3coor ( r1, r2, r3 )
    !
    !**********************************************************************
    !   Calculates D3H symmetry coordinates and derivatives
    !**********************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, arg, az2, bet0, bet1, bet2, beta, &
        & beta1, beta2, beta3, cart, cartnu, cd0, cd1, chh, ck0, ck1, corrs, corrt, cphi3, damp, dampa, dampb, dbt1, &
        & dbt2, dbt3, dcorrs, dcorrt, dcphi1, dcphi2, dcphi3, ddamp, deesdr, degsdr, deijdr, descart, dgm1, dgm2, dgm3, &
        & dgscart, dijcart, dphi, dq1, dq2, dq3, drho1, drho2, drho3, ds1, ds2, ds3, dsphi3, engyes, engygs, engyij, &
        & ezero, gam2, gamma, h2r0, h2rm, h2rmt, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, pengyes, &
        & pengygs, pengyij, per, per2, phi, phi3, qcoord, r, r1, r12, r2, r22, r3, r32, rho, rhol, rtan, s, s2, sphi3, &
        & sqrt3, t, t1, t12, t2, t3, theta, xpar
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    common /h3crcm/ corrs(3), dcorrs(3), &
        & corrt(3), dcorrt(3), &
        & damp(3, 3), ddamp(3, 3)
    !
    common/coupcm/phi, dphi(3), sphi3, dsphi3(3)
    !
    per  = r1+r2+r3
    per2 = per*per
    !
    !   QCOORD and its derivatives
    !
    r12    = r1*r1
    r22    = r2*r2
    r32    = r3*r3
    qcoord = r12+r22+r32
    dq1    = 2.0d0*r1
    dq2    = 2.0d0*r2
    dq3    = 2.0d0*r3
    !
    !   RHO and its derivatives
    !
    rho   = sqrt(qcoord/3.0d0)
    t     = 1.0d0/(6.0d0*rho)
    drho1 = t*dq1
    drho2 = t*dq2
    drho3 = t*dq3
    !
    !   S, CPHI3 (cos(phi3)), and their derivatives
    !
    gamma = 2.0d0*r12-r22-r32
    gam2  = gamma*gamma
    dgm1  = 2.0d0*dq1
    dgm2  = -dq2
    dgm3  = -dq3
    beta  = sqrt3*(r22-r32)
    bet2  = beta*beta
    dbt1  = 0.0d0
    dbt2  = sqrt3*dq2
    dbt3  = -sqrt3*dq3
    t12   = bet2 + gam2
    t1    = sqrt(t12)
    s     = t1/qcoord
    s2    = s*s
    if (s .eq. 0.0d0) then
        ds1 = 0.0d0
        ds2 = 0.0d0
        ds3 = 0.0d0
        !
        !   For S=0, CPHI3 and its derivative should not be
        !   used anywhere but set to zero anyway.
        !
        cphi3  = 0.0d0
        dcphi1 = 0.0d0
        dcphi2 = 0.0d0
        dcphi3 = 0.0d0
    else
        ds1    = s*((beta*dbt1+gamma*dgm1)/t12 - dq1/qcoord)
        ds2    = s*((beta*dbt2+gamma*dgm2)/t12 - dq2/qcoord)
        ds3    = s*((beta*dbt3+gamma*dgm3)/t12 - dq3/qcoord)
        t2     = 1.0d0/(t1*t12)
        cphi3  = gamma*(3.0d0*bet2 - gam2)*t2
        t3     = 3.0d0*beta*(3.0d0*gam2 - bet2)*t2/t12
        dcphi1 = t3*(gamma*dbt1 - beta*dgm1)
        dcphi2 = t3*(gamma*dbt2 - beta*dgm2)
        dcphi3 = t3*(gamma*dbt3 - beta*dgm3)
    end if
    if (nasurf(2, 2).eq.0.and.nflag(4).eq.0) go to 10
    arg = gamma/t1
    if(arg.le.-1.d0)arg=-1.d0
    if(arg.ge.1.d0)arg=1.d0
    theta     = acos(arg)
    phi       = pi-theta
    phi3      = 3.0d0*phi
    sphi3     = sin(phi3)
    rtan      = cphi3/sphi3
    dsphi3(1) = -rtan*dcphi1
    dsphi3(2) = -rtan*dcphi2
    dsphi3(3) = -rtan*dcphi3
    dphi(1)   = -dcphi1/(sphi3*3.d0)
    dphi(2)   = -dcphi2/(sphi3*3.d0)
    dphi(3)   = -dcphi3/(sphi3*3.d0)
    10 continue
end subroutine h3coor

!
!> @brief Evaluate H3COR2 for the H3_GEN_DMBE_1987 potential.
subroutine h3cor2
    !
    !***********************************************************************
    !     Calculates 2-body correlation energies for singlet and triplet
    !     states of H2, the damping factors for the dispersion terms, and
    !     their 1st derivatives
    !***********************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, az2, bet0, bet1, beta1, beta2, beta3, &
        & cart, cartnu, cd0, cd1, chh, ck0, ck1, corrs, corrt, cphi3, d, damp, dampa, dampb, dcorrs, dcorrt, dcphi1, &
        & dcphi2, dcphi3, dd, ddamp, deesdr, degsdr, deijdr, denom, descart, dgscart, dijcart, dpol, dq1, dq2, dq3, &
        & drho1, drho2, drho3, ds1, ds2, ds3, engyes, engygs, engyij, ezero, fnexp, h2r0, h2rm, h2rmt, h2ta, h2tb1, &
        & h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, pengyes, pengygs, pengyij, per, per2, pol, qcoord, r, r12, &
        & r22, r32, rho, rhol, rsave, s, s2, sqrt3, t, t1, t2, tt, tt1, tt2, x, xpar
    integer :: icartr, id, indexes, irctnt, j, mder, msurf, nasurf, &
        & natoms, nder, nexp, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !      DIMENSION RR(3)
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COR2)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    common /h3crcm/ corrs(3), dcorrs(3), &
        & corrt(3), dcorrt(3), &
        & damp(3, 3), ddamp(3, 3)
    !
    !
    !   Loop over three coordinates
    !
    do 10 j = 1, 3
        rsave     = r(j)
        corrs(j)  = 0.0d0
        dcorrs(j) = 0.0d0
        corrt(j)  = 0.0d0
        dcorrt(j) = 0.0d0
        t1        = 1.0d0/rsave
        t         = t1*t1
        t2        = t*t
        t1        = t2*t1
        nexp      = 4
        !
        !      Loop over terms in dispersion expansion
        !
        do 1 id = 1, 3
            nexp  = nexp + 2
            fnexp = dble(nexp)
            t2    = t2*t
            t1    = t1*t
            !
            !     singlet
            !
            !        Calculate damping function for the nth dispersion coefficient
            !        and its 1st derivative
            !
            denom = h2rm + 2.5d0*h2r0
            x     = 2.0d0*rsave/denom
            tt    = dampb(id)*x
            tt1   = dampa(id) + tt
            pol   = x*tt1
            dpol  = tt1 + tt
            tt    = exp(-pol)
            tt1   = 1.0d0 - tt
            tt2   = tt1**(nexp-1)
            d     = tt1*tt2
            dd    = 2.0d0*fnexp*dpol*tt*tt2/denom
            !
            !     store damping factors (including CHH coeff and 1/R**NEXP)
            !     FOR later use in computing the 3-body correlation energy.
            !
            damp(id, j)  = chh(id)*d*t2
            ddamp(id, j) = chh(id)*(dd*t2 - fnexp*d*t1)
            corrs(j)     = corrs(j) - damp(id, j)
            dcorrs(j)    = dcorrs(j) - ddamp(id, j)
            !
            !     triplet
            !        Calculate damping function for the nth dispersion
            !        coefficient and its 1st derivative
            !
            denom     = h2rmt + 2.5d0*h2r0
            x         = 2.0d0*rsave/denom
            tt        = dampb(id)*x
            tt1       = dampa(id) + tt
            pol       = x*tt1
            dpol      = tt1 + tt
            tt        = exp(-pol)
            tt1       = 1.0d0 - tt
            tt2       = tt1**(nexp-1)
            d         = tt1*tt2
            dd        = 2.0d0*fnexp*dpol*tt*tt2/denom
            corrt(j)  = corrt(j) - chh(id)*d*t2
            dcorrt(j) = dcorrt(j) - chh(id)*(dd*t2 - fnexp*d*t1)
        1 continue
    10 continue
end subroutine h3cor2

!
!> @brief Evaluate H3COR3 for the H3_GEN_DMBE_1987 potential.
subroutine h3cor3 ( ce3, dc3 )
    !
    !***********************************************************************
    !     Calculates 3-body correlation energy and its 1st derivatives
    !***********************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, az2, bet0, bet1, beta1, beta2, beta3, &
        & cart, cartnu, cd0, cd1, ce3, chh, ck0, ck1, corrs, corrt, cphi3, damp, dampa, dampb, dc3, dcorrs, dcorrt, &
        & dcphi1, dcphi2, dcphi3, ddamp, deesdr, degsdr, deijdr, descart, dgscart, dijcart, dq1, dq2, dq3, drho1, drho2, &
        & drho3, ds1, ds2, ds3, engyes, engygs, engyij, ezero, g, gd, h, h2r0, h2rm, h2rmt, h2ta, h2tb1, h2tb2, h2tb3, &
        & h2tb4, hd, hfa1, hfa2, hfa3, hfd, hfgam, hysec, hytan, pengyes, pengygs, pengyij, per, per2, qcoord, r, r12, &
        & r22, r32, rho, rhol, rr, s, s2, sgnt, sqrt3, t, t1, t2, xpar
    integer :: i, i2, i3, icartr, id, indexes, irctnt, mder, msurf, nasurf, natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    dimension dc3(3), g(3), gd(3), h(3), hd(3)
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    common /h3crcm/ corrs(3), dcorrs(3), &
        & corrt(3), dcorrt(3), &
        & damp(3, 3), ddamp(3, 3)
    !
    ce3    = 0.0d0
    dc3(1) = 0.0d0
    dc3(2) = 0.0d0
    dc3(3) = 0.0d0
    !
    !   Loop over terms in dispersion expansion; dispersion damping
    !   FACTORS are computed in H3COR2 and passed through COMMON /H3CRCM/.
    !
    do 30 id = 1, 3
        !
        !   Loop over 3 coordinates
        !
        do 10 i = 1, 3
            rr = r(i)
            !
            !     Calculate G function and its 1st derivative
            !
            t     = ck0(id)*exp(-ck1(id)*(rr-h2rm))
            g(i)  = 1.0d0 + t
            gd(i) = -ck1(id)*t
            !
            !     Calculate H function and its 1st derivative.
            !
            t    = ck1(id)*rr
            sgnt = 1.0d0
            if (t .lt. 0.0d0) then
                t    = -t
                sgnt = -1.0d0
            end if
            t     = exp(-t)
            t2    = t*t
            t1    = 1.0d0/(1.0d0+t2)
            hysec = 2.0d0*t*t1
            hytan = sgnt*(1.0d0-t2)*t1
            t1    = hytan**5
            h(i)  = hytan*t1
            hd(i) = 6.0d0*ck1(id)*t1*hysec*hysec
        10 continue
        do 20 i = 1, 3
            i2      = mod(i, 3) + 1
            i3      = mod(i+1, 3) + 1
            t       = 1.0d0 - 0.5d0*(g(i2)*h(i3) + g(i3)*h(i2))
            ce3     = ce3 + t*damp(id, i)
            dc3(i)  = dc3(i) + t*ddamp(id, i)
            dc3(i2) = dc3(i2) - 0.5d0*(gd(i2)*h(i3)+g(i3)*hd(i2))*damp(id, i)
            dc3(i3) = dc3(i3) - 0.5d0*(g(i2)*hd(i3)+gd(i3)*h(i2))*damp(id, i)
        20 continue
    30 continue
end subroutine h3cor3

!
!> @brief Evaluate H3LEPS for the H3_GEN_DMBE_1987 potential.
subroutine h3leps ( vleps, vleps2, dlep, dlep2 )
    !
    !***********************************************************************
    !     Calculates 3-body extended-Hartree-Fock energy defined by a
    !     LEPS-type function.
    !***********************************************************************
    !
    !   VLEPS  IS LEPS LOWER SURFACE.
    !   VLEPS2 IS LEPS UPPER SURFACE.
    !   DLEP(3) ARE LEPS LOWER SURFACE DERIVATIVES
    !   DLEP2(3) "   "   UPPER    "       "
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, at, az2, bet0, bet1, beta1, beta2, &
        & beta3, cart, cartnu, cd0, cd1, chh, ck0, ck1, corrs, corrt, cphi3, damp, dampa, dampb, dat, dcorrs, dcorrt, &
        & dcphi1, dcphi2, dcphi3, ddamp, deesdr, degsdr, deijdr, descart, desdr, df, dgscart, dijcart, dlep, dlep2, dq, &
        & dq1, dq2, dq3, dr, drho1, drho2, drho3, ds1, ds2, ds3, dwe, dxj, engyes, engygs, engyij, es, et, exch, ezero, &
        & f, f1, f2, f3, fterm, h2r0, h2rm, h2rmt, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, omf, &
        & pengyes, pengygs, pengyij, per, per2, qcoord, qsum, r, r12, r22, r32, rho, rhol, rr, rri, s, s2, sqrt3, t, t1, &
        & t2, tdc, tdq, tds, texch, vleps, vleps2, we, xj, xpar, xterm, xterm1, xterm2, xterm3
    integer :: i, icartr, indexes, irctnt, j, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    common /h3crcm/ corrs(3), dcorrs(3), &
        & corrt(3), dcorrt(3), &
        & damp(3, 3), ddamp(3, 3)
    !
    dimension dlep(3), dlep2(3)
    dimension df(3), dq(3), xj(3), &
        & dxj(3, 3)
    !
    !  Calculate switching term for the Hartree-Fock component of the
    !  diatomic triplet state function.  The notation of Thompson et al.
    !  (J. Chem. Phys., 82,5597,1985; and references therein) is used
    !  throughout.
    !
    t1 = -az2*qcoord
    if (s .eq. 0.0d0) then
        f     = exp(t1*qcoord)
        t1    = 2.0d0*t1*f
        df(1) = t1*dq1
        df(2) = t1*dq2
        df(3) = t1*dq3
    else
        t2    = 1.0d0 + s*s2*cphi3
        f     = exp(t1*qcoord*t2)
        t1    = t1*f
        tdq   = 2.0d0*t2
        t2    = s2*qcoord
        tds   = 3.0d0*cphi3*t2
        tdc   = s*t2
        df(1) = t1*(tdq*dq1 + tds*ds1 + tdc*dcphi1)
        df(2) = t1*(tdq*dq2 + tds*ds2 + tdc*dcphi2)
        df(3) = t1*(tdq*dq3 + tds*ds3 + tdc*dcphi3)
    end if
    omf = 1.0d0 - f
    !
    qsum = 0.0d0
    do 10 i = 1, 3
        dq(i) = 0.0d0
    10 continue
    !
    !  Loop over 3 coordinates
    !
    do 50 i = 1, 3
        rr  = r(i)
        rri = 1.0d0/rr
        !
        !     Calculate extended-Hartree-Fock curve for ground singlet state of
        !     H2 [A.J.C. Varandas & J.D. Silva, J. Chem. Soc. Faraday II
        !     (submitted)] and 1st derivative of 2-body extended-Hartree-Fock
        !     curve for the ground-singlet state of H2.
        !
        dr    = rr - h2rm
        t     = exp(-hfgam*dr)
        es    = -hfd*(1.0d0 + dr*(hfa1 + dr*(hfa2 + dr*hfa3)))*t
        desdr = -hfgam*es - hfd*(hfa1 + dr*(2.0d0*hfa2 + 3.0d0*dr*hfa3))*t
        !
        !     Compute the HFACE (Hartree-Fock-approximate correlation energy)
        !     potential for the lowest triplet state of H2 [A.J.C. Varandas &
        !     J. Brandao Mol. Phys.,45,1982,857] without the 2-body correlation
        !     energy.
        !
        t   = rr*(h2tb1 + rr*(h2tb2 + rr*(h2tb3 + rr*h2tb4)))
        at  = h2ta*exp(-t)*rri
        t   = h2tb1 + rr*(2.0d0*h2tb2 + rr*(3.0d0*h2tb3 + rr*4.0d0*h2tb4))
        dat = -at*(t + rri)
        !
        !     Add in triplet and subtract singlet 2-body correlation terms
        !
        at  = at + corrt(i) - corrs(i)
        dat = dat + dcorrt(i) - dcorrs(i)
        !
        !     Calculate effective diatomic triplet state curve.
        !
        t   = exp(-al3*rr)
        we  = (al0 + rr*(al1 + rr*al2))*t
        dwe = (al1 + 2.0d0*al2*rr)*t - al3*we
        !
        !     Triplet energy
        !
        et = f*we + omf*at
        !
        !     Contribution to coulomb term
        !
        qsum = qsum + 0.5d0*(es + et)
        !
        !     Contribution to exchange term
        !
        xj(i) = 0.5d0*(es - et)
        !
        !     Derivatives of coulomb and exchange parts
        !
        fterm = we - at
        do 20 j = 1, 3
            t         = 0.5d0*fterm*df(j)
            dq(j)     = dq(j) + t
            dxj(i, j) = -t
        20 continue
        t         = f*dwe + omf*dat
        dq(i)     = dq(i) + 0.5d0*(desdr + t)
        dxj(i, i) = dxj(i, i) + 0.5d0*(desdr - t)
    50 continue
    !
    f1     = (xj(1) - xj(2))**2
    f2     = (xj(2) - xj(3))**2
    f3     = (xj(3) - xj(1))**2
    exch   = sqrt(0.5d0*(f1 + f2 + f3))
    vleps  = qsum - exch
    vleps2 = qsum + exch
    xterm1 = 2.0d0*xj(1) - xj(2) - xj(3)
    xterm2 = 2.0d0*xj(2) - xj(1) - xj(3)
    xterm3 = 2.0d0*xj(3) - xj(1) - xj(2)
    do 80 i = 1, 3
        texch = exch
        if(texch.le.0.0d0)texch = 1.0d0
        xterm    = 0.5d0/texch*(xterm1*dxj(1, i)  + xterm2*dxj(2, i) + xterm3*dxj(3, i))
        dlep(i)  = dq(i) - xterm
        dlep2(i) = dq(i) + xterm
    80 continue
end subroutine h3leps

!
!> @brief Evaluate H3VA for the H3_GEN_DMBE_1987 potential.
subroutine h3va ( r1, r2, r3, va, dva )
    !
    !***********************************************************************
    !     Calculates Va correction energy and its 1st derivative.
    !***********************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, az2, bet0, bet1, beta1, beta2, beta3, &
        & cart, cartnu, cd0, cd1, chh, ck0, ck1, cphi3, dampa, dampb, dcphi1, dcphi2, dcphi3, deesdr, degsdr, deijdr, &
        & descart, dgscart, dijcart, dq1, dq2, dq3, drho1, drho2, drho3, ds1, ds2, ds3, dv, dva, engyes, engygs, engyij, &
        & expva, ezero, h2r0, h2rm, h2rmt, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, pengyes, &
        & pengygs, pengyij, per, per2, qcoord, r, r1, r12, r2, r22, r3, r32, rho, rhol, s, s2, sqrt3, t, t1, t2, t3, t4, &
        & va, xpar
    integer :: icartr, indexes, irctnt, j, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    dimension dva(3)
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
    !     +                CORRT(3),DCORRT(3),
    !     +                DAMP(3,3),DDAMP(3,3)
    !
    t     = alpha5*per2
    expva = exp(-t*per)
    if (expva .eq. 0.0d0) then
        va     = 0.0d0
        dva(1) = 0.0d0
        dva(2) = 0.0d0
        dva(3) = 0.0d0
    else
        va = 0.0d0
        dv = 0.0d0
        t3 = (r1-r2)*(r2-r3)*(r3-r1)
        t1 = t3*t3
        t2 = 1.0d0
        t4 = 2.0d0
        do 1 j=1, 4
            t2 = t2*t1
            va = va + xpar(j)*t2
            dv = dv + t4*xpar(j)*t3
            t3 = t3*t1
            t4 = t4 + 2.0d0
        1 continue
        va     = va*expva
        t1     = dv*expva
        t2     = 3.0d0*t*va
        dva(1) = t1*(r2-r3)*(r3+r2-2.0d0*r1) - t2
        dva(2) = t1*(r3-r1)*(r1+r3-2.0d0*r2) - t2
        dva(3) = t1*(r1-r2)*(r1+r2-2.0d0*r3) - t2
    end if
end subroutine h3va

!
!> @brief Evaluate H3VII for the H3_GEN_DMBE_1987 potential.
subroutine h3vii ( r1, r2, r3, e, dv2 )
    !
    !***********************************************************************
    !     Calculates VII correction energy and its 1st derivative.
    !***********************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, az2, bet0, bet1, beta1, beta2, beta3, &
        & cart, cartnu, cd0, cd1, chh, ck0, ck1, cos1, cos2, cos3, cphi3, dampa, dampb, dcphi1, dcphi2, dcphi3, deesdr, &
        & degsdr, deijdr, descart, dgscart, dijcart, dq1, dq2, dq3, dr12, dr23, dr31, drho1, drho2, drho3, ds1, ds2, ds3, &
        & dv2, e, eb1, eb1p1, eb1p2, eb1p3, eb1pr, eb1pwb, eb1t, eb2, eb2p1, eb2p2, eb2p3, eb2pr, eb2pwb, eb3t, eb4, &
        & eb4a, eb4ap1, eb4ap2, eb4ap3, eb4apr, eb4apw, eb4b, eb4bp1, eb4bp2, eb4bp3, eb4bpr, eb4bpw, eb4p1, eb4p2, &
        & eb4p3, engyes, engygs, engyij, eq, eqp1, eqp2, eqp3, exp1, exp2, exp3, ezero, h2r0, h2rm, h2rmt, h2ta, h2tb1, &
        & h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, pengyes, pengygs, pengyij, per, per2, qcoord, r, r1, r12, &
        & r1i, r2, r22, r2i, r3, r32, r3i, rho, rhol, ri, s, s2, sqrt3, t1, t2, t3, wb, wb1p, wb2, wb2p, wb3p, xpar
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    dimension dv2(3)
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
    !     +                CORRT(3),DCORRT(3),
    !     +                DAMP(3,3),DDAMP(3,3)
    !
    !   Compute B1 function.
    !
    r1i  = 1.0d0/r1
    r2i  = 1.0d0/r2
    r3i  = 1.0d0/r3
    cos1 = 0.5d0*r2i*r3i*(r12-r22-r32)
    cos2 = 0.5d0*r1i*r3i*(r22-r12-r32)
    cos3 = 0.5d0*r1i*r2i*(r32-r12-r22)
    wb   = 1.0d0+cos1+cos2+cos3
    wb2  = wb*wb
    !
    !   WB derivatives
    !
    wb1p = (r1*r3i - 1.0d0)*r2i - r3i - (cos2+cos3)*r1i
    wb2p = (r2*r3i - 1.0d0)*r1i - r3i - (cos1+cos3)*r2i
    wb3p = (r3*r2i - 1.0d0)*r1i - r2i - (cos1+cos2)*r3i
    !
    !   EB1 term
    !
    exp1 = exp(-beta1*per)
    exp3 = exp(-beta3*per)
    eb1t = (xpar(5) + xpar(6)*per)*exp1
    eb3t = (xpar(14) + xpar(15)*per2)*exp3
    eb1  = wb*(eb1t + eb3t)
    !
    !   EB1 derivatives
    !
    eb1pr  = wb*(-beta1*eb1t - beta3*eb3t + xpar(6)*exp1 + 2.0d0*per*xpar(15)*exp3)
    eb1pwb = eb1t + eb3t
    eb1p1  = eb1pwb*wb1p + eb1pr
    eb1p2  = eb1pwb*wb2p + eb1pr
    eb1p3  = eb1pwb*wb3p + eb1pr
    !
    !   EB2 term
    !
    t1   = beta2*per
    exp2 = exp(-t1*per)
    eb2  = wb2*(xpar(7) + wb*(xpar(8) + wb*xpar(9)))*exp2
    !
    !   EB2 derivatives
    !
    eb2pwb = wb*(2.0d0*xpar(7) + wb*(3.0d0*xpar(8) + wb*4.0d0*xpar(9)) )*exp2
    eb2pr  = -2.0d0*t1*eb2
    eb2p1  = eb2pwb*wb1p + eb2pr
    eb2p2  = eb2pwb*wb2p + eb2pr
    eb2p3  = eb2pwb*wb3p + eb2pr
    !
    !   EB4 term
    !      EB4A
    !
    t2   = xpar(10)*exp1
    t3   = wb*xpar(11)*exp2
    eb4a = wb*(t2 + t3)
    !
    !      EB4A derivatives
    !
    eb4apw = t2 + 2.0d0*t3
    eb4apr = -wb*(beta1*t2 + 2.0d0*t1*t3)
    eb4ap1 = eb4apw*wb1p + eb4apr
    eb4ap2 = eb4apw*wb2p + eb4apr
    eb4ap3 = eb4apw*wb3p + eb4apr
    !
    !      EB4B
    !
    t2   = xpar(12)*exp1
    t3   = xpar(13)*exp2
    eb4b = wb*(t2 + t3)
    !
    !      EB4B derivatives
    !
    eb4bpw = t2 + t3
    eb4bpr = -wb*(beta1*t2 + 2.0d0*t1*t3)
    eb4bp1 = eb4bpw*wb1p + eb4bpr
    eb4bp2 = eb4bpw*wb2p + eb4bpr
    eb4bp3 = eb4bpw*wb3p + eb4bpr
    dr12   = r1 - r2
    dr23   = r2 - r3
    dr31   = r3 - r1
    eq     = dr12*dr12 + dr23*dr23 + dr31*dr31
    eqp1   = 2.0d0*(dr12 - dr31)
    eqp2   = 2.0d0*(-dr12 + dr23)
    eqp3   = 2.0d0*(-dr23 + dr31)
    ri     = r1i + r2i + r3i
    eb4    = eb4a*ri + eb4b*eq
    !
    !   EB4 derivatives
    !
    eb4p1  = eb4ap1*ri - eb4a/r12 + eb4bp1*eq + eb4b*eqp1
    eb4p2  = eb4ap2*ri - eb4a/r22 + eb4bp2*eq + eb4b*eqp2
    eb4p3  = eb4ap3*ri - eb4a/r32 + eb4bp3*eq + eb4b*eqp3
    e      = eb1 + eb2 + eb4
    dv2(1) = eb1p1 + eb2p1 + eb4p1
    dv2(2) = eb1p2 + eb2p2 + eb4p2
    dv2(3) = eb1p3 + eb2p3 + eb4p3
end subroutine h3vii

!
!> @brief Evaluate H3VIII for the H3_GEN_DMBE_1987 potential.
subroutine h3viii ( viii, dv3 )
    !
    !**********************************************************************
    !     Calculates VIII correction energy and it 1st derivatives
    !***********************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, az2, bet0, bet1, beta1, beta2, beta3, &
        & cart, cartnu, cd0, cd1, chh, ck0, ck1, cphi3, dampa, dampb, dcphi1, dcphi2, dcphi3, deesdr, degsdr, deijdr, &
        & descart, dgscart, dijcart, dq1, dq2, dq3, drho1, drho2, drho3, ds1, ds2, ds3, dv3, engyes, engygs, engyij, &
        & expv3, ezero, h2r0, h2rm, h2rmt, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, pengyes, &
        & pengygs, pengyij, per, per2, qcoord, r, r12, r22, r32, rho, rhol, s, s2, sqrt3, t1, t2, t3, t4, t5, tdcph, &
        & tdrho, tds, viii, xpar
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    dimension dv3(3)
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
    !     +                CORRT(3),DCORRT(3),
    !     +                DAMP(3,3),DDAMP(3,3)
    !
    if (s .eq. 0.0d0) then
        viii   = 0.0d0
        dv3(1) = 0.0d0
        dv3(2) = 0.0d0
        dv3(3) = 0.0d0
    else
        t1    = rho - rhol
        t5    = alph2*t1
        expv3 = exp(-t5*t1)
        if (expv3 .eq. 0.0d0) then
            viii   = 0.0d0
            dv3(1) = 0.0d0
            dv3(2) = 0.0d0
            dv3(3) = 0.0d0
        else
            t1     = s2*cphi3
            t2     = 1.0d0 + s*t1
            t3     = s2*t2
            t4     = cd0 + cd1*rho
            viii   = t3*t4*expv3
            tds    = (2.0d0*s*t2 + 3.0d0*s2*t1)*t4
            tdcph  = s2*s*s2*t4
            tdrho  = t3*(cd1 - 2.0d0*t5*t4)
            dv3(1) = (tds*ds1 + tdcph*dcphi1 + tdrho*drho1)*expv3
            dv3(2) = (tds*ds2 + tdcph*dcphi2 + tdrho*drho2)*expv3
            dv3(3) = (tds*ds3 + tdcph*dcphi3 + tdrho*drho3)*expv3
        end if
    end if
end subroutine h3viii

!
!> @brief Evaluate CAPF for the H3_GEN_DMBE_1987 potential.
subroutine capf( r1, r2, r3 )
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, arg, az2, bet0, bet1, beta1, beta2, &
        & beta3, brak, capk, cart, cartnu, cd0, cd1, chh, ck0, ck1, cphi3, dampa, dampb, darg, dcapk, dcphi, dcphi1, &
        & dcphi2, dcphi3, deesdr, degsdr, deijdr, denom, descart, dfu, dfz, dgscart, dgz, dijcart, dphi, dq1, dq2, dq3, &
        & drho1, drho2, drho3, ds, ds1, ds2, ds3, dsphi3, dx, engyes, engygs, engyij, ezero, f, fu, fz, gz, h2r0, h2rm, &
        & h2rmt, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, pengyes, pengygs, pengyij, per, per2, &
        & phi, qcoord, r, r1, r12, r2, r22, r3, r32, rcapf, rho, rhol, s, s2, sphi3, sqrt3, x, x2, xi, xi3, xpar
    integer :: i, icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
    !     +                CORRT(3),DCORRT(3),
    !     +                DAMP(3,3),DDAMP(3,3)
    !
    common/coupcm/phi, dphi(3), sphi3, dsphi3(3)
    dimension dx(3), dcapk(4), dfz(3), &
        & dfu(3), dgz(3), f(3, jsurf), rcapf(3), &
        & ds(3), dcphi(3)
    rcapf(1) = r1
    rcapf(2) = r2
    rcapf(3) = r3
    call h3coor(rcapf(1), rcapf(2), rcapf(3))
    !
    ds(1)    = ds1
    ds(2)    = ds2
    ds(3)    = ds3
    dcphi(1) = dcphi1
    dcphi(2) = dcphi2
    dcphi(3) = dcphi3
    !
    x   = sqrt(qcoord/3.0d0)
    xi  = 1.0d0/x
    xi3 = x/3.0d0
    x2  = x*x
    call kappa(x, capk, dcapk)
    fz    = 0.375d0*x*dcapk(1)
    fu    = 0.03515625d0*(x*dcapk(1)-x2*dcapk(2)+x2*xi3*dcapk(3))
    gz    = 0.046875d0*(-x*dcapk(1)+x2*dcapk(2))
    denom = fz+gz*cphi3*s+fu*s2
    arg   = gz*sphi3*s/denom
    !     FCT=PHI-ATAN(ARG)
    do 10 i=1, 3
        dx(i)     = rcapf(i)*xi
        dfz(i)    = 0.125d0*dx(i)*(dcapk(1)+x*dcapk(2))
        dfu(i)    = 0.01171875d0*dx(i)*(dcapk(1)-x*dcapk(2)+x2*xi3*dcapk(4))
        dgz(i)    = 0.015625d0*dx(i)*(-dcapk(1)+x*dcapk(2)+x2*dcapk(3))
        brak      = s*dgz(i)+gz*ds(i)
        darg      = -arg*(dfz(i)+cphi3*brak+s*(gz*dcphi(i)+s*dfu(i)+ 2.0d0*fu*ds(i)))
        darg      = (darg+gz*s*dsphi3(i)+sphi3*brak)/denom
        engyij(i) = dphi(i)-1.0d0/(1.0d0+arg**2)*darg
    10 continue
end subroutine capf

!
!> @brief Evaluate H2COR2 for the H3_GEN_DMBE_1987 potential.
subroutine h2cor2 ( rr )
    !
    !*****************************************************************
    !     CALCULATES 1-4 DERIVATIVES OF 2-BODY CORRELATION TERMS FOR
    !     BOTH SINGLET AND TRIPLET STATES AND RETURNS THEIR DIFFERENCE
    !******************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, az2, bet0, bet1, beta1, beta2, beta3, &
        & cart, cartnu, cd0, cd1, chh, ck0, ck1, corrs, corrt, cphi3, d, d2, d3, d4, dampa, dampb, dcor2, dcorrs, dcorrt, &
        & dcphi1, dcphi2, dcphi3, dd, deesdr, degsdr, deijdr, denom, descart, dgscart, dijcart, dpol, dq1, dq2, dq3, &
        & drho1, drho2, drho3, ds1, ds2, ds3, engyes, engygs, engyij, ezero, fnexp, h2r0, h2rm, h2rmt, h2ta, h2tb1, &
        & h2tb2, h2tb3, h2tb4, hfa1, hfa2, hfa3, hfd, hfgam, parb, pengyes, pengygs, pengyij, per, per2, pol, qcoord, r, &
        & r12, r22, r32, rho, rhol, rr, rsave, s, s2, sqrt3, t, t1, t2, tem2, tem3, tem4, tt, tt1, tt2, x, xpar, y
    integer :: icartr, id, indexes, irctnt, j, mder, msurf, nasurf, &
        & natoms, nder, nexp, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
    !     +                CORRT(3),DCORRT(3),
    !     +                DAMP(3,3),DDAMP(3,3)
    !
    common/dmpcm/dpol, y, tt, tt1, parb
    common/cor2cm/corrs, corrt, dcor2(4)
    rsave  = rr
    corrs  = 0.0d0
    dcorrs = 0.0d0
    corrt  = 0.0d0
    dcorrt = 0.0d0
    do 11 j=1, 4
    11 dcor2(j)=0.0d0
    t1   = 1.0d0/rsave
    t    = t1*t1
    t2   = t*t
    t1   = t2*t1
    nexp = 4
    !      Loop over terms in dispersion expansion
    do 1 id = 1, 3
        nexp  = nexp + 2
        fnexp = dble(nexp)
        t2    = t2*t
        t1    = t1*t
        !
        !******SINGLET
        !
        denom  = h2rm + 2.5d0*h2r0
        y      = 2.0d0/denom
        x      = rsave*y
        parb   = dampb(id)
        tt     = dampb(id)*x
        tt1    = dampa(id) + tt
        pol    = x*tt1
        dpol   = tt1 + tt
        tt     = exp(-pol)
        tt1    = 1.0d0 - tt
        tt2    = tt1**(nexp-1)
        d      = tt1*tt2
        dd     = 2.0d0*fnexp*dpol*tt*tt2/denom
        corrs  = corrs - chh(id)*d*t2
        dcorrs = dcorrs - chh(id)*(dd*t2-fnexp*d*t1)
        call damp(rsave, nexp, d2, d3, d4)
        tem2 = d2
        tem3 = d3
        tem4 = d4
        !
        !*****TRIPLET
        !
        denom  = h2rmt + 2.5d0*h2r0
        y      = 2.0d0/denom
        x      = rsave*y
        tt     = dampb(id)*x
        tt1    = dampa(id) + tt
        pol    = x*tt1
        dpol   = tt1 + tt
        tt     = exp(-pol)
        tt1    = 1.0d0 - tt
        tt2    = tt1**(nexp-1)
        d      = tt1*tt2
        dd     = 2.0d0*fnexp*dpol*tt*tt2/denom
        corrt  = corrt - chh(id)*d*t2
        dcorrt = dcorrt - chh(id)*(dd*t2 - fnexp*d*t1)
        call damp(rsave, nexp, d2, d3, d4)
        dcor2(2) = dcor2(2)+chh(id)*(d2-tem2)
        dcor2(3) = dcor2(3)+chh(id)*(d3-tem3)
        dcor2(4) = dcor2(4)+chh(id)*(d4-tem4)
        3 continue
    1 continue
    dcor2(1) = dcorrt-dcorrs
    10 continue
end subroutine h2cor2

!
!> @brief Evaluate DAMP for the H3_GEN_DMBE_1987 potential.
subroutine damp( rmp, n, d2, d3, d4 )
    !
    !**********************************************************************
    !      CALCULATES 2ND,3RD AND 4TH DERIVATIVES OF 2-BODY CORRELATION
    !      ENERGIES FOR SINGLET AND TRIPLET STATES FOR USE IN NON-ADIABATIC
    !      COUPLINGS OF H3 SURFACES 1 AND 2. SENT TO SUBROUTINE H3COR2
    !**********************************************************************
    !
    implicit none
    double precision :: anuzero, b, bm2, bm3, bm4, cart, cartnu, d1b, d2, d2b, d3, d3b, d4, d4b, dampb, deesdr, degsdr, &
        & deijdr, descart, dgscart, dijcart, dpol, dpol2, engyes, engygs, engyij, ezero, fr, pengyes, pengygs, pengyij, &
        & r, rmp, rr, rrn, t, texp, th, tm, tm2, tm3, tp, tp2, tp3, tw, y
    integer :: icartr, indexes, irctnt, mder, msurf, n, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    !      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,
    !     +                AL0,AL1,AL2,AL3,AZ2,
    !     +                BETA1,BETA2,BETA3,BET0,BET1,
    !     +                CD0,CD1,CHH(3),CK0(3),CK1(3),
    !     +                DAMPA(3),DAMPB(3),
    !     +                HFD,HFA1,HFA2,HFA3,HFGAM,
    !     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,
    !     +                RHOL,SQRT3,XPAR(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    !      COMMON /H3COCM/ PER,PER2,R12,R22,R32,
    !     +                QCOORD,DQ1,DQ2,DQ3,
    !     +                RHO,DRHO1,DRHO2,DRHO3,
    !     +                S,S2,DS1,DS2,DS3,
    !     +                CPHI3,DCPHI1,DCPHI2,DCPHI3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
    !     +                CORRT(3),DCORRT(3),
    !     +                DAMP(3,3),DDAMP(3,3)
    !
    common/dmpcm/dpol, y, texp, b, dampb
    dpol2 = dpol**2
    d1b   = y*dpol*texp
    d2b   = y**2*(2.0d0*dampb-dpol2)*texp
    d3b   = y**3*dpol*(dpol2-6.0d0*dampb)*texp
    d4b   = y**4*(-12.d0*dampb*dampb+dpol2*(12.d0*dampb-dpol2))*texp
    t     = dble(n)
    tm    = dble(n-1)
    tm2   = dble(n-2)
    tm3   = dble(n-3)
    tp    = dble(n+1)
    tp2   = dble(n+2)
    tp3   = dble(n+3)
    th    = 3.0d0
    tw    = 2.0d0
    fr    = 4.0d0
    rr    = 1.0d0/rmp
    rrn   = rr**n
    bm4   = b**(n-4)
    bm3   = b*bm4
    bm2   = b*bm3
    d2    = t*bm2*(-tp*(b*rr)**2+tw*t*b*rr*d1b-tm*d1b**2-b*d2b)*rrn
    d3=t*bm3*(tp*tp2*(b*rr)**3-th*t*tp*(b*rr)**2*d1b+th*tm*t*b*d1b**2* rr-tm2*tm*d1b**3+th*t*rr*b**2*d2b-th*tm*b*d1b* &
        & d2b-b**2*d3b)*rrn
    d4=t*bm4*(-tp*tp2*tp3*(b*rr)**4+fr*t*tp*tp2*(b*rr)**3*d1b-tw*th*tm *t*tp*(b*rr*d1b)**2+fr*tm2*tm*t*b*rr*d1b**3-tw* &
        & th*t*tp*b*(b*rr)** 2*d2b+th*fr*tm*t*b*b*rr*d1b*d2b-tw*th*tm2*tm*b*d1b*d1b*d2b+fr*t* rr*b**3*d3b-tm3*tm2*tm* &
        & d1b**4-th*tm*(b*d2b)**2-fr*tm*b*b*d1b *d3b-b**3*d4b)*rrn
end subroutine damp

!
!> @brief Evaluate KAPPA for the H3_GEN_DMBE_1987 potential.
subroutine kappa( qlc, capk, dcapk )
    !
    !*****************************************************************
    !     CALCULATES FUNCTION CAPK AND ITS FIRST ROUR DERIVATIVES
    !*****************************************************************
    !
    implicit none
    double precision :: al0, al1, al2, al3, alph0, alph1, alph2, alpha5, anuzero, at, az2, bet0, bet1, beta1, beta2, &
        & beta3, capk, cart, cartnu, cd0, cd1, chh, ck0, ck1, corrs, corrt, cphi3, d1uehf, d1vhf, d1w, d2uehf, d2vhf, &
        & d2w, d3uehf, d3vhf, d3w, d4uehf, d4vhf, d4w, dampa, dampb, dat, dcapk, dcor2, dcphi1, dcphi2, dcphi3, deesdr, &
        & degsdr, deijdr, descart, desdr, dgscart, dijcart, dq1, dq2, dq3, dr, drho1, drho2, drho3, ds1, ds2, ds3, dwe, &
        & engyes, engygs, engyij, es, ezero, fzero, h2r0, h2rm, h2rmt, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, hfa1, hfa2, &
        & hfa3, hfd, hfgam, o1, o2, omfz, p1, p2, p3, pengyes, pengygs, pengyij, per, per2, q1, q2, q3, qcoord, qlc, r, &
        & r12, r22, r32, rho, rhol, rr, rri, s, s2, sqrt3, t, t1, we, xpar
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    !***********************************************************************
    !
    !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
    !
    !***********************************************************************
    !
    common /h3dmcm/ alph2, alpha5, alph0, alph1, &
        & al0, al1, al2, al3, az2, &
        & beta1, beta2, beta3, bet0, bet1, &
        & cd0, cd1, chh(3), ck0(3), ck1(3), &
        & dampa(3), dampb(3), &
        & hfd, hfa1, hfa2, hfa3, hfgam, &
        & h2rm, h2rmt, h2r0, h2ta, h2tb1, h2tb2, h2tb3, h2tb4, &
        & rhol, sqrt3, xpar(15)
    !
    !******************************************************
    !                                                     *
    !  Common block for coordinates (computed in H3COOR)  *
    !                                                     *
    !******************************************************
    !
    common /h3cocm/ per, per2, r12, r22, r32, &
        & qcoord, dq1, dq2, dq3, &
        & rho, drho1, drho2, drho3, &
        & s, s2, ds1, ds2, ds3, &
        & cphi3, dcphi1, dcphi2, dcphi3
    !
    !**************************************************
    !                                                 *
    !  Common block for 2-body correlation energies   *
    !  and damping factors (computed in H3COR2)       *
    !                                                 *
    !**************************************************
    !
    !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
    !     +                CORRT(3),DCORRT(3),
    !     +                DAMP(3,3),DDAMP(3,3)
    !
    !
    !
    common/cor2cm/corrs, corrt, dcor2(4)
    dimension dcapk(4)
    call h2cor2(qlc)
    t1    = -az2*qcoord
    fzero = exp(t1*qcoord)
    rr    = qlc
    rri   = 1.0d0/rr
    !
    !*****CALCULATE 1-4 DERIVATIVES OF VEHF
    !
    dr       = rr - h2rm
    t        = exp(-hfgam*dr)
    es       = -hfd*(1.0d0 + dr*(hfa1 + dr*(hfa2 + dr*hfa3)))*t
    desdr    = -hfgam*es - hfd*(hfa1 + dr*(2.0d0*hfa2 + 3.0d0*dr*hfa3))*t
    p1       = hfa1+2.0d0*hfa2*dr+3.0d0*hfa3*dr*dr
    p3       = 6.0d0*hfa3
    p2       = 2.0d0*hfa2 +p3*dr
    d1vhf    = desdr
    d2vhf    = -hfgam*d1vhf-hfd*t*(p2-hfgam*p1)
    d3vhf    = -hfgam*d2vhf-hfd*t*(p3-2.d0*hfgam*p2+hfgam*hfgam*p1)
    d4vhf    = -hfgam*d3vhf+hfd*t*(3.d0*hfgam*p3-3.d0*hfgam*hfgam*p2+ hfgam*hfgam*hfgam*p1)
    dcapk(1) = d1vhf
    dcapk(2) = d2vhf
    dcapk(3) = d3vhf
    dcapk(4) = d4vhf
    !
    !*****CALCULATE 1-4 DERIVATIVES OF UEHF
    !
    t      = rr*(h2tb1 + rr*(h2tb2 + rr*(h2tb3 + rr*h2tb4)))
    at     = h2ta*exp(-t)*rri
    t      = h2tb1 + rr*(2.0d0*h2tb2 + rr*(3.0d0*h2tb3 + rr*4.0d0*h2tb4))+rri
    dat    = -at*t
    q1     = -rri**2+2.d0*h2tb2+6.d0*h2tb3*rr+12.d0*h2tb4*rr**2
    q2     = 2.d0*rri**3+6.d0*h2tb3+24.d0*h2tb4*rr
    q3     = -6.d0*rri**4+24.d0*h2tb4
    d1uehf = dat
    d2uehf = -t*d1uehf-at*q1
    d3uehf = -t*d2uehf-2.d0*q1*d1uehf-at*q2
    d4uehf = -t*d3uehf-3.d0*q1*d2uehf-3.d0*q2*d1uehf-at*q3
    !
    !     Add in triplet and subtract singlet 2-body correlation terms
    !
    at = at + corrt - corrs
    !
    !*****INCLUDE CORRELATION PART OF CAPK,(USUBC-VSUBC)*(1-FZERO)
    !
    omfz     = 1.0d0-fzero
    dcapk(1) = dcapk(1)-(d1uehf+dcor2(1))*omfz
    dcapk(2) = dcapk(2)-(d2uehf+dcor2(2))*omfz
    dcapk(3) = dcapk(3)-(d3uehf+dcor2(3))*omfz
    dcapk(4) = dcapk(4)-(d4uehf+dcor2(4))*omfz
    !
    !*****CALCULATE 1-4 DERIVATIVES OF W-EFFECTIVE
    !
    t        = exp(-al3*rr)
    we       = (al0 + rr*(al1 + rr*al2))*t
    dwe      = (al1 + 2.0d0*al2*rr)*t - al3*we
    o2       = 2.0d0*al2
    o1       = al1+o2*rr
    d1w      = dwe
    d2w      = -al3*d1w+t*(o2-al3*o1)
    d3w      = -al3*d2w+t*(al3*o1-2.0d0*o2)*al3
    d4w      = -al3*d3w+t*(3.0d0*o2-al3*o1)*al3**2
    dcapk(1) = dcapk(1)-fzero*d1w
    dcapk(2) = dcapk(2)-fzero*d2w
    dcapk(3) = dcapk(3)-fzero*d3w
    dcapk(4) = dcapk(4)-fzero*d4w
    capk     = es-fzero*we-omfz*at
end subroutine kappa

!
!> @brief Evaluate POTINFO for the H3_GEN_DMBE_1987 potential.
subroutine potinfo
    implicit none
    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, deesdr, degsdr, deijdr, descart, descartnu, &
        & dgscart, dgscartnu, dijcart, dijcartnu, engyes, engygs, engyij, ezero, pengyes, pengygs, pengyij, r
    integer :: i, icartr, indexes, irctnt, ireorder, j, k, kediag, keoffd, kmax, ksdiag, ksoffd, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    character(len=75)           :: ref(5)
    integer, parameter          :: n3atom = 75
    integer, parameter          :: natom  = 25
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    common/utilcm/ dgscartnu(natom, 3), descartnu(natom, 3, isurf), &
        & dijcartnu(natom, 3, jsurf), cnvrtd, cnvrte, &
        & cnvrtde, ireorder, ksdiag, kediag, ksoffd, keoffd
    write(nflag(18), '(/)')
    do i =1, 5
        write(nflag(18), '(2x,a75)') ref(i)
    end do
    write(nflag(18), '(/)')
    kmax = 0
    do i = 1, isurf+1
        do j = 1, isurf+1
            if(nasurf(i, j).ne.0.and.kmax.lt.max(i, j)) kmax = max(i, j)
        end do
    end do
    write(nflag(18), '(2x,'' MAX. AND ACTUAL NO. OF EXCITED SURFACES: '',I3,5x,I3)') msurf, kmax-1
    if(kmax-1.gt.msurf) then
        write(6, *) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
        stop
    end if
    ksdiag = 0
    kediag = 0
    do i = 2, isurf+1
        if(nasurf(i, i).ne.0) then
            kediag = i-1
            if(ksdiag.eq.0) ksdiag = i-1
        end if
    end do
    ksoffd = 0
    keoffd = 0
    k      = 0
    do i = 1, isurf
        do j = i+1, isurf+1
            k = k+1
            if(nasurf(i, j)+nasurf(j, i).ne.0) then
                keoffd = k
                if(ksoffd.eq.0) ksoffd = k
            end if
        end do
    end do
    write(nflag(18), '(2x,'' MAX. AND ACTUAL ORDER OF DERIVATIVES:    '',I3,5x,I3)') mder, nder
    if(nder.gt.mder) then
        write(6, *) ' WRONG INPUT ON ORDER OF DERIVATIVES'
        stop
    end if
    if(nflag(19).eq.1) then
        write(nflag(18), '(/)')
        write(nflag(18), '(2x,''Cartesian coordinates are supplied by'',/, 2x,''the user in the array CART.'',//)')
        write(nflag(18), '(2x,''Provide cartesian coordinates in the'',/, 2x,''following order using the array CART'',' // &
            '//, 2x,'' CART(1,1)...CART(1,3)   => ATOM 1'',/, 2x,'' CART(2,1)...CART(2,3)   => AT' // &
            'OM 2'',/, 2x,'' CART(3,1)...CART(3,3)   => ATOM 3'',/, 2x,'' CART(N,1)...CART(N,3)  ' // &
            ' => ATOM N'',/, 2x,''CART(25,1)...CART(25,3)  => ATOM 25'',/)')
        write(nflag(18), '(2x,''If the user wishes to relabel the atoms,'',/, 2x,''set the variable IREORDER equal to ' // &
            '1'',/, 2x,''in the PARAMETER statement.  The user'',/, 2x,''must also specify the ne' // &
            'w labeling'',/, 2x,''scheme.  This is done using the array'',/, 2x,''NULBL in the fo' // &
            'llowing manner:'',//, 2x,''NULBL(i) = j'',/, 2x,''where:  i => old index'',/, 2x,'' ' // &
            '       j => new index'',//)')
        write(nflag(18), '(2x,''Cartesian coordinates can be provided to'',/, 2x,''the potential routine in a variety ' // &
            'of units.'',/, 2x,''The input units will be converted to Bohr'',/, 2x,''based on the' // &
            ' following values of the NFLAG'',/, 2x,''variable:'',//, 2x,''NFLAG(1)  =  1  =>  CA' // &
            'RTESIANS IN BOHR (no'',/, 2x,''                    conversion required)'',/, 2x,''NF' // &
            'LAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS'',//)')
        write(nflag(18), '(2x,''The value of the energy and derivatives'',/, 2x,''(if computed) can be reported in a v' // &
            'ariety'',/, 2x,''units.  A units conversion will take place'',/, 2x,''as specified b' // &
            'y the following values of the'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(2) = 1 =>  ' // &
            'ENERGIES REPORTED IN HARTEEE'',/, 2x,''NFLAG(2) = 2 =>  ENERGIES REPORTED IN mHARTRE' // &
            'E'',/, 2x,''NFLAG(2) = 3 =>  ENERGIES REPORTED IN eV'',/, 2x,''NFLAG(2) = 4 =>  ENER' // &
            'GIES REPORTED IN kcal/mol'',/, 2x,''NFLAG(2) = 5 =>  ENERGIES REPORTED IN cm**-1'',/' // &
            '/)')
        write(nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
            'he'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(1)=1 & NFLAG(2)=1 => DERIVATIVES REPOR' // &
            'TED IN'',/, 2x,''                           HARTEEE/BOHR'',/, 2x,''NFLAG(1)=1 & NFLA' // &
            'G(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/BOHR' // &
            ''',/, 2x,''NFLAG(1)=1 & NFLAG(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''            ' // &
            '               eV/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=4 => DERIVATIVES REPORTED IN' // &
            ''',/, 2x,''                           kcal/mol/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=' // &
            '5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm**-1/BOHR'',//)')
        write(nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
            'he'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(1)=2 & NFLAG(2)=1 => DERIVATIVES REPOR' // &
            'TED IN'',/, 2x,''                           HARTEEE/ANGSTROM'',/, 2x,''NFLAG(1)=2 & ' // &
            'NFLAG(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/' // &
            'ANGSTROM'',/, 2x,''NFLAG(1)=2 & NFLAG(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''    ' // &
            '                       eV/ANGSTROM'',/, 2x,''NFLAG(1)=2 & NFLAG(2)=4 => DERIVATIVES ' // &
            'REPORTED IN'',/, 2x,''                           kcal/mol/ANGSTROM'',/, 2x,''NFLAG(1' // &
            ')=2 & NFLAG(2)=5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm*' // &
            '*-1/ANGSTROM'',//)')
    end if
end subroutine potinfo

!> @brief Evaluate ANCVRT for the H3_GEN_DMBE_1987 potential.
subroutine ancvrt
    implicit none
    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, deesdr, degsdr, deijdr, descart, descartnu, &
        & dgscart, dgscartnu, dijcart, dijcartnu, engyes, engygs, engyij, ezero, pengyes, pengygs, pengyij, r
    integer :: i, ianum, icartr, inc1, inc2, ind, ind2, indexes, inu, irctnt, ireorder, isave, isum, j, jsave, kediag, &
        & keoffd, ksdiag, ksoffd, mder, msurf, nasurf, natoms, nder, nflag, nulbl
    character(len=75)  :: ref(5)
    character(len=3)   :: periodic_1(7,32)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    character(len=2)   :: name1(natom)
    character(len=2)   :: name2(natom)
    character(len=1)   :: iblank
    character(len=20)  :: distance
    character(len=20)  :: units
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    common/utilcm/ dgscartnu(natom, 3), descartnu(natom, 3, isurf), &
        & dijcartnu(natom, 3, jsurf), cnvrtd, cnvrte, &
        & cnvrtde, ireorder, ksdiag, kediag, ksoffd, keoffd
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    dimension ianum(7,32)
    dimension isave(natom), jsave(natom)
    double precision, parameter :: pi         = 3.141592653589793d0
    double precision, parameter :: clight     = 2.99792458d08
    double precision, parameter :: cmu_0      = 4.0d0*pi*1.0d-07
    double precision, parameter :: cepsilon_0 = 1.0d0/(cmu_0*clight**2)
    double precision, parameter :: ce         = 1.602176462d-19
    double precision, parameter :: cplanck    = 6.62606876d-34
    double precision, parameter :: cm_e       = 9.10938188d-31
    double precision, parameter :: cang       = 1.0d-10
    double precision, parameter :: cavogadro  = 6.02214199d23
    double precision, parameter :: ckcal      = 4.184d10
    double precision, parameter :: htomillh   = 1000.d0
    double precision, parameter :: htoev      = 27.2113834d0
    double precision, parameter :: htokcal    = 627.509470d0
    double precision, parameter :: htowave    = 219474.631d0
    double precision, parameter :: htokj      = 2625.49962d0
    double precision, parameter :: bohr_a     = .5291772083d0
    do i=1, 7
        do j=1, 32
            ianum(i,j)      = 0
            periodic_1(i,j) = ' '
        end do
    end do
    distance         = 'BOHR                '
    units            = 'HARTREE             '
    ianum(1,1)       = 1
    ianum(1,32)      = 2
    ianum(2,1)       = 3
    ianum(2,2)       = 4
    ianum(2,27)      = 5
    ianum(2,28)      = 6
    ianum(2,29)      = 7
    ianum(2,30)      = 8
    ianum(2,31)      = 9
    ianum(2,32)      = 10
    ianum(3,1)       = 11
    ianum(3,2)       = 12
    ianum(3,27)      = 13
    ianum(3,28)      = 14
    ianum(3,29)      = 15
    ianum(3,30)      = 16
    ianum(3,31)      = 17
    ianum(3,32)      = 18
    ianum(4,1)       = 19
    ianum(4,2)       = 20
    ianum(4,17)      = 21
    ianum(4,18)      = 22
    ianum(4,19)      = 23
    ianum(4,20)      = 24
    ianum(4,21)      = 25
    ianum(4,22)      = 26
    ianum(4,23)      = 27
    ianum(4,24)      = 28
    ianum(4,25)      = 29
    ianum(4,26)      = 30
    ianum(4,27)      = 31
    ianum(4,28)      = 32
    ianum(4,29)      = 33
    ianum(4,30)      = 34
    ianum(4,31)      = 35
    ianum(4,32)      = 36
    ianum(5,1)       = 37
    ianum(5,2)       = 38
    ianum(5,17)      = 39
    ianum(5,18)      = 40
    ianum(5,19)      = 41
    ianum(5,20)      = 42
    ianum(5,21)      = 43
    ianum(5,22)      = 44
    ianum(5,23)      = 45
    ianum(5,24)      = 46
    ianum(5,25)      = 47
    ianum(5,26)      = 48
    ianum(5,27)      = 49
    ianum(5,28)      = 50
    ianum(5,29)      = 51
    ianum(5,30)      = 52
    ianum(5,31)      = 53
    ianum(5,32)      = 54
    ianum(6,1)       = 55
    ianum(6,2)       = 56
    ianum(6,3)       = 57
    ianum(6,4)       = 58
    ianum(6,5)       = 59
    ianum(6,6)       = 60
    ianum(6,7)       = 61
    ianum(6,8)       = 62
    ianum(6,9)       = 63
    ianum(6,10)      = 64
    ianum(6,11)      = 65
    ianum(6,12)      = 66
    ianum(6,13)      = 67
    ianum(6,14)      = 68
    ianum(6,15)      = 69
    ianum(6,16)      = 70
    ianum(6,17)      = 71
    ianum(6,18)      = 72
    ianum(6,19)      = 73
    ianum(6,20)      = 74
    ianum(6,21)      = 75
    ianum(6,22)      = 76
    ianum(6,23)      = 77
    ianum(6,24)      = 78
    ianum(6,25)      = 79
    ianum(6,26)      = 80
    ianum(6,27)      = 81
    ianum(6,28)      = 82
    ianum(6,29)      = 83
    ianum(6,30)      = 84
    ianum(6,31)      = 85
    ianum(6,32)      = 86
    ianum(7,1)       = 87
    ianum(7,2)       = 88
    ianum(7,3)       = 89
    ianum(7,4)       = 90
    ianum(7,5)       = 91
    ianum(7,6)       = 92
    ianum(7,7)       = 93
    ianum(7,8)       = 94
    ianum(7,9)       = 95
    ianum(7,10)      = 96
    ianum(7,11)      = 97
    ianum(7,12)      = 98
    ianum(7,13)      = 99
    ianum(7,14)      = 100
    ianum(7,15)      = 101
    ianum(7,16)      = 102
    ianum(7,17)      = 103
    ianum(7,18)      = 104
    ianum(7,19)      = 105
    ianum(7,20)      = 106
    ianum(7,21)      = 107
    ianum(7,22)      = 108
    ianum(7,23)      = 109
    ianum(7,24)      = 110
    ianum(7,25)      = 111
    ianum(7,26)      = 112
    ianum(7,27)      = 113
    ianum(7,28)      = 114
    ianum(7,29)      = 115
    ianum(7,30)      = 116
    ianum(7,31)      = 117
    ianum(7,32)      = 120
    periodic_1(1,1)  = 'H  '
    periodic_1(1,32) = 'He '
    periodic_1(2,1)  = 'Li '
    periodic_1(2,2)  = 'Be '
    periodic_1(2,27) = 'B  '
    periodic_1(2,28) = 'C  '
    periodic_1(2,29) = 'N  '
    periodic_1(2,30) = 'O  '
    periodic_1(2,31) = 'F  '
    periodic_1(2,32) = 'Ne '
    periodic_1(3,1)  = 'Na '
    periodic_1(3,2)  = 'Mg '
    periodic_1(3,27) = 'Al '
    periodic_1(3,28) = 'Si '
    periodic_1(3,29) = 'P  '
    periodic_1(3,30) = 'S  '
    periodic_1(3,31) = 'Cl '
    periodic_1(3,32) = 'Ar '
    periodic_1(4,1)  = 'K  '
    periodic_1(4,2)  = 'Ca '
    periodic_1(4,17) = 'Sc '
    periodic_1(4,18) = 'Ti '
    periodic_1(4,19) = 'V  '
    periodic_1(4,20) = 'Cr '
    periodic_1(4,21) = 'Mn '
    periodic_1(4,22) = 'Fe '
    periodic_1(4,23) = 'Co '
    periodic_1(4,24) = 'Ni '
    periodic_1(4,25) = 'Cu '
    periodic_1(4,26) = 'Zn '
    periodic_1(4,27) = 'Ga '
    periodic_1(4,28) = 'Ge '
    periodic_1(4,29) = 'As '
    periodic_1(4,30) = 'Se '
    periodic_1(4,31) = 'Br '
    periodic_1(4,32) = 'Kr '
    periodic_1(5,1)  = 'Rb '
    periodic_1(5,2)  = 'Sr '
    periodic_1(5,17) = 'Y  '
    periodic_1(5,18) = 'Zr '
    periodic_1(5,19) = 'Nb '
    periodic_1(5,20) = 'Mo '
    periodic_1(5,21) = 'Tc '
    periodic_1(5,22) = 'Ru '
    periodic_1(5,23) = 'Rh '
    periodic_1(5,24) = 'Pd '
    periodic_1(5,25) = 'Ag '
    periodic_1(5,26) = 'Cd '
    periodic_1(5,27) = 'In '
    periodic_1(5,28) = 'Sn '
    periodic_1(5,29) = 'Sb '
    periodic_1(5,30) = 'Te '
    periodic_1(5,31) = 'I  '
    periodic_1(5,32) = 'Xe '
    periodic_1(5,32) = 'Xe '
    do i=1, natoms
        isave(i) = 0
        jsave(i) = 0
        name1(i) = '  '
        name2(i) = '  '
    end do
    iblank = ' '
    do ind=1, natoms
        do i=1, 7
            do j=1, 32
                if(indexes(ind).eq.ianum(i,j)) then
                    isave(ind) = i
                    jsave(ind) = j
                end if
            end do
        end do
    end do

    do ind=1, natoms
        ind2 = nulbl(ind)
        if(ind2.eq.0) ind2=ind
    end do
    inc1 = 0
    do ind=1, irctnt-1
        inc1        = inc1+1
        name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
    end do
    inc2 = 0
    do ind=irctnt, natoms
        inc2        = inc2+1
        name2(inc2) = periodic_1(isave(ind),jsave(ind))(:2)
    end do
    if(nflag(1).eq.2) distance = 'ANGSTROMS           '
    if(nflag(2).eq.2) then
        units = 'MILLIHARTREE        '
    elseif(nflag(2).eq.3) then
        units = 'EV                  '
    elseif(nflag(2).eq.4) then
        units = 'KCAL PER MOLE       '
    elseif(nflag(2).eq.5) then
        units = 'WAVENUMBERS         '
    elseif(nflag(2).eq.6) then
        units = 'KILOJOULES PER MOLE '
    end if
    cnvrtd  = 1.d0
    cnvrte  = 1.d0
    cnvrtde = 1.d0
    if(nflag(1).eq.2) cnvrtd = bohr_a
    if(nflag(2).eq.2) then
        cnvrte = cnvrte*htomillh
    elseif(nflag(2).eq.3) then
        cnvrte = cnvrte*htoev
    elseif(nflag(2).eq.4) then
        cnvrte = cnvrte*htokcal
    elseif(nflag(2).eq.5) then
        cnvrte = cnvrte*htowave
    elseif(nflag(2).eq.6) then
        cnvrte = cnvrte*htokj
    end if
    cnvrtde = cnvrte/cnvrtd
    isum    = 0
    do inu=1, 25
        isum = isum + nulbl(inu)
    end do
    ireorder = 0
    if(isum.ne.0) ireorder = 1
end subroutine ancvrt

!> @brief Evaluate CARTOU for the H3_GEN_DMBE_1987 potential.
subroutine cartou
    use chempotpy_legacy_coordinates, only: convert_cartesian_input
    implicit none
    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, descartnu, dgscartnu, &
        & dijcartnu
    integer :: i, icartr, indexes, irctnt, ireorder, j, kediag, keoffd, ksdiag, ksoffd, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    common/utilcm/ dgscartnu(natom, 3), descartnu(natom, 3, isurf), &
        & dijcartnu(natom, 3, jsurf), cnvrtd, cnvrte, &
        & cnvrtde, ireorder, ksdiag, kediag, ksoffd, keoffd
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    call convert_cartesian_input( cart, nulbl, natoms, cnvrtd, ireorder == 1, cartnu )
end subroutine cartou

!> @brief Evaluate CARTTOR for the H3_GEN_DMBE_1987 potential.
subroutine carttor
    use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
    implicit none
    integer          :: coordinate_status
    double precision :: anuzero, cart, cartnu, cmm, cphi, degsdr, engygs, flm, hhd, hym, phi, pi, q, &
        & q1, q2, q3, r, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, &
        & zcm1, zcm2, zcm3, zrm1, zrm2
    integer :: i, icartr, ind, indexes, irctnt, k, l, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    call cartesian_to_potlib_coordinates( cartnu, natoms, icartr, r, coordinate_status )
    if (coordinate_status /= 0) then
        write(nflag(18), '(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') icartr
        stop
    end if
end subroutine carttor

!> @brief Evaluate EUNITZERO for the H3_GEN_DMBE_1987 potential.
subroutine eunitzero
    implicit none
    double precision :: anuzero, cart, cnvrtd, cnvrtde, cnvrte, deesdr, degsdr, deijdr, descart, &
        & descartnu, dgscart, dgscartnu, dijcart, dijcartnu, engyes, engygs, engyij, ezero, pengyes, &
        & pengygs, pengyij, r
    integer :: i, ireorder, j, kediag, keoffd, ksdiag, ksoffd, nasurf, &
        & nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    common/utilcm/ dgscartnu(natom, 3), descartnu(natom, 3, isurf), &
        & dijcartnu(natom, 3, jsurf), cnvrtd, cnvrte, &
        & cnvrtde, ireorder, ksdiag, kediag, ksoffd, keoffd
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    pengygs = engygs * cnvrte - anuzero
    if(ksdiag.ne.0) then
        do i=ksdiag, kediag
            pengyes(i) = engyes(i) * cnvrte - anuzero
        end do
    end if
    if(ksoffd.ne.0) then
        do j=ksoffd, keoffd
            pengyij(j) = engyij(j) * cnvrte
        end do
    end if
end subroutine eunitzero

!> @brief Evaluate RTOCART for the H3_GEN_DMBE_1987 potential.
subroutine rtocart
    implicit none
    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, deesdr, degsdr, deijdr, &
        & descartnu, dgscartnu, dijcartnu, engyes, engygs, engyij, eps, epsp, ezero, phi, r, r1, r2, &
        & rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, termy, termz, theta1, &
        & theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes, ygs, yij, za, zb, zc, &
        & zd
    integer :: i, icartr, ind, indexes, irctnt, ireorder, j, j1, j2, k, kediag, keoffd, ksdiag, ksoffd, m1, mder, msurf, &
        & nasurf, natoms, nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/  engyij(jsurf), deijdr(n3atom, jsurf)
    common /utilcm/ dgscartnu(natom, 3), descartnu(natom, 3, isurf), &
        & dijcartnu(natom, 3, jsurf), cnvrtd, cnvrte, cnvrtde, &
        & ireorder, ksdiag, kediag, ksoffd, keoffd
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    dimension ygs(n3atom), yes(n3atom,isurf), yij(n3atom,jsurf)
    if(icartr.eq.1) then
        do i = 1, natoms
            ind             = 3*i-2
            dgscartnu(i, 1) = degsdr(ind)
            dgscartnu(i, 2) = degsdr(ind+1)
            dgscartnu(i, 3) = degsdr(ind+2)
            if(ksdiag.ne.0) then
                do j = ksdiag, kediag
                    descartnu(i, 1, j) = deesdr(ind, j)
                    descartnu(i, 2, j) = deesdr(ind+1, j)
                    descartnu(i, 3, j) = deesdr(ind+2, j)
                end do
            end if
            if(keoffd.ne.0) then
                do k = ksoffd, keoffd
                    dijcartnu(i, 1, k) = deijdr(ind, k)
                    dijcartnu(i, 2, k) = deijdr(ind+1, k)
                    dijcartnu(i, 3, k) = deijdr(ind+2, k)
                end do
            end if
        end do
    elseif(icartr.eq.2) then
        do i = 1, natoms
            dgscartnu(i, 1) = 0.d0
            dgscartnu(i, 2) = 0.d0
            dgscartnu(i, 3) = 0.d0
            if(ksdiag.ne.0) then
                do j1=ksdiag, kediag
                    descartnu(i, 1, j1) = 0.d0
                    descartnu(i, 2, j1) = 0.d0
                    descartnu(i, 3, j1) = 0.d0
                end do
            end if
            if(ksoffd.ne.0) then
                do j2=ksoffd, keoffd
                    dijcartnu(i, 1, j2) = 0.d0
                    dijcartnu(i, 2, j2) = 0.d0
                    dijcartnu(i, 3, j2) = 0.d0
                end do
            end if
            do j = 1, natoms
                if(j.lt.i) then
                    m1 = natoms*(j-1) - (j*(j-1))/2 + i-j
                elseif(j.gt.i) then
                    m1 = natoms*(i-1) - (i*(i-1))/2 + j-i
                else
                    go to 20
                end if
                y               = degsdr(m1)
                termx           = (cartnu(i, 1)-cartnu(j, 1))/r(m1)
                termy           = (cartnu(i, 2)-cartnu(j, 2))/r(m1)
                termz           = (cartnu(i, 3)-cartnu(j, 3))/r(m1)
                dgscartnu(i, 1) = dgscartnu(i, 1) + termx*y
                dgscartnu(i, 2) = dgscartnu(i, 2) + termy*y
                dgscartnu(i, 3) = dgscartnu(i, 3) + termz*y
                if(ksdiag.gt.0) then
                    y = deesdr(m1, j1)
                    do j1=ksdiag, kediag
                        descartnu(i, 1, j1) = descartnu(i, 1, j1) + termx*y
                        descartnu(i, 2, j1) = descartnu(i, 2, j1) + termy*y
                        descartnu(i, 3, j1) = descartnu(i, 3, j1) + termz*y
                    end do
                elseif(ksoffd.gt.0) then
                    do j2=ksoffd, keoffd
                        y                   = deijdr(m1, j2)
                        dijcartnu(i, 1, j2) = dijcartnu(i, 1, j2) + termx*y
                        dijcartnu(i, 2, j2) = dijcartnu(i, 2, j2) + termy*y
                        dijcartnu(i, 3, j2) = dijcartnu(i, 3, j2) + termz*y
                    end do
                end if
                20 continue
            end do
        end do
    elseif(icartr.eq.3) then
        do i = 1, natoms
            ygs(i) = degsdr(i)/r(i)
            if(ksdiag.ne.0) then
                do j=ksdiag, kediag
                    yes(i,j) = deesdr(i, j)/r(i)
                end do
            end if
            if(ksoffd.ne.0) then
                do k=ksoffd, keoffd
                    yij(i,k) = deijdr(i, k)/r(i)
                end do
            end if
        end do
        do k = 1, 3
            term12          = cartnu(1, k)-cartnu(2, k)
            term23          = cartnu(2, k)-cartnu(3, k)
            term13          = cartnu(1, k)-cartnu(3, k)
            dgscartnu(1, k) = term12*ygs(1) + term13*ygs(3)
            dgscartnu(2, k) = -term12*ygs(1) + term23*ygs(2)
            dgscartnu(3, k) = -term13*ygs(3) - term23*ygs(2)
            if(ksdiag.ne.0) then
                do j1=ksdiag, kediag
                    descartnu(1, k, j1) = term12*yes(1,j1) + term13*yes(3,j1)
                    descartnu(2, k, j1) = -term12*yes(1,j1) + term23*yes(2,j1)
                    descartnu(3, k, j1) = -term13*yes(3,j1) - term23*yes(2,j1)
                end do
            end if
            if(ksoffd.ne.0) then
                do j2=ksoffd, keoffd
                    dijcartnu(1, k, j2) = term12*yij(1,j2) + term13*yij(3,j2)
                    dijcartnu(2, k, j2) = -term12*yij(1,j2) + term23*yij(2,j2)
                    dijcartnu(3, k, j2) = -term13*yij(3,j2) - term23*yij(2,j2)
                end do
            end if
        end do
    elseif(icartr.eq.4) then
        wh   = 1.007825d0
        wf   = 18.99840d0
        sum  = wh+wf
        eps  = wf/sum
        epsp = wh/sum
        u1   = cos(theta1)
        u2   = cos(theta2)
        u3   = cos(phi)
        ss1  = sin(theta1)
        ss2  = sin(theta2)
        ss3  = sin(phi)
        ya   = 0.0d0
        yb   = 0.0d0
        t0   = r1*u1
        za   = -epsp*t0
        zb   = eps*t0
        t0   = r1*ss1
        xa   = -epsp*t0
        xbb  = eps*t0
        t0   = r2*ss2
        t1   = t0*u3
        xc   = -epsp*t1
        xd   = eps*t1
        t1   = t0*ss3
        yc   = -epsp*t1
        yd   = eps*t1
        t0   = r2*u2
        zc   = -epsp*t0+rcm
        zd   = eps*t0+rcm
        rff  = sqrt((xa-xc)**2+yc**2+(za-zc)**2)
    else
        write(nflag(18), '(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') icartr
        stop
    end if
end subroutine rtocart

!> @brief Evaluate DEDCOU for the H3_GEN_DMBE_1987 potential.
subroutine dedcou
    implicit none
    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, descart, descartnu, dgscart, &
        & dgscartnu, dijcart, dijcartnu, pengyes, pengygs, pengyij
    integer :: i, icartr, indexes, irctnt, ireorder, j, k, kediag, keoffd, ksdiag, ksoffd, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    common /utilcm/ dgscartnu(natom, 3), descartnu(natom, 3, isurf), &
        & dijcartnu(natom, 3, jsurf), cnvrtd, cnvrte, cnvrtde, &
        & ireorder, ksdiag, kediag, ksoffd, keoffd
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    if (ireorder.eq.1) then
        do i = 1, natoms
            dgscart(i, 1) = dgscartnu(nulbl(i), 1) * cnvrtde
            dgscart(i, 2) = dgscartnu(nulbl(i), 2) * cnvrtde
            dgscart(i, 3) = dgscartnu(nulbl(i), 3) * cnvrtde
            if(ksdiag.ne.0) then
                do j=ksdiag, kediag
                    descart(i, 1, j) = descartnu(nulbl(i), 1, j) * cnvrtde
                    descart(i, 2, j) = descartnu(nulbl(i), 2, j) * cnvrtde
                    descart(i, 3, j) = descartnu(nulbl(i), 3, j) * cnvrtde
                end do
            end if
            if(ksoffd.ne.0) then
                do k=ksoffd, keoffd
                    dijcart(i, 1, k) = dijcartnu(nulbl(i), 1, k) * cnvrtde
                    dijcart(i, 2, k) = dijcartnu(nulbl(i), 2, k) * cnvrtde
                    dijcart(i, 3, k) = dijcartnu(nulbl(i), 3, k) * cnvrtde
                end do
            end if
        end do
    else
        do i = 1, natoms
            dgscart(i, 1) = dgscartnu(i, 1) * cnvrtde
            dgscart(i, 2) = dgscartnu(i, 2) * cnvrtde
            dgscart(i, 3) = dgscartnu(i, 3) * cnvrtde
            if(ksdiag.ne.0) then
                do j=ksdiag, kediag
                    descart(i, 1, j) = descartnu(i, 1, j) * cnvrtde
                    descart(i, 2, j) = descartnu(i, 2, j) * cnvrtde
                    descart(i, 3, j) = descartnu(i, 3, j) * cnvrtde
                end do
            end if
            if(ksoffd.ne.0) then
                do k=ksoffd, keoffd
                    dijcart(i, 1, k) = dijcartnu(i, 1, k) * cnvrtde
                    dijcart(i, 2, k) = dijcartnu(i, 2, k) * cnvrtde
                    dijcart(i, 3, k) = dijcartnu(i, 3, k) * cnvrtde
                end do
            end if
        end do
    end if
end subroutine dedcou
