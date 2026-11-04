!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the BrH2_LEPS_LTBZ1_1995 potential.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
    implicit none
    ! number of electronic state
    integer, parameter            :: nstates = 1
    integer, parameter            :: natoms  = 3
    integer, parameter            :: natom   = 25
    integer, parameter            :: isurf   = 5
    integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    double precision              :: anuzero, cart, descart, dgscart, dijcart, pengyes, pengygs, pengyij, v
    integer                       :: iatom, idir, istate, nasurf, nder, nflag, nulbl
    common/usrocm/ pengygs, pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom, 3), descart(natom, 3, isurf), &
        & dijcart(natom, 3, jsurf)
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    logical, save :: first_time_data=.true.
    !initialize
    v    = 0.d0; g = 0.d0; d = 0.d0
    cart = 0.d0
    do iatom=1, natoms
        do idir=1, 3
            cart(iatom, idir) = x(iatom,idir)/legacy_bohr_to_angstrom
        enddo
    enddo
    nder = igrad
    if( first_time_data )then
        call prepot
        first_time_data = .false.
    endif
    call pot
    if( igrad == 0 )then
        do istate=1, nstates
            p(istate) = pengygs*legacy_hartree_to_electronvolt
        enddo
    else if( igrad == 1 )then
        do istate=1, nstates
            p(istate) = pengygs*legacy_hartree_to_electronvolt
        enddo
        do iatom=1, natoms
            do idir=1, 3
                g(1,iatom,idir) = dgscart(iatom, idir)*51.422067
            enddo
        enddo
    else if( igrad == 2 )then
        write (*, *) 'Only energy and gradient are available'
    endif
end subroutine pes

!
!> @brief Evaluate PREPOT for the BrH2_LEPS_LTBZ1_1995 potential.
subroutine prepot
    !
    !   System:          BrH2
    !   Functional form: extended-LEPS plus three-center term
    !   Common name:     BrH2SEC
    !   Interface:       potlib2001
    !   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g.
    !                    J. Phys. Chem. 99, 207 (1995)
    !
    !   Number of bodies:              3
    !   Number of derivatives:         1
    !   Number of electronic surfaces: 1
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
    !      coordinates - bohr
    !      derivatives - hartrees/bohr
    !
    !   Surfaces:
    !      ground electronic state
    !
    !   Zero of energy:
    !      The classical potential energy is set equal to zero for the Br
    !      infinitely far from the H2 diatomic and R(H2) set equal to the
    !      H2 equilibrium diatomic value.
    !
    !   Parameters:
    !      Set in the BLOCK DATA subprogram PTPACM
    !
    !   Coordinates:
    !      Internal, Definition: R(1) = R(Br-H)
    !                            R(2) = R(H-H)
    !                            R(3) = R(Br-H)
    !
    !   Common Blocks (used between the calling program and this potential):
    !        passes the coordinates, ground state electronic energy, and
    !        derivatives of the ground electronic state energy with respect
    !        to the coordinates.
    !        passes the control flags where
    !          NSURF  - not used
    !        NDER  = 0 => no derivatives are computed
    !              = 1 => derivatives of the energy for the ground electroni
    !                     state with respect to the coordinates are computed
    !          NFLAG  - Control flags
    !      NFLAG(18-20)
    !        passes the FORTRAN unit number used for potential output
    !      /ASYCM/ EASYAB, EASYBC, EASYAC
    !        passes the energy in the three asymptotic valleys for an A + BC
    !        The energy in the AB valley, EASYAB, is equal to the energy of
    !        C atom "infinitely" far from the AB diatomic and R(AB) set equa
    !        Re(AB), the equilibrium bond length for the AB diatomic.
    !        In this potential the AB valley represents H infinitely far fro
    !        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the t
    !        EASYBC and EASYAC represent the energies in the H2 and the BrH
    !        valleys, respectively.
    !
    !   Default Parameter Values:
    !      Variable      Default value
    !      NDER             1
    !      NFLAG(18)        6
    !
    !*****
    !
    implicit none
    double precision :: anuzero, beta, cart, cartnu, coul, de, deq, derad, do4z, dv3c, dxdr, dz, easyab, easyac, easybc, &
        & exch, ezero, op3z, rad, re, top3z, tzp3, v3c, v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, x, z, zhbrp1, zhbrp2, &
        & zhhp1, zhhp2, zp3, zpo
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
    common /pt3cm/ ezero(isurf+1)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    common /asycm/ easyab, easybc, easyac
    !
    double precision, parameter :: cangau = 0.52917706d0
    double precision, parameter :: charkc = 627.5095d0
    double precision, parameter :: r2     = 1.41421356d0
    !
    common /satocm/ zhhp1, zhhp2, zhbrp1, zhbrp2, &
        & z(3), dz(3, 3)
    common /v3ccm/ v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, &
        & v3c, dv3c(3)
    common /lengcm/ x(3), coul(3), exch(3), rad
    common /lsatcm/ zpo(3), op3z(3), zp3(3), tzp3(3), top3z(3), &
        & do4z(3), beta(3)
    common /ldercm/ dxdr(3, 3), deq(3), derad(3)
    common /parmcm/ de(3), re(3)
    !
    save
    !
    if( natoms .gt. 25 )then
        write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    endif
    !
    !
    !   Convert the constants to atomic units
    !
    do 10 i = 1, 3
        de(i)   = de(i) / charkc
        re(i)   = re(i) / cangau
        beta(i) = beta(i) * cangau
    10 continue
    !
    !    Initialize the energy in the three asymptotic valleys
    !
    easyab = de(1)
    easybc = de(2)
    easyac = de(3)
    !
    !
    ezero(1) = de(2)
    !
    do i=1, 5
        ref(i) = ' '
    enddo
    !
    ref(1) = 'G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g. Zhao'
    ref(2) = 'J. Phys. Chem. 99, 207 (1995)'
    !
    indexes(1) = 35
    indexes(2) = 1
    indexes(3) = 1
    !
    !
    !
    irctnt = 2
    !
    !CALL POTINFO
    !
    call ancvrt
    !
end subroutine prepot

!
!> @brief Evaluate POT for the BrH2_LEPS_LTBZ1_1995 potential.
subroutine pot
    !
    !   The potential parameters must be passed through the COMMON BLOCKS
    !   PT1CM, ASYCM, PT2CM,  All information passed through the
    !   COMMON BLOCKS PT1CM, ASYCM must be in Hartree atomic units.
    !
    !   Coordinates:
    !      Internal, Definition: R(1) = R(Br-H)
    !                            R(2) = R(H-H)
    !                            R(3) = R(Br-H)
    !
    !        The energy in the AB valley, EASYAB, is equal to the energy of
    !        C atom "infinitely" far from the AB diatomic and R(AB) set equa
    !        Re(AB), the equilibrium bond length for the AB diatomic.
    !
    !        In this potential the AB valley represents H infinitely far fro
    !        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the t
    !        EASYBC and EASYAC represent the energies in the H2 and the HBr
    !        valleys, respectively.
    !
    !      ENTRY POT
    implicit none
    double precision :: anuzero, beta, cart, cartnu, coul, de, deesdr, degsdr, deijdr, deq, derad, descart, dgscart, &
        & dijcart, do4z, dv3c, dxdr, dz, easyab, easyac, easybc, engyes, engygs, engyij, exch, ezero, op3z, pengyes, &
        & pengygs, pengyij, r, rad, re, top3z, tzp3, v3c, v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, x, z, zhbrp1, zhbrp2, &
        & zhhp1, zhhp2, zp3, zpo
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
    common /pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/ ezero(isurf+1)
    common /pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /asycm/ easyab, easybc, easyac
    !
    double precision, parameter :: cangau = 0.52917706d0
    double precision, parameter :: charkc = 627.5095d0
    double precision, parameter :: r2     = 1.41421356d0
    !
    common /satocm/ zhhp1, zhhp2, zhbrp1, zhbrp2, &
        & z(3), dz(3, 3)
    common /v3ccm/ v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, &
        & v3c, dv3c(3)
    common /lengcm/ x(3), coul(3), exch(3), rad
    common /lsatcm/ zpo(3), op3z(3), zp3(3), tzp3(3), top3z(3), &
        & do4z(3), beta(3)
    common /ldercm/ dxdr(3, 3), deq(3), derad(3)
    common /parmcm/ de(3), re(3)
    !
    call cartou
    call carttor
    !
    !   Check the value of NDER
    !
    !   Compute bond dependent Sato parameters
    !
    call sato
    !
    !   Compute the LEPS constants
    !
    do 20 i = 1, 3
        zpo(i)   = 1.0d0+z(i)
        op3z(i)  = 1.0d0+3.0d0*z(i)
        top3z(i) = 2.0d0*op3z(i)
        zp3(i)   = z(i)+3.0d0
        tzp3(i)  = 2.0d0*zp3(i)
        do4z(i)  = de(i)/4.0d0/zpo(i)
    20 continue
    !
    engygs = 0.d0
    !
    do 30 i = 1, 3
        x(i)    = exp(-beta(i)*(r(i)-re(i)))
        coul(i) = do4z(i) * (zp3(i)*x(i) - top3z(i)) * x(i)
        exch(i) = do4z(i) * (op3z(i)*x(i) - tzp3(i)) * x(i)
        engygs  = engygs + coul(i)
    30 continue
    !
    rad = sqrt((exch(1) - exch(2)) ** 2 + (exch(2) - exch(3)) ** 2 + (exch(3) - exch(1)) ** 2 )
    !
    engygs = engygs - (rad / r2) + ezero(1)
    !
    !   Compute the three-center term
    call v3cfcn
    !
    !   Sum all the terms that make up the energy
    !
    engygs = engygs + v3c
    !
    !   Compute the derivatives if NDER = 1
    !
    if( nder .ne. 1 ) go to 700
    !
    !   DCOUL computes the derivatives of the coulombic terms w.r.t. R(I)
    !   DEXCH computes the derivatives of the exchange terms w.r.t. R(I) and
    !         also computes the derivative of the term RAD w.r.t. R(i)
    !
    !   DCOUL must be called before DEXCH becaues this subprogram computes
    !   the DXDR terms needed by both subprograms.
    !
    !   Compute DEQ = dE/dQ * [dQ/dX * dX/dR + dQ/dZ * dZ/dR]
    call dcoul
    !   Compute DERAD = dRAD/dJ * [dJ/dX * dX/dR + dJ/dZ * dZ/dR]
    call dexch
    !
    !   Put all the partial derivatives together:
    !      DEGSDR  = DEQ + DERAD
    !      DEQ   = dE/dQ * DQ/DR
    !      DERAD = dE/dRAD * DRAD/DR
    !
    do 40 i = 1, 3
        degsdr(i) = deq(i) - derad(i)/r2
    40 continue
    !
    !   Add in the derivatives of V3C w.r.t. R(I)
    !
    do 50 i = 1, 3
        degsdr(i) = degsdr(i) + dv3c(i)
    50 continue
    !
    700 continue
    !
    !
    call eunitzero
    if( nder .ne. 0 )then
        call rtocart
        call dedcou
    endif
    !
end subroutine pot

!
!****
!
!> @brief Evaluate SATO for the BrH2_LEPS_LTBZ1_1995 potential.
subroutine sato
    !
    !   This subroutine determines the Sato parameters as a function of the
    !   bond lengths.
    !
    implicit none
    double precision :: anuzero, capr, capr2, capri, cart, cartnu, coschi, dcapr1, dcapr2, dcapr3, dchi1, dchi2, dchi3, &
        & deesdr, degsdr, deijdr, descart, dgscart, dijcart, dz, dzchic, dzchir, easyab, easyac, easybc, engyes, engygs, &
        & engyij, ezero, pengyes, pengygs, pengyij, r, r1p2, r2i, r2p2, r3p2, schip2, tmp, tmp1, tmp2, z, zchi, zhbrp1, &
        & zhbrp2, zhhp1, zhhp2
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
    common /pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/ ezero(isurf+1)
    common /pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /asycm/ easyab, easybc, easyac
    !
    common /satocm/ zhhp1, zhhp2, zhbrp1, zhbrp2, &
        & z(3), dz(3, 3)
    save
    !
    do 10 i = 1, 3
        do 10 j = 1, 3
    10 dz(i, j) = 0.0d0
    !
    r1p2 = r(1) * r(1)
    r2p2 = r(2) * r(2)
    r3p2 = r(3) * r(3)
    r2i  = 1.0d0/r(2)
    zchi = 0.0d0
    !
    !   CAPR IS THE Br TO H2 DISTANCE
    !
    capr2 = 0.5d0 * (r1p2 + r3p2 - 0.5d0 * r2p2)
    if( capr2 .le. 0.d0 )then
        capr2 = 0.0d0
        capr  = 0.0d0
    else
        capr = sqrt(capr2)
        !
        !   Compute Chi, the angle between CAPR and RHH
        !
        coschi = 0.5d0*(r3p2 - r1p2)/(r(2)*capr)
        schip2 = 1.0d0 - coschi*coschi
        zchi   = capr2*(zhbrp1*schip2 + zhbrp2*schip2*schip2)
    endif
    !
    z(1) = 0.1675d0 + zchi
    z(3) = z(1)
    z(2) = zhhp1*r2p2/(zhhp2+r2p2)
    !
    !   Compute the derivative of the Sato parameters with respect to
    !   bond lengths.
    !   The HH Sato parameter depends on R(2) only
    !
    if( nder .ne. 1 ) go to 700
    tmp      = r(2)/(zhhp2 + r2p2)
    dz(2, 2) = 2.0d0*z(2)*(r2i - tmp)
    !
    !   Compute the derivative of ZHBr with respect R(i) if CHI is
    !   is between 0 and 90 degrees.
    !
    if( capr .eq. 0.d0 ) go to 100
    capri = 1.0d0/capr
    !
    !   Compute the derivative of R(Br-HH) with respect to R(i)
    !
    dcapr1 = 0.5d0*r(1)*capri
    dcapr2 = -0.25d0*r(2)*capri
    dcapr3 = 0.5d0*r(3)*capri
    !
    !   Compute the derivative of CHI with respect to R(i)
    !   These derivatives a actually sin(CHI)*dCHI/dR(i)
    !
    tmp1  = coschi*capri
    tmp2  = 1.0d0/(r(2)*capr)
    dchi1 = r(1)*tmp2 + tmp1*dcapr1
    dchi2 = tmp1*(capr + r(2)*dcapr2)*r2i
    dchi3 = tmp1*dcapr3 - r(3)*tmp2
    !
    !   Compute the derivative of ZHBr with respect to CAPR and CHI
    !
    dzchir = 2.0d0*zchi*capri
    dzchic = 2.0d0*capr2*coschi*(zhbrp1 + 2.0d0*zhbrp2*schip2)
    !
    !   Compute the derivative of ZHBr with respect to R(i)
    !
    dz(1, 1) = dzchir*dcapr1 + dzchic*dchi1
    dz(1, 2) = dzchir*dcapr2 + dzchic*dchi2
    dz(1, 3) = dzchir*dcapr3 + dzchic*dchi3
    dz(3, 1) = dz(1, 1)
    dz(3, 2) = dz(1, 2)
    dz(3, 3) = dz(1, 3)
    !
    100 continue
    700 continue
    !
end subroutine sato

!
!*****
!
!> @brief Evaluate V3CFCN for the BrH2_LEPS_LTBZ1_1995 potential.
subroutine v3cfcn
    !
    !   This subroutine computes the energy and the derivatives of the
    !   three-center term for the H-Br-H barrier as a function of the geomet
    !
    implicit none
    double precision :: a1, a2, a3, a4, a5, a6, anuzero, cart, cartnu, costh, deesdr, degsdr, deijdr, descart, dfcn1, &
        & dfcn2, dfcn3, dfcn4, dgscart, dijcart, dist, dthdr1, dthdr2, dthdr3, dv3c, easyab, easyac, easybc, engyes, &
        & engygs, engyij, ezero, fcn1, fcn2, fcn3, fcn4, pengyes, pengygs, pengyij, r, r132p2, r1m3, r1m3p2, r1p3, &
        & r1p3m2, sintp2, sintp6, sintp8, v3c, v3cdth, v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6
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
    common /pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/ ezero(isurf+1)
    common /pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /asycm/ easyab, easybc, easyac
    !
    common /v3ccm/ v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, &
        & v3c, dv3c(3)
    save
    !
    a1 = v3cp1
    a2 = v3cp2
    a3 = v3cp3
    a4 = v3cp4
    a5 = v3cp5
    a6 = v3cp6
    !
    r1p3   = r(1) + r(3)
    r1m3   = r(1) - r(3)
    r1p3m2 = r1p3 - r(2)
    r1m3p2 = r1m3*r1m3
    r132p2 = r1p3m2*r1p3m2
    !
    !   Compute the bond angle from the three internal coordinates
    !
    dist   = r(1)*r(1) + r(3)*r(3) - r(2)*r(2)
    costh  = dist/(2.0d0*r(1)*r(3))
    sintp2 = 1.0d0 - costh*costh
    sintp6 = sintp2*sintp2*sintp2
    sintp8 = sintp6*sintp2
    !
    fcn1 = exp(-a2*r1p3)
    fcn2 = exp(-a3*r1m3p2)
    fcn3 = exp(-a4*r132p2)
    fcn4 = 1.d0 + a5*sintp2 + a6*sintp8
    !
    v3c = a1*fcn1*fcn2*fcn3*fcn4
    !
    !   Compute the derivatives of V3C with respect to the internal coordina
    !   if NDER = 1
    !
    if( nder .ne. 1 ) go to 700
    !
    !   Compute the partial derivatives of the FCN terms
    !
    dfcn1 = -a2
    dfcn2 = -2.0d0*a3*r1m3
    dfcn3 = -2.0d0*a4*r1p3m2
    dfcn4 = 2.0d0*a5*costh + 8.0d0*a6*costh*sintp6
    !
    !   Compute the partial derivative of THETA w.r.t. R(i)
    !   DTHDR(i) = SIN(THETA)*DTHDR(I)
    !
    dthdr1 = costh/r(1) - 1.0d0/r(3)
    dthdr2 = r(2)/(r(1)*r(3))
    dthdr3 = costh/r(3) - 1.0d0/r(1)
    v3cdth = a1*fcn1*fcn2*fcn3*dfcn4
    !
    dv3c(1) = v3c * (dfcn1 + dfcn2 + dfcn3) + v3cdth*dthdr1
    dv3c(2) = -v3c * dfcn3 + v3cdth*dthdr2
    dv3c(3) = v3c * (dfcn1 - dfcn2 + dfcn3) + v3cdth*dthdr3
    !
    700 continue
    !
end subroutine v3cfcn

!
!*****
!
!> @brief Evaluate DCOUL for the BrH2_LEPS_LTBZ1_1995 potential.
subroutine dcoul
    !
    !    This subprogram computes the contribution to the derivative of
    !    the energy with respect to R(i) from the coulombic terms.
    !    That is:
    !            DE/DR(j) = DE/DQ(i) * DQ(i)/DR(j)
    !
    !    In this subprogram the following convention is assumed for the
    !    numbering of the internal coordinates:
    !    For Br + AB -> ABr + B
    !                   R(1) = R(Br-A)
    !                   R(2) = R(A-B)
    !                   R(3) = R(Br-B)
    !
    implicit none
    double precision :: anuzero, beta, cart, cartnu, coul, deesdr, degsdr, deijdr, deq, derad, descart, dgscart, dijcart, &
        & do4z, dqdx, dqdz, dxdr, dz, easyab, easyac, easybc, engyes, engygs, engyij, exch, ezero, op3z, pengyes, &
        & pengygs, pengyij, r, rad, top3z, tzp3, x, z, zhbrp1, zhbrp2, zhhp1, zhhp2, zp3, zpo
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
    common /pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/ ezero(isurf+1)
    common /pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /asycm/ easyab, easybc, easyac
    !
    common /satocm/ zhhp1, zhhp2, zhbrp1, zhbrp2, &
        & z(3), dz(3, 3)
    common /lengcm/ x(3), coul(3), exch(3), rad
    common /lsatcm/ zpo(3), op3z(3), zp3(3), tzp3(3), top3z(3), &
        & do4z(3), beta(3)
    common /ldercm/ dxdr(3, 3), deq(3), derad(3)
    !
    dimension dqdx(3), dqdz(3)
    !
    save
    !
    do 10 i = 1, 3
        deq(i)  = 0.d0
        dqdx(i) = 0.d0
        dqdz(i) = 0.d0
        do 10 j = 1, 3
            dxdr(i, j) = 0.d0
    10 continue
    !
    !   First compute the partial derivatives dQ/dX, dQ/dZ, and
    !   the derivative of X w.r.t. R(i), DX/DR.
    !
    do 20 i = 1, 3
        dqdx(i)    = do4z(i)*(tzp3(i)*x(i) - top3z(i))
        dqdz(i)    = do4z(i)*x(i)*(x(i) - 6.0d0) - coul(i)/zpo(i)
        dxdr(i, i) = -beta(i)*x(i)
    20 continue
    !
    !   Compute the contribution to the derivate of E w.r.t. R(i) from
    !   the coulombic terms, i.e., dE/dR(i) = dQ/dX * DX/DR(i) + dQ/dZ *DZ/D
    !
    do 30 i = 1, 3
        do 30 j = 1, 3
            deq(i) = deq(i) + dqdx(j)*dxdr(j, i) + dqdz(j)*dz(j, i)
    30 continue
    !
end subroutine dcoul

!
!****
!
!> @brief Evaluate DEXCH for the BrH2_LEPS_LTBZ1_1995 potential.
subroutine dexch
    !
    !    This subprogram computes the contribution to the derivative of
    !    the energy with respect to R(i) from the exchange terms.
    !    That is:
    !            DE/DR(j) = DE/DJ(i) * DJ(i)/DR(j)
    !    where
    !            E(LEPS)  = Q1 + Q2 + Q3 - RAD/SQRT(2)
    !
    !    In this subprogram the following convention is assumed for the
    !    numbering of the internal coordinates:
    !    For Br + AB -> ABr + B
    !                   R(1) = R(Br-A)
    !                   R(2) = R(A-B)
    !                   R(3) = R(Br-B)
    !
    implicit none
    double precision :: anuzero, beta, cart, cartnu, coul, deesdr, degsdr, deijdr, deq, derad, descart, dgscart, dijcart, &
        & djdr, djdx, djdz, do4z, draddj, dxdr, dz, easyab, easyac, easybc, engyes, engygs, engyij, exch, ezero, op3z, &
        & pengyes, pengygs, pengyij, r, rad, sumj, top3z, tzp3, x, z, zhbrp1, zhbrp2, zhhp1, zhhp2, zp3, zpo
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
    common /pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/ ezero(isurf+1)
    common /pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common /pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /asycm/ easyab, easybc, easyac
    !
    common /satocm/ zhhp1, zhhp2, zhbrp1, zhbrp2, &
        & z(3), dz(3, 3)
    common /lengcm/ x(3), coul(3), exch(3), rad
    common /lsatcm/ zpo(3), op3z(3), zp3(3), tzp3(3), top3z(3), &
        & do4z(3), beta(3)
    common /ldercm/ dxdr(3, 3), deq(3), derad(3)
    !
    dimension djdx(3), djdz(3), djdr(3,3), draddj(3)
    !
    save
    !
    do 10 i = 1, 3
        derad(i) = 0.d0
        djdx(i)  = 0.d0
        djdz(i)  = 0.d0
    10 continue
    !
    sumj = 0.0d0
    !
    !   First compute the partial derivatives dJ/dX, dJ/dZ
    !
    do 20 i = 1, 3
        sumj    = sumj + exch(i)
        djdx(i) = do4z(i)*(top3z(i)*x(i) - tzp3(i))
        djdz(i) = do4z(i)*x(i)*(3.0d0*x(i) - 2.0d0) - exch(i)/zpo(i)
    20 continue
    !
    !   Compute the partial derivatives of the exchange term, J, w.r.t. R(i)
    !   The outer loop, I, is over the bond lengths, the inner loop is over
    !   each X and Z term.
    !
    do 30 i = 1, 3
        do 30 j = 1, 3
            djdr(j,i) = djdx(j)*dxdr(j, i) + djdz(j)*dz(j, i)
    30 continue
    !
    !   Compute the partial derivatives of the term RAD w.r.t. J(I).
    !
    do 40 i = 1, 3
        draddj(i) = (3.0d0*exch(i) - sumj)/rad
    40 continue
    !
    do 60 i = 1, 3
        do 70 j = 1, 3
            derad(i) = derad(i) + draddj(j)*djdr(j,i)
        70 continue
    60 continue
    !
end subroutine dexch

!
!*****
!
block data ptpacm
    implicit none
    double precision :: anuzero, beta, cart, cartnu, de, do4z, dv3c, dz, easyab, easyac, easybc, ezero, op3z, re, top3z, &
        & tzp3, v3c, v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, z, zhbrp1, zhbrp2, zhhp1, zhhp2, zp3, zpo
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
    common /pt3cm/ ezero(isurf+1)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    common /asycm/ easyab, easybc, easyac
    !
    common /satocm/ zhhp1, zhhp2, zhbrp1, zhbrp2, &
        & z(3), dz(3, 3)
    common /v3ccm/ v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, &
        & v3c, dv3c(3)
    common /lsatcm/ zpo(3), op3z(3), zp3(3), tzp3(3), top3z(3), &
        & do4z(3), beta(3)
    common /parmcm/ de(3), re(3)
    !
    save
    !
    !   Initialize the flags and the I/O unit numbers for the potential
    !
    data nasurf /1, 35*0/
    data nder /0/
    data nflag /1, 1, 15*0, 6, 0, 0/
    !
    data anuzero /0.0d0/
    data icartr, msurf, mder/3, 0, 1/
    data nulbl /25*0/
    data natoms /3/
    !
    !   Initialize the potential parameters
    !
    data de/90.4473d0, 109.559d0, 90.4473d0/
    data re/1.41443d0, 0.74144d0, 1.41443d0/
    data beta/1.81088d0, 1.94198d0, 1.81088d0/
    data zhhp1, zhhp2 /2.91230499d-1, 4.99160974d0/
    data zhbrp1, zhbrp2 /3.84262775d-3, 2.87642885d-3/
    data v3cp1, v3cp2 /1.67457884d0, 7.09518875d-1/
    data v3cp3, v3cp4 /1.87582679d-1, 2.43517074d-1/
    data v3cp5, v3cp6 /2.04145206d-1, 2.57833022d0/
    !
end

!
!*****

!> @brief Evaluate POTINFO for the BrH2_LEPS_LTBZ1_1995 potential.
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
    enddo
    write(nflag(18), '(/)')
    kmax = 0
    do i = 1, isurf+1
        do j = 1, isurf+1
            if( nasurf(i, j) .ne. 0.and.kmax .lt. max(i, j) ) kmax = max(i, j)
        enddo
    enddo
    write(nflag(18), '(2x,'' MAX. AND ACTUAL NO. OF EXCITED SURFACES: '',I3,5x,I3)') msurf, kmax-1
    if( kmax-1 .gt. msurf )then
        write(6, *) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
        stop
    endif
    ksdiag = 0
    kediag = 0
    do i = 2, isurf+1
        if( nasurf(i, i) .ne. 0 )then
            kediag = i-1
            if( ksdiag .eq. 0 ) ksdiag = i-1
        endif
    enddo
    ksoffd = 0
    keoffd = 0
    k      = 0
    do i = 1, isurf
        do j = i+1, isurf+1
            k = k+1
            if( nasurf(i, j)+nasurf(j, i) .ne. 0 )then
                keoffd = k
                if( ksoffd .eq. 0 ) ksoffd = k
            endif
        enddo
    enddo
    write(nflag(18), '(2x,'' MAX. AND ACTUAL ORDER OF DERIVATIVES:    '',I3,5x,I3)') mder, nder
    if( nder .gt. mder )then
        write(6, *) ' WRONG INPUT ON ORDER OF DERIVATIVES'
        stop
    endif
    if( nflag(19) .eq. 1 )then
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
    endif
end subroutine potinfo

!> @brief Evaluate ANCVRT for the BrH2_LEPS_LTBZ1_1995 potential.
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
        enddo
    enddo
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
    enddo
    iblank = ' '
    do ind=1, natoms
        do i=1, 7
            do j=1, 32
                if( indexes(ind) .eq. ianum(i,j) )then
                    isave(ind) = i
                    jsave(ind) = j
                endif
            enddo
        enddo
    enddo
    do ind=1, natoms
        ind2 = nulbl(ind)
        if( ind2 .eq. 0 ) ind2=ind
    enddo
    inc1 = 0
    do ind=1, irctnt-1
        inc1        = inc1+1
        name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
    enddo
    inc2 = 0
    do ind=irctnt, natoms
        inc2        = inc2+1
        name2(inc2) = periodic_1(isave(ind),jsave(ind))(:2)
    enddo
    if( nflag(1) .eq. 2 ) distance = 'ANGSTROMS           '
    if( nflag(2) .eq. 2 )then
        units = 'MILLIHARTREE        '
    else if( nflag(2) .eq. 3 )then
        units = 'EV                  '
    else if( nflag(2) .eq. 4 )then
        units = 'KCAL PER MOLE       '
    else if( nflag(2) .eq. 5 )then
        units = 'WAVENUMBERS         '
    else if( nflag(2) .eq. 6 )then
        units = 'KILOJOULES PER MOLE '
    endif
    cnvrtd  = 1.d0
    cnvrte  = 1.d0
    cnvrtde = 1.d0
    if( nflag(1) .eq. 2 ) cnvrtd = bohr_a
    if( nflag(2) .eq. 2 )then
        cnvrte = cnvrte*htomillh
    else if( nflag(2) .eq. 3 )then
        cnvrte = cnvrte*htoev
    else if( nflag(2) .eq. 4 )then
        cnvrte = cnvrte*htokcal
    else if( nflag(2) .eq. 5 )then
        cnvrte = cnvrte*htowave
    else if( nflag(2) .eq. 6 )then
        cnvrte = cnvrte*htokj
    endif
    cnvrtde = cnvrte/cnvrtd
    isum    = 0
    do inu=1, 25
        isum = isum + nulbl(inu)
    enddo
    ireorder = 0
    if( isum .ne. 0 ) ireorder = 1
end subroutine ancvrt

!> @brief Evaluate CARTOU for the BrH2_LEPS_LTBZ1_1995 potential.
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
    call convert_cartesian_input(cart, nulbl, natoms, cnvrtd, ireorder == 1, cartnu)
end subroutine cartou

!> @brief Evaluate CARTTOR for the BrH2_LEPS_LTBZ1_1995 potential.
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
    call cartesian_to_potlib_coordinates(cartnu, natoms, icartr, r, coordinate_status)
    if( coordinate_status /= 0 )then
        write(nflag(18), '(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') icartr
        stop
    endif
end subroutine carttor

!> @brief Evaluate EUNITZERO for the BrH2_LEPS_LTBZ1_1995 potential.
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
    if( ksdiag .ne. 0 )then
        do i=ksdiag, kediag
            pengyes(i) = engyes(i) * cnvrte - anuzero
        enddo
    endif
    if( ksoffd .ne. 0 )then
        do j=ksoffd, keoffd
            pengyij(j) = engyij(j) * cnvrte
        enddo
    endif
end subroutine eunitzero

!> @brief Evaluate RTOCART for the BrH2_LEPS_LTBZ1_1995 potential.
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
    if( icartr .eq. 1 )then
        do i = 1, natoms
            ind             = 3*i-2
            dgscartnu(i, 1) = degsdr(ind)
            dgscartnu(i, 2) = degsdr(ind+1)
            dgscartnu(i, 3) = degsdr(ind+2)
            if( ksdiag .ne. 0 )then
                do j = ksdiag, kediag
                    descartnu(i, 1, j) = deesdr(ind, j)
                    descartnu(i, 2, j) = deesdr(ind+1, j)
                    descartnu(i, 3, j) = deesdr(ind+2, j)
                enddo
            endif
            if( keoffd .ne. 0 )then
                do k = ksoffd, keoffd
                    dijcartnu(i, 1, k) = deijdr(ind, k)
                    dijcartnu(i, 2, k) = deijdr(ind+1, k)
                    dijcartnu(i, 3, k) = deijdr(ind+2, k)
                enddo
            endif
        enddo
    else if( icartr .eq. 2 )then
        do i = 1, natoms
            dgscartnu(i, 1) = 0.d0
            dgscartnu(i, 2) = 0.d0
            dgscartnu(i, 3) = 0.d0
            if( ksdiag .ne. 0 )then
                do j1=ksdiag, kediag
                    descartnu(i, 1, j1) = 0.d0
                    descartnu(i, 2, j1) = 0.d0
                    descartnu(i, 3, j1) = 0.d0
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do j2=ksoffd, keoffd
                    dijcartnu(i, 1, j2) = 0.d0
                    dijcartnu(i, 2, j2) = 0.d0
                    dijcartnu(i, 3, j2) = 0.d0
                enddo
            endif
            do j = 1, natoms
                if( j .lt. i )then
                    m1 = natoms*(j-1) - (j*(j-1))/2 + i-j
                else if( j .gt. i )then
                    m1 = natoms*(i-1) - (i*(i-1))/2 + j-i
                else
                    go to 20
                endif
                y               = degsdr(m1)
                termx           = (cartnu(i, 1)-cartnu(j, 1))/r(m1)
                termy           = (cartnu(i, 2)-cartnu(j, 2))/r(m1)
                termz           = (cartnu(i, 3)-cartnu(j, 3))/r(m1)
                dgscartnu(i, 1) = dgscartnu(i, 1) + termx*y
                dgscartnu(i, 2) = dgscartnu(i, 2) + termy*y
                dgscartnu(i, 3) = dgscartnu(i, 3) + termz*y
                if( ksdiag .gt. 0 )then
                    y = deesdr(m1, j1)
                    do j1=ksdiag, kediag
                        descartnu(i, 1, j1) = descartnu(i, 1, j1) + termx*y
                        descartnu(i, 2, j1) = descartnu(i, 2, j1) + termy*y
                        descartnu(i, 3, j1) = descartnu(i, 3, j1) + termz*y
                    enddo
                else if( ksoffd .gt. 0 )then
                    do j2=ksoffd, keoffd
                        y                   = deijdr(m1, j2)
                        dijcartnu(i, 1, j2) = dijcartnu(i, 1, j2) + termx*y
                        dijcartnu(i, 2, j2) = dijcartnu(i, 2, j2) + termy*y
                        dijcartnu(i, 3, j2) = dijcartnu(i, 3, j2) + termz*y
                    enddo
                endif
                20 continue
            enddo
        enddo
    else if( icartr .eq. 3 )then
        do i = 1, natoms
            ygs(i) = degsdr(i)/r(i)
            if( ksdiag .ne. 0 )then
                do j=ksdiag, kediag
                    yes(i,j) = deesdr(i, j)/r(i)
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do k=ksoffd, keoffd
                    yij(i,k) = deijdr(i, k)/r(i)
                enddo
            endif
        enddo
        do k = 1, 3
            term12          = cartnu(1, k)-cartnu(2, k)
            term23          = cartnu(2, k)-cartnu(3, k)
            term13          = cartnu(1, k)-cartnu(3, k)
            dgscartnu(1, k) = term12*ygs(1) + term13*ygs(3)
            dgscartnu(2, k) = -term12*ygs(1) + term23*ygs(2)
            dgscartnu(3, k) = -term13*ygs(3) - term23*ygs(2)
            if( ksdiag .ne. 0 )then
                do j1=ksdiag, kediag
                    descartnu(1, k, j1) = term12*yes(1,j1) + term13*yes(3,j1)
                    descartnu(2, k, j1) = -term12*yes(1,j1) + term23*yes(2,j1)
                    descartnu(3, k, j1) = -term13*yes(3,j1) - term23*yes(2,j1)
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do j2=ksoffd, keoffd
                    dijcartnu(1, k, j2) = term12*yij(1,j2) + term13*yij(3,j2)
                    dijcartnu(2, k, j2) = -term12*yij(1,j2) + term23*yij(2,j2)
                    dijcartnu(3, k, j2) = -term13*yij(3,j2) - term23*yij(2,j2)
                enddo
            endif
        enddo
    else if( icartr .eq. 4 )then
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
    endif
end subroutine rtocart

!> @brief Evaluate DEDCOU for the BrH2_LEPS_LTBZ1_1995 potential.
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
    if( ireorder .eq. 1 )then
        do i = 1, natoms
            dgscart(i, 1) = dgscartnu(nulbl(i), 1) * cnvrtde
            dgscart(i, 2) = dgscartnu(nulbl(i), 2) * cnvrtde
            dgscart(i, 3) = dgscartnu(nulbl(i), 3) * cnvrtde
            if( ksdiag .ne. 0 )then
                do j=ksdiag, kediag
                    descart(i, 1, j) = descartnu(nulbl(i), 1, j) * cnvrtde
                    descart(i, 2, j) = descartnu(nulbl(i), 2, j) * cnvrtde
                    descart(i, 3, j) = descartnu(nulbl(i), 3, j) * cnvrtde
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do k=ksoffd, keoffd
                    dijcart(i, 1, k) = dijcartnu(nulbl(i), 1, k) * cnvrtde
                    dijcart(i, 2, k) = dijcartnu(nulbl(i), 2, k) * cnvrtde
                    dijcart(i, 3, k) = dijcartnu(nulbl(i), 3, k) * cnvrtde
                enddo
            endif
        enddo
    else
        do i = 1, natoms
            dgscart(i, 1) = dgscartnu(i, 1) * cnvrtde
            dgscart(i, 2) = dgscartnu(i, 2) * cnvrtde
            dgscart(i, 3) = dgscartnu(i, 3) * cnvrtde
            if( ksdiag .ne. 0 )then
                do j=ksdiag, kediag
                    descart(i, 1, j) = descartnu(i, 1, j) * cnvrtde
                    descart(i, 2, j) = descartnu(i, 2, j) * cnvrtde
                    descart(i, 3, j) = descartnu(i, 3, j) * cnvrtde
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do k=ksoffd, keoffd
                    dijcart(i, 1, k) = dijcartnu(i, 1, k) * cnvrtde
                    dijcart(i, 2, k) = dijcartnu(i, 2, k) * cnvrtde
                    dijcart(i, 3, k) = dijcartnu(i, 3, k) * cnvrtde
                enddo
            endif
        enddo
    endif
end subroutine dedcou
