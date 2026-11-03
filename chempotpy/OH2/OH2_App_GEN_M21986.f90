!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the OH2_App_GEN_M21986 potential.
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
    v = 0.d0; g = 0.d0; d = 0.d0

    nder   = igrad
    nasurf = 1
    cart   = 0.d0
    do iatom=1, natoms
        do idir=1, 3
            cart(iatom, idir) = x(iatom,idir)/legacy_bohr_to_angstrom
        end do
    end do

    if(first_time_data) then
        call prepot
        first_time_data = .false.
    end if

    call pot

    if (igrad==0) then
        do istate=1, nstates
            p(istate) = pengygs*legacy_hartree_to_electronvolt
        end do
    else if (igrad==1) then
        do istate=1, nstates
            p(istate) = pengygs*legacy_hartree_to_electronvolt
        end do
        do iatom=1, natoms
            do idir=1, 3
                g(1,iatom,idir) = dgscart(iatom, idir)*51.422067
            end do
        end do
    else if (igrad==2) then
        write (*, *) 'Only energy and gradient are available'
    end if

end subroutine pes

!
!> @brief Evaluate PREPOT for the OH2_App_GEN_M21986 potential.
subroutine prepot
    !
    !   System:          OH2  A'' surface
    !   Common name:     M2
    !   Functional form: Rotated-Morse-oscillator-spline (RMOS) potential
    !                    for collinear geometeries plus anti-Morse
    !                    bending potential.
    !   References:      B. C. Garrett and D. G. Truhlar
    !                    Int. J. Quantum Chem. 29, 1463 (1986)
    !
    !   Interface:                     potlib2001
    !   Number of bodies:              3
    !   Number of derivatives:         1
    !   Number of electronic surfaces: 1
    !
    !   PREPOT must be called once before any calls to POT.
    !   The control flags for the potential are contained in the block data
    !   subprogram PTPACM.
    !   Coordinates, potential energy, and derivatives are passed
    !   The potential energy in the three asymptotic valleys are
    !   stored in the common block ASYCM:
    !                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC
    !   The potential energy in the AB valley, EASYAB, is equal to the poten
    !   energy of the H "infinitely" far from the OH diatomic, with the
    !   OH diatomic at its equilibrium configuration.  Similarly, the terms
    !   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
    !   respectively.
    !   All the information passed through the common blocks PT1CM and ASYCM
    !   is in hartree atomic units.
    !
    !        This potential is written such that:
    !                       R(1) = R(O-H)
    !                       R(2) = R(H-H)
    !                       R(3) = R(H-O)
    !   The classical potential energy is set equal to zero for the O
    !   infinitely far from the equilibrium H2 diatomic.
    !
    !   The flags that indicate what calculations should be carried out in
    !   the potential routine are passed through the common block PT2CM:
    !   where:
    !        NASURF - which electronic states availalble
    !                 (1,1) = 1 as only gs state available
    !        NDER  = 0 => no derivatives should be calculated
    !        NDER  = 1 => calculate first derivatives
    !        NFLAG - these integer values can be used to flag options
    !                within the potential;
    !
    !
    !   Potential parameters' default settings
    !                  Variable            Default value
    !                  NDER                1
    !                  NFLAG(18)           6
    !
    !
    implicit none
    double precision :: a, anuzero, beta, cart, cartnu, d, dbcoor, de, easyab, easyac, easybc, ezero, gamp, rbb, req, &
        & rhh
    integer :: i, icartr, indexes, irctnt, it, mder, msurf, nasurf, &
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
    double precision :: nab, nab1, nbc, nbc1
    !
    !   Conversion for Angstroms to Bohr (Angstroms/Bohr)
    !
    double precision, parameter :: bohr = .52917706d0
    !
    !   Conversion for kcal to Hartree (kcal/Hartree)
    !
    double precision, parameter :: hartre = 627.5095d0
    !
    !  Vaule of tolerance
    !
    double precision, parameter :: tol = 1.d-05
    !
    common /ambpcm/ d(3), beta, req(3), a, gamp(3), rbb, rhh, de
    common /precm/  dbcoor(2)
    !
    !   Echo the name of the potential energy surface
    !
    if(natoms.gt.25) then
        write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    end if
    !
    !         WRITE (NFLAG(18), 600)
    !
    !   Assign the parameters for the collinear part of the potential
    !
    call prepot2
    !
    !   Echo the parameters for the bending potential
    !
    !      WRITE (NFLAG(18), 620) REQ,D
    !      WRITE (NFLAG(18), 630) BETA,A,RBB,GAMP
    !
    !       CONVERT TO ATOMIC UNITS
    !
    do 50 it = 1, 3
        req(it) = req(it)/bohr
        d(it)   = d(it)/hartre
    50 continue
    !
    beta = beta * bohr
    de   = d(1)
    a    = a/bohr
    rhh  = 0.74144d0 / bohr
    rbb  = rbb /bohr
    !
    !   Initialize the energy in the asymptotic valleys
    !
    easyab = d(1)
    easybc = d(2)
    easyac = d(3)
    !
    !
    !
    ezero(1) = d(2)
    !
    do i=1, 5
        ref(i) = ' '
    end do
    !
    ref(1) = 'B. C. Garrett, D. G. Truhlar'
    ref(2) = 'Int. J. Quantum Chem. 29, 1463(1986)'
    !
    indexes(1) = 8
    indexes(2) = 1
    indexes(3) = 1
    !
    !
    !
    irctnt = 2
    !
    !      CALL POTINFO
    !
    call ancvrt
    !
end subroutine prepot

!
!> @brief Evaluate POT for the OH2_App_GEN_M21986 potential.
subroutine pot
    !
    !   The potential energy in the AB valley, EASYAB, is equal to the poten
    !   energy of the H "infinitely" far from the OH diatomic, with the
    !   OH diatomic at its equilibrium configuration.  Similarly, the terms
    !   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
    !   respectively.
    !
    !        This potential is written such that:
    !                       R(1) = R(O-H)
    !                       R(2) = R(H-H)
    !                       R(3) = R(H-O)
    !   The classical potential energy is set equal to zero for the O
    !   infinitely far from the equilibrium H2 diatomic.
    !
    !      ENTRY POT
    implicit none
    double precision :: a, anuzero, bcoor, beta, c, cart, cartnu, d, dac, dbcoor, de, deesdr, degsdr, deijdr, descart, &
        & dex1, dex2, dgam, dgscart, dijcart, easyab, easyac, easybc, engyes, engygs, engyij, ex, ex1, ex2, ezero, gamma, &
        & gamp, pengyes, pengygs, pengyij, r, raccol, ractmp, rbb, req, rhh, rob, stuff, t
    integer :: icartr, indexes, irctnt, it, mder, msurf, nasurf, natoms, &
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
    double precision :: nab, nab1, nbc, nbc1
    !
    !   Conversion for Angstroms to Bohr (Angstroms/Bohr)
    !
    double precision, parameter :: bohr = .52917706d0
    !
    !   Conversion for kcal to Hartree (kcal/Hartree)
    !
    double precision, parameter :: hartre = 627.5095d0
    !
    !  Vaule of tolerance
    !
    double precision, parameter :: tol = 1.d-05
    !
    common /ambpcm/ d(3), beta, req(3), a, gamp(3), rbb, rhh, de
    common /precm/  dbcoor(2)
    !
    call cartou
    call carttor
    !
    !   Check the values of NASURF and NDER for validity.
    !
    !      IF (NASURF(1,1) .EQ. 0) THEN
    !         WRITE(NFLAG(18), 900) NASURF(1,1)
    !         STOP
    !      ENDIF
    !         IF (NDER .GT. 1) THEN
    !             WRITE (NFLAG(18), 910) NDER
    !             STOP 'POT 2'
    !         ENDIF
    !
    if (r(1) .gt. tol .and. r(2) .gt. tol .and. r(3) .gt. tol) then
        raccol = r(1) + r(2)
        bcoor  = 0.0d0
        if (nder .eq. 1) then
            dbcoor(1) = 0.0d0
            dbcoor(2) = 0.0d0
            dac       = 0.0d0
        end if
        rob = 1.0d0
        if(abs(rhh-rbb).gt.1.d-10) then
            !
            !     FIRST CALCULATE SOME USEFUL NUMBERS
            !     NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
            !     AB AND BC RESPECTIVELY. (EQUATION 2.2)
            !
            nab = exp((req(1) - r(1))/a)
            if (nab.lt.1.d-15) nab = 1.d-15
            nbc = exp((req(2) - r(2))/a)
            if (nbc.lt.1.d-15) nbc = 1.d-15
            !
            !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
            !       SEE EQUATION 2.7A, 2.7B
            !
            c = 1.d0 - nab/nbc
            if (abs(c) .lt. 1.d-14) then
                nab1 = .5d0
                nbc1 = .5d0
            else
                c    = c/nab
                nab1 = (2.d0 + c - sqrt(4.d0+c*c))/(2.d0*c)
                nbc1 = 1.d0 - nab1
            end if
            !
            !          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
            !          FIRST TRAP OUT ANY ZERO ARGUMENTS
            !
            if (nab1*nbc1 .le. 0.0d0) then
                rob = 1.d0
            else
                stuff = a * log(nab1*nbc1)
                t     = rhh - stuff
                rob   = 1.0d0
                if(abs(t).gt.1.d-15) rob = (rbb - stuff)/t
            end if
        end if
        !
        !        CALCULATE BENDING CORRECTION
        !        EVALUATE GAMMA AND DERIVATIVE
        !
        ex    = exp(raccol*(gamp(2) + gamp(3)*raccol))
        gamma = gamp(1) + ex
        if (nder .eq. 1) dgam = (gamp(2) + 2.d0*gamp(3)*raccol)*ex
        gamma = gamma*de
        if (nder .eq. 1) dgam = dgam*de
        !
        !           CALCULATE V(R(3)) TERM
        !
        ex  = exp(-beta*(r(3)-req(3))/rob)
        ex1 = ex*(1.d0 + 0.5d0*ex)
        if (nder .eq. 1) dex1 = ex*(1.d0 + ex)
        !
        !           CALCULATE V(R(1)+R(2)) TERM
        !
        ex  = exp(beta*(req(3)-r(2)-r(1))/rob)
        ex2 = ex*(1.d0 + 0.5d0*ex)
        if (nder .eq. 1) dex2 = ex*(1.d0 + ex)
        !
        !        HERE IS THE BENDING CORRECTION
        !
        bcoor = gamma*(ex1 - ex2)
        !
        !        WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
        !        OF BCOOR WITH RESPECT TO R(1)  R(2) AND R(3).
        !        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
        !
        if (nder .eq. 1) then
            gamma = gamma*beta/rob
            dac   = -gamma*dex1
            !
            !          CALCULATE CORRECTIONS TO DAB, DBC
            !
            dbcoor(1) = gamma*dex2 + dgam*(ex1-ex2)
            dbcoor(2) = dbcoor(1)
        end if
        !
        !       NOW CALCULATE THE COLLINEAR ENERGY
        !       STORE R(3) IN RACTMP THEN SET R(3) TO
        !       THE COLLINEAR GEOMETRY AND CALL POT2
        !
        ractmp = r(3)
        r(3)   = r(2) + r(1)
        !
        call pot2
        !
        r(3) = ractmp
        !
        engygs = engygs + bcoor
        !
        if (nder .eq. 1) then
            do 300 it = 1, 2
                degsdr(it) = dbcoor(it) + degsdr(it) + degsdr(3)
            300 continue
            degsdr(3) = dac
        end if
    else
        engygs = 1.d30
        if (nder .eq. 1) then
            degsdr(1) = 0.0d0
            if (r(1) .le. tol) degsdr(1) = -engygs
            degsdr(2) = 0.0d0
            if (r(1) .le. tol) degsdr(2) = -engygs
            degsdr(3) = 0.0d0
            if (r(1) .le. tol) degsdr(3) = -engygs
        end if
    end if
    !
    call eunitzero
    if(nder.ne.0) then
        call rtocart
        call dedcou
    end if
    !
    !
end subroutine pot

!
!> @brief Evaluate PREPOT2 for the OH2_App_GEN_M21986 potential.
subroutine prepot2
    !
    !  COLLINEAR A + BC POTENTIAL FUNCTION
    !  RMO-SPLINE TYPE FIT
    !  COORDINATES AND ENERGIES PASSED TO POTENTIAL IN ATOMIC UNITS
    !
    implicit none
    double precision :: a1dif, a1inf, a2dif, a2inf, aer, alph1, alph2, anuzero, b1dif, b1inf, b2dif, b2inf, bet, bsp, &
        & cart, cartnu, csp, d1dif, d1inf, d2dif, d2inf, dasy, dis, easyab, easyac, easybc, ezero, phi, r1dif, r1inf, &
        & r1s, r2dif, r2inf, r2s, req2, t, t1, vsp
    integer :: i, iayes, icartr, ii, indexes, irctnt, j, mder, msurf, nasurf, natoms, nder, nflag, nnn, nphi, nphip2, &
        & nulbl
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
    double precision, parameter :: detrad = 0.0174532d0
    double precision, parameter :: ckcal  = 627.5095d0
    double precision, parameter :: bohr   = 0.529177d0
    !
    common /pre2cm/ d1inf, r1inf, b1inf, d2inf, &
        & r2inf, b2inf, a1inf, a2inf, &
        & r1s, r2s, vsp, alph1, alph2, &
        & phi(20), bsp(20, 3), csp(22, 3), &
        & d1dif, r1dif, b1dif, a1dif, &
        & d2dif, r2dif, b2dif, a2dif, &
        & dasy, iayes, nphi
    !
    character(len=20) :: title
    !
    !  TITLE
    !
    title = ' O+H2 MAB 9-12-79 '
    !
    !  ASSIGN COLLINEAR PARAMETERS (UNITS OF KCAL, ANGSTROM, DEGREES)
    !  IN BLOCK DATA SUBROUTINE
    !
    a1inf = 0.0d0
    a2inf = 0.0d0
    iayes = 0
    !
    !  ASSIGN PARAMETERS FOR THE ASYMPTOTIC REACTANT AND PRODUCT REGIONS
    !  IN BLOCK DATA SUBROUTINE
    !
    if(a2inf .eq.0.0d0) a2inf = d1inf-d2inf
    !
    !  ASSIGN THE PIVOT POINT, POTENTIAL AT THE PIVOT,
    !  FLAG, AND ASYMPTOTIC EXPONENTIAL PARAMETERS
    !  IN BLOCK DATA SUBROUTINE.  ECHO.
    !
    !      WRITE (NFLAG(18), 600) TITLE
    !      WRITE (NFLAG(18), 601) D1INF,R1INF,B1INF,D2INF,R2INF,B2INF,A1INF,
    !     2INF
    !      WRITE (NFLAG(18), 602) R1S,R2S,VSP,ALPH1,ALPH2
    !
    !      IF(IAYES.NE.0) WRITE (NFLAG(18), 603)
    !      IF(IAYES.EQ.0) WRITE (NFLAG(18), 604)
    !
    nphip2 = nphi + 2
    !
    !  ECHO POTENTIAL PARAMTERS ASSIGNED BY BLOCK DATA SUBROUTINE
    !
    !      WRITE (NFLAG(18), 605) NPHI
    !      WRITE (NFLAG(18), 610)
    !
    !      DO 30 I = 1,NPHI
    !         WRITE (NFLAG(18), 606) I,PHI(I),(BSP(I,J),CSP(I,J),J=1,3)
    !   30 CONTINUE
    !
    !      WRITE (NFLAG(18), 611)
    !      DO 40 I=1,2
    !         K = I + NPHI
    !         WRITE (NFLAG(18), 607) K,(CSP(K,J),J=1,3)
    !   40 CONTINUE
    !
    !  CONVERT TO ATOMIC UNITS
    !
    d1inf = d1inf/ckcal
    a1inf = a1inf/ckcal
    r1inf = r1inf/bohr
    b1inf = b1inf*bohr
    d2inf = d2inf/ckcal
    a2inf = a2inf/ckcal
    r2inf = r2inf/bohr
    b2inf = b2inf*bohr
    r1s   = r1s/bohr
    r2s   = r2s/bohr
    vsp   = vsp/ckcal
    alph1 = alph1*bohr
    alph2 = alph2*bohr
    do 50 i = 1, nphi
        bsp(i, 1) = bsp(i, 1)/ckcal
        bsp(i, 2) = bsp(i, 2)/bohr
        bsp(i, 3) = bsp(i, 3)*bohr
    50 continue
    nnn = nphi + 2
    do 60 i = 1, nnn
        csp(i, 1) = csp(i, 1)/ckcal
        csp(i, 2) = csp(i, 2)/bohr
        csp(i, 3) = csp(i, 3)*bohr
    60 continue
    !
    !  CONVERT FROM DEGREES TO RADIANS
    !
    t = 1.0d0/detrad
    do 70 j = 1, 3
        csp(2, j) = csp(2, j)*t
    70 continue
    t = t**3
    do 80 ii = 1, nphi
        phi(ii) = phi(ii)*detrad
        i       = ii + 2
        do 80 j = 1, 3
            csp(i, j) = csp(i, j)*t
    80 continue
    !
    !  SET UP CONSTANTS
    !
    !     REACTANT ASYMPTOTE
    !
    dis   = bsp(1, 1)
    req2  = bsp(1, 2)
    bet   = bsp(1, 3)
    t1    = 1.0d0 - exp(-bet*req2)
    aer   = vsp - dis*t1*t1
    req2  = r2s - req2
    d1dif = dis - d1inf
    r1dif = req2 - r1inf
    b1dif = bet - b1inf
    a1dif = aer - a1inf
    !
    !     PRODUCT ASYMPTOTE
    !
    dis   = bsp(nphi, 1)
    req2  = bsp(nphi, 2)
    bet   = bsp(nphi, 3)
    t1    = 1.0d0 - exp(-bet*req2)
    aer   = vsp - dis*t1*t1
    req2  = r1s - req2
    d2dif = dis - d2inf
    r2dif = req2 - r2inf
    b2dif = bet - b2inf
    a2dif = aer - a2inf
    !
    !  DEFINE ZERO OF ENERGY
    !
    dasy = d1inf
    !
    !
end subroutine prepot2

!
!> @brief Evaluate POT2 for the OH2_App_GEN_M21986 potential.
subroutine pot2
    !
    !      ENTRY POT2
    implicit none
    double precision :: a1dif, a1inf, a2dif, a2inf, aer, alph1, alph2, angle, anuzero, b1dif, b1inf, b2dif, b2inf, bet, &
        & bsp, cart, cartnu, csp, d1dif, d1inf, d2dif, d2inf, daer, dasy, dbet, ddis, deesdr, degsdr, deijdr, delr, &
        & descart, dgscart, dijcart, dis, dphd1, dphd2, dr1, dr2, dreq, drrd1, drrd2, dvdp, easyab, easyac, easybc, &
        & engyes, engygs, engyij, ezero, pengyes, pengygs, pengyij, phi, r, r1, r1dif, r1inf, r1s, r2, r2dif, r2inf, r2s, &
        & r3, req2, rr, scrap, spline, t1, t2, t3, vsp
    integer :: iayes, icartr, indexes, irctnt, mder, msurf, nasurf, &
        & natoms, nder, nflag, nphi, nulbl
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
    double precision, parameter :: detrad = 0.0174532d0
    double precision, parameter :: ckcal  = 627.5095d0
    double precision, parameter :: bohr   = 0.529177d0
    !
    common /pre2cm/ d1inf, r1inf, b1inf, d2inf, &
        & r2inf, b2inf, a1inf, a2inf, &
        & r1s, r2s, vsp, alph1, alph2, &
        & phi(20), bsp(20, 3), csp(22, 3), &
        & d1dif, r1dif, b1dif, a1dif, &
        & d2dif, r2dif, b2dif, a2dif, &
        & dasy, iayes, nphi
    !
    character(len=20) :: title
    !
    r1 = r(1)
    r2 = r(2)
    r3 = r(3)
    if (r1 .gt. r1s) then
        if (r2 .gt. r2s) then
            !
            !   3-ATOM BREAK UP REGION
            !
            engygs = dasy
            if (nder .eq. 1) then
                degsdr(1) = 0.0d0
                degsdr(2) = 0.0d0
            end if
        else
            !
            !   REACTANT ASYMPTOTIC REGION
            !
            dr1 = r1 - r1s
            if (dr1*alph1 .lt. -70.0d0) write (nflag(18), '(/,2X,T5,''R1, R2, R3, DELR, ALPH1:'',/,2X,T10,5(1PE13.5,1' // &
                'X))') r1, r2, &
                !     R3,
                & dr1, alph1
            t1 = exp(-dr1*alph1)
            t2 = -alph1*t1
            !
            !.....INTERPOLATE FOR MORSE PARAMETERS
            !
            dis  = d1inf + d1dif*t1
            ddis = d1dif*t2
            req2 = r1inf + r1dif*t1
            dreq = r1dif*t2
            bet  = b1inf + b1dif*t1
            dbet = b1dif*t2
            if(iayes.eq.0) then
                aer  = vsp - dis
                daer = -ddis
            else
                aer  = a1inf + a1dif*t1
                daer = a1dif*t2
            end if
            delr = r2 - req2
            if (delr*bet .lt.  -70.0d0) &
                & write (nflag(18), '(/,2X,T5,''R1, R2, R3, DR2, BET:'',/,2X,T10,5(1PE13.5,1X))') r1, r2, r3, delr, bet
            t1     = exp(-bet*delr)
            t2     = 1.0d0 - t1
            engygs = aer + dis*t2*t2
            if (nder .eq. 1) then
                t3        = 2.0d0*dis*t1
                degsdr(1) = daer + t2*(ddis*t2 - t3*(bet*dreq - dbet*delr))
                degsdr(2) = t3*bet*t2
            end if
        end if
    else if (r2 .gt. r2s) then
        !
        !   PRODUCT ASYMPTOTIC REGION
        !
        dr2 = r2 - r2s
        if (dr2*alph2 .lt.  -70.d0) &
            & write (nflag(18), '(/,2X,T5,''R1, R2, R3, DELR, ALPH2:'',/,2X,T10,5(1PE13.5,1X))') r1, r2, r3, dr2, alph2
        !
        !.....INTERPOLATE FOR MORSE PARAMETERS
        !
        t1   = exp(-dr2*alph2)
        t2   = -alph2*t1
        dis  = d2inf + d2dif*t1
        ddis = d2dif*t2
        req2 = r2inf + r2dif*t1
        dreq = r2dif*t2
        bet  = b2inf + b2dif*t1
        dbet = b2dif*t2
        if(iayes.eq.0) then
            aer  = vsp - dis
            daer = -ddis
        else
            aer  = a2inf + a2dif*t1
            daer = a2dif*t2
        end if
        delr = r1 - req2
        if (delr*bet .lt.  -70.0d0) &
            & write (nflag(18), '(/,2X,T5,''R1, R2, R3, DR2, BET:'',/,2X,T10,5(1PE13.5,1X))') r1, r2, r3, delr, bet
        t1     = exp(-bet*delr)
        t2     = 1.0d0 - t1
        engygs = aer + dis*t2*t2
        if (nder .eq. 1) then
            t3        = 2.0d0*dis*t1
            degsdr(2) = daer + t2*(ddis*t2 - t3*(bet*dreq - dbet*delr))
            !
            degsdr(1) = t3*bet*t2
        end if
    else
        !
        !   SPLINE INTERPOLATION REGION
        !
        dr1   = r1s - r1
        dr2   = r2s - r2
        t1    = dr1*dr1 + dr2*dr2
        rr    = sqrt(t1)
        drrd1 = -dr1/rr
        drrd2 = -dr2/rr
        angle = atan(dr1/dr2)
        dphd1 = -dr2/t1
        dphd2 = dr1/t1
        dis   = spline(3, nphi, phi, bsp(1, 1), csp(1, 1), scrap, angle)
        ddis  = spline(4, nphi, phi, bsp(1, 1), csp(1, 1), scrap, angle)
        req2  = spline(3, nphi, phi, bsp(1, 2), csp(1, 2), scrap, angle)
        dreq  = spline(4, nphi, phi, bsp(1, 2), csp(1, 2), scrap, angle)
        bet   = spline(3, nphi, phi, bsp(1, 3), csp(1, 3), scrap, angle)
        dbet  = spline(4, nphi, phi, bsp(1, 3), csp(1, 3), scrap, angle)
        if(iayes.eq.0) then
            aer  = vsp - dis
            daer = -ddis
        else
            t1   = exp(-bet*req2)
            t2   = 1.0d0 - t1
            aer  = vsp - dis*t2*t2
            daer = -t2*(ddis*t2 + 2.0d0*dis*t1*(bet*dreq + req2*dbet))
        end if
        if (bet*(rr-req2) .gt.  70) &
            & write (nflag(18), '(/,2X,T5,''R1, R2, R3, RR, ANGLE, DIS, REQ2, BET:'', /,(2X,T5,1PE20.10,T30,1PE20.10,' // &
            'T55,1PE20.10))') r1, r2, r3, rr, angle, dis, req2, bet
        delr   = rr - req2
        t1     = exp(bet*delr)
        t2     = 1.0d0 - t1
        engygs = aer + dis*t2*t2
        if(nder.eq.1) then
            t3        = 2.0d0*dis*t1
            dvdp      = daer + t2*(ddis*t2 - t3*(dbet*delr - bet*dreq))
            t3        = t3*t2*bet
            degsdr(1) = dvdp*dphd1 - t3*drrd1
            degsdr(2) = dvdp*dphd2 - t3*drrd2
        end if
    end if
    if(nder.eq.1) then
        degsdr(1) = degsdr(1)
        degsdr(2) = degsdr(2)
        degsdr(3) = 0.0d0
    end if
    !
    !
end subroutine pot2

!
!> @brief Evaluate SPLINE for the OH2_App_GEN_M21986 potential.
function spline( isw, nn, x, y, c, d, xpt )
    implicit none
    double precision :: anuzero, c, cart, cartnu, cfact, d, deesdr, degsdr, deijdr, deriv, descart, dfact, dgscart, &
        & dijcart, easyab, easyac, easybc, engyes, engygs, engyij, ezero, pengyes, pengygs, pengyij, pivot, r, spline, &
        & supd, v, v2, x, xpt, y, z
    integer :: i, icartr, indexes, irctnt, isw, j, mder, msurf, n, nasurf, natoms, nder, nflag, nn, np1, np2, nulbl
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
    !     PARAMETER (IPRT = 6)
    dimension x(nn), y(nn), c(nn+2), d(nn+2)
    !
    !  THIS IS A SUBROUTINE FOR FITTING DATA WITH A CUBIC SPLINE
    !  POLYNOMIAL AND EVALUATING THAT POLYNOMIAL AT A GIVEN POINT
    !  OR ITS DERIVATIVE AT A GIVEN POINT
    !
    !  CALLING SEQUENCE .......
    !     ISW ... CONTROL OPTION
    !         ISW=1  IF A CUBIC SPLINE IS TO BE FITTED TO THE SET OF KNOTS
    !                DEFINED BY THE ARRAYS X AND Y.  THE SPLINE COEFFICIENTS
    !                ARE STORED IN THE ARRAY C.
    !         ISW=2  IF THE SPLINE DEFINED BY THE COEFFICIENT ARRAY 'C' IS
    !                TO BE EVALUATED (INTERPOLATED) AT THE POINT DEFINED BY
    !                THE PARAMETER 'XPT'.
    !         ISW=3  AS IN ISW=2, ONLY THE DERIVATIVE/3.D0 IS ALSO CALCULATE
    !         ISW=4  THE DERIVATIVE CALCULATED BY THE LAST USE OF SPLINE WIT
    !                IS RETURNED.
    !
    !     NN ... THE NUMBER OF KNOTS (DATA POINTS) TO WHICH THE SPLINE IS TO
    !            BE FITTED
    !
    !     X,Y ... THE ARRAYS DEFINING THE KNOTS.  THE X-VALUES MUST BE IN
    !             INCREASING ORDER.  THE ARRAYS MUST BE DIMENSIONED AT LEAST
    !             NN.
    !
    !     C ... THE ARRAY THAT CONTAINS THE CUBIC SPLINE COEFFICIENTS.
    !           MUST BE DIMENSIONED AT LEAST NN+2 .
    !
    !     D ... A WORK SPACE.  MUST BE DIMENSIONED AT LEAST NN+2 .
    !
    !     XPT ... THE POINT AT WHICH THE INTERPOLATION IS DESIRED (IF ISW IS
    !              SET TO 2).  THE VALUE OF SPLINE IS SET TO THE
    !              INTERPOLATED VALUE.
    !
    !
    ! *****  USER NOTES  *****
    !
    !     INTERPOLATION INVOLVES AT LEAST TWO STEPS .......
    !
    !       A.  CALL SPLINE WITH THE KNOTS.  THIS SETS UP THE
    !           COEFFICIENT ARRAY C.
    !           EG.  DUMY=SPLINE(1,NN,X,Y,C,D,XPT)
    !
    !       B.  CALL SPLINE WITH THE ARRAY C WHICH WAS DEFINED BY THE
    !           PREVIOUS CALL AND WILL BE USED TO FIND THE VALUE AT THE
    !           POINT 'XPT' .
    !           EG.   VALUE=SPLINE(2,NN,X,Y,C,D,XPT)
    !
    !     STEP 'A' NEED BE EXECUTED ONLY ONCE FOR A GIVEN SET OF KNOTS.
    !     STEP B MAY BE EXECUTED AS MANY TIMES AS NECESSARY.
    !
    !     Output from this file is written to UNIT NFLAG(18) which is
    !     defined in a parameter statement.
    !
    n   = nn
    np1 = n+1
    np2 = n+2
    z   = xpt
    if (isw .eq. 1) then
        c(1)   = y(1)
        d(1)   = 1.0d0
        c(np1) = 0.0d0
        d(np1) = 0.0d0
        c(np2) = 0.0d0
        d(np2) = 0.0d0
        do 41 i = 2, n
            c(i) = y(i)-y(1)
        41 d(i) = x(i)-x(1)
        do 410 i = 3, np2
            if(d(i-1).eq.0.0d0) then
                write(nflag(18), '(1X,5X,'' ***** ERROR IN SPLINE ... DIVIDE FAULT *****'')')
                stop 'SPLINE 1'
            end if
            pivot = 1.0d0/d(i-1)
            if(i.lt.np2) then
                supd = x(i-1)-x(i-2)
                if(supd.lt.0.0d0) then
                    write(nflag(18), '(1X,5X,''***** ERROR IN SPLINE ... UNORDERED X-VALUES *****'')')
                    stop 'SPLINE 2'
                end if
                supd = supd*supd*supd
            else
                supd = 1.0d0
            end if
            dfact = supd*pivot
            cfact = c(i-1)*pivot
            if(i.le.n) then
                do 47 j = i, n
                    v    = x(j)-x(i-2)
                    c(j) = c(j)-d(j)*cfact
                47 d(j) = v*v*v-d(j)*dfact
            end if
            if(i.lt.np2) then
                c(np1) = c(np1)-d(np1)*cfact
                d(np1) = 1.0d0-d(np1)*dfact
            end if
            c(np2) = c(np2)-d(np2)*cfact
            d(np2) = x(i-2)-d(np2)*dfact
        410 continue
        do 411 i = 1, n
            j = np2-i
            if(j.eq.np1) then
                v = 1.0d0
            else
                v = x(j)-x(j-1)
                v = v*v*v
            end if
            if(d(j+1).eq.0.0d0) then
                write(nflag(18), '(1X,5X,'' ***** ERROR IN SPLINE ... DIVIDE FAULT *****'')')
                stop 'SPLINE 3'
            end if
            c(j+1) = c(j+1)/d(j+1)
            c(j)   = c(j)-c(j+1)*v
        411 continue
        if(d(2).eq.0.0d0) then
            write(nflag(18), '(1X,5X,'' ***** ERROR IN SPLINE ... DIVIDE FAULT *****'')')
            stop 'SPLINE 4'
        end if
        c(2) = c(2)/d(2)
    else if (isw .eq. 2) then
        spline = c(1) + c(2)*(z-x(1))
        do 51 i = 1, n
            v = z-x(i)
            if(v.le.0) return
            spline = spline + c(i+2)*v**3
        51 continue
    else if (isw .eq. 3) then
        spline = c(1) + c(2)*(z-x(1))
        deriv  = c(2)/3.d0
        do 53 i = 1, n
            v = z-x(i)
            if(v.le.0) return
            v2     = v*v
            spline = spline + c(i+2)*v2*v
            deriv  = deriv + c(i+2)*v2
        53 continue
    else if (isw .eq. 4) then
        spline = 3.d0*deriv
    end if
end function spline

!
block data ptpacm
    implicit none
    double precision :: a, a1dif, a1inf, a2dif, a2inf, alph1, alph2, anuzero, b1dif, b1inf, b2dif, b2inf, beta, bsp, &
        & cart, cartnu, csp, d, d1dif, d1inf, d2dif, d2inf, dasy, de, easyab, easyac, easybc, ezero, gamp, phi, r1dif, &
        & r1inf, r1s, r2dif, r2inf, r2s, rbb, req, rhh, vsp
    integer :: i, iayes, icartr, indexes, irctnt, mder, msurf, nasurf, natoms, nder, nflag, nphi, nulbl
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
    common /ambpcm/ d(3), beta, req(3), a, gamp(3), rbb, rhh, de
    common /pre2cm/ d1inf, r1inf, b1inf, d2inf, &
        & r2inf, b2inf, a1inf, a2inf, &
        & r1s, r2s, vsp, alph1, alph2, &
        & phi(20), bsp(20, 3), csp(22, 3), &
        & d1dif, r1dif, b1dif, a1dif, &
        & d2dif, r2dif, b2dif, a2dif, &
        & dasy, iayes, nphi
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
    !--- bend potential parameters for A'' surface of O + H2
    !   Assign the parameters for the non-collinear part of the potential
    !
    !   Dissociation energy in kcal/mol
    !
    data d /106.56d0, 109.472d0, 106.56d0/
    !
    !   Morse beta in reciprocal Angstroms
    !
    data beta /2.07942d0/
    !
    !   Equilibrium bond lengths in Angstroms
    !
    data req /0.96966d0, 0.74144d0, 0.96966d0/
    !
    !   Pauling parameters in Angstroms
    !
    data a /0.26d0/
    !
    !   Gamma parameters for the bend correction
    !
    data gamp /0.210820d0, 5.42027d0, -1.44235d0/
    !
    !   RBB in Angstroms - this number is used in
    !   a scaling calculation for the bending correction
    !
    data rbb /0.74144d0/
    !
    data d1inf /  95.3d0    /
    data r1inf /    .76d0   /
    data b1inf /   2.0202d0 /
    data d2inf /  92.4d0    /
    data r2inf /    .98d0   /
    data b2inf /   2.41105d0/
    data a1inf /   0.0d0    /
    data a2inf /   0.0d0    /
    data r1s   /  2.21208d0 /
    data r2s   /  2.21531d0 /
    data vsp   / 95.3d0     /
    data iayes /  0         /
    data alph1 /  1.6d0     /
    data alph2 /  1.6d0     /
    data nphi / 11 /
    data (phi(i), i=1, 11) &
        & / 0.0000000d+00, 0.1000000d+02, &
        & 0.2000000d+02, 0.3000000d+02, &
        & 0.3733000d+02, 0.4000000d+02, &
        & 0.5000000d+02, 0.6000000d+02, &
        & 0.7000000d+02, 0.8000000d+02, &
        & 0.9000000d+02/
    data (bsp(i, 1), i=1, 11) &
        & / 0.9477920d+02, 0.9316200d+02, &
        & 0.9035320d+02, 0.8517960d+02, &
        & 0.8272720d+02, 0.8304320d+02, &
        & 0.8855580d+02, 0.9003100d+02, &
        & 0.9129800d+02, 0.9184900d+02, &
        & 0.9209100d+02/
    data (bsp(i, 2), i=1, 11) &
        & / 0.1456500d+01, 0.1475500d+01, &
        & 0.1530500d+01, 0.1606400d+01, &
        & 0.1634900d+01, 0.1632300d+01, &
        & 0.1534600d+01, 0.1396600d+01, &
        & 0.1299800d+01, 0.1243800d+01, &
        & 0.1225600d+01/
    data (bsp(i, 3), i=1, 11) &
        & / 0.2041900d+01, 0.1987500d+01, &
        & 0.1905400d+01, 0.1672700d+01, &
        & 0.1554000d+01, 0.1539400d+01, &
        & 0.1749300d+01, 0.2014500d+01, &
        & 0.2227100d+01, 0.2348800d+01, &
        & 0.2383800d+01/
    data (csp(i, 1), i=1, 13) &
        & / 0.9477920d+02, -0.1478277d+00, &
        & -0.1389235d-03, -0.3580591d-03, &
        & 0.1450636d-02, 0.1148717d-02, &
        & -0.9458414d-03, -0.4983968d-02, &
        & 0.5921762d-02, -0.2814975d-02, &
        & 0.1001140d-02, -0.2747850d-03, &
        & -0.5702495d-05/
    data (csp(i, 2), i=1, 13) &
        & / 0.1456500d+01, 0.1137130d-02, &
        & 0.7628700d-05, -0.9772201d-05, &
        & -0.1201120d-04, 0.2800151d-05, &
        & -0.4114549d-04, 0.7866376d-04, &
        & -0.1122307d-04, -0.1936703d-04, &
        & 0.6791209d-05, -0.1039780d-04, &
        & 0.8032967d-05/
    data (csp(i, 3), i=1, 13) &
        & / 0.2041900d+01, -0.5887558d-02, &
        & 0.4475579d-05, -0.5455348d-04, &
        & 0.1230139d-03, -0.6997501d-04, &
        & 0.2822753d-03, -0.4019424d-03, &
        & 0.1215636d-03, -0.1548133d-04, &
        & 0.9961687d-05, 0.1813458d-04, &
        & -0.1747243d-04/
end

!
!*****

!> @brief Evaluate POTINFO for the OH2_App_GEN_M21986 potential.
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

!> @brief Evaluate ANCVRT for the OH2_App_GEN_M21986 potential.
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

!> @brief Evaluate CARTOU for the OH2_App_GEN_M21986 potential.
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

!> @brief Evaluate CARTTOR for the OH2_App_GEN_M21986 potential.
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

!> @brief Evaluate EUNITZERO for the OH2_App_GEN_M21986 potential.
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

!> @brief Evaluate RTOCART for the OH2_App_GEN_M21986 potential.
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

!> @brief Evaluate DEDCOU for the OH2_App_GEN_M21986 potential.
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
