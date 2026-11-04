!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the OH2_GEN_DIM1981 potential.
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
    v      = 0.d0; g = 0.d0; d = 0.d0
    nder   = igrad
    nasurf = 1
    cart   = 0.d0
    do iatom=1, natoms
        do idir=1, 3
            cart(iatom, idir) = x(iatom,idir)/legacy_bohr_to_angstrom
        enddo
    enddo
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
!> @brief Evaluate PREPOT for the OH2_GEN_DIM1981 potential.
subroutine prepot
    !
    !
    !   System:          OH2
    !   Functional form: Anti-Morse Bends plus Bowman's collinear spline
    !   Common name:     OH2DIM
    !   Interface:       potlib2001
    !   Reference:
    !
    !   Number of bodies:              3
    !   Number of derivatives:         1
    !   Number of electronic surfaces: 1
    !
    !   PREPOT must be called once before any calls to POT.
    !   The control flags for the potential, such as IPRT, are
    !   initialized in the block data subprogram PTPACM.
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
    !   is in Hartree atomic units.
    !
    !   This potential is written such that:
    !                  R(1) = R(O-H)
    !                  R(2) = R(H-H)
    !                  R(3) = R(H-O)
    !   The zero of energy is defined at O "infinitely" far from the H2 diat
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
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !          GENERAL INFORMATION
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !   POTA.FOR HAS ANTI-MORSE BEND PLUS BOWMAN'S COLLINEAR SPLINE ROUTINE.
    !   ANTI-MORSE BENDING POTENTIAL IS ADDED TO COLLINEAR POTENTIAL.
    !   THIS SUBROUTINE CAN BE USED TO ADD A BENDING CORRECTION TO
    !   ANY COLLINEAR POTENTIAL SUBROUTINE.
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCC
    !          DEFINE VARIABLES
    !CCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    implicit none
    double precision :: a, aab, abc, anuzero, beta, bsp, cab, cart, cartnu, cbc, dbcoor, de, deesdr, degsdr, deijdr, &
        & descart, dgscart, dijcart, easyab, easyac, easybc, engyes, engygs, engyij, ezero, gamma, pengyes, pengygs, &
        & pengyij, phi, r, r1s, r2s, rbb, req, rhh, step, vsp, xab, xbc
    integer :: i, icartr, indexes, irctnt, it, mder, msurf, nasurf, &
        & natoms, nder, nflag, nphi, nulbl
    !
    double precision :: nab, nab1, nbc, nbc1, junk
    dimension junk(4)
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
    ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
    !
    double precision, parameter :: bohr = .52917706d0
    !
    ! Conversion constant from Hartree to kcal (kcal/Hartree)
    !
    double precision, parameter :: hrtree = 627.5095d0
    !
    common /nstep/  step
    common /potcm/  xbc(3), xab(3), cbc(3), abc(3), cab(3), &
        & aab(3), phi(40), bsp(40, 3), &
        & de, beta, req(3), a, gamma, rbb, rhh, &
        & r1s, r2s, vsp, dbcoor(2), nphi
    !
    if( natoms .gt. 25 )then
        write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    endif
    !
    !      WRITE(NFLAG(18), 600)
    !
    call prepot2
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !      BENDING POTENTIAL PARAMETERS
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !      WRITE(NFLAG(18),3)
    !
    !   Echo the potential parameters to the file linked to UNIT NFLAG(18)
    !
    !      WRITE (NFLAG(18),25) (REQ(IT), IT = 1,3)
    !      WRITE(NFLAG(18),10) DE
    !      WRITE(NFLAG(18),21) BETA
    !      WRITE (NFLAG(18),30) A
    !      WRITE (NFLAG(18),31) GAMMA
    !      WRITE(NFLAG(18),35) RBB
    !
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !       CONVERT TO ATOMIC UNITS
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    do 50 it = 1, 3
        req(it) = req(it)/bohr
    50 continue
    beta = beta * bohr
    de   = de/hrtree
    a    = a/bohr
    rhh  = 1.4007977d0
    rbb  = rbb /bohr
    !
    !CCCCCCCCCCC
    !
    !
    ezero(1) = vsp
    !
    do i=1, 5
        ref(i) = ' '
    enddo
    !
    ref(1) = 'No Reference'
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
!CCCCCCCCCCC
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!        MOVE TO POT SUBROUTINE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!> @brief Evaluate POT for the OH2_GEN_DIM1981 potential.
subroutine pot
    !
    !   The potential energy in the AB valley, EASYAB, is equal to the poten
    !   energy of the H "infinitely" far from the OH diatomic, with the
    !   OH diatomic at its equilibrium configuration.  Similarly, the terms
    !   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
    !   respectively.
    !
    !   This potential is written such that:
    !                  R(1) = R(O-H)
    !                  R(2) = R(H-H)
    !                  R(3) = R(H-O)
    !   The zero of energy is defined at O "infinitely" far from the H2 diat
    !
    !CCCCCCCCCCCCCCC
    !      ENTRY POT
    !CCCCCCCCCCCCCCC
    !
    implicit none
    double precision :: a, aab, abc, anuzero, bcoor, beta, bsp, c, cab, cart, cartnu, cbc, dac, dbcoor, de, deesdr, &
        & degsdr, deijdr, descart, dgscart, dijcart, easyab, easyac, easybc, engyes, engygs, engyij, ezero, gamma, &
        & pengyes, pengygs, pengyij, phi, r, r1s, r2s, ractmp, rbb, req, rhh, rob, step, stuff, vabbc, vac, vsp, xab, &
        & xbc
    integer :: icartr, indexes, irctnt, it, mder, msurf, nasurf, natoms, &
        & nder, nflag, nphi, nulbl
    !
    double precision :: nab, nab1, nbc, nbc1, junk
    dimension junk(4)
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
    ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
    !
    double precision, parameter :: bohr = .52917706d0
    !
    ! Conversion constant from Hartree to kcal (kcal/Hartree)
    !
    double precision, parameter :: hrtree = 627.5095d0
    !
    common /nstep/  step
    common /potcm/  xbc(3), xab(3), cbc(3), abc(3), cab(3), &
        & aab(3), phi(40), bsp(40, 3), &
        & de, beta, req(3), a, gamma, rbb, rhh, &
        & r1s, r2s, vsp, dbcoor(2), nphi
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
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !        FIRST CALCULATE SOME USEFUL NUMBERS
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
    !        AB AND BC RESPECTIVELY. (EQUATION 2.2)
    !
    nab = exp((req(1) - r(1))/a)
    !
    if( nab .lt. 1.d-15 ) nab = 1.d-15
    !
    nbc = exp((req(2) - r(2))/a)
    !
    if( nbc .lt. 1.d-15 ) nbc = 1.d-15
    !
    !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
    !       SEE EQUATION 2.7A, 2.7B
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
    !          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
    !
    !          FIRST TRAP OUT ANY ZERO ARGUEMENTS
    !
    95 if (nab1*nbc1.gt.0.0d0) go to 103
    rob = 1.d0
    go to 107
    !
    103 stuff = a * log(nab1*nbc1)
    rob = (rbb - stuff)/(rhh - stuff)
    !
    107 continue
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !        CALCULATE BENDING CORRECTION
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !        CALCULATE V(R(3))
    !
    junk(3) = exp(-beta*(r(3)-req(3))/rob)
    vac     = gamma *1.0d-14 * de * junk(3)* (1.d0 + 0.5d0*junk(3))
    !
    !        CALCULATE V(R(1)+R(2))
    !
    junk(4) = exp(beta*(req(3)-r(2)-r(1))/rob)
    vabbc   = gamma *1.0d-14 * de * junk(4) * (1.d0 + 0.5d0*junk(4))
    !
    !        HERE IS THE BENDING CORRECTION
    !
    bcoor = (vac - vabbc)*1.0d14
    !
    !          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
    !          OF BCOOR WITH RESPECT TO R(1), R(2), AND R(3).
    !
    !        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
    !
    if( nder .eq. 1 )then
        dac       = -gamma*de*beta*junk(3)*(1.d0+junk(3))/rob
        dbcoor(1) = gamma * de * beta * junk(4) * (1.d0 + junk(4))/rob
        dbcoor(2) = dbcoor(1)
    endif
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !        NOW CALCULATE THE COLLINEAR ENERGY
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR
    !          GEOMETRY AND CALL POT2
    !
    ractmp = r(3)
    r(3)   = r(2) + r(1)
    !
    call pot2
    !
    if( nder .eq. 1 )then
        degsdr(1) = degsdr(1) + degsdr(3)
        degsdr(2) = degsdr(2) + degsdr(3)
        degsdr(3) = dac
    endif
    r(3) = ractmp
    !
    engygs = engygs + bcoor
    !
    if( nder .eq. 1 )then
        do 300 it = 1, 2
            degsdr(it) = dbcoor(it) + degsdr(it)
        300 continue
    endif
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
!> @brief Evaluate PREPOT2 for the OH2_GEN_DIM1981 potential.
subroutine prepot2
    !
    !   POTENTIAL FUNCTION FOR RATE1D SERIES
    !   RMO-SPLINE TYPE FIT WITH BEBO TYPE CONTINUATION INTO THE ASYMBTOTIC
    !   POTENTIAL ROUTINE INTERACTS WITH REST OF PROGRAM WITH UNITS OF EV/BO
    !   POTENTIAL ROUTINE ASSUMES SPLINE INFO IS IN UNITS OF KCAL/ANGSTROM
    !
    implicit none
    double precision :: a, aab, abc, anuzero, beta, bsp, cab, cart, cartnu, cbc, csp, dbcoor, de, easyab, easyac, easybc, &
        & ezero, gamma, phi, r1s, r2s, rbb, req, rhh, scrap, step, v, vsp, xab, xbc
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nphi, nulbl
    !
    double precision :: nab, nab1, nbc, nbc1, junk
    dimension junk(4)
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
    ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
    !
    double precision, parameter :: bohr = .52917706d0
    !
    ! Conversion constant from Hartree to kcal (kcal/Hartree)
    !
    double precision, parameter :: hrtree = 627.5095d0
    !
    double precision, parameter :: detrad = .0174532d0
    !
    common /nstep/  step
    common /potcm/  xbc(3), xab(3), cbc(3), abc(3), cab(3), &
        & aab(3), phi(40), bsp(40, 3), &
        & de, beta, req(3), a, gamma, rbb, rhh, &
        & r1s, r2s, vsp, dbcoor(2), nphi
    !
    dimension       v(4), csp(48,3), scrap(48)
    !
    character(len=80) :: title
    !
    title = ' RMOFIT to DIM surface with BEBO continuations '
    !
    !      WRITE(NFLAG(18),845)TITLE
    !      WRITE(NFLAG(18),125)STEP
    !      WRITE(NFLAG(18),111)XAB(1), XBC(1), XAB(1),
    !     *               XAB(2), XBC(2), XAB(2),
    !     *               XAB(3), XBC(3), XAB(3)
    !      WRITE(NFLAG(18),105)R1S,R2S,VSP
    !      WRITE(NFLAG(18),113)(PHI(I),(BSP(I,J),J = 1,3),I = 1,NPHI)
    !      WRITE(NFLAG(18),117)(CBC(I),ABC(I),I = 1,3),
    !     +                    (CAB(I),AAB(I),I = 1,3)
    !
    !   Initialize the energy in the asymptotic valleys
    !
    easyab = xab(1)/hrtree
    easybc = xbc(1)/hrtree
    easyac = easyab
    !
    !
end subroutine prepot2

!
!> @brief Evaluate POT2 for the OH2_GEN_DIM1981 potential.
subroutine pot2
    !
    ! COCK SPLINES
    !
    !      ENTRY POT2
    implicit none
    double precision :: a, aab, abc, angle, anuzero, bet, beta, betpr, betrg, bsp, cab, cart, cartnu, cbc, csp, dbc, &
        & dbcoor, de, deesdr, degsdr, deijdr, depr, derg, descart, dgscart, dijcart, dis, dr1, dr2, easyab, easyac, &
        & easybc, engyes, engygs, engyij, ezero, gamma, pengyes, pengygs, pengyij, phi, ppp, r, r1, r1s, r2, r2s, r3, &
        & rbb, repr, req, req2, rerg, rhh, rr, scrap, sep, spline, step, v, vsp, x, xab, xbc
    integer :: i, ic, icartr, indexes, irctnt, j, k, mder, msurf, nasurf, natoms, nder, nflag, nphi, nulbl
    !
    double precision :: nab, nab1, nbc, nbc1, junk, bohr
    dimension junk(4)
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
    double precision, parameter :: hrtree = 627.5095d0
    double precision, parameter :: detrad = .0174532d0
    !
    common /nstep/  step
    common /potcm/  xbc(3), xab(3), cbc(3), abc(3), cab(3), &
        & aab(3), phi(40), bsp(40, 3), &
        & de, beta, req(3), a, gamma, rbb, rhh, &
        & r1s, r2s, vsp, dbcoor(2), nphi
    !
    dimension v(4), csp(48,3), scrap(48)
    !
    ic = 1
    10 go to (81, 82, 83, 84), ic
    81 r1=r(1)
    r2 = r(2)
    r3 = r(3)
    go to 85
    82 r1=r(1) + step
    r2 = r(2)
    r3 = r(3)
    go to 85
    83 r1=r(1)
    r2 = r(2) + step
    r3 = r(3)
    go to 85
    84 r1=r(1)
    r2 = r(2)
    r3 = r(3) + step
    85 r1=r1*0.529177d0
    r2    = r2*0.529177d0
    r3    = r3*0.529177d0
    dr1   = r1s - r1
    dr2   = r2s - r2
    angle = atan(dr1/dr2)
    angle = angle/detrad
    do 500 k=1, 3
        ppp = spline(1, nphi, phi, bsp(1, k), csp(1,k), scrap, angle)
    500 continue
    !
    ! SETUP ANY CONSTANTS
    !
    dbc   = xbc(1)
    derg  = xbc(1)/23.061d0
    rerg  = (r2s-xbc(2))/.529177d0
    betrg = xbc(3)*.529177d0
    depr  = xab(1)/23.061d0
    repr  = (r1s-xab(2))/.529177d0
    betpr = xab(3)*.529177d0
    !
    ! ENTRY INTO POTENTIAL CALCULATOR
    !
    if( r1 .gt. r1s ) go to 100
    if( r2 .gt. r2s ) go to 102
    !
    ! INTERACTION REGION OF POLAR COORDINATES
    !
    dr1   = r1s-r1
    dr2   = r2s-r2
    angle = atan(dr1/dr2)
    angle = angle/detrad
    rr    = sqrt(dr1*dr1+dr2*dr2)
    dis   = spline(2, nphi, phi, bsp(1, 1), csp(1,1), scrap, angle)
    req2  = spline(2, nphi, phi, bsp(1, 2), csp(1,2), scrap, angle)
    bet   = spline(2, nphi, phi, bsp(1, 3), csp(1,3), scrap, angle)
    !      D(1)=SPLINE(3,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)
    !      D(2)=SPLINE(3,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)
    !30      D(1)=D(1)/627.5095D0*0.529177D0
    !      D(2)=D(2)/627.5095D0*0.529177D0
    !
    ! 30   X = DIS*((1.D0-EXP(BET*(RR-REQ2)))**2.0D0 -1.D0) + VSP
    30 x = dis*((1.d0-exp(bet*(rr-req2)))**2.0d0 -1.d0) +&
        &ezero(1) - (anuzero*627.5095d0)
    !       1
    v(ic) = x
    ic    = ic+1
    if( ic .lt. 5 ) go to 10
    do 333 i=1, 3
    333 v(i)=v(i)/627.5095d0
    if( nder .eq. 1 )then
        degsdr(1) = (v(2)-v(1))/step
        degsdr(2) = (v(3)-v(1))/step
        degsdr(3) = 0.0d0
    endif
    engygs = v(1)
    93 return
    !
    ! LARGE R1 ASYMPTOTIC REGION
    !
    100 if(r2.gt.r2s) go to 104
    rr   = r2s-r2
    sep  = r1-r1s
    dis  = xbc(1) + cbc(1)*exp(-min(abc(1)*sep, 25.d0))
    req2 = xbc(2) + cbc(2)*exp(-min(abc(2)*sep, 25.d0))
    bet  = xbc(3) + cbc(3)*exp(-min(abc(3)*sep, 25.d0))
    !      IF(R1.GT.10.D0)D(1)=1.0D1
    !      IF(R1.GT.10.D0)D(2)=1.0D-8
    !      IF(R1.GT.10.D0)GO TO 30
    !      D(1)=-CBC(1)*ABC(1)*((1.D0-EXP(-BET*(RR-REQ2)))**2-1.D0)-
    !     +     DIS*2.D0*(1.D0-EXP(-BET*(RR-REQ2)))*EXP(-BET*(RR-REQ2))*
    !     +     ((-CBC(3)*ABC(3)*EXP(-ABC(3)*SEP)*(RR-REQ2) -
    !     +     CBC(2)*ABC(2)*EXP(ABC(2)*SEP)*BET))*BET*(RR-REQ2)
    !      D(2)=2.D0*DIS*(1.D0-EXP(-BET*(RR-REQ2)))*BET*(RR-REQ2)*
    !     +     EXP(-BET*(RR-REQ2))*BET
    go to 30
    !
    ! LARGE R2 ASYMPTOTIC REGION
    !
    102 continue
    rr   = r1s-r1
    sep  = r2-r2s
    dis  = xab(1) + cab(1)*exp(-min(aab(1)*sep, 25.d0))
    req2 = xab(2) + cab(2)*exp(-min(aab(2)*sep, 25.d0))
    bet  = xab(3) + cab(3)*exp(-min(aab(3)*sep, 25.d0))
    !      IF(R2.GT.10.D0)D(2)=1.0D1
    !      IF(R2.GT.10.D0)D(1)=1.0D-8
    !      IF(R2.GT.10.D0)GO TO 30
    !      D(1)=-CAB(1)*AAB(1)*((1.D0-EXP(-BET*(RR-REQ2)))**2-1.D0)-
    !     +     DIS*2.D0*(1.D0-EXP(-BET*(RR-REQ2)))*EXP(-BET*(RR-REQ2))*
    !     +     ((-CAB(3)*AAB(3)*EXP(-AAC(3)*SEP)*(RR-REQ2) -
    !     +     CAB(2)*AAB(2)*EXP(AAB(2)*SEP)*BET))*BET*(RR-REQ2)
    !      D(2)=2.D0*DIS*(1.D0-EXP(-BET*(RR-REQ2)))*BET*(RR-REQ2)*
    !     +     EXP(-BET*(RR-REQ2))*BET
    go to 30
    !
    ! LARGE R1,R2 3 BODY BREAK-UP REGION
    !
    104 continue
    j      = j+1
    engygs = xbc(1)/627.5095d0
    if( nder .eq. 1 )then
        degsdr(1) = 1.0d-9
        degsdr(2) = 1.0d-9
        degsdr(3) = 0.0d0
    endif
    go to 93
end subroutine pot2

!
!> @brief Evaluate SPLINE for the OH2_GEN_DIM1981 potential.
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
    !      DIMENSION X(*),Y(*),C(*),D(*)
    dimension x(nn+2), y(nn+2), c(nn+2), d(nn+2)
    !
    !
    !  THIS IS A SUBROUTINE FOR FITTING DATA WITH A CUBIC SPLINE
    !  POLYNOMIAL AND EVALUATING THAT POLYNOMIAL AT A GIVEN POINT
    !  OR ITS DERIVATIVE AT A GIVEN POINT
    !
    !  CALLING SEQUENCE .......
    !
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
    !     Output from this subprogram is written to unit 6.
    !
    2 n=nn
    np1 = n+1
    np2 = n+2
    z   = xpt
    24 go to (4, 5, 6, 7), isw
    4 c(1)=y(1)
    d(1)   = 1.0d0
    c(np1) = 0.0d0
    d(np1) = 0.0d0
    c(np2) = 0.0d0
    d(np2) = 0.0d0
    do 41 i=2, n
        c(i) = y(i)-y(1)
    41 d(i)=x(i)-x(1)
    do 410 i=3, np2
        if( d(i-1) .ne. 0 ) go to 43
        !      WRITE(NFLAG(18),1001)
        stop 'SPLINE 1'
        43 pivot=1.0d0/d(i-1)
        if( i .ge. np2 ) go to 45
        supd = x(i-1)-x(i-2)
        if( supd .ge. 0 ) go to 44
        !      WRITE(NFLAG(18),1000)
        stop 'SPLINE 2'
        44 supd=supd*supd*supd
        go to 46
        45 supd=1.0d0
        46 dfact=supd*pivot
        cfact = c(i-1)*pivot
        if( i .gt. n ) go to 48
        do 47 j=i, n
            v    = x(j)-x(i-2)
            c(j) = c(j)-d(j)*cfact
        47 d(j)=v*v*v-d(j)*dfact
        48 continue
        if( i .ge. np2 ) go to 49
        c(np1) = c(np1)-d(np1)*cfact
        d(np1) = 1.0d0-d(np1)*dfact
        49 c(np2)=c(np2)-d(np2)*cfact
    410 d(np2)=x(i-2)-d(np2)*dfact
    do 411 i=1, n
        j = np2-i
        if( j .ne. np1 ) go to 413
        v = 1.0d0
        go to 414
        413 v=x(j)-x(j-1)
        v = v*v*v
        414 if(d(j+1).ne.0)go to 415
        !      WRITE(NFLAG(18),1001)
        stop 'SPLINE 3'
        415 c(j+1)=c(j+1)/d(j+1)
    411 c(j)=c(j)-c(j+1)*v
    if( d(2) .ne. 0 ) go to 416
    !      WRITE(NFLAG(18),1001)
    stop 'SPLINE 4'
    416 c(2)=c(2)/d(2)
    return
    5 spline=c(1)+c(2)*(z-x(1))
    do 51 i=1, n
        v = z-x(i)
        if( v .gt. 0 ) go to 51
        return
    51 spline=spline+c(i+2)*v*v*v
    return
    6 continue
    spline = c(1)+c(2)*(z-x(1))
    deriv  = c(2)/3.d0
    do 53 i = 1, n
        v = z-x(i)
        if( v .le. 0 ) return
        v2     = v*v
        spline = spline + c(i+2)*v2*v
        deriv  = deriv + c(i+2)*v2
    53 continue
    return
    7 continue
    spline = 3.d0*deriv
end function spline

!
!*****
!
block data ptpacm
    implicit none
    double precision :: a, aab, abc, anuzero, beta, bsp, cab, cart, cartnu, cbc, dbcoor, de, easyab, easyac, easybc, &
        & ezero, gamma, phi, r1s, r2s, rbb, req, rhh, step, vsp, xab, xbc
    integer :: i, icartr, indexes, irctnt, mder, msurf, nasurf, &
        & natoms, nder, nflag, nphi, nulbl
    !
    !         REAL*8 NAB,NAB1,NBC,NBC1,JUNK,HRTREE,BOHR
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
    common /nstep/  step
    common /potcm/  xbc(3), xab(3), cbc(3), abc(3), cab(3), &
        & aab(3), phi(40), bsp(40, 3), &
        & de, beta, req(3), a, gamma, rbb, rhh, &
        & r1s, r2s, vsp, dbcoor(2), nphi
    !
    !   Initialize the flags and the I/O unit numbers for the potential
    !
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
    data step /0.000000008d0/
    data nphi /17/
    data r1s /0.2212079d+01/
    data r2s /0.2215309d+01/
    data vsp /0.1094890d+03/
    data xbc / 0.1094890d+03, 0.1473409d+01, 0.1941998d+01/
    data xab / 0.1056400d+03, 0.1241479d+01, 0.2304798d+01/
    data cbc /-0.4815674d-01, 0.3814697d-05, 0.4415512d-03/
    data abc / 0.7873216d+01, 0.7873216d+01, 0.7873216d+01/
    data cab /-0.1482804d+01, -0.2084732d-02, -0.3727913d-02/
    data aab / 0.1932070d+01, 0.1932070d+01, 0.1932070d+01/
    data (phi(i), i=1, 17) &
        & / 0.0000000d+00, 0.2000000d+01, &
        & 0.4000000d+01, 0.6000000d+01, &
        & 0.1800000d+02, 0.2400000d+02, &
        & 0.3000000d+02, 0.3400000d+02, &
        & 0.3800000d+02, 0.4400000d+02, &
        & 0.4800000d+02, 0.5400000d+02, &
        & 0.6000000d+02, 0.7600000d+02, &
        & 0.8600000d+02, 0.8800000d+02, &
        & 0.9000000d+02/
    data (bsp(i, 1), i=1, 17) &
        & / 0.1094409d+03, 0.1094212d+03, &
        & 0.1094004d+03, 0.1093775d+03, &
        & 0.1088412d+03, 0.1077417d+03, &
        & 0.1055117d+03, 0.1030257d+03, &
        & 0.9966266d+02, 0.9630693d+02, &
        & 0.9637480d+02, 0.9811765d+02, &
        & 0.1000138d+03, 0.1029822d+03, &
        & 0.1038929d+03, 0.1040309d+03, &
        & 0.1041572d+03/
    data (bsp(i, 2), i=1, 17) &
        & / 0.1473413d+01, 0.1474257d+01, &
        & 0.1476847d+01, 0.1481213d+01, &
        & 0.1545280d+01, 0.1599182d+01, &
        & 0.1655956d+01, 0.1681611d+01, &
        & 0.1680401d+01, 0.1626349d+01, &
        & 0.1572269d+01, 0.1486160d+01, &
        & 0.1408802d+01, 0.1273705d+01, &
        & 0.1241724d+01, 0.1239825d+01, &
        & 0.1239394d+01/
    data (bsp(i, 3), i=1, 17) &
        & / 0.1942440d+01, 0.1940545d+01, &
        & 0.1936212d+01, 0.1929558d+01, &
        & 0.1825750d+01, 0.1719964d+01, &
        & 0.1603891d+01, 0.1538857d+01, &
        & 0.1490752d+01, 0.1525683d+01, &
        & 0.1620584d+01, 0.1788983d+01, &
        & 0.1942842d+01, 0.2222638d+01, &
        & 0.2294156d+01, 0.2299114d+01, &
        & 0.2301070d+01/
    !
    !       DISSOCIATION ENERGY IN KCAL/MOLE
    !
    data de /106.56d0/
    !
    !       MORSE BETAS--RECIPROCAL ANGSTROMS
    !
    data beta /2.07942d0/
    !
    !       THE EQUILIBRIUM DISTANCES IN ANGSTROMS
    !
    data req /0.96966d0, 0.74144d0, 0.96966d0/
    !
    !       PAULING PARAMETER IN ANGSTROMS
    !
    data a /0.26d0/
    !
    !       GAMMA FOR THE BEND CORRECTION
    !       GAMMA IS DIMENSIONLESS AND RANGES FROM .5 TO .65
    !
    data gamma / 0.4131d0/
    !
    !       RBB IN ANGSTROMS--THIS NUMBER IS USED IN A SCALING
    !       CALCULATION FOR THE BENDING CORRECTION
    !
    data rbb / 0.74144d0/
    !
end

!
!*****
!> @brief Evaluate POTINFO for the OH2_GEN_DIM1981 potential.
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

!> @brief Evaluate ANCVRT for the OH2_GEN_DIM1981 potential.
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

!> @brief Evaluate CARTOU for the OH2_GEN_DIM1981 potential.
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

!> @brief Evaluate CARTTOR for the OH2_GEN_DIM1981 potential.
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

!> @brief Evaluate EUNITZERO for the OH2_GEN_DIM1981 potential.
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

!> @brief Evaluate RTOCART for the OH2_GEN_DIM1981 potential.
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

!> @brief Evaluate DEDCOU for the OH2_GEN_DIM1981 potential.
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
