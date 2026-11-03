!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the FH2_LEPS_NO4_1984 potential.
subroutine pes( x, igrad, p, g, d )

    implicit none

    double precision :: anuzero, cart, descart, dgscart, dijcart, pengyes, pengygs, pengyij, v
    integer          :: iatom, idir, istate, nasurf, nder, nflag, nulbl
    ! number of electronic state
    integer, parameter            :: nstates = 1
    integer, parameter            :: natoms  = 3
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)

    integer, parameter :: natom = 25
    integer, parameter :: isurf = 5
    integer, parameter :: jsurf = int(isurf*(isurf+1)/2)

    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    logical, save :: first_time_data=.true.

    !initialize
    v = 0.d0
    g = 0.d0
    d = 0.d0

    cart = 0.d0
    do iatom=1,natoms
        do idir=1,3
            cart(iatom,idir) = x(iatom,idir)/0.529177211
        enddo
    enddo

    nder = igrad

    if(first_time_data) then
        call prepot
        first_time_data = .false.
    endif

    call pot

    if (igrad==0) then
        do istate=1,nstates
            p(istate) = pengygs*27.211386
        enddo
    else if (igrad==1) then
        do istate=1,nstates
            p(istate) = pengygs*27.211386
        enddo
        do iatom=1,natoms
            do idir=1,3
                g(1,iatom,idir) = dgscart(iatom,idir)*51.422067
            enddo
        enddo
    else if (igrad==2) then
        write (*,*) 'Only energy and gradient are available'
    endif

end subroutine pes
!
!> @brief Evaluate PREPOT for the FH2_LEPS_NO4_1984 potential.
subroutine prepot
    !
    !   System:          FH2
    !   Functional form: eLEPS (extended-London-Eyring-Polanyi-Sato) plus
    !                    a three-centered term
    !   Common name:     no. 4
    !   Reference:       F. B. Brown, R. Steckler, D. W. Schwenke,
    !                    D. G. Truhlar, and B. C. Garrett,
    !                    J. Chem. Phys. 82, 188-201 (1985).
    !
    !   Notes:           In this version of the FH2 potential the AB and AC
    !                    Sato parameters are functions of PHI, where PHI is
    !                    the angle of F-to-BC vector with the BC axis.
    !                    The BC Sato parameter is a function of R(BC), the
    !                    BC bond length, and PHI.
    !
    ! Number of bodies: 3
    ! Interface: potlib2001
    ! Number of electronic states: 1
    ! Number of derivatives: 1
    !
    !   PREPOT must be called once before any calls to POT.
    !   The potential parameters are included in DATA statements
    !   and in the block data subprogram PTPACM.
    !   Coordinates, potential energy, and derivatives are passed
    !   The potential energy in the three asymptotic valleys are
    !   stored in the common block ASYCM:
    !                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC
    !   The potential energy in the AB valley, EASYAB, is equal to the poten
    !   energy of the H "infinitely" far from the FH diatomic, with the
    !   FH diatomic at its equilibrium configuration.  Similarly, the terms
    !   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
    !   respectively.
    !   All the information passed through the common blocks PT1CM and ASYCM
    !   is in Hartree atomic units.
    !
    !   This potential is written such that:
    !                  R(1) = R(F-H)
    !                  R(2) = R(H-H)
    !                  R(3) = R(H-F)
    !   The zero of energy is defined at F "infinitely" far from the H2 diat
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
    implicit none

    double precision :: a, a1, a2, alf, alfp, alft, anuzero, b, bet, beta, betp, c, c1, c2, c3, c3c, c4, c5, c6, cart, &
        & cartnu, d, deg1, deg2, deg3, deg3p, deg4, dj, dm, dq, dz, easyab, easyac, easybc, ezero, f, g, g2, ha2, omp3c, &
        & omq3c, p3c, q, q3c, re, rr, xj, z, zfh0, zhh0, zhhdif, zhhst
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
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /asycm/ easyab,easybc,easyac
    !
    double precision, parameter :: ckcau  = 627.5095d0
    double precision, parameter :: cangau = 0.52917706d0
    double precision, parameter :: tdeg   = 57.29577951d0
    !
    common /satocm/ dm(3), re(3), beta(3)
    !
    !      LOGICAL LCOL, LDERIV, LZERO
    common /precm/ z(3),q(3),xj(3),dq(3,3),dj(3,3),dz(3,3), &
        & rr(3),a1,a2,a,b,c,d,f,g,bet,betp,zfh0, &
        & zhh0,zhhst,deg1,deg2,deg3,deg4,deg3p, &
        & c1,c2,c3,c4,c5,c6,c3c,alf,alfp,alft,p3c, &
        & q3c,zhhdif,omp3c,omq3c,ha2,g2
    !
    if(natoms.gt.25) then
        write(nflag(18),'(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    end if
    !
    !
    !
    do 10 i = 1,3
        dm(i)   = dm(i) / ckcau
        re(i)   = re(i) / cangau
        beta(i) = beta(i) * cangau
    10 continue
    !
    !   USEFUL CONSTANTS
    !
    zhhdif = zhhst - zhh0
    omp3c  = 1.d0 - p3c
    omq3c  = 1.d0 - q3c
    ha2    = 0.5d0 * a2
    g2     = g * g
    !
    !    Set the values of the classical energy in the three asymptotic vall
    !
    easyab = dm(1)
    easybc = dm(2)
    easyac = dm(3)
    !
    !
    ezero(1) = dm(2)
    !
    do i=1,5
        ref(i) = ' '
    end do
    !
    ref(1) = 'F. B. Brown, R. Steckler, D. W. Schwenke,'
    ref(2) = 'D. G. Truhlar, B. C. Garrett,'
    ref(3) = 'J. Chem. Phys. 82, 188(1985)'
    !
    indexes(1) = 9
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
!> @brief Evaluate POT for the FH2_LEPS_NO4_1984 potential.
subroutine pot
    !
    !   The potential energy in the AB valley, EASYAB, is equal to the poten
    !   energy of the H "infinitely" far from the FH diatomic, with the
    !   FH diatomic at its equilibrium configuration.  Similarly, the terms
    !   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
    !   respectively.
    !
    !
    !   This potential is written such that:
    !                  R(1) = R(F-H)
    !                  R(2) = R(H-H)
    !                  R(3) = R(H-F)
    !   The zero of energy is defined at F "infinitely" far from the H2 diat
    !
    !      ENTRY POT
    implicit none

    double precision :: a, a1, a2, alf, alfp, alft, anuzero, b, bet, beta, betp, c, c1, c2, c3, c3c, c4, c5, c6, capr, &
        & capr2, cart, cartnu, cosphi, costh, cphi, d, dcos1, dcos2, dcos3, dcphi, dd2, de3c1, de3c2, de3c3, deesdr, deg, &
        & deg1, deg2, deg3, deg3p, deg4, degdif, degsdr, deijdr, descart, dgscart, dijcart, dj, dm, dphi1, dphi2, dphi3, &
        & dq, dz, dzphi, e3c, easyab, easyac, easybc, engyes, engygs, engyij, ex, ex1, ex2, ex3, exx, ezero, f, g, g2, &
        & ha2, hsech, omp3c, omq3c, op3z, p3c, pengyes, pengygs, pengyij, phi, q, q3c, r, r12, r1i, r22, r32, r3i, rdif, &
        & rdif2, re, rr, rsum, scale1, scale2, sgnphi, sinp2, sinphi, t, t1, t2, t3, t4, t5, xj, xjx, z, zfh0, zhh, zhh0, &
        & zhhdif, zhhst, zp1, zp3, zphi
    integer :: i, i1, i2, i3, icartr, indexes, irctnt, j, mder, msurf, nasurf, natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter          :: n3atom = 75
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    !
    common /pt1cm/ r(n3atom),engygs,degsdr(n3atom)
    common /pt3cm/ ezero(isurf+1)
    common /pt4cm/ engyes(isurf),deesdr(n3atom,isurf)
    common /pt5cm/ engyij(jsurf),deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /asycm/ easyab,easybc,easyac
    !
    double precision, parameter :: ckcau  = 627.5095d0
    double precision, parameter :: cangau = 0.52917706d0
    double precision, parameter :: tdeg   = 57.29577951d0
    double precision, parameter :: exlarg = 80.d0
    !
    !   EXLARG IS LN(LARGEST FLOATING POINT NUMBER ALLOWED)
    !
    common /satocm/ dm(3), re(3), beta(3)
    !
    logical :: lcol, lderiv, lzero
    !
    common /precm/ z(3),q(3),xj(3),dq(3,3),dj(3,3),dz(3,3), &
        & rr(3),a1,a2,a,b,c,d,f,g,bet,betp,zfh0, &
        & zhh0,zhhst,deg1,deg2,deg3,deg4,deg3p, &
        & c1,c2,c3,c4,c5,c6,c3c,alf,alfp,alft,p3c, &
        & q3c,zhhdif,omp3c,omq3c,ha2,g2
    !
    call cartou
    call carttor
    !
    !   Check the values of NASURF and NDER for validity.
    !
    !   Set the logical flag LDERIV which is the variable that actually
    !   controls the calculation of the derivatives.
    !
    if (nder .eq. 1) then
        lderiv = .true.
    else
        lderiv = .false.
    endif
    !
    !   TEST R'S FOR PHYSICAL VALUES
    !
    do 20 i = 1,3
        rr(i) = r(i)
    20 continue
    do 30 i1 = 1,3
        i2 = mod(i1,3) + 1
        i3 = mod(i2,3) + 1
        t  = rr(i2) + rr(i3)
        if (rr(i1) .gt. t) rr(i1) = t
    30 continue
    r12 = rr(1) * rr(1)
    r22 = rr(2) * rr(2)
    r32 = rr(3) * rr(3)
    !
    !   CAPR IS THE F TO H2 DISTANCE
    !
    capr2 = 0.5d0 * (r12 + r32 - 0.5d0 * r22)
    if (capr2 .lt. 0.d0) capr2 = 0.0d0
    capr = sqrt(capr2)
    !
    !   PHI IS THE ANGLE BETWEEN CAPR AND RHH
    !   PHI IS RESTRICTED TO BE BETWEEN O AND 90 DEGREES
    !
    t = 0.5d0 * (r32 - r12)
    if (abs(t) .gt. 1.d-7) t = t / (rr(2) * capr)
    t2     = 1.0d0
    sgnphi = sign(t2,t)
    cosphi = sgnphi * t
    cosphi = min(cosphi,t2)
    phi    = acos(cosphi)
    sinp2  = 1.d0 - cosphi * cosphi
    if(sinp2 .lt. 0.d0) sinp2 = 0.d0
    sinphi = sqrt(sinp2)
    !
    !   LCOL IS SET TRUE IF COLLINEAR GEOMETRY
    !
    lcol = abs(sinphi) .lt. 1.d-6 .or. abs(capr) .lt. 1.d-6
    !
    !   COMPUTE COLLINEAR SATO PARAMETERS
    !
    z(1) = zfh0
    t    = a * (rr(2) - b)
    if (.not.(t .lt. exlarg)) go to 35
    t     = exp(t)
    t1    = 1.d0 / t
    hsech = 1.d0/(t + t1)
    zhh   = a1 + ha2 * ((t - t1) * hsech + 1.d0)
    go to 36
    35 continue
    hsech = 0.d0
    zhh   = a1 + a2
    36 continue
    t    = rr(2) - f
    exx  = c * exp(d * t * t)
    zhh  = zhh - exx
    z(2) = zhh
    if (.not.lderiv) go to 50
    do 40 i = 1,3
        do 40 j = 1,3
            dz(i,j) = 0.d0
    40 continue
    !
    !  ONLY THE HH SATO PARAMETER HAS A NONZERO DERIVATIVE - W.R.T. R2
    !
    dz(2,2) = 2.d0 * a * a2 * hsech * hsech - 2.d0 * d * t * exx
    50 continue
    !
    !   COMPUTE NONCOLLINEAR CONTRIBUTIONS TO SATO PARAMETERS
    !
    if (lcol) go to 160
    !
    !   COMPUTE Z(PHI) FOR FH SATO PARAMETER
    !   DZPHI IS THE DERIVATIVE OF Z(PHI) W.R.T. PHI
    !
    deg = tdeg * phi
    if (.not.(deg.lt.0. .or. deg.gt.90.d0)) go to 60
    write(nflag(18),'(/,2X,T5,''In the FH2 surface no. 4, PHI = '', F10.5, ''degrees'', /,2X,T5,''but PHI must be gre' // &
        'ater than zero.'')') deg
    stop 'POT 1'
    60 continue
    if (.not.(deg .le. deg1)) go to 70
    zphi  = 0.d0
    dzphi = 0.d0
    go to 140
    70 continue
    if (.not.(deg .le. deg2)) go to 80
    degdif = deg - deg1
    dd2    = degdif * degdif
    zphi   = dd2 * (c1 + c2 * degdif)
    if (lderiv) dzphi = tdeg * &
        & degdif * (2.d0 * c1 + 3.d0 * c2 * degdif)
    go to 130
    80 continue
    if (.not.(deg .le. deg3)) go to 90
    degdif = deg - deg3p
    zphi   = c3 + c4 * degdif
    if (lderiv) dzphi = c4 * tdeg
    go to 120
    90 continue
    if (.not.(deg .le. deg4)) go to 100
    zphi = c5 + c6 * sinphi * sinphi
    if (lderiv) dzphi = 2.d0 * c6 * sinphi * cosphi
    go to 110
    100 continue
    write (nflag(18),'(/,2X,T5,''In the FH2 surface no. 4, PHI = '', F10.5, ''degrees'', /,2X,T5,''but PHI must be le' // &
        'ss than 90 degrees.'')') deg
    stop 'POT 2'
    110 continue
    120 continue
    130 continue
    140 continue
    !
    !   COMPUTE SCALE FACTOR TO DAMP OUT CORRECTION AT LARGE DISTANCE
    !
    scale1 = capr2 / (g2 + capr2)
    z(1)   = z(1) + scale1 * zphi
    !
    !   COMPUTE SCALE FACTOR FOR HH SATO PARAMETER
    !
    scale2 = (zhh - zhh0) / zhhdif
    !
    !   COMPUTE PHI DEPENDENT PART OF HH SATO PARAMETER
    !
    cphi = (bet + betp * sinp2) * sinp2
    !
    !   COMPUTE HH SATO PARAMTER
    !
    z(2) = z(2) + scale1 * scale2 * cphi
    if (.not.lderiv) go to 150
    !
    !   COMPUTE DERIVATIVES
    !   FIRST, DERIVATIVES OF PHI W.R.T. R1, R2, R3
    !
    t     = 1.d0 / capr
    t1    = 0.5d0 * cosphi * t
    t2    = sgnphi / rr(2)
    dphi1 = rr(1) * (t1 + t2)
    dphi3 = rr(3) * (t1 - t2)
    t1    = 1.d0 / sinphi
    t2    = t1 * t
    dphi1 = dphi1 * t2
    dphi3 = dphi3 * t2
    dphi2 = rr(2) * cosphi * t1 * (1.d0 / r22 - 0.25d0 / capr2)
    t1    = (1.d0 - scale1) / capr2
    t2    = t1 * zphi
    !
    !   NEXT, DERIVATIVES OF FH SATO PARAMETERS
    !
    dz(1,1) = scale1 * (t2 * rr(1) + dzphi * dphi1)
    dz(1,2) = scale1 * (-0.5d0 * t2 * rr(2) + dzphi * dphi2)
    dz(1,3) = scale1 * (t2 * rr(3) + dzphi * dphi3)
    t1      = t1 * cphi
    !
    !   DERIVATIVE OF THE PHI DEPENDENT PART OF HH SATO W.R.T. PHI
    !
    dcphi = 2.d0 * (bet + 2.d0 * betp * sinp2) * sinphi * cosphi
    t2    = scale1 * scale2
    !
    !   LAST, DERIVATIVES OF HH SATO PARAMETER
    !
    dz(2,1) = t2 * (t1 * rr(1) + dcphi * dphi1)
    dz(2,2) = dz(2,2) * (1.d0 + scale1 * cphi / zhhdif) + t2 * (-0.5d0 * t1 * rr(2) + dcphi * dphi2)
    dz(2,3) = t2 * (t1 * rr(3) + dcphi * dphi3)
    150 continue
    160 continue
    !
    !   COMPUTE POTENTIAL FROM LEPS FORM
    !
    engygs = 0.d0
    z(3)   = z(1)
    if (.not.lderiv) go to 170
    do 165 i = 1,3
        dz(3,i)   = dz(1,i)
        degsdr(i) = 0.d0
        q(i)      = 0.d0
        xj(i)     = 0.d0
        do 162 j = 1,3
            dq(i,j) = 0.d0
            dj(i,j) = 0.d0
        162 continue
    165 continue
    170 continue
    lzero = .true.
    do 200 i = 1,3
        ex    = exp(beta(i) * (re(i) - rr(i)))
        lzero = lzero .and. ex.eq.0.0d0
        if (ex.eq.0.0d0) go to 195
        zp1  = z(i) + 1.d0
        zp3  = zp1 + 2.d0
        op3z = 1.d0 + 3.d0 * z(i)
        t    = 0.25d0 * dm(i) * ex / zp1
        t1   = zp3 * ex
        t2   = op3z * ex
        !
        !   Q AND XJ ARE THE COULUMB AND EXCHANGE PARTS, INDEX I RUNS OVER THE
        !   COORDINATES OF WHICH THEY ARE EXPLICIT FUNCTIONS
        !
        q(i)  = t * (t1 - 2.d0  * op3z)
        xj(i) = t * (t2 - 2.d0  * zp3)
        !
        !   PUT SUM OF Q'S IN ENERGY
        !
        engygs = engygs + q(i)
        if (.not.lderiv) go to 190
        !
        !   DERIVATIVE OF Q AND ZJ W.R.T. R1, R2, R3
        !   ELEMENT I,J IS THE DERIVATIVE OF Q(I) W.R.T. RJ
        !
        t3 = 2.d0 * t * (ex + 2.d0) / zp1
        do 180 j = 1,3
            dq(i,j) = -t3 * dz(i,j)
            dj(i,j) = -dq(i,j)
        180 continue
        t       = 2.d0*t*beta(i)
        dq(i,i) = dq(i,i) - t * (t1 - op3z)
        dj(i,i) = dj(i,i) - t * (t2 - zp3)
        190 continue
        195 continue
    200 continue
    if (lzero) go to 270
    xjx = 0.d0
    !
    !   COMPUTE THE TOTAL EXCHANGE TERM
    !
    do 230 i1 = 1,3
        i2  = mod(i1,3) + 1
        t   = xj(i1) - xj(i2)
        xjx = xjx + t * t
        if (.not.lderiv) go to 220
        !
        !   PUT THE DERIVATIVES OF THE TOTAL EXCHANGE TERM IN DV
        !
        do 210 j = 1,3
            degsdr(j) = degsdr(j) + t * (dj(i1,j) - dj(i2,j))
        210 continue
        220 continue
    230 continue
    xjx = sqrt(0.5d0 * xjx)
    !
    !   COMPUTE THE LEPS PART OF THE POTENTIAL ENERGY
    !
    engygs = engygs - xjx
    if (.not.lderiv) go to 260
    !
    !   COMPUTE DERIVATIVES OF THE LEPS PART
    !
    if(abs(xjx) .gt. 1.d-14) t = 0.5d0 / xjx
    do 250 j = 1,3
        !
        !   THE DERIVATIVE OF THE EXCHANGE PART MUST BE DIVIDED BY 2*XJX
        !
        degsdr(j) = -t * degsdr(j)
        do 240  i = 1,3
            !
            !   THEN ADD IN THE DERIVATIVE OF Q
            !
            degsdr(j) = degsdr(j) + dq(i,j)
        240 continue
    250 continue
    260 continue
    270 continue
    !
    !   EZERO ADDED TO PUT ZERO AT REQ FOR REACTANTS
    !
    engygs = engygs + ezero(1)
    !
    !   COMPUTE THE 3-CENTER TERM
    !
    rsum  = rr(1) + rr(3)
    rdif  = rr(1) - rr(3)
    rdif2 = rdif * rdif
    t1    = -alfp * rsum
    ex1   = c3c * exp(t1 * rsum)
    t2    = -alf * rdif
    ex2   = omp3c * exp(t2 * rdif)
    t3    = -alft * rdif2
    ex3   = p3c * exp(t3 * rdif2)
    r1i   = 1.d0 / rr(1)
    r3i   = 1.d0 / rr(3)
    t4    = r1i*r3i
    !
    !   COSTH IS THE ANGLE BETWEEN R1 AND R3
    !
    costh = 0.5d0 * (r12 + r32 - r22) * t4
    t5    = omq3c * costh
    t     = q3c + t5 * costh
    !
    !   E3C CONTAINS THE THREE CENTER CORRECTION TO THE ENERGY
    !
    e3c    = (ex2 + ex3) * t * ex1
    engygs = engygs + e3c
    if (.not.lderiv) go to 280
    !
    !   COMPUTE DERIVATIVE OF THE 3C TERM
    !   FIRST, DERIVATIVE OF COSTH
    !
    dcos1 = r3i - costh*r1i
    dcos2 = - rr(2)*t4
    dcos3 = r1i - costh*r3i
    t4    = ex2 + ex3
    t5    = t4 * t5
    !
    !   NOW, DERIVATIVES OF E3C
    !
    de3c1     = 2.d0*((t2 * ex2 + 2.d0 * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos1) * ex1
    de3c2     = 2.d0 * t5 * dcos2 * ex1
    de3c3     = 2.d0*((-t2 * ex2 - 2.d0 * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos3) * ex1
    degsdr(1) = degsdr(1) + de3c1
    degsdr(2) = degsdr(2) + de3c2
    degsdr(3) = degsdr(3) + de3c3
    280 continue
    !
    !
    call eunitzero
    if(nder.ne.0) then
        call rtocart
        call dedcou
    endif
    !
    !
end subroutine pot
!
!*****
!
block data ptpacm
    implicit none

    double precision :: a, a1, a2, alf, alfp, alft, anuzero, b, bet, beta, betp, c, c1, c2, c3, c3c, c4, c5, c6, cart, &
        & cartnu, d, deg1, deg2, deg3, deg3p, deg4, dj, dm, dq, dz, easyab, easyac, easybc, ezero, f, g, g2, ha2, omp3c, &
        & omq3c, p3c, q, q3c, re, rr, xj, z, zfh0, zhh0, zhhdif, zhhst
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
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /asycm/ easyab,easybc,easyac
    !
    common /satocm/ dm(3), re(3), beta(3)
    common /precm/ z(3),q(3),xj(3),dq(3,3),dj(3,3),dz(3,3), &
        & rr(3),a1,a2,a,b,c,d,f,g,bet,betp,zfh0, &
        & zhh0,zhhst,deg1,deg2,deg3,deg4,deg3p, &
        & c1,c2,c3,c4,c5,c6,c3c,alf,alfp,alft,p3c, &
        & q3c,zhhdif,omp3c,omq3c,ha2,g2
    !
    !   Initialize the flags and the I/O unit numbers for the potential
    !
    data nasurf /1,35*0/
    data nder /0/
    data nflag /1,1,15*0,6,0,0/
    !
    data anuzero /0.0d0/
    data icartr,msurf,mder/3,0,1/
    data nulbl /25*0/
    data natoms /3/
    !
    !   Initialize the potential parameters; the energy parameters are in
    !   kcal/mol, and the lengths are in Angstroms.
    !
    data dm /141.196d0, 109.449d0, 141.196d0/
    data re /0.9170d0, 0.7419d0, 0.9170d0/
    data beta /2.2187d0, 1.9420d0, 2.2187d0/
    !
    data a1,a2,a,b,c,d,f,g &
        & / 0.094d0,0.116d0,0.87d0,2.28d0,0.01056d0, &
        & -13.6221d0,2.70d0,0.2d0/
    data bet,betp,zfh0,zhh0,zhhst &
        & /-0.09627d0,-0.18075d0,0.175d0, &
        & 0.1146d0,0.167065d0/
    data deg1,deg2,deg3,deg4,deg3p &
        & /10.0d0,20.0d0,60.0d0,90.0d0,41.0d0/
    data c1,c2,c3,c4,c5,c6 &
        & /1.15d-4,-6.0d-6,0.016d0,0.0005d0,0.00069d0, &
        & 3.308d-2/
    data c3c,alf,alfp,alft,p3c,q3c &
        & /0.268340d0,0.68d0,0.078d0,2.5d0, &
        & 0.85d0,0.64d0/
    !
end
!
!*****




!> @brief Evaluate POTINFO for the FH2_LEPS_NO4_1984 potential.
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
    common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
    common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte, &
        & cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    write(nflag(18),'(/)')
    do i =1,5
        write(nflag(18),'(2x,a75)') ref(i)
    end do
    write(nflag(18),'(/)')
    kmax = 0
    do i = 1,isurf+1
        do j = 1,isurf+1
            if(nasurf(i,j).ne.0.and.kmax.lt.max(i,j)) kmax = max(i,j)
        enddo
    enddo
    write(nflag(18),'(2x,'' MAX. AND ACTUAL NO. OF EXCITED SURFACES: '',I3,5x,I3)') msurf,kmax-1
    if(kmax-1.gt.msurf) then
        write(6,*) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
        stop
    endif
    ksdiag = 0
    kediag = 0
    do i = 2,isurf+1
        if(nasurf(i,i).ne.0) then
            kediag = i-1
            if(ksdiag.eq.0) ksdiag = i-1
        endif
    enddo
    ksoffd = 0
    keoffd = 0
    k      = 0
    do i = 1,isurf
        do j = i+1,isurf+1
            k = k+1
            if(nasurf(i,j)+nasurf(j,i).ne.0) then
                keoffd = k
                if(ksoffd.eq.0) ksoffd = k
            endif
        enddo
    enddo
    write(nflag(18),'(2x,'' MAX. AND ACTUAL ORDER OF DERIVATIVES:    '',I3,5x,I3)') mder,nder
    if(nder.gt.mder) then
        write(6,*) ' WRONG INPUT ON ORDER OF DERIVATIVES'
        stop
    endif
    if(nflag(19).eq.1) then
        write(nflag(18),'(/)')
        write(nflag(18),'(2x,''Cartesian coordinates are supplied by'',/, 2x,''the user in the array CART.'',//)')
        write(nflag(18),'(2x,''Provide cartesian coordinates in the'',/, 2x,''following order using the array CART'',' // &
            '//, 2x,'' CART(1,1)...CART(1,3)   => ATOM 1'',/, 2x,'' CART(2,1)...CART(2,3)   => AT' // &
            'OM 2'',/, 2x,'' CART(3,1)...CART(3,3)   => ATOM 3'',/, 2x,'' CART(N,1)...CART(N,3)  ' // &
            ' => ATOM N'',/, 2x,''CART(25,1)...CART(25,3)  => ATOM 25'',/)')
        write(nflag(18),'(2x,''If the user wishes to relabel the atoms,'',/, 2x,''set the variable IREORDER equal to ' // &
            '1'',/, 2x,''in the PARAMETER statement.  The user'',/, 2x,''must also specify the ne' // &
            'w labeling'',/, 2x,''scheme.  This is done using the array'',/, 2x,''NULBL in the fo' // &
            'llowing manner:'',//, 2x,''NULBL(i) = j'',/, 2x,''where:  i => old index'',/, 2x,'' ' // &
            '       j => new index'',//)')
        write(nflag(18),'(2x,''Cartesian coordinates can be provided to'',/, 2x,''the potential routine in a variety ' // &
            'of units.'',/, 2x,''The input units will be converted to Bohr'',/, 2x,''based on the' // &
            ' following values of the NFLAG'',/, 2x,''variable:'',//, 2x,''NFLAG(1)  =  1  =>  CA' // &
            'RTESIANS IN BOHR (no'',/, 2x,''                    conversion required)'',/, 2x,''NF' // &
            'LAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS'',//)')
        write(nflag(18),'(2x,''The value of the energy and derivatives'',/, 2x,''(if computed) can be reported in a v' // &
            'ariety'',/, 2x,''units.  A units conversion will take place'',/, 2x,''as specified b' // &
            'y the following values of the'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(2) = 1 =>  ' // &
            'ENERGIES REPORTED IN HARTEEE'',/, 2x,''NFLAG(2) = 2 =>  ENERGIES REPORTED IN mHARTRE' // &
            'E'',/, 2x,''NFLAG(2) = 3 =>  ENERGIES REPORTED IN eV'',/, 2x,''NFLAG(2) = 4 =>  ENER' // &
            'GIES REPORTED IN kcal/mol'',/, 2x,''NFLAG(2) = 5 =>  ENERGIES REPORTED IN cm**-1'',/' // &
            '/)')
        write(nflag(18),'(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
            'he'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(1)=1 & NFLAG(2)=1 => DERIVATIVES REPOR' // &
            'TED IN'',/, 2x,''                           HARTEEE/BOHR'',/, 2x,''NFLAG(1)=1 & NFLA' // &
            'G(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/BOHR' // &
            ''',/, 2x,''NFLAG(1)=1 & NFLAG(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''            ' // &
            '               eV/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=4 => DERIVATIVES REPORTED IN' // &
            ''',/, 2x,''                           kcal/mol/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=' // &
            '5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm**-1/BOHR'',//)')
        write(nflag(18),'(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
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


!> @brief Evaluate ANCVRT for the FH2_LEPS_NO4_1984 potential.
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
    common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
    common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte, &
        & cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    dimension ianum(7,32)
    dimension isave(natom),jsave(natom)
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
    do i=1,7
        do j=1,32
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
    do i=1,natoms
        isave(i) = 0
        jsave(i) = 0
        name1(i) = '  '
        name2(i) = '  '
    end do
    iblank = ' '
    do ind=1,natoms
        do i=1,7
            do j=1,32
                if(indexes(ind).eq.ianum(i,j)) then
                    isave(ind) = i
                    jsave(ind) = j
                end if
            end do
        end do
    end do

    do ind=1,natoms
        ind2 = nulbl(ind)
        if(ind2.eq.0) ind2=ind
    end do
    inc1 = 0
    do ind=1,irctnt-1
        inc1        = inc1+1
        name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
    end do
    inc2 = 0
    do ind=irctnt,natoms
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
    endif
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
    endif
    cnvrtde = cnvrte/cnvrtd
    isum    = 0
    do inu=1,25
        isum = isum + nulbl(inu)
    end do
    ireorder = 0
    if(isum.ne.0) ireorder = 1
end subroutine ancvrt
!> @brief Evaluate CARTOU for the FH2_LEPS_NO4_1984 potential.
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
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte, &
        & cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    call convert_cartesian_input( cart, nulbl, natoms, cnvrtd, ireorder == 1, cartnu )
end subroutine cartou

!> @brief Evaluate CARTTOR for the FH2_LEPS_NO4_1984 potential.
subroutine carttor
    use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
    implicit none
    integer :: coordinate_status

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
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    call cartesian_to_potlib_coordinates( cartnu, natoms, icartr, r, coordinate_status )
    if (coordinate_status /= 0) then
        write(nflag(18),'(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') icartr
        stop
    end if
end subroutine carttor


!> @brief Evaluate EUNITZERO for the FH2_LEPS_NO4_1984 potential.
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
    common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
    common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
    common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte, &
        & cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    pengygs = engygs * cnvrte - anuzero
    if(ksdiag.ne.0) then
        do i=ksdiag,kediag
            pengyes(i) = engyes(i) * cnvrte - anuzero
        end do
    endif
    if(ksoffd.ne.0) then
        do j=ksoffd,keoffd
            pengyij(j) = engyij(j) * cnvrte
        end do
    endif
end subroutine eunitzero

!> @brief Evaluate RTOCART for the FH2_LEPS_NO4_1984 potential.
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
    common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
    common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
    common /utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte,cnvrtde, &
        & ireorder,ksdiag,kediag,ksoffd,keoffd
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    dimension ygs(n3atom),yes(n3atom,isurf),yij(n3atom,jsurf)
    if(icartr.eq.1) then
        do i = 1, natoms
            ind            = 3*i-2
            dgscartnu(i,1) = degsdr(ind)
            dgscartnu(i,2) = degsdr(ind+1)
            dgscartnu(i,3) = degsdr(ind+2)
            if(ksdiag.ne.0) then
                do j = ksdiag,kediag
                    descartnu(i,1,j) = deesdr(ind,j)
                    descartnu(i,2,j) = deesdr(ind+1,j)
                    descartnu(i,3,j) = deesdr(ind+2,j)
                end do
            endif
            if(keoffd.ne.0) then
                do k = ksoffd,keoffd
                    dijcartnu(i,1,k) = deijdr(ind,k)
                    dijcartnu(i,2,k) = deijdr(ind+1,k)
                    dijcartnu(i,3,k) = deijdr(ind+2,k)
                end do
            endif
        end do
    elseif(icartr.eq.2) then
        do i = 1, natoms
            dgscartnu(i,1) = 0.d0
            dgscartnu(i,2) = 0.d0
            dgscartnu(i,3) = 0.d0
            if(ksdiag.ne.0) then
                do j1=ksdiag,kediag
                    descartnu(i,1,j1) = 0.d0
                    descartnu(i,2,j1) = 0.d0
                    descartnu(i,3,j1) = 0.d0
                enddo
            endif
            if(ksoffd.ne.0) then
                do j2=ksoffd,keoffd
                    dijcartnu(i,1,j2) = 0.d0
                    dijcartnu(i,2,j2) = 0.d0
                    dijcartnu(i,3,j2) = 0.d0
                enddo
            endif
            do j = 1,natoms
                if(j.lt.i) then
                    m1 = natoms*(j-1) - (j*(j-1))/2 + i-j
                elseif(j.gt.i) then
                    m1 = natoms*(i-1) - (i*(i-1))/2 + j-i
                else
                    go to 20
                endif
                y              = degsdr(m1)
                termx          = (cartnu(i,1)-cartnu(j,1))/r(m1)
                termy          = (cartnu(i,2)-cartnu(j,2))/r(m1)
                termz          = (cartnu(i,3)-cartnu(j,3))/r(m1)
                dgscartnu(i,1) = dgscartnu(i,1) + termx*y
                dgscartnu(i,2) = dgscartnu(i,2) + termy*y
                dgscartnu(i,3) = dgscartnu(i,3) + termz*y
                if(ksdiag.gt.0) then
                    y = deesdr(m1,j1)
                    do j1=ksdiag,kediag
                        descartnu(i,1,j1) = descartnu(i,1,j1) + termx*y
                        descartnu(i,2,j1) = descartnu(i,2,j1) + termy*y
                        descartnu(i,3,j1) = descartnu(i,3,j1) + termz*y
                    enddo
                elseif(ksoffd.gt.0) then
                    do j2=ksoffd,keoffd
                        y                 = deijdr(m1,j2)
                        dijcartnu(i,1,j2) = dijcartnu(i,1,j2) + termx*y
                        dijcartnu(i,2,j2) = dijcartnu(i,2,j2) + termy*y
                        dijcartnu(i,3,j2) = dijcartnu(i,3,j2) + termz*y
                    enddo
                endif
                20 continue
            enddo
        enddo
    elseif(icartr.eq.3) then
        do i = 1, natoms
            ygs(i) = degsdr(i)/r(i)
            if(ksdiag.ne.0) then
                do j=ksdiag,kediag
                    yes(i,j) = deesdr(i,j)/r(i)
                enddo
            endif
            if(ksoffd.ne.0) then
                do k=ksoffd,keoffd
                    yij(i,k) = deijdr(i,k)/r(i)
                enddo
            endif
        enddo
        do k = 1,3
            term12         = cartnu(1,k)-cartnu(2,k)
            term23         = cartnu(2,k)-cartnu(3,k)
            term13         = cartnu(1,k)-cartnu(3,k)
            dgscartnu(1,k) = term12*ygs(1) + term13*ygs(3)
            dgscartnu(2,k) = -term12*ygs(1) + term23*ygs(2)
            dgscartnu(3,k) = -term13*ygs(3) - term23*ygs(2)
            if(ksdiag.ne.0) then
                do j1=ksdiag,kediag
                    descartnu(1,k,j1) = term12*yes(1,j1) + term13*yes(3,j1)
                    descartnu(2,k,j1) = -term12*yes(1,j1) + term23*yes(2,j1)
                    descartnu(3,k,j1) = -term13*yes(3,j1) - term23*yes(2,j1)
                enddo
            endif
            if(ksoffd.ne.0) then
                do j2=ksoffd,keoffd
                    dijcartnu(1,k,j2) = term12*yij(1,j2) + term13*yij(3,j2)
                    dijcartnu(2,k,j2) = -term12*yij(1,j2) + term23*yij(2,j2)
                    dijcartnu(3,k,j2) = -term13*yij(3,j2) - term23*yij(2,j2)
                enddo
            endif
        enddo
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
        write(nflag(18),'(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') icartr
        stop
    endif
end subroutine rtocart

!> @brief Evaluate DEDCOU for the FH2_LEPS_NO4_1984 potential.
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
    common /utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte,cnvrtde, &
        & ireorder,ksdiag,kediag,ksoffd,keoffd
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    if (ireorder.eq.1) then
        do i = 1, natoms
            dgscart(i,1) = dgscartnu(nulbl(i),1) * cnvrtde
            dgscart(i,2) = dgscartnu(nulbl(i),2) * cnvrtde
            dgscart(i,3) = dgscartnu(nulbl(i),3) * cnvrtde
            if(ksdiag.ne.0) then
                do j=ksdiag,kediag
                    descart(i,1,j) = descartnu(nulbl(i),1,j) * cnvrtde
                    descart(i,2,j) = descartnu(nulbl(i),2,j) * cnvrtde
                    descart(i,3,j) = descartnu(nulbl(i),3,j) * cnvrtde
                end do
            endif
            if(ksoffd.ne.0) then
                do k=ksoffd,keoffd
                    dijcart(i,1,k) = dijcartnu(nulbl(i),1,k) * cnvrtde
                    dijcart(i,2,k) = dijcartnu(nulbl(i),2,k) * cnvrtde
                    dijcart(i,3,k) = dijcartnu(nulbl(i),3,k) * cnvrtde
                end do
            endif
        end do
    else
        do i = 1, natoms
            dgscart(i,1) = dgscartnu(i,1) * cnvrtde
            dgscart(i,2) = dgscartnu(i,2) * cnvrtde
            dgscart(i,3) = dgscartnu(i,3) * cnvrtde
            if(ksdiag.ne.0) then
                do j=ksdiag,kediag
                    descart(i,1,j) = descartnu(i,1,j) * cnvrtde
                    descart(i,2,j) = descartnu(i,2,j) * cnvrtde
                    descart(i,3,j) = descartnu(i,3,j) * cnvrtde
                end do
            endif
            if(ksoffd.ne.0) then
                do k=ksoffd,keoffd
                    dijcart(i,1,k) = dijcartnu(i,1,k) * cnvrtde
                    dijcart(i,2,k) = dijcartnu(i,2,k) * cnvrtde
                    dijcart(i,3,k) = dijcartnu(i,3,k) * cnvrtde
                end do
            endif
        end do
    endif
end subroutine dedcou
