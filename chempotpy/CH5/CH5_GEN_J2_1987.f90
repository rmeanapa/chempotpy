!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the CH5_GEN_J2_1987 potential.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
    implicit none
    ! number of electronic state
    integer, parameter            :: nstates = 1
    integer, parameter            :: natoms  = 6
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
! System:          CH5
! Functional form:
! Common name:     J2
! Reference  :     T. Joseph, R. Steckler, and D. G. Truhlar
!                  J. Chem. Phys. Vol. 87, p. 7036, 1987
!
! Number of bodies: 6
! Number of electronic states: 1
! Number of derivatives: 0
! Interface: 6-1S
! Notes: This routine requires the utility 'utility.f' located in the
!        'utilities' section of POTLIB-online.
!
!***********************************************************************
!> @brief Evaluate POT for the CH5_GEN_J2_1987 potential.
subroutine pot
    !
    !***********************************************************************
    !     SUBROUTINE TO CALCULATE THE CH4-H SURFACE
    !     ALSO COMPUTES DE/DR FOR ALL R AND STORES RESULTS IN ARRAY DER
    !     DEFINITION OF CARTESIAN COORDINATES X(I), Y(I), AND Z(I)
    !      I = 1, 2, 3, 4   REFERS TO H1, H2, H3, H4
    !      I = 5            REFERS TO C
    !      I = 6            REFERS TO H'
    !
    !     DEFINITION OF THE R-ARRAY  AR(I)
    !     H1 IS THE ATOM THAT WILL FORM H-H WITH H5
    !      I = 1, 2, 3, 4   REFERS TO C - H(I)
    !      I = 5            REFERS TO C - H'
    !      I = 6, 7, 8, 9   REFERS TO H' - H(I)
    !      I = 10, 11, 12   REFERS TO H1 - H(2,3,4)
    !      I = 13, 14       REFERS TO H2 - H(3,4)
    !      I = 15           REFERS TO H3 - H4
    !***********************************************************************
    !
    implicit none
    double precision :: a11, a12, a13, a14, a21, a22, a23, a24, aa, acs, ag, alph, anuzero, aparm, ar, asb, asf, ass, b1, &
        & beta, c1, cart, cartnu, cc, coord, cp, d1, d3, daa, dadr1, dalp, dasb, dasf, deccr1, deccr5, deccr6, dedx, &
        & deesdr, degsdr, deijdr, der, der2, der3, der4, descart, df, dfr, dgscart, dijcart, dk, dr2, drdx, dtemp, dtot, &
        & dvoop, dxx, e, ee, energy, engyes, engygs, engyij, etem, evcon, exa11, exa12, exa13, exa14, exa21, exa22, &
        & exa23, exa24, exf1, exf2, exf3, exf4, ezero, f, o, pengyes, pengygs, pengyij, r, r2, re, refv, rhmax, rtem, &
        & stp, sum, tau, v, voop, w, x, xa, xa2, xb, xb2, xc, xc2, xd, xd2, xxa, xxb, xxc, xxd, y, y1, y2, z
    integer :: i, ib1, ib2, ib3, icartr, if1, if2, if3, indexes, ipt, irctnt, klm, kmax, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    dimension dfr(4), dvoop(18), y1(6), y2(6), f(4), df(4), dk(6)
    !
    dimension coord(18), dxx(18)
    !
    common /icorcm/ ar(15), der(15)
    !
    common /potxcm/ x(6), y(6), z(6), energy, drdx(15,18), dedx(18)
    common /potrcm/ ag(6), r2(15), acs(6), ass(6)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    common /dervcp/ deccr1, deccr5, deccr6
    common /temcom/ rtem(9), dtemp(9), etem
    common /alpdcm/ dalp
    data evcon / 2.721161d+1 /
    !
    call cartou
    call carttor
    !
    !     PUT COORDINATES IN PROPER ARRAYS
    !
    do 1 i=1, 18
        coord(i) = r(i)
    1 continue
    !
    ipt = 1
    do 5 i = 1, 6
        x(i) = coord(ipt)
        y(i) = coord(ipt+1)
        z(i) = coord(ipt+2)
        ipt  = ipt+3
    5 continue
    !
    ! COMPUTE THE INTERPARTICLE DISTANCES R
    !
    call dista
    !
    !***********************************************************************
    ! CHECK FOR THE FIFTH ATOM--I.E.--CH4---(H)
    ! BY DEFINITION THE FIFTH H ATOM IS THE ONE WHOSE C-H DISTANCE IS MAXIM
    !***********************************************************************
    !
    rhmax = ar(1)
    kmax  = 1
    if( ar(2) .lt. rhmax ) go to 10
    rhmax = ar(2)
    kmax  = 2
    10 if (ar(3).lt.rhmax) go to 20
    rhmax = ar(3)
    kmax  = 3
    20 if (ar(4).lt.rhmax) go to 30
    rhmax = ar(4)
    kmax  = 4
    30 if (ar(5).lt.rhmax) go to 40
    rhmax = ar(5)
    kmax  = 5
    40 call switch (kmax, 1, v, coord, dxx)
    !
    ! COMPUTE THE HCH ANGLES AND STORE RESULTS IN ARRAY AG
    !
    call angle
    !
    do 50 i = 1, 15
        der(i) = 0.0d0
    50 continue
    !
    !  COMPUTE THE MORSE RANGE PARAMETERS ALPH:
    !
    call calpha
    !
    ! COMPUTE THE FOUR TRIATOMIC TERMS
    !
    call tri (ar(1), ar(5), ar(6), der(1), der(5), der(6), 1, 3, 2, ee, dtot)
    der(1) = der(1)+dtot
    der2   = dtot
    der3   = dtot
    der4   = dtot
    call tri (ar(2), ar(5), ar(7), der(2), dr2, der(7), 1, 3, 2, e, dtot)
    ee     = ee+e
    der(2) = der(2)+der2+dtot
    der(1) = der(1)+dtot
    der3   = der3+dtot
    der4   = der4+dtot
    der(5) = der(5)+dr2
    call tri (ar(3), ar(5), ar(8), der(3), dr2, der(8), 1, 3, 2, e, dtot)
    ee     = ee+e
    der(3) = der(3)+der3+dtot
    der(1) = der(1)+dtot
    der(2) = der(2)+dtot
    der4   = der4+dtot
    der(5) = der(5)+dr2
    call tri (ar(4), ar(5), ar(9), der(4), dr2, der(9), 1, 3, 2, e, dtot)
    ee     = ee+e
    der(4) = der(4)+der4+dtot
    der(1) = der(1)+dtot
    der(2) = der(2)+dtot
    der(3) = der(3)+dtot
    der(5) = der(5)+dr2
    !
    ! COMPUTE FORCE CONSTANTS
    ! COMPUTE ATTENUATION TERMS
    !
    xa    = ar(1)-re(1)
    xb    = ar(2)-re(1)
    xc    = ar(3)-re(1)
    xd    = ar(4)-re(1)
    xa2   = xa*xa
    xb2   = xb*xb
    xc2   = xc*xc
    xd2   = xd*xd
    xxa   = ar(6)-re(2)
    xxb   = ar(7)-re(2)
    xxc   = ar(8)-re(2)
    xxd   = ar(9)-re(2)
    exa11 = exp(-aparm(2)*r2(6))
    exa12 = exp(-aparm(2)*r2(7))
    exa13 = exp(-aparm(2)*r2(8))
    exa14 = exp(-aparm(2)*r2(9))
    exa21 = exp(-aparm(5)*xxa*xxa)
    exa22 = exp(-aparm(5)*xxb*xxb)
    exa23 = exp(-aparm(5)*xxc*xxc)
    exa24 = exp(-aparm(5)*xxd*xxd)
    a11   = 1.0d0-exa11
    a12   = 1.0d0-exa12
    a13   = 1.0d0-exa13
    a14   = 1.0d0-exa14
    a21   = aparm(3)+aparm(4)*exa21
    a22   = aparm(3)+aparm(4)*exa22
    a23   = aparm(3)+aparm(4)*exa23
    a24   = aparm(3)+aparm(4)*exa24
    exf1  = exp(-a21*xa2)
    exf2  = exp(-a22*xb2)
    exf3  = exp(-a23*xc2)
    exf4  = exp(-a24*xd2)
    f(1)  = a11*exf1
    f(2)  = a12*exf2
    f(3)  = a13*exf3
    f(4)  = a14*exf4
    !
    ! DK'S ARE THE FORCE CONSTANTS
    !     FIRST COMPUTE THE K0 SINCE IT IS NOW A FUNCTION OF AR(1)
    !
    call force (dadr1)
    dk(1) = aparm(1)*f(1)*f(2)
    dk(2) = aparm(1)*f(1)*f(3)
    dk(3) = aparm(1)*f(1)*f(4)
    dk(4) = aparm(1)*f(2)*f(3)
    dk(5) = aparm(1)*f(2)*f(4)
    dk(6) = aparm(1)*f(3)*f(4)
    !
    ! COMPUTE THE EQUILIBRIUM ANGLES--THETA AND D(THETA)/DR
    !
    call theta (asf, dasf, asb, dasb)
    daa = 0.0d0
    !
    ! COMPUTE CONTRIBUTION TO THE TOTAL ENERGY FROM THE ANGLE TERMS
    !
    sum = 0.0d0
    if1 = 1
    if2 = 2
    if3 = 3
    ib1 = 4
    ib2 = 5
    ib3 = 6
    do 60 i = 1, 3
        o     = ag(i)-asf
        y1(i) = o
        y2(i) = o*o
        sum   = sum+dk(i)*y2(i)
    60 continue
    do 70 i = 4, 6
        y1(i) = ag(i)-asb
        y2(i) = y1(i)*y1(i)
        sum   = sum+dk(i)*y2(i)
    70 continue
    ee = ee+0.5d0*sum
    !
    ! COMPUTE DERIVATIVES OF THE ANGLE ATTENUATION TERMS WITH RESPECT TO R6,
    !                       AND R9
    !
    dfr(1) = 2.0d0*aparm(2)*ar(6)*exa11*exf1+a11*exf1*xa2*2.0d0* aparm(4)*aparm(5)*xxa*exa21
    dfr(2) = 2.0d0*aparm(2)*ar(7)*exa12*exf2+a12*exf2*xb2*2.0d0* aparm(4)*aparm(5)*xxb*exa22
    dfr(3) = 2.0d0*aparm(2)*ar(8)*exa13*exf3+a13*exf3*xc2*2.0d0* aparm(4)*aparm(5)*xxc*exa23
    dfr(4) = 2.0d0*aparm(2)*ar(9)*exa14*exf4+a14*exf4*xd2*2.0d0* aparm(4)*aparm(5)*xxd*exa24
    !
    ! COMPUTE DERIVATIVES OF ATTENUATION TERMS WITH RESPECT TO R1,R2,R3, AND
    !
    df(1) = -2.0d0*a21*xa*f(1)
    df(2) = -2.0d0*a22*xb*f(2)
    df(3) = -2.0d0*a23*xc*f(3)
    df(4) = -2.0d0*a24*xd*f(4)
    !
    ! ADD IN CONTRIBUTION TO DER FROM DK/DR ANGLE TERMS AND D(THETA)/DR TERM
    ! FIRST THREE ANGLE TERMS
    !
    do 80 i = 1, 3
        stp = 0.5d0*dk(i)*y2(i)
        w   = dk(i)*y1(i)
        if( abs(f(1)) .gt. 1.0d-10 ) der(1)=der(1)+stp*(df(1)/f(1)+dadr1/ aparm(1))+w*(-1.0d0/(ar(i+1)*ass(i))+acs(i)/ &
            & (ar(1)*ass(i)))
        if( abs(f(i+1)) .gt. 1.0d-10 ) der(i+1)=der(i+1)+stp*df(i+1)/f(i+ 1)+w*(-1.0d0/(ar(1)*ass(i))+acs(i)/(ass(i)* &
            & ar(i+1)))
        if( abs(f(1)) .gt. 1.0d-10 ) der(6)     = der(6)+stp*dfr(1)/f(1)
        if( abs(f(i+1)) .gt. 1.0d-10 ) der(i+6) = der(i+6)+stp*dfr(i+1)/f(i +1)
        !
        der(i+9) = der(i+9)+w*ar(i+9)/(ar(1)*ar(i+1)*ass(i))
    80 continue
    !
    ! FOURTH AND FIFTH ANGLE TERMS
    !
    do 90 i = 4, 5
        stp    = 0.5d0*dk(i)*y2(i)
        w      = dk(i)*y1(i)
        der(1) = der(1)+stp*dadr1/aparm(1)
        der(2) = der(2)+stp*df(2)/f(2)+w*(-1.0d0/(ass(i)*ar(i-1))+ acs(i)/(ass(i)*ar(2)))
        der(7) = der(7)+stp*dfr(2)/f(2)
        if( abs(f(i-1)) .gt. 1.0d-10 ) der(i-1)=der(i-1)+stp*df(i-1)/f(i- 1)+w*(-1.0d0/(ass(i)*ar(2))+acs(i)/(ass(i)* &
            & ar(i-1)))
        !
        if( abs(f(i-1)) .gt. 1.0d-10 ) der(i+4) = der(i+4)+stp*dfr(i-1)/f(i -1)
        !
        der(i+9) = der(i+9)+w*ar(i+9)/(ar(2)*ar(i-1)*ass(i))
    90 continue
    !
    ! CONTRIBUTION FROM SIXTH ANGLE TERM
    !
    stp     = 0.5d0*dk(6)*y2(6)
    w       = dk(6)*y1(6)
    der(1)  = der(1)+stp*dadr1/aparm(1)
    der(3)  = der(3)+stp*df(3)/f(3)+w*(-1.0d0/(ass(6)*ar(4))+acs(6)/ (ass(6)*ar(3)))
    der(8)  = der(8)+stp*dfr(3)/f(3)
    der(4)  = der(4)+stp*df(4)/f(4)+w*(-1.0d0/(ass(6)*ar(3))+acs(6)/ (ass(6)*ar(4)))
    der(9)  = der(9)+stp*dfr(4)/f(4)
    der(15) = der(15)+w*ar(15)/(ass(6)*ar(3)*ar(4))
    !
    ! CONTRIBUTION FROM EQ. ANGLE VARIATION TO D(ANGLE)/DR
    !
    do 100 klm = 1, 4
        der(klm) = der(klm)-(dk(if1)*y1(if1)+dk(if2)*y1(if2)+dk(if3)*y1 (if3))*dasf-(dk(ib1)*y1(ib1)+dk(ib2)*y1(ib2)+ &
            & dk(ib3)*y1(ib3) )*dasb
    100 continue
    der(5) = der(5)-(dk(if1)*y1(if1)+dk(if2)*y1(if2)+dk(if3)*y1(if3))* daa-(dk(ib1)*y1(ib1)+dk(ib2)*y1(ib2)+dk(ib3)* &
        & y1(ib3))*3.0d0* sin(asf)*cos(asf)*daa/sin(asb)
    !
    call switch (kmax, 2, v, coord, dxx)
    !
    energy = ee/evcon+refv
    do 110 i = 1, 15
        der(i) = der(i)/evcon
    110 continue
    !
    !      CALCULATE DRDX ARRAY
    !
    call drdxc
    !
    !      CALCULATE DEDX = DEDR(DER) * DRDX
    !
    call mult (der, 1, 15, drdx, 15, 18, dedx)
    ipt = 1
    do 120 i = 1, 6
        dxx(ipt)   = dedx(i)
        dxx(ipt+1) = dedx(i+6)
        dxx(ipt+2) = dedx(i+12)
        ipt        = ipt+3
    120 continue
    !
    !     NOW ADD IN THE OUT OF PLANE BEND TERM
    !
    call oopb (voop, dvoop, v, coord, dxx)
    energy = energy+voop
    do 130 i = 1, 18
        dxx(i) = dxx(i)+dvoop(i)
    130 continue
    !
    engygs = energy
    call eunitzero
    if( nder .ne. 0 )then
        do i=1, 18
            degsdr(i) = dxx(i)
        enddo
        call rtocart
        call dedcou
    endif
    !
end subroutine pot

!
!***********************************************************************
!
!> @brief Evaluate ANGLE for the CH5_GEN_J2_1987 potential.
subroutine angle
    !
    ! ROUTINE TO COMPUTE CH4 ANGLES FROM THE INTERPARTICLE DISTANCES
    !***********************************************************************
    !***********************************************************************
    !
    implicit none
    double precision :: aa, acs, ag, alph, anuzero, aparm, ar, ass, b1, beta, c1, cart, cartnu, cc, cp, d1, d3, dedx, &
        & deesdr, degsdr, deijdr, der, descart, dgscart, dijcart, drdx, energy, engyes, engygs, engyij, ezero, pengyes, &
        & pengygs, pengyij, r, r2, re, refv, tau, x, xx, y, z
    integer :: i, icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /icorcm/ ar(15), der(15)
    !
    common /potxcm/ xx(6), y(6), z(6), energy, drdx(15,18), dedx(18)
    common /potrcm/ ag(6), r2(15), acs(6), ass(6)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    do 10 i = 1, 3
        x      = (r2(1)+r2(i+1)-r2(i+9))/(2.0d0*ar(1)*ar(i+1))
        acs(i) = x
        ass(i) = sqrt(1.0d0-x*x)
        ag(i)  = acos(x)
    10 continue
    do 20 i = 1, 2
        x        = (r2(2)+r2(i+2)-r2(i+12))/(2.0d0*ar(2)*ar(i+2))
        acs(i+3) = x
        ass(i+3) = sqrt(1.0d0-x*x)
        ag(i+3)  = acos(x)
    20 continue
    x      = (r2(3)+r2(4)-r2(15))/(2.0d0*ar(3)*ar(4))
    acs(6) = x
    ass(6) = sqrt(1.0d0-x*x)
    ag(6)  = acos(x)
end subroutine angle

!
!
!> @brief Evaluate CALPHA for the CH5_GEN_J2_1987 potential.
subroutine calpha
    !
    !     COMPUTE THE C-H AND C-H' MORSE RANGE PARAMETERS FOR THE SINGLET
    !     AND TRIPLET ENERGY CURVES.  THESE WILL BE A FUNCTION OF THE
    !     AVERAGE C-H BOND LENGTH IN METHANE.
    !
    implicit none
    double precision :: aa, alph, anuzero, aparm, ar, b1, beta, c1, ca, cart, cartnu, cb, cc, cp, d1, d3, dalp, deesdr, &
        & degsdr, deijdr, der, descart, dgscart, dijcart, engyes, engygs, engyij, ezero, fa, fk0, htan, pengyes, pengygs, &
        & pengyij, r, ravg, re, refv, tau
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /icorcm/ ar(15), der(15)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    common /pdt3cm/ fk0, fa, ca, cb
    common /alpdcm/ dalp
    !
    !     FIRST COMPUTE THE AVERAGE METHANE BOND LENGTH
    !
    ravg = (ar(1)+ar(2)+ar(3)+ar(4))/4.0d0
    !
    htan    = tanh(c1*(ravg-re(1)))
    alph(1) = ca+cb*((htan+1.0d0)/2.0d0)
    alph(3) = alph(1)
    beta(1) = alph(1)
    beta(3) = alph(1)
    !
    !     NOW THE DERIVATIVES---DALPH WILL BE THE DERIVATIVE OF ALPHA WITH
    !     RESPECT TO R1,R2,R3 OR R4
    !
    dalp = (c1*cb/8.0d0)*(1.0d0-htan**2)
    !
end subroutine calpha

!
!***********************************************************************
!***********************************************************************
!
!> @brief Evaluate DISTA for the CH5_GEN_J2_1987 potential.
subroutine dista
    !
    ! ROUTINE TO COMPUTE THE INTERPARTICLE DISTANCES FROM THE CARTESIAN
    !     COORDINATES OF THE ATOMS
    !***********************************************************************
    !***********************************************************************
    !
    implicit none
    double precision :: aa, acs, ag, alph, anuzero, aparm, ar, ass, b1, beta, c1, cart, cartnu, cc, cp, d1, d3, dedx, &
        & deesdr, degsdr, deijdr, der, descart, df, dgscart, dijcart, drdx, energy, engyes, engygs, engyij, ezero, &
        & pengyes, pengygs, pengyij, r, r2, re, refv, tau, x, y, z
    integer :: i, icartr, indexes, irctnt, j, k, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    dimension df(15,3)
    !
    common /icorcm/ ar(15), der(15)
    !
    common /potxcm/ x(6), y(6), z(6), energy, drdx(15,18), dedx(18)
    common /potrcm/ ag(6), r2(15), acs(6), ass(6)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    do 10 i = 1, 4
        df(i,1) = x(i)-x(5)
        df(i,2) = y(i)-y(5)
        df(i,3) = z(i)-z(5)
        r2(i)   = df(i,1)*df(i,1)+df(i,2)*df(i,2)+df(i,3)*df(i,3)
        ar(i)   = sqrt(r2(i))
    10 continue
    df(5,1) = x(5)-x(6)
    df(5,2) = y(5)-y(6)
    df(5,3) = z(5)-z(6)
    r2(5)   = df(5,1)*df(5,1)+df(5,2)*df(5,2)+df(5,3)*df(5,3)
    ar(5)   = sqrt(r2(5))
    do 20 i = 1, 4
        j       = i+5
        df(j,1) = x(i)-x(6)
        df(j,2) = y(i)-y(6)
        df(j,3) = z(i)-z(6)
        r2(j)   = df(j,1)*df(j,1)+df(j,2)*df(j,2)+df(j,3)*df(j,3)
        ar(j)   = sqrt(r2(j))
    20 continue
    do 30 i = 1, 3
        j       = i+1
        k       = i+9
        df(k,1) = x(1)-x(j)
        df(k,2) = y(1)-y(j)
        df(k,3) = z(1)-z(j)
        r2(k)   = df(k,1)*df(k,1)+df(k,2)*df(k,2)+df(k,3)*df(k,3)
        ar(k)   = sqrt(r2(k))
    30 continue
    do 40 i = 1, 2
        j       = i+2
        k       = i+12
        df(k,1) = x(2)-x(j)
        df(k,2) = y(2)-y(j)
        df(k,3) = z(2)-z(j)
        r2(k)   = df(k,1)*df(k,1)+df(k,2)*df(k,2)+df(k,3)*df(k,3)
        ar(k)   = sqrt(r2(k))
    40 continue
    df(15,1) = x(3)-x(4)
    df(15,2) = y(3)-y(4)
    df(15,3) = z(3)-z(4)
    r2(15)   = df(15,1)*df(15,1)+df(15,2)*df(15,2)+df(15,3)*df(15,3)
    ar(15)   = sqrt(r2(15))
end subroutine dista

!> @brief Evaluate DRDXC for the CH5_GEN_J2_1987 potential.
subroutine drdxc
    !
    ! THIS SUBROUTINE CALCULATES DRDX
    !
    !     ADDED BY FBB AT THE U. OF MINN.  1/24/84
    !
    implicit none
    double precision :: aa, acs, ag, alph, anuzero, aparm, ar, ass, b1, beta, c1, cart, cartnu, cc, cp, d1, d3, dedx, &
        & deesdr, degsdr, deijdr, der, descart, dgscart, dijcart, drdx, energy, engyes, engygs, engyij, ezero, pengyes, &
        & pengygs, pengyij, r, r2, re, refv, tau, x, y, z
    integer :: i, icartr, indexes, irctnt, j, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /icorcm/ ar(15), der(15)
    !
    common /potxcm/ x(6), y(6), z(6), energy, drdx(15,18), dedx(18)
    common /potrcm/ ag(6), r2(15), acs(6), ass(6)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    !
    !
    ! ZERO OUT DRDX
    !
    do 10 i = 1, 18
        do 10 j = 1, 15
            drdx(j,i) = 0.0d0
    10 continue
    !
    ! CALCULATE DRDX FOR C - H(1,2,3,4) BONDS
    !
    do 20 i = 1, 4
        drdx(i,i)    = (x(i)-x(5))/ar(i)
        drdx(i,5)    = -drdx(i,i)
        drdx(i,i+6)  = (y(i)-y(5))/ar(i)
        drdx(i,11)   = -drdx(i,i+6)
        drdx(i,i+12) = (z(i)-z(5))/ar(i)
        drdx(i,17)   = -drdx(i,i+12)
    20 continue
    !
    ! CALCULATE DRDX FOR C - H' BOND
    !
    drdx(5,5)  = (x(5)-x(6))/ar(5)
    drdx(5,6)  = -drdx(5,5)
    drdx(5,11) = (y(5)-y(6))/ar(5)
    drdx(5,12) = -drdx(5,11)
    drdx(5,17) = (z(5)-z(6))/ar(5)
    drdx(5,18) = -drdx(5,17)
    !
    ! CALCULATE DRDX FOR H - H' BONDS
    !
    do 30 i = 1, 4
        drdx(5+i,i)    = (x(i)-x(6))/ar(5+i)
        drdx(5+i,6)    = -drdx(5+i,i)
        drdx(5+i,i+6)  = (y(i)-y(6))/ar(5+i)
        drdx(5+i,12)   = -drdx(5+i,i+6)
        drdx(5+i,i+12) = (z(i)-z(6))/ar(5+i)
        drdx(5+i,18)   = -drdx(5+i,i+12)
    30 continue
    !
    ! CALCULATE DRDX FOR H1 - H(2,3,4) BONDS
    !
    do 40 i = 1, 3
        drdx(9+i,i+1)  = -(x(1)-x(i+1))/ar(9+i)
        drdx(9+i,1)    = -drdx(9+i,i+1)
        drdx(9+i,i+7)  = -(y(1)-y(i+1))/ar(9+i)
        drdx(9+i,7)    = -drdx(9+i,i+7)
        drdx(9+i,i+13) = -(z(1)-z(i+1))/ar(9+i)
        drdx(9+i,13)   = -drdx(9+i,i+13)
    40 continue
    !
    ! CALCULATE DRDX FOR H2 - H(3,4)
    !
    do 50 i = 1, 2
        drdx(12+i,i+2)  = -(x(2)-x(i+2))/ar(12+i)
        drdx(12+i,2)    = -drdx(12+i,i+2)
        drdx(12+i,i+8)  = -(y(2)-y(i+2))/ar(12+i)
        drdx(12+i,8)    = -drdx(12+i,i+8)
        drdx(12+i,i+14) = -(z(2)-z(i+2))/ar(12+i)
        drdx(12+i,14)   = -drdx(12+i,i+14)
    50 continue
    !
    ! CALCULATE DRDX FOR H3 - H4
    !
    drdx(15,4)  = -(x(3)-x(4))/ar(15)
    drdx(15,3)  = -drdx(15,4)
    drdx(15,10) = -(y(3)-y(4))/ar(15)
    drdx(15,9)  = -drdx(15,10)
    drdx(15,16) = -(z(3)-z(4))/ar(15)
    drdx(15,15) = -drdx(15,16)
    !
end subroutine drdxc

!
!
!> @brief Evaluate FORCE for the CH5_GEN_J2_1987 potential.
subroutine force ( dadr1 )
    !
    !     THIS SUBROUTINE COMPUTES THE BENDING FORCE CONSTANT AS A FUNCTION
    !     OF AR(1). (APARM(1))   IT ALSO COMPUTES THE DERIVATIVE OF THIS
    !     TERM WITH RESPECT TO AR(1).
    !
    implicit none
    double precision :: aa, alph, anuzero, aparm, ar, b1, beta, c1, ca, cart, cartnu, cb, cc, cp, d1, d3, dadr1, deesdr, &
        & degsdr, deijdr, der, descart, dgscart, dijcart, engyes, engygs, engyij, ezero, fa, fex, fk0, pengyes, pengygs, &
        & pengyij, r, re, refv, tau
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    !
    common /icorcm/ ar(15), der(15)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    common /pdt3cm/ fk0, fa, ca, cb
    !
    fex      = exp(-b1*(ar(1)-re(1))**2)
    aparm(1) = fk0+fa*fex
    dadr1    = -2.0d0*fa*b1*(ar(1)-re(1))*fex
end subroutine force

!> @brief Evaluate MULT for the CH5_GEN_J2_1987 potential.
subroutine mult ( a, na, ma, b, nb, mb, c )
    implicit none
    double precision :: a, b, c, sum
    integer          :: i, j, k, ma, mb, na, nb
    !
    ! THIS SUBROUTINE PERFORMS THE MATRIX MULTIPLICATION :
    !
    !         C = A X B
    !
    !
    dimension a(na,ma), b(nb,mb), c(na,mb)
    !
    do 30 i = 1, na
        !
        do 20 j = 1, mb
            !
            sum = 0.0d0
            do 10 k = 1, ma
                !
                sum = sum+a(i,k)*b(k,j)
            10 continue
            !
            c(i,j) = sum
        20 continue
        !
    30 continue
end subroutine mult

!> @brief Evaluate OOPB for the CH5_GEN_J2_1987 potential.
subroutine oopb ( voop, dvoop, v, x, dx )
    !
    !     THIS SUBROUTINE CALCULTES THE OUT OF PLANE BEND POTENTIAL THAT
    !     RESULTS FROM THE METHYL.  THE FORM USED IS THE SAME AS THAT
    !     GIVEN BY DUCHOVIC ET. AL. (JPC 89, 1339, (1984)).
    !     THE DERIVATIVES WITH RESPECT TO THE CARTESIAN COORDINATES (X)
    !     ARE ALSO COMPUTED AND STORED IN DVOOP.
    !
    implicit none
    double precision :: a0, a3, anuzero, aphi, b0, b3, bphi, c0, cart, cartnu, cphi, crossx, crossy, crossz, cst, d, &
        & dcross, ddelta, deesdr, degsdr, deijdr, delt, delta, descart, dfds, dfdx, dgamma, dgcros, dgscart, dhds, dhdx, &
        & dij, dijcart, dphi, drdx, drv, dsdr, dtop, dvoop, dx, dxmagr, engyes, engygs, engyij, erest, ezero, f, fch3, &
        & gamma, gcross, h, hch3, pch4, pengyes, pengygs, pengyij, phinot, r, rest, rex, rnot, rstar, rv, s3, sech, sphi, &
        & term1, term2, term3, term4, texp, top, v, voop, x, xdif, xmagr, ydif, zdif
    integer :: i, icartr, indexes, irctnt, j, k, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    dimension dvoop(18), xdif(2), ydif(2), zdif(2), rv(3,3), drv(3,3,18)
    dimension xmagr(3), dxmagr(3,18), delta(3), ddelta(3,18)
    dimension gamma(3), dgamma(3,18), top(3), dtop(3,18), dcross(3,18)
    dimension dgcros(18), dphi(18), dfdx(18), dhdx(18), drdx(18)
    !
    !
    dimension x(18), dx(18)
    !
    common /testcm/ gamma, dgamma
    common /pdt2cm/ aphi, bphi, cphi, pch4, fch3, hch3, rnot, a3, b3, delt, &
        & dij, a0, b0, c0, d, rex
    !
    !     START BY COMPUTING THE NECESSARY DIFFERENCES IN THE R1,R2 AND
    !     R3 VECTORS.
    !
    xdif(1) = x(7)-x(4)
    ydif(1) = x(8)-x(5)
    zdif(1) = x(9)-x(6)
    xdif(2) = x(10)-x(4)
    ydif(2) = x(11)-x(5)
    zdif(2) = x(12)-x(6)
    !
    !     NOW COMPUTE THE INDIVIDUAL COMPONENTS OF THE R1,R2 AND R3 VECTORS
    !
    j = 1
    do 10 i = 1, 3
        rv(i,1) = x(j+3)-x(13)
        rv(i,2) = x(j+4)-x(14)
        rv(i,3) = x(j+5)-x(15)
        j       = j+3
    10 continue
    !
    !     RSTAR IS THE C-H4 BOND LENGTH
    !
    rstar = sqrt((x(1)-x(13))**2+(x(2)-x(14))**2+(x(3)-x(15))**2)
    !
    !     COMPUTE THE CROSS PRODUCT OF THESE DIFFERENCES
    !
    crossx = ydif(1)*zdif(2)-ydif(2)*zdif(1)
    crossy = zdif(1)*xdif(2)-xdif(1)*zdif(2)
    crossz = xdif(1)*ydif(2)-xdif(2)*ydif(1)
    gcross = sqrt(crossx**2+crossy**2+crossz**2)
    !
    !     COMPUTE THE THREE TRIPLE PRODUCTS--WITH R1,R2 AND R3 (TOP)
    !     AND THEN COMPUTE THE ANGLE GAMMA USED IN DELTA.
    !
    do 20 i = 1, 3
        xmagr(i) = sqrt(rv(i,1)**2+rv(i,2)**2+rv(i,3)**2)
        top(i)   = crossx*rv(i,1)+crossy*rv(i,2)+crossz*rv(i,3)
        gamma(i) = +top(i)/(gcross*xmagr(i))
    20 continue
    !
    !     NOW THE SECOND TERM IN DELTA--PHINOT
    !     FIRST THE SWITCHING FUNCTION SPHI THEN PHINOT
    !
    sphi = 0.0d0
    rest = bphi*(rstar-cphi)**3
    if( rest .lt. 70.0d0 )then
        erest = exp(rest)
        sphi  = 1.0d0-tanh(aphi*(rstar-rnot)*erest)
    endif
    phinot = pch4+(pch4-pi/2.0d0)*(sphi-1.0d0)
    !
    !     NOW PUT THESE TWO TERMS TOGETHER TO GET EACH DELTA
    !
    do 30 i = 1, 3
        delta(i) = acos(gamma(i))-phinot
    30 continue
    !
    !     NEXT COMPUTE THE PARTIALS OF DELTA W/R TO EACH X(I)
    !      IN DOING THIS COMPUTATION FIRST FOR EACH TERM CALCULATED
    !      ABOVE THE PARTIAL WILL BE COMPUTED AND A STORED USING THE
    !      SAME VARIABLE WITH A 'D' IN FRONT.
    !
    !
    !     START WITH THE PARTIAL OF CROSS WITH RESPECT TO EACH X(I)
    !
    do 60 i = 1, 3
        do 40 j = 1, 3
            dcross(i,j) = 0.0d0
        40 continue
        do 50 j = 13, 18
            dcross(i,j) = 0.0d0
        50 continue
        if( i .eq. 1 )then
            dcross(i,4)  = 0.0d0
            dcross(i,5)  = zdif(1)-zdif(2)
            dcross(i,6)  = ydif(2)-ydif(1)
            dcross(i,7)  = 0.0d0
            dcross(i,8)  = zdif(2)
            dcross(i,9)  = -ydif(2)
            dcross(i,10) = 0.0d0
            dcross(i,11) = -zdif(1)
            dcross(i,12) = ydif(1)
        else if( i .eq. 2 )then
            dcross(i,4)  = zdif(2)-zdif(1)
            dcross(i,5)  = 0.0d0
            dcross(i,6)  = xdif(1)-xdif(2)
            dcross(i,7)  = -zdif(2)
            dcross(i,8)  = 0.0d0
            dcross(i,9)  = xdif(2)
            dcross(i,10) = zdif(1)
            dcross(i,11) = 0.0d0
            dcross(i,12) = -xdif(1)
        else if( i .eq. 3 )then
            dcross(i,4)  = ydif(1)-ydif(2)
            dcross(i,5)  = xdif(2)-xdif(1)
            dcross(i,6)  = 0.0d0
            dcross(i,7)  = ydif(2)
            dcross(i,8)  = -xdif(2)
            dcross(i,9)  = 0.0d0
            dcross(i,10) = -ydif(1)
            dcross(i,11) = xdif(1)
            dcross(i,12) = 0.0d0
        endif
    60 continue
    !
    !     NOW COMPUTE THE PARTIALS FOR THE THREE R VECTORS
    !
    do 70 i = 1, 3
        do 70 j = 1, 3
            !            DO 70 K = 1, 3
            ! BUGFIX:  DRV was not fully initialized to zero causing trouble on some
            !          platforms.
            do 70 k = 1, 18
                drv(i,j,k) = 0.0d0
    70 continue
    do 80 k = 4, 6
        drv(1,k-3,k) = 1.0d0
    80 continue
    do 90 k = 7, 9
        drv(2,k-6,k) = 1.0d0
    90 continue
    do 100 k = 10, 12
        drv(3,k-9,k) = 1.0d0
    100 continue
    do 110 i = 1, 3
        do 110 j = 1, 3
            drv(i,j,j+12) = -1.0d0
    110 continue
    !
    !     NOW COMPUTE THE PARTIALS FOR THE TRIPLE PRODUCT  (TOP)
    !
    do 120 i = 1, 18
        do 120 j = 1, 3
            dtop(j,i) = dcross(1,i)*rv(j,1)+crossx*drv(j,1,i)+dcross(2,i )*rv(j,2)+crossy*drv(j,2,i)+ dcross(3,i)* &
                & rv(j,3)+crossz* drv(j,3,i)
    120 continue
    !
    !     PARTIALS FOR THE MAGNITUDE OF THE CROSS PRODUCT  (GCROSS)
    !
    do 130 i = 1, 18
        dgcros(i) = (1.0d0/gcross)*(dcross(1,i)*crossx+dcross(2,i)* crossy+dcross(3,i)*crossz)
    130 continue
    !
    !     PARTIALS FOR THE MAGNITUDE OF R1, R2 AND R3
    !
    do 150 i = 1, 3
        do 150 j = 1, 18
            dxmagr(i,j) = 0.0d0
            do 140 k = 1, 3
                dxmagr(i,j) = dxmagr(i,j)+(1.0d0/xmagr(i))*(rv(i,k)*drv(i,k,j))
            140 continue
    150 continue
    !
    !     NOW WE'RE READY FOR THE PARTIALS OF GAMMA
    !
    do 160 i = 1, 3
        do 160 j = 1, 18
            dgamma(i,j) = +gcross*xmagr(i)*dtop(i,j)-top(i)*(xmagr(i)* dgcros(j)+gcross*dxmagr(i,j))
            dgamma(i,j) = dgamma(i,j)/(gcross*xmagr(i))**2
            dgamma(i,j) = dgamma(i,j)*(-1.0d0)/(sqrt(1.0d0-gamma(i)**2) )
    160 continue
    !
    !     NOW THE PARTIAL FOR THE SECOND TERM IN DELTA -- PHINOT
    !       FIRST COMPUTE THE PARTIAL OF S(SWITCHING FUNCTION) WITH
    !       RESPECT TO RSTAR, THE THE PARTIAL OF RSTAR W/R X(I) AND
    !       FINALLY THE PARTIAL OF PHINOT W/R TO S.
    !
    dsdr = 0.0d0
    rest = bphi*(rstar-cphi)**3
    if( rest .lt. 70.0d0 )then
        texp = exp(rest)
        dsdr = -((sech(aphi*(rstar-rnot)*texp))**2)
        dsdr = dsdr*aphi*texp*(1.0d0+(rstar-rnot)*bphi*3.0d0*(rstar- cphi)**2)
    endif
    do 170 i = 1, 18
        drdx(i) = 0.0d0
        dphi(i) = 0.0d0
    170 continue
    !
    cst = 1.0d0/(2.0d0*rstar)
    do 180 i = 1, 3
        drdx(i)    = cst*2.0d0*(x(i)-x(i+12))
        drdx(i+12) = -drdx(i)
        dphi(i)    = drdx(i)*dsdr*(pch4-pi/2.0d0)
        dphi(i+12) = -dphi(i)
    180 continue
    !
    !     NOW PUT BOTH TERMS TOGETHER TO GET THE PARTIALS OF DELTA
    !
    do 190 i = 1, 3
        do 190 j = 1, 18
            ddelta(i,j) = dgamma(i,j)-dphi(j)
    190 continue
    !
    !     WE HAVE NOW COMPUTED ALL THE DELTA TERMS.  NEXT COMPUTE THE
    !       FORCE CONSTANTS (F AND H) AND THEIR DERIVATIVES.  THE
    !       OUT-OF-PLANE BEND ENERGY WILL ALSO BE COMPUTED HERE (VOOP)
    !
    !     FIRST THE SWITCHING FUNCTION
    !
    s3   = 1.0d0-tanh(a3*(rstar-rnot)*(rstar-b3)**2)
    f    = (1.0d0-s3)*fch3
    h    = (1.0d0-s3)*hch3
    voop = 0.0d0
    do 200 i = 1, 3
        voop = voop+f*(delta(i)**2)+h*(delta(i)**4)
    200 continue
    !
    !     PARTIAL OF THE SWITCHING FUNCTION W/R TO RSTAR
    !
    dsdr = -(sech(a3*(rstar-rnot)*(rstar-b3)**2)**2)*a3*((rstar-b3)**2 +(rstar-rnot)*2*(rstar-b3))
    dfds = -fch3
    dhds = -hch3
    !
    !     PARTIAL OF EACH FORCE CONSTANT
    !
    do 210 i = 1, 18
        dfdx(i) = 0.0d0
        dhdx(i) = 0.0d0
    210 continue
    do 220 i = 1, 3
        dfdx(i)    = dfds*dsdr*drdx(i)
        dhdx(i)    = dhds*dsdr*drdx(i)
        dfdx(i+12) = -dfdx(i)
        dhdx(i+12) = -dhdx(i)
    220 continue
    !
    !    NOW WE ARE READY FOR THE FULL PARTIAL DERIVATIVES OF VOOP
    !
    do 250 i = 1, 18
        term1 = 0.0d0
        do 230 j = 1, 3
            term1 = term1+dfdx(i)*delta(j)**2
        230 continue
        term2 = 2.0d0*delta(1)*ddelta(1,i)+2.0d0*delta(2)*ddelta(2,i) +2.0d0*delta(3)*ddelta(3,i)
        term2 = term2*f
        term3 = 0.0d0
        do 240 j = 1, 3
            term3 = term3+dhdx(i)*delta(j)**4
        240 continue
        term4    = 4.0d0*delta(1)**3*ddelta(1,i)+4.0d0*delta(2)**3* ddelta(2,i)+4.0d0*delta(3)**3*ddelta(3,i)
        term4    = term4*h
        dvoop(i) = term1+term2+term3+term4
    250 continue
end subroutine oopb

!
!> @brief Evaluate SECH for the CH5_GEN_J2_1987 potential.
function sech ( x )
    implicit none
    double precision :: sech, x
    if( x .gt. 70.0d0 )then
        sech = 0.0d0
    else
        sech = 1.0d0/cosh(x)
    endif
end function sech

!> @brief Evaluate SWITCH for the CH5_GEN_J2_1987 potential.
subroutine switch ( kmax, lmg, v, coord, dxx )
    !
    ! **********************************************************************
    ! THIS ROUTINE SWITCHES MATRICES R,R2,DER TO ALLOW FOR DIFFERENT H ATOMS
    !      BEING OFF THE CH4 MOLECULE.
    ! LMG IS A CONTROL PARAMETER
    ! KMAX IS THE LABEL OF THE H ATOM WITH MAXIMUM C-H DISTANCE.
    ! IF LMG=1, R AND R2 MATRICES ARE SWITCHED APPROPRIATELY.
    ! IF LMG=2, R AND R2 MATRICES ARE SWITCHED BACK AND MATRIX DER IS SWITCH
    ! **********************************************************************
    !
    implicit none
    double precision :: aa, acs, ag, alph, anuzero, aparm, ar, ass, b1, beta, c1, cart, cartnu, cc, coord, cp, d1, d3, &
        & dedx, deesdr, degsdr, deijdr, der, descart, dgscart, dijcart, drdx, dxx, energy, engyes, engygs, engyij, ezero, &
        & pengyes, pengygs, pengyij, r, r2, re, refv, store, tau, v, x, y, z
    integer :: icartr, indexes, irctnt, k, kk, kmax, l, ll, lmg, m, mder, mm, msurf, n, nasurf, natoms, nder, nflag, nn, &
        & nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    dimension coord(18), dxx(18)
    !
    common /icorcm/ ar(15), der(15)
    !
    common /potxcm/ x(6), y(6), z(6), energy, drdx(15,18), dedx(18)
    common /potrcm/ ag(6), r2(15), acs(6), ass(6)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    go to (10, 20, 30, 40, 50), kmax
    !
    ! H1 IS MAXIMUM
    !
    10 k = 1
    kk = 5
    l  = 7
    ll = 10
    m  = 8
    mm = 11
    n  = 9
    nn = 12
    go to 60
    !
    ! H2 IS MAXIMUM
    !
    20 k = 2
    kk = 5
    l  = 6
    ll = 10
    m  = 8
    mm = 13
    n  = 9
    nn = 14
    go to 60
    !
    ! H3 IS MAXIMUM
    !
    30 k = 3
    kk = 5
    l  = 6
    ll = 11
    m  = 7
    mm = 13
    n  = 9
    nn = 15
    go to 60
    !
    ! H4 IS MAXIMUM
    !
    40 k = 4
    kk = 5
    l  = 6
    ll = 12
    m  = 7
    mm = 14
    n  = 8
    nn = 15
    go to 60
    !
    ! H5 IS MAXIMUM
    ! NO SWITCHES ARE REQUIRED
    !
    50 go to 80
    !
    ! SWITCH R AND R2 MATRICES
    !
    60 store = ar(k)
    ar(k)  = ar(kk)
    ar(kk) = store
    store  = ar(l)
    ar(l)  = ar(ll)
    ar(ll) = store
    store  = ar(m)
    ar(m)  = ar(mm)
    ar(mm) = store
    store  = ar(n)
    ar(n)  = ar(nn)
    ar(nn) = store
    r2(k)  = ar(k)*ar(k)
    r2(kk) = ar(kk)*ar(kk)
    r2(l)  = ar(l)*ar(l)
    r2(ll) = ar(ll)*ar(ll)
    r2(m)  = ar(m)*ar(m)
    r2(mm) = ar(mm)*ar(mm)
    r2(n)  = ar(n)*ar(n)
    r2(nn) = ar(nn)*ar(nn)
    if( lmg-1 ) 80, 80, 70
    !
    ! SWITCH DER MATRIX *******
    !
    70 store = der(k)
    der(k)  = der(kk)
    der(kk) = store
    store   = der(l)
    der(l)  = der(ll)
    der(ll) = store
    store   = der(m)
    der(m)  = der(mm)
    der(mm) = store
    store   = der(n)
    der(n)  = der(nn)
    der(nn) = store
    80 return
end subroutine switch

!> @brief Evaluate THETA for the CH5_GEN_J2_1987 potential.
subroutine theta ( asf, dasf, asb, dasb )
    !
    !     THIS SUBROUTINE COMPUTES THE EQUILIBRIUM THETA AS A FUNCTION OF
    !     R (THE AVERAGE C-H BOND LENGTH IN METHANE).
    !     THE DERIVATIVE OF THETA WITH RESPECT TO AR(I) IS ALSO COMPUTED.
    !     NOTE HERE THAT BECAUSE OF SYMMETRY EACH OF THESE PARTIALS IS
    !     THE SAME.       ADDED 11/22/85
    !
    implicit none
    double precision :: a, a0, a3, aa, alph, anuzero, aparm, aphi, ar, asb, asf, b, b0, b1, b3, beta, bphi, c, c0, c1, &
        & cart, cartnu, cc, cp, cphi, ctau, ctau2, d0, d1, d3, dasb, dasf, deesdr, degsdr, deijdr, delta, delta2, der, &
        & descart, dgscart, dij, dijcart, engyes, engygs, engyij, ezero, f, fch3, fk1, fk2, hch3, pch4, pengyes, pengygs, &
        & pengyij, r, ravg, rdif, rdif2, rdif3, rdif4, rdif5, re, refv, rex, rnot, tau, xyz
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /icorcm/ ar(15), der(15)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    common /pdt2cm/ aphi, bphi, cphi, pch4, fch3, hch3, rnot, a3, b3, delta, &
        & dij, a0, b0, c0, d0, rex
    !
    ravg = (ar(1)+ar(2)+ar(3)+ar(4))/4.0d0
    !
    !     COMPUTE THETA
    !
    if( ravg .gt. dij )then
        asf  = pi/2.0d0
        dasf = 0.0d0
    else if( ravg .le. re(1) )then
        asf  = tau
        dasf = 0.0d0
    else if( ravg .gt. re(1).and.ravg .le. re(1)+delta )then
        delta2 = delta*delta
        fk1    = -4.0d0*cp/(delta2)
        a      = 3.0d0*fk1/(delta2)
        b      = -(5.0d0*a*delta2+fk1)/(2.0d0*delta)
        c      = -(20.0d0*a*delta2+12.0d0*b*delta)/6.0d0
        f      = tau
        rdif   = ravg-re(1)
        rdif2  = rdif*rdif
        rdif3  = rdif2*rdif
        rdif4  = rdif3*rdif
        rdif5  = rdif4*rdif
        asf    = a*rdif5+b*rdif4+c*rdif3+f
        dasf   = 5.0d0*a*rdif4+4.0d0*b*rdif3+3.0d0*c*rdif2
    else if( ravg .le. dij-delta )then
        asf  = tau-4.0d0*cp*(ravg-re(1))
        dasf = -4.0d0*cp
    else if( ravg .le. dij.and.ravg .gt. dij-delta )then
        fk1   = (tau-4.0d0*cp*(dij-delta-re(1))-f)/(delta**3)
        fk2   = (-4.0d0*cp)/(delta**2)
        a     = -3.0d0*(2.0d0*fk1+fk2)/(delta**2)
        b     = (fk2+5.0d0*a*delta**2)/(2.0d0*delta)
        c     = (12.0d0*b*delta-20.0d0*a*delta**2)/6.0d0
        f     = pi/2.0d0
        rdif  = ravg-dij
        rdif3 = rdif*rdif*rdif
        rdif4 = rdif3*rdif
        rdif5 = rdif4*rdif
        asf   = a*rdif5+b*rdif4+c*rdif3+f
        dasf  = 5.0d0*a*rdif4+4.0d0*b*rdif3+3.0d0*c*rdif*rdif
    endif
    !
    ! COMPUTE BACKSIDE EQUILIBRIUM ANGLE
    !
    xyz   = cos(asf)
    ctau  = cos(tau)
    ctau2 = ctau*ctau
    a     = (ctau+0.5d0)/ctau2
    asb   = acos(a*xyz*xyz-0.5d0)
    dasb  = 2.0d0*a*xyz*sin(asf)*dasf/sin(asb)
    !
    ! NOW COMPUTE THE DER W/R TO RI INSTEAD OF RAVG
    !
    dasf = dasf/4.0d0
    dasb = dasb/4.0d0
end subroutine theta

!
!***********************************************************************
!
!> @brief Evaluate TRI for the CH5_GEN_J2_1987 potential.
subroutine tri ( rr1, rr2, rr3, dr1, dr2, dr3, i, j, k, e, dtot )
    !
    ! ROUTINE TO CALCULATE THREE-BODY ENERGY FOR ABC SYSTEM
    ! RR1,RR2,RR3 ARE THE THREE INTERPARTICLE DISTANCES
    ! DR1,DR2,DR3 ARE THE THREE DERIVATIVES WITH RESPECT TO RR1,RR2, AND RR3
    ! I,J, AND K GIVE THE ARRAY NUMBER FOR THE APPROPRIATE POTENTIAL PARAMET
    ! E IS THE FINAL THREE-BODY ENERGY
    !***********************************************************************
    !***********************************************************************
    !
    implicit none
    double precision :: a, a1, a2, aa, aa1, aa2, acs, ag, alab, alac, albc, alph, anuzero, aparm, ar, ass, b, b1, b2, bb, &
        & bb1, bb2, beta, bp1, c, c1, c2, cart, cartnu, cc, cc1, cc2, cp, cp1, d1, d3, dab1, dab3, dac1, dac3, dalp, &
        & de1ab, de1ac, de1bc, de3ab, de3ac, de3bc, deccr1, deccr5, deccr6, dedx, deesdr, degsdr, deijdr, der, descart, &
        & dgscart, dijcart, dqr, dr1, dr2, dr3, drdx, dtot, dur, dux, duy, duz, dxr, dyr, dzr, e, e1ab, e1ac, e1bc, e3ab, &
        & e3ac, e3bc, energy, engyes, engygs, engyij, ezero, pengyes, pengygs, pengyij, qab, qac, qbc, r, r2, re, refv, &
        & rr1, rr2, rr3, tau, u, x, xx, y, yy, z, zz
    integer :: i, icartr, indexes, irctnt, j, k, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /icorcm/ ar(15), der(15)
    !
    common /potxcm/ x(6), y(6), z(6), energy, drdx(15,18), dedx(18)
    common /potrcm/ ag(6), r2(15), acs(6), ass(6)
    common /dervcp/ deccr1, deccr5, deccr6
    common /alpdcm/ dalp
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, bp1, cp1
    a     = rr1-re(i)
    b     = rr2-re(j)
    c     = rr3-re(k)
    a1    = exp(-alph(i)*a)
    b1    = exp(-alph(j)*b)
    c1    = exp(-alph(k)*c)
    a2    = a1*a1
    b2    = b1*b1
    c2    = c1*c1
    e1ab  = d1(i)*(a2-2.0d0*a1)
    e1ac  = d1(j)*(b2-2.0d0*b1)
    e1bc  = d1(k)*(c2-2.0d0*c1)
    de1ab = 2.0d0*d1(i)*(a1-a2)*alph(i)
    dab1  = 2.0d0*d1(i)*(a1-a2)*(a*dalp)
    de1ac = 2.0d0*d1(j)*(b1-b2)*alph(j)
    dac1  = 2.0d0*d1(j)*(b1-b2)*(b*dalp)
    de1bc = 2.0d0*alph(k)*d1(k)*(c1-c2)
    aa1   = exp(-beta(i)*a)
    aa2   = aa1*aa1
    e3ab  = d3(i)*(aa2+2.0d0*aa1)
    de3ab = -2.0d0*d3(i)*(aa1+aa2)*beta(i)
    dab3  = -2.0d0*d3(i)*(aa1+aa2)*a*dalp
    bb1   = exp(-beta(j)*b)
    bb2   = bb1*bb1
    e3ac  = d3(j)*(bb2+2.0d0*bb1)
    de3ac = -2.0d0*d3(j)*(bb1+bb2)*beta(j)
    dac3  = -2.0d0*d3(j)*(bb1+bb2)*b*dalp
    cc1   = exp(-beta(k)*c)
    cc2   = cc1*cc1
    e3bc  = d3(k)*(cc2+2.0d0*cc1)
    de3bc = -2.0d0*beta(k)*d3(k)*(cc1+cc2)
    qab   = (e1ab+e3ab)/2.0d0
    qac   = (e1ac+e3ac)/2.0d0
    qbc   = (e1bc+e3bc)/2.0d0
    alab  = (e1ab-e3ab)/2.0d0
    alac  = (e1ac-e3ac)/2.0d0
    albc  = (e1bc-e3bc)/2.0d0
    xx    = alab-albc
    yy    = albc-alac
    zz    = alac-alab
    u     = (xx*xx+yy*yy+zz*zz)/2.0d0
    u     = sqrt(u)
    e     = qab+qbc+qac-u
    bb    = 4.0d0*u
    dxr   = (dab1-dab3)/2.0d0
    dyr   = (dac3-dac1)/2.0d0
    dzr   = -(dxr+dyr)
    dux   = xx/(2.0d0*u)
    duy   = yy/(2.0d0*u)
    duz   = zz/(2.0d0*u)
    dur   = dux*dxr+duy*dyr+duz*dzr
    dqr   = (dab1+dac1+dab3+dac3)/2.0d0
    dtot  = dqr-dur
    dr1   = (de1ab+de3ab)/2.0d0-(2.0d0*alab-albc-alac)*(de1ab-de3ab)/bb
    dr2   = (de1ac+de3ac)/2.0d0-(2.0d0*alac-alab-albc)*(de1ac-de3ac)/bb
    dr3   = (de1bc+de3bc)/2.0d0-(2.0d0*albc-alac-alab)*(de1bc-de3bc)/bb
end subroutine tri

!> @brief Evaluate TRIPLE for the CH5_GEN_J2_1987 potential.
subroutine triple
    !
    !     THIS SUBROUTINE COMPUTES THE C-H' SATO PARAMETER WHICH IS A
    !     FUNCTION OF AR(H'-H) AND THE C-H-H' BEND ANGLE.
    !
    implicit none
    double precision :: a, a3, aa, alph, alpha, anuzero, aparm, aphi, ar, b, b1, b3, beta, bphi, c0, c1, cart, cartnu, &
        & cc, check, cp, cp1, cphi, d, d1, d3, dalr1, dalr5, dalr6, dc1r6, dcyx, deccr1, deccr5, deccr6, deesdr, deg20, &
        & degsdr, deijdr, delt, der, descart, dgscart, dij, dijcart, dyxr, engyes, engygs, engyij, ezero, fch3, hch3, &
        & pch4, pengyes, pengygs, pengyij, r, r6, re, refv, rex, rnot, tau, tem, tem2, tem22, tembot, u, ui, yx
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom, isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom, jsurf)
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
    common /icorcm/ ar(15), der(15)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, cp1
    common /pdt2cm/ aphi, bphi, cphi, pch4, fch3, hch3, rnot, a3, b3, delt, &
        & dij, a, b, c0, d, rex
    common /dervcp/ deccr1, deccr5, deccr6
    !
    !     COMPUTE C1--IT IS A FUNCTION OF R6--THE H-H' BOND DISTANCE
    !
    if( ar(6) .gt. 15.0d0 )then
        r6 = ar(6)/3.0d0
        yx = (exp(2.0d0*pi*(r6-rex)/d))**3
    else
        yx = exp(2.0d0*pi*(ar(6)-rex)/d)
    endif
    c1    = -(a*yx/(yx+1.0d0)+b/yx)-c0
    dyxr  = (2.0d0*pi/d)*yx
    dcyx  = b/yx-(a*yx/(1+yx))/(1+yx)
    dc1r6 = dcyx*dyxr/yx
    !
    !     COMPUTE THE C-H4-H5 ANGLE
    !
    tem    = ar(5)**2-ar(1)**2-ar(6)**2
    tembot = -2.0d0*ar(1)*ar(6)
    check  = tem/tembot
    if( check .le. -1.0d0 )then
        alpha  = 0.0d0
        cc(3)  = c0
        deccr1 = 0.0d0
        deccr5 = 0.0d0
        deccr6 = 0.0d0
        return
    else
        alpha = acos(tem/(-2.0d0*ar(1)*ar(6)))
        alpha = pi-alpha
    endif
    !
    !     NOW COMPUTE THE ANGLE DERIVATIVES (DERIV. ALPHA W.R. R1,R5,R6)
    !
    u     = tem/(-2.0d0*ar(1)*ar(6))
    ui    = -1.0d0/(sqrt(1.0d0-u**2))
    dalr5 = -ui*ar(5)/(ar(1)*ar(6))
    tem2  = 1.0d0/(-2.0d0*ar(1)*ar(6))
    tem22 = tem2*tem2
    dalr1 = -2.0d0*ui*(ar(1)*tem2-ar(6)*tem*tem22)
    dalr6 = -2.0d0*ui*(ar(6)*tem2-ar(1)*tem*tem22)
    !
    !     COMBINE THE RADIAL AND ANGLE COMPONENT TO COMPUTE THE TRIPLET
    !     PARAMETER CC(3) AND THE DERIVATIVES.  HERE TEM IS JUST THE PARTIAL
    !     OF CC W.R.T. ALPHA.
    !
    deg20 = pi/9.0d0
    if( alpha .le. deg20 )then
        cc(3)  = c0+c1*(sin(9.0d0*alpha/2.0d0))**4
        tem    = 18.0d0*(sin(9.0d0*alpha/2.0d0))**3*cos(9.0d0*alpha/ 2.0d0)
        deccr1 = c1*tem*dalr1
        deccr5 = c1*tem*dalr5
        deccr6 = dc1r6*(sin(9.0d0*alpha/2.0d0))**4+c1*tem*dalr6
    else
        cc(3)  = c0+c1
        deccr1 = 0.0d0
        deccr5 = 0.0d0
        deccr6 = dc1r6
    endif
end subroutine triple

!
!> @brief Evaluate PREPOT for the CH5_GEN_J2_1987 potential.
subroutine prepot
    !
    !   N3TMMN = 3 * NATOMS
    !   NATOMS = the number of atoms represented by this potential function
    !
    !   The variable N3TMMN is the minimum value of N3TM allowed to be
    !   passed by the calling routine for the number of cartesian
    !   coordinates needed to represent the full system represented by this
    !   potential energy surface routine.
    !   N3TM must be greater than or equal to N3TMMN.
    !
    implicit none
    double precision :: aa, acs, ag, alph, anuzero, aparm, ar, ass, b1, beta, c1, cart, cartnu, cc, coord, cp, d1, d3, &
        & dalp, deccr1, deccr5, deccr6, dedx, der, df, dfr, dk, drdx, dtemp, dvoop, dxx, energy, etem, evcon, ezero, f, &
        & r2, re, refv, rtem, tau, x, y, y1, y2, z
    integer :: i, icartr, indexes, irctnt, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    integer, parameter          :: n3tmmn = 18
    !
    common/pt3cm/ ezero(isurf+1)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    dimension dfr(4), dvoop(18), y1(6), y2(6), f(4), df(4), dk(6)
    !
    dimension coord(18), dxx(18)
    !
    common /icorcm/ ar(15), der(15)
    !
    common /potxcm/ x(6), y(6), z(6), energy, drdx(15,18), dedx(18)
    common /potrcm/ ag(6), r2(15), acs(6), ass(6)
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    common /dervcp/ deccr1, deccr5, deccr6
    common /temcom/ rtem(9), dtemp(9), etem
    common /alpdcm/ dalp
    data evcon / 2.721161d+1 /
    !
    !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
    !
    if( natoms .gt. 25 )then
        write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    endif
    !
    do i=1, 5
        ref(i) = ' '
    enddo
    !
    ref(1) = 'T. Joseph, R. Steckler, and D. G. Truhlar'
    ref(2) = 'J. Chem. Phys. Vol. 87, p. 7036, 1987'
    !
    indexes(1) = 1
    indexes(2) = 1
    indexes(3) = 1
    indexes(4) = 1
    indexes(5) = 6
    indexes(6) = 1
    !
    !
    !
    irctnt = 6
    !
    !      CALL POTINFO
    !
    call ancvrt
    !
end subroutine prepot

!
block data ptpacm
    !
    implicit none
    double precision :: a, a3, aa, alph, anuzero, aparm, aphi, b, b1, b3, beta, bphi, c0, c1, ca, cart, cartnu, cb, cc, &
        & cp, cphi, d, d1, d3, delt, dij, ezero, fa, fch3, fk0, hch3, pch4, re, refv, rex, rnot, tau
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt3cm/ ezero(isurf+1)
    !
    common/infocm/ cartnu(natom, 3), indexes(natom), &
        & irctnt, natoms, icartr, mder, msurf, ref
    !
    !
    common/usricm/ cart(natom, 3), anuzero, &
        & nulbl(natom), nflag(20), &
        & nasurf(isurf+1, isurf+1), nder
    !
    common /pdatcm/ d1(3), d3(3), alph(3), re(3), beta(3), cc(3), aa(3), &
        & aparm(5), refv, tau, cp, b1, c1
    common /pdt2cm/ aphi, bphi, cphi, pch4, fch3, hch3, rnot, a3, b3, delt, &
        & dij, a, b, c0, d, rex
    common /pdt3cm/ fk0, fa, ca, cb
    !
    data nasurf /1, 35*0/
    data nder /1/
    data nflag /1, 1, 15*0, 6, 0, 0/
    !
    data anuzero /0.0d0/
    data icartr, msurf, mder/1, 0, 1/
    data nulbl /25*0/
    data natoms /6/
    !
    data d1 /  4.8668d0, 4.74660d0, 1.1452d0 /
    data d3 /  1.6840d0, 1.72000d0, 0.8700d0 /
    data b1 / 3.00d0 /
    data c1 / 3.50d0 /
    data alph / 0.939d0, 1.029955d0, 0.939d0 /
    data re / 2.0673d0, 1.402d0, 2.0673d0 /
    data beta / 0.939d0, 1.029955d0, 0.939d0 /
    data aa / 0.0d0, 0.0d0, -2.3603d0 /
    data refv / 0.71098329d0 /
    data tau / 1.9105588d0 /
    data cp / 0.080768d0 /
    data aparm / 4.0d0, 0.8999975d0, 0.448035d0, 0.606528d0, 3.239921d0 /
    !
    data aphi / 2.7982366d-1 /
    data bphi / 5.937217d-2 /
    data cphi / 3.6301534d0 /
    data pch4 / 1.9106d0 /
    data fch3 / 2.19d-2 /
    data hch3 / 4.38d-02 /
    data rnot / 2.052d0 /
    data a3 / 2.102956d-2 /
    data b3 / 5.798533d-1 /
    data delt / 0.1d0 /
    data dij / 3.1190d0 /
    data a / -1.56944d+3 /
    data b / -1.12008d+3 /
    data c0 / 1.5786151d+3 /
    data d / 2.67146d0 /
    data ca / 0.903d0 /
    data cb / 0.072d0 /
    data rex / 1.18584d0 /
    data fk0 / 2.55d0 /
    data fa / 0.8d0 /
end

!> @brief Evaluate POTINFO for the CH5_GEN_J2_1987 potential.
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

!> @brief Evaluate ANCVRT for the CH5_GEN_J2_1987 potential.
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

!> @brief Evaluate CARTOU for the CH5_GEN_J2_1987 potential.
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

!> @brief Evaluate CARTTOR for the CH5_GEN_J2_1987 potential.
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

!> @brief Evaluate EUNITZERO for the CH5_GEN_J2_1987 potential.
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

!> @brief Evaluate RTOCART for the CH5_GEN_J2_1987 potential.
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

!> @brief Evaluate DEDCOU for the CH5_GEN_J2_1987 potential.
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
