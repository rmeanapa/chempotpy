!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the FH2_GEN_6SEC_1993 potential.
subroutine pes( x, igrad, p, g, d )

    implicit none

    ! number of electronic state
    integer, parameter :: nstates = 1
    integer, parameter :: natoms  = 3
    integer, parameter :: natom   = 25
    integer, parameter :: isurf   = 5
    integer, parameter :: jsurf   = int(isurf*(isurf+1)/2)

    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)

    double precision :: anuzero, cart, descart, dgscart, dijcart, pengyes, pengygs, pengyij, v
    integer          :: iatom, idir, istate, nasurf, nder, nflag, nulbl

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
!> @brief Evaluate PREPOT for the FH2_GEN_6SEC_1993 potential.
subroutine prepot
    !
    !   System:             FH2
    !   Functional form:    Double many-body expansion
    !   Common name:        6SEC
    !   Reference:          S. L. Mielke, G. C. Lynch,
    !                       D. G. Truhlar, and D. W. Schwenke
    !                       Chem. Phys. Lett. 216, 441-446 (1993)
    !
    ! Number of bodies: 3
    ! Interface: potlib2001
    ! Number of electronic states: 1
    ! Number of derivatives: 1
    !
    !
    !   PREPOT must be called once before any call to POT.
    !   The potential parameters are included through the subprogram
    !   BLOCK DATA PTPACM.
    !   The coordinates, potential energy for the ground electronic state, a
    !   the derivatives of the potential energy for the ground electronic st
    !   with respect to the coordinates are passed through the common block
    !   The potential energy in the three asymptotic valleys are
    !   stored in the common block ASYCM:
    !                  COMMON /PT5COM/ EASYAB, EASYBC, EASYAC
    !   The potential energy in the AB valley, EASYAB, is equal to the poten
    !   energy of the H "infinitely" far from the FH diatomic, with the
    !   FH diatomic at its equilibrium configuration.  Similarly, the terms
    !   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
    !   respectively.
    !   All the information passed through the common blocks PT1CM and ASYCM
    !   is in Hartree atomic units.
    !
    !        This potential is written such that:
    !                       R(1) = R(F-H)
    !                       R(2) = R(H-H)
    !                       R(3) = R(H-F)
    !
    !   The flags that indicate what calculations should be carried out in
    !   the potential routine are passed through the common block PT2CM:
    !   where:
    !        NASURF - which electronic states availalble
    !                 (1,1) = 1 as only gs state available
    !        NDER  = 0 => no derivatives should be calculated
    !        NDER  = 1 => calculate first derivatives
    !        NFLAG - these integer values can be used to flag options
    !                 within the potential;
    !
    !
    !   Potential parameters' default settings
    !                  Variable            Default value
    !                  NDER                1
    !                  NFLAG(18)           6
    !
    !   ENERGY  = E(EHF) + E(CORR) + De(H2)
    !   E(EHF)  = E(London) + E(3-body)
    !   E(CORR) = ECORR(2-body) + ECORR(3-body)
    !
    !   Note:  The zero of energy for this surface is F infinitely far
    !          from H2 at R(H2) = Re(H2).
    !          To get the zero of energy R1 must be at least 50 a.u.
    !
    implicit none

    double precision :: a1, a2, a3, ad1, ad2, ad3, alf, alfp, alft, alph0, alph1, anuzero, bd1, bd2, bd3, bet0, bet1, &
        & betp, bhf, bhh, c3c, cart, cartnu, cfhh, chf, chh, chhf, ck1, ck2, ck3, ckp1, ckp2, ckp3, d, dcorr2, dcorr3, &
        & ddp1, ddp2, ddp3, de3c, dec, disse, dldr, dp1, dp2, dp3, drtanh, dsw, e3c, easyab, easyac, easybc, ecor2b, &
        & ecor3b, elond, eta1, eta2, eta3, etap1, etap2, etap3, ezero, fcsi, fcsi1, fcsi2, gamma, p, p3c, q3c, r0, rb, &
        & rm, rmt, ro, rtanh, rtmp, sw, vdum, xhf, xhh
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
    common /potdcm/ vdum(10)
    common /decm/  disse(3)
    !
    dimension       rtmp(3)
    common /mint/   rmt(3)
    common /dist/   rm(3),r0(3)
    common /ref/    ro(3)
    common /tpar/   xhh(3),xhf(3),p(4)
    common /atat/   chh(10),chf(10)
    common /coeff/  cfhh,chhf
    common /trial/  eta1(10),eta2(10),eta3(10), &
        & etap1(10),etap2(10),etap3(10)
    common /kkp/    ck1(10),ck2(10),ck3(10),ckp1(10),ckp2(10),ckp3(10)
    common /diat/   d(3),a1(3),a2(3),a3(3),gamma(3)
    common /swpar/  dec(3),rb(3),sw(3),dsw(3,3)
    common /csi/    fcsi,fcsi1,fcsi2
    common /local/  c3c,alf,alfp,alft,p3c,q3c,betp
    common /dispc/  alph0,alph1,bet0,bet1
    !
    common /dampc/  ad1,ad2,ad3,bd1,bd2,bd3
    common /dpc/    dp1(3),dp2(3),dp3(3),ddp1(3),ddp2(3),ddp3(3)
    common /degsdrcm/ dcorr2(3),dcorr3(3),dldr(3),de3c(3)
    common /engycm/ ecor2b, ecor3b, elond, e3c
    common /par/    bhh(5),bhf(5)
    common /vtanh/  rtanh(3,10),drtanh(3,10)
    !
    double precision, parameter :: deh2 = 0.174472692d0
    double precision, parameter :: dehf = 0.224989333d0
    !
    if(natoms.gt.25) then
        write(nflag(18),'(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    end if
    !
    !     CONSTANTS NEEDED FOR THE DAMPING FUNCTION
    !     A=ALPH0/DBLE(N)**ALPH1        N=6,8,10
    !     B=BET0*EXP(-BET1*DBLE(N))     N=6,8,10
    !
    ad1 = alph0/6.0d0**alph1
    ad2 = alph0/8.0d0**alph1
    ad3 = alph0/10.0d0**alph1
    bd1 = bet0*exp(-bet1*6.0d0)
    bd2 = bet0*exp(-bet1*8.0d0)
    bd3 = bet0*exp(-bet1*10.0d0)
    !
    !   Initialize the energy in the asymptotic valleys
    !
    easyab = disse(1)
    easybc = disse(2)
    easyac = disse(3)
    !
    !
    !
    !
    ezero(1) = disse(2)
    !
    do i=1,5
        ref(i) = ' '
    end do
    !   Reference:
    !                       D. G. Truhlar, and D. W. Schwenke
    !
    !
    ref(1) = 'S. L. Mielke, G. C. Lynch,'
    ref(2) = 'D. G. Truhlar, D. W. Schwenke,'
    ref(3) = 'Chem. Phys. Lett. 216, 441(1993)'
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
!=======================================================================
!   POT: Subroutine for computing the energy and the derivatives.
!=======================================================================
!
!> @brief Evaluate POT for the FH2_GEN_6SEC_1993 potential.
subroutine pot
    !
    !        This potential is written such that:
    !                       R(1) = R(F-H)
    !                       R(2) = R(H-H)
    !                       R(3) = R(H-F)
    !
    !   The potential energy in the AB valley, EASYAB, is equal to the poten
    !   energy of the H "infinitely" far from the FH diatomic, with the
    !   FH diatomic at its equilibrium configuration.  Similarly, the terms
    !   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
    !   respectively.
    !
    !   Note:  The zero of energy for this surface is F infinitely far
    !          from H2 at R(H2) = Re(H2).
    !          To get the zero of energy R1 must be at least 50 a.u.
    !
    !      ENTRY POT
    implicit none

    double precision :: a1, a2, a3, ad1, ad2, ad3, alf, alfp, alft, alph0, alph1, anuzero, bd1, bd2, bd3, bet0, bet1, &
        & betp, bhf, bhh, c3c, cart, cartnu, cfhh, chf, chh, chhf, ck1, ck2, ck3, ckp1, ckp2, ckp3, d, dcorr2, dcorr3, &
        & ddp1, ddp2, ddp3, de3c, dec, deesdr, degsdr, deijdr, descart, dgscart, dijcart, disse, dldr, dp1, dp2, dp3, &
        & drtanh, dsw, e3c, easyab, easyac, easybc, ecor2b, ecor3b, elond, engyes, engygs, engyij, eta1, eta2, eta3, &
        & etap1, etap2, etap3, ezero, fcsi, fcsi1, fcsi2, gamma, p, p3c, pengyes, pengygs, pengyij, q3c, r, r0, rb, rm, &
        & rmt, ro, rtanh, rtmp, sw, t, vdum, xhf, xhh
    integer :: i, i1, i2, i3, icartr, indexes, irctnt, mder, msurf, nasurf, natoms, nder, nflag, nulbl
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
    common /potdcm/ vdum(10)
    common /decm/  disse(3)
    !
    dimension       rtmp(3)
    common /mint/   rmt(3)
    common /dist/   rm(3),r0(3)
    common /ref/    ro(3)
    common /tpar/   xhh(3),xhf(3),p(4)
    common /atat/   chh(10),chf(10)
    common /coeff/  cfhh,chhf
    common /trial/  eta1(10),eta2(10),eta3(10), &
        & etap1(10),etap2(10),etap3(10)
    common /kkp/    ck1(10),ck2(10),ck3(10),ckp1(10),ckp2(10),ckp3(10)
    common /diat/   d(3),a1(3),a2(3),a3(3),gamma(3)
    common /swpar/  dec(3),rb(3),sw(3),dsw(3,3)
    common /csi/    fcsi,fcsi1,fcsi2
    common /local/  c3c,alf,alfp,alft,p3c,q3c,betp
    common /dispc/  alph0,alph1,bet0,bet1
    !
    common /dampc/  ad1,ad2,ad3,bd1,bd2,bd3
    common /dpc/    dp1(3),dp2(3),dp3(3),ddp1(3),ddp2(3),ddp3(3)
    common /degsdrcm/ dcorr2(3),dcorr3(3),dldr(3),de3c(3)
    common /engycm/ ecor2b, ecor3b, elond, e3c
    common /par/    bhh(5),bhf(5)
    common /vtanh/  rtanh(3,10),drtanh(3,10)
    !
    double precision, parameter :: deh2 = 0.174472692d0
    double precision, parameter :: dehf = 0.224989333d0
    !
    call cartou
    call carttor
    !
    !   Check the values of NASURF and NDER for validity.
    !
    !
    !
    !   Test the values of Ri
    !   If the values are unrealistic then print Ri and exit
    !
    do 11 i=1,3
        rtmp(i) = r(i)
    11 continue
    do 12 i1=1,3
        i2 = mod(i1,3)+1
        i3 = mod(i2,3)+1
        t  = r(i2) + r(i3)
        if(r(i1) .gt. t)then
            write(nflag(18),'(//,10X,''R1 = '',1PE20.10,5X,''R2 = '',1PE20.10,5X, ''R3 = '',1PE20.10)') r(i1),r(i2),r(i3)
            stop
        endif
    12 continue
    !
    !   I: Compute the damping function for the N th. dispersion coefficient
    !      and, if NDER = 1, compute the derivative of the damping function
    !      w.r.t. R(i)
    !
    call damp(r,rm)
    !
    !   II: Compute the two-body dynamical correlation term and, if NDER = 1
    !       compute the derivative of the two-body dynamical correlation ter
    !       w.r.t. R(i)
    !
    call corr2
    !
    !   III: Compute the auxiliary function, hn(R), and, if NDER = 1,
    !            compute the derivative of hn(R) w.r.t. R
    !
    call xtanh(1,eta1)
    call xtanh(2,eta2)
    call xtanh(3,eta3)
    !
    !   IV: Compute the three-body dynamical correlation term and, if NDER =
    !       compute the derivative of the three-body dynamical correlation t
    !       w.r.t. R(i)
    !
    call corr3
    !
    !   V: Compute the switching function, S(R), which is used in the comput
    !      of the H2 triplet state curves.  If NDER = 1, compute the derivat
    !      of S w.r.t. R
    !
    call switch
    !
    !   VI: Compute the lower eigenvalue of the London eq. and, if NDER = 1,
    !       compute the derivatives of the London energy w.r.t. R.
    !
    call vlond
    !
    !   VII: Compute the three-body term and, if NDER = 1, compute the
    !        the derivatives of the three-body energy w.r.t. R.
    !
    call thcent
    !
    !   IX: Compute the total energy which is a sum of the extended-Hartree-
    !       and the dynamical correlation terms.
    !
    engygs = ecor2b + ecor3b + elond + e3c + ezero(1) - anuzero
    !
    !   Compute the derivatives of the energy with respect to the coordinate
    !
    if (nder .eq. 1) then
        do 10 i = 1, 3
            degsdr(i) = dcorr2(i) + dcorr3(i) + dldr(i) + de3c(i)
        10 continue
    endif
    !
    !
    call eunitzero
    if(nder.ne.0) then
        call rtocart
        call dedcou
    endif
    !
end subroutine pot
!
!*****
!
!> @brief Evaluate CORR2 for the FH2_GEN_6SEC_1993 potential.
subroutine corr2
    !=======================================================================
    !   Compute the two-body dynamical energies and derivatives.
    !=======================================================================
    !
    implicit none

    double precision :: anuzero, cart, cartnu, chf, chh, cor2r1, cor2r2, cor2r3, dcorr2, dcorr3, ddp1, ddp2, ddp3, de3c, &
        & deesdr, degsdr, deijdr, descart, dgscart, dijcart, dldr, dp1, dp2, dp3, e3c, easyab, easyac, easybc, ecor2b, &
        & ecor3b, elond, engyes, engygs, engyij, ezero, pengyes, pengygs, pengyij, r, r1p10, r1p6, r1p8, r2p10, r2p6, &
        & r2p8, r3p10, r3p6, r3p8, vdum
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
    common /potdcm/ vdum(10)
    common /atat/chh(10),chf(10)
    common /degsdrcm/dcorr2(3),dcorr3(3),dldr(3),de3c(3)
    common /dpc/dp1(3),dp2(3),dp3(3),ddp1(3),ddp2(3),ddp3(3)
    common /engycm/ecor2b, ecor3b, elond, e3c
    !
    r1p6  = r(1)**6
    r1p8  = r(1)**8
    r1p10 = r(1)**10
    r2p6  = r(2)**6
    r2p8  = r(2)**8
    r2p10 = r(2)**10
    r3p6  = r(3)**6
    r3p8  = r(3)**8
    r3p10 = r(3)**10
    !
    cor2r1 = -dp1(1)*chf(6)/r1p6 -dp1(2)*chf(8)/r1p8 -dp1(3)*chf(10)/r1p10
    !
    cor2r2 = -dp2(1)*chh(6)/r2p6 -dp2(2)*chh(8)/r2p8 -dp2(3)*chh(10)/r2p10
    !
    cor2r3 = -dp3(1)*chf(6)/r3p6 -dp3(2)*chf(8)/r3p8 -dp3(3)*chf(10)/r3p10
    !
    ecor2b = cor2r1 + cor2r2 + cor2r3
    !
    !   Compute the derivatives w.r.t. R(i) if NDER = 1
    !
    if (nder .ne. 1) go to 100
    !
    dcorr2(1) = -chf(6)*(ddp1(1)-6.0d0*dp1(1)/r(1))/r1p6 -chf(8)*(ddp1(2)-8.0d0*dp1(2)/r(1))/r1p8 -chf(10)*(ddp1(3)- &
        & 10.0d0*dp1(3)/r(1))/r1p10
    !
    dcorr2(2) = -chh(6)*(ddp2(1)-6.0d0*dp2(1)/r(2))/r2p6 -chh(8)*(ddp2(2)-8.0d0*dp2(2)/r(2))/r2p8 -chh(10)*(ddp2(3)- &
        & 10.0d0*dp2(3)/r(2))/r2p10
    !
    dcorr2(3) = -chf(6)*(ddp3(1)-6.0d0*dp3(1)/r(3))/r3p6 -chf(8)*(ddp3(2)-8.0d0*dp3(2)/r(3))/r3p8 -chf(10)*(ddp3(3)- &
        & 10.0d0*dp3(3)/r(3))/r3p10
    !
    100 continue
    !
end subroutine corr2
!
!****
!
!> @brief Evaluate DAMP for the FH2_GEN_6SEC_1993 potential.
subroutine damp( rt1, rt2 )
    !
    !=======================================================================
    !   Compute the damping function for the Nth. dispersion coefficient and
    !   derivative w.r.t. R(i)
    !=======================================================================
    !
    implicit none

    double precision :: ad1, ad2, ad3, anuzero, bd1, bd2, bd3, cart, cartnu, ddp1, ddp2, ddp3, deesdr, degsdr, deijdr, &
        & descart, dgscart, dijcart, dp1, dp2, dp3, dpol1, dpol2, dpol3, dx, easyab, easyac, easybc, engyes, engygs, &
        & engyij, ezero, pengyes, pengygs, pengyij, pol1, pol2, pol3, r, r0, rm, rt1, rt2, vdum, x
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
    common /potdcm/ vdum(10)
    common /dampc/ ad1,ad2,ad3,bd1,bd2,bd3
    common /dpc/ dp1(3),dp2(3),dp3(3),ddp1(3),ddp2(3),ddp3(3)
    common /dist/rm(3),r0(3)
    dimension rt1(3),rt2(3),x(3),dx(3),pol1(3),pol2(3),pol3(3)
    dimension dpol1(3),dpol2(3),dpol3(3)
    !
    x(1) = 2.0d0*rt1(1)/(rt2(1)+2.5d0*r0(1))
    x(2) = 2.0d0*rt1(2)/(rt2(2)+2.5d0*r0(2))
    x(3) = 2.0d0*rt1(3)/(rt2(3)+2.5d0*r0(3))
    !
    pol1(1) = ad1*x(1)+bd1*x(1)**2
    pol2(1) = ad2*x(1)+bd2*x(1)**2
    pol3(1) = ad3*x(1)+bd3*x(1)**2
    !
    pol1(2) = ad1*x(2)+bd1*x(2)**2
    pol2(2) = ad2*x(2)+bd2*x(2)**2
    pol3(2) = ad3*x(2)+bd3*x(2)**2
    !
    pol1(3) = ad1*x(3)+bd1*x(3)**2
    pol2(3) = ad2*x(3)+bd2*x(3)**2
    pol3(3) = ad3*x(3)+bd3*x(3)**2
    !
    dp1(1) = (1.0d0-exp(-pol1(1)))**6
    dp1(2) = (1.0d0-exp(-pol2(1)))**8
    dp1(3) = (1.0d0-exp(-pol3(1)))**10
    !
    dp2(1) = (1.0d0-exp(-pol1(2)))**6
    dp2(2) = (1.0d0-exp(-pol2(2)))**8
    dp2(3) = (1.0d0-exp(-pol3(2)))**10
    !
    dp3(1) = (1.0d0-exp(-pol1(3)))**6
    dp3(2) = (1.0d0-exp(-pol2(3)))**8
    dp3(3) = (1.0d0-exp(-pol3(3)))**10
    !
    if (nder .eq. 1) then
        dx(1) = 2.0d0/(rt2(1)+2.5d0*r0(1))
        dx(2) = 2.0d0/(rt2(2)+2.5d0*r0(2))
        dx(3) = 2.0d0/(rt2(3)+2.5d0*r0(3))
        !
        dpol1(1) = ad1+2.0d0*bd1*x(1)
        dpol2(1) = ad2+2.0d0*bd2*x(1)
        dpol3(1) = ad3+2.0d0*bd3*x(1)
        !
        dpol1(2) = ad1+2.0d0*bd1*x(2)
        dpol2(2) = ad2+2.0d0*bd2*x(2)
        dpol3(2) = ad3+2.0d0*bd3*x(2)
        !
        dpol1(3) = ad1+2.0d0*bd1*x(3)
        dpol2(3) = ad2+2.0d0*bd2*x(3)
        dpol3(3) = ad3+2.0d0*bd3*x(3)
        !
        ddp1(1) = 6.0d0*dp1(1)*exp(-pol1(1))*dpol1(1)*dx(1)/ (1.0d0-exp(-pol1(1)))
        ddp1(2) = 8.0d0*dp1(2)*exp(-pol2(1))*dpol2(1)*dx(1)/ (1.0d0-exp(-pol2(1)))
        ddp1(3) = 10.0d0*dp1(3)*exp(-pol3(1))*dpol3(1)*dx(1)/ (1.0d0-exp(-pol3(1)))
        !
        ddp2(1) = 6.0d0*dp2(1)*exp(-pol1(2))*dpol1(2)*dx(2)/ (1.0d0-exp(-pol1(2)))
        ddp2(2) = 8.0d0*dp2(2)*exp(-pol2(2))*dpol2(2)*dx(2)/ (1.0d0-exp(-pol2(2)))
        ddp2(3) = 10.0d0*dp2(3)*exp(-pol3(2))*dpol3(2)*dx(2)/ (1.0d0-exp(-pol3(2)))
        !
        ddp3(1) = 6.0d0*dp3(1)*exp(-pol1(3))*dpol1(3)*dx(3)/ (1.0d0-exp(-pol1(3)))
        ddp3(2) = 8.0d0*dp3(2)*exp(-pol2(3))*dpol2(3)*dx(3)/ (1.0d0-exp(-pol2(3)))
        ddp3(3) = 10.0d0*dp3(3)*exp(-pol3(3))*dpol3(3)*dx(3)/ (1.0d0-exp(-pol3(3)))
    endif
    !
end subroutine damp
!
!*****
!
!> @brief Evaluate CORR3 for the FH2_GEN_6SEC_1993 potential.
subroutine corr3
    !
    !=======================================================================
    !   Compute the three-body dynamical energies and derivatives.
    !=======================================================================
    !
    implicit none

    double precision :: anuzero, cart, cartnu, chf, chh, ck1, ck2, ck3, ckp1, ckp2, ckp3, co3, dc3dr1, dc3dr2, dc3dr3, &
        & dcorr2, dcorr3, ddp1, ddp2, ddp3, de3c, deesdr, degsdr, deijdr, descart, dg1, dg2, dg3, dgscart, dh1, dh2, dh3, &
        & dijcart, dldr, dp1, dp2, dp3, drtanh, dt1dr1, dt1dr2, dt1dr3, dt2dr1, dt2dr2, dt2dr3, dt3dr1, dt3dr2, dt3dr3, &
        & e3c, easyab, easyac, easybc, ecor2b, ecor3b, elond, engyes, engygs, engyij, eta1, eta2, eta3, etap1, etap2, &
        & etap3, ezero, fi, g1, g2, g3, h1, h2, h3, pengyes, pengygs, pengyij, r, r0, r1pi, r2pi, r3pi, rm, ro, rtanh, &
        & t1, t2, t3, vdum
    integer :: i, icartr, indexes, irctnt, jj, mder, msurf, nasurf, &
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
    common /potdcm/ vdum(10)
    common /atat/chh(10),chf(10)
    common /degsdrcm/dcorr2(3),dcorr3(3),dldr(3),de3c(3)
    common /engycm/ecor2b, ecor3b, elond, e3c
    common /dpc/dp1(3),dp2(3),dp3(3),ddp1(3),ddp2(3),ddp3(3)
    common /dist/rm(3),r0(3)
    common /kkp/ck1(10),ck2(10),ck3(10),ckp1(10),ckp2(10),ckp3(10)
    common /ref/ro(3)
    common /trial/eta1(10),eta2(10),eta3(10), &
        & etap1(10),etap2(10),etap3(10)
    common /vtanh/ rtanh(3,10),drtanh(3,10)
    dimension co3(3),dc3dr1(3),dc3dr2(3),dc3dr3(3)
    !
    jj = 0
    !
    do 10 i=6,10,2
        jj   = jj+1
        fi   = dble(i)
        r1pi = r(1)**i
        r2pi = r(2)**i
        r3pi = r(3)**i
        !
        !   Compute the auxiliary function gn(R)
        !
        g1 = 1.0d0+ck1(i)*exp(-ckp1(i)*(r(1)-ro(1)))
        g2 = 1.0d0+ck2(i)*exp(-ckp2(i)*(r(2)-ro(2)))
        g3 = 1.0d0+ck3(i)*exp(-ckp3(i)*(r(3)-ro(3)))
        !
        !   Compute the auxiliary function hn(R)
        !
        h1 = rtanh(1,i)**etap1(i)
        h2 = rtanh(2,i)**etap2(i)
        h3 = rtanh(3,i)**etap3(i)
        !
        t1 = chf(i)*(1.0d0-0.5d0*(g2*h3+g3*h2)) *dp1(jj)/r1pi
        t2 = chh(i)*(1.0d0-0.5d0*(g3*h1+g1*h3)) *dp2(jj)/r2pi
        t3 = chf(i)*(1.0d0-0.5d0*(g1*h2+g2*h1)) *dp3(jj)/r3pi
        !
        co3(jj) = t1+t2+t3
        !
        if (nder .eq. 1) then
            !
            !   Compute the partial derivatives w.r.t. R(i)
            !
            dg1 = -ckp1(i)*(g1-1.0d0)
            dg2 = -ckp2(i)*(g2-1.0d0)
            dg3 = -ckp3(i)*(g3-1.0d0)
            !
            dh1 = etap1(i)*(rtanh(1,i)**(etap1(i)-1))*drtanh(1,i)
            dh2 = etap2(i)*(rtanh(2,i)**(etap2(i)-1))*drtanh(2,i)
            dh3 = etap3(i)*(rtanh(3,i)**(etap3(i)-1))*drtanh(3,i)
            !
            dt1dr1 = t1*(ddp1(jj)/dp1(jj)-fi/r(1))
            dt1dr2 = -0.5d0*chf(i)*dp1(jj)*(dg2*h3+g3*dh2)/r1pi
            dt1dr3 = -0.5d0*chf(i)*dp1(jj)*(g2*dh3+dg3*h2)/r1pi
            !
            dt2dr1 = -0.5d0*chh(i)*dp2(jj)*(g3*dh1+dg1*h3)/r2pi
            dt2dr2 = t2*ddp2(jj)/dp2(jj)-t2*fi/r(2)
            dt2dr3 = -0.5d0*chh(i)*dp2(jj)*(dg3*h1+g1*dh3)/r2pi
            !
            dt3dr1 = -0.5d0*chf(i)*dp3(jj)*(dg1*h2+g2*dh1)/r3pi
            dt3dr2 = -0.5d0*chf(i)*dp3(jj)*(g1*dh2+dg2*h1)/r3pi
            dt3dr3 = t3*(ddp3(jj)/dp3(jj)-fi/r(3))
            !
            dc3dr1(jj) = dt1dr1+dt2dr1+dt3dr1
            dc3dr2(jj) = dt1dr2+dt2dr2+dt3dr2
            dc3dr3(jj) = dt1dr3+dt2dr3+dt3dr3
        endif
    10 continue
    !
    !   Compute the three-body dynamical correlation term by summing over al
    !   C6, C8, and C10 terms.
    !
    ecor3b = co3(1)+co3(2)+co3(3)
    !
    if (nder .eq. 1) then
        !
        !   Sum all the partial derivatives.
        !
        dcorr3(1) = dc3dr1(1)+dc3dr1(2)+dc3dr1(3)
        dcorr3(2) = dc3dr2(1)+dc3dr2(2)+dc3dr2(3)
        dcorr3(3) = dc3dr3(1)+dc3dr3(2)+dc3dr3(3)
    endif
    !
end subroutine corr3
!
!*****
!
!> @brief Evaluate VLOND for the FH2_GEN_6SEC_1993 potential.
subroutine vlond
    !
    !=======================================================================
    !   Compute the energy expressed by the lower eignenvalue of the London
    !=======================================================================
    !
    implicit none

    double precision :: a1, a2, a3, ad1, ad2, ad3, alf, alfp, alft, anuzero, arg, asc1, asc2, asc3, bd1, bd2, bd3, betp, &
        & bhf, bhh, c3c, cart, cartnu, chf, chh, csi, csip2, csip3, csip4, d, darg, dcorr2, dcorr3, dcsi, ddp1, ddp2, &
        & ddp3, de3c, dec, deesdr, degsdr, dehf, dehft, deijdr, descart, detrip, dex, dfcsi1, dfcsi2, dgscart, dijcart, &
        & dldr, dp1, dp2, dp3, dq, dr, dri, drp2, dscale, dsumq, dsw, dti, dtscle, dx, e3c, easyab, easyac, easybc, &
        & ecor2b, ecor3b, ehf, ehft, elond, engyes, engygs, engyij, etrip, ex, exch, ezero, fcsi, fcsi1, fcsi1t, fcsi2, &
        & fcsi2t, ff1, ff2, ff3, gamma, p, p3c, pengyes, pengygs, pengyij, q, q3c, r, r0, r1p2, r2p2, r2p3, r2p4, r3p2, &
        & rb, rm, rmt, scale, sw, texch, ti, tmp1, tmp1p2, tmp2, tscale, vdum, x, xhf, xhh, xterm, xterm1, xterm2, &
        & xterm3
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
    common /potdcm/ vdum(10)
    !
    common /atat/chh(10),chf(10)
    common /csi/fcsi,fcsi1,fcsi2
    common /dampc/ad1,ad2,ad3,bd1,bd2,bd3
    common /degsdrcm/dcorr2(3),dcorr3(3),dldr(3),de3c(3)
    common /diat/d,a1,a2,a3,gamma
    common /dist/rm,r0
    common /dpc/dp1(3),dp2(3),dp3(3),ddp1(3),ddp2(3),ddp3(3)
    common /engycm/ecor2b, ecor3b, elond, e3c
    common /local/c3c,alf,alfp,alft,p3c,q3c,betp
    common /mint/rmt
    common /par/bhh(5),bhf(5)
    common /swpar/ dec(3),rb(3),sw(3),dsw(3,3)
    common /tpar/xhh(3),xhf(3),p(4)
    !
    dimension rm(3),r0(3)
    dimension d(3),a1(3),a2(3),a3(3),gamma(3)
    dimension rmt(3)
    dimension dcsi(3),dfcsi1(3),dfcsi2(3),ehf(3),dehf(3)
    dimension ehft(3),dehft(3),dsumq(3),ti(3),dti(3,3)
    dimension etrip(3),detrip(3,3),dx(3),dri(3),x(3)
    dimension q(3),dq(3,3),ex(3),dex(3,3)
    dimension tscale(3), dtscle(3,3)
    !
    !   Initialize the arrays that will be used to store the partial
    !   derivatives w.r.t. R(i)
    !
    do 11 i = 1, 3
        do 12 j = 1, 3
            detrip(i,j) = 0.0d0
            dtscle(i,j) = 0.0d0
        12 continue
    11 continue
    !
    !   Compute the pseudo-angular coordinate
    !
    csi   = (r(1)-r(3))/r(2)
    csip2 = csi**2
    csip3 = csi**3
    csip4 = csi**4
    !
    !    Compute the F1 and F2 terms
    !
    fcsi1t = 1.0d0+p(1)*(1.0d0-csip2)+p(2)*(1.0d0-csip4)
    fcsi2t = 1.0d0+p(3)*(1.0d0-csip2)+p(4)*(1.0d0-csip4)
    !
    fcsi1 = fcsi1t**2
    fcsi2 = fcsi2t**2
    !
    if (nder .eq. 1) then
        !
        !   Compute the derivatives w.r.t. R of the terms computed to this point
        !
        dcsi(1) = 1.0d0/r(2)
        dcsi(2) = -csi/r(2)
        dcsi(3) = -1.0d0/r(2)
        !
        tmp1 = 2.0d0*fcsi1t*(-2.0d0*p(1)*csi-4.0d0*p(2)*csip3)
        tmp2 = 2.0d0*fcsi2t*(-2.0d0*p(3)*csi-4.0d0*p(4)*csip3)
        !
        do 10 i=1,3
            dfcsi1(i) = tmp1*dcsi(i)
            dfcsi2(i) = tmp2*dcsi(i)
        10 continue
    endif
    !
    r1p2 = r(1)**2
    r2p2 = r(2)**2
    r2p3 = r(2)**3
    r2p4 = r(2)**4
    r3p2 = r(3)**2
    !
    !   Compute the extended-Hartree-Fock curves and the derivatives w.r.t.
    !   for the singlet state of H2 and HF.
    !   Ref: A. J. C. Varandas and J. D. Silva, J. Chem. Soc. Faraday Trans.
    !        82, 593 (1986).
    !
    do 20 i=1,3
        dr     = r(i)-rm(i)
        drp2   = dr**2
        ehf(i) = -d(i)*(1.0d0+a1(i)*dr+a2(i)*drp2+a3(i)*drp2*dr) *exp(-gamma(i)*dr)
        if (nder .eq. 1) &
            & dehf(i) = -d(i)*(a1(i) + 2.0d0*a2(i)*dr +&
            &3.0d0*a3(i)*drp2)*exp(-gamma(i)*dr) -&
            &gamma(i)*ehf(i)
        !
    20 continue
    !
    !   Compute the triplet state curve and the derivative of the triplet st
    !   curve w.r.t. R(i) for the H2 diatomic
    !   Ref.: A. J. C. Varandas and J. Brando
    !         Mol. Phys. 45, 857 (1982)
    !
    scale = 1.0d0+xhh(1)*r(2)+xhh(2)*r2p2+xhh(3)*r2p3
    !
    arg = bhh(2)*r(2)+bhh(3)*r2p2+bhh(4)*r2p3+bhh(5)*r2p4
    !
    ehft(2) = bhh(1)*exp(-arg)/r(2)
    !
    etrip(2) = scale*ehft(2)
    !
    if (nder .eq. 1) then
        !
        !   Compute the derivatives of the H2 triplet curve w.r.t. R(i)
        !
        dscale   = xhh(1)+2.0d0*xhh(2)*r(2)+3.0d0*xhh(3)*r2p2
        darg     = bhh(2)+2.0d0*bhh(3)*r(2)+3.0d0*bhh(4)*r2p2+ 4.0d0*bhh(5)*r2p3
        dehft(2) = ehft(2)*(-darg-1.0d0/r(2))
        !
        detrip(2,2) = dscale*ehft(2)+scale*dehft(2)
        detrip(2,1) = 0.0d0
        detrip(2,3) = 0.0d0
    endif
    !
    !   Compute the triplet state curve and the derivatives of the triplet
    !   state curve w.r.t. R(i) for the HF diatomic.
    !
    do 30 i=1, 3, 2
        arg      = bhf(2)*r(i)+bhf(3)*r(i)**2
        ehft(i)  = bhf(1)*exp(-arg)/r(i)
        x(i)     = exp(-xhf(2)*r(i))
        etrip(i) = xhf(1)*(fcsi1*x(i)*x(i)+fcsi2*xhf(3)*x(i))/r(i)
    30 continue
    !
    if (nder .eq. 1) then
        !
        !   Compute the derivatives of the HF triplet state curve w.r.t. R(i)
        !
        do 50 i = 1, 3, 2
            darg     = bhf(2)+2.0d0*bhf(3)*r(i)
            dehft(i) = ehft(i)*(-1.0d0/r(i)-darg)
            dx(i)    = -xhf(2)*x(i)
            dri(i)   = 1.0d0
            do 40 j = 1, 3
                if (j .ne. i) then
                    dx(j)  = 0.0d0
                    dri(j) = 0.0d0
                endif
                detrip(i,j) = xhf(1) * (dfcsi1(j)*x(i)*x(i) + 2.0d0*fcsi1*x(i)*dx(j) + dfcsi2(j)*xhf(3)*x(i) + fcsi2* &
                    & xhf(3)*dx(j))/r(i) - etrip(i)*dri(j)/r(i)
            40 continue
        50 continue
    endif
    !
    !   Compute the scaling factor for the triplet state curves
    !
    do 51 i = 1, 3
        tscale(i) = 1.0d0
    51 continue
    asc1      = 1.0d0
    asc2      = 0.103d0
    tmp1      = 1.0d0/0.38d0
    asc3      = tmp1*tmp1*tmp1*tmp1
    tmp1      = (r(2) - 2.d0)
    tmp1p2    = tmp1 * tmp1
    tmp2      = exp(-asc3*tmp1p2*tmp1p2)
    tscale(2) = asc1 + asc2*tmp2
    if (nder .eq. 1) &
        & dtscle(2,2) = asc2*tmp2*(-asc3*4.0d0*tmp1p2*tmp1)
    !
    !   Compute the triplet state curves and combine the triplet and the
    !   singlet state curves to determine the eigenvalue of the
    !   London equation.
    !
    do 60 i = 1, 3
        !
        !   Compute the triplet state curves
        !
        ti(i) = (1.0d0-sw(i))*etrip(i)+sw(i)*ehft(i)
        !
        !   Scale the triplet state curves
        !
        ti(i) = ti(i)*tscale(i)
        !
        !   Compute the coulomb integrals
        !
        q(i) = 0.5d0*(ehf(i)+ti(i))
        !
        !   Compute the exchange integrals
        !
        ex(i) = 0.5d0*(ehf(i)-ti(i))
    60 continue
    !
    !  Diagnostics
    !
    vdum(1) = tscale(2)
    vdum(2) = ti(2)
    !
    ff1 = (ex(1) - ex(2)) ** 2
    ff2 = (ex(2) - ex(3)) ** 2
    ff3 = (ex(3) - ex(1)) ** 2
    !
    exch = sqrt(0.5d0*(ff1+ff2+ff3))
    !
    elond = q(1) + q(2) + q(3) - exch
    !
    if (nder .eq. 1) then
        !
        !   Compute the derivatives of the lowest eignenvalue of the London equa
        !   w.r.t. R(i)
        !
        do 70 i = 1, 3
            do 71 j = 1, 3
                dti(i,j) = dsw(i,j)*(ehft(i)-etrip(i)) + detrip(i,j)*(1.0d0 - sw(i))
                if (j .eq. i) dti(i,j) = dti(i,j) + sw(i)*dehft(i)
                dti(i,j) = dti(i,j)*tscale(i) + ti(i)*dtscle(i,j)/tscale(i)
            71 continue
        70 continue
        !
        do 80 i = 1, 3
            do 81 j = 1, 3
                if(j .eq. i)then
                    dq(i,j)  = 0.5d0*(dehf(j)+dti(i,j))
                    dex(i,j) = 0.5d0*(dehf(j)-dti(i,j))
                else
                    dq(i,j)  = 0.5d0*dti(i,j)
                    dex(i,j) = -dq(i,j)
                endif
            81 continue
        80 continue
        !
        dsumq(1) = dq(1,1)+dq(2,1)+dq(3,1)
        dsumq(2) = dq(1,2)+dq(2,2)+dq(3,2)
        dsumq(3) = dq(1,3)+dq(2,3)+dq(3,3)
        !
        xterm1 = 2.0d0*ex(1)-ex(2)-ex(3)
        xterm2 = 2.0d0*ex(2)-ex(1)-ex(3)
        xterm3 = 2.0d0*ex(3)-ex(1)-ex(2)
        !
        do 90 i = 1, 3
            texch = exch
            if (texch .eq. 0.0d0)texch = 1.0d0
            xterm   = 0.5d0*(xterm1*dex(1,i) + xterm2*dex(2,i)+ xterm3*dex(3,i))/texch
            dldr(i) = dsumq(i)-xterm
        90 continue
    endif
    !
end subroutine vlond
!
!****
!
!> @brief Evaluate THCENT for the FH2_GEN_6SEC_1993 potential.
subroutine thcent
    !
    !=======================================================================
    !   Compute the localized three-body term and the derivatives w.r.t. R(i
    !=======================================================================
    !
    implicit none

    double precision :: alf, alfp, alft, anuzero, betp, c3c, cart, cartnu, costh, dcorr2, dcorr3, dcos1, dcos2, dcos3, &
        & de3c, deesdr, degsdr, deijdr, descart, dgscart, dijcart, dldr, e3c, easyab, easyac, easybc, ecor2b, ecor3b, &
        & elond, engyes, engygs, engyij, ex1, ex2, ex3, ezero, f3, omp3c, omq3c, p3c, pengyes, pengygs, pengyij, q3c, r, &
        & r12, r1i, r22, r32, r3i, rdif, rdif2, rsum, t, t1, t2, t3, t4, t5, vdum
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
    common /potdcm/ vdum(10)
    common /degsdrcm/dcorr2(3),dcorr3(3),dldr(3),de3c(3)
    common /engycm/ecor2b, ecor3b, elond, e3c
    common /local/c3c,alf,alfp,alft,p3c,q3c,betp
    !
    omp3c = 1.d0 - p3c
    omq3c = 1.d0 - q3c
    r12   = r(1) * r(1)
    r22   = r(2) * r(2)
    r32   = r(3) * r(3)
    rsum  = r(1) + r(3)
    rdif  = r(1) - r(3)
    rdif2 = rdif * rdif
    t1    = -alfp * rsum
    ex1   = c3c * exp(t1 * rsum)
    t2    = -alf * rdif
    ex2   = omp3c * exp(t2 * rdif)
    t3    = -alft * rdif2
    ex3   = p3c * exp(t3 * rdif2)
    r1i   = 1.d0 / r(1)
    r3i   = 1.d0 / r(3)
    t4    = r1i*r3i
    !
    !   COSTH IS THE ANGLE BETWEEN R1 AND R3
    !
    costh = 0.5d0 * (r12 + r32 - r22) * t4
    t5    = omq3c * costh
    !
    !   SET MODIFICATION ...
    !
    t = q3c + t5 * costh
    !
    !   E3C CONTAINS THE THREE CENTER CORRECTION TO V
    !
    e3c = (ex2 + ex3) * t * ex1
    f3  = exp(-betp*(r(1)+r(3)-r(2))**2)
    e3c = e3c*f3
    if (nder .ne. 1) go to 280

    !   COMPUTE DERIVATIVE OF THE 3C TERM
    !   FIRST, DERIVATIVE OF COSTH
    !
    dcos1 = r3i - costh*r1i
    dcos2 = - r(2)*t4
    dcos3 = r1i - costh*r3i
    t4    = ex2 + ex3
    t5    = t4 * t5
    !
    !   NOW, DERIVATIVES OF E3C
    !
    de3c(1) = 2.d0*((t2 * ex2 + 2.d0 * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos1) * ex1*f3+ 2.0d0*e3c*(-betp*(r(1)+ &
        & r(3)-r(2)))
    de3c(2) = 2.d0 * t5 * dcos2 * ex1*f3+ 2.0d0*e3c*(betp*(r(1)+r(3)-r(2)))
    de3c(3) = 2.d0*((-t2 * ex2 - 2.d0 * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos3) * ex1*f3+ 2.0d0*e3c*(-betp*(r(1)+ &
        & r(3)-r(2)))
    280 continue
end subroutine thcent
!
!*****
!
!> @brief Evaluate XTANH for the FH2_GEN_6SEC_1993 potential.
subroutine xtanh( i, a )
    !
    !=======================================================================
    !   Compute the auxiliary function, hn, and the derivatives of hn w.r.t.
    !=======================================================================
    !
    implicit none

    double precision :: a, anuzero, cart, cartnu, deesdr, degsdr, deijdr, descart, dgscart, dijcart, drtanh, easyab, &
        & easyac, easybc, engyes, engygs, engyij, ex, ezero, pengyes, pengygs, pengyij, r, rtanh, rval, sgn, t, th, vdum
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
    common /potdcm/ vdum(10)
    common /vtanh/ rtanh(3,10),drtanh(3,10)
    dimension a(10)
    !
    do 10 j = 6, 10, 2
        rval = r(i) * a(j)
        if(rval .gt. 0.0d0) then
            ex  = exp( -2.0d0 * rval)
            sgn = 1.0d0
        else
            ex  = exp( 2.0d0 * rval)
            sgn = -1.0d0
        endif
        t  = 1.0d0 + ex
        th = 1.0d0/t
        if (nder .eq. 1) drtanh(i,j) = 4.0d0 * a(j) * ex * th * th
        rtanh(i,j) = sgn * th * (1.0d0 - ex)
    10 continue
end subroutine xtanh
!
!*****
!
!> @brief Evaluate SWITCH for the FH2_GEN_6SEC_1993 potential.
subroutine switch
    !
    !=======================================================================
    !   Compute the switching function, S, and the derivatives of S w.r.t. R
    !=======================================================================
    !
    implicit none

    double precision :: anuzero, cart, cartnu, dec, deesdr, degsdr, deijdr, descart, dgscart, dijcart, dsw, easyab, &
        & easyac, easybc, engyes, engygs, engyij, ex2u, exu, ezero, hysec, hytan, pengyes, pengygs, pengyij, r, rb, rmrb, &
        & rval, sgn, sw, temp1, temp2, vdum
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
    common /potdcm/ vdum(10)
    common /swpar/ dec(3),rb(3),sw(3),dsw(3,3)
    !
    !   Initialize the array that will be used to store the derivatives of
    !   the switching function w.r.t. R(i)
    !
    do 10 i = 1, 3
        do 10 j = 1, 3
            dsw(i,j) = 0.0d0
    10 continue
    !
    !   Compute the switching function.
    !
    do 20 i=1,3
        rmrb = r(i)-rb(i)
        rval = rmrb*dec(i)
        if (rval .gt. 0.0d0)then
            ex2u = exp(-2.0d0*rval)
            exu  = exp(-rval)
            sgn  = 1.0d0
        else
            ex2u = exp(2.0d0*rval)
            exu  = exp(rval)
            sgn  = -1.0d0
        endif
        temp1 = 1.0d0 + ex2u
        hytan = sgn*(1.0d0 - ex2u)/temp1
        sw(i) = 0.5d0*(1.0d0 + hytan)
    20 continue
    if (nder .ne. 1) go to 100
    do 30 i = 1, 3
        rmrb = r(i)-rb(i)
        rval = rmrb*dec(i)
        if (rval .gt. 0.0d0)then
            exu = exp(-rval)
            sgn = 1.0d0
        else
            exu = exp(rval)
            sgn = -1.0d0
        endif
        temp2    = 1.0d0 + exu*exu
        hysec    = sgn*2.0d0*exu/temp2
        dsw(i,i) = 0.5d0*dec(i)*hysec*hysec
    30 continue
    !
    100 continue
end subroutine switch
!
!*****
!=======================================================================
block data ptpacm
    !=======================================================================
    !     INPUT DATA FOR FH2
    !
    !=======================================================================
    implicit none

    double precision :: a1, a2, a3, alf, alfp, alft, alph0, alph1, anuzero, bet0, bet1, betp, bhf, bhh, c3c, cart, &
        & cartnu, chf, chh, ck1, ck2, ck3, ckp1, ckp2, ckp3, d, dec, disse, dsw, easyab, easyac, easybc, eta1, eta2, &
        & eta3, etap1, etap2, etap3, ezero, gamma, p, p3c, q3c, r0, rb, rm, rmt, ro, sw, xhf, xhh
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
    common /decm/ disse(3)
    common /atat/chh(10),chf(10)
    common /diat/d(3),a1(3),a2(3),a3(3),gamma(3)
    common /dist/rm(3),r0(3)
    common /dispc/alph0,alph1,bet0,bet1
    common /kkp/ck1(10),ck2(10),ck3(10), &
        & ckp1(10),ckp2(10),ckp3(10)
    common /local/c3c,alf,alfp,alft,p3c,q3c,betp
    common /mint/rmt(3)
    common /par/bhh(5),bhf(5)
    common /ref/ro(3)
    common /swpar/ dec(3),rb(3),sw(3),dsw(3,3)
    common /tpar/xhh(3),xhf(3),p(4)
    common /trial/eta1(10),eta2(10),eta3(10), &
        & etap1(10),etap2(10),etap3(10)
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
    data disse /0.224989333d0, 0.174472692d0, 0.224989333d0/
    !
    data xhf /2.633449d0, 0.417523d0, -0.090289d0/
    data xhh /-0.393828d0, 0.167810d0, -0.017081d0/
    data p  /-0.071547d0, -0.108550d0, 0.050461d0, -0.370255d0/
    data c3c, alfp /0.296668d0, 0.078606d0/
    !
    data alf,alft/0.18d0,2.14d0/
    data p3c,q3c,betp/0.95d0,1.0d0,0.18d0/
    data dec/2.0d0,2.0d0,2.0d0/
    data rb/6.80d0,5.00d0,6.80d0/
    data bhh/0.448467d0,-0.056687d0,0.246831d0,-0.018419d0,0.000598d0/
    data bhf/0.9149555d0,0.3222197d0,0.1273001d0,0.0d0,0.0d0/
    data chh(6),chh(8),chh(10)/6.499027d0,1.243991d2,3.2858d3/
    data chf(6),chf(8),chf(10)/6.68799d0,1.0289812d2,2.07391451d3/
    data rm(1),rm(2),rm(3)/1.7329d0,1.4010d0,1.7329d0/
    data r0(1),r0(2),r0(3)/5.9d0,6.9283d0,5.9d0/
    data ro(1),ro(2),ro(3)/1.7329d0,1.449d0,1.7329d0/
    data alph0,alph1,bet0,bet1/2.59528d1,1.1868d0,1.57381d1,9.729d-2/
    data eta1(6),eta1(8),eta1(10)/0.9438421608d0,0.9774440533d0, &
        & 0.9452240321d0/
    data eta2(6),eta2(8),eta2(10)/0.1085550150d1,0.1117269260d1, &
        & 0.1138536961d1/
    data eta3(6),eta3(8),eta3(10)/0.9438421608d0,0.9774440533d0, &
        & 0.9452240321d0/
    data etap1(6),etap1(8),etap1(10)/6.0d0,6.0d0,6.0d0/
    data etap2(6),etap2(8),etap2(10)/6.0d0,6.0d0,6.0d0/
    data etap3(6),etap3(8),etap3(10)/6.0d0,6.0d0,6.0d0/
    data ck1(6),ck1(8),ck1(10)/-0.2610767389d-1,-0.3428701964d-1, &
        & -0.4858536634d-1/
    data ck2(6),ck2(8),ck2(10)/-0.9709847270d-1,-0.1248186655d0, &
        & -0.1426680807d0/
    data ck3(6),ck3(8),ck3(10)/-0.2610767389d-1,-0.3428701964d-1, &
        & -0.4858536634d-1/
    data ckp1(6),ckp1(8),ckp1(10)/0.9438421608d0,0.9774440533d0, &
        & 0.9452240321d0/
    data ckp2(6),ckp2(8),ckp2(10)/0.1085550150d1,0.1117269260d1, &
        & 0.1138536961d1/
    data ckp3(6),ckp3(8),ckp3(10)/0.9438421608d0,0.9774440533d0, &
        & 0.9452240321d0/
    data d(1),a1(1),a2(1),a3(1),gamma(1)/0.19383609d0,2.4822787d0, &
        & 1.5435337d0,0.83093855d0,2.3992999d0/
    data d(2),a1(2),a2(2),a3(2),gamma(2)/0.15796326d0,2.1977034d0, &
        & 1.2932502d0,0.64375666d0,2.1835071d0/
    data d(3),a1(3),a2(3),a3(3),gamma(3)/0.19383609d0,2.4822787d0, &
        & 1.5435337d0,0.83093855d0,2.3992999d0/
    data rmt/1.0d5,1.0d5,1.0d5/
end
!=======================================================================
!> @brief Evaluate POTINFO for the FH2_GEN_6SEC_1993 potential.
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

!> @brief Evaluate ANCVRT for the FH2_GEN_6SEC_1993 potential.
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
!> @brief Evaluate CARTOU for the FH2_GEN_6SEC_1993 potential.
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

!> @brief Evaluate CARTTOR for the FH2_GEN_6SEC_1993 potential.
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

!> @brief Evaluate EUNITZERO for the FH2_GEN_6SEC_1993 potential.
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

!> @brief Evaluate RTOCART for the FH2_GEN_6SEC_1993 potential.
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

!> @brief Evaluate DEDCOU for the FH2_GEN_6SEC_1993 potential.
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
