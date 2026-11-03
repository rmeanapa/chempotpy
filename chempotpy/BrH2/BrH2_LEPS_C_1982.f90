!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!
!> @brief Evaluate pes for the BrH2_LEPS_C_1982 potential.
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

!> @brief Evaluate PREPOT for the BrH2_LEPS_C_1982 potential.
subroutine prepot
    !
    !   System:          BrH2
    !   Functional form: Diatomics-in-molecules plus three-center term
    !   Common name:     BrH2 DIM-3C
    !   Interface:       potlib2001
    !   References:      D. C. Clary
    !                    Chem. Phys. 71, 117 (1982)
    !                    I. Last and M. Baer
    !                    in Potential Energy Surfaces and Dynamics
    !                    edited by D. G. Truhlar (Plenum, New York, 1981)
    !                    p. 519
    !
    !   Number of bodies: 3
    !   Number of derivatives: 1
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
    !      Set in data statements
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
    !      NFLAG(18-20)
    !        passes the FORTRAN unit number used for output from the potenti
    !      /ASYCM/ EASYAB, EASYBC, EASYAC
    !        passes the energy in the three asymptotic valleys for an A + BC
    !        The energy in the AB valley, EASYAB, is equal to the energy of
    !        C atom "infinitely" far from the AB diatomic and R(AB) set equa
    !        Re(AB), the equilibrium bond length for the AB diatomic.
    !        In this potential the AB valley represents H infinitely far fro
    !        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the t
    !        EASYBC and EASYAC represent the energies in the H2 and the HBr
    !        valleys, respectively.
    !
    !   Default Parameter Values:
    !      Variable      Default value
    !      NFLAG(18)        6
    !
    !*****
    !
    implicit none

    double precision :: ahh, ahx, alfw, anuzero, bhh, bhx, c3, cart, cartnu, denom, dh1, dh2, dh3, dhh, dhx, e, easyab, &
        & easyac, easybc, eps, eta1p, eta3p, eta3s, etahh, ezero, g, h, one, rhh, rhx, rt34, rt38, scr1, scr2, t, u, xib, &
        & xih, xix
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
    common /precm/ h(10),dh1(10),dh2(10),dh3(10),u(4,4), &
        & e(4),scr1(4),scr2(4), &
        & xib,c3,denom,rt34,rt38
    common /potcm/ eps,one,rhh,ahh,dhh,bhh,xih,xix,g,alfw,rhx,ahx, &
        & dhx,bhx,etahh,eta3s,eta1p,eta3p
    !
    !      LOGICAL LCOL
    !      DIMENSION H(10),DH1(10),DH2(10),DH3(10),U(4,4),E(4),SCR1(4),
    !     1          SCR2(4)
    !      COMMON /POTCOM/ EPS,ONE,RHH,AHH,DHH,BHH,XIH,XIX,G,ALFW,RHX,AHX,
    !     &                DHX,BHX,ETAHH,ETA3S,ETA1P,ETA3P
    !
    !   Echo potential parameters to UNIT NFLAG(18)
    !
    if(natoms.gt.25) then
        write(nflag(18),'(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    end if
    !
    easyab = dhx
    easybc = dhh
    easyac = dhx
    !
    xib   = 0.5d0*(xih + xix)
    t     = xib*rhx
    c3    = 1.d0/3.d0
    denom = 2.d0*rhx*(1.d0+t*(1.d0+c3*t))*exp(-t)
    g     = g/rhx
    rt34  = 0.25d0*sqrt(3.d0)
    rt38  = 0.5d0*rt34
    !
    !
    !
    ezero(1) = dhh
    !
    do i=1,5
        ref(i) = ' '
    end do
    !
    ref(1) = 'D. C. Clary, Chem. Phys. 71, 117 (1982)'
    ref(2) = 'I. Last and M. Baer in'
    ref(3) = 'Potential Energy Surfaces and Dynamics'
    ref(4) = 'D. G. Truhlar (Ed.) (Plenum, New York, 1981), p. 519'
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
!> @brief Evaluate POT for the BrH2_LEPS_C_1982 potential.
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

    double precision :: ahh, ahx, alfw, anuzero, bhh, bhx, c3, cart, cartnu, csa, csa2, csb, csb2, csg2, d1, d2, d3, &
        & dcsa21, dcsa22, dcsa23, dcsb21, dcsb22, dcsb23, dcsg21, dcsg22, dcsg23, deesdr, degsdr, deijdr, denom, descart, &
        & dgscart, dh1, dh2, dh3, dhh, dhx, dijcart, dp11, dp12, dp21, dp22, dp31, dp32, ds11, ds12, ds21, ds22, ds31, &
        & ds32, dshh, dshx1, dshx2, dsn2a1, dsn2a2, dsn2a3, dsn2b1, dsn2b2, dsn2b3, dv1hh, dv1p1, dv1p2, dv1s1, dv1s2, &
        & dv3hh, dv3p1, dv3p2, dv3s1, dv3s2, e, easyab, easyac, easybc, engyes, engygs, engyij, eps, eta1p, eta3p, eta3s, &
        & etahh, ex1, ex2, ezero, g, h, one, p11, p12, p21, p22, p31, p32, pengyes, pengygs, pengyij, r, r1, r12, r13, &
        & r1s, r2, r23, r2s, r3, r3s, rdif, rdif2, rhh, rhx, rt34, rt38, s11, s12, s21, s22, s31, s32, scr1, scr2, shh, &
        & shx1, shx2, sn2a, sn2b, sna, sna2, snb, snb2, t, t1, t11, t12, t2, t21, t22, t3, t31, t32, u, v1hh, v1p1, v1p2, &
        & v1s1, v1s2, v3hh, v3p1, v3p2, v3s1, v3s2, xib, xih, xix
    integer :: icartr, ierr, index, indexes, irctnt, k, kk, l, ll, mder, msurf, nasurf, natoms, nder, nflag, nulbl
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
    logical :: lcol
    common /precm/ h(10),dh1(10),dh2(10),dh3(10),u(4,4), &
        & e(4),scr1(4),scr2(4), &
        & xib,c3,denom,rt34,rt38
    common /potcm/ eps,one,rhh,ahh,dhh,bhh,xih,xix,g,alfw,rhx,ahx, &
        & dhx,bhx,etahh,eta3s,eta1p,eta3p
    !
    call cartou
    call carttor
    !
    !   The potential routines define R1=R(Br-H), R2=R(H-Br), and R3=R(H-H);
    !   rearrange the input coordinates to match these definitions
    !
    r1 = r(1)
    r2 = r(3)
    r3 = r(2)
    !
    if(r1.gt.r2+r3) r1 = r2+r3
    if(r2.gt.r1+r3) r2 = r1+r3
    if(r3.gt.r1+r2) r3 = r1+r2
    r1s = r1*r1
    r2s = r2*r2
    r3s = r3*r3
    r12 = r1*r2
    r13 = r1*r3
    r23 = r2*r3
    t3  = (r1s+r2s-r3s)*0.5d0
    t2  = (r1s+r3s-r2s)*0.5d0
    t1  = (r2s+r3s-r1s)*0.5d0
    csa = t2/r13
    if(abs(csa) .gt. one) csa = sign(one,csa)
    csa2 = csa*csa
    sna2 = 1.d0-csa2
    lcol = sna2 .lt. eps
    sna  = sqrt(sna2)
    t22  = 2.d0*csa
    t11  = 0.d0
    if(.not.lcol) t11 = 2.d0*(sna2-csa2)/sna
    sn2a   = t22*sna
    t      = t3/(r1*r13)
    dcsa21 = t22*t
    dsn2a1 = t11*t
    t      = -r2/r13
    dcsa22 = t22*t
    dsn2a2 = t11*t
    t      = t1/(r13*r3)
    dcsa23 = t22*t
    dsn2a3 = t11*t
    csb    = t1/r23
    if(abs(csb).gt.one) csb = sign(one,csb)
    csb2 = csb*csb
    snb2 = 1.d0-csb2
    snb  = sqrt(snb2)
    t11  = 0.d0
    if(.not.lcol) t11 = 2.d0*(snb2-csb2)/snb
    t22    = 2.d0*csb
    sn2b   = t22*snb
    t      = -r1/r23
    dcsb21 = t22*t
    dsn2b1 = t11*t
    t      = t3/(r2*r23)
    dcsb22 = t22*t
    dsn2b2 = t11*t
    t      = t2/(r23*r3)
    dcsb23 = t22*t
    dsn2b3 = t11*t
    t      = t3/r12
    csg2   = t*t
    t      = 2.d0*t
    dcsg21 = t*t2/(r1*r12)
    dcsg22 = t*t1/(r12*r2)
    dcsg23 = -t*r3/r12
    !
    !  DIATOMIC CURVES
    !
    !    HH
    !
    rdif  = r3 - rhh
    ex1   = exp(-ahh*rdif)
    rdif2 = rdif*rdif
    ex2   = exp(-bhh*rdif*rdif2)
    t1    = dhh*ex1*ex2
    v1hh  = t1*(ex1-2.d0)
    v3hh  = etahh*t1*(ex1+2.d0)
    t1    = 2.d0*ahh*t1
    t2    = 3.d0*bhh*rdif2
    dv1hh = t1*(1.d0-ex1) - t2*v1hh
    dv3hh = -t1*etahh*(1.d0+ex1) - t2*v3hh
    !
    !    HX
    !
    call vhx(r1,v1s1,v3s1,v1p1,v3p1,dv1s1,dv3s1,dv1p1,dv3p1)
    call vhx(r2,v1s2,v3s2,v1p2,v3p2,dv1s2,dv3s2,dv1p2,dv3p2)
    !
    s11 = v1s1 + 3.d0*v3s1
    s12 = v1s2 + 3.d0*v3s2
    s21 = 3.d0*v1s1 + v3s1
    s22 = 3.d0*v1s2 + v3s2
    s31 = v1s1 - v3s1
    s32 = v1s2 - v3s2
    p11 = v1p1 + 3.d0*v3p1
    p12 = v1p2 + 3.d0*v3p2
    p21 = 3.d0*v1p1 + v3p1
    p22 = 3.d0*v1p2 + v3p2
    p31 = v1p1 - v3p1
    p32 = v1p2 - v3p2
    !
    ds11 = dv1s1 + 3.d0*dv3s1
    ds12 = dv1s2 + 3.d0*dv3s2
    ds21 = 3.d0*dv1s1 + dv3s1
    ds22 = 3.d0*dv1s2 + dv3s2
    ds31 = dv1s1 - dv3s1
    ds32 = dv1s2 - dv3s2
    dp11 = dv1p1 + 3.d0*dv3p1
    dp12 = dv1p2 + 3.d0*dv3p2
    dp21 = 3.d0*dv1p1 + dv3p1
    dp22 = 3.d0*dv1p2 + dv3p2
    dp31 = dv1p1 - dv3p1
    dp32 = dv1p2 - dv3p2
    !
    !    CONSTRUCT 4X4 HAMILTONIAN PUT IN PACKED ARRAY
    !
    h(1)               = v1hh + 0.25d0*(s11*csa2 + p11*sna2 + s12*csb2 + p12*snb2)
    h(3)               = v1hh + 0.25d0*(s11*sna2 + p11*csa2 + s12*snb2 + p12*csb2)
    h(6)               = v3hh + 0.25d0*(s21*csa2 + p21*sna2 + s22*csb2 + p22*snb2)
    h(10)              = v3hh + 0.25d0*(s21*sna2 + p21*csa2 + s22*snb2 + p22*csb2)
    t11                = s11 - p11
    t12                = s12 - p12
    h(2)               = 0.d0
    if(.not.lcol) h(2) = 0.125d0*(t11*sn2a - t12*sn2b)
    h(4)               = rt34*(s31*csa2 + p31*sna2 - s32*csb2 - p32*snb2)
    t31                = s31 - p31
    t32                = s32 - p32
    h(5)               = 0.d0
    if(.not.lcol) h(5) = rt38*(t31*sn2a + t32*sn2b)
    h(7)               = h(5)
    h(8)               = rt34*(s31*sna2 + p31*csa2 - s32*snb2 - p32*csb2)
    t21                = s21 - p21
    t22                = s22 - p22
    h(9)               = 0.d0
    if(.not.lcol) h(9) = 0.125d0*(t21*sn2a - t22*sn2b)
    !
    !   CONSTRUCT DERIVATIVE OF HAMILTONIAN MATRIX
    !
    t       = t11*dcsa21 + t12*dcsb21
    dh1(1)  = 0.25d0*(ds11*csa2 + dp11*sna2 + t)
    dh1(2)  = 0.125d0*((ds11-dp11)*sn2a + t11*dsn2a1 - t12*dsn2b1)
    dh1(3)  = 0.25d0*(ds11*sna2 + dp11*csa2 - t)
    t       = t11*dcsa22 + t12*dcsb22
    dh2(1)  = 0.25d0*(ds12*csb2 + dp12*snb2 + t)
    dh2(2)  = 0.125d0*(-(ds12-dp12)*sn2b + t11*dsn2a2 - t12*dsn2b2)
    dh2(3)  = 0.25d0*(ds12*snb2 + dp12*csb2 - t)
    t       = 0.25d0*(t11*dcsa23 + t12*dcsb23)
    dh3(1)  = dv1hh + t
    dh3(2)  = 0.125d0*(t11*dsn2a3 - t12*dsn2b3)
    dh3(3)  = dv1hh - t
    t       = t21*dcsa21 + t22*dcsb21
    dh1(6)  = 0.25d0*(ds21*csa2 + dp21*sna2 + t)
    dh1(9)  = 0.125d0*((ds21-dp21)*sn2a + t21*dsn2a1 - t22*dsn2b1)
    dh1(10) = 0.25d0*(ds21*sna2 + dp21*csa2 - t)
    t       = t21*dcsa22 + t22*dcsb22
    dh2(6)  = 0.25d0*(ds22*csb2 + dp22*snb2 + t)
    dh2(9)  = 0.125d0*((dp22-ds22)*sn2b + t21*dsn2a2 - t22*dsn2b2)
    dh2(10) = 0.25d0*(ds22*snb2 + dp22*csb2 - t)
    t       = 0.25d0*(t21*dcsa23 + t22*dcsb23)
    dh3(6)  = dv3hh + t
    dh3(9)  = 0.125d0*(t21*dsn2a3 - t22*dsn2b3)
    dh3(10) = dv3hh - t
    t       = t31*dcsa21 - t32*dcsb21
    dh1(4)  = rt34*(ds31*csa2 + dp31*sna2 + t)
    dh1(5)  = rt38*((ds31-dp31)*sn2a + t31*dsn2a1 + t32*dsn2b1)
    dh1(7)  = dh1(5)
    dh1(8)  = rt34*(ds31*sna2 + dp31*csa2 - t)
    t       = t31*dcsa22 - t32*dcsb22
    dh2(4)  = rt34*(-ds32*csb2 - dp32*snb2 + t)
    dh2(5)  = rt38*((ds32-dp32)*sn2b + t31*dsn2a2 + t32*dsn2b2)
    dh2(7)  = dh2(5)
    dh2(8)  = rt34*(-ds32*snb2 - dp32*csb2 - t)
    t       = t31*dcsa23 - t32*dcsb23
    dh3(4)  = rt34*t
    dh3(5)  = rt38*(t31*dsn2a3 + t32*dsn2b3)
    dh3(7)  = dh3(5)
    dh3(8)  = -dh3(4)
    !
    !    DIAGONALIZE H PUT EIGENVECTORS INTO U, EIGENVALUES INTO E
    !    USE EISPACK ROUTINE RSP
    !
    call rsp(4,4,10,h,e,1,u,scr1,scr2,ierr)
    if(ierr.ne.0) then
        write(nflag(18),'(/,2X,''Error returned by the subprogram RSP in POT'', /,2X,''R = '',T10,3(1PE13.5,1X), /,2X' // &
            ',''H = '',(T10,4(1PE13.5,1X)))') r1,r2,r3,h
        write(nflag(18),'(2X,''V1HH,V3HH,V1S1,V3S1,V1S2,V3S2,V1P1,V3P1,V1P2,V3P2='', /,2X,(T10,4(1PE13.5,1X)))') v1hh,v3hh,v1s1,v3s1,v1s2, &
            & v3s2,v1p1,v3p1,v1p2,v3p2
        stop 'POT 3'
    end if
    !
    !    EVALUATE DERIVATIVES OF LOWEST ROOT
    !
    d1 = 0.d0
    d2 = 0.d0
    d3 = 0.d0
    do 20 l = 1,4
        t1 = u(l,1)
        ll = (l*(l-1))/2
        do 20 k = 1,4
            t2 = t1*u(k,1)
            kk = (k*(k-1))/2
            if(l.ge.k) index = ll + k
            if(k.gt.l) index = kk + l
            d1 = d1 + t2*dh1(index)
            d2 = d2 + t2*dh2(index)
            d3 = d3 + t2*dh3(index)
    20 continue
    !
    !   SUPPLEMENTARY TERM
    !
    t      = xih*r3
    ex1    = exp(-t)
    shh    = (1.d0+t*(1.d0+c3*t))*ex1
    dshh   = -xih*t*(1.d0+t)*c3*ex1
    t      = xib*r1
    ex1    = exp(-t)
    shx1   = r1*(1.d0+t*(1.d0+c3*t))*ex1/denom
    dshx1  = (1.d0+t*(1.d0-c3*t*t))*ex1/denom
    t      = xib*r2
    ex1    = exp(-t)
    shx2   = r2*(1.d0+t*(1.d0+c3*t))*ex1/denom
    dshx2  = (1.d0+t*(1.d0-c3*t*t))*ex1/denom
    rdif   = r1 - r2
    t1     = g*exp(-alfw*rdif*rdif)
    t2     = shh*(shx1+shx2) + shx1*shx2
    engygs = e(1) + t2*t1*csg2 + ezero(1)
    !
    call eunitzero
    if(nder.ne.0) then
        t3        = 2.d0*alfw*rdif*csg2
        degsdr(1) = d1 + (dshx1*(shh+shx2)*csg2 + t2*(dcsg21-t3))*t1
        degsdr(3) = d2 + (dshx2*(shh+shx1)*csg2 + t2*(dcsg22+t3))*t1
        degsdr(2) = d3 + (dshh*(shx1+shx2)*csg2 + t2*dcsg23)*t1
        call rtocart
        call dedcou
    endif
    !
    !
end subroutine pot
!
!
!> @brief Evaluate VHX for the BrH2_LEPS_C_1982 potential.
subroutine vhx( rr, v1s, v3s, v1p, v3p, dv1s, dv3s, dv1p, dv3p )
    implicit none

    double precision :: ahh, ahx, alfw, anuzero, bhh, bhx, cart, cartnu, deesdr, degsdr, deijdr, descart, dgscart, dhh, &
        & dhx, dijcart, dv1, dv1p, dv1s, dv3p, dv3s, easyab, easyac, easybc, engyes, engygs, engyij, eps, eta1p, eta3p, &
        & eta3s, etahh, ex1, ex2, ezero, g, one, pengyes, pengygs, pengyij, r, rdif, rdif2, rhh, rhx, rr, t1, t2, v1, &
        & v1p, v1s, v3p, v3s, xih, xix
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
    common /potcm/ eps,one,rhh,ahh,dhh,bhh,xih,xix,g,alfw,rhx,ahx, &
        & dhx,bhx,etahh,eta3s,eta1p,eta3p
    rdif  = rr - rhx
    rdif2 = rdif*rdif
    ex1   = exp(-ahx*rdif)
    ex2   = exp(-bhx*rdif2*rdif)
    t1    = dhx*ex1*ex2
    v1s   = t1*(ex1-2.d0)
    v1    = t1*(ex1+2.d0)
    v3s   = eta3s*v1
    v1p   = eta1p*v1
    v3p   = eta3p*v1
    t1    = 2.d0*ahx*t1
    t2    = 3.d0*bhx*rdif2
    dv1s  = t1*(1.d0-ex1) - t2*v1s
    dv1   = -t1*(1.d0+ex1) - t2*v1
    dv3s  = eta3s*dv1
    dv1p  = eta1p*dv1
    dv3p  = eta3p*dv1
end subroutine vhx
!
block data ptpacm
    implicit none

    double precision :: ahh, ahx, alfw, anuzero, bhh, bhx, cart, cartnu, dhh, dhx, easyab, easyac, easybc, eps, eta1p, &
        & eta3p, eta3s, etahh, ezero, g, one, rhh, rhx, xih, xix
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
    common /potcm/ eps,one,rhh,ahh,dhh,bhh,xih,xix,g,alfw,rhx,ahx, &
        & dhx,bhx,etahh,eta3s,eta1p,eta3p
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
    data eps /1.d-6/, one /1.0d0/
    data rhh,ahh,dhh,bhh /1.4016d0, 1.0291d0, 0.17447d0, 0.018d0/
    data xih,xix,g,alfw /1.d0, 1.6d0, 0.22d0, 1.0d0/
    !
    data rhx /2.673d0/
    data ahx /0.957d0/
    data dhx /0.1439d0/
    data bhx /0.012d0/
    data etahh /0.393764d0/
    data eta3s /0.322d0/
    data eta1p /0.20286d0/
    data eta3p /0.1771d0/
    !
end
!
!> @brief Evaluate RSP for the BrH2_LEPS_C_1982 potential.
subroutine rsp( nm, n, nv, a, w, matz, z, fv1, fv2, ierr )
    implicit none

    double precision :: anuzero, cart, cartnu, deesdr, degsdr, deijdr, descart, dgscart, dijcart, easyab, easyac, easybc, &
        & engyes, engygs, engyij, ezero, one, pengyes, pengygs, pengyij, r, zero
    integer :: i, icartr, ierr, indexes, irctnt, j, matz, mder, msurf, n, nasurf, natoms, nder, nflag, nm, nulbl, nv
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
    double precision :: a(nv),w(n),z(nm,n),fv1(n),fv2(n)
    !
    !     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
    !     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
    !     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
    !     OF A REAL SYMMETRIC PACKED MATRIX.
    !
    !     ON INPUT-
    !
    !        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
    !        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
    !        DIMENSION STATEMENT,
    !
    !        N  IS THE ORDER OF THE MATRIX  A,
    !
    !        NV  IS AN INTEGER VARIABLE SET EQUAL TO THE
    !        DIMENSION OF THE ARRAY  A  AS SPECIFIED FOR
    !        A  IN THE CALLING PROGRAM.  NV  MUST NOT BE
    !        LESS THAN  N*(N+1)/2,
    !
    !        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
    !        PACKED MATRIX STORED ROW-WISE,
    !
    !        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
    !        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO
    !        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
    !
    !     ON OUTPUT-
    !
    !        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,
    !
    !        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,
    !
    !        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN
    !        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE
    !        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,
    !
    !        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
    !
    !     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
    !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
    !
    !     ------------------------------------------------------------------
    !
    data zero,one/0.0d0,1.0d0/
    if (n .le. nm) go to 5
    ierr = 10 * n
    go to 50
    5 if (nv .ge. (n * (n + 1)) / 2) go to 10
    ierr = 20 * n
    go to 50
    !
    10 call  tred3(n,nv,a,w,fv1,fv2)
    if (matz .ne. 0) go to 20
    !
    !     ********** FIND EIGENVALUES ONLY **********
    !
    call  tqlrat(n,w,fv2,ierr)
    go to 50
    !
    !     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
    !
    20 do 40 i = 1, n
        !
        do 30 j = 1, n
            z(j,i) = zero
        30 continue
        !
        z(i,i) = one
    40 continue
    !
    call  tql2(nm,n,w,fv1,z,ierr)
    if (ierr .ne. 0) go to 50
    call  trbak3(nm,n,nv,a,n,z)
    50 return
    !
    !     ********** LAST CARD OF RSP **********
    !
end subroutine rsp
!
!> @brief Evaluate TRED3 for the BrH2_LEPS_C_1982 potential.
subroutine tred3( n, nv, a, d, e, e2 )
    !
    implicit none

    double precision :: anuzero, cart, cartnu, deesdr, degsdr, deijdr, descart, dgscart, dijcart, easyab, easyac, easybc, &
        & engyes, engygs, engyij, ezero, f, g, h, hh, pengyes, pengygs, pengyij, r, scale, zero
    integer :: i, icartr, ii, indexes, irctnt, iz, j, jk, k, l, mder, msurf, n, nasurf, natoms, nder, nflag, nulbl, nv
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
    double precision :: a(nv),d(n),e(n),e2(n)
    !
    !     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
    !     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
    !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
    !
    !     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
    !     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
    !     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
    !
    !     ON INPUT-
    !
    !        N IS THE ORDER OF THE MATRIX,
    !
    !        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
    !          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
    !
    !        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
    !          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
    !          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
    !
    !     ON OUTPUT-
    !
    !        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
    !          TRANSFORMATIONS USED IN THE REDUCTION,
    !
    !        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
    !
    !        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
    !          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
    !
    !        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
    !          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
    !
    !     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
    !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
    !
    !     ------------------------------------------------------------------
    !
    !     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
    !
    data zero/0.0d0/
    do  300 ii = 1, n
        i     = n + 1 - ii
        l     = i - 1
        iz    = (i * l) / 2
        h     = zero
        scale = zero
        if (l .lt. 1) go to 130
        !
        !     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
        !
        do 120 k = 1, l
            iz    = iz + 1
            d(k)  = a(iz)
            scale = scale + abs(d(k))
        120 continue
        !
        if (scale .ne. zero) go to 140
        130 e(i) = zero
        e2(i) = zero
        go to 290
        !
        140 do 150 k = 1, l
            d(k) = d(k) / scale
            h    = h + d(k) * d(k)
        150 continue
        !
        e2(i) = scale * scale * h
        f     = d(l)
        g     = -sign(sqrt(h),f)
        e(i)  = scale * g
        h     = h - f * g
        d(l)  = f - g
        a(iz) = scale * d(l)
        if (l .eq. 1) go to 290
        f = zero
        !
        do 240 j = 1, l
            g  = zero
            jk = (j * (j-1)) / 2
            !
            !     ********** FORM ELEMENT OF A*U **********
            !
            do 180 k = 1, l
                jk = jk + 1
                if (k .gt. j) jk = jk + k - 2
                g = g + a(jk) * d(k)
            180 continue
            !
            !     ********** FORM ELEMENT OF P **********
            !
            e(j) = g / h
            f    = f + e(j) * d(j)
        240 continue
        !
        hh = f / (h + h)
        jk = 0
        !
        !     ********** FORM REDUCED A **********
        !
        do 260 j = 1, l
            f    = d(j)
            g    = e(j) - hh * f
            e(j) = g
            !
            do 260 k = 1, j
                jk    = jk + 1
                a(jk) = a(jk) - f * e(k) - g * d(k)
        260 continue
        !
        290 d(i) = a(iz+1)
        a(iz+1) = scale * sqrt(h)
    300 continue
    !
    !
    !     ********** LAST CARD OF TRED3 **********
    !
end subroutine tred3
!
!> @brief Evaluate TQLRAT for the BrH2_LEPS_C_1982 potential.
subroutine tqlrat( n, d, e2, ierr )
    !
    implicit none

    double precision :: anuzero, b, c, cart, cartnu, d, deesdr, degsdr, deijdr, descart, dgscart, dijcart, e2, easyab, &
        & easyac, easybc, engyes, engygs, engyij, ezero, f, g, h, one, p, pengyes, pengygs, pengyij, r, rr, s, two, zero
    integer :: i, icartr, ierr, ii, indexes, irctnt, j, l, l1, m, mder, mml, msurf, n, nasurf, natoms, nder, nflag, &
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
    dimension d(n),e2(n)
    double precision :: machep
    !
    !     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
    !     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
    !
    !     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
    !     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
    !
    !     ON INPUT-
    !
    !        N IS THE ORDER OF THE MATRIX,
    !
    !        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
    !
    !        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
    !          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
    !
    !      ON OUTPUT-
    !
    !        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
    !          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
    !          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
    !          THE SMALLEST EIGENVALUES,
    !
    !        E2 HAS BEEN DESTROYED,
    !
    !        IERR IS SET TO
    !          ZERO       FOR NORMAL RETURN,
    !          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
    !                     DETERMINED AFTER 30 ITERATIONS.
    !
    !     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
    !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
    !
    !     ------------------------------------------------------------------
    !
    !                **********
    !
    data zero,one,two/0.0d0,1.0d0,2.0d0/
    !
    !     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
    !                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
    !
    machep = two**(-37)
    !
    ierr = 0
    if (n .eq. 1) go to 1001
    !
    do 100 i = 2, n
    100 e2(i-1) = e2(i)
    !
    f     = zero
    b     = zero
    e2(n) = zero
    !
    do 290 l = 1, n
        j = 0
        h = machep * (abs(d(l)) + sqrt(e2(l)))
        if (b .gt. h) go to 105
        b = h
        c = b * b
        !
        !     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
        !
        105 do 110 m = l, n
            if (e2(m) .le. c) go to 120
            !
            !     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
            !                THROUGH THE BOTTOM OF THE LOOP **********
            !
        110 continue
        !
        120 if (m .eq. l) go to 210
        130 if (j .eq. 30) go to 1000
        j = j + 1
        !
        !     ********** FORM SHIFT **********
        !
        l1   = l + 1
        s    = sqrt(e2(l))
        g    = d(l)
        p    = (d(l1) - g) / (two * s)
        rr   = sqrt(p*p+one)
        d(l) = s / (p + sign(rr,p))
        h    = g - d(l)
        !
        do 140 i = l1, n
        140 d(i) = d(i) - h
        !
        f = f + h
        !
        !     ********** RATIONAL QL TRANSFORMATION **********
        !
        g = d(m)
        if (g .eq. zero) g = b
        h   = g
        s   = zero
        mml = m - l
        !
        !     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
        !
        do 200 ii = 1, mml
            i       = m - ii
            p       = g * h
            rr      = p + e2(i)
            e2(i+1) = s * rr
            s       = e2(i) / rr
            d(i+1)  = h + s * (h + d(i))
            g       = d(i) - e2(i) / g
            if (g .eq. zero) g = b
            h = g * p / rr
        200 continue
        !
        e2(l) = s * g
        d(l)  = h
        !
        !     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********
        !
        if (h .eq. zero) go to 210
        if (abs(e2(l)) .le. abs(c/h)) go to 210
        e2(l) = h * e2(l)
        if (e2(l) .ne. zero) go to 130
        210 p = d(l) + f
        !
        !     ********** ORDER EIGENVALUES **********
        !
        if (l .eq. 1) go to 250
        !
        !     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
        !
        do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
        230 continue
        !
        250 i = 1
        270 d(i) = p
    290 continue
    !
    go to 1001
    !
    !     ********** SET ERROR -- NO CONVERGENCE TO AN
    !                EIGENVALUE AFTER 30 ITERATIONS **********
    !
    1000 ierr = l
    1001 return
    !
    !     ********** LAST CARD OF TQLRAT **********
    !
end subroutine tqlrat
!
!> @brief Evaluate TQL2 for the BrH2_LEPS_C_1982 potential.
subroutine tql2( nm, n, d, e, z, ierr )
    !
    implicit none

    double precision :: anuzero, b, c, cart, cartnu, d, deesdr, degsdr, deijdr, descart, dgscart, dijcart, e, easyab, &
        & easyac, easybc, engyes, engygs, engyij, ezero, f, g, h, one, p, pengyes, pengygs, pengyij, r, rr, s, two, z, &
        & zero
    integer :: i, icartr, ierr, ii, indexes, irctnt, j, k, l, l1, m, mder, mml, msurf, n, nasurf, natoms, nder, nflag, &
        & nm, nulbl
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
    dimension d(n),e(n),z(nm,n)
    double precision :: machep
    !
    !     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
    !     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
    !     WILKINSON.
    !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
    !
    !     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
    !     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
    !     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
    !     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
    !     FULL MATRIX TO TRIDIAGONAL FORM.
    !
    !     ON INPUT-
    !
    !        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
    !          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
    !          DIMENSION STATEMENT,
    !
    !        N IS THE ORDER OF THE MATRIX,
    !
    !        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
    !
    !        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
    !          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
    !
    !        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
    !          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
    !          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
    !          THE IDENTITY MATRIX.
    !
    !      ON OUTPUT-
    !
    !        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
    !          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
    !          UNORDERED FOR INDICES 1,2,...,IERR-1,
    !
    !        E HAS BEEN DESTROYED,
    !
    !        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
    !          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
    !          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
    !          EIGENVALUES,
    !
    !        IERR IS SET TO
    !          ZERO       FOR NORMAL RETURN,
    !          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
    !                     DETERMINED AFTER 30 ITERATIONS.
    !
    !     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
    !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
    !
    !     ------------------------------------------------------------------
    !
    !                **********
    !
    data zero,one,two/0.0d0,1.0d0,2.0d0/
    !
    !     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
    !                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
    !
    machep = two**(-37)
    !
    ierr = 0
    if (n .eq. 1) go to 1001
    !
    do 100 i = 2, n
    100 e(i-1) = e(i)
    !
    f    = zero
    b    = zero
    e(n) = zero
    !
    do 240 l = 1, n
        j = 0
        h = machep * (abs(d(l)) + abs(e(l)))
        if (b .lt. h) b = h
        !
        !     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
        !
        do 110 m = l, n
            if (abs(e(m)) .le. b) go to 120
            !
            !     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
            !                THROUGH THE BOTTOM OF THE LOOP **********
            !
        110 continue
        !
        120 if (m .eq. l) go to 220
        130 if (j .eq. 30) go to 1000
        j = j + 1
        !
        !     ********** FORM SHIFT **********
        !
        l1   = l + 1
        g    = d(l)
        p    = (d(l1) - g) / (two * e(l))
        rr   = sqrt(p*p+one)
        d(l) = e(l) / (p + sign(rr,p))
        h    = g - d(l)
        !
        do 140 i = l1, n
        140 d(i) = d(i) - h
        !
        f = f + h
        !
        !     ********** QL TRANSFORMATION **********
        !
        p   = d(m)
        c   = one
        s   = zero
        mml = m - l
        !
        !     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
        !
        do 200 ii = 1, mml
            i = m - ii
            g = c * e(i)
            h = c * p
            if (abs(p) .lt. abs(e(i))) go to 150
            c      = e(i) / p
            rr     = sqrt(c*c+one)
            e(i+1) = s * p * rr
            s      = c / rr
            c      = one / rr
            go to 160
            150 c = p / e(i)
            rr     = sqrt(c*c+one)
            e(i+1) = s * e(i) * rr
            s      = one / rr
            c      = c * s
            160 p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
            !
            !     ********** FORM VECTOR **********
            !
            do 180 k = 1, n
                h        = z(k,i+1)
                z(k,i+1) = s * z(k,i) + c * h
                z(k,i)   = c * z(k,i) - s * h
            180 continue
            !
        200 continue
        !
        e(l) = s * p
        d(l) = c * p
        if (abs(e(l)) .gt. b) go to 130
        220 d(l) = d(l) + f
    240 continue
    !
    !     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
    !
    do 300 ii = 2, n
        i = ii - 1
        k = i
        p = d(i)
        !
        do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
        260 continue
        !
        if (k .eq. i) go to 300
        d(k) = d(i)
        d(i) = p
        !
        do 280 j = 1, n
            p      = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
        280 continue
        !
    300 continue
    !
    go to 1001
    !
    !     ********** SET ERROR -- NO CONVERGENCE TO AN
    !                EIGENVALUE AFTER 30 ITERATIONS **********
    !
    1000 ierr = l
    1001 return
    !
    !     ********** LAST CARD OF TQL2 **********
    !
end subroutine tql2
!
!> @brief Evaluate TRBAK3 for the BrH2_LEPS_C_1982 potential.
subroutine trbak3( nm, n, nv, a, m, z )
    !
    implicit none

    double precision :: a, anuzero, cart, cartnu, deesdr, degsdr, deijdr, descart, dgscart, dijcart, easyab, easyac, &
        & easybc, engyes, engygs, engyij, ezero, h, pengyes, pengygs, pengyij, r, s, z, zero
    integer :: i, icartr, ik, indexes, irctnt, iz, j, k, l, m, mder, msurf, n, nasurf, natoms, nder, nflag, nm, nulbl, &
        & nv
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
    dimension a(nv),z(nm,m)
    !
    !     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
    !     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
    !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
    !
    !     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
    !     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
    !     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.
    !
    !     ON INPUT-
    !
    !        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
    !          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
    !          DIMENSION STATEMENT,
    !
    !        N IS THE ORDER OF THE MATRIX,
    !
    !        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
    !          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
    !
    !        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
    !          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST
    !          N*(N+1)/2 POSITIONS,
    !
    !        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
    !
    !        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
    !          IN ITS FIRST M COLUMNS.
    !
    !     ON OUTPUT-
    !
    !        Z CONTAINS THE TRANSFORMED EIGENVECTORS
    !          IN ITS FIRST M COLUMNS.
    !
    !     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
    !
    !     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
    !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
    !
    !     ------------------------------------------------------------------
    !
    data zero/0.0d0/
    if (m .eq. 0) go to 200
    if (n .eq. 1) go to 200
    !
    do 140 i = 2, n
        l  = i - 1
        iz = (i * l) / 2
        ik = iz + i
        h  = a(ik)
        if (h .eq. zero) go to 140
        !
        do 130 j = 1, m
            s  = zero
            ik = iz
            !
            do 110 k = 1, l
                ik = ik + 1
                s  = s + a(ik) * z(k,j)
            110 continue
            !
            !     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********
            !
            s  = (s / h) / h
            ik = iz
            !
            do 120 k = 1, l
                ik     = ik + 1
                z(k,j) = z(k,j) - s * a(ik)
            120 continue
            !
        130 continue
        !
    140 continue
    !
    200 return
    !
    !     ********** LAST CARD OF TRBAK3 **********
    !
end subroutine trbak3

!> @brief Evaluate POTINFO for the BrH2_LEPS_C_1982 potential.
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


!> @brief Evaluate ANCVRT for the BrH2_LEPS_C_1982 potential.
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
!> @brief Evaluate CARTOU for the BrH2_LEPS_C_1982 potential.
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

!> @brief Evaluate CARTTOR for the BrH2_LEPS_C_1982 potential.
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


!> @brief Evaluate EUNITZERO for the BrH2_LEPS_C_1982 potential.
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

!> @brief Evaluate RTOCART for the BrH2_LEPS_C_1982 potential.
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

!> @brief Evaluate DEDCOU for the BrH2_LEPS_C_1982 potential.
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
