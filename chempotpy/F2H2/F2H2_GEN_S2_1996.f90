!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the F2H2_GEN_S2_1996 potential.
subroutine pes( x, igrad, p, g, d )

    implicit none

    double precision :: anuzero, cart, descart, dgscart, dijcart, pengyes, pengygs, pengyij, v
    integer          :: iatom, idir, istate, nasurf, nder, nflag, nulbl
    ! number of electronic state
    integer, parameter            :: nstates = 1
    integer, parameter            :: natoms  = 4
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

    if(first_time_data) then
        call prepot
        first_time_data = .false.
    endif

    call pot

    if (igrad==0) then
        do istate=1,nstates
            p(istate) = pengygs/8065.7112013
        enddo
    else
        write (*,*) 'Only energy is available'
    endif

end subroutine pes


!
!> @brief Evaluate prepot for the F2H2_GEN_S2_1996 potential.
subroutine prepot
    !
    !   System:          F2H2
    !   Common name:     S2
    !   Number of derivatives: 0
    !   Number of bodies: 4
    !   Number of electronic surfaces: 1
    !
    !   Interface:       potlib2001
    !
    !   Note:            This surface is based on the SQSBDE surface
    !                    by M.A. Suhm.  It was refit to give better
    !                    tunneling splitting for the ground state.
    !
    !   Protocol:
    !
    !      PREPOT - initializes the potential's variables
    !               must be called once before any calls to POT
    !      POT    - driver for the evaluation of the energy
    !
    !   Units:
    !
    !      energies    - cm-1
    !      coordinates - distances: bohr; angles: degrees
    !
    !   Surfaces:
    !
    !      ground electronic state
    !
    ! M. A. SUHM  SEPTEMBER 1990
    !
    ! INPUT : SIX INTERNAL COORDINATES
    ! --R1,R2,R IN ATOMIC UNITS (BOHR)
    ! --TH1,TH2,TAU IN RADIANS
    ! DEFINITION AS IN QUACK+SUHM MOL.PHYS. 69 (1990) P791
    ! OUTPUT : HF DIMER (SQSBDE) POTENTIAL ENERGY IN HARTREE
    !
    ! THE PARAMETERS OF THE POTENTIAL ARE READ IN FROM UNIT 1
    ! (FILENAME GA.IN)
    !
    ! THE 6-D POTENTIAL IS BASED ON 1070 AB INITIO POINTS FROM
    ! KOFRANEK ET AL (CHEM PHYS 121(1988) 137) UP TO 42000 CM-1
    ! AND A PERTURBATION TREATMENT BY RIJKS AND WORMER
    ! (JCP 90(1989) 6507) AS WELL AS MONOMER We,WeXe FROM
    ! HUBER/HERZBERG. 29 FREE, 15 MONOMER OR AB INITIO
    ! BASED AND 18 CONSTRAINED PARAMETERS ARE USED.
    ! SQSBDE IS AN EMPIRICALLY REFINED VERSION WHICH HAS THE
    ! CORRECT EXPERIMENTAL HYDROGEN BOND WAVENUMBER OF
    ! 1062CM-1 (MILLER, ACC. CHEM. RES. 23 (1990) P10) AND
    ! B ROTATIONAL CONSTANT. FOR THIS PURPOSE, THE AB INITIO
    ! RAB COORDINATES WERE SCALED BY 1/1.035.
    ! THE RMS OF THE WEIGHTED FIT (W=1 FOR R1,2=/=1.7374 OR
    ! R1,2=1.7374 AND E<2220CM-1, W=0.1 FOR R1,2=1.7374 AND
    ! 2220<E<3210CM-1, W=0.001 FOR R1,2=1.7374 AND E>3210CM-1)
    ! IS 32.8CM-1.
    !
    ! E(1.7445,1.7404,5.1440,9.0007,64.1361,180)=-1.307CM-1 IS
    ! THE ABSOLUTE MINIMUM OF THE POTENTIAL
    !
    ! M. A. SUHM  SEPTEMBER 1989
    !
    ! preparatory calculations for gssqs-routine
    !
    implicit none

    double precision :: a, ac, anuzero, ar, cart, cartnu, cg3j, da, deesdr, degsdr, deijdr, descart, dgscart, dijcart, e, &
        & efct, engyes, engygs, engyij, esum, ezero, f, fac, fct, g, pengyes, pengygs, pengyij, pf, r, sum, x
    integer :: i, i1, i2, ia, ib, ic, icartr, id, ie, if, indexes, irctnt, j, j1, j2, j3, l, m, m1, m2, m3, mder, msurf, &
        & nasurf, natoms, nder, nflag, nmax, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    integer, parameter          :: ma    = 135
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
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
    common/potcm/  a(ma),da(ma),fac(0:2*4), &
        & fct(401),efct(401), &
        & ar(0:4,0:4,0:2*4), &
        & ac(3,0:4,0:4,0:2*4)

    common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)
    !
    do i=1,5
        ref(i) = ' '
    end do
    !
    !   References:      W.C. Necoechea, D.G. Truhlar Chem. Phys. Lett. 248,
    !                    (1995) 182
    !                    M. Quack, M. Suhm Mol. Phys. 69 (1990) P791
    !                    Kofranek et al Chem. Phys. 121 (1988) 137
    !                    Rijks and Wormer J. Chem. Phys. 90(1989) 6507
    !
    ref(1) = 'W.C. Necoechea, D.G. Truhlar,'
    ref(2) = 'Chem. Phys. Lett. 248, 182(1995)'
    ref(3) = 'M. Quack, M. Suhm, Mol. Phys. 69, P791(1990)'
    ref(4) = 'Kofranek et al., Chem. Phys. 121, 137(1988)'
    ref(5) = 'Rijks and Wormer, J. Chem. Phys. 90, 6507(1989)'
    !
    indexes(1) = 1
    indexes(2) = 9
    indexes(3) = 1
    indexes(4) = 9
    !
    irctnt = 3
    !
    !      CALL POTINFO
    !
    call ancvrt
    !
    ! prefactors for polynomials
    !
    fac(0) = 1.
    do 101 i=1,2*4
        fac(i) = fac(i-1)*i
    101 continue
    do 102 l=0,4
        do 103 m=0,l
            pf(l,m) = (-1)**m*sqrt((2*l+1)*fac(l-m)/(4*pi*fac(l+m)))
        103 continue
    102 continue
    !
    ! 3j-symbols for potential energy surface
    ! (abc/d-d0)=g(a,b,c,d)
    !
    !
    fct(1)  = 1.0d0
    efct(1) = 0.0d0
    do 104 i=1,400
        x       = i
        j       = i+1
        efct(j) = efct(i)
        fct(j)  = fct(i)*x
        1 if(fct(j).ge.10.0d0) then
            fct(j)  = 0.1d0*fct(j)
            efct(j) = efct(j)+1.0d0
            goto 1
        endif
    104 continue
    !
    !
    m2 = 0
    m3 = 0
    do 108 j1=0,4*2,2
        do 107 j2=0,4*2,2
            do 106 j3=iabs(j1-j2),iabs(j1+j2),4
                do 105 m1=0,min(j1,j2),2
                    m2   = -m1
                    i1   = max0(0,(j2-j3-m1)/2,(j1-j3+m2)/2)
                    i2   = min0((j1+j2-j3)/2,(j1-m1)/2,(j2+m2)/2)
                    ia   = i1
                    ib   = (j1+j2-j3)/2-i1
                    ic   = (j1-m1)/2-i1
                    id   = (j2+m2)/2-i1
                    ie   = (j3-j2+m1)/2+i1
                    if   = (j3-j1-m2)/2+i1
                    sum  = (-1)**i1/(fct(ia+1)*fct(ib+1)*fct(ic+1)*fct(id+1)* fct(ie+1)*fct(if+1))
                    esum = -efct(ia+1)-efct(ib+1)-efct(ic+1)-efct(id+1)- efct(ie+1)-efct(if+1)
                    nmax = i2-i1
                    f    = 1.0d0
                    if(nmax.gt.0) then
                        do 109 i=nmax,1,-1
                            x = (ib-i+1)*(ic-i+1)*(id-i+1)
                            x = x/((ia+i)*(ie+i)*(if+i))
                            f = 1.0d0-x*f
                        109 continue
                    endif
                    sum = sum*f
                    f=fct((j1+j2-j3)/2+1)*fct((j1-j2+j3)/2+1)* fct((j2-j1+j3)/2+1)*fct((j1+m1)/2+1)* fct((j1-m1)/2+1)* &
                        & fct((j2+m2)/2+1)*fct((j2-m2)/2+1)*fct((j3+m3)/2+1)* fct((j3-m3)/2+1)/fct((j1+j2+j3)/2+2)
                    e=efct((j1+j2-j3)/2+1)+efct((j1-j2+j3)/2+1)+ efct((j2-j1+j3)/2+1)+efct((j1+m1)/2+1)+ efct((j1-m1)/ &
                        & 2+ 1)+ efct((j2+m2)/2+1)+efct((j2-m2)/2+1)+ efct((j3+m3)/2+1)+ efct((j3-m3)/2+1)-efct((j1+ &
                        & j2+j3)/ 2+2)
                    cg3j                   = sum*dsqrt(f)*10.0d0**(esum+0.5d0*e)* (-1)**((j2-j1-m3)/2)
                    g(j1/2,j2/2,j3/2,m1/2) = cg3j
                105 continue
            106 continue
        107 continue
    108 continue
    !
    !
end subroutine prepot
!
! M. A. SUHM  SEPTEMBER 1989
!
!> @brief Evaluate gsrsqs for the F2H2_GEN_S2_1996 potential.
subroutine gsrsqs
    implicit none

    double precision :: a, ac, anuzero, ar, cart, cartnu, da, dab1, dab2, dab3, deesdr, degsdr, deijdr, descart, dgscart, &
        & dijcart, efct, engyes, engygs, engyij, ezero, fac, fct, g, pengyes, pengygs, pengyij, pf, r
    integer :: icartr, indexes, irctnt, j1, j2, j3, mder, msurf, nasurf, natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    integer, parameter :: ndmax  = 1100
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    integer, parameter          :: ma    = 135
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
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
    common/potcm/  a(ma),da(ma),fac(0:2*4), &
        & fct(401),efct(401), &
        & ar(0:4,0:4,0:2*4), &
        & ac(3,0:4,0:4,0:2*4)
    common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)
    !
    !      dab1=rcm1-3.0d0
    !      dab2=rcm1-4.0d0
    !      dab3=rcm1-4.9d0
    !
    dab1 = r(1)-3.0d0
    dab2 = r(1)-4.0d0
    dab3 = r(1)-4.09d0
    do 301 j1=0,4,1
        do 302 j2=j1,4,1
            do 303 j3=(j2-j1),(j1+j2),2
                ar(j1,j2,j3)=ac(1,j1,j2,j3)*exp(-0.20d0*dab1**2)+ ac(2,j1,j2,j3)*exp(-0.35d0*dab2**2)+ ac(3,j1,j2,j3)* &
                    & exp(-0.80d0*dab3**2)
                ar(j2,j1,j3) = (-1)**(j1+j2)*ar(j1,j2,j3)
            303 continue
        302 continue
    301 continue
    !
    !
end subroutine gsrsqs
!
! M. A. SUHM  SEPTEMBER 1989
!
!> @brief Evaluate pot for the F2H2_GEN_S2_1996 potential.
subroutine pot
    !
    implicit none

    double precision :: a, ac, amf, amh, anuzero, ar, c, c2, c3, c4, c5, cart, cartnu, cs, da, deesdr, degsdr, deijdr, &
        & descart, dgscart, dijcart, e, efct, engyes, engygs, engyij, ezero, fac, fac2, fct, g, p1, p2, pengyes, pengygs, &
        & pengyij, pf, phicon, pif, ptntl, r, r1, r12, r14, r2, r23, r34, rcm, s, s2, s3, s4, s5, tau, th1, th2, vp, x, &
        & y, z
    integer :: icartr, icu, indexes, irctnt, j, j1, j2, j3, m, mder, msurf, nasurf, natoms, nder, nflag, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    integer, parameter          :: ma    = 135
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
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
    common/potcm/  a(ma),da(ma),fac(0:2*4), &
        & fct(401),efct(401), &
        & ar(0:4,0:4,0:2*4), &
        & ac(3,0:4,0:4,0:2*4)
    common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)
    !
    dimension x(4),y(4),z(4)
    dimension p1(0:4,0:4),p2(0:4,0:4),cs(4)
    dimension vp(0:4,0:4,0:2*4)
    !
    data amh/1837.152756d0/amf/34631.9678d0/
    !
    !     PUT COORDINATES IN PROPER ARRAYS
    !
    call cartou
    call carttor
    !
    !*********************************
    !                                *
    !   Assign entries from R-array  *
    !   to variables used internally *
    !   by the POT subroutine        *
    !                                *
    !*********************************
    !
    rcm = r(1)
    r1  = r(2)
    r2  = r(3)
    th1 = r(4)
    th2 = r(5)
    tau = r(6)
    !         rcm=2.5d0
    !         r1=0.9d0
    !         r2=0.9d0
    !         th1=1.04719755120d0
    !         th2=1.04719755120d0
    !         tau=0.78539816340d0
    phicon = 0.01745329252d0
    pif    = 1.0d0/(2.0d0*sqrt(pi))
    fac2   = amh/(amh+amf)
    x(3)   = 0.0d0-r1*cos(th1)*fac2
    y(3)   = 0.0d0-r1*sin(th1)*fac2
    z(3)   = 0.0d0
    x(1)   = r1*cos(th1)*(1.0d0-fac2)
    y(1)   = r1*sin(th1)*(1.0d0-fac2)
    z(1)   = 0.0d0
    x(4)   = rcm-r2*cos(th2)*fac2
    y(4)   = 0.0d0-r2*sin(th2)*cos(tau)*fac2
    z(4)   = 0.0d0-r2*sin(th2)*sin(tau)*fac2
    x(2)   = rcm+r2*cos(th2)*(1.0d0-fac2)
    y(2)   = 0.0d0+r2*sin(th2)*cos(tau)*(1.0d0-fac2)
    z(2)   = 0.0d0+r2*sin(th2)*sin(tau)*(1.0d0-fac2)
    r12    = sqrt((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
    r23    = sqrt((x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2)
    r14    = sqrt((x(4)-x(1))**2+(y(4)-y(1))**2+(z(4)-z(1))**2)
    r34    = sqrt((x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2)
    do 101 j=1,4
        cs(j) = cos(j*tau)
    101 continue
    s       = sin(th1)
    c       = cos(th1)
    s2      = s*s
    s3      = s2*s
    s4      = s3*s
    s5      = s4*s
    c2      = c*c
    c3      = c2*c
    c4      = c3*c
    c5      = c4*c
    p1(0,0) = pf(0,0)
    p1(1,0) = pf(1,0)*c
    p1(2,0) = pf(2,0)*0.5d0*(3.0d0*c2-1.0d0)
    p1(3,0) = pf(3,0)*0.5d0*(5.0d0*c3-3.0d0*c)
    p1(4,0) = pf(4,0)*0.125d0*(35.0d0*c4-30.0d0*c2+3.0d0)
    p1(1,1) = pf(1,1)*s
    p1(2,1) = pf(2,1)*3.0d0*s*c
    p1(3,1) = pf(3,1)*1.5d0*s*(5.0d0*c2-1.0d0)
    p1(4,1) = pf(4,1)*2.5d0*s*(7.0d0*c3-3.0d0*c)
    p1(2,2) = pf(2,2)*3.0d0*s2
    p1(3,2) = pf(3,2)*15.0d0*s2*c
    p1(4,2) = pf(4,2)*7.5d0*s2*(7.0d0*c2-1.0d0)
    p1(3,3) = pf(3,3)*15.0d0*s3
    p1(4,3) = pf(4,3)*105.0d0*s3*c
    p1(4,4) = pf(4,4)*105.0d0*s4
    s       = sin(th2)
    c       = cos(th2)
    s2      = s*s
    s3      = s2*s
    s4      = s3*s
    s5      = s4*s
    c2      = c*c
    c3      = c2*c
    c4      = c3*c
    c5      = c4*c
    p2(0,0) = pf(0,0)
    p2(1,0) = pf(1,0)*c
    p2(2,0) = pf(2,0)*0.5d0*(3.0d0*c2-1.0d0)
    p2(3,0) = pf(3,0)*0.5d0*(5.0d0*c3-3.0d0*c)
    p2(4,0) = pf(4,0)*0.125d0*(35.0d0*c4-30.0d0*c2+3.0d0)
    p2(1,1) = pf(1,1)*s
    p2(2,1) = pf(2,1)*3.0d0*s*c
    p2(3,1) = pf(3,1)*1.5d0*s*(5.0d0*c2-1.0d0)
    p2(4,1) = pf(4,1)*2.5d0*s*(7.0d0*c3-3.0d0*c)
    p2(2,2) = pf(2,2)*3.0d0*s2
    p2(3,2) = pf(3,2)*15.0d0*s2*c
    p2(4,2) = pf(4,2)*7.5d0*s2*(7.0d0*c2-1.0d0)
    p2(3,3) = pf(3,3)*15.0d0*s3
    p2(4,3) = pf(4,3)*105.0d0*s3*c
    p2(4,4) = pf(4,4)*105.0d0*s4
    icu     = 31
    do 102 j1=0,4,1
        do 103 j2=j1,4,1
            do 104 j3=(j2-j1),(j1+j2),2
                ac(1,j1,j2,j3) = a(icu)
                ac(2,j1,j2,j3) = a(icu+35)
                ac(3,j1,j2,j3) = a(icu+70)
                icu            = icu+1
            104 continue
        103 continue
    102 continue
    !
    call gsrsqs
    !
    do 105 j1=0,4,1
        do 106 j2=0,4,1
            do 107 j3=abs(j2-j1),(j2+j1),2
                vp(j1,j2,j3) = g(j1,j2,j3,0)*p1(j1,0)*p2(j2,0)
                do 108 m=1,min(j1,j2),1
                    vp(j1,j2,j3) = vp(j1,j2,j3)+(-1)**m*2*g(j1,j2,j3,m)* cs(m)*p1(j1,m)*p2(j2,m)
                108 continue
                vp(j1,j2,j3) = vp(j1,j2,j3)*(-1)**(j1-j2)*(2*j3+1)*pif
            107 continue
        106 continue
    105 continue
    !
    e = 0.0d0
    !
    do 109 j1=0,4,1
        do 110 j2=0,4,1
            do 111 j3=abs(j2-j1),(j2+j1),2
                e = e+vp(j1,j2,j3)*ar(j1,j2,j3)
            111 continue
        110 continue
    109 continue
    !
    !     e=e+a(1)
    !
    if(rcm.lt.7.2d0)then
        e = e+vp(0,0,0)*exp(-((7.2d0-rcm)/rcm)**2)*a(2)/rcm**6
    else
        e = e+vp(0,0,0)*a(2)/rcm**6
    endif
    if(rcm.lt.8.4d0)then
        e = e+vp(0,0,0)*exp(-((8.4d0-rcm)/rcm)**2)*a(3)/rcm**8
    else
        e = e+vp(0,0,0)*a(3)/rcm**8
    endif
    if(rcm.lt.9.6d0)then
        e = e+vp(0,0,0)*exp(-((9.6d0-rcm)/rcm)**2)*a(4)/rcm**10
    else
        e = e+vp(0,0,0)*a(4)/rcm**10
    endif
    if(rcm.lt.9.5d0) then
        e = e+vp(1,1,2)*a(5)/rcm**3*exp(-((9.5d0-rcm)/rcm)**2)
    else
        e = e+vp(1,1,2)*a(5)/rcm**3
    endif
    if(rcm.lt.10.5d0) then
        e = e+(vp(2,1,3)-vp(1,2,3))*a(6)/rcm**4* exp(-((10.5d0-rcm)/rcm)**2)
    else
        e = e+(vp(2,1,3)-vp(1,2,3))*a(6)/rcm**4
    endif
    if(rcm.lt.11.5d0) then
        e = e+vp(2,2,4)*a(7)/rcm**5*exp(-((11.5d0-rcm)/rcm)**2)
        e = e+(vp(3,1,4)+vp(1,3,4))*a(8)/rcm**5* exp(-((11.5d0-rcm)/rcm)**2)
    else
        e = e+vp(2,2,4)*a(7)/rcm**5
        e = e+(vp(3,1,4)+vp(1,3,4))*a(8)/rcm**5
    endif
    e = e+a(9)*exp(-0.5d0*r12)
    e = e+a(10)*exp(-1.5d0*r12)
    e = e+a(13)*exp(-1.5d0*r34)
    e = e+a(17)*(exp(-10.0d0*r14)+exp(-10.0d0*r23))
    if(rcm.lt.4.5d0)then
        e = e+a(18)/rcm**6*((3.0d0*(cos(th1))**2+1.0d0)+ (3.0d0*(cos(th2))**2+1.0d0))* exp(-((4.5d0-rcm)/rcm)**2)
    else
        e = e+a(18)/rcm**6*((3.0d0*(cos(th1))**2+1.0d0)+ (3.0d0*(cos(th2))**2+1.0d0))
    endif
    !
    ! monomer terms
    !
    e=e+47635.0d0*(1.0d0-a(19)*exp(-0.1d0*(r23-3.0d0)**2))* (1.0d0-exp(-1.1953d0* (1.0d0-a(26)*(1.0d0-tanh(a(27)*(r14- &
        & 3.0d0))) )* (r1-1.7374d0*(1.0d0 -a(20)* (1.0d0-tanh(a(21)*(r14-3.0d0))) ))))**2
    e=e+47635.0d0*(1.0d0-a(19)*exp(-0.1d0*(r14-3.0d0)**2))* (1.0d0-exp(-1.1953d0* (1.0d0-a(26)*(1.0d0-tanh(a(27)*(r23- &
        & 3.0d0))) )* (r2-1.7374d0*(1.0d0 -a(20)* (1.0d0-tanh(a(21)*(r23-3.0d0))) ))))**2
    !
    !     next line corrects for asymptotic diatomic potential which is
    !     added in separately.
    !
    e = e-(ptntl(r1-1.7374d0)+ptntl(r2-1.7374d0))*219474.63067d0
    e = e+a(28)*(exp(-4.0d0*(r14-r1)**2)+exp(-4.0d0*(r23-r2)**2))
    e = e+a(22)*exp(-4.0d0*r14)*exp(-4.0d0*r23)
    e = e+a(23)
    e = e+a(24)*(exp(-r14)+exp(-r23))* (exp(-15.0d0*r1)+exp(-15.0d0*r2))
    !
    !   Units conversion:  Conversion to atomic units
    !
    e = e/219474.63067d0
    !
    engygs = e
    !
    !
    call eunitzero
    !
    if(nder.ne.0) then
        call rtocart
        call dedcou
    endif
    !
end subroutine pot
!
!> @brief Evaluate PTNTL for the F2H2_GEN_S2_1996 potential.
function ptntl( x )
    !
    !     HF DIATOMIC POTENTIAL CURVE FOR SQSBDE SURFACE, CALCULATED USING
    !     PARAMETERS GIVEN IN JCP, VOL. 95(5), P. 31.
    !     NOTICE THAT THE PARAMETER DM IS GIVEN IN ATOMIC UNITS HERE AND
    !     IN CM-1 IN THE JCP REFERENCE.
    !
    implicit none

    double precision :: a, ac, anuzero, ar, cart, cartnu, da, deesdr, degsdr, deijdr, descart, dgscart, dijcart, efct, &
        & engyes, engygs, engyij, ezero, fac, fct, g, pengyes, pengygs, pengyij, pf, ptntl, r, x
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
    integer, parameter          :: ma    = 135
    double precision, parameter :: am    = 1.1953d0
    double precision, parameter :: dm    = 0.21704d0
    double precision, parameter :: one   = 1.0d0
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
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
    common/potcm/  a(ma),da(ma),fac(0:2*4), &
        & fct(401),efct(401), &
        & ar(0:4,0:4,0:2*4), &
        & ac(3,0:4,0:4,0:2*4)

    common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)
    !
    !      SAVE
    !
    ptntl = dm*(one-exp(-am*x))**2
    !
end function ptntl
!
block data ptpacm
    !
    implicit none

    double precision :: a, ac, anuzero, ar, cart, cartnu, da, efct, ezero, fac, fct, g, pf
    integer          :: icartr, indexes, irctnt, mder, msurf, nasurf, &
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
    integer, parameter          :: ma    = 135
    !
    common/pt3cm/ ezero(isurf+1)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common/potcm/  a(ma),da(ma),fac(0:2*4), &
        & fct(401),efct(401), &
        & ar(0:4,0:4,0:2*4), &
        & ac(3,0:4,0:4,0:2*4)
    common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)
    !
    data nasurf /1,35*0/
    data nder /0/
    data nflag /1,1,15*0,6,0,0/
    !
    data anuzero /0.0d0/
    data icartr,msurf,mder/4,0,0/
    data nulbl /25*0/
    data natoms /4/
    !
    data a(1),da(1)  / 1558.0d0,  0.0d0/
    data a(2),da(2)  / -202900000.0d0,  0.0d0/
    data a(3),da(3)  / -3703000000.0d0,  0.0d0/
    data a(4),da(4)  / -11350000000.0d0,  0.0d0/
    data a(5),da(5)  / -1736693.094999999d0,  0.0d0/
    data a(6),da(6)  / -4416487.722000003d0,  0.0d0/
    data a(7),da(7)  / 16012389.50999999d0,  0.0d0/
    data a(8),da(8)  / -6389444.854999989d0,  0.0d0/
    data a(9),da(9)  / -866.86737555133d0,  117.5897899275601d0/
    data a(10),da(10)  / 129789.3634701688d0,  1937.563168749279d0/
    data a(11),da(11)  / 0.0d0,  0.0d0/
    data a(12),da(12)  / 0.0d0,  0.0d0/
    data a(13),da(13)  / 2622184.756658763d0,  29726.63014604442d0/
    data a(14),da(14)  / 0.0d0,  0.0d0/
    data a(15),da(15)  / 0.0d0,  0.0d0/
    data a(16),da(16)  / 0.0d0,  0.0d0/
    data a(17),da(17)  / 184254701843.3105d0,  17553092934.44629d0/
    data a(18),da(18)  / -320000.0d0,  0.0d0/
    data a(19),da(19)  / 3.326438943860243d-2,  1.353645657747024d-3/
    data a(20),da(20)  / -4.576887206887037d-3,  1.865017706495818d-4/
    data a(21),da(21)  / 0.4097306660396534d0,  2.777020547964315d-2/
    data a(22),da(22)  / 32623551943283.75d0,  3630087056915.047d0/
    data a(23),da(23)  / 0.0d0,  0.0d0/
    data a(24),da(24)  / 71754497894202.0d0,  2525041881368.016d0/
    data a(25),da(25)  / 0.0d0,  0.0d0/
    data a(26),da(26)  / 5.032387872479127d-2,  1.860431838521345d-3/
    data a(27),da(27)  / 1.46832057915281d0,  9.71324803167759d-2/
    data a(28),da(28)  / 3384.631559286528d0,  85.55678229458908d0/
    data a(29),da(29)  / 0.0d0,  0.0d0/
    data a(30),da(30)  / 0.0d0,  0.0d0/
    data a(31),da(31)  / -160267.101272529d0,  4055.338776717705d0/
    data a(32),da(32)  / 0.0d0,  0.0d0/
    data a(33),da(33)  / 0.0d0,  0.0d0/
    data a(34),da(34)  / 0.0d0,  0.0d0/
    data a(35),da(35)  / 0.0d0,  0.0d0/
    data a(36),da(36)  / -16268.89465981041d0,  541.1396773442175d0/
    data a(37),da(37)  / 0.0d0,  0.0d0/
    data a(38),da(38)  / 0.0d0,  0.0d0/
    data a(39),da(39)  / 2802.576101972474d0,  56.80382831682482d0/
    data a(40),da(40)  / 0.0d0,  0.0d0/
    data a(41),da(41)  / 0.0d0,  0.0d0/
    data a(42),da(42)  / 0.0d0,  0.0d0/
    data a(43),da(43)  / 0.0d0,  0.0d0/
    data a(44),da(44)  / 0.0d0,  0.0d0/
    data a(45),da(45)  / 0.0d0,  0.0d0/
    data a(46),da(46)  / 2860.063943395988d0,  44.36333848084973d0/
    data a(47),da(47)  / 0.0d0,  0.0d0/
    data a(48),da(48)  / 0.0d0,  0.0d0/
    data a(49),da(49)  / 0.0d0,  0.0d0/
    data a(50),da(50)  / 0.0d0,  0.0d0/
    data a(51),da(51)  / 0.0d0,  0.0d0/
    data a(52),da(52)  / 509.1020221188392d0,  18.73758800462167d0/
    data a(53),da(53)  / 0.0d0,  0.0d0/
    data a(54),da(54)  / 0.0d0,  0.0d0/
    data a(55),da(55)  / 0.0d0,  0.0d0/
    data a(56),da(56)  / 0.0d0,  0.0d0/
    data a(57),da(57)  / 0.0d0,  0.0d0/
    data a(58),da(58)  / 0.0d0,  0.0d0/
    data a(59),da(59)  / 0.0d0,  0.0d0/
    data a(60),da(60)  / 0.0d0,  0.0d0/
    data a(61),da(61)  / 0.0d0,  0.0d0/
    data a(62),da(62)  / 0.0d0,  0.0d0/
    data a(63),da(63)  / 0.0d0,  0.0d0/
    data a(64),da(64)  / 0.0d0,  0.0d0/
    data a(65),da(65)  / 0.0d0,  0.0d0/
    data a(66),da(66)  / 0.0d0,  0.0d0/
    data a(67),da(67)  / 0.0d0,  0.0d0/
    data a(68),da(68)  / -3467.234170178781d0,  52.91056956195189d0/
    data a(69),da(69)  / 0.0d0,  0.0d0/
    data a(70),da(70)  / 0.0d0,  0.0d0/
    data a(71),da(71)  / 0.0d0,  0.0d0/
    data a(72),da(72)  / 0.0d0,  0.0d0/
    data a(73),da(73)  / 2304.632858038298d0,  44.02327215328341d0/
    data a(74),da(74)  / 0.0d0,  0.0d0/
    data a(75),da(75)  / 0.0d0,  0.0d0/
    data a(76),da(76)  / 0.0d0,  0.0d0/
    data a(77),da(77)  / 0.0d0,  0.0d0/
    data a(78),da(78)  / 0.0d0,  0.0d0/
    data a(79),da(79)  / 0.0d0,  0.0d0/
    data a(80),da(80)  / 0.0d0,  0.0d0/
    data a(81),da(81)  / 0.0d0,  0.0d0/
    data a(82),da(82)  / 0.0d0,  0.0d0/
    data a(83),da(83)  / 0.0d0,  0.0d0/
    data a(84),da(84)  / -644.8150830210325d0,  18.1151965381963d0/
    data a(85),da(85)  / 0.0d0,  0.0d0/
    data a(86),da(86)  / 0.0d0,  0.0d0/
    data a(87),da(87)  / 0.0d0,  0.0d0/
    data a(88),da(88)  / 0.0d0,  0.0d0/
    data a(89),da(89)  / 0.0d0,  0.0d0/
    data a(90),da(90)  / 0.0d0,  0.0d0/
    data a(91),da(91)  / 0.0d0,  0.0d0/
    data a(92),da(92)  / 0.0d0,  0.0d0/
    data a(93),da(93)  / 0.0d0,  0.0d0/
    data a(94),da(94)  / 0.0d0,  0.0d0/
    data a(95),da(95)  / 0.0d0,  0.0d0/
    data a(96),da(96)  / 0.0d0,  0.0d0/
    data a(97),da(97)  / 0.0d0,  0.0d0/
    data a(98),da(98)  / 0.0d0,  0.0d0/
    data a(99),da(99)  / 0.0d0,  0.0d0/
    data a(100),da(100)  / 0.0d0,  0.0d0/
    data a(101),da(101)  / 21489.15561258642d0,  989.768118754484d0/
    data a(102),da(102)  / 1895.153684423094d0,  111.7819362295832d0/
    data a(103),da(103)  / 0.0d0,  0.0d0/
    data a(104),da(104)  / 768.3954452456346d0,  18.43253058647917d0/
    data a(105),da(105)  / -395.2301759253078d0,  12.84348126946628d0/
    data a(106),da(106)  / 2041.09421622329d0,  246.1094667134039d0/
    data a(107),da(107)  / -2422.952146150754d0,  70.83017811025411d0/
    data a(108),da(108)  / 0.0d0,  0.0d0/
    data a(109),da(109)  / 0.0d0,  0.0d0/
    data a(110),da(110)  / -297.3729718893483d0,  19.35802742071735d0/
    data a(111),da(111)  / 0.0d0,  0.0d0/
    data a(112),da(112)  / 0.0d0,  0.0d0/
    data a(113),da(113)  / 39.46815694764132d0,  11.34523118392423d0/
    data a(114),da(114)  / 0.0d0,  0.0d0/
    data a(115),da(115)  / 720.5552230122485d0,  33.69084291038598d0/
    data a(116),da(116)  / 0.0d0,  0.0d0/
    data a(117),da(117)  / 0.0d0,  0.0d0/
    data a(118),da(118)  / 0.0d0,  0.0d0/
    data a(119),da(119)  / 0.0d0,  0.0d0/
    data a(120),da(120)  / 0.0d0,  0.0d0/
    data a(121),da(121)  / 0.0d0,  0.0d0/
    data a(122),da(122)  / 0.0d0,  0.0d0/
    data a(123),da(123)  / 0.0d0,  0.0d0/
    data a(124),da(124)  / 0.0d0,  0.0d0/
    data a(125),da(125)  / 0.0d0,  0.0d0/
    data a(126),da(126)  / 0.0d0,  0.0d0/
    data a(127),da(127)  / 0.0d0,  0.0d0/
    data a(128),da(128)  / 0.0d0,  0.0d0/
    data a(129),da(129)  / 0.0d0,  0.0d0/
    data a(130),da(130)  / 0.0d0,  0.0d0/
    data a(131),da(131)  / 0.0d0,  0.0d0/
    data a(132),da(132)  / 0.0d0,  0.0d0/
    data a(133),da(133)  / 0.0d0,  0.0d0/
    data a(134),da(134)  / 0.0d0,  0.0d0/
    data a(135),da(135)  / 0.0d0,  0.0d0/
    !
end
!> @brief Evaluate POTINFO for the F2H2_GEN_S2_1996 potential.
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


!> @brief Evaluate ANCVRT for the F2H2_GEN_S2_1996 potential.
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
!> @brief Evaluate CARTOU for the F2H2_GEN_S2_1996 potential.
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

!> @brief Evaluate CARTTOR for the F2H2_GEN_S2_1996 potential.
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


!> @brief Evaluate EUNITZERO for the F2H2_GEN_S2_1996 potential.
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

!> @brief Evaluate RTOCART for the F2H2_GEN_S2_1996 potential.
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

!> @brief Evaluate DEDCOU for the F2H2_GEN_S2_1996 potential.
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
