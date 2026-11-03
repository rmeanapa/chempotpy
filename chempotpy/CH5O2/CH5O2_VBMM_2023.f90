!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the CH5O2_VBMM_2023 potential.
subroutine pes( x, igrad, p, g, d )

    implicit none

    ! number of electronic state
    integer, parameter :: nstates = 1
    integer, parameter :: natoms  = 8
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
            cart(iatom,idir) = x(iatom,idir)
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
    else
        write (*,*) 'Only energy is available'
    endif

end subroutine pes

!:**********************************************************************
!     SUBROUTINE SURF(V, COORD, DX, N3TM)
!> @brief Evaluate POT for the CH5O2_VBMM_2023 potential.
subroutine pot
    !
    !   System:    CH3OH + OH --> CH2OH + H2O.
    !                        --> CH3O + H2O
    !
    ! Como una primera prueba SOLO considero el canal 1, que es el
    ! mayoritario
    !              Joaquin Espinosa-Garcia Junio 2022
    ! Esta version inlcuye el potencial de van der Waals
    !
    !
    !   All the information passed to and from the potential energy surface
    !   routine is in hartree atomic units.
    !
    !        This potential is written such that:
    !                       X(1)  - X(3)  : X, Y, Z for C
    !                       X(4)  - X(6)  : X, Y, Z for O
    !                       X(7)  - X(9)  : X, Y, Z for H1
    !                       X(10) - X(12) : X, Y, Z for H2
    !                       X(13) - X(15) : X, Y, Z for H3
    !                       X(16) - X(18) : X, Y, Z for H4
    !                       X(19) - X(21) : X, Y, Z for O
    !                       X(22) - X(24) : X, Y, Z for H(O)
    !
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, bphi1, bphi2, &
        & btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, d1c1h, d1c2b, d1c2h, &
        & d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, dedho, dedo, deesdr, &
        & degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, en, engyes, engygs, engyij, evdw, ezero, f1, fch3, &
        & fdelta, fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0c1b, r0c1h, &
        & r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, rb2h, rc1b, rc1h, &
        & rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t2, tb1h, tb2h, tc1b, tc1h, tc2b, tc2h, tcc, theta0, &
        & thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, vip, vop, vstr, vtor, w1, w10, w2, w3, w4, w5, w6, w7, w8, &
        & w9
    integer :: i, icartr, indexes, irctnt, mder, msurf, n, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, &
        & nflag, nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !
    !     DIMENSION COORD(24),DX(24)
    !     DIMENSION dxvstr(21),dxvtor(21),dxvop(21),dxvip(21)
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
    !  include file for CH3OH+OH PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    ! Start of changes by EMN.
    ! The lines below do not seem to be needed because the Cartesian
    ! coordinates are in common qpdot.
    ! Besides, BOHR_A in utility.f does not match the value below.
    ! As a consequence, q(i) changes and energy is not conserved
    !c call de venus
    call cartou
    call carttor
    !
    !     PUT COORDINATES IN PROPER ARRAYS
    !
    ! Inicialization de variables
    !
    do i=1,27
        q(i)    = r(i)
        pdot(i) = 0.d0
    enddo
    !
    ! Changing units from angstrom to bohr
    !
    !      DO  I = 1, 27
    !c         q(I) = COORD(I)*0.52918d0
    !        q(I) = R(I)*0.52918d0
    !      enddo
    !
    ! End of corrections by EMN
    !
    !  calculate relative coordinates and bond lengths
    !
    en = 0.0d0
    call coorden
    !
    call torsion(vtor)
    !
    !  calculate switching functions
    !
    call switchf
    !
    !  calculate reference angles and their derivatives
    !
    call refangles
    !
    !  calculate stretching potential
    !
    call stretch(vstr)
    !
    !  calculate out of plane bending potential
    !
    call opbend(vop)
    !
    !  jeg:calculate van der Waals term, vdw
    !
    call vander(evdw)
    !
    !  calculate in plane bending potential
    !
    call ipbend(vip)
    !
    !  total potential energy is vstr+vop+vip
    !
    en = vstr+vop+vip+vtor+evdw
    !
    !  convert from 10(5) j/mol to au
    !
    en = en*0.03812d0
    !      V = en
    engygs = en
    !
    call eunitzero
    if(nder.ne.0) then
        !
        !  Initialize derivatives
        !
        do n=1,24
            !          DX(n)=0.0d0
            degsdr(i) = pdot(i)*0.0201723d0
        enddo
        call rtocart
        call dedcou
    endif
    !
end subroutine pot
!******************************************************
!
!------------------------------------------------------
!> @brief Evaluate coorden for the CH5O2_VBMM_2023 potential.
subroutine coorden
    !------------------------------------------------------
    !
    !  calculates relative coordinates and bond lengths
    !
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, argmax, argp5, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, &
        & bphi1, bphi2, btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, d1c1h, &
        & d1c2b, d1c2h, d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, dedho, &
        & dedo, deesdr, degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, engyes, engygs, engyij, ezero, f1, fch3, &
        & fdelta, fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, p5, pdot, pengyes, pengygs, pengyij, q, r, r0c1b, &
        & r0c1h, r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rav4, rb1h, rb2h, &
        & rc1b, rc1h, rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t2, t5tmp, tb1h, tb2h, tc1b, tc1h, tc2b, tc2h, &
        & tcc, theta0, thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, w4, w5, w6, w7, w8, w9
    integer :: i, icartr, ind, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, &
        & nflag, nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !      include 'ch3oh_h.inc'
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    !
    !  calculate relative coordinates
    !
    do ind=1,3
        tc1b(ind) = q(nc1(ind))-q(nhb(ind))
        tc2b(ind) = q(nc2(ind))-q(nhb(ind))
        tcc(ind)  = q(nc1(ind))-q(nc2(ind))
        tno(ind)  = q(no(ind))-q(nhb(ind))
        ! jeg vdw
        toh(ind) = q(no(ind))-q(nc2(ind))
    enddo
    !
    do ind=1,3
        do i=1,3
            tc1h(i,ind) = q(nc1(ind))-q(n1h(i,ind))
        enddo
    enddo
    !
    do ind=1,3
        tc2h(1,ind) = q(nc2(ind))-q(n2h(1,ind))
    enddo
    !
    do ind=1,3
        do i=1,3
            tb1h(i,ind) = q(nhb(ind))-q(n1h(i,ind))
        enddo
        tb2h(1,ind) = q(nhb(ind))-q(n2h(1,ind))
    enddo
    !
    rc1b = sqrt(tc1b(1)*tc1b(1)+tc1b(2)*tc1b(2)+tc1b(3)*tc1b(3))
    rc2b = sqrt(tc2b(1)*tc2b(1)+tc2b(2)*tc2b(2)+tc2b(3)*tc2b(3))
    rcc  = sqrt(tcc(1)*tcc(1)+tcc(2)*tcc(2)+tcc(3)*tcc(3))
    rno  = sqrt(tno(1)*tno(1)+tno(2)*tno(2)+tno(3)*tno(3))
    ! jeg vdw
    roh = sqrt(toh(1)*toh(1)+toh(2)*toh(2)+toh(3)*toh(3))
    !
    do i=1,3
        rc1h(i) = sqrt(tc1h(i,1)*tc1h(i,1)+tc1h(i,2)*tc1h(i,2)+ tc1h(i,3)*tc1h(i,3))
    enddo
    !
    rc2h(1) = sqrt(tc2h(1,1)*tc2h(1,1)+tc2h(1,2)*tc2h(1,2)+ tc2h(1,3)*tc2h(1,3))
    !
    do i=1,3
        rb1h(i) = sqrt(tb1h(i,1)*tb1h(i,1)+tb1h(i,2)*tb1h(i,2)+ tb1h(i,3)*tb1h(i,3))
    enddo
    !
    rb2h(1) = sqrt(tb2h(1,1)*tb2h(1,1)+tb2h(1,2)*tb2h(1,2)+ tb2h(1,3)*tb2h(1,3))
    !
    ! crd 2013 modificaciones de las distancias
    !
    !
    ! CRD 2021 modificacion para hacer solo un parametro que
    ! modifique todos los r0c en funcion de rav
    !
    argmax = 19.d0
    rav4   = (rc1h(1)+rc1h(2)+rc1h(3)+rc2h(1))/4.0d0
    p5     = 1.d0
    argp5  = (w9*(rav4-w10))
    if (argp5.lt.argmax) then
        t5tmp = 1-tanh(argp5)
    else
        t5tmp = 0.d0
    endif
    p5    = p5*t5tmp
    r0cc  = p5*r0ccr+(1-p5)*r0ccp
    r0c1h = p5*r0c1hr+(1-p5)*r0c1hp
    r0c2h = p5*r0c2hr+(1-p5)*r0c2hp
    !c
end subroutine coorden
!
!------------------------------------------------------
! JEG 2022 subrutina para calcular el potencial de van der Waals
!> @brief Evaluate vander for the CH5O2_VBMM_2023 potential.
subroutine vander( evdw )
    !------------------------------------------------------
    !
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, bphi1, bphi2, &
        & btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, d1c1h, d1c2b, d1c2h, &
        & d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, dedho, dedo, deesdr, &
        & degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, engyes, engygs, engyij, evdw, ezero, f1, fch3, fdelta, &
        & fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, p1vdw, p2vdw, pdot, pengyes, pengygs, pengyij, q, r, r0c1b, &
        & r0c1h, r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, rb2h, rc1b, &
        & rc1h, rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t1, t2, tb1h, tb2h, tc1b, tc1h, tc2b, tc2h, tcc, &
        & theta0, thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, w4, w5, w6, w7, w8, w9
    integer :: icartr, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, &
        & natoms, nc, nc1, nc1h, nc2, nc2h, nder, nflag, nhb, nnb, nnc1, nnc2, nno, no, nulbl
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
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    dimension t1(3)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    !
    p1vdw = vdwc1*(vdwr/roh)**6.d0
    p2vdw = vdwc2*exp(-12.d0*(roh/vdwr))
    evdw  = vdwe/0.03812d0/627.52d0*(p1vdw+p2vdw)
    !
end subroutine vander
!
!
! CRD 2011 incluimos un apartado para calcular los angulos dihedro
! para poder calcular despues vtorsion
!------------------------------------------------------
! CRD 2011 subrutina para calcular el potencial de torsion
!> @brief Evaluate torsion for the CH5O2_VBMM_2023 potential.
subroutine torsion( vtor )
    !------------------------------------------------------
    !
    implicit none

    double precision :: a, a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, &
        & aa22, aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, &
        & aphi1, aphi2, arg1, arg2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, &
        & bphi, bphi1, bphi2, btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, &
        & d1c1h, d1c2b, d1c2h, d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, &
        & dedho, dedo, deesdr, degsdr, deijdr, descart, dgscart, dihed, dijcart, dkdr, dkdx, engyes, engygs, engyij, &
        & ezero, f1, fch3, fdelta, fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, &
        & r, r0c1b, r0c1h, r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, &
        & rb2h, rc1b, rc1h, rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t1, t2, tb1h, tb2h, tc1b, tc1h, tc2b, &
        & tc2h, tcc, theta0, thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, vtor, w1, w10, w2, w3, w4, w5, w6, w7, &
        & w8, w9
    integer :: i, icartr, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, nflag, &
        & nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    dimension t1(3)
    !      double precision pi
    !      pi=4.0d0*atan(1.0d0)
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    !
    ! CRD 2011 incluimos un apartado para calcular los angulos dihedro
    ! para poder calcular despues vtorsion
    !
    call dihedro(3,1,2,6,dihed)
    thetahcch(1) = dihed*pi/180.0d0
    call dihedro(4,1,2,6,dihed)
    thetahcch(2) = dihed*pi/180.0d0
    call dihedro(5,1,2,6,dihed)
    thetahcch(3) = dihed*pi/180.0d0
    !
    ! CRD 2011 calculamos los terminos switching t1 y t2
    !
    do i=1,3
        !       t1(i)=0.5d0*(1-tanh(w1*(rc1h(i)-w2)))
        !       t1(i)=(1-tanh(w1*(rc1h(i)-w2)))
        t1(i) = (1-tanh(w1*(rc1h(i)-r0c1h)))
    enddo
    !       t2(1)=0.5d0*(1-tanh(w7*(rc2h(1)-w8)))
    !       t2(1)=(1-tanh(w7*(rc2h(1)-w8)))
    t2 = (1-tanh(w7*(rc2h(1)-r0c2h)))
    !
    !     write(*,*) t1(1)*t2(1),t1(2)*t2(1),t1(3)*t2(1)
    !
    vtor = 0.0d0
    !
    ! CRD 2021 incluimos terminos para pasar de 3 a 2 minimos
    !
    !     arg1=(t1(1)*t1(2)*t1(3)*t2)**3.0d0
    arg1 = (rb1h(1)-r0h1h)*(rb1h(2)-r0h1h)*(rb1h(3)-r0h1h)* (rb2h(1)-r0h2h)
    arg2 = tanh(arg1)
    if (arg2.lt.0.0) then
        a = 2.0d0
    else
        a = 2.0d0+1.d0*arg2
    endif
    !     write(45,*) arg1, arg2, A
    ! CRD linea para eliminar la correccion
    a = 3.0d0
    !
    !      vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(1)))*t1(1)*t2
    !      vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(2)))*t1(2)*t2
    !      vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(3)))*t1(3)*t2
    vtor = vtor+(v3/3.0d0)*(1+cos(a*thetahcch(1)))*t1(1)*t2
    vtor = vtor+(v3/3.0d0)*(1+cos(a*thetahcch(2)))*t1(2)*t2
    vtor = vtor+(v3/3.0d0)*(1+cos(a*thetahcch(3)))*t1(3)*t2
end subroutine torsion
!
!------------------------------------------------------
!      double precision function dihed(i,j,k,l)
!> @brief Evaluate dihedro for the CH5O2_VBMM_2023 potential.
subroutine dihedro( i, j, k, l, dihed )
    !------------------------------------------------------
    !
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, bphi1, bphi2, &
        & btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, d1c1h, d1c2b, d1c2h, &
        & d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, dedho, dedo, deesdr, &
        & degsdr, deijdr, descart, dgscart, dihed, dijcart, dkdr, dkdx, engyes, engygs, engyij, ezero, f1, fch3, fdelta, &
        & fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0c1b, r0c1h, r0c1hp, &
        & r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, rb2h, rc1b, rc1h, rc2b, &
        & rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t2, tb1h, tb2h, tc1b, tc1h, tc2b, tc2h, tcc, theta0, thetahcch, &
        & tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, w4, w5, w6, w7, w8, w9
    integer :: icartr, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, &
        & natoms, nc, nc1, nc1h, nc2, nc2h, nder, nflag, nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !        include 'ch3oh_h.inc'
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    integer :: i, j, k,l
    integer :: m
    !        double precision x(300)
    !        double precision dista
    !        double precision d12,d13,d14,d23,d24,d34
    !        double precision p,q
    !        double precision p
    !        double precision pi
    !
    double precision :: anorm
    double precision :: dotprod
    double precision :: rij(3),rjk(3),rkj(3),rkl(3)
    double precision :: zijjk(3),zijkj(3), zkjkl(3), zjkkl(3)
    double precision :: tmpsign(3), tmp, arg1,arg2
    !
    !        pi=4.0d0*atan(1.0d0)
    !
    do m = 1, 3
        rij(m) = q((i-1)*3+m)-q((j-1)*3+m)
        rjk(m) = q((j-1)*3+m)-q((k-1)*3+m)
        rkl(m) = q((k-1)*3+m)-q((l-1)*3+m)
    enddo
    call crosprod(rjk,rkl,zjkkl)
    call crosprod(rij,rjk,zijjk)
    arg2 = anorm(rjk)*dotprod(rij,zjkkl)
    arg1 = dotprod(zijjk,zjkkl)
    !
    !
    if (abs(arg1).lt.1.d-6) then
        if (arg2.gt.0.d0) dihed = 90.d0
        if (arg2.lt.0.d0) dihed = -90.d0
        return
    endif
    if (arg1.gt.0.d0) then
        dihed = atan(arg2/arg1) *180.d0/pi
        return
    endif
    if (arg2.ge.0.d0) then
        dihed = 180.d0 + atan(arg2/arg1) *180.d0/pi
        return
    endif
    if (arg2.lt.0.d0) then
        dihed = -180.d0 + atan(arg2/arg1) *180.d0/pi
        return
    endif
    !
    dihed = dihed*pi/180.d0
    !
end subroutine dihedro
!
!------------------------------------------------------
!> @brief Evaluate crosprod for the CH5O2_VBMM_2023 potential.
subroutine crosprod( x, y, z )
    !------------------------------------------------------
    !
    implicit none

    double precision :: x(3),y(3),z(3)
    !
    z(1) = x(2)*y(3)-x(3)*y(2)
    z(2) = x(3)*y(1)-x(1)*y(3)
    z(3) = x(1)*y(2)-x(2)*y(1)
    !
end subroutine crosprod
!
!------------------------------------------------------
!> @brief Evaluate dotprod for the CH5O2_VBMM_2023 potential.
double precision function dotprod( x, y )
    !------------------------------------------------------
    !
    implicit none

    integer          :: k
    double precision :: x(3),y(3), tmp
    !
    tmp = 0.d0
    do k = 1, 3
        tmp = tmp + x(k)*y(k)
    enddo
    dotprod = tmp
    !
end function dotprod
!
!------------------------------------------------------
!> @brief Evaluate anorm for the CH5O2_VBMM_2023 potential.
double precision function anorm( x )
    !------------------------------------------------------
    implicit none

    integer          :: k
    double precision :: x(3), tmp
    !
    tmp = 0.d0
    do k = 1, 3
        tmp = tmp + x(k)*x(k)
    enddo
    anorm = dsqrt(tmp)
    !
end function anorm
!
!------------------------------------------------------
!> @brief Evaluate refangles for the CH5O2_VBMM_2023 potential.
subroutine refangles
    !------------------------------------------------------
    !
    !  subroutine calculates reference angles for the "in-plane" potential
    !
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, bphi1, bphi2, &
        & btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, d1c1h, d1c2b, d1c2h, &
        & d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, dedho, dedo, deesdr, &
        & degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, &
        & fkh2o, fkh2oeq, fkinf1, fkinf2, halfpi, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0c1b, r0c1h, &
        & r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, rb2h, rc1b, rc1h, &
        & rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t2, tau, tauco, tauoc, tb1h, tb2h, tc1b, tc1h, tc2b, tc2h, &
        & tcc, theta0, thetahcch, tno, toh, twopi, v3, vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, w4, w5, w6, w7, w8, w9
    integer :: i, icartr, indexes, irctnt, j, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, &
        & nflag, nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+OH PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    !      dimension sumd2(4),sumd4(4),ddr(4,4)
    tau = acos(-1.0d0/3.0d0)
    !
    ! CRD 2012 taucc para el angulo H-C-C
    !
    tauco = acos(-0.3616d0)
    tauoc = acos(-0.4558d0)
    !
    !      pi=4.0d0*atan(1.0d0)
    !      halfpi=0.5d0*pi
    halfpi = 0.6555d0*pi
    twopi  = 2.0d0*pi
    ! jcc-2010
    ! valor para ajustar el angulo
    !      anghch=0.649d0*pi
    !      tausih = pi - asin ( sin(anghch/2.d0) / sin(pi/3.d0) )
    !      tausih=acos(-0.4648d0)
    ! jcc-2010
    !
    !
    !  set diagonal elements to zero
    !
    do i=1,6
        theta0(i,i) = 0.0d0
    enddo
    !
    !  calculate reference angles
    !
    ! el atomo 4 es el c2
    !
    theta0(1,2) = tau+(tau-halfpi)*(sphi(1)*sphi(2)-1.0d0) +(tau-twopi/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
    theta0(1,3) = tau+(tau-halfpi)*(sphi(1)*sphi(3)-1.0d0) +(tau-twopi/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
    theta0(1,4) = tauco+(tau-halfpi)*(sphi(1)*sphi(4)-1.0d0) +(tau-twopi/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
    theta0(2,3) = tau+(tau-halfpi)*(sphi(2)*sphi(3)-1.0d0) +(tau-twopi/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
    theta0(2,4) = tauco+(tau-halfpi)*(sphi(2)*sphi(4)-1.0d0) +(tau-twopi/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
    theta0(3,4) = tauco+(tau-halfpi)*(sphi(3)*sphi(4)-1.0d0) +(tau-twopi/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
    !
    !  el atomo 6 es el c1
    !
    theta0(5,6) = tauoc+(tauoc-halfpi)*(sphi(5)*sphi(6)-1.0d0)
    !
    !  fill in the other half of the matrix
    !
    do i=1,3
        do j=i+1,4
            theta0(j,i) = theta0(i,j)
        enddo
    enddo
    theta0(6,5) = theta0(5,6)
    !
    !      do i=1,3
    !         do j=i+1,4
    !          write(44,*) 'theta',i,j,'=', theta0(i,j)*180/pi
    !         enddo
    !      enddo
    !          write(44,*) 'theta(5,6)=', theta0(5,6)*180/pi
end subroutine refangles
!
!******************************************************
!
!
!------------------------------------------------------
!> @brief Evaluate stretch for the CH5O2_VBMM_2023 potential.
subroutine stretch( vstr )
    !------------------------------------------------------
    !
    !  subroutine to calculate leps-type stretching potential and its
    !  derivatives
    !
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ach1, ach2, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, &
        & aphi, aphi1, aphi2, arga1, arga2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, &
        & b3scc, bphi, bphi1, bphi2, btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, &
        & d1c1b, d1c1h, d1c2b, d1c2h, d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, &
        & dedhi, dedho, dedo, deesdr, degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, dt, e1, e3, engyes, engygs, &
        & engyij, expterm, ezero, f1, fch3, fdelta, fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, &
        & pengygs, pengyij, q, r, r0c1b, r0c1h, r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, &
        & r0h2h, r0hh, rav1, rav2, rb1h, rb2h, rc1b, rc1h, rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t2, tb1h, &
        & tb2h, tc1b, tc1h, tc2b, tc2h, tcc, theta0, thetahcch, tno, toh, v3, vcc, vdwc1, vdwc2, vdwe, vdwr, vj, vjb1h, &
        & vjb2h, vjc1b, vjc1h, vjc2b, vjc2h, vno, vq, vqb1h, vqb2h, vqc1b, vqc1h, vqc2b, vqc2h, vstr, w1, w10, w2, w3, &
        & w4, w5, w6, w7, w8, w9
    integer :: i, icartr, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, nflag, &
        & nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !      include 'ch3oh_h.inc'
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    dimension vqb1h(3),vjb1h(3),vqb2h(1),vjb2h(1),vq(4),vj(4), &
        & vqc1h(3), vqc2h(1), vjc1h(3), vjc2h(1)
    !
    !  calculate avergage bond length for the methane moiety
    !
    rav1 = (rc1h(1)+rc1h(2)+rc1h(3))/3.0d0
    rav2 = rc2h(1)
    !      write(44,*) rav1, rav2
    !
    !  initialise:
    !
    vstr = 0.0d0
    !
    !  ach:
    !
    !  nb: in double precision tanh(19.0d0)=1.0d0 and we put the if statemen
    !  in to avoid overflow/underflow errors
    !
    arga1 = c1c1h*(rav1-r0c1h)
    arga2 = c1c2h*(rav2-r0c2h)
    if(arga1.lt.19.0d0)then
        ach1 = a1c1h+b1c1h*(tanh(arga1)+1.0d0)*0.5d0
    else
        ach1 = a1c1h+b1c1h
    endif
    if(arga2.lt.19.0d0)then
        ach2 = a1c2h+b1c2h*(tanh(arga2)+1.0d0)*0.5d0
    else
        ach2 = a1c2h+b1c2h
    endif
    !      write(44,*) 'r0c1h',r0c1h,'rav1',rav1
    !      write(44,*) 'r0c2h',r0c2h, 'rav2',rav2
    !      write(44,*) 'ach1',ach1,'ach2', ach2
    !
    !  calculate singlet: e1, triplet: e3 energies and vq and vj
    !  terms for each bond
    !
    !      write(44,*) d1c1b,d3c2b,ac1b,rc1b,r0c1b
    !      write(44,*) d1c2b,d3c2b,ac2b,rc2b,r0c2b
    e1    = d1c1b*(exp(-2.0d0*ac1b*(rc1b-r0c1b)) -2.0d0*exp(-ac1b*(rc1b-r0c1b)))
    e3    = d3c1b*(exp(-2.0d0*ac1b*(rc1b-r0c1b)) +2.0d0*exp(-ac1b*(rc1b-r0c1b)))
    vqc1b = (e1+e3)*0.5d0
    vjc1b = (e1-e3)*0.5d0
    !      write(44,*) vqc1b, vjc1b
    !
    e1    = d1c2b*(exp(-2.0d0*ac2b*(rc2b-r0c2b)) -2.0d0*exp(-ac2b*(rc2b-r0c2b)))
    e3    = d3c2b*(exp(-2.0d0*ac2b*(rc2b-r0c2b)) +2.0d0*exp(-ac2b*(rc2b-r0c2b)))
    vqc2b = (e1+e3)*0.5d0
    vjc2b = (e1-e3)*0.5d0
    !     write(44,*) vqc2b, vjc2b
    ! O-H  new term (simple morse term)
    !
    dt      = (rno-r0hh)
    expterm = exp(-ahh*dt)
    vno     = d1hh*(1.d0-expterm)**2.d0
    !
    do i=1,3
        !      write(44,*) 'ach1',ach1,'rc1h',rc1h(i),r0c1h
        e1       = d1c1h*(exp(-2.0d0*ach1*(rc1h(i)-r0c1h)) -2.0d0*exp(-ach1*(rc1h(i)-r0c1h)))
        e3       = d3c1h*(exp(-2.0d0*ach1*(rc1h(i)-r0c1h)) +2.0d0*exp(-ach1*(rc1h(i)-r0c1h)))
        vqc1h(i) = (e1+e3)*0.5d0
        vjc1h(i) = (e1-e3)*0.5d0
        !      write(44,*) 'vqc1h',i, vqc1h(i), 'vjc1h',i, vjc1h(i)
    enddo
    !
    e1       = d1c2h*(exp(-2.0d0*ach2*(rc2h(1)-r0c2h)) -2.0d0*exp(-ach2*(rc2h(1)-r0c2h)))
    e3       = d3c2h*(exp(-2.0d0*ach2*(rc2h(1)-r0c2h)) +2.0d0*exp(-ach2*(rc2h(1)-r0c2h)))
    vqc2h(1) = (e1+e3)*0.5d0
    vjc2h(1) = (e1-e3)*0.5d0
    !      write(44,*) 'vqc2h', vqc2h(1), 'vjc2h',vjc2h(1)
    !
    do i=1,3
        e1       = d1h1h*(exp(-2.0d0*ah1h*(rb1h(i)-r0h1h)) -2.0d0*exp(-ah1h*(rb1h(i)-r0h1h)))
        e3       = d3h1h*(exp(-2.0d0*ah1h*(rb1h(i)-r0h1h)) +2.0d0*exp(-ah1h*(rb1h(i)-r0h1h)))
        vqb1h(i) = (e1+e3)*0.5d0
        vjb1h(i) = (e1-e3)*0.5d0
        !      write(44,*) 'vqb1h',i, vqb1h(i), 'vjb1h',i, vjb1h(i)
    enddo
    e1       = d1h2h*(exp(-2.0d0*ah2h*(rb2h(1)-r0h2h)) -2.0d0*exp(-ah2h*(rb2h(1)-r0h2h)))
    e3       = d3h2h*(exp(-2.0d0*ah2h*(rb2h(1)-r0h2h)) +2.0d0*exp(-ah2h*(rb2h(1)-r0h2h)))
    vqb2h(1) = (e1+e3)*0.5d0
    vjb2h(1) = (e1-e3)*0.5d0
    !      write(44,*) 'vqb2h', vqb2h(1), 'vji2h',vjb2h(1)
    !
    !  calculate 3 body potential
    !
    do i=1,3
        !       write(44,*) i,'vqc1h',vqc1h(i),'vqc1b',vqc1b, 'vqb1h',vqb1h(i)
        vq(i) = vqc1h(i)+vqc1b+vqb1h(i)
        vj(i) = -sqrt(((vjc1h(i)-vjc1b)**2+(vjc1b-vjb1h(i))**2 +(vjb1h(i)-vjc1h(i))**2)*0.5d0)
        vstr  = vstr+vq(i)+vj(i)
    enddo
    !
    !
    vq(4) = vqc2h(1)+vqc2b+vqb2h(1)
    vj(4) = -sqrt(((vjc2h(1)-vjc2b)**2+(vjc2b-vjb2h(1))**2 +(vjb2h(1)-vjc2h(1))**2)*0.5d0)
    vstr  = vstr+vq(4)+vj(4)
    !
    ! CRD 2012 potencial morse para el enlace C-C
    !
    dt      = (rcc-r0cc)
    expterm = exp(-acc*dt)
    vcc     = d1cc*(1.0d0-expterm)**2.0d0
    vstr    = vstr + vcc + vno
    !      write(44,*) 'rcc',rcc,'r0cc',r0cc,'acc',acc
    !
end subroutine stretch
!
!
!------------------------------------------------------
!> @brief Evaluate opbend for the CH5O2_VBMM_2023 potential.
subroutine opbend( vop )
    !------------------------------------------------------
    !
    !  subroutine calculates symmetrized vop potential and derivatives
    !
    implicit none

    double precision :: a, a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, &
        & aa22, aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, &
        & aphi1, aphi2, argd, atheta, athetacc, axb, b, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, &
        & bphi, bphi1, bphi2, btheta, bthetacc, c, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, &
        & d1c1b, d1c1h, d1c2b, d1c2h, d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, &
        & dedhi, dedho, dedo, deesdr, degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, engyes, engygs, engyij, &
        & ezero, f1, fch3, fdelta, fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, &
        & r, r0c1b, r0c1h, r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, &
        & rb2h, rc1b, rc1h, rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, sum2, sum4, sumd2, sumd4, t2, tb1h, &
        & tb2h, tc1b, tc1h, tc2b, tc2h, tcc, theta0, thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, vop, w1, w10, w2, &
        & w3, w4, w5, w6, w7, w8, w9
    integer :: i, icartr, in, indexes, irctnt, j, k, l, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, &
        & nder, nflag, nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    double precision :: norma
    dimension sumd2(4),sumd4(4)
    dimension in(3),a(3),b(3),axb(3),c(4,3),argd(4)
    !      dimension sumd2(8),sumd4(8)
    !      dimension in(3),a(3),b(3),axb(3),c(8,3),argd(8)
    !
    !
    vop = 0.0d0
    !
    !  calculate force constants
    !
    call opforce
    !
    !  calculate out-of-plane angle
    !
    ! CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C2
    !
    do i=1,4
        j = i+1
        if(j.gt.4)j=j-4
        k = j+1
        if(k.gt.4)k=k-4
        l = k+1
        if(l.gt.4)l=l-4
        !
        call calcdelta1(i,j,k,l,sum2,sum4)
        sumd2(i) = sum2
        sumd4(i) = sum4
        vop      = vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
        !      write(44,*) i,j,k,l,sumd2(i)*fdelta(i),sumd4(i)*hdelta(i),vop
    enddo
    !
end subroutine opbend
!
!------------------------------------------------------
!> @brief Evaluate calcdelta1 for the CH5O2_VBMM_2023 potential.
subroutine calcdelta1( i, j, k, l, sum2, sum4 )
    !------------------------------------------------------
    !
    !  subroutine calculates out of plane angle delta, loops
    !  through delta(i,j), delta(i,k), delta(i,l)
    !
    !   also calculates the derivatives wrt delta
    !
    implicit none

    double precision :: a, a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, &
        & aa22, aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, &
        & aphi1, aphi2, argd, atemp2, atheta, athetacc, axb, b, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, &
        & b3scc, bphi, bphi1, bphi2, btheta, bthetacc, c, c1c1h, c1c2h, cart, cartnu, cdot, cphi, cphi1, cphi2, ctheta, &
        & cthetacc, d1c1b, d1c1h, d1c2b, d1c2h, d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, &
        & daxb, ded, dedhi, dedho, dedo, deesdr, degsdr, deijdr, delta, descart, dgscart, dijcart, dkdr, dkdx, engyes, &
        & engygs, engyij, ezero, f1, fch3, fdelta, fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, &
        & pengygs, pengyij, q, r, r0c1b, r0c1h, r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, &
        & r0h2h, r0hh, rb1h, rb2h, rc1b, rc1h, rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, sum2, sum4, t2, tb1h, &
        & tb2h, tc1b, tc1h, tc2b, tc2h, tcc, theta0, thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, &
        & w4, w5, w6, w7, w8, w9
    integer :: i, icartr, ii, in, ind, indexes, irctnt, j, k, l, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, &
        & nc2, nc2h, nder, nflag, nhb, nnb, nnc1, nnc2, nno, no, nulbl
    double precision :: norma
    !      include 'ch3oh_h.inc'
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !

    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    dimension  delta(8),in(3),a(3),b(3),axb(3),c(8,3),argd(8), &
        & daxb(8,3,3),cdot(8,3,3),atemp2(3)
    !
    !  initialise
    !
    sum2 = 0.0d0
    sum4 = 0.0d0
    !
    !  set j,k,l indices
    !
    in(1) = j
    in(2) = k
    in(3) = l
    !
    ! CRD 2005 incluimos bucles if para asegurar que siempre
    ! que una de los valores j, k o l sea 4 las coordenadas
    ! corresponden al atomo nuevo.
    !
    !  vector a is rk-rj, vector b is rl-rj
    !
    if(j.eq.4)then
        do ind=1,3
            a(ind) = q(n1h(k,ind))-q(nc2(ind))
            b(ind) = q(n1h(l,ind))-q(nc2(ind))
        enddo
        axb(1) = a(2)*b(3)-a(3)*b(2)
        axb(2) = a(3)*b(1)-a(1)*b(3)
        axb(3) = a(1)*b(2)-a(2)*b(1)
        norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
        norma  = sqrt(norma)
    else
    endif
    !
    if(k.eq.4)then
        do ind=1,3
            a(ind) = q(nc2(ind))-q(n1h(j,ind))
            b(ind) = q(n1h(l,ind))-q(n1h(j,ind))
        enddo
        axb(1) = a(2)*b(3)-a(3)*b(2)
        axb(2) = a(3)*b(1)-a(1)*b(3)
        axb(3) = a(1)*b(2)-a(2)*b(1)
        norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
        norma  = sqrt(norma)
    else
    endif
    !
    if(l.eq.4)then
        do ind=1,3
            a(ind) = q(n1h(k,ind))-q(n1h(j,ind))
            b(ind) = q(nc2(ind))-q(n1h(j,ind))
        enddo
        axb(1) = a(2)*b(3)-a(3)*b(2)
        axb(2) = a(3)*b(1)-a(1)*b(3)
        axb(3) = a(1)*b(2)-a(2)*b(1)
        !
        ! modificacion para que hiciese vxu
        !      axb(1)=b(2)*a(3)-b(3)*a(2)
        !      axb(2)=a(1)*b(3)-b(1)*a(3)
        !      axb(3)=b(1)*a(2)-a(1)*b(2)
        !
        norma = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
        norma = sqrt(norma)
    else
    endif
    !
    ! CRD 2012 eliminamos porque seria la salida del C2
    !
    if(k.lt.4.and.j.lt.4.and.l.lt.4)then
        do ind=1,3
            a(ind) = q(n1h(k,ind))-q(n1h(j,ind))
            b(ind) = q(n1h(l,ind))-q(n1h(j,ind))
        enddo
        axb(1) = a(2)*b(3)-a(3)*b(2)
        axb(2) = a(3)*b(1)-a(1)*b(3)
        axb(3) = a(1)*b(2)-a(2)*b(1)
        norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
        norma  = sqrt(norma)
    else
    endif
    !
    !  c is position vector of h(ii): calculate c(j),c(k),c(l)
    !
    do ii=1,3
        if(in(ii).eq.4) then
            do ind=1,3
                c(in(ii),ind) = -tcc(ind)/rcc
            enddo
        else
            do ind=1,3
                c(in(ii),ind) = -tc1h(in(ii),ind)/rc1h(in(ii))
            enddo
        endif
    enddo
    !
    !  argd is the dot product axb dot c
    !
    do ii=1,3
        argd(in(ii))  = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
        argd(in(ii))  = argd(in(ii))/norma
        delta(in(ii)) = acos(argd(in(ii)))-theta0(i,in(ii))
        sum2          = sum2+delta(in(ii))**2
        sum4          = sum4+delta(in(ii))**4
    enddo
    !
end subroutine calcdelta1
!
!------------------------------------------------------
!> @brief Evaluate opforce for the CH5O2_VBMM_2023 potential.
subroutine opforce
    !------------------------------------------------------
    !
    !  calculates the out-of-plane bending force constants
    !  and their derivatives
    !
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, bphi1, bphi2, &
        & btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, d1c1h, d1c2b, d1c2h, &
        & d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, dedho, dedo, deesdr, &
        & degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, &
        & fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0c1b, r0c1h, r0c1hp, &
        & r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, rb2h, rc1b, rc1h, rc2b, &
        & rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, switch, t2, tb1h, tb2h, tc1b, tc1h, tc2b, tc2h, tcc, theta0, &
        & thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, w4, w5, w6, w7, w8, w9
    integer :: i, icartr, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, nflag, &
        & nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !      include 'ch3oh_h.inc'
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    dimension switch(4)
    !      dimension switch(8)
    !
    !  calculate switching functions:
    !
    switch(1) = (1.0d0-s3(1))*s3(2)*s3(3)*s3(4)
    switch(2) = (1.0d0-s3(2))*s3(3)*s3(4)*s3(1)
    switch(3) = (1.0d0-s3(3))*s3(4)*s3(1)*s3(2)
    switch(4) = (1.0d0-s3(4))*s3(1)*s3(2)*s3(3)
    !      switch(5)=(1.0d0-s3(5))*s3(6)*s3(7)*s3(8)
    !      switch(6)=(1.0d0-s3(6))*s3(7)*s3(8)*s3(5)
    !      switch(7)=(1.0d0-s3(7))*s3(8)*s3(5)*s3(6)
    !      switch(8)=(1.0d0-s3(8))*s3(5)*s3(6)*s3(7)
    !      do i=1,4
    !       write(44,*) 'switch',i,'=',switch(i)
    !      enddo
    !       write(*,*) 'pasa'
    !
    !  calculate the force constants
    !
    !      do i=1,8
    do i=1,4
        fdelta(i) = switch(i)*fch3
        hdelta(i) = switch(i)*hch3
    enddo
    !
end subroutine opforce
!
!
!------------------------------------------------------
!> @brief Evaluate ipbend for the CH5O2_VBMM_2023 potential.
subroutine ipbend( vip )
    !------------------------------------------------------
    !
    !  subroutine calculates symmetrised in plane bend term
    !
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, arga, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, bphi1, &
        & bphi2, btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cosine, costh, cphi, cphi1, cphi2, ctheta, cthetacc, &
        & d1c1b, d1c1h, d1c2b, d1c2h, d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, dang, &
        & ded, dedhi, dedho, dedo, deesdr, degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, dot, dth, engyes, &
        & engygs, engyij, ezero, f1, fch3, fdelta, fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, &
        & pengygs, pengyij, q, r, r0c1b, r0c1h, r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, &
        & r0h2h, r0hh, rb1h, rb2h, rc1b, rc1h, rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t2, tb1h, tb2h, tc1b, &
        & tc1h, tc2b, tc2h, tcc, theta, theta0, thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, vip, w1, w10, w2, w3, &
        & w4, w5, w6, w7, w8, w9
    integer :: i, icartr, indexes, irctnt, j, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, &
        & nflag, nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    dimension costh(6,6),theta(6,6),dth(6,6)
    !
    !  initialise
    !
    vip = 0.0d0
    !
    !  calculate force constants: fk0(i,j), f1(i)
    !  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
    !
    call ipforce
    !
    !  calculate theta(i,j) and in plane bend potential
    !
    ! c1
    !
    do i=1,3
        do j=i+1,4
            if(j.lt.4.0) then
                costh(i,j) = tc1h(i,1)*tc1h(j,1)+tc1h(i,2)*tc1h(j,2) +tc1h(i,3)*tc1h(j,3)
                costh(i,j) = costh(i,j)/rc1h(i)/rc1h(j)
                theta(i,j) = acos(costh(i,j))
                dth(i,j)   = theta(i,j)-theta0(i,j)
                !      write(44,*) 'f1 i',i,f1(i),'f1 j',j,f1(j)
                vip = vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
            else
                costh(i,j) = tc1h(i,1)*tcc(1)+tc1h(i,2)*tcc(2) +tc1h(i,3)*tcc(3)
                costh(i,j) = costh(i,j)/rc1h(i)/rcc
                theta(i,j) = acos(costh(i,j))
                dth(i,j)   = theta(i,j)-theta0(i,j)
                vip        = vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
                !     write(44,*)'fk0i,j',fk0(i,j),'f1i',f1(i),'f1j',f1(j),'dt',dth(i,j)
            endif
        enddo
    enddo
    !     write(44,*) 'vip do 1', vip
    !
    ! c2
    !
    costh(5,6) = tc2h(1,1)*(-tcc(1))+tc2h(1,2)*(-tcc(2)) +tc2h(1,3)*(-tcc(3))
    costh(5,6) = costh(5,6)/rc2h(1)/rcc
    theta(5,6) = acos(costh(5,6))
    dth(5,6)   = theta(5,6)-theta0(5,6)
    vip        = vip+0.5d0*fk0(5,6)*f1(5)*f1(6)*dth(5,6)**2
    !     write(44,*) 'theta 5,6', theta(5,6), theta0(5,6)
    !     write(44,*)'fk05,6',fk0(5,6),'f15',f1(5),'f16',f1(6),'dt',dth(5,6)
    !     write(44,*) 'vip do 2', vip
    !
    !  Now, calculate h2o bending energies.
    !
    !  First, calculate the angles:
    ! CRD 2022 INICIO
    do i = 1,3
        dot = 0.d0
        do j = 1,3
            dot = dot - tno(j)*tb1h(i,j)
        enddo
        cosine    = dot / (rno*rb1h(i))
        cosine    = min(1.d0, max(-1.d0, cosine))
        angh2o(i) = acos(cosine)
    enddo
    dot = 0.d0
    do j = 1,3
        dot = dot - tno(j)*tb2h(1,j)
    enddo
    cosine    = dot / (rno*rb2h(1))
    cosine    = min(1.d0, max(-1.d0, cosine))
    angh2o(4) = acos(cosine)
    ! CRD 2022 FIN
    !  Now, calculate each force constant fkh2o(i) as a function of the
    !  O-H(I)
    !  distance, rbh(i)
    ! CRD 2022 INICIO
    do i=1,3
        arga = alph2o * (rb1h(i) - r0hh)
        if(arga.lt.19.0d0)then
            fkh2o(i) = fkh2oeq * (1 - tanh(arga))
        else
            fkh2o(i) = 0.0d0
        endif
    enddo
    arga = alph2o * (rb2h(1) - r0hh)
    if(arga.lt.19.0d0)then
        fkh2o(4) = fkh2oeq * (1 - tanh(arga))
    else
        fkh2o(4) = 0.0d0
    endif
    ! CRD 2022 FIN
    !  Now calculate the contribution to the energy for each H(I)-O-H(O)
    !  harmonic bending term
    !
    do i=1,4
        dang = (angh2o(i) - angh2oeq)
        vip  = vip+0.5d0* fkh2o(i) * dang * dang
    enddo
    !
end subroutine ipbend
!
!------------------------------------------------------
!> @brief Evaluate ipforce for the CH5O2_VBMM_2023 potential.
subroutine ipforce
    !------------------------------------------------------
    !
    !  calculates the symmetrised in plane bend force constants and
    !  all partial derivatives involving them
    !
    implicit none

    double precision :: a1, a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, &
        & aa21, aa22, aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, &
        & aphi, aphi1, aphi2, arga1, arga2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, &
        & b3scc, bphi, bphi1, bphi2, btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, &
        & d1c1b, d1c1h, d1c2b, d1c2h, d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, &
        & dedhi, dedho, dedo, deesdr, degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, engyes, engygs, engyij, &
        & ezero, f0, f1, f2, fch3, fdelta, fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, pengygs, &
        & pengyij, q, r, r0c1b, r0c1h, r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, &
        & r0hh, rb1h, rb2h, rc1b, rc1h, rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t2, tb1h, tb2h, tc1b, tc1h, &
        & tc2b, tc2h, tcc, theta0, thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, w4, w5, w6, w7, &
        & w8, w9
    integer :: i, icartr, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, nflag, &
        & nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !      include 'ch3oh_h.inc'
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    !  set force constant at asymptotes
    !
    f0 = fkinf1+ak1
    f2 = fkinf1
    !
    fk0(1,2) = f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)*s2(4)-1.0d0)
    fk0(1,3) = f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)*s2(4)-1.0d0)
    fk0(1,4) = f0+f0*(s1(1)*s1(4)-1.0d0)+(f0-f2)*(s2(2)*s2(3)-1.0d0)
    fk0(2,3) = f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)*s2(4)-1.0d0)
    fk0(2,4) = f0+f0*(s1(2)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(3)-1.0d0)
    fk0(3,4) = f0+f0*(s1(3)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(2)-1.0d0)
    !
    !      fk0(5,6)=fkinf2*(s1(5)*s1(6))
    !      fk0(5,6)=fkinf2*s1(5)
    fk0(5,6) = f0+f0*(s1(5)*s1(6)-1.0d0)+(f0-f2)*(s2(5)*s2(6)-1.0d0)
    !
    ! f1 para los h del c1
    !
    do i=1,3
        !
        !  calc derivatives of fk0 wrt each of the rch(i) bonds
        !  calculate the terms f1(i)
        !
        !       write(44,*) 'rb1h',i, rb1h(i)
        arga1 = aa11*rb1h(i)*rb1h(i)
        arga2 = aa41*(rb1h(i)-r0hh)*(rb1h(i)-r0hh)
        a1    = 1.0d0-exp(-arga1)
        a2    = aa21+aa31*exp(-arga2)
        !       write(44,*) 'a1',i, a1,'a2',i,a2,'aa11',aa11
        f1(i) = a1*exp(-a2*(rc1h(i)-r0c1h)**2)
    enddo
    !
    ! f1 para el c2
    !
    ! CRD 2012 eliminamos el calculo de a1c y a2c pq es constante
    ! lo incluimos en el CONST
    !
    f1(4) = a1c*exp(-a2c*(rcc-r0cc)**2)
    f1(6) = f1(4)
    !
    ! f1 para los h del c2
    !
    !
    !  calc derivatives of fk0 wrt each of the rch(i) bonds
    !  calculate the terms f1(i)
    !
    arga1 = aa12*rb2h(1)*rb2h(1)
    arga2 = aa42*(rb2h(1)-r0hh)*(rb2h(1)-r0hh)
    a1    = 1.0d0-exp(-arga1)
    a2    = aa22+aa32*exp(-arga2)
    f1(5) = a1*exp(-a2*(rc2h(1)-r0c2h)**2)
    !
end subroutine ipforce
!
!******************************************************
!
!------------------------------------------------------
!> @brief Evaluate switchf for the CH5O2_VBMM_2023 potential.
subroutine switchf
    !------------------------------------------------------
    !
    !  calculates switching functions: s3,sphi,stheta
    !
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, argmax, args1, args2, args3, argsc1, argsc2, argsc3, argsphi, argsphic, argstheta, argsthetac, atheta, &
        & athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, bphi1, bphi2, btheta, bthetacc, &
        & c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, d1c1h, d1c2b, d1c2h, d1cc, d1h1h, &
        & d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, dedho, dedo, deesdr, degsdr, deijdr, &
        & descart, dgscart, dijcart, dkdr, dkdx, engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, fkh2o, fkh2oeq, &
        & fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0c1b, r0c1h, r0c1hp, r0c1hr, r0c2b, &
        & r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, rb2h, rc1b, rc1h, rc2b, rc2h, rcc, rno, &
        & roh, s1, s2, s3, sphi, stheta, t2, tb1h, tb2h, tc1b, tc1h, tc2b, tc2h, tcc, theta0, thetahcch, tno, toh, v3, &
        & vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, w4, w5, w6, w7, w8, w9
    integer :: i, icartr, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, nflag, &
        & nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !      include 'ch3oh_h.inc'
    !C INCLUDE DEL VENUS
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    !  nb remember that integration units are:
    !  energy in   1.0d+05 j/mol
    !  time in     1.0d-14 s
    !
    ! original     a1s=1.5313681d-7
    !      a1s=1.5313681d-8
    !      b1s=-4.6696246d0
    ! original      a2s=1.0147402d-7
    !      a2s=1.0147402d-8
    !      b2s=-12.363798d0
    !
    !  use double precision criterion:
    !
    !  tanh(19.0d0)=1.0d0
    !
    argmax = 19.0d0
    !
    !  calculate s1
    !
    do i=1,3
        args1 = a1s1*(rc1h(i)-r0c1h)*(rc1h(i)-b1s1)**8
        if(args1.lt.argmax)then
            s1(i) = 1.0d0-tanh(args1)
        else
            s1(i) = 0.0d0
        endif
    enddo
    !
    args1 = a2s1*(rc2h(1)-r0c2h)*(rc2h(1)-b2s1)**8
    if(args1.lt.argmax)then
        s1(5) = 1.0d0-tanh(args1)
    else
        s1(5) = 0.0d0
    endif
    !
    ! CRD 2012 todos los S de los carbonos son 1
    !
    argsc1 = a1scc*(rcc-r0cc)*(rcc-b1scc)**8
    if(argsc1.lt.argmax)then
        s1(4) = 1.0d0-tanh(argsc1)
        s1(6) = s1(4)
    else
        s1(4) = 0.0d0
        s1(6) = s1(4)
    endif
    !
    ! CRD correccion para dejar C-C fijo
    !
    s1(4) = 1.0d0
    s1(6) = 1.0d0
    !
    !  calculate s2 and ds2
    !
    do i=1,3
        args2 = a1s2*(rc1h(i)-r0c1h)*(rc1h(i)-b1s2)**6
        if(args2.lt.argmax)then
            s2(i) = 1.0d0-tanh(args2)
        else
            s2(i) = 0.0d0
        endif
    enddo
    !
    args2 = a2s2*(rc2h(1)-r0c2h)*(rc2h(1)-b2s2)**6
    if(args2.lt.argmax)then
        s2(5) = 1.0d0-tanh(args2)
    else
        s2(5) = 0.0d0
    endif
    !
    argsc2 = a2scc*(rcc-r0cc)*(rcc-b2scc)**6
    if(argsc2.lt.argmax)then
        s2(4) = 1.0d0-tanh(argsc2)
        s2(6) = s2(4)
    else
        s2(4) = 0.0d0
        s2(6) = s2(4)
    endif
    !
    ! CRD correccion para dejar C-C fijo
    !
    s2(4) = 1.0d0
    s2(6) = 1.0d0
    !
    !  calculate s3 and ds3
    !
    do i=1,3
        args3 = a3s*(rc1h(i)-r0c1h)*(rc1h(i)-b3s)**2
        if (args3.lt.argmax)then
            s3(i) = 1.0d0-tanh(args3)
        else
            s3(i) = 0.0d0
        endif
    enddo
    !
    argsc3 = a3scc*(rcc-r0cc)*(rcc-b3scc)**2
    if (argsc3.lt.argmax)then
        s3(4) = 1.0d0-tanh(argsc3)
    else
        s3(4) = 0.0d0
    endif
    !
    ! CRD correccion para dejar C-C fijo
    !
    s3(4) = 1.0d0
    !
    !  calculate sphi
    !
    !  condition here is on the bondlength rch(i)
    !  st argsphi is lt approx 19.0d0
    !
    do i=1,3
        if(rc1h(i).lt.3.8d0)then
            argsphi = aphi1*(rc1h(i)-r0c1h)*exp(bphi1*(rc1h(i)-cphi1)**3)
            sphi(i) = 1.0d0-tanh(argsphi)
        else
            sphi(i) = 0.0d0
        endif
    enddo
    !
    if(rc2h(1).lt.3.8d0)then
        argsphi = aphi2*(rc2h(1)-r0c2h)*exp(bphi2*(rc2h(1)-cphi2)**3)
        sphi(5) = 1.0d0-tanh(argsphi)
    else
        sphi(5) = 0.0d0
    endif
    !
    if(rcc.lt.3.8d0)then
        argsphic = aphi*(rcc-r0cc)*exp(bphi*(rcc-cphi)**3)
        sphi(4)  = 1.0d0-tanh(argsphic)
        sphi(6)  = sphi(4)
    else
        sphi(4) = 0.0d0
        sphi(6) = sphi(4)
    endif
    !
    !          sphi(4)=1.0d0
    !          sphi(6)=1.0d0
    !
    !      do i=1,6
    !          write(44,*) 'sphi',i,sphi(i)
    !      enddo
    !
    !  calculate stheta and dstheta
    !
    do i=1,3
        if(rc1h(i).lt.3.8d0)then
            argstheta = atheta*(rc1h(i)-r0c1h)*exp(btheta*(rc1h(i)-ctheta)**3)
            stheta(i) = 1.0d0-tanh(argstheta)
        else
            stheta(i) = 0.0d0
        endif
    enddo
    !
    !      do i=1,3
    !        if(rc2h(i).lt.3.8d0)then
    !      argstheta=atheta*(rc2h(i)-r0c2h)*exp(btheta*(rc2h(i)-ctheta)**3)
    !          stheta(i+4)=1.0d0-tanh(argstheta)
    !        else
    !          stheta(i+4)=0.0d0
    !        endif
    !      enddo
    !
    ! CRD 2012 todos los S de los carbonos son 1
    !
    if(rcc.lt.3.8d0)then
        argsthetac = athetacc*(rcc-r0cc)*exp(bthetacc*(rcc-cthetacc)**3)
        stheta(4)  = 1.0d0-tanh(argsthetac)
        stheta(6)  = stheta(4)
    else
        stheta(4) = 0.0d0
        stheta(6) = stheta(4)
    endif
    !
    !          stheta(4)=1.0d0
    !          stheta(6)=1.0d0
    !
end subroutine switchf
!
!
!> @brief Evaluate PREPOT for the CH5O2_VBMM_2023 potential.
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
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, bphi1, bphi2, &
        & btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, d1c1h, d1c2b, d1c2h, &
        & d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, dedho, dedo, deesdr, &
        & degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, engyes, engygs, engyij, ezero, f1, fact1, fact2, fact3, &
        & fch3, fdelta, fk0, fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0c1b, &
        & r0c1h, r0c1hp, r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, rb2h, rc1b, &
        & rc1h, rc2b, rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t2, tb1h, tb2h, tc1b, tc1h, tc2b, tc2h, tcc, theta0, &
        & thetahcch, tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, w4, w5, w6, w7, w8, w9
    integer :: i, icartr, icount, ind, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, &
        & nder, nflag, nhb, nnb, nnc1, nnc2, nno, no, nulbl
    !C INCLUDE DEL VENUS
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    integer, parameter          :: n3tmmn = 27
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2
    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    !c call de venus
    !
    !
    !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
    !
    if(natoms.gt.25) then
        write(nflag(18),'(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    end if
    !
    do i=1,5
        ref(i) = ' '
    end do
    !
    !      REF(1)='J. Espinosa-Garcia and J. C. Corchado'
    !      REF(2)='J. Chem. Phys., Vol. 112, p. 5731, 2000'
    !
    indexes(1) = 6
    indexes(2) = 8
    indexes(3) = 1
    indexes(4) = 1
    indexes(5) = 1
    indexes(6) = 1
    indexes(7) = 8
    indexes(8) = 1
    !
    !     IRCTNT=9
    irctnt = 7
    !
    !      CALL POTINFO
    !
    call ancvrt
    !
    !  calculate indexes for coordinates
    !
    do ind=1,3
        icount   = ind-3
        nc1(ind) = 3*nnc1+icount
        nc2(ind) = 3*nnc2+icount
        nhb(ind) = 3*nnb+icount
        no(ind)  = 3*nno+icount
        !
        do i=1,3
            n1h(i,ind) = 3*nc1h(i)+icount
        enddo
        n2h(1,ind) = 3*nc2h(1)+icount
    enddo
    !
    !
    !  convert to appropriate units:
    !
    !  From:
    !
    !  kcal/mol
    !  mdyn A-1        -> 1.0d+05 j/mol...
    !  mdyn A rad-1
    !
    !  To
    !
    !  energy   in 1.0d+05 j/mol
    !  distance in 1.0d-10 m
    !  angles   in radians
    !  mass     in amu
    !
    fact1 = 0.041840d0
    fact2 = 6.022045d0
    fact3 = 2.d0*3.1415926d0/360.d0
    !
    d1c1h = d1c1h*fact1
    d3c1h = d3c1h*fact1
    !
    d1c2h = d1c2h*fact1
    d3c2h = d3c2h*fact1
    !
    d1cc = d1cc*fact1
    !
    d1c1b = d1c1b*fact1
    d3c1b = d3c1b*fact1
    !
    d1c2b = d1c2b*fact1
    d3c2b = d3c2b*fact1
    !
    d1h1h = d1h1h*fact1
    d3h1h = d3h1h*fact1
    !
    d1h2h = d1h2h*fact1
    d3h2h = d3h2h*fact1
    !
    d1hh = d1hh*fact1
    d3hh = d3hh*fact1

    fkinf1   = fkinf1*fact2
    ak1      = ak1*fact2
    fkinf2   = fkinf2*fact2
    fkh2oeq  = fkh2oeq*fact2
    angh2oeq = angh2oeq*fact3
    !
    a1s1 = a1s1*1.d-8
    a2s1 = a2s1*1.d-8
    !
    a1s2 = a1s2*1.d-8
    a2s2 = a2s2*1.d-8
    !
    a1scc = a1scc*1.d-8
    a2scc = a2scc*1.d-8
    !
end subroutine prepot
!
block data ptpacm
    !
    !C INCLUDE DEL VENUS
    implicit none

    double precision :: a1c, a1c1h, a1c2h, a1s1, a1s2, a1scc, a2c, a2s1, a2s2, a2scc, a3s, a3scc, aa11, aa12, aa21, aa22, &
        & aa31, aa32, aa41, aa42, ac1b, ac2b, acc, ah1h, ah2h, ahh, ak1, alph2o, angh2o, angh2oeq, anuzero, aphi, aphi1, &
        & aphi2, atheta, athetacc, b1c1h, b1c2h, b1s1, b1s2, b1scc, b2s1, b2s2, b2scc, b3s, b3scc, bphi, bphi1, bphi2, &
        & btheta, bthetacc, c1c1h, c1c2h, cart, cartnu, cphi, cphi1, cphi2, ctheta, cthetacc, d1c1b, d1c1h, d1c2b, d1c2h, &
        & d1cc, d1h1h, d1h2h, d1hh, d3c1b, d3c1h, d3c2b, d3c2h, d3h1h, d3h2h, d3hh, ded, dedhi, dedho, dedo, deesdr, &
        & degsdr, deijdr, descart, dgscart, dijcart, dkdr, dkdx, engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, &
        & fkh2o, fkh2oeq, fkinf1, fkinf2, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0c1b, r0c1h, r0c1hp, &
        & r0c1hr, r0c2b, r0c2h, r0c2hp, r0c2hr, r0cc, r0ccp, r0ccr, r0h1h, r0h2h, r0hh, rb1h, rb2h, rc1b, rc1h, rc2b, &
        & rc2h, rcc, rno, roh, s1, s2, s3, sphi, stheta, t2, tb1h, tb2h, tc1b, tc1h, tc2b, tc2h, tcc, theta0, thetahcch, &
        & tno, toh, v3, vdwc1, vdwc2, vdwe, vdwr, w1, w10, w2, w3, w4, w5, w6, w7, w8, w9
    integer :: icartr, indexes, irctnt, mder, msurf, n1h, n2h, nasurf, natoms, nc, nc1, nc1h, nc2, nc2h, nder, nflag, &
        & nhb, nnb, nnc1, nnc2, nno, no, nulbl
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
    !  include file for CH3OH+H PES
    !
    !  CONST
    !
    !     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
    common /potcm/ nnc1,nnc2,nnb, &
        & r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h, &
        & a1c1h,b1c1h,c1c1h, &
        & r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h, &
        & a1c2h,b1c2h,c1c2h, &
        & r0ccr,r0ccp,w9,w10,d1cc,acc, &
        & a1c,a2c, &
        & r0h1h,d1h1h,d3h1h,ah1h, &
        & r0h2h,d1h2h,d3h2h,ah2h, &
        & r0c1b,d1c1b,d3c1b,ac1b, &
        & r0c2b,d1c2b,d3c2b,ac2b, &
        & aphi,bphi,cphi, &
        & aphi1,bphi1,cphi1, &
        & aphi2,bphi2,cphi2, &
        & a3s,b3s,a3scc,b3scc, &
        & atheta,btheta,ctheta, &
        & athetacc,bthetacc,cthetacc, &
        & fch3,hch3, &
        & fkinf1,ak1,aa11,aa21,aa31,aa41, &
        & fkinf2,aa12,aa22,aa32,aa42, &
        & v3,w1,w2,w7,w8, &
        & a1s1,a2s1,b1s1,b2s1, &
        & a1s2,a2s2,b1s2,b2s2, &
        & a1scc,a2scc,b1scc,b2scc, &
        & vdwe,vdwr,vdwc1,vdwc2

    !
    ! Subroutine Coorden
    !
    common /potcm2/ nc1h(3),nc2h(1)
    common /qpdot/ q(24),pdot(24)
    common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
    common /coord3/ tc1b(3),tc2b(3),tcc(3)
    common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
    common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
    common /coord6/ r0c1h,r0c2h,r0cc
    !
    ! Subroutine torsion
    !
    common /torsion1/ thetahcch(3),t2
    !
    ! Subroutine refangles
    !
    common /refangles1/ theta0(6,6),sphi(6),stheta(6)
    !
    ! Subroutine vop
    !
    common /vop6/ fdelta(4),hdelta(4)
    !
    ! Subroutine ipforce
    !
    common /ipforce1/ fk0(6,6),f1(6)
    !
    ! Subroutine Switch
    !
    common /switch1/ s1(6),s2(6),s3(4)
    !
    ! Add Newh2o
    !
    common /newh2o/nno, no(3), tno(3), rno, &
        & ded(3),angh2o(4),fkh2oeq,fkh2o(4), &
        & alph2o,angh2oeq,dkdr(4),dkdx(4,3), &
        & dedo(3),dedhi(3),dedho(3), &
        & ahh, d3hh, d1hh, r0hh
    !
    ! jeg vdw
    !
    common /newvdw/roh, toh(3)
    !
    data nasurf /1,35*0/
    data nder /1/
    data nflag /1,1,15*0,6,0,0/
    !
    data anuzero /0.0d0/
    data icartr,msurf,mder/1,0,1/
    data nulbl /25*0/
    data natoms /8/
    !
    ! CONST-2-O
    !
    data nnc1 /1/
    data nnc2 /2/
    data nnb /7/
    data nc1h /3,4,5/
    data nc2h /6/
    data nno /8/
    data r0c1hr /1.09047d0/
    data r0c1hp /1.08177d0/
    !     data w3 /1.0000d0/
    !     data w4 /1.09047d0/
    data d1c1h /106.23000d0/
    data d3c1h /37.00400d0/
    data a1c1h /1.60300d0/
    data b1c1h /0.52500d0/
    data c1c1h /0.55404d0/
    data r0c2hr /0.95897d0/
    data r0c2hp /0.96077d0/
    !     data w5 /1.000d0/
    !     data w6 /0.96897d0/
    data d1c2h /117.23000d0/
    data d3c2h /42.00400d0/
    data a1c2h /2.05300d0/
    data b1c2h /0.42500d0/
    data c1c2h /0.51404d0/
    data r0hh / 0.9719d0/
    data d1hh/ 130.458d0/
    data d3hh / 13.264d0/
    data ahh /2.0557d0/
    data r0ccr /1.42007d0/
    data r0ccp /1.37500d0/
    data w9 /1.5000d0/
    data w10 /1.08507d0/
    data d1cc /109.94900d0/
    data acc /1.934d0/
    data a1c /1.20997d0/
    data a2c /1.949d0/
    data r0h1h /0.97191d0/
    data d1h1h /130.458d0/
    data d3h1h /33.264d0/
    data ah1h /2.0557d0/
    data r0h2h /0.97191d0/
    data d1h2h /130.458d0/
    data d3h2h /33.264d0/
    data ah2h /2.0557d0/
    data r0c1b /1.69597d0/
    data d1c1b /45.009d0/
    data d3c1b /22.200d0/
    data ac1b /1.6390285d0/
    data r0c2b /1.85897d0/
    data d1c2b /32.609d0/
    data d3c2b /22.400d0/
    data ac2b /2.1170285d0/
    data aphi /0.528790d0/
    data bphi /  1.200664d0/
    data cphi /  1.420074d0/
    data aphi1 /   0.19700d0/
    data bphi1 /   1.20600d0/
    data cphi1 /   1.08200d0/
    data aphi2 /   0.19700d0/
    data bphi2 /   1.20600d0/
    data cphi2 /   0.96800d0/
    data a3s /   0.12770d0/
    data b3s /   1.08200d0/
    data a3scc /   0.297000d0/
    data b3scc /  1.420664d0/
    data atheta /   1.65781d0/
    data btheta /   0.950600d0/
    data ctheta /   1.08201d0/
    data athetacc /   1.45781d0/
    data bthetacc /   0.950600d0/
    data cthetacc /   1.42001d0/
    data fch3 /   0.19770d0/
    data hch3 /   0.09600d0/
    data fkinf1 /   0.44770d0/
    data ak1 /   0.12600d0/
    data aa11 /   1.49001d0/
    data aa21 /   1.59995d0/
    data aa31 /   2.16996d0/
    data aa41 /  11.56595d0/
    data fkinf2 /   0.45770d0/
    data aa12 /   0.15001d0/
    data aa22 /   1.30095d0/
    data aa32 /   1.49996d0/
    data aa42 /  11.46595d0/
    data fkh2oeq / 0.6700d0/
    data alph2o / 3.70800d0/
    data angh2oeq / 104.7132d0/
    data v3 /   0.03000d0/
    data w1 /   0.11000d0/
    data w2 /   1.08697d0/
    data w7 /   0.11000d0/
    data w8 /   1.08367d0/
    data a1s1 /   1.53137d0/
    data b1s1 /   1.09463d0/
    data a2s1 /   1.01474d0/
    data b2s1 /   0.96898d0/
    data a1s2 /   1.53137d0/
    data b1s2 /   1.09463d0/
    data a2s2 /   1.51474d0/
    data b2s2 /   0.96898d0/
    data a1scc /   1.53137d0/
    data b1scc /   1.42003d0/
    data a2scc /   1.51474d0/
    data b2scc /   1.42008d0/
    data vdwe /    0.500d0/
    data vdwr /    3.7100d0/
    data vdwc1/    -6.7000d0/
    data vdwc2 /  225000.0d0/
end

!> @brief Evaluate POTINFO for the CH5O2_VBMM_2023 potential.
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

!> @brief Evaluate ANCVRT for the CH5O2_VBMM_2023 potential.
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
!> @brief Evaluate CARTOU for the CH5O2_VBMM_2023 potential.
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

!> @brief Evaluate CARTTOR for the CH5O2_VBMM_2023 potential.
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

!> @brief Evaluate EUNITZERO for the CH5O2_VBMM_2023 potential.
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

!> @brief Evaluate RTOCART for the CH5O2_VBMM_2023 potential.
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

!> @brief Evaluate DEDCOU for the CH5O2_VBMM_2023 potential.
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
