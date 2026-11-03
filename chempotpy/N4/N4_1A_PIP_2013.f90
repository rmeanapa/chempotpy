!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!   System:                     N4
!   Functional form:            permutation-invariant polynomials
!   Common name:                N4(adiabatic ground state)
!   Number of derivatives:      1
!   Number of bodies:           4
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References:: Yuliya, Paukku, Ke R. Yang, Zoltan Varga,
!   and Donald G. Truhlar J. Chem. Phys. 139, 044309 (2013).
!
!   Notes:    PES of N4 with special emphasize for
!             N2 + N2 --> N2 + N + N
!             (Refit with the corrected data points on 12/05/2013)
!             (named n4pes-gpip-version2)
!
!     N1--N2
!
!     N3--N4
!
!
!> @brief Evaluate pes for the N4_1A_PIP_2013 potential.
subroutine pes( x, igrad, p, g, d )

    implicit none

    ! number of electronic state
    integer, parameter            :: nstates = 1
    integer, parameter            :: natoms  = 4
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    double precision, parameter   :: eref = -218.40801323d0

    double precision :: v, tx(12), tg(12)
    integer          :: iatom, idir, j, istate
    !initialize
    v = 0.d0
    g = 0.d0
    d = 0.d0

    do iatom=1,natoms
        do idir=1,3
            j     = 3*(iatom-1)+idir
            tx(j) = x(iatom, idir)
        enddo
    enddo
    tg = 0.d0

    call n4pes(tx,v,tg,igrad)

    v  = v/23.0609 !+Eref*27.211386d0
    tg = tg/23.0609

    do istate=1,nstates
        p(istate) = v
    enddo

    do istate=1,nstates
        do iatom=1,natoms
            do idir=1,3
                j                    = 3*(iatom-1)+idir
                g(istate,iatom,idir) = tg(j)
            enddo
        enddo
    enddo

    d = 0.d0

end subroutine pes

!
! Local variables used in the N4 PES subroutine
! input coordinate matrix: X(12)in Ang
!                          N1: X(1),X(2),X(3)
!                          N2: X(4),X(5),X(6)
!                          N3: X(7),X(8),X(9)
!                          N4: X(10),X(11),X(12)
! input flag: igrad     igrad=0 energy-only calculation
!                       igrad=1 energy + gradient
! output potential energy:      v    in kcal/mol
! output gradient:              dVdX in kcal/(mol*Ang)
!
!> @brief Evaluate n4pes for the N4_1A_PIP_2013 potential.
subroutine n4pes( x, v, dvdx, igrad )
    !**********************************************************************
    ! Subroutine to calculate the potential energy V and gradient dVdX
    ! for given Cartesian coordinates X(12)
    ! R:            Interatomic bond distance (6)
    ! V:            Calculated potential energy
    ! dVdX:         The derivative of V w.r.t. X, dim(12)
    ! dVdR:         The derivative of V w.r.t. R, dim(6)
    ! dPdR:         The derivative of basis functions w.r.t. R
    !               dim(6*306)
    ! dMdR:         The derivative of monomials w.r.t. R
    !               dim(6*112)
    ! dRdX:         The derivative of R w.r.t. X, dim(6*12)
    !**********************************************************************

    implicit none


    integer          :: i,igrad,j,nob,k
    double precision :: v
    double precision :: dvdx(12),x(12)

    !      common /coord/    R(6)
    !      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
    !      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
    !     $                  dRdX(6,12),dBdR(6,276)
    !      common /msprmt/   a,re

    ! Read Cartesian coordinate from input file
    call coord_convt(x)

    if (igrad .le. 1) then
        ! Call subroutine Evv to evaluate potential energy V
        call evv(v)

        if (igrad .eq. 1) then
            ! Call EvdVdX to evaluate the derivatives of V w.r.t. X
            call evdvdx(x,dvdx)
        endif
    else
        write (*,*) 'Only energy and gradient are available'
    endif

end subroutine n4pes
!> @brief Evaluate coord_convt for the N4_1A_PIP_2013 potential.
subroutine coord_convt( x )
    !**********************************************************************
    !  Program to calculate the six interatomic distance
    !  by reading XYZ coordinate
    !**********************************************************************
    implicit none

    double precision :: r

    integer          :: i
    double precision :: x(12)

    common /coord/    r(6)

    !**********************************************************************
    !  Now, calculate the inter-atomic distance
    !  r1 = r(N1N2)    r2 = r(N1N3)
    !  r3 = r(N1N4)    r4 = r(N2N3)
    !  r5 = r(N2N4)    r6 = r(N3N4)
    !**********************************************************************

    r(1) = sqrt((x(4)-x(1))**2 + (x(5)-x(2))**2 + (x(6)-x(3))**2)
    r(2) = sqrt((x(7)-x(1))**2 + (x(8)-x(2))**2 + (x(9)-x(3))**2)
    r(3) = sqrt((x(10)-x(1))**2 + (x(11)-x(2))**2 + (x(12)-x(3))**2)
    r(4) = sqrt((x(4)-x(7))**2 + (x(5)-x(8))**2 + (x(6)-x(9))**2)
    r(5) = sqrt((x(4)-x(10))**2 + (x(5)-x(11))**2 + (x(6)-x(12))**2)
    r(6) = sqrt((x(7)-x(10))**2 + (x(8)-x(11))**2 + (x(9)-x(12))**2)


end subroutine coord_convt
!> @brief Evaluate EvV for the N4_1A_PIP_2013 potential.
subroutine evv( v )
    !**********************************************************************
    ! Subroutine to evaluate V for given R
    ! V(R) = C*P
    ! C:            Coefficients, stored in 'dim.inc'
    ! P:            Basis functions evaluated for given R
    ! rMs:          rMs(6), six Morse terms
    ! a:            Nonlinear parameters in Morse terrms(Angstrom)
    ! re:           Equilibrium bond length(Angstrom)
    ! nop:          number of points
    ! nom:          number of monomials
    ! nob:          number of basis functions(polynomials)
    ! rM(0:111):    Array to store monomials
    ! P(0:305):     Array to store polynomials
    ! B(1:276):     Array to store basis functions
    !**********************************************************************

    implicit none

    double precision :: b, c, p, r, rm, rms

    integer          :: i,j,k
    double precision :: dist,dv2dr,v,v2

    common /coord/    r(6)
    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)

    ! Calculate the six Morse terms for each point
    call evmorse

    ! Calculate the monomials for each point by using six Morse terms
    call evmono

    ! Calculate the polynomials (basis functions) by using monomials
    call evpoly

    ! Calculate the basis functions by removing unconnected and 2-body terms
    call evbas

    ! Initialized v to be 2De 2*228.7 kcal/mol
    v = 457.4d0
    ! Evaluate 2-body interactions
    do i=1,6
        dist = r(i)
        call ev2gm2(dist,v2,dv2dr,1,0)
        v = v+v2
    enddo

    ! Evaluate V by taken the product of C and Basis function array
    do i=1,276
        v = v + c(i)*b(i)
    enddo

    !      Write(*,9999) V
    ! 9999 Format('The potential energy is ',F20.14,' kcal/mol')


end subroutine evv
!> @brief Evaluate EvdVdX for the N4_1A_PIP_2013 potential.
subroutine evdvdx( x, dvdx )
    !**********************************************************************
    ! Subroutine to evaluate dRdX for given R and X
    ! R:            R(6), 6 bond lengths
    ! X:            X(12), 12 Cartesian coordinates
    ! rM(0:111):    Array to store monomials
    ! P(0:305):     Array to store polynomials
    ! dVdX:         dVdX(12), derivatives of V w.r.t. Cartesian coordinates
    ! dVdR:         dVdR(6), derivatives of V w.r.t.6 bond length
    ! dRdX:         dRdX(6,12), derivatives of R(6) w.r.t. 12
    !               Cartesian coordinates
    !**********************************************************************

    implicit none

    double precision :: b, c, dbdr, dmdr, dmsdr, dpdr, drdx, dvdr, p, r, rm, rms

    integer          :: i,j
    double precision :: dvdx(12),x(12)

    common /coord/    r(6)
    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)
    common /gradt/    dmsdr(6,6),dmdr(6,0:111),dpdr(6,0:305),dvdr(6), &
        & drdx(6,12),dbdr(6,276)

    ! Initialize dVdX
    do i=1,12
        dvdx(i) = 0.0d0
    enddo

    ! Call EvdVdR to evaluate dVdR(6)
    call evdvdr

    ! Call EvdRdX to evaluate dRdX(6,12)
    call evdrdx(x)

    ! Calculate dVdX by using chain rule: dV/dXi=(dV/dRj)*(dRj/dXi), j=1 to
    do i=1,12
        do j=1,6
            dvdx(i) = dvdx(i) + dvdr(j)*drdx(j,i)
        enddo
    enddo

    !      write(*,*) 'The 12 dVdX are:'
    !      Write(*,9999) (dVdX(i),i=1,12)
    ! 9999 Format(1x,3F15.8)

end subroutine evdvdx
!> @brief Evaluate EvMorse for the N4_1A_PIP_2013 potential.
subroutine evmorse

    implicit none

    double precision :: a, b, c, p, r, re, rm, rms

    integer :: i

    common /coord/    r(6)
    common /msprmt/   a,re
    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)

    ! Morse term ms = exp(-(r-re)/a)
    ! re:   equlibrium bond length(1.098A)
    ! a:    nonlinear parameter optimized for diatomic molecules, unit Anstr
    do i=1,6
        rms(i) = exp(-(r(i)-re)/a)
    enddo

end subroutine evmorse
!> @brief Evaluate EvMono for the N4_1A_PIP_2013 potential.
subroutine evmono
    !**********************************************************************
    !  The subroutine reads six Morse variables(X) and calculates the
    !  monomials(M) that do not have usable decomposition.
    !  For A4 with max. degree 9, the number of monomials is nom.
    !**********************************************************************

    implicit none

    double precision :: b, c, p, rm, rms

    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)

    rm(0)   = 1.0d0
    rm(1)   = rms(6)
    rm(2)   = rms(5)
    rm(3)   = rms(4)
    rm(4)   = rms(3)
    rm(5)   = rms(2)
    rm(6)   = rms(1)
    rm(7)   = rm(3)*rm(4)
    rm(8)   = rm(2)*rm(5)
    rm(9)   = rm(1)*rm(6)
    rm(10)  = rm(1)*rm(2)
    rm(11)  = rm(1)*rm(3)
    rm(12)  = rm(2)*rm(3)
    rm(13)  = rm(1)*rm(4)
    rm(14)  = rm(2)*rm(4)
    rm(15)  = rm(1)*rm(5)
    rm(16)  = rm(3)*rm(5)
    rm(17)  = rm(4)*rm(5)
    rm(18)  = rm(2)*rm(6)
    rm(19)  = rm(3)*rm(6)
    rm(20)  = rm(4)*rm(6)
    rm(21)  = rm(5)*rm(6)
    rm(22)  = rm(1)*rm(7)
    rm(23)  = rm(2)*rm(7)
    rm(24)  = rm(1)*rm(8)
    rm(25)  = rm(2)*rm(16)
    rm(26)  = rm(2)*rm(17)
    rm(27)  = rm(3)*rm(17)
    rm(28)  = rm(1)*rm(18)
    rm(29)  = rm(1)*rm(19)
    rm(30)  = rm(1)*rm(20)
    rm(31)  = rm(3)*rm(20)
    rm(32)  = rm(1)*rm(21)
    rm(33)  = rm(2)*rm(21)
    rm(34)  = rm(1)*rm(12)
    rm(35)  = rm(1)*rm(17)
    rm(36)  = rm(2)*rm(20)
    rm(37)  = rm(3)*rm(21)
    rm(38)  = rm(1)*rm(14)
    rm(39)  = rm(1)*rm(16)
    rm(40)  = rm(2)*rm(19)
    rm(41)  = rm(4)*rm(21)
    rm(42)  = rm(2)*rm(27)
    rm(43)  = rm(1)*rm(31)
    rm(44)  = rm(1)*rm(33)
    rm(45)  = rm(1)*rm(23)
    rm(46)  = rm(1)*rm(25)
    rm(47)  = rm(1)*rm(26)
    rm(48)  = rm(1)*rm(27)
    rm(49)  = rm(1)*rm(40)
    rm(50)  = rm(1)*rm(36)
    rm(51)  = rm(2)*rm(31)
    rm(52)  = rm(1)*rm(37)
    rm(53)  = rm(2)*rm(37)
    rm(54)  = rm(1)*rm(41)
    rm(55)  = rm(2)*rm(41)
    rm(56)  = rm(3)*rm(41)
    rm(57)  = rm(1)*rm(42)
    rm(58)  = rm(1)*rm(51)
    rm(59)  = rm(1)*rm(53)
    rm(60)  = rm(1)*rm(55)
    rm(61)  = rm(1)*rm(56)
    rm(62)  = rm(2)*rm(56)
    rm(63)  = rm(1)*rm(62)
    rm(64)  = rm(2)*rm(57)
    rm(65)  = rm(3)*rm(57)
    rm(66)  = rm(4)*rm(57)
    rm(67)  = rm(5)*rm(57)
    rm(68)  = rm(1)*rm(58)
    rm(69)  = rm(3)*rm(58)
    rm(70)  = rm(4)*rm(58)
    rm(71)  = rm(1)*rm(59)
    rm(72)  = rm(2)*rm(59)
    rm(73)  = rm(1)*rm(60)
    rm(74)  = rm(2)*rm(60)
    rm(75)  = rm(1)*rm(61)
    rm(76)  = rm(2)*rm(62)
    rm(77)  = rm(3)*rm(61)
    rm(78)  = rm(3)*rm(62)
    rm(79)  = rm(4)*rm(61)
    rm(80)  = rm(4)*rm(62)
    rm(81)  = rm(5)*rm(59)
    rm(82)  = rm(5)*rm(60)
    rm(83)  = rm(5)*rm(62)
    rm(84)  = rm(6)*rm(58)
    rm(85)  = rm(6)*rm(59)
    rm(86)  = rm(6)*rm(60)
    rm(87)  = rm(6)*rm(61)
    rm(88)  = rm(2)*rm(64)
    rm(89)  = rm(3)*rm(65)
    rm(90)  = rm(4)*rm(66)
    rm(91)  = rm(5)*rm(67)
    rm(92)  = rm(1)*rm(68)
    rm(93)  = rm(3)*rm(69)
    rm(94)  = rm(4)*rm(70)
    rm(95)  = rm(1)*rm(71)
    rm(96)  = rm(2)*rm(72)
    rm(97)  = rm(1)*rm(73)
    rm(98)  = rm(2)*rm(74)
    rm(99)  = rm(1)*rm(75)
    rm(100) = rm(2)*rm(76)
    rm(101) = rm(3)*rm(77)
    rm(102) = rm(3)*rm(78)
    rm(103) = rm(4)*rm(79)
    rm(104) = rm(4)*rm(80)
    rm(105) = rm(5)*rm(81)
    rm(106) = rm(5)*rm(82)
    rm(107) = rm(5)*rm(83)
    rm(108) = rm(6)*rm(84)
    rm(109) = rm(6)*rm(85)
    rm(110) = rm(6)*rm(86)
    rm(111) = rm(6)*rm(87)



end subroutine evmono
!> @brief Evaluate EvPoly for the N4_1A_PIP_2013 potential.
subroutine evpoly
    !**********************************************************************
    !  The subroutine reads monomials(m) and calculates the
    !  permutation invariant polynomials(p)
    !  For A4 with max. degree 9, the number of polynomials is nob.
    !**********************************************************************

    implicit none

    double precision :: b, c, p, rm, rms

    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)

    p(0)  = rm(0)
    p(1)  = rm(1) + rm(2) + rm(3) + rm(4) + rm(5) + rm(6)
    p(2)  = rm(7) + rm(8) + rm(9)
    p(3)  = rm(10) + rm(11) + rm(12) + rm(13) + rm(14) + rm(15) + rm(16) + rm(17) + rm(18) + rm(19) + rm(20) + rm(21)
    p(4)  = p(1)*p(1) - p(3) - p(2) - p(3) - p(2)
    p(5)  = rm(22) + rm(23) + rm(24) + rm(25) + rm(26) + rm(27) + rm(28) + rm(29) + rm(30) + rm(31) + rm(32) + rm(33)
    p(6)  = rm(34) + rm(35) + rm(36) + rm(37)
    p(7)  = rm(38) + rm(39) + rm(40) + rm(41)
    p(8)  = p(1)*p(2) - p(5)
    p(9)  = p(1)*p(3) - p(6) - p(7) - p(5) - p(6) - p(7) - p(5) - p(6) - p(7)
    p(10) = p(1)*p(4) - p(9) - p(8)
    p(11) = rm(42) + rm(43) + rm(44)
    p(12) = rm(45) + rm(46) + rm(47) + rm(48) + rm(49) + rm(50) + rm(51) + rm(52) + rm(53) + rm(54) + rm(55) + rm(56)
    p(13) = p(2)*p(3) - p(12)
    p(14) = p(1)*p(5) - p(12) - p(11) - p(13) - p(12) - p(11) - p(11) - p(11)
    p(15) = p(1)*p(6) - p(12)
    p(16) = p(1)*p(7) - p(12)
    p(17) = p(2)*p(2) - p(11) - p(11)
    p(18) = p(3)*p(3) - p(12) - p(11) - p(15) - p(16) - p(14) - p(12) - p(11) - p(15) - p(16) - p(14) - p(12) - &
        & p(11) - p(12) - p(11)
    p(19) = p(2)*p(4) - p(14)
    p(20) = p(3)*p(4) - p(15) - p(16) - p(13)
    p(21) = p(1)*p(10) - p(20) - p(19)
    p(22) = rm(57) + rm(58) + rm(59) + rm(60) + rm(61) + rm(62)
    p(23) = p(1)*p(11) - p(22)
    p(24) = p(2)*p(6)
    p(25) = p(2)*p(7)
    p(26) = p(1)*p(12) - p(22) - p(24) - p(25) - p(22) - p(22) - p(22)
    p(27) = p(2)*p(5) - p(22) - p(23) - p(22)
    p(28) = p(3)*p(5) - p(22) - p(26) - p(24) - p(25) - p(23) - p(22) - p(24) - p(25) - p(23) - p(22) - p(22)
    p(29) = p(3)*p(6) - p(22) - p(26) - p(22)
    p(30) = p(3)*p(7) - p(22) - p(26) - p(22)
    p(31) = p(2)*p(9) - p(26) - p(28)
    p(32) = p(1)*p(14) - p(26) - p(23) - p(28)
    p(33) = p(4)*p(6) - p(25)
    p(34) = p(4)*p(7) - p(24)
    p(35) = p(1)*p(17) - p(27)
    p(36) = p(1)*p(18) - p(29) - p(30) - p(28)
    p(37) = p(2)*p(10) - p(32)
    p(38) = p(3)*p(10) - p(33) - p(34) - p(31)
    p(39) = p(1)*p(21) - p(38) - p(37)
    p(40) = rm(63)
    p(41) = rm(64) + rm(65) + rm(66) + rm(67) + rm(68) + rm(69) + rm(70) + rm(71) + rm(72) + rm(73) + rm(74) + &
        & rm(75) + rm(76) + rm(77) + rm(78) + rm(79) + rm(80) + rm(81) + rm(82) + rm(83) + rm(84) + rm(85) + rm(86) + &
        & rm(87)
    p(42) = p(1)*p(22) - p(40) - p(41) - p(40) - p(40) - p(40) - p(40) - p(40)
    p(43) = p(2)*p(11) - p(40) - p(40) - p(40)
    p(44) = p(2)*p(12) - p(41)
    p(45) = p(3)*p(11) - p(41)
    p(46) = p(5)*p(6) - p(41)
    p(47) = p(5)*p(7) - p(41)
    p(48) = p(6)*p(7) - p(40) - p(40) - p(40) - p(40)
    p(49) = p(4)*p(11) - p(42)
    p(50) = p(2)*p(15) - p(46)
    p(51) = p(2)*p(16) - p(47)
    p(52) = p(4)*p(12) - p(41) - p(50) - p(51)
    p(53) = p(2)*p(14) - p(42) - p(49) - p(42)
    p(54) = p(6)*p(6) - p(42) - p(42)
    p(55) = p(7)*p(7) - p(42) - p(42)
    p(56) = p(3)*p(17) - p(44)
    p(57) = p(2)*p(18) - p(48)
    p(58) = p(3)*p(14) - p(41) - p(52) - p(46) - p(47) - p(45) - p(45)
    p(59) = p(6)*p(9) - p(41) - p(52) - p(47)
    p(60) = p(7)*p(9) - p(41) - p(52) - p(46)
    p(61) = p(2)*p(20) - p(52) - p(58)
    p(62) = p(1)*p(32) - p(52) - p(49) - p(58)
    p(63) = p(6)*p(10) - p(51)
    p(64) = p(7)*p(10) - p(50)
    p(65) = p(2)*p(17) - p(43)
    p(66) = p(3)*p(18) - p(46) - p(47) - p(45) - p(59) - p(60) - p(58)
    p(67) = p(2)*p(19) - p(49)
    p(68) = p(1)*p(36) - p(59) - p(60) - p(58) - p(57) - p(66) - p(66)
    p(69) = p(2)*p(21) - p(62)
    p(70) = p(3)*p(21) - p(63) - p(64) - p(61)
    p(71) = p(1)*p(39) - p(70) - p(69)
    p(72) = p(40)*p(1)
    p(73) = p(2)*p(22) - p(72)
    p(74) = p(6)*p(11)
    p(75) = p(7)*p(11)
    p(76) = p(3)*p(22) - p(72) - p(74) - p(75) - p(72) - p(72) - p(72)
    p(77) = rm(88) + rm(89) + rm(90) + rm(91) + rm(92) + rm(93) + rm(94) + rm(95) + rm(96) + rm(97) + rm(98) + &
        & rm(99) + rm(100) + rm(101) + rm(102) + rm(103) + rm(104) + rm(105) + rm(106) + rm(107) + rm(108) + rm(109) + &
        & rm(110) + rm(111)
    p(78)  = p(1)*p(42) - p(72) - p(76)
    p(79)  = p(5)*p(11) - p(72) - p(73) - p(72)
    p(80)  = p(2)*p(26) - p(76) - p(77)
    p(81)  = p(6)*p(12) - p(72) - p(76) - p(72)
    p(82)  = p(7)*p(12) - p(72) - p(76) - p(72)
    p(83)  = p(8)*p(11) - p(72)
    p(84)  = p(6)*p(17)
    p(85)  = p(7)*p(17)
    p(86)  = p(9)*p(11) - p(76) - p(77)
    p(87)  = p(2)*p(29) - p(81)
    p(88)  = p(6)*p(14) - p(75) - p(75)
    p(89)  = p(2)*p(30) - p(82)
    p(90)  = p(7)*p(14) - p(74) - p(74)
    p(91)  = p(1)*p(48) - p(76) - p(81) - p(82)
    p(92)  = p(10)*p(11) - p(78)
    p(93)  = p(2)*p(33) - p(88)
    p(94)  = p(2)*p(34) - p(90)
    p(95)  = p(10)*p(12) - p(77) - p(93) - p(94)
    p(96)  = p(2)*p(27) - p(73) - p(79)
    p(97)  = p(2)*p(28) - p(76) - p(86)
    p(98)  = p(1)*p(53) - p(80) - p(79) - p(97)
    p(99)  = p(1)*p(54) - p(81)
    p(100) = p(1)*p(55) - p(82)
    p(101) = p(5)*p(18) - p(76) - p(91) - p(87) - p(89) - p(86)
    p(102) = p(6)*p(18) - p(74) - p(90)
    p(103) = p(7)*p(18) - p(75) - p(88)
    p(104) = p(2)*p(31) - p(77) - p(86)
    p(105) = p(2)*p(36) - p(91) - p(101)
    p(106) = p(3)*p(32) - p(77) - p(95) - p(88) - p(90) - p(86)
    p(107) = p(4)*p(29) - p(82) - p(80) - p(99)
    p(108) = p(4)*p(30) - p(81) - p(80) - p(100)
    p(109) = p(2)*p(38) - p(95) - p(106)
    p(110) = p(1)*p(62) - p(95) - p(92) - p(106)
    p(111) = p(6)*p(21) - p(94)
    p(112) = p(7)*p(21) - p(93)
    p(113) = p(1)*p(65) - p(96)
    p(114) = p(1)*p(66) - p(102) - p(103) - p(101)
    p(115) = p(2)*p(37) - p(92)
    p(116) = p(10)*p(18) - p(99) - p(100) - p(97)
    p(117) = p(2)*p(39) - p(110)
    p(118) = p(3)*p(39) - p(111) - p(112) - p(109)
    p(119) = p(1)*p(71) - p(118) - p(117)
    p(120) = p(40)*p(2)
    p(121) = p(40)*p(3)
    p(122) = p(40)*p(4)
    p(123) = p(11)*p(12) - p(121)
    p(124) = p(2)*p(42) - p(122)
    p(125) = p(6)*p(22) - p(121)
    p(126) = p(7)*p(22) - p(121)
    p(127) = p(2)*p(41) - p(121) - p(123) - p(121)
    p(128) = p(6)*p(23) - p(123)
    p(129) = p(7)*p(23) - p(123)
    p(130) = p(3)*p(41) - p(122) - p(121) - p(120) - p(125) - p(128) - p(126) - p(129) - p(124) - p(123) - p(122) - &
        & p(121) - p(120) - p(125) - p(126) - p(124) - p(123) - p(122) - p(121) - p(120) - p(122) - p(121) - p(120) - &
        & p(120) - p(120) - p(120) - p(120)
    p(131) = p(3)*p(42) - p(121) - p(125) - p(126) - p(121)
    p(132) = p(4)*p(41) - p(121) - p(131) - p(128) - p(129) - p(127) - p(121)
    p(133) = p(1)*p(78) - p(122) - p(131)
    p(134) = p(11)*p(11) - p(120) - p(120)
    p(135) = p(2)*p(48) - p(130)
    p(136) = p(11)*p(17) - p(120)
    p(137) = p(2)*p(44) - p(123)
    p(138) = p(2)*p(45) - p(121)
    p(139) = p(11)*p(14) - p(122) - p(124) - p(122)
    p(140) = p(6)*p(27) - p(127)
    p(141) = p(2)*p(54)
    p(142) = p(7)*p(27) - p(127)
    p(143) = p(2)*p(52) - p(131) - p(132)
    p(144) = p(1)*p(81) - p(125) - p(141) - p(135) - p(125)
    p(145) = p(2)*p(55)
    p(146) = p(1)*p(82) - p(126) - p(135) - p(145) - p(126)
    p(147) = p(11)*p(18) - p(130)
    p(148) = p(6)*p(28) - p(129) - p(123) - p(143)
    p(149) = p(7)*p(28) - p(128) - p(123) - p(143)
    p(150) = p(6)*p(30) - p(121) - p(146)
    p(151) = p(11)*p(19) - p(122)
    p(152) = p(2)*p(50) - p(128)
    p(153) = p(2)*p(51) - p(129)
    p(154) = p(11)*p(20) - p(131) - p(132)
    p(155) = p(2)*p(59) - p(144) - p(148)
    p(156) = p(6)*p(32) - p(129)
    p(157) = p(2)*p(60) - p(146) - p(149)
    p(158) = p(7)*p(32) - p(128)
    p(159) = p(6)*p(34) - p(122) - p(145) - p(122)
    p(160) = p(11)*p(21) - p(133)
    p(161) = p(2)*p(63) - p(156)
    p(162) = p(2)*p(64) - p(158)
    p(163) = p(12)*p(21) - p(132) - p(161) - p(162)
    p(164) = p(2)*p(53) - p(124) - p(139)
    p(165) = p(2)*p(58) - p(131) - p(154)
    p(166) = p(3)*p(54) - p(125) - p(144)
    p(167) = p(3)*p(55) - p(126) - p(146)
    p(168) = p(3)*p(65) - p(137)
    p(169) = p(17)*p(18) - p(135)
    p(170) = p(1)*p(98) - p(143) - p(139) - p(165)
    p(171) = p(4)*p(54) - p(135)
    p(172) = p(4)*p(55) - p(135)
    p(173) = p(2)*p(66) - p(150)
    p(174) = p(1)*p(101) - p(148) - p(149) - p(147) - p(173) - p(165) - p(147)
    p(175) = p(1)*p(102) - p(150) - p(148) - p(166)
    p(176) = p(1)*p(103) - p(150) - p(149) - p(167)
    p(177) = p(2)*p(61) - p(132) - p(154)
    p(178) = p(2)*p(68) - p(159) - p(174)
    p(179) = p(3)*p(62) - p(132) - p(163) - p(156) - p(158) - p(154)
    p(180) = p(6)*p(38) - p(132) - p(163) - p(157)
    p(181) = p(7)*p(38) - p(132) - p(163) - p(155)
    p(182) = p(2)*p(70) - p(163) - p(179)
    p(183) = p(1)*p(110) - p(163) - p(160) - p(179)
    p(184) = p(6)*p(39) - p(162)
    p(185) = p(7)*p(39) - p(161)
    p(186) = p(2)*p(65) - p(136)
    p(187) = p(3)*p(66) - p(148) - p(149) - p(147) - p(175) - p(176) - p(174)
    p(188) = p(2)*p(67) - p(151)
    p(189) = p(4)*p(66) - p(166) - p(167) - p(165)
    p(190) = p(2)*p(69) - p(160)
    p(191) = p(18)*p(21) - p(171) - p(172) - p(169)
    p(192) = p(2)*p(71) - p(183)
    p(193) = p(3)*p(71) - p(184) - p(185) - p(182)
    p(194) = p(1)*p(119) - p(193) - p(192)
    p(195) = p(40)*p(5)
    p(196) = p(40)*p(6)
    p(197) = p(40)*p(7)
    p(198) = p(40)*p(8)
    p(199) = p(40)*p(9)
    p(200) = p(40)*p(10)
    p(201) = p(11)*p(22) - p(195)
    p(202) = p(12)*p(22) - p(196) - p(197) - p(195) - p(196) - p(197) - p(195) - p(196) - p(197)
    p(203) = p(17)*p(22) - p(198)
    p(204) = p(6)*p(43)
    p(205) = p(7)*p(43)
    p(206) = p(11)*p(26) - p(199) - p(202)
    p(207) = p(2)*p(76) - p(199) - p(202)
    p(208) = p(2)*p(78) - p(200)
    p(209) = p(6)*p(41) - p(199) - p(195) - p(202) - p(195)
    p(210) = p(6)*p(42) - p(197) - p(197) - p(197)
    p(211) = p(7)*p(41) - p(199) - p(195) - p(202) - p(195)
    p(212) = p(7)*p(42) - p(196) - p(196) - p(196)
    p(213) = p(11)*p(29) - p(209)
    p(214) = p(11)*p(30) - p(211)
    p(215) = p(18)*p(22) - p(199) - p(213) - p(214)
    p(216) = p(2)*p(77) - p(199) - p(206)
    p(217) = p(6)*p(49) - p(205)
    p(218) = p(7)*p(49) - p(204)
    p(219) = p(3)*p(77) - p(200) - p(199) - p(198) - p(209) - p(217) - p(211) - p(218) - p(207) - p(204) - p(205) - &
        & p(200) - p(199) - p(198) - p(200) - p(198) - p(200) - p(198)
    p(220) = p(3)*p(78) - p(199) - p(210) - p(212)
    p(221) = p(10)*p(41) - p(199) - p(220) - p(217) - p(218) - p(216)
    p(222) = p(1)*p(133) - p(200) - p(220)
    p(223) = p(11)*p(27) - p(195) - p(203)
    p(224) = p(1)*p(134) - p(201)
    p(225) = p(2)*p(81) - p(209)
    p(226) = p(2)*p(80) - p(202) - p(206)
    p(227) = p(2)*p(82) - p(211)
    p(228) = p(2)*p(91) - p(215) - p(219)
    p(229) = p(11)*p(28) - p(199) - p(207)
    p(230) = p(6)*p(53) - p(205)
    p(231) = p(5)*p(54) - p(209)
    p(232) = p(7)*p(53) - p(204)
    p(233) = p(7)*p(54) - p(196)
    p(234) = p(5)*p(55) - p(211)
    p(235) = p(6)*p(55) - p(197)
    p(236) = p(11)*p(35) - p(198)
    p(237) = p(6)*p(65)
    p(238) = p(7)*p(65)
    p(239) = p(2)*p(86) - p(199) - p(229)
    p(240) = p(1)*p(139) - p(206) - p(229) - p(224)
    p(241) = p(17)*p(29) - p(225)
    p(242) = p(2)*p(99) - p(231)
    p(243) = p(17)*p(30) - p(227)
    p(244) = p(2)*p(95) - p(220) - p(221)
    p(245) = p(4)*p(81) - p(202) - p(242) - p(227)
    p(246) = p(2)*p(100) - p(234)
    p(247) = p(4)*p(82) - p(202) - p(225) - p(246)
    p(248) = p(11)*p(36) - p(215) - p(219)
    p(249) = p(2)*p(102) - p(233)
    p(250) = p(6)*p(58) - p(206) - p(214) - p(244) - p(214)
    p(251) = p(2)*p(103) - p(235)
    p(252) = p(7)*p(58) - p(213) - p(206) - p(244) - p(213)
    p(253) = p(1)*p(150) - p(215) - p(233) - p(235)
    p(254) = p(11)*p(37) - p(200)
    p(255) = p(2)*p(93) - p(217)
    p(256) = p(2)*p(94) - p(218)
    p(257) = p(11)*p(38) - p(220) - p(221)
    p(258) = p(2)*p(107) - p(245) - p(250)
    p(259) = p(6)*p(62) - p(218)
    p(260) = p(2)*p(108) - p(247) - p(252)
    p(261) = p(7)*p(62) - p(217)
    p(262) = p(6)*p(64) - p(200) - p(246) - p(200)
    p(263) = p(11)*p(39) - p(222)
    p(264) = p(2)*p(111) - p(259)
    p(265) = p(2)*p(112) - p(261)
    p(266) = p(12)*p(39) - p(221) - p(264) - p(265)
    p(267) = p(2)*p(98) - p(208) - p(240)
    p(268) = p(6)*p(54) - p(210)
    p(269) = p(7)*p(55) - p(212)
    p(270) = p(2)*p(96) - p(203) - p(223)
    p(271) = p(2)*p(97) - p(207) - p(229)
    p(272) = p(2)*p(101) - p(215) - p(248)
    p(273) = p(2)*p(106) - p(220) - p(257)
    p(274) = p(6)*p(59) - p(220) - p(215) - p(209)
    p(275) = p(7)*p(60) - p(220) - p(215) - p(211)
    p(276) = p(5)*p(66) - p(215) - p(253) - p(249) - p(251) - p(248)
    p(277) = p(6)*p(66) - p(213) - p(252)
    p(278) = p(7)*p(66) - p(214) - p(250)
    p(279) = p(2)*p(104) - p(216) - p(239)
    p(280) = p(2)*p(105) - p(219) - p(248)
    p(281) = p(1)*p(170) - p(244) - p(240) - p(273)
    p(282) = p(10)*p(54) - p(227)
    p(283) = p(10)*p(55) - p(225)
    p(284) = p(2)*p(114) - p(253) - p(276)
    p(285) = p(1)*p(174) - p(250) - p(252) - p(248) - p(276) - p(273)
    p(286) = p(6)*p(68) - p(217) - p(261) - p(251)
    p(287) = p(7)*p(68) - p(218) - p(259) - p(249)
    p(288) = p(2)*p(109) - p(221) - p(257)
    p(289) = p(2)*p(116) - p(262) - p(285)
    p(290) = p(3)*p(110) - p(221) - p(266) - p(259) - p(261) - p(257)
    p(291) = p(6)*p(70) - p(221) - p(266) - p(260)
    p(292) = p(7)*p(70) - p(221) - p(266) - p(258)
    p(293) = p(2)*p(118) - p(266) - p(290)
    p(294) = p(1)*p(183) - p(266) - p(263) - p(290)
    p(295) = p(6)*p(71) - p(265)
    p(296) = p(7)*p(71) - p(264)
    p(297) = p(1)*p(186) - p(270)
    p(298) = p(1)*p(187) - p(277) - p(278) - p(276)
    p(299) = p(2)*p(115) - p(254)
    p(300) = p(1)*p(189) - p(286) - p(287) - p(285) - p(284) - p(298)
    p(301) = p(2)*p(117) - p(263)
    p(302) = p(18)*p(39) - p(282) - p(283) - p(280)
    p(303) = p(2)*p(119) - p(294)
    p(304) = p(3)*p(119) - p(295) - p(296) - p(293)
    p(305) = p(1)*p(194) - p(304) - p(303)


end subroutine evpoly
!> @brief Evaluate evbas for the N4_1A_PIP_2013 potential.
subroutine evbas
    !**********************************************************************
    !  The subroutine eliminate the 2-body terms in Bowman's approach
    !**********************************************************************

    implicit none

    double precision :: b, c, p, rm, rms

    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)

    integer          :: i
    double precision :: b1(306)

    ! Pass P(0:305) to BM1(1:306)
    do i=1,306
        b1(i) = p(i-1)
    enddo


    ! Remove unconnected terms and 2-body terms and pass to B(1:276)
    b(1) = b1(4)

    do i=2,4
        b(i) = b1(i+4)
    enddo

    b(5) = b1(10)

    do i=6,11
        b(i) = b1(i+6)
    enddo

    b(12) = b1(19)
    b(13) = b1(21)

    do i=14,26
        b(i) = b1(i+9)
    enddo

    b(27) = b1(37)
    b(28) = b1(39)

    do i=29,53
        b(i) = b1(i+12)
    enddo

    b(54) = b1(67)
    b(55) = b1(69)
    b(56) = b1(71)

    do i=57,97
        b(i) = b1(i+16)
    enddo

    b(98)  = b1(115)
    b(99)  = b1(117)
    b(100) = b1(119)

    do i=101,166
        b(i) = b1(i+20)
    enddo

    b(167) = b1(188)
    b(168) = b1(190)
    b(169) = b1(192)
    b(170) = b1(194)

    do i=171,272
        b(i) = b1(i+25)
    enddo

    b(273) = b1(299)
    b(274) = b1(301)
    b(275) = b1(303)
    b(276) = b1(305)


end subroutine evbas
!> @brief Evaluate ev2gm2 for the N4_1A_PIP_2013 potential.
subroutine ev2gm2( r, v2, dv2dr, imol, igrad )
    !**********************************************************************
    !
    ! Subroutine to evaluate the 2-body potential energy and gradient
    ! for given r with generalized Morse potential with Scheme 2.
    ! V(r) = 0 for r -> inifity
    ! V(r) = De*(1-exp(-f(y)*(r-re)))^2 - De
    ! f(y) = c0 + c1*y + c2*y^2 + c3*y^3 + c4*y^4 + c5*y^5 + c6*y^6
    ! y = (r^4 - re^4)/(r^4 + re^4)
    !
    ! imol .eq. 1   Paramemters for N2 dissociation without SEC
    !               (for Yuliya's N4 data)
    ! imol .eq. 2   Paramemters for N2 dissociation with SEC
    !               (for Zoltan's N2O2 data)
    ! imol .eq. 3   Paramemters for NO dissociation with SEC
    !               (for Zoltan's N2O2 data)
    ! imol .eq. 4   Paramemters for O2 dissociation with SEC
    !               (for Zoltan's N2O2 data and future O4 data)
    !**********************************************************************

    implicit none

    double precision :: dfdy, v

    integer          :: igrad,imol
    double precision :: c(0:10)
    double precision :: re,de
    double precision :: fy,dfdr,r,u,v2,dv2dr,y,dydr

    if (imol .eq. 1) then
        ! Parameter for N2 dissociation without SEC (for Yuliya's N4 data)
        re   = 1.098d0
        de   = 228.7d0
        c(0) = 2.70963254293d0
        c(1) = 1.32620177271d-1
        c(2) = 2.96757048793d-1
        c(3) = 1.97112432229d-1
        c(4) = -5.02002309588d-1
        c(5) = 3.80734244606d-1
        c(6) = 1.21001628750d0
    else if (imol .eq. 2) then
        ! Parameter for N2 dissociation with SEC (for Zoltan's N2O2)
        re   = 1.098d0
        de   = 228.4d0
        c(0) = 2.71405774451d0
        c(1) = 1.32757649829d-1
        c(2) = 2.66756890408d-1
        c(3) = 1.95350725241d-1
        c(4) = -4.08663480982d-1
        c(5) = 3.92451705557d-1
        c(6) = 1.13006674877d0
    else if (imol .eq. 3) then
        ! Parameter for NO dissociation with SEC (for Zoltan's N2O2)
        re   = 1.1508d0
        de   = 149.9d0
        c(0) = 2.81134495569d0
        c(1) = 1.43241169611d-1
        c(2) = 1.35760038863d-2
        c(3) = 3.92892178507d-1
        c(4) = 9.29495534058d-1
        c(5) = 2.66966672332d-1
        c(6) = -3.68118714223d-1
    else if (imol .eq. 4) then
        ! Parameter for O2 dissociation with SEC (for Zoltan's N2O2 and KRY's O4
        re   = 1.208d0
        de   = 120.243d0
        c(0) = 2.69132890094d0
        c(1) = 3.39550045614d-1
        c(2) = 3.46401777195d-1
        c(3) = -7.76983671636d-1
        c(4) = -3.29632972405d-1
        c(5) = 2.40883331247d0
        c(6) = 2.09264029009d0
    endif

    y = (r*r*r*r - re*re*re*re)/(r*r*r*r + re*re*re*re)

    fy = c(0) + c(1)*y + c(2)*y*y + c(3)*y*y*y + c(4)*y*y*y*y + c(5)*y*y*y*y*y + c(6)*y*y*y*y*y*y

    u = exp(-fy*(r-re))

    v2 = de*(1.0d0-u)*(1.0d0-u)-de

    if (igrad .eq. 1) then
        dfdy  = c(1) + 2.0d0*c(2)*y + 3.0d0*c(3)*y*y + 4.0d0*c(4)*y*y*y + 5.0d0*c(5)*y*y*y*y + 6.0d0*c(6)*y*y*y*y*y
        v     = r*r*r*r + re*re*re*re
        dydr  = 8.0d0*r*r*r*re*re*re*re/(v*v)
        dfdr  = dfdy*dydr
        dv2dr = 2.0d0*de*(1-u)*u*(dfdr*(r-re)+fy)
    endif

end subroutine ev2gm2
!> @brief Evaluate EvdVdR for the N4_1A_PIP_2013 potential.
subroutine evdvdr
    !**********************************************************************
    ! Subroutine to evaluate dVdR for given R
    ! dVdR = dV2dR + C*dBdR
    ! C:            Coefficients, stored in 'dim.inc'
    ! P:            Basis functions evaluated for given R
    ! M:            Monomials evaluated for given R
    ! dV2dR:        Gradient of 2-body interactions
    ! dMsdR:        dMsdR(6,6), 6 Morse terms w.r.t. 6 bond length
    ! dMdR:         dMdR(6,nom), nom monomials w.r.t.6 bond length
    ! dPdR:         dPdR(6,nob), nop polynomial basis functions
    !               w.r.t. 6 bond length
    ! nom:          number of monomials
    ! nob:          number of basis functions(polynomials)
    ! M(nom):       Array to store monomials
    ! P(nob):       Array to store polynomials
    !**********************************************************************

    implicit none

    double precision :: b, c, dbdr, dmdr, dmsdr, dpdr, drdx, dvdr, p, r, rm, rms

    integer          :: i,j
    double precision :: dist,v2,dv2dr

    common /coord/    r(6)
    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)
    common /gradt/    dmsdr(6,6),dmdr(6,0:111),dpdr(6,0:305),dvdr(6), &
        & drdx(6,12),dbdr(6,276)

    ! Initialize dVdR(6)
    do i=1,6
        dvdr(i) = 0.0d0
    enddo

    ! Add dV2dR(i) to dVdR
    do i=1,6
        dist = r(i)
        call ev2gm2(dist,v2,dv2dr,1,1)
        dvdr(i) = dv2dr
    enddo

    ! Calculate dMorse/dr(6,6) for given R(6)
    call evdmsdr

    ! Calculate the monomials for each point by using six Morse terms
    call evdmdr

    ! Calculate the polynomials by using monomials
    call evdpdr

    ! Remove 2-body interactions and unconnected terms from polynomials
    call evdbdr

    ! Evaluate dVdR(6) by taken the product of C(j) and dPdR(i,j)
    do i=1,6
        do j=1,276
            dvdr(i) = dvdr(i) + c(j)*dbdr(i,j)
        enddo
    enddo

end subroutine evdvdr
!> @brief Evaluate EvdRdX for the N4_1A_PIP_2013 potential.
subroutine evdrdx( x )
    !**********************************************************************
    ! Subroutine to evaluate dRdX for given R and X
    ! R:            R(6), 6 bond lengths
    ! X:            X(12), 12 Cartesian coordinates
    !
    ! dMdR:         dMdR(6,nom), nom monomials w.r.t.6 bond length
    ! dPdR:         dPdR(6,nob), nop polynomial basis functions
    !               w.r.t. 6 bond length
    ! M(nom):       Array to store monomials
    ! P(nob):       Array to store polynomials
    !**********************************************************************

    implicit none

    double precision :: dbdr, dmdr, dmsdr, dpdr, drdx, dvdr, r

    integer          :: i,j
    double precision :: x(12)

    common /coord/    r(6)
    common /gradt/    dmsdr(6,6),dmdr(6,0:111),dpdr(6,0:305),dvdr(6), &
        & drdx(6,12),dbdr(6,276)

    ! Initialize dRdX(6,12)
    do i=1,6
        do j=1,12
            drdx(i,j) = 0.0d0
        enddo
    enddo

    ! Start to calculate the non-zero dRdX
    ! dr1dx
    drdx(1,1) = (x(1)-x(4))/r(1)
    drdx(1,2) = (x(2)-x(5))/r(1)
    drdx(1,3) = (x(3)-x(6))/r(1)
    drdx(1,4) = -drdx(1,1)
    drdx(1,5) = -drdx(1,2)
    drdx(1,6) = -drdx(1,3)

    ! dr2dx
    drdx(2,1) = (x(1)-x(7))/r(2)
    drdx(2,2) = (x(2)-x(8))/r(2)
    drdx(2,3) = (x(3)-x(9))/r(2)
    drdx(2,7) = -drdx(2,1)
    drdx(2,8) = -drdx(2,2)
    drdx(2,9) = -drdx(2,3)

    ! dr3dx
    drdx(3,1)  = (x(1)-x(10))/r(3)
    drdx(3,2)  = (x(2)-x(11))/r(3)
    drdx(3,3)  = (x(3)-x(12))/r(3)
    drdx(3,10) = -drdx(3,1)
    drdx(3,11) = -drdx(3,2)
    drdx(3,12) = -drdx(3,3)

    ! dr4dx
    drdx(4,4) = (x(4)-x(7))/r(4)
    drdx(4,5) = (x(5)-x(8))/r(4)
    drdx(4,6) = (x(6)-x(9))/r(4)
    drdx(4,7) = -drdx(4,4)
    drdx(4,8) = -drdx(4,5)
    drdx(4,9) = -drdx(4,6)

    ! dr5dx
    drdx(5,4)  = (x(4)-x(10))/r(5)
    drdx(5,5)  = (x(5)-x(11))/r(5)
    drdx(5,6)  = (x(6)-x(12))/r(5)
    drdx(5,10) = -drdx(5,4)
    drdx(5,11) = -drdx(5,5)
    drdx(5,12) = -drdx(5,6)

    ! dr6dx
    drdx(6,7)  = (x(7)-x(10))/r(6)
    drdx(6,8)  = (x(8)-x(11))/r(6)
    drdx(6,9)  = (x(9)-x(12))/r(6)
    drdx(6,10) = -drdx(6,7)
    drdx(6,11) = -drdx(6,8)
    drdx(6,12) = -drdx(6,9)
    ! Finish the calculation of non-zero dRdX


end subroutine evdrdx
!> @brief Evaluate EvdMsdR for the N4_1A_PIP_2013 potential.
subroutine evdmsdr
    !**********************************************************************
    ! Subroutine to evaluate the derivatives of Morse term X
    ! w.r.t. interatomic distance R(6)
    ! dmsdR:        Local variables, dirm(6,6)
    ! a:            Nonlinear pamameter(Angstrom)
    ! re:           equilibrium bond length(Angstrom)
    !**********************************************************************

    implicit none

    double precision :: a, dbdr, dmdr, dmsdr, dpdr, drdx, dvdr, r, re

    integer :: i,j

    common /coord/    r(6)
    common /gradt/    dmsdr(6,6),dmdr(6,0:111),dpdr(6,0:305),dvdr(6), &
        & drdx(6,12),dbdr(6,276)
    common /msprmt/   a,re

    ! Initialize dmsdr
    do i=1,6
        do j=1,6
            dmsdr(i,j) = 0.0d0
        enddo
    enddo
    !
    ! Morse term dmsdr = exp(-(r-re)/a)
    ! dmsdr(i,j)=0  i!=j
    ! dmsdr(i,i)= -(1/a)*Exp(-r(i)-re)/a)
    !
    do i=1,6
        dmsdr(i,i) = -(1/a)*exp(-(r(i)-re)/a)
    enddo


end subroutine evdmsdr
!> @brief Evaluate EvdMdR for the N4_1A_PIP_2013 potential.
subroutine evdmdr
    !**********************************************************************
    !  The subroutine reads M(nom) and dMSdR(6,6) and calculates the
    !  dMdR(6,nom) that do not have usable decomposition.
    !  For A4 with max. degree 9, the number of monomials is nom.
    !**********************************************************************

    implicit none

    double precision :: b, c, dbdr, dmdr, dmsdr, dpdr, drdx, dvdr, p, rm, rms

    integer :: i

    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)
    common /gradt/    dmsdr(6,6),dmdr(6,0:111),dpdr(6,0:305),dvdr(6), &
        & drdx(6,12),dbdr(6,276)

    do i=1,6
        dmdr(i,0)   = 0.0d0
        dmdr(i,1)   = dmsdr(i,6)
        dmdr(i,2)   = dmsdr(i,5)
        dmdr(i,3)   = dmsdr(i,4)
        dmdr(i,4)   = dmsdr(i,3)
        dmdr(i,5)   = dmsdr(i,2)
        dmdr(i,6)   = dmsdr(i,1)
        dmdr(i,7)   = dmdr(i,3)*rm(4) + rm(3)*dmdr(i,4)
        dmdr(i,8)   = dmdr(i,2)*rm(5) + rm(2)*dmdr(i,5)
        dmdr(i,9)   = dmdr(i,1)*rm(6) + rm(1)*dmdr(i,6)
        dmdr(i,10)  = dmdr(i,1)*rm(2) + rm(1)*dmdr(i,2)
        dmdr(i,11)  = dmdr(i,1)*rm(3) + rm(1)*dmdr(i,3)
        dmdr(i,12)  = dmdr(i,2)*rm(3) + rm(2)*dmdr(i,3)
        dmdr(i,13)  = dmdr(i,1)*rm(4) + rm(1)*dmdr(i,4)
        dmdr(i,14)  = dmdr(i,2)*rm(4) + rm(2)*dmdr(i,4)
        dmdr(i,15)  = dmdr(i,1)*rm(5) + rm(1)*dmdr(i,5)
        dmdr(i,16)  = dmdr(i,3)*rm(5) + rm(3)*dmdr(i,5)
        dmdr(i,17)  = dmdr(i,4)*rm(5) + rm(4)*dmdr(i,5)
        dmdr(i,18)  = dmdr(i,2)*rm(6) + rm(2)*dmdr(i,6)
        dmdr(i,19)  = dmdr(i,3)*rm(6) + rm(3)*dmdr(i,6)
        dmdr(i,20)  = dmdr(i,4)*rm(6) + rm(4)*dmdr(i,6)
        dmdr(i,21)  = dmdr(i,5)*rm(6) + rm(5)*dmdr(i,6)
        dmdr(i,22)  = dmdr(i,1)*rm(7) + rm(1)*dmdr(i,7)
        dmdr(i,23)  = dmdr(i,2)*rm(7) + rm(2)*dmdr(i,7)
        dmdr(i,24)  = dmdr(i,1)*rm(8) + rm(1)*dmdr(i,8)
        dmdr(i,25)  = dmdr(i,2)*rm(16) + rm(2)*dmdr(i,16)
        dmdr(i,26)  = dmdr(i,2)*rm(17) + rm(2)*dmdr(i,17)
        dmdr(i,27)  = dmdr(i,3)*rm(17) + rm(3)*dmdr(i,17)
        dmdr(i,28)  = dmdr(i,1)*rm(18) + rm(1)*dmdr(i,18)
        dmdr(i,29)  = dmdr(i,1)*rm(19) + rm(1)*dmdr(i,19)
        dmdr(i,30)  = dmdr(i,1)*rm(20) + rm(1)*dmdr(i,20)
        dmdr(i,31)  = dmdr(i,3)*rm(20) + rm(3)*dmdr(i,20)
        dmdr(i,32)  = dmdr(i,1)*rm(21) + rm(1)*dmdr(i,21)
        dmdr(i,33)  = dmdr(i,2)*rm(21) + rm(2)*dmdr(i,21)
        dmdr(i,34)  = dmdr(i,1)*rm(12) + rm(1)*dmdr(i,12)
        dmdr(i,35)  = dmdr(i,1)*rm(17) + rm(1)*dmdr(i,17)
        dmdr(i,36)  = dmdr(i,2)*rm(20) + rm(2)*dmdr(i,20)
        dmdr(i,37)  = dmdr(i,3)*rm(21) + rm(3)*dmdr(i,21)
        dmdr(i,38)  = dmdr(i,1)*rm(14) + rm(1)*dmdr(i,14)
        dmdr(i,39)  = dmdr(i,1)*rm(16) + rm(1)*dmdr(i,16)
        dmdr(i,40)  = dmdr(i,2)*rm(19) + rm(2)*dmdr(i,19)
        dmdr(i,41)  = dmdr(i,4)*rm(21) + rm(4)*dmdr(i,21)
        dmdr(i,42)  = dmdr(i,2)*rm(27) + rm(2)*dmdr(i,27)
        dmdr(i,43)  = dmdr(i,1)*rm(31) + rm(1)*dmdr(i,31)
        dmdr(i,44)  = dmdr(i,1)*rm(33) + rm(1)*dmdr(i,33)
        dmdr(i,45)  = dmdr(i,1)*rm(23) + rm(1)*dmdr(i,23)
        dmdr(i,46)  = dmdr(i,1)*rm(25) + rm(1)*dmdr(i,25)
        dmdr(i,47)  = dmdr(i,1)*rm(26) + rm(1)*dmdr(i,26)
        dmdr(i,48)  = dmdr(i,1)*rm(27) + rm(1)*dmdr(i,27)
        dmdr(i,49)  = dmdr(i,1)*rm(40) + rm(1)*dmdr(i,40)
        dmdr(i,50)  = dmdr(i,1)*rm(36) + rm(1)*dmdr(i,36)
        dmdr(i,51)  = dmdr(i,2)*rm(31) + rm(2)*dmdr(i,31)
        dmdr(i,52)  = dmdr(i,1)*rm(37) + rm(1)*dmdr(i,37)
        dmdr(i,53)  = dmdr(i,2)*rm(37) + rm(2)*dmdr(i,37)
        dmdr(i,54)  = dmdr(i,1)*rm(41) + rm(1)*dmdr(i,41)
        dmdr(i,55)  = dmdr(i,2)*rm(41) + rm(2)*dmdr(i,41)
        dmdr(i,56)  = dmdr(i,3)*rm(41) + rm(3)*dmdr(i,41)
        dmdr(i,57)  = dmdr(i,1)*rm(42) + rm(1)*dmdr(i,42)
        dmdr(i,58)  = dmdr(i,1)*rm(51) + rm(1)*dmdr(i,51)
        dmdr(i,59)  = dmdr(i,1)*rm(53) + rm(1)*dmdr(i,53)
        dmdr(i,60)  = dmdr(i,1)*rm(55) + rm(1)*dmdr(i,55)
        dmdr(i,61)  = dmdr(i,1)*rm(56) + rm(1)*dmdr(i,56)
        dmdr(i,62)  = dmdr(i,2)*rm(56) + rm(2)*dmdr(i,56)
        dmdr(i,63)  = dmdr(i,1)*rm(62) + rm(1)*dmdr(i,62)
        dmdr(i,64)  = dmdr(i,2)*rm(57) + rm(2)*dmdr(i,57)
        dmdr(i,65)  = dmdr(i,3)*rm(57) + rm(3)*dmdr(i,57)
        dmdr(i,66)  = dmdr(i,4)*rm(57) + rm(4)*dmdr(i,57)
        dmdr(i,67)  = dmdr(i,5)*rm(57) + rm(5)*dmdr(i,57)
        dmdr(i,68)  = dmdr(i,1)*rm(58) + rm(1)*dmdr(i,58)
        dmdr(i,69)  = dmdr(i,3)*rm(58) + rm(3)*dmdr(i,58)
        dmdr(i,70)  = dmdr(i,4)*rm(58) + rm(4)*dmdr(i,58)
        dmdr(i,71)  = dmdr(i,1)*rm(59) + rm(1)*dmdr(i,59)
        dmdr(i,72)  = dmdr(i,2)*rm(59) + rm(2)*dmdr(i,59)
        dmdr(i,73)  = dmdr(i,1)*rm(60) + rm(1)*dmdr(i,60)
        dmdr(i,74)  = dmdr(i,2)*rm(60) + rm(2)*dmdr(i,60)
        dmdr(i,75)  = dmdr(i,1)*rm(61) + rm(1)*dmdr(i,61)
        dmdr(i,76)  = dmdr(i,2)*rm(62) + rm(2)*dmdr(i,62)
        dmdr(i,77)  = dmdr(i,3)*rm(61) + rm(3)*dmdr(i,61)
        dmdr(i,78)  = dmdr(i,3)*rm(62) + rm(3)*dmdr(i,62)
        dmdr(i,79)  = dmdr(i,4)*rm(61) + rm(4)*dmdr(i,61)
        dmdr(i,80)  = dmdr(i,4)*rm(62) + rm(4)*dmdr(i,62)
        dmdr(i,81)  = dmdr(i,5)*rm(59) + rm(5)*dmdr(i,59)
        dmdr(i,82)  = dmdr(i,5)*rm(60) + rm(5)*dmdr(i,60)
        dmdr(i,83)  = dmdr(i,5)*rm(62) + rm(5)*dmdr(i,62)
        dmdr(i,84)  = dmdr(i,6)*rm(58) + rm(6)*dmdr(i,58)
        dmdr(i,85)  = dmdr(i,6)*rm(59) + rm(6)*dmdr(i,59)
        dmdr(i,86)  = dmdr(i,6)*rm(60) + rm(6)*dmdr(i,60)
        dmdr(i,87)  = dmdr(i,6)*rm(61) + rm(6)*dmdr(i,61)
        dmdr(i,88)  = dmdr(i,2)*rm(64) + rm(2)*dmdr(i,64)
        dmdr(i,89)  = dmdr(i,3)*rm(65) + rm(3)*dmdr(i,65)
        dmdr(i,90)  = dmdr(i,4)*rm(66) + rm(4)*dmdr(i,66)
        dmdr(i,91)  = dmdr(i,5)*rm(67) + rm(5)*dmdr(i,67)
        dmdr(i,92)  = dmdr(i,1)*rm(68) + rm(1)*dmdr(i,68)
        dmdr(i,93)  = dmdr(i,3)*rm(69) + rm(3)*dmdr(i,69)
        dmdr(i,94)  = dmdr(i,4)*rm(70) + rm(4)*dmdr(i,70)
        dmdr(i,95)  = dmdr(i,1)*rm(71) + rm(1)*dmdr(i,71)
        dmdr(i,96)  = dmdr(i,2)*rm(72) + rm(2)*dmdr(i,72)
        dmdr(i,97)  = dmdr(i,1)*rm(73) + rm(1)*dmdr(i,73)
        dmdr(i,98)  = dmdr(i,2)*rm(74) + rm(2)*dmdr(i,74)
        dmdr(i,99)  = dmdr(i,1)*rm(75) + rm(1)*dmdr(i,75)
        dmdr(i,100) = dmdr(i,2)*rm(76) + rm(2)*dmdr(i,76)
        dmdr(i,101) = dmdr(i,3)*rm(77) + rm(3)*dmdr(i,77)
        dmdr(i,102) = dmdr(i,3)*rm(78) + rm(3)*dmdr(i,78)
        dmdr(i,103) = dmdr(i,4)*rm(79) + rm(4)*dmdr(i,79)
        dmdr(i,104) = dmdr(i,4)*rm(80) + rm(4)*dmdr(i,80)
        dmdr(i,105) = dmdr(i,5)*rm(81) + rm(5)*dmdr(i,81)
        dmdr(i,106) = dmdr(i,5)*rm(82) + rm(5)*dmdr(i,82)
        dmdr(i,107) = dmdr(i,5)*rm(83) + rm(5)*dmdr(i,83)
        dmdr(i,108) = dmdr(i,6)*rm(84) + rm(6)*dmdr(i,84)
        dmdr(i,109) = dmdr(i,6)*rm(85) + rm(6)*dmdr(i,85)
        dmdr(i,110) = dmdr(i,6)*rm(86) + rm(6)*dmdr(i,86)
        dmdr(i,111) = dmdr(i,6)*rm(87) + rm(6)*dmdr(i,87)
    enddo


end subroutine evdmdr
!> @brief Evaluate EvdPdr for the N4_1A_PIP_2013 potential.
subroutine evdpdr
    !**********************************************************************
    !  The subroutine reads monomials(m) and calculates the
    !  permutation invariant polynomials(p)
    !  For A4 with max. degree 9, the number of polynomials is nob.
    !**********************************************************************

    implicit none

    double precision :: b, c, dbdr, dmdr, dmsdr, dpdr, drdx, dvdr, p, rm, rms

    integer :: i

    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)
    common /gradt/    dmsdr(6,6),dmdr(6,0:111),dpdr(6,0:305),dvdr(6), &
        & drdx(6,12),dbdr(6,276)

    do i=1,6
        dpdr(i,0) = dmdr(i,0)
        dpdr(i,1) = dmdr(i,1) + dmdr(i,2) + dmdr(i,3) + dmdr(i,4) + dmdr(i,5) + dmdr(i,6)
        dpdr(i,2) = dmdr(i,7) + dmdr(i,8) + dmdr(i,9)
        dpdr(i,3) = dmdr(i,10) + dmdr(i,11) + dmdr(i,12) + dmdr(i,13) + dmdr(i,14) + dmdr(i,15) + dmdr(i,16) + dmdr(i, &
            & 17) + dmdr(i,18) + dmdr(i,19) + dmdr(i,20) + dmdr(i,21)
        dpdr(i,4) = dpdr(i,1)*p(1) + p(1)*dpdr(i,1) - dpdr(i,3) - dpdr(i,2) - dpdr(i,3) - dpdr(i,2)
        dpdr(i,5) = dmdr(i,22) + dmdr(i,23) + dmdr(i,24) + dmdr(i,25) + dmdr(i,26) + dmdr(i,27) + dmdr(i,28) + dmdr(i, &
            & 29) + dmdr(i,30) + dmdr(i,31) + dmdr(i,32) + dmdr(i,33)
        dpdr(i,6) = dmdr(i,34) + dmdr(i,35) + dmdr(i,36) + dmdr(i,37)
        dpdr(i,7) = dmdr(i,38) + dmdr(i,39) + dmdr(i,40) + dmdr(i,41)
        dpdr(i,8) = dpdr(i,1)*p(2) + p(1)*dpdr(i,2) - dpdr(i,5)
        dpdr(i,9) = dpdr(i,1)*p(3) + p(1)*dpdr(i,3) - dpdr(i,6) - dpdr(i,7) - dpdr(i,5) - dpdr(i,6) - dpdr(i,7) - &
            & dpdr(i,5) - dpdr(i,6) - dpdr(i,7)
        dpdr(i,10) = dpdr(i,1)*p(4) + p(1)*dpdr(i,4) - dpdr(i,9) - dpdr(i,8)
        dpdr(i,11) = dmdr(i,42) + dmdr(i,43) + dmdr(i,44)
        dpdr(i,12) = dmdr(i,45) + dmdr(i,46) + dmdr(i,47) + dmdr(i,48) + dmdr(i,49) + dmdr(i,50) + dmdr(i,51) + &
            & dmdr(i, 52) + dmdr(i,53) + dmdr(i,54) + dmdr(i,55) + dmdr(i,56)
        dpdr(i,13) = dpdr(i,2)*p(3) + p(2)*dpdr(i,3) - dpdr(i,12)
        dpdr(i,14) = dpdr(i,1)*p(5) + p(1)*dpdr(i,5) - dpdr(i,12) - dpdr(i,11) - dpdr(i,13) - dpdr(i,12) - dpdr(i, &
            & 11) - dpdr(i,11) - dpdr(i,11)
        dpdr(i,15) = dpdr(i,1)*p(6) + p(1)*dpdr(i,6) - dpdr(i,12)
        dpdr(i,16) = dpdr(i,1)*p(7) + p(1)*dpdr(i,7) - dpdr(i,12)
        dpdr(i,17) = dpdr(i,2)*p(2) + p(2)*dpdr(i,2) - dpdr(i,11) - dpdr(i,11)
        dpdr(i,18) = dpdr(i,3)*p(3) + p(3)*dpdr(i,3) - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) - dpdr(i,16) - dpdr(i, &
            & 14) - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) - dpdr(i,16) - dpdr(i,14) - dpdr(i,12) - dpdr(i,11) - dpdr(i, &
            & 12) - dpdr(i,11)
        dpdr(i,19) = dpdr(i,2)*p(4) + p(2)*dpdr(i,4) - dpdr(i,14)
        dpdr(i,20) = dpdr(i,3)*p(4) + p(3)*dpdr(i,4) - dpdr(i,15) - dpdr(i,16) - dpdr(i,13)
        dpdr(i,21) = dpdr(i,1)*p(10) + p(1)*dpdr(i,10) - dpdr(i,20) - dpdr(i,19)
        dpdr(i,22) = dmdr(i,57) + dmdr(i,58) + dmdr(i,59) + dmdr(i,60) + dmdr(i,61) + dmdr(i,62)
        dpdr(i,23) = dpdr(i,1)*p(11) + p(1)*dpdr(i,11) - dpdr(i,22)
        dpdr(i,24) = dpdr(i,2)*p(6) + p(2)*dpdr(i,6)
        dpdr(i,25) = dpdr(i,2)*p(7) + p(2)*dpdr(i,7)
        dpdr(i,26) = dpdr(i,1)*p(12) + p(1)*dpdr(i,12) - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) - dpdr(i,22) - dpdr(i, &
            & 22) - dpdr(i,22)
        dpdr(i,27) = dpdr(i,2)*p(5) + p(2)*dpdr(i,5) - dpdr(i,22) - dpdr(i,23) - dpdr(i,22)
        dpdr(i,28) = dpdr(i,3)*p(5) + p(3)*dpdr(i,5) - dpdr(i,22) - dpdr(i,26) - dpdr(i,24) - dpdr(i,25) - dpdr(i, &
            & 23) - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) - dpdr(i,23) - dpdr(i,22) - dpdr(i,22)
        dpdr(i,29) = dpdr(i,3)*p(6) + p(3)*dpdr(i,6) - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
        dpdr(i,30) = dpdr(i,3)*p(7) + p(3)*dpdr(i,7) - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
        dpdr(i,31) = dpdr(i,2)*p(9) + p(2)*dpdr(i,9) - dpdr(i,26) - dpdr(i,28)
        dpdr(i,32) = dpdr(i,1)*p(14) + p(1)*dpdr(i,14) - dpdr(i,26) - dpdr(i,23) - dpdr(i,28)
        dpdr(i,33) = dpdr(i,4)*p(6) + p(4)*dpdr(i,6) - dpdr(i,25)
        dpdr(i,34) = dpdr(i,4)*p(7) + p(4)*dpdr(i,7) - dpdr(i,24)
        dpdr(i,35) = dpdr(i,1)*p(17) + p(1)*dpdr(i,17) - dpdr(i,27)
        dpdr(i,36) = dpdr(i,1)*p(18) + p(1)*dpdr(i,18) - dpdr(i,29) - dpdr(i,30) - dpdr(i,28)
        dpdr(i,37) = dpdr(i,2)*p(10) + p(2)*dpdr(i,10) - dpdr(i,32)
        dpdr(i,38) = dpdr(i,3)*p(10) + p(3)*dpdr(i,10) - dpdr(i,33) - dpdr(i,34) - dpdr(i,31)
        dpdr(i,39) = dpdr(i,1)*p(21) + p(1)*dpdr(i,21) - dpdr(i,38) - dpdr(i,37)
        dpdr(i,40) = dmdr(i,63)
        dpdr(i,41) = dmdr(i,64) + dmdr(i,65) + dmdr(i,66) + dmdr(i,67) + dmdr(i,68) + dmdr(i,69) + dmdr(i,70) + &
            & dmdr(i, 71) + dmdr(i,72) + dmdr(i,73) + dmdr(i,74) + dmdr(i,75) + dmdr(i,76) + dmdr(i,77) + dmdr(i,78) + &
            & dmdr(i, 79) + dmdr(i,80) + dmdr(i,81) + dmdr(i,82) + dmdr(i,83) + dmdr(i,84) + dmdr(i,85) + dmdr(i,86) + &
            & dmdr(i, 87)
        dpdr(i,42) = dpdr(i,1)*p(22) + p(1)*dpdr(i,22) - dpdr(i,40) - dpdr(i,41) - dpdr(i,40) - dpdr(i,40) - dpdr(i, &
            & 40) - dpdr(i,40) - dpdr(i,40)
        dpdr(i,43) = dpdr(i,2)*p(11) + p(2)*dpdr(i,11) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
        dpdr(i,44) = dpdr(i,2)*p(12) + p(2)*dpdr(i,12) - dpdr(i,41)
        dpdr(i,45) = dpdr(i,3)*p(11) + p(3)*dpdr(i,11) - dpdr(i,41)
        dpdr(i,46) = dpdr(i,5)*p(6) + p(5)*dpdr(i,6) - dpdr(i,41)
        dpdr(i,47) = dpdr(i,5)*p(7) + p(5)*dpdr(i,7) - dpdr(i,41)
        dpdr(i,48) = dpdr(i,6)*p(7) + p(6)*dpdr(i,7) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
        dpdr(i,49) = dpdr(i,4)*p(11) + p(4)*dpdr(i,11) - dpdr(i,42)
        dpdr(i,50) = dpdr(i,2)*p(15) + p(2)*dpdr(i,15) - dpdr(i,46)
        dpdr(i,51) = dpdr(i,2)*p(16) + p(2)*dpdr(i,16) - dpdr(i,47)
        dpdr(i,52) = dpdr(i,4)*p(12) + p(4)*dpdr(i,12) - dpdr(i,41) - dpdr(i,50) - dpdr(i,51)
        dpdr(i,53) = dpdr(i,2)*p(14) + p(2)*dpdr(i,14) - dpdr(i,42) - dpdr(i,49) - dpdr(i,42)
        dpdr(i,54) = dpdr(i,6)*p(6) + p(6)*dpdr(i,6) - dpdr(i,42) - dpdr(i,42)
        dpdr(i,55) = dpdr(i,7)*p(7) + p(7)*dpdr(i,7) - dpdr(i,42) - dpdr(i,42)
        dpdr(i,56) = dpdr(i,3)*p(17) + p(3)*dpdr(i,17) - dpdr(i,44)
        dpdr(i,57) = dpdr(i,2)*p(18) + p(2)*dpdr(i,18) - dpdr(i,48)
        dpdr(i,58) = dpdr(i,3)*p(14) + p(3)*dpdr(i,14) - dpdr(i,41) - dpdr(i,52) - dpdr(i,46) - dpdr(i,47) - dpdr(i, &
            & 45) - dpdr(i,45)
        dpdr(i,59) = dpdr(i,6)*p(9) + p(6)*dpdr(i,9) - dpdr(i,41) - dpdr(i,52) - dpdr(i,47)
        dpdr(i,60) = dpdr(i,7)*p(9) + p(7)*dpdr(i,9) - dpdr(i,41) - dpdr(i,52) - dpdr(i,46)
        dpdr(i,61) = dpdr(i,2)*p(20) + p(2)*dpdr(i,20) - dpdr(i,52) - dpdr(i,58)
        dpdr(i,62) = dpdr(i,1)*p(32) + p(1)*dpdr(i,32) - dpdr(i,52) - dpdr(i,49) - dpdr(i,58)
        dpdr(i,63) = dpdr(i,6)*p(10) + p(6)*dpdr(i,10) - dpdr(i,51)
        dpdr(i,64) = dpdr(i,7)*p(10) + p(7)*dpdr(i,10) - dpdr(i,50)
        dpdr(i,65) = dpdr(i,2)*p(17) + p(2)*dpdr(i,17) - dpdr(i,43)
        dpdr(i,66) = dpdr(i,3)*p(18) + p(3)*dpdr(i,18) - dpdr(i,46) - dpdr(i,47) - dpdr(i,45) - dpdr(i,59) - dpdr(i, &
            & 60) - dpdr(i,58)
        dpdr(i,67) = dpdr(i,2)*p(19) + p(2)*dpdr(i,19) - dpdr(i,49)
        dpdr(i,68) = dpdr(i,1)*p(36) + p(1)*dpdr(i,36) - dpdr(i,59) - dpdr(i,60) - dpdr(i,58) - dpdr(i,57) - dpdr(i, &
            & 66) - dpdr(i,66)
        dpdr(i,69) = dpdr(i,2)*p(21) + p(2)*dpdr(i,21) - dpdr(i,62)
        dpdr(i,70) = dpdr(i,3)*p(21) + p(3)*dpdr(i,21) - dpdr(i,63) - dpdr(i,64) - dpdr(i,61)
        dpdr(i,71) = dpdr(i,1)*p(39) + p(1)*dpdr(i,39) - dpdr(i,70) - dpdr(i,69)
        dpdr(i,72) = dpdr(i,40)*p(1) + p(40)*dpdr(i,1)
        dpdr(i,73) = dpdr(i,2)*p(22) + p(2)*dpdr(i,22) - dpdr(i,72)
        dpdr(i,74) = dpdr(i,6)*p(11) + p(6)*dpdr(i,11)
        dpdr(i,75) = dpdr(i,7)*p(11) + p(7)*dpdr(i,11)
        dpdr(i,76) = dpdr(i,3)*p(22) + p(3)*dpdr(i,22) - dpdr(i,72) - dpdr(i,74) - dpdr(i,75) - dpdr(i,72) - dpdr(i, &
            & 72) - dpdr(i,72)
        dpdr(i,77) = dmdr(i,88) + dmdr(i,89) + dmdr(i,90) + dmdr(i,91) + dmdr(i,92) + dmdr(i,93) + dmdr(i,94) + &
            & dmdr(i, 95) + dmdr(i,96) + dmdr(i,97) + dmdr(i,98) + dmdr(i,99) + dmdr(i,100) + dmdr(i,101) + dmdr(i, &
            & 102) + dmdr(i,103) + dmdr(i,104) + dmdr(i,105) + dmdr(i,106) + dmdr(i,107) + dmdr(i,108) + dmdr(i,109) + &
            & dmdr(i, 110) + dmdr(i,111)
        dpdr(i,78)  = dpdr(i,1)*p(42) + p(1)*dpdr(i,42) - dpdr(i,72) - dpdr(i,76)
        dpdr(i,79)  = dpdr(i,5)*p(11) + p(5)*dpdr(i,11) - dpdr(i,72) - dpdr(i,73) - dpdr(i,72)
        dpdr(i,80)  = dpdr(i,2)*p(26) + p(2)*dpdr(i,26) - dpdr(i,76) - dpdr(i,77)
        dpdr(i,81)  = dpdr(i,6)*p(12) + p(6)*dpdr(i,12) - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
        dpdr(i,82)  = dpdr(i,7)*p(12) + p(7)*dpdr(i,12) - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
        dpdr(i,83)  = dpdr(i,8)*p(11) + p(8)*dpdr(i,11) - dpdr(i,72)
        dpdr(i,84)  = dpdr(i,6)*p(17) + p(6)*dpdr(i,17)
        dpdr(i,85)  = dpdr(i,7)*p(17) + p(7)*dpdr(i,17)
        dpdr(i,86)  = dpdr(i,9)*p(11) + p(9)*dpdr(i,11) - dpdr(i,76) - dpdr(i,77)
        dpdr(i,87)  = dpdr(i,2)*p(29) + p(2)*dpdr(i,29) - dpdr(i,81)
        dpdr(i,88)  = dpdr(i,6)*p(14) + p(6)*dpdr(i,14) - dpdr(i,75) - dpdr(i,75)
        dpdr(i,89)  = dpdr(i,2)*p(30) + p(2)*dpdr(i,30) - dpdr(i,82)
        dpdr(i,90)  = dpdr(i,7)*p(14) + p(7)*dpdr(i,14) - dpdr(i,74) - dpdr(i,74)
        dpdr(i,91)  = dpdr(i,1)*p(48) + p(1)*dpdr(i,48) - dpdr(i,76) - dpdr(i,81) - dpdr(i,82)
        dpdr(i,92)  = dpdr(i,10)*p(11) + p(10)*dpdr(i,11) - dpdr(i,78)
        dpdr(i,93)  = dpdr(i,2)*p(33) + p(2)*dpdr(i,33) - dpdr(i,88)
        dpdr(i,94)  = dpdr(i,2)*p(34) + p(2)*dpdr(i,34) - dpdr(i,90)
        dpdr(i,95)  = dpdr(i,10)*p(12) + p(10)*dpdr(i,12) - dpdr(i,77) - dpdr(i,93) - dpdr(i,94)
        dpdr(i,96)  = dpdr(i,2)*p(27) + p(2)*dpdr(i,27) - dpdr(i,73) - dpdr(i,79)
        dpdr(i,97)  = dpdr(i,2)*p(28) + p(2)*dpdr(i,28) - dpdr(i,76) - dpdr(i,86)
        dpdr(i,98)  = dpdr(i,1)*p(53) + p(1)*dpdr(i,53) - dpdr(i,80) - dpdr(i,79) - dpdr(i,97)
        dpdr(i,99)  = dpdr(i,1)*p(54) + p(1)*dpdr(i,54) - dpdr(i,81)
        dpdr(i,100) = dpdr(i,1)*p(55) + p(1)*dpdr(i,55) - dpdr(i,82)
        dpdr(i,101) = dpdr(i,5)*p(18) + p(5)*dpdr(i,18) - dpdr(i,76) - dpdr(i,91) - dpdr(i,87) - dpdr(i,89) - dpdr(i,86)
        dpdr(i,102) = dpdr(i,6)*p(18) + p(6)*dpdr(i,18) - dpdr(i,74) - dpdr(i,90)
        dpdr(i,103) = dpdr(i,7)*p(18) + p(7)*dpdr(i,18) - dpdr(i,75) - dpdr(i,88)
        dpdr(i,104) = dpdr(i,2)*p(31) + p(2)*dpdr(i,31) - dpdr(i,77) - dpdr(i,86)
        dpdr(i,105) = dpdr(i,2)*p(36) + p(2)*dpdr(i,36) - dpdr(i,91) - dpdr(i,101)
        dpdr(i,106) = dpdr(i,3)*p(32) + p(3)*dpdr(i,32) - dpdr(i,77) - dpdr(i,95) - dpdr(i,88) - dpdr(i,90) - dpdr(i,86)
        dpdr(i,107) = dpdr(i,4)*p(29) + p(4)*dpdr(i,29) - dpdr(i,82) - dpdr(i,80) - dpdr(i,99)
        dpdr(i,108) = dpdr(i,4)*p(30) + p(4)*dpdr(i,30) - dpdr(i,81) - dpdr(i,80) - dpdr(i,100)
        dpdr(i,109) = dpdr(i,2)*p(38) + p(2)*dpdr(i,38) - dpdr(i,95) - dpdr(i,106)
        dpdr(i,110) = dpdr(i,1)*p(62) + p(1)*dpdr(i,62) - dpdr(i,95) - dpdr(i,92) - dpdr(i,106)
        dpdr(i,111) = dpdr(i,6)*p(21) + p(6)*dpdr(i,21) - dpdr(i,94)
        dpdr(i,112) = dpdr(i,7)*p(21) + p(7)*dpdr(i,21) - dpdr(i,93)
        dpdr(i,113) = dpdr(i,1)*p(65) + p(1)*dpdr(i,65) - dpdr(i,96)
        dpdr(i,114) = dpdr(i,1)*p(66) + p(1)*dpdr(i,66) - dpdr(i,102) - dpdr(i,103) - dpdr(i,101)
        dpdr(i,115) = dpdr(i,2)*p(37) + p(2)*dpdr(i,37) - dpdr(i,92)
        dpdr(i,116) = dpdr(i,10)*p(18) + p(10)*dpdr(i,18) - dpdr(i,99) - dpdr(i,100) - dpdr(i,97)
        dpdr(i,117) = dpdr(i,2)*p(39) + p(2)*dpdr(i,39) - dpdr(i,110)
        dpdr(i,118) = dpdr(i,3)*p(39) + p(3)*dpdr(i,39) - dpdr(i,111) - dpdr(i,112) - dpdr(i,109)
        dpdr(i,119) = dpdr(i,1)*p(71) + p(1)*dpdr(i,71) - dpdr(i,118) - dpdr(i,117)
        dpdr(i,120) = dpdr(i,40)*p(2) + p(40)*dpdr(i,2)
        dpdr(i,121) = dpdr(i,40)*p(3) + p(40)*dpdr(i,3)
        dpdr(i,122) = dpdr(i,40)*p(4) + p(40)*dpdr(i,4)
        dpdr(i,123) = dpdr(i,11)*p(12) + p(11)*dpdr(i,12) - dpdr(i,121)
        dpdr(i,124) = dpdr(i,2)*p(42) + p(2)*dpdr(i,42) - dpdr(i,122)
        dpdr(i,125) = dpdr(i,6)*p(22) + p(6)*dpdr(i,22) - dpdr(i,121)
        dpdr(i,126) = dpdr(i,7)*p(22) + p(7)*dpdr(i,22) - dpdr(i,121)
        dpdr(i,127) = dpdr(i,2)*p(41) + p(2)*dpdr(i,41) - dpdr(i,121) - dpdr(i,123) - dpdr(i,121)
        dpdr(i,128) = dpdr(i,6)*p(23) + p(6)*dpdr(i,23) - dpdr(i,123)
        dpdr(i,129) = dpdr(i,7)*p(23) + p(7)*dpdr(i,23) - dpdr(i,123)
        dpdr(i,130) = dpdr(i,3)*p(41) + p(3)*dpdr(i,41) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,125) - &
            & dpdr(i,128) - dpdr(i,126) - dpdr(i,129) - dpdr(i,124) - dpdr(i,123) - dpdr(i,122) - dpdr(i,121) - &
            & dpdr(i, 120) - dpdr(i,125) - dpdr(i,126) - dpdr(i,124) - dpdr(i,123) - dpdr(i,122) - dpdr(i,121) - &
            & dpdr(i,120) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120) - &
            & dpdr(i,120)
        dpdr(i,131) = dpdr(i,3)*p(42) + p(3)*dpdr(i,42) - dpdr(i,121) - dpdr(i,125) - dpdr(i,126) - dpdr(i,121)
        dpdr(i,132) = dpdr(i,4)*p(41) + p(4)*dpdr(i,41) - dpdr(i,121) - dpdr(i,131) - dpdr(i,128) - dpdr(i,129) - &
            & dpdr(i,127) - dpdr(i,121)
        dpdr(i,133) = dpdr(i,1)*p(78) + p(1)*dpdr(i,78) - dpdr(i,122) - dpdr(i,131)
        dpdr(i,134) = dpdr(i,11)*p(11) + p(11)*dpdr(i,11) - dpdr(i,120) - dpdr(i,120)
        dpdr(i,135) = dpdr(i,2)*p(48) + p(2)*dpdr(i,48) - dpdr(i,130)
        dpdr(i,136) = dpdr(i,11)*p(17) + p(11)*dpdr(i,17) - dpdr(i,120)
        dpdr(i,137) = dpdr(i,2)*p(44) + p(2)*dpdr(i,44) - dpdr(i,123)
        dpdr(i,138) = dpdr(i,2)*p(45) + p(2)*dpdr(i,45) - dpdr(i,121)
        dpdr(i,139) = dpdr(i,11)*p(14) + p(11)*dpdr(i,14) - dpdr(i,122) - dpdr(i,124) - dpdr(i,122)
        dpdr(i,140) = dpdr(i,6)*p(27) + p(6)*dpdr(i,27) - dpdr(i,127)
        dpdr(i,141) = dpdr(i,2)*p(54) + p(2)*dpdr(i,54)
        dpdr(i,142) = dpdr(i,7)*p(27) + p(7)*dpdr(i,27) - dpdr(i,127)
        dpdr(i,143) = dpdr(i,2)*p(52) + p(2)*dpdr(i,52) - dpdr(i,131) - dpdr(i,132)
        dpdr(i,144) = dpdr(i,1)*p(81) + p(1)*dpdr(i,81) - dpdr(i,125) - dpdr(i,141) - dpdr(i,135) - dpdr(i,125)
        dpdr(i,145) = dpdr(i,2)*p(55) + p(2)*dpdr(i,55)
        dpdr(i,146) = dpdr(i,1)*p(82) + p(1)*dpdr(i,82) - dpdr(i,126) - dpdr(i,135) - dpdr(i,145) - dpdr(i,126)
        dpdr(i,147) = dpdr(i,11)*p(18) + p(11)*dpdr(i,18) - dpdr(i,130)
        dpdr(i,148) = dpdr(i,6)*p(28) + p(6)*dpdr(i,28) - dpdr(i,129) - dpdr(i,123) - dpdr(i,143)
        dpdr(i,149) = dpdr(i,7)*p(28) + p(7)*dpdr(i,28) - dpdr(i,128) - dpdr(i,123) - dpdr(i,143)
        dpdr(i,150) = dpdr(i,6)*p(30) + p(6)*dpdr(i,30) - dpdr(i,121) - dpdr(i,146)
        dpdr(i,151) = dpdr(i,11)*p(19) + p(11)*dpdr(i,19) - dpdr(i,122)
        dpdr(i,152) = dpdr(i,2)*p(50) + p(2)*dpdr(i,50) - dpdr(i,128)
        dpdr(i,153) = dpdr(i,2)*p(51) + p(2)*dpdr(i,51) - dpdr(i,129)
        dpdr(i,154) = dpdr(i,11)*p(20) + p(11)*dpdr(i,20) - dpdr(i,131) - dpdr(i,132)
        dpdr(i,155) = dpdr(i,2)*p(59) + p(2)*dpdr(i,59) - dpdr(i,144) - dpdr(i,148)
        dpdr(i,156) = dpdr(i,6)*p(32) + p(6)*dpdr(i,32) - dpdr(i,129)
        dpdr(i,157) = dpdr(i,2)*p(60) + p(2)*dpdr(i,60) - dpdr(i,146) - dpdr(i,149)
        dpdr(i,158) = dpdr(i,7)*p(32) + p(7)*dpdr(i,32) - dpdr(i,128)
        dpdr(i,159) = dpdr(i,6)*p(34) + p(6)*dpdr(i,34) - dpdr(i,122) - dpdr(i,145) - dpdr(i,122)
        dpdr(i,160) = dpdr(i,11)*p(21) + p(11)*dpdr(i,21) - dpdr(i,133)
        dpdr(i,161) = dpdr(i,2)*p(63) + p(2)*dpdr(i,63) - dpdr(i,156)
        dpdr(i,162) = dpdr(i,2)*p(64) + p(2)*dpdr(i,64) - dpdr(i,158)
        dpdr(i,163) = dpdr(i,12)*p(21) + p(12)*dpdr(i,21) - dpdr(i,132) - dpdr(i,161) - dpdr(i,162)
        dpdr(i,164) = dpdr(i,2)*p(53) + p(2)*dpdr(i,53) - dpdr(i,124) - dpdr(i,139)
        dpdr(i,165) = dpdr(i,2)*p(58) + p(2)*dpdr(i,58) - dpdr(i,131) - dpdr(i,154)
        dpdr(i,166) = dpdr(i,3)*p(54) + p(3)*dpdr(i,54) - dpdr(i,125) - dpdr(i,144)
        dpdr(i,167) = dpdr(i,3)*p(55) + p(3)*dpdr(i,55) - dpdr(i,126) - dpdr(i,146)
        dpdr(i,168) = dpdr(i,3)*p(65) + p(3)*dpdr(i,65) - dpdr(i,137)
        dpdr(i,169) = dpdr(i,17)*p(18) + p(17)*dpdr(i,18) - dpdr(i,135)
        dpdr(i,170) = dpdr(i,1)*p(98) + p(1)*dpdr(i,98) - dpdr(i,143) - dpdr(i,139) - dpdr(i,165)
        dpdr(i,171) = dpdr(i,4)*p(54) + p(4)*dpdr(i,54) - dpdr(i,135)
        dpdr(i,172) = dpdr(i,4)*p(55) + p(4)*dpdr(i,55) - dpdr(i,135)
        dpdr(i,173) = dpdr(i,2)*p(66) + p(2)*dpdr(i,66) - dpdr(i,150)
        dpdr(i,174) = dpdr(i,1)*p(101) + p(1)*dpdr(i,101) - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,173) - &
            & dpdr(i,165) - dpdr(i,147)
        dpdr(i,175) = dpdr(i,1)*p(102) + p(1)*dpdr(i,102) - dpdr(i,150) - dpdr(i,148) - dpdr(i,166)
        dpdr(i,176) = dpdr(i,1)*p(103) + p(1)*dpdr(i,103) - dpdr(i,150) - dpdr(i,149) - dpdr(i,167)
        dpdr(i,177) = dpdr(i,2)*p(61) + p(2)*dpdr(i,61) - dpdr(i,132) - dpdr(i,154)
        dpdr(i,178) = dpdr(i,2)*p(68) + p(2)*dpdr(i,68) - dpdr(i,159) - dpdr(i,174)
        dpdr(i,179) = dpdr(i,3)*p(62) + p(3)*dpdr(i,62) - dpdr(i,132) - dpdr(i,163) - dpdr(i,156) - dpdr(i,158) - &
            & dpdr(i,154)
        dpdr(i,180) = dpdr(i,6)*p(38) + p(6)*dpdr(i,38) - dpdr(i,132) - dpdr(i,163) - dpdr(i,157)
        dpdr(i,181) = dpdr(i,7)*p(38) + p(7)*dpdr(i,38) - dpdr(i,132) - dpdr(i,163) - dpdr(i,155)
        dpdr(i,182) = dpdr(i,2)*p(70) + p(2)*dpdr(i,70) - dpdr(i,163) - dpdr(i,179)
        dpdr(i,183) = dpdr(i,1)*p(110) + p(1)*dpdr(i,110) - dpdr(i,163) - dpdr(i,160) - dpdr(i,179)
        dpdr(i,184) = dpdr(i,6)*p(39) + p(6)*dpdr(i,39) - dpdr(i,162)
        dpdr(i,185) = dpdr(i,7)*p(39) + p(7)*dpdr(i,39) - dpdr(i,161)
        dpdr(i,186) = dpdr(i,2)*p(65) + p(2)*dpdr(i,65) - dpdr(i,136)
        dpdr(i,187) = dpdr(i,3)*p(66) + p(3)*dpdr(i,66) - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,175) - &
            & dpdr(i,176) - dpdr(i,174)
        dpdr(i,188) = dpdr(i,2)*p(67) + p(2)*dpdr(i,67) - dpdr(i,151)
        dpdr(i,189) = dpdr(i,4)*p(66) + p(4)*dpdr(i,66) - dpdr(i,166) - dpdr(i,167) - dpdr(i,165)
        dpdr(i,190) = dpdr(i,2)*p(69) + p(2)*dpdr(i,69) - dpdr(i,160)
        dpdr(i,191) = dpdr(i,18)*p(21) + p(18)*dpdr(i,21) - dpdr(i,171) - dpdr(i,172) - dpdr(i,169)
        dpdr(i,192) = dpdr(i,2)*p(71) + p(2)*dpdr(i,71) - dpdr(i,183)
        dpdr(i,193) = dpdr(i,3)*p(71) + p(3)*dpdr(i,71) - dpdr(i,184) - dpdr(i,185) - dpdr(i,182)
        dpdr(i,194) = dpdr(i,1)*p(119) + p(1)*dpdr(i,119) - dpdr(i,193) - dpdr(i,192)
        dpdr(i,195) = dpdr(i,40)*p(5) + p(40)*dpdr(i,5)
        dpdr(i,196) = dpdr(i,40)*p(6) + p(40)*dpdr(i,6)
        dpdr(i,197) = dpdr(i,40)*p(7) + p(40)*dpdr(i,7)
        dpdr(i,198) = dpdr(i,40)*p(8) + p(40)*dpdr(i,8)
        dpdr(i,199) = dpdr(i,40)*p(9) + p(40)*dpdr(i,9)
        dpdr(i,200) = dpdr(i,40)*p(10) + p(40)*dpdr(i,10)
        dpdr(i,201) = dpdr(i,11)*p(22) + p(11)*dpdr(i,22) - dpdr(i,195)
        dpdr(i,202) = dpdr(i,12)*p(22) + p(12)*dpdr(i,22) - dpdr(i,196) - dpdr(i,197) - dpdr(i,195) - dpdr(i,196) - &
            & dpdr(i,197) - dpdr(i,195) - dpdr(i,196) - dpdr(i,197)
        dpdr(i,203) = dpdr(i,17)*p(22) + p(17)*dpdr(i,22) - dpdr(i,198)
        dpdr(i,204) = dpdr(i,6)*p(43) + p(6)*dpdr(i,43)
        dpdr(i,205) = dpdr(i,7)*p(43) + p(7)*dpdr(i,43)
        dpdr(i,206) = dpdr(i,11)*p(26) + p(11)*dpdr(i,26) - dpdr(i,199) - dpdr(i,202)
        dpdr(i,207) = dpdr(i,2)*p(76) + p(2)*dpdr(i,76) - dpdr(i,199) - dpdr(i,202)
        dpdr(i,208) = dpdr(i,2)*p(78) + p(2)*dpdr(i,78) - dpdr(i,200)
        dpdr(i,209) = dpdr(i,6)*p(41) + p(6)*dpdr(i,41) - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
        dpdr(i,210) = dpdr(i,6)*p(42) + p(6)*dpdr(i,42) - dpdr(i,197) - dpdr(i,197) - dpdr(i,197)
        dpdr(i,211) = dpdr(i,7)*p(41) + p(7)*dpdr(i,41) - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
        dpdr(i,212) = dpdr(i,7)*p(42) + p(7)*dpdr(i,42) - dpdr(i,196) - dpdr(i,196) - dpdr(i,196)
        dpdr(i,213) = dpdr(i,11)*p(29) + p(11)*dpdr(i,29) - dpdr(i,209)
        dpdr(i,214) = dpdr(i,11)*p(30) + p(11)*dpdr(i,30) - dpdr(i,211)
        dpdr(i,215) = dpdr(i,18)*p(22) + p(18)*dpdr(i,22) - dpdr(i,199) - dpdr(i,213) - dpdr(i,214)
        dpdr(i,216) = dpdr(i,2)*p(77) + p(2)*dpdr(i,77) - dpdr(i,199) - dpdr(i,206)
        dpdr(i,217) = dpdr(i,6)*p(49) + p(6)*dpdr(i,49) - dpdr(i,205)
        dpdr(i,218) = dpdr(i,7)*p(49) + p(7)*dpdr(i,49) - dpdr(i,204)
        dpdr(i,219) = dpdr(i,3)*p(77) + p(3)*dpdr(i,77) - dpdr(i,200) - dpdr(i,199) - dpdr(i,198) - dpdr(i,209) - &
            & dpdr(i,217) - dpdr(i,211) - dpdr(i,218) - dpdr(i,207) - dpdr(i,204) - dpdr(i,205) - dpdr(i,200) - &
            & dpdr(i, 199) - dpdr(i,198) - dpdr(i,200) - dpdr(i,198) - dpdr(i,200) - dpdr(i,198)
        dpdr(i,220) = dpdr(i,3)*p(78) + p(3)*dpdr(i,78) - dpdr(i,199) - dpdr(i,210) - dpdr(i,212)
        dpdr(i,221) = dpdr(i,10)*p(41) + p(10)*dpdr(i,41) - dpdr(i,199) - dpdr(i,220) - dpdr(i,217) - dpdr(i,218) - &
            & dpdr(i,216)
        dpdr(i,222) = dpdr(i,1)*p(133) + p(1)*dpdr(i,133) - dpdr(i,200) - dpdr(i,220)
        dpdr(i,223) = dpdr(i,11)*p(27) + p(11)*dpdr(i,27) - dpdr(i,195) - dpdr(i,203)
        dpdr(i,224) = dpdr(i,1)*p(134) + p(1)*dpdr(i,134) - dpdr(i,201)
        dpdr(i,225) = dpdr(i,2)*p(81) + p(2)*dpdr(i,81) - dpdr(i,209)
        dpdr(i,226) = dpdr(i,2)*p(80) + p(2)*dpdr(i,80) - dpdr(i,202) - dpdr(i,206)
        dpdr(i,227) = dpdr(i,2)*p(82) + p(2)*dpdr(i,82) - dpdr(i,211)
        dpdr(i,228) = dpdr(i,2)*p(91) + p(2)*dpdr(i,91) - dpdr(i,215) - dpdr(i,219)
        dpdr(i,229) = dpdr(i,11)*p(28) + p(11)*dpdr(i,28) - dpdr(i,199) - dpdr(i,207)
        dpdr(i,230) = dpdr(i,6)*p(53) + p(6)*dpdr(i,53) - dpdr(i,205)
        dpdr(i,231) = dpdr(i,5)*p(54) + p(5)*dpdr(i,54) - dpdr(i,209)
        dpdr(i,232) = dpdr(i,7)*p(53) + p(7)*dpdr(i,53) - dpdr(i,204)
        dpdr(i,233) = dpdr(i,7)*p(54) + p(7)*dpdr(i,54) - dpdr(i,196)
        dpdr(i,234) = dpdr(i,5)*p(55) + p(5)*dpdr(i,55) - dpdr(i,211)
        dpdr(i,235) = dpdr(i,6)*p(55) + p(6)*dpdr(i,55) - dpdr(i,197)
        dpdr(i,236) = dpdr(i,11)*p(35) + p(11)*dpdr(i,35) - dpdr(i,198)
        dpdr(i,237) = dpdr(i,6)*p(65) + p(6)*dpdr(i,65)
        dpdr(i,238) = dpdr(i,7)*p(65) + p(7)*dpdr(i,65)
        dpdr(i,239) = dpdr(i,2)*p(86) + p(2)*dpdr(i,86) - dpdr(i,199) - dpdr(i,229)
        dpdr(i,240) = dpdr(i,1)*p(139) + p(1)*dpdr(i,139) - dpdr(i,206) - dpdr(i,229) - dpdr(i,224)
        dpdr(i,241) = dpdr(i,17)*p(29) + p(17)*dpdr(i,29) - dpdr(i,225)
        dpdr(i,242) = dpdr(i,2)*p(99) + p(2)*dpdr(i,99) - dpdr(i,231)
        dpdr(i,243) = dpdr(i,17)*p(30) + p(17)*dpdr(i,30) - dpdr(i,227)
        dpdr(i,244) = dpdr(i,2)*p(95) + p(2)*dpdr(i,95) - dpdr(i,220) - dpdr(i,221)
        dpdr(i,245) = dpdr(i,4)*p(81) + p(4)*dpdr(i,81) - dpdr(i,202) - dpdr(i,242) - dpdr(i,227)
        dpdr(i,246) = dpdr(i,2)*p(100) + p(2)*dpdr(i,100) - dpdr(i,234)
        dpdr(i,247) = dpdr(i,4)*p(82) + p(4)*dpdr(i,82) - dpdr(i,202) - dpdr(i,225) - dpdr(i,246)
        dpdr(i,248) = dpdr(i,11)*p(36) + p(11)*dpdr(i,36) - dpdr(i,215) - dpdr(i,219)
        dpdr(i,249) = dpdr(i,2)*p(102) + p(2)*dpdr(i,102) - dpdr(i,233)
        dpdr(i,250) = dpdr(i,6)*p(58) + p(6)*dpdr(i,58) - dpdr(i,206) - dpdr(i,214) - dpdr(i,244) - dpdr(i,214)
        dpdr(i,251) = dpdr(i,2)*p(103) + p(2)*dpdr(i,103) - dpdr(i,235)
        dpdr(i,252) = dpdr(i,7)*p(58) + p(7)*dpdr(i,58) - dpdr(i,213) - dpdr(i,206) - dpdr(i,244) - dpdr(i,213)
        dpdr(i,253) = dpdr(i,1)*p(150) + p(1)*dpdr(i,150) - dpdr(i,215) - dpdr(i,233) - dpdr(i,235)
        dpdr(i,254) = dpdr(i,11)*p(37) + p(11)*dpdr(i,37) - dpdr(i,200)
        dpdr(i,255) = dpdr(i,2)*p(93) + p(2)*dpdr(i,93) - dpdr(i,217)
        dpdr(i,256) = dpdr(i,2)*p(94) + p(2)*dpdr(i,94) - dpdr(i,218)
        dpdr(i,257) = dpdr(i,11)*p(38) + p(11)*dpdr(i,38) - dpdr(i,220) - dpdr(i,221)
        dpdr(i,258) = dpdr(i,2)*p(107) + p(2)*dpdr(i,107) - dpdr(i,245) - dpdr(i,250)
        dpdr(i,259) = dpdr(i,6)*p(62) + p(6)*dpdr(i,62) - dpdr(i,218)
        dpdr(i,260) = dpdr(i,2)*p(108) + p(2)*dpdr(i,108) - dpdr(i,247) - dpdr(i,252)
        dpdr(i,261) = dpdr(i,7)*p(62) + p(7)*dpdr(i,62) - dpdr(i,217)
        dpdr(i,262) = dpdr(i,6)*p(64) + p(6)*dpdr(i,64) - dpdr(i,200) - dpdr(i,246) - dpdr(i,200)
        dpdr(i,263) = dpdr(i,11)*p(39) + p(11)*dpdr(i,39) - dpdr(i,222)
        dpdr(i,264) = dpdr(i,2)*p(111) + p(2)*dpdr(i,111) - dpdr(i,259)
        dpdr(i,265) = dpdr(i,2)*p(112) + p(2)*dpdr(i,112) - dpdr(i,261)
        dpdr(i,266) = dpdr(i,12)*p(39) + p(12)*dpdr(i,39) - dpdr(i,221) - dpdr(i,264) - dpdr(i,265)
        dpdr(i,267) = dpdr(i,2)*p(98) + p(2)*dpdr(i,98) - dpdr(i,208) - dpdr(i,240)
        dpdr(i,268) = dpdr(i,6)*p(54) + p(6)*dpdr(i,54) - dpdr(i,210)
        dpdr(i,269) = dpdr(i,7)*p(55) + p(7)*dpdr(i,55) - dpdr(i,212)
        dpdr(i,270) = dpdr(i,2)*p(96) + p(2)*dpdr(i,96) - dpdr(i,203) - dpdr(i,223)
        dpdr(i,271) = dpdr(i,2)*p(97) + p(2)*dpdr(i,97) - dpdr(i,207) - dpdr(i,229)
        dpdr(i,272) = dpdr(i,2)*p(101) + p(2)*dpdr(i,101) - dpdr(i,215) - dpdr(i,248)
        dpdr(i,273) = dpdr(i,2)*p(106) + p(2)*dpdr(i,106) - dpdr(i,220) - dpdr(i,257)
        dpdr(i,274) = dpdr(i,6)*p(59) + p(6)*dpdr(i,59) - dpdr(i,220) - dpdr(i,215) - dpdr(i,209)
        dpdr(i,275) = dpdr(i,7)*p(60) + p(7)*dpdr(i,60) - dpdr(i,220) - dpdr(i,215) - dpdr(i,211)
        dpdr(i,276) = dpdr(i,5)*p(66) + p(5)*dpdr(i,66) - dpdr(i,215) - dpdr(i,253) - dpdr(i,249) - dpdr(i,251) - &
            & dpdr(i,248)
        dpdr(i,277) = dpdr(i,6)*p(66) + p(6)*dpdr(i,66) - dpdr(i,213) - dpdr(i,252)
        dpdr(i,278) = dpdr(i,7)*p(66) + p(7)*dpdr(i,66) - dpdr(i,214) - dpdr(i,250)
        dpdr(i,279) = dpdr(i,2)*p(104) + p(2)*dpdr(i,104) - dpdr(i,216) - dpdr(i,239)
        dpdr(i,280) = dpdr(i,2)*p(105) + p(2)*dpdr(i,105) - dpdr(i,219) - dpdr(i,248)
        dpdr(i,281) = dpdr(i,1)*p(170) + p(1)*dpdr(i,170) - dpdr(i,244) - dpdr(i,240) - dpdr(i,273)
        dpdr(i,282) = dpdr(i,10)*p(54) + p(10)*dpdr(i,54) - dpdr(i,227)
        dpdr(i,283) = dpdr(i,10)*p(55) + p(10)*dpdr(i,55) - dpdr(i,225)
        dpdr(i,284) = dpdr(i,2)*p(114) + p(2)*dpdr(i,114) - dpdr(i,253) - dpdr(i,276)
        dpdr(i,285) = dpdr(i,1)*p(174) + p(1)*dpdr(i,174) - dpdr(i,250) - dpdr(i,252) - dpdr(i,248) - dpdr(i,276) - &
            & dpdr(i,273)
        dpdr(i,286) = dpdr(i,6)*p(68) + p(6)*dpdr(i,68) - dpdr(i,217) - dpdr(i,261) - dpdr(i,251)
        dpdr(i,287) = dpdr(i,7)*p(68) + p(7)*dpdr(i,68) - dpdr(i,218) - dpdr(i,259) - dpdr(i,249)
        dpdr(i,288) = dpdr(i,2)*p(109) + p(2)*dpdr(i,109) - dpdr(i,221) - dpdr(i,257)
        dpdr(i,289) = dpdr(i,2)*p(116) + p(2)*dpdr(i,116) - dpdr(i,262) - dpdr(i,285)
        dpdr(i,290) = dpdr(i,3)*p(110) + p(3)*dpdr(i,110) - dpdr(i,221) - dpdr(i,266) - dpdr(i,259) - dpdr(i,261) - &
            & dpdr(i,257)
        dpdr(i,291) = dpdr(i,6)*p(70) + p(6)*dpdr(i,70) - dpdr(i,221) - dpdr(i,266) - dpdr(i,260)
        dpdr(i,292) = dpdr(i,7)*p(70) + p(7)*dpdr(i,70) - dpdr(i,221) - dpdr(i,266) - dpdr(i,258)
        dpdr(i,293) = dpdr(i,2)*p(118) + p(2)*dpdr(i,118) - dpdr(i,266) - dpdr(i,290)
        dpdr(i,294) = dpdr(i,1)*p(183) + p(1)*dpdr(i,183) - dpdr(i,266) - dpdr(i,263) - dpdr(i,290)
        dpdr(i,295) = dpdr(i,6)*p(71) + p(6)*dpdr(i,71) - dpdr(i,265)
        dpdr(i,296) = dpdr(i,7)*p(71) + p(7)*dpdr(i,71) - dpdr(i,264)
        dpdr(i,297) = dpdr(i,1)*p(186) + p(1)*dpdr(i,186) - dpdr(i,270)
        dpdr(i,298) = dpdr(i,1)*p(187) + p(1)*dpdr(i,187) - dpdr(i,277) - dpdr(i,278) - dpdr(i,276)
        dpdr(i,299) = dpdr(i,2)*p(115) + p(2)*dpdr(i,115) - dpdr(i,254)
        dpdr(i,300) = dpdr(i,1)*p(189) + p(1)*dpdr(i,189) - dpdr(i,286) - dpdr(i,287) - dpdr(i,285) - dpdr(i,284) - &
            & dpdr(i,298)
        dpdr(i,301) = dpdr(i,2)*p(117) + p(2)*dpdr(i,117) - dpdr(i,263)
        dpdr(i,302) = dpdr(i,18)*p(39) + p(18)*dpdr(i,39) - dpdr(i,282) - dpdr(i,283) - dpdr(i,280)
        dpdr(i,303) = dpdr(i,2)*p(119) + p(2)*dpdr(i,119) - dpdr(i,294)
        dpdr(i,304) = dpdr(i,3)*p(119) + p(3)*dpdr(i,119) - dpdr(i,295) - dpdr(i,296) - dpdr(i,293)
        dpdr(i,305) = dpdr(i,1)*p(194) + p(1)*dpdr(i,194) - dpdr(i,304) - dpdr(i,303)
    enddo


end subroutine evdpdr
!> @brief Evaluate evdbdr for the N4_1A_PIP_2013 potential.
subroutine evdbdr
    !**********************************************************************
    !  The subroutine eliminate the 2-body terms in Bowman's approach
    !**********************************************************************

    implicit none

    double precision :: dbdr, dmdr, dmsdr, dpdr, drdx, dvdr
    integer          :: j

    common /gradt/    dmsdr(6,6),dmdr(6,0:111),dpdr(6,0:305),dvdr(6), &
        & drdx(6,12),dbdr(6,276)

    integer          :: i
    double precision :: db1dr(6,306)

    ! Pass P(0:305) to BM1(1:306)
    do j=1,6
        do i=1,306
            db1dr(j,i) = dpdr(j,i-1)
        enddo
    enddo

    ! Remove unconnected terms and 2-body terms and pass to B(1:276)
    do j=1,6
        dbdr(j,1) = db1dr(j,4)

        do i=2,4
            dbdr(j,i) = db1dr(j,i+4)
        enddo

        dbdr(j,5) = db1dr(j,10)

        do i=6,11
            dbdr(j,i) = db1dr(j,i+6)
        enddo

        dbdr(j,12) = db1dr(j,19)
        dbdr(j,13) = db1dr(j,21)

        do i=14,26
            dbdr(j,i) = db1dr(j,i+9)
        enddo

        dbdr(j,27) = db1dr(j,37)
        dbdr(j,28) = db1dr(j,39)

        do i=29,53
            dbdr(j,i) = db1dr(j,i+12)
        enddo

        dbdr(j,54) = db1dr(j,67)
        dbdr(j,55) = db1dr(j,69)
        dbdr(j,56) = db1dr(j,71)

        do i=57,97
            dbdr(j,i) = db1dr(j,i+16)
        enddo

        dbdr(j,98)  = db1dr(j,115)
        dbdr(j,99)  = db1dr(j,117)
        dbdr(j,100) = db1dr(j,119)

        do i=101,166
            dbdr(j,i) = db1dr(j,i+20)
        enddo

        dbdr(j,167) = db1dr(j,188)
        dbdr(j,168) = db1dr(j,190)
        dbdr(j,169) = db1dr(j,192)
        dbdr(j,170) = db1dr(j,194)

        do i=171,272
            dbdr(j,i) = db1dr(j,i+25)
        enddo

        dbdr(j,273) = db1dr(j,299)
        dbdr(j,274) = db1dr(j,301)
        dbdr(j,275) = db1dr(j,303)
        dbdr(j,276) = db1dr(j,305)
    enddo


end subroutine evdbdr
! Begin
block data prmt

    implicit none

    double precision :: a, b, c, dbdr, dmdr, dmsdr, dpdr, drdx, dvdr, p, r, re, rm, rms

    common /coord/    r(6)
    common /epot/     rms(6),rm(0:111),p(0:305),c(276),b(276)
    common /gradt/    dmsdr(6,6),dmdr(6,0:111),dpdr(6,0:305),dvdr(6), &
        & drdx(6,12),dbdr(6,276)
    common /msprmt/   a,re

    ! Nonlinear parameters of a(0.9A) and re (1.098A)
    data a    /0.9d0/
    data re   /1.098d0/

    ! Linear parameters optimized by the weighted-least square fitting
    data c(  1)   /  0.102403978367d+04/
    data c(  2)   / -0.329677065461d+04/
    data c(  3)   /  0.563720430512d+04/
    data c(  4)   / -0.398579928313d+04/
    data c(  5)   / -0.385965542114d+04/
    data c(  6)   /  0.388707170536d+04/
    data c(  7)   /  0.313073995764d+04/
    data c(  8)   /  0.448635322841d+04/
    data c(  9)   /  0.145177478364d+05/
    data c( 10)   / -0.520482498458d+04/
    data c( 11)   /  0.927577190398d+04/
    data c( 12)   /  0.577514178728d+04/
    data c( 13)   /  0.549924812055d+04/
    data c( 14)   /  0.128789377708d+05/
    data c( 15)   / -0.123915284843d+05/
    data c( 16)   / -0.152516198607d+05/
    data c( 17)   / -0.311653729683d+05/
    data c( 18)   /  0.241802095549d+04/
    data c( 19)   /  0.552823062458d+04/
    data c( 20)   / -0.157606516638d+05/
    data c( 21)   /  0.183999978622d+05/
    data c( 22)   / -0.730065813635d+04/
    data c( 23)   / -0.154254898066d+04/
    data c( 24)   / -0.270137001447d+05/
    data c( 25)   / -0.712833478072d+04/
    data c( 26)   / -0.157507777456d+05/
    data c( 27)   /  0.150777397542d+03/
    data c( 28)   / -0.460118218767d+04/
    data c( 29)   /  0.621757808397d+05/
    data c( 30)   / -0.387325947409d+05/
    data c( 31)   /  0.956015975839d+05/
    data c( 32)   /  0.598731286868d+05/
    data c( 33)   /  0.184378520851d+05/
    data c( 34)   /  0.118525258622d+05/
    data c( 35)   / -0.192034792924d+05/
    data c( 36)   /  0.102203620381d+05/
    data c( 37)   / -0.358701539635d+04/
    data c( 38)   /  0.174513120739d+05/
    data c( 39)   /  0.620986763680d+05/
    data c( 40)   /  0.382171290501d+05/
    data c( 41)   /  0.372876912748d+04/
    data c( 42)   /  0.458049260553d+04/
    data c( 43)   /  0.974878688873d+05/
    data c( 44)   /  0.182211600789d+05/
    data c( 45)   / -0.117103856781d+05/
    data c( 46)   /  0.889576596375d+04/
    data c( 47)   /  0.176322778899d+05/
    data c( 48)   / -0.292083982135d+05/
    data c( 49)   /  0.394476552928d+04/
    data c( 50)   / -0.490908973766d+04/
    data c( 51)   /  0.341307888032d+05/
    data c( 52)   / -0.186725715169d+05/
    data c( 53)   /  0.147822633435d+05/
    data c( 54)   / -0.160310707626d+04/
    data c( 55)   /  0.529073124595d+03/
    data c( 56)   /  0.436589884373d+04/
    data c( 57)   / -0.275977248391d+05/
    data c( 58)   /  0.916659840404d+05/
    data c( 59)   /  0.797381187300d+04/
    data c( 60)   / -0.266127694396d+05/
    data c( 61)   / -0.366120071607d+05/
    data c( 62)   /  0.302613610927d+05/
    data c( 63)   / -0.472403494286d+05/
    data c( 64)   / -0.437423663348d+04/
    data c( 65)   / -0.299548954583d+05/
    data c( 66)   / -0.650332702893d+05/
    data c( 67)   /  0.710278406294d+05/
    data c( 68)   / -0.754784583336d+05/
    data c( 69)   /  0.230385041533d+05/
    data c( 70)   / -0.371044004395d+05/
    data c( 71)   / -0.106153162568d+04/
    data c( 72)   /  0.415993233011d+05/
    data c( 73)   /  0.180415447082d+05/
    data c( 74)   / -0.181546212206d+05/
    data c( 75)   /  0.133142608504d+04/
    data c( 76)   /  0.190256765666d+05/
    data c( 77)   / -0.150450493013d+05/
    data c( 78)   / -0.920409716245d+05/
    data c( 79)   / -0.216561262057d+04/
    data c( 80)   / -0.171625055996d+05/
    data c( 81)   /  0.153782282485d+05/
    data c( 82)   /  0.158621441332d+04/
    data c( 83)   / -0.785674943671d+04/
    data c( 84)   / -0.209994189557d+05/
    data c( 85)   / -0.333435521700d+05/
    data c( 86)   / -0.224157355756d+05/
    data c( 87)   /  0.182745060537d+05/
    data c( 88)   / -0.501427327882d+04/
    data c( 89)   /  0.735584333269d+04/
    data c( 90)   /  0.413504676578d+03/
    data c( 91)   /  0.109674543527d+04/
    data c( 92)   / -0.186740710649d+04/
    data c( 93)   /  0.359599299302d+04/
    data c( 94)   /  0.721485422686d+04/
    data c( 95)   / -0.330055110073d+05/
    data c( 96)   /  0.730992533723d+05/
    data c( 97)   / -0.799352299005d+04/
    data c( 98)   /  0.112075241137d+05/
    data c( 99)   / -0.156320449278d+05/
    data c(100)   / -0.408715780841d+04/
    data c(101)   /  0.268188694785d+06/
    data c(102)   / -0.102983762167d+06/
    data c(103)   /  0.583699882510d+05/
    data c(104)   /  0.260451292047d+05/
    data c(105)   / -0.177541637884d+05/
    data c(106)   / -0.203359237620d+05/
    data c(107)   /  0.480058134079d+05/
    data c(108)   / -0.885961013155d+05/
    data c(109)   /  0.409859892393d+05/
    data c(110)   /  0.325921040295d+05/
    data c(111)   / -0.260852766104d+05/
    data c(112)   /  0.370242414223d+05/
    data c(113)   / -0.664758905617d+04/
    data c(114)   / -0.275618526134d+05/
    data c(115)   / -0.107805359174d+06/
    data c(116)   /  0.287387764086d+05/
    data c(117)   /  0.666908326522d+05/
    data c(118)   / -0.156541620151d+05/
    data c(119)   /  0.336154826289d+05/
    data c(120)   / -0.414314869617d+05/
    data c(121)   /  0.176679423717d+04/
    data c(122)   /  0.577289611931d+05/
    data c(123)   /  0.302830113555d+05/
    data c(124)   / -0.807830452831d+04/
    data c(125)   /  0.219029154199d+03/
    data c(126)   / -0.413421879006d+05/
    data c(127)   / -0.451799646223d+05/
    data c(128)   / -0.159295124884d+05/
    data c(129)   / -0.424907280397d+05/
    data c(130)   / -0.140611779798d+04/
    data c(131)   /  0.568547939978d+03/
    data c(132)   /  0.264987028493d+05/
    data c(133)   / -0.843157148611d+03/
    data c(134)   /  0.275201914701d+05/
    data c(135)   / -0.105393958865d+05/
    data c(136)   / -0.169105249253d+05/
    data c(137)   /  0.104932391725d+05/
    data c(138)   /  0.119833997133d+05/
    data c(139)   /  0.139893856185d+05/
    data c(140)   / -0.181001837340d+03/
    data c(141)   /  0.153415246742d+05/
    data c(142)   /  0.424481598925d+05/
    data c(143)   / -0.197444461986d+05/
    data c(144)   /  0.490620372926d+04/
    data c(145)   / -0.148185485287d+04/
    data c(146)   /  0.624665998166d+04/
    data c(147)   /  0.245340563235d+05/
    data c(148)   /  0.305015590467d+05/
    data c(149)   / -0.611332793065d+04/
    data c(150)   / -0.555190868114d+04/
    data c(151)   /  0.130155396245d+04/
    data c(152)   / -0.455810219138d+05/
    data c(153)   /  0.115493696632d+05/
    data c(154)   /  0.888993953795d+04/
    data c(155)   /  0.103815845402d+05/
    data c(156)   /  0.562575679728d+04/
    data c(157)   / -0.346648114071d+04/
    data c(158)   / -0.240765749446d+04/
    data c(159)   /  0.206433652571d+04/
    data c(160)   / -0.758049531330d+04/
    data c(161)   /  0.171657661399d+05/
    data c(162)   / -0.260208297109d+04/
    data c(163)   / -0.372114666463d+04/
    data c(164)   /  0.192848747262d+05/
    data c(165)   / -0.614612031387d+05/
    data c(166)   /  0.284938866770d+04/
    data c(167)   /  0.580215877101d+03/
    data c(168)   / -0.142884180958d+05/
    data c(169)   /  0.164345541019d+05/
    data c(170)   /  0.204700489332d+04/
    data c(171)   /  0.232847535995d+06/
    data c(172)   / -0.215988386103d+06/
    data c(173)   / -0.543268059091d+06/
    data c(174)   / -0.291087990782d+06/
    data c(175)   /  0.993671706214d+05/
    data c(176)   / -0.324811444690d+05/
    data c(177)   / -0.165451102346d+06/
    data c(178)   /  0.503736859146d+05/
    data c(179)   /  0.441027718212d+05/
    data c(180)   /  0.427514011421d+05/
    data c(181)   /  0.316984455059d+05/
    data c(182)   / -0.916644985736d+05/
    data c(183)   / -0.497678170749d+05/
    data c(184)   / -0.570481120151d+04/
    data c(185)   /  0.392692240737d+05/
    data c(186)   / -0.636234723414d+05/
    data c(187)   /  0.605343078929d+05/
    data c(188)   / -0.504945371764d+05/
    data c(189)   / -0.394790420088d+05/
    data c(190)   / -0.296012156034d+05/
    data c(191)   / -0.189032641400d+05/
    data c(192)   /  0.315573783036d+05/
    data c(193)   / -0.123145782845d+05/
    data c(194)   / -0.183417710262d+05/
    data c(195)   /  0.137514063372d+04/
    data c(196)   /  0.282584683995d+05/
    data c(197)   /  0.255488952485d+03/
    data c(198)   / -0.214177073532d+05/
    data c(199)   / -0.329389386184d+05/
    data c(200)   /  0.326844593866d+05/
    data c(201)   / -0.268690289613d+05/
    data c(202)   /  0.688737993203d+03/
    data c(203)   / -0.174784217735d+05/
    data c(204)   /  0.592388749822d+04/
    data c(205)   /  0.782006070657d+04/
    data c(206)   /  0.182368489092d+05/
    data c(207)   / -0.195025918004d+05/
    data c(208)   /  0.140155756116d+05/
    data c(209)   /  0.411097619946d+05/
    data c(210)   /  0.363810881096d+03/
    data c(211)   / -0.280890295507d+04/
    data c(212)   / -0.105569797133d+05/
    data c(213)   /  0.380998787463d+04/
    data c(214)   / -0.360964430447d+04/
    data c(215)   / -0.986412341243d+04/
    data c(216)   /  0.226662646459d+05/
    data c(217)   / -0.376753011539d+04/
    data c(218)   /  0.578159426611d+03/
    data c(219)   / -0.121449139228d+05/
    data c(220)   / -0.750765164465d+04/
    data c(221)   / -0.140227007314d+05/
    data c(222)   /  0.203120231276d+04/
    data c(223)   /  0.118811735987d+05/
    data c(224)   /  0.519651005890d+04/
    data c(225)   /  0.947714025376d+04/
    data c(226)   /  0.179018909494d+05/
    data c(227)   /  0.121879281662d+05/
    data c(228)   / -0.126436586574d+05/
    data c(229)   / -0.452354749524d+04/
    data c(230)   / -0.684935604354d+04/
    data c(231)   / -0.174197586324d+04/
    data c(232)   / -0.238302572429d+04/
    data c(233)   /  0.537465564700d+04/
    data c(234)   / -0.157114135550d+03/
    data c(235)   / -0.652324696358d+04/
    data c(236)   / -0.323096366831d+04/
    data c(237)   /  0.296083943992d+04/
    data c(238)   / -0.122452541606d+04/
    data c(239)   / -0.593094058930d+04/
    data c(240)   / -0.379005180042d+04/
    data c(241)   /  0.777078060474d+04/
    data c(242)   /  0.193407389171d+04/
    data c(243)   / -0.251175830196d+04/
    data c(244)   / -0.342507311002d+05/
    data c(245)   / -0.817218571265d+04/
    data c(246)   /  0.167151796682d+04/
    data c(247)   /  0.140032630687d+04/
    data c(248)   / -0.108222853767d+04/
    data c(249)   / -0.289156752055d+04/
    data c(250)   /  0.671923615356d+04/
    data c(251)   / -0.186517211143d+04/
    data c(252)   / -0.608928447848d+04/
    data c(253)   / -0.320885381491d+04/
    data c(254)   /  0.168580430914d+03/
    data c(255)   /  0.561938258944d+03/
    data c(256)   /  0.178839575085d+04/
    data c(257)   / -0.127093039569d+02/
    data c(258)   /  0.158954991882d+05/
    data c(259)   / -0.728671162222d+04/
    data c(260)   /  0.273935428053d+03/
    data c(261)   /  0.210205486354d+03/
    data c(262)   / -0.726829028390d+04/
    data c(263)   /  0.234342832052d+04/
    data c(264)   /  0.450306951911d+03/
    data c(265)   / -0.225605345434d+04/
    data c(266)   /  0.146430478883d+04/
    data c(267)   / -0.502512737258d+04/
    data c(268)   /  0.529554234299d+03/
    data c(269)   /  0.665683392527d+03/
    data c(270)   / -0.446965631596d+04/
    data c(271)   /  0.155396579789d+05/
    data c(272)   / -0.690692642506d+03/
    data c(273)   /  0.128917160396d+04/
    data c(274)   /  0.373090368655d+04/
    data c(275)   / -0.454736355855d+04/
    data c(276)   / -0.377918615671d+03/

end
