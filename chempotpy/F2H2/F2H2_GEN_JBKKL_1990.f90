!> @modern CHEMPOTPY modern Fortran style is enforced for this file.


!> @brief Evaluate pes for the F2H2_GEN_JBKKL_1990 potential.
subroutine pes( x, igrad, path, p, g, d )

    implicit none

    integer :: iatom, idir, istate
    ! number of electronic state
    integer, parameter              :: nstates = 1
    integer, parameter              :: natoms  = 4
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    double precision, intent(in)    :: x(natoms,3)
    double precision, intent(out)   :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out)   :: d(nstates,nstates,natoms,3)

    double precision :: cartnu(natoms,3), r(6), v
    double precision :: rcm, r1, r2
    logical, save    :: first_time_data=.true.

    !initialize
    v = 0.d0
    g = 0.d0
    d = 0.d0

    cartnu = 0.d0
    do iatom=1,natoms
        do idir=1,3
            cartnu(iatom,idir) = x(iatom,idir)/0.529177211
        enddo
    enddo

    if(first_time_data) then
        call bjinit(path)
        first_time_data = .false.
    endif


    call carttor(cartnu,r)
    rcm = r(1)
    r1  = r(2)
    r2  = r(3)

    call trbunk(rcm,r1,r2,v,4,0,0)

    if (igrad==0) then
        do istate=1,nstates
            p(istate) = v*27.211386
        enddo
    else
        write (*,*) 'Only energy is available'
    endif

end subroutine pes

!> @brief Evaluate CARTTOR for the F2H2_GEN_JBKKL_1990 potential.
subroutine carttor( cartnu, r )
    implicit none

    double precision :: cmm, cphi, flm, hhd, hym, phi, pi, q, q1, q2, q3, theta1, theta2, xcm1, xcm2, &
        & xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
    double precision :: cartnu(4,3), r(6)
    flm    = 18.99840d0
    hym    = 1.007825d0
    xcm1   = (hym*cartnu(1,1)+flm*cartnu(2,1))/(flm+hym)
    ycm1   = (hym*cartnu(1,2)+flm*cartnu(2,2))/(flm+hym)
    zcm1   = (hym*cartnu(1,3)+flm*cartnu(2,3))/(flm+hym)
    xcm2   = (hym*cartnu(3,1)+flm*cartnu(4,1))/(flm+hym)
    ycm2   = (hym*cartnu(3,2)+flm*cartnu(4,2))/(flm+hym)
    zcm2   = (hym*cartnu(3,3)+flm*cartnu(4,3))/(flm+hym)
    xcm3   = xcm2-xcm1
    ycm3   = ycm2-ycm1
    zcm3   = zcm2-zcm1
    xrm1   = cartnu(1,1)-xcm1
    yrm1   = cartnu(1,2)-ycm1
    zrm1   = cartnu(1,3)-zcm1
    theta1 = (xrm1*xcm3+yrm1*ycm3+zrm1*zcm3)
    theta1 = theta1/(sqrt(xrm1**2+yrm1**2+zrm1**2))
    theta1 = theta1/(sqrt(xcm3**2+ycm3**2+zcm3**2))
    if(theta1.gt.1.0d0)theta1=1.0d0
    if(theta1.lt.-1.0d0)theta1=-1.0d0
    theta1 = acos(theta1)
    xrm2   = cartnu(3,1)-xcm2
    yrm2   = cartnu(3,2)-ycm2
    zrm2   = cartnu(3,3)-zcm2
    theta2 = (xrm2*(-xcm3)+yrm2*(-ycm3)+zrm2*(-zcm3))
    theta2 = theta2/(sqrt(xrm2**2+yrm2**2+zrm2**2))
    theta2 = theta2/(sqrt(xcm3**2+ycm3**2+zcm3**2))
    if(theta2.gt.1.0d0)theta2=1.0d0
    if(theta2.lt.-1.0d0)theta2=-1.0d0
    theta2 = acos(theta2)
    pi     = acos(-1.0d0)
    theta2 = pi-theta2
    q1     = sqrt(xrm1**2+yrm1**2+zrm1**2)
    q2     = sqrt(xrm2**2+yrm2**2+zrm2**2)
    cmm    = (xcm3**2+ycm3**2+zcm3**2)
    cmm    = sqrt(cmm)
    hhd    = (cartnu(1,1)-cartnu(3,1))**2 + (cartnu(1,2)-cartnu(3,2))**2 + (cartnu(1,3)-cartnu(3,3))**2
    hhd    = sqrt(hhd)
    q      = cmm-q1*cos(theta1)+q2*cos(theta2)
    q3     = sqrt(abs(hhd**2-q**2))
    q1     = q1*sin(theta1)
    q2     = q2*sin(theta2)
    cphi   = (q1**2+q2**2-q3**2)/(2.*q1*q2)
    if(cphi.lt.-1.0d0)cphi=-1.0d0
    if(cphi.gt.1.0d0)cphi=1.0d0
    phi  = acos(cphi)
    r(1) = sqrt(xcm3**2+ycm3**2+zcm3**2)
    r(2) = (sqrt(xrm1**2+yrm1**2+zrm1**2))*(flm+hym)/flm
    r(3) = (sqrt(xrm2**2+yrm2**2+zrm2**2))*(flm+hym)/flm
    r(4) = theta1
    r(5) = theta2
    r(6) = phi
endsubroutine carttor

! bjinit reads in parameters for BJKKL potential from file bjkkl.par
!> @brief Evaluate bjinit for the F2H2_GEN_JBKKL_1990 potential.
subroutine bjinit( path )
    implicit none
    integer :: file_unit_11

    double precision                :: c
    integer                         :: i
    character(len=1024), intent(in) :: path
    character(len=1024)             :: file_path1
    dimension c(359)
    common/params/c


    file_path1 = trim(path)//"/F2H2/bjkkl.par"
    open(newunit=file_unit_11,file=file_path1)
    read(file_unit_11,'(5F14.10)')(c(i),i=1,359)
    close(file_unit_11)

end subroutine bjinit
!
!
!
!  bjinit should be called before the first call to trbunk
!
!> @brief Evaluate trbunk for the F2H2_GEN_JBKKL_1990 potential.
subroutine trbunk( r, r1, r2, v, iq1, iq2, imu )

    !     subroutine trbunk transforms the surface of Bunker, Jensen, Kofran
    !     Lischka and Karpfen(the angular part only; no two-body terms) to t
    !     representation of Launay et al.

    !     references

    !     Bunker potential: P.R. Bunker, P. Jensen, A. Karpfen, M. Kofranek,
    !     H. Lischka, J. Chem. Phys. 92, 7432 (1990).
    !     Launay expansion: J.M. Launay, J. Phys. B 10, 3665 (1977).

    implicit none

    double precision :: a, bo, c, cdenom, const, const1, const2, const3, ha, pc, phas, pi, r, r1, r2, threej, v, vl, xl, &
        & xl1, xl2, xm
    integer :: imu, iq1, iq2, l, ll, lm, lu
    save
    integer, parameter          :: lmax    = 4
    integer, parameter          :: lsummax = 9
    integer, parameter          :: ipara   = 1
    double precision, parameter :: zero    = 0.0d0
    double precision, parameter :: one     = 1.d0
    double precision, parameter :: two     = 2.d0
    double precision, parameter :: four    = 4.d0
    integer, parameter          :: izero   = 0
    !     dimension a(0:4,0:4,0:8),v(ipara)
    !     dimension iq1(ipara),iq2(ipara),imu(ipara)
    dimension a(0:4,0:4,0:8)
    dimension c(359)
    common/acoef/a
    common/params/c
    common/consts/pi,bo,ha,vl,pc

    bo = 0.529177249d0
    ha = 219474.63067d0
    pi = 3.1415926535898d0
    vl = 2.99792458d0
    pc = 6.6260755d0

    !     next get the a's from subroutine vcalc

    call vcalc(r1,r2,r)

    xl1 = float(iq1)
    xl2 = float(iq2)
    xm  = float(imu)
    lm  = min(iq1,iq2)
    ll  = abs(iq1-iq2)
    lu  = iq1+iq2
    v   = 0.0d0

    do 20 l=ll,lu,2
        xl     = float(l)
        phas   = (-1)**(iq1-iq2)
        const1 = one
        if(imu.eq.0)const1=two
        const1 = one/dsqrt(const1)
        const2 = two*xl+one
        cdenom = (four*pi)**3
        const3 = dsqrt(two/cdenom)
        const  = phas*const1*const2*const3
        v      = v+a(iq1,iq2,l)*const*threej(xl1,xl2,xl,xm,-xm,0.d0)
    20 continue

    !     divide next two terms by three because there is a bug in Bunker's
    !     function subprogram "P" which computes the associated Legendre
    !     function. (this is a very small correction).

    if(iq1.eq.4.and.iq2.eq.2.and.imu.eq.2)v=v/3.d0
    if(iq1.eq.2.and.iq2.eq.4.and.imu.eq.2)v=v/3.d0

    !     vpot=zero
    !     kount=kount-1
    !     do 50 j=1,kount
    !     write(65,*)v(j)
    !50   continue

end subroutine trbunk
!> @brief Evaluate THREEJ for the F2H2_GEN_JBKKL_1990 potential.
function threej( a, b, c, al, be, ga )
    !
    !     CALCULATES WIGNER 3-J SYMBOLS USING FORMULA OF P34 OF BRINK &
    !     SATCHLER IN "ANGULAT MOMENTUM"
    !     NOTE FAKT(A)=FACT(A)/10**A
    !     A=J1,B=J2,C=J,AL=M1,BE=M2,GA=-M
    !
    !
    implicit none

    double precision :: a, al, b, be, c, clebsh, delta, dnew, fakt, ga, summ, term, term1, term2, &
        & term3, term4, term5, term6, threej
    integer :: iphase, new, newmax, newmin
    save
    !
    threej = 0.d0
    !
    !     (J1+J2).GE.J AND J.GE.ABS(A-B)    -M=M1+M2    J1,J2,J.GE.0
    !     ABS(M1).LE.J1    ABS(M2).LE.J2   ABS(M).LE.J
    !
    if(c.gt.a+b)return
    if(c.lt.dabs(a-b))return
    if(a.lt.0.d0.or.b.lt.0.d0.or.c.lt.0.d0)return
    if(a.lt.dabs(al).or.b.lt.dabs(be).or.c.lt.dabs(ga))return
    if(-1.d0*ga.ne.al+be)return
    !
    !
    !     COMPUTE DELTA(ABC)
    !
    delta = dsqrt(fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.d0))
    !
    !
    term1 = fakt(a+al)*fakt(a-al)
    term2 = fakt(b-be)*fakt(b+be)
    term3 = fakt(c+ga)*fakt(c-ga)
    term  = dsqrt((2.d0*c+1.d0)*term1*term2*term3)
    !
    !
    !     NOW COMPUTE SUMMATION TERM
    !
    !     SUM TO GET SUMMATION IN EQ(2.34) OF BRINK AND SATCHLER.  SUM UNTIL
    !     A TERM INSIDE FACTORIAL GOES NEGATIVE.  NEW IS INDEX FOR SUMMATION
    !     .  NOW FIND WHAT THE RANGE OF NEW IS.
    !
    !
    newmin = idnint(dmax1((a+be-c),(b-c-al),0.d0))
    newmax = idnint(dmin1((a-al),(b+be),(a+b-c)))
    !
    !
    summ = 0.d0
    !
    !
    do 10 new=newmin,newmax
        dnew  = dfloat(new)
        term4 = fakt(a-al-dnew)*fakt(c-b+al+dnew)
        term5 = fakt(b+be-dnew)*fakt(c-a-be+dnew)
        term6 = fakt(dnew)*fakt(a+b-c-dnew)
        summ  = summ+(-1.d0)**new/(term4*term5*term6)
    10 continue
    !
    !     SO CLEBSCH-GORDEN <J1J2M1M2LJM> IS CLEBSH
    !
    clebsh = delta*term*summ/dsqrt(10.d0)
    !
    !     CONVERT CLEBSCH-GORDEN TO THREEJ
    !
    iphase = idnint(a-b-ga)
    threej = (-1.d0)**(iphase)*clebsh/dsqrt(2.d0*c+1.d0)
    !
    !
end function threej
!> @brief Evaluate VCALC for the F2H2_GEN_JBKKL_1990 potential.
subroutine vcalc( r1, r2, rab )
    !
    !  CALCULATES V AT POINT (R1,R2,RAB,THETA1,THETA2,TAU)
    !  USING THE CURRENT VALUES OF THE PARAMETERS.
    !    V CALCULATED USING GREEN'S EXPRESSION; SEE JCP 62,2271 (1975).
    !       G IS THE SCRIPT I USED BY GREEN IN HIS EQ.(A9)
    !
    implicit none

    double precision :: a, bo, c, cd, cdq, cq, d1, d2, d3, de1, de2, de3, de4, g, ha, pc, pi, r1, r2, &
        & r6, rab, req, rin01, rin10, rmiddle, rout10, t1, t1deg, t2, t2deg, tau, taudeg, vl, y1, &
        & y1py2, y1y2, y1y22, y1y23, y1y24, y2, y3
    integer :: i, il, iu, j, k
    save
    integer :: out
    dimension c(359),g(0:4,0:4,0:8),a(0:4,0:4,0:8),req(6)
    common/params/c
    common/consts/pi,bo,ha,vl,pc
    common/acoef/a
    !
    out = 6
    !
    !     R1=REQ(1)
    !     R2=REQ(2)
    !     RAB=REQ(3)
    !     T1DEG=REQ(4)
    !     T2DEG=REQ(5)
    !     TAUDEG=REQ(6)
    !     T1=T1DEG*PI/180.D0
    !     T2=T2DEG*PI/180.D0
    !     TAU=TAUDEG*PI/180.D0
    !
    !
    !     NEXT CALL SUBROUTINE HF2INT WHICH CALCULATES THE RELEVANT ATOM-ATO
    !     DISTANCES NEEDED TO COMPUTE THE TWO-BODY TERMS. THE ANGLES MUST
    !     BE IN DEGREES SO CONVERT THEM HERE.                           WCN
    t1deg  = t1*180.d0/pi
    t2deg  = t2*180.d0/pi
    taudeg = tau*180.d0/pi
    !     CALL HF2INT(R1,RAB,R2,T1DEG,T2DEG,TAUDEG,RHH,RHF1,RHF2,RFF)
    !
    !     DO 100 L1=0,4
    !     DO 100 L2=0,4
    !     LL=IABS(L1-L2)
    !     LU=L1+L2
    !     LM=L2
    !     IF(L1.LT.L2)LM=L1
    !     DO 100 L=LL,LU,2
    !     XL=DFLOAT(L)
    !     ZZ=((2.D0*XL+1.D0)/2.D0)/DSQRT(PI)
    !     XL1=DFLOAT(L1)
    !     XL2=DFLOAT(L2)
    !     COEFF=((-1)**(L1-L2))*ZZ
    !     G(L1,L2,L)=THREEJ(XL1,XL2,XL,0.D0,0.D0,0.D0)*P(L1,0,T1)*
    !    X           P(L2,0,T2)
    !     IF(LM.EQ.0)GO TO 100
    !       DO 90 M=1,LM
    !       XM=DFLOAT(M)
    ! 90    G(L1,L2,L)=G(L1,L2,L)+
    !    X             ((-1)**M)*2.0D0*THREEJ(XL1,XL2,XL,XM,-XM,0.D0)*
    !    X             P(L1,M,T1)*P(L2,M,T2)*DCOS(XM*TAU)
    !100  G(L1,L2,L)=COEFF*G(L1,L2,L)
    !
    do 20 i=0,4
        do 20 j=0,4

            do 20 k=0,8
    20 a(i,j,k)=0.0d0
    !
    d1  = r1-c(1)
    d2  = r2-c(1)
    d3  = rab-c(2)
    y1  = 1.d0-dexp(-c(3)*d1)
    y2  = 1.d0-dexp(-c(3)*d2)
    y3  = 1.d0-dexp(-d3)
    de1 = dexp(-d3)
    de2 = dexp(-2.0d0*d3)
    de3 = dexp(-3.0d0*d3)
    de4 = dexp(-4.0d0*d3)
    r6  = rab**6
    !
    !  DAMPING FACTORS AROUND RIN AND ROUT
    !
    rout10  = 0.5d0*(1.d0-dtanh(0.5d0*(rab-7.5d0)))
    rin10   = 0.5d0*(1.d0-dtanh(2.d0*(rab-4.0d0)))
    rin01   = 0.5d0*(1.d0-dtanh(2.d0*(4.0d0-rab)))
    rmiddle = rin01*rout10
    !
    y1y22 = y1**2+y2**2
    y1y23 = y1**3+y2**3
    y1y24 = y1**4+y2**4
    y1py2 = y1+y2
    y1y2  = y1*y2
    !
    !    CD, CQ AND CMQ ARE CONVERSION FACTORS FOR MULTIPOLE MOMENTS
    !
    cd  = (bo**3)*pc*vl*ha/100000.d0
    cq  = cd*bo*bo
    cdq = cd*bo
    !
    !    C(4) TO C(9), C(91) AND C(111) DESCRIBE THE DISPERSION INTERACTION
    !    C(9) IS IN A(2,2,4).
    !
    !    C(41) IS THE HF DIPOLE MOMENT IN DEBYE
    !    C(101) IS THE HF QUADRUPOLE MOMENT IN DEBYE.ANGSTROM
    !    C(131) IS THE HF OCTUPOLE MOMENT IN DEBYE.ANGSTROM**2
    !
    !
    !     the zero of energy is chosen as the two monomers infinitely separa
    !     at r1=r2=1.744 a.u. thus c(11) has been omitted from a(0,0,0).
    !

    a(0,0,0)=c(357)*de2*rin10 +c(20)*y1y2*rout10 +(c(12)*de1+c(13)*de2+c(14)*de3+c(15)*de4)*rmiddle +(c(4)/r6 + c(5)/ &
        & (rab**7) + c(6)/(rab**8) +c(7)/(rab**9)+c(8)/(rab**10))*rin01 +(c(157)*(y1+y2)*y3 +c(158)*y1*y2*y3 +c(166)* &
        & y1* y2*y3**2 +c(167)*y1*y2*y3**3)*rmiddle

    a(1,0,1)=(c(22)*de1+c(23)*de2 +c(24)*de3+c(25)*de4 +c(26)*y1 +c(27)*y2 +c(21)*y1*y3 +c(156)*y2*y3 +c(28)*y1**2 + &
        & c(29)*y2**2 +c(30)*y1*y2 +c(116)*y1**3 +c(117)*y2**3 +c(118)*y1*y2**2 +c(119)*y2*y1**2 +c(159)*y1*y2*y3 + &
        & c(169)*y1*y2*y3**2 +c(179)*y1*y2**2*y3)*rmiddle
    a(0,1,1)=-(c(22)*de1+c(23)*de2 +c(24)*de3+c(25)*de4 +c(26)*y2 +c(27)*y1 +c(21)*y2*y3 +c(156)*y1*y3 +c(28)*y2**2 + &
        & c(29)*y1**2 +c(30)*y1*y2 +c(116)*y2**3 +c(117)*y1**3 +c(118)*y2*y1**2 +c(119)*y1*y2**2 +c(159)*y1*y2*y3 + &
        & c(169)*y1*y2*y3**2 +c(179)*y1**2*y1*y3)*rmiddle
    !
    a(1,1,0)=(c(32)*de1+c(33)*de2 +c(34)*de3+c(35)*de4 +c(36)*y1y22 +c(37)*y1y23 +c(38)*y1py2 +c(39)*y1py2*y1y2 + &
        & c(40)* y1y2)*rmiddle
    a(1,1,2)=(c(42)*de1+c(43)*de2 +c(44)*de3+c(45)*de4)*rmiddle -((8.d0*pi**1.5d0*dsqrt(2.d0/15.d0))*(c(41)**2)/(cd* &
        & rab**3))*rin01 +(c(46)*y1y22 +c(47)*y1y23 +c(48)*y1py2 +c(49)*y1py2*y1y2 +c(50)*y1y2)*rmiddle
    !
    a(2,0,2)=c(356)/r6+(c(52)*de1+c(53)*de2 +c(54)*de3+c(55)*de4 +c(56)*y1 +c(57)*y2 +c(58)*y1**2 +c(59)*y2**2 +c(60)* &
        & y1*y2 +c(126)*y1**3 +c(127)*y2**3 +c(128)*y1*y2**2 +c(129)*y2*y1**2)*rmiddle
    a(0,2,2)=c(356)/r6+(c(52)*de1+c(53)*de2 +c(54)*de3+c(55)*de4 +c(56)*y2 +c(57)*y1 +c(58)*y2**2 +c(59)*y1**2 +c(60)* &
        & y1*y2 +c(126)*y2**3 +c(127)*y1**3 +c(128)*y2*y1**2 +c(129)*y1*y2**2)*rmiddle
    !
    a(2,1,1)=(c(62)*de1+c(63)*de2 +c(64)*de3+c(65)*de4 +c(66)*y1 +c(67)*y2 +c(68)*y1**2 +c(69)*y2**2 +c(70)*y1*y2 + &
        & c(136)*y1**3 +c(137)*y2**3 +c(138)*y1*y2**2 +c(139)*y2*y1**2)*rmiddle
    a(1,2,1)=-(c(62)*de1+c(63)*de2 +c(64)*de3+c(65)*de4 +c(66)*y2 +c(67)*y1 +c(68)*y2**2 +c(69)*y1**2 +c(70)*y1*y2 + &
        & c(136)*y2**3 +c(137)*y1**3 +c(138)*y2*y1**2 +c(139)*y1*y2**2)*rmiddle
    a(2,1,3)=(c(72)*de1+c(73)*de2 +c(74)*de3+c(75)*de4)*rmiddle -rin01*(8.d0*pi**1.5d0/dsqrt(7.d0))*c(41)*c(101)/(cdq* &
        & rab**4) +(c(76)*y1 +c(77)*y2 +c(78)*y1**2 +c(79)*y2**2 +c(80)*y1*y2 +c(146)*y1**3 +c(147)*y2**3 +c(148)*y1* &
        & y2**2 +c(149)*y2*y1**2)*rmiddle
    a(1,2,3)=-(c(72)*de1+c(73)*de2 +c(74)*de3+c(75)*de4)*rmiddle +rin01*(8.d0*pi**1.5d0/dsqrt(7.d0))*c(41)*c(101)/ &
        & (cdq* rab**4) -(c(76)*y2 +c(77)*y1 +c(78)*y2**2 +c(79)*y1**2 +c(80)*y1*y2 +c(146)*y2**3 +c(147)*y1**3 + &
        & c(148)*y2* y1**2 +c(149)*y1*y2**2)*rmiddle
    !
    a(2,2,0)=(c(82)*de1+c(83)*de2 +c(84)*de3+c(85)*de4 +c(86)*y1y22 +c(87)*y1y23 +c(88)*y1py2 +c(89)*y1py2*y1y2 + &
        & c(90)* y1y2)*rmiddle
    a(2,2,2)=c(91)/r6 +(c(92)*de1+c(93)*de2 +c(94)*de3+c(95)*de4 +c(96)*y1y22 +c(97)*y1y23 +c(98)*y1py2 +c(99)*y1py2* &
        & y1y2 +c(100)*y1y2)*rmiddle
    a(2,2,4)=(c(102)*de1+c(103)*de2 +c(104)*de3+c(105)*de4)*rmiddle +(8.d0*pi**1.5d0/3.d0)*dsqrt(14.d0/5.d0)*rin01* &
        & (c(101)**2)/(cq*rab**5) +(c(106)*y1y22 +c(107)*y1y23 +c(108)*y1py2 +c(109)*y1py2*y1y2 +c(110)*y1y2)*rmiddle
    !
    a(3,0,3) = c(111)/(rab**7) +(c(112)*de1+c(113)*de2 +c(114)*de3+c(115)*de4)*rmiddle
    a(0,3,3) = -a(3,0,3)
    a(3,1,2) = (c(122)*de1+c(123)*de2 +c(124)*de3+c(125)*de4)*rmiddle
    a(1,3,2) = a(3,1,2)
    a(3,1,4)=(c(132)*de1+c(133)*de2 +c(134)*de3+c(135)*de4)*rmiddle -rin01*(16.d0*((pi/3.d0)**1.5d0))*c(41)*c(131)/ &
        & (cq* rab**5)
    a(1,3,4) = a(3,1,4)
    a(3,2,1) = (c(142)*de1+c(143)*de2 +c(144)*de3+c(145)*de4)*rmiddle
    a(2,3,1) = -a(3,2,1)
    a(3,2,3) = (c(152)*de1+c(153)*de2 +c(154)*de3+c(155)*de4)*rmiddle
    a(2,3,3) = -a(3,2,3)
    a(3,2,5) = (c(162)*dexp(-0.5d0*d3)+c(163)*dexp(-0.75d0*d3) +c(164)*de1+c(165)*de4)*rmiddle
    a(2,3,5) = -a(3,2,5)
    a(3,3,0) = (c(172)*de1+c(173)*de2 +c(174)*de3+c(175)*de4 +c(176)*y1y22+c(177)*y1y2)*rmiddle
    a(3,3,2) = (c(182)*de1+c(183)*de2 +c(184)*de3+c(185)*de4 +c(186)*y1y22+c(187)*y1y2)*rmiddle
    a(3,3,4) = (c(192)*de1+c(193)*de2 +c(194)*de3+c(195)*de4 +c(196)*y1y22+c(197)*y1y2)*rmiddle
    a(3,3,6) = (c(202)*de1+c(203)*de2 +c(204)*de3+c(205)*de4 +c(206)*y1y22+c(207)*y1y2)*rmiddle
    !
    a(4,0,4) = c(212)*de1+c(213)*de2 +c(214)*de3+c(215)*de4
    a(0,4,4) = a(4,0,4)
    a(4,1,3) = c(222)*de1+c(223)*de2 +c(224)*de3+c(225)*de4
    a(1,4,3) = -a(4,1,3)
    a(4,1,5) = c(232)*dexp(-0.5d0*d3)+c(233)*dexp(-0.75d0*d3) +c(234)*de1+c(235)*de2
    a(1,4,5) = -a(4,1,5)
    a(4,2,2) = c(242)*de1+c(243)*de2 +c(244)*de3+c(245)*de4
    a(2,4,2) = a(4,2,2)
    a(4,2,4) = c(252)*de1+c(253)*de2 +c(254)*de3+c(255)*de4
    a(2,4,4) = a(4,2,4)
    a(4,2,6) = c(262)*dexp(-0.5d0*d3)+c(263)*dexp(-0.75d0*d3) +c(264)*de1+c(265)*de2
    a(2,4,6) = a(4,2,6)
    a(4,3,1) = c(272)*de1+c(273)*de2 +c(274)*de3+c(275)*de4
    a(3,4,1) = -a(4,3,1)
    a(4,3,3) = c(282)*de1+c(283)*de2 +c(284)*de3+c(285)*de4
    a(3,4,3) = -a(4,3,3)
    a(4,3,5) = c(292)*de1+c(293)*de2 +c(294)*de3+c(295)*de4
    a(3,4,5) = -a(4,3,5)
    a(4,3,7) = c(302)*de1+c(303)*de2 +c(304)*de3+c(305)*de4
    a(3,4,7) = -a(4,3,7)
    a(4,4,0) = c(312)*de1+c(313)*de2 +c(314)*de3+c(315)*de4 +(c(316)*y1y22+c(317)*y1y2)
    a(4,4,2) = c(322)*de1+c(323)*de2 +c(324)*de3+c(325)*de4
    a(4,4,4) = c(332)*de1+c(333)*de2 +c(334)*de3+c(335)*de4
    a(4,4,6) = c(342)*de1+c(343)*de2 +c(344)*de3+c(345)*de4
    a(4,4,8) = c(352)*de1+c(353)*de2 +c(354)*de3+c(355)*de4
    !
    do 992 i=0,4
        do 992 j=0,4
            il = iabs(i-j)
            iu = i+j
            do 992 k=il,iu,2
                if(i.ne.4.and.j.ne.4)go to 992
                a(i,j,k) = a(i,j,k)*rmiddle
    992 continue
    !
    !     V = 0.0D0
    !     DO 2 I=0,4
    !     DO 2 J=0,4
    !     IL=IABS(I-J)
    !     IU=I+J
    !     DO 2 K=IL,IU,2
    !2    V = V + A(I,J,K)*G(I,J,K)
    !
    !  HH, FF, HF1 AND HF2 2-BODY TERMS COMBINED IN TERM VSRG.
    !  VSRG IS DAMPED TO ZERO IF RHF.GT.3.0
    !  IN ORDER TO PREVENT SPURIOUS 'HOLES' IN V.
    !  DAMPED TO ZERO IF RHH.GT.2.3, RHF1/2.GT.2.3, AND RFF.GT.2.3
    !  SO THAT BOUND STATES CAN BE CALCULATED WITH THE NEGLECT
    !  OF VSRG AND CLOSE COUPLING PROGRAM 'BOUND' USED WITH V.
    !  eliminate two-body terms completely in this version wcn 3/10/92
    !
    !     YRHF=0.25D0*(1.D0-DTANH(7.D0*(R1-3.0D0)))
    !    X           *(1.D0-DTANH(7.D0*(R2-3.0D0)))
    !     YRHH=0.50D0*(1.D0-DTANH(7.D0*(RHH-2.3D0)))
    !     YRFF=0.50D0*(1.D0-DTANH(7.D0*(RFF-2.3D0)))
    !     YRHF1=0.50D0*(1.D0-DTANH(7.D0*(RHF1-2.3D0)))
    !     YRHF2=0.50D0*(1.D0-DTANH(7.D0*(RHF2-2.3D0)))
    !
    !     VHH=C(161)*DEXP(-C(171)*RHH)+C(181)*DEXP(-C(191)*RHH)
    !     VHF1=C(201)*DEXP(-C(211)*RHF1)+C(221)*DEXP(-C(231)*RHF1)
    !     VHF2=C(201)*DEXP(-C(211)*RHF2)+C(221)*DEXP(-C(231)*RHF2)
    !     VFF=C(241)*DEXP(-C(251)*RFF)+C(261)*DEXP(-C(271)*RFF)
    !
    !     VSRG=(VHH*YRHH+VHF1*YRHF1+VHF2*YRHF2+VFF*YRFF)*YRHF
    !
    !     V=V+VSRG
    !
end subroutine vcalc
!> @brief Evaluate FAKT for the F2H2_GEN_JBKKL_1990 potential.
function fakt( a )
    !     THIS FUNCTION COMPUTES N!/(10**N).  STOP 101 IMPLIES TRIED
    !     TO TAKE FACTORIAL OF A NEGATIVE NO.
    !
    implicit none

    double precision :: a, ax, fakt
    integer          :: i, ic
    save
    !
    ax   = a
    fakt = 1.d0
    if(ax.eq.0.d0)return
    fakt = .1d0
    if(ax.lt.0.d0)stop 101
    !
    ic   = idnint(ax)
    ax   = ax/10.d0
    fakt = ax
    do 10  i=1,ic-1
        fakt = fakt*(ax-dfloat(i)*.1d0)
    10 continue
end function fakt
