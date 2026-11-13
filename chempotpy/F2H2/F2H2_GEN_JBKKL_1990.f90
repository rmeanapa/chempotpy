!> @brief Implement F2H2_GEN_JBKKL_1990, a potential energy surface for F2H2.
module chempotpy_f2h2_gen_jbkkl_1990_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the F2H2_GEN_JBKKL_1990 potential energy surface for F2H2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_f2h2_jbkkl_parameters, only: params_c => jbkkl_coefficients
        use chempotpy_units,                 only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds,                 only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 4
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: acoef_a(0:4,0:4,0:8)
        real(wp)              :: consts_bo
        real(wp)              :: consts_ha
        real(wp)              :: consts_pc
        real(wp)              :: consts_pi
        real(wp)              :: consts_vl
        integer               :: iatom, idir, istate
        real(wp)              :: cartnu(natoms,3), r(6), v
        real(wp)              :: rcm, r1, r2
        !initialize
        acoef_a = 0.e0_wp; consts_bo = 0.e0_wp; consts_ha = 0.e0_wp; consts_pc = 0.e0_wp; consts_pi = 0.e0_wp; consts_vl = 0.e0_wp
        v       = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        cartnu  = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                cartnu(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        call carttor(cartnu, r)
        rcm = r(1)
        r1  = r(2)
        r2  = r(3)
        call trbunk(rcm, r1, r2, v, 4, 0, 0)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = v*hartree_to_electronvolt
            enddo
        else
            write (*, *) 'Only energy is available'
        endif

        contains

        !> @brief Evaluate CARTTOR for the F2H2_GEN_JBKKL_1990 potential.
        subroutine carttor( cartnu, r )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: cmm, cphi, flm, hhd, hym, phi, pi, q, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
            real(wp) :: cartnu(4,3), r(6)
            flm    = 18.99840e0_wp
            hym    = 1.007825e0_wp
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
            if( theta1 .gt. 1.0e0_wp ) theta1=1.0e0_wp
            if( theta1 .lt. -1.0e0_wp ) theta1=-1.0e0_wp
            theta1 = acos(theta1)
            xrm2   = cartnu(3,1)-xcm2
            yrm2   = cartnu(3,2)-ycm2
            zrm2   = cartnu(3,3)-zcm2
            theta2 = (xrm2*(-xcm3)+yrm2*(-ycm3)+zrm2*(-zcm3))
            theta2 = theta2/(sqrt(xrm2**2+yrm2**2+zrm2**2))
            theta2 = theta2/(sqrt(xcm3**2+ycm3**2+zcm3**2))
            if( theta2 .gt. 1.0e0_wp ) theta2=1.0e0_wp
            if( theta2 .lt. -1.0e0_wp ) theta2=-1.0e0_wp
            theta2 = acos(theta2)
            pi     = acos(-1.0e0_wp)
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
            if( cphi .lt. -1.0e0_wp ) cphi=-1.0e0_wp
            if( cphi .gt. 1.0e0_wp ) cphi=1.0e0_wp
            phi  = acos(cphi)
            r(1) = sqrt(xcm3**2+ycm3**2+zcm3**2)
            r(2) = (sqrt(xrm1**2+yrm1**2+zrm1**2))*(flm+hym)/flm
            r(3) = (sqrt(xrm2**2+yrm2**2+zrm2**2))*(flm+hym)/flm
            r(4) = theta1
            r(5) = theta2
            r(6) = phi
        end subroutine carttor

        !> @brief Evaluate trbunk for the F2H2_GEN_JBKKL_1990 potential.
        subroutine trbunk( r, r1, r2, v, iq1, iq2, imu )
            !     subroutine trbunk transforms the surface of Bunker, Jensen, Kofran
            !     Lischka and Karpfen(the angular part only; no two-body terms) to t
            !     representation of Launay et al.
            !     references
            !     Bunker potential: P.R. Bunker, P. Jensen, A. Karpfen, M. Kofranek,
            !     H. Lischka, J. Chem. Phys. 92, 7432 (1990).
            !     Launay expansion: J.M. Launay, J. Phys. B 10, 3665 (1977).
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: cdenom, const, const1, const2, const3, phas, r, r1, r2, v, xl, xl1, xl2, xm
            integer             :: imu, iq1, iq2, l, ll, lm, lu
            integer, parameter  :: lmax    = 4
            integer, parameter  :: lsummax = 9
            integer, parameter  :: ipara   = 1
            real(wp), parameter :: zero    = 0.0e0_wp
            real(wp), parameter :: one     = 1.e0_wp
            real(wp), parameter :: two     = 2.e0_wp
            real(wp), parameter :: four    = 4.e0_wp
            integer, parameter  :: izero   = 0
            !     dimension a(0:4,0:4,0:8),v(ipara)
            !     dimension iq1(ipara),iq2(ipara),imu(ipara)
            consts_bo = 0.529177249e0_wp
            consts_ha = 219474.63067e0_wp
            consts_pi = 3.1415926535898e0_wp
            consts_vl = 2.99792458e0_wp
            consts_pc = 6.6260755e0_wp
            !     next get the a's from subroutine vcalc
            call vcalc(r1, r2, r)
            xl1 = real(iq1, wp)
            xl2 = real(iq2, wp)
            xm  = real(imu, wp)
            lm  = min(iq1, iq2)
            ll  = abs(iq1-iq2)
            lu  = iq1+iq2
            v   = 0.0e0_wp
            do l=ll, lu, 2
                xl     = real(l, wp)
                phas   = (-1)**(iq1-iq2)
                const1 = one
                if( imu .eq. 0 ) const1=two
                const1 = one/sqrt(const1)
                const2 = two*xl+one
                cdenom = (four*consts_pi)**3
                const3 = sqrt(two/cdenom)
                const  = phas*const1*const2*const3
                v      = v+acoef_a(iq1,iq2,l)*const*threej(xl1, xl2, xl, xm, -xm, 0.e0_wp)
            enddo
            !     divide next two terms by three because there is a bug in Bunker's
            !     function subprogram "P" which computes the associated Legendre
            !     function. (this is a very small correction).
            if( iq1 .eq. 4.and.iq2 .eq. 2.and.imu .eq. 2 ) v=v/3.e0_wp
            if( iq1 .eq. 2.and.iq2 .eq. 4.and.imu .eq. 2 ) v=v/3.e0_wp
            !     vpot=zero
            !     kount=kount-1
            !     do 50 j=1,kount
            !     write(65,*)v(j)
            !50   continue
        end subroutine trbunk

        !> @brief Evaluate THREEJ for the F2H2_GEN_JBKKL_1990 potential.
        real(wp) function threej( a, b, c, al, be, ga )
            !     CALCULATES WIGNER 3-J SYMBOLS USING FORMULA OF P34 OF BRINK &
            !     SATCHLER IN "ANGULAT MOMENTUM"
            !     NOTE FAKT(A)=FACT(A)/10**A
            !     A=J1,B=J2,C=J,AL=M1,BE=M2,GA=-M
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a, al, b, be, c, clebsh, delta, dnew, ga, summ, term, term1, term2, term3, term4, term5, term6
            integer  :: iphase, new, newmax, newmin
            threej = 0.e0_wp
            !     (J1+J2).GE.J AND J.GE.ABS(A-B)    -M=M1+M2    J1,J2,J.GE.0
            !     ABS(M1).LE.J1    ABS(M2).LE.J2   ABS(M).LE.J
            if( c .gt. a+b ) return
            if( c .lt. abs(a-b) ) return
            if( a .lt. 0.e0_wp.or.b .lt. 0.e0_wp.or.c .lt. 0.e0_wp ) return
            if( a .lt. abs(al).or.b .lt. abs(be).or.c .lt. abs(ga) ) return
            if( -1.e0_wp*ga .ne. al+be ) return
            !     COMPUTE DELTA(ABC)
            delta = sqrt(fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.e0_wp))
            term1 = fakt(a+al)*fakt(a-al)
            term2 = fakt(b-be)*fakt(b+be)
            term3 = fakt(c+ga)*fakt(c-ga)
            term  = sqrt((2.e0_wp*c+1.e0_wp)*term1*term2*term3)
            !     NOW COMPUTE SUMMATION TERM
            !     SUM TO GET SUMMATION IN EQ(2.34) OF BRINK AND SATCHLER.  SUM UNTIL
            !     A TERM INSIDE FACTORIAL GOES NEGATIVE.  NEW IS INDEX FOR SUMMATION
            !     .  NOW FIND WHAT THE RANGE OF NEW IS.
            newmin = idnint(max((a+be-c), (b-c-al), 0.e0_wp))
            newmax = idnint(min((a-al), (b+be), (a+b-c)))
            summ   = 0.e0_wp
            do new=newmin, newmax
                dnew  = dfloat(new)
                term4 = fakt(a-al-dnew)*fakt(c-b+al+dnew)
                term5 = fakt(b+be-dnew)*fakt(c-a-be+dnew)
                term6 = fakt(dnew)*fakt(a+b-c-dnew)
                summ  = summ+(-1.e0_wp)**new/(term4*term5*term6)
            enddo
            !     SO CLEBSCH-GORDEN <J1J2M1M2LJM> IS CLEBSH
            clebsh = delta*term*summ/sqrt(10.e0_wp)
            !     CONVERT CLEBSCH-GORDEN TO THREEJ
            iphase = idnint(a-b-ga)
            threej = (-1.e0_wp)**(iphase)*clebsh/sqrt(2.e0_wp*c+1.e0_wp)
        end function threej

        !> @brief Evaluate VCALC for the F2H2_GEN_JBKKL_1990 potential.
        subroutine vcalc( r1, r2, rab )
            !  CALCULATES V AT POINT (R1,R2,RAB,THETA1,THETA2,TAU)
            !  USING THE CURRENT VALUES OF THE PARAMETERS.
            !    V CALCULATED USING GREEN'S EXPRESSION; SEE JCP 62,2271 (1975).
            !       G IS THE SCRIPT I USED BY GREEN IN HIS EQ.(A9)
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: cd, cdq, cq, d1, d2, d3, de1, de2, de3, de4, g(0:4,0:4,0:8), r1, r2, r6, rab, req(6), rin01, rin10, &
                & rmiddle, rout10, t1, t1deg, t2, t2deg, tau, taudeg, y1, y1py2, y1y2, y1y22, y1y23, y1y24, y2, y3
            integer :: i, il, iu, j, k
            integer :: out
            out = 6
            !     R1=REQ(1)
            !     R2=REQ(2)
            !     RAB=REQ(3)
            !     T1DEG=REQ(4)
            !     T2DEG=REQ(5)
            !     TAUDEG=REQ(6)
            !     T1=T1DEG*PI/180.D0
            !     T2=T2DEG*PI/180.D0
            !     TAU=TAUDEG*PI/180.D0
            !     NEXT CALL SUBROUTINE HF2INT WHICH CALCULATES THE RELEVANT ATOM-ATO
            !     DISTANCES NEEDED TO COMPUTE THE TWO-BODY TERMS. THE ANGLES MUST
            !     BE IN DEGREES SO CONVERT THEM HERE.                           WCN
            t1deg  = t1*180.e0_wp/consts_pi
            t2deg  = t2*180.e0_wp/consts_pi
            taudeg = tau*180.e0_wp/consts_pi
            !     CALL HF2INT(R1,RAB,R2,T1DEG,T2DEG,TAUDEG,RHH,RHF1,RHF2,RFF)
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
            do i=0, 4
                do j=0, 4
                    do k=0, 8
                        acoef_a(i,j,k) = 0.0e0_wp
                    enddo
                enddo
            enddo
            d1  = r1-params_c(1)
            d2  = r2-params_c(1)
            d3  = rab-params_c(2)
            y1  = 1.e0_wp-exp(-params_c(3)*d1)
            y2  = 1.e0_wp-exp(-params_c(3)*d2)
            y3  = 1.e0_wp-exp(-d3)
            de1 = exp(-d3)
            de2 = exp(-2.0e0_wp*d3)
            de3 = exp(-3.0e0_wp*d3)
            de4 = exp(-4.0e0_wp*d3)
            r6  = rab**6
            !  DAMPING FACTORS AROUND RIN AND ROUT
            rout10  = 0.5e0_wp*(1.e0_wp-tanh(0.5e0_wp*(rab-7.5e0_wp)))
            rin10   = 0.5e0_wp*(1.e0_wp-tanh(2.e0_wp*(rab-4.0e0_wp)))
            rin01   = 0.5e0_wp*(1.e0_wp-tanh(2.e0_wp*(4.0e0_wp-rab)))
            rmiddle = rin01*rout10
            y1y22   = y1**2+y2**2
            y1y23   = y1**3+y2**3
            y1y24   = y1**4+y2**4
            y1py2   = y1+y2
            y1y2    = y1*y2
            !    CD, CQ AND CMQ ARE CONVERSION FACTORS FOR MULTIPOLE MOMENTS
            cd  = (consts_bo**3)*consts_pc*consts_vl*consts_ha/100000.e0_wp
            cq  = cd*consts_bo*consts_bo
            cdq = cd*consts_bo
            !    C(4) TO C(9), C(91) AND C(111) DESCRIBE THE DISPERSION INTERACTION
            !    C(9) IS IN A(2,2,4).
            !    C(41) IS THE HF DIPOLE MOMENT IN DEBYE
            !    C(101) IS THE HF QUADRUPOLE MOMENT IN DEBYE.ANGSTROM
            !    C(131) IS THE HF OCTUPOLE MOMENT IN DEBYE.ANGSTROM**2
            !     the zero of energy is chosen as the two monomers infinitely separa
            !     at r1=r2=1.744 a.u. thus c(11) has been omitted from a(0,0,0).
            acoef_a(0,0,0)=params_c(357)*de2*rin10 +params_c(20)*y1y2*rout10 +(params_c(12)*de1+ params_c(13)*de2+ &
                & params_c(14)*de3+params_c(15)*de4)*rmiddle +(params_c(4)/r6 + params_c(5)/ (rab**7) + params_c(6)/ &
                & (rab**8) +params_c(7)/(rab**9)+params_c(8)/(rab** 10))*rin01 +(params_c(157)*(y1+y2)* y3 + &
                & params_c(158)*y1*y2* y3 + params_c(166)* y1* y2*y3**2 +params_c(167)*y1*y2*y3**3)*rmiddle
            acoef_a(1,0,1)=(params_c(22)*de1+params_c(23)*de2 +params_c(24)*de3+params_c(25)*de4 + params_c(26)* y1 + &
                & params_c(27)*y2 +params_c(21)*y1*y3 +params_c(156)*y2*y3 + params_c(28)*y1** 2 + params_c(29)*y2** &
                & 2 + params_c(30)*y1*y2 +params_c(116)*y1**3 + params_c(117)*y2**3 +params_c(118)* y1*y2**2 + &
                & params_c(119)*y2*y1**2 +params_c(159)*y1* y2*y3 + params_c(169)*y1*y2*y3**2 + params_c(179)*y1*y2** &
                & 2*y3)*rmiddle
            acoef_a(0,1,1)=-(params_c(22)*de1+params_c(23)*de2 +params_c(24)*de3+params_c(25)*de4 + params_c(26)* y2 + &
                & params_c(27)*y1 +params_c(21)*y2*y3 +params_c(156)*y1*y3 + params_c(28)*y2** 2 + params_c(29)*y1** &
                & 2 + params_c(30)*y1*y2 +params_c(116)*y2**3 + params_c(117)*y1**3 +params_c(118)* y2*y1**2 + &
                & params_c(119)*y1*y2**2 +params_c(159)*y1* y2* y3 + params_c(169)*y1*y2*y3**2 + params_c(179)*y1**2* &
                & y1*y3)*rmiddle
            acoef_a(1,1,0)=(params_c(32)*de1+params_c(33)*de2 +params_c(34)*de3+params_c(35)*de4 + params_c(36)* &
                & y1y22 +params_c(37)*y1y23 +params_c(38)*y1py2 +params_c(39)*y1py2*y1y2 + params_c(40)* y1y2)* rmiddle
            acoef_a(1,1,2)=(params_c(42)*de1+params_c(43)*de2 +params_c(44)*de3+params_c(45)*de4)* rmiddle -((8.e0_wp* &
                & consts_pi**1.5e0_wp*sqrt(2.e0_wp/15.e0_wp))*(params_c(41)**2)/ (cd* rab**3))* rin01 +(params_c(46)* &
                & y1y22 + params_c(47)*y1y23 +params_c(48)*y1py2 +params_c(49)* y1py2*y1y2 +params_c(50)*y1y2)*rmiddle
            acoef_a(2,0,2)=params_c(356)/r6+(params_c(52)*de1+params_c(53)*de2 +params_c(54)*de3+ params_c(55)* de4 + &
                & params_c(56)*y1 +params_c(57)*y2 +params_c(58)*y1**2 +params_c(59)* y2**2 + params_c(60)* y1*y2 + &
                & params_c(126)*y1**3 +params_c(127)*y2**3 +params_c(128)* y1*y2**2 +params_c(129)*y2*y1** 2)*rmiddle
            acoef_a(0,2,2)=params_c(356)/r6+(params_c(52)*de1+params_c(53)*de2 +params_c(54)*de3+ params_c(55)* de4 + &
                & params_c(56)*y2 +params_c(57)*y1 +params_c(58)*y2**2 +params_c(59)* y1**2 + params_c(60)* y1*y2 + &
                & params_c(126)*y2**3 +params_c(127)*y1**3 +params_c(128)* y2*y1**2 +params_c(129)*y1*y2** 2)*rmiddle
            acoef_a(2,1,1)=(params_c(62)*de1+params_c(63)*de2 +params_c(64)*de3+params_c(65)*de4 + params_c(66)* y1 + &
                & params_c(67)*y2 +params_c(68)*y1**2 +params_c(69)*y2**2 + params_c(70)*y1* y2 + params_c(136)*y1** &
                & 3 + params_c(137)*y2**3 +params_c(138)*y1*y2** 2 +params_c(139)*y2*y1**2)*rmiddle
            acoef_a(1,2,1)=-(params_c(62)*de1+params_c(63)*de2 +params_c(64)*de3+params_c(65)*de4 + params_c(66)* y2 + &
                & params_c(67)*y1 +params_c(68)*y2**2 +params_c(69)*y1**2 + params_c(70)*y1* y2 + params_c(136)*y2** &
                & 3 + params_c(137)*y1**3 +params_c(138)*y2*y1** 2 +params_c(139)*y1*y2**2)*rmiddle
            acoef_a(2,1,3)=(params_c(72)*de1+params_c(73)*de2 +params_c(74)*de3+params_c(75)*de4)* rmiddle -rin01* &
                & (8.e0_wp*consts_pi**1.5e0_wp/sqrt(7.e0_wp))*params_c(41)*params_c(101)/ (cdq* rab**4) + &
                & (params_c(76)* y1 + params_c(77)*y2 +params_c(78)*y1**2 +params_c(79)*y2**2 + params_c(80)*y1*y2 + &
                & params_c(146)* y1**3 + params_c(147)*y2**3 + params_c(148)*y1* y2** 2 +params_c(149)*y2*y1**2)*rmiddle
            acoef_a(1,2,3)=-(params_c(72)*de1+params_c(73)*de2 +params_c(74)*de3+params_c(75)*de4)* rmiddle +rin01* &
                & (8.e0_wp*consts_pi**1.5e0_wp/sqrt(7.e0_wp))*params_c(41)*params_c(101)/ (cdq* rab**4) - &
                & (params_c(76)* y2 + params_c(77)*y1 +params_c(78)*y2**2 +params_c(79)*y1**2 + params_c(80)*y1*y2 + &
                & params_c(146)* y2**3 + params_c(147)*y1**3 + params_c(148)*y2* y1** 2 +params_c(149)*y1*y2**2)*rmiddle
            acoef_a(2,2,0)=(params_c(82)*de1+params_c(83)*de2 +params_c(84)*de3+params_c(85)*de4 + params_c(86)* &
                & y1y22 +params_c(87)*y1y23 +params_c(88)*y1py2 +params_c(89)*y1py2*y1y2 + params_c(90)* y1y2)* rmiddle
            acoef_a(2,2,2)=params_c(91)/r6 +(params_c(92)*de1+params_c(93)*de2 +params_c(94)*de3+ params_c(95)* de4 + &
                & params_c(96)*y1y22 +params_c(97)*y1y23 +params_c(98)*y1py2 + params_c(99)* y1py2* y1y2 + &
                & params_c(100)*y1y2)*rmiddle
            acoef_a(2,2,4)=(params_c(102)*de1+params_c(103)*de2 +params_c(104)*de3+params_c(105)*de4)* rmiddle + &
                & (8.e0_wp* consts_pi**1.5e0_wp/3.e0_wp)*sqrt(14.e0_wp/5.e0_wp)* rin01* (params_c(101)**2)/ (cq*rab** &
                & 5) + (params_c(106)* y1y22 +params_c(107)*y1y23 +params_c(108)*y1py2 + params_c(109)*y1py2*y1y2 + &
                & params_c(110)* y1y2)* rmiddle
            acoef_a(3,0,3) = params_c(111)/(rab**7) +(params_c(112)*de1+params_c(113)*de2 +params_c(114)* de3+ &
                & params_c(115)*de4)*rmiddle
            acoef_a(0,3,3) = -acoef_a(3,0,3)
            acoef_a(3,1,2) = (params_c(122)*de1+params_c(123)*de2 +params_c(124)*de3+params_c(125)*de4)* rmiddle
            acoef_a(1,3,2) = acoef_a(3,1,2)
            acoef_a(3,1,4)=(params_c(132)*de1+params_c(133)*de2 +params_c(134)*de3+params_c(135)*de4)* rmiddle - &
                & rin01* (16.e0_wp*((consts_pi/3.e0_wp)**1.5e0_wp))*params_c(41)* params_c(131)/ (cq* rab**5)
            acoef_a(1,3,4) = acoef_a(3,1,4)
            acoef_a(3,2,1) = (params_c(142)*de1+params_c(143)*de2 +params_c(144)*de3+params_c(145)*de4)* rmiddle
            acoef_a(2,3,1) = -acoef_a(3,2,1)
            acoef_a(3,2,3) = (params_c(152)*de1+params_c(153)*de2 +params_c(154)*de3+params_c(155)*de4)* rmiddle
            acoef_a(2,3,3) = -acoef_a(3,2,3)
            acoef_a(3,2,5) = (params_c(162)*exp(-0.5e0_wp*d3)+params_c(163)*exp(-0.75e0_wp*d3) +params_c(164)* de1+ &
                & params_c(165)*de4)*rmiddle
            acoef_a(2,3,5) = -acoef_a(3,2,5)
            acoef_a(3,3,0) = (params_c(172)*de1+params_c(173)*de2 +params_c(174)*de3+params_c(175)*de4 + &
                & params_c(176)*y1y22+params_c(177)*y1y2)*rmiddle
            acoef_a(3,3,2) = (params_c(182)*de1+params_c(183)*de2 +params_c(184)*de3+params_c(185)*de4 + &
                & params_c(186)*y1y22+params_c(187)*y1y2)*rmiddle
            acoef_a(3,3,4) = (params_c(192)*de1+params_c(193)*de2 +params_c(194)*de3+params_c(195)*de4 + &
                & params_c(196)*y1y22+params_c(197)*y1y2)*rmiddle
            acoef_a(3,3,6) = (params_c(202)*de1+params_c(203)*de2 +params_c(204)*de3+params_c(205)*de4 + &
                & params_c(206)*y1y22+params_c(207)*y1y2)*rmiddle
            acoef_a(4,0,4) = params_c(212)*de1+params_c(213)*de2 +params_c(214)*de3+params_c(215)*de4
            acoef_a(0,4,4) = acoef_a(4,0,4)
            acoef_a(4,1,3) = params_c(222)*de1+params_c(223)*de2 +params_c(224)*de3+params_c(225)*de4
            acoef_a(1,4,3) = -acoef_a(4,1,3)
            acoef_a(4,1,5) = params_c(232)*exp(-0.5e0_wp*d3)+params_c(233)*exp(-0.75e0_wp*d3) +params_c(234)* de1+ &
                & params_c(235)*de2
            acoef_a(1,4,5) = -acoef_a(4,1,5)
            acoef_a(4,2,2) = params_c(242)*de1+params_c(243)*de2 +params_c(244)*de3+params_c(245)*de4
            acoef_a(2,4,2) = acoef_a(4,2,2)
            acoef_a(4,2,4) = params_c(252)*de1+params_c(253)*de2 +params_c(254)*de3+params_c(255)*de4
            acoef_a(2,4,4) = acoef_a(4,2,4)
            acoef_a(4,2,6) = params_c(262)*exp(-0.5e0_wp*d3)+params_c(263)*exp(-0.75e0_wp*d3) +params_c(264)* de1+ &
                & params_c(265)*de2
            acoef_a(2,4,6) = acoef_a(4,2,6)
            acoef_a(4,3,1) = params_c(272)*de1+params_c(273)*de2 +params_c(274)*de3+params_c(275)*de4
            acoef_a(3,4,1) = -acoef_a(4,3,1)
            acoef_a(4,3,3) = params_c(282)*de1+params_c(283)*de2 +params_c(284)*de3+params_c(285)*de4
            acoef_a(3,4,3) = -acoef_a(4,3,3)
            acoef_a(4,3,5) = params_c(292)*de1+params_c(293)*de2 +params_c(294)*de3+params_c(295)*de4
            acoef_a(3,4,5) = -acoef_a(4,3,5)
            acoef_a(4,3,7) = params_c(302)*de1+params_c(303)*de2 +params_c(304)*de3+params_c(305)*de4
            acoef_a(3,4,7) = -acoef_a(4,3,7)
            acoef_a(4,4,0) = params_c(312)*de1+params_c(313)*de2 +params_c(314)*de3+params_c(315)*de4 + &
                & (params_c(316)*y1y22+params_c(317)*y1y2)
            acoef_a(4,4,2) = params_c(322)*de1+params_c(323)*de2 +params_c(324)*de3+params_c(325)*de4
            acoef_a(4,4,4) = params_c(332)*de1+params_c(333)*de2 +params_c(334)*de3+params_c(335)*de4
            acoef_a(4,4,6) = params_c(342)*de1+params_c(343)*de2 +params_c(344)*de3+params_c(345)*de4
            acoef_a(4,4,8) = params_c(352)*de1+params_c(353)*de2 +params_c(354)*de3+params_c(355)*de4
            do i=0, 4
                do j=0, 4
                    il = iabs(i-j)
                    iu = i+j
                    do k=il, iu, 2
                        if( i .ne. 4.and.j .ne. 4 ) cycle
                        acoef_a(i,j,k) = acoef_a(i,j,k)*rmiddle
                    enddo
                enddo
            enddo
            !     V = 0.0D0
            !     DO 2 I=0,4
            !     DO 2 J=0,4
            !     IL=IABS(I-J)
            !     IU=I+J
            !     DO 2 K=IL,IU,2
            !2    V = V + A(I,J,K)*G(I,J,K)
            !  HH, FF, HF1 AND HF2 2-BODY TERMS COMBINED IN TERM VSRG.
            !  VSRG IS DAMPED TO ZERO IF RHF.GT.3.0
            !  IN ORDER TO PREVENT SPURIOUS 'HOLES' IN V.
            !  DAMPED TO ZERO IF RHH.GT.2.3, RHF1/2.GT.2.3, AND RFF.GT.2.3
            !  SO THAT BOUND STATES CAN BE CALCULATED WITH THE NEGLECT
            !  OF VSRG AND CLOSE COUPLING PROGRAM 'BOUND' USED WITH V.
            !  eliminate two-body terms completely in this version wcn 3/10/92
            !     YRHF=0.25D0*(1.D0-DTANH(7.D0*(R1-3.0D0)))
            !    X           *(1.D0-DTANH(7.D0*(R2-3.0D0)))
            !     YRHH=0.50D0*(1.D0-DTANH(7.D0*(RHH-2.3D0)))
            !     YRFF=0.50D0*(1.D0-DTANH(7.D0*(RFF-2.3D0)))
            !     YRHF1=0.50D0*(1.D0-DTANH(7.D0*(RHF1-2.3D0)))
            !     YRHF2=0.50D0*(1.D0-DTANH(7.D0*(RHF2-2.3D0)))
            !     VHH=C(161)*DEXP(-C(171)*RHH)+C(181)*DEXP(-C(191)*RHH)
            !     VHF1=C(201)*DEXP(-C(211)*RHF1)+C(221)*DEXP(-C(231)*RHF1)
            !     VHF2=C(201)*DEXP(-C(211)*RHF2)+C(221)*DEXP(-C(231)*RHF2)
            !     VFF=C(241)*DEXP(-C(251)*RFF)+C(261)*DEXP(-C(271)*RFF)
            !     VSRG=(VHH*YRHH+VHF1*YRHF1+VHF2*YRHF2+VFF*YRFF)*YRHF
            !     V=V+VSRG
        end subroutine vcalc

        !> @brief Evaluate FAKT for the F2H2_GEN_JBKKL_1990 potential.
        real(wp) function fakt( a )
            !     THIS FUNCTION COMPUTES N!/(10**N).  STOP 101 IMPLIES TRIED
            !     TO TAKE FACTORIAL OF A NEGATIVE NO.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a, ax
            integer  :: i, ic
            ax   = a
            fakt = 1.e0_wp
            if( ax .eq. 0.e0_wp ) return
            fakt = .1e0_wp
            if( ax .lt. 0.e0_wp ) stop 101
            ic   = idnint(ax)
            ax   = ax/10.e0_wp
            fakt = ax
            do i=1, ic-1
                fakt = fakt*(ax-dfloat(i)*.1e0_wp)
            enddo
        end function fakt
    end subroutine evaluate
end module chempotpy_f2h2_gen_jbkkl_1990_pes

!> @brief Expose F2H2_GEN_JBKKL_1990 through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_f2h2_gen_jbkkl_1990_pes, only: evaluate
    use chempotpy_kinds,                   only: wp
    implicit none
    integer, parameter              :: natoms  = 4
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
