module zmasses
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: zmass1, zmass2
    real(wp) :: zmass1, zmass2
end module

!> @brief Evaluate pot.
subroutine pot( rvp, evp, dvp )
    use zmasses,         only: zmass1, zmass2
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), parameter   :: ccishift = 0.1744759989649372e0_wp
    real(wp), parameter   :: towave   = 219474.7e0_wp
    real(wp), parameter   :: xmass1   = 1836.15264e0_wp, xmass2 = 1836.15264e0_wp, xmass3 = 1836.15264e0_wp
    real(wp), intent(in)  :: rvp(3)
    real(wp), intent(out) :: evp, dvp(3)
    integer               :: icase, info, i
    real(wp)              :: r1x, r2x, r3x, ecci
    real(wp)              :: dcci(3), d1(3), d2(3), d3(3)
    real(wp)              :: xmh, xmd, xmminv(3,3), amat(3), t
    real(wp)              :: e1, e2, e3, emat(3), ezero, hatom_bodc, xmin
    real(wp)              :: ezero2
    !     write(6,*)present(newmasses)
    !     if(present(newmasses))then
    !       xmass1=newmasses(1)
    !       xmass2=newmasses(2)
    !       xmass3=newmasses(3)
    !     endif
    !    We need to zero the PES with respect to one arrangement
    !    If only one mass is distinct, choose the target mass as that, otherwise choose
    !    the lightest particle.  If another choice is desired, add appropriate assignment
    !    statements here
    if( xmass1 .eq. xmass2 )then
        zmass1 = xmass1
        zmass2 = xmass2
    elseif( xmass1 .eq. xmass3 )then
        zmass1 = xmass1
        zmass2 = xmass3
    elseif( xmass2 .eq. xmass3 )then
        zmass1 = xmass2
        zmass2 = xmass3
    else
        zmass1 = max(xmass1, xmass2, xmass3)
        if( xmass1 .lt. zmass1 )then
            zmass2 = max(xmass1, min(xmass2, xmass3))
        else
            zmass2 = max(xmass2, xmass3)
        endif
    endif
    xmh        = 1836.15264e0_wp   ! H nuclear mass (after substracting 1 for electron mass)
    xmd        = 3670.48293e0_wp   ! D nuclear mass (after substracting 1 for electron mass)
    hatom_bodc = 59.6818e0_wp  ! the BODC in cm-1 of H with an aug-cc-pVTZ basis set
    ezero      = hatom_bodc*xmh*(1e0_wp/xmass1+1e0_wp/xmass2+1e0_wp/xmass3)
    call golden(1.3e0_wp, 1.4e0_wp, 1.5e0_wp, 1.e-14_wp, xmin, ezero2)
    ezero = ezero/towave+ezero2+ccishift
    if( (xmass1.eq.xmd) .and. (xmass2.eq.xmh) &
        &      .and. (xmass3.eq.xmh) )then
    icase = 1  ! We are solving for D + H2 so we don't need to solve for amat
else
    icase       = -1 ! We are not solving for D + H2 so we need 3 times as many evaluations so we can determine amat
    t           = (2.0e0_wp*xmd+xmh)*(xmd-xmh)
    xmminv(:,:) = xmh*xmd*xmd/t
    t           = -xmh*xmd*(xmd+xmh)/t
    xmminv(1,1) = t
    xmminv(2,2) = t
    xmminv(3,3) = t
    endif
    ! R1 is the distance between mass1 and mass2, R2 is the distance between mass2 and mass3,
    ! and R3 is the distance between mass1 and mass3.  Note that the D+H2 code requires the homonuclear
    ! diatom distance to be R3, so before calling all the various subroutines we need to switch R2 and R3
    r1x = rvp(1)
    r3x = rvp(2)
    r2x = rvp(3)
    call potcbs(r1x, r2x, r3x, ecci, dcci)
    call potbodcdh2(r1x, r2x, r3x, e1, d1)
    if( icase .eq. 1 )then
        evp    = e1 - ezero
        dvp(:) = d1(:)
    else
        call potbodcdh2(r1x, r3x, r2x, e2, d2)
        t     = d2(2)
        d2(2) = d2(3)
        d2(3) = t
        call potbodcdh2(r3x, r2x, r1x, e3, d3)
        t       = d3(1)
        d3(1)   = d3(3)
        d3(3)   = t
        emat(1) = e1
        emat(2) = e2
        emat(3) = e3
        amat(:) = matmul(xmminv, emat)
        evp     = amat(1)/xmass1+amat(2)/xmass2+amat(3)/xmass3 - ezero
        do i =1, 3
            emat(1) = d1(i)
            emat(2) = d2(i)
            emat(3) = d3(i)
            amat(:) = matmul(xmminv, emat)
            dvp(i)  = amat(1)/xmass1+amat(2)/xmass2+amat(3)/xmass3
        enddo
    endif
    evp    = evp/towave+ecci-ezero
    dvp(1) = dvp(1)/towave+dcci(1)
    t      = dvp(2)/towave+dcci(2)
    dvp(2) = dvp(3)/towave+dcci(3)
    dvp(3) = t
end

module chempotpy_h3_gen_bh_2009_potbodcdh2_legacy_state_b8d1df42
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r1x, r2x, r3x, eee, dv )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r1x, r2x, r3x, eee, dv
        integer, parameter :: n2 = 5000
        integer            :: iind(n2), jind(n2), kind(n2), ipar(n2)
        integer            :: ii, it, maxi, maxi2, maxi2a, i
        real(wp)           :: r1x, r2x, r3x, rho1, rho2, rho3, eee, dv(:)
        real(wp)           :: drho1, drho2, drho3
        real(wp)           :: e1body, xmd, xmh, hatom_bodc, bodc2body
        real(wp)           :: round, ald, betad, beta, deldamp, rhodamp
        real(wp)           :: bb2, bb3, bb4, bb5, v, rnought, towave
        real(wp)           :: xma, xmb, xmc, v1, v2, gtot, rsum, &
            &                    betadsum, rdampsum, v3c
        real(wp) :: sum1, sum2, sum3, sum4, sum5, rprod
        real(wp) :: bigq, slit, slitsq, damper, rpass(3), &
            &                    sdamp, rdamp1, rdamp2, rdamp3
        real(wp) :: slit1, slit2, slit3, ddamper(1:3), dsdamp
        real(wp) :: dgtot(1:3), dbb2(1:3), dbb3(1:3)
        real(wp) :: dv1, dv2, dv3
        real(wp) :: dsum1(1:3), dsum2(1:3), dsum3(1:3)
        real(wp) :: rho1t(0:13), rho2t(0:13), rho3t(0:13)
        real(wp) :: drho1t(0:13), drho2t(0:13), drho3t(0:13)
        real(wp) :: drprod1, drprod2, drprod3
        real(wp) :: xlin1(180), xlin2(180), xlin3(180)
        ald        = 0.e0_wp
        bb2        = 0.e0_wp
        bb3        = 0.e0_wp
        bb4        = 0.e0_wp
        bb5        = 0.e0_wp
        beta       = 0.e0_wp
        betad      = 0.e0_wp
        betadsum   = 0.e0_wp
        bigq       = 0.e0_wp
        damper     = 0.e0_wp
        dbb2       = 0.e0_wp
        dbb3       = 0.e0_wp
        ddamper    = 0.e0_wp
        deldamp    = 0.e0_wp
        dgtot      = 0.e0_wp
        drho1      = 0.e0_wp
        drho1t     = 0.e0_wp
        drho2      = 0.e0_wp
        drho2t     = 0.e0_wp
        drho3      = 0.e0_wp
        drho3t     = 0.e0_wp
        drprod1    = 0.e0_wp
        drprod2    = 0.e0_wp
        drprod3    = 0.e0_wp
        dsdamp     = 0.e0_wp
        dsum1      = 0.e0_wp
        dsum2      = 0.e0_wp
        dsum3      = 0.e0_wp
        dv1        = 0.e0_wp
        dv2        = 0.e0_wp
        dv3        = 0.e0_wp
        e1body     = 0.e0_wp
        gtot       = 0.e0_wp
        hatom_bodc = 0.e0_wp
        rdamp1     = 0.e0_wp
        rdamp2     = 0.e0_wp
        rdamp3     = 0.e0_wp
        rdampsum   = 0.e0_wp
        rho1       = 0.e0_wp
        rho1t      = 0.e0_wp
        rho2       = 0.e0_wp
        rho2t      = 0.e0_wp
        rho3       = 0.e0_wp
        rho3t      = 0.e0_wp
        rhodamp    = 0.e0_wp
        rnought    = 0.e0_wp
        round      = 0.e0_wp
        rpass      = 0.e0_wp
        rprod      = 0.e0_wp
        rsum       = 0.e0_wp
        sdamp      = 0.e0_wp
        slit       = 0.e0_wp
        slit1      = 0.e0_wp
        slit2      = 0.e0_wp
        slit3      = 0.e0_wp
        slitsq     = 0.e0_wp
        sum1       = 0.e0_wp
        sum2       = 0.e0_wp
        sum3       = 0.e0_wp
        sum4       = 0.e0_wp
        sum5       = 0.e0_wp
        towave     = 0.e0_wp
        v          = 0.e0_wp
        v1         = 0.e0_wp
        v2         = 0.e0_wp
        v3c        = 0.e0_wp
        xlin1      = 0.e0_wp
        xlin2      = 0.e0_wp
        xlin3      = 0.e0_wp
        xma        = 0.e0_wp
        xmb        = 0.e0_wp
        xmc        = 0.e0_wp
        xmd        = 0.e0_wp
        xmh        = 0.e0_wp
        xlin1(1)   = 2.451839720665176e+02_wp
        xlin1(2)   = -3.023514674602468e+02_wp
        xlin1(3)   = -3.643179260300581e+02_wp
        xlin1(4)   = 1.033283164173005e+03_wp
        xlin1(5)   = -9.649202796347731e+02_wp
        xlin1(6)   = 4.723249503697386e+02_wp
        xlin1(7)   = -1.312024315948558e+02_wp
        xlin1(8)   = 2.053242829600808e+01_wp
        xlin1(9)   = -1.663116010000327e+00_wp
        xlin1(10)  = 5.317955346689944e-02_wp
        xlin1(11)  = 3.693835631800750e+02_wp
        xlin1(12)  = 9.643897029819937e+01_wp
        xlin1(13)  = -3.794338103064540e+02_wp
        xlin1(14)  = 2.184383291329870e+02_wp
        xlin1(15)  = -3.787026565563110e+01_wp
        xlin1(16)  = -3.960932914729867e+00_wp
        xlin1(17)  = 1.733261249407424e+00_wp
        xlin1(18)  = -1.233088325428533e-01_wp
        xlin1(19)  = 7.634113468160084e+02_wp
        xlin1(20)  = -6.089390321694615e+02_wp
        xlin1(21)  = 2.145121698941533e+02_wp
        xlin1(22)  = -4.303659182480859e+01_wp
        xlin1(23)  = 4.793100434454963e+00_wp
        xlin1(24)  = -2.362075732031712e-01_wp
        xlin1(25)  = 5.161649762181137e+02_wp
        xlin1(26)  = -1.379218427454547e+02_wp
        xlin1(27)  = 1.222214446630326e+01_wp
        xlin1(28)  = -1.700714347704069e-01_wp
        xlin1(29)  = 2.590291319801722e+01_wp
        xlin1(30)  = -1.187599837006957e+00_wp
        xlin1(31)  = 2.949665817120014e+02_wp
        xlin1(32)  = -4.423967394072385e+02_wp
        xlin1(33)  = -1.499691959357717e+02_wp
        xlin1(34)  = 8.615959233030945e+02_wp
        xlin1(35)  = -9.063992109731515e+02_wp
        xlin1(36)  = 4.800278837250469e+02_wp
        xlin1(37)  = -1.444131347943261e+02_wp
        xlin1(38)  = 2.489091872075205e+01_wp
        xlin1(39)  = -2.288354420765698e+00_wp
        xlin1(40)  = 8.706064188372947e-02_wp
        xlin1(41)  = -1.561653042358118e+03_wp
        xlin1(42)  = 2.795630390969673e+03_wp
        xlin1(43)  = -4.442497231862677e+03_wp
        xlin1(44)  = 2.653340911293203e+03_wp
        xlin1(45)  = -7.274300003583212e+02_wp
        xlin1(46)  = 1.223835539827257e+02_wp
        xlin1(47)  = -1.767422177777780e+01_wp
        xlin1(48)  = 1.080509007725006e+00_wp
        xlin1(49)  = 8.007412890755861e-02_wp
        xlin1(50)  = -6.866852114392939e+02_wp
        xlin1(51)  = 1.965474846931427e+03_wp
        xlin1(52)  = -8.327550428277026e+02_wp
        xlin1(53)  = 2.921969416962440e+00_wp
        xlin1(54)  = 5.934336236873139e+01_wp
        xlin1(55)  = -1.157242040507434e+01_wp
        xlin1(56)  = 6.777877552427620e-01_wp
        xlin1(57)  = 6.809701703798451e+02_wp
        xlin1(58)  = -6.460103807549953e+02_wp
        xlin1(59)  = 1.080711387406044e+02_wp
        xlin1(60)  = -2.581355973993933e+00_wp
        xlin1(61)  = -2.702457089997727e-01_wp
        xlin1(62)  = 3.373040929271790e+02_wp
        xlin1(63)  = -3.319320528532528e+01_wp
        xlin1(64)  = -8.252443297418132e-01_wp
        xlin1(65)  = 2.242813487013438e+00_wp
        xlin1(66)  = -4.786033544785251e+02_wp
        xlin1(67)  = 7.379465395733814e+02_wp
        xlin1(68)  = -2.997606193372937e+02_wp
        xlin1(69)  = -1.444004106661323e+02_wp
        xlin1(70)  = 1.160834932238466e+02_wp
        xlin1(71)  = -1.055041157003008e+01_wp
        xlin1(72)  = -6.084935176995607e+00_wp
        xlin1(73)  = 1.255348914553237e+00_wp
        xlin1(74)  = -4.967414449352265e-02_wp
        xlin1(75)  = 2.869805514805477e+03_wp
        xlin1(76)  = -2.901158537312932e+02_wp
        xlin1(77)  = 2.593678864031728e+03_wp
        xlin1(78)  = -1.474577262776355e+03_wp
        xlin1(79)  = 1.739661050571814e+02_wp
        xlin1(80)  = 4.108033331334562e+01_wp
        xlin1(81)  = -9.031479084332030e+00_wp
        xlin1(82)  = 3.265756035116496e-01_wp
        xlin1(83)  = -9.394389739950026e+03_wp
        xlin1(84)  = 1.281273536289386e+03_wp
        xlin1(85)  = 9.186863730169819e+02_wp
        xlin1(86)  = -2.224206878555418e+02_wp
        xlin1(87)  = 2.852507226985430e+01_wp
        xlin1(88)  = -2.975868406356395e+00_wp
        xlin1(89)  = -2.887753159216629e+02_wp
        xlin1(90)  = -6.232442825355130e+01_wp
        xlin1(91)  = -2.080038525889382e+01_wp
        xlin1(92)  = 5.846990317948459e+00_wp
        xlin1(93)  = 1.494611359198648e+01_wp
        xlin1(94)  = 1.430793238496208e+00_wp
        xlin1(95)  = -1.161743802353405e+02_wp
        xlin1(96)  = -4.435211957384280e+02_wp
        xlin1(97)  = 7.949179810694009e+02_wp
        xlin1(98)  = -4.621878870840682e+02_wp
        xlin1(99)  = 1.829448335674127e+02_wp
        xlin1(100) = -4.859832146779939e+01_wp
        xlin1(101) = 6.467586280687250e+00_wp
        xlin1(102) = -3.014121429251370e-01_wp
        xlin1(103) = -4.727165211936207e+03_wp
        xlin1(104) = 1.681945178807448e+03_wp
        xlin1(105) = 6.001879815346814e+02_wp
        xlin1(106) = -5.284816509082527e+02_wp
        xlin1(107) = 9.509901063900256e+01_wp
        xlin1(108) = -8.355776044216967e+00_wp
        xlin1(109) = 4.667856406455721e-01_wp
        xlin1(110) = 5.648107714717939e+02_wp
        xlin1(111) = -2.902912088460810e+02_wp
        xlin1(112) = -9.407789952750028e+01_wp
        xlin1(113) = -1.631433254869286e+00_wp
        xlin1(114) = 5.242491229250340e+00_wp
        xlin1(115) = 3.214464002696617e+01_wp
        xlin1(116) = 2.913754508725327e+01_wp
        xlin1(117) = -1.568292873005413e+01_wp
        xlin1(118) = 5.462484623657479e+00_wp
        xlin1(119) = 9.083635945770520e+02_wp
        xlin1(120) = 9.799890386730289e+01_wp
        xlin1(121) = -6.411746299039472e+02_wp
        xlin1(122) = 4.298716219171475e+02_wp
        xlin1(123) = -1.205657322203302e+02_wp
        xlin1(124) = 1.315814978109771e+01_wp
        xlin1(125) = -3.975978827444066e-01_wp
        xlin1(126) = 3.180124260914723e+03_wp
        xlin1(127) = -5.872907359889355e+02_wp
        xlin1(128) = -6.259492461216178e+02_wp
        xlin1(129) = 3.434767513553659e+02_wp
        xlin1(130) = -3.640741128738960e+01_wp
        xlin1(131) = -9.236985549909328e-01_wp
        xlin1(132) = 1.155215128322806e+03_wp
        xlin1(133) = -9.784160598115349e+01_wp
        xlin1(134) = 1.777745091732398e+01_wp
        xlin1(135) = 2.485226250947354e-01_wp
        xlin1(136) = 1.521755613790718e+01_wp
        xlin1(137) = 7.156889276881836e+00_wp
        xlin1(138) = -1.016142747224836e+03_wp
        xlin1(139) = -5.200008131651686e+01_wp
        xlin1(140) = 2.946853566528521e+02_wp
        xlin1(141) = -1.316509239684804e+02_wp
        xlin1(142) = 2.366837195107599e+01_wp
        xlin1(143) = -1.174923797943324e+00_wp
        xlin1(144) = -1.159057124509925e+03_wp
        xlin1(145) = -1.521893512764090e+02_wp
        xlin1(146) = 1.104768477989274e+02_wp
        xlin1(147) = -3.656291398178829e+01_wp
        xlin1(148) = 3.488333934608939e+00_wp
        xlin1(149) = -2.286294041696426e+02_wp
        xlin1(150) = -7.953044862181258e+00_wp
        xlin1(151) = -2.999788001060466e-02_wp
        xlin1(152) = -1.546533722850765e+01_wp
        xlin1(153) = 5.670669717065498e+02_wp
        xlin1(154) = 5.041037222513963e+01_wp
        xlin1(155) = -7.690800748197464e+01_wp
        xlin1(156) = 1.401968534663515e+01_wp
        xlin1(157) = -1.189136309799305e+00_wp
        xlin1(158) = 2.944715172811239e+02_wp
        xlin1(159) = 1.011155288939896e+02_wp
        xlin1(160) = -5.420733436869213e+00_wp
        xlin1(161) = -8.900491551445843e-01_wp
        xlin1(162) = 2.785072259806540e+01_wp
        xlin1(163) = 6.092864636803840e+00_wp
        xlin1(164) = -1.799957321893311e+02_wp
        xlin1(165) = -1.892894466520360e+01_wp
        xlin1(166) = 1.036160665441249e+01_wp
        xlin1(167) = -3.761882545122726e-01_wp
        xlin1(168) = -5.436396122931865e+01_wp
        xlin1(169) = -1.702038609057893e+01_wp
        xlin1(170) = 4.969068890501028e-03_wp
        xlin1(171) = -3.519006024232350e+00_wp
        xlin1(172) = 3.284752951627190e+01_wp
        xlin1(173) = 2.783127622573391e+00_wp
        xlin1(174) = -5.624072554917517e-01_wp
        xlin1(175) = 5.370237229208703e+00_wp
        xlin1(176) = 9.768755790120242e-01_wp
        xlin1(177) = -3.206502231938864e+00_wp
        xlin1(178) = -1.267702449737564e-01_wp
        xlin1(179) = -1.430186406716300e-01_wp
        xlin1(180) = 1.296229968602821e-01_wp
        xlin2(1)   = 8.145631479853698e-01_wp
        xlin2(2)   = 1.770730113974052e+01_wp
        xlin2(3)   = -2.173959339189124e+01_wp
        xlin2(4)   = 1.961131724453491e+01_wp
        xlin2(5)   = -4.104999264119586e+01_wp
        xlin2(6)   = 2.337793724047610e+01_wp
        xlin2(7)   = -5.449082199443552e+00_wp
        xlin2(8)   = 6.201541014865459e-01_wp
        xlin2(9)   = -4.728365816655994e-02_wp
        xlin2(10)  = 3.321586948535380e-03_wp
        xlin2(11)  = -6.652255446067896e+01_wp
        xlin2(12)  = 5.962375314289732e+01_wp
        xlin2(13)  = 1.469998744381738e+01_wp
        xlin2(14)  = 1.153758314494926e+00_wp
        xlin2(15)  = -7.591897754372836e+00_wp
        xlin2(16)  = 2.103874186533575e+00_wp
        xlin2(17)  = -1.702965824754904e-01_wp
        xlin2(18)  = 1.346878024842888e-04_wp
        xlin2(19)  = -1.178260275080853e+02_wp
        xlin2(20)  = 3.489556535478315e+01_wp
        xlin2(21)  = -1.251243111695901e+00_wp
        xlin2(22)  = 2.471519774632120e-01_wp
        xlin2(23)  = -2.474793708434775e-01_wp
        xlin2(24)  = 2.687242988983517e-02_wp
        xlin2(25)  = -2.026760191030008e+01_wp
        xlin2(26)  = 2.756556273693077e+00_wp
        xlin2(27)  = 8.357078625422904e-02_wp
        xlin2(28)  = -1.640242746249622e-02_wp
        xlin2(29)  = -6.876043147213573e-01_wp
        xlin2(30)  = 2.609337260717558e-02_wp
        xlin2(31)  = 3.002603577505805e+00_wp
        xlin2(32)  = 1.498828705957069e+01_wp
        xlin2(33)  = -2.086306249380632e+01_wp
        xlin2(34)  = 2.081499551542078e+01_wp
        xlin2(35)  = -3.171758822478321e+01_wp
        xlin2(36)  = 1.723622638287100e+01_wp
        xlin2(37)  = -4.028169514216751e+00_wp
        xlin2(38)  = 4.671596528203490e-01_wp
        xlin2(39)  = -3.670331530469922e-02_wp
        xlin2(40)  = 2.629265006741015e-03_wp
        xlin2(41)  = -8.927371067901906e+01_wp
        xlin2(42)  = 8.867171901574056e+01_wp
        xlin2(43)  = -1.133027069526800e+02_wp
        xlin2(44)  = 1.602734361579233e+02_wp
        xlin2(45)  = -1.023650687832928e+02_wp
        xlin2(46)  = 2.471829619377287e+01_wp
        xlin2(47)  = -2.300433949989794e+00_wp
        xlin2(48)  = 2.592872235150836e-01_wp
        xlin2(49)  = -4.379218716825992e-02_wp
        xlin2(50)  = -1.231186961261672e+01_wp
        xlin2(51)  = 1.420549968400411e+01_wp
        xlin2(52)  = 1.232687494983029e+00_wp
        xlin2(53)  = 1.841045801282512e+01_wp
        xlin2(54)  = -6.141054432632874e+00_wp
        xlin2(55)  = 2.034847145705890e-01_wp
        xlin2(56)  = 8.860867509074010e-02_wp
        xlin2(57)  = -9.644503856306309e+01_wp
        xlin2(58)  = 2.392130829529128e+01_wp
        xlin2(59)  = -3.146354866365471e+00_wp
        xlin2(60)  = 9.562642394369082e-01_wp
        xlin2(61)  = -2.269962411748502e-01_wp
        xlin2(62)  = -7.627358272973419e+00_wp
        xlin2(63)  = 4.570619067451160e-01_wp
        xlin2(64)  = 2.237747868045461e-01_wp
        xlin2(65)  = -3.169396242810660e-01_wp
        xlin2(66)  = 1.231481888044332e+01_wp
        xlin2(67)  = -5.861438736123795e+01_wp
        xlin2(68)  = 5.294643266448595e+01_wp
        xlin2(69)  = -4.088837194494829e+00_wp
        xlin2(70)  = 5.495540512301918e+00_wp
        xlin2(71)  = -6.300138454926564e+00_wp
        xlin2(72)  = 1.765098522094756e+00_wp
        xlin2(73)  = -1.780074643720451e-01_wp
        xlin2(74)  = 4.114464551103923e-03_wp
        xlin2(75)  = 8.917318947290573e+01_wp
        xlin2(76)  = -1.770180683611587e+01_wp
        xlin2(77)  = -5.602678714811219e+00_wp
        xlin2(78)  = 2.973755372611144e+01_wp
        xlin2(79)  = 5.605823953149696e+00_wp
        xlin2(80)  = -4.949227356210042e+00_wp
        xlin2(81)  = 2.420977258547833e-01_wp
        xlin2(82)  = 8.167453315236323e-02_wp
        xlin2(83)  = -1.135834657778073e+02_wp
        xlin2(84)  = 8.290397544044251e+01_wp
        xlin2(85)  = -7.296380654367582e+01_wp
        xlin2(86)  = 1.347009356803447e+01_wp
        xlin2(87)  = -2.500870639262219e-01_wp
        xlin2(88)  = -1.272556171177672e-01_wp
        xlin2(89)  = 5.669201615035288e+00_wp
        xlin2(90)  = 4.449130216981585e+00_wp
        xlin2(91)  = -1.272014612117903e+00_wp
        xlin2(92)  = 2.139168528456700e-01_wp
        xlin2(93)  = -4.777213854750650e-02_wp
        xlin2(94)  = -1.275863250945433e-02_wp
        xlin2(95)  = -1.642934979019225e+01_wp
        xlin2(96)  = 5.594030291218376e+01_wp
        xlin2(97)  = -8.305010271115899e+01_wp
        xlin2(98)  = 3.140219381847346e+01_wp
        xlin2(99)  = -2.534163029719606e+00_wp
        xlin2(100) = 1.142479699516450e-01_wp
        xlin2(101) = -1.263578178350560e-01_wp
        xlin2(102) = 1.629095589290311e-02_wp
        xlin2(103) = -1.136985820244078e+02_wp
        xlin2(104) = 3.262096657063820e+01_wp
        xlin2(105) = -1.087121813266270e+02_wp
        xlin2(106) = 2.673092202326776e+01_wp
        xlin2(107) = -5.444160084123269e-01_wp
        xlin2(108) = 6.138628773543448e-01_wp
        xlin2(109) = -2.227712501792724e-01_wp
        xlin2(110) = 5.512545979043862e+01_wp
        xlin2(111) = 2.559974311865491e+01_wp
        xlin2(112) = 5.216798906963263e-01_wp
        xlin2(113) = -1.496009262806457e+00_wp
        xlin2(114) = 2.638941938405173e-01_wp
        xlin2(115) = -9.723650088603947e+00_wp
        xlin2(116) = 6.868905539663074e-01_wp
        xlin2(117) = -4.177757286042882e-01_wp
        xlin2(118) = 2.781285313472535e-01_wp
        xlin2(119) = 1.563554076363040e+01_wp
        xlin2(120) = -6.542137033344432e+00_wp
        xlin2(121) = 2.193791357437970e+01_wp
        xlin2(122) = -1.576158545476312e+01_wp
        xlin2(123) = 2.587945022189958e+00_wp
        xlin2(124) = -3.893301430803128e-02_wp
        xlin2(125) = -5.190707866429820e-03_wp
        xlin2(126) = 1.367205185988990e+02_wp
        xlin2(127) = 1.188103897312893e+01_wp
        xlin2(128) = 1.606164459906229e+01_wp
        xlin2(129) = -9.262867579726493e+00_wp
        xlin2(130) = 6.271894032241021e-01_wp
        xlin2(131) = 2.285179285451767e-01_wp
        xlin2(132) = -7.103397905922931e+01_wp
        xlin2(133) = 3.121061890947882e+00_wp
        xlin2(134) = 6.354928967335751e-01_wp
        xlin2(135) = -4.163368478300122e-02_wp
        xlin2(136) = 1.020366559611126e+00_wp
        xlin2(137) = 1.837284088182978e-01_wp
        xlin2(138) = -2.832630010141147e+01_wp
        xlin2(139) = 1.275654589067907e+01_wp
        xlin2(140) = -4.603547842975937e-01_wp
        xlin2(141) = 1.995827205765568e+00_wp
        xlin2(142) = -5.232664167379651e-01_wp
        xlin2(143) = 1.710301952383997e-02_wp
        xlin2(144) = -8.590072407323460e+01_wp
        xlin2(145) = 1.674416687913584e+01_wp
        xlin2(146) = 1.953210020027019e-02_wp
        xlin2(147) = 2.520910415921704e-01_wp
        xlin2(148) = -3.423339521867145e-01_wp
        xlin2(149) = 1.169816725323726e+01_wp
        xlin2(150) = -1.761561723165250e+00_wp
        xlin2(151) = 3.013931319831423e-02_wp
        xlin2(152) = -4.165165498393316e-01_wp
        xlin2(153) = 1.386201564606759e+01_wp
        xlin2(154) = -9.091983474080589e+00_wp
        xlin2(155) = 2.768047770222391e-01_wp
        xlin2(156) = 5.763259032772141e-02_wp
        xlin2(157) = 3.029683243309628e-02_wp
        xlin2(158) = 1.740587394078216e+01_wp
        xlin2(159) = -6.199982916108434e+00_wp
        xlin2(160) = 9.489522270750911e-01_wp
        xlin2(161) = 2.535531600295191e-01_wp
        xlin2(162) = 1.698179626605738e-01_wp
        xlin2(163) = 2.305666871236848e-01_wp
        xlin2(164) = -2.509462487790768e+00_wp
        xlin2(165) = 2.075944787757823e+00_wp
        xlin2(166) = -2.275822797369400e-01_wp
        xlin2(167) = -2.065529414930606e-02_wp
        xlin2(168) = -1.074236086914455e+00_wp
        xlin2(169) = 6.803905517919883e-02_wp
        xlin2(170) = -2.615499139877834e-01_wp
        xlin2(171) = -1.335426493176628e-01_wp
        xlin2(172) = 1.840487213722926e-01_wp
        xlin2(173) = -1.569973762062688e-01_wp
        xlin2(174) = 3.050419120227986e-02_wp
        xlin2(175) = 2.653719569186909e-01_wp
        xlin2(176) = 1.094965343816732e-01_wp
        xlin2(177) = -2.187339773412837e-02_wp
        xlin2(178) = -1.873849846674096e-03_wp
        xlin2(179) = -5.603548924512692e-02_wp
        xlin2(180) = 3.558364313199398e-03_wp
        xlin3(1)   = -7.638557823809794e+01_wp
        xlin3(2)   = 1.703707427310439e+02_wp
        xlin3(3)   = -1.025601533591727e+02_wp
        xlin3(4)   = -5.296578652461962e+01_wp
        xlin3(5)   = 1.093595183012879e+02_wp
        xlin3(6)   = -6.588531176536077e+01_wp
        xlin3(7)   = 2.040793884905717e+01_wp
        xlin3(8)   = -3.479631418111812e+00_wp
        xlin3(9)   = 3.088169472711390e-01_wp
        xlin3(10)  = -1.111590540608995e-02_wp
        xlin3(11)  = -3.361962954982916e+02_wp
        xlin3(12)  = 2.446590342879648e+02_wp
        xlin3(13)  = -6.622902002288096e+01_wp
        xlin3(14)  = 1.191674304065106e+01_wp
        xlin3(15)  = -1.201463499030519e+01_wp
        xlin3(16)  = 6.122631636838099e+00_wp
        xlin3(17)  = -1.187879235377073e+00_wp
        xlin3(18)  = 7.926806077855969e-02_wp
        xlin3(19)  = -5.725326950390917e+01_wp
        xlin3(20)  = 2.868589802572544e+01_wp
        xlin3(21)  = -2.971320475494586e+01_wp
        xlin3(22)  = 1.008864240170261e+01_wp
        xlin3(23)  = -1.365325434437133e+00_wp
        xlin3(24)  = 6.554733610556563e-02_wp
        xlin3(25)  = -4.045186273038085e+01_wp
        xlin3(26)  = 1.874648094047193e+01_wp
        xlin3(27)  = -2.454684686634876e+00_wp
        xlin3(28)  = 6.765838416051127e-02_wp
        xlin3(29)  = -4.059257274704468e+00_wp
        xlin3(30)  = 2.165639317481252e-01_wp
        xlin3(31)  = -7.534329078823575e+01_wp
        xlin3(32)  = 1.796639204536162e+02_wp
        xlin3(33)  = -1.564733447866279e+02_wp
        xlin3(34)  = 3.775541937573816e+01_wp
        xlin3(35)  = 3.375720237731807e+01_wp
        xlin3(36)  = -3.005181622792514e+01_wp
        xlin3(37)  = 1.029486036217851e+01_wp
        xlin3(38)  = -1.791612509905809e+00_wp
        xlin3(39)  = 1.538953705090221e-01_wp
        xlin3(40)  = -5.019936328398713e-03_wp
        xlin3(41)  = -4.774634443417375e+01_wp
        xlin3(42)  = 3.863906310035196e+02_wp
        xlin3(43)  = 1.909693601494610e+00_wp
        xlin3(44)  = 7.896215873407198e+01_wp
        xlin3(45)  = -1.295233067421028e+02_wp
        xlin3(46)  = 4.285454822181649e+01_wp
        xlin3(47)  = -3.947563931965661e+00_wp
        xlin3(48)  = -2.453031697183805e-02_wp
        xlin3(49)  = -4.072260116080702e-03_wp
        xlin3(50)  = -1.512501065486055e+03_wp
        xlin3(51)  = 9.514665008404671e+02_wp
        xlin3(52)  = -6.253456734186879e+02_wp
        xlin3(53)  = 2.584487207114875e+02_wp
        xlin3(54)  = -5.147932463240105e+01_wp
        xlin3(55)  = 4.978444306783625e+00_wp
        xlin3(56)  = -2.242912909135229e-01_wp
        xlin3(57)  = -5.814167788685131e+02_wp
        xlin3(58)  = 2.286668120955779e+02_wp
        xlin3(59)  = -6.032832571577611e+01_wp
        xlin3(60)  = 8.065674990175529e+00_wp
        xlin3(61)  = -3.404150604733470e-01_wp
        xlin3(62)  = -3.061436485306763e+01_wp
        xlin3(63)  = 2.200637753633254e+00_wp
        xlin3(64)  = 1.057897733012469e-01_wp
        xlin3(65)  = 7.289136325949297e-03_wp
        xlin3(66)  = 1.896105176076998e+02_wp
        xlin3(67)  = -3.648470785265107e+02_wp
        xlin3(68)  = 2.595841400715867e+02_wp
        xlin3(69)  = -6.319152600268254e+01_wp
        xlin3(70)  = 8.494896397888191e+00_wp
        xlin3(71)  = -1.053255418482307e+01_wp
        xlin3(72)  = 5.467310061893341e+00_wp
        xlin3(73)  = -1.031401599859863e+00_wp
        xlin3(74)  = 6.598583555227435e-02_wp
        xlin3(75)  = 4.919698991437962e+02_wp
        xlin3(76)  = -1.523587257667702e+03_wp
        xlin3(77)  = 1.014167387686189e+03_wp
        xlin3(78)  = -6.508692392742124e+02_wp
        xlin3(79)  = 2.713378987371789e+02_wp
        xlin3(80)  = -5.735977155603892e+01_wp
        xlin3(81)  = 5.830582807234884e+00_wp
        xlin3(82)  = -2.440949592570337e-01_wp
        xlin3(83)  = -6.735007851702401e+02_wp
        xlin3(84)  = 7.901400404656945e+02_wp
        xlin3(85)  = -1.229867349843578e+02_wp
        xlin3(86)  = -3.009484917235691e+01_wp
        xlin3(87)  = 2.526743872042880e+00_wp
        xlin3(88)  = 5.970110295569673e-01_wp
        xlin3(89)  = -3.532469213851834e+02_wp
        xlin3(90)  = 4.005928931048262e+01_wp
        xlin3(91)  = 5.976702173098687e+00_wp
        xlin3(92)  = -1.614616921246194e+00_wp
        xlin3(93)  = -9.871235251782917e+00_wp
        xlin3(94)  = 3.290059681876524e-01_wp
        xlin3(95)  = -1.813749415305320e+02_wp
        xlin3(96)  = 3.109697012254384e+02_wp
        xlin3(97)  = -3.205683056506578e+01_wp
        xlin3(98)  = -3.657410087464648e+00_wp
        xlin3(99)  = -2.624482850677591e+01_wp
        xlin3(100) = 1.173325195979073e+01_wp
        xlin3(101) = -1.590994311096534e+00_wp
        xlin3(102) = 5.486195466674193e-02_wp
        xlin3(103) = -2.440497562681578e+02_wp
        xlin3(104) = 1.038109037697877e+03_wp
        xlin3(105) = -7.152807104622796e+02_wp
        xlin3(106) = 2.980486621854296e+02_wp
        xlin3(107) = -7.381873725420134e+01_wp
        xlin3(108) = 9.536084102100837e+00_wp
        xlin3(109) = -4.034487399636460e-01_wp
        xlin3(110) = 6.749359685299933e+02_wp
        xlin3(111) = -3.958160095292775e+02_wp
        xlin3(112) = 7.043339305226182e+01_wp
        xlin3(113) = 4.280031318910807e+00_wp
        xlin3(114) = -1.799431596685355e+00_wp
        xlin3(115) = 1.572684765203245e+02_wp
        xlin3(116) = -1.497852043495127e+01_wp
        xlin3(117) = 3.450483851677588e+00_wp
        xlin3(118) = -1.651601456234870e+00_wp
        xlin3(119) = 6.267522188135509e+01_wp
        xlin3(120) = -1.490721159969828e+02_wp
        xlin3(121) = 7.203117438959657e+01_wp
        xlin3(122) = -4.910294960766687e+01_wp
        xlin3(123) = 2.126436044051092e+01_wp
        xlin3(124) = -3.286263182633637e+00_wp
        xlin3(125) = 1.432920523463851e-01_wp
        xlin3(126) = 3.482161123451701e+02_wp
        xlin3(127) = -6.634624347379950e+02_wp
        xlin3(128) = 2.978996288910582e+02_wp
        xlin3(129) = -6.102181144139994e+01_wp
        xlin3(130) = 4.913896432859977e+00_wp
        xlin3(131) = 9.220877926952298e-02_wp
        xlin3(132) = -9.040330715822043e+01_wp
        xlin3(133) = 3.763878231734758e+01_wp
        xlin3(134) = -1.466110603121045e+01_wp
        xlin3(135) = 4.286621953679466e-01_wp
        xlin3(136) = -1.171184819157715e+01_wp
        xlin3(137) = -8.282936408674097e-01_wp
        xlin3(138) = 2.329172506772613e+01_wp
        xlin3(139) = 7.205432276781195e+01_wp
        xlin3(140) = -7.134506264990834e+01_wp
        xlin3(141) = 2.778750878256592e+01_wp
        xlin3(142) = -4.744233867644683e+00_wp
        xlin3(143) = 2.489849788505943e-01_wp
        xlin3(144) = -3.035952094511171e+02_wp
        xlin3(145) = 2.774672648114976e+02_wp
        xlin3(146) = -7.050154798680367e+01_wp
        xlin3(147) = 7.641437548457593e+00_wp
        xlin3(148) = -1.880492973015346e-01_wp
        xlin3(149) = -1.984086397967394e+01_wp
        xlin3(150) = 7.604978254357780e+00_wp
        xlin3(151) = 2.014923057356265e-01_wp
        xlin3(152) = 3.308108128880802e+00_wp
        xlin3(153) = -2.969993736353019e+01_wp
        xlin3(154) = -3.682715000434093e+01_wp
        xlin3(155) = 2.210331765988686e+01_wp
        xlin3(156) = -4.301393898599036e+00_wp
        xlin3(157) = 2.671370386619945e-01_wp
        xlin3(158) = 1.124784775141184e+02_wp
        xlin3(159) = -6.041199428997342e+01_wp
        xlin3(160) = 6.988933983482634e+00_wp
        xlin3(161) = -1.362655996129701e-01_wp
        xlin3(162) = -2.997626634973072e-02_wp
        xlin3(163) = -1.846855278627695e+00_wp
        xlin3(164) = 1.154319604248443e+01_wp
        xlin3(165) = 1.192127579658776e+01_wp
        xlin3(166) = -2.699670781911932e+00_wp
        xlin3(167) = 1.920744268281359e-01_wp
        xlin3(168) = -1.990490891239322e+01_wp
        xlin3(169) = 7.011127199328610e+00_wp
        xlin3(170) = -9.502550138734989e-02_wp
        xlin3(171) = 7.777433674389070e-01_wp
        xlin3(172) = -2.235484930636709e+00_wp
        xlin3(173) = -1.909433565086661e+00_wp
        xlin3(174) = 1.029497435932162e-01_wp
        xlin3(175) = 1.784299385800847e+00_wp
        xlin3(176) = -3.937636905107279e-01_wp
        xlin3(177) = 2.159590396303356e-01_wp
        xlin3(178) = 1.161199705458321e-01_wp
        xlin3(179) = -7.827610505544675e-02_wp
        xlin3(180) = -8.180052967228596e-03_wp
        !     distances are in Bohr, energies are in cm-1
        towave     = 219474.7e0_wp  ! hartree to cm-1 conversion
        xmh        = 1836.15264e0_wp   ! H nuclear mass (after substracting 1 for electron mass)
        xmd        = 3670.48293e0_wp   ! D nuclear mass (after substracting 1 for electron mass)
        hatom_bodc = 59.6818e0_wp  ! the BODC in cm-1 of H with an aug-cc-pVTZ basis set
        xma        = xmd
        xmb        = xmh
        xmc        = xmh
        round      = 1.e-2_wp
        ald        = 0.02e0_wp
        betad      = 0.72e0_wp
        betadsum   = 3.3e0_wp
        deldamp    = 0.25e0_wp
        rhodamp    = 0.01e0_wp
        beta       = 1.414254132358239e0_wp
        rnought    = 2.0e0_wp
        maxi2      = 0
        maxi       = 0
        call indexa2(11, 11, iind, jind, kind, ipar, maxi, maxi2)
        maxi2a = maxi2
        rsum   = r1x+r2x+r3x
        if((r1x.le.betad).or.(r2x.le.betad).or.(r3x.le.betad)&
            &    .or.(rsum.le.betadsum) )then
        damper     = 0.0e0_wp
        ddamper(:) = 0.0e0_wp
    else
        bigq  = r1x*r1x+r2x*r2x+r3x*r3x
        slit1 = 2.0e0_wp*r1x**2-r2x**2-r3x**2
        slit2 = r2x**2-r3x**2
        slit3 = slit1*slit1+3.0e0_wp*slit2*slit2
        slit  = sqrt(slit3)/bigq
        if( slit .le. deldamp )then
            damper     = 0.0e0_wp
            ddamper(:) = 0.0e0_wp
        else
            sdamp    = rhodamp/(deldamp-slit)
            rdampsum = ald/(betadsum-rsum)
            rdamp1   = ald/(betad-r1x)
            rdamp2   = ald/(betad-r2x)
            rdamp3   = ald/(betad-r3x)
            damper   = exp(rdampsum+rdamp1+rdamp2+rdamp3+sdamp)
            dsdamp   = sdamp/(deldamp-slit)
            ddamper(1)= damper*(dsdamp*slit*r1x* (4.0e0_wp*slit1/slit3 - 2.0e0_wp/bigq) +rdampsum/(betadsum-rsum)+ &
                & rdamp1/ (betad- r1x))
            ddamper(2)=damper*(-dsdamp*slit*r2x* (4.0e0_wp*(r1x*r1x-2.0e0_wp*r2x*r2x+r3x*r3x)/slit3+2.0e0_wp/bigq) + &
                & rdampsum/ (betadsum-rsum)+rdamp2/(betad-r2x))
            ddamper(3)=damper*(-dsdamp*slit*r3x* (4.0e0_wp*(r1x*r1x+r2x*r2x-2.0e0_wp*r3x*r3x)/slit3+2.0e0_wp/bigq) + &
                & rdampsum/ (betadsum-rsum)+rdamp3/(betad-r3x))
        endif
        endif
        rpass(1) = r1x
        rpass(3) = r2x
        rpass(2) = r3x
        call leps(rpass, v1, v2, xma, xmb, xmc, gtot, dgtot)
        sum1     = dgtot(2)
        dgtot(2) = dgtot(3)
        dgtot(3) = sum1
        bb2      = towave*gtot
        do i = 1, 3
            dbb2(i) = towave*dgtot(i)
        enddo
        bb3      = sqrt(round+(r1x-r3x)**2+(r1x-r2x)**2+(r2x-r3x)**2)
        dbb3(1)  = (2.0e0_wp*r1x-r2x-r3x)/bb3
        dbb3(2)  = (2.0e0_wp*r2x-r1x-r3x)/bb3
        dbb3(3)  = (2.0e0_wp*r3x-r1x-r2x)/bb3
        e1body   = (2.e0_wp + xmh/xmd) * hatom_bodc
        v        = e1body + ((xmd+xmh)/(2.e0_wp*xmd))*(bodc2body(r1x, dv1) +bodc2body(r2x, dv2)) + bodc2body(r3x, dv3)
        dv(1)    = ((xmd+xmh)/(2.e0_wp*xmd))*dv1
        dv(2)    = ((xmd+xmh)/(2.e0_wp*xmd))*dv2
        dv(3)    = dv3
        sum1     = 0.e0_wp
        sum2     = 0.e0_wp
        sum3     = 0.e0_wp
        dsum1(:) = 0.0e0_wp
        dsum2(:) = 0.0e0_wp
        dsum3(:) = 0.0e0_wp
        rho1     = exp(-beta*(r1x-rnought))
        drho1    = (1.0-beta*r1x)*rho1
        rho1     = r1x*rho1
        rho2     = exp(-beta*(r2x-rnought))
        drho2    = (1.0-beta*r2x)*rho2
        rho2     = r2x*rho2
        rho3     = exp(-beta*(r3x-rnought))
        drho3    = (1.0-beta*r3x)*rho3
        rho3     = r3x*rho3
        rho1t(0) = 1.e0_wp
        rho2t(0) = 1.e0_wp
        rho3t(0) = 1.e0_wp
        do it=1, 13
            rho1t(it)  = rho1t(it-1)*rho1
            rho2t(it)  = rho2t(it-1)*rho2
            rho3t(it)  = rho3t(it-1)*rho3
            drho1t(it) = it*rho1t(it-1)*drho1
            drho2t(it) = it*rho2t(it-1)*drho2
            drho3t(it) = it*rho3t(it-1)*drho3
        enddo
        do ii=1, maxi2a
            rprod    = rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
            sum1     = sum1+xlin1(ipar(ii))*rprod
            sum2     = sum2+xlin2(ipar(ii))*rprod
            sum3     = sum3+xlin3(ipar(ii))*rprod
            drprod1  = drho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
            drprod2  = rho1t(kind(ii))*drho2t(jind(ii))*rho3t(iind(ii))
            drprod3  = rho1t(kind(ii))*rho2t(jind(ii))*drho3t(iind(ii))
            dsum1(1) = dsum1(1)+xlin1(ipar(ii))*drprod1
            dsum1(2) = dsum1(2)+xlin1(ipar(ii))*drprod2
            dsum1(3) = dsum1(3)+xlin1(ipar(ii))*drprod3
            dsum2(1) = dsum2(1)+xlin2(ipar(ii))*drprod1
            dsum2(2) = dsum2(2)+xlin2(ipar(ii))*drprod2
            dsum2(3) = dsum2(3)+xlin2(ipar(ii))*drprod3
            dsum3(1) = dsum3(1)+xlin3(ipar(ii))*drprod1
            dsum3(2) = dsum3(2)+xlin3(ipar(ii))*drprod2
            dsum3(3) = dsum3(3)+xlin3(ipar(ii))*drprod3
        enddo
        v3c = sum1+sum2*bb2+sum3*bb3
        eee = v+damper*v3c
        do i = 1, 3
            dv(i) = dv(i) +ddamper(i)*v3c +damper*(sum2*dbb2(i)+sum3*dbb3(i) +dsum1(i)+dsum2(i)*bb2+dsum3(i)*bb3)
        enddo
    end subroutine dispatch
end module chempotpy_h3_gen_bh_2009_potbodcdh2_legacy_state_b8d1df42

!> @brief Compatibility interface for potbodcdh2.
subroutine potbodcdh2( r1x, r2x, r3x, eee, dv )
    use chempotpy_h3_gen_bh_2009_potbodcdh2_legacy_state_b8d1df42, only: dispatch
    use chempotpy_kinds,                                           only: wp
    implicit none
    real(wp) :: r1x
    real(wp) :: r2x
    real(wp) :: r3x
    real(wp) :: eee
    real(wp) :: dv(1:3)
    call dispatch(r1x=r1x, r2x=r2x, r3x=r3x, eee=eee, dv=dv)
end subroutine potbodcdh2

!> @brief Evaluate bodc2body.
real(wp) function bodc2body( r1x, dv )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp) :: xlp(14), beta, re, v, r1x, dv, t
    integer  :: j, indx
    !     distances are in Bohr, energies are in cm-1
    !     The global minimum is at: R=2.3178250879779552  E= -17.5988842162493135 cm-1
    beta    = 3.4059169e0_wp
    xlp( 1) = -4.942622126e0_wp
    xlp( 2) = -4.742197464e1_wp
    xlp( 3) = -1.115075510e2_wp
    xlp( 4) = -1.419033104e2_wp
    xlp( 5) = -1.165969889e2_wp
    xlp( 6) = -6.279822284e1_wp
    xlp( 7) = -2.626186059e1_wp
    xlp( 8) = -1.795750569e1_wp
    xlp( 9) = -6.826792412e0_wp
    xlp(10) = 7.563692060e0_wp
    xlp(11) = 2.630227612e0_wp
    xlp(12) = -2.595885106e0_wp
    xlp(13) = 5.647782847e-1_wp
    xlp(14) = -4.042171086e-2_wp
    re      = 1.401e0_wp
    v       = xlp(14)
    dv      = -beta*xlp(14)
    indx    = 14
    do j=2, 14
        indx = indx-1
        v    = v*(r1x-re)+xlp(indx)
        dv   = dv*(r1x-re)+indx*xlp(indx+1)-beta*xlp(indx)
    enddo
    t         = exp(-beta*(r1x-re))
    v         = v*t
    dv        = dv*t
    bodc2body = v
end

!> @brief Evaluate indexa2.
subroutine indexa2( mbig, mtop, iind, jind, kind, ipar, maxi, maxi2 )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)    :: mbig, mtop
    integer, intent(inout) :: iind(5000), jind(5000), kind(5000), ipar(5000), maxi, maxi2
    integer                :: l, l2, i, j, k, isum
    l  = 0
    l2 = 0
    do i=0, min(mtop, mbig)
        do j=0, min(mtop, mbig)
            do k=j, min(mtop, mbig)
                isum = i+j+k
                if( isum .le. mbig )then
                    if(((isum-i).gt.0).and.((isum-j).gt.0).and.&
                        &        ((isum-k).gt.0))then
                    l2      = l2+1
                    l       = l+1
                    ipar(l) = l2
                    iind(l) = i
                    jind(l) = j
                    kind(l) = k
                    if( j .ne. k )then
                        l       = l+1
                        ipar(l) = l2
                        iind(l) = i
                        jind(l) = k
                        kind(l) = j
                    endif
                endif
                endif
            enddo
        enddo
    enddo
    maxi  = l2
    maxi2 = l
end

!> @brief Evaluate leps.
subroutine leps ( r, v1, v2, xma, xmb, xmc, gtot, dgtot )
    !   r1 is the distance between A and B
    !   r2 is the distance between B and C
    !   r3 is the distance between A and C
    !   v1 is the lower energy root
    !   v2 is the upper energy root
    !   xma, xmb, xmc are the atomic masses
    !   ga, gb, gc are the contributions to the nonadiabatic corrections
    !      from each atom
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: i, j
    real(wp) :: r(3), v1, v2, xma, xmb, xmc, gtot, dgtot(3)
    real(wp) :: ex, vs(3), vt(3), dvs(3), dvt(3)
    real(wp) :: h11, h12, h22, hsum, hdif, dhdif(3), rad, delh
    real(wp) :: f12(3), cos1, cos2, cos3, ga, gb, gc
    real(wp) :: d2vs(3), d2vt(3), dh11(3), dh12(3), dh22(3)
    real(wp) :: d2h12(3), d2hdif(3), drad(3), df12(3,3)
    real(wp) :: dcos1(3), dcos2(3), dcos3(3)
    real(wp) :: d(3), re(3), beta(3), sato(3), rt3
    d    = 109.4583e0_wp/627.509552e0_wp
    re   = 0.7412871e0_wp/0.52917725e0_wp
    beta = 1.973536e0_wp*0.52917725e0_wp
    sato = 0.132e0_wp
    rt3  = sqrt(3.e0_wp)
    ! Compute single and triplet diatomic potentials
    do i = 1, 3
        ex      = exp(-beta(i)*(r(i)-re(i)))
        vs(i)   = d(i)*(ex - 2.0e0_wp)*ex
        dvs(i)  = -2.0e0_wp*d(i)*beta(i)*(ex-1.0e0_wp)*ex
        d2vs(i) = 2.0e0_wp*d(i)*beta(i)*beta(i)*(2.e0_wp*ex-1.0e0_wp)*ex
        vt(i)   = 0.5e0_wp*d(i)*(ex + 2.0)*ex*(1.0-sato(i))/(1.0+sato(i))
        dvt(i)  = -d(i)*beta(i)*(ex+1.0)*ex*(1.0-sato(i))/(1.0+sato(i))
        d2vt(i) = d(i)*beta(i)*beta(i)*(2.0e0_wp*ex+1.0)*ex *(1.0-sato(i))/(1.0+sato(i))
    enddo
    ! Construct 2x2 DIM Hamiltonian
    h11 = d(2) + vs(2) + 0.25e0_wp*(vs(1)+vs(3)) + 0.75e0_wp*(vt(1)+vt(3))
    h12 = 0.25e0_wp*rt3*(vs(1) - vs(3) - vt(1) + vt(3))
    h22 = d(2) + vt(2) + 0.75*(vs(1)+vs(3)) + 0.25e0_wp*(vt(1)+vt(3))
    ! Calculate eigenvalues (V1 and V2) and quantities needed for matrix elements
    !   of d/dRi
    hsum     = h11 + h22
    hdif     = h22 - h11
    dhdif(1) = 0.5e0_wp*(dvs(1)-dvt(1))
    dhdif(2) = -(dvs(2) - dvt(2))
    dhdif(3) = 0.5e0_wp*(dvs(3)-dvt(3))
    rad      = 0.25e0_wp*hdif*hdif + h12*h12
    delh     = sqrt(rad)
    v1       = 0.5e0_wp*hsum - delh
    v2       = 0.5e0_wp*hsum + delh
    ! derivatives
    dh11(1)   = 0.25e0_wp*dvs(1)+0.75e0_wp*dvt(1)
    dh11(2)   = dvs(2)
    dh11(3)   = 0.25e0_wp*dvs(3)+0.75e0_wp*dvt(3)
    dh12(1)   = 0.25e0_wp*rt3*(dvs(1)-dvt(1))
    dh12(2)   = 0.0e0_wp
    dh12(3)   = -0.25e0_wp*rt3*(dvs(3)-dvt(3))
    dh22(1)   = 0.75e0_wp*dvs(1)+0.25e0_wp*dvt(1)
    dh22(2)   = dvt(2)
    dh22(3)   = 0.75e0_wp*dvs(3)+0.25e0_wp*dvt(3)
    d2h12(1)  = 0.25e0_wp*rt3*(d2vs(1)-d2vt(1))
    d2h12(2)  = 0.0e0_wp
    d2h12(3)  = -0.25e0_wp*rt3*(d2vs(3)-d2vt(3))
    d2hdif(1) = 0.5e0_wp*(d2vs(1) - d2vt(1))
    d2hdif(2) = -(d2vs(2) - d2vt(2))
    d2hdif(3) = 0.5e0_wp*(d2vs(3) - d2vt(3))
    drad(1)   = 0.5e0_wp*hdif*dhdif(1) +2.0e0_wp*h12*dh12(1)
    drad(2)   = 0.5e0_wp*hdif*dhdif(2) +2.0e0_wp*h12*dh12(2)
    drad(3)   = 0.5e0_wp*hdif*dhdif(3) +2.0e0_wp*h12*dh12(3)
    ! Calculate matrix elements of d/dRi; note that at conical intersection f12
    !   diverges, but is set to zero in the code.
    do i = 1, 3
        f12(i) = 0.25e0_wp*(h12*dhdif(i) - hdif*dh12(i))
        if( rad .ne. 0.0e0_wp )then
            f12(i) = f12(i)/rad
            do j = 1, 3
                df12(i,j) = (0.25e0_wp*(dh12(j)*dhdif(i)-dhdif(j)*dh12(i)) -f12(i)*drad(j))/rad
            enddo
            df12(i,i) = df12(i,i) + 0.25e0_wp*(h12*d2hdif(i)-hdif*d2h12(i))/rad
        endif
    enddo
    ! Calculate diagonal matrix element of kinetic energy operator
    cos1     = (r(2)*r(2) + r(3)*r(3) - r(1)*r(1))/(2.0*r(2)*r(3))
    cos2     = (r(1)*r(1) + r(3)*r(3) - r(2)*r(2))/(2.0*r(1)*r(3))
    cos3     = (r(1)*r(1) + r(2)*r(2) - r(3)*r(3))/(2.0*r(1)*r(2))
    ga       = (f12(1)*(f12(1) + 2.0*cos2*f12(3)) + f12(3)*f12(3))/(2.e0_wp*xma)
    gb       = (f12(1)*(f12(1) + 2.0*cos3*f12(2)) + f12(2)*f12(2))/(2.e0_wp*xmb)
    gc       = (f12(2)*(f12(2) + 2.0*cos1*f12(3)) + f12(3)*f12(3))/(2.e0_wp*xmc)
    gtot     = ga+gb+gc
    dcos1(1) = -r(1)/(r(2)*r(3))
    dcos1(2) = 1.0e0_wp/r(3) - cos1/r(2)
    dcos1(3) = 1.0e0_wp/r(2) - cos1/r(3)
    dcos2(1) = 1.0e0_wp/r(3) - cos2/r(1)
    dcos2(2) = -r(2)/(r(1)*r(3))
    dcos2(3) = 1.0e0_wp/r(1) - cos2/r(3)
    dcos3(1) = 1.0e0_wp/r(2) - cos3/r(1)
    dcos3(2) = 1.0e0_wp/r(1) - cos3/r(2)
    dcos3(3) = -r(3)/(r(1)*r(2))
    do i = 1, 3
        ! contribution from ga
        dgtot(i) = (f12(1)*(df12(1,i) + dcos2(i)*f12(3) + cos2*df12(3,i)) + cos2*df12(1,i)*f12(3) + f12(3)*df12(3, &
            & i))/ xma
        ! contribution from ga
        dgtot(i) = dgtot(i) + (f12(1)*(df12(1,i) + dcos3(i)*f12(2) + cos3*df12(2,i)) + cos3*df12(1,i)*f12(2) + f12(2)* &
            & df12(2,i))/xmb
        ! contribution from ga
        dgtot(i) = dgtot(i) + (f12(2)*(df12(2,i) + dcos1(i)*f12(3) + cos1*df12(3,i)) + cos1*df12(2,i)*f12(3) + f12(3)* &
            & df12(3,i))/xmc
    enddo
end

module chempotpy_h3_gen_bh_2009_potcbs_legacy_state_b8d1df42
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r1x, r2x, r3x, eee, dv )
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter  :: n2      = 5000
        real(wp), parameter :: eshift  = 0.1744759989649372e0_wp
        real(wp), parameter :: beta    = 0.9246250191375499e0_wp
        real(wp), parameter :: rnought = 2.e0_wp
        optional :: r1x, r2x, r3x, eee, dv
        integer  :: iind(n2), jind(n2), kind(n2), ipar(n2), maxi, maxi2
        integer  :: maxi2a, maxi2b, ii, i, it
        real(wp) :: epslon, epscem, ald, betad
        real(wp) :: r1x, r2x, r3x, eee, dv(:)
        real(wp) :: etrip1, etrip2, etrip3, detrip1, detrip2, detrip3
        real(wp) :: es1, es2, es3, des1, des2, des3
        real(wp) :: q1, q2, q3, xj1, xj2, xj3, xj, vlondon
        real(wp) :: dq1, dq2, dq3, dxj1, dxj2, dxj3, dvlon(3)
        real(wp) :: damper, ddamp(3), rdamp1, rdamp2, rdamp3
        real(wp) :: bb, dbb(3), warp, dwarp(3), rho1, rho2, rho3
        real(wp) :: rho1t(0:12), rho2t(0:12), rho3t(0:12)
        real(wp) :: sum1, sum2, sum3, sum4
        real(wp) :: dsum1(3), dsum2(3), dsum3(3), dsum4(3)
        real(wp) :: rprod, t1, v3c
        real(wp) :: xlin1(89), xlin2(89), xlin3(71), xlin4(71)
        ald       = 0.e0_wp
        bb        = 0.e0_wp
        betad     = 0.e0_wp
        damper    = 0.e0_wp
        dbb       = 0.e0_wp
        ddamp     = 0.e0_wp
        des1      = 0.e0_wp
        des2      = 0.e0_wp
        des3      = 0.e0_wp
        detrip1   = 0.e0_wp
        detrip2   = 0.e0_wp
        detrip3   = 0.e0_wp
        dq1       = 0.e0_wp
        dq2       = 0.e0_wp
        dq3       = 0.e0_wp
        dsum1     = 0.e0_wp
        dsum2     = 0.e0_wp
        dsum3     = 0.e0_wp
        dsum4     = 0.e0_wp
        dvlon     = 0.e0_wp
        dwarp     = 0.e0_wp
        dxj1      = 0.e0_wp
        dxj2      = 0.e0_wp
        dxj3      = 0.e0_wp
        epscem    = 0.e0_wp
        epslon    = 0.e0_wp
        es1       = 0.e0_wp
        es2       = 0.e0_wp
        es3       = 0.e0_wp
        etrip1    = 0.e0_wp
        etrip2    = 0.e0_wp
        etrip3    = 0.e0_wp
        q1        = 0.e0_wp
        q2        = 0.e0_wp
        q3        = 0.e0_wp
        rdamp1    = 0.e0_wp
        rdamp2    = 0.e0_wp
        rdamp3    = 0.e0_wp
        rho1      = 0.e0_wp
        rho1t     = 0.e0_wp
        rho2      = 0.e0_wp
        rho2t     = 0.e0_wp
        rho3      = 0.e0_wp
        rho3t     = 0.e0_wp
        rprod     = 0.e0_wp
        sum1      = 0.e0_wp
        sum2      = 0.e0_wp
        sum3      = 0.e0_wp
        sum4      = 0.e0_wp
        t1        = 0.e0_wp
        v3c       = 0.e0_wp
        vlondon   = 0.e0_wp
        warp      = 0.e0_wp
        xj        = 0.e0_wp
        xj1       = 0.e0_wp
        xj2       = 0.e0_wp
        xj3       = 0.e0_wp
        xlin1     = 0.e0_wp
        xlin2     = 0.e0_wp
        xlin3     = 0.e0_wp
        xlin4     = 0.e0_wp
        xlin1(1)  = -0.1682626669809615e0_wp
        xlin1(2)  = 1.697952225005269e0_wp
        xlin1(3)  = 3.959367441888326e0_wp
        xlin1(4)  = -7.159854944676235e0_wp
        xlin1(5)  = -2.462488887191634e0_wp
        xlin1(6)  = 7.998434264390882e0_wp
        xlin1(7)  = -3.665712084901418e0_wp
        xlin1(8)  = -0.2155463656483861e0_wp
        xlin1(9)  = 0.5620418518154215e0_wp
        xlin1(10) = -0.1369812453279634e0_wp
        xlin1(11) = 1.0063614579308723e-2_wp
        xlin1(12) = 41.07358514471853e0_wp
        xlin1(13) = 67.02217912830267e0_wp
        xlin1(14) = -100.0014758720550e0_wp
        xlin1(15) = 2.481875721304641e0_wp
        xlin1(16) = 32.46244740624601e0_wp
        xlin1(17) = -15.32460432587821e0_wp
        xlin1(18) = 3.814631583515083e0_wp
        xlin1(19) = -0.3155985331636701e0_wp
        xlin1(20) = -4.2198516063816239e-3_wp
        xlin1(21) = -157.7203098415152e0_wp
        xlin1(22) = 17.11141015275416e0_wp
        xlin1(23) = -61.48638359701442e0_wp
        xlin1(24) = 4.244841908917134e0_wp
        xlin1(25) = 5.300393723900312e0_wp
        xlin1(26) = -0.6170558358973661e0_wp
        xlin1(27) = -6.2143674437616274e-3_wp
        xlin1(28) = 5.785766830955630e0_wp
        xlin1(29) = -0.1133784523542150e0_wp
        xlin1(30) = -0.2539145873400511e0_wp
        xlin1(31) = -4.3406559735048131e-2_wp
        xlin1(32) = 2.4548851428702549e-2_wp
        xlin1(33) = -1.092841685240530e0_wp
        xlin1(34) = 0.8131178021738225e0_wp
        xlin1(35) = -7.7933784638615652e-2_wp
        xlin1(36) = -0.1484455564689179e0_wp
        xlin1(37) = -7.351078313940695e0_wp
        xlin1(38) = 128.7480567271992e0_wp
        xlin1(39) = 199.6058760392980e0_wp
        xlin1(40) = -50.18444648711906e0_wp
        xlin1(41) = -178.2184990028168e0_wp
        xlin1(42) = 129.6832433574917e0_wp
        xlin1(43) = -37.57036725437851e0_wp
        xlin1(44) = 5.747838042162251e0_wp
        xlin1(45) = 0.1434298979682676e0_wp
        xlin1(46) = -8.3110343245055174e-2_wp
        xlin1(47) = 428.7635440340811e0_wp
        xlin1(48) = 495.9564113906064e0_wp
        xlin1(49) = -1742.123837168904e0_wp
        xlin1(50) = -259.0564024322437e0_wp
        xlin1(51) = 250.2352551214438e0_wp
        xlin1(52) = -44.15232947384926e0_wp
        xlin1(53) = 0.9963955068997190e0_wp
        xlin1(54) = 0.3136096547683901e0_wp
        xlin1(55) = -1032.921631265462e0_wp
        xlin1(56) = 199.9708129167335e0_wp
        xlin1(57) = -135.4683186684037e0_wp
        xlin1(58) = -14.50534672808297e0_wp
        xlin1(59) = 5.650612567917722e0_wp
        xlin1(60) = -0.1473535040913008e0_wp
        xlin1(61) = 1.351710421967908e0_wp
        xlin1(62) = 29.60972395839821e0_wp
        xlin1(63) = -7.224575522867537e0_wp
        xlin1(64) = 0.3582134126835411e0_wp
        xlin1(65) = 0.2589923362572861e0_wp
        xlin1(66) = 0.1819402887127863e0_wp
        xlin1(67) = 1312.633253770123e0_wp
        xlin1(68) = 5189.281073844522e0_wp
        xlin1(69) = -2727.020597875299e0_wp
        xlin1(70) = -122.2193440638810e0_wp
        xlin1(71) = 152.4859537357477e0_wp
        xlin1(72) = -5.562127541501408e0_wp
        xlin1(73) = -1.549418235204051e0_wp
        xlin1(74) = 2647.710764337011e0_wp
        xlin1(75) = 502.1991025406292e0_wp
        xlin1(76) = -152.3069488960497e0_wp
        xlin1(77) = -5.740041181528469e0_wp
        xlin1(78) = 1.362209926084833e0_wp
        xlin1(79) = 17.96504193139269e0_wp
        xlin1(80) = 10.13930464010770e0_wp
        xlin1(81) = -2.287463377365149e0_wp
        xlin1(82) = 1.932492687974671e0_wp
        xlin1(83) = -3337.181754475553e0_wp
        xlin1(84) = 182.9671187483376e0_wp
        xlin1(85) = 20.29606065582956e0_wp
        xlin1(86) = 1.192500349837268e0_wp
        xlin1(87) = -22.69079850990068e0_wp
        xlin1(88) = -0.5967747962499963e0_wp
        xlin1(89) = 1.091111370927579e0_wp
        xlin2(1)  = 3.2289532152874986e-2_wp
        xlin2(2)  = -8.6686506736620789e-2_wp
        xlin2(3)  = -0.3266024449171150e0_wp
        xlin2(4)  = -0.7972798579742405e0_wp
        xlin2(5)  = 2.812854025206053e0_wp
        xlin2(6)  = -2.374871328104482e0_wp
        xlin2(7)  = 0.6089164699935086e0_wp
        xlin2(8)  = 0.1635658322754857e0_wp
        xlin2(9)  = -0.1263080841790984e0_wp
        xlin2(10) = 2.5221027611065153e-2_wp
        xlin2(11) = -1.6677760633386995e-3_wp
        xlin2(12) = -0.3372069966080046e0_wp
        xlin2(13) = -8.377537389356892e0_wp
        xlin2(14) = 3.694668699246826e0_wp
        xlin2(15) = 11.36868357221171e0_wp
        xlin2(16) = -11.88889428297186e0_wp
        xlin2(17) = 4.270226804928180e0_wp
        xlin2(18) = -0.4526863369836274e0_wp
        xlin2(19) = -5.2518975083660345e-2_wp
        xlin2(20) = 9.0223742463199624e-3_wp
        xlin2(21) = 38.03127084604298e0_wp
        xlin2(22) = -12.82065963469141e0_wp
        xlin2(23) = -5.663705994784566e0_wp
        xlin2(24) = 7.300080908068341e0_wp
        xlin2(25) = -0.5157232521089076e0_wp
        xlin2(26) = -0.2200954472805637e0_wp
        xlin2(27) = 2.0686337667851093e-2_wp
        xlin2(28) = 14.19614073502447e0_wp
        xlin2(29) = 1.726255538125508e0_wp
        xlin2(30) = -2.526391344985463e0_wp
        xlin2(31) = 6.0598667264186397e-2_wp
        xlin2(32) = 3.8560006141770319e-2_wp
        xlin2(33) = 2.306725001955670e0_wp
        xlin2(34) = -0.3197660708712612e0_wp
        xlin2(35) = -1.9982819587147187e-2_wp
        xlin2(36) = 0.1036787525350017e0_wp
        xlin2(37) = -2.030231243513758e0_wp
        xlin2(38) = 25.16474735933574e0_wp
        xlin2(39) = -14.13898387128746e0_wp
        xlin2(40) = -62.01219930765056e0_wp
        xlin2(41) = 75.56471352787402e0_wp
        xlin2(42) = -30.25154027581496e0_wp
        xlin2(43) = 4.257663341071418e0_wp
        xlin2(44) = 0.6142213624561507e0_wp
        xlin2(45) = -0.2458345005137628e0_wp
        xlin2(46) = 1.9135101060142316e-2_wp
        xlin2(47) = 168.2482649287792e0_wp
        xlin2(48) = -270.3587962365308e0_wp
        xlin2(49) = 72.46467150857366e0_wp
        xlin2(50) = 32.12278113625301e0_wp
        xlin2(51) = 11.89435258275814e0_wp
        xlin2(52) = -5.521268211450904e0_wp
        xlin2(53) = 0.8934653854937415e0_wp
        xlin2(54) = -7.6805320500971247e-2_wp
        xlin2(55) = -422.2672942154591e0_wp
        xlin2(56) = -35.92489169679082e0_wp
        xlin2(57) = 14.74665733902546e0_wp
        xlin2(58) = 6.020706922914545e0_wp
        xlin2(59) = 2.4538599994799737e-3_wp
        xlin2(60) = -0.1370624616096560e0_wp
        xlin2(61) = 44.47654636843973e0_wp
        xlin2(62) = -11.18166310447657e0_wp
        xlin2(63) = 3.6110526441458540e-2_wp
        xlin2(64) = 0.1539535582849851e0_wp
        xlin2(65) = 2.065981296572923e0_wp
        xlin2(66) = -0.1624439939014451e0_wp
        xlin2(67) = -658.3723053145102e0_wp
        xlin2(68) = -883.4947620258348e0_wp
        xlin2(69) = -263.8479489759550e0_wp
        xlin2(70) = 44.68359598895833e0_wp
        xlin2(71) = 29.88970883430300e0_wp
        xlin2(72) = -5.663152586850906e0_wp
        xlin2(73) = 0.2413124485342995e0_wp
        xlin2(74) = 656.9813540494517e0_wp
        xlin2(75) = -91.56710742448996e0_wp
        xlin2(76) = -30.38175110815412e0_wp
        xlin2(77) = 2.066013906655882e0_wp
        xlin2(78) = 0.5314632792073760e0_wp
        xlin2(79) = 54.99422880821747e0_wp
        xlin2(80) = -0.7887249891070589e0_wp
        xlin2(81) = -0.9164068632403487e0_wp
        xlin2(82) = 0.6576346543435602e0_wp
        xlin2(83) = 292.2665706493461e0_wp
        xlin2(84) = -32.04152511166628e0_wp
        xlin2(85) = -2.104943871500433e0_wp
        xlin2(86) = -0.7251129002870176e0_wp
        xlin2(87) = 4.709122797258667e0_wp
        xlin2(88) = 1.119971732384046e0_wp
        xlin2(89) = -4.066519840824071e0_wp
        xlin3(1)  = 6.3998636332917941e-2_wp
        xlin3(2)  = -0.5826083150212633e0_wp
        xlin3(3)  = -3.294968443823759e0_wp
        xlin3(4)  = 5.751817019955082e0_wp
        xlin3(5)  = -1.486176533347650e0_wp
        xlin3(6)  = -0.6845411439568367e0_wp
        xlin3(7)  = -0.9083627566591157e0_wp
        xlin3(8)  = 1.305678536333854e0_wp
        xlin3(9)  = -0.4821780067143284e0_wp
        xlin3(10) = 5.5182580699196751e-2_wp
        xlin3(11) = -22.57213192201081e0_wp
        xlin3(12) = -61.37436407473290e0_wp
        xlin3(13) = 52.64511100643983e0_wp
        xlin3(14) = 29.99521825023724e0_wp
        xlin3(15) = -38.31341936937510e0_wp
        xlin3(16) = 14.86681566870566e0_wp
        xlin3(17) = -3.460352146284503e0_wp
        xlin3(18) = 0.1343449644795501e0_wp
        xlin3(19) = 53.21459427503272e0_wp
        xlin3(20) = 64.60699757442232e0_wp
        xlin3(21) = 30.89466488573955e0_wp
        xlin3(22) = 8.412732049211240e0_wp
        xlin3(23) = -4.274560634850636e0_wp
        xlin3(24) = -0.1788031871862021e0_wp
        xlin3(25) = -46.47154313325355e0_wp
        xlin3(26) = 24.88995914944731e0_wp
        xlin3(27) = -2.757433396373881e0_wp
        xlin3(28) = -0.3635530438220282e0_wp
        xlin3(29) = 1.049428472748707e0_wp
        xlin3(30) = -0.9909460260393127e0_wp
        xlin3(31) = 7.623381850055868e0_wp
        xlin3(32) = -75.76387696233419e0_wp
        xlin3(33) = -162.8879929453359e0_wp
        xlin3(34) = -29.98317728892977e0_wp
        xlin3(35) = 168.7972614497797e0_wp
        xlin3(36) = -102.0864921135087e0_wp
        xlin3(37) = 31.38363768116054e0_wp
        xlin3(38) = -6.084136544334489e0_wp
        xlin3(39) = -2.9400115678577254e-2_wp
        xlin3(40) = -417.4808742201730e0_wp
        xlin3(41) = -601.3246381855070e0_wp
        xlin3(42) = 1111.017825017578e0_wp
        xlin3(43) = 749.0346863822200e0_wp
        xlin3(44) = -254.9798923160399e0_wp
        xlin3(45) = 37.02896288679148e0_wp
        xlin3(46) = -0.7010161055223685e0_wp
        xlin3(47) = 928.6100310886904e0_wp
        xlin3(48) = 513.3554862860957e0_wp
        xlin3(49) = 138.4248255502364e0_wp
        xlin3(50) = -0.9347938923678116e0_wp
        xlin3(51) = 3.867048618649862e0_wp
        xlin3(52) = -109.6794304656507e0_wp
        xlin3(53) = 10.30110710850976e0_wp
        xlin3(54) = 0.5118243051520034e0_wp
        xlin3(55) = 0.2069543381912151e0_wp
        xlin3(56) = -796.4010440692600e0_wp
        xlin3(57) = -4693.456970734772e0_wp
        xlin3(58) = 2597.246850944603e0_wp
        xlin3(59) = 630.4845628327919e0_wp
        xlin3(60) = -221.0729177223342e0_wp
        xlin3(61) = 9.839396652075672e0_wp
        xlin3(62) = -5227.790162633692e0_wp
        xlin3(63) = -529.4879950032479e0_wp
        xlin3(64) = 30.53701394166368e0_wp
        xlin3(65) = -14.22864486129778e0_wp
        xlin3(66) = -35.14731691232372e0_wp
        xlin3(67) = 2.140914121987492e0_wp
        xlin3(68) = 2250.367423678786e0_wp
        xlin3(69) = 138.5994085512685e0_wp
        xlin3(70) = 70.30416282155014e0_wp
        xlin3(71) = -62.74844809925895e0_wp
        xlin4(1)  = -1.2612782515309912e-2_wp
        xlin4(2)  = 1.5724446750473936e-2_wp
        xlin4(3)  = 0.2653887372495336e0_wp
        xlin4(4)  = 0.1216019078797281e0_wp
        xlin4(5)  = -0.6102844930670878e0_wp
        xlin4(6)  = 4.6420446436398755e-3_wp
        xlin4(7)  = 0.5563651571664677e0_wp
        xlin4(8)  = -0.3890198684837317e0_wp
        xlin4(9)  = 0.1052569138798340e0_wp
        xlin4(10) = -1.0080347996402828e-2_wp
        xlin4(11) = 0.2409741384095758e0_wp
        xlin4(12) = 5.745106082351191e0_wp
        xlin4(13) = 1.047295525278513e0_wp
        xlin4(14) = -12.03507114814484e0_wp
        xlin4(15) = 9.035186991905373e0_wp
        xlin4(16) = -2.267198360038944e0_wp
        xlin4(17) = 2.0836099694222385e-2_wp
        xlin4(18) = 4.7445938282087222e-2_wp
        xlin4(19) = -18.26537244396756e0_wp
        xlin4(20) = -4.809068965727743e0_wp
        xlin4(21) = 5.658647992475390e0_wp
        xlin4(22) = -1.284093508478602e0_wp
        xlin4(23) = -2.164733987255563e0_wp
        xlin4(24) = 0.3630073042947094e0_wp
        xlin4(25) = 3.760490137844927e0_wp
        xlin4(26) = -5.927286374087972e0_wp
        xlin4(27) = 0.1278839922444694e0_wp
        xlin4(28) = 0.5381095462111573e0_wp
        xlin4(29) = -2.179370013756247e0_wp
        xlin4(30) = 0.4127566904703107e0_wp
        xlin4(31) = 0.7168269469196163e0_wp
        xlin4(32) = -9.502943939736065e0_wp
        xlin4(33) = -9.597638311349975e0_wp
        xlin4(34) = 65.60860814240448e0_wp
        xlin4(35) = -56.04003254284753e0_wp
        xlin4(36) = 15.38423202570607e0_wp
        xlin4(37) = -0.3167914651468724e0_wp
        xlin4(38) = -0.7643266917502949e0_wp
        xlin4(39) = 0.1110705439245336e0_wp
        xlin4(40) = -146.5638710976931e0_wp
        xlin4(41) = 151.8079548350792e0_wp
        xlin4(42) = 15.86486188082157e0_wp
        xlin4(43) = -65.56160562330304e0_wp
        xlin4(44) = -1.126379669058077e0_wp
        xlin4(45) = -1.759058937005427e0_wp
        xlin4(46) = -3.1582037336166588e-4_wp
        xlin4(47) = 410.9371301349821e0_wp
        xlin4(48) = 142.2842073695382e0_wp
        xlin4(49) = -15.14809406881419e0_wp
        xlin4(50) = -9.858662824266736e0_wp
        xlin4(51) = 0.2622049408092991e0_wp
        xlin4(52) = -72.67984057864001e0_wp
        xlin4(53) = 11.28350883179339e0_wp
        xlin4(54) = -0.5323962789000201e0_wp
        xlin4(55) = -2.447540583751350e0_wp
        xlin4(56) = 148.9575982135905e0_wp
        xlin4(57) = 1249.939675860206e0_wp
        xlin4(58) = 268.6478693121586e0_wp
        xlin4(59) = 33.67689052007242e0_wp
        xlin4(60) = -35.65384689269180e0_wp
        xlin4(61) = 3.464021415311834e0_wp
        xlin4(62) = -103.4190791219546e0_wp
        xlin4(63) = 76.76961997409407e0_wp
        xlin4(64) = 11.59182442022548e0_wp
        xlin4(65) = -0.5254501536114633e0_wp
        xlin4(66) = -32.07157614897808e0_wp
        xlin4(67) = 1.102243506822152e0_wp
        xlin4(68) = -952.9757997893194e0_wp
        xlin4(69) = 52.11461967121677e0_wp
        xlin4(70) = 6.929202177070909e0_wp
        xlin4(71) = -6.322557284152793e0_wp
        epslon    = 1.e-12_wp
        epscem    = 1.e-12_wp
        ald       = 0.02e0_wp
        betad     = 0.72e0_wp
        maxi2     = 0
        maxi      = 0
        call indexa3(12, 12, iind, jind, kind, ipar, maxi, maxi2)
        maxi2a = maxi2
        maxi   = 0
        call indexa3(11, 11, iind, jind, kind, ipar, maxi, maxi2)
        maxi2b = maxi2-maxi2a
        call trip(r1x, etrip1, detrip1)
        call trip(r2x, etrip2, detrip2)
        call trip(r3x, etrip3, detrip3)
        call singlet(r1x, es1, des1)
        call singlet(r2x, es2, des2)
        call singlet(r3x, es3, des3)
        q1       = 0.5e0_wp*(es1+etrip1)
        q2       = 0.5e0_wp*(es2+etrip2)
        q3       = 0.5e0_wp*(es3+etrip3)
        xj1      = 0.5e0_wp*(es1-etrip1)
        xj2      = 0.5e0_wp*(es2-etrip2)
        xj3      = 0.5e0_wp*(es3-etrip3)
        xj       = sqrt(epslon+0.5e0_wp*((xj3-xj1)**2+(xj2-xj1)**2+(xj3-xj2)**2))
        vlondon  = q1+q2+q3-xj
        dq1      = 0.5e0_wp*(des1+detrip1)
        dq2      = 0.5e0_wp*(des2+detrip2)
        dq3      = 0.5e0_wp*(des3+detrip3)
        dxj1     = 0.5e0_wp*(des1-detrip1)
        dxj2     = 0.5e0_wp*(des2-detrip2)
        dxj3     = 0.5e0_wp*(des3-detrip3)
        dvlon(1) = dq1-0.5e0_wp*dxj1*(2.0e0_wp*xj1-xj2-xj3)/xj
        dvlon(2) = dq2-0.5e0_wp*dxj2*(2.0e0_wp*xj2-xj1-xj3)/xj
        dvlon(3) = dq3-0.5e0_wp*dxj3*(2.0e0_wp*xj3-xj1-xj2)/xj
        if( (r1x .le. betad).or.(r2x .le. betad).or.(r3x .le. betad) )then
            damper   = 0.e0_wp
            ddamp(1) = 0.e0_wp
            ddamp(2) = 0.e0_wp
            ddamp(3) = 0.e0_wp
        else
            rdamp1   = ald/(betad-r1x)
            rdamp2   = ald/(betad-r2x)
            rdamp3   = ald/(betad-r3x)
            damper   = exp(rdamp1+rdamp2+rdamp3)
            ddamp(1) = rdamp1/(betad-r1x)
            ddamp(2) = rdamp2/(betad-r2x)
            ddamp(3) = rdamp3/(betad-r3x)
        endif
        bb       = sqrt(epscem+(r1x-r3x)**2+(r1x-r2x)**2+(r2x-r3x)**2)
        dbb(1)   = (2.e0_wp*r1x-r2x-r3x)/bb
        dbb(2)   = (2.e0_wp*r2x-r1x-r3x)/bb
        dbb(3)   = (2.e0_wp*r3x-r1x-r2x)/bb
        warp     = 1.e0_wp/(1.e0_wp/r1x+1.e0_wp/r2x+1.e0_wp/r3x)
        dwarp(1) = (warp/r1x)**2
        dwarp(2) = (warp/r2x)**2
        dwarp(3) = (warp/r3x)**2
        rho1     = exp(-beta*(r1x-rnought))
        rho2     = exp(-beta*(r2x-rnought))
        rho3     = exp(-beta*(r3x-rnought))
        rho1t(0) = 1.e0_wp
        rho2t(0) = 1.e0_wp
        rho3t(0) = 1.e0_wp
        do it=1, 12
            rho1t(it) = rho1t(it-1)*rho1
            rho2t(it) = rho2t(it-1)*rho2
            rho3t(it) = rho3t(it-1)*rho3
        enddo
        sum1 = 0.e0_wp
        sum2 = 0.e0_wp
        sum3 = 0.e0_wp
        sum4 = 0.e0_wp
        do i = 1, 3
            dsum1(i) = 0.e0_wp
            dsum2(i) = 0.e0_wp
            dsum3(i) = 0.e0_wp
            dsum4(i) = 0.e0_wp
        enddo
        do ii=1, maxi2a
            rprod    = rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
            t1       = xlin1(ipar(ii))*rprod
            sum1     = sum1+t1
            dsum1(1) = dsum1(1)-beta*kind(ii)*t1
            dsum1(2) = dsum1(2)-beta*jind(ii)*t1
            dsum1(3) = dsum1(3)-beta*iind(ii)*t1
            t1       = xlin2(ipar(ii))*rprod
            sum2     = sum2+t1
            dsum2(1) = dsum2(1)-beta*kind(ii)*t1
            dsum2(2) = dsum2(2)-beta*jind(ii)*t1
            dsum2(3) = dsum2(3)-beta*iind(ii)*t1
        enddo
        do ii=maxi2a+1, maxi2a+maxi2b
            rprod    = rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
            t1       = xlin3(ipar(ii))*rprod
            sum3     = sum3+t1
            dsum3(1) = dsum3(1)-beta*kind(ii)*t1
            dsum3(2) = dsum3(2)-beta*jind(ii)*t1
            dsum3(3) = dsum3(3)-beta*iind(ii)*t1
            t1       = xlin4(ipar(ii))*rprod
            sum4     = sum4+t1
            dsum4(1) = dsum4(1)-beta*kind(ii)*t1
            dsum4(2) = dsum4(2)-beta*jind(ii)*t1
            dsum4(3) = dsum4(3)-beta*iind(ii)*t1
        enddo
        v3c = damper*(sum1+sum2*bb+warp*(sum3+sum4*bb))
        eee = vlondon+v3c+eshift
        do i = 1, 3
            dv(i)=dvlon(i) +ddamp(i)*v3c +damper*(dsum1(i)+dbb(i)*(sum2+warp*sum4)+bb*dsum2(i) +dwarp(i)*(sum3+bb* &
                & sum4)+ warp*(dsum3(i)+bb*dsum4(i)))
        enddo
    end subroutine dispatch
end module chempotpy_h3_gen_bh_2009_potcbs_legacy_state_b8d1df42

!> @brief Compatibility interface for potcbs.
subroutine potcbs( r1x, r2x, r3x, eee, dv )
    use chempotpy_h3_gen_bh_2009_potcbs_legacy_state_b8d1df42, only: dispatch
    use chempotpy_kinds,                                       only: wp
    implicit none
    real(wp) :: r1x
    real(wp) :: r2x
    real(wp) :: r3x
    real(wp) :: eee
    real(wp) :: dv(3)
    call dispatch(r1x=r1x, r2x=r2x, r3x=r3x, eee=eee, dv=dv)
end subroutine potcbs

!> @brief Evaluate indexa3.
subroutine indexa3( mbig, mtop, iind, jind, kind, ipar, maxi, maxi2 )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)    :: mbig, mtop
    integer, intent(inout) :: iind(5000), jind(5000), kind(5000), ipar(5000), maxi, maxi2
    integer                :: l, l2, i, j, k, isum
    l  = maxi2
    l2 = maxi
    do  i=0, min(mtop, mbig)
        do  j=i, min(mtop, mbig)
            do  k=j, min(mtop, mbig)
                isum = i+j+k
                if( isum .le. mbig )then
                    if(((isum-i).gt.0).and.((isum-j).gt.0).and.&
                        &        ((isum-k).gt.0))then
                    l2      = l2+1
                    l       = l+1
                    ipar(l) = l2
                    iind(l) = i
                    jind(l) = j
                    kind(l) = k
                    if( i .ne. j )then
                        if( j .ne. k )then
                            l       = l+1
                            ipar(l) = l2
                            iind(l) = i
                            jind(l) = k
                            kind(l) = j
                            l       = l+1
                            ipar(l) = l2
                            iind(l) = j
                            jind(l) = i
                            kind(l) = k
                            l       = l+1
                            ipar(l) = l2
                            iind(l) = j
                            jind(l) = k
                            kind(l) = i
                            l       = l+1
                            ipar(l) = l2
                            iind(l) = k
                            jind(l) = i
                            kind(l) = j
                            l       = l+1
                            ipar(l) = l2
                            iind(l) = k
                            jind(l) = j
                            kind(l) = i
                        else
                            l       = l+1
                            ipar(l) = l2
                            iind(l) = j
                            jind(l) = i
                            kind(l) = k
                            l       = l+1
                            ipar(l) = l2
                            iind(l) = j
                            jind(l) = k
                            kind(l) = i
                        endif
                    elseif( j .ne. k )then
                        l       = l+1
                        ipar(l) = l2
                        iind(l) = j
                        jind(l) = k
                        kind(l) = i
                        l       = l+1
                        ipar(l) = l2
                        iind(l) = k
                        jind(l) = j
                        kind(l) = i
                    endif
                endif
                endif
            enddo
        enddo
    enddo
    maxi  = l2
    maxi2 = l
end

!> @brief Evaluate trip.
subroutine trip( r, v, dv )
    use chempotpy_kinds, only: wp
    implicit none
    integer             :: j
    real(wp)            :: r, v, dv, t0, t1, t2, vlr, dvlr, vsr, dvsr
    real(wp)            :: damp, xd, xd2, prefac
    real(wp), parameter :: c6    = -6.499027e0_wp
    real(wp), parameter :: c8    = -124.3991e0_wp
    real(wp), parameter :: c10   = -3285.828e0_wp
    real(wp), parameter :: beta  = 0.2e0_wp
    real(wp), parameter :: re    = 1.401e0_wp
    real(wp), parameter :: alpha = 4.342902497495857e0_wp
    real(wp)            :: xlp(17)
    xlp(1)  = 2.190254292077946e-1_wp
    xlp(1)  = 2.190254292077946e-1_wp
    xlp(1)  = 2.190254292077946e-1_wp
    xlp(2)  = 6.903970886899540e-1_wp
    xlp(3)  = 1.163805164813647e+0_wp
    xlp(4)  = 1.337274224278977e+0_wp
    xlp(5)  = 1.176131501286896e+0_wp
    xlp(6)  = 9.093025316067311e-1_wp
    xlp(7)  = 5.908490740317524e-1_wp
    xlp(8)  = 1.467271248985162e-1_wp
    xlp(9)  = 6.218667627370255e-2_wp
    xlp(10) = 1.105747186790804e-1_wp
    xlp(11) = -7.112594221275334e-2_wp
    xlp(12) = -8.866929977503379e-3_wp
    xlp(13) = 3.823523348718375e-2_wp
    xlp(14) = -2.055915301253834e-2_wp
    xlp(15) = 5.419433399532862e-3_wp
    xlp(16) = -7.263483303294653e-4_wp
    xlp(17) = 4.154351972583585e-5_wp
    t0      = r-re
    prefac  = exp(-alpha*t0)
    if( (r-re) .ne. 0.e0_wp )then
        vsr  = 0.0
        dvsr = 0.0
        t1   = 1.0e0_wp
        t2   = 0.0e0_wp
        do j=1, 17
            vsr  = vsr+xlp(j)*t1*prefac
            dvsr = dvsr+(j-1)*xlp(j)*t2*prefac
            t2   = t1
            t1   = t1*t0
        enddo
    else
        vsr  = xlp(1)*prefac
        dvsr = xlp(2)*prefac
    endif
    dvsr = dvsr-alpha*vsr
    t1   = r*r
    t2   = exp(-beta*t1)
    damp = 1.e0_wp-t2
    xd   = damp/r
    xd2  = xd*xd
    vlr  = (c6+xd2*(c8+xd2*c10))*xd2**3
    dvlr = -(6.0e0_wp*c6+xd2*(8.0e0_wp*c8+xd2*10.0e0_wp*c10))* ((1/r)-2.0e0_wp*beta*t2/xd)*xd2**3
    v    = vlr+vsr
    dv   = dvlr+dvsr
end

!> @brief Evaluate singlet.
subroutine singlet( r, v, dv )
    use chempotpy_kinds, only: wp
    implicit none
    integer             :: j
    real(wp)            :: r, v, dv, t0, t1, t2, vsr, damp, xd, xd2, prefac
    real(wp), parameter :: c6    = -6.499027e0_wp
    real(wp), parameter :: c8    = -124.3991e0_wp
    real(wp), parameter :: c10   = -3285.828e0_wp
    real(wp), parameter :: beta  = 0.2e0_wp
    real(wp), parameter :: re    = 1.401e0_wp
    real(wp), parameter :: alpha = 3.980917850296971e0_wp
    real(wp)            :: xlp(17)
    xlp(1)  = -1.709663318869429e-1_wp
    xlp(2)  = -6.675286769482015e-1_wp
    xlp(3)  = -1.101876055072129e+0_wp
    xlp(4)  = -1.106460095658282e+0_wp
    xlp(5)  = -7.414724284525231e-1_wp
    xlp(6)  = -3.487882731923523e-1_wp
    xlp(7)  = -1.276736255756598e-1_wp
    xlp(8)  = -5.875965709151867e-2_wp
    xlp(9)  = -4.030128017933840e-2_wp
    xlp(10) = -2.038653237221007e-2_wp
    xlp(11) = -7.639198558462706e-4_wp
    xlp(12) = 2.912954920483885e-3_wp
    xlp(13) = -2.628273116815280e-4_wp
    xlp(14) = -4.622088855684211e-4_wp
    xlp(15) = 1.278468948126147e-4_wp
    xlp(16) = -1.157434070240206e-5_wp
    xlp(17) = -2.609691840882097e-12_wp
    t1      = r*r
    t2      = exp(-beta*t1)
    damp    = 1.e0_wp-t2
    xd      = damp/r
    xd2     = xd*xd
    v       = (c6+xd2*(c8+xd2*c10))*xd2**3
    dv      = -(6.0e0_wp*c6+xd2*(8.0e0_wp*c8+xd2*10.0e0_wp*c10))* ((1/r)-2.0e0_wp*beta*t2/xd)*xd2**3
    t0      = r-re
    prefac  = exp(-alpha*t0)
    vsr     = 0.0
    if( (r-re) .ne. 0.e0_wp )then
        t1 = 1.0e0_wp
        t2 = 0.0e0_wp
        do j=1, 17
            vsr = vsr+xlp(j)*t1*prefac
            dv  = dv+(j-1)*xlp(j)*t2*prefac
            t2  = t1
            t1  = t1*t0
        enddo
    else
        vsr = xlp(1)*prefac
        dv  = dv+xlp(2)*prefac
    endif
    v  = v+vsr
    dv = dv-alpha*vsr
end

!> @brief Evaluate golden.
subroutine golden( ax, bx, cx, tol, xmin, fval )
    use chempotpy_kinds, only: wp
    implicit none
    external fcn
    real(wp), parameter :: r = 0.61803399e0_wp, c = 1e0_wp-r
    real(wp)            :: f1, f2, x0, x1, x2, x3, tol
    real(wp)            :: fval, xmin, bx, ax, cx, fcn
    real(wp)            :: fa, fb, fc
    integer             :: icount
    integer, parameter  :: maxcount = 5000
    fa = fcn(ax)
    fb = fcn(bx)
    fc = fcn(cx)
    if( (fb .ge. fa).or.(fb .ge. fc) )then
        write(6, *)'bracketing error in golden'
        write(6, *)ax, fa
        write(6, *)bx, fb
        write(6, *)cx, fc
        stop
    endif
    x0 = ax
    x3 = cx
    if( abs(cx-bx) .gt. abs(bx-ax) )then
        x1 = bx
        f1 = fb
        x2 = bx+c*(cx-bx)
        f2 = fcn(x2)
    else
        x2 = bx
        f2 = fb
        x1 = bx-c*(bx-ax)
        f1 = fcn(x1)
    endif
    icount = 0
    do
        if( abs(x3-x0) .le. tol*(abs(x1)+abs(x2)) ) exit
        icount = icount+1
        if( icount .ge. maxcount )then
            write(6, *)'too many evaluations in golden'
            stop
        endif
        if( f2 .lt. f1 )then
            x0 = x1
            x1 = x2
            x2 = r*x1+c*x3
            f1 = f2
            f2 = fcn(x2)
        else
            x3 = x2
            x2 = x1
            x1 = r*x2+c*x0
            f2 = f1
            f1 = fcn(x1)
        endif
    enddo
    if( f1 .lt. f2 )then
        fval = f1
        xmin = x1
    else
        fval = f2
        xmin = x2
    endif
end

!> @brief Evaluate fcn.
real(wp) function fcn( x )
    use zmasses,         only: zmass1, zmass2
    use chempotpy_kinds, only: wp
    implicit none
    real(wp) :: x, e, bodc2body, v, dv, towave, xmh, massfact
    !      Evaluates the Born-Huang energy of an isomer of H2 with masses of zmass1 and zmass2
    xmh      = 1836.15264e0_wp   ! H nuclear mass (after substracting 1 for electron mass)
    massfact = 0.5e0_wp*xmh*(zmass1+zmass2)/(zmass1*zmass2)
    towave   = 219474.7e0_wp
    call singlet(x, v, dv)
    fcn = massfact*bodc2body(x, dv)/towave+v
end

!> @brief Implement H3_GEN_BH_2009, a potential energy surface for H3.
module chempotpy_h3_gen_bh_2009_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_BH_2009 potential energy surface for H3.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,              only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_units,              only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: triatomic_distances_to_bohr
        use chempotpy_kinds,              only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: r(3), dvp(3), v
        real(wp)              :: tx(9), r2(3), drdx(3,9)
        real(wp)              :: gx(nstates,9)
        integer               :: iatom, idir, i, j
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is HHH
        call triatomic_distances_to_bohr(x, r)
        r2 = r*bohr_to_angstrom
        call pot(r, v, dvp)
        p   = v*hartree_to_electronvolt
        dvp = dvp*hartree_per_bohr_to_ev_per_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(dvp, drdx, gx(1,:))
        do iatom=1, natoms
            do idir=1, 3
                j               = (iatom-1)*3+idir
                g(:,iatom,idir) = gx(:,j)
            enddo
        enddo
        if( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
    end subroutine evaluate
end module chempotpy_h3_gen_bh_2009_pes

!> @brief Expose H3_GEN_BH_2009 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_bh_2009_pes, only: evaluate
    use chempotpy_kinds,              only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    use chempotpy_kinds,    only: wp
    implicit none
    integer, parameter    :: first_atoms(3)  = [1, 2, 1]
    integer, parameter    :: second_atoms(3) = [2, 3, 3]
    real(wp), intent(in)  :: x(9), r(3)
    real(wp), intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx
