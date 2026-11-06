module zmasses
    implicit none
    private
    public :: zmass1, zmass2
    double precision :: zmass1, zmass2
end module

!> @brief Evaluate pot.
subroutine pot( rvp, evp, dvp )
    use zmasses, only: zmass1, zmass2
    implicit none
    double precision, parameter   :: ccishift = 0.1744759989649372d0
    double precision, parameter   :: towave   = 219474.7d0
    double precision, parameter   :: xmass1   = 1836.15264d0, xmass2 = 1836.15264d0, xmass3 = 1836.15264d0
    double precision, intent(in)  :: rvp(3)
    double precision, intent(out) :: evp, dvp(3)
    integer                       :: icase, info, i
    double precision              :: r1x, r2x, r3x, ecci
    double precision              :: dcci(3), d1(3), d2(3), d3(3)
    double precision              :: xmh, xmd, xmminv(3,3), amat(3), t
    double precision              :: e1, e2, e3, emat(3), ezero, hatom_bodc, xmin
    double precision              :: ezero2
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
    else if( xmass1 .eq. xmass3 )then
        zmass1 = xmass1
        zmass2 = xmass3
    else if( xmass2 .eq. xmass3 )then
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
    xmh        = 1836.15264d0   ! H nuclear mass (after substracting 1 for electron mass)
    xmd        = 3670.48293d0   ! D nuclear mass (after substracting 1 for electron mass)
    hatom_bodc = 59.6818d0  ! the BODC in cm-1 of H with an aug-cc-pVTZ basis set
    ezero      = hatom_bodc*xmh*(1d0/xmass1+1d0/xmass2+1d0/xmass3)
    call golden(1.3d0, 1.4d0, 1.5d0, 1.d-14, xmin, ezero2)
    ezero = ezero/towave+ezero2+ccishift
    if( (xmass1.eq.xmd) .and. (xmass2.eq.xmh) &
        &      .and. (xmass3.eq.xmh) )then
    icase = 1  ! We are solving for D + H2 so we don't need to solve for amat
else
    icase       = -1 ! We are not solving for D + H2 so we need 3 times as many evaluations so we can determine amat
    t           = (2.0d0*xmd+xmh)*(xmd-xmh)
    xmminv(:,:) = xmh*xmd*xmd/t
    t           = -xmh*xmd*(xmd+xmh)/t
    xmminv(1,1) = t
    xmminv(2,2) = t
    xmminv(3,3) = t
    endif
    ! R1 is the distance between mass1 and mass2, R2 is the distance between mass2 and mass3,
    ! and R3 is the distance between mass1 and mass3.  Note that the D+H2 code requires the homonuclear
    ! diatom distance to be R3, so before calling all the various subroutines we need to switch R2 and R3
    !
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
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r1x, r2x, r3x, eee, dv )
        implicit none
        optional :: r1x, r2x, r3x, eee, dv
        integer, parameter :: n2 = 5000
        dimension iind(n2), jind(n2), kind(n2), ipar(n2)
        integer          :: iind, jind, kind, ipar  ! arrays of length n2
        integer          :: ii, it, maxi, maxi2, maxi2a, i
        double precision :: r1x, r2x, r3x, rho1, rho2, rho3, eee, dv(:)
        double precision :: drho1, drho2, drho3
        double precision :: e1body, xmd, xmh, hatom_bodc, bodc2body
        double precision :: round, ald, betad, beta, deldamp, rhodamp
        double precision :: bb2, bb3, bb4, bb5, v, rnought, towave
        double precision :: xma, xmb, xmc, v1, v2, gtot, rsum, &
            &                    betadsum, rdampsum, v3c
        double precision :: sum1, sum2, sum3, sum4, sum5, rprod
        double precision :: bigq, slit, slitsq, damper, rpass(3), &
            &                    sdamp, rdamp1, rdamp2, rdamp3
        double precision :: slit1, slit2, slit3, ddamper(1:3), dsdamp
        double precision :: dgtot(1:3), dbb2(1:3), dbb3(1:3)
        double precision :: dv1, dv2, dv3
        double precision :: dsum1(1:3), dsum2(1:3), dsum3(1:3)
        double precision :: rho1t(0:13), rho2t(0:13), rho3t(0:13)
        double precision :: drho1t(0:13), drho2t(0:13), drho3t(0:13)
        double precision :: drprod1, drprod2, drprod3
        double precision :: xlin1(180), xlin2(180), xlin3(180)
        ald        = 0.d0
        bb2        = 0.d0
        bb3        = 0.d0
        bb4        = 0.d0
        bb5        = 0.d0
        beta       = 0.d0
        betad      = 0.d0
        betadsum   = 0.d0
        bigq       = 0.d0
        damper     = 0.d0
        dbb2       = 0.d0
        dbb3       = 0.d0
        ddamper    = 0.d0
        deldamp    = 0.d0
        dgtot      = 0.d0
        drho1      = 0.d0
        drho1t     = 0.d0
        drho2      = 0.d0
        drho2t     = 0.d0
        drho3      = 0.d0
        drho3t     = 0.d0
        drprod1    = 0.d0
        drprod2    = 0.d0
        drprod3    = 0.d0
        dsdamp     = 0.d0
        dsum1      = 0.d0
        dsum2      = 0.d0
        dsum3      = 0.d0
        dv1        = 0.d0
        dv2        = 0.d0
        dv3        = 0.d0
        e1body     = 0.d0
        gtot       = 0.d0
        hatom_bodc = 0.d0
        rdamp1     = 0.d0
        rdamp2     = 0.d0
        rdamp3     = 0.d0
        rdampsum   = 0.d0
        rho1       = 0.d0
        rho1t      = 0.d0
        rho2       = 0.d0
        rho2t      = 0.d0
        rho3       = 0.d0
        rho3t      = 0.d0
        rhodamp    = 0.d0
        rnought    = 0.d0
        round      = 0.d0
        rpass      = 0.d0
        rprod      = 0.d0
        rsum       = 0.d0
        sdamp      = 0.d0
        slit       = 0.d0
        slit1      = 0.d0
        slit2      = 0.d0
        slit3      = 0.d0
        slitsq     = 0.d0
        sum1       = 0.d0
        sum2       = 0.d0
        sum3       = 0.d0
        sum4       = 0.d0
        sum5       = 0.d0
        towave     = 0.d0
        v          = 0.d0
        v1         = 0.d0
        v2         = 0.d0
        v3c        = 0.d0
        xlin1      = 0.d0
        xlin2      = 0.d0
        xlin3      = 0.d0
        xma        = 0.d0
        xmb        = 0.d0
        xmc        = 0.d0
        xmd        = 0.d0
        xmh        = 0.d0
        xlin1(1)   = 2.451839720665176d+02
        xlin1(2)   = -3.023514674602468d+02
        xlin1(3)   = -3.643179260300581d+02
        xlin1(4)   = 1.033283164173005d+03
        xlin1(5)   = -9.649202796347731d+02
        xlin1(6)   = 4.723249503697386d+02
        xlin1(7)   = -1.312024315948558d+02
        xlin1(8)   = 2.053242829600808d+01
        xlin1(9)   = -1.663116010000327d+00
        xlin1(10)  = 5.317955346689944d-02
        xlin1(11)  = 3.693835631800750d+02
        xlin1(12)  = 9.643897029819937d+01
        xlin1(13)  = -3.794338103064540d+02
        xlin1(14)  = 2.184383291329870d+02
        xlin1(15)  = -3.787026565563110d+01
        xlin1(16)  = -3.960932914729867d+00
        xlin1(17)  = 1.733261249407424d+00
        xlin1(18)  = -1.233088325428533d-01
        xlin1(19)  = 7.634113468160084d+02
        xlin1(20)  = -6.089390321694615d+02
        xlin1(21)  = 2.145121698941533d+02
        xlin1(22)  = -4.303659182480859d+01
        xlin1(23)  = 4.793100434454963d+00
        xlin1(24)  = -2.362075732031712d-01
        xlin1(25)  = 5.161649762181137d+02
        xlin1(26)  = -1.379218427454547d+02
        xlin1(27)  = 1.222214446630326d+01
        xlin1(28)  = -1.700714347704069d-01
        xlin1(29)  = 2.590291319801722d+01
        xlin1(30)  = -1.187599837006957d+00
        xlin1(31)  = 2.949665817120014d+02
        xlin1(32)  = -4.423967394072385d+02
        xlin1(33)  = -1.499691959357717d+02
        xlin1(34)  = 8.615959233030945d+02
        xlin1(35)  = -9.063992109731515d+02
        xlin1(36)  = 4.800278837250469d+02
        xlin1(37)  = -1.444131347943261d+02
        xlin1(38)  = 2.489091872075205d+01
        xlin1(39)  = -2.288354420765698d+00
        xlin1(40)  = 8.706064188372947d-02
        xlin1(41)  = -1.561653042358118d+03
        xlin1(42)  = 2.795630390969673d+03
        xlin1(43)  = -4.442497231862677d+03
        xlin1(44)  = 2.653340911293203d+03
        xlin1(45)  = -7.274300003583212d+02
        xlin1(46)  = 1.223835539827257d+02
        xlin1(47)  = -1.767422177777780d+01
        xlin1(48)  = 1.080509007725006d+00
        xlin1(49)  = 8.007412890755861d-02
        xlin1(50)  = -6.866852114392939d+02
        xlin1(51)  = 1.965474846931427d+03
        xlin1(52)  = -8.327550428277026d+02
        xlin1(53)  = 2.921969416962440d+00
        xlin1(54)  = 5.934336236873139d+01
        xlin1(55)  = -1.157242040507434d+01
        xlin1(56)  = 6.777877552427620d-01
        xlin1(57)  = 6.809701703798451d+02
        xlin1(58)  = -6.460103807549953d+02
        xlin1(59)  = 1.080711387406044d+02
        xlin1(60)  = -2.581355973993933d+00
        xlin1(61)  = -2.702457089997727d-01
        xlin1(62)  = 3.373040929271790d+02
        xlin1(63)  = -3.319320528532528d+01
        xlin1(64)  = -8.252443297418132d-01
        xlin1(65)  = 2.242813487013438d+00
        xlin1(66)  = -4.786033544785251d+02
        xlin1(67)  = 7.379465395733814d+02
        xlin1(68)  = -2.997606193372937d+02
        xlin1(69)  = -1.444004106661323d+02
        xlin1(70)  = 1.160834932238466d+02
        xlin1(71)  = -1.055041157003008d+01
        xlin1(72)  = -6.084935176995607d+00
        xlin1(73)  = 1.255348914553237d+00
        xlin1(74)  = -4.967414449352265d-02
        xlin1(75)  = 2.869805514805477d+03
        xlin1(76)  = -2.901158537312932d+02
        xlin1(77)  = 2.593678864031728d+03
        xlin1(78)  = -1.474577262776355d+03
        xlin1(79)  = 1.739661050571814d+02
        xlin1(80)  = 4.108033331334562d+01
        xlin1(81)  = -9.031479084332030d+00
        xlin1(82)  = 3.265756035116496d-01
        xlin1(83)  = -9.394389739950026d+03
        xlin1(84)  = 1.281273536289386d+03
        xlin1(85)  = 9.186863730169819d+02
        xlin1(86)  = -2.224206878555418d+02
        xlin1(87)  = 2.852507226985430d+01
        xlin1(88)  = -2.975868406356395d+00
        xlin1(89)  = -2.887753159216629d+02
        xlin1(90)  = -6.232442825355130d+01
        xlin1(91)  = -2.080038525889382d+01
        xlin1(92)  = 5.846990317948459d+00
        xlin1(93)  = 1.494611359198648d+01
        xlin1(94)  = 1.430793238496208d+00
        xlin1(95)  = -1.161743802353405d+02
        xlin1(96)  = -4.435211957384280d+02
        xlin1(97)  = 7.949179810694009d+02
        xlin1(98)  = -4.621878870840682d+02
        xlin1(99)  = 1.829448335674127d+02
        xlin1(100) = -4.859832146779939d+01
        xlin1(101) = 6.467586280687250d+00
        xlin1(102) = -3.014121429251370d-01
        xlin1(103) = -4.727165211936207d+03
        xlin1(104) = 1.681945178807448d+03
        xlin1(105) = 6.001879815346814d+02
        xlin1(106) = -5.284816509082527d+02
        xlin1(107) = 9.509901063900256d+01
        xlin1(108) = -8.355776044216967d+00
        xlin1(109) = 4.667856406455721d-01
        xlin1(110) = 5.648107714717939d+02
        xlin1(111) = -2.902912088460810d+02
        xlin1(112) = -9.407789952750028d+01
        xlin1(113) = -1.631433254869286d+00
        xlin1(114) = 5.242491229250340d+00
        xlin1(115) = 3.214464002696617d+01
        xlin1(116) = 2.913754508725327d+01
        xlin1(117) = -1.568292873005413d+01
        xlin1(118) = 5.462484623657479d+00
        xlin1(119) = 9.083635945770520d+02
        xlin1(120) = 9.799890386730289d+01
        xlin1(121) = -6.411746299039472d+02
        xlin1(122) = 4.298716219171475d+02
        xlin1(123) = -1.205657322203302d+02
        xlin1(124) = 1.315814978109771d+01
        xlin1(125) = -3.975978827444066d-01
        xlin1(126) = 3.180124260914723d+03
        xlin1(127) = -5.872907359889355d+02
        xlin1(128) = -6.259492461216178d+02
        xlin1(129) = 3.434767513553659d+02
        xlin1(130) = -3.640741128738960d+01
        xlin1(131) = -9.236985549909328d-01
        xlin1(132) = 1.155215128322806d+03
        xlin1(133) = -9.784160598115349d+01
        xlin1(134) = 1.777745091732398d+01
        xlin1(135) = 2.485226250947354d-01
        xlin1(136) = 1.521755613790718d+01
        xlin1(137) = 7.156889276881836d+00
        xlin1(138) = -1.016142747224836d+03
        xlin1(139) = -5.200008131651686d+01
        xlin1(140) = 2.946853566528521d+02
        xlin1(141) = -1.316509239684804d+02
        xlin1(142) = 2.366837195107599d+01
        xlin1(143) = -1.174923797943324d+00
        xlin1(144) = -1.159057124509925d+03
        xlin1(145) = -1.521893512764090d+02
        xlin1(146) = 1.104768477989274d+02
        xlin1(147) = -3.656291398178829d+01
        xlin1(148) = 3.488333934608939d+00
        xlin1(149) = -2.286294041696426d+02
        xlin1(150) = -7.953044862181258d+00
        xlin1(151) = -2.999788001060466d-02
        xlin1(152) = -1.546533722850765d+01
        xlin1(153) = 5.670669717065498d+02
        xlin1(154) = 5.041037222513963d+01
        xlin1(155) = -7.690800748197464d+01
        xlin1(156) = 1.401968534663515d+01
        xlin1(157) = -1.189136309799305d+00
        xlin1(158) = 2.944715172811239d+02
        xlin1(159) = 1.011155288939896d+02
        xlin1(160) = -5.420733436869213d+00
        xlin1(161) = -8.900491551445843d-01
        xlin1(162) = 2.785072259806540d+01
        xlin1(163) = 6.092864636803840d+00
        xlin1(164) = -1.799957321893311d+02
        xlin1(165) = -1.892894466520360d+01
        xlin1(166) = 1.036160665441249d+01
        xlin1(167) = -3.761882545122726d-01
        xlin1(168) = -5.436396122931865d+01
        xlin1(169) = -1.702038609057893d+01
        xlin1(170) = 4.969068890501028d-03
        xlin1(171) = -3.519006024232350d+00
        xlin1(172) = 3.284752951627190d+01
        xlin1(173) = 2.783127622573391d+00
        xlin1(174) = -5.624072554917517d-01
        xlin1(175) = 5.370237229208703d+00
        xlin1(176) = 9.768755790120242d-01
        xlin1(177) = -3.206502231938864d+00
        xlin1(178) = -1.267702449737564d-01
        xlin1(179) = -1.430186406716300d-01
        xlin1(180) = 1.296229968602821d-01
        xlin2(1)   = 8.145631479853698d-01
        xlin2(2)   = 1.770730113974052d+01
        xlin2(3)   = -2.173959339189124d+01
        xlin2(4)   = 1.961131724453491d+01
        xlin2(5)   = -4.104999264119586d+01
        xlin2(6)   = 2.337793724047610d+01
        xlin2(7)   = -5.449082199443552d+00
        xlin2(8)   = 6.201541014865459d-01
        xlin2(9)   = -4.728365816655994d-02
        xlin2(10)  = 3.321586948535380d-03
        xlin2(11)  = -6.652255446067896d+01
        xlin2(12)  = 5.962375314289732d+01
        xlin2(13)  = 1.469998744381738d+01
        xlin2(14)  = 1.153758314494926d+00
        xlin2(15)  = -7.591897754372836d+00
        xlin2(16)  = 2.103874186533575d+00
        xlin2(17)  = -1.702965824754904d-01
        xlin2(18)  = 1.346878024842888d-04
        xlin2(19)  = -1.178260275080853d+02
        xlin2(20)  = 3.489556535478315d+01
        xlin2(21)  = -1.251243111695901d+00
        xlin2(22)  = 2.471519774632120d-01
        xlin2(23)  = -2.474793708434775d-01
        xlin2(24)  = 2.687242988983517d-02
        xlin2(25)  = -2.026760191030008d+01
        xlin2(26)  = 2.756556273693077d+00
        xlin2(27)  = 8.357078625422904d-02
        xlin2(28)  = -1.640242746249622d-02
        xlin2(29)  = -6.876043147213573d-01
        xlin2(30)  = 2.609337260717558d-02
        xlin2(31)  = 3.002603577505805d+00
        xlin2(32)  = 1.498828705957069d+01
        xlin2(33)  = -2.086306249380632d+01
        xlin2(34)  = 2.081499551542078d+01
        xlin2(35)  = -3.171758822478321d+01
        xlin2(36)  = 1.723622638287100d+01
        xlin2(37)  = -4.028169514216751d+00
        xlin2(38)  = 4.671596528203490d-01
        xlin2(39)  = -3.670331530469922d-02
        xlin2(40)  = 2.629265006741015d-03
        xlin2(41)  = -8.927371067901906d+01
        xlin2(42)  = 8.867171901574056d+01
        xlin2(43)  = -1.133027069526800d+02
        xlin2(44)  = 1.602734361579233d+02
        xlin2(45)  = -1.023650687832928d+02
        xlin2(46)  = 2.471829619377287d+01
        xlin2(47)  = -2.300433949989794d+00
        xlin2(48)  = 2.592872235150836d-01
        xlin2(49)  = -4.379218716825992d-02
        xlin2(50)  = -1.231186961261672d+01
        xlin2(51)  = 1.420549968400411d+01
        xlin2(52)  = 1.232687494983029d+00
        xlin2(53)  = 1.841045801282512d+01
        xlin2(54)  = -6.141054432632874d+00
        xlin2(55)  = 2.034847145705890d-01
        xlin2(56)  = 8.860867509074010d-02
        xlin2(57)  = -9.644503856306309d+01
        xlin2(58)  = 2.392130829529128d+01
        xlin2(59)  = -3.146354866365471d+00
        xlin2(60)  = 9.562642394369082d-01
        xlin2(61)  = -2.269962411748502d-01
        xlin2(62)  = -7.627358272973419d+00
        xlin2(63)  = 4.570619067451160d-01
        xlin2(64)  = 2.237747868045461d-01
        xlin2(65)  = -3.169396242810660d-01
        xlin2(66)  = 1.231481888044332d+01
        xlin2(67)  = -5.861438736123795d+01
        xlin2(68)  = 5.294643266448595d+01
        xlin2(69)  = -4.088837194494829d+00
        xlin2(70)  = 5.495540512301918d+00
        xlin2(71)  = -6.300138454926564d+00
        xlin2(72)  = 1.765098522094756d+00
        xlin2(73)  = -1.780074643720451d-01
        xlin2(74)  = 4.114464551103923d-03
        xlin2(75)  = 8.917318947290573d+01
        xlin2(76)  = -1.770180683611587d+01
        xlin2(77)  = -5.602678714811219d+00
        xlin2(78)  = 2.973755372611144d+01
        xlin2(79)  = 5.605823953149696d+00
        xlin2(80)  = -4.949227356210042d+00
        xlin2(81)  = 2.420977258547833d-01
        xlin2(82)  = 8.167453315236323d-02
        xlin2(83)  = -1.135834657778073d+02
        xlin2(84)  = 8.290397544044251d+01
        xlin2(85)  = -7.296380654367582d+01
        xlin2(86)  = 1.347009356803447d+01
        xlin2(87)  = -2.500870639262219d-01
        xlin2(88)  = -1.272556171177672d-01
        xlin2(89)  = 5.669201615035288d+00
        xlin2(90)  = 4.449130216981585d+00
        xlin2(91)  = -1.272014612117903d+00
        xlin2(92)  = 2.139168528456700d-01
        xlin2(93)  = -4.777213854750650d-02
        xlin2(94)  = -1.275863250945433d-02
        xlin2(95)  = -1.642934979019225d+01
        xlin2(96)  = 5.594030291218376d+01
        xlin2(97)  = -8.305010271115899d+01
        xlin2(98)  = 3.140219381847346d+01
        xlin2(99)  = -2.534163029719606d+00
        xlin2(100) = 1.142479699516450d-01
        xlin2(101) = -1.263578178350560d-01
        xlin2(102) = 1.629095589290311d-02
        xlin2(103) = -1.136985820244078d+02
        xlin2(104) = 3.262096657063820d+01
        xlin2(105) = -1.087121813266270d+02
        xlin2(106) = 2.673092202326776d+01
        xlin2(107) = -5.444160084123269d-01
        xlin2(108) = 6.138628773543448d-01
        xlin2(109) = -2.227712501792724d-01
        xlin2(110) = 5.512545979043862d+01
        xlin2(111) = 2.559974311865491d+01
        xlin2(112) = 5.216798906963263d-01
        xlin2(113) = -1.496009262806457d+00
        xlin2(114) = 2.638941938405173d-01
        xlin2(115) = -9.723650088603947d+00
        xlin2(116) = 6.868905539663074d-01
        xlin2(117) = -4.177757286042882d-01
        xlin2(118) = 2.781285313472535d-01
        xlin2(119) = 1.563554076363040d+01
        xlin2(120) = -6.542137033344432d+00
        xlin2(121) = 2.193791357437970d+01
        xlin2(122) = -1.576158545476312d+01
        xlin2(123) = 2.587945022189958d+00
        xlin2(124) = -3.893301430803128d-02
        xlin2(125) = -5.190707866429820d-03
        xlin2(126) = 1.367205185988990d+02
        xlin2(127) = 1.188103897312893d+01
        xlin2(128) = 1.606164459906229d+01
        xlin2(129) = -9.262867579726493d+00
        xlin2(130) = 6.271894032241021d-01
        xlin2(131) = 2.285179285451767d-01
        xlin2(132) = -7.103397905922931d+01
        xlin2(133) = 3.121061890947882d+00
        xlin2(134) = 6.354928967335751d-01
        xlin2(135) = -4.163368478300122d-02
        xlin2(136) = 1.020366559611126d+00
        xlin2(137) = 1.837284088182978d-01
        xlin2(138) = -2.832630010141147d+01
        xlin2(139) = 1.275654589067907d+01
        xlin2(140) = -4.603547842975937d-01
        xlin2(141) = 1.995827205765568d+00
        xlin2(142) = -5.232664167379651d-01
        xlin2(143) = 1.710301952383997d-02
        xlin2(144) = -8.590072407323460d+01
        xlin2(145) = 1.674416687913584d+01
        xlin2(146) = 1.953210020027019d-02
        xlin2(147) = 2.520910415921704d-01
        xlin2(148) = -3.423339521867145d-01
        xlin2(149) = 1.169816725323726d+01
        xlin2(150) = -1.761561723165250d+00
        xlin2(151) = 3.013931319831423d-02
        xlin2(152) = -4.165165498393316d-01
        xlin2(153) = 1.386201564606759d+01
        xlin2(154) = -9.091983474080589d+00
        xlin2(155) = 2.768047770222391d-01
        xlin2(156) = 5.763259032772141d-02
        xlin2(157) = 3.029683243309628d-02
        xlin2(158) = 1.740587394078216d+01
        xlin2(159) = -6.199982916108434d+00
        xlin2(160) = 9.489522270750911d-01
        xlin2(161) = 2.535531600295191d-01
        xlin2(162) = 1.698179626605738d-01
        xlin2(163) = 2.305666871236848d-01
        xlin2(164) = -2.509462487790768d+00
        xlin2(165) = 2.075944787757823d+00
        xlin2(166) = -2.275822797369400d-01
        xlin2(167) = -2.065529414930606d-02
        xlin2(168) = -1.074236086914455d+00
        xlin2(169) = 6.803905517919883d-02
        xlin2(170) = -2.615499139877834d-01
        xlin2(171) = -1.335426493176628d-01
        xlin2(172) = 1.840487213722926d-01
        xlin2(173) = -1.569973762062688d-01
        xlin2(174) = 3.050419120227986d-02
        xlin2(175) = 2.653719569186909d-01
        xlin2(176) = 1.094965343816732d-01
        xlin2(177) = -2.187339773412837d-02
        xlin2(178) = -1.873849846674096d-03
        xlin2(179) = -5.603548924512692d-02
        xlin2(180) = 3.558364313199398d-03
        xlin3(1)   = -7.638557823809794d+01
        xlin3(2)   = 1.703707427310439d+02
        xlin3(3)   = -1.025601533591727d+02
        xlin3(4)   = -5.296578652461962d+01
        xlin3(5)   = 1.093595183012879d+02
        xlin3(6)   = -6.588531176536077d+01
        xlin3(7)   = 2.040793884905717d+01
        xlin3(8)   = -3.479631418111812d+00
        xlin3(9)   = 3.088169472711390d-01
        xlin3(10)  = -1.111590540608995d-02
        xlin3(11)  = -3.361962954982916d+02
        xlin3(12)  = 2.446590342879648d+02
        xlin3(13)  = -6.622902002288096d+01
        xlin3(14)  = 1.191674304065106d+01
        xlin3(15)  = -1.201463499030519d+01
        xlin3(16)  = 6.122631636838099d+00
        xlin3(17)  = -1.187879235377073d+00
        xlin3(18)  = 7.926806077855969d-02
        xlin3(19)  = -5.725326950390917d+01
        xlin3(20)  = 2.868589802572544d+01
        xlin3(21)  = -2.971320475494586d+01
        xlin3(22)  = 1.008864240170261d+01
        xlin3(23)  = -1.365325434437133d+00
        xlin3(24)  = 6.554733610556563d-02
        xlin3(25)  = -4.045186273038085d+01
        xlin3(26)  = 1.874648094047193d+01
        xlin3(27)  = -2.454684686634876d+00
        xlin3(28)  = 6.765838416051127d-02
        xlin3(29)  = -4.059257274704468d+00
        xlin3(30)  = 2.165639317481252d-01
        xlin3(31)  = -7.534329078823575d+01
        xlin3(32)  = 1.796639204536162d+02
        xlin3(33)  = -1.564733447866279d+02
        xlin3(34)  = 3.775541937573816d+01
        xlin3(35)  = 3.375720237731807d+01
        xlin3(36)  = -3.005181622792514d+01
        xlin3(37)  = 1.029486036217851d+01
        xlin3(38)  = -1.791612509905809d+00
        xlin3(39)  = 1.538953705090221d-01
        xlin3(40)  = -5.019936328398713d-03
        xlin3(41)  = -4.774634443417375d+01
        xlin3(42)  = 3.863906310035196d+02
        xlin3(43)  = 1.909693601494610d+00
        xlin3(44)  = 7.896215873407198d+01
        xlin3(45)  = -1.295233067421028d+02
        xlin3(46)  = 4.285454822181649d+01
        xlin3(47)  = -3.947563931965661d+00
        xlin3(48)  = -2.453031697183805d-02
        xlin3(49)  = -4.072260116080702d-03
        xlin3(50)  = -1.512501065486055d+03
        xlin3(51)  = 9.514665008404671d+02
        xlin3(52)  = -6.253456734186879d+02
        xlin3(53)  = 2.584487207114875d+02
        xlin3(54)  = -5.147932463240105d+01
        xlin3(55)  = 4.978444306783625d+00
        xlin3(56)  = -2.242912909135229d-01
        xlin3(57)  = -5.814167788685131d+02
        xlin3(58)  = 2.286668120955779d+02
        xlin3(59)  = -6.032832571577611d+01
        xlin3(60)  = 8.065674990175529d+00
        xlin3(61)  = -3.404150604733470d-01
        xlin3(62)  = -3.061436485306763d+01
        xlin3(63)  = 2.200637753633254d+00
        xlin3(64)  = 1.057897733012469d-01
        xlin3(65)  = 7.289136325949297d-03
        xlin3(66)  = 1.896105176076998d+02
        xlin3(67)  = -3.648470785265107d+02
        xlin3(68)  = 2.595841400715867d+02
        xlin3(69)  = -6.319152600268254d+01
        xlin3(70)  = 8.494896397888191d+00
        xlin3(71)  = -1.053255418482307d+01
        xlin3(72)  = 5.467310061893341d+00
        xlin3(73)  = -1.031401599859863d+00
        xlin3(74)  = 6.598583555227435d-02
        xlin3(75)  = 4.919698991437962d+02
        xlin3(76)  = -1.523587257667702d+03
        xlin3(77)  = 1.014167387686189d+03
        xlin3(78)  = -6.508692392742124d+02
        xlin3(79)  = 2.713378987371789d+02
        xlin3(80)  = -5.735977155603892d+01
        xlin3(81)  = 5.830582807234884d+00
        xlin3(82)  = -2.440949592570337d-01
        xlin3(83)  = -6.735007851702401d+02
        xlin3(84)  = 7.901400404656945d+02
        xlin3(85)  = -1.229867349843578d+02
        xlin3(86)  = -3.009484917235691d+01
        xlin3(87)  = 2.526743872042880d+00
        xlin3(88)  = 5.970110295569673d-01
        xlin3(89)  = -3.532469213851834d+02
        xlin3(90)  = 4.005928931048262d+01
        xlin3(91)  = 5.976702173098687d+00
        xlin3(92)  = -1.614616921246194d+00
        xlin3(93)  = -9.871235251782917d+00
        xlin3(94)  = 3.290059681876524d-01
        xlin3(95)  = -1.813749415305320d+02
        xlin3(96)  = 3.109697012254384d+02
        xlin3(97)  = -3.205683056506578d+01
        xlin3(98)  = -3.657410087464648d+00
        xlin3(99)  = -2.624482850677591d+01
        xlin3(100) = 1.173325195979073d+01
        xlin3(101) = -1.590994311096534d+00
        xlin3(102) = 5.486195466674193d-02
        xlin3(103) = -2.440497562681578d+02
        xlin3(104) = 1.038109037697877d+03
        xlin3(105) = -7.152807104622796d+02
        xlin3(106) = 2.980486621854296d+02
        xlin3(107) = -7.381873725420134d+01
        xlin3(108) = 9.536084102100837d+00
        xlin3(109) = -4.034487399636460d-01
        xlin3(110) = 6.749359685299933d+02
        xlin3(111) = -3.958160095292775d+02
        xlin3(112) = 7.043339305226182d+01
        xlin3(113) = 4.280031318910807d+00
        xlin3(114) = -1.799431596685355d+00
        xlin3(115) = 1.572684765203245d+02
        xlin3(116) = -1.497852043495127d+01
        xlin3(117) = 3.450483851677588d+00
        xlin3(118) = -1.651601456234870d+00
        xlin3(119) = 6.267522188135509d+01
        xlin3(120) = -1.490721159969828d+02
        xlin3(121) = 7.203117438959657d+01
        xlin3(122) = -4.910294960766687d+01
        xlin3(123) = 2.126436044051092d+01
        xlin3(124) = -3.286263182633637d+00
        xlin3(125) = 1.432920523463851d-01
        xlin3(126) = 3.482161123451701d+02
        xlin3(127) = -6.634624347379950d+02
        xlin3(128) = 2.978996288910582d+02
        xlin3(129) = -6.102181144139994d+01
        xlin3(130) = 4.913896432859977d+00
        xlin3(131) = 9.220877926952298d-02
        xlin3(132) = -9.040330715822043d+01
        xlin3(133) = 3.763878231734758d+01
        xlin3(134) = -1.466110603121045d+01
        xlin3(135) = 4.286621953679466d-01
        xlin3(136) = -1.171184819157715d+01
        xlin3(137) = -8.282936408674097d-01
        xlin3(138) = 2.329172506772613d+01
        xlin3(139) = 7.205432276781195d+01
        xlin3(140) = -7.134506264990834d+01
        xlin3(141) = 2.778750878256592d+01
        xlin3(142) = -4.744233867644683d+00
        xlin3(143) = 2.489849788505943d-01
        xlin3(144) = -3.035952094511171d+02
        xlin3(145) = 2.774672648114976d+02
        xlin3(146) = -7.050154798680367d+01
        xlin3(147) = 7.641437548457593d+00
        xlin3(148) = -1.880492973015346d-01
        xlin3(149) = -1.984086397967394d+01
        xlin3(150) = 7.604978254357780d+00
        xlin3(151) = 2.014923057356265d-01
        xlin3(152) = 3.308108128880802d+00
        xlin3(153) = -2.969993736353019d+01
        xlin3(154) = -3.682715000434093d+01
        xlin3(155) = 2.210331765988686d+01
        xlin3(156) = -4.301393898599036d+00
        xlin3(157) = 2.671370386619945d-01
        xlin3(158) = 1.124784775141184d+02
        xlin3(159) = -6.041199428997342d+01
        xlin3(160) = 6.988933983482634d+00
        xlin3(161) = -1.362655996129701d-01
        xlin3(162) = -2.997626634973072d-02
        xlin3(163) = -1.846855278627695d+00
        xlin3(164) = 1.154319604248443d+01
        xlin3(165) = 1.192127579658776d+01
        xlin3(166) = -2.699670781911932d+00
        xlin3(167) = 1.920744268281359d-01
        xlin3(168) = -1.990490891239322d+01
        xlin3(169) = 7.011127199328610d+00
        xlin3(170) = -9.502550138734989d-02
        xlin3(171) = 7.777433674389070d-01
        xlin3(172) = -2.235484930636709d+00
        xlin3(173) = -1.909433565086661d+00
        xlin3(174) = 1.029497435932162d-01
        xlin3(175) = 1.784299385800847d+00
        xlin3(176) = -3.937636905107279d-01
        xlin3(177) = 2.159590396303356d-01
        xlin3(178) = 1.161199705458321d-01
        xlin3(179) = -7.827610505544675d-02
        xlin3(180) = -8.180052967228596d-03
        !     distances are in Bohr, energies are in cm-1
        towave     = 219474.7d0  ! hartree to cm-1 conversion
        xmh        = 1836.15264d0   ! H nuclear mass (after substracting 1 for electron mass)
        xmd        = 3670.48293d0   ! D nuclear mass (after substracting 1 for electron mass)
        hatom_bodc = 59.6818d0  ! the BODC in cm-1 of H with an aug-cc-pVTZ basis set
        xma        = xmd
        xmb        = xmh
        xmc        = xmh
        round      = 1.d-2
        ald        = 0.02d0
        betad      = 0.72d0
        betadsum   = 3.3d0
        deldamp    = 0.25d0
        rhodamp    = 0.01d0
        beta       = 1.414254132358239d0
        rnought    = 2.0d0
        maxi2      = 0
        maxi       = 0
        call indexa2(11, 11, iind, jind, kind, ipar, maxi, maxi2)
        maxi2a = maxi2
        rsum   = r1x+r2x+r3x
        if((r1x.le.betad).or.(r2x.le.betad).or.(r3x.le.betad)&
            &    .or.(rsum.le.betadsum) )then
        damper     = 0.0d0
        ddamper(:) = 0.0d0
    else
        bigq  = r1x*r1x+r2x*r2x+r3x*r3x
        slit1 = 2.0d0*r1x**2-r2x**2-r3x**2
        slit2 = r2x**2-r3x**2
        slit3 = slit1*slit1+3.0d0*slit2*slit2
        slit  = sqrt(slit3)/bigq
        if( slit .le. deldamp )then
            damper     = 0.0d0
            ddamper(:) = 0.0d0
        else
            sdamp    = rhodamp/(deldamp-slit)
            rdampsum = ald/(betadsum-rsum)
            rdamp1   = ald/(betad-r1x)
            rdamp2   = ald/(betad-r2x)
            rdamp3   = ald/(betad-r3x)
            damper   = exp(rdampsum+rdamp1+rdamp2+rdamp3+sdamp)
            dsdamp   = sdamp/(deldamp-slit)
            ddamper(1)= damper*(dsdamp*slit*r1x* (4.0d0*slit1/slit3 - 2.0d0/bigq) +rdampsum/(betadsum-rsum)+ rdamp1/ &
                & (betad- r1x))
            ddamper(2)=damper*(-dsdamp*slit*r2x* (4.0d0*(r1x*r1x-2.0d0*r2x*r2x+r3x*r3x)/slit3+2.0d0/bigq) + rdampsum/ &
                & (betadsum-rsum)+rdamp2/(betad-r2x))
            ddamper(3)=damper*(-dsdamp*slit*r3x* (4.0d0*(r1x*r1x+r2x*r2x-2.0d0*r3x*r3x)/slit3+2.0d0/bigq) + rdampsum/ &
                & (betadsum-rsum)+rdamp3/(betad-r3x))
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
        dbb3(1)  = (2.0d0*r1x-r2x-r3x)/bb3
        dbb3(2)  = (2.0d0*r2x-r1x-r3x)/bb3
        dbb3(3)  = (2.0d0*r3x-r1x-r2x)/bb3
        e1body   = (2.d0 + xmh/xmd) * hatom_bodc
        v        = e1body + ((xmd+xmh)/(2.d0*xmd))*(bodc2body(r1x, dv1) +bodc2body(r2x, dv2)) + bodc2body(r3x, dv3)
        dv(1)    = ((xmd+xmh)/(2.d0*xmd))*dv1
        dv(2)    = ((xmd+xmh)/(2.d0*xmd))*dv2
        dv(3)    = dv3
        sum1     = 0.d0
        sum2     = 0.d0
        sum3     = 0.d0
        dsum1(:) = 0.0d0
        dsum2(:) = 0.0d0
        dsum3(:) = 0.0d0
        rho1     = exp(-beta*(r1x-rnought))
        drho1    = (1.0-beta*r1x)*rho1
        rho1     = r1x*rho1
        rho2     = exp(-beta*(r2x-rnought))
        drho2    = (1.0-beta*r2x)*rho2
        rho2     = r2x*rho2
        rho3     = exp(-beta*(r3x-rnought))
        drho3    = (1.0-beta*r3x)*rho3
        rho3     = r3x*rho3
        rho1t(0) = 1.d0
        rho2t(0) = 1.d0
        rho3t(0) = 1.d0
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
    implicit none
    double precision :: r1x
    double precision :: r2x
    double precision :: r3x
    double precision :: eee
    double precision :: dv(1:3)
    call dispatch(r1x=r1x, r2x=r2x, r3x=r3x, eee=eee, dv=dv)
end subroutine potbodcdh2

!> @brief Evaluate bodc2body.
double precision function bodc2body( r1x, dv )
    implicit none
    double precision :: xlp(14), beta, re, v, r1x, dv, t
    integer          :: j, indx
    !
    !     distances are in Bohr, energies are in cm-1
    !     The global minimum is at: R=2.3178250879779552  E= -17.5988842162493135 cm-1
    !
    beta    = 3.4059169d0
    xlp( 1) = -4.942622126d0
    xlp( 2) = -4.742197464d1
    xlp( 3) = -1.115075510d2
    xlp( 4) = -1.419033104d2
    xlp( 5) = -1.165969889d2
    xlp( 6) = -6.279822284d1
    xlp( 7) = -2.626186059d1
    xlp( 8) = -1.795750569d1
    xlp( 9) = -6.826792412d0
    xlp(10) = 7.563692060d0
    xlp(11) = 2.630227612d0
    xlp(12) = -2.595885106d0
    xlp(13) = 5.647782847d-1
    xlp(14) = -4.042171086d-2
    re      = 1.401d0
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
    implicit none
    integer, intent(in)    :: mbig, mtop
    integer, intent(inout) :: iind(5000), jind(5000), kind(5000), ipar(5000), maxi, maxi2
    integer                :: l, l2, i, j, k, isum
    l  = 0
    l2 = 0
    do 10 i=0, min(mtop, mbig)
        do 20 j=0, min(mtop, mbig)
            do 30 k=j, min(mtop, mbig)
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
            30 continue
        20 continue
    10 continue
    maxi  = l2
    maxi2 = l
end

!> @brief Evaluate leps.
subroutine leps ( r, v1, v2, xma, xmb, xmc, gtot, dgtot )
    !
    !   r1 is the distance between A and B
    !   r2 is the distance between B and C
    !   r3 is the distance between A and C
    !   v1 is the lower energy root
    !   v2 is the upper energy root
    !   xma, xmb, xmc are the atomic masses
    !   ga, gb, gc are the contributions to the nonadiabatic corrections
    !      from each atom
    !
    implicit none
    integer          :: i, j
    double precision :: r(3), v1, v2, xma, xmb, xmc, gtot, dgtot(3)
    double precision :: ex, vs(3), vt(3), dvs(3), dvt(3)
    double precision :: h11, h12, h22, hsum, hdif, dhdif(3), rad, delh
    double precision :: f12(3), cos1, cos2, cos3, ga, gb, gc
    double precision :: d2vs(3), d2vt(3), dh11(3), dh12(3), dh22(3)
    double precision :: d2h12(3), d2hdif(3), drad(3), df12(3,3)
    double precision :: dcos1(3), dcos2(3), dcos3(3)
    double precision :: d(3), re(3), beta(3), sato(3), rt3
    d    = 109.4583d0/627.509552d0
    re   = 0.7412871d0/0.52917725d0
    beta = 1.973536d0*0.52917725d0
    sato = 0.132d0
    rt3  = sqrt(3.d0)
    !
    ! Compute single and triplet diatomic potentials
    do i = 1, 3
        ex      = exp(-beta(i)*(r(i)-re(i)))
        vs(i)   = d(i)*(ex - 2.0d0)*ex
        dvs(i)  = -2.0d0*d(i)*beta(i)*(ex-1.0d0)*ex
        d2vs(i) = 2.0d0*d(i)*beta(i)*beta(i)*(2.d0*ex-1.0d0)*ex
        vt(i)   = 0.5d0*d(i)*(ex + 2.0)*ex*(1.0-sato(i))/(1.0+sato(i))
        dvt(i)  = -d(i)*beta(i)*(ex+1.0)*ex*(1.0-sato(i))/(1.0+sato(i))
        d2vt(i) = d(i)*beta(i)*beta(i)*(2.0d0*ex+1.0)*ex *(1.0-sato(i))/(1.0+sato(i))
    enddo
    !
    ! Construct 2x2 DIM Hamiltonian
    h11 = d(2) + vs(2) + 0.25d0*(vs(1)+vs(3)) + 0.75d0*(vt(1)+vt(3))
    h12 = 0.25d0*rt3*(vs(1) - vs(3) - vt(1) + vt(3))
    h22 = d(2) + vt(2) + 0.75*(vs(1)+vs(3)) + 0.25d0*(vt(1)+vt(3))
    !
    ! Calculate eigenvalues (V1 and V2) and quantities needed for matrix elements
    !   of d/dRi
    !
    hsum     = h11 + h22
    hdif     = h22 - h11
    dhdif(1) = 0.5d0*(dvs(1)-dvt(1))
    dhdif(2) = -(dvs(2) - dvt(2))
    dhdif(3) = 0.5d0*(dvs(3)-dvt(3))
    rad      = 0.25d0*hdif*hdif + h12*h12
    delh     = sqrt(rad)
    v1       = 0.5d0*hsum - delh
    v2       = 0.5d0*hsum + delh
    ! derivatives
    !
    dh11(1)   = 0.25d0*dvs(1)+0.75d0*dvt(1)
    dh11(2)   = dvs(2)
    dh11(3)   = 0.25d0*dvs(3)+0.75d0*dvt(3)
    dh12(1)   = 0.25d0*rt3*(dvs(1)-dvt(1))
    dh12(2)   = 0.0d0
    dh12(3)   = -0.25d0*rt3*(dvs(3)-dvt(3))
    dh22(1)   = 0.75d0*dvs(1)+0.25d0*dvt(1)
    dh22(2)   = dvt(2)
    dh22(3)   = 0.75d0*dvs(3)+0.25d0*dvt(3)
    d2h12(1)  = 0.25d0*rt3*(d2vs(1)-d2vt(1))
    d2h12(2)  = 0.0d0
    d2h12(3)  = -0.25d0*rt3*(d2vs(3)-d2vt(3))
    d2hdif(1) = 0.5d0*(d2vs(1) - d2vt(1))
    d2hdif(2) = -(d2vs(2) - d2vt(2))
    d2hdif(3) = 0.5d0*(d2vs(3) - d2vt(3))
    drad(1)   = 0.5d0*hdif*dhdif(1) +2.0d0*h12*dh12(1)
    drad(2)   = 0.5d0*hdif*dhdif(2) +2.0d0*h12*dh12(2)
    drad(3)   = 0.5d0*hdif*dhdif(3) +2.0d0*h12*dh12(3)
    ! Calculate matrix elements of d/dRi; note that at conical intersection f12
    !   diverges, but is set to zero in the code.
    do i = 1, 3
        f12(i) = 0.25d0*(h12*dhdif(i) - hdif*dh12(i))
        if( rad .ne. 0.0d0 )then
            f12(i) = f12(i)/rad
            do j = 1, 3
                df12(i,j) = (0.25d0*(dh12(j)*dhdif(i)-dhdif(j)*dh12(i)) -f12(i)*drad(j))/rad
            enddo
            df12(i,i) = df12(i,i) + 0.25d0*(h12*d2hdif(i)-hdif*d2h12(i))/rad
        endif
    enddo
    !
    ! Calculate diagonal matrix element of kinetic energy operator
    !
    cos1     = (r(2)*r(2) + r(3)*r(3) - r(1)*r(1))/(2.0*r(2)*r(3))
    cos2     = (r(1)*r(1) + r(3)*r(3) - r(2)*r(2))/(2.0*r(1)*r(3))
    cos3     = (r(1)*r(1) + r(2)*r(2) - r(3)*r(3))/(2.0*r(1)*r(2))
    ga       = (f12(1)*(f12(1) + 2.0*cos2*f12(3)) + f12(3)*f12(3))/(2.d0*xma)
    gb       = (f12(1)*(f12(1) + 2.0*cos3*f12(2)) + f12(2)*f12(2))/(2.d0*xmb)
    gc       = (f12(2)*(f12(2) + 2.0*cos1*f12(3)) + f12(3)*f12(3))/(2.d0*xmc)
    gtot     = ga+gb+gc
    dcos1(1) = -r(1)/(r(2)*r(3))
    dcos1(2) = 1.0d0/r(3) - cos1/r(2)
    dcos1(3) = 1.0d0/r(2) - cos1/r(3)
    dcos2(1) = 1.0d0/r(3) - cos2/r(1)
    dcos2(2) = -r(2)/(r(1)*r(3))
    dcos2(3) = 1.0d0/r(1) - cos2/r(3)
    dcos3(1) = 1.0d0/r(2) - cos3/r(1)
    dcos3(2) = 1.0d0/r(1) - cos3/r(2)
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
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r1x, r2x, r3x, eee, dv )
        implicit none
        optional :: r1x, r2x, r3x, eee, dv
        integer            :: iind, jind, kind, ipar, maxi, maxi2
        integer            :: maxi2a, maxi2b, ii, i, it
        integer, parameter :: n2 = 5000
        dimension iind(n2), jind(n2), kind(n2), ipar(n2)
        double precision            :: epslon, epscem, ald, betad
        double precision            :: r1x, r2x, r3x, eee, dv(:)
        double precision            :: etrip1, etrip2, etrip3, detrip1, detrip2, detrip3
        double precision            :: es1, es2, es3, des1, des2, des3
        double precision            :: q1, q2, q3, xj1, xj2, xj3, xj, vlondon
        double precision            :: dq1, dq2, dq3, dxj1, dxj2, dxj3, dvlon(3)
        double precision            :: damper, ddamp(3), rdamp1, rdamp2, rdamp3
        double precision            :: bb, dbb(3), warp, dwarp(3), rho1, rho2, rho3
        double precision            :: rho1t(0:12), rho2t(0:12), rho3t(0:12)
        double precision            :: sum1, sum2, sum3, sum4
        double precision            :: dsum1(3), dsum2(3), dsum3(3), dsum4(3)
        double precision            :: rprod, t1, v3c
        double precision, parameter :: eshift  = 0.1744759989649372d0
        double precision, parameter :: beta    = 0.9246250191375499d0
        double precision, parameter :: rnought = 2.d0
        double precision            :: xlin1(89), xlin2(89), xlin3(71), xlin4(71)
        ald       = 0.d0
        bb        = 0.d0
        betad     = 0.d0
        damper    = 0.d0
        dbb       = 0.d0
        ddamp     = 0.d0
        des1      = 0.d0
        des2      = 0.d0
        des3      = 0.d0
        detrip1   = 0.d0
        detrip2   = 0.d0
        detrip3   = 0.d0
        dq1       = 0.d0
        dq2       = 0.d0
        dq3       = 0.d0
        dsum1     = 0.d0
        dsum2     = 0.d0
        dsum3     = 0.d0
        dsum4     = 0.d0
        dvlon     = 0.d0
        dwarp     = 0.d0
        dxj1      = 0.d0
        dxj2      = 0.d0
        dxj3      = 0.d0
        epscem    = 0.d0
        epslon    = 0.d0
        es1       = 0.d0
        es2       = 0.d0
        es3       = 0.d0
        etrip1    = 0.d0
        etrip2    = 0.d0
        etrip3    = 0.d0
        q1        = 0.d0
        q2        = 0.d0
        q3        = 0.d0
        rdamp1    = 0.d0
        rdamp2    = 0.d0
        rdamp3    = 0.d0
        rho1      = 0.d0
        rho1t     = 0.d0
        rho2      = 0.d0
        rho2t     = 0.d0
        rho3      = 0.d0
        rho3t     = 0.d0
        rprod     = 0.d0
        sum1      = 0.d0
        sum2      = 0.d0
        sum3      = 0.d0
        sum4      = 0.d0
        t1        = 0.d0
        v3c       = 0.d0
        vlondon   = 0.d0
        warp      = 0.d0
        xj        = 0.d0
        xj1       = 0.d0
        xj2       = 0.d0
        xj3       = 0.d0
        xlin1     = 0.d0
        xlin2     = 0.d0
        xlin3     = 0.d0
        xlin4     = 0.d0
        xlin1(1)  = -0.1682626669809615d0
        xlin1(2)  = 1.697952225005269d0
        xlin1(3)  = 3.959367441888326d0
        xlin1(4)  = -7.159854944676235d0
        xlin1(5)  = -2.462488887191634d0
        xlin1(6)  = 7.998434264390882d0
        xlin1(7)  = -3.665712084901418d0
        xlin1(8)  = -0.2155463656483861d0
        xlin1(9)  = 0.5620418518154215d0
        xlin1(10) = -0.1369812453279634d0
        xlin1(11) = 1.0063614579308723d-2
        xlin1(12) = 41.07358514471853d0
        xlin1(13) = 67.02217912830267d0
        xlin1(14) = -100.0014758720550d0
        xlin1(15) = 2.481875721304641d0
        xlin1(16) = 32.46244740624601d0
        xlin1(17) = -15.32460432587821d0
        xlin1(18) = 3.814631583515083d0
        xlin1(19) = -0.3155985331636701d0
        xlin1(20) = -4.2198516063816239d-3
        xlin1(21) = -157.7203098415152d0
        xlin1(22) = 17.11141015275416d0
        xlin1(23) = -61.48638359701442d0
        xlin1(24) = 4.244841908917134d0
        xlin1(25) = 5.300393723900312d0
        xlin1(26) = -0.6170558358973661d0
        xlin1(27) = -6.2143674437616274d-3
        xlin1(28) = 5.785766830955630d0
        xlin1(29) = -0.1133784523542150d0
        xlin1(30) = -0.2539145873400511d0
        xlin1(31) = -4.3406559735048131d-2
        xlin1(32) = 2.4548851428702549d-2
        xlin1(33) = -1.092841685240530d0
        xlin1(34) = 0.8131178021738225d0
        xlin1(35) = -7.7933784638615652d-2
        xlin1(36) = -0.1484455564689179d0
        xlin1(37) = -7.351078313940695d0
        xlin1(38) = 128.7480567271992d0
        xlin1(39) = 199.6058760392980d0
        xlin1(40) = -50.18444648711906d0
        xlin1(41) = -178.2184990028168d0
        xlin1(42) = 129.6832433574917d0
        xlin1(43) = -37.57036725437851d0
        xlin1(44) = 5.747838042162251d0
        xlin1(45) = 0.1434298979682676d0
        xlin1(46) = -8.3110343245055174d-2
        xlin1(47) = 428.7635440340811d0
        xlin1(48) = 495.9564113906064d0
        xlin1(49) = -1742.123837168904d0
        xlin1(50) = -259.0564024322437d0
        xlin1(51) = 250.2352551214438d0
        xlin1(52) = -44.15232947384926d0
        xlin1(53) = 0.9963955068997190d0
        xlin1(54) = 0.3136096547683901d0
        xlin1(55) = -1032.921631265462d0
        xlin1(56) = 199.9708129167335d0
        xlin1(57) = -135.4683186684037d0
        xlin1(58) = -14.50534672808297d0
        xlin1(59) = 5.650612567917722d0
        xlin1(60) = -0.1473535040913008d0
        xlin1(61) = 1.351710421967908d0
        xlin1(62) = 29.60972395839821d0
        xlin1(63) = -7.224575522867537d0
        xlin1(64) = 0.3582134126835411d0
        xlin1(65) = 0.2589923362572861d0
        xlin1(66) = 0.1819402887127863d0
        xlin1(67) = 1312.633253770123d0
        xlin1(68) = 5189.281073844522d0
        xlin1(69) = -2727.020597875299d0
        xlin1(70) = -122.2193440638810d0
        xlin1(71) = 152.4859537357477d0
        xlin1(72) = -5.562127541501408d0
        xlin1(73) = -1.549418235204051d0
        xlin1(74) = 2647.710764337011d0
        xlin1(75) = 502.1991025406292d0
        xlin1(76) = -152.3069488960497d0
        xlin1(77) = -5.740041181528469d0
        xlin1(78) = 1.362209926084833d0
        xlin1(79) = 17.96504193139269d0
        xlin1(80) = 10.13930464010770d0
        xlin1(81) = -2.287463377365149d0
        xlin1(82) = 1.932492687974671d0
        xlin1(83) = -3337.181754475553d0
        xlin1(84) = 182.9671187483376d0
        xlin1(85) = 20.29606065582956d0
        xlin1(86) = 1.192500349837268d0
        xlin1(87) = -22.69079850990068d0
        xlin1(88) = -0.5967747962499963d0
        xlin1(89) = 1.091111370927579d0
        xlin2(1)  = 3.2289532152874986d-2
        xlin2(2)  = -8.6686506736620789d-2
        xlin2(3)  = -0.3266024449171150d0
        xlin2(4)  = -0.7972798579742405d0
        xlin2(5)  = 2.812854025206053d0
        xlin2(6)  = -2.374871328104482d0
        xlin2(7)  = 0.6089164699935086d0
        xlin2(8)  = 0.1635658322754857d0
        xlin2(9)  = -0.1263080841790984d0
        xlin2(10) = 2.5221027611065153d-2
        xlin2(11) = -1.6677760633386995d-3
        xlin2(12) = -0.3372069966080046d0
        xlin2(13) = -8.377537389356892d0
        xlin2(14) = 3.694668699246826d0
        xlin2(15) = 11.36868357221171d0
        xlin2(16) = -11.88889428297186d0
        xlin2(17) = 4.270226804928180d0
        xlin2(18) = -0.4526863369836274d0
        xlin2(19) = -5.2518975083660345d-2
        xlin2(20) = 9.0223742463199624d-3
        xlin2(21) = 38.03127084604298d0
        xlin2(22) = -12.82065963469141d0
        xlin2(23) = -5.663705994784566d0
        xlin2(24) = 7.300080908068341d0
        xlin2(25) = -0.5157232521089076d0
        xlin2(26) = -0.2200954472805637d0
        xlin2(27) = 2.0686337667851093d-2
        xlin2(28) = 14.19614073502447d0
        xlin2(29) = 1.726255538125508d0
        xlin2(30) = -2.526391344985463d0
        xlin2(31) = 6.0598667264186397d-2
        xlin2(32) = 3.8560006141770319d-2
        xlin2(33) = 2.306725001955670d0
        xlin2(34) = -0.3197660708712612d0
        xlin2(35) = -1.9982819587147187d-2
        xlin2(36) = 0.1036787525350017d0
        xlin2(37) = -2.030231243513758d0
        xlin2(38) = 25.16474735933574d0
        xlin2(39) = -14.13898387128746d0
        xlin2(40) = -62.01219930765056d0
        xlin2(41) = 75.56471352787402d0
        xlin2(42) = -30.25154027581496d0
        xlin2(43) = 4.257663341071418d0
        xlin2(44) = 0.6142213624561507d0
        xlin2(45) = -0.2458345005137628d0
        xlin2(46) = 1.9135101060142316d-2
        xlin2(47) = 168.2482649287792d0
        xlin2(48) = -270.3587962365308d0
        xlin2(49) = 72.46467150857366d0
        xlin2(50) = 32.12278113625301d0
        xlin2(51) = 11.89435258275814d0
        xlin2(52) = -5.521268211450904d0
        xlin2(53) = 0.8934653854937415d0
        xlin2(54) = -7.6805320500971247d-2
        xlin2(55) = -422.2672942154591d0
        xlin2(56) = -35.92489169679082d0
        xlin2(57) = 14.74665733902546d0
        xlin2(58) = 6.020706922914545d0
        xlin2(59) = 2.4538599994799737d-3
        xlin2(60) = -0.1370624616096560d0
        xlin2(61) = 44.47654636843973d0
        xlin2(62) = -11.18166310447657d0
        xlin2(63) = 3.6110526441458540d-2
        xlin2(64) = 0.1539535582849851d0
        xlin2(65) = 2.065981296572923d0
        xlin2(66) = -0.1624439939014451d0
        xlin2(67) = -658.3723053145102d0
        xlin2(68) = -883.4947620258348d0
        xlin2(69) = -263.8479489759550d0
        xlin2(70) = 44.68359598895833d0
        xlin2(71) = 29.88970883430300d0
        xlin2(72) = -5.663152586850906d0
        xlin2(73) = 0.2413124485342995d0
        xlin2(74) = 656.9813540494517d0
        xlin2(75) = -91.56710742448996d0
        xlin2(76) = -30.38175110815412d0
        xlin2(77) = 2.066013906655882d0
        xlin2(78) = 0.5314632792073760d0
        xlin2(79) = 54.99422880821747d0
        xlin2(80) = -0.7887249891070589d0
        xlin2(81) = -0.9164068632403487d0
        xlin2(82) = 0.6576346543435602d0
        xlin2(83) = 292.2665706493461d0
        xlin2(84) = -32.04152511166628d0
        xlin2(85) = -2.104943871500433d0
        xlin2(86) = -0.7251129002870176d0
        xlin2(87) = 4.709122797258667d0
        xlin2(88) = 1.119971732384046d0
        xlin2(89) = -4.066519840824071d0
        xlin3(1)  = 6.3998636332917941d-2
        xlin3(2)  = -0.5826083150212633d0
        xlin3(3)  = -3.294968443823759d0
        xlin3(4)  = 5.751817019955082d0
        xlin3(5)  = -1.486176533347650d0
        xlin3(6)  = -0.6845411439568367d0
        xlin3(7)  = -0.9083627566591157d0
        xlin3(8)  = 1.305678536333854d0
        xlin3(9)  = -0.4821780067143284d0
        xlin3(10) = 5.5182580699196751d-2
        xlin3(11) = -22.57213192201081d0
        xlin3(12) = -61.37436407473290d0
        xlin3(13) = 52.64511100643983d0
        xlin3(14) = 29.99521825023724d0
        xlin3(15) = -38.31341936937510d0
        xlin3(16) = 14.86681566870566d0
        xlin3(17) = -3.460352146284503d0
        xlin3(18) = 0.1343449644795501d0
        xlin3(19) = 53.21459427503272d0
        xlin3(20) = 64.60699757442232d0
        xlin3(21) = 30.89466488573955d0
        xlin3(22) = 8.412732049211240d0
        xlin3(23) = -4.274560634850636d0
        xlin3(24) = -0.1788031871862021d0
        xlin3(25) = -46.47154313325355d0
        xlin3(26) = 24.88995914944731d0
        xlin3(27) = -2.757433396373881d0
        xlin3(28) = -0.3635530438220282d0
        xlin3(29) = 1.049428472748707d0
        xlin3(30) = -0.9909460260393127d0
        xlin3(31) = 7.623381850055868d0
        xlin3(32) = -75.76387696233419d0
        xlin3(33) = -162.8879929453359d0
        xlin3(34) = -29.98317728892977d0
        xlin3(35) = 168.7972614497797d0
        xlin3(36) = -102.0864921135087d0
        xlin3(37) = 31.38363768116054d0
        xlin3(38) = -6.084136544334489d0
        xlin3(39) = -2.9400115678577254d-2
        xlin3(40) = -417.4808742201730d0
        xlin3(41) = -601.3246381855070d0
        xlin3(42) = 1111.017825017578d0
        xlin3(43) = 749.0346863822200d0
        xlin3(44) = -254.9798923160399d0
        xlin3(45) = 37.02896288679148d0
        xlin3(46) = -0.7010161055223685d0
        xlin3(47) = 928.6100310886904d0
        xlin3(48) = 513.3554862860957d0
        xlin3(49) = 138.4248255502364d0
        xlin3(50) = -0.9347938923678116d0
        xlin3(51) = 3.867048618649862d0
        xlin3(52) = -109.6794304656507d0
        xlin3(53) = 10.30110710850976d0
        xlin3(54) = 0.5118243051520034d0
        xlin3(55) = 0.2069543381912151d0
        xlin3(56) = -796.4010440692600d0
        xlin3(57) = -4693.456970734772d0
        xlin3(58) = 2597.246850944603d0
        xlin3(59) = 630.4845628327919d0
        xlin3(60) = -221.0729177223342d0
        xlin3(61) = 9.839396652075672d0
        xlin3(62) = -5227.790162633692d0
        xlin3(63) = -529.4879950032479d0
        xlin3(64) = 30.53701394166368d0
        xlin3(65) = -14.22864486129778d0
        xlin3(66) = -35.14731691232372d0
        xlin3(67) = 2.140914121987492d0
        xlin3(68) = 2250.367423678786d0
        xlin3(69) = 138.5994085512685d0
        xlin3(70) = 70.30416282155014d0
        xlin3(71) = -62.74844809925895d0
        xlin4(1)  = -1.2612782515309912d-2
        xlin4(2)  = 1.5724446750473936d-2
        xlin4(3)  = 0.2653887372495336d0
        xlin4(4)  = 0.1216019078797281d0
        xlin4(5)  = -0.6102844930670878d0
        xlin4(6)  = 4.6420446436398755d-3
        xlin4(7)  = 0.5563651571664677d0
        xlin4(8)  = -0.3890198684837317d0
        xlin4(9)  = 0.1052569138798340d0
        xlin4(10) = -1.0080347996402828d-2
        xlin4(11) = 0.2409741384095758d0
        xlin4(12) = 5.745106082351191d0
        xlin4(13) = 1.047295525278513d0
        xlin4(14) = -12.03507114814484d0
        xlin4(15) = 9.035186991905373d0
        xlin4(16) = -2.267198360038944d0
        xlin4(17) = 2.0836099694222385d-2
        xlin4(18) = 4.7445938282087222d-2
        xlin4(19) = -18.26537244396756d0
        xlin4(20) = -4.809068965727743d0
        xlin4(21) = 5.658647992475390d0
        xlin4(22) = -1.284093508478602d0
        xlin4(23) = -2.164733987255563d0
        xlin4(24) = 0.3630073042947094d0
        xlin4(25) = 3.760490137844927d0
        xlin4(26) = -5.927286374087972d0
        xlin4(27) = 0.1278839922444694d0
        xlin4(28) = 0.5381095462111573d0
        xlin4(29) = -2.179370013756247d0
        xlin4(30) = 0.4127566904703107d0
        xlin4(31) = 0.7168269469196163d0
        xlin4(32) = -9.502943939736065d0
        xlin4(33) = -9.597638311349975d0
        xlin4(34) = 65.60860814240448d0
        xlin4(35) = -56.04003254284753d0
        xlin4(36) = 15.38423202570607d0
        xlin4(37) = -0.3167914651468724d0
        xlin4(38) = -0.7643266917502949d0
        xlin4(39) = 0.1110705439245336d0
        xlin4(40) = -146.5638710976931d0
        xlin4(41) = 151.8079548350792d0
        xlin4(42) = 15.86486188082157d0
        xlin4(43) = -65.56160562330304d0
        xlin4(44) = -1.126379669058077d0
        xlin4(45) = -1.759058937005427d0
        xlin4(46) = -3.1582037336166588d-4
        xlin4(47) = 410.9371301349821d0
        xlin4(48) = 142.2842073695382d0
        xlin4(49) = -15.14809406881419d0
        xlin4(50) = -9.858662824266736d0
        xlin4(51) = 0.2622049408092991d0
        xlin4(52) = -72.67984057864001d0
        xlin4(53) = 11.28350883179339d0
        xlin4(54) = -0.5323962789000201d0
        xlin4(55) = -2.447540583751350d0
        xlin4(56) = 148.9575982135905d0
        xlin4(57) = 1249.939675860206d0
        xlin4(58) = 268.6478693121586d0
        xlin4(59) = 33.67689052007242d0
        xlin4(60) = -35.65384689269180d0
        xlin4(61) = 3.464021415311834d0
        xlin4(62) = -103.4190791219546d0
        xlin4(63) = 76.76961997409407d0
        xlin4(64) = 11.59182442022548d0
        xlin4(65) = -0.5254501536114633d0
        xlin4(66) = -32.07157614897808d0
        xlin4(67) = 1.102243506822152d0
        xlin4(68) = -952.9757997893194d0
        xlin4(69) = 52.11461967121677d0
        xlin4(70) = 6.929202177070909d0
        xlin4(71) = -6.322557284152793d0
        epslon    = 1.d-12
        epscem    = 1.d-12
        ald       = 0.02d0
        betad     = 0.72d0
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
        q1       = 0.5d0*(es1+etrip1)
        q2       = 0.5d0*(es2+etrip2)
        q3       = 0.5d0*(es3+etrip3)
        xj1      = 0.5d0*(es1-etrip1)
        xj2      = 0.5d0*(es2-etrip2)
        xj3      = 0.5d0*(es3-etrip3)
        xj       = sqrt(epslon+0.5d0*((xj3-xj1)**2+(xj2-xj1)**2+(xj3-xj2)**2))
        vlondon  = q1+q2+q3-xj
        dq1      = 0.5d0*(des1+detrip1)
        dq2      = 0.5d0*(des2+detrip2)
        dq3      = 0.5d0*(des3+detrip3)
        dxj1     = 0.5d0*(des1-detrip1)
        dxj2     = 0.5d0*(des2-detrip2)
        dxj3     = 0.5d0*(des3-detrip3)
        dvlon(1) = dq1-0.5d0*dxj1*(2.0d0*xj1-xj2-xj3)/xj
        dvlon(2) = dq2-0.5d0*dxj2*(2.0d0*xj2-xj1-xj3)/xj
        dvlon(3) = dq3-0.5d0*dxj3*(2.0d0*xj3-xj1-xj2)/xj
        if( (r1x .le. betad).or.(r2x .le. betad).or.(r3x .le. betad) )then
            damper   = 0.d0
            ddamp(1) = 0.d0
            ddamp(2) = 0.d0
            ddamp(3) = 0.d0
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
        dbb(1)   = (2.d0*r1x-r2x-r3x)/bb
        dbb(2)   = (2.d0*r2x-r1x-r3x)/bb
        dbb(3)   = (2.d0*r3x-r1x-r2x)/bb
        warp     = 1.d0/(1.d0/r1x+1.d0/r2x+1.d0/r3x)
        dwarp(1) = (warp/r1x)**2
        dwarp(2) = (warp/r2x)**2
        dwarp(3) = (warp/r3x)**2
        rho1     = exp(-beta*(r1x-rnought))
        rho2     = exp(-beta*(r2x-rnought))
        rho3     = exp(-beta*(r3x-rnought))
        rho1t(0) = 1.d0
        rho2t(0) = 1.d0
        rho3t(0) = 1.d0
        do it=1, 12
            rho1t(it) = rho1t(it-1)*rho1
            rho2t(it) = rho2t(it-1)*rho2
            rho3t(it) = rho3t(it-1)*rho3
        enddo
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0
        sum4 = 0.d0
        do i = 1, 3
            dsum1(i) = 0.d0
            dsum2(i) = 0.d0
            dsum3(i) = 0.d0
            dsum4(i) = 0.d0
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
    implicit none
    double precision :: r1x
    double precision :: r2x
    double precision :: r3x
    double precision :: eee
    double precision :: dv(3)
    call dispatch(r1x=r1x, r2x=r2x, r3x=r3x, eee=eee, dv=dv)
end subroutine potcbs

!> @brief Evaluate indexa3.
subroutine indexa3( mbig, mtop, iind, jind, kind, ipar, maxi, maxi2 )
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
                    else if( j .ne. k )then
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
    implicit none
    integer                     :: j
    double precision            :: r, v, dv, t0, t1, t2, vlr, dvlr, vsr, dvsr
    double precision            :: damp, xd, xd2, prefac
    double precision, parameter :: c6    = -6.499027d0
    double precision, parameter :: c8    = -124.3991d0
    double precision, parameter :: c10   = -3285.828d0
    double precision, parameter :: beta  = 0.2d0
    double precision, parameter :: re    = 1.401d0
    double precision, parameter :: alpha = 4.342902497495857d0
    double precision            :: xlp(17)
    xlp(1)  = 2.190254292077946d-1
    xlp(1)  = 2.190254292077946d-1
    xlp(1)  = 2.190254292077946d-1
    xlp(2)  = 6.903970886899540d-1
    xlp(3)  = 1.163805164813647d+0
    xlp(4)  = 1.337274224278977d+0
    xlp(5)  = 1.176131501286896d+0
    xlp(6)  = 9.093025316067311d-1
    xlp(7)  = 5.908490740317524d-1
    xlp(8)  = 1.467271248985162d-1
    xlp(9)  = 6.218667627370255d-2
    xlp(10) = 1.105747186790804d-1
    xlp(11) = -7.112594221275334d-2
    xlp(12) = -8.866929977503379d-3
    xlp(13) = 3.823523348718375d-2
    xlp(14) = -2.055915301253834d-2
    xlp(15) = 5.419433399532862d-3
    xlp(16) = -7.263483303294653d-4
    xlp(17) = 4.154351972583585d-5
    t0      = r-re
    prefac  = exp(-alpha*t0)
    if( (r-re) .ne. 0.d0 )then
        vsr  = 0.0
        dvsr = 0.0
        t1   = 1.0d0
        t2   = 0.0d0
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
    damp = 1.d0-t2
    xd   = damp/r
    xd2  = xd*xd
    vlr  = (c6+xd2*(c8+xd2*c10))*xd2**3
    dvlr = -(6.0d0*c6+xd2*(8.0d0*c8+xd2*10.0d0*c10))* ((1/r)-2.0d0*beta*t2/xd)*xd2**3
    v    = vlr+vsr
    dv   = dvlr+dvsr
end

!> @brief Evaluate singlet.
subroutine singlet( r, v, dv )
    implicit none
    integer                     :: j
    double precision            :: r, v, dv, t0, t1, t2, vsr, damp, xd, xd2, prefac
    double precision, parameter :: c6    = -6.499027d0
    double precision, parameter :: c8    = -124.3991d0
    double precision, parameter :: c10   = -3285.828d0
    double precision, parameter :: beta  = 0.2d0
    double precision, parameter :: re    = 1.401d0
    double precision, parameter :: alpha = 3.980917850296971d0
    double precision            :: xlp(17)
    xlp(1)  = -1.709663318869429d-1
    xlp(2)  = -6.675286769482015d-1
    xlp(3)  = -1.101876055072129d+0
    xlp(4)  = -1.106460095658282d+0
    xlp(5)  = -7.414724284525231d-1
    xlp(6)  = -3.487882731923523d-1
    xlp(7)  = -1.276736255756598d-1
    xlp(8)  = -5.875965709151867d-2
    xlp(9)  = -4.030128017933840d-2
    xlp(10) = -2.038653237221007d-2
    xlp(11) = -7.639198558462706d-4
    xlp(12) = 2.912954920483885d-3
    xlp(13) = -2.628273116815280d-4
    xlp(14) = -4.622088855684211d-4
    xlp(15) = 1.278468948126147d-4
    xlp(16) = -1.157434070240206d-5
    xlp(17) = -2.609691840882097d-12
    t1      = r*r
    t2      = exp(-beta*t1)
    damp    = 1.d0-t2
    xd      = damp/r
    xd2     = xd*xd
    v       = (c6+xd2*(c8+xd2*c10))*xd2**3
    dv      = -(6.0d0*c6+xd2*(8.0d0*c8+xd2*10.0d0*c10))* ((1/r)-2.0d0*beta*t2/xd)*xd2**3
    t0      = r-re
    prefac  = exp(-alpha*t0)
    vsr     = 0.0
    if( (r-re) .ne. 0.d0 )then
        t1 = 1.0d0
        t2 = 0.0d0
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
    implicit none
    external fcn
    double precision, parameter :: r = 0.61803399d0, c = 1d0-r
    double precision            :: f1, f2, x0, x1, x2, x3, tol
    double precision            :: fval, xmin, bx, ax, cx, fcn
    double precision            :: fa, fb, fc
    integer                     :: icount
    integer, parameter          :: maxcount = 5000
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
double precision function fcn( x )
    use zmasses, only: zmass1, zmass2
    implicit none
    double precision :: x, e, bodc2body, v, dv, towave, xmh, massfact
    !      Evaluates the Born-Huang energy of an isomer of H2 with masses of zmass1 and zmass2
    xmh      = 1836.15264d0   ! H nuclear mass (after substracting 1 for electron mass)
    massfact = 0.5d0*xmh*(zmass1+zmass2)/(zmass1*zmass2)
    towave   = 219474.7d0
    call singlet(x, v, dv)
    fcn = massfact*bodc2body(x, dv)/towave+v
end

!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_h3_gen_bh_2009_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_legacy_units,       only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: legacy_triatomic_distances_to_bohr
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: r(3), dvp(3), v
        double precision              :: tx(9), r2(3), drdx(3,9)
        double precision              :: gx(nstates,9)
        integer                       :: iatom, idir, i, j
        !initialize
        v = 0.d0; g = 0.d0; d = 0.d0
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is HHH
        call legacy_triatomic_distances_to_bohr(x, r)
        r2 = r*legacy_bohr_to_angstrom
        call pot(r, v, dvp)
        p   = v*legacy_hartree_to_electronvolt
        dvp = dvp*51.422067
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

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_bh_2009_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX.
subroutine evdrdx( x, r, drdx )
    implicit none
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    integer                       :: i, j
    ! Initialize dRdX(3,9)
    do i=1, 3
        do j=1, 9
            drdx(i,j) = 0.0d0
        enddo
    enddo
    drdx(1,1) = (x(1)-x(4))/r(1)
    drdx(1,2) = (x(2)-x(5))/r(1)
    drdx(1,3) = (x(3)-x(6))/r(1)
    drdx(1,4) = -drdx(1,1)
    drdx(1,5) = -drdx(1,2)
    drdx(1,6) = -drdx(1,3)
    drdx(2,4) = (x(4)-x(7))/r(2)
    drdx(2,5) = (x(5)-x(8))/r(2)
    drdx(2,6) = (x(6)-x(9))/r(2)
    drdx(2,7) = -drdx(2,4)
    drdx(2,8) = -drdx(2,5)
    drdx(2,9) = -drdx(2,6)
    drdx(3,1) = (x(1)-x(7))/r(3)
    drdx(3,2) = (x(2)-x(8))/r(3)
    drdx(3,3) = (x(3)-x(9))/r(3)
    drdx(3,7) = -drdx(3,1)
    drdx(3,8) = -drdx(3,2)
    drdx(3,9) = -drdx(3,3)
end subroutine
