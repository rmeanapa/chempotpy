!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!   System:                        N4
!   Functional form:               permutation-invariant polynomials
!   Common name:                   N4(adiabatic ground state)
!   Number of derivatives:         1
!   Number of bodies:              4
!   Number of electronic surfaces: 1
!   Interface:                     Section-2
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
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_n4_1a_pip_2013_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, p, g, d )
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 4
        double precision, parameter   :: eref    = -218.40801323d0
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: coord_r(6)
        double precision              :: epot_b(276)
        double precision              :: epot_c(276)
        double precision              :: epot_p(0:305)
        double precision              :: epot_rm(0:111)
        double precision              :: epot_rms(6)
        double precision              :: gradt_dbdr(6,276)
        double precision              :: gradt_dmdr(6,0:111)
        double precision              :: gradt_dmsdr(6,6)
        double precision              :: gradt_dpdr(6,0:305)
        double precision              :: gradt_drdx(6,12)
        double precision              :: gradt_dvdr(6)
        double precision              :: msprmt_a
        double precision              :: msprmt_re
        double precision              :: v, tx(12), tg(12)
        integer                       :: iatom, idir, j, istate
        !initialize
        coord_r     = 0.d0; epot_b = 0.d0; epot_c = 0.d0; epot_p = 0.d0; epot_rm = 0.d0; epot_rms = 0.d0; gradt_dbdr = 0.d0
        gradt_dmdr  = 0.d0; gradt_dmsdr = 0.d0; gradt_dpdr = 0.d0; gradt_drdx = 0.d0; gradt_dvdr = 0.d0; msprmt_a = 0.d0
        msprmt_re   = 0.d0; msprmt_a = 0.9d0; msprmt_re = 1.098d0
        epot_c(  1) = 0.102403978367d+04
        epot_c(  2) = -0.329677065461d+04
        epot_c(  3) = 0.563720430512d+04
        epot_c(  4) = -0.398579928313d+04
        epot_c(  5) = -0.385965542114d+04
        epot_c(  6) = 0.388707170536d+04
        epot_c(  7) = 0.313073995764d+04
        epot_c(  8) = 0.448635322841d+04
        epot_c(  9) = 0.145177478364d+05
        epot_c( 10) = -0.520482498458d+04
        epot_c( 11) = 0.927577190398d+04
        epot_c( 12) = 0.577514178728d+04
        epot_c( 13) = 0.549924812055d+04
        epot_c( 14) = 0.128789377708d+05
        epot_c( 15) = -0.123915284843d+05
        epot_c( 16) = -0.152516198607d+05
        epot_c( 17) = -0.311653729683d+05
        epot_c( 18) = 0.241802095549d+04
        epot_c( 19) = 0.552823062458d+04
        epot_c( 20) = -0.157606516638d+05
        epot_c( 21) = 0.183999978622d+05
        epot_c( 22) = -0.730065813635d+04
        epot_c( 23) = -0.154254898066d+04
        epot_c( 24) = -0.270137001447d+05
        epot_c( 25) = -0.712833478072d+04
        epot_c( 26) = -0.157507777456d+05
        epot_c( 27) = 0.150777397542d+03
        epot_c( 28) = -0.460118218767d+04
        epot_c( 29) = 0.621757808397d+05
        epot_c( 30) = -0.387325947409d+05
        epot_c( 31) = 0.956015975839d+05
        epot_c( 32) = 0.598731286868d+05
        epot_c( 33) = 0.184378520851d+05
        epot_c( 34) = 0.118525258622d+05
        epot_c( 35) = -0.192034792924d+05
        epot_c( 36) = 0.102203620381d+05
        epot_c( 37) = -0.358701539635d+04
        epot_c( 38) = 0.174513120739d+05
        epot_c( 39) = 0.620986763680d+05
        epot_c( 40) = 0.382171290501d+05
        epot_c( 41) = 0.372876912748d+04
        epot_c( 42) = 0.458049260553d+04
        epot_c( 43) = 0.974878688873d+05
        epot_c( 44) = 0.182211600789d+05
        epot_c( 45) = -0.117103856781d+05
        epot_c( 46) = 0.889576596375d+04
        epot_c( 47) = 0.176322778899d+05
        epot_c( 48) = -0.292083982135d+05
        epot_c( 49) = 0.394476552928d+04
        epot_c( 50) = -0.490908973766d+04
        epot_c( 51) = 0.341307888032d+05
        epot_c( 52) = -0.186725715169d+05
        epot_c( 53) = 0.147822633435d+05
        epot_c( 54) = -0.160310707626d+04
        epot_c( 55) = 0.529073124595d+03
        epot_c( 56) = 0.436589884373d+04
        epot_c( 57) = -0.275977248391d+05
        epot_c( 58) = 0.916659840404d+05
        epot_c( 59) = 0.797381187300d+04
        epot_c( 60) = -0.266127694396d+05
        epot_c( 61) = -0.366120071607d+05
        epot_c( 62) = 0.302613610927d+05
        epot_c( 63) = -0.472403494286d+05
        epot_c( 64) = -0.437423663348d+04
        epot_c( 65) = -0.299548954583d+05
        epot_c( 66) = -0.650332702893d+05
        epot_c( 67) = 0.710278406294d+05
        epot_c( 68) = -0.754784583336d+05
        epot_c( 69) = 0.230385041533d+05
        epot_c( 70) = -0.371044004395d+05
        epot_c( 71) = -0.106153162568d+04
        epot_c( 72) = 0.415993233011d+05
        epot_c( 73) = 0.180415447082d+05
        epot_c( 74) = -0.181546212206d+05
        epot_c( 75) = 0.133142608504d+04
        epot_c( 76) = 0.190256765666d+05
        epot_c( 77) = -0.150450493013d+05
        epot_c( 78) = -0.920409716245d+05
        epot_c( 79) = -0.216561262057d+04
        epot_c( 80) = -0.171625055996d+05
        epot_c( 81) = 0.153782282485d+05
        epot_c( 82) = 0.158621441332d+04
        epot_c( 83) = -0.785674943671d+04
        epot_c( 84) = -0.209994189557d+05
        epot_c( 85) = -0.333435521700d+05
        epot_c( 86) = -0.224157355756d+05
        epot_c( 87) = 0.182745060537d+05
        epot_c( 88) = -0.501427327882d+04
        epot_c( 89) = 0.735584333269d+04
        epot_c( 90) = 0.413504676578d+03
        epot_c( 91) = 0.109674543527d+04
        epot_c( 92) = -0.186740710649d+04
        epot_c( 93) = 0.359599299302d+04
        epot_c( 94) = 0.721485422686d+04
        epot_c( 95) = -0.330055110073d+05
        epot_c( 96) = 0.730992533723d+05
        epot_c( 97) = -0.799352299005d+04
        epot_c( 98) = 0.112075241137d+05
        epot_c( 99) = -0.156320449278d+05
        epot_c(100) = -0.408715780841d+04
        epot_c(101) = 0.268188694785d+06
        epot_c(102) = -0.102983762167d+06
        epot_c(103) = 0.583699882510d+05
        epot_c(104) = 0.260451292047d+05
        epot_c(105) = -0.177541637884d+05
        epot_c(106) = -0.203359237620d+05
        epot_c(107) = 0.480058134079d+05
        epot_c(108) = -0.885961013155d+05
        epot_c(109) = 0.409859892393d+05
        epot_c(110) = 0.325921040295d+05
        epot_c(111) = -0.260852766104d+05
        epot_c(112) = 0.370242414223d+05
        epot_c(113) = -0.664758905617d+04
        epot_c(114) = -0.275618526134d+05
        epot_c(115) = -0.107805359174d+06
        epot_c(116) = 0.287387764086d+05
        epot_c(117) = 0.666908326522d+05
        epot_c(118) = -0.156541620151d+05
        epot_c(119) = 0.336154826289d+05
        epot_c(120) = -0.414314869617d+05
        epot_c(121) = 0.176679423717d+04
        epot_c(122) = 0.577289611931d+05
        epot_c(123) = 0.302830113555d+05
        epot_c(124) = -0.807830452831d+04
        epot_c(125) = 0.219029154199d+03
        epot_c(126) = -0.413421879006d+05
        epot_c(127) = -0.451799646223d+05
        epot_c(128) = -0.159295124884d+05
        epot_c(129) = -0.424907280397d+05
        epot_c(130) = -0.140611779798d+04
        epot_c(131) = 0.568547939978d+03
        epot_c(132) = 0.264987028493d+05
        epot_c(133) = -0.843157148611d+03
        epot_c(134) = 0.275201914701d+05
        epot_c(135) = -0.105393958865d+05
        epot_c(136) = -0.169105249253d+05
        epot_c(137) = 0.104932391725d+05
        epot_c(138) = 0.119833997133d+05
        epot_c(139) = 0.139893856185d+05
        epot_c(140) = -0.181001837340d+03
        epot_c(141) = 0.153415246742d+05
        epot_c(142) = 0.424481598925d+05
        epot_c(143) = -0.197444461986d+05
        epot_c(144) = 0.490620372926d+04
        epot_c(145) = -0.148185485287d+04
        epot_c(146) = 0.624665998166d+04
        epot_c(147) = 0.245340563235d+05
        epot_c(148) = 0.305015590467d+05
        epot_c(149) = -0.611332793065d+04
        epot_c(150) = -0.555190868114d+04
        epot_c(151) = 0.130155396245d+04
        epot_c(152) = -0.455810219138d+05
        epot_c(153) = 0.115493696632d+05
        epot_c(154) = 0.888993953795d+04
        epot_c(155) = 0.103815845402d+05
        epot_c(156) = 0.562575679728d+04
        epot_c(157) = -0.346648114071d+04
        epot_c(158) = -0.240765749446d+04
        epot_c(159) = 0.206433652571d+04
        epot_c(160) = -0.758049531330d+04
        epot_c(161) = 0.171657661399d+05
        epot_c(162) = -0.260208297109d+04
        epot_c(163) = -0.372114666463d+04
        epot_c(164) = 0.192848747262d+05
        epot_c(165) = -0.614612031387d+05
        epot_c(166) = 0.284938866770d+04
        epot_c(167) = 0.580215877101d+03
        epot_c(168) = -0.142884180958d+05
        epot_c(169) = 0.164345541019d+05
        epot_c(170) = 0.204700489332d+04
        epot_c(171) = 0.232847535995d+06
        epot_c(172) = -0.215988386103d+06
        epot_c(173) = -0.543268059091d+06
        epot_c(174) = -0.291087990782d+06
        epot_c(175) = 0.993671706214d+05
        epot_c(176) = -0.324811444690d+05
        epot_c(177) = -0.165451102346d+06
        epot_c(178) = 0.503736859146d+05
        epot_c(179) = 0.441027718212d+05
        epot_c(180) = 0.427514011421d+05
        epot_c(181) = 0.316984455059d+05
        epot_c(182) = -0.916644985736d+05
        epot_c(183) = -0.497678170749d+05
        epot_c(184) = -0.570481120151d+04
        epot_c(185) = 0.392692240737d+05
        epot_c(186) = -0.636234723414d+05
        epot_c(187) = 0.605343078929d+05
        epot_c(188) = -0.504945371764d+05
        epot_c(189) = -0.394790420088d+05
        epot_c(190) = -0.296012156034d+05
        epot_c(191) = -0.189032641400d+05
        epot_c(192) = 0.315573783036d+05
        epot_c(193) = -0.123145782845d+05
        epot_c(194) = -0.183417710262d+05
        epot_c(195) = 0.137514063372d+04
        epot_c(196) = 0.282584683995d+05
        epot_c(197) = 0.255488952485d+03
        epot_c(198) = -0.214177073532d+05
        epot_c(199) = -0.329389386184d+05
        epot_c(200) = 0.326844593866d+05
        epot_c(201) = -0.268690289613d+05
        epot_c(202) = 0.688737993203d+03
        epot_c(203) = -0.174784217735d+05
        epot_c(204) = 0.592388749822d+04
        epot_c(205) = 0.782006070657d+04
        epot_c(206) = 0.182368489092d+05
        epot_c(207) = -0.195025918004d+05
        epot_c(208) = 0.140155756116d+05
        epot_c(209) = 0.411097619946d+05
        epot_c(210) = 0.363810881096d+03
        epot_c(211) = -0.280890295507d+04
        epot_c(212) = -0.105569797133d+05
        epot_c(213) = 0.380998787463d+04
        epot_c(214) = -0.360964430447d+04
        epot_c(215) = -0.986412341243d+04
        epot_c(216) = 0.226662646459d+05
        epot_c(217) = -0.376753011539d+04
        epot_c(218) = 0.578159426611d+03
        epot_c(219) = -0.121449139228d+05
        epot_c(220) = -0.750765164465d+04
        epot_c(221) = -0.140227007314d+05
        epot_c(222) = 0.203120231276d+04
        epot_c(223) = 0.118811735987d+05
        epot_c(224) = 0.519651005890d+04
        epot_c(225) = 0.947714025376d+04
        epot_c(226) = 0.179018909494d+05
        epot_c(227) = 0.121879281662d+05
        epot_c(228) = -0.126436586574d+05
        epot_c(229) = -0.452354749524d+04
        epot_c(230) = -0.684935604354d+04
        epot_c(231) = -0.174197586324d+04
        epot_c(232) = -0.238302572429d+04
        epot_c(233) = 0.537465564700d+04
        epot_c(234) = -0.157114135550d+03
        epot_c(235) = -0.652324696358d+04
        epot_c(236) = -0.323096366831d+04
        epot_c(237) = 0.296083943992d+04
        epot_c(238) = -0.122452541606d+04
        epot_c(239) = -0.593094058930d+04
        epot_c(240) = -0.379005180042d+04
        epot_c(241) = 0.777078060474d+04
        epot_c(242) = 0.193407389171d+04
        epot_c(243) = -0.251175830196d+04
        epot_c(244) = -0.342507311002d+05
        epot_c(245) = -0.817218571265d+04
        epot_c(246) = 0.167151796682d+04
        epot_c(247) = 0.140032630687d+04
        epot_c(248) = -0.108222853767d+04
        epot_c(249) = -0.289156752055d+04
        epot_c(250) = 0.671923615356d+04
        epot_c(251) = -0.186517211143d+04
        epot_c(252) = -0.608928447848d+04
        epot_c(253) = -0.320885381491d+04
        epot_c(254) = 0.168580430914d+03
        epot_c(255) = 0.561938258944d+03
        epot_c(256) = 0.178839575085d+04
        epot_c(257) = -0.127093039569d+02
        epot_c(258) = 0.158954991882d+05
        epot_c(259) = -0.728671162222d+04
        epot_c(260) = 0.273935428053d+03
        epot_c(261) = 0.210205486354d+03
        epot_c(262) = -0.726829028390d+04
        epot_c(263) = 0.234342832052d+04
        epot_c(264) = 0.450306951911d+03
        epot_c(265) = -0.225605345434d+04
        epot_c(266) = 0.146430478883d+04
        epot_c(267) = -0.502512737258d+04
        epot_c(268) = 0.529554234299d+03
        epot_c(269) = 0.665683392527d+03
        epot_c(270) = -0.446965631596d+04
        epot_c(271) = 0.155396579789d+05
        epot_c(272) = -0.690692642506d+03
        epot_c(273) = 0.128917160396d+04
        epot_c(274) = 0.373090368655d+04
        epot_c(275) = -0.454736355855d+04
        epot_c(276) = -0.377918615671d+03
        v           = 0.d0; g = 0.d0; d = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                j     = 3*(iatom-1)+idir
                tx(j) = x(iatom,idir)
            enddo
        enddo
        tg = 0.d0
        call n4pes(tx, v, tg, igrad)
        v  = v/23.0609 !+Eref*27.211386d0
        tg = tg/23.0609
        do istate=1, nstates
            p(istate) = v
        enddo
        do istate=1, nstates
            do iatom=1, natoms
                do idir=1, 3
                    j                    = 3*(iatom-1)+idir
                    g(istate,iatom,idir) = tg(j)
                enddo
            enddo
        enddo
        d = 0.d0
        contains
        !> @brief Evaluate n4pes for the N4_1A_PIP_2013 potential.
        subroutine n4pes( x, v, dvdx, igrad )
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
            implicit none
            integer          :: i, igrad, j, nob, k
            double precision :: v
            double precision :: dvdx(12), x(12)
            !      common /coord/    R(6)
            !      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
            !      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
            !     $                  dRdX(6,12),dBdR(6,276)
            !      common /msprmt/   a,re
            ! Read Cartesian coordinate from input file
            call coord_convt(x)
            if( igrad .le. 1 )then
                ! Call subroutine Evv to evaluate potential energy V
                call evv(v)
                if( igrad .eq. 1 )then
                    ! Call EvdVdX to evaluate the derivatives of V w.r.t. X
                    call evdvdx(x, dvdx)
                endif
            else
                write (*, *) 'Only energy and gradient are available'
            endif
        end subroutine n4pes
        !> @brief Evaluate coord_convt for the N4_1A_PIP_2013 potential.
        subroutine coord_convt( x )
            !  Program to calculate the six interatomic distance
            !  by reading XYZ coordinate
            implicit none
            integer          :: i
            double precision :: x(12)
            !  Now, calculate the inter-atomic distance
            !  r1 = r(N1N2)    r2 = r(N1N3)
            !  r3 = r(N1N4)    r4 = r(N2N3)
            !  r5 = r(N2N4)    r6 = r(N3N4)
            coord_r(1) = sqrt((x(4)-x(1))**2 + (x(5)-x(2))**2 + (x(6)-x(3))**2)
            coord_r(2) = sqrt((x(7)-x(1))**2 + (x(8)-x(2))**2 + (x(9)-x(3))**2)
            coord_r(3) = sqrt((x(10)-x(1))**2 + (x(11)-x(2))**2 + (x(12)-x(3))**2)
            coord_r(4) = sqrt((x(4)-x(7))**2 + (x(5)-x(8))**2 + (x(6)-x(9))**2)
            coord_r(5) = sqrt((x(4)-x(10))**2 + (x(5)-x(11))**2 + (x(6)-x(12))**2)
            coord_r(6) = sqrt((x(7)-x(10))**2 + (x(8)-x(11))**2 + (x(9)-x(12))**2)
        end subroutine coord_convt
        !> @brief Evaluate EvV for the N4_1A_PIP_2013 potential.
        subroutine evv( v )
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
            implicit none
            integer          :: i, j, k
            double precision :: dist, dv2dr, v, v2
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
            do i=1, 6
                dist = coord_r(i)
                call ev2gm2(dist, v2, dv2dr, 1, 0)
                v = v+v2
            enddo
            ! Evaluate V by taken the product of C and Basis function array
            do i=1, 276
                v = v + epot_c(i)*epot_b(i)
            enddo
            !      Write(*,9999) V
            ! 9999 Format('The potential energy is ',F20.14,' kcal/mol')
        end subroutine evv
        !> @brief Evaluate EvdVdX for the N4_1A_PIP_2013 potential.
        subroutine evdvdx( x, dvdx )
            ! Subroutine to evaluate dRdX for given R and X
            ! R:            R(6), 6 bond lengths
            ! X:            X(12), 12 Cartesian coordinates
            ! rM(0:111):    Array to store monomials
            ! P(0:305):     Array to store polynomials
            ! dVdX:         dVdX(12), derivatives of V w.r.t. Cartesian coordinates
            ! dVdR:         dVdR(6), derivatives of V w.r.t.6 bond length
            ! dRdX:         dRdX(6,12), derivatives of R(6) w.r.t. 12
            !               Cartesian coordinates
            implicit none
            integer          :: i, j
            double precision :: dvdx(12), x(12)
            ! Initialize dVdX
            do i=1, 12
                dvdx(i) = 0.0d0
            enddo
            ! Call EvdVdR to evaluate dVdR(6)
            call evdvdr
            ! Call EvdRdX to evaluate dRdX(6,12)
            call evdrdx(x)
            ! Calculate dVdX by using chain rule: dV/dXi=(dV/dRj)*(dRj/dXi), j=1 to
            do i=1, 12
                do j=1, 6
                    dvdx(i) = dvdx(i) + gradt_dvdr(j)*gradt_drdx(j,i)
                enddo
            enddo
            !      write(*,*) 'The 12 dVdX are:'
            !      Write(*,9999) (dVdX(i),i=1,12)
            ! 9999 Format(1x,3F15.8)
        end subroutine evdvdx
        !> @brief Evaluate EvMorse for the N4_1A_PIP_2013 potential.
        subroutine evmorse
            implicit none
            integer :: i
            ! Morse term ms = exp(-(r-re)/a)
            ! re:   equlibrium bond length(1.098A)
            ! a:    nonlinear parameter optimized for diatomic molecules, unit Anstr
            do i=1, 6
                epot_rms(i) = exp(-(coord_r(i)-msprmt_re)/msprmt_a)
            enddo
        end subroutine evmorse
        !> @brief Evaluate EvMono for the N4_1A_PIP_2013 potential.
        subroutine evmono
            !  The subroutine reads six Morse variables(X) and calculates the
            !  monomials(M) that do not have usable decomposition.
            !  For A4 with max. degree 9, the number of monomials is nom.
            implicit none
            epot_rm(0)   = 1.0d0
            epot_rm(1)   = epot_rms(6)
            epot_rm(2)   = epot_rms(5)
            epot_rm(3)   = epot_rms(4)
            epot_rm(4)   = epot_rms(3)
            epot_rm(5)   = epot_rms(2)
            epot_rm(6)   = epot_rms(1)
            epot_rm(7)   = epot_rm(3)*epot_rm(4)
            epot_rm(8)   = epot_rm(2)*epot_rm(5)
            epot_rm(9)   = epot_rm(1)*epot_rm(6)
            epot_rm(10)  = epot_rm(1)*epot_rm(2)
            epot_rm(11)  = epot_rm(1)*epot_rm(3)
            epot_rm(12)  = epot_rm(2)*epot_rm(3)
            epot_rm(13)  = epot_rm(1)*epot_rm(4)
            epot_rm(14)  = epot_rm(2)*epot_rm(4)
            epot_rm(15)  = epot_rm(1)*epot_rm(5)
            epot_rm(16)  = epot_rm(3)*epot_rm(5)
            epot_rm(17)  = epot_rm(4)*epot_rm(5)
            epot_rm(18)  = epot_rm(2)*epot_rm(6)
            epot_rm(19)  = epot_rm(3)*epot_rm(6)
            epot_rm(20)  = epot_rm(4)*epot_rm(6)
            epot_rm(21)  = epot_rm(5)*epot_rm(6)
            epot_rm(22)  = epot_rm(1)*epot_rm(7)
            epot_rm(23)  = epot_rm(2)*epot_rm(7)
            epot_rm(24)  = epot_rm(1)*epot_rm(8)
            epot_rm(25)  = epot_rm(2)*epot_rm(16)
            epot_rm(26)  = epot_rm(2)*epot_rm(17)
            epot_rm(27)  = epot_rm(3)*epot_rm(17)
            epot_rm(28)  = epot_rm(1)*epot_rm(18)
            epot_rm(29)  = epot_rm(1)*epot_rm(19)
            epot_rm(30)  = epot_rm(1)*epot_rm(20)
            epot_rm(31)  = epot_rm(3)*epot_rm(20)
            epot_rm(32)  = epot_rm(1)*epot_rm(21)
            epot_rm(33)  = epot_rm(2)*epot_rm(21)
            epot_rm(34)  = epot_rm(1)*epot_rm(12)
            epot_rm(35)  = epot_rm(1)*epot_rm(17)
            epot_rm(36)  = epot_rm(2)*epot_rm(20)
            epot_rm(37)  = epot_rm(3)*epot_rm(21)
            epot_rm(38)  = epot_rm(1)*epot_rm(14)
            epot_rm(39)  = epot_rm(1)*epot_rm(16)
            epot_rm(40)  = epot_rm(2)*epot_rm(19)
            epot_rm(41)  = epot_rm(4)*epot_rm(21)
            epot_rm(42)  = epot_rm(2)*epot_rm(27)
            epot_rm(43)  = epot_rm(1)*epot_rm(31)
            epot_rm(44)  = epot_rm(1)*epot_rm(33)
            epot_rm(45)  = epot_rm(1)*epot_rm(23)
            epot_rm(46)  = epot_rm(1)*epot_rm(25)
            epot_rm(47)  = epot_rm(1)*epot_rm(26)
            epot_rm(48)  = epot_rm(1)*epot_rm(27)
            epot_rm(49)  = epot_rm(1)*epot_rm(40)
            epot_rm(50)  = epot_rm(1)*epot_rm(36)
            epot_rm(51)  = epot_rm(2)*epot_rm(31)
            epot_rm(52)  = epot_rm(1)*epot_rm(37)
            epot_rm(53)  = epot_rm(2)*epot_rm(37)
            epot_rm(54)  = epot_rm(1)*epot_rm(41)
            epot_rm(55)  = epot_rm(2)*epot_rm(41)
            epot_rm(56)  = epot_rm(3)*epot_rm(41)
            epot_rm(57)  = epot_rm(1)*epot_rm(42)
            epot_rm(58)  = epot_rm(1)*epot_rm(51)
            epot_rm(59)  = epot_rm(1)*epot_rm(53)
            epot_rm(60)  = epot_rm(1)*epot_rm(55)
            epot_rm(61)  = epot_rm(1)*epot_rm(56)
            epot_rm(62)  = epot_rm(2)*epot_rm(56)
            epot_rm(63)  = epot_rm(1)*epot_rm(62)
            epot_rm(64)  = epot_rm(2)*epot_rm(57)
            epot_rm(65)  = epot_rm(3)*epot_rm(57)
            epot_rm(66)  = epot_rm(4)*epot_rm(57)
            epot_rm(67)  = epot_rm(5)*epot_rm(57)
            epot_rm(68)  = epot_rm(1)*epot_rm(58)
            epot_rm(69)  = epot_rm(3)*epot_rm(58)
            epot_rm(70)  = epot_rm(4)*epot_rm(58)
            epot_rm(71)  = epot_rm(1)*epot_rm(59)
            epot_rm(72)  = epot_rm(2)*epot_rm(59)
            epot_rm(73)  = epot_rm(1)*epot_rm(60)
            epot_rm(74)  = epot_rm(2)*epot_rm(60)
            epot_rm(75)  = epot_rm(1)*epot_rm(61)
            epot_rm(76)  = epot_rm(2)*epot_rm(62)
            epot_rm(77)  = epot_rm(3)*epot_rm(61)
            epot_rm(78)  = epot_rm(3)*epot_rm(62)
            epot_rm(79)  = epot_rm(4)*epot_rm(61)
            epot_rm(80)  = epot_rm(4)*epot_rm(62)
            epot_rm(81)  = epot_rm(5)*epot_rm(59)
            epot_rm(82)  = epot_rm(5)*epot_rm(60)
            epot_rm(83)  = epot_rm(5)*epot_rm(62)
            epot_rm(84)  = epot_rm(6)*epot_rm(58)
            epot_rm(85)  = epot_rm(6)*epot_rm(59)
            epot_rm(86)  = epot_rm(6)*epot_rm(60)
            epot_rm(87)  = epot_rm(6)*epot_rm(61)
            epot_rm(88)  = epot_rm(2)*epot_rm(64)
            epot_rm(89)  = epot_rm(3)*epot_rm(65)
            epot_rm(90)  = epot_rm(4)*epot_rm(66)
            epot_rm(91)  = epot_rm(5)*epot_rm(67)
            epot_rm(92)  = epot_rm(1)*epot_rm(68)
            epot_rm(93)  = epot_rm(3)*epot_rm(69)
            epot_rm(94)  = epot_rm(4)*epot_rm(70)
            epot_rm(95)  = epot_rm(1)*epot_rm(71)
            epot_rm(96)  = epot_rm(2)*epot_rm(72)
            epot_rm(97)  = epot_rm(1)*epot_rm(73)
            epot_rm(98)  = epot_rm(2)*epot_rm(74)
            epot_rm(99)  = epot_rm(1)*epot_rm(75)
            epot_rm(100) = epot_rm(2)*epot_rm(76)
            epot_rm(101) = epot_rm(3)*epot_rm(77)
            epot_rm(102) = epot_rm(3)*epot_rm(78)
            epot_rm(103) = epot_rm(4)*epot_rm(79)
            epot_rm(104) = epot_rm(4)*epot_rm(80)
            epot_rm(105) = epot_rm(5)*epot_rm(81)
            epot_rm(106) = epot_rm(5)*epot_rm(82)
            epot_rm(107) = epot_rm(5)*epot_rm(83)
            epot_rm(108) = epot_rm(6)*epot_rm(84)
            epot_rm(109) = epot_rm(6)*epot_rm(85)
            epot_rm(110) = epot_rm(6)*epot_rm(86)
            epot_rm(111) = epot_rm(6)*epot_rm(87)
        end subroutine evmono
        !> @brief Evaluate EvPoly for the N4_1A_PIP_2013 potential.
        subroutine evpoly
            !  The subroutine reads monomials(m) and calculates the
            !  permutation invariant polynomials(p)
            !  For A4 with max. degree 9, the number of polynomials is nob.
            implicit none
            epot_p(0) = epot_rm(0)
            epot_p(1) = epot_rm(1) + epot_rm(2) + epot_rm(3) + epot_rm(4) + epot_rm(5) + epot_rm(6)
            epot_p(2) = epot_rm(7) + epot_rm(8) + epot_rm(9)
            epot_p(3)  = epot_rm(10) + epot_rm(11) + epot_rm(12) + epot_rm(13) + epot_rm(14) + epot_rm(15) + &
                & epot_rm(16) + epot_rm(17) + epot_rm(18) + epot_rm(19) + epot_rm(20) + epot_rm(21)
            epot_p(4) = epot_p(1)*epot_p(1) - epot_p(3) - epot_p(2) - epot_p(3) - epot_p(2)
            epot_p(5)  = epot_rm(22) + epot_rm(23) + epot_rm(24) + epot_rm(25) + epot_rm(26) + epot_rm(27) + &
                & epot_rm(28) + epot_rm(29) + epot_rm(30) + epot_rm(31) + epot_rm(32) + epot_rm(33)
            epot_p(6) = epot_rm(34) + epot_rm(35) + epot_rm(36) + epot_rm(37)
            epot_p(7) = epot_rm(38) + epot_rm(39) + epot_rm(40) + epot_rm(41)
            epot_p(8) = epot_p(1)*epot_p(2) - epot_p(5)
            epot_p(9)  = epot_p(1)*epot_p(3) - epot_p(6) - epot_p(7) - epot_p(5) - epot_p(6) - epot_p(7) - epot_p(5) - &
                & epot_p(6) - epot_p(7)
            epot_p(10) = epot_p(1)*epot_p(4) - epot_p(9) - epot_p(8)
            epot_p(11) = epot_rm(42) + epot_rm(43) + epot_rm(44)
            epot_p(12) = epot_rm(45) + epot_rm(46) + epot_rm(47) + epot_rm(48) + epot_rm(49) + epot_rm(50) + &
                & epot_rm(51) + epot_rm(52) + epot_rm(53) + epot_rm(54) + epot_rm(55) + epot_rm(56)
            epot_p(13) = epot_p(2)*epot_p(3) - epot_p(12)
            epot_p(14) = epot_p(1)*epot_p(5) - epot_p(12) - epot_p(11) - epot_p(13) - epot_p(12) - epot_p(11) - &
                & epot_p(11) - epot_p(11)
            epot_p(15) = epot_p(1)*epot_p(6) - epot_p(12)
            epot_p(16) = epot_p(1)*epot_p(7) - epot_p(12)
            epot_p(17) = epot_p(2)*epot_p(2) - epot_p(11) - epot_p(11)
            epot_p(18) = epot_p(3)*epot_p(3) - epot_p(12) - epot_p(11) - epot_p(15) - epot_p(16) - epot_p(14) - &
                & epot_p(12) - epot_p(11) - epot_p(15) - epot_p(16) - epot_p(14) - epot_p(12) - epot_p(11) - &
                & epot_p(12) - epot_p(11)
            epot_p(19) = epot_p(2)*epot_p(4) - epot_p(14)
            epot_p(20) = epot_p(3)*epot_p(4) - epot_p(15) - epot_p(16) - epot_p(13)
            epot_p(21) = epot_p(1)*epot_p(10) - epot_p(20) - epot_p(19)
            epot_p(22) = epot_rm(57) + epot_rm(58) + epot_rm(59) + epot_rm(60) + epot_rm(61) + epot_rm(62)
            epot_p(23) = epot_p(1)*epot_p(11) - epot_p(22)
            epot_p(24) = epot_p(2)*epot_p(6)
            epot_p(25) = epot_p(2)*epot_p(7)
            epot_p(26) = epot_p(1)*epot_p(12) - epot_p(22) - epot_p(24) - epot_p(25) - epot_p(22) - epot_p(22) - &
                & epot_p(22)
            epot_p(27) = epot_p(2)*epot_p(5) - epot_p(22) - epot_p(23) - epot_p(22)
            epot_p(28) = epot_p(3)*epot_p(5) - epot_p(22) - epot_p(26) - epot_p(24) - epot_p(25) - epot_p(23) - &
                & epot_p(22) - epot_p(24) - epot_p(25) - epot_p(23) - epot_p(22) - epot_p(22)
            epot_p(29) = epot_p(3)*epot_p(6) - epot_p(22) - epot_p(26) - epot_p(22)
            epot_p(30) = epot_p(3)*epot_p(7) - epot_p(22) - epot_p(26) - epot_p(22)
            epot_p(31) = epot_p(2)*epot_p(9) - epot_p(26) - epot_p(28)
            epot_p(32) = epot_p(1)*epot_p(14) - epot_p(26) - epot_p(23) - epot_p(28)
            epot_p(33) = epot_p(4)*epot_p(6) - epot_p(25)
            epot_p(34) = epot_p(4)*epot_p(7) - epot_p(24)
            epot_p(35) = epot_p(1)*epot_p(17) - epot_p(27)
            epot_p(36) = epot_p(1)*epot_p(18) - epot_p(29) - epot_p(30) - epot_p(28)
            epot_p(37) = epot_p(2)*epot_p(10) - epot_p(32)
            epot_p(38) = epot_p(3)*epot_p(10) - epot_p(33) - epot_p(34) - epot_p(31)
            epot_p(39) = epot_p(1)*epot_p(21) - epot_p(38) - epot_p(37)
            epot_p(40) = epot_rm(63)
            epot_p(41) = epot_rm(64) + epot_rm(65) + epot_rm(66) + epot_rm(67) + epot_rm(68) + epot_rm(69) + &
                & epot_rm(70) + epot_rm(71) + epot_rm(72) + epot_rm(73) + epot_rm(74) + epot_rm(75) + epot_rm(76) + &
                & epot_rm(77) + epot_rm(78) + epot_rm(79) + epot_rm(80) + epot_rm(81) + epot_rm(82) + epot_rm(83) + &
                & epot_rm(84) + epot_rm(85) + epot_rm(86) + epot_rm(87)
            epot_p(42) = epot_p(1)*epot_p(22) - epot_p(40) - epot_p(41) - epot_p(40) - epot_p(40) - epot_p(40) - &
                & epot_p(40) - epot_p(40)
            epot_p(43) = epot_p(2)*epot_p(11) - epot_p(40) - epot_p(40) - epot_p(40)
            epot_p(44) = epot_p(2)*epot_p(12) - epot_p(41)
            epot_p(45) = epot_p(3)*epot_p(11) - epot_p(41)
            epot_p(46) = epot_p(5)*epot_p(6) - epot_p(41)
            epot_p(47) = epot_p(5)*epot_p(7) - epot_p(41)
            epot_p(48) = epot_p(6)*epot_p(7) - epot_p(40) - epot_p(40) - epot_p(40) - epot_p(40)
            epot_p(49) = epot_p(4)*epot_p(11) - epot_p(42)
            epot_p(50) = epot_p(2)*epot_p(15) - epot_p(46)
            epot_p(51) = epot_p(2)*epot_p(16) - epot_p(47)
            epot_p(52) = epot_p(4)*epot_p(12) - epot_p(41) - epot_p(50) - epot_p(51)
            epot_p(53) = epot_p(2)*epot_p(14) - epot_p(42) - epot_p(49) - epot_p(42)
            epot_p(54) = epot_p(6)*epot_p(6) - epot_p(42) - epot_p(42)
            epot_p(55) = epot_p(7)*epot_p(7) - epot_p(42) - epot_p(42)
            epot_p(56) = epot_p(3)*epot_p(17) - epot_p(44)
            epot_p(57) = epot_p(2)*epot_p(18) - epot_p(48)
            epot_p(58) = epot_p(3)*epot_p(14) - epot_p(41) - epot_p(52) - epot_p(46) - epot_p(47) - epot_p(45) - &
                & epot_p(45)
            epot_p(59) = epot_p(6)*epot_p(9) - epot_p(41) - epot_p(52) - epot_p(47)
            epot_p(60) = epot_p(7)*epot_p(9) - epot_p(41) - epot_p(52) - epot_p(46)
            epot_p(61) = epot_p(2)*epot_p(20) - epot_p(52) - epot_p(58)
            epot_p(62) = epot_p(1)*epot_p(32) - epot_p(52) - epot_p(49) - epot_p(58)
            epot_p(63) = epot_p(6)*epot_p(10) - epot_p(51)
            epot_p(64) = epot_p(7)*epot_p(10) - epot_p(50)
            epot_p(65) = epot_p(2)*epot_p(17) - epot_p(43)
            epot_p(66) = epot_p(3)*epot_p(18) - epot_p(46) - epot_p(47) - epot_p(45) - epot_p(59) - epot_p(60) - &
                & epot_p(58)
            epot_p(67) = epot_p(2)*epot_p(19) - epot_p(49)
            epot_p(68) = epot_p(1)*epot_p(36) - epot_p(59) - epot_p(60) - epot_p(58) - epot_p(57) - epot_p(66) - &
                & epot_p(66)
            epot_p(69) = epot_p(2)*epot_p(21) - epot_p(62)
            epot_p(70) = epot_p(3)*epot_p(21) - epot_p(63) - epot_p(64) - epot_p(61)
            epot_p(71) = epot_p(1)*epot_p(39) - epot_p(70) - epot_p(69)
            epot_p(72) = epot_p(40)*epot_p(1)
            epot_p(73) = epot_p(2)*epot_p(22) - epot_p(72)
            epot_p(74) = epot_p(6)*epot_p(11)
            epot_p(75) = epot_p(7)*epot_p(11)
            epot_p(76) = epot_p(3)*epot_p(22) - epot_p(72) - epot_p(74) - epot_p(75) - epot_p(72) - epot_p(72) - &
                & epot_p(72)
            epot_p(77) = epot_rm(88) + epot_rm(89) + epot_rm(90) + epot_rm(91) + epot_rm(92) + epot_rm(93) + &
                & epot_rm(94) + epot_rm(95) + epot_rm(96) + epot_rm(97) + epot_rm(98) + epot_rm(99) + epot_rm(100) + &
                & epot_rm(101) + epot_rm(102) + epot_rm(103) + epot_rm(104) + epot_rm(105) + epot_rm(106) + &
                & epot_rm(107) + epot_rm(108) + epot_rm(109) + epot_rm(110) + epot_rm(111)
            epot_p(78)  = epot_p(1)*epot_p(42) - epot_p(72) - epot_p(76)
            epot_p(79)  = epot_p(5)*epot_p(11) - epot_p(72) - epot_p(73) - epot_p(72)
            epot_p(80)  = epot_p(2)*epot_p(26) - epot_p(76) - epot_p(77)
            epot_p(81)  = epot_p(6)*epot_p(12) - epot_p(72) - epot_p(76) - epot_p(72)
            epot_p(82)  = epot_p(7)*epot_p(12) - epot_p(72) - epot_p(76) - epot_p(72)
            epot_p(83)  = epot_p(8)*epot_p(11) - epot_p(72)
            epot_p(84)  = epot_p(6)*epot_p(17)
            epot_p(85)  = epot_p(7)*epot_p(17)
            epot_p(86)  = epot_p(9)*epot_p(11) - epot_p(76) - epot_p(77)
            epot_p(87)  = epot_p(2)*epot_p(29) - epot_p(81)
            epot_p(88)  = epot_p(6)*epot_p(14) - epot_p(75) - epot_p(75)
            epot_p(89)  = epot_p(2)*epot_p(30) - epot_p(82)
            epot_p(90)  = epot_p(7)*epot_p(14) - epot_p(74) - epot_p(74)
            epot_p(91)  = epot_p(1)*epot_p(48) - epot_p(76) - epot_p(81) - epot_p(82)
            epot_p(92)  = epot_p(10)*epot_p(11) - epot_p(78)
            epot_p(93)  = epot_p(2)*epot_p(33) - epot_p(88)
            epot_p(94)  = epot_p(2)*epot_p(34) - epot_p(90)
            epot_p(95)  = epot_p(10)*epot_p(12) - epot_p(77) - epot_p(93) - epot_p(94)
            epot_p(96)  = epot_p(2)*epot_p(27) - epot_p(73) - epot_p(79)
            epot_p(97)  = epot_p(2)*epot_p(28) - epot_p(76) - epot_p(86)
            epot_p(98)  = epot_p(1)*epot_p(53) - epot_p(80) - epot_p(79) - epot_p(97)
            epot_p(99)  = epot_p(1)*epot_p(54) - epot_p(81)
            epot_p(100) = epot_p(1)*epot_p(55) - epot_p(82)
            epot_p(101) = epot_p(5)*epot_p(18) - epot_p(76) - epot_p(91) - epot_p(87) - epot_p(89) - epot_p(86)
            epot_p(102) = epot_p(6)*epot_p(18) - epot_p(74) - epot_p(90)
            epot_p(103) = epot_p(7)*epot_p(18) - epot_p(75) - epot_p(88)
            epot_p(104) = epot_p(2)*epot_p(31) - epot_p(77) - epot_p(86)
            epot_p(105) = epot_p(2)*epot_p(36) - epot_p(91) - epot_p(101)
            epot_p(106) = epot_p(3)*epot_p(32) - epot_p(77) - epot_p(95) - epot_p(88) - epot_p(90) - epot_p(86)
            epot_p(107) = epot_p(4)*epot_p(29) - epot_p(82) - epot_p(80) - epot_p(99)
            epot_p(108) = epot_p(4)*epot_p(30) - epot_p(81) - epot_p(80) - epot_p(100)
            epot_p(109) = epot_p(2)*epot_p(38) - epot_p(95) - epot_p(106)
            epot_p(110) = epot_p(1)*epot_p(62) - epot_p(95) - epot_p(92) - epot_p(106)
            epot_p(111) = epot_p(6)*epot_p(21) - epot_p(94)
            epot_p(112) = epot_p(7)*epot_p(21) - epot_p(93)
            epot_p(113) = epot_p(1)*epot_p(65) - epot_p(96)
            epot_p(114) = epot_p(1)*epot_p(66) - epot_p(102) - epot_p(103) - epot_p(101)
            epot_p(115) = epot_p(2)*epot_p(37) - epot_p(92)
            epot_p(116) = epot_p(10)*epot_p(18) - epot_p(99) - epot_p(100) - epot_p(97)
            epot_p(117) = epot_p(2)*epot_p(39) - epot_p(110)
            epot_p(118) = epot_p(3)*epot_p(39) - epot_p(111) - epot_p(112) - epot_p(109)
            epot_p(119) = epot_p(1)*epot_p(71) - epot_p(118) - epot_p(117)
            epot_p(120) = epot_p(40)*epot_p(2)
            epot_p(121) = epot_p(40)*epot_p(3)
            epot_p(122) = epot_p(40)*epot_p(4)
            epot_p(123) = epot_p(11)*epot_p(12) - epot_p(121)
            epot_p(124) = epot_p(2)*epot_p(42) - epot_p(122)
            epot_p(125) = epot_p(6)*epot_p(22) - epot_p(121)
            epot_p(126) = epot_p(7)*epot_p(22) - epot_p(121)
            epot_p(127) = epot_p(2)*epot_p(41) - epot_p(121) - epot_p(123) - epot_p(121)
            epot_p(128) = epot_p(6)*epot_p(23) - epot_p(123)
            epot_p(129) = epot_p(7)*epot_p(23) - epot_p(123)
            epot_p(130) = epot_p(3)*epot_p(41) - epot_p(122) - epot_p(121) - epot_p(120) - epot_p(125) - epot_p(128) - &
                & epot_p(126) - epot_p(129) - epot_p(124) - epot_p(123) - epot_p(122) - epot_p(121) - epot_p(120) - &
                & epot_p(125) - epot_p(126) - epot_p(124) - epot_p(123) - epot_p(122) - epot_p(121) - epot_p(120) - &
                & epot_p(122) - epot_p(121) - epot_p(120) - epot_p(120) - epot_p(120) - epot_p(120) - epot_p(120)
            epot_p(131) = epot_p(3)*epot_p(42) - epot_p(121) - epot_p(125) - epot_p(126) - epot_p(121)
            epot_p(132) = epot_p(4)*epot_p(41) - epot_p(121) - epot_p(131) - epot_p(128) - epot_p(129) - epot_p(127) - &
                & epot_p(121)
            epot_p(133) = epot_p(1)*epot_p(78) - epot_p(122) - epot_p(131)
            epot_p(134) = epot_p(11)*epot_p(11) - epot_p(120) - epot_p(120)
            epot_p(135) = epot_p(2)*epot_p(48) - epot_p(130)
            epot_p(136) = epot_p(11)*epot_p(17) - epot_p(120)
            epot_p(137) = epot_p(2)*epot_p(44) - epot_p(123)
            epot_p(138) = epot_p(2)*epot_p(45) - epot_p(121)
            epot_p(139) = epot_p(11)*epot_p(14) - epot_p(122) - epot_p(124) - epot_p(122)
            epot_p(140) = epot_p(6)*epot_p(27) - epot_p(127)
            epot_p(141) = epot_p(2)*epot_p(54)
            epot_p(142) = epot_p(7)*epot_p(27) - epot_p(127)
            epot_p(143) = epot_p(2)*epot_p(52) - epot_p(131) - epot_p(132)
            epot_p(144) = epot_p(1)*epot_p(81) - epot_p(125) - epot_p(141) - epot_p(135) - epot_p(125)
            epot_p(145) = epot_p(2)*epot_p(55)
            epot_p(146) = epot_p(1)*epot_p(82) - epot_p(126) - epot_p(135) - epot_p(145) - epot_p(126)
            epot_p(147) = epot_p(11)*epot_p(18) - epot_p(130)
            epot_p(148) = epot_p(6)*epot_p(28) - epot_p(129) - epot_p(123) - epot_p(143)
            epot_p(149) = epot_p(7)*epot_p(28) - epot_p(128) - epot_p(123) - epot_p(143)
            epot_p(150) = epot_p(6)*epot_p(30) - epot_p(121) - epot_p(146)
            epot_p(151) = epot_p(11)*epot_p(19) - epot_p(122)
            epot_p(152) = epot_p(2)*epot_p(50) - epot_p(128)
            epot_p(153) = epot_p(2)*epot_p(51) - epot_p(129)
            epot_p(154) = epot_p(11)*epot_p(20) - epot_p(131) - epot_p(132)
            epot_p(155) = epot_p(2)*epot_p(59) - epot_p(144) - epot_p(148)
            epot_p(156) = epot_p(6)*epot_p(32) - epot_p(129)
            epot_p(157) = epot_p(2)*epot_p(60) - epot_p(146) - epot_p(149)
            epot_p(158) = epot_p(7)*epot_p(32) - epot_p(128)
            epot_p(159) = epot_p(6)*epot_p(34) - epot_p(122) - epot_p(145) - epot_p(122)
            epot_p(160) = epot_p(11)*epot_p(21) - epot_p(133)
            epot_p(161) = epot_p(2)*epot_p(63) - epot_p(156)
            epot_p(162) = epot_p(2)*epot_p(64) - epot_p(158)
            epot_p(163) = epot_p(12)*epot_p(21) - epot_p(132) - epot_p(161) - epot_p(162)
            epot_p(164) = epot_p(2)*epot_p(53) - epot_p(124) - epot_p(139)
            epot_p(165) = epot_p(2)*epot_p(58) - epot_p(131) - epot_p(154)
            epot_p(166) = epot_p(3)*epot_p(54) - epot_p(125) - epot_p(144)
            epot_p(167) = epot_p(3)*epot_p(55) - epot_p(126) - epot_p(146)
            epot_p(168) = epot_p(3)*epot_p(65) - epot_p(137)
            epot_p(169) = epot_p(17)*epot_p(18) - epot_p(135)
            epot_p(170) = epot_p(1)*epot_p(98) - epot_p(143) - epot_p(139) - epot_p(165)
            epot_p(171) = epot_p(4)*epot_p(54) - epot_p(135)
            epot_p(172) = epot_p(4)*epot_p(55) - epot_p(135)
            epot_p(173) = epot_p(2)*epot_p(66) - epot_p(150)
            epot_p(174) = epot_p(1)*epot_p(101) - epot_p(148) - epot_p(149) - epot_p(147) - epot_p(173) - &
                & epot_p(165) - epot_p(147)
            epot_p(175) = epot_p(1)*epot_p(102) - epot_p(150) - epot_p(148) - epot_p(166)
            epot_p(176) = epot_p(1)*epot_p(103) - epot_p(150) - epot_p(149) - epot_p(167)
            epot_p(177) = epot_p(2)*epot_p(61) - epot_p(132) - epot_p(154)
            epot_p(178) = epot_p(2)*epot_p(68) - epot_p(159) - epot_p(174)
            epot_p(179) = epot_p(3)*epot_p(62) - epot_p(132) - epot_p(163) - epot_p(156) - epot_p(158) - epot_p(154)
            epot_p(180) = epot_p(6)*epot_p(38) - epot_p(132) - epot_p(163) - epot_p(157)
            epot_p(181) = epot_p(7)*epot_p(38) - epot_p(132) - epot_p(163) - epot_p(155)
            epot_p(182) = epot_p(2)*epot_p(70) - epot_p(163) - epot_p(179)
            epot_p(183) = epot_p(1)*epot_p(110) - epot_p(163) - epot_p(160) - epot_p(179)
            epot_p(184) = epot_p(6)*epot_p(39) - epot_p(162)
            epot_p(185) = epot_p(7)*epot_p(39) - epot_p(161)
            epot_p(186) = epot_p(2)*epot_p(65) - epot_p(136)
            epot_p(187) = epot_p(3)*epot_p(66) - epot_p(148) - epot_p(149) - epot_p(147) - epot_p(175) - epot_p(176) - &
                & epot_p(174)
            epot_p(188) = epot_p(2)*epot_p(67) - epot_p(151)
            epot_p(189) = epot_p(4)*epot_p(66) - epot_p(166) - epot_p(167) - epot_p(165)
            epot_p(190) = epot_p(2)*epot_p(69) - epot_p(160)
            epot_p(191) = epot_p(18)*epot_p(21) - epot_p(171) - epot_p(172) - epot_p(169)
            epot_p(192) = epot_p(2)*epot_p(71) - epot_p(183)
            epot_p(193) = epot_p(3)*epot_p(71) - epot_p(184) - epot_p(185) - epot_p(182)
            epot_p(194) = epot_p(1)*epot_p(119) - epot_p(193) - epot_p(192)
            epot_p(195) = epot_p(40)*epot_p(5)
            epot_p(196) = epot_p(40)*epot_p(6)
            epot_p(197) = epot_p(40)*epot_p(7)
            epot_p(198) = epot_p(40)*epot_p(8)
            epot_p(199) = epot_p(40)*epot_p(9)
            epot_p(200) = epot_p(40)*epot_p(10)
            epot_p(201) = epot_p(11)*epot_p(22) - epot_p(195)
            epot_p(202) = epot_p(12)*epot_p(22) - epot_p(196) - epot_p(197) - epot_p(195) - epot_p(196) - &
                & epot_p(197) - epot_p(195) - epot_p(196) - epot_p(197)
            epot_p(203) = epot_p(17)*epot_p(22) - epot_p(198)
            epot_p(204) = epot_p(6)*epot_p(43)
            epot_p(205) = epot_p(7)*epot_p(43)
            epot_p(206) = epot_p(11)*epot_p(26) - epot_p(199) - epot_p(202)
            epot_p(207) = epot_p(2)*epot_p(76) - epot_p(199) - epot_p(202)
            epot_p(208) = epot_p(2)*epot_p(78) - epot_p(200)
            epot_p(209) = epot_p(6)*epot_p(41) - epot_p(199) - epot_p(195) - epot_p(202) - epot_p(195)
            epot_p(210) = epot_p(6)*epot_p(42) - epot_p(197) - epot_p(197) - epot_p(197)
            epot_p(211) = epot_p(7)*epot_p(41) - epot_p(199) - epot_p(195) - epot_p(202) - epot_p(195)
            epot_p(212) = epot_p(7)*epot_p(42) - epot_p(196) - epot_p(196) - epot_p(196)
            epot_p(213) = epot_p(11)*epot_p(29) - epot_p(209)
            epot_p(214) = epot_p(11)*epot_p(30) - epot_p(211)
            epot_p(215) = epot_p(18)*epot_p(22) - epot_p(199) - epot_p(213) - epot_p(214)
            epot_p(216) = epot_p(2)*epot_p(77) - epot_p(199) - epot_p(206)
            epot_p(217) = epot_p(6)*epot_p(49) - epot_p(205)
            epot_p(218) = epot_p(7)*epot_p(49) - epot_p(204)
            epot_p(219) = epot_p(3)*epot_p(77) - epot_p(200) - epot_p(199) - epot_p(198) - epot_p(209) - epot_p(217) - &
                & epot_p(211) - epot_p(218) - epot_p(207) - epot_p(204) - epot_p(205) - epot_p(200) - epot_p(199) - &
                & epot_p(198) - epot_p(200) - epot_p(198) - epot_p(200) - epot_p(198)
            epot_p(220) = epot_p(3)*epot_p(78) - epot_p(199) - epot_p(210) - epot_p(212)
            epot_p(221) = epot_p(10)*epot_p(41) - epot_p(199) - epot_p(220) - epot_p(217) - epot_p(218) - epot_p(216)
            epot_p(222) = epot_p(1)*epot_p(133) - epot_p(200) - epot_p(220)
            epot_p(223) = epot_p(11)*epot_p(27) - epot_p(195) - epot_p(203)
            epot_p(224) = epot_p(1)*epot_p(134) - epot_p(201)
            epot_p(225) = epot_p(2)*epot_p(81) - epot_p(209)
            epot_p(226) = epot_p(2)*epot_p(80) - epot_p(202) - epot_p(206)
            epot_p(227) = epot_p(2)*epot_p(82) - epot_p(211)
            epot_p(228) = epot_p(2)*epot_p(91) - epot_p(215) - epot_p(219)
            epot_p(229) = epot_p(11)*epot_p(28) - epot_p(199) - epot_p(207)
            epot_p(230) = epot_p(6)*epot_p(53) - epot_p(205)
            epot_p(231) = epot_p(5)*epot_p(54) - epot_p(209)
            epot_p(232) = epot_p(7)*epot_p(53) - epot_p(204)
            epot_p(233) = epot_p(7)*epot_p(54) - epot_p(196)
            epot_p(234) = epot_p(5)*epot_p(55) - epot_p(211)
            epot_p(235) = epot_p(6)*epot_p(55) - epot_p(197)
            epot_p(236) = epot_p(11)*epot_p(35) - epot_p(198)
            epot_p(237) = epot_p(6)*epot_p(65)
            epot_p(238) = epot_p(7)*epot_p(65)
            epot_p(239) = epot_p(2)*epot_p(86) - epot_p(199) - epot_p(229)
            epot_p(240) = epot_p(1)*epot_p(139) - epot_p(206) - epot_p(229) - epot_p(224)
            epot_p(241) = epot_p(17)*epot_p(29) - epot_p(225)
            epot_p(242) = epot_p(2)*epot_p(99) - epot_p(231)
            epot_p(243) = epot_p(17)*epot_p(30) - epot_p(227)
            epot_p(244) = epot_p(2)*epot_p(95) - epot_p(220) - epot_p(221)
            epot_p(245) = epot_p(4)*epot_p(81) - epot_p(202) - epot_p(242) - epot_p(227)
            epot_p(246) = epot_p(2)*epot_p(100) - epot_p(234)
            epot_p(247) = epot_p(4)*epot_p(82) - epot_p(202) - epot_p(225) - epot_p(246)
            epot_p(248) = epot_p(11)*epot_p(36) - epot_p(215) - epot_p(219)
            epot_p(249) = epot_p(2)*epot_p(102) - epot_p(233)
            epot_p(250) = epot_p(6)*epot_p(58) - epot_p(206) - epot_p(214) - epot_p(244) - epot_p(214)
            epot_p(251) = epot_p(2)*epot_p(103) - epot_p(235)
            epot_p(252) = epot_p(7)*epot_p(58) - epot_p(213) - epot_p(206) - epot_p(244) - epot_p(213)
            epot_p(253) = epot_p(1)*epot_p(150) - epot_p(215) - epot_p(233) - epot_p(235)
            epot_p(254) = epot_p(11)*epot_p(37) - epot_p(200)
            epot_p(255) = epot_p(2)*epot_p(93) - epot_p(217)
            epot_p(256) = epot_p(2)*epot_p(94) - epot_p(218)
            epot_p(257) = epot_p(11)*epot_p(38) - epot_p(220) - epot_p(221)
            epot_p(258) = epot_p(2)*epot_p(107) - epot_p(245) - epot_p(250)
            epot_p(259) = epot_p(6)*epot_p(62) - epot_p(218)
            epot_p(260) = epot_p(2)*epot_p(108) - epot_p(247) - epot_p(252)
            epot_p(261) = epot_p(7)*epot_p(62) - epot_p(217)
            epot_p(262) = epot_p(6)*epot_p(64) - epot_p(200) - epot_p(246) - epot_p(200)
            epot_p(263) = epot_p(11)*epot_p(39) - epot_p(222)
            epot_p(264) = epot_p(2)*epot_p(111) - epot_p(259)
            epot_p(265) = epot_p(2)*epot_p(112) - epot_p(261)
            epot_p(266) = epot_p(12)*epot_p(39) - epot_p(221) - epot_p(264) - epot_p(265)
            epot_p(267) = epot_p(2)*epot_p(98) - epot_p(208) - epot_p(240)
            epot_p(268) = epot_p(6)*epot_p(54) - epot_p(210)
            epot_p(269) = epot_p(7)*epot_p(55) - epot_p(212)
            epot_p(270) = epot_p(2)*epot_p(96) - epot_p(203) - epot_p(223)
            epot_p(271) = epot_p(2)*epot_p(97) - epot_p(207) - epot_p(229)
            epot_p(272) = epot_p(2)*epot_p(101) - epot_p(215) - epot_p(248)
            epot_p(273) = epot_p(2)*epot_p(106) - epot_p(220) - epot_p(257)
            epot_p(274) = epot_p(6)*epot_p(59) - epot_p(220) - epot_p(215) - epot_p(209)
            epot_p(275) = epot_p(7)*epot_p(60) - epot_p(220) - epot_p(215) - epot_p(211)
            epot_p(276) = epot_p(5)*epot_p(66) - epot_p(215) - epot_p(253) - epot_p(249) - epot_p(251) - epot_p(248)
            epot_p(277) = epot_p(6)*epot_p(66) - epot_p(213) - epot_p(252)
            epot_p(278) = epot_p(7)*epot_p(66) - epot_p(214) - epot_p(250)
            epot_p(279) = epot_p(2)*epot_p(104) - epot_p(216) - epot_p(239)
            epot_p(280) = epot_p(2)*epot_p(105) - epot_p(219) - epot_p(248)
            epot_p(281) = epot_p(1)*epot_p(170) - epot_p(244) - epot_p(240) - epot_p(273)
            epot_p(282) = epot_p(10)*epot_p(54) - epot_p(227)
            epot_p(283) = epot_p(10)*epot_p(55) - epot_p(225)
            epot_p(284) = epot_p(2)*epot_p(114) - epot_p(253) - epot_p(276)
            epot_p(285) = epot_p(1)*epot_p(174) - epot_p(250) - epot_p(252) - epot_p(248) - epot_p(276) - epot_p(273)
            epot_p(286) = epot_p(6)*epot_p(68) - epot_p(217) - epot_p(261) - epot_p(251)
            epot_p(287) = epot_p(7)*epot_p(68) - epot_p(218) - epot_p(259) - epot_p(249)
            epot_p(288) = epot_p(2)*epot_p(109) - epot_p(221) - epot_p(257)
            epot_p(289) = epot_p(2)*epot_p(116) - epot_p(262) - epot_p(285)
            epot_p(290) = epot_p(3)*epot_p(110) - epot_p(221) - epot_p(266) - epot_p(259) - epot_p(261) - epot_p(257)
            epot_p(291) = epot_p(6)*epot_p(70) - epot_p(221) - epot_p(266) - epot_p(260)
            epot_p(292) = epot_p(7)*epot_p(70) - epot_p(221) - epot_p(266) - epot_p(258)
            epot_p(293) = epot_p(2)*epot_p(118) - epot_p(266) - epot_p(290)
            epot_p(294) = epot_p(1)*epot_p(183) - epot_p(266) - epot_p(263) - epot_p(290)
            epot_p(295) = epot_p(6)*epot_p(71) - epot_p(265)
            epot_p(296) = epot_p(7)*epot_p(71) - epot_p(264)
            epot_p(297) = epot_p(1)*epot_p(186) - epot_p(270)
            epot_p(298) = epot_p(1)*epot_p(187) - epot_p(277) - epot_p(278) - epot_p(276)
            epot_p(299) = epot_p(2)*epot_p(115) - epot_p(254)
            epot_p(300) = epot_p(1)*epot_p(189) - epot_p(286) - epot_p(287) - epot_p(285) - epot_p(284) - epot_p(298)
            epot_p(301) = epot_p(2)*epot_p(117) - epot_p(263)
            epot_p(302) = epot_p(18)*epot_p(39) - epot_p(282) - epot_p(283) - epot_p(280)
            epot_p(303) = epot_p(2)*epot_p(119) - epot_p(294)
            epot_p(304) = epot_p(3)*epot_p(119) - epot_p(295) - epot_p(296) - epot_p(293)
            epot_p(305) = epot_p(1)*epot_p(194) - epot_p(304) - epot_p(303)
        end subroutine evpoly
        !> @brief Evaluate evbas for the N4_1A_PIP_2013 potential.
        subroutine evbas
            !  The subroutine eliminate the 2-body terms in Bowman's approach
            implicit none
            integer          :: i
            double precision :: b1(306)
            ! Pass P(0:305) to BM1(1:306)
            do i=1, 306
                b1(i) = epot_p(i-1)
            enddo
            ! Remove unconnected terms and 2-body terms and pass to B(1:276)
            epot_b(1) = b1(4)
            do i=2, 4
                epot_b(i) = b1(i+4)
            enddo
            epot_b(5) = b1(10)
            do i=6, 11
                epot_b(i) = b1(i+6)
            enddo
            epot_b(12) = b1(19)
            epot_b(13) = b1(21)
            do i=14, 26
                epot_b(i) = b1(i+9)
            enddo
            epot_b(27) = b1(37)
            epot_b(28) = b1(39)
            do i=29, 53
                epot_b(i) = b1(i+12)
            enddo
            epot_b(54) = b1(67)
            epot_b(55) = b1(69)
            epot_b(56) = b1(71)
            do i=57, 97
                epot_b(i) = b1(i+16)
            enddo
            epot_b(98)  = b1(115)
            epot_b(99)  = b1(117)
            epot_b(100) = b1(119)
            do i=101, 166
                epot_b(i) = b1(i+20)
            enddo
            epot_b(167) = b1(188)
            epot_b(168) = b1(190)
            epot_b(169) = b1(192)
            epot_b(170) = b1(194)
            do i=171, 272
                epot_b(i) = b1(i+25)
            enddo
            epot_b(273) = b1(299)
            epot_b(274) = b1(301)
            epot_b(275) = b1(303)
            epot_b(276) = b1(305)
        end subroutine evbas
        !> @brief Evaluate ev2gm2 for the N4_1A_PIP_2013 potential.
        subroutine ev2gm2( r, v2, dv2dr, imol, igrad )
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
            implicit none
            double precision :: dfdy, v
            integer          :: igrad, imol
            double precision :: c(0:10)
            double precision :: re, de
            double precision :: fy, dfdr, r, u, v2, dv2dr, y, dydr
            if( imol .eq. 1 )then
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
            else if( imol .eq. 2 )then
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
            else if( imol .eq. 3 )then
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
            else if( imol .eq. 4 )then
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
            y  = (r*r*r*r - re*re*re*re)/(r*r*r*r + re*re*re*re)
            fy = c(0) + c(1)*y + c(2)*y*y + c(3)*y*y*y + c(4)*y*y*y*y + c(5)*y*y*y*y*y + c(6)*y*y*y*y*y*y
            u  = exp(-fy*(r-re))
            v2 = de*(1.0d0-u)*(1.0d0-u)-de
            if( igrad .eq. 1 )then
                dfdy  = c(1) + 2.0d0*c(2)*y + 3.0d0*c(3)*y*y + 4.0d0*c(4)*y*y*y + 5.0d0*c(5)*y*y*y*y + 6.0d0*c(6)*y*y* &
                    & y*y*y
                v     = r*r*r*r + re*re*re*re
                dydr  = 8.0d0*r*r*r*re*re*re*re/(v*v)
                dfdr  = dfdy*dydr
                dv2dr = 2.0d0*de*(1-u)*u*(dfdr*(r-re)+fy)
            endif
        end subroutine ev2gm2
        !> @brief Evaluate EvdVdR for the N4_1A_PIP_2013 potential.
        subroutine evdvdr
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
            implicit none
            integer          :: i, j
            double precision :: dist, v2, dv2dr
            ! Initialize dVdR(6)
            do i=1, 6
                gradt_dvdr(i) = 0.0d0
            enddo
            ! Add dV2dR(i) to dVdR
            do i=1, 6
                dist = coord_r(i)
                call ev2gm2(dist, v2, dv2dr, 1, 1)
                gradt_dvdr(i) = dv2dr
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
            do i=1, 6
                do j=1, 276
                    gradt_dvdr(i) = gradt_dvdr(i) + epot_c(j)*gradt_dbdr(i,j)
                enddo
            enddo
        end subroutine evdvdr
        !> @brief Evaluate EvdRdX for the N4_1A_PIP_2013 potential.
        subroutine evdrdx( x )
            ! Subroutine to evaluate dRdX for given R and X
            ! R:            R(6), 6 bond lengths
            ! X:            X(12), 12 Cartesian coordinates
            !
            ! dMdR:         dMdR(6,nom), nom monomials w.r.t.6 bond length
            ! dPdR:         dPdR(6,nob), nop polynomial basis functions
            !               w.r.t. 6 bond length
            ! M(nom):       Array to store monomials
            ! P(nob):       Array to store polynomials
            implicit none
            integer          :: i, j
            double precision :: x(12)
            ! Initialize dRdX(6,12)
            do i=1, 6
                do j=1, 12
                    gradt_drdx(i,j) = 0.0d0
                enddo
            enddo
            ! Start to calculate the non-zero dRdX
            ! dr1dx
            gradt_drdx(1,1) = (x(1)-x(4))/coord_r(1)
            gradt_drdx(1,2) = (x(2)-x(5))/coord_r(1)
            gradt_drdx(1,3) = (x(3)-x(6))/coord_r(1)
            gradt_drdx(1,4) = -gradt_drdx(1,1)
            gradt_drdx(1,5) = -gradt_drdx(1,2)
            gradt_drdx(1,6) = -gradt_drdx(1,3)
            ! dr2dx
            gradt_drdx(2,1) = (x(1)-x(7))/coord_r(2)
            gradt_drdx(2,2) = (x(2)-x(8))/coord_r(2)
            gradt_drdx(2,3) = (x(3)-x(9))/coord_r(2)
            gradt_drdx(2,7) = -gradt_drdx(2,1)
            gradt_drdx(2,8) = -gradt_drdx(2,2)
            gradt_drdx(2,9) = -gradt_drdx(2,3)
            ! dr3dx
            gradt_drdx(3,1)  = (x(1)-x(10))/coord_r(3)
            gradt_drdx(3,2)  = (x(2)-x(11))/coord_r(3)
            gradt_drdx(3,3)  = (x(3)-x(12))/coord_r(3)
            gradt_drdx(3,10) = -gradt_drdx(3,1)
            gradt_drdx(3,11) = -gradt_drdx(3,2)
            gradt_drdx(3,12) = -gradt_drdx(3,3)
            ! dr4dx
            gradt_drdx(4,4) = (x(4)-x(7))/coord_r(4)
            gradt_drdx(4,5) = (x(5)-x(8))/coord_r(4)
            gradt_drdx(4,6) = (x(6)-x(9))/coord_r(4)
            gradt_drdx(4,7) = -gradt_drdx(4,4)
            gradt_drdx(4,8) = -gradt_drdx(4,5)
            gradt_drdx(4,9) = -gradt_drdx(4,6)
            ! dr5dx
            gradt_drdx(5,4)  = (x(4)-x(10))/coord_r(5)
            gradt_drdx(5,5)  = (x(5)-x(11))/coord_r(5)
            gradt_drdx(5,6)  = (x(6)-x(12))/coord_r(5)
            gradt_drdx(5,10) = -gradt_drdx(5,4)
            gradt_drdx(5,11) = -gradt_drdx(5,5)
            gradt_drdx(5,12) = -gradt_drdx(5,6)
            ! dr6dx
            gradt_drdx(6,7)  = (x(7)-x(10))/coord_r(6)
            gradt_drdx(6,8)  = (x(8)-x(11))/coord_r(6)
            gradt_drdx(6,9)  = (x(9)-x(12))/coord_r(6)
            gradt_drdx(6,10) = -gradt_drdx(6,7)
            gradt_drdx(6,11) = -gradt_drdx(6,8)
            gradt_drdx(6,12) = -gradt_drdx(6,9)
            ! Finish the calculation of non-zero dRdX
        end subroutine evdrdx
        !> @brief Evaluate EvdMsdR for the N4_1A_PIP_2013 potential.
        subroutine evdmsdr
            ! Subroutine to evaluate the derivatives of Morse term X
            ! w.r.t. interatomic distance R(6)
            ! dmsdR:        Local variables, dirm(6,6)
            ! a:            Nonlinear pamameter(Angstrom)
            ! re:           equilibrium bond length(Angstrom)
            implicit none
            integer :: i, j
            ! Initialize dmsdr
            do i=1, 6
                do j=1, 6
                    gradt_dmsdr(i,j) = 0.0d0
                enddo
            enddo
            !
            ! Morse term dmsdr = exp(-(r-re)/a)
            ! dmsdr(i,j)=0  i!=j
            ! dmsdr(i,i)= -(1/a)*Exp(-r(i)-re)/a)
            !
            do i=1, 6
                gradt_dmsdr(i,i) = -(1/msprmt_a)*exp(-(coord_r(i)-msprmt_re)/msprmt_a)
            enddo
        end subroutine evdmsdr
        !> @brief Evaluate EvdMdR for the N4_1A_PIP_2013 potential.
        subroutine evdmdr
            !  The subroutine reads M(nom) and dMSdR(6,6) and calculates the
            !  dMdR(6,nom) that do not have usable decomposition.
            !  For A4 with max. degree 9, the number of monomials is nom.
            implicit none
            integer :: i
            do i=1, 6
                gradt_dmdr(i,0)   = 0.0d0
                gradt_dmdr(i,1)   = gradt_dmsdr(i,6)
                gradt_dmdr(i,2)   = gradt_dmsdr(i,5)
                gradt_dmdr(i,3)   = gradt_dmsdr(i,4)
                gradt_dmdr(i,4)   = gradt_dmsdr(i,3)
                gradt_dmdr(i,5)   = gradt_dmsdr(i,2)
                gradt_dmdr(i,6)   = gradt_dmsdr(i,1)
                gradt_dmdr(i,7)   = gradt_dmdr(i,3)*epot_rm(4) + epot_rm(3)*gradt_dmdr(i,4)
                gradt_dmdr(i,8)   = gradt_dmdr(i,2)*epot_rm(5) + epot_rm(2)*gradt_dmdr(i,5)
                gradt_dmdr(i,9)   = gradt_dmdr(i,1)*epot_rm(6) + epot_rm(1)*gradt_dmdr(i,6)
                gradt_dmdr(i,10)  = gradt_dmdr(i,1)*epot_rm(2) + epot_rm(1)*gradt_dmdr(i,2)
                gradt_dmdr(i,11)  = gradt_dmdr(i,1)*epot_rm(3) + epot_rm(1)*gradt_dmdr(i,3)
                gradt_dmdr(i,12)  = gradt_dmdr(i,2)*epot_rm(3) + epot_rm(2)*gradt_dmdr(i,3)
                gradt_dmdr(i,13)  = gradt_dmdr(i,1)*epot_rm(4) + epot_rm(1)*gradt_dmdr(i,4)
                gradt_dmdr(i,14)  = gradt_dmdr(i,2)*epot_rm(4) + epot_rm(2)*gradt_dmdr(i,4)
                gradt_dmdr(i,15)  = gradt_dmdr(i,1)*epot_rm(5) + epot_rm(1)*gradt_dmdr(i,5)
                gradt_dmdr(i,16)  = gradt_dmdr(i,3)*epot_rm(5) + epot_rm(3)*gradt_dmdr(i,5)
                gradt_dmdr(i,17)  = gradt_dmdr(i,4)*epot_rm(5) + epot_rm(4)*gradt_dmdr(i,5)
                gradt_dmdr(i,18)  = gradt_dmdr(i,2)*epot_rm(6) + epot_rm(2)*gradt_dmdr(i,6)
                gradt_dmdr(i,19)  = gradt_dmdr(i,3)*epot_rm(6) + epot_rm(3)*gradt_dmdr(i,6)
                gradt_dmdr(i,20)  = gradt_dmdr(i,4)*epot_rm(6) + epot_rm(4)*gradt_dmdr(i,6)
                gradt_dmdr(i,21)  = gradt_dmdr(i,5)*epot_rm(6) + epot_rm(5)*gradt_dmdr(i,6)
                gradt_dmdr(i,22)  = gradt_dmdr(i,1)*epot_rm(7) + epot_rm(1)*gradt_dmdr(i,7)
                gradt_dmdr(i,23)  = gradt_dmdr(i,2)*epot_rm(7) + epot_rm(2)*gradt_dmdr(i,7)
                gradt_dmdr(i,24)  = gradt_dmdr(i,1)*epot_rm(8) + epot_rm(1)*gradt_dmdr(i,8)
                gradt_dmdr(i,25)  = gradt_dmdr(i,2)*epot_rm(16) + epot_rm(2)*gradt_dmdr(i,16)
                gradt_dmdr(i,26)  = gradt_dmdr(i,2)*epot_rm(17) + epot_rm(2)*gradt_dmdr(i,17)
                gradt_dmdr(i,27)  = gradt_dmdr(i,3)*epot_rm(17) + epot_rm(3)*gradt_dmdr(i,17)
                gradt_dmdr(i,28)  = gradt_dmdr(i,1)*epot_rm(18) + epot_rm(1)*gradt_dmdr(i,18)
                gradt_dmdr(i,29)  = gradt_dmdr(i,1)*epot_rm(19) + epot_rm(1)*gradt_dmdr(i,19)
                gradt_dmdr(i,30)  = gradt_dmdr(i,1)*epot_rm(20) + epot_rm(1)*gradt_dmdr(i,20)
                gradt_dmdr(i,31)  = gradt_dmdr(i,3)*epot_rm(20) + epot_rm(3)*gradt_dmdr(i,20)
                gradt_dmdr(i,32)  = gradt_dmdr(i,1)*epot_rm(21) + epot_rm(1)*gradt_dmdr(i,21)
                gradt_dmdr(i,33)  = gradt_dmdr(i,2)*epot_rm(21) + epot_rm(2)*gradt_dmdr(i,21)
                gradt_dmdr(i,34)  = gradt_dmdr(i,1)*epot_rm(12) + epot_rm(1)*gradt_dmdr(i,12)
                gradt_dmdr(i,35)  = gradt_dmdr(i,1)*epot_rm(17) + epot_rm(1)*gradt_dmdr(i,17)
                gradt_dmdr(i,36)  = gradt_dmdr(i,2)*epot_rm(20) + epot_rm(2)*gradt_dmdr(i,20)
                gradt_dmdr(i,37)  = gradt_dmdr(i,3)*epot_rm(21) + epot_rm(3)*gradt_dmdr(i,21)
                gradt_dmdr(i,38)  = gradt_dmdr(i,1)*epot_rm(14) + epot_rm(1)*gradt_dmdr(i,14)
                gradt_dmdr(i,39)  = gradt_dmdr(i,1)*epot_rm(16) + epot_rm(1)*gradt_dmdr(i,16)
                gradt_dmdr(i,40)  = gradt_dmdr(i,2)*epot_rm(19) + epot_rm(2)*gradt_dmdr(i,19)
                gradt_dmdr(i,41)  = gradt_dmdr(i,4)*epot_rm(21) + epot_rm(4)*gradt_dmdr(i,21)
                gradt_dmdr(i,42)  = gradt_dmdr(i,2)*epot_rm(27) + epot_rm(2)*gradt_dmdr(i,27)
                gradt_dmdr(i,43)  = gradt_dmdr(i,1)*epot_rm(31) + epot_rm(1)*gradt_dmdr(i,31)
                gradt_dmdr(i,44)  = gradt_dmdr(i,1)*epot_rm(33) + epot_rm(1)*gradt_dmdr(i,33)
                gradt_dmdr(i,45)  = gradt_dmdr(i,1)*epot_rm(23) + epot_rm(1)*gradt_dmdr(i,23)
                gradt_dmdr(i,46)  = gradt_dmdr(i,1)*epot_rm(25) + epot_rm(1)*gradt_dmdr(i,25)
                gradt_dmdr(i,47)  = gradt_dmdr(i,1)*epot_rm(26) + epot_rm(1)*gradt_dmdr(i,26)
                gradt_dmdr(i,48)  = gradt_dmdr(i,1)*epot_rm(27) + epot_rm(1)*gradt_dmdr(i,27)
                gradt_dmdr(i,49)  = gradt_dmdr(i,1)*epot_rm(40) + epot_rm(1)*gradt_dmdr(i,40)
                gradt_dmdr(i,50)  = gradt_dmdr(i,1)*epot_rm(36) + epot_rm(1)*gradt_dmdr(i,36)
                gradt_dmdr(i,51)  = gradt_dmdr(i,2)*epot_rm(31) + epot_rm(2)*gradt_dmdr(i,31)
                gradt_dmdr(i,52)  = gradt_dmdr(i,1)*epot_rm(37) + epot_rm(1)*gradt_dmdr(i,37)
                gradt_dmdr(i,53)  = gradt_dmdr(i,2)*epot_rm(37) + epot_rm(2)*gradt_dmdr(i,37)
                gradt_dmdr(i,54)  = gradt_dmdr(i,1)*epot_rm(41) + epot_rm(1)*gradt_dmdr(i,41)
                gradt_dmdr(i,55)  = gradt_dmdr(i,2)*epot_rm(41) + epot_rm(2)*gradt_dmdr(i,41)
                gradt_dmdr(i,56)  = gradt_dmdr(i,3)*epot_rm(41) + epot_rm(3)*gradt_dmdr(i,41)
                gradt_dmdr(i,57)  = gradt_dmdr(i,1)*epot_rm(42) + epot_rm(1)*gradt_dmdr(i,42)
                gradt_dmdr(i,58)  = gradt_dmdr(i,1)*epot_rm(51) + epot_rm(1)*gradt_dmdr(i,51)
                gradt_dmdr(i,59)  = gradt_dmdr(i,1)*epot_rm(53) + epot_rm(1)*gradt_dmdr(i,53)
                gradt_dmdr(i,60)  = gradt_dmdr(i,1)*epot_rm(55) + epot_rm(1)*gradt_dmdr(i,55)
                gradt_dmdr(i,61)  = gradt_dmdr(i,1)*epot_rm(56) + epot_rm(1)*gradt_dmdr(i,56)
                gradt_dmdr(i,62)  = gradt_dmdr(i,2)*epot_rm(56) + epot_rm(2)*gradt_dmdr(i,56)
                gradt_dmdr(i,63)  = gradt_dmdr(i,1)*epot_rm(62) + epot_rm(1)*gradt_dmdr(i,62)
                gradt_dmdr(i,64)  = gradt_dmdr(i,2)*epot_rm(57) + epot_rm(2)*gradt_dmdr(i,57)
                gradt_dmdr(i,65)  = gradt_dmdr(i,3)*epot_rm(57) + epot_rm(3)*gradt_dmdr(i,57)
                gradt_dmdr(i,66)  = gradt_dmdr(i,4)*epot_rm(57) + epot_rm(4)*gradt_dmdr(i,57)
                gradt_dmdr(i,67)  = gradt_dmdr(i,5)*epot_rm(57) + epot_rm(5)*gradt_dmdr(i,57)
                gradt_dmdr(i,68)  = gradt_dmdr(i,1)*epot_rm(58) + epot_rm(1)*gradt_dmdr(i,58)
                gradt_dmdr(i,69)  = gradt_dmdr(i,3)*epot_rm(58) + epot_rm(3)*gradt_dmdr(i,58)
                gradt_dmdr(i,70)  = gradt_dmdr(i,4)*epot_rm(58) + epot_rm(4)*gradt_dmdr(i,58)
                gradt_dmdr(i,71)  = gradt_dmdr(i,1)*epot_rm(59) + epot_rm(1)*gradt_dmdr(i,59)
                gradt_dmdr(i,72)  = gradt_dmdr(i,2)*epot_rm(59) + epot_rm(2)*gradt_dmdr(i,59)
                gradt_dmdr(i,73)  = gradt_dmdr(i,1)*epot_rm(60) + epot_rm(1)*gradt_dmdr(i,60)
                gradt_dmdr(i,74)  = gradt_dmdr(i,2)*epot_rm(60) + epot_rm(2)*gradt_dmdr(i,60)
                gradt_dmdr(i,75)  = gradt_dmdr(i,1)*epot_rm(61) + epot_rm(1)*gradt_dmdr(i,61)
                gradt_dmdr(i,76)  = gradt_dmdr(i,2)*epot_rm(62) + epot_rm(2)*gradt_dmdr(i,62)
                gradt_dmdr(i,77)  = gradt_dmdr(i,3)*epot_rm(61) + epot_rm(3)*gradt_dmdr(i,61)
                gradt_dmdr(i,78)  = gradt_dmdr(i,3)*epot_rm(62) + epot_rm(3)*gradt_dmdr(i,62)
                gradt_dmdr(i,79)  = gradt_dmdr(i,4)*epot_rm(61) + epot_rm(4)*gradt_dmdr(i,61)
                gradt_dmdr(i,80)  = gradt_dmdr(i,4)*epot_rm(62) + epot_rm(4)*gradt_dmdr(i,62)
                gradt_dmdr(i,81)  = gradt_dmdr(i,5)*epot_rm(59) + epot_rm(5)*gradt_dmdr(i,59)
                gradt_dmdr(i,82)  = gradt_dmdr(i,5)*epot_rm(60) + epot_rm(5)*gradt_dmdr(i,60)
                gradt_dmdr(i,83)  = gradt_dmdr(i,5)*epot_rm(62) + epot_rm(5)*gradt_dmdr(i,62)
                gradt_dmdr(i,84)  = gradt_dmdr(i,6)*epot_rm(58) + epot_rm(6)*gradt_dmdr(i,58)
                gradt_dmdr(i,85)  = gradt_dmdr(i,6)*epot_rm(59) + epot_rm(6)*gradt_dmdr(i,59)
                gradt_dmdr(i,86)  = gradt_dmdr(i,6)*epot_rm(60) + epot_rm(6)*gradt_dmdr(i,60)
                gradt_dmdr(i,87)  = gradt_dmdr(i,6)*epot_rm(61) + epot_rm(6)*gradt_dmdr(i,61)
                gradt_dmdr(i,88)  = gradt_dmdr(i,2)*epot_rm(64) + epot_rm(2)*gradt_dmdr(i,64)
                gradt_dmdr(i,89)  = gradt_dmdr(i,3)*epot_rm(65) + epot_rm(3)*gradt_dmdr(i,65)
                gradt_dmdr(i,90)  = gradt_dmdr(i,4)*epot_rm(66) + epot_rm(4)*gradt_dmdr(i,66)
                gradt_dmdr(i,91)  = gradt_dmdr(i,5)*epot_rm(67) + epot_rm(5)*gradt_dmdr(i,67)
                gradt_dmdr(i,92)  = gradt_dmdr(i,1)*epot_rm(68) + epot_rm(1)*gradt_dmdr(i,68)
                gradt_dmdr(i,93)  = gradt_dmdr(i,3)*epot_rm(69) + epot_rm(3)*gradt_dmdr(i,69)
                gradt_dmdr(i,94)  = gradt_dmdr(i,4)*epot_rm(70) + epot_rm(4)*gradt_dmdr(i,70)
                gradt_dmdr(i,95)  = gradt_dmdr(i,1)*epot_rm(71) + epot_rm(1)*gradt_dmdr(i,71)
                gradt_dmdr(i,96)  = gradt_dmdr(i,2)*epot_rm(72) + epot_rm(2)*gradt_dmdr(i,72)
                gradt_dmdr(i,97)  = gradt_dmdr(i,1)*epot_rm(73) + epot_rm(1)*gradt_dmdr(i,73)
                gradt_dmdr(i,98)  = gradt_dmdr(i,2)*epot_rm(74) + epot_rm(2)*gradt_dmdr(i,74)
                gradt_dmdr(i,99)  = gradt_dmdr(i,1)*epot_rm(75) + epot_rm(1)*gradt_dmdr(i,75)
                gradt_dmdr(i,100) = gradt_dmdr(i,2)*epot_rm(76) + epot_rm(2)*gradt_dmdr(i,76)
                gradt_dmdr(i,101) = gradt_dmdr(i,3)*epot_rm(77) + epot_rm(3)*gradt_dmdr(i,77)
                gradt_dmdr(i,102) = gradt_dmdr(i,3)*epot_rm(78) + epot_rm(3)*gradt_dmdr(i,78)
                gradt_dmdr(i,103) = gradt_dmdr(i,4)*epot_rm(79) + epot_rm(4)*gradt_dmdr(i,79)
                gradt_dmdr(i,104) = gradt_dmdr(i,4)*epot_rm(80) + epot_rm(4)*gradt_dmdr(i,80)
                gradt_dmdr(i,105) = gradt_dmdr(i,5)*epot_rm(81) + epot_rm(5)*gradt_dmdr(i,81)
                gradt_dmdr(i,106) = gradt_dmdr(i,5)*epot_rm(82) + epot_rm(5)*gradt_dmdr(i,82)
                gradt_dmdr(i,107) = gradt_dmdr(i,5)*epot_rm(83) + epot_rm(5)*gradt_dmdr(i,83)
                gradt_dmdr(i,108) = gradt_dmdr(i,6)*epot_rm(84) + epot_rm(6)*gradt_dmdr(i,84)
                gradt_dmdr(i,109) = gradt_dmdr(i,6)*epot_rm(85) + epot_rm(6)*gradt_dmdr(i,85)
                gradt_dmdr(i,110) = gradt_dmdr(i,6)*epot_rm(86) + epot_rm(6)*gradt_dmdr(i,86)
                gradt_dmdr(i,111) = gradt_dmdr(i,6)*epot_rm(87) + epot_rm(6)*gradt_dmdr(i,87)
            enddo
        end subroutine evdmdr
        !> @brief Evaluate EvdPdr for the N4_1A_PIP_2013 potential.
        subroutine evdpdr
            !  The subroutine reads monomials(m) and calculates the
            !  permutation invariant polynomials(p)
            !  For A4 with max. degree 9, the number of polynomials is nob.
            implicit none
            integer :: i
            do i=1, 6
                gradt_dpdr(i,0) = gradt_dmdr(i,0)
                gradt_dpdr(i,1) = gradt_dmdr(i,1) + gradt_dmdr(i,2) + gradt_dmdr(i,3) + gradt_dmdr(i,4) + &
                    & gradt_dmdr(i,5) + gradt_dmdr(i,6)
                gradt_dpdr(i,2) = gradt_dmdr(i,7) + gradt_dmdr(i,8) + gradt_dmdr(i,9)
                gradt_dpdr(i,3) = gradt_dmdr(i,10) + gradt_dmdr(i,11) + gradt_dmdr(i,12) + gradt_dmdr(i,13) + &
                    & gradt_dmdr(i,14) + gradt_dmdr(i,15) + gradt_dmdr(i,16) + gradt_dmdr(i,17) + gradt_dmdr(i,18) + &
                    & gradt_dmdr(i,19) + gradt_dmdr(i,20) + gradt_dmdr(i,21)
                gradt_dpdr(i,4) = gradt_dpdr(i,1)*epot_p(1) + epot_p(1)*gradt_dpdr(i,1) - gradt_dpdr(i,3) - &
                    & gradt_dpdr(i,2) - gradt_dpdr(i,3) - gradt_dpdr(i,2)
                gradt_dpdr(i,5) = gradt_dmdr(i,22) + gradt_dmdr(i,23) + gradt_dmdr(i,24) + gradt_dmdr(i,25) + &
                    & gradt_dmdr(i,26) + gradt_dmdr(i,27) + gradt_dmdr(i,28) + gradt_dmdr(i,29) + gradt_dmdr(i,30) + &
                    & gradt_dmdr(i,31) + gradt_dmdr(i,32) + gradt_dmdr(i,33)
                gradt_dpdr(i,6) = gradt_dmdr(i,34) + gradt_dmdr(i,35) + gradt_dmdr(i,36) + gradt_dmdr(i,37)
                gradt_dpdr(i,7) = gradt_dmdr(i,38) + gradt_dmdr(i,39) + gradt_dmdr(i,40) + gradt_dmdr(i,41)
                gradt_dpdr(i,8) = gradt_dpdr(i,1)*epot_p(2) + epot_p(1)*gradt_dpdr(i,2) - gradt_dpdr(i,5)
                gradt_dpdr(i,9) = gradt_dpdr(i,1)*epot_p(3) + epot_p(1)*gradt_dpdr(i,3) - gradt_dpdr(i,6) - &
                    & gradt_dpdr(i,7) - gradt_dpdr(i,5) - gradt_dpdr(i,6) - gradt_dpdr(i,7) - gradt_dpdr(i,5) - &
                    & gradt_dpdr(i,6) - gradt_dpdr(i,7)
                gradt_dpdr(i,10) = gradt_dpdr(i,1)*epot_p(4) + epot_p(1)*gradt_dpdr(i,4) - gradt_dpdr(i,9) - &
                    & gradt_dpdr(i,8)
                gradt_dpdr(i,11) = gradt_dmdr(i,42) + gradt_dmdr(i,43) + gradt_dmdr(i,44)
                gradt_dpdr(i,12) = gradt_dmdr(i,45) + gradt_dmdr(i,46) + gradt_dmdr(i,47) + gradt_dmdr(i,48) + &
                    & gradt_dmdr(i,49) + gradt_dmdr(i,50) + gradt_dmdr(i,51) + gradt_dmdr(i,52) + gradt_dmdr(i,53) + &
                    & gradt_dmdr(i,54) + gradt_dmdr(i,55) + gradt_dmdr(i,56)
                gradt_dpdr(i,13) = gradt_dpdr(i,2)*epot_p(3) + epot_p(2)*gradt_dpdr(i,3) - gradt_dpdr(i,12)
                gradt_dpdr(i,14) = gradt_dpdr(i,1)*epot_p(5) + epot_p(1)*gradt_dpdr(i,5) - gradt_dpdr(i,12) - &
                    & gradt_dpdr(i,11) - gradt_dpdr(i,13) - gradt_dpdr(i,12) - gradt_dpdr(i,11) - gradt_dpdr(i,11) - &
                    & gradt_dpdr(i,11)
                gradt_dpdr(i,15) = gradt_dpdr(i,1)*epot_p(6) + epot_p(1)*gradt_dpdr(i,6) - gradt_dpdr(i,12)
                gradt_dpdr(i,16) = gradt_dpdr(i,1)*epot_p(7) + epot_p(1)*gradt_dpdr(i,7) - gradt_dpdr(i,12)
                gradt_dpdr(i,17) = gradt_dpdr(i,2)*epot_p(2) + epot_p(2)*gradt_dpdr(i,2) - gradt_dpdr(i,11) - &
                    & gradt_dpdr(i,11)
                gradt_dpdr(i,18) = gradt_dpdr(i,3)*epot_p(3) + epot_p(3)*gradt_dpdr(i,3) - gradt_dpdr(i,12) - &
                    & gradt_dpdr(i,11) - gradt_dpdr(i,15) - gradt_dpdr(i,16) - gradt_dpdr(i,14) - gradt_dpdr(i,12) - &
                    & gradt_dpdr(i,11) - gradt_dpdr(i,15) - gradt_dpdr(i,16) - gradt_dpdr(i,14) - gradt_dpdr(i,12) - &
                    & gradt_dpdr(i,11) - gradt_dpdr(i,12) - gradt_dpdr(i,11)
                gradt_dpdr(i,19) = gradt_dpdr(i,2)*epot_p(4) + epot_p(2)*gradt_dpdr(i,4) - gradt_dpdr(i,14)
                gradt_dpdr(i,20) = gradt_dpdr(i,3)*epot_p(4) + epot_p(3)*gradt_dpdr(i,4) - gradt_dpdr(i,15) - &
                    & gradt_dpdr(i,16) - gradt_dpdr(i,13)
                gradt_dpdr(i,21) = gradt_dpdr(i,1)*epot_p(10) + epot_p(1)*gradt_dpdr(i,10) - gradt_dpdr(i,20) - &
                    & gradt_dpdr(i,19)
                gradt_dpdr(i,22) = gradt_dmdr(i,57) + gradt_dmdr(i,58) + gradt_dmdr(i,59) + gradt_dmdr(i,60) + &
                    & gradt_dmdr(i,61) + gradt_dmdr(i,62)
                gradt_dpdr(i,23) = gradt_dpdr(i,1)*epot_p(11) + epot_p(1)*gradt_dpdr(i,11) - gradt_dpdr(i,22)
                gradt_dpdr(i,24) = gradt_dpdr(i,2)*epot_p(6) + epot_p(2)*gradt_dpdr(i,6)
                gradt_dpdr(i,25) = gradt_dpdr(i,2)*epot_p(7) + epot_p(2)*gradt_dpdr(i,7)
                gradt_dpdr(i,26) = gradt_dpdr(i,1)*epot_p(12) + epot_p(1)*gradt_dpdr(i,12) - gradt_dpdr(i,22) - &
                    & gradt_dpdr(i,24) - gradt_dpdr(i,25) - gradt_dpdr(i,22) - gradt_dpdr(i,22) - gradt_dpdr(i,22)
                gradt_dpdr(i,27) = gradt_dpdr(i,2)*epot_p(5) + epot_p(2)*gradt_dpdr(i,5) - gradt_dpdr(i,22) - &
                    & gradt_dpdr(i,23) - gradt_dpdr(i,22)
                gradt_dpdr(i,28) = gradt_dpdr(i,3)*epot_p(5) + epot_p(3)*gradt_dpdr(i,5) - gradt_dpdr(i,22) - &
                    & gradt_dpdr(i,26) - gradt_dpdr(i,24) - gradt_dpdr(i,25) - gradt_dpdr(i,23) - gradt_dpdr(i,22) - &
                    & gradt_dpdr(i,24) - gradt_dpdr(i,25) - gradt_dpdr(i,23) - gradt_dpdr(i,22) - gradt_dpdr(i,22)
                gradt_dpdr(i,29) = gradt_dpdr(i,3)*epot_p(6) + epot_p(3)*gradt_dpdr(i,6) - gradt_dpdr(i,22) - &
                    & gradt_dpdr(i,26) - gradt_dpdr(i,22)
                gradt_dpdr(i,30) = gradt_dpdr(i,3)*epot_p(7) + epot_p(3)*gradt_dpdr(i,7) - gradt_dpdr(i,22) - &
                    & gradt_dpdr(i,26) - gradt_dpdr(i,22)
                gradt_dpdr(i,31) = gradt_dpdr(i,2)*epot_p(9) + epot_p(2)*gradt_dpdr(i,9) - gradt_dpdr(i,26) - &
                    & gradt_dpdr(i,28)
                gradt_dpdr(i,32) = gradt_dpdr(i,1)*epot_p(14) + epot_p(1)*gradt_dpdr(i,14) - gradt_dpdr(i,26) - &
                    & gradt_dpdr(i,23) - gradt_dpdr(i,28)
                gradt_dpdr(i,33) = gradt_dpdr(i,4)*epot_p(6) + epot_p(4)*gradt_dpdr(i,6) - gradt_dpdr(i,25)
                gradt_dpdr(i,34) = gradt_dpdr(i,4)*epot_p(7) + epot_p(4)*gradt_dpdr(i,7) - gradt_dpdr(i,24)
                gradt_dpdr(i,35) = gradt_dpdr(i,1)*epot_p(17) + epot_p(1)*gradt_dpdr(i,17) - gradt_dpdr(i,27)
                gradt_dpdr(i,36) = gradt_dpdr(i,1)*epot_p(18) + epot_p(1)*gradt_dpdr(i,18) - gradt_dpdr(i,29) - &
                    & gradt_dpdr(i,30) - gradt_dpdr(i,28)
                gradt_dpdr(i,37) = gradt_dpdr(i,2)*epot_p(10) + epot_p(2)*gradt_dpdr(i,10) - gradt_dpdr(i,32)
                gradt_dpdr(i,38) = gradt_dpdr(i,3)*epot_p(10) + epot_p(3)*gradt_dpdr(i,10) - gradt_dpdr(i,33) - &
                    & gradt_dpdr(i,34) - gradt_dpdr(i,31)
                gradt_dpdr(i,39) = gradt_dpdr(i,1)*epot_p(21) + epot_p(1)*gradt_dpdr(i,21) - gradt_dpdr(i,38) - &
                    & gradt_dpdr(i,37)
                gradt_dpdr(i,40) = gradt_dmdr(i,63)
                gradt_dpdr(i,41) = gradt_dmdr(i,64) + gradt_dmdr(i,65) + gradt_dmdr(i,66) + gradt_dmdr(i,67) + &
                    & gradt_dmdr(i,68) + gradt_dmdr(i,69) + gradt_dmdr(i,70) + gradt_dmdr(i,71) + gradt_dmdr(i,72) + &
                    & gradt_dmdr(i,73) + gradt_dmdr(i,74) + gradt_dmdr(i,75) + gradt_dmdr(i,76) + gradt_dmdr(i,77) + &
                    & gradt_dmdr(i,78) + gradt_dmdr(i,79) + gradt_dmdr(i,80) + gradt_dmdr(i,81) + gradt_dmdr(i,82) + &
                    & gradt_dmdr(i,83) + gradt_dmdr(i,84) + gradt_dmdr(i,85) + gradt_dmdr(i,86) + gradt_dmdr(i,87)
                gradt_dpdr(i,42) = gradt_dpdr(i,1)*epot_p(22) + epot_p(1)*gradt_dpdr(i,22) - gradt_dpdr(i,40) - &
                    & gradt_dpdr(i,41) - gradt_dpdr(i,40) - gradt_dpdr(i,40) - gradt_dpdr(i,40) - gradt_dpdr(i,40) - &
                    & gradt_dpdr(i,40)
                gradt_dpdr(i,43) = gradt_dpdr(i,2)*epot_p(11) + epot_p(2)*gradt_dpdr(i,11) - gradt_dpdr(i,40) - &
                    & gradt_dpdr(i,40) - gradt_dpdr(i,40)
                gradt_dpdr(i,44) = gradt_dpdr(i,2)*epot_p(12) + epot_p(2)*gradt_dpdr(i,12) - gradt_dpdr(i,41)
                gradt_dpdr(i,45) = gradt_dpdr(i,3)*epot_p(11) + epot_p(3)*gradt_dpdr(i,11) - gradt_dpdr(i,41)
                gradt_dpdr(i,46) = gradt_dpdr(i,5)*epot_p(6) + epot_p(5)*gradt_dpdr(i,6) - gradt_dpdr(i,41)
                gradt_dpdr(i,47) = gradt_dpdr(i,5)*epot_p(7) + epot_p(5)*gradt_dpdr(i,7) - gradt_dpdr(i,41)
                gradt_dpdr(i,48) = gradt_dpdr(i,6)*epot_p(7) + epot_p(6)*gradt_dpdr(i,7) - gradt_dpdr(i,40) - &
                    & gradt_dpdr(i,40) - gradt_dpdr(i,40) - gradt_dpdr(i,40)
                gradt_dpdr(i,49) = gradt_dpdr(i,4)*epot_p(11) + epot_p(4)*gradt_dpdr(i,11) - gradt_dpdr(i,42)
                gradt_dpdr(i,50) = gradt_dpdr(i,2)*epot_p(15) + epot_p(2)*gradt_dpdr(i,15) - gradt_dpdr(i,46)
                gradt_dpdr(i,51) = gradt_dpdr(i,2)*epot_p(16) + epot_p(2)*gradt_dpdr(i,16) - gradt_dpdr(i,47)
                gradt_dpdr(i,52) = gradt_dpdr(i,4)*epot_p(12) + epot_p(4)*gradt_dpdr(i,12) - gradt_dpdr(i,41) - &
                    & gradt_dpdr(i,50) - gradt_dpdr(i,51)
                gradt_dpdr(i,53) = gradt_dpdr(i,2)*epot_p(14) + epot_p(2)*gradt_dpdr(i,14) - gradt_dpdr(i,42) - &
                    & gradt_dpdr(i,49) - gradt_dpdr(i,42)
                gradt_dpdr(i,54) = gradt_dpdr(i,6)*epot_p(6) + epot_p(6)*gradt_dpdr(i,6) - gradt_dpdr(i,42) - &
                    & gradt_dpdr(i,42)
                gradt_dpdr(i,55) = gradt_dpdr(i,7)*epot_p(7) + epot_p(7)*gradt_dpdr(i,7) - gradt_dpdr(i,42) - &
                    & gradt_dpdr(i,42)
                gradt_dpdr(i,56) = gradt_dpdr(i,3)*epot_p(17) + epot_p(3)*gradt_dpdr(i,17) - gradt_dpdr(i,44)
                gradt_dpdr(i,57) = gradt_dpdr(i,2)*epot_p(18) + epot_p(2)*gradt_dpdr(i,18) - gradt_dpdr(i,48)
                gradt_dpdr(i,58) = gradt_dpdr(i,3)*epot_p(14) + epot_p(3)*gradt_dpdr(i,14) - gradt_dpdr(i,41) - &
                    & gradt_dpdr(i,52) - gradt_dpdr(i,46) - gradt_dpdr(i,47) - gradt_dpdr(i,45) - gradt_dpdr(i,45)
                gradt_dpdr(i,59) = gradt_dpdr(i,6)*epot_p(9) + epot_p(6)*gradt_dpdr(i,9) - gradt_dpdr(i,41) - &
                    & gradt_dpdr(i,52) - gradt_dpdr(i,47)
                gradt_dpdr(i,60) = gradt_dpdr(i,7)*epot_p(9) + epot_p(7)*gradt_dpdr(i,9) - gradt_dpdr(i,41) - &
                    & gradt_dpdr(i,52) - gradt_dpdr(i,46)
                gradt_dpdr(i,61) = gradt_dpdr(i,2)*epot_p(20) + epot_p(2)*gradt_dpdr(i,20) - gradt_dpdr(i,52) - &
                    & gradt_dpdr(i,58)
                gradt_dpdr(i,62) = gradt_dpdr(i,1)*epot_p(32) + epot_p(1)*gradt_dpdr(i,32) - gradt_dpdr(i,52) - &
                    & gradt_dpdr(i,49) - gradt_dpdr(i,58)
                gradt_dpdr(i,63) = gradt_dpdr(i,6)*epot_p(10) + epot_p(6)*gradt_dpdr(i,10) - gradt_dpdr(i,51)
                gradt_dpdr(i,64) = gradt_dpdr(i,7)*epot_p(10) + epot_p(7)*gradt_dpdr(i,10) - gradt_dpdr(i,50)
                gradt_dpdr(i,65) = gradt_dpdr(i,2)*epot_p(17) + epot_p(2)*gradt_dpdr(i,17) - gradt_dpdr(i,43)
                gradt_dpdr(i,66) = gradt_dpdr(i,3)*epot_p(18) + epot_p(3)*gradt_dpdr(i,18) - gradt_dpdr(i,46) - &
                    & gradt_dpdr(i,47) - gradt_dpdr(i,45) - gradt_dpdr(i,59) - gradt_dpdr(i,60) - gradt_dpdr(i,58)
                gradt_dpdr(i,67) = gradt_dpdr(i,2)*epot_p(19) + epot_p(2)*gradt_dpdr(i,19) - gradt_dpdr(i,49)
                gradt_dpdr(i,68) = gradt_dpdr(i,1)*epot_p(36) + epot_p(1)*gradt_dpdr(i,36) - gradt_dpdr(i,59) - &
                    & gradt_dpdr(i,60) - gradt_dpdr(i,58) - gradt_dpdr(i,57) - gradt_dpdr(i,66) - gradt_dpdr(i,66)
                gradt_dpdr(i,69) = gradt_dpdr(i,2)*epot_p(21) + epot_p(2)*gradt_dpdr(i,21) - gradt_dpdr(i,62)
                gradt_dpdr(i,70) = gradt_dpdr(i,3)*epot_p(21) + epot_p(3)*gradt_dpdr(i,21) - gradt_dpdr(i,63) - &
                    & gradt_dpdr(i,64) - gradt_dpdr(i,61)
                gradt_dpdr(i,71) = gradt_dpdr(i,1)*epot_p(39) + epot_p(1)*gradt_dpdr(i,39) - gradt_dpdr(i,70) - &
                    & gradt_dpdr(i,69)
                gradt_dpdr(i,72) = gradt_dpdr(i,40)*epot_p(1) + epot_p(40)*gradt_dpdr(i,1)
                gradt_dpdr(i,73) = gradt_dpdr(i,2)*epot_p(22) + epot_p(2)*gradt_dpdr(i,22) - gradt_dpdr(i,72)
                gradt_dpdr(i,74) = gradt_dpdr(i,6)*epot_p(11) + epot_p(6)*gradt_dpdr(i,11)
                gradt_dpdr(i,75) = gradt_dpdr(i,7)*epot_p(11) + epot_p(7)*gradt_dpdr(i,11)
                gradt_dpdr(i,76) = gradt_dpdr(i,3)*epot_p(22) + epot_p(3)*gradt_dpdr(i,22) - gradt_dpdr(i,72) - &
                    & gradt_dpdr(i,74) - gradt_dpdr(i,75) - gradt_dpdr(i,72) - gradt_dpdr(i,72) - gradt_dpdr(i,72)
                gradt_dpdr(i,77) = gradt_dmdr(i,88) + gradt_dmdr(i,89) + gradt_dmdr(i,90) + gradt_dmdr(i,91) + &
                    & gradt_dmdr(i,92) + gradt_dmdr(i,93) + gradt_dmdr(i,94) + gradt_dmdr(i,95) + gradt_dmdr(i,96) + &
                    & gradt_dmdr(i,97) + gradt_dmdr(i,98) + gradt_dmdr(i,99) + gradt_dmdr(i,100) + gradt_dmdr(i, &
                    & 101) + gradt_dmdr(i,102) + gradt_dmdr(i,103) + gradt_dmdr(i,104) + gradt_dmdr(i,105) + &
                    & gradt_dmdr(i,106) + gradt_dmdr(i,107) + gradt_dmdr(i,108) + gradt_dmdr(i,109) + gradt_dmdr(i, &
                    & 110) + gradt_dmdr(i,111)
                gradt_dpdr(i,78)  = gradt_dpdr(i,1)*epot_p(42) + epot_p(1)*gradt_dpdr(i,42) - gradt_dpdr(i,72) - &
                    & gradt_dpdr(i,76)
                gradt_dpdr(i,79)  = gradt_dpdr(i,5)*epot_p(11) + epot_p(5)*gradt_dpdr(i,11) - gradt_dpdr(i,72) - &
                    & gradt_dpdr(i,73) - gradt_dpdr(i,72)
                gradt_dpdr(i,80)  = gradt_dpdr(i,2)*epot_p(26) + epot_p(2)*gradt_dpdr(i,26) - gradt_dpdr(i,76) - &
                    & gradt_dpdr(i,77)
                gradt_dpdr(i,81)  = gradt_dpdr(i,6)*epot_p(12) + epot_p(6)*gradt_dpdr(i,12) - gradt_dpdr(i,72) - &
                    & gradt_dpdr(i,76) - gradt_dpdr(i,72)
                gradt_dpdr(i,82)  = gradt_dpdr(i,7)*epot_p(12) + epot_p(7)*gradt_dpdr(i,12) - gradt_dpdr(i,72) - &
                    & gradt_dpdr(i,76) - gradt_dpdr(i,72)
                gradt_dpdr(i,83) = gradt_dpdr(i,8)*epot_p(11) + epot_p(8)*gradt_dpdr(i,11) - gradt_dpdr(i,72)
                gradt_dpdr(i,84) = gradt_dpdr(i,6)*epot_p(17) + epot_p(6)*gradt_dpdr(i,17)
                gradt_dpdr(i,85) = gradt_dpdr(i,7)*epot_p(17) + epot_p(7)*gradt_dpdr(i,17)
                gradt_dpdr(i,86)  = gradt_dpdr(i,9)*epot_p(11) + epot_p(9)*gradt_dpdr(i,11) - gradt_dpdr(i,76) - &
                    & gradt_dpdr(i,77)
                gradt_dpdr(i,87) = gradt_dpdr(i,2)*epot_p(29) + epot_p(2)*gradt_dpdr(i,29) - gradt_dpdr(i,81)
                gradt_dpdr(i,88)  = gradt_dpdr(i,6)*epot_p(14) + epot_p(6)*gradt_dpdr(i,14) - gradt_dpdr(i,75) - &
                    & gradt_dpdr(i,75)
                gradt_dpdr(i,89) = gradt_dpdr(i,2)*epot_p(30) + epot_p(2)*gradt_dpdr(i,30) - gradt_dpdr(i,82)
                gradt_dpdr(i,90)  = gradt_dpdr(i,7)*epot_p(14) + epot_p(7)*gradt_dpdr(i,14) - gradt_dpdr(i,74) - &
                    & gradt_dpdr(i,74)
                gradt_dpdr(i,91)  = gradt_dpdr(i,1)*epot_p(48) + epot_p(1)*gradt_dpdr(i,48) - gradt_dpdr(i,76) - &
                    & gradt_dpdr(i,81) - gradt_dpdr(i,82)
                gradt_dpdr(i,92) = gradt_dpdr(i,10)*epot_p(11) + epot_p(10)*gradt_dpdr(i,11) - gradt_dpdr(i,78)
                gradt_dpdr(i,93) = gradt_dpdr(i,2)*epot_p(33) + epot_p(2)*gradt_dpdr(i,33) - gradt_dpdr(i,88)
                gradt_dpdr(i,94) = gradt_dpdr(i,2)*epot_p(34) + epot_p(2)*gradt_dpdr(i,34) - gradt_dpdr(i,90)
                gradt_dpdr(i,95)  = gradt_dpdr(i,10)*epot_p(12) + epot_p(10)*gradt_dpdr(i,12) - gradt_dpdr(i,77) - &
                    & gradt_dpdr(i,93) - gradt_dpdr(i,94)
                gradt_dpdr(i,96)  = gradt_dpdr(i,2)*epot_p(27) + epot_p(2)*gradt_dpdr(i,27) - gradt_dpdr(i,73) - &
                    & gradt_dpdr(i,79)
                gradt_dpdr(i,97)  = gradt_dpdr(i,2)*epot_p(28) + epot_p(2)*gradt_dpdr(i,28) - gradt_dpdr(i,76) - &
                    & gradt_dpdr(i,86)
                gradt_dpdr(i,98)  = gradt_dpdr(i,1)*epot_p(53) + epot_p(1)*gradt_dpdr(i,53) - gradt_dpdr(i,80) - &
                    & gradt_dpdr(i,79) - gradt_dpdr(i,97)
                gradt_dpdr(i,99)  = gradt_dpdr(i,1)*epot_p(54) + epot_p(1)*gradt_dpdr(i,54) - gradt_dpdr(i,81)
                gradt_dpdr(i,100) = gradt_dpdr(i,1)*epot_p(55) + epot_p(1)*gradt_dpdr(i,55) - gradt_dpdr(i,82)
                gradt_dpdr(i,101) = gradt_dpdr(i,5)*epot_p(18) + epot_p(5)*gradt_dpdr(i,18) - gradt_dpdr(i,76) - &
                    & gradt_dpdr(i,91) - gradt_dpdr(i,87) - gradt_dpdr(i,89) - gradt_dpdr(i,86)
                gradt_dpdr(i,102) = gradt_dpdr(i,6)*epot_p(18) + epot_p(6)*gradt_dpdr(i,18) - gradt_dpdr(i,74) - &
                    & gradt_dpdr(i,90)
                gradt_dpdr(i,103) = gradt_dpdr(i,7)*epot_p(18) + epot_p(7)*gradt_dpdr(i,18) - gradt_dpdr(i,75) - &
                    & gradt_dpdr(i,88)
                gradt_dpdr(i,104) = gradt_dpdr(i,2)*epot_p(31) + epot_p(2)*gradt_dpdr(i,31) - gradt_dpdr(i,77) - &
                    & gradt_dpdr(i,86)
                gradt_dpdr(i,105) = gradt_dpdr(i,2)*epot_p(36) + epot_p(2)*gradt_dpdr(i,36) - gradt_dpdr(i,91) - &
                    & gradt_dpdr(i,101)
                gradt_dpdr(i,106) = gradt_dpdr(i,3)*epot_p(32) + epot_p(3)*gradt_dpdr(i,32) - gradt_dpdr(i,77) - &
                    & gradt_dpdr(i,95) - gradt_dpdr(i,88) - gradt_dpdr(i,90) - gradt_dpdr(i,86)
                gradt_dpdr(i,107) = gradt_dpdr(i,4)*epot_p(29) + epot_p(4)*gradt_dpdr(i,29) - gradt_dpdr(i,82) - &
                    & gradt_dpdr(i,80) - gradt_dpdr(i,99)
                gradt_dpdr(i,108) = gradt_dpdr(i,4)*epot_p(30) + epot_p(4)*gradt_dpdr(i,30) - gradt_dpdr(i,81) - &
                    & gradt_dpdr(i,80) - gradt_dpdr(i,100)
                gradt_dpdr(i,109) = gradt_dpdr(i,2)*epot_p(38) + epot_p(2)*gradt_dpdr(i,38) - gradt_dpdr(i,95) - &
                    & gradt_dpdr(i,106)
                gradt_dpdr(i,110) = gradt_dpdr(i,1)*epot_p(62) + epot_p(1)*gradt_dpdr(i,62) - gradt_dpdr(i,95) - &
                    & gradt_dpdr(i,92) - gradt_dpdr(i,106)
                gradt_dpdr(i,111) = gradt_dpdr(i,6)*epot_p(21) + epot_p(6)*gradt_dpdr(i,21) - gradt_dpdr(i,94)
                gradt_dpdr(i,112) = gradt_dpdr(i,7)*epot_p(21) + epot_p(7)*gradt_dpdr(i,21) - gradt_dpdr(i,93)
                gradt_dpdr(i,113) = gradt_dpdr(i,1)*epot_p(65) + epot_p(1)*gradt_dpdr(i,65) - gradt_dpdr(i,96)
                gradt_dpdr(i,114) = gradt_dpdr(i,1)*epot_p(66) + epot_p(1)*gradt_dpdr(i,66) - gradt_dpdr(i,102) - &
                    & gradt_dpdr(i,103) - gradt_dpdr(i,101)
                gradt_dpdr(i,115) = gradt_dpdr(i,2)*epot_p(37) + epot_p(2)*gradt_dpdr(i,37) - gradt_dpdr(i,92)
                gradt_dpdr(i,116) = gradt_dpdr(i,10)*epot_p(18) + epot_p(10)*gradt_dpdr(i,18) - gradt_dpdr(i,99) - &
                    & gradt_dpdr(i,100) - gradt_dpdr(i,97)
                gradt_dpdr(i,117) = gradt_dpdr(i,2)*epot_p(39) + epot_p(2)*gradt_dpdr(i,39) - gradt_dpdr(i,110)
                gradt_dpdr(i,118) = gradt_dpdr(i,3)*epot_p(39) + epot_p(3)*gradt_dpdr(i,39) - gradt_dpdr(i,111) - &
                    & gradt_dpdr(i,112) - gradt_dpdr(i,109)
                gradt_dpdr(i,119) = gradt_dpdr(i,1)*epot_p(71) + epot_p(1)*gradt_dpdr(i,71) - gradt_dpdr(i,118) - &
                    & gradt_dpdr(i,117)
                gradt_dpdr(i,120) = gradt_dpdr(i,40)*epot_p(2) + epot_p(40)*gradt_dpdr(i,2)
                gradt_dpdr(i,121) = gradt_dpdr(i,40)*epot_p(3) + epot_p(40)*gradt_dpdr(i,3)
                gradt_dpdr(i,122) = gradt_dpdr(i,40)*epot_p(4) + epot_p(40)*gradt_dpdr(i,4)
                gradt_dpdr(i,123) = gradt_dpdr(i,11)*epot_p(12) + epot_p(11)*gradt_dpdr(i,12) - gradt_dpdr(i,121)
                gradt_dpdr(i,124) = gradt_dpdr(i,2)*epot_p(42) + epot_p(2)*gradt_dpdr(i,42) - gradt_dpdr(i,122)
                gradt_dpdr(i,125) = gradt_dpdr(i,6)*epot_p(22) + epot_p(6)*gradt_dpdr(i,22) - gradt_dpdr(i,121)
                gradt_dpdr(i,126) = gradt_dpdr(i,7)*epot_p(22) + epot_p(7)*gradt_dpdr(i,22) - gradt_dpdr(i,121)
                gradt_dpdr(i,127) = gradt_dpdr(i,2)*epot_p(41) + epot_p(2)*gradt_dpdr(i,41) - gradt_dpdr(i,121) - &
                    & gradt_dpdr(i,123) - gradt_dpdr(i,121)
                gradt_dpdr(i,128) = gradt_dpdr(i,6)*epot_p(23) + epot_p(6)*gradt_dpdr(i,23) - gradt_dpdr(i,123)
                gradt_dpdr(i,129) = gradt_dpdr(i,7)*epot_p(23) + epot_p(7)*gradt_dpdr(i,23) - gradt_dpdr(i,123)
                gradt_dpdr(i,130) = gradt_dpdr(i,3)*epot_p(41) + epot_p(3)*gradt_dpdr(i,41) - gradt_dpdr(i,122) - &
                    & gradt_dpdr(i,121) - gradt_dpdr(i,120) - gradt_dpdr(i,125) - gradt_dpdr(i,128) - gradt_dpdr(i, &
                    & 126) - gradt_dpdr(i,129) - gradt_dpdr(i,124) - gradt_dpdr(i,123) - gradt_dpdr(i,122) - &
                    & gradt_dpdr(i,121) - gradt_dpdr(i,120) - gradt_dpdr(i,125) - gradt_dpdr(i,126) - gradt_dpdr(i, &
                    & 124) - gradt_dpdr(i,123) - gradt_dpdr(i,122) - gradt_dpdr(i,121) - gradt_dpdr(i,120) - &
                    & gradt_dpdr(i,122) - gradt_dpdr(i,121) - gradt_dpdr(i,120) - gradt_dpdr(i,120) - gradt_dpdr(i, &
                    & 120) - gradt_dpdr(i,120) - gradt_dpdr(i,120)
                gradt_dpdr(i,131) = gradt_dpdr(i,3)*epot_p(42) + epot_p(3)*gradt_dpdr(i,42) - gradt_dpdr(i,121) - &
                    & gradt_dpdr(i,125) - gradt_dpdr(i,126) - gradt_dpdr(i,121)
                gradt_dpdr(i,132) = gradt_dpdr(i,4)*epot_p(41) + epot_p(4)*gradt_dpdr(i,41) - gradt_dpdr(i,121) - &
                    & gradt_dpdr(i,131) - gradt_dpdr(i,128) - gradt_dpdr(i,129) - gradt_dpdr(i,127) - gradt_dpdr(i,121)
                gradt_dpdr(i,133) = gradt_dpdr(i,1)*epot_p(78) + epot_p(1)*gradt_dpdr(i,78) - gradt_dpdr(i,122) - &
                    & gradt_dpdr(i,131)
                gradt_dpdr(i,134) = gradt_dpdr(i,11)*epot_p(11) + epot_p(11)*gradt_dpdr(i,11) - gradt_dpdr(i,120) - &
                    & gradt_dpdr(i,120)
                gradt_dpdr(i,135) = gradt_dpdr(i,2)*epot_p(48) + epot_p(2)*gradt_dpdr(i,48) - gradt_dpdr(i,130)
                gradt_dpdr(i,136) = gradt_dpdr(i,11)*epot_p(17) + epot_p(11)*gradt_dpdr(i,17) - gradt_dpdr(i,120)
                gradt_dpdr(i,137) = gradt_dpdr(i,2)*epot_p(44) + epot_p(2)*gradt_dpdr(i,44) - gradt_dpdr(i,123)
                gradt_dpdr(i,138) = gradt_dpdr(i,2)*epot_p(45) + epot_p(2)*gradt_dpdr(i,45) - gradt_dpdr(i,121)
                gradt_dpdr(i,139) = gradt_dpdr(i,11)*epot_p(14) + epot_p(11)*gradt_dpdr(i,14) - gradt_dpdr(i,122) - &
                    & gradt_dpdr(i,124) - gradt_dpdr(i,122)
                gradt_dpdr(i,140) = gradt_dpdr(i,6)*epot_p(27) + epot_p(6)*gradt_dpdr(i,27) - gradt_dpdr(i,127)
                gradt_dpdr(i,141) = gradt_dpdr(i,2)*epot_p(54) + epot_p(2)*gradt_dpdr(i,54)
                gradt_dpdr(i,142) = gradt_dpdr(i,7)*epot_p(27) + epot_p(7)*gradt_dpdr(i,27) - gradt_dpdr(i,127)
                gradt_dpdr(i,143) = gradt_dpdr(i,2)*epot_p(52) + epot_p(2)*gradt_dpdr(i,52) - gradt_dpdr(i,131) - &
                    & gradt_dpdr(i,132)
                gradt_dpdr(i,144) = gradt_dpdr(i,1)*epot_p(81) + epot_p(1)*gradt_dpdr(i,81) - gradt_dpdr(i,125) - &
                    & gradt_dpdr(i,141) - gradt_dpdr(i,135) - gradt_dpdr(i,125)
                gradt_dpdr(i,145) = gradt_dpdr(i,2)*epot_p(55) + epot_p(2)*gradt_dpdr(i,55)
                gradt_dpdr(i,146) = gradt_dpdr(i,1)*epot_p(82) + epot_p(1)*gradt_dpdr(i,82) - gradt_dpdr(i,126) - &
                    & gradt_dpdr(i,135) - gradt_dpdr(i,145) - gradt_dpdr(i,126)
                gradt_dpdr(i,147) = gradt_dpdr(i,11)*epot_p(18) + epot_p(11)*gradt_dpdr(i,18) - gradt_dpdr(i,130)
                gradt_dpdr(i,148) = gradt_dpdr(i,6)*epot_p(28) + epot_p(6)*gradt_dpdr(i,28) - gradt_dpdr(i,129) - &
                    & gradt_dpdr(i,123) - gradt_dpdr(i,143)
                gradt_dpdr(i,149) = gradt_dpdr(i,7)*epot_p(28) + epot_p(7)*gradt_dpdr(i,28) - gradt_dpdr(i,128) - &
                    & gradt_dpdr(i,123) - gradt_dpdr(i,143)
                gradt_dpdr(i,150) = gradt_dpdr(i,6)*epot_p(30) + epot_p(6)*gradt_dpdr(i,30) - gradt_dpdr(i,121) - &
                    & gradt_dpdr(i,146)
                gradt_dpdr(i,151) = gradt_dpdr(i,11)*epot_p(19) + epot_p(11)*gradt_dpdr(i,19) - gradt_dpdr(i,122)
                gradt_dpdr(i,152) = gradt_dpdr(i,2)*epot_p(50) + epot_p(2)*gradt_dpdr(i,50) - gradt_dpdr(i,128)
                gradt_dpdr(i,153) = gradt_dpdr(i,2)*epot_p(51) + epot_p(2)*gradt_dpdr(i,51) - gradt_dpdr(i,129)
                gradt_dpdr(i,154) = gradt_dpdr(i,11)*epot_p(20) + epot_p(11)*gradt_dpdr(i,20) - gradt_dpdr(i,131) - &
                    & gradt_dpdr(i,132)
                gradt_dpdr(i,155) = gradt_dpdr(i,2)*epot_p(59) + epot_p(2)*gradt_dpdr(i,59) - gradt_dpdr(i,144) - &
                    & gradt_dpdr(i,148)
                gradt_dpdr(i,156) = gradt_dpdr(i,6)*epot_p(32) + epot_p(6)*gradt_dpdr(i,32) - gradt_dpdr(i,129)
                gradt_dpdr(i,157) = gradt_dpdr(i,2)*epot_p(60) + epot_p(2)*gradt_dpdr(i,60) - gradt_dpdr(i,146) - &
                    & gradt_dpdr(i,149)
                gradt_dpdr(i,158) = gradt_dpdr(i,7)*epot_p(32) + epot_p(7)*gradt_dpdr(i,32) - gradt_dpdr(i,128)
                gradt_dpdr(i,159) = gradt_dpdr(i,6)*epot_p(34) + epot_p(6)*gradt_dpdr(i,34) - gradt_dpdr(i,122) - &
                    & gradt_dpdr(i,145) - gradt_dpdr(i,122)
                gradt_dpdr(i,160) = gradt_dpdr(i,11)*epot_p(21) + epot_p(11)*gradt_dpdr(i,21) - gradt_dpdr(i,133)
                gradt_dpdr(i,161) = gradt_dpdr(i,2)*epot_p(63) + epot_p(2)*gradt_dpdr(i,63) - gradt_dpdr(i,156)
                gradt_dpdr(i,162) = gradt_dpdr(i,2)*epot_p(64) + epot_p(2)*gradt_dpdr(i,64) - gradt_dpdr(i,158)
                gradt_dpdr(i,163) = gradt_dpdr(i,12)*epot_p(21) + epot_p(12)*gradt_dpdr(i,21) - gradt_dpdr(i,132) - &
                    & gradt_dpdr(i,161) - gradt_dpdr(i,162)
                gradt_dpdr(i,164) = gradt_dpdr(i,2)*epot_p(53) + epot_p(2)*gradt_dpdr(i,53) - gradt_dpdr(i,124) - &
                    & gradt_dpdr(i,139)
                gradt_dpdr(i,165) = gradt_dpdr(i,2)*epot_p(58) + epot_p(2)*gradt_dpdr(i,58) - gradt_dpdr(i,131) - &
                    & gradt_dpdr(i,154)
                gradt_dpdr(i,166) = gradt_dpdr(i,3)*epot_p(54) + epot_p(3)*gradt_dpdr(i,54) - gradt_dpdr(i,125) - &
                    & gradt_dpdr(i,144)
                gradt_dpdr(i,167) = gradt_dpdr(i,3)*epot_p(55) + epot_p(3)*gradt_dpdr(i,55) - gradt_dpdr(i,126) - &
                    & gradt_dpdr(i,146)
                gradt_dpdr(i,168) = gradt_dpdr(i,3)*epot_p(65) + epot_p(3)*gradt_dpdr(i,65) - gradt_dpdr(i,137)
                gradt_dpdr(i,169) = gradt_dpdr(i,17)*epot_p(18) + epot_p(17)*gradt_dpdr(i,18) - gradt_dpdr(i,135)
                gradt_dpdr(i,170) = gradt_dpdr(i,1)*epot_p(98) + epot_p(1)*gradt_dpdr(i,98) - gradt_dpdr(i,143) - &
                    & gradt_dpdr(i,139) - gradt_dpdr(i,165)
                gradt_dpdr(i,171) = gradt_dpdr(i,4)*epot_p(54) + epot_p(4)*gradt_dpdr(i,54) - gradt_dpdr(i,135)
                gradt_dpdr(i,172) = gradt_dpdr(i,4)*epot_p(55) + epot_p(4)*gradt_dpdr(i,55) - gradt_dpdr(i,135)
                gradt_dpdr(i,173) = gradt_dpdr(i,2)*epot_p(66) + epot_p(2)*gradt_dpdr(i,66) - gradt_dpdr(i,150)
                gradt_dpdr(i,174) = gradt_dpdr(i,1)*epot_p(101) + epot_p(1)*gradt_dpdr(i,101) - gradt_dpdr(i,148) - &
                    & gradt_dpdr(i,149) - gradt_dpdr(i,147) - gradt_dpdr(i,173) - gradt_dpdr(i,165) - gradt_dpdr(i,147)
                gradt_dpdr(i,175) = gradt_dpdr(i,1)*epot_p(102) + epot_p(1)*gradt_dpdr(i,102) - gradt_dpdr(i,150) - &
                    & gradt_dpdr(i,148) - gradt_dpdr(i,166)
                gradt_dpdr(i,176) = gradt_dpdr(i,1)*epot_p(103) + epot_p(1)*gradt_dpdr(i,103) - gradt_dpdr(i,150) - &
                    & gradt_dpdr(i,149) - gradt_dpdr(i,167)
                gradt_dpdr(i,177) = gradt_dpdr(i,2)*epot_p(61) + epot_p(2)*gradt_dpdr(i,61) - gradt_dpdr(i,132) - &
                    & gradt_dpdr(i,154)
                gradt_dpdr(i,178) = gradt_dpdr(i,2)*epot_p(68) + epot_p(2)*gradt_dpdr(i,68) - gradt_dpdr(i,159) - &
                    & gradt_dpdr(i,174)
                gradt_dpdr(i,179) = gradt_dpdr(i,3)*epot_p(62) + epot_p(3)*gradt_dpdr(i,62) - gradt_dpdr(i,132) - &
                    & gradt_dpdr(i,163) - gradt_dpdr(i,156) - gradt_dpdr(i,158) - gradt_dpdr(i,154)
                gradt_dpdr(i,180) = gradt_dpdr(i,6)*epot_p(38) + epot_p(6)*gradt_dpdr(i,38) - gradt_dpdr(i,132) - &
                    & gradt_dpdr(i,163) - gradt_dpdr(i,157)
                gradt_dpdr(i,181) = gradt_dpdr(i,7)*epot_p(38) + epot_p(7)*gradt_dpdr(i,38) - gradt_dpdr(i,132) - &
                    & gradt_dpdr(i,163) - gradt_dpdr(i,155)
                gradt_dpdr(i,182) = gradt_dpdr(i,2)*epot_p(70) + epot_p(2)*gradt_dpdr(i,70) - gradt_dpdr(i,163) - &
                    & gradt_dpdr(i,179)
                gradt_dpdr(i,183) = gradt_dpdr(i,1)*epot_p(110) + epot_p(1)*gradt_dpdr(i,110) - gradt_dpdr(i,163) - &
                    & gradt_dpdr(i,160) - gradt_dpdr(i,179)
                gradt_dpdr(i,184) = gradt_dpdr(i,6)*epot_p(39) + epot_p(6)*gradt_dpdr(i,39) - gradt_dpdr(i,162)
                gradt_dpdr(i,185) = gradt_dpdr(i,7)*epot_p(39) + epot_p(7)*gradt_dpdr(i,39) - gradt_dpdr(i,161)
                gradt_dpdr(i,186) = gradt_dpdr(i,2)*epot_p(65) + epot_p(2)*gradt_dpdr(i,65) - gradt_dpdr(i,136)
                gradt_dpdr(i,187) = gradt_dpdr(i,3)*epot_p(66) + epot_p(3)*gradt_dpdr(i,66) - gradt_dpdr(i,148) - &
                    & gradt_dpdr(i,149) - gradt_dpdr(i,147) - gradt_dpdr(i,175) - gradt_dpdr(i,176) - gradt_dpdr(i,174)
                gradt_dpdr(i,188) = gradt_dpdr(i,2)*epot_p(67) + epot_p(2)*gradt_dpdr(i,67) - gradt_dpdr(i,151)
                gradt_dpdr(i,189) = gradt_dpdr(i,4)*epot_p(66) + epot_p(4)*gradt_dpdr(i,66) - gradt_dpdr(i,166) - &
                    & gradt_dpdr(i,167) - gradt_dpdr(i,165)
                gradt_dpdr(i,190) = gradt_dpdr(i,2)*epot_p(69) + epot_p(2)*gradt_dpdr(i,69) - gradt_dpdr(i,160)
                gradt_dpdr(i,191) = gradt_dpdr(i,18)*epot_p(21) + epot_p(18)*gradt_dpdr(i,21) - gradt_dpdr(i,171) - &
                    & gradt_dpdr(i,172) - gradt_dpdr(i,169)
                gradt_dpdr(i,192) = gradt_dpdr(i,2)*epot_p(71) + epot_p(2)*gradt_dpdr(i,71) - gradt_dpdr(i,183)
                gradt_dpdr(i,193) = gradt_dpdr(i,3)*epot_p(71) + epot_p(3)*gradt_dpdr(i,71) - gradt_dpdr(i,184) - &
                    & gradt_dpdr(i,185) - gradt_dpdr(i,182)
                gradt_dpdr(i,194) = gradt_dpdr(i,1)*epot_p(119) + epot_p(1)*gradt_dpdr(i,119) - gradt_dpdr(i,193) - &
                    & gradt_dpdr(i,192)
                gradt_dpdr(i,195) = gradt_dpdr(i,40)*epot_p(5) + epot_p(40)*gradt_dpdr(i,5)
                gradt_dpdr(i,196) = gradt_dpdr(i,40)*epot_p(6) + epot_p(40)*gradt_dpdr(i,6)
                gradt_dpdr(i,197) = gradt_dpdr(i,40)*epot_p(7) + epot_p(40)*gradt_dpdr(i,7)
                gradt_dpdr(i,198) = gradt_dpdr(i,40)*epot_p(8) + epot_p(40)*gradt_dpdr(i,8)
                gradt_dpdr(i,199) = gradt_dpdr(i,40)*epot_p(9) + epot_p(40)*gradt_dpdr(i,9)
                gradt_dpdr(i,200) = gradt_dpdr(i,40)*epot_p(10) + epot_p(40)*gradt_dpdr(i,10)
                gradt_dpdr(i,201) = gradt_dpdr(i,11)*epot_p(22) + epot_p(11)*gradt_dpdr(i,22) - gradt_dpdr(i,195)
                gradt_dpdr(i,202) = gradt_dpdr(i,12)*epot_p(22) + epot_p(12)*gradt_dpdr(i,22) - gradt_dpdr(i,196) - &
                    & gradt_dpdr(i,197) - gradt_dpdr(i,195) - gradt_dpdr(i,196) - gradt_dpdr(i,197) - gradt_dpdr(i, &
                    & 195) - gradt_dpdr(i,196) - gradt_dpdr(i,197)
                gradt_dpdr(i,203) = gradt_dpdr(i,17)*epot_p(22) + epot_p(17)*gradt_dpdr(i,22) - gradt_dpdr(i,198)
                gradt_dpdr(i,204) = gradt_dpdr(i,6)*epot_p(43) + epot_p(6)*gradt_dpdr(i,43)
                gradt_dpdr(i,205) = gradt_dpdr(i,7)*epot_p(43) + epot_p(7)*gradt_dpdr(i,43)
                gradt_dpdr(i,206) = gradt_dpdr(i,11)*epot_p(26) + epot_p(11)*gradt_dpdr(i,26) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,202)
                gradt_dpdr(i,207) = gradt_dpdr(i,2)*epot_p(76) + epot_p(2)*gradt_dpdr(i,76) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,202)
                gradt_dpdr(i,208) = gradt_dpdr(i,2)*epot_p(78) + epot_p(2)*gradt_dpdr(i,78) - gradt_dpdr(i,200)
                gradt_dpdr(i,209) = gradt_dpdr(i,6)*epot_p(41) + epot_p(6)*gradt_dpdr(i,41) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,195) - gradt_dpdr(i,202) - gradt_dpdr(i,195)
                gradt_dpdr(i,210) = gradt_dpdr(i,6)*epot_p(42) + epot_p(6)*gradt_dpdr(i,42) - gradt_dpdr(i,197) - &
                    & gradt_dpdr(i,197) - gradt_dpdr(i,197)
                gradt_dpdr(i,211) = gradt_dpdr(i,7)*epot_p(41) + epot_p(7)*gradt_dpdr(i,41) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,195) - gradt_dpdr(i,202) - gradt_dpdr(i,195)
                gradt_dpdr(i,212) = gradt_dpdr(i,7)*epot_p(42) + epot_p(7)*gradt_dpdr(i,42) - gradt_dpdr(i,196) - &
                    & gradt_dpdr(i,196) - gradt_dpdr(i,196)
                gradt_dpdr(i,213) = gradt_dpdr(i,11)*epot_p(29) + epot_p(11)*gradt_dpdr(i,29) - gradt_dpdr(i,209)
                gradt_dpdr(i,214) = gradt_dpdr(i,11)*epot_p(30) + epot_p(11)*gradt_dpdr(i,30) - gradt_dpdr(i,211)
                gradt_dpdr(i,215) = gradt_dpdr(i,18)*epot_p(22) + epot_p(18)*gradt_dpdr(i,22) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,213) - gradt_dpdr(i,214)
                gradt_dpdr(i,216) = gradt_dpdr(i,2)*epot_p(77) + epot_p(2)*gradt_dpdr(i,77) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,206)
                gradt_dpdr(i,217) = gradt_dpdr(i,6)*epot_p(49) + epot_p(6)*gradt_dpdr(i,49) - gradt_dpdr(i,205)
                gradt_dpdr(i,218) = gradt_dpdr(i,7)*epot_p(49) + epot_p(7)*gradt_dpdr(i,49) - gradt_dpdr(i,204)
                gradt_dpdr(i,219) = gradt_dpdr(i,3)*epot_p(77) + epot_p(3)*gradt_dpdr(i,77) - gradt_dpdr(i,200) - &
                    & gradt_dpdr(i,199) - gradt_dpdr(i,198) - gradt_dpdr(i,209) - gradt_dpdr(i,217) - gradt_dpdr(i, &
                    & 211) - gradt_dpdr(i,218) - gradt_dpdr(i,207) - gradt_dpdr(i,204) - gradt_dpdr(i,205) - &
                    & gradt_dpdr(i,200) - gradt_dpdr(i,199) - gradt_dpdr(i,198) - gradt_dpdr(i,200) - gradt_dpdr(i, &
                    & 198) - gradt_dpdr(i,200) - gradt_dpdr(i,198)
                gradt_dpdr(i,220) = gradt_dpdr(i,3)*epot_p(78) + epot_p(3)*gradt_dpdr(i,78) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,210) - gradt_dpdr(i,212)
                gradt_dpdr(i,221) = gradt_dpdr(i,10)*epot_p(41) + epot_p(10)*gradt_dpdr(i,41) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,220) - gradt_dpdr(i,217) - gradt_dpdr(i,218) - gradt_dpdr(i,216)
                gradt_dpdr(i,222) = gradt_dpdr(i,1)*epot_p(133) + epot_p(1)*gradt_dpdr(i,133) - gradt_dpdr(i,200) - &
                    & gradt_dpdr(i,220)
                gradt_dpdr(i,223) = gradt_dpdr(i,11)*epot_p(27) + epot_p(11)*gradt_dpdr(i,27) - gradt_dpdr(i,195) - &
                    & gradt_dpdr(i,203)
                gradt_dpdr(i,224) = gradt_dpdr(i,1)*epot_p(134) + epot_p(1)*gradt_dpdr(i,134) - gradt_dpdr(i,201)
                gradt_dpdr(i,225) = gradt_dpdr(i,2)*epot_p(81) + epot_p(2)*gradt_dpdr(i,81) - gradt_dpdr(i,209)
                gradt_dpdr(i,226) = gradt_dpdr(i,2)*epot_p(80) + epot_p(2)*gradt_dpdr(i,80) - gradt_dpdr(i,202) - &
                    & gradt_dpdr(i,206)
                gradt_dpdr(i,227) = gradt_dpdr(i,2)*epot_p(82) + epot_p(2)*gradt_dpdr(i,82) - gradt_dpdr(i,211)
                gradt_dpdr(i,228) = gradt_dpdr(i,2)*epot_p(91) + epot_p(2)*gradt_dpdr(i,91) - gradt_dpdr(i,215) - &
                    & gradt_dpdr(i,219)
                gradt_dpdr(i,229) = gradt_dpdr(i,11)*epot_p(28) + epot_p(11)*gradt_dpdr(i,28) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,207)
                gradt_dpdr(i,230) = gradt_dpdr(i,6)*epot_p(53) + epot_p(6)*gradt_dpdr(i,53) - gradt_dpdr(i,205)
                gradt_dpdr(i,231) = gradt_dpdr(i,5)*epot_p(54) + epot_p(5)*gradt_dpdr(i,54) - gradt_dpdr(i,209)
                gradt_dpdr(i,232) = gradt_dpdr(i,7)*epot_p(53) + epot_p(7)*gradt_dpdr(i,53) - gradt_dpdr(i,204)
                gradt_dpdr(i,233) = gradt_dpdr(i,7)*epot_p(54) + epot_p(7)*gradt_dpdr(i,54) - gradt_dpdr(i,196)
                gradt_dpdr(i,234) = gradt_dpdr(i,5)*epot_p(55) + epot_p(5)*gradt_dpdr(i,55) - gradt_dpdr(i,211)
                gradt_dpdr(i,235) = gradt_dpdr(i,6)*epot_p(55) + epot_p(6)*gradt_dpdr(i,55) - gradt_dpdr(i,197)
                gradt_dpdr(i,236) = gradt_dpdr(i,11)*epot_p(35) + epot_p(11)*gradt_dpdr(i,35) - gradt_dpdr(i,198)
                gradt_dpdr(i,237) = gradt_dpdr(i,6)*epot_p(65) + epot_p(6)*gradt_dpdr(i,65)
                gradt_dpdr(i,238) = gradt_dpdr(i,7)*epot_p(65) + epot_p(7)*gradt_dpdr(i,65)
                gradt_dpdr(i,239) = gradt_dpdr(i,2)*epot_p(86) + epot_p(2)*gradt_dpdr(i,86) - gradt_dpdr(i,199) - &
                    & gradt_dpdr(i,229)
                gradt_dpdr(i,240) = gradt_dpdr(i,1)*epot_p(139) + epot_p(1)*gradt_dpdr(i,139) - gradt_dpdr(i,206) - &
                    & gradt_dpdr(i,229) - gradt_dpdr(i,224)
                gradt_dpdr(i,241) = gradt_dpdr(i,17)*epot_p(29) + epot_p(17)*gradt_dpdr(i,29) - gradt_dpdr(i,225)
                gradt_dpdr(i,242) = gradt_dpdr(i,2)*epot_p(99) + epot_p(2)*gradt_dpdr(i,99) - gradt_dpdr(i,231)
                gradt_dpdr(i,243) = gradt_dpdr(i,17)*epot_p(30) + epot_p(17)*gradt_dpdr(i,30) - gradt_dpdr(i,227)
                gradt_dpdr(i,244) = gradt_dpdr(i,2)*epot_p(95) + epot_p(2)*gradt_dpdr(i,95) - gradt_dpdr(i,220) - &
                    & gradt_dpdr(i,221)
                gradt_dpdr(i,245) = gradt_dpdr(i,4)*epot_p(81) + epot_p(4)*gradt_dpdr(i,81) - gradt_dpdr(i,202) - &
                    & gradt_dpdr(i,242) - gradt_dpdr(i,227)
                gradt_dpdr(i,246) = gradt_dpdr(i,2)*epot_p(100) + epot_p(2)*gradt_dpdr(i,100) - gradt_dpdr(i,234)
                gradt_dpdr(i,247) = gradt_dpdr(i,4)*epot_p(82) + epot_p(4)*gradt_dpdr(i,82) - gradt_dpdr(i,202) - &
                    & gradt_dpdr(i,225) - gradt_dpdr(i,246)
                gradt_dpdr(i,248) = gradt_dpdr(i,11)*epot_p(36) + epot_p(11)*gradt_dpdr(i,36) - gradt_dpdr(i,215) - &
                    & gradt_dpdr(i,219)
                gradt_dpdr(i,249) = gradt_dpdr(i,2)*epot_p(102) + epot_p(2)*gradt_dpdr(i,102) - gradt_dpdr(i,233)
                gradt_dpdr(i,250) = gradt_dpdr(i,6)*epot_p(58) + epot_p(6)*gradt_dpdr(i,58) - gradt_dpdr(i,206) - &
                    & gradt_dpdr(i,214) - gradt_dpdr(i,244) - gradt_dpdr(i,214)
                gradt_dpdr(i,251) = gradt_dpdr(i,2)*epot_p(103) + epot_p(2)*gradt_dpdr(i,103) - gradt_dpdr(i,235)
                gradt_dpdr(i,252) = gradt_dpdr(i,7)*epot_p(58) + epot_p(7)*gradt_dpdr(i,58) - gradt_dpdr(i,213) - &
                    & gradt_dpdr(i,206) - gradt_dpdr(i,244) - gradt_dpdr(i,213)
                gradt_dpdr(i,253) = gradt_dpdr(i,1)*epot_p(150) + epot_p(1)*gradt_dpdr(i,150) - gradt_dpdr(i,215) - &
                    & gradt_dpdr(i,233) - gradt_dpdr(i,235)
                gradt_dpdr(i,254) = gradt_dpdr(i,11)*epot_p(37) + epot_p(11)*gradt_dpdr(i,37) - gradt_dpdr(i,200)
                gradt_dpdr(i,255) = gradt_dpdr(i,2)*epot_p(93) + epot_p(2)*gradt_dpdr(i,93) - gradt_dpdr(i,217)
                gradt_dpdr(i,256) = gradt_dpdr(i,2)*epot_p(94) + epot_p(2)*gradt_dpdr(i,94) - gradt_dpdr(i,218)
                gradt_dpdr(i,257) = gradt_dpdr(i,11)*epot_p(38) + epot_p(11)*gradt_dpdr(i,38) - gradt_dpdr(i,220) - &
                    & gradt_dpdr(i,221)
                gradt_dpdr(i,258) = gradt_dpdr(i,2)*epot_p(107) + epot_p(2)*gradt_dpdr(i,107) - gradt_dpdr(i,245) - &
                    & gradt_dpdr(i,250)
                gradt_dpdr(i,259) = gradt_dpdr(i,6)*epot_p(62) + epot_p(6)*gradt_dpdr(i,62) - gradt_dpdr(i,218)
                gradt_dpdr(i,260) = gradt_dpdr(i,2)*epot_p(108) + epot_p(2)*gradt_dpdr(i,108) - gradt_dpdr(i,247) - &
                    & gradt_dpdr(i,252)
                gradt_dpdr(i,261) = gradt_dpdr(i,7)*epot_p(62) + epot_p(7)*gradt_dpdr(i,62) - gradt_dpdr(i,217)
                gradt_dpdr(i,262) = gradt_dpdr(i,6)*epot_p(64) + epot_p(6)*gradt_dpdr(i,64) - gradt_dpdr(i,200) - &
                    & gradt_dpdr(i,246) - gradt_dpdr(i,200)
                gradt_dpdr(i,263) = gradt_dpdr(i,11)*epot_p(39) + epot_p(11)*gradt_dpdr(i,39) - gradt_dpdr(i,222)
                gradt_dpdr(i,264) = gradt_dpdr(i,2)*epot_p(111) + epot_p(2)*gradt_dpdr(i,111) - gradt_dpdr(i,259)
                gradt_dpdr(i,265) = gradt_dpdr(i,2)*epot_p(112) + epot_p(2)*gradt_dpdr(i,112) - gradt_dpdr(i,261)
                gradt_dpdr(i,266) = gradt_dpdr(i,12)*epot_p(39) + epot_p(12)*gradt_dpdr(i,39) - gradt_dpdr(i,221) - &
                    & gradt_dpdr(i,264) - gradt_dpdr(i,265)
                gradt_dpdr(i,267) = gradt_dpdr(i,2)*epot_p(98) + epot_p(2)*gradt_dpdr(i,98) - gradt_dpdr(i,208) - &
                    & gradt_dpdr(i,240)
                gradt_dpdr(i,268) = gradt_dpdr(i,6)*epot_p(54) + epot_p(6)*gradt_dpdr(i,54) - gradt_dpdr(i,210)
                gradt_dpdr(i,269) = gradt_dpdr(i,7)*epot_p(55) + epot_p(7)*gradt_dpdr(i,55) - gradt_dpdr(i,212)
                gradt_dpdr(i,270) = gradt_dpdr(i,2)*epot_p(96) + epot_p(2)*gradt_dpdr(i,96) - gradt_dpdr(i,203) - &
                    & gradt_dpdr(i,223)
                gradt_dpdr(i,271) = gradt_dpdr(i,2)*epot_p(97) + epot_p(2)*gradt_dpdr(i,97) - gradt_dpdr(i,207) - &
                    & gradt_dpdr(i,229)
                gradt_dpdr(i,272) = gradt_dpdr(i,2)*epot_p(101) + epot_p(2)*gradt_dpdr(i,101) - gradt_dpdr(i,215) - &
                    & gradt_dpdr(i,248)
                gradt_dpdr(i,273) = gradt_dpdr(i,2)*epot_p(106) + epot_p(2)*gradt_dpdr(i,106) - gradt_dpdr(i,220) - &
                    & gradt_dpdr(i,257)
                gradt_dpdr(i,274) = gradt_dpdr(i,6)*epot_p(59) + epot_p(6)*gradt_dpdr(i,59) - gradt_dpdr(i,220) - &
                    & gradt_dpdr(i,215) - gradt_dpdr(i,209)
                gradt_dpdr(i,275) = gradt_dpdr(i,7)*epot_p(60) + epot_p(7)*gradt_dpdr(i,60) - gradt_dpdr(i,220) - &
                    & gradt_dpdr(i,215) - gradt_dpdr(i,211)
                gradt_dpdr(i,276) = gradt_dpdr(i,5)*epot_p(66) + epot_p(5)*gradt_dpdr(i,66) - gradt_dpdr(i,215) - &
                    & gradt_dpdr(i,253) - gradt_dpdr(i,249) - gradt_dpdr(i,251) - gradt_dpdr(i,248)
                gradt_dpdr(i,277) = gradt_dpdr(i,6)*epot_p(66) + epot_p(6)*gradt_dpdr(i,66) - gradt_dpdr(i,213) - &
                    & gradt_dpdr(i,252)
                gradt_dpdr(i,278) = gradt_dpdr(i,7)*epot_p(66) + epot_p(7)*gradt_dpdr(i,66) - gradt_dpdr(i,214) - &
                    & gradt_dpdr(i,250)
                gradt_dpdr(i,279) = gradt_dpdr(i,2)*epot_p(104) + epot_p(2)*gradt_dpdr(i,104) - gradt_dpdr(i,216) - &
                    & gradt_dpdr(i,239)
                gradt_dpdr(i,280) = gradt_dpdr(i,2)*epot_p(105) + epot_p(2)*gradt_dpdr(i,105) - gradt_dpdr(i,219) - &
                    & gradt_dpdr(i,248)
                gradt_dpdr(i,281) = gradt_dpdr(i,1)*epot_p(170) + epot_p(1)*gradt_dpdr(i,170) - gradt_dpdr(i,244) - &
                    & gradt_dpdr(i,240) - gradt_dpdr(i,273)
                gradt_dpdr(i,282) = gradt_dpdr(i,10)*epot_p(54) + epot_p(10)*gradt_dpdr(i,54) - gradt_dpdr(i,227)
                gradt_dpdr(i,283) = gradt_dpdr(i,10)*epot_p(55) + epot_p(10)*gradt_dpdr(i,55) - gradt_dpdr(i,225)
                gradt_dpdr(i,284) = gradt_dpdr(i,2)*epot_p(114) + epot_p(2)*gradt_dpdr(i,114) - gradt_dpdr(i,253) - &
                    & gradt_dpdr(i,276)
                gradt_dpdr(i,285) = gradt_dpdr(i,1)*epot_p(174) + epot_p(1)*gradt_dpdr(i,174) - gradt_dpdr(i,250) - &
                    & gradt_dpdr(i,252) - gradt_dpdr(i,248) - gradt_dpdr(i,276) - gradt_dpdr(i,273)
                gradt_dpdr(i,286) = gradt_dpdr(i,6)*epot_p(68) + epot_p(6)*gradt_dpdr(i,68) - gradt_dpdr(i,217) - &
                    & gradt_dpdr(i,261) - gradt_dpdr(i,251)
                gradt_dpdr(i,287) = gradt_dpdr(i,7)*epot_p(68) + epot_p(7)*gradt_dpdr(i,68) - gradt_dpdr(i,218) - &
                    & gradt_dpdr(i,259) - gradt_dpdr(i,249)
                gradt_dpdr(i,288) = gradt_dpdr(i,2)*epot_p(109) + epot_p(2)*gradt_dpdr(i,109) - gradt_dpdr(i,221) - &
                    & gradt_dpdr(i,257)
                gradt_dpdr(i,289) = gradt_dpdr(i,2)*epot_p(116) + epot_p(2)*gradt_dpdr(i,116) - gradt_dpdr(i,262) - &
                    & gradt_dpdr(i,285)
                gradt_dpdr(i,290) = gradt_dpdr(i,3)*epot_p(110) + epot_p(3)*gradt_dpdr(i,110) - gradt_dpdr(i,221) - &
                    & gradt_dpdr(i,266) - gradt_dpdr(i,259) - gradt_dpdr(i,261) - gradt_dpdr(i,257)
                gradt_dpdr(i,291) = gradt_dpdr(i,6)*epot_p(70) + epot_p(6)*gradt_dpdr(i,70) - gradt_dpdr(i,221) - &
                    & gradt_dpdr(i,266) - gradt_dpdr(i,260)
                gradt_dpdr(i,292) = gradt_dpdr(i,7)*epot_p(70) + epot_p(7)*gradt_dpdr(i,70) - gradt_dpdr(i,221) - &
                    & gradt_dpdr(i,266) - gradt_dpdr(i,258)
                gradt_dpdr(i,293) = gradt_dpdr(i,2)*epot_p(118) + epot_p(2)*gradt_dpdr(i,118) - gradt_dpdr(i,266) - &
                    & gradt_dpdr(i,290)
                gradt_dpdr(i,294) = gradt_dpdr(i,1)*epot_p(183) + epot_p(1)*gradt_dpdr(i,183) - gradt_dpdr(i,266) - &
                    & gradt_dpdr(i,263) - gradt_dpdr(i,290)
                gradt_dpdr(i,295) = gradt_dpdr(i,6)*epot_p(71) + epot_p(6)*gradt_dpdr(i,71) - gradt_dpdr(i,265)
                gradt_dpdr(i,296) = gradt_dpdr(i,7)*epot_p(71) + epot_p(7)*gradt_dpdr(i,71) - gradt_dpdr(i,264)
                gradt_dpdr(i,297) = gradt_dpdr(i,1)*epot_p(186) + epot_p(1)*gradt_dpdr(i,186) - gradt_dpdr(i,270)
                gradt_dpdr(i,298) = gradt_dpdr(i,1)*epot_p(187) + epot_p(1)*gradt_dpdr(i,187) - gradt_dpdr(i,277) - &
                    & gradt_dpdr(i,278) - gradt_dpdr(i,276)
                gradt_dpdr(i,299) = gradt_dpdr(i,2)*epot_p(115) + epot_p(2)*gradt_dpdr(i,115) - gradt_dpdr(i,254)
                gradt_dpdr(i,300) = gradt_dpdr(i,1)*epot_p(189) + epot_p(1)*gradt_dpdr(i,189) - gradt_dpdr(i,286) - &
                    & gradt_dpdr(i,287) - gradt_dpdr(i,285) - gradt_dpdr(i,284) - gradt_dpdr(i,298)
                gradt_dpdr(i,301) = gradt_dpdr(i,2)*epot_p(117) + epot_p(2)*gradt_dpdr(i,117) - gradt_dpdr(i,263)
                gradt_dpdr(i,302) = gradt_dpdr(i,18)*epot_p(39) + epot_p(18)*gradt_dpdr(i,39) - gradt_dpdr(i,282) - &
                    & gradt_dpdr(i,283) - gradt_dpdr(i,280)
                gradt_dpdr(i,303) = gradt_dpdr(i,2)*epot_p(119) + epot_p(2)*gradt_dpdr(i,119) - gradt_dpdr(i,294)
                gradt_dpdr(i,304) = gradt_dpdr(i,3)*epot_p(119) + epot_p(3)*gradt_dpdr(i,119) - gradt_dpdr(i,295) - &
                    & gradt_dpdr(i,296) - gradt_dpdr(i,293)
                gradt_dpdr(i,305) = gradt_dpdr(i,1)*epot_p(194) + epot_p(1)*gradt_dpdr(i,194) - gradt_dpdr(i,304) - &
                    & gradt_dpdr(i,303)
            enddo
        end subroutine evdpdr
        !> @brief Evaluate evdbdr for the N4_1A_PIP_2013 potential.
        subroutine evdbdr
            !  The subroutine eliminate the 2-body terms in Bowman's approach
            implicit none
            integer          :: j
            integer          :: i
            double precision :: db1dr(6,306)
            ! Pass P(0:305) to BM1(1:306)
            do j=1, 6
                do i=1, 306
                    db1dr(j,i) = gradt_dpdr(j,i-1)
                enddo
            enddo
            ! Remove unconnected terms and 2-body terms and pass to B(1:276)
            do j=1, 6
                gradt_dbdr(j,1) = db1dr(j,4)
                do i=2, 4
                    gradt_dbdr(j,i) = db1dr(j,i+4)
                enddo
                gradt_dbdr(j,5) = db1dr(j,10)
                do i=6, 11
                    gradt_dbdr(j,i) = db1dr(j,i+6)
                enddo
                gradt_dbdr(j,12) = db1dr(j,19)
                gradt_dbdr(j,13) = db1dr(j,21)
                do i=14, 26
                    gradt_dbdr(j,i) = db1dr(j,i+9)
                enddo
                gradt_dbdr(j,27) = db1dr(j,37)
                gradt_dbdr(j,28) = db1dr(j,39)
                do i=29, 53
                    gradt_dbdr(j,i) = db1dr(j,i+12)
                enddo
                gradt_dbdr(j,54) = db1dr(j,67)
                gradt_dbdr(j,55) = db1dr(j,69)
                gradt_dbdr(j,56) = db1dr(j,71)
                do i=57, 97
                    gradt_dbdr(j,i) = db1dr(j,i+16)
                enddo
                gradt_dbdr(j,98)  = db1dr(j,115)
                gradt_dbdr(j,99)  = db1dr(j,117)
                gradt_dbdr(j,100) = db1dr(j,119)
                do i=101, 166
                    gradt_dbdr(j,i) = db1dr(j,i+20)
                enddo
                gradt_dbdr(j,167) = db1dr(j,188)
                gradt_dbdr(j,168) = db1dr(j,190)
                gradt_dbdr(j,169) = db1dr(j,192)
                gradt_dbdr(j,170) = db1dr(j,194)
                do i=171, 272
                    gradt_dbdr(j,i) = db1dr(j,i+25)
                enddo
                gradt_dbdr(j,273) = db1dr(j,299)
                gradt_dbdr(j,274) = db1dr(j,301)
                gradt_dbdr(j,275) = db1dr(j,303)
                gradt_dbdr(j,276) = db1dr(j,305)
            enddo
        end subroutine evdbdr
    end subroutine evaluate
end module chempotpy_n4_1a_pip_2013_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_n4_1a_pip_2013_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 4
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
