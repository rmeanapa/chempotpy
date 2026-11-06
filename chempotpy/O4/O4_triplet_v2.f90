!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!   System:                        O4
!   Functional form:               permutation-invariant polynomials
!   Common name:                   O4 triplet(adiabatic ground state)
!   Number of derivatives:         1
!   Number of bodies:              4
!   Number of electronic surfaces: 1
!   Interface:                     Section-2
!
!   Reference - "Potential energy surface of triplet O4," Yuliya
!   Paukku, Zoltan Varga, and Donald G. Truhlar, to be published
!   (Feb. 09, 2018)
!
!   Notes:    PES of triplet O4 with special emphasis on
!                  O2 + O2 -> O2 + O + O
!            - Fit is based on a dataset of 10180 points
!            - Mixed-exponential-Gaussian (MEG)
!               variables are applied to describe
!               the long-range interactions
!            - D3 dispersion correction is implemented
!               for pairwise interactions
!
!     O1--O2
!
!     O3--O4
!
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_o4_triplet_v2_pes
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
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: coord_r(6)
        double precision              :: epot_b(430)
        double precision              :: epot_c(430)
        double precision              :: epot_p(0:465)
        double precision              :: epot_rm(0:111)
        double precision              :: epot_rms(6)
        double precision              :: gradt_dbdr(6,430)
        double precision              :: gradt_dmdr(6,0:111)
        double precision              :: gradt_dmsdr(6,6)
        double precision              :: gradt_dpdr(6,0:465)
        double precision              :: gradt_drdx(6,12)
        double precision              :: gradt_dvdr(6)
        double precision              :: msprmt_a
        double precision              :: msprmt_ab
        double precision              :: msprmt_ra
        double precision              :: msprmt_rb
        double precision              :: v, tx(12), tg(12)
        integer                       :: iatom, idir, j, istate
        !initialize
        coord_r      = 0.d0; epot_b = 0.d0; epot_c = 0.d0; epot_p = 0.d0; epot_rm = 0.d0; epot_rms = 0.d0; gradt_dbdr = 0.d0
        gradt_dmdr   = 0.d0; gradt_dmsdr = 0.d0; gradt_dpdr = 0.d0; gradt_drdx = 0.d0; gradt_dvdr = 0.d0; msprmt_a = 0.d0
        msprmt_ab    = 0.d0; msprmt_ra = 0.d0; msprmt_rb = 0.d0; msprmt_a = 0.9d0; msprmt_ab = 4.5d0; msprmt_ra = 1.3d0
        msprmt_rb    = 1.1d0
        epot_c(   1) = 0.434720431798d+02
        epot_c(   2) = -0.554668579283d+03
        epot_c(   3) = 0.257445901353d+04
        epot_c(   4) = -0.213807688288d+04
        epot_c(   5) = -0.248886942161d+03
        epot_c(   6) = 0.719084654648d+03
        epot_c(   7) = -0.785231200932d+03
        epot_c(   8) = 0.268844019872d+04
        epot_c(   9) = -0.894783391739d+03
        epot_c(  10) = -0.357957764652d+04
        epot_c(  11) = 0.362675454825d+04
        epot_c(  12) = 0.102629622676d+04
        epot_c(  13) = 0.795258994983d+03
        epot_c(  14) = -0.946936811683d+04
        epot_c(  15) = 0.571579383749d+04
        epot_c(  16) = 0.934372178824d+03
        epot_c(  17) = -0.292732443657d+04
        epot_c(  18) = 0.219656326428d+04
        epot_c(  19) = -0.676509263255d+04
        epot_c(  20) = -0.282505207332d+04
        epot_c(  21) = 0.104924792198d+05
        epot_c(  22) = -0.658586513846d+04
        epot_c(  23) = -0.432459822098d+04
        epot_c(  24) = 0.316198143890d+04
        epot_c(  25) = -0.337163445639d+04
        epot_c(  26) = -0.713571319552d+03
        epot_c(  27) = -0.134555448022d+04
        epot_c(  28) = -0.150137169111d+04
        epot_c(  29) = -0.358000122745d+04
        epot_c(  30) = 0.140590867542d+05
        epot_c(  31) = -0.459125633418d+04
        epot_c(  32) = -0.550769788422d+04
        epot_c(  33) = 0.217230838750d+04
        epot_c(  34) = -0.632774344127d+03
        epot_c(  35) = -0.100999591686d+05
        epot_c(  36) = 0.169165667764d+04
        epot_c(  37) = 0.167842440492d+04
        epot_c(  38) = -0.257911596712d+05
        epot_c(  39) = 0.127815645611d+05
        epot_c(  40) = -0.537883442310d+03
        epot_c(  41) = -0.133764280242d+04
        epot_c(  42) = -0.212370345146d+04
        epot_c(  43) = 0.161754184625d+05
        epot_c(  44) = 0.355956338618d+04
        epot_c(  45) = 0.913378256342d+04
        epot_c(  46) = 0.610950125310d+04
        epot_c(  47) = 0.515181004953d+04
        epot_c(  48) = -0.149547604838d+05
        epot_c(  49) = 0.447834665993d+04
        epot_c(  50) = 0.289335443481d+04
        epot_c(  51) = -0.430236737245d+04
        epot_c(  52) = 0.115862975164d+05
        epot_c(  53) = -0.223649282872d+04
        epot_c(  54) = -0.424281221393d+04
        epot_c(  55) = 0.607690277921d+04
        epot_c(  56) = 0.682307953653d+03
        epot_c(  57) = -0.146643434299d+05
        epot_c(  58) = -0.995963785829d+04
        epot_c(  59) = 0.151656352036d+05
        epot_c(  60) = -0.713080950000d+04
        epot_c(  61) = 0.992264149235d+04
        epot_c(  62) = -0.225192674918d+05
        epot_c(  63) = -0.227521700630d+05
        epot_c(  64) = 0.356736751019d+03
        epot_c(  65) = -0.149009596280d+05
        epot_c(  66) = 0.547202199914d+04
        epot_c(  67) = -0.562297176466d+04
        epot_c(  68) = 0.533482100629d+04
        epot_c(  69) = 0.173076951226d+05
        epot_c(  70) = -0.470318017955d+04
        epot_c(  71) = 0.252423648786d+04
        epot_c(  72) = 0.397204506011d+04
        epot_c(  73) = 0.149604469867d+05
        epot_c(  74) = 0.261468562194d+04
        epot_c(  75) = 0.148832634817d+04
        epot_c(  76) = -0.525471073088d+04
        epot_c(  77) = 0.393212581141d+05
        epot_c(  78) = -0.306109210725d+05
        epot_c(  79) = 0.721202884138d+04
        epot_c(  80) = 0.352945175290d+04
        epot_c(  81) = -0.784511969410d+04
        epot_c(  82) = 0.201473287242d+03
        epot_c(  83) = 0.124609448656d+05
        epot_c(  84) = -0.103099257764d+05
        epot_c(  85) = 0.123788976352d+04
        epot_c(  86) = -0.845440125549d+03
        epot_c(  87) = -0.652328346462d+04
        epot_c(  88) = -0.397970507273d+04
        epot_c(  89) = -0.813327575032d+04
        epot_c(  90) = -0.854960333368d+04
        epot_c(  91) = -0.146818322655d+05
        epot_c(  92) = 0.250112509218d+05
        epot_c(  93) = -0.120475653865d+04
        epot_c(  94) = 0.894965121433d+03
        epot_c(  95) = 0.687552457372d+04
        epot_c(  96) = -0.155946821955d+05
        epot_c(  97) = -0.556039535294d+03
        epot_c(  98) = 0.412239459401d+04
        epot_c(  99) = -0.103729683198d+05
        epot_c( 100) = 0.120021137587d+04
        epot_c( 101) = -0.106685046856d+05
        epot_c( 102) = 0.580461159121d+04
        epot_c( 103) = 0.323628855252d+05
        epot_c( 104) = 0.801279078904d+03
        epot_c( 105) = -0.144743305379d+05
        epot_c( 106) = 0.245891704416d+05
        epot_c( 107) = 0.734963482362d+04
        epot_c( 108) = 0.579008306105d+04
        epot_c( 109) = -0.209003437783d+05
        epot_c( 110) = 0.788390162743d+03
        epot_c( 111) = 0.795154017688d+03
        epot_c( 112) = -0.218632472140d+05
        epot_c( 113) = 0.172512443045d+05
        epot_c( 114) = 0.540100143425d+05
        epot_c( 115) = 0.338674945681d+04
        epot_c( 116) = -0.386001569801d+05
        epot_c( 117) = 0.239712386089d+04
        epot_c( 118) = -0.613824367400d+04
        epot_c( 119) = -0.488209303900d+04
        epot_c( 120) = 0.701192258148d+04
        epot_c( 121) = 0.167984032510d+05
        epot_c( 122) = -0.259650106547d+05
        epot_c( 123) = 0.177986468669d+05
        epot_c( 124) = 0.150944386869d+05
        epot_c( 125) = 0.121802007803d+05
        epot_c( 126) = -0.978234467399d+04
        epot_c( 127) = 0.184330692467d+05
        epot_c( 128) = 0.921531962740d+04
        epot_c( 129) = 0.183336479381d+05
        epot_c( 130) = 0.587334293511d+02
        epot_c( 131) = -0.167218462260d+05
        epot_c( 132) = -0.154061569920d+04
        epot_c( 133) = -0.198628866071d+05
        epot_c( 134) = -0.595045726927d+04
        epot_c( 135) = -0.614631655966d+04
        epot_c( 136) = -0.643860187595d+04
        epot_c( 137) = -0.328961702219d+05
        epot_c( 138) = -0.312220303325d+04
        epot_c( 139) = -0.117173908650d+05
        epot_c( 140) = 0.102108801243d+05
        epot_c( 141) = -0.296838514570d+05
        epot_c( 142) = 0.316743949442d+05
        epot_c( 143) = -0.547675942794d+04
        epot_c( 144) = -0.428760549879d+03
        epot_c( 145) = -0.281670316262d+03
        epot_c( 146) = -0.150670294917d+05
        epot_c( 147) = -0.251378269785d+05
        epot_c( 148) = -0.109642831245d+05
        epot_c( 149) = 0.452857325822d+04
        epot_c( 150) = 0.282458960852d+04
        epot_c( 151) = 0.177576149909d+04
        epot_c( 152) = 0.461479349791d+05
        epot_c( 153) = 0.677295737371d+04
        epot_c( 154) = -0.544615362371d+04
        epot_c( 155) = 0.797407925067d+04
        epot_c( 156) = 0.337166459111d+04
        epot_c( 157) = 0.531810770795d+04
        epot_c( 158) = 0.415106816445d+04
        epot_c( 159) = 0.842950842709d+04
        epot_c( 160) = 0.131227142369d+05
        epot_c( 161) = -0.278620017629d+05
        epot_c( 162) = -0.467327134275d+04
        epot_c( 163) = -0.274452920819d+04
        epot_c( 164) = -0.732936170936d+04
        epot_c( 165) = 0.128933968259d+05
        epot_c( 166) = 0.500270166783d+04
        epot_c( 167) = -0.832810100467d+04
        epot_c( 168) = 0.231068303833d+04
        epot_c( 169) = 0.736066804875d+04
        epot_c( 170) = -0.167133346329d+04
        epot_c( 171) = 0.271382328852d+05
        epot_c( 172) = -0.299226497133d+04
        epot_c( 173) = -0.300405606047d+05
        epot_c( 174) = -0.417234719382d+05
        epot_c( 175) = -0.925938352157d+04
        epot_c( 176) = 0.150738290565d+05
        epot_c( 177) = 0.621392381200d+04
        epot_c( 178) = -0.250466052381d+05
        epot_c( 179) = 0.189034139443d+05
        epot_c( 180) = 0.645847033559d+04
        epot_c( 181) = -0.122509319704d+04
        epot_c( 182) = 0.674497164564d+04
        epot_c( 183) = 0.165909951935d+04
        epot_c( 184) = 0.314085526715d+05
        epot_c( 185) = -0.567978068986d+04
        epot_c( 186) = 0.804344254374d+04
        epot_c( 187) = 0.718296953579d+04
        epot_c( 188) = 0.805081988188d+04
        epot_c( 189) = 0.959057200505d+04
        epot_c( 190) = 0.173175419540d+05
        epot_c( 191) = 0.192224951993d+04
        epot_c( 192) = -0.897479623343d+04
        epot_c( 193) = 0.489185253916d+04
        epot_c( 194) = -0.176425900434d+05
        epot_c( 195) = 0.373806080219d+04
        epot_c( 196) = -0.499566648122d+04
        epot_c( 197) = -0.453676607117d+04
        epot_c( 198) = -0.214019574689d+05
        epot_c( 199) = -0.760378970714d+04
        epot_c( 200) = -0.291704914069d+04
        epot_c( 201) = 0.237429310741d+05
        epot_c( 202) = -0.111950792650d+05
        epot_c( 203) = 0.293253496618d+05
        epot_c( 204) = 0.789023058401d+03
        epot_c( 205) = -0.424411435768d+04
        epot_c( 206) = -0.117087625656d+05
        epot_c( 207) = 0.694946236291d+04
        epot_c( 208) = -0.149281145977d+05
        epot_c( 209) = -0.136623142423d+04
        epot_c( 210) = -0.141052998243d+05
        epot_c( 211) = -0.102173374749d+05
        epot_c( 212) = -0.566906448414d+03
        epot_c( 213) = 0.789959020997d+04
        epot_c( 214) = -0.391936833293d+04
        epot_c( 215) = 0.999203207828d+04
        epot_c( 216) = 0.565812874635d+03
        epot_c( 217) = -0.882979001006d+04
        epot_c( 218) = 0.421173930005d+04
        epot_c( 219) = -0.625164802951d+04
        epot_c( 220) = 0.142249936194d+04
        epot_c( 221) = -0.122348255800d+05
        epot_c( 222) = 0.144730122892d+05
        epot_c( 223) = 0.716888617746d+02
        epot_c( 224) = -0.296205282215d+04
        epot_c( 225) = -0.135157525654d+05
        epot_c( 226) = 0.282282739575d+04
        epot_c( 227) = -0.427724924831d+04
        epot_c( 228) = 0.112055827751d+05
        epot_c( 229) = 0.570635552195d+04
        epot_c( 230) = -0.360280649223d+04
        epot_c( 231) = 0.102208442500d+05
        epot_c( 232) = 0.102712354122d+05
        epot_c( 233) = 0.261880233029d+04
        epot_c( 234) = 0.756346425878d+04
        epot_c( 235) = 0.164971200139d+05
        epot_c( 236) = 0.145688784480d+04
        epot_c( 237) = 0.296086198502d+04
        epot_c( 238) = -0.312536908014d+04
        epot_c( 239) = 0.121650939853d+05
        epot_c( 240) = -0.166454887214d+05
        epot_c( 241) = 0.812412831947d+02
        epot_c( 242) = -0.136497078414d+04
        epot_c( 243) = 0.870474131812d+04
        epot_c( 244) = -0.368363201039d+05
        epot_c( 245) = 0.255189655169d+05
        epot_c( 246) = 0.206043389737d+02
        epot_c( 247) = 0.488181964434d+03
        epot_c( 248) = 0.822653317496d+04
        epot_c( 249) = 0.115473149939d+04
        epot_c( 250) = 0.180422606051d+05
        epot_c( 251) = -0.187109086247d+04
        epot_c( 252) = 0.200831958953d+04
        epot_c( 253) = -0.144908090443d+05
        epot_c( 254) = 0.369931716826d+04
        epot_c( 255) = -0.251207665591d+04
        epot_c( 256) = -0.378840937076d+04
        epot_c( 257) = -0.251019035579d+04
        epot_c( 258) = -0.352218080889d+05
        epot_c( 259) = -0.157646059642d+04
        epot_c( 260) = 0.135390050658d+04
        epot_c( 261) = -0.673592136935d+04
        epot_c( 262) = 0.513499801117d+04
        epot_c( 263) = -0.597893982397d+04
        epot_c( 264) = -0.629318425103d+03
        epot_c( 265) = -0.327904556063d+04
        epot_c( 266) = -0.425819219862d+04
        epot_c( 267) = 0.138231691435d+05
        epot_c( 268) = 0.442022386425d+04
        epot_c( 269) = 0.149938258302d+04
        epot_c( 270) = 0.348504332442d+04
        epot_c( 271) = -0.561845437718d+04
        epot_c( 272) = -0.371238889983d+04
        epot_c( 273) = 0.266879264931d+04
        epot_c( 274) = -0.241086213543d+04
        epot_c( 275) = -0.256487559688d+04
        epot_c( 276) = 0.807623594819d+03
        epot_c( 277) = -0.146790064277d+05
        epot_c( 278) = 0.555226660757d+04
        epot_c( 279) = 0.235811102582d+05
        epot_c( 280) = -0.139884844872d+05
        epot_c( 281) = 0.329907188988d+05
        epot_c( 282) = -0.805451894503d+04
        epot_c( 283) = -0.103162828961d+05
        epot_c( 284) = 0.161264355953d+05
        epot_c( 285) = 0.114990797542d+05
        epot_c( 286) = -0.462324369583d+04
        epot_c( 287) = -0.368515824108d+04
        epot_c( 288) = -0.316340302073d+04
        epot_c( 289) = 0.515823486141d+03
        epot_c( 290) = -0.145836944077d+04
        epot_c( 291) = -0.244069295350d+04
        epot_c( 292) = 0.486854300243d+04
        epot_c( 293) = 0.840847301500d+04
        epot_c( 294) = 0.104526133121d+04
        epot_c( 295) = 0.884893182738d+03
        epot_c( 296) = -0.344498293927d+04
        epot_c( 297) = -0.208868965702d+05
        epot_c( 298) = -0.420935027569d+04
        epot_c( 299) = -0.838893564378d+04
        epot_c( 300) = 0.744963159083d+04
        epot_c( 301) = -0.542353889751d+04
        epot_c( 302) = -0.195893556763d+04
        epot_c( 303) = -0.466724683081d+04
        epot_c( 304) = 0.275184249484d+04
        epot_c( 305) = 0.277522069152d+04
        epot_c( 306) = -0.662767970556d+04
        epot_c( 307) = -0.319233977664d+03
        epot_c( 308) = -0.896417468066d+03
        epot_c( 309) = 0.123175057570d+05
        epot_c( 310) = -0.139755586021d+05
        epot_c( 311) = 0.506568047186d+04
        epot_c( 312) = -0.258952868793d+04
        epot_c( 313) = 0.264676924023d+04
        epot_c( 314) = -0.195443921477d+04
        epot_c( 315) = -0.127269946878d+05
        epot_c( 316) = 0.828536232661d+04
        epot_c( 317) = 0.214896348017d+04
        epot_c( 318) = -0.483068390482d+04
        epot_c( 319) = 0.565494771451d+04
        epot_c( 320) = 0.446561697944d+03
        epot_c( 321) = 0.364782524646d+04
        epot_c( 322) = 0.376605752474d+03
        epot_c( 323) = 0.579822188327d+03
        epot_c( 324) = -0.211839304989d+04
        epot_c( 325) = -0.472264011837d+04
        epot_c( 326) = 0.345025317650d+04
        epot_c( 327) = -0.676130963175d+01
        epot_c( 328) = 0.269359548788d+04
        epot_c( 329) = 0.276503090233d+04
        epot_c( 330) = 0.177885809699d+05
        epot_c( 331) = -0.586550985360d+03
        epot_c( 332) = -0.542475290546d+04
        epot_c( 333) = 0.509567352435d+04
        epot_c( 334) = 0.113613052156d+04
        epot_c( 335) = 0.170984112536d+03
        epot_c( 336) = 0.150487958530d+04
        epot_c( 337) = -0.753940462606d+03
        epot_c( 338) = 0.169631995550d+04
        epot_c( 339) = -0.121674596114d+05
        epot_c( 340) = 0.285415113429d+04
        epot_c( 341) = -0.922400820349d+04
        epot_c( 342) = 0.204864026100d+04
        epot_c( 343) = -0.871423146423d+03
        epot_c( 344) = 0.248567624983d+04
        epot_c( 345) = 0.614342185794d+04
        epot_c( 346) = -0.187662157284d+04
        epot_c( 347) = -0.857978013673d+04
        epot_c( 348) = -0.414781715613d+04
        epot_c( 349) = 0.357346302603d+04
        epot_c( 350) = 0.214164315004d+04
        epot_c( 351) = -0.590074931550d+04
        epot_c( 352) = 0.120405323962d+04
        epot_c( 353) = 0.645612000769d+03
        epot_c( 354) = 0.136516824296d+04
        epot_c( 355) = -0.169012829652d+04
        epot_c( 356) = 0.800622583324d+04
        epot_c( 357) = 0.992615757371d+03
        epot_c( 358) = 0.290405325075d+04
        epot_c( 359) = -0.598813578191d+03
        epot_c( 360) = -0.161195626593d+04
        epot_c( 361) = 0.545457796504d+03
        epot_c( 362) = -0.340800286247d+04
        epot_c( 363) = -0.377880997467d+03
        epot_c( 364) = 0.198426234996d+04
        epot_c( 365) = 0.705249584289d+04
        epot_c( 366) = 0.909789266216d+03
        epot_c( 367) = -0.119710848734d+04
        epot_c( 368) = 0.662721663144d+04
        epot_c( 369) = -0.388094987001d+04
        epot_c( 370) = -0.120538116588d+04
        epot_c( 371) = 0.121038056363d+04
        epot_c( 372) = -0.108996864000d+04
        epot_c( 373) = -0.554436730351d+04
        epot_c( 374) = 0.136234762955d+04
        epot_c( 375) = -0.474165011293d+04
        epot_c( 376) = -0.106383826047d+04
        epot_c( 377) = 0.212555444736d+04
        epot_c( 378) = -0.236013258835d+04
        epot_c( 379) = -0.301482092993d+04
        epot_c( 380) = 0.264793530210d+02
        epot_c( 381) = -0.205207581287d+04
        epot_c( 382) = -0.173304927521d+04
        epot_c( 383) = -0.581724015357d+03
        epot_c( 384) = 0.736106916388d+03
        epot_c( 385) = -0.862989104295d+03
        epot_c( 386) = -0.240480120333d+04
        epot_c( 387) = 0.344863969249d+04
        epot_c( 388) = 0.673930809805d+03
        epot_c( 389) = 0.538758861001d+03
        epot_c( 390) = 0.147967846147d+03
        epot_c( 391) = -0.189845771577d+04
        epot_c( 392) = -0.343381266242d+03
        epot_c( 393) = 0.226151971199d+04
        epot_c( 394) = -0.348165990383d+04
        epot_c( 395) = -0.137886410095d+04
        epot_c( 396) = -0.452426797538d+04
        epot_c( 397) = -0.546311644413d+02
        epot_c( 398) = -0.480113520257d+02
        epot_c( 399) = -0.260536444352d+03
        epot_c( 400) = -0.133327190832d+04
        epot_c( 401) = 0.632765868126d+03
        epot_c( 402) = -0.661117799374d+03
        epot_c( 403) = 0.172506945261d+04
        epot_c( 404) = -0.615072273819d+03
        epot_c( 405) = -0.383815665429d+03
        epot_c( 406) = 0.169025584445d+04
        epot_c( 407) = -0.808965608264d+03
        epot_c( 408) = 0.593728899296d+03
        epot_c( 409) = 0.119569641767d+04
        epot_c( 410) = 0.430656118447d+03
        epot_c( 411) = 0.587203407297d+04
        epot_c( 412) = -0.773418173876d+03
        epot_c( 413) = 0.204092911587d+03
        epot_c( 414) = 0.198285065744d+04
        epot_c( 415) = -0.145843638158d+04
        epot_c( 416) = 0.154180304340d+04
        epot_c( 417) = -0.113615015119d+03
        epot_c( 418) = 0.341159952346d+03
        epot_c( 419) = 0.128014772782d+03
        epot_c( 420) = -0.241151083483d+04
        epot_c( 421) = -0.888711069985d+03
        epot_c( 422) = -0.255786858400d+03
        epot_c( 423) = -0.516753605238d+03
        epot_c( 424) = 0.952501144399d+03
        epot_c( 425) = 0.734163721480d+03
        epot_c( 426) = -0.103271529318d+04
        epot_c( 427) = 0.212633590214d+03
        epot_c( 428) = 0.303411042713d+03
        epot_c( 429) = 0.409859079006d+03
        epot_c( 430) = -0.144512479307d+03
        v            = 0.d0; g = 0.d0; d = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                j     = 3*(iatom-1)+idir
                tx(j) = x(iatom,idir)
            enddo
        enddo
        tg = 0.d0
        call o4pes(tx, v, tg, igrad)
        v  = v/23.0609
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
        !> @brief Evaluate o4pes for the O4_triplet_v2 potential.
        subroutine o4pes( x, v, dvdx, igrad )
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
            !      common /epot/     rMs(6),rM(0:111),P(0:305),C(430),B(430)
            !      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
            !     $                  dRdX(6,12),dBdR(6,276)
            !      common /msprmt/   a,ab,ra,rb
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
        end subroutine o4pes
        !> @brief Evaluate coord_convt for the O4_triplet_v2 potential.
        subroutine coord_convt( x )
            !  Program to calculate the six interatomic distance
            !  by reading XYZ coordinate
            implicit none
            integer          :: i
            double precision :: x(12)
            !  Now, calculate the interatomic distances
            !  r1 = r(O1O2)    r2 = r(O1O3)
            !  r3 = r(O1O4)    r4 = r(O2O3)
            !  r5 = r(O2O4)    r6 = r(O3O4)
            coord_r(1) = sqrt((x(4)-x(1))**2 + (x(5)-x(2))**2 + (x(6)-x(3))**2)
            coord_r(2) = sqrt((x(7)-x(1))**2 + (x(8)-x(2))**2 + (x(9)-x(3))**2)
            coord_r(3) = sqrt((x(10)-x(1))**2 + (x(11)-x(2))**2 + (x(12)-x(3))**2)
            coord_r(4) = sqrt((x(4)-x(7))**2 + (x(5)-x(8))**2 + (x(6)-x(9))**2)
            coord_r(5) = sqrt((x(4)-x(10))**2 + (x(5)-x(11))**2 + (x(6)-x(12))**2)
            coord_r(6) = sqrt((x(7)-x(10))**2 + (x(8)-x(11))**2 + (x(9)-x(12))**2)
        end subroutine coord_convt
        !> @brief Evaluate EvV for the O4_triplet_v2 potential.
        subroutine evv( v )
            ! Subroutine to evaluate V for given R
            ! V(R) = C*P
            ! C:            Coefficients, stored in 'dim.inc'
            ! P:            Basis functions evaluated for given R
            ! rMs:          rMs(6), six mixed exponetial-Gaussian terms (MEG)
            ! a:            Nonlinear parameters in Morse terms (Angstrom)
            ! ab:           Nonlinear parameters in Gauss terms (Angstrom^2)
            ! re:           Equilibrium bond length(Angstrom)
            ! nop:          number of points
            ! nom:          number of monomials
            ! nob:          number of basis functions (polynomials)
            ! rM(0:111):    Array to store monomials
            ! P(0:305):     Array to store polynomials
            ! B(1:276):     Array to store basis functions
            implicit none
            integer          :: i, j, k
            double precision :: dist, dv2dr, v, v2, disp, dispdr(6)
            ! Calculate the six MEG terms for each point
            call evmorse
            ! Calculate the monomials for each point by using six MEG terms
            call evmono
            ! Calculate the polynomials (basis functions) by using monomials
            call evpoly
            ! Calculate the basis functions by removing unconnected and 2-body terms
            call evbas
            ! Initialized v to be 2De 2*120.243 kcal/mol
            v = 240.486d0
            ! Evaluate 2-body interactions
            do i=1, 6
                dist = coord_r(i)
                call ev2gm2(dist, v2, dv2dr, 4, 0)
                v = v+v2
            enddo
            ! Add D3 dispersion correction
            call d3disp(coord_r, disp, dispdr, 0)
            v = v+disp
            ! Evaluate V by taken the product of C and Basis function array
            do i=1, 430
                v = v + epot_c(i)*epot_b(i)
            enddo
            !      Write(*,9999) V
            ! 9999 Format('The potential energy is ',F20.14,' kcal/mol')
        end subroutine evv
        !> @brief Evaluate EvdVdX for the O4_triplet_v2 potential.
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
            ! 6
            do i=1, 12
                do j=1, 6
                    dvdx(i) = dvdx(i) + gradt_dvdr(j)*gradt_drdx(j,i)
                enddo
            enddo
            !      write(*,*) 'The 12 dVdX are:'
            !      Write(*,9999) (dVdX(i),i=1,12)
            ! 9999 Format(1x,3F15.8)
        end subroutine evdvdx
        !> @brief Evaluate EvMorse for the O4_triplet_v2 potential.
        subroutine evmorse
            implicit none
            integer :: i
            ! mixed exponential-Gaussian terms = exp(-(r-re)/a-(r-re)^2/ab)
            ! re:   equlibrium bond length(1.208A)
            ! a:    nonlinear parameter optimized for diatomic molecules, unit
            ! Anstrom
            ! ab:   nonlinear parameter optimized for diatomic molecules, unit
            ! Anstrom^2
            do i=1, 6
                epot_rms(i) = exp(-(coord_r(i)-msprmt_ra)/msprmt_a-((coord_r(i)-msprmt_rb)**2.0d0)/msprmt_ab)
            enddo
        end subroutine evmorse
        !> @brief Evaluate EvMono for the O4_triplet_v2 potential.
        subroutine evmono
            !  The subroutine reads six MEG variables(X) and calculates the
            !  monomials (M), the polynomial is built from the monomials.
            !  For A4 with max. degree 10, the number of monomials is nom.
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
        !> @brief Evaluate EvPoly for the O4_triplet_v2 potential.
        subroutine evpoly
            !  The subroutine reads monomials (m) and calculates the
            !  permutation invariant polynomials (p).
            !  For A4 with max. degree 10, the number of polynomials is nob.
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
            epot_p(306) = epot_p(40)*epot_p(11)
            epot_p(307) = epot_p(40)*epot_p(12)
            epot_p(308) = epot_p(40)*epot_p(17)
            epot_p(309) = epot_p(40)*epot_p(13)
            epot_p(310) = epot_p(40)*epot_p(14)
            epot_p(311) = epot_p(40)*epot_p(15)
            epot_p(312) = epot_p(40)*epot_p(16)
            epot_p(313) = epot_p(40)*epot_p(18)
            epot_p(314) = epot_p(40)*epot_p(19)
            epot_p(315) = epot_p(40)*epot_p(20)
            epot_p(316) = epot_p(40)*epot_p(21)
            epot_p(317) = epot_p(11)*epot_p(42) - epot_p(310)
            epot_p(318) = epot_p(11)*epot_p(44) - epot_p(307)
            epot_p(319) = epot_p(12)*epot_p(43) - epot_p(309) - epot_p(318)
            epot_p(320) = epot_p(17)*epot_p(42) - epot_p(314)
            epot_p(321) = epot_p(2)*epot_p(125) - epot_p(311)
            epot_p(322) = epot_p(2)*epot_p(126) - epot_p(312)
            epot_p(323) = epot_p(11)*epot_p(48) - epot_p(313)
            epot_p(324) = epot_p(12)*epot_p(42) - epot_p(311) - epot_p(312) - epot_p(307) - epot_p(307)
            epot_p(325) = epot_p(6)*epot_p(79) - epot_p(319)
            epot_p(326) = epot_p(11)*epot_p(54)
            epot_p(327) = epot_p(7)*epot_p(79) - epot_p(319)
            epot_p(328) = epot_p(2)*epot_p(131) - epot_p(315) - epot_p(324)
            epot_p(329) = epot_p(22)*epot_p(29) - epot_p(313) - epot_p(311) - epot_p(326) - epot_p(311)
            epot_p(330) = epot_p(11)*epot_p(55)
            epot_p(331) = epot_p(22)*epot_p(30) - epot_p(313) - epot_p(312) - epot_p(330) - epot_p(312)
            epot_p(332) = epot_p(2)*epot_p(127) - epot_p(309) - epot_p(319)
            epot_p(333) = epot_p(6)*epot_p(83) - epot_p(318)
            epot_p(334) = epot_p(7)*epot_p(83) - epot_p(318)
            epot_p(335) = epot_p(11)*epot_p(52) - epot_p(315) - epot_p(324)
            epot_p(336) = epot_p(2)*epot_p(130) - epot_p(313) - epot_p(323) - epot_p(313)
            epot_p(337) = epot_p(2)*epot_p(133) - epot_p(316)
            epot_p(338) = epot_p(6)*epot_p(77) - epot_p(315) - epot_p(309) - epot_p(322)
            epot_p(339) = epot_p(6)*epot_p(78) - epot_p(312)
            epot_p(340) = epot_p(7)*epot_p(77) - epot_p(315) - epot_p(309) - epot_p(321)
            epot_p(341) = epot_p(7)*epot_p(78) - epot_p(311)
            epot_p(342) = epot_p(11)*epot_p(59) - epot_p(329) - epot_p(338)
            epot_p(343) = epot_p(11)*epot_p(60) - epot_p(331) - epot_p(340)
            epot_p(344) = epot_p(3)*epot_p(130) - epot_p(315) - epot_p(311) - epot_p(312) - epot_p(309) - &
                & epot_p(329) - epot_p(338) - epot_p(331) - epot_p(340) - epot_p(328) - epot_p(321) - epot_p(322) - &
                & epot_p(311) - epot_p(312)
            epot_p(345) = epot_p(18)*epot_p(42) - epot_p(310) - epot_p(326) - epot_p(330) - epot_p(310)
            epot_p(346) = epot_p(2)*epot_p(132) - epot_p(315) - epot_p(335)
            epot_p(347) = epot_p(6)*epot_p(92) - epot_p(334)
            epot_p(348) = epot_p(7)*epot_p(92) - epot_p(333)
            epot_p(349) = epot_p(3)*epot_p(132) - epot_p(316) - epot_p(315) - epot_p(314) - epot_p(338) - &
                & epot_p(347) - epot_p(340) - epot_p(348) - epot_p(336) - epot_p(333) - epot_p(334) - epot_p(316) - &
                & epot_p(315) - epot_p(314) - epot_p(316) - epot_p(314) - epot_p(316) - epot_p(314)
            epot_p(350) = epot_p(3)*epot_p(133) - epot_p(315) - epot_p(339) - epot_p(341)
            epot_p(351) = epot_p(4)*epot_p(132) - epot_p(315) - epot_p(344) - epot_p(342) - epot_p(343) - epot_p(332)
            epot_p(352) = epot_p(1)*epot_p(222) - epot_p(316) - epot_p(350)
            epot_p(353) = epot_p(2)*epot_p(134) - epot_p(306)
            epot_p(354) = epot_p(2)*epot_p(135) - epot_p(323)
            epot_p(355) = epot_p(3)*epot_p(134) - epot_p(319)
            epot_p(356) = epot_p(11)*epot_p(53) - epot_p(310) - epot_p(320)
            epot_p(357) = epot_p(2)*epot_p(144) - epot_p(329) - epot_p(338)
            epot_p(358) = epot_p(2)*epot_p(143) - epot_p(324) - epot_p(335)
            epot_p(359) = epot_p(6)*epot_p(81) - epot_p(311) - epot_p(324)
            epot_p(360) = epot_p(2)*epot_p(146) - epot_p(331) - epot_p(340)
            epot_p(361) = epot_p(2)*epot_p(150) - epot_p(344)
            epot_p(362) = epot_p(7)*epot_p(82) - epot_p(312) - epot_p(324)
            epot_p(363) = epot_p(11)*epot_p(65) - epot_p(308)
            epot_p(364) = epot_p(2)*epot_p(137) - epot_p(318)
            epot_p(365) = epot_p(17)*epot_p(45) - epot_p(307)
            epot_p(366) = epot_p(4)*epot_p(134) - epot_p(317)
            epot_p(367) = epot_p(6)*epot_p(96) - epot_p(332)
            epot_p(368) = epot_p(17)*epot_p(54)
            epot_p(369) = epot_p(7)*epot_p(96) - epot_p(332)
            epot_p(370) = epot_p(17)*epot_p(55)
            epot_p(371) = epot_p(2)*epot_p(159) - epot_p(345) - epot_p(349)
            epot_p(372) = epot_p(2)*epot_p(147) - epot_p(313)
            epot_p(373) = epot_p(11)*epot_p(58) - epot_p(315) - epot_p(328)
            epot_p(374) = epot_p(2)*epot_p(148) - epot_p(329) - epot_p(342)
            epot_p(375) = epot_p(6)*epot_p(98) - epot_p(327)
            epot_p(376) = epot_p(2)*epot_p(166) - epot_p(359)
            epot_p(377) = epot_p(14)*epot_p(54) - epot_p(323)
            epot_p(378) = epot_p(2)*epot_p(149) - epot_p(331) - epot_p(343)
            epot_p(379) = epot_p(7)*epot_p(98) - epot_p(325)
            epot_p(380) = epot_p(7)*epot_p(99) - epot_p(311) - epot_p(359)
            epot_p(381) = epot_p(2)*epot_p(167) - epot_p(362)
            epot_p(382) = epot_p(14)*epot_p(55) - epot_p(323)
            epot_p(383) = epot_p(6)*epot_p(100) - epot_p(312) - epot_p(362)
            epot_p(384) = epot_p(11)*epot_p(66) - epot_p(344)
            epot_p(385) = epot_p(6)*epot_p(101) - epot_p(325) - epot_p(343) - epot_p(379)
            epot_p(386) = epot_p(7)*epot_p(101) - epot_p(342) - epot_p(327) - epot_p(375)
            epot_p(387) = epot_p(6)*epot_p(103) - epot_p(313) - epot_p(382)
            epot_p(388) = epot_p(11)*epot_p(67) - epot_p(314)
            epot_p(389) = epot_p(2)*epot_p(152) - epot_p(333)
            epot_p(390) = epot_p(2)*epot_p(153) - epot_p(334)
            epot_p(391) = epot_p(2)*epot_p(154) - epot_p(315) - epot_p(373)
            epot_p(392) = epot_p(1)*epot_p(240) - epot_p(335) - epot_p(373) - epot_p(366)
            epot_p(393) = epot_p(2)*epot_p(155) - epot_p(338) - epot_p(342)
            epot_p(394) = epot_p(2)*epot_p(171) - epot_p(377)
            epot_p(395) = epot_p(2)*epot_p(157) - epot_p(340) - epot_p(343)
            epot_p(396) = epot_p(2)*epot_p(163) - epot_p(350) - epot_p(351)
            epot_p(397) = epot_p(6)*epot_p(95) - epot_p(315) - epot_p(350) - epot_p(340)
            epot_p(398) = epot_p(2)*epot_p(172) - epot_p(382)
            epot_p(399) = epot_p(7)*epot_p(95) - epot_p(315) - epot_p(350) - epot_p(338)
            epot_p(400) = epot_p(11)*epot_p(68) - epot_p(345) - epot_p(349)
            epot_p(401) = epot_p(2)*epot_p(175) - epot_p(380) - epot_p(385)
            epot_p(402) = epot_p(6)*epot_p(106) - epot_p(335) - epot_p(343) - epot_p(396)
            epot_p(403) = epot_p(2)*epot_p(176) - epot_p(383) - epot_p(386)
            epot_p(404) = epot_p(7)*epot_p(106) - epot_p(342) - epot_p(335) - epot_p(396)
            epot_p(405) = epot_p(4)*epot_p(150) - epot_p(328) - epot_p(359) - epot_p(362)
            epot_p(406) = epot_p(11)*epot_p(69) - epot_p(316)
            epot_p(407) = epot_p(2)*epot_p(161) - epot_p(347)
            epot_p(408) = epot_p(2)*epot_p(162) - epot_p(348)
            epot_p(409) = epot_p(11)*epot_p(70) - epot_p(350) - epot_p(351)
            epot_p(410) = epot_p(2)*epot_p(180) - epot_p(397) - epot_p(402)
            epot_p(411) = epot_p(6)*epot_p(110) - epot_p(348)
            epot_p(412) = epot_p(2)*epot_p(181) - epot_p(399) - epot_p(404)
            epot_p(413) = epot_p(7)*epot_p(110) - epot_p(347)
            epot_p(414) = epot_p(6)*epot_p(112) - epot_p(316) - epot_p(398) - epot_p(316)
            epot_p(415) = epot_p(11)*epot_p(71) - epot_p(352)
            epot_p(416) = epot_p(2)*epot_p(184) - epot_p(411)
            epot_p(417) = epot_p(2)*epot_p(185) - epot_p(413)
            epot_p(418) = epot_p(12)*epot_p(71) - epot_p(351) - epot_p(416) - epot_p(417)
            epot_p(419) = epot_p(2)*epot_p(164) - epot_p(320) - epot_p(356)
            epot_p(420) = epot_p(2)*epot_p(165) - epot_p(328) - epot_p(373)
            epot_p(421) = epot_p(2)*epot_p(170) - epot_p(337) - epot_p(392)
            epot_p(422) = epot_p(1)*epot_p(268) - epot_p(359)
            epot_p(423) = epot_p(1)*epot_p(269) - epot_p(362)
            epot_p(424) = epot_p(2)*epot_p(174) - epot_p(345) - epot_p(400)
            epot_p(425) = epot_p(6)*epot_p(102) - epot_p(345) - epot_p(326)
            epot_p(426) = epot_p(7)*epot_p(103) - epot_p(345) - epot_p(330)
            epot_p(427) = epot_p(3)*epot_p(186) - epot_p(364)
            epot_p(428) = epot_p(18)*epot_p(65) - epot_p(354)
            epot_p(429) = epot_p(17)*epot_p(66) - epot_p(361)
            epot_p(430) = epot_p(2)*epot_p(179) - epot_p(350) - epot_p(409)
            epot_p(431) = epot_p(4)*epot_p(166) - epot_p(361) - epot_p(357) - epot_p(422)
            epot_p(432) = epot_p(4)*epot_p(167) - epot_p(361) - epot_p(360) - epot_p(423)
            epot_p(433) = epot_p(2)*epot_p(187) - epot_p(387)
            epot_p(434) = epot_p(14)*epot_p(66) - epot_p(328) - epot_p(405) - epot_p(376) - epot_p(381) - epot_p(373)
            epot_p(435) = epot_p(1)*epot_p(277) - epot_p(387) - epot_p(385) - epot_p(425)
            epot_p(436) = epot_p(1)*epot_p(278) - epot_p(387) - epot_p(386) - epot_p(426)
            epot_p(437) = epot_p(2)*epot_p(177) - epot_p(346) - epot_p(391)
            epot_p(438) = epot_p(2)*epot_p(178) - epot_p(349) - epot_p(400)
            epot_p(439) = epot_p(1)*epot_p(281) - epot_p(396) - epot_p(392) - epot_p(430)
            epot_p(440) = epot_p(21)*epot_p(54) - epot_p(370)
            epot_p(441) = epot_p(21)*epot_p(55) - epot_p(368)
            epot_p(442) = epot_p(2)*epot_p(189) - epot_p(405) - epot_p(434)
            epot_p(443) = epot_p(1)*epot_p(285) - epot_p(402) - epot_p(404) - epot_p(400) - epot_p(434) - epot_p(430)
            epot_p(444) = epot_p(6)*epot_p(116) - epot_p(347) - epot_p(413) - epot_p(403)
            epot_p(445) = epot_p(7)*epot_p(116) - epot_p(348) - epot_p(411) - epot_p(401)
            epot_p(446) = epot_p(2)*epot_p(182) - epot_p(351) - epot_p(409)
            epot_p(447) = epot_p(2)*epot_p(191) - epot_p(414) - epot_p(443)
            epot_p(448) = epot_p(3)*epot_p(183) - epot_p(351) - epot_p(418) - epot_p(411) - epot_p(413) - epot_p(409)
            epot_p(449) = epot_p(6)*epot_p(118) - epot_p(351) - epot_p(418) - epot_p(412)
            epot_p(450) = epot_p(7)*epot_p(118) - epot_p(351) - epot_p(418) - epot_p(410)
            epot_p(451) = epot_p(2)*epot_p(193) - epot_p(418) - epot_p(448)
            epot_p(452) = epot_p(1)*epot_p(294) - epot_p(418) - epot_p(415) - epot_p(448)
            epot_p(453) = epot_p(6)*epot_p(119) - epot_p(417)
            epot_p(454) = epot_p(7)*epot_p(119) - epot_p(416)
            epot_p(455) = epot_p(2)*epot_p(186) - epot_p(363)
            epot_p(456) = epot_p(3)*epot_p(187) - epot_p(385) - epot_p(386) - epot_p(384) - epot_p(435) - &
                & epot_p(436) - epot_p(434)
            epot_p(457) = epot_p(2)*epot_p(188) - epot_p(388)
            epot_p(458) = epot_p(4)*epot_p(187) - epot_p(425) - epot_p(426) - epot_p(424)
            epot_p(459) = epot_p(2)*epot_p(190) - epot_p(406)
            epot_p(460) = epot_p(21)*epot_p(66) - epot_p(422) - epot_p(423) - epot_p(420)
            epot_p(461) = epot_p(2)*epot_p(192) - epot_p(415)
            epot_p(462) = epot_p(18)*epot_p(71) - epot_p(440) - epot_p(441) - epot_p(438)
            epot_p(463) = epot_p(2)*epot_p(194) - epot_p(452)
            epot_p(464) = epot_p(3)*epot_p(194) - epot_p(453) - epot_p(454) - epot_p(451)
            epot_p(465) = epot_p(1)*epot_p(305) - epot_p(464) - epot_p(463)
        end subroutine evpoly
        !> @brief Evaluate evbas for the O4_triplet_v2 potential.
        subroutine evbas
            !  This subroutine eliminates the 2-body terms in Bowman's approach.
            implicit none
            integer          :: i
            double precision :: b1(466)
            ! Pass P(0:465) to BM1(1:466)
            do i=1, 466
                b1(i) = epot_p(i-1)
            enddo
            ! Remove unconnected terms and 2-body terms and pass to B(1:430)
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
            do i=277, 425
                epot_b(i) = b1(i+30)
            enddo
            epot_b(426) = b1(457)
            epot_b(427) = b1(459)
            epot_b(428) = b1(461)
            epot_b(429) = b1(463)
            epot_b(430) = b1(465)
        end subroutine evbas
        !> @brief Evaluate ev2gm2 for the O4_triplet_v2 potential.
        subroutine ev2gm2( r, v, grad, imol, igrad )
            !
            ! Compute the diatomic potential of ground-state triplet O2
            !
            ! References: J. Chem. Phys. 132, 074307 (2010)
            !
            ! Input:  r      interatomic distance in Angstrom
            ! Output: V      potential in kcal/mol
            !         grad   gradient (kcal/mol)/Angstrom
            !
            implicit none
            integer          :: imol
            double precision :: r
            double precision :: v, grad
            ! Parameters of analytical even-tempered Gaussian expansions for the
            ! ground
            ! state potential energy curve of O2 CBS+SR+SO+CV. Units: alpha in
            ! Anstromgs^-2, beta=dimensionless, a_k in milihartree.
            double precision :: alpha, beta, a(0:7)
            integer          :: k, igrad
            ! Original parameters
            !      alpha = 0.785d0
            !      beta = 1.307d0
            !      a(0) = -2388.5641690d0
            !      a(1) = 18086.977116d0
            !      a(2) = -71760.197585d0
            !      a(3) = 154738.09175d0
            !      a(4) = -215074.85646d0
            !      a(5) = 214799.54567d0
            !      a(6) = -148395.42850d0
            !      a(7) = 73310.781453d0
            ! Modified parameters for D3(BJ)
            alpha = 9.439784362354936d-1
            beta  = 1.262242998506810d0
            a(0)  = -1.488979427684798d3
            a(1)  = 1.881435846488955d4
            a(2)  = -1.053475425838226d5
            a(3)  = 2.755135591229064d5
            a(4)  = -4.277588997761775d5
            a(5)  = 4.404104009614092d5
            a(6)  = -2.946204062950765d5
            a(7)  = 1.176861219078620d5
            v     = 0.d0
            do k=0, 7
                v = v + a(k)*dexp(-alpha*beta**k*r**2)
            enddo
            ! From milihartree to kcal/mol
            v = v*627.509523475149d-3
            ! Compute the gradient if needed
            if( igrad .eq. 1 )then
                grad = 0.d0
                do k=0, 7
                    grad = grad-2.d0*a(k)*alpha*beta**k*r*dexp(-alpha*beta**k*r**2)
                enddo
                ! Convert from milihartree/A to i(kcal/mol)/A
                grad = grad*627.509523475149d-3
            endif
        end subroutine ev2gm2
        !> @brief Evaluate EvdVdR for the O4_triplet_v2 potential.
        subroutine evdvdr
            ! Subroutine to evaluate dVdR for given R
            ! dVdR = dV2dR + C*dBdR
            ! C:            Coefficients, stored in 'dim.inc'
            ! P:            Basis functions evaluated for given R
            ! M:            Monomials evaluated for given R
            ! dV2dR:        Gradient of 2-body interactions
            ! dMsdR:        dMsdR(6,6), 6 MEG terms w.r.t. 6 bond lengths
            ! dMdR:         dMdR(6,nom), nom monomials w.r.t. 6 bond lengths
            ! dPdR:         dPdR(6,nob), nop polynomial basis functions
            !               w.r.t. 6 bond length
            ! nom:          number of monomials
            ! nob:          number of basis functions (polynomials)
            ! M(nom):       Array to store monomials
            ! P(nob):       Array to store polynomials
            implicit none
            integer          :: i, j
            double precision :: dist, v2, dv2dr, disp, dispdr(6)
            ! Initialize dVdR(6)
            do i=1, 6
                gradt_dvdr(i) = 0.0d0
            enddo
            ! Add dV2dR(i) to dVdR
            do i=1, 6
                dist = coord_r(i)
                call ev2gm2(dist, v2, dv2dr, 4, 1)
                gradt_dvdr(i) = dv2dr
            enddo
            ! Add numerical gradient of D3 dispersion correction
            call d3disp(coord_r, disp, dispdr, 1)
            do i=1, 6
                gradt_dvdr(i) = gradt_dvdr(i) + dispdr(i)
            enddo
            ! Calculate dMEG/dr(6,6) for given R(6)
            call evdmsdr
            ! Calculate the monomials for each point by using six MEG terms
            call evdmdr
            ! Calculate the polynomials by using monomials
            call evdpdr
            ! Remove 2-body interactions and unconnected terms from polynomials
            call evdbdr
            ! Evaluate dVdR(6) by taken the product of C(j) and dPdR(i,j)
            do i=1, 6
                do j=1, 430
                    gradt_dvdr(i) = gradt_dvdr(i) + epot_c(j)*gradt_dbdr(i,j)
                enddo
            enddo
        end subroutine evdvdr
        !> @brief Evaluate EvdRdX for the O4_triplet_v2 potential.
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
        !> @brief Evaluate EvdMsdR for the O4_triplet_v2 potential.
        subroutine evdmsdr
            ! Subroutine to evaluate the derivatives of MEG term X
            ! w.r.t. interatomic distance R(6)
            ! dmsdR:        Local variables, dirm(6,6)
            ! a:            Nonlinear pamameter(Angstrom)
            ! ab:       Nonlinear pamameter(Angstrom^2)
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
            ! MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
            ! dmsdr(i,j)=0  i!=j
            !
            do i=1, 6
                gradt_dmsdr(i,i) = (-2.0d0*(coord_r(i)-msprmt_rb)/msprmt_ab-1/msprmt_a)* exp(-(coord_r(i)-msprmt_ra)/ &
                    & msprmt_a-((coord_r(i)-msprmt_rb)**2.0d0)/msprmt_ab)
            enddo
        end subroutine evdmsdr
        !> @brief Evaluate EvdMdR for the O4_triplet_v2 potential.
        subroutine evdmdr
            !  The subroutine reads M(nom) and dMSdR(6,6) and calculates the
            !  dMdR(6,nom), the polynomial is built from the monomials.
            !  For A4 with max. degree 10, the number of monomials is nom.
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
        !> @brief Evaluate EvdPdr for the O4_triplet_v2 potential.
        subroutine evdpdr
            !  The subroutine reads monomials (m) and calculates the
            !  permutation invariant polynomials (p)
            !  For A4 with max. degree 10, the number of polynomials is nob.
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
                gradt_dpdr(i,306) = gradt_dpdr(i,40)*epot_p(11) + epot_p(40)*gradt_dpdr(i,11)
                gradt_dpdr(i,307) = gradt_dpdr(i,40)*epot_p(12) + epot_p(40)*gradt_dpdr(i,12)
                gradt_dpdr(i,308) = gradt_dpdr(i,40)*epot_p(17) + epot_p(40)*gradt_dpdr(i,17)
                gradt_dpdr(i,309) = gradt_dpdr(i,40)*epot_p(13) + epot_p(40)*gradt_dpdr(i,13)
                gradt_dpdr(i,310) = gradt_dpdr(i,40)*epot_p(14) + epot_p(40)*gradt_dpdr(i,14)
                gradt_dpdr(i,311) = gradt_dpdr(i,40)*epot_p(15) + epot_p(40)*gradt_dpdr(i,15)
                gradt_dpdr(i,312) = gradt_dpdr(i,40)*epot_p(16) + epot_p(40)*gradt_dpdr(i,16)
                gradt_dpdr(i,313) = gradt_dpdr(i,40)*epot_p(18) + epot_p(40)*gradt_dpdr(i,18)
                gradt_dpdr(i,314) = gradt_dpdr(i,40)*epot_p(19) + epot_p(40)*gradt_dpdr(i,19)
                gradt_dpdr(i,315) = gradt_dpdr(i,40)*epot_p(20) + epot_p(40)*gradt_dpdr(i,20)
                gradt_dpdr(i,316) = gradt_dpdr(i,40)*epot_p(21) + epot_p(40)*gradt_dpdr(i,21)
                gradt_dpdr(i,317) = gradt_dpdr(i,11)*epot_p(42) + epot_p(11)*gradt_dpdr(i,42) - gradt_dpdr(i,310)
                gradt_dpdr(i,318) = gradt_dpdr(i,11)*epot_p(44) + epot_p(11)*gradt_dpdr(i,44) - gradt_dpdr(i,307)
                gradt_dpdr(i,319) = gradt_dpdr(i,12)*epot_p(43) + epot_p(12)*gradt_dpdr(i,43) - gradt_dpdr(i,309) - &
                    & gradt_dpdr(i,318)
                gradt_dpdr(i,320) = gradt_dpdr(i,17)*epot_p(42) + epot_p(17)*gradt_dpdr(i,42) - gradt_dpdr(i,314)
                gradt_dpdr(i,321) = gradt_dpdr(i,2)*epot_p(125) + epot_p(2)*gradt_dpdr(i,125) - gradt_dpdr(i,311)
                gradt_dpdr(i,322) = gradt_dpdr(i,2)*epot_p(126) + epot_p(2)*gradt_dpdr(i,126) - gradt_dpdr(i,312)
                gradt_dpdr(i,323) = gradt_dpdr(i,11)*epot_p(48) + epot_p(11)*gradt_dpdr(i,48) - gradt_dpdr(i,313)
                gradt_dpdr(i,324) = gradt_dpdr(i,12)*epot_p(42) + epot_p(12)*gradt_dpdr(i,42) - gradt_dpdr(i,311) - &
                    & gradt_dpdr(i,312) - gradt_dpdr(i,307) - gradt_dpdr(i,307)
                gradt_dpdr(i,325) = gradt_dpdr(i,6)*epot_p(79)  + epot_p(6)*gradt_dpdr(i,79) - gradt_dpdr(i,319)
                gradt_dpdr(i,326) = gradt_dpdr(i,11)*epot_p(54) + epot_p(11)*gradt_dpdr(i,54)
                gradt_dpdr(i,327) = gradt_dpdr(i,7)*epot_p(79)  + epot_p(7)*gradt_dpdr(i,79) - gradt_dpdr(i,319)
                gradt_dpdr(i,328) = gradt_dpdr(i,2)*epot_p(131) + epot_p(2)*gradt_dpdr(i,131) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,324)
                gradt_dpdr(i,329) = gradt_dpdr(i,22)*epot_p(29) + epot_p(22)*gradt_dpdr(i,29) - gradt_dpdr(i,313) - &
                    & gradt_dpdr(i,311) - gradt_dpdr(i,326) - gradt_dpdr(i,311)
                gradt_dpdr(i,330) = gradt_dpdr(i,11)*epot_p(55) + epot_p(11)*gradt_dpdr(i,55)
                gradt_dpdr(i,331) = gradt_dpdr(i,22)*epot_p(30) + epot_p(22)*gradt_dpdr(i,30) - gradt_dpdr(i,313) - &
                    & gradt_dpdr(i,312) - gradt_dpdr(i,330) - gradt_dpdr(i,312)
                gradt_dpdr(i,332) = gradt_dpdr(i,2)*epot_p(127) + epot_p(2)*gradt_dpdr(i,127) - gradt_dpdr(i,309) - &
                    & gradt_dpdr(i,319)
                gradt_dpdr(i,333) = gradt_dpdr(i,6)*epot_p(83)  + epot_p(6)*gradt_dpdr(i,83) - gradt_dpdr(i,318)
                gradt_dpdr(i,334) = gradt_dpdr(i,7)*epot_p(83)  + epot_p(7)*gradt_dpdr(i,83) - gradt_dpdr(i,318)
                gradt_dpdr(i,335) = gradt_dpdr(i,11)*epot_p(52) + epot_p(11)*gradt_dpdr(i,52) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,324)
                gradt_dpdr(i,336) = gradt_dpdr(i,2)*epot_p(130) + epot_p(2)*gradt_dpdr(i,130) - gradt_dpdr(i,313) - &
                    & gradt_dpdr(i,323) - gradt_dpdr(i,313)
                gradt_dpdr(i,337) = gradt_dpdr(i,2)*epot_p(133) + epot_p(2)*gradt_dpdr(i,133) - gradt_dpdr(i,316)
                gradt_dpdr(i,338) = gradt_dpdr(i,6)*epot_p(77)  + epot_p(6)*gradt_dpdr(i,77)  - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,309) - gradt_dpdr(i,322)
                gradt_dpdr(i,339) = gradt_dpdr(i,6)*epot_p(78)  + epot_p(6)*gradt_dpdr(i,78) - gradt_dpdr(i,312)
                gradt_dpdr(i,340) = gradt_dpdr(i,7)*epot_p(77)  + epot_p(7)*gradt_dpdr(i,77) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,309) - gradt_dpdr(i,321)
                gradt_dpdr(i,341) = gradt_dpdr(i,7)*epot_p(78)  + epot_p(7)*gradt_dpdr(i,78) - gradt_dpdr(i,311)
                gradt_dpdr(i,342) = gradt_dpdr(i,11)*epot_p(59) + epot_p(11)*gradt_dpdr(i,59) - gradt_dpdr(i,329) - &
                    & gradt_dpdr(i,338)
                gradt_dpdr(i,343) = gradt_dpdr(i,11)*epot_p(60) + epot_p(11)*gradt_dpdr(i,60) - gradt_dpdr(i,331) - &
                    & gradt_dpdr(i,340)
                gradt_dpdr(i,344) = gradt_dpdr(i,3)*epot_p(130) + epot_p(3)*gradt_dpdr(i,130) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,311) - gradt_dpdr(i,312) - gradt_dpdr(i,309) - gradt_dpdr(i,329) - gradt_dpdr(i, &
                    & 338) - gradt_dpdr(i,331) - gradt_dpdr(i,340) - gradt_dpdr(i,328) - gradt_dpdr(i,321) - &
                    & gradt_dpdr(i,322) - gradt_dpdr(i,311) - gradt_dpdr(i,312)
                gradt_dpdr(i,345) = gradt_dpdr(i,18)*epot_p(42) + epot_p(18)*gradt_dpdr(i,42) - gradt_dpdr(i,310) - &
                    & gradt_dpdr(i,326) - gradt_dpdr(i,330) - gradt_dpdr(i,310)
                gradt_dpdr(i,346) = gradt_dpdr(i,2)*epot_p(132) + epot_p(2)*gradt_dpdr(i,132) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,335)
                gradt_dpdr(i,347) = gradt_dpdr(i,6)*epot_p(92)  + epot_p(6)*gradt_dpdr(i,92) - gradt_dpdr(i,334)
                gradt_dpdr(i,348) = gradt_dpdr(i,7)*epot_p(92)  + epot_p(7)*gradt_dpdr(i,92) - gradt_dpdr(i,333)
                gradt_dpdr(i,349) = gradt_dpdr(i,3)*epot_p(132) + epot_p(3)*gradt_dpdr(i,132) - gradt_dpdr(i,316) - &
                    & gradt_dpdr(i,315) - gradt_dpdr(i,314) - gradt_dpdr(i,338) - gradt_dpdr(i,347) - gradt_dpdr(i, &
                    & 340) - gradt_dpdr(i,348) - gradt_dpdr(i,336) - gradt_dpdr(i,333) - gradt_dpdr(i,334) - &
                    & gradt_dpdr(i,316) - gradt_dpdr(i,315) - gradt_dpdr(i,314) - gradt_dpdr(i,316) - gradt_dpdr(i, &
                    & 314) - gradt_dpdr(i,316) - gradt_dpdr(i,314)
                gradt_dpdr(i,350) = gradt_dpdr(i,3)*epot_p(133) + epot_p(3)*gradt_dpdr(i,133) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,339) - gradt_dpdr(i,341)
                gradt_dpdr(i,351) = gradt_dpdr(i,4)*epot_p(132) + epot_p(4)*gradt_dpdr(i,132) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,344) - gradt_dpdr(i,342) - gradt_dpdr(i,343) - gradt_dpdr(i,332)
                gradt_dpdr(i,352) = gradt_dpdr(i,1)*epot_p(222) + epot_p(1)*gradt_dpdr(i,222) - gradt_dpdr(i,316) - &
                    & gradt_dpdr(i,350)
                gradt_dpdr(i,353) = gradt_dpdr(i,2)*epot_p(134) + epot_p(2)*gradt_dpdr(i,134) - gradt_dpdr(i,306)
                gradt_dpdr(i,354) = gradt_dpdr(i,2)*epot_p(135) + epot_p(2)*gradt_dpdr(i,135) - gradt_dpdr(i,323)
                gradt_dpdr(i,355) = gradt_dpdr(i,3)*epot_p(134) + epot_p(3)*gradt_dpdr(i,134) - gradt_dpdr(i,319)
                gradt_dpdr(i,356) = gradt_dpdr(i,11)*epot_p(53) + epot_p(11)*gradt_dpdr(i,53) - gradt_dpdr(i,310) - &
                    & gradt_dpdr(i,320)
                gradt_dpdr(i,357) = gradt_dpdr(i,2)*epot_p(144) + epot_p(2)*gradt_dpdr(i,144) - gradt_dpdr(i,329) - &
                    & gradt_dpdr(i,338)
                gradt_dpdr(i,358) = gradt_dpdr(i,2)*epot_p(143) + epot_p(2)*gradt_dpdr(i,143) - gradt_dpdr(i,324) - &
                    & gradt_dpdr(i,335)
                gradt_dpdr(i,359) = gradt_dpdr(i,6)*epot_p(81)  + epot_p(6)*gradt_dpdr(i,81) - gradt_dpdr(i,311) - &
                    & gradt_dpdr(i,324)
                gradt_dpdr(i,360) = gradt_dpdr(i,2)*epot_p(146) + epot_p(2)*gradt_dpdr(i,146) - gradt_dpdr(i,331) - &
                    & gradt_dpdr(i,340)
                gradt_dpdr(i,361) = gradt_dpdr(i,2)*epot_p(150) + epot_p(2)*gradt_dpdr(i,150) - gradt_dpdr(i,344)
                gradt_dpdr(i,362) = gradt_dpdr(i,7)*epot_p(82)  + epot_p(7)*gradt_dpdr(i,82) - gradt_dpdr(i,312) - &
                    & gradt_dpdr(i,324)
                gradt_dpdr(i,363) = gradt_dpdr(i,11)*epot_p(65) + epot_p(11)*gradt_dpdr(i,65) - gradt_dpdr(i,308)
                gradt_dpdr(i,364) = gradt_dpdr(i,2)*epot_p(137) + epot_p(2)*gradt_dpdr(i,137) - gradt_dpdr(i,318)
                gradt_dpdr(i,365) = gradt_dpdr(i,17)*epot_p(45) + epot_p(17)*gradt_dpdr(i,45) - gradt_dpdr(i,307)
                gradt_dpdr(i,366) = gradt_dpdr(i,4)*epot_p(134) + epot_p(4)*gradt_dpdr(i,134) - gradt_dpdr(i,317)
                gradt_dpdr(i,367) = gradt_dpdr(i,6)*epot_p(96)  + epot_p(6)*gradt_dpdr(i,96) - gradt_dpdr(i,332)
                gradt_dpdr(i,368) = gradt_dpdr(i,17)*epot_p(54) + epot_p(17)*gradt_dpdr(i,54)
                gradt_dpdr(i,369) = gradt_dpdr(i,7)*epot_p(96)  + epot_p(7)*gradt_dpdr(i,96) - gradt_dpdr(i,332)
                gradt_dpdr(i,370) = gradt_dpdr(i,17)*epot_p(55) + epot_p(17)*gradt_dpdr(i,55)
                gradt_dpdr(i,371) = gradt_dpdr(i,2)*epot_p(159) + epot_p(2)*gradt_dpdr(i,159) - gradt_dpdr(i,345) - &
                    & gradt_dpdr(i,349)
                gradt_dpdr(i,372) = gradt_dpdr(i,2)*epot_p(147) + epot_p(2)*gradt_dpdr(i,147) - gradt_dpdr(i,313)
                gradt_dpdr(i,373) = gradt_dpdr(i,11)*epot_p(58) + epot_p(11)*gradt_dpdr(i,58) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,328)
                gradt_dpdr(i,374) = gradt_dpdr(i,2)*epot_p(148) + epot_p(2)*gradt_dpdr(i,148) - gradt_dpdr(i,329) - &
                    & gradt_dpdr(i,342)
                gradt_dpdr(i,375) = gradt_dpdr(i,6)*epot_p(98)  + epot_p(6)*gradt_dpdr(i,98) - gradt_dpdr(i,327)
                gradt_dpdr(i,376) = gradt_dpdr(i,2)*epot_p(166) + epot_p(2)*gradt_dpdr(i,166) - gradt_dpdr(i,359)
                gradt_dpdr(i,377) = gradt_dpdr(i,14)*epot_p(54) + epot_p(14)*gradt_dpdr(i,54) - gradt_dpdr(i,323)
                gradt_dpdr(i,378) = gradt_dpdr(i,2)*epot_p(149) + epot_p(2)*gradt_dpdr(i,149) - gradt_dpdr(i,331) - &
                    & gradt_dpdr(i,343)
                gradt_dpdr(i,379) = gradt_dpdr(i,7)*epot_p(98)  + epot_p(7)*gradt_dpdr(i,98) - gradt_dpdr(i,325)
                gradt_dpdr(i,380) = gradt_dpdr(i,7)*epot_p(99)  + epot_p(7)*gradt_dpdr(i,99) - gradt_dpdr(i,311) - &
                    & gradt_dpdr(i,359)
                gradt_dpdr(i,381) = gradt_dpdr(i,2)*epot_p(167) + epot_p(2)*gradt_dpdr(i,167) - gradt_dpdr(i,362)
                gradt_dpdr(i,382) = gradt_dpdr(i,14)*epot_p(55) + epot_p(14)*gradt_dpdr(i,55) - gradt_dpdr(i,323)
                gradt_dpdr(i,383) = gradt_dpdr(i,6)*epot_p(100) + epot_p(6)*gradt_dpdr(i,100) - gradt_dpdr(i,312) - &
                    & gradt_dpdr(i,362)
                gradt_dpdr(i,384) = gradt_dpdr(i,11)*epot_p(66) + epot_p(11)*gradt_dpdr(i,66) - gradt_dpdr(i,344)
                gradt_dpdr(i,385) = gradt_dpdr(i,6)*epot_p(101) + epot_p(6)*gradt_dpdr(i,101) - gradt_dpdr(i,325) - &
                    & gradt_dpdr(i,343) - gradt_dpdr(i,379)
                gradt_dpdr(i,386) = gradt_dpdr(i,7)*epot_p(101) + epot_p(7)*gradt_dpdr(i,101) - gradt_dpdr(i,342) - &
                    & gradt_dpdr(i,327) - gradt_dpdr(i,375)
                gradt_dpdr(i,387) = gradt_dpdr(i,6)*epot_p(103) + epot_p(6)*gradt_dpdr(i,103) - gradt_dpdr(i,313) - &
                    & gradt_dpdr(i,382)
                gradt_dpdr(i,388) = gradt_dpdr(i,11)*epot_p(67) + epot_p(11)*gradt_dpdr(i,67) - gradt_dpdr(i,314)
                gradt_dpdr(i,389) = gradt_dpdr(i,2)*epot_p(152) + epot_p(2)*gradt_dpdr(i,152) - gradt_dpdr(i,333)
                gradt_dpdr(i,390) = gradt_dpdr(i,2)*epot_p(153) + epot_p(2)*gradt_dpdr(i,153) - gradt_dpdr(i,334)
                gradt_dpdr(i,391) = gradt_dpdr(i,2)*epot_p(154) + epot_p(2)*gradt_dpdr(i,154) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,373)
                gradt_dpdr(i,392) = gradt_dpdr(i,1)*epot_p(240) + epot_p(1)*gradt_dpdr(i,240) - gradt_dpdr(i,335) - &
                    & gradt_dpdr(i,373) - gradt_dpdr(i,366)
                gradt_dpdr(i,393) = gradt_dpdr(i,2)*epot_p(155) + epot_p(2)*gradt_dpdr(i,155) - gradt_dpdr(i,338) - &
                    & gradt_dpdr(i,342)
                gradt_dpdr(i,394) = gradt_dpdr(i,2)*epot_p(171) + epot_p(2)*gradt_dpdr(i,171) - gradt_dpdr(i,377)
                gradt_dpdr(i,395) = gradt_dpdr(i,2)*epot_p(157) + epot_p(2)*gradt_dpdr(i,157) - gradt_dpdr(i,340) - &
                    & gradt_dpdr(i,343)
                gradt_dpdr(i,396) = gradt_dpdr(i,2)*epot_p(163) + epot_p(2)*gradt_dpdr(i,163) - gradt_dpdr(i,350) - &
                    & gradt_dpdr(i,351)
                gradt_dpdr(i,397) = gradt_dpdr(i,6)*epot_p(95)  + epot_p(6)*gradt_dpdr(i,95) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,350) - gradt_dpdr(i,340)
                gradt_dpdr(i,398) = gradt_dpdr(i,2)*epot_p(172) + epot_p(2)*gradt_dpdr(i,172) - gradt_dpdr(i,382)
                gradt_dpdr(i,399) = gradt_dpdr(i,7)*epot_p(95)  + epot_p(7)*gradt_dpdr(i,95) - gradt_dpdr(i,315) - &
                    & gradt_dpdr(i,350) - gradt_dpdr(i,338)
                gradt_dpdr(i,400) = gradt_dpdr(i,11)*epot_p(68) + epot_p(11)*gradt_dpdr(i,68) - gradt_dpdr(i,345) - &
                    & gradt_dpdr(i,349)
                gradt_dpdr(i,401) = gradt_dpdr(i,2)*epot_p(175) + epot_p(2)*gradt_dpdr(i,175) - gradt_dpdr(i,380) - &
                    & gradt_dpdr(i,385)
                gradt_dpdr(i,402) = gradt_dpdr(i,6)*epot_p(106) + epot_p(6)*gradt_dpdr(i,106) - gradt_dpdr(i,335) - &
                    & gradt_dpdr(i,343) - gradt_dpdr(i,396)
                gradt_dpdr(i,403) = gradt_dpdr(i,2)*epot_p(176) + epot_p(2)*gradt_dpdr(i,176) - gradt_dpdr(i,383) - &
                    & gradt_dpdr(i,386)
                gradt_dpdr(i,404) = gradt_dpdr(i,7)*epot_p(106) + epot_p(7)*gradt_dpdr(i,106) - gradt_dpdr(i,342) - &
                    & gradt_dpdr(i,335) - gradt_dpdr(i,396)
                gradt_dpdr(i,405) = gradt_dpdr(i,4)*epot_p(150) + epot_p(4)*gradt_dpdr(i,150) - gradt_dpdr(i,328) - &
                    & gradt_dpdr(i,359) - gradt_dpdr(i,362)
                gradt_dpdr(i,406) = gradt_dpdr(i,11)*epot_p(69) + epot_p(11)*gradt_dpdr(i,69) - gradt_dpdr(i,316)
                gradt_dpdr(i,407) = gradt_dpdr(i,2)*epot_p(161) + epot_p(2)*gradt_dpdr(i,161) - gradt_dpdr(i,347)
                gradt_dpdr(i,408) = gradt_dpdr(i,2)*epot_p(162) + epot_p(2)*gradt_dpdr(i,162) - gradt_dpdr(i,348)
                gradt_dpdr(i,409) = gradt_dpdr(i,11)*epot_p(70) + epot_p(11)*gradt_dpdr(i,70) - gradt_dpdr(i,350) - &
                    & gradt_dpdr(i,351)
                gradt_dpdr(i,410) = gradt_dpdr(i,2)*epot_p(180) + epot_p(2)*gradt_dpdr(i,180) - gradt_dpdr(i,397) - &
                    & gradt_dpdr(i,402)
                gradt_dpdr(i,411) = gradt_dpdr(i,6)*epot_p(110) + epot_p(6)*gradt_dpdr(i,110) - gradt_dpdr(i,348)
                gradt_dpdr(i,412) = gradt_dpdr(i,2)*epot_p(181) + epot_p(2)*gradt_dpdr(i,181) - gradt_dpdr(i,399) - &
                    & gradt_dpdr(i,404)
                gradt_dpdr(i,413) = gradt_dpdr(i,7)*epot_p(110) + epot_p(7)*gradt_dpdr(i,110) - gradt_dpdr(i,347)
                gradt_dpdr(i,414) = gradt_dpdr(i,6)*epot_p(112) + epot_p(6)*gradt_dpdr(i,112) - gradt_dpdr(i,316) - &
                    & gradt_dpdr(i,398) - gradt_dpdr(i,316)
                gradt_dpdr(i,415) = gradt_dpdr(i,11)*epot_p(71) + epot_p(11)*gradt_dpdr(i,71) - gradt_dpdr(i,352)
                gradt_dpdr(i,416) = gradt_dpdr(i,2)*epot_p(184) + epot_p(2)*gradt_dpdr(i,184) - gradt_dpdr(i,411)
                gradt_dpdr(i,417) = gradt_dpdr(i,2)*epot_p(185) + epot_p(2)*gradt_dpdr(i,185) - gradt_dpdr(i,413)
                gradt_dpdr(i,418) = gradt_dpdr(i,12)*epot_p(71) + epot_p(12)*gradt_dpdr(i,71) - gradt_dpdr(i,351) - &
                    & gradt_dpdr(i,416) - gradt_dpdr(i,417)
                gradt_dpdr(i,419) = gradt_dpdr(i,2)*epot_p(164) + epot_p(2)*gradt_dpdr(i,164) - gradt_dpdr(i,320) - &
                    & gradt_dpdr(i,356)
                gradt_dpdr(i,420) = gradt_dpdr(i,2)*epot_p(165) + epot_p(2)*gradt_dpdr(i,165) - gradt_dpdr(i,328) - &
                    & gradt_dpdr(i,373)
                gradt_dpdr(i,421) = gradt_dpdr(i,2)*epot_p(170) + epot_p(2)*gradt_dpdr(i,170) - gradt_dpdr(i,337) - &
                    & gradt_dpdr(i,392)
                gradt_dpdr(i,422) = gradt_dpdr(i,1)*epot_p(268) + epot_p(1)*gradt_dpdr(i,268) - gradt_dpdr(i,359)
                gradt_dpdr(i,423) = gradt_dpdr(i,1)*epot_p(269) + epot_p(1)*gradt_dpdr(i,269) - gradt_dpdr(i,362)
                gradt_dpdr(i,424) = gradt_dpdr(i,2)*epot_p(174) + epot_p(2)*gradt_dpdr(i,174) - gradt_dpdr(i,345) - &
                    & gradt_dpdr(i,400)
                gradt_dpdr(i,425) = gradt_dpdr(i,6)*epot_p(102) + epot_p(6)*gradt_dpdr(i,102) - gradt_dpdr(i,345) - &
                    & gradt_dpdr(i,326)
                gradt_dpdr(i,426) = gradt_dpdr(i,7)*epot_p(103) + epot_p(7)*gradt_dpdr(i,103) - gradt_dpdr(i,345) - &
                    & gradt_dpdr(i,330)
                gradt_dpdr(i,427) = gradt_dpdr(i,3)*epot_p(186) + epot_p(3)*gradt_dpdr(i,186) - gradt_dpdr(i,364)
                gradt_dpdr(i,428) = gradt_dpdr(i,18)*epot_p(65) + epot_p(18)*gradt_dpdr(i,65) - gradt_dpdr(i,354)
                gradt_dpdr(i,429) = gradt_dpdr(i,17)*epot_p(66) + epot_p(17)*gradt_dpdr(i,66) - gradt_dpdr(i,361)
                gradt_dpdr(i,430) = gradt_dpdr(i,2)*epot_p(179) + epot_p(2)*gradt_dpdr(i,179) - gradt_dpdr(i,350) - &
                    & gradt_dpdr(i,409)
                gradt_dpdr(i,431) = gradt_dpdr(i,4)*epot_p(166) + epot_p(4)*gradt_dpdr(i,166) - gradt_dpdr(i,361) - &
                    & gradt_dpdr(i,357) - gradt_dpdr(i,422)
                gradt_dpdr(i,432) = gradt_dpdr(i,4)*epot_p(167) + epot_p(4)*gradt_dpdr(i,167) - gradt_dpdr(i,361) - &
                    & gradt_dpdr(i,360) - gradt_dpdr(i,423)
                gradt_dpdr(i,433) = gradt_dpdr(i,2)*epot_p(187) + epot_p(2)*gradt_dpdr(i,187) - gradt_dpdr(i,387)
                gradt_dpdr(i,434) = gradt_dpdr(i,14)*epot_p(66) + epot_p(14)*gradt_dpdr(i,66) - gradt_dpdr(i,328) - &
                    & gradt_dpdr(i,405) - gradt_dpdr(i,376) - gradt_dpdr(i,381) - gradt_dpdr(i,373)
                gradt_dpdr(i,435) = gradt_dpdr(i,1)*epot_p(277) + epot_p(1)*gradt_dpdr(i,277) - gradt_dpdr(i,387) - &
                    & gradt_dpdr(i,385) - gradt_dpdr(i,425)
                gradt_dpdr(i,436) = gradt_dpdr(i,1)*epot_p(278) + epot_p(1)*gradt_dpdr(i,278) - gradt_dpdr(i,387) - &
                    & gradt_dpdr(i,386) - gradt_dpdr(i,426)
                gradt_dpdr(i,437) = gradt_dpdr(i,2)*epot_p(177) + epot_p(2)*gradt_dpdr(i,177) - gradt_dpdr(i,346) - &
                    & gradt_dpdr(i,391)
                gradt_dpdr(i,438) = gradt_dpdr(i,2)*epot_p(178) + epot_p(2)*gradt_dpdr(i,178) - gradt_dpdr(i,349) - &
                    & gradt_dpdr(i,400)
                gradt_dpdr(i,439) = gradt_dpdr(i,1)*epot_p(281) + epot_p(1)*gradt_dpdr(i,281) - gradt_dpdr(i,396) - &
                    & gradt_dpdr(i,392) - gradt_dpdr(i,430)
                gradt_dpdr(i,440) = gradt_dpdr(i,21)*epot_p(54) + epot_p(21)*gradt_dpdr(i,54) - gradt_dpdr(i,370)
                gradt_dpdr(i,441) = gradt_dpdr(i,21)*epot_p(55) + epot_p(21)*gradt_dpdr(i,55) - gradt_dpdr(i,368)
                gradt_dpdr(i,442) = gradt_dpdr(i,2)*epot_p(189) + epot_p(2)*gradt_dpdr(i,189) - gradt_dpdr(i,405) - &
                    & gradt_dpdr(i,434)
                gradt_dpdr(i,443) = gradt_dpdr(i,1)*epot_p(285) + epot_p(1)*gradt_dpdr(i,285) - gradt_dpdr(i,402) - &
                    & gradt_dpdr(i,404) - gradt_dpdr(i,400) - gradt_dpdr(i,434)  - gradt_dpdr(i,430)
                gradt_dpdr(i,444) = gradt_dpdr(i,6)*epot_p(116) + epot_p(6)*gradt_dpdr(i,116) - gradt_dpdr(i,347) - &
                    & gradt_dpdr(i,413) - gradt_dpdr(i,403)
                gradt_dpdr(i,445) = gradt_dpdr(i,7)*epot_p(116) + epot_p(7)*gradt_dpdr(i,116) - gradt_dpdr(i,348) - &
                    & gradt_dpdr(i,411) - gradt_dpdr(i,401)
                gradt_dpdr(i,446) = gradt_dpdr(i,2)*epot_p(182) + epot_p(2)*gradt_dpdr(i,182) - gradt_dpdr(i,351) - &
                    & gradt_dpdr(i,409)
                gradt_dpdr(i,447) = gradt_dpdr(i,2)*epot_p(191) + epot_p(2)*gradt_dpdr(i,191) - gradt_dpdr(i,414) - &
                    & gradt_dpdr(i,443)
                gradt_dpdr(i,448) = gradt_dpdr(i,3)*epot_p(183) + epot_p(3)*gradt_dpdr(i,183) - gradt_dpdr(i,351) - &
                    & gradt_dpdr(i,418) - gradt_dpdr(i,411) - gradt_dpdr(i,413) - gradt_dpdr(i,409)
                gradt_dpdr(i,449) = gradt_dpdr(i,6)*epot_p(118) + epot_p(6)*gradt_dpdr(i,118) - gradt_dpdr(i,351) - &
                    & gradt_dpdr(i,418) - gradt_dpdr(i,412)
                gradt_dpdr(i,450) = gradt_dpdr(i,7)*epot_p(118) + epot_p(7)*gradt_dpdr(i,118) - gradt_dpdr(i,351) - &
                    & gradt_dpdr(i,418) - gradt_dpdr(i,410)
                gradt_dpdr(i,451) = gradt_dpdr(i,2)*epot_p(193) + epot_p(2)*gradt_dpdr(i,193) - gradt_dpdr(i,418) - &
                    & gradt_dpdr(i,448)
                gradt_dpdr(i,452) = gradt_dpdr(i,1)*epot_p(294) + epot_p(1)*gradt_dpdr(i,294) - gradt_dpdr(i,418) - &
                    & gradt_dpdr(i,415) - gradt_dpdr(i,448)
                gradt_dpdr(i,453) = gradt_dpdr(i,6)*epot_p(119) + epot_p(6)*gradt_dpdr(i,119) - gradt_dpdr(i,417)
                gradt_dpdr(i,454) = gradt_dpdr(i,7)*epot_p(119) + epot_p(7)*gradt_dpdr(i,119) - gradt_dpdr(i,416)
                gradt_dpdr(i,455) = gradt_dpdr(i,2)*epot_p(186) + epot_p(2)*gradt_dpdr(i,186) - gradt_dpdr(i,363)
                gradt_dpdr(i,456) = gradt_dpdr(i,3)*epot_p(187) + epot_p(3)*gradt_dpdr(i,187) - gradt_dpdr(i,385) - &
                    & gradt_dpdr(i,386) - gradt_dpdr(i,384) - gradt_dpdr(i,435) - gradt_dpdr(i,436) - gradt_dpdr(i,434)
                gradt_dpdr(i,457) = gradt_dpdr(i,2)*epot_p(188) + epot_p(2)*gradt_dpdr(i,188) - gradt_dpdr(i,388)
                gradt_dpdr(i,458) = gradt_dpdr(i,4)*epot_p(187) + epot_p(4)*gradt_dpdr(i,187) - gradt_dpdr(i,425) - &
                    & gradt_dpdr(i,426) - gradt_dpdr(i,424)
                gradt_dpdr(i,459) = gradt_dpdr(i,2)*epot_p(190) + epot_p(2)*gradt_dpdr(i,190) - gradt_dpdr(i,406)
                gradt_dpdr(i,460) = gradt_dpdr(i,21)*epot_p(66) + epot_p(21)*gradt_dpdr(i,66) - gradt_dpdr(i,422) - &
                    & gradt_dpdr(i,423) - gradt_dpdr(i,420)
                gradt_dpdr(i,461) = gradt_dpdr(i,2)*epot_p(192) + epot_p(2)*gradt_dpdr(i,192) - gradt_dpdr(i,415)
                gradt_dpdr(i,462) = gradt_dpdr(i,18)*epot_p(71) + epot_p(18)*gradt_dpdr(i,71) - gradt_dpdr(i,440) - &
                    & gradt_dpdr(i,441) - gradt_dpdr(i,438)
                gradt_dpdr(i,463) = gradt_dpdr(i,2)*epot_p(194) + epot_p(2)*gradt_dpdr(i,194) - gradt_dpdr(i,452)
                gradt_dpdr(i,464) = gradt_dpdr(i,3)*epot_p(194) + epot_p(3)*gradt_dpdr(i,194) - gradt_dpdr(i,453) - &
                    & gradt_dpdr(i,454) - gradt_dpdr(i,451)
                gradt_dpdr(i,465) = gradt_dpdr(i,1)*epot_p(305) + epot_p(1)*gradt_dpdr(i,305) - gradt_dpdr(i,464) - &
                    & gradt_dpdr(i,463)
            enddo
        end subroutine evdpdr
        !> @brief Evaluate evdbdr for the O4_triplet_v2 potential.
        subroutine evdbdr
            !  This subroutine eliminates the 2-body terms in Bowman's approach.
            implicit none
            integer          :: j
            integer          :: i
            double precision :: db1dr(6,466)
            ! Pass P(0:465) to BM1(1:466)
            do j=1, 6
                do i=1, 466
                    db1dr(j,i) = gradt_dpdr(j,i-1)
                enddo
            enddo
            ! Remove unconnected terms and 2-body terms and pass to B(1:430)
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
                do i=277, 425
                    gradt_dbdr(j,i) = db1dr(j,i+30)
                enddo
                gradt_dbdr(j,426) = db1dr(j,457)
                gradt_dbdr(j,427) = db1dr(j,459)
                gradt_dbdr(j,428) = db1dr(j,461)
                gradt_dbdr(j,429) = db1dr(j,463)
                gradt_dbdr(j,430) = db1dr(j,465)
            enddo
        end subroutine evdbdr
        !> @brief Evaluate d3disp for the O4_triplet_v2 potential.
        subroutine d3disp( dist, disp, dispdr, igrad )
            implicit none
            double precision :: cn(4), s6, s8, rs6, rs8
            double precision :: dist(6), e6(6), e8(6), disp, dispdr(6), c6(6)
            double precision :: e6dr(6), e8dr(6)
            integer          :: iz(4), mxc(94), i, j, igrad
            double precision :: c6ab(94,94,5,5,3)
            double precision :: r2r4(94)
            double precision :: autoang, autokcal
            autoang  = 0.52917726d0
            autokcal = 627.509541d0
            ! Generalized parameters for BJ damping
            s6  = 1.0d0
            s8  = 2.0d0
            rs6 = 0.5299d0
            rs8 = 2.20d0
            do i=1, 6
                dist(i) = dist(i)/autoang
            enddo
            ! iz for O4 system
            iz(1) = 8
            iz(2) = 8
            iz(3) = 8
            iz(4) = 8
            ! C6 for O4 system
            c6(1) = 12.8d0
            c6(2) = 12.8d0
            c6(3) = 12.8d0
            c6(4) = 12.8d0
            c6(5) = 12.8d0
            c6(6) = 12.8d0
            ! Calculate dispersion correction
            call edisp(94, 5, 4, dist, iz, mxc, &
                & rs6, rs8, e6, e8, e6dr, e8dr, c6, 0)
            disp = 0.0d0
            do i=1, 6
                disp = disp + (-s6*e6(i)-s8*e8(i))*autokcal
            enddo
            if( igrad .eq. 1 )then
                call edisp(94, 5, 4, dist, iz, mxc, &
                    & rs6, rs8, e6, e8, e6dr, e8dr, c6, 1)
                dispdr(:) = 0.0d0
                do i=1, 6
                    dispdr(i) = dispdr(i) + (-s6*e6dr(i)-s8*e8dr(i))*autokcal/autoang
                enddo
            endif
            do i=1, 6
                dist(i) = dist(i)*autoang
            enddo
        end subroutine d3disp
        !> @brief Evaluate edisp for the O4_triplet_v2 potential.
        subroutine edisp( max_elem, maxc, n, dist, iz, mxc, &
                & rs6, rs8, e6, e8, e6dr, e8dr, c6a, igrad )
            implicit none
            integer          :: n, iz(4), max_elem, maxc, mxc(max_elem)
            double precision :: dist(6), r2r4(max_elem), r0ab(max_elem,max_elem)
            double precision :: rs6, rs8, rcov(max_elem)
            double precision :: c6ab(max_elem,max_elem,maxc,maxc,3)
            double precision :: e6(6), e8(6), c6a(6), e6dr(6), e8dr(6)
            integer          :: iat, jat, igrad
            double precision :: r, tmp, c6, c8, a1, a2
            double precision :: damp6, damp8
            double precision :: cn(n)
            double precision :: r2ab(n*n), cc6ab(n*n), dmp(n*n)
            integer          :: step
            e6(:)   = 0.0d0
            e8(:)   = 0.0d0
            e6dr(:) = 0.0d0
            e8dr(:) = 0.0d0
            a1      = rs6
            a2      = rs8
            !  r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
            r2r4(1:94)=(/ 2.00734898, 1.56637132, 5.01986934, 3.85379032, 3.64446594, 3.10492822, 2.71175247, &
                & 2.59361680, 2.38825250, 2.21522516, 6.58585536, 5.46295967, 5.65216669, 4.88284902, 4.29727576, &
                & 4.04108902, 3.72932356, 3.44677275, 7.97762753, 7.07623947, 6.60844053, 6.28791364, 6.07728703, &
                & 5.54643096, 5.80491167, 5.58415602, 5.41374528, 5.28497229, 5.22592821, 5.09817141, 6.12149689, &
                & 5.54083734, 5.06696878, 4.87005108, 4.59089647, 4.31176304, 9.55461698, 8.67396077, 7.97210197, &
                & 7.43439917, 6.58711862, 6.19536215, 6.01517290, 5.81623410, 5.65710424, 5.52640661, 5.44263305, &
                & 5.58285373, 7.02081898, 6.46815523, 5.98089120, 5.81686657, 5.53321815, 5.25477007, 11.02204549, &
                & 10.15679528, 9.35167836, 9.06926079, 8.97241155, 8.90092807, 8.85984840, 8.81736827, 8.79317710, &
                & 7.89969626, 8.80588454, 8.42439218, 8.54289262, 8.47583370, 8.45090888, 8.47339339, 7.83525634, &
                & 8.20702843, 7.70559063, 7.32755997, 7.03887381, 6.68978720, 6.05450052, 5.88752022, 5.70661499, &
                & 5.78450695, 7.79780729, 7.26443867, 6.78151984, 6.67883169, 6.39024318, 6.09527958, 11.79156076, &
                & 11.10997644, 9.51377795, 8.67197068, 8.77140725, 8.65402716, 8.53923501, 8.85024712 /)
            ! these new data are scaled with k2=4./3. and converted to a_0 via
            ! autoang=0.52917726d0
            rcov(1:94)=(/ 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, 1.88972601, 1.78894056, &
                & 1.58736983, 1.61256616, 1.68815527, 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
                & 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, 3.35111422, 3.07395437, 3.04875805, &
                & 2.77159820, 2.69600923, 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, 2.82199085, &
                & 2.74640188, 2.89757982, 2.77159820, 2.87238349, 2.94797246, 4.76210950, 4.20778980, 3.70386304, &
                & 3.50229216, 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, 2.72120556, 2.89757982, &
                & 3.09915070, 3.22513231, 3.17473967, 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
                & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, 3.93062995, 3.90543362, 3.80464833, &
                & 3.82984466, 3.80464833, 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, 3.67866672, &
                & 3.45189952, 3.30072128, 3.09915070, 2.97316878, 2.92277614, 2.79679452, 2.82199085, 2.84718717, &
                & 3.32591790, 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, 3.57788113, 5.06446567, &
                & 4.56053862, 4.20778980, 3.98102289, 3.82984466, 3.85504098, 3.88023730, 3.90543362 /)
            ! DFT-D3
            step = 0
            do iat=1, n-1
                do jat=iat+1, n
                    step = step+1
                    r    = dist(step)
                    c6   = c6a(step)
                    ! r2r4 stored in main as sqrt
                    c8 = 3.0d0*c6*r2r4(iz(iat))*r2r4(iz(jat))
                    ! energy for BJ damping
                    tmp      = sqrt(c8/c6)
                    e6(step) = c6/(r**6+(a1*tmp+a2)**6)
                    e8(step) = c8/(r**8+(a1*tmp+a2)**8)
                    ! calculate gradients
                    if( igrad .eq. 1 )then
                        ! grad for BJ damping
                        e6dr(step) = c6*(-6*r**5)/(r**6+(a1*tmp+a2)**6)**2
                        e8dr(step) = c8*(-8*r**7)/(r**8+(a1*tmp+a2)**8)**2
                    endif
                enddo
            enddo
        end subroutine edisp
    end subroutine evaluate
end module chempotpy_o4_triplet_v2_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_o4_triplet_v2_pes, only: evaluate
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
