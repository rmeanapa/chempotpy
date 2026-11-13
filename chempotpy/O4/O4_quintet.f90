!   System:                        O4
!   Functional form:               permutation-invariant polynomials
!   Common name:                   O4 quintet (adiabatic ground state)
!   Number of derivatives:         1
!   Number of bodies:              4
!   Number of electronic surfaces: 1
!   Interface:                     Section-2
!   References:: Y. Paukku, K. R. Yang, Z. Varga, G. Song, J. D. Bender,
!                D. G. Truhlar, J. Chem. Phys. 147, 034301/1-11 (2017)
!   Notes:    PES of quintet O4 with special emphasize for
!             O2 + O2 --> O2 + O + O
!            - Fit is based on extended dataset
!              containing 10543 points
!            - Mixed-exponential- gaussian (MEG)
!              variables are applied to describe the
!              long-range interactions
!            - D3 dispersion correction is implemented for diatoms
!     O1--O2
!     O3--O4

!> @brief Implement O4_quintet, a potential energy surface for O4.
module chempotpy_o4_quintet_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the O4_quintet potential energy surface for O4.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units, only: electronvolt_to_kilocalorie_per_mole
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 4
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: coord_r(6)
        real(wp)              :: epot_b(430)
        real(wp)              :: epot_c(430)
        real(wp)              :: epot_p(0:465)
        real(wp)              :: epot_rm(0:111)
        real(wp)              :: epot_rms(6)
        real(wp)              :: gradt_dbdr(6,430)
        real(wp)              :: gradt_dmdr(6,0:111)
        real(wp)              :: gradt_dmsdr(6,6)
        real(wp)              :: gradt_dpdr(6,0:465)
        real(wp)              :: gradt_drdx(6,12)
        real(wp)              :: gradt_dvdr(6)
        real(wp)              :: msprmt_a
        real(wp)              :: msprmt_ab
        real(wp)              :: msprmt_re
        real(wp)              :: v, tx(12), tg(12)
        integer               :: iatom, idir, j, istate
        !initialize
        coord_r      = 0.e0_wp; epot_b = 0.e0_wp; epot_c = 0.e0_wp; epot_p = 0.e0_wp; epot_rm = 0.e0_wp; epot_rms = 0.e0_wp; gradt_dbdr = 0.e0_wp
        gradt_dmdr   = 0.e0_wp; gradt_dmsdr = 0.e0_wp; gradt_dpdr = 0.e0_wp; gradt_drdx = 0.e0_wp; gradt_dvdr = 0.e0_wp; msprmt_a = 0.e0_wp
        msprmt_ab    = 0.e0_wp; msprmt_re = 0.e0_wp; msprmt_a = 1.1e0_wp; msprmt_ab = 1.5e0_wp; msprmt_re = 1.208e0_wp
        epot_c(   1) = -0.510307106137e+02_wp
        epot_c(   2) = 0.152013704450e+04_wp
        epot_c(   3) = 0.702040380028e+03_wp
        epot_c(   4) = -0.189549071273e+04_wp
        epot_c(   5) = 0.360584974224e+03_wp
        epot_c(   6) = 0.308127051642e+04_wp
        epot_c(   7) = 0.522152294002e+04_wp
        epot_c(   8) = -0.551032010095e+04_wp
        epot_c(   9) = -0.760116812667e+04_wp
        epot_c(  10) = -0.489616491044e+04_wp
        epot_c(  11) = 0.114503095229e+04_wp
        epot_c(  12) = 0.694353589977e+04_wp
        epot_c(  13) = -0.944972881209e+03_wp
        epot_c(  14) = -0.448350361660e+04_wp
        epot_c(  15) = 0.744328826942e+03_wp
        epot_c(  16) = -0.166069418117e+05_wp
        epot_c(  17) = -0.350500168935e+05_wp
        epot_c(  18) = 0.119751621658e+05_wp
        epot_c(  19) = 0.151254753150e+05_wp
        epot_c(  20) = 0.127366352123e+05_wp
        epot_c(  21) = 0.115026737908e+05_wp
        epot_c(  22) = -0.120144773118e+05_wp
        epot_c(  23) = 0.108996977211e+05_wp
        epot_c(  24) = 0.190252591139e+05_wp
        epot_c(  25) = 0.464172231847e+04_wp
        epot_c(  26) = 0.100731093634e+05_wp
        epot_c(  27) = -0.152784530393e+05_wp
        epot_c(  28) = -0.932727432671e+03_wp
        epot_c(  29) = 0.192822098873e+06_wp
        epot_c(  30) = -0.261777445045e+05_wp
        epot_c(  31) = 0.673970374508e+05_wp
        epot_c(  32) = -0.135291303011e+05_wp
        epot_c(  33) = 0.768377903110e+05_wp
        epot_c(  34) = 0.348971795640e+05_wp
        epot_c(  35) = -0.299247786372e+05_wp
        epot_c(  36) = 0.114396343819e+05_wp
        epot_c(  37) = -0.349772875041e+05_wp
        epot_c(  38) = -0.374183986788e+05_wp
        epot_c(  39) = 0.409507131624e+05_wp
        epot_c(  40) = 0.646522945728e+05_wp
        epot_c(  41) = -0.233491096554e+05_wp
        epot_c(  42) = -0.481404447423e+05_wp
        epot_c(  43) = 0.191243065502e+05_wp
        epot_c(  44) = 0.143418012720e+05_wp
        epot_c(  45) = -0.193846545774e+05_wp
        epot_c(  46) = -0.208233298400e+05_wp
        epot_c(  47) = 0.413795523782e+04_wp
        epot_c(  48) = -0.178794922080e+05_wp
        epot_c(  49) = 0.217365997113e+05_wp
        epot_c(  50) = -0.779025994199e+04_wp
        epot_c(  51) = -0.478679487639e+05_wp
        epot_c(  52) = 0.141953259269e+05_wp
        epot_c(  53) = -0.325087673057e+05_wp
        epot_c(  54) = 0.384547469912e+05_wp
        epot_c(  55) = 0.162777907591e+05_wp
        epot_c(  56) = 0.309637649927e+04_wp
        epot_c(  57) = -0.131766366127e+06_wp
        epot_c(  58) = 0.697147369273e+05_wp
        epot_c(  59) = 0.424260868457e+05_wp
        epot_c(  60) = 0.714505919981e+05_wp
        epot_c(  61) = -0.462441555428e+05_wp
        epot_c(  62) = 0.131868582560e+05_wp
        epot_c(  63) = -0.593828467900e+05_wp
        epot_c(  64) = -0.512285289738e+05_wp
        epot_c(  65) = -0.267396244032e+05_wp
        epot_c(  66) = -0.293523105895e+05_wp
        epot_c(  67) = -0.327071965736e+04_wp
        epot_c(  68) = 0.419582761819e+05_wp
        epot_c(  69) = -0.727168820408e+05_wp
        epot_c(  70) = -0.111558301125e+06_wp
        epot_c(  71) = -0.204223217692e+05_wp
        epot_c(  72) = 0.281162250737e+05_wp
        epot_c(  73) = 0.717115867098e+05_wp
        epot_c(  74) = -0.373288985801e+05_wp
        epot_c(  75) = 0.297368128373e+05_wp
        epot_c(  76) = 0.786320441836e+05_wp
        epot_c(  77) = 0.634406834466e+05_wp
        epot_c(  78) = -0.398191838378e+05_wp
        epot_c(  79) = -0.410848233354e+05_wp
        epot_c(  80) = -0.236967649327e+05_wp
        epot_c(  81) = 0.869012967583e+04_wp
        epot_c(  82) = 0.608216006948e+05_wp
        epot_c(  83) = 0.544073252898e+04_wp
        epot_c(  84) = -0.163481530703e+05_wp
        epot_c(  85) = -0.174872643458e+05_wp
        epot_c(  86) = -0.174933092212e+05_wp
        epot_c(  87) = -0.416656147984e+05_wp
        epot_c(  88) = -0.705209143942e+05_wp
        epot_c(  89) = 0.145961809449e+05_wp
        epot_c(  90) = 0.169336076298e+05_wp
        epot_c(  91) = 0.560117299053e+04_wp
        epot_c(  92) = 0.470027890016e+05_wp
        epot_c(  93) = -0.107340638355e+05_wp
        epot_c(  94) = -0.590457881328e+04_wp
        epot_c(  95) = 0.695182719331e+05_wp
        epot_c(  96) = -0.466772368967e+05_wp
        epot_c(  97) = 0.557575819406e+05_wp
        epot_c(  98) = -0.162153283191e+05_wp
        epot_c(  99) = -0.409572300690e+05_wp
        epot_c( 100) = 0.657176261329e+04_wp
        epot_c( 101) = 0.365869515549e+05_wp
        epot_c( 102) = 0.714740505571e+05_wp
        epot_c( 103) = 0.304351919989e+05_wp
        epot_c( 104) = 0.134107887391e+05_wp
        epot_c( 105) = -0.905525657746e+05_wp
        epot_c( 106) = -0.469306631319e+05_wp
        epot_c( 107) = 0.485210273391e+05_wp
        epot_c( 108) = -0.493741519903e+05_wp
        epot_c( 109) = -0.326582246691e+05_wp
        epot_c( 110) = -0.791203184349e+05_wp
        epot_c( 111) = 0.615395268964e+05_wp
        epot_c( 112) = 0.526454043626e+05_wp
        epot_c( 113) = 0.196253010880e+05_wp
        epot_c( 114) = 0.157012027807e+05_wp
        epot_c( 115) = 0.388731994724e+04_wp
        epot_c( 116) = 0.543270505041e+05_wp
        epot_c( 117) = -0.478098147077e+05_wp
        epot_c( 118) = 0.742896921260e+05_wp
        epot_c( 119) = 0.228879220180e+05_wp
        epot_c( 120) = 0.194548206486e+05_wp
        epot_c( 121) = 0.137663654388e+05_wp
        epot_c( 122) = 0.232524805536e+05_wp
        epot_c( 123) = 0.323648720523e+05_wp
        epot_c( 124) = -0.127768713938e+04_wp
        epot_c( 125) = -0.568425575517e+02_wp
        epot_c( 126) = 0.222977436754e+05_wp
        epot_c( 127) = -0.419461937707e+05_wp
        epot_c( 128) = 0.448197349719e+05_wp
        epot_c( 129) = -0.705712320284e+05_wp
        epot_c( 130) = -0.139508690794e+05_wp
        epot_c( 131) = -0.754023766262e+05_wp
        epot_c( 132) = -0.210133695356e+05_wp
        epot_c( 133) = 0.260690789849e+05_wp
        epot_c( 134) = 0.628402824141e+05_wp
        epot_c( 135) = -0.166785947231e+05_wp
        epot_c( 136) = -0.196387600130e+04_wp
        epot_c( 137) = -0.320126646380e+05_wp
        epot_c( 138) = 0.364726282275e+05_wp
        epot_c( 139) = -0.201126324358e+05_wp
        epot_c( 140) = -0.637492206610e+05_wp
        epot_c( 141) = -0.470165005984e+05_wp
        epot_c( 142) = 0.427448027960e+04_wp
        epot_c( 143) = -0.421886803752e+04_wp
        epot_c( 144) = 0.596967070038e+05_wp
        epot_c( 145) = -0.242932629820e+05_wp
        epot_c( 146) = -0.265855057514e+05_wp
        epot_c( 147) = -0.122495388487e+05_wp
        epot_c( 148) = 0.437525439748e+05_wp
        epot_c( 149) = -0.422669718398e+04_wp
        epot_c( 150) = -0.416601179186e+05_wp
        epot_c( 151) = 0.335811865831e+05_wp
        epot_c( 152) = 0.567559120175e+05_wp
        epot_c( 153) = -0.143996770828e+05_wp
        epot_c( 154) = 0.394197970540e+04_wp
        epot_c( 155) = 0.272627603959e+05_wp
        epot_c( 156) = 0.341689603878e+05_wp
        epot_c( 157) = 0.536544688440e+05_wp
        epot_c( 158) = -0.918899354247e+03_wp
        epot_c( 159) = 0.378103419779e+04_wp
        epot_c( 160) = -0.330065395607e+05_wp
        epot_c( 161) = -0.754760171191e+05_wp
        epot_c( 162) = 0.144422257472e+04_wp
        epot_c( 163) = 0.114743451986e+05_wp
        epot_c( 164) = -0.519428980791e+05_wp
        epot_c( 165) = 0.645816290971e+05_wp
        epot_c( 166) = -0.559562864270e+05_wp
        epot_c( 167) = 0.351434441601e+05_wp
        epot_c( 168) = -0.222876206069e+05_wp
        epot_c( 169) = 0.645932799169e+05_wp
        epot_c( 170) = -0.188756450257e+05_wp
        epot_c( 171) = -0.411557452775e+05_wp
        epot_c( 172) = 0.871164609151e+05_wp
        epot_c( 173) = 0.109427636168e+06_wp
        epot_c( 174) = -0.214205078353e+05_wp
        epot_c( 175) = -0.536046993412e+05_wp
        epot_c( 176) = 0.689673470448e+05_wp
        epot_c( 177) = -0.875029866881e+05_wp
        epot_c( 178) = 0.248733257308e+05_wp
        epot_c( 179) = 0.284065112563e+05_wp
        epot_c( 180) = 0.527907246958e+05_wp
        epot_c( 181) = 0.400805997696e+05_wp
        epot_c( 182) = -0.748862499635e+05_wp
        epot_c( 183) = -0.257023688521e+05_wp
        epot_c( 184) = 0.933105553485e+05_wp
        epot_c( 185) = 0.192805527414e+05_wp
        epot_c( 186) = 0.133533877002e+05_wp
        epot_c( 187) = 0.273476233435e+05_wp
        epot_c( 188) = -0.110970161308e+06_wp
        epot_c( 189) = 0.452363098065e+05_wp
        epot_c( 190) = 0.317357392909e+05_wp
        epot_c( 191) = -0.500130837286e+05_wp
        epot_c( 192) = 0.974073881056e+04_wp
        epot_c( 193) = -0.135260998634e+05_wp
        epot_c( 194) = 0.316465658569e+05_wp
        epot_c( 195) = -0.167969409973e+05_wp
        epot_c( 196) = -0.122333689905e+05_wp
        epot_c( 197) = -0.236503786581e+05_wp
        epot_c( 198) = 0.491950077131e+04_wp
        epot_c( 199) = -0.375060734982e+05_wp
        epot_c( 200) = 0.119839064062e+05_wp
        epot_c( 201) = -0.632683840358e+05_wp
        epot_c( 202) = 0.782970258033e+04_wp
        epot_c( 203) = -0.616381849464e+05_wp
        epot_c( 204) = 0.877762760481e+04_wp
        epot_c( 205) = -0.890632421155e+04_wp
        epot_c( 206) = 0.399574195744e+05_wp
        epot_c( 207) = 0.133345057057e+05_wp
        epot_c( 208) = 0.429744773402e+03_wp
        epot_c( 209) = -0.510067816928e+05_wp
        epot_c( 210) = 0.210420794625e+05_wp
        epot_c( 211) = 0.116771418497e+05_wp
        epot_c( 212) = 0.276416046642e+05_wp
        epot_c( 213) = -0.191126669127e+05_wp
        epot_c( 214) = -0.497838484987e+05_wp
        epot_c( 215) = 0.227434935870e+04_wp
        epot_c( 216) = 0.204783224371e+05_wp
        epot_c( 217) = -0.228229605635e+05_wp
        epot_c( 218) = -0.169282578950e+05_wp
        epot_c( 219) = -0.118567829402e+05_wp
        epot_c( 220) = -0.184635101201e+05_wp
        epot_c( 221) = 0.217971884985e+05_wp
        epot_c( 222) = -0.335686945048e+05_wp
        epot_c( 223) = 0.251965657435e+05_wp
        epot_c( 224) = -0.209516813660e+05_wp
        epot_c( 225) = 0.283837054949e+05_wp
        epot_c( 226) = 0.138575344046e+05_wp
        epot_c( 227) = 0.356662443318e+04_wp
        epot_c( 228) = 0.145872016971e+05_wp
        epot_c( 229) = 0.404623669906e+05_wp
        epot_c( 230) = -0.154213084865e+05_wp
        epot_c( 231) = 0.150964575232e+05_wp
        epot_c( 232) = 0.547052418883e+04_wp
        epot_c( 233) = 0.248812674183e+05_wp
        epot_c( 234) = 0.575781808631e+03_wp
        epot_c( 235) = -0.252983623231e+04_wp
        epot_c( 236) = -0.145760802893e+05_wp
        epot_c( 237) = -0.264655411633e+04_wp
        epot_c( 238) = 0.268789876250e+05_wp
        epot_c( 239) = 0.196337921266e+05_wp
        epot_c( 240) = 0.599487068073e+04_wp
        epot_c( 241) = 0.907170535070e+04_wp
        epot_c( 242) = -0.385787216363e+05_wp
        epot_c( 243) = 0.127360926778e+05_wp
        epot_c( 244) = -0.333260849874e+05_wp
        epot_c( 245) = -0.276667324002e+05_wp
        epot_c( 246) = -0.690896140070e+04_wp
        epot_c( 247) = 0.531155525030e+04_wp
        epot_c( 248) = 0.218346541254e+05_wp
        epot_c( 249) = -0.163468389759e+05_wp
        epot_c( 250) = 0.198206951580e+05_wp
        epot_c( 251) = -0.160613533852e+05_wp
        epot_c( 252) = -0.262658343058e+05_wp
        epot_c( 253) = -0.337292252138e+05_wp
        epot_c( 254) = -0.273555473972e+05_wp
        epot_c( 255) = 0.780883098252e+04_wp
        epot_c( 256) = 0.103119987171e+05_wp
        epot_c( 257) = -0.588945308438e+04_wp
        epot_c( 258) = -0.820258852264e+05_wp
        epot_c( 259) = 0.231203712109e+05_wp
        epot_c( 260) = 0.714408377876e+04_wp
        epot_c( 261) = 0.124519779737e+05_wp
        epot_c( 262) = -0.320954486494e+03_wp
        epot_c( 263) = -0.239663176424e+05_wp
        epot_c( 264) = -0.755342506767e+04_wp
        epot_c( 265) = -0.116758179678e+05_wp
        epot_c( 266) = 0.168358732882e+05_wp
        epot_c( 267) = 0.596355046660e+05_wp
        epot_c( 268) = 0.197655665318e+04_wp
        epot_c( 269) = -0.449289741911e+04_wp
        epot_c( 270) = 0.223113843909e+05_wp
        epot_c( 271) = -0.453737763679e+05_wp
        epot_c( 272) = 0.280595124075e+05_wp
        epot_c( 273) = -0.124784725736e+04_wp
        epot_c( 274) = 0.446159786020e+04_wp
        epot_c( 275) = -0.374658266723e+05_wp
        epot_c( 276) = 0.141748554587e+05_wp
        epot_c( 277) = -0.273907047888e+06_wp
        epot_c( 278) = 0.122807601440e+06_wp
        epot_c( 279) = 0.588161701690e+05_wp
        epot_c( 280) = 0.594703606379e+05_wp
        epot_c( 281) = -0.200711338048e+06_wp
        epot_c( 282) = -0.155086098516e+05_wp
        epot_c( 283) = -0.325318142519e+04_wp
        epot_c( 284) = -0.264114724806e+05_wp
        epot_c( 285) = -0.580448831341e+05_wp
        epot_c( 286) = 0.660415466876e+05_wp
        epot_c( 287) = -0.993752404420e+05_wp
        epot_c( 288) = -0.184756777849e+06_wp
        epot_c( 289) = -0.626635853007e+05_wp
        epot_c( 290) = 0.698227922811e+05_wp
        epot_c( 291) = 0.260001259487e+05_wp
        epot_c( 292) = -0.585081386306e+05_wp
        epot_c( 293) = 0.114414673982e+04_wp
        epot_c( 294) = 0.370039135760e+04_wp
        epot_c( 295) = 0.755172183095e+05_wp
        epot_c( 296) = 0.230703799826e+05_wp
        epot_c( 297) = 0.478360609088e+05_wp
        epot_c( 298) = 0.262454477132e+05_wp
        epot_c( 299) = 0.132189293223e+05_wp
        epot_c( 300) = -0.158263431578e+05_wp
        epot_c( 301) = 0.192050804662e+05_wp
        epot_c( 302) = -0.156209756331e+05_wp
        epot_c( 303) = -0.826771163057e+02_wp
        epot_c( 304) = -0.740582883450e+04_wp
        epot_c( 305) = -0.145528651276e+05_wp
        epot_c( 306) = 0.158481752716e+05_wp
        epot_c( 307) = 0.166142844683e+05_wp
        epot_c( 308) = -0.816222634838e+05_wp
        epot_c( 309) = -0.304694875842e+05_wp
        epot_c( 310) = -0.878539688044e+04_wp
        epot_c( 311) = -0.413657632383e+05_wp
        epot_c( 312) = 0.249072629792e+05_wp
        epot_c( 313) = -0.305664150039e+05_wp
        epot_c( 314) = -0.230920204604e+05_wp
        epot_c( 315) = 0.615068223835e+05_wp
        epot_c( 316) = -0.120889953297e+05_wp
        epot_c( 317) = 0.164637448911e+04_wp
        epot_c( 318) = 0.414748789671e+05_wp
        epot_c( 319) = 0.253719412740e+05_wp
        epot_c( 320) = -0.300346454853e+05_wp
        epot_c( 321) = 0.159591233040e+05_wp
        epot_c( 322) = 0.652824812464e+04_wp
        epot_c( 323) = -0.136383102684e+05_wp
        epot_c( 324) = 0.242513795330e+05_wp
        epot_c( 325) = 0.572302226251e+05_wp
        epot_c( 326) = -0.337878542901e+05_wp
        epot_c( 327) = 0.373430894053e+04_wp
        epot_c( 328) = 0.183651694937e+05_wp
        epot_c( 329) = -0.150571498364e+05_wp
        epot_c( 330) = -0.171299354624e+05_wp
        epot_c( 331) = -0.268642857746e+05_wp
        epot_c( 332) = -0.178660784070e+05_wp
        epot_c( 333) = 0.133616164573e+05_wp
        epot_c( 334) = -0.221998771924e+05_wp
        epot_c( 335) = 0.967011806075e+04_wp
        epot_c( 336) = 0.120624474979e+05_wp
        epot_c( 337) = -0.782822493834e+04_wp
        epot_c( 338) = -0.479813127413e+04_wp
        epot_c( 339) = -0.946901219142e+04_wp
        epot_c( 340) = -0.646428111867e+04_wp
        epot_c( 341) = 0.242503818114e+05_wp
        epot_c( 342) = 0.926875675589e+03_wp
        epot_c( 343) = 0.672611252496e+03_wp
        epot_c( 344) = 0.764618624199e+04_wp
        epot_c( 345) = -0.700849579197e+04_wp
        epot_c( 346) = -0.117338336192e+05_wp
        epot_c( 347) = 0.143351002006e+05_wp
        epot_c( 348) = -0.179222755118e+05_wp
        epot_c( 349) = -0.866977720025e+04_wp
        epot_c( 350) = 0.182152735461e+05_wp
        epot_c( 351) = 0.248547542680e+05_wp
        epot_c( 352) = 0.479361755609e+04_wp
        epot_c( 353) = -0.547845404358e+04_wp
        epot_c( 354) = -0.789526384742e+02_wp
        epot_c( 355) = 0.184400748048e+05_wp
        epot_c( 356) = -0.143635356336e+05_wp
        epot_c( 357) = -0.129060844688e+05_wp
        epot_c( 358) = -0.405557959228e+05_wp
        epot_c( 359) = 0.393018703427e+04_wp
        epot_c( 360) = 0.734848646986e+03_wp
        epot_c( 361) = 0.104858540824e+05_wp
        epot_c( 362) = -0.606179549280e+04_wp
        epot_c( 363) = -0.239677089837e+05_wp
        epot_c( 364) = 0.124449293377e+05_wp
        epot_c( 365) = -0.565010495917e+04_wp
        epot_c( 366) = 0.429338995993e+04_wp
        epot_c( 367) = 0.523148566869e+04_wp
        epot_c( 368) = -0.911368721817e+04_wp
        epot_c( 369) = 0.564870179683e+04_wp
        epot_c( 370) = -0.502152609750e+04_wp
        epot_c( 371) = -0.465251526831e+04_wp
        epot_c( 372) = 0.222237537150e+04_wp
        epot_c( 373) = 0.112406687182e+05_wp
        epot_c( 374) = 0.539796908287e+04_wp
        epot_c( 375) = -0.690742754175e+03_wp
        epot_c( 376) = 0.603749679769e+04_wp
        epot_c( 377) = 0.799555982423e+04_wp
        epot_c( 378) = -0.113923260696e+05_wp
        epot_c( 379) = -0.112708572512e+05_wp
        epot_c( 380) = -0.805117918051e+04_wp
        epot_c( 381) = -0.390384250914e+04_wp
        epot_c( 382) = 0.193037973836e+04_wp
        epot_c( 383) = 0.334200703653e+03_wp
        epot_c( 384) = 0.314461704840e+04_wp
        epot_c( 385) = -0.104231726717e+05_wp
        epot_c( 386) = -0.415096092629e+04_wp
        epot_c( 387) = 0.615322194725e+03_wp
        epot_c( 388) = 0.193752788883e+03_wp
        epot_c( 389) = 0.972094902709e+04_wp
        epot_c( 390) = 0.457165203225e+04_wp
        epot_c( 391) = 0.200579473774e+04_wp
        epot_c( 392) = -0.570363045192e+04_wp
        epot_c( 393) = -0.189376029202e+04_wp
        epot_c( 394) = 0.141074483583e+04_wp
        epot_c( 395) = 0.569963181391e+04_wp
        epot_c( 396) = 0.335923390845e+05_wp
        epot_c( 397) = 0.980132445079e+04_wp
        epot_c( 398) = 0.115856281358e+04_wp
        epot_c( 399) = -0.471651079258e+04_wp
        epot_c( 400) = -0.104793710728e+05_wp
        epot_c( 401) = 0.552426474937e+04_wp
        epot_c( 402) = -0.262317260809e+05_wp
        epot_c( 403) = -0.234962207650e+04_wp
        epot_c( 404) = 0.101852776750e+04_wp
        epot_c( 405) = 0.804852499577e+04_wp
        epot_c( 406) = 0.454207213061e+04_wp
        epot_c( 407) = 0.636981421673e+04_wp
        epot_c( 408) = -0.283308858173e+04_wp
        epot_c( 409) = 0.320137033865e+04_wp
        epot_c( 410) = -0.372629403204e+04_wp
        epot_c( 411) = 0.456011869848e+05_wp
        epot_c( 412) = -0.563675461730e+04_wp
        epot_c( 413) = -0.530456535380e+03_wp
        epot_c( 414) = -0.118525831719e+05_wp
        epot_c( 415) = -0.380166689616e+04_wp
        epot_c( 416) = 0.516078752085e+04_wp
        epot_c( 417) = 0.337551470491e+04_wp
        epot_c( 418) = 0.327969347317e+04_wp
        epot_c( 419) = 0.169364835464e+03_wp
        epot_c( 420) = -0.187461896052e+05_wp
        epot_c( 421) = -0.130075907733e+04_wp
        epot_c( 422) = 0.913355294232e+02_wp
        epot_c( 423) = -0.521028908314e+04_wp
        epot_c( 424) = 0.126326318708e+05_wp
        epot_c( 425) = -0.523334675988e+04_wp
        epot_c( 426) = 0.130911936078e+05_wp
        epot_c( 427) = -0.120879055239e+05_wp
        epot_c( 428) = 0.761622252765e+04_wp
        epot_c( 429) = 0.603116293566e+04_wp
        epot_c( 430) = -0.339051917129e+04_wp
        v            = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                j     = 3*(iatom-1)+idir
                tx(j) = x(iatom,idir)
            enddo
        enddo
        tg = 0.e0_wp
        call o4pes(tx, v, tg, igrad)
        v  = v/electronvolt_to_kilocalorie_per_mole
        tg = tg/electronvolt_to_kilocalorie_per_mole
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
        d = 0.e0_wp

        contains

        !> @brief Evaluate o4pes for the O4_quintet potential.
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
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i, igrad, j, nob, k
            real(wp) :: v
            real(wp) :: dvdx(12), x(12)
            !      common /coord/    R(6)
            !      common /epot/     rMs(6),rM(0:111),P(0:305),C(430),B(430)
            !      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
            !     $                  dRdX(6,12),dBdR(6,276)
            !      common /msprmt/   a,ab,re
            ! Convert the Cartesian coordinate array to pair distances
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

        !> @brief Evaluate coord_convt for the O4_quintet potential.
        subroutine coord_convt( x )
            !  Program to calculate the six interatomic distance
            !  from the Cartesian coordinate array
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i
            real(wp) :: x(12)
            !  Now, calculate the inter-atomic distance
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

        !> @brief Evaluate EvV for the O4_quintet potential.
        subroutine evv( v )
            ! Subroutine to evaluate V for given R
            ! V(R) = C*P
            ! C:            Coefficients, stored in 'dim.inc'
            ! P:            Basis functions evaluated for given R
            ! rMs:          rMs(6), six mixed exponential-Gaussian terms (MEG)
            ! a:            Nonlinear parameters in Morse terms(Angstrom)
            ! ab:           Nonlinear parameters in Gauss terms(Angstrom^2)
            ! re:           Equilibrium bond length(Angstrom)
            ! nop:          number of points
            ! nom:          number of monomials
            ! nob:          number of basis functions(polynomials)
            ! rM(0:111):    Array to store monomials
            ! P(0:305):     Array to store polynomials
            ! B(1:276):     Array to store basis functions
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i, j, k
            real(wp) :: dist, dv2dr, v, v2, disp, dispdr(6)
            ! Calculate the six MEG terms for each point
            call evmorse
            ! Calculate the monomials for each point by using six MEG terms
            call evmono
            ! Calculate the polynomials (basis functions) by using monomials
            call evpoly
            ! Calculate the basis functions by removing unconnected and 2-body terms
            call evbas
            ! Initialized v to be 2De 2*120.243 kcal/mol
            v = 240.486e0_wp
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

        !> @brief Evaluate EvdVdX for the O4_quintet potential.
        subroutine evdvdx( x, dvdx )
            use chempotpy_geometry, only: transform_internal_gradient_to_cartesian
            ! Subroutine to evaluate dRdX for given R and X
            ! R:            R(6), 6 bond lengths
            ! X:            X(12), 12 Cartesian coordinates
            ! rM(0:111):    Array to store monomials
            ! P(0:305):     Array to store polynomials
            ! dVdX:         dVdX(12), derivatives of V w.r.t. Cartesian coordinates
            ! dVdR:         dVdR(6), derivatives of V w.r.t.6 bond length
            ! dRdX:         dRdX(6,12), derivatives of R(6) w.r.t. 12
            !               Cartesian coordinates
            use chempotpy_kinds, only: wp
            implicit none
            real(wp), intent(in)  :: x(12)
            real(wp), intent(out) :: dvdx(12)
            call evdvdr
            call evdrdx(x)
            call transform_internal_gradient_to_cartesian(gradt_dvdr, gradt_drdx, dvdx)
        end subroutine evdvdx

        !> @brief Evaluate EvMorse for the O4_quintet potential.
        subroutine evmorse
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i
            ! mixed exponential-Gaussian terms = exp(-(r-re)/a-(r-re)^2/ab)
            ! re:   equlibrium bond length(1.208A)
            ! a:    nonlinear parameter optimized for diatomic molecules, unit
            ! Anstrom
            ! ab:   nonlinear parameter optimized for diatomic molecules, unit
            ! Anstrom^2
            do i=1, 6
                epot_rms(i) = exp(-(coord_r(i)-msprmt_re)/msprmt_a-((coord_r(i)-msprmt_re)**2.0e0_wp)/msprmt_ab)
            enddo
        end subroutine evmorse

        !> @brief Evaluate EvMono for the O4_quintet potential.
        subroutine evmono
            !  The subroutine accepts six MEG variables (X) and calculates the
            !  monomials(M) that do not have usable decomposition.
            !  For A4 with max. degree 10, the number of monomials is nom.
            use chempotpy_kinds, only: wp
            implicit none
            epot_rm(0)   = 1.0e0_wp
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

        !> @brief Evaluate EvPoly for the O4_quintet potential.
        subroutine evpoly
            !  The subroutine accepts monomials (m) and calculates the
            !  permutation invariant polynomials(p)
            !  For A4 with max. degree 10, the number of polynomials is nob.
            use chempotpy_kinds, only: wp
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

        !> @brief Evaluate evbas for the O4_quintet potential.
        subroutine evbas
            !  The subroutine eliminate the 2-body terms in Bowman's approach
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i
            real(wp) :: b1(466)
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

        !> @brief Evaluate ev2gm2 for the O4_quintet potential.
        subroutine ev2gm2( r, v, grad, imol, igrad )
            ! Compute the diatomic potential of ground-state triplet O2
            ! References: J. Chem. Phys. 132, 074307 (2010)
            ! Input:  r      interatomic distance in Angstrom
            ! Output: V      potential in kcal/mol
            !         grad   gradient (kcal/mol)/Angstrom
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: imol
            real(wp) :: r
            real(wp) :: v, grad
            ! Parameters of analytical even-tempered Gaussian expansions for the
            ! ground
            ! state potential energy curve of O2 CBS+SR+SO+CV. Units: alpha in
            ! Anstromgs^-2, beta=dimensionless, a_k in millihartree.
            real(wp) :: alpha, beta, a(0:7)
            integer  :: k, igrad
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
            alpha = 9.439784362354936e-1_wp
            beta  = 1.262242998506810e0_wp
            a(0)  = -1.488979427684798e3_wp
            a(1)  = 1.881435846488955e4_wp
            a(2)  = -1.053475425838226e5_wp
            a(3)  = 2.755135591229064e5_wp
            a(4)  = -4.277588997761775e5_wp
            a(5)  = 4.404104009614092e5_wp
            a(6)  = -2.946204062950765e5_wp
            a(7)  = 1.176861219078620e5_wp
            v     = 0.e0_wp
            do k=0, 7
                v = v + a(k)*exp(-alpha*beta**k*r**2)
            enddo
            ! From millihartree to kcal/mol
            v = v*627.509523475149e-3_wp
            ! Compute the gradient if needed
            if( igrad .eq. 1 )then
                grad = 0.e0_wp
                do k=0, 7
                    grad = grad-2.e0_wp*a(k)*alpha*beta**k*r*exp(-alpha*beta**k*r**2)
                enddo
                ! Convert from millihartree/A to i(kcal/mol)/A
                grad = grad*627.509523475149e-3_wp
            endif
        end subroutine ev2gm2

        !> @brief Evaluate EvdVdR for the O4_quintet potential.
        subroutine evdvdr
            ! Subroutine to evaluate dVdR for given R
            ! dVdR = dV2dR + C*dBdR
            ! C:            Coefficients, stored in 'dim.inc'
            ! P:            Basis functions evaluated for given R
            ! M:            Monomials evaluated for given R
            ! dV2dR:        Gradient of 2-body interactions
            ! dMsdR:        dMsdR(6,6), 6 MEG terms w.r.t. 6 bond length
            ! dMdR:         dMdR(6,nom), nom monomials w.r.t.6 bond length
            ! dPdR:         dPdR(6,nob), nop polynomial basis functions
            !               w.r.t. 6 bond length
            ! nom:          number of monomials
            ! nob:          number of basis functions(polynomials)
            ! M(nom):       Array to store monomials
            ! P(nob):       Array to store polynomials
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i, j
            real(wp) :: dist, v2, dv2dr, disp, dispdr(6)
            ! Initialize dVdR(6)
            do i=1, 6
                gradt_dvdr(i) = 0.0e0_wp
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

        !> @brief Evaluate EvdRdX for the O4_quintet potential.
        subroutine evdrdx( x )
            use chempotpy_geometry, only: all_pair_distance_jacobian
            use chempotpy_kinds,    only: wp
            implicit none
            real(wp), intent(in) :: x(12)
            call all_pair_distance_jacobian(x, coord_r, gradt_drdx)
        end subroutine evdrdx

        !> @brief Evaluate EvdMsdR for the O4_quintet potential.
        subroutine evdmsdr
            ! Subroutine to evaluate the derivatives of MEG term X
            ! w.r.t. interatomic distance R(6)
            ! dmsdR:        Local variables, dirm(6,6)
            ! a:            Nonlinear pamameter(Angstrom)
            ! ab:       Nonlinear pamameter(Angstrom^2)
            ! re:           equilibrium bond length(Angstrom)
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, j
            ! Initialize dmsdr
            do i=1, 6
                do j=1, 6
                    gradt_dmsdr(i,j) = 0.0e0_wp
                enddo
            enddo
            ! MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
            ! dmsdr(i,j)=0  i!=j
            do i=1, 6
                gradt_dmsdr(i,i) = (-2.0e0_wp*(coord_r(i)-msprmt_re)/msprmt_ab-1/msprmt_a)* exp(-(coord_r(i)- &
                    & msprmt_re)/ msprmt_a-((coord_r(i)-msprmt_re)**2.0e0_wp)/msprmt_ab)
            enddo
        end subroutine evdmsdr

        !> @brief Evaluate EvdMdR for the O4_quintet potential.
        subroutine evdmdr
            !  The subroutine accepts M(nom) and dMSdR(6,6) and calculates the
            !  dMdR(6,nom) that do not have usable decomposition.
            !  For A4 with max. degree 10, the number of monomials is nom.
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i
            do i=1, 6
                gradt_dmdr(i,0)   = 0.0e0_wp
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

        !> @brief Evaluate EvdPdr for the O4_quintet potential.
        subroutine evdpdr
            !  The subroutine accepts monomials (m) and calculates the
            !  permutation invariant polynomials(p)
            !  For A4 with max. degree 10, the number of polynomials is nob.
            use chempotpy_kinds, only: wp
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

        !> @brief Evaluate evdbdr for the O4_quintet potential.
        subroutine evdbdr
            !  The subroutine eliminate the 2-body terms in Bowman's approach
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: j
            integer  :: i
            real(wp) :: db1dr(6,466)
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

        !> @brief Evaluate d3disp for the O4_quintet potential.
        subroutine d3disp( dist, disp, dispdr, igrad )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: cn(4), s6, s8, rs6, rs8
            real(wp) :: dist(6), e6(6), e8(6), disp, dispdr(6), c6(6)
            real(wp) :: e6dr(6), e8dr(6)
            integer  :: iz(4), mxc(94), i, j, igrad
            real(wp) :: c6ab(94,94,5,5,3)
            real(wp) :: r2r4(94)
            real(wp) :: autoang, autokcal
            autoang  = 0.52917726e0_wp
            autokcal = 627.509541e0_wp
            ! Pragya's generalized parameters for BJ damping
            s6  = 1.0e0_wp
            s8  = 2.0e0_wp
            rs6 = 0.5299e0_wp
            rs8 = 2.20e0_wp
            do i=1, 6
                dist(i) = dist(i)/autoang
            enddo
            ! iz for O4 system
            iz(1) = 8
            iz(2) = 8
            iz(3) = 8
            iz(4) = 8
            ! C6 for O4 system
            c6(1) = 12.8e0_wp
            c6(2) = 12.8e0_wp
            c6(3) = 12.8e0_wp
            c6(4) = 12.8e0_wp
            c6(5) = 12.8e0_wp
            c6(6) = 12.8e0_wp
            ! Calculate dispersion correction
            call edisp(94, 5, 4, dist, iz, mxc, &
                & rs6, rs8, e6, e8, e6dr, e8dr, c6, 0)
            disp = 0.0e0_wp
            do i=1, 6
                disp = disp + (-s6*e6(i)-s8*e8(i))*autokcal
            enddo
            if( igrad .eq. 1 )then
                call edisp(94, 5, 4, dist, iz, mxc, &
                    & rs6, rs8, e6, e8, e6dr, e8dr, c6, 1)
                dispdr(:) = 0.0e0_wp
                do i=1, 6
                    dispdr(i) = dispdr(i) + (-s6*e6dr(i)-s8*e8dr(i))*autokcal/autoang
                enddo
            endif
            do i=1, 6
                dist(i) = dist(i)*autoang
            enddo
        end subroutine d3disp

        !> @brief Evaluate edisp for the O4_quintet potential.
        subroutine edisp( max_elem, maxc, n, dist, iz, mxc, &
                & rs6, rs8, e6, e8, e6dr, e8dr, c6a, igrad )
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: n, iz(4), max_elem, maxc, mxc(max_elem)
            real(wp) :: dist(6), r2r4(max_elem), r0ab(max_elem,max_elem)
            real(wp) :: rs6, rs8, rcov(max_elem)
            real(wp) :: c6ab(max_elem,max_elem,maxc,maxc,3)
            real(wp) :: e6(6), e8(6), c6a(6), e6dr(6), e8dr(6)
            integer  :: iat, jat, igrad
            real(wp) :: r, tmp, c6, c8, a1, a2
            real(wp) :: damp6, damp8
            real(wp) :: cn(n)
            real(wp) :: r2ab(n*n), cc6ab(n*n), dmp(n*n)
            integer  :: step
            e6(:)   = 0.0e0_wp
            e8(:)   = 0.0e0_wp
            e6dr(:) = 0.0e0_wp
            e8dr(:) = 0.0e0_wp
            a1      = rs6
            a2      = rs8
            !  r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
            !  the large number of digits is just to keep the results consistent
            !  with older versions. They should not imply any higher accuracy than
            !  the old values
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
                    c8 = 3.0e0_wp*c6*r2r4(iz(iat))*r2r4(iz(jat))
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
end module chempotpy_o4_quintet_pes

!> @brief Expose O4_quintet through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_o4_quintet_pes, only: evaluate
    use chempotpy_kinds,          only: wp
    implicit none
    integer, parameter    :: natoms  = 4
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
