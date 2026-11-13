!   System:                        H3
!   Functional form:
!   Common name:                   A3
!   Number of derivatives:         0
!   Number of bodies:              3
!   Number of electronic surfaces: 1
!   Interface:                     3-1V
!   Data file:
!   References:      SL Mielke, BC Garrett, and KA Peterson, J. Chem. Ph
!   Notes:    3 H-H distances passed in bohr.
!             Energy returned in Hatrees relative to the classical minim

!> @brief Implement H3_GEN_A3_2002, a potential energy surface for H3.
module chempotpy_h3_gen_a3_2002_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_A3_2002 potential energy surface for H3.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,              only: energy_to_electronvolt
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
        integer               :: istate
        real(wp)              :: e(1), r(1,3), v
        integer               :: iatom, idir, i, j
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        ! input cartesian is HHH
        call triatomic_distances_to_bohr(x, r(1,:))
        call pot(r, e, 1)
        v = energy_to_electronvolt(e(1))
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = v
            enddo
        else
            write (*, *) 'Only energy is available'
        endif
    end subroutine evaluate
end module chempotpy_h3_gen_a3_2002_pes

!> @brief Expose H3_GEN_A3_2002 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_a3_2002_pes, only: evaluate
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

module a3data
    use chempotpy_h3_2002_basis, only: primary_term_count, secondary_term_count
    use chempotpy_kinds,         only: wp
    implicit none
    private
    public :: ald, beta, betad, epscem, epslon, eshift, maxi2a, maxi2b, rnought
    public :: xlin1, xlin2, xlin3, xlin4
    real(wp), parameter     :: ald     = 0.02e0_wp
    real(wp), parameter     :: betad   = 0.72e0_wp
    real(wp), parameter     :: epscem  = 1.e-12_wp
    real(wp), parameter     :: epslon  = 1.e-12_wp
    integer, parameter      :: maxi2a  = primary_term_count
    integer, parameter      :: maxi2b  = secondary_term_count
    real(wp), parameter     :: eshift  = 0.1729928728793815e0_wp
    real(wp), parameter     :: beta    = 0.9291969140819419e0_wp
    real(wp), parameter     :: rnought = 2.e0_wp
    real(wp), dimension(89) :: xlin1=(/-0.1523291433408970e0_wp, &
        & 0.8814508465229537e0_wp, 6.537053395327828e0_wp, &
        & -11.34209491672260e0_wp, 2.504271304238934e0_wp, &
        & 4.326668451461330e0_wp, -2.239606610799320e0_wp, &
        & -0.4128984024042961e0_wp, 0.5186887309355593e0_wp, &
        & -0.1196607910913041e0_wp, 8.6003994509478154e-3_wp, &
        & 41.15414291845340e0_wp, 86.45559795459729e0_wp, &
        & -116.3713616382037e0_wp, 12.58972597209381e0_wp, &
        & 33.81929163465899e0_wp, -17.15593990627275e0_wp, &
        & 4.196478522411963e0_wp, -0.3411919636902630e0_wp, &
        & -3.4206236345359939e-3_wp, -127.6850556787161e0_wp, &
        & 12.44205218121803e0_wp, -46.65469830727286e0_wp, &
        & 2.132190766110241e0_wp, 5.029968337261709e0_wp, &
        & -0.6386010902854538e0_wp, 5.0081780443521871e-3_wp, &
        & 32.66461240039222e0_wp, -12.69768482032793e0_wp, &
        & 0.8275601413171687e0_wp, 3.9468614677319450e-2_wp, &
        & 1.5477194647200569e-2_wp, 1.089948589115662e0_wp, &
        & 0.7206200046869276e0_wp, -7.8425404380267669e-2_wp, &
        & -0.1603878592389938e0_wp, -9.724970811924141e0_wp, &
        & 135.5222874536157e0_wp, 270.0021723964958e0_wp, &
        & -175.4825227217210e0_wp, -58.58906289532513e0_wp, &
        & 100.5430564127615e0_wp, -37.00315174873428e0_wp, &
        & 6.447449239215345e0_wp, 8.2570662069006429e-3_wp, &
        & -7.1098331653255367e-2_wp, 715.3552761780156e0_wp, &
        & 703.9102525501766e0_wp, -1861.078911265242e0_wp, &
        & -125.9587491921284e0_wp, 225.8731414579132e0_wp, &
        & -46.07417109551250e0_wp, 1.490246977236569e0_wp, &
        & 0.2982638988343370e0_wp, -999.4775422429434e0_wp, &
        & 151.7862715652945e0_wp, -133.0715723394372e0_wp, &
        & -12.30680794373364e0_wp, 5.825782052452050e0_wp, &
        & -0.2439712905709224e0_wp, 85.16524762327094e0_wp, &
        & 16.46204547051492e0_wp, -7.133858657001830e0_wp, &
        & 0.4030617759680488e0_wp, 0.6909270689548164e0_wp, &
        & 0.3097014450084296e0_wp, 2556.785909933932e0_wp, &
        & 5103.315112906826e0_wp, -2912.860957206383e0_wp, &
        & -160.8805563304034e0_wp, 173.7542607534209e0_wp, &
        & -7.478371849288498e0_wp, -1.555852990003698e0_wp, &
        & 2468.715016248928e0_wp, 529.4607240124748e0_wp, &
        & -162.6956551745768e0_wp, -4.113114408714173e0_wp, &
        & 1.611041727204823e0_wp, 56.28848656173622e0_wp, &
        & 9.735269264701914e0_wp, -3.377129269530819e0_wp, &
        & 2.600888740672099e0_wp, -3246.118792587758e0_wp, &
        & 161.4154798439835e0_wp, 14.65512118402112e0_wp, &
        & 1.684352485791373e0_wp, -19.56031027772914e0_wp, &
        & -8.9132714864690765e-2_wp, -0.6919443797821230e0_wp/)
    real(wp), dimension(89) :: xlin2=(/2.6159085145117750e-2_wp, &
        & 1.2241854330557397e-2_wp, -0.6763682764316550e0_wp, &
        & -0.1562836764561804e0_wp, 1.996594730253928e0_wp, &
        & -1.707513510844818e0_wp, 0.2797466754569802e0_wp, &
        & 0.2576622535966453e0_wp, -0.1404344725580116e0_wp, &
        & 2.6052658461384465e-2_wp, -1.6693444722953859e-3_wp, &
        & -0.4833014947701142e0_wp, -11.17866172577437e0_wp, &
        & 7.607092641774579e0_wp, 9.519929742191135e0_wp, &
        & -12.37037283407236e0_wp, 4.707623579732490e0_wp, &
        & -0.5730444600084309e0_wp, -3.4653598431171197e-2_wp, &
        & 7.8187537452263211e-3_wp, 38.72674320125338e0_wp, &
        & -7.269232169491375e0_wp, -9.879886268080352e0_wp, &
        & 8.517458481793820e0_wp, -0.8199496374616673e0_wp, &
        & -0.1559314546076156e0_wp, 1.5401566524384906e-2_wp, &
        & 9.554678972364389e0_wp, 3.201632506258832e0_wp, &
        & -2.321941631684996e0_wp, 3.6930680231572893e-2_wp, &
        & 3.2742627057071150e-2_wp, 1.900429427415868e0_wp, &
        & -0.3327189910649105e0_wp, -7.7188825506846991e-3_wp, &
        & 9.3760580872625787e-2_wp, -1.014102607634918e0_wp, &
        & 22.02809848224747e0_wp, -24.76473768346871e0_wp, &
        & -37.02983718150585e0_wp, 58.80030375084864e0_wp, &
        & -28.47495308394478e0_wp, 5.757365917529928e0_wp, &
        & 4.3246175682466526e-2_wp, -0.1676517786924019e0_wp, &
        & 1.5097351831581923e-2_wp, 135.0679946550904e0_wp, &
        & -213.7000109608329e0_wp, 95.89118538821953e0_wp, &
        & 25.40918010829430e0_wp, 9.782035588385490e0_wp, &
        & -4.168396476529826e0_wp, 0.6883577777053669e0_wp, &
        & -6.9367683522399748e-2_wp, -262.5899913042082e0_wp, &
        & -53.22977186145765e0_wp, 6.121343202849875e0_wp, &
        & 7.823177388578925e0_wp, -0.2711421320419219e0_wp, &
        & -0.1043327957891995e0_wp, 53.77598190726489e0_wp, &
        & -12.01741605017266e0_wp, -1.4784757119912589e-2_wp, &
        & 0.1488112444684154e0_wp, 2.732678494168803e0_wp, &
        & -0.1970131259131940e0_wp, -467.4973230906062e0_wp, &
        & -444.5887032753010e0_wp, -315.4834996683269e0_wp, &
        & 48.67913126874383e0_wp, 24.84508509899326e0_wp, &
        & -4.729140122994481e0_wp, 0.2069558855324623e0_wp, &
        & 781.9447144772034e0_wp, -143.3043797424338e0_wp, &
        & -28.68693796003896e0_wp, 2.344077954027283e0_wp, &
        & 0.4918333516655770e0_wp, 72.15250134833730e0_wp, &
        & -2.105811641390979e0_wp, -0.8842156350073892e0_wp, &
        & 0.7824630993546668e0_wp, 303.3073086352659e0_wp, &
        & -32.34761963775958e0_wp, -2.157141352189826e0_wp, &
        & -0.8000872307506535e0_wp, 5.610906605940581e0_wp, &
        & 1.199667001896809e0_wp, -4.678892380269374e0_wp/)
    real(wp), dimension(71) :: xlin3=(/5.3353477682006475e-2_wp, &
        & -0.1599144960966421e0_wp, -4.312643718679910e0_wp, &
        & 6.814315854917867e0_wp, -2.205456828188728e0_wp, &
        & -0.7143749272312716e0_wp, -0.4114395291526620e0_wp, &
        & 0.9520707372222905e0_wp, -0.3745853003354750e0_wp, &
        & 4.2811012132827980e-2_wp, -21.36048998715569e0_wp, &
        & -74.47880934502300e0_wp, 59.54491453184337e0_wp, &
        & 26.38093015799090e0_wp, -40.49367250732830e0_wp, &
        & 15.47625249844092e0_wp, -3.376506055551386e0_wp, &
        & 0.1062027161738525e0_wp, 18.81967136271011e0_wp, &
        & 71.16187431300146e0_wp, 14.71335233143087e0_wp, &
        & 6.721378533743240e0_wp, -3.339779478411798e0_wp, &
        & -0.2339343234862702e0_wp, -72.80454829701864e0_wp, &
        & 30.21027407834825e0_wp, -1.905603230726993e0_wp, &
        & -0.3996844993912013e0_wp, 2.076285530599899e0_wp, &
        & -1.070380536340887e0_wp, 8.598507152662780e0_wp, &
        & -73.61703213544804e0_wp, -219.8778318595011e0_wp, &
        & 48.95327228138328e0_wp, 109.3741906195230e0_wp, &
        & -104.0415972823556e0_wp, 35.11211365107425e0_wp, &
        & -6.534860269104033e0_wp, -8.0422603679405902e-3_wp, &
        & -600.3417941867114e0_wp, -878.4454407436597e0_wp, &
        & 1243.686176699445e0_wp, 591.3052045918021e0_wp, &
        & -243.0727894093502e0_wp, 39.64018613563970e0_wp, &
        & -0.8090321241234495e0_wp, 753.6354030176653e0_wp, &
        & 592.5340756904949e0_wp, 109.1519077002292e0_wp, &
        & 3.933746884407120e0_wp, 3.506148188290761e0_wp, &
        & -189.3502887433559e0_wp, 13.94268601946511e0_wp, &
        & 0.6424807674312493e0_wp, 1.419771258623844e0_wp, &
        & -2050.644716844225e0_wp, -5080.398907091004e0_wp, &
        & 2822.341530025133e0_wp, 638.1501839778731e0_wp, &
        & -217.6432578060511e0_wp, 9.608610580227475e0_wp, &
        & -4953.965720751209e0_wp, -460.4461493597798e0_wp, &
        & 40.61534110405616e0_wp, -18.26218983170062e0_wp, &
        & -113.4926963135811e0_wp, 3.382375857465263e0_wp, &
        & 2169.169337193152e0_wp, 154.5462714936868e0_wp, &
        & 84.83291904037159e0_wp, -71.12679286953023e0_wp/)
    real(wp), dimension(71) :: xlin4=(/-9.4876856817530143e-3_wp, &
        & -3.6483362893877616e-2_wp, 0.4205707974558235e0_wp, &
        & -0.1059572058101154e0_wp, -0.3707179832723506e0_wp, &
        & -0.1445341683023923e0_wp, 0.5958052755316040e0_wp, &
        & -0.3854655745032128e0_wp, 0.1010458670628552e0_wp, &
        & -9.4255692193045623e-3_wp, 0.1509885666028168e0_wp, &
        & 7.682780170568845e0_wp, -0.6840277749251158e0_wp, &
        & -11.98139522620621e0_wp, 9.946678161114002e0_wp, &
        & -2.590064247380303e0_wp, 7.6410732889741748e-2_wp, &
        & 4.3207950587928659e-2_wp, -15.82021851395566e0_wp, &
        & -11.55029253337579e0_wp, 9.357730866212668e0_wp, &
        & -1.817290225301895e0_wp, -2.051176763692164e0_wp, &
        & 0.3422211063996889e0_wp, 5.935560909793717e0_wp, &
        & -5.899083969710368e0_wp, -0.3644702094335688e0_wp, &
        & 0.5301830356858073e0_wp, -2.119504022080795e0_wp, &
        & 0.4144161928859211e0_wp, 0.2636535169916772e0_wp, &
        & -9.256727433362473e0_wp, 1.010605951519258e0_wp, &
        & 48.54597955205001e0_wp, -48.87059486081539e0_wp, &
        & 17.98526173595513e0_wp, -2.237557599107648e0_wp, &
        & -0.3789995719107901e0_wp, 8.6650123871130774e-2_wp, &
        & -119.5839332405538e0_wp, 133.3734002478926e0_wp, &
        & -28.89035273972550e0_wp, -56.60471823462028e0_wp, &
        & -9.1459472557995802e-2_wp, -2.754447481273484e0_wp, &
        & 9.6332311722992170e-2_wp, 253.4638121163468e0_wp, &
        & 112.6549501199819e0_wp, -2.630187545871118e0_wp, &
        & -11.07273872864548e0_wp, 0.3646509060305644e0_wp, &
        & -68.65504615698455e0_wp, 10.86627678597302e0_wp, &
        & -0.6358917122975185e0_wp, -2.352160870433293e0_wp, &
        & 84.30979852066358e0_wp, 797.2622595018731e0_wp, &
        & 214.5994942620052e0_wp, 30.59962513407079e0_wp, &
        & -31.32236829922999e0_wp, 2.833941401102188e0_wp, &
        & -476.2563698124825e0_wp, 122.6747471786113e0_wp, &
        & 17.34552327162847e0_wp, -0.8126973285879644e0_wp, &
        & -38.30796451594332e0_wp, 1.661856665041447e0_wp, &
        & -1019.748561415903e0_wp, 64.84858716655161e0_wp, &
        & 5.856777175810681e0_wp, -8.373832093263392e0_wp/)
end module

!> @brief Evaluate a batch of H3 A3 geometries.
subroutine pot( rvp, evp, nt )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)   :: nt
    real(wp), intent(in)  :: rvp(nt,3)
    real(wp), intent(out) :: evp(nt)
    integer               :: ivp
    real(wp)              :: r1x, r2x, r3x
    do ivp = 1, nt
        r1x = rvp(ivp,1)
        r2x = rvp(ivp,2)
        r3x = rvp(ivp,3)
        call pota3(r1x, r2x, r3x, evp(ivp))
    enddo
end subroutine pot

!> @brief Evaluate one H3 A3 geometry.
subroutine pota3( r1x, r2x, r3x, eee )
    use a3data,                  only: ald, beta, betad, epscem, epslon, eshift, maxi2a, maxi2b, rnought, xlin1, xlin2, xlin3, xlin4
    use chempotpy_h3_2002_basis, only: get_h3_basis_term
    use chempotpy_kinds,         only: wp
    implicit none
    real(wp), intent(in)  :: r1x, r2x, r3x
    real(wp), intent(out) :: eee
    real(wp)              :: bb, damper, es1, es2, es3, etrip1, etrip2, etrip3, q1, q2, q3, rdamp1, rdamp2, &
        & rdamp3, rho1, rho2, rho3, rprod, sum1, sum2, sum3, sum4, v, vlondon, warp, xj, xj1, xj2, xj3, xjs
    real(wp) :: rho1t(0:12), rho2t(0:12), rho3t(0:12)
    integer  :: coefficient_index, first_power, ii, it, second_power, third_power
    call trip(r1x, etrip1)
    call trip(r2x, etrip2)
    call trip(r3x, etrip3)
    call singlet(r1x, es1)
    call singlet(r2x, es2)
    call singlet(r3x, es3)
    q1      = 0.5e0_wp*(es1+etrip1)
    q2      = 0.5e0_wp*(es2+etrip2)
    q3      = 0.5e0_wp*(es3+etrip3)
    xj1     = 0.5e0_wp*(es1-etrip1)
    xj2     = 0.5e0_wp*(es2-etrip2)
    xj3     = 0.5e0_wp*(es3-etrip3)
    xj      = sqrt(epslon+0.5e0_wp*((xj3-xj1)**2+(xj2-xj1)**2+(xj3-xj2)**2))
    xjs     = 0.25e0_wp*(xj1+xj2+xj3)
    vlondon = q1+q2+q3-xj
    rho1    = exp(-beta*(r1x-rnought))
    rho2    = exp(-beta*(r2x-rnought))
    rho3    = exp(-beta*(r3x-rnought))
    if( (r1x .le. betad).or.(r2x .le. betad).or.(r3x .le. betad) )then
        damper = 0.e0_wp
    else
        rdamp1 = ald/(betad-r1x)
        rdamp2 = ald/(betad-r2x)
        rdamp3 = ald/(betad-r3x)
        damper = exp(rdamp1+rdamp2+rdamp3)
    endif
    bb       = sqrt(epscem+(r1x-r3x)**2+(r1x-r2x)**2+(r2x-r3x)**2)
    warp     = 1.e0_wp/(1.e0_wp/r1x+1.e0_wp/r2x+1.e0_wp/r3x)
    v        = eshift+vlondon
    sum1     = 0.e0_wp
    sum2     = 0.e0_wp
    sum3     = 0.e0_wp
    sum4     = 0.e0_wp
    rho1t(0) = 1.e0_wp
    rho2t(0) = 1.e0_wp
    rho3t(0) = 1.e0_wp
    do it=1, 12
        rho1t(it) = rho1t(it-1)*rho1
        rho2t(it) = rho2t(it-1)*rho2
        rho3t(it) = rho3t(it-1)*rho3
    enddo
    do ii=1, maxi2a
        call get_h3_basis_term(ii, first_power, second_power, third_power, coefficient_index)
        rprod = rho1t(third_power)*rho2t(second_power)*rho3t(first_power)
        sum1  = sum1+xlin1(coefficient_index)*rprod
        sum2  = sum2+xlin2(coefficient_index)*rprod
    enddo
    do ii=maxi2a+1, maxi2a+maxi2b
        call get_h3_basis_term(ii, first_power, second_power, third_power, coefficient_index)
        rprod = rho1t(third_power)*rho2t(second_power)*rho3t(first_power)
        sum3  = sum3+xlin3(coefficient_index)*rprod
        sum4  = sum4+xlin4(coefficient_index)*rprod
    enddo
    eee = v+damper*(sum1+sum2*bb+sum3*warp+sum4*bb*warp)
end subroutine pota3

!> @brief Evaluate trip for the H3_GEN_A3_2002 potential.
subroutine trip( r, v )
    !     H2 triplet curve for FCI/aug-cc-pVTZ
    use chempotpy_kinds, only: wp
    implicit none
    integer                 :: j
    real(wp)                :: r, v, damp, xd, prefac
    real(wp)                :: c6=-6.499027e0_wp
    real(wp)                :: c8=-124.3991e0_wp
    real(wp)                :: c10=-3285.828e0_wp
    real(wp)                :: beta=0.2e0_wp
    real(wp)                :: re=1.401e0_wp
    real(wp)                :: alpha=3.8716418156250960e0_wp
    real(wp), dimension(17) :: xlp=(/2.1948631629921860e-1_wp, &
        & 5.8892583097268760e-1_wp, 8.6235600154107670e-1_wp, &
        & 8.6194702003339070e-1_wp, 6.8604820043881530e-1_wp, &
        & 4.5436692359636230e-1_wp, 2.3371369821182240e-1_wp, &
        & 6.9843424078195320e-2_wp, -2.6494129220450960e-3_wp, &
        & -8.5475462878786820e-3_wp, 3.0498779190429310e-3_wp, &
        & 3.0448279356991980e-3_wp, -5.6431515849655190e-4_wp, &
        & -4.0765667872606770e-4_wp, 1.9838486248362640e-4_wp, &
        & -3.2129748801006750e-5_wp, 1.9377307867602390e-6_wp/)
    damp   = 1.e0_wp-exp(-beta*r**2)
    xd     = damp/r
    v      = c6*xd**6+c8*xd**8+c10*xd**10
    prefac = exp(- alpha*(r-re))
    if( (r-re) .ne. 0.e0_wp )then
        do j=1, 17
            v = v+xlp(j)*(r-re)**(j-1)*prefac
        enddo
    else
        v = v+xlp(1)*prefac
    endif
end subroutine trip

!> @brief Evaluate singlet for the H3_GEN_A3_2002 potential.
subroutine singlet( r, v )
    !     H2 singlet potential for FCI/aug-cc-pVTZ
    use chempotpy_kinds, only: wp
    implicit none
    integer                 :: j
    real(wp)                :: r, v, damp, xd, prefac
    real(wp)                :: c6=-6.499027e0_wp
    real(wp)                :: c8=-124.3991e0_wp
    real(wp)                :: c10=-3285.828e0_wp
    real(wp)                :: beta=0.2e0_wp
    real(wp)                :: re=1.401e0_wp
    real(wp)                :: alpha=  3.9259981290585980e0_wp
    real(wp), dimension(17) :: xlp= (/-1.6948148539493060e-1_wp, &
        & -6.5341157509085330e-1_wp, -1.0568246834574850e0_wp, &
        & -1.0380175566115990e0_wp, -6.7657643586978360e-1_wp, &
        & -3.0434043205949370e-1_wp, -1.0305186076403050e-1_wp, &
        & -4.9128168199975780e-2_wp, -3.7985765587970510e-2_wp, &
        & -1.8881471279090120e-2_wp, 4.4201285958012090e-4_wp, &
        & 3.2421117771421200e-3_wp, -3.9290175622769260e-4_wp, &
        & -5.7935611043841150e-4_wp, 2.0970556978282260e-4_wp, &
        & -2.8399222117118230e-5_wp, 1.2735256346272810e-6_wp/)
    damp   = 1.e0_wp-exp(-beta*r**2)
    xd     = damp/r
    prefac = exp(- alpha*(r-re))
    v      = c6*xd**6+c8*xd**8+c10*xd**10
    if( (r-re) .ne. 0.e0_wp )then
        do j=1, 17
            v = v+xlp(j)*(r-re)**(j-1)*prefac
        enddo
    else
        v = v+xlp(1)*prefac
    endif
end subroutine singlet
