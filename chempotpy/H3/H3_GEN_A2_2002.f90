!   System:                        H3
!   Functional form:
!   Common name:                   a2
!   Number of derivatives:         0
!   Number of bodies:              3
!   Number of electronic surfaces: 1
!   Interface:                     3-1V
!   Data file:
!   References:      SL Mielke, BC Garrett, and KA Peterson, J. Chem. Ph
!   Notes:    3 H-H distances passed in bohr.
!             Energy returned in Hatrees relative to the classical minim

!> @brief Implement H3_GEN_A2_2002, a potential energy surface for H3.
module chempotpy_h3_gen_a2_2002_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_A2_2002 potential energy surface for H3.
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
end module chempotpy_h3_gen_a2_2002_pes

!> @brief Expose H3_GEN_A2_2002 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_a2_2002_pes, only: evaluate
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

module a2data
    use chempotpy_h3_2002_basis, only: primary_term_count, secondary_term_count
    use chempotpy_kinds,         only: wp
    implicit none
    private
    public :: ald, beta, betad, epscem, epslon, eshift, maxi2a, maxi2b, rnought
    public :: xlin1, xlin2, xlin3, xlin4
    real(wp), parameter :: ald     = 0.02e0_wp
    real(wp), parameter :: betad   = 0.72e0_wp
    real(wp), parameter :: epscem  = 1.e-12_wp
    real(wp), parameter :: epslon  = 1.e-12_wp
    integer, parameter  :: maxi2a  = primary_term_count
    integer, parameter  :: maxi2b  = secondary_term_count
    real(wp), parameter :: eshift  = 0.1662307680057728e0_wp
    real(wp), parameter :: beta    = 0.900e0_wp
    real(wp), parameter :: rnought = 2.e0_wp
    real(wp)            :: xlin1(89)=(/0.3122096268380428e0_wp, &
        & -0.5091675305404870e0_wp, 17.09937094766591e0_wp, &
        & -47.55205835431389e0_wp, 58.80549483490758e0_wp, &
        & -38.83145155729753e0_wp, 11.88747689358562e0_wp, &
        & 0.2617252059369060e0_wp, -1.142054614272080e0_wp, &
        & 0.2728967965513714e0_wp, -1.9323083032991791e-2_wp, &
        & 57.04567478144305e0_wp, 148.5636591372696e0_wp, &
        & -318.6179220492799e0_wp, 182.3742028409573e0_wp, &
        & -144.0665991614994e0_wp, 62.04257754042266e0_wp, &
        & -5.877331776092552e0_wp, 0.4288520217861720e0_wp, &
        & -0.1036073697125561e0_wp, 120.0475914921467e0_wp, &
        & 263.6164106497357e0_wp, -462.6236113095622e0_wp, &
        & 116.0608660637904e0_wp, -19.91101916036278e0_wp, &
        & 1.193693981250914e0_wp, 0.1422177054289994e0_wp, &
        & -160.3114246747076e0_wp, 89.72253917752164e0_wp, &
        & -9.071109638571208e0_wp, -2.0135523158362437e-2_wp, &
        & 0.1243557501187861e0_wp, -12.96200461060273e0_wp, &
        & -2.347887744857734e0_wp, 0.4369816290614699e0_wp, &
        & 0.1861905441025299e0_wp, -8.284586975531951e0_wp, &
        & -20.42241348042063e0_wp, 958.1242541736171e0_wp, &
        & -1928.583970392355e0_wp, 1119.594819408116e0_wp, &
        & -696.0896180464565e0_wp, 272.2524874095656e0_wp, &
        & -17.44951726252673e0_wp, -1.588843767916973e0_wp, &
        & -8.3560790211881619e-2_wp, 2386.104049771720e0_wp, &
        & 3941.711173912721e0_wp, -5088.409610323098e0_wp, &
        & 2147.437212950090e0_wp, -1209.070332780972e0_wp, &
        & 228.5650959646700e0_wp, -18.75623216539055e0_wp, &
        & 1.331518943744410e0_wp, 1323.582901230102e0_wp, &
        & 1366.507728689508e0_wp, -470.6376874755080e0_wp, &
        & -4.840542534010448e0_wp, 4.418210908054986e0_wp, &
        & 0.2615031877550023e0_wp, 365.2205662123347e0_wp, &
        & 56.60696871480308e0_wp, 20.20287252717788e0_wp, &
        & -8.939307433105093e0_wp, -7.579362430843242e0_wp, &
        & 1.945102121992439e0_wp, 23893.31732811452e0_wp, &
        & 14785.75023538304e0_wp, -5590.061708848788e0_wp, &
        & 2659.366722765638e0_wp, -1103.540790332176e0_wp, &
        & 133.2401706836243e0_wp, -2.779352658653005e0_wp, &
        & -2487.131159960726e0_wp, 223.7986868185370e0_wp, &
        & 940.7101507094950e0_wp, -119.8581541479964e0_wp, &
        & -5.864853168420393e0_wp, -1471.392587738388e0_wp, &
        & 20.46733470843232e0_wp, 44.81137223542368e0_wp, &
        & -37.13875660767655e0_wp, -4114.467297104680e0_wp, &
        & 435.2902050180833e0_wp, 110.1713440583027e0_wp, &
        & -6.086172632709467e0_wp, -58.60516318227123e0_wp, &
        & -30.37217546242515e0_wp, 114.1717244633075e0_wp/)
    real(wp) :: xlin2(89)=(/-2.3691917237142486e-2_wp, &
        & 0.6433593796852592e0_wp, -3.191080817325861e0_wp, &
        & 5.082622923008778e0_wp, -4.303785253384310e0_wp, &
        & 2.677917664616706e0_wp, -1.482618068308200e0_wp, &
        & 0.6784912558836848e0_wp, -0.2049921467505755e0_wp, &
        & 3.1560425722410100e-2_wp, -1.9786106444049946e-3_wp, &
        & -1.510131880107122e0_wp, -28.60321443064080e0_wp, &
        & 53.22691305234915e0_wp, -40.11411939395515e0_wp, &
        & 20.81861842742526e0_wp, -4.447977855272196e0_wp, &
        & 0.2301695113250538e0_wp, -0.1295406658545045e0_wp, &
        & 2.4496342125756746e-2_wp, 23.93815251101609e0_wp, &
        & -15.97857074145742e0_wp, -19.75274167695325e0_wp, &
        & 15.86955351115030e0_wp, 0.4841116429447808e0_wp, &
        & -1.3958200902611101e-2_wp, -6.8392975615942639e-2_wp, &
        & -6.396616748692148e0_wp, -0.5798087307661497e0_wp, &
        & 1.564469079111058e0_wp, -9.8289900209482819e-2_wp, &
        & -4.9186020637759537e-2_wp, 16.99519683721401e0_wp, &
        & -2.466567635268690e0_wp, 1.0845662786683728e-2_wp, &
        & 0.2941408959952037e0_wp, 1.777059184168285e0_wp, &
        & 46.60356704672536e0_wp, -149.4923339330999e0_wp, &
        & 153.0203582840253e0_wp, -129.5745301688355e0_wp, &
        & 83.01390198308063e0_wp, -20.50411666044072e0_wp, &
        & 1.329739978962026e0_wp, 8.1037167237220931e-2_wp, &
        & -5.6360423765908704e-3_wp, -224.7409562531643e0_wp, &
        & -445.7702964503590e0_wp, -182.9799944501726e0_wp, &
        & -133.0432327364785e0_wp, 227.3077678979116e0_wp, &
        & -34.61637061865006e0_wp, -8.9508185399590101e-2_wp, &
        & 7.1955361459616729e-2_wp, 442.7525395172289e0_wp, &
        & 2.458625503871524e0_wp, -202.8428507849074e0_wp, &
        & 32.31023899069420e0_wp, -2.282723397214718e0_wp, &
        & 0.3055028229001248e0_wp, 181.8487130166304e0_wp, &
        & -69.54292969634861e0_wp, -6.889401006201262e0_wp, &
        & 2.218301674443548e0_wp, 32.61450081881122e0_wp, &
        & -1.829579173006778e0_wp, -3053.146122439747e0_wp, &
        & 2039.754739773428e0_wp, -1925.218456128402e0_wp, &
        & 149.3339608613740e0_wp, 28.43501684561451e0_wp, &
        & 0.3254831648262675e0_wp, -0.3970879353555506e0_wp, &
        & 3299.500750582070e0_wp, -594.7502330927563e0_wp, &
        & -65.96404318767495e0_wp, 28.40550550560758e0_wp, &
        & -1.865692389204451e0_wp, 448.8499645057911e0_wp, &
        & -30.33721822994903e0_wp, -3.221780195307217e0_wp, &
        & 3.085642390961864e0_wp, 498.9201490733377e0_wp, &
        & -151.2816443510426e0_wp, -76.13951656899020e0_wp, &
        & 3.665378992703659e0_wp, 72.17259482865131e0_wp, &
        & 6.152241900437081e0_wp, -31.61340871373989e0_wp/)
    real(wp) :: xlin3(71)=(/-0.1344994646896329e0_wp, &
        & -3.2141410415166161e-2_wp, -7.517359288762402e0_wp, &
        & 18.44051637133763e0_wp, -14.31587410374696e0_wp, &
        & -2.834545980827915e0_wp, 12.12504845332679e0_wp, &
        & -8.087779341389851e0_wp, 2.266439914858519e0_wp, &
        & -0.2458975078101396e0_wp, -24.78954125497926e0_wp, &
        & -130.9696482818210e0_wp, 208.5056810553900e0_wp, &
        & -102.7798259674824e0_wp, 94.50458094288865e0_wp, &
        & -26.66242091873472e0_wp, -3.796095988884443e0_wp, &
        & -0.2135440885386056e0_wp, -101.9329951911512e0_wp, &
        & -237.0984558426526e0_wp, 340.3675666528260e0_wp, &
        & -5.600055272704649e0_wp, -6.626620457234992e0_wp, &
        & 2.250749424450074e0_wp, 72.61834696241013e0_wp, &
        & 28.58895062725391e0_wp, -25.96494974227168e0_wp, &
        & 1.775974157654598e0_wp, -10.06825158236352e0_wp, &
        & 7.123101687824712e0_wp, 2.773473562095701e0_wp, &
        & 34.25830296691608e0_wp, -613.3162919232790e0_wp, &
        & 1003.318114503314e0_wp, -293.8447014270225e0_wp, &
        & 238.0037203719179e0_wp, -45.54821241600535e0_wp, &
        & -38.15460316807260e0_wp, 2.152193113454373e0_wp, &
        & -1265.105998633106e0_wp, -4182.774652447684e0_wp, &
        & 4323.977859927758e0_wp, -1511.567227580535e0_wp, &
        & 1013.128818248927e0_wp, -145.0171911214698e0_wp, &
        & 17.32391332654659e0_wp, -2672.569843994831e0_wp, &
        & -494.6368772637625e0_wp, -176.6910389060321e0_wp, &
        & 333.1709756767322e0_wp, -44.43847568245701e0_wp, &
        & -715.3808858788194e0_wp, -86.17184561227721e0_wp, &
        & -19.32116072356918e0_wp, -16.23431382053576e0_wp, &
        & -17732.22229184206e0_wp, -21836.37425726624e0_wp, &
        & 6027.463948377972e0_wp, -3413.560409537637e0_wp, &
        & 1462.120516808686e0_wp, -132.4258810442149e0_wp, &
        & -5461.158978972337e0_wp, 1036.230524299047e0_wp, &
        & -1687.326937421790e0_wp, 136.8271126588730e0_wp, &
        & 1563.706329850235e0_wp, 55.53727444282163e0_wp, &
        & 6546.413397438990e0_wp, 273.3238154461579e0_wp, &
        & -123.8048659677774e0_wp, -69.99455204122073e0_wp/)
    real(wp) :: xlin4(71)=(/1.2472601237000719e-2_wp, &
        & -0.2563237394846267e0_wp, 1.141626638949766e0_wp, &
        & -1.268789070877211e0_wp, 0.2669120728039442e0_wp, &
        & 0.4795337070923109e0_wp, -0.4870693693992560e0_wp, &
        & 0.2180413266051965e0_wp, -5.4940116209431157e-2_wp, &
        & 8.1898224212963791e-3_wp, -0.9232881113457609e0_wp, &
        & 22.37373927583910e0_wp, -31.43191815577488e0_wp, &
        & 15.83511958871624e0_wp, -4.774566352751313e0_wp, &
        & -1.809376059889746e0_wp, 0.8755187623690739e0_wp, &
        & 3.5850954944800173e-2_wp, -7.088532812068095e0_wp, &
        & 5.757713051107022e0_wp, 10.86426353414864e0_wp, &
        & -5.823082328751028e0_wp, -3.528581474490192e0_wp, &
        & -0.1229368602139166e0_wp, -10.79753576660658e0_wp, &
        & 24.25244965677145e0_wp, -8.953867875122164e0_wp, &
        & 0.4615248048240077e0_wp, -23.16730300552956e0_wp, &
        & 0.3228037267248082e0_wp, -0.4960917683113312e0_wp, &
        & -29.27293769041133e0_wp, 77.13556531104766e0_wp, &
        & -65.28389682361755e0_wp, 48.11222198276875e0_wp, &
        & -29.05873944480265e0_wp, -0.1063025548510765e0_wp, &
        & 2.448839826203129e0_wp, -0.2802338471724751e0_wp, &
        & 51.92481774609280e0_wp, 530.2135642610900e0_wp, &
        & -9.028255840313982e0_wp, 346.5780338983027e0_wp, &
        & -240.1913305251466e0_wp, 8.295317216645913e0_wp, &
        & -0.2022074738812011e0_wp, -181.7378257010903e0_wp, &
        & 56.03276885734908e0_wp, 155.3619996347792e0_wp, &
        & -45.12505868315263e0_wp, 10.96008703215353e0_wp, &
        & -242.7878951933611e0_wp, 101.8270659259155e0_wp, &
        & -0.4864797895058910e0_wp, -7.874126046335711e0_wp, &
        & 2460.968493617616e0_wp, -65.55982339327284e0_wp, &
        & 1112.649677115970e0_wp, 300.2383741022979e0_wp, &
        & -150.4025249750507e0_wp, 23.70714634112154e0_wp, &
        & -4536.070861321965e0_wp, 451.5619825013044e0_wp, &
        & 167.3937904161444e0_wp, -51.79228062736297e0_wp, &
        & -288.7604569805079e0_wp, 10.58119538655392e0_wp, &
        & -2266.622734963641e0_wp, 164.3048434953309e0_wp, &
        & 102.2707553778119e0_wp, -71.61259882125256e0_wp/)
end module

!> @brief Evaluate a batch of H3 A2 geometries.
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
        call pota2(r1x, r2x, r3x, evp(ivp))
    enddo
end subroutine pot

!> @brief Evaluate one H3 A2 geometry.
subroutine pota2( r1x, r2x, r3x, eee )
    use a2data,                  only: ald, beta, betad, epscem, epslon, eshift, maxi2a, maxi2b, rnought, xlin1, xlin2, xlin3, xlin4
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
end subroutine pota2

!> @brief Evaluate trip for the H3_GEN_A2_2002 potential.
subroutine trip( r, v )
    !     H2 triplet curve for FCI/aug-cc-pVDZ
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: j
    real(wp) :: r, v, damp, xd, prefac
    real(wp) :: c6=-6.499027e0_wp
    real(wp) :: c8=-124.3991e0_wp
    real(wp) :: c10=-3285.828e0_wp
    real(wp) :: beta=0.2e0_wp
    real(wp) :: re=1.401e0_wp
    real(wp) :: alpha=3.9139661662872520e0_wp
    real(wp) :: xlp(17)=(/2.2349220271817220e-01_wp, &
        & 6.0356997147268380e-01_wp, 8.9383829921944970e-01_wp, &
        & 9.0317491929985030e-01_wp, 7.2775729411232600e-01_wp, &
        & 4.9697572693224850e-01_wp, 2.5866961299917790e-01_wp, &
        & 6.8824590522579050e-02_wp, 4.8273710815396330e-03_wp, &
        & 2.8416581979833760e-03_wp, -5.7032138025732140e-03_wp, &
        & 2.3004916075318040e-03_wp, 3.1740661346757730e-03_wp, &
        & -2.3028140070272670e-03_wp, 6.4778830165694920e-04_wp, &
        & -8.5334197479958010e-05_wp, 4.5529457450819180e-06_wp/)
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

!> @brief Evaluate singlet for the H3_GEN_A2_2002 potential.
subroutine singlet( r, v )
    !     H2 singlet potential for FCI/aug-cc-pVDZ
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: j
    real(wp) :: r, v, damp, xd, prefac
    real(wp) :: c6=-6.499027e0_wp
    real(wp) :: c8=-124.3991e0_wp
    real(wp) :: c10=-3285.828e0_wp
    real(wp) :: beta=0.2e0_wp
    real(wp) :: re=1.401e0_wp
    real(wp) :: alpha=4.0578706387229880e0_wp
    real(wp) :: xlp(17)= (/-1.6244445106378250e-1_wp, &
        & -6.6081156795821980e-01_wp, -1.1236729331029230e+00_wp, &
        & -1.1564282634392240e+00_wp, -7.9003598933823580e-01_wp, &
        & -3.7515947232144770e-01_wp, -1.3997385158833170e-01_wp, &
        & -6.5955219641305790e-02_wp, -4.0788777170244650e-02_wp, &
        & -1.8988748368877420e-02_wp, -2.2257283819080430e-03_wp, &
        & 2.5357629098565200e-03_wp, 1.5612692551680850e-05_wp, &
        & -6.9340042488041930e-04_wp, 2.3175123448115890e-04_wp, &
        & -3.3226980442299310e-05_wp, 1.6041565480651940e-06_wp/)
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
