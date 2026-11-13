!   System:                        H3
!   Functional form:
!   Common name:                   A4
!   Number of derivatives:         0
!   Number of bodies:              3
!   Number of electronic surfaces: 1
!   Interface:                     3-1V
!   Data file:
!   References:      SL Mielke, BC Garrett, and KA Peterson, J. Chem. Ph
!   Notes:    3 H-H distances passed in bohr.
!             Energy returned in Hatrees relative to the classical minim

!> @brief Implement H3_GEN_A4_2002, a potential energy surface for H3.
module chempotpy_h3_gen_a4_2002_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_A4_2002 potential energy surface for H3.
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
end module chempotpy_h3_gen_a4_2002_pes

!> @brief Expose H3_GEN_A4_2002 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_a4_2002_pes, only: evaluate
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

module a4data
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
    real(wp), parameter :: eshift  = 0.1739709107962980e0_wp
    real(wp), parameter :: beta    = 0.9201192312615871e0_wp
    real(wp), parameter :: rnought = 2.e0_wp
    real(wp)            :: xlin1(89)=(/-0.1747684349138167e0_wp, &
        & 1.277318541919902e0_wp, 4.865120178227985e0_wp, &
        & -8.469658615576670e0_wp, -0.5743180823752541e0_wp, &
        & 6.240393306718566e0_wp, -2.708137131341480e0_wp, &
        & -0.5692746751976894e0_wp, 0.6546217258826842e0_wp, &
        & -0.1509110106019215e0_wp, 1.0897871662071377e-2_wp, &
        & 40.08907725517741e0_wp, 78.14109140356737e0_wp, &
        & -114.9895503001828e0_wp, 14.58341862464093e0_wp, &
        & 32.62546038173688e0_wp, -17.51073267443663e0_wp, &
        & 4.374012288933156e0_wp, -0.3395820839216175e0_wp, &
        & -7.0640001462913965e-3_wp, -151.5692355055670e0_wp, &
        & 21.01354666590584e0_wp, -55.11286823230184e0_wp, &
        & 2.770716740885964e0_wp, 5.009995499090712e0_wp, &
        & -0.6464346348303938e0_wp, 8.3910658414713523e-3_wp, &
        & 49.25955250094846e0_wp, -9.603115825941249e0_wp, &
        & 4.8272796970700860e-3_wp, 6.5210196404216089e-2_wp, &
        & 1.2052781423687001e-2_wp, 1.153704488560043e0_wp, &
        & 0.6779128659000427e0_wp, -8.7223899834716512e-2_wp, &
        & -0.1303377419808928e0_wp, -10.90700005830876e0_wp, &
        & 145.7535863536604e0_wp, 222.5532957606632e0_wp, &
        & -153.3725179930681e0_wp, -65.61649004631528e0_wp, &
        & 98.77749621876600e0_wp, -37.38968164171358e0_wp, &
        & 6.448018323092472e0_wp, 0.1406206598020832e0_wp, &
        & -9.1763752206667951e-2_wp, 577.9387530303504e0_wp, &
        & 321.2258522767713e0_wp, -1956.057673174400e0_wp, &
        & -75.51415928657248e0_wp, 207.0415614864125e0_wp, &
        & -44.79804897124119e0_wp, 1.427806116568303e0_wp, &
        & 0.3525406620212973e0_wp, -1255.200747980566e0_wp, &
        & 196.2632339741470e0_wp, -138.6365615558797e0_wp, &
        & -10.64943395699545e0_wp, 6.425094067862800e0_wp, &
        & -0.3444399405796634e0_wp, 90.36916178542219e0_wp, &
        & 16.37161644951711e0_wp, -7.480040057276788e0_wp, &
        & 0.4011296752128877e0_wp, 1.298019306755686e0_wp, &
        & 0.2891734971459664e0_wp, 1285.602684915883e0_wp, &
        & 4565.798100354337e0_wp, -3155.673834080026e0_wp, &
        & -44.92705598116694e0_wp, 168.9260873699247e0_wp, &
        & -7.844938785348392e0_wp, -1.666140934947715e0_wp, &
        & 2779.863836531601e0_wp, 537.1244197181024e0_wp, &
        & -170.2582980443451e0_wp, -5.290922323422062e0_wp, &
        & 1.770830652534150e0_wp, 52.20593427757114e0_wp, &
        & 9.347628365683427e0_wp, -3.092797446921548e0_wp, &
        & 2.379156980285720e0_wp, -3334.415420363777e0_wp, &
        & 161.8383070634118e0_wp, 19.96531657285394e0_wp, &
        & 1.327410294913859e0_wp, -21.45733100652736e0_wp, &
        & -0.1662058917002993e0_wp, -0.2962531422178993e0_wp/)
    real(wp) :: xlin2(89)=(/3.6653811266866398e-2_wp, &
        & -6.0901908406171948e-2_wp, -0.4296166129329135e0_wp, &
        & -0.5999553551949833e0_wp, 2.516208126297607e0_wp, &
        & -2.048432510828575e0_wp, 0.3663315280773198e0_wp, &
        & 0.2774060109660687e0_wp, -0.1575944525599335e0_wp, &
        & 2.9695991705353031e-2_wp, -1.9203688037688026e-3_wp, &
        & -0.5834482485892581e0_wp, -10.26004544230231e0_wp, &
        & 7.260618335575135e0_wp, 9.670699153937617e0_wp, &
        & -12.63290323745150e0_wp, 4.813462603326546e0_wp, &
        & -0.4949553387010042e0_wp, -6.7585812514542612e-2_wp, &
        & 1.0929473861305191e-2_wp, 38.39652796798982e0_wp, &
        & -8.046863706246480e0_wp, -9.322531512667055e0_wp, &
        & 8.319538557744806e0_wp, -0.5844887973497459e0_wp, &
        & -0.2072100517931691e0_wp, 1.7662601024712416e-2_wp, &
        & 10.84542783635474e0_wp, 2.719206724534208e0_wp, &
        & -2.661681915907651e0_wp, 5.4511820020259769e-2_wp, &
        & 4.0917094713795558e-2_wp, 2.310535942300969e0_wp, &
        & -0.3313532425660889e0_wp, -1.9282783753774072e-2_wp, &
        & 0.1076070240762595e0_wp, -1.535785256390712e0_wp, &
        & 21.55814329603826e0_wp, -20.67412669364389e0_wp, &
        & -42.22253006381749e0_wp, 63.24620365892324e0_wp, &
        & -29.72684500701624e0_wp, 5.845208967229377e0_wp, &
        & 0.1479499544852192e0_wp, -0.1927854989069523e0_wp, &
        & 1.6801714288723839e-2_wp, 137.1873413615389e0_wp, &
        & -257.7551343332649e0_wp, 135.0209376494197e0_wp, &
        & 13.64539095240564e0_wp, 12.09164291321869e0_wp, &
        & -5.169736649305458e0_wp, 0.8698478087116743e0_wp, &
        & -8.0291725093975599e-2_wp, -347.8994751159709e0_wp, &
        & -29.90289317716435e0_wp, 2.681813980780868e0_wp, &
        & 8.722110659596989e0_wp, -0.2420626678309532e0_wp, &
        & -0.1274172885241889e0_wp, 55.06868860080869e0_wp, &
        & -12.74535530370019e0_wp, -0.2131001164206556e0_wp, &
        & 0.1837891447312188e0_wp, 3.041924830651179e0_wp, &
        & -0.2126407877934366e0_wp, -742.3965481348824e0_wp, &
        & -784.5779176647559e0_wp, -204.0838570881213e0_wp, &
        & 27.41420672246095e0_wp, 27.81184144136420e0_wp, &
        & -5.107775576751197e0_wp, 0.2426879842278437e0_wp, &
        & 673.8570609208289e0_wp, -110.1522754394637e0_wp, &
        & -30.37003163984198e0_wp, 2.612700034489588e0_wp, &
        & 0.4985395581182200e0_wp, 67.07263299708157e0_wp, &
        & -2.044069233731526e0_wp, -1.022956710686919e0_wp, &
        & 0.9605395343509358e0_wp, 302.0563415645296e0_wp, &
        & -34.61819182491933e0_wp, -2.686757365342857e0_wp, &
        & -0.7220696373666605e0_wp, 6.015885938467894e0_wp, &
        & 1.207183412067369e0_wp, -4.708046093073749e0_wp/)
    real(wp) :: xlin3(71)=(/6.2633244088892953e-2_wp, &
        & -0.3594764662369648e0_wp, -3.600210064900942e0_wp, &
        & 6.077024302044581e0_wp, -2.103247648143528e0_wp, &
        & -7.5336710734737355e-2_wp, -1.214802246636612e0_wp, &
        & 1.424478564483824e0_wp, -0.5152381133998951e0_wp, &
        & 5.8690423728471254e-2_wp, -21.07550548519476e0_wp, &
        & -68.04052224304678e0_wp, 57.93242378809885e0_wp, &
        & 27.37591131596061e0_wp, -42.01788015786860e0_wp, &
        & 17.01097780464082e0_wp, -3.810546557911935e0_wp, &
        & 0.1341788526198947e0_wp, 39.21740752426907e0_wp, &
        & 69.74078714802111e0_wp, 19.18393280026359e0_wp, &
        & 8.117596367547538e0_wp, -3.122024525347344e0_wp, &
        & -0.2382844161241255e0_wp, -89.25083728779680e0_wp, &
        & 25.77704137475598e0_wp, -1.969043985146592e0_wp, &
        & -0.4855106421075169e0_wp, 0.7855919428604119e0_wp, &
        & -0.8441424883696753e0_wp, 9.448681574551424e0_wp, &
        & -79.69875321240659e0_wp, -193.2070036563253e0_wp, &
        & 42.92224561943964e0_wp, 111.7341031859155e0_wp, &
        & -104.4278320210573e0_wp, 37.98166156896379e0_wp, &
        & -7.211132085860729e0_wp, -4.6723492922660093e-3_wp, &
        & -538.7716000645521e0_wp, -545.3793327560588e0_wp, &
        & 1369.806489379943e0_wp, 584.9126962811058e0_wp, &
        & -235.2608524242808e0_wp, 44.32276578556720e0_wp, &
        & -1.162019393434994e0_wp, 1190.467988321617e0_wp, &
        & 620.3261156354067e0_wp, 108.9394942519943e0_wp, &
        & 6.060056785130230e0_wp, 2.581991625425053e0_wp, &
        & -218.0491226601461e0_wp, 13.36602093278375e0_wp, &
        & 0.7718484142277729e0_wp, 0.8430743860674987e0_wp, &
        & -957.4048659566447e0_wp, -4033.862537342876e0_wp, &
        & 3191.552283364393e0_wp, 574.5261882010860e0_wp, &
        & -220.9679992384400e0_wp, 9.589933194356654e0_wp, &
        & -4986.101521958573e0_wp, -483.1000657856682e0_wp, &
        & 16.05429422497049e0_wp, -15.62404181993110e0_wp, &
        & -90.73046075338512e0_wp, 3.835596218545594e0_wp, &
        & 2040.882666165481e0_wp, 168.7823751649685e0_wp, &
        & 84.39701838464045e0_wp, -72.00299079742861e0_wp/)
    real(wp) :: xlin4(71)=(/-1.3789322831534202e-2_wp, &
        & -5.3276062951444891e-3_wp, 0.3280620147342902e0_wp, &
        & 2.7087208142318464e-3_wp, -0.4121991913062390e0_wp, &
        & -0.2194705118222190e0_wp, 0.7128681500449559e0_wp, &
        & -0.4529979958049952e0_wp, 0.1189948966982933e0_wp, &
        & -1.1240573716569334e-2_wp, 0.2786903571008327e0_wp, &
        & 6.967391109018378e0_wp, -0.3729050190579263e0_wp, &
        & -12.24220706099335e0_wp, 10.19906396811849e0_wp, &
        & -2.638016620492531e0_wp, 9.1394461539987537e-3_wp, &
        & 5.5889568628281538e-2_wp, -16.32524207951339e0_wp, &
        & -10.27043474888980e0_wp, 8.195640957287562e0_wp, &
        & -1.026303873216019e0_wp, -2.533027615178429e0_wp, &
        & 0.3931306992906919e0_wp, 5.700236993312759e0_wp, &
        & -5.964165300509528e0_wp, -4.2236603189111335e-2_wp, &
        & 0.5679808795239981e0_wp, -2.155880855745120e0_wp, &
        & 0.3863665921458664e0_wp, 0.4200525718978940e0_wp, &
        & -8.018391528470787e0_wp, -3.363559188467546e0_wp, &
        & 53.69613583466297e0_wp, -52.75603349194013e0_wp, &
        & 19.06668600609882e0_wp, -2.550942982823412e0_wp, &
        & -0.3464476386899833e0_wp, 7.6947106143076857e-2_wp, &
        & -122.7689330815471e0_wp, 158.8043086748957e0_wp, &
        & -45.97927653390762e0_wp, -56.47466715674170e0_wp, &
        & 0.6113508220195394e0_wp, -2.641958330496622e0_wp, &
        & 0.1064836300639488e0_wp, 343.3356986270791e0_wp, &
        & 102.5423266427451e0_wp, -0.4240044167349897e0_wp, &
        & -11.09858306988299e0_wp, 0.2559185246181059e0_wp, &
        & -81.76030598060535e0_wp, 11.92782620599065e0_wp, &
        & -0.6869535394737264e0_wp, -2.253291923744354e0_wp, &
        & 245.0489145059600e0_wp, 1182.806070145991e0_wp, &
        & 152.7550042593918e0_wp, 45.24212872272277e0_wp, &
        & -30.12697910280025e0_wp, 2.918896114976833e0_wp, &
        & -184.2429751946113e0_wp, 83.71871876387934e0_wp, &
        & 19.99821094524841e0_wp, -1.471571531480851e0_wp, &
        & -47.00968833983665e0_wp, 2.161626248435158e0_wp, &
        & -971.4885779654726e0_wp, 58.68272991198457e0_wp, &
        & 7.609238909723311e0_wp, -8.492629420669424e0_wp/)
end module

!> @brief Evaluate a batch of H3 A4 geometries.
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
        call pota4(r1x, r2x, r3x, evp(ivp))
    enddo
end subroutine pot

!> @brief Evaluate one H3 A4 geometry.
subroutine pota4( r1x, r2x, r3x, eee )
    use a4data,                  only: ald, beta, betad, epscem, epslon, eshift, maxi2a, maxi2b, rnought, xlin1, xlin2, xlin3, xlin4
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
end subroutine pota4

!> @brief Evaluate trip for the H3_GEN_A4_2002 potential.
subroutine trip( r, v )
    !     H2 triplet curve for FCI/aug-cc-pVQZ
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: j
    real(wp) :: r, v, damp, xd, prefac
    real(wp) :: c6=-6.499027e0_wp
    real(wp) :: c8=-124.3991e0_wp
    real(wp) :: c10=-3285.828e0_wp
    real(wp) :: beta=0.2e0_wp
    real(wp) :: re=1.401e0_wp
    real(wp) :: alpha= 3.8542531266774340e0_wp
    real(wp) :: xlp(17)=(/2.1913359720738420e-1_wp, &
        & 5.8422835449484010e-1_wp, 8.5061017374705010e-1_wp, &
        & 8.4614320758188910e-1_wp, 6.7185038364340950e-1_wp, &
        & 4.4194539374617390e-1_wp, 2.2339961981172340e-1_wp, &
        & 6.5347598849164070e-2_wp, -3.1377420901498310e-3_wp, &
        & -8.3737736675708600e-3_wp, 3.4989034143871880e-3_wp, &
        & 2.6754734558345340e-3_wp, -7.7370541086400900e-4_wp, &
        & -1.6880774643325540e-4_wp, 1.1493550202349050e-4_wp, &
        & -1.9309505011102930e-5_wp, 1.1841183250969540e-6_wp/)
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

!> @brief Evaluate singlet for the H3_GEN_A4_2002 potential.
subroutine singlet( r, v )
    !     H2 singlet potential for FCI/aug-cc-pVQZ
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: j
    real(wp) :: r, v, damp, xd, prefac
    real(wp) :: c6=-6.499027e0_wp
    real(wp) :: c8=-124.3991e0_wp
    real(wp) :: c10=-3285.828e0_wp
    real(wp) :: beta=0.2e0_wp
    real(wp) :: re=1.401e0_wp
    real(wp) :: alpha= 4.0590126171616700e0_wp
    real(wp) :: xlp(17)= (/-1.7046099537113510e-1_wp, &
        & -6.7922917067867770e-1_wp, -1.1517142440791580e0_wp, &
        & -1.1910821267977240e0_wp, -8.2784442225164410e-1_wp, &
        & -4.0807832245803330e-1_wp, -1.5672842672450390e-1_wp, &
        & -6.8869155703744850e-2_wp, -4.3662333239840700e-2_wp, &
        & -2.3870880222600530e-2_wp, -3.4451065483225290e-3_wp, &
        & 3.2887970319820380e-3_wp, 1.9772166585174660e-4_wp, &
        & -7.4858600680843820e-4_wp, 1.8385885697114190e-4_wp, &
        & -1.6249086298603110e-5_wp, 3.6338862344567620e-10_wp/)
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
