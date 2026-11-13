!   System:                        H3
!   Functional form:
!   Common name:                   CCI
!   Number of derivatives:         0
!   Number of bodies:              3
!   Number of electronic surfaces: 1
!   Interface:                     3-1V
!   Data file:
!   References:      SL Mielke, BC Garrett, and KA Peterson, J. Chem. Ph
!   Notes:    3 H-H distances passed in bohr.
!             Energy returned in Hatrees relative to the classical minim

!> @brief Implement H3_GEN_CCI_2002, a potential energy surface for H3.
module chempotpy_h3_gen_cci_2002_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_CCI_2002 potential energy surface for H3.
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
end module chempotpy_h3_gen_cci_2002_pes

!> @brief Expose H3_GEN_CCI_2002 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_cci_2002_pes, only: evaluate
    use chempotpy_kinds,               only: wp
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

module ccidata
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
    real(wp), parameter :: eshift  = 0.1744759989649372e0_wp
    real(wp), parameter :: beta    = 0.9246250191375499e0_wp
    real(wp), parameter :: rnought = 2.e0_wp
    real(wp)            :: xlin1(89)=(/-0.1682626669809615e0_wp, &
        & 1.697952225005269e0_wp, 3.959367441888326e0_wp, &
        & -7.159854944676235e0_wp, -2.462488887191634e0_wp, &
        & 7.998434264390882e0_wp, -3.665712084901418e0_wp, &
        & -0.2155463656483861e0_wp, 0.5620418518154215e0_wp, &
        & -0.1369812453279634e0_wp, 1.0063614579308723e-2_wp, &
        & 41.07358514471853e0_wp, 67.02217912830267e0_wp, &
        & -100.0014758720550e0_wp, 2.481875721304641e0_wp, &
        & 32.46244740624601e0_wp, -15.32460432587821e0_wp, &
        & 3.814631583515083e0_wp, -0.3155985331636701e0_wp, &
        & -4.2198516063816239e-3_wp, -157.7203098415152e0_wp, &
        & 17.11141015275416e0_wp, -61.48638359701442e0_wp, &
        & 4.244841908917134e0_wp, 5.300393723900312e0_wp, &
        & -0.6170558358973661e0_wp, -6.2143674437616274e-3_wp, &
        & 5.785766830955630e0_wp, -0.1133784523542150e0_wp, &
        & -0.2539145873400511e0_wp, -4.3406559735048131e-2_wp, &
        & 2.4548851428702549e-2_wp, -1.092841685240530e0_wp, &
        & 0.8131178021738225e0_wp, -7.7933784638615652e-2_wp, &
        & -0.1484455564689179e0_wp, -7.351078313940695e0_wp, &
        & 128.7480567271992e0_wp, 199.6058760392980e0_wp, &
        & -50.18444648711906e0_wp, -178.2184990028168e0_wp, &
        & 129.6832433574917e0_wp, -37.57036725437851e0_wp, &
        & 5.747838042162251e0_wp, 0.1434298979682676e0_wp, &
        & -8.3110343245055174e-2_wp, 428.7635440340811e0_wp, &
        & 495.9564113906064e0_wp, -1742.123837168904e0_wp, &
        & -259.0564024322437e0_wp, 250.2352551214438e0_wp, &
        & -44.15232947384926e0_wp, 0.9963955068997190e0_wp, &
        & 0.3136096547683901e0_wp, -1032.921631265462e0_wp, &
        & 199.9708129167335e0_wp, -135.4683186684037e0_wp, &
        & -14.50534672808297e0_wp, 5.650612567917722e0_wp, &
        & -0.1473535040913008e0_wp, 1.351710421967908e0_wp, &
        & 29.60972395839821e0_wp, -7.224575522867537e0_wp, &
        & 0.3582134126835411e0_wp, 0.2589923362572861e0_wp, &
        & 0.1819402887127863e0_wp, 1312.633253770123e0_wp, &
        & 5189.281073844522e0_wp, -2727.020597875299e0_wp, &
        & -122.2193440638810e0_wp, 152.4859537357477e0_wp, &
        & -5.562127541501408e0_wp, -1.549418235204051e0_wp, &
        & 2647.710764337011e0_wp, 502.1991025406292e0_wp, &
        & -152.3069488960497e0_wp, -5.740041181528469e0_wp, &
        & 1.362209926084833e0_wp, 17.96504193139269e0_wp, &
        & 10.13930464010770e0_wp, -2.287463377365149e0_wp, &
        & 1.932492687974671e0_wp, -3337.181754475553e0_wp, &
        & 182.9671187483376e0_wp, 20.29606065582956e0_wp, &
        & 1.192500349837268e0_wp, -22.69079850990068e0_wp, &
        & -0.5967747962499963e0_wp, 1.091111370927579e0_wp/)
    real(wp) :: xlin2(89)=(/3.2289532152874986e-2_wp, &
        & -8.6686506736620789e-2_wp, -0.3266024449171150e0_wp, &
        & -0.7972798579742405e0_wp, 2.812854025206053e0_wp, &
        & -2.374871328104482e0_wp, 0.6089164699935086e0_wp, &
        & 0.1635658322754857e0_wp, -0.1263080841790984e0_wp, &
        & 2.5221027611065153e-2_wp, -1.6677760633386995e-3_wp, &
        & -0.3372069966080046e0_wp, -8.377537389356892e0_wp, &
        & 3.694668699246826e0_wp, 11.36868357221171e0_wp, &
        & -11.88889428297186e0_wp, 4.270226804928180e0_wp, &
        & -0.4526863369836274e0_wp, -5.2518975083660345e-2_wp, &
        & 9.0223742463199624e-3_wp, 38.03127084604298e0_wp, &
        & -12.82065963469141e0_wp, -5.663705994784566e0_wp, &
        & 7.300080908068341e0_wp, -0.5157232521089076e0_wp, &
        & -0.2200954472805637e0_wp, 2.0686337667851093e-2_wp, &
        & 14.19614073502447e0_wp, 1.726255538125508e0_wp, &
        & -2.526391344985463e0_wp, 6.0598667264186397e-2_wp, &
        & 3.8560006141770319e-2_wp, 2.306725001955670e0_wp, &
        & -0.3197660708712612e0_wp, -1.9982819587147187e-2_wp, &
        & 0.1036787525350017e0_wp, -2.030231243513758e0_wp, &
        & 25.16474735933574e0_wp, -14.13898387128746e0_wp, &
        & -62.01219930765056e0_wp, 75.56471352787402e0_wp, &
        & -30.25154027581496e0_wp, 4.257663341071418e0_wp, &
        & 0.6142213624561507e0_wp, -0.2458345005137628e0_wp, &
        & 1.9135101060142316e-2_wp, 168.2482649287792e0_wp, &
        & -270.3587962365308e0_wp, 72.46467150857366e0_wp, &
        & 32.12278113625301e0_wp, 11.89435258275814e0_wp, &
        & -5.521268211450904e0_wp, 0.8934653854937415e0_wp, &
        & -7.6805320500971247e-2_wp, -422.2672942154591e0_wp, &
        & -35.92489169679082e0_wp, 14.74665733902546e0_wp, &
        & 6.020706922914545e0_wp, 2.4538599994799737e-3_wp, &
        & -0.1370624616096560e0_wp, 44.47654636843973e0_wp, &
        & -11.18166310447657e0_wp, 3.6110526441458540e-2_wp, &
        & 0.1539535582849851e0_wp, 2.065981296572923e0_wp, &
        & -0.1624439939014451e0_wp, -658.3723053145102e0_wp, &
        & -883.4947620258348e0_wp, -263.8479489759550e0_wp, &
        & 44.68359598895833e0_wp, 29.88970883430300e0_wp, &
        & -5.663152586850906e0_wp, 0.2413124485342995e0_wp, &
        & 656.9813540494517e0_wp, -91.56710742448996e0_wp, &
        & -30.38175110815412e0_wp, 2.066013906655882e0_wp, &
        & 0.5314632792073760e0_wp, 54.99422880821747e0_wp, &
        & -0.7887249891070589e0_wp, -0.9164068632403487e0_wp, &
        & 0.6576346543435602e0_wp, 292.2665706493461e0_wp, &
        & -32.04152511166628e0_wp, -2.104943871500433e0_wp, &
        & -0.7251129002870176e0_wp, 4.709122797258667e0_wp, &
        & 1.119971732384046e0_wp, -4.066519840824071e0_wp/)
    real(wp) :: xlin3(71)=(/6.3998636332917941e-2_wp, &
        & -0.5826083150212633e0_wp, -3.294968443823759e0_wp, &
        & 5.751817019955082e0_wp, -1.486176533347650e0_wp, &
        & -0.6845411439568367e0_wp, -0.9083627566591157e0_wp, &
        & 1.305678536333854e0_wp, -0.4821780067143284e0_wp, &
        & 5.5182580699196751e-2_wp, -22.57213192201081e0_wp, &
        & -61.37436407473290e0_wp, 52.64511100643983e0_wp, &
        & 29.99521825023724e0_wp, -38.31341936937510e0_wp, &
        & 14.86681566870566e0_wp, -3.460352146284503e0_wp, &
        & 0.1343449644795501e0_wp, 53.21459427503272e0_wp, &
        & 64.60699757442232e0_wp, 30.89466488573955e0_wp, &
        & 8.412732049211240e0_wp, -4.274560634850636e0_wp, &
        & -0.1788031871862021e0_wp, -46.47154313325355e0_wp, &
        & 24.88995914944731e0_wp, -2.757433396373881e0_wp, &
        & -0.3635530438220282e0_wp, 1.049428472748707e0_wp, &
        & -0.9909460260393127e0_wp, 7.623381850055868e0_wp, &
        & -75.76387696233419e0_wp, -162.8879929453359e0_wp, &
        & -29.98317728892977e0_wp, 168.7972614497797e0_wp, &
        & -102.0864921135087e0_wp, 31.38363768116054e0_wp, &
        & -6.084136544334489e0_wp, -2.9400115678577254e-2_wp, &
        & -417.4808742201730e0_wp, -601.3246381855070e0_wp, &
        & 1111.017825017578e0_wp, 749.0346863822200e0_wp, &
        & -254.9798923160399e0_wp, 37.02896288679148e0_wp, &
        & -0.7010161055223685e0_wp, 928.6100310886904e0_wp, &
        & 513.3554862860957e0_wp, 138.4248255502364e0_wp, &
        & -0.9347938923678116e0_wp, 3.867048618649862e0_wp, &
        & -109.6794304656507e0_wp, 10.30110710850976e0_wp, &
        & 0.5118243051520034e0_wp, 0.2069543381912151e0_wp, &
        & -796.4010440692600e0_wp, -4693.456970734772e0_wp, &
        & 2597.246850944603e0_wp, 630.4845628327919e0_wp, &
        & -221.0729177223342e0_wp, 9.839396652075672e0_wp, &
        & -5227.790162633692e0_wp, -529.4879950032479e0_wp, &
        & 30.53701394166368e0_wp, -14.22864486129778e0_wp, &
        & -35.14731691232372e0_wp, 2.140914121987492e0_wp, &
        & 2250.367423678786e0_wp, 138.5994085512685e0_wp, &
        & 70.30416282155014e0_wp, -62.74844809925895e0_wp/)
    real(wp) :: xlin4(71)=(/-1.2612782515309912e-2_wp, &
        & 1.5724446750473936e-2_wp, 0.2653887372495336e0_wp, &
        & 0.1216019078797281e0_wp, -0.6102844930670878e0_wp, &
        & 4.6420446436398755e-3_wp, 0.5563651571664677e0_wp, &
        & -0.3890198684837317e0_wp, 0.1052569138798340e0_wp, &
        & -1.0080347996402828e-2_wp, 0.2409741384095758e0_wp, &
        & 5.745106082351191e0_wp, 1.047295525278513e0_wp, &
        & -12.03507114814484e0_wp, 9.035186991905373e0_wp, &
        & -2.267198360038944e0_wp, 2.0836099694222385e-2_wp, &
        & 4.7445938282087222e-2_wp, -18.26537244396756e0_wp, &
        & -4.809068965727743e0_wp, 5.658647992475390e0_wp, &
        & -1.284093508478602e0_wp, -2.164733987255563e0_wp, &
        & 0.3630073042947094e0_wp, 3.760490137844927e0_wp, &
        & -5.927286374087972e0_wp, 0.1278839922444694e0_wp, &
        & 0.5381095462111573e0_wp, -2.179370013756247e0_wp, &
        & 0.4127566904703107e0_wp, 0.7168269469196163e0_wp, &
        & -9.502943939736065e0_wp, -9.597638311349975e0_wp, &
        & 65.60860814240448e0_wp, -56.04003254284753e0_wp, &
        & 15.38423202570607e0_wp, -0.3167914651468724e0_wp, &
        & -0.7643266917502949e0_wp, 0.1110705439245336e0_wp, &
        & -146.5638710976931e0_wp, 151.8079548350792e0_wp, &
        & 15.86486188082157e0_wp, -65.56160562330304e0_wp, &
        & -1.126379669058077e0_wp, -1.759058937005427e0_wp, &
        & -3.1582037336166588e-4_wp, 410.9371301349821e0_wp, &
        & 142.2842073695382e0_wp, -15.14809406881419e0_wp, &
        & -9.858662824266736e0_wp, 0.2622049408092991e0_wp, &
        & -72.67984057864001e0_wp, 11.28350883179339e0_wp, &
        & -0.5323962789000201e0_wp, -2.447540583751350e0_wp, &
        & 148.9575982135905e0_wp, 1249.939675860206e0_wp, &
        & 268.6478693121586e0_wp, 33.67689052007242e0_wp, &
        & -35.65384689269180e0_wp, 3.464021415311834e0_wp, &
        & -103.4190791219546e0_wp, 76.76961997409407e0_wp, &
        & 11.59182442022548e0_wp, -0.5254501536114633e0_wp, &
        & -32.07157614897808e0_wp, 1.102243506822152e0_wp, &
        & -952.9757997893194e0_wp, 52.11461967121677e0_wp, &
        & 6.929202177070909e0_wp, -6.322557284152793e0_wp/)
end module

!> @brief Evaluate a batch of H3 CCI geometries.
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
        call potcbs(r1x, r2x, r3x, evp(ivp))
    enddo
end subroutine pot

!> @brief Evaluate one H3 CCI geometry.
subroutine potcbs( r1x, r2x, r3x, eee )
    use ccidata,                 only: ald, beta, betad, epscem, epslon, eshift, maxi2a, maxi2b, rnought, xlin1, xlin2, xlin3, xlin4
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
end subroutine potcbs

!> @brief Evaluate trip for the H3_GEN_CCI_2002 potential.
subroutine trip( r, v )
    !     H2 triplet curve for FCI/CBS
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: j
    real(wp) :: r, v, damp, xd, prefac
    real(wp) :: c6=-6.499027e0_wp
    real(wp) :: c8=-124.3991e0_wp
    real(wp) :: c10=-3285.828e0_wp
    real(wp) :: beta=0.2e0_wp
    real(wp) :: re=1.401e0_wp
    real(wp) :: alpha=4.342902497495857e0_wp
    real(wp) :: xlp(17)=(/2.190254292077946e-1_wp, &
        & 6.903970886899540e-1_wp, 1.163805164813647e+0_wp, &
        & 1.337274224278977e+0_wp, 1.176131501286896e+0_wp, &
        & 9.093025316067311e-1_wp, 5.908490740317524e-1_wp, &
        & 1.467271248985162e-1_wp, 6.218667627370255e-2_wp, &
        & 1.105747186790804e-1_wp, -7.112594221275334e-2_wp, &
        & -8.866929977503379e-3_wp, 3.823523348718375e-2_wp, &
        & -2.055915301253834e-2_wp, 5.419433399532862e-3_wp, &
        & -7.263483303294653e-4_wp, 4.154351972583585e-5_wp/)
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

!> @brief Evaluate singlet for the H3_GEN_CCI_2002 potential.
subroutine singlet( r, v )
    !     H2 singlet potential for FCI/CBS
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: j
    real(wp) :: r, v, damp, xd, prefac
    real(wp) :: c6=-6.499027e0_wp
    real(wp) :: c8=-124.3991e0_wp
    real(wp) :: c10=-3285.828e0_wp
    real(wp) :: beta=0.2e0_wp
    real(wp) :: re=1.401e0_wp
    real(wp) :: alpha=3.980917850296971e0_wp
    real(wp) :: xlp(17)= (/ -1.709663318869429e-1_wp, &
        & -6.675286769482015e-1_wp, -1.101876055072129e+0_wp, &
        & -1.106460095658282e+0_wp, -7.414724284525231e-1_wp, &
        & -3.487882731923523e-1_wp, -1.276736255756598e-1_wp, &
        & -5.875965709151867e-2_wp, -4.030128017933840e-2_wp, &
        & -2.038653237221007e-2_wp, -7.639198558462706e-4_wp, &
        & 2.912954920483885e-3_wp, -2.628273116815280e-4_wp, &
        & -4.622088855684211e-4_wp, 1.278468948126147e-4_wp, &
        & -1.157434070240206e-5_wp, -2.609691840882097e-12_wp/)
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
