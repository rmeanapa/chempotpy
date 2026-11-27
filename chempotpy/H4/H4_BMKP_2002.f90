! System:                        H4
! Functional form:               Global analytic many-body expansion with Schwenke pair decompositions
! Common name:                   BMKP H4
! Reference:                     A. I. Boothroyd, P. G. Martin, W. J. Keogh, and M. R. Peterson, J. Chem. Phys. 116, 666-689 (2002), doi:10.1063/1.1405008
! Number of derivatives:         1
! Number of bodies:              4
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable released fit coefficients are embedded
! Notes:                         Full-dimensional H2 + H2 surface; total energy is relative to four separated H atoms.

!> @brief Implement the full-dimensional BMKP potential for the H4 system.
module chempotpy_h4_bmkp_2002_pes
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    use chempotpy_status,    only: chempotpy_status_t, status_derivative_unavailable
    use chempotpy_units,     only: bohr_to_angstrom, hartree_per_bohr_to_ev_per_angstrom, hartree_to_electronvolt
    implicit none
    private
    public :: evaluate
    type :: h4_workspace_t
    real(wp) :: ap_vpcalc_2
    real(wp) :: ap_vpcalc_3
    real(wp) :: aselect(8)
    real(wp) :: beta(8)
    real(wp) :: beta_j
    real(wp) :: beta_p
    real(wp) :: coef(2918)
    real(wp) :: dcomdcc(6,3,12)
    real(wp) :: dfsadcc(3,12)
    real(wp) :: dfsadri(3,6)
    real(wp) :: dgeodcc(3,6,12)
    real(wp) :: dgeodri(3,6,6)
    real(wp) :: dhsf_lond
    real(wp) :: dhsfinv_lond
    real(wp) :: dridcc(6,12)
    real(wp) :: drvdcc(6,3,12)
    real(wp) :: dsumvpdri(3,6)
    real(wp) :: dyydgeo(3,6,6)
    real(wp) :: eps_lond
    real(wp) :: fsa(3)
    real(wp) :: geo(3,6)
    real(wp) :: geo_sin(3,0:5)
    integer  :: i1sh(4)
    integer  :: i2sh(4)
    integer  :: icomp
    integer  :: ieqp(8)
    integer  :: igen1
    integer  :: igen2
    integer  :: ikmax_hi
    integer  :: ikmin_lo
    integer  :: imax(8)
    integer  :: imin(8)
    integer  :: indmax
    integer  :: inonlinh3
    integer  :: ipart1(6)
    integer  :: ipart2(6)
    integer  :: ir_ij(4,4)
    integer  :: isqtet
    integer  :: iswtyp_lond
    integer  :: iv(8)
    integer  :: jco_geo(6)
    integer  :: jco_ri(6,3)
    integer  :: k1_hi
    integer  :: k1_lo
    integer  :: k2_hi(12)
    integer  :: k2_lo(12)
    integer  :: k3_hi(144)
    integer  :: k3_lo(144)
    integer  :: k4_hi(1728)
    integer  :: k4_lo(1728)
    integer  :: k5_hi(3000)
    integer  :: k5_lo(3000)
    integer  :: k6_hi(3000)
    integer  :: k6_lo(3000)
    integer  :: kmax(8)
    integer  :: kmin(8)
    integer  :: mbe_perm(6,24)
    integer  :: ncol(8)
    integer  :: njac
    integer  :: nkeep(8)
    integer  :: nll(8)
    integer  :: nll4(8)
    real(wp) :: ordcoef(2048)
    real(wp) :: pt(2,3)
    real(wp) :: ptm(2)
    real(wp) :: rdelta_lond
    real(wp) :: rt_zfl0m0
    real(wp) :: rt_zfm0
    real(wp) :: rt_zfm1_3
    real(wp) :: rt_zfm2_3
    real(wp) :: s0_nsb_1
    real(wp) :: s0_nsb_2
    real(wp) :: sfhi_lond
    real(wp) :: sflo_lond
    real(wp) :: sfmid_lond
    real(wp) :: sh_nsb_1
    real(wp) :: sh_nsb_2
    real(wp) :: st_nsb_1
    real(wp) :: st_nsb_2
    real(wp) :: sumvp(3)
    real(wp) :: vlinear(3)
    real(wp) :: vlinmbe
    real(wp) :: yy(3,6)
    real(wp) :: z4pi
    real(wp) :: z4pi_rt2
    end type h4_workspace_t

    contains

    !> @brief Evaluate BMKP H4 energies and analytic Cartesian gradients.
    subroutine evaluate( x, igrad, p, g, d, status )
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms = 4, nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        type(h4_workspace_t), allocatable :: state
        real(wp) :: coordinates(natoms,3), distances(6), distance_gradient(6)
        real(wp) :: derivative_error, native_energy, native_gradient(natoms,3)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad < 0.or.igrad > 1 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'H4_BMKP_2002 provides energies and first derivatives')
            return
        endif
        allocate(state)
        call prepare_workspace(state)
        coordinates = x/bohr_to_angstrom
        call h4bmkp(state, coordinates, native_energy, distances, native_gradient, &
            distance_gradient, derivative_error, igrad)
        p(1) = native_energy*hartree_to_electronvolt
        if( igrad == 1 ) g(1,:,:) = native_gradient*hartree_per_bohr_to_ev_per_angstrom
    end subroutine evaluate

    !> @brief Prepare one evaluation-local BMKP workspace from embedded fit data.
    subroutine prepare_workspace( state )
        type(h4_workspace_t), intent(out) :: state
        state%ap_vpcalc_2  = 0.0_wp
        state%ap_vpcalc_3  = 0.0_wp
        state%aselect      = 0.0_wp
        state%beta         = 0.0_wp
        state%beta_j       = 0.0_wp
        state%beta_p       = 0.0_wp
        state%coef         = 0.0_wp
        state%dcomdcc      = 0.0_wp
        state%dfsadcc      = 0.0_wp
        state%dfsadri      = 0.0_wp
        state%dgeodcc      = 0.0_wp
        state%dgeodri      = 0.0_wp
        state%dhsf_lond    = 0.0_wp
        state%dhsfinv_lond = 0.0_wp
        state%dridcc       = 0.0_wp
        state%drvdcc       = 0.0_wp
        state%dsumvpdri    = 0.0_wp
        state%dyydgeo      = 0.0_wp
        state%eps_lond     = 0.0_wp
        state%fsa          = 0.0_wp
        state%geo          = 0.0_wp
        state%geo_sin      = 0.0_wp
        state%i1sh         = 0
        state%i2sh         = 0
        state%icomp        = 0
        state%ieqp         = 0
        state%igen1        = 0
        state%igen2        = 0
        state%ikmax_hi     = 0
        state%ikmin_lo     = 0
        state%imax         = 0
        state%imin         = 0
        state%indmax       = 0
        state%inonlinh3    = 0
        state%ipart1       = 0
        state%ipart2       = 0
        state%ir_ij        = 0
        state%isqtet       = 0
        state%iswtyp_lond  = 0
        state%iv           = 0
        state%jco_geo      = 0
        state%jco_ri       = 0
        state%k1_hi        = 0
        state%k1_lo        = 0
        state%k2_hi        = 0
        state%k2_lo        = 0
        state%k3_hi        = 0
        state%k3_lo        = 0
        state%k4_hi        = 0
        state%k4_lo        = 0
        state%k5_hi        = 0
        state%k5_lo        = 0
        state%k6_hi        = 0
        state%k6_lo        = 0
        state%kmax         = 0
        state%kmin         = 0
        state%mbe_perm     = 0
        state%ncol         = 0
        state%njac         = 0
        state%nkeep        = 0
        state%nll          = 0
        state%nll4         = 0
        state%ordcoef      = 0.0_wp
        state%pt           = 0.0_wp
        state%ptm          = 0.0_wp
        state%rdelta_lond  = 0.0_wp
        state%rt_zfl0m0    = 0.0_wp
        state%rt_zfm0      = 0.0_wp
        state%rt_zfm1_3    = 0.0_wp
        state%rt_zfm2_3    = 0.0_wp
        state%s0_nsb_1     = 0.0_wp
        state%s0_nsb_2     = 0.0_wp
        state%sfhi_lond    = 0.0_wp
        state%sflo_lond    = 0.0_wp
        state%sfmid_lond   = 0.0_wp
        state%sh_nsb_1     = 0.0_wp
        state%sh_nsb_2     = 0.0_wp
        state%st_nsb_1     = 0.0_wp
        state%st_nsb_2     = 0.0_wp
        state%sumvp        = 0.0_wp
        state%vlinear      = 0.0_wp
        state%vlinmbe      = 0.0_wp
        state%yy           = 0.0_wp
        state%z4pi         = 0.0_wp
        state%z4pi_rt2     = 0.0_wp
        state%coef(1:30) = [11.41656861561020e0_wp, 2.461981010111942e0_wp, -1.467350524873491e0_wp, 0.e0_wp, &
            & 0.4320536169168858e0_wp, -.2425867014300068e0_wp, 1.369331303380642e0_wp, -1.433069016328931e0_wp, &
            & 0.3878029623800905e0_wp, 0.2243777140165565e0_wp, 0.e0_wp, 0.e0_wp, 0.9565185801891209e0_wp, 0.e0_wp, &
            & 0.1679487106395962e0_wp, -.8043324172825887e-1_wp, -.5553208620290703e-1_wp, 0.4021797399144415e-1_wp, &
            & 0.e0_wp, 0.e0_wp, 0.1407145922487997e0_wp, 0.e0_wp, -.5119002797896408e-1_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, 0.2732447662063222e-1_wp, 0.e0_wp]
        state%coef(31:60) = [0.2569544783161231e-2_wp, 0.e0_wp, -.3351870046177481e-1_wp, 0.e0_wp, 0.e0_wp, - &
            & .2326424596451241e-2_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -.5881369223191887e-2_wp, 0.e0_wp, - &
            & .4674389804646140e-1_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.2235313229563889e-1_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, -.4422143271573974e-2_wp, 0.2313604270707399e-2_wp, 0.e0_wp, -.2588537141170540e-1_wp, &
            & 0.1826575181931046e-1_wp, 0.e0_wp, 0.e0_wp, -.1007716530727486e-1_wp, -.3307300279251005e-2_wp]
        state%coef(61:90) = [0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -.4247747992767682e-2_wp, &
            & 0.e0_wp, 0.2044170535757541e-2_wp, 0.e0_wp, 0.e0_wp, 0.1522380389286762e0_wp, 0.e0_wp, &
            & 0.2193168336511549e-1_wp, -.3864155496673826e0_wp, 0.2952821198469801e-1_wp, 0.e0_wp, - &
            & .5122925515081258e-1_wp, -.7933584692639094e-1_wp, 0.e0_wp, -.9640444763980278e-1_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, -.1766881903778219e-1_wp, -.7235812426737219e-1_wp, 0.1185383058919827e0_wp, - &
            & .1391655300384850e-1_wp]
        state%coef(91:120) = [0.e0_wp, 0.6580256079215258e-2_wp, -.6442806985257051e-2_wp, 0.e0_wp, 0.e0_wp, &
            & 0.1061238256124183e-1_wp, -.7080023799253410e-2_wp, -.9177354648147840e-2_wp, 0.e0_wp, - &
            & .1189054490030095e0_wp, 0.8131113974522991e-1_wp, -.6028550495422950e-1_wp, -.1340099467289728e-1_wp, &
            & 0.1188970731299674e-1_wp, 0.e0_wp, 0.1216426793454426e0_wp, -.8154964357971532e-1_wp, &
            & 0.6642703603111079e-1_wp, 0.1318918383813553e-1_wp, -.1392311809243365e-1_wp, 0.2084072455893501e-2_wp, &
            & - .2821790556887770e-1_wp, 0.1836263349587651e-1_wp, -.1461994658707290e-1_wp, -.2933600332178104e-2_wp, &
            & 0.3279655300405622e-2_wp, -.7536040194164841e-3_wp, 0.9662904089135244e-3_wp, 0.e0_wp, - &
            & .2353041625429331e-2_wp]
        state%coef(121:150) = [0.e0_wp, 0.5611934951252549e-3_wp, 0.9564226935654422e-3_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, -.6829845855150903e-3_wp, 0.e0_wp, 0.e0_wp, 0.1903275078627081e-3_wp, - &
            & .3308243596276508e-4_wp, 0.e0_wp, 0.7794519273561473e-4_wp, 0.e0_wp, 0.e0_wp, 0.2101122114040021e-2_wp, &
            & 0.e0_wp, -.4040602236056038e-3_wp, -.2587008024119695e-3_wp, -.6364645397740569e-3_wp, 0.e0_wp, - &
            & .9353055974199075e-3_wp, 0.e0_wp, 0.e0_wp, 0.4711266249527347e-3_wp, 0.e0_wp, 0.1569194222582758e-3_wp, &
            & 0.e0_wp]
        state%coef(151:180) = [0.e0_wp, 0.e0_wp, -.7820219543883007e-4_wp, 0.e0_wp, 0.e0_wp, 0.2074127775895557e-3_wp, &
            & 0.e0_wp, 0.e0_wp, 0.3488603138430314e-3_wp, 0.e0_wp, 0.e0_wp, -.4895717143534099e-3_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -.1371570251840013e-4_wp, 0.7871566414221949e-4_wp, - &
            & .8492275975618101e-4_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, - &
            & .1708344657789834e-2_wp]
        state%coef(181:210) = [0.e0_wp, 0.3323527691292804e-3_wp, -.7333214741878197e-4_wp, 0.e0_wp, &
            & 0.6460682421301904e-3_wp, 0.e0_wp, -.2277390762352338e-3_wp, 0.e0_wp, 0.9128512806546217e-3_wp, - &
            & .6958228400160269e-3_wp, 0.2660140098402603e-3_wp, 0.e0_wp, 0.5383770268573472e-3_wp, - &
            & .5860104570730753e-3_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -.7050220161376720e-4_wp, 0.e0_wp, - &
            & .1849834003028463e0_wp, 0.e0_wp, 0.4620812091093298e0_wp, -.1757213672839463e0_wp, 0.e0_wp, 0.e0_wp, - &
            & .4057065045231522e0_wp, 0.5750547530386685e0_wp, -.8603309454113299e-1_wp, 0.e0_wp, 0.e0_wp]
        state%coef(211:240) = [0.2344020773325997e0_wp, 0.e0_wp, -.8795349352684374e-1_wp, 0.1264880308849537e-1_wp, &
            & 0.e0_wp, 0.2668482770587186e-2_wp, 0.e0_wp, 0.e0_wp, -.1100610087674357e0_wp, 0.e0_wp, &
            & 0.4198142714188278e-1_wp, -.5309915669636895e-1_wp, 0.e0_wp, 0.8721717815328730e-2_wp, &
            & 0.4286969491463773e-1_wp, 0.e0_wp, -.8749194176570741e-2_wp, 0.2024906272110951e-1_wp, &
            & 0.9083193758180202e-2_wp, -.7903699985451101e-2_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -.2249087223226226e-2_wp, &
            & 0.e0_wp, -.7488035985516799e-1_wp, -.7884539373833625e-1_wp, 0.4784488215937220e-1_wp, 0.e0_wp, &
            & 0.7757934134790705e-2_wp]
        state%coef(241:270) = [0.5405421234966801e-1_wp, 0.8802835957276156e-1_wp, 0.e0_wp, -.4030577591563147e-1_wp, &
            & 0.e0_wp, 0.e0_wp, -.2984565152037004e-1_wp, -.1727219033185700e-1_wp, 0.1182144991915245e-1_wp, &
            & 0.7938201173244584e-2_wp, -.1949620877956762e-2_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -.3235243307954152e-3_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp]
        state%coef(271:300) = [-.7620268765924905e-1_wp, 0.2615432676679663e0_wp, 0.e0_wp, -.3733918059753137e0_wp, &
            & 0.e0_wp, 0.e0_wp, 0.2108908100274858e0_wp, -.3908497880553245e-1_wp, -.4135382830384930e-1_wp, &
            & 0.1522872914884577e0_wp, -.2855267631463690e0_wp, 0.1187218798350170e-1_wp, 0.e0_wp, &
            & 0.1720929635828025e0_wp, 0.e0_wp, -.4762896947490703e-1_wp, 0.e0_wp, 0.1296479727681136e-1_wp, 0.e0_wp, &
            & 0.3875683711725621e-1_wp, 0.e0_wp, 0.e0_wp, -.2728020475068825e-1_wp, -.3568966505779430e-2_wp, &
            & 0.2391400007431497e-2_wp, 0.e0_wp, 0.e0_wp, 0.7092715349016440e-1_wp, -.2171083138201483e-1_wp, 0.e0_wp]
        state%coef(301:330) = [0.e0_wp, 0.3703995087126884e-2_wp, -.2434327063040666e-2_wp, -.9427841028470206e-1_wp, &
            & 0.e0_wp, 0.e0_wp, 0.1294971819195222e-2_wp, 0.9682682112195861e-3_wp, 0.e0_wp, 0.1631275239132733e-1_wp, &
            & 0.3636984154327500e-2_wp, -.1340984765576260e-2_wp, 0.e0_wp, -.5956593779898524e-3_wp, &
            & 0.2314568357119641e-3_wp, 0.e0_wp, 0.e0_wp, -.1992197266813824e-2_wp, 0.e0_wp, 0.e0_wp, - &
            & .1694357150309008e-3_wp, 0.e0_wp, 0.7029485453338565e-3_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.1070975861800831e-3_wp, -.2113079739386820e-3_wp, 0.e0_wp, 0.e0_wp]
        state%coef(331:360) = [-.2539283910284728e-4_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.1099815327496916e-2_wp, - &
            & .1765748911285442e-2_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.1998984184634926e-4_wp, &
            & 0.e0_wp, 0.e0_wp, 0.e0_wp, -.6705593097363792e-4_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, -.7370107278577147e-4_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp]
        state%coef(361:390) = [0.e0_wp, 0.e0_wp, -.1226471127884504e-3_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp, 0.3602375326308723e-4_wp, -.3519096757861601e-2_wp, -.1411948724683564e-2_wp, 0.e0_wp, - &
            & .1573773984487967e-1_wp, 0.1185480315550147e-1_wp, 0.e0_wp, 0.5103394787098075e-2_wp, - &
            & .3614550200861926e-2_wp, -.1770324228084872e-3_wp, 0.e0_wp, -.1992171779083550e-2_wp, 0.e0_wp, &
            & 0.1270692393828489e-1_wp, -.4921205014438497e-2_wp, 0.2064885831382025e-3_wp, -.2060673869866738e-2_wp, &
            & 0.1474585053339496e-2_wp, 0.e0_wp, 0.e0_wp, 0.9978301872727278e-3_wp, -.3964660536460676e-4_wp]
        state%coef(391:420) = [-.1773396097672859e-2_wp, 0.e0_wp, 0.e0_wp, 0.1172919080826781e-3_wp, - &
            & .9737749811915252e-4_wp, 0.e0_wp, -45.26288823530060e0_wp, -8.088211185381054e0_wp, &
            & 3.668642918880288e0_wp, 0.4786804751301441e0_wp, -1.707321870517821e0_wp, 1.088841015252361e0_wp, - &
            & 3.303208576830779e0_wp, 3.936877551205625e0_wp, 0.e0_wp, 0.4179110229268445e0_wp, 0.e0_wp, - &
            & .3959245897503227e0_wp, 3.770841516697761e0_wp, 0.e0_wp, -1.114631944449590e0_wp, - &
            & .8088285725825271e-1_wp, 0.1566691026629507e0_wp, -.1855320817030686e-1_wp, 0.e0_wp, - &
            & .6410649862995392e-1_wp, -.3266963544501503e0_wp, 0.4155996356413736e-1_wp, 0.e0_wp, 0.e0_wp]
        state%coef(421:450) = [0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, - &
            & .5501384523353580e-2_wp, 0.1060007442378859e0_wp, -.1608713465509136e-1_wp, 0.1437949646266004e-1_wp, &
            & 0.e0_wp, 0.e0_wp, 0.3394636523584611e-1_wp, -.1134732985014914e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.1277036538813291e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.8296264006392355e-2_wp, -.3992780588062202e-1_wp, &
            & -.2201005099651295e-1_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.1166889972496094e-1_wp, - &
            & .5709122448787395e-2_wp]
        state%coef(451:480) = [0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.5276267894590526e-1_wp, &
            & 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.5828792976399929e-2_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, - .6901495268242744e0_wp, 0.e0_wp, 0.e0_wp, 1.463682894109018e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.6284874210232054e0_wp, 0.7177356034107231e-1_wp, 0.1114026925935851e0_wp, 0.7639309350060595e0_wp, &
            & 0.e0_wp, 0.e0_wp]
        state%coef(481:510) = [0.e0_wp, 0.e0_wp, 0.e0_wp, -.3105415152600133e0_wp, -.1846415240159753e0_wp, - &
            & .1857352512650449e-1_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -.1474876842579986e0_wp, -.3788771252443011e-1_wp, &
            & 0.1371285727211139e-1_wp, 0.e0_wp, 0.6625356371830053e-1_wp, -.1026112057095907e-1_wp, 0.e0_wp, - &
            & .1105251492579490e-5_wp, 0.e0_wp, 0.e0_wp, -.6418485229163887e-6_wp, 0.1132970198243758e-5_wp, &
            & 0.3088336285550473e-5_wp, 0.6006136636933950e-6_wp, 0.3965148926799415e-6_wp, 0.e0_wp, 0.e0_wp, - &
            & .4447322393731460e-6_wp, -.1699716222759334e-5_wp, 0.e0_wp, 0.e0_wp]
        state%coef(511:540) = [0.4528105444195479e-7_wp, 0.e0_wp, 0.4873262657912689e-7_wp, 0.e0_wp, - &
            & .1069438713799491e-5_wp, 0.1111275907252112e-5_wp, 0.e0_wp, 0.6396482750065100e-7_wp, 0.e0_wp, &
            & 0.2620093358105838e-5_wp, 0.e0_wp, -.6800012598833259e-6_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, - &
            & .1859877217865141e-5_wp, 0.5267767186786455e-6_wp, 0.e0_wp, -.4490791393108958e-7_wp, 0.e0_wp, 0.e0_wp, &
            & - .1808544844425166e-5_wp, 0.e0_wp, 0.1503227604390245e-5_wp, -.8450361230700443e-7_wp, 0.e0_wp, &
            & 0.e0_wp, 0.1594595406051139e-5_wp, 0.e0_wp, -.1658370074628046e-5_wp]
        state%coef(541:570) = [0.e0_wp, 0.1504941865060415e-6_wp, -.1372299792933966e-7_wp, -.8036673503006943e-6_wp, &
            & 0.4019859576399040e-6_wp, 0.e0_wp, -.4451491359093753e-7_wp, 0.e0_wp, 0.e0_wp, 0.4390555540312022e-5_wp, &
            & - .1988249891489124e-5_wp, 0.e0_wp, -.2560677863118263e-6_wp, 0.9815971961161848e-6_wp, - &
            & .8183070027992053e-6_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.1870275113972953e-6_wp, 0.e0_wp, &
            & 0.e0_wp, 0.e0_wp, 0.1899881667341878e-7_wp, -.4642199767652417e-7_wp, 0.e0_wp, -.3818998927536044e-5_wp, &
            & 0.e0_wp, 0.e0_wp]
        state%coef(571:600) = [0.4643005856103144e-5_wp, 0.e0_wp, 0.e0_wp, -.6998860287758503e-6_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp, 0.1159361237893339e-5_wp, 0.e0_wp, -.3532821151131939e-5_wp, 0.e0_wp, -.1119561115279527e-6_wp, &
            & 0.6730732445924999e-6_wp, 0.e0_wp, 0.e0_wp, 0.1743906906514087e-5_wp, -.1011848690997574e-5_wp, &
            & 0.8010838221017618e-7_wp, 0.1751392240339859e-6_wp, 0.1887909194821102e-6_wp, 0.1006090339987867e-7_wp, &
            & - .1179071341543131e-6_wp, 0.e0_wp, 0.e0_wp, -.2577636134332770e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 0.2950066706001413e-2_wp, -.8436249003701826e-3_wp]
        state%coef(601:630) = [0.9650392433268715e-1_wp, 0.1180865411962495e-1_wp, -.1890768945728852e-1_wp, 0.e0_wp, &
            & 0.e0_wp, 0.1396056667617217e-3_wp, -.4784573309839698e-2_wp, -.2307816824880044e-2_wp, &
            & 0.2236248930560251e-2_wp, 0.2741664444144839e-3_wp, -.1050250271971434e-3_wp, 0.e0_wp, 0.e0_wp, - &
            & .1319897018744042e-2_wp, -.4284037762263289e-2_wp, 0.9989080207036391e-3_wp, 0.e0_wp, - &
            & .1884447664581848e-3_wp, 0.1495754666678066e-3_wp, 0.e0_wp, 0.4216827533402400e-2_wp, - &
            & .5873115173081413e-3_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -.4041229946271942e-3_wp, &
            & 0.6443324536316071e-4_wp, -.8167457890020569e-5_wp, 0.1056967699024861e-4_wp]
        state%coef(631:660) = [0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.2803923282671573e-3_wp, &
            & 0.1255823209120723e-3_wp, 0.e0_wp, 0.e0_wp, 0.3631834621471845e-4_wp, -.1237743740137913e-3_wp, - &
            & .7679069449277003e-4_wp, -.5544353440897096e-4_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.2226610345487547e-4_wp, &
            & 0.e0_wp, 0.1102858209688644e-2_wp, -.9707421030189208e-4_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, - &
            & .1788348599144874e-4_wp, -.8680038386835745e-4_wp, 0.e0_wp, 0.7502974676591654e-4_wp, 0.e0_wp, 0.e0_wp, &
            & 0.e0_wp]
        state%coef(661:690) = [0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -.1140689346615696e-5_wp, &
            & 0.1338360967887129e-1_wp, -.4734314575965839e-2_wp, 0.e0_wp, 0.1164975543347115e-1_wp, &
            & 0.2368768109003181e-2_wp, 0.2801552358660103e-3_wp, 0.7028567446182806e-3_wp, -.3124827789848925e-3_wp, &
            & - .3151243066290188e-3_wp, 0.e0_wp, 0.5887365420884373e-3_wp, 0.e0_wp, -.8018026247461792e-2_wp, &
            & 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.8600408000548635e-4_wp, -.2211835739470488e-3_wp, 0.e0_wp, &
            & 0.3075426445130659e-4_wp, 0.8401140913442456e-3_wp, 0.e0_wp, -.4243877152533952e-4_wp]
        state%coef(691:720) = [0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 4.570473204155594e0_wp, -.2798626320494389e0_wp, - &
            & 6.940100755011418e0_wp, 0.e0_wp, 64.98997989641670e0_wp, 1.715566496848222e0_wp, - &
            & 7.759401344839803e0_wp, - 7.307106314539124e0_wp, -11.18609692968977e0_wp, -92.20698555542199e0_wp, &
            & 29.72626438699660e0_wp, 11.98208769527558e0_wp, -39.24233472335177e0_wp, 137.9163648926951e0_wp, &
            & 65.69977054586563e0_wp, 0.e0_wp, 214.0380086983486e0_wp, 55.46915239130455e0_wp, - &
            & 166.2744780226556e0_wp, 0.e0_wp, 40.81904694776546e0_wp, 0.e0_wp, -62.76548193762972e0_wp, - &
            & 37.71415377213584e0_wp, -167.4967713102476e0_wp, -345.0583498103054e0_wp]
        state%coef(721:750) = [-192.0902057087892e0_wp, -79.86702739644112e0_wp, 0.e0_wp, 0.e0_wp, &
            & 413.0483295509064e0_wp, 0.e0_wp, 39.64091918839881e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 96.66100850781979e0_wp, - 169.3760884098744e0_wp, 0.e0_wp, 573.3282420428860e0_wp, - &
            & 3.671947638444790e0_wp, 25.53422541109841e0_wp, - 19.08478197765944e0_wp, 14.20835333585394e0_wp, &
            & 0.e0_wp, 283.3923273325448e0_wp, -217.6861835311140e0_wp, - 62.28512734587679e0_wp, - &
            & 62.71380848477281e0_wp, -21.74157883889574e0_wp, 86.04408325527467e0_wp, - 323.5339117089853e0_wp, &
            & 260.0626158828048e0_wp, -103.0644269330105e0_wp]
        state%coef(751:780) = [-271.2087580376460e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, &
            & 69.07644925729612e0_wp, -52.23384556701277e0_wp, 0.e0_wp, 228.6493095830203e0_wp, &
            & 124.2085335019347e0_wp, - 85.21025019384091e0_wp, 496.0374997733861e0_wp, -684.4131378969931e0_wp, &
            & 415.4879765735129e0_wp, 890.6952128611791e0_wp, -387.9510910233527e0_wp, 230.2546692076336e0_wp, - &
            & 101.1696229503429e0_wp, 0.e0_wp, - 118.3505341147770e0_wp, 0.e0_wp, -734.6888505805208e0_wp, 0.e0_wp, - &
            & 387.8540736616000e0_wp, 0.e0_wp, 0.e0_wp, -249.9596065837551e0_wp, 0.e0_wp, 0.e0_wp]
        state%coef(781:788) = [-99.09862868012513e0_wp, 0.e0_wp, 0.e0_wp, 0.e0_wp, -155.3886754928112e0_wp, &
            & 563.8830985327528e0_wp, -94.35505172463631e0_wp, 189.7518257337113e0_wp]
        state%beta   = [0.006e0_wp, 0.006e0_wp, 0.006e0_wp, 0.006e0_wp, 0.006e0_wp, 0.003e0_wp, 0.006e0_wp, 0.0e0_wp]
        state%iv     = [2, 3, 3, 3, 3, 3, 3, 100]
        state%imin   = 0
        state%imax   = [2, 2, 2, 2, 2, 2, 2, 3]
        state%kmin   = 0
        state%kmax   = [2, 2, 2, 2, 2, 2, 2, 8]
        state%njac   = 7
        state%beta_p = 0.8e0_wp
        state%mbe_perm = reshape([1, 2, 3, 4, 5, 6, 1, 3, 2, 5, 4, 6, 2, 1, 3, 4, 6, 5, 2, 3, 1, 6, 4, 5, 3, 1, 2, 5, &
            & 6, 4, 3, 2, 1, 6, 5, 4, 1, 4, 5, 2, 3, 6, 1, 5, 4, 3, 2, 6, 4, 1, 5, 2, 6, 3, 4, 5, 1, 6, 2, 3, 5, 1, 4, &
            & 3, 6, 2, 5, 4, 1, 6, 3, 2, 2, 4, 6, 1, 3, 5, 2, 6, 4, 3, 1, 5, 4, 2, 6, 1, 5, 3, 4, 6, 2, 5, 1, 3, 6, 2, &
            & 4, 3, 5, 1, 6, 4, 2, 5, 3, 1, 3, 5, 6, 1, 2, 4, 3, 6, 5, 2, 1, 4, 5, 3, 6, 1, 4, 2, 5, 6, 3, 4, 1, 2, 6, &
            & 3, 5, 2, 4, 1, 6, 5, 3, 4, 2, 1], shape(state%mbe_perm))
        state%indmax      = 3
        state%k1_lo       = 1
        state%k1_hi       = 3
        state%k2_lo       = 1
        state%k2_hi(1:12) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        state%k3_lo       = 0
        state%k3_hi(1:6)  = [1, 1, 2, 1, 2, 2]
        state%k4_lo       = 0
        state%k4_hi(1:15) = [0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 0]
        state%k5_lo(1:32) = [1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, &
            & 0, 0]
        state%k5_hi(1:32) = [1, 0, 1, 1, 0, 0, 1, 2, 1, 2, 2, 1, 0, 1, 0, 1, 0, 0, 1, 2, 1, 0, 2, 1, 0, 0, 0, 2, 1, 1, &
            & 0, 0]
        state%k6_lo(1:52) = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, &
            & 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        state%k6_hi(1:52) = [1, 0, 0, 1, 2, 2, 1, 2, 2, 2, 1, 2, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 3, 2, 3, 2, &
            & 1, 3, 2, 1, 2, 1, 1, 2, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0]
        state%beta_j       = 10.e0_wp
        state%eps_lond     = 0.775e0_wp
        state%sflo_lond    = 0.e0_wp
        state%sfhi_lond    = 7.5e0_wp
        state%sfmid_lond   = 3.75e0_wp
        state%dhsf_lond    = 3.75e0_wp
        state%dhsfinv_lond = (1.0_wp/3.75_wp)
        state%rdelta_lond  = 3.e0_wp
        state%iswtyp_lond  = 4
        state%s0_nsb_1     = 10.04e0_wp
        state%st_nsb_1     = 0.312e0_wp
        state%sh_nsb_1     = 1.61e0_wp
        state%s0_nsb_2     = 4.03e0_wp
        state%st_nsb_2     = 0.53e0_wp
        state%sh_nsb_2     = 1.079e0_wp
        state%ap_vpcalc_2  = 0.2593e0_wp
        state%ap_vpcalc_3  = 1.6813e0_wp
        state%ipart1       = [1, 1, 1, 2, 2, 3]
        state%ipart2       = [2, 3, 4, 3, 4, 4]
        state%ir_ij        = reshape([0, 1, 2, 3, 1, 0, 4, 5, 2, 4, 0, 6, 3, 5, 6, 0], shape(state%ir_ij))
        state%geo          = 0.e0_wp
        state%dridcc       = 0.e0_wp
        state%dgeodcc      = 0.e0_wp
        state%fsa          = 0.e0_wp
        state%sumvp        = 0.e0_wp
        state%dfsadcc      = 0.e0_wp
        state%dfsadri      = 0.e0_wp
        state%dsumvpdri    = 0.e0_wp
        state%yy           = 0.e0_wp
        state%dyydgeo      = 0.e0_wp
        state%dgeodri      = 0.e0_wp
        state%jco_ri       = reshape([1, 2, 3, 4, 5, 6, 2, 1, 3, 4, 6, 5, 3, 1, 2, 5, 6, 4], shape(state%jco_ri))
        state%jco_geo      = [3, 1, 2, 4, 5, 6]
        state%geo_sin      = 0.e0_wp
        state%drvdcc       = 0.e0_wp
        state%dcomdcc      = 0.e0_wp
        state%isqtet       = 3
        state%igen1        = 5
        state%icomp        = 0
        state%inonlinh3    = 6
        state%igen2        = 8
        state%i1sh         = [1, 1, 2, 2]
        state%i2sh         = [3, 4, 3, 4]
        call build_polynomial_tables(state)
    end subroutine prepare_workspace

    !> @brief Evaluate the released BMKP energy and Cartesian or distance derivatives.
    subroutine h4bmkp( state, cc, vfit, r, dvdcc, dvdr, derr, ideriv )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: cc(4,3)
        real(wp), intent(inout) :: vfit
        real(wp), intent(inout) :: r(6)
        real(wp), intent(inout) :: dvdcc(4,3)
        real(wp), intent(inout) :: dvdr(6)
        real(wp), intent(inout) :: derr
        integer, intent(in)     :: ideriv
        integer                 :: ii, j, jj, k
        real(wp)                :: dcc
        real(wp), parameter     :: tol_derr = 1.e-6_wp
        derr = 0.e0_wp
        call surftry_h4bmkp(state, cc, vfit, r, dvdcc, dvdr, ideriv)
        if( vfit .gt. 98.e0_wp )then
            derr = 99.e0_wp
            return
        endif
        if( ideriv .gt. 1 )then
            k = 0
            do jj = 1, 3
                do ii = 1, 4
                    k   = k + 1
                    dcc = 0.e0_wp
                    do j = 1, 6
                        dcc = dcc + dvdr(j) * state%dridcc(j,k)
                    enddo
                    derr = max( derr, abs( dvdcc(ii,jj) - dcc ) )
                enddo
            enddo
            if( derr .gt. tol_derr ) call fixdvdr_h4bmkp(state, cc, r, dvdcc, dvdr, derr, ideriv)
        endif
    end subroutine h4bmkp

    !> @brief Stabilize distance derivatives near geometrically singular configurations.
    subroutine fixdvdr_h4bmkp( state, cc, r, dvdcc, dvdr, derr, ideriv )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: cc(4,3)
        real(wp), intent(in)    :: r(6)
        real(wp), intent(in)    :: dvdcc(4,3)
        real(wp), intent(inout) :: dvdr(6)
        real(wp), intent(inout) :: derr
        integer, intent(in)     :: ideriv
        integer                 :: i, i1, i2, ia(4), iden, ii, ipt, j, j1, j2, j3, jj, k, k12, k13, k24, k34, klong, kloop, kloop_hi, knot, nok_1, nok_2, nuse
        real(wp)                :: ccsh(4,3), cosphi, d, dcc, dcci(3), dersh(9,3), dershmin_1, dershmin_2, dlim, dpdrrv(3), dpdrvp(2,3), dridcc_0(6,12), dsh_1, dsh_2, dsh_at, dsh_hi, dvdcc_0(12), dvdccsh(12), dvdrsh(6,9,3), esh, f_1, f_2, rmax, rrv(3), rsh(6), rvp(2,3)
        real(wp), parameter     :: v_big    = 1.e+36_wp
        real(wp), parameter     :: v_bigger = 1.01e+36_wp
        real(wp), parameter     :: tol_derr = 1.e-6_wp
        real(wp), parameter     :: tol_pt   = 1.e-7_wp
        real(wp), parameter     :: dsh_lo   = 1.e-7_wp
        rmax  = r(1)
        klong = 1
        do k = 2, 6
            if( r(k) .gt. r(klong) ) klong = k
        enddo
        ia(1) = state%ipart1(klong)
        ia(4) = state%ipart2(klong)
        knot  = 7 - klong
        ia(2) = state%ipart1(knot)
        ia(3) = state%ipart2(knot)
        k12   = state%ir_ij(ia(1),ia(2))
        k13   = state%ir_ij(ia(1),ia(3))
        k24   = state%ir_ij(ia(2),ia(4))
        k34   = state%ir_ij(ia(3),ia(4))
        if( r(k13) .le. r(k12) .and. r(k24) .le. r(k34) )then
            ia(3) = state%ipart1(knot)
            ia(2) = state%ipart2(knot)
            k12   = state%ir_ij(ia(1),ia(2))
            k13   = state%ir_ij(ia(1),ia(3))
            k24   = state%ir_ij(ia(2),ia(4))
            k34   = state%ir_ij(ia(3),ia(4))
        endif
        do i = 1, 3
            rrv(i)   = cc(ia(4),i) - cc(ia(1),i)
            rvp(1,i) = cc(ia(2),i) - cc(ia(1),i)
            rvp(2,i) = cc(ia(4),i) - cc(ia(3),i)
        enddo
        call cross_h4bmkp(state, rvp, rrv, cosphi, dpdrvp, dpdrrv, 0, iden)
        if( state%ptm(1) .ge. state%ptm(2) )then
            ipt = 1
        else
            ipt = 2
        endif
        if( state%ptm(ipt) .ge. tol_pt )then
            do i = 1, 3
                dcci(i) = state%pt(ipt,i) / state%ptm(ipt)
            enddo
        else
            if( abs(rrv(1)) .ge. max( abs(rrv(2)), abs(rrv(3)) ) )then
                j1 = 1
                if( abs(rrv(2)) .ge. abs(rrv(3)) )then
                    j2 = 2
                    j3 = 3
                else
                    j2 = 3
                    j3 = 2
                endif
            elseif( abs(rrv(2)) .ge. abs(rrv(3)) )then
                j1 = 2
                if( abs(rrv(1)) .ge. abs(rrv(3)) )then
                    j2 = 1
                    j3 = 3
                else
                    j2 = 3
                    j3 = 1
                endif
            else
                j1 = 3
                if( abs(rrv(1)) .ge. abs(rrv(2)) )then
                    j2 = 1
                    j3 = 2
                else
                    j2 = 2
                    j3 = 1
                endif
            endif
            d        = sqrt( rrv(j1)**2 + rrv(j2)**2 )
            dcci(j1) = - rrv(j2) / d
            dcci(j2) = rrv(j1) / d
            dcci(j3) = 0.e0_wp
        endif
        j = 0
        do jj = 1, 3
            do ii = 1, 4
                j           = j + 1
                dvdcc_0(j)  = dvdcc(ii,jj)
                ccsh(ii,jj) = cc(ii,jj)
            enddo
        enddo
        do i = 1, 6
            do j = 1, 12
                dridcc_0(i,j) = state%dridcc(i,j)
            enddo
        enddo
        dsh_at     = dsh_lo
        dsh_1      = 0.e0_wp
        dsh_2      = dsh_lo
        dershmin_1 = v_big
        dershmin_2 = v_bigger
        dsh_hi     = min( 1.e0_wp, 0.001e0_wp * r(klong), 0.1e0_wp * min( r(knot), r(k12), r(k13), r(k24), r(k34) ) )
        nok_1      = 0
        nok_2      = 0
        kloop      = 0
        kloop_hi   = 30
        do while ( dsh_at .le. dsh_hi .and. ( nok_1 .eq. 0 .or. nok_2 .eq. 0 .or. dershmin_2 .lt. dershmin_1 ) .and. ( kloop .lt. kloop_hi .or. dershmin_1 .gt. 1.e-3_wp ) )
            kloop = kloop + 1
            if( kloop .eq. 2 .or. dershmin_2 .lt. dershmin_1 )then
                dsh_1      = dsh_2
                dershmin_1 = dershmin_2
                nok_1      = nok_2
                do k = 1, 9
                    dersh(k,1) = dersh(k,2)
                    do j = 1, 6
                        dvdrsh(j,k,1) = dvdrsh(j,k,2)
                    enddo
                enddo
                if( dershmin_1 .lt. 1.e-5_wp )then
                    kloop_hi = kloop + 5
                elseif( dershmin_1 .lt. 1.e-4_wp )then
                    kloop_hi = kloop + 7
                elseif( dershmin_1 .lt. 1.e-3_wp )then
                    kloop_hi = kloop + 10
                else
                    kloop_hi = 30
                endif
            endif
            nok_2 = 0
            dsh_2 = dsh_at
            do k = 1, 8
                ccsh = cc
                select case (k)
                case (1:4)
                    i1 = k
                    i2 = 0
                case default
                    i1 = ia(state%i1sh(max(1, k - 4)))
                    i2 = ia(state%i2sh(max(1, k - 4)))
                end select
                ccsh(i1,:)              = cc(i1,:) + dcci(:) * dsh_at
                if( i2 > 0 ) ccsh(i2,:) = cc(i2,:) + dcci(:) * dsh_at
                call surftry_h4bmkp(state, ccsh, esh, rsh, dvdccsh, dvdrsh(1,k,2), 2)
                dersh(k,2) = 0.e0_wp
                do j = 1, 12
                    dcc = 0.e0_wp
                    do i = 1, 6
                        dcc = dcc + dvdrsh(i,k,2) * dridcc_0(i,j)
                    enddo
                    dersh(k,2) = max( dersh(k,2), abs( dvdcc_0(j) - dcc ) )
                enddo
                if( k .eq. 1 )then
                    dershmin_2 = dersh(1,2)
                else
                    dershmin_2 = min( dershmin_2, dersh(k,2) )
                endif
                if( dersh(k,2) .le. tol_derr ) nok_2 = nok_2 + 1
            enddo
            do i = 1, 6
                dvdrsh(i,9,2) = 0.e0_wp
            enddo
            nuse = 0
            dlim = 10.e0_wp * dershmin_2
            do k = 1, 8
                if( dersh(k,2) .le. dlim )then
                    nuse = nuse + 1
                    do i = 1, 6
                        dvdrsh(i,9,2) = dvdrsh(i,9,2) + dvdrsh(i,k,2)
                    enddo
                endif
            enddo
            do i = 1, 6
                dvdrsh(i,9,2) = dvdrsh(i,9,2) / nuse
            enddo
            dersh(9,2) = 0.e0_wp
            do j = 1, 12
                dcc = 0.e0_wp
                do i = 1, 6
                    dcc = dcc + dvdrsh(i,9,2) * dridcc_0(i,j)
                enddo
                dersh(9,2) = max( dersh(9,2), abs( dvdcc_0(j) - dcc ) )
            enddo
            if( dersh(9,2) .le. tol_derr ) nok_2 = nok_2 + 1
            dsh_at = dsh_at * 2.e0_wp
        enddo
        f_2 = - dsh_1 / ( dsh_2 - dsh_1 )
        f_1 = 1.e0_wp - f_2
        do k = 1, 9
            do i = 1, 6
                dvdrsh(i,k,3) = f_1 * dvdrsh(i,k,1) + f_2 * dvdrsh(i,k,2)
            enddo
            dersh(k,3) = 0.e0_wp
            do j = 1, 12
                dcc = 0.e0_wp
                do i = 1, 6
                    dcc = dcc + dvdrsh(i,k,3) * dridcc_0(i,j)
                enddo
                dersh(k,3) = max( dersh(k,3), abs( dvdcc_0(j) - dcc ) )
            enddo
            do j = 1, 3
                if( dersh(k,j) .lt. derr )then
                    derr = dersh(k,j)
                    do i = 1, 6
                        dvdr(i) = dvdrsh(i,k,j)
                    enddo
                endif
            enddo
        enddo
        if( ideriv .gt. 2 )then
            call surftry_h4bmkp(state, cc, esh, rsh, dvdccsh, dvdrsh(1,1,1), 2)
        endif
    end subroutine fixdvdr_h4bmkp

    !> @brief Assemble pair, three-body, linear, and four-body BMKP contributions.
    subroutine surftry_h4bmkp( state, cc, vfit, r, dvdcc, dvdr, ideriv )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: cc(4,3)
        real(wp), intent(inout) :: vfit
        real(wp), intent(inout) :: r(6)
        real(wp), intent(inout) :: dvdcc(4,3)
        real(wp), intent(inout) :: dvdr(6)
        integer, intent(in)     :: ideriv
        integer                 :: i, ierr, ig, ii, j, jj, k
        real(wp)                :: dcc, dvlindgeo(3,6), dvlindri(3,6), dvlondfsa(3), dvlondri(6), dvmbedri(6), dvnsbdri(6), vc_tmp, vlonh4, vnsb, vr_tmp
        do i = 1, 4
            dvdcc(i,1) = 0.e0_wp
            dvdcc(i,2) = 0.e0_wp
            dvdcc(i,3) = 0.e0_wp
        enddo
        do i = 1, 6
            dvdr(i) = 0.e0_wp
        enddo
        call geomet_h4bmkp(state, r, cc, ideriv, ierr)
        if( ierr .gt. 0 )then
            vfit = 99.e0_wp
            return
        endif
        call vcorrec_h4bmkp(state, ideriv)
        call vpcalc_h4bmkp(state, r, ideriv)
        call getyy_h4bmkp(state, ideriv)
        call comlon_h4bmkp(state, r, ideriv, vlonh4, dvlondri, dvlondfsa)
        call vnsb_h4bmkp(state, r, ideriv, vnsb, dvnsbdri)
        call vlinear_h4bmkp(state, r, ideriv, state%vlinear, state%vlinmbe, dvlindri, dvlindgeo, dvmbedri)
        vfit = vlonh4 + vnsb + state%fsa(1) * state%vlinear(1) + state%fsa(2) * state%vlinear(2) + state%fsa(3) * &
            & state%vlinear(3) + state%vlinmbe
        if( ideriv .gt. 0 )then
            j = 0
            do jj = 1, 3
                do ii = 1, 4
                    j   = j + 1
                    dcc = 0.e0_wp
                    do ig = 1, 3
                        if( state%fsa(ig) .gt. 0.e0_wp )then
                            dcc = dcc + ( state%vlinear(ig) + dvlondfsa(ig) ) * state%dfsadcc(ig,j)
                            do k = 1, 6
                                dcc = dcc + ( dvlindgeo(ig,k) * state%dgeodcc(ig,k,j) + dvlindri(ig,k) * &
                                    & state%dridcc(k,j) ) * state%fsa(ig)
                            enddo
                        endif
                    enddo
                    do k = 1, 6
                        dcc = dcc + ( dvlondri(k) + dvnsbdri(k) + dvmbedri(k) ) * state%dridcc(k,j)
                    enddo
                    dvdcc(ii,jj) = dcc
                enddo
            enddo
            if( ideriv .gt. 1 )then
                do j = 1, 6
                    dcc = 0.e0_wp
                    do ig = 1, 3
                        if( state%fsa(ig) .gt. 0.e0_wp )then
                            dcc = dcc + ( state%vlinear(ig) + dvlondfsa(ig) ) * state%dfsadri(ig,j) + dvlindri(ig,j) * &
                                & state%fsa(ig)
                            do k = 1, 6
                                dcc = dcc + dvlindgeo(ig,k) * state%dgeodri(ig,k,j) * state%fsa(ig)
                            enddo
                        endif
                    enddo
                    dvdr(j) = dcc + dvlondri(j) + dvnsbdri(j) + dvmbedri(j)
                enddo
            endif
        endif
        vr_tmp = dvdr(1)
        vc_tmp = dvdcc(1,1)
    end subroutine surftry_h4bmkp

    !> @brief Derive permutation and polynomial tables from immutable fit parameters.
    subroutine build_polynomial_tables( state )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        integer             :: i, ij, ikeep, j, ji, k, kshift, l, lshift, mshift, n, n_terms, ni, nisq, nk, nprev, nqqm, nshift
        real(wp), parameter :: z3     = 3.e0_wp
        real(wp), parameter :: z4     = 4.e0_wp
        real(wp), parameter :: z5     = 5.e0_wp
        real(wp), parameter :: z24    = 24.e0_wp
        real(wp), parameter :: z96    = 96.e0_wp
        real(wp), parameter :: zfm0   = z5/(z4*PI)
        real(wp), parameter :: zfm1   = z5/(z24*PI)
        real(wp), parameter :: zfm2   = z5/(z96*PI)
        real(wp), parameter :: zfl0m0 = 1.e0_wp/(4.e0_wp*PI)
        real(wp), parameter :: z2     = 2.e0_wp
        integer, parameter  :: nterms = 8
        state%z4pi      = z4 * PI
        state%z4pi_rt2  = state%z4pi * sqrt( z2 )
        state%rt_zfl0m0 = sqrt( zfl0m0 )
        state%rt_zfm0   = sqrt( zfm0 )
        state%rt_zfm1_3 = sqrt( zfm1 ) * z3
        state%rt_zfm2_3 = sqrt( zfm2 ) * z3
        if( state%iv(8) .le. 3 )then
            state%njac = nterms
        else
            state%njac = nterms - 1
        endif
        if( state%isqtet .gt. 0 )then
            if( state%iv(state%isqtet) .le. 0 ) state%isqtet = 0
        endif
        if( state%igen1 .gt. 0 )then
            if( state%iv(state%igen1) .le. 0 ) state%igen1 = 0
        endif
        if( state%icomp .gt. 0 )then
            if( state%iv(state%icomp) .le. 0 ) state%icomp = 0
        endif
        if( state%inonlinh3 .gt. 0 )then
            if( state%iv(state%inonlinh3) .le. 0 ) state%inonlinh3 = 0
        endif
        if( state%igen2 .gt. 0 )then
            if( state%iv(state%igen2) .le. 0 ) state%igen2 = 0
        endif
        state%ikmax_hi = -9
        state%ikmin_lo = 9
        nshift         = 0
        mshift         = 0
        nprev          = 0
        n_terms        = 0
        do n = 1, state%njac
            if( state%iv(n) .le. 0 )then
                state%ncol(n) = 0
                state%nll4(n) = 0
                state%ieqp(n) = 0
            else
                state%ieqp(n) = 1
                if( nprev .gt. 0 )then
                    if( state%imin(n) .ne. state%imin(nprev) .or. state%imax(n) .ne. state%imax(nprev) .or. state%kmin(n) .ne. state%kmin(nprev) .or. state%kmax(n) .ne. state%kmax(nprev) ) state%ieqp(n) = 0
                endif
                state%ikmin_lo = min( state%ikmin_lo, state%imin(n), state%kmin(n) )
                state%ikmax_hi = max( state%ikmax_hi, state%imax(n), state%kmax(n) )
                nk             = state%kmax(n) - state%kmin(n) + 1
                ni             = state%imax(n) - state%imin(n) + 1
                nisq           = ni * ni
                state%nll(n)   = nk * nisq
                state%nll4(n)  = 4 * state%nll(n)
                nqqm           = ( nk * ni * ( ni + 1 ) ) / 2
                state%nkeep(n) = nqqm * 4
                state%ncol(n)  = state%nkeep(n) + state%nll(n)
                n_terms        = n_terms + state%ncol(n)
                ikeep          = nshift
                lshift         = 0
                do l = 1, 4
                    kshift = lshift
                    do k = 1, nk
                        do i = 1, ni
                            do j = 1, i
                                ij                       = (i-1)*ni + j + kshift
                                ikeep                    = ikeep + 1
                                state%ordcoef(ij+mshift) = state%coef(ikeep)
                                if( i .ne. j )then
                                    ji                       = (j-1)*ni + i + kshift
                                    state%ordcoef(ji+mshift) = state%coef(ikeep)
                                endif
                            enddo
                        enddo
                        kshift = kshift + nisq
                    enddo
                    lshift = lshift + state%nll(n)
                enddo
                mshift = mshift + state%nll4(n)
                nshift = nshift + state%ncol(n)
                nprev  = n
            endif
        enddo
        if( state%iswtyp_lond .gt. 0 )then
            state%sfmid_lond   = ( state%sflo_lond + state%sfhi_lond ) * 0.5e0_wp
            state%dhsf_lond    = ( state%sfhi_lond - state%sflo_lond ) * 0.5e0_wp
            state%dhsfinv_lond = 1.e0_wp / state%dhsf_lond
        endif
    end subroutine build_polynomial_tables

    !> @brief Build many-body-expansion index bounds and equivalent permutations.
    !> @brief Transform six H--H distances to one H2--H2 Jacobi coordinate set.
    subroutine jcoord_h4bmkp( state, rij, q, derr )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: rij(6)
        real(wp), intent(inout) :: q(6)
        real(wp), intent(inout) :: derr(6,6)
        integer                 :: i, j, k
        real(wp)                :: cc, ccs, cm1x, cm1y, cm1z, cm2x, cm2y, cm2z, dd(9,6), der(6), r1z, r2x, r2y, r2z, rij2(6), riji(6), rx, ry, rz, x(3,4), xd(3,4,6)
        real(wp), parameter     :: zero   = 0.e0_wp
        real(wp), parameter     :: half   = 0.5e0_wp
        real(wp), parameter     :: one    = 1.e0_wp
        real(wp), parameter     :: epsilo = 1.e-30_wp
        do i=1, 6
            riji(i) = 1.e0_wp/rij(i)
            rij2(i) = rij(i)*rij(i)
        enddo
        do i=1, 4
            do j=1, 3
                x(j,i) = zero
            enddo
        enddo
        do i=1, 6
            do j=1, 4
                do k=1, 3
                    xd(k,j,i) = zero
                enddo
            enddo
        enddo
        x(3,2)    = rij(1)
        xd(3,2,1) = one
        cc        = -(rij2(4)-rij2(1)-rij2(2))*half *(riji(1)*riji(2))
        der(1)    = (riji(2))-cc   *riji(1)
        der(2)    = (riji(1))-cc   *riji(2)
        der(3)    = zero
        der(4)    = -rij(4)*(riji(1)*riji(2))
        der(5)    = zero
        der(6)    = zero
        ccs       = one-cc   *cc
        x(3,3)    = rij(2)*cc
        do i=1, 6
            xd(3,3,i) = rij(2)*der(i)
        enddo
        xd(3,3,2) = xd(3,3,2)+cc
        ccs       = sqrt(abs(ccs   ))
        cc        = rij(2)*cc   /(ccs   +epsilo)
        x(2,3)    = rij(2)*ccs
        do i=1, 6
            xd(2,3,i) = -cc   *der(i)
        enddo
        xd(2,3,2) = xd(2,3,2)+ccs
        x(3,4)    = (rij2(3)+rij2(1)-rij2(5))*half*riji(1)
        xd(3,4,1) = one-x(3,4)*riji(1)
        xd(3,4,3) = rij(3)*riji(1)
        xd(3,4,5) = -rij(5)*riji(1)
        der(1)    = one/(x(2,3)+epsilo)
        x(2,4)    = (rij2(3)+x(2,3)**2+x(3,3)**2-2e0*x(3,4)* x(3,3)-rij2(6))*half*der(1)
        do i=1, 6
            xd(2,4,i) = (x(2,3)*xd(2,3,i)+x(3,3)*(xd(3,3,i) -xd(3,4,i))-x(3,4)*xd(3,3,i)-x(2,4) *xd(2,3,i))*der(1)
        enddo
        xd(2,4,3) = xd(2,4,3)+rij(3)*der(1)
        xd(2,4,6) = xd(2,4,6)-rij(6)*der(1)
        x(1,4)    = (rij2(3)-x(2,4)**2-x(3,4)**2)
        cc        = one/(sqrt(abs(x(1,4)))+epsilo)
        do i=1, 6
            xd(1,4,i) = -(x(2,4)*xd(2,4,i)+x(3,4)*xd(3,4,i))*cc
        enddo
        xd(1,4,3) = xd(1,4,3)+rij(3)*cc
        x(1,4)    = one/(cc   +epsilo)
        rx        = half*(x(1,3)+x(1,4))
        ry        = half*(x(2,3)+x(2,4))
        rz        = half*(x(3,3)+x(3,4)-rij(1))
        do i=1, 6
            dd(1,i) = half*(xd(1,3,i)+xd(1,4,i))
            dd(2,i) = half*(xd(2,3,i)+xd(2,4,i))
            dd(3,i) = half*(xd(3,3,i)+xd(3,4,i))
        enddo
        dd(3,1) = dd(3,1)-half
        q(1)    = sqrt(abs(rx*rx + ry*ry + rz*rz))
        cc      = one/(q(1)+epsilo)
        do i=1, 6
            derr(1,i) = (dd(1,i)*rx   +dd(2,i)*ry +dd(3,i)*rz   )*cc
        enddo
        r1z  = rij(1)
        q(2) = rij(1)
        do i=1, 6
            derr(2,i) = zero
        enddo
        derr(2,1) = one
        r2x       = x(1,4)-x(1,3)
        r2y       = x(2,4)-x(2,3)
        r2z       = x(3,4)-x(3,3)
        do i=1, 6
            dd(4,i) = xd(1,4,i)-xd(1,3,i)
            dd(5,i) = xd(2,4,i)-xd(2,3,i)
            dd(6,i) = xd(3,4,i)-xd(3,3,i)
        enddo
        q(3) = rij(6)
        ccs  = one/(q(3)+epsilo)
        do i=1, 6
            derr(3,i) = (dd(4,i)*r2x   +dd(5,i)*r2y +dd(6,i)*r2z   )*ccs
        enddo
        q(4) = rz*cc
        do i=1, 6
            derr(4,i) = (dd(3,i)-rz   *derr(1,i)*cc   )*cc
        enddo
        ccs  = q(1)*q(3)
        ccs  = one/(ccs   +epsilo)
        q(5) = (r2x*rx + r2y*ry + r2z*rz)*ccs
        do i=1, 6
            derr(5,i) = (dd(4,i)*rx + dd(1,i)*r2x + dd(5,i)*ry +dd(2,i)*r2y+ dd(6,i)*rz  + dd(3,i)*r2z -(q(1)*derr(3, &
                & i) + derr(1,i)*q(3))*q(5)) *ccs
        enddo
        cm1x = -r1z   *ry
        cm1y = r1z   *rx
        cm1z = zero
        cm2x = r2y *rz - r2z *ry
        cm2y = r2z *rx - r2x *rz
        cm2z = r2x *ry - r2y *rx
        do i=1, 6
            dd(7,i) = r2y *dd(3,i) + dd(5,i)*rz - r2z *dd(2,i) - dd(6,i)*ry
            dd(8,i) = r2z *dd(1,i) + dd(6,i)*rx - r2x *dd(3,i) - dd(4,i)*rz
            dd(9,i) = r2x *dd(2,i) + dd(4,i)*ry - r2y *dd(1,i) - dd(5,i)*rx
            dd(4,i) = -r1z *dd(2,i)
            dd(5,i) = r1z *dd(1,i)
            dd(6,i) = zero
        enddo
        dd(4,1) = dd(4,1)-ry
        dd(5,1) = dd(5,1)+rx
        r2x     = sqrt(cm1x**2 + cm1y**2 + cm1z**2)
        r2y     = sqrt(cm2x**2 + cm2y**2 + cm2z**2)
        cc      = one/(r2x+epsilo)
        ccs     = one/(r2y+epsilo)
        do i=1, 6
            dd(1,i) = (cm1x*dd(4,i) + cm1y*dd(5,i) + cm1z*dd(6,i))*cc
            dd(2,i) = (cm2x*dd(7,i) + cm2y*dd(8,i) + cm2z*dd(9,i))*ccs
        enddo
        rz   = cc *ccs
        q(6) = ( cm1x*cm2x + cm1y*cm2y + cm1z*cm2z )*rz
        do i=1, 6
            derr(6,i) = (cm1x*dd(7,i) + dd(4,i)*cm2x +  cm1y*dd(8,i) + dd(5,i)*cm2y +  cm1z*dd(9,i) + dd(6,i)*cm2z - &
                & (r2x*dd(2,i)  + dd(1,i)*r2y )*q(6))*rz
        enddo
    end subroutine jcoord_h4bmkp

    !> @brief Build three H2--H2 arrangements and their Cartesian Jacobians.
    subroutine geomet_h4bmkp( state, r, cc, ideriv, ierr )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(inout) :: r(6)
        real(wp), intent(in)    :: cc(4,3)
        integer, intent(in)     :: ideriv
        integer, intent(inout)  :: ierr
        integer                 :: i, i1, i2, iden, ig, isin1, isin2, j, k, l1, l2
        real(wp)                :: com(6,3), cosphi, costh1, costh2, dpdrrv(3), dpdrvp(2,3), dqdrjac(6,6), drrv1, drrv2, drrv3, qjac(6), r1, r2, r3, rjac(6), rrv(3), rv(6,3), rvp(2,3)
        real(wp), parameter     :: epsilon  = 1.e-06_wp
        real(wp), parameter     :: zh       = 0.5e0_wp
        real(wp), parameter     :: rmin_err = 0.301e0_wp
        real(wp), parameter     :: tolsin   = 1.e-7_wp
        do k = 1, 6
            i1      = state%ipart1(k)
            i2      = state%ipart2(k)
            rv(k,1) = cc(i2,1) - cc(i1,1)
            rv(k,2) = cc(i2,2) - cc(i1,2)
            rv(k,3) = cc(i2,3) - cc(i1,3)
            r(k)    = sqrt( rv(k,1)**2 + rv(k,2)**2 + rv(k,3)**2 )
        enddo
        ierr = 0
        do k = 1, 6
            if( r(k) .lt. rmin_err )then
                ierr = 1
                return
            endif
        enddo
        do i = 1, 3
            com(1,i) = zh*rv(1,i)
            com(2,i) = zh*rv(2,i)
            com(3,i) = zh*rv(3,i)
            com(4,i) = rv(1,i) + zh*rv(4,i)
            com(5,i) = rv(1,i) + zh*rv(5,i)
            com(6,i) = rv(2,i) + zh*rv(6,i)
        enddo
        if( ideriv .gt. 0 )then
            do k = 1, 6
                i1                     = state%ipart1(k)
                i2                     = state%ipart2(k)
                state%dridcc(k,i2)     = rv(k,1) / r(k)
                state%dridcc(k,i1)     = - state%dridcc(k,i2)
                state%dridcc(k,i2+4)   = rv(k,2) / r(k)
                state%dridcc(k,i1+4)   = - state%dridcc(k,i2+4)
                state%dridcc(k,i2+8)   = rv(k,3) / r(k)
                state%dridcc(k,i1+8)   = - state%dridcc(k,i2+8)
                state%drvdcc(k,1,i1)   = -1.e0_wp
                state%drvdcc(k,1,i2)   = 1.e0_wp
                state%drvdcc(k,2,i1+4) = -1.e0_wp
                state%drvdcc(k,2,i2+4) = 1.e0_wp
                state%drvdcc(k,3,i1+8) = -1.e0_wp
                state%drvdcc(k,3,i2+8) = 1.e0_wp
            enddo
            do i = 1, 3
                do j = 1, 12
                    state%dcomdcc(1,i,j) = zh * state%drvdcc(1,i,j)
                    state%dcomdcc(2,i,j) = zh * state%drvdcc(2,i,j)
                    state%dcomdcc(3,i,j) = zh * state%drvdcc(3,i,j)
                    state%dcomdcc(4,i,j) = state%drvdcc(1,i,j) + zh * state%drvdcc(4,i,j)
                    state%dcomdcc(5,i,j) = state%drvdcc(1,i,j) + zh * state%drvdcc(5,i,j)
                    state%dcomdcc(6,i,j) = state%drvdcc(2,i,j) + zh * state%drvdcc(6,i,j)
                enddo
            enddo
        endif
        do ig = 1, 3
            l1 = ig
            l2 = 7 - ig
            do i = 1, 3
                rrv(i)   = com(l2,i) - com(l1,i)
                rvp(1,i) = rv(l1,i)
                rvp(2,i) = rv(l2,i)
            enddo
            r1              = r(l1)
            r2              = r(l2)
            r3              = sqrt( rrv(1)**2 + rrv(2)**2 + rrv(3)**2 )
            state%geo(ig,1) = r1
            state%geo(ig,2) = r2
            state%geo(ig,3) = r3
            if( r3 .le. epsilon )then
                state%geo(ig,4) = 0.e0_wp
                state%geo(ig,5) = 0.e0_wp
                state%geo(ig,6) = 0.e0_wp
                do i = 0, 5
                    state%geo_sin(ig,i) = 0.e0_wp
                enddo
                if( ideriv .gt. 0 )then
                    do i = 1, 6
                        do j = 1, 12
                            state%dgeodcc(ig,i,j) = 0.e0_wp
                        enddo
                        do j = 1, 6
                            state%dgeodri(ig,i,j) = 0.e0_wp
                        enddo
                    enddo
                endif
            else
                costh1 = ( rv(l1,1)*rrv(1) + rv(l1,2)*rrv(2) + rv(l1,3)*rrv(3) ) / ( r3 * r1 )
                costh2 = ( rv(l2,1)*rrv(1) + rv(l2,2)*rrv(2) + rv(l2,3)*rrv(3) ) / ( r3 * r2 )
                call cross_h4bmkp(state, rvp, rrv, cosphi, dpdrvp, dpdrrv, ideriv, iden)
                isin1 = 1
                if( abs(costh1) .lt. 0.7e0_wp )then
                    state%geo_sin(ig,0) = costh1**2
                    state%geo_sin(ig,2) = 1.e0_wp - state%geo_sin(ig,0)
                    state%geo_sin(ig,4) = sqrt( state%geo_sin(ig,2) )
                else
                    state%geo_sin(ig,4) = state%ptm(1) / ( r1 * r3 )
                    state%geo_sin(ig,2) = state%geo_sin(ig,4)**2
                    if( state%geo_sin(ig,4) .lt. tolsin ) isin1 = 0
                    state%geo_sin(ig,0) = 1.e0_wp - state%geo_sin(ig,2)
                    if( costh1 .ge. 0.e0_wp )then
                        costh1 = sqrt( state%geo_sin(ig,0) )
                    else
                        costh1 = - sqrt( state%geo_sin(ig,0) )
                    endif
                endif
                isin2 = 1
                if( abs(costh2) .lt. 0.7e0_wp )then
                    state%geo_sin(ig,1) = costh2**2
                    state%geo_sin(ig,3) = 1.e0_wp - state%geo_sin(ig,1)
                    state%geo_sin(ig,5) = sqrt( state%geo_sin(ig,3) )
                else
                    state%geo_sin(ig,5) = state%ptm(2) / ( r2 * r3 )
                    state%geo_sin(ig,3) = state%geo_sin(ig,5)**2
                    if( state%geo_sin(ig,5) .lt. tolsin ) isin2 = 0
                    state%geo_sin(ig,1) = 1.e0_wp - state%geo_sin(ig,5)**2
                    if( costh2 .ge. 0.e0_wp )then
                        costh2 = sqrt( state%geo_sin(ig,1) )
                    else
                        costh2 = - sqrt( state%geo_sin(ig,1) )
                    endif
                endif
                state%geo(ig,4) = costh1
                state%geo(ig,5) = costh2
                if( isin1 .eq. 0 .or. isin2 .eq. 0 .or. iden .eq. 0 )then
                    cosphi = 0.e0_wp
                    iden   = 0
                endif
                state%geo(ig,6) = cosphi
                if( ideriv .gt. 0 )then
                    do j = 1, 12
                        drrv1                 = state%dcomdcc(l2,1,j) - state%dcomdcc(l1,1,j)
                        drrv2                 = state%dcomdcc(l2,2,j) - state%dcomdcc(l1,2,j)
                        drrv3                 = state%dcomdcc(l2,3,j) - state%dcomdcc(l1,3,j)
                        state%dgeodcc(ig,1,j) = state%dridcc(l1,j)
                        state%dgeodcc(ig,2,j) = state%dridcc(l2,j)
                        state%dgeodcc(ig,3,j) = ( rrv(1) * drrv1 + rrv(2) * drrv2 + rrv(3) * drrv3 ) / r3
                        if( isin1 .eq. 0 )then
                            state%dgeodcc(ig,4,j) = 0.e0_wp
                        else
                            state%dgeodcc(ig,4,j) = ( state%drvdcc(l1,1,j) * rrv(1) + rv(l1,1) * drrv1 + &
                                & state%drvdcc(l1,2,j) * rrv(2) + rv(l1,2) * drrv2 + state%drvdcc(l1,3,j) * rrv(3) + &
                                & rv(l1,3) * drrv3 ) / ( r3 * r1 ) - costh1 * ( state%dgeodcc(ig,3,j) / r3 + &
                                & state%dridcc(l1,j) / r1 )
                        endif
                        if( isin2 .eq. 0 )then
                            state%dgeodcc(ig,5,j) = 0.e0_wp
                        else
                            state%dgeodcc(ig,5,j) = ( state%drvdcc(l2,1,j) * rrv(1) + rv(l2,1) * drrv1 + &
                                & state%drvdcc(l2,2,j) * rrv(2) + rv(l2,2) * drrv2 + state%drvdcc(l2,3,j) * rrv(3) + &
                                & rv(l2,3) * drrv3 ) / ( r3 * r2 ) - costh2 * ( state%dgeodcc(ig,3,j) / r3 + &
                                & state%dridcc(l2,j) / r2 )
                        endif
                        if( iden .eq. 0 )then
                            state%dgeodcc(ig,6,j) = 0.e0_wp
                        else
                            state%dgeodcc(ig,6,j) = dpdrvp(1,1) * state%drvdcc(l1,1,j) + dpdrvp(1,2) * &
                                & state%drvdcc(l1,2,j) + dpdrvp(1,3) * state%drvdcc(l1,3,j) + dpdrvp(2,1) * &
                                & state%drvdcc(l2,1,j) + dpdrvp(2,2) * state%drvdcc(l2,2,j) + dpdrvp(2,3) * &
                                & state%drvdcc(l2,3,j) + dpdrrv(1) * drrv1 + dpdrrv(2) * drrv2 + dpdrrv(3) * drrv3
                        endif
                    enddo
                    ! derivatives w.r.t distances r():
                    if( ideriv .gt. 1 )then
                        do j = 1, 6
                            rjac(j) = r( state%jco_ri(j,ig) )
                        enddo
                        call jcoord_h4bmkp(state, rjac, qjac, dqdrjac)
                        do i = 1, 6
                            do j = 1, 6
                                state%dgeodri( ig,state%jco_geo(i),state%jco_ri(j,ig) ) = dqdrjac(i,j)
                            enddo
                        enddo
                        if( iden .eq. 0 )then
                            do j = 1, 6
                                state%dgeodri(ig,6,j) = 0.e0_wp
                            enddo
                        endif
                    endif
                endif
            endif
        enddo
    end subroutine geomet_h4bmkp

    !> @brief Normalize arrangement weights and differentiate them when requested.
    subroutine vcorrec_h4bmkp( state, ideriv )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        integer, intent(in) :: ideriv
        integer             :: ig, j
        real(wp)            :: dsadcc(3,12), dsadgeo_1, dsadgeo_2, dsadgeo_3, dsadri(3,6), eta(3), etarg(3), sa(3), sasum, sasumtol
        real(wp), parameter :: epsilon = 1.e-06_wp
        real(wp), parameter :: tolfsa  = 1.e-18_wp
        do ig = 1, 3
            if( state%geo(ig,3) .le. epsilon )then
                sa(ig) = 0.e0_wp
            else
                etarg(ig) = 1.e0_wp / state%geo(ig,1)**2 + 1.e0_wp / state%geo(ig,2)**2
                eta(ig)   = 1.e0_wp / ( state%geo(ig,3) * sqrt( etarg(ig) ) )
                sa(ig)    = exp( -state%beta_j * eta(ig) )
            endif
        enddo
        sasumtol = ( sa(1) + sa(2) + sa(3) ) * tolfsa
        do ig = 1, 3
            if( sa(ig) .lt. sasumtol ) sa(ig) = 0.e0_wp
        enddo
        sasum = sa(1) + sa(2) + sa(3)
        do ig = 1, 3
            state%fsa(ig) = sa(ig) / sasum
        enddo
        if( ideriv .gt. 0 )then
            do ig = 1, 3
                if( sa(ig) .eq. 0.e0_wp )then
                    do j = 1, 12
                        state%dfsadcc(ig,j) = 0.e0_wp
                        dsadcc(ig,j)        = 0.e0_wp
                    enddo
                    do j = 1, 6
                        state%dfsadri(ig,j) = 0.e0_wp
                        dsadri(ig,j)        = 0.e0_wp
                    enddo
                else
                    dsadgeo_1 = - state%beta_j * sa(ig) * eta(ig) / ( etarg(ig) * state%geo(ig,1)**3 )
                    dsadgeo_2 = - state%beta_j * sa(ig) * eta(ig) / ( etarg(ig) * state%geo(ig,2)**3 )
                    dsadgeo_3 = state%beta_j * sa(ig) * eta(ig) / state%geo(ig,3)
                    do j = 1, 12
                        dsadcc(ig,j) = dsadgeo_1 * state%dgeodcc(ig,1,j) + dsadgeo_2 * state%dgeodcc(ig,2,j) + &
                            & dsadgeo_3 * state%dgeodcc(ig,3,j)
                    enddo
                    if( ideriv .gt. 1 )then
                        do j = 1, 6
                            dsadri(ig,j) = dsadgeo_1 * state%dgeodri(ig,1,j) + dsadgeo_2 * state%dgeodri(ig,2,j) + &
                                & dsadgeo_3 * state%dgeodri(ig,3,j)
                        enddo
                    endif
                endif
            enddo
            do ig = 1, 3
                if( sa(ig) .gt. 0.e0_wp )then
                    do j = 1, 12
                        state%dfsadcc(ig,j) = ( dsadcc(ig,j) - state%fsa(ig) * ( dsadcc(1,j) + dsadcc(2,j) + dsadcc(3, &
                            & j) ) ) / sasum
                    enddo
                    if( ideriv .gt. 1 )then
                        do j = 1, 6
                            state%dfsadri(ig,j) = ( dsadri(ig,j) - state%fsa(ig) * ( dsadri(1,j) + dsadri(2,j) + &
                                & dsadri(3,j) ) ) / sasum
                        enddo
                    endif
                endif
            enddo
        endif
    end subroutine vcorrec_h4bmkp

    !> @brief Compute the dihedral cosine and derivatives from Jacobi vectors.
    subroutine cross_h4bmkp( state, rvp, rrv, cosphi, dpdrvp, dpdrrv, ideriv, iden )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: rvp(2,3)
        real(wp), intent(in)    :: rrv(3)
        real(wp), intent(inout) :: cosphi
        real(wp), intent(inout) :: dpdrvp(2,3)
        real(wp), intent(inout) :: dpdrrv(3)
        integer, intent(in)     :: ideriv
        integer, intent(inout)  :: iden
        real(wp)                :: denom, ptm1sq, ptm2sq
        real(wp), parameter     :: toldenom = 1.e-15_wp
        state%pt(1,1) = rvp(1,2)*rrv(3) - rrv(2)*rvp(1,3)
        state%pt(1,2) = rrv(1)*rvp(1,3) - rvp(1,1)*rrv(3)
        state%pt(1,3) = rvp(1,1)*rrv(2) - rrv(1)*rvp(1,2)
        ptm1sq        = state%pt(1,1)**2 + state%pt(1,2)**2 + state%pt(1,3)**2
        state%ptm(1)  = sqrt( ptm1sq )
        state%pt(2,1) = rvp(2,2)*rrv(3) - rrv(2)*rvp(2,3)
        state%pt(2,2) = rrv(1)*rvp(2,3) - rvp(2,1)*rrv(3)
        state%pt(2,3) = rvp(2,1)*rrv(2) - rrv(1)*rvp(2,2)
        ptm2sq        = state%pt(2,1)**2 + state%pt(2,2)**2 + state%pt(2,3)**2
        state%ptm(2)  = sqrt( ptm2sq )
        denom         = state%ptm(1) * state%ptm(2)
        if( min( state%ptm(1), state%ptm(2), denom ) .gt. toldenom )then
            cosphi = ( state%pt(1,1)*state%pt(2,1) + state%pt(1,2)*state%pt(2,2) + state%pt(1,3)*state%pt(2,3) ) / denom
            iden   = 1
            if( ideriv .gt. 0 )then
                dpdrvp(1,1) = ( rrv(2) * state%pt(2,3) - rrv(3) * state%pt(2,2) ) / denom - cosphi * ( rrv(2) * &
                    & state%pt(1,3) - rrv(3) * state%pt(1,2) ) / ptm1sq
                dpdrvp(1,2) = ( rrv(3) * state%pt(2,1) - rrv(1) * state%pt(2,3) ) / denom - cosphi * ( rrv(3) * &
                    & state%pt(1,1) - rrv(1) * state%pt(1,3) ) / ptm1sq
                dpdrvp(1,3) = ( rrv(1) * state%pt(2,2) - rrv(2) * state%pt(2,1) ) / denom - cosphi * ( rrv(1) * &
                    & state%pt(1,2) - rrv(2) * state%pt(1,1) ) / ptm1sq
                dpdrvp(2,1) = ( rrv(2) * state%pt(1,3) - rrv(3) * state%pt(1,2) ) / denom - cosphi * ( rrv(2) * &
                    & state%pt(2,3) - rrv(3) * state%pt(2,2) ) / ptm2sq
                dpdrvp(2,2) = ( rrv(3) * state%pt(1,1) - rrv(1) * state%pt(1,3) ) / denom - cosphi * ( rrv(3) * &
                    & state%pt(2,1) - rrv(1) * state%pt(2,3) ) / ptm2sq
                dpdrvp(2,3) = ( rrv(1) * state%pt(1,2) - rrv(2) * state%pt(1,1) ) / denom - cosphi * ( rrv(1) * &
                    & state%pt(2,2) - rrv(2) * state%pt(2,1) ) / ptm2sq
                dpdrrv(1) = ( rvp(1,3) * state%pt(2,2) - rvp(1,2) * state%pt(2,3) + rvp(2,3) * state%pt(1,2) - rvp(2, &
                    & 2) * state%pt(1,3) ) / denom - cosphi * ( ( rvp(1,3) * state%pt(1,2) - rvp(1,2) * state%pt(1, &
                    & 3) ) / ptm1sq + ( rvp(2,3) * state%pt(2,2) - rvp(2,2) * state%pt(2,3) ) / ptm2sq )
                dpdrrv(2) = ( rvp(1,1) * state%pt(2,3) - rvp(1,3) * state%pt(2,1) + rvp(2,1) * state%pt(1,3) - rvp(2, &
                    & 3) * state%pt(1,1) ) / denom - cosphi * ( ( rvp(1,1) * state%pt(1,3) - rvp(1,3) * state%pt(1, &
                    & 1) ) / ptm1sq + ( rvp(2,1) * state%pt(2,3) - rvp(2,3) * state%pt(2,1) ) / ptm2sq )
                dpdrrv(3) = ( rvp(1,2) * state%pt(2,1) - rvp(1,1) * state%pt(2,2) + rvp(2,2) * state%pt(1,1) - rvp(2, &
                    & 1) * state%pt(1,2) ) / denom - cosphi * ( ( rvp(1,2) * state%pt(1,1) - rvp(1,1) * state%pt(1, &
                    & 2) ) / ptm1sq + ( rvp(2,2) * state%pt(2,1) - rvp(2,1) * state%pt(2,2) ) / ptm2sq )
            endif
        else
            cosphi = 0.e0_wp
            iden   = 0
        endif
    end subroutine cross_h4bmkp

    !> @brief Sum the four embedded BKMP2 H3 nonseparable contributions.
    subroutine vnsb_h4bmkp( state, r, ideriv, vnsb, dvnsbdri )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r(6)
        integer, intent(in)     :: ideriv
        real(wp), intent(inout) :: vnsb
        real(wp), intent(inout) :: dvnsbdri(6)
        real(wp)                :: dvabc(6), dvabd(6), dvacd(6), dvbcd(6), r1, r2, r3, r4, r5, r6, vabc, vabd, vacd, vbcd
        r1 = r(1)
        r2 = r(2)
        r3 = r(3)
        r4 = r(4)
        r5 = r(5)
        r6 = r(6)
        call vxyzcal_h4bmkp(state, r1, r4, r2, r3, r5, r6, ideriv, vabc, dvabc)
        call vxyzcal_h4bmkp(state, r1, r5, r3, r2, r4, r6, ideriv, vabd, dvabd)
        call vxyzcal_h4bmkp(state, r2, r6, r3, r1, r4, r5, ideriv, vacd, dvacd)
        call vxyzcal_h4bmkp(state, r4, r6, r5, r1, r2, r3, ideriv, vbcd, dvbcd)
        vnsb = vabc + vabd + vacd + vbcd
        if( ideriv .gt. 0 )then
            dvnsbdri(1) = dvabc(1) + dvabd(1) + dvacd(4) + dvbcd(4)
            dvnsbdri(2) = dvabc(3) + dvabd(4) + dvacd(1) + dvbcd(5)
            dvnsbdri(3) = dvabc(4) + dvabd(3) + dvacd(3) + dvbcd(6)
            dvnsbdri(4) = dvabc(2) + dvabd(5) + dvacd(5) + dvbcd(1)
            dvnsbdri(5) = dvabc(5) + dvabd(2) + dvacd(6) + dvbcd(3)
            dvnsbdri(6) = dvabc(6) + dvabd(6) + dvacd(2) + dvbcd(2)
        endif
    end subroutine vnsb_h4bmkp

    !> @brief Evaluate geometry-selector powers and their distance derivatives.
    subroutine vpcalc_h4bmkp( state, r, ideriv )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in) :: r(6)
        integer, intent(in)  :: ideriv
        integer              :: j
        real(wp)             :: dvpdri, vp(6)
        do j = 1, 6
            vp(j) = exp( -state%ap_vpcalc_2 * r(j)**state%ap_vpcalc_3 )
        enddo
        state%sumvp(1) = vp(2) + vp(5) + vp(3) + vp(4)
        state%sumvp(2) = vp(1) + vp(6) + vp(3) + vp(4)
        state%sumvp(3) = vp(1) + vp(6) + vp(2) + vp(5)
        if( ideriv .gt. 0 )then
            do j = 1, 6
                dvpdri                                             = -state%ap_vpcalc_2 * state%ap_vpcalc_3 * r(j)** &
                    & (state%ap_vpcalc_3-1.e0_wp) * vp(j)
                if( j .ne. 1 .and. j .ne. 6 ) state%dsumvpdri(1,j) = dvpdri
                if( j .ne. 2 .and. j .ne. 5 ) state%dsumvpdri(2,j) = dvpdri
                if( j .ne. 3 .and. j .ne. 4 ) state%dsumvpdri(3,j) = dvpdri
            enddo
        endif
    end subroutine vpcalc_h4bmkp

    !> @brief Evaluate the composite London H2--H2 interaction contribution.
    subroutine comlon_h4bmkp( state, r, ideriv, vlon, dvlondri, dvlondfsa )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r(6)
        integer, intent(in)     :: ideriv
        real(wp), intent(inout) :: vlon
        real(wp), intent(inout) :: dvlondri(6)
        real(wp), intent(inout) :: dvlondfsa(3)
        integer                 :: i
        real(wp)                :: deps2dfsa(3), deps2dri(6), derfac, djdri(6), dqdri(6), dsumjhalfdri(6), eps, eps2, esing(3), etrip(3), o_msw, o_mxd, o_pxd, rpd(6), s_f, sfterm1, sfterm2, sfterm3, sumjhalf, sumq, sw_use, termj1, termj2, termj3, termjsub, vj(6), vq(6)
        real(wp), parameter     :: f_1o4      = 1.e0_wp / 4.e0_wp
        real(wp), parameter     :: eps_h3lond = 1.e-6_wp
        real(wp), parameter     :: half       = 0.5e0_wp
        do i = 1, 6
            rpd(i)      = r(i) + state%rdelta_lond
            deps2dri(i) = 0.e0_wp
        enddo
        deps2dfsa(1) = 0.e0_wp
        deps2dfsa(2) = 0.e0_wp
        deps2dfsa(3) = 0.e0_wp
        sfterm1      = rpd(2) * rpd(3) * rpd(4) * rpd(5)
        sfterm2      = rpd(1) * rpd(3) * rpd(4) * rpd(6)
        sfterm3      = rpd(1) * rpd(2) * rpd(5) * rpd(6)
        if( state%iswtyp_lond .eq. 0 )then
            eps = eps_h3lond
        else
            s_f = ( state%fsa(1) * sfterm1 + state%fsa(2) * sfterm2 + state%fsa(3) * sfterm3 )**f_1o4
            if( s_f .le. state%sflo_lond )then
                eps = state%eps_lond
            elseif( s_f .ge. state%sfhi_lond )then
                eps = eps_h3lond
            else
                o_mxd = 1.e0_wp - ( s_f - state%sfmid_lond ) * state%dhsfinv_lond
                o_pxd = 1.e0_wp + ( s_f - state%sfmid_lond ) * state%dhsfinv_lond
                if( s_f .le. state%sfmid_lond )then
                    sw_use = o_pxd**3 * o_mxd * 0.5e0_wp
                    o_msw  = 1.e0_wp - sw_use
                    if( ideriv .gt. 0 ) derfac = 0.5e0_wp * state%dhsfinv_lond * o_pxd**2 * ( 3.e0_wp * o_mxd - o_pxd ) * ( eps_h3lond - state%eps_lond )
                else
                    o_msw  = o_mxd**3 * o_pxd * 0.5e0_wp
                    sw_use = 1.e0_wp - o_msw
                    if( ideriv .gt. 0 ) derfac = 0.5e0_wp * state%dhsfinv_lond * o_mxd**2 * ( o_mxd - 3.e0_wp * o_pxd ) * ( state%eps_lond - eps_h3lond )
                endif
                eps = o_msw * state%eps_lond + sw_use * eps_h3lond
                if( ideriv .gt. 0 )then
                    derfac       = derfac * 2.e0_wp * eps * f_1o4 / s_f**3
                    deps2dri(1)  = ( state%fsa(2) * rpd(3) * rpd(4) + state%fsa(3) * rpd(2) * rpd(5) ) * rpd(6) * derfac
                    deps2dri(2)  = ( state%fsa(1) * rpd(3) * rpd(4) + state%fsa(3) * rpd(1) * rpd(6) ) * rpd(5) * derfac
                    deps2dri(3)  = ( state%fsa(1) * rpd(2) * rpd(5) + state%fsa(2) * rpd(1) * rpd(6) ) * rpd(4) * derfac
                    deps2dri(4)  = ( state%fsa(1) * rpd(2) * rpd(5) + state%fsa(2) * rpd(1) * rpd(6) ) * rpd(3) * derfac
                    deps2dri(5)  = ( state%fsa(1) * rpd(3) * rpd(4) + state%fsa(3) * rpd(1) * rpd(6) ) * rpd(2) * derfac
                    deps2dri(6)  = ( state%fsa(2) * rpd(3) * rpd(4) + state%fsa(3) * rpd(2) * rpd(5) ) * rpd(1) * derfac
                    deps2dfsa(1) = sfterm1 * derfac
                    deps2dfsa(2) = sfterm2 * derfac
                    deps2dfsa(3) = sfterm3 * derfac
                endif
            endif
        endif
        eps2 = eps**2
        do i = 1, 6
            call vh2opt95_h4bmkp(state, r(i), esing, ideriv)
            call triplet95_h4bmkp(state, r(i), etrip, ideriv)
            vq(i) = half * ( esing(1) + etrip(1) )
            vj(i) = half * ( esing(1) - etrip(1) )
            if( ideriv .gt. 0 )then
                dqdri(i) = half * ( esing(2) + etrip(2) )
                djdri(i) = half * ( esing(2) - etrip(2) )
            endif
        enddo
        sumq     = vq(1) + vq(2) + vq(3) + vq(4) + vq(5) + vq(6)
        termj1   = vj(1) + vj(6) - vj(2) - vj(5)
        termj2   = vj(2) + vj(5) - vj(3) - vj(4)
        termj3   = vj(3) + vj(4) - vj(1) - vj(6)
        sumjhalf = half * ( termj1**2 + termj2**2 + termj3**2 )
        termjsub = sqrt( eps2 + sumjhalf )
        vlon     = sumq - termjsub
        if( ideriv .gt. 0 )then
            dsumjhalfdri(1) = ( termj1 - termj3 ) * djdri(1)
            dsumjhalfdri(2) = ( termj2 - termj1 ) * djdri(2)
            dsumjhalfdri(3) = ( termj3 - termj2 ) * djdri(3)
            dsumjhalfdri(4) = ( termj3 - termj2 ) * djdri(4)
            dsumjhalfdri(5) = ( termj2 - termj1 ) * djdri(5)
            dsumjhalfdri(6) = ( termj1 - termj3 ) * djdri(6)
            derfac          = -0.5e0_wp / termjsub
            do i = 1, 6
                dvlondri(i) = dqdri(i) + ( deps2dri(i) + dsumjhalfdri(i) ) * derfac
            enddo
            dvlondfsa(1) = deps2dfsa(1) * derfac
            dvlondfsa(2) = deps2dfsa(2) * derfac
            dvlondfsa(3) = deps2dfsa(3) * derfac
        endif
    end subroutine comlon_h4bmkp

    !> @brief Form angular basis functions and their geometry derivatives.
    subroutine getyy_h4bmkp( state, ideriv )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        integer, intent(in) :: ideriv
        integer             :: ig
        real(wp)            :: c2phi, cphi, ct1, ct1sq, ct2, ct2sq, st1, st2, y00, y20_1, y20_2, y21_1, y21_2, y22_1, y22_2
        real(wp), parameter :: zh     = 0.5e0_wp
        real(wp), parameter :: z3h    = 1.5e0_wp
        real(wp), parameter :: twoz3h = 3.e0_wp
        real(wp), parameter :: tol    = 1.e-7_wp
        do ig = 1, 3
            if( state%fsa(ig) .gt. 0.e0_wp )then
                ct1            = state%geo(ig,4)
                ct2            = state%geo(ig,5)
                ct1sq          = state%geo_sin(ig,0)
                ct2sq          = state%geo_sin(ig,1)
                cphi           = state%geo(ig,6)
                st1            = state%geo_sin(ig,4)
                st2            = state%geo_sin(ig,5)
                c2phi          = 2.e0_wp * cphi**2 - 1.e0_wp
                y00            = state%rt_zfl0m0
                y20_1          = state%rt_zfm0 * ( z3h * ct1sq - zh )
                y20_2          = state%rt_zfm0 * ( z3h * ct2sq - zh )
                y21_1          = state%rt_zfm1_3 * ct1 * st1
                y21_2          = state%rt_zfm1_3 * ct2 * st2
                y22_1          = state%rt_zfm2_3 * state%geo_sin(ig,2)
                y22_2          = state%rt_zfm2_3 * state%geo_sin(ig,3)
                state%yy(ig,1) = state%z4pi * y00**2
                state%yy(ig,2) = state%z4pi * y20_1 * y20_2
                state%yy(ig,3) = state%z4pi_rt2 * y21_1 * y21_2 * cphi
                state%yy(ig,4) = state%z4pi_rt2 * y22_1 * y22_2 * c2phi
                state%yy(ig,5) = state%z4pi * y20_1 * y00
                state%yy(ig,6) = state%z4pi * y00 * y20_2
                if( ideriv .gt. 0 )then
                    state%dyydgeo(ig,2,4) = state%z4pi * state%rt_zfm0 * twoz3h * ct1 * y20_2
                    state%dyydgeo(ig,2,5) = state%z4pi * y20_1 * state%rt_zfm0 * twoz3h * ct2
                    if( st1 .gt. tol )then
                        state%dyydgeo(ig,3,4) = state%z4pi_rt2 * state%rt_zfm1_3 * ( st1 - ct1sq / st1 ) * y21_2 * cphi
                    else
                        state%dyydgeo(ig,3,4) = 0.e0_wp
                    endif
                    if( st2 .gt. tol )then
                        state%dyydgeo(ig,3,5) = state%z4pi_rt2 * y21_1 * state%rt_zfm1_3 * ( st2 - ct2sq / st2 ) * cphi
                    else
                        state%dyydgeo(ig,3,5) = 0.e0_wp
                    endif
                    state%dyydgeo(ig,3,6) = state%z4pi_rt2 * y21_1 * y21_2
                    state%dyydgeo(ig,4,4) = -2.e0_wp * state%z4pi_rt2 * state%rt_zfm2_3 * ct1 * y22_2 * c2phi
                    state%dyydgeo(ig,4,5) = -2.e0_wp * state%z4pi_rt2 * y22_1 * state%rt_zfm2_3 * ct2 * c2phi
                    state%dyydgeo(ig,4,6) = state%z4pi_rt2 * y22_1 * y22_2 * 4.e0_wp * cphi
                    state%dyydgeo(ig,5,4) = state%z4pi * state%rt_zfm0 * twoz3h * ct1 * y00
                    state%dyydgeo(ig,6,5) = state%z4pi * y00 * state%rt_zfm0 * twoz3h * ct2
                endif
            endif
        enddo
    end subroutine getyy_h4bmkp

    !> @brief Evaluate linear-geometry and many-body polynomial corrections.
    subroutine vlinear_h4bmkp( state, r, ideriv, vlinear, vlinmbe, dvlindri, dvlindgeo, dvmbedri )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r(6)
        integer, intent(in)     :: ideriv
        real(wp), intent(inout) :: vlinear(3)
        real(wp), intent(inout) :: vlinmbe
        real(wp), intent(inout) :: dvlindri(3,6)
        real(wp), intent(inout) :: dvlindgeo(3,6)
        real(wp), intent(inout) :: dvmbedri(6)
        integer                 :: i, ig, ii, iklo, j, jj, jm1, jp1, k, kk, kt, nd0, nd1, nd2, nd3, nd4, ndone
        real(wp)                :: dadri(8,6), dds, detdri(8,6), et(8), r_3, r_a(-9:9), r_b(-9:9), r_r(-9:9), rbar, rho, tcda(64), tcdb(64), tcdr(64), tcor(64), teda(64), tedb(64), tedr(64), temp(64), v_s, vmbe, vs(6), vsa(6), vsb(6), vsr(6)
        integer, parameter      :: nterms  = 8
        real(wp), parameter     :: r3toler = 1.e-6_wp
        vmbe             = 0.e0_wp
        state%vlinear(1) = 0.e0_wp
        state%vlinear(2) = 0.e0_wp
        state%vlinear(3) = 0.e0_wp
        call asqsub_h4bmkp(state, r, ideriv, rbar, rho, state%aselect, dadri)
        et(1) = exp( -state%beta(1) * rho**3 )
        if( ideriv .gt. 0 )then
            dds = -3.e0_wp * state%beta(1) * rho * et(1)
            do i = 1, 6
                detdri(1,i)    = dds * r(i)
                dvlindri(1,i)  = 0.e0_wp
                dvlindri(2,i)  = 0.e0_wp
                dvlindri(3,i)  = 0.e0_wp
                dvlindgeo(1,i) = 0.e0_wp
                dvlindgeo(2,i) = 0.e0_wp
                dvlindgeo(3,i) = 0.e0_wp
            enddo
        endif
        do kt = 2, 8
            if( state%beta(kt) .le. 0.e0_wp )then
                et(kt) = 1.e0_wp
                do i = 1, 6
                    detdri(kt,i) = 0.e0_wp
                enddo
            elseif( state%beta(kt) .eq. state%beta(1) )then
                et(kt) = et(1)
                do i = 1, 6
                    detdri(kt,i) = detdri(1,i)
                enddo
            else
                et(kt) = exp( -state%beta(kt) * rho**3 )
                if( ideriv .gt. 0 )then
                    dds = -3.e0_wp * state%beta(kt) * rho * et(kt)
                    do i = 1, 6
                        detdri(kt,i) = dds * r(i)
                    enddo
                endif
            endif
        enddo
        if( state%ikmin_lo .eq. 0 .and. ideriv .gt. 0 )then
            iklo    = 0
            r_a(-1) = 1.e0_wp
            r_b(-1) = 1.e0_wp
            r_b(-1) = 1.e0_wp
        else
            iklo = state%ikmin_lo - min( ideriv, 1 )
        endif
        r_a(0) = 1.e0_wp
        r_b(0) = 1.e0_wp
        r_r(0) = 1.e0_wp
        do ig = 1, 3
            if( state%fsa(ig) .gt. 0.e0_wp )then
                r_a(1) = state%geo(ig,1)
                r_b(1) = state%geo(ig,2)
                r_r(1) = state%geo(ig,3)
                jm1    = 1
                do j = 2, state%ikmax_hi
                    r_a(j) = r_a(jm1) * r_a(1)
                    r_b(j) = r_b(jm1) * r_b(1)
                    r_r(j) = r_r(jm1) * r_r(1)
                    jm1    = j
                enddo
                if( iklo .lt. 0 )then
                    jp1 = 0
                    r_3 = max( r_r(1), r3toler )
                    do j = -1, iklo, -1
                        r_a(j) = r_a(jp1) / r_a(1)
                        r_b(j) = r_b(jp1) / r_b(1)
                        r_r(j) = r_r(jp1) / r_3
                        jp1    = j
                    enddo
                endif
                nd0   = 0
                ndone = 0
                if( state%iv(1) .gt. 0 )then
                    j = 0
                    do kk = state%kmin(1), state%kmax(1)
                        do ii = state%imin(1), state%imax(1)
                            do jj = state%imin(1), state%imax(1)
                                j       = j + 1
                                tcor(j) = r_a(ii) * r_b(jj) * r_r(kk)
                                temp(j) = r_a(jj) * r_b(ii) * r_r(kk)
                                if( ideriv .gt. 0 )then
                                    tcda(j) = ii * r_a(ii-1) * r_b(jj) * r_r(kk)
                                    tcdb(j) = jj * r_a(ii) * r_b(jj-1) * r_r(kk)
                                    tcdr(j) = kk * r_a(ii) * r_b(jj) * r_r(kk-1)
                                    teda(j) = jj * r_a(jj-1) * r_b(ii) * r_r(kk)
                                    tedb(j) = ii * r_a(jj) * r_b(ii-1) * r_r(kk)
                                    tedr(j) = kk * r_a(jj) * r_b(ii) * r_r(kk-1)
                                endif
                            enddo
                        enddo
                    enddo
                    nd1 = nd0 + state%nll(1)
                    nd2 = nd1 + state%nll(1)
                    nd3 = nd2 + state%nll(1)
                    nd4 = state%nkeep(1)
                    do i = 1, 6
                        vs(i) = 0.e0_wp
                    enddo
                    do j = 1, state%nll(1)
                        vs(1) = vs(1) + tcor(j) * state%ordcoef(j)
                        vs(2) = vs(2) + tcor(j) * state%ordcoef(j+nd1)
                        vs(3) = vs(3) + tcor(j) * state%ordcoef(j+nd2)
                        vs(4) = vs(4) + tcor(j) * state%ordcoef(j+nd3)
                        vs(5) = vs(5) + tcor(j) * state%coef(j+nd4)
                        vs(6) = vs(6) + temp(j) * state%coef(j+nd4)
                    enddo
                    v_s = vs(1) * state%yy(ig,1) + vs(2) * state%yy(ig,2) + vs(3) * state%yy(ig,3) + vs(4) * &
                        & state%yy(ig,4) + vs(5) * state%yy(ig,5) + vs(6) * state%yy(ig,6)
                    state%vlinear(ig) = state%sumvp(ig) * v_s * et(1)
                    if( ideriv .gt. 0 )then
                        do i = 1, 6
                            vsa(i) = 0.e0_wp
                            vsb(i) = 0.e0_wp
                            vsr(i) = 0.e0_wp
                        enddo
                        do j = 1, state%nll(1)
                            vsa(1) = vsa(1) + tcda(j) * state%ordcoef(j)
                            vsa(2) = vsa(2) + tcda(j) * state%ordcoef(j+nd1)
                            vsa(3) = vsa(3) + tcda(j) * state%ordcoef(j+nd2)
                            vsa(4) = vsa(4) + tcda(j) * state%ordcoef(j+nd3)
                            vsa(5) = vsa(5) + tcda(j) * state%coef(j+nd4)
                            vsa(6) = vsa(6) + teda(j) * state%coef(j+nd4)
                            vsb(1) = vsb(1) + tcdb(j) * state%ordcoef(j)
                            vsb(2) = vsb(2) + tcdb(j) * state%ordcoef(j+nd1)
                            vsb(3) = vsb(3) + tcdb(j) * state%ordcoef(j+nd2)
                            vsb(4) = vsb(4) + tcdb(j) * state%ordcoef(j+nd3)
                            vsb(5) = vsb(5) + tcdb(j) * state%coef(j+nd4)
                            vsb(6) = vsb(6) + tedb(j) * state%coef(j+nd4)
                            vsr(1) = vsr(1) + tcdr(j) * state%ordcoef(j)
                            vsr(2) = vsr(2) + tcdr(j) * state%ordcoef(j+nd1)
                            vsr(3) = vsr(3) + tcdr(j) * state%ordcoef(j+nd2)
                            vsr(4) = vsr(4) + tcdr(j) * state%ordcoef(j+nd3)
                            vsr(5) = vsr(5) + tcdr(j) * state%coef(j+nd4)
                            vsr(6) = vsr(6) + tedr(j) * state%coef(j+nd4)
                        enddo
                        do j = 1, 6
                            dvlindri(ig,j) = ( state%sumvp(ig) * detdri(1,j) + state%dsumvpdri(ig,j) * et(1) ) * v_s
                        enddo
                        dvlindgeo(ig,1) = state%sumvp(ig) * et(1) * ( vsa(1) * state%yy(ig,1) + vsa(2) * state%yy(ig, &
                            & 2) + vsa(3) * state%yy(ig,3) + vsa(4) * state%yy(ig,4) + vsa(5) * state%yy(ig,5) + &
                            & vsa(6) * state%yy(ig,6) )
                        dvlindgeo(ig,2) = state%sumvp(ig) * et(1) * ( vsb(1) * state%yy(ig,1) + vsb(2) * state%yy(ig, &
                            & 2) + vsb(3) * state%yy(ig,3) + vsb(4) * state%yy(ig,4) + vsb(5) * state%yy(ig,5) + &
                            & vsb(6) * state%yy(ig,6) )
                        dvlindgeo(ig,3) = state%sumvp(ig) * et(1) * ( vsr(1) * state%yy(ig,1) + vsr(2) * state%yy(ig, &
                            & 2) + vsr(3) * state%yy(ig,3) + vsr(4) * state%yy(ig,4) + vsr(5) * state%yy(ig,5) + &
                            & vsr(6) * state%yy(ig,6) )
                        dvlindgeo(ig,4) = ( vs(2) * state%dyydgeo(ig,2,4) + vs(3) * state%dyydgeo(ig,3,4) + vs(4) * &
                            & state%dyydgeo(ig,4,4) + vs(5) * state%dyydgeo(ig,5,4) ) * state%sumvp(ig) * et(1)
                        dvlindgeo(ig,5) = ( vs(2) * state%dyydgeo(ig,2,5) + vs(3) * state%dyydgeo(ig,3,5) + vs(4) * &
                            & state%dyydgeo(ig,4,5) + vs(6) * state%dyydgeo(ig,6,5) ) * state%sumvp(ig) * et(1)
                        dvlindgeo(ig,6) = ( vs(3) * state%dyydgeo(ig,3,6) + vs(4) * state%dyydgeo(ig,4,6) ) * &
                            & state%sumvp(ig) * et(1)
                    endif
                    ndone = state%ncol(1)
                    nd0   = state%nll4(1)
                endif
                do k = 2, state%njac
                    if( state%iv(k) .gt. 0 )then
                        if( state%ieqp(k) .eq. 0 )then
                            j = 0
                            do kk = state%kmin(k), state%kmax(k)
                                do ii = state%imin(k), state%imax(k)
                                    do jj = state%imin(k), state%imax(k)
                                        j       = j + 1
                                        tcor(j) = r_a(ii) * r_b(jj) * r_r(kk)
                                        temp(j) = r_a(jj) * r_b(ii) * r_r(kk)
                                        if( ideriv .gt. 0 )then
                                            tcda(j) = ii * r_a(ii-1) * r_b(jj) * r_r(kk)
                                            tcdb(j) = jj * r_a(ii) * r_b(jj-1) * r_r(kk)
                                            tcdr(j) = kk * r_a(ii) * r_b(jj) * r_r(kk-1)
                                            teda(j) = jj * r_a(jj-1) * r_b(ii) * r_r(kk)
                                            tedb(j) = ii * r_a(jj) * r_b(ii-1) * r_r(kk)
                                            tedr(j) = kk * r_a(jj) * r_b(ii) * r_r(kk-1)
                                        endif
                                    enddo
                                enddo
                            enddo
                        endif
                        nd1 = nd0 + state%nll(k)
                        nd2 = nd1 + state%nll(k)
                        nd3 = nd2 + state%nll(k)
                        nd4 = ndone + state%nkeep(k)
                        do i = 1, 6
                            vs(i) = 0.e0_wp
                        enddo
                        do j = 1, state%nll(k)
                            vs(1) = vs(1) + tcor(j) * state%ordcoef(j+nd0)
                            vs(2) = vs(2) + tcor(j) * state%ordcoef(j+nd1)
                            vs(3) = vs(3) + tcor(j) * state%ordcoef(j+nd2)
                            vs(4) = vs(4) + tcor(j) * state%ordcoef(j+nd3)
                            vs(5) = vs(5) + tcor(j) * state%coef(j+nd4)
                            vs(6) = vs(6) + temp(j) * state%coef(j+nd4)
                        enddo
                        v_s = vs(1) * state%yy(ig,1) + vs(2) * state%yy(ig,2) + vs(3) * state%yy(ig,3) + vs(4) * &
                            & state%yy(ig,4) + vs(5) * state%yy(ig,5) + vs(6) * state%yy(ig,6)
                        state%vlinear(ig) = state%vlinear(ig) + v_s * state%aselect(k) * et(k)
                        if( ideriv .gt. 0 )then
                            do i = 1, 6
                                vsa(i) = 0.e0_wp
                                vsb(i) = 0.e0_wp
                                vsr(i) = 0.e0_wp
                            enddo
                            do j = 1, state%nll(k)
                                vsa(1) = vsa(1) + tcda(j) * state%ordcoef(j+nd0)
                                vsa(2) = vsa(2) + tcda(j) * state%ordcoef(j+nd1)
                                vsa(3) = vsa(3) + tcda(j) * state%ordcoef(j+nd2)
                                vsa(4) = vsa(4) + tcda(j) * state%ordcoef(j+nd3)
                                vsa(5) = vsa(5) + tcda(j) * state%coef(j+nd4)
                                vsa(6) = vsa(6) + teda(j) * state%coef(j+nd4)
                                vsb(1) = vsb(1) + tcdb(j) * state%ordcoef(j+nd0)
                                vsb(2) = vsb(2) + tcdb(j) * state%ordcoef(j+nd1)
                                vsb(3) = vsb(3) + tcdb(j) * state%ordcoef(j+nd2)
                                vsb(4) = vsb(4) + tcdb(j) * state%ordcoef(j+nd3)
                                vsb(5) = vsb(5) + tcdb(j) * state%coef(j+nd4)
                                vsb(6) = vsb(6) + tedb(j) * state%coef(j+nd4)
                                vsr(1) = vsr(1) + tcdr(j) * state%ordcoef(j+nd0)
                                vsr(2) = vsr(2) + tcdr(j) * state%ordcoef(j+nd1)
                                vsr(3) = vsr(3) + tcdr(j) * state%ordcoef(j+nd2)
                                vsr(4) = vsr(4) + tcdr(j) * state%ordcoef(j+nd3)
                                vsr(5) = vsr(5) + tcdr(j) * state%coef(j+nd4)
                                vsr(6) = vsr(6) + tedr(j) * state%coef(j+nd4)
                            enddo
                            do j = 1, 6
                                dvlindri(ig,j) = dvlindri(ig,j) + ( state%aselect(k) * detdri(k,j) + dadri(k,j) * &
                                    & et(k) ) * v_s
                            enddo
                            dvlindgeo(ig,1) = dvlindgeo(ig,1) + ( vsa(1) * state%yy(ig,1) + vsa(2) * state%yy(ig,2) + &
                                & vsa(3) * state%yy(ig,3) + vsa(4) * state%yy(ig,4) + vsa(5) * state%yy(ig,5) + &
                                & vsa(6) * state%yy(ig,6) ) * state%aselect(k) * et(k)
                            dvlindgeo(ig,2) = dvlindgeo(ig,2) + ( vsb(1) * state%yy(ig,1) + vsb(2) * state%yy(ig,2) + &
                                & vsb(3) * state%yy(ig,3) + vsb(4) * state%yy(ig,4) + vsb(5) * state%yy(ig,5) + &
                                & vsb(6) * state%yy(ig,6) ) * state%aselect(k) * et(k)
                            dvlindgeo(ig,3) = dvlindgeo(ig,3) + ( vsr(1) * state%yy(ig,1) + vsr(2) * state%yy(ig,2) + &
                                & vsr(3) * state%yy(ig,3) + vsr(4) * state%yy(ig,4) + vsr(5) * state%yy(ig,5) + &
                                & vsr(6) * state%yy(ig,6) ) * state%aselect(k) * et(k)
                            dvlindgeo(ig,4) = dvlindgeo(ig,4) + ( vs(2) * state%dyydgeo(ig,2,4) + vs(3) * &
                                & state%dyydgeo(ig,3,4) + vs(4) * state%dyydgeo(ig,4,4) + vs(5) * state%dyydgeo(ig,5, &
                                & 4) ) * state%aselect(k) * et(k)
                            dvlindgeo(ig,5) = dvlindgeo(ig,5) + ( vs(2) * state%dyydgeo(ig,2,5) + vs(3) * &
                                & state%dyydgeo(ig,3,5) + vs(4) * state%dyydgeo(ig,4,5) + vs(6) * state%dyydgeo(ig,6, &
                                & 5) ) * state%aselect(k) * et(k)
                            dvlindgeo(ig,6) = dvlindgeo(ig,6) + ( vs(3) * state%dyydgeo(ig,3,6) + vs(4) * &
                                & state%dyydgeo(ig,4,6) ) * state%aselect(k) * et(k)
                        endif
                        ndone = ndone + state%ncol(k)
                        nd0   = nd0 + state%nll4(k)
                    endif
                enddo
            endif
        enddo
        if( state%njac .lt. nterms .and. state%iv(nterms) .gt. 0 )then
            call mbe_h4bmkp(state, r, ideriv, ndone, vmbe, dvmbedri)
            state%vlinmbe = vmbe * state%aselect(nterms) * et(nterms)
            if( ideriv .gt. 0 )then
                do i = 1, 6
                    dvmbedri(i) = dvmbedri(i) * state%aselect(nterms) * et(nterms) + ( dadri(nterms,i) * et(nterms) + &
                        & state%aselect(nterms) * detdri(nterms,i) ) * vmbe
                enddo
            endif
        endif
    end subroutine vlinear_h4bmkp

    !> @brief Construct geometry-selection functions for the polynomial expansion.
    subroutine asqsub_h4bmkp( state, r, ideriv, rbar, rho, aselect, dadri )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r(6)
        integer, intent(in)     :: ideriv
        real(wp), intent(inout) :: rbar
        real(wp), intent(inout) :: rho
        real(wp), intent(inout) :: aselect(8)
        real(wp), intent(inout) :: dadri(8,6)
        integer                 :: i, j
        real(wp)                :: acomp, agen, ah2h2, ah3h, akite, anonlinh3h, asq0, asq2, asq2x, asqden, atet, atetden, atetx, ddd, ddr1, ddr2, ddr3, dds, r1, r1256, r1346, r1pr6, r2, r2345, r2pr5, r3, r3pr4, r4, r5, r6, rat16sq, rat25sq, rat34sq, rat43sq, rat52sq, rat61sq, rbar16, rbar2, rbar4, rbar8, ri, ta, ta1d, ta2d, tb, tb1d, tb2d, tc, tc1d, tc2d, td, td1d, td2d, te, te1d, te2d, tf, tf1d, tf2d, v16d, v16s, v25d, v25s, v34d, v34s, w1, w2, w3, w4, wn1, wn1_w1_11, wn1a, wn1b, wn1c, wn2, wn2_w2_11, wn2a, wn2b, wn2c, wn3, wn3_w3_11, wn3a, wn3b, wn3c, wn4, wn4_w4_11, wn4a, wn4b, wn4c
        real(wp), parameter     :: fside      = 0.87867965644035752e0_wp
        real(wp), parameter     :: fs2        = 2.e0_wp*fside
        real(wp), parameter     :: fdiag      = 1.2426406871192852e0_wp
        real(wp), parameter     :: fd2        = 2.e0_wp*fdiag
        real(wp), parameter     :: rr         = 1.15e0_wp
        real(wp), parameter     :: rt         = 1.25e0_wp
        real(wp), parameter     :: delta      = 0.1e0_wp
        real(wp), parameter     :: f_1_3      = 1.e0_wp / 3.e0_wp
        real(wp), parameter     :: f_2_3      = 2.e0_wp / 3.e0_wp
        real(wp), parameter     :: f13fs2     = f_1_3 * fs2
        real(wp), parameter     :: twomf13fs2 = 2.e0_wp - f13fs2
        real(wp), parameter     :: f13fd2     = f_1_3 * fd2
        real(wp), parameter     :: twomf13fd2 = 2.e0_wp - f13fd2
        real(wp), parameter     :: f_8_3      = 8.e0_wp / 3.e0_wp
        r1     = r(1)
        r2     = r(2)
        r3     = r(3)
        r4     = r(4)
        r5     = r(5)
        r6     = r(6)
        rbar   = ( r1+r2+r3+r4+r5+r6 ) / 6.e0_wp
        rbar2  = rbar * rbar
        rbar4  = rbar2 * rbar2
        rbar8  = rbar4 * rbar4
        rbar16 = rbar8 * rbar8
        rho    = sqrt( r1**2 + r2**2 + r3**2 + r4**2 + r5**2 + r6**2 )
        if( state%isqtet .gt. 0 )then
            atet                        = (r1-rbar)**2 + (r2-rbar)**2 + (r3-rbar)**2 + (r4-rbar)**2 + (r5-rbar)**2 + &
                & (r6-rbar)**2
            atetden                     = delta * rbar2 + atet
            atetx                       = rbar4 * delta / atetden
            v16s                        = r1 + r6 - fs2 * rbar
            v16d                        = r1 + r6 - fd2 * rbar
            v25s                        = r2 + r5 - fs2 * rbar
            v25d                        = r2 + r5 - fd2 * rbar
            v34s                        = r3 + r4 - fs2 * rbar
            v34d                        = r3 + r4 - fd2 * rbar
            asq0                        = ( ( v16s * v16d )**2 + ( v25s * v25d )**2 + ( v34s * v34d )**2 ) / rbar2
            asq2                        = asq0 + ( r1-r6 )**2 + ( r2-r5 )**2 + ( r3-r4 )**2
            asqden                      = delta * rbar2 + asq2
            asq2x                       = rbar4 * delta / asqden
            state%aselect(state%isqtet) = atetx + asq2x
            if( ideriv .gt. 0 )then
                ddr1 = ( twomf13fs2 * v16s * v16d**2 + twomf13fd2 * v16s**2 * v16d - f13fs2 * ( v25s * v25d**2 + &
                    & v34s * v34d**2 ) - f13fd2 * ( v25s**2 * v25d + v34s**2 * v34d ) ) / rbar2
                ddr2 = ( twomf13fs2 * v25s * v25d**2 + twomf13fd2 * v25s**2 * v25d - f13fs2 * ( v16s * v16d**2 + &
                    & v34s * v34d**2 ) - f13fd2 * ( v16s**2 * v16d + v34s**2 * v34d ) ) / rbar2
                ddr3 = ( twomf13fs2 * v34s * v34d**2 + twomf13fd2 * v34s**2 * v34d - f13fs2 * ( v25s * v25d**2 + &
                    & v16s * v16d**2 ) - f13fd2 * ( v25s**2 * v25d + v16s**2 * v16d ) ) / rbar2
                dds                   = - asq2x / asqden
                dadri(state%isqtet,4) = ( ddr3 + 2.e0_wp * ( r4 - r3 ) ) * dds
                dadri(state%isqtet,5) = ( ddr2 + 2.e0_wp * ( r5 - r2 ) ) * dds
                dadri(state%isqtet,6) = ( ddr1 + 2.e0_wp * ( r6 - r1 ) ) * dds
                dadri(state%isqtet,1) = ( ddr1 + 2.e0_wp * ( r1 - r6 ) ) * dds
                dadri(state%isqtet,2) = ( ddr2 + 2.e0_wp * ( r2 - r5 ) ) * dds
                dadri(state%isqtet,3) = ( ddr3 + 2.e0_wp * ( r3 - r4 ) ) * dds
                ddd = ( f_2_3 / rbar + f_1_3 * ( asq0 / rbar - delta * rbar ) / asqden ) * asq2x + ( f_2_3 / rbar - &
                    & delta * f_1_3 * rbar / atetden ) * atetx
                dds = 2.e0_wp * atetx / atetden
                do j = 1, 6
                    dadri(state%isqtet,j) = dadri(state%isqtet,j) + ddd - ( r(j) - rbar ) * dds
                enddo
            endif
        endif
        if( state%inonlinh3 .gt. 0 .or. state%iv(4) .gt. 0 )then
            w1               = ( r1 + r2 + r4 ) - ( r3 + r5 + r6 )
            w2               = ( r1 + r3 + r5 ) - ( r2 + r4 + r6 )
            w3               = ( r2 + r3 + r6 ) - ( r1 + r4 + r5 )
            w4               = ( r4 + r5 + r6 ) - ( r1 + r2 + r3 )
            ah3h             = w1 * w1 + w2 * w2 + w3 * w3 + w4 * w4
            state%aselect(4) = ah3h
            if( ideriv .gt. 0 )then
                dadri(4,1) = 2.e0_wp * ( w1 + w2 - w3 - w4 )
                dadri(4,2) = 2.e0_wp * ( w1 - w2 + w3 - w4 )
                dadri(4,3) = 2.e0_wp * ( w2 - w1 + w3 - w4 )
                dadri(4,4) = 2.e0_wp * ( w1 - w2 - w3 + w4 )
                dadri(4,5) = 2.e0_wp * ( w2 - w1 - w3 + w4 )
                dadri(4,6) = 2.e0_wp * ( w3 + w4 - w1 - w2 )
            endif
        endif
        if( state%inonlinh3 .gt. 0 )then
            wn1a = r1+r2-r4
            wn1b = r1+r4-r2
            wn1c = r2+r4-r1
            wn1  = wn1a * wn1b * wn1c
            wn2a = r1+r3-r5
            wn2b = r1+r5-r3
            wn2c = r3+r5-r1
            wn2  = wn2a * wn2b * wn2c
            wn3a = r2+r3-r6
            wn3b = r2+r6-r3
            wn3c = r3+r6-r2
            wn3  = wn3a * wn3b * wn3c
            wn4a = r4+r5-r6
            wn4b = r4+r6-r5
            wn4c = r5+r6-r4
            wn4  = wn4a * wn4b * wn4c
            anonlinh3h                     = ( wn1**2 * w1**12 + wn2**2 * w2**12 + wn3**2 * w3**12 + wn4**2 * w4** &
                & 12 ) / rbar16
            state%aselect(state%inonlinh3) = anonlinh3h
            if( ideriv .gt. 0 )then
                wn1_w1_11 = wn1 * w1**11
                wn2_w2_11 = wn2 * w2**11
                wn3_w3_11 = wn3 * w3**11
                wn4_w4_11 = wn4 * w4**11
                dadri(state%inonlinh3,1) = ( 2.e0_wp * ( wn1b * wn1c + wn1a * wn1c - wn1a * wn1b ) * w1 + 12.e0_wp * &
                    & wn1 ) * wn1_w1_11 + ( 2.e0_wp * ( wn2b * wn2c + wn2a * wn2c - wn2a * wn2b ) * w2 + 12.e0_wp * &
                    & wn2 ) * wn2_w2_11 - 12.e0_wp * ( wn3 * wn3_w3_11 + wn4 * wn4_w4_11 )
                dadri(state%inonlinh3,2) = ( 2.e0_wp * ( wn1b * wn1c - wn1a * wn1c + wn1a * wn1b ) * w1 + 12.e0_wp * &
                    & wn1 ) * wn1_w1_11 + ( 2.e0_wp * ( wn3b * wn3c + wn3a * wn3c - wn3a * wn3b ) * w3 + 12.e0_wp * &
                    & wn3 ) * wn3_w3_11 - 12.e0_wp * ( wn2 * wn2_w2_11 + wn4 * wn4_w4_11 )
                dadri(state%inonlinh3,3) = ( 2.e0_wp * ( wn2b * wn2c - wn2a * wn2c + wn2a * wn2b ) * w2 + 12.e0_wp * &
                    & wn2 ) * wn2_w2_11 + ( 2.e0_wp * ( wn3b * wn3c - wn3a * wn3c + wn3a * wn3b ) * w3 + 12.e0_wp * &
                    & wn3 ) * wn3_w3_11 - 12.e0_wp * ( wn1 * wn1_w1_11 + wn4 * wn4_w4_11 )
                dadri(state%inonlinh3,4) = ( 2.e0_wp * ( wn1a * wn1b + wn1a * wn1c - wn1b * wn1c ) * w1 + 12.e0_wp * &
                    & wn1 ) * wn1_w1_11 + ( 2.e0_wp * ( wn4b * wn4c + wn4a * wn4c - wn4a * wn4b ) * w4 + 12.e0_wp * &
                    & wn4 ) * wn4_w4_11 - 12.e0_wp * ( wn2 * wn2_w2_11 + wn3 * wn3_w3_11 )
                dadri(state%inonlinh3,5) = ( 2.e0_wp * ( wn2a * wn2b + wn2a * wn2c - wn2b * wn2c ) * w2 + 12.e0_wp * &
                    & wn2 ) * wn2_w2_11 + ( 2.e0_wp * ( wn4b * wn4c - wn4a * wn4c + wn4a * wn4b ) * w4 + 12.e0_wp * &
                    & wn4 ) * wn4_w4_11 - 12.e0_wp * ( wn1 * wn1_w1_11 + wn3 * wn3_w3_11 )
                dadri(state%inonlinh3,6) = ( 2.e0_wp * ( wn3a * wn3b + wn3a * wn3c - wn3b * wn3c ) * w3 + 12.e0_wp * &
                    & wn3 ) * wn3_w3_11 + ( 2.e0_wp * ( wn4a * wn4b + wn4a * wn4c - wn4b * wn4c ) * w4 + 12.e0_wp * &
                    & wn4 ) * wn4_w4_11 - 12.e0_wp * ( wn1 * wn1_w1_11 + wn2 * wn2_w2_11 )
                ddd = f_8_3 * anonlinh3h / rbar
                do j = 1, 6
                    dadri(state%inonlinh3,j) = dadri(state%inonlinh3,j) / rbar16 - ddd
                enddo
            endif
        endif
        if( state%iv(7) .gt. 0 )then
            r1pr6            = r1+r6
            r2pr5            = r2+r5
            r3pr4            = r3+r4
            r2345            = r2+r3+r4+r5
            r1346            = r1+r3+r4+r6
            r1256            = r1+r2+r5+r6
            ah2h2            = r1pr6**2 - r2345**2 + r2pr5**2 - r1346**2 + r3pr4**2 - r1256**2
            state%aselect(7) = ah2h2
            if( ideriv .gt. 0 )then
                dadri(7,1) = 2.e0_wp * ( r1pr6 - r1346 - r1256 )
                dadri(7,2) = 2.e0_wp * ( r2pr5 - r2345 - r1256 )
                dadri(7,3) = 2.e0_wp * ( r3pr4 - r2345 - r1346 )
                dadri(7,4) = 2.e0_wp * ( r3pr4 - r2345 - r1346 )
                dadri(7,5) = 2.e0_wp * ( r2pr5 - r2345 - r1256 )
                dadri(7,6) = 2.e0_wp * ( r1pr6 - r1346 - r1256 )
            endif
        endif
        if( state%igen1 .gt. 0 )then
            agen                       = 1.e0_wp
            state%aselect(state%igen1) = agen
            if( ideriv .gt. 0 )then
                do j = 1, 6
                    dadri(state%igen1,j) = 0.e0_wp
                enddo
            endif
        endif
        if( state%iv(2) .gt. 0 )then
            ta1d             = (r1-r5)**2 + rbar2
            ta2d             = (r2-r6)**2 + rbar2
            rat43sq          = (r4/r3)**2
            ta               = rat43sq * ( rbar2 / ta1d + rbar2 / ta2d )
            tb1d             = (r1-r2)**2 + rbar2
            tb2d             = (r5-r6)**2 + rbar2
            rat34sq          = (r3/r4)**2
            tb               = rat34sq * ( rbar2 / tb1d + rbar2 / tb2d )
            tc1d             = (r2-r4)**2 + rbar2
            tc2d             = (r3-r5)**2 + rbar2
            rat61sq          = (r6/r1)**2
            tc               = rat61sq * ( rbar2 / tc1d + rbar2 / tc2d )
            td1d             = (r2-r3)**2 + rbar2
            td2d             = (r4-r5)**2 + rbar2
            rat16sq          = (r1/r6)**2
            td               = rat16sq * ( rbar2 / td1d + rbar2 / td2d )
            te1d             = (r1-r4)**2 + rbar2
            te2d             = (r3-r6)**2 + rbar2
            rat52sq          = (r5/r2)**2
            te               = rat52sq * ( rbar2 / te1d + rbar2 / te2d )
            tf1d             = (r1-r3)**2 + rbar2
            tf2d             = (r4-r6)**2 + rbar2
            rat25sq          = (r2/r5)**2
            tf               = rat25sq * ( rbar2 / tf1d + rbar2 / tf2d )
            akite            = ta + tb + tc + td + te + tf
            state%aselect(2) = akite
            if( ideriv .gt. 0 )then
                ddd = f_1_3 * akite / rbar - rbar2 * f_1_3 * rbar * ( rat43sq / ta1d**2 + rat43sq / ta2d**2 + &
                    & rat34sq / tb1d**2 + rat34sq / tb2d**2 + rat61sq / tc1d**2 + rat61sq / tc2d**2 + rat16sq / td1d** &
                    & 2 + rat16sq / td2d**2 + rat52sq / te1d**2 + rat52sq / te2d**2 + rat25sq / tf1d**2 + rat25sq / &
                    & tf2d**2 )
                dadri(2,1) = ddd + 2.e0_wp * ( ( td - tc ) / r1 - ( rat43sq * (r1-r5) / ta1d**2 + rat34sq * (r1-r2) / &
                    & tb1d**2 + rat52sq * (r1-r4) / te1d**2 + rat25sq * (r1-r3) / tf1d**2 ) * rbar2 )
                dadri(2,2) = ddd + 2.e0_wp * ( ( tf - te ) / r2 - ( rat43sq * (r2-r6) / ta2d**2 - rat34sq * (r1-r2) / &
                    & tb1d**2 + rat61sq * (r2-r4) / tc1d**2 + rat16sq * (r2-r3) / td1d**2 ) * rbar2 )
                dadri(2,3) = ddd + 2.e0_wp * ( ( tb - ta ) / r3 - ( rat61sq * (r3-r5) / tc2d**2 - rat16sq * (r2-r3) / &
                    & td1d**2 + rat52sq * (r3-r6) / te2d**2 - rat25sq * (r1-r3) / tf1d**2 ) * rbar2 )
                dadri(2,4) = ddd + 2.e0_wp * ( ( ta - tb ) / r4 - ( rat16sq * (r4-r5) / td2d**2 - rat61sq * (r2-r4) / &
                    & tc1d**2 - rat52sq * (r1-r4) / te1d**2 + rat25sq * (r4-r6) / tf2d**2 ) * rbar2 )
                dadri(2,5) = ddd + 2.e0_wp * ( ( te - tf ) / r5 - ( rat34sq * (r5-r6) / tb2d**2 - rat43sq * (r1-r5) / &
                    & ta1d**2 - rat61sq * (r3-r5) / tc2d**2 - rat16sq * (r4-r5) / td2d**2 ) * rbar2 )
                dadri(2,6) = ddd + 2.e0_wp * ( ( tc - td ) / r6 - ( rat43sq * (r6-r2) / ta2d**2 - rat34sq * (r5-r6) / &
                    & tb2d**2 - rat52sq * (r3-r6) / te2d**2 - rat25sq * (r4-r6) / tf2d**2 ) * rbar2 )
            endif
        endif
        if( state%icomp .gt. 0 )then
            acomp = 0.e0_wp
            do i = 1, 6
                dadri(state%icomp,i) = 0.e0_wp
                if( r(i) .lt. rr )then
                    ri                                       = r(i)
                    acomp                                    = acomp + ( rr-ri )**3 / (rt-ri)
                    if( ideriv .gt. 0 ) dadri(state%icomp,i) = (rr-ri)**3 / (rt-ri)**2 - 3.e0_wp * (rr-ri)**2 / (rt-ri)
                endif
            enddo
            state%aselect(state%icomp) = acomp
        endif
        if( state%igen2 .gt. 0 )then
            state%aselect(state%igen2) = 1.e0_wp
            if( ideriv .gt. 0 )then
                do j = 1, 6
                    dadri(state%igen2,j) = 0.e0_wp
                enddo
            endif
        endif
    end subroutine asqsub_h4bmkp

    !> @brief Evaluate the permutationally symmetric four-body expansion.
    subroutine mbe_h4bmkp( state, r, ideriv, ndone, vmbe, dvmbedri )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r(6)
        integer, intent(in)     :: ideriv
        integer, intent(in)     :: ndone
        real(wp), intent(inout) :: vmbe
        real(wp), intent(inout) :: dvmbedri(6)
        integer                 :: i, i1, i2, i3, i4, i5, i6, j, kcoef, kind3, kind4, kind5, kind6
        real(wp)                :: dlnpdri(6), dlnvdri(24,6), pi(6), ppow(0:12,6), val_1(24), val_2(24), val_3(24), val_4(24), val_5(24), val_6(24), value
        do i = 1, 6
            pi(i)     = r(i) * exp( -state%beta_p * r(i) )
            ppow(0,i) = 1.e0_wp
            ppow(1,i) = pi(i)
            do j = 2, state%indmax
                ppow(j,i) = ppow(j-1,i) * pi(i)
            enddo
            dvmbedri(i) = 0.e0_wp
            dlnpdri(i)  = 1.e0_wp / r(i) - state%beta_p
        enddo
        vmbe  = 0.e0_wp
        kcoef = ndone
        kind3 = 0
        kind4 = 0
        kind5 = 0
        kind6 = 0
        do i1 = state%k1_lo, state%k1_hi
            do i = 1, 24
                val_1(i) = ppow( i1,state%mbe_perm(1,i) )
            enddo
            if( ideriv .gt. 0 )then
                if( i1 .eq. 0 )then
                    do i = 1, 24
                        dlnvdri( i,state%mbe_perm(1,i) ) = 0.e0_wp
                    enddo
                else
                    do i = 1, 24
                        dlnvdri( i,state%mbe_perm(1,i) ) = i1 * dlnpdri( state%mbe_perm(1,i) )
                    enddo
                endif
            endif
            do i2 = state%k2_lo(i1), state%k2_hi(i1)
                do i = 1, 24
                    val_2(i) = val_1(i) * ppow( i2,state%mbe_perm(2,i) )
                enddo
                if( ideriv .gt. 0 )then
                    if( i2 .eq. 0 )then
                        do i = 1, 24
                            dlnvdri( i,state%mbe_perm(2,i) ) = 0.e0_wp
                        enddo
                    else
                        do i = 1, 24
                            dlnvdri( i,state%mbe_perm(2,i) ) = i2 * dlnpdri( state%mbe_perm(2,i) )
                        enddo
                    endif
                endif
                kind3 = kind3 + 1
                do i3 = state%k3_lo(kind3), state%k3_hi(kind3)
                    do i = 1, 24
                        val_3(i) = val_2(i) * ppow( i3,state%mbe_perm(3,i) )
                    enddo
                    if( ideriv .gt. 0 )then
                        if( i3 .eq. 0 )then
                            do i = 1, 24
                                dlnvdri( i,state%mbe_perm(3,i) ) = 0.e0_wp
                            enddo
                        else
                            do i = 1, 24
                                dlnvdri( i,state%mbe_perm(3,i) ) = i3 * dlnpdri( state%mbe_perm(3,i) )
                            enddo
                        endif
                    endif
                    kind4 = kind4 + 1
                    do i4 = state%k4_lo(kind4), state%k4_hi(kind4)
                        do i = 1, 24
                            val_4(i) = val_3(i) * ppow( i4,state%mbe_perm(4,i) )
                        enddo
                        if( ideriv .gt. 0 )then
                            if( i4 .eq. 0 )then
                                do i = 1, 24
                                    dlnvdri( i,state%mbe_perm(4,i) ) = 0.e0_wp
                                enddo
                            else
                                do i = 1, 24
                                    dlnvdri( i,state%mbe_perm(4,i) ) = i4 * dlnpdri( state%mbe_perm(4,i) )
                                enddo
                            endif
                        endif
                        kind5 = kind5 + 1
                        do i5 = state%k5_lo(kind5), state%k5_hi(kind5)
                            do i = 1, 24
                                val_5(i) = val_4(i) * ppow( i5,state%mbe_perm(5,i) )
                            enddo
                            if( ideriv .gt. 0 )then
                                if( i5 .eq. 0 )then
                                    do i = 1, 24
                                        dlnvdri( i,state%mbe_perm(5,i) ) = 0.e0_wp
                                    enddo
                                else
                                    do i = 1, 24
                                        dlnvdri( i,state%mbe_perm(5,i) ) = i5 * dlnpdri( state%mbe_perm(5,i) )
                                    enddo
                                endif
                            endif
                            kind6 = kind6 + 1
                            do i6 = state%k6_lo(kind6), state%k6_hi(kind6)
                                kcoef = kcoef + 1
                                if( state%coef(kcoef) .ne. 0.e0_wp )then
                                    value = 0.e0_wp
                                    do i = 1, 24
                                        val_6(i) = val_5(i) * ppow( i6,state%mbe_perm(6,i) )
                                        value    = value + val_6(i)
                                    enddo
                                    vmbe = vmbe + value * state%coef(kcoef)
                                    if( ideriv .gt. 0 )then
                                        if( i6 .eq. 0 )then
                                            do i = 1, 24
                                                dlnvdri( i,state%mbe_perm(6,i) ) = 0.e0_wp
                                            enddo
                                        else
                                            do i = 1, 24
                                                dlnvdri( i,state%mbe_perm(6,i) ) = i6 * dlnpdri( state%mbe_perm(6,i) )
                                            enddo
                                        endif
                                        do j = 1, 6
                                            value = 0.e0_wp
                                            do i = 1, 24
                                                value = value + val_6(i) * dlnvdri(i,j)
                                            enddo
                                            dvmbedri(j) = dvmbedri(j) + value * state%coef(kcoef)
                                        enddo
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine mbe_h4bmkp

    !> @brief Evaluate one H3 BKMP2 building block and its derivatives.
    subroutine vxyzcal_h4bmkp( state, rt1, rt2, rt3, rh1, rh2, rh3, ideriv, vxyz, dv )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: rt1
        real(wp), intent(in)    :: rt2
        real(wp), intent(in)    :: rt3
        real(wp), intent(in)    :: rh1
        real(wp), intent(in)    :: rh2
        real(wp), intent(in)    :: rh3
        integer, intent(in)     :: ideriv
        real(wp), intent(inout) :: vxyz
        real(wp), intent(inout) :: dv(6)
        real(wp)                :: calasbe, dcalasbe(3), derfac, dvasbe(3), es1, es2, f1, f2, roh, rohsq, rot, rotsq, rt(3), s_1, s_2, vasbe
        rotsq = 1.e0_wp/rt1**2 + 1.e0_wp/rt2**2 + 1.e0_wp/rt3**2
        rohsq = 1.e0_wp/rh1**2 + 1.e0_wp/rh2**2 + 1.e0_wp/rh3**2
        rot   = sqrt( rotsq )
        roh   = sqrt( rohsq )
        es1   = state%s0_nsb_1 * exp( state%st_nsb_1 / rot - state%sh_nsb_1 / roh )
        es2   = state%s0_nsb_2 * exp( state%st_nsb_2 / rot - state%sh_nsb_2 / roh )
        s_1   = 1.e0_wp / ( 1.e0_wp + es1 )
        s_2   = 1.e0_wp / ( 1.e0_wp + es2 )
        rt(1) = rt1
        rt(2) = rt2
        rt(3) = rt3
        call h3bend_h4bmkp(state, rt, ideriv, vasbe, dvasbe, calasbe, dcalasbe)
        vxyz = s_1 * vasbe + s_2 * calasbe
        if( ideriv .gt. 0 )then
            f1     = vasbe * s_1**2 * es1
            f2     = calasbe * s_2**2 * es2
            derfac = ( f1 * state%st_nsb_1 + f2 * state%st_nsb_2 ) / ( rot * rotsq )
            dv(1)  = s_1 * dvasbe(1) + s_2 * dcalasbe(1) - derfac / rt1**3
            dv(2)  = s_1 * dvasbe(2) + s_2 * dcalasbe(2) - derfac / rt2**3
            dv(3)  = s_1 * dvasbe(3) + s_2 * dcalasbe(3) - derfac / rt3**3
            derfac = ( f1 * state%sh_nsb_1 + f2 * state%sh_nsb_2 ) / ( roh * rohsq )
            dv(4)  = derfac / rh1**3
            dv(5)  = derfac / rh2**3
            dv(6)  = derfac / rh3**3
        endif
    end subroutine vxyzcal_h4bmkp

    !> @brief Combine H3 antisymmetric, compact, and bending corrections.
    subroutine h3bend_h4bmkp( state, r, ideriv, vasbe, dvasbe, calasbe, dcalasbe )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r(3)
        integer, intent(in)     :: ideriv
        real(wp), intent(inout) :: vasbe
        real(wp), intent(inout) :: dvasbe(3)
        real(wp), intent(inout) :: calasbe
        real(wp), intent(inout) :: dcalasbe(3)
        integer                 :: i, icompc
        real(wp)                :: acalc, cal, cas, cbnda, cbndb, dacalc(3), dcal(3), dcas(3), dcbnda(3), dcbndb(3), dt(3), dvas(3), dvbnda(3), dvbndb(3), dvlon(3), t(3), vas, vbnda, vbndb, vlon
        vlon  = 0.e0_wp
        vas   = 0.e0_wp
        vbnda = 0.e0_wp
        vbndb = 0.e0_wp
        cal   = 0.e0_wp
        cas   = 0.e0_wp
        cbnda = 0.e0_wp
        cbndb = 0.e0_wp
        do i = 1, 3
            dvlon(i)    = 0.e0_wp
            dvas(i)     = 0.e0_wp
            dvbnda(i)   = 0.e0_wp
            dvbndb(i)   = 0.e0_wp
            dcal(i)     = 0.e0_wp
            dcas(i)     = 0.e0_wp
            dcbnda(i)   = 0.e0_wp
            dcbndb(i)   = 0.e0_wp
            dvasbe(i)   = 0.e0_wp
            dcalasbe(i) = 0.e0_wp
            dacalc(i)   = 0.e0_wp
        enddo
        call vascal95_h4bmkp(state, r, vas, dvas, ideriv, acalc, dacalc)
        call compac95_h4bmkp(state, r, icompc, t, dt)
        if( icompc .gt. 0 )then
            call csym95_h4bmkp(state, r, cal, dcal)
            call casym95_h4bmkp(state, r, cas, dcas, ideriv, t, dt, acalc, dacalc)
        endif
        call vbcb95_h4bmkp(state, r, icompc, t, dt, ideriv, vbnda, vbndb, dvbnda, dvbndb, cbnda, cbndb, dcbnda, dcbndb)
        vasbe   = vas + vbnda + vbndb
        calasbe = cal + cas + cbnda + cbndb
        if( ideriv .gt. 0 )then
            dvasbe(1)   = dvas(1) + dvbnda(1) + dvbndb(1)
            dvasbe(2)   = dvas(2) + dvbnda(2) + dvbndb(2)
            dvasbe(3)   = dvas(3) + dvbnda(3) + dvbndb(3)
            dcalasbe(1) = dcal(1) + dcas(1) + dcbnda(1) + dcbndb(1)
            dcalasbe(2) = dcal(2) + dcas(2) + dcbnda(2) + dcbndb(2)
            dcalasbe(3) = dcal(3) + dcas(3) + dcbnda(3) + dcbndb(3)
        endif
    end subroutine h3bend_h4bmkp

    !> @brief Evaluate the 1995 triplet H2 curve with its short-range splice.
    subroutine triplet95_h4bmkp( state, r, e3, ideriv )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r
        real(wp), intent(inout) :: e3(3)
        integer, intent(in)     :: ideriv
        real(wp)                :: dr, e1(3), exdr, ra6, ra61, rsq
        real(wp), parameter     :: rl = 0.95e0_wp
        real(wp), parameter     :: rr = 1.15e0_wp
        real(wp), parameter     :: z1 = 1.e0_wp
        real(wp), parameter     :: z2 = 2.e0_wp
        real(wp), parameter     :: a1 = -0.0298546962e0_wp
        real(wp), parameter     :: a2 = -23.9604445036e0_wp
        real(wp), parameter     :: a3 = -42.5185569474e0_wp
        real(wp), parameter     :: a4 = 2.0382390988e0_wp
        real(wp), parameter     :: a5 = -11.5214861455e0_wp
        real(wp), parameter     :: a6 = 1.5309487826e0_wp
        real(wp), parameter     :: c1 = -0.4106358351531854e0_wp
        real(wp), parameter     :: c2 = -0.0770355790707090e0_wp
        real(wp), parameter     :: c3 = 0.4303193846943223e0_wp
        e3(3) = 0.e0_wp
        e3(2) = 0.e0_wp
        if( r .ge. rr )then
            exdr  = exp( -a4*r )
            rsq   = r*r
            ra6   = r**(-a6)
            e3(1) = a1* ( a2 + r + a3*rsq + a5*ra6 )*exdr
            if( ideriv .gt. 0 )then
                ra61  = r**(-a6-z1)
                e3(2) = a1*exdr * ( z1 -a2*a4 +(z2*a3-a4)*r -a3*a4*rsq -a5*a6*ra61 -a4*a5*ra6 )
            endif
        else
            dr = r - rl
            call vh2opt95_h4bmkp(state, r, e1, ideriv)
            if( r .le. rl )then
                e3(1)                     = e1(1) + c2*dr + c3
                if( ideriv .gt. 0 ) e3(2) = e1(2) + c2
            else
                e3(1)                     = e1(1) + c1*dr*dr*dr + c2*dr + c3
                if( ideriv .gt. 0 ) e3(2) = e1(2) + 3.e0_wp*c1*dr*dr + c2
            endif
        endif
    end subroutine triplet95_h4bmkp

    !> @brief Evaluate the optimized 1995 singlet H2 reference potential.
    subroutine vh2opt95_h4bmkp( state, r, e, ideriv )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r
        real(wp), intent(inout) :: e(3)
        integer, intent(in)     :: ideriv
        real(wp)                :: alphar, dalphr, dvlr, dvsr, exalph, r02, r04, r06, r10, r11, r12, r13, r14, r15, r2, r3, r4, r5, r6, r7, r8, r9, rr2, rr25, rr26, rr4, rr43, rr6, vlr, vsr
        real(wp), parameter     :: a0  = 0.03537359271649620e0_wp
        real(wp), parameter     :: a1  = 2.013977588700072e0_wp
        real(wp), parameter     :: a2  = -2.827452449964767e0_wp
        real(wp), parameter     :: a3  = 2.713257715593500e0_wp
        real(wp), parameter     :: a4  = -2.792039234205731e0_wp
        real(wp), parameter     :: a5  = 2.166542078766724e0_wp
        real(wp), parameter     :: a6  = -1.272679684173909e0_wp
        real(wp), parameter     :: a7  = 0.5630423099212294e0_wp
        real(wp), parameter     :: a8  = -0.1879397372273814e0_wp
        real(wp), parameter     :: a9  = 0.04719891893374140e0_wp
        real(wp), parameter     :: a10 = -0.008851622656489644e0_wp
        real(wp), parameter     :: a11 = 0.001224998776243630e0_wp
        real(wp), parameter     :: a12 = -1.227820520228028e-04_wp
        real(wp), parameter     :: a13 = 8.638783190083473e-06_wp
        real(wp), parameter     :: a14 = -4.036967926499151e-07_wp
        real(wp), parameter     :: a15 = 1.123286608335365e-08_wp
        real(wp), parameter     :: a16 = -1.406619156782167e-10_wp
        real(wp), parameter     :: r0  = 3.5284882e0_wp
        real(wp), parameter     :: dd  = 0.160979391e0_wp
        real(wp), parameter     :: c6  = 6.499027e0_wp
        real(wp), parameter     :: c8  = 124.3991e0_wp
        real(wp), parameter     :: c10 = 3285.828e0_wp
        e(2) = 0.e0_wp
        e(3) = 0.e0_wp
        r2   = r  *r
        r3   = r2 *r
        r4   = r3 *r
        r5   = r4 *r
        r6   = r5 *r
        r7   = r6 *r
        r8   = r7 *r
        r9   = r8 *r
        r10  = r9 *r
        r11  = r10*r
        r12  = r11*r
        r13  = r12*r
        r14  = r13*r
        r15  = r14*r
        r02  = r0*r0
        r04  = r02*r02
        r06  = r04*r02
        rr2  = r2 + r02
        rr4  = r4 + r04
        rr6  = r6 + r06
        rr25 = rr2*rr2*rr2*rr2*rr2
        alphar =  a0/r + a1 + a2 *r   + a3 *r2  + a4 *r3  + a5 *r4  + a6 *r5 + a7 *r6  + a8 *r7  + a9 *r8  + a10*r9  + &
            & a11* r10 + a12*r11 + a13*r12 + a14*r13 + a15*r14 + a16*r15
        exalph = exp(alphar)
        vsr    = dd*(exalph-1.e0_wp)*(exalph-1.e0_wp) - dd
        vlr    = -c6/rr6   -c8/(rr4*rr4)   -c10/rr25
        e(1)   = vsr +  vlr
        if( ideriv .gt. 0 )then
            rr26 = rr25*rr2
            rr43 = rr4*rr4*rr4
            dalphr = -a0/r2 + a2 + 2.e0_wp*a3*r + 3.e0_wp*a4*r2 + 4.e0_wp*a5*r3 + 5.e0_wp*a6*r4 + 6.e0_wp*a7*r5 + &
                & 7.e0_wp*a8*r6 + 8.e0_wp*a9*r7 + 9.e0_wp*a10*r8 + 10.e0_wp*a11*r9 + 11.e0_wp*a12*r10 + 12.e0_wp*a13* &
                & r11 + 13.e0_wp*a14*r12 + 14.e0_wp*a15*r13 + 15.e0_wp*a16*r14
            dvsr = 2.e0_wp*dd*(exalph - 1.e0_wp)*exalph*dalphr
            dvlr = 6.e0_wp*c6*r5/(rr6*rr6) + 8.e0_wp*c8*r3/rr43 + 10.e0_wp*c10*r/rr26
            e(2) = dvsr + dvlr
        endif
    end subroutine vh2opt95_h4bmkp

    !> @brief Evaluate the antisymmetric H3 correction and analytic derivative.
    subroutine vascal95_h4bmkp( state, rpass, vas, dvas, ideriv, a, da )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: rpass(3)
        real(wp), intent(inout) :: vas
        real(wp), intent(inout) :: dvas(3)
        integer, intent(in)     :: ideriv
        real(wp), intent(inout) :: a
        real(wp), intent(inout) :: da(3)
        integer                 :: isgn_a
        real(wp)                :: a2, a3, a4, a5, ds(3), exp1, exp6, r, r1, r2, r3, rcu, rsq, s
        real(wp), parameter     :: aa1 = 0.3788951192e-02_wp
        real(wp), parameter     :: aa2 = 0.1478100901e-02_wp
        real(wp), parameter     :: aa3 = -.1848513849e-03_wp
        real(wp), parameter     :: aa4 = 0.9230803609e-05_wp
        real(wp), parameter     :: aa5 = -.1293180255e-06_wp
        real(wp), parameter     :: aa6 = 0.5237179303e+00_wp
        real(wp), parameter     :: aa7 = -.1112326215e-02_wp
        r1  = rpass(1)
        r2  = rpass(2)
        r3  = rpass(3)
        r   = r1 + r2 + r3
        rsq = r*r
        rcu = rsq*r
        a   = (r1-r2)*(r2-r3)*(r3-r1)
        if( a .lt. 0.e0_wp )then
            isgn_a = -1
            a      = -a
        else
            isgn_a = 1
        endif
        a2   = a *a
        a3   = a2*a
        a4   = a3*a
        a5   = a4*a
        exp1 = exp(-aa1*rcu)
        exp6 = exp(-aa6*r)
        s    = aa2*a2 + aa3*a3 + aa4*a4 + aa5*a5
        vas  = s*exp1 +  aa7*a2 *exp6 / r
        if( ideriv .gt. 0 )then
            if( isgn_a .gt. 0 )then
                da(1) = ( r2 + r3 - 2.e0_wp*r1 ) * (r2-r3)
                da(2) = ( r3 + r1 - 2.e0_wp*r2 ) * (r3-r1)
                da(3) = ( r1 + r2 - 2.e0_wp*r3 ) * (r1-r2)
            else
                da(1) = ( 2.e0_wp*r1 - r2 - r3 ) * (r2-r3)
                da(2) = ( 2.e0_wp*r2 - r3 - r1 ) * (r3-r1)
                da(3) = ( 2.e0_wp*r3 - r1 - r2 ) * (r1-r2)
            endif
            ds(1) = ( 2.e0_wp*aa2*a  +3.e0_wp*aa3*a2 +4.e0_wp*aa4*a3 +5.e0_wp*aa5*a4) * da(1)
            dvas(1) = -3.e0_wp*aa1*rsq*s*exp1 + ds(1)*exp1 -aa7*a2*exp6/rsq  + 2.e0_wp*aa7*a*da(1)*exp6/r -aa6*aa7*a2* &
                & exp6/r
            ds(2) = ( 2.e0_wp*aa2*a  +3.e0_wp*aa3*a2 +4.e0_wp*aa4*a3 +5.e0_wp*aa5*a4) * da(2)
            dvas(2) = -3.e0_wp*aa1*rsq*s*exp1 + ds(2)*exp1 -aa7*a2*exp6/rsq  + 2.e0_wp*aa7*a*da(2)*exp6/r -aa6*aa7*a2* &
                & exp6/r
            ds(3) = ( 2.e0_wp*aa2*a  +3.e0_wp*aa3*a2 +4.e0_wp*aa4*a3 +5.e0_wp*aa5*a4) * da(3)
            dvas(3) = -3.e0_wp*aa1*rsq*s*exp1 + ds(3)*exp1 -aa7*a2*exp6/rsq  + 2.e0_wp*aa7*a*da(3)*exp6/r -aa6*aa7*a2* &
                & exp6/r
        endif
    end subroutine vascal95_h4bmkp

    !> @brief Identify compact H3 geometries and compute compactness weights.
    subroutine compac95_h4bmkp( state, r, icompc, t, dt )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r(3)
        integer, intent(inout)  :: icompc
        real(wp), intent(inout) :: t(3)
        real(wp), intent(inout) :: dt(3)
        integer                 :: i
        real(wp)                :: bot, bot2, top, top2, top3
        real(wp), parameter     :: rr = 1.15e0_wp
        real(wp), parameter     :: rp = 1.25e0_wp
        icompc = 0
        do i = 1, 3
            if( r(i) .lt. rr )then
                icompc = icompc + 1
                top    = rr-r(i)
                bot    = rp-r(i)
                top2   = top  * top
                top3   = top2 * top
                bot2   = bot  * bot
                t(i)   = top3 / bot
                dt(i)  = ( top / bot - 3.e0_wp ) * top2 / bot
            else
                t(i)  = 0.e0_wp
                dt(i) = 0.e0_wp
            endif
        enddo
    end subroutine compac95_h4bmkp

    !> @brief Evaluate the compact antisymmetric H3 correction.
    subroutine casym95_h4bmkp( state, r, cas, dcas, ideriv, t, dt, a, da )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r(3)
        real(wp), intent(inout) :: cas
        real(wp), intent(inout) :: dcas(3)
        integer, intent(in)     :: ideriv
        real(wp), intent(in)    :: t(3)
        real(wp), intent(in)    :: dt(3)
        real(wp), intent(in)    :: a
        real(wp), intent(in)    :: da(3)
        integer                 :: i
        real(wp)                :: a2, determ(3), dpr(3), dseries(3), dsumt(3), dterm1(3), eterm, pr, pr2, pr3, series, sr, sr2, sr3, sumt, term1
        real(wp), parameter     :: u1  = 0.2210243144e+00_wp
        real(wp), parameter     :: u2  = 0.4367417579e+00_wp
        real(wp), parameter     :: u3  = 0.6994985432e-02_wp
        real(wp), parameter     :: u4  = 0.1491096501e+01_wp
        real(wp), parameter     :: u5  = 0.1602896673e+01_wp
        real(wp), parameter     :: u6  = -.2821747323e+01_wp
        real(wp), parameter     :: u7  = 0.4948310833e+00_wp
        real(wp), parameter     :: u8  = -.3540394679e-01_wp
        real(wp), parameter     :: u9  = -.3305809954e+01_wp
        real(wp), parameter     :: u10 = 0.3644382172e+01_wp
        real(wp), parameter     :: u11 = -.9997570970e+00_wp
        real(wp), parameter     :: u12 = 0.7989919534e-01_wp
        real(wp), parameter     :: u13 = -.1075807322e-02_wp
        cas     = 0.e0_wp
        dcas(1) = 0.e0_wp
        dcas(2) = 0.e0_wp
        dcas(3) = 0.e0_wp
        a2      = a*a
        sumt    = t(1) + t(2) + t(3)
        sr      = r(1) + r(2) + r(3)
        pr      = r(1) * r(2) * r(3)
        sr2     = sr*sr
        sr3     = sr2*sr
        pr2     = pr*pr
        pr3     = pr2*pr
        series  = 1.e0_wp + u4/pr2 +  u5/pr + u6 +   u7*pr +  u8*pr2 + a*(u9/pr2 + u10/pr + u11 + u12*pr + u13*pr2)
        term1   = u1/pr**u2
        eterm   = exp(-u3*sr3)
        cas     = sumt * a2 * term1 * series * eterm
        if( ideriv .gt. 0 )then
            dpr(1) = r(2)*r(3)
            dpr(2) = r(3)*r(1)
            dpr(3) = r(1)*r(2)
            do i = 1, 3
                dsumt(i)  = dt(i)
                dterm1(i) = -u1*u2*pr**(-u2 - 1.e0_wp)*dpr(i)
                determ(i) = eterm*(-3.e0_wp*u3*sr2)
                dseries(i) = dpr(i)*(-2.e0_wp*u4/pr3 - u5/pr2 + u7 + 2.e0_wp*u8*pr) + da(i)*(u9/pr2 + u10/pr + u11 + &
                    & u12*pr + u13*pr2) + a*dpr(i)*(-2.e0_wp*u9/pr3 - u10/pr2 + u12 + 2.e0_wp*u13*pr)
                dcas(i) = dsumt(i)*a2*term1*series*eterm + 2.e0_wp*a*da(i)*sumt*term1*series*eterm + dterm1(i)*sumt* &
                    & a2*series*eterm + dseries(i)*sumt*a2*term1*eterm + determ(i)*sumt*a2*term1*series
            enddo
        endif
    end subroutine casym95_h4bmkp

    !> @brief Evaluate the compact symmetric H3 correction.
    subroutine csym95_h4bmkp( state, r, cal, dcal )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: r(3)
        real(wp), intent(inout) :: cal
        real(wp), intent(inout) :: dcal(3)
        integer                 :: i
        real(wp)                :: dexp3, dg(3), exp3, g(3), ri, rpri, rpri2, rrri, rrri2, rrri3, sr, sr2, sr3, sumg, sumv(3)
        real(wp), parameter     :: rr = 1.15e0_wp
        real(wp), parameter     :: rp = 1.25e0_wp
        real(wp), parameter     :: v1 = -.2071708868e+00_wp
        real(wp), parameter     :: v2 = -.5672350377e+00_wp
        real(wp), parameter     :: v3 = 0.9058780367e-02_wp
        sr    = r(1)+r(2)+r(3)
        sr2   = sr*sr
        sr3   = sr*sr2
        exp3  = exp( -v3*sr3 )
        dexp3 = -3.e0_wp*v3*sr2*exp3
        do i = 1, 3
            ri = r(i)
            if( ri .lt. rr )then
                rrri    = rr-ri
                rrri2   = rrri*rrri
                rrri3   = rrri*rrri2
                rpri    = rp-ri
                rpri2   = rpri*rpri
                sumv(i) = v1 + v1*v2*ri
                g(i)    = (rrri3/rpri)*sumv(i)
                dg(i)   = (rrri3/rpri2)*sumv(i) -3.e0_wp*(rrri2/rpri)*sumv(i) + (rrri3/rpri)*v1*v2
            else
                g(i)  = 0.e0_wp
                dg(i) = 0.e0_wp
            endif
        enddo
        sumg    = g(1) + g(2) + g(3)
        cal     = sumg*exp3
        dcal(1) = dg(1)*exp3 + sumg*dexp3
        dcal(2) = dg(2)*exp3 + sumg*dexp3
        dcal(3) = dg(3)*exp3 + sumg*dexp3
    end subroutine csym95_h4bmkp

    !> @brief Evaluate H3 bending corrections for both adiabatic components.
    subroutine vbcb95_h4bmkp( state, rpass, icompc, t, dt, ideriv, vbnda, vbndb, dvbnda, dvbndb, cbnda, cbndb, dcbnda, dcbndb )
        implicit none
        type(h4_workspace_t), intent(inout) :: state
        real(wp), intent(in)    :: rpass(3)
        integer, intent(in)     :: icompc
        real(wp), intent(in)    :: t(3)
        real(wp), intent(in)    :: dt(3)
        integer, intent(in)     :: ideriv
        real(wp), intent(inout) :: vbnda
        real(wp), intent(inout) :: vbndb
        real(wp), intent(inout) :: dvbnda(3)
        real(wp), intent(inout) :: dvbndb(3)
        real(wp), intent(inout) :: cbnda
        real(wp), intent(inout) :: cbndb
        real(wp), intent(inout) :: dcbnda(3)
        real(wp), intent(inout) :: dcbndb(3)
        real(wp)                :: aasum, asum, b1, b12, b13, b14, b15, b1a, b1b, b2, b3, b3b, bsum, c1, c1cube, c2, c2cube, c3, c3cube, cb(2,25), cbnds(2), cos3t1, cos3t2, cos3t3, csum, cx1, d1, d2, d3, daasum, dasum, db1a(3), db1b(3), db2(3), db3(6), dbsum, dc1dr1, dc1dr2, dc1dr3, dc2dr1, dc2dr2, dc2dr3, dc3dr1, dc3dr2, dc3dr3, dcb1(3), dcb2(3), dcb3(3), dcb4(3), dcb5(3), dcb6(3), dcb7(3), dcb8(3), dcbnds(2,3), ddsum1, ddsum2, dexp1, dexp2, dexp7, dffsum, dfsum, dgsum, dp(3), dsum1, dsum2, dvb1(3), dvb2(3), dvb3(3), dvb4(3), dvb5(3), dx1, eps2, exp1, exp2, exp7, ffsum, fsum, gsum, p, pcu, psq, r, r1, r2, r3, rcu, rsq, sum, sumb, sumt, t1, t2, t3, vb(2,25), vbnd(2)
        real(wp), parameter     :: z58   = 0.625e0_wp
        real(wp), parameter     :: z38   = 0.375e0_wp
        real(wp), parameter     :: beta1 = 0.52e0_wp
        real(wp), parameter     :: beta2 = 0.052e0_wp
        real(wp), parameter     :: beta3 = 0.79e0_wp
        real(wp), parameter     :: a11   = -.1838073394e+03_wp
        real(wp), parameter     :: a12   = 0.1334593242e+02_wp
        real(wp), parameter     :: a13   = -.2358129537e+00_wp
        real(wp), parameter     :: a21   = -.4668193478e+01_wp
        real(wp), parameter     :: a22   = 0.7197506670e+01_wp
        real(wp), parameter     :: a23   = 0.2162004275e+02_wp
        real(wp), parameter     :: a24   = 0.2106294028e+02_wp
        real(wp), parameter     :: a31   = 0.4242962586e+01_wp
        real(wp), parameter     :: a32   = 0.4453505045e+01_wp
        real(wp), parameter     :: a41   = -.1456918088e+00_wp
        real(wp), parameter     :: a42   = -.1692657366e-01_wp
        real(wp), parameter     :: a43   = 0.1279520698e+01_wp
        real(wp), parameter     :: a44   = -.4898940075e+00_wp
        real(wp), parameter     :: a51   = 0.1742295219e+03_wp
        real(wp), parameter     :: a52   = 0.3142175348e+02_wp
        real(wp), parameter     :: a53   = 0.5152903406e+01_wp
        real(wp), parameter     :: g11   = -.4765732725e+02_wp
        real(wp), parameter     :: g12   = 0.3648933563e+01_wp
        real(wp), parameter     :: g13   = -.7141145244e-01_wp
        real(wp), parameter     :: g21   = 0.1002349176e-01_wp
        real(wp), parameter     :: g22   = 0.9989856329e-02_wp
        real(wp), parameter     :: g23   = -.4161953634e-02_wp
        real(wp), parameter     :: g24   = 0.9075807910e-03_wp
        real(wp), parameter     :: g31   = -.2693628729e+00_wp
        real(wp), parameter     :: g32   = -.1399065763e-01_wp
        real(wp), parameter     :: g41   = -.1417634346e-01_wp
        real(wp), parameter     :: g42   = -.4870024792e-03_wp
        real(wp), parameter     :: g43   = 0.1312231847e+00_wp
        real(wp), parameter     :: g44   = -.4409850519e-01_wp
        real(wp), parameter     :: g51   = 0.5382970863e+02_wp
        real(wp), parameter     :: g52   = 0.4587102824e+01_wp
        real(wp), parameter     :: g53   = 0.1768550515e+01_wp
        real(wp), parameter     :: c11   = 0.1860299931e+04_wp
        real(wp), parameter     :: c12   = -.6134458037e+03_wp
        real(wp), parameter     :: c13   = 0.7337207161e+02_wp
        real(wp), parameter     :: c14   = -.2676717625e+04_wp
        real(wp), parameter     :: c15   = 0.1344099415e+04_wp
        real(wp), parameter     :: c21   = 0.1538913137e+03_wp
        real(wp), parameter     :: c22   = 0.4348007369e+02_wp
        real(wp), parameter     :: c23   = 0.1719720677e+03_wp
        real(wp), parameter     :: c24   = 0.2115963042e+03_wp
        real(wp), parameter     :: c31   = -.7026089414e+02_wp
        real(wp), parameter     :: c32   = -.1300938992e+03_wp
        real(wp), parameter     :: c41   = 0.1310273564e+01_wp
        real(wp), parameter     :: c42   = -.6175149574e+00_wp
        real(wp), parameter     :: c43   = -.2679089358e+02_wp
        real(wp), parameter     :: c44   = 0.5577477171e+01_wp
        real(wp), parameter     :: c51   = -.3543353539e+04_wp
        real(wp), parameter     :: c52   = -.3740709591e+03_wp
        real(wp), parameter     :: c53   = 0.7979303144e+02_wp
        real(wp), parameter     :: c61   = -.1104230585e+04_wp
        real(wp), parameter     :: c62   = 0.4603572025e+04_wp
        real(wp), parameter     :: c63   = -.5593496634e+04_wp
        real(wp), parameter     :: c71   = -.1069406434e+02_wp
        real(wp), parameter     :: c72   = 0.1021807153e+01_wp
        real(wp), parameter     :: c73   = 0.6669828341e-01_wp
        real(wp), parameter     :: c74   = 0.4168542348e+02_wp
        real(wp), parameter     :: c75   = 0.1751608567e+02_wp
        real(wp), parameter     :: c81   = 0.9486883238e+02_wp
        real(wp), parameter     :: c82   = -.1519334221e+02_wp
        real(wp), parameter     :: c83   = 0.4024697252e+04_wp
        real(wp), parameter     :: c84   = -.2225159395e+02_wp
        real(wp), parameter     :: d11   = 0.4203543357e+03_wp
        real(wp), parameter     :: d12   = -.4922474096e+02_wp
        real(wp), parameter     :: d13   = 0.3362942544e+00_wp
        real(wp), parameter     :: d14   = -.3827423082e+03_wp
        real(wp), parameter     :: d15   = 0.1746726001e+03_wp
        real(wp), parameter     :: d21   = 0.1699995737e-01_wp
        real(wp), parameter     :: d22   = 0.1513036778e-01_wp
        real(wp), parameter     :: d23   = 0.2659119354e-01_wp
        real(wp), parameter     :: d24   = -.5760387483e-02_wp
        real(wp), parameter     :: d31   = 0.1020622621e+02_wp
        real(wp), parameter     :: d32   = 0.1050536271e-01_wp
        real(wp), parameter     :: d41   = 0.6836172780e+00_wp
        real(wp), parameter     :: d42   = -.1627858240e+00_wp
        real(wp), parameter     :: d43   = -.6925485045e+01_wp
        real(wp), parameter     :: d44   = 0.1632567385e+01_wp
        real(wp), parameter     :: d51   = 0.1083595009e+04_wp
        real(wp), parameter     :: d52   = 0.4641431791e+01_wp
        real(wp), parameter     :: d53   = -.8233144461e+00_wp
        real(wp), parameter     :: d61   = -.6157225942e+02_wp
        real(wp), parameter     :: d62   = 0.3094361471e+03_wp
        real(wp), parameter     :: d63   = -.3299631143e+03_wp
        real(wp), parameter     :: d71   = 0.8866227120e+01_wp
        real(wp), parameter     :: d72   = -.1382126854e+01_wp
        real(wp), parameter     :: d73   = 0.7620770145e-01_wp
        real(wp), parameter     :: d74   = -.5145757859e+02_wp
        real(wp), parameter     :: d75   = 0.2046097265e+01_wp
        real(wp), parameter     :: d81   = 0.2540775558e+01_wp
        real(wp), parameter     :: d82   = -.4889246569e+00_wp
        real(wp), parameter     :: d83   = -.1127439280e+04_wp
        real(wp), parameter     :: d84   = -.2269932295e+01_wp
        cx1    = c51 + c83
        dx1    = d51 + d83
        r1     = rpass(1)
        r2     = rpass(2)
        r3     = rpass(3)
        t1     = r1*r1 -r2*r2 -r3*r3
        t2     = r2*r2 -r3*r3 -r1*r1
        t3     = r3*r3 -r1*r1 -r2*r2
        c1     = t1 / (-2.e0_wp*r2*r3)
        c2     = t2 / (-2.e0_wp*r3*r1)
        c3     = t3 / (-2.e0_wp*r1*r2)
        sum    = c1 + c2 + c3
        b1a    = 1.e0_wp - sum
        c1cube = c1*c1*c1
        c2cube = c2*c2*c2
        c3cube = c3*c3*c3
        cos3t1 = 4.e0_wp*c1cube - 3.e0_wp*c1
        cos3t2 = 4.e0_wp*c2cube - 3.e0_wp*c2
        cos3t3 = 4.e0_wp*c3cube - 3.e0_wp*c3
        sumb   = cos3t1 + cos3t2 + cos3t3
        b1b    = 1.e0_wp - ( z58*sumb + z38*sum )
        if( ideriv .gt. 0 )then
            dc1dr1  = -r1/(r2*r3)
            dc2dr2  = -r2/(r1*r3)
            dc3dr3  = -r3/(r1*r2)
            dc1dr2  = ( t1/(r2*r2) + 2.e0_wp )/(2.e0_wp*r3)
            dc1dr3  = ( t1/(r3*r3) + 2.e0_wp )/(2.e0_wp*r2)
            dc2dr1  = ( t2/(r1*r1) + 2.e0_wp )/(2.e0_wp*r3)
            dc2dr3  = ( t2/(r3*r3) + 2.e0_wp )/(2.e0_wp*r1)
            dc3dr1  = ( t3/(r1*r1) + 2.e0_wp )/(2.e0_wp*r2)
            dc3dr2  = ( t3/(r2*r2) + 2.e0_wp )/(2.e0_wp*r1)
            db1a(1) = -1.e0_wp*( dc1dr1 + dc2dr1 + dc3dr1 )
            db1a(2) = -1.e0_wp*( dc1dr2 + dc2dr2 + dc3dr2 )
            db1a(3) = -1.e0_wp*( dc1dr3 + dc2dr3 + dc3dr3 )
            d1      = 12.e0_wp*c1*c1 - 3.e0_wp
            d2      = 12.e0_wp*c2*c2 - 3.e0_wp
            d3      = 12.e0_wp*c3*c3 - 3.e0_wp
            db1b(1) = -z58*( d1*dc1dr1 + d2*dc2dr1 + d3*dc3dr1 ) -z38*(    dc1dr1 +    dc2dr1 +    dc3dr1 )
            db1b(2) = -z58*( d1*dc1dr2 + d2*dc2dr2 + d3*dc3dr2 ) -z38*(    dc1dr2 +    dc2dr2 +    dc3dr2 )
            db1b(3) = -z58*( d1*dc1dr3 + d2*dc2dr3 + d3*dc3dr3 ) -z38*(    dc1dr3 +    dc2dr3 +    dc3dr3 )
        endif
        r       = r1 + r2 + r3
        rsq     = r*r
        b2      = 1.e0_wp/r1 + 1.e0_wp/r2 + 1.e0_wp/r3
        b3      = (r2-r1)*(r2-r1) + (r3-r2)*(r3-r2) + (r1-r3)*(r1-r3)
        eps2    = 1.e-12_wp
        b3b     = sqrt( b3 + eps2 )
        exp1    = exp(-beta1*r)
        exp2    = exp(-beta2*rsq)
        exp7    = exp(-beta3*r)
        dexp1   = -beta1*exp1
        dexp2   = -2.e0_wp*beta2*r*exp2
        dexp7   = -beta3*exp7
        b1      = b1a
        b12     = b1*b1
        b13     = b12*b1
        b14     = b13*b1
        b15     = b14*b1
        asum    = a11 +  a12*r +  a13*rsq
        bsum    = a21*b12 +  a22*b13 +  a23*b14 +  a24*b15
        csum    = a31*b1*exp1 +  a32*b12*exp2
        dsum1   = a41*exp1 +  a42*exp2
        dsum2   = a43*exp1 +  a44*exp2
        fsum    = a51 +  a52*r +  a53*rsq
        vb(1,1) = b1*asum*exp1
        vb(1,2) = bsum * exp2
        vb(1,3) = b2*csum
        vb(1,4) = b1*b3*dsum1 + b1*b3b*dsum2
        vb(1,5) = b1* fsum  *exp7
        vbnd(1) = vb(1,1)+vb(1,2)+vb(1,3)+vb(1,4)+vb(1,5)
        if( ideriv .gt. 0 )then
            db2(1)  = -1.e0_wp/(r1*r1)
            db2(2)  = -1.e0_wp/(r2*r2)
            db2(3)  = -1.e0_wp/(r3*r3)
            db3(1)  = 4.e0_wp*r1 - 2.e0_wp*r2 - 2.e0_wp*r3
            db3(2)  = 4.e0_wp*r2 - 2.e0_wp*r3 - 2.e0_wp*r1
            db3(3)  = 4.e0_wp*r3 - 2.e0_wp*r1 - 2.e0_wp*r2
            db3(4)  = 0.5e0_wp*db3(1)/b3b
            db3(5)  = 0.5e0_wp*db3(2)/b3b
            db3(6)  = 0.5e0_wp*db3(3)/b3b
            dasum   = a12 + 2.e0_wp*a13*r
            dbsum   = 2.e0_wp*a21*b1 + 3.e0_wp*a22*b12 + 4.e0_wp*a23*b13 + 5.e0_wp*a24*b14
            ddsum1  = a41*dexp1 + a42*dexp2
            ddsum2  = a43*dexp1 + a44*dexp2
            dfsum   = a52 + 2.e0_wp*a53*r
            dvb1(1) = db1a(1)*asum*exp1 + b1*dasum*exp1 + b1*asum*dexp1
            dvb1(2) = db1a(2)*asum*exp1 + b1*dasum*exp1 + b1*asum*dexp1
            dvb1(3) = db1a(3)*asum*exp1 + b1*dasum*exp1 + b1*asum*dexp1
            dvb2(1) = dbsum*db1a(1)*exp2 + bsum*dexp2
            dvb2(2) = dbsum*db1a(2)*exp2 + bsum*dexp2
            dvb2(3) = dbsum*db1a(3)*exp2 + bsum*dexp2
            dvb3(1) = db2(1)*csum + b2*(a31*db1a(1)*exp1 + a31*b1*dexp1 + 2.e0_wp*a32*b1*db1a(1)*exp2 + a32*b12*dexp2)
            dvb3(2) = db2(2)*csum + b2*(a31*db1a(2)*exp1 + a31*b1*dexp1 + 2.e0_wp*a32*b1*db1a(2)*exp2 + a32*b12*dexp2)
            dvb3(3) = db2(3)*csum + b2*(a31*db1a(3)*exp1 + a31*b1*dexp1 + 2.e0_wp*a32*b1*db1a(3)*exp2 + a32*b12*dexp2)
            dvb4(1) = db1a(1)*b3*dsum1 + b1*db3(1)*dsum1 + b1*b3*ddsum1 + db1a(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1* &
                & b3b*ddsum2
            dvb4(2) = db1a(2)*b3*dsum1 + b1*db3(2)*dsum1 + b1*b3*ddsum1 + db1a(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1* &
                & b3b*ddsum2
            dvb4(3) = db1a(3)*b3*dsum1 + b1*db3(3)*dsum1 + b1*b3*ddsum1 + db1a(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1* &
                & b3b*ddsum2
            dvb5(1)   = db1a(1)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
            dvb5(2)   = db1a(2)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
            dvb5(3)   = db1a(3)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
            dvbnda(1) = dvb1(1) + dvb2(1) + dvb3(1) + dvb4(1) + dvb5(1)
            dvbnda(2) = dvb1(2) + dvb2(2) + dvb3(2) + dvb4(2) + dvb5(2)
            dvbnda(3) = dvb1(3) + dvb2(3) + dvb3(3) + dvb4(3) + dvb5(3)
        endif
        vbnda   = vbnd(1)
        b1      = b1b
        b12     = b1*b1
        b13     = b12*b1
        b14     = b13*b1
        b15     = b14*b1
        asum    = g11 +  g12*r +  g13*rsq
        bsum    = g21*b12 +  g22*b13 +  g23*b14 +  g24*b15
        csum    = g31*b1*exp1 +  g32*b12*exp2
        dsum1   = g41*exp1 +  g42*exp2
        dsum2   = g43*exp1 +  g44*exp2
        fsum    = g51 +  g52*r +  g53*rsq
        vb(2,1) = b1*asum*exp1
        vb(2,2) = bsum * exp2
        vb(2,3) = b2*csum
        vb(2,4) = b1*b3*dsum1 + b1*b3b*dsum2
        vb(2,5) = b1* fsum  *exp7
        vbnd(2) = vb(2,1)+vb(2,2)+vb(2,3)+vb(2,4)+vb(2,5)
        if( ideriv .gt. 0 )then
            dasum   = g12 + 2.e0_wp*g13*r
            dbsum   = 2.e0_wp*g21*b1 + 3.e0_wp*g22*b12 + 4.e0_wp*g23*b13 + 5.e0_wp*g24*b14
            ddsum1  = g41*dexp1 + g42*dexp2
            ddsum2  = g43*dexp1 + g44*dexp2
            dfsum   = g52 + 2.e0_wp*g53*r
            dvb1(1) = db1b(1)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
            dvb1(2) = db1b(2)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
            dvb1(3) = db1b(3)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
            dvb2(1) = dbsum*db1b(1)*exp2 + bsum*dexp2
            dvb2(2) = dbsum*db1b(2)*exp2 + bsum*dexp2
            dvb2(3) = dbsum*db1b(3)*exp2 + bsum*dexp2
            dvb3(1) = db2(1)*csum + b2*(  g31*db1b(1)*exp1 +  g31*b1*dexp1 +2.e0_wp* g32*b1*db1b(1)*exp2 +  g32*b12* &
                & dexp2)
            dvb3(2) = db2(2)*csum + b2*(  g31*db1b(2)*exp1 +  g31*b1*dexp1 +2.e0_wp* g32*b1*db1b(2)*exp2 +  g32*b12* &
                & dexp2)
            dvb3(3) = db2(3)*csum + b2*(  g31*db1b(3)*exp1 +  g31*b1*dexp1 +2.e0_wp* g32*b1*db1b(3)*exp2 +  g32*b12* &
                & dexp2)
            dvb4(1) = db1b(1)*b3 *dsum1+  b1*db3(1)*dsum1 + b1*b3 *ddsum1 + db1b(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1* &
                & b3b*ddsum2
            dvb4(2) = db1b(2)*b3 *dsum1+  b1*db3(2)*dsum1 + b1*b3 *ddsum1 + db1b(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1* &
                & b3b*ddsum2
            dvb4(3)= db1b(3)*b3*dsum1 + b1*db3(3)*dsum1 + b1*b3*ddsum1 + db1b(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1*b3b* &
                & ddsum2
            dvb5(1)   = db1b(1)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
            dvb5(2)   = db1b(2)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
            dvb5(3)   = db1b(3)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
            dvbndb(1) = dvb1(1)+dvb2(1)+dvb3(1)+dvb4(1)+dvb5(1)
            dvbndb(2) = dvb1(2)+dvb2(2)+dvb3(2)+dvb4(2)+dvb5(2)
            dvbndb(3) = dvb1(3)+dvb2(3)+dvb3(3)+dvb4(3)+dvb5(3)
        endif
        vbndb = vbnd(2)
        if( icompc .eq. 0 ) return
        sumt      = t(1) + t(2) + t(3)
        rcu       = rsq*r
        p         = r1*r2*r3
        psq       = p*p
        pcu       = psq*p
        dp(1)     = r2*r3
        dp(2)     = r3*r1
        dp(3)     = r1*r2
        cbnds(1)  = 0.e0_wp
        cbnds(2)  = 0.e0_wp
        dcbnda(1) = 0.e0_wp
        dcbnda(2) = 0.e0_wp
        dcbnda(3) = 0.e0_wp
        dcbndb(1) = 0.e0_wp
        dcbndb(2) = 0.e0_wp
        dcbndb(3) = 0.e0_wp
        exp7      = exp( -beta2*rcu )
        dexp7     = -3.e0_wp*beta2*rsq*exp7
        b1        = b1a
        b12       = b1*b1
        b13       = b12*b1
        b14       = b13*b1
        b15       = b14*b1
        asum      = c11 + c12*r + c13*rsq + c14/r + c15/rsq
        bsum      = c21*b12 + c22*b13 + c23*b14 + c24*b15
        csum      = c31*b1 *exp1  + c32*b12*exp2
        dsum1     = c41*exp1 + c42*exp2
        dsum2     = c43*exp1 + c44*exp2
        fsum      = cx1 + c52*r + c53*rsq
        gsum      = c61 + c62/r + c63/rsq
        aasum     = c71 + c72*p + c73*psq + c74/p + c75/psq
        ffsum     = c81*p + c82*psq + c84/psq
        dasum     = c12 + 2.e0_wp*c13*r -c14/rsq -2.e0_wp*c15/rcu
        dbsum     = 2.e0_wp*c21*b1 +3.e0_wp*c22*b12 +4.e0_wp*c23*b13 +5.e0_wp*c24*b14
        ddsum1    = c41*dexp1 + c42*dexp2
        ddsum2    = c43*dexp1 + c44*dexp2
        dfsum     = c52 + 2.e0_wp*c53*r
        dgsum     = -c62/rsq - 2.e0_wp*c63/rcu
        daasum    = c72 +2.e0_wp*c73*p -c74/psq -2.e0_wp*c75/pcu
        dffsum    = c81 + 2.e0_wp*c82*p -2.e0_wp*c84/pcu
        cb(1,1)   = b1*asum * exp1  /p
        cb(1,2)   = bsum * exp2
        cb(1,3)   = b2*csum
        cb(1,4)   = b1*b3*dsum1 + b1*b3b*dsum2
        cb(1,5)   = b1 * fsum * exp7 /p
        cb(1,6)   = b1 * gsum * exp7
        cb(1,7)   = b1 * aasum * exp2
        cb(1,8)   = b1 * ffsum * exp7
        cbnds(1)  = cb(1,1)+cb(1,2)+cb(1,3)+cb(1,4)+cb(1,5) +cb(1,6)+cb(1,7)+cb(1,8)
        if( ideriv .gt. 0 )then
            dcb1(1) = db1a(1)*asum*exp1/p + b1*dasum*exp1/p + b1*asum*dexp1/p - b1*asum*exp1*dp(1)/psq
            dcb1(2) = db1a(2)*asum*exp1/p + b1*dasum*exp1/p + b1*asum*dexp1/p - b1*asum*exp1*dp(2)/psq
            dcb1(3) = db1a(3)*asum*exp1/p + b1*dasum*exp1/p + b1*asum*dexp1/p - b1*asum*exp1*dp(3)/psq
            dcb2(1) = dbsum*db1a(1)*exp2 + bsum*dexp2
            dcb2(2) = dbsum*db1a(2)*exp2 + bsum*dexp2
            dcb2(3) = dbsum*db1a(3)*exp2 + bsum*dexp2
            dcb3(1) = db2(1)*csum + b2*(c31*db1a(1)*exp1 + c31*b1*dexp1 + 2.e0_wp*c32*b1*db1a(1)*exp2 + c32*b12*dexp2)
            dcb3(2) = db2(2)*csum + b2*(c31*db1a(2)*exp1 + c31*b1*dexp1 + 2.e0_wp*c32*b1*db1a(2)*exp2 + c32*b12*dexp2)
            dcb3(3) = db2(3)*csum + b2*(c31*db1a(3)*exp1 + c31*b1*dexp1 + 2.e0_wp*c32*b1*db1a(3)*exp2 + c32*b12*dexp2)
            dcb4(1) = db1a(1)*b3*dsum1 + b1*db3(1)*dsum1 + b1*b3*ddsum1 + db1a(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1* &
                & b3b*ddsum2
            dcb4(2) = db1a(2)*b3*dsum1 + b1*db3(2)*dsum1 + b1*b3*ddsum1 + db1a(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1* &
                & b3b*ddsum2
            dcb4(3) = db1a(3)*b3*dsum1 + b1*db3(3)*dsum1 + b1*b3*ddsum1 + db1a(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1* &
                & b3b*ddsum2
            dcb5(1)     = db1a(1)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(1)/psq
            dcb5(2)     = db1a(2)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(2)/psq
            dcb5(3)     = db1a(3)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(3)/psq
            dcb6(1)     = db1a(1)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
            dcb6(2)     = db1a(2)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
            dcb6(3)     = db1a(3)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
            dcb7(1)     = db1a(1)*aasum*exp2 + b1*daasum*dp(1)*exp2 + b1*aasum*dexp2
            dcb7(2)     = db1a(2)*aasum*exp2 + b1*daasum*dp(2)*exp2 + b1*aasum*dexp2
            dcb7(3)     = db1a(3)*aasum*exp2 + b1*daasum*dp(3)*exp2 + b1*aasum*dexp2
            dcb8(1)     = db1a(1)*ffsum*exp7 + b1*dffsum*dp(1)*exp7 + b1*ffsum*dexp7
            dcb8(2)     = db1a(2)*ffsum*exp7 + b1*dffsum*dp(2)*exp7 + b1*ffsum*dexp7
            dcb8(3)     = db1a(3)*ffsum*exp7 + b1*dffsum*dp(3)*exp7 + b1*ffsum*dexp7
            dcbnds(1,1) = dcb1(1) +dcb2(1) +dcb3(1) +dcb4(1) +dcb5(1) +dcb6(1) +dcb7(1) +dcb8(1)
            dcbnds(1,2) = dcb1(2) +dcb2(2) +dcb3(2) +dcb4(2) +dcb5(2) +dcb6(2) +dcb7(2) +dcb8(2)
            dcbnds(1,3) = dcb1(3) +dcb2(3) +dcb3(3) +dcb4(3) +dcb5(3) +dcb6(3) +dcb7(3) +dcb8(3)
        endif
        b1       = b1b
        b12      = b1*b1
        b13      = b12*b1
        b14      = b13*b1
        b15      = b14*b1
        asum     = d11 + d12*r + d13*rsq + d14/r + d15/rsq
        bsum     = d21*b12 + d22*b13 + d23*b14 + d24*b15
        csum     = d31*b1 *exp1  + d32*b12*exp2
        dsum1    = d41*exp1 + d42*exp2
        dsum2    = d43*exp1 + d44*exp2
        fsum     = dx1 + d52*r + d53*rsq
        gsum     = d61 + d62/r + d63/rsq
        aasum    = d71 + d72*p + d73*psq + d74/p + d75/psq
        ffsum    = d81*p + d82*psq + d84/psq
        dasum    = d12 + 2.e0_wp*d13*r -d14/rsq -2.e0_wp*d15/rcu
        dbsum    = 2.e0_wp*d21*b1 +3.e0_wp*d22*b12 +4.e0_wp*d23*b13 +5.e0_wp*d24*b14
        ddsum1   = d41*dexp1 + d42*dexp2
        ddsum2   = d43*dexp1 + d44*dexp2
        dfsum    = d52 + 2.e0_wp*d53*r
        dgsum    = -d62/rsq - 2.e0_wp*d63/rcu
        daasum   = d72 +2.e0_wp*d73*p -d74/psq -2.e0_wp*d75/pcu
        dffsum   = d81 + 2.e0_wp*d82*p -2.e0_wp*d84/pcu
        cb(2,1)  = b1*asum * exp1  /p
        cb(2,2)  = bsum * exp2
        cb(2,3)  = b2*csum
        cb(2,4)  = b1*b3*dsum1 + b1*b3b*dsum2
        cb(2,5)  = b1 * fsum * exp7 /p
        cb(2,6)  = b1 * gsum * exp7
        cb(2,7)  = b1 * aasum * exp2
        cb(2,8)  = b1 * ffsum * exp7
        cbnds(2) = cb(2,1)+cb(2,2)+cb(2,3)+cb(2,4)+cb(2,5) +cb(2,6)+cb(2,7)+cb(2,8)
        if( ideriv .gt. 0 )then
            dcb1(1) = db1b(1)*asum*exp1/p + b1*dasum*exp1/p + b1*asum*dexp1/p - b1*asum*exp1*dp(1)/psq
            dcb1(2) = db1b(2)*asum*exp1/p + b1*dasum*exp1/p + b1*asum*dexp1/p - b1*asum*exp1*dp(2)/psq
            dcb1(3) = db1b(3)*asum*exp1/p + b1*dasum*exp1/p + b1*asum*dexp1/p - b1*asum*exp1*dp(3)/psq
            dcb2(1) = dbsum*db1b(1)*exp2 + bsum*dexp2
            dcb2(2) = dbsum*db1b(2)*exp2 + bsum*dexp2
            dcb2(3) = dbsum*db1b(3)*exp2 + bsum*dexp2
            dcb3(1) = db2(1)*csum + b2*(d31*db1b(1)*exp1 + d31*b1*dexp1 + 2.e0_wp*d32*b1*db1b(1)*exp2 + d32*b12*dexp2)
            dcb3(2) = db2(2)*csum + b2*(d31*db1b(2)*exp1 + d31*b1*dexp1 + 2.e0_wp*d32*b1*db1b(2)*exp2 + d32*b12*dexp2)
            dcb3(3) = db2(3)*csum + b2*(d31*db1b(3)*exp1 + d31*b1*dexp1 + 2.e0_wp*d32*b1*db1b(3)*exp2 + d32*b12*dexp2)
            dcb4(1) = db1b(1)*b3*dsum1 + b1*db3(1)*dsum1 + b1*b3*ddsum1 + db1b(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1* &
                & b3b*ddsum2
            dcb4(2) = db1b(2)*b3*dsum1 + b1*db3(2)*dsum1 + b1*b3*ddsum1 + db1b(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1* &
                & b3b*ddsum2
            dcb4(3) = db1b(3)*b3*dsum1 + b1*db3(3)*dsum1 + b1*b3*ddsum1 + db1b(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1* &
                & b3b*ddsum2
            dcb5(1)     = db1b(1)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(1)/psq
            dcb5(2)     = db1b(2)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(2)/psq
            dcb5(3)     = db1b(3)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(3)/psq
            dcb6(1)     = db1b(1)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
            dcb6(2)     = db1b(2)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
            dcb6(3)     = db1b(3)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
            dcb7(1)     = db1b(1)*aasum*exp2 + b1*daasum*dp(1)*exp2 + b1*aasum*dexp2
            dcb7(2)     = db1b(2)*aasum*exp2 + b1*daasum*dp(2)*exp2 + b1*aasum*dexp2
            dcb7(3)     = db1b(3)*aasum*exp2 + b1*daasum*dp(3)*exp2 + b1*aasum*dexp2
            dcb8(1)     = db1b(1)*ffsum*exp7 + b1*dffsum*dp(1)*exp7 + b1*ffsum*dexp7
            dcb8(2)     = db1b(2)*ffsum*exp7 + b1*dffsum*dp(2)*exp7 + b1*ffsum*dexp7
            dcb8(3)     = db1b(3)*ffsum*exp7 + b1*dffsum*dp(3)*exp7 + b1*ffsum*dexp7
            dcbnds(2,1) = dcb1(1) +dcb2(1) +dcb3(1) +dcb4(1) +dcb5(1) +dcb6(1) +dcb7(1) +dcb8(1)
            dcbnds(2,2) = dcb1(2) +dcb2(2) +dcb3(2) +dcb4(2) +dcb5(2) +dcb6(2) +dcb7(2) +dcb8(2)
            dcbnds(2,3) = dcb1(3) +dcb2(3) +dcb3(3) +dcb4(3) +dcb5(3) +dcb6(3) +dcb7(3) +dcb8(3)
            dcbnda(1)   = dt(1) * cbnds(1)    + sumt * dcbnds(1,1)
            dcbndb(1)   = dt(1) * cbnds(2)    + sumt * dcbnds(2,1)
            dcbnda(2)   = dt(2) * cbnds(1)    + sumt * dcbnds(1,2)
            dcbndb(2)   = dt(2) * cbnds(2)    + sumt * dcbnds(2,2)
            dcbnda(3)   = dt(3) * cbnds(1)    + sumt * dcbnds(1,3)
            dcbndb(3)   = dt(3) * cbnds(2)    + sumt * dcbnds(2,3)
        endif
        cbnda = sumt*cbnds(1)
        cbndb = sumt*cbnds(2)
    end subroutine vbcb95_h4bmkp
end module chempotpy_h4_bmkp_2002_pes

!> @brief Expose H4_BMKP_2002 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h4_bmkp_2002_pes, only: evaluate
    use chempotpy_kinds,            only: wp
    implicit none
    integer, parameter    :: natoms = 4, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
