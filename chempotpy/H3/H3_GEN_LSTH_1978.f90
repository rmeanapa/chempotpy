! System:                        H3
! Functional form:               Liu-Siegbahn-Truhlar-Horowitz
! Common name:                   LSTH
! Reference:                     B. Liu, J. Chem. Phys. 58, 1925 (1973); B. Liu and P. Siegbahn,
!                                J. Chem. Phys. 68, 2457 (1978); D. G. Truhlar and C. J. Horowitz,
!                                J. Chem. Phys. 68, 2466 (1978); 71, 1514(E) (1979).
! Number of derivatives:         1 (centered finite difference)
! Number of bodies:              3
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; the released H2 spline is embedded below
! Notes:                         Ground-state H + H2 surface with atom order H, H, H.

!> @brief Evaluate the LSTH ground-state potential for H + H2.
module chempotpy_h3_gen_lsth_1978_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    real(wp), parameter :: c6        = 6.89992032_wp, c8 = 219.9997304_wp
    real(wp), parameter :: london_c  = -1.2148730613_wp, london_a = -1.514663474_wp
    real(wp), parameter :: london_a1 = -1.46_wp, london_f = 2.088442_wp
    real(wp), parameter :: fns       = 0.0035_wp
    real(wp), parameter :: an1       = 0.0012646477_wp, an2 = -0.0001585792_wp
    real(wp), parameter :: an3       = 0.0000079707_wp, an4 = -0.0000001151_wp
    real(wp), parameter :: f1        = 0.52_wp, b1 = 3.0231771503_wp, b2 = -1.08935219_wp
    real(wp), parameter :: f2        = 0.052_wp, w1 = 1.7732141742_wp
    real(wp), parameter :: w2        = -2.0979468223_wp, w3 = -3.978850217_wp
    real(wp), parameter :: d1        = 0.4908116374_wp, d2 = -0.8718696387_wp
    real(wp), parameter :: d3        = 0.1612118092_wp, d4 = -0.1273731045_wp
    real(wp), parameter :: f3        = 0.79_wp, xl2 = 0.9877930913_wp, xl1 = -13.3599568553_wp
    real(wp), parameter :: rkw(87)   = [ &
        .400000000e+00_wp, .450000000e+00_wp, .500000000e+00_wp, .550000000e+00_wp, &
        .600000000e+00_wp, .650000000e+00_wp, .700000000e+00_wp, .750000000e+00_wp, &
        .800000000e+00_wp, .900000000e+00_wp, .100000000e+01_wp, .110000000e+01_wp, &
        .120000000e+01_wp, .130000000e+01_wp, .135000000e+01_wp, .139000000e+01_wp, &
        .140000000e+01_wp, .140100001e+01_wp, .140109999e+01_wp, .141000000e+01_wp, &
        .145000000e+01_wp, .150000000e+01_wp, .160000000e+01_wp, .170000000e+01_wp, &
        .180000000e+01_wp, .190000000e+01_wp, .200000000e+01_wp, .210000000e+01_wp, &
        .220000000e+01_wp, .230000000e+01_wp, .240000000e+01_wp, .250000000e+01_wp, &
        .260000000e+01_wp, .270000000e+01_wp, .280000000e+01_wp, .290000000e+01_wp, &
        .300000000e+01_wp, .310000000e+01_wp, .320000000e+01_wp, .330000000e+01_wp, &
        .340000000e+01_wp, .350000000e+01_wp, .360000000e+01_wp, .370000000e+01_wp, &
        .380000000e+01_wp, .390000000e+01_wp, .400000000e+01_wp, .410000000e+01_wp, &
        .420000000e+01_wp, .430000000e+01_wp, .440000000e+01_wp, .450000000e+01_wp, &
        .460000000e+01_wp, .470000000e+01_wp, .480000000e+01_wp, .490000000e+01_wp, &
        .500000000e+01_wp, .510000000e+01_wp, .520000000e+01_wp, .530000000e+01_wp, &
        .540000000e+01_wp, .550000000e+01_wp, .560000000e+01_wp, .570000000e+01_wp, &
        .580000000e+01_wp, .590000000e+01_wp, .600000000e+01_wp, .610000000e+01_wp, &
        .620000000e+01_wp, .630000000e+01_wp, .640000000e+01_wp, .650000000e+01_wp, &
        .660000000e+01_wp, .670000000e+01_wp, .680000000e+01_wp, .690000000e+01_wp, &
        .700000000e+01_wp, .720000000e+01_wp, .740000000e+01_wp, .760000000e+01_wp, &
        .780000000e+01_wp, .800000000e+01_wp, .824999991e+01_wp, .850000000e+01_wp, &
        .900000000e+01_wp, .950000000e+01_wp, .100000000e+02_wp &
        ]
    real(wp), parameter :: ekw(87) = [ &
        .879796188e+00_wp, .649071056e+00_wp, .473372447e+00_wp, .337228924e+00_wp, &
        .230365628e+00_wp, .145638432e+00_wp, .779738117e-01_wp, .236642733e-01_wp, &
        -.200555771e-01_wp, -.836421044e-01_wp, -.124538356e+00_wp, -.150056027e+00_wp, &
        -.164934012e+00_wp, -.172345701e+00_wp, -.173962500e+00_wp, -.174451499e+00_wp, &
        -.174474200e+00_wp, -.174474400e+00_wp, -.174474400e+00_wp, -.174459699e+00_wp, &
        -.174055600e+00_wp, -.172853502e+00_wp, -.168579707e+00_wp, -.162456813e+00_wp, &
        -.155066822e+00_wp, -.146849432e+00_wp, -.138131041e+00_wp, -.129156051e+00_wp, &
        -.120123163e+00_wp, -.111172372e+00_wp, -.102412583e+00_wp, -.939271927e-01_wp, &
        -.857809026e-01_wp, -.780163108e-01_wp, -.706699181e-01_wp, -.637640270e-01_wp, &
        -.573117349e-01_wp, -.513184414e-01_wp, -.457831464e-01_wp, -.407002530e-01_wp, &
        -.360577581e-01_wp, -.318401624e-01_wp, -.280271683e-01_wp, -.245977718e-01_wp, &
        -.215296753e-01_wp, -.187966785e-01_wp, -.163688812e-01_wp, -.142246837e-01_wp, &
        -.123370858e-01_wp, -.106809878e-01_wp, -.923028934e-02_wp, -.796819096e-02_wp, &
        -.687029215e-02_wp, -.591779314e-02_wp, -.509229414e-02_wp, -.437819496e-02_wp, &
        -.376259562e-02_wp, -.323089623e-02_wp, -.277399691e-02_wp, -.237999732e-02_wp, &
        -.204229767e-02_wp, -.175209799e-02_wp, -.150299828e-02_wp, -.128989853e-02_wp, &
        -.110689874e-02_wp, -.949798920e-03_wp, -.814999069e-03_wp, -.700199190e-03_wp, &
        -.602999302e-03_wp, -.516199400e-03_wp, -.446599479e-03_wp, -.386399548e-03_wp, &
        -.332799617e-03_wp, -.290599668e-03_wp, -.246599722e-03_wp, -.215399753e-03_wp, &
        -.188899784e-03_wp, -.143399836e-03_wp, -.108599875e-03_wp, -.867998994e-04_wp, &
        -.681999214e-04_wp, -.527999393e-04_wp, -.403999540e-04_wp, -.313999636e-04_wp, &
        -.184999787e-04_wp, -.120999861e-04_wp, -.909998949e-05_wp &
        ]
    real(wp), parameter :: wkw(87) = [ &
        .308019605e+02_wp, .214419954e+02_wp, .154937452e+02_wp, .115151545e+02_wp, &
        .871827707e+01_wp, .673831756e+01_wp, .527864661e+01_wp, .419929947e+01_wp, &
        .333940643e+01_wp, .219403463e+01_wp, .149861953e+01_wp, .103863661e+01_wp, &
        .730647471e+00_wp, .518552387e+00_wp, .441110777e+00_wp, .383461006e+00_wp, &
        .373946396e+00_wp, .358559402e+00_wp, .372215569e+00_wp, .356670198e+00_wp, &
        .312744133e+00_wp, .261523038e+00_wp, .180817537e+00_wp, .124665543e+00_wp, &
        .807794104e-01_wp, .486562494e-01_wp, .251952492e-01_wp, .452257820e-02_wp, &
        -.854560161e-02_wp, -.196001146e-01_wp, -.276538076e-01_wp, -.344244662e-01_wp, &
        -.381080935e-01_wp, -.421628973e-01_wp, -.441600287e-01_wp, -.454966841e-01_wp, &
        -.460129217e-01_wp, -.458513118e-01_wp, -.453815149e-01_wp, -.440623159e-01_wp, &
        -.426089183e-01_wp, -.404417185e-01_wp, -.383839285e-01_wp, -.361823035e-01_wp, &
        -.336666088e-01_wp, -.302110314e-01_wp, -.286090554e-01_wp, -.255125522e-01_wp, &
        -.233005599e-01_wp, -.201850499e-01_wp, -.191990995e-01_wp, -.161784216e-01_wp, &
        -.146071006e-01_wp, -.126330766e-01_wp, -.110605069e-01_wp, -.996481997e-02_wp, &
        -.818014482e-02_wp, -.765454189e-02_wp, -.608163613e-02_wp, -.575887028e-02_wp, &
        -.466284400e-02_wp, -.408972107e-02_wp, -.363824334e-02_wp, -.295728079e-02_wp, &
        -.259261281e-02_wp, -.221225014e-02_wp, -.193837141e-02_wp, -.203425060e-02_wp, &
        -.484614204e-03_wp, -.226728547e-02_wp, -.766232140e-03_wp, -.307779418e-03_wp, &
        -.196264565e-02_wp, .131836977e-02_wp, -.223083472e-02_wp, -.750220030e-04_wp, &
        -.289074004e-03_wp, -.220265690e-03_wp, -.434861384e-03_wp, .971346041e-05_wp, &
        -.839919101e-04_wp, -.153745275e-03_wp, -.369227366e-04_wp, -.249634065e-04_wp, &
        -.290482724e-04_wp, -.148433244e-04_wp, .682166282e-05_wp &
        ]

    contains

    !> @brief Return the LSTH energy and optional numerical Cartesian gradient.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_finite_difference, only: central_finite_difference_gradient
        use chempotpy_status,            only: chempotpy_status_t, status_derivative_unavailable
        integer, parameter    :: natoms = 3, nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        type(chempotpy_status_t), intent(out), optional :: status
        type(chempotpy_status_t) :: evaluation_status
        real(wp), parameter :: step = 1.0e-5_wp
        real(wp)            :: flat_gradient(1,9)
        integer             :: atom, direction
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        call evaluation_status%clear()
        if( igrad == 0 )then
            call evaluate_energy_ev(x, p, evaluation_status)
        elseif( igrad == 1 )then
            call central_finite_difference_gradient(evaluate_energy_ev, x, step, p, flat_gradient, evaluation_status)
            if( .not. evaluation_status%failed() )then
                do atom=1, natoms
                    do direction=1, 3
                        g(1,atom,direction) = flat_gradient(1,3*(atom-1)+direction)
                    enddo
                enddo
            endif
        else
            call evaluation_status%set(status_derivative_unavailable, &
                & 'H3_GEN_LSTH_1978 provides energy and a numerical first derivative')
        endif
        if( present(status) ) status = evaluation_status
    end subroutine evaluate

    !> @brief Evaluate the LSTH energy in electronvolt for one Cartesian geometry.
    pure subroutine evaluate_energy_ev( coordinates, values, status )
        use chempotpy_internal_coordinates, only: prepare_legacy_triatomic_coordinates
        use chempotpy_status,               only: chempotpy_status_t
        use chempotpy_units,                only: hartree_to_electronvolt
        real(wp), intent(in)  :: coordinates(:,:)
        real(wp), intent(out) :: values(:)
        type(chempotpy_status_t), intent(inout) :: status
        real(wp) :: distance_jacobian(3,9), distances(3)
        if( status%failed() ) return
        call prepare_legacy_triatomic_coordinates(coordinates, distances, distance_jacobian)
        values(1) = lsth_energy_hartree(distances)*hartree_to_electronvolt
    end subroutine evaluate_energy_ev

    !> @brief Evaluate the released LSTH internal-distance energy in hartree.
    pure function lsth_energy_hartree( x ) result( energy )
        real(wp), intent(in) :: x(3)
        real(wp)             :: energy, s(3,3), t(3), xq_terms(3), xj_terms(3)
        real(wp)             :: x2(3), xq, xj, elond, radius, radius2, radius3
        real(wp)             :: exns, wnt, wn, wb, wb2, wb3, wb4, cos_sum, cosines(3)
        real(wp)             :: exf1, exf2, exf3, eb1t, eb3t, eb1, eb2, eq, reciprocal_sum
        real(wp)             :: eb4a, eb4b, eb4, ens, minimum_energy, minimum_radius
        integer              :: i, iteration
        minimum_radius = 1.40105_wp
        do iteration=1, 30
            call spline_evaluate(minimum_radius, s(:,1))
            if( abs(s(2,1)) < 1.0e-10_wp ) exit
            minimum_radius = minimum_radius-s(2,1)/s(3,1)
        enddo
        minimum_energy = s(1,1)
        do i=1, 3
            t(i) = london_c*(london_a+x(i)+london_a1*x(i)*x(i))*exp(-london_f*x(i))
            call h2_curve(x(i), s(:,i))
            xq_terms(i) = s(1,i)+t(i)
            xj_terms(i) = s(1,i)-t(i)
        enddo
        xq     = sum(xq_terms)/2.0_wp
        xj     = sqrt(((xj_terms(1)-xj_terms(2))**2+(xj_terms(2)-xj_terms(3))**2+ (xj_terms(3)-xj_terms(1))**2)/8.0_wp)
        elond  = xq-xj
        radius = sum(x); radius2 = radius*radius; radius3 = radius2*radius
        exns   = exp(-fns*radius3); ens = 0.0_wp
        if( exns /= 0.0_wp )then
            wnt = (x(1)-x(2))*(x(2)-x(3))*(x(3)-x(1)); wn = abs(wnt)
            ens = (an1*wn**2+an2*wn**3+an3*wn**4+an4*wn**5)*exns
        endif
        x2             = x*x; cos_sum = sum(x2)/2.0_wp
        cosines(1)     = (x2(1)-cos_sum)/(x(2)*x(3))
        cosines(2)     = (x2(2)-cos_sum)/(x(1)*x(3))
        cosines(3)     = (x2(3)-cos_sum)/(x(1)*x(2))
        wb             = sum(cosines)+1.0_wp; wb2 = wb*wb; wb3 = wb2*wb; wb4 = wb3*wb
        exf1           = exp(-f1*radius); exf2 = exp(-f2*radius2); exf3 = exp(-f3*radius)
        eb1t           = (b1+b2*radius)*exf1; eb3t = (xl1+xl2*radius2)*exf3
        eb1            = wb*(eb1t+eb3t); eb2 = (wb2*w1+wb3*w2+wb4*w3)*exf2
        eq             = (x(1)-x(2))**2+(x(2)-x(3))**2+(x(3)-x(1))**2
        reciprocal_sum = sum(1.0_wp/x)
        eb4a           = wb*d1*exf1+wb2*d2*exf2; eb4b = d3*exf1+d4*exf2
        eb4            = eb4a*reciprocal_sum+eb4b*wb*eq
        energy         = elond+ens+eb1+eb2+eb4-minimum_energy
    end function lsth_energy_hartree

    !> @brief Evaluate the embedded H2 curve and its first two derivatives.
    pure subroutine h2_curve( distance, values )
        real(wp), intent(in)  :: distance
        real(wp), intent(out) :: values(3)
        if( distance > 10.0_wp )then
            values(1) = -(c6+c8/(distance*distance))/distance**6
            values(2) = (6.0_wp*c6+8.0_wp*c8/(distance*distance))/distance**7
            values(3) = 0.0_wp
        else
            call spline_evaluate(distance, values)
        endif
    end subroutine h2_curve

    !> @brief Evaluate the released cubic spline and its first two derivatives.
    pure subroutine spline_evaluate( coordinate, values )
        real(wp), intent(in)  :: coordinate
        real(wp), intent(out) :: values(3)
        real(wp)              :: left_weight, right_weight, interval
        integer               :: i, k
        if( coordinate <= rkw(1) )then
            i = 1
        elseif( coordinate >= rkw(87) )then
            i = 86
        else
            i = 1
            do k=2, 87
                if( rkw(k) > coordinate )then
                    i = k-1
                    exit
                endif
            enddo
        endif
        interval     = rkw(i+1)-rkw(i)
        left_weight  = rkw(i+1)-coordinate
        right_weight = coordinate-rkw(i)
        values(1) = (wkw(i)*left_weight**3+wkw(i+1)*right_weight**3)/(6.0_wp*interval) + (ekw(i+1)/interval-wkw(i+1)* &
            & interval/6.0_wp)*right_weight + (ekw(i)/interval-interval*wkw(i)/6.0_wp)*left_weight
        values(2) = (wkw(i+1)*right_weight**2-wkw(i)*left_weight**2)/(2.0_wp*interval) + (ekw(i+1)-ekw(i))/interval+ &
            & interval*(wkw(i)-wkw(i+1))/6.0_wp
        values(3) = (wkw(i)*left_weight+wkw(i+1)*right_weight)/interval
    end subroutine spline_evaluate
end module chempotpy_h3_gen_lsth_1978_pes

!> @brief Expose H3_GEN_LSTH_1978 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_lsth_1978_pes, only: evaluate
    use chempotpy_kinds,                only: wp
    implicit none
    integer, parameter    :: natoms = 3, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
