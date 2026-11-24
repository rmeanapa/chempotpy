! System:                        N4
! Functional form:               many-body permutationally invariant polynomial with MEG variables
! Common name:                   N4 MEG2
! Reference:                     J. Li et al., J. Chem. Theory Comput. 16, 4822 (2020)
! Number of derivatives:         1
! Number of bodies:              4
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; 276 fit coefficients and pair parameters are embedded
! Notes:                         The released MEG2 refit uses the unified N2 pair potential.

module n4_1a_meg2_2020_par
    ! This code is based on n4pes-gpip-meg.f (FORTRAN77 version) where some
    ! keywords were repleced by more modern ones.
    ! The surface itself was not changed.
    ! Conversion factors
    ! Cconv: bohr to Angstrom
    !        1 bohr = 0.52917721092 angstrom
    ! Econv: kcal/mol to hartree
    !        1 kcal/mol = 0.159360144 * 10^-2 hartree
    ! Gconv: kcal/(mol*Angstrom) to hartree/bohr
    !        1 kcal mol^-1 angstrom^-1 = 0.843297564 * 10^-3 hartree/bohr
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: a, ab, b, c, cconv, dbdr, dmdr, dmsdr, dpdr, drdx, dvdr, econv, eref, gconv, p, r, re, &
        & rm, rms, totdiss
    real(wp), parameter :: cconv = 0.52917721092e0_wp
    real(wp), parameter :: econv = 0.159360144e-2_wp
    real(wp), parameter :: gconv = 0.843297564e-3_wp
    ! Common variables
    ! R(6):         Interatomic bond distance
    ! rMs(6):       Array to store base terms
    ! rM(0:111):    Array to store monomials
    ! P(0:305):     Array to store polynomials
    ! B(1:276):     Array to store basis functions
    ! dMsdR(6,6):   The derivative of base terms w.r.t. R
    ! dMdR(6,112):  The derivative of monomials w.r.t. R
    ! dPdR(6,306):  The derivative of basis functions w.r.t. R
    ! dVdR(6):      The derivative of V w.r.t. R
    ! dRdX(6,12):   The derivative of R w.r.t. X
    ! dBdR(6,276)   The derivative of B w.r.t. R
    real(wp) :: r(6)
    real(wp) :: rms(6), rm(0:111), p(0:305), b(276)
    real(wp) :: dmsdr(6,6), dmdr(6,0:111), dpdr(6,0:305)
    real(wp) :: dvdr(6), drdx(6,12), dbdr(6,276)
    ! Nonlinear parameters:
    ! a(in Ang)
    ! ab (in Ang^2)
    ! re (in Ang)
    real(wp), parameter :: a  = 1.0e0_wp
    real(wp), parameter :: ab = 1.5e0_wp
    real(wp), parameter :: re = 1.098e0_wp
    ! Reference energy of infinitely separated N2 + N2 in hartree
    real(wp), parameter :: eref = &
        -218.40801323e0_wp
    ! For N2 + N2 framework total diss. energy is 2De 2*228.4 kcal/mol
    real(wp), parameter :: totdiss = 457.4e0_wp
    ! Linear parameters optimized by the weighted-least square fitting
    real(wp), parameter :: c(276) = [ &
        0.200045033574e+03_wp, -0.339587904853e+03_wp, 0.118281385756e+04_wp, -0.384724185628e+04_wp, &
        0.319289278735e+03_wp, -0.115519723268e+05_wp, 0.252027171567e+04_wp, -0.289208359130e+04_wp, &
        0.349443534212e+04_wp, -0.198436047133e+04_wp, 0.687053262948e+04_wp, 0.753563625768e+04_wp, &
        -0.545260141586e+04_wp, -0.114630276188e+05_wp, 0.236563605925e+05_wp, -0.119192368852e+05_wp, &
        -0.912576972886e+04_wp, -0.229169453850e+04_wp, 0.570788086262e+04_wp, -0.945021515234e+04_wp, &
        0.337994426027e+05_wp, -0.631904260249e+04_wp, 0.153707864612e+05_wp, 0.213441251389e+04_wp, &
        -0.273771704033e+05_wp, -0.797540663873e+04_wp, -0.138428253765e+05_wp, 0.170856161290e+05_wp, &
        0.252881781428e+05_wp, 0.492377722884e+04_wp, 0.222301272797e+05_wp, -0.214551231633e+05_wp, &
        0.133890649506e+05_wp, -0.160799276021e+05_wp, -0.230501473789e+05_wp, 0.436369395726e+04_wp, &
        0.929665682096e+04_wp, -0.325460562103e+05_wp, 0.540428295327e+05_wp, 0.129689259134e+05_wp, &
        0.124152282316e+05_wp, 0.650532720470e+04_wp, 0.367541239704e+05_wp, 0.691794919054e+04_wp, &
        -0.111691063525e+05_wp, 0.358989802803e+04_wp, 0.143766721085e+05_wp, -0.510379825278e+05_wp, &
        -0.130755495968e+04_wp, -0.296853652857e+05_wp, -0.149540023081e+05_wp, 0.664896590914e+05_wp, &
        0.984963447172e+04_wp, 0.397977493345e+05_wp, 0.685598820256e+04_wp, -0.239542285316e+05_wp, &
        -0.238883260009e+05_wp, 0.369027824454e+05_wp, 0.146405783403e+05_wp, -0.120639356045e+05_wp, &
        -0.212764090382e+05_wp, -0.996760085438e+04_wp, -0.114203668448e+05_wp, 0.860487170791e+03_wp, &
        0.179521990451e+05_wp, -0.248280953228e+05_wp, 0.128125732573e+05_wp, 0.146892384576e+05_wp, &
        -0.310561860468e+05_wp, -0.300235800681e+05_wp, 0.118752525706e+05_wp, 0.256220118404e+05_wp, &
        0.166614401550e+05_wp, -0.428097738687e+04_wp, -0.968743662259e+04_wp, -0.262825854360e+04_wp, &
        0.209936422565e+05_wp, -0.686431237995e+05_wp, 0.466643315622e+04_wp, -0.144423013899e+05_wp, &
        0.912682289747e+04_wp, 0.673295123156e+04_wp, -0.253958354980e+05_wp, -0.124444752729e+05_wp, &
        -0.361575451759e+04_wp, -0.199420375442e+05_wp, 0.595889037797e+05_wp, 0.164930100978e+04_wp, &
        0.104744114739e+05_wp, 0.392551298439e+04_wp, 0.458355853446e+04_wp, 0.233986502329e+05_wp, &
        0.435411612641e+04_wp, 0.289173081157e+05_wp, 0.112615108708e+05_wp, -0.544574890050e+05_wp, &
        -0.692008718449e+04_wp, -0.297316156508e+05_wp, 0.376765466370e+03_wp, 0.166728805447e+05_wp, &
        0.263521218630e+06_wp, -0.700952905121e+05_wp, 0.432573346813e+05_wp, -0.209541622096e+05_wp, &
        0.296541834142e+05_wp, -0.387836993321e+04_wp, 0.206649206972e+05_wp, -0.449279701024e+05_wp, &
        0.237081890198e+04_wp, 0.153525123905e+05_wp, 0.124677379142e+05_wp, 0.994696583885e+04_wp, &
        0.130309257880e+05_wp, -0.106975335877e+05_wp, -0.312119308966e+04_wp, 0.112687719797e+05_wp, &
        0.105122140236e+05_wp, 0.184728205663e+05_wp, 0.935043882376e+04_wp, -0.684240183902e+04_wp, &
        -0.838802788225e+04_wp, 0.105939545947e+05_wp, -0.242397105178e+04_wp, -0.802728240144e+04_wp, &
        0.689599683244e+04_wp, 0.256851073635e+04_wp, -0.204914317179e+05_wp, 0.371343313984e+04_wp, &
        -0.200053798516e+05_wp, 0.627015233273e+03_wp, 0.506980614817e+04_wp, -0.120770736358e+05_wp, &
        0.210950295548e+05_wp, 0.214398871094e+05_wp, -0.134638828136e+05_wp, -0.735553099954e+04_wp, &
        -0.634045384878e+04_wp, -0.189914284611e+04_wp, 0.890459213525e+04_wp, 0.426629792878e+04_wp, &
        -0.736114543350e+03_wp, 0.341463139053e+05_wp, -0.141453866024e+05_wp, 0.765679590097e+04_wp, &
        -0.625611915152e+04_wp, 0.756120993705e+04_wp, 0.168630249290e+05_wp, 0.560152312341e+04_wp, &
        -0.433667515287e+04_wp, -0.714554027895e+04_wp, 0.121970702972e+05_wp, -0.167351518220e+05_wp, &
        0.317753037485e+04_wp, 0.986911427341e+04_wp, 0.839370989093e+04_wp, -0.265555241547e+05_wp, &
        -0.478926268856e+04_wp, -0.454264288636e+04_wp, -0.349253518537e+04_wp, -0.108097959996e+05_wp, &
        0.437224905573e+04_wp, -0.839827549948e+02_wp, -0.141822952576e+05_wp, -0.694085197424e+02_wp, &
        0.145143436119e+05_wp, 0.814238046137e+02_wp, 0.136289737004e+05_wp, 0.116193607256e+05_wp, &
        -0.130457227675e+04_wp, -0.524171145511e+04_wp, 0.244027092661e+05_wp, -0.370798015239e+05_wp, &
        -0.585765270588e+05_wp, -0.858818617973e+05_wp, 0.252580104931e+05_wp, -0.117094851663e+05_wp, &
        0.573346366430e+04_wp, 0.193944188486e+04_wp, 0.208539287977e+05_wp, 0.915612220009e+04_wp, &
        0.117836454931e+04_wp, -0.582869284469e+04_wp, -0.948498016240e+04_wp, 0.278406842692e+04_wp, &
        -0.162758137182e+04_wp, 0.994561527718e+04_wp, 0.802259126226e+04_wp, -0.204178391709e+05_wp, &
        -0.189478538194e+03_wp, -0.729758019808e+04_wp, -0.447258528240e+04_wp, 0.108022979288e+05_wp, &
        -0.152796723044e+04_wp, -0.443217892640e+04_wp, -0.471318020185e+04_wp, -0.323805553264e+04_wp, &
        -0.474853294031e+04_wp, 0.727870654658e+04_wp, -0.480003990295e+04_wp, -0.901692402681e+03_wp, &
        0.361531617058e+04_wp, -0.129458775524e+04_wp, 0.239872383321e+04_wp, -0.613392039351e+04_wp, &
        -0.909009369483e+03_wp, 0.400238255588e+04_wp, -0.150326338488e+05_wp, 0.218914472743e+04_wp, &
        0.239895086789e+05_wp, 0.114385023318e+04_wp, 0.713054834829e+04_wp, -0.399896928899e+04_wp, &
        -0.471596056754e+04_wp, -0.242573924903e+04_wp, -0.281473702631e+03_wp, 0.359512322699e+04_wp, &
        -0.795359212080e+03_wp, 0.103710393997e+05_wp, -0.114623852531e+04_wp, -0.578299696350e+03_wp, &
        -0.536080316130e+04_wp, -0.372551232771e+04_wp, 0.804771343275e+04_wp, -0.750205359657e+03_wp, &
        0.557336024002e+04_wp, 0.841003412805e+04_wp, 0.227272169978e+04_wp, -0.215450813740e+04_wp, &
        -0.707931699033e+04_wp, 0.361887573377e+04_wp, -0.395197029618e+04_wp, -0.581254265074e+04_wp, &
        0.418895474492e+04_wp, -0.205947976269e+04_wp, 0.111660262741e+04_wp, 0.208384577312e+04_wp, &
        -0.855630467033e+03_wp, 0.323076283511e+03_wp, -0.248160287086e+04_wp, -0.532414460924e+04_wp, &
        0.545271789083e+04_wp, -0.215264356789e+04_wp, 0.364258865680e+03_wp, 0.362514360746e+04_wp, &
        -0.274969056573e+04_wp, 0.655738634214e+03_wp, 0.133875558063e+04_wp, -0.357414511891e+03_wp, &
        -0.243707590740e+04_wp, -0.354551353111e+04_wp, -0.980037368257e+03_wp, -0.203741247593e+04_wp, &
        0.457683411673e+04_wp, 0.264339195897e+03_wp, 0.878262929685e+03_wp, 0.186639375953e+04_wp, &
        -0.190585874076e+04_wp, 0.605305696451e+04_wp, -0.316825009188e+04_wp, -0.242062478394e+04_wp, &
        -0.195288006588e+04_wp, 0.579975066295e+04_wp, 0.342138743584e+04_wp, 0.793139343330e+03_wp, &
        0.949678467586e+03_wp, 0.389147092370e+04_wp, -0.404975571529e+04_wp, -0.122288585004e+04_wp, &
        0.274193648357e+04_wp, -0.156596842910e+04_wp, 0.309688435256e+03_wp, 0.121880743886e+04_wp, &
        -0.243236889006e+04_wp, -0.215167746180e+04_wp, 0.336174192236e+03_wp, 0.472454712468e+03_wp &
        ]
end module n4_1a_meg2_2020_par

!> @brief Implement N4_1A_MEG2_2020, a potential energy surface for N4.
module chempotpy_n4_1a_meg2_2020_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the N4_1A_MEG2_2020 potential energy surface for N4.
    subroutine evaluate( x, igrad, potential, gradient, dvec, status )
        use chempotpy_geometry, only: flatten_cartesian_coordinates, reshape_cartesian_coordinates
        use chempotpy_status,   only: chempotpy_status_t
        use chempotpy_units,    only: electronvolt_to_kilocalorie_per_mole
        implicit none
        interface
            !> @brief Evaluate the released MEG2 kernel in kcal/mol and kcal/(mol angstrom).
            subroutine n4pes( x, v, dvdx, igrad )
                use chempotpy_kinds, only: wp
                implicit none
                integer, intent(in)   :: igrad
                real(wp), intent(in)  :: x(12)
                real(wp), intent(out) :: v, dvdx(12)
            end subroutine n4pes
        end interface
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: natoms  = 4
        integer, parameter    :: nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: potential(nstates)
        real(wp), intent(out) :: gradient(nstates,natoms,3)
        real(wp), intent(out) :: dvec(nstates,nstates,natoms,3)
        real(wp)              :: v, tx(12), tg(12)
        if( present(status) ) call status%clear()
        potential = 0.e0_wp; gradient = 0.e0_wp; dvec = 0.e0_wp
        call flatten_cartesian_coordinates(x, tx)
        tg = 0.e0_wp
        call n4pes(tx, v, tg, igrad)
        v            = v/electronvolt_to_kilocalorie_per_mole
        tg           = tg/electronvolt_to_kilocalorie_per_mole
        potential(1) = v
        call reshape_cartesian_coordinates(tg, gradient(1,:,:))
    end subroutine evaluate
end module chempotpy_n4_1a_meg2_2020_pes

!> @brief Expose N4_1A_MEG2_2020 through the Python-facing PES entry point.
subroutine pes( x, igrad, potential, gradient, dvec )
    use chempotpy_n4_1a_meg2_2020_pes, only: evaluate
    use chempotpy_kinds,               only: wp
    implicit none
    integer, parameter    :: natoms  = 4
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: potential(nstates)
    real(wp), intent(out) :: gradient(nstates,natoms,3)
    real(wp), intent(out) :: dvec(nstates,nstates,natoms,3)
    call evaluate(x, igrad, potential, gradient, dvec)
end subroutine pes

!> @brief Evaluate n4pes.
subroutine n4pes( x, v, dvdx, igrad )
    ! Subroutine to calculate the potential energy V and gradient dVdX
    ! for given Cartesian coordinates X(12)
    ! R:        Interatomic bond distance (6)
    ! V:        Calculated potential energy
    ! dVdX:     The derivative of V w.r.t. X, dim(12)
    ! dVdR:     The derivative of V w.r.t. R, dim(6)
    ! dPdR:     The derivative of basis functions w.r.t. R
    !       dim(6*306)
    ! dMdR:     The derivative of monomials w.r.t. R
    !       dim(6*112)
    ! dRdX:     The derivative of R w.r.t. X, dim(6*12)
    use chempotpy_kinds, only: wp
    implicit none
    interface
        !> @brief Evaluate coord_convt.
        subroutine coord_convt( x )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp), intent(in) :: x(12)
        end subroutine coord_convt

        !> @brief Evaluate evdvdx.
        subroutine evdvdx( x, dvdx )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp), intent(in)  :: x(12)
            real(wp), intent(out) :: dvdx(12)
        end subroutine evdvdx

        !> @brief Evaluate evv.
        subroutine evv( v )
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i, j, k
            real(wp) :: dist, dv2dr, v, v2
        end subroutine evv
    end interface
    integer  :: igrad
    real(wp) :: v
    real(wp) :: dvdx(12), x(12)
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
    endif
end subroutine n4pes

!> @brief Evaluate coord_convt.
subroutine coord_convt( x )
    use chempotpy_geometry,  only: all_pair_distances
    use n4_1a_meg2_2020_par, only: r
    !  Program to calculate the six interatomic distance
    !  from the Cartesian coordinate array
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), intent(in) :: x(12)
    call all_pair_distances(x, r)
end subroutine coord_convt

!> @brief Evaluate EvV.
subroutine evv( v )
    use n4_1a_meg2_2020_par, only: b, c, r, totdiss
    ! Subroutine to evaluate V for given R
    ! V(R) = C*P
    ! C:        Coefficients, stored in 'dim.inc'
    ! P:        Basis functions evaluated for given R
    ! rMs:      rMs(6), six mixed exponential gaussian terms (MEG)
    ! a:        Nonlinear parameters in Morse terms(Angstrom)
    ! ab:       Nonlinear parameters in Gauss terms(Angstrom^2)
    ! re:       Equilibrium bond length(Angstrom)
    ! nop:      number of points
    ! nom:      number of monomials
    ! nob:      number of basis functions(polynomials)
    ! rM(0:111):    Array to store monomials
    ! P(0:305): Array to store polynomials
    ! B(1:276):     Array to store basis functions
    use chempotpy_kinds, only: wp
    implicit none
    interface
        !> @brief Evaluate ev2gm2.
        subroutine ev2gm2( r, v, grad, imol, igrad )
            use chempotpy_kinds, only: wp
            implicit none
            integer, intent(in) :: imol, igrad
            real(wp)            :: r
            real(wp)            :: v, grad
            real(wp)            :: re, de, as(0:6)
            real(wp)            :: y, fy, u
            real(wp)            :: dfdy, rr, dydr, dfdr
            integer             :: k
            real(wp)            :: dist(1), disp, dispdr(1)
        end subroutine ev2gm2

        !> @brief Evaluate evbas.
        subroutine evbas
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i
            real(wp) :: b1(306)
        end subroutine evbas

        !> @brief Evaluate evmono.
        subroutine evmono
            use chempotpy_kinds, only: wp
            implicit none
        end subroutine evmono

        !> @brief Evaluate evmorse.
        subroutine evmorse
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i
        end subroutine evmorse

        !> @brief Evaluate evpoly.
        subroutine evpoly
            use chempotpy_kinds, only: wp
            implicit none
        end subroutine evpoly
    end interface
    integer  :: i
    real(wp) :: dist, dv2dr, v, v2
    ! Calculate the six MEG terms for each point
    call evmorse
    ! Calculate the monomials for each point by using six MEG terms
    call evmono
    ! Calculate the polynomials (basis functions) by using monomials
    call evpoly
    ! Calculate the basis functions by removing unconnected and 2-body terms
    call evbas
    ! Initialized v to be totdiss
    v = totdiss
    ! Evaluate 2-body interactions
    do i=1, 6
        dist = r(i)
        call ev2gm2(dist, v2, dv2dr, 1, 0)
        v = v+v2
    enddo
    ! Evaluate V by taken the product of C and Basis function array
    do i=1, 276
        v = v + c(i)*b(i)
    enddo
end subroutine evv

!> @brief Evaluate EvdVdX.
subroutine evdvdx( x, dvdx )
    use chempotpy_geometry,  only: transform_internal_gradient_to_cartesian
    use n4_1a_meg2_2020_par, only: drdx, dvdr
    ! Subroutine to evaluate dRdX for given R and X
    ! R:        R(6), 6 bond lengths
    ! X:        X(12), 12 Cartesian coordinates
    ! rM(0:111):    Array to store monomials
    ! P(0:305):     Array to store polynomials
    ! dVdX:     dVdX(12), derivatives of V w.r.t. Cartesian coordinates
    ! dVdR:     dVdR(6), derivatives of V w.r.t.6 bond length
    ! dRdX:     dRdX(6,12), derivatives of R(6) w.r.t. 12
    !       Cartesian coordinates
    use chempotpy_kinds, only: wp
    implicit none
    interface
        !> @brief Evaluate evdrdx.
        subroutine evdrdx( x )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp), intent(in) :: x(12)
        end subroutine evdrdx

        !> @brief Evaluate evdvdr.
        subroutine evdvdr
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i, j
            real(wp) :: dist, v2, dv2dr
        end subroutine evdvdr
    end interface
    real(wp), intent(in)  :: x(12)
    real(wp), intent(out) :: dvdx(12)
    call evdvdr
    call evdrdx(x)
    call transform_internal_gradient_to_cartesian(dvdr, drdx, dvdx)
end subroutine evdvdx

!> @brief Evaluate EvMorse.
subroutine evmorse
    use n4_1a_meg2_2020_par, only: a, ab, r, re, rms
    ! mixed exponential gaussian term rms = exp(-(r-re)/a-(r-re)^2/ab)
    ! re:   equlibrium bond length
    ! a:    nonlinear parameter, unit Ang
    ! ab:   nonlinear parameter, unit Ang^2
    use chempotpy_kinds, only: wp
    implicit none
    integer :: i
    do i=1, 6
        rms(i) = exp(-(r(i)-re)/a-((r(i)-re)**2.0e0_wp)/ab)
    enddo
end subroutine evmorse

!> @brief Evaluate EvMono.
subroutine evmono
    use n4_1a_meg2_2020_par, only: rm, rms
    !  The subroutine accepts six MEG variables (X) and calculates the
    !  monomials(M) that do not have usable decomposition.
    !  For A4 with max. degree 9, the number of monomials is nom.
    use chempotpy_kinds, only: wp
    implicit none
    rm(0)   = 1.0e0_wp
    rm(1)   = rms(6)
    rm(2)   = rms(5)
    rm(3)   = rms(4)
    rm(4)   = rms(3)
    rm(5)   = rms(2)
    rm(6)   = rms(1)
    rm(7)   = rm(3)*rm(4)
    rm(8)   = rm(2)*rm(5)
    rm(9)   = rm(1)*rm(6)
    rm(10)  = rm(1)*rm(2)
    rm(11)  = rm(1)*rm(3)
    rm(12)  = rm(2)*rm(3)
    rm(13)  = rm(1)*rm(4)
    rm(14)  = rm(2)*rm(4)
    rm(15)  = rm(1)*rm(5)
    rm(16)  = rm(3)*rm(5)
    rm(17)  = rm(4)*rm(5)
    rm(18)  = rm(2)*rm(6)
    rm(19)  = rm(3)*rm(6)
    rm(20)  = rm(4)*rm(6)
    rm(21)  = rm(5)*rm(6)
    rm(22)  = rm(1)*rm(7)
    rm(23)  = rm(2)*rm(7)
    rm(24)  = rm(1)*rm(8)
    rm(25)  = rm(2)*rm(16)
    rm(26)  = rm(2)*rm(17)
    rm(27)  = rm(3)*rm(17)
    rm(28)  = rm(1)*rm(18)
    rm(29)  = rm(1)*rm(19)
    rm(30)  = rm(1)*rm(20)
    rm(31)  = rm(3)*rm(20)
    rm(32)  = rm(1)*rm(21)
    rm(33)  = rm(2)*rm(21)
    rm(34)  = rm(1)*rm(12)
    rm(35)  = rm(1)*rm(17)
    rm(36)  = rm(2)*rm(20)
    rm(37)  = rm(3)*rm(21)
    rm(38)  = rm(1)*rm(14)
    rm(39)  = rm(1)*rm(16)
    rm(40)  = rm(2)*rm(19)
    rm(41)  = rm(4)*rm(21)
    rm(42)  = rm(2)*rm(27)
    rm(43)  = rm(1)*rm(31)
    rm(44)  = rm(1)*rm(33)
    rm(45)  = rm(1)*rm(23)
    rm(46)  = rm(1)*rm(25)
    rm(47)  = rm(1)*rm(26)
    rm(48)  = rm(1)*rm(27)
    rm(49)  = rm(1)*rm(40)
    rm(50)  = rm(1)*rm(36)
    rm(51)  = rm(2)*rm(31)
    rm(52)  = rm(1)*rm(37)
    rm(53)  = rm(2)*rm(37)
    rm(54)  = rm(1)*rm(41)
    rm(55)  = rm(2)*rm(41)
    rm(56)  = rm(3)*rm(41)
    rm(57)  = rm(1)*rm(42)
    rm(58)  = rm(1)*rm(51)
    rm(59)  = rm(1)*rm(53)
    rm(60)  = rm(1)*rm(55)
    rm(61)  = rm(1)*rm(56)
    rm(62)  = rm(2)*rm(56)
    rm(63)  = rm(1)*rm(62)
    rm(64)  = rm(2)*rm(57)
    rm(65)  = rm(3)*rm(57)
    rm(66)  = rm(4)*rm(57)
    rm(67)  = rm(5)*rm(57)
    rm(68)  = rm(1)*rm(58)
    rm(69)  = rm(3)*rm(58)
    rm(70)  = rm(4)*rm(58)
    rm(71)  = rm(1)*rm(59)
    rm(72)  = rm(2)*rm(59)
    rm(73)  = rm(1)*rm(60)
    rm(74)  = rm(2)*rm(60)
    rm(75)  = rm(1)*rm(61)
    rm(76)  = rm(2)*rm(62)
    rm(77)  = rm(3)*rm(61)
    rm(78)  = rm(3)*rm(62)
    rm(79)  = rm(4)*rm(61)
    rm(80)  = rm(4)*rm(62)
    rm(81)  = rm(5)*rm(59)
    rm(82)  = rm(5)*rm(60)
    rm(83)  = rm(5)*rm(62)
    rm(84)  = rm(6)*rm(58)
    rm(85)  = rm(6)*rm(59)
    rm(86)  = rm(6)*rm(60)
    rm(87)  = rm(6)*rm(61)
    rm(88)  = rm(2)*rm(64)
    rm(89)  = rm(3)*rm(65)
    rm(90)  = rm(4)*rm(66)
    rm(91)  = rm(5)*rm(67)
    rm(92)  = rm(1)*rm(68)
    rm(93)  = rm(3)*rm(69)
    rm(94)  = rm(4)*rm(70)
    rm(95)  = rm(1)*rm(71)
    rm(96)  = rm(2)*rm(72)
    rm(97)  = rm(1)*rm(73)
    rm(98)  = rm(2)*rm(74)
    rm(99)  = rm(1)*rm(75)
    rm(100) = rm(2)*rm(76)
    rm(101) = rm(3)*rm(77)
    rm(102) = rm(3)*rm(78)
    rm(103) = rm(4)*rm(79)
    rm(104) = rm(4)*rm(80)
    rm(105) = rm(5)*rm(81)
    rm(106) = rm(5)*rm(82)
    rm(107) = rm(5)*rm(83)
    rm(108) = rm(6)*rm(84)
    rm(109) = rm(6)*rm(85)
    rm(110) = rm(6)*rm(86)
    rm(111) = rm(6)*rm(87)
end subroutine evmono

!> @brief Evaluate EvPoly.
subroutine evpoly
    use n4_1a_meg2_2020_par, only: p, rm
    !  The subroutine accepts monomials (m) and calculates the
    !  permutation-invariant polynomials(p)
    !  For A4 with max. degree 9, the number of polynomials is nob.
    implicit none
    p(0)  = rm(0)
    p(1)  = rm(1) + rm(2) + rm(3) + rm(4) + rm(5) + rm(6)
    p(2)  = rm(7) + rm(8) + rm(9)
    p(3)  = rm(10) + rm(11) + rm(12) + rm(13) + rm(14) + rm(15) + rm(16) + rm(17) + rm(18) + rm(19) + rm(20) + rm(21)
    p(4)  = p(1)*p(1) - p(3) - p(2) - p(3) - p(2)
    p(5)  = rm(22) + rm(23) + rm(24) + rm(25) + rm(26) + rm(27) + rm(28) + rm(29) + rm(30) + rm(31) + rm(32) + rm(33)
    p(6)  = rm(34) + rm(35) + rm(36) + rm(37)
    p(7)  = rm(38) + rm(39) + rm(40) + rm(41)
    p(8)  = p(1)*p(2) - p(5)
    p(9)  = p(1)*p(3) - p(6) - p(7) - p(5) - p(6) - p(7) - p(5) - p(6) - p(7)
    p(10) = p(1)*p(4) - p(9) - p(8)
    p(11) = rm(42) + rm(43) + rm(44)
    p(12) = rm(45) + rm(46) + rm(47) + rm(48) + rm(49) + rm(50) + rm(51) + rm(52) + rm(53) + rm(54) + rm(55) + rm(56)
    p(13) = p(2)*p(3) - p(12)
    p(14) = p(1)*p(5) - p(12) - p(11) - p(13) - p(12) - p(11) - p(11) - p(11)
    p(15) = p(1)*p(6) - p(12)
    p(16) = p(1)*p(7) - p(12)
    p(17) = p(2)*p(2) - p(11) - p(11)
    p(18) = p(3)*p(3) - p(12) - p(11) - p(15) - p(16) - p(14) - p(12) - p(11) - p(15) - p(16) - p(14) - p(12) - &
        & p(11) - p(12) - p(11)
    p(19) = p(2)*p(4) - p(14)
    p(20) = p(3)*p(4) - p(15) - p(16) - p(13)
    p(21) = p(1)*p(10) - p(20) - p(19)
    p(22) = rm(57) + rm(58) + rm(59) + rm(60) + rm(61) + rm(62)
    p(23) = p(1)*p(11) - p(22)
    p(24) = p(2)*p(6)
    p(25) = p(2)*p(7)
    p(26) = p(1)*p(12) - p(22) - p(24) - p(25) - p(22) - p(22) - p(22)
    p(27) = p(2)*p(5) - p(22) - p(23) - p(22)
    p(28) = p(3)*p(5) - p(22) - p(26) - p(24) - p(25) - p(23) - p(22) - p(24) - p(25) - p(23) - p(22) - p(22)
    p(29) = p(3)*p(6) - p(22) - p(26) - p(22)
    p(30) = p(3)*p(7) - p(22) - p(26) - p(22)
    p(31) = p(2)*p(9) - p(26) - p(28)
    p(32) = p(1)*p(14) - p(26) - p(23) - p(28)
    p(33) = p(4)*p(6) - p(25)
    p(34) = p(4)*p(7) - p(24)
    p(35) = p(1)*p(17) - p(27)
    p(36) = p(1)*p(18) - p(29) - p(30) - p(28)
    p(37) = p(2)*p(10) - p(32)
    p(38) = p(3)*p(10) - p(33) - p(34) - p(31)
    p(39) = p(1)*p(21) - p(38) - p(37)
    p(40) = rm(63)
    p(41) = rm(64) + rm(65) + rm(66) + rm(67) + rm(68) + rm(69) + rm(70) + rm(71) + rm(72) + rm(73) + rm(74) + &
        & rm(75) + rm(76) + rm(77) + rm(78) + rm(79) + rm(80) + rm(81) + rm(82) + rm(83) + rm(84) + rm(85) + rm(86) + &
        & rm(87)
    p(42) = p(1)*p(22) - p(40) - p(41) - p(40) - p(40) - p(40) - p(40) - p(40)
    p(43) = p(2)*p(11) - p(40) - p(40) - p(40)
    p(44) = p(2)*p(12) - p(41)
    p(45) = p(3)*p(11) - p(41)
    p(46) = p(5)*p(6) - p(41)
    p(47) = p(5)*p(7) - p(41)
    p(48) = p(6)*p(7) - p(40) - p(40) - p(40) - p(40)
    p(49) = p(4)*p(11) - p(42)
    p(50) = p(2)*p(15) - p(46)
    p(51) = p(2)*p(16) - p(47)
    p(52) = p(4)*p(12) - p(41) - p(50) - p(51)
    p(53) = p(2)*p(14) - p(42) - p(49) - p(42)
    p(54) = p(6)*p(6) - p(42) - p(42)
    p(55) = p(7)*p(7) - p(42) - p(42)
    p(56) = p(3)*p(17) - p(44)
    p(57) = p(2)*p(18) - p(48)
    p(58) = p(3)*p(14) - p(41) - p(52) - p(46) - p(47) - p(45) - p(45)
    p(59) = p(6)*p(9) - p(41) - p(52) - p(47)
    p(60) = p(7)*p(9) - p(41) - p(52) - p(46)
    p(61) = p(2)*p(20) - p(52) - p(58)
    p(62) = p(1)*p(32) - p(52) - p(49) - p(58)
    p(63) = p(6)*p(10) - p(51)
    p(64) = p(7)*p(10) - p(50)
    p(65) = p(2)*p(17) - p(43)
    p(66) = p(3)*p(18) - p(46) - p(47) - p(45) - p(59) - p(60) - p(58)
    p(67) = p(2)*p(19) - p(49)
    p(68) = p(1)*p(36) - p(59) - p(60) - p(58) - p(57) - p(66) - p(66)
    p(69) = p(2)*p(21) - p(62)
    p(70) = p(3)*p(21) - p(63) - p(64) - p(61)
    p(71) = p(1)*p(39) - p(70) - p(69)
    p(72) = p(40)*p(1)
    p(73) = p(2)*p(22) - p(72)
    p(74) = p(6)*p(11)
    p(75) = p(7)*p(11)
    p(76) = p(3)*p(22) - p(72) - p(74) - p(75) - p(72) - p(72) - p(72)
    p(77) = rm(88) + rm(89) + rm(90) + rm(91) + rm(92) + rm(93) + rm(94) + rm(95) + rm(96) + rm(97) + rm(98) + &
        & rm(99) + rm(100) + rm(101) + rm(102) + rm(103) + rm(104) + rm(105) + rm(106) + rm(107) + rm(108) + rm(109) + &
        & rm(110) + rm(111)
    p(78)  = p(1)*p(42) - p(72) - p(76)
    p(79)  = p(5)*p(11) - p(72) - p(73) - p(72)
    p(80)  = p(2)*p(26) - p(76) - p(77)
    p(81)  = p(6)*p(12) - p(72) - p(76) - p(72)
    p(82)  = p(7)*p(12) - p(72) - p(76) - p(72)
    p(83)  = p(8)*p(11) - p(72)
    p(84)  = p(6)*p(17)
    p(85)  = p(7)*p(17)
    p(86)  = p(9)*p(11) - p(76) - p(77)
    p(87)  = p(2)*p(29) - p(81)
    p(88)  = p(6)*p(14) - p(75) - p(75)
    p(89)  = p(2)*p(30) - p(82)
    p(90)  = p(7)*p(14) - p(74) - p(74)
    p(91)  = p(1)*p(48) - p(76) - p(81) - p(82)
    p(92)  = p(10)*p(11) - p(78)
    p(93)  = p(2)*p(33) - p(88)
    p(94)  = p(2)*p(34) - p(90)
    p(95)  = p(10)*p(12) - p(77) - p(93) - p(94)
    p(96)  = p(2)*p(27) - p(73) - p(79)
    p(97)  = p(2)*p(28) - p(76) - p(86)
    p(98)  = p(1)*p(53) - p(80) - p(79) - p(97)
    p(99)  = p(1)*p(54) - p(81)
    p(100) = p(1)*p(55) - p(82)
    p(101) = p(5)*p(18) - p(76) - p(91) - p(87) - p(89) - p(86)
    p(102) = p(6)*p(18) - p(74) - p(90)
    p(103) = p(7)*p(18) - p(75) - p(88)
    p(104) = p(2)*p(31) - p(77) - p(86)
    p(105) = p(2)*p(36) - p(91) - p(101)
    p(106) = p(3)*p(32) - p(77) - p(95) - p(88) - p(90) - p(86)
    p(107) = p(4)*p(29) - p(82) - p(80) - p(99)
    p(108) = p(4)*p(30) - p(81) - p(80) - p(100)
    p(109) = p(2)*p(38) - p(95) - p(106)
    p(110) = p(1)*p(62) - p(95) - p(92) - p(106)
    p(111) = p(6)*p(21) - p(94)
    p(112) = p(7)*p(21) - p(93)
    p(113) = p(1)*p(65) - p(96)
    p(114) = p(1)*p(66) - p(102) - p(103) - p(101)
    p(115) = p(2)*p(37) - p(92)
    p(116) = p(10)*p(18) - p(99) - p(100) - p(97)
    p(117) = p(2)*p(39) - p(110)
    p(118) = p(3)*p(39) - p(111) - p(112) - p(109)
    p(119) = p(1)*p(71) - p(118) - p(117)
    p(120) = p(40)*p(2)
    p(121) = p(40)*p(3)
    p(122) = p(40)*p(4)
    p(123) = p(11)*p(12) - p(121)
    p(124) = p(2)*p(42) - p(122)
    p(125) = p(6)*p(22) - p(121)
    p(126) = p(7)*p(22) - p(121)
    p(127) = p(2)*p(41) - p(121) - p(123) - p(121)
    p(128) = p(6)*p(23) - p(123)
    p(129) = p(7)*p(23) - p(123)
    p(130) = p(3)*p(41) - p(122) - p(121) - p(120) - p(125) - p(128) - p(126) - p(129) - p(124) - p(123) - p(122) - &
        & p(121) - p(120) - p(125) - p(126) - p(124) - p(123) - p(122) - p(121) - p(120) - p(122) - p(121) - p(120) - &
        & p(120) - p(120) - p(120) - p(120)
    p(131) = p(3)*p(42) - p(121) - p(125) - p(126) - p(121)
    p(132) = p(4)*p(41) - p(121) - p(131) - p(128) - p(129) - p(127) - p(121)
    p(133) = p(1)*p(78) - p(122) - p(131)
    p(134) = p(11)*p(11) - p(120) - p(120)
    p(135) = p(2)*p(48) - p(130)
    p(136) = p(11)*p(17) - p(120)
    p(137) = p(2)*p(44) - p(123)
    p(138) = p(2)*p(45) - p(121)
    p(139) = p(11)*p(14) - p(122) - p(124) - p(122)
    p(140) = p(6)*p(27) - p(127)
    p(141) = p(2)*p(54)
    p(142) = p(7)*p(27) - p(127)
    p(143) = p(2)*p(52) - p(131) - p(132)
    p(144) = p(1)*p(81) - p(125) - p(141) - p(135) - p(125)
    p(145) = p(2)*p(55)
    p(146) = p(1)*p(82) - p(126) - p(135) - p(145) - p(126)
    p(147) = p(11)*p(18) - p(130)
    p(148) = p(6)*p(28) - p(129) - p(123) - p(143)
    p(149) = p(7)*p(28) - p(128) - p(123) - p(143)
    p(150) = p(6)*p(30) - p(121) - p(146)
    p(151) = p(11)*p(19) - p(122)
    p(152) = p(2)*p(50) - p(128)
    p(153) = p(2)*p(51) - p(129)
    p(154) = p(11)*p(20) - p(131) - p(132)
    p(155) = p(2)*p(59) - p(144) - p(148)
    p(156) = p(6)*p(32) - p(129)
    p(157) = p(2)*p(60) - p(146) - p(149)
    p(158) = p(7)*p(32) - p(128)
    p(159) = p(6)*p(34) - p(122) - p(145) - p(122)
    p(160) = p(11)*p(21) - p(133)
    p(161) = p(2)*p(63) - p(156)
    p(162) = p(2)*p(64) - p(158)
    p(163) = p(12)*p(21) - p(132) - p(161) - p(162)
    p(164) = p(2)*p(53) - p(124) - p(139)
    p(165) = p(2)*p(58) - p(131) - p(154)
    p(166) = p(3)*p(54) - p(125) - p(144)
    p(167) = p(3)*p(55) - p(126) - p(146)
    p(168) = p(3)*p(65) - p(137)
    p(169) = p(17)*p(18) - p(135)
    p(170) = p(1)*p(98) - p(143) - p(139) - p(165)
    p(171) = p(4)*p(54) - p(135)
    p(172) = p(4)*p(55) - p(135)
    p(173) = p(2)*p(66) - p(150)
    p(174) = p(1)*p(101) - p(148) - p(149) - p(147) - p(173) - p(165) - p(147)
    p(175) = p(1)*p(102) - p(150) - p(148) - p(166)
    p(176) = p(1)*p(103) - p(150) - p(149) - p(167)
    p(177) = p(2)*p(61) - p(132) - p(154)
    p(178) = p(2)*p(68) - p(159) - p(174)
    p(179) = p(3)*p(62) - p(132) - p(163) - p(156) - p(158) - p(154)
    p(180) = p(6)*p(38) - p(132) - p(163) - p(157)
    p(181) = p(7)*p(38) - p(132) - p(163) - p(155)
    p(182) = p(2)*p(70) - p(163) - p(179)
    p(183) = p(1)*p(110) - p(163) - p(160) - p(179)
    p(184) = p(6)*p(39) - p(162)
    p(185) = p(7)*p(39) - p(161)
    p(186) = p(2)*p(65) - p(136)
    p(187) = p(3)*p(66) - p(148) - p(149) - p(147) - p(175) - p(176) - p(174)
    p(188) = p(2)*p(67) - p(151)
    p(189) = p(4)*p(66) - p(166) - p(167) - p(165)
    p(190) = p(2)*p(69) - p(160)
    p(191) = p(18)*p(21) - p(171) - p(172) - p(169)
    p(192) = p(2)*p(71) - p(183)
    p(193) = p(3)*p(71) - p(184) - p(185) - p(182)
    p(194) = p(1)*p(119) - p(193) - p(192)
    p(195) = p(40)*p(5)
    p(196) = p(40)*p(6)
    p(197) = p(40)*p(7)
    p(198) = p(40)*p(8)
    p(199) = p(40)*p(9)
    p(200) = p(40)*p(10)
    p(201) = p(11)*p(22) - p(195)
    p(202) = p(12)*p(22) - p(196) - p(197) - p(195) - p(196) - p(197) - p(195) - p(196) - p(197)
    p(203) = p(17)*p(22) - p(198)
    p(204) = p(6)*p(43)
    p(205) = p(7)*p(43)
    p(206) = p(11)*p(26) - p(199) - p(202)
    p(207) = p(2)*p(76) - p(199) - p(202)
    p(208) = p(2)*p(78) - p(200)
    p(209) = p(6)*p(41) - p(199) - p(195) - p(202) - p(195)
    p(210) = p(6)*p(42) - p(197) - p(197) - p(197)
    p(211) = p(7)*p(41) - p(199) - p(195) - p(202) - p(195)
    p(212) = p(7)*p(42) - p(196) - p(196) - p(196)
    p(213) = p(11)*p(29) - p(209)
    p(214) = p(11)*p(30) - p(211)
    p(215) = p(18)*p(22) - p(199) - p(213) - p(214)
    p(216) = p(2)*p(77) - p(199) - p(206)
    p(217) = p(6)*p(49) - p(205)
    p(218) = p(7)*p(49) - p(204)
    p(219) = p(3)*p(77) - p(200) - p(199) - p(198) - p(209) - p(217) - p(211) - p(218) - p(207) - p(204) - p(205) - &
        & p(200) - p(199) - p(198) - p(200) - p(198) - p(200) - p(198)
    p(220) = p(3)*p(78) - p(199) - p(210) - p(212)
    p(221) = p(10)*p(41) - p(199) - p(220) - p(217) - p(218) - p(216)
    p(222) = p(1)*p(133) - p(200) - p(220)
    p(223) = p(11)*p(27) - p(195) - p(203)
    p(224) = p(1)*p(134) - p(201)
    p(225) = p(2)*p(81) - p(209)
    p(226) = p(2)*p(80) - p(202) - p(206)
    p(227) = p(2)*p(82) - p(211)
    p(228) = p(2)*p(91) - p(215) - p(219)
    p(229) = p(11)*p(28) - p(199) - p(207)
    p(230) = p(6)*p(53) - p(205)
    p(231) = p(5)*p(54) - p(209)
    p(232) = p(7)*p(53) - p(204)
    p(233) = p(7)*p(54) - p(196)
    p(234) = p(5)*p(55) - p(211)
    p(235) = p(6)*p(55) - p(197)
    p(236) = p(11)*p(35) - p(198)
    p(237) = p(6)*p(65)
    p(238) = p(7)*p(65)
    p(239) = p(2)*p(86) - p(199) - p(229)
    p(240) = p(1)*p(139) - p(206) - p(229) - p(224)
    p(241) = p(17)*p(29) - p(225)
    p(242) = p(2)*p(99) - p(231)
    p(243) = p(17)*p(30) - p(227)
    p(244) = p(2)*p(95) - p(220) - p(221)
    p(245) = p(4)*p(81) - p(202) - p(242) - p(227)
    p(246) = p(2)*p(100) - p(234)
    p(247) = p(4)*p(82) - p(202) - p(225) - p(246)
    p(248) = p(11)*p(36) - p(215) - p(219)
    p(249) = p(2)*p(102) - p(233)
    p(250) = p(6)*p(58) - p(206) - p(214) - p(244) - p(214)
    p(251) = p(2)*p(103) - p(235)
    p(252) = p(7)*p(58) - p(213) - p(206) - p(244) - p(213)
    p(253) = p(1)*p(150) - p(215) - p(233) - p(235)
    p(254) = p(11)*p(37) - p(200)
    p(255) = p(2)*p(93) - p(217)
    p(256) = p(2)*p(94) - p(218)
    p(257) = p(11)*p(38) - p(220) - p(221)
    p(258) = p(2)*p(107) - p(245) - p(250)
    p(259) = p(6)*p(62) - p(218)
    p(260) = p(2)*p(108) - p(247) - p(252)
    p(261) = p(7)*p(62) - p(217)
    p(262) = p(6)*p(64) - p(200) - p(246) - p(200)
    p(263) = p(11)*p(39) - p(222)
    p(264) = p(2)*p(111) - p(259)
    p(265) = p(2)*p(112) - p(261)
    p(266) = p(12)*p(39) - p(221) - p(264) - p(265)
    p(267) = p(2)*p(98) - p(208) - p(240)
    p(268) = p(6)*p(54) - p(210)
    p(269) = p(7)*p(55) - p(212)
    p(270) = p(2)*p(96) - p(203) - p(223)
    p(271) = p(2)*p(97) - p(207) - p(229)
    p(272) = p(2)*p(101) - p(215) - p(248)
    p(273) = p(2)*p(106) - p(220) - p(257)
    p(274) = p(6)*p(59) - p(220) - p(215) - p(209)
    p(275) = p(7)*p(60) - p(220) - p(215) - p(211)
    p(276) = p(5)*p(66) - p(215) - p(253) - p(249) - p(251) - p(248)
    p(277) = p(6)*p(66) - p(213) - p(252)
    p(278) = p(7)*p(66) - p(214) - p(250)
    p(279) = p(2)*p(104) - p(216) - p(239)
    p(280) = p(2)*p(105) - p(219) - p(248)
    p(281) = p(1)*p(170) - p(244) - p(240) - p(273)
    p(282) = p(10)*p(54) - p(227)
    p(283) = p(10)*p(55) - p(225)
    p(284) = p(2)*p(114) - p(253) - p(276)
    p(285) = p(1)*p(174) - p(250) - p(252) - p(248) - p(276) - p(273)
    p(286) = p(6)*p(68) - p(217) - p(261) - p(251)
    p(287) = p(7)*p(68) - p(218) - p(259) - p(249)
    p(288) = p(2)*p(109) - p(221) - p(257)
    p(289) = p(2)*p(116) - p(262) - p(285)
    p(290) = p(3)*p(110) - p(221) - p(266) - p(259) - p(261) - p(257)
    p(291) = p(6)*p(70) - p(221) - p(266) - p(260)
    p(292) = p(7)*p(70) - p(221) - p(266) - p(258)
    p(293) = p(2)*p(118) - p(266) - p(290)
    p(294) = p(1)*p(183) - p(266) - p(263) - p(290)
    p(295) = p(6)*p(71) - p(265)
    p(296) = p(7)*p(71) - p(264)
    p(297) = p(1)*p(186) - p(270)
    p(298) = p(1)*p(187) - p(277) - p(278) - p(276)
    p(299) = p(2)*p(115) - p(254)
    p(300) = p(1)*p(189) - p(286) - p(287) - p(285) - p(284) - p(298)
    p(301) = p(2)*p(117) - p(263)
    p(302) = p(18)*p(39) - p(282) - p(283) - p(280)
    p(303) = p(2)*p(119) - p(294)
    p(304) = p(3)*p(119) - p(295) - p(296) - p(293)
    p(305) = p(1)*p(194) - p(304) - p(303)
end subroutine evpoly

!> @brief Evaluate evbas.
subroutine evbas
    use n4_1a_meg2_2020_par, only: b, p
    !  The subroutine eliminates the 2-body terms in Bowman's approach
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: i
    real(wp) :: b1(306)
    !     Pass P(0:305) to BM1(1:306)
    do i=1, 306
        b1(i) = p(i-1)
    enddo
    ! Remove unconnected terms and 2-body terms and pass to B(1:276)
    b(1) = b1(4)
    do i=2, 4
        b(i) = b1(i+4)
    enddo
    b(5) = b1(10)
    do i=6, 11
        b(i) = b1(i+6)
    enddo
    b(12) = b1(19)
    b(13) = b1(21)
    do i=14, 26
        b(i) = b1(i+9)
    enddo
    b(27) = b1(37)
    b(28) = b1(39)
    do i=29, 53
        b(i) = b1(i+12)
    enddo
    b(54) = b1(67)
    b(55) = b1(69)
    b(56) = b1(71)
    do i=57, 97
        b(i) = b1(i+16)
    enddo
    b(98)  = b1(115)
    b(99)  = b1(117)
    b(100) = b1(119)
    do i=101, 166
        b(i) = b1(i+20)
    enddo
    b(167) = b1(188)
    b(168) = b1(190)
    b(169) = b1(192)
    b(170) = b1(194)
    do i=171, 272
        b(i) = b1(i+25)
    enddo
    b(273) = b1(299)
    b(274) = b1(301)
    b(275) = b1(303)
    b(276) = b1(305)
end subroutine evbas

!> @brief Evaluate ev2gm2.
subroutine ev2gm2( r, v, grad, imol, igrad )
    ! Compute the diatomic potential of ground-state N2 (singlet)
    ! Input:  r      interatomic distance in Angstrom
    !         imol   1 - stands for N2 for N4 PES
    !         igrad  1 - gradient is calculated
    !                0 - gradient is not calculated
    ! Output: V      potential in kcal/mol
    !         grad   gradient (kcal/mol)/Angstrom
    use chempotpy_kinds, only: wp
    implicit none
    interface
        !> @brief Evaluate d3disp.
        subroutine d3disp( dist, disp, dispdr, igrad, imol )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: cn(2), s6, s8, rs6, rs8
            real(wp) :: dist(1), e6(1), e8(1), disp, dispdr(1), c6(2)
            real(wp) :: e6dr(1), e8dr(1)
            integer  :: iz(2), mxc(94), i, j, igrad
            real(wp) :: c6ab(94,94,5,5,3)
            real(wp) :: r2r4(94)
            real(wp) :: autoang, autokcal
            integer  :: imol
        end subroutine d3disp
    end interface
    integer, intent(in) :: imol, igrad
    real(wp)            :: r
    real(wp)            :: v, grad
    real(wp)            :: re, de, as(0:6)
    real(wp)            :: y, fy, u
    real(wp)            :: dfdy, rr, dydr, dfdr
    ! Dispersion variables
    real(wp) :: dist(1), disp, dispdr(1)
    ! Modified parameters for D3(BJ)
    re    = 1.098e0_wp
    de    = 225.213e0_wp ! 228.4d0
    as(0) = 2.7475450369759e0_wp  ! 2.71405774451d0
    as(1) = 0.218868498415108e0_wp  ! 1.32757649829d-1
    as(2) = 0.248885765371433e0_wp  ! 2.66756890408d-1
    as(3) = -0.229295466336412e0_wp  ! 1.95350725241d-1
    as(4) = -0.653389048592838e0_wp  !-4.08663480982d-1
    as(5) = 1.03611964035396e0_wp  ! 3.92451705557d-1
    as(6) = 1.71287482791961e0_wp  ! 1.13006674877d0
    v     = 0.e0_wp
    y     = (r*r*r*r - re*re*re*re)/(r*r*r*r + re*re*re*re)
    fy    = as(0) + as(1)*y + as(2)*y*y + as(3)*y*y*y + as(4)*y*y*y*y + as(5)*y*y*y*y*y + as(6)*y*y*y*y*y*y
    u     = exp(-fy*(r-re))
    v     = de*(1.0e0_wp-u)*(1.0e0_wp-u)-de
    ! Add D3 dispersion correction
    dist(1) = r
    call d3disp(dist, disp, dispdr, 0, imol)
    v = v+disp
    ! Compute the gradient if needed
    if( igrad .eq. 1 )then
        grad = 0.e0_wp
        dfdy = as(1) + 2.0e0_wp*as(2)*y + 3.0e0_wp*as(3)*y*y + 4.0e0_wp*as(4)*y*y*y + 5.0e0_wp*as(5)*y*y*y*y + &
            & 6.0e0_wp*as(6)*y*y*y*y*y
        rr   = r*r*r*r + re*re*re*re
        dydr = 8.0e0_wp*r*r*r*re*re*re*re/(rr*rr)
        dfdr = dfdy*dydr
        grad = 2.0e0_wp*de*(1-u)*u*(dfdr*(r-re)+fy)
        ! Add analytical gradient of D3 dispersion correction
        call d3disp(dist, disp, dispdr, 1, imol)
        grad = grad + dispdr(1)
    endif
end subroutine ev2gm2

!> @brief Evaluate EvdVdR.
subroutine evdvdr
    use n4_1a_meg2_2020_par, only: c, dbdr, dvdr, r
    ! Subroutine to evaluate dVdR for given R
    ! dVdR = dV2dR + C*dBdR
    ! C:        Coefficients, stored in 'dim.inc'
    ! P:        Basis functions evaluated for given R
    ! M:        Monomials evaluated for given R
    ! dV2dR:        Gradient of 2-body interactions
    ! dMsdR:    dMsdR(6,6), 6 MEG terms w.r.t. 6 bond length
    ! dMdR:     dMdR(6,nom), nom monomials w.r.t.6 bond length
    ! dPdR:     dPdR(6,nob), nop polynomial basis functions
    !       w.r.t. 6 bond length
    ! nom:      number of monomials
    ! nob:      number of basis functions(polynomials)
    ! M(nom):   Array to store monomials
    ! P(nob):   Array to store polynomials
    use chempotpy_kinds, only: wp
    implicit none
    interface
        !> @brief Evaluate ev2gm2.
        subroutine ev2gm2( r, v, grad, imol, igrad )
            use chempotpy_kinds, only: wp
            implicit none
            integer, intent(in) :: imol, igrad
            real(wp)            :: r
            real(wp)            :: v, grad
            real(wp)            :: re, de, as(0:6)
            real(wp)            :: y, fy, u
            real(wp)            :: dfdy, rr, dydr, dfdr
            integer             :: k
            real(wp)            :: dist(1), disp, dispdr(1)
        end subroutine ev2gm2

        !> @brief Evaluate evdbdr.
        subroutine evdbdr
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i, j
            real(wp) :: db1dr(6,306)
        end subroutine evdbdr

        !> @brief Evaluate evdmdr.
        subroutine evdmdr
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i
        end subroutine evdmdr

        !> @brief Evaluate evdmsdr.
        subroutine evdmsdr
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, j
        end subroutine evdmsdr

        !> @brief Evaluate evdpdr.
        subroutine evdpdr
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i
        end subroutine evdpdr
    end interface
    integer  :: i, j
    real(wp) :: dist, v2, dv2dr
    ! Initialize dVdR(6)
    do i=1, 6
        dvdr(i) = 0.0e0_wp
    enddo
    ! Add dV2dR(i) to dVdR
    do i=1, 6
        dist = r(i)
        call ev2gm2(dist, v2, dv2dr, 1, 1)
        dvdr(i) = dv2dr
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
        do j=1, 276
            dvdr(i) = dvdr(i) + c(j)*dbdr(i,j)
        enddo
    enddo
end subroutine evdvdr

!> @brief Evaluate EvdRdX.
subroutine evdrdx( x )
    use chempotpy_geometry,  only: all_pair_distance_jacobian
    use n4_1a_meg2_2020_par, only: drdx, r
    use chempotpy_kinds,     only: wp
    implicit none
    real(wp), intent(in) :: x(12)
    call all_pair_distance_jacobian(x, r, drdx)
end subroutine evdrdx

!> @brief Evaluate EvdMsdR.
subroutine evdmsdr
    use n4_1a_meg2_2020_par, only: a, ab, dmsdr, r, re
    ! Subroutine to evaluate the derivatives of MEG term X
    ! w.r.t. interatomic distance R(6)
    ! dmsdR:    Local variables, dirm(6,6)
    ! a:        Nonlinear parameter(Angstrom)
    ! ab:       Nonlinear parameter(Angstrom^2)
    ! re:       equilibrium bond length(Angstrom)
    use chempotpy_kinds, only: wp
    implicit none
    integer :: i, j
    ! Initialize dmsdr
    do i=1, 6
        do j=1, 6
            dmsdr(i,j) = 0.0e0_wp
        enddo
    enddo
    ! MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
    ! dmsdr(i,j)=0  i!=j
    do i=1, 6
        dmsdr(i,i) = (-2.0e0_wp*(r(i)-re)/ab-1/a) * exp(-(r(i)-re)/a-((r(i)-re)**2.0e0_wp)/ab)
    enddo
end subroutine evdmsdr

!> @brief Evaluate EvdMdR.
subroutine evdmdr
    use n4_1a_meg2_2020_par, only: dmdr, dmsdr, rm
    !  The subroutine accepts M(nom) and dMSdR(6,6) and calculates the
    !  dMdR(6,nom) that do not have usable decomposition.
    !  For A4 with max. degree 9, the number of monomials is nom.
    use chempotpy_kinds, only: wp
    implicit none
    integer :: i
    do i=1, 6
        dmdr(i,0)   = 0.0e0_wp
        dmdr(i,1)   = dmsdr(i,6)
        dmdr(i,2)   = dmsdr(i,5)
        dmdr(i,3)   = dmsdr(i,4)
        dmdr(i,4)   = dmsdr(i,3)
        dmdr(i,5)   = dmsdr(i,2)
        dmdr(i,6)   = dmsdr(i,1)
        dmdr(i,7)   = dmdr(i,3)*rm(4) + rm(3)*dmdr(i,4)
        dmdr(i,8)   = dmdr(i,2)*rm(5) + rm(2)*dmdr(i,5)
        dmdr(i,9)   = dmdr(i,1)*rm(6) + rm(1)*dmdr(i,6)
        dmdr(i,10)  = dmdr(i,1)*rm(2) + rm(1)*dmdr(i,2)
        dmdr(i,11)  = dmdr(i,1)*rm(3) + rm(1)*dmdr(i,3)
        dmdr(i,12)  = dmdr(i,2)*rm(3) + rm(2)*dmdr(i,3)
        dmdr(i,13)  = dmdr(i,1)*rm(4) + rm(1)*dmdr(i,4)
        dmdr(i,14)  = dmdr(i,2)*rm(4) + rm(2)*dmdr(i,4)
        dmdr(i,15)  = dmdr(i,1)*rm(5) + rm(1)*dmdr(i,5)
        dmdr(i,16)  = dmdr(i,3)*rm(5) + rm(3)*dmdr(i,5)
        dmdr(i,17)  = dmdr(i,4)*rm(5) + rm(4)*dmdr(i,5)
        dmdr(i,18)  = dmdr(i,2)*rm(6) + rm(2)*dmdr(i,6)
        dmdr(i,19)  = dmdr(i,3)*rm(6) + rm(3)*dmdr(i,6)
        dmdr(i,20)  = dmdr(i,4)*rm(6) + rm(4)*dmdr(i,6)
        dmdr(i,21)  = dmdr(i,5)*rm(6) + rm(5)*dmdr(i,6)
        dmdr(i,22)  = dmdr(i,1)*rm(7) + rm(1)*dmdr(i,7)
        dmdr(i,23)  = dmdr(i,2)*rm(7) + rm(2)*dmdr(i,7)
        dmdr(i,24)  = dmdr(i,1)*rm(8) + rm(1)*dmdr(i,8)
        dmdr(i,25)  = dmdr(i,2)*rm(16) + rm(2)*dmdr(i,16)
        dmdr(i,26)  = dmdr(i,2)*rm(17) + rm(2)*dmdr(i,17)
        dmdr(i,27)  = dmdr(i,3)*rm(17) + rm(3)*dmdr(i,17)
        dmdr(i,28)  = dmdr(i,1)*rm(18) + rm(1)*dmdr(i,18)
        dmdr(i,29)  = dmdr(i,1)*rm(19) + rm(1)*dmdr(i,19)
        dmdr(i,30)  = dmdr(i,1)*rm(20) + rm(1)*dmdr(i,20)
        dmdr(i,31)  = dmdr(i,3)*rm(20) + rm(3)*dmdr(i,20)
        dmdr(i,32)  = dmdr(i,1)*rm(21) + rm(1)*dmdr(i,21)
        dmdr(i,33)  = dmdr(i,2)*rm(21) + rm(2)*dmdr(i,21)
        dmdr(i,34)  = dmdr(i,1)*rm(12) + rm(1)*dmdr(i,12)
        dmdr(i,35)  = dmdr(i,1)*rm(17) + rm(1)*dmdr(i,17)
        dmdr(i,36)  = dmdr(i,2)*rm(20) + rm(2)*dmdr(i,20)
        dmdr(i,37)  = dmdr(i,3)*rm(21) + rm(3)*dmdr(i,21)
        dmdr(i,38)  = dmdr(i,1)*rm(14) + rm(1)*dmdr(i,14)
        dmdr(i,39)  = dmdr(i,1)*rm(16) + rm(1)*dmdr(i,16)
        dmdr(i,40)  = dmdr(i,2)*rm(19) + rm(2)*dmdr(i,19)
        dmdr(i,41)  = dmdr(i,4)*rm(21) + rm(4)*dmdr(i,21)
        dmdr(i,42)  = dmdr(i,2)*rm(27) + rm(2)*dmdr(i,27)
        dmdr(i,43)  = dmdr(i,1)*rm(31) + rm(1)*dmdr(i,31)
        dmdr(i,44)  = dmdr(i,1)*rm(33) + rm(1)*dmdr(i,33)
        dmdr(i,45)  = dmdr(i,1)*rm(23) + rm(1)*dmdr(i,23)
        dmdr(i,46)  = dmdr(i,1)*rm(25) + rm(1)*dmdr(i,25)
        dmdr(i,47)  = dmdr(i,1)*rm(26) + rm(1)*dmdr(i,26)
        dmdr(i,48)  = dmdr(i,1)*rm(27) + rm(1)*dmdr(i,27)
        dmdr(i,49)  = dmdr(i,1)*rm(40) + rm(1)*dmdr(i,40)
        dmdr(i,50)  = dmdr(i,1)*rm(36) + rm(1)*dmdr(i,36)
        dmdr(i,51)  = dmdr(i,2)*rm(31) + rm(2)*dmdr(i,31)
        dmdr(i,52)  = dmdr(i,1)*rm(37) + rm(1)*dmdr(i,37)
        dmdr(i,53)  = dmdr(i,2)*rm(37) + rm(2)*dmdr(i,37)
        dmdr(i,54)  = dmdr(i,1)*rm(41) + rm(1)*dmdr(i,41)
        dmdr(i,55)  = dmdr(i,2)*rm(41) + rm(2)*dmdr(i,41)
        dmdr(i,56)  = dmdr(i,3)*rm(41) + rm(3)*dmdr(i,41)
        dmdr(i,57)  = dmdr(i,1)*rm(42) + rm(1)*dmdr(i,42)
        dmdr(i,58)  = dmdr(i,1)*rm(51) + rm(1)*dmdr(i,51)
        dmdr(i,59)  = dmdr(i,1)*rm(53) + rm(1)*dmdr(i,53)
        dmdr(i,60)  = dmdr(i,1)*rm(55) + rm(1)*dmdr(i,55)
        dmdr(i,61)  = dmdr(i,1)*rm(56) + rm(1)*dmdr(i,56)
        dmdr(i,62)  = dmdr(i,2)*rm(56) + rm(2)*dmdr(i,56)
        dmdr(i,63)  = dmdr(i,1)*rm(62) + rm(1)*dmdr(i,62)
        dmdr(i,64)  = dmdr(i,2)*rm(57) + rm(2)*dmdr(i,57)
        dmdr(i,65)  = dmdr(i,3)*rm(57) + rm(3)*dmdr(i,57)
        dmdr(i,66)  = dmdr(i,4)*rm(57) + rm(4)*dmdr(i,57)
        dmdr(i,67)  = dmdr(i,5)*rm(57) + rm(5)*dmdr(i,57)
        dmdr(i,68)  = dmdr(i,1)*rm(58) + rm(1)*dmdr(i,58)
        dmdr(i,69)  = dmdr(i,3)*rm(58) + rm(3)*dmdr(i,58)
        dmdr(i,70)  = dmdr(i,4)*rm(58) + rm(4)*dmdr(i,58)
        dmdr(i,71)  = dmdr(i,1)*rm(59) + rm(1)*dmdr(i,59)
        dmdr(i,72)  = dmdr(i,2)*rm(59) + rm(2)*dmdr(i,59)
        dmdr(i,73)  = dmdr(i,1)*rm(60) + rm(1)*dmdr(i,60)
        dmdr(i,74)  = dmdr(i,2)*rm(60) + rm(2)*dmdr(i,60)
        dmdr(i,75)  = dmdr(i,1)*rm(61) + rm(1)*dmdr(i,61)
        dmdr(i,76)  = dmdr(i,2)*rm(62) + rm(2)*dmdr(i,62)
        dmdr(i,77)  = dmdr(i,3)*rm(61) + rm(3)*dmdr(i,61)
        dmdr(i,78)  = dmdr(i,3)*rm(62) + rm(3)*dmdr(i,62)
        dmdr(i,79)  = dmdr(i,4)*rm(61) + rm(4)*dmdr(i,61)
        dmdr(i,80)  = dmdr(i,4)*rm(62) + rm(4)*dmdr(i,62)
        dmdr(i,81)  = dmdr(i,5)*rm(59) + rm(5)*dmdr(i,59)
        dmdr(i,82)  = dmdr(i,5)*rm(60) + rm(5)*dmdr(i,60)
        dmdr(i,83)  = dmdr(i,5)*rm(62) + rm(5)*dmdr(i,62)
        dmdr(i,84)  = dmdr(i,6)*rm(58) + rm(6)*dmdr(i,58)
        dmdr(i,85)  = dmdr(i,6)*rm(59) + rm(6)*dmdr(i,59)
        dmdr(i,86)  = dmdr(i,6)*rm(60) + rm(6)*dmdr(i,60)
        dmdr(i,87)  = dmdr(i,6)*rm(61) + rm(6)*dmdr(i,61)
        dmdr(i,88)  = dmdr(i,2)*rm(64) + rm(2)*dmdr(i,64)
        dmdr(i,89)  = dmdr(i,3)*rm(65) + rm(3)*dmdr(i,65)
        dmdr(i,90)  = dmdr(i,4)*rm(66) + rm(4)*dmdr(i,66)
        dmdr(i,91)  = dmdr(i,5)*rm(67) + rm(5)*dmdr(i,67)
        dmdr(i,92)  = dmdr(i,1)*rm(68) + rm(1)*dmdr(i,68)
        dmdr(i,93)  = dmdr(i,3)*rm(69) + rm(3)*dmdr(i,69)
        dmdr(i,94)  = dmdr(i,4)*rm(70) + rm(4)*dmdr(i,70)
        dmdr(i,95)  = dmdr(i,1)*rm(71) + rm(1)*dmdr(i,71)
        dmdr(i,96)  = dmdr(i,2)*rm(72) + rm(2)*dmdr(i,72)
        dmdr(i,97)  = dmdr(i,1)*rm(73) + rm(1)*dmdr(i,73)
        dmdr(i,98)  = dmdr(i,2)*rm(74) + rm(2)*dmdr(i,74)
        dmdr(i,99)  = dmdr(i,1)*rm(75) + rm(1)*dmdr(i,75)
        dmdr(i,100) = dmdr(i,2)*rm(76) + rm(2)*dmdr(i,76)
        dmdr(i,101) = dmdr(i,3)*rm(77) + rm(3)*dmdr(i,77)
        dmdr(i,102) = dmdr(i,3)*rm(78) + rm(3)*dmdr(i,78)
        dmdr(i,103) = dmdr(i,4)*rm(79) + rm(4)*dmdr(i,79)
        dmdr(i,104) = dmdr(i,4)*rm(80) + rm(4)*dmdr(i,80)
        dmdr(i,105) = dmdr(i,5)*rm(81) + rm(5)*dmdr(i,81)
        dmdr(i,106) = dmdr(i,5)*rm(82) + rm(5)*dmdr(i,82)
        dmdr(i,107) = dmdr(i,5)*rm(83) + rm(5)*dmdr(i,83)
        dmdr(i,108) = dmdr(i,6)*rm(84) + rm(6)*dmdr(i,84)
        dmdr(i,109) = dmdr(i,6)*rm(85) + rm(6)*dmdr(i,85)
        dmdr(i,110) = dmdr(i,6)*rm(86) + rm(6)*dmdr(i,86)
        dmdr(i,111) = dmdr(i,6)*rm(87) + rm(6)*dmdr(i,87)
    enddo
end subroutine evdmdr

!> @brief Evaluate EvdPdr.
subroutine evdpdr
    use n4_1a_meg2_2020_par, only: dmdr, dpdr, p
    !  The subroutine accepts monomials (m) and calculates the
    !  permutation-invariant polynomials(p)
    !  For A4 with max. degree 9, the number of polynomials is nob.
    implicit none
    integer :: i
    do i=1, 6
        dpdr(i,0) = dmdr(i,0)
        dpdr(i,1) = dmdr(i,1) + dmdr(i,2) + dmdr(i,3) + dmdr(i,4) + dmdr(i,5) + dmdr(i,6)
        dpdr(i,2) = dmdr(i,7) + dmdr(i,8) + dmdr(i,9)
        dpdr(i,3) = dmdr(i,10) + dmdr(i,11) + dmdr(i,12) + dmdr(i,13) + dmdr(i,14) + dmdr(i,15) + dmdr(i,16) + dmdr(i, &
            & 17) + dmdr(i,18) + dmdr(i,19) + dmdr(i,20) + dmdr(i,21)
        dpdr(i,4) = dpdr(i,1)*p(1) + p(1)*dpdr(i,1) - dpdr(i,3) - dpdr(i,2) - dpdr(i,3) - dpdr(i,2)
        dpdr(i,5) = dmdr(i,22) + dmdr(i,23) + dmdr(i,24) + dmdr(i,25) + dmdr(i,26) + dmdr(i,27) + dmdr(i,28) + dmdr(i, &
            & 29) + dmdr(i,30) + dmdr(i,31) + dmdr(i,32) + dmdr(i,33)
        dpdr(i,6) = dmdr(i,34) + dmdr(i,35) + dmdr(i,36) + dmdr(i,37)
        dpdr(i,7) = dmdr(i,38) + dmdr(i,39) + dmdr(i,40) + dmdr(i,41)
        dpdr(i,8) = dpdr(i,1)*p(2) + p(1)*dpdr(i,2) - dpdr(i,5)
        dpdr(i,9) = dpdr(i,1)*p(3) + p(1)*dpdr(i,3) - dpdr(i,6) - dpdr(i,7) - dpdr(i,5) - dpdr(i,6) - dpdr(i,7) - &
            & dpdr(i,5) - dpdr(i,6) - dpdr(i,7)
        dpdr(i,10) = dpdr(i,1)*p(4) + p(1)*dpdr(i,4) - dpdr(i,9) - dpdr(i,8)
        dpdr(i,11) = dmdr(i,42) + dmdr(i,43) + dmdr(i,44)
        dpdr(i,12) = dmdr(i,45) + dmdr(i,46) + dmdr(i,47) + dmdr(i,48) + dmdr(i,49) + dmdr(i,50) + dmdr(i,51) + &
            & dmdr(i,52) + dmdr(i,53) + dmdr(i,54) + dmdr(i,55) + dmdr(i,56)
        dpdr(i,13) = dpdr(i,2)*p(3) + p(2)*dpdr(i,3) - dpdr(i,12)
        dpdr(i,14) = dpdr(i,1)*p(5) + p(1)*dpdr(i,5) - dpdr(i,12) - dpdr(i,11) - dpdr(i,13) - dpdr(i,12) - dpdr(i, &
            & 11) - dpdr(i,11) - dpdr(i,11)
        dpdr(i,15) = dpdr(i,1)*p(6) + p(1)*dpdr(i,6) - dpdr(i,12)
        dpdr(i,16) = dpdr(i,1)*p(7) + p(1)*dpdr(i,7) - dpdr(i,12)
        dpdr(i,17) = dpdr(i,2)*p(2) + p(2)*dpdr(i,2) - dpdr(i,11) - dpdr(i,11)
        dpdr(i,18) = dpdr(i,3)*p(3) + p(3)*dpdr(i,3) - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) - dpdr(i,16) - dpdr(i, &
            & 14) - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) - dpdr(i,16) - dpdr(i,14) - dpdr(i,12) - dpdr(i,11) - dpdr(i, &
            & 12) - dpdr(i,11)
        dpdr(i,19) = dpdr(i,2)*p(4) + p(2)*dpdr(i,4) - dpdr(i,14)
        dpdr(i,20) = dpdr(i,3)*p(4) + p(3)*dpdr(i,4) - dpdr(i,15) - dpdr(i,16) - dpdr(i,13)
        dpdr(i,21) = dpdr(i,1)*p(10) + p(1)*dpdr(i,10) - dpdr(i,20) - dpdr(i,19)
        dpdr(i,22) = dmdr(i,57) + dmdr(i,58) + dmdr(i,59) + dmdr(i,60) + dmdr(i,61) + dmdr(i,62)
        dpdr(i,23) = dpdr(i,1)*p(11) + p(1)*dpdr(i,11) - dpdr(i,22)
        dpdr(i,24) = dpdr(i,2)*p(6) + p(2)*dpdr(i,6)
        dpdr(i,25) = dpdr(i,2)*p(7) + p(2)*dpdr(i,7)
        dpdr(i,26) = dpdr(i,1)*p(12) + p(1)*dpdr(i,12) - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) - dpdr(i,22) - dpdr(i, &
            & 22) - dpdr(i,22)
        dpdr(i,27) = dpdr(i,2)*p(5) + p(2)*dpdr(i,5) - dpdr(i,22) - dpdr(i,23) - dpdr(i,22)
        dpdr(i,28) = dpdr(i,3)*p(5) + p(3)*dpdr(i,5) - dpdr(i,22) - dpdr(i,26) - dpdr(i,24) - dpdr(i,25) - dpdr(i, &
            & 23) - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) - dpdr(i,23) - dpdr(i,22) - dpdr(i,22)
        dpdr(i,29) = dpdr(i,3)*p(6) + p(3)*dpdr(i,6) - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
        dpdr(i,30) = dpdr(i,3)*p(7) + p(3)*dpdr(i,7) - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
        dpdr(i,31) = dpdr(i,2)*p(9) + p(2)*dpdr(i,9) - dpdr(i,26) - dpdr(i,28)
        dpdr(i,32) = dpdr(i,1)*p(14) + p(1)*dpdr(i,14) - dpdr(i,26) - dpdr(i,23) - dpdr(i,28)
        dpdr(i,33) = dpdr(i,4)*p(6) + p(4)*dpdr(i,6) - dpdr(i,25)
        dpdr(i,34) = dpdr(i,4)*p(7) + p(4)*dpdr(i,7) - dpdr(i,24)
        dpdr(i,35) = dpdr(i,1)*p(17) + p(1)*dpdr(i,17) - dpdr(i,27)
        dpdr(i,36) = dpdr(i,1)*p(18) + p(1)*dpdr(i,18) - dpdr(i,29) - dpdr(i,30) - dpdr(i,28)
        dpdr(i,37) = dpdr(i,2)*p(10) + p(2)*dpdr(i,10) - dpdr(i,32)
        dpdr(i,38) = dpdr(i,3)*p(10) + p(3)*dpdr(i,10) - dpdr(i,33) - dpdr(i,34) - dpdr(i,31)
        dpdr(i,39) = dpdr(i,1)*p(21) + p(1)*dpdr(i,21) - dpdr(i,38) - dpdr(i,37)
        dpdr(i,40) = dmdr(i,63)
        dpdr(i,41) = dmdr(i,64) + dmdr(i,65) + dmdr(i,66) + dmdr(i,67) + dmdr(i,68) + dmdr(i,69) + dmdr(i,70) + &
            & dmdr(i,71) + dmdr(i,72) + dmdr(i,73) + dmdr(i,74) + dmdr(i,75) + dmdr(i,76) + dmdr(i,77) + dmdr(i,78) + &
            & dmdr(i,79) + dmdr(i,80) + dmdr(i,81) + dmdr(i,82) + dmdr(i,83) + dmdr(i,84) + dmdr(i,85) + dmdr(i,86) + &
            & dmdr(i,87)
        dpdr(i,42) = dpdr(i,1)*p(22) + p(1)*dpdr(i,22) - dpdr(i,40) - dpdr(i,41) - dpdr(i,40) - dpdr(i,40) - dpdr(i, &
            & 40) - dpdr(i,40) - dpdr(i,40)
        dpdr(i,43) = dpdr(i,2)*p(11) + p(2)*dpdr(i,11) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
        dpdr(i,44) = dpdr(i,2)*p(12) + p(2)*dpdr(i,12) - dpdr(i,41)
        dpdr(i,45) = dpdr(i,3)*p(11) + p(3)*dpdr(i,11) - dpdr(i,41)
        dpdr(i,46) = dpdr(i,5)*p(6) + p(5)*dpdr(i,6) - dpdr(i,41)
        dpdr(i,47) = dpdr(i,5)*p(7) + p(5)*dpdr(i,7) - dpdr(i,41)
        dpdr(i,48) = dpdr(i,6)*p(7) + p(6)*dpdr(i,7) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
        dpdr(i,49) = dpdr(i,4)*p(11) + p(4)*dpdr(i,11) - dpdr(i,42)
        dpdr(i,50) = dpdr(i,2)*p(15) + p(2)*dpdr(i,15) - dpdr(i,46)
        dpdr(i,51) = dpdr(i,2)*p(16) + p(2)*dpdr(i,16) - dpdr(i,47)
        dpdr(i,52) = dpdr(i,4)*p(12) + p(4)*dpdr(i,12) - dpdr(i,41) - dpdr(i,50) - dpdr(i,51)
        dpdr(i,53) = dpdr(i,2)*p(14) + p(2)*dpdr(i,14) - dpdr(i,42) - dpdr(i,49) - dpdr(i,42)
        dpdr(i,54) = dpdr(i,6)*p(6) + p(6)*dpdr(i,6) - dpdr(i,42) - dpdr(i,42)
        dpdr(i,55) = dpdr(i,7)*p(7) + p(7)*dpdr(i,7) - dpdr(i,42) - dpdr(i,42)
        dpdr(i,56) = dpdr(i,3)*p(17) + p(3)*dpdr(i,17) - dpdr(i,44)
        dpdr(i,57) = dpdr(i,2)*p(18) + p(2)*dpdr(i,18) - dpdr(i,48)
        dpdr(i,58) = dpdr(i,3)*p(14) + p(3)*dpdr(i,14) - dpdr(i,41) - dpdr(i,52) - dpdr(i,46) - dpdr(i,47) - dpdr(i, &
            & 45) - dpdr(i,45)
        dpdr(i,59) = dpdr(i,6)*p(9) + p(6)*dpdr(i,9) - dpdr(i,41) - dpdr(i,52) - dpdr(i,47)
        dpdr(i,60) = dpdr(i,7)*p(9) + p(7)*dpdr(i,9) - dpdr(i,41) - dpdr(i,52) - dpdr(i,46)
        dpdr(i,61) = dpdr(i,2)*p(20) + p(2)*dpdr(i,20) - dpdr(i,52) - dpdr(i,58)
        dpdr(i,62) = dpdr(i,1)*p(32) + p(1)*dpdr(i,32) - dpdr(i,52) - dpdr(i,49) - dpdr(i,58)
        dpdr(i,63) = dpdr(i,6)*p(10) + p(6)*dpdr(i,10) - dpdr(i,51)
        dpdr(i,64) = dpdr(i,7)*p(10) + p(7)*dpdr(i,10) - dpdr(i,50)
        dpdr(i,65) = dpdr(i,2)*p(17) + p(2)*dpdr(i,17) - dpdr(i,43)
        dpdr(i,66) = dpdr(i,3)*p(18) + p(3)*dpdr(i,18) - dpdr(i,46) - dpdr(i,47) - dpdr(i,45) - dpdr(i,59) - dpdr(i, &
            & 60) - dpdr(i,58)
        dpdr(i,67) = dpdr(i,2)*p(19) + p(2)*dpdr(i,19) - dpdr(i,49)
        dpdr(i,68) = dpdr(i,1)*p(36) + p(1)*dpdr(i,36) - dpdr(i,59) - dpdr(i,60) - dpdr(i,58) - dpdr(i,57) - dpdr(i, &
            & 66) - dpdr(i,66)
        dpdr(i,69) = dpdr(i,2)*p(21) + p(2)*dpdr(i,21) - dpdr(i,62)
        dpdr(i,70) = dpdr(i,3)*p(21) + p(3)*dpdr(i,21) - dpdr(i,63) - dpdr(i,64) - dpdr(i,61)
        dpdr(i,71) = dpdr(i,1)*p(39) + p(1)*dpdr(i,39) - dpdr(i,70) - dpdr(i,69)
        dpdr(i,72) = dpdr(i,40)*p(1) + p(40)*dpdr(i,1)
        dpdr(i,73) = dpdr(i,2)*p(22) + p(2)*dpdr(i,22) - dpdr(i,72)
        dpdr(i,74) = dpdr(i,6)*p(11) + p(6)*dpdr(i,11)
        dpdr(i,75) = dpdr(i,7)*p(11) + p(7)*dpdr(i,11)
        dpdr(i,76) = dpdr(i,3)*p(22) + p(3)*dpdr(i,22) - dpdr(i,72) - dpdr(i,74) - dpdr(i,75) - dpdr(i,72) - dpdr(i, &
            & 72) - dpdr(i,72)
        dpdr(i,77) = dmdr(i,88) + dmdr(i,89) + dmdr(i,90) + dmdr(i,91) + dmdr(i,92) + dmdr(i,93) + dmdr(i,94) + &
            & dmdr(i,95) + dmdr(i,96) + dmdr(i,97) + dmdr(i,98) + dmdr(i,99) + dmdr(i,100) + dmdr(i,101) + dmdr(i, &
            & 102) + dmdr(i,103) + dmdr(i,104) + dmdr(i,105) + dmdr(i,106) + dmdr(i,107) + dmdr(i,108) + dmdr(i, &
            & 109) + dmdr(i,110) + dmdr(i,111)
        dpdr(i,78)  = dpdr(i,1)*p(42) + p(1)*dpdr(i,42) - dpdr(i,72) - dpdr(i,76)
        dpdr(i,79)  = dpdr(i,5)*p(11) + p(5)*dpdr(i,11) - dpdr(i,72) - dpdr(i,73) - dpdr(i,72)
        dpdr(i,80)  = dpdr(i,2)*p(26) + p(2)*dpdr(i,26) - dpdr(i,76) - dpdr(i,77)
        dpdr(i,81)  = dpdr(i,6)*p(12) + p(6)*dpdr(i,12) - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
        dpdr(i,82)  = dpdr(i,7)*p(12) + p(7)*dpdr(i,12) - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
        dpdr(i,83)  = dpdr(i,8)*p(11) + p(8)*dpdr(i,11) - dpdr(i,72)
        dpdr(i,84)  = dpdr(i,6)*p(17) + p(6)*dpdr(i,17)
        dpdr(i,85)  = dpdr(i,7)*p(17) + p(7)*dpdr(i,17)
        dpdr(i,86)  = dpdr(i,9)*p(11) + p(9)*dpdr(i,11) - dpdr(i,76) - dpdr(i,77)
        dpdr(i,87)  = dpdr(i,2)*p(29) + p(2)*dpdr(i,29) - dpdr(i,81)
        dpdr(i,88)  = dpdr(i,6)*p(14) + p(6)*dpdr(i,14) - dpdr(i,75) - dpdr(i,75)
        dpdr(i,89)  = dpdr(i,2)*p(30) + p(2)*dpdr(i,30) - dpdr(i,82)
        dpdr(i,90)  = dpdr(i,7)*p(14) + p(7)*dpdr(i,14) - dpdr(i,74) - dpdr(i,74)
        dpdr(i,91)  = dpdr(i,1)*p(48) + p(1)*dpdr(i,48) - dpdr(i,76) - dpdr(i,81) - dpdr(i,82)
        dpdr(i,92)  = dpdr(i,10)*p(11) + p(10)*dpdr(i,11) - dpdr(i,78)
        dpdr(i,93)  = dpdr(i,2)*p(33) + p(2)*dpdr(i,33) - dpdr(i,88)
        dpdr(i,94)  = dpdr(i,2)*p(34) + p(2)*dpdr(i,34) - dpdr(i,90)
        dpdr(i,95)  = dpdr(i,10)*p(12) + p(10)*dpdr(i,12) - dpdr(i,77) - dpdr(i,93) - dpdr(i,94)
        dpdr(i,96)  = dpdr(i,2)*p(27) + p(2)*dpdr(i,27) - dpdr(i,73) - dpdr(i,79)
        dpdr(i,97)  = dpdr(i,2)*p(28) + p(2)*dpdr(i,28) - dpdr(i,76) - dpdr(i,86)
        dpdr(i,98)  = dpdr(i,1)*p(53) + p(1)*dpdr(i,53) - dpdr(i,80) - dpdr(i,79) - dpdr(i,97)
        dpdr(i,99)  = dpdr(i,1)*p(54) + p(1)*dpdr(i,54) - dpdr(i,81)
        dpdr(i,100) = dpdr(i,1)*p(55) + p(1)*dpdr(i,55) - dpdr(i,82)
        dpdr(i,101) = dpdr(i,5)*p(18) + p(5)*dpdr(i,18) - dpdr(i,76) - dpdr(i,91) - dpdr(i,87) - dpdr(i,89) - dpdr(i,86)
        dpdr(i,102) = dpdr(i,6)*p(18) + p(6)*dpdr(i,18) - dpdr(i,74) - dpdr(i,90)
        dpdr(i,103) = dpdr(i,7)*p(18) + p(7)*dpdr(i,18) - dpdr(i,75) - dpdr(i,88)
        dpdr(i,104) = dpdr(i,2)*p(31) + p(2)*dpdr(i,31) - dpdr(i,77) - dpdr(i,86)
        dpdr(i,105) = dpdr(i,2)*p(36) + p(2)*dpdr(i,36) - dpdr(i,91) - dpdr(i,101)
        dpdr(i,106) = dpdr(i,3)*p(32) + p(3)*dpdr(i,32) - dpdr(i,77) - dpdr(i,95) - dpdr(i,88) - dpdr(i,90) - dpdr(i,86)
        dpdr(i,107) = dpdr(i,4)*p(29) + p(4)*dpdr(i,29) - dpdr(i,82) - dpdr(i,80) - dpdr(i,99)
        dpdr(i,108) = dpdr(i,4)*p(30) + p(4)*dpdr(i,30) - dpdr(i,81) - dpdr(i,80) - dpdr(i,100)
        dpdr(i,109) = dpdr(i,2)*p(38) + p(2)*dpdr(i,38) - dpdr(i,95) - dpdr(i,106)
        dpdr(i,110) = dpdr(i,1)*p(62) + p(1)*dpdr(i,62) - dpdr(i,95) - dpdr(i,92) - dpdr(i,106)
        dpdr(i,111) = dpdr(i,6)*p(21) + p(6)*dpdr(i,21) - dpdr(i,94)
        dpdr(i,112) = dpdr(i,7)*p(21) + p(7)*dpdr(i,21) - dpdr(i,93)
        dpdr(i,113) = dpdr(i,1)*p(65) + p(1)*dpdr(i,65) - dpdr(i,96)
        dpdr(i,114) = dpdr(i,1)*p(66) + p(1)*dpdr(i,66) - dpdr(i,102) - dpdr(i,103) - dpdr(i,101)
        dpdr(i,115) = dpdr(i,2)*p(37) + p(2)*dpdr(i,37) - dpdr(i,92)
        dpdr(i,116) = dpdr(i,10)*p(18) + p(10)*dpdr(i,18) - dpdr(i,99) - dpdr(i,100) - dpdr(i,97)
        dpdr(i,117) = dpdr(i,2)*p(39) + p(2)*dpdr(i,39) - dpdr(i,110)
        dpdr(i,118) = dpdr(i,3)*p(39) + p(3)*dpdr(i,39) - dpdr(i,111) - dpdr(i,112) - dpdr(i,109)
        dpdr(i,119) = dpdr(i,1)*p(71) + p(1)*dpdr(i,71) - dpdr(i,118) - dpdr(i,117)
        dpdr(i,120) = dpdr(i,40)*p(2) + p(40)*dpdr(i,2)
        dpdr(i,121) = dpdr(i,40)*p(3) + p(40)*dpdr(i,3)
        dpdr(i,122) = dpdr(i,40)*p(4) + p(40)*dpdr(i,4)
        dpdr(i,123) = dpdr(i,11)*p(12) + p(11)*dpdr(i,12) - dpdr(i,121)
        dpdr(i,124) = dpdr(i,2)*p(42) + p(2)*dpdr(i,42) - dpdr(i,122)
        dpdr(i,125) = dpdr(i,6)*p(22) + p(6)*dpdr(i,22) - dpdr(i,121)
        dpdr(i,126) = dpdr(i,7)*p(22) + p(7)*dpdr(i,22) - dpdr(i,121)
        dpdr(i,127) = dpdr(i,2)*p(41) + p(2)*dpdr(i,41) - dpdr(i,121) - dpdr(i,123) - dpdr(i,121)
        dpdr(i,128) = dpdr(i,6)*p(23) + p(6)*dpdr(i,23) - dpdr(i,123)
        dpdr(i,129) = dpdr(i,7)*p(23) + p(7)*dpdr(i,23) - dpdr(i,123)
        dpdr(i,130) = dpdr(i,3)*p(41) + p(3)*dpdr(i,41) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,125) - &
            & dpdr(i,128) - dpdr(i,126) - dpdr(i,129) - dpdr(i,124) - dpdr(i,123) - dpdr(i,122) - dpdr(i,121) - &
            & dpdr(i,120) - dpdr(i,125) - dpdr(i,126) - dpdr(i,124) - dpdr(i,123) - dpdr(i,122) - dpdr(i,121) - &
            & dpdr(i,120) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120) - &
            & dpdr(i,120)
        dpdr(i,131) = dpdr(i,3)*p(42) + p(3)*dpdr(i,42) - dpdr(i,121) - dpdr(i,125) - dpdr(i,126) - dpdr(i,121)
        dpdr(i,132) = dpdr(i,4)*p(41) + p(4)*dpdr(i,41) - dpdr(i,121) - dpdr(i,131) - dpdr(i,128) - dpdr(i,129) - &
            & dpdr(i,127) - dpdr(i,121)
        dpdr(i,133) = dpdr(i,1)*p(78) + p(1)*dpdr(i,78) - dpdr(i,122) - dpdr(i,131)
        dpdr(i,134) = dpdr(i,11)*p(11) + p(11)*dpdr(i,11) - dpdr(i,120) - dpdr(i,120)
        dpdr(i,135) = dpdr(i,2)*p(48) + p(2)*dpdr(i,48) - dpdr(i,130)
        dpdr(i,136) = dpdr(i,11)*p(17) + p(11)*dpdr(i,17) - dpdr(i,120)
        dpdr(i,137) = dpdr(i,2)*p(44) + p(2)*dpdr(i,44) - dpdr(i,123)
        dpdr(i,138) = dpdr(i,2)*p(45) + p(2)*dpdr(i,45) - dpdr(i,121)
        dpdr(i,139) = dpdr(i,11)*p(14) + p(11)*dpdr(i,14) - dpdr(i,122) - dpdr(i,124) - dpdr(i,122)
        dpdr(i,140) = dpdr(i,6)*p(27) + p(6)*dpdr(i,27) - dpdr(i,127)
        dpdr(i,141) = dpdr(i,2)*p(54) + p(2)*dpdr(i,54)
        dpdr(i,142) = dpdr(i,7)*p(27) + p(7)*dpdr(i,27) - dpdr(i,127)
        dpdr(i,143) = dpdr(i,2)*p(52) + p(2)*dpdr(i,52) - dpdr(i,131) - dpdr(i,132)
        dpdr(i,144) = dpdr(i,1)*p(81) + p(1)*dpdr(i,81) - dpdr(i,125) - dpdr(i,141) - dpdr(i,135) - dpdr(i,125)
        dpdr(i,145) = dpdr(i,2)*p(55) + p(2)*dpdr(i,55)
        dpdr(i,146) = dpdr(i,1)*p(82) + p(1)*dpdr(i,82) - dpdr(i,126) - dpdr(i,135) - dpdr(i,145) - dpdr(i,126)
        dpdr(i,147) = dpdr(i,11)*p(18) + p(11)*dpdr(i,18) - dpdr(i,130)
        dpdr(i,148) = dpdr(i,6)*p(28) + p(6)*dpdr(i,28) - dpdr(i,129) - dpdr(i,123) - dpdr(i,143)
        dpdr(i,149) = dpdr(i,7)*p(28) + p(7)*dpdr(i,28) - dpdr(i,128) - dpdr(i,123) - dpdr(i,143)
        dpdr(i,150) = dpdr(i,6)*p(30) + p(6)*dpdr(i,30) - dpdr(i,121) - dpdr(i,146)
        dpdr(i,151) = dpdr(i,11)*p(19) + p(11)*dpdr(i,19) - dpdr(i,122)
        dpdr(i,152) = dpdr(i,2)*p(50) + p(2)*dpdr(i,50) - dpdr(i,128)
        dpdr(i,153) = dpdr(i,2)*p(51) + p(2)*dpdr(i,51) - dpdr(i,129)
        dpdr(i,154) = dpdr(i,11)*p(20) + p(11)*dpdr(i,20) - dpdr(i,131) - dpdr(i,132)
        dpdr(i,155) = dpdr(i,2)*p(59) + p(2)*dpdr(i,59) - dpdr(i,144) - dpdr(i,148)
        dpdr(i,156) = dpdr(i,6)*p(32) + p(6)*dpdr(i,32) - dpdr(i,129)
        dpdr(i,157) = dpdr(i,2)*p(60) + p(2)*dpdr(i,60) - dpdr(i,146) - dpdr(i,149)
        dpdr(i,158) = dpdr(i,7)*p(32) + p(7)*dpdr(i,32) - dpdr(i,128)
        dpdr(i,159) = dpdr(i,6)*p(34) + p(6)*dpdr(i,34) - dpdr(i,122) - dpdr(i,145) - dpdr(i,122)
        dpdr(i,160) = dpdr(i,11)*p(21) + p(11)*dpdr(i,21) - dpdr(i,133)
        dpdr(i,161) = dpdr(i,2)*p(63) + p(2)*dpdr(i,63) - dpdr(i,156)
        dpdr(i,162) = dpdr(i,2)*p(64) + p(2)*dpdr(i,64) - dpdr(i,158)
        dpdr(i,163) = dpdr(i,12)*p(21) + p(12)*dpdr(i,21) - dpdr(i,132) - dpdr(i,161) - dpdr(i,162)
        dpdr(i,164) = dpdr(i,2)*p(53) + p(2)*dpdr(i,53) - dpdr(i,124) - dpdr(i,139)
        dpdr(i,165) = dpdr(i,2)*p(58) + p(2)*dpdr(i,58) - dpdr(i,131) - dpdr(i,154)
        dpdr(i,166) = dpdr(i,3)*p(54) + p(3)*dpdr(i,54) - dpdr(i,125) - dpdr(i,144)
        dpdr(i,167) = dpdr(i,3)*p(55) + p(3)*dpdr(i,55) - dpdr(i,126) - dpdr(i,146)
        dpdr(i,168) = dpdr(i,3)*p(65) + p(3)*dpdr(i,65) - dpdr(i,137)
        dpdr(i,169) = dpdr(i,17)*p(18) + p(17)*dpdr(i,18) - dpdr(i,135)
        dpdr(i,170) = dpdr(i,1)*p(98) + p(1)*dpdr(i,98) - dpdr(i,143) - dpdr(i,139) - dpdr(i,165)
        dpdr(i,171) = dpdr(i,4)*p(54) + p(4)*dpdr(i,54) - dpdr(i,135)
        dpdr(i,172) = dpdr(i,4)*p(55) + p(4)*dpdr(i,55) - dpdr(i,135)
        dpdr(i,173) = dpdr(i,2)*p(66) + p(2)*dpdr(i,66) - dpdr(i,150)
        dpdr(i,174) = dpdr(i,1)*p(101) + p(1)*dpdr(i,101) - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,173) - &
            & dpdr(i,165) - dpdr(i,147)
        dpdr(i,175) = dpdr(i,1)*p(102) + p(1)*dpdr(i,102) - dpdr(i,150) - dpdr(i,148) - dpdr(i,166)
        dpdr(i,176) = dpdr(i,1)*p(103) + p(1)*dpdr(i,103) - dpdr(i,150) - dpdr(i,149) - dpdr(i,167)
        dpdr(i,177) = dpdr(i,2)*p(61) + p(2)*dpdr(i,61) - dpdr(i,132) - dpdr(i,154)
        dpdr(i,178) = dpdr(i,2)*p(68) + p(2)*dpdr(i,68) - dpdr(i,159) - dpdr(i,174)
        dpdr(i,179) = dpdr(i,3)*p(62) + p(3)*dpdr(i,62) - dpdr(i,132) - dpdr(i,163) - dpdr(i,156) - dpdr(i,158) - &
            & dpdr(i,154)
        dpdr(i,180) = dpdr(i,6)*p(38) + p(6)*dpdr(i,38) - dpdr(i,132) - dpdr(i,163) - dpdr(i,157)
        dpdr(i,181) = dpdr(i,7)*p(38) + p(7)*dpdr(i,38) - dpdr(i,132) - dpdr(i,163) - dpdr(i,155)
        dpdr(i,182) = dpdr(i,2)*p(70) + p(2)*dpdr(i,70) - dpdr(i,163) - dpdr(i,179)
        dpdr(i,183) = dpdr(i,1)*p(110) + p(1)*dpdr(i,110) - dpdr(i,163) - dpdr(i,160) - dpdr(i,179)
        dpdr(i,184) = dpdr(i,6)*p(39) + p(6)*dpdr(i,39) - dpdr(i,162)
        dpdr(i,185) = dpdr(i,7)*p(39) + p(7)*dpdr(i,39) - dpdr(i,161)
        dpdr(i,186) = dpdr(i,2)*p(65) + p(2)*dpdr(i,65) - dpdr(i,136)
        dpdr(i,187) = dpdr(i,3)*p(66) + p(3)*dpdr(i,66) - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,175) - &
            & dpdr(i,176) - dpdr(i,174)
        dpdr(i,188) = dpdr(i,2)*p(67) + p(2)*dpdr(i,67) - dpdr(i,151)
        dpdr(i,189) = dpdr(i,4)*p(66) + p(4)*dpdr(i,66) - dpdr(i,166) - dpdr(i,167) - dpdr(i,165)
        dpdr(i,190) = dpdr(i,2)*p(69) + p(2)*dpdr(i,69) - dpdr(i,160)
        dpdr(i,191) = dpdr(i,18)*p(21) + p(18)*dpdr(i,21) - dpdr(i,171) - dpdr(i,172) - dpdr(i,169)
        dpdr(i,192) = dpdr(i,2)*p(71) + p(2)*dpdr(i,71) - dpdr(i,183)
        dpdr(i,193) = dpdr(i,3)*p(71) + p(3)*dpdr(i,71) - dpdr(i,184) - dpdr(i,185) - dpdr(i,182)
        dpdr(i,194) = dpdr(i,1)*p(119) + p(1)*dpdr(i,119) - dpdr(i,193) - dpdr(i,192)
        dpdr(i,195) = dpdr(i,40)*p(5) + p(40)*dpdr(i,5)
        dpdr(i,196) = dpdr(i,40)*p(6) + p(40)*dpdr(i,6)
        dpdr(i,197) = dpdr(i,40)*p(7) + p(40)*dpdr(i,7)
        dpdr(i,198) = dpdr(i,40)*p(8) + p(40)*dpdr(i,8)
        dpdr(i,199) = dpdr(i,40)*p(9) + p(40)*dpdr(i,9)
        dpdr(i,200) = dpdr(i,40)*p(10) + p(40)*dpdr(i,10)
        dpdr(i,201) = dpdr(i,11)*p(22) + p(11)*dpdr(i,22) - dpdr(i,195)
        dpdr(i,202) = dpdr(i,12)*p(22) + p(12)*dpdr(i,22) - dpdr(i,196) - dpdr(i,197) - dpdr(i,195) - dpdr(i,196) - &
            & dpdr(i,197) - dpdr(i,195) - dpdr(i,196) - dpdr(i,197)
        dpdr(i,203) = dpdr(i,17)*p(22) + p(17)*dpdr(i,22) - dpdr(i,198)
        dpdr(i,204) = dpdr(i,6)*p(43) + p(6)*dpdr(i,43)
        dpdr(i,205) = dpdr(i,7)*p(43) + p(7)*dpdr(i,43)
        dpdr(i,206) = dpdr(i,11)*p(26) + p(11)*dpdr(i,26) - dpdr(i,199) - dpdr(i,202)
        dpdr(i,207) = dpdr(i,2)*p(76) + p(2)*dpdr(i,76) - dpdr(i,199) - dpdr(i,202)
        dpdr(i,208) = dpdr(i,2)*p(78) + p(2)*dpdr(i,78) - dpdr(i,200)
        dpdr(i,209) = dpdr(i,6)*p(41) + p(6)*dpdr(i,41) - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
        dpdr(i,210) = dpdr(i,6)*p(42) + p(6)*dpdr(i,42) - dpdr(i,197) - dpdr(i,197) - dpdr(i,197)
        dpdr(i,211) = dpdr(i,7)*p(41) + p(7)*dpdr(i,41) - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
        dpdr(i,212) = dpdr(i,7)*p(42) + p(7)*dpdr(i,42) - dpdr(i,196) - dpdr(i,196) - dpdr(i,196)
        dpdr(i,213) = dpdr(i,11)*p(29) + p(11)*dpdr(i,29) - dpdr(i,209)
        dpdr(i,214) = dpdr(i,11)*p(30) + p(11)*dpdr(i,30) - dpdr(i,211)
        dpdr(i,215) = dpdr(i,18)*p(22) + p(18)*dpdr(i,22) - dpdr(i,199) - dpdr(i,213) - dpdr(i,214)
        dpdr(i,216) = dpdr(i,2)*p(77) + p(2)*dpdr(i,77) - dpdr(i,199) - dpdr(i,206)
        dpdr(i,217) = dpdr(i,6)*p(49) + p(6)*dpdr(i,49) - dpdr(i,205)
        dpdr(i,218) = dpdr(i,7)*p(49) + p(7)*dpdr(i,49) - dpdr(i,204)
        dpdr(i,219) = dpdr(i,3)*p(77) + p(3)*dpdr(i,77) - dpdr(i,200) - dpdr(i,199) - dpdr(i,198) - dpdr(i,209) - &
            & dpdr(i,217) - dpdr(i,211) - dpdr(i,218) - dpdr(i,207) - dpdr(i,204) - dpdr(i,205) - dpdr(i,200) - &
            & dpdr(i,199) - dpdr(i,198) - dpdr(i,200) - dpdr(i,198) - dpdr(i,200) - dpdr(i,198)
        dpdr(i,220) = dpdr(i,3)*p(78) + p(3)*dpdr(i,78) - dpdr(i,199) - dpdr(i,210) - dpdr(i,212)
        dpdr(i,221) = dpdr(i,10)*p(41) + p(10)*dpdr(i,41) - dpdr(i,199) - dpdr(i,220) - dpdr(i,217) - dpdr(i,218) - &
            & dpdr(i,216)
        dpdr(i,222) = dpdr(i,1)*p(133) + p(1)*dpdr(i,133) - dpdr(i,200) - dpdr(i,220)
        dpdr(i,223) = dpdr(i,11)*p(27) + p(11)*dpdr(i,27) - dpdr(i,195) - dpdr(i,203)
        dpdr(i,224) = dpdr(i,1)*p(134) + p(1)*dpdr(i,134) - dpdr(i,201)
        dpdr(i,225) = dpdr(i,2)*p(81) + p(2)*dpdr(i,81) - dpdr(i,209)
        dpdr(i,226) = dpdr(i,2)*p(80) + p(2)*dpdr(i,80) - dpdr(i,202) - dpdr(i,206)
        dpdr(i,227) = dpdr(i,2)*p(82) + p(2)*dpdr(i,82) - dpdr(i,211)
        dpdr(i,228) = dpdr(i,2)*p(91) + p(2)*dpdr(i,91) - dpdr(i,215) - dpdr(i,219)
        dpdr(i,229) = dpdr(i,11)*p(28) + p(11)*dpdr(i,28) - dpdr(i,199) - dpdr(i,207)
        dpdr(i,230) = dpdr(i,6)*p(53) + p(6)*dpdr(i,53) - dpdr(i,205)
        dpdr(i,231) = dpdr(i,5)*p(54) + p(5)*dpdr(i,54) - dpdr(i,209)
        dpdr(i,232) = dpdr(i,7)*p(53) + p(7)*dpdr(i,53) - dpdr(i,204)
        dpdr(i,233) = dpdr(i,7)*p(54) + p(7)*dpdr(i,54) - dpdr(i,196)
        dpdr(i,234) = dpdr(i,5)*p(55) + p(5)*dpdr(i,55) - dpdr(i,211)
        dpdr(i,235) = dpdr(i,6)*p(55) + p(6)*dpdr(i,55) - dpdr(i,197)
        dpdr(i,236) = dpdr(i,11)*p(35) + p(11)*dpdr(i,35) - dpdr(i,198)
        dpdr(i,237) = dpdr(i,6)*p(65) + p(6)*dpdr(i,65)
        dpdr(i,238) = dpdr(i,7)*p(65) + p(7)*dpdr(i,65)
        dpdr(i,239) = dpdr(i,2)*p(86) + p(2)*dpdr(i,86) - dpdr(i,199) - dpdr(i,229)
        dpdr(i,240) = dpdr(i,1)*p(139) + p(1)*dpdr(i,139) - dpdr(i,206) - dpdr(i,229) - dpdr(i,224)
        dpdr(i,241) = dpdr(i,17)*p(29) + p(17)*dpdr(i,29) - dpdr(i,225)
        dpdr(i,242) = dpdr(i,2)*p(99) + p(2)*dpdr(i,99) - dpdr(i,231)
        dpdr(i,243) = dpdr(i,17)*p(30) + p(17)*dpdr(i,30) - dpdr(i,227)
        dpdr(i,244) = dpdr(i,2)*p(95) + p(2)*dpdr(i,95) - dpdr(i,220) - dpdr(i,221)
        dpdr(i,245) = dpdr(i,4)*p(81) + p(4)*dpdr(i,81) - dpdr(i,202) - dpdr(i,242) - dpdr(i,227)
        dpdr(i,246) = dpdr(i,2)*p(100) + p(2)*dpdr(i,100) - dpdr(i,234)
        dpdr(i,247) = dpdr(i,4)*p(82) + p(4)*dpdr(i,82) - dpdr(i,202) - dpdr(i,225) - dpdr(i,246)
        dpdr(i,248) = dpdr(i,11)*p(36) + p(11)*dpdr(i,36) - dpdr(i,215) - dpdr(i,219)
        dpdr(i,249) = dpdr(i,2)*p(102) + p(2)*dpdr(i,102) - dpdr(i,233)
        dpdr(i,250) = dpdr(i,6)*p(58) + p(6)*dpdr(i,58) - dpdr(i,206) - dpdr(i,214) - dpdr(i,244) - dpdr(i,214)
        dpdr(i,251) = dpdr(i,2)*p(103) + p(2)*dpdr(i,103) - dpdr(i,235)
        dpdr(i,252) = dpdr(i,7)*p(58) + p(7)*dpdr(i,58) - dpdr(i,213) - dpdr(i,206) - dpdr(i,244) - dpdr(i,213)
        dpdr(i,253) = dpdr(i,1)*p(150) + p(1)*dpdr(i,150) - dpdr(i,215) - dpdr(i,233) - dpdr(i,235)
        dpdr(i,254) = dpdr(i,11)*p(37) + p(11)*dpdr(i,37) - dpdr(i,200)
        dpdr(i,255) = dpdr(i,2)*p(93) + p(2)*dpdr(i,93) - dpdr(i,217)
        dpdr(i,256) = dpdr(i,2)*p(94) + p(2)*dpdr(i,94) - dpdr(i,218)
        dpdr(i,257) = dpdr(i,11)*p(38) + p(11)*dpdr(i,38) - dpdr(i,220) - dpdr(i,221)
        dpdr(i,258) = dpdr(i,2)*p(107) + p(2)*dpdr(i,107) - dpdr(i,245) - dpdr(i,250)
        dpdr(i,259) = dpdr(i,6)*p(62) + p(6)*dpdr(i,62) - dpdr(i,218)
        dpdr(i,260) = dpdr(i,2)*p(108) + p(2)*dpdr(i,108) - dpdr(i,247) - dpdr(i,252)
        dpdr(i,261) = dpdr(i,7)*p(62) + p(7)*dpdr(i,62) - dpdr(i,217)
        dpdr(i,262) = dpdr(i,6)*p(64) + p(6)*dpdr(i,64) - dpdr(i,200) - dpdr(i,246) - dpdr(i,200)
        dpdr(i,263) = dpdr(i,11)*p(39) + p(11)*dpdr(i,39) - dpdr(i,222)
        dpdr(i,264) = dpdr(i,2)*p(111) + p(2)*dpdr(i,111) - dpdr(i,259)
        dpdr(i,265) = dpdr(i,2)*p(112) + p(2)*dpdr(i,112) - dpdr(i,261)
        dpdr(i,266) = dpdr(i,12)*p(39) + p(12)*dpdr(i,39) - dpdr(i,221) - dpdr(i,264) - dpdr(i,265)
        dpdr(i,267) = dpdr(i,2)*p(98) + p(2)*dpdr(i,98) - dpdr(i,208) - dpdr(i,240)
        dpdr(i,268) = dpdr(i,6)*p(54) + p(6)*dpdr(i,54) - dpdr(i,210)
        dpdr(i,269) = dpdr(i,7)*p(55) + p(7)*dpdr(i,55) - dpdr(i,212)
        dpdr(i,270) = dpdr(i,2)*p(96) + p(2)*dpdr(i,96) - dpdr(i,203) - dpdr(i,223)
        dpdr(i,271) = dpdr(i,2)*p(97) + p(2)*dpdr(i,97) - dpdr(i,207) - dpdr(i,229)
        dpdr(i,272) = dpdr(i,2)*p(101) + p(2)*dpdr(i,101) - dpdr(i,215) - dpdr(i,248)
        dpdr(i,273) = dpdr(i,2)*p(106) + p(2)*dpdr(i,106) - dpdr(i,220) - dpdr(i,257)
        dpdr(i,274) = dpdr(i,6)*p(59) + p(6)*dpdr(i,59) - dpdr(i,220) - dpdr(i,215) - dpdr(i,209)
        dpdr(i,275) = dpdr(i,7)*p(60) + p(7)*dpdr(i,60) - dpdr(i,220) - dpdr(i,215) - dpdr(i,211)
        dpdr(i,276) = dpdr(i,5)*p(66) + p(5)*dpdr(i,66) - dpdr(i,215) - dpdr(i,253) - dpdr(i,249) - dpdr(i,251) - &
            & dpdr(i,248)
        dpdr(i,277) = dpdr(i,6)*p(66) + p(6)*dpdr(i,66) - dpdr(i,213) - dpdr(i,252)
        dpdr(i,278) = dpdr(i,7)*p(66) + p(7)*dpdr(i,66) - dpdr(i,214) - dpdr(i,250)
        dpdr(i,279) = dpdr(i,2)*p(104) + p(2)*dpdr(i,104) - dpdr(i,216) - dpdr(i,239)
        dpdr(i,280) = dpdr(i,2)*p(105) + p(2)*dpdr(i,105) - dpdr(i,219) - dpdr(i,248)
        dpdr(i,281) = dpdr(i,1)*p(170) + p(1)*dpdr(i,170) - dpdr(i,244) - dpdr(i,240) - dpdr(i,273)
        dpdr(i,282) = dpdr(i,10)*p(54) + p(10)*dpdr(i,54) - dpdr(i,227)
        dpdr(i,283) = dpdr(i,10)*p(55) + p(10)*dpdr(i,55) - dpdr(i,225)
        dpdr(i,284) = dpdr(i,2)*p(114) + p(2)*dpdr(i,114) - dpdr(i,253) - dpdr(i,276)
        dpdr(i,285) = dpdr(i,1)*p(174) + p(1)*dpdr(i,174) - dpdr(i,250) - dpdr(i,252) - dpdr(i,248) - dpdr(i,276) - &
            & dpdr(i,273)
        dpdr(i,286) = dpdr(i,6)*p(68) + p(6)*dpdr(i,68) - dpdr(i,217) - dpdr(i,261) - dpdr(i,251)
        dpdr(i,287) = dpdr(i,7)*p(68) + p(7)*dpdr(i,68) - dpdr(i,218) - dpdr(i,259) - dpdr(i,249)
        dpdr(i,288) = dpdr(i,2)*p(109) + p(2)*dpdr(i,109) - dpdr(i,221) - dpdr(i,257)
        dpdr(i,289) = dpdr(i,2)*p(116) + p(2)*dpdr(i,116) - dpdr(i,262) - dpdr(i,285)
        dpdr(i,290) = dpdr(i,3)*p(110) + p(3)*dpdr(i,110) - dpdr(i,221) - dpdr(i,266) - dpdr(i,259) - dpdr(i,261) - &
            & dpdr(i,257)
        dpdr(i,291) = dpdr(i,6)*p(70) + p(6)*dpdr(i,70) - dpdr(i,221) - dpdr(i,266) - dpdr(i,260)
        dpdr(i,292) = dpdr(i,7)*p(70) + p(7)*dpdr(i,70) - dpdr(i,221) - dpdr(i,266) - dpdr(i,258)
        dpdr(i,293) = dpdr(i,2)*p(118) + p(2)*dpdr(i,118) - dpdr(i,266) - dpdr(i,290)
        dpdr(i,294) = dpdr(i,1)*p(183) + p(1)*dpdr(i,183) - dpdr(i,266) - dpdr(i,263) - dpdr(i,290)
        dpdr(i,295) = dpdr(i,6)*p(71) + p(6)*dpdr(i,71) - dpdr(i,265)
        dpdr(i,296) = dpdr(i,7)*p(71) + p(7)*dpdr(i,71) - dpdr(i,264)
        dpdr(i,297) = dpdr(i,1)*p(186) + p(1)*dpdr(i,186) - dpdr(i,270)
        dpdr(i,298) = dpdr(i,1)*p(187) + p(1)*dpdr(i,187) - dpdr(i,277) - dpdr(i,278) - dpdr(i,276)
        dpdr(i,299) = dpdr(i,2)*p(115) + p(2)*dpdr(i,115) - dpdr(i,254)
        dpdr(i,300) = dpdr(i,1)*p(189) + p(1)*dpdr(i,189) - dpdr(i,286) - dpdr(i,287) - dpdr(i,285) - dpdr(i,284) - &
            & dpdr(i,298)
        dpdr(i,301) = dpdr(i,2)*p(117) + p(2)*dpdr(i,117) - dpdr(i,263)
        dpdr(i,302) = dpdr(i,18)*p(39) + p(18)*dpdr(i,39) - dpdr(i,282) - dpdr(i,283) - dpdr(i,280)
        dpdr(i,303) = dpdr(i,2)*p(119) + p(2)*dpdr(i,119) - dpdr(i,294)
        dpdr(i,304) = dpdr(i,3)*p(119) + p(3)*dpdr(i,119) - dpdr(i,295) - dpdr(i,296) - dpdr(i,293)
        dpdr(i,305) = dpdr(i,1)*p(194) + p(1)*dpdr(i,194) - dpdr(i,304) - dpdr(i,303)
    enddo
end subroutine evdpdr

!> @brief Evaluate evdbdr.
subroutine evdbdr
    use n4_1a_meg2_2020_par, only: dbdr, dpdr
    !  The subroutine eliminates the 2-body terms in Bowman's approach
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: i, j
    real(wp) :: db1dr(6,306)
    ! Pass P(0:305) to BM1(1:306)
    do j=1, 6
        do i=1, 306
            db1dr(j,i) = dpdr(j,i-1)
        enddo
    enddo
    ! Remove unconnected terms and 2-body terms and pass to B(1:276)
    do j=1, 6
        dbdr(j,1) = db1dr(j,4)
        do i=2, 4
            dbdr(j,i) = db1dr(j,i+4)
        enddo
        dbdr(j,5) = db1dr(j,10)
        do i=6, 11
            dbdr(j,i) = db1dr(j,i+6)
        enddo
        dbdr(j,12) = db1dr(j,19)
        dbdr(j,13) = db1dr(j,21)
        do i=14, 26
            dbdr(j,i) = db1dr(j,i+9)
        enddo
        dbdr(j,27) = db1dr(j,37)
        dbdr(j,28) = db1dr(j,39)
        do i=29, 53
            dbdr(j,i) = db1dr(j,i+12)
        enddo
        dbdr(j,54) = db1dr(j,67)
        dbdr(j,55) = db1dr(j,69)
        dbdr(j,56) = db1dr(j,71)
        do i=57, 97
            dbdr(j,i) = db1dr(j,i+16)
        enddo
        dbdr(j,98)  = db1dr(j,115)
        dbdr(j,99)  = db1dr(j,117)
        dbdr(j,100) = db1dr(j,119)
        do i=101, 166
            dbdr(j,i) = db1dr(j,i+20)
        enddo
        dbdr(j,167) = db1dr(j,188)
        dbdr(j,168) = db1dr(j,190)
        dbdr(j,169) = db1dr(j,192)
        dbdr(j,170) = db1dr(j,194)
        do i=171, 272
            dbdr(j,i) = db1dr(j,i+25)
        enddo
        dbdr(j,273) = db1dr(j,299)
        dbdr(j,274) = db1dr(j,301)
        dbdr(j,275) = db1dr(j,303)
        dbdr(j,276) = db1dr(j,305)
    enddo
end subroutine evdbdr

!> @brief Evaluate d3disp.
subroutine d3disp( dist, disp, dispdr, igrad, imol )
    use chempotpy_dispersion, only: evaluate_legacy_d3_diatomic
    use chempotpy_kinds,      only: wp
    implicit none
    real(wp) :: dist(1)
    real(wp) :: disp, dispdr(1)
    integer  :: igrad
    integer  :: imol
    call evaluate_legacy_d3_diatomic(dist, disp, dispdr, igrad, imol)
end subroutine d3disp

! compute energy
!> @brief Evaluate edisp.
