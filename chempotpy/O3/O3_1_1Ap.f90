!   System:                        O3
!   Functional form:               permutation-invariant polynomials
!   Common name:                   O3 1 1Ap (1st adiab. singlet A' state)
!   Number of derivatives:         1
!   Number of bodies:              3
!   Number of electronic surfaces: 1
!   Interface:                     Section-2
!   References:: Z. Varga, Y. Paukku, and D. G. Truhlar,
!     "Potential energy surfaces for O + O2 collisions"
!      J. Chem. Phys. 147, 154312/1-17 (2017).
!   Notes:    -Mixed-exponential-Gaussian (MEG)
!              variables are applied
!             -Diatomic potential contains
!              dispersion correction

!> @brief Implement O3_1_1Ap, a potential energy surface for O3.
module chempotpy_o3_1_1ap_pes
    use chempotpy_units, only: hartree_to_kilocalorie_per_mole
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the O3_1_1Ap potential energy surface for O3.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t
        use chempotpy_units,  only: electronvolt_to_kilocalorie_per_mole
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: coord_r(3)
        real(wp)              :: epot_b(56)
        real(wp)              :: epot_c(56)
        real(wp)              :: epot_p(0:66)
        real(wp)              :: epot_rm(0:7)
        real(wp)              :: epot_rms(3)
        real(wp)              :: gradt_dbdr(3,56)
        real(wp)              :: gradt_dmdr(3,0:7)
        real(wp)              :: gradt_dmsdr(3,3)
        real(wp)              :: gradt_dpdr(3,0:66)
        real(wp)              :: gradt_drdx(3,9)
        real(wp)              :: gradt_dvdr(3)
        real(wp)              :: msprmt_a
        real(wp)              :: msprmt_ab
        real(wp)              :: msprmt_ra
        real(wp)              :: msprmt_rb
        real(wp)              :: v, tx(9), tg(9)
        integer               :: iatom, idir, j, istate
        !initialize
        coord_r      = 0.e0_wp; epot_b = 0.e0_wp; epot_c = 0.e0_wp; epot_p = 0.e0_wp; epot_rm = 0.e0_wp; epot_rms = 0.e0_wp; gradt_dbdr = 0.e0_wp
        gradt_dmdr   = 0.e0_wp; gradt_dmsdr = 0.e0_wp; gradt_dpdr = 0.e0_wp; gradt_drdx = 0.e0_wp; gradt_dvdr = 0.e0_wp; msprmt_a = 0.e0_wp
        msprmt_ab    = 0.e0_wp; msprmt_ra = 0.e0_wp; msprmt_rb = 0.e0_wp; msprmt_a = 1.24e0_wp; msprmt_ab = 1.79e0_wp; msprmt_ra = 1.18e0_wp
        msprmt_rb    = 1.11e0_wp
        epot_c(   1) = -0.101038755197e+03_wp
        epot_c(   2) = -0.173846060523e+04_wp
        epot_c(   3) = 0.182561485864e+04_wp
        epot_c(   4) = 0.128228081972e+05_wp
        epot_c(   5) = 0.131286683935e+04_wp
        epot_c(   6) = -0.180148145918e+05_wp
        epot_c(   7) = 0.917232862414e+04_wp
        epot_c(   8) = -0.898548422285e+05_wp
        epot_c(   9) = 0.236443117836e+04_wp
        epot_c(  10) = 0.743556379026e+05_wp
        epot_c(  11) = 0.169485321184e+06_wp
        epot_c(  12) = -0.658211759764e+05_wp
        epot_c(  13) = 0.293154900522e+06_wp
        epot_c(  14) = 0.310179195421e+05_wp
        epot_c(  15) = -0.174367178858e+05_wp
        epot_c(  16) = -0.171278144114e+06_wp
        epot_c(  17) = -0.196708992219e+06_wp
        epot_c(  18) = 0.135370348197e+06_wp
        epot_c(  19) = 0.126786481833e+06_wp
        epot_c(  20) = -0.483041955100e+06_wp
        epot_c(  21) = -0.421201296092e+05_wp
        epot_c(  22) = 0.133523045430e+05_wp
        epot_c(  23) = 0.237659552123e+06_wp
        epot_c(  24) = -0.351906377925e+05_wp
        epot_c(  25) = 0.376297007075e+06_wp
        epot_c(  26) = -0.149922827153e+06_wp
        epot_c(  27) = -0.141832776431e+06_wp
        epot_c(  28) = 0.428298060452e+06_wp
        epot_c(  29) = 0.338495737904e+05_wp
        epot_c(  30) = 0.486620284500e+05_wp
        epot_c(  31) = 0.247392770653e+04_wp
        epot_c(  32) = -0.195510275257e+06_wp
        epot_c(  33) = -0.507721539659e+05_wp
        epot_c(  34) = 0.126435459710e+05_wp
        epot_c(  35) = -0.227222097390e+05_wp
        epot_c(  36) = -0.255415567632e+06_wp
        epot_c(  37) = 0.145033029297e+06_wp
        epot_c(  38) = 0.640437783008e+05_wp
        epot_c(  39) = -0.193559554130e+06_wp
        epot_c(  40) = -0.142880796980e+05_wp
        epot_c(  41) = -0.368068499717e+05_wp
        epot_c(  42) = -0.446680420671e+04_wp
        epot_c(  43) = 0.873374592726e+05_wp
        epot_c(  44) = 0.239902768734e+05_wp
        epot_c(  45) = -0.108935637797e+06_wp
        epot_c(  46) = 0.564785581551e+05_wp
        epot_c(  47) = 0.124730127388e+05_wp
        epot_c(  48) = 0.288994701384e+05_wp
        epot_c(  49) = -0.542994079358e+05_wp
        epot_c(  50) = -0.387502623052e+04_wp
        epot_c(  51) = 0.341436319238e+05_wp
        epot_c(  52) = -0.358667851341e+04_wp
        epot_c(  53) = 0.690127419637e+04_wp
        epot_c(  54) = 0.966576122637e+04_wp
        epot_c(  55) = 0.842319423210e+03_wp
        epot_c(  56) = -0.162693146596e+05_wp
        v            = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                j     = 3*(iatom-1)+idir
                tx(j) = x(iatom,idir)
            enddo
        enddo
        tg = 0.e0_wp
        call o3pes(tx, v, tg, igrad)
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

        !> @brief Evaluate o3pes for the O3_1_1Ap potential.
        subroutine o3pes( x, v, dvdx, igrad )
            ! Subroutine to calculate the potential energy V and gradient dVdX
            ! for given Cartesian coordinates X(9)
            ! R:            Interatomic bond distance (3)
            ! V:            Calculated potential energy
            ! dVdX:         The derivative of V w.r.t. X, dim(9)
            ! dVdR:         The derivative of V w.r.t. R, dim(3)
            ! dPdR:         The derivative of basis functions w.r.t. R
            !               dim(3*67)
            ! dMdR:         The derivative of monomials w.r.t. R
            !               dim(3*8)
            ! dRdX:         The derivative of R w.r.t. X, dim(3*9)
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: igrad
            real(wp) :: v
            real(wp) :: dvdx(9), x(9)
            !     $                  dRdX(3,9),dBdR(3,56)
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
        end subroutine o3pes

        !> @brief Evaluate coord_convt for the O3_1_1Ap potential.
        subroutine coord_convt( x )
            !  Program to calculate the six interatomic distance
            !  from the Cartesian coordinate array
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: x(9)
            !  Now, calculate the inter-atomic distance
            !  r1 = r(O1O2)
            !  r2 = r(O1O3)
            !  r3 = r(O2O3)
            coord_r(1) = sqrt((x(4)-x(1))**2 + (x(5)-x(2))**2 + (x(6)-x(3))**2)
            coord_r(2) = sqrt((x(7)-x(1))**2 + (x(8)-x(2))**2 + (x(9)-x(3))**2)
            coord_r(3) = sqrt((x(4)-x(7))**2 + (x(5)-x(8))**2 + (x(6)-x(9))**2)
        end subroutine coord_convt

        !> @brief Evaluate EvV for the O3_1_1Ap potential.
        subroutine evv( v )
            ! Subroutine to evaluate V for given R
            ! V(R) = C*P
            ! C:            Coefficients, stored in 'dim.inc'
            ! P:            Basis functions evaluated for given R
            ! rMs:          rMs(3), six mixed exponential-Gaussian terms (MEG)
            ! a:            Nonlinear parameters in Morse terms(Angstrom)
            ! ab:           Nonlinear parameters in Gauss terms(Angstrom^2)
            ! re:           Equilibrium bond length(Angstrom)
            ! nop:          number of points
            ! nom:          number of monomials
            ! nob:          number of basis functions(polynomials)
            ! rM(0:7):      Array to store monomials
            ! P(0:66):      Array to store polynomials
            ! B(1:56):      Array to store basis functions
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i
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
            do i=1, 3
                dist = coord_r(i)
                call ev2gm2(dist, v2, dv2dr, 4, 0)
                v = v+v2
            enddo
            ! Add D3 dispersion correction
            call d3disp(coord_r, disp, dispdr, 0)
            v = v+disp
            ! Evaluate V by taken the product of C and Basis function array
            do i=1, 56
                v = v + epot_c(i)*epot_b(i)
            enddo
        end subroutine evv

        !> @brief Evaluate EvdVdX for the O3_1_1Ap potential.
        subroutine evdvdx( x, dvdx )
            use chempotpy_geometry, only: transform_internal_gradient_to_cartesian
            ! Subroutine to evaluate dRdX for given R and X
            ! R:            R(3), 3 bond lengths
            ! X:            X(9), 9 Cartesian coordinates
            ! rM(0:7):      Array to store monomials
            ! P(0:66):      Array to store polynomials
            ! dVdX:         dVdX(9), derivatives of V w.r.t. Cartesian coordinates
            ! dVdR:         dVdR(3), derivatives of V w.r.t. 3 bond lengths
            ! dRdX:         dRdX(3,9), derivatives of R(3) w.r.t. 9
            !               Cartesian coordinates
            use chempotpy_kinds, only: wp
            implicit none
            real(wp), intent(in)  :: x(9)
            real(wp), intent(out) :: dvdx(9)
            call evdvdr
            call evdrdx(x)
            call transform_internal_gradient_to_cartesian(gradt_dvdr, gradt_drdx, dvdx)
        end subroutine evdvdx

        !> @brief Evaluate EvMorse for the O3_1_1Ap potential.
        subroutine evmorse
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i
            ! mixed exponential-Gaussian terms = exp(-(r-ra)/a-(r-rb)^2/ab)
            ! ra:   reference bond length
            ! rb:   reference bond length
            ! a:    nonlinear parameter, unit Anstrom
            ! ab:   nonlinear parameter, unit Anstrom^2
            do i=1, 3
                epot_rms(i) = exp(-(coord_r(i)-msprmt_ra)/msprmt_a-((coord_r(i)-msprmt_rb)**2.0e0_wp)/msprmt_ab)
            enddo
        end subroutine evmorse

        !> @brief Evaluate EvMono for the O3_1_1Ap potential.
        subroutine evmono
            !  The subroutine accepts six MEG variables (X) and calculates the
            !  monomials(M) that do not have usable decomposition.
            !  For A4 with max. degree 10, the number of monomials is nom.
            use chempotpy_kinds, only: wp
            implicit none
            epot_rm(0) = 1.0e0_wp
            epot_rm(1) = epot_rms(3)
            epot_rm(2) = epot_rms(2)
            epot_rm(3) = epot_rms(1)
            epot_rm(4) = epot_rm(1)*epot_rm(2)
            epot_rm(5) = epot_rm(1)*epot_rm(3)
            epot_rm(6) = epot_rm(2)*epot_rm(3)
            epot_rm(7) = epot_rm(1)*epot_rm(6)
        end subroutine evmono

        !> @brief Evaluate EvPoly for the O3_1_1Ap potential.
        subroutine evpoly
            !  The subroutine accepts monomials (m) and calculates the
            !  permutation invariant polynomials(p)
            !  For A4 with max. degree 10, the number of polynomials is nob.
            implicit none
            epot_p( 0) = epot_rm(0)
            epot_p( 1) = epot_rm(1) + epot_rm(2) + epot_rm(3)
            epot_p( 2) = epot_rm(4) + epot_rm(5) + epot_rm(6)
            epot_p( 3) = epot_p( 1)*epot_p( 1) - epot_p( 2) - epot_p( 2)
            epot_p( 4) = epot_rm(7)
            epot_p( 5) = epot_p( 1)*epot_p( 2) - epot_p( 4) - epot_p( 4) - epot_p( 4)
            epot_p( 6) = epot_p( 1)*epot_p( 3) - epot_p( 5)
            epot_p( 7) = epot_p( 1)*epot_p( 4)
            epot_p( 8) = epot_p( 2)*epot_p( 2) - epot_p( 7) - epot_p( 7)
            epot_p( 9) = epot_p( 2)*epot_p( 3) - epot_p( 7)
            epot_p(10) = epot_p( 1)*epot_p( 6) - epot_p( 9)
            epot_p(11) = epot_p( 2)*epot_p( 4)
            epot_p(12) = epot_p( 3)*epot_p( 4)
            epot_p(13) = epot_p( 1)*epot_p( 8) - epot_p(11)
            epot_p(14) = epot_p( 2)*epot_p( 6) - epot_p(12)
            epot_p(15) = epot_p( 1)*epot_p(10) - epot_p(14)
            epot_p(16) = epot_p( 4)*epot_p( 4)
            epot_p(17) = epot_p( 4)*epot_p( 5)
            epot_p(18) = epot_p( 4)*epot_p( 6)
            epot_p(19) = epot_p( 2)*epot_p( 8) - epot_p(17)
            epot_p(20) = epot_p( 1)*epot_p(13) - epot_p(17) - epot_p(19) - epot_p(19)
            epot_p(21) = epot_p( 2)*epot_p(10) - epot_p(18)
            epot_p(22) = epot_p( 1)*epot_p(15) - epot_p(21)
            epot_p(23) = epot_p( 1)*epot_p(16)
            epot_p(24) = epot_p( 4)*epot_p( 8)
            epot_p(25) = epot_p( 3)*epot_p(11) - epot_p(23)
            epot_p(26) = epot_p( 4)*epot_p(10)
            epot_p(27) = epot_p( 1)*epot_p(19) - epot_p(24)
            epot_p(28) = epot_p( 6)*epot_p( 8) - epot_p(23)
            epot_p(29) = epot_p( 2)*epot_p(15) - epot_p(26)
            epot_p(30) = epot_p( 1)*epot_p(22) - epot_p(29)
            epot_p(31) = epot_p( 2)*epot_p(16)
            epot_p(32) = epot_p( 3)*epot_p(16)
            epot_p(33) = epot_p( 1)*epot_p(24) - epot_p(31)
            epot_p(34) = epot_p( 4)*epot_p(14)
            epot_p(35) = epot_p( 4)*epot_p(15)
            epot_p(36) = epot_p( 2)*epot_p(19) - epot_p(33)
            epot_p(37) = epot_p( 3)*epot_p(19) - epot_p(31)
            epot_p(38) = epot_p( 8)*epot_p(10) - epot_p(32)
            epot_p(39) = epot_p( 2)*epot_p(22) - epot_p(35)
            epot_p(40) = epot_p( 1)*epot_p(30) - epot_p(39)
            epot_p(41) = epot_p( 4)*epot_p(16)
            epot_p(42) = epot_p( 4)*epot_p(17)
            epot_p(43) = epot_p( 4)*epot_p(19)
            epot_p(44) = epot_p( 6)*epot_p(16)
            epot_p(45) = epot_p( 4)*epot_p(20)
            epot_p(46) = epot_p( 4)*epot_p(21)
            epot_p(47) = epot_p( 4)*epot_p(22)
            epot_p(48) = epot_p( 1)*epot_p(36) - epot_p(43)
            epot_p(49) = epot_p( 1)*epot_p(37) - epot_p(45) - epot_p(48)
            epot_p(50) = epot_p( 8)*epot_p(15) - epot_p(44)
            epot_p(51) = epot_p( 2)*epot_p(30) - epot_p(47)
            epot_p(52) = epot_p( 1)*epot_p(40) - epot_p(51)
            epot_p(53) = epot_p( 1)*epot_p(41)
            epot_p(54) = epot_p( 4)*epot_p(24)
            epot_p(55) = epot_p( 3)*epot_p(31) - epot_p(53)
            epot_p(56) = epot_p( 1)*epot_p(43) - epot_p(54)
            epot_p(57) = epot_p(10)*epot_p(16)
            epot_p(58) = epot_p( 4)*epot_p(28)
            epot_p(59) = epot_p( 4)*epot_p(29)
            epot_p(60) = epot_p( 4)*epot_p(30)
            epot_p(61) = epot_p( 2)*epot_p(36) - epot_p(56)
            epot_p(62) = epot_p( 3)*epot_p(36) - epot_p(54)
            epot_p(63) = epot_p(10)*epot_p(19) - epot_p(53)
            epot_p(64) = epot_p( 8)*epot_p(22) - epot_p(57)
            epot_p(65) = epot_p( 2)*epot_p(40) - epot_p(60)
            epot_p(66) = epot_p( 1)*epot_p(52) - epot_p(65)
        end subroutine evpoly

        !> @brief Evaluate evbas for the O3_1_1Ap potential.
        subroutine evbas
            !  The subroutine eliminate the 2-body terms in Bowman's approach
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i
            real(wp) :: b1(67)
            ! Pass P(0:66) to BM1(1:67)
            do i=1, 67
                b1(i) = epot_p(i-1)
            enddo
            ! Remove unconnected terms and 2-body terms and pass to B(1:430)
            epot_b(1) = b1(3)
            do i=2, 3
                epot_b(i) = b1(i+3)
            enddo
            do i=4, 6
                epot_b(i) = b1(i+4)
            enddo
            do i=7, 10
                epot_b(i) = b1(i+5)
            enddo
            do i=11, 16
                epot_b(i) = b1(i+6)
            enddo
            do i=17, 23
                epot_b(i) = b1(i+7)
            enddo
            do i=24, 32
                epot_b(i) = b1(i+8)
            enddo
            do i=33, 43
                epot_b(i) = b1(i+9)
            enddo
            do i=44, 56
                epot_b(i) = b1(i+10)
            enddo
        end subroutine evbas

        !> @brief Evaluate ev2gm2 for the O3_1_1Ap potential.
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
            ! Anstromgs^-2, beta=dimensionless, a_k in milihartree.
            real(wp) :: alpha, beta, a(0:7)
            integer  :: k, igrad
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
            ! From milihartree to kcal/mol
            v = v*hartree_to_kilocalorie_per_mole*1.0e-3_wp
            ! Compute the gradient if needed
            if( igrad .eq. 1 )then
                grad = 0.e0_wp
                do k=0, 7
                    grad = grad-2.e0_wp*a(k)*alpha*beta**k*r*exp(-alpha*beta**k*r**2)
                enddo
                ! Convert from milihartree/A to i(kcal/mol)/A
                grad = grad*hartree_to_kilocalorie_per_mole*1.0e-3_wp
            endif
        end subroutine ev2gm2

        !> @brief Evaluate EvdVdR for the O3_1_1Ap potential.
        subroutine evdvdr
            ! Subroutine to evaluate dVdR for given R
            ! dVdR = dV2dR + C*dBdR
            ! C:            Coefficients, stored in 'dim.inc'
            ! P:            Basis functions evaluated for given R
            ! M:            Monomials evaluated for given R
            ! dV2dR:        Gradient of 2-body interactions
            ! dMsdR:        dMsdR(3,3), 3 MEG terms w.r.t. 3 bond lengths
            ! dMdR:         dMdR(3,nom), nom monomials w.r.t.3 bond length
            ! dPdR:         dPdR(3,nob), nop polynomial basis functions
            !               w.r.t. 3 bond lengths
            ! nom:          number of monomials
            ! nob:          number of basis functions(polynomials)
            ! M(nom):       Array to store monomials
            ! P(nob):       Array to store polynomials
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: i, j
            real(wp) :: dist, v2, dv2dr, disp, dispdr(3)
            ! Initialize dVdR(3)
            do i=1, 3
                gradt_dvdr(i) = 0.0e0_wp
            enddo
            ! Add dV2dR(i) to dVdR
            do i=1, 3
                dist = coord_r(i)
                call ev2gm2(dist, v2, dv2dr, 4, 1)
                gradt_dvdr(i) = dv2dr
            enddo
            ! Add numerical gradient of D3 dispersion correction
            call d3disp(coord_r, disp, dispdr, 1)
            do i=1, 3
                gradt_dvdr(i) = gradt_dvdr(i) + dispdr(i)
            enddo
            ! Calculate dMEG/dr(3,3) for given R(3)
            call evdmsdr
            ! Calculate the monomials for each point by using six MEG terms
            call evdmdr
            ! Calculate the polynomials by using monomials
            call evdpdr
            ! Remove 2-body interactions and unconnected terms from polynomials
            call evdbdr
            ! Evaluate dVdR(3) by taken the product of C(j) and dPdR(i,j)
            do i=1, 3
                do j=1, 56
                    gradt_dvdr(i) = gradt_dvdr(i) + epot_c(j)*gradt_dbdr(i,j)
                enddo
            enddo
        end subroutine evdvdr

        !> @brief Evaluate EvdRdX for the O3_1_1Ap potential.
        subroutine evdrdx( x )
            use chempotpy_geometry, only: ordered_pair_distance_jacobian
            ! Subroutine to evaluate dRdX for given R and X
            ! R:            R(3), 3 bond lengths
            ! X:            X(9), 9 Cartesian coordinates
            ! dMdR:         dMdR(3,nom), nom monomials w.r.t.3 bond length
            ! dPdR:         dPdR(3,nob), nop polynomial basis functions
            !               w.r.t. 3 bond lengths
            ! M(nom):       Array to store monomials
            ! P(nob):       Array to store polynomials
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter   :: first_atoms(3)  = [1, 1, 2]
            integer, parameter   :: second_atoms(3) = [2, 3, 3]
            real(wp), intent(in) :: x(9)
            call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, coord_r, gradt_drdx)
        end subroutine evdrdx

        !> @brief Evaluate EvdMsdR for the O3_1_1Ap potential.
        subroutine evdmsdr
            ! Subroutine to evaluate the derivatives of MEG term X
            ! w.r.t. interatomic distance R(3)
            ! dmsdR:        Local variables, dirm(3,3)
            ! a:            Nonlinear pamameter(Angstrom)
            ! ab:           Nonlinear pamameter(Angstrom^2)
            ! ra:           reference bond length(Angstrom)
            ! rb:           reference bond length(Angstrom)
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, j
            ! Initialize dmsdr
            do i=1, 3
                do j=1, 3
                    gradt_dmsdr(i,j) = 0.0e0_wp
                enddo
            enddo
            ! MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
            ! dmsdr(i,j)=0  i!=j
            do i=1, 3
                gradt_dmsdr(i,i) = (-2.0e0_wp*(coord_r(i)-msprmt_rb)/msprmt_ab-1/msprmt_a)* exp(-(coord_r(i)- &
                    & msprmt_ra)/ msprmt_a-((coord_r(i)-msprmt_rb)**2.0e0_wp)/msprmt_ab)
            enddo
        end subroutine evdmsdr

        !> @brief Evaluate EvdMdR for the O3_1_1Ap potential.
        subroutine evdmdr
            !  The subroutine accepts M(nom) and dMSdR(3,3) and calculates the
            !  dMdR(3,nom) that do not have usable decomposition.
            !  For A4 with max. degree 10, the number of monomials is nom.
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i
            do i=1, 3
                gradt_dmdr(i,0) = 0.0e0_wp
                gradt_dmdr(i,1) = gradt_dmsdr(i,3)
                gradt_dmdr(i,2) = gradt_dmsdr(i,2)
                gradt_dmdr(i,3) = gradt_dmsdr(i,1)
                gradt_dmdr(i,4) = gradt_dmdr(i,1)*epot_rm(2) + epot_rm(1)*gradt_dmdr(i,2)
                gradt_dmdr(i,5) = gradt_dmdr(i,1)*epot_rm(3) + epot_rm(1)*gradt_dmdr(i,3)
                gradt_dmdr(i,6) = gradt_dmdr(i,2)*epot_rm(3) + epot_rm(2)*gradt_dmdr(i,3)
                gradt_dmdr(i,7) = gradt_dmdr(i,1)*epot_rm(6) + epot_rm(1)*gradt_dmdr(i,6)
            enddo
        end subroutine evdmdr

        !> @brief Evaluate EvdPdr for the O3_1_1Ap potential.
        subroutine evdpdr
            !  The subroutine accepts monomials (m) and calculates the
            !  permutation invariant polynomials(p)
            !  For A4 with max. degree 10, the number of polynomials is nob.
            implicit none
            integer :: i
            do i=1, 3
                gradt_dpdr(i,0) = gradt_dmdr(i,0)
                gradt_dpdr(i,1) = gradt_dmdr(i,1) + gradt_dmdr(i,2) + gradt_dmdr(i,3)
                gradt_dpdr(i,2) = gradt_dmdr(i,4) + gradt_dmdr(i,5) + gradt_dmdr(i,6)
                gradt_dpdr(i,3)  = gradt_dpdr(i,1)*epot_p( 1) + epot_p( 1)*gradt_dpdr(i,1) - gradt_dpdr(i,2) - &
                    & gradt_dpdr(i,2)
                gradt_dpdr(i,4) = gradt_dmdr(i,7)
                gradt_dpdr(i,5)  = gradt_dpdr(i,1)*epot_p( 2) + epot_p( 1)*gradt_dpdr(i,2) - gradt_dpdr(i,4) - &
                    & gradt_dpdr(i,4) - gradt_dpdr(i,4)
                gradt_dpdr(i,6) = gradt_dpdr(i,1)*epot_p( 3) + epot_p( 1)*gradt_dpdr(i,3) - gradt_dpdr(i,5)
                gradt_dpdr(i,7) = gradt_dpdr(i,1)*epot_p( 4) + epot_p( 1)*gradt_dpdr(i,4)
                gradt_dpdr(i,8)  = gradt_dpdr(i,2)*epot_p( 2) + epot_p( 2)*gradt_dpdr(i,2) - gradt_dpdr(i,7) - &
                    & gradt_dpdr(i,7)
                gradt_dpdr(i,9)  = gradt_dpdr(i,2)*epot_p( 3) + epot_p( 2)*gradt_dpdr(i,3) - gradt_dpdr(i,7)
                gradt_dpdr(i,10) = gradt_dpdr(i,1)*epot_p( 6) + epot_p( 1)*gradt_dpdr(i,6) - gradt_dpdr(i,9)
                gradt_dpdr(i,11) = gradt_dpdr(i,2)*epot_p( 4) + epot_p( 2)*gradt_dpdr(i,4)
                gradt_dpdr(i,12) = gradt_dpdr(i,3)*epot_p( 4) + epot_p( 3)*gradt_dpdr(i,4)
                gradt_dpdr(i,13) = gradt_dpdr(i,1)*epot_p( 8) + epot_p( 1)*gradt_dpdr(i,8) - gradt_dpdr(i,11)
                gradt_dpdr(i,14) = gradt_dpdr(i,2)*epot_p( 6) + epot_p( 2)*gradt_dpdr(i,6) - gradt_dpdr(i,12)
                gradt_dpdr(i,15) = gradt_dpdr(i,1)*epot_p(10) + epot_p( 1)*gradt_dpdr(i,10) - gradt_dpdr(i,14)
                gradt_dpdr(i,16) = gradt_dpdr(i,4)*epot_p( 4) + epot_p( 4)*gradt_dpdr(i,4)
                gradt_dpdr(i,17) = gradt_dpdr(i,4)*epot_p( 5) + epot_p( 4)*gradt_dpdr(i,5)
                gradt_dpdr(i,18) = gradt_dpdr(i,4)*epot_p( 6) + epot_p( 4)*gradt_dpdr(i,6)
                gradt_dpdr(i,19) = gradt_dpdr(i,2)*epot_p( 8) + epot_p( 2)*gradt_dpdr(i,8) - gradt_dpdr(i,17)
                gradt_dpdr(i,20) = gradt_dpdr(i,1)*epot_p(13) + epot_p( 1)*gradt_dpdr(i,13) - gradt_dpdr(i,17) - &
                    & gradt_dpdr(i,19) - gradt_dpdr(i,19)
                gradt_dpdr(i,21) = gradt_dpdr(i,2)*epot_p(10) + epot_p( 2)*gradt_dpdr(i,10) - gradt_dpdr(i,18)
                gradt_dpdr(i,22) = gradt_dpdr(i,1)*epot_p(15) + epot_p( 1)*gradt_dpdr(i,15) - gradt_dpdr(i,21)
                gradt_dpdr(i,23) = gradt_dpdr(i,1)*epot_p(16) + epot_p( 1)*gradt_dpdr(i,16)
                gradt_dpdr(i,24) = gradt_dpdr(i,4)*epot_p( 8) + epot_p( 4)*gradt_dpdr(i,8)
                gradt_dpdr(i,25) = gradt_dpdr(i,3)*epot_p(11) + epot_p( 3)*gradt_dpdr(i,11) - gradt_dpdr(i,23)
                gradt_dpdr(i,26) = gradt_dpdr(i,4)*epot_p(10) + epot_p( 4)*gradt_dpdr(i,10)
                gradt_dpdr(i,27) = gradt_dpdr(i,1)*epot_p(19) + epot_p( 1)*gradt_dpdr(i,19) - gradt_dpdr(i,24)
                gradt_dpdr(i,28) = gradt_dpdr(i,6)*epot_p( 8) + epot_p( 6)*gradt_dpdr(i,8) - gradt_dpdr(i,23)
                gradt_dpdr(i,29) = gradt_dpdr(i,2)*epot_p(15) + epot_p( 2)*gradt_dpdr(i,15) - gradt_dpdr(i,26)
                gradt_dpdr(i,30) = gradt_dpdr(i,1)*epot_p(22) + epot_p( 1)*gradt_dpdr(i,22) - gradt_dpdr(i,29)
                gradt_dpdr(i,31) = gradt_dpdr(i,2)*epot_p(16) + epot_p( 2)*gradt_dpdr(i,16)
                gradt_dpdr(i,32) = gradt_dpdr(i,3)*epot_p(16) + epot_p( 3)*gradt_dpdr(i,16)
                gradt_dpdr(i,33) = gradt_dpdr(i,1)*epot_p(24) + epot_p( 1)*gradt_dpdr(i,24) - gradt_dpdr(i,31)
                gradt_dpdr(i,34) = gradt_dpdr(i,4)*epot_p(14) + epot_p( 4)*gradt_dpdr(i,14)
                gradt_dpdr(i,35) = gradt_dpdr(i,4)*epot_p(15) + epot_p( 4)*gradt_dpdr(i,15)
                gradt_dpdr(i,36) = gradt_dpdr(i,2)*epot_p(19) + epot_p( 2)*gradt_dpdr(i,19) - gradt_dpdr(i,33)
                gradt_dpdr(i,37) = gradt_dpdr(i,3)*epot_p(19) + epot_p( 3)*gradt_dpdr(i,19) - gradt_dpdr(i,31)
                gradt_dpdr(i,38) = gradt_dpdr(i,8)*epot_p(10) + epot_p( 8)*gradt_dpdr(i,10) - gradt_dpdr(i,32)
                gradt_dpdr(i,39) = gradt_dpdr(i,2)*epot_p(22) + epot_p( 2)*gradt_dpdr(i,22) - gradt_dpdr(i,35)
                gradt_dpdr(i,40) = gradt_dpdr(i,1)*epot_p(30) + epot_p( 1)*gradt_dpdr(i,30) - gradt_dpdr(i,39)
                gradt_dpdr(i,41) = gradt_dpdr(i,4)*epot_p(16) + epot_p( 4)*gradt_dpdr(i,16)
                gradt_dpdr(i,42) = gradt_dpdr(i,4)*epot_p(17) + epot_p( 4)*gradt_dpdr(i,17)
                gradt_dpdr(i,43) = gradt_dpdr(i,4)*epot_p(19) + epot_p( 4)*gradt_dpdr(i,19)
                gradt_dpdr(i,44) = gradt_dpdr(i,6)*epot_p(16) + epot_p( 6)*gradt_dpdr(i,16)
                gradt_dpdr(i,45) = gradt_dpdr(i,4)*epot_p(20) + epot_p( 4)*gradt_dpdr(i,20)
                gradt_dpdr(i,46) = gradt_dpdr(i,4)*epot_p(21) + epot_p( 4)*gradt_dpdr(i,21)
                gradt_dpdr(i,47) = gradt_dpdr(i,4)*epot_p(22) + epot_p( 4)*gradt_dpdr(i,22)
                gradt_dpdr(i,48) = gradt_dpdr(i,1)*epot_p(36) + epot_p( 1)*gradt_dpdr(i,36) - gradt_dpdr(i,43)
                gradt_dpdr(i,49) = gradt_dpdr(i,1)*epot_p(37) + epot_p( 1)*gradt_dpdr(i,37) - gradt_dpdr(i,45) - &
                    & gradt_dpdr(i,48)
                gradt_dpdr(i,50) = gradt_dpdr(i,8)*epot_p(15) + epot_p( 8)*gradt_dpdr(i,15) - gradt_dpdr(i,44)
                gradt_dpdr(i,51) = gradt_dpdr(i,2)*epot_p(30) + epot_p( 2)*gradt_dpdr(i,30) - gradt_dpdr(i,47)
                gradt_dpdr(i,52) = gradt_dpdr(i,1)*epot_p(40) + epot_p( 1)*gradt_dpdr(i,40) - gradt_dpdr(i,51)
                gradt_dpdr(i,53) = gradt_dpdr(i,1)*epot_p(41) + epot_p( 1)*gradt_dpdr(i,41)
                gradt_dpdr(i,54) = gradt_dpdr(i,4)*epot_p(24) + epot_p( 4)*gradt_dpdr(i,24)
                gradt_dpdr(i,55) = gradt_dpdr(i,3)*epot_p(31) + epot_p( 3)*gradt_dpdr(i,31) - gradt_dpdr(i,53)
                gradt_dpdr(i,56) = gradt_dpdr(i,1)*epot_p(43) + epot_p( 1)*gradt_dpdr(i,43) - gradt_dpdr(i,54)
                gradt_dpdr(i,57) = gradt_dpdr(i,10)*epot_p(16) + epot_p(10)*gradt_dpdr(i,16)
                gradt_dpdr(i,58) = gradt_dpdr(i,4)*epot_p(28) + epot_p( 4)*gradt_dpdr(i,28)
                gradt_dpdr(i,59) = gradt_dpdr(i,4)*epot_p(29) + epot_p( 4)*gradt_dpdr(i,29)
                gradt_dpdr(i,60) = gradt_dpdr(i,4)*epot_p(30) + epot_p( 4)*gradt_dpdr(i,30)
                gradt_dpdr(i,61) = gradt_dpdr(i,2)*epot_p(36) + epot_p( 2)*gradt_dpdr(i,36) - gradt_dpdr(i,56)
                gradt_dpdr(i,62) = gradt_dpdr(i,3)*epot_p(36) + epot_p( 3)*gradt_dpdr(i,36) - gradt_dpdr(i,54)
                gradt_dpdr(i,63) = gradt_dpdr(i,10)*epot_p(19) + epot_p(10)*gradt_dpdr(i,19) - gradt_dpdr(i,53)
                gradt_dpdr(i,64) = gradt_dpdr(i,8)*epot_p(22) + epot_p( 8)*gradt_dpdr(i,22) - gradt_dpdr(i,57)
                gradt_dpdr(i,65) = gradt_dpdr(i,2)*epot_p(40) + epot_p( 2)*gradt_dpdr(i,40) - gradt_dpdr(i,60)
                gradt_dpdr(i,66) = gradt_dpdr(i,1)*epot_p(52) + epot_p( 1)*gradt_dpdr(i,52) - gradt_dpdr(i,65)
            enddo
        end subroutine evdpdr

        !> @brief Evaluate evdbdr for the O3_1_1Ap potential.
        subroutine evdbdr
            !  The subroutine eliminate the 2-body terms in Bowman's approach
            use chempotpy_kinds, only: wp
            implicit none
            integer  :: j
            integer  :: i
            real(wp) :: db1dr(3,67)
            ! Pass P(3,0:66) to BM1(3,1:67)
            do j=1, 3
                do i=1, 67
                    db1dr(j,i) = gradt_dpdr(j,i-1)
                enddo
            enddo
            ! Remove unconnected terms and 2-body terms and pass to B(1:56)
            do j=1, 3
                gradt_dbdr(j,1) = db1dr(j,3)
                do i=2, 3
                    gradt_dbdr(j,i) = db1dr(j,i+3)
                enddo
                do i=4, 6
                    gradt_dbdr(j,i) = db1dr(j,i+4)
                enddo
                do i=7, 10
                    gradt_dbdr(j,i) = db1dr(j,i+5)
                enddo
                do i=11, 16
                    gradt_dbdr(j,i) = db1dr(j,i+6)
                enddo
                do i=17, 23
                    gradt_dbdr(j,i) = db1dr(j,i+7)
                enddo
                do i=24, 32
                    gradt_dbdr(j,i) = db1dr(j,i+8)
                enddo
                do i=33, 43
                    gradt_dbdr(j,i) = db1dr(j,i+9)
                enddo
                do i=44, 56
                    gradt_dbdr(j,i) = db1dr(j,i+10)
                enddo
            enddo
        end subroutine evdbdr

        !> @brief Evaluate d3disp for the O3_1_1Ap potential.
        subroutine d3disp( dist, disp, dispdr, igrad )
            use chempotpy_dispersion, only: evaluate_legacy_d3_o3
            use chempotpy_kinds,      only: wp
            implicit none
            real(wp) :: dist(3)
            real(wp) :: disp, dispdr(3)
            integer  :: igrad
            call evaluate_legacy_d3_o3(dist, disp, dispdr, igrad)
        end subroutine d3disp
        !> @brief Evaluate edisp for the O3_1_1Ap potential.
    end subroutine evaluate
end module chempotpy_o3_1_1ap_pes

!> @brief Expose O3_1_1Ap through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_o3_1_1ap_pes, only: evaluate
    use chempotpy_kinds,        only: wp
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
