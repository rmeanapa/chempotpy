!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!   System:                        O3
!   Functional form:               permutation-invariant polynomials
!   Common name:                   O3 1 1Ap (1st adiab. singlet A' state)
!   Number of derivatives:         1
!   Number of bodies:              3
!   Number of electronic surfaces: 1
!   Interface:                     Section-2
!
!   References:: Z. Varga, Y. Paukku, and D. G. Truhlar,
!     "Potential energy surfaces for O + O2 collisions"
!      J. Chem. Phys. 147, 154312/1-17 (2017).
!
!   Notes:    -Mixed-exponential-Gaussian (MEG)
!              variables are applied
!             -Diatomic potential contains
!              dispersion correction
!

!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_o3_1_1ap_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, p, g, d )
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: coord_r(3)
        double precision              :: epot_b(56)
        double precision              :: epot_c(56)
        double precision              :: epot_p(0:66)
        double precision              :: epot_rm(0:7)
        double precision              :: epot_rms(3)
        double precision              :: gradt_dbdr(3,56)
        double precision              :: gradt_dmdr(3,0:7)
        double precision              :: gradt_dmsdr(3,3)
        double precision              :: gradt_dpdr(3,0:66)
        double precision              :: gradt_drdx(3,9)
        double precision              :: gradt_dvdr(3)
        double precision              :: msprmt_a
        double precision              :: msprmt_ab
        double precision              :: msprmt_ra
        double precision              :: msprmt_rb
        double precision              :: v, tx(9), tg(9)
        integer                       :: iatom, idir, j, istate
        !initialize
        coord_r      = 0.d0; epot_b = 0.d0; epot_c = 0.d0; epot_p = 0.d0; epot_rm = 0.d0; epot_rms = 0.d0; gradt_dbdr = 0.d0
        gradt_dmdr   = 0.d0; gradt_dmsdr = 0.d0; gradt_dpdr = 0.d0; gradt_drdx = 0.d0; gradt_dvdr = 0.d0; msprmt_a = 0.d0
        msprmt_ab    = 0.d0; msprmt_ra = 0.d0; msprmt_rb = 0.d0; msprmt_a = 1.24d0; msprmt_ab = 1.79d0; msprmt_ra = 1.18d0
        msprmt_rb    = 1.11d0
        epot_c(   1) = -0.101038755197d+03
        epot_c(   2) = -0.173846060523d+04
        epot_c(   3) = 0.182561485864d+04
        epot_c(   4) = 0.128228081972d+05
        epot_c(   5) = 0.131286683935d+04
        epot_c(   6) = -0.180148145918d+05
        epot_c(   7) = 0.917232862414d+04
        epot_c(   8) = -0.898548422285d+05
        epot_c(   9) = 0.236443117836d+04
        epot_c(  10) = 0.743556379026d+05
        epot_c(  11) = 0.169485321184d+06
        epot_c(  12) = -0.658211759764d+05
        epot_c(  13) = 0.293154900522d+06
        epot_c(  14) = 0.310179195421d+05
        epot_c(  15) = -0.174367178858d+05
        epot_c(  16) = -0.171278144114d+06
        epot_c(  17) = -0.196708992219d+06
        epot_c(  18) = 0.135370348197d+06
        epot_c(  19) = 0.126786481833d+06
        epot_c(  20) = -0.483041955100d+06
        epot_c(  21) = -0.421201296092d+05
        epot_c(  22) = 0.133523045430d+05
        epot_c(  23) = 0.237659552123d+06
        epot_c(  24) = -0.351906377925d+05
        epot_c(  25) = 0.376297007075d+06
        epot_c(  26) = -0.149922827153d+06
        epot_c(  27) = -0.141832776431d+06
        epot_c(  28) = 0.428298060452d+06
        epot_c(  29) = 0.338495737904d+05
        epot_c(  30) = 0.486620284500d+05
        epot_c(  31) = 0.247392770653d+04
        epot_c(  32) = -0.195510275257d+06
        epot_c(  33) = -0.507721539659d+05
        epot_c(  34) = 0.126435459710d+05
        epot_c(  35) = -0.227222097390d+05
        epot_c(  36) = -0.255415567632d+06
        epot_c(  37) = 0.145033029297d+06
        epot_c(  38) = 0.640437783008d+05
        epot_c(  39) = -0.193559554130d+06
        epot_c(  40) = -0.142880796980d+05
        epot_c(  41) = -0.368068499717d+05
        epot_c(  42) = -0.446680420671d+04
        epot_c(  43) = 0.873374592726d+05
        epot_c(  44) = 0.239902768734d+05
        epot_c(  45) = -0.108935637797d+06
        epot_c(  46) = 0.564785581551d+05
        epot_c(  47) = 0.124730127388d+05
        epot_c(  48) = 0.288994701384d+05
        epot_c(  49) = -0.542994079358d+05
        epot_c(  50) = -0.387502623052d+04
        epot_c(  51) = 0.341436319238d+05
        epot_c(  52) = -0.358667851341d+04
        epot_c(  53) = 0.690127419637d+04
        epot_c(  54) = 0.966576122637d+04
        epot_c(  55) = 0.842319423210d+03
        epot_c(  56) = -0.162693146596d+05
        v            = 0.d0; g = 0.d0; d = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                j     = 3*(iatom-1)+idir
                tx(j) = x(iatom,idir)
            enddo
        enddo
        tg = 0.d0
        call o3pes(tx, v, tg, igrad)
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
            implicit none
            integer          :: i, igrad, j, nob, k
            double precision :: v
            double precision :: dvdx(9), x(9)
            !      common /coord/    R(3)
            !      common /epot/     rMs(3),rM(0:7),P(0:66),C(56),B(56)
            !      common /gradt/    dMsdR(3,3),dMdR(3,0:7),dPdR(3,0:66),dVdR(3),
            !     $                  dRdX(3,9),dBdR(3,56)
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
        end subroutine o3pes
        !> @brief Evaluate coord_convt for the O3_1_1Ap potential.
        subroutine coord_convt( x )
            !  Program to calculate the six interatomic distance
            !  by reading XYZ coordinate
            implicit none
            integer          :: i
            double precision :: x(9)
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
            !      Write(*,9999) V
            ! 9999 Format('The potential energy is ',F20.14,' kcal/mol')
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
            implicit none
            double precision, intent(in)  :: x(9)
            double precision, intent(out) :: dvdx(9)
            call evdvdr
            call evdrdx(x)
            call transform_internal_gradient_to_cartesian(gradt_dvdr, gradt_drdx, dvdx)
        end subroutine evdvdx
        !> @brief Evaluate EvMorse for the O3_1_1Ap potential.
        subroutine evmorse
            implicit none
            integer :: i
            ! mixed exponential-Gaussian terms = exp(-(r-ra)/a-(r-rb)^2/ab)
            ! ra:   reference bond length
            ! rb:   reference bond length
            ! a:    nonlinear parameter, unit Anstrom
            ! ab:   nonlinear parameter, unit Anstrom^2
            do i=1, 3
                epot_rms(i) = exp(-(coord_r(i)-msprmt_ra)/msprmt_a-((coord_r(i)-msprmt_rb)**2.0d0)/msprmt_ab)
            enddo
        end subroutine evmorse
        !> @brief Evaluate EvMono for the O3_1_1Ap potential.
        subroutine evmono
            !  The subroutine reads six MEG variables(X) and calculates the
            !  monomials(M) that do not have usable decomposition.
            !  For A4 with max. degree 10, the number of monomials is nom.
            implicit none
            epot_rm(0) = 1.0d0
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
            !  The subroutine reads monomials(m) and calculates the
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
            implicit none
            integer          :: i
            double precision :: b1(67)
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
            implicit none
            integer          :: i, j
            double precision :: dist, v2, dv2dr, disp, dispdr(3)
            ! Initialize dVdR(3)
            do i=1, 3
                gradt_dvdr(i) = 0.0d0
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
            !
            ! dMdR:         dMdR(3,nom), nom monomials w.r.t.3 bond length
            ! dPdR:         dPdR(3,nob), nop polynomial basis functions
            !               w.r.t. 3 bond lengths
            ! M(nom):       Array to store monomials
            ! P(nob):       Array to store polynomials
            implicit none
            integer, parameter           :: first_atoms(3)  = [1, 1, 2]
            integer, parameter           :: second_atoms(3) = [2, 3, 3]
            double precision, intent(in) :: x(9)
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
            implicit none
            integer :: i, j
            ! Initialize dmsdr
            do i=1, 3
                do j=1, 3
                    gradt_dmsdr(i,j) = 0.0d0
                enddo
            enddo
            !
            ! MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
            ! dmsdr(i,j)=0  i!=j
            !
            do i=1, 3
                gradt_dmsdr(i,i) = (-2.0d0*(coord_r(i)-msprmt_rb)/msprmt_ab-1/msprmt_a)* exp(-(coord_r(i)-msprmt_ra)/ &
                    & msprmt_a-((coord_r(i)-msprmt_rb)**2.0d0)/msprmt_ab)
            enddo
        end subroutine evdmsdr
        !> @brief Evaluate EvdMdR for the O3_1_1Ap potential.
        subroutine evdmdr
            !  The subroutine reads M(nom) and dMSdR(3,3) and calculates the
            !  dMdR(3,nom) that do not have usable decomposition.
            !  For A4 with max. degree 10, the number of monomials is nom.
            implicit none
            integer :: i
            do i=1, 3
                gradt_dmdr(i,0) = 0.0d0
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
            !  The subroutine reads monomials(m) and calculates the
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
            implicit none
            integer          :: j
            integer          :: i
            double precision :: db1dr(3,67)
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
            implicit none
            double precision :: cn(3), s6, s8, rs6, rs8
            double precision :: dist(3), e6(3), e8(3), disp, dispdr(3), c6(3)
            double precision :: e6dr(3), e8dr(3)
            integer          :: iz(3), mxc(94), i, j, igrad
            double precision :: c6ab(94,94,5,5,3)
            double precision :: r2r4(94)
            double precision :: autoang, autokcal
            autoang  = 0.52917726d0
            autokcal = 627.509541d0
            ! Generalized parameters for BJ damping from P. Verma, B. Wang,
            ! L. E. Fernandez, and D. G. Truhlar, J. Phys. Chem. A 121, 2855 (2017).
            s6  = 1.0d0
            s8  = 2.0d0
            rs6 = 0.5299d0
            rs8 = 2.20d0
            do i=1, 3
                dist(i) = dist(i)/autoang
            enddo
            ! iz for O4 system
            iz(1) = 8
            iz(2) = 8
            iz(3) = 8
            ! C6 for O4 system
            c6(1) = 12.8d0
            c6(2) = 12.8d0
            c6(3) = 12.8d0
            ! Calculate dispersion correction
            call edisp(94, 5, 3, dist, iz, mxc, &
                & rs6, rs8, e6, e8, e6dr, e8dr, c6, 0)
            disp = 0.0d0
            do i=1, 3
                disp = disp + (-s6*e6(i)-s8*e8(i))*autokcal
            enddo
            if( igrad .eq. 1 )then
                call edisp(94, 5, 3, dist, iz, mxc, &
                    & rs6, rs8, e6, e8, e6dr, e8dr, c6, 1)
                dispdr(:) = 0.0d0
                do i=1, 3
                    dispdr(i) = dispdr(i) + (-s6*e6dr(i)-s8*e8dr(i))*autokcal/autoang
                enddo
            endif
            do i=1, 3
                dist(i) = dist(i)*autoang
            enddo
        end subroutine d3disp
        !> @brief Evaluate edisp for the O3_1_1Ap potential.
        subroutine edisp( max_elem, maxc, n, dist, iz, mxc, &
                & rs6, rs8, e6, e8, e6dr, e8dr, c6a, igrad )
            implicit none
            integer          :: n, iz(3), max_elem, maxc, mxc(max_elem)
            double precision :: dist(3), r2r4(max_elem), r0ab(max_elem,max_elem)
            double precision :: rs6, rs8, rcov(max_elem)
            double precision :: c6ab(max_elem,max_elem,maxc,maxc,3)
            double precision :: e6(3), e8(3), c6a(3), e6dr(3), e8dr(3)
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
end module chempotpy_o3_1_1ap_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_o3_1_1ap_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
