!   System:          FH2
!   Functional form: eLEPS (extended-London-Eyring-Polanyi-Sato)
!                    plus E3C (a three-center term).
!   Common name: no. 5A
!   Reference:       R. Steckler, D. G. Truhlar, and B. C. Garrett
!                    J. Chem. Phys. 82, 5499-5505 (1985).
!   Notes:           In this version of the potential the AB and AC
!                    Sato parameters are functions of PHI, where PHI is
!                    the angle of BC-to-A vector with the BC axis.
!                    The BC Sato parameter is a function of BC bond
!                    length and the angle PHI.  A three-center term is
!                    also added to correct the H-F-H barrier height.
!   prepare_local_context must be called once before any calls to POT.
!   The potential parameters are included in DATA statements
!   and in the block data subprogram PTPACM.
!   Coordinates, potential energy, and derivatives are passed
!   The potential energy in the three asymptotic valleys are
!   stored in the common block ASYCM:
!   The potential energy in the AB valley, EASYAB, is equal to the poten
!   energy of the H "infinitely" far from the FH diatomic, with the
!   FH diatomic at its equilibrium configuration.  Similarly, the terms
!   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
!   respectively.
!   All the information passed through the common blocks PT1CM and ASYCM
!   is in Hartree atomic units.
!   This potential is written such that:
!                  R(1) = R(F-H)
!                  R(2) = R(H-H)
!                  R(3) = R(H-F)
!   The zero of energy is defined at F "infinitely" far from the H2 diat
!   The flags that indicate what calculations should be carried out in
!   the potential routine are passed through the common block PT2CM:
!   where:
!        surface_requests - which electronic states availalble
!                 (1,1) = 1 as only gs state available
!        requested_derivative_order  = 0 => no derivatives should be calculated
!        requested_derivative_order  = 1 => calculate first derivatives
!        NFLAG - these integer values can be used to flag options
!                within the potential;
!   Potential parameters' default settings
!                  Variable            Default value
!                  requested_derivative_order                1
!                  NFLAG(18)           6

!> @brief Implement FH2_LEPS_NO5a_1985, a potential energy surface for FH2.
module chempotpy_fh2_leps_no5a_1985_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the FH2_LEPS_NO5a_1985 potential energy surface for FH2.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units,  only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter    :: n3atom  = 75
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: a
        real(wp)              :: a1
        real(wp)              :: a2
        real(wp)              :: alf
        real(wp)              :: alfp
        real(wp)              :: alft
        real(wp)              :: b
        real(wp)              :: bet
        real(wp)              :: beta(3)
        real(wp)              :: betp
        real(wp)              :: c
        real(wp)              :: c1
        real(wp)              :: c2
        real(wp)              :: c3
        real(wp)              :: c3c
        real(wp)              :: c4
        real(wp)              :: c5
        real(wp)              :: c6
        real(wp)              :: input_cartesian(natom,3)
        real(wp)              :: cartesian_coordinates(natom,3)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: deg1
        real(wp)              :: deg2
        real(wp)              :: deg3
        real(wp)              :: deg3p
        real(wp)              :: deg4
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        real(wp)              :: dj(3,3)
        real(wp)              :: dm(3)
        real(wp)              :: dq(3,3)
        real(wp)              :: dz(3,3)
        real(wp)              :: easyab
        real(wp)              :: easyac
        real(wp)              :: easybc
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: ground_state_energy
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: ezero(isurf+1)
        real(wp)              :: f
        real(wp)              :: g2
        real(wp)              :: ha2
        integer               :: coordinate_representation
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        real(wp)              :: hh_sato_exponent
        real(wp)              :: noncollinear_width
        integer               :: atom_count
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: omp3c
        real(wp)              :: omq3c
        real(wp)              :: p3c
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: q(3)
        real(wp)              :: q3c
        real(wp)              :: r(n3atom)
        real(wp)              :: re(3)
        character(len=75)     :: ref(5)
        real(wp)              :: rr(3)
        real(wp)              :: xj(3)
        real(wp)              :: z(3)
        real(wp)              :: zfh0
        real(wp)              :: zhh0
        real(wp)              :: zhhdif
        real(wp)              :: zhhst
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        a                                 = 0.e0_wp; a1 = 0.e0_wp; a2 = 0.e0_wp; alf = 0.e0_wp; alfp = 0.e0_wp; alft = 0.e0_wp; b = 0.e0_wp; bet = 0.e0_wp
        beta                              = 0.e0_wp; betp = 0.e0_wp; c = 0.e0_wp; c1 = 0.e0_wp; c2 = 0.e0_wp; c3 = 0.e0_wp; c3c = 0.e0_wp; c4 = 0.e0_wp; c5 = 0.e0_wp; c6 = 0.e0_wp
        input_cartesian                   = 0.e0_wp; cartesian_coordinates = 0.e0_wp; excited_state_internal_gradients = 0.e0_wp; deg1 = 0.e0_wp
        deg2                              = 0.e0_wp; deg3 = 0.e0_wp; deg3p = 0.e0_wp; deg4 = 0.e0_wp; ground_state_internal_gradient = 0.e0_wp; diabatic_internal_gradients = 0.e0_wp; output_excited_state_gradients = 0.e0_wp
        excited_state_cartesian_gradients = 0.e0_wp; output_ground_state_gradient = 0.e0_wp; ground_state_cartesian_gradient = 0.e0_wp; output_diabatic_gradients = 0.e0_wp; diabatic_cartesian_gradients = 0.e0_wp; dj = 0.e0_wp; dm = 0.e0_wp
        dq                                = 0.e0_wp; dz = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp; excited_state_energies = 0.e0_wp; ground_state_energy = 0.e0_wp; diabatic_couplings = 0.e0_wp
        ezero                             = 0.e0_wp; f = 0.e0_wp; g2 = 0.e0_wp; ha2 = 0.e0_wp; coordinate_representation = 0; last_requested_excited_state = 0
        last_requested_coupling           = 0; first_requested_excited_state = 0; first_requested_coupling = 0; hh_sato_exponent = 0.e0_wp; noncollinear_width = 0.e0_wp; atom_count = 0; maximum_derivative_order = 0; maximum_excited_state = 0
        surface_requests                  = 0; requested_derivative_order = 0; omp3c = 0.e0_wp; omq3c = 0.e0_wp; p3c = 0.e0_wp; output_excited_state_energies = 0.e0_wp
        output_ground_state_energy        = 0.e0_wp; output_diabatic_couplings = 0.e0_wp; q = 0.e0_wp; q3c = 0.e0_wp; r = 0.e0_wp; re = 0.e0_wp
        ref                               = ' '
        rr                                = 0.e0_wp
        xj                                = 0.e0_wp
        z                                 = 0.e0_wp
        zfh0                              = 0.e0_wp
        zhh0                              = 0.e0_wp
        zhhdif                            = 0.e0_wp
        zhhst                             = 0.e0_wp
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 0
        coordinate_representation  = 3
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 3
        dm                         = [141.196e0_wp, 109.449e0_wp, 141.196e0_wp]
        re                         = [0.9170e0_wp, 0.7419e0_wp, 0.9170e0_wp]
        beta                       = [2.2187e0_wp, 1.9420e0_wp, 2.2187e0_wp]
        a1                         = 0.0395e0_wp
        a2                         = 0.201e0_wp
        a                          = 1.26e0_wp
        b                          = 1.60e0_wp
        c                          = 0.0e0_wp
        hh_sato_exponent           = 0.0e0_wp
        f                          = 0.0e0_wp
        noncollinear_width         = 0.2e0_wp
        bet                        = -0.101168e0_wp
        betp                       = -0.46183e0_wp
        zfh0                       = 0.170e0_wp
        zhh0                       = 0.120e0_wp
        zhhst                      = 0.23107e0_wp
        deg1                       = 10.0e0_wp
        deg2                       = 20.0e0_wp
        deg3                       = 60.0e0_wp
        deg4                       = 90.0e0_wp
        deg3p                      = 41.0e0_wp
        c1                         = 1.15e-4_wp
        c2                         = -6.0e-6_wp
        c3                         = 0.016e0_wp
        c4                         = 0.0005e0_wp
        c5                         = 0.00069e0_wp
        c6                         = 3.308e-2_wp
        c3c                        = 0.270063e0_wp
        alf                        = 0.18e0_wp
        alfp                       = 0.078e0_wp
        alft                       = 2.14e0_wp
        p3c                        = 0.95e0_wp
        q3c                        = 0.48e0_wp
        v                          = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        input_cartesian            = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                input_cartesian(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        requested_derivative_order = igrad
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = output_ground_state_energy*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = output_ground_state_energy*hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = output_ground_state_gradient(iatom,idir)*hartree_per_bohr_to_ev_per_angstrom
                enddo
            enddo
        elseif( igrad == 2 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives above first order are unavailable')
        endif

        contains

        !> @brief Evaluate prepare_local_context for the FH2_LEPS_NO5a_1985 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            real(wp), parameter :: ckcau  = hartree_to_kilocalorie_per_mole
            real(wp), parameter :: cangau = 0.52917706e0_wp
            !   EXLARG IS LN(LARGEST FLOATING POINT NUMBER ALLOWED)
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            do i = 1, 3
                dm(i)   = dm(i) / ckcau
                re(i)   = re(i) / cangau
                beta(i) = beta(i) * cangau
            enddo
            !   Set the energy of the potential in the each of the three asymptotic
            easyab = dm(1)
            easybc = dm(2)
            easyac = dm(3)
            !   USEFUL CONSTANTS
            zhhdif   = zhhst - zhh0
            omp3c    = 1.e0_wp - p3c
            omq3c    = 1.e0_wp - q3c
            ha2      = 0.5e0_wp * a2
            g2       = noncollinear_width * noncollinear_width
            ezero(1) = dm(2)
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1) = 'R. Steckler, D. G. Truhlar, B. C. Garrett,'
            ref(2) = 'J. Chem. Phys. 82, 5499(1985)'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the FH2_LEPS_NO5a_1985 potential.
        subroutine pot
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the FH diatomic, with the
            !   FH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
            !   respectively.
            !   This potential is written such that:
            !                  R(1) = R(F-H)
            !                  R(2) = R(H-H)
            !                  R(3) = R(H-F)
            !   The zero of energy is defined at F "infinitely" far from the H2 diat
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: capr, capr2, cosphi, costh, cphi, dcos1, dcos2, dcos3, dcphi, dd2, de3c1, de3c2, de3c3, deg, degdif, dphi1, dphi2, dphi3, dzphi, e3c, ex, ex1, ex2, ex3, exx, hsech, op3z, phi, r12, r1i, r22, r32, r3i, rdif, rdif2, rsum, scale1, scale2, sgnphi, sinp2, sinphi, t, t1, t2, t3, t4, t5, xjx, zhh, zp1, zp3, zphi
            integer             :: i, i1, i2, i3, j
            real(wp), parameter :: tdeg   = 57.29577951e0_wp
            real(wp), parameter :: exlarg = 80.e0_wp
            logical             :: lcol, lderiv, lzero
            lderiv                                = .false.
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   Check the values of surface_requests and requested_derivative_order for validity.
            !   Check the value of requested_derivative_order
            !   The logical LDERIV is used to flag the derivative calculations
            !   in the code; LDERIV is set equal to .TRUE. if requested_derivative_order = 1.
            if( requested_derivative_order .eq. 1 ) lderiv = .true.
            !   TEST R'S FOR PHYSICAL VALUES
            do i = 1, 3
                rr(i) = r(i)
            enddo
            do i1 = 1, 3
                i2                         = mod(i1, 3) + 1
                i3                         = mod(i2, 3) + 1
                t                          = rr(i2) + rr(i3)
                if( rr(i1) .gt. t ) rr(i1) = t
            enddo
            r12 = rr(1) * rr(1)
            r22 = rr(2) * rr(2)
            r32 = rr(3) * rr(3)
            !   CAPR IS THE F TO H2 DISTANCE
            capr2 = 0.5e0_wp * (r12 + r32 - 0.5e0_wp * r22)
            if( capr2 .lt. 0.e0_wp ) capr2 = 0.0e0_wp
            capr = sqrt(capr2)
            !   PHI IS THE ANGLE BETWEEN CAPR AND RHH
            !   PHI IS RESTRICTED TO BE BETWEEN O AND 90 DEGREES
            t = 0.5e0_wp * (r32 - r12)
            if( abs(t) .gt. 1.e-7_wp ) t = t / (rr(2) * capr)
            t2     = 1.0e0_wp
            sgnphi = sign(t2, t)
            cosphi = sgnphi * t
            cosphi = min(cosphi, t2)
            phi    = acos(cosphi)
            sinp2  = 1.e0_wp - cosphi * cosphi
            if( sinp2 .lt. 0.e0_wp ) sinp2 = 0.e0_wp
            sinphi = sqrt(sinp2)
            !   LCOL IS SET TRUE IF COLLINEAR GEOMETRY
            lcol = abs(sinphi) .lt. 1.e-6_wp .or. abs(capr) .lt. 1.e-6_wp
            !   COMPUTE COLLINEAR SATO PARAMETERS
            z(1) = zfh0
            t    = a * (rr(2) - b)
            if( t < exlarg )then
                t     = exp(t)
                t1    = 1.0e0_wp/t
                hsech = 1.0e0_wp/(t+t1)
                zhh   = a1+ha2*((t-t1)*hsech+1.0e0_wp)
            else
                hsech = 0.0e0_wp
                zhh   = a1+a2
            endif
            t    = rr(2) - f
            exx  = c * exp(hh_sato_exponent * t * t)
            zhh  = zhh - exx
            z(2) = zhh
            goto_block_50_451: block
                if( .not.lderiv ) exit goto_block_50_451
                do i = 1, 3
                    do j = 1, 3
                        dz(i,j) = 0.e0_wp
                    enddo
                enddo
                !  ONLY THE HH SATO PARAMETER HAS A NONZERO DERIVATIVE - W.R.T. R2
                dz(2,2) = 2.e0_wp * a * a2 * hsech * hsech - 2.e0_wp * hh_sato_exponent * t * exx
            end block goto_block_50_451
            !   COMPUTE NONCOLLINEAR CONTRIBUTIONS TO SATO PARAMETERS
            goto_block_160_550: block
                if( lcol ) exit goto_block_160_550
                !   COMPUTE Z(PHI) FOR FH SATO PARAMETER
                !   DZPHI IS THE DERIVATIVE OF Z(PHI) W.R.T. PHI
                deg = tdeg * phi
                goto_block_60_465: block
                    if( .not.(deg .lt. 0.e0_wp .or. deg .gt. 90.e0_wp) ) exit goto_block_60_465
                    if( present(status) ) call status%set(status_invalid_input, &
                        'POT 3')
                    return
                end block goto_block_60_465
                if( deg <= deg1 )then
                    zphi  = 0.0e0_wp
                    dzphi = 0.0e0_wp
                elseif( deg <= deg2 )then
                    degdif = deg-deg1
                    dd2    = degdif*degdif
                    zphi   = dd2*(c1+c2*degdif)
                    if( lderiv ) dzphi = tdeg*degdif*(2.0e0_wp*c1+3.0e0_wp*c2*degdif)
                elseif( deg <= deg3 )then
                    degdif = deg-deg3p
                    zphi   = c3+c4*degdif
                    if( lderiv ) dzphi = c4*tdeg
                elseif( deg <= deg4 )then
                    zphi = c5+c6*sinphi*sinphi
                    if( lderiv ) dzphi = 2.0e0_wp*c6*sinphi*cosphi
                else
                    if( present(status) ) call status%set(status_invalid_input, &
                        'POT 4')
                    return
                endif
                !   COMPUTE SCALE FACTOR TO DAMP OUT CORRECTION AT LARGE DISTANCE
                scale1 = capr2 / (g2 + capr2)
                z(1)   = z(1) + scale1 * zphi
                !   COMPUTE SCALE FACTOR FOR HH SATO PARAMETER
                scale2 = (zhh - zhh0) / zhhdif
                !   COMPUTE PHI DEPENDENT PART OF HH SATO PARAMETER
                cphi = (bet + betp * sinp2) * sinp2
                !   COMPUTE HH SATO PARAMTER
                z(2) = z(2) + scale1 * scale2 * cphi
                goto_block_150_549: block
                    if( .not.lderiv ) exit goto_block_150_549
                    !   COMPUTE DERIVATIVES
                    !   FIRST, DERIVATIVES OF PHI W.R.T. R1, R2, R3
                    t     = 1.e0_wp / capr
                    t1    = 0.5e0_wp * cosphi * t
                    t2    = sgnphi / rr(2)
                    dphi1 = rr(1) * (t1 + t2)
                    dphi3 = rr(3) * (t1 - t2)
                    t1    = 1.e0_wp / sinphi
                    t2    = t1 * t
                    dphi1 = dphi1 * t2
                    dphi3 = dphi3 * t2
                    dphi2 = rr(2) * cosphi * t1 * (1.e0_wp / r22 - 0.25e0_wp / capr2)
                    t1    = (1.e0_wp - scale1) / capr2
                    t2    = t1 * zphi
                    !   NEXT, DERIVATIVES OF FH SATO PARAMETERS
                    dz(1,1) = scale1 * (t2 * rr(1) + dzphi * dphi1)
                    dz(1,2) = scale1 * (-0.5e0_wp * t2 * rr(2) + dzphi * dphi2)
                    dz(1,3) = scale1 * (t2 * rr(3) + dzphi * dphi3)
                    t1      = t1 * cphi
                    !   DERIVATIVE OF THE PHI DEPENDENT PART OF HH SATO W.R.T. PHI
                    dcphi = 2.e0_wp * (bet + 2.e0_wp * betp * sinp2) * sinphi * cosphi
                    t2    = scale1 * scale2
                    !   LAST, DERIVATIVES OF HH SATO PARAMETER
                    dz(2,1) = t2 * (t1 * rr(1) + dcphi * dphi1)
                    dz(2,2) = dz(2,2) * (1.e0_wp + scale1 * cphi / zhhdif) + t2 * (-0.5e0_wp * t1 * rr(2) + dcphi * &
                        & dphi2)
                    dz(2,3) = t2 * (t1 * rr(3) + dcphi * dphi3)
                end block goto_block_150_549
            end block goto_block_160_550
            !   COMPUTE POTENTIAL FROM LEPS FORM
            ground_state_energy = 0.e0_wp
            z(3)                = z(1)
            goto_block_170_567: block
                if( .not.lderiv ) exit goto_block_170_567
                do i = 1, 3
                    dz(3,i)                           = dz(1,i)
                    ground_state_internal_gradient(i) = 0.e0_wp
                    q(i)                              = 0.e0_wp
                    xj(i)                             = 0.e0_wp
                    do j = 1, 3
                        dq(i,j) = 0.e0_wp
                        dj(i,j) = 0.e0_wp
                    enddo
                enddo
            end block goto_block_170_567
            lzero = .true.
            do i = 1, 3
                ex    = exp(beta(i) * (re(i) - rr(i)))
                lzero = lzero .and. ex.eq.0.0e0_wp
                if( ex .eq. 0.0e0_wp ) cycle
                zp1  = z(i) + 1.e0_wp
                zp3  = zp1 + 2.e0_wp
                op3z = 1.e0_wp + 3.e0_wp * z(i)
                t    = 0.25e0_wp * dm(i) * ex / zp1
                t1   = zp3 * ex
                t2   = op3z * ex
                !   Q AND XJ ARE THE COULUMB AND EXCHANGE PARTS, INDEX I RUNS OVER THE
                !   COORDINATES OF WHICH THEY ARE EXPLICIT FUNCTIONS
                q(i)  = t * (t1 - 2.e0_wp  * op3z)
                xj(i) = t * (t2 - 2.e0_wp  * zp3)
                !   PUT SUM OF Q'S IN ENERGY
                ground_state_energy = ground_state_energy + q(i)
                goto_block_190_602: block
                    if( .not.lderiv ) exit goto_block_190_602
                    !   DERIVATIVE OF Q AND ZJ W.R.T. R1, R2, R3
                    !   ELEMENT I,J IS THE DERIVATIVE OF Q(I) W.R.T. RJ
                    t3 = 2.e0_wp * t * (ex + 2.e0_wp) / zp1
                    do j = 1, 3
                        dq(i,j) = -t3 * dz(i,j)
                        dj(i,j) = -dq(i,j)
                    enddo
                    t       = 2.e0_wp*t*beta(i)
                    dq(i,i) = dq(i,i) - t * (t1 - op3z)
                    dj(i,i) = dj(i,i) - t * (t2 - zp3)
                end block goto_block_190_602
            enddo
            goto_block_270_646: block
                if( lzero ) exit goto_block_270_646
                xjx = 0.e0_wp
                !   COMPUTE THE TOTAL EXCHANGE TERM
                do i1 = 1, 3
                    i2  = mod(i1, 3) + 1
                    t   = xj(i1) - xj(i2)
                    xjx = xjx + t * t
                    if( .not.lderiv ) cycle
                    !   PUT THE DERIVATIVES OF THE TOTAL EXCHANGE TERM IN ground_state_internal_gradient
                    do j = 1, 3
                        ground_state_internal_gradient(j) = ground_state_internal_gradient(j) + t * (dj(i1,j) - dj(i2, &
                            & j))
                    enddo
                enddo
                xjx = sqrt(0.5e0_wp * xjx)
                !   COMPUTE THE LEPS PART OF THE POTENTIAL ENERGY
                ground_state_energy = ground_state_energy - xjx
                goto_block_260_645: block
                    if( .not.lderiv ) exit goto_block_260_645
                    !   COMPUTE DERIVATIVES OF THE LEPS PART
                    if( abs(xjx) .gt. 1.e-14_wp ) t = 0.5e0_wp / xjx
                    do j = 1, 3
                        !   THE DERIVATIVE OF THE EXCHANGE PART MUST BE DIVIDED BY 2*XJX
                        ground_state_internal_gradient(j) = -t *ground_state_internal_gradient(j)
                        do i = 1, 3
                            !   THEN ADD IN THE DERIVATIVE OF Q
                            ground_state_internal_gradient(j) = ground_state_internal_gradient(j) + dq(i,j)
                        enddo
                    enddo
                end block goto_block_260_645
            end block goto_block_270_646
            !   EZERO ADDED TO PUT ZERO AT REQ FOR REACTANTS
            ground_state_energy = ground_state_energy + ezero(1)
            !   COMPUTE THE 3-CENTER TERM
            rsum  = rr(1) + rr(3)
            rdif  = rr(1) - rr(3)
            rdif2 = rdif * rdif
            t1    = -alfp * rsum
            ex1   = c3c * exp(t1 * rsum)
            t2    = -alf * rdif
            ex2   = omp3c * exp(t2 * rdif)
            t3    = -alft * rdif2
            ex3   = p3c * exp(t3 * rdif2)
            r1i   = 1.e0_wp / rr(1)
            r3i   = 1.e0_wp / rr(3)
            t4    = r1i*r3i
            !   COSTH IS THE ANGLE BETWEEN R1 AND R3
            costh = 0.5e0_wp * (r12 + r32 - r22) * t4
            t5    = omq3c * costh
            t     = q3c + t5 * costh
            !   E3C CONTAINS THE THREE CENTER CORRECTION TO THE ENERGY
            e3c                 = (ex2 + ex3) * t * ex1
            ground_state_energy = ground_state_energy + e3c
            goto_block_280_696: block
                if( .not.lderiv ) exit goto_block_280_696
                !   COMPUTE DERIVATIVE OF THE 3C TERM
                !   FIRST, DERIVATIVE OF COSTH
                dcos1 = r3i - costh*r1i
                dcos2 = - rr(2)*t4
                dcos3 = r1i - costh*r3i
                t4    = ex2 + ex3
                t5    = t4 * t5
                !   NOW, DERIVATIVES OF E3C
                de3c1                             = 2.e0_wp*((t2 * ex2 + 2.e0_wp * t3 * rdif * ex3 + t1 * t4) * t + &
                    & t5 * dcos1) * ex1
                de3c2 = 2.e0_wp * t5 * dcos2 * ex1
                de3c3                             = 2.e0_wp*((-t2 * ex2 - 2.e0_wp * t3 * rdif * ex3 + t1 * t4) * t + &
                    & t5 * dcos3) * ex1
                ground_state_internal_gradient(1) = ground_state_internal_gradient(1) + de3c1
                ground_state_internal_gradient(2) = ground_state_internal_gradient(2) + de3c2
                ground_state_internal_gradient(3) = ground_state_internal_gradient(3) + de3c3
            end block goto_block_280_696
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate CARTTOR for the FH2_LEPS_NO5a_1985 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(cartesian_coordinates, atom_count, coordinate_representation, r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the FH2_LEPS_NO5a_1985 potential.
        subroutine eunitzero
            implicit none
            integer :: i, j
            output_ground_state_energy = ground_state_energy
            if( first_requested_excited_state .ne. 0 )then
                do i=first_requested_excited_state, last_requested_excited_state
                    output_excited_state_energies(i) = excited_state_energies(i)
                enddo
            endif
            if( first_requested_coupling .ne. 0 )then
                do j=first_requested_coupling, last_requested_coupling
                    output_diabatic_couplings(j) = diabatic_couplings(j)
                enddo
            endif
        end subroutine eunitzero

        !> @brief Evaluate RTOCART for the FH2_LEPS_NO5a_1985 potential.
        subroutine rtocart
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            real(wp)           :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, &
                & termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes(n3atom, isurf), &
                & ygs(n3atom), yij(n3atom, jsurf), za, zb, zc, zd
            integer :: i, ind, j, j1, j2, k, m1
            if( coordinate_representation .eq. 1 )then
                do i = 1, atom_count
                    ind                                  = 3*i-2
                    ground_state_cartesian_gradient(i,1) = ground_state_internal_gradient(ind)
                    ground_state_cartesian_gradient(i,2) = ground_state_internal_gradient(ind+1)
                    ground_state_cartesian_gradient(i,3) = ground_state_internal_gradient(ind+2)
                    if( first_requested_excited_state .ne. 0 )then
                        do j = first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(i,1,j) = excited_state_internal_gradients(ind,j)
                            excited_state_cartesian_gradients(i,2,j) = excited_state_internal_gradients(ind+1,j)
                            excited_state_cartesian_gradients(i,3,j) = excited_state_internal_gradients(ind+2,j)
                        enddo
                    endif
                    if( last_requested_coupling .ne. 0 )then
                        do k = first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(i,1,k) = diabatic_internal_gradients(ind,k)
                            diabatic_cartesian_gradients(i,2,k) = diabatic_internal_gradients(ind+1,k)
                            diabatic_cartesian_gradients(i,3,k) = diabatic_internal_gradients(ind+2,k)
                        enddo
                    endif
                enddo
            elseif( coordinate_representation .eq. 2 )then
                do i = 1, atom_count
                    ground_state_cartesian_gradient(i,1) = 0.e0_wp
                    ground_state_cartesian_gradient(i,2) = 0.e0_wp
                    ground_state_cartesian_gradient(i,3) = 0.e0_wp
                    if( first_requested_excited_state .ne. 0 )then
                        do j1=first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(i,1,j1) = 0.e0_wp
                            excited_state_cartesian_gradients(i,2,j1) = 0.e0_wp
                            excited_state_cartesian_gradients(i,3,j1) = 0.e0_wp
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do j2=first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(i,1,j2) = 0.e0_wp
                            diabatic_cartesian_gradients(i,2,j2) = 0.e0_wp
                            diabatic_cartesian_gradients(i,3,j2) = 0.e0_wp
                        enddo
                    endif
                    do j = 1, atom_count
                        if( j .lt. i )then
                            m1 = atom_count*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
                            m1 = atom_count*(i-1) - (i*(i-1))/2 + j-i
                        else
                            cycle
                        endif
                        y = ground_state_internal_gradient(m1)
                        termx                                = (cartesian_coordinates(i,1)-cartesian_coordinates(j, &
                            & 1))/r(m1)
                        termy                                = (cartesian_coordinates(i,2)-cartesian_coordinates(j, &
                            & 2))/r(m1)
                        termz                                = (cartesian_coordinates(i,3)-cartesian_coordinates(j, &
                            & 3))/r(m1)
                        ground_state_cartesian_gradient(i,1) = ground_state_cartesian_gradient(i,1) + termx*y
                        ground_state_cartesian_gradient(i,2) = ground_state_cartesian_gradient(i,2) + termy*y
                        ground_state_cartesian_gradient(i,3) = ground_state_cartesian_gradient(i,3) + termz*y
                        if( first_requested_excited_state .gt. 0 )then
                            y = excited_state_internal_gradients(m1,j1)
                            do j1=first_requested_excited_state, last_requested_excited_state
                                excited_state_cartesian_gradients(i,1,j1) = excited_state_cartesian_gradients(i,1, &
                                    & j1) + termx*y
                                excited_state_cartesian_gradients(i,2,j1) = excited_state_cartesian_gradients(i,2, &
                                    & j1) + termy*y
                                excited_state_cartesian_gradients(i,3,j1) = excited_state_cartesian_gradients(i,3, &
                                    & j1) + termz*y
                            enddo
                        elseif( first_requested_coupling .gt. 0 )then
                            do j2=first_requested_coupling, last_requested_coupling
                                y                                    = diabatic_internal_gradients(m1,j2)
                                diabatic_cartesian_gradients(i,1,j2) = diabatic_cartesian_gradients(i,1,j2) + termx*y
                                diabatic_cartesian_gradients(i,2,j2) = diabatic_cartesian_gradients(i,2,j2) + termy*y
                                diabatic_cartesian_gradients(i,3,j2) = diabatic_cartesian_gradients(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            elseif( coordinate_representation .eq. 3 )then
                do i = 1, atom_count
                    ygs(i) = ground_state_internal_gradient(i)/r(i)
                    if( first_requested_excited_state .ne. 0 )then
                        do j=first_requested_excited_state, last_requested_excited_state
                            yes(i, j) = excited_state_internal_gradients(i,j)/r(i)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do k=first_requested_coupling, last_requested_coupling
                            yij(i, k) = diabatic_internal_gradients(i,k)/r(i)
                        enddo
                    endif
                enddo
                do k = 1, 3
                    term12                               = cartesian_coordinates(1,k)-cartesian_coordinates(2,k)
                    term23                               = cartesian_coordinates(2,k)-cartesian_coordinates(3,k)
                    term13                               = cartesian_coordinates(1,k)-cartesian_coordinates(3,k)
                    ground_state_cartesian_gradient(1,k) = term12*ygs(1) + term13*ygs(3)
                    ground_state_cartesian_gradient(2,k) = -term12*ygs(1) + term23*ygs(2)
                    ground_state_cartesian_gradient(3,k) = -term13*ygs(3) - term23*ygs(2)
                    if( first_requested_excited_state .ne. 0 )then
                        do j1=first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(1,k,j1) = term12*yes(1, j1) + term13*yes(3, j1)
                            excited_state_cartesian_gradients(2,k,j1) = -term12*yes(1, j1) + term23*yes(2, j1)
                            excited_state_cartesian_gradients(3,k,j1) = -term13*yes(3, j1) - term23*yes(2, j1)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do j2=first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(1,k,j2) = term12*yij(1, j2) + term13*yij(3, j2)
                            diabatic_cartesian_gradients(2,k,j2) = -term12*yij(1, j2) + term23*yij(2, j2)
                            diabatic_cartesian_gradients(3,k,j2) = -term13*yij(3, j2) - term23*yij(2, j2)
                        enddo
                    endif
                enddo
            elseif( coordinate_representation .eq. 4 )then
                wh   = 1.007825e0_wp
                wf   = 18.99840e0_wp
                sum  = wh+wf
                eps  = wf/sum
                epsp = wh/sum
                u1   = cos(theta1)
                u2   = cos(theta2)
                u3   = cos(phi)
                ss1  = sin(theta1)
                ss2  = sin(theta2)
                ss3  = sin(phi)
                ya   = 0.0e0_wp
                yb   = 0.0e0_wp
                t0   = r1*u1
                za   = -epsp*t0
                zb   = eps*t0
                t0   = r1*ss1
                xa   = -epsp*t0
                xbb  = eps*t0
                t0   = r2*ss2
                t1   = t0*u3
                xc   = -epsp*t1
                xd   = eps*t1
                t1   = t0*ss3
                yc   = -epsp*t1
                yd   = eps*t1
                t0   = r2*u2
                zc   = -epsp*t0+rcm
                zd   = eps*t0+rcm
                rff  = sqrt((xa-xc)**2+yc**2+(za-zc)**2)
            else
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the FH2_LEPS_NO5a_1985 potential.
        subroutine dedcou
            implicit none
            integer :: i, j, k
            do i = 1, atom_count
                output_ground_state_gradient(i,1) = ground_state_cartesian_gradient(i,1)
                output_ground_state_gradient(i,2) = ground_state_cartesian_gradient(i,2)
                output_ground_state_gradient(i,3) = ground_state_cartesian_gradient(i,3)
                if( first_requested_excited_state .ne. 0 )then
                    do j=first_requested_excited_state, last_requested_excited_state
                        output_excited_state_gradients(i,1,j) = excited_state_cartesian_gradients(i,1,j)
                        output_excited_state_gradients(i,2,j) = excited_state_cartesian_gradients(i,2,j)
                        output_excited_state_gradients(i,3,j) = excited_state_cartesian_gradients(i,3,j)
                    enddo
                endif
                if( first_requested_coupling .ne. 0 )then
                    do k=first_requested_coupling, last_requested_coupling
                        output_diabatic_gradients(i,1,k) = diabatic_cartesian_gradients(i,1,k)
                        output_diabatic_gradients(i,2,k) = diabatic_cartesian_gradients(i,2,k)
                        output_diabatic_gradients(i,3,k) = diabatic_cartesian_gradients(i,3,k)
                    enddo
                endif
            enddo
        end subroutine dedcou
    end subroutine evaluate
end module chempotpy_fh2_leps_no5a_1985_pes

!> @brief Expose FH2_LEPS_NO5a_1985 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_fh2_leps_no5a_1985_pes, only: evaluate
    use chempotpy_kinds,                  only: wp
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
