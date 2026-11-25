!   System:          OH2
!   Functional form: Anti-Morse Bends plus Bowman's collinear spline
!   Common name:     OH2DIM
!   Interface:       potlib2001
!   Reference:
!   Number of bodies:              3
!   Number of derivatives:         1
!   Number of electronic surfaces: 1
!   prepare_local_context must be called once before any calls to POT.
!   The control flags for the potential, such as IPRT, are
!   initialized in the block data subprogram PTPACM.
!   Coordinates, potential energy, and derivatives are passed
!   The potential energy in the three asymptotic valleys are
!   stored in the common block ASYCM:
!   The potential energy in the AB valley, EASYAB, is equal to the poten
!   energy of the H "infinitely" far from the OH diatomic, with the
!   OH diatomic at its equilibrium configuration.  Similarly, the terms
!   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
!   respectively.
!   All the information passed through the common blocks PT1CM and ASYCM
!   is in Hartree atomic units.
!   This potential is written such that:
!                  R(1) = R(O-H)
!                  R(2) = R(H-H)
!                  R(3) = R(H-O)
!   The zero of energy is defined at O "infinitely" far from the H2 diat
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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          GENERAL INFORMATION
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   POTA.FOR HAS ANTI-MORSE BEND PLUS BOWMAN'S COLLINEAR SPLINE ROUTINE.
!   ANTI-MORSE BENDING POTENTIAL IS ADDED TO COLLINEAR POTENTIAL.
!   THIS SUBROUTINE CAN BE USED TO ADD A BENDING CORRECTION TO
!   ANY COLLINEAR POTENTIAL SUBROUTINE.
!CCCCCCCCCCCCCCCCCCCCCCCCCC
!          DEFINE VARIABLES
!CCCCCCCCCCCCCCCCCCCCCCCCCC

!> @brief Implement OH2_GEN_DIM1981, a potential energy surface for OH2.
module chempotpy_oh2_gen_dim1981_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the OH2_GEN_DIM1981 potential energy surface for OH2.
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
        real(wp)              :: asycm_easyab
        real(wp)              :: asycm_easyac
        real(wp)              :: asycm_easybc
        real(wp)              :: cartesian_coordinates(natom,3)
        integer               :: coordinate_representation
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: atom_count
        character(len=75)     :: infocm_ref(5)
        real(wp)              :: nstep_step
        real(wp)              :: potcm_a
        real(wp)              :: potcm_aab(3)
        real(wp)              :: potcm_abc(3)
        real(wp)              :: potcm_beta
        real(wp)              :: potcm_bsp(40,3)
        real(wp)              :: potcm_cab(3)
        real(wp)              :: potcm_cbc(3)
        real(wp)              :: potcm_dbcoor(2)
        real(wp)              :: potcm_de
        real(wp)              :: potcm_gamma
        integer               :: potcm_nphi
        real(wp)              :: potcm_phi(40)
        real(wp)              :: potcm_r1s
        real(wp)              :: potcm_r2s
        real(wp)              :: potcm_rbb
        real(wp)              :: potcm_req(3)
        real(wp)              :: potcm_rhh
        real(wp)              :: potcm_vsp
        real(wp)              :: potcm_xab(3)
        real(wp)              :: potcm_xbc(3)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: ground_state_energy
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: input_cartesian(natom,3)
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        asycm_easyab                      = 0.e0_wp; asycm_easyac = 0.e0_wp; asycm_easybc = 0.e0_wp; cartesian_coordinates = 0.e0_wp; coordinate_representation = 0
        maximum_derivative_order          = 0; maximum_excited_state = 0; atom_count = 0
        infocm_ref                        = ' '
        nstep_step                        = 0.e0_wp
        potcm_a                           = 0.e0_wp
        potcm_aab                         = 0.e0_wp
        potcm_abc                         = 0.e0_wp
        potcm_beta                        = 0.e0_wp
        potcm_bsp                         = 0.e0_wp
        potcm_cab                         = 0.e0_wp
        potcm_cbc                         = 0.e0_wp
        potcm_dbcoor                      = 0.e0_wp
        potcm_de                          = 0.e0_wp
        potcm_gamma                       = 0.e0_wp
        potcm_nphi                        = 0
        potcm_phi                         = 0.e0_wp
        potcm_r1s                         = 0.e0_wp
        potcm_r2s                         = 0.e0_wp
        potcm_rbb                         = 0.e0_wp
        potcm_req                         = 0.e0_wp
        potcm_rhh                         = 0.e0_wp
        potcm_vsp                         = 0.e0_wp
        potcm_xab                         = 0.e0_wp
        potcm_xbc                         = 0.e0_wp
        ground_state_internal_gradient    = 0.e0_wp
        ground_state_energy               = 0.e0_wp
        pt1cm_r                           = 0.e0_wp
        pt3cm_ezero                       = 0.e0_wp
        excited_state_internal_gradients  = 0.e0_wp
        excited_state_energies            = 0.e0_wp
        diabatic_internal_gradients       = 0.e0_wp
        diabatic_couplings                = 0.e0_wp
        input_cartesian                   = 0.e0_wp
        surface_requests                  = 0
        requested_derivative_order        = 0
        output_excited_state_gradients    = 0.e0_wp
        output_ground_state_gradient      = 0.e0_wp
        output_diabatic_gradients         = 0.e0_wp
        output_excited_state_energies     = 0.e0_wp
        output_ground_state_energy        = 0.e0_wp
        output_diabatic_couplings         = 0.e0_wp
        excited_state_cartesian_gradients = 0.e0_wp
        ground_state_cartesian_gradient   = 0.e0_wp
        diabatic_cartesian_gradients      = 0.e0_wp
        last_requested_excited_state      = 0
        last_requested_coupling           = 0
        first_requested_excited_state     = 0
        first_requested_coupling          = 0
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 0
        coordinate_representation  = 3
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 3
        nstep_step                 = 0.000000008e0_wp
        potcm_nphi                 = 17
        potcm_r1s                  = 0.2212079e+01_wp
        potcm_r2s                  = 0.2215309e+01_wp
        potcm_vsp                  = 0.1094890e+03_wp
        potcm_xbc                  = [0.1094890e+03_wp, 0.1473409e+01_wp, 0.1941998e+01_wp]
        potcm_xab                  = [0.1056400e+03_wp, 0.1241479e+01_wp, 0.2304798e+01_wp]
        potcm_cbc                  = [-0.4815674e-01_wp, 0.3814697e-05_wp, 0.4415512e-03_wp]
        potcm_abc                  = [0.7873216e+01_wp, 0.7873216e+01_wp, 0.7873216e+01_wp]
        potcm_cab                  = [-0.1482804e+01_wp, -0.2084732e-02_wp, -0.3727913e-02_wp]
        potcm_aab                  = [0.1932070e+01_wp, 0.1932070e+01_wp, 0.1932070e+01_wp]
        potcm_phi(1:17) = [0.0000000e+00_wp, 0.2000000e+01_wp, 0.4000000e+01_wp, 0.6000000e+01_wp, 0.1800000e+02_wp, &
            & 0.2400000e+02_wp, 0.3000000e+02_wp, 0.3400000e+02_wp, 0.3800000e+02_wp, 0.4400000e+02_wp, &
            & 0.4800000e+02_wp, 0.5400000e+02_wp, 0.6000000e+02_wp, 0.7600000e+02_wp, 0.8600000e+02_wp, &
            & 0.8800000e+02_wp, 0.9000000e+02_wp]
        potcm_bsp(1:17,1) = [0.1094409e+03_wp, 0.1094212e+03_wp, 0.1094004e+03_wp, 0.1093775e+03_wp, 0.1088412e+03_wp, &
            & 0.1077417e+03_wp, 0.1055117e+03_wp, 0.1030257e+03_wp, 0.9966266e+02_wp, 0.9630693e+02_wp, &
            & 0.9637480e+02_wp, 0.9811765e+02_wp, 0.1000138e+03_wp, 0.1029822e+03_wp, 0.1038929e+03_wp, &
            & 0.1040309e+03_wp, 0.1041572e+03_wp]
        potcm_bsp(1:17,2) = [0.1473413e+01_wp, 0.1474257e+01_wp, 0.1476847e+01_wp, 0.1481213e+01_wp, 0.1545280e+01_wp, &
            & 0.1599182e+01_wp, 0.1655956e+01_wp, 0.1681611e+01_wp, 0.1680401e+01_wp, 0.1626349e+01_wp, &
            & 0.1572269e+01_wp, 0.1486160e+01_wp, 0.1408802e+01_wp, 0.1273705e+01_wp, 0.1241724e+01_wp, &
            & 0.1239825e+01_wp, 0.1239394e+01_wp]
        potcm_bsp(1:17,3) = [0.1942440e+01_wp, 0.1940545e+01_wp, 0.1936212e+01_wp, 0.1929558e+01_wp, 0.1825750e+01_wp, &
            & 0.1719964e+01_wp, 0.1603891e+01_wp, 0.1538857e+01_wp, 0.1490752e+01_wp, 0.1525683e+01_wp, &
            & 0.1620584e+01_wp, 0.1788983e+01_wp, 0.1942842e+01_wp, 0.2222638e+01_wp, 0.2294156e+01_wp, &
            & 0.2299114e+01_wp, 0.2301070e+01_wp]
        potcm_de                   = 106.56e0_wp
        potcm_beta                 = 2.07942e0_wp
        potcm_req                  = [0.96966e0_wp, 0.74144e0_wp, 0.96966e0_wp]
        potcm_a                    = 0.26e0_wp
        potcm_gamma                = 0.4131e0_wp
        potcm_rbb                  = 0.74144e0_wp
        v                          = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        requested_derivative_order = igrad
        surface_requests           = 1
        input_cartesian            = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                input_cartesian(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
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

        !> @brief Evaluate prepare_local_context for the OH2_GEN_DIM1981 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, it
            ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
            real(wp), parameter :: bohr = .52917706e0_wp
            ! Conversion constant from Hartree to kcal (kcal/Hartree)
            real(wp), parameter :: hrtree = hartree_to_kilocalorie_per_mole
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            call prepare_local_context_2
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !      BENDING POTENTIAL PARAMETERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !   Echo the potential parameters to the file linked to UNIT NFLAG(18)
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !       CONVERT TO ATOMIC UNITS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            do it = 1, 3
                potcm_req(it) = potcm_req(it)/bohr
            enddo
            potcm_beta = potcm_beta * bohr
            potcm_de   = potcm_de/hrtree
            potcm_a    = potcm_a/bohr
            potcm_rhh  = 1.4007977e0_wp
            potcm_rbb  = potcm_rbb /bohr
            !CCCCCCCCCCC
            pt3cm_ezero(1) = potcm_vsp
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            infocm_ref(1) = 'No Reference'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the OH2_GEN_DIM1981 potential.
        subroutine pot
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the OH diatomic, with the
            !   OH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
            !   respectively.
            !   This potential is written such that:
            !                  R(1) = R(O-H)
            !                  R(2) = R(H-H)
            !                  R(3) = R(H-O)
            !   The zero of energy is defined at O "infinitely" far from the H2 diat
            !CCCCCCCCCCCCCCC
            !CCCCCCCCCCCCCCC
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: bcoor, c, dac, ractmp, rob, stuff, vabbc, vac
            integer  :: it
            real(wp) :: nab, nab1, nbc, nbc1, junk(4)
            ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
            ! Conversion constant from Hartree to kcal (kcal/Hartree)
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   Check the values of surface_requests and requested_derivative_order for validity.
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        FIRST CALCULATE SOME USEFUL NUMBERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
            !        AB AND BC RESPECTIVELY. (EQUATION 2.2)
            nab = exp((potcm_req(1) - pt1cm_r(1))/potcm_a)
            if( nab .lt. 1.e-15_wp ) nab = 1.e-15_wp
            nbc = exp((potcm_req(2) - pt1cm_r(2))/potcm_a)
            if( nbc .lt. 1.e-15_wp ) nbc = 1.e-15_wp
            !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
            !       SEE EQUATION 2.7A, 2.7B
            c = 1.0e0_wp-nab/nbc
            if( abs(c) <= 1.0e-14_wp )then
                nab1 = 0.5e0_wp
                nbc1 = 0.5e0_wp
            else
                c    = c/nab
                nab1 = (2.0e0_wp+c-sqrt(4.0e0_wp+c*c))/(2.0e0_wp*c)
                nbc1 = 1.0e0_wp-nab1
            endif
            !          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
            !          FIRST TRAP OUT ANY ZERO ARGUEMENTS
            if( nab1*nbc1 > 0.0e0_wp )then
                stuff = potcm_a * log(nab1*nbc1)
                rob   = (potcm_rbb - stuff)/(potcm_rhh - stuff)
            else
                rob = 1.0e0_wp
            endif
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        CALCULATE BENDING CORRECTION
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        CALCULATE V(R(3))
            junk(3) = exp(-potcm_beta*(pt1cm_r(3)-potcm_req(3))/rob)
            vac     = potcm_gamma *1.0e-14_wp * potcm_de * junk(3)* (1.e0_wp + 0.5e0_wp*junk(3))
            !        CALCULATE V(R(1)+R(2))
            junk(4) = exp(potcm_beta*(potcm_req(3)-pt1cm_r(2)-pt1cm_r(1))/rob)
            vabbc   = potcm_gamma *1.0e-14_wp * potcm_de * junk(4) * (1.e0_wp + 0.5e0_wp*junk(4))
            !        HERE IS THE BENDING CORRECTION
            bcoor = (vac - vabbc)*1.0e14_wp
            !          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
            !          OF BCOOR WITH RESPECT TO R(1), R(2), AND R(3).
            !        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
            if( requested_derivative_order .eq. 1 )then
                dac             = -potcm_gamma*potcm_de*potcm_beta*junk(3)*(1.e0_wp+junk(3))/rob
                potcm_dbcoor(1) = potcm_gamma * potcm_de * potcm_beta * junk(4) * (1.e0_wp + junk(4))/rob
                potcm_dbcoor(2) = potcm_dbcoor(1)
            endif
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        NOW CALCULATE THE COLLINEAR ENERGY
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR
            !          GEOMETRY AND CALL POT2
            ractmp     = pt1cm_r(3)
            pt1cm_r(3) = pt1cm_r(2) + pt1cm_r(1)
            call pot2
            if( requested_derivative_order .eq. 1 )then
                ground_state_internal_gradient(1) = ground_state_internal_gradient(1) + &
                    & ground_state_internal_gradient(3)
                ground_state_internal_gradient(2) = ground_state_internal_gradient(2) + &
                    & ground_state_internal_gradient(3)
                ground_state_internal_gradient(3) = dac
            endif
            pt1cm_r(3)          = ractmp
            ground_state_energy = ground_state_energy + bcoor
            if( requested_derivative_order .eq. 1 )then
                do it = 1, 2
                    ground_state_internal_gradient(it) = potcm_dbcoor(it) + ground_state_internal_gradient(it)
                enddo
            endif
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate prepare_local_context_2 for the OH2_GEN_DIM1981 potential.
        subroutine prepare_local_context_2
            !   POTENTIAL FUNCTION FOR RATE1D SERIES
            !   RMO-SPLINE TYPE FIT WITH BEBO TYPE CONTINUATION INTO THE ASYMBTOTIC
            !   POTENTIAL ROUTINE INTERACTS WITH REST OF PROGRAM WITH UNITS OF EV/BO
            !   POTENTIAL ROUTINE ASSUMES SPLINE INFO IS IN UNITS OF KCAL/ANGSTROM
            use chempotpy_kinds, only: wp
            implicit none
            ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
            ! Conversion constant from Hartree to kcal (kcal/Hartree)
            real(wp), parameter :: hrtree = hartree_to_kilocalorie_per_mole
            character(len=80)   :: title
            title = ' RMOFIT to DIM surface with BEBO continuations '
            !     *               XAB(2), XBC(2), XAB(2),
            !     *               XAB(3), XBC(3), XAB(3)
            !   Initialize the energy in the asymptotic valleys
            asycm_easyab = potcm_xab(1)/hrtree
            asycm_easybc = potcm_xbc(1)/hrtree
            asycm_easyac = asycm_easyab
        end subroutine prepare_local_context_2

        !> @brief Evaluate POT2 for the OH2_GEN_DIM1981 potential.
        subroutine pot2
            ! COCK SPLINES
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: angle, bet, betpr, betrg, csp(48,3), dbc, depr, derg, dis, dr1, dr2, ppp, r1, r2, r3, repr, req2, &
                & rerg, rr, scrap(48), sep, v(4), x
            integer             :: ic, k
            real(wp), parameter :: detrad = .0174532e0_wp
            do ic = 1, 4
                select case( ic )
                case( 1 )
                    r1 = pt1cm_r(1)
                    r2 = pt1cm_r(2)
                    r3 = pt1cm_r(3)
                case( 2 )
                    r1 = pt1cm_r(1)+nstep_step
                    r2 = pt1cm_r(2)
                    r3 = pt1cm_r(3)
                case( 3 )
                    r1 = pt1cm_r(1)
                    r2 = pt1cm_r(2)+nstep_step
                    r3 = pt1cm_r(3)
                case( 4 )
                    r1 = pt1cm_r(1)
                    r2 = pt1cm_r(2)
                    r3 = pt1cm_r(3)+nstep_step
                end select
                r1    = r1*0.529177e0_wp
                r2    = r2*0.529177e0_wp
                r3    = r3*0.529177e0_wp
                dr1   = potcm_r1s-r1
                dr2   = potcm_r2s-r2
                angle = atan(dr1/dr2)/detrad
                do k = 1, 3
                    ppp = spline(1, potcm_nphi, potcm_phi, potcm_bsp(1,k), csp(1,k), scrap, angle)
                enddo
                dbc   = potcm_xbc(1)
                derg  = potcm_xbc(1)/23.061e0_wp
                rerg  = (potcm_r2s-potcm_xbc(2))/.529177e0_wp
                betrg = potcm_xbc(3)*.529177e0_wp
                depr  = potcm_xab(1)/23.061e0_wp
                repr  = (potcm_r1s-potcm_xab(2))/.529177e0_wp
                betpr = potcm_xab(3)*.529177e0_wp
                if( r1 > potcm_r1s .and. r2 > potcm_r2s )then
                    ground_state_energy = potcm_xbc(1)/hartree_to_kilocalorie_per_mole
                    if( requested_derivative_order == 1 )then
                        ground_state_internal_gradient(1) = 1.0e-9_wp
                        ground_state_internal_gradient(2) = 1.0e-9_wp
                        ground_state_internal_gradient(3) = 0.0e0_wp
                    endif
                    return
                elseif( r1 > potcm_r1s )then
                    rr   = potcm_r2s-r2
                    sep  = r1-potcm_r1s
                    dis  = potcm_xbc(1)+potcm_cbc(1)*exp(-min(potcm_abc(1)*sep, 25.0e0_wp))
                    req2 = potcm_xbc(2)+potcm_cbc(2)*exp(-min(potcm_abc(2)*sep, 25.0e0_wp))
                    bet  = potcm_xbc(3)+potcm_cbc(3)*exp(-min(potcm_abc(3)*sep, 25.0e0_wp))
                elseif( r2 > potcm_r2s )then
                    rr   = potcm_r1s-r1
                    sep  = r2-potcm_r2s
                    dis  = potcm_xab(1)+potcm_cab(1)*exp(-min(potcm_aab(1)*sep, 25.0e0_wp))
                    req2 = potcm_xab(2)+potcm_cab(2)*exp(-min(potcm_aab(2)*sep, 25.0e0_wp))
                    bet  = potcm_xab(3)+potcm_cab(3)*exp(-min(potcm_aab(3)*sep, 25.0e0_wp))
                else
                    dr1   = potcm_r1s-r1
                    dr2   = potcm_r2s-r2
                    angle = atan(dr1/dr2)/detrad
                    rr    = sqrt(dr1*dr1+dr2*dr2)
                    dis   = spline(2, potcm_nphi, potcm_phi, potcm_bsp(1,1), csp(1,1), scrap, angle)
                    req2  = spline(2, potcm_nphi, potcm_phi, potcm_bsp(1,2), csp(1,2), scrap, angle)
                    bet   = spline(2, potcm_nphi, potcm_phi, potcm_bsp(1,3), csp(1,3), scrap, angle)
                endif
                x     = dis*((1.0e0_wp-exp(bet*(rr-req2)))**2-1.0e0_wp)+pt3cm_ezero(1)
                v(ic) = x
            enddo
            v(1:3)              = v(1:3)/hartree_to_kilocalorie_per_mole
            ground_state_energy = v(1)
            if( requested_derivative_order == 1 )then
                ground_state_internal_gradient(1) = (v(2)-v(1))/nstep_step
                ground_state_internal_gradient(2) = (v(3)-v(1))/nstep_step
                ground_state_internal_gradient(3) = 0.0e0_wp
            endif
        end subroutine pot2

        !> @brief Fit or evaluate the legacy cubic-spline representation.
        real(wp) function spline( isw, nn, x, y, c, d, xpt )
            use chempotpy_kinds, only: wp
            implicit none
            integer, intent(in)     :: isw, nn
            real(wp), intent(in)    :: x(nn), y(nn), xpt
            real(wp), intent(inout) :: c(nn+2), d(nn+2)
            integer                 :: i, j, n, np1, np2
            real(wp)                :: cfact, dfact, pivot, supd, v, v2, z
            real(wp)                :: deriv = 0.0e0_wp
            n   = nn
            np1 = n+1
            np2 = n+2
            z   = xpt
            select case( isw )
            case( 1 )
                c(1)   = y(1)
                d(1)   = 1.0e0_wp
                c(np1) = 0.0e0_wp
                d(np1) = 0.0e0_wp
                c(np2) = 0.0e0_wp
                d(np2) = 0.0e0_wp
                do i = 2, n
                    c(i) = y(i)-y(1)
                    d(i) = x(i)-x(1)
                enddo
                do i = 3, np2
                    if( d(i-1) == 0.0e0_wp )then
                        if( present(status) ) call status%set(status_invalid_input, &
                            'SPLINE: singular knot system')
                        return
                    endif
                    pivot = 1.0e0_wp/d(i-1)
                    if( i < np2 )then
                        supd = x(i-1)-x(i-2)
                        if( supd < 0.0e0_wp )then
                            if( present(status) ) call status%set(status_invalid_input, &
                                'SPLINE: knots must be ordered')
                            return
                        endif
                        supd = supd**3
                    else
                        supd = 1.0e0_wp
                    endif
                    dfact = supd*pivot
                    cfact = c(i-1)*pivot
                    if( i <= n )then
                        do j = i, n
                            v    = x(j)-x(i-2)
                            c(j) = c(j)-d(j)*cfact
                            d(j) = v**3-d(j)*dfact
                        enddo
                    endif
                    if( i < np2 )then
                        c(np1) = c(np1)-d(np1)*cfact
                        d(np1) = 1.0e0_wp-d(np1)*dfact
                    endif
                    c(np2) = c(np2)-d(np2)*cfact
                    d(np2) = x(i-2)-d(np2)*dfact
                enddo
                do i = 1, n
                    j = np2-i
                    if( j == np1 )then
                        v = 1.0e0_wp
                    else
                        v = (x(j)-x(j-1))**3
                    endif
                    if( d(j+1) == 0.0e0_wp )then
                        if( present(status) ) call status%set(status_invalid_input, &
                            'SPLINE: singular back substitution')
                        return
                    endif
                    c(j+1) = c(j+1)/d(j+1)
                    c(j)   = c(j)-c(j+1)*v
                enddo
                if( d(2) == 0.0e0_wp )then
                    if( present(status) ) call status%set(status_invalid_input, &
                        'SPLINE: singular leading coefficient')
                    return
                endif
                c(2)   = c(2)/d(2)
                spline = 0.0e0_wp
            case( 2 )
                spline = c(1)+c(2)*(z-x(1))
                do i = 1, n
                    v = z-x(i)
                    if( v <= 0.0e0_wp ) exit
                    spline = spline+c(i+2)*v**3
                enddo
            case( 3 )
                spline = c(1)+c(2)*(z-x(1))
                deriv  = c(2)/3.0e0_wp
                do i = 1, n
                    v = z-x(i)
                    if( v <= 0.0e0_wp ) exit
                    v2     = v*v
                    spline = spline+c(i+2)*v2*v
                    deriv  = deriv+c(i+2)*v2
                enddo
            case( 4 )
                spline = 3.0e0_wp*deriv
            case default
                if( present(status) ) call status%set(status_invalid_input, &
                    'SPLINE: invalid operation selector')
                return
            end select
        end function spline

        !> @brief Evaluate POTINFO for the OH2_GEN_DIM1981 potential.
        subroutine potinfo
            implicit none
            integer            :: i, j, k, kmax
            integer, parameter :: isurf = 5
            do i =1, 5
            enddo
            kmax = 0
            do i = 1, isurf+1
                do j = 1, isurf+1
                    if( surface_requests(i,j) .ne. 0.and.kmax .lt. max(i, j) ) kmax = max(i, j)
                enddo
            enddo
            if( kmax-1 .gt. maximum_excited_state )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            first_requested_excited_state = 0
            last_requested_excited_state  = 0
            do i = 2, isurf+1
                if( surface_requests(i,i) .ne. 0 )then
                    last_requested_excited_state = i-1
                    if( first_requested_excited_state .eq. 0 ) first_requested_excited_state = i-1
                endif
            enddo
            first_requested_coupling = 0
            last_requested_coupling  = 0
            k                        = 0
            do i = 1, isurf
                do j = i+1, isurf+1
                    k = k+1
                    if( surface_requests(i,j)+surface_requests(j,i) .ne. 0 )then
                        last_requested_coupling = k
                        if( first_requested_coupling .eq. 0 ) first_requested_coupling = k
                    endif
                enddo
            enddo
            if( requested_derivative_order .gt. maximum_derivative_order )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine potinfo

        !> @brief Evaluate CARTTOR for the OH2_GEN_DIM1981 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(cartesian_coordinates, atom_count, coordinate_representation, pt1cm_r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the OH2_GEN_DIM1981 potential.
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

        !> @brief Evaluate RTOCART for the OH2_GEN_DIM1981 potential.
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
                            & 1))/pt1cm_r(m1)
                        termy                                = (cartesian_coordinates(i,2)-cartesian_coordinates(j, &
                            & 2))/pt1cm_r(m1)
                        termz                                = (cartesian_coordinates(i,3)-cartesian_coordinates(j, &
                            & 3))/pt1cm_r(m1)
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
                    ygs(i) = ground_state_internal_gradient(i)/pt1cm_r(i)
                    if( first_requested_excited_state .ne. 0 )then
                        do j=first_requested_excited_state, last_requested_excited_state
                            yes(i, j) = excited_state_internal_gradients(i,j)/pt1cm_r(i)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do k=first_requested_coupling, last_requested_coupling
                            yij(i, k) = diabatic_internal_gradients(i,k)/pt1cm_r(i)
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

        !> @brief Evaluate DEDCOU for the OH2_GEN_DIM1981 potential.
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
end module chempotpy_oh2_gen_dim1981_pes

!> @brief Expose OH2_GEN_DIM1981 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh2_gen_dim1981_pes, only: evaluate
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
