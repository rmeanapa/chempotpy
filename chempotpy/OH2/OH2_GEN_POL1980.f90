!   System:          OH2
!   Functional form: RMO-Spline fit
!   Common name:     POL
!   Interface:       potlib2001
!   Reference:
!   Number of bodies:              3
!   Number of derivatives:         1
!   Number of electronic surfaces: 1
!   prepare_local_context must be called once before any calls to POT.
!   The potential parameters are assigned in the Block Data Subprogram
!   PTPACM.  The control flags for the potential also are initialized
!   in the Block Data Subprogram PPTARM.
!   Coordinates, potential energy, and derivatives are passed
!   All the information passed through the common block PT1CM
!   is in Hartree atomic units.
!        This potential is written such that:
!                       R(1) = R(O-H)
!                       R(2) = R(H-H)
!                       R(3) = R(H-O)
!   The classical potential energy is set equal to zero for the O
!   infinitely far from the equilibrium H2 diatomic.
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
!          ANTI-MORSE BENDING POTENTIAL
!          SUBROUTINE BEND WILL ADD A BENDING CORRECTION TO ANY
!          COLLINEAR POTENTIAL SUBROUTINE.
!          PREBEND MUST BE CALLED BEFORE ANY CALLS TO BEND
!          CAN BE MADE.  THE POTENTIAL PARAMETERS FOR THE COL-
!          LINEAR POTENTIALS ARE PROVIDED BY prepare_local_context
!          (CONTAINED IN PREBEND).  THEN THE BENDING PARAMETERS ARE
!          COORDINATES, POTENTIAL ENERGY, AND DERIVATIVES ARE
!          PASSED THROUGH COMMON /PT1CM/ R(3), ground_state_energy, ground_state_internal_gradient(3).
!CCCCCCCCCCCCCCCCCCCCCCCCCC
!          DEFINE VARIABLES
!CCCCCCCCCCCCCCCCCCCCCCCCCC

!> @brief Implement OH2_GEN_POL1980, a potential energy surface for OH2.
module chempotpy_oh2_gen_pol1980_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the OH2_GEN_POL1980 potential energy surface for OH2.
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
        real(wp)              :: ambpcm_a
        real(wp)              :: ambpcm_beta
        real(wp)              :: ambpcm_de
        real(wp)              :: ambpcm_gamma
        real(wp)              :: ambpcm_rbb
        real(wp)              :: ambpcm_req(3)
        real(wp)              :: ambpcm_rhh
        real(wp)              :: asycm_easyab
        real(wp)              :: asycm_easyac
        real(wp)              :: asycm_easybc
        real(wp)              :: cartesian_coordinates(natom,3)
        integer               :: coordinate_representation
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: atom_count
        character(len=75)     :: infocm_ref(5)
        real(wp)              :: pre2cm_a1
        real(wp)              :: pre2cm_a1dif
        real(wp)              :: pre2cm_a1inf
        real(wp)              :: pre2cm_a2
        real(wp)              :: pre2cm_a2dif
        real(wp)              :: pre2cm_a2inf
        real(wp)              :: pre2cm_alph1
        real(wp)              :: pre2cm_alph2
        real(wp)              :: pre2cm_b1
        real(wp)              :: pre2cm_b1dif
        real(wp)              :: pre2cm_b1inf
        real(wp)              :: pre2cm_b2
        real(wp)              :: pre2cm_b2dif
        real(wp)              :: pre2cm_b2inf
        real(wp)              :: pre2cm_bsp(20,5)
        real(wp)              :: pre2cm_csp(22,5)
        real(wp)              :: pre2cm_d1
        real(wp)              :: pre2cm_d1dif
        real(wp)              :: pre2cm_d1inf
        real(wp)              :: pre2cm_d2
        real(wp)              :: pre2cm_d2dif
        real(wp)              :: pre2cm_d2inf
        integer               :: pre2cm_iayes
        integer               :: pre2cm_nphi
        integer               :: pre2cm_nv
        real(wp)              :: pre2cm_phi(20)
        real(wp)              :: pre2cm_r1dif
        real(wp)              :: pre2cm_r1e
        real(wp)              :: pre2cm_r1einf
        real(wp)              :: pre2cm_r1s
        real(wp)              :: pre2cm_r2dif
        real(wp)              :: pre2cm_r2e
        real(wp)              :: pre2cm_r2inf
        real(wp)              :: pre2cm_r2s
        real(wp)              :: pre2cm_step
        real(wp)              :: pre2cm_v(4)
        real(wp)              :: pre2cm_vsp
        real(wp)              :: precm_bcoor
        real(wp)              :: precm_dbcoor(2)
        real(wp)              :: precm_junk(4)
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
        ambpcm_a                          = 0.e0_wp; ambpcm_beta = 0.e0_wp; ambpcm_de = 0.e0_wp; ambpcm_gamma = 0.e0_wp; ambpcm_rbb = 0.e0_wp; ambpcm_req = 0.e0_wp
        ambpcm_rhh                        = 0.e0_wp; asycm_easyab = 0.e0_wp; asycm_easyac = 0.e0_wp; asycm_easybc = 0.e0_wp; cartesian_coordinates = 0.e0_wp
        coordinate_representation         = 0; maximum_derivative_order = 0; maximum_excited_state = 0; atom_count = 0
        infocm_ref                        = ' '
        pre2cm_a1                         = 0.e0_wp
        pre2cm_a1dif                      = 0.e0_wp
        pre2cm_a1inf                      = 0.e0_wp
        pre2cm_a2                         = 0.e0_wp
        pre2cm_a2dif                      = 0.e0_wp
        pre2cm_a2inf                      = 0.e0_wp
        pre2cm_alph1                      = 0.e0_wp
        pre2cm_alph2                      = 0.e0_wp
        pre2cm_b1                         = 0.e0_wp
        pre2cm_b1dif                      = 0.e0_wp
        pre2cm_b1inf                      = 0.e0_wp
        pre2cm_b2                         = 0.e0_wp
        pre2cm_b2dif                      = 0.e0_wp
        pre2cm_b2inf                      = 0.e0_wp
        pre2cm_bsp                        = 0.e0_wp
        pre2cm_csp                        = 0.e0_wp
        pre2cm_d1                         = 0.e0_wp
        pre2cm_d1dif                      = 0.e0_wp
        pre2cm_d1inf                      = 0.e0_wp
        pre2cm_d2                         = 0.e0_wp
        pre2cm_d2dif                      = 0.e0_wp
        pre2cm_d2inf                      = 0.e0_wp
        pre2cm_iayes                      = 0
        pre2cm_nphi                       = 0
        pre2cm_nv                         = 0
        pre2cm_phi                        = 0.e0_wp
        pre2cm_r1dif                      = 0.e0_wp
        pre2cm_r1e                        = 0.e0_wp
        pre2cm_r1einf                     = 0.e0_wp
        pre2cm_r1s                        = 0.e0_wp
        pre2cm_r2dif                      = 0.e0_wp
        pre2cm_r2e                        = 0.e0_wp
        pre2cm_r2inf                      = 0.e0_wp
        pre2cm_r2s                        = 0.e0_wp
        pre2cm_step                       = 0.e0_wp
        pre2cm_v                          = 0.e0_wp
        pre2cm_vsp                        = 0.e0_wp
        precm_bcoor                       = 0.e0_wp
        precm_dbcoor                      = 0.e0_wp
        precm_junk                        = 0.e0_wp
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
        ambpcm_de                  = 106.56e0_wp
        ambpcm_beta                = 2.07942e0_wp
        ambpcm_req                 = [0.96966e0_wp, 0.74144e0_wp, 0.96966e0_wp]
        ambpcm_a                   = 0.26e0_wp
        ambpcm_gamma               = 0.4131e0_wp
        ambpcm_rbb                 = 0.74144e0_wp
        pre2cm_d1inf               = 97.23777e0_wp
        pre2cm_r1einf              = .75933e0_wp
        pre2cm_b1inf               = 2.01627e0_wp
        pre2cm_d2inf               = 98.81423e0_wp
        pre2cm_r2inf               = .98712e0_wp
        pre2cm_b2inf               = 2.30153e0_wp
        pre2cm_a1inf               = 0.0e0_wp
        pre2cm_a2inf               = -0.1e0_wp
        pre2cm_r1s                 = 2.21208e0_wp
        pre2cm_r2s                 = 2.21531e0_wp
        pre2cm_vsp                 = 87.65999e0_wp
        pre2cm_iayes               = 1
        pre2cm_alph1               = 1.0e0_wp
        pre2cm_alph2               = 1.0e0_wp
        pre2cm_step                = 0.000000010e0_wp
        pre2cm_nphi                = 11
        pre2cm_nv                  = 3
        pre2cm_phi(1:11) = [0.0000000e+00_wp, 0.1000000e+02_wp, 0.2000000e+02_wp, 0.3000000e+02_wp, 0.3732999e+02_wp, &
            & 0.4000000e+02_wp, 0.5000000e+02_wp, 0.6000000e+02_wp, 0.7000000e+02_wp, 0.8000000e+02_wp, &
            & 0.9000000e+02_wp]
        pre2cm_bsp(1:11,1) = [0.9723777e+02_wp, 0.9597664e+02_wp, 0.9292151e+02_wp, 0.9047624e+02_wp, &
            & 0.9036362e+02_wp, 0.9132085e+02_wp, 0.9458319e+02_wp, 0.9728180e+02_wp, 0.9809436e+02_wp, &
            & 0.9854700e+02_wp, 0.9881423e+02_wp]
        pre2cm_bsp(1:11,2) = [0.1455982e+01_wp, 0.1474749e+01_wp, 0.1529910e+01_wp, 0.1605620e+01_wp, &
            & 0.1633813e+01_wp, 0.1631247e+01_wp, 0.1533492e+01_wp, 0.1395812e+01_wp, 0.1298821e+01_wp, &
            & 0.1242818e+01_wp, 0.1224964e+01_wp]
        pre2cm_bsp(1:11,3) = [0.2016270e+01_wp, 0.1955010e+01_wp, 0.1879873e+01_wp, 0.1619415e+01_wp, &
            & 0.1484336e+01_wp, 0.1467517e+01_wp, 0.1689022e+01_wp, 0.1935001e+01_wp, 0.2141747e+01_wp, &
            & 0.2261585e+01_wp, 0.2301529e+01_wp]
        pre2cm_csp(1:13,1) = [0.9723777e+02_wp, -0.7558417e-01_wp, -0.5052923e-03_wp, 0.1237764e-02_wp, - &
            & 0.7532287e-03_wp, 0.2249910e-02_wp, -0.1043168e-01_wp, 0.8624236e-02_wp, -0.9957934e-03_wp, &
            & 0.1126430e-02_wp, -0.6615264e-03_wp, 0.1681010e-03_wp, -0.5891778e-04_wp]
        pre2cm_csp(1:13,2) = [0.1455982e+01_wp, 0.1101632e-02_wp, 0.7749966e-05_wp, -0.1010505e-04_wp, - &
            & 0.1182017e-04_wp, 0.3430826e-05_wp, -0.4369002e-04_wp, 0.8096568e-04_wp, -0.1197233e-04_wp, - &
            & 0.1871181e-04_wp, 0.6506310e-05_wp, - 0.1045387e-04_wp, 0.8100479e-05_wp]
        pre2cm_csp(1:13,3) = [0.2016270e+01_wp, -0.7200198e-02_wp, 0.1074172e-04_wp, -0.7832711e-04_wp, &
            & 0.1557412e-03_wp, - 0.8529179e-04_wp, 0.3140934e-03_wp, -0.4585322e-03_wp, 0.1656866e-03_wp, - &
            & 0.4269112e-04_wp, 0.2110969e-04_wp, 0.1294235e-04_wp, -0.1547269e-04_wp]
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

        !> @brief Evaluate prepare_local_context for the OH2_GEN_POL1980 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, it
            !   Conversion for Angstroms to Bohr (Angstrom/Bohr)
            real(wp), parameter :: bohr = .52917706e0_wp
            !   Conversion for kcal to Hartree (kcal/Hartree)
            real(wp), parameter :: hrtree = hartree_to_kilocalorie_per_mole
            !   Echo the name of the potential routine to UNIT NFLAG(18)
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !      ASSIGN COLLINEAR POTENTIAL PARAMETERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            call prepare_local_context_2
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !      ASSIGN BENDING POTENTIAL PARAMETERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !       CONVERT TO ATOMIC UNITS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            do it = 1, 3
                ambpcm_req(it) = ambpcm_req(it)/bohr
            enddo
            ambpcm_beta = ambpcm_beta * bohr
            ambpcm_de   = ambpcm_de/hrtree
            ambpcm_a    = ambpcm_a/bohr
            ambpcm_rhh  = 1.4007977e0_wp
            ambpcm_rbb  = ambpcm_rbb /bohr
            !    Echo the potential parameters to the file linked to UNIT NFLAG(18)
            !CCCCCCCCCCC
            pt3cm_ezero(1) = 0.0e0_wp
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            infocm_ref(1) = 'No Reference'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the OH2_GEN_POL1980 potential.
        subroutine pot
            !        This potential is written such that:
            !                       R(1) = R(O-H)
            !                       R(2) = R(H-H)
            !                       R(3) = R(H-O)
            !   The classical potential energy is set equal to zero for the O
            !   infinitely far from the equilibrium H2 diatomic.
            !CCCCCCCCCCCCCCC
            !CCCCCCCCCCCCCCC
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: c, dac, ractmp, rob, stuff, vabbc, vac
            integer  :: it
            real(wp) :: nab, nab1, nbc, nbc1
            !   Conversion for Angstroms to Bohr (Angstrom/Bohr)
            !   Conversion for kcal to Hartree (kcal/Hartree)
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   Check the values of surface_requests and requested_derivative_order for validity.
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        FIRST CALCULATE SOME USEFUL NUMBERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
            !        AB AND BC RESPECTIVELY. (EQUATION 2.2)
            nab = exp((ambpcm_req(1) - pt1cm_r(1))/ambpcm_a)
            if( nab .lt. 1.e-15_wp ) nab = 1.e-15_wp
            nbc = exp((ambpcm_req(2) - pt1cm_r(2))/ambpcm_a)
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
                stuff = ambpcm_a * log(nab1*nbc1)
                rob   = (ambpcm_rbb - stuff)/(ambpcm_rhh - stuff)
            else
                rob = 1.0e0_wp
            endif
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        CALCULATE BENDING CORRECTION
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        CALCULATE V(R(3))
            precm_junk(3) = exp(-ambpcm_beta*(pt1cm_r(3)-ambpcm_req(3))/rob)
            vac           = ambpcm_gamma * ambpcm_de * precm_junk(3)* (1.e0_wp + 0.5e0_wp*precm_junk(3))
            !        CALCULATE V(R(1)+R(2))
            precm_junk(4) = exp(ambpcm_beta*(ambpcm_req(3)-pt1cm_r(2)-pt1cm_r(1))/rob)
            vabbc         = ambpcm_gamma * ambpcm_de * precm_junk(4) * (1.e0_wp + 0.5e0_wp*precm_junk(4))
            !        HERE IS THE BENDING CORRECTION
            precm_bcoor = vac - vabbc
            !          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
            !          OF BCOOR WITH RESPECT TO R(1)  R(2) AND R(3).
            !        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
            if( requested_derivative_order .eq. 1 )then
                dac             = -ambpcm_gamma*ambpcm_de*ambpcm_beta*precm_junk(3)*(1.e0_wp+precm_junk(3))/rob
                precm_dbcoor(1) = ambpcm_gamma * ambpcm_de * ambpcm_beta * precm_junk(4) * (1.e0_wp + precm_junk(4))/rob
                precm_dbcoor(2) = precm_dbcoor(1)
            endif
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        NOW CALCULATE THE COLLINEAR ENERGY
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR
            !          GEOMETRY AND CALL POT
            ractmp     = pt1cm_r(3)
            pt1cm_r(3) = pt1cm_r(2) + pt1cm_r(1)
            call pot2
            !   Add the energy terms together
            ground_state_energy = ground_state_energy + precm_bcoor
            !   If requested_derivative_order = 1, compute the derivatives
            if( requested_derivative_order .ne. 1 )then
                call eunitzero
                if( requested_derivative_order .ne. 0 )then
                    call rtocart
                    call dedcou
                endif
                return
            endif
            ground_state_internal_gradient(1) = ground_state_internal_gradient(1) + ground_state_internal_gradient(3)
            ground_state_internal_gradient(2) = ground_state_internal_gradient(2) + ground_state_internal_gradient(3)
            ground_state_internal_gradient(3) = dac
            pt1cm_r(3)                        = ractmp
            do it = 1, 2
                ground_state_internal_gradient(it) = precm_dbcoor(it) + ground_state_internal_gradient(it)
            enddo
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate prepare_local_context_2 for the OH2_GEN_POL1980 potential.
        subroutine prepare_local_context_2
            !   POTENTIAL FUNCTION FOR RATEID SERIES
            !   RMO-SPLINE TYPE FIT
            !   POTENTIAL ROUTINE INTERACTS WITH REST OF PROGRAM
            !   IN UNITS OF HARTREE/BOHR
            !   POTENTIAL ROUTINE ASSUMES SPLINE INFO IN UNITS OF KCAL/ANGSTROM
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)          :: betpr, betrg, depr, derg, repr, rerg
            integer           :: nphip2
            character(len=20) :: title
            !  TITLE
            title = '   O+H2 AB 9-12-79 '
            !  ASSIGN COLLINEAR PARAMETERS (UNITS OF KCAL, ANGSTROM, DEGREES)
            !  IN BLOCK DATA SUBROUTINE
            !  ASSIGN PARAMETERS FOR THE ASYMPTOTIC REACTANT AND PRODUCT REGIONS
            !  IN BLOCK DATA SUBROUTINE
            if( pre2cm_a2inf .eq. 0.0e0_wp ) pre2cm_a2inf= pre2cm_d1inf-pre2cm_d2inf
            !  ASSIGN THE PIVOT POINT, POTENTIAL AT THE PIVOT,
            !  FLAG, AND ASYMPTOTIC EXPONENTIAL PARAMETERS
            !  IN BLOCK DATA SUBROUTINE.  ECHO.
            !     NF
            nphip2 = pre2cm_nphi + 2
            !  ECHO POTENTIAL PARAMTERS ASSIGNED BY BLOCK DATA SUBROUTINE
            !  20  CONTINUE
            !  22  CONTINUE
            pre2cm_d1  = pre2cm_bsp(1,1)
            pre2cm_d2  = pre2cm_bsp(pre2cm_nphi,1)
            pre2cm_r1e = pre2cm_bsp(1,2)
            pre2cm_r2e = pre2cm_bsp(pre2cm_nphi,2)
            pre2cm_b1  = pre2cm_bsp(1,3)
            pre2cm_b2  = pre2cm_bsp(pre2cm_nphi,3)
            pre2cm_a1  = pre2cm_vsp-pre2cm_d1*(1.0e0_wp-exp(-pre2cm_b1*pre2cm_r1e))**2.0e0_wp
            pre2cm_a2  = pre2cm_vsp-pre2cm_d2*(1.0e0_wp-exp(-pre2cm_b2*pre2cm_r2e))**2.0e0_wp
            pre2cm_r1e = pre2cm_r2s -pre2cm_r1e
            pre2cm_r2e = pre2cm_r1s -pre2cm_r2e
            depr       = pre2cm_d2inf/23.061e0_wp
            betpr      = pre2cm_b2inf*0.529177e0_wp
            repr       = pre2cm_r2inf/0.529177e0_wp
            derg       = pre2cm_d1inf/23.061e0_wp
            betrg      = pre2cm_b1inf*0.529177e0_wp
            rerg       = pre2cm_r1einf/0.529177e0_wp
        end subroutine prepare_local_context_2

        !> @brief Evaluate POT2 for the OH2_GEN_POL1980 potential.
        subroutine pot2
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: angle, anu, ber, bet, da, dbe, delr, der, dis, dr, dr1, dr2, dre, factr, potsp, r1, r2, &
                & r3, req2, rer, rr, scrap(20)
            integer             :: ic
            real(wp), parameter :: detrad = 0.0174532e0_wp
            real(wp), parameter :: ckcal  = hartree_to_kilocalorie_per_mole
            real(wp), parameter :: bohr   = 0.529177e0_wp
            do ic = 1, 3
                select case( ic )
                case( 1 )
                    r1 = pt1cm_r(1)
                    r2 = pt1cm_r(2)
                    r3 = pt1cm_r(3)
                case( 2 )
                    r1 = pt1cm_r(1)+pre2cm_step
                    r2 = pt1cm_r(2)
                    r3 = pt1cm_r(3)
                case( 3 )
                    r1 = pt1cm_r(1)
                    r2 = pt1cm_r(2)+pre2cm_step
                    r3 = pt1cm_r(3)
                end select
                r1 = r1*bohr
                r2 = r2*bohr
                r3 = r3*bohr
                if( r1 > pre2cm_r1s .and. r2 >= pre2cm_r2s )then
                    ground_state_energy               = pre2cm_vsp/ckcal
                    ground_state_internal_gradient(1) = 1.0e-9_wp
                    ground_state_internal_gradient(2) = 1.0e-9_wp
                    ground_state_internal_gradient(3) = 0.0e0_wp
                    return
                elseif( r1 > pre2cm_r1s )then
                    dr    = pre2cm_d1-pre2cm_d1inf
                    dre   = pre2cm_r1e-pre2cm_r1einf
                    dbe   = pre2cm_b1-pre2cm_b1inf
                    da    = pre2cm_a1-pre2cm_a1inf
                    delr  = r1-pre2cm_r1s
                    factr = exp(-delr*pre2cm_alph1)
                    der   = pre2cm_d1inf+dr*factr
                    rer   = pre2cm_r1einf+dre*factr
                    ber   = pre2cm_b1inf+dbe*factr
                    if( pre2cm_iayes == 0 ) anu = pre2cm_vsp-der
                    if( pre2cm_iayes /= 0 ) anu = pre2cm_a1inf+da*factr
                    dr2   = r2-rer
                    potsp = der*(1.0e0_wp-exp(-ber*dr2))**2+anu
                elseif( r2 > pre2cm_r2s )then
                    dr    = pre2cm_d2-pre2cm_d2inf
                    dre   = pre2cm_r2e-pre2cm_r2inf
                    dbe   = pre2cm_b2-pre2cm_b2inf
                    da    = pre2cm_a2-pre2cm_a2inf
                    delr  = r2-pre2cm_r2s
                    factr = exp(-delr*pre2cm_alph2)
                    der   = pre2cm_d2inf+dr*factr
                    rer   = pre2cm_r2inf+dre*factr
                    ber   = pre2cm_b2inf+dbe*factr
                    if( pre2cm_iayes == 0 ) anu = pre2cm_vsp-der
                    if( pre2cm_iayes /= 0 ) anu = pre2cm_a2inf+da*factr
                    dr1   = r1-rer
                    potsp = der*(1.0e0_wp-exp(-ber*dr1))**2+anu
                else
                    dr1   = pre2cm_r1s-r1
                    dr2   = pre2cm_r2s-r2
                    angle = atan(dr1/dr2)/detrad
                    rr    = sqrt(dr1*dr1+dr2*dr2)
                    dis   = spline(2, pre2cm_nphi, pre2cm_phi, pre2cm_bsp(1,1), pre2cm_csp(1,1), scrap, angle)
                    req2  = spline(2, pre2cm_nphi, pre2cm_phi, pre2cm_bsp(1,2), pre2cm_csp(1,2), scrap, angle)
                    bet   = spline(2, pre2cm_nphi, pre2cm_phi, pre2cm_bsp(1,3), pre2cm_csp(1,3), scrap, angle)
                    if( pre2cm_iayes == 0 ) anu = pre2cm_vsp-dis
                    if( pre2cm_iayes /= 0 ) anu = pre2cm_vsp-dis*(1.0e0_wp-exp(-bet*req2))**2
                    potsp = dis*(1.0e0_wp-exp(bet*(rr-req2)))**2+anu
                endif
                pre2cm_v(ic) = potsp
            enddo
            pre2cm_v(1:3)       = pre2cm_v(1:3)/ckcal
            ground_state_energy = pre2cm_v(1)
            if( requested_derivative_order == 1 )then
                ground_state_internal_gradient(1) = (pre2cm_v(2)-pre2cm_v(1))/pre2cm_step
                ground_state_internal_gradient(2) = (pre2cm_v(3)-pre2cm_v(1))/pre2cm_step
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

        !> @brief Evaluate POTINFO for the OH2_GEN_POL1980 potential.
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

        !> @brief Evaluate CARTTOR for the OH2_GEN_POL1980 potential.
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

        !> @brief Evaluate EUNITZERO for the OH2_GEN_POL1980 potential.
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

        !> @brief Evaluate RTOCART for the OH2_GEN_POL1980 potential.
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

        !> @brief Evaluate DEDCOU for the OH2_GEN_POL1980 potential.
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
end module chempotpy_oh2_gen_pol1980_pes

!> @brief Expose OH2_GEN_POL1980 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh2_gen_pol1980_pes, only: evaluate
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
