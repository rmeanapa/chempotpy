!   System:          OH2  A'' surface
!   Common name:     M2
!   Functional form: Rotated-Morse-oscillator-spline (RMOS) potential
!                    for collinear geometeries plus anti-Morse
!                    bending potential.
!   References:      B. C. Garrett and D. G. Truhlar
!                    Int. J. Quantum Chem. 29, 1463 (1986)
!   Interface:                     potlib2001
!   Number of bodies:              3
!   Number of derivatives:         1
!   Number of electronic surfaces: 1
!   prepare_local_context must be called once before any calls to POT.
!   The control flags for the potential are contained in the block data
!   subprogram PTPACM.
!   Coordinates, potential energy, and derivatives are passed
!   The potential energy in the three asymptotic valleys are
!   stored in the common block ASYCM:
!   The potential energy in the AB valley, EASYAB, is equal to the poten
!   energy of the H "infinitely" far from the OH diatomic, with the
!   OH diatomic at its equilibrium configuration.  Similarly, the terms
!   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
!   respectively.
!   All the information passed through the common blocks PT1CM and ASYCM
!   is in hartree atomic units.
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

!> @brief Implement OH2_App_GEN_M21986, a potential energy surface for OH2.
module chempotpy_oh2_app_gen_m21986_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the OH2_App_GEN_M21986 potential energy surface for OH2.
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
        real(wp)              :: ambpcm_d(3)
        real(wp)              :: ambpcm_de
        real(wp)              :: ambpcm_gamp(3)
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
        real(wp)              :: pre2cm_a1dif
        real(wp)              :: pre2cm_a1inf
        real(wp)              :: pre2cm_a2dif
        real(wp)              :: pre2cm_a2inf
        real(wp)              :: pre2cm_alph1
        real(wp)              :: pre2cm_alph2
        real(wp)              :: pre2cm_b1dif
        real(wp)              :: pre2cm_b1inf
        real(wp)              :: pre2cm_b2dif
        real(wp)              :: pre2cm_b2inf
        real(wp)              :: pre2cm_bsp(20,3)
        real(wp)              :: pre2cm_csp(22,3)
        real(wp)              :: pre2cm_d1dif
        real(wp)              :: pre2cm_d1inf
        real(wp)              :: pre2cm_d2dif
        real(wp)              :: pre2cm_d2inf
        real(wp)              :: pre2cm_dasy
        integer               :: pre2cm_iayes
        integer               :: pre2cm_nphi
        real(wp)              :: pre2cm_phi(20)
        real(wp)              :: pre2cm_r1dif
        real(wp)              :: pre2cm_r1inf
        real(wp)              :: pre2cm_r1s
        real(wp)              :: pre2cm_r2dif
        real(wp)              :: pre2cm_r2inf
        real(wp)              :: pre2cm_r2s
        real(wp)              :: pre2cm_vsp
        real(wp)              :: precm_dbcoor(2)
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
        ambpcm_a                          = 0.e0_wp; ambpcm_beta = 0.e0_wp; ambpcm_d = 0.e0_wp; ambpcm_de = 0.e0_wp; ambpcm_gamp = 0.e0_wp; ambpcm_rbb = 0.e0_wp
        ambpcm_req                        = 0.e0_wp; ambpcm_rhh = 0.e0_wp; asycm_easyab = 0.e0_wp; asycm_easyac = 0.e0_wp; asycm_easybc = 0.e0_wp
        cartesian_coordinates             = 0.e0_wp; coordinate_representation = 0; maximum_derivative_order = 0
        maximum_excited_state             = 0; atom_count = 0
        infocm_ref                        = ' '
        pre2cm_a1dif                      = 0.e0_wp
        pre2cm_a1inf                      = 0.e0_wp
        pre2cm_a2dif                      = 0.e0_wp
        pre2cm_a2inf                      = 0.e0_wp
        pre2cm_alph1                      = 0.e0_wp
        pre2cm_alph2                      = 0.e0_wp
        pre2cm_b1dif                      = 0.e0_wp
        pre2cm_b1inf                      = 0.e0_wp
        pre2cm_b2dif                      = 0.e0_wp
        pre2cm_b2inf                      = 0.e0_wp
        pre2cm_bsp                        = 0.e0_wp
        pre2cm_csp                        = 0.e0_wp
        pre2cm_d1dif                      = 0.e0_wp
        pre2cm_d1inf                      = 0.e0_wp
        pre2cm_d2dif                      = 0.e0_wp
        pre2cm_d2inf                      = 0.e0_wp
        pre2cm_dasy                       = 0.e0_wp
        pre2cm_iayes                      = 0
        pre2cm_nphi                       = 0
        pre2cm_phi                        = 0.e0_wp
        pre2cm_r1dif                      = 0.e0_wp
        pre2cm_r1inf                      = 0.e0_wp
        pre2cm_r1s                        = 0.e0_wp
        pre2cm_r2dif                      = 0.e0_wp
        pre2cm_r2inf                      = 0.e0_wp
        pre2cm_r2s                        = 0.e0_wp
        pre2cm_vsp                        = 0.e0_wp
        precm_dbcoor                      = 0.e0_wp
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
        ambpcm_d                   = [106.56e0_wp, 109.472e0_wp, 106.56e0_wp]
        ambpcm_beta                = 2.07942e0_wp
        ambpcm_req                 = [0.96966e0_wp, 0.74144e0_wp, 0.96966e0_wp]
        ambpcm_a                   = 0.26e0_wp
        ambpcm_gamp                = [0.210820e0_wp, 5.42027e0_wp, -1.44235e0_wp]
        ambpcm_rbb                 = 0.74144e0_wp
        pre2cm_d1inf               = 95.3e0_wp
        pre2cm_r1inf               = .76e0_wp
        pre2cm_b1inf               = 2.0202e0_wp
        pre2cm_d2inf               = 92.4e0_wp
        pre2cm_r2inf               = .98e0_wp
        pre2cm_b2inf               = 2.41105e0_wp
        pre2cm_a1inf               = 0.0e0_wp
        pre2cm_a2inf               = 0.0e0_wp
        pre2cm_r1s                 = 2.21208e0_wp
        pre2cm_r2s                 = 2.21531e0_wp
        pre2cm_vsp                 = 95.3e0_wp
        pre2cm_iayes               = 0
        pre2cm_alph1               = 1.6e0_wp
        pre2cm_alph2               = 1.6e0_wp
        pre2cm_nphi                = 11
        pre2cm_phi(1:11) = [0.0000000e+00_wp, 0.1000000e+02_wp, 0.2000000e+02_wp, 0.3000000e+02_wp, 0.3733000e+02_wp, &
            & 0.4000000e+02_wp, 0.5000000e+02_wp, 0.6000000e+02_wp, 0.7000000e+02_wp, 0.8000000e+02_wp, &
            & 0.9000000e+02_wp]
        pre2cm_bsp(1:11,1) = [0.9477920e+02_wp, 0.9316200e+02_wp, 0.9035320e+02_wp, 0.8517960e+02_wp, &
            & 0.8272720e+02_wp, 0.8304320e+02_wp, 0.8855580e+02_wp, 0.9003100e+02_wp, 0.9129800e+02_wp, &
            & 0.9184900e+02_wp, 0.9209100e+02_wp]
        pre2cm_bsp(1:11,2) = [0.1456500e+01_wp, 0.1475500e+01_wp, 0.1530500e+01_wp, 0.1606400e+01_wp, &
            & 0.1634900e+01_wp, 0.1632300e+01_wp, 0.1534600e+01_wp, 0.1396600e+01_wp, 0.1299800e+01_wp, &
            & 0.1243800e+01_wp, 0.1225600e+01_wp]
        pre2cm_bsp(1:11,3) = [0.2041900e+01_wp, 0.1987500e+01_wp, 0.1905400e+01_wp, 0.1672700e+01_wp, &
            & 0.1554000e+01_wp, 0.1539400e+01_wp, 0.1749300e+01_wp, 0.2014500e+01_wp, 0.2227100e+01_wp, &
            & 0.2348800e+01_wp, 0.2383800e+01_wp]
        pre2cm_csp(1:13,1) = [0.9477920e+02_wp, -0.1478277e+00_wp, -0.1389235e-03_wp, -0.3580591e-03_wp, &
            & 0.1450636e-02_wp, 0.1148717e-02_wp, -0.9458414e-03_wp, -0.4983968e-02_wp, 0.5921762e-02_wp, - &
            & 0.2814975e-02_wp, 0.1001140e-02_wp, - 0.2747850e-03_wp, -0.5702495e-05_wp]
        pre2cm_csp(1:13,2) = [0.1456500e+01_wp, 0.1137130e-02_wp, 0.7628700e-05_wp, -0.9772201e-05_wp, - &
            & 0.1201120e-04_wp, 0.2800151e-05_wp, -0.4114549e-04_wp, 0.7866376e-04_wp, -0.1122307e-04_wp, - &
            & 0.1936703e-04_wp, 0.6791209e-05_wp, - 0.1039780e-04_wp, 0.8032967e-05_wp]
        pre2cm_csp(1:13,3) = [0.2041900e+01_wp, -0.5887558e-02_wp, 0.4475579e-05_wp, -0.5455348e-04_wp, &
            & 0.1230139e-03_wp, - 0.6997501e-04_wp, 0.2822753e-03_wp, -0.4019424e-03_wp, 0.1215636e-03_wp, - &
            & 0.1548133e-04_wp, 0.9961687e-05_wp, 0.1813458e-04_wp, -0.1747243e-04_wp]
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

        !> @brief Evaluate prepare_local_context for the OH2_App_GEN_M21986 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, it
            !   Conversion for Angstroms to Bohr (Angstroms/Bohr)
            real(wp), parameter :: bohr = .52917706e0_wp
            !   Conversion for kcal to Hartree (kcal/Hartree)
            real(wp), parameter :: hartre = hartree_to_kilocalorie_per_mole
            !  Vaule of tolerance
            !   Echo the name of the potential energy surface
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            !   Assign the parameters for the collinear part of the potential
            call prepare_local_context_2
            !   Echo the parameters for the bending potential
            !       CONVERT TO ATOMIC UNITS
            do it = 1, 3
                ambpcm_req(it) = ambpcm_req(it)/bohr
                ambpcm_d(it)   = ambpcm_d(it)/hartre
            enddo
            ambpcm_beta = ambpcm_beta * bohr
            ambpcm_de   = ambpcm_d(1)
            ambpcm_a    = ambpcm_a/bohr
            ambpcm_rhh  = 0.74144e0_wp / bohr
            ambpcm_rbb  = ambpcm_rbb /bohr
            !   Initialize the energy in the asymptotic valleys
            asycm_easyab   = ambpcm_d(1)
            asycm_easybc   = ambpcm_d(2)
            asycm_easyac   = ambpcm_d(3)
            pt3cm_ezero(1) = ambpcm_d(2)
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            infocm_ref(1) = 'B. C. Garrett, D. G. Truhlar'
            infocm_ref(2) = 'Int. J. Quantum Chem. 29, 1463(1986)'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the OH2_App_GEN_M21986 potential.
        subroutine pot
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the OH diatomic, with the
            !   OH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
            !   respectively.
            !        This potential is written such that:
            !                       R(1) = R(O-H)
            !                       R(2) = R(H-H)
            !                       R(3) = R(H-O)
            !   The classical potential energy is set equal to zero for the O
            !   infinitely far from the equilibrium H2 diatomic.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: bcoor, c, dac, dex1, dex2, dgam, ex, ex1, ex2, gamma, raccol, ractmp, rob, stuff, t
            integer  :: it
            real(wp) :: nab, nab1, nbc, nbc1
            !   Conversion for Angstroms to Bohr (Angstroms/Bohr)
            !   Conversion for kcal to Hartree (kcal/Hartree)
            !  Vaule of tolerance
            real(wp), parameter :: tol = 1.e-05_wp
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   Check the values of surface_requests and requested_derivative_order for validity.
            if( pt1cm_r(1) .gt. tol .and. pt1cm_r(2) .gt. tol .and. pt1cm_r(3) .gt. tol )then
                raccol = pt1cm_r(1) + pt1cm_r(2)
                bcoor  = 0.0e0_wp
                if( requested_derivative_order .eq. 1 )then
                    precm_dbcoor(1) = 0.0e0_wp
                    precm_dbcoor(2) = 0.0e0_wp
                    dac             = 0.0e0_wp
                endif
                rob = 1.0e0_wp
                if( abs(ambpcm_rhh-ambpcm_rbb) .gt. 1.e-10_wp )then
                    !     FIRST CALCULATE SOME USEFUL NUMBERS
                    !     NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
                    !     AB AND BC RESPECTIVELY. (EQUATION 2.2)
                    nab = exp((ambpcm_req(1) - pt1cm_r(1))/ambpcm_a)
                    if( nab .lt. 1.e-15_wp ) nab = 1.e-15_wp
                    nbc = exp((ambpcm_req(2) - pt1cm_r(2))/ambpcm_a)
                    if( nbc .lt. 1.e-15_wp ) nbc = 1.e-15_wp
                    !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
                    !       SEE EQUATION 2.7A, 2.7B
                    c = 1.e0_wp - nab/nbc
                    if( abs(c) .lt. 1.e-14_wp )then
                        nab1 = .5e0_wp
                        nbc1 = .5e0_wp
                    else
                        c    = c/nab
                        nab1 = (2.e0_wp + c - sqrt(4.e0_wp+c*c))/(2.e0_wp*c)
                        nbc1 = 1.e0_wp - nab1
                    endif
                    !          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
                    !          FIRST TRAP OUT ANY ZERO ARGUMENTS
                    if( nab1*nbc1 .le. 0.0e0_wp )then
                        rob = 1.e0_wp
                    else
                        stuff = ambpcm_a * log(nab1*nbc1)
                        t     = ambpcm_rhh - stuff
                        rob   = 1.0e0_wp
                        if( abs(t) .gt. 1.e-15_wp ) rob = (ambpcm_rbb - stuff)/t
                    endif
                endif
                !        CALCULATE BENDING CORRECTION
                !        EVALUATE GAMMA AND DERIVATIVE
                ex    = exp(raccol*(ambpcm_gamp(2) + ambpcm_gamp(3)*raccol))
                gamma = ambpcm_gamp(1) + ex
                if( requested_derivative_order .eq. 1 ) dgam = (ambpcm_gamp(2) + 2.e0_wp*ambpcm_gamp(3)*raccol)*ex
                gamma = gamma*ambpcm_de
                if( requested_derivative_order .eq. 1 ) dgam = dgam*ambpcm_de
                !           CALCULATE V(R(3)) TERM
                ex  = exp(-ambpcm_beta*(pt1cm_r(3)-ambpcm_req(3))/rob)
                ex1 = ex*(1.e0_wp + 0.5e0_wp*ex)
                if( requested_derivative_order .eq. 1 ) dex1 = ex*(1.e0_wp + ex)
                !           CALCULATE V(R(1)+R(2)) TERM
                ex  = exp(ambpcm_beta*(ambpcm_req(3)-pt1cm_r(2)-pt1cm_r(1))/rob)
                ex2 = ex*(1.e0_wp + 0.5e0_wp*ex)
                if( requested_derivative_order .eq. 1 ) dex2 = ex*(1.e0_wp + ex)
                !        HERE IS THE BENDING CORRECTION
                bcoor = gamma*(ex1 - ex2)
                !        WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
                !        OF BCOOR WITH RESPECT TO R(1)  R(2) AND R(3).
                !        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
                if( requested_derivative_order .eq. 1 )then
                    gamma = gamma*ambpcm_beta/rob
                    dac   = -gamma*dex1
                    !          CALCULATE CORRECTIONS TO DAB, DBC
                    precm_dbcoor(1) = gamma*dex2 + dgam*(ex1-ex2)
                    precm_dbcoor(2) = precm_dbcoor(1)
                endif
                !       NOW CALCULATE THE COLLINEAR ENERGY
                !       STORE R(3) IN RACTMP THEN SET R(3) TO
                !       THE COLLINEAR GEOMETRY AND CALL POT2
                ractmp     = pt1cm_r(3)
                pt1cm_r(3) = pt1cm_r(2) + pt1cm_r(1)
                call pot2
                pt1cm_r(3)          = ractmp
                ground_state_energy = ground_state_energy + bcoor
                if( requested_derivative_order .eq. 1 )then
                    do it = 1, 2
                        ground_state_internal_gradient(it) = precm_dbcoor(it) + ground_state_internal_gradient(it) + &
                            & ground_state_internal_gradient(3)
                    enddo
                    ground_state_internal_gradient(3) = dac
                endif
            else
                ground_state_energy = 1.e30_wp
                if( requested_derivative_order .eq. 1 )then
                    ground_state_internal_gradient(1)                           = 0.0e0_wp
                    if( pt1cm_r(1) .le. tol ) ground_state_internal_gradient(1) = -ground_state_energy
                    ground_state_internal_gradient(2)                           = 0.0e0_wp
                    if( pt1cm_r(1) .le. tol ) ground_state_internal_gradient(2) = -ground_state_energy
                    ground_state_internal_gradient(3)                           = 0.0e0_wp
                    if( pt1cm_r(1) .le. tol ) ground_state_internal_gradient(3) = -ground_state_energy
                endif
            endif
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate prepare_local_context_2 for the OH2_App_GEN_M21986 potential.
        subroutine prepare_local_context_2
            !  COLLINEAR A + BC POTENTIAL FUNCTION
            !  RMO-SPLINE TYPE FIT
            !  COORDINATES AND ENERGIES PASSED TO POTENTIAL IN ATOMIC UNITS
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: aer, bet, dis, req2, t, t1
            integer             :: i, ii, j, nnn, nphip2
            real(wp), parameter :: detrad = 0.0174532e0_wp
            real(wp), parameter :: ckcal  = hartree_to_kilocalorie_per_mole
            real(wp), parameter :: bohr   = 0.529177e0_wp
            character(len=20)   :: title
            !  TITLE
            title = ' O+H2 MAB 9-12-79 '
            !  ASSIGN COLLINEAR PARAMETERS (UNITS OF KCAL, ANGSTROM, DEGREES)
            !  IN BLOCK DATA SUBROUTINE
            pre2cm_a1inf = 0.0e0_wp
            pre2cm_a2inf = 0.0e0_wp
            pre2cm_iayes = 0
            !  ASSIGN PARAMETERS FOR THE ASYMPTOTIC REACTANT AND PRODUCT REGIONS
            !  IN BLOCK DATA SUBROUTINE
            if( pre2cm_a2inf .eq. 0.0e0_wp ) pre2cm_a2inf = pre2cm_d1inf-pre2cm_d2inf
            !  ASSIGN THE PIVOT POINT, POTENTIAL AT THE PIVOT,
            !  FLAG, AND ASYMPTOTIC EXPONENTIAL PARAMETERS
            !  IN BLOCK DATA SUBROUTINE.  ECHO.
            !     2INF
            nphip2 = pre2cm_nphi + 2
            !  ECHO POTENTIAL PARAMTERS ASSIGNED BY BLOCK DATA SUBROUTINE
            !   30 CONTINUE
            !   40 CONTINUE
            !  CONVERT TO ATOMIC UNITS
            pre2cm_d1inf = pre2cm_d1inf/ckcal
            pre2cm_a1inf = pre2cm_a1inf/ckcal
            pre2cm_r1inf = pre2cm_r1inf/bohr
            pre2cm_b1inf = pre2cm_b1inf*bohr
            pre2cm_d2inf = pre2cm_d2inf/ckcal
            pre2cm_a2inf = pre2cm_a2inf/ckcal
            pre2cm_r2inf = pre2cm_r2inf/bohr
            pre2cm_b2inf = pre2cm_b2inf*bohr
            pre2cm_r1s   = pre2cm_r1s/bohr
            pre2cm_r2s   = pre2cm_r2s/bohr
            pre2cm_vsp   = pre2cm_vsp/ckcal
            pre2cm_alph1 = pre2cm_alph1*bohr
            pre2cm_alph2 = pre2cm_alph2*bohr
            do i = 1, pre2cm_nphi
                pre2cm_bsp(i,1) = pre2cm_bsp(i,1)/ckcal
                pre2cm_bsp(i,2) = pre2cm_bsp(i,2)/bohr
                pre2cm_bsp(i,3) = pre2cm_bsp(i,3)*bohr
            enddo
            nnn = pre2cm_nphi + 2
            do i = 1, nnn
                pre2cm_csp(i,1) = pre2cm_csp(i,1)/ckcal
                pre2cm_csp(i,2) = pre2cm_csp(i,2)/bohr
                pre2cm_csp(i,3) = pre2cm_csp(i,3)*bohr
            enddo
            !  CONVERT FROM DEGREES TO RADIANS
            t = 1.0e0_wp/detrad
            do j = 1, 3
                pre2cm_csp(2,j) = pre2cm_csp(2,j)*t
            enddo
            t = t**3
            do ii = 1, pre2cm_nphi
                pre2cm_phi(ii) = pre2cm_phi(ii)*detrad
                i              = ii + 2
                do j = 1, 3
                    pre2cm_csp(i,j) = pre2cm_csp(i,j)*t
                enddo
            enddo
            !  SET UP CONSTANTS
            !     REACTANT ASYMPTOTE
            dis          = pre2cm_bsp(1,1)
            req2         = pre2cm_bsp(1,2)
            bet          = pre2cm_bsp(1,3)
            t1           = 1.0e0_wp - exp(-bet*req2)
            aer          = pre2cm_vsp - dis*t1*t1
            req2         = pre2cm_r2s - req2
            pre2cm_d1dif = dis - pre2cm_d1inf
            pre2cm_r1dif = req2 - pre2cm_r1inf
            pre2cm_b1dif = bet - pre2cm_b1inf
            pre2cm_a1dif = aer - pre2cm_a1inf
            !     PRODUCT ASYMPTOTE
            dis          = pre2cm_bsp(pre2cm_nphi,1)
            req2         = pre2cm_bsp(pre2cm_nphi,2)
            bet          = pre2cm_bsp(pre2cm_nphi,3)
            t1           = 1.0e0_wp - exp(-bet*req2)
            aer          = pre2cm_vsp - dis*t1*t1
            req2         = pre2cm_r1s - req2
            pre2cm_d2dif = dis - pre2cm_d2inf
            pre2cm_r2dif = req2 - pre2cm_r2inf
            pre2cm_b2dif = bet - pre2cm_b2inf
            pre2cm_a2dif = aer - pre2cm_a2inf
            !  DEFINE ZERO OF ENERGY
            pre2cm_dasy = pre2cm_d1inf
        end subroutine prepare_local_context_2

        !> @brief Evaluate POT2 for the OH2_App_GEN_M21986 potential.
        subroutine pot2
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: aer, angle, bet, daer, dbet, ddis, delr, dis, dphd1, dphd2, dr1, dr2, dreq, drrd1, drrd2, dvdp, r1, r2, r3, req2, rr, scrap(pre2cm_nphi+2), t1, t2, t3
            r1 = pt1cm_r(1)
            r2 = pt1cm_r(2)
            r3 = pt1cm_r(3)
            if( r1 .gt. pre2cm_r1s )then
                if( r2 .gt. pre2cm_r2s )then
                    !   3-ATOM BREAK UP REGION
                    ground_state_energy = pre2cm_dasy
                    if( requested_derivative_order .eq. 1 )then
                        ground_state_internal_gradient(1) = 0.0e0_wp
                        ground_state_internal_gradient(2) = 0.0e0_wp
                    endif
                else
                    !   REACTANT ASYMPTOTIC REGION
                    dr1 = r1 - pre2cm_r1s
                    t1  = exp(-dr1*pre2cm_alph1)
                    t2  = -pre2cm_alph1*t1
                    !.....INTERPOLATE FOR MORSE PARAMETERS
                    dis  = pre2cm_d1inf + pre2cm_d1dif*t1
                    ddis = pre2cm_d1dif*t2
                    req2 = pre2cm_r1inf + pre2cm_r1dif*t1
                    dreq = pre2cm_r1dif*t2
                    bet  = pre2cm_b1inf + pre2cm_b1dif*t1
                    dbet = pre2cm_b1dif*t2
                    if( pre2cm_iayes .eq. 0 )then
                        aer  = pre2cm_vsp - dis
                        daer = -ddis
                    else
                        aer  = pre2cm_a1inf + pre2cm_a1dif*t1
                        daer = pre2cm_a1dif*t2
                    endif
                    delr                = r2 - req2
                    t1                  = exp(-bet*delr)
                    t2                  = 1.0e0_wp - t1
                    ground_state_energy = aer + dis*t2*t2
                    if( requested_derivative_order .eq. 1 )then
                        t3                                = 2.0e0_wp*dis*t1
                        ground_state_internal_gradient(1) = daer + t2*(ddis*t2 - t3*(bet*dreq - dbet*delr))
                        ground_state_internal_gradient(2) = t3*bet*t2
                    endif
                endif
            elseif( r2 .gt. pre2cm_r2s )then
                !   PRODUCT ASYMPTOTIC REGION
                dr2 = r2 - pre2cm_r2s
                !.....INTERPOLATE FOR MORSE PARAMETERS
                t1   = exp(-dr2*pre2cm_alph2)
                t2   = -pre2cm_alph2*t1
                dis  = pre2cm_d2inf + pre2cm_d2dif*t1
                ddis = pre2cm_d2dif*t2
                req2 = pre2cm_r2inf + pre2cm_r2dif*t1
                dreq = pre2cm_r2dif*t2
                bet  = pre2cm_b2inf + pre2cm_b2dif*t1
                dbet = pre2cm_b2dif*t2
                if( pre2cm_iayes .eq. 0 )then
                    aer  = pre2cm_vsp - dis
                    daer = -ddis
                else
                    aer  = pre2cm_a2inf + pre2cm_a2dif*t1
                    daer = pre2cm_a2dif*t2
                endif
                delr                = r1 - req2
                t1                  = exp(-bet*delr)
                t2                  = 1.0e0_wp - t1
                ground_state_energy = aer + dis*t2*t2
                if( requested_derivative_order .eq. 1 )then
                    t3                                = 2.0e0_wp*dis*t1
                    ground_state_internal_gradient(2) = daer + t2*(ddis*t2 - t3*(bet*dreq - dbet*delr))
                    ground_state_internal_gradient(1) = t3*bet*t2
                endif
            else
                !   SPLINE INTERPOLATION REGION
                dr1   = pre2cm_r1s - r1
                dr2   = pre2cm_r2s - r2
                t1    = dr1*dr1 + dr2*dr2
                rr    = sqrt(t1)
                drrd1 = -dr1/rr
                drrd2 = -dr2/rr
                angle = atan(dr1/dr2)
                dphd1 = -dr2/t1
                dphd2 = dr1/t1
                dis   = spline(3, pre2cm_nphi, pre2cm_phi, pre2cm_bsp(1,1), pre2cm_csp(1,1), scrap, angle)
                ddis  = spline(4, pre2cm_nphi, pre2cm_phi, pre2cm_bsp(1,1), pre2cm_csp(1,1), scrap, angle)
                req2  = spline(3, pre2cm_nphi, pre2cm_phi, pre2cm_bsp(1,2), pre2cm_csp(1,2), scrap, angle)
                dreq  = spline(4, pre2cm_nphi, pre2cm_phi, pre2cm_bsp(1,2), pre2cm_csp(1,2), scrap, angle)
                bet   = spline(3, pre2cm_nphi, pre2cm_phi, pre2cm_bsp(1,3), pre2cm_csp(1,3), scrap, angle)
                dbet  = spline(4, pre2cm_nphi, pre2cm_phi, pre2cm_bsp(1,3), pre2cm_csp(1,3), scrap, angle)
                if( pre2cm_iayes .eq. 0 )then
                    aer  = pre2cm_vsp - dis
                    daer = -ddis
                else
                    t1   = exp(-bet*req2)
                    t2   = 1.0e0_wp - t1
                    aer  = pre2cm_vsp - dis*t2*t2
                    daer = -t2*(ddis*t2 + 2.0e0_wp*dis*t1*(bet*dreq + req2*dbet))
                endif
                delr                = rr - req2
                t1                  = exp(bet*delr)
                t2                  = 1.0e0_wp - t1
                ground_state_energy = aer + dis*t2*t2
                if( requested_derivative_order .eq. 1 )then
                    t3                                = 2.0e0_wp*dis*t1
                    dvdp                              = daer + t2*(ddis*t2 - t3*(dbet*delr - bet*dreq))
                    t3                                = t3*t2*bet
                    ground_state_internal_gradient(1) = dvdp*dphd1 - t3*drrd1
                    ground_state_internal_gradient(2) = dvdp*dphd2 - t3*drrd2
                endif
            endif
            if( requested_derivative_order .eq. 1 )then
                ground_state_internal_gradient(1) = ground_state_internal_gradient(1)
                ground_state_internal_gradient(2) = ground_state_internal_gradient(2)
                ground_state_internal_gradient(3) = 0.0e0_wp
            endif
        end subroutine pot2

        !> @brief Evaluate SPLINE for the OH2_App_GEN_M21986 potential.
        function spline( isw, nn, x, y, c, d, xpt )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: c(nn+2), cfact, d(nn+2), deriv, dfact, pivot, spline, supd, v, v2, x(nn), xpt, y(nn), z
            integer  :: i, isw, j, n, nn, np1, np2
            !  THIS IS A SUBROUTINE FOR FITTING DATA WITH A CUBIC SPLINE
            !  POLYNOMIAL AND EVALUATING THAT POLYNOMIAL AT A GIVEN POINT
            !  OR ITS DERIVATIVE AT A GIVEN POINT
            !  CALLING SEQUENCE .......
            !     ISW ... CONTROL OPTION
            !         ISW=1  IF A CUBIC SPLINE IS TO BE FITTED TO THE SET OF KNOTS
            !                DEFINED BY THE ARRAYS X AND Y.  THE SPLINE COEFFICIENTS
            !                ARE STORED IN THE ARRAY C.
            !         ISW=2  IF THE SPLINE DEFINED BY THE COEFFICIENT ARRAY 'C' IS
            !                TO BE EVALUATED (INTERPOLATED) AT THE POINT DEFINED BY
            !                THE PARAMETER 'XPT'.
            !         ISW=3  AS IN ISW=2, ONLY THE DERIVATIVE/3.D0 IS ALSO CALCULATE
            !         ISW=4  THE DERIVATIVE CALCULATED BY THE LAST USE OF SPLINE WIT
            !                IS RETURNED.
            !     NN ... THE NUMBER OF KNOTS (DATA POINTS) TO WHICH THE SPLINE IS TO
            !            BE FITTED
            !     X,Y ... THE ARRAYS DEFINING THE KNOTS.  THE X-VALUES MUST BE IN
            !             INCREASING ORDER.  THE ARRAYS MUST BE DIMENSIONED AT LEAST
            !             NN.
            !     C ... THE ARRAY THAT CONTAINS THE CUBIC SPLINE COEFFICIENTS.
            !           MUST BE DIMENSIONED AT LEAST NN+2 .
            !     D ... A WORK SPACE.  MUST BE DIMENSIONED AT LEAST NN+2 .
            !     XPT ... THE POINT AT WHICH THE INTERPOLATION IS DESIRED (IF ISW IS
            !              SET TO 2).  THE VALUE OF SPLINE IS SET TO THE
            !              INTERPOLATED VALUE.
            ! *****  USER NOTES  *****
            !     INTERPOLATION INVOLVES AT LEAST TWO STEPS .......
            !       A.  CALL SPLINE WITH THE KNOTS.  THIS SETS UP THE
            !           COEFFICIENT ARRAY C.
            !           EG.  DUMY=SPLINE(1,NN,X,Y,C,D,XPT)
            !       B.  CALL SPLINE WITH THE ARRAY C WHICH WAS DEFINED BY THE
            !           PREVIOUS CALL AND WILL BE USED TO FIND THE VALUE AT THE
            !           POINT 'XPT' .
            !           EG.   VALUE=SPLINE(2,NN,X,Y,C,D,XPT)
            !     STEP 'A' NEED BE EXECUTED ONLY ONCE FOR A GIVEN SET OF KNOTS.
            !     STEP B MAY BE EXECUTED AS MANY TIMES AS NECESSARY.
            !     Output from this file is written to UNIT NFLAG(18) which is
            !     defined in a parameter statement.
            n   = nn
            np1 = n+1
            np2 = n+2
            z   = xpt
            if( isw .eq. 1 )then
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
                    if( d(i-1) .eq. 0.0e0_wp )then
                        if( present(status) ) call status%set(status_invalid_input, &
                            'SPLINE 1')
                        return
                    endif
                    pivot = 1.0e0_wp/d(i-1)
                    if( i .lt. np2 )then
                        supd = x(i-1)-x(i-2)
                        if( supd .lt. 0.0e0_wp )then
                            if( present(status) ) call status%set(status_invalid_input, &
                                'SPLINE 2')
                            return
                        endif
                        supd = supd*supd*supd
                    else
                        supd = 1.0e0_wp
                    endif
                    dfact = supd*pivot
                    cfact = c(i-1)*pivot
                    if( i .le. n )then
                        do j = i, n
                            v    = x(j)-x(i-2)
                            c(j) = c(j)-d(j)*cfact
                            d(j) = v*v*v-d(j)*dfact
                        enddo
                    endif
                    if( i .lt. np2 )then
                        c(np1) = c(np1)-d(np1)*cfact
                        d(np1) = 1.0e0_wp-d(np1)*dfact
                    endif
                    c(np2) = c(np2)-d(np2)*cfact
                    d(np2) = x(i-2)-d(np2)*dfact
                enddo
                do i = 1, n
                    j = np2-i
                    if( j .eq. np1 )then
                        v = 1.0e0_wp
                    else
                        v = x(j)-x(j-1)
                        v = v*v*v
                    endif
                    if( d(j+1) .eq. 0.0e0_wp )then
                        if( present(status) ) call status%set(status_invalid_input, &
                            'SPLINE 3')
                        return
                    endif
                    c(j+1) = c(j+1)/d(j+1)
                    c(j)   = c(j)-c(j+1)*v
                enddo
                if( d(2) .eq. 0.0e0_wp )then
                    if( present(status) ) call status%set(status_invalid_input, &
                        'SPLINE 4')
                    return
                endif
                c(2) = c(2)/d(2)
            elseif( isw .eq. 2 )then
                spline = c(1) + c(2)*(z-x(1))
                do i = 1, n
                    v = z-x(i)
                    if( v .le. 0 ) return
                    spline = spline + c(i+2)*v**3
                enddo
            elseif( isw .eq. 3 )then
                spline = c(1) + c(2)*(z-x(1))
                deriv  = c(2)/3.e0_wp
                do i = 1, n
                    v = z-x(i)
                    if( v .le. 0 ) return
                    v2     = v*v
                    spline = spline + c(i+2)*v2*v
                    deriv  = deriv + c(i+2)*v2
                enddo
            elseif( isw .eq. 4 )then
                spline = 3.e0_wp*deriv
            endif
        end function spline

        !> @brief Evaluate POTINFO for the OH2_App_GEN_M21986 potential.
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

        !> @brief Evaluate CARTTOR for the OH2_App_GEN_M21986 potential.
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

        !> @brief Evaluate EUNITZERO for the OH2_App_GEN_M21986 potential.
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

        !> @brief Evaluate RTOCART for the OH2_App_GEN_M21986 potential.
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

        !> @brief Evaluate DEDCOU for the OH2_App_GEN_M21986 potential.
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
end module chempotpy_oh2_app_gen_m21986_pes

!> @brief Expose OH2_App_GEN_M21986 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh2_app_gen_m21986_pes, only: evaluate
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
