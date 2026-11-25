!   System:          FH2
!   Functional form: DMBE (double-many-body expansion)
!   Common name:     5SEC
!   References:      G. C. Lynch, R. Steckler, D.W. Schwenke, A.J.C. Var
!                    D.G. Truhlar, and B.C. Garrett,
!                    J. Chem. Phys. 94, 7136-7149 (1991)
! Number of bodies: 3
! Interface: potlib2001
! Number of electronic states: 1
! Number of derivatives: 1
!   Protocol:
!      prepare_local_context - initializes the potential's variables
!               must be called once before any calls to POT
!      POT    - driver for the evaluation of the energy and the derivati
!               of the energy with respect to the coordinates for a give
!               geometric configuration
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      ground electronic state
!   Zero of energy:
!      The classical potential energy is set equal to zero for the F
!      infinitely far from the H2 diatomic and R(H2) set equal to the
!      H2 equilibrium diatomic value.
!   Parameters:
!      Set in the BLOCK DATA subprogram PTPACM
!   Coordinates:
!      Internal, Definition: R(1) = R(F-H)
!                            R(2) = R(H-H)
!                            R(3) = R(H-F)
!   Common Blocks (used between the calling program and this potential):
!        passes the coordinates, energy, and derivatives of
!        the energy with respect to the coordinates.
!        passes the control flags where
!        requested_derivative_order  = 0 => no derivatives are calculated
!        requested_derivative_order  = 1 => calculate first derivatives of the energy for the
!                     ground electronic state with respect to the coordi
!        NFLAG  - Control flags
!      NFLAG(18-20)
!        passes the FORTRAN unit number used for potential output
!      /ASYCM/ EASYAB, EASYBC, EASYAC
!        passes the energy in the three asymptotic valleys for an A+BC s
!        The energy in the AB valley, EASYAB, is equal to the energy
!        of the C atom "infinitely" far from the AB diatomic and R(AB) s
!        equal to Re(AB), the equilibrium bond length for the AB diatomi
!        In this potential the AB valley represents H infinitely far fro
!        the FH diatomic and R(FH) equal to Re(FH).
!        Similarly, the terms EASYBC and EASYAC represent the energies i
!        H2 and the other HF asymptotic valleys, respectively.
!   Default Parameter Values:
!      Variable      Default value
!      requested_derivative_order             1
!      NFLAG(18)        6
!       First:  Calculate the dynamical correlation energy e(corr)
!       Second: Calculate the Extended-Hartree-Fock type energy given
!               by the LEPS function calibrated from SEC energies e(ehf)
!               etot    = e(ehf) + e(corr) + deh2
!               e(ehf)  = e(leps) + e(3c)
!               e(corr) = e(corr2) + e(corr3)
!       Note:   The zero of energy for this surface is F infinitely far
!               from H2 at r(H2) = re(H2) "deh2"
!               To get the zero of energy r1 must be at least 50 a.u.

!> @brief Implement FH2_GEN_5SEC_1991, a potential energy surface for FH2.
module chempotpy_fh2_gen_5sec_1991_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the FH2_GEN_5SEC_1991 potential energy surface for FH2.
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
        real(wp)              :: a1(3)
        real(wp)              :: a2(3)
        real(wp)              :: a3(3)
        real(wp)              :: ad1
        real(wp)              :: ad2
        real(wp)              :: ad3
        real(wp)              :: alf
        real(wp)              :: alfp
        real(wp)              :: alft
        real(wp)              :: alph0
        real(wp)              :: alph1
        real(wp)              :: bd1
        real(wp)              :: bd2
        real(wp)              :: bd3
        real(wp)              :: bet0
        real(wp)              :: bet1
        real(wp)              :: betp
        real(wp)              :: bhf(5)
        real(wp)              :: bhh(5)
        real(wp)              :: c3c
        real(wp)              :: input_cartesian(natom,3)
        real(wp)              :: cartesian_coordinates(natom,3)
        real(wp)              :: cfhh
        real(wp)              :: chf(10)
        real(wp)              :: chh(10)
        real(wp)              :: chhf
        real(wp)              :: ck1(10)
        real(wp)              :: ck2(10)
        real(wp)              :: ck3(10)
        real(wp)              :: ckp1(10)
        real(wp)              :: ckp2(10)
        real(wp)              :: ckp3(10)
        real(wp)              :: dcorr2(3)
        real(wp)              :: dcorr3(3)
        real(wp)              :: ddp1(3)
        real(wp)              :: ddp2(3)
        real(wp)              :: ddp3(3)
        real(wp)              :: dec(3)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        real(wp)              :: dleps(3)
        real(wp)              :: dp1(3)
        real(wp)              :: dp2(3)
        real(wp)              :: dp3(3)
        real(wp)              :: drtanh(3,10)
        real(wp)              :: dsw(3)
        real(wp)              :: e3c
        real(wp)              :: easyab
        real(wp)              :: easyac
        real(wp)              :: easybc
        real(wp)              :: ecorr2(3)
        real(wp)              :: eleps
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: ground_state_energy
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: eta1(10)
        real(wp)              :: eta2(10)
        real(wp)              :: eta3(10)
        real(wp)              :: etap1(10)
        real(wp)              :: etap2(10)
        real(wp)              :: etap3(10)
        real(wp)              :: ezero(isurf+1)
        real(wp)              :: fcsi
        real(wp)              :: fcsi1
        real(wp)              :: fcsi2
        real(wp)              :: gamma(3)
        integer               :: coordinate_representation
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        real(wp)              :: dissociation_energies(3)
        integer               :: atom_count
        real(wp)              :: angular_correction_coefficients(4)
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: p3c
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: q3c
        real(wp)              :: r(n3atom)
        real(wp)              :: r0(3)
        real(wp)              :: rb(3)
        character(len=75)     :: ref(5)
        real(wp)              :: rm(3)
        real(wp)              :: rmt(3)
        real(wp)              :: ro(3)
        real(wp)              :: rtanh(3,10)
        real(wp)              :: sw(3)
        real(wp)              :: thbc
        real(wp)              :: twobc
        real(wp)              :: xhf(3)
        real(wp)              :: xhh(3)
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        a1                               = 0.e0_wp; a2 = 0.e0_wp; a3 = 0.e0_wp; ad1 = 0.e0_wp; ad2 = 0.e0_wp; ad3 = 0.e0_wp; alf = 0.e0_wp; alfp = 0.e0_wp; alft = 0.e0_wp
        alph0                            = 0.e0_wp; alph1 = 0.e0_wp; bd1 = 0.e0_wp; bd2 = 0.e0_wp; bd3 = 0.e0_wp; bet0 = 0.e0_wp; bet1 = 0.e0_wp
        betp                             = 0.e0_wp; bhf = 0.e0_wp; bhh = 0.e0_wp; c3c = 0.e0_wp; input_cartesian = 0.e0_wp; cartesian_coordinates = 0.e0_wp; cfhh = 0.e0_wp; chf = 0.e0_wp; chh = 0.e0_wp
        chhf                             = 0.e0_wp; ck1 = 0.e0_wp; ck2 = 0.e0_wp; ck3 = 0.e0_wp; ckp1 = 0.e0_wp; ckp2 = 0.e0_wp; ckp3 = 0.e0_wp
        dcorr2                           = 0.e0_wp; dcorr3 = 0.e0_wp; ddp1 = 0.e0_wp; ddp2 = 0.e0_wp; ddp3 = 0.e0_wp; dec = 0.e0_wp
        excited_state_internal_gradients = 0.e0_wp; ground_state_internal_gradient = 0.e0_wp; diabatic_internal_gradients = 0.e0_wp; output_excited_state_gradients = 0.e0_wp; excited_state_cartesian_gradients = 0.e0_wp; output_ground_state_gradient = 0.e0_wp; ground_state_cartesian_gradient = 0.e0_wp
        output_diabatic_gradients        = 0.e0_wp; diabatic_cartesian_gradients = 0.e0_wp; dleps = 0.e0_wp; dp1 = 0.e0_wp; dp2 = 0.e0_wp; dp3 = 0.e0_wp; drtanh = 0.e0_wp; dsw = 0.e0_wp
        e3c                              = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp; ecorr2 = 0.e0_wp; eleps = 0.e0_wp; excited_state_energies = 0.e0_wp
        ground_state_energy              = 0.e0_wp; diabatic_couplings = 0.e0_wp; eta1 = 0.e0_wp; eta2 = 0.e0_wp; eta3 = 0.e0_wp; etap1 = 0.e0_wp; etap2 = 0.e0_wp; etap3 = 0.e0_wp
        ezero                            = 0.e0_wp; fcsi = 0.e0_wp; fcsi1 = 0.e0_wp; fcsi2 = 0.e0_wp; gamma = 0.e0_wp; coordinate_representation = 0
        last_requested_excited_state     = 0; last_requested_coupling = 0; first_requested_excited_state = 0; first_requested_coupling = 0; dissociation_energies = 0.e0_wp; atom_count = 0
        angular_correction_coefficients  = 0.e0_wp; maximum_derivative_order = 0; maximum_excited_state = 0; surface_requests = 0; requested_derivative_order = 0; p3c = 0.e0_wp; output_excited_state_energies = 0.e0_wp
        output_ground_state_energy       = 0.e0_wp; output_diabatic_couplings = 0.e0_wp; q3c = 0.e0_wp; r = 0.e0_wp; r0 = 0.e0_wp; rb = 0.e0_wp
        ref                              = ' '
        rm                               = 0.e0_wp
        rmt                              = 0.e0_wp
        ro                               = 0.e0_wp
        rtanh                            = 0.e0_wp
        sw                               = 0.e0_wp
        thbc                             = 0.e0_wp
        twobc                            = 0.e0_wp
        xhf                              = 0.e0_wp
        xhh                              = 0.e0_wp
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 0
        coordinate_representation  = 3
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 3
        xhf                        = [2.633449e0_wp, 0.417523e0_wp, -0.090289e0_wp]
        xhh                        = [-0.393828e0_wp, 0.167810e0_wp, -0.017081e0_wp]
        angular_correction_coefficients                   = [-0.071547e0_wp, -0.108550e0_wp, 0.050461e0_wp, - &
            & 0.370255e0_wp]
        c3c                      = 0.296668e0_wp
        alfp                     = 0.078606e0_wp
        alf                      = 0.18e0_wp
        alft                     = 2.14e0_wp
        p3c                      = 0.95e0_wp
        q3c                      = 1.0e0_wp
        betp                     = 0.18e0_wp
        dec                      = [2.0e0_wp, 2.0e0_wp, 2.0e0_wp]
        rb                       = [6.80e0_wp, 5.00e0_wp, 6.80e0_wp]
        bhh                      = [0.448467e0_wp, -0.056687e0_wp, 0.246831e0_wp, -0.018419e0_wp, 0.000598e0_wp]
        bhf                      = [0.9149555e0_wp, 0.3222197e0_wp, 0.1273001e0_wp, 0.0e0_wp, 0.0e0_wp]
        chh(6)                   = 6.499027e0_wp
        chh(8)                   = 1.243991e2_wp
        chh(10)                  = 3.2858e3_wp
        chf(6)                   = 6.68799e0_wp
        chf(8)                   = 1.0289812e2_wp
        chf(10)                  = 2.07391451e3_wp
        rm(1)                    = 1.7329e0_wp
        rm(2)                    = 1.4010e0_wp
        rm(3)                    = 1.7329e0_wp
        r0(1)                    = 5.9e0_wp
        r0(2)                    = 6.9283e0_wp
        r0(3)                    = 5.9e0_wp
        ro(1)                    = 1.7329e0_wp
        ro(2)                    = 1.449e0_wp
        ro(3)                    = 1.7329e0_wp
        alph0                    = 2.59528e1_wp
        alph1                    = 1.1868e0_wp
        bet0                     = 1.57381e1_wp
        bet1                     = 9.729e-2_wp
        eta1(6)                  = 0.9438421608e0_wp
        eta1(8)                  = 0.9774440533e0_wp
        eta1(10)                 = 0.9452240321e0_wp
        eta2(6)                  = 0.1085550150e1_wp
        eta2(8)                  = 0.1117269260e1_wp
        eta2(10)                 = 0.1138536961e1_wp
        eta3(6)                  = 0.9438421608e0_wp
        eta3(8)                  = 0.9774440533e0_wp
        eta3(10)                 = 0.9452240321e0_wp
        etap1(6)                 = 6.0e0_wp
        etap1(8)                 = 6.0e0_wp
        etap1(10)                = 6.0e0_wp
        etap2(6)                 = 6.0e0_wp
        etap2(8)                 = 6.0e0_wp
        etap2(10)                = 6.0e0_wp
        etap3(6)                 = 6.0e0_wp
        etap3(8)                 = 6.0e0_wp
        etap3(10)                = 6.0e0_wp
        ck1(6)                   = -0.2610767389e-1_wp
        ck1(8)                   = -0.3428701964e-1_wp
        ck1(10)                  = -0.4858536634e-1_wp
        ck2(6)                   = -0.9709847270e-1_wp
        ck2(8)                   = -0.1248186655e0_wp
        ck2(10)                  = -0.1426680807e0_wp
        ck3(6)                   = -0.2610767389e-1_wp
        ck3(8)                   = -0.3428701964e-1_wp
        ck3(10)                  = -0.4858536634e-1_wp
        ckp1(6)                  = 0.9438421608e0_wp
        ckp1(8)                  = 0.9774440533e0_wp
        ckp1(10)                 = 0.9452240321e0_wp
        ckp2(6)                  = 0.1085550150e1_wp
        ckp2(8)                  = 0.1117269260e1_wp
        ckp2(10)                 = 0.1138536961e1_wp
        ckp3(6)                  = 0.9438421608e0_wp
        ckp3(8)                  = 0.9774440533e0_wp
        ckp3(10)                 = 0.9452240321e0_wp
        dissociation_energies(1) = 0.19383609e0_wp
        a1(1)                    = 2.4822787e0_wp
        a2(1)                    = 1.5435337e0_wp
        a3(1)                    = 0.83093855e0_wp
        gamma(1)                 = 2.3992999e0_wp
        dissociation_energies(2) = 0.15796326e0_wp
        a1(2)                    = 2.1977034e0_wp
        a2(2)                    = 1.2932502e0_wp
        a3(2)                    = 0.64375666e0_wp
        gamma(2)                 = 2.1835071e0_wp
        dissociation_energies(3) = 0.19383609e0_wp
        a1(3)                    = 2.4822787e0_wp
        a2(3)                    = 1.5435337e0_wp
        a3(3)                    = 0.83093855e0_wp
        gamma(3)                 = 2.3992999e0_wp
        rmt                      = [1.0e5_wp, 1.0e5_wp, 1.0e5_wp]
        v                        = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        input_cartesian          = 0.e0_wp
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

        !> @brief Evaluate prepare_local_context for the FH2_GEN_5SEC_1991 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            real(wp), parameter :: deh2 = 0.174472692e0_wp
            real(wp), parameter :: dehf = 0.224989333e0_wp
            !      SAVE
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            !     CONSTANTS NEEDED FOR THE DAMPING FUNCTION
            !     A=ALPH0/FLOAT(N)**ALPH1        N=6,8,10
            !     B=BET0*EXP(-BET1*FLOAT(N))     N=6,8,10
            ad1 = alph0/6.0e0_wp**alph1
            ad2 = alph0/8.0e0_wp**alph1
            ad3 = alph0/10.0e0_wp**alph1
            bd1 = bet0*exp(-bet1*6.0e0_wp)
            bd2 = bet0*exp(-bet1*8.0e0_wp)
            bd3 = bet0*exp(-bet1*10.0e0_wp)
            !   Initialize the energy in the asymptotic valleys
            easyab   = dehf
            easybc   = deh2
            easyac   = dehf
            ezero(1) = deh2
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1) = 'G. C. Lynch, R. Steckler, D.W. Schwenke,'
            ref(2) = 'A.J.C. Varandas, D.G. Truhlar, and B.C. Garrett,'
            ref(3) = 'J. Chem. Phys. 94, 7136(1991)'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the FH2_GEN_5SEC_1991 potential.
        subroutine pot
            !   Surfaces:
            !      ground electronic state
            !   Zero of energy:
            !      The classical potential energy is set equal to zero for the F
            !      infinitely far from the H2 diatomic and R(H2) set equal to the
            !      H2 equilibrium diatomic value.
            !   Coordinates:
            !      Internal, Definition: R(1) = R(F-H)
            !                            R(2) = R(H-H)
            !                            R(3) = R(H-F)
            !        The energy in the AB valley, EASYAB, is equal to the energy
            !        of the C atom "infinitely" far from the AB diatomic and R(AB) s
            !        equal to Re(AB), the equilibrium bond length for the AB diatomi
            !        In this potential the AB valley represents H infinitely far fro
            !        the FH diatomic and R(FH) equal to Re(FH).
            !        Similarly, the terms EASYBC and EASYAC represent the energies i
            !        H2 and the other HF asymptotic valleys, respectively.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dec1, dec2, dec3, ec
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   Check the value of requested_derivative_order
            !  PART IIa.  CALCULATE THE DYNAMICAL CORRELATION ENERGY AND ITS DERIVAT
            !     FIRST: CALL THE SUBROUTINE DAMP WHICH CALCULATES THE DAMPING FUNCT
            !            FOR THE NTH DISPERSION COEFFICIENT AND ITS DERIVATIVES
            call damp(r, rm)
            !     SECOND: CALL THE SUBROUTINE CORR2 WHICH CALCULATES THE TWO BODY
            !             CORRELATION ENERGY AND ITS DERIVATIVES
            call corr2
            !     THIRD:  CALL XTANH TO COMPUTE HYPERBOLIC TANGENT VALUES NEEDED
            !             IN THE SUBROUTINE CORR3, THEN
            !             CALL THE SUBROUTINE CORR3 WHICH CALCULATES THE THREE BODY
            !             CORRELATION ENERGY AND ITS DERIVATIVES
            call xtanh(1, eta1)
            call xtanh(2, eta2)
            call xtanh(3, eta3)
            call corr3
            !     FINALLY: SUM THE TWO AND THE THREE BODY CORRELATION TERMS TO FORM
            !              AND DE(CORR)- THE DERIVATIVE OF E(CORR)
            ec = twobc + thbc
            if( requested_derivative_order .eq. 1 )then
                dec1 = dcorr2(1) + dcorr3(1)
                dec2 = dcorr2(2) + dcorr3(2)
                dec3 = dcorr2(3) + dcorr3(3)
            endif
            !       CALCULATE THE SWITCHING FUNCTION VALUES AND
            !       THE CORRESPONDING DERIVATIVES
            call switch
            call vleps
            ground_state_energy = ec + eleps + ezero(1)
            if( requested_derivative_order .ne. 1 )then
                call eunitzero
                if( requested_derivative_order .ne. 0 )then
                    call rtocart
                    call dedcou
                endif
                return
            endif
            ground_state_internal_gradient(1) = dec1 + dleps(1)
            ground_state_internal_gradient(2) = dec2 + dleps(2)
            ground_state_internal_gradient(3) = dec3 + dleps(3)
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate CORR2 for the FH2_GEN_5SEC_1991 potential.
        subroutine corr2
            !     CALCULATES TWO-BODY CORRELATION ENERGY.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: r1p10, r1p6, r1p8, r2p10, r2p6, r2p8, r3p10, r3p6, r3p8
            r1p6      = r(1)**6
            r1p8      = r(1)**8
            r1p10     = r(1)**10
            r2p6      = r(2)**6
            r2p8      = r(2)**8
            r2p10     = r(2)**10
            r3p6      = r(3)**6
            r3p8      = r(3)**8
            r3p10     = r(3)**10
            ecorr2(1) = -dp1(1)*chf(6)/r1p6 -dp1(2)*chf(8)/r1p8 -dp1(3)*chf(10)/r1p10
            ecorr2(2) = -dp2(1)*chh(6)/r2p6 -dp2(2)*chh(8)/r2p8 -dp2(3)*chh(10)/r2p10
            ecorr2(3) = -dp3(1)*chf(6)/r3p6 -dp3(2)*chf(8)/r3p8 -dp3(3)*chf(10)/r3p10
            twobc     = ecorr2(1)+ecorr2(2)+ecorr2(3)
            if( requested_derivative_order .ne. 1 ) return
            dcorr2(1) = -chf(6)*(ddp1(1)-6.0e0_wp*dp1(1)/r(1))/r1p6 -chf(8)*(ddp1(2)-8.0e0_wp*dp1(2)/r(1))/r1p8 - &
                & chf(10)* (ddp1(3)- 10.0e0_wp*dp1(3)/r(1))/r1p10
            dcorr2(2) = -chh(6)*(ddp2(1)-6.0e0_wp*dp2(1)/r(2))/r2p6 -chh(8)*(ddp2(2)-8.0e0_wp*dp2(2)/r(2))/r2p8 - &
                & chh(10)* (ddp2(3)- 10.0e0_wp*dp2(3)/r(2))/r2p10
            dcorr2(3) = -chf(6)*(ddp3(1)-6.0e0_wp*dp3(1)/r(3))/r3p6 -chf(8)*(ddp3(2)-8.0e0_wp*dp3(2)/r(3))/r3p8 - &
                & chf(10)* (ddp3(3)- 10.0e0_wp*dp3(3)/r(3))/r3p10
        end subroutine corr2

        !> @brief Evaluate DAMP for the FH2_GEN_5SEC_1991 potential.
        subroutine damp( rt1, rt2 )
            !     CALCULATES DAMPING FUNCTION FOR THE NTH DISPERSION COEFFICIENT
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dpol1(3), dpol2(3), dpol3(3), dx(3), pol1(3), pol2(3), pol3(3), rt1(3), rt2(3), x(3)
            x(1) = 2.0e0_wp*rt1(1)/(rt2(1)+2.5e0_wp*r0(1))
            x(2) = 2.0e0_wp*rt1(2)/(rt2(2)+2.5e0_wp*r0(2))
            x(3) = 2.0e0_wp*rt1(3)/(rt2(3)+2.5e0_wp*r0(3))
            if( requested_derivative_order .eq. 1 )then
                dx(1) = 2.0e0_wp/(rt2(1)+2.5e0_wp*r0(1))
                dx(2) = 2.0e0_wp/(rt2(2)+2.5e0_wp*r0(2))
                dx(3) = 2.0e0_wp/(rt2(3)+2.5e0_wp*r0(3))
            endif
            pol1(1) = ad1*x(1)+bd1*x(1)**2
            pol2(1) = ad2*x(1)+bd2*x(1)**2
            pol3(1) = ad3*x(1)+bd3*x(1)**2
            pol1(2) = ad1*x(2)+bd1*x(2)**2
            pol2(2) = ad2*x(2)+bd2*x(2)**2
            pol3(2) = ad3*x(2)+bd3*x(2)**2
            pol1(3) = ad1*x(3)+bd1*x(3)**2
            pol2(3) = ad2*x(3)+bd2*x(3)**2
            pol3(3) = ad3*x(3)+bd3*x(3)**2
            if( requested_derivative_order .eq. 1 )then
                dpol1(1) = ad1+2.0e0_wp*bd1*x(1)
                dpol2(1) = ad2+2.0e0_wp*bd2*x(1)
                dpol3(1) = ad3+2.0e0_wp*bd3*x(1)
                dpol1(2) = ad1+2.0e0_wp*bd1*x(2)
                dpol2(2) = ad2+2.0e0_wp*bd2*x(2)
                dpol3(2) = ad3+2.0e0_wp*bd3*x(2)
                dpol1(3) = ad1+2.0e0_wp*bd1*x(3)
                dpol2(3) = ad2+2.0e0_wp*bd2*x(3)
                dpol3(3) = ad3+2.0e0_wp*bd3*x(3)
            endif
            dp1(1) = (1.0e0_wp-exp(-pol1(1)))**6
            dp1(2) = (1.0e0_wp-exp(-pol2(1)))**8
            dp1(3) = (1.0e0_wp-exp(-pol3(1)))**10
            dp2(1) = (1.0e0_wp-exp(-pol1(2)))**6
            dp2(2) = (1.0e0_wp-exp(-pol2(2)))**8
            dp2(3) = (1.0e0_wp-exp(-pol3(2)))**10
            dp3(1) = (1.0e0_wp-exp(-pol1(3)))**6
            dp3(2) = (1.0e0_wp-exp(-pol2(3)))**8
            dp3(3) = (1.0e0_wp-exp(-pol3(3)))**10
            if( requested_derivative_order .ne. 1 ) return
            ddp1(1) = 6.0e0_wp*dp1(1)*exp(-pol1(1))*dpol1(1)*dx(1)/ (1.0e0_wp-exp(-pol1(1)))
            ddp1(2) = 8.0e0_wp*dp1(2)*exp(-pol2(1))*dpol2(1)*dx(1)/ (1.0e0_wp-exp(-pol2(1)))
            ddp1(3) = 10.0e0_wp*dp1(3)*exp(-pol3(1))*dpol3(1)*dx(1)/ (1.0e0_wp-exp(-pol3(1)))
            ddp2(1) = 6.0e0_wp*dp2(1)*exp(-pol1(2))*dpol1(2)*dx(2)/ (1.0e0_wp-exp(-pol1(2)))
            ddp2(2) = 8.0e0_wp*dp2(2)*exp(-pol2(2))*dpol2(2)*dx(2)/ (1.0e0_wp-exp(-pol2(2)))
            ddp2(3) = 10.0e0_wp*dp2(3)*exp(-pol3(2))*dpol3(2)*dx(2)/ (1.0e0_wp-exp(-pol3(2)))
            ddp3(1) = 6.0e0_wp*dp3(1)*exp(-pol1(3))*dpol1(3)*dx(3)/ (1.0e0_wp-exp(-pol1(3)))
            ddp3(2) = 8.0e0_wp*dp3(2)*exp(-pol2(3))*dpol2(3)*dx(3)/ (1.0e0_wp-exp(-pol2(3)))
            ddp3(3) = 10.0e0_wp*dp3(3)*exp(-pol3(3))*dpol3(3)*dx(3)/ (1.0e0_wp-exp(-pol3(3)))
        end subroutine damp

        !> @brief Evaluate CORR3 for the FH2_GEN_5SEC_1991 potential.
        subroutine corr3
            !     CALCULATES THREE-BODY CORRELATION ENERGY.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: co3(3), dc3dr1(3), dc3dr2(3), dc3dr3(3), dg1, dg2, dg3, dh1, dh2, dh3, dt1dr1, dt1dr2, &
                & dt1dr3, dt2dr1, dt2dr2, dt2dr3, dt3dr1, dt3dr2, dt3dr3, fi, g1, g2, g3, h1, h2, h3, r1pi, r2pi, r3pi, t1, t2, t3
            integer :: i, jj
            jj = 0
            do i=6, 10, 2
                jj   = jj+1
                fi   = real(i, wp)
                r1pi = r(1)**i
                r2pi = r(2)**i
                r3pi = r(3)**i
                g1   = 1.0e0_wp+ck1(i)*exp(-ckp1(i)*(r(1)-ro(1)))
                g2   = 1.0e0_wp+ck2(i)*exp(-ckp2(i)*(r(2)-ro(2)))
                g3   = 1.0e0_wp+ck3(i)*exp(-ckp3(i)*(r(3)-ro(3)))
                if( requested_derivative_order .eq. 1 )then
                    dg1 = -ckp1(i)*(g1-1.0e0_wp)
                    dg2 = -ckp2(i)*(g2-1.0e0_wp)
                    dg3 = -ckp3(i)*(g3-1.0e0_wp)
                endif
                h1 = rtanh(1,i)**etap1(i)
                h2 = rtanh(2,i)**etap2(i)
                h3 = rtanh(3,i)**etap3(i)
                if( requested_derivative_order .eq. 1 )then
                    dh1 = etap1(i)*(rtanh(1,i)**(etap1(i)-1))*drtanh(1,i)
                    dh2 = etap2(i)*(rtanh(2,i)**(etap2(i)-1))*drtanh(2,i)
                    dh3 = etap3(i)*(rtanh(3,i)**(etap3(i)-1))*drtanh(3,i)
                endif
                t1 = chf(i)*(1.0e0_wp-0.5e0_wp*(g2*h3+g3*h2)) *dp1(jj)/r1pi
                t2 = chh(i)*(1.0e0_wp-0.5e0_wp*(g3*h1+g1*h3)) *dp2(jj)/r2pi
                t3 = chf(i)*(1.0e0_wp-0.5e0_wp*(g1*h2+g2*h1)) *dp3(jj)/r3pi
                if( requested_derivative_order .eq. 1 )then
                    dt1dr1 = t1*(ddp1(jj)/dp1(jj)-fi/r(1))
                    dt1dr2 = -0.5e0_wp*chf(i)*dp1(jj)*(dg2*h3+g3*dh2)/r1pi
                    dt1dr3 = -0.5e0_wp*chf(i)*dp1(jj)*(g2*dh3+dg3*h2)/r1pi
                    dt2dr1 = -0.5e0_wp*chh(i)*dp2(jj)*(g3*dh1+dg1*h3)/r2pi
                    dt2dr2 = t2*ddp2(jj)/dp2(jj)-t2*fi/r(2)
                    dt2dr3 = -0.5e0_wp*chh(i)*dp2(jj)*(dg3*h1+g1*dh3)/r2pi
                    dt3dr1 = -0.5e0_wp*chf(i)*dp3(jj)*(dg1*h2+g2*dh1)/r3pi
                    dt3dr2 = -0.5e0_wp*chf(i)*dp3(jj)*(g1*dh2+dg2*h1)/r3pi
                    dt3dr3 = t3*(ddp3(jj)/dp3(jj)-fi/r(3))
                endif
                co3(jj) = t1+t2+t3
                if( requested_derivative_order .eq. 1 )then
                    dc3dr1(jj) = dt1dr1+dt2dr1+dt3dr1
                    dc3dr2(jj) = dt1dr2+dt2dr2+dt3dr2
                    dc3dr3(jj) = dt1dr3+dt2dr3+dt3dr3
                endif
            enddo
            thbc = co3(1)+co3(2)+co3(3)
            if( requested_derivative_order .ne. 1 ) return
            dcorr3(1) = dc3dr1(1)+dc3dr1(2)+dc3dr1(3)
            dcorr3(2) = dc3dr2(1)+dc3dr2(2)+dc3dr2(3)
            dcorr3(3) = dc3dr3(1)+dc3dr3(2)+dc3dr3(3)
        end subroutine corr3

        !> @brief Evaluate VLEPS for the FH2_GEN_5SEC_1991 potential.
        subroutine vleps
            !     CALCULATES THREE-BODY EXTENDED-HARTREE-FOCK ENERGY DEFINED BY A
            !     LEPS-TYPE FUNCTION
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: arg, csi, csip2, csip3, csip4, darg, dcsi(3), de3c1, de3c2, de3c3, dehf(3), dehft(3), &
                & detrip(3, 3), dex(3, 3), dfcsi1(3), dfcsi2(3), dq(3, 3), dr, dri(3), drp2, dscale, dsumq(3), dti(3, 3), dx(3), ehf(3), &
                & ehft(3), etrip(3), ex(3), exch, fcsi1t, fcsi2t, ff1, ff2, ff3, q(3), r1p2, r2p2, r2p3, r2p4, r3p2, scale, sumq, &
                & texch, ti(3), tmp1, tmp2, x(3), xip2, xterm, xterm1, xterm2, xterm3
            integer :: i, j
            !     CSI IS SIMILAR TO A COSINE FUNCTION
            csi   = (r(1)-r(3))/r(2)
            csip2 = csi**2
            csip3 = csi**3
            csip4 = csi**4
            if( requested_derivative_order .eq. 1 )then
                dcsi(1) = 1.0e0_wp/r(2)
                dcsi(2) = -csi/r(2)
                dcsi(3) = -1.0e0_wp/r(2)
            endif
            fcsi1t = 1.0e0_wp+angular_correction_coefficients(1)*(1.0e0_wp-csip2)+angular_correction_coefficients(2)* &
                & (1.0e0_wp-csip4)
            fcsi2t = 1.0e0_wp+angular_correction_coefficients(3)*(1.0e0_wp-csip2)+angular_correction_coefficients(4)* &
                & (1.0e0_wp-csip4)
            fcsi1 = fcsi1t**2
            fcsi2 = fcsi2t**2
            if( requested_derivative_order .eq. 1 )then
                tmp1 = 2.0e0_wp*fcsi1t*(-2.0e0_wp*angular_correction_coefficients(1)*csi-4.0e0_wp* &
                    & angular_correction_coefficients(2)*csip3)
                tmp2 = 2.0e0_wp*fcsi2t*(-2.0e0_wp*angular_correction_coefficients(3)*csi-4.0e0_wp* &
                    & angular_correction_coefficients(4)*csip3)
                do i=1, 3
                    dfcsi1(i) = tmp1*dcsi(i)
                    dfcsi2(i) = tmp2*dcsi(i)
                enddo
            endif
            r1p2 = r(1)**2
            r2p2 = r(2)**2
            r2p3 = r(2)**3
            r2p4 = r(2)**4
            r3p2 = r(3)**2
            !     CALCULATE THE SWITCHING FUNCTION FOR THE DIATOMIC TRIPLET STATE
            !     CURVES AND ITS DERIVATIVES : SW(I),DSW(I)
            !     CALCULATE THE EXTENDED-HARTREE-FOCK CURVE FOR THE GROUND SINGLET
            !     STATE OF H2 OR HF AND ITS DERIVATIVES : EHF(I),DEHF(I)
            !     I => R(J) J=1,2,3
            !                         A.J.C. VARANDAS & J.D. SILVA,
            !                         J. CHEM. SOC. FARADAY TRANS. II,
            !                         82,593 (1986)
            !     FOR AN UPDATE FOR C8 AND C10 DISPERSION COEFFICIENTS OF HF :
            !                         A.J.C. VARANDAS
            !                         ADV. CHEM. PHYS., 1988 IN PRESS
            do i=1, 3
                dr     = r(i)-rm(i)
                drp2   = dr**2
                ehf(i) = -dissociation_energies(i)*(1.0e0_wp+a1(i)*dr+a2(i)*drp2+a3(i)*drp2*dr) *exp(-gamma(i)*dr)
                if( requested_derivative_order .eq. 1 ) dehf(i) = -dissociation_energies(i)*(a1(i)+2.0e0_wp*a2(i)*dr+ &
                    & 3.0e0_wp* a3(i)*drp2)*exp(-gamma(i)* dr)- gamma(i)*ehf(i)
            enddo
            !     CALCULATE THE DIATOMIC TRIPLET STATE CURVES ETRIP(I)
            !     FIRST: CALCULATE ETRIP(I) FOR I=2 => FOR H2
            !            SCALE IS THE SCALING FUNCTION FOR THE DIATOMIC
            !            TRIPLET STATE CURVES USED FOR H2 ONLY
            !            EHFT(I) IS THE POTENTIAL FOR THE LOWEST
            !            TRIPLET STATE CURVE
            !                 FOR H2 : A.J.C. VARANDAS & J. BARANDO
            !                          MOL. PHYS. 45,857 (1982)
            !            H2 => R(2)
            scale = 1.0e0_wp+xhh(1)*r(2)+xhh(2)*r2p2+xhh(3)*r2p3
            if( requested_derivative_order .eq. 1 ) dscale = xhh(1)+2.0e0_wp*xhh(2)*r(2)+&
                &3.0e0_wp*xhh(3)*r2p2
            arg = bhh(2)*r(2)+bhh(3)*r2p2+bhh(4)*r2p3+bhh(5)*r2p4
            if( requested_derivative_order .eq. 1 ) darg = bhh(2)+2.0e0_wp*bhh(3)*r(2)+&
                &3.0e0_wp*bhh(4)*r2p2+4.0e0_wp*bhh(5)*r2p3
            ehft(2)                                          = bhh(1)*exp(-arg)/r(2)
            if( requested_derivative_order .eq. 1 ) dehft(2) = ehft(2)*(-darg-1.0e0_wp/r(2))
            etrip(2)                                         = scale*ehft(2)
            if( requested_derivative_order .eq. 1 )then
                detrip(2, 2) = dscale*ehft(2)+scale*dehft(2)
                detrip(2, 1) = 0.0e0_wp
                detrip(2, 3) = 0.0e0_wp
            endif
            !     NEXT:
            !            CALCULATE THE EHF POTENTIAL FOR THE LOWEST SIGMA TRIPLET
            !            STATE OF HF;
            !            FIT TO ab initio DATA FROM
            !                                      T.H. DUNNING
            !                                      J. CHEM. PHYS.  65,3854 (1976)
            !            THEN CALCULATE THE TRIPLET STATE CURVES FOR HF
            do i=1, 3, 2
                arg = bhf(2)*r(i)+bhf(3)*r(i)**2
                if( requested_derivative_order .eq. 1 ) darg = bhf(2)+2.0e0_wp*bhf(3)*r(i)
                ehft(i)                                          = bhf(1)*exp(-arg)/r(i)
                if( requested_derivative_order .eq. 1 ) dehft(i) = ehft(i)*(-1.0e0_wp/r(i)-darg)
                x(i)                                             = exp(-xhf(2)*r(i))
                if( requested_derivative_order .eq. 1 ) dx(i)    = -xhf(2)*x(i)
                dri(i)                                           = 1.0e0_wp
                xip2                                             = x(i)**2
                etrip(i)                                         = xhf(1)*(fcsi1*xip2+fcsi2*xhf(3)*x(i))/r(i)
                if( requested_derivative_order .ne. 1 ) cycle
                do j=1, 3
                    if( j .ne. i )then
                        dx(j)  = 0.0e0_wp
                        dri(j) = 0.0e0_wp
                    endif
                    detrip(i, j) = xhf(1)*(dfcsi1(j)*xip2+2.0e0_wp*fcsi1*x(i)*dx(j)+ dfcsi2(j)*xhf(3)*x(i)+fcsi2* &
                        & xhf(3)* dx(j))/ r(i)- etrip(i)*dri(j)/r(i)
                enddo
            enddo
            !     FINALLY:
            !             COMBINE THE THREE-BODY EXTENDED-HARTREE-FOCK ENERGY BY A
            !             LEPS TYPE FUNCTION
            sumq = 0.0e0_wp
            if( requested_derivative_order .eq. 1 )then
                dsumq(1) = 0.0e0_wp
                dsumq(2) = 0.0e0_wp
                dsumq(3) = 0.0e0_wp
            endif
            do i=1, 3
                ti(i) = (1.0e0_wp-sw(i))*etrip(i)+sw(i)*ehft(i)
                q(i)  = 0.5e0_wp*(ehf(i)+ti(i))
                ex(i) = 0.5e0_wp*(ehf(i)-ti(i))
                sumq  = sumq+q(i)
                if( requested_derivative_order .ne. 1 ) cycle
                do j=1, 3
                    if( j .eq. i )then
                        dti(i, j) = (1.0e0_wp-sw(i))*detrip(i, j)-etrip(i)*dsw(j)+ ehft(i)*dsw(j)+sw(i)*dehft(j)
                        dq(i, j)  = 0.5e0_wp*(dehf(j)+dti(i, j))
                        dex(i, j) = 0.5e0_wp*(dehf(j)-dti(i, j))
                    else
                        dti(i, j) = (1.0e0_wp-sw(i))*detrip(i, j)
                        dq(i, j)  = 0.5e0_wp*dti(i, j)
                        dex(i, j) = -dq(i, j)
                    endif
                enddo
            enddo
            if( requested_derivative_order .eq. 1 )then
                dsumq(1) = dq(1, 1)+dq(2, 1)+dq(3, 1)
                dsumq(2) = dq(1, 2)+dq(2, 2)+dq(3, 2)
                dsumq(3) = dq(1, 3)+dq(2, 3)+dq(3, 3)
            endif
            ff1   = (ex(1) - ex(2)) ** 2
            ff2   = (ex(2) - ex(3)) ** 2
            ff3   = (ex(3) - ex(1)) ** 2
            exch  = sqrt(0.5e0_wp*(ff1+ff2+ff3))
            eleps = sumq - exch
            goto_block_700_1004: block
                if( requested_derivative_order .ne. 1 ) exit goto_block_700_1004
                xterm1 = 2.0e0_wp*ex(1)-ex(2)-ex(3)
                xterm2 = 2.0e0_wp*ex(2)-ex(1)-ex(3)
                xterm3 = 2.0e0_wp*ex(3)-ex(1)-ex(2)
                do i=1, 3
                    texch = exch
                    if( texch .eq. 0.0e0_wp ) texch = 1.0e0_wp
                    xterm    = 0.5e0_wp*(xterm1*dex(1, i)+xterm2*dex(2, i)+ xterm3*dex(3, i))/texch
                    dleps(i) = dsumq(i)-xterm
                enddo
            end block goto_block_700_1004
            call thcent(de3c1, de3c2, de3c3)
            eleps = eleps + e3c
            if( requested_derivative_order .ne. 1 ) return
            dleps(1) = dleps(1) + de3c1
            dleps(2) = dleps(2) + de3c2
            dleps(3) = dleps(3) + de3c3
        end subroutine vleps

        !> @brief Evaluate THCENT for the FH2_GEN_5SEC_1991 potential.
        subroutine thcent( de3c1, de3c2, de3c3 )
            !  LOCAL THREE-BODY EHF-TYPE TERM; ADAPTED FROM BROWN ET AL, J. CHEM. PH
            !  82(1985)188 WITH MODIFICATIONS
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: costh, dcos1, dcos2, dcos3, de3c1, de3c2, de3c3, ex1, ex2, ex3, f3, omp3c, omq3c, r12, r1i, r22, r32, r3i, rdif, rdif2, rsum, t, t1, t2, t3, t4, t5
            logical  :: lderiv
            lderiv = .false.
            if( requested_derivative_order .eq. 1 ) lderiv = .true.
            omp3c = 1.e0_wp - p3c
            omq3c = 1.e0_wp - q3c
            r12   = r(1) * r(1)
            r22   = r(2) * r(2)
            r32   = r(3) * r(3)
            rsum  = r(1) + r(3)
            rdif  = r(1) - r(3)
            rdif2 = rdif * rdif
            t1    = -alfp * rsum
            ex1   = c3c * exp(t1 * rsum)
            t2    = -alf * rdif
            ex2   = omp3c * exp(t2 * rdif)
            t3    = -alft * rdif2
            ex3   = p3c * exp(t3 * rdif2)
            r1i   = 1.e0_wp / r(1)
            r3i   = 1.e0_wp / r(3)
            t4    = r1i*r3i
            !   COSTH IS THE ANGLE BETWEEN R1 AND R3
            costh = 0.5e0_wp * (r12 + r32 - r22) * t4
            t5    = omq3c * costh
            !   SET MODIFICATION ...
            t = q3c + t5 * costh
            !   E3C CONTAINS THE THREE CENTER CORRECTION TO V
            e3c = (ex2 + ex3) * t * ex1
            f3  = exp(-betp*(r(1)+r(3)-r(2))**2)
            e3c = e3c*f3
            goto_block_280_1096: block
                if( .not.lderiv ) exit goto_block_280_1096
                !   COMPUTE DERIVATIVE OF THE 3C TERM
                !   FIRST, DERIVATIVE OF COSTH
                dcos1 = r3i - costh*r1i
                dcos2 = - r(2)*t4
                dcos3 = r1i - costh*r3i
                t4    = ex2 + ex3
                t5    = t4 * t5
                !   NOW, DERIVATIVES OF E3C
                de3c1 = 2.e0_wp*((t2 * ex2 + 2.e0_wp * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos1) * ex1*f3+ &
                    & 2.0e0_wp*e3c*(- betp* (r(1)+ r(3)- r(2)))
                de3c2 = 2.e0_wp * t5 * dcos2 * ex1*f3+ 2.0e0_wp*e3c*(betp*(r(1)+r(3)-r(2)))
                de3c3 = 2.e0_wp*((-t2 * ex2 - 2.e0_wp * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos3) * ex1*f3+ &
                    & 2.0e0_wp*e3c*(- betp* (r(1)+ r(3)-r(2)))
            end block goto_block_280_1096
        end subroutine thcent

        !> @brief Evaluate XTANH for the FH2_GEN_5SEC_1991 potential.
        subroutine xtanh( i, a )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(*), ex, rval, sgn, t, th
            integer  :: i, j
            do j=6, 10, 2
                rval = r(i) * a(j)
                if( rval .gt. 0.0e0_wp )then
                    ex  = exp( -2.0e0_wp * rval)
                    sgn = 1.0e0_wp
                else
                    ex  = exp( 2.0e0_wp * rval)
                    sgn = -1.0e0_wp
                endif
                t                                                   = 1.0e0_wp + ex
                th                                                  = 1.0e0_wp/t
                if( requested_derivative_order .eq. 1 ) drtanh(i,j) = 4.0e0_wp*a(j)*ex*th*th
                rtanh(i,j)                                          = sgn * th * (1.0e0_wp - ex)
            enddo
        end subroutine xtanh

        !> @brief Evaluate SWITCH for the FH2_GEN_5SEC_1991 potential.
        subroutine switch
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dth, ex, rmrb, rval, sgn, t, th
            integer  :: i
            do i=1, 3
                rmrb = r(i)-rb(i)
                rval = rmrb*dec(i)
                if( rval .gt. 0.0e0_wp )then
                    ex  = exp(-2.0e0_wp*rval)
                    sgn = 1.0e0_wp
                else
                    ex  = exp(2.0e0_wp*rval)
                    sgn = -1.0e0_wp
                endif
                t  = 1.0e0_wp + ex
                th = 1.0e0_wp/t
                if( requested_derivative_order .eq. 1 ) dth = 4.0e0_wp*dec(i)*ex*th*th
                th                                             = sgn*th*(1.0e0_wp-ex)
                sw(i)                                          = 0.5e0_wp*(1.0e0_wp+th)
                if( requested_derivative_order .eq. 1 ) dsw(i) = 0.5e0_wp*dth
            enddo
        end subroutine switch

        !> @brief Evaluate POTINFO for the FH2_GEN_5SEC_1991 potential.
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

        !> @brief Evaluate CARTTOR for the FH2_GEN_5SEC_1991 potential.
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

        !> @brief Evaluate EUNITZERO for the FH2_GEN_5SEC_1991 potential.
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

        !> @brief Evaluate RTOCART for the FH2_GEN_5SEC_1991 potential.
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

        !> @brief Evaluate DEDCOU for the FH2_GEN_5SEC_1991 potential.
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
end module chempotpy_fh2_gen_5sec_1991_pes

!> @brief Expose FH2_GEN_5SEC_1991 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_fh2_gen_5sec_1991_pes, only: evaluate
    use chempotpy_kinds,                 only: wp
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
