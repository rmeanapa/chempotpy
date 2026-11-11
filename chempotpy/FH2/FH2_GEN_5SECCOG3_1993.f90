!   System:          FH2
!   Functional form: Double many-body expansion
!   Common name:     5SEC-C2
!   Reference:          unpublished
!   Cross reference:    S. L. Mielke, G. C. Lynch,
!                       diat_d. G. Truhlar, and diat_d. W. Schwenke
!                       Chem. Phys. Lett. 213, 10-16 (1993)
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
!      N.B.: R(1),[R(FH)], must be at least 50 a.u. for the energy to be
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
!        usricm_nder  = 0 => no derivatives are calculated
!        usricm_nder  = 1 => calculate first derivatives of the energy for the
!                     ground electronic state with respect to the coordi
!        usricm_nflag  - Control flags
!      usricm_nflag(18-20)
!        passes the FORTRAN unit number used for potential output
!      /ASYCM/ asycm_easyab, asycm_easybc, asycm_easyac
!        passes the energy in the three asymptotic valleys for an A+BC s
!        The energy in the AB valley, asycm_easyab, is equal to the energy
!        of the C atom "infinitely" far from the AB diatomic and R(AB) s
!        equal to Re(AB), the equilibrium bond length for the AB diatomi
!        In this potential the AB valley represents H infinitely far fro
!        the FH diatomic and R(FH) equal to Re(FH).
!        Similarly, the terms asycm_easybc and asycm_easyac represent the energies i
!        H2 and the other HF asymptotic valleys, respectively.
!   Default Parameter Values:
!      Variable      Default value
!      usricm_nder             1
!      usricm_nflag(18)        6
!   ENERGY  = E(EHF) + E(CORR) + De(H2)
!   E(EHF)  = E(London) + E(3-body)
!   E(CORR) = ECORR(2-body) + ECORR(3-body)

!> @brief Implement FH2_GEN_5SECCOG3_1993, a potential energy surface for FH2.
module chempotpy_fh2_gen_5seccog3_1993_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the FH2_GEN_5SECCOG3_1993 potential energy surface for FH2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 3
        integer, parameter            :: natom   = 25
        integer, parameter            :: isurf   = 5
        integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter            :: n3atom  = 75
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: asycm_easyab
        double precision              :: asycm_easyac
        double precision              :: asycm_easybc
        double precision              :: atat_chf(10)
        double precision              :: atat_chh(10)
        double precision              :: dampc_ad1
        double precision              :: dampc_ad2
        double precision              :: dampc_ad3
        double precision              :: dampc_bd1
        double precision              :: dampc_bd2
        double precision              :: dampc_bd3
        double precision              :: decm_disse(3)
        double precision              :: degsdrcm_dcorr2(3)
        double precision              :: degsdrcm_dcorr3(3)
        double precision              :: degsdrcm_de3c(3)
        double precision              :: degsdrcm_dldr(3)
        double precision              :: diat_a1(3)
        double precision              :: diat_a2(3)
        double precision              :: diat_a3(3)
        double precision              :: diat_d(3)
        double precision              :: diat_gamma(3)
        double precision              :: dispc_alph0
        double precision              :: dispc_alph1
        double precision              :: dispc_bet0
        double precision              :: dispc_bet1
        double precision              :: dist_r0(3)
        double precision              :: dist_rm(3)
        double precision              :: dpc_ddp1(3)
        double precision              :: dpc_ddp2(3)
        double precision              :: dpc_ddp3(3)
        double precision              :: dpc_dp1(3)
        double precision              :: dpc_dp2(3)
        double precision              :: dpc_dp3(3)
        double precision              :: engycm_e3c
        double precision              :: engycm_ecor2b
        double precision              :: engycm_ecor3b
        double precision              :: engycm_elond
        double precision              :: infocm_cartnu(natom,3)
        integer                       :: infocm_icartr
        integer                       :: infocm_indexes(natom)
        integer                       :: infocm_irctnt
        integer                       :: infocm_mder
        integer                       :: infocm_msurf
        integer                       :: infocm_natoms
        character                     :: infocm_ref(5)
        double precision              :: kkp_ck1(10)
        double precision              :: kkp_ck2(10)
        double precision              :: kkp_ck3(10)
        double precision              :: kkp_ckp1(10)
        double precision              :: kkp_ckp2(10)
        double precision              :: kkp_ckp3(10)
        double precision              :: local_alf
        double precision              :: local_alfp
        double precision              :: local_alft
        double precision              :: local_betp
        double precision              :: local_c3c
        double precision              :: local_p3c
        double precision              :: local_q3c
        double precision              :: par_bhf(5)
        double precision              :: par_bhh(5)
        double precision              :: pt1cm_degsdr(n3atom)
        double precision              :: pt1cm_engygs
        double precision              :: pt1cm_r(n3atom)
        double precision              :: pt3cm_ezero(isurf+1)
        double precision              :: pt4cm_deesdr(n3atom,isurf)
        double precision              :: pt4cm_engyes(isurf)
        double precision              :: pt5cm_deijdr(n3atom,jsurf)
        double precision              :: pt5cm_engyij(jsurf)
        double precision              :: ref_ro(3)
        double precision              :: swpar_dec(3)
        double precision              :: swpar_dsw(3,3)
        double precision              :: swpar_rb(3)
        double precision              :: swpar_sw(3)
        double precision              :: tpar_p(4)
        double precision              :: tpar_xhf(3)
        double precision              :: tpar_xhh(3)
        double precision              :: trial_eta1(10)
        double precision              :: trial_eta2(10)
        double precision              :: trial_eta3(10)
        double precision              :: trial_etap1(10)
        double precision              :: trial_etap2(10)
        double precision              :: trial_etap3(10)
        double precision              :: usricm_anuzero
        double precision              :: usricm_cart(natom,3)
        integer                       :: usricm_nasurf(isurf+1,isurf+1)
        integer                       :: usricm_nder
        integer                       :: usricm_nflag(20)
        integer                       :: usricm_nulbl(natom)
        double precision              :: usrocm_descart(natom,3,isurf)
        double precision              :: usrocm_dgscart(natom,3)
        double precision              :: usrocm_dijcart(natom,3,jsurf)
        double precision              :: usrocm_pengyes(isurf)
        double precision              :: usrocm_pengygs
        double precision              :: usrocm_pengyij(jsurf)
        double precision              :: utilcm_cnvrtd
        double precision              :: utilcm_cnvrtde
        double precision              :: utilcm_cnvrte
        double precision              :: utilcm_descartnu(natom,3,isurf)
        double precision              :: utilcm_dgscartnu(natom,3)
        double precision              :: utilcm_dijcartnu(natom,3,jsurf)
        integer                       :: utilcm_ireorder
        integer                       :: utilcm_kediag
        integer                       :: utilcm_keoffd
        integer                       :: utilcm_ksdiag
        integer                       :: utilcm_ksoffd
        double precision              :: vtanh_drtanh(3,10)
        double precision              :: vtanh_rtanh(3,10)
        double precision              :: v
        integer                       :: iatom, idir, istate
        !initialize
        asycm_easyab     = 0.d0; asycm_easyac = 0.d0; asycm_easybc = 0.d0; atat_chf = 0.d0; atat_chh = 0.d0
        dampc_ad1        = 0.d0; dampc_ad2 = 0.d0; dampc_ad3 = 0.d0; dampc_bd1 = 0.d0; dampc_bd2 = 0.d0; dampc_bd3 = 0.d0
        decm_disse       = 0.d0; degsdrcm_dcorr2 = 0.d0; degsdrcm_dcorr3 = 0.d0; degsdrcm_de3c = 0.d0; degsdrcm_dldr = 0.d0
        diat_a1          = 0.d0; diat_a2 = 0.d0; diat_a3 = 0.d0; diat_d = 0.d0; diat_gamma = 0.d0; dispc_alph0 = 0.d0
        dispc_alph1      = 0.d0; dispc_bet0 = 0.d0; dispc_bet1 = 0.d0; dist_r0 = 0.d0; dist_rm = 0.d0; dpc_ddp1 = 0.d0
        dpc_ddp2         = 0.d0; dpc_ddp3 = 0.d0; dpc_dp1 = 0.d0; dpc_dp2 = 0.d0; dpc_dp3 = 0.d0; engycm_e3c = 0.d0
        engycm_ecor2b    = 0.d0; engycm_ecor3b = 0.d0; engycm_elond = 0.d0; infocm_cartnu = 0.d0; infocm_icartr = 0
        infocm_indexes   = 0; infocm_irctnt = 0; infocm_mder = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref       = ' '
        kkp_ck1          = 0.d0
        kkp_ck2          = 0.d0
        kkp_ck3          = 0.d0
        kkp_ckp1         = 0.d0
        kkp_ckp2         = 0.d0
        kkp_ckp3         = 0.d0
        local_alf        = 0.d0
        local_alfp       = 0.d0
        local_alft       = 0.d0
        local_betp       = 0.d0
        local_c3c        = 0.d0
        local_p3c        = 0.d0
        local_q3c        = 0.d0
        par_bhf          = 0.d0
        par_bhh          = 0.d0
        pt1cm_degsdr     = 0.d0
        pt1cm_engygs     = 0.d0
        pt1cm_r          = 0.d0
        pt3cm_ezero      = 0.d0
        pt4cm_deesdr     = 0.d0
        pt4cm_engyes     = 0.d0
        pt5cm_deijdr     = 0.d0
        pt5cm_engyij     = 0.d0
        ref_ro           = 0.d0
        swpar_dec        = 0.d0
        swpar_dsw        = 0.d0
        swpar_rb         = 0.d0
        swpar_sw         = 0.d0
        tpar_p           = 0.d0
        tpar_xhf         = 0.d0
        tpar_xhh         = 0.d0
        trial_eta1       = 0.d0
        trial_eta2       = 0.d0
        trial_eta3       = 0.d0
        trial_etap1      = 0.d0
        trial_etap2      = 0.d0
        trial_etap3      = 0.d0
        usricm_anuzero   = 0.d0
        usricm_cart      = 0.d0
        usricm_nasurf    = 0
        usricm_nder      = 0
        usricm_nflag     = 0
        usricm_nulbl     = 0
        usrocm_descart   = 0.d0
        usrocm_dgscart   = 0.d0
        usrocm_dijcart   = 0.d0
        usrocm_pengyes   = 0.d0
        usrocm_pengygs   = 0.d0
        usrocm_pengyij   = 0.d0
        utilcm_cnvrtd    = 0.d0
        utilcm_cnvrtde   = 0.d0
        utilcm_cnvrte    = 0.d0
        utilcm_descartnu = 0.d0
        utilcm_dgscartnu = 0.d0
        utilcm_dijcartnu = 0.d0
        utilcm_ireorder  = 0
        utilcm_kediag    = 0
        utilcm_keoffd    = 0
        utilcm_ksdiag    = 0
        utilcm_ksoffd    = 0
        vtanh_drtanh     = 0.d0
        vtanh_rtanh      = 0.d0
        usricm_nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0], shape(usricm_nasurf))
        usricm_nder     = 0
        usricm_nflag    = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero  = 0.0d0
        infocm_icartr   = 3
        infocm_msurf    = 0
        infocm_mder     = 1
        usricm_nulbl    = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms   = 3
        decm_disse      = [0.224989335d0, 0.174472692d0, 0.224989355d0]
        tpar_xhf        = [2.633449d0, 0.417523d0, -0.090289d0]
        tpar_xhh        = [-0.393828d0, 0.167810d0, -0.017081d0]
        tpar_p          = [-0.071547d0, -0.108550d0, 0.050461d0, -0.370255d0]
        local_c3c       = 0.296668d0
        local_alfp      = 0.078606d0
        local_alf       = 0.18d0
        local_alft      = 2.14d0
        local_p3c       = 0.95d0
        local_q3c       = 1.0d0
        local_betp      = 0.18d0
        swpar_dec       = [2.0d0, 2.0d0, 2.0d0]
        swpar_rb        = [6.80d0, 5.00d0, 6.80d0]
        par_bhh         = [0.448467d0, -0.056687d0, 0.246831d0, -0.018419d0, 0.000598d0]
        par_bhf         = [0.9149555d0, 0.3222197d0, 0.1273001d0, 0.0d0, 0.0d0]
        atat_chh(6)     = 6.499027d0
        atat_chh(8)     = 1.243991d2
        atat_chh(10)    = 3.2858d3
        atat_chf(6)     = 6.68799d0
        atat_chf(8)     = 1.0289812d2
        atat_chf(10)    = 2.07391451d3
        dist_rm(1)      = 1.7329d0
        dist_rm(2)      = 1.4010d0
        dist_rm(3)      = 1.7329d0
        dist_r0(1)      = 5.9d0
        dist_r0(2)      = 6.9283d0
        dist_r0(3)      = 5.9d0
        ref_ro(1)       = 1.7329d0
        ref_ro(2)       = 1.449d0
        ref_ro(3)       = 1.7329d0
        dispc_alph0     = 2.59528d1
        dispc_alph1     = 1.1868d0
        dispc_bet0      = 1.57381d1
        dispc_bet1      = 9.729d-2
        trial_eta1(6)   = 0.9438421608d0
        trial_eta1(8)   = 0.9774440533d0
        trial_eta1(10)  = 0.9452240321d0
        trial_eta2(6)   = 0.1085550150d1
        trial_eta2(8)   = 0.1117269260d1
        trial_eta2(10)  = 0.1138536961d1
        trial_eta3(6)   = 0.9438421608d0
        trial_eta3(8)   = 0.9774440533d0
        trial_eta3(10)  = 0.9452240321d0
        trial_etap1(6)  = 6.0d0
        trial_etap1(8)  = 6.0d0
        trial_etap1(10) = 6.0d0
        trial_etap2(6)  = 6.0d0
        trial_etap2(8)  = 6.0d0
        trial_etap2(10) = 6.0d0
        trial_etap3(6)  = 6.0d0
        trial_etap3(8)  = 6.0d0
        trial_etap3(10) = 6.0d0
        kkp_ck1(6)      = -0.2610767389d-1
        kkp_ck1(8)      = -0.3428701964d-1
        kkp_ck1(10)     = -0.4858536634d-1
        kkp_ck2(6)      = -0.9709847270d-1
        kkp_ck2(8)      = -0.1248186655d0
        kkp_ck2(10)     = -0.1426680807d0
        kkp_ck3(6)      = -0.2610767389d-1
        kkp_ck3(8)      = -0.3428701964d-1
        kkp_ck3(10)     = -0.4858536634d-1
        kkp_ckp1(6)     = 0.9438421608d0
        kkp_ckp1(8)     = 0.9774440533d0
        kkp_ckp1(10)    = 0.9452240321d0
        kkp_ckp2(6)     = 0.1085550150d1
        kkp_ckp2(8)     = 0.1117269260d1
        kkp_ckp2(10)    = 0.1138536961d1
        kkp_ckp3(6)     = 0.9438421608d0
        kkp_ckp3(8)     = 0.9774440533d0
        kkp_ckp3(10)    = 0.9452240321d0
        diat_d(1)       = 0.19383609d0
        diat_a1(1)      = 2.4822787d0
        diat_a2(1)      = 1.5435337d0
        diat_a3(1)      = 0.83093855d0
        diat_gamma(1)   = 2.3992999d0
        diat_d(2)       = 0.15796326d0
        diat_a1(2)      = 2.1977034d0
        diat_a2(2)      = 1.2932502d0
        diat_a3(2)      = 0.64375666d0
        diat_gamma(2)   = 2.1835071d0
        diat_d(3)       = 0.19383609d0
        diat_a1(3)      = 2.4822787d0
        diat_a2(3)      = 1.5435337d0
        diat_a3(3)      = 0.83093855d0
        diat_gamma(3)   = 2.3992999d0
        v               = 0.d0; g = 0.d0; d = 0.d0
        usricm_cart     = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                usricm_cart(iatom,idir) = x(iatom,idir)/legacy_bohr_to_angstrom
            enddo
        enddo
        usricm_nder = igrad
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*legacy_hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*legacy_hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = usrocm_dgscart(iatom,idir)*51.422067
                enddo
            enddo
        elseif( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif

        contains

        !> @brief Evaluate prepare_local_context for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine prepare_local_context
            implicit none
            integer                     :: i
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !   Echo the potential name and parameters to the
            !   file linked to UNIT usricm_nflag(18)
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !     CONSTANTS NEEDED FOR THE DAMPING FUNCTION
            !     A=dispc_alph0/FLOAT(N)**dispc_alph1        N=6,8,10
            !     B=dispc_bet0*EXP(-dispc_bet1*FLOAT(N))     N=6,8,10
            dampc_ad1 = dispc_alph0/6.0d0**dispc_alph1
            dampc_ad2 = dispc_alph0/8.0d0**dispc_alph1
            dampc_ad3 = dispc_alph0/10.0d0**dispc_alph1
            dampc_bd1 = dispc_bet0*exp(-dispc_bet1*6.0d0)
            dampc_bd2 = dispc_bet0*exp(-dispc_bet1*8.0d0)
            dampc_bd3 = dispc_bet0*exp(-dispc_bet1*10.0d0)
            !    Initialize the energy in the three asymptotic valleys
            asycm_easyab   = decm_disse(1)
            asycm_easybc   = decm_disse(2)
            asycm_easyac   = decm_disse(3)
            pt3cm_ezero(1) = decm_disse(2)
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            infocm_ref(1)     = 'S. L. Mielke, G. C. Lynch,'
            infocm_ref(2)     = 'diat_d. G. Truhlar, diat_d. W. Schwenke,'
            infocm_ref(3)     = 'Chem. Phys. Lett. 213, 10(1993)'
            infocm_indexes(1) = 9
            infocm_indexes(2) = 1
            infocm_indexes(3) = 1
            infocm_irctnt     = 2
            !      CALL POTINFO
            call ancvrt
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine pot
            !   Surfaces:
            !      ground electronic state
            !   Zero of energy:
            !      The classical potential energy is set equal to zero for the F
            !      infinitely far from the H2 diatomic and pt1cm_r(H2) set equal to the
            !      H2 equilibrium diatomic value.
            !      N.B.: pt1cm_r(1),[pt1cm_r(FH)], must be at least 50 a.u. for the energy to be
            !   Coordinates:
            !      Internal, Definition: pt1cm_r(1) = pt1cm_r(F-H)
            !                            pt1cm_r(2) = pt1cm_r(H-H)
            !                            pt1cm_r(3) = pt1cm_r(H-F)
            !        The energy in the AB valley, asycm_easyab, is equal to the energy
            !        of the C atom "infinitely" far from the AB diatomic and pt1cm_r(AB) s
            !        equal to Re(AB), the equilibrium bond length for the AB diatomi
            !        In this potential the AB valley represents H infinitely far fro
            !        the FH diatomic and pt1cm_r(FH) equal to Re(FH).
            !        Similarly, the terms asycm_easybc and asycm_easyac represent the energies i
            !        H2 and the other HF asymptotic valleys, respectively.
            implicit none
            integer                     :: i
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            call cartou
            call carttor
            !   Check the value of usricm_nder
            !   I: Compute the damping function for the N th. dispersion coefficient
            !      and, if usricm_nder = 1, compute the derivative of the damping function
            !      w.pt1cm_r.t. pt1cm_r(i)
            call damp(pt1cm_r, dist_rm)
            !   II: Compute the two-body dynamical correlation term and, if usricm_nder = 1
            !       compute the derivative of the two-body dynamical correlation ter
            !       w.pt1cm_r.t. pt1cm_r(i)
            call corr2
            !   III: Compute the auxiliary function, hn(pt1cm_r), and, if usricm_nder = 1,
            !            compute the derivative of hn(pt1cm_r) w.pt1cm_r.t. pt1cm_r
            call xtanh(1, trial_eta1)
            call xtanh(2, trial_eta2)
            call xtanh(3, trial_eta3)
            !   IV: Compute the three-body dynamical correlation term and, if usricm_nder =
            !       compute the derivative of the three-body dynamical correlation t
            !       w.pt1cm_r.t. pt1cm_r(i)
            call corr3
            !   V: Compute the switching function, S(pt1cm_r), which is used in the comput
            !      of the H2 triplet state curves.  If usricm_nder = 1, compute the derivat
            !      of S w.pt1cm_r.t. pt1cm_r
            call switch
            !   VI: Compute the lower eigenvalue of the London eq. and, if usricm_nder = 1,
            !       compute the derivatives of the London energy w.pt1cm_r.t. pt1cm_r.
            call vlond
            !   VII: Compute the three-body term and, if usricm_nder = 1, compute the
            !        the derivatives of the three-body energy w.pt1cm_r.t. pt1cm_r.
            call thcent
            !   IX: Compute the total energy which is a sum of the extended-Hartree-
            !       and the dynamical correlation terms.
            pt1cm_engygs = engycm_ecor2b + engycm_ecor3b + engycm_elond + engycm_e3c + pt3cm_ezero(1) - usricm_anuzero
            !   Compute the derivatives of the energy with respect to the coordinate
            if( usricm_nder .eq. 1 )then
                do i = 1, 3
                    pt1cm_degsdr(i) = degsdrcm_dcorr2(i) + degsdrcm_dcorr3(i) + degsdrcm_dldr(i) + degsdrcm_de3c(i)
                enddo
            endif
            call eunitzero
            if( usricm_nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate CORR2 for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine corr2
            !   Compute the two-body dynamical energies and derivatives.
            implicit none
            double precision            :: cor2r1, cor2r2, cor2r3, r1p10, r1p6, r1p8, r2p10, r2p6, r2p8, r3p10, r3p6, r3p8
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            r1p6          = pt1cm_r(1)**6
            r1p8          = pt1cm_r(1)**8
            r1p10         = pt1cm_r(1)**10
            r2p6          = pt1cm_r(2)**6
            r2p8          = pt1cm_r(2)**8
            r2p10         = pt1cm_r(2)**10
            r3p6          = pt1cm_r(3)**6
            r3p8          = pt1cm_r(3)**8
            r3p10         = pt1cm_r(3)**10
            cor2r1        = -dpc_dp1(1)*atat_chf(6)/r1p6 -dpc_dp1(2)*atat_chf(8)/r1p8 -dpc_dp1(3)*atat_chf(10)/r1p10
            cor2r2        = -dpc_dp2(1)*atat_chh(6)/r2p6 -dpc_dp2(2)*atat_chh(8)/r2p8 -dpc_dp2(3)*atat_chh(10)/r2p10
            cor2r3        = -dpc_dp3(1)*atat_chf(6)/r3p6 -dpc_dp3(2)*atat_chf(8)/r3p8 -dpc_dp3(3)*atat_chf(10)/r3p10
            engycm_ecor2b = cor2r1 + cor2r2 + cor2r3
            !   Compute the derivatives w.pt1cm_r.t. pt1cm_r(i) if usricm_nder = 1
            goto_block_100_640: block
                if( usricm_nder .ne. 1 ) exit goto_block_100_640
                degsdrcm_dcorr2(1) = -atat_chf(6)*(dpc_ddp1(1)-6.0d0*dpc_dp1(1)/pt1cm_r(1))/r1p6 -atat_chf(8)* &
                    & (dpc_ddp1(2)-8.0d0*dpc_dp1(2)/pt1cm_r(1))/r1p8 -atat_chf(10)*(dpc_ddp1(3)- 10.0d0*dpc_dp1(3)/ &
                    & pt1cm_r(1))/r1p10
                degsdrcm_dcorr2(2) = -atat_chh(6)*(dpc_ddp2(1)-6.0d0*dpc_dp2(1)/pt1cm_r(2))/r2p6 -atat_chh(8)* &
                    & (dpc_ddp2(2)-8.0d0*dpc_dp2(2)/pt1cm_r(2))/r2p8 -atat_chh(10)*(dpc_ddp2(3)- 10.0d0*dpc_dp2(3)/ &
                    & pt1cm_r(2))/r2p10
                degsdrcm_dcorr2(3) = -atat_chf(6)*(dpc_ddp3(1)-6.0d0*dpc_dp3(1)/pt1cm_r(3))/r3p6 -atat_chf(8)* &
                    & (dpc_ddp3(2)-8.0d0*dpc_dp3(2)/pt1cm_r(3))/r3p8 -atat_chf(10)*(dpc_ddp3(3)- 10.0d0*dpc_dp3(3)/ &
                    & pt1cm_r(3))/r3p10
            end block goto_block_100_640
        end subroutine corr2

        !> @brief Evaluate DAMP for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine damp( rt1, rt2 )
            !   Compute the damping function for the Nth. dispersion coefficient and
            !   derivative w.pt1cm_r.t. pt1cm_r(i)
            implicit none
            double precision            :: dpol1(3), dpol2(3), dpol3(3), dx(3), pol1(3), pol2(3), pol3(3), rt1(3), rt2(3), x(3)
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            x(1)       = 2.0d0*rt1(1)/(rt2(1)+2.5d0*dist_r0(1))
            x(2)       = 2.0d0*rt1(2)/(rt2(2)+2.5d0*dist_r0(2))
            x(3)       = 2.0d0*rt1(3)/(rt2(3)+2.5d0*dist_r0(3))
            pol1(1)    = dampc_ad1*x(1)+dampc_bd1*x(1)**2
            pol2(1)    = dampc_ad2*x(1)+dampc_bd2*x(1)**2
            pol3(1)    = dampc_ad3*x(1)+dampc_bd3*x(1)**2
            pol1(2)    = dampc_ad1*x(2)+dampc_bd1*x(2)**2
            pol2(2)    = dampc_ad2*x(2)+dampc_bd2*x(2)**2
            pol3(2)    = dampc_ad3*x(2)+dampc_bd3*x(2)**2
            pol1(3)    = dampc_ad1*x(3)+dampc_bd1*x(3)**2
            pol2(3)    = dampc_ad2*x(3)+dampc_bd2*x(3)**2
            pol3(3)    = dampc_ad3*x(3)+dampc_bd3*x(3)**2
            dpc_dp1(1) = (1.0d0-exp(-pol1(1)))**6
            dpc_dp1(2) = (1.0d0-exp(-pol2(1)))**8
            dpc_dp1(3) = (1.0d0-exp(-pol3(1)))**10
            dpc_dp2(1) = (1.0d0-exp(-pol1(2)))**6
            dpc_dp2(2) = (1.0d0-exp(-pol2(2)))**8
            dpc_dp2(3) = (1.0d0-exp(-pol3(2)))**10
            dpc_dp3(1) = (1.0d0-exp(-pol1(3)))**6
            dpc_dp3(2) = (1.0d0-exp(-pol2(3)))**8
            dpc_dp3(3) = (1.0d0-exp(-pol3(3)))**10
            if( usricm_nder .eq. 1 )then
                dx(1)       = 2.0d0/(rt2(1)+2.5d0*dist_r0(1))
                dx(2)       = 2.0d0/(rt2(2)+2.5d0*dist_r0(2))
                dx(3)       = 2.0d0/(rt2(3)+2.5d0*dist_r0(3))
                dpol1(1)    = dampc_ad1+2.0d0*dampc_bd1*x(1)
                dpol2(1)    = dampc_ad2+2.0d0*dampc_bd2*x(1)
                dpol3(1)    = dampc_ad3+2.0d0*dampc_bd3*x(1)
                dpol1(2)    = dampc_ad1+2.0d0*dampc_bd1*x(2)
                dpol2(2)    = dampc_ad2+2.0d0*dampc_bd2*x(2)
                dpol3(2)    = dampc_ad3+2.0d0*dampc_bd3*x(2)
                dpol1(3)    = dampc_ad1+2.0d0*dampc_bd1*x(3)
                dpol2(3)    = dampc_ad2+2.0d0*dampc_bd2*x(3)
                dpol3(3)    = dampc_ad3+2.0d0*dampc_bd3*x(3)
                dpc_ddp1(1) = 6.0d0*dpc_dp1(1)*exp(-pol1(1))*dpol1(1)*dx(1)/ (1.0d0-exp(-pol1(1)))
                dpc_ddp1(2) = 8.0d0*dpc_dp1(2)*exp(-pol2(1))*dpol2(1)*dx(1)/ (1.0d0-exp(-pol2(1)))
                dpc_ddp1(3) = 10.0d0*dpc_dp1(3)*exp(-pol3(1))*dpol3(1)*dx(1)/ (1.0d0-exp(-pol3(1)))
                dpc_ddp2(1) = 6.0d0*dpc_dp2(1)*exp(-pol1(2))*dpol1(2)*dx(2)/ (1.0d0-exp(-pol1(2)))
                dpc_ddp2(2) = 8.0d0*dpc_dp2(2)*exp(-pol2(2))*dpol2(2)*dx(2)/ (1.0d0-exp(-pol2(2)))
                dpc_ddp2(3) = 10.0d0*dpc_dp2(3)*exp(-pol3(2))*dpol3(2)*dx(2)/ (1.0d0-exp(-pol3(2)))
                dpc_ddp3(1) = 6.0d0*dpc_dp3(1)*exp(-pol1(3))*dpol1(3)*dx(3)/ (1.0d0-exp(-pol1(3)))
                dpc_ddp3(2) = 8.0d0*dpc_dp3(2)*exp(-pol2(3))*dpol2(3)*dx(3)/ (1.0d0-exp(-pol2(3)))
                dpc_ddp3(3) = 10.0d0*dpc_dp3(3)*exp(-pol3(3))*dpol3(3)*dx(3)/ (1.0d0-exp(-pol3(3)))
            endif
        end subroutine damp

        !> @brief Evaluate CORR3 for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine corr3
            !   Compute the three-body dynamical energies and derivatives.
            implicit none
            double precision            :: co3(3), dc3dr1(3), dc3dr2(3), dc3dr3(3), dg1, dg2, dg3, dh1, dh2, dh3, dt1dr1, dt1dr2, &
                & dt1dr3, dt2dr1, dt2dr2, dt2dr3, dt3dr1, dt3dr2, dt3dr3, fi, g1, g2, g3, h1, h2, h3, r1pi, r2pi, r3pi, t1, t2, t3
            integer                     :: i, jj
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            jj = 0
            do i=6, 10, 2
                jj   = jj+1
                fi   = dble(i)
                r1pi = pt1cm_r(1)**i
                r2pi = pt1cm_r(2)**i
                r3pi = pt1cm_r(3)**i
                !   Compute the auxiliary function gn(pt1cm_r)
                g1 = 1.0d0+kkp_ck1(i)*exp(-kkp_ckp1(i)*(pt1cm_r(1)-ref_ro(1)))
                g2 = 1.0d0+kkp_ck2(i)*exp(-kkp_ckp2(i)*(pt1cm_r(2)-ref_ro(2)))
                g3 = 1.0d0+kkp_ck3(i)*exp(-kkp_ckp3(i)*(pt1cm_r(3)-ref_ro(3)))
                !   Compute the auxiliary function hn(pt1cm_r)
                h1      = vtanh_rtanh(1,i)**trial_etap1(i)
                h2      = vtanh_rtanh(2,i)**trial_etap2(i)
                h3      = vtanh_rtanh(3,i)**trial_etap3(i)
                t1      = atat_chf(i)*(1.0d0-0.5d0*(g2*h3+g3*h2)) *dpc_dp1(jj)/r1pi
                t2      = atat_chh(i)*(1.0d0-0.5d0*(g3*h1+g1*h3)) *dpc_dp2(jj)/r2pi
                t3      = atat_chf(i)*(1.0d0-0.5d0*(g1*h2+g2*h1)) *dpc_dp3(jj)/r3pi
                co3(jj) = t1+t2+t3
                if( usricm_nder .eq. 1 )then
                    !   Compute the partial derivatives w.pt1cm_r.t. pt1cm_r(i)
                    dg1        = -kkp_ckp1(i)*(g1-1.0d0)
                    dg2        = -kkp_ckp2(i)*(g2-1.0d0)
                    dg3        = -kkp_ckp3(i)*(g3-1.0d0)
                    dh1        = trial_etap1(i)*(vtanh_rtanh(1,i)**(trial_etap1(i)-1))*vtanh_drtanh(1,i)
                    dh2        = trial_etap2(i)*(vtanh_rtanh(2,i)**(trial_etap2(i)-1))*vtanh_drtanh(2,i)
                    dh3        = trial_etap3(i)*(vtanh_rtanh(3,i)**(trial_etap3(i)-1))*vtanh_drtanh(3,i)
                    dt1dr1     = t1*(dpc_ddp1(jj)/dpc_dp1(jj)-fi/pt1cm_r(1))
                    dt1dr2     = -0.5d0*atat_chf(i)*dpc_dp1(jj)*(dg2*h3+g3*dh2)/r1pi
                    dt1dr3     = -0.5d0*atat_chf(i)*dpc_dp1(jj)*(g2*dh3+dg3*h2)/r1pi
                    dt2dr1     = -0.5d0*atat_chh(i)*dpc_dp2(jj)*(g3*dh1+dg1*h3)/r2pi
                    dt2dr2     = t2*dpc_ddp2(jj)/dpc_dp2(jj)-t2*fi/pt1cm_r(2)
                    dt2dr3     = -0.5d0*atat_chh(i)*dpc_dp2(jj)*(dg3*h1+g1*dh3)/r2pi
                    dt3dr1     = -0.5d0*atat_chf(i)*dpc_dp3(jj)*(dg1*h2+g2*dh1)/r3pi
                    dt3dr2     = -0.5d0*atat_chf(i)*dpc_dp3(jj)*(g1*dh2+dg2*h1)/r3pi
                    dt3dr3     = t3*(dpc_ddp3(jj)/dpc_dp3(jj)-fi/pt1cm_r(3))
                    dc3dr1(jj) = dt1dr1+dt2dr1+dt3dr1
                    dc3dr2(jj) = dt1dr2+dt2dr2+dt3dr2
                    dc3dr3(jj) = dt1dr3+dt2dr3+dt3dr3
                endif
            enddo
            !   Compute the three-body dynamical correlation term by summing over al
            !   C6, C8, and C10 terms.
            engycm_ecor3b = co3(1)+co3(2)+co3(3)
            if( usricm_nder .eq. 1 )then
                !   Sum all the partial derivatives.
                degsdrcm_dcorr3(1) = dc3dr1(1)+dc3dr1(2)+dc3dr1(3)
                degsdrcm_dcorr3(2) = dc3dr2(1)+dc3dr2(2)+dc3dr2(3)
                degsdrcm_dcorr3(3) = dc3dr3(1)+dc3dr3(2)+dc3dr3(3)
            endif
        end subroutine corr3

        !> @brief Evaluate VLOND for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine vlond
            !   Compute the energy expressed by the lower eignenvalue of the London
            implicit none
            double precision            :: arg, asc1, asc2, csi, csip2, csip3, csip4, darg, dcsi(3), dehf(3), dehft(3), detrip(3,3), &
                & dex(3,3), dfcsi1(3), dfcsi2(3), dq(3,3), dr, dri(3), drp2, dsc, dscale, dsumq(3), dti(3,3), dtscle(3,3), dx(3), &
                & ehf(3), ehft(3), etrip(3), ex(3), exch, fcsi1, fcsi1t, fcsi2, fcsi2t, ff1, ff2, ff3, q(3), r2p2, r2p3, r2p4, rsc, &
                & scale, texch, ti(3), tmp1, tmp2, tscale(3), x(3), xterm, xterm1, xterm2, xterm3
            integer                     :: i, j
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !   Initialize the arrays that will be used to store the partial
            !   derivatives w.pt1cm_r.t. pt1cm_r(i)
            do i = 1, 3
                do j = 1, 3
                    detrip(i,j) = 0.0d0
                    dtscle(i,j) = 0.0d0
                enddo
            enddo
            !   Compute the pseudo-angular coordinate
            csi   = (pt1cm_r(1)-pt1cm_r(3))/pt1cm_r(2)
            csip2 = csi**2
            csip3 = csi**3
            csip4 = csi**4
            !    Compute the F1 and F2 terms
            fcsi1t = 1.0d0+tpar_p(1)*(1.0d0-csip2)+tpar_p(2)*(1.0d0-csip4)
            fcsi2t = 1.0d0+tpar_p(3)*(1.0d0-csip2)+tpar_p(4)*(1.0d0-csip4)
            fcsi1  = fcsi1t**2
            fcsi2  = fcsi2t**2
            if( usricm_nder .eq. 1 )then
                !   Compute the derivatives w.pt1cm_r.t. pt1cm_r of the terms computed to this point
                dcsi(1) = 1.0d0/pt1cm_r(2)
                dcsi(2) = -csi/pt1cm_r(2)
                dcsi(3) = -1.0d0/pt1cm_r(2)
                tmp1    = 2.0d0*fcsi1t*(-2.0d0*tpar_p(1)*csi-4.0d0*tpar_p(2)*csip3)
                tmp2    = 2.0d0*fcsi2t*(-2.0d0*tpar_p(3)*csi-4.0d0*tpar_p(4)*csip3)
                do i=1, 3
                    dfcsi1(i) = tmp1*dcsi(i)
                    dfcsi2(i) = tmp2*dcsi(i)
                enddo
            endif
            r2p2 = pt1cm_r(2)**2
            r2p3 = pt1cm_r(2)**3
            r2p4 = pt1cm_r(2)**4
            !   Compute the extended-Hartree-Fock curves and the derivatives w.pt1cm_r.t.
            !   for the singlet state of H2 and HF.
            !   infocm_ref: A. J. C. Varandas and J. diat_d. Silva, J. Chem. Soc. Faraday Trans.
            !        82, 593 (1986).
            do i=1, 3
                dr     = pt1cm_r(i)-dist_rm(i)
                drp2   = dr**2
                ehf(i) = -diat_d(i)*(1.0d0+diat_a1(i)*dr+diat_a2(i)*drp2+diat_a3(i)*drp2*dr) *exp(-diat_gamma(i)*dr)
                if( usricm_nder .eq. 1 ) &
                    & dehf(i) = -diat_d(i)*(diat_a1(i) + 2.0d0*diat_a2(i)*dr +&
                    &3.0d0*diat_a3(i)*drp2)*exp(-diat_gamma(i)*dr) -&
                    &diat_gamma(i)*ehf(i)
            enddo
            !   Compute the triplet state curve and the derivative of the triplet st
            !   curve w.pt1cm_r.t. pt1cm_r(i) for the H2 diatomic
            !   infocm_ref.: A. J. C. Varandas and J. Brando
            !         Mol. Phys. 45, 857 (1982)
            scale    = 1.0d0+tpar_xhh(1)*pt1cm_r(2)+tpar_xhh(2)*r2p2+tpar_xhh(3)*r2p3
            arg      = par_bhh(2)*pt1cm_r(2)+par_bhh(3)*r2p2+par_bhh(4)*r2p3+par_bhh(5)*r2p4
            ehft(2)  = par_bhh(1)*exp(-arg)/pt1cm_r(2)
            etrip(2) = scale*ehft(2)
            if( usricm_nder .eq. 1 )then
                !   Compute the derivatives of the H2 triplet curve w.pt1cm_r.t. pt1cm_r(i)
                dscale      = tpar_xhh(1)+2.0d0*tpar_xhh(2)*pt1cm_r(2)+3.0d0*tpar_xhh(3)*r2p2
                darg        = par_bhh(2)+2.0d0*par_bhh(3)*pt1cm_r(2)+3.0d0*par_bhh(4)*r2p2+ 4.0d0*par_bhh(5)*r2p3
                dehft(2)    = ehft(2)*(-darg-1.0d0/pt1cm_r(2))
                detrip(2,2) = dscale*ehft(2)+scale*dehft(2)
                detrip(2,1) = 0.0d0
                detrip(2,3) = 0.0d0
            endif
            !   Compute the triplet state curve and the derivatives of the triplet
            !   state curve w.pt1cm_r.t. pt1cm_r(i) for the HF diatomic.
            do i=1, 3, 2
                arg      = par_bhf(2)*pt1cm_r(i)+par_bhf(3)*pt1cm_r(i)**2
                ehft(i)  = par_bhf(1)*exp(-arg)/pt1cm_r(i)
                x(i)     = exp(-tpar_xhf(2)*pt1cm_r(i))
                etrip(i) = tpar_xhf(1)*(fcsi1*x(i)*x(i)+fcsi2*tpar_xhf(3)*x(i))/pt1cm_r(i)
            enddo
            if( usricm_nder .eq. 1 )then
                !   Compute the derivatives of the HF triplet state curve w.pt1cm_r.t. pt1cm_r(i)
                do i = 1, 3, 2
                    darg     = par_bhf(2)+2.0d0*par_bhf(3)*pt1cm_r(i)
                    dehft(i) = ehft(i)*(-1.0d0/pt1cm_r(i)-darg)
                    dx(i)    = -tpar_xhf(2)*x(i)
                    dri(i)   = 1.0d0
                    do j = 1, 3
                        if( j .ne. i )then
                            dx(j)  = 0.0d0
                            dri(j) = 0.0d0
                        endif
                        detrip(i,j) = tpar_xhf(1) * (dfcsi1(j)*x(i)*x(i) + 2.0d0*fcsi1*x(i)*dx(j) + dfcsi2(j)* &
                            & tpar_xhf(3)*x(i) + fcsi2* tpar_xhf(3)*dx(j))/pt1cm_r(i) - etrip(i)*dri(j)/pt1cm_r(i)
                    enddo
                enddo
            endif
            !   Compute the scaling factor for the triplet state curves
            do i = 1, 3
                tscale(i) = 1.0d0
            enddo
            asc1 = 0.09d0
            asc2 = 0.1d0
            rsc  = (6.0d0 - pt1cm_r(2)) / 4.0d0
            dsc  = 1.0d0 - rsc*rsc
            if( pt1cm_r(2) .gt. 2.0d0 .and. pt1cm_r(2) .lt. 10.0d0 )then
                tscale(2) = 1.0d0+asc1*exp(-asc2/dsc)
                if( usricm_nder .eq. 1 ) &
                    & dtscle(2,2) = asc1*exp(-asc2/dsc)* &
                    & (asc2*2.0d0*rsc*0.25d0/(dsc*dsc))
            endif
            !   Compute the triplet state curves and combine the triplet and the
            !   singlet state curves to determine the eigenvalue of the
            !   London equation.
            do i = 1, 3
                !   Compute the triplet state curves
                ti(i) = (1.0d0-swpar_sw(i))*etrip(i)+swpar_sw(i)*ehft(i)
                !   Scale the triplet state curves
                ti(i) = ti(i)*tscale(i)
                !   Compute the coulomb integrals
                q(i) = 0.5d0*(ehf(i)+ti(i))
                !   Compute the exchange integrals
                ex(i) = 0.5d0*(ehf(i)-ti(i))
            enddo
            ff1          = (ex(1) - ex(2)) ** 2
            ff2          = (ex(2) - ex(3)) ** 2
            ff3          = (ex(3) - ex(1)) ** 2
            exch         = sqrt(0.5d0*(ff1+ff2+ff3))
            engycm_elond = q(1) + q(2) + q(3) - exch
            if( usricm_nder .eq. 1 )then
                !   Compute the derivatives of the lowest eignenvalue of the London equa
                !   w.pt1cm_r.t. pt1cm_r(i)
                do i = 1, 3
                    do j = 1, 3
                        dti(i,j)                = swpar_dsw(i,j)*(ehft(i)-etrip(i)) + detrip(i,j)*(1.0d0 - swpar_sw(i))
                        if( j .eq. i ) dti(i,j) = dti(i,j) + swpar_sw(i)*dehft(i)
                        dti(i,j)                = dti(i,j)*tscale(i) + ti(i)*dtscle(i,j)/tscale(i)
                    enddo
                enddo
                do i = 1, 3
                    do j = 1, 3
                        if( j .eq. i )then
                            dq(i,j)  = 0.5d0*(dehf(j)+dti(i,j))
                            dex(i,j) = 0.5d0*(dehf(j)-dti(i,j))
                        else
                            dq(i,j)  = 0.5d0*dti(i,j)
                            dex(i,j) = -dq(i,j)
                        endif
                    enddo
                enddo
                dsumq(1) = dq(1,1)+dq(2,1)+dq(3,1)
                dsumq(2) = dq(1,2)+dq(2,2)+dq(3,2)
                dsumq(3) = dq(1,3)+dq(2,3)+dq(3,3)
                xterm1   = 2.0d0*ex(1)-ex(2)-ex(3)
                xterm2   = 2.0d0*ex(2)-ex(1)-ex(3)
                xterm3   = 2.0d0*ex(3)-ex(1)-ex(2)
                do i = 1, 3
                    texch = exch
                    if( texch .eq. 0.0d0 ) texch = 1.0d0
                    xterm            = 0.5d0*(xterm1*dex(1,i) + xterm2*dex(2,i)+ xterm3*dex(3,i))/texch
                    degsdrcm_dldr(i) = dsumq(i)-xterm
                enddo
            endif
        end subroutine vlond

        !> @brief Evaluate THCENT for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine thcent
            !   Compute the localized three-body term and the derivatives w.pt1cm_r.t. pt1cm_r(i
            implicit none
            double precision            :: costh, dcos1, dcos2, dcos3, ex1, ex2, ex3, f3, omp3c, omq3c, r12, r1i, r22, r32, r3i, rdif, rdif2, rsum, t, t1, t2, t3, t4, t5
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            omp3c = 1.d0 - local_p3c
            omq3c = 1.d0 - local_q3c
            r12   = pt1cm_r(1) * pt1cm_r(1)
            r22   = pt1cm_r(2) * pt1cm_r(2)
            r32   = pt1cm_r(3) * pt1cm_r(3)
            rsum  = pt1cm_r(1) + pt1cm_r(3)
            rdif  = pt1cm_r(1) - pt1cm_r(3)
            rdif2 = rdif * rdif
            t1    = -local_alfp * rsum
            ex1   = local_c3c * exp(t1 * rsum)
            t2    = -local_alf * rdif
            ex2   = omp3c * exp(t2 * rdif)
            t3    = -local_alft * rdif2
            ex3   = local_p3c * exp(t3 * rdif2)
            r1i   = 1.d0 / pt1cm_r(1)
            r3i   = 1.d0 / pt1cm_r(3)
            t4    = r1i*r3i
            !   COSTH IS THE ANGLE BETWEEN R1 AND R3
            costh = 0.5d0 * (r12 + r32 - r22) * t4
            t5    = omq3c * costh
            !   SET MODIFICATION ...
            t = local_q3c + t5 * costh
            !   engycm_e3c CONTAINS THE THREE CENTER CORRECTION TO V
            engycm_e3c = (ex2 + ex3) * t * ex1
            f3         = exp(-local_betp*(pt1cm_r(1)+pt1cm_r(3)-pt1cm_r(2))**2)
            engycm_e3c = engycm_e3c*f3
            goto_block_280_1132: block
                if( usricm_nder .ne. 1 ) exit goto_block_280_1132
                !   COMPUTE DERIVATIVE OF THE 3C TERM
                !   FIRST, DERIVATIVE OF COSTH
                dcos1 = r3i - costh*r1i
                dcos2 = - pt1cm_r(2)*t4
                dcos3 = r1i - costh*r3i
                t4    = ex2 + ex3
                t5    = t4 * t5
                !   NOW, DERIVATIVES OF engycm_e3c
                degsdrcm_de3c(1) = 2.d0*((t2 * ex2 + 2.d0 * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos1) * ex1*f3+ &
                    & 2.0d0* engycm_e3c*(-local_betp*(pt1cm_r(1)+ pt1cm_r(3)-pt1cm_r(2)))
                degsdrcm_de3c(2) = 2.d0 * t5 * dcos2 * ex1*f3+ 2.0d0*engycm_e3c*(local_betp*(pt1cm_r(1)+pt1cm_r(3)- &
                    & pt1cm_r(2)))
                degsdrcm_de3c(3) = 2.d0*((-t2 * ex2 - 2.d0 * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos3) * ex1*f3+ &
                    & 2.0d0* engycm_e3c*(-local_betp*(pt1cm_r(1)+ pt1cm_r(3)-pt1cm_r(2)))
            end block goto_block_280_1132
        end subroutine thcent

        !> @brief Evaluate XTANH for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine xtanh( i, a )
            !   Compute the auxiliary function, hn, and the derivatives of hn w.pt1cm_r.t.
            implicit none
            double precision            :: a(10), ex, rval, sgn, t, th
            integer                     :: i, j
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            do j = 6, 10, 2
                rval = pt1cm_r(i) * a(j)
                if( rval .gt. 0.0d0 )then
                    ex  = exp( -2.0d0 * rval)
                    sgn = 1.0d0
                else
                    ex  = exp( 2.0d0 * rval)
                    sgn = -1.0d0
                endif
                t                                          = 1.0d0 + ex
                th                                         = 1.0d0/t
                if( usricm_nder .eq. 1 ) vtanh_drtanh(i,j) = 4.0d0 * a(j) * ex * th * th
                vtanh_rtanh(i,j)                           = sgn * th * (1.0d0 - ex)
            enddo
        end subroutine xtanh

        !> @brief Evaluate SWITCH for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine switch
            !   Compute the switching function, S, and the derivatives of S w.pt1cm_r.t. pt1cm_r
            implicit none
            double precision            :: ex2u, exu, hysec, hytan, rmrb, rval, sgn, temp1, temp2
            integer                     :: i, j
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !   Initialize the array that will be used to store the derivatives of
            !   the switching function w.pt1cm_r.t. pt1cm_r(i)
            do i = 1, 3
                do j = 1, 3
                    swpar_dsw(i,j) = 0.0d0
                enddo
            enddo
            !   Compute the switching function.
            do i=1, 3
                rmrb = pt1cm_r(i)-swpar_rb(i)
                rval = rmrb*swpar_dec(i)
                if( rval .gt. 0.0d0 )then
                    ex2u = exp(-2.0d0*rval)
                    exu  = exp(-rval)
                    sgn  = 1.0d0
                else
                    ex2u = exp(2.0d0*rval)
                    exu  = exp(rval)
                    sgn  = -1.0d0
                endif
                temp1       = 1.0d0 + ex2u
                hytan       = sgn*(1.0d0 - ex2u)/temp1
                swpar_sw(i) = 0.5d0*(1.0d0 + hytan)
            enddo
            goto_block_100_1239: block
                if( usricm_nder .ne. 1 ) exit goto_block_100_1239
                do i = 1, 3
                    rmrb = pt1cm_r(i)-swpar_rb(i)
                    rval = rmrb*swpar_dec(i)
                    if( rval .gt. 0.0d0 )then
                        exu = exp(-rval)
                        sgn = 1.0d0
                    else
                        exu = exp(rval)
                        sgn = -1.0d0
                    endif
                    temp2          = 1.0d0 + exu*exu
                    hysec          = sgn*2.0d0*exu/temp2
                    swpar_dsw(i,i) = 0.5d0*swpar_dec(i)*hysec*hysec
                enddo
            end block goto_block_100_1239
        end subroutine switch

        !> @brief Evaluate POTINFO for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine potinfo
            implicit none
            integer                     :: i, j, k, kmax
            integer, parameter          :: n3atom = 75
            integer, parameter          :: natom  = 25
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            write(usricm_nflag(18), '(/)')
            do i =1, 5
                write(usricm_nflag(18), '(2x,a75)') infocm_ref(i)
            enddo
            write(usricm_nflag(18), '(/)')
            kmax = 0
            do i = 1, isurf+1
                do j = 1, isurf+1
                    if( usricm_nasurf(i,j) .ne. 0.and.kmax .lt. max(i, j) ) kmax = max(i, j)
                enddo
            enddo
            write(usricm_nflag(18), '(2x,'' MAX. AND ACTUAL NO. OF EXCITED SURFACES: '',I3,5x,I3)') infocm_msurf, kmax-1
            if( kmax-1 .gt. infocm_msurf )then
                write(6, *) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
                stop
            endif
            utilcm_ksdiag = 0
            utilcm_kediag = 0
            do i = 2, isurf+1
                if( usricm_nasurf(i,i) .ne. 0 )then
                    utilcm_kediag = i-1
                    if( utilcm_ksdiag .eq. 0 ) utilcm_ksdiag = i-1
                endif
            enddo
            utilcm_ksoffd = 0
            utilcm_keoffd = 0
            k             = 0
            do i = 1, isurf
                do j = i+1, isurf+1
                    k = k+1
                    if( usricm_nasurf(i,j)+usricm_nasurf(j,i) .ne. 0 )then
                        utilcm_keoffd = k
                        if( utilcm_ksoffd .eq. 0 ) utilcm_ksoffd = k
                    endif
                enddo
            enddo
            write(usricm_nflag(18), '(2x,'' MAX. AND ACTUAL ORDER OF DERIVATIVES:    '',I3,5x,I3)') infocm_mder, usricm_nder
            if( usricm_nder .gt. infocm_mder )then
                write(6, *) ' WRONG INPUT ON ORDER OF DERIVATIVES'
                stop
            endif
            if( usricm_nflag(19) .eq. 1 )then
                write(usricm_nflag(18), '(/)')
                write(usricm_nflag(18), '(2x,''Cartesian coordinates are supplied by'',/, 2x,''the user in the array usricm_cart.'',//)')
                write(usricm_nflag(18), '(2x,''Provide cartesian coordinates in the'',/, 2x,''following order using the array usricm_cart'',' // &
                    '//, 2x,'' usricm_cart(1,1)...usricm_cart(1,3)   => ATOM 1'',/, 2x,'' usricm_cart(2,1)...usricm_cart(2,3)   => AT' // &
                    'OM 2'',/, 2x,'' usricm_cart(3,1)...usricm_cart(3,3)   => ATOM 3'',/, 2x,'' usricm_cart(N,1)...usricm_cart(N,3)  ' // &
                    ' => ATOM N'',/, 2x,''usricm_cart(25,1)...usricm_cart(25,3)  => ATOM 25'',/)')
                write(usricm_nflag(18), '(2x,''If the user wishes to relabel the atoms,'',/, 2x,''set the variable utilcm_ireorder equal to ' // &
                    '1'',/, 2x,''in the PARAMETER statement.  The user'',/, 2x,''must also specify the ne' // &
                    'w labeling'',/, 2x,''scheme.  This is done using the array'',/, 2x,''usricm_nulbl in the fo' // &
                    'llowing manner:'',//, 2x,''usricm_nulbl(i) = j'',/, 2x,''where:  i => old index'',/, 2x,'' ' // &
                    '       j => new index'',//)')
                write(usricm_nflag(18), '(2x,''Cartesian coordinates can be provided to'',/, 2x,''the potential routine in a variety ' // &
                    'of units.'',/, 2x,''The input units will be converted to Bohr'',/, 2x,''based on the' // &
                    ' following values of the usricm_nflag'',/, 2x,''variable:'',//, 2x,''usricm_nflag(1)  =  1  =>  CA' // &
                    'RTESIANS IN BOHR (no'',/, 2x,''                    conversion required)'',/, 2x,''NF' // &
                    'LAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS'',//)')
                write(usricm_nflag(18), '(2x,''The value of the energy and derivatives'',/, 2x,''(if computed) can be reported in a v' // &
                    'ariety'',/, 2x,''units.  A units conversion will take place'',/, 2x,''as specified b' // &
                    'y the following values of the'',/, 2x,''usricm_nflag variable:'',//, 2x,''usricm_nflag(2) = 1 =>  ' // &
                    'ENERGIES REPORTED IN HARTEEE'',/, 2x,''usricm_nflag(2) = 2 =>  ENERGIES REPORTED IN mHARTRE' // &
                    'E'',/, 2x,''usricm_nflag(2) = 3 =>  ENERGIES REPORTED IN eV'',/, 2x,''usricm_nflag(2) = 4 =>  ENER' // &
                    'GIES REPORTED IN kcal/mol'',/, 2x,''usricm_nflag(2) = 5 =>  ENERGIES REPORTED IN cm**-1'',/' // &
                    '/)')
                write(usricm_nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
                    'he'',/, 2x,''usricm_nflag variable:'',//, 2x,''usricm_nflag(1)=1 & usricm_nflag(2)=1 => DERIVATIVES REPOR' // &
                    'TED IN'',/, 2x,''                           HARTEEE/BOHR'',/, 2x,''usricm_nflag(1)=1 & NFLA' // &
                    'G(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/BOHR' // &
                    ''',/, 2x,''usricm_nflag(1)=1 & usricm_nflag(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''            ' // &
                    '               eV/BOHR'',/, 2x,''usricm_nflag(1)=1 & usricm_nflag(2)=4 => DERIVATIVES REPORTED IN' // &
                    ''',/, 2x,''                           kcal/mol/BOHR'',/, 2x,''usricm_nflag(1)=1 & usricm_nflag(2)=' // &
                    '5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm**-1/BOHR'',//)')
                write(usricm_nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
                    'he'',/, 2x,''usricm_nflag variable:'',//, 2x,''usricm_nflag(1)=2 & usricm_nflag(2)=1 => DERIVATIVES REPOR' // &
                    'TED IN'',/, 2x,''                           HARTEEE/ANGSTROM'',/, 2x,''usricm_nflag(1)=2 & ' // &
                    'usricm_nflag(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/' // &
                    'ANGSTROM'',/, 2x,''usricm_nflag(1)=2 & usricm_nflag(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''    ' // &
                    '                       eV/ANGSTROM'',/, 2x,''usricm_nflag(1)=2 & usricm_nflag(2)=4 => DERIVATIVES ' // &
                    'REPORTED IN'',/, 2x,''                           kcal/mol/ANGSTROM'',/, 2x,''usricm_nflag(1' // &
                    ')=2 & usricm_nflag(2)=5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm*' // &
                    '*-1/ANGSTROM'',//)')
            endif
        end subroutine potinfo

        !> @brief Evaluate ANCVRT for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine ancvrt
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi         = 3.141592653589793d0
            double precision, parameter :: clight     = 2.99792458d08
            double precision, parameter :: cmu_0      = 4.0d0*pi*1.0d-07
            double precision, parameter :: cepsilon_0 = 1.0d0/(cmu_0*clight**2)
            double precision, parameter :: ce         = 1.602176462d-19
            double precision, parameter :: cplanck    = 6.62606876d-34
            double precision, parameter :: cm_e       = 9.10938188d-31
            double precision, parameter :: cang       = 1.0d-10
            double precision, parameter :: cavogadro  = 6.02214199d23
            double precision, parameter :: ckcal      = 4.184d10
            double precision, parameter :: htomillh   = 1000.d0
            double precision, parameter :: htoev      = 27.2113834d0
            double precision, parameter :: htokcal    = 627.509470d0
            double precision, parameter :: htowave    = 219474.631d0
            double precision, parameter :: htokj      = 2625.49962d0
            double precision, parameter :: bohr_a     = .5291772083d0
            integer            :: i, ianum(7,32), inc1, inc2, ind, ind2, inu, isave(natom), isum, j, jsave(natom)
            character(len=3)   :: periodic_1(7,32)
            character(len=2)   :: name1(natom)
            character(len=2)   :: name2(natom)
            character(len=1)   :: iblank
            character(len=20)  :: distance
            character(len=20)  :: units
            do i=1, 7
                do j=1, 32
                    ianum(i,j)      = 0
                    periodic_1(i,j) = ' '
                enddo
            enddo
            distance         = 'BOHR                '
            units            = 'HARTREE             '
            ianum(1,1)       = 1
            ianum(1,32)      = 2
            ianum(2,1)       = 3
            ianum(2,2)       = 4
            ianum(2,27)      = 5
            ianum(2,28)      = 6
            ianum(2,29)      = 7
            ianum(2,30)      = 8
            ianum(2,31)      = 9
            ianum(2,32)      = 10
            ianum(3,1)       = 11
            ianum(3,2)       = 12
            ianum(3,27)      = 13
            ianum(3,28)      = 14
            ianum(3,29)      = 15
            ianum(3,30)      = 16
            ianum(3,31)      = 17
            ianum(3,32)      = 18
            ianum(4,1)       = 19
            ianum(4,2)       = 20
            ianum(4,17)      = 21
            ianum(4,18)      = 22
            ianum(4,19)      = 23
            ianum(4,20)      = 24
            ianum(4,21)      = 25
            ianum(4,22)      = 26
            ianum(4,23)      = 27
            ianum(4,24)      = 28
            ianum(4,25)      = 29
            ianum(4,26)      = 30
            ianum(4,27)      = 31
            ianum(4,28)      = 32
            ianum(4,29)      = 33
            ianum(4,30)      = 34
            ianum(4,31)      = 35
            ianum(4,32)      = 36
            ianum(5,1)       = 37
            ianum(5,2)       = 38
            ianum(5,17)      = 39
            ianum(5,18)      = 40
            ianum(5,19)      = 41
            ianum(5,20)      = 42
            ianum(5,21)      = 43
            ianum(5,22)      = 44
            ianum(5,23)      = 45
            ianum(5,24)      = 46
            ianum(5,25)      = 47
            ianum(5,26)      = 48
            ianum(5,27)      = 49
            ianum(5,28)      = 50
            ianum(5,29)      = 51
            ianum(5,30)      = 52
            ianum(5,31)      = 53
            ianum(5,32)      = 54
            ianum(6,1)       = 55
            ianum(6,2)       = 56
            ianum(6,3)       = 57
            ianum(6,4)       = 58
            ianum(6,5)       = 59
            ianum(6,6)       = 60
            ianum(6,7)       = 61
            ianum(6,8)       = 62
            ianum(6,9)       = 63
            ianum(6,10)      = 64
            ianum(6,11)      = 65
            ianum(6,12)      = 66
            ianum(6,13)      = 67
            ianum(6,14)      = 68
            ianum(6,15)      = 69
            ianum(6,16)      = 70
            ianum(6,17)      = 71
            ianum(6,18)      = 72
            ianum(6,19)      = 73
            ianum(6,20)      = 74
            ianum(6,21)      = 75
            ianum(6,22)      = 76
            ianum(6,23)      = 77
            ianum(6,24)      = 78
            ianum(6,25)      = 79
            ianum(6,26)      = 80
            ianum(6,27)      = 81
            ianum(6,28)      = 82
            ianum(6,29)      = 83
            ianum(6,30)      = 84
            ianum(6,31)      = 85
            ianum(6,32)      = 86
            ianum(7,1)       = 87
            ianum(7,2)       = 88
            ianum(7,3)       = 89
            ianum(7,4)       = 90
            ianum(7,5)       = 91
            ianum(7,6)       = 92
            ianum(7,7)       = 93
            ianum(7,8)       = 94
            ianum(7,9)       = 95
            ianum(7,10)      = 96
            ianum(7,11)      = 97
            ianum(7,12)      = 98
            ianum(7,13)      = 99
            ianum(7,14)      = 100
            ianum(7,15)      = 101
            ianum(7,16)      = 102
            ianum(7,17)      = 103
            ianum(7,18)      = 104
            ianum(7,19)      = 105
            ianum(7,20)      = 106
            ianum(7,21)      = 107
            ianum(7,22)      = 108
            ianum(7,23)      = 109
            ianum(7,24)      = 110
            ianum(7,25)      = 111
            ianum(7,26)      = 112
            ianum(7,27)      = 113
            ianum(7,28)      = 114
            ianum(7,29)      = 115
            ianum(7,30)      = 116
            ianum(7,31)      = 117
            ianum(7,32)      = 120
            periodic_1(1,1)  = 'H  '
            periodic_1(1,32) = 'He '
            periodic_1(2,1)  = 'Li '
            periodic_1(2,2)  = 'Be '
            periodic_1(2,27) = 'B  '
            periodic_1(2,28) = 'C  '
            periodic_1(2,29) = 'N  '
            periodic_1(2,30) = 'O  '
            periodic_1(2,31) = 'F  '
            periodic_1(2,32) = 'Ne '
            periodic_1(3,1)  = 'Na '
            periodic_1(3,2)  = 'Mg '
            periodic_1(3,27) = 'Al '
            periodic_1(3,28) = 'Si '
            periodic_1(3,29) = 'P  '
            periodic_1(3,30) = 'S  '
            periodic_1(3,31) = 'Cl '
            periodic_1(3,32) = 'Ar '
            periodic_1(4,1)  = 'K  '
            periodic_1(4,2)  = 'Ca '
            periodic_1(4,17) = 'Sc '
            periodic_1(4,18) = 'Ti '
            periodic_1(4,19) = 'V  '
            periodic_1(4,20) = 'Cr '
            periodic_1(4,21) = 'Mn '
            periodic_1(4,22) = 'Fe '
            periodic_1(4,23) = 'Co '
            periodic_1(4,24) = 'Ni '
            periodic_1(4,25) = 'Cu '
            periodic_1(4,26) = 'Zn '
            periodic_1(4,27) = 'Ga '
            periodic_1(4,28) = 'Ge '
            periodic_1(4,29) = 'As '
            periodic_1(4,30) = 'Se '
            periodic_1(4,31) = 'Br '
            periodic_1(4,32) = 'Kr '
            periodic_1(5,1)  = 'Rb '
            periodic_1(5,2)  = 'Sr '
            periodic_1(5,17) = 'Y  '
            periodic_1(5,18) = 'Zr '
            periodic_1(5,19) = 'Nb '
            periodic_1(5,20) = 'Mo '
            periodic_1(5,21) = 'Tc '
            periodic_1(5,22) = 'Ru '
            periodic_1(5,23) = 'Rh '
            periodic_1(5,24) = 'Pd '
            periodic_1(5,25) = 'Ag '
            periodic_1(5,26) = 'Cd '
            periodic_1(5,27) = 'In '
            periodic_1(5,28) = 'Sn '
            periodic_1(5,29) = 'Sb '
            periodic_1(5,30) = 'Te '
            periodic_1(5,31) = 'I  '
            periodic_1(5,32) = 'Xe '
            periodic_1(5,32) = 'Xe '
            do i=1, infocm_natoms
                isave(i) = 0
                jsave(i) = 0
                name1(i) = '  '
                name2(i) = '  '
            enddo
            iblank = ' '
            do ind=1, infocm_natoms
                do i=1, 7
                    do j=1, 32
                        if( infocm_indexes(ind) .eq. ianum(i,j) )then
                            isave(ind) = i
                            jsave(ind) = j
                        endif
                    enddo
                enddo
            enddo
            do ind=1, infocm_natoms
                ind2 = usricm_nulbl(ind)
                if( ind2 .eq. 0 ) ind2=ind
            enddo
            inc1 = 0
            do ind=1, infocm_irctnt-1
                inc1        = inc1+1
                name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
            enddo
            inc2 = 0
            do ind=infocm_irctnt, infocm_natoms
                inc2        = inc2+1
                name2(inc2) = periodic_1(isave(ind),jsave(ind))(:2)
            enddo
            if( usricm_nflag(1) .eq. 2 ) distance = 'ANGSTROMS           '
            if( usricm_nflag(2) .eq. 2 )then
                units = 'MILLIHARTREE        '
            elseif( usricm_nflag(2) .eq. 3 )then
                units = 'EV                  '
            elseif( usricm_nflag(2) .eq. 4 )then
                units = 'KCAL PER MOLE       '
            elseif( usricm_nflag(2) .eq. 5 )then
                units = 'WAVENUMBERS         '
            elseif( usricm_nflag(2) .eq. 6 )then
                units = 'KILOJOULES PER MOLE '
            endif
            utilcm_cnvrtd  = 1.d0
            utilcm_cnvrte  = 1.d0
            utilcm_cnvrtde = 1.d0
            if( usricm_nflag(1) .eq. 2 ) utilcm_cnvrtd = bohr_a
            if( usricm_nflag(2) .eq. 2 )then
                utilcm_cnvrte = utilcm_cnvrte*htomillh
            elseif( usricm_nflag(2) .eq. 3 )then
                utilcm_cnvrte = utilcm_cnvrte*htoev
            elseif( usricm_nflag(2) .eq. 4 )then
                utilcm_cnvrte = utilcm_cnvrte*htokcal
            elseif( usricm_nflag(2) .eq. 5 )then
                utilcm_cnvrte = utilcm_cnvrte*htowave
            elseif( usricm_nflag(2) .eq. 6 )then
                utilcm_cnvrte = utilcm_cnvrte*htokj
            endif
            utilcm_cnvrtde = utilcm_cnvrte/utilcm_cnvrtd
            isum           = 0
            do inu=1, 25
                isum = isum + usricm_nulbl(inu)
            enddo
            utilcm_ireorder = 0
            if( isum .ne. 0 ) utilcm_ireorder = 1
        end subroutine ancvrt

        !> @brief Evaluate CARTOU for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            implicit none
            integer            :: i, j
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            call convert_cartesian_input(usricm_cart, usricm_nulbl, infocm_natoms, utilcm_cnvrtd, utilcm_ireorder == 1, infocm_cartnu)
        end subroutine cartou

        !> @brief Evaluate CARTTOR for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer            :: coordinate_status
            double precision   :: cmm, cphi, flm, hhd, hym, phi, pi, q, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
            integer            :: i, ind, k, l
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            call cartesian_to_potlib_coordinates(infocm_cartnu, infocm_natoms, infocm_icartr, pt1cm_r, coordinate_status)
            if( coordinate_status /= 0 )then
                write(usricm_nflag(18), '(2x,''WRONG infocm_icartr FOR infocm_cartnu; infocm_icartr ='',i5//)') infocm_icartr
                stop
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine eunitzero
            implicit none
            integer            :: i, j
            character(len=75)  :: ref(5)
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            usrocm_pengygs = pt1cm_engygs * utilcm_cnvrte - usricm_anuzero
            if( utilcm_ksdiag .ne. 0 )then
                do i=utilcm_ksdiag, utilcm_kediag
                    usrocm_pengyes(i) = pt4cm_engyes(i) * utilcm_cnvrte - usricm_anuzero
                enddo
            endif
            if( utilcm_ksoffd .ne. 0 )then
                do j=utilcm_ksoffd, utilcm_keoffd
                    usrocm_pengyij(j) = pt5cm_engyij(j) * utilcm_cnvrte
                enddo
            endif
        end subroutine eunitzero

        !> @brief Evaluate RTOCART for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine rtocart
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            double precision   :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, &
                & termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes(n3atom,isurf), &
                & ygs(n3atom), yij(n3atom,jsurf), za, zb, zc, zd
            integer            :: i, ind, j, j1, j2, k, m1
            if( infocm_icartr .eq. 1 )then
                do i = 1, infocm_natoms
                    ind                   = 3*i-2
                    utilcm_dgscartnu(i,1) = pt1cm_degsdr(ind)
                    utilcm_dgscartnu(i,2) = pt1cm_degsdr(ind+1)
                    utilcm_dgscartnu(i,3) = pt1cm_degsdr(ind+2)
                    if( utilcm_ksdiag .ne. 0 )then
                        do j = utilcm_ksdiag, utilcm_kediag
                            utilcm_descartnu(i,1,j) = pt4cm_deesdr(ind,j)
                            utilcm_descartnu(i,2,j) = pt4cm_deesdr(ind+1,j)
                            utilcm_descartnu(i,3,j) = pt4cm_deesdr(ind+2,j)
                        enddo
                    endif
                    if( utilcm_keoffd .ne. 0 )then
                        do k = utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(i,1,k) = pt5cm_deijdr(ind,k)
                            utilcm_dijcartnu(i,2,k) = pt5cm_deijdr(ind+1,k)
                            utilcm_dijcartnu(i,3,k) = pt5cm_deijdr(ind+2,k)
                        enddo
                    endif
                enddo
            elseif( infocm_icartr .eq. 2 )then
                do i = 1, infocm_natoms
                    utilcm_dgscartnu(i,1) = 0.d0
                    utilcm_dgscartnu(i,2) = 0.d0
                    utilcm_dgscartnu(i,3) = 0.d0
                    if( utilcm_ksdiag .ne. 0 )then
                        do j1=utilcm_ksdiag, utilcm_kediag
                            utilcm_descartnu(i,1,j1) = 0.d0
                            utilcm_descartnu(i,2,j1) = 0.d0
                            utilcm_descartnu(i,3,j1) = 0.d0
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do j2=utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(i,1,j2) = 0.d0
                            utilcm_dijcartnu(i,2,j2) = 0.d0
                            utilcm_dijcartnu(i,3,j2) = 0.d0
                        enddo
                    endif
                    do j = 1, infocm_natoms
                        if( j .lt. i )then
                            m1 = infocm_natoms*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
                            m1 = infocm_natoms*(i-1) - (i*(i-1))/2 + j-i
                        else
                            cycle
                        endif
                        y                     = pt1cm_degsdr(m1)
                        termx                 = (infocm_cartnu(i,1)-infocm_cartnu(j,1))/pt1cm_r(m1)
                        termy                 = (infocm_cartnu(i,2)-infocm_cartnu(j,2))/pt1cm_r(m1)
                        termz                 = (infocm_cartnu(i,3)-infocm_cartnu(j,3))/pt1cm_r(m1)
                        utilcm_dgscartnu(i,1) = utilcm_dgscartnu(i,1) + termx*y
                        utilcm_dgscartnu(i,2) = utilcm_dgscartnu(i,2) + termy*y
                        utilcm_dgscartnu(i,3) = utilcm_dgscartnu(i,3) + termz*y
                        if( utilcm_ksdiag .gt. 0 )then
                            y = pt4cm_deesdr(m1,j1)
                            do j1=utilcm_ksdiag, utilcm_kediag
                                utilcm_descartnu(i,1,j1) = utilcm_descartnu(i,1,j1) + termx*y
                                utilcm_descartnu(i,2,j1) = utilcm_descartnu(i,2,j1) + termy*y
                                utilcm_descartnu(i,3,j1) = utilcm_descartnu(i,3,j1) + termz*y
                            enddo
                        elseif( utilcm_ksoffd .gt. 0 )then
                            do j2=utilcm_ksoffd, utilcm_keoffd
                                y                        = pt5cm_deijdr(m1,j2)
                                utilcm_dijcartnu(i,1,j2) = utilcm_dijcartnu(i,1,j2) + termx*y
                                utilcm_dijcartnu(i,2,j2) = utilcm_dijcartnu(i,2,j2) + termy*y
                                utilcm_dijcartnu(i,3,j2) = utilcm_dijcartnu(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            elseif( infocm_icartr .eq. 3 )then
                do i = 1, infocm_natoms
                    ygs(i) = pt1cm_degsdr(i)/pt1cm_r(i)
                    if( utilcm_ksdiag .ne. 0 )then
                        do j=utilcm_ksdiag, utilcm_kediag
                            yes(i,j) = pt4cm_deesdr(i,j)/pt1cm_r(i)
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do k=utilcm_ksoffd, utilcm_keoffd
                            yij(i,k) = pt5cm_deijdr(i,k)/pt1cm_r(i)
                        enddo
                    endif
                enddo
                do k = 1, 3
                    term12                = infocm_cartnu(1,k)-infocm_cartnu(2,k)
                    term23                = infocm_cartnu(2,k)-infocm_cartnu(3,k)
                    term13                = infocm_cartnu(1,k)-infocm_cartnu(3,k)
                    utilcm_dgscartnu(1,k) = term12*ygs(1) + term13*ygs(3)
                    utilcm_dgscartnu(2,k) = -term12*ygs(1) + term23*ygs(2)
                    utilcm_dgscartnu(3,k) = -term13*ygs(3) - term23*ygs(2)
                    if( utilcm_ksdiag .ne. 0 )then
                        do j1=utilcm_ksdiag, utilcm_kediag
                            utilcm_descartnu(1,k,j1) = term12*yes(1,j1) + term13*yes(3,j1)
                            utilcm_descartnu(2,k,j1) = -term12*yes(1,j1) + term23*yes(2,j1)
                            utilcm_descartnu(3,k,j1) = -term13*yes(3,j1) - term23*yes(2,j1)
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do j2=utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(1,k,j2) = term12*yij(1,j2) + term13*yij(3,j2)
                            utilcm_dijcartnu(2,k,j2) = -term12*yij(1,j2) + term23*yij(2,j2)
                            utilcm_dijcartnu(3,k,j2) = -term13*yij(3,j2) - term23*yij(2,j2)
                        enddo
                    endif
                enddo
            elseif( infocm_icartr .eq. 4 )then
                wh   = 1.007825d0
                wf   = 18.99840d0
                sum  = wh+wf
                eps  = wf/sum
                epsp = wh/sum
                u1   = cos(theta1)
                u2   = cos(theta2)
                u3   = cos(phi)
                ss1  = sin(theta1)
                ss2  = sin(theta2)
                ss3  = sin(phi)
                ya   = 0.0d0
                yb   = 0.0d0
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
                write(usricm_nflag(18), '(2X,'' WRONG infocm_icartr FOR DERIVATIVE; infocm_icartr ='',I5//)') infocm_icartr
                stop
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the FH2_GEN_5SECCOG3_1993 potential.
        subroutine dedcou
            implicit none
            integer            :: i, j, k
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            if( utilcm_ireorder .eq. 1 )then
                do i = 1, infocm_natoms
                    usrocm_dgscart(i,1) = utilcm_dgscartnu(usricm_nulbl(i),1) * utilcm_cnvrtde
                    usrocm_dgscart(i,2) = utilcm_dgscartnu(usricm_nulbl(i),2) * utilcm_cnvrtde
                    usrocm_dgscart(i,3) = utilcm_dgscartnu(usricm_nulbl(i),3) * utilcm_cnvrtde
                    if( utilcm_ksdiag .ne. 0 )then
                        do j=utilcm_ksdiag, utilcm_kediag
                            usrocm_descart(i,1,j) = utilcm_descartnu(usricm_nulbl(i),1,j) * utilcm_cnvrtde
                            usrocm_descart(i,2,j) = utilcm_descartnu(usricm_nulbl(i),2,j) * utilcm_cnvrtde
                            usrocm_descart(i,3,j) = utilcm_descartnu(usricm_nulbl(i),3,j) * utilcm_cnvrtde
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do k=utilcm_ksoffd, utilcm_keoffd
                            usrocm_dijcart(i,1,k) = utilcm_dijcartnu(usricm_nulbl(i),1,k) * utilcm_cnvrtde
                            usrocm_dijcart(i,2,k) = utilcm_dijcartnu(usricm_nulbl(i),2,k) * utilcm_cnvrtde
                            usrocm_dijcart(i,3,k) = utilcm_dijcartnu(usricm_nulbl(i),3,k) * utilcm_cnvrtde
                        enddo
                    endif
                enddo
            else
                do i = 1, infocm_natoms
                    usrocm_dgscart(i,1) = utilcm_dgscartnu(i,1) * utilcm_cnvrtde
                    usrocm_dgscart(i,2) = utilcm_dgscartnu(i,2) * utilcm_cnvrtde
                    usrocm_dgscart(i,3) = utilcm_dgscartnu(i,3) * utilcm_cnvrtde
                    if( utilcm_ksdiag .ne. 0 )then
                        do j=utilcm_ksdiag, utilcm_kediag
                            usrocm_descart(i,1,j) = utilcm_descartnu(i,1,j) * utilcm_cnvrtde
                            usrocm_descart(i,2,j) = utilcm_descartnu(i,2,j) * utilcm_cnvrtde
                            usrocm_descart(i,3,j) = utilcm_descartnu(i,3,j) * utilcm_cnvrtde
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do k=utilcm_ksoffd, utilcm_keoffd
                            usrocm_dijcart(i,1,k) = utilcm_dijcartnu(i,1,k) * utilcm_cnvrtde
                            usrocm_dijcart(i,2,k) = utilcm_dijcartnu(i,2,k) * utilcm_cnvrtde
                            usrocm_dijcart(i,3,k) = utilcm_dijcartnu(i,3,k) * utilcm_cnvrtde
                        enddo
                    endif
                enddo
            endif
        end subroutine dedcou
    end subroutine evaluate
end module chempotpy_fh2_gen_5seccog3_1993_pes

!> @brief Expose FH2_GEN_5SECCOG3_1993 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_fh2_gen_5seccog3_1993_pes, only: evaluate
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
