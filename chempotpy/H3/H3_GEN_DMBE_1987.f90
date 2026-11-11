!   System:          H3
!   Functional Form: Double many-body expansion
!   Common Name:     DMBE
!   Reference:       A. J. C. Varandas, F. B. Brown, C. A. Mead,
!                    D. G. Truhlar, and N. C. Blais
!                    J. Chem. Phys. 86, 6258 (1987).
!   Interface: potlib2001
! Number of bodies: 3
! Number of electronic states: 1
! Number of derivatives: 1
!   Calling Sequence:
!      prepare_local_context - initializes the potential'h3cocm_s variables and
!               must be called once before any calls to POT
!      POT    - driver for the evaluation of the energy and the derivati
!               of the energy with respect to the coordinates for a give
!               geometric configuration
!   Units:
!      energies    - hartrees
!      coordinates - bohrs
!      derivatives - hartrees/bohr
!   Surfaces:
!      ground electronic state
!   Zero of energy:
!      The classical potential energy is set equal to zero for the H
!      infinitely far from the H2 diatomic and R(H2) set equal to the
!      H2 equilibrium diatomic value.
!   Parameters:
!      Set in the BLOCK DATA subprogram PTPACM
!   Coordinates:
!      Internal, Definition: R(1) = R( first H -- second H)
!                            R(2) = R(second H -- third H )
!                            R(3) = R( third H -- first H )
!   Common Blocks (used between the calling program and this potential):
!        passes the coordinates, ground state electronic energy, and
!        derivatives of the ground electronic state energy with respect
!        to the coordinates.
!        passes the control flags where
!              = 1, first excited electronic state energy and derivative
!              = 2, ground and first excited electronic state energies a
!                   derivatives
!        usricm_nder  - not used
!        usricm_nflag - Control flags
!              derivatives for the excited electronic states
!              derivatives for the nonadiabatic coupling between ground
!              excited states.
!   Default Parameter Values:
!      Variable      Default value
!      usricm_nflag(18)        6
!     Calculates double many-body expansion of the H3 potential.
!        ENGYGS  - energy of surface 1 (Ground State).        [ARRAY]
!        DEGSDR  - derivatives for surface 1 (Ground State).  [ARRAY]
!        ENGYES  - energy of surface 2 (Excited State).       [ARRAY]
!        DEESDR  - derivatives for surface 2 (Excited State). [ARRAY]
!        ENGYIJ  - nonadiabatic coupling energy.              [ARRAY]
!        CAPF   - nonadiabatic coupling.  [SUBROUTINE TO CALCULATE ENGYI

!> @brief Implement H3_GEN_DMBE_1987, a potential energy surface for H3.
module chempotpy_h3_gen_dmbe_1987_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_DMBE_1987 potential energy surface for H3.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_geometry,     only: cartesian_distance, flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 2
        integer, parameter            :: natoms  = 3
        integer, parameter            :: n3atom  = 75
        integer, parameter            :: natom   = 25
        integer, parameter            :: isurf   = 5
        integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: cor2cm_corrs
        double precision              :: cor2cm_corrt
        double precision              :: cor2cm_dcor2(4)
        double precision              :: coupcm_dphi(3)
        double precision              :: coupcm_dsphi3(3)
        double precision              :: coupcm_phi
        double precision              :: coupcm_sphi3
        double precision              :: dmpcm_b
        double precision              :: dmpcm_dampb
        double precision              :: dmpcm_dpol
        double precision              :: dmpcm_parb
        double precision              :: dmpcm_texp
        double precision              :: dmpcm_tt
        double precision              :: dmpcm_tt1
        double precision              :: dmpcm_y
        double precision              :: h3cocm_cphi3
        double precision              :: h3cocm_dcphi1
        double precision              :: h3cocm_dcphi2
        double precision              :: h3cocm_dcphi3
        double precision              :: h3cocm_dq1
        double precision              :: h3cocm_dq2
        double precision              :: h3cocm_dq3
        double precision              :: h3cocm_drho1
        double precision              :: h3cocm_drho2
        double precision              :: h3cocm_drho3
        double precision              :: h3cocm_ds1
        double precision              :: h3cocm_ds2
        double precision              :: h3cocm_ds3
        double precision              :: h3cocm_per
        double precision              :: h3cocm_per2
        double precision              :: h3cocm_qcoord
        double precision              :: h3cocm_r12
        double precision              :: h3cocm_r22
        double precision              :: h3cocm_r32
        double precision              :: h3cocm_rho
        double precision              :: h3cocm_s
        double precision              :: h3cocm_s2
        double precision              :: h3crcm_corrs(3)
        double precision              :: h3crcm_corrt(3)
        double precision              :: h3crcm_damp(3,3)
        double precision              :: h3crcm_dcorrs(3)
        double precision              :: h3crcm_dcorrt(3)
        double precision              :: h3crcm_ddamp(3,3)
        double precision              :: h3dmcm_al0
        double precision              :: h3dmcm_al1
        double precision              :: h3dmcm_al2
        double precision              :: h3dmcm_al3
        double precision              :: h3dmcm_alph0
        double precision              :: h3dmcm_alph1
        double precision              :: h3dmcm_alph2
        double precision              :: h3dmcm_alpha5
        double precision              :: h3dmcm_az2
        double precision              :: h3dmcm_bet0
        double precision              :: h3dmcm_bet1
        double precision              :: h3dmcm_beta1
        double precision              :: h3dmcm_beta2
        double precision              :: h3dmcm_beta3
        double precision              :: h3dmcm_cd0
        double precision              :: h3dmcm_cd1
        double precision              :: h3dmcm_chh(3)
        double precision              :: h3dmcm_ck0(3)
        double precision              :: h3dmcm_ck1(3)
        double precision              :: h3dmcm_dampa(3)
        double precision              :: h3dmcm_dampb(3)
        double precision              :: h3dmcm_h2r0
        double precision              :: h3dmcm_h2rm
        double precision              :: h3dmcm_h2rmt
        double precision              :: h3dmcm_h2ta
        double precision              :: h3dmcm_h2tb1
        double precision              :: h3dmcm_h2tb2
        double precision              :: h3dmcm_h2tb3
        double precision              :: h3dmcm_h2tb4
        double precision              :: h3dmcm_hfa1
        double precision              :: h3dmcm_hfa2
        double precision              :: h3dmcm_hfa3
        double precision              :: h3dmcm_hfd
        double precision              :: h3dmcm_hfgam
        double precision              :: h3dmcm_rhol
        double precision              :: h3dmcm_sqrt3
        double precision              :: h3dmcm_xpar(15)
        double precision              :: ijconvert_cnvrtchi
        double precision              :: ijconvert_cnvrtq
        double precision              :: ijconvert_cnvrts
        integer                       :: ijconvert_icouple
        double precision              :: infocm_cartnu(natom,3)
        integer                       :: infocm_icartr
        integer                       :: infocm_indexes(natom)
        integer                       :: infocm_irctnt
        integer                       :: infocm_mder
        integer                       :: infocm_msurf
        integer                       :: infocm_natoms
        character                     :: infocm_ref(5)
        double precision              :: pt1cm_degsdr(n3atom)
        double precision              :: pt1cm_engygs
        double precision              :: pt1cm_r(n3atom)
        double precision              :: pt3cm_ezero(isurf+1)
        double precision              :: pt4cm_deesdr(n3atom,isurf)
        double precision              :: pt4cm_engyes(isurf)
        double precision              :: pt5cm_deijdr(n3atom,jsurf)
        double precision              :: pt5cm_engyij(jsurf)
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
        double precision              :: v
        integer                       :: iatom
        double precision              :: tx(9), r1(3), drdx(3,9)
        double precision              :: hr(nstates,nstates,3), gr(nstates,3)
        double precision              :: hx(nstates,nstates,9), gx(nstates,9)
        integer                       :: istate, idir, i, j, k, l
        !initialize
        cor2cm_corrs     = 0.d0; cor2cm_corrt = 0.d0; cor2cm_dcor2 = 0.d0; coupcm_dphi = 0.d0; coupcm_dsphi3 = 0.d0
        coupcm_phi       = 0.d0; coupcm_sphi3 = 0.d0; dmpcm_b = 0.d0; dmpcm_dampb = 0.d0; dmpcm_dpol = 0.d0; dmpcm_parb = 0.d0
        dmpcm_texp       = 0.d0; dmpcm_tt = 0.d0; dmpcm_tt1 = 0.d0; dmpcm_y = 0.d0; h3cocm_cphi3 = 0.d0; h3cocm_dcphi1 = 0.d0
        h3cocm_dcphi2    = 0.d0; h3cocm_dcphi3 = 0.d0; h3cocm_dq1 = 0.d0; h3cocm_dq2 = 0.d0; h3cocm_dq3 = 0.d0
        h3cocm_drho1     = 0.d0; h3cocm_drho2 = 0.d0; h3cocm_drho3 = 0.d0; h3cocm_ds1 = 0.d0; h3cocm_ds2 = 0.d0
        h3cocm_ds3       = 0.d0; h3cocm_per = 0.d0; h3cocm_per2 = 0.d0; h3cocm_qcoord = 0.d0; h3cocm_r12 = 0.d0
        h3cocm_r22       = 0.d0; h3cocm_r32 = 0.d0; h3cocm_rho = 0.d0; h3cocm_s = 0.d0; h3cocm_s2 = 0.d0; h3crcm_corrs = 0.d0
        h3crcm_corrt     = 0.d0; h3crcm_damp = 0.d0; h3crcm_dcorrs = 0.d0; h3crcm_dcorrt = 0.d0; h3crcm_ddamp = 0.d0
        h3dmcm_al0       = 0.d0; h3dmcm_al1 = 0.d0; h3dmcm_al2 = 0.d0; h3dmcm_al3 = 0.d0; h3dmcm_alph0 = 0.d0
        h3dmcm_alph1     = 0.d0; h3dmcm_alph2 = 0.d0; h3dmcm_alpha5 = 0.d0; h3dmcm_az2 = 0.d0; h3dmcm_bet0 = 0.d0
        h3dmcm_bet1      = 0.d0; h3dmcm_beta1 = 0.d0; h3dmcm_beta2 = 0.d0; h3dmcm_beta3 = 0.d0; h3dmcm_cd0 = 0.d0
        h3dmcm_cd1       = 0.d0; h3dmcm_chh = 0.d0; h3dmcm_ck0 = 0.d0; h3dmcm_ck1 = 0.d0; h3dmcm_dampa = 0.d0
        h3dmcm_dampb     = 0.d0; h3dmcm_h2r0 = 0.d0; h3dmcm_h2rm = 0.d0; h3dmcm_h2rmt = 0.d0; h3dmcm_h2ta = 0.d0
        h3dmcm_h2tb1     = 0.d0; h3dmcm_h2tb2 = 0.d0; h3dmcm_h2tb3 = 0.d0; h3dmcm_h2tb4 = 0.d0; h3dmcm_hfa1 = 0.d0
        h3dmcm_hfa2      = 0.d0; h3dmcm_hfa3 = 0.d0; h3dmcm_hfd = 0.d0; h3dmcm_hfgam = 0.d0; h3dmcm_rhol = 0.d0
        h3dmcm_sqrt3     = 0.d0; h3dmcm_xpar = 0.d0; ijconvert_cnvrtchi = 0.d0; ijconvert_cnvrtq = 0.d0
        ijconvert_cnvrts = 0.d0; ijconvert_icouple = 0; infocm_cartnu = 0.d0; infocm_icartr = 0; infocm_indexes = 0
        infocm_irctnt    = 0; infocm_mder = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref       = ' '
        pt1cm_degsdr     = 0.d0
        pt1cm_engygs     = 0.d0
        pt1cm_r          = 0.d0
        pt3cm_ezero      = 0.d0
        pt4cm_deesdr     = 0.d0
        pt4cm_engyes     = 0.d0
        pt5cm_deijdr     = 0.d0
        pt5cm_engyij     = 0.d0
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
        pt3cm_ezero(1)   = 0.174474112d0
        pt3cm_ezero(2)   = 0.174474112d0
        usricm_nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0], shape(usricm_nasurf))
        usricm_nder    = 0
        usricm_nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero = 0.0d0
        infocm_icartr  = 3
        infocm_msurf   = 1
        infocm_mder    = 1
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 3
        h3dmcm_rhol    = 2.4848d0
        h3dmcm_alph0   = 25.9528d0
        h3dmcm_alph1   = 1.1868d0
        h3dmcm_bet0    = 15.7381d0
        h3dmcm_bet1    = 0.09729d0
        h3dmcm_alpha5  = 8.2433033d-03
        h3dmcm_beta1   = 0.53302897d0
        h3dmcm_beta2   = 0.39156612d-01
        h3dmcm_beta3   = 0.69996945d0
        h3dmcm_alph2   = 4.735364d-1
        h3dmcm_al0     = 0.45024472d+01
        h3dmcm_al1     = -0.62467617d+01
        h3dmcm_al2     = 0.40966542d+01
        h3dmcm_al3     = 0.21813012d+01
        h3dmcm_az2     = 4.453649d-4
        h3dmcm_cd0     = 6.333404d-03
        h3dmcm_cd1     = -1.726839d-03
        h3dmcm_chh     = [6.499027d0, 1.243991d+02, 3285.8d0]
        h3dmcm_ck0     = [-1.372843d-01, -1.638459d-01, -1.973814d-01]
        h3dmcm_ck1     = [1.011204d0, 9.988099d-01, 9.399411d-01]
        h3dmcm_hfd     = 0.15796326d0
        h3dmcm_hfa1    = 2.1977034d0
        h3dmcm_hfa2    = 1.2932502d0
        h3dmcm_hfa3    = 0.64375666d0
        h3dmcm_hfgam   = 2.1835071d0
        h3dmcm_h2rm    = 1.401d0
        h3dmcm_h2r0    = 6.928203d0
        h3dmcm_h2rmt   = 7.82d0
        h3dmcm_h2ta    = 0.448467d0
        h3dmcm_h2tb1   = -0.056687d0
        h3dmcm_h2tb2   = 0.246831d0
        h3dmcm_h2tb3   = -0.018419d0
        h3dmcm_h2tb4   = 0.000598d0
        h3dmcm_sqrt3   = 1.73205080756887d0
        h3dmcm_xpar = [-0.9286579d-02, 0.2811592d-03, -0.4665659d-05, 0.2069800d-07, 0.2903613d+02, -0.2934824d+01, &
            & 0.7181886d0, -0.3753218d0, -0.1114538d+01, 0.2134221d+01, -0.4164343d0, 0.2022584d0, -0.4662687d-01, - &
            & 0.4818623d+02, 0.2988468d0]
        v           = 0.d0; g = 0.d0; d = 0.d0
        usricm_cart = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                usricm_cart(iatom,idir) = x(iatom,idir)/legacy_bohr_to_angstrom
            enddo
        enddo
        if( igrad .eq. 0 )then
            usricm_nder = 0
        elseif( igrad .eq. 1 .or. igrad .eq. 2 )then
            usricm_nder = 1
        endif
        usricm_nasurf = 1
        call prepare_local_context
        call pot
        p(1) = pt1cm_engygs*legacy_hartree_to_electronvolt
        do istate=2, nstates
            i         = istate-1
            p(istate) = pt4cm_engyes(i)*legacy_hartree_to_electronvolt
        enddo
        ! now compute gradients
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is HHH
        r1(1) = cartesian_distance(x(1,:), x(2,:))
        r1(2) = cartesian_distance(x(2,:), x(3,:))
        r1(3) = cartesian_distance(x(1,:), x(3,:))
        call evdrdx(tx, r1, drdx)
        hr = 0.d0
        gr = 0.d0
        do i=1, 3
            hr(1,2,i) = usrocm_pengyij(i)/legacy_bohr_to_angstrom
            hr(2,1,i) = -usrocm_pengyij(i)/legacy_bohr_to_angstrom
        enddo
        do i=1, 3
            gr(1,i) = pt1cm_degsdr(i)*51.422067
            gr(2,i) = pt4cm_deesdr(i,1)*51.422067
        enddo
        call transform_internal_gradient_to_cartesian(hr, drdx, hx)
        call transform_internal_gradient_to_cartesian(gr, drdx, gx)
        do iatom=1, natoms
            do idir=1, 3
                j                 = (iatom-1)*3+idir
                g(:,iatom,idir)   = gx(:,j)
                d(:,:,iatom,idir) = hx(:,:,j)
            enddo
        enddo

        contains

        !> @brief Evaluate EvdRdX for the H3_GEN_DMBE_1987 potential.
        subroutine evdrdx( x, r, drdx )
            use chempotpy_geometry, only: ordered_pair_distance_jacobian
            implicit none
            integer, parameter            :: first_atoms(3)  = [1, 2, 1]
            integer, parameter            :: second_atoms(3) = [2, 3, 3]
            double precision, intent(in)  :: x(9), r(3)
            double precision, intent(out) :: drdx(3,9)
            call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
        end subroutine evdrdx

        !> @brief Evaluate prepare_local_context for the H3_GEN_DMBE_1987 potential.
        subroutine prepare_local_context
            implicit none
            double precision            :: cnvrtd, cnvrte
            integer                     :: i, id, nexp
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !                                                 *
            !  units conversion for anomalous adiabatic       *
            !  coupling pengyij                               *
            !                                                 *
            !     DIMENSION DLEP(3), DLEP2(3), DV2(3), DV3(3), DVA(3), DC3(3)
            !     DIMENSION R(3)
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            ! NON STANDARD COUPLING BETWEEN SURFACES
            ijconvert_icouple  = 1
            usricm_nasurf(1,2) = 0
            usricm_nasurf(2,1) = 0
            cnvrtd             = 1.d0
            if( usricm_nflag(1) .eq. 2 ) cnvrtd = .529177249d0
            cnvrte = 1.d0
            if( usricm_nflag(2) .eq. 2 ) cnvrte = 1000.d0
            if( usricm_nflag(2) .eq. 3 ) cnvrte = 27.211395d0
            if( usricm_nflag(2) .eq. 4 ) cnvrte = 627.5095d0
            if( usricm_nflag(2) .eq. 5 ) cnvrte = 219474.6d0
            if( usricm_nflag(2) .eq. 6 ) cnvrte = 2625.500d0
            ijconvert_cnvrtq   = cnvrte/cnvrtd
            ijconvert_cnvrts   = cnvrte/cnvrtd
            ijconvert_cnvrtchi = cnvrte
            nexp               = 4
            do id = 1, 3
                nexp             = nexp + 2
                h3dmcm_dampa(id) = h3dmcm_alph0/dble(nexp)**h3dmcm_alph1
                h3dmcm_dampb(id) = h3dmcm_bet0*exp(-h3dmcm_bet1*dble(nexp))
            enddo
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            infocm_ref(1)     = 'A. J. C. Varandas, F. B. Brown, C. A. Mead,'
            infocm_ref(2)     = 'D. G. Truhlar, N. C. Blais,'
            infocm_ref(3)     = 'J. Chem. Phys. 86, 6258(1987)'
            infocm_indexes(1) = 1
            infocm_indexes(2) = 1
            infocm_indexes(3) = 1
            infocm_irctnt     = 2
            !      CALL POTINFO
            call ancvrt
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the H3_GEN_DMBE_1987 potential.
        subroutine pot
            !   Surfaces:
            !      ground electronic state
            !   Zero of energy:
            !      The classical potential energy is set equal to zero for the H
            !      infinitely far from the H2 diatomic and pt1cm_r(H2) set equal to the
            !      H2 equilibrium diatomic value.
            !   Coordinates:
            !      Internal, Definition: pt1cm_r(1) = pt1cm_r( first H -- second H)
            !                            pt1cm_r(2) = pt1cm_r(second H -- third H )
            !                            pt1cm_r(3) = pt1cm_r( third H -- first H )
            implicit none
            double precision            :: ce2, ce3, dc3(3), dlep(3), dlep2(3), dv2(3), dv3(3), dva(3), t, va, vc, vii, viii, vleps, &
                & vleps2
            integer                     :: i, indes
            integer, parameter          :: n3atom = 75
            integer, parameter          :: natom  = 25
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !                                                 *
            !  units conversion for anomalous adiabatic       *
            !  coupling usrocm_pengyij                               *
            !                                                 *
            !     DIMENSION pt1cm_r(3)
            call cartou
            call carttor
            ! Set up symmetry coordinates and their derivatives
            call h3coor (pt1cm_r(1), pt1cm_r(2), pt1cm_r(3))
            ! Set up 2-body correlation energies, dispersion damping terms, and
            !    their derivatives
            call h3cor2
            ! Get Leps potentials and derivatives
            call h3leps (vleps, vleps2, dlep, dlep2)
            ! Get VA term and its derivatives
            call h3va (pt1cm_r(1), pt1cm_r(2), pt1cm_r(3), va, dva)
            ! Get VII term and its derivatives
            call h3vii (pt1cm_r(1), pt1cm_r(2), pt1cm_r(3), vii, dv2)
            ! Get VIII term and its derivatives
            call h3viii (viii, dv3)
            ! Get 3-body correlation and its derivative
            call h3cor3 (ce3, dc3)
            ! 2-body correlation term
            ce2 = h3crcm_corrs(1) + h3crcm_corrs(2) + h3crcm_corrs(3)
            ! ground state
            if( usricm_nasurf(1,1) .gt. 0 )then
                vc           = va + vii + viii + ce2 + ce3
                pt1cm_engygs = vleps + vc + pt3cm_ezero(1)
                if( usricm_nder .eq. 1 )then
                    do i = 1, 3
                        t               = dva(i) + dv2(i) + dv3(i) + h3crcm_dcorrs(i) + dc3(i)
                        pt1cm_degsdr(i) = ( dlep(i) + t)
                    enddo
                endif
            endif
            ! first electronic state
            if( usricm_nasurf(2,2) .gt. 0 )then
                indes               = 1
                vc                  = va + vii + viii + ce2 + ce3
                pt4cm_engyes(indes) = vleps2 + vc + pt3cm_ezero(2)
                if( usricm_nder .eq. 1 )then
                    do i = 1, 3
                        t                     = dva(i) + dv2(i) + dv3(i) + h3crcm_dcorrs(i) + dc3(i)
                        pt4cm_deesdr(i,indes) = ( dlep2(i) + t)
                    enddo
                endif
            endif
            ! nonadiabatic coupling in an anomalous form where usrocm_pengyij is treated as
            ! -  a vector with <psi|dpsi/dQ>,<psi|dpsi/ds>,<psi|dpsi/dchi>
            ! -  already transformed to the user'h3cocm_s units
            !    (thus no RTOCART, DEDCOU calls via setting is usricm_nasurf(i,j))
            ! -  having no coded derivatives (no matter what usricm_nder is set to )
            if( ijconvert_icouple .gt. 0 )then
                call capf(pt1cm_r(1), pt1cm_r(2), pt1cm_r(3))
                usrocm_pengyij(1) = pt5cm_engyij(1) * ijconvert_cnvrtq
                usrocm_pengyij(2) = pt5cm_engyij(2) * ijconvert_cnvrts
                usrocm_pengyij(3) = pt5cm_engyij(3) * ijconvert_cnvrtchi
            endif
            call eunitzero
            if( usricm_nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate H3COOR for the H3_GEN_DMBE_1987 potential.
        subroutine h3coor ( r1, r2, r3 )
            !   Calculates D3H symmetry coordinates and derivatives
            implicit none
            double precision            :: arg, bet2, beta, dbt1, dbt2, dbt3, dgm1, dgm2, dgm3, gam2, gamma, phi3, r1, r2, r3, rtan, t, t1, t12, t2, t3, theta
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            h3cocm_per  = r1+r2+r3
            h3cocm_per2 = h3cocm_per*h3cocm_per
            !   h3cocm_qcoord and its derivatives
            h3cocm_r12    = r1*r1
            h3cocm_r22    = r2*r2
            h3cocm_r32    = r3*r3
            h3cocm_qcoord = h3cocm_r12+h3cocm_r22+h3cocm_r32
            h3cocm_dq1    = 2.0d0*r1
            h3cocm_dq2    = 2.0d0*r2
            h3cocm_dq3    = 2.0d0*r3
            !   h3cocm_rho and its derivatives
            h3cocm_rho   = sqrt(h3cocm_qcoord/3.0d0)
            t            = 1.0d0/(6.0d0*h3cocm_rho)
            h3cocm_drho1 = t*h3cocm_dq1
            h3cocm_drho2 = t*h3cocm_dq2
            h3cocm_drho3 = t*h3cocm_dq3
            !   h3cocm_s, h3cocm_cphi3 (cos(phi3)), and their derivatives
            gamma     = 2.0d0*h3cocm_r12-h3cocm_r22-h3cocm_r32
            gam2      = gamma*gamma
            dgm1      = 2.0d0*h3cocm_dq1
            dgm2      = -h3cocm_dq2
            dgm3      = -h3cocm_dq3
            beta      = h3dmcm_sqrt3*(h3cocm_r22-h3cocm_r32)
            bet2      = beta*beta
            dbt1      = 0.0d0
            dbt2      = h3dmcm_sqrt3*h3cocm_dq2
            dbt3      = -h3dmcm_sqrt3*h3cocm_dq3
            t12       = bet2 + gam2
            t1        = sqrt(t12)
            h3cocm_s  = t1/h3cocm_qcoord
            h3cocm_s2 = h3cocm_s*h3cocm_s
            if( h3cocm_s .eq. 0.0d0 )then
                h3cocm_ds1 = 0.0d0
                h3cocm_ds2 = 0.0d0
                h3cocm_ds3 = 0.0d0
                !   For h3cocm_s=0, h3cocm_cphi3 and its derivative should not be
                !   used anywhere but set to zero anyway.
                h3cocm_cphi3  = 0.0d0
                h3cocm_dcphi1 = 0.0d0
                h3cocm_dcphi2 = 0.0d0
                h3cocm_dcphi3 = 0.0d0
            else
                h3cocm_ds1    = h3cocm_s*((beta*dbt1+gamma*dgm1)/t12 - h3cocm_dq1/h3cocm_qcoord)
                h3cocm_ds2    = h3cocm_s*((beta*dbt2+gamma*dgm2)/t12 - h3cocm_dq2/h3cocm_qcoord)
                h3cocm_ds3    = h3cocm_s*((beta*dbt3+gamma*dgm3)/t12 - h3cocm_dq3/h3cocm_qcoord)
                t2            = 1.0d0/(t1*t12)
                h3cocm_cphi3  = gamma*(3.0d0*bet2 - gam2)*t2
                t3            = 3.0d0*beta*(3.0d0*gam2 - bet2)*t2/t12
                h3cocm_dcphi1 = t3*(gamma*dbt1 - beta*dgm1)
                h3cocm_dcphi2 = t3*(gamma*dbt2 - beta*dgm2)
                h3cocm_dcphi3 = t3*(gamma*dbt3 - beta*dgm3)
            endif
            goto_block_10_741: block
                if( usricm_nasurf(2,2) .eq. 0.and.usricm_nflag(4) .eq. 0 ) exit goto_block_10_741
                arg = gamma/t1
                if( arg .le. -1.d0 ) arg=-1.d0
                if( arg .ge. 1.d0 ) arg=1.d0
                theta            = acos(arg)
                coupcm_phi       = pi-theta
                phi3             = 3.0d0*coupcm_phi
                coupcm_sphi3     = sin(phi3)
                rtan             = h3cocm_cphi3/coupcm_sphi3
                coupcm_dsphi3(1) = -rtan*h3cocm_dcphi1
                coupcm_dsphi3(2) = -rtan*h3cocm_dcphi2
                coupcm_dsphi3(3) = -rtan*h3cocm_dcphi3
                coupcm_dphi(1)   = -h3cocm_dcphi1/(coupcm_sphi3*3.d0)
                coupcm_dphi(2)   = -h3cocm_dcphi2/(coupcm_sphi3*3.d0)
                coupcm_dphi(3)   = -h3cocm_dcphi3/(coupcm_sphi3*3.d0)
            end block goto_block_10_741
        end subroutine h3coor

        !> @brief Evaluate H3COR2 for the H3_GEN_DMBE_1987 potential.
        subroutine h3cor2
            !     Calculates 2-body correlation energies for singlet and triplet
            !     states of H2, the damping factors for the dispersion terms, and
            !     their 1st derivatives
            implicit none
            double precision            :: d, dd, denom, dpol, fnexp, pol, rsave, t, t1, t2, tt, tt1, tt2, x
            integer                     :: id, j, nexp
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !      DIMENSION RR(3)
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COR2)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !   Loop over three coordinates
            do j = 1, 3
                rsave            = pt1cm_r(j)
                h3crcm_corrs(j)  = 0.0d0
                h3crcm_dcorrs(j) = 0.0d0
                h3crcm_corrt(j)  = 0.0d0
                h3crcm_dcorrt(j) = 0.0d0
                t1               = 1.0d0/rsave
                t                = t1*t1
                t2               = t*t
                t1               = t2*t1
                nexp             = 4
                !      Loop over terms in dispersion expansion
                do id = 1, 3
                    nexp  = nexp + 2
                    fnexp = dble(nexp)
                    t2    = t2*t
                    t1    = t1*t
                    !     singlet
                    !        Calculate damping function for the nth dispersion coefficient
                    !        and its 1st derivative
                    denom = h3dmcm_h2rm + 2.5d0*h3dmcm_h2r0
                    x     = 2.0d0*rsave/denom
                    tt    = h3dmcm_dampb(id)*x
                    tt1   = h3dmcm_dampa(id) + tt
                    pol   = x*tt1
                    dpol  = tt1 + tt
                    tt    = exp(-pol)
                    tt1   = 1.0d0 - tt
                    tt2   = tt1**(nexp-1)
                    d     = tt1*tt2
                    dd    = 2.0d0*fnexp*dpol*tt*tt2/denom
                    !     store damping factors (including h3dmcm_chh coeff and 1/pt1cm_r**NEXP)
                    !     FOR later use in computing the 3-body correlation energy.
                    h3crcm_damp(id,j)  = h3dmcm_chh(id)*d*t2
                    h3crcm_ddamp(id,j) = h3dmcm_chh(id)*(dd*t2 - fnexp*d*t1)
                    h3crcm_corrs(j)    = h3crcm_corrs(j) - h3crcm_damp(id,j)
                    h3crcm_dcorrs(j)   = h3crcm_dcorrs(j) - h3crcm_ddamp(id,j)
                    !     triplet
                    !        Calculate damping function for the nth dispersion
                    !        coefficient and its 1st derivative
                    denom            = h3dmcm_h2rmt + 2.5d0*h3dmcm_h2r0
                    x                = 2.0d0*rsave/denom
                    tt               = h3dmcm_dampb(id)*x
                    tt1              = h3dmcm_dampa(id) + tt
                    pol              = x*tt1
                    dpol             = tt1 + tt
                    tt               = exp(-pol)
                    tt1              = 1.0d0 - tt
                    tt2              = tt1**(nexp-1)
                    d                = tt1*tt2
                    dd               = 2.0d0*fnexp*dpol*tt*tt2/denom
                    h3crcm_corrt(j)  = h3crcm_corrt(j) - h3dmcm_chh(id)*d*t2
                    h3crcm_dcorrt(j) = h3crcm_dcorrt(j) - h3dmcm_chh(id)*(dd*t2 - fnexp*d*t1)
                enddo
            enddo
        end subroutine h3cor2

        !> @brief Evaluate H3COR3 for the H3_GEN_DMBE_1987 potential.
        subroutine h3cor3 ( ce3, dc3 )
            !     Calculates 3-body correlation energy and its 1st derivatives
            implicit none
            double precision            :: ce3, dc3(3), g(3), gd(3), h(3), hd(3), hysec, hytan, rr, sgnt, t, t1, t2
            integer                     :: i, i2, i3, id
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            ce3    = 0.0d0
            dc3(1) = 0.0d0
            dc3(2) = 0.0d0
            dc3(3) = 0.0d0
            !   Loop over terms in dispersion expansion; dispersion damping
            !   FACTORS are computed in H3COR2 and passed through COMMON /H3CRCM/.
            do id = 1, 3
                !   Loop over 3 coordinates
                do i = 1, 3
                    rr = pt1cm_r(i)
                    !     Calculate G function and its 1st derivative
                    t     = h3dmcm_ck0(id)*exp(-h3dmcm_ck1(id)*(rr-h3dmcm_h2rm))
                    g(i)  = 1.0d0 + t
                    gd(i) = -h3dmcm_ck1(id)*t
                    !     Calculate H function and its 1st derivative.
                    t    = h3dmcm_ck1(id)*rr
                    sgnt = 1.0d0
                    if( t .lt. 0.0d0 )then
                        t    = -t
                        sgnt = -1.0d0
                    endif
                    t     = exp(-t)
                    t2    = t*t
                    t1    = 1.0d0/(1.0d0+t2)
                    hysec = 2.0d0*t*t1
                    hytan = sgnt*(1.0d0-t2)*t1
                    t1    = hytan**5
                    h(i)  = hytan*t1
                    hd(i) = 6.0d0*h3dmcm_ck1(id)*t1*hysec*hysec
                enddo
                do i = 1, 3
                    i2      = mod(i, 3) + 1
                    i3      = mod(i+1, 3) + 1
                    t       = 1.0d0 - 0.5d0*(g(i2)*h(i3) + g(i3)*h(i2))
                    ce3     = ce3 + t*h3crcm_damp(id,i)
                    dc3(i)  = dc3(i) + t*h3crcm_ddamp(id,i)
                    dc3(i2) = dc3(i2) - 0.5d0*(gd(i2)*h(i3)+g(i3)*hd(i2))*h3crcm_damp(id,i)
                    dc3(i3) = dc3(i3) - 0.5d0*(g(i2)*hd(i3)+gd(i3)*h(i2))*h3crcm_damp(id,i)
                enddo
            enddo
        end subroutine h3cor3

        !> @brief Evaluate H3LEPS for the H3_GEN_DMBE_1987 potential.
        subroutine h3leps ( vleps, vleps2, dlep, dlep2 )
            !     Calculates 3-body extended-Hartree-Fock energy defined by a
            !     LEPS-type function.
            !   VLEPS  IS LEPS LOWER SURFACE.
            !   VLEPS2 IS LEPS UPPER SURFACE.
            !   DLEP(3) ARE LEPS LOWER SURFACE DERIVATIVES
            !   DLEP2(3) "   "   UPPER    "       "
            implicit none
            double precision            :: at, dat, desdr, df(3), dlep(3), dlep2(3), dq(3), dr, dwe, dxj(3,3), es, et, exch, f, f1, &
                & f2, f3, fterm, omf, qsum, rr, rri, t, t1, t2, tdc, tdq, tds, texch, vleps, vleps2, we, xj(3), xterm, xterm1, &
                & xterm2, xterm3
            integer                     :: i, j
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !  Calculate switching term for the Hartree-Fock component of the
            !  diatomic triplet state function.  The notation of Thompson et al.
            !  (J. Chem. Phys., 82,5597,1985; and references therein) is used
            !  throughout.
            t1 = -h3dmcm_az2*h3cocm_qcoord
            if( h3cocm_s .eq. 0.0d0 )then
                f     = exp(t1*h3cocm_qcoord)
                t1    = 2.0d0*t1*f
                df(1) = t1*h3cocm_dq1
                df(2) = t1*h3cocm_dq2
                df(3) = t1*h3cocm_dq3
            else
                t2    = 1.0d0 + h3cocm_s*h3cocm_s2*h3cocm_cphi3
                f     = exp(t1*h3cocm_qcoord*t2)
                t1    = t1*f
                tdq   = 2.0d0*t2
                t2    = h3cocm_s2*h3cocm_qcoord
                tds   = 3.0d0*h3cocm_cphi3*t2
                tdc   = h3cocm_s*t2
                df(1) = t1*(tdq*h3cocm_dq1 + tds*h3cocm_ds1 + tdc*h3cocm_dcphi1)
                df(2) = t1*(tdq*h3cocm_dq2 + tds*h3cocm_ds2 + tdc*h3cocm_dcphi2)
                df(3) = t1*(tdq*h3cocm_dq3 + tds*h3cocm_ds3 + tdc*h3cocm_dcphi3)
            endif
            omf  = 1.0d0 - f
            qsum = 0.0d0
            do i = 1, 3
                dq(i) = 0.0d0
            enddo
            !  Loop over 3 coordinates
            do i = 1, 3
                rr  = pt1cm_r(i)
                rri = 1.0d0/rr
                !     Calculate extended-Hartree-Fock curve for ground singlet state of
                !     H2 [A.J.C. Varandas & J.D. Silva, J. Chem. Soc. Faraday II
                !     (submitted)] and 1st derivative of 2-body extended-Hartree-Fock
                !     curve for the ground-singlet state of H2.
                dr    = rr - h3dmcm_h2rm
                t     = exp(-h3dmcm_hfgam*dr)
                es    = -h3dmcm_hfd*(1.0d0 + dr*(h3dmcm_hfa1 + dr*(h3dmcm_hfa2 + dr*h3dmcm_hfa3)))*t
                desdr = -h3dmcm_hfgam*es - h3dmcm_hfd*(h3dmcm_hfa1 + dr*(2.0d0*h3dmcm_hfa2 + 3.0d0*dr*h3dmcm_hfa3))*t
                !     Compute the HFACE (Hartree-Fock-approximate correlation energy)
                !     potential for the lowest triplet state of H2 [A.J.C. Varandas &
                !     J. Brandao Mol. Phys.,45,1982,857] without the 2-body correlation
                !     energy.
                t   = rr*(h3dmcm_h2tb1 + rr*(h3dmcm_h2tb2 + rr*(h3dmcm_h2tb3 + rr*h3dmcm_h2tb4)))
                at  = h3dmcm_h2ta*exp(-t)*rri
                t   = h3dmcm_h2tb1 + rr*(2.0d0*h3dmcm_h2tb2 + rr*(3.0d0*h3dmcm_h2tb3 + rr*4.0d0*h3dmcm_h2tb4))
                dat = -at*(t + rri)
                !     Add in triplet and subtract singlet 2-body correlation terms
                at  = at + h3crcm_corrt(i) - h3crcm_corrs(i)
                dat = dat + h3crcm_dcorrt(i) - h3crcm_dcorrs(i)
                !     Calculate effective diatomic triplet state curve.
                t   = exp(-h3dmcm_al3*rr)
                we  = (h3dmcm_al0 + rr*(h3dmcm_al1 + rr*h3dmcm_al2))*t
                dwe = (h3dmcm_al1 + 2.0d0*h3dmcm_al2*rr)*t - h3dmcm_al3*we
                !     Triplet energy
                et = f*we + omf*at
                !     Contribution to coulomb term
                qsum = qsum + 0.5d0*(es + et)
                !     Contribution to exchange term
                xj(i) = 0.5d0*(es - et)
                !     Derivatives of coulomb and exchange parts
                fterm = we - at
                do j = 1, 3
                    t        = 0.5d0*fterm*df(j)
                    dq(j)    = dq(j) + t
                    dxj(i,j) = -t
                enddo
                t        = f*dwe + omf*dat
                dq(i)    = dq(i) + 0.5d0*(desdr + t)
                dxj(i,i) = dxj(i,i) + 0.5d0*(desdr - t)
            enddo
            f1     = (xj(1) - xj(2))**2
            f2     = (xj(2) - xj(3))**2
            f3     = (xj(3) - xj(1))**2
            exch   = sqrt(0.5d0*(f1 + f2 + f3))
            vleps  = qsum - exch
            vleps2 = qsum + exch
            xterm1 = 2.0d0*xj(1) - xj(2) - xj(3)
            xterm2 = 2.0d0*xj(2) - xj(1) - xj(3)
            xterm3 = 2.0d0*xj(3) - xj(1) - xj(2)
            do i = 1, 3
                texch = exch
                if( texch .le. 0.0d0 ) texch = 1.0d0
                xterm    = 0.5d0/texch*(xterm1*dxj(1,i)  + xterm2*dxj(2,i) + xterm3*dxj(3,i))
                dlep(i)  = dq(i) - xterm
                dlep2(i) = dq(i) + xterm
            enddo
        end subroutine h3leps

        !> @brief Evaluate H3VA for the H3_GEN_DMBE_1987 potential.
        subroutine h3va ( r1, r2, r3, va, dva )
            !     Calculates Va correction energy and its 1st derivative.
            implicit none
            double precision            :: dv, dva(3), expva, r1, r2, r3, t, t1, t2, t3, t4, va
            integer                     :: j
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
            !     +                CORRT(3),DCORRT(3),
            !     +                DAMP(3,3),DDAMP(3,3)
            t     = h3dmcm_alpha5*h3cocm_per2
            expva = exp(-t*h3cocm_per)
            if( expva .eq. 0.0d0 )then
                va     = 0.0d0
                dva(1) = 0.0d0
                dva(2) = 0.0d0
                dva(3) = 0.0d0
            else
                va = 0.0d0
                dv = 0.0d0
                t3 = (r1-r2)*(r2-r3)*(r3-r1)
                t1 = t3*t3
                t2 = 1.0d0
                t4 = 2.0d0
                do j=1, 4
                    t2 = t2*t1
                    va = va + h3dmcm_xpar(j)*t2
                    dv = dv + t4*h3dmcm_xpar(j)*t3
                    t3 = t3*t1
                    t4 = t4 + 2.0d0
                enddo
                va     = va*expva
                t1     = dv*expva
                t2     = 3.0d0*t*va
                dva(1) = t1*(r2-r3)*(r3+r2-2.0d0*r1) - t2
                dva(2) = t1*(r3-r1)*(r1+r3-2.0d0*r2) - t2
                dva(3) = t1*(r1-r2)*(r1+r2-2.0d0*r3) - t2
            endif
        end subroutine h3va

        !> @brief Evaluate H3VII for the H3_GEN_DMBE_1987 potential.
        subroutine h3vii ( r1, r2, r3, e, dv2 )
            !     Calculates VII correction energy and its 1st derivative.
            implicit none
            double precision            :: cos1, cos2, cos3, dr12, dr23, dr31, dv2(3), e, eb1, eb1p1, eb1p2, eb1p3, eb1pr, eb1pwb, &
                & eb1t, eb2, eb2p1, eb2p2, eb2p3, eb2pr, eb2pwb, eb3t, eb4, eb4a, eb4ap1, eb4ap2, eb4ap3, eb4apr, eb4apw, eb4b, &
                & eb4bp1, eb4bp2, eb4bp3, eb4bpr, eb4bpw, eb4p1, eb4p2, eb4p3, eq, eqp1, eqp2, eqp3, exp1, exp2, exp3, r1, r1i, r2, &
                & r2i, r3, r3i, ri, t1, t2, t3, wb, wb1p, wb2, wb2p, wb3p
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
            !     +                CORRT(3),DCORRT(3),
            !     +                DAMP(3,3),DDAMP(3,3)
            !   Compute B1 function.
            r1i  = 1.0d0/r1
            r2i  = 1.0d0/r2
            r3i  = 1.0d0/r3
            cos1 = 0.5d0*r2i*r3i*(h3cocm_r12-h3cocm_r22-h3cocm_r32)
            cos2 = 0.5d0*r1i*r3i*(h3cocm_r22-h3cocm_r12-h3cocm_r32)
            cos3 = 0.5d0*r1i*r2i*(h3cocm_r32-h3cocm_r12-h3cocm_r22)
            wb   = 1.0d0+cos1+cos2+cos3
            wb2  = wb*wb
            !   WB derivatives
            wb1p = (r1*r3i - 1.0d0)*r2i - r3i - (cos2+cos3)*r1i
            wb2p = (r2*r3i - 1.0d0)*r1i - r3i - (cos1+cos3)*r2i
            wb3p = (r3*r2i - 1.0d0)*r1i - r2i - (cos1+cos2)*r3i
            !   EB1 term
            exp1 = exp(-h3dmcm_beta1*h3cocm_per)
            exp3 = exp(-h3dmcm_beta3*h3cocm_per)
            eb1t = (h3dmcm_xpar(5) + h3dmcm_xpar(6)*h3cocm_per)*exp1
            eb3t = (h3dmcm_xpar(14) + h3dmcm_xpar(15)*h3cocm_per2)*exp3
            eb1  = wb*(eb1t + eb3t)
            !   EB1 derivatives
            eb1pr  = wb*(-h3dmcm_beta1*eb1t - h3dmcm_beta3*eb3t + h3dmcm_xpar(6)*exp1 + 2.0d0*h3cocm_per* &
                & h3dmcm_xpar(15)*exp3)
            eb1pwb = eb1t + eb3t
            eb1p1  = eb1pwb*wb1p + eb1pr
            eb1p2  = eb1pwb*wb2p + eb1pr
            eb1p3  = eb1pwb*wb3p + eb1pr
            !   EB2 term
            t1   = h3dmcm_beta2*h3cocm_per
            exp2 = exp(-t1*h3cocm_per)
            eb2  = wb2*(h3dmcm_xpar(7) + wb*(h3dmcm_xpar(8) + wb*h3dmcm_xpar(9)))*exp2
            !   EB2 derivatives
            eb2pwb = wb*(2.0d0*h3dmcm_xpar(7) + wb*(3.0d0*h3dmcm_xpar(8) + wb*4.0d0*h3dmcm_xpar(9)) )*exp2
            eb2pr  = -2.0d0*t1*eb2
            eb2p1  = eb2pwb*wb1p + eb2pr
            eb2p2  = eb2pwb*wb2p + eb2pr
            eb2p3  = eb2pwb*wb3p + eb2pr
            !   EB4 term
            !      EB4A
            t2   = h3dmcm_xpar(10)*exp1
            t3   = wb*h3dmcm_xpar(11)*exp2
            eb4a = wb*(t2 + t3)
            !      EB4A derivatives
            eb4apw = t2 + 2.0d0*t3
            eb4apr = -wb*(h3dmcm_beta1*t2 + 2.0d0*t1*t3)
            eb4ap1 = eb4apw*wb1p + eb4apr
            eb4ap2 = eb4apw*wb2p + eb4apr
            eb4ap3 = eb4apw*wb3p + eb4apr
            !      EB4B
            t2   = h3dmcm_xpar(12)*exp1
            t3   = h3dmcm_xpar(13)*exp2
            eb4b = wb*(t2 + t3)
            !      EB4B derivatives
            eb4bpw = t2 + t3
            eb4bpr = -wb*(h3dmcm_beta1*t2 + 2.0d0*t1*t3)
            eb4bp1 = eb4bpw*wb1p + eb4bpr
            eb4bp2 = eb4bpw*wb2p + eb4bpr
            eb4bp3 = eb4bpw*wb3p + eb4bpr
            dr12   = r1 - r2
            dr23   = r2 - r3
            dr31   = r3 - r1
            eq     = dr12*dr12 + dr23*dr23 + dr31*dr31
            eqp1   = 2.0d0*(dr12 - dr31)
            eqp2   = 2.0d0*(-dr12 + dr23)
            eqp3   = 2.0d0*(-dr23 + dr31)
            ri     = r1i + r2i + r3i
            eb4    = eb4a*ri + eb4b*eq
            !   EB4 derivatives
            eb4p1  = eb4ap1*ri - eb4a/h3cocm_r12 + eb4bp1*eq + eb4b*eqp1
            eb4p2  = eb4ap2*ri - eb4a/h3cocm_r22 + eb4bp2*eq + eb4b*eqp2
            eb4p3  = eb4ap3*ri - eb4a/h3cocm_r32 + eb4bp3*eq + eb4b*eqp3
            e      = eb1 + eb2 + eb4
            dv2(1) = eb1p1 + eb2p1 + eb4p1
            dv2(2) = eb1p2 + eb2p2 + eb4p2
            dv2(3) = eb1p3 + eb2p3 + eb4p3
        end subroutine h3vii

        !> @brief Evaluate H3VIII for the H3_GEN_DMBE_1987 potential.
        subroutine h3viii ( viii, dv3 )
            !     Calculates VIII correction energy and it 1st derivatives
            implicit none
            double precision            :: dv3(3), expv3, t1, t2, t3, t4, t5, tdcph, tdrho, tds, viii
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
            !     +                CORRT(3),DCORRT(3),
            !     +                DAMP(3,3),DDAMP(3,3)
            if( h3cocm_s .eq. 0.0d0 )then
                viii   = 0.0d0
                dv3(1) = 0.0d0
                dv3(2) = 0.0d0
                dv3(3) = 0.0d0
            else
                t1    = h3cocm_rho - h3dmcm_rhol
                t5    = h3dmcm_alph2*t1
                expv3 = exp(-t5*t1)
                if( expv3 .eq. 0.0d0 )then
                    viii   = 0.0d0
                    dv3(1) = 0.0d0
                    dv3(2) = 0.0d0
                    dv3(3) = 0.0d0
                else
                    t1     = h3cocm_s2*h3cocm_cphi3
                    t2     = 1.0d0 + h3cocm_s*t1
                    t3     = h3cocm_s2*t2
                    t4     = h3dmcm_cd0 + h3dmcm_cd1*h3cocm_rho
                    viii   = t3*t4*expv3
                    tds    = (2.0d0*h3cocm_s*t2 + 3.0d0*h3cocm_s2*t1)*t4
                    tdcph  = h3cocm_s2*h3cocm_s*h3cocm_s2*t4
                    tdrho  = t3*(h3dmcm_cd1 - 2.0d0*t5*t4)
                    dv3(1) = (tds*h3cocm_ds1 + tdcph*h3cocm_dcphi1 + tdrho*h3cocm_drho1)*expv3
                    dv3(2) = (tds*h3cocm_ds2 + tdcph*h3cocm_dcphi2 + tdrho*h3cocm_drho2)*expv3
                    dv3(3) = (tds*h3cocm_ds3 + tdcph*h3cocm_dcphi3 + tdrho*h3cocm_drho3)*expv3
                endif
            endif
        end subroutine h3viii

        !> @brief Evaluate CAPF for the H3_GEN_DMBE_1987 potential.
        subroutine capf( r1, r2, r3 )
            implicit none
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            double precision            :: arg, brak, capk, darg, dcapk(4), dcphi(3), denom, dfu(3), dfz(3), dgz(3), ds(3), dx(3), &
                & f(3,jsurf), fu, fz, gz, r1, r2, r3, rcapf(3), x, x2, xi, xi3
            integer                     :: i
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
            !     +                CORRT(3),DCORRT(3),
            !     +                DAMP(3,3),DDAMP(3,3)
            rcapf(1) = r1
            rcapf(2) = r2
            rcapf(3) = r3
            call h3coor(rcapf(1), rcapf(2), rcapf(3))
            ds(1)    = h3cocm_ds1
            ds(2)    = h3cocm_ds2
            ds(3)    = h3cocm_ds3
            dcphi(1) = h3cocm_dcphi1
            dcphi(2) = h3cocm_dcphi2
            dcphi(3) = h3cocm_dcphi3
            x        = sqrt(h3cocm_qcoord/3.0d0)
            xi       = 1.0d0/x
            xi3      = x/3.0d0
            x2       = x*x
            call kappa(x, capk, dcapk)
            fz    = 0.375d0*x*dcapk(1)
            fu    = 0.03515625d0*(x*dcapk(1)-x2*dcapk(2)+x2*xi3*dcapk(3))
            gz    = 0.046875d0*(-x*dcapk(1)+x2*dcapk(2))
            denom = fz+gz*h3cocm_cphi3*h3cocm_s+fu*h3cocm_s2
            arg   = gz*coupcm_sphi3*h3cocm_s/denom
            !     FCT=coupcm_phi-ATAN(ARG)
            do i=1, 3
                dx(i)           = rcapf(i)*xi
                dfz(i)          = 0.125d0*dx(i)*(dcapk(1)+x*dcapk(2))
                dfu(i)          = 0.01171875d0*dx(i)*(dcapk(1)-x*dcapk(2)+x2*xi3*dcapk(4))
                dgz(i)          = 0.015625d0*dx(i)*(-dcapk(1)+x*dcapk(2)+x2*dcapk(3))
                brak            = h3cocm_s*dgz(i)+gz*ds(i)
                darg            = -arg*(dfz(i)+h3cocm_cphi3*brak+h3cocm_s*(gz*dcphi(i)+h3cocm_s*dfu(i)+ 2.0d0*fu*ds(i)))
                darg            = (darg+gz*h3cocm_s*coupcm_dsphi3(i)+coupcm_sphi3*brak)/denom
                pt5cm_engyij(i) = coupcm_dphi(i)-1.0d0/(1.0d0+arg**2)*darg
            enddo
        end subroutine capf

        !> @brief Evaluate H2COR2 for the H3_GEN_DMBE_1987 potential.
        subroutine h2cor2 ( rr )
            !     CALCULATES 1-4 DERIVATIVES OF 2-BODY CORRELATION TERMS FOR
            !     BOTH SINGLET AND TRIPLET STATES AND RETURNS THEIR DIFFERENCE
            implicit none
            double precision            :: d, d2, d3, d4, dcorrs, dcorrt, dd, denom, fnexp, pol, rr, rsave, t, t1, t2, tem2, tem3, tem4, tt2, x
            integer                     :: id, j, nexp
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !      COMMON /H3CRCM/ cor2cm_corrs(3),DCORRS(3),
            !     +                cor2cm_corrt(3),DCORRT(3),
            !     +                DAMP(3,3),DDAMP(3,3)
            rsave        = rr
            cor2cm_corrs = 0.0d0
            dcorrs       = 0.0d0
            cor2cm_corrt = 0.0d0
            dcorrt       = 0.0d0
            do j=1, 4
                cor2cm_dcor2(j) = 0.0d0
            enddo
            t1   = 1.0d0/rsave
            t    = t1*t1
            t2   = t*t
            t1   = t2*t1
            nexp = 4
            !      Loop over terms in dispersion expansion
            do id = 1, 3
                nexp  = nexp + 2
                fnexp = dble(nexp)
                t2    = t2*t
                t1    = t1*t
                !******SINGLET
                denom        = h3dmcm_h2rm + 2.5d0*h3dmcm_h2r0
                dmpcm_y      = 2.0d0/denom
                x            = rsave*dmpcm_y
                dmpcm_parb   = h3dmcm_dampb(id)
                dmpcm_tt     = h3dmcm_dampb(id)*x
                dmpcm_tt1    = h3dmcm_dampa(id) + dmpcm_tt
                pol          = x*dmpcm_tt1
                dmpcm_dpol   = dmpcm_tt1 + dmpcm_tt
                dmpcm_tt     = exp(-pol)
                dmpcm_tt1    = 1.0d0 - dmpcm_tt
                tt2          = dmpcm_tt1**(nexp-1)
                d            = dmpcm_tt1*tt2
                dd           = 2.0d0*fnexp*dmpcm_dpol*dmpcm_tt*tt2/denom
                cor2cm_corrs = cor2cm_corrs - h3dmcm_chh(id)*d*t2
                dcorrs       = dcorrs - h3dmcm_chh(id)*(dd*t2-fnexp*d*t1)
                call damp(rsave, nexp, d2, d3, d4)
                tem2 = d2
                tem3 = d3
                tem4 = d4
                !*****TRIPLET
                denom        = h3dmcm_h2rmt + 2.5d0*h3dmcm_h2r0
                dmpcm_y      = 2.0d0/denom
                x            = rsave*dmpcm_y
                dmpcm_tt     = h3dmcm_dampb(id)*x
                dmpcm_tt1    = h3dmcm_dampa(id) + dmpcm_tt
                pol          = x*dmpcm_tt1
                dmpcm_dpol   = dmpcm_tt1 + dmpcm_tt
                dmpcm_tt     = exp(-pol)
                dmpcm_tt1    = 1.0d0 - dmpcm_tt
                tt2          = dmpcm_tt1**(nexp-1)
                d            = dmpcm_tt1*tt2
                dd           = 2.0d0*fnexp*dmpcm_dpol*dmpcm_tt*tt2/denom
                cor2cm_corrt = cor2cm_corrt - h3dmcm_chh(id)*d*t2
                dcorrt       = dcorrt - h3dmcm_chh(id)*(dd*t2 - fnexp*d*t1)
                call damp(rsave, nexp, d2, d3, d4)
                cor2cm_dcor2(2) = cor2cm_dcor2(2)+h3dmcm_chh(id)*(d2-tem2)
                cor2cm_dcor2(3) = cor2cm_dcor2(3)+h3dmcm_chh(id)*(d3-tem3)
                cor2cm_dcor2(4) = cor2cm_dcor2(4)+h3dmcm_chh(id)*(d4-tem4)
            enddo
            cor2cm_dcor2(1) = dcorrt-dcorrs
        end subroutine h2cor2

        !> @brief Evaluate DAMP for the H3_GEN_DMBE_1987 potential.
        subroutine damp( rmp, n, d2, d3, d4 )
            !      CALCULATES 2ND,3RD AND 4TH DERIVATIVES OF 2-BODY CORRELATION
            !      ENERGIES FOR SINGLET AND TRIPLET STATES FOR USE IN NON-ADIABATIC
            !      COUPLINGS OF H3 SURFACES 1 AND 2. SENT TO SUBROUTINE H3COR2
            implicit none
            double precision            :: bm2, bm3, bm4, d1b, d2, d2b, d3, d3b, d4, d4b, dpol2, fr, rmp, rr, rrn, t, th, tm, tm2, tm3, tp, tp2, tp3, tw
            integer                     :: n
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,
            !     +                AL0,AL1,AL2,AL3,AZ2,
            !     +                BETA1,BETA2,BETA3,BET0,BET1,
            !     +                CD0,CD1,CHH(3),CK0(3),CK1(3),
            !     +                DAMPA(3),dmpcm_dampb(3),
            !     +                HFD,HFA1,HFA2,HFA3,HFGAM,
            !     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,
            !     +                RHOL,SQRT3,XPAR(15)
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !      COMMON /H3COCM/ PER,PER2,R12,R22,R32,
            !     +                QCOORD,DQ1,DQ2,DQ3,
            !     +                RHO,DRHO1,DRHO2,DRHO3,
            !     +                S,S2,DS1,DS2,DS3,
            !     +                CPHI3,DCPHI1,DCPHI2,DCPHI3
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !      COMMON /H3CRCM/ CORRS(3),DCORRS(3),
            !     +                CORRT(3),DCORRT(3),
            !     +                DAMP(3,3),DDAMP(3,3)
            dpol2 = dmpcm_dpol**2
            d1b   = dmpcm_y*dmpcm_dpol*dmpcm_texp
            d2b   = dmpcm_y**2*(2.0d0*dmpcm_dampb-dpol2)*dmpcm_texp
            d3b   = dmpcm_y**3*dmpcm_dpol*(dpol2-6.0d0*dmpcm_dampb)*dmpcm_texp
            d4b   = dmpcm_y**4*(-12.d0*dmpcm_dampb*dmpcm_dampb+dpol2*(12.d0*dmpcm_dampb-dpol2))*dmpcm_texp
            t     = dble(n)
            tm    = dble(n-1)
            tm2   = dble(n-2)
            tm3   = dble(n-3)
            tp    = dble(n+1)
            tp2   = dble(n+2)
            tp3   = dble(n+3)
            th    = 3.0d0
            tw    = 2.0d0
            fr    = 4.0d0
            rr    = 1.0d0/rmp
            rrn   = rr**n
            bm4   = dmpcm_b**(n-4)
            bm3   = dmpcm_b*bm4
            bm2   = dmpcm_b*bm3
            d2    = t*bm2*(-tp*(dmpcm_b*rr)**2+tw*t*dmpcm_b*rr*d1b-tm*d1b**2-dmpcm_b*d2b)*rrn
            d3=t*bm3*(tp*tp2*(dmpcm_b*rr)**3-th*t*tp*(dmpcm_b*rr)**2*d1b+th*tm*t*dmpcm_b*d1b**2* rr-tm2*tm*d1b**3+th* &
                & t*rr*dmpcm_b**2*d2b-th*tm*dmpcm_b*d1b* d2b-dmpcm_b**2*d3b)*rrn
            d4=t*bm4*(-tp*tp2*tp3*(dmpcm_b*rr)**4+fr*t*tp*tp2*(dmpcm_b*rr)**3*d1b-tw*th*tm *t*tp*(dmpcm_b*rr*d1b)**2+ &
                & fr*tm2*tm*t*dmpcm_b*rr*d1b**3-tw* th*t*tp*dmpcm_b*(dmpcm_b*rr)** 2*d2b+th*fr*tm*t*dmpcm_b*dmpcm_b* &
                & rr*d1b*d2b-tw*th*tm2*tm*dmpcm_b*d1b*d1b*d2b+fr*t* rr*dmpcm_b**3*d3b-tm3*tm2*tm* d1b**4-th*tm* &
                & (dmpcm_b*d2b)**2-fr*tm*dmpcm_b*dmpcm_b*d1b *d3b-dmpcm_b**3*d4b)*rrn
        end subroutine damp

        !> @brief Evaluate KAPPA for the H3_GEN_DMBE_1987 potential.
        subroutine kappa( qlc, capk, dcapk )
            !     CALCULATES FUNCTION CAPK AND ITS FIRST ROUR DERIVATIVES
            implicit none
            double precision            :: at, capk, d1uehf, d1vhf, d1w, d2uehf, d2vhf, d2w, d3uehf, d3vhf, d3w, d4uehf, d4vhf, d4w, &
                & dat, dcapk(4), desdr, dr, dwe, es, fzero, o1, o2, omfz, p1, p2, p3, q1, q2, q3, qlc, rr, rri, t, t1, we
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !  Common block for H3 potential parameters, set in BLOCK DATA PTPACM
            !                                                     *
            !  Common block for coordinates (computed in H3COOR)  *
            !                                                     *
            !                                                 *
            !  Common block for 2-body correlation energies   *
            !  and damping factors (computed in H3COR2)       *
            !                                                 *
            !      COMMON /H3CRCM/ cor2cm_corrs(3),DCORRS(3),
            !     +                cor2cm_corrt(3),DCORRT(3),
            !     +                DAMP(3,3),DDAMP(3,3)
            call h2cor2(qlc)
            t1    = -h3dmcm_az2*h3cocm_qcoord
            fzero = exp(t1*h3cocm_qcoord)
            rr    = qlc
            rri   = 1.0d0/rr
            !*****CALCULATE 1-4 DERIVATIVES OF VEHF
            dr    = rr - h3dmcm_h2rm
            t     = exp(-h3dmcm_hfgam*dr)
            es    = -h3dmcm_hfd*(1.0d0 + dr*(h3dmcm_hfa1 + dr*(h3dmcm_hfa2 + dr*h3dmcm_hfa3)))*t
            desdr = -h3dmcm_hfgam*es - h3dmcm_hfd*(h3dmcm_hfa1 + dr*(2.0d0*h3dmcm_hfa2 + 3.0d0*dr*h3dmcm_hfa3))*t
            p1    = h3dmcm_hfa1+2.0d0*h3dmcm_hfa2*dr+3.0d0*h3dmcm_hfa3*dr*dr
            p3    = 6.0d0*h3dmcm_hfa3
            p2    = 2.0d0*h3dmcm_hfa2 +p3*dr
            d1vhf = desdr
            d2vhf = -h3dmcm_hfgam*d1vhf-h3dmcm_hfd*t*(p2-h3dmcm_hfgam*p1)
            d3vhf = -h3dmcm_hfgam*d2vhf-h3dmcm_hfd*t*(p3-2.d0*h3dmcm_hfgam*p2+h3dmcm_hfgam*h3dmcm_hfgam*p1)
            d4vhf    = -h3dmcm_hfgam*d3vhf+h3dmcm_hfd*t*(3.d0*h3dmcm_hfgam*p3-3.d0*h3dmcm_hfgam*h3dmcm_hfgam*p2+ &
                & h3dmcm_hfgam*h3dmcm_hfgam*h3dmcm_hfgam*p1)
            dcapk(1) = d1vhf
            dcapk(2) = d2vhf
            dcapk(3) = d3vhf
            dcapk(4) = d4vhf
            !*****CALCULATE 1-4 DERIVATIVES OF UEHF
            t      = rr*(h3dmcm_h2tb1 + rr*(h3dmcm_h2tb2 + rr*(h3dmcm_h2tb3 + rr*h3dmcm_h2tb4)))
            at     = h3dmcm_h2ta*exp(-t)*rri
            t      = h3dmcm_h2tb1 + rr*(2.0d0*h3dmcm_h2tb2 + rr*(3.0d0*h3dmcm_h2tb3 + rr*4.0d0*h3dmcm_h2tb4))+rri
            dat    = -at*t
            q1     = -rri**2+2.d0*h3dmcm_h2tb2+6.d0*h3dmcm_h2tb3*rr+12.d0*h3dmcm_h2tb4*rr**2
            q2     = 2.d0*rri**3+6.d0*h3dmcm_h2tb3+24.d0*h3dmcm_h2tb4*rr
            q3     = -6.d0*rri**4+24.d0*h3dmcm_h2tb4
            d1uehf = dat
            d2uehf = -t*d1uehf-at*q1
            d3uehf = -t*d2uehf-2.d0*q1*d1uehf-at*q2
            d4uehf = -t*d3uehf-3.d0*q1*d2uehf-3.d0*q2*d1uehf-at*q3
            !     Add in triplet and subtract singlet 2-body correlation terms
            at = at + cor2cm_corrt - cor2cm_corrs
            !*****INCLUDE CORRELATION PART OF CAPK,(USUBC-VSUBC)*(1-FZERO)
            omfz     = 1.0d0-fzero
            dcapk(1) = dcapk(1)-(d1uehf+cor2cm_dcor2(1))*omfz
            dcapk(2) = dcapk(2)-(d2uehf+cor2cm_dcor2(2))*omfz
            dcapk(3) = dcapk(3)-(d3uehf+cor2cm_dcor2(3))*omfz
            dcapk(4) = dcapk(4)-(d4uehf+cor2cm_dcor2(4))*omfz
            !*****CALCULATE 1-4 DERIVATIVES OF W-EFFECTIVE
            t        = exp(-h3dmcm_al3*rr)
            we       = (h3dmcm_al0 + rr*(h3dmcm_al1 + rr*h3dmcm_al2))*t
            dwe      = (h3dmcm_al1 + 2.0d0*h3dmcm_al2*rr)*t - h3dmcm_al3*we
            o2       = 2.0d0*h3dmcm_al2
            o1       = h3dmcm_al1+o2*rr
            d1w      = dwe
            d2w      = -h3dmcm_al3*d1w+t*(o2-h3dmcm_al3*o1)
            d3w      = -h3dmcm_al3*d2w+t*(h3dmcm_al3*o1-2.0d0*o2)*h3dmcm_al3
            d4w      = -h3dmcm_al3*d3w+t*(3.0d0*o2-h3dmcm_al3*o1)*h3dmcm_al3**2
            dcapk(1) = dcapk(1)-fzero*d1w
            dcapk(2) = dcapk(2)-fzero*d2w
            dcapk(3) = dcapk(3)-fzero*d3w
            dcapk(4) = dcapk(4)-fzero*d4w
            capk     = es-fzero*we-omfz*at
        end subroutine kappa

        !> @brief Evaluate POTINFO for the H3_GEN_DMBE_1987 potential.
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

        !> @brief Evaluate ANCVRT for the H3_GEN_DMBE_1987 potential.
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

        !> @brief Evaluate CARTOU for the H3_GEN_DMBE_1987 potential.
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

        !> @brief Evaluate CARTTOR for the H3_GEN_DMBE_1987 potential.
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

        !> @brief Evaluate EUNITZERO for the H3_GEN_DMBE_1987 potential.
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

        !> @brief Evaluate RTOCART for the H3_GEN_DMBE_1987 potential.
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

        !> @brief Evaluate DEDCOU for the H3_GEN_DMBE_1987 potential.
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
end module chempotpy_h3_gen_dmbe_1987_pes

!> @brief Expose H3_GEN_DMBE_1987 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_dmbe_1987_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 2
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
