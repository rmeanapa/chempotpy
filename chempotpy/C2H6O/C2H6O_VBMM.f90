!   System:                        C2H6O
!   Common name:                   C2H6+O
!   Number of derivatives:         1
!   Number of bodies:              9
!   Number of electronic surfaces: 1
!   Interface:                     potlib2001
!   References:: J. Espinosa-Garcia, C. Rangel, J.C. Corchado and
!   M. Garcia-Chamorro, PCCP, 22, 22591 (2020)

!> @brief Implement C2H6O_VBMM, a potential energy surface for C2H6O.
module chempotpy_c2h6o_vbmm_pes
    use chempotpy_geometry,  only: cross_product_3d
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the C2H6O_VBMM potential energy surface for C2H6O.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: hartree_to_electronvolt
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 9
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter    :: n3atom  = 75
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: angles_dtheta0(8,8,8)
        real(wp)              :: angles_theta0(8,8)
        real(wp)              :: bonds2_rbh(6)
        real(wp)              :: bonds2_rc1h(3)
        real(wp)              :: bonds2_rc2h(3)
        real(wp)              :: bonds2_rhh(3,6)
        real(wp)              :: bonds3_rc1b
        real(wp)              :: bonds3_rc2b
        real(wp)              :: bonds_r0c1h
        real(wp)              :: bonds_r0c2h
        real(wp)              :: coords1_tbh(6,3)
        real(wp)              :: coords1_tc1h(3,3)
        real(wp)              :: coords1_tc2h(3,3)
        real(wp)              :: coords1_tcb(3)
        real(wp)              :: coords2_tcc(3)
        real(wp)              :: coords2_thh(3,6,3)
        real(wp)              :: coords3_tc1b(3)
        real(wp)              :: coords3_tc2b(3)
        real(wp)              :: delta1_fdelta(8)
        real(wp)              :: delta1_hdelta(8)
        real(wp)              :: delta2_dfdelta(8,8)
        real(wp)              :: delta2_dhdelta(8,8)
        real(wp)              :: force1_dfdc(8,8,8)
        real(wp)              :: force1_dfdh(8,8,8)
        real(wp)              :: force1_f1(8)
        real(wp)              :: force1_fk0(8,8)
        real(wp)              :: infocm_cartnu(natom,3)
        integer               :: infocm_icartr
        integer               :: infocm_indexes(natom)
        integer               :: infocm_irctnt
        integer               :: infocm_mder
        integer               :: infocm_msurf
        integer               :: infocm_natoms
        character(len=75)     :: infocm_ref(5)
        real(wp)              :: ip1_ds1(8)
        real(wp)              :: ip1_ds2(8)
        real(wp)              :: ip1_s1(8)
        real(wp)              :: ip1_s2(8)
        integer               :: ndx2_nc1(3)
        integer               :: ndx2_nc2(3)
        integer               :: ndx2_nnh(6)
        integer               :: ndx_nc(3)
        integer               :: ndx_nh(8,3)
        integer               :: ndx_nhb(3)
        real(wp)              :: op1_ds3(8)
        real(wp)              :: op1_s3(8)
        real(wp)              :: potcm_a1c
        real(wp)              :: potcm_a1c1h
        real(wp)              :: potcm_a1c2h
        real(wp)              :: potcm_a1s
        real(wp)              :: potcm_a1scc
        real(wp)              :: potcm_a2c
        real(wp)              :: potcm_a2s
        real(wp)              :: potcm_a2scc
        real(wp)              :: potcm_a3s
        real(wp)              :: potcm_a3scc
        real(wp)              :: potcm_aa1
        real(wp)              :: potcm_aa2
        real(wp)              :: potcm_aa3
        real(wp)              :: potcm_aa4
        real(wp)              :: potcm_acb
        real(wp)              :: potcm_acc
        real(wp)              :: potcm_ahh
        real(wp)              :: potcm_ak
        real(wp)              :: potcm_aphi
        real(wp)              :: potcm_atheta
        real(wp)              :: potcm_b1c1h
        real(wp)              :: potcm_b1c2h
        real(wp)              :: potcm_b1s
        real(wp)              :: potcm_b1scc
        real(wp)              :: potcm_b2s
        real(wp)              :: potcm_b2scc
        real(wp)              :: potcm_b3s
        real(wp)              :: potcm_b3scc
        real(wp)              :: potcm_bk
        real(wp)              :: potcm_bphi
        real(wp)              :: potcm_btheta
        real(wp)              :: potcm_c1c1h
        real(wp)              :: potcm_c1c2h
        real(wp)              :: potcm_cphi
        real(wp)              :: potcm_ctheta
        real(wp)              :: potcm_d1c1h
        real(wp)              :: potcm_d1c2h
        real(wp)              :: potcm_d1cb
        real(wp)              :: potcm_d1cc
        real(wp)              :: potcm_d1hh
        real(wp)              :: potcm_d3c1h
        real(wp)              :: potcm_d3c2h
        real(wp)              :: potcm_d3cb
        real(wp)              :: potcm_d3hh
        real(wp)              :: potcm_fch3
        real(wp)              :: potcm_fkinf
        real(wp)              :: potcm_hch3
        integer               :: potcm_nnb
        integer               :: potcm_nnc1
        integer               :: potcm_nnc2
        real(wp)              :: potcm_r0c1hp
        real(wp)              :: potcm_r0c1hr
        real(wp)              :: potcm_r0c2hp
        real(wp)              :: potcm_r0c2hr
        real(wp)              :: potcm_r0cb
        real(wp)              :: potcm_r0cc
        real(wp)              :: potcm_r0hh
        real(wp)              :: potcm_rcc
        real(wp)              :: potcm_v3
        real(wp)              :: potcm_w1
        real(wp)              :: potcm_w2
        real(wp)              :: potcm_w3
        real(wp)              :: potcm_w4
        real(wp)              :: potcm_w5
        real(wp)              :: potcm_w6
        real(wp)              :: pt1cm_degsdr(n3atom)
        real(wp)              :: pt1cm_engygs
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: pt4cm_deesdr(n3atom,isurf)
        real(wp)              :: pt4cm_engyes(isurf)
        real(wp)              :: pt5cm_deijdr(n3atom,jsurf)
        real(wp)              :: pt5cm_engyij(jsurf)
        real(wp)              :: qpdot_pdot(27)
        real(wp)              :: qpdot_q(27)
        real(wp)              :: stret1_a1cc
        real(wp)              :: stret1_b1cc
        real(wp)              :: stret1_c1cc
        real(wp)              :: stret1_d1ch
        real(wp)              :: stret1_d3ch
        real(wp)              :: switch1_dsphi(8)
        real(wp)              :: switch1_dstheta(8)
        real(wp)              :: switch1_sphi(8)
        real(wp)              :: switch1_stheta(8)
        real(wp)              :: thetahcch_thetahcch(9)
        real(wp)              :: usricm_anuzero
        real(wp)              :: usricm_cart(natom,3)
        integer               :: usricm_nasurf(isurf+1,isurf+1)
        integer               :: usricm_nder
        integer               :: usricm_nflag(20)
        integer               :: usricm_nulbl(natom)
        real(wp)              :: usrocm_descart(natom,3,isurf)
        real(wp)              :: usrocm_dgscart(natom,3)
        real(wp)              :: usrocm_dijcart(natom,3,jsurf)
        real(wp)              :: usrocm_pengyes(isurf)
        real(wp)              :: usrocm_pengygs
        real(wp)              :: usrocm_pengyij(jsurf)
        real(wp)              :: utilcm_cnvrtd
        real(wp)              :: utilcm_cnvrtde
        real(wp)              :: utilcm_cnvrte
        real(wp)              :: utilcm_descartnu(natom,3,isurf)
        real(wp)              :: utilcm_dgscartnu(natom,3)
        real(wp)              :: utilcm_dijcartnu(natom,3,jsurf)
        integer               :: utilcm_ireorder
        integer               :: utilcm_kediag
        integer               :: utilcm_keoffd
        integer               :: utilcm_ksdiag
        integer               :: utilcm_ksoffd
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        angles_dtheta0      = 0.e0_wp; angles_theta0 = 0.e0_wp; bonds2_rbh = 0.e0_wp; bonds2_rc1h = 0.e0_wp; bonds2_rc2h = 0.e0_wp
        bonds2_rhh          = 0.e0_wp; bonds3_rc1b = 0.e0_wp; bonds3_rc2b = 0.e0_wp; bonds_r0c1h = 0.e0_wp; bonds_r0c2h = 0.e0_wp
        coords1_tbh         = 0.e0_wp; coords1_tc1h = 0.e0_wp; coords1_tc2h = 0.e0_wp; coords1_tcb = 0.e0_wp; coords2_tcc = 0.e0_wp
        coords2_thh         = 0.e0_wp; coords3_tc1b = 0.e0_wp; coords3_tc2b = 0.e0_wp; delta1_fdelta = 0.e0_wp; delta1_hdelta = 0.e0_wp
        delta2_dfdelta      = 0.e0_wp; delta2_dhdelta = 0.e0_wp; force1_dfdc = 0.e0_wp; force1_dfdh = 0.e0_wp; force1_f1 = 0.e0_wp
        force1_fk0          = 0.e0_wp; infocm_cartnu = 0.e0_wp; infocm_icartr = 0; infocm_indexes = 0; infocm_irctnt = 0
        infocm_mder         = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref          = ' '
        ip1_ds1             = 0.e0_wp
        ip1_ds2             = 0.e0_wp
        ip1_s1              = 0.e0_wp
        ip1_s2              = 0.e0_wp
        ndx2_nc1            = 0
        ndx2_nc2            = 0
        ndx2_nnh            = 0
        ndx_nc              = 0
        ndx_nh              = 0
        ndx_nhb             = 0
        op1_ds3             = 0.e0_wp
        op1_s3              = 0.e0_wp
        potcm_a1c           = 0.e0_wp
        potcm_a1c1h         = 0.e0_wp
        potcm_a1c2h         = 0.e0_wp
        potcm_a1s           = 0.e0_wp
        potcm_a1scc         = 0.e0_wp
        potcm_a2c           = 0.e0_wp
        potcm_a2s           = 0.e0_wp
        potcm_a2scc         = 0.e0_wp
        potcm_a3s           = 0.e0_wp
        potcm_a3scc         = 0.e0_wp
        potcm_aa1           = 0.e0_wp
        potcm_aa2           = 0.e0_wp
        potcm_aa3           = 0.e0_wp
        potcm_aa4           = 0.e0_wp
        potcm_acb           = 0.e0_wp
        potcm_acc           = 0.e0_wp
        potcm_ahh           = 0.e0_wp
        potcm_ak            = 0.e0_wp
        potcm_aphi          = 0.e0_wp
        potcm_atheta        = 0.e0_wp
        potcm_b1c1h         = 0.e0_wp
        potcm_b1c2h         = 0.e0_wp
        potcm_b1s           = 0.e0_wp
        potcm_b1scc         = 0.e0_wp
        potcm_b2s           = 0.e0_wp
        potcm_b2scc         = 0.e0_wp
        potcm_b3s           = 0.e0_wp
        potcm_b3scc         = 0.e0_wp
        potcm_bk            = 0.e0_wp
        potcm_bphi          = 0.e0_wp
        potcm_btheta        = 0.e0_wp
        potcm_c1c1h         = 0.e0_wp
        potcm_c1c2h         = 0.e0_wp
        potcm_cphi          = 0.e0_wp
        potcm_ctheta        = 0.e0_wp
        potcm_d1c1h         = 0.e0_wp
        potcm_d1c2h         = 0.e0_wp
        potcm_d1cb          = 0.e0_wp
        potcm_d1cc          = 0.e0_wp
        potcm_d1hh          = 0.e0_wp
        potcm_d3c1h         = 0.e0_wp
        potcm_d3c2h         = 0.e0_wp
        potcm_d3cb          = 0.e0_wp
        potcm_d3hh          = 0.e0_wp
        potcm_fch3          = 0.e0_wp
        potcm_fkinf         = 0.e0_wp
        potcm_hch3          = 0.e0_wp
        potcm_nnb           = 0
        potcm_nnc1          = 0
        potcm_nnc2          = 0
        potcm_r0c1hp        = 0.e0_wp
        potcm_r0c1hr        = 0.e0_wp
        potcm_r0c2hp        = 0.e0_wp
        potcm_r0c2hr        = 0.e0_wp
        potcm_r0cb          = 0.e0_wp
        potcm_r0cc          = 0.e0_wp
        potcm_r0hh          = 0.e0_wp
        potcm_rcc           = 0.e0_wp
        potcm_v3            = 0.e0_wp
        potcm_w1            = 0.e0_wp
        potcm_w2            = 0.e0_wp
        potcm_w3            = 0.e0_wp
        potcm_w4            = 0.e0_wp
        potcm_w5            = 0.e0_wp
        potcm_w6            = 0.e0_wp
        pt1cm_degsdr        = 0.e0_wp
        pt1cm_engygs        = 0.e0_wp
        pt1cm_r             = 0.e0_wp
        pt3cm_ezero         = 0.e0_wp
        pt4cm_deesdr        = 0.e0_wp
        pt4cm_engyes        = 0.e0_wp
        pt5cm_deijdr        = 0.e0_wp
        pt5cm_engyij        = 0.e0_wp
        qpdot_pdot          = 0.e0_wp
        qpdot_q             = 0.e0_wp
        stret1_a1cc         = 0.e0_wp
        stret1_b1cc         = 0.e0_wp
        stret1_c1cc         = 0.e0_wp
        stret1_d1ch         = 0.e0_wp
        stret1_d3ch         = 0.e0_wp
        switch1_dsphi       = 0.e0_wp
        switch1_dstheta     = 0.e0_wp
        switch1_sphi        = 0.e0_wp
        switch1_stheta      = 0.e0_wp
        thetahcch_thetahcch = 0.e0_wp
        usricm_anuzero      = 0.e0_wp
        usricm_cart         = 0.e0_wp
        usricm_nasurf       = 0
        usricm_nder         = 0
        usricm_nflag        = 0
        usricm_nulbl        = 0
        usrocm_descart      = 0.e0_wp
        usrocm_dgscart      = 0.e0_wp
        usrocm_dijcart      = 0.e0_wp
        usrocm_pengyes      = 0.e0_wp
        usrocm_pengygs      = 0.e0_wp
        usrocm_pengyij      = 0.e0_wp
        utilcm_cnvrtd       = 0.e0_wp
        utilcm_cnvrtde      = 0.e0_wp
        utilcm_cnvrte       = 0.e0_wp
        utilcm_descartnu    = 0.e0_wp
        utilcm_dgscartnu    = 0.e0_wp
        utilcm_dijcartnu    = 0.e0_wp
        utilcm_ireorder     = 0
        utilcm_kediag       = 0
        utilcm_keoffd       = 0
        utilcm_ksdiag       = 0
        utilcm_ksoffd       = 0
        usricm_nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0], shape(usricm_nasurf))
        usricm_nder    = 1
        usricm_nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero = 0.0e0_wp
        infocm_icartr  = 1
        infocm_msurf   = 0
        infocm_mder    = 1
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 9
        potcm_nnc1     = 1
        potcm_nnc2     = 2
        potcm_nnb      = 9
        ndx2_nnh       = [3, 4, 5, 6, 7, 8]
        potcm_r0c1hr   = 1.08897e0_wp
        potcm_r0c1hp   = 1.07577e0_wp
        potcm_w3       = 1.00000e0_wp
        potcm_w4       = 1.08897e0_wp
        potcm_d1c1h    = 129.23000e0_wp
        potcm_d3c1h    = 36.00400e0_wp
        potcm_a1c1h    = 1.60300e0_wp
        potcm_b1c1h    = 0.12500e0_wp
        potcm_c1c1h    = 3.51404e0_wp
        potcm_r0c2hr   = 1.08897e0_wp
        potcm_r0c2hp   = 1.07577e0_wp
        potcm_w5       = 1.00000e0_wp
        potcm_w6       = 1.08897e0_wp
        potcm_d1c2h    = 129.23000e0_wp
        potcm_d3c2h    = 36.00400e0_wp
        potcm_a1c2h    = 1.60300e0_wp
        potcm_b1c2h    = 0.12500e0_wp
        potcm_c1c2h    = 3.51404e0_wp
        potcm_r0cc     = 1.50007e0_wp
        potcm_d1cc     = 117.94900e0_wp
        potcm_acc      = 1.93400e0_wp
        potcm_a1c      = 0.90007e0_wp
        potcm_a2c      = 0.94900e0_wp
        potcm_r0hh     = 0.97101e0_wp
        potcm_d1hh     = 127.45800e0_wp
        potcm_d3hh     = 40.36400e0_wp
        potcm_ahh      = 2.10400e0_wp
        potcm_r0cb     = 1.77097e0_wp
        potcm_d1cb     = 48.20900e0_wp
        potcm_d3cb     = 15.66000e0_wp
        potcm_acb      = 1.3290285e0_wp
        potcm_a3s      = 1.8419147e0_wp
        potcm_b3s      = 1.0889300e0_wp
        potcm_aphi     = 1.0087903e0_wp
        potcm_bphi     = 1.2006638e0_wp
        potcm_cphi     = 1.0889937e0_wp
        potcm_atheta   = 1.2078714e0_wp
        potcm_btheta   = 0.3548859e0_wp
        potcm_ctheta   = 1.8915497e0_wp
        potcm_fch3     = 0.1357500e0_wp
        potcm_hch3     = 0.2615000e0_wp
        potcm_fkinf    = 0.2977000e0_wp
        potcm_ak       = 0.3060000e0_wp
        potcm_bk       = 0.2032000e0_wp
        potcm_aa1      = 0.500952e0_wp
        potcm_aa2      = 0.499963e0_wp
        potcm_aa3      = 6.465953e0_wp
        potcm_aa4      = 1.569977e0_wp
        potcm_v3       = 13.730000e0_wp
        potcm_w1       = 9.3548859e0_wp
        potcm_w2       = 0.913970e0_wp
        potcm_a1s      = 2.5313681e0_wp
        potcm_b1s      = 1.088625e0_wp
        potcm_a2s      = 2.0147402e0_wp
        potcm_b2s      = 1.088980e0_wp
        potcm_a1scc    = 2.5313681e0_wp
        potcm_b1scc    = 1.550070e0_wp
        potcm_a2scc    = 2.0147402e0_wp
        potcm_b2scc    = 1.550070e0_wp
        potcm_a3scc    = 0.531368e0_wp
        potcm_b3scc    = 1.550070e0_wp
        v              = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        usricm_cart    = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                usricm_cart(iatom,idir) = x(iatom,idir)
            enddo
        enddo
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives are unavailable')
        elseif( igrad == 2 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives are unavailable')
        endif

        contains

        !> @brief Evaluate POT for the C2H6O_VBMM potential.
        subroutine pot
            !   System: C2H6 + O --> C2H5 + HO.
            !              Joaquin Espinosa-Garcia Junio 2020
            !              Basada en la PES H + C2H6 REv-2019 para RPMD.
            !   All the information passed to and from the potential energy surface
            !   routine is in hartree atomic units.
            !        This potential is written such that:
            !                       X(1)  - X(3)  : X, Y, Z for C1
            !                       X(4)  - X(6)  : X, Y, Z for C2
            !                       X(7)  - X(9)  : X, Y, Z for H1
            !                       X(10) - X(12) : X, Y, Z for H2
            !                       X(13) - X(15) : X, Y, Z for H3
            !                       X(16) - X(18) : X, Y, Z for H4
            !                       X(19) - X(21) : X, Y, Z for H5
            !                       X(22) - X(24) : X, Y, Z for H6
            !                       X(25) - X(27) : X, Y, Z for Hb
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: en, vip, vop, vstr, vtor
            integer  :: i, n
            ! Start of changes by EMN.
            ! The lines below do not seem to be needed because the Cartesian
            ! coordinates are in common qpdot.
            ! Besides, BOHR_A in utility.f does not match the value below.
            ! As a consequence, qpdot_q(i) changes and energy is not conserved
            call cartou
            call carttor
            !     PUT COORDINATES IN PROPER ARRAYS
            ! Inicialization de variables
            do i=1, 27
                qpdot_q(i)    = pt1cm_r(i)
                qpdot_pdot(i) = 0.e0_wp
            enddo
            ! Changing units from angstrom to bohr
            ! End of corrections by EMN
            en = 0.0e0_wp
            call coorden
            !  calculate switching functions
            call switchf
            !  calculate reference angles and their derivatives
            call refangles
            ! CRD 2011 incluimos los terminos de potencial de torsion
            call torsion(vtor)
            !  calculate stretching potential
            call stretch(vstr)
            !  calculate out of plane bending potential
            call opbend(vop)
            !  calculate in plane bending potential
            call ipbend(vip)
            en = vstr+vop+vip+vtor
            !  convert from 10(5) j/mol to au
            en           = en*0.03812e0_wp
            pt1cm_engygs = en
            call eunitzero
            if( usricm_nder .ne. 0 )then
                !  Initialize derivatives
                do n=1, 27
                    pt1cm_degsdr(n) = qpdot_pdot(n)*0.0201723e0_wp
                enddo
                ! Start of changes by EMN.
                ! in venuspotlib_3.f. Anyway usricm_nder=1
                ! Convert from 10^5 J/mol to hartree.
                ! 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
                ! qpdot_pdot sale de la derivada en hartrees/bohr por eso no lo multiplico
                ! End of corrections by EMN
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate coorden for the C2H6O_VBMM potential.
        subroutine coorden
            !  calculates relative coordinates and bond lengths
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: argmax, argp3, argp4, p3, p4, t3tmp, t4tmp
            integer  :: i, ind
            !  calculate relative coordinates
            do ind=1, 3
                coords3_tc1b(ind) = qpdot_q(ndx2_nc1(ind))-qpdot_q(ndx_nhb(ind))
                coords3_tc2b(ind) = qpdot_q(ndx2_nc2(ind))-qpdot_q(ndx_nhb(ind))
                coords2_tcc(ind)  = qpdot_q(ndx2_nc1(ind))-qpdot_q(ndx2_nc2(ind))
            enddo
            do ind=1, 3
                do i=1, 3
                    coords1_tc1h(i,ind)  = qpdot_q(ndx2_nc1(ind))-qpdot_q(ndx_nh(i,ind))
                    coords1_tc2h(i,ind)  = qpdot_q(ndx2_nc2(ind))-qpdot_q(ndx_nh(i+3,ind))
                    coords1_tbh(i,ind)   = qpdot_q(ndx_nhb(ind))-qpdot_q(ndx_nh(i,ind))
                    coords1_tbh(i+3,ind) = qpdot_q(ndx_nhb(ind))-qpdot_q(ndx_nh(i+3,ind))
                enddo
            enddo
            ! coordenadas entre los hidrogenos de los dos carbonos
            !  calculate bond lengths
            bonds3_rc1b = sqrt(coords3_tc1b(1)*coords3_tc1b(1)+coords3_tc1b(2)*coords3_tc1b(2)+coords3_tc1b(3)* &
                & coords3_tc1b(3))
            bonds3_rc2b = sqrt(coords3_tc2b(1)*coords3_tc2b(1)+coords3_tc2b(2)*coords3_tc2b(2)+coords3_tc2b(3)* &
                & coords3_tc2b(3))
            potcm_rcc = sqrt(coords2_tcc(1)*coords2_tcc(1)+coords2_tcc(2)*coords2_tcc(2)+coords2_tcc(3)*coords2_tcc(3))
            do i=1, 3
                bonds2_rc1h(i) = sqrt(coords1_tc1h(i,1)*coords1_tc1h(i,1)+coords1_tc1h(i,2)*coords1_tc1h(i,2)+ &
                    & coords1_tc1h(i,3)*coords1_tc1h(i,3))
                ! redondeo para corregir las frecuencias
                bonds2_rc2h(i) = sqrt(coords1_tc2h(i,1)*coords1_tc2h(i,1)+coords1_tc2h(i,2)*coords1_tc2h(i,2)+ &
                    & coords1_tc2h(i,3)*coords1_tc2h(i,3))
                ! redondeo para corregir las frecuencias
            enddo
            do i=1, 6
                bonds2_rbh(i) = sqrt(coords1_tbh(i,1)*coords1_tbh(i,1)+coords1_tbh(i,2)*coords1_tbh(i,2)+ &
                    & coords1_tbh(i,3)*coords1_tbh(i,3))
            enddo
            ! calculamos las distancias H-H
            !    *                coords2_thh(i,j,3)*coords2_thh(i,j,3))
            ! crd 2013 modificaciones de las distancias
            argmax = 19.e0_wp
            p3     = 1.e0_wp
            do i=1, 3
                argp3 = (potcm_w3*(bonds2_rc1h(i)-potcm_w4))
                if( argp3 .lt. argmax )then
                    t3tmp = 1-tanh(argp3)
                else
                    t3tmp = 0.e0_wp
                endif
                p3 = p3*t3tmp
            enddo
            bonds_r0c1h = p3*potcm_r0c1hr+(1-p3)*potcm_r0c1hp
            p4          = 1.e0_wp
            do i=1, 3
                argp4 = (potcm_w5*(bonds2_rc2h(i)-potcm_w6))
                if( argp4 .lt. argmax )then
                    t4tmp = 1-tanh(argp4)
                else
                    t4tmp = 0.e0_wp
                endif
                p4 = p4*t4tmp
            enddo
            bonds_r0c2h = p4*potcm_r0c2hr+(1-p4)*potcm_r0c2hp
        end subroutine coorden

        !> @brief Evaluate torsion for the C2H6O_VBMM potential.
        subroutine torsion( vtor )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dihed, t1(3), t2(3), vtor
            integer  :: i
            !      pi=4.0d0*atan(1.0d0)
            ! CRD 2011 incluimos un apartado para calcular los angulos dihedro
            ! para poder calcular despues vtorsion
            call dihedro(3, 1, 2, 6, dihed)
            thetahcch_thetahcch(1) = dihed*PI/180.0e0_wp
            call dihedro(3, 1, 2, 7, dihed)
            thetahcch_thetahcch(2) = dihed*PI/180.0e0_wp
            call dihedro(3, 1, 2, 8, dihed)
            thetahcch_thetahcch(3) = dihed*PI/180.0e0_wp
            call dihedro(4, 1, 2, 6, dihed)
            thetahcch_thetahcch(4) = dihed*PI/180.0e0_wp
            call dihedro(4, 1, 2, 7, dihed)
            thetahcch_thetahcch(5) = dihed*PI/180.0e0_wp
            call dihedro(4, 1, 2, 8, dihed)
            thetahcch_thetahcch(6) = dihed*PI/180.0e0_wp
            call dihedro(5, 1, 2, 6, dihed)
            thetahcch_thetahcch(7) = dihed*PI/180.0e0_wp
            call dihedro(5, 1, 2, 7, dihed)
            thetahcch_thetahcch(8) = dihed*PI/180.0e0_wp
            call dihedro(5, 1, 2, 8, dihed)
            thetahcch_thetahcch(9) = dihed*PI/180.0e0_wp
            ! CRD 2011 calculamos los terminos switching t1 y t2
            do i=1, 3
                t1(i) = 0.5e0_wp*(1-tanh(potcm_w1*(bonds2_rc1h(i)-potcm_w2)))
                t2(i) = 0.5e0_wp*(1-tanh(potcm_w1*(bonds2_rc2h(i)-potcm_w2)))
            enddo
            ! CRD 2011 inicializamos vtor
            vtor = 0.0e0_wp
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(3*thetahcch_thetahcch(1)))*t1(1)*t2(1)
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(3*thetahcch_thetahcch(2)))*t1(1)*t2(2)
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(3*thetahcch_thetahcch(3)))*t1(1)*t2(3)
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(3*thetahcch_thetahcch(4)))*t1(2)*t2(1)
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(3*thetahcch_thetahcch(5)))*t1(2)*t2(2)
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(3*thetahcch_thetahcch(6)))*t1(2)*t2(3)
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(3*thetahcch_thetahcch(7)))*t1(3)*t2(1)
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(3*thetahcch_thetahcch(8)))*t1(3)*t2(2)
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(3*thetahcch_thetahcch(9)))*t1(3)*t2(3)
        end subroutine torsion

        !> @brief Evaluate dihedro for the C2H6O_VBMM potential.
        subroutine dihedro( i, j, k, l, dihed )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dihed
            integer  :: i, j, k, l
            integer  :: m
            real(wp) :: rij(3), rjk(3), rkl(3)
            real(wp) :: zijjk(3), zjkkl(3)
            real(wp) :: arg1, arg2
            !        pi=4.0d0*atan(1.0d0)
            do m = 1, 3
                rij(m) = qpdot_q((i-1)*3+m)-qpdot_q((j-1)*3+m)
                rjk(m) = qpdot_q((j-1)*3+m)-qpdot_q((k-1)*3+m)
                rkl(m) = qpdot_q((k-1)*3+m)-qpdot_q((l-1)*3+m)
            enddo
            zjkkl = cross_product_3d(rjk, rkl)
            zijjk = cross_product_3d(rij, rjk)
            arg2  = norm2(rjk)*dot_product(rij, zjkkl)
            arg1  = dot_product(zijjk, zjkkl)
            if( abs(arg1) .lt. 1.e-6_wp )then
                if( arg2 .gt. 0.e0_wp ) dihed = 90.e0_wp
                if( arg2 .lt. 0.e0_wp ) dihed = -90.e0_wp
                return
            endif
            if( arg1 .gt. 0.e0_wp )then
                dihed = atan(arg2/arg1) *180.e0_wp/PI
                return
            endif
            if( arg2 .ge. 0.e0_wp )then
                dihed = 180.e0_wp + atan(arg2/arg1) *180.e0_wp/PI
                return
            endif
            if( arg2 .lt. 0.e0_wp )then
                dihed = -180.e0_wp + atan(arg2/arg1) *180.e0_wp/PI
                return
            endif
        end subroutine dihedro

        !> @brief Evaluate refangles for the C2H6O_VBMM potential.
        subroutine refangles
            !  subroutine calculates reference angles for the "in-plane" potential
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: halfpi, tau, taucc, twopi
            integer  :: i, j
            tau = acos(-1.0e0_wp/3.0e0_wp)
            ! CRD 2012 taucc para el angulo H-C-C
            taucc  = acos(-0.3616e0_wp)
            halfpi = 0.5e0_wp*PI
            twopi  = 2.0e0_wp*PI
            ! jcc-2010
            ! valor para ajustar el angulo
            ! jcc-2010
            !  set diagonal elements to zero
            do i=1, 8
                angles_theta0(i,i) = 0.0e0_wp
            enddo
            !  calculate reference angles
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(3)*switch1_stheta(4)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(2)*switch1_stheta(4)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(2)*switch1_stheta(3)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(1)*switch1_stheta(4)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(1)*switch1_stheta(3)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(1)*switch1_stheta(2)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(7)*switch1_stheta(8)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(6)*switch1_stheta(8)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(6)*switch1_stheta(7)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(5)*switch1_stheta(8)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(5)*switch1_stheta(7)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(5)*switch1_stheta(6)-1.0d0)
            ! el atomo 4 es el c2
            angles_theta0(1,2) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(2)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(3)*switch1_stheta(4)-1.0e0_wp)
            angles_theta0(1,3) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(3)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(2)*switch1_stheta(4)-1.0e0_wp)
            angles_theta0(1,4) = taucc+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(4)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(2)*switch1_stheta(3)-1.0e0_wp)
            angles_theta0(2,3) = tau+(tau-halfpi)*(switch1_sphi(2)*switch1_sphi(3)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(1)*switch1_stheta(4)-1.0e0_wp)
            angles_theta0(2,4) = taucc+(tau-halfpi)*(switch1_sphi(2)*switch1_sphi(4)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(1)*switch1_stheta(3)-1.0e0_wp)
            angles_theta0(3,4) = taucc+(tau-halfpi)*(switch1_sphi(3)*switch1_sphi(4)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(1)*switch1_stheta(2)-1.0e0_wp)
            !  el atomo 8 es el c1
            angles_theta0(5,6) = tau+(tau-halfpi)*(switch1_sphi(5)*switch1_sphi(6)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(7)*switch1_stheta(8)-1.0e0_wp)
            angles_theta0(5,7) = tau+(tau-halfpi)*(switch1_sphi(5)*switch1_sphi(7)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(6)*switch1_stheta(8)-1.0e0_wp)
            angles_theta0(5,8) = taucc+(tau-halfpi)*(switch1_sphi(5)*switch1_sphi(8)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(6)*switch1_stheta(7)-1.0e0_wp)
            angles_theta0(6,7) = tau+(tau-halfpi)*(switch1_sphi(6)*switch1_sphi(7)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(5)*switch1_stheta(8)-1.0e0_wp)
            angles_theta0(6,8) = taucc+(tau-halfpi)*(switch1_sphi(6)*switch1_sphi(8)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(5)*switch1_stheta(7)-1.0e0_wp)
            angles_theta0(7,8) = taucc+(tau-halfpi)*(switch1_sphi(7)*switch1_sphi(8)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(5)*switch1_stheta(6)-1.0e0_wp)
            !  fill in the other half of the matrix
            do i=1, 3
                do j=i+1, 4
                    angles_theta0(j,i) = angles_theta0(i,j)
                enddo
            enddo
            do i=5, 7
                do j=i+1, 8
                    angles_theta0(j,i) = angles_theta0(i,j)
                enddo
            enddo
        end subroutine refangles

        !> @brief Evaluate stretch for the C2H6O_VBMM potential.
        subroutine stretch( vstr )
            !  subroutine to calculate leps-type stretching potential and its
            !  derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: ach1, ach2, arga1, arga2, dt, e1, e3, expterm, rav1, rav2, vcc, vj(3), vjbh(6), vjc1b, vjc1h(3), vjc2b, vjc2h(3), vq(3), vqbh(6), vqc1b, vqc1h(3), vqc2b, vqc2h(3), vstr
            integer  :: i
            !  calculate avergage bond length for the methane moiety
            rav1 = (bonds2_rc1h(1)+bonds2_rc1h(2)+bonds2_rc1h(3))/3.0e0_wp
            rav2 = (bonds2_rc2h(1)+bonds2_rc2h(2)+bonds2_rc2h(3))/3.0e0_wp
            !  initialise:
            vstr = 0.0e0_wp
            !  ach:
            !  nb: in double precision tanh(19.0d0)=1.0d0 and we put the if statemen
            !  in to avoid overflow/underflow errors
            arga1 = potcm_c1c1h*(rav1-bonds_r0c1h)
            arga2 = potcm_c1c2h*(rav2-bonds_r0c2h)
            if( arga1 .lt. 19.0e0_wp )then
                ach1 = potcm_a1c1h+potcm_b1c1h*(tanh(arga1)+1.0e0_wp)*0.5e0_wp
            else
                ach1 = potcm_a1c1h+potcm_b1c1h
            endif
            if( arga2 .lt. 19.0e0_wp )then
                ach2 = potcm_a1c2h+potcm_b1c2h*(tanh(arga2)+1.0e0_wp)*0.5e0_wp
            else
                ach2 = potcm_a1c2h+potcm_b1c2h
            endif
            !  calculate singlet: e1, triplet: e3 energies and vq and vj
            !  terms for each bond
            e1    = potcm_d1cb*(exp(-2.0e0_wp*potcm_acb*(bonds3_rc1b-potcm_r0cb))-2.0e0_wp*exp(-potcm_acb* &
                & (bonds3_rc1b- potcm_r0cb)))
            e3    = potcm_d3cb*(exp(-2.0e0_wp*potcm_acb*(bonds3_rc1b-potcm_r0cb))+2.0e0_wp*exp(-potcm_acb* &
                & (bonds3_rc1b- potcm_r0cb)))
            vqc1b = (e1+e3)*0.5e0_wp
            vjc1b = (e1-e3)*0.5e0_wp
            e1    = potcm_d1cb*(exp(-2.0e0_wp*potcm_acb*(bonds3_rc2b-potcm_r0cb))-2.0e0_wp*exp(-potcm_acb* &
                & (bonds3_rc2b- potcm_r0cb)))
            e3    = potcm_d3cb*(exp(-2.0e0_wp*potcm_acb*(bonds3_rc2b-potcm_r0cb))+2.0e0_wp*exp(-potcm_acb* &
                & (bonds3_rc2b- potcm_r0cb)))
            vqc2b = (e1+e3)*0.5e0_wp
            vjc2b = (e1-e3)*0.5e0_wp
            do i=1, 3
                e1       = potcm_d1c1h*(exp(-2.0e0_wp*ach1*(bonds2_rc1h(i)-bonds_r0c1h)) -2.0e0_wp*exp(-ach1* &
                    & (bonds2_rc1h(i)-bonds_r0c1h)))
                e3       = potcm_d3c1h*(exp(-2.0e0_wp*ach1*(bonds2_rc1h(i)-bonds_r0c1h)) +2.0e0_wp*exp(-ach1* &
                    & (bonds2_rc1h(i)-bonds_r0c1h)))
                vqc1h(i) = (e1+e3)*0.5e0_wp
                vjc1h(i) = (e1-e3)*0.5e0_wp
                e1       = potcm_d1c2h*(exp(-2.0e0_wp*ach2*(bonds2_rc2h(i)-bonds_r0c2h)) -2.0e0_wp*exp(-ach2* &
                    & (bonds2_rc2h(i)-bonds_r0c2h)))
                e3       = potcm_d3c2h*(exp(-2.0e0_wp*ach2*(bonds2_rc2h(i)-bonds_r0c2h)) +2.0e0_wp*exp(-ach2* &
                    & (bonds2_rc2h(i)-bonds_r0c2h)))
                vqc2h(i) = (e1+e3)*0.5e0_wp
                vjc2h(i) = (e1-e3)*0.5e0_wp
            enddo
            do i=1, 6
                e1      = potcm_d1hh*(exp(-2.0e0_wp*potcm_ahh*(bonds2_rbh(i)-potcm_r0hh)) -2.0e0_wp*exp(-potcm_ahh* &
                    & (bonds2_rbh(i)-potcm_r0hh)))
                e3      = potcm_d3hh*(exp(-2.0e0_wp*potcm_ahh*(bonds2_rbh(i)-potcm_r0hh)) +2.0e0_wp*exp(-potcm_ahh* &
                    & (bonds2_rbh(i)-potcm_r0hh)))
                vqbh(i) = (e1+e3)*0.5e0_wp
                vjbh(i) = (e1-e3)*0.5e0_wp
            enddo
            !  calculate 3 body potential
            do i=1, 3
                vq(i) = vqc1h(i)+vqc1b+vqbh(i)
                vj(i) = -sqrt(((vjc1h(i)-vjc1b)**2+(vjc1b-vjbh(i))**2 +(vjbh(i)-vjc1h(i))**2)*0.5e0_wp)
                vstr  = vstr+vq(i)+vj(i)
            enddo
            do i=1, 3
                vq(i) = vqc2h(i)+vqc2b+vqbh(i+3)
                vj(i) = -sqrt(((vjc2h(i)-vjc2b)**2+(vjc2b-vjbh(i+3))**2 +(vjbh(i+3)-vjc2h(i))**2)*0.5e0_wp)
                vstr  = vstr+vq(i)+vj(i)
            enddo
            ! CRD 2012 potencial morse para el enlace C-C
            dt      = (potcm_rcc-potcm_r0cc)
            expterm = exp(-potcm_acc*dt)
            vcc     = potcm_d1cc*(1.0e0_wp-expterm)**2.0e0_wp
            vstr    = vstr + vcc
        end subroutine stretch

        !> @brief Evaluate opbend for the C2H6O_VBMM potential.
        subroutine opbend( vop )
            !  subroutine calculates symmetrized vop potential and derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: sum2, sum4, sumd2(8), sumd4(8), vop
            integer  :: i, j, k, l
            vop = 0.0e0_wp
            !  calculate force constants and their derivatives
            call opforce
            !  calculate out-of-plane angle and derivatives
            ! CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C2
            i = 4
            j = i+1
            if( j .gt. 4 ) j=j-4
            k = j+1
            if( k .gt. 4 ) k=k-4
            l = k+1
            if( l .gt. 4 ) l=l-4
            call calcdelta1(i, j, k, l, sum2, sum4)
            sumd2(i) = sum2
            sumd4(i) = sum4
            vop      = vop+delta1_fdelta(i)*sumd2(i)+delta1_hdelta(i)*sumd4(i)
            ! CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C1
            i = 8
            j = i+1
            if( j .gt. 8 ) j=j-4
            k = j+1
            if( k .gt. 8 ) k=k-4
            l = k+1
            if( l .gt. 8 ) l=l-4
            call calcdelta2(i, j, k, l, sum2, sum4)
            sumd2(i) = sum2
            sumd4(i) = sum4
            vop      = vop+delta1_fdelta(i)*sumd2(i)+delta1_hdelta(i)*sumd4(i)
        end subroutine opbend

        !> @brief Evaluate ipbend for the C2H6O_VBMM potential.
        subroutine ipbend( vip )
            !  subroutine calculates symmetrised in plane bend term
            !  and its derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: costh(8,8), dth(8,8), theta(8,8), vip
            integer  :: i, j
            !  initialise
            vip = 0.0e0_wp
            !  calculate force constants: force1_fk0(i,j), force1_f1(i)
            !  and derivatives wrt rch(k) and bonds2_rbh(k): force1_dfdc(i,j,k), force1_dfdh(i,j,k)
            call ipforce
            !  calculate theta(i,j) and in plane bend potential
            ! bucle para el c1
            do i=1, 3
                do j=i+1, 3
                    costh(i,j) = coords1_tc1h(i,1)*coords1_tc1h(j,1)+coords1_tc1h(i,2)*coords1_tc1h(j,2) + &
                        & coords1_tc1h(i,3)*coords1_tc1h(j,3)
                    costh(i,j) = costh(i,j)/bonds2_rc1h(i)/bonds2_rc1h(j)
                    theta(i,j) = acos(costh(i,j))
                    dth(i,j)   = theta(i,j)-angles_theta0(i,j)
                    vip        = vip+0.5e0_wp*force1_fk0(i,j)*force1_f1(i)*force1_f1(j)*dth(i,j)**2
                enddo
                j = 4
                costh(i,j) = coords1_tc1h(i,1)*coords2_tcc(1)+coords1_tc1h(i,2)*coords2_tcc(2) + coords1_tc1h(i,3)* &
                    & coords2_tcc(3)
                costh(i,j) = costh(i,j)/bonds2_rc1h(i)/potcm_rcc
                theta(i,j) = acos(costh(i,j))
                dth(i,j)   = theta(i,j)-angles_theta0(i,j)
                vip        = vip+0.5e0_wp*force1_fk0(i,j)*force1_f1(i)*force1_f1(j)*dth(i,j)**2
            enddo
            ! bucle para el c2 para los theta 5,6,7,8
            do i=1, 3
                do j=i+1, 3
                    costh(i+4,j+4) = coords1_tc2h(i,1)*coords1_tc2h(j,1)+coords1_tc2h(i,2)*coords1_tc2h(j,2) + &
                        & coords1_tc2h(i,3)*coords1_tc2h(j,3)
                    costh(i+4,j+4) = costh(i+4,j+4)/bonds2_rc2h(i)/bonds2_rc2h(j)
                    theta(i+4,j+4) = acos(costh(i+4,j+4))
                    dth(i+4,j+4)   = theta(i+4,j+4)-angles_theta0(i+4,j+4)
                    vip            = vip+0.5e0_wp*force1_fk0(i+4,j+4)*force1_f1(i+4)*force1_f1(j+4)*dth(i+4,j+4)**2
                enddo
                j = 4
                ! modificacion para corregir los angulos del c2
                !    *                       +coords1_tc2h(i,3)*coords2_tcc(3)
                costh(i+4,j+4) = coords1_tc2h(i,1)*(-coords2_tcc(1))+coords1_tc2h(i,2)*(-coords2_tcc(2)) + &
                    & coords1_tc2h(i,3)*(-coords2_tcc(3))
                costh(i+4,j+4) = costh(i+4,j+4)/bonds2_rc2h(i)/potcm_rcc
                theta(i+4,j+4) = acos(costh(i+4,j+4))
                dth(i+4,j+4)   = theta(i+4,j+4)-angles_theta0(i+4,j+4)
                vip            = vip+0.5e0_wp*force1_fk0(i+4,j+4)*force1_f1(i+4)*force1_f1(j+4)*dth(i+4,j+4)**2
            enddo
        end subroutine ipbend

        !> @brief Evaluate calcdelta1 for the C2H6O_VBMM potential.
        subroutine calcdelta1( i, j, k, l, sum2, sum4 )
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3), argd(8), axb(3), b(3), c(8,3), delta(8), sum2, sum4
            integer  :: i, ii, in(3), ind, j, k, l
            real(wp) :: norma
            !  initialise
            sum2 = 0.0e0_wp; sum4 = 0.0e0_wp
            !  set j,k,l indices
            in(1) = j
            in(2) = k
            in(3) = l
            ! CRD 2005 incluimos bucles if para asegurar que siempre
            ! que una de los valores j, k o l sea 4 las coordenadas
            ! corresponden al atomo nuevo.
            !  vector a is rk-rj, vector b is rl-rj
            if( j .eq. 4 )then
                do ind=1, 3
                    a(ind) = qpdot_q(ndx_nh(k,ind))-qpdot_q(ndx2_nc2(ind))
                    b(ind) = qpdot_q(ndx_nh(l,ind))-qpdot_q(ndx2_nc2(ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            if( k .eq. 4 )then
                do ind=1, 3
                    a(ind) = qpdot_q(ndx2_nc2(ind))-qpdot_q(ndx_nh(j,ind))
                    b(ind) = qpdot_q(ndx_nh(l,ind))-qpdot_q(ndx_nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            if( l .eq. 4 )then
                do ind=1, 3
                    a(ind) = qpdot_q(ndx_nh(k,ind))-qpdot_q(ndx_nh(j,ind))
                    b(ind) = qpdot_q(ndx2_nc2(ind))-qpdot_q(ndx_nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                ! modificacion para que hiciese vxu
                norma = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma = sqrt(norma)
            else
            endif
            ! CRD 2012 eliminamos porque seria la salida del C2
            if( k .lt. 4.and.j .lt. 4.and.l .lt. 4 )then
                do ind=1, 3
                    a(ind) = qpdot_q(ndx_nh(k,ind))-qpdot_q(ndx_nh(j,ind))
                    b(ind) = qpdot_q(ndx_nh(l,ind))-qpdot_q(ndx_nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            !  c is position vector of h(ii): calculate c(j),c(k),c(l)
            do ii=1, 3
                if( in(ii) .eq. 4 )then
                    do ind=1, 3
                        c(in(ii),ind) = -coords2_tcc(ind)/potcm_rcc
                    enddo
                else
                    do ind=1, 3
                        c(in(ii),ind) = -coords1_tc1h(in(ii),ind)/bonds2_rc1h(in(ii))
                    enddo
                endif
            enddo
            !  argd is the dot product axb dot c
            do ii=1, 3
                argd(in(ii)) = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                argd(in(ii)) = argd(in(ii))/norma
                ! Start of corrections by EMN
                ! The two lines below are added to avoid problems when acos is
                ! calculated
                if( argd(in(ii)) > 1.e0_wp ) argd(in(ii))  = 1
                if( argd(in(ii)) < -1.e0_wp ) argd(in(ii)) = -1
                ! End of corrections by EMN
                delta(in(ii)) = acos(argd(in(ii)))-angles_theta0(i,in(ii))
                sum2          = sum2+delta(in(ii))**2
                sum4          = sum4+delta(in(ii))**4
            enddo
        end subroutine calcdelta1

        !> @brief Evaluate calcdelta2 for the C2H6O_VBMM potential.
        subroutine calcdelta2( i, j, k, l, sum2, sum4 )
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3), argd(8), axb(3), b(3), c(8,3), delta(8), sum2, sum4
            integer  :: i, ii, in(3), ind, j, k, l
            real(wp) :: norma
            !  initialise
            sum2 = 0.0e0_wp; sum4 = 0.0e0_wp
            !  set j,k,l indices
            in(1) = j
            in(2) = k
            in(3) = l
            ! CRD 2005 incluimos bucles if para asegurar que siempre
            ! que una de los valores j, k o l sea 4 las coordenadas
            ! corresponden al atomo nuevo.
            !  vector a is rk-rj, vector b is rl-rj
            if( j .eq. 8 )then
                do ind=1, 3
                    a(ind) = qpdot_q(ndx_nh(k-1,ind))-qpdot_q(ndx2_nc1(ind))
                    b(ind) = qpdot_q(ndx_nh(l-1,ind))-qpdot_q(ndx2_nc1(ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            if( k .eq. 8 )then
                do ind=1, 3
                    a(ind) = qpdot_q(ndx2_nc1(ind))-qpdot_q(ndx_nh(j-1,ind))
                    b(ind) = qpdot_q(ndx_nh(l-1,ind))-qpdot_q(ndx_nh(j-1,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            if( l .eq. 8 )then
                do ind=1, 3
                    a(ind) = qpdot_q(ndx_nh(k-1,ind))-qpdot_q(ndx_nh(j-1,ind))
                    b(ind) = qpdot_q(ndx2_nc1(ind))-qpdot_q(ndx_nh(j-1,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            ! CRD 2012 eliminamos porque seria la salida del C1
            if( k .lt. 8.and.j .lt. 8.and.l .lt. 8 )then
                do ind=1, 3
                    a(ind) = qpdot_q(ndx_nh(k-1,ind))-qpdot_q(ndx_nh(j-1,ind))
                    b(ind) = qpdot_q(ndx_nh(l-1,ind))-qpdot_q(ndx_nh(j-1,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            !  c is position vector of h(ii): calculate c(j),c(k),c(l)
            do ii=1, 3
                if( in(ii) .eq. 4 )then
                    do ind=1, 3
                        c(in(ii),ind) = -coords2_tcc(ind)/potcm_rcc
                    enddo
                else
                    do ind=1, 3
                        c(in(ii),ind) = -coords1_tc2h(in(ii)-4,ind)/bonds2_rc2h(in(ii)-4)
                    enddo
                endif
            enddo
            !  argd is the dot product axb dot c
            do ii=1, 3
                argd(in(ii)) = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                argd(in(ii)) = argd(in(ii))/norma
                ! Start of corrections by EMN
                ! The two lines below are added to avoid problems when acos is
                ! calculated
                if( argd(in(ii)) > 1.e0_wp ) argd(in(ii))  = 1
                if( argd(in(ii)) < -1.e0_wp ) argd(in(ii)) = -1
                ! End of corrections by EMN
                delta(in(ii)) = acos(argd(in(ii)))-angles_theta0(i,in(ii))
                sum2          = sum2+delta(in(ii))**2
                sum4          = sum4+delta(in(ii))**4
            enddo
        end subroutine calcdelta2

        !> @brief Evaluate opforce for the C2H6O_VBMM potential.
        subroutine opforce
            !  calculates the out-of-plane bending force constants
            !  and their derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: switch(8)
            integer  :: i
            !  calculate switching functions:
            switch(1) = (1.0e0_wp-op1_s3(1))*op1_s3(2)*op1_s3(3)*op1_s3(4)
            switch(2) = (1.0e0_wp-op1_s3(2))*op1_s3(3)*op1_s3(4)*op1_s3(1)
            switch(3) = (1.0e0_wp-op1_s3(3))*op1_s3(4)*op1_s3(1)*op1_s3(2)
            switch(4) = (1.0e0_wp-op1_s3(4))*op1_s3(1)*op1_s3(2)*op1_s3(3)
            switch(5) = (1.0e0_wp-op1_s3(5))*op1_s3(6)*op1_s3(7)*op1_s3(8)
            switch(6) = (1.0e0_wp-op1_s3(6))*op1_s3(7)*op1_s3(8)*op1_s3(5)
            switch(7) = (1.0e0_wp-op1_s3(7))*op1_s3(8)*op1_s3(5)*op1_s3(6)
            switch(8) = (1.0e0_wp-op1_s3(8))*op1_s3(5)*op1_s3(6)*op1_s3(7)
            !  calculate the force constants
            do i=1, 8
                delta1_fdelta(i) = switch(i)*potcm_fch3
                delta1_hdelta(i) = switch(i)*potcm_hch3
            enddo
        end subroutine opforce

        !> @brief Evaluate ipforce for the C2H6O_VBMM potential.
        subroutine ipforce
            !  calculates the symmetrised in plane bend force constants and
            !  all partial derivatives involving them
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a1, a2, arga1, arga2, f0, f2
            integer  :: i
            !  set force constant at asymptotes
            f0              = potcm_fkinf+potcm_ak
            f2              = potcm_fkinf
            force1_fk0(1,2) = f0+f0*(ip1_s1(1)*ip1_s1(2)-1.0e0_wp)+(f0-f2)*(ip1_s2(3)*ip1_s2(4)-1.0e0_wp)
            force1_fk0(1,3) = f0+f0*(ip1_s1(1)*ip1_s1(3)-1.0e0_wp)+(f0-f2)*(ip1_s2(2)*ip1_s2(4)-1.0e0_wp)
            force1_fk0(1,4) = f0+f0*(ip1_s1(1)*ip1_s1(4)-1.0e0_wp)+(f0-f2)*(ip1_s2(2)*ip1_s2(3)-1.0e0_wp)
            force1_fk0(2,3) = f0+f0*(ip1_s1(2)*ip1_s1(3)-1.0e0_wp)+(f0-f2)*(ip1_s2(1)*ip1_s2(4)-1.0e0_wp)
            force1_fk0(2,4) = f0+f0*(ip1_s1(2)*ip1_s1(4)-1.0e0_wp)+(f0-f2)*(ip1_s2(1)*ip1_s2(3)-1.0e0_wp)
            force1_fk0(3,4) = f0+f0*(ip1_s1(3)*ip1_s1(4)-1.0e0_wp)+(f0-f2)*(ip1_s2(1)*ip1_s2(2)-1.0e0_wp)
            force1_fk0(5,6) = f0+f0*(ip1_s1(5)*ip1_s1(6)-1.0e0_wp)+(f0-f2)*(ip1_s2(7)*ip1_s2(8)-1.0e0_wp)
            force1_fk0(5,7) = f0+f0*(ip1_s1(5)*ip1_s1(7)-1.0e0_wp)+(f0-f2)*(ip1_s2(6)*ip1_s2(8)-1.0e0_wp)
            force1_fk0(5,8) = f0+f0*(ip1_s1(5)*ip1_s1(8)-1.0e0_wp)+(f0-f2)*(ip1_s2(6)*ip1_s2(7)-1.0e0_wp)
            force1_fk0(6,7) = f0+f0*(ip1_s1(6)*ip1_s1(7)-1.0e0_wp)+(f0-f2)*(ip1_s2(5)*ip1_s2(8)-1.0e0_wp)
            force1_fk0(6,8) = f0+f0*(ip1_s1(6)*ip1_s1(8)-1.0e0_wp)+(f0-f2)*(ip1_s2(5)*ip1_s2(7)-1.0e0_wp)
            force1_fk0(7,8) = f0+f0*(ip1_s1(7)*ip1_s1(8)-1.0e0_wp)+(f0-f2)*(ip1_s2(5)*ip1_s2(6)-1.0e0_wp)
            ! force1_f1 para los h del c1
            do i=1, 3
                !  calc derivatives of force1_fk0 wrt each of the rch(i) bonds
                !  calculate the terms force1_f1(i)
                arga1        = potcm_aa1*bonds2_rbh(i)*bonds2_rbh(i)
                arga2        = potcm_aa4*(bonds2_rbh(i)-potcm_r0hh)*(bonds2_rbh(i)-potcm_r0hh)
                a1           = 1.0e0_wp-exp(-arga1)
                a2           = potcm_aa2+potcm_aa3*exp(-arga2)
                force1_f1(i) = a1*exp(-a2*(bonds2_rc1h(i)-bonds_r0c1h)**2)
            enddo
            ! force1_f1 para el c2
            ! CRD 2012 eliminamos el calculo de potcm_a1c y potcm_a2c pq es constante
            ! lo incluimos en el CONST
            force1_f1(4) = potcm_a1c*exp(-potcm_a2c*(potcm_rcc-potcm_r0cc)**2)
            force1_f1(8) = force1_f1(4)
            ! force1_f1 para los h del c2
            do i=5, 7
                !  calc derivatives of force1_fk0 wrt each of the rch(i) bonds
                !  calculate the terms force1_f1(i)
                arga1        = potcm_aa1*bonds2_rbh(i-1)*bonds2_rbh(i-1)
                arga2        = potcm_aa4*(bonds2_rbh(i-1)-potcm_r0hh)*(bonds2_rbh(i-1)-potcm_r0hh)
                a1           = 1.0e0_wp-exp(-arga1)
                a2           = potcm_aa2+potcm_aa3*exp(-arga2)
                force1_f1(i) = a1*exp(-a2*(bonds2_rc2h(i-4)-bonds_r0c2h)**2)
            enddo
            ! force1_f1 para el c2
            ! CRD 2012 eliminamos el calculo de potcm_a1c y potcm_a2c pq es constante
            ! lo incluimos en el CONST
        end subroutine ipforce

        !> @brief Evaluate switchf for the C2H6O_VBMM potential.
        subroutine switchf
            !  calculates switching functions: op1_s3,switch1_sphi,switch1_stheta
            !  and their derivatives op1_ds3,switch1_dsphi,switch1_dstheta
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: argmax, args1, args2, args3, argsc1, argsc2, argsc3, argsphi, argsphic, argstheta, argsthetac
            integer  :: i
            !  nb remember that integration units are:
            !  energy in   1.0d+05 j/mol
            !  time in     1.0d-14 s
            ! original     potcm_a1s=1.5313681d-7
            !      potcm_a1s=1.5313681d-8
            !      potcm_b1s=-4.6696246d0
            ! original      potcm_a2s=1.0147402d-7
            !      potcm_a2s=1.0147402d-8
            !      potcm_b2s=-12.363798d0
            !  use double precision criterion:
            !  tanh(19.0d0)=1.0d0
            argmax = 19.0e0_wp
            !  calculate ip1_s1 and ip1_ds1
            do i=1, 3
                args1 = potcm_a1s*(bonds2_rc1h(i)-bonds_r0c1h)*(bonds2_rc1h(i)-potcm_b1s)**8
                if( args1 .lt. argmax )then
                    ip1_s1(i) = 1.0e0_wp-tanh(args1)
                else
                    ip1_s1(i) = 0.0e0_wp
                endif
            enddo
            do i=1, 3
                args1 = potcm_a1s*(bonds2_rc2h(i)-bonds_r0c2h)*(bonds2_rc2h(i)-potcm_b1s)**8
                if( args1 .lt. argmax )then
                    ip1_s1(i+4) = 1.0e0_wp-tanh(args1)
                else
                    ip1_s1(i+4) = 0.0e0_wp
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            argsc1 = potcm_a1scc*(potcm_rcc-potcm_r0cc)*(potcm_rcc-potcm_b1scc)**8
            if( argsc1 .lt. argmax )then
                ip1_s1(4) = 1.0e0_wp-tanh(argsc1)
                ip1_s1(8) = ip1_s1(4)
            else
                ip1_s1(4) = 0.0e0_wp
                ip1_s1(8) = ip1_s1(4)
            endif
            !  calculate ip1_s2 and ip1_ds2
            do i=1, 3
                args2 = potcm_a2s*(bonds2_rc1h(i)-bonds_r0c1h)*(bonds2_rc1h(i)-potcm_b2s)**6
                if( args2 .lt. argmax )then
                    ip1_s2(i) = 1.0e0_wp-tanh(args2)
                else
                    ip1_s2(i) = 0.0e0_wp
                endif
            enddo
            do i=1, 3
                args2 = potcm_a2s*(bonds2_rc2h(i)-bonds_r0c2h)*(bonds2_rc2h(i)-potcm_b2s)**6
                if( args2 .lt. argmax )then
                    ip1_s2(i+4) = 1.0e0_wp-tanh(args2)
                else
                    ip1_s2(i+4) = 0.0e0_wp
                endif
            enddo
            ! CRD 2005 ip1_s2(4)for the new atom
            ! CRD 2012 todos los S de los carbonos son 1
            argsc2 = potcm_a2scc*(potcm_rcc-potcm_r0cc)*(potcm_rcc-potcm_b2scc)**6
            if( argsc2 .lt. argmax )then
                ip1_s2(4) = 1.0e0_wp-tanh(argsc2)
                ip1_s2(8) = ip1_s2(4)
            else
                ip1_s2(4) = 0.0e0_wp
                ip1_s2(8) = ip1_s2(4)
            endif
            !  calculate op1_s3 and op1_ds3
            do i=1, 3
                args3 = potcm_a3s*(bonds2_rc1h(i)-bonds_r0c1h)*(bonds2_rc1h(i)-potcm_b3s)**2
                if( args3 .lt. argmax )then
                    op1_s3(i) = 1.0e0_wp-tanh(args3)
                else
                    op1_s3(i) = 0.0e0_wp
                endif
            enddo
            do i=1, 3
                args3 = potcm_a3s*(bonds2_rc2h(i)-bonds_r0c2h)*(bonds2_rc2h(i)-potcm_b3s)**2
                if( args3 .lt. argmax )then
                    op1_s3(i+4) = 1.0e0_wp-tanh(args3)
                else
                    op1_s3(i+4) = 0.0e0_wp
                endif
            enddo
            ! CRD 2005 op1_s3(4) for the new atom
            ! CRD 2012 todos los S de los carbonos son 1
            argsc3 = potcm_a3scc*(potcm_rcc-potcm_r0cc)*(potcm_rcc-potcm_b3scc)**2
            if( argsc3 .lt. argmax )then
                op1_s3(4) = 1.0e0_wp-tanh(argsc3)
                op1_s3(8) = op1_s3(4)
            else
                op1_s3(4) = 0.0e0_wp
                op1_s3(8) = op1_s3(4)
            endif
            ! JEG 2019. cambios por problemas en RPMD de F,Cl,H + C2H6!!
            op1_s3(4) = 1.0e0_wp
            op1_s3(8) = 1.0e0_wp
            ! JEG 2019
            !  calculate switch1_sphi and switch1_dsphi
            !  condition here is on the bondlength rch(i)
            !  st argsphi is lt approx 19.0d0
            do i=1, 3
                if( bonds2_rc1h(i) .lt. 3.8e0_wp )then
                    argsphi         = potcm_aphi*(bonds2_rc1h(i)-bonds_r0c1h)*exp(potcm_bphi*(bonds2_rc1h(i)- &
                        & potcm_cphi)**3)
                    switch1_sphi(i) = 1.0e0_wp-tanh(argsphi)
                else
                    switch1_sphi(i) = 0.0e0_wp
                endif
            enddo
            do i=1, 3
                if( bonds2_rc2h(i) .lt. 3.8e0_wp )then
                    argsphi           = potcm_aphi*(bonds2_rc2h(i)-bonds_r0c2h)*exp(potcm_bphi*(bonds2_rc2h(i)- &
                        & potcm_cphi)**3)
                    switch1_sphi(i+4) = 1.0e0_wp-tanh(argsphi)
                else
                    switch1_sphi(i+4) = 0.0e0_wp
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            if( potcm_rcc .lt. 3.8e0_wp )then
                argsphic        = potcm_aphi*(potcm_rcc-potcm_r0cc)*exp(potcm_bphi*(potcm_rcc-potcm_cphi)**3)
                switch1_sphi(4) = 1.0e0_wp-tanh(argsphic)
                switch1_sphi(8) = switch1_sphi(4)
            else
                switch1_sphi(4) = 0.0e0_wp
                switch1_sphi(8) = switch1_sphi(4)
            endif
            !  calculate switch1_stheta and switch1_dstheta
            do i=1, 3
                if( bonds2_rc1h(i) .lt. 3.8e0_wp )then
                    argstheta = potcm_atheta*(bonds2_rc1h(i)-bonds_r0c1h)*exp(potcm_btheta*(bonds2_rc1h(i)- &
                        & potcm_ctheta)**3)
                    switch1_stheta(i) = 1.0e0_wp-tanh(argstheta)
                else
                    switch1_stheta(i) = 0.0e0_wp
                endif
            enddo
            do i=1, 3
                if( bonds2_rc2h(i) .lt. 3.8e0_wp )then
                    argstheta   = potcm_atheta*(bonds2_rc2h(i)-bonds_r0c2h)*exp(potcm_btheta*(bonds2_rc2h(i)- &
                        & potcm_ctheta)**3)
                    switch1_stheta(i+4) = 1.0e0_wp-tanh(argstheta)
                else
                    switch1_stheta(i+4) = 0.0e0_wp
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            if( potcm_rcc .lt. 3.8e0_wp )then
                argsthetac        = potcm_atheta*(potcm_rcc-potcm_r0cc)*exp(potcm_btheta*(potcm_rcc-potcm_ctheta)**3)
                switch1_stheta(4) = 1.0e0_wp-tanh(argsthetac)
                switch1_stheta(8) = switch1_stheta(4)
            else
                switch1_stheta(4) = 0.0e0_wp
                switch1_stheta(8) = switch1_stheta(4)
            endif
        end subroutine switchf

        !> @brief Evaluate prepare_local_context for the C2H6O_VBMM potential.
        subroutine prepare_local_context
            !   N3TMMN = 3 * infocm_natoms
            !   infocm_natoms = the number of atoms represented by this potential function
            !   The variable N3TMMN is the minimum value of N3TM allowed to be
            !   passed by the calling routine for the number of cartesian
            !   coordinates needed to represent the full system represented by this
            !   potential energy surface routine.
            !   N3TM must be greater than or equal to N3TMMN.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: d1bc, d3bc, fact1, fact2, vdwc2
            integer  :: i, icount, ind
            !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
            if( infocm_natoms .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !      infocm_ref(1)='J. Espinosa-Garcia and J. C. Corchado'
            infocm_indexes(1) = 6
            infocm_indexes(2) = 6
            infocm_indexes(3) = 1
            infocm_indexes(4) = 1
            infocm_indexes(5) = 1
            infocm_indexes(6) = 1
            infocm_indexes(7) = 1
            infocm_indexes(8) = 1
            infocm_indexes(9) = 17
            infocm_irctnt     = 9
            call ancvrt
            !  calculate infocm_indexes for coordinates
            do ind=1, 3
                icount        = ind-3
                ndx2_nc1(ind) = 3*potcm_nnc1+icount
                ndx2_nc2(ind) = 3*potcm_nnc2+icount
                ndx_nhb(ind)  = 3*potcm_nnb+icount
                do i=1, 6
                    ndx_nh(i,ind) = 3*ndx2_nnh(i)+icount
                enddo
            enddo
            !  convert to appropriate units:
            !  From:
            !  kcal/mol
            !  mdyn A-1        -> 1.0d+05 j/mol...
            !  mdyn A rad-1
            !  To
            !  energy   in 1.0d+05 j/mol
            !  distance in 1.0d-10 m
            !  angles   in radians
            !  mass     in amu
            fact1       = 0.041840e0_wp
            fact2       = 6.022045e0_wp
            stret1_d1ch = stret1_d1ch*fact1
            stret1_d3ch = stret1_d3ch*fact1
            potcm_d1c1h = potcm_d1c1h*fact1
            potcm_d3c1h = potcm_d3c1h*fact1
            potcm_d1c2h = potcm_d1c2h*fact1
            potcm_d3c2h = potcm_d3c2h*fact1
            potcm_d1cc  = potcm_d1cc*fact1
            d1bc        = d1bc*fact1
            d3bc        = d3bc*fact1
            potcm_d1cb  = potcm_d1cb*fact1
            potcm_d3cb  = potcm_d3cb*fact1
            potcm_d1hh  = potcm_d1hh*fact1
            potcm_d3hh  = potcm_d3hh*fact1
            potcm_fch3  = potcm_fch3*fact2
            potcm_hch3  = potcm_hch3*fact2
            potcm_fkinf = potcm_fkinf*fact2
            potcm_ak    = potcm_ak*fact2
            potcm_a1s   = potcm_a1s*1.e-8_wp
            potcm_a2s   = potcm_a2s*1.e-8_wp
            potcm_a1scc = potcm_a1scc*1.e-8_wp
            potcm_a2scc = potcm_a2scc*1.e-8_wp
            vdwc2       = vdwc2*1.e+4_wp
        end subroutine prepare_local_context

        !> @brief Evaluate POTINFO for the C2H6O_VBMM potential.
        subroutine potinfo
            implicit none
            integer            :: i, j, k, kmax
            integer, parameter :: isurf = 5
            do i =1, 5
            enddo
            kmax = 0
            do i = 1, isurf+1
                do j = 1, isurf+1
                    if( usricm_nasurf(i,j) .ne. 0.and.kmax .lt. max(i, j) ) kmax = max(i, j)
                enddo
            enddo
            if( kmax-1 .gt. infocm_msurf )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
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
            if( usricm_nder .gt. infocm_mder )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            if( usricm_nflag(19) .eq. 1 )then
            endif
        end subroutine potinfo

        !> @brief Evaluate ANCVRT for the C2H6O_VBMM potential.
        subroutine ancvrt
            use chempotpy_atomic_elements, only: all_atomic_numbers_valid
            use chempotpy_kinds,           only: wp
            implicit none
            real(wp), parameter :: htomillh = 1000.e0_wp
            real(wp), parameter :: htoev    = 27.2113834e0_wp
            real(wp), parameter :: htokcal  = hartree_to_kilocalorie_per_mole
            real(wp), parameter :: htowave  = 219474.631e0_wp
            real(wp), parameter :: htokj    = 2625.49962e0_wp
            real(wp), parameter :: bohr_a   = .5291772083e0_wp
            integer             :: inu, isum
            if( .not. all_atomic_numbers_valid(infocm_indexes(1:infocm_natoms)) )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'potential declares an invalid atomic number')
                return
            endif
            utilcm_cnvrtd  = 1.e0_wp
            utilcm_cnvrte  = 1.e0_wp
            utilcm_cnvrtde = 1.e0_wp
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
            do inu=1, size(usricm_nulbl)
                isum = isum+usricm_nulbl(inu)
            enddo
            utilcm_ireorder = 0
            if( isum .ne. 0 ) utilcm_ireorder = 1
        end subroutine ancvrt

        !> @brief Evaluate CARTOU for the C2H6O_VBMM potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            implicit none
            call convert_cartesian_input(usricm_cart, usricm_nulbl, infocm_natoms, utilcm_cnvrtd, utilcm_ireorder == 1, infocm_cartnu)
        end subroutine cartou

        !> @brief Evaluate CARTTOR for the C2H6O_VBMM potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(infocm_cartnu, infocm_natoms, infocm_icartr, pt1cm_r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the C2H6O_VBMM potential.
        subroutine eunitzero
            implicit none
            integer :: i, j
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

        !> @brief Evaluate RTOCART for the C2H6O_VBMM potential.
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
                    utilcm_dgscartnu(i,1) = 0.e0_wp
                    utilcm_dgscartnu(i,2) = 0.e0_wp
                    utilcm_dgscartnu(i,3) = 0.e0_wp
                    if( utilcm_ksdiag .ne. 0 )then
                        do j1=utilcm_ksdiag, utilcm_kediag
                            utilcm_descartnu(i,1,j1) = 0.e0_wp
                            utilcm_descartnu(i,2,j1) = 0.e0_wp
                            utilcm_descartnu(i,3,j1) = 0.e0_wp
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do j2=utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(i,1,j2) = 0.e0_wp
                            utilcm_dijcartnu(i,2,j2) = 0.e0_wp
                            utilcm_dijcartnu(i,3,j2) = 0.e0_wp
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
                            yes(i, j) = pt4cm_deesdr(i,j)/pt1cm_r(i)
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do k=utilcm_ksoffd, utilcm_keoffd
                            yij(i, k) = pt5cm_deijdr(i,k)/pt1cm_r(i)
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
                            utilcm_descartnu(1,k,j1) = term12*yes(1, j1) + term13*yes(3, j1)
                            utilcm_descartnu(2,k,j1) = -term12*yes(1, j1) + term23*yes(2, j1)
                            utilcm_descartnu(3,k,j1) = -term13*yes(3, j1) - term23*yes(2, j1)
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do j2=utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(1,k,j2) = term12*yij(1, j2) + term13*yij(3, j2)
                            utilcm_dijcartnu(2,k,j2) = -term12*yij(1, j2) + term23*yij(2, j2)
                            utilcm_dijcartnu(3,k,j2) = -term13*yij(3, j2) - term23*yij(2, j2)
                        enddo
                    endif
                enddo
            elseif( infocm_icartr .eq. 4 )then
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

        !> @brief Evaluate DEDCOU for the C2H6O_VBMM potential.
        subroutine dedcou
            implicit none
            integer :: i, j, k
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
end module chempotpy_c2h6o_vbmm_pes

!> @brief Expose C2H6O_VBMM through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_c2h6o_vbmm_pes, only: evaluate
    use chempotpy_kinds,          only: wp
    implicit none
    integer, parameter    :: natoms  = 9
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
