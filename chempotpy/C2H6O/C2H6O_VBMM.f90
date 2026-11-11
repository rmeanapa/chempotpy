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
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the C2H6O_VBMM potential energy surface for C2H6O.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_legacy_units, only: legacy_hartree_to_electronvolt
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 9
        integer, parameter            :: natom   = 25
        integer, parameter            :: isurf   = 5
        integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter            :: n3atom  = 75
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: angles_dtheta0(8,8,8)
        double precision              :: angles_theta0(8,8)
        double precision              :: bonds2_rbh(6)
        double precision              :: bonds2_rc1h(3)
        double precision              :: bonds2_rc2h(3)
        double precision              :: bonds2_rhh(3,6)
        double precision              :: bonds3_rc1b
        double precision              :: bonds3_rc2b
        double precision              :: bonds_r0c1h
        double precision              :: bonds_r0c2h
        double precision              :: coords1_tbh(6,3)
        double precision              :: coords1_tc1h(3,3)
        double precision              :: coords1_tc2h(3,3)
        double precision              :: coords1_tcb(3)
        double precision              :: coords2_tcc(3)
        double precision              :: coords2_thh(3,6,3)
        double precision              :: coords3_tc1b(3)
        double precision              :: coords3_tc2b(3)
        double precision              :: delta1_fdelta(8)
        double precision              :: delta1_hdelta(8)
        double precision              :: delta2_dfdelta(8,8)
        double precision              :: delta2_dhdelta(8,8)
        double precision              :: force1_dfdc(8,8,8)
        double precision              :: force1_dfdh(8,8,8)
        double precision              :: force1_f1(8)
        double precision              :: force1_fk0(8,8)
        double precision              :: infocm_cartnu(natom,3)
        integer                       :: infocm_icartr
        integer                       :: infocm_indexes(natom)
        integer                       :: infocm_irctnt
        integer                       :: infocm_mder
        integer                       :: infocm_msurf
        integer                       :: infocm_natoms
        character                     :: infocm_ref(5)
        double precision              :: ip1_ds1(8)
        double precision              :: ip1_ds2(8)
        double precision              :: ip1_s1(8)
        double precision              :: ip1_s2(8)
        integer                       :: ndx2_nc1(3)
        integer                       :: ndx2_nc2(3)
        integer                       :: ndx2_nnh(6)
        integer                       :: ndx_nc(3)
        integer                       :: ndx_nh(8,3)
        integer                       :: ndx_nhb(3)
        double precision              :: op1_ds3(8)
        double precision              :: op1_s3(8)
        double precision              :: potcm_a1c
        double precision              :: potcm_a1c1h
        double precision              :: potcm_a1c2h
        double precision              :: potcm_a1s
        double precision              :: potcm_a1scc
        double precision              :: potcm_a2c
        double precision              :: potcm_a2s
        double precision              :: potcm_a2scc
        double precision              :: potcm_a3s
        double precision              :: potcm_a3scc
        double precision              :: potcm_aa1
        double precision              :: potcm_aa2
        double precision              :: potcm_aa3
        double precision              :: potcm_aa4
        double precision              :: potcm_acb
        double precision              :: potcm_acc
        double precision              :: potcm_ahh
        double precision              :: potcm_ak
        double precision              :: potcm_aphi
        double precision              :: potcm_atheta
        double precision              :: potcm_b1c1h
        double precision              :: potcm_b1c2h
        double precision              :: potcm_b1s
        double precision              :: potcm_b1scc
        double precision              :: potcm_b2s
        double precision              :: potcm_b2scc
        double precision              :: potcm_b3s
        double precision              :: potcm_b3scc
        double precision              :: potcm_bk
        double precision              :: potcm_bphi
        double precision              :: potcm_btheta
        double precision              :: potcm_c1c1h
        double precision              :: potcm_c1c2h
        double precision              :: potcm_cphi
        double precision              :: potcm_ctheta
        double precision              :: potcm_d1c1h
        double precision              :: potcm_d1c2h
        double precision              :: potcm_d1cb
        double precision              :: potcm_d1cc
        double precision              :: potcm_d1hh
        double precision              :: potcm_d3c1h
        double precision              :: potcm_d3c2h
        double precision              :: potcm_d3cb
        double precision              :: potcm_d3hh
        double precision              :: potcm_fch3
        double precision              :: potcm_fkinf
        double precision              :: potcm_hch3
        integer                       :: potcm_nnb
        integer                       :: potcm_nnc1
        integer                       :: potcm_nnc2
        double precision              :: potcm_r0c1hp
        double precision              :: potcm_r0c1hr
        double precision              :: potcm_r0c2hp
        double precision              :: potcm_r0c2hr
        double precision              :: potcm_r0cb
        double precision              :: potcm_r0cc
        double precision              :: potcm_r0hh
        double precision              :: potcm_rcc
        double precision              :: potcm_v3
        double precision              :: potcm_w1
        double precision              :: potcm_w2
        double precision              :: potcm_w3
        double precision              :: potcm_w4
        double precision              :: potcm_w5
        double precision              :: potcm_w6
        double precision              :: pt1cm_degsdr(n3atom)
        double precision              :: pt1cm_engygs
        double precision              :: pt1cm_r(n3atom)
        double precision              :: pt3cm_ezero(isurf+1)
        double precision              :: pt4cm_deesdr(n3atom,isurf)
        double precision              :: pt4cm_engyes(isurf)
        double precision              :: pt5cm_deijdr(n3atom,jsurf)
        double precision              :: pt5cm_engyij(jsurf)
        double precision              :: qpdot_pdot(27)
        double precision              :: qpdot_q(27)
        double precision              :: stret1_a1cc
        double precision              :: stret1_b1cc
        double precision              :: stret1_c1cc
        double precision              :: stret1_d1ch
        double precision              :: stret1_d3ch
        double precision              :: switch1_dsphi(8)
        double precision              :: switch1_dstheta(8)
        double precision              :: switch1_sphi(8)
        double precision              :: switch1_stheta(8)
        double precision              :: thetahcch_thetahcch(9)
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
        integer                       :: iatom, idir, istate
        !initialize
        angles_dtheta0      = 0.d0; angles_theta0 = 0.d0; bonds2_rbh = 0.d0; bonds2_rc1h = 0.d0; bonds2_rc2h = 0.d0
        bonds2_rhh          = 0.d0; bonds3_rc1b = 0.d0; bonds3_rc2b = 0.d0; bonds_r0c1h = 0.d0; bonds_r0c2h = 0.d0
        coords1_tbh         = 0.d0; coords1_tc1h = 0.d0; coords1_tc2h = 0.d0; coords1_tcb = 0.d0; coords2_tcc = 0.d0
        coords2_thh         = 0.d0; coords3_tc1b = 0.d0; coords3_tc2b = 0.d0; delta1_fdelta = 0.d0; delta1_hdelta = 0.d0
        delta2_dfdelta      = 0.d0; delta2_dhdelta = 0.d0; force1_dfdc = 0.d0; force1_dfdh = 0.d0; force1_f1 = 0.d0
        force1_fk0          = 0.d0; infocm_cartnu = 0.d0; infocm_icartr = 0; infocm_indexes = 0; infocm_irctnt = 0
        infocm_mder         = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref          = ' '
        ip1_ds1             = 0.d0
        ip1_ds2             = 0.d0
        ip1_s1              = 0.d0
        ip1_s2              = 0.d0
        ndx2_nc1            = 0
        ndx2_nc2            = 0
        ndx2_nnh            = 0
        ndx_nc              = 0
        ndx_nh              = 0
        ndx_nhb             = 0
        op1_ds3             = 0.d0
        op1_s3              = 0.d0
        potcm_a1c           = 0.d0
        potcm_a1c1h         = 0.d0
        potcm_a1c2h         = 0.d0
        potcm_a1s           = 0.d0
        potcm_a1scc         = 0.d0
        potcm_a2c           = 0.d0
        potcm_a2s           = 0.d0
        potcm_a2scc         = 0.d0
        potcm_a3s           = 0.d0
        potcm_a3scc         = 0.d0
        potcm_aa1           = 0.d0
        potcm_aa2           = 0.d0
        potcm_aa3           = 0.d0
        potcm_aa4           = 0.d0
        potcm_acb           = 0.d0
        potcm_acc           = 0.d0
        potcm_ahh           = 0.d0
        potcm_ak            = 0.d0
        potcm_aphi          = 0.d0
        potcm_atheta        = 0.d0
        potcm_b1c1h         = 0.d0
        potcm_b1c2h         = 0.d0
        potcm_b1s           = 0.d0
        potcm_b1scc         = 0.d0
        potcm_b2s           = 0.d0
        potcm_b2scc         = 0.d0
        potcm_b3s           = 0.d0
        potcm_b3scc         = 0.d0
        potcm_bk            = 0.d0
        potcm_bphi          = 0.d0
        potcm_btheta        = 0.d0
        potcm_c1c1h         = 0.d0
        potcm_c1c2h         = 0.d0
        potcm_cphi          = 0.d0
        potcm_ctheta        = 0.d0
        potcm_d1c1h         = 0.d0
        potcm_d1c2h         = 0.d0
        potcm_d1cb          = 0.d0
        potcm_d1cc          = 0.d0
        potcm_d1hh          = 0.d0
        potcm_d3c1h         = 0.d0
        potcm_d3c2h         = 0.d0
        potcm_d3cb          = 0.d0
        potcm_d3hh          = 0.d0
        potcm_fch3          = 0.d0
        potcm_fkinf         = 0.d0
        potcm_hch3          = 0.d0
        potcm_nnb           = 0
        potcm_nnc1          = 0
        potcm_nnc2          = 0
        potcm_r0c1hp        = 0.d0
        potcm_r0c1hr        = 0.d0
        potcm_r0c2hp        = 0.d0
        potcm_r0c2hr        = 0.d0
        potcm_r0cb          = 0.d0
        potcm_r0cc          = 0.d0
        potcm_r0hh          = 0.d0
        potcm_rcc           = 0.d0
        potcm_v3            = 0.d0
        potcm_w1            = 0.d0
        potcm_w2            = 0.d0
        potcm_w3            = 0.d0
        potcm_w4            = 0.d0
        potcm_w5            = 0.d0
        potcm_w6            = 0.d0
        pt1cm_degsdr        = 0.d0
        pt1cm_engygs        = 0.d0
        pt1cm_r             = 0.d0
        pt3cm_ezero         = 0.d0
        pt4cm_deesdr        = 0.d0
        pt4cm_engyes        = 0.d0
        pt5cm_deijdr        = 0.d0
        pt5cm_engyij        = 0.d0
        qpdot_pdot          = 0.d0
        qpdot_q             = 0.d0
        stret1_a1cc         = 0.d0
        stret1_b1cc         = 0.d0
        stret1_c1cc         = 0.d0
        stret1_d1ch         = 0.d0
        stret1_d3ch         = 0.d0
        switch1_dsphi       = 0.d0
        switch1_dstheta     = 0.d0
        switch1_sphi        = 0.d0
        switch1_stheta      = 0.d0
        thetahcch_thetahcch = 0.d0
        usricm_anuzero      = 0.d0
        usricm_cart         = 0.d0
        usricm_nasurf       = 0
        usricm_nder         = 0
        usricm_nflag        = 0
        usricm_nulbl        = 0
        usrocm_descart      = 0.d0
        usrocm_dgscart      = 0.d0
        usrocm_dijcart      = 0.d0
        usrocm_pengyes      = 0.d0
        usrocm_pengygs      = 0.d0
        usrocm_pengyij      = 0.d0
        utilcm_cnvrtd       = 0.d0
        utilcm_cnvrtde      = 0.d0
        utilcm_cnvrte       = 0.d0
        utilcm_descartnu    = 0.d0
        utilcm_dgscartnu    = 0.d0
        utilcm_dijcartnu    = 0.d0
        utilcm_ireorder     = 0
        utilcm_kediag       = 0
        utilcm_keoffd       = 0
        utilcm_ksdiag       = 0
        utilcm_ksoffd       = 0
        usricm_nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0], shape(usricm_nasurf))
        usricm_nder    = 1
        usricm_nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero = 0.0d0
        infocm_icartr  = 1
        infocm_msurf   = 0
        infocm_mder    = 1
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 9
        potcm_nnc1     = 1
        potcm_nnc2     = 2
        potcm_nnb      = 9
        ndx2_nnh       = [3, 4, 5, 6, 7, 8]
        potcm_r0c1hr   = 1.08897d0
        potcm_r0c1hp   = 1.07577d0
        potcm_w3       = 1.00000d0
        potcm_w4       = 1.08897d0
        potcm_d1c1h    = 129.23000d0
        potcm_d3c1h    = 36.00400d0
        potcm_a1c1h    = 1.60300d0
        potcm_b1c1h    = 0.12500d0
        potcm_c1c1h    = 3.51404d0
        potcm_r0c2hr   = 1.08897d0
        potcm_r0c2hp   = 1.07577d0
        potcm_w5       = 1.00000d0
        potcm_w6       = 1.08897d0
        potcm_d1c2h    = 129.23000d0
        potcm_d3c2h    = 36.00400d0
        potcm_a1c2h    = 1.60300d0
        potcm_b1c2h    = 0.12500d0
        potcm_c1c2h    = 3.51404d0
        potcm_r0cc     = 1.50007d0
        potcm_d1cc     = 117.94900d0
        potcm_acc      = 1.93400d0
        potcm_a1c      = 0.90007d0
        potcm_a2c      = 0.94900d0
        potcm_r0hh     = 0.97101d0
        potcm_d1hh     = 127.45800d0
        potcm_d3hh     = 40.36400d0
        potcm_ahh      = 2.10400d0
        potcm_r0cb     = 1.77097d0
        potcm_d1cb     = 48.20900d0
        potcm_d3cb     = 15.66000d0
        potcm_acb      = 1.3290285d0
        potcm_a3s      = 1.8419147d0
        potcm_b3s      = 1.0889300d0
        potcm_aphi     = 1.0087903d0
        potcm_bphi     = 1.2006638d0
        potcm_cphi     = 1.0889937d0
        potcm_atheta   = 1.2078714d0
        potcm_btheta   = 0.3548859d0
        potcm_ctheta   = 1.8915497d0
        potcm_fch3     = 0.1357500d0
        potcm_hch3     = 0.2615000d0
        potcm_fkinf    = 0.2977000d0
        potcm_ak       = 0.3060000d0
        potcm_bk       = 0.2032000d0
        potcm_aa1      = 0.500952d0
        potcm_aa2      = 0.499963d0
        potcm_aa3      = 6.465953d0
        potcm_aa4      = 1.569977d0
        potcm_v3       = 13.730000d0
        potcm_w1       = 9.3548859d0
        potcm_w2       = 0.913970d0
        potcm_a1s      = 2.5313681d0
        potcm_b1s      = 1.088625d0
        potcm_a2s      = 2.0147402d0
        potcm_b2s      = 1.088980d0
        potcm_a1scc    = 2.5313681d0
        potcm_b1scc    = 1.550070d0
        potcm_a2scc    = 2.0147402d0
        potcm_b2scc    = 1.550070d0
        potcm_a3scc    = 0.531368d0
        potcm_b3scc    = 1.550070d0
        v              = 0.d0; g = 0.d0; d = 0.d0
        usricm_cart    = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                usricm_cart(iatom,idir) = x(iatom,idir)
            enddo
        enddo
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*legacy_hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            write (*, *) 'Only energy is available'
            !       do istate=1,nstates
            !         p(istate)=usrocm_pengygs*legacy_hartree_to_electronvolt
            !       enddo
            !       do iatom=1,natoms
            !       do idir=1,3
            !         g(1,iatom,idir)=usrocm_dgscart(iatom,idir)*51.422067
            !       enddo
            !       enddo
        elseif( igrad == 2 )then
            !       write (*,*) 'Only energy and gradient are available'
            write (*, *) 'Only energy is available'
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
            implicit none
            double precision :: coord(27), dx(27), en, vip, vop, vstr, vtor
            integer          :: i, n
            !      DIMENSION COORD(27),DX(27)
            !      DIMENSION dxvstr(27),dxvtor(27),dxvop(27),dxvip(27)
            !C      include 'c2h6cl.inc'
            !     write(97,*) 'pasa'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            ! Start of changes by EMN.
            ! The lines below do not seem to be needed because the Cartesian
            ! coordinates are in common qpdot.
            ! Besides, BOHR_A in utility.f does not match the value below.
            ! As a consequence, qpdot_q(i) changes and energy is not conserved
            !c call de venus
            call cartou
            call carttor
            !     PUT COORDINATES IN PROPER ARRAYS
            ! Inicialization de variables
            do i=1, 27
                qpdot_q(i) = pt1cm_r(i)
                !          qpdot_q(I) = pt1cm_r(I)*0.52918d0
                qpdot_pdot(i) = 0.d0
            enddo
            ! Changing units from angstrom to bohr
            !      DO  I = 1, 27
            !c         qpdot_q(I) = COORD(I)*0.52918d0
            !        qpdot_q(I) = pt1cm_r(I)*0.52918d0
            !      enddo
            ! End of corrections by EMN
            en = 0.0d0
            call coorden
            !  calculate switching functions
            call switchf
            !  calculate reference angles and their derivatives
            call refangles
            ! CRD 2011 incluimos los terminos de potencial de torsion
            call torsion(vtor)
            !     write(*,*) 'vtor', vtor*0.03812D0
            !  calculate stretching potential
            call stretch(vstr)
            !     write(*,*) 'vstr', vstr*0.03812D0
            !  calculate out of plane bending potential
            call opbend(vop)
            !     write(*,*) 'vop', vop*0.03812D0
            !  calculate in plane bending potential
            call ipbend(vip)
            !     write(*,*) 'vip', vip*0.03812D0
            en = vstr+vop+vip+vtor
            !  convert from 10(5) j/mol to au
            en = en*0.03812d0
            !      V = en
            pt1cm_engygs = en
            call eunitzero
            if( usricm_nder .ne. 0 )then
                !  Initialize derivatives
                do n=1, 27
                    !          DX(n)=0.0d0
                    pt1cm_degsdr(i) = qpdot_pdot(i)*0.0201723d0
                enddo
                ! Start of changes by EMN.
                ! The lines below are not needed because the derivatives are calculated
                ! in venuspotlib_3.f. Anyway usricm_nder=1
                !      PASO= 1.0D-8
                !     PASO= 1.0D-6
                !      do I=1,27
                !        qpdot_q(I)=qpdot_q(I) + PASO
                !        call coorden
                !        call switchf
                !        call refangles
                !        call torsion(vtor)
                !        call stretch(vstr)
                !        call opbend(vop)
                !        call ipbend(vip)
                !        en=vstr+vop+vip+vtor
                !C
                !c  convert from 10(5) j/mol to au
                !c
                !        en=en*0.03812D0
                !c
                !C       DX(I)=(en-V)/PASO
                !C       DX(I)=DX(I)*0.52918d0
                !C       qpdot_pdot(I)=DX(I)
                !        pt1cm_degsdr(I)=(en-pt1cm_engygs)/PASO
                !        pt1cm_degsdr(I)=pt1cm_degsdr(I)*0.52918d0
                !        qpdot_q(I)=qpdot_q(I)-PASO
                !      enddo
                !c
                ! 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
                ! qpdot_pdot sale de la derivada en hartrees/bohr por eso no lo multiplico
                !        do i=1,27
                !           pt1cm_degsdr(i)=qpdot_pdot(i)*0.0201723d0
                !           pt1cm_degsdr(i)=qpdot_pdot(i)
                !        enddo
                ! End of corrections by EMN
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate coorden for the C2H6O_VBMM potential.
        subroutine coorden
            !  calculates relative coordinates and bond lengths
            implicit none
            double precision :: argmax, argp3, argp4, p3, p4, t3tmp, t4tmp
            integer          :: i, ind
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
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
            !      do i=1,3
            !         do j=4,6
            !            do ind=1,3
            !          coords2_thh(i,j,ind)=qpdot_q(ndx_nh(j,ind))-qpdot_q(ndx_nh(i,ind))
            !            enddo
            !         enddo
            !      enddo
            !  calculate bond lengths
            bonds3_rc1b = sqrt(coords3_tc1b(1)*coords3_tc1b(1)+coords3_tc1b(2)*coords3_tc1b(2)+coords3_tc1b(3)* &
                & coords3_tc1b(3))
            bonds3_rc2b = sqrt(coords3_tc2b(1)*coords3_tc2b(1)+coords3_tc2b(2)*coords3_tc2b(2)+coords3_tc2b(3)* &
                & coords3_tc2b(3))
            potcm_rcc = sqrt(coords2_tcc(1)*coords2_tcc(1)+coords2_tcc(2)*coords2_tcc(2)+coords2_tcc(3)*coords2_tcc(3))
            !      write(*,*) bonds3_rc1b, bonds3_rc2b
            do i=1, 3
                bonds2_rc1h(i) = sqrt(coords1_tc1h(i,1)*coords1_tc1h(i,1)+coords1_tc1h(i,2)*coords1_tc1h(i,2)+ &
                    & coords1_tc1h(i,3)*coords1_tc1h(i,3))
                ! redondeo para corregir las frecuencias
                !       bonds2_rc1h(i)=(anint(bonds2_rc1h(i)*100000))/100000
                bonds2_rc2h(i) = sqrt(coords1_tc2h(i,1)*coords1_tc2h(i,1)+coords1_tc2h(i,2)*coords1_tc2h(i,2)+ &
                    & coords1_tc2h(i,3)*coords1_tc2h(i,3))
                ! redondeo para corregir las frecuencias
                !       bonds2_rc2h(i)=(anint(bonds2_rc2h(i)*100000))/100000
            enddo
            !      write(*,*) 'bonds2_rc1h(1)',bonds2_rc1h(1),'bonds2_rc1h(2)',bonds2_rc1h(2),'bonds2_rc1h(3)',bonds2_rc1h(3)
            !      write(*,*) 'bonds2_rc2h(1)',bonds2_rc2h(1),'bonds2_rc2h(2)',bonds2_rc2h(2),'bonds2_rc2h(3)',bonds2_rc2h(3)
            do i=1, 6
                bonds2_rbh(i) = sqrt(coords1_tbh(i,1)*coords1_tbh(i,1)+coords1_tbh(i,2)*coords1_tbh(i,2)+ &
                    & coords1_tbh(i,3)*coords1_tbh(i,3))
            enddo
            ! calculamos las distancias H-H
            !      do i=1,3
            !         do j=4,6
            !        bonds2_rhh(i,j)=sqrt(coords2_thh(i,j,1)*coords2_thh(i,j,1)+coords2_thh(i,j,2)*coords2_thh(i,j,2)+
            !    *                coords2_thh(i,j,3)*coords2_thh(i,j,3))
            !      write(55,*) 'bonds2_rhh', i, j, '=', bonds2_rhh(i,j)
            !         enddo
            !      enddo
            ! crd 2013 modificaciones de las distancias
            argmax = 19.d0
            p3     = 1.d0
            do i=1, 3
                argp3 = (potcm_w3*(bonds2_rc1h(i)-potcm_w4))
                if( argp3 .lt. argmax )then
                    t3tmp = 1-tanh(argp3)
                else
                    t3tmp = 0.d0
                endif
                p3 = p3*t3tmp
            enddo
            bonds_r0c1h = p3*potcm_r0c1hr+(1-p3)*potcm_r0c1hp
            !      bonds_r0c1h=potcm_r0c1hr
            p4 = 1.d0
            do i=1, 3
                argp4 = (potcm_w5*(bonds2_rc2h(i)-potcm_w6))
                if( argp4 .lt. argmax )then
                    t4tmp = 1-tanh(argp4)
                else
                    t4tmp = 0.d0
                endif
                p4 = p4*t4tmp
            enddo
            bonds_r0c2h = p4*potcm_r0c2hr+(1-p4)*potcm_r0c2hp
            !      bonds_r0c2h=potcm_r0c2hr
            !      write(*,*) 'bonds_r0c1h', bonds_r0c1h, 'bonds_r0c2h',bonds_r0c2h
        end subroutine coorden

        !> @brief Evaluate torsion for the C2H6O_VBMM potential.
        subroutine torsion( vtor )
            implicit none
            double precision :: dihed, t1(3), t2(3), vtor
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !      double precision pi
            !      pi=4.0d0*atan(1.0d0)
            ! CRD 2011 incluimos un apartado para calcular los angulos dihedro
            ! para poder calcular despues vtorsion
            call dihedro(3, 1, 2, 6, dihed)
            thetahcch_thetahcch(1) = dihed*pi/180.0d0
            call dihedro(3, 1, 2, 7, dihed)
            thetahcch_thetahcch(2) = dihed*pi/180.0d0
            call dihedro(3, 1, 2, 8, dihed)
            thetahcch_thetahcch(3) = dihed*pi/180.0d0
            call dihedro(4, 1, 2, 6, dihed)
            thetahcch_thetahcch(4) = dihed*pi/180.0d0
            call dihedro(4, 1, 2, 7, dihed)
            thetahcch_thetahcch(5) = dihed*pi/180.0d0
            call dihedro(4, 1, 2, 8, dihed)
            thetahcch_thetahcch(6) = dihed*pi/180.0d0
            call dihedro(5, 1, 2, 6, dihed)
            thetahcch_thetahcch(7) = dihed*pi/180.0d0
            call dihedro(5, 1, 2, 7, dihed)
            thetahcch_thetahcch(8) = dihed*pi/180.0d0
            call dihedro(5, 1, 2, 8, dihed)
            thetahcch_thetahcch(9) = dihed*pi/180.0d0
            ! CRD 2011 calculamos los terminos switching t1 y t2
            do i=1, 3
                t1(i) = 0.5d0*(1-tanh(potcm_w1*(bonds2_rc1h(i)-potcm_w2)))
                t2(i) = 0.5d0*(1-tanh(potcm_w1*(bonds2_rc2h(i)-potcm_w2)))
                !       write(*,*) 't1',i,t1(i),'t2',i,t2(i)
                !       t1(i)=1-tanh(potcm_w1*(bonds2_rc1h(i)-potcm_w2))
                !       t2(i)=1-tanh(potcm_w1*(bonds2_rc2h(i)-potcm_w2))
            enddo
            ! CRD 2011 inicializamos vtor
            vtor = 0.0d0
            vtor = vtor+(potcm_v3/3.0d0)*(1+cos(3*thetahcch_thetahcch(1)))*t1(1)*t2(1)
            vtor = vtor+(potcm_v3/3.0d0)*(1+cos(3*thetahcch_thetahcch(2)))*t1(1)*t2(2)
            vtor = vtor+(potcm_v3/3.0d0)*(1+cos(3*thetahcch_thetahcch(3)))*t1(1)*t2(3)
            vtor = vtor+(potcm_v3/3.0d0)*(1+cos(3*thetahcch_thetahcch(4)))*t1(2)*t2(1)
            vtor = vtor+(potcm_v3/3.0d0)*(1+cos(3*thetahcch_thetahcch(5)))*t1(2)*t2(2)
            vtor = vtor+(potcm_v3/3.0d0)*(1+cos(3*thetahcch_thetahcch(6)))*t1(2)*t2(3)
            vtor = vtor+(potcm_v3/3.0d0)*(1+cos(3*thetahcch_thetahcch(7)))*t1(3)*t2(1)
            vtor = vtor+(potcm_v3/3.0d0)*(1+cos(3*thetahcch_thetahcch(8)))*t1(3)*t2(2)
            vtor = vtor+(potcm_v3/3.0d0)*(1+cos(3*thetahcch_thetahcch(9)))*t1(3)*t2(3)
        end subroutine torsion

        !> @brief Evaluate dihedro for the C2H6O_VBMM potential.
        subroutine dihedro( i, j, k, l, dihed )
            implicit none
            double precision :: dihed
            !        include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            integer :: i, j, k, l
            integer :: m
            !        double precision x(300)
            !        double precision dista
            !        double precision d12,d13,d14,d23,d24,d34
            !        double precision p,qpdot_q
            !        double precision p
            !        double precision pi
            double precision :: rij(3), rjk(3), rkj(3), rkl(3)
            double precision :: zijjk(3), zijkj(3), zkjkl(3), zjkkl(3)
            double precision :: tmpsign(3), tmp, arg1, arg2
            !        pi=4.0d0*atan(1.0d0)
            do m = 1, 3
                rij(m) = qpdot_q((i-1)*3+m)-qpdot_q((j-1)*3+m)
                rjk(m) = qpdot_q((j-1)*3+m)-qpdot_q((k-1)*3+m)
                rkl(m) = qpdot_q((k-1)*3+m)-qpdot_q((l-1)*3+m)
            enddo
            call crosprod(rjk, rkl, zjkkl)
            call crosprod(rij, rjk, zijjk)
            arg2 = anorm(rjk)*dotprod(rij, zjkkl)
            arg1 = dotprod(zijjk, zjkkl)
            if( abs(arg1) .lt. 1.d-6 )then
                if( arg2 .gt. 0.d0 ) dihed = 90.d0
                if( arg2 .lt. 0.d0 ) dihed = -90.d0
                return
            endif
            if( arg1 .gt. 0.d0 )then
                dihed = atan(arg2/arg1) *180.d0/pi
                return
            endif
            if( arg2 .ge. 0.d0 )then
                dihed = 180.d0 + atan(arg2/arg1) *180.d0/pi
                return
            endif
            if( arg2 .lt. 0.d0 )then
                dihed = -180.d0 + atan(arg2/arg1) *180.d0/pi
                return
            endif
        end subroutine dihedro

        !> @brief Evaluate crosprod for the C2H6O_VBMM potential.
        subroutine crosprod( x, y, z )
            implicit none
            double precision :: x(3), y(3), z(3)
            z(1) = x(2)*y(3)-x(3)*y(2)
            z(2) = x(3)*y(1)-x(1)*y(3)
            z(3) = x(1)*y(2)-x(2)*y(1)
        end subroutine crosprod

        !> @brief Evaluate dotprod for the C2H6O_VBMM potential.
        double precision function dotprod( x, y )
            implicit none
            integer          :: k
            double precision :: x(3), y(3), tmp
            tmp = 0.d0
            do k = 1, 3
                tmp = tmp + x(k)*y(k)
            enddo
            dotprod = tmp
        end function dotprod

        !> @brief Evaluate anorm for the C2H6O_VBMM potential.
        double precision function anorm( x )
            implicit none
            integer          :: k
            double precision :: x(3), tmp
            tmp = 0.d0
            do k = 1, 3
                tmp = tmp + x(k)*x(k)
            enddo
            anorm = dsqrt(tmp)
        end function anorm

        !> @brief Evaluate refangles for the C2H6O_VBMM potential.
        subroutine refangles
            !  subroutine calculates reference angles for the "in-plane" potential
            implicit none
            double precision :: ddr(4,4), halfpi, sumd2(4), sumd4(4), tau, taucc, twopi
            integer          :: i, j
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            tau = acos(-1.0d0/3.0d0)
            ! CRD 2012 taucc para el angulo H-C-C
            taucc = acos(-0.3616d0)
            !      write(*,*) 'tau', tau, 'acostau', tau*57.2957795
            !      write(*,*) 'taucc', taucc, 'acostaucc', taucc*57.2957795
            !      pi=4.0d0*atan(1.0d0)
            halfpi = 0.5d0*pi
            twopi  = 2.0d0*pi
            !      write(38,*) tau, pi, halfpi, twopi
            ! jcc-2010
            ! valor para ajustar el angulo
            !      anghch=0.649d0*pi
            !      tausih = pi - asin ( sin(anghch/2.d0) / sin(pi/3.d0) )
            !      tausih=acos(-0.4648d0)
            ! jcc-2010
            !  set diagonal elements to zero
            do i=1, 8
                angles_theta0(i,i) = 0.0d0
            enddo
            !  calculate reference angles
            !      angles_theta0(1,2)=tau+(tau-tausih)*(switch1_sphi(1)*switch1_sphi(2)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(3)*switch1_stheta(4)-1.0d0)
            !      angles_theta0(1,3)=tau+(tau-tausih)*(switch1_sphi(1)*switch1_sphi(3)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(2)*switch1_stheta(4)-1.0d0)
            !      angles_theta0(1,4)=taucc+(tau-tausih)*(switch1_sphi(1)*switch1_sphi(4)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(2)*switch1_stheta(3)-1.0d0)
            !      angles_theta0(2,3)=tau+(tau-tausih)*(switch1_sphi(2)*switch1_sphi(3)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(1)*switch1_stheta(4)-1.0d0)
            !      angles_theta0(2,4)=taucc+(tau-tausih)*(switch1_sphi(2)*switch1_sphi(4)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(1)*switch1_stheta(3)-1.0d0)
            !      angles_theta0(3,4)=taucc+(tau-tausih)*(switch1_sphi(3)*switch1_sphi(4)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(1)*switch1_stheta(2)-1.0d0)
            !      angles_theta0(5,6)=tau+(tau-tausih)*(switch1_sphi(5)*switch1_sphi(6)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(7)*switch1_stheta(8)-1.0d0)
            !      angles_theta0(5,7)=tau+(tau-tausih)*(switch1_sphi(5)*switch1_sphi(7)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(6)*switch1_stheta(8)-1.0d0)
            !      angles_theta0(5,8)=taucc+(tau-tausih)*(switch1_sphi(5)*switch1_sphi(8)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(6)*switch1_stheta(7)-1.0d0)
            !      angles_theta0(6,7)=tau+(tau-tausih)*(switch1_sphi(6)*switch1_sphi(7)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(5)*switch1_stheta(8)-1.0d0)
            !      angles_theta0(6,8)=taucc+(tau-tausih)*(switch1_sphi(6)*switch1_sphi(8)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(5)*switch1_stheta(7)-1.0d0)
            !      angles_theta0(7,8)=taucc+(tau-tausih)*(switch1_sphi(7)*switch1_sphi(8)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(switch1_stheta(5)*switch1_stheta(6)-1.0d0)
            ! el atomo 4 es el c2
            angles_theta0(1,2) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(2)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(3)*switch1_stheta(4)-1.0d0)
            angles_theta0(1,3) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(3)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(2)*switch1_stheta(4)-1.0d0)
            angles_theta0(1,4) = taucc+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(4)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(2)*switch1_stheta(3)-1.0d0)
            angles_theta0(2,3) = tau+(tau-halfpi)*(switch1_sphi(2)*switch1_sphi(3)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(1)*switch1_stheta(4)-1.0d0)
            angles_theta0(2,4) = taucc+(tau-halfpi)*(switch1_sphi(2)*switch1_sphi(4)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(1)*switch1_stheta(3)-1.0d0)
            angles_theta0(3,4) = taucc+(tau-halfpi)*(switch1_sphi(3)*switch1_sphi(4)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(1)*switch1_stheta(2)-1.0d0)
            !  el atomo 8 es el c1
            angles_theta0(5,6) = tau+(tau-halfpi)*(switch1_sphi(5)*switch1_sphi(6)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(7)*switch1_stheta(8)-1.0d0)
            angles_theta0(5,7) = tau+(tau-halfpi)*(switch1_sphi(5)*switch1_sphi(7)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(6)*switch1_stheta(8)-1.0d0)
            angles_theta0(5,8) = taucc+(tau-halfpi)*(switch1_sphi(5)*switch1_sphi(8)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(6)*switch1_stheta(7)-1.0d0)
            angles_theta0(6,7) = tau+(tau-halfpi)*(switch1_sphi(6)*switch1_sphi(7)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(5)*switch1_stheta(8)-1.0d0)
            angles_theta0(6,8) = taucc+(tau-halfpi)*(switch1_sphi(6)*switch1_sphi(8)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(5)*switch1_stheta(7)-1.0d0)
            angles_theta0(7,8) = taucc+(tau-halfpi)*(switch1_sphi(7)*switch1_sphi(8)-1.0d0) +(tau-twopi/3.0d0)* &
                & (switch1_stheta(5)*switch1_stheta(6)-1.0d0)
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
            !      do i=1,3
            !         do j=i+1,4
            !          write(*,*) 'theta',i,j,'=', angles_theta0(i,j)*180/pi
            !         enddo
            !      enddo
            !      do i=5,7
            !         do j=i+1,8
            !          write(*,*) 'theta',i,j,'=', angles_theta0(i,j)*180/pi
            !         enddo
            !      enddo
        end subroutine refangles

        !> @brief Evaluate stretch for the C2H6O_VBMM potential.
        subroutine stretch( vstr )
            !  subroutine to calculate leps-type stretching potential and its
            !  derivatives
            implicit none
            double precision :: ach1, ach2, achdc(3), achdh(4,3), achdx(3), arga1, arga2, dt, e1, e3, expterm, factj(4), rav1, rav2, &
                & vcc, vj(3), vjbh(6), vjc1b, vjc1h(3), vjc2b, vjc2h(3), vq(3), vqbh(6), vqc1b, vqc1h(3), vqc2b, vqc2h(3), vstr
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !      write(*,*) potcm_nnc1,potcm_nnc2,potcm_nnb
            !      write(*,*) potcm_r0c1hr,potcm_r0c1hp,potcm_w3,potcm_w4,potcm_d1c1h,potcm_d3c1h
            !      write(*,*) potcm_a1c1h,potcm_b1c1h,potcm_c1c1h
            !      write(*,*) potcm_r0c2hr,potcm_r0c2hp,potcm_w5,potcm_w6,potcm_d1c2h,potcm_d3c2h
            !      write(*,*) potcm_a1c2h,potcm_b1c2h,potcm_c1c2h
            !      write(*,*) potcm_r0cc,potcm_rcc,potcm_d1cc,potcm_acc
            !      write(*,*) potcm_a1c,potcm_a2c
            !      write(*,*) potcm_r0hh,potcm_d1hh,potcm_d3hh,potcm_ahh
            !      write(*,*) potcm_r0cb,potcm_d1cb,potcm_d3cb,potcm_acb
            !      write(*,*) potcm_a3s,potcm_b3s,potcm_aphi,potcm_bphi,potcm_cphi
            !      write(*,*) potcm_atheta,potcm_btheta,potcm_ctheta
            !      write(*,*) potcm_fch3,potcm_hch3
            !      write(*,*) potcm_fkinf,potcm_ak,potcm_bk,potcm_aa1,potcm_aa2,potcm_aa3,potcm_aa4
            !      write(*,*) potcm_v3,potcm_w1,potcm_w2
            !      write(*,*) potcm_a1s,potcm_b1s,potcm_a2s,potcm_b2s
            !      write(*,*) potcm_a1scc,potcm_b1scc,potcm_a2scc,potcm_b2scc,potcm_a3scc,potcm_b3scc
            !  calculate avergage bond length for the methane moiety
            !C       rav=(bonds2_rc1h(1)+bonds2_rc1h(2)+bonds2_rc1h(3)+bonds2_rc2h(1)+bonds2_rc2h(2)+bonds2_rc2h(3))/6.0d0
            rav1 = (bonds2_rc1h(1)+bonds2_rc1h(2)+bonds2_rc1h(3))/3.0d0
            rav2 = (bonds2_rc2h(1)+bonds2_rc2h(2)+bonds2_rc2h(3))/3.0d0
            !      write(*,*) 'rav1',rav1,'rav2', rav2
            !  initialise:
            vstr = 0.0d0
            !  ach:
            !  nb: in double precision tanh(19.0d0)=1.0d0 and we put the if statemen
            !  in to avoid overflow/underflow errors
            !     write(*,*) 'bonds_r0c1h',bonds_r0c1h,'bonds_r0c2h',bonds_r0c2h
            arga1 = potcm_c1c1h*(rav1-bonds_r0c1h)
            arga2 = potcm_c1c2h*(rav2-bonds_r0c2h)
            !      write(*,*) 'arga1', arga1,'arga2',arga2
            if( arga1 .lt. 19.0d0 )then
                ach1 = potcm_a1c1h+potcm_b1c1h*(tanh(arga1)+1.0d0)*0.5d0
            else
                ach1 = potcm_a1c1h+potcm_b1c1h
            endif
            if( arga2 .lt. 19.0d0 )then
                ach2 = potcm_a1c2h+potcm_b1c2h*(tanh(arga2)+1.0d0)*0.5d0
            else
                ach2 = potcm_a1c2h+potcm_b1c2h
            endif
            !      write(6,*) bonds_r0c1h, rav1, bonds_r0c2h, rav2, ach1, ach2
            !      argac=stret1_c1cc*(potcm_rcc-potcm_r0cc)
            !      if(argac.lt.19.0d0)then
            !        potcm_acc=stret1_a1cc+stret1_b1cc*(tanh(argac)+1.0d0)*0.5d0
            !      else
            !        potcm_acc=stret1_a1cc+stret1_b1cc
            !      endif
            !  calculate singlet: e1, triplet: e3 energies and vq and vj
            !  terms for each bond
            e1    = potcm_d1cb*(exp(-2.0d0*potcm_acb*(bonds3_rc1b-potcm_r0cb))-2.0d0*exp(-potcm_acb*(bonds3_rc1b- &
                & potcm_r0cb)))
            e3    = potcm_d3cb*(exp(-2.0d0*potcm_acb*(bonds3_rc1b-potcm_r0cb))+2.0d0*exp(-potcm_acb*(bonds3_rc1b- &
                & potcm_r0cb)))
            vqc1b = (e1+e3)*0.5d0
            vjc1b = (e1-e3)*0.5d0
            e1    = potcm_d1cb*(exp(-2.0d0*potcm_acb*(bonds3_rc2b-potcm_r0cb))-2.0d0*exp(-potcm_acb*(bonds3_rc2b- &
                & potcm_r0cb)))
            e3    = potcm_d3cb*(exp(-2.0d0*potcm_acb*(bonds3_rc2b-potcm_r0cb))+2.0d0*exp(-potcm_acb*(bonds3_rc2b- &
                & potcm_r0cb)))
            vqc2b = (e1+e3)*0.5d0
            vjc2b = (e1-e3)*0.5d0
            !     write(*,*) 'vqc1b',vqc1b,'vqc2b',vqc2b,'vjc1b',vjc1b,'vjc2b',vjc2b
            !      e1=potcm_d1cc*(exp(-2.0d0*potcm_acc*(potcm_rcc-potcm_r0cc))-2.0d0*exp(-potcm_acc*(potcm_rcc-potcm_r0cc)))
            !      e3=d3cc*(exp(-2.0d0*potcm_acc*(potcm_rcc-potcm_r0cc))+2.0d0*exp(-potcm_acc*(potcm_rcc-potcm_r0cc)))
            !      vqcc=(e1+e3)*0.5d0
            !      vjcc=(e1-e3)*0.5d0
            do i=1, 3
                e1       = potcm_d1c1h*(exp(-2.0d0*ach1*(bonds2_rc1h(i)-bonds_r0c1h)) -2.0d0*exp(-ach1* &
                    & (bonds2_rc1h(i)-bonds_r0c1h)))
                e3       = potcm_d3c1h*(exp(-2.0d0*ach1*(bonds2_rc1h(i)-bonds_r0c1h)) +2.0d0*exp(-ach1* &
                    & (bonds2_rc1h(i)-bonds_r0c1h)))
                vqc1h(i) = (e1+e3)*0.5d0
                vjc1h(i) = (e1-e3)*0.5d0
                e1       = potcm_d1c2h*(exp(-2.0d0*ach2*(bonds2_rc2h(i)-bonds_r0c2h)) -2.0d0*exp(-ach2* &
                    & (bonds2_rc2h(i)-bonds_r0c2h)))
                e3       = potcm_d3c2h*(exp(-2.0d0*ach2*(bonds2_rc2h(i)-bonds_r0c2h)) +2.0d0*exp(-ach2* &
                    & (bonds2_rc2h(i)-bonds_r0c2h)))
                vqc2h(i) = (e1+e3)*0.5d0
                vjc2h(i) = (e1-e3)*0.5d0
            enddo
            do i=1, 6
                e1      = potcm_d1hh*(exp(-2.0d0*potcm_ahh*(bonds2_rbh(i)-potcm_r0hh)) -2.0d0*exp(-potcm_ahh* &
                    & (bonds2_rbh(i)-potcm_r0hh)))
                e3      = potcm_d3hh*(exp(-2.0d0*potcm_ahh*(bonds2_rbh(i)-potcm_r0hh)) +2.0d0*exp(-potcm_ahh* &
                    & (bonds2_rbh(i)-potcm_r0hh)))
                vqbh(i) = (e1+e3)*0.5d0
                vjbh(i) = (e1-e3)*0.5d0
            enddo
            !  calculate 3 body potential
            do i=1, 3
                vq(i) = vqc1h(i)+vqc1b+vqbh(i)
                vj(i) = -sqrt(((vjc1h(i)-vjc1b)**2+(vjc1b-vjbh(i))**2 +(vjbh(i)-vjc1h(i))**2)*0.5d0)
                vstr  = vstr+vq(i)+vj(i)
            enddo
            do i=1, 3
                vq(i) = vqc2h(i)+vqc2b+vqbh(i+3)
                vj(i) = -sqrt(((vjc2h(i)-vjc2b)**2+(vjc2b-vjbh(i+3))**2 +(vjbh(i+3)-vjc2h(i))**2)*0.5d0)
                vstr  = vstr+vq(i)+vj(i)
            enddo
            ! CRD 2012 potencial morse para el enlace C-C
            dt      = (potcm_rcc-potcm_r0cc)
            expterm = exp(-potcm_acc*dt)
            vcc     = potcm_d1cc*(1.0d0-expterm)**2.0d0
            vstr    = vstr + vcc
        end subroutine stretch

        !> @brief Evaluate opbend for the C2H6O_VBMM potential.
        subroutine opbend( vop )
            !  subroutine calculates symmetrized vop potential and derivatives
            implicit none
            double precision :: a(3), argd(8), axb(3), b(3), c(8,3), sum2, sum4, sumd2(8), sumd4(8), vop
            integer          :: i, in(3), j, k, l
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            double precision :: norma
            vop = 0.0d0
            !  calculate force constants and their derivatives
            call opforce
            !  calculate out-of-plane angle and derivatives
            ! CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C2
            !      do i=1,4
            i = 4
            !      do i=1,3
            j = i+1
            if( j .gt. 4 ) j=j-4
            k = j+1
            if( k .gt. 4 ) k=k-4
            l = k+1
            if( l .gt. 4 ) l=l-4
            call calcdelta1(i, j, k, l, sum2, sum4)
            sumd2(i) = sum2
            sumd4(i) = sum4
            !      write(*,*) 'sum2',sum2,'sum4',sum4
            !      write(*,*) 'sumd2',sumd2(i),'sumd4',sumd4(i)
            vop = vop+delta1_fdelta(i)*sumd2(i)+delta1_hdelta(i)*sumd4(i)
            !      write(*,*) i,j,k,l,sumd2(i),delta1_fdelta(i),sumd4(i),delta1_hdelta(i)
            !      enddo
            ! CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C1
            !      do i=5,8
            i = 8
            !      do i=4,6
            j = i+1
            if( j .gt. 8 ) j=j-4
            k = j+1
            if( k .gt. 8 ) k=k-4
            l = k+1
            if( l .gt. 8 ) l=l-4
            call calcdelta2(i, j, k, l, sum2, sum4)
            sumd2(i) = sum2
            sumd4(i) = sum4
            !      write(*,*) 'sum2',sum2,'sum4',sum4
            !      write(*,*) 'sumd2',sumd2(i),'sumd4',sumd4(i)
            vop = vop+delta1_fdelta(i)*sumd2(i)+delta1_hdelta(i)*sumd4(i)
            !      write(*,*) i,j,k,l,sumd2(i),delta1_fdelta(i),sumd4(i),delta1_hdelta(i)
            !      enddo
        end subroutine opbend

        !> @brief Evaluate ipbend for the C2H6O_VBMM potential.
        subroutine ipbend( vip )
            !  subroutine calculates symmetrised in plane bend term
            !  and its derivatives
            implicit none
            double precision :: costh(8,8), dth(8,8), theta(8,8), vip
            integer          :: i, j
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  initialise
            vip = 0.0d0
            !  calculate force constants: force1_fk0(i,j), force1_f1(i)
            !  and derivatives wrt rch(k) and bonds2_rbh(k): force1_dfdc(i,j,k), force1_dfdh(i,j,k)
            call ipforce
            !      do i=1,8
            !      write(*,*) force1_f1(i)
            !      enddo
            !      write(*,*) 'fin force1_f1(i)'
            !  calculate theta(i,j) and in plane bend potential
            ! bucle para el c1
            do i=1, 3
                do j=i+1, 4
                    if( j .lt. 4.0 )then
                        costh(i,j) = coords1_tc1h(i,1)*coords1_tc1h(j,1)+coords1_tc1h(i,2)*coords1_tc1h(j,2) + &
                            & coords1_tc1h(i,3)*coords1_tc1h(j,3)
                        costh(i,j) = costh(i,j)/bonds2_rc1h(i)/bonds2_rc1h(j)
                        theta(i,j) = acos(costh(i,j))
                        dth(i,j)   = theta(i,j)-angles_theta0(i,j)
                        vip        = vip+0.5d0*force1_fk0(i,j)*force1_f1(i)*force1_f1(j)*dth(i,j)**2
                    else
                        costh(i,j) = coords1_tc1h(i,1)*coords2_tcc(1)+coords1_tc1h(i,2)*coords2_tcc(2) + &
                            & coords1_tc1h(i,3)*coords2_tcc(3)
                        !      write(92,*) 'bonds2_rc1h',i,'=', bonds2_rc1h(i), 'potcm_rcc=', potcm_rcc
                        costh(i,j) = costh(i,j)/bonds2_rc1h(i)/potcm_rcc
                        theta(i,j) = acos(costh(i,j))
                        dth(i,j)   = theta(i,j)-angles_theta0(i,j)
                        vip        = vip+0.5d0*force1_fk0(i,j)*force1_f1(i)*force1_f1(j)*dth(i,j)**2
                    endif
                    !      write(*,*) 'force1_fk0',i,j,'=',force1_fk0(i,j)
                    !      write(*,*) 'force1_f1',i,'=',force1_f1(i)
                    !      write(*,*) 'force1_f1',j,'=',force1_f1(j)
                    !      write(*,*) 'potcm_rcc','=',potcm_rcc
                    !      write(*,*) 'costh',i,j,'=',costh(i,j)
                    !      write(*,*) 'theta',i,j,'=',theta(i,j)
                    !      write(*,*) 'angles_theta0',i,j,'=',angles_theta0(i,j)
                    !      write(*,*) 'dth',i,j,'=',dth(i,j)
                    !      write(*,*) 'costh',i,j,'=',costh(i,j)
                enddo
            enddo
            ! bucle para el c2 para los theta 5,6,7,8
            do i=1, 3
                do j=i+1, 4
                    if( j .lt. 4.0 )then
                        costh(i+4,j+4) = coords1_tc2h(i,1)*coords1_tc2h(j,1)+coords1_tc2h(i,2)*coords1_tc2h(j,2) + &
                            & coords1_tc2h(i,3)*coords1_tc2h(j,3)
                        costh(i+4,j+4) = costh(i+4,j+4)/bonds2_rc2h(i)/bonds2_rc2h(j)
                        theta(i+4,j+4) = acos(costh(i+4,j+4))
                        dth(i+4,j+4)   = theta(i+4,j+4)-angles_theta0(i+4,j+4)
                        vip            = vip+0.5d0*force1_fk0(i+4,j+4)*force1_f1(i+4)*force1_f1(j+4)*dth(i+4,j+4)**2
                    else
                        ! modificacion para corregir los angulos del c2
                        !          costh(i,j)=coords1_tc2h(i,1)*coords2_tcc(1)+coords1_tc2h(i,2)*coords2_tcc(2)
                        !    *                       +coords1_tc2h(i,3)*coords2_tcc(3)
                        costh(i+4,j+4) = coords1_tc2h(i,1)*(-coords2_tcc(1))+coords1_tc2h(i,2)*(-coords2_tcc(2)) + &
                            & coords1_tc2h(i,3)*(-coords2_tcc(3))
                        costh(i+4,j+4) = costh(i+4,j+4)/bonds2_rc2h(i)/potcm_rcc
                        theta(i+4,j+4) = acos(costh(i+4,j+4))
                        dth(i+4,j+4)   = theta(i+4,j+4)-angles_theta0(i+4,j+4)
                        vip            = vip+0.5d0*force1_fk0(i+4,j+4)*force1_f1(i+4)*force1_f1(j+4)*dth(i+4,j+4)**2
                    endif
                enddo
            enddo
        end subroutine ipbend

        !> @brief Evaluate calcdelta1 for the C2H6O_VBMM potential.
        subroutine calcdelta1( i, j, k, l, sum2, sum4 )
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            implicit none
            double precision :: a(3), argd(8), atemp2(3), axb(3), b(3), c(8,3), cdot(8,3,3), daxb(8,3,3), delta(8), sum2, sum4
            integer          :: i, ii, in(3), ind, j, k, l
            double precision :: norma
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  initialise
            sum2 = 0.0d0; sum4 = 0.0d0
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
                !      axb(1)=b(2)*a(3)-b(3)*a(2)
                !      axb(2)=a(1)*b(3)-b(1)*a(3)
                !      axb(3)=b(1)*a(2)-a(1)*b(2)
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
                !      write(*,*) 'norma=',norma
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
                        !      write(*,*) 'calcdelta1 in(ii)', in(ii),c(in(ii),ind)
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
                if( argd(in(ii)) > 1.d0 ) argd(in(ii))  = 1
                if( argd(in(ii)) < -1.d0 ) argd(in(ii)) = -1
                ! End of corrections by EMN
                delta(in(ii)) = acos(argd(in(ii)))-angles_theta0(i,in(ii))
                sum2          = sum2+delta(in(ii))**2
                sum4          = sum4+delta(in(ii))**4
                !      write(*,*) 'calcdelta1 sum2', sum2
                !      write(*,*) 'calcdelta1 sum4', sum4
            enddo
        end subroutine calcdelta1

        !> @brief Evaluate calcdelta2 for the C2H6O_VBMM potential.
        subroutine calcdelta2( i, j, k, l, sum2, sum4 )
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            implicit none
            double precision :: a(3), argd(8), atemp2(3), axb(3), b(3), c(8,3), cdot(8,3,3), daxb(8,3,3), delta(8), sum2, sum4
            integer          :: i, ii, in(3), ind, j, k, l
            double precision :: norma
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  initialise
            sum2 = 0.0d0; sum4 = 0.0d0
            !  set j,k,l indices
            !      in(1)=j
            !      in(2)=k
            !      in(3)=l
            in(1) = 4
            in(2) = 5
            in(3) = 6
            ! CRD 2005 incluimos bucles if para asegurar que siempre
            ! que una de los valores j, k o l sea 4 las coordenadas
            ! corresponden al atomo nuevo.
            !  vector a is rk-rj, vector b is rl-rj
            if( j .eq. 8 )then
                do ind=1, 3
                    a(ind) = qpdot_q(ndx_nh(k,ind))-qpdot_q(ndx2_nc1(ind))
                    b(ind) = qpdot_q(ndx_nh(l,ind))-qpdot_q(ndx2_nc1(ind))
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
                    a(ind) = qpdot_q(ndx2_nc1(ind))-qpdot_q(ndx_nh(j,ind))
                    b(ind) = qpdot_q(ndx_nh(l,ind))-qpdot_q(ndx_nh(j,ind))
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
                    a(ind) = qpdot_q(ndx_nh(k,ind))-qpdot_q(ndx_nh(j,ind))
                    b(ind) = qpdot_q(ndx2_nc1(ind))-qpdot_q(ndx_nh(j,ind))
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
                    a(ind) = qpdot_q(ndx_nh(k,ind))-qpdot_q(ndx_nh(j,ind))
                    b(ind) = qpdot_q(ndx_nh(l,ind))-qpdot_q(ndx_nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
                !      write(*,*) 'calcdelta2 norma',norma
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
                        c(in(ii),ind) = -coords1_tc2h(in(ii),ind)/bonds2_rc2h(in(ii))
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
                if( argd(in(ii)) > 1.d0 ) argd(in(ii))  = 1
                if( argd(in(ii)) < -1.d0 ) argd(in(ii)) = -1
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
            implicit none
            double precision :: switch(8)
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  calculate switching functions:
            switch(1) = (1.0d0-op1_s3(1))*op1_s3(2)*op1_s3(3)*op1_s3(4)
            switch(2) = (1.0d0-op1_s3(2))*op1_s3(3)*op1_s3(4)*op1_s3(1)
            switch(3) = (1.0d0-op1_s3(3))*op1_s3(4)*op1_s3(1)*op1_s3(2)
            switch(4) = (1.0d0-op1_s3(4))*op1_s3(1)*op1_s3(2)*op1_s3(3)
            switch(5) = (1.0d0-op1_s3(5))*op1_s3(6)*op1_s3(7)*op1_s3(8)
            switch(6) = (1.0d0-op1_s3(6))*op1_s3(7)*op1_s3(8)*op1_s3(5)
            switch(7) = (1.0d0-op1_s3(7))*op1_s3(8)*op1_s3(5)*op1_s3(6)
            switch(8) = (1.0d0-op1_s3(8))*op1_s3(5)*op1_s3(6)*op1_s3(7)
            !      do i=1,8
            !       write(*,*) 'switch',i,'=',switch(i)
            !      enddo
            !       write(*,*) 'pasa'
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
            implicit none
            double precision :: a1, a2, arga1, arga2, f0, f2
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  set force constant at asymptotes
            f0              = potcm_fkinf+potcm_ak
            f2              = potcm_fkinf
            force1_fk0(1,2) = f0+f0*(ip1_s1(1)*ip1_s1(2)-1.0d0)+(f0-f2)*(ip1_s2(3)*ip1_s2(4)-1.0d0)
            force1_fk0(1,3) = f0+f0*(ip1_s1(1)*ip1_s1(3)-1.0d0)+(f0-f2)*(ip1_s2(2)*ip1_s2(4)-1.0d0)
            force1_fk0(1,4) = f0+f0*(ip1_s1(1)*ip1_s1(4)-1.0d0)+(f0-f2)*(ip1_s2(2)*ip1_s2(3)-1.0d0)
            force1_fk0(2,3) = f0+f0*(ip1_s1(2)*ip1_s1(3)-1.0d0)+(f0-f2)*(ip1_s2(1)*ip1_s2(4)-1.0d0)
            force1_fk0(2,4) = f0+f0*(ip1_s1(2)*ip1_s1(4)-1.0d0)+(f0-f2)*(ip1_s2(1)*ip1_s2(3)-1.0d0)
            force1_fk0(3,4) = f0+f0*(ip1_s1(3)*ip1_s1(4)-1.0d0)+(f0-f2)*(ip1_s2(1)*ip1_s2(2)-1.0d0)
            force1_fk0(5,6) = f0+f0*(ip1_s1(5)*ip1_s1(6)-1.0d0)+(f0-f2)*(ip1_s2(7)*ip1_s2(8)-1.0d0)
            force1_fk0(5,7) = f0+f0*(ip1_s1(5)*ip1_s1(7)-1.0d0)+(f0-f2)*(ip1_s2(6)*ip1_s2(8)-1.0d0)
            force1_fk0(5,8) = f0+f0*(ip1_s1(5)*ip1_s1(8)-1.0d0)+(f0-f2)*(ip1_s2(6)*ip1_s2(7)-1.0d0)
            force1_fk0(6,7) = f0+f0*(ip1_s1(6)*ip1_s1(7)-1.0d0)+(f0-f2)*(ip1_s2(5)*ip1_s2(8)-1.0d0)
            force1_fk0(6,8) = f0+f0*(ip1_s1(6)*ip1_s1(8)-1.0d0)+(f0-f2)*(ip1_s2(5)*ip1_s2(7)-1.0d0)
            force1_fk0(7,8) = f0+f0*(ip1_s1(7)*ip1_s1(8)-1.0d0)+(f0-f2)*(ip1_s2(5)*ip1_s2(6)-1.0d0)
            ! force1_f1 para los h del c1
            do i=1, 3
                !  calc derivatives of force1_fk0 wrt each of the rch(i) bonds
                !  calculate the terms force1_f1(i)
                arga1        = potcm_aa1*bonds2_rbh(i)*bonds2_rbh(i)
                arga2        = potcm_aa4*(bonds2_rbh(i)-potcm_r0hh)*(bonds2_rbh(i)-potcm_r0hh)
                a1           = 1.0d0-exp(-arga1)
                a2           = potcm_aa2+potcm_aa3*exp(-arga2)
                force1_f1(i) = a1*exp(-a2*(bonds2_rc1h(i)-bonds_r0c1h)**2)
            enddo
            ! force1_f1 para el c2
            ! CRD 2012 eliminamos el calculo de potcm_a1c y potcm_a2c pq es constante
            ! lo incluimos en el CONST
            !        arga1c=potcm_aa1*rbc*rbc
            !        arga2c=potcm_aa4*(rbc-potcm_r0cb)*(rbc-potcm_r0cb)
            !        potcm_a1c=1.0d0-exp(-arga1c)
            !        potcm_a2c=potcm_aa2+potcm_aa3*exp(-arga2c)
            force1_f1(4) = potcm_a1c*exp(-potcm_a2c*(potcm_rcc-potcm_r0cc)**2)
            force1_f1(8) = force1_f1(4)
            ! force1_f1 para los h del c2
            do i=5, 7
                !  calc derivatives of force1_fk0 wrt each of the rch(i) bonds
                !  calculate the terms force1_f1(i)
                arga1        = potcm_aa1*bonds2_rbh(i-1)*bonds2_rbh(i-1)
                arga2        = potcm_aa4*(bonds2_rbh(i-1)-potcm_r0hh)*(bonds2_rbh(i-1)-potcm_r0hh)
                a1           = 1.0d0-exp(-arga1)
                a2           = potcm_aa2+potcm_aa3*exp(-arga2)
                force1_f1(i) = a1*exp(-a2*(bonds2_rc2h(i-4)-bonds_r0c2h)**2)
            enddo
            ! force1_f1 para el c2
            ! CRD 2012 eliminamos el calculo de potcm_a1c y potcm_a2c pq es constante
            ! lo incluimos en el CONST
            !        arga1c=potcm_aa1*rcb*rcb
            !        arga2c=potcm_aa4*(rcb-potcm_r0cb)*(rcb-potcm_r0cb)
            !        potcm_a1c=1.0d0-exp(-arga1c)
            !        potcm_a2c=potcm_aa2+potcm_aa3*exp(-arga2c)
            !        force1_f1(8)=potcm_a1c*exp(-potcm_a2c*(potcm_rcc-potcm_r0cc)**2)
        end subroutine ipforce

        !> @brief Evaluate switchf for the C2H6O_VBMM potential.
        subroutine switchf
            !  calculates switching functions: op1_s3,switch1_sphi,switch1_stheta
            !  and their derivatives op1_ds3,switch1_dsphi,switch1_dstheta
            implicit none
            double precision :: argmax, args1, args2, args3, argsc1, argsc2, argsc3, argsphi, argsphic, argstheta, argsthetac
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
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
            argmax = 19.0d0
            !  calculate ip1_s1 and ip1_ds1
            do i=1, 3
                args1 = potcm_a1s*(bonds2_rc1h(i)-bonds_r0c1h)*(bonds2_rc1h(i)-potcm_b1s)**8
                if( args1 .lt. argmax )then
                    ip1_s1(i) = 1.0d0-tanh(args1)
                else
                    ip1_s1(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                args1 = potcm_a1s*(bonds2_rc2h(i)-bonds_r0c2h)*(bonds2_rc2h(i)-potcm_b1s)**8
                if( args1 .lt. argmax )then
                    ip1_s1(i+4) = 1.0d0-tanh(args1)
                else
                    ip1_s1(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            !      ip1_s1(4)=1.0d0
            !      ip1_s1(8)=1.0d0
            argsc1 = potcm_a1scc*(potcm_rcc-potcm_r0cc)*(potcm_rcc-potcm_b1scc)**8
            if( argsc1 .lt. argmax )then
                ip1_s1(4) = 1.0d0-tanh(argsc1)
                ip1_s1(8) = ip1_s1(4)
            else
                ip1_s1(4) = 0.0d0
                ip1_s1(8) = ip1_s1(4)
            endif
            !  calculate ip1_s2 and ip1_ds2
            do i=1, 3
                args2 = potcm_a2s*(bonds2_rc1h(i)-bonds_r0c1h)*(bonds2_rc1h(i)-potcm_b2s)**6
                if( args2 .lt. argmax )then
                    ip1_s2(i) = 1.0d0-tanh(args2)
                else
                    ip1_s2(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                args2 = potcm_a2s*(bonds2_rc2h(i)-bonds_r0c2h)*(bonds2_rc2h(i)-potcm_b2s)**6
                if( args2 .lt. argmax )then
                    ip1_s2(i+4) = 1.0d0-tanh(args2)
                else
                    ip1_s2(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2005 ip1_s2(4)for the new atom
            ! CRD 2012 todos los S de los carbonos son 1
            !      ip1_s2(4)=1.0d0
            !      ip1_s2(8)=1.0d0
            argsc2 = potcm_a2scc*(potcm_rcc-potcm_r0cc)*(potcm_rcc-potcm_b2scc)**6
            if( argsc2 .lt. argmax )then
                ip1_s2(4) = 1.0d0-tanh(argsc2)
                ip1_s2(8) = ip1_s2(4)
            else
                ip1_s2(4) = 0.0d0
                ip1_s2(8) = ip1_s2(4)
            endif
            !  calculate op1_s3 and op1_ds3
            do i=1, 3
                args3 = potcm_a3s*(bonds2_rc1h(i)-bonds_r0c1h)*(bonds2_rc1h(i)-potcm_b3s)**2
                if( args3 .lt. argmax )then
                    op1_s3(i) = 1.0d0-tanh(args3)
                else
                    op1_s3(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                args3 = potcm_a3s*(bonds2_rc2h(i)-bonds_r0c2h)*(bonds2_rc2h(i)-potcm_b3s)**2
                if( args3 .lt. argmax )then
                    op1_s3(i+4) = 1.0d0-tanh(args3)
                else
                    op1_s3(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2005 op1_s3(4) for the new atom
            ! CRD 2012 todos los S de los carbonos son 1
            !      op1_s3(4)=1.0d0
            !      op1_s3(8)=1.0d0
            argsc3 = potcm_a3scc*(potcm_rcc-potcm_r0cc)*(potcm_rcc-potcm_b3scc)**2
            if( argsc3 .lt. argmax )then
                op1_s3(4) = 1.0d0-tanh(argsc3)
                op1_s3(8) = op1_s3(4)
            else
                op1_s3(4) = 0.0d0
                op1_s3(8) = op1_s3(4)
            endif
            ! JEG 2019. cambios por problemas en RPMD de F,Cl,H + C2H6!!
            op1_s3(4) = 1.0d0
            op1_s3(8) = 1.0d0
            ! JEG 2019
            !  calculate switch1_sphi and switch1_dsphi
            !  condition here is on the bondlength rch(i)
            !  st argsphi is lt approx 19.0d0
            do i=1, 3
                if( bonds2_rc1h(i) .lt. 3.8d0 )then
                    argsphi         = potcm_aphi*(bonds2_rc1h(i)-bonds_r0c1h)*exp(potcm_bphi*(bonds2_rc1h(i)- &
                        & potcm_cphi)**3)
                    switch1_sphi(i) = 1.0d0-tanh(argsphi)
                else
                    switch1_sphi(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                if( bonds2_rc2h(i) .lt. 3.8d0 )then
                    argsphi           = potcm_aphi*(bonds2_rc2h(i)-bonds_r0c2h)*exp(potcm_bphi*(bonds2_rc2h(i)- &
                        & potcm_cphi)**3)
                    switch1_sphi(i+4) = 1.0d0-tanh(argsphi)
                else
                    switch1_sphi(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            !      switch1_sphi(4)=1.0d0
            !      switch1_sphi(8)=1.0d0
            if( potcm_rcc .lt. 3.8d0 )then
                argsphic        = potcm_aphi*(potcm_rcc-potcm_r0cc)*exp(potcm_bphi*(potcm_rcc-potcm_cphi)**3)
                switch1_sphi(4) = 1.0d0-tanh(argsphic)
                switch1_sphi(8) = switch1_sphi(4)
            else
                switch1_sphi(4) = 0.0d0
                switch1_sphi(8) = switch1_sphi(4)
            endif
            !  calculate switch1_stheta and switch1_dstheta
            do i=1, 3
                if( bonds2_rc1h(i) .lt. 3.8d0 )then
                    argstheta = potcm_atheta*(bonds2_rc1h(i)-bonds_r0c1h)*exp(potcm_btheta*(bonds2_rc1h(i)- &
                        & potcm_ctheta)**3)
                    switch1_stheta(i) = 1.0d0-tanh(argstheta)
                else
                    switch1_stheta(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                if( bonds2_rc2h(i) .lt. 3.8d0 )then
                    argstheta   = potcm_atheta*(bonds2_rc2h(i)-bonds_r0c2h)*exp(potcm_btheta*(bonds2_rc2h(i)- &
                        & potcm_ctheta)**3)
                    switch1_stheta(i+4) = 1.0d0-tanh(argstheta)
                else
                    switch1_stheta(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            !      switch1_stheta(4)=1.0d0
            !      switch1_stheta(8)=1.0d0
            if( potcm_rcc .lt. 3.8d0 )then
                argsthetac        = potcm_atheta*(potcm_rcc-potcm_r0cc)*exp(potcm_btheta*(potcm_rcc-potcm_ctheta)**3)
                switch1_stheta(4) = 1.0d0-tanh(argsthetac)
                switch1_stheta(8) = switch1_stheta(4)
            else
                switch1_stheta(4) = 0.0d0
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
            implicit none
            double precision :: d1bc, d3bc, fact1, fact2, vdwc2
            integer          :: i, icount, ind
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            integer, parameter          :: n3tmmn = 27
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !      infocm_ref(1)='J. Espinosa-Garcia and J. C. Corchado'
            !      infocm_ref(2)='J. Chem. Phys., Vol. 112, p. 5731, 2000'
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
            !CALL POTINFO
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
            fact1       = 0.041840d0
            fact2       = 6.022045d0
            stret1_d1ch = stret1_d1ch*fact1
            stret1_d3ch = stret1_d3ch*fact1
            potcm_d1c1h = potcm_d1c1h*fact1
            potcm_d3c1h = potcm_d3c1h*fact1
            potcm_d1c2h = potcm_d1c2h*fact1
            potcm_d3c2h = potcm_d3c2h*fact1
            potcm_d1cc  = potcm_d1cc*fact1
            !      d3cc=d3cc*fact1
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
            potcm_a1s   = potcm_a1s*1.d-8
            potcm_a2s   = potcm_a2s*1.d-8
            potcm_a1scc = potcm_a1scc*1.d-8
            potcm_a2scc = potcm_a2scc*1.d-8
            vdwc2       = vdwc2*1.d+4
        end subroutine prepare_local_context

        !> @brief Evaluate POTINFO for the C2H6O_VBMM potential.
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

        !> @brief Evaluate ANCVRT for the C2H6O_VBMM potential.
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

        !> @brief Evaluate CARTOU for the C2H6O_VBMM potential.
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

        !> @brief Evaluate CARTTOR for the C2H6O_VBMM potential.
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

        !> @brief Evaluate EUNITZERO for the C2H6O_VBMM potential.
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

        !> @brief Evaluate RTOCART for the C2H6O_VBMM potential.
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

        !> @brief Evaluate DEDCOU for the C2H6O_VBMM potential.
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
end module chempotpy_c2h6o_vbmm_pes

!> @brief Expose C2H6O_VBMM through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_c2h6o_vbmm_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 9
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
