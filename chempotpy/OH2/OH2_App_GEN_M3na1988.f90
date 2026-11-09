!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_oh2_app_gen_m3na1988_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
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
        double precision              :: ambpcm_a
        double precision              :: ambpcm_beta
        double precision              :: ambpcm_d(3)
        double precision              :: ambpcm_de
        double precision              :: ambpcm_gamp(3)
        double precision              :: ambpcm_rbb
        double precision              :: ambpcm_req(3)
        double precision              :: ambpcm_rhh
        double precision              :: asycm_easyab
        double precision              :: asycm_easyac
        double precision              :: asycm_easybc
        double precision              :: infocm_cartnu(natom,3)
        integer                       :: infocm_icartr
        integer                       :: infocm_indexes(natom)
        integer                       :: infocm_irctnt
        integer                       :: infocm_mder
        integer                       :: infocm_msurf
        integer                       :: infocm_natoms
        character                     :: infocm_ref(5)
        double precision              :: pre22cm_a1dif
        double precision              :: pre22cm_a1inf
        double precision              :: pre22cm_a2dif
        double precision              :: pre22cm_a2inf
        double precision              :: pre22cm_alph1
        double precision              :: pre22cm_alph2
        double precision              :: pre22cm_b1dif
        double precision              :: pre22cm_b1inf
        double precision              :: pre22cm_b2dif
        double precision              :: pre22cm_b2inf
        double precision              :: pre22cm_bsp(20,3)
        double precision              :: pre22cm_csp(22,3)
        double precision              :: pre22cm_d1dif
        double precision              :: pre22cm_d1inf
        double precision              :: pre22cm_d2dif
        double precision              :: pre22cm_d2inf
        double precision              :: pre22cm_dasy
        integer                       :: pre22cm_iayes
        integer                       :: pre22cm_nphi
        double precision              :: pre22cm_phi(20)
        double precision              :: pre22cm_r1dif
        double precision              :: pre22cm_r1inf
        double precision              :: pre22cm_r1s
        double precision              :: pre22cm_r2dif
        double precision              :: pre22cm_r2inf
        double precision              :: pre22cm_r2s
        double precision              :: pre22cm_vsp
        double precision              :: pre2cm_dbcoor(2)
        double precision              :: precm_a2
        double precision              :: precm_asq
        double precision              :: precm_dcosdr(3)
        double precision              :: precm_de1(3)
        double precision              :: precm_de2(3)
        double precision              :: precm_dfdr(3)
        integer                       :: precm_iopsw
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
        integer                       :: iatom, idir, istate
        !initialize
        ambpcm_a         = 0.d0; ambpcm_beta = 0.d0; ambpcm_d = 0.d0; ambpcm_de = 0.d0; ambpcm_gamp = 0.d0; ambpcm_rbb = 0.d0
        ambpcm_req       = 0.d0; ambpcm_rhh = 0.d0; asycm_easyab = 0.d0; asycm_easyac = 0.d0; asycm_easybc = 0.d0
        infocm_cartnu    = 0.d0; infocm_icartr = 0; infocm_indexes = 0; infocm_irctnt = 0; infocm_mder = 0
        infocm_msurf     = 0; infocm_natoms = 0
        infocm_ref       = ' '
        pre22cm_a1dif    = 0.d0
        pre22cm_a1inf    = 0.d0
        pre22cm_a2dif    = 0.d0
        pre22cm_a2inf    = 0.d0
        pre22cm_alph1    = 0.d0
        pre22cm_alph2    = 0.d0
        pre22cm_b1dif    = 0.d0
        pre22cm_b1inf    = 0.d0
        pre22cm_b2dif    = 0.d0
        pre22cm_b2inf    = 0.d0
        pre22cm_bsp      = 0.d0
        pre22cm_csp      = 0.d0
        pre22cm_d1dif    = 0.d0
        pre22cm_d1inf    = 0.d0
        pre22cm_d2dif    = 0.d0
        pre22cm_d2inf    = 0.d0
        pre22cm_dasy     = 0.d0
        pre22cm_iayes    = 0
        pre22cm_nphi     = 0
        pre22cm_phi      = 0.d0
        pre22cm_r1dif    = 0.d0
        pre22cm_r1inf    = 0.d0
        pre22cm_r1s      = 0.d0
        pre22cm_r2dif    = 0.d0
        pre22cm_r2inf    = 0.d0
        pre22cm_r2s      = 0.d0
        pre22cm_vsp      = 0.d0
        pre2cm_dbcoor    = 0.d0
        precm_a2         = 0.d0
        precm_asq        = 0.d0
        precm_dcosdr     = 0.d0
        precm_de1        = 0.d0
        precm_de2        = 0.d0
        precm_dfdr       = 0.d0
        precm_iopsw      = 0
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
        usricm_nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0], shape(usricm_nasurf))
        usricm_nder    = 0
        usricm_nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero = 0.0d0
        infocm_icartr  = 3
        infocm_msurf   = 0
        infocm_mder    = 1
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 3
        ambpcm_d       = [106.56d0, 109.472d0, 106.56d0]
        ambpcm_beta    = 2.07942d0
        ambpcm_req     = [0.96966d0, 0.74144d0, 0.96966d0]
        ambpcm_a       = 0.26d0
        ambpcm_gamp    = [0.210820d0, 5.42027d0, -1.44235d0]
        ambpcm_rbb     = 0.74144d0
        precm_iopsw    = 2
        precm_a2       = 0.4d0
        pre22cm_d1inf  = 95.3d0
        pre22cm_r1inf  = .76d0
        pre22cm_b1inf  = 2.0202d0
        pre22cm_d2inf  = 92.4d0
        pre22cm_r2inf  = .98d0
        pre22cm_b2inf  = 2.41105d0
        pre22cm_a1inf  = 0.0d0
        pre22cm_a2inf  = 0.0d0
        pre22cm_r1s    = 2.21208d0
        pre22cm_r2s    = 2.21531d0
        pre22cm_vsp    = 95.3d0
        pre22cm_iayes  = 0
        pre22cm_alph1  = 1.6d0
        pre22cm_alph2  = 1.6d0
        pre22cm_nphi   = 11
        pre22cm_phi(1:11) = [0.0000000d+00, 0.1000000d+02, 0.2000000d+02, 0.3000000d+02, 0.3733000d+02, 0.4000000d+02, &
            & 0.5000000d+02, 0.6000000d+02, 0.7000000d+02, 0.8000000d+02, 0.9000000d+02]
        pre22cm_bsp(1:11,1) = [0.9477920d+02, 0.9316200d+02, 0.9035320d+02, 0.8517960d+02, 0.8272720d+02, &
            & 0.8304320d+02, 0.8855580d+02, 0.9003100d+02, 0.9129800d+02, 0.9184900d+02, 0.9209100d+02]
        pre22cm_bsp(1:11,2) = [0.1456500d+01, 0.1475500d+01, 0.1530500d+01, 0.1606400d+01, 0.1634900d+01, &
            & 0.1632300d+01, 0.1534600d+01, 0.1396600d+01, 0.1299800d+01, 0.1243800d+01, 0.1225600d+01]
        pre22cm_bsp(1:11,3) = [0.2041900d+01, 0.1987500d+01, 0.1905400d+01, 0.1672700d+01, 0.1554000d+01, &
            & 0.1539400d+01, 0.1749300d+01, 0.2014500d+01, 0.2227100d+01, 0.2348800d+01, 0.2383800d+01]
        pre22cm_csp(1:13,1) = [0.9477920d+02, -0.1478277d+00, -0.1389235d-03, -0.3580591d-03, 0.1450636d-02, &
            & 0.1148717d-02, -0.9458414d-03, -0.4983968d-02, 0.5921762d-02, -0.2814975d-02, 0.1001140d-02, - &
            & 0.2747850d-03, -0.5702495d-05]
        pre22cm_csp(1:13,2) = [0.1456500d+01, 0.1137130d-02, 0.7628700d-05, -0.9772201d-05, -0.1201120d-04, &
            & 0.2800151d-05, -0.4114549d-04, 0.7866376d-04, -0.1122307d-04, -0.1936703d-04, 0.6791209d-05, - &
            & 0.1039780d-04, 0.8032967d-05]
        pre22cm_csp(1:13,3) = [0.2041900d+01, -0.5887558d-02, 0.4475579d-05, -0.5455348d-04, 0.1230139d-03, - &
            & 0.6997501d-04, 0.2822753d-03, -0.4019424d-03, 0.1215636d-03, -0.1548133d-04, 0.9961687d-05, &
            & 0.1813458d-04, -0.1747243d-04]
        v             = 0.d0; g = 0.d0; d = 0.d0
        usricm_nder   = igrad
        usricm_nasurf = 1
        usricm_cart   = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                usricm_cart(iatom,idir) = x(iatom,idir)/legacy_bohr_to_angstrom
            enddo
        enddo
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*legacy_hartree_to_electronvolt
            enddo
        else if( igrad == 1 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*legacy_hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = usrocm_dgscart(iatom,idir)*51.422067
                enddo
            enddo
        else if( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
        contains
        !> @brief Evaluate prepare_local_context for the OH2_App_GEN_M3na1988 potential.
        subroutine prepare_local_context
            !
            !   System:          OH2   A'' surface
            !   Common name:     M3na, symmetrized M2 potential
            !   Functional form: Rotated-Morse-oscillator-spline (RMOS) potential
            !   Reference:       unpublished
            !   Cross reference: T. Joseph, D. G. Truhlar, and B. C. Garrett
            !                    J. Chem. Phys. 88, 6982-6990 (1988).
            !
            !   Interface:                     potlib2001
            !   Number of bodies:              3
            !   Number of derivatives:         1
            !   Number of electronic surfaces: 1
            !
            !   prepare_local_context must be called once before any calls to POT.
            !   The potential parameters are included in the block data subprogram P
            !   Coordinates, potential energy, and derivatives are passed
            !   The potential energy in the three asymptotic valleys are
            !   stored in the common block ASYCM:
            !                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the OH diatomic, with the
            !   OH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
            !   respectively.
            !   All the information passed through the common blocks PT1CM and ASYCM
            !   is in hartree atomic units.
            !
            !        This potential is written such that:
            !                       R(1) = R(O-H)
            !                       R(2) = R(H-H)
            !                       R(3) = R(H-O)
            !   The classical potential energy is set equal to zero for the O
            !   infinitely far from the equilibrium H2 diatomic.
            !
            !   The flags that indicate what calculations should be carried out in
            !   the potential routine are passed through the common block PT2CM:
            !   where:
            !        NASURF - which electronic states availalble
            !                 (1,1) = 1 as only gs state available
            !        NDER  = 0 => no derivatives should be calculated
            !        NDER  = 1 => calculate first derivatives
            !        NFLAG - these integer values can be used to flag options
            !                within the potential;
            !
            !
            !   Potential parameters' default settings
            !                  Variable            Default value
            !                  NDER                1
            !                  NFLAG(18)           6
            !
            !
            implicit none
            double precision :: amn, hpi
            integer          :: i
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            !
            !
            double precision, parameter :: one = 1.0d0
            !
            !
            !   Echo the name of the potential energy surface
            !
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !
            !      WRITE (NFLAG(18), 600)
            !
            hpi = 2.0d0*atan(one)
            !
            precm_iopsw = max(0, precm_iopsw)
            amn         = 0.0d0
            precm_a2    = max(amn, precm_a2)
            !
            !      WRITE (NFLAG(18),601) IOPSW, A2
            !
            !      CALL PREM2
            call prepare_local_context_2
            !
            precm_asq = precm_a2*precm_a2
            !
            !
            pt3cm_ezero(1) = ambpcm_d(2)
            !
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !
            infocm_ref(1) = 'T. Joseph, D. G. Truhlar, B. C. Garrett,'
            infocm_ref(2) = 'J. Chem. Phys. 88, 6982(1988)'
            !
            infocm_indexes(1) = 8
            infocm_indexes(2) = 1
            infocm_indexes(3) = 1
            !
            !
            !
            infocm_irctnt = 2
            !
            !      CALL POTINFO
            !
            call ancvrt
            !
        end subroutine prepare_local_context
        !> @brief Evaluate POT for the OH2_App_GEN_M3na1988 potential.
        subroutine pot
            !
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the OH diatomic, with the
            !   OH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
            !   respectively.
            !
            !        This potential is written such that:
            !                       R(1) = R(O-H)
            !                       R(2) = R(H-H)
            !                       R(3) = R(H-O)
            !   The classical potential energy is set equal to zero for the O
            !   infinitely far from the equilibrium H2 diatomic.
            !
            implicit none
            double precision :: bigr2, coschi, delv, dfdx, dgdx, dr13, f, frho2, g, hpi, omf, r12, r22, r32, rho, rho2, t, v1, v2, zero
            integer          :: i
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            !   Conversion for Angstroms to Bohr (Angstroms/Bohr)
            !
            double precision, parameter :: bohr = .52917706d0
            !
            !   Conversion for kcal to Hartree (kcal/Hartree)
            !
            double precision, parameter :: hartre = 627.5095d0
            !
            double precision, parameter :: tol = 1.d-5
            double precision, parameter :: one = 1.0d0
            !
            !
            call cartou
            call carttor
            !
            hpi = 2.0d0*atan(one)
            !
            !   Check the values of NASURF and NDER for validity.
            !
            !      IF (NASURF(1,1) .EQ. 0) THEN
            !         WRITE(NFLAG(18), 900) NASURF(1,1)
            !         STOP
            !      ENDIF
            !         IF (NDER .GT. 1) THEN
            !             WRITE (NFLAG(18), 910) NDER
            !             STOP 'POT 2'
            !         ENDIF
            !
            !  CHECK FOR SMALL VALUES OF R
            !
            if (pt1cm_r(1) .lt. 0.001d0 .or. pt1cm_r(2).lt. 0.001d0 .or.&
                &pt1cm_r(3) .lt. 0.001d0) then
            pt1cm_engygs = 1.d10
            if( usricm_nder .ne. 1 )then
                call eunitzero
                if( usricm_nder .ne. 0 )then
                    call rtocart
                    call dedcou
                endif
                return
            endif
            pt1cm_degsdr(1) = -pt1cm_engygs
            pt1cm_degsdr(2) = -pt1cm_engygs
            pt1cm_degsdr(3) = -pt1cm_engygs
            !
            call eunitzero
            if( usricm_nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
            !
            return
            endif
            !
            !    EVALUATE SWITCHING FUNCTION AND ITS DERIVATIVE
            !
            r12    = pt1cm_r(1)*pt1cm_r(1)
            r22    = pt1cm_r(2)*pt1cm_r(2)
            r32    = pt1cm_r(3)*pt1cm_r(3)
            bigr2  = 0.5d0*(r12 - 0.5d0*r22 + r32)
            bigr2  = max(zero, bigr2)
            rho2   = bigr2 + precm_asq
            frho2  = 4.0d0*rho2
            rho    = sqrt(rho2)
            dr13   = r12 - r32
            coschi = dr13/(2.0d0*pt1cm_r(2)*rho)
            !
            !    DERIVATIVES OF COS
            !
            if( usricm_nder .eq. 1 )then
                t               = 0.25d0/(pt1cm_r(2)*rho2)
                precm_dcosdr(2) = -coschi*(frho2 - r22)*t
                t               = t/rho
                precm_dcosdr(1) = pt1cm_r(1)*(frho2 - dr13)*t
                precm_dcosdr(3) = -pt1cm_r(3)*(frho2 + dr13)*t
            endif
            !
            !    EVALUATE SWITCHING FUNCTION
            !
            g    = coschi
            dgdx = 1.0d0
            if( precm_iopsw .gt. 0 )then
                do i = 1, precm_iopsw
                    t    = hpi*g
                    g    = sin(t)
                    dgdx = hpi*cos(t)*dgdx
                enddo
            endif
            f = 0.5d0*(1.0d0 - g)
            if( usricm_nder .eq. 1 )then
                dfdx = -0.5d0*dgdx
                do i = 1, 3
                    precm_dfdr(i) = dfdx*precm_dcosdr(i)
                enddo
            endif
            omf = 1.0d0 - f
            if( f .gt. 1.d-10 )then
                !
                !    EVALUATE POTENTIAL FOR O + H2 DIRECTION
                !
                !         CALL POTM2
                call pot2
                v1 = pt1cm_engygs
                if( usricm_nder .eq. 1 )then
                    do i = 1, 3
                        precm_de1(i) = pt1cm_degsdr(i)
                    enddo
                endif
            else
                v1 = 0.0d0
                if( usricm_nder .eq. 1 )then
                    do i = 1, 3
                        precm_de1(i) = 0.0d0
                    enddo
                endif
            endif
            if( abs(f-0.5d0) .lt. 1.d-10 )then
                v2 = pt1cm_engygs
                if( usricm_nder .eq. 1 )then
                    precm_de2(1) = pt1cm_degsdr(3)
                    precm_de2(2) = pt1cm_degsdr(2)
                    precm_de2(3) = pt1cm_degsdr(1)
                endif
            else if( omf .gt. 1.d-10 )then
                !
                !    INTERCHANGE R(1) AND R(3) AND RECOMPUTE POTENTIAL
                !
                t          = pt1cm_r(1)
                pt1cm_r(1) = pt1cm_r(3)
                pt1cm_r(3) = t
                !         CALL POTM2
                call pot2
                v2 = pt1cm_engygs
                if( usricm_nder .eq. 1 )then
                    precm_de2(1) = pt1cm_degsdr(3)
                    precm_de2(2) = pt1cm_degsdr(2)
                    precm_de2(3) = pt1cm_degsdr(1)
                endif
                t          = pt1cm_r(1)
                pt1cm_r(1) = pt1cm_r(3)
                pt1cm_r(3) = t
            else
                v2 = 0.0d0
                if( usricm_nder .eq. 1 )then
                    do i = 1, 3
                        precm_de2(i) = 0.0d0
                    enddo
                endif
            endif
            !
            pt1cm_engygs = f*v1 + omf*v2
            if( usricm_nder .ne. 1 )then
                call eunitzero
                if( usricm_nder .ne. 0 )then
                    call rtocart
                    call dedcou
                endif
                return
            endif
            delv = v1 - v2
            do i = 1, 3
                pt1cm_degsdr(i) = f*precm_de1(i) + omf*precm_de2(i) + delv*precm_dfdr(i)
            enddo
            !
            !
            call eunitzero
            if( usricm_nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
            !
        end subroutine pot
        !> @brief Evaluate prepare_local_context_2 for the OH2_App_GEN_M3na1988 potential.
        subroutine prepare_local_context_2
            !
            !   The subprogram prepare_local_context_2 determines the bending correction.
            !
            !   The potential parameters for the collinear potential are
            !   initialized in the block data subprogram PTPACM.  The
            !   potential parameters for the bending potential are assigned
            !   in the subprogram prepare_local_context_22.  The subprograms prepare_local_context_2 and
            !   prepare_local_context_22 do not need to be called in the MAIN program.
            !
            implicit none
            integer :: it
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            double precision :: nab, nab1, nbc, nbc1
            !
            !   Conversion for Angstroms to Bohr (Angstroms/Bohr)
            !
            double precision, parameter :: bohr = .52917706d0
            !
            !   Conversion for kcal to Hartree (kcal/Hartree)
            !
            double precision, parameter :: hartre = 627.5095d0
            !
            !  Vaule of tolerance
            !
            double precision, parameter :: tol = 1.d-5
            !
            !
            !   Assign parameters for the bending potential
            !
            !      CALL PREM22
            !
            call prepare_local_context_22
            !
            !   Echo the potential parameters for the bending correction.
            !
            !      WRITE (NFLAG(18),600) REQ,D
            !      WRITE (NFLAG(18),601) BETA,A,RBB,GAMP
            !
            !       CONVERT TO ATOMIC UNITS
            !
            do it = 1, 3
                ambpcm_req(it) = ambpcm_req(it)/bohr
                ambpcm_d(it)   = ambpcm_d(it)/hartre
            enddo
            ambpcm_beta = ambpcm_beta * bohr
            ambpcm_de   = ambpcm_d(1)
            ambpcm_a    = ambpcm_a/bohr
            ambpcm_rhh  = 0.74144d0 / bohr
            ambpcm_rbb  = ambpcm_rbb /bohr
            !
            !   Initialize the energy in the asymptotic valleys
            !
            asycm_easyab = ambpcm_d(1)
            asycm_easybc = ambpcm_d(2)
            asycm_easyac = ambpcm_d(3)
            !
            !
        end subroutine prepare_local_context_2
        !> @brief Evaluate POT2 for the OH2_App_GEN_M3na1988 potential.
        subroutine pot2
            !
            implicit none
            double precision :: bcoor, c, dac, dex1, dex2, dgam, ex, ex1, ex2, gamma, raccol, ractmp, rob, stuff, t
            integer          :: it
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            double precision :: nab, nab1, nbc, nbc1
            !
            !   Conversion for Angstroms to Bohr (Angstroms/Bohr)
            !
            double precision, parameter :: bohr = .52917706d0
            !
            !   Conversion for kcal to Hartree (kcal/Hartree)
            !
            double precision, parameter :: hartre = 627.5095d0
            !
            !  Vaule of tolerance
            !
            double precision, parameter :: tol = 1.d-5
            !
            !
            if( pt1cm_r(1) .gt. tol .and. pt1cm_r(2) .gt. tol .and. pt1cm_r(3) .gt. tol )then
                raccol = pt1cm_r(1) + pt1cm_r(2)
                bcoor  = 0.0d0
                if( usricm_nder .eq. 1 )then
                    pre2cm_dbcoor(1) = 0.0d0
                    pre2cm_dbcoor(2) = 0.0d0
                endif
                dac = 0.0d0
                rob = 1.0d0
                if( abs(ambpcm_rhh-ambpcm_rbb) .gt. 1.d-10 )then
                    !
                    !     FIRST CALCULATE SOME USEFUL NUMBERS
                    !     NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
                    !     AB AND BC RESPECTIVELY. (EQUATION 2.2)
                    !
                    nab = exp((ambpcm_req(1) - pt1cm_r(1))/ambpcm_a)
                    if( nab .lt. 1.d-15 ) nab = 1.d-15
                    nbc = exp((ambpcm_req(2) - pt1cm_r(2))/ambpcm_a)
                    if( nbc .lt. 1.d-15 ) nbc = 1.d-15
                    !
                    !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
                    !       SEE EQUATION 2.7A, 2.7B
                    !
                    c = 1.d0 - nab/nbc
                    if( abs(c) .lt. 1.d-14 )then
                        nab1 = .5d0
                        nbc1 = .5d0
                    else
                        c    = c/nab
                        nab1 = (2.d0 + c - sqrt(4.d0+c*c))/(2.d0*c)
                        nbc1 = 1.d0 - nab1
                    endif
                    !
                    !          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
                    !          FIRST TRAP OUT ANY ZERO ARGUMENTS
                    !
                    if( nab1*nbc1 .le. 0.0d0 )then
                        rob = 1.
                    else
                        stuff = ambpcm_a * log(nab1*nbc1)
                        t     = ambpcm_rhh - stuff
                        rob   = 1.0d0
                        if( abs(t) .gt. 1.d-15 ) rob = (ambpcm_rbb - stuff)/t
                    endif
                endif
                !
                !        CALCULATE BENDING CORRECTION
                !        EVALUATE GAMMA AND DERIVATIVE
                !
                ex    = exp(raccol*(ambpcm_gamp(2) + ambpcm_gamp(3)*raccol))
                gamma = ambpcm_gamp(1) + ex
                if( usricm_nder .eq. 1 ) dgam = (ambpcm_gamp(2) + 2.d0*ambpcm_gamp(3)*raccol)*ex
                gamma = gamma*ambpcm_de
                if( usricm_nder .eq. 1 ) dgam = dgam*ambpcm_de
                !
                !           CALCULATE V(R(3)) TERM
                !
                ex  = exp(-ambpcm_beta*(pt1cm_r(3)-ambpcm_req(3))/rob)
                ex1 = ex*(1.d0 + 0.5d0*ex)
                if( usricm_nder .eq. 1 ) dex1 = ex*(1.d0 + ex)
                !
                !           CALCULATE V(R(1)+R(2)) TERM
                !
                ex  = exp(ambpcm_beta*(ambpcm_req(3)-pt1cm_r(2)-pt1cm_r(1))/rob)
                ex2 = ex*(1.d0 + 0.5d0*ex)
                if( usricm_nder .eq. 1 ) dex2 = ex*(1.d0 + ex)
                !
                !        HERE IS THE BENDING CORRECTION
                !
                bcoor = gamma*(ex1 - ex2)
                !
                !        WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
                !        OF BCOOR WITH RESPECT TO R(1), R(2), AND R(3).
                !        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
                !
                gamma = gamma*ambpcm_beta/rob
                if( usricm_nder .eq. 1 )then
                    dac = -gamma*dex1
                    !          CALCULATE CORRECTIONS TO DAB, DBC
                    pre2cm_dbcoor(1) = gamma*dex2 + dgam*(ex1-ex2)
                    pre2cm_dbcoor(2) = pre2cm_dbcoor(1)
                endif
                !
                !       NOW CALCULATE THE COLLINEAR ENERGY
                !       STORE R(3) IN RACTMP THEN SET R(3) TO
                !       THE COLLINEAR GEOMETRY AND CALL POTM22
                !
                ractmp     = pt1cm_r(3)
                pt1cm_r(3) = pt1cm_r(2) + pt1cm_r(1)
                !
                !         CALL POTM22
                call pot22
                !
                pt1cm_r(3) = ractmp
                !
                pt1cm_engygs = pt1cm_engygs + bcoor
                !
                if( usricm_nder .eq. 1 )then
                    do it = 1, 2
                        pt1cm_degsdr(it) = pre2cm_dbcoor(it) + pt1cm_degsdr(it) + pt1cm_degsdr(3)
                    enddo
                    pt1cm_degsdr(3) = dac
                endif
            else
                pt1cm_engygs = 1.d30
                if( usricm_nder .ne. 1 ) return
                pt1cm_degsdr(1)                           = 0.0d0
                if( pt1cm_r(1) .le. tol ) pt1cm_degsdr(1) = -pt1cm_engygs
                pt1cm_degsdr(2)                           = 0.0d0
                if( pt1cm_r(1) .le. tol ) pt1cm_degsdr(2) = -pt1cm_engygs
                pt1cm_degsdr(3)                           = 0.0d0
                if( pt1cm_r(1) .le. tol ) pt1cm_degsdr(3) = -pt1cm_engygs
            endif
            !
        end subroutine pot2
        !> @brief Evaluate prepare_local_context_22 for the OH2_App_GEN_M3na1988 potential.
        subroutine prepare_local_context_22
            !
            !     COLLINEAR A + BC POTENTIAL FUNCTION
            !     RMO-SPLINE TYPE FIT
            !     COORDINATES AND ENERGIES PASSED TO POTENTIAL IN ATOMIC UNITS
            !
            implicit none
            double precision :: aer, bet, dis, req2, t, t1
            integer          :: i, ii, j, nnn, nphip2
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            double precision, parameter :: detrad = 0.0174532d0
            double precision, parameter :: ckcal  = 627.5095d0
            double precision, parameter :: bohr   = 0.529177d0
            !
            !
            character(len=20) :: title
            !
            !  TITLE
            !
            title = '   OH2 MAB 9-12-79 '
            !
            !      WRITE (NFLAG(18), 599)
            !
            !  ASSIGN COLLINEAR PARAMETERS (UNITS OF KCAL, ANGSTROM, DEGREES)
            !  IN BLOCK DATA SUBROUTINE
            !
            pre22cm_a1inf = 0.0d0
            pre22cm_a2inf = 0.0d0
            pre22cm_iayes = 0
            !
            !  ASSIGN PARAMETERS FOR THE ASYMPTOTIC REACTANT AND PRODUCT REGIONS
            !  IN BLOCK DATA SUBROUTINE
            !
            if( pre22cm_a2inf .eq. 0.0d0 ) pre22cm_a2inf = pre22cm_d1inf-pre22cm_d2inf
            !
            !  ASSIGN THE PIVOT POINT, POTENTIAL AT THE PIVOT,
            !  FLAG, AND ASYMPTOTIC EXPONENTIAL PARAMETERS
            !  IN BLOCK DATA SUBROUTINE.  ECHO.
            !
            !      WRITE (NFLAG(18), 600) TITLE
            !      WRITE (NFLAG(18), 601) D1INF,R1INF,B1INF,D2INF,R2INF,B2INF,A1INF,
            !     2INF
            !      WRITE (NFLAG(18), 602) R1S,R2S,VSP,ALPH1,ALPH2
            !
            !      IF(IAYES.NE.0) WRITE (NFLAG(18), 603)
            !      IF(IAYES.EQ.0) WRITE (NFLAG(18), 604)
            !
            nphip2 = pre22cm_nphi + 2
            !
            !  ECHO POTENTIAL PARAMTERS ASSIGNED BY BLOCK DATA SUBROUTINE
            !
            !      WRITE (NFLAG(18), 605) NPHI
            !      WRITE (NFLAG(18), 610)
            !
            !      DO 30 I = 1,NPHI
            !         WRITE (NFLAG(18), 606) I,PHI(I),(BSP(I,J),CSP(I,J),J=1,3)
            !   30 CONTINUE
            !
            !      WRITE (NFLAG(18), 611)
            !      DO 40 I=1,2
            !         K = I + NPHI
            !         WRITE (NFLAG(18), 607) K,(CSP(K,J),J=1,3)
            !   40 CONTINUE
            !
            !  CONVERT TO ATOMIC UNITS
            !
            pre22cm_d1inf = pre22cm_d1inf/ckcal
            pre22cm_a1inf = pre22cm_a1inf/ckcal
            pre22cm_r1inf = pre22cm_r1inf/bohr
            pre22cm_b1inf = pre22cm_b1inf*bohr
            pre22cm_d2inf = pre22cm_d2inf/ckcal
            pre22cm_a2inf = pre22cm_a2inf/ckcal
            pre22cm_r2inf = pre22cm_r2inf/bohr
            pre22cm_b2inf = pre22cm_b2inf*bohr
            pre22cm_r1s   = pre22cm_r1s/bohr
            pre22cm_r2s   = pre22cm_r2s/bohr
            pre22cm_vsp   = pre22cm_vsp/ckcal
            pre22cm_alph1 = pre22cm_alph1*bohr
            pre22cm_alph2 = pre22cm_alph2*bohr
            do i = 1, pre22cm_nphi
                pre22cm_bsp(i,1) = pre22cm_bsp(i,1)/ckcal
                pre22cm_bsp(i,2) = pre22cm_bsp(i,2)/bohr
                pre22cm_bsp(i,3) = pre22cm_bsp(i,3)*bohr
            enddo
            nnn = pre22cm_nphi + 2
            do i = 1, nnn
                pre22cm_csp(i,1) = pre22cm_csp(i,1)/ckcal
                pre22cm_csp(i,2) = pre22cm_csp(i,2)/bohr
                pre22cm_csp(i,3) = pre22cm_csp(i,3)*bohr
            enddo
            !
            !  CONVERT FROM DEGREES TO RADIANS
            !
            t = 1.0d0/detrad
            do j = 1, 3
                pre22cm_csp(2,j) = pre22cm_csp(2,j)*t
            enddo
            t = t**3
            do ii = 1, pre22cm_nphi
                pre22cm_phi(ii) = pre22cm_phi(ii)*detrad
                i               = ii + 2
                do j = 1, 3
                    pre22cm_csp(i,j) = pre22cm_csp(i,j)*t
                enddo
            enddo
            !
            !  SET UP CONSTANTS
            !
            !     REACTANT ASYMPTOTE
            !
            dis           = pre22cm_bsp(1,1)
            req2          = pre22cm_bsp(1,2)
            bet           = pre22cm_bsp(1,3)
            t1            = 1.0d0 - exp(-bet*req2)
            aer           = pre22cm_vsp - dis*t1*t1
            req2          = pre22cm_r2s - req2
            pre22cm_d1dif = dis - pre22cm_d1inf
            pre22cm_r1dif = req2 - pre22cm_r1inf
            pre22cm_b1dif = bet - pre22cm_b1inf
            pre22cm_a1dif = aer - pre22cm_a1inf
            !
            !     PRODUCT ASYMPTOTE
            !
            dis           = pre22cm_bsp(pre22cm_nphi,1)
            req2          = pre22cm_bsp(pre22cm_nphi,2)
            bet           = pre22cm_bsp(pre22cm_nphi,3)
            t1            = 1.0d0 - exp(-bet*req2)
            aer           = pre22cm_vsp - dis*t1*t1
            req2          = pre22cm_r1s - req2
            pre22cm_d2dif = dis - pre22cm_d2inf
            pre22cm_r2dif = req2 - pre22cm_r2inf
            pre22cm_b2dif = bet - pre22cm_b2inf
            pre22cm_a2dif = aer - pre22cm_a2inf
            !
            !  DEFINE ASYMPTOTIC ENERGY AND ZERO OF ENERGY
            !
            pre22cm_dasy = pre22cm_d1inf
            !
        end subroutine prepare_local_context_22
        !> @brief Evaluate POT22 for the OH2_App_GEN_M3na1988 potential.
        subroutine pot22
            !
            implicit none
            double precision :: aer, angle, bet, daer, dbet, ddis, delr, dis, dphd1, dphd2, dr1, dr2, dreq, drrd1, drrd2, dvdp, req2, rr, scrap(pre22cm_nphi+2), t1, t2, t3
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            double precision, parameter :: detrad = 0.0174532d0
            double precision, parameter :: ckcal  = 627.5095d0
            double precision, parameter :: bohr   = 0.529177d0
            !
            !
            character(len=20) :: title
            !
            !      R1 = R(1)
            !      R2 = R(2)
            !      R3 = R(3)
            if( pt1cm_r(1) .gt. pre22cm_r1s )then
                if( pt1cm_r(2) .gt. pre22cm_r2s )then
                    !
                    !   3-ATOM BREAK UP REGION
                    !
                    pt1cm_engygs = pre22cm_dasy
                    if( usricm_nder .eq. 1 )then
                        pt1cm_degsdr(1) = 0.0d0
                        pt1cm_degsdr(2) = 0.0d0
                    endif
                else
                    !
                    !   REACTANT ASYMPTOTIC REGION
                    !
                    dr1 = pt1cm_r(1) - pre22cm_r1s
                    if( dr1*pre22cm_alph1 .lt. -70.0d0 ) &
                        & write (usricm_nflag(18), '(/,2X,T5,''R(1), R(2), R(3), DELR, ALPH1:'', /,2X,T10,5(1PE13.5,1X))') pt1cm_r(1), pt1cm_r(2), pt1cm_r(3), dr1, pre22cm_alph1
                    t1 = exp(-dr1*pre22cm_alph1)
                    t2 = -pre22cm_alph1*t1
                    !
                    !.....INTERPOLATE FOR MORSE PARAMETERS
                    !
                    dis  = pre22cm_d1inf + pre22cm_d1dif*t1
                    ddis = pre22cm_d1dif*t2
                    req2 = pre22cm_r1inf + pre22cm_r1dif*t1
                    dreq = pre22cm_r1dif*t2
                    bet  = pre22cm_b1inf + pre22cm_b1dif*t1
                    dbet = pre22cm_b1dif*t2
                    if( pre22cm_iayes .eq. 0 )then
                        aer  = pre22cm_vsp - dis
                        daer = -ddis
                    else
                        aer  = pre22cm_a1inf + pre22cm_a1dif*t1
                        daer = pre22cm_a1dif*t2
                    endif
                    delr = pt1cm_r(2) - req2
                    if( delr*bet .lt. -70.0d0 ) &
                        & write (usricm_nflag(18), '(/,2X,T5,''R(1), R(2), R(3), DR2, BET:'', /,2X,T10,5(1PE13.5,1X))') pt1cm_r(1), pt1cm_r(2), pt1cm_r(3), delr, bet
                    t1           = exp(-bet*delr)
                    t2           = 1.0d0 - t1
                    pt1cm_engygs = aer + dis*t2*t2
                    if( usricm_nder .eq. 1 )then
                        t3              = 2.0d0*dis*t1
                        pt1cm_degsdr(1) = daer + t2*(ddis*t2 - t3*(bet*dreq - dbet*delr))
                        pt1cm_degsdr(2) = t3*bet*t2
                    endif
                endif
            else if( pt1cm_r(2) .gt. pre22cm_r2s )then
                !
                !   PRODUCT ASYMPTOTIC REGION
                !
                dr2 = pt1cm_r(2) - pre22cm_r2s
                if( dr2*pre22cm_alph2 .lt. -70.d0 ) &
                    & write (usricm_nflag(18), '(/,2X,T5,''R(1), R(2), R(3), DELR, ALPH2:'', /,2X,T10,5(1PE13.5,1X))') pt1cm_r(1), pt1cm_r(2), pt1cm_r(3), dr2, pre22cm_alph2
                !
                !.....INTERPOLATE FOR MORSE PARAMETERS
                !
                t1   = exp(-dr2*pre22cm_alph2)
                t2   = -pre22cm_alph2*t1
                dis  = pre22cm_d2inf + pre22cm_d2dif*t1
                ddis = pre22cm_d2dif*t2
                req2 = pre22cm_r2inf + pre22cm_r2dif*t1
                dreq = pre22cm_r2dif*t2
                bet  = pre22cm_b2inf + pre22cm_b2dif*t1
                dbet = pre22cm_b2dif*t2
                if( pre22cm_iayes .eq. 0 )then
                    aer  = pre22cm_vsp - dis
                    daer = -ddis
                else
                    aer  = pre22cm_a2inf + pre22cm_a2dif*t1
                    daer = pre22cm_a2dif*t2
                endif
                delr = pt1cm_r(1) - req2
                if( delr*bet .lt. -70.0d0 ) &
                    & write (usricm_nflag(18), '(/,2X,T5,''R(1), R(2), R(3), DR2, BET:'', /,2X,T10,5(1PE13.5,1X))') pt1cm_r(1), pt1cm_r(2), pt1cm_r(3), delr, bet
                t1           = exp(-bet*delr)
                t2           = 1.0d0 - t1
                pt1cm_engygs = aer + dis*t2*t2
                if( usricm_nder .eq. 1 )then
                    t3              = 2.0d0*dis*t1
                    pt1cm_degsdr(2) = daer + t2*(ddis*t2 - t3*(bet*dreq - dbet*delr))
                    !
                    pt1cm_degsdr(1) = t3*bet*t2
                endif
            else
                !
                !   SPLINE INTERPOLATION REGION
                !
                dr1   = pre22cm_r1s - pt1cm_r(1)
                dr2   = pre22cm_r2s - pt1cm_r(2)
                t1    = dr1*dr1 + dr2*dr2
                rr    = sqrt(t1)
                drrd1 = -dr1/rr
                drrd2 = -dr2/rr
                angle = atan(dr1/dr2)
                dphd1 = -dr2/t1
                dphd2 = dr1/t1
                dis   = potspl(3, pre22cm_nphi, pre22cm_phi, pre22cm_bsp(1,1), pre22cm_csp(1,1), scrap, angle)
                ddis  = potspl(4, pre22cm_nphi, pre22cm_phi, pre22cm_bsp(1,1), pre22cm_csp(1,1), scrap, angle)
                req2  = potspl(3, pre22cm_nphi, pre22cm_phi, pre22cm_bsp(1,2), pre22cm_csp(1,2), scrap, angle)
                dreq  = potspl(4, pre22cm_nphi, pre22cm_phi, pre22cm_bsp(1,2), pre22cm_csp(1,2), scrap, angle)
                bet   = potspl(3, pre22cm_nphi, pre22cm_phi, pre22cm_bsp(1,3), pre22cm_csp(1,3), scrap, angle)
                dbet  = potspl(4, pre22cm_nphi, pre22cm_phi, pre22cm_bsp(1,3), pre22cm_csp(1,3), scrap, angle)
                if( pre22cm_iayes .eq. 0 )then
                    aer  = pre22cm_vsp - dis
                    daer = -ddis
                else
                    t1   = exp(-bet*req2)
                    t2   = 1.0d0 - t1
                    aer  = pre22cm_vsp - dis*t2*t2
                    daer = -t2*(ddis*t2 + 2.0d0*dis*t1*(bet*dreq + req2*dbet))
                endif
                if( bet*(rr-req2) .gt. 70 ) &
                    & write (usricm_nflag(18), '(/,2X,T5,''R(1), R(2), R(3), RR, ANGLE, DIS, REQ2, BET:'', /,(2X,T5,1PE20.10,T30,1PE' // &
                    '20.10,T55,1PE20.10))') pt1cm_r(1), pt1cm_r(2), pt1cm_r(3), rr, angle, dis, &
                    & req2, bet
                delr         = rr - req2
                t1           = exp(bet*delr)
                t2           = 1.0d0 - t1
                pt1cm_engygs = aer + dis*t2*t2
                if( usricm_nder .eq. 1 )then
                    t3              = 2.0d0*dis*t1
                    dvdp            = daer + t2*(ddis*t2 - t3*(dbet*delr - bet*dreq))
                    t3              = t3*t2*bet
                    pt1cm_degsdr(1) = dvdp*dphd1 - t3*drrd1
                    pt1cm_degsdr(2) = dvdp*dphd2 - t3*drrd2
                endif
            endif
            pt1cm_engygs = pt1cm_engygs
            if( usricm_nder .ne. 1 ) return
            pt1cm_degsdr(1) = pt1cm_degsdr(1)
            pt1cm_degsdr(2) = pt1cm_degsdr(2)
            pt1cm_degsdr(3) = 0.0d0
            !
            !
        end subroutine pot22
        !> @brief Evaluate POTSPL for the OH2_App_GEN_M3na1988 potential.
        function potspl( isw, nn, x, y, c, d, xpt )
            implicit none
            double precision :: c, cfact, d, deriv, dfact, pivot, potspl, supd, v, v2, x, xpt, y, z
            integer          :: i, isw, j, n, nn, np1, np2
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            dimension x(nn), y(nn), c(nn+2), d(nn+2)
            !
            !  THIS IS A SUBROUTINE FOR FITTING DATA WITH A CUBIC SPLINE
            !  POLYNOMIAL AND EVALUATING THAT POLYNOMIAL AT A GIVEN POINT
            !  OR ITS DERIVATIVE AT A GIVEN POINT
            !
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
            !
            !     NN ... THE NUMBER OF KNOTS (DATA POINTS) TO WHICH THE SPLINE IS TO
            !            BE FITTED
            !
            !     X,Y ... THE ARRAYS DEFINING THE KNOTS.  THE X-VALUES MUST BE IN
            !             INCREASING ORDER.  THE ARRAYS MUST BE DIMENSIONED AT LEAST
            !             NN.
            !
            !     C ... THE ARRAY THAT CONTAINS THE CUBIC SPLINE COEFFICIENTS.
            !           MUST BE DIMENSIONED AT LEAST NN+2 .
            !
            !     D ... A WORK SPACE.  MUST BE DIMENSIONED AT LEAST NN+2 .
            !
            !     XPT ... THE POINT AT WHICH THE INTERPOLATION IS DESIRED (IF ISW IS
            !              SET TO 2).  THE VALUE OF POTSPL IS SET TO THE
            !              INTERPOLATED VALUE.
            !
            !   Output from this subprogram is written to UNIT NFLAG(18)
            !
            ! *****  USER NOTES  *****
            !
            !     INTERPOLATION INVOLVES AT LEAST TWO STEPS .......
            !
            !       A.  CALL POTSPL WITH THE KNOTS.  THIS SETS UP THE
            !           COEFFICIENT ARRAY C.
            !           EG.  DUMY=POTSPL(1,NN,X,Y,C,D,XPT)
            !
            !       B.  CALL POTSPL WITH THE ARRAY C WHICH WAS DEFINED BY THE
            !           PREVIOUS CALL AND WILL BE USED TO FIND THE VALUE AT THE
            !           POINT 'XPT' .
            !           EG.   VALUE=POTSPL(2,NN,X,Y,C,D,XPT)
            !
            !     STEP 'A' NEED BE EXECUTED ONLY ONCE FOR A GIVEN SET OF KNOTS.
            !     STEP B MAY BE EXECUTED AS MANY TIMES AS NECESSARY.
            !
            n   = nn
            np1 = n+1
            np2 = n+2
            z   = xpt
            if( isw .eq. 1 )then
                c(1)   = y(1)
                d(1)   = 1.0d0
                c(np1) = 0.0d0
                d(np1) = 0.0d0
                c(np2) = 0.0d0
                d(np2) = 0.0d0
                do i = 2, n
                    c(i) = y(i)-y(1)
                    d(i) = x(i)-x(1)
                enddo
                do i = 3, np2
                    if( d(i-1) .eq. 0.0d0 )then
                        write(usricm_nflag(18), '(1X,5X,'' ***** ERROR IN POTSPL ... DIVIDE FAULT *****'')')
                        stop 'POTSPL 1'
                    endif
                    pivot = 1.0d0/d(i-1)
                    if( i .lt. np2 )then
                        supd = x(i-1)-x(i-2)
                        if( supd .lt. 0.0d0 )then
                            write(usricm_nflag(18), '(1X,5X,''***** ERROR IN POTSPL ... UNORDERED X-VALUES *****'')')
                            stop 'POTSPL 2'
                        endif
                        supd = supd*supd*supd
                    else
                        supd = 1.0d0
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
                        d(np1) = 1.0d0-d(np1)*dfact
                    endif
                    c(np2) = c(np2)-d(np2)*cfact
                    d(np2) = x(i-2)-d(np2)*dfact
                enddo
                do i = 1, n
                    j = np2-i
                    if( j .eq. np1 )then
                        v = 1.0d0
                    else
                        v = x(j)-x(j-1)
                        v = v*v*v
                    endif
                    if( d(j+1) .eq. 0.0d0 )then
                        write(usricm_nflag(18), '(1X,5X,'' ***** ERROR IN POTSPL ... DIVIDE FAULT *****'')')
                        stop 'POTSPL 3'
                    endif
                    c(j+1) = c(j+1)/d(j+1)
                    c(j)   = c(j)-c(j+1)*v
                enddo
                if( d(2) .eq. 0.0d0 )then
                    write(usricm_nflag(18), '(1X,5X,'' ***** ERROR IN POTSPL ... DIVIDE FAULT *****'')')
                    stop 'POTSPL 4'
                endif
                c(2) = c(2)/d(2)
            else if( isw .eq. 2 )then
                potspl = c(1) + c(2)*(z-x(1))
                do i = 1, n
                    v = z-x(i)
                    if( v .le. 0 ) return
                    potspl = potspl + c(i+2)*v**3
                enddo
            else if( isw .eq. 3 )then
                potspl = c(1) + c(2)*(z-x(1))
                deriv  = c(2)/3.d0
                do i = 1, n
                    v = z-x(i)
                    if( v .le. 0 ) return
                    v2     = v*v
                    potspl = potspl + c(i+2)*v2*v
                    deriv  = deriv + c(i+2)*v2
                enddo
            else if( isw .eq. 4 )then
                potspl = 3.d0*deriv
            endif
        end function potspl
        !> @brief Evaluate POTINFO for the OH2_App_GEN_M3na1988 potential.
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
                write(usricm_nflag(18), '(2x,''Cartesian coordinates are supplied by'',/, 2x,''the user in the array CART.'',//)')
                write(usricm_nflag(18), '(2x,''Provide cartesian coordinates in the'',/, 2x,''following order using the array CART'',' // &
                    '//, 2x,'' CART(1,1)...CART(1,3)   => ATOM 1'',/, 2x,'' CART(2,1)...CART(2,3)   => AT' // &
                    'OM 2'',/, 2x,'' CART(3,1)...CART(3,3)   => ATOM 3'',/, 2x,'' CART(N,1)...CART(N,3)  ' // &
                    ' => ATOM N'',/, 2x,''CART(25,1)...CART(25,3)  => ATOM 25'',/)')
                write(usricm_nflag(18), '(2x,''If the user wishes to relabel the atoms,'',/, 2x,''set the variable IREORDER equal to ' // &
                    '1'',/, 2x,''in the PARAMETER statement.  The user'',/, 2x,''must also specify the ne' // &
                    'w labeling'',/, 2x,''scheme.  This is done using the array'',/, 2x,''NULBL in the fo' // &
                    'llowing manner:'',//, 2x,''NULBL(i) = j'',/, 2x,''where:  i => old index'',/, 2x,'' ' // &
                    '       j => new index'',//)')
                write(usricm_nflag(18), '(2x,''Cartesian coordinates can be provided to'',/, 2x,''the potential routine in a variety ' // &
                    'of units.'',/, 2x,''The input units will be converted to Bohr'',/, 2x,''based on the' // &
                    ' following values of the NFLAG'',/, 2x,''variable:'',//, 2x,''NFLAG(1)  =  1  =>  CA' // &
                    'RTESIANS IN BOHR (no'',/, 2x,''                    conversion required)'',/, 2x,''NF' // &
                    'LAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS'',//)')
                write(usricm_nflag(18), '(2x,''The value of the energy and derivatives'',/, 2x,''(if computed) can be reported in a v' // &
                    'ariety'',/, 2x,''units.  A units conversion will take place'',/, 2x,''as specified b' // &
                    'y the following values of the'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(2) = 1 =>  ' // &
                    'ENERGIES REPORTED IN HARTEEE'',/, 2x,''NFLAG(2) = 2 =>  ENERGIES REPORTED IN mHARTRE' // &
                    'E'',/, 2x,''NFLAG(2) = 3 =>  ENERGIES REPORTED IN eV'',/, 2x,''NFLAG(2) = 4 =>  ENER' // &
                    'GIES REPORTED IN kcal/mol'',/, 2x,''NFLAG(2) = 5 =>  ENERGIES REPORTED IN cm**-1'',/' // &
                    '/)')
                write(usricm_nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
                    'he'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(1)=1 & NFLAG(2)=1 => DERIVATIVES REPOR' // &
                    'TED IN'',/, 2x,''                           HARTEEE/BOHR'',/, 2x,''NFLAG(1)=1 & NFLA' // &
                    'G(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/BOHR' // &
                    ''',/, 2x,''NFLAG(1)=1 & NFLAG(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''            ' // &
                    '               eV/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=4 => DERIVATIVES REPORTED IN' // &
                    ''',/, 2x,''                           kcal/mol/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=' // &
                    '5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm**-1/BOHR'',//)')
                write(usricm_nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
                    'he'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(1)=2 & NFLAG(2)=1 => DERIVATIVES REPOR' // &
                    'TED IN'',/, 2x,''                           HARTEEE/ANGSTROM'',/, 2x,''NFLAG(1)=2 & ' // &
                    'NFLAG(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/' // &
                    'ANGSTROM'',/, 2x,''NFLAG(1)=2 & NFLAG(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''    ' // &
                    '                       eV/ANGSTROM'',/, 2x,''NFLAG(1)=2 & NFLAG(2)=4 => DERIVATIVES ' // &
                    'REPORTED IN'',/, 2x,''                           kcal/mol/ANGSTROM'',/, 2x,''NFLAG(1' // &
                    ')=2 & NFLAG(2)=5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm*' // &
                    '*-1/ANGSTROM'',//)')
            endif
        end subroutine potinfo
        !> @brief Evaluate ANCVRT for the OH2_App_GEN_M3na1988 potential.
        subroutine ancvrt
            implicit none
            integer            :: i, ianum, inc1, inc2, ind, ind2, inu, isave, isum, j, jsave
            character(len=3)   :: periodic_1(7,32)
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            character(len=2)   :: name1(natom)
            character(len=2)   :: name2(natom)
            character(len=1)   :: iblank
            character(len=20)  :: distance
            character(len=20)  :: units
            dimension ianum(7,32)
            dimension isave(natom), jsave(natom)
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
            else if( usricm_nflag(2) .eq. 3 )then
                units = 'EV                  '
            else if( usricm_nflag(2) .eq. 4 )then
                units = 'KCAL PER MOLE       '
            else if( usricm_nflag(2) .eq. 5 )then
                units = 'WAVENUMBERS         '
            else if( usricm_nflag(2) .eq. 6 )then
                units = 'KILOJOULES PER MOLE '
            endif
            utilcm_cnvrtd  = 1.d0
            utilcm_cnvrte  = 1.d0
            utilcm_cnvrtde = 1.d0
            if( usricm_nflag(1) .eq. 2 ) utilcm_cnvrtd = bohr_a
            if( usricm_nflag(2) .eq. 2 )then
                utilcm_cnvrte = utilcm_cnvrte*htomillh
            else if( usricm_nflag(2) .eq. 3 )then
                utilcm_cnvrte = utilcm_cnvrte*htoev
            else if( usricm_nflag(2) .eq. 4 )then
                utilcm_cnvrte = utilcm_cnvrte*htokcal
            else if( usricm_nflag(2) .eq. 5 )then
                utilcm_cnvrte = utilcm_cnvrte*htowave
            else if( usricm_nflag(2) .eq. 6 )then
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
        !> @brief Evaluate CARTOU for the OH2_App_GEN_M3na1988 potential.
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
        !> @brief Evaluate CARTTOR for the OH2_App_GEN_M3na1988 potential.
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
                write(usricm_nflag(18), '(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') infocm_icartr
                stop
            endif
        end subroutine carttor
        !> @brief Evaluate EUNITZERO for the OH2_App_GEN_M3na1988 potential.
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
        !> @brief Evaluate RTOCART for the OH2_App_GEN_M3na1988 potential.
        subroutine rtocart
            implicit none
            double precision   :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes, ygs, yij, za, zb, zc, zd
            integer            :: i, ind, j, j1, j2, k, m1
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            dimension ygs(n3atom), yes(n3atom,isurf), yij(n3atom,jsurf)
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
            else if( infocm_icartr .eq. 2 )then
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
                        else if( j .gt. i )then
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
                        else if( utilcm_ksoffd .gt. 0 )then
                            do j2=utilcm_ksoffd, utilcm_keoffd
                                y                        = pt5cm_deijdr(m1,j2)
                                utilcm_dijcartnu(i,1,j2) = utilcm_dijcartnu(i,1,j2) + termx*y
                                utilcm_dijcartnu(i,2,j2) = utilcm_dijcartnu(i,2,j2) + termy*y
                                utilcm_dijcartnu(i,3,j2) = utilcm_dijcartnu(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            else if( infocm_icartr .eq. 3 )then
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
            else if( infocm_icartr .eq. 4 )then
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
                write(usricm_nflag(18), '(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') infocm_icartr
                stop
            endif
        end subroutine rtocart
        !> @brief Evaluate DEDCOU for the OH2_App_GEN_M3na1988 potential.
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
end module chempotpy_oh2_app_gen_m3na1988_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh2_app_gen_m3na1988_pes, only: evaluate
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
