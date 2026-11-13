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
!                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC
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
!        NASURF - which electronic states availalble
!                 (1,1) = 1 as only gs state available
!        NDER  = 0 => no derivatives should be calculated
!        NDER  = 1 => calculate first derivatives
!        NFLAG - these integer values can be used to flag options
!                within the potential;
!   Potential parameters' default settings
!                  Variable            Default value
!                  NDER                1
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
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the OH2_GEN_DIM1981 potential energy surface for OH2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units, only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
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
        real(wp)              :: infocm_cartnu(natom,3)
        integer               :: infocm_icartr
        integer               :: infocm_indexes(natom)
        integer               :: infocm_irctnt
        integer               :: infocm_mder
        integer               :: infocm_msurf
        integer               :: infocm_natoms
        character             :: infocm_ref(5)
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
        real(wp)              :: pt1cm_degsdr(n3atom)
        real(wp)              :: pt1cm_engygs
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: pt4cm_deesdr(n3atom,isurf)
        real(wp)              :: pt4cm_engyes(isurf)
        real(wp)              :: pt5cm_deijdr(n3atom,jsurf)
        real(wp)              :: pt5cm_engyij(jsurf)
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
        asycm_easyab     = 0.e0_wp; asycm_easyac = 0.e0_wp; asycm_easybc = 0.e0_wp; infocm_cartnu = 0.e0_wp; infocm_icartr = 0
        infocm_indexes   = 0; infocm_irctnt = 0; infocm_mder = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref       = ' '
        nstep_step       = 0.e0_wp
        potcm_a          = 0.e0_wp
        potcm_aab        = 0.e0_wp
        potcm_abc        = 0.e0_wp
        potcm_beta       = 0.e0_wp
        potcm_bsp        = 0.e0_wp
        potcm_cab        = 0.e0_wp
        potcm_cbc        = 0.e0_wp
        potcm_dbcoor     = 0.e0_wp
        potcm_de         = 0.e0_wp
        potcm_gamma      = 0.e0_wp
        potcm_nphi       = 0
        potcm_phi        = 0.e0_wp
        potcm_r1s        = 0.e0_wp
        potcm_r2s        = 0.e0_wp
        potcm_rbb        = 0.e0_wp
        potcm_req        = 0.e0_wp
        potcm_rhh        = 0.e0_wp
        potcm_vsp        = 0.e0_wp
        potcm_xab        = 0.e0_wp
        potcm_xbc        = 0.e0_wp
        pt1cm_degsdr     = 0.e0_wp
        pt1cm_engygs     = 0.e0_wp
        pt1cm_r          = 0.e0_wp
        pt3cm_ezero      = 0.e0_wp
        pt4cm_deesdr     = 0.e0_wp
        pt4cm_engyes     = 0.e0_wp
        pt5cm_deijdr     = 0.e0_wp
        pt5cm_engyij     = 0.e0_wp
        usricm_anuzero   = 0.e0_wp
        usricm_cart      = 0.e0_wp
        usricm_nasurf    = 0
        usricm_nder      = 0
        usricm_nflag     = 0
        usricm_nulbl     = 0
        usrocm_descart   = 0.e0_wp
        usrocm_dgscart   = 0.e0_wp
        usrocm_dijcart   = 0.e0_wp
        usrocm_pengyes   = 0.e0_wp
        usrocm_pengygs   = 0.e0_wp
        usrocm_pengyij   = 0.e0_wp
        utilcm_cnvrtd    = 0.e0_wp
        utilcm_cnvrtde   = 0.e0_wp
        utilcm_cnvrte    = 0.e0_wp
        utilcm_descartnu = 0.e0_wp
        utilcm_dgscartnu = 0.e0_wp
        utilcm_dijcartnu = 0.e0_wp
        utilcm_ireorder  = 0
        utilcm_kediag    = 0
        utilcm_keoffd    = 0
        utilcm_ksdiag    = 0
        utilcm_ksoffd    = 0
        usricm_nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0], shape(usricm_nasurf))
        usricm_nder    = 0
        usricm_nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero = 0.0e0_wp
        infocm_icartr  = 3
        infocm_msurf   = 0
        infocm_mder    = 1
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 3
        nstep_step     = 0.000000008e0_wp
        potcm_nphi     = 17
        potcm_r1s      = 0.2212079e+01_wp
        potcm_r2s      = 0.2215309e+01_wp
        potcm_vsp      = 0.1094890e+03_wp
        potcm_xbc      = [0.1094890e+03_wp, 0.1473409e+01_wp, 0.1941998e+01_wp]
        potcm_xab      = [0.1056400e+03_wp, 0.1241479e+01_wp, 0.2304798e+01_wp]
        potcm_cbc      = [-0.4815674e-01_wp, 0.3814697e-05_wp, 0.4415512e-03_wp]
        potcm_abc      = [0.7873216e+01_wp, 0.7873216e+01_wp, 0.7873216e+01_wp]
        potcm_cab      = [-0.1482804e+01_wp, -0.2084732e-02_wp, -0.3727913e-02_wp]
        potcm_aab      = [0.1932070e+01_wp, 0.1932070e+01_wp, 0.1932070e+01_wp]
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
        potcm_de      = 106.56e0_wp
        potcm_beta    = 2.07942e0_wp
        potcm_req     = [0.96966e0_wp, 0.74144e0_wp, 0.96966e0_wp]
        potcm_a       = 0.26e0_wp
        potcm_gamma   = 0.4131e0_wp
        potcm_rbb     = 0.74144e0_wp
        v             = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        usricm_nder   = igrad
        usricm_nasurf = 1
        usricm_cart   = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                usricm_cart(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = usrocm_dgscart(iatom,idir)*hartree_per_bohr_to_ev_per_angstrom
                enddo
            enddo
        elseif( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif

        contains

        !> @brief Evaluate prepare_local_context for the OH2_GEN_DIM1981 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i, it
            real(wp)            :: nab, nab1, nbc, nbc1, junk(4)
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
            real(wp), parameter :: bohr = .52917706e0_wp
            ! Conversion constant from Hartree to kcal (kcal/Hartree)
            real(wp), parameter :: hrtree = 627.5095e0_wp
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !      WRITE(NFLAG(18), 600)
            call prepare_local_context_2
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !      BENDING POTENTIAL PARAMETERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !      WRITE(NFLAG(18),3)
            !   Echo the potential parameters to the file linked to UNIT NFLAG(18)
            !      WRITE (NFLAG(18),25) (REQ(IT), IT = 1,3)
            !      WRITE(NFLAG(18),10) DE
            !      WRITE(NFLAG(18),21) BETA
            !      WRITE (NFLAG(18),30) A
            !      WRITE (NFLAG(18),31) GAMMA
            !      WRITE(NFLAG(18),35) RBB
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
            infocm_ref(1)     = 'No Reference'
            infocm_indexes(1) = 8
            infocm_indexes(2) = 1
            infocm_indexes(3) = 1
            infocm_irctnt     = 2
            !      CALL POTINFO
            call ancvrt
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
            real(wp)            :: bcoor, c, dac, ractmp, rob, stuff, vabbc, vac
            integer             :: it
            real(wp)            :: nab, nab1, nbc, nbc1, junk(4)
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
            real(wp), parameter :: bohr = .52917706e0_wp
            ! Conversion constant from Hartree to kcal (kcal/Hartree)
            real(wp), parameter :: hrtree = 627.5095e0_wp
            call cartou
            call carttor
            !   Check the values of NASURF and NDER for validity.
            !      IF (NASURF(1,1) .EQ. 0) THEN
            !         WRITE(NFLAG(18), 900) NASURF(1,1)
            !         STOP
            !      ENDIF
            !         IF (NDER .GT. 1) THEN
            !             WRITE (NFLAG(18), 910) NDER
            !             STOP 'POT 2'
            !         ENDIF
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
            if( usricm_nder .eq. 1 )then
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
            if( usricm_nder .eq. 1 )then
                pt1cm_degsdr(1) = pt1cm_degsdr(1) + pt1cm_degsdr(3)
                pt1cm_degsdr(2) = pt1cm_degsdr(2) + pt1cm_degsdr(3)
                pt1cm_degsdr(3) = dac
            endif
            pt1cm_r(3)   = ractmp
            pt1cm_engygs = pt1cm_engygs + bcoor
            if( usricm_nder .eq. 1 )then
                do it = 1, 2
                    pt1cm_degsdr(it) = potcm_dbcoor(it) + pt1cm_degsdr(it)
                enddo
            endif
            call eunitzero
            if( usricm_nder .ne. 0 )then
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
            real(wp)            :: csp(48,3), scrap(48), v(4)
            real(wp)            :: nab, nab1, nbc, nbc1, junk(4)
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
            real(wp), parameter :: bohr = .52917706e0_wp
            ! Conversion constant from Hartree to kcal (kcal/Hartree)
            real(wp), parameter :: hrtree = 627.5095e0_wp
            real(wp), parameter :: detrad = .0174532e0_wp
            character(len=80)   :: title
            title = ' RMOFIT to DIM surface with BEBO continuations '
            !      WRITE(NFLAG(18),845)TITLE
            !      WRITE(NFLAG(18),125)STEP
            !      WRITE(NFLAG(18),111)XAB(1), XBC(1), XAB(1),
            !     *               XAB(2), XBC(2), XAB(2),
            !     *               XAB(3), XBC(3), XAB(3)
            !      WRITE(NFLAG(18),105)R1S,R2S,VSP
            !      WRITE(NFLAG(18),113)(PHI(I),(BSP(I,J),J = 1,3),I = 1,NPHI)
            !      WRITE(NFLAG(18),117)(CBC(I),ABC(I),I = 1,3),
            !     +                    (CAB(I),AAB(I),I = 1,3)
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
            integer             :: i, ic, j, k
            real(wp)            :: nab, nab1, nbc, nbc1, junk(4), bohr
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            real(wp), parameter :: hrtree = 627.5095e0_wp
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
                    pt1cm_engygs = potcm_xbc(1)/627.5095e0_wp
                    if( usricm_nder == 1 )then
                        pt1cm_degsdr(1) = 1.0e-9_wp
                        pt1cm_degsdr(2) = 1.0e-9_wp
                        pt1cm_degsdr(3) = 0.0e0_wp
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
                x     = dis*((1.0e0_wp-exp(bet*(rr-req2)))**2-1.0e0_wp)+pt3cm_ezero(1)-usricm_anuzero*627.5095e0_wp
                v(ic) = x
            enddo
            v(1:3)       = v(1:3)/627.5095e0_wp
            pt1cm_engygs = v(1)
            if( usricm_nder == 1 )then
                pt1cm_degsdr(1) = (v(2)-v(1))/nstep_step
                pt1cm_degsdr(2) = (v(3)-v(1))/nstep_step
                pt1cm_degsdr(3) = 0.0e0_wp
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
                    if( d(i-1) == 0.0e0_wp ) error stop 'SPLINE: singular knot system'
                    pivot = 1.0e0_wp/d(i-1)
                    if( i < np2 )then
                        supd = x(i-1)-x(i-2)
                        if( supd < 0.0e0_wp ) error stop 'SPLINE: knots must be ordered'
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
                    if( d(j+1) == 0.0e0_wp ) error stop 'SPLINE: singular back substitution'
                    c(j+1) = c(j+1)/d(j+1)
                    c(j)   = c(j)-c(j+1)*v
                enddo
                if( d(2) == 0.0e0_wp ) error stop 'SPLINE: singular leading coefficient'
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
                error stop 'SPLINE: invalid operation selector'
            end select
        end function spline

        !> @brief Evaluate POTINFO for the OH2_GEN_DIM1981 potential.
        subroutine potinfo
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i, j, k, kmax
            integer, parameter  :: n3atom = 75
            integer, parameter  :: natom  = 25
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
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

        !> @brief Evaluate ANCVRT for the OH2_GEN_DIM1981 potential.
        subroutine ancvrt
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter  :: n3atom     = 75
            integer, parameter  :: natom      = 25
            integer, parameter  :: isurf      = 5
            integer, parameter  :: jsurf      = isurf*(isurf+1)/2
            real(wp), parameter :: pi         = 3.141592653589793e0_wp
            real(wp), parameter :: clight     = 2.99792458e08_wp
            real(wp), parameter :: cmu_0      = 4.0e0_wp*pi*1.0e-07_wp
            real(wp), parameter :: cepsilon_0 = 1.0e0_wp/(cmu_0*clight**2)
            real(wp), parameter :: ce         = 1.602176462e-19_wp
            real(wp), parameter :: cplanck    = 6.62606876e-34_wp
            real(wp), parameter :: cm_e       = 9.10938188e-31_wp
            real(wp), parameter :: cang       = 1.0e-10_wp
            real(wp), parameter :: cavogadro  = 6.02214199e23_wp
            real(wp), parameter :: ckcal      = 4.184e10_wp
            real(wp), parameter :: htomillh   = 1000.e0_wp
            real(wp), parameter :: htoev      = 27.2113834e0_wp
            real(wp), parameter :: htokcal    = 627.509470e0_wp
            real(wp), parameter :: htowave    = 219474.631e0_wp
            real(wp), parameter :: htokj      = 2625.49962e0_wp
            real(wp), parameter :: bohr_a     = .5291772083e0_wp
            integer             :: i, ianum(7,32), inc1, inc2, ind, ind2, inu, isave(natom), isum, j, jsave(natom)
            character(len=3)    :: periodic_1(7,32)
            character(len=2)    :: name1(natom)
            character(len=2)    :: name2(natom)
            character(len=1)    :: iblank
            character(len=20)   :: distance
            character(len=20)   :: units
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
            do inu=1, 25
                isum = isum + usricm_nulbl(inu)
            enddo
            utilcm_ireorder = 0
            if( isum .ne. 0 ) utilcm_ireorder = 1
        end subroutine ancvrt

        !> @brief Evaluate CARTOU for the OH2_GEN_DIM1981 potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            use chempotpy_kinds,              only: wp
            implicit none
            integer            :: i, j
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            call convert_cartesian_input(usricm_cart, usricm_nulbl, infocm_natoms, utilcm_cnvrtd, utilcm_ireorder == 1, infocm_cartnu)
        end subroutine cartou

        !> @brief Evaluate CARTTOR for the OH2_GEN_DIM1981 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            use chempotpy_kinds,              only: wp
            implicit none
            integer            :: coordinate_status
            real(wp)           :: cmm, cphi, flm, hhd, hym, phi, pi, q, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
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

        !> @brief Evaluate EUNITZERO for the OH2_GEN_DIM1981 potential.
        subroutine eunitzero
            use chempotpy_kinds, only: wp
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

        !> @brief Evaluate RTOCART for the OH2_GEN_DIM1981 potential.
        subroutine rtocart
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
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
                write(usricm_nflag(18), '(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') infocm_icartr
                stop
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the OH2_GEN_DIM1981 potential.
        subroutine dedcou
            use chempotpy_kinds, only: wp
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
