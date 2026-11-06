!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_oh2_gen_dim1981_pes
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
        double precision              :: nstep_step
        double precision              :: potcm_a
        double precision              :: potcm_aab(3)
        double precision              :: potcm_abc(3)
        double precision              :: potcm_beta
        double precision              :: potcm_bsp(40,3)
        double precision              :: potcm_cab(3)
        double precision              :: potcm_cbc(3)
        double precision              :: potcm_dbcoor(2)
        double precision              :: potcm_de
        double precision              :: potcm_gamma
        integer                       :: potcm_nphi
        double precision              :: potcm_phi(40)
        double precision              :: potcm_r1s
        double precision              :: potcm_r2s
        double precision              :: potcm_rbb
        double precision              :: potcm_req(3)
        double precision              :: potcm_rhh
        double precision              :: potcm_vsp
        double precision              :: potcm_xab(3)
        double precision              :: potcm_xbc(3)
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
        asycm_easyab     = 0.d0; asycm_easyac = 0.d0; asycm_easybc = 0.d0; infocm_cartnu = 0.d0; infocm_icartr = 0
        infocm_indexes   = 0; infocm_irctnt = 0; infocm_mder = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref       = ' '
        nstep_step       = 0.d0
        potcm_a          = 0.d0
        potcm_aab        = 0.d0
        potcm_abc        = 0.d0
        potcm_beta       = 0.d0
        potcm_bsp        = 0.d0
        potcm_cab        = 0.d0
        potcm_cbc        = 0.d0
        potcm_dbcoor     = 0.d0
        potcm_de         = 0.d0
        potcm_gamma      = 0.d0
        potcm_nphi       = 0
        potcm_phi        = 0.d0
        potcm_r1s        = 0.d0
        potcm_r2s        = 0.d0
        potcm_rbb        = 0.d0
        potcm_req        = 0.d0
        potcm_rhh        = 0.d0
        potcm_vsp        = 0.d0
        potcm_xab        = 0.d0
        potcm_xbc        = 0.d0
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
        nstep_step     = 0.000000008d0
        potcm_nphi     = 17
        potcm_r1s      = 0.2212079d+01
        potcm_r2s      = 0.2215309d+01
        potcm_vsp      = 0.1094890d+03
        potcm_xbc      = [0.1094890d+03, 0.1473409d+01, 0.1941998d+01]
        potcm_xab      = [0.1056400d+03, 0.1241479d+01, 0.2304798d+01]
        potcm_cbc      = [-0.4815674d-01, 0.3814697d-05, 0.4415512d-03]
        potcm_abc      = [0.7873216d+01, 0.7873216d+01, 0.7873216d+01]
        potcm_cab      = [-0.1482804d+01, -0.2084732d-02, -0.3727913d-02]
        potcm_aab      = [0.1932070d+01, 0.1932070d+01, 0.1932070d+01]
        potcm_phi(1:17) = [0.0000000d+00, 0.2000000d+01, 0.4000000d+01, 0.6000000d+01, 0.1800000d+02, 0.2400000d+02, &
            & 0.3000000d+02, 0.3400000d+02, 0.3800000d+02, 0.4400000d+02, 0.4800000d+02, 0.5400000d+02, 0.6000000d+02, &
            & 0.7600000d+02, 0.8600000d+02, 0.8800000d+02, 0.9000000d+02]
        potcm_bsp(1:17,1) = [0.1094409d+03, 0.1094212d+03, 0.1094004d+03, 0.1093775d+03, 0.1088412d+03, 0.1077417d+03, &
            & 0.1055117d+03, 0.1030257d+03, 0.9966266d+02, 0.9630693d+02, 0.9637480d+02, 0.9811765d+02, 0.1000138d+03, &
            & 0.1029822d+03, 0.1038929d+03, 0.1040309d+03, 0.1041572d+03]
        potcm_bsp(1:17,2) = [0.1473413d+01, 0.1474257d+01, 0.1476847d+01, 0.1481213d+01, 0.1545280d+01, 0.1599182d+01, &
            & 0.1655956d+01, 0.1681611d+01, 0.1680401d+01, 0.1626349d+01, 0.1572269d+01, 0.1486160d+01, 0.1408802d+01, &
            & 0.1273705d+01, 0.1241724d+01, 0.1239825d+01, 0.1239394d+01]
        potcm_bsp(1:17,3) = [0.1942440d+01, 0.1940545d+01, 0.1936212d+01, 0.1929558d+01, 0.1825750d+01, 0.1719964d+01, &
            & 0.1603891d+01, 0.1538857d+01, 0.1490752d+01, 0.1525683d+01, 0.1620584d+01, 0.1788983d+01, 0.1942842d+01, &
            & 0.2222638d+01, 0.2294156d+01, 0.2299114d+01, 0.2301070d+01]
        potcm_de      = 106.56d0
        potcm_beta    = 2.07942d0
        potcm_req     = [0.96966d0, 0.74144d0, 0.96966d0]
        potcm_a       = 0.26d0
        potcm_gamma   = 0.4131d0
        potcm_rbb     = 0.74144d0
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
        !> @brief Evaluate prepare_local_context for the OH2_GEN_DIM1981 potential.
        subroutine prepare_local_context
            !
            !
            !   System:          OH2
            !   Functional form: Anti-Morse Bends plus Bowman's collinear spline
            !   Common name:     OH2DIM
            !   Interface:       potlib2001
            !   Reference:
            !
            !   Number of bodies:              3
            !   Number of derivatives:         1
            !   Number of electronic surfaces: 1
            !
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
            !
            !   This potential is written such that:
            !                  R(1) = R(O-H)
            !                  R(2) = R(H-H)
            !                  R(3) = R(H-O)
            !   The zero of energy is defined at O "infinitely" far from the H2 diat
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
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !          GENERAL INFORMATION
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !
            !   POTA.FOR HAS ANTI-MORSE BEND PLUS BOWMAN'S COLLINEAR SPLINE ROUTINE.
            !   ANTI-MORSE BENDING POTENTIAL IS ADDED TO COLLINEAR POTENTIAL.
            !   THIS SUBROUTINE CAN BE USED TO ADD A BENDING CORRECTION TO
            !   ANY COLLINEAR POTENTIAL SUBROUTINE.
            !
            !CCCCCCCCCCCCCCCCCCCCCCCCCC
            !          DEFINE VARIABLES
            !CCCCCCCCCCCCCCCCCCCCCCCCCC
            !
            implicit none
            integer :: i, it
            !
            double precision :: nab, nab1, nbc, nbc1, junk
            dimension junk(4)
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
            ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
            !
            double precision, parameter :: bohr = .52917706d0
            !
            ! Conversion constant from Hartree to kcal (kcal/Hartree)
            !
            double precision, parameter :: hrtree = 627.5095d0
            !
            !
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !
            !      WRITE(NFLAG(18), 600)
            !
            call prepare_local_context_2
            !
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !      BENDING POTENTIAL PARAMETERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !
            !      WRITE(NFLAG(18),3)
            !
            !   Echo the potential parameters to the file linked to UNIT NFLAG(18)
            !
            !      WRITE (NFLAG(18),25) (REQ(IT), IT = 1,3)
            !      WRITE(NFLAG(18),10) DE
            !      WRITE(NFLAG(18),21) BETA
            !      WRITE (NFLAG(18),30) A
            !      WRITE (NFLAG(18),31) GAMMA
            !      WRITE(NFLAG(18),35) RBB
            !
            !
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !       CONVERT TO ATOMIC UNITS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !
            do it = 1, 3
                potcm_req(it) = potcm_req(it)/bohr
            50 enddo
            potcm_beta = potcm_beta * bohr
            potcm_de   = potcm_de/hrtree
            potcm_a    = potcm_a/bohr
            potcm_rhh  = 1.4007977d0
            potcm_rbb  = potcm_rbb /bohr
            !
            !CCCCCCCCCCC
            !
            !
            pt3cm_ezero(1) = potcm_vsp
            !
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !
            infocm_ref(1) = 'No Reference'
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
        !> @brief Evaluate POT for the OH2_GEN_DIM1981 potential.
        subroutine pot
            !
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the OH diatomic, with the
            !   OH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
            !   respectively.
            !
            !   This potential is written such that:
            !                  R(1) = R(O-H)
            !                  R(2) = R(H-H)
            !                  R(3) = R(H-O)
            !   The zero of energy is defined at O "infinitely" far from the H2 diat
            !
            !CCCCCCCCCCCCCCC
            !CCCCCCCCCCCCCCC
            !
            implicit none
            double precision :: bcoor, c, dac, ractmp, rob, stuff, vabbc, vac
            integer          :: it
            !
            double precision :: nab, nab1, nbc, nbc1, junk
            dimension junk(4)
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
            ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
            !
            double precision, parameter :: bohr = .52917706d0
            !
            ! Conversion constant from Hartree to kcal (kcal/Hartree)
            !
            double precision, parameter :: hrtree = 627.5095d0
            !
            !
            call cartou
            call carttor
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
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        FIRST CALCULATE SOME USEFUL NUMBERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !
            !        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
            !        AB AND BC RESPECTIVELY. (EQUATION 2.2)
            !
            nab = exp((potcm_req(1) - pt1cm_r(1))/potcm_a)
            !
            if( nab .lt. 1.d-15 ) nab = 1.d-15
            !
            nbc = exp((potcm_req(2) - pt1cm_r(2))/potcm_a)
            !
            if( nbc .lt. 1.d-15 ) nbc = 1.d-15
            !
            !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
            !       SEE EQUATION 2.7A, 2.7B
            !
            c = 1.0d0-nab/nbc
            if( abs(c) <= 1.0d-14 )then
                nab1 = 0.5d0
                nbc1 = 0.5d0
            else
                c    = c/nab
                nab1 = (2.0d0+c-sqrt(4.0d0+c*c))/(2.0d0*c)
                nbc1 = 1.0d0-nab1
            endif
            !
            !          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
            !
            !          FIRST TRAP OUT ANY ZERO ARGUEMENTS
            !
            if( nab1*nbc1 > 0.0d0 )then
                stuff = potcm_a * log(nab1*nbc1)
                rob   = (potcm_rbb - stuff)/(potcm_rhh - stuff)
            else
                rob = 1.0d0
            endif
            !
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        CALCULATE BENDING CORRECTION
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !
            !        CALCULATE V(R(3))
            !
            junk(3) = exp(-potcm_beta*(pt1cm_r(3)-potcm_req(3))/rob)
            vac     = potcm_gamma *1.0d-14 * potcm_de * junk(3)* (1.d0 + 0.5d0*junk(3))
            !
            !        CALCULATE V(R(1)+R(2))
            !
            junk(4) = exp(potcm_beta*(potcm_req(3)-pt1cm_r(2)-pt1cm_r(1))/rob)
            vabbc   = potcm_gamma *1.0d-14 * potcm_de * junk(4) * (1.d0 + 0.5d0*junk(4))
            !
            !        HERE IS THE BENDING CORRECTION
            !
            bcoor = (vac - vabbc)*1.0d14
            !
            !          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
            !          OF BCOOR WITH RESPECT TO R(1), R(2), AND R(3).
            !
            !        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
            !
            if( usricm_nder .eq. 1 )then
                dac             = -potcm_gamma*potcm_de*potcm_beta*junk(3)*(1.d0+junk(3))/rob
                potcm_dbcoor(1) = potcm_gamma * potcm_de * potcm_beta * junk(4) * (1.d0 + junk(4))/rob
                potcm_dbcoor(2) = potcm_dbcoor(1)
            endif
            !
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        NOW CALCULATE THE COLLINEAR ENERGY
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !
            !          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR
            !          GEOMETRY AND CALL POT2
            !
            ractmp     = pt1cm_r(3)
            pt1cm_r(3) = pt1cm_r(2) + pt1cm_r(1)
            !
            call pot2
            !
            if( usricm_nder .eq. 1 )then
                pt1cm_degsdr(1) = pt1cm_degsdr(1) + pt1cm_degsdr(3)
                pt1cm_degsdr(2) = pt1cm_degsdr(2) + pt1cm_degsdr(3)
                pt1cm_degsdr(3) = dac
            endif
            pt1cm_r(3) = ractmp
            !
            pt1cm_engygs = pt1cm_engygs + bcoor
            !
            if( usricm_nder .eq. 1 )then
                do it = 1, 2
                    pt1cm_degsdr(it) = potcm_dbcoor(it) + pt1cm_degsdr(it)
                300 enddo
            endif
            !
            !
            call eunitzero
            if( usricm_nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
            !
        end subroutine pot
        !> @brief Evaluate prepare_local_context_2 for the OH2_GEN_DIM1981 potential.
        subroutine prepare_local_context_2
            !
            !   POTENTIAL FUNCTION FOR RATE1D SERIES
            !   RMO-SPLINE TYPE FIT WITH BEBO TYPE CONTINUATION INTO THE ASYMBTOTIC
            !   POTENTIAL ROUTINE INTERACTS WITH REST OF PROGRAM WITH UNITS OF EV/BO
            !   POTENTIAL ROUTINE ASSUMES SPLINE INFO IS IN UNITS OF KCAL/ANGSTROM
            !
            implicit none
            double precision :: csp, scrap, v
            !
            double precision :: nab, nab1, nbc, nbc1, junk
            dimension junk(4)
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
            ! Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
            !
            double precision, parameter :: bohr = .52917706d0
            !
            ! Conversion constant from Hartree to kcal (kcal/Hartree)
            !
            double precision, parameter :: hrtree = 627.5095d0
            !
            double precision, parameter :: detrad = .0174532d0
            !
            !
            dimension v(4), csp(48,3), scrap(48)
            !
            character(len=80) :: title
            !
            title = ' RMOFIT to DIM surface with BEBO continuations '
            !
            !      WRITE(NFLAG(18),845)TITLE
            !      WRITE(NFLAG(18),125)STEP
            !      WRITE(NFLAG(18),111)XAB(1), XBC(1), XAB(1),
            !     *               XAB(2), XBC(2), XAB(2),
            !     *               XAB(3), XBC(3), XAB(3)
            !      WRITE(NFLAG(18),105)R1S,R2S,VSP
            !      WRITE(NFLAG(18),113)(PHI(I),(BSP(I,J),J = 1,3),I = 1,NPHI)
            !      WRITE(NFLAG(18),117)(CBC(I),ABC(I),I = 1,3),
            !     +                    (CAB(I),AAB(I),I = 1,3)
            !
            !   Initialize the energy in the asymptotic valleys
            !
            asycm_easyab = potcm_xab(1)/hrtree
            asycm_easybc = potcm_xbc(1)/hrtree
            asycm_easyac = asycm_easyab
            !
            !
        end subroutine prepare_local_context_2
        !> @brief Evaluate POT2 for the OH2_GEN_DIM1981 potential.
        subroutine pot2
            !
            ! COCK SPLINES
            !
            implicit none
            double precision :: angle, bet, betpr, betrg, csp, dbc, depr, derg, dis, dr1, dr2, ppp, r1, r2, r3, repr, req2, rerg, rr, scrap, sep, v, x
            integer          :: i, ic, j, k
            !
            double precision :: nab, nab1, nbc, nbc1, junk, bohr
            dimension junk(4)
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
            double precision, parameter :: hrtree = 627.5095d0
            double precision, parameter :: detrad = .0174532d0
            !
            !
            dimension v(4), csp(48,3), scrap(48)
            !
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
                r1    = r1*0.529177d0
                r2    = r2*0.529177d0
                r3    = r3*0.529177d0
                dr1   = potcm_r1s-r1
                dr2   = potcm_r2s-r2
                angle = atan(dr1/dr2)/detrad
                do k = 1, 3
                    ppp = spline(1, potcm_nphi, potcm_phi, potcm_bsp(1,k), csp(1,k), scrap, angle)
                enddo
                dbc   = potcm_xbc(1)
                derg  = potcm_xbc(1)/23.061d0
                rerg  = (potcm_r2s-potcm_xbc(2))/.529177d0
                betrg = potcm_xbc(3)*.529177d0
                depr  = potcm_xab(1)/23.061d0
                repr  = (potcm_r1s-potcm_xab(2))/.529177d0
                betpr = potcm_xab(3)*.529177d0
                if( r1 > potcm_r1s .and. r2 > potcm_r2s )then
                    pt1cm_engygs = potcm_xbc(1)/627.5095d0
                    if( usricm_nder == 1 )then
                        pt1cm_degsdr(1) = 1.0d-9
                        pt1cm_degsdr(2) = 1.0d-9
                        pt1cm_degsdr(3) = 0.0d0
                    endif
                    return
                else if( r1 > potcm_r1s )then
                    rr   = potcm_r2s-r2
                    sep  = r1-potcm_r1s
                    dis  = potcm_xbc(1)+potcm_cbc(1)*exp(-min(potcm_abc(1)*sep, 25.0d0))
                    req2 = potcm_xbc(2)+potcm_cbc(2)*exp(-min(potcm_abc(2)*sep, 25.0d0))
                    bet  = potcm_xbc(3)+potcm_cbc(3)*exp(-min(potcm_abc(3)*sep, 25.0d0))
                else if( r2 > potcm_r2s )then
                    rr   = potcm_r1s-r1
                    sep  = r2-potcm_r2s
                    dis  = potcm_xab(1)+potcm_cab(1)*exp(-min(potcm_aab(1)*sep, 25.0d0))
                    req2 = potcm_xab(2)+potcm_cab(2)*exp(-min(potcm_aab(2)*sep, 25.0d0))
                    bet  = potcm_xab(3)+potcm_cab(3)*exp(-min(potcm_aab(3)*sep, 25.0d0))
                else
                    dr1   = potcm_r1s-r1
                    dr2   = potcm_r2s-r2
                    angle = atan(dr1/dr2)/detrad
                    rr    = sqrt(dr1*dr1+dr2*dr2)
                    dis   = spline(2, potcm_nphi, potcm_phi, potcm_bsp(1,1), csp(1,1), scrap, angle)
                    req2  = spline(2, potcm_nphi, potcm_phi, potcm_bsp(1,2), csp(1,2), scrap, angle)
                    bet   = spline(2, potcm_nphi, potcm_phi, potcm_bsp(1,3), csp(1,3), scrap, angle)
                endif
                x     = dis*((1.0d0-exp(bet*(rr-req2)))**2-1.0d0)+pt3cm_ezero(1)-usricm_anuzero*627.5095d0
                v(ic) = x
            enddo
            v(1:3)       = v(1:3)/627.5095d0
            pt1cm_engygs = v(1)
            if( usricm_nder == 1 )then
                pt1cm_degsdr(1) = (v(2)-v(1))/nstep_step
                pt1cm_degsdr(2) = (v(3)-v(1))/nstep_step
                pt1cm_degsdr(3) = 0.0d0
            endif
        end subroutine pot2
        !> @brief Fit or evaluate the legacy cubic-spline representation.
        double precision function spline( isw, nn, x, y, c, d, xpt )
            implicit none
            integer, intent(in)             :: isw, nn
            double precision, intent(in)    :: x(nn), y(nn), xpt
            double precision, intent(inout) :: c(nn+2), d(nn+2)
            integer                         :: i, j, n, np1, np2
            double precision                :: cfact, dfact, pivot, supd, v, v2, z
            double precision, save          :: deriv = 0.0d0
            n   = nn
            np1 = n+1
            np2 = n+2
            z   = xpt
            select case( isw )
            case( 1 )
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
                    if( d(i-1) == 0.0d0 ) error stop 'SPLINE: singular knot system'
                    pivot = 1.0d0/d(i-1)
                    if( i < np2 )then
                        supd = x(i-1)-x(i-2)
                        if( supd < 0.0d0 ) error stop 'SPLINE: knots must be ordered'
                        supd = supd**3
                    else
                        supd = 1.0d0
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
                        d(np1) = 1.0d0-d(np1)*dfact
                    endif
                    c(np2) = c(np2)-d(np2)*cfact
                    d(np2) = x(i-2)-d(np2)*dfact
                enddo
                do i = 1, n
                    j = np2-i
                    if( j == np1 )then
                        v = 1.0d0
                    else
                        v = (x(j)-x(j-1))**3
                    endif
                    if( d(j+1) == 0.0d0 ) error stop 'SPLINE: singular back substitution'
                    c(j+1) = c(j+1)/d(j+1)
                    c(j)   = c(j)-c(j+1)*v
                enddo
                if( d(2) == 0.0d0 ) error stop 'SPLINE: singular leading coefficient'
                c(2)   = c(2)/d(2)
                spline = 0.0d0
            case( 2 )
                spline = c(1)+c(2)*(z-x(1))
                do i = 1, n
                    v = z-x(i)
                    if( v <= 0.0d0 ) exit
                    spline = spline+c(i+2)*v**3
                enddo
            case( 3 )
                spline = c(1)+c(2)*(z-x(1))
                deriv  = c(2)/3.0d0
                do i = 1, n
                    v = z-x(i)
                    if( v <= 0.0d0 ) exit
                    v2     = v*v
                    spline = spline+c(i+2)*v2*v
                    deriv  = deriv+c(i+2)*v2
                enddo
            case( 4 )
                spline = 3.0d0*deriv
            case default
                error stop 'SPLINE: invalid operation selector'
            end select
        end function spline
        !> @brief Evaluate POTINFO for the OH2_GEN_DIM1981 potential.
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
        !> @brief Evaluate ANCVRT for the OH2_GEN_DIM1981 potential.
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
        !> @brief Evaluate CARTOU for the OH2_GEN_DIM1981 potential.
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
        !> @brief Evaluate CARTTOR for the OH2_GEN_DIM1981 potential.
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
        !> @brief Evaluate EUNITZERO for the OH2_GEN_DIM1981 potential.
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
        !> @brief Evaluate RTOCART for the OH2_GEN_DIM1981 potential.
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
                        20 continue
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
        !> @brief Evaluate DEDCOU for the OH2_GEN_DIM1981 potential.
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
end module chempotpy_oh2_gen_dim1981_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh2_gen_dim1981_pes, only: evaluate
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
