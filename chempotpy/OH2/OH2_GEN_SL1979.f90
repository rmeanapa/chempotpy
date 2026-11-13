!   System:          OH2
!   Functional form:
!   Common name:     OH2SL
!   Reference:
!   Interface:                     potlib2001
!   Number of bodies:              3
!   Number of derivatives:         1
!   Number of electronic surfaces: 1
!   Notes: This is a collinear potential energy surface, i.e., the poten
!          energy is not defined for non-linear geometries.
!          The derivatives of the potential energy with respect to the
!          coordinates are determined numerically using a step size STEP
!          if R(1) and R(2) are less than 50 bohr.
!          The variable STEP is initialized in the block data subprogram
!          PTPACM.
!   prepare_local_context must be called once before any calls to POT.
!   The potential parameters are assigned in the Block Data Subprogram P
!   Coordinates, potential energy, and derivatives are passed
!   through the common block PT1CM:
!                  COMMON /PT1CM/ R(3), ENGYGS, DEGSDR(3)
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
!   The flags that indicate which calculations should be carried out in
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

!> @brief Implement OH2_GEN_SL1979, a potential energy surface for OH2.
module chempotpy_oh2_gen_sl1979_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the OH2_GEN_SL1979 potential energy surface for OH2.
    subroutine evaluate( x, igrad, p, g, d )
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
        real(wp)              :: ndercm_step
        real(wp)              :: precm_a
        real(wp)              :: precm_al(3)
        real(wp)              :: precm_alf
        real(wp)              :: precm_b
        real(wp)              :: precm_c(50)
        real(wp)              :: precm_chi(3)
        real(wp)              :: precm_gam
        real(wp)              :: precm_s(3)
        real(wp)              :: precm_v(4)
        real(wp)              :: precm_vm(3)
        real(wp)              :: precm_x(3)
        real(wp)              :: pt1cm_degsdr(n3atom)
        real(wp)              :: pt1cm_engygs
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: pt4cm_deesdr(n3atom,isurf)
        real(wp)              :: pt4cm_engyes(isurf)
        real(wp)              :: pt5cm_deijdr(n3atom,jsurf)
        real(wp)              :: pt5cm_engyij(jsurf)
        real(wp)              :: satocm_bet(3)
        real(wp)              :: satocm_de(3)
        real(wp)              :: satocm_xre(3)
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
        ndercm_step      = 0.e0_wp
        precm_a          = 0.e0_wp
        precm_al         = 0.e0_wp
        precm_alf        = 0.e0_wp
        precm_b          = 0.e0_wp
        precm_c          = 0.e0_wp
        precm_chi        = 0.e0_wp
        precm_gam        = 0.e0_wp
        precm_s          = 0.e0_wp
        precm_v          = 0.e0_wp
        precm_vm         = 0.e0_wp
        precm_x          = 0.e0_wp
        pt1cm_degsdr     = 0.e0_wp
        pt1cm_engygs     = 0.e0_wp
        pt1cm_r          = 0.e0_wp
        pt3cm_ezero      = 0.e0_wp
        pt4cm_deesdr     = 0.e0_wp
        pt4cm_engyes     = 0.e0_wp
        pt5cm_deijdr     = 0.e0_wp
        pt5cm_engyij     = 0.e0_wp
        satocm_bet       = 0.e0_wp
        satocm_de        = 0.e0_wp
        satocm_xre       = 0.e0_wp
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
        satocm_de      = [4.035e0_wp, 4.1662e0_wp, 4.035e0_wp]
        satocm_xre     = [0.9923e0_wp, 0.7590e0_wp, 0.9923e0_wp]
        satocm_bet     = [2.4113e0_wp, 2.0682e0_wp, 2.4113e0_wp]
        ndercm_step    = 5.0e-7_wp
        precm_a        = 0.157002e-03_wp
        precm_b        = -0.162301e-01_wp
        precm_alf      = 0.9e0_wp
        precm_gam      = 0.045e0_wp
        precm_al       = [1.0e0_wp, 0.8e0_wp, 1.0e0_wp]
        precm_c = [0.354195e+01_wp, 0.302780e+01_wp, 0.577312e+00_wp, 0.539877e+01_wp, 0.233027e+01_wp, &
            & 0.157276e+02_wp, 0.701851e+01_wp, - 0.478535e+02_wp, -0.106915e+01_wp, 0.268137e+02_wp, - &
            & 0.104806e+02_wp, -0.182660e+02_wp, 0.833941e+01_wp, 0.872497e+02_wp, 0.917870e+01_wp, -0.697766e+02_wp, &
            & 0.434862e+02_wp, -0.156633e+02_wp, 0.380162e+02_wp, -0.796542e+01_wp, 0.602854e+02_wp, 0.500800e+02_wp, &
            & -0.545761e+02_wp, 0.410210e+02_wp, -0.890739e+01_wp, 0.191135e+02_wp, -0.512653e+02_wp, 0.365096e+02_wp, &
            & - 0.235972e+02_wp, -0.102322e+03_wp, 0.837770e+00_wp, 0.258200e+02_wp, -0.198105e+01_wp, - &
            & 0.484189e+01_wp, 0.105377e+02_wp, 0.372352e+01_wp, -0.207237e+02_wp, 0.223342e+02_wp, 0.244787e+01_wp, &
            & 0.979572e+01_wp, -0.135612e+01_wp, -0.724496e+01_wp, 0.234984e+02_wp, 0.136364e+03_wp, 0.253799e+01_wp, &
            & -0.203821e+02_wp, 0.699846e+01_wp, 0.198529e+02_wp, -0.475902e+01_wp, 0.236328e+01_wp]
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
                    g(1,iatom,idir) = usrocm_dgscart(iatom,idir)*hartree_to_electronvolt
                enddo
            enddo
        elseif( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif

        contains

        !> @brief Evaluate prepare_local_context for the OH2_GEN_SL1979 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, j
            !        IMPLICIT REAL*16(L)
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            real(wp), parameter :: bohr   = 0.529177e0_wp
            real(wp), parameter :: hart   = 27.21161e0_wp
            !      DIMENSION C(50),AL(3),X(3),S(3),VM(3),R(3),CHI(3)
            !   Echo potential data to the file linked to UNIT NFLAG(18)
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !      WRITE (NFLAG(18),600)
            !      WRITE (NFLAG(18),601) STEP
            !      WRITE (NFLAG(18),602)
            j = 1
            !    Set the values of the classical energy in the three asymptotic vall
            asycm_easyab   = satocm_de(1)/hart
            asycm_easybc   = satocm_de(2)/hart
            asycm_easyac   = satocm_de(3)/hart
            pt3cm_ezero(1) = satocm_de(2)
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

        !> @brief Evaluate POT for the OH2_GEN_SL1979 potential.
        subroutine pot
            !     CALCULATE ENERGY AND DERIVATIVE
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
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: ap, bp, check, d1d, d2d, d3d, da, da2, da3, dap1, dap2, dap3, db, db1, db2, db3, dc, dc2, dc3, dchi1, dchi2, dchi3, dcsr1, dcsr2, dcsr3, dd, dd2, dd3, df, df2, df3, dg, dq1, dq2, dq3, dvr1, dvr2, dvr3, gg, p, psi13, q, soh, sr2, t, temp, tt, u, w, ww, xter, y, z, zz
            integer  :: i, ic, j
            !        IMPLICIT REAL*16(L)
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            real(wp), parameter :: bohr   = 0.529177e0_wp
            real(wp), parameter :: hart   = 27.21161e0_wp
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
            ic = 4
            if( pt1cm_r(1) > 30.0e0_wp .or. pt1cm_r(2) > 30.0e0_wp ) ic = 1
            evaluation_loop: do
                select case( ic )
                case( 1 )
                    sr2        = pt1cm_r(2)
                    soh        = pt1cm_r(1)
                    pt1cm_r(3) = pt1cm_r(3)+ndercm_step
                case( 2 )
                    sr2 = pt1cm_r(2)+ndercm_step
                    soh = pt1cm_r(1)
                case( 3 )
                    sr2 = pt1cm_r(2)
                    soh = pt1cm_r(1)+ndercm_step
                case( 4 )
                    soh = pt1cm_r(1)
                    sr2 = pt1cm_r(2)
                case default
                    error stop 'OH2 SL: invalid finite-difference state'
                end select
                sr2        = sr2*bohr
                soh        = soh*bohr
                pt1cm_r(3) = pt1cm_r(3)*bohr
                pt1cm_r(1) = soh
                pt1cm_r(2) = sr2
                q          = 1.0e0_wp
                do i=1, 3
                    precm_x(i)   = pt1cm_r(i)
                    precm_chi(i) = (1.0e0_wp - exp(-precm_alf*precm_x(i)**2))
                    precm_s(i)   = precm_x(i)-1.e0_wp
                    precm_vm(i)  = satocm_de(i)*(1.0e0_wp-exp(-satocm_bet(i)*(precm_x(i)-satocm_xre(i))))**2- &
                        & satocm_de(i)
                    q = q*(1.0e0_wp-tanh(precm_al(i)*precm_s(i)))
                enddo
                u = precm_vm(1)*precm_chi(2)*precm_chi(3)+precm_vm(2)*precm_chi(1)*precm_chi(3) + precm_vm(3)* &
                    & precm_chi(1)*precm_chi(2)
                bp    = precm_b*(1.0e0_wp/precm_x(1)**4+1.0e0_wp/precm_x(3)**4)
                check = (precm_x(1)*precm_x(3)*(precm_x(2)-0.2e0_wp))
                if( check > 1000.0e0_wp )then
                    ap = 0.0e0_wp
                else
                    ap = precm_a/(precm_x(1)*precm_x(3)*(precm_x(2)-0.2e0_wp))**12
                endif
                u = u+ap+bp
                t = precm_c(1) + precm_c(2)*(precm_s(1)+precm_s(3)) + precm_c(3)*precm_s(2) + precm_c(4)*(precm_s(1)** &
                    & 2 + precm_s(3)**2) + precm_c(5)*precm_s(2)**2 + precm_c(6)*precm_s(2)*(precm_s(1) + &
                    & precm_s(3)) + precm_c(7)* precm_s(1)*precm_s(3) + precm_c(8)*(precm_s(1)**3 + precm_s(3)**3) + &
                    & precm_c(9)* precm_s(2)**3 + precm_c(10)*precm_s(1)*precm_s(3)*(precm_s(1) + precm_s(3))
                tt = precm_c(11)*precm_s(2)**2*(precm_s(1) + precm_s(3)) + precm_c(12)*precm_s(2)*(precm_s(1)**2 + &
                    & precm_s(3)**2) + precm_c(13)*precm_s(1)*precm_s(2)*precm_s(3) + precm_c(14)*(precm_s(1)**4 + &
                    & precm_s(3)** 4) + precm_c(15)*precm_s(2)**4 + precm_c(16)*(precm_s(1)**3*precm_s(3) + &
                    & precm_s(1)* precm_s(3)**3)
                ww = precm_c(17)*precm_s(1)**2*precm_s(3)**2 + precm_c(18)*precm_s(2)*(precm_s(1)**3 + precm_s(3)** &
                    & 3) + precm_c(19)*precm_s(2)**2*(precm_s(1)**2 + precm_s(3)**2) + precm_c(20)*precm_s(2)**3* &
                    & (precm_s(1) + precm_s(3)) + precm_c(21)*precm_s(2)**2*precm_s(1)*precm_s(3) + precm_c(22)* &
                    & precm_s(2)*(precm_s(1)** 2*precm_s(3) + precm_s(1)*precm_s(3)**2)
                w = precm_c(23)*(precm_s(1)**5 + precm_s(3)**5) + precm_c(24)*(precm_s(1)**4*precm_s(3) + precm_s(1)* &
                    & precm_s(3)**4) + precm_c(25)*(precm_s(1)**3*precm_s(3)**2 + precm_s(1)**2*precm_s(3)**3) + &
                    & precm_c(26)*precm_s(2)*(precm_s(1)**4 + precm_s(3)**4) + precm_c(27)*precm_s(2)*(precm_s(1)**3* &
                    & precm_s(3) + precm_s(1)*precm_s(3)**3)
                gg = precm_c(28)*precm_s(2)*precm_s(1)**2*precm_s(3)**2 + precm_c(29)*precm_s(2)**2*(precm_s(1)**3 + &
                    & precm_s(3)**3) + precm_c(30)*precm_s(2)**2*(precm_s(1)**2*precm_s(3) + precm_s(1)*precm_s(3)** &
                    & 2) + precm_c(31)*precm_s(2)**3*(precm_s(1)**2 + precm_s(3)**2)
                y = precm_c(32)*precm_s(2)**3*precm_s(1)*precm_s(3) + precm_c(33)*precm_s(2)**4*(precm_s(1) + &
                    & precm_s(3)) + precm_c(34)*precm_s(2)**5 + precm_c(35)*(precm_s(1)**6 + precm_s(3)**6) + &
                    & precm_c(36)* (precm_s(1)**5*precm_s(3) + precm_s(1)*precm_s(3)**5) + precm_c(37)*(precm_s(1)**4* &
                    & precm_s(3)**2 + precm_s(1)**2*precm_s(3)**4) + precm_c(38)*precm_s(1)**3*precm_s(3)**3
                z = precm_c(39)*precm_s(2)*(precm_s(1)**5 + precm_s(3)**5) + precm_c(40)*precm_s(2)*(precm_s(1)**4* &
                    & precm_s(3) + precm_s(1)*precm_s(3)**4) + precm_c(41)*precm_s(2)*(precm_s(1)**3*precm_s(3)**2 + &
                    & precm_s(1)** 2*precm_s(3)**3) + precm_c(42)*precm_s(2)**2*(precm_s(1)**4 + precm_s(3)**4) + &
                    & precm_c(43)*precm_s(2)**2*(precm_s(1)**3*precm_s(3) + precm_s(1)*precm_s(3)**3) + precm_c(44)* &
                    & precm_s(1)**2* precm_s(2)**2*precm_s(3)**2
                zz = precm_c(45)*precm_s(2)**3*(precm_s(1)**3 + precm_s(3)**3) + precm_c(46)*precm_s(2)**3* &
                    & (precm_s(1)**2* precm_s(3) + precm_s(1)*precm_s(3)**2) + precm_c(47)*precm_s(2)**4*(precm_s(1)** &
                    & 2 + precm_s(3)**2) + precm_c(48)*precm_s(2)**4*precm_s(1)*precm_s(3) + precm_c(49)*precm_s(2)** &
                    & 5*(precm_s(1) + precm_s(3)) + precm_c(50)*precm_s(2)**6
                p           = t + tt + ww + w + gg + y + z + zz
                temp        = exp(-precm_gam*(precm_x(1) + precm_x(3))**2)
                u           = u + p*q*exp(-precm_gam*(precm_x(1) + precm_x(3))**2)
                xter        = u + pt3cm_ezero(1)
                precm_v(ic) = xter/hart
                if( ic < 4 )then
                    ic = ic+1
                    cycle evaluation_loop
                endif
                if( pt1cm_r(1) > 30.0e0_wp .or. pt1cm_r(2) >= 30.0e0_wp )then
                    if( usricm_nder == 1 )then
                        pt1cm_degsdr(1) = (precm_v(3)-precm_v(4))/ndercm_step
                        pt1cm_degsdr(2) = (precm_v(2)-precm_v(4))/ndercm_step
                        pt1cm_degsdr(3) = (precm_v(1)-precm_v(4))/ndercm_step
                    endif
                else
                    ! Analytic derivatives in the interaction region.
                    if( usricm_nder .eq. 1 )then
                        da = precm_c(2) + 2.0e0_wp*precm_s(1)*precm_c(4) + precm_s(2)*precm_c(6) + precm_s(3)* &
                            & precm_c(7) + 3.0e0_wp* precm_s(1)**2*precm_c(8) + (precm_s(3)**2 + 2.0e0_wp*precm_s(1)* &
                            & precm_s(3))*precm_c(10) + precm_s(2)**2*precm_c(11) + 2.0e0_wp*precm_s(1)*precm_s(2)* &
                            & precm_c(12) + precm_s(2)*precm_s(3)* precm_c(13) + 4.0e0_wp*precm_s(1)**3*precm_c(14)
                        db = (3.0e0_wp*precm_s(1)**2*precm_s(3) + precm_s(3)**3)*precm_c(16) + 2.0e0_wp*precm_s(1)* &
                            & precm_s(3)**2* precm_c(17) + 3.0e0_wp*precm_s(1)**2*precm_s(2)*precm_c(18) + 2.0e0_wp* &
                            & precm_s(1)* precm_s(2)**2* precm_c(19) + precm_s(2)**3*precm_c(20) + precm_s(2)**2* &
                            & precm_s(3)*precm_c(21) + (2.0e0_wp* precm_s(1)*precm_s(2)*precm_s(3) + precm_s(2)* &
                            & precm_s(3)**2)*precm_c(22)
                        dc = 5.0e0_wp*precm_s(1)**4*precm_c(23) + (4.0e0_wp*precm_s(1)**3*precm_s(3) + precm_s(3)**4)* &
                            & precm_c(24) + (3.0e0_wp*precm_s(1)**2*precm_s(3)**2 + 2.0e0_wp*precm_s(1)*precm_s(3)** &
                            & 3)* precm_c(25) + 4.0e0_wp* precm_s(1)**3*precm_s(2)*precm_c(26) + (3.0e0_wp*precm_s(2)* &
                            & precm_s(1)**2*precm_s(3) + precm_s(2)* precm_s(3)**3)*precm_c(27)
                        dd = 2.0e0_wp*precm_s(1)*precm_s(2)*precm_s(3)**2*precm_c(28) + 3.0e0_wp*precm_s(1)**2* &
                            & precm_s(2)** 2* precm_c(29) + (2.0e0_wp*precm_s(1)*precm_s(2)**2*precm_s(3) + &
                            & precm_s(2)**2*precm_s(3)** 2)* precm_c(30) + 2.0e0_wp*precm_s(1)*precm_s(2)**3* &
                            & precm_c(31) + precm_s(2)**3*precm_s(3)* precm_c(32) + precm_s(2)**4*precm_c(33) + &
                            & 6.0e0_wp*precm_s(1)**5*precm_c(35) + (5.0e0_wp* precm_s(1)** 4*precm_s(3) + precm_s(3)** &
                            & 5)*precm_c(36)
                        df = (4.0e0_wp*precm_s(1)**3*precm_s(3)**2 + 2.0e0_wp*precm_s(1)*precm_s(3)**4)*precm_c(37) + &
                            & 3.0e0_wp* precm_s(1)**2*precm_s(3)**3*precm_c(38) + 5.0e0_wp*precm_s(1)**4*precm_s(2)* &
                            & precm_c(39) + (4.0e0_wp* precm_s(2)*precm_s(1)**3*precm_s(3) + precm_s(2)*precm_s(3)** &
                            & 4)*precm_c(40) + (3.0e0_wp* precm_s(1)**2*precm_s(2)*precm_s(3)**2 + 2.0e0_wp* &
                            & precm_s(1)*precm_s(2)* precm_s(3)**3)* precm_c(41) + 4.0e0_wp*precm_s(1)**3*precm_s(2)** &
                            & 2*precm_c(42)
                        dg = (3.0e0_wp*precm_s(1)**2*precm_s(2)**2*precm_s(3) + precm_s(2)**2*precm_s(3)**3)* &
                            & precm_c(43) + 2.0e0_wp* precm_s(1)*precm_s(2)**2*precm_s(3)**2*precm_c(44) + 3.0e0_wp* &
                            & precm_s(1)**2* precm_s(2)**3* precm_c(45) + (2.0e0_wp*precm_s(1)*precm_s(2)**3* &
                            & precm_s(3) + precm_s(2)**3*precm_s(3)**2)* precm_c(46) + 2.0e0_wp*precm_s(1)* &
                            & precm_s(2)**4*precm_c(47) + precm_s(2)**4* precm_s(3)* precm_c(48) + precm_s(2)**5* &
                            & precm_c(49)
                        d1d = da + db + dc + dd + df + dg
                        da2 = precm_c(3) + 2.0e0_wp*precm_s(2)*precm_c(5) + (precm_s(1) + precm_s(3))*precm_c(6) + &
                            & 3.0e0_wp* precm_s(2)**2*precm_c(9) + 2.0e0_wp*precm_s(2)*(precm_s(1) + precm_s(3))* &
                            & precm_c(11) + (precm_s(1)**2 + precm_s(3)**2)*precm_c(12) + precm_s(1)*precm_s(3)* &
                            & precm_c(13) + 4.0e0_wp* precm_s(2)**3*precm_c(15) + (precm_s(1)**3 + precm_s(3)**3)* &
                            & precm_c(18) + 2.0e0_wp* precm_s(2)* (precm_s(1)**2 + precm_s(3)**2)*precm_c(19)
                        db2 = 3.0e0_wp*precm_s(2)**2*(precm_s(1) + precm_s(3))*precm_c(20) + 2.0e0_wp*precm_s(1)* &
                            & precm_s(2)* precm_s(3)*precm_c(21) + (precm_s(1)**2*precm_s(3) + precm_s(1)*precm_s(3)** &
                            & 2)* precm_c(22) + (precm_s(1)**4 + precm_s(3)**4)*precm_c(26) + (precm_s(1)**3* &
                            & precm_s(3) + precm_s(1)* precm_s(3)**3)*precm_c(27) + precm_s(1)**2*precm_s(3)**2* &
                            & precm_c(28) + 2.0e0_wp*precm_s(2)* (precm_s(1)**3 + precm_s(3)**3)*precm_c(29)
                        dc2 = 2.0e0_wp*precm_s(2)*(precm_s(1)**2*precm_s(3) + precm_s(1)*precm_s(3)**2)*precm_c(30) + &
                            & 3.0e0_wp* precm_s(2)**2*(precm_s(1)**2 + precm_s(3)**2)*precm_c(31) + 3.0e0_wp* &
                            & precm_s(2)** 2* precm_s(1)* precm_s(3)*precm_c(32) + 4.0e0_wp*precm_s(2)**3* &
                            & (precm_s(1) + precm_s(3))* precm_c(33) + 5.0e0_wp* precm_s(2)**4*precm_c(34)
                        dd2 = (precm_s(1)**5 + precm_s(3)**5)*precm_c(39) + (precm_s(1)**4*precm_s(3) + precm_s(1)* &
                            & precm_s(3)**4)*precm_c(40) + (precm_s(1)**3*precm_s(3)**2 + precm_s(1)**2*precm_s(3)** &
                            & 3)* precm_c(41) + 2.0e0_wp*precm_s(2)*(precm_s(1)**4 + precm_s(3)**4)*precm_c(42) + &
                            & 2.0e0_wp* precm_s(2)* (precm_s(1)**3*precm_s(3) + precm_s(1)*precm_s(3)**3)* &
                            & precm_c(43) + 2.0e0_wp* precm_s(1)**2* precm_s(2)*precm_s(3)**2*precm_c(44)
                        df2 = 3.0e0_wp*precm_s(2)**2*(precm_s(1)**3 + precm_s(3)**3)*precm_c(45) + 3.0e0_wp* &
                            & precm_s(2)**2* (precm_s(1)**2*precm_s(3) + precm_s(1)*precm_s(3)**2)*precm_c(46) + &
                            & 4.0e0_wp* precm_s(2)** 3* (precm_s(1)**2 + precm_s(3)**2)*precm_c(47) + 4.0e0_wp* &
                            & precm_s(2)**3*precm_s(1)* precm_s(3)* precm_c(48) + 5.0e0_wp*precm_s(2)**4*(precm_s(1) + &
                            & precm_s(3))*precm_c(49) + 6.0e0_wp* precm_s(2)**5* precm_c(50)
                        d2d = da2 + db2 + dc2 + dd2 + df2
                        da3 = precm_c(2) + 2.0e0_wp*precm_s(3)*precm_c(4) + precm_s(2)*precm_c(6) + precm_s(1)* &
                            & precm_c(7) + 3.0e0_wp*precm_s(3)**2*precm_c(8) + (2.0e0_wp*precm_s(1)*precm_s(3) + &
                            & precm_s(1)**2)* precm_c(10) + precm_s(2)**2*precm_c(11) + 2.0e0_wp*precm_s(2)* &
                            & precm_s(3)* precm_c(12) + precm_s(1)*precm_s(2)* precm_c(13) + 4.0e0_wp*precm_s(3)**3* &
                            & precm_c(14)
                        db3 = (precm_s(1)**3 + 3.0e0_wp*precm_s(1)*precm_s(3)**2)*precm_c(16) + 2.0e0_wp*precm_s(1)** &
                            & 2* precm_s(3)* precm_c(17) + 3.0e0_wp*precm_s(2)*precm_s(3)**2*precm_c(18) + 2.0e0_wp* &
                            & precm_s(2)** 2*precm_s(3)* precm_c(19) + precm_s(2)**3*precm_c(20) + precm_s(2)**2* &
                            & precm_s(1)*precm_c(21) + precm_s(2)* (precm_s(1)**2 + 2.0e0_wp*precm_s(1)*precm_s(3))* &
                            & precm_c(22) + 5.0e0_wp*precm_s(3)** 4*precm_c(23) + (precm_s(1)**4 + 4.0e0_wp* &
                            & precm_s(1)* precm_s(3)**3)*precm_c(24)
                        dc3 = (2.0e0_wp*precm_s(1)**3*precm_s(3) + 3.0e0_wp*precm_s(1)**2*precm_s(3)**2)*precm_c(25) + &
                            & 4.0e0_wp* precm_s(2)*precm_s(3)**3*precm_c(26) + precm_s(2)*(precm_s(1)**3 + 3.0e0_wp* &
                            & precm_s(1)* precm_s(3)**2)*precm_c(27) + 2.0e0_wp*precm_s(2)*precm_s(1)**2*precm_s(3)* &
                            & precm_c(28) + 3.0e0_wp* precm_s(2)**2*precm_s(3)**2*precm_c(29)
                        dd3 = precm_s(2)**2*(precm_s(1)**2 + 2.0e0_wp*precm_s(1)*precm_s(3))*precm_c(30) + 2.0e0_wp* &
                            & precm_s(2)**3* precm_s(3)*precm_c(31) + precm_s(2)**3*precm_s(1)*precm_c(32) + &
                            & precm_s(2)**4* precm_c(33) + 6.0e0_wp*precm_s(3)**5*precm_c(35) + (precm_s(1)**5 + &
                            & 5.0e0_wp* precm_s(1)*precm_s(3)**4)*precm_c(36) + (2.0e0_wp*precm_s(1)**4*precm_s(3) + &
                            & 4.0e0_wp* precm_s(1)** 2* precm_s(3)**3)*precm_c(37) + 3.0e0_wp* precm_s(1)**3* &
                            & precm_s(3)**2* precm_c(38)
                        df3 = 5.0e0_wp*precm_s(2)*precm_s(3)**4*precm_c(39) + precm_s(2)*(precm_s(1)**4 + 4.0e0_wp* &
                            & precm_s(1)* precm_s(3)**3)*precm_c(40) + precm_s(2)*(2.0e0_wp*precm_s(1)**3*precm_s(3) + &
                            & 3.0e0_wp* precm_s(1)**2* precm_s(3)**2)*precm_c(41) + 4.0e0_wp*precm_s(2)**2* &
                            & precm_s(3)**3* precm_c(42) + precm_s(2)**2* (precm_s(1)**3 + 3.0e0_wp*precm_s(1)* &
                            & precm_s(3)**2)* precm_c(43) + 2.0e0_wp*precm_s(1)**2*precm_s(2)** 2*precm_s(3)* &
                            & precm_c(44) + 3.0e0_wp* precm_s(2)**3*precm_s(3)**2*precm_c(45) + precm_s(2)**3* &
                            & (precm_s(1)**2 + 2.0e0_wp* precm_s(1)*precm_s(3))* precm_c(46) + 2.0e0_wp*precm_s(2)**4* &
                            & precm_s(3)* precm_c(47) + precm_s(2)**4*precm_s(1)*precm_c(48) + precm_s(2)**5*precm_c(49)
                        d3d  = da3 + db3 + dc3 + dd3 + df3
                        dq1  = (1.0e0_wp - tanh(precm_al(2)*precm_s(2)))*(1.0e0_wp-tanh(precm_al(3)*precm_s(3)))
                        dq2  = (1.0e0_wp - tanh(precm_al(1)*precm_s(1)))*(1.0e0_wp-tanh(precm_al(3)*precm_s(3)))
                        dq3  = (1.0e0_wp - tanh(precm_al(1)*precm_s(1)))*(1.0e0_wp-tanh(precm_al(2)*precm_s(2)))
                        dap1 = -12.0e0_wp*ap/pt1cm_r(1)
                        dap2 = -12.0e0_wp*ap/(pt1cm_r(2) - 0.2e0_wp)
                        dap3 = -12.0e0_wp*ap/pt1cm_r(3)
                        db1  = -4.0e0_wp*precm_b/pt1cm_r(1)**5
                        db3  = -4.0e0_wp*precm_b/pt1cm_r(3)**5
                        !      IF(DCOSH(AL(1)*S(1)).LT.1.0D-15)DCRS1=1.0D30
                        !      IF(DCOSH(AL(1)*S(1)).LT.1.0D-15)GO TO 144
                        dcsr1 = -precm_al(1)/(cosh(precm_al(1)*precm_s(1))**2)
                        !144   IF(COSH(AL(2)*S(2)).LT.1.0D-15)DCRS2=1.0D30
                        !      IF(COSH(AL(2)*S(2)).LT.1.0D-15)GO TO 145
                        dcsr2 = -precm_al(2)/(cosh(precm_al(2)*precm_s(2))**2)
                        !145   IF(COSH(AL(3)*S(3)).LT.1.0D-15)DCRS3=1.0D30
                        !      IF(COSH(AL(3)*S(3)).LT.1.0D-15)GO TO 146
                        dcsr3 = -precm_al(3)/(cosh(precm_al(3)*precm_s(3))**2)
                        psi13 = 2.0e0_wp*(-precm_gam)*temp*(pt1cm_r(1) + pt1cm_r(3))
                        dchi1 = 2.0e0_wp*precm_alf*pt1cm_r(1)*exp(-precm_alf*(pt1cm_r(1)**2))
                        dchi2 = 2.0e0_wp*precm_alf*pt1cm_r(2)*exp(-precm_alf*(pt1cm_r(2)**2))
                        dchi3 = 2.0e0_wp*precm_alf*pt1cm_r(3)*exp(-precm_alf*(pt1cm_r(3)**2))
                        !      IF(VM(1).GT.1.0D20)DVR1=1.0D30
                        !      IF(VM(1).GT.1.0D20)GO TO 153
                        dvr1 = 2.0e0_wp*satocm_de(1)*(1.0e0_wp - exp(-satocm_bet(1)*(pt1cm_r(1)-satocm_xre(1))))* &
                            & satocm_bet(1)* exp(-satocm_bet(1)*(pt1cm_r(1)-satocm_xre(1)))
                        ! 153  IF(VM(2).GT.1.0D20)DVR2=1.0D30
                        !      IF(VM(2).GT.1.0D20)GO TO 154
                        dvr2 = 2.0e0_wp*satocm_de(2)*(1.0e0_wp - exp(-satocm_bet(2)*(pt1cm_r(2)-satocm_xre(2))))* &
                            & satocm_bet(2)* exp(-satocm_bet(2)*(pt1cm_r(2)-satocm_xre(2)))
                        ! 154  IF(VM(3).GT.1.0D20)DVR3=1.0D30
                        !      IF(VM(3).GT.1.0D20)GO TO 155
                        dvr3 = 2.0e0_wp*satocm_de(3)*(1.0e0_wp - exp(-satocm_bet(3)*(pt1cm_r(3)-satocm_xre(3))))* &
                            & satocm_bet(3)* exp(-satocm_bet(3)*(pt1cm_r(3)-satocm_xre(3)))
                        pt1cm_degsdr(1) = precm_vm(2)*precm_chi(3)*dchi1 + precm_chi(2)*precm_chi(3)*dvr1 + &
                            & precm_vm(3)*precm_chi(2)*dchi1 + d1d*q*temp + p*dq1*temp*dcsr1 + dap1 + db1 + p*q*psi13
                        pt1cm_degsdr(2) = dvr2*precm_chi(1)*precm_chi(3) + precm_vm(1)*precm_chi(3)*dchi2 + &
                            & precm_vm(3)* precm_chi(1)*dchi2 + d2d*q*temp + p*dq2*temp*dcsr2 + dap2
                        pt1cm_degsdr(3) = precm_vm(2)*precm_chi(1)*dchi3 + precm_vm(1)*precm_chi(2)*dchi3 + dvr3* &
                            & precm_chi(1)* precm_chi(2) + d3d*q*temp + p*dq3*temp*dcsr3 + p* q* psi13 + dap3 + db3
                        do i=1, 3
                            pt1cm_degsdr(i) = pt1cm_degsdr(i)/hart*bohr
                        enddo
                        ! END OF NDER SWITCH AROUND DERIVATIVES
                    endif
                endif
                if( pt1cm_r(2) >= 20.0e0_wp )then
                    precm_v(4) = int(precm_v(4)*1.0e8_wp)*1.0e-8_wp
                endif
                exit evaluation_loop
            enddo evaluation_loop
            pt1cm_engygs = precm_v(4)
            call eunitzero
            if( usricm_nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
            return
            stop
        end subroutine pot

        !> @brief Evaluate POTINFO for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate ANCVRT for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate CARTOU for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate CARTTOR for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate EUNITZERO for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate RTOCART for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate DEDCOU for the OH2_GEN_SL1979 potential.
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
end module chempotpy_oh2_gen_sl1979_pes

!> @brief Expose OH2_GEN_SL1979 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh2_gen_sl1979_pes, only: evaluate
    use chempotpy_kinds,              only: wp
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
