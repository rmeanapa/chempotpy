!   System:          BrH2
!   Functional form: extended LEPS plus three-center term
!   Common name:     BrH2SO
!   Interface:       potlib2001
!   Reference:        G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g.
!                     J. Phys. Chem. 99, 207 (1995)
!   Number of bodies:              3
!   Number of derivatives:         0
!   Number of electronic surfaces: 2

!> @brief Implement BrH2_LEPS_LTBZSO_1995, a potential energy surface for BrH2.
module chempotpy_brh2_leps_ltbzso_1995_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the BrH2_LEPS_LTBZSO_1995 potential energy surface for BrH2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 2
        integer, parameter    :: natoms  = 3
        integer, parameter    :: n3atom  = 75
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: anuzero
        real(wp)              :: beta11(3)
        real(wp)              :: beta22(3)
        real(wp)              :: betaa11(3)
        real(wp)              :: betaa22(3)
        real(wp)              :: cart(natom,3)
        real(wp)              :: cartnu(natom,3)
        real(wp)              :: ccon
        real(wp)              :: cnvrtd
        real(wp)              :: cnvrtde
        real(wp)              :: cnvrte
        real(wp)              :: coul(3)
        real(wp)              :: d11(3)
        real(wp)              :: d22(3)
        real(wp)              :: de11(3)
        real(wp)              :: de22(3)
        real(wp)              :: deesdr(n3atom,isurf)
        real(wp)              :: degsdr(n3atom)
        real(wp)              :: deijdr(n3atom,jsurf)
        real(wp)              :: deq(3)
        real(wp)              :: derad(3)
        real(wp)              :: descart(natom,3,isurf)
        real(wp)              :: descartnu(natom,3,isurf)
        real(wp)              :: dgscart(natom,3)
        real(wp)              :: dgscartnu(natom,3)
        real(wp)              :: dijcart(natom,3,jsurf)
        real(wp)              :: dijcartnu(natom,3,jsurf)
        real(wp)              :: do4z(3)
        real(wp)              :: dv3c(3)
        real(wp)              :: dxdr(3,3)
        real(wp)              :: dz(3,3)
        real(wp)              :: easyac
        real(wp)              :: easyba
        real(wp)              :: easybc
        real(wp)              :: engyes(isurf)
        real(wp)              :: engygs
        real(wp)              :: engyij(jsurf)
        real(wp)              :: exch(3)
        real(wp)              :: ezero(isurf+1)
        integer               :: icartr
        integer               :: indexes(natom)
        integer               :: irctnt
        integer               :: ireorder
        integer               :: kediag
        integer               :: keoffd
        integer               :: ksdiag
        integer               :: ksoffd
        integer               :: legacy_natoms
        real(wp)              :: legacy_x(3)
        integer               :: mder
        integer               :: msurf
        integer               :: nasurf(isurf+1,isurf+1)
        integer               :: nder
        integer               :: nflag(20)
        integer               :: nulbl(natom)
        real(wp)              :: op3z(3)
        real(wp)              :: pengyes(isurf)
        real(wp)              :: pengygs
        real(wp)              :: pengyij(jsurf)
        real(wp)              :: r(n3atom)
        real(wp)              :: r211
        real(wp)              :: r222
        real(wp)              :: rad
        real(wp)              :: re11(3)
        real(wp)              :: re22(3)
        character             :: ref(5)
        real(wp)              :: top3z(3)
        real(wp)              :: tzp3(3)
        real(wp)              :: v3c
        real(wp)              :: v3cp1
        real(wp)              :: v3cp2
        real(wp)              :: v3cp3
        real(wp)              :: v3cp4
        real(wp)              :: v3cp5
        real(wp)              :: v3cp6
        real(wp)              :: z(3)
        real(wp)              :: zhbrp1
        real(wp)              :: zhbrp2
        real(wp)              :: zhhp1
        real(wp)              :: zhhp2
        real(wp)              :: zp3(3)
        real(wp)              :: zpo(3)
        integer               :: iatom, idir
        real(wp)              :: u(nstates,nstates)
        real(wp)              :: t(nstates,nstates)
        !initialize
        anuzero   = 0.e0_wp; beta11 = 0.e0_wp; beta22 = 0.e0_wp; betaa11 = 0.e0_wp; betaa22 = 0.e0_wp; cart = 0.e0_wp; cartnu = 0.e0_wp
        ccon      = 0.e0_wp; cnvrtd = 0.e0_wp; cnvrtde = 0.e0_wp; cnvrte = 0.e0_wp; coul = 0.e0_wp; d11 = 0.e0_wp; d22 = 0.e0_wp; de11 = 0.e0_wp
        de22      = 0.e0_wp; deesdr = 0.e0_wp; degsdr = 0.e0_wp; deijdr = 0.e0_wp; deq = 0.e0_wp; derad = 0.e0_wp; descart = 0.e0_wp
        descartnu = 0.e0_wp; dgscart = 0.e0_wp; dgscartnu = 0.e0_wp; dijcart = 0.e0_wp; dijcartnu = 0.e0_wp; do4z = 0.e0_wp; dv3c = 0.e0_wp
        dxdr      = 0.e0_wp; dz = 0.e0_wp; easyac = 0.e0_wp; easyba = 0.e0_wp; easybc = 0.e0_wp; engyes = 0.e0_wp; engygs = 0.e0_wp; engyij = 0.e0_wp
        exch      = 0.e0_wp; ezero = 0.e0_wp; icartr = 0; indexes = 0; irctnt = 0; ireorder = 0; kediag = 0; keoffd = 0; ksdiag = 0
        ksoffd    = 0; legacy_natoms = 0; legacy_x = 0.e0_wp; mder = 0; msurf = 0; nasurf = 0; nder = 0; nflag = 0; nulbl = 0
        op3z      = 0.e0_wp; pengyes = 0.e0_wp; pengygs = 0.e0_wp; pengyij = 0.e0_wp; r = 0.e0_wp; r211 = 0.e0_wp; r222 = 0.e0_wp; rad = 0.e0_wp
        re11      = 0.e0_wp; re22 = 0.e0_wp
        ref       = ' '
        top3z     = 0.e0_wp
        tzp3      = 0.e0_wp
        v3c       = 0.e0_wp
        v3cp1     = 0.e0_wp
        v3cp2     = 0.e0_wp
        v3cp3     = 0.e0_wp
        v3cp4     = 0.e0_wp
        v3cp5     = 0.e0_wp
        v3cp6     = 0.e0_wp
        z         = 0.e0_wp
        zhbrp1    = 0.e0_wp
        zhbrp2    = 0.e0_wp
        zhhp1     = 0.e0_wp
        zhhp2     = 0.e0_wp
        zp3       = 0.e0_wp
        zpo       = 0.e0_wp
        nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0], shape(nasurf))
        nder          = 0
        nflag         = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        anuzero       = 0.0e0_wp
        icartr        = 3
        msurf         = 1
        mder          = 0
        nulbl         = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        legacy_natoms = 3
        d11           = [90.4473e0_wp, 109.559e0_wp, 90.4473e0_wp]
        re11          = [1.41443e0_wp, 0.74144e0_wp, 1.41443e0_wp]
        betaa11       = [1.81088e0_wp, 1.94198e0_wp, 1.81088e0_wp]
        r211          = 1.41421356e0_wp
        d22           = [90.4473e0_wp, 109.559e0_wp, 90.4473e0_wp]
        re22          = [1.41443e0_wp, 0.74144e0_wp, 1.41443e0_wp]
        betaa22       = [1.81088e0_wp, 1.94198e0_wp, 1.81088e0_wp]
        r222          = 1.41421356e0_wp
        ccon          = 0.5e0_wp
        zhhp1         = 2.91230499e-1_wp
        zhhp2         = 4.99160974e0_wp
        zhbrp1        = 3.84262775e-3_wp
        zhbrp2        = 2.87642885e-3_wp
        v3cp1         = 1.67457884e0_wp
        v3cp2         = 7.09518875e-1_wp
        v3cp3         = 1.87582679e-1_wp
        v3cp4         = 2.43517074e-1_wp
        v3cp5         = 2.04145206e-1_wp
        v3cp6         = 2.57833022e0_wp
        u             = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        nasurf        = 1
        cart          = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                cart(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        nder = igrad
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            u(1,1) = engygs*hartree_to_electronvolt
            u(2,2) = engyes(1)*hartree_to_electronvolt
            u(1,2) = engyij(1)*hartree_to_electronvolt
            u(2,1) = engyij(1)*hartree_to_electronvolt
            call diagonalize(nstates, u, p, t)
        else
            write (*, *) 'Only energy is available'
        endif

        contains

        !> @brief Evaluate diagonalize for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine diagonalize( n, a_ss, ev_s, u_ss )
            use chempotpy_kinds, only: wp
            implicit none
            integer, intent(in)   :: n
            real(wp), intent(in)  :: a_ss(n,n)
            real(wp), intent(out) :: u_ss(n,n)
            real(wp), intent(out) :: ev_s(n)
            integer               :: io, i
            real(wp), allocatable :: work_d(:)
            integer               :: lwork_d
            u_ss    = a_ss
            lwork_d = -1
            allocate( work_d(2) )
            call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
            lwork_d = int(work_d(1))
            deallocate ( work_d )
            allocate( work_d(lwork_d) )
            call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
        end subroutine diagonalize

        !> @brief Evaluate prepare_local_context for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            if( legacy_natoms .gt. 25 )then
                write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            call prepare_local_context_11
            call prepare_local_context_12
            call prepare_local_context_22
            ezero(1) = de11(2)
            ezero(2) = de22(2)
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1)     = 'No Reference'
            indexes(1) = 35
            indexes(2) = 1
            indexes(3) = 1
            irctnt     = 2
            call ancvrt
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine pot
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            call cartou
            call carttor
            if( nasurf(1,1) .ne. 0 ) call pot11
            if( nasurf(2,2) .ne. 0 ) call pot22
            if( nasurf(1,2)+nasurf(2,1) .ne. 0 ) call pot12
            call eunitzero
            if( nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate prepare_local_context_11 for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine prepare_local_context_11
            !   System:          BrH2
            !   Functional form: extended-LEPS plus three-center term
            !   Common name:     BrH2SEC
            !   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, et. al.
            !                    to be published
            !   prepare_local_context must be called once before any calls to POT.
            !   The potential parameters are included in DATA statements.
            !   Coordinates, potential energy, and derivatives are passed
            !   through the common block POTCM:
            !   The potential energy in the three asymptotic valleys are
            !   stored in the common block ASYCOM:
            !   For the this potential energy surface the potential energy in the
            !   BC valley is equal to the energy of the HH diatomic.
            !   The energy in the AC valley is equal to the energy of the HBr diatom
            !   which is also the energy value in the BA valley.
            !   All the information passed through the common blocks PV1CM and ASYCM
            !   is in Hartree atomic units.
            !   This potential is written such that:
            !                  R(1) = R(Br-H)
            !                  R(2) = R(H-H)
            !                  R(3) = R(H-Br)
            !   The zero of energy is defined at Br "infinitely" far from the H2 dia
            !   The flags that indicate what calculations should be carried out in
            !   the potential routine are passed through the common block PT2CM:
            !   where:
            !        NASURF - which electronic states are available
            !        NDER  = 0 => no derivatives should be calculated
            !        NDER  = 1 => calculate first derivatives
            !        NFLAG  - these 20 integer values can be used to flag options
            !                within the potential;
            !   IPRT - FORTRAN unit number to which all output will be written
            !   Potential parameters' default settings
            !                  Variable            Default value
            !                  NASURF(1,1)         1
            !                  NDER                1
            !                  NFLAG(18)           6
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            real(wp), parameter :: cangau = 0.52917706e0_wp
            real(wp), parameter :: charkc = 627.5095e0_wp
            !         A1  =  V3CP1
            !         A2  =  V3CP2
            !         A3  =  V3CP3
            !         A4  =  V3CP4
            !         A5  =  V3CP5
            !         A6  =  V3CP6
            !   Convert the constants to atomic units
            do i = 1, 3
                de11(i)   = d11(i) / charkc
                re11(i)   = re11(i) / cangau
                beta11(i) = betaa11(i) * cangau
            enddo
            !    Initialize the asyptotic energy terms
            easybc = de11(2)
            easyac = de11(3)
            easyba = de11(1)
        end subroutine prepare_local_context_11

        !> @brief Evaluate POT11 for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine pot11
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a1, a2, a3, a4, a5, a6, b(3), capr, capr2, coschi, costh, dist, fcn4, pe, r132p2, r1m3, &
                & r1m3p2, r1p2, r1p3, r1p3m2, r2i, r2p2, r3p2, rr(3), schip2, sintp2, sintp6, sintp8, zchi
            integer             :: i
            integer, parameter  :: n3atom  = 75
            integer, parameter  :: isurf   = 5
            integer, parameter  :: jsurf   = isurf*(isurf+1)/2
            real(wp), parameter :: pi      = 3.141592653589793e0_wp
            integer, parameter  :: natom   = 25
            real(wp), parameter :: cangau  = 0.52917706e0_wp
            real(wp), parameter :: charkc  = 627.5095e0_wp
            real(wp), parameter :: xlamb11 = 0.4568e0_wp/27.21161e0_wp
            rr(1) = r(1)
            rr(2) = r(2)
            rr(3) = r(3)
            !   Compute bond dependent Sato parameters
            r1p2 = rr(1) * rr(1)
            r2p2 = rr(2) * rr(2)
            r3p2 = rr(3) * rr(3)
            r2i  = 1.0e0_wp/rr(2)
            zchi = 0.0e0_wp
            !   CAPR IS THE Br TO H2 DISTANCE
            capr2 = 0.5e0_wp * (r1p2 + r3p2 - 0.5e0_wp * r2p2)
            if( capr2 .le. 0.e0_wp )then
                capr2 = 0.0e0_wp
                capr  = 0.0e0_wp
            else
                capr = sqrt(capr2)
                !   Compute Chi, the angle between CAPR and RHH
                coschi = 0.5e0_wp*(r3p2 - r1p2)/(rr(2)*capr)
                schip2 = 1.0e0_wp - coschi*coschi
                zchi   = capr2*(zhbrp1*schip2 + zhbrp2*schip2*schip2)
            endif
            z(1) = 0.1675e0_wp + zchi
            z(3) = z(1)
            z(2) = zhhp1*r2p2/(zhhp2+r2p2)
            !   Compute the LEPS constants
            do i = 1, 3
                zpo(i)   = 1.0e0_wp+z(i)
                op3z(i)  = 1.0e0_wp+3.0e0_wp*z(i)
                top3z(i) = 2.0e0_wp*op3z(i)
                zp3(i)   = z(i)+3.0e0_wp
                tzp3(i)  = 2.0e0_wp*zp3(i)
                do4z(i)  = de11(i)/4.0e0_wp/zpo(i)
                b(i)     = beta11(i)*do4z(i)*2.0e0_wp
            enddo
            pe = 0.e0_wp
            do i = 1, 3
                legacy_x(i) = exp(-beta11(i)*(rr(i)-re11(i)))
                coul(i)     = do4z(i) * (zp3(i)*legacy_x(i) - top3z(i)) * legacy_x(i)
                exch(i)     = do4z(i) * (op3z(i)*legacy_x(i) - tzp3(i)) * legacy_x(i)
                pe          = pe + coul(i)
            enddo
            rad = sqrt((exch(1) - exch(2)) ** 2 + (exch(2) - exch(3)) ** 2 + (exch(3) - exch(1)) ** 2 )
            pe  = pe - (rad / r211) + ezero(1)
            !   Compute the three-center term
            r1p3   = rr(1) + rr(3)
            a1     = v3cp1
            a2     = v3cp2
            a3     = v3cp3
            a4     = v3cp4
            a5     = v3cp5
            a6     = v3cp6
            r1m3   = rr(1) - rr(3)
            r1p3m2 = r1p3 - rr(2)
            r1m3p2 = r1m3*r1m3
            r132p2 = r1p3m2*r1p3m2
            !   Compute the bond angle from the three internal coordinates
            dist   = rr(1)*rr(1) + rr(3)*rr(3) - rr(2)*rr(2)
            costh  = dist/(2.0e0_wp*rr(1)*rr(3))
            sintp2 = 1.0e0_wp - costh*costh
            sintp6 = sintp2*sintp2*sintp2
            sintp8 = sintp6*sintp2
            fcn4   = 1.e0_wp + a5*sintp2 + a6*sintp8
            v3c    = a1*fcn4*exp(-(a2*r1p3+a3*r1m3p2+a4*r132p2))
            !   Sum all the terms that make up the energy
            pe     = pe + v3c
            engygs = pe+xlamb11/3.e0_wp
        end subroutine pot11

        !> @brief Evaluate prepare_local_context_12 for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine prepare_local_context_12
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            real(wp), parameter :: cangau = 0.52917706e0_wp
            real(wp), parameter :: charkc = 627.5095e0_wp
        end subroutine prepare_local_context_12

        !> @brief Evaluate POT12 for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine pot12
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter  :: n3atom  = 75
            integer, parameter  :: isurf   = 5
            integer, parameter  :: jsurf   = isurf*(isurf+1)/2
            real(wp), parameter :: pi      = 3.141592653589793e0_wp
            integer, parameter  :: natom   = 25
            real(wp), parameter :: cangau  = 0.52917706e0_wp
            real(wp), parameter :: charkc  = 627.5095e0_wp
            real(wp), parameter :: xlamb12 = 0.4568e0_wp/27.21161e0_wp
            engyij(1) = -sqrt(2.e0_wp)*xlamb12/3.e0_wp
        end subroutine pot12

        !> @brief Evaluate prepare_local_context_22 for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine prepare_local_context_22
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            real(wp), parameter :: cangau = 0.52917706e0_wp
            real(wp), parameter :: charkc = 627.5095e0_wp
            !   Convert the constants to atomic units
            do i = 1, 3
                de22(i)   = d22(i) / charkc
                re22(i)   = re22(i) / cangau
                beta22(i) = betaa22(i) * cangau
            enddo
            !    Initialize the asyptotic energy terms
            easybc = de22(2)
            easyac = de22(3)
            easyba = de22(1)
        end subroutine prepare_local_context_22

        !> @brief Evaluate POT22 for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine pot22
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: pe, rr(3), vanti1, vanti3, vmorse, y1, y2, y3
            integer, parameter  :: n3atom  = 75
            integer, parameter  :: isurf   = 5
            integer, parameter  :: jsurf   = isurf*(isurf+1)/2
            real(wp), parameter :: pi      = 3.141592653589793e0_wp
            integer, parameter  :: natom   = 25
            real(wp), parameter :: cangau  = 0.52917706e0_wp
            real(wp), parameter :: charkc  = 627.5095e0_wp
            real(wp), parameter :: xlamb22 = 0.4568e0_wp/27.21161e0_wp
            rr(1)     = r(1)
            rr(2)     = r(2)
            rr(3)     = r(3)
            y2        = exp(-beta22(2)*(rr(2)-re22(2)))
            y1        = exp(-beta22(1)*(rr(1)-re22(1)))
            y3        = exp(-beta22(3)*(rr(3)-re22(3)))
            vmorse    = de22(2)*(y2**2-2.e0_wp*y2)
            vanti1    = de22(1)*(y1**2+2.e0_wp*y1)
            vanti3    = de22(3)*(y3**2+2.e0_wp*y3)
            pe        = ezero(2)  + vmorse+(vanti1+vanti3)*ccon+2.e0_wp*xlamb22/3.e0_wp
            engyes(1) = pe
        end subroutine pot22

        !> @brief Evaluate ANCVRT for the BrH2_LEPS_LTBZSO_1995 potential.
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
            do i=1, legacy_natoms
                isave(i) = 0
                jsave(i) = 0
                name1(i) = '  '
                name2(i) = '  '
            enddo
            iblank = ' '
            do ind=1, legacy_natoms
                do i=1, 7
                    do j=1, 32
                        if( indexes(ind) .eq. ianum(i,j) )then
                            isave(ind) = i
                            jsave(ind) = j
                        endif
                    enddo
                enddo
            enddo
            do ind=1, legacy_natoms
                ind2 = nulbl(ind)
                if( ind2 .eq. 0 ) ind2=ind
            enddo
            inc1 = 0
            do ind=1, irctnt-1
                inc1        = inc1+1
                name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
            enddo
            inc2 = 0
            do ind=irctnt, legacy_natoms
                inc2        = inc2+1
                name2(inc2) = periodic_1(isave(ind),jsave(ind))(:2)
            enddo
            if( nflag(1) .eq. 2 ) distance = 'ANGSTROMS           '
            if( nflag(2) .eq. 2 )then
                units = 'MILLIHARTREE        '
            elseif( nflag(2) .eq. 3 )then
                units = 'EV                  '
            elseif( nflag(2) .eq. 4 )then
                units = 'KCAL PER MOLE       '
            elseif( nflag(2) .eq. 5 )then
                units = 'WAVENUMBERS         '
            elseif( nflag(2) .eq. 6 )then
                units = 'KILOJOULES PER MOLE '
            endif
            cnvrtd  = 1.e0_wp
            cnvrte  = 1.e0_wp
            cnvrtde = 1.e0_wp
            if( nflag(1) .eq. 2 ) cnvrtd = bohr_a
            if( nflag(2) .eq. 2 )then
                cnvrte = cnvrte*htomillh
            elseif( nflag(2) .eq. 3 )then
                cnvrte = cnvrte*htoev
            elseif( nflag(2) .eq. 4 )then
                cnvrte = cnvrte*htokcal
            elseif( nflag(2) .eq. 5 )then
                cnvrte = cnvrte*htowave
            elseif( nflag(2) .eq. 6 )then
                cnvrte = cnvrte*htokj
            endif
            cnvrtde = cnvrte/cnvrtd
            isum    = 0
            do inu=1, 25
                isum = isum + nulbl(inu)
            enddo
            ireorder = 0
            if( isum .ne. 0 ) ireorder = 1
        end subroutine ancvrt

        !> @brief Evaluate CARTOU for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            use chempotpy_kinds,              only: wp
            implicit none
            integer            :: i, j
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            call convert_cartesian_input(cart, nulbl, legacy_natoms, cnvrtd, ireorder == 1, cartnu)
        end subroutine cartou

        !> @brief Evaluate CARTTOR for the BrH2_LEPS_LTBZSO_1995 potential.
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
            call cartesian_to_potlib_coordinates(cartnu, legacy_natoms, icartr, r, coordinate_status)
            if( coordinate_status /= 0 )then
                write(nflag(18), '(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') icartr
                stop
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine eunitzero
            use chempotpy_kinds, only: wp
            implicit none
            integer            :: i, j
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            pengygs = engygs * cnvrte - anuzero
            if( ksdiag .ne. 0 )then
                do i=ksdiag, kediag
                    pengyes(i) = engyes(i) * cnvrte - anuzero
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do j=ksoffd, keoffd
                    pengyij(j) = engyij(j) * cnvrte
                enddo
            endif
        end subroutine eunitzero

        !> @brief Evaluate RTOCART for the BrH2_LEPS_LTBZSO_1995 potential.
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
            if( icartr .eq. 1 )then
                do i = 1, legacy_natoms
                    ind            = 3*i-2
                    dgscartnu(i,1) = degsdr(ind)
                    dgscartnu(i,2) = degsdr(ind+1)
                    dgscartnu(i,3) = degsdr(ind+2)
                    if( ksdiag .ne. 0 )then
                        do j = ksdiag, kediag
                            descartnu(i,1,j) = deesdr(ind,j)
                            descartnu(i,2,j) = deesdr(ind+1,j)
                            descartnu(i,3,j) = deesdr(ind+2,j)
                        enddo
                    endif
                    if( keoffd .ne. 0 )then
                        do k = ksoffd, keoffd
                            dijcartnu(i,1,k) = deijdr(ind,k)
                            dijcartnu(i,2,k) = deijdr(ind+1,k)
                            dijcartnu(i,3,k) = deijdr(ind+2,k)
                        enddo
                    endif
                enddo
            elseif( icartr .eq. 2 )then
                do i = 1, legacy_natoms
                    dgscartnu(i,1) = 0.e0_wp
                    dgscartnu(i,2) = 0.e0_wp
                    dgscartnu(i,3) = 0.e0_wp
                    if( ksdiag .ne. 0 )then
                        do j1=ksdiag, kediag
                            descartnu(i,1,j1) = 0.e0_wp
                            descartnu(i,2,j1) = 0.e0_wp
                            descartnu(i,3,j1) = 0.e0_wp
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do j2=ksoffd, keoffd
                            dijcartnu(i,1,j2) = 0.e0_wp
                            dijcartnu(i,2,j2) = 0.e0_wp
                            dijcartnu(i,3,j2) = 0.e0_wp
                        enddo
                    endif
                    do j = 1, legacy_natoms
                        if( j .lt. i )then
                            m1 = legacy_natoms*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
                            m1 = legacy_natoms*(i-1) - (i*(i-1))/2 + j-i
                        else
                            cycle
                        endif
                        y              = degsdr(m1)
                        termx          = (cartnu(i,1)-cartnu(j,1))/r(m1)
                        termy          = (cartnu(i,2)-cartnu(j,2))/r(m1)
                        termz          = (cartnu(i,3)-cartnu(j,3))/r(m1)
                        dgscartnu(i,1) = dgscartnu(i,1) + termx*y
                        dgscartnu(i,2) = dgscartnu(i,2) + termy*y
                        dgscartnu(i,3) = dgscartnu(i,3) + termz*y
                        if( ksdiag .gt. 0 )then
                            y = deesdr(m1,j1)
                            do j1=ksdiag, kediag
                                descartnu(i,1,j1) = descartnu(i,1,j1) + termx*y
                                descartnu(i,2,j1) = descartnu(i,2,j1) + termy*y
                                descartnu(i,3,j1) = descartnu(i,3,j1) + termz*y
                            enddo
                        elseif( ksoffd .gt. 0 )then
                            do j2=ksoffd, keoffd
                                y                 = deijdr(m1,j2)
                                dijcartnu(i,1,j2) = dijcartnu(i,1,j2) + termx*y
                                dijcartnu(i,2,j2) = dijcartnu(i,2,j2) + termy*y
                                dijcartnu(i,3,j2) = dijcartnu(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            elseif( icartr .eq. 3 )then
                do i = 1, legacy_natoms
                    ygs(i) = degsdr(i)/r(i)
                    if( ksdiag .ne. 0 )then
                        do j=ksdiag, kediag
                            yes(i, j) = deesdr(i,j)/r(i)
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do k=ksoffd, keoffd
                            yij(i, k) = deijdr(i,k)/r(i)
                        enddo
                    endif
                enddo
                do k = 1, 3
                    term12         = cartnu(1,k)-cartnu(2,k)
                    term23         = cartnu(2,k)-cartnu(3,k)
                    term13         = cartnu(1,k)-cartnu(3,k)
                    dgscartnu(1,k) = term12*ygs(1) + term13*ygs(3)
                    dgscartnu(2,k) = -term12*ygs(1) + term23*ygs(2)
                    dgscartnu(3,k) = -term13*ygs(3) - term23*ygs(2)
                    if( ksdiag .ne. 0 )then
                        do j1=ksdiag, kediag
                            descartnu(1,k,j1) = term12*yes(1, j1) + term13*yes(3, j1)
                            descartnu(2,k,j1) = -term12*yes(1, j1) + term23*yes(2, j1)
                            descartnu(3,k,j1) = -term13*yes(3, j1) - term23*yes(2, j1)
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do j2=ksoffd, keoffd
                            dijcartnu(1,k,j2) = term12*yij(1, j2) + term13*yij(3, j2)
                            dijcartnu(2,k,j2) = -term12*yij(1, j2) + term23*yij(2, j2)
                            dijcartnu(3,k,j2) = -term13*yij(3, j2) - term23*yij(2, j2)
                        enddo
                    endif
                enddo
            elseif( icartr .eq. 4 )then
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
                write(nflag(18), '(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') icartr
                stop
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine dedcou
            use chempotpy_kinds, only: wp
            implicit none
            integer            :: i, j, k
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            if( ireorder .eq. 1 )then
                do i = 1, legacy_natoms
                    dgscart(i,1) = dgscartnu(nulbl(i),1) * cnvrtde
                    dgscart(i,2) = dgscartnu(nulbl(i),2) * cnvrtde
                    dgscart(i,3) = dgscartnu(nulbl(i),3) * cnvrtde
                    if( ksdiag .ne. 0 )then
                        do j=ksdiag, kediag
                            descart(i,1,j) = descartnu(nulbl(i),1,j) * cnvrtde
                            descart(i,2,j) = descartnu(nulbl(i),2,j) * cnvrtde
                            descart(i,3,j) = descartnu(nulbl(i),3,j) * cnvrtde
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do k=ksoffd, keoffd
                            dijcart(i,1,k) = dijcartnu(nulbl(i),1,k) * cnvrtde
                            dijcart(i,2,k) = dijcartnu(nulbl(i),2,k) * cnvrtde
                            dijcart(i,3,k) = dijcartnu(nulbl(i),3,k) * cnvrtde
                        enddo
                    endif
                enddo
            else
                do i = 1, legacy_natoms
                    dgscart(i,1) = dgscartnu(i,1) * cnvrtde
                    dgscart(i,2) = dgscartnu(i,2) * cnvrtde
                    dgscart(i,3) = dgscartnu(i,3) * cnvrtde
                    if( ksdiag .ne. 0 )then
                        do j=ksdiag, kediag
                            descart(i,1,j) = descartnu(i,1,j) * cnvrtde
                            descart(i,2,j) = descartnu(i,2,j) * cnvrtde
                            descart(i,3,j) = descartnu(i,3,j) * cnvrtde
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do k=ksoffd, keoffd
                            dijcart(i,1,k) = dijcartnu(i,1,k) * cnvrtde
                            dijcart(i,2,k) = dijcartnu(i,2,k) * cnvrtde
                            dijcart(i,3,k) = dijcartnu(i,3,k) * cnvrtde
                        enddo
                    endif
                enddo
            endif
        end subroutine dedcou
    end subroutine evaluate
end module chempotpy_brh2_leps_ltbzso_1995_pes

!> @brief Expose BrH2_LEPS_LTBZSO_1995 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_brh2_leps_ltbzso_1995_pes, only: evaluate
    use chempotpy_kinds,                     only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
