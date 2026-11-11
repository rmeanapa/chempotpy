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
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the BrH2_LEPS_LTBZSO_1995 potential energy surface for BrH2.
    subroutine evaluate( x, igrad, p, g, d )
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
        double precision              :: anuzero
        double precision              :: beta11(3)
        double precision              :: beta22(3)
        double precision              :: betaa11(3)
        double precision              :: betaa22(3)
        double precision              :: cart(natom,3)
        double precision              :: cartnu(natom,3)
        double precision              :: ccon
        double precision              :: cnvrtd
        double precision              :: cnvrtde
        double precision              :: cnvrte
        double precision              :: coul(3)
        double precision              :: d11(3)
        double precision              :: d22(3)
        double precision              :: de11(3)
        double precision              :: de22(3)
        double precision              :: deesdr(n3atom,isurf)
        double precision              :: degsdr(n3atom)
        double precision              :: deijdr(n3atom,jsurf)
        double precision              :: deq(3)
        double precision              :: derad(3)
        double precision              :: descart(natom,3,isurf)
        double precision              :: descartnu(natom,3,isurf)
        double precision              :: dgscart(natom,3)
        double precision              :: dgscartnu(natom,3)
        double precision              :: dijcart(natom,3,jsurf)
        double precision              :: dijcartnu(natom,3,jsurf)
        double precision              :: do4z(3)
        double precision              :: dv3c(3)
        double precision              :: dxdr(3,3)
        double precision              :: dz(3,3)
        double precision              :: easyac
        double precision              :: easyba
        double precision              :: easybc
        double precision              :: engyes(isurf)
        double precision              :: engygs
        double precision              :: engyij(jsurf)
        double precision              :: exch(3)
        double precision              :: ezero(isurf+1)
        integer                       :: icartr
        integer                       :: indexes(natom)
        integer                       :: irctnt
        integer                       :: ireorder
        integer                       :: kediag
        integer                       :: keoffd
        integer                       :: ksdiag
        integer                       :: ksoffd
        integer                       :: legacy_natoms
        double precision              :: legacy_x(3)
        integer                       :: mder
        integer                       :: msurf
        integer                       :: nasurf(isurf+1,isurf+1)
        integer                       :: nder
        integer                       :: nflag(20)
        integer                       :: nulbl(natom)
        double precision              :: op3z(3)
        double precision              :: pengyes(isurf)
        double precision              :: pengygs
        double precision              :: pengyij(jsurf)
        double precision              :: r(n3atom)
        double precision              :: r211
        double precision              :: r222
        double precision              :: rad
        double precision              :: re11(3)
        double precision              :: re22(3)
        character                     :: ref(5)
        double precision              :: top3z(3)
        double precision              :: tzp3(3)
        double precision              :: v3c
        double precision              :: v3cp1
        double precision              :: v3cp2
        double precision              :: v3cp3
        double precision              :: v3cp4
        double precision              :: v3cp5
        double precision              :: v3cp6
        double precision              :: z(3)
        double precision              :: zhbrp1
        double precision              :: zhbrp2
        double precision              :: zhhp1
        double precision              :: zhhp2
        double precision              :: zp3(3)
        double precision              :: zpo(3)
        integer                       :: iatom, idir
        double precision              :: u(nstates,nstates)
        double precision              :: t(nstates,nstates)
        !initialize
        anuzero   = 0.d0; beta11 = 0.d0; beta22 = 0.d0; betaa11 = 0.d0; betaa22 = 0.d0; cart = 0.d0; cartnu = 0.d0
        ccon      = 0.d0; cnvrtd = 0.d0; cnvrtde = 0.d0; cnvrte = 0.d0; coul = 0.d0; d11 = 0.d0; d22 = 0.d0; de11 = 0.d0
        de22      = 0.d0; deesdr = 0.d0; degsdr = 0.d0; deijdr = 0.d0; deq = 0.d0; derad = 0.d0; descart = 0.d0
        descartnu = 0.d0; dgscart = 0.d0; dgscartnu = 0.d0; dijcart = 0.d0; dijcartnu = 0.d0; do4z = 0.d0; dv3c = 0.d0
        dxdr      = 0.d0; dz = 0.d0; easyac = 0.d0; easyba = 0.d0; easybc = 0.d0; engyes = 0.d0; engygs = 0.d0; engyij = 0.d0
        exch      = 0.d0; ezero = 0.d0; icartr = 0; indexes = 0; irctnt = 0; ireorder = 0; kediag = 0; keoffd = 0; ksdiag = 0
        ksoffd    = 0; legacy_natoms = 0; legacy_x = 0.d0; mder = 0; msurf = 0; nasurf = 0; nder = 0; nflag = 0; nulbl = 0
        op3z      = 0.d0; pengyes = 0.d0; pengygs = 0.d0; pengyij = 0.d0; r = 0.d0; r211 = 0.d0; r222 = 0.d0; rad = 0.d0
        re11      = 0.d0; re22 = 0.d0
        ref       = ' '
        top3z     = 0.d0
        tzp3      = 0.d0
        v3c       = 0.d0
        v3cp1     = 0.d0
        v3cp2     = 0.d0
        v3cp3     = 0.d0
        v3cp4     = 0.d0
        v3cp5     = 0.d0
        v3cp6     = 0.d0
        z         = 0.d0
        zhbrp1    = 0.d0
        zhbrp2    = 0.d0
        zhhp1     = 0.d0
        zhhp2     = 0.d0
        zp3       = 0.d0
        zpo       = 0.d0
        nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0], shape(nasurf))
        nder          = 0
        nflag         = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        anuzero       = 0.0d0
        icartr        = 3
        msurf         = 1
        mder          = 0
        nulbl         = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        legacy_natoms = 3
        d11           = [90.4473d0, 109.559d0, 90.4473d0]
        re11          = [1.41443d0, 0.74144d0, 1.41443d0]
        betaa11       = [1.81088d0, 1.94198d0, 1.81088d0]
        r211          = 1.41421356d0
        d22           = [90.4473d0, 109.559d0, 90.4473d0]
        re22          = [1.41443d0, 0.74144d0, 1.41443d0]
        betaa22       = [1.81088d0, 1.94198d0, 1.81088d0]
        r222          = 1.41421356d0
        ccon          = 0.5d0
        zhhp1         = 2.91230499d-1
        zhhp2         = 4.99160974d0
        zhbrp1        = 3.84262775d-3
        zhbrp2        = 2.87642885d-3
        v3cp1         = 1.67457884d0
        v3cp2         = 7.09518875d-1
        v3cp3         = 1.87582679d-1
        v3cp4         = 2.43517074d-1
        v3cp5         = 2.04145206d-1
        v3cp6         = 2.57833022d0
        u             = 0.d0; g = 0.d0; d = 0.d0
        nasurf        = 1
        cart          = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                cart(iatom,idir) = x(iatom,idir)/legacy_bohr_to_angstrom
            enddo
        enddo
        nder = igrad
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            u(1,1) = engygs*legacy_hartree_to_electronvolt
            u(2,2) = engyes(1)*legacy_hartree_to_electronvolt
            u(1,2) = engyij(1)*legacy_hartree_to_electronvolt
            u(2,1) = engyij(1)*legacy_hartree_to_electronvolt
            call diagonalize(nstates, u, p, t)
        else
            write (*, *) 'Only energy is available'
        endif

        contains

        !> @brief Evaluate diagonalize for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine diagonalize( n, a_ss, ev_s, u_ss )
            implicit none
            integer, intent(in)           :: n
            double precision, intent(in)  :: a_ss(n,n)
            double precision, intent(out) :: u_ss(n,n)
            double precision, intent(out) :: ev_s(n)
            integer                       :: io, i
            double precision, allocatable :: work_d(:)
            integer                       :: lwork_d
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
            implicit none
            integer                     :: i
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
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
            !CALL POTINFO
            call ancvrt
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine pot
            implicit none
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
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
            !                  COMMON /ASYCM/ EASYBC, EASYAC, EASYBA
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
            implicit none
            integer                     :: i
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            double precision, parameter :: cangau = 0.52917706d0
            double precision, parameter :: charkc = 627.5095d0
            !         A1  =  V3CP1
            !         A2  =  V3CP2
            !         A3  =  V3CP3
            !         A4  =  V3CP4
            !         A5  =  V3CP5
            !         A6  =  V3CP6
            !         WRITE (NFLAG(18),100) D11, RE11, BETAA11,
            !     +                    ZHHP1, ZHHP2, ZHBRP1, ZHBRP2,
            !     +                    V3CP1, V3CP2, V3CP3, V3CP4, V3CP5
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
            implicit none
            double precision            :: a1, a2, a3, a4, a5, a6, b(3), capr, capr2, coschi, costh, dist, fcn4, pe, r132p2, r1m3, &
                & r1m3p2, r1p2, r1p3, r1p3m2, r2i, r2p2, r3p2, rr(3), schip2, sintp2, sintp6, sintp8, zchi
            integer                     :: i
            integer, parameter          :: n3atom  = 75
            integer, parameter          :: isurf   = 5
            integer, parameter          :: jsurf   = isurf*(isurf+1)/2
            double precision, parameter :: pi      = 3.141592653589793d0
            integer, parameter          :: natom   = 25
            double precision, parameter :: cangau  = 0.52917706d0
            double precision, parameter :: charkc  = 627.5095d0
            double precision, parameter :: xlamb11 = 0.4568d0/27.21161d0
            rr(1) = r(1)
            rr(2) = r(2)
            rr(3) = r(3)
            !   Compute bond dependent Sato parameters
            r1p2 = rr(1) * rr(1)
            r2p2 = rr(2) * rr(2)
            r3p2 = rr(3) * rr(3)
            r2i  = 1.0d0/rr(2)
            zchi = 0.0d0
            !   CAPR IS THE Br TO H2 DISTANCE
            capr2 = 0.5d0 * (r1p2 + r3p2 - 0.5d0 * r2p2)
            if( capr2 .le. 0.d0 )then
                capr2 = 0.0d0
                capr  = 0.0d0
            else
                capr = sqrt(capr2)
                !   Compute Chi, the angle between CAPR and RHH
                coschi = 0.5d0*(r3p2 - r1p2)/(rr(2)*capr)
                schip2 = 1.0d0 - coschi*coschi
                zchi   = capr2*(zhbrp1*schip2 + zhbrp2*schip2*schip2)
            endif
            z(1) = 0.1675d0 + zchi
            z(3) = z(1)
            z(2) = zhhp1*r2p2/(zhhp2+r2p2)
            !   Compute the LEPS constants
            do i = 1, 3
                zpo(i)   = 1.0d0+z(i)
                op3z(i)  = 1.0d0+3.0d0*z(i)
                top3z(i) = 2.0d0*op3z(i)
                zp3(i)   = z(i)+3.0d0
                tzp3(i)  = 2.0d0*zp3(i)
                do4z(i)  = de11(i)/4.0d0/zpo(i)
                b(i)     = beta11(i)*do4z(i)*2.0d0
            enddo
            pe = 0.d0
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
            costh  = dist/(2.0d0*rr(1)*rr(3))
            sintp2 = 1.0d0 - costh*costh
            sintp6 = sintp2*sintp2*sintp2
            sintp8 = sintp6*sintp2
            fcn4   = 1.d0 + a5*sintp2 + a6*sintp8
            v3c    = a1*fcn4*exp(-(a2*r1p3+a3*r1m3p2+a4*r132p2))
            !   Sum all the terms that make up the energy
            pe     = pe + v3c
            engygs = pe+xlamb11/3.d0
        end subroutine pot11

        !> @brief Evaluate prepare_local_context_12 for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine prepare_local_context_12
            implicit none
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            double precision, parameter :: cangau = 0.52917706d0
            double precision, parameter :: charkc = 627.5095d0
        end subroutine prepare_local_context_12

        !> @brief Evaluate POT12 for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine pot12
            implicit none
            integer, parameter          :: n3atom  = 75
            integer, parameter          :: isurf   = 5
            integer, parameter          :: jsurf   = isurf*(isurf+1)/2
            double precision, parameter :: pi      = 3.141592653589793d0
            integer, parameter          :: natom   = 25
            double precision, parameter :: cangau  = 0.52917706d0
            double precision, parameter :: charkc  = 627.5095d0
            double precision, parameter :: xlamb12 = 0.4568d0/27.21161d0
            engyij(1) = -dsqrt(2.d0)*xlamb12/3.d0
        end subroutine pot12

        !> @brief Evaluate prepare_local_context_22 for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine prepare_local_context_22
            implicit none
            integer                     :: i
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            double precision, parameter :: cangau = 0.52917706d0
            double precision, parameter :: charkc = 627.5095d0
            !         WRITE (NFLAG(18),100) D22, RE22, BETAA22,
            !     +                 ZHHP1, ZHHP2, ZHBRP1, ZHBRP2,
            !     +                 V3CP1, V3CP2, V3CP3, V3CP4, V3CP5
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
            implicit none
            double precision            :: pe, rr(3), vanti1, vanti3, vmorse, y1, y2, y3
            integer, parameter          :: n3atom  = 75
            integer, parameter          :: isurf   = 5
            integer, parameter          :: jsurf   = isurf*(isurf+1)/2
            double precision, parameter :: pi      = 3.141592653589793d0
            integer, parameter          :: natom   = 25
            double precision, parameter :: cangau  = 0.52917706d0
            double precision, parameter :: charkc  = 627.5095d0
            double precision, parameter :: xlamb22 = 0.4568d0/27.21161d0
            rr(1)     = r(1)
            rr(2)     = r(2)
            rr(3)     = r(3)
            y2        = exp(-beta22(2)*(rr(2)-re22(2)))
            y1        = exp(-beta22(1)*(rr(1)-re22(1)))
            y3        = exp(-beta22(3)*(rr(3)-re22(3)))
            vmorse    = de22(2)*(y2**2-2.d0*y2)
            vanti1    = de22(1)*(y1**2+2.d0*y1)
            vanti3    = de22(3)*(y3**2+2.d0*y3)
            pe        = ezero(2)  + vmorse+(vanti1+vanti3)*ccon+2.d0*xlamb22/3.d0
            engyes(1) = pe
        end subroutine pot22

        !> @brief Evaluate ANCVRT for the BrH2_LEPS_LTBZSO_1995 potential.
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
            cnvrtd  = 1.d0
            cnvrte  = 1.d0
            cnvrtde = 1.d0
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
            implicit none
            integer            :: coordinate_status
            double precision   :: cmm, cphi, flm, hhd, hym, phi, pi, q, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
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
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            double precision   :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, &
                & termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes(n3atom,isurf), &
                & ygs(n3atom), yij(n3atom,jsurf), za, zb, zc, zd
            integer            :: i, ind, j, j1, j2, k, m1
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
                    dgscartnu(i,1) = 0.d0
                    dgscartnu(i,2) = 0.d0
                    dgscartnu(i,3) = 0.d0
                    if( ksdiag .ne. 0 )then
                        do j1=ksdiag, kediag
                            descartnu(i,1,j1) = 0.d0
                            descartnu(i,2,j1) = 0.d0
                            descartnu(i,3,j1) = 0.d0
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do j2=ksoffd, keoffd
                            dijcartnu(i,1,j2) = 0.d0
                            dijcartnu(i,2,j2) = 0.d0
                            dijcartnu(i,3,j2) = 0.d0
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
                            yes(i,j) = deesdr(i,j)/r(i)
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do k=ksoffd, keoffd
                            yij(i,k) = deijdr(i,k)/r(i)
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
                            descartnu(1,k,j1) = term12*yes(1,j1) + term13*yes(3,j1)
                            descartnu(2,k,j1) = -term12*yes(1,j1) + term23*yes(2,j1)
                            descartnu(3,k,j1) = -term13*yes(3,j1) - term23*yes(2,j1)
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do j2=ksoffd, keoffd
                            dijcartnu(1,k,j2) = term12*yij(1,j2) + term13*yij(3,j2)
                            dijcartnu(2,k,j2) = -term12*yij(1,j2) + term23*yij(2,j2)
                            dijcartnu(3,k,j2) = -term13*yij(3,j2) - term23*yij(2,j2)
                        enddo
                    endif
                enddo
            elseif( icartr .eq. 4 )then
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
                write(nflag(18), '(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') icartr
                stop
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the BrH2_LEPS_LTBZSO_1995 potential.
        subroutine dedcou
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
