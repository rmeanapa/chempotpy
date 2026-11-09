!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_i2h_leps_kk_1981_pes
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
        double precision              :: anuzero
        double precision              :: b(3)
        double precision              :: cart(natom,3)
        double precision              :: cartnu(natom,3)
        double precision              :: cnvrtd
        double precision              :: cnvrtde
        double precision              :: cnvrte
        double precision              :: de(3)
        double precision              :: deesdr(n3atom,isurf)
        double precision              :: degsdr(n3atom)
        double precision              :: deijdr(n3atom,jsurf)
        double precision              :: descart(natom,3,isurf)
        double precision              :: descartnu(natom,3,isurf)
        double precision              :: dgscart(natom,3)
        double precision              :: dgscartnu(natom,3)
        double precision              :: dijcart(natom,3,jsurf)
        double precision              :: dijcartnu(natom,3,jsurf)
        double precision              :: easyab
        double precision              :: easyac
        double precision              :: easybc
        double precision              :: engyes(isurf)
        double precision              :: engygs
        double precision              :: engyij(jsurf)
        double precision              :: ezero(isurf+1)
        double precision              :: hde(3)
        double precision              :: hpf(3)
        integer                       :: icartr
        integer                       :: indexes(natom)
        integer                       :: irctnt
        integer                       :: ireorder
        integer                       :: kediag
        integer                       :: keoffd
        integer                       :: ksdiag
        integer                       :: ksoffd
        integer                       :: legacy_natoms
        integer                       :: mder
        integer                       :: msurf
        integer                       :: nasurf(isurf+1,isurf+1)
        integer                       :: nder
        integer                       :: nflag(20)
        integer                       :: nulbl(natom)
        double precision              :: pengyes(isurf)
        double precision              :: pengygs
        double precision              :: pengyij(jsurf)
        double precision              :: pf
        double precision              :: pfs(3)
        double precision              :: pft(3)
        double precision              :: r(n3atom)
        double precision              :: re(3)
        character                     :: ref(5)
        double precision              :: s(3)
        double precision              :: v
        integer                       :: iatom, idir, istate
        !initialize
        anuzero = 0.d0; b = 0.d0; cart = 0.d0; cartnu = 0.d0; cnvrtd = 0.d0; cnvrtde = 0.d0; cnvrte = 0.d0; de = 0.d0
        deesdr  = 0.d0; degsdr = 0.d0; deijdr = 0.d0; descart = 0.d0; descartnu = 0.d0; dgscart = 0.d0; dgscartnu = 0.d0
        dijcart = 0.d0; dijcartnu = 0.d0; easyab = 0.d0; easyac = 0.d0; easybc = 0.d0; engyes = 0.d0; engygs = 0.d0
        engyij  = 0.d0; ezero = 0.d0; hde = 0.d0; hpf = 0.d0; icartr = 0; indexes = 0; irctnt = 0; ireorder = 0
        kediag  = 0; keoffd = 0; ksdiag = 0; ksoffd = 0; legacy_natoms = 0; mder = 0; msurf = 0; nasurf = 0; nder = 0
        nflag   = 0; nulbl = 0; pengyes = 0.d0; pengygs = 0.d0; pengyij = 0.d0; pf = 0.d0; pfs = 0.d0; pft = 0.d0
        r       = 0.d0; re = 0.d0
        ref     = ' '
        s       = 0.d0
        de      = [3.3303d0, 3.3303d0, 1.5567d0]
        re      = [3.0236d0, 3.0236d0, 5.0457d0]
        b       = [0.9260d0, 0.9260d0, 0.9843d0]
        s       = [0.2000d0, 0.2000d0, 0.1250d0]
        nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0], shape(nasurf))
        nder          = 0
        nflag         = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        anuzero       = 0.0d0
        icartr        = 3
        msurf         = 0
        mder          = 1
        nulbl         = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        legacy_natoms = 3
        v             = 0.d0; g = 0.d0; d = 0.d0
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
            do istate=1, nstates
                p(istate) = pengygs*legacy_hartree_to_electronvolt
            enddo
        else if( igrad == 1 )then
            do istate=1, nstates
                p(istate) = pengygs*legacy_hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = dgscart(iatom,idir)*51.422067
                enddo
            enddo
        else if( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
        contains
        !> @brief Evaluate prepare_local_context for the I2H_LEPS_KK_1981 potential.
        subroutine prepare_local_context
            !
            ! System:          HI2
            ! Functional form: eLEPS (extended London-Eyring-Polanyi-Sato)
            ! Common name:     KK
            ! Reference:       J. A. Kaye and A. Kuppermann, Chem. Phys. Lett.
            !                  77, 573 (1981).
            !
            ! Number of bodies: 3
            ! Interface: potlib2001
            ! Number of electronic states: 1
            ! Number of derivatives: 1
            !
            ! prepare_local_context must be called once before any calls to POT are made.  It shoul
            ! never be called more than once.
            ! The potential parameters are included in DATA statements.
            ! Coordinates, potential energy, derivatives, and other information
            ! is passed through three common blocks:
            !           COMMON /PT1CM/ R(3), ENGYGS, DEGSDR(3)
            !           COMMON /ASYCM/ EASYAB, EASYBC, EASYAC
            ! The common block PT1CM contains the internal coordinates, the potentia
            ! energy, and the derivatives with respect to internal coordinates.  The
            ! common block PT2CM contains the a several variables which control
            ! various aspects of the calculation as outlined below.  The last common
            ! block, ASYCM, contains information about the differences in the
            ! dissociation energies which is useful for certain calculations.
            !
            ! Internuclear distances should be expressed in bohr and energies, etc.
            ! returned in Hartree atomic units.
            !
            ! The coordinates are defined as follows:
            !           R(1) = R(I-H)
            !           R(2) = R(H-I')
            !           R(3) = R(I-I')
            !
            ! The zero of energy for this potential occurs when the I atom is
            ! "infinitely" far from the HI molecule, and the HI molecule is at
            ! its equilibrium separation.
            !
            ! The array DEGSDR contains the analytic first derivatives of the potent
            ! energy with respect to the internal coordinates if the variable NDER i
            ! equal to 1 (see not below).  The derivatives are defined as follows:
            !           DEGSDR(1) = dV(R)/dR(1)
            !           DEGSDR(2) = dV(R)/dR(2)
            !           DEGSDR(3) = dV(R)/dR(3)
            ! where V(R) = V { R(1), R(2), R(3) }.
            !
            ! The variable NDER controls whether derivatives are calculated.
            ! If NDER = 1 first derivatives are calculated, otherwise no derivatives
            ! are calculated.  Note:  the default value of NDER is 1.
            !
            ! The array NFLAG is unused in this potential surface.
            !
            ! The common block ASYCM is defined as follows:
            !           EASYAB = DE(I-H)
            !           EASYBC = DE(H-I')
            !           EASYAC = DE(I-I')
            ! where DE is the dissociation energy in atomic units.
            !
            !   Potential parameters' default settings
            !                  Variable            Default value
            !                  NDER                1
            !                  NFLAG(18)           6
            !
            implicit none
            integer :: i
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
            double precision, parameter :: cevau = 27.2113961d0
            !
            !
            !
            !      save de,re,b,s,cevau
            !
            ! Write header potential parameters to UNIT NFLAG(18).
            !
            if( legacy_natoms .gt. 25 )then
                write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !
            ! Calculate first derivatives (default).  If no first derivatives are
            ! desired, nder should be set equal to zero after prepare_local_context is called.
            !
            !      nder=1
            !
            ! Convert to atomic units and calculate some useful constants.
            !
            do i=1, 3
                de(i)  = de(i)/cevau
                pf     = 0.5d0*de(i)*((1.d0-s(i))/(1.d0+s(i)))
                hde(i) = 0.5d0*de(i)
                hpf(i) = 0.5d0*pf
                pfs(i) = -de(i)*b(i)
                pft(i) = -pf*b(i)
            enddo
            easyab = de(1)
            easybc = de(2)
            easyac = de(3)
            !
            ezero(1) = de(2)
            !
            do i=1, 5
                ref(i) = ' '
            enddo
            !
            ref(1) = 'J. A. Kaye and A. Kuppermann,'
            ref(2) = 'Chem. Phys. Lett. 77, 573(1981)'
            !
            indexes(1) = 53
            indexes(2) = 1
            indexes(3) = 53
            !
            !
            !
            irctnt = 2
            !
            !CALL POTINFO
            !
            call ancvrt
            !
        end subroutine prepare_local_context
        !> @brief Evaluate POT for the I2H_LEPS_KK_1981 potential.
        subroutine pot
            !
            !
            ! The coordinates are defined as follows:
            !           R(1) = R(I-H)
            !           R(2) = R(H-I')
            !           R(3) = R(I-I')
            !
            ! The zero of energy for this potential occurs when the I atom is
            ! "infinitely" far from the HI molecule, and the HI molecule is at
            ! its equilibrium separation.
            !
            implicit none
            double precision :: dj1dr1, dj2dr2, dj3dr3, dq1dr1, dq2dr2, dq3dr3, hds1dr1, hds2dr2, hds3dr3, hdt1dr1, hdt2dr2, hdt3dr3, hradi, hs1, hs2, hs3, ht1, ht2, ht3, rad, radi, wj1, wj2, wj3, wq1, wq2, wq3, y1, y2, y3
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
            double precision, parameter :: cevau = 27.2113961d0
            !
            !
            !
            !      save de,re,b,s,cevau
            !
            call cartou
            call carttor
            !
            ! compute coulomb and exchange terms
            !
            y1  = dexp(b(1)*(re(1)-r(1)))
            hs1 = hde(1)*(y1-2.d0)*y1
            ht1 = hpf(1)*(y1+2.d0)*y1
            wq1 = hs1+ht1
            wj1 = hs1-ht1
            y2  = dexp(b(2)*(re(2)-r(2)))
            hs2 = hde(2)*(y2-2.d0)*y2
            ht2 = hpf(2)*(y2+2.d0)*y2
            wq2 = hs2+ht2
            wj2 = hs2-ht2
            y3  = dexp(b(3)*(re(3)-r(3)))
            hs3 = hde(3)*(y3-2.d0)*y3
            ht3 = hpf(3)*(y3+2.d0)*y3
            wq3 = hs3+ht3
            wj3 = hs3-ht3
            !
            ! compute potential energy
            !
            rad    = dsqrt(wj1*(wj1-wj2)+wj2*(wj2-wj3)+wj3*(wj3-wj1))
            engygs = ezero(1)  +wq1+wq2+wq3-rad
            !
            !  compute first derivatives if necessary
            !
            if( nder .eq. 1 )then
                radi      = 1.d0/rad
                hradi     = 0.5d0*radi
                hds1dr1   = pfs(1)*(y1-1.d0)*y1
                hds2dr2   = pfs(2)*(y2-1.d0)*y2
                hds3dr3   = pfs(3)*(y3-1.d0)*y3
                hdt1dr1   = pft(1)*(y1+1.d0)*y1
                hdt2dr2   = pft(2)*(y2+1.d0)*y2
                hdt3dr3   = pft(3)*(y3+1.d0)*y3
                dq1dr1    = hds1dr1+hdt1dr1
                dq2dr2    = hds2dr2+hdt2dr2
                dq3dr3    = hds3dr3+hdt3dr3
                dj1dr1    = hds1dr1-hdt1dr1
                dj2dr2    = hds2dr2-hdt2dr2
                dj3dr3    = hds3dr3-hdt3dr3
                degsdr(1) = dq1dr1-hradi*dj1dr1*(2.d0*wj1-wj2-wj3)
                degsdr(2) = dq2dr2-hradi*dj2dr2*(2.d0*wj2-wj1-wj3)
                degsdr(3) = dq3dr3-hradi*dj3dr3*(2.d0*wj3-wj1-wj2)
            endif
            call eunitzero
            if( nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
            !
        end subroutine pot
        !> @brief Evaluate ANCVRT for the I2H_LEPS_KK_1981 potential.
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
            else if( nflag(2) .eq. 3 )then
                units = 'EV                  '
            else if( nflag(2) .eq. 4 )then
                units = 'KCAL PER MOLE       '
            else if( nflag(2) .eq. 5 )then
                units = 'WAVENUMBERS         '
            else if( nflag(2) .eq. 6 )then
                units = 'KILOJOULES PER MOLE '
            endif
            cnvrtd  = 1.d0
            cnvrte  = 1.d0
            cnvrtde = 1.d0
            if( nflag(1) .eq. 2 ) cnvrtd = bohr_a
            if( nflag(2) .eq. 2 )then
                cnvrte = cnvrte*htomillh
            else if( nflag(2) .eq. 3 )then
                cnvrte = cnvrte*htoev
            else if( nflag(2) .eq. 4 )then
                cnvrte = cnvrte*htokcal
            else if( nflag(2) .eq. 5 )then
                cnvrte = cnvrte*htowave
            else if( nflag(2) .eq. 6 )then
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
        !> @brief Evaluate CARTOU for the I2H_LEPS_KK_1981 potential.
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
        !> @brief Evaluate CARTTOR for the I2H_LEPS_KK_1981 potential.
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
        !> @brief Evaluate EUNITZERO for the I2H_LEPS_KK_1981 potential.
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
        !> @brief Evaluate RTOCART for the I2H_LEPS_KK_1981 potential.
        subroutine rtocart
            implicit none
            double precision   :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes, ygs, yij, za, zb, zc, zd
            integer            :: i, ind, j, j1, j2, k, m1
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            dimension ygs(n3atom), yes(n3atom,isurf), yij(n3atom,jsurf)
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
            else if( icartr .eq. 2 )then
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
                        else if( j .gt. i )then
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
                        else if( ksoffd .gt. 0 )then
                            do j2=ksoffd, keoffd
                                y                 = deijdr(m1,j2)
                                dijcartnu(i,1,j2) = dijcartnu(i,1,j2) + termx*y
                                dijcartnu(i,2,j2) = dijcartnu(i,2,j2) + termy*y
                                dijcartnu(i,3,j2) = dijcartnu(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            else if( icartr .eq. 3 )then
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
            else if( icartr .eq. 4 )then
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
        !> @brief Evaluate DEDCOU for the I2H_LEPS_KK_1981 potential.
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
end module chempotpy_i2h_leps_kk_1981_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_i2h_leps_kk_1981_pes, only: evaluate
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
