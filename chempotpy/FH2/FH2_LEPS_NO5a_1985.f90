!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_fh2_leps_no5a_1985_pes
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
        double precision              :: a
        double precision              :: a1
        double precision              :: a2
        double precision              :: alf
        double precision              :: alfp
        double precision              :: alft
        double precision              :: anuzero
        double precision              :: b
        double precision              :: bet
        double precision              :: beta(3)
        double precision              :: betp
        double precision              :: c
        double precision              :: c1
        double precision              :: c2
        double precision              :: c3
        double precision              :: c3c
        double precision              :: c4
        double precision              :: c5
        double precision              :: c6
        double precision              :: cart(natom,3)
        double precision              :: cartnu(natom,3)
        double precision              :: cnvrtd
        double precision              :: cnvrtde
        double precision              :: cnvrte
        double precision              :: deesdr(n3atom,isurf)
        double precision              :: deg1
        double precision              :: deg2
        double precision              :: deg3
        double precision              :: deg3p
        double precision              :: deg4
        double precision              :: degsdr(n3atom)
        double precision              :: deijdr(n3atom,jsurf)
        double precision              :: descart(natom,3,isurf)
        double precision              :: descartnu(natom,3,isurf)
        double precision              :: dgscart(natom,3)
        double precision              :: dgscartnu(natom,3)
        double precision              :: dijcart(natom,3,jsurf)
        double precision              :: dijcartnu(natom,3,jsurf)
        double precision              :: dj(3,3)
        double precision              :: dm(3)
        double precision              :: dq(3,3)
        double precision              :: dz(3,3)
        double precision              :: easyab
        double precision              :: easyac
        double precision              :: easybc
        double precision              :: engyes(isurf)
        double precision              :: engygs
        double precision              :: engyij(jsurf)
        double precision              :: ezero(isurf+1)
        double precision              :: f
        double precision              :: g2
        double precision              :: ha2
        integer                       :: icartr
        integer                       :: indexes(natom)
        integer                       :: irctnt
        integer                       :: ireorder
        integer                       :: kediag
        integer                       :: keoffd
        integer                       :: ksdiag
        integer                       :: ksoffd
        double precision              :: legacy_d
        double precision              :: legacy_g
        integer                       :: legacy_natoms
        integer                       :: mder
        integer                       :: msurf
        integer                       :: nasurf(isurf+1,isurf+1)
        integer                       :: nder
        integer                       :: nflag(20)
        integer                       :: nulbl(natom)
        double precision              :: omp3c
        double precision              :: omq3c
        double precision              :: p3c
        double precision              :: pengyes(isurf)
        double precision              :: pengygs
        double precision              :: pengyij(jsurf)
        double precision              :: q(3)
        double precision              :: q3c
        double precision              :: r(n3atom)
        double precision              :: re(3)
        character                     :: ref(5)
        double precision              :: rr(3)
        double precision              :: xj(3)
        double precision              :: z(3)
        double precision              :: zfh0
        double precision              :: zhh0
        double precision              :: zhhdif
        double precision              :: zhhst
        double precision              :: v
        integer                       :: iatom, idir, istate
        !initialize
        a         = 0.d0; a1 = 0.d0; a2 = 0.d0; alf = 0.d0; alfp = 0.d0; alft = 0.d0; anuzero = 0.d0; b = 0.d0; bet = 0.d0
        beta      = 0.d0; betp = 0.d0; c = 0.d0; c1 = 0.d0; c2 = 0.d0; c3 = 0.d0; c3c = 0.d0; c4 = 0.d0; c5 = 0.d0; c6 = 0.d0
        cart      = 0.d0; cartnu = 0.d0; cnvrtd = 0.d0; cnvrtde = 0.d0; cnvrte = 0.d0; deesdr = 0.d0; deg1 = 0.d0
        deg2      = 0.d0; deg3 = 0.d0; deg3p = 0.d0; deg4 = 0.d0; degsdr = 0.d0; deijdr = 0.d0; descart = 0.d0
        descartnu = 0.d0; dgscart = 0.d0; dgscartnu = 0.d0; dijcart = 0.d0; dijcartnu = 0.d0; dj = 0.d0; dm = 0.d0
        dq        = 0.d0; dz = 0.d0; easyab = 0.d0; easyac = 0.d0; easybc = 0.d0; engyes = 0.d0; engygs = 0.d0; engyij = 0.d0
        ezero     = 0.d0; f = 0.d0; g2 = 0.d0; ha2 = 0.d0; icartr = 0; indexes = 0; irctnt = 0; ireorder = 0; kediag = 0
        keoffd    = 0; ksdiag = 0; ksoffd = 0; legacy_d = 0.d0; legacy_g = 0.d0; legacy_natoms = 0; mder = 0; msurf = 0
        nasurf    = 0; nder = 0; nflag = 0; nulbl = 0; omp3c = 0.d0; omq3c = 0.d0; p3c = 0.d0; pengyes = 0.d0
        pengygs   = 0.d0; pengyij = 0.d0; q = 0.d0; q3c = 0.d0; r = 0.d0; re = 0.d0
        ref       = ' '
        rr        = 0.d0
        xj        = 0.d0
        z         = 0.d0
        zfh0      = 0.d0
        zhh0      = 0.d0
        zhhdif    = 0.d0
        zhhst     = 0.d0
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
        dm            = [141.196d0, 109.449d0, 141.196d0]
        re            = [0.9170d0, 0.7419d0, 0.9170d0]
        beta          = [2.2187d0, 1.9420d0, 2.2187d0]
        a1            = 0.0395d0
        a2            = 0.201d0
        a             = 1.26d0
        b             = 1.60d0
        c             = 0.0d0
        legacy_d      = 0.0d0
        f             = 0.0d0
        legacy_g      = 0.2d0
        bet           = -0.101168d0
        betp          = -0.46183d0
        zfh0          = 0.170d0
        zhh0          = 0.120d0
        zhhst         = 0.23107d0
        deg1          = 10.0d0
        deg2          = 20.0d0
        deg3          = 60.0d0
        deg4          = 90.0d0
        deg3p         = 41.0d0
        c1            = 1.15d-4
        c2            = -6.0d-6
        c3            = 0.016d0
        c4            = 0.0005d0
        c5            = 0.00069d0
        c6            = 3.308d-2
        c3c           = 0.270063d0
        alf           = 0.18d0
        alfp          = 0.078d0
        alft          = 2.14d0
        p3c           = 0.95d0
        q3c           = 0.48d0
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
        !> @brief Evaluate prepare_local_context for the FH2_LEPS_NO5a_1985 potential.
        subroutine prepare_local_context
            !
            !   System:          FH2
            !   Functional form: eLEPS (extended-London-Eyring-Polanyi-Sato)
            !                    plus E3C (a three-center term).
            !   Common name: no. 5A
            !   Reference:       R. Steckler, legacy_d. legacy_g. Truhlar, and B. C. Garrett
            !                    J. Chem. Phys. 82, 5499-5505 (1985).
            !
            !   Notes:           In this version of the potential the AB and AC
            !                    Sato parameters are functions of PHI, where PHI is
            !                    the angle of BC-to-A vector with the BC axis.
            !                    The BC Sato parameter is a function of BC bond
            !                    length and the angle PHI.  A three-center term is
            !                    also added to correct the H-F-H barrier height.
            !
            !   prepare_local_context must be called once before any calls to POT.
            !   The potential parameters are included in DATA statements
            !   and in the block data subprogram PTPACM.
            !   Coordinates, potential energy, and derivatives are passed
            !   The potential energy in the three asymptotic valleys are
            !   stored in the common block ASYCM:
            !                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the FH diatomic, with the
            !   FH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
            !   respectively.
            !   All the information passed through the common blocks PT1CM and ASYCM
            !   is in Hartree atomic units.
            !
            !   This potential is written such that:
            !                  R(1) = R(F-H)
            !                  R(2) = R(H-H)
            !                  R(3) = R(H-F)
            !   The zero of energy is defined at F "infinitely" far from the H2 diat
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
            double precision, parameter :: ckcau  = 627.5095d0
            double precision, parameter :: cangau = 0.52917706d0
            !
            !   EXLARG IS LN(LARGEST FLOATING POINT NUMBER ALLOWED)
            !
            !
            !      LOGICAL LCOL, LDERIV, LZERO
            !
            !
            if( legacy_natoms .gt. 25 )then
                write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !
            !
            !
            do i = 1, 3
                dm(i)   = dm(i) / ckcau
                re(i)   = re(i) / cangau
                beta(i) = beta(i) * cangau
            10 enddo
            !
            !   Set the energy of the potential in the each of the three asymptotic
            !
            easyab = dm(1)
            easybc = dm(2)
            easyac = dm(3)
            !
            !   USEFUL CONSTANTS
            !
            zhhdif = zhhst - zhh0
            omp3c  = 1.d0 - p3c
            omq3c  = 1.d0 - q3c
            ha2    = 0.5d0 * a2
            g2     = legacy_g * legacy_g
            !
            !
            ezero(1) = dm(2)
            !
            do i=1, 5
                ref(i) = ' '
            enddo
            !
            ref(1) = 'R. Steckler, legacy_d. legacy_g. Truhlar, B. C. Garrett,'
            ref(2) = 'J. Chem. Phys. 82, 5499(1985)'
            !
            indexes(1) = 9
            indexes(2) = 1
            indexes(3) = 1
            !
            !
            !
            irctnt = 2
            !
            !      CALL POTINFO
            !
            call ancvrt
            !
        end subroutine prepare_local_context
        !> @brief Evaluate POT for the FH2_LEPS_NO5a_1985 potential.
        subroutine pot
            !
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the FH diatomic, with the
            !   FH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
            !   respectively.
            !
            !   This potential is written such that:
            !                  R(1) = R(F-H)
            !                  R(2) = R(H-H)
            !                  R(3) = R(H-F)
            !   The zero of energy is defined at F "infinitely" far from the H2 diat
            !
            implicit none
            double precision :: capr, capr2, cosphi, costh, cphi, dcos1, dcos2, dcos3, dcphi, dd2, de3c1, de3c2, de3c3, deg, degdif, dphi1, dphi2, dphi3, dzphi, e3c, ex, ex1, ex2, ex3, exx, hsech, op3z, phi, r12, r1i, r22, r32, r3i, rdif, rdif2, rsum, scale1, scale2, sgnphi, sinp2, sinphi, t, t1, t2, t3, t4, t5, xjx, zhh, zp1, zp3, zphi
            integer          :: i, i1, i2, i3, j
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
            double precision, parameter :: ckcau  = 627.5095d0
            double precision, parameter :: cangau = 0.52917706d0
            double precision, parameter :: tdeg   = 57.29577951d0
            double precision, parameter :: exlarg = 80.d0
            !
            !
            logical :: lcol, lderiv, lzero
            !
            data lderiv /.false./
            !
            !
            call cartou
            call carttor
            !
            !   Check the values of NASURF and NDER for validity.
            !
            !   Check the value of NDER
            !
            !
            !   The logical LDERIV is used to flag the derivative calculations
            !   in the code; LDERIV is set equal to .TRUE. if NDER = 1.
            !
            if( nder .eq. 1 ) lderiv = .true.
            !
            !   TEST R'S FOR PHYSICAL VALUES
            !
            do i = 1, 3
                rr(i) = r(i)
            20 enddo
            do i1 = 1, 3
                i2                         = mod(i1, 3) + 1
                i3                         = mod(i2, 3) + 1
                t                          = rr(i2) + rr(i3)
                if( rr(i1) .gt. t ) rr(i1) = t
            30 enddo
            r12 = rr(1) * rr(1)
            r22 = rr(2) * rr(2)
            r32 = rr(3) * rr(3)
            !
            !   CAPR IS THE F TO H2 DISTANCE
            !
            capr2 = 0.5d0 * (r12 + r32 - 0.5d0 * r22)
            if( capr2 .lt. 0.d0 ) capr2 = 0.0d0
            capr = sqrt(capr2)
            !
            !   PHI IS THE ANGLE BETWEEN CAPR AND RHH
            !   PHI IS RESTRICTED TO BE BETWEEN O AND 90 DEGREES
            !
            t = 0.5d0 * (r32 - r12)
            if( abs(t) .gt. 1.d-7 ) t = t / (rr(2) * capr)
            t2     = 1.0d0
            sgnphi = sign(t2, t)
            cosphi = sgnphi * t
            cosphi = min(cosphi, t2)
            phi    = acos(cosphi)
            sinp2  = 1.d0 - cosphi * cosphi
            if( sinp2 .lt. 0.d0 ) sinp2 = 0.d0
            sinphi = sqrt(sinp2)
            !
            !   LCOL IS SET TRUE IF COLLINEAR GEOMETRY
            !
            lcol = abs(sinphi) .lt. 1.d-6 .or. abs(capr) .lt. 1.d-6
            !
            !   COMPUTE COLLINEAR SATO PARAMETERS
            !
            z(1) = zfh0
            t    = a * (rr(2) - b)
            if( t < exlarg )then
                t     = exp(t)
                t1    = 1.0d0/t
                hsech = 1.0d0/(t+t1)
                zhh   = a1+ha2*((t-t1)*hsech+1.0d0)
            else
                hsech = 0.0d0
                zhh   = a1+a2
            endif
            t    = rr(2) - f
            exx  = c * exp(legacy_d * t * t)
            zhh  = zhh - exx
            z(2) = zhh
            goto_block_50_451: block
                if( .not.lderiv ) exit goto_block_50_451
                do i = 1, 3
                    do j = 1, 3
                        dz(i,j) = 0.d0
                    40 enddo
                enddo
                !
                !  ONLY THE HH SATO PARAMETER HAS A NONZERO DERIVATIVE - W.R.T. R2
                !
                dz(2,2) = 2.d0 * a * a2 * hsech * hsech - 2.d0 * legacy_d * t * exx
            end block goto_block_50_451
            50 continue
            !
            !   COMPUTE NONCOLLINEAR CONTRIBUTIONS TO SATO PARAMETERS
            !
            goto_block_160_550: block
                if( lcol ) exit goto_block_160_550
                !
                !   COMPUTE Z(PHI) FOR FH SATO PARAMETER
                !   DZPHI IS THE DERIVATIVE OF Z(PHI) W.R.T. PHI
                !
                deg = tdeg * phi
                goto_block_60_465: block
                    if( .not.(deg .lt. 0.d0 .or. deg .gt. 90.d0) ) exit goto_block_60_465
                    write(nflag(18), '(/,2X,T5,''In the FH2 surface no. 5A, PHI = '', F10.5, ''degrees'', /,2X,T5,''but PHI must be gr' // &
                        'eater than zero.'')') deg
                    stop 'POT 3'
                end block goto_block_60_465
                60 continue
                if( deg <= deg1 )then
                    zphi  = 0.0d0
                    dzphi = 0.0d0
                else if( deg <= deg2 )then
                    degdif = deg-deg1
                    dd2    = degdif*degdif
                    zphi   = dd2*(c1+c2*degdif)
                    if( lderiv ) dzphi = tdeg*degdif*(2.0d0*c1+3.0d0*c2*degdif)
                else if( deg <= deg3 )then
                    degdif = deg-deg3p
                    zphi   = c3+c4*degdif
                    if( lderiv ) dzphi = c4*tdeg
                else if( deg <= deg4 )then
                    zphi = c5+c6*sinphi*sinphi
                    if( lderiv ) dzphi = 2.0d0*c6*sinphi*cosphi
                else
                    write (nflag(18), '(/,2X,T5,''In the FH2 surface no. 5A, PHI = '', F10.5, ''degrees'', /,2X,T5,''but PHI must be l' // &
                        'ess than 90 degrees.'')') deg
                    stop 'POT 4'
                endif
                !
                !   COMPUTE SCALE FACTOR TO DAMP OUT CORRECTION AT LARGE DISTANCE
                !
                scale1 = capr2 / (g2 + capr2)
                z(1)   = z(1) + scale1 * zphi
                !
                !   COMPUTE SCALE FACTOR FOR HH SATO PARAMETER
                !
                scale2 = (zhh - zhh0) / zhhdif
                !
                !   COMPUTE PHI DEPENDENT PART OF HH SATO PARAMETER
                !
                cphi = (bet + betp * sinp2) * sinp2
                !
                !   COMPUTE HH SATO PARAMTER
                !
                z(2) = z(2) + scale1 * scale2 * cphi
                goto_block_150_549: block
                    if( .not.lderiv ) exit goto_block_150_549
                    !
                    !   COMPUTE DERIVATIVES
                    !   FIRST, DERIVATIVES OF PHI W.R.T. R1, R2, R3
                    !
                    t     = 1.d0 / capr
                    t1    = 0.5d0 * cosphi * t
                    t2    = sgnphi / rr(2)
                    dphi1 = rr(1) * (t1 + t2)
                    dphi3 = rr(3) * (t1 - t2)
                    t1    = 1.d0 / sinphi
                    t2    = t1 * t
                    dphi1 = dphi1 * t2
                    dphi3 = dphi3 * t2
                    dphi2 = rr(2) * cosphi * t1 * (1.d0 / r22 - 0.25d0 / capr2)
                    t1    = (1.d0 - scale1) / capr2
                    t2    = t1 * zphi
                    !
                    !   NEXT, DERIVATIVES OF FH SATO PARAMETERS
                    !
                    dz(1,1) = scale1 * (t2 * rr(1) + dzphi * dphi1)
                    dz(1,2) = scale1 * (-0.5d0 * t2 * rr(2) + dzphi * dphi2)
                    dz(1,3) = scale1 * (t2 * rr(3) + dzphi * dphi3)
                    t1      = t1 * cphi
                    !
                    !   DERIVATIVE OF THE PHI DEPENDENT PART OF HH SATO W.R.T. PHI
                    !
                    dcphi = 2.d0 * (bet + 2.d0 * betp * sinp2) * sinphi * cosphi
                    t2    = scale1 * scale2
                    !
                    !   LAST, DERIVATIVES OF HH SATO PARAMETER
                    !
                    dz(2,1) = t2 * (t1 * rr(1) + dcphi * dphi1)
                    dz(2,2) = dz(2,2) * (1.d0 + scale1 * cphi / zhhdif) + t2 * (-0.5d0 * t1 * rr(2) + dcphi * dphi2)
                    dz(2,3) = t2 * (t1 * rr(3) + dcphi * dphi3)
                end block goto_block_150_549
                150 continue
            end block goto_block_160_550
            160 continue
            !
            !   COMPUTE POTENTIAL FROM LEPS FORM
            !
            engygs = 0.d0
            z(3)   = z(1)
            goto_block_170_567: block
                if( .not.lderiv ) exit goto_block_170_567
                do i = 1, 3
                    dz(3,i)   = dz(1,i)
                    degsdr(i) = 0.d0
                    q(i)      = 0.d0
                    xj(i)     = 0.d0
                    do j = 1, 3
                        dq(i,j) = 0.d0
                        dj(i,j) = 0.d0
                    162 enddo
                165 enddo
            end block goto_block_170_567
            170 continue
            lzero = .true.
            do i = 1, 3
                ex    = exp(beta(i) * (re(i) - rr(i)))
                lzero = lzero .and. ex.eq.0.0d0
                if( ex .eq. 0.0d0 ) cycle
                zp1  = z(i) + 1.d0
                zp3  = zp1 + 2.d0
                op3z = 1.d0 + 3.d0 * z(i)
                t    = 0.25d0 * dm(i) * ex / zp1
                t1   = zp3 * ex
                t2   = op3z * ex
                !
                !   Q AND XJ ARE THE COULUMB AND EXCHANGE PARTS, INDEX I RUNS OVER THE
                !   COORDINATES OF WHICH THEY ARE EXPLICIT FUNCTIONS
                !
                q(i)  = t * (t1 - 2.d0  * op3z)
                xj(i) = t * (t2 - 2.d0  * zp3)
                !
                !   PUT SUM OF Q'S IN ENERGY
                !
                engygs = engygs + q(i)
                goto_block_190_602: block
                    if( .not.lderiv ) exit goto_block_190_602
                    !
                    !   DERIVATIVE OF Q AND ZJ W.R.T. R1, R2, R3
                    !   ELEMENT I,J IS THE DERIVATIVE OF Q(I) W.R.T. RJ
                    !
                    t3 = 2.d0 * t * (ex + 2.d0) / zp1
                    do j = 1, 3
                        dq(i,j) = -t3 * dz(i,j)
                        dj(i,j) = -dq(i,j)
                    180 enddo
                    t       = 2.d0*t*beta(i)
                    dq(i,i) = dq(i,i) - t * (t1 - op3z)
                    dj(i,i) = dj(i,i) - t * (t2 - zp3)
                end block goto_block_190_602
                190 continue
                195 continue
            200 enddo
            goto_block_270_646: block
                if( lzero ) exit goto_block_270_646
                xjx = 0.d0
                !
                !   COMPUTE THE TOTAL EXCHANGE TERM
                !
                do i1 = 1, 3
                    i2  = mod(i1, 3) + 1
                    t   = xj(i1) - xj(i2)
                    xjx = xjx + t * t
                    if( .not.lderiv ) cycle
                    !
                    !   PUT THE DERIVATIVES OF THE TOTAL EXCHANGE TERM IN DEGSDR
                    !
                    do j = 1, 3
                        degsdr(j) = degsdr(j) + t * (dj(i1,j) - dj(i2,j))
                    210 enddo
                    220 continue
                230 enddo
                xjx = sqrt(0.5d0 * xjx)
                !
                !   COMPUTE THE LEPS PART OF THE POTENTIAL ENERGY
                !
                engygs = engygs - xjx
                goto_block_260_645: block
                    if( .not.lderiv ) exit goto_block_260_645
                    !
                    !   COMPUTE DERIVATIVES OF THE LEPS PART
                    !
                    if( abs(xjx) .gt. 1.d-14 ) t = 0.5d0 / xjx
                    do j = 1, 3
                        !
                        !   THE DERIVATIVE OF THE EXCHANGE PART MUST BE DIVIDED BY 2*XJX
                        !
                        degsdr(j) = -t *degsdr(j)
                        do i = 1, 3
                            !
                            !   THEN ADD IN THE DERIVATIVE OF Q
                            !
                            degsdr(j) = degsdr(j) + dq(i,j)
                        240 enddo
                    250 enddo
                end block goto_block_260_645
                260 continue
            end block goto_block_270_646
            270 continue
            !
            !   EZERO ADDED TO PUT ZERO AT REQ FOR REACTANTS
            !
            engygs = engygs + ezero(1)
            !
            !   COMPUTE THE 3-CENTER TERM
            !
            rsum  = rr(1) + rr(3)
            rdif  = rr(1) - rr(3)
            rdif2 = rdif * rdif
            t1    = -alfp * rsum
            ex1   = c3c * exp(t1 * rsum)
            t2    = -alf * rdif
            ex2   = omp3c * exp(t2 * rdif)
            t3    = -alft * rdif2
            ex3   = p3c * exp(t3 * rdif2)
            r1i   = 1.d0 / rr(1)
            r3i   = 1.d0 / rr(3)
            t4    = r1i*r3i
            !
            !   COSTH IS THE ANGLE BETWEEN R1 AND R3
            !
            costh = 0.5d0 * (r12 + r32 - r22) * t4
            t5    = omq3c * costh
            t     = q3c + t5 * costh
            !
            !   E3C CONTAINS THE THREE CENTER CORRECTION TO THE ENERGY
            !
            e3c    = (ex2 + ex3) * t * ex1
            engygs = engygs + e3c
            goto_block_280_696: block
                if( .not.lderiv ) exit goto_block_280_696
                !
                !   COMPUTE DERIVATIVE OF THE 3C TERM
                !   FIRST, DERIVATIVE OF COSTH
                !
                dcos1 = r3i - costh*r1i
                dcos2 = - rr(2)*t4
                dcos3 = r1i - costh*r3i
                t4    = ex2 + ex3
                t5    = t4 * t5
                !
                !   NOW, DERIVATIVES OF E3C
                !
                de3c1     = 2.d0*((t2 * ex2 + 2.d0 * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos1) * ex1
                de3c2     = 2.d0 * t5 * dcos2 * ex1
                de3c3     = 2.d0*((-t2 * ex2 - 2.d0 * t3 * rdif * ex3 + t1 * t4) * t + t5 * dcos3) * ex1
                degsdr(1) = degsdr(1) + de3c1
                degsdr(2) = degsdr(2) + de3c2
                degsdr(3) = degsdr(3) + de3c3
            end block goto_block_280_696
            280 continue
            !
            call eunitzero
            if( nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
            !
        end subroutine pot
        !> @brief Evaluate ANCVRT for the FH2_LEPS_NO5a_1985 potential.
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
        !> @brief Evaluate CARTOU for the FH2_LEPS_NO5a_1985 potential.
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
        !> @brief Evaluate CARTTOR for the FH2_LEPS_NO5a_1985 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer            :: coordinate_status
            double precision   :: cmm, cphi, flm, hhd, hym, phi, pi, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
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
        !> @brief Evaluate EUNITZERO for the FH2_LEPS_NO5a_1985 potential.
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
        !> @brief Evaluate RTOCART for the FH2_LEPS_NO5a_1985 potential.
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
                        20 continue
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
        !> @brief Evaluate DEDCOU for the FH2_LEPS_NO5a_1985 potential.
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
end module chempotpy_fh2_leps_no5a_1985_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_fh2_leps_no5a_1985_pes, only: evaluate
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
