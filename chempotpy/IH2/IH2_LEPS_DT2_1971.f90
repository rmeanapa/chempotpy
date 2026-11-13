!   System:          IH2
!   Functional form:
!   Common name:     RMC
!   References:
!   Note: This is the three-dimensional version of the IH2
!         RMC potential energy surface; the potential energy is
!         defined as E = Ecollinear + Ebend.
! Number of bodies: 3
! Interface: potlib2001
! Number of electronic states: 1
! Number of derivatives: 1
!   prepare_local_context must be called once before any calls to POT.
!   The potential parameters are included in the block data subprogram P
!   Coordinates, potential energy, and derivatives are passed
!   The potential energy in the three asymptotic valleys are
!   stored in the common block ASYCM:
!   The potential energy in the AB valley, EASYAB, is equal to the poten
!   energy of the H "infinitely" far from the IH diatomic, with the
!   IH diatomic at its equilibrium configuration.  Similarly, the terms
!   EASYBC and EASYAC represent the H2 and the IH asymptotic valleys,
!   respectively.
!   All the information passed through the common blocks PT1CM and ASYCM
!   is in Hartree atomic units.
!        This potential is written such that:
!                       R(1) = R(I-H)
!                       R(2) = R(H-H)
!                       R(3) = R(H-I)
!   The classical potential energy is set equal to zero for the I
!   infinitely far from the equilibrium H2 diatomic.
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

!> @brief Implement IH2_LEPS_DT2_1971, a potential energy surface for IH2.
module chempotpy_ih2_leps_dt2_1971_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the IH2_LEPS_DT2_1971 potential energy surface for IH2.
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
        real(wp)              :: a
        real(wp)              :: a2(2)
        real(wp)              :: alfd
        real(wp)              :: alpha(2)
        real(wp)              :: anuzero
        real(wp)              :: beta
        real(wp)              :: c(2,2)
        real(wp)              :: cart(natom,3)
        real(wp)              :: cartnu(natom,3)
        real(wp)              :: clbh
        real(wp)              :: cnvrtd
        real(wp)              :: cnvrtde
        real(wp)              :: cnvrte
        real(wp)              :: dadp(2)
        real(wp)              :: dbcoor(2)
        real(wp)              :: de
        real(wp)              :: deesdr(n3atom,isurf)
        real(wp)              :: degsdr(n3atom)
        real(wp)              :: deijdr(n3atom,jsurf)
        real(wp)              :: descart(natom,3,isurf)
        real(wp)              :: descartnu(natom,3,isurf)
        real(wp)              :: dgscart(natom,3)
        real(wp)              :: dgscartnu(natom,3)
        real(wp)              :: dijcart(natom,3,jsurf)
        real(wp)              :: dijcartnu(natom,3,jsurf)
        real(wp)              :: dmdn
        real(wp)              :: dpdr(2)
        real(wp)              :: easyab
        real(wp)              :: easyac
        real(wp)              :: easybc
        real(wp)              :: engyes(isurf)
        real(wp)              :: engygs
        real(wp)              :: engyij(jsurf)
        real(wp)              :: ezero(isurf+1)
        real(wp)              :: gamma
        integer               :: icartr
        integer               :: indexes(natom)
        integer               :: irctnt
        integer               :: ireorder
        integer               :: irxn
        real(wp)              :: junk(4)
        integer               :: kediag
        integer               :: keoffd
        integer               :: ksdiag
        integer               :: ksoffd
        integer               :: l(2)
        real(wp)              :: legacy_d(2)
        integer               :: legacy_natoms
        integer               :: mder
        integer               :: msurf
        integer               :: nasurf(isurf+1,isurf+1)
        integer               :: nder
        integer               :: nflag(20)
        integer               :: nulbl(natom)
        real(wp)              :: pengyes(isurf)
        real(wp)              :: pengygs
        real(wp)              :: pengyij(jsurf)
        real(wp)              :: pid2
        real(wp)              :: pid4
        real(wp)              :: pre(2)
        real(wp)              :: r(n3atom)
        real(wp)              :: rbb
        real(wp)              :: re(2)
        character             :: ref(5)
        real(wp)              :: req(3)
        real(wp)              :: rhh
        real(wp)              :: rr(2)
        real(wp)              :: rr1ms1
        real(wp)              :: rr2ms2
        real(wp)              :: rrmr(2)
        real(wp)              :: rrmre(2)
        real(wp)              :: rrs1d2
        real(wp)              :: s(2)
        real(wp)              :: sr1r2
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        a             = 0.e0_wp; a2 = 0.e0_wp; alfd = 0.e0_wp; alpha = 0.e0_wp; anuzero = 0.e0_wp; beta = 0.e0_wp; c = 0.e0_wp; cart = 0.e0_wp
        cartnu        = 0.e0_wp; clbh = 0.e0_wp; cnvrtd = 0.e0_wp; cnvrtde = 0.e0_wp; cnvrte = 0.e0_wp; dadp = 0.e0_wp; dbcoor = 0.e0_wp; de = 0.e0_wp
        deesdr        = 0.e0_wp; degsdr = 0.e0_wp; deijdr = 0.e0_wp; descart = 0.e0_wp; descartnu = 0.e0_wp; dgscart = 0.e0_wp; dgscartnu = 0.e0_wp
        dijcart       = 0.e0_wp; dijcartnu = 0.e0_wp; dmdn = 0.e0_wp; dpdr = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp
        engyes        = 0.e0_wp; engygs = 0.e0_wp; engyij = 0.e0_wp; ezero = 0.e0_wp; gamma = 0.e0_wp; icartr = 0; indexes = 0; irctnt = 0
        ireorder      = 0; irxn = 0; junk = 0.e0_wp; kediag = 0; keoffd = 0; ksdiag = 0; ksoffd = 0; l = 0; legacy_d = 0.e0_wp
        legacy_natoms = 0; mder = 0; msurf = 0; nasurf = 0; nder = 0; nflag = 0; nulbl = 0; pengyes = 0.e0_wp
        pengygs       = 0.e0_wp; pengyij = 0.e0_wp; pid2 = 0.e0_wp; pid4 = 0.e0_wp; pre = 0.e0_wp; r = 0.e0_wp; rbb = 0.e0_wp; re = 0.e0_wp
        ref           = ' '
        req           = 0.e0_wp
        rhh           = 0.e0_wp
        rr            = 0.e0_wp
        rr1ms1        = 0.e0_wp
        rr2ms2        = 0.e0_wp
        rrmr          = 0.e0_wp
        rrmre         = 0.e0_wp
        rrs1d2        = 0.e0_wp
        s             = 0.e0_wp
        sr1r2         = 0.e0_wp
        nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0], shape(nasurf))
        nder          = 0
        nflag         = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        anuzero       = 0.0e0_wp
        icartr        = 3
        msurf         = 0
        mder          = 1
        nulbl         = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        legacy_natoms = 3
        alpha         = [1.789e0_wp, 1.974e0_wp]
        legacy_d      = [73.682e0_wp, 109.499e0_wp]
        re            = [1.604e0_wp, 0.742e0_wp]
        s             = [1.6783e0_wp, 1.2801e0_wp]
        rr            = [2.750e0_wp, 1.750e0_wp]
        alfd          = 1.500e0_wp
        clbh          = 35.879e0_wp
        l             = [12, 6]
        irxn          = 1
        de            = 73.682e0_wp
        beta          = 1.789e0_wp
        req           = [1.604e0_wp, 0.742e0_wp, 1.604e0_wp]
        a             = 0.26e0_wp
        gamma         = 0.50e0_wp
        rbb           = 0.7412871e0_wp
        v             = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
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
            do istate=1, nstates
                p(istate) = pengygs*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = pengygs*hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = dgscart(iatom,idir)*hartree_per_bohr_to_ev_per_angstrom
                enddo
            enddo
        elseif( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif

        contains

        !> @brief Evaluate prepare_local_context for the IH2_LEPS_DT2_1971 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i, it
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            real(wp)            :: nab, nab1, nbc, nbc1
            !   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr
            real(wp), parameter :: ckcal = 627.5095e0_wp
            real(wp), parameter :: cangs = 0.52917706e0_wp
            if( legacy_natoms .gt. 25 )then
                write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !      COLLINEAR POTENTIAL PARAMETERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !   Initialize the parameters for the collinear part of the potential
            call prepare_local_context_2
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !      BENDING POTENTIAL PARAMETERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !   Echo the parameters for the bend potential to the
            !   file linked to UNIT NFLAG(18)
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !       CONVERT TO ATOMIC UNITS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            do it = 1, 3
                req(it) = req(it)/cangs
            enddo
            beta = beta * cangs
            de   = de/ckcal
            a    = a/cangs
            rhh  = 1.4007977e0_wp
            rbb  = rbb /cangs
            !CCCCCCCCCCC
            ezero(1) = legacy_d(2)
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1)     = 'No Reference'
            indexes(1) = 53
            indexes(2) = 1
            indexes(3) = 1
            irctnt     = 2
            call ancvrt
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the IH2_LEPS_DT2_1971 potential.
        subroutine pot
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the IH diatomic, with the
            !   IH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the IH asymptotic valleys,
            !   respectively.
            !        This potential is written such that:
            !                       R(1) = R(I-H)
            !                       R(2) = R(H-H)
            !                       R(3) = R(H-I)
            !   The classical potential energy is set equal to zero for the I
            !   infinitely far from the equilibrium H2 diatomic.
            !CCCCCCCCCCCCCCC
            !CCCCCCCCCCCCCCC
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: bcoor, c2, dac, ractmp, rob, stuff, vabbc, vac
            integer             :: it
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            real(wp)            :: nab, nab1, nbc, nbc1
            !   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr
            real(wp), parameter :: ckcal = 627.5095e0_wp
            real(wp), parameter :: cangs = 0.52917706e0_wp
            call cartou
            call carttor
            !   Check the values of NASURF and NDER for validity.
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        FIRST CALCULATE SOME USEFUL NUMBERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
            !        AB AND BC RESPECTIVELY. (EQUATION 2.2)
            nab = exp((req(1) - r(1))/a)
            if( nab .lt. 1.e-15_wp ) nab = 1.e-15_wp
            nbc = exp((req(2) - r(2))/a)
            if( nbc .lt. 1.e-15_wp ) nbc = 1.e-15_wp
            !       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
            !       SEE EQUATION 2.7A, 2.7B
            c2 = 1.0e0_wp-nab/nbc
            if( abs(c2) <= 1.0e-14_wp )then
                nab1 = 0.5e0_wp
                nbc1 = 0.5e0_wp
            else
                c2   = c2/nab
                nab1 = (2.0e0_wp+c2-sqrt(4.0e0_wp+c2*c2))/(2.0e0_wp*c2)
                nbc1 = 1.0e0_wp-nab1
            endif
            !          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
            !          FIRST TRAP OUT ANY ZERO ARGUEMENTS
            if( nab1*nbc1 > 0.0e0_wp )then
                stuff = a * log(nab1*nbc1)
                rob   = (rbb - stuff)/(rhh - stuff)
            else
                rob = 1.0e0_wp
            endif
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        CALCULATE BENDING CORRECTION
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        CALCULATE V(R(3))
            junk(3) = exp(-beta*(r(3)-req(3))/rob)
            vac     = gamma * de * junk(3)* (1.e0_wp + 0.5e0_wp*junk(3))
            !        CALCULATE V(R(1)+R(2))
            junk(4) = exp(beta*(req(3)-r(2)-r(1))/rob)
            vabbc   = gamma * de * junk(4) * (1.e0_wp + 0.5e0_wp*junk(4))
            !        HERE IS THE BENDING CORRECTION
            bcoor = vac - vabbc
            !          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
            !          OF BCOOR WITH RESPECT TO R(1), R(2), AND R(3).
            !        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
            if( nder .eq. 1 )then
                dac       = -gamma*de*beta*junk(3)*(1.e0_wp+junk(3))/rob
                dbcoor(1) = gamma * de * beta * junk(4) * (1.e0_wp + junk(4))/rob
                dbcoor(2) = dbcoor(1)
            endif
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        NOW CALCULATE THE COLLINEAR ENERGY
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR
            !          GEOMETRY AND CALL POT
            ractmp = r(3)
            r(3)   = r(2) + r(1)
            call pot2
            !   Reset R(3) the its correct value
            r(3) = ractmp
            !   Add the energy terms together, i.e., E = Ecollinear + Ebend
            engygs = engygs + bcoor
            if( nder .ne. 1 )then
                call eunitzero
                if( nder .ne. 0 )then
                    call rtocart
                    call dedcou
                endif
                return
            endif
            !   If NDER = 1, then compute the derivatives of the energy with respect
            !   the internal coordinates, i.e., dE/dRi = dEcollinear/dRi + dEbend/dR
            degsdr(1) = degsdr(1) + degsdr(3)
            degsdr(2) = degsdr(2) + degsdr(3)
            degsdr(3) = dac
            do it = 1, 2
                degsdr(it) = dbcoor(it) + degsdr(it)
            enddo
            call eunitzero
            if( nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate prepare_local_context_2 for the IH2_LEPS_DT2_1971 potential.
        subroutine prepare_local_context_2
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i, m, n
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr
            real(wp), parameter :: ckcal = 627.5095e0_wp
            real(wp), parameter :: cangs = 0.52917706e0_wp
            !   This subprogram initializes the constants
            !   for the collinear part of the potential.
            !   prepare_local_context_2 must be called once before any calls to POT2.
            pid4 = atan(1.e0_wp)
            pid2 = 2.e0_wp*pid4
            m    = irxn
            n    = 1
            if( irxn .eq. 1 ) n = 2
            c(m,1)                   = clbh/legacy_d(m)
            if( irxn .eq. 1 ) c(m,1) = c(m,1) + (legacy_d(m)-legacy_d(n))/legacy_d(m)
            c(m,2)                   = 0.0e0_wp
            dmdn                     = legacy_d(m)/legacy_d(n)
            c(n,2)                   = 2.0e0_wp*dmdn*c(m,1)*real(l(n)-l(m), wp) + 2.0e0_wp*real(l(n), wp)*(1.0e0_wp- &
                & dmdn)
            c(n,1) = 1.0e0_wp-dmdn*(1.0e0_wp-c(m,1))+pid4*pid4*c(n,2)
            !   Echo the potential parameters to the file linked to UNIT NFLAG(18).
            !   CONVERT TO ATOMIC UNITS
            do i = 1, 2
                rr(i)       = rr(i)/cangs
                s(i)        = s(i)/cangs
                re(i)       = re(i)/cangs
                alpha(i)    = alpha(i)*cangs
                legacy_d(i) = legacy_d(i)/ckcal
            enddo
            alfd = alfd*cangs
            !   CALCULATE USEFUL CONSTANTS
            rr1ms1 = rr(1)-s(1)
            rr2ms2 = rr(2)-s(2)
            rrs1d2 = rr1ms1/rr2ms2
            sr1r2  = sqrt(rr1ms1**2+rr2ms2**2)
            do n=1, 2
                rrmre(n) = rr(n)-re(n)
                pre(n)   = sr1r2-rrmre(n)
            enddo
            !   Initialize the asymptotic energy terms
            easyab = legacy_d(1)
            easybc = legacy_d(2)
            easyac = legacy_d(1)
        end subroutine prepare_local_context_2

        !> @brief Evaluate POT2 for the IH2_LEPS_DT2_1971 potential.
        subroutine pot2
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: alf, cs2pp, daldr, dgadr, ds4dr, dsldr, dxidr, fac, gam, phip, quan, quan2, rr12, sn2pp, sn2pp3, sn2pp4, sn2ppl, tn2pp, xi
            integer             :: i, j, lm3
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr
            real(wp), parameter :: ckcal = 627.5095e0_wp
            real(wp), parameter :: cangs = 0.52917706e0_wp
            !   This subprogram initializes the constants
            !   for the collinear part of the potential.
            !   PREPT2 must be called once before any calls to POT2.
            logical :: lrr1, lrr2
            !   This subprogram computes the collinear part of the potential energy
            rrmr(1) = rr(1)-r(1)
            rrmr(2) = rr(2)-r(2)
            lrr1    = rrmr(1).le.0.0e0_wp
            lrr2    = rrmr(2).le.0.0e0_wp
            if( nder == 1 ) degsdr = 0.0e0_wp
            engygs = ezero(1)
            if( lrr1 .and. lrr2 ) return
            if( lrr1 .or. lrr2 )then
                j = 1
                if( lrr1 ) j              = 2
                quan   = 1.0e0_wp-exp(-alpha(j)*(r(j)-re(j)))
                engygs = legacy_d(j)*(quan**2-1.0e0_wp)+ezero(1)
                if( nder == 1 ) degsdr(j) = 2.0e0_wp*legacy_d(j)*alpha(j)*quan*(1.0e0_wp-quan)
            else
                phip   = atan(rrs1d2*rrmr(2)/rrmr(1))
                sn2pp  = sin(2.0e0_wp*phip)
                cs2pp  = cos(2.0e0_wp*phip)
                tn2pp  = sn2pp/cs2pp
                sn2pp3 = sn2pp**3
                sn2pp4 = sn2pp3*sn2pp
                i      = 1
                if( phip > pid4 ) i = 2
                a2(i)   = c(i,1)+c(i,2)*(phip-pid2)*phip
                dadp(i) = c(i,2)*(2.0e0_wp*phip-pid2)
                rr12    = sqrt(rrmr(1)**2+rrmr(2)**2)
                lm3     = l(i)-3
                sn2ppl  = sn2pp3*sn2pp**lm3
                alf     = alpha(i)+(alfd-alpha(i))*sn2pp4
                gam     = legacy_d(i)*(1.0e0_wp-a2(i)*sn2ppl)
                xi      = pre(i)*sn2pp4+rrmre(i)-rr12
                quan    = 1.0e0_wp-exp(-alf*xi)
                quan2   = quan*quan
                engygs  = gam*(quan2-1.0e0_wp)+ezero(1)
                if( nder == 1 )then
                    fac     = rrs1d2/(rrmr(1)**2+(rrs1d2*rrmr(2))**2)
                    dpdr(1) = fac*rrmr(2)
                    dpdr(2) = -fac*rrmr(1)
                    do j = 1, 2
                        ds4dr     = 8.0e0_wp*sn2pp4/tn2pp*dpdr(j)
                        dsldr     = 2.0e0_wp*real(l(i), wp)*sn2ppl/tn2pp*dpdr(j)
                        daldr     = (alfd-alpha(i))*ds4dr
                        dgadr     = -legacy_d(i)*a2(i)*dsldr-legacy_d(i)*dadp(i)*sn2ppl*dpdr(j)
                        dxidr     = pre(i)*ds4dr+rrmr(j)/rr12
                        degsdr(j) = dgadr*(quan2-1.0e0_wp)+2.0e0_wp*gam*quan*(1.0e0_wp-quan)*(alf*dxidr+xi*daldr)
                    enddo
                endif
            endif
        end subroutine pot2

        !> @brief Evaluate ANCVRT for the IH2_LEPS_DT2_1971 potential.
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

        !> @brief Evaluate CARTOU for the IH2_LEPS_DT2_1971 potential.
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

        !> @brief Evaluate CARTTOR for the IH2_LEPS_DT2_1971 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            use chempotpy_kinds,              only: wp
            implicit none
            integer            :: coordinate_status
            real(wp)           :: cmm, cphi, flm, hhd, hym, phi, pi, q, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
            integer            :: i, ind, k
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            call cartesian_to_potlib_coordinates(cartnu, legacy_natoms, icartr, r, coordinate_status)
            if( coordinate_status /= 0 )then
                write(nflag(18), '(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') icartr
                stop
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the IH2_LEPS_DT2_1971 potential.
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

        !> @brief Evaluate RTOCART for the IH2_LEPS_DT2_1971 potential.
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

        !> @brief Evaluate DEDCOU for the IH2_LEPS_DT2_1971 potential.
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
end module chempotpy_ih2_leps_dt2_1971_pes

!> @brief Expose IH2_LEPS_DT2_1971 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_ih2_leps_dt2_1971_pes, only: evaluate
    use chempotpy_kinds,                 only: wp
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
