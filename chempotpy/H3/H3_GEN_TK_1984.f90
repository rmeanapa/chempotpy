!   System:      H3
!   Common name: TK
!   Reference:       legacy_d. G. Truhlar and A. Kuppermann
!                    J. Chem. Phys. 56, 2232 (1984).
! Number of bodies: 3
! Interface: potlib2001
! Number of electronic states: 1
! Number of derivatives: 1
!   prepare_local_context must be called once before any calls to POT.
!   The potential parameters are included in DATA statements.
!   and in the block data subprogram PTPACM.
!   Coordinates, potential energy, and derivatives are passed
!   COORDINATES:
!            R(1) = RAB
!            R(2) = RBC
!            R(3) = RAC
!   The potential energy in the three asymptotic valleys are
!   stored in the common block ASYCM:
!   The potential energy in the AB valley, EASYAB, is equal to the poten
!   energy of the H "infinitely" far from the H2 diatomic, with the
!   H2 diatomic at its equilibrium configuration.  For this potential en
!   surface the potential energy in the three asymptotic valleys are equ
!   All the information passed through the common blocks PT1CM and ASYCM
!   is in Hartree atomic units.
!   The potential energy is defined to be equal to zero when one H is
!   "infinitely" far from the H2 diatomic, and the H2 diatomic bond
!   length is equal to the H2 equilibrium bond length.
!   The flags that indicate what calculations should be carried out in
!   the potential routine are passed through the common block PT2CM:
!   where:
!        NASURF - which electronic states availalble
!                 (1,1) = 1 as only gs state available
!        NDER  - the order of the derivatives of the energy with respect
!                the coordinates to be calculates.
!        NDER  = 1 => calculate first derivatives
!                This is the only option supported in this version
!                of the potential.
!        NFLAG - these integer values can be used to flag options
!                within the potential;
!   Potential parameters' default settings
!                  Variable            Default value
!                  NDER                1
!                  NFLAG(18)           6

!> @brief Implement H3_GEN_TK_1984, a potential energy surface for H3.
module chempotpy_h3_gen_tk_1984_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_TK_1984 potential energy surface for H3.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units,  only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
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
        real(wp)              :: a2
        real(wp)              :: alfd
        real(wp)              :: alfdif
        real(wp)              :: alfn
        real(wp)              :: anuzero
        real(wp)              :: beta
        real(wp)              :: cart(natom,3)
        real(wp)              :: cartnu(natom,3)
        real(wp)              :: cnvrtd
        real(wp)              :: cnvrtde
        real(wp)              :: cnvrte
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
        real(wp)              :: junk(4)
        integer               :: kediag
        integer               :: keoffd
        integer               :: ksdiag
        integer               :: ksoffd
        integer               :: l
        real(wp)              :: legacy_d
        integer               :: legacy_natoms
        integer               :: mder
        integer               :: msurf
        real(wp)              :: nab
        real(wp)              :: nab1
        integer               :: nasurf(isurf+1,isurf+1)
        real(wp)              :: nbc
        real(wp)              :: nbc1
        integer               :: nder
        integer               :: nflag(20)
        integer               :: nulbl(natom)
        real(wp)              :: pengyes(isurf)
        real(wp)              :: pengygs
        real(wp)              :: pengyij(jsurf)
        real(wp)              :: r(n3atom)
        real(wp)              :: r2t
        real(wp)              :: rbb
        character(len=75)     :: ref(5)
        real(wp)              :: req(3)
        real(wp)              :: rhh
        real(wp)              :: rn
        real(wp)              :: rnmsn
        real(wp)              :: rnmxn
        real(wp)              :: sn
        real(wp)              :: xn
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        a         = 0.e0_wp; a2 = 0.e0_wp; alfd = 0.e0_wp; alfdif = 0.e0_wp; alfn = 0.e0_wp; anuzero = 0.e0_wp; beta = 0.e0_wp; cart = 0.e0_wp
        cartnu    = 0.e0_wp; cnvrtd = 0.e0_wp; cnvrtde = 0.e0_wp; cnvrte = 0.e0_wp; dbcoor = 0.e0_wp; de = 0.e0_wp; deesdr = 0.e0_wp
        degsdr    = 0.e0_wp; deijdr = 0.e0_wp; descart = 0.e0_wp; descartnu = 0.e0_wp; dgscart = 0.e0_wp; dgscartnu = 0.e0_wp; dijcart = 0.e0_wp
        dijcartnu = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp; engyes = 0.e0_wp; engygs = 0.e0_wp; engyij = 0.e0_wp
        ezero     = 0.e0_wp; gamma = 0.e0_wp; icartr = 0; indexes = 0; irctnt = 0; ireorder = 0; junk = 0.e0_wp; kediag = 0
        keoffd    = 0; ksdiag = 0; ksoffd = 0; l = 0; legacy_d = 0.e0_wp; legacy_natoms = 0; mder = 0; msurf = 0; nab = 0.e0_wp
        nab1      = 0.e0_wp; nasurf = 0; nbc = 0.e0_wp; nbc1 = 0.e0_wp; nder = 0; nflag = 0; nulbl = 0; pengyes = 0.e0_wp; pengygs = 0.e0_wp
        pengyij   = 0.e0_wp; r = 0.e0_wp; r2t = 0.e0_wp; rbb = 0.e0_wp
        ref       = ' '
        req       = 0.e0_wp
        rhh       = 0.e0_wp
        rn        = 0.e0_wp
        rnmsn     = 0.e0_wp
        rnmxn     = 0.e0_wp
        sn        = 0.e0_wp
        xn        = 0.e0_wp
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
        de            = 109.47e0_wp
        beta          = 1.97354
        req           = [0.741287e0_wp, 0.741287e0_wp, 0.741287e0_wp]
        a             = 0.26e0_wp
        gamma         = 0.5e0_wp
        rbb           = 0.741270e0_wp
        alfd          = 0.70152e0_wp
        a2            = 0.08939e0_wp
        sn            = 1.765e0_wp
        legacy_d      = 109.47e0_wp
        alfn          = 1.04435e0_wp
        xn            = 1.40083e0_wp
        l             = 3
        rnmsn         = 1.85630e0_wp
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
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives above first order are unavailable')
        endif

        contains

        !> @brief Evaluate prepare_local_context for the H3_GEN_TK_1984 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, it
            !   Conversion constants
            !   BOHR (Angstrom/bohr) converts Angstroms to Bohr
            !   HARTRE (kcal/Hartree) converts from Hartree to kcal/mol
            real(wp), parameter :: bohr   = 0.52917706e0_wp
            real(wp), parameter :: hartre = hartree_to_kilocalorie_per_mole
            if( legacy_natoms .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !      SET COLLINEAR POTENTIAL PARAMETERS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            call prepare_local_context_2
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !   Echo the bending potential parameters
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !   DE    = dissociation energy in kcal/mol
            !   BETA  = Morse beta in reciprocal Angstroms
            !   A     = Pauling parameter in Angstroms
            !   GAMMA = gamma for the bend correction
            !           gamma is dimensionless and ranges from 0.5 to 0.65
            !   RBB   = parameter used in the scaling calculation for the
            !           bending correction
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !       CONVERT TO ATOMIC UNITS
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            do it = 1, 3
                req(it) = req(it)/bohr
            enddo
            beta = beta * bohr
            de   = de/hartre
            a    = a/bohr
            rhh  = 1.4007977e0_wp
            rbb  = rbb /bohr
            !CCCCCCCCCC
            !                      *
            !   The vale of legacy_d is   *
            !   set in prepare_local_context_2     *
            !                      *
            ezero(1) = legacy_d
            do i=1, 5
                ref(i) = ' '
            enddo
            !   Reference:
            ref(1)     = 'legacy_d. G. Truhlar and A. Kuppermann'
            ref(2)     = 'J. Chem. Phys. 56, 2232(1984)'
            indexes(1) = 1
            indexes(2) = 1
            indexes(3) = 1
            irctnt     = 2
            call ancvrt
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the H3_GEN_TK_1984 potential.
        subroutine pot
            !   COORDINATES:
            !            R(1) = RAB
            !            R(2) = RBC
            !            R(3) = RAC
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the H2 diatomic, with the
            !   H2 diatomic at its equilibrium configuration.  For this potential en
            !   surface the potential energy in the three asymptotic valleys are equ
            !   All the information passed through the common blocks PT1CM and ASYCM
            !   is in Hartree atomic units.
            !   The potential energy is defined to be equal to zero when one H is
            !   "infinitely" far from the H2 diatomic, and the H2 diatomic bond
            !   length is equal to the H2 equilibrium bond length.
            !CCCCCCCCCCCCCCC
            !CCCCCCCCCCCCCCC
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: bcoor, c, dac, ractmp, rob, stuff, vabbc, vac
            integer  :: it
            !   Conversion constants
            !   BOHR (Angstrom/bohr) converts Angstroms to Bohr
            !   HARTRE (kcal/Hartree) converts from Hartree to kcal/mol
            call cartou
            call carttor
            !   Check the values of NASURF and NDER for validity.
            if( nasurf(1,1) .eq. 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            if( nder .gt. 1 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'POT 2')
                return
            endif
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
            !          OF BCOOR WITH RESPECT TO R(1)  R(2) AND R(3).
            !        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))
            dac       = -gamma*de*beta*junk(3)*(1.e0_wp+junk(3))/rob
            dbcoor(1) = gamma * de * beta * junk(4) * (1.e0_wp + junk(4))/rob
            dbcoor(2) = dbcoor(1)
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !        NOW CALCULATE THE COLLINEAR ENERGY
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR
            !          GEOMETRY AND CALL POT
            ractmp = r(3)
            r(3)   = r(2) + r(1)
            call pot2
            r(3)   = ractmp
            engygs = engygs + bcoor
            call eunitzero
            if( nder .ne. 0 )then
                degsdr(1) = degsdr(1) + degsdr(3)
                degsdr(2) = degsdr(2) + degsdr(3)
                degsdr(3) = dac
                do it = 1, 2
                    degsdr(it) = dbcoor(it) + degsdr(it)
                enddo
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate prepare_local_context_2 for the H3_GEN_TK_1984 potential.
        subroutine prepare_local_context_2
            use chempotpy_kinds, only: wp
            implicit none
            real(wp), parameter :: hartre = hartree_to_kilocalorie_per_mole
            !         COORDINATES,POTENTIAL ENERGY, AND DERIVATIVES ARE
            !  THIS SUBPROGRAM SHOULD BE CALLED
            !  ONCE BEFORE THE FIRST CALL TO POT
            rn       = rnmsn+sn
            rnmxn    = rn-xn
            alfdif   = alfd-alfn
            r2t      = 1.4142136e0_wp*rnmsn
            legacy_d = legacy_d /hartre
            !    Set the values of the classical energy in the three asymptotic vall
            easyab = legacy_d
            easybc = easyab
            easyac = easyab
        end subroutine prepare_local_context_2

        !> @brief Evaluate POT2 for the H3_GEN_TK_1984 potential.
        subroutine pot2
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: alf, brack, ct, ctt, dalfx, dalfy, dgamx, dgamy, dsqigx, dsqigy, gam, quan, rcrc, rnmx, rnmy, squig, st, stt, stt4, sttl, t1, t2, theta
            integer  :: lm4
            logical  :: rabl, rbcl
            !         COORDINATES,POTENTIAL ENERGY, AND DERIVATIVES ARE
            !         PASSED THROUGH COMMON /PT1CM/ R(3), ENGYGS, DEGSDR(3)
            ! THIS VERSION WORKS ONLY FOR L .LT. 4
            rabl = r(1).gt.rn
            rbcl = r(2).gt.rn
            if( rabl .or. rbcl )then
                gam = legacy_d
                alf = alfn
                if( rabl ) squig = r(2)-xn
                if( rbcl ) squig = r(1)-xn
                if( rabl .and. rbcl ) squig = squig+r(2)-rn
            else
                rnmx  = rn-r(1)
                rnmy  = rn - r(2)
                theta = atan(rnmy/rnmx)
                st    = sin(theta)
                ct    = cos(theta)
                stt   = 2.e0_wp*st*ct
                ctt   = ct*ct - st*st
                lm4   = l-4
                sttl  = stt**l
                stt4  = sttl*stt**(-lm4)
                gam   = legacy_d*(1.e0_wp-a2*sttl)
                alf   = alfn+alfdif*stt4
                brack = r2t-rnmxn
                rcrc  = 2.e0_wp*rn*rn+r(1)*r(1)+r(2)*r(2)-2.e0_wp*rn*(r(1)+r(2))
                squig = brack*stt4+rnmxn-sqrt(rcrc)
                ! ENERGY AND DERIVATIVE INFORMATION
                t1     = sttl*ctt*3.e0_wp*a2*legacy_d
                dgamx  = -t1/rnmx
                dgamy  = t1/rnmy
                t2     = stt4*ctt*4.e0_wp
                t1     = t2*alfdif
                dalfx  = t1/rnmx
                dalfy  = -t1/rnmy
                t1     = t2*brack
                dsqigx = ct + t1/rnmx
                dsqigy = st - t1/rnmy
            endif
            t1     = exp(-alf*squig)
            quan   = 1.e0_wp-t1
            t2     = quan*quan - 1.e0_wp
            engygs = gam*t2 + ezero(1)
            if( nder .eq. 0 ) return
            ! EXCLUSIVELY DERIVATIVES
            if( rabl .or. rbcl )then
                t2                   = 2.0e0_wp*alf*gam*t1*quan
                degsdr(1)            = 0.0e0_wp
                degsdr(2)            = 0.0e0_wp
                if( rabl ) degsdr(2) = t2
                if( rbcl ) degsdr(1) = t2
                degsdr(3)            = 0.0e0_wp
            else
                degsdr(1) = dgamx*t2
                degsdr(2) = dgamy * t2
                t2        = t1*quan*2.e0_wp*gam
                degsdr(1) = degsdr(1) + t2*(alf*dsqigx + squig*dalfx)
                degsdr(2) = degsdr(2) + t2*(alf*dsqigy + squig*dalfy)
                degsdr(3) = 0.e0_wp
            endif
        end subroutine pot2

        !> @brief Evaluate ANCVRT for the H3_GEN_TK_1984 potential.
        subroutine ancvrt
            use chempotpy_atomic_elements, only: all_atomic_numbers_valid
            use chempotpy_kinds,           only: wp
            implicit none
            real(wp), parameter :: htomillh = 1000.e0_wp
            real(wp), parameter :: htoev    = 27.2113834e0_wp
            real(wp), parameter :: htokcal  = hartree_to_kilocalorie_per_mole
            real(wp), parameter :: htowave  = 219474.631e0_wp
            real(wp), parameter :: htokj    = 2625.49962e0_wp
            real(wp), parameter :: bohr_a   = .5291772083e0_wp
            integer             :: inu, isum
            if( .not. all_atomic_numbers_valid(indexes(1:legacy_natoms)) )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'potential declares an invalid atomic number')
                return
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
            do inu=1, size(nulbl)
                isum = isum+nulbl(inu)
            enddo
            ireorder = 0
            if( isum .ne. 0 ) ireorder = 1
        end subroutine ancvrt

        !> @brief Evaluate CARTOU for the H3_GEN_TK_1984 potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            implicit none
            call convert_cartesian_input(cart, nulbl, legacy_natoms, cnvrtd, ireorder == 1, cartnu)
        end subroutine cartou

        !> @brief Evaluate CARTTOR for the H3_GEN_TK_1984 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(cartnu, legacy_natoms, icartr, r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the H3_GEN_TK_1984 potential.
        subroutine eunitzero
            implicit none
            integer :: i, j
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

        !> @brief Evaluate RTOCART for the H3_GEN_TK_1984 potential.
        subroutine rtocart
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: n3atom = 75
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
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the H3_GEN_TK_1984 potential.
        subroutine dedcou
            implicit none
            integer :: i, j, k
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
end module chempotpy_h3_gen_tk_1984_pes

!> @brief Expose H3_GEN_TK_1984 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_tk_1984_pes, only: evaluate
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
