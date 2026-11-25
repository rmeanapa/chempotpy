! System:          F2H
! Functional form: eLEPS (extended London-Eyring-Polanyi-Sato)
! Common name:     JOT-II
! Reference:       N. Jonathon, S. Okuda, and D. Timlin, Mol. Phys.
!                  24, 1143 (1972).
! Number of bodies: 3
! Interface: potlib2001
! Number of electronic states: 1
! Number of derivatives: 1
! prepare_local_context must be called once before any calls to POT are made.  It shoul
! never be called more than once.
! The potential parameters are included in DATA statements.
! Coordinates, potential energy, derivatives, and other information
! is passed through three common blocks:
! The common block PT1CM contains the internal coordinates, the potentia
! energy, and the derivatives with respect to internal coordinates.  The
! common block PT2CM contains the a several variables which control
! various aspects of the calculation as outlined below.  The last common
! block, ASYCM, contains information about the differences in the
! dissociation energies which is useful for certain calculations.
! Internuclear distances should be expressed in bohr and energies, etc.
! returned in Hartree atomic units.
! The coordinates are defined as follows:
!           R(1) = R(H-F)
!           R(2) = R(F-F')
!           R(3) = R(H-F')
! The zero of energy for this potential occurs when the H atom is
! "infinitely" far from the F2 diatom, and the F2 diatom is at
! its equilibrium separation.
! The array DEGSDR contains the analytic first derivatives of the potent
! energy with respect to the internal coordinates if the variable NDER i
! equal to 1 (see not below).  The derivatives are defined as follows:
!           DEGSDR(1) = dV(R)/dR(1)
!           DEGSDR(2) = dV(R)/dR(2)
!           DEGSDR(3) = dV(R)/dR(3)
! where V(R) = V { R(1), R(2), R(3) }.
! The variable NDER controls whether derivatives are calculated.
! If NDER = 1 first derivatives are calculated, otherwise no derivatives
! are calculated.  Note:  the default value of NDER is 1.
! The array NFLAG is unused in this potential surface.
! The common block ASYCM is defined as follows:
!           EASYAB = DE(H-F)
!           EASYBC = DE(F-F')
!           EASYAC = DE(H-F')
! where DE is the dissociation energy in atomic units.
!   Potential parameters' default settings
!                  Variable            Default value
!                  NDER                1
!                  NFLAG(18)           6

!> @brief Implement F2H_LEPS_JOT_1972, a potential energy surface for F2H.
module chempotpy_f2h_leps_jot_1972_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the F2H_LEPS_JOT_1972 potential energy surface for F2H.
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
        real(wp)              :: anuzero
        real(wp)              :: b(3)
        real(wp)              :: canau
        real(wp)              :: cart(natom,3)
        real(wp)              :: cartnu(natom,3)
        real(wp)              :: cevau
        real(wp)              :: cnvrtd
        real(wp)              :: cnvrtde
        real(wp)              :: cnvrte
        real(wp)              :: de(3)
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
        real(wp)              :: hde(3)
        real(wp)              :: hpf(3)
        integer               :: icartr
        integer               :: indexes(natom)
        integer               :: irctnt
        integer               :: ireorder
        integer               :: kediag
        integer               :: keoffd
        integer               :: ksdiag
        integer               :: ksoffd
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
        real(wp)              :: pfs(3)
        real(wp)              :: pft(3)
        real(wp)              :: r(n3atom)
        real(wp)              :: re(3)
        character(len=75)     :: ref(5)
        real(wp)              :: s(3)
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        anuzero = 0.e0_wp; b = 0.e0_wp; canau = 0.e0_wp; cart = 0.e0_wp; cartnu = 0.e0_wp; cevau = 0.e0_wp; cnvrtd = 0.e0_wp; cnvrtde = 0.e0_wp
        cnvrte  = 0.e0_wp; de = 0.e0_wp; deesdr = 0.e0_wp; degsdr = 0.e0_wp; deijdr = 0.e0_wp; descart = 0.e0_wp; descartnu = 0.e0_wp
        dgscart = 0.e0_wp; dgscartnu = 0.e0_wp; dijcart = 0.e0_wp; dijcartnu = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp
        engyes  = 0.e0_wp; engygs = 0.e0_wp; engyij = 0.e0_wp; ezero = 0.e0_wp; hde = 0.e0_wp; hpf = 0.e0_wp; icartr = 0; indexes = 0
        irctnt  = 0; ireorder = 0; kediag = 0; keoffd = 0; ksdiag = 0; ksoffd = 0; legacy_natoms = 0; mder = 0; msurf = 0
        nasurf  = 0; nder = 0; nflag = 0; nulbl = 0; pengyes = 0.e0_wp; pengygs = 0.e0_wp; pengyij = 0.e0_wp; pfs = 0.e0_wp
        pft     = 0.e0_wp; r = 0.e0_wp; re = 0.e0_wp
        ref     = ' '
        s       = 0.e0_wp
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
        de            = [6.12e0_wp, 1.63e0_wp, 6.12e0_wp]
        re            = [0.9171e0_wp, 1.418e0_wp, 0.9171e0_wp]
        b             = [2.219e0_wp, 2.920e0_wp, 2.219e0_wp]
        s             = [0.e0_wp, -0.35e0_wp, 0.e0_wp]
        cevau         = 27.2113961e0_wp
        canau         = 0.529177249e0_wp
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

        !> @brief Evaluate prepare_local_context for the F2H_LEPS_JOT_1972 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: pf
            integer  :: i
            ! Write header potential parameters to UNIT NFLAG(18).
            if( legacy_natoms .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            ! Calculate first derivatives (default).  If no first derivatives are
            ! desired, nder should be set equal to zero after prepare_local_context is called.
            nder = 1
            ! Convert to atomic units and calculate some useful constants.
            do i=1, 3
                de(i)  = de(i)/cevau
                re(i)  = re(i)/canau
                b(i)   = b(i)*canau
                pf     = 0.5e0_wp*de(i)*((1.e0_wp-s(i))/(1.e0_wp+s(i)))
                hde(i) = 0.5e0_wp*de(i)
                hpf(i) = 0.5e0_wp*pf
                pfs(i) = -de(i)*b(i)
                pft(i) = -pf*b(i)
            enddo
            easyab   = de(1)
            easybc   = de(2)
            easyac   = de(3)
            ezero(1) = de(2)
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1)     = 'N. Jonathon, S. Okuda, D. Timlin,'
            ref(2)     = 'Mol. Phys. 24, 1143(1972)'
            indexes(1) = 1
            indexes(2) = 9
            indexes(3) = 9
            irctnt     = 2
            call ancvrt
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the F2H_LEPS_JOT_1972 potential.
        subroutine pot
            ! The coordinates are defined as follows:
            !           R(1) = R(H-F)
            !           R(2) = R(F-F')
            !           R(3) = R(H-F')
            ! The zero of energy for this potential occurs when the H atom is
            ! "infinitely" far from the F2 diatom, and the F2 diatom is at
            ! its equilibrium separation.
            ! The common block ASYCM is defined as follows:
            !           EASYAB = DE(H-F)
            !           EASYBC = DE(F-F')
            !           EASYAC = DE(H-F')
            ! where DE is the dissociation energy in atomic units.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dj1dr1, dj2dr2, dj3dr3, dq1dr1, dq2dr2, dq3dr3, hds1dr1, hds2dr2, hds3dr3, hdt1dr1, hdt2dr2, hdt3dr3, hradi, hs1, hs2, hs3, ht1, ht2, ht3, rad, radi, wj1, wj2, wj3, wq1, wq2, wq3, y1, y2, y3
            call cartou
            call carttor
            ! compute coulomb and exchange terms
            y1  = exp(b(1)*(re(1)-r(1)))
            hs1 = hde(1)*(y1-2.e0_wp)*y1
            ht1 = hpf(1)*(y1+2.e0_wp)*y1
            wq1 = hs1+ht1
            wj1 = hs1-ht1
            y2  = exp(b(2)*(re(2)-r(2)))
            hs2 = hde(2)*(y2-2.e0_wp)*y2
            ht2 = hpf(2)*(y2+2.e0_wp)*y2
            wq2 = hs2+ht2
            wj2 = hs2-ht2
            y3  = exp(b(3)*(re(3)-r(3)))
            hs3 = hde(3)*(y3-2.e0_wp)*y3
            ht3 = hpf(3)*(y3+2.e0_wp)*y3
            wq3 = hs3+ht3
            wj3 = hs3-ht3
            ! compute potential energy
            rad    = sqrt(wj1*(wj1-wj2)+wj2*(wj2-wj3)+wj3*(wj3-wj1))
            engygs = ezero(1)  + wq1+wq2+wq3-rad
            !  compute first derivatives if necessary
            if( nder .eq. 1 )then
                radi      = 1.e0_wp/rad
                hradi     = 0.5e0_wp*radi
                hds1dr1   = pfs(1)*(y1-1.e0_wp)*y1
                hds2dr2   = pfs(2)*(y2-1.e0_wp)*y2
                hds3dr3   = pfs(3)*(y3-1.e0_wp)*y3
                hdt1dr1   = pft(1)*(y1+1.e0_wp)*y1
                hdt2dr2   = pft(2)*(y2+1.e0_wp)*y2
                hdt3dr3   = pft(3)*(y3+1.e0_wp)*y3
                dq1dr1    = hds1dr1+hdt1dr1
                dq2dr2    = hds2dr2+hdt2dr2
                dq3dr3    = hds3dr3+hdt3dr3
                dj1dr1    = hds1dr1-hdt1dr1
                dj2dr2    = hds2dr2-hdt2dr2
                dj3dr3    = hds3dr3-hdt3dr3
                degsdr(1) = dq1dr1-hradi*dj1dr1*(2.e0_wp*wj1-wj2-wj3)
                degsdr(2) = dq2dr2-hradi*dj2dr2*(2.e0_wp*wj2-wj1-wj3)
                degsdr(3) = dq3dr3-hradi*dj3dr3*(2.e0_wp*wj3-wj1-wj2)
            endif
            call eunitzero
            if( nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate ANCVRT for the F2H_LEPS_JOT_1972 potential.
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

        !> @brief Evaluate CARTOU for the F2H_LEPS_JOT_1972 potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            implicit none
            call convert_cartesian_input(cart, nulbl, legacy_natoms, cnvrtd, ireorder == 1, cartnu)
        end subroutine cartou

        !> @brief Evaluate CARTTOR for the F2H_LEPS_JOT_1972 potential.
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

        !> @brief Evaluate EUNITZERO for the F2H_LEPS_JOT_1972 potential.
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

        !> @brief Evaluate RTOCART for the F2H_LEPS_JOT_1972 potential.
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

        !> @brief Evaluate DEDCOU for the F2H_LEPS_JOT_1972 potential.
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
end module chempotpy_f2h_leps_jot_1972_pes

!> @brief Expose F2H_LEPS_JOT_1972 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_f2h_leps_jot_1972_pes, only: evaluate
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
