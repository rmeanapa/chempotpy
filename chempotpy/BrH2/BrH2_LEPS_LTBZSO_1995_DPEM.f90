!   System:          BrH2
!   Functional form: extended LEPS plus three-center term
!   Common name:     BrH2SO
!   Interface:       potlib2001
!   Reference:        G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g.
!                     J. Phys. Chem. 99, 207 (1995)
!   Number of bodies:              3
!   Number of derivatives:         0
!   Number of electronic surfaces: 2

!> @brief Implement BrH2_LEPS_LTBZSO_1995, a diabatic potential matrix for BrH2.
module chempotpy_brh2_leps_ltbzso_1995_dpem_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the BrH2_LEPS_LTBZSO_1995 diabatic potential matrix for BrH2.
    subroutine evaluate( x, igrad, u, ug, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 2
        integer, parameter    :: natoms  = 3
        integer, parameter    :: n3atom  = 75
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: u(nstates,nstates)
        real(wp), intent(out) :: ug(nstates,nstates,natoms,3)
        real(wp)              :: beta11(3)
        real(wp)              :: beta22(3)
        real(wp)              :: betaa11(3)
        real(wp)              :: betaa22(3)
        real(wp)              :: input_cartesian(natom,3)
        real(wp)              :: cartesian_coordinates(natom,3)
        real(wp)              :: ccon
        real(wp)              :: coul(3)
        real(wp)              :: d11(3)
        real(wp)              :: d22(3)
        real(wp)              :: de11(3)
        real(wp)              :: de22(3)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: deq(3)
        real(wp)              :: derad(3)
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        real(wp)              :: do4z(3)
        real(wp)              :: dv3c(3)
        real(wp)              :: dxdr(3,3)
        real(wp)              :: dz(3,3)
        real(wp)              :: easyac
        real(wp)              :: easyba
        real(wp)              :: easybc
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: ground_state_energy
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: exch(3)
        real(wp)              :: ezero(isurf+1)
        integer               :: coordinate_representation
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        integer               :: atom_count
        real(wp)              :: morse_exponentials(3)
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: op3z(3)
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: r(n3atom)
        real(wp)              :: r211
        real(wp)              :: r222
        real(wp)              :: rad
        real(wp)              :: re11(3)
        real(wp)              :: re22(3)
        character(len=75)     :: ref(5)
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
        !initialize
        beta11                            = 0.e0_wp; beta22 = 0.e0_wp; betaa11 = 0.e0_wp; betaa22 = 0.e0_wp; input_cartesian = 0.e0_wp; cartesian_coordinates = 0.e0_wp
        ccon                              = 0.e0_wp; coul = 0.e0_wp; d11 = 0.e0_wp; d22 = 0.e0_wp; de11 = 0.e0_wp
        de22                              = 0.e0_wp; excited_state_internal_gradients = 0.e0_wp; ground_state_internal_gradient = 0.e0_wp; diabatic_internal_gradients = 0.e0_wp; deq = 0.e0_wp; derad = 0.e0_wp; output_excited_state_gradients = 0.e0_wp
        excited_state_cartesian_gradients = 0.e0_wp; output_ground_state_gradient = 0.e0_wp; ground_state_cartesian_gradient = 0.e0_wp; output_diabatic_gradients = 0.e0_wp; diabatic_cartesian_gradients = 0.e0_wp; do4z = 0.e0_wp; dv3c = 0.e0_wp
        dxdr                              = 0.e0_wp; dz = 0.e0_wp; easyac = 0.e0_wp; easyba = 0.e0_wp; easybc = 0.e0_wp; excited_state_energies = 0.e0_wp; ground_state_energy = 0.e0_wp; diabatic_couplings = 0.e0_wp
        exch                              = 0.e0_wp; ezero = 0.e0_wp; coordinate_representation = 0; last_requested_excited_state = 0; last_requested_coupling = 0; first_requested_excited_state = 0
        first_requested_coupling          = 0; atom_count = 0; morse_exponentials = 0.e0_wp; maximum_derivative_order = 0; maximum_excited_state = 0; surface_requests = 0; requested_derivative_order = 0
        op3z                              = 0.e0_wp; output_excited_state_energies = 0.e0_wp; output_ground_state_energy = 0.e0_wp; output_diabatic_couplings = 0.e0_wp; r = 0.e0_wp; r211 = 0.e0_wp; r222 = 0.e0_wp; rad = 0.e0_wp
        re11                              = 0.e0_wp; re22 = 0.e0_wp
        ref                               = ' '
        top3z                             = 0.e0_wp
        tzp3                              = 0.e0_wp
        v3c                               = 0.e0_wp
        v3cp1                             = 0.e0_wp
        v3cp2                             = 0.e0_wp
        v3cp3                             = 0.e0_wp
        v3cp4                             = 0.e0_wp
        v3cp5                             = 0.e0_wp
        v3cp6                             = 0.e0_wp
        z                                 = 0.e0_wp
        zhbrp1                            = 0.e0_wp
        zhbrp2                            = 0.e0_wp
        zhhp1                             = 0.e0_wp
        zhhp2                             = 0.e0_wp
        zp3                               = 0.e0_wp
        zpo                               = 0.e0_wp
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 0
        coordinate_representation  = 3
        maximum_excited_state      = 1
        maximum_derivative_order   = 0
        atom_count                 = 3
        d11                        = [90.4473e0_wp, 109.559e0_wp, 90.4473e0_wp]
        re11                       = [1.41443e0_wp, 0.74144e0_wp, 1.41443e0_wp]
        betaa11                    = [1.81088e0_wp, 1.94198e0_wp, 1.81088e0_wp]
        r211                       = 1.41421356e0_wp
        d22                        = [90.4473e0_wp, 109.559e0_wp, 90.4473e0_wp]
        re22                       = [1.41443e0_wp, 0.74144e0_wp, 1.41443e0_wp]
        betaa22                    = [1.81088e0_wp, 1.94198e0_wp, 1.81088e0_wp]
        r222                       = 1.41421356e0_wp
        ccon                       = 0.5e0_wp
        zhhp1                      = 2.91230499e-1_wp
        zhhp2                      = 4.99160974e0_wp
        zhbrp1                     = 3.84262775e-3_wp
        zhbrp2                     = 2.87642885e-3_wp
        v3cp1                      = 1.67457884e0_wp
        v3cp2                      = 7.09518875e-1_wp
        v3cp3                      = 1.87582679e-1_wp
        v3cp4                      = 2.43517074e-1_wp
        v3cp5                      = 2.04145206e-1_wp
        v3cp6                      = 2.57833022e0_wp
        u                          = 0.e0_wp; ug = 0.e0_wp
        surface_requests           = 1
        input_cartesian            = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                input_cartesian(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        requested_derivative_order = igrad
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            u(1,1) = ground_state_energy*hartree_to_electronvolt
            u(2,2) = excited_state_energies(1)*hartree_to_electronvolt
            u(1,2) = diabatic_couplings(1)*hartree_to_electronvolt
            u(2,1) = diabatic_couplings(1)*hartree_to_electronvolt
        else
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives are unavailable')
        endif

        contains

        !> @brief Evaluate prepare_local_context for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine prepare_local_context
            implicit none
            integer :: i
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            call prepare_local_context_11
            call prepare_local_context_12
            call prepare_local_context_22
            ezero(1) = de11(2)
            ezero(2) = de22(2)
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1) = 'No Reference'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine pot
            implicit none
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            if( surface_requests(1,1) .ne. 0 ) call pot11
            if( surface_requests(2,2) .ne. 0 ) call pot22
            if( surface_requests(1,2)+surface_requests(2,1) .ne. 0 ) call pot12
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate prepare_local_context_11 for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
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
            !        surface_requests - which electronic states are available
            !        requested_derivative_order  = 0 => no derivatives should be calculated
            !        requested_derivative_order  = 1 => calculate first derivatives
            !        NFLAG  - these 20 integer values can be used to flag options
            !                within the potential;
            !   IPRT - FORTRAN unit number to which all output will be written
            !   Potential parameters' default settings
            !                  Variable            Default value
            !                  surface_requests(1,1)         1
            !                  requested_derivative_order                1
            !                  NFLAG(18)           6
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            real(wp), parameter :: cangau = 0.52917706e0_wp
            real(wp), parameter :: charkc = hartree_to_kilocalorie_per_mole
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

        !> @brief Evaluate POT11 for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine pot11
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a1, a2, a3, a4, a5, a6, b(3), capr, capr2, coschi, costh, dist, fcn4, pe, r132p2, r1m3, &
                & r1m3p2, r1p2, r1p3, r1p3m2, r2i, r2p2, r3p2, rr(3), schip2, sintp2, sintp6, sintp8, zchi
            integer             :: i
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
                morse_exponentials(i) = exp(-beta11(i)*(rr(i)-re11(i)))
                coul(i)               = do4z(i) * (zp3(i)*morse_exponentials(i) - top3z(i)) * morse_exponentials(i)
                exch(i)               = do4z(i) * (op3z(i)*morse_exponentials(i) - tzp3(i)) * morse_exponentials(i)
                pe                    = pe + coul(i)
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
            pe                  = pe + v3c
            ground_state_energy = pe+xlamb11/3.e0_wp
        end subroutine pot11

        !> @brief Evaluate prepare_local_context_12 for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine prepare_local_context_12
            implicit none
        end subroutine prepare_local_context_12

        !> @brief Evaluate POT12 for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine pot12
            use chempotpy_kinds, only: wp
            implicit none
            real(wp), parameter :: xlamb12 = 0.4568e0_wp/27.21161e0_wp
            diabatic_couplings(1) = -sqrt(2.e0_wp)*xlamb12/3.e0_wp
        end subroutine pot12

        !> @brief Evaluate prepare_local_context_22 for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine prepare_local_context_22
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            real(wp), parameter :: cangau = 0.52917706e0_wp
            real(wp), parameter :: charkc = hartree_to_kilocalorie_per_mole
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

        !> @brief Evaluate POT22 for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine pot22
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: pe, rr(3), vanti1, vanti3, vmorse, y1, y2, y3
            real(wp), parameter :: xlamb22 = 0.4568e0_wp/27.21161e0_wp
            rr(1)                     = r(1)
            rr(2)                     = r(2)
            rr(3)                     = r(3)
            y2                        = exp(-beta22(2)*(rr(2)-re22(2)))
            y1                        = exp(-beta22(1)*(rr(1)-re22(1)))
            y3                        = exp(-beta22(3)*(rr(3)-re22(3)))
            vmorse                    = de22(2)*(y2**2-2.e0_wp*y2)
            vanti1                    = de22(1)*(y1**2+2.e0_wp*y1)
            vanti3                    = de22(3)*(y3**2+2.e0_wp*y3)
            pe                        = ezero(2)  + vmorse+(vanti1+vanti3)*ccon+2.e0_wp*xlamb22/3.e0_wp
            excited_state_energies(1) = pe
        end subroutine pot22

        !> @brief Evaluate CARTTOR for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(cartesian_coordinates, atom_count, coordinate_representation, r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine eunitzero
            implicit none
            integer :: i, j
            output_ground_state_energy = ground_state_energy
            if( first_requested_excited_state .ne. 0 )then
                do i=first_requested_excited_state, last_requested_excited_state
                    output_excited_state_energies(i) = excited_state_energies(i)
                enddo
            endif
            if( first_requested_coupling .ne. 0 )then
                do j=first_requested_coupling, last_requested_coupling
                    output_diabatic_couplings(j) = diabatic_couplings(j)
                enddo
            endif
        end subroutine eunitzero

        !> @brief Evaluate RTOCART for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
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
            if( coordinate_representation .eq. 1 )then
                do i = 1, atom_count
                    ind                                  = 3*i-2
                    ground_state_cartesian_gradient(i,1) = ground_state_internal_gradient(ind)
                    ground_state_cartesian_gradient(i,2) = ground_state_internal_gradient(ind+1)
                    ground_state_cartesian_gradient(i,3) = ground_state_internal_gradient(ind+2)
                    if( first_requested_excited_state .ne. 0 )then
                        do j = first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(i,1,j) = excited_state_internal_gradients(ind,j)
                            excited_state_cartesian_gradients(i,2,j) = excited_state_internal_gradients(ind+1,j)
                            excited_state_cartesian_gradients(i,3,j) = excited_state_internal_gradients(ind+2,j)
                        enddo
                    endif
                    if( last_requested_coupling .ne. 0 )then
                        do k = first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(i,1,k) = diabatic_internal_gradients(ind,k)
                            diabatic_cartesian_gradients(i,2,k) = diabatic_internal_gradients(ind+1,k)
                            diabatic_cartesian_gradients(i,3,k) = diabatic_internal_gradients(ind+2,k)
                        enddo
                    endif
                enddo
            elseif( coordinate_representation .eq. 2 )then
                do i = 1, atom_count
                    ground_state_cartesian_gradient(i,1) = 0.e0_wp
                    ground_state_cartesian_gradient(i,2) = 0.e0_wp
                    ground_state_cartesian_gradient(i,3) = 0.e0_wp
                    if( first_requested_excited_state .ne. 0 )then
                        do j1=first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(i,1,j1) = 0.e0_wp
                            excited_state_cartesian_gradients(i,2,j1) = 0.e0_wp
                            excited_state_cartesian_gradients(i,3,j1) = 0.e0_wp
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do j2=first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(i,1,j2) = 0.e0_wp
                            diabatic_cartesian_gradients(i,2,j2) = 0.e0_wp
                            diabatic_cartesian_gradients(i,3,j2) = 0.e0_wp
                        enddo
                    endif
                    do j = 1, atom_count
                        if( j .lt. i )then
                            m1 = atom_count*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
                            m1 = atom_count*(i-1) - (i*(i-1))/2 + j-i
                        else
                            cycle
                        endif
                        y = ground_state_internal_gradient(m1)
                        termx                                = (cartesian_coordinates(i,1)-cartesian_coordinates(j, &
                            & 1))/r(m1)
                        termy                                = (cartesian_coordinates(i,2)-cartesian_coordinates(j, &
                            & 2))/r(m1)
                        termz                                = (cartesian_coordinates(i,3)-cartesian_coordinates(j, &
                            & 3))/r(m1)
                        ground_state_cartesian_gradient(i,1) = ground_state_cartesian_gradient(i,1) + termx*y
                        ground_state_cartesian_gradient(i,2) = ground_state_cartesian_gradient(i,2) + termy*y
                        ground_state_cartesian_gradient(i,3) = ground_state_cartesian_gradient(i,3) + termz*y
                        if( first_requested_excited_state .gt. 0 )then
                            y = excited_state_internal_gradients(m1,j1)
                            do j1=first_requested_excited_state, last_requested_excited_state
                                excited_state_cartesian_gradients(i,1,j1) = excited_state_cartesian_gradients(i,1, &
                                    & j1) + termx*y
                                excited_state_cartesian_gradients(i,2,j1) = excited_state_cartesian_gradients(i,2, &
                                    & j1) + termy*y
                                excited_state_cartesian_gradients(i,3,j1) = excited_state_cartesian_gradients(i,3, &
                                    & j1) + termz*y
                            enddo
                        elseif( first_requested_coupling .gt. 0 )then
                            do j2=first_requested_coupling, last_requested_coupling
                                y                                    = diabatic_internal_gradients(m1,j2)
                                diabatic_cartesian_gradients(i,1,j2) = diabatic_cartesian_gradients(i,1,j2) + termx*y
                                diabatic_cartesian_gradients(i,2,j2) = diabatic_cartesian_gradients(i,2,j2) + termy*y
                                diabatic_cartesian_gradients(i,3,j2) = diabatic_cartesian_gradients(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            elseif( coordinate_representation .eq. 3 )then
                do i = 1, atom_count
                    ygs(i) = ground_state_internal_gradient(i)/r(i)
                    if( first_requested_excited_state .ne. 0 )then
                        do j=first_requested_excited_state, last_requested_excited_state
                            yes(i, j) = excited_state_internal_gradients(i,j)/r(i)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do k=first_requested_coupling, last_requested_coupling
                            yij(i, k) = diabatic_internal_gradients(i,k)/r(i)
                        enddo
                    endif
                enddo
                do k = 1, 3
                    term12                               = cartesian_coordinates(1,k)-cartesian_coordinates(2,k)
                    term23                               = cartesian_coordinates(2,k)-cartesian_coordinates(3,k)
                    term13                               = cartesian_coordinates(1,k)-cartesian_coordinates(3,k)
                    ground_state_cartesian_gradient(1,k) = term12*ygs(1) + term13*ygs(3)
                    ground_state_cartesian_gradient(2,k) = -term12*ygs(1) + term23*ygs(2)
                    ground_state_cartesian_gradient(3,k) = -term13*ygs(3) - term23*ygs(2)
                    if( first_requested_excited_state .ne. 0 )then
                        do j1=first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(1,k,j1) = term12*yes(1, j1) + term13*yes(3, j1)
                            excited_state_cartesian_gradients(2,k,j1) = -term12*yes(1, j1) + term23*yes(2, j1)
                            excited_state_cartesian_gradients(3,k,j1) = -term13*yes(3, j1) - term23*yes(2, j1)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do j2=first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(1,k,j2) = term12*yij(1, j2) + term13*yij(3, j2)
                            diabatic_cartesian_gradients(2,k,j2) = -term12*yij(1, j2) + term23*yij(2, j2)
                            diabatic_cartesian_gradients(3,k,j2) = -term13*yij(3, j2) - term23*yij(2, j2)
                        enddo
                    endif
                enddo
            elseif( coordinate_representation .eq. 4 )then
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

        !> @brief Evaluate DEDCOU for the BrH2_LEPS_LTBZSO_1995_DPEM potential.
        subroutine dedcou
            implicit none
            integer :: i, j, k
            do i = 1, atom_count
                output_ground_state_gradient(i,1) = ground_state_cartesian_gradient(i,1)
                output_ground_state_gradient(i,2) = ground_state_cartesian_gradient(i,2)
                output_ground_state_gradient(i,3) = ground_state_cartesian_gradient(i,3)
                if( first_requested_excited_state .ne. 0 )then
                    do j=first_requested_excited_state, last_requested_excited_state
                        output_excited_state_gradients(i,1,j) = excited_state_cartesian_gradients(i,1,j)
                        output_excited_state_gradients(i,2,j) = excited_state_cartesian_gradients(i,2,j)
                        output_excited_state_gradients(i,3,j) = excited_state_cartesian_gradients(i,3,j)
                    enddo
                endif
                if( first_requested_coupling .ne. 0 )then
                    do k=first_requested_coupling, last_requested_coupling
                        output_diabatic_gradients(i,1,k) = diabatic_cartesian_gradients(i,1,k)
                        output_diabatic_gradients(i,2,k) = diabatic_cartesian_gradients(i,2,k)
                        output_diabatic_gradients(i,3,k) = diabatic_cartesian_gradients(i,3,k)
                    enddo
                endif
            enddo
        end subroutine dedcou
    end subroutine evaluate
end module chempotpy_brh2_leps_ltbzso_1995_dpem_pes

!> @brief Expose BrH2_LEPS_LTBZSO_1995 through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_brh2_leps_ltbzso_1995_dpem_pes, only: evaluate
    use chempotpy_kinds,                          only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: u(nstates,nstates)
    real(wp), intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem
