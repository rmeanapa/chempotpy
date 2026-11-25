! System:          HI2
! Functional form: eLEPS (extended London-Eyring-Polanyi-Sato)
! Common name:     KK
! Reference:       J. A. Kaye and A. Kuppermann, Chem. Phys. Lett.
!                  77, 573 (1981).
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
!           R(1) = R(I-H)
!           R(2) = R(H-I')
!           R(3) = R(I-I')
! The zero of energy for this potential occurs when the I atom is
! "infinitely" far from the HI molecule, and the HI molecule is at
! its equilibrium separation.
! The array ground_state_internal_gradient contains the analytic first derivatives of the potent
! energy with respect to the internal coordinates if the variable requested_derivative_order i
! equal to 1 (see not below).  The derivatives are defined as follows:
!           ground_state_internal_gradient(1) = dV(R)/dR(1)
!           ground_state_internal_gradient(2) = dV(R)/dR(2)
!           ground_state_internal_gradient(3) = dV(R)/dR(3)
! where V(R) = V { R(1), R(2), R(3) }.
! The variable requested_derivative_order controls whether derivatives are calculated.
! If requested_derivative_order = 1 first derivatives are calculated, otherwise no derivatives
! are calculated.  Note:  the default value of requested_derivative_order is 1.
! The array NFLAG is unused in this potential surface.
! The common block ASYCM is defined as follows:
!           EASYAB = DE(I-H)
!           EASYBC = DE(H-I')
!           EASYAC = DE(I-I')
! where DE is the dissociation energy in atomic units.
!   Potential parameters' default settings
!                  Variable            Default value
!                  requested_derivative_order                1
!                  NFLAG(18)           6

!> @brief Implement I2H_LEPS_KK_1981, a potential energy surface for I2H.
module chempotpy_i2h_leps_kk_1981_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the I2H_LEPS_KK_1981 potential energy surface for I2H.
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
        real(wp)              :: b(3)
        real(wp)              :: input_cartesian(natom,3)
        real(wp)              :: cartesian_coordinates(natom,3)
        real(wp)              :: de(3)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        real(wp)              :: easyab
        real(wp)              :: easyac
        real(wp)              :: easybc
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: ground_state_energy
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: ezero(isurf+1)
        real(wp)              :: hde(3)
        real(wp)              :: hpf(3)
        integer               :: coordinate_representation
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        integer               :: atom_count
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: pf
        real(wp)              :: pfs(3)
        real(wp)              :: pft(3)
        real(wp)              :: r(n3atom)
        real(wp)              :: re(3)
        character(len=75)     :: ref(5)
        real(wp)              :: s(3)
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        b                                = 0.e0_wp; input_cartesian = 0.e0_wp; cartesian_coordinates = 0.e0_wp; de = 0.e0_wp
        excited_state_internal_gradients = 0.e0_wp; ground_state_internal_gradient = 0.e0_wp; diabatic_internal_gradients = 0.e0_wp; output_excited_state_gradients = 0.e0_wp; excited_state_cartesian_gradients = 0.e0_wp; output_ground_state_gradient = 0.e0_wp; ground_state_cartesian_gradient = 0.e0_wp
        output_diabatic_gradients        = 0.e0_wp; diabatic_cartesian_gradients = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp; excited_state_energies = 0.e0_wp; ground_state_energy = 0.e0_wp
        diabatic_couplings               = 0.e0_wp; ezero = 0.e0_wp; hde = 0.e0_wp; hpf = 0.e0_wp; coordinate_representation = 0
        last_requested_excited_state     = 0; last_requested_coupling = 0; first_requested_excited_state = 0; first_requested_coupling = 0; atom_count = 0; maximum_derivative_order = 0; maximum_excited_state = 0; surface_requests = 0; requested_derivative_order = 0
        output_excited_state_energies    = 0.e0_wp; output_ground_state_energy = 0.e0_wp; output_diabatic_couplings = 0.e0_wp; pf = 0.e0_wp; pfs = 0.e0_wp; pft = 0.e0_wp
        r                                = 0.e0_wp; re = 0.e0_wp
        ref                              = ' '
        s                                = 0.e0_wp
        de                               = [3.3303e0_wp, 3.3303e0_wp, 1.5567e0_wp]
        re                               = [3.0236e0_wp, 3.0236e0_wp, 5.0457e0_wp]
        b                                = [0.9260e0_wp, 0.9260e0_wp, 0.9843e0_wp]
        s                                = [0.2000e0_wp, 0.2000e0_wp, 0.1250e0_wp]
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 0
        coordinate_representation  = 3
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 3
        v                          = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
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
            do istate=1, nstates
                p(istate) = output_ground_state_energy*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = output_ground_state_energy*hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = output_ground_state_gradient(iatom,idir)*hartree_per_bohr_to_ev_per_angstrom
                enddo
            enddo
        elseif( igrad == 2 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives above first order are unavailable')
        endif

        contains

        !> @brief Evaluate prepare_local_context for the I2H_LEPS_KK_1981 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            real(wp), parameter :: cevau = 27.2113961e0_wp
            ! Write header potential parameters to UNIT NFLAG(18).
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            ! Calculate first derivatives (default).  If no first derivatives are
            ! desired, requested_derivative_order should be set equal to zero after prepare_local_context is called.
            ! Convert to atomic units and calculate some useful constants.
            do i=1, 3
                de(i)  = de(i)/cevau
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
            ref(1) = 'J. A. Kaye and A. Kuppermann,'
            ref(2) = 'Chem. Phys. Lett. 77, 573(1981)'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the I2H_LEPS_KK_1981 potential.
        subroutine pot
            ! The coordinates are defined as follows:
            !           R(1) = R(I-H)
            !           R(2) = R(H-I')
            !           R(3) = R(I-I')
            ! The zero of energy for this potential occurs when the I atom is
            ! "infinitely" far from the HI molecule, and the HI molecule is at
            ! its equilibrium separation.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dj1dr1, dj2dr2, dj3dr3, dq1dr1, dq2dr2, dq3dr3, hds1dr1, hds2dr2, hds3dr3, hdt1dr1, hdt2dr2, hdt3dr3, hradi, hs1, hs2, hs3, ht1, ht2, ht3, rad, radi, wj1, wj2, wj3, wq1, wq2, wq3, y1, y2, y3
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
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
            rad                 = sqrt(wj1*(wj1-wj2)+wj2*(wj2-wj3)+wj3*(wj3-wj1))
            ground_state_energy = ezero(1)  +wq1+wq2+wq3-rad
            !  compute first derivatives if necessary
            if( requested_derivative_order .eq. 1 )then
                radi                              = 1.e0_wp/rad
                hradi                             = 0.5e0_wp*radi
                hds1dr1                           = pfs(1)*(y1-1.e0_wp)*y1
                hds2dr2                           = pfs(2)*(y2-1.e0_wp)*y2
                hds3dr3                           = pfs(3)*(y3-1.e0_wp)*y3
                hdt1dr1                           = pft(1)*(y1+1.e0_wp)*y1
                hdt2dr2                           = pft(2)*(y2+1.e0_wp)*y2
                hdt3dr3                           = pft(3)*(y3+1.e0_wp)*y3
                dq1dr1                            = hds1dr1+hdt1dr1
                dq2dr2                            = hds2dr2+hdt2dr2
                dq3dr3                            = hds3dr3+hdt3dr3
                dj1dr1                            = hds1dr1-hdt1dr1
                dj2dr2                            = hds2dr2-hdt2dr2
                dj3dr3                            = hds3dr3-hdt3dr3
                ground_state_internal_gradient(1) = dq1dr1-hradi*dj1dr1*(2.e0_wp*wj1-wj2-wj3)
                ground_state_internal_gradient(2) = dq2dr2-hradi*dj2dr2*(2.e0_wp*wj2-wj1-wj3)
                ground_state_internal_gradient(3) = dq3dr3-hradi*dj3dr3*(2.e0_wp*wj3-wj1-wj2)
            endif
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate CARTTOR for the I2H_LEPS_KK_1981 potential.
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

        !> @brief Evaluate EUNITZERO for the I2H_LEPS_KK_1981 potential.
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

        !> @brief Evaluate RTOCART for the I2H_LEPS_KK_1981 potential.
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

        !> @brief Evaluate DEDCOU for the I2H_LEPS_KK_1981 potential.
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
end module chempotpy_i2h_leps_kk_1981_pes

!> @brief Expose I2H_LEPS_KK_1981 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_i2h_leps_kk_1981_pes, only: evaluate
    use chempotpy_kinds,                only: wp
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
