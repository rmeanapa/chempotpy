!   System:          BrH2
!   Functional form: extended-LEPS plus three-center term
!   Common name:     BrH2SEC
!   Interface:       potlib2001
!   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g.
!                    J. Phys. Chem. 99, 207 (1995)
!   Number of bodies:              3
!   Number of derivatives:         1
!   Number of electronic surfaces: 1
!   Calling Sequence:
!      prepare_local_context - initializes the potential's variables and
!               must be called once before any calls to POT
!      POT    - driver for the evaluation of the energy and the derivati
!               of the energy with respect to the coordinates for a give
!               geometric configuration
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      ground electronic state
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Br
!      infinitely far from the H2 diatomic and R(H2) set equal to the
!      H2 equilibrium diatomic value.
!   Parameters:
!      Set in the BLOCK DATA subprogram PTPACM
!   Coordinates:
!      Internal, Definition: R(1) = R(Br-H)
!                            R(2) = R(H-H)
!                            R(3) = R(Br-H)
!   Common Blocks (used between the calling program and this potential):
!        passes the coordinates, ground state electronic energy, and
!        derivatives of the ground electronic state energy with respect
!        to the coordinates.
!        passes the control flags where
!          NSURF  - not used
!        requested_derivative_order  = 0 => no derivatives are computed
!              = 1 => derivatives of the energy for the ground electroni
!                     state with respect to the coordinates are computed
!          NFLAG  - Control flags
!      NFLAG(18-20)
!        passes the FORTRAN unit number used for potential output
!      /ASYCM/ EASYAB, EASYBC, EASYAC
!        passes the energy in the three asymptotic valleys for an A + BC
!        The energy in the AB valley, EASYAB, is equal to the energy of
!        C atom "infinitely" far from the AB diatomic and R(AB) set equa
!        Re(AB), the equilibrium bond length for the AB diatomic.
!        In this potential the AB valley represents H infinitely far fro
!        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the t
!        EASYBC and EASYAC represent the energies in the H2 and the BrH
!        valleys, respectively.
!   Default Parameter Values:
!      Variable      Default value
!      requested_derivative_order             1
!      NFLAG(18)        6

!> @brief Implement BrH2_LEPS_LTBZ1_1995, a potential energy surface for BrH2.
module chempotpy_brh2_leps_ltbz1_1995_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the BrH2_LEPS_LTBZ1_1995 potential energy surface for BrH2.
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
        real(wp)              :: beta(3)
        real(wp)              :: input_cartesian(natom,3)
        real(wp)              :: cartesian_coordinates(natom,3)
        real(wp)              :: coul(3)
        real(wp)              :: de(3)
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
        real(wp)              :: easyab
        real(wp)              :: easyac
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
        real(wp)              :: rad
        real(wp)              :: re(3)
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
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        beta                              = 0.e0_wp; input_cartesian = 0.e0_wp; cartesian_coordinates = 0.e0_wp
        coul                              = 0.e0_wp; de = 0.e0_wp; excited_state_internal_gradients = 0.e0_wp; ground_state_internal_gradient = 0.e0_wp; diabatic_internal_gradients = 0.e0_wp; deq = 0.e0_wp; derad = 0.e0_wp; output_excited_state_gradients = 0.e0_wp
        excited_state_cartesian_gradients = 0.e0_wp; output_ground_state_gradient = 0.e0_wp; ground_state_cartesian_gradient = 0.e0_wp; output_diabatic_gradients = 0.e0_wp; diabatic_cartesian_gradients = 0.e0_wp; do4z = 0.e0_wp; dv3c = 0.e0_wp
        dxdr                              = 0.e0_wp; dz = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp; excited_state_energies = 0.e0_wp; ground_state_energy = 0.e0_wp; diabatic_couplings = 0.e0_wp
        exch                              = 0.e0_wp; ezero = 0.e0_wp; coordinate_representation = 0; last_requested_excited_state = 0; last_requested_coupling = 0; first_requested_excited_state = 0
        first_requested_coupling          = 0; atom_count = 0; morse_exponentials = 0.e0_wp; maximum_derivative_order = 0; maximum_excited_state = 0; surface_requests = 0; requested_derivative_order = 0
        op3z                              = 0.e0_wp; output_excited_state_energies = 0.e0_wp; output_ground_state_energy = 0.e0_wp; output_diabatic_couplings = 0.e0_wp; r = 0.e0_wp; rad = 0.e0_wp; re = 0.e0_wp
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
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 3
        de                         = [90.4473e0_wp, 109.559e0_wp, 90.4473e0_wp]
        re                         = [1.41443e0_wp, 0.74144e0_wp, 1.41443e0_wp]
        beta                       = [1.81088e0_wp, 1.94198e0_wp, 1.81088e0_wp]
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

        !> @brief Evaluate prepare_local_context for the BrH2_LEPS_LTBZ1_1995 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            real(wp), parameter :: cangau = 0.52917706e0_wp
            real(wp), parameter :: charkc = hartree_to_kilocalorie_per_mole
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            !   Convert the constants to atomic units
            do i = 1, 3
                de(i)   = de(i) / charkc
                re(i)   = re(i) / cangau
                beta(i) = beta(i) * cangau
            enddo
            !    Initialize the energy in the three asymptotic valleys
            easyab   = de(1)
            easybc   = de(2)
            easyac   = de(3)
            ezero(1) = de(2)
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1) = 'G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g. Zhao'
            ref(2) = 'J. Phys. Chem. 99, 207 (1995)'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the BrH2_LEPS_LTBZ1_1995 potential.
        subroutine pot
            !   The potential parameters must be passed through the COMMON BLOCKS
            !   PT1CM, ASYCM, PT2CM,  All information passed through the
            !   COMMON BLOCKS PT1CM, ASYCM must be in Hartree atomic units.
            !   Coordinates:
            !      Internal, Definition: R(1) = R(Br-H)
            !                            R(2) = R(H-H)
            !                            R(3) = R(Br-H)
            !        The energy in the AB valley, EASYAB, is equal to the energy of
            !        C atom "infinitely" far from the AB diatomic and R(AB) set equa
            !        Re(AB), the equilibrium bond length for the AB diatomic.
            !        In this potential the AB valley represents H infinitely far fro
            !        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the t
            !        EASYBC and EASYAC represent the energies in the H2 and the HBr
            !        valleys, respectively.
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            real(wp), parameter :: r2 = 1.41421356e0_wp
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   Check the value of requested_derivative_order
            !   Compute bond dependent Sato parameters
            call sato
            !   Compute the LEPS constants
            do i = 1, 3
                zpo(i)   = 1.0e0_wp+z(i)
                op3z(i)  = 1.0e0_wp+3.0e0_wp*z(i)
                top3z(i) = 2.0e0_wp*op3z(i)
                zp3(i)   = z(i)+3.0e0_wp
                tzp3(i)  = 2.0e0_wp*zp3(i)
                do4z(i)  = de(i)/4.0e0_wp/zpo(i)
            enddo
            ground_state_energy = 0.e0_wp
            do i = 1, 3
                morse_exponentials(i) = exp(-beta(i)*(r(i)-re(i)))
                coul(i)               = do4z(i) * (zp3(i)*morse_exponentials(i) - top3z(i)) * morse_exponentials(i)
                exch(i)               = do4z(i) * (op3z(i)*morse_exponentials(i) - tzp3(i)) * morse_exponentials(i)
                ground_state_energy   = ground_state_energy + coul(i)
            enddo
            rad                 = sqrt((exch(1) - exch(2)) ** 2 + (exch(2) - exch(3)) ** 2 + (exch(3) - exch(1)) ** 2 )
            ground_state_energy = ground_state_energy - (rad / r2) + ezero(1)
            !   Compute the three-center term
            call v3cfcn
            !   Sum all the terms that make up the energy
            ground_state_energy = ground_state_energy + v3c
            !   Compute the derivatives if requested_derivative_order = 1
            goto_block_700_419: block
                if( requested_derivative_order .ne. 1 ) exit goto_block_700_419
                !   DCOUL computes the derivatives of the coulombic terms w.r.t. R(I)
                !   DEXCH computes the derivatives of the exchange terms w.r.t. R(I) and
                !         also computes the derivative of the term RAD w.r.t. R(i)
                !   DCOUL must be called before DEXCH becaues this subprogram computes
                !   the DXDR terms needed by both subprograms.
                !   Compute DEQ = dE/dQ * [dQ/dX * dX/dR + dQ/dZ * dZ/dR]
                call dcoul
                !   Compute DERAD = dRAD/dJ * [dJ/dX * dX/dR + dJ/dZ * dZ/dR]
                call dexch
                !   Put all the partial derivatives together:
                do i = 1, 3
                    ground_state_internal_gradient(i) = deq(i) - derad(i)/r2
                enddo
                !   Add in the derivatives of V3C w.r.t. R(I)
                do i = 1, 3
                    ground_state_internal_gradient(i) = ground_state_internal_gradient(i) + dv3c(i)
                enddo
            end block goto_block_700_419
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate SATO for the BrH2_LEPS_LTBZ1_1995 potential.
        subroutine sato
            !   This subroutine determines the Sato parameters as a function of the
            !   bond lengths.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: capr, capr2, capri, coschi, dcapr1, dcapr2, dcapr3, dchi1, dchi2, dchi3, dzchic, dzchir, r1p2, r2i, r2p2, r3p2, schip2, tmp, tmp1, tmp2, zchi
            integer  :: i, j
            do i = 1, 3
                do j = 1, 3
                    dz(i,j) = 0.0e0_wp
                enddo
            enddo
            r1p2 = r(1) * r(1)
            r2p2 = r(2) * r(2)
            r3p2 = r(3) * r(3)
            r2i  = 1.0e0_wp/r(2)
            zchi = 0.0e0_wp
            !   CAPR IS THE Br TO H2 DISTANCE
            capr2 = 0.5e0_wp * (r1p2 + r3p2 - 0.5e0_wp * r2p2)
            if( capr2 .le. 0.e0_wp )then
                capr2 = 0.0e0_wp
                capr  = 0.0e0_wp
            else
                capr = sqrt(capr2)
                !   Compute Chi, the angle between CAPR and RHH
                coschi = 0.5e0_wp*(r3p2 - r1p2)/(r(2)*capr)
                schip2 = 1.0e0_wp - coschi*coschi
                zchi   = capr2*(zhbrp1*schip2 + zhbrp2*schip2*schip2)
            endif
            z(1) = 0.1675e0_wp + zchi
            z(3) = z(1)
            z(2) = zhhp1*r2p2/(zhhp2+r2p2)
            !   Compute the derivative of the Sato parameters with respect to
            !   bond lengths.
            !   The HH Sato parameter depends on R(2) only
            goto_block_700_527: block
                if( requested_derivative_order .ne. 1 ) exit goto_block_700_527
                tmp     = r(2)/(zhhp2 + r2p2)
                dz(2,2) = 2.0e0_wp*z(2)*(r2i - tmp)
                !   Compute the derivative of ZHBr with respect R(i) if CHI is
                !   is between 0 and 90 degrees.
                goto_block_100_526: block
                    if( capr .eq. 0.e0_wp ) exit goto_block_100_526
                    capri = 1.0e0_wp/capr
                    !   Compute the derivative of R(Br-HH) with respect to R(i)
                    dcapr1 = 0.5e0_wp*r(1)*capri
                    dcapr2 = -0.25e0_wp*r(2)*capri
                    dcapr3 = 0.5e0_wp*r(3)*capri
                    !   Compute the derivative of CHI with respect to R(i)
                    !   These derivatives a actually sin(CHI)*dCHI/dR(i)
                    tmp1  = coschi*capri
                    tmp2  = 1.0e0_wp/(r(2)*capr)
                    dchi1 = r(1)*tmp2 + tmp1*dcapr1
                    dchi2 = tmp1*(capr + r(2)*dcapr2)*r2i
                    dchi3 = tmp1*dcapr3 - r(3)*tmp2
                    !   Compute the derivative of ZHBr with respect to CAPR and CHI
                    dzchir = 2.0e0_wp*zchi*capri
                    dzchic = 2.0e0_wp*capr2*coschi*(zhbrp1 + 2.0e0_wp*zhbrp2*schip2)
                    !   Compute the derivative of ZHBr with respect to R(i)
                    dz(1,1) = dzchir*dcapr1 + dzchic*dchi1
                    dz(1,2) = dzchir*dcapr2 + dzchic*dchi2
                    dz(1,3) = dzchir*dcapr3 + dzchic*dchi3
                    dz(3,1) = dz(1,1)
                    dz(3,2) = dz(1,2)
                    dz(3,3) = dz(1,3)
                end block goto_block_100_526
            end block goto_block_700_527
        end subroutine sato

        !> @brief Evaluate V3CFCN for the BrH2_LEPS_LTBZ1_1995 potential.
        subroutine v3cfcn
            !   This subroutine computes the energy and the derivatives of the
            !   three-center term for the H-Br-H barrier as a function of the geomet
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a1, a2, a3, a4, a5, a6, costh, dfcn1, dfcn2, dfcn3, dfcn4, dist, dthdr1, dthdr2, dthdr3, fcn1, fcn2, fcn3, fcn4, r132p2, r1m3, r1m3p2, r1p3, r1p3m2, sintp2, sintp6, sintp8, v3cdth
            a1     = v3cp1
            a2     = v3cp2
            a3     = v3cp3
            a4     = v3cp4
            a5     = v3cp5
            a6     = v3cp6
            r1p3   = r(1) + r(3)
            r1m3   = r(1) - r(3)
            r1p3m2 = r1p3 - r(2)
            r1m3p2 = r1m3*r1m3
            r132p2 = r1p3m2*r1p3m2
            !   Compute the bond angle from the three internal coordinates
            dist   = r(1)*r(1) + r(3)*r(3) - r(2)*r(2)
            costh  = dist/(2.0e0_wp*r(1)*r(3))
            sintp2 = 1.0e0_wp - costh*costh
            sintp6 = sintp2*sintp2*sintp2
            sintp8 = sintp6*sintp2
            fcn1   = exp(-a2*r1p3)
            fcn2   = exp(-a3*r1m3p2)
            fcn3   = exp(-a4*r132p2)
            fcn4   = 1.e0_wp + a5*sintp2 + a6*sintp8
            v3c    = a1*fcn1*fcn2*fcn3*fcn4
            !   Compute the derivatives of V3C with respect to the internal coordina
            !   if requested_derivative_order = 1
            goto_block_700_605: block
                if( requested_derivative_order .ne. 1 ) exit goto_block_700_605
                !   Compute the partial derivatives of the FCN terms
                dfcn1 = -a2
                dfcn2 = -2.0e0_wp*a3*r1m3
                dfcn3 = -2.0e0_wp*a4*r1p3m2
                dfcn4 = 2.0e0_wp*a5*costh + 8.0e0_wp*a6*costh*sintp6
                !   Compute the partial derivative of THETA w.r.t. R(i)
                !   DTHDR(i) = SIN(THETA)*DTHDR(I)
                dthdr1  = costh/r(1) - 1.0e0_wp/r(3)
                dthdr2  = r(2)/(r(1)*r(3))
                dthdr3  = costh/r(3) - 1.0e0_wp/r(1)
                v3cdth  = a1*fcn1*fcn2*fcn3*dfcn4
                dv3c(1) = v3c * (dfcn1 + dfcn2 + dfcn3) + v3cdth*dthdr1
                dv3c(2) = -v3c * dfcn3 + v3cdth*dthdr2
                dv3c(3) = v3c * (dfcn1 - dfcn2 + dfcn3) + v3cdth*dthdr3
            end block goto_block_700_605
        end subroutine v3cfcn

        !> @brief Evaluate DCOUL for the BrH2_LEPS_LTBZ1_1995 potential.
        subroutine dcoul
            !    This subprogram computes the contribution to the derivative of
            !    the energy with respect to R(i) from the coulombic terms.
            !    That is:
            !            DE/DR(j) = DE/DQ(i) * DQ(i)/DR(j)
            !    In this subprogram the following convention is assumed for the
            !    numbering of the internal coordinates:
            !    For Br + AB -> ABr + B
            !                   R(1) = R(Br-A)
            !                   R(2) = R(A-B)
            !                   R(3) = R(Br-B)
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dqdx(3), dqdz(3)
            integer  :: i, j
            do i = 1, 3
                deq(i)  = 0.e0_wp
                dqdx(i) = 0.e0_wp
                dqdz(i) = 0.e0_wp
                do j = 1, 3
                    dxdr(i,j) = 0.e0_wp
                enddo
            enddo
            !   First compute the partial derivatives dQ/dX, dQ/dZ, and
            !   the derivative of morse_exponentials w.r.t. R(i), DX/DR.
            do i = 1, 3
                dqdx(i)   = do4z(i)*(tzp3(i)*morse_exponentials(i) - top3z(i))
                dqdz(i)   = do4z(i)*morse_exponentials(i)*(morse_exponentials(i) - 6.0e0_wp) - coul(i)/zpo(i)
                dxdr(i,i) = -beta(i)*morse_exponentials(i)
            enddo
            !   Compute the contribution to the derivate of E w.r.t. R(i) from
            !   the coulombic terms, i.e., dE/dR(i) = dQ/dX * DX/DR(i) + dQ/dZ *DZ/D
            do i = 1, 3
                do j = 1, 3
                    deq(i) = deq(i) + dqdx(j)*dxdr(j,i) + dqdz(j)*dz(j,i)
                enddo
            enddo
        end subroutine dcoul

        !> @brief Evaluate DEXCH for the BrH2_LEPS_LTBZ1_1995 potential.
        subroutine dexch
            !    This subprogram computes the contribution to the derivative of
            !    the energy with respect to R(i) from the exchange terms.
            !    That is:
            !            DE/DR(j) = DE/DJ(i) * DJ(i)/DR(j)
            !    where
            !            E(LEPS)  = Q1 + Q2 + Q3 - RAD/SQRT(2)
            !    In this subprogram the following convention is assumed for the
            !    numbering of the internal coordinates:
            !    For Br + AB -> ABr + B
            !                   R(1) = R(Br-A)
            !                   R(2) = R(A-B)
            !                   R(3) = R(Br-B)
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: djdr(3,3), djdx(3), djdz(3), draddj(3), sumj
            integer  :: i, j
            do i = 1, 3
                derad(i) = 0.e0_wp
                djdx(i)  = 0.e0_wp
                djdz(i)  = 0.e0_wp
            enddo
            sumj = 0.0e0_wp
            !   First compute the partial derivatives dJ/dX, dJ/dZ
            do i = 1, 3
                sumj    = sumj + exch(i)
                djdx(i) = do4z(i)*(top3z(i)*morse_exponentials(i) - tzp3(i))
                djdz(i) = do4z(i)*morse_exponentials(i)*(3.0e0_wp*morse_exponentials(i) - 2.0e0_wp) - exch(i)/zpo(i)
            enddo
            !   Compute the partial derivatives of the exchange term, J, w.r.t. R(i)
            !   The outer loop, I, is over the bond lengths, the inner loop is over
            !   each morse_exponentials and Z term.
            do i = 1, 3
                do j = 1, 3
                    djdr(j,i) = djdx(j)*dxdr(j,i) + djdz(j)*dz(j,i)
                enddo
            enddo
            !   Compute the partial derivatives of the term RAD w.r.t. J(I).
            do i = 1, 3
                draddj(i) = (3.0e0_wp*exch(i) - sumj)/rad
            enddo
            do i = 1, 3
                do j = 1, 3
                    derad(i) = derad(i) + draddj(j)*djdr(j,i)
                enddo
            enddo
        end subroutine dexch

        !> @brief Evaluate CARTTOR for the BrH2_LEPS_LTBZ1_1995 potential.
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

        !> @brief Evaluate EUNITZERO for the BrH2_LEPS_LTBZ1_1995 potential.
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

        !> @brief Evaluate RTOCART for the BrH2_LEPS_LTBZ1_1995 potential.
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

        !> @brief Evaluate DEDCOU for the BrH2_LEPS_LTBZ1_1995 potential.
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
end module chempotpy_brh2_leps_ltbz1_1995_pes

!> @brief Expose BrH2_LEPS_LTBZ1_1995 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_brh2_leps_ltbz1_1995_pes, only: evaluate
    use chempotpy_kinds,                    only: wp
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
