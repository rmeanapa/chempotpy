!   System:          BrH2
!   Functional form: Diatomics-in-molecules plus three-center term
!   Common name:     BrH2 DIM-3C
!   Interface:       potlib2001
!   References:      D. C. Clary
!                    Chem. Phys. 71, 117 (1982)
!                    I. Last and M. Baer
!                    in Potential Energy Surfaces and Dynamics
!                    edited by D. G. Truhlar (Plenum, New York, 1981)
!                    p. 519
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
!      Set in data statements
!   Coordinates:
!      Internal, Definition: R(1) = R(Br-H)
!                            R(2) = R(H-H)
!                            R(3) = R(Br-H)
!   Common Blocks (used between the calling program and this potential):
!        passes the coordinates, ground state electronic energy, and
!        derivatives of the ground electronic state energy with respect
!        to the coordinates.
!      NFLAG(18-20)
!        passes the FORTRAN unit number used for output from the potenti
!      /ASYCM/ EASYAB, EASYBC, EASYAC
!        passes the energy in the three asymptotic valleys for an A + BC
!        The energy in the AB valley, EASYAB, is equal to the energy of
!        C atom "infinitely" far from the AB diatomic and R(AB) set equa
!        Re(AB), the equilibrium bond length for the AB diatomic.
!        In this potential the AB valley represents H infinitely far fro
!        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the t
!        EASYBC and EASYAC represent the energies in the H2 and the HBr
!        valleys, respectively.
!   Default Parameter Values:
!      Variable      Default value
!      NFLAG(18)        6

!> @brief Implement BrH2_LEPS_C_1982, a potential energy surface for BrH2.
module chempotpy_brh2_leps_c_1982_pes
    use chempotpy_linear_algebra, only: symmetric_eigendecomposition
    use chempotpy_units,          only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants,      only: PI
    use chempotpy_kinds,          only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the BrH2_LEPS_C_1982 potential energy surface for BrH2.
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
        real(wp)              :: ahh
        real(wp)              :: ahx
        real(wp)              :: alfw
        real(wp)              :: bhh
        real(wp)              :: bhx
        real(wp)              :: c3
        real(wp)              :: input_cartesian(natom,3)
        real(wp)              :: cartesian_coordinates(natom,3)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: denom
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: dh1(10)
        real(wp)              :: dh2(10)
        real(wp)              :: dh3(10)
        real(wp)              :: dhh
        real(wp)              :: dhx
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        real(wp)              :: e(4)
        real(wp)              :: easyab
        real(wp)              :: easyac
        real(wp)              :: easybc
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: ground_state_energy
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: eps
        real(wp)              :: eta1p
        real(wp)              :: eta3p
        real(wp)              :: eta3s
        real(wp)              :: etahh
        real(wp)              :: ezero(isurf+1)
        real(wp)              :: h(10)
        integer               :: coordinate_representation
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        real(wp)              :: three_center_decay
        integer               :: atom_count
        real(wp)              :: eigenvectors(4,4)
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: one
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: r(n3atom)
        character(len=75)     :: ref(5)
        real(wp)              :: rhh
        real(wp)              :: rhx
        real(wp)              :: rt34
        real(wp)              :: rt38
        real(wp)              :: scr1(4)
        real(wp)              :: scr2(4)
        real(wp)              :: xib
        real(wp)              :: xih
        real(wp)              :: xix
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        ahh                           = 0.e0_wp; ahx = 0.e0_wp; alfw = 0.e0_wp; bhh = 0.e0_wp; bhx = 0.e0_wp; c3 = 0.e0_wp; input_cartesian = 0.e0_wp
        cartesian_coordinates         = 0.e0_wp; excited_state_internal_gradients = 0.e0_wp; ground_state_internal_gradient = 0.e0_wp; diabatic_internal_gradients = 0.e0_wp
        denom                         = 0.e0_wp; output_excited_state_gradients = 0.e0_wp; excited_state_cartesian_gradients = 0.e0_wp; output_ground_state_gradient = 0.e0_wp; ground_state_cartesian_gradient = 0.e0_wp; dh1 = 0.e0_wp; dh2 = 0.e0_wp
        dh3                           = 0.e0_wp; dhh = 0.e0_wp; dhx = 0.e0_wp; output_diabatic_gradients = 0.e0_wp; diabatic_cartesian_gradients = 0.e0_wp; e = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp
        easybc                        = 0.e0_wp; excited_state_energies = 0.e0_wp; ground_state_energy = 0.e0_wp; diabatic_couplings = 0.e0_wp; eps = 0.e0_wp; eta1p = 0.e0_wp; eta3p = 0.e0_wp; eta3s = 0.e0_wp
        etahh                         = 0.e0_wp; ezero = 0.e0_wp; h = 0.e0_wp; coordinate_representation = 0; last_requested_excited_state = 0; last_requested_coupling = 0
        first_requested_excited_state = 0; first_requested_coupling = 0; three_center_decay = 0.e0_wp; atom_count = 0; eigenvectors = 0.e0_wp; maximum_derivative_order = 0; maximum_excited_state = 0; surface_requests = 0
        requested_derivative_order    = 0; one = 0.e0_wp; output_excited_state_energies = 0.e0_wp; output_ground_state_energy = 0.e0_wp; output_diabatic_couplings = 0.e0_wp; r = 0.e0_wp
        ref                           = ' '
        rhh                           = 0.e0_wp
        rhx                           = 0.e0_wp
        rt34                          = 0.e0_wp
        rt38                          = 0.e0_wp
        scr1                          = 0.e0_wp
        scr2                          = 0.e0_wp
        xib                           = 0.e0_wp
        xih                           = 0.e0_wp
        xix                           = 0.e0_wp
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 0
        coordinate_representation  = 3
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 3
        eps                        = 1.e-6_wp
        one                        = 1.0e0_wp
        rhh                        = 1.4016e0_wp
        ahh                        = 1.0291e0_wp
        dhh                        = 0.17447e0_wp
        bhh                        = 0.018e0_wp
        xih                        = 1.e0_wp
        xix                        = 1.6e0_wp
        three_center_decay         = 0.22e0_wp
        alfw                       = 1.0e0_wp
        rhx                        = 2.673e0_wp
        ahx                        = 0.957e0_wp
        dhx                        = 0.1439e0_wp
        bhx                        = 0.012e0_wp
        etahh                      = 0.393764e0_wp
        eta3s                      = 0.322e0_wp
        eta1p                      = 0.20286e0_wp
        eta3p                      = 0.1771e0_wp
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

        !> @brief Evaluate prepare_local_context for the BrH2_LEPS_C_1982 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: t
            integer  :: i
            !     1          SCR2(4)
            !   Echo potential parameters to UNIT NFLAG(18)
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            easyab             = dhx
            easybc             = dhh
            easyac             = dhx
            xib                = 0.5e0_wp*(xih + xix)
            t                  = xib*rhx
            c3                 = 1.e0_wp/3.e0_wp
            denom              = 2.e0_wp*rhx*(1.e0_wp+t*(1.e0_wp+c3*t))*exp(-t)
            three_center_decay = three_center_decay/rhx
            rt34               = 0.25e0_wp*sqrt(3.e0_wp)
            rt38               = 0.5e0_wp*rt34
            ezero(1)           = dhh
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1) = 'D. C. Clary, Chem. Phys. 71, 117 (1982)'
            ref(2) = 'I. Last and M. Baer in'
            ref(3) = 'Potential Energy Surfaces and Dynamics'
            ref(4) = 'D. G. Truhlar (Ed.) (Plenum, New York, 1981), p. 519'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the BrH2_LEPS_C_1982 potential.
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
            real(wp) :: csa, csa2, csb, csb2, csg2, d1, d2, d3, dcsa21, dcsa22, dcsa23, dcsb21, dcsb22, dcsb23, dcsg21, dcsg22, dcsg23, dp11, dp12, dp21, dp22, dp31, dp32, ds11, ds12, ds21, ds22, ds31, ds32, dshh, dshx1, dshx2, dsn2a1, dsn2a2, dsn2a3, dsn2b1, dsn2b2, dsn2b3, dv1hh, dv1p1, dv1p2, dv1s1, dv1s2, dv3hh, dv3p1, dv3p2, dv3s1, dv3s2, ex1, ex2, p11, p12, p21, p22, p31, p32, r1, r12, r13, r1s, r2, r23, r2s, r3, r3s, rdif, rdif2, s11, s12, s21, s22, s31, s32, shh, shx1, shx2, sn2a, sn2b, sna, sna2, snb, snb2, t, t1, t11, t12, t2, t21, t22, t3, t31, t32, v1hh, v1p1, v1p2, v1s1, v1s2, v3hh, v3p1, v3p2, v3s1, v3s2
            integer  :: ierr, index, k, kk, l, ll
            logical  :: lcol
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   The potential routines define R1=R(Br-H), R2=R(H-Br), and R3=R(H-H);
            !   rearrange the input coordinates to match these definitions
            r1 = r(1)
            r2 = r(3)
            r3 = r(2)
            if( r1 .gt. r2+r3 ) r1 = r2+r3
            if( r2 .gt. r1+r3 ) r2 = r1+r3
            if( r3 .gt. r1+r2 ) r3 = r1+r2
            r1s = r1*r1
            r2s = r2*r2
            r3s = r3*r3
            r12 = r1*r2
            r13 = r1*r3
            r23 = r2*r3
            t3  = (r1s+r2s-r3s)*0.5e0_wp
            t2  = (r1s+r3s-r2s)*0.5e0_wp
            t1  = (r2s+r3s-r1s)*0.5e0_wp
            csa = t2/r13
            if( abs(csa) .gt. one ) csa = sign(one, csa)
            csa2 = csa*csa
            sna2 = 1.e0_wp-csa2
            lcol = sna2 .lt. eps
            sna  = sqrt(sna2)
            t22  = 2.e0_wp*csa
            t11  = 0.e0_wp
            if( .not.lcol ) t11 = 2.e0_wp*(sna2-csa2)/sna
            sn2a   = t22*sna
            t      = t3/(r1*r13)
            dcsa21 = t22*t
            dsn2a1 = t11*t
            t      = -r2/r13
            dcsa22 = t22*t
            dsn2a2 = t11*t
            t      = t1/(r13*r3)
            dcsa23 = t22*t
            dsn2a3 = t11*t
            csb    = t1/r23
            if( abs(csb) .gt. one ) csb = sign(one, csb)
            csb2 = csb*csb
            snb2 = 1.e0_wp-csb2
            snb  = sqrt(snb2)
            t11  = 0.e0_wp
            if( .not.lcol ) t11 = 2.e0_wp*(snb2-csb2)/snb
            t22    = 2.e0_wp*csb
            sn2b   = t22*snb
            t      = -r1/r23
            dcsb21 = t22*t
            dsn2b1 = t11*t
            t      = t3/(r2*r23)
            dcsb22 = t22*t
            dsn2b2 = t11*t
            t      = t2/(r23*r3)
            dcsb23 = t22*t
            dsn2b3 = t11*t
            t      = t3/r12
            csg2   = t*t
            t      = 2.e0_wp*t
            dcsg21 = t*t2/(r1*r12)
            dcsg22 = t*t1/(r12*r2)
            dcsg23 = -t*r3/r12
            !  DIATOMIC CURVES
            !    HH
            rdif  = r3 - rhh
            ex1   = exp(-ahh*rdif)
            rdif2 = rdif*rdif
            ex2   = exp(-bhh*rdif*rdif2)
            t1    = dhh*ex1*ex2
            v1hh  = t1*(ex1-2.e0_wp)
            v3hh  = etahh*t1*(ex1+2.e0_wp)
            t1    = 2.e0_wp*ahh*t1
            t2    = 3.e0_wp*bhh*rdif2
            dv1hh = t1*(1.e0_wp-ex1) - t2*v1hh
            dv3hh = -t1*etahh*(1.e0_wp+ex1) - t2*v3hh
            !    HX
            call vhx(r1, v1s1, v3s1, v1p1, v3p1, dv1s1, dv3s1, dv1p1, dv3p1)
            call vhx(r2, v1s2, v3s2, v1p2, v3p2, dv1s2, dv3s2, dv1p2, dv3p2)
            s11  = v1s1 + 3.e0_wp*v3s1
            s12  = v1s2 + 3.e0_wp*v3s2
            s21  = 3.e0_wp*v1s1 + v3s1
            s22  = 3.e0_wp*v1s2 + v3s2
            s31  = v1s1 - v3s1
            s32  = v1s2 - v3s2
            p11  = v1p1 + 3.e0_wp*v3p1
            p12  = v1p2 + 3.e0_wp*v3p2
            p21  = 3.e0_wp*v1p1 + v3p1
            p22  = 3.e0_wp*v1p2 + v3p2
            p31  = v1p1 - v3p1
            p32  = v1p2 - v3p2
            ds11 = dv1s1 + 3.e0_wp*dv3s1
            ds12 = dv1s2 + 3.e0_wp*dv3s2
            ds21 = 3.e0_wp*dv1s1 + dv3s1
            ds22 = 3.e0_wp*dv1s2 + dv3s2
            ds31 = dv1s1 - dv3s1
            ds32 = dv1s2 - dv3s2
            dp11 = dv1p1 + 3.e0_wp*dv3p1
            dp12 = dv1p2 + 3.e0_wp*dv3p2
            dp21 = 3.e0_wp*dv1p1 + dv3p1
            dp22 = 3.e0_wp*dv1p2 + dv3p2
            dp31 = dv1p1 - dv3p1
            dp32 = dv1p2 - dv3p2
            !    CONSTRUCT 4X4 HAMILTONIAN PUT IN PACKED ARRAY
            h(1)                 = v1hh + 0.25e0_wp*(s11*csa2 + p11*sna2 + s12*csb2 + p12*snb2)
            h(3)                 = v1hh + 0.25e0_wp*(s11*sna2 + p11*csa2 + s12*snb2 + p12*csb2)
            h(6)                 = v3hh + 0.25e0_wp*(s21*csa2 + p21*sna2 + s22*csb2 + p22*snb2)
            h(10)                = v3hh + 0.25e0_wp*(s21*sna2 + p21*csa2 + s22*snb2 + p22*csb2)
            t11                  = s11 - p11
            t12                  = s12 - p12
            h(2)                 = 0.e0_wp
            if( .not.lcol ) h(2) = 0.125e0_wp*(t11*sn2a - t12*sn2b)
            h(4)                 = rt34*(s31*csa2 + p31*sna2 - s32*csb2 - p32*snb2)
            t31                  = s31 - p31
            t32                  = s32 - p32
            h(5)                 = 0.e0_wp
            if( .not.lcol ) h(5) = rt38*(t31*sn2a + t32*sn2b)
            h(7)                 = h(5)
            h(8)                 = rt34*(s31*sna2 + p31*csa2 - s32*snb2 - p32*csb2)
            t21                  = s21 - p21
            t22                  = s22 - p22
            h(9)                 = 0.e0_wp
            if( .not.lcol ) h(9) = 0.125e0_wp*(t21*sn2a - t22*sn2b)
            !   CONSTRUCT DERIVATIVE OF HAMILTONIAN MATRIX
            t       = t11*dcsa21 + t12*dcsb21
            dh1(1)  = 0.25e0_wp*(ds11*csa2 + dp11*sna2 + t)
            dh1(2)  = 0.125e0_wp*((ds11-dp11)*sn2a + t11*dsn2a1 - t12*dsn2b1)
            dh1(3)  = 0.25e0_wp*(ds11*sna2 + dp11*csa2 - t)
            t       = t11*dcsa22 + t12*dcsb22
            dh2(1)  = 0.25e0_wp*(ds12*csb2 + dp12*snb2 + t)
            dh2(2)  = 0.125e0_wp*(-(ds12-dp12)*sn2b + t11*dsn2a2 - t12*dsn2b2)
            dh2(3)  = 0.25e0_wp*(ds12*snb2 + dp12*csb2 - t)
            t       = 0.25e0_wp*(t11*dcsa23 + t12*dcsb23)
            dh3(1)  = dv1hh + t
            dh3(2)  = 0.125e0_wp*(t11*dsn2a3 - t12*dsn2b3)
            dh3(3)  = dv1hh - t
            t       = t21*dcsa21 + t22*dcsb21
            dh1(6)  = 0.25e0_wp*(ds21*csa2 + dp21*sna2 + t)
            dh1(9)  = 0.125e0_wp*((ds21-dp21)*sn2a + t21*dsn2a1 - t22*dsn2b1)
            dh1(10) = 0.25e0_wp*(ds21*sna2 + dp21*csa2 - t)
            t       = t21*dcsa22 + t22*dcsb22
            dh2(6)  = 0.25e0_wp*(ds22*csb2 + dp22*snb2 + t)
            dh2(9)  = 0.125e0_wp*((dp22-ds22)*sn2b + t21*dsn2a2 - t22*dsn2b2)
            dh2(10) = 0.25e0_wp*(ds22*snb2 + dp22*csb2 - t)
            t       = 0.25e0_wp*(t21*dcsa23 + t22*dcsb23)
            dh3(6)  = dv3hh + t
            dh3(9)  = 0.125e0_wp*(t21*dsn2a3 - t22*dsn2b3)
            dh3(10) = dv3hh - t
            t       = t31*dcsa21 - t32*dcsb21
            dh1(4)  = rt34*(ds31*csa2 + dp31*sna2 + t)
            dh1(5)  = rt38*((ds31-dp31)*sn2a + t31*dsn2a1 + t32*dsn2b1)
            dh1(7)  = dh1(5)
            dh1(8)  = rt34*(ds31*sna2 + dp31*csa2 - t)
            t       = t31*dcsa22 - t32*dcsb22
            dh2(4)  = rt34*(-ds32*csb2 - dp32*snb2 + t)
            dh2(5)  = rt38*((ds32-dp32)*sn2b + t31*dsn2a2 + t32*dsn2b2)
            dh2(7)  = dh2(5)
            dh2(8)  = rt34*(-ds32*snb2 - dp32*csb2 - t)
            t       = t31*dcsa23 - t32*dcsb23
            dh3(4)  = rt34*t
            dh3(5)  = rt38*(t31*dsn2a3 + t32*dsn2b3)
            dh3(7)  = dh3(5)
            dh3(8)  = -dh3(4)
            !    DIAGONALIZE H PUT EIGENVECTORS INTO eigenvectors, EIGENVALUES INTO E
            !    USE EISPACK ROUTINE RSP
            call rsp(4, 4, 10, h, e, 1, eigenvectors, scr1, scr2, ierr)
            if( ierr .ne. 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'POT 3')
                return
            endif
            !    EVALUATE DERIVATIVES OF LOWEST ROOT
            d1 = 0.e0_wp
            d2 = 0.e0_wp
            d3 = 0.e0_wp
            do l = 1, 4
                t1 = eigenvectors(l,1)
                ll = (l*(l-1))/2
                do k = 1, 4
                    t2 = t1*eigenvectors(k,1)
                    kk = (k*(k-1))/2
                    if( l .ge. k ) index = ll + k
                    if( k .gt. l ) index = kk + l
                    d1 = d1 + t2*dh1(index)
                    d2 = d2 + t2*dh2(index)
                    d3 = d3 + t2*dh3(index)
                enddo
            enddo
            !   SUPPLEMENTARY TERM
            t                   = xih*r3
            ex1                 = exp(-t)
            shh                 = (1.e0_wp+t*(1.e0_wp+c3*t))*ex1
            dshh                = -xih*t*(1.e0_wp+t)*c3*ex1
            t                   = xib*r1
            ex1                 = exp(-t)
            shx1                = r1*(1.e0_wp+t*(1.e0_wp+c3*t))*ex1/denom
            dshx1               = (1.e0_wp+t*(1.e0_wp-c3*t*t))*ex1/denom
            t                   = xib*r2
            ex1                 = exp(-t)
            shx2                = r2*(1.e0_wp+t*(1.e0_wp+c3*t))*ex1/denom
            dshx2               = (1.e0_wp+t*(1.e0_wp-c3*t*t))*ex1/denom
            rdif                = r1 - r2
            t1                  = three_center_decay*exp(-alfw*rdif*rdif)
            t2                  = shh*(shx1+shx2) + shx1*shx2
            ground_state_energy = e(1) + t2*t1*csg2 + ezero(1)
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                t3                                = 2.e0_wp*alfw*rdif*csg2
                ground_state_internal_gradient(1) = d1 + (dshx1*(shh+shx2)*csg2 + t2*(dcsg21-t3))*t1
                ground_state_internal_gradient(3) = d2 + (dshx2*(shh+shx1)*csg2 + t2*(dcsg22+t3))*t1
                ground_state_internal_gradient(2) = d3 + (dshh*(shx1+shx2)*csg2 + t2*dcsg23)*t1
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate VHX for the BrH2_LEPS_C_1982 potential.
        subroutine vhx( rr, v1s, v3s, v1p, v3p, dv1s, dv3s, dv1p, dv3p )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dv1, dv1p, dv1s, dv3p, dv3s, ex1, ex2, rdif, rdif2, rr, t1, t2, v1, v1p, v1s, v3p, v3s
            rdif  = rr - rhx
            rdif2 = rdif*rdif
            ex1   = exp(-ahx*rdif)
            ex2   = exp(-bhx*rdif2*rdif)
            t1    = dhx*ex1*ex2
            v1s   = t1*(ex1-2.e0_wp)
            v1    = t1*(ex1+2.e0_wp)
            v3s   = eta3s*v1
            v1p   = eta1p*v1
            v3p   = eta3p*v1
            t1    = 2.e0_wp*ahx*t1
            t2    = 3.e0_wp*bhx*rdif2
            dv1s  = t1*(1.e0_wp-ex1) - t2*v1s
            dv1   = -t1*(1.e0_wp+ex1) - t2*v1
            dv3s  = eta3s*dv1
            dv1p  = eta1p*dv1
            dv3p  = eta3p*dv1
        end subroutine vhx

        !> @brief Evaluate RSP for the BrH2_LEPS_C_1982 potential.
        !> @brief Diagonalize a real symmetric packed matrix with the shared core eigensolver.
        subroutine rsp( nm, n, nv, a, w, matz, z, fv1, fv2, ierr )
            use chempotpy_kinds, only: wp
            implicit none
            integer, intent(in)   :: nm, n, nv, matz
            real(wp), intent(in)  :: a(nv)
            real(wp), intent(out) :: w(n), z(nm,n), fv1(n), fv2(n)
            integer, intent(out)  :: ierr
            integer               :: i, j, packed_index
            real(wp), allocatable :: eigenvectors(:,:), matrix(:,:)
            ierr = 0
            z    = 0.0e0_wp
            fv1  = 0.0e0_wp
            fv2  = 0.0e0_wp
            if( n > nm )then
                ierr = 10*n
                return
            endif
            if( nv < n*(n+1)/2 )then
                ierr = 20*n
                return
            endif
            allocate(eigenvectors(n,n), matrix(n,n))
            matrix       = 0.0e0_wp
            packed_index = 0
            do i = 1, n
                do j = 1, i
                    packed_index = packed_index+1
                    matrix(i,j)  = a(packed_index)
                    matrix(j,i)  = a(packed_index)
                enddo
            enddo
            if( matz == 0 )then
                call symmetric_eigendecomposition(matrix, w, info=ierr)
            else
                call symmetric_eigendecomposition(matrix, w, eigenvectors, ierr)
                if( ierr == 0 ) z(1:n,1:n) = eigenvectors
            endif
            deallocate(eigenvectors, matrix)
        end subroutine rsp

        !> @brief Evaluate POTINFO for the BrH2_LEPS_C_1982 potential.
        subroutine potinfo
            implicit none
            integer            :: i, j, k, kmax
            integer, parameter :: isurf = 5
            do i =1, 5
            enddo
            kmax = 0
            do i = 1, isurf+1
                do j = 1, isurf+1
                    if( surface_requests(i,j) .ne. 0.and.kmax .lt. max(i, j) ) kmax = max(i, j)
                enddo
            enddo
            if( kmax-1 .gt. maximum_excited_state )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            first_requested_excited_state = 0
            last_requested_excited_state  = 0
            do i = 2, isurf+1
                if( surface_requests(i,i) .ne. 0 )then
                    last_requested_excited_state = i-1
                    if( first_requested_excited_state .eq. 0 ) first_requested_excited_state = i-1
                endif
            enddo
            first_requested_coupling = 0
            last_requested_coupling  = 0
            k                        = 0
            do i = 1, isurf
                do j = i+1, isurf+1
                    k = k+1
                    if( surface_requests(i,j)+surface_requests(j,i) .ne. 0 )then
                        last_requested_coupling = k
                        if( first_requested_coupling .eq. 0 ) first_requested_coupling = k
                    endif
                enddo
            enddo
            if( requested_derivative_order .gt. maximum_derivative_order )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine potinfo

        !> @brief Evaluate CARTTOR for the BrH2_LEPS_C_1982 potential.
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

        !> @brief Evaluate EUNITZERO for the BrH2_LEPS_C_1982 potential.
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

        !> @brief Evaluate RTOCART for the BrH2_LEPS_C_1982 potential.
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

        !> @brief Evaluate DEDCOU for the BrH2_LEPS_C_1982 potential.
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
end module chempotpy_brh2_leps_c_1982_pes

!> @brief Expose BrH2_LEPS_C_1982 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_brh2_leps_c_1982_pes, only: evaluate
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
