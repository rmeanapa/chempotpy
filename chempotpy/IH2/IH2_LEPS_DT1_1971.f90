!   System:          IH2
!   Functional form:
!   Common name:     RMC
!   References:
! Number of bodies: 3
! Interface: potlib2001
! Number of electronic states: 1
! Number of derivatives: 1
!   Note: This is the one-dimensional (collinear) version of the IH2
!         RMC potential energy surface; there are no angle dependent
!         terms in this version of the potential.
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
!        surface_requests - which electronic states availalble
!                 (1,1) = 1 as only gs state available
!        requested_derivative_order  = 0 => no derivatives should be calculated
!        requested_derivative_order  = 1 => calculate first derivatives
!        NFLAG - these integer values can be used to flag options
!                within the potential;
!   Potential parameters' default settings
!                  Variable            Default value
!                  requested_derivative_order                1
!                  NFLAG(18)           6

!> @brief Implement IH2_LEPS_DT1_1971, a potential energy surface for IH2.
module chempotpy_ih2_leps_dt1_1971_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the IH2_LEPS_DT1_1971 potential energy surface for IH2.
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
        real(wp)              :: a(2)
        real(wp)              :: alfd
        real(wp)              :: alpha(2)
        real(wp)              :: c(2,2)
        real(wp)              :: input_cartesian(natom,3)
        real(wp)              :: cartesian_coordinates(natom,3)
        real(wp)              :: clbh
        real(wp)              :: dadp(2)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        real(wp)              :: dmdn
        real(wp)              :: dpdr(2)
        real(wp)              :: easyab
        real(wp)              :: easyac
        real(wp)              :: easybc
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: ground_state_energy
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: ezero(isurf+1)
        integer               :: coordinate_representation
        integer               :: irxn
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        integer               :: l(2)
        real(wp)              :: dissociation_energies(2)
        integer               :: atom_count
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: pid2
        real(wp)              :: pid4
        real(wp)              :: pre(2)
        real(wp)              :: r(n3atom)
        real(wp)              :: re(2)
        character(len=75)     :: ref(5)
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
        a                              = 0.e0_wp; alfd = 0.e0_wp; alpha = 0.e0_wp; c = 0.e0_wp; input_cartesian = 0.e0_wp; cartesian_coordinates = 0.e0_wp; clbh = 0.e0_wp
        dadp                           = 0.e0_wp; excited_state_internal_gradients = 0.e0_wp; ground_state_internal_gradient = 0.e0_wp; diabatic_internal_gradients = 0.e0_wp
        output_excited_state_gradients = 0.e0_wp; excited_state_cartesian_gradients = 0.e0_wp; output_ground_state_gradient = 0.e0_wp; ground_state_cartesian_gradient = 0.e0_wp; output_diabatic_gradients = 0.e0_wp; diabatic_cartesian_gradients = 0.e0_wp
        dmdn                           = 0.e0_wp; dpdr = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp; excited_state_energies = 0.e0_wp; ground_state_energy = 0.e0_wp
        diabatic_couplings             = 0.e0_wp; ezero = 0.e0_wp; coordinate_representation = 0; irxn = 0; last_requested_excited_state = 0; last_requested_coupling = 0
        first_requested_excited_state  = 0; first_requested_coupling = 0; l = 0; dissociation_energies = 0.e0_wp; atom_count = 0; maximum_derivative_order = 0; maximum_excited_state = 0; surface_requests = 0; requested_derivative_order = 0
        output_excited_state_energies  = 0.e0_wp; output_ground_state_energy = 0.e0_wp; output_diabatic_couplings = 0.e0_wp; pid2 = 0.e0_wp; pid4 = 0.e0_wp; pre = 0.e0_wp
        r                              = 0.e0_wp; re = 0.e0_wp
        ref                            = ' '
        rr                             = 0.e0_wp
        rr1ms1                         = 0.e0_wp
        rr2ms2                         = 0.e0_wp
        rrmr                           = 0.e0_wp
        rrmre                          = 0.e0_wp
        rrs1d2                         = 0.e0_wp
        s                              = 0.e0_wp
        sr1r2                          = 0.e0_wp
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 0
        coordinate_representation  = 3
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 3
        alpha                      = [1.789e0_wp, 1.974e0_wp]
        dissociation_energies      = [73.682e0_wp, 109.499e0_wp]
        re                         = [1.604e0_wp, 0.742e0_wp]
        s                          = [1.6783e0_wp, 1.2801e0_wp]
        rr                         = [2.750e0_wp, 1.750e0_wp]
        alfd                       = 1.500e0_wp
        clbh                       = 35.879e0_wp
        l                          = [12, 6]
        irxn                       = 1
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

        !> @brief Evaluate prepare_local_context for the IH2_LEPS_DT1_1971 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, m, n
            !   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr
            real(wp), parameter :: ckcal = hartree_to_kilocalorie_per_mole
            real(wp), parameter :: cangs = 0.52917706e0_wp
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            pid4 = atan(1.e0_wp)
            pid2 = 2.e0_wp*pid4
            m    = irxn
            n    = 1
            if( irxn .eq. 1 ) n = 2
            c(m,1) = clbh/dissociation_energies(m)
            if( irxn .eq. 1 ) c(m,1) = c(m,1) + (dissociation_energies(m)-dissociation_energies(n))/ &
                & dissociation_energies(m)
            c(m,2) = 0.0e0_wp
            dmdn   = dissociation_energies(m)/dissociation_energies(n)
            c(n,2)                   = 2.0e0_wp*dmdn*c(m,1)*real(l(n)-l(m), wp) + 2.0e0_wp*real(l(n), wp)*(1.0e0_wp- &
                & dmdn)
            c(n,1) = 1.0e0_wp-dmdn*(1.0e0_wp-c(m,1))+pid4*pid4*c(n,2)
            !   Echo the potential parameters to the file linked to UNIT NFLAG(18).
            !   CONVERT TO ATOMIC UNITS
            do i = 1, 2
                rr(i)                    = rr(i)/cangs
                s(i)                     = s(i)/cangs
                re(i)                    = re(i)/cangs
                alpha(i)                 = alpha(i)*cangs
                dissociation_energies(i) = dissociation_energies(i)/ckcal
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
            easyab   = dissociation_energies(1)
            easybc   = dissociation_energies(2)
            easyac   = dissociation_energies(1)
            ezero(1) = dissociation_energies(2)
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1) = 'No Reference'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the IH2_LEPS_DT1_1971 potential.
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
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: alf, cs2pp, daldr, dgadr, ds4dr, dsldr, dxidr, fac, gam, phip, quan, quan2, rr12, sn2pp, sn2pp3, sn2pp4, sn2ppl, tn2pp, xi
            integer  :: i, j, lm3
            !   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr
            logical :: lrr1, lrr2
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   Check the values of surface_requests and requested_derivative_order for validity.
            if( surface_requests(1,1) .eq. 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            if( requested_derivative_order .gt. 1 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'POT 2')
                return
            endif
            rrmr(1) = rr(1)-r(1)
            rrmr(2) = rr(2)-r(2)
            lrr1    = rrmr(1).le.0.0e0_wp
            lrr2    = rrmr(2).le.0.0e0_wp
            if( requested_derivative_order == 1 ) ground_state_internal_gradient = 0.0e0_wp
            ground_state_energy = ezero(1)
            if( lrr1 .and. lrr2 )then
                call eunitzero
                if( requested_derivative_order /= 0 )then
                    call rtocart
                    call dedcou
                endif
                return
            endif
            if( lrr1 .or. lrr2 )then
                j = 1
                if( lrr1 ) j                = 2
                quan                = 1.0e0_wp-exp(-alpha(j)*(r(j)-re(j)))
                ground_state_energy = dissociation_energies(j)*(quan**2-1.0e0_wp)+ezero(1)
                if( requested_derivative_order == 1 ) ground_state_internal_gradient(j)   = 2.0e0_wp*dissociation_energies(j)*alpha(j)*quan*(1.0e0_wp-quan)
            else
                phip   = atan(rrs1d2*rrmr(2)/rrmr(1))
                sn2pp  = sin(2.0e0_wp*phip)
                cs2pp  = cos(2.0e0_wp*phip)
                tn2pp  = sn2pp/cs2pp
                sn2pp3 = sn2pp**3
                sn2pp4 = sn2pp3*sn2pp
                i      = 1
                if( phip > pid4 ) i = 2
                a(i) = c(i,1)+c(i,2)*(phip-pid2)*phip
                if( requested_derivative_order == 1 ) dadp(i)  = c(i,2)*(2.0e0_wp*phip-pid2)
                rr12                = sqrt(rrmr(1)**2+rrmr(2)**2)
                lm3                 = l(i)-3
                sn2ppl              = sn2pp3*sn2pp**lm3
                alf                 = alpha(i)+(alfd-alpha(i))*sn2pp4
                gam                 = dissociation_energies(i)*(1.0e0_wp-a(i)*sn2ppl)
                xi                  = pre(i)*sn2pp4+rrmre(i)-rr12
                quan                = 1.0e0_wp-exp(-alf*xi)
                quan2               = quan*quan
                ground_state_energy = gam*(quan2-1.0e0_wp)+dissociation_energies(2)
                if( requested_derivative_order == 1 )then
                    fac     = rrs1d2/(rrmr(1)**2+(rrs1d2*rrmr(2))**2)
                    dpdr(1) = fac*rrmr(2)
                    dpdr(2) = -fac*rrmr(1)
                    do j = 1, 2
                        ds4dr = 8.0e0_wp*sn2pp4/tn2pp*dpdr(j)
                        dsldr = 2.0e0_wp*real(l(i), wp)*sn2ppl/tn2pp*dpdr(j)
                        daldr = (alfd-alpha(i))*ds4dr
                        dgadr = -dissociation_energies(i)*a(i)*dsldr-dissociation_energies(i)*dadp(i)*sn2ppl*dpdr(j)
                        dxidr = pre(i)*ds4dr+rrmr(j)/rr12
                        ground_state_internal_gradient(j) = dgadr*(quan2-1.0e0_wp)+2.0e0_wp*gam*quan*(1.0e0_wp-quan)* &
                            & (alf*dxidr+xi*daldr)
                    enddo
                endif
            endif
            call eunitzero
            if( requested_derivative_order /= 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate CARTTOR for the IH2_LEPS_DT1_1971 potential.
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

        !> @brief Evaluate EUNITZERO for the IH2_LEPS_DT1_1971 potential.
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

        !> @brief Evaluate RTOCART for the IH2_LEPS_DT1_1971 potential.
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

        !> @brief Evaluate DEDCOU for the IH2_LEPS_DT1_1971 potential.
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
end module chempotpy_ih2_leps_dt1_1971_pes

!> @brief Expose IH2_LEPS_DT1_1971 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_ih2_leps_dt1_1971_pes, only: evaluate
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
