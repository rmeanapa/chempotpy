!> @brief Implement H3_GEN_PK2_1964, a potential energy surface for H3.
module chempotpy_h3_gen_pk2_1964_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_PK2_1964 potential energy surface for H3.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units,  only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: nstates               = 1
        integer, parameter    :: number_of_atoms       = 3
        integer, parameter    :: natom                 = 25
        integer, parameter    :: n3atom                = 75
        integer, parameter    :: isurf                 = 5
        integer, parameter    :: jsurf                 = isurf*(isurf+1)/2
        real(wp), parameter   :: legacy_ev_per_hartree = 27.21161e0_wp
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(number_of_atoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,number_of_atoms,3)
        real(wp), intent(out) :: d(nstates,nstates,number_of_atoms,3)
        real(wp)              :: alph, beta, d1, d3, del, easyab
        real(wp)              :: easyac, easybc, ground_state_energy, eps, fkap, gam, output_ground_state_energy, re, v
        real(wp)              :: input_cartesian(natom,3), cartesian_coordinates(natom,3), excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: ground_state_internal_gradient(n3atom), diabatic_internal_gradients(n3atom,jsurf), output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf), output_ground_state_gradient(natom,3), ground_state_cartesian_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf), diabatic_cartesian_gradients(natom,3,jsurf), excited_state_energies(isurf)
        real(wp)              :: diabatic_couplings(jsurf), ezero(isurf+1), output_excited_state_energies(isurf), output_diabatic_couplings(jsurf)
        real(wp)              :: r(n3atom)
        integer               :: iatom, coordinate_representation, idir, istate
        integer               :: last_requested_excited_state, last_requested_coupling, first_requested_excited_state, first_requested_coupling, maximum_derivative_order, maximum_excited_state, natoms, requested_derivative_order
        integer               :: surface_requests(isurf+1,isurf+1)
        character(len=75)     :: ref(5)
        v                              = 0.e0_wp; p = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        alph                           = 0.e0_wp; beta = 0.e0_wp
        d1                             = 0.e0_wp; d3 = 0.e0_wp; del = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp
        ground_state_energy            = 0.e0_wp; eps = 0.e0_wp; fkap = 0.e0_wp; gam = 0.e0_wp; output_ground_state_energy = 0.e0_wp; re = 0.e0_wp
        input_cartesian                = 0.e0_wp; cartesian_coordinates = 0.e0_wp; excited_state_internal_gradients = 0.e0_wp; ground_state_internal_gradient = 0.e0_wp; diabatic_internal_gradients = 0.e0_wp
        output_excited_state_gradients = 0.e0_wp; excited_state_cartesian_gradients = 0.e0_wp; output_ground_state_gradient = 0.e0_wp; ground_state_cartesian_gradient = 0.e0_wp
        output_diabatic_gradients      = 0.e0_wp; diabatic_cartesian_gradients = 0.e0_wp; excited_state_energies = 0.e0_wp; diabatic_couplings = 0.e0_wp; ezero = 0.e0_wp
        output_excited_state_energies  = 0.e0_wp; output_diabatic_couplings = 0.e0_wp; r = 0.e0_wp
        last_requested_excited_state   = 0; last_requested_coupling = 0; first_requested_excited_state = 0; first_requested_coupling = 0
        maximum_derivative_order    = 1; maximum_excited_state = 0; natoms = number_of_atoms; requested_derivative_order = 0; surface_requests = 0; surface_requests(1, &
            & &
            & 1) = 1
        coordinate_representation = 3; ref = ' '
        d1                        = 4.7466e0_wp; d3 = 1.9668e0_wp; re = 1.40083e0_wp; alph = 1.04435e0_wp; beta = 1.000122e0_wp
        eps                       = -17.5e0_wp; del = 28.2e0_wp; fkap = 0.60e0_wp; gam = 0.65e0_wp
        do iatom=1, number_of_atoms
            do idir=1, 3
                input_cartesian(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        requested_derivative_order = igrad
        d1                         = d1/legacy_ev_per_hartree; d3 = d3/legacy_ev_per_hartree
        eps                        = eps/legacy_ev_per_hartree; del = del/legacy_ev_per_hartree
        easyab                     = d1; easybc = easyab; easyac = easyab; ezero(1) = d1
        ref(1)                     = 'R. N. Porter and M. Karplus'
        ref(2)                     = 'J. Chem. Phys. 40, 1105(1964)'
        call evaluate_kernel
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = output_ground_state_energy*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = output_ground_state_energy*hartree_to_electronvolt
            enddo
            do iatom=1, number_of_atoms
                do idir=1, 3
                    g(1,iatom,idir) = output_ground_state_gradient(iatom,idir)*hartree_per_bohr_to_ev_per_angstrom
                enddo
            enddo
        elseif( igrad == 2 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives above first order are unavailable')
        endif

        contains

        !> @brief Evaluate the local PK2 potential kernel.
        subroutine evaluate_kernel
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the H2 diatomic, with the
            !   H2 diatomic at its equilibrium configuration.  For this potential en
            !   surface the potential energy in the three asymptotic valleys are equ
            !   The potential energy is defined to be equal to zero when one H is
            !   "infinitely" far from the H2 diatomic, and the H2 diatomic bond
            !   length is equal to the H2 equilibrium bond length.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a, a1, a2, a3, b, b1, b2, b3, bac, c, c1, c2, c3, cs123, de1(3), de3(3), delj(3), &
                & dife(3), dj(3), dj112, dj113, dj12, dj123, dj13, dj21, dj212, dj213, dj223, dj23, dj31, dj312, dj313, dj32, dj323, &
                & dq(3), ds(3), e1(3), e3(3), engygs0, ex(3), exa, exb, expg, expr, f, fj(3), fj12, fj13, fj23, gamr, q(3), qt, &
                & rbac, rel, rho, rr(3), s(3), s12, s123, s13, s23, ss(3), sume(3), t, zeta
            integer :: i
            !   Conversion constants
            !   CEV (eV) converts from Hartree to kcal/mol
            cartesian_coordinates(1:natoms,:) = input_cartesian(1:natoms,:)
            call carttor
            !   Check the values of surface_requests and requested_derivative_order for validity.
            do i = 1, 3
                gamr    = gam*r(i)
                expg    = exp(-gamr)
                zeta    = 1.0e0_wp + fkap*expg
                rho     = r(i)*zeta
                rr(i)   = 1.0e0_wp/r(i)
                ex(i)   = del*exp(-r(i)-r(i))
                expr    = exp(-rho)
                s(i)    = (1.0e0_wp+rho+(rho*rho)/3.0e0_wp)*expr
                ss(i)   = s(i)*s(i)
                rel     = re-r(i)
                exa     = exp(alph*rel)
                exb     = exp(beta*rel)
                e1(i)   = d1*(1.0e0_wp-exa)**2-d1
                e3(i)   = d3*(1.0e0_wp+exb)**2-d3
                sume(i) = e1(i)+e3(i)
                dife(i) = e1(i)-e3(i)
                q(i)    = 0.5e0_wp*(sume(i)+ss(i)*dife(i))
                fj(i)   = 0.5e0_wp*(dife(i)+ss(i)*sume(i))
                de1(i)  = 2.0e0_wp*d1*alph*exa*(1.0e0_wp-exa)
                de3(i)  = 2.0e0_wp*d3*beta*exb*(-1.0e0_wp-exb)
                ds(i)   = -rho*(1.0e0_wp+rho)*expr*(zeta+gamr*(1.0e0_wp-zeta))/3.0e0_wp
                dj(i)   = s(i)*sume(i)*ds(i)+0.5e0_wp*((1.0e0_wp+ss(i))*de1(i)-(1.0e0_wp-ss(i ))* de3(i))
                delj(i) = (1.0e0_wp+rr(i))*ex(i)
                dq(i)   = s(i)*dife(i)*ds(i)+0.5e0_wp*((1.0e0_wp+ss(i))*de1(i)+(1.0e0_wp-ss(i ))* de3(i))
            enddo
            dj(1)   = dj(1) + 2.0e0_wp*s(1)*ds(1)*(delj(2)+delj(3))
            dj(2)   = dj(2) + 2.0e0_wp*s(2)*ds(2)*(delj(1)+delj(3))
            dj(3)   = dj(3) + 2.0e0_wp*s(3)*ds(3)*(delj(1)+delj(2))
            fj(1)   = fj(1) + ss(1)*(delj(2)+delj(3))
            fj(2)   = fj(2) + ss(2)*(delj(1)+delj(3))
            fj(3)   = fj(3) + ss(3)*(delj(1)+delj(2))
            s12     = ss(1)-ss(2)
            s23     = ss(2)-ss(3)
            s13     = s12+s23
            fj12    = fj(1)-fj(2)
            fj23    = fj(2)-fj(3)
            fj13    = fj12+fj23
            dj12    = -ss(2)*ex(1)*(2.0e0_wp+2.0e0_wp*rr(1)+rr(1)*rr(1))
            dj13    = -ss(3)*ex(1)*(2.0e0_wp+2.0e0_wp*rr(1)+rr(1)*rr(1))
            dj21    = -ss(1)*ex(2)*(2.0e0_wp+2.0e0_wp*rr(2)+rr(2)*rr(2))
            dj23    = -ss(3)*ex(2)*(2.0e0_wp+2.0e0_wp*rr(2)+rr(2)*rr(2))
            dj31    = -ss(1)*ex(3)*(2.0e0_wp+2.0e0_wp*rr(3)+rr(3)*rr(3))
            dj32    = -ss(2)*ex(3)*(2.0e0_wp+2.0e0_wp*rr(3)+rr(3)*rr(3))
            dj112   = dj(1)-dj12
            dj113   = dj(1)-dj13
            dj212   = dj21-dj(2)
            dj223   = dj(2)-dj23
            dj313   = dj31-dj(3)
            dj323   = dj32-dj(3)
            dj123   = dj12-dj13
            dj213   = dj21-dj23
            dj312   = dj31-dj32
            s123    = s(1)*s(2)*s(3)
            qt      = q(1)+q(2)+q(3)-eps*s123
            cs123   = 1.0e0_wp-s123
            a       = cs123*cs123-0.5e0_wp*(s12*s12+s23*s23+s13*s13)
            b       = -qt*cs123+0.5e0_wp*(fj12*s12+fj23*s23+fj13*s13)
            c       = qt*qt-0.5e0_wp*(fj12*fj12+fj23*fj23+fj13*fj13)
            a1      = -(ds(1)+ds(1))*(cs123*s(2)*s(3)+s(1)*(s12+s13))
            a2      = -(ds(2)+ds(2))*(cs123*s(1)*s(3)+s(2)*(s23-s12))
            a3      = -(ds(3)+ds(3))*(cs123*s(1)*s(2)-s(3)*(s13+s23))
            f       = cs123*eps+qt
            b1      = -cs123*dq(1)+ds(1)*(s(2)*s(3)*f+s(1)*(fj12+fj13)) +0.5e0_wp*(s12*dj112+s23*dj123+s13*dj113)
            b2      = -cs123*dq(2)+ds(2)*(s(1)*s(3)*f+s(2)*(fj23-fj12)) +0.5e0_wp*(s12*dj212+s23*dj223+s13*dj213)
            b3      = -cs123*dq(3)+ds(3)*(s(1)*s(2)*f-s(3)*(fj13+fj23)) +0.5e0_wp*(s12*dj312+s23*dj323+s13*dj313)
            f       = qt+qt
            c1      = f*(dq(1)-eps*s(2)*s(3)*ds(1))-(fj12*dj112+fj23*dj123 +fj13*dj113)
            c2      = f*(dq(2)-eps*s(1)*s(3)*ds(2))-(fj12*dj212+fj23*dj223 +fj13*dj213)
            c3      = f*(dq(3)-eps*s(2)*s(1)*ds(3))-(fj12*dj312+fj23*dj323 +fj13*dj313)
            bac     = b*b-a*c
            rbac    = sqrt(bac)
            engygs0 = -b - rbac
            if( engygs0 .ne. 0.e0_wp ) engygs0 = engygs0/a
            ground_state_energy = engygs0 + ezero(1)
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                t = 2.0e0_wp*rbac
                ground_state_internal_gradient(1)                              = engygs0*engygs0*a1+(engygs0+engygs0)* &
                    & b1+c1
                if( ground_state_internal_gradient(1) .ne. 0.e0_wp ) ground_state_internal_gradient(1) = ground_state_internal_gradient(1)/ &
                    & &
                    & t
                ground_state_internal_gradient(2)                              = engygs0*engygs0*a2+(engygs0+engygs0)* &
                    & b2+c2
                if( ground_state_internal_gradient(2) .ne. 0.e0_wp ) ground_state_internal_gradient(2) = ground_state_internal_gradient(2)/ &
                    & &
                    & t
                ground_state_internal_gradient(3)                              = engygs0*engygs0*a3+(engygs0+engygs0)* &
                    & b3+c3
                if( ground_state_internal_gradient(3) .ne. 0.e0_wp ) ground_state_internal_gradient(3) = ground_state_internal_gradient(3)/ &
                    & &
                    & t
                call rtocart
                call dedcou
            endif
        end subroutine evaluate_kernel

        !> @brief Evaluate CARTTOR for the H3_GEN_PK2_1964 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(cartesian_coordinates, natoms, coordinate_representation, r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the H3_GEN_PK2_1964 potential.
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

        !> @brief Evaluate RTOCART for the H3_GEN_PK2_1964 potential.
        subroutine rtocart
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            real(wp)           :: epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, termy, &
                & termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes(n3atom, isurf), ygs(n3atom), &
                & yij(n3atom, jsurf), za, zb, zc, zd
            integer :: i, ind, j, j1, j2, k, m1
            if( coordinate_representation .eq. 1 )then
                do i = 1, natoms
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
                do i = 1, natoms
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
                    do j = 1, natoms
                        if( j .lt. i )then
                            m1 = natoms*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
                            m1 = natoms*(i-1) - (i*(i-1))/2 + j-i
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
                do i = 1, natoms
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

        !> @brief Evaluate DEDCOU for the H3_GEN_PK2_1964 potential.
        subroutine dedcou
            implicit none
            integer :: i, j, k
            do i = 1, natoms
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
end module chempotpy_h3_gen_pk2_1964_pes

!> @brief Expose H3_GEN_PK2_1964 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_pk2_1964_pes, only: evaluate
    use chempotpy_kinds,               only: wp
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
