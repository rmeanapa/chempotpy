!   System:                        OH3
!   Functional form:
!   Common name:
!   Number of derivatives:         1
!   Number of bodies:              4
!   Number of electronic surfaces: 1
!   Interface:                     4-XS
!   Reference: G. C. Schatz and H. Elgersma
!              Chem. Phys. Lett. 73, 21 (1980).

!> @brief Implement OH3_LEPS_SE_1986, a potential energy surface for OH3.
module chempotpy_oh3_leps_se_1986_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the OH3_LEPS_SE_1986 potential energy surface for OH3.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units,  only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 4
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter    :: n3atom  = 75
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: concom_acon(2)
        real(wp)              :: concom_alp(4)
        real(wp)              :: concom_beta(3)
        real(wp)              :: concom_clam(4)
        real(wp)              :: concom_con(7)
        real(wp)              :: concom_de(3)
        real(wp)              :: concom_dum(22)
        real(wp)              :: concom_dum1(10)
        real(wp)              :: concom_dum2(10)
        real(wp)              :: concom_gam(3)
        real(wp)              :: concom_re(3)
        real(wp)              :: concom_rehh
        real(wp)              :: concom_reoh
        real(wp)              :: concom_sato
        real(wp)              :: cartesian_coordinates(natom,3)
        integer               :: coordinate_representation
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: atom_count
        character(len=75)     :: infocm_ref(5)
        real(wp)              :: pot2cm_dedr(4)
        real(wp)              :: pot2cm_energy
        real(wp)              :: pot2cm_r(4)
        real(wp)              :: potcm_dvdr(6)
        real(wp)              :: potcm_r(6)
        real(wp)              :: potcm_vtot
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: ground_state_energy
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: input_cartesian(natom,3)
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        concom_acon                       = 0.e0_wp; concom_alp = 0.e0_wp; concom_beta = 0.e0_wp; concom_clam = 0.e0_wp; concom_con = 0.e0_wp
        concom_de                         = 0.e0_wp; concom_dum = 0.e0_wp; concom_dum1 = 0.e0_wp; concom_dum2 = 0.e0_wp; concom_gam = 0.e0_wp; concom_re = 0.e0_wp
        concom_rehh                       = 0.e0_wp; concom_reoh = 0.e0_wp; concom_sato = 0.e0_wp; cartesian_coordinates = 0.e0_wp; coordinate_representation = 0
        maximum_derivative_order          = 0; maximum_excited_state = 0; atom_count = 0
        infocm_ref                        = ' '
        pot2cm_dedr                       = 0.e0_wp
        pot2cm_energy                     = 0.e0_wp
        pot2cm_r                          = 0.e0_wp
        potcm_dvdr                        = 0.e0_wp
        potcm_r                           = 0.e0_wp
        potcm_vtot                        = 0.e0_wp
        ground_state_internal_gradient    = 0.e0_wp
        ground_state_energy               = 0.e0_wp
        pt1cm_r                           = 0.e0_wp
        pt3cm_ezero                       = 0.e0_wp
        excited_state_internal_gradients  = 0.e0_wp
        excited_state_energies            = 0.e0_wp
        diabatic_internal_gradients       = 0.e0_wp
        diabatic_couplings                = 0.e0_wp
        input_cartesian                   = 0.e0_wp
        surface_requests                  = 0
        requested_derivative_order        = 0
        output_excited_state_gradients    = 0.e0_wp
        output_ground_state_gradient      = 0.e0_wp
        output_diabatic_gradients         = 0.e0_wp
        output_excited_state_energies     = 0.e0_wp
        output_ground_state_energy        = 0.e0_wp
        output_diabatic_couplings         = 0.e0_wp
        excited_state_cartesian_gradients = 0.e0_wp
        ground_state_cartesian_gradient   = 0.e0_wp
        diabatic_cartesian_gradients      = 0.e0_wp
        last_requested_excited_state      = 0
        last_requested_coupling           = 0
        first_requested_excited_state     = 0
        first_requested_coupling          = 0
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 1
        coordinate_representation  = 2
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 4
        concom_de                  = [0.148201e0_wp, 0.0275690e0_wp, 0.151548e0_wp]
        concom_beta                = [1.260580e0_wp, 0.924180e0_wp, 1.068620e0_wp]
        concom_re                  = [1.863300e0_wp, 2.907700e0_wp, 1.428600e0_wp]
        concom_sato                = 0.10e0_wp
        concom_gam                 = [2.399700e0_wp, 1.058350e0_wp, 2.399700e0_wp]
        concom_reoh                = 1.808090e0_wp
        concom_rehh                = 2.861590e0_wp
        concom_con     = [-.0015920e0_wp, 0.026963e0_wp, 0.0014689e0_wp, 0.080011e0_wp, 0.085816e0_wp, -0.063179e0_wp, &
            & 0.101380e0_wp]
        concom_alp                 = [4.773e0_wp, 7.14e0_wp, 2.938e0_wp, 5.28e0_wp]
        concom_clam                = [0.10e0_wp, 0.10e0_wp, 0.20e0_wp, 0.03e0_wp]
        concom_acon                = [0.10e0_wp, 0.009e0_wp]
        v                          = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        requested_derivative_order = igrad
        surface_requests           = 1
        input_cartesian            = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                input_cartesian(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
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

        !> @brief Evaluate prepare_local_context for the OH3_LEPS_SE_1986 potential.
        subroutine prepare_local_context
            implicit none
            integer :: i
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            infocm_ref(1) = 'OH + H2 potential energy function'
            infocm_ref(2) = 'G. C. Schatz and H. Elgersma'
            infocm_ref(3) = 'Chem. Phys. Lett. 73, 21 (1980)'
        end subroutine prepare_local_context

        !> @brief Evaluate pot for the OH3_LEPS_SE_1986 potential.
        subroutine pot
            !   This subprogram evaluates the OH + H2 potential energy and derivativ
            !   of the potential with respect to the internal coordinates.
            !   The subprogram PREPEF must be called once before any calls to this
            !   subprogram.
            !   All calculations in this subprogram are in hartree atomic units.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: b, d, dvmor, rr, t, vmor
            integer  :: i
            !   Define the statement functions VMOR and DVMOR which evaluate the Mor
            !   diatomic potential energy  and the derivatives of the Morse potentia
            vmor(d, b, t, rr)                     = d*(1.0e0_wp-exp(-b*(rr-t)))**2
            dvmor(d, b, t, rr)                    = 2.0e0_wp*b*d*(1.0e0_wp-exp(-b*(rr-t)))*exp(-b*(rr-t))
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            do i = 1, 6
                potcm_r(i) = pt1cm_r(i)
            enddo
            !   Zero the array which will contain the derivatives of the energy
            !   with respect to the internal coordinates.
            if( requested_derivative_order .eq. 1 )then
                do i = 1, 6
                    potcm_dvdr(i) = 0.0e0_wp
                enddo
            endif
            !   Calculate the Morse part of the potential energy.
            potcm_vtot = vmor(concom_de(1), concom_beta(1), concom_re(1), potcm_r(1))+vmor(concom_de(2), &
                & concom_beta(2), concom_re(2), potcm_r(4 ))+vmor(concom_de(2), concom_beta(2), concom_re(2), &
                & potcm_r(5))
            !   Calculate the derivatives of the Morse part of the potential energy.
            if( requested_derivative_order .eq. 1 )then
                potcm_dvdr(1) = potcm_dvdr(1)+dvmor(concom_de(1), concom_beta(1), concom_re(1), potcm_r(1))
                potcm_dvdr(4) = potcm_dvdr(4)+dvmor(concom_de(2), concom_beta(2), concom_re(2), potcm_r(4))
                potcm_dvdr(5) = potcm_dvdr(5)+dvmor(concom_de(2), concom_beta(2), concom_re(2), potcm_r(5))
            endif
            !   Initialize the coordinates for the three-body LEPS part of the poten
            pot2cm_r(1) = potcm_r(2)
            pot2cm_r(2) = potcm_r(3)
            pot2cm_r(3) = potcm_r(6)
            !   Calculate the three-body LEPS portion of the potential and update th
            !   energy term.
            call v3pot
            potcm_vtot = potcm_vtot+pot2cm_energy
            !   Update the array containing the derivatives with the LEPS derivative
            if( requested_derivative_order .eq. 1 )then
                potcm_dvdr(2) = potcm_dvdr(2)+pot2cm_dedr(1)
                potcm_dvdr(3) = potcm_dvdr(3)+pot2cm_dedr(2)
                potcm_dvdr(6) = potcm_dvdr(6)+pot2cm_dedr(3)
            endif
            !   Initialize the coordinates for the H2O part of the potential for H1
            pot2cm_r(1) = potcm_r(1)
            pot2cm_r(2) = potcm_r(2)
            pot2cm_r(3) = potcm_r(4)
            !   Calculate the H2O part of the potential and update the energy term.
            call vh2o
            potcm_vtot = potcm_vtot+pot2cm_energy
            !   Update the array containing the derivatives with the H2O derivatives
            if( requested_derivative_order .eq. 1 )then
                potcm_dvdr(1) = potcm_dvdr(1)+pot2cm_dedr(1)
                potcm_dvdr(2) = potcm_dvdr(2)+pot2cm_dedr(2)
                potcm_dvdr(4) = potcm_dvdr(4)+pot2cm_dedr(3)
            endif
            !   Initialize the coordinates for the H2O part of the potential for H1
            pot2cm_r(1) = potcm_r(1)
            pot2cm_r(2) = potcm_r(3)
            pot2cm_r(3) = potcm_r(5)
            !   Calculate the H2O part of the potential and update the energy term.
            call vh2o
            potcm_vtot = potcm_vtot+pot2cm_energy
            !   Update the array containing the derivatives with the H2O derivatives
            if( requested_derivative_order .eq. 1 )then
                potcm_dvdr(1) = potcm_dvdr(1)+pot2cm_dedr(1)
                potcm_dvdr(3) = potcm_dvdr(3)+pot2cm_dedr(2)
                potcm_dvdr(5) = potcm_dvdr(5)+pot2cm_dedr(3)
            endif
            !   Initialize the coordinates for the four-body part of the potential.
            pot2cm_r(1) = potcm_r(2)
            pot2cm_r(2) = potcm_r(3)
            pot2cm_r(3) = potcm_r(4)
            pot2cm_r(4) = potcm_r(5)
            !   Calculate the four-body part of the potential and update the energy
            call v4pot
            potcm_vtot = potcm_vtot+pot2cm_energy
            !   Update the array containing the derivatives.
            if( requested_derivative_order .eq. 1 )then
                potcm_dvdr(2) = potcm_dvdr(2)+pot2cm_dedr(1)
                potcm_dvdr(3) = potcm_dvdr(3)+pot2cm_dedr(2)
                potcm_dvdr(4) = potcm_dvdr(4)+pot2cm_dedr(3)
                potcm_dvdr(5) = potcm_dvdr(5)+pot2cm_dedr(4)
            endif
            !   Adjust the potential to the correct zero, which corresponds to
            !   OH(R=RE) and H2(R=RE) at infinity.
            potcm_vtot          = potcm_vtot-2.0e0_wp*concom_de(2)+concom_de(3)
            ground_state_energy = potcm_vtot
            do i = 1, 6
                ground_state_internal_gradient(i) = potcm_dvdr(i)
            enddo
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate V3POT for the OH3_LEPS_SE_1986 potential.
        subroutine v3pot
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: b(3), beta(3), coul(3), de(3), do4z(3), exch(3), op3z(3), rad, re(3), s, top3z(3), &
                & tzp3(3), x(3), z(3), zp3(3), zpo(3)
            integer :: i
            !   Set up the values of DE, BETA, and RE for the three-body LEPS potent
            de(1)   = concom_de(1)
            beta(1) = concom_beta(1)
            re(1)   = concom_re(1)
            de(2)   = concom_de(1)
            beta(2) = concom_beta(1)
            re(2)   = concom_re(1)
            de(3)   = concom_de(3)
            beta(3) = concom_beta(3)
            re(3)   = concom_re(3)
            !   Compute useful constants.
            do i = 1, 3
                z(i)     = concom_sato
                zpo(i)   = 1.0e0_wp+z(i)
                op3z(i)  = 1.0e0_wp+3.0e0_wp*z(i)
                top3z(i) = 2.0e0_wp*op3z(i)
                zp3(i)   = z(i)+3.0e0_wp
                tzp3(i)  = 2.0e0_wp*zp3(i)
                do4z(i)  = de(i)/4.0e0_wp/zpo(i)
                b(i)     = beta(i)*do4z(i)*2.0e0_wp
            enddo
            !   Initialize the variable used for storing the energy.
            pot2cm_energy = 0.e0_wp
            do i = 1, 3
                x(i)          = exp(-beta(i)*(pot2cm_r(i)-re(i)))
                coul(i)       = do4z(i)*(zp3(i)*x(i)-top3z(i))*x(i)
                exch(i)       = do4z(i)*(op3z(i)*x(i)-tzp3(i))*x(i)
                pot2cm_energy = pot2cm_energy + coul(i)
            enddo
            rad           = sqrt((exch(1)-exch(2))**2+(exch(2)-exch(3))**2+(exch(3)-exch (1))**2)
            pot2cm_energy = pot2cm_energy - rad/sqrt(2.e0_wp)
            !   Compute the derivatives of the energy with respect to the internal
            !   coordinates.
            if( requested_derivative_order .eq. 1 )then
                s = exch(1) + exch(2) + exch(3)
                do i = 1, 3
                    pot2cm_dedr(i) = b(i)*x(i)*((3.0e0_wp*exch(i)-s)/sqrt(2.e0_wp)* (op3z(i)*x(i)-zp3(i))/rad- zp3(i)* &
                        & x(i)+ op3z(i))
                enddo
            endif
        end subroutine v3pot

        !> @brief Evaluate V4POT for the OH3_LEPS_SE_1986 potential.
        subroutine v4pot
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: t1, t2
            t1 = exp(-concom_clam(1)*(pot2cm_r(1)-concom_alp(1))**2-concom_clam(1)*(pot2cm_r(2)-concom_alp(1))**2- &
                & concom_clam(3)*(pot2cm_r(3)-concom_alp(3)) **2-concom_clam(3)*(pot2cm_r(4)-concom_alp(3))**2)* &
                & concom_acon(1)
            t2 = exp(-concom_clam(2)*(pot2cm_r(1)-concom_alp(2))**2-concom_clam(2)*(pot2cm_r(2)-concom_alp(2))**2- &
                & concom_clam(4)*(pot2cm_r(3)-concom_alp(4)) **2-concom_clam(4)*(pot2cm_r(4)-concom_alp(4))**2)* &
                & concom_acon(2)
            pot2cm_energy = t1+t2
            if( requested_derivative_order .eq. 1 )then
                pot2cm_dedr(1) = -2.0e0_wp*(t1*concom_clam(1)*(pot2cm_r(1)-concom_alp(1))+t2*concom_clam(2)* &
                    & (pot2cm_r(1)- concom_alp(2)))
                pot2cm_dedr(2) = -2.0e0_wp*(t1*concom_clam(1)*(pot2cm_r(2)-concom_alp(1))+t2*concom_clam(2)* &
                    & (pot2cm_r(2)- concom_alp(2)))
                pot2cm_dedr(3) = -2.0e0_wp*(t1*concom_clam(3)*(pot2cm_r(3)-concom_alp(3))+t2*concom_clam(4)* &
                    & (pot2cm_r(3)- concom_alp(4)))
                pot2cm_dedr(4) = -2.0e0_wp*(t1*concom_clam(3)*(pot2cm_r(4)-concom_alp(3))+t2*concom_clam(4)* &
                    & (pot2cm_r(4)- concom_alp(4)))
            endif
        end subroutine v4pot

        !> @brief Evaluate VH2O for the OH3_LEPS_SE_1986 potential.
        subroutine vh2o
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dp(3), dq(3), p, q(3), s(3), tanlg, temp, trm1, x(3), xmax1, xmax2, xx
            integer  :: i
            !     XMAX1 FOR WHEN TANH SET=1.0 ON VAX
            !     XMAX2 FOR PREVENTING OVERFLOWS ON VAX
            !     STATEMENT FUNCTION
            tanlg(xx) = 2.0e0_wp/(1.0e0_wp+exp(2.0e0_wp*xx))
            s(1)      = pot2cm_r(1)-concom_reoh
            s(3)      = pot2cm_r(2)-concom_reoh
            s(2)      = pot2cm_r(3)-concom_rehh
            xmax1     = 15.0e0_wp; xmax2 = 43.0e0_wp
            do i = 1, 3
                x(i) = 0.5e0_wp*concom_gam(i)*s(i)
                q(i) = 1.0e0_wp-tanh(x(i))
                if( x(i) < xmax1 )then
                    dq(i) = -0.5e0_wp*concom_gam(i)/cosh(x(i))**2
                elseif( x(i) < xmax2 )then
                    q(i)  = tanlg(x(i))
                    dq(i) = -0.5e0_wp*concom_gam(i)/cosh(x(i))**2
                else
                    q(i)  = 0.0e0_wp
                    dq(i) = 0.0e0_wp
                endif
            enddo
            p = concom_con(1)+concom_con(2)*(s(1)+s(3))+concom_con(3)*s(2)+0.5e0_wp*concom_con(4)*(s(1)*s(1)+s(3)* &
                & s( 3))+ 0.5e0_wp*concom_con(5)*s(2)*s(2)+concom_con(6)*s(2)*(s(1)+s(3))+ concom_con(7)*s(1)*s(3)
            pot2cm_energy = q(1)*q(2)*q(3)*p
            if( requested_derivative_order .eq. 1 )then
                dp(1) = concom_con(2)+concom_con(4)*s(1)+concom_con(6)*s(2)+concom_con(7)*s(3)
                dp(2) = concom_con(3)+concom_con(5)*s(2)+concom_con(6)*(s(1)+s(3))
                dp(3) = concom_con(2)+concom_con(4)*s(3)+concom_con(6)*s(2)+concom_con(7)*s(1)
                do i = 1, 3
                    trm1 = 0.0e0_wp
                    if( q(i) /= 0.0e0_wp ) trm1 = dq(i)/q(i)
                    pot2cm_dedr(i) = pot2cm_energy*(trm1+dp(i)/p)
                enddo
                temp           = pot2cm_dedr(2)
                pot2cm_dedr(2) = pot2cm_dedr(3)
                pot2cm_dedr(3) = temp
            endif
        end subroutine vh2o

        !> @brief Evaluate POTINFO for the OH3_LEPS_SE_1986 potential.
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
            if( kmax-1 .lt. maximum_excited_state )then
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

        !> @brief Evaluate CARTTOR for the OH3_LEPS_SE_1986 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(cartesian_coordinates, atom_count, coordinate_representation, pt1cm_r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the OH3_LEPS_SE_1986 potential.
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

        !> @brief Evaluate RTOCART for the OH3_LEPS_SE_1986 potential.
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
                            & 1))/pt1cm_r(m1)
                        termy                                = (cartesian_coordinates(i,2)-cartesian_coordinates(j, &
                            & 2))/pt1cm_r(m1)
                        termz                                = (cartesian_coordinates(i,3)-cartesian_coordinates(j, &
                            & 3))/pt1cm_r(m1)
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
                    ygs(i) = ground_state_internal_gradient(i)/pt1cm_r(i)
                    if( first_requested_excited_state .ne. 0 )then
                        do j=first_requested_excited_state, last_requested_excited_state
                            yes(i, j) = excited_state_internal_gradients(i,j)/pt1cm_r(i)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do k=first_requested_coupling, last_requested_coupling
                            yij(i, k) = diabatic_internal_gradients(i,k)/pt1cm_r(i)
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

        !> @brief Evaluate DEDCOU for the OH3_LEPS_SE_1986 potential.
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
end module chempotpy_oh3_leps_se_1986_pes

!> @brief Expose OH3_LEPS_SE_1986 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh3_leps_se_1986_pes, only: evaluate
    use chempotpy_kinds,                only: wp
    implicit none
    integer, parameter    :: natoms  = 4
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
