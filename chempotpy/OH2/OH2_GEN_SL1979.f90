!   System:          OH2
!   Functional form:
!   Common name:     OH2SL
!   Reference:
!   Interface:                     potlib2001
!   Number of bodies:              3
!   Number of derivatives:         1
!   Number of electronic surfaces: 1
!   Notes: This is a collinear potential energy surface, i.e., the poten
!          energy is not defined for non-linear geometries.
!          The derivatives of the potential energy with respect to the
!          coordinates are determined numerically using a step size STEP
!          if R(1) and R(2) are less than 50 bohr.
!          The variable STEP is initialized in the block data subprogram
!          PTPACM.
!   prepare_local_context must be called once before any calls to POT.
!   The potential parameters are assigned in the Block Data Subprogram P
!   Coordinates, potential energy, and derivatives are passed
!   through the common block PT1CM:
!   The potential energy in the three asymptotic valleys are
!   stored in the common block ASYCM:
!   The potential energy in the AB valley, EASYAB, is equal to the poten
!   energy of the H "infinitely" far from the OH diatomic, with the
!   OH diatomic at its equilibrium configuration.  Similarly, the terms
!   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
!   respectively.
!   All the information passed through the common blocks PT1CM and ASYCM
!   is in Hartree atomic units.
!   This potential is written such that:
!                  R(1) = R(O-H)
!                  R(2) = R(H-H)
!                  R(3) = R(H-O)
!   The zero of energy is defined at O "infinitely" far from the H2 diat
!   The flags that indicate which calculations should be carried out in
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

!> @brief Implement OH2_GEN_SL1979, a potential energy surface for OH2.
module chempotpy_oh2_gen_sl1979_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the OH2_GEN_SL1979 potential energy surface for OH2.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
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
        real(wp)              :: asycm_easyab
        real(wp)              :: asycm_easyac
        real(wp)              :: asycm_easybc
        real(wp)              :: cartesian_coordinates(natom,3)
        integer               :: coordinate_representation
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: atom_count
        character(len=75)     :: infocm_ref(5)
        real(wp)              :: ndercm_step
        real(wp)              :: precm_a
        real(wp)              :: precm_al(3)
        real(wp)              :: precm_alf
        real(wp)              :: precm_b
        real(wp)              :: precm_c(50)
        real(wp)              :: precm_chi(3)
        real(wp)              :: precm_gam
        real(wp)              :: precm_s(3)
        real(wp)              :: precm_v(4)
        real(wp)              :: precm_vm(3)
        real(wp)              :: precm_x(3)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: ground_state_energy
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: satocm_bet(3)
        real(wp)              :: satocm_de(3)
        real(wp)              :: satocm_xre(3)
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
        asycm_easyab                      = 0.e0_wp; asycm_easyac = 0.e0_wp; asycm_easybc = 0.e0_wp; cartesian_coordinates = 0.e0_wp; coordinate_representation = 0
        maximum_derivative_order          = 0; maximum_excited_state = 0; atom_count = 0
        infocm_ref                        = ' '
        ndercm_step                       = 0.e0_wp
        precm_a                           = 0.e0_wp
        precm_al                          = 0.e0_wp
        precm_alf                         = 0.e0_wp
        precm_b                           = 0.e0_wp
        precm_c                           = 0.e0_wp
        precm_chi                         = 0.e0_wp
        precm_gam                         = 0.e0_wp
        precm_s                           = 0.e0_wp
        precm_v                           = 0.e0_wp
        precm_vm                          = 0.e0_wp
        precm_x                           = 0.e0_wp
        ground_state_internal_gradient    = 0.e0_wp
        ground_state_energy               = 0.e0_wp
        pt1cm_r                           = 0.e0_wp
        pt3cm_ezero                       = 0.e0_wp
        excited_state_internal_gradients  = 0.e0_wp
        excited_state_energies            = 0.e0_wp
        diabatic_internal_gradients       = 0.e0_wp
        diabatic_couplings                = 0.e0_wp
        satocm_bet                        = 0.e0_wp
        satocm_de                         = 0.e0_wp
        satocm_xre                        = 0.e0_wp
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
        requested_derivative_order = 0
        coordinate_representation  = 3
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 3
        satocm_de                  = [4.035e0_wp, 4.1662e0_wp, 4.035e0_wp]
        satocm_xre                 = [0.9923e0_wp, 0.7590e0_wp, 0.9923e0_wp]
        satocm_bet                 = [2.4113e0_wp, 2.0682e0_wp, 2.4113e0_wp]
        ndercm_step                = 5.0e-7_wp
        precm_a                    = 0.157002e-03_wp
        precm_b                    = -0.162301e-01_wp
        precm_alf                  = 0.9e0_wp
        precm_gam                  = 0.045e0_wp
        precm_al                   = [1.0e0_wp, 0.8e0_wp, 1.0e0_wp]
        precm_c = [0.354195e+01_wp, 0.302780e+01_wp, 0.577312e+00_wp, 0.539877e+01_wp, 0.233027e+01_wp, &
            & 0.157276e+02_wp, 0.701851e+01_wp, - 0.478535e+02_wp, -0.106915e+01_wp, 0.268137e+02_wp, - &
            & 0.104806e+02_wp, -0.182660e+02_wp, 0.833941e+01_wp, 0.872497e+02_wp, 0.917870e+01_wp, -0.697766e+02_wp, &
            & 0.434862e+02_wp, -0.156633e+02_wp, 0.380162e+02_wp, -0.796542e+01_wp, 0.602854e+02_wp, 0.500800e+02_wp, &
            & -0.545761e+02_wp, 0.410210e+02_wp, -0.890739e+01_wp, 0.191135e+02_wp, -0.512653e+02_wp, 0.365096e+02_wp, &
            & - 0.235972e+02_wp, -0.102322e+03_wp, 0.837770e+00_wp, 0.258200e+02_wp, -0.198105e+01_wp, - &
            & 0.484189e+01_wp, 0.105377e+02_wp, 0.372352e+01_wp, -0.207237e+02_wp, 0.223342e+02_wp, 0.244787e+01_wp, &
            & 0.979572e+01_wp, -0.135612e+01_wp, -0.724496e+01_wp, 0.234984e+02_wp, 0.136364e+03_wp, 0.253799e+01_wp, &
            & -0.203821e+02_wp, 0.699846e+01_wp, 0.198529e+02_wp, -0.475902e+01_wp, 0.236328e+01_wp]
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
                    g(1,iatom,idir) = output_ground_state_gradient(iatom,idir)*hartree_to_electronvolt
                enddo
            enddo
        elseif( igrad == 2 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives above first order are unavailable')
        endif

        contains

        !> @brief Evaluate prepare_local_context for the OH2_GEN_SL1979 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i, j
            real(wp), parameter :: hart = 27.21161e0_wp
            !   Echo potential data to the file linked to UNIT NFLAG(18)
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            j = 1
            !    Set the values of the classical energy in the three asymptotic vall
            asycm_easyab   = satocm_de(1)/hart
            asycm_easybc   = satocm_de(2)/hart
            asycm_easyac   = satocm_de(3)/hart
            pt3cm_ezero(1) = satocm_de(2)
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            infocm_ref(1) = 'No Reference'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the OH2_GEN_SL1979 potential.
        subroutine pot
            !     CALCULATE ENERGY AND DERIVATIVE
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the OH diatomic, with the
            !   OH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,
            !   respectively.
            !   This potential is written such that:
            !                  R(1) = R(O-H)
            !                  R(2) = R(H-H)
            !                  R(3) = R(H-O)
            !   The zero of energy is defined at O "infinitely" far from the H2 diat
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: ap, bp, check, d1d, d2d, d3d, da, da2, da3, dap1, dap2, dap3, db, db1, db2, db3, dc, dc2, dc3, dchi1, dchi2, dchi3, dcsr1, dcsr2, dcsr3, dd, dd2, dd3, df, df2, df3, dg, dq1, dq2, dq3, dvr1, dvr2, dvr3, gg, p, psi13, q, soh, sr2, t, temp, tt, u, w, ww, xter, y, z, zz
            integer             :: i, ic
            real(wp), parameter :: bohr = 0.529177e0_wp
            real(wp), parameter :: hart = 27.21161e0_wp
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   Check the values of surface_requests and requested_derivative_order for validity.
            ic = 4
            if( pt1cm_r(1) > 30.0e0_wp .or. pt1cm_r(2) > 30.0e0_wp ) ic = 1
            evaluation_loop: do
                select case( ic )
                case( 1 )
                    sr2        = pt1cm_r(2)
                    soh        = pt1cm_r(1)
                    pt1cm_r(3) = pt1cm_r(3)+ndercm_step
                case( 2 )
                    sr2 = pt1cm_r(2)+ndercm_step
                    soh = pt1cm_r(1)
                case( 3 )
                    sr2 = pt1cm_r(2)
                    soh = pt1cm_r(1)+ndercm_step
                case( 4 )
                    soh = pt1cm_r(1)
                    sr2 = pt1cm_r(2)
                case default
                    if( present(status) ) call status%set(status_invalid_input, &
                        'OH2 SL: invalid finite-difference state')
                    return
                end select
                sr2        = sr2*bohr
                soh        = soh*bohr
                pt1cm_r(3) = pt1cm_r(3)*bohr
                pt1cm_r(1) = soh
                pt1cm_r(2) = sr2
                q          = 1.0e0_wp
                do i=1, 3
                    precm_x(i)   = pt1cm_r(i)
                    precm_chi(i) = (1.0e0_wp - exp(-precm_alf*precm_x(i)**2))
                    precm_s(i)   = precm_x(i)-1.e0_wp
                    precm_vm(i)  = satocm_de(i)*(1.0e0_wp-exp(-satocm_bet(i)*(precm_x(i)-satocm_xre(i))))**2- &
                        & satocm_de(i)
                    q = q*(1.0e0_wp-tanh(precm_al(i)*precm_s(i)))
                enddo
                u = precm_vm(1)*precm_chi(2)*precm_chi(3)+precm_vm(2)*precm_chi(1)*precm_chi(3) + precm_vm(3)* &
                    & precm_chi(1)*precm_chi(2)
                bp    = precm_b*(1.0e0_wp/precm_x(1)**4+1.0e0_wp/precm_x(3)**4)
                check = (precm_x(1)*precm_x(3)*(precm_x(2)-0.2e0_wp))
                if( check > 1000.0e0_wp )then
                    ap = 0.0e0_wp
                else
                    ap = precm_a/(precm_x(1)*precm_x(3)*(precm_x(2)-0.2e0_wp))**12
                endif
                u = u+ap+bp
                t = precm_c(1) + precm_c(2)*(precm_s(1)+precm_s(3)) + precm_c(3)*precm_s(2) + precm_c(4)*(precm_s(1)** &
                    & 2 + precm_s(3)**2) + precm_c(5)*precm_s(2)**2 + precm_c(6)*precm_s(2)*(precm_s(1) + &
                    & precm_s(3)) + precm_c(7)* precm_s(1)*precm_s(3) + precm_c(8)*(precm_s(1)**3 + precm_s(3)**3) + &
                    & precm_c(9)* precm_s(2)**3 + precm_c(10)*precm_s(1)*precm_s(3)*(precm_s(1) + precm_s(3))
                tt = precm_c(11)*precm_s(2)**2*(precm_s(1) + precm_s(3)) + precm_c(12)*precm_s(2)*(precm_s(1)**2 + &
                    & precm_s(3)**2) + precm_c(13)*precm_s(1)*precm_s(2)*precm_s(3) + precm_c(14)*(precm_s(1)**4 + &
                    & precm_s(3)** 4) + precm_c(15)*precm_s(2)**4 + precm_c(16)*(precm_s(1)**3*precm_s(3) + &
                    & precm_s(1)* precm_s(3)**3)
                ww = precm_c(17)*precm_s(1)**2*precm_s(3)**2 + precm_c(18)*precm_s(2)*(precm_s(1)**3 + precm_s(3)** &
                    & 3) + precm_c(19)*precm_s(2)**2*(precm_s(1)**2 + precm_s(3)**2) + precm_c(20)*precm_s(2)**3* &
                    & (precm_s(1) + precm_s(3)) + precm_c(21)*precm_s(2)**2*precm_s(1)*precm_s(3) + precm_c(22)* &
                    & precm_s(2)*(precm_s(1)** 2*precm_s(3) + precm_s(1)*precm_s(3)**2)
                w = precm_c(23)*(precm_s(1)**5 + precm_s(3)**5) + precm_c(24)*(precm_s(1)**4*precm_s(3) + precm_s(1)* &
                    & precm_s(3)**4) + precm_c(25)*(precm_s(1)**3*precm_s(3)**2 + precm_s(1)**2*precm_s(3)**3) + &
                    & precm_c(26)*precm_s(2)*(precm_s(1)**4 + precm_s(3)**4) + precm_c(27)*precm_s(2)*(precm_s(1)**3* &
                    & precm_s(3) + precm_s(1)*precm_s(3)**3)
                gg = precm_c(28)*precm_s(2)*precm_s(1)**2*precm_s(3)**2 + precm_c(29)*precm_s(2)**2*(precm_s(1)**3 + &
                    & precm_s(3)**3) + precm_c(30)*precm_s(2)**2*(precm_s(1)**2*precm_s(3) + precm_s(1)*precm_s(3)** &
                    & 2) + precm_c(31)*precm_s(2)**3*(precm_s(1)**2 + precm_s(3)**2)
                y = precm_c(32)*precm_s(2)**3*precm_s(1)*precm_s(3) + precm_c(33)*precm_s(2)**4*(precm_s(1) + &
                    & precm_s(3)) + precm_c(34)*precm_s(2)**5 + precm_c(35)*(precm_s(1)**6 + precm_s(3)**6) + &
                    & precm_c(36)* (precm_s(1)**5*precm_s(3) + precm_s(1)*precm_s(3)**5) + precm_c(37)*(precm_s(1)**4* &
                    & precm_s(3)**2 + precm_s(1)**2*precm_s(3)**4) + precm_c(38)*precm_s(1)**3*precm_s(3)**3
                z = precm_c(39)*precm_s(2)*(precm_s(1)**5 + precm_s(3)**5) + precm_c(40)*precm_s(2)*(precm_s(1)**4* &
                    & precm_s(3) + precm_s(1)*precm_s(3)**4) + precm_c(41)*precm_s(2)*(precm_s(1)**3*precm_s(3)**2 + &
                    & precm_s(1)** 2*precm_s(3)**3) + precm_c(42)*precm_s(2)**2*(precm_s(1)**4 + precm_s(3)**4) + &
                    & precm_c(43)*precm_s(2)**2*(precm_s(1)**3*precm_s(3) + precm_s(1)*precm_s(3)**3) + precm_c(44)* &
                    & precm_s(1)**2* precm_s(2)**2*precm_s(3)**2
                zz = precm_c(45)*precm_s(2)**3*(precm_s(1)**3 + precm_s(3)**3) + precm_c(46)*precm_s(2)**3* &
                    & (precm_s(1)**2* precm_s(3) + precm_s(1)*precm_s(3)**2) + precm_c(47)*precm_s(2)**4*(precm_s(1)** &
                    & 2 + precm_s(3)**2) + precm_c(48)*precm_s(2)**4*precm_s(1)*precm_s(3) + precm_c(49)*precm_s(2)** &
                    & 5*(precm_s(1) + precm_s(3)) + precm_c(50)*precm_s(2)**6
                p           = t + tt + ww + w + gg + y + z + zz
                temp        = exp(-precm_gam*(precm_x(1) + precm_x(3))**2)
                u           = u + p*q*exp(-precm_gam*(precm_x(1) + precm_x(3))**2)
                xter        = u + pt3cm_ezero(1)
                precm_v(ic) = xter/hart
                if( ic < 4 )then
                    ic = ic+1
                    cycle evaluation_loop
                endif
                if( pt1cm_r(1) > 30.0e0_wp .or. pt1cm_r(2) >= 30.0e0_wp )then
                    if( requested_derivative_order == 1 )then
                        ground_state_internal_gradient(1) = (precm_v(3)-precm_v(4))/ndercm_step
                        ground_state_internal_gradient(2) = (precm_v(2)-precm_v(4))/ndercm_step
                        ground_state_internal_gradient(3) = (precm_v(1)-precm_v(4))/ndercm_step
                    endif
                else
                    ! Analytic derivatives in the interaction region.
                    if( requested_derivative_order .eq. 1 )then
                        da = precm_c(2) + 2.0e0_wp*precm_s(1)*precm_c(4) + precm_s(2)*precm_c(6) + precm_s(3)* &
                            & precm_c(7) + 3.0e0_wp* precm_s(1)**2*precm_c(8) + (precm_s(3)**2 + 2.0e0_wp*precm_s(1)* &
                            & precm_s(3))*precm_c(10) + precm_s(2)**2*precm_c(11) + 2.0e0_wp*precm_s(1)*precm_s(2)* &
                            & precm_c(12) + precm_s(2)*precm_s(3)* precm_c(13) + 4.0e0_wp*precm_s(1)**3*precm_c(14)
                        db = (3.0e0_wp*precm_s(1)**2*precm_s(3) + precm_s(3)**3)*precm_c(16) + 2.0e0_wp*precm_s(1)* &
                            & precm_s(3)**2* precm_c(17) + 3.0e0_wp*precm_s(1)**2*precm_s(2)*precm_c(18) + 2.0e0_wp* &
                            & precm_s(1)* precm_s(2)**2* precm_c(19) + precm_s(2)**3*precm_c(20) + precm_s(2)**2* &
                            & precm_s(3)*precm_c(21) + (2.0e0_wp* precm_s(1)*precm_s(2)*precm_s(3) + precm_s(2)* &
                            & precm_s(3)**2)*precm_c(22)
                        dc = 5.0e0_wp*precm_s(1)**4*precm_c(23) + (4.0e0_wp*precm_s(1)**3*precm_s(3) + precm_s(3)**4)* &
                            & precm_c(24) + (3.0e0_wp*precm_s(1)**2*precm_s(3)**2 + 2.0e0_wp*precm_s(1)*precm_s(3)** &
                            & 3)* precm_c(25) + 4.0e0_wp* precm_s(1)**3*precm_s(2)*precm_c(26) + (3.0e0_wp*precm_s(2)* &
                            & precm_s(1)**2*precm_s(3) + precm_s(2)* precm_s(3)**3)*precm_c(27)
                        dd = 2.0e0_wp*precm_s(1)*precm_s(2)*precm_s(3)**2*precm_c(28) + 3.0e0_wp*precm_s(1)**2* &
                            & precm_s(2)** 2* precm_c(29) + (2.0e0_wp*precm_s(1)*precm_s(2)**2*precm_s(3) + &
                            & precm_s(2)**2*precm_s(3)** 2)* precm_c(30) + 2.0e0_wp*precm_s(1)*precm_s(2)**3* &
                            & precm_c(31) + precm_s(2)**3*precm_s(3)* precm_c(32) + precm_s(2)**4*precm_c(33) + &
                            & 6.0e0_wp*precm_s(1)**5*precm_c(35) + (5.0e0_wp* precm_s(1)** 4*precm_s(3) + precm_s(3)** &
                            & 5)*precm_c(36)
                        df = (4.0e0_wp*precm_s(1)**3*precm_s(3)**2 + 2.0e0_wp*precm_s(1)*precm_s(3)**4)*precm_c(37) + &
                            & 3.0e0_wp* precm_s(1)**2*precm_s(3)**3*precm_c(38) + 5.0e0_wp*precm_s(1)**4*precm_s(2)* &
                            & precm_c(39) + (4.0e0_wp* precm_s(2)*precm_s(1)**3*precm_s(3) + precm_s(2)*precm_s(3)** &
                            & 4)*precm_c(40) + (3.0e0_wp* precm_s(1)**2*precm_s(2)*precm_s(3)**2 + 2.0e0_wp* &
                            & precm_s(1)*precm_s(2)* precm_s(3)**3)* precm_c(41) + 4.0e0_wp*precm_s(1)**3*precm_s(2)** &
                            & 2*precm_c(42)
                        dg = (3.0e0_wp*precm_s(1)**2*precm_s(2)**2*precm_s(3) + precm_s(2)**2*precm_s(3)**3)* &
                            & precm_c(43) + 2.0e0_wp* precm_s(1)*precm_s(2)**2*precm_s(3)**2*precm_c(44) + 3.0e0_wp* &
                            & precm_s(1)**2* precm_s(2)**3* precm_c(45) + (2.0e0_wp*precm_s(1)*precm_s(2)**3* &
                            & precm_s(3) + precm_s(2)**3*precm_s(3)**2)* precm_c(46) + 2.0e0_wp*precm_s(1)* &
                            & precm_s(2)**4*precm_c(47) + precm_s(2)**4* precm_s(3)* precm_c(48) + precm_s(2)**5* &
                            & precm_c(49)
                        d1d = da + db + dc + dd + df + dg
                        da2 = precm_c(3) + 2.0e0_wp*precm_s(2)*precm_c(5) + (precm_s(1) + precm_s(3))*precm_c(6) + &
                            & 3.0e0_wp* precm_s(2)**2*precm_c(9) + 2.0e0_wp*precm_s(2)*(precm_s(1) + precm_s(3))* &
                            & precm_c(11) + (precm_s(1)**2 + precm_s(3)**2)*precm_c(12) + precm_s(1)*precm_s(3)* &
                            & precm_c(13) + 4.0e0_wp* precm_s(2)**3*precm_c(15) + (precm_s(1)**3 + precm_s(3)**3)* &
                            & precm_c(18) + 2.0e0_wp* precm_s(2)* (precm_s(1)**2 + precm_s(3)**2)*precm_c(19)
                        db2 = 3.0e0_wp*precm_s(2)**2*(precm_s(1) + precm_s(3))*precm_c(20) + 2.0e0_wp*precm_s(1)* &
                            & precm_s(2)* precm_s(3)*precm_c(21) + (precm_s(1)**2*precm_s(3) + precm_s(1)*precm_s(3)** &
                            & 2)* precm_c(22) + (precm_s(1)**4 + precm_s(3)**4)*precm_c(26) + (precm_s(1)**3* &
                            & precm_s(3) + precm_s(1)* precm_s(3)**3)*precm_c(27) + precm_s(1)**2*precm_s(3)**2* &
                            & precm_c(28) + 2.0e0_wp*precm_s(2)* (precm_s(1)**3 + precm_s(3)**3)*precm_c(29)
                        dc2 = 2.0e0_wp*precm_s(2)*(precm_s(1)**2*precm_s(3) + precm_s(1)*precm_s(3)**2)*precm_c(30) + &
                            & 3.0e0_wp* precm_s(2)**2*(precm_s(1)**2 + precm_s(3)**2)*precm_c(31) + 3.0e0_wp* &
                            & precm_s(2)** 2* precm_s(1)* precm_s(3)*precm_c(32) + 4.0e0_wp*precm_s(2)**3* &
                            & (precm_s(1) + precm_s(3))* precm_c(33) + 5.0e0_wp* precm_s(2)**4*precm_c(34)
                        dd2 = (precm_s(1)**5 + precm_s(3)**5)*precm_c(39) + (precm_s(1)**4*precm_s(3) + precm_s(1)* &
                            & precm_s(3)**4)*precm_c(40) + (precm_s(1)**3*precm_s(3)**2 + precm_s(1)**2*precm_s(3)** &
                            & 3)* precm_c(41) + 2.0e0_wp*precm_s(2)*(precm_s(1)**4 + precm_s(3)**4)*precm_c(42) + &
                            & 2.0e0_wp* precm_s(2)* (precm_s(1)**3*precm_s(3) + precm_s(1)*precm_s(3)**3)* &
                            & precm_c(43) + 2.0e0_wp* precm_s(1)**2* precm_s(2)*precm_s(3)**2*precm_c(44)
                        df2 = 3.0e0_wp*precm_s(2)**2*(precm_s(1)**3 + precm_s(3)**3)*precm_c(45) + 3.0e0_wp* &
                            & precm_s(2)**2* (precm_s(1)**2*precm_s(3) + precm_s(1)*precm_s(3)**2)*precm_c(46) + &
                            & 4.0e0_wp* precm_s(2)** 3* (precm_s(1)**2 + precm_s(3)**2)*precm_c(47) + 4.0e0_wp* &
                            & precm_s(2)**3*precm_s(1)* precm_s(3)* precm_c(48) + 5.0e0_wp*precm_s(2)**4*(precm_s(1) + &
                            & precm_s(3))*precm_c(49) + 6.0e0_wp* precm_s(2)**5* precm_c(50)
                        d2d = da2 + db2 + dc2 + dd2 + df2
                        da3 = precm_c(2) + 2.0e0_wp*precm_s(3)*precm_c(4) + precm_s(2)*precm_c(6) + precm_s(1)* &
                            & precm_c(7) + 3.0e0_wp*precm_s(3)**2*precm_c(8) + (2.0e0_wp*precm_s(1)*precm_s(3) + &
                            & precm_s(1)**2)* precm_c(10) + precm_s(2)**2*precm_c(11) + 2.0e0_wp*precm_s(2)* &
                            & precm_s(3)* precm_c(12) + precm_s(1)*precm_s(2)* precm_c(13) + 4.0e0_wp*precm_s(3)**3* &
                            & precm_c(14)
                        db3 = (precm_s(1)**3 + 3.0e0_wp*precm_s(1)*precm_s(3)**2)*precm_c(16) + 2.0e0_wp*precm_s(1)** &
                            & 2* precm_s(3)* precm_c(17) + 3.0e0_wp*precm_s(2)*precm_s(3)**2*precm_c(18) + 2.0e0_wp* &
                            & precm_s(2)** 2*precm_s(3)* precm_c(19) + precm_s(2)**3*precm_c(20) + precm_s(2)**2* &
                            & precm_s(1)*precm_c(21) + precm_s(2)* (precm_s(1)**2 + 2.0e0_wp*precm_s(1)*precm_s(3))* &
                            & precm_c(22) + 5.0e0_wp*precm_s(3)** 4*precm_c(23) + (precm_s(1)**4 + 4.0e0_wp* &
                            & precm_s(1)* precm_s(3)**3)*precm_c(24)
                        dc3 = (2.0e0_wp*precm_s(1)**3*precm_s(3) + 3.0e0_wp*precm_s(1)**2*precm_s(3)**2)*precm_c(25) + &
                            & 4.0e0_wp* precm_s(2)*precm_s(3)**3*precm_c(26) + precm_s(2)*(precm_s(1)**3 + 3.0e0_wp* &
                            & precm_s(1)* precm_s(3)**2)*precm_c(27) + 2.0e0_wp*precm_s(2)*precm_s(1)**2*precm_s(3)* &
                            & precm_c(28) + 3.0e0_wp* precm_s(2)**2*precm_s(3)**2*precm_c(29)
                        dd3 = precm_s(2)**2*(precm_s(1)**2 + 2.0e0_wp*precm_s(1)*precm_s(3))*precm_c(30) + 2.0e0_wp* &
                            & precm_s(2)**3* precm_s(3)*precm_c(31) + precm_s(2)**3*precm_s(1)*precm_c(32) + &
                            & precm_s(2)**4* precm_c(33) + 6.0e0_wp*precm_s(3)**5*precm_c(35) + (precm_s(1)**5 + &
                            & 5.0e0_wp* precm_s(1)*precm_s(3)**4)*precm_c(36) + (2.0e0_wp*precm_s(1)**4*precm_s(3) + &
                            & 4.0e0_wp* precm_s(1)** 2* precm_s(3)**3)*precm_c(37) + 3.0e0_wp* precm_s(1)**3* &
                            & precm_s(3)**2* precm_c(38)
                        df3 = 5.0e0_wp*precm_s(2)*precm_s(3)**4*precm_c(39) + precm_s(2)*(precm_s(1)**4 + 4.0e0_wp* &
                            & precm_s(1)* precm_s(3)**3)*precm_c(40) + precm_s(2)*(2.0e0_wp*precm_s(1)**3*precm_s(3) + &
                            & 3.0e0_wp* precm_s(1)**2* precm_s(3)**2)*precm_c(41) + 4.0e0_wp*precm_s(2)**2* &
                            & precm_s(3)**3* precm_c(42) + precm_s(2)**2* (precm_s(1)**3 + 3.0e0_wp*precm_s(1)* &
                            & precm_s(3)**2)* precm_c(43) + 2.0e0_wp*precm_s(1)**2*precm_s(2)** 2*precm_s(3)* &
                            & precm_c(44) + 3.0e0_wp* precm_s(2)**3*precm_s(3)**2*precm_c(45) + precm_s(2)**3* &
                            & (precm_s(1)**2 + 2.0e0_wp* precm_s(1)*precm_s(3))* precm_c(46) + 2.0e0_wp*precm_s(2)**4* &
                            & precm_s(3)* precm_c(47) + precm_s(2)**4*precm_s(1)*precm_c(48) + precm_s(2)**5*precm_c(49)
                        d3d   = da3 + db3 + dc3 + dd3 + df3
                        dq1   = (1.0e0_wp - tanh(precm_al(2)*precm_s(2)))*(1.0e0_wp-tanh(precm_al(3)*precm_s(3)))
                        dq2   = (1.0e0_wp - tanh(precm_al(1)*precm_s(1)))*(1.0e0_wp-tanh(precm_al(3)*precm_s(3)))
                        dq3   = (1.0e0_wp - tanh(precm_al(1)*precm_s(1)))*(1.0e0_wp-tanh(precm_al(2)*precm_s(2)))
                        dap1  = -12.0e0_wp*ap/pt1cm_r(1)
                        dap2  = -12.0e0_wp*ap/(pt1cm_r(2) - 0.2e0_wp)
                        dap3  = -12.0e0_wp*ap/pt1cm_r(3)
                        db1   = -4.0e0_wp*precm_b/pt1cm_r(1)**5
                        db3   = -4.0e0_wp*precm_b/pt1cm_r(3)**5
                        dcsr1 = -precm_al(1)/(cosh(precm_al(1)*precm_s(1))**2)
                        !144   IF(COSH(AL(2)*S(2)).LT.1.0D-15)DCRS2=1.0D30
                        dcsr2 = -precm_al(2)/(cosh(precm_al(2)*precm_s(2))**2)
                        !145   IF(COSH(AL(3)*S(3)).LT.1.0D-15)DCRS3=1.0D30
                        dcsr3 = -precm_al(3)/(cosh(precm_al(3)*precm_s(3))**2)
                        psi13 = 2.0e0_wp*(-precm_gam)*temp*(pt1cm_r(1) + pt1cm_r(3))
                        dchi1 = 2.0e0_wp*precm_alf*pt1cm_r(1)*exp(-precm_alf*(pt1cm_r(1)**2))
                        dchi2 = 2.0e0_wp*precm_alf*pt1cm_r(2)*exp(-precm_alf*(pt1cm_r(2)**2))
                        dchi3 = 2.0e0_wp*precm_alf*pt1cm_r(3)*exp(-precm_alf*(pt1cm_r(3)**2))
                        dvr1 = 2.0e0_wp*satocm_de(1)*(1.0e0_wp - exp(-satocm_bet(1)*(pt1cm_r(1)-satocm_xre(1))))* &
                            & satocm_bet(1)* exp(-satocm_bet(1)*(pt1cm_r(1)-satocm_xre(1)))
                        ! 153  IF(VM(2).GT.1.0D20)DVR2=1.0D30
                        dvr2 = 2.0e0_wp*satocm_de(2)*(1.0e0_wp - exp(-satocm_bet(2)*(pt1cm_r(2)-satocm_xre(2))))* &
                            & satocm_bet(2)* exp(-satocm_bet(2)*(pt1cm_r(2)-satocm_xre(2)))
                        ! 154  IF(VM(3).GT.1.0D20)DVR3=1.0D30
                        dvr3 = 2.0e0_wp*satocm_de(3)*(1.0e0_wp - exp(-satocm_bet(3)*(pt1cm_r(3)-satocm_xre(3))))* &
                            & satocm_bet(3)* exp(-satocm_bet(3)*(pt1cm_r(3)-satocm_xre(3)))
                        ground_state_internal_gradient(1) = precm_vm(2)*precm_chi(3)*dchi1 + precm_chi(2)* &
                            & precm_chi(3)*dvr1 + precm_vm(3)*precm_chi(2)*dchi1 + d1d*q*temp + p*dq1*temp*dcsr1 + &
                            & dap1 + db1 + p*q*psi13
                        ground_state_internal_gradient(2) = dvr2*precm_chi(1)*precm_chi(3) + precm_vm(1)*precm_chi(3)* &
                            & dchi2 + precm_vm(3)* precm_chi(1)*dchi2 + d2d*q*temp + p*dq2*temp*dcsr2 + dap2
                        ground_state_internal_gradient(3) = precm_vm(2)*precm_chi(1)*dchi3 + precm_vm(1)*precm_chi(2)* &
                            & dchi3 + dvr3* precm_chi(1)* precm_chi(2) + d3d*q*temp + p*dq3*temp*dcsr3 + p* q* psi13 + &
                            & dap3 + db3
                        do i=1, 3
                            ground_state_internal_gradient(i) = ground_state_internal_gradient(i)/hart*bohr
                        enddo
                        ! END OF requested_derivative_order SWITCH AROUND DERIVATIVES
                    endif
                endif
                if( pt1cm_r(2) >= 20.0e0_wp )then
                    precm_v(4) = int(precm_v(4)*1.0e8_wp)*1.0e-8_wp
                endif
                exit evaluation_loop
            enddo evaluation_loop
            ground_state_energy = precm_v(4)
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
            return
            if( present(status) ) call status%set(status_invalid_input, &
                'legacy kernel rejected the evaluation request')
        end subroutine pot

        !> @brief Evaluate POTINFO for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate CARTTOR for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate EUNITZERO for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate RTOCART for the OH2_GEN_SL1979 potential.
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

        !> @brief Evaluate DEDCOU for the OH2_GEN_SL1979 potential.
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
end module chempotpy_oh2_gen_sl1979_pes

!> @brief Expose OH2_GEN_SL1979 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh2_gen_sl1979_pes, only: evaluate
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
