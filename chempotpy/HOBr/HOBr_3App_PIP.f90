!   System:                        OHBr_3App
!   Functional form:               permutationally invariant polynomials
!   Common name:                   OHBr (3A")
!   Number of derivatives:         1
!   Number of bodies:              3
!   Number of electronic surfaces: 1
!   Interface:                     Section-2
!   References: A. G. S. de Oliveira-Filho, F. R. Ornellas and K. A.
!   Peterson  "Accurate ab initio potential energy surfaces for the
!   3A" and 3A' electronic states of the O(3P)+HBr system"
!   J. Chem. Phys. 136, 174316 (2012).
!   Notes:
!      r1 = R(O-H)
!      r2 = R(H-Br)
!      r3 = R(Br-O)
!   Input: r1,r2,r3                     in units of bohr
!   Output: value                       in units of hartree
!   Output: dr1,dr2,dr3                 in units of hartree/bohr

!> @brief Implement HOBr_3App_PIP, a potential energy surface for HOBr.
module chempotpy_hobr_3app_pip_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the HOBr_3App_PIP potential energy surface for HOBr.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,    only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_geometry, only: cartesian_distance, flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_units,    only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds,    only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: v
        real(wp)              :: r_oh, r_hbr, r_obr
        real(wp)              :: dr_oh, dr_hbr, dr_obr
        real(wp)              :: r(3), tx(9)
        real(wp)              :: dr(3), drdx(3,9), dx(9)
        integer               :: iatom, idir, i, j, istate
        v     = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        r_obr = cartesian_distance(x(2,:), x(3,:))
        r_oh  = cartesian_distance(x(1,:), x(2,:))
        r_hbr = cartesian_distance(x(1,:), x(3,:))
        call flatten_cartesian_coordinates(x, tx)
        r = [r_oh, r_hbr, r_obr]
        call mypot(r_oh, r_hbr, r_obr, v, dr_oh, dr_hbr, dr_obr)
        v     = v*hartree_to_electronvolt
        dr(1) = dr_oh*hartree_to_electronvolt
        dr(2) = dr_hbr*hartree_to_electronvolt
        dr(3) = dr_obr*hartree_to_electronvolt
        !      dr(1)=dr_OH*hartree_per_bohr_to_ev_per_angstrom
        !      dr(2)=dr_HBr*hartree_per_bohr_to_ev_per_angstrom
        !      dr(3)=dr_OBr*hartree_per_bohr_to_ev_per_angstrom
        call evdrdx(tx, r, drdx)
        call transform_internal_gradient_to_cartesian(dr, drdx, dx)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = v
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = v
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    j               = (iatom-1)*3+idir
                    g(1,iatom,idir) = dx(j)
                enddo
            enddo
        elseif( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
    end subroutine evaluate
end module chempotpy_hobr_3app_pip_pes

!> @brief Expose HOBr_3App_PIP through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_hobr_3app_pip_pes, only: evaluate
    use chempotpy_kinds,             only: wp
    implicit none
    integer, parameter              :: natoms  = 3
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    use chempotpy_kinds,    only: wp
    implicit none
    integer, parameter    :: first_atoms(3)  = [1, 1, 2]
    integer, parameter    :: second_atoms(3) = [2, 3, 3]
    real(wp), intent(in)  :: x(9), r(3)
    real(wp), intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!> @brief Evaluate myprepot.
module chempotpy_hobr_3app_pip_mypot_legacy_state_fe33e804
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r1, r2, r3, value, dr1, dr2, dr3, initialize_only )
        use chempotpy_hobr_3app_parameters, only: amn => rkhs_coefficients, phi_hbr, phi_oh, x_1 => rkhs_x_1, x_2 => rkhs_x_2, y_1 => rkhs_y_1
        !..
        !.. The O(3P) + HBr --> OH + Br MRCI+Q/CBS(aug-cc-pVnZ(-PP); n = Q,5)+SO 3A"
        !.. potential surface of A. G. S. de Oliveira-Filho, F. R. Ornellas
        !.. and K. A. Peterson
        !..
        !.. Reference:
        !.. Antonio G. S. de Oliveira-Filho, Fernando R. Ornellas and Kirk A. Peterson
        !.. J. Chem. Phys. 136, 174316 (2012).
        !..
        !.. The Reproducing Kernel Hilbert Space method of Ho and Rabitz is used
        !.. to interpolate the surface consisting of 1110 geometries spanning
        !.. O-H-Br angles of 60-180 deg.  This surface does not contain the
        !.. H + BrO arrangement and no claims of accuracy are made for O-H-Br
        !.. angles smaller than 60 deg.
        !  USAGE:
        !  On Input:
        !      r1 = R(O-H)
        !      r2 = R(H-Br)
        !      r3 = R(Br-O).
        ! On Output:
        !      value = Energy, relative to the asymptotic reactants
        !              valley (O + HBr(r_e)), in hartree
        !      dr1 = Derivative of the potential with respect to r1 (R(O-H))
        !            in hartree/bohr
        !      dr2 = Derivative of the potential with respect to r2 (R(H-Br))
        !            in hartree/bohr
        !      dr3 = Derivative of the potential with respect to r3 (R(O-Br))
        !            in hartree/bohr
        !  NOTE:
        !    Before any actual potential energy calculations are made, a single
        !    call to myprepot must be made:
        !      call myprepot
        !    Later, the potential energy is computed by calling mypot:
        !      call mypot(r1,r2,r3,value,dr1,dr2,dr3)
        !   The released o_hbr_a_pp.par parameters are embedded in the parameter module
        !                              Saddle point
        !        r1                          r2                        r3
        !   2.6289099999999999        2.8486750000000001        5.1064215904894388
        !        energy
        !   7.9866812204093529E-003
        !        dr1                         dr2                       dr3
        !   1.5320730578638475E-007   1.5109533511581397E-006  -1.5577335148869720E-007
        !                             Reactants side vdW well
        !        r1                          r2                        r3
        !   4.3746320000000001        2.6855400000000000        7.0601719999996142
        !        energy
        !  -2.5958471320488563E-003
        !        dr1                         dr2                       dr3
        !   1.4650519388234480E-003   1.4648796113941348E-003  -1.4650322069311608E-003
        !                             Products  side vdW well
        !        r1                          r2                        r3
        !   1.8372599999999999        4.8222319999999996        4.4678939513580840
        !        energy
        !  -3.3549004949264938E-002
        !        dr1                         dr2                       dr3
        !   6.5041529057150577E-007  -4.8876198921465885E-010   6.6615884752874166E-009
        use chempotpy_kinds, only: wp
        implicit none
        logical, intent(in) :: initialize_only
        optional :: r1, r2, r3, value, dr1, dr2, dr3
        real(wp) :: r1, r2, r3, value, dr1, dr2, dr3
        real(wp) :: r_oh, r_hbr, r_bro
        real(wp) :: sum1, sum2, dzdr, dfdr
        real(wp) :: dbetadz, dedr_oh, dedr_hbr
        real(wp) :: v_oh, v_hbr, v_bro, v3
        real(wp) :: th_y, theta
        real(wp) :: z, re, de, phi_inf, alpha_s, r_s, z_i, phi_z, f_s
        real(wp) :: x_1_max, x_2_max, y_1_max
        real(wp) :: x_1_min, x_2_min, y_1_min
        real(wp) :: q_1, q_2, q_3, dq_1, dq_2, dq_3
        integer  :: i
        if( .not. initialize_only )then
            r_oh    = r1
            r_hbr   = r2
            r_bro   = r3
            re      = phi_oh(1)
            de      = phi_oh(2)
            phi_inf = phi_oh(3)
            alpha_s = phi_oh(4)
            r_s     = phi_oh(5)
            z       = (r_oh-re)/(r_oh+re)
            phi_z   = -phi_inf
            z_i     = 1.0_wp
            sum2    = 0.0_wp
            do i=6, 14, 1
                phi_z = phi_z+phi_oh(i)*z_i
                sum2  = sum2+(i-6)*phi_oh(i)*z_i/z
                z_i   = z_i*z
            enddo
            sum1    = phi_z
            f_s     = 1.0_wp/(exp(alpha_s*(r_oh-r_s))+1.0_wp)
            phi_z   = f_s*phi_z+phi_inf
            v_oh    = de*((1.0_wp-(re/r_oh)**6*exp(-phi_z*z))**2-1.0_wp)
            dzdr    = 2.0_wp*re/((r_oh+re)**2)
            dfdr    = -alpha_s*exp(alpha_s*(r_oh-r_s))*f_s*f_s
            dbetadz = phi_z+((1/dzdr)*dfdr*sum1+f_s*sum2)*z
            dedr_oh = 2.0_wp*de*(1.0_wp-(re/r_oh)**6*exp(-phi_z*z))*((re/r_oh)**6)*exp(-phi_z*z)*(dzdr*dbetadz+6.0_wp/ &
                & r_oh)
            re      = phi_hbr(1)
            de      = phi_hbr(2)
            phi_inf = phi_hbr(3)
            alpha_s = phi_hbr(4)
            r_s     = phi_hbr(5)
            z       = (r_hbr-re)/(r_hbr+re)
            phi_z   = -phi_inf
            z_i     = 1.0_wp
            sum2    = 0.0_wp
            do i=6, 14, 1
                phi_z = phi_z+phi_hbr(i)*z_i
                sum2  = sum2+(i-6)*phi_hbr(i)*z_i/z
                z_i   = z_i*z
            enddo
            sum1    = phi_z
            f_s     = 1.0_wp/(exp(alpha_s*(r_hbr-r_s))+1.0_wp)
            phi_z   = f_s*phi_z+phi_inf
            v_hbr   = de*((1.0_wp-(re/r_hbr)**6*exp(-phi_z*z))**2-1.0_wp)
            dzdr    = 2.0_wp*re/((r_hbr+re)**2)
            dfdr    = -alpha_s*exp(alpha_s*(r_hbr-r_s))*f_s*f_s
            dbetadz = phi_z+((1/dzdr)*dfdr*sum1+f_s*sum2)*z
            dedr_hbr = 2.0_wp*de*(1.0_wp-(re/r_hbr)**6*exp(-phi_z*z))*((re/r_hbr)**6)*exp(-phi_z*z)*(dzdr*dbetadz+ &
                & 6.0_wp/r_hbr)
            v_bro = 0.0868406_wp*exp(-1.91_wp*(r_bro-3.2509_wp))
            th_y  = (1.0_wp-((r_oh*r_oh+r_hbr*r_hbr-r_bro*r_bro)/(2.0_wp*r_oh*r_hbr)))/2.0_wp
            v3    = 0.0_wp
            dr1   = 0.0_wp
            dr2   = 0.0_wp
            dr3   = 0.0_wp
            do i=1, 1110, 1
                if( r_oh .ge. x_1(i) )then
                    q_1  = (1.0_wp/14.0_wp)*(r_oh**(-7))*(1.0_wp-(7.0_wp/9.0_wp)*(x_1(i)/r_oh))
                    dq_1 = 0.5_wp*(r_oh**(-8))*((8.0_wp/9.0_wp)*(x_1(i)/r_oh)-1.0_wp)
                else
                    q_1  = (1.0_wp/14.0_wp)*(x_1(i)**(-7))*(1.0_wp-(7.0_wp/9.0_wp)*(r_oh/x_1(i)))
                    dq_1 = -x_1(i)**(-8)/18.0_wp
                endif
                if( r_hbr .ge. x_2(i) )then
                    q_2  = (1.0_wp/14.0_wp)*(r_hbr**(-7))*(1.0_wp-(7.0_wp/9.0_wp)*(x_2(i)/r_hbr))
                    dq_2 = 0.5_wp*(r_hbr**(-8))*((8.0_wp/9.0_wp)*(x_2(i)/r_hbr)-1.0_wp)
                else
                    q_2  = (1.0_wp/14.0_wp)*(x_2(i)**(-7))*(1.0_wp-(7.0_wp/9.0_wp)*(r_hbr/x_2(i)))
                    dq_2 = -x_2(i)**(-8)/18.0_wp
                endif
                if( th_y .ge. y_1(i) )then
                    q_3  = 1.0_wp+th_y*y_1(i)+2.0_wp*y_1(i)**2*th_y*(1.0_wp-(1.0_wp/3.0_wp)*(y_1(i)/th_y))
                    dq_3 = y_1(i)+2.0_wp*y_1(i)**2
                else
                    q_3  = 1.0_wp+th_y*y_1(i)+2.0_wp*th_y**2*y_1(i)*(1.0_wp-(1.0_wp/3.0_wp)*(th_y/y_1(i)))
                    dq_3 = y_1(i)+4.0_wp*y_1(i)*th_y-2.0_wp*th_y**2
                endif
                v3  = v3+amn(i)*q_1*q_2*q_3
                dr1 = dr1+amn(i)*dq_1*q_2*q_3
                dr2 = dr2+amn(i)*q_1*dq_2*q_3
                dr3 = dr3+amn(i)*q_1*q_2*dq_3
            enddo
            value = 0.97714_wp*v3+v_oh+v_hbr+v_bro+0.14242006_wp
            dr3   = 0.97714_wp*dr3*r3/(2.0_wp*r1*r2)
            dr1   = 0.97714_wp*dr1-dr3*((r1-r2*((r_oh**2+r_hbr**2-r_bro**2)/(2*r_oh*r_hbr)))/r3)
            dr2   = 0.97714_wp*dr2-dr3*((r2-r1*((r_oh**2+r_hbr**2-r_bro**2)/(2*r_oh*r_hbr)))/r3)
            dr1   = dedr_oh+dr1
            dr2   = dedr_hbr+dr2
            dr3   = dr3-1.91_wp*0.0868406_wp*exp(-1.91_wp*(r_bro-3.2509_wp))
        endif
    end subroutine dispatch
end module chempotpy_hobr_3app_pip_mypot_legacy_state_fe33e804

!> @brief Compatibility interface for myprepot.

!> @brief Compatibility interface for mypot.
subroutine mypot( r1, r2, r3, value, dr1, dr2, dr3 )
    use chempotpy_hobr_3app_pip_mypot_legacy_state_fe33e804, only: dispatch
    use chempotpy_kinds,                                     only: wp
    implicit none
    real(wp) :: r1
    real(wp) :: r2
    real(wp) :: r3
    real(wp) :: value
    real(wp) :: dr1
    real(wp) :: dr2
    real(wp) :: dr3
    call dispatch(r1=r1, r2=r2, r3=r3, value=value, dr1=dr1, dr2=dr2, dr3=dr3, initialize_only=.false.)
end subroutine mypot
