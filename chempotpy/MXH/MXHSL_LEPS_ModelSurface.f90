!   System:                        MXH
!   Common name:                   MXH-SL
!   Functional form:               Modified extended LEPS
!   Number of derivatives:         1
!   Number of electronic surfaces: 2
!   Interface:                     3-2V
!   Reference:  Y. L. Volobuev, M. D. Hack, M. S. Topaler, and D. G. Tru
!               J. Chem. Phys. 112 9716 (2000).
!   Notes:
!   Protocol: POT evaluates the stateless diabatic surfaces.
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      ground and first excited electronic states, diabatic representati
!   Zero of energy:
!      The classical potential energy is set equal to zero for M
!      infinitely far from the HX diatomic and R(HX) set equal to the
!      HX equilibrium diatomic value.
!   Parameters:
!   Coordinates:
!      Internal, Definition: R(1) = R(M-H)
!                            R(2) = R(H-X)
!                            R(3) = R(M-X))

!> @brief Implement MXHSL_LEPS_ModelSurface, a potential energy surface for MXH.
module chempotpy_mxhsl_leps_modelsurface_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the MXHSL_LEPS_ModelSurface potential energy surface for MXH.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,              only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_units,              only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: triatomic_distances_to_bohr
        use chempotpy_kinds,              only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 2
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        integer               :: iatom, idir, nder
        real(wp)              :: nt, r(1,3), r2(3), v(1)
        real(wp)              :: u11(1), u12(1), u22(1)
        real(wp)              :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
        real(wp)              :: u(nstates,nstates), dudr(nstates,nstates,3)
        real(wp)              :: t(nstates,nstates)
        real(wp)              :: tmpmat(nstates,nstates)
        real(wp)              :: hr(nstates,nstates,3), gr(nstates,3)
        real(wp)              :: tx(9), drdx(3,9)
        real(wp)              :: hx(nstates,nstates,9), gx(nstates,9)
        integer               :: i, j, k, l
        !initialize
        u  = 0.e0_wp; p = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        nt = 1
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is ClHH
        call triatomic_distances_to_bohr(x, r(1,:))
        nder = igrad
        call pot(r, u11, de_u11, 1, 1)
        call pot(r, u12, de_u12, 1, 2)
        call pot(r, u22, de_u22, 1, 3)
        u(1,1)      = u11(1)*hartree_to_electronvolt
        u(1,2)      = u12(1)*hartree_to_electronvolt
        u(2,1)      = u12(1)*hartree_to_electronvolt
        u(2,2)      = u22(1)*hartree_to_electronvolt
        dudr(1,1,:) = de_u11(:,1)*hartree_per_bohr_to_ev_per_angstrom
        dudr(1,2,:) = de_u12(:,1)*hartree_per_bohr_to_ev_per_angstrom
        dudr(2,1,:) = de_u12(:,1)*hartree_per_bohr_to_ev_per_angstrom
        dudr(2,2,:) = de_u22(:,1)*hartree_per_bohr_to_ev_per_angstrom
        call diagonalize(nstates, u, p, t)
        do i=1, 3
            tmpmat(:,:) = dudr(:,:,i)
            tmpmat      = matmul(transpose(t), matmul(tmpmat, t))
            do j=1, nstates
                do k=j+1, nstates
                    if( (p(k)-p(j)) > 1.e-8_wp )then
                        hr(j,k,i) = tmpmat(j,k)/(p(k)-p(j))
                    else
                        hr(j,k,i) = tmpmat(j,k)/1.e-8_wp
                    endif
                    hr(k,j,i) = -hr(j,k,i)
                enddo
            enddo
            do j=1, nstates
                gr(j,i) = tmpmat(j,j)
            enddo
        enddo
        r2 = r(1,:)*bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(hr, drdx, hx)
        call transform_internal_gradient_to_cartesian(gr, drdx, gx)
        do iatom=1, natoms
            do idir=1, 3
                j                 = (iatom-1)*3+idir
                g(:,iatom,idir)   = gx(:,j)
                d(:,:,iatom,idir) = hx(:,:,j)
            enddo
        enddo
    end subroutine evaluate
end module chempotpy_mxhsl_leps_modelsurface_pes

!> @brief Expose MXHSL_LEPS_ModelSurface through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_mxhsl_leps_modelsurface_pes, only: evaluate
    use chempotpy_kinds,                       only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX for the MXHSL_LEPS_ModelSurface potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    use chempotpy_kinds,    only: wp
    implicit none
    integer, parameter    :: first_atoms(3)  = [1, 2, 1]
    integer, parameter    :: second_atoms(3) = [2, 3, 3]
    real(wp), intent(in)  :: x(9), r(3)
    real(wp), intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!> @brief Evaluate diagonalize for the MXHSL_LEPS_ModelSurface potential.
subroutine diagonalize( n, a_ss, ev_s, u_ss )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)   :: n
    real(wp), intent(in)  :: a_ss(n,n)
    real(wp), intent(out) :: u_ss(n,n)
    real(wp), intent(out) :: ev_s(n)
    integer               :: io, i
    real(wp), allocatable :: work_d(:)
    integer               :: lwork_d
    u_ss    = a_ss
    lwork_d = -1
    allocate( work_d(1) )
    call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
    lwork_d = int(work_d(1))
    deallocate( work_d )
    allocate( work_d(lwork_d) )
    call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
end subroutine diagonalize

!> @brief Evaluate pot for the MXHSL_LEPS_ModelSurface potential.
subroutine pot( r, e, de, nt, nsurf )
    use chempotpy_kinds, only: wp
    implicit none
    integer, parameter :: g_p_  = 3
    integer, parameter :: ldg_r = 3
    integer, parameter :: ldg_e = 3
    real(wp)           :: ddis(3), de(3,nt), de1(3,1), de2(3,1), de3(3,1), e(nt), g_e(ldg_e,nt), g_r(ldg_r,nt,3), r(nt,3)
    integer            :: i, j, k, nsurf, nt
    ! initialize seed matrix here
    do i=1, nt
        do j=1, 3
            do k=1, 3
                if( j .eq. k )then
                    g_r(j,i,k) = 1.0e0_wp
                else
                    g_r(j,i,k) = 0.0e0_wp
                endif
            enddo
        enddo
    enddo
    if( nsurf .eq. 1 )then
        call g_pot11(g_p_, r, g_r, ldg_r, e, de, ldg_e, nt)
    elseif( nsurf .eq. 2 )then
        call g_pot12(g_p_, r, g_r, ldg_r, e, de, ldg_e, nt)
    elseif( nsurf .eq. 3 )then
        call g_pot22(g_p_, r, g_r, ldg_r, e, de, ldg_e, nt)
    else
        write(*, *)'Error!'
        stop
    endif
end subroutine pot

!  U11 surface - lower diabatic surface
!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_mxhsl_leps_modelsurface_g_pot11_legacy_state_62de5a21
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt
        real(wp) :: r(:,:), e(:), des1, bes1, res1, rehf, dehf, des3, bes3, res3
        real(wp) :: s1, s2, s3, xes1, xehf, xes3, yhf, behf
        real(wp) :: cevau, cosa
        integer  :: i, nt
        ! MH singlet
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d1_p, d5_b, d2_v, d3_v, d4_b, d4_v, d1_w, g_d1_w(g_pmax_), g_r(:,:,:), g_xes1(g_pmax_)
        real(wp)           :: g_s1(g_pmax_), g_yhf(g_pmax_), g_behf(g_pmax_), &
            & g_xehf(g_pmax_), g_s2(g_pmax_), g_xes3(g_pmax_), g_s3(g_pmax_), g_e(:,:)
        integer :: g_ehfid
        ! HF singlet
        ! MF singlet
        !  misc constants
        ! Common blocks for the diatomic potential
        behf    = 0.e0_wp
        bes1    = 0.e0_wp
        bes3    = 0.e0_wp
        cevau   = 0.e0_wp
        cosa    = 0.e0_wp
        d1_p    = 0.e0_wp
        d1_w    = 0.e0_wp
        d2_v    = 0.e0_wp
        d3_v    = 0.e0_wp
        d4_b    = 0.e0_wp
        d4_v    = 0.e0_wp
        d5_b    = 0.e0_wp
        dehf    = 0.e0_wp
        des1    = 0.e0_wp
        des3    = 0.e0_wp
        g_behf  = 0.e0_wp
        g_d1_w  = 0.e0_wp
        g_s1    = 0.e0_wp
        g_s2    = 0.e0_wp
        g_s3    = 0.e0_wp
        g_xehf  = 0.e0_wp
        g_xes1  = 0.e0_wp
        g_xes3  = 0.e0_wp
        g_yhf   = 0.e0_wp
        rehf    = 0.e0_wp
        res1    = 0.e0_wp
        res3    = 0.e0_wp
        s1      = 0.e0_wp
        s2      = 0.e0_wp
        s3      = 0.e0_wp
        xehf    = 0.e0_wp
        xes1    = 0.e0_wp
        xes3    = 0.e0_wp
        yhf     = 0.e0_wp
        des1    = 1.67e0_wp
        bes1    = 1.7e0_wp
        res1    = 1.5e0_wp
        rehf    = 1.733e0_wp
        dehf    = 6.122e0_wp
        des3    = 0.1e0_wp
        bes3    = 1.2e0_wp
        res3    = 4.0e0_wp
        g_ehfid = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        cevau = 1.e0_wp / 27.2113961e0_wp
        do i = 1, nt
            !     *(r(i, 2) * r(i, 3))
            !==============================================================k
            !This section contains the diatomic curves
            ! The  MH singlet
            res1 = 2.0e0_wp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bes1) * g_r(g_i_,i,1)
            enddo
            d1_w = (-bes1) * (r(i,1) - res1)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xes1(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xes1 = d2_v
            d2_v = des1 * xes1
            d3_v = xes1 + 2.0e0_wp
            d4_b = d2_v + d3_v * des1
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d4_b * g_xes1(g_i_)
            enddo
            s1 = d2_v * d3_v
            ! The HF singlet
            rehf = 1.733e0_wp
            dehf = 6.122e0_wp
            do g_i_ = 1, g_p_
                g_yhf(g_i_) = g_r(g_i_,i,2)
            enddo
            yhf  = r(i,2) - 2.1042e0_wp
            d4_v = yhf * yhf
            d1_p = 2.0e0_wp * yhf
            d5_b = 0.059062e0_wp * d1_p + (-0.025647e0_wp)
            do g_i_ = 1, g_p_
                g_behf(g_i_) = d5_b * g_yhf(g_i_)
            enddo
            behf = 1.1622e0_wp - 0.025647e0_wp * yhf + 0.059062e0_wp * d4_v
            d4_v = r(i,2) - rehf
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-behf) * g_r(g_i_,i,2) + (-d4_v) * g_behf( g_i_)
            enddo
            d1_w = (-behf) * d4_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xehf(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xehf = d2_v
            d2_v = dehf * xehf
            d3_v = xehf - 2.0e0_wp
            d4_b = d2_v + d3_v * dehf
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d4_b * g_xehf(g_i_)
            enddo
            s2 = d2_v * d3_v
            ! The MF singlet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bes3) * g_r(g_i_,i,3)
            enddo
            d1_w = (-bes3) * (r(i,3) - res3)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xes3(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xes3 = d2_v
            des3 = 0.060
            d2_v = des3 * xes3
            d3_v = xes3 - 2.0e0_wp
            d4_b = d2_v + d3_v * des3
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d4_b * g_xes3(g_i_)
            enddo
            s3 = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = g_s3(g_i_) + g_s2(g_i_) + g_s1(g_i_)
            enddo
            e(i) = s1 + s2 + s3 + dehf
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
        enddo
    end subroutine dispatch
end module chempotpy_mxhsl_leps_modelsurface_g_pot11_legacy_state_62de5a21

!> @brief Compatibility interface for g_pot11.
subroutine g_pot11( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_mxhsl_leps_modelsurface_g_pot11_legacy_state_62de5a21, only: dispatch
    use chempotpy_kinds,                                                 only: wp
    implicit none
    integer  :: g_p_
    real(wp) :: r(nt,3)
    real(wp) :: g_r(ldg_r,nt,3)
    integer  :: ldg_r
    real(wp) :: e(nt)
    real(wp) :: g_e(ldg_e,nt)
    integer  :: ldg_e
    integer  :: nt
    call dispatch(g_p_=g_p_, r=r, g_r=g_r, ldg_r=ldg_r, e=e, g_e=g_e, ldg_e=ldg_e, nt=nt)
end subroutine g_pot11

! U12 coupling surface
!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_mxhsl_leps_modelsurface_g_pot12_legacy_state_62de5a21
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt
        real(wp)           :: ap, alp, bet, rd0, r20, theta
        real(wp)           :: gc, gc2, rnag, r1s, r2s, r3s
        real(wp)           :: cevau
        integer            :: nt, i
        real(wp)           :: r(:,:), e(:)
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d2_p, d2_w, d8_b, d7_b, d6_b, d5_b, d4_b, d2_v, &
            & d4_v, d5_v
        real(wp) :: d1_p, d1_w, d3_v, g_r1s(g_pmax_), g_r(:,:,:), g_r2s(g_pmax_), g_r3s(g_pmax_), g_d1_w(g_pmax_), g_rnag(g_pmax_), g_gc(g_pmax_)
        real(wp) :: g_gc2(g_pmax_), g_d2_w(g_pmax_), g_e(:,:)
        integer  :: g_ehfid
        alp     = 0.e0_wp
        ap      = 0.e0_wp
        bet     = 0.e0_wp
        cevau   = 0.e0_wp
        d1_p    = 0.e0_wp
        d1_w    = 0.e0_wp
        d2_p    = 0.e0_wp
        d2_v    = 0.e0_wp
        d2_w    = 0.e0_wp
        d3_v    = 0.e0_wp
        d4_b    = 0.e0_wp
        d4_v    = 0.e0_wp
        d5_b    = 0.e0_wp
        d5_v    = 0.e0_wp
        d6_b    = 0.e0_wp
        d7_b    = 0.e0_wp
        d8_b    = 0.e0_wp
        g_d1_w  = 0.e0_wp
        g_d2_w  = 0.e0_wp
        g_gc    = 0.e0_wp
        g_gc2   = 0.e0_wp
        g_r1s   = 0.e0_wp
        g_r2s   = 0.e0_wp
        g_r3s   = 0.e0_wp
        g_rnag  = 0.e0_wp
        gc      = 0.e0_wp
        gc2     = 0.e0_wp
        r1s     = 0.e0_wp
        r20     = 0.e0_wp
        r2s     = 0.e0_wp
        r3s     = 0.e0_wp
        rd0     = 0.e0_wp
        rnag    = 0.e0_wp
        theta   = 0.e0_wp
        ap      = 0.3017e0_wp
        alp     = 3.0e0_wp
        rd0     = 3.0e0_wp
        bet     = 0.5e0_wp
        r20     = -1.3e0_wp
        g_ehfid = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        theta = 15.6e0_wp * acos(-1.0e0_wp) / 180.0e0_wp
        cevau = 1.e0_wp / 27.2113961e0_wp
        do i = 1, nt
            d2_v = r(i,1) * r(i,1)
            d1_p = 2.0e0_wp * r(i,1)
            do g_i_ = 1, g_p_
                g_r1s(g_i_) = d1_p * g_r(g_i_,i,1)
            enddo
            r1s  = d2_v
            d2_v = r(i,2) * r(i,2)
            d1_p = 2.0e0_wp * r(i,2)
            do g_i_ = 1, g_p_
                g_r2s(g_i_) = d1_p * g_r(g_i_,i,2)
            enddo
            r2s  = d2_v
            d2_v = r(i,3) * r(i,3)
            d1_p = 2.0e0_wp * r(i,3)
            do g_i_ = 1, g_p_
                g_r3s(g_i_) = d1_p * g_r(g_i_,i,3)
            enddo
            r3s  = d2_v
            d4_b = -(2.0e0_wp / 9.0e0_wp)
            d7_b = 1.0e0_wp / 3.0e0_wp
            d8_b = 2.0e0_wp / 3.0e0_wp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_r2s(g_i_) + d7_b * g_r1s(g_i_) + d8_b * g_r3s(g_i_)
            enddo
            d1_w = 2.0e0_wp / 3.0e0_wp * r3s + 1.0e0_wp / 3.0e0_wp * r1s - 2.0e0_wp / 9.0e0_wp * r2s
            d2_v = sqrt(d1_w)
            if( d1_w .gt. 0.0e0_wp )then
                d1_p = 1.0e0_wp / (2.0e0_wp *  d2_v)
            else
                d1_p = 0.0e0_wp
            endif
            do g_i_ = 1, g_p_
                g_rnag(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            rnag = d2_v
            d4_b = sin(theta)
            d5_b = cos(theta)
            do g_i_ = 1, g_p_
                g_gc(g_i_) = d4_b * g_r(g_i_,i,2) + d5_b * g_rnag(g_i_)
            enddo
            gc   = cos(theta) * rnag + sin(theta) * r(i,2)
            d4_b = -cos(theta)
            d5_b = sin(theta)
            do g_i_ = 1, g_p_
                g_gc2(g_i_) = d4_b * g_r(g_i_,i,2) + d5_b * g_rnag(g_i_)
            enddo
            gc2  = sin(theta) * rnag - cos(theta) * r(i,2)
            d4_v = ((gc - rd0) / alp) ** ( 4 - 2)
            d4_v = d4_v * ((gc - rd0) / alp)
            d1_p = 4 *  d4_v
            d4_v = d4_v * ((gc - rd0) / alp)
            d5_b = (-d1_p) * (1.0e0_wp / alp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_gc(g_i_)
            enddo
            d1_w = -d4_v
            d4_v = ((gc2 - r20) / bet) * ((gc2 - r20) / bet)
            d1_p = 2.0e0_wp * ((gc2 - r20) / bet)
            d5_b = (-d1_p) * (1.0e0_wp / bet)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = d5_b * g_gc2(g_i_)
            enddo
            d2_w = -d4_v
            d2_v = exp(d1_w)
            d2_p = d2_v
            d3_v = ap * d2_v
            d5_v = exp(d2_w)
            d1_p = d5_v
            d4_b = d3_v * d1_p
            d6_b = d5_v * ap * d2_p
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d4_b * g_d2_w(g_i_) + d6_b * g_d1_w(g_i_)
            enddo
            e(i) = d3_v * d5_v
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
        enddo
    end subroutine dispatch
end module chempotpy_mxhsl_leps_modelsurface_g_pot12_legacy_state_62de5a21

!> @brief Compatibility interface for g_pot12.
subroutine g_pot12( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_mxhsl_leps_modelsurface_g_pot12_legacy_state_62de5a21, only: dispatch
    use chempotpy_kinds,                                                 only: wp
    implicit none
    integer  :: g_p_
    real(wp) :: r(nt,3)
    real(wp) :: g_r(ldg_r,nt,3)
    integer  :: ldg_r
    real(wp) :: e(nt)
    real(wp) :: g_e(ldg_e,nt)
    integer  :: ldg_e
    integer  :: nt
    call dispatch(g_p_=g_p_, r=r, g_r=g_r, ldg_r=ldg_r, e=e, g_e=g_e, ldg_e=ldg_e, nt=nt)
end subroutine g_pot12

! the upper diabatic surface  U22
!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_mxhsl_leps_modelsurface_g_pot22_legacy_state_62de5a21
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt
        real(wp) :: r(:,:), e(:), denah, renah, benah, rehf, dehf, demf, remf, bemf
        real(wp) :: c2a, c2b, c2c, um2p, cevau, onethd, twothd
        real(wp) :: s1, s2, s3, t1, t2, t3, s2a, s2b, s2m
        real(wp) :: coul1, coul2, coul3, exch1, exch2, exch3, cplg2, w
        real(wp) :: xenah, xehfa, xehfb, xemf, sdehf
        real(wp) :: yhf, behf
        real(wp) :: cosa, rave
        integer  :: i, j, k, nt
        ! MH singlet & triplet
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d3_p, d2_p, d10_v, d10_b, d7_v, d5_b, d1_w, d7_b, d2_v, d3_v
        real(wp)           :: d6_b, d2_b, d3_b, d1_p, d4_v, d4_b, g_d1_w(g_pmax_), g_r(:,:,:), g_xenah(g_pmax_), g_s1(g_pmax_)
        real(wp)           :: g_t1(g_pmax_), g_yhf(g_pmax_), g_behf(g_pmax_), &
            & g_xehfa(g_pmax_), g_s2(g_pmax_), g_t2(g_pmax_), g_xemf(g_pmax_), &
            & g_t3(g_pmax_), g_s3(g_pmax_), g_coul1(g_pmax_)
        real(wp) :: g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_), g_exch3(g_pmax_), g_w(g_pmax_), g_cplg2(g_pmax_), g_e(:,:)
        integer  :: g_ehfid
        ! HF singlet & triplet
        ! LiF singlet & triplet
        ! H2 singlets
        ! misc constants
        behf    = 0.e0_wp
        bemf    = 0.e0_wp
        benah   = 0.e0_wp
        c2a     = 0.e0_wp
        c2b     = 0.e0_wp
        c2c     = 0.e0_wp
        cevau   = 0.e0_wp
        cosa    = 0.e0_wp
        coul1   = 0.e0_wp
        coul2   = 0.e0_wp
        coul3   = 0.e0_wp
        cplg2   = 0.e0_wp
        d10_b   = 0.e0_wp
        d10_v   = 0.e0_wp
        d1_p    = 0.e0_wp
        d1_w    = 0.e0_wp
        d2_b    = 0.e0_wp
        d2_p    = 0.e0_wp
        d2_v    = 0.e0_wp
        d3_b    = 0.e0_wp
        d3_p    = 0.e0_wp
        d3_v    = 0.e0_wp
        d4_b    = 0.e0_wp
        d4_v    = 0.e0_wp
        d5_b    = 0.e0_wp
        d6_b    = 0.e0_wp
        d7_b    = 0.e0_wp
        d7_v    = 0.e0_wp
        dehf    = 0.e0_wp
        demf    = 0.e0_wp
        denah   = 0.e0_wp
        exch1   = 0.e0_wp
        exch2   = 0.e0_wp
        exch3   = 0.e0_wp
        g_behf  = 0.e0_wp
        g_coul1 = 0.e0_wp
        g_coul2 = 0.e0_wp
        g_coul3 = 0.e0_wp
        g_cplg2 = 0.e0_wp
        g_d1_w  = 0.e0_wp
        g_exch1 = 0.e0_wp
        g_exch2 = 0.e0_wp
        g_exch3 = 0.e0_wp
        g_s1    = 0.e0_wp
        g_s2    = 0.e0_wp
        g_s3    = 0.e0_wp
        g_t1    = 0.e0_wp
        g_t2    = 0.e0_wp
        g_t3    = 0.e0_wp
        g_w     = 0.e0_wp
        g_xehfa = 0.e0_wp
        g_xemf  = 0.e0_wp
        g_xenah = 0.e0_wp
        g_yhf   = 0.e0_wp
        onethd  = 0.e0_wp
        rave    = 0.e0_wp
        rehf    = 0.e0_wp
        remf    = 0.e0_wp
        renah   = 0.e0_wp
        s1      = 0.e0_wp
        s2      = 0.e0_wp
        s2a     = 0.e0_wp
        s2b     = 0.e0_wp
        s2m     = 0.e0_wp
        s3      = 0.e0_wp
        sdehf   = 0.e0_wp
        t1      = 0.e0_wp
        t2      = 0.e0_wp
        t3      = 0.e0_wp
        twothd  = 0.e0_wp
        um2p    = 0.e0_wp
        w       = 0.e0_wp
        xehfa   = 0.e0_wp
        xehfb   = 0.e0_wp
        xemf    = 0.e0_wp
        xenah   = 0.e0_wp
        yhf     = 0.e0_wp
        denah   = 0.30e0_wp
        renah   = 3.566044e0_wp
        benah   = 2.1e0_wp
        rehf    = 1.733e0_wp
        dehf    = 6.122e0_wp
        demf    = 5.452e0_wp
        remf    = 2.96e0_wp
        bemf    = 0.85e0_wp
        c2a     = 1.38661e0_wp
        c2b     = 0.27362e0_wp
        c2c     = 0.15e0_wp
        um2p    = 0.76e0_wp
        g_ehfid = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        cevau  = 1.e0_wp / 27.2113961e0_wp
        onethd = 1.0e0_wp / 3.0e0_wp
        twothd = 2.0e0_wp / 3.0e0_wp
        do i = 1, nt
            cosa = 0.5e0_wp * (r(i,2) ** 2 + r(i,3) ** 2 - r(i,1) ** 2) / (r(i,2) * r(i,3))
            rave = 0.5e0_wp * (r(i,1) + r(i,3))
            ! The diatomic curves . . .
            ! The MH singlet & triplet
            benah = 2.0e0_wp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-benah) * g_r(g_i_,i,1)
            enddo
            d1_w = (-benah) * (r(i,1) - renah)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xenah(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xenah = d2_v
            d2_v  = denah * xenah
            d3_v  = xenah - 2.0e0_wp
            d4_b  = d2_v + d3_v * denah
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d4_b * g_xenah(g_i_)
            enddo
            s1   = d2_v * d3_v
            d2_v = 0.05e0_wp * denah * xenah
            d3_v = xenah + 2.0e0_wp
            d4_b = d2_v + d3_v * (0.05e0_wp * denah)
            do g_i_ = 1, g_p_
                g_t1(g_i_) = d4_b * g_xenah(g_i_)
            enddo
            t1 = d2_v * d3_v
            ! The HF singlet and triplet
            do g_i_ = 1, g_p_
                g_yhf(g_i_) = g_r(g_i_,i,2)
            enddo
            yhf  = r(i,2) - 2.1042e0_wp
            d4_v = yhf * yhf
            d1_p = 2.0e0_wp * yhf
            d5_b = 0.059062e0_wp * d1_p + (-0.025647e0_wp)
            do g_i_ = 1, g_p_
                g_behf(g_i_) = d5_b * g_yhf(g_i_)
            enddo
            behf = 1.1622e0_wp - 0.025647e0_wp * yhf + 0.059062e0_wp * d4_v
            d4_v = r(i,2) - rehf
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-behf) * g_r(g_i_,i,2) + (-d4_v) * g_behf( g_i_)
            enddo
            d1_w = (-behf) * d4_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xehfa(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xehfa = d2_v
            d2_v  = dehf * xehfa
            d3_v  = xehfa - 2.0e0_wp
            d4_b  = d2_v + d3_v * dehf
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d4_b * g_xehfa(g_i_)
            enddo
            s2   = d2_v * d3_v
            d2_v = xehfa * xehfa
            d1_p = 2.0e0_wp * xehfa
            d2_b = 0.5e0_wp * dehf
            d5_b = d2_b * 1.1e0_wp + d2_b * d1_p
            do g_i_ = 1, g_p_
                g_t2(g_i_) = d5_b * g_xehfa(g_i_)
            enddo
            t2 = 0.5e0_wp * dehf * (d2_v + 1.1e0_wp * xehfa)
            ! The Modified HF singlet
            ! The  LiF singlet & triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bemf) * g_r(g_i_,i,3)
            enddo
            d1_w = (-bemf) * (r(i,3) - remf)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xemf(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xemf = d2_v
            d2_v = xemf * xemf
            d1_p = 2.0e0_wp * xemf
            d2_b = 0.5e0_wp * (demf + um2p)
            d5_b = d2_b * 1.3e0_wp + d2_b * d1_p
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d5_b * g_xemf(g_i_)
            enddo
            t3   = 0.5e0_wp * (demf + um2p) * (d2_v + 1.3e0_wp * xemf)
            d2_v = (demf + um2p) * xemf
            d3_v = xemf - 2.0e0_wp
            d4_b = d2_v + d3_v * (demf + um2p)
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d4_b * g_xemf(g_i_)
            enddo
            s3 = d2_v * d3_v
            ! The final LEPS form
            do g_i_ = 1, g_p_
                g_coul1(g_i_) = 0.5e0_wp * g_t1(g_i_) + 0.5e0_wp * g_s1(g_i_)
            enddo
            coul1 = 0.5e0_wp * (s1 + t1)
            do g_i_ = 1, g_p_
                g_coul2(g_i_) = 0.5e0_wp * g_t2(g_i_) + 0.5e0_wp * g_s2(g_i_)
            enddo
            coul2 = 0.5e0_wp * (s2 + t2)
            do g_i_ = 1, g_p_
                g_coul3(g_i_) = 0.5e0_wp * g_t3(g_i_) + 0.5e0_wp * g_s3(g_i_)
            enddo
            coul3 = 0.5e0_wp * (s3 + t3)
            do g_i_ = 1, g_p_
                g_exch1(g_i_) = (-0.5e0_wp) * g_t1(g_i_) + 0.5e0_wp * g_s1(g_i_)
            enddo
            exch1 = 0.5e0_wp * (s1 - t1)
            do g_i_ = 1, g_p_
                g_exch2(g_i_) = (-0.5e0_wp) * g_t2(g_i_) + 0.5e0_wp * g_s2(g_i_)
            enddo
            exch2 = 0.5e0_wp * (s2 - t2)
            do g_i_ = 1, g_p_
                g_exch3(g_i_) = (-0.5e0_wp) * g_t3(g_i_) + 0.5e0_wp * g_s3(g_i_)
            enddo
            exch3 = 0.5e0_wp * (s3 - t3)
            d4_v  = (exch1 - exch2) * (exch1 - exch2)
            d3_p  = 2.0e0_wp * (exch1 - exch2)
            d7_v  = (exch2 - exch3) * (exch2 - exch3)
            d2_p  = 2.0e0_wp * (exch2 - exch3)
            d10_v = (exch3 - exch1) * (exch3 - exch1)
            d1_p  = 2.0e0_wp * (exch3 - exch1)
            d5_b  = d1_p + (-d2_p)
            d6_b  = -d1_p + d3_p
            d10_b = d2_p + (-d3_p)
            do g_i_ = 1, g_p_
                g_w(g_i_) = d5_b * g_exch3(g_i_) + d10_b * g_exch2(g_i_) + d6_b * g_exch1(g_i_)
            enddo
            w = d4_v + d7_v + d10_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-c2c) * g_r(g_i_,i,3) + (-c2c) * g_r(g_i_,i,2) + (-c2c) * g_r(g_i_,i,1) + (-c2b) * &
                    & g_w(g_i_)
            enddo
            d1_w = (-c2b) * w - c2c * (r(i,1) + r(i,2) + r(i,3))
            d2_v = exp(d1_w)
            d1_p = d2_v
            d3_b = c2a * d1_p
            do g_i_ = 1, g_p_
                g_cplg2(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            cplg2 = c2a * d2_v
            d4_b  = cplg2 + cplg2
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_cplg2(g_i_) + g_w(g_i_)
            enddo
            d1_w = w + cplg2 * cplg2
            d7_v = sqrt(d1_w)
            if( d1_w .gt. 0.0e0_wp )then
                d1_p = 1.0e0_wp / (2.0e0_wp *  d7_v)
            else
                d1_p = 0.0e0_wp
            endif
            d7_b = (-(1.0e0_wp / sqrt(2.0e0_wp))) * d1_p
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d7_b * g_d1_w(g_i_) + g_coul3(g_i_) + g_coul2 (g_i_) + g_coul1(g_i_)
            enddo
            e(i) = coul1 + coul2 + coul3 - d7_v / sqrt(2.0e0_wp) + um2p + dehf
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
        enddo
    end subroutine dispatch
end module chempotpy_mxhsl_leps_modelsurface_g_pot22_legacy_state_62de5a21

!> @brief Compatibility interface for g_pot22.
subroutine g_pot22( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_mxhsl_leps_modelsurface_g_pot22_legacy_state_62de5a21, only: dispatch
    use chempotpy_kinds,                                                 only: wp
    implicit none
    integer  :: g_p_
    real(wp) :: r(nt,3)
    real(wp) :: g_r(ldg_r,nt,3)
    integer  :: ldg_r
    real(wp) :: e(nt)
    real(wp) :: g_e(ldg_e,nt)
    integer  :: ldg_e
    integer  :: nt
    call dispatch(g_p_=g_p_, r=r, g_r=g_r, ldg_r=ldg_r, e=e, g_e=g_e, ldg_e=ldg_e, nt=nt)
end subroutine g_pot22
