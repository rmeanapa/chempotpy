!   System:                        YRH
!   Common name:                   YRH
!   Functional form:               Modified Extended LEPS (2x2 diabatic fit)
!   Number of derivatives:         1
!   Number of electronic surfaces: 2
!   Interface:                     3-2V
!   References:   A. W. Jasper, M. D. Hack, and D. G. Truhlar, J. Chem.
!                 Vol. 115, p. 1804, 2001.
!   Notes:  This is a model surface set where the diabatic surfaces are
!           weakly coupled and do not cross.  The parameter SCALE determ
!           the strength of the diabatic coupling.
!   Protocol: POT evaluates the stateless diabatic surfaces.
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Y
!      infinitely far from the RH diatomic and R(RH) set equal to the
!      RH equilibrium diatomic value.
!   Parameters:
!   Coordinates:
!      Internal, Definition: R(1) = R(Y-R)
!                            R(2) = R(R-H)
!                            R(3) = R(H-Y)

!> @brief Implement YRH_LEPS_ModelSurface, a diabatic potential matrix for YRH.
module chempotpy_yrh_leps_modelsurface_dpem_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the YRH_LEPS_ModelSurface diabatic potential matrix for YRH.
    subroutine evaluate( x, igrad, u, ug )
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
        real(wp), intent(out) :: u(nstates,nstates)
        real(wp), intent(out) :: ug(nstates,nstates,natoms,3)
        real(wp)              :: d, g
        integer               :: iatom, idir, nder
        real(wp)              :: nt, r(1,3), r2(3), v(1)
        real(wp)              :: p(nstates)
        real(wp)              :: u11(1), u12(1), u22(1)
        real(wp)              :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
        real(wp)              :: dudr(nstates,nstates,3)
        real(wp)              :: dudx(nstates,nstates,9)
        real(wp)              :: t(nstates,nstates)
        real(wp)              :: tmpmat(nstates,nstates)
        real(wp)              :: tx(9), drdx(3,9)
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
        r2 = r(1,:)*bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(dudr, drdx, dudx)
        do iatom=1, natoms
            do idir=1, 3
                j                  = (iatom-1)*3+idir
                ug(:,:,iatom,idir) = dudx(:,:,j)
            enddo
        enddo
    end subroutine evaluate
end module chempotpy_yrh_leps_modelsurface_dpem_pes

!> @brief Expose YRH_LEPS_ModelSurface through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_yrh_leps_modelsurface_dpem_pes, only: evaluate
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

!> @brief Evaluate EvdRdX for the YRH_LEPS_ModelSurface_DPEM potential.
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

!> @brief Evaluate diagonalize for the YRH_LEPS_ModelSurface_DPEM potential.
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

!> @brief Evaluate pot for the YRH_LEPS_ModelSurface_DPEM potential.
subroutine pot( ro, e, de, nt, nsurf )
    use chempotpy_kinds, only: wp
    implicit none
    integer, parameter :: g_p_  = 3
    integer, parameter :: ldg_r = 3
    integer, parameter :: ldg_e = 3
    real(wp)           :: de(3,nt), deo(3,nt), e(nt), g_e(ldg_e,nt), g_r(ldg_r,nt,3), r(nt,3), ro(nt,3)
    integer            :: i, j, k, nsurf, nt
    !     permute the bonds
    !     want r1 = YR
    !          r2 = RH
    !          r3 = YH
    !     so arr #1 = Y + RH
    !            #2 = R + YH
    !            #3 = H + YR
    !     in code:  r1 = RH
    !               r2 = YH
    !               r3 = YR
    !     so...
    do i=1, nt
        r(i,1) = ro(i,2)
        r(i,2) = ro(i,3)
        r(i,3) = ro(i,1)
    enddo
    !     seed matrix
    do i=1, nt
        do j=1, 3
            do k=1, 3
                if( j .eq. k ) g_r(j,i,k) = 1.0e0_wp
                if( j .ne. k ) g_r(j,i,k) = 0.0e0_wp
            enddo
        enddo
    enddo
    if( nsurf .eq. 1 )then
        call g_pot11(g_p_, r, g_r, ldg_r, e, deo, ldg_e, nt)
    elseif( nsurf .eq. 2 )then
        call g_pot12(g_p_, r, g_r, ldg_r, e, deo, ldg_e, nt)
    else
        call g_pot22(g_p_, r, g_r, ldg_r, e, deo, ldg_e, nt)
    endif
    do i=1, nt
        de(1,i) = deo(3,i)
        de(2,i) = deo(1,i)
        de(3,i) = deo(2,i)
    enddo
end subroutine pot

! U11- lower diabat
!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_yrh_leps_modelsurface_dpem_g_pot11_entry_st_d6ebbc37
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
        real(wp) :: beta1_11, beta2_11, beta3_11, c1, ca_11, cb_11, cc_11, cevau, coul1, coul2, coul3, de1_11, de2_11, &
            & de3_11, e(:), exch1, exch2, exch3, r(:,:), r1, r2, r3, re1_11, re2_11, re3_11, rootw, s1, s2, s3, sp1_11, sp2_11, &
            & sp3_11, t1, t2, t3, w
        integer :: i, nt
        !     data is in bohr, and eV units                                    x
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d3_p, d7_b, d10_b, d7_v, d2_p, d5_b, d1_w, d10_v, d2_v, d6_b
        real(wp)           :: d1_p, d2_b, d3_b, d4_v, g_r1(g_pmax_), g_r(:,:,:), g_r2(g_pmax_), g_r3(g_pmax_), g_d1_w(g_pmax_), g_s1(g_pmax_)
        real(wp)           :: g_s2(g_pmax_), g_s3(g_pmax_), g_t1(g_pmax_), g_t2(g_pmax_), g_t3(g_pmax_), g_coul1(g_pmax_), g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_)
        real(wp)           :: g_exch3(g_pmax_), g_w(g_pmax_), g_c1(g_pmax_), &
            & g_rootw(g_pmax_), g_e(:,:)
        integer :: g_ehfid
        de1_11   = 3.9e0_wp
        de2_11   = 4.3e0_wp
        de3_11   = 0.4e0_wp
        re1_11   = 2.1e0_wp
        re2_11   = 2.1e0_wp
        re3_11   = 2.5e0_wp
        beta1_11 = 1.0e0_wp
        beta2_11 = 1.0e0_wp
        beta3_11 = 1.5e0_wp
        ca_11    = 1.5e0_wp
        cb_11    = 1.0e0_wp
        cc_11    = 0.15e0_wp
        sp1_11   = 0.9e0_wp
        sp2_11   = 0.9e0_wp
        sp3_11   = 0.2e0_wp
        g_ehfid  = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        cevau = 1.e0_wp / 27.2113961e0_wp
        do i = 1, nt
            do g_i_ = 1, g_p_
                g_r1(g_i_) = g_r(g_i_,i,1)
            enddo
            r1 = r(i,1)
            do g_i_ = 1, g_p_
                g_r2(g_i_) = g_r(g_i_,i,2)
            enddo
            r2 = r(i,2)
            do g_i_ = 1, g_p_
                g_r3(g_i_) = g_r(g_i_,i,3)
            enddo
            r3 = r(i,3)
            !   SINGLETS (s1, s2, s3)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta1_11) * g_r1(g_i_)
            enddo
            d1_w = (-beta1_11) * (r1 - re1_11)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.e0_wp - d2_v) * (1.e0_wp - d2_v)
            d1_p = 2.0e0_wp * (1.e0_wp - d2_v)
            d6_b = (-(de1_11 * d1_p)) * d2_p
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d6_b * g_d1_w(g_i_)
            enddo
            s1 = de1_11 * d4_v - de1_11
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta2_11) * g_r2(g_i_)
            enddo
            d1_w = (-beta2_11) * (r2 - re2_11)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.e0_wp - d2_v) * (1.e0_wp - d2_v)
            d1_p = 2.0e0_wp * (1.e0_wp - d2_v)
            d6_b = (-(de2_11 * d1_p)) * d2_p
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d6_b * g_d1_w(g_i_)
            enddo
            s2 = de2_11 * d4_v - de2_11
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta3_11) * g_r3(g_i_)
            enddo
            d1_w = (-beta3_11) * (r3 - re3_11)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.e0_wp - d2_v) * (1.e0_wp - d2_v)
            d1_p = 2.0e0_wp * (1.e0_wp - d2_v)
            d6_b = (-(de3_11 * d1_p)) * d2_p
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d6_b * g_d1_w(g_i_)
            enddo
            s3 = de3_11 * d4_v - de3_11
            !   TRIPLETS (t1, t2, t3)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta1_11) * g_r1(g_i_)
            enddo
            d1_w = (-beta1_11) * (r1 - re1_11)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.e0_wp + d2_v) * (1.e0_wp + d2_v)
            d1_p = 2.0e0_wp * (1.e0_wp + d2_v)
            d7_b = 0.5e0_wp * de1_11 * d1_p * d2_p
            do g_i_ = 1, g_p_
                g_t1(g_i_) = d7_b * g_d1_w(g_i_)
            enddo
            t1 = 0.5e0_wp * (de1_11 * d4_v - de1_11)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta2_11) * g_r2(g_i_)
            enddo
            d1_w = (-beta2_11) * (r2 - re2_11)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.e0_wp + d2_v) * (1.e0_wp + d2_v)
            d1_p = 2.0e0_wp * (1.e0_wp + d2_v)
            d7_b = 0.5e0_wp * de2_11 * d1_p * d2_p
            do g_i_ = 1, g_p_
                g_t2(g_i_) = d7_b * g_d1_w(g_i_)
            enddo
            t2 = 0.5e0_wp * (de2_11 * d4_v - de2_11)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta3_11) * g_r3(g_i_)
            enddo
            d1_w = (-beta3_11) * (r3 - re3_11)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.e0_wp + d2_v) * (1.e0_wp + d2_v)
            d1_p = 2.0e0_wp * (1.e0_wp + d2_v)
            d7_b = 0.5e0_wp * de3_11 * d1_p * d2_p
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d7_b * g_d1_w(g_i_)
            enddo
            t3 = 0.5e0_wp * (de3_11 * d4_v - de3_11)
            do g_i_ = 1, g_p_
                g_t1(g_i_) = sp1_11 * g_t1(g_i_)
            enddo
            t1 = t1 * sp1_11
            do g_i_ = 1, g_p_
                g_t2(g_i_) = sp2_11 * g_t2(g_i_)
            enddo
            t2 = t2 * sp2_11
            do g_i_ = 1, g_p_
                g_t3(g_i_) = sp3_11 * g_t3(g_i_)
            enddo
            t3 = t3 * sp3_11
            !  EXCHANGE AND COULOMB  (exch1, coul1, etc.)
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
                g_d1_w(g_i_) = (-cc_11) * g_r3(g_i_) + (-cc_11) * g_r2(g_i_) + (-cc_11) * g_r1(g_i_) + (-cb_11) * &
                    & g_w(g_i_)
            enddo
            d1_w = (-cb_11) * w - cc_11 * (r1 + r2 + r3)
            d2_v = exp(d1_w)
            d1_p = d2_v
            d3_b = ca_11 * d1_p
            do g_i_ = 1, g_p_
                g_c1(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            c1   = ca_11 * d2_v
            d2_b = 1.0e0_wp / 2.e0_wp
            d5_b = d2_b * c1 + d2_b * c1
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_c1(g_i_) + d2_b * g_w(g_i_)
            enddo
            d1_w = (w + c1 * c1) / 2.e0_wp
            d2_v = sqrt(d1_w)
            if( d1_w .gt. 0.0e0_wp )then
                d1_p = 1.0e0_wp / (2.0e0_wp *  d2_v)
            else
                d1_p = 0.e0_wp
            endif
            do g_i_ = 1, g_p_
                g_rootw(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            rootw = d2_v
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = -g_rootw(g_i_) + g_coul3(g_i_) + g_coul2(g_i_ ) + g_coul1(g_i_)
            enddo
            e(i) = coul1 + coul2 + coul3 - rootw + de2_11
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
            !  end the nt loop here:
        enddo
    end subroutine dispatch
end module chempotpy_yrh_leps_modelsurface_dpem_g_pot11_entry_st_d6ebbc37

!> @brief Compatibility interface for g_pot11.
subroutine g_pot11( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_yrh_leps_modelsurface_dpem_g_pot11_entry_st_d6ebbc37, only: dispatch
    use chempotpy_kinds,                                                only: wp
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

!  U22- upper diabat
!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_yrh_leps_modelsurface_dpem_g_pot22_entry_st_d6ebbc37
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
        real(wp) :: beta1_22, beta2_22, beta3_22, cevau, de1_22, de2_11, de2_22, de3_22, e(:), excite, r(:,:), r1, r2, &
            & r3, re1_22, re2_22, re3_22, s1, t2, t3, wgt
        integer :: i, nt
        !     data is in bohr and eV units
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d6_b, d2_v, d7_b, d4_v, d2_p, d1_p, d1_w, g_r1( &
            & g_pmax_), g_r(:,:,:), g_r2(g_pmax_)
        real(wp) :: g_r3(g_pmax_), g_d1_w(g_pmax_), g_s1(g_pmax_), &
            & g_t2(g_pmax_), g_t3(g_pmax_), g_e(:,:)
        integer :: g_ehfid
        de2_11   = 4.3e0_wp
        excite   = 0.36e0_wp
        wgt      = 0.2e0_wp
        de1_22   = 3.9e0_wp
        de2_22   = 4.3e0_wp
        de3_22   = 0.4e0_wp
        re1_22   = 2.1e0_wp
        re2_22   = 2.1e0_wp
        re3_22   = 2.5e0_wp
        beta1_22 = 1.e0_wp
        beta2_22 = 1.0e0_wp
        beta3_22 = 1.5e0_wp
        g_ehfid  = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        cevau = 1.e0_wp / 27.2113961e0_wp
        do i = 1, nt
            do g_i_ = 1, g_p_
                g_r1(g_i_) = g_r(g_i_,i,1)
            enddo
            r1 = r(i,1)
            do g_i_ = 1, g_p_
                g_r2(g_i_) = g_r(g_i_,i,2)
            enddo
            r2 = r(i,2)
            do g_i_ = 1, g_p_
                g_r3(g_i_) = g_r(g_i_,i,3)
            enddo
            r3 = r(i,3)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta1_22) * g_r1(g_i_)
            enddo
            d1_w = (-beta1_22) * (r1 - re1_22)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.e0_wp - d2_v) * (1.e0_wp - d2_v)
            d1_p = 2.0e0_wp * (1.e0_wp - d2_v)
            d6_b = (-(de1_22 * d1_p)) * d2_p
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d6_b * g_d1_w(g_i_)
            enddo
            s1 = de1_22 * d4_v - de1_22
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta2_22) * g_r2(g_i_)
            enddo
            d1_w = (-beta2_22) * (r2 - re2_22)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.e0_wp + d2_v) * (1.e0_wp + d2_v)
            d1_p = 2.0e0_wp * (1.e0_wp + d2_v)
            d7_b = 0.5e0_wp * de2_22 * d1_p * d2_p
            do g_i_ = 1, g_p_
                g_t2(g_i_) = d7_b * g_d1_w(g_i_)
            enddo
            t2 = 0.5e0_wp * (de2_22 * d4_v - de2_22)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta3_22) * g_r3(g_i_)
            enddo
            d1_w = (-beta3_22) * (r3 - re3_22)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.e0_wp + d2_v) * (1.e0_wp + d2_v)
            d1_p = 2.0e0_wp * (1.e0_wp + d2_v)
            d7_b = 0.5e0_wp * de3_22 * d1_p * d2_p
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d7_b * g_d1_w(g_i_)
            enddo
            t3 = 0.5e0_wp * (de3_22 * d4_v - de3_22)
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = wgt * g_t3(g_i_) + wgt * g_t2(g_i_) + g_s1(g_i_)
            enddo
            e(i) = s1 + (t2 + t3) * wgt + excite + de2_11
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
            !  end the nt loop here:
        enddo
    end subroutine dispatch
end module chempotpy_yrh_leps_modelsurface_dpem_g_pot22_entry_st_d6ebbc37

!> @brief Compatibility interface for g_pot22.
subroutine g_pot22( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_yrh_leps_modelsurface_dpem_g_pot22_entry_st_d6ebbc37, only: dispatch
    use chempotpy_kinds,                                                only: wp
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

!  U12- diabatic coupling
!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_yrh_leps_modelsurface_dpem_g_pot12_entry_st_d6ebbc37
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
        real(wp) :: a1, a2, cevau, csth, e(:), r(:,:), r1, r2, r3, re1_12, re2_12, s1, s2, scale, sith, thrad, thshift, &
            & tmp1, tmp2
        integer            :: i, nt
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d1_p, d1_w, d4_b, d3_b, d2_v, g_r1(g_pmax_), g_r(:,:,:), g_r2(g_pmax_), g_s1(g_pmax_), g_s2(g_pmax_)
        real(wp)           :: g_d1_w(g_pmax_), g_tmp1(g_pmax_), g_tmp2(g_pmax_), g_e(:,:)
        integer            :: g_ehfid
        a1      = 0.8e0_wp
        a2      = 1.0e0_wp
        re1_12  = 2.2e0_wp
        re2_12  = 4.2e0_wp
        thshift = 45.e0_wp
        scale   = 0.2e0_wp
        g_ehfid = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        cevau = 1.e0_wp / 27.2113961e0_wp
        do i = 1, nt
            do g_i_ = 1, g_p_
                g_r1(g_i_) = g_r(g_i_,i,1)
            enddo
            r1 = r(i,1)
            do g_i_ = 1, g_p_
                g_r2(g_i_) = g_r(g_i_,i,2)
            enddo
            r2    = r(i,2)
            r3    = r(i,3)
            thrad = acos(-1.e0_wp) / 180.e0_wp * thshift
            csth  = cos(thrad)
            sith  = sin(thrad)
            do g_i_ = 1, g_p_
                g_s1(g_i_) = sith * g_r2(g_i_) + csth * g_r1(g_i_)
            enddo
            s1 = (r1 - re1_12) * csth + (r2 - re2_12) * sith
            do g_i_ = 1, g_p_
                g_s2(g_i_) = (-csth) * g_r2(g_i_) + sith * g_r1(g_i_)
            enddo
            s2   = (r1 - re1_12) * sith - (r2 - re2_12) * csth
            d2_v = s1 * s1
            d1_p = 2.0e0_wp * s1
            d3_b = (-a1) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_s1(g_i_)
            enddo
            d1_w = (-a1) * d2_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_tmp1(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            tmp1 = d2_v
            d2_v = s2 * s2
            d1_p = 2.0e0_wp * s2
            d3_b = (-a2) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_s2(g_i_)
            enddo
            d1_w = (-a2) * d2_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_tmp2(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            tmp2 = d2_v
            d3_b = scale * tmp2
            d4_b = scale * tmp1
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d4_b * g_tmp2(g_i_) + d3_b * g_tmp1(g_i_)
            enddo
            e(i) = tmp1 * tmp2 * scale
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
            !   end the nt loop here:
        enddo
    end subroutine dispatch
end module chempotpy_yrh_leps_modelsurface_dpem_g_pot12_entry_st_d6ebbc37

!> @brief Compatibility interface for g_pot12.
subroutine g_pot12( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_yrh_leps_modelsurface_dpem_g_pot12_entry_st_d6ebbc37, only: dispatch
    use chempotpy_kinds,                                                only: wp
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
