!   System:          BrH2
!   Functional form: extended-LEPS plus three-center term with
!                    internal-angle-dependent Sato parameters
!   Common name: BrH2SECD
!   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, et. al.
!                    J. Phys. Chem. 99, 207-225 (1995)
!   Analytic derivatives added, Yuri Volobuev, July-Aug 1998
!   The potential parameters are included in DATA statements.
!   For the this potential energy surface the potential energy in the
!   BC valley is equal to the energy of the HH diatomic.
!   The energy in the AC valley is equal to the energy of the HBr diatom
!   which is also the energy value in the BA valley.
!   This potential is written such that:
!                  R(1) = R(Br-H)
!                  R(2) = R(H-H)
!                  R(3) = R(H-Br)
!   The zero of energy is defined at Br "infinitely" far from the H2 dia

!> @brief Implement BrH2_LEPS_LTBZ2_1995, a diabatic potential matrix for BrH2.
module chempotpy_brh2_leps_ltbz2_1995_dpem_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the BrH2_LEPS_LTBZ2_1995 diabatic potential matrix for BrH2.
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
        real(wp)              :: d, g, p
        integer               :: iatom, idir, nder
        real(wp)              :: nt, r(1,3), r2(3), v(1)
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
        r2          = r(1,:)*bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(dudr, drdx, dudx)
        do iatom=1, natoms
            do idir=1, 3
                j                  = (iatom-1)*3+idir
                ug(:,:,iatom,idir) = dudx(:,:,j)
            enddo
        enddo
    end subroutine evaluate
end module chempotpy_brh2_leps_ltbz2_1995_dpem_pes

!> @brief Expose BrH2_LEPS_LTBZ2_1995 through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_brh2_leps_ltbz2_1995_dpem_pes, only: evaluate
    use chempotpy_kinds,                         only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: u(nstates,nstates)
    real(wp), intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem

!> @brief Evaluate EvdRdX for the BrH2_LEPS_LTBZ2_1995_DPEM potential.
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

!> @brief Evaluate diagonalize for the BrH2_LEPS_LTBZ2_1995_DPEM potential.
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

!> @brief Evaluate pot for the BrH2_LEPS_LTBZ2_1995_DPEM potential.
subroutine pot( r, e, de, nt, nsurf )
    use chempotpy_kinds, only: wp
    implicit none
    integer, parameter :: g_p_  = 3
    integer, parameter :: ldg_e = 3
    real(wp)           :: de(3,nt), e(nt), g_e(ldg_e,nt), r(nt,3)
    integer            :: i, j, k, ldg_r, nsurf, nt
    if( nsurf .eq. 1 )then
        call g_pot11(g_p_, r, e, de, ldg_e, nt)
    elseif( nsurf .eq. 2 )then
        call pot12(r, e, de, nt)
    elseif( nsurf .eq. 3 )then
        call pot22(r, e, de, nt)
    endif
end subroutine pot

!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_brh2_leps_ltbz2_1995_dpem_g_pot11_entry_sta_432f3192
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, rr, e, g_e, ldg_e, n )
        !   this is a vectorized version of the BrH2SEC potential with derivativ
        use chempotpy_kinds, only: wp
        implicit none
        optional :: g_p_, rr, e, g_e, ldg_e, n
        real(wp) :: a1, a2, a3, a4, a5, a6, b(3), beta(3), betaa(3), capr, capr2, coschi, costh, coul(3), d(3), de(3), &
            & deq(3), derad(3), dist, do4z(3), dv3c(3), dxdr(3, 3), dz(3,3), e(:), exch(3), fcn4, op3z(3), pe, r132p2, r1m3, &
            & r1m3p2, r1p2, r1p3, r1p3m2, r2, r2i, r2p2, r3p2, rad, re(3), rr(:,:), schip2, sintp2, sintp6, sintp8, top3z(3), &
            & tzp3(3), v3c, v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, x(3), xlamb, z(3), zchi, zhbrp1, zhbrp2, zhhp1, zhhp2, &
            & zp3(3), zpo(3)
        integer             :: i, k, n
        real(wp)            :: r(3)
        real(wp), parameter :: cangau  = 0.52917706e0_wp
        real(wp), parameter :: charkc  = 627.5095e0_wp
        integer, parameter  :: g_pmax_ = 3
        integer             :: g_i_, g_p_, ldg_rr, ldg_e
        real(wp)            :: d3_p, d2_p, d2_v, d10_v, d2_b, d3_v, d4_v, &
            & d5_v, d6_v, d7_v, d3_b, d4_b, d5_b, d6_b, d7_b, d1_p, d8_v, &
            & d8_b, d1_w, d10_b, g_r(g_pmax_, 3), g_r1p2(g_pmax_), &
            & g_r2p2(g_pmax_), g_r3p2(g_pmax_), g_zchi(g_pmax_), &
            & g_capr2(g_pmax_), g_capr(g_pmax_), g_coschi(g_pmax_), &
            & g_schip2(g_pmax_), g_z(g_pmax_, 3), g_zpo(g_pmax_, 3), &
            & g_op3z(g_pmax_, 3), g_top3z(g_pmax_, 3), g_zp3(g_pmax_, 3), &
            & g_tzp3(g_pmax_, 3), g_do4z(g_pmax_, 3), g_pe(g_pmax_), &
            & g_d1_w(g_pmax_), g_x(g_pmax_, 3), g_coul(g_pmax_, 3), &
            & g_exch(g_pmax_, 3), g_rad(g_pmax_), g_r1p3(g_pmax_), &
            & g_r1m3(g_pmax_), g_r1p3m2(g_pmax_), g_r1m3p2(g_pmax_), &
            & g_r132p2(g_pmax_), g_dist(g_pmax_), g_costh(g_pmax_), &
            & g_sintp2(g_pmax_), g_sintp6(g_pmax_), g_sintp8(g_pmax_), &
            & g_fcn4(g_pmax_), g_v3c(g_pmax_), g_e(:,:), g_zhhp1(g_pmax_), &
            & g_zhhp2(g_pmax_), g_zhbrp1(g_pmax_), g_zhbrp2(g_pmax_), &
            & g_dz(g_pmax_, 3, 3), g_beta(g_pmax_, 3), g_v3cp1(g_pmax_), &
            & g_v3cp2(g_pmax_), g_v3cp3(g_pmax_), g_v3cp4(g_pmax_), &
            & g_v3cp5(g_pmax_), g_v3cp6(g_pmax_), g_dv3c(g_pmax_, 3)
        integer :: g_myrank, g_nprocs
        a1       = 0.e0_wp
        a2       = 0.e0_wp
        a3       = 0.e0_wp
        a4       = 0.e0_wp
        a5       = 0.e0_wp
        a6       = 0.e0_wp
        b        = 0.e0_wp
        beta     = 0.e0_wp
        betaa    = 0.e0_wp
        capr     = 0.e0_wp
        capr2    = 0.e0_wp
        coschi   = 0.e0_wp
        costh    = 0.e0_wp
        coul     = 0.e0_wp
        d        = 0.e0_wp
        d10_b    = 0.e0_wp
        d10_v    = 0.e0_wp
        d1_p     = 0.e0_wp
        d1_w     = 0.e0_wp
        d2_b     = 0.e0_wp
        d2_p     = 0.e0_wp
        d2_v     = 0.e0_wp
        d3_b     = 0.e0_wp
        d3_p     = 0.e0_wp
        d3_v     = 0.e0_wp
        d4_b     = 0.e0_wp
        d4_v     = 0.e0_wp
        d5_b     = 0.e0_wp
        d5_v     = 0.e0_wp
        d6_b     = 0.e0_wp
        d6_v     = 0.e0_wp
        d7_b     = 0.e0_wp
        d7_v     = 0.e0_wp
        d8_b     = 0.e0_wp
        d8_v     = 0.e0_wp
        de       = 0.e0_wp
        deq      = 0.e0_wp
        derad    = 0.e0_wp
        dist     = 0.e0_wp
        do4z     = 0.e0_wp
        dv3c     = 0.e0_wp
        dxdr     = 0.e0_wp
        dz       = 0.e0_wp
        exch     = 0.e0_wp
        fcn4     = 0.e0_wp
        g_beta   = 0.e0_wp
        g_capr   = 0.e0_wp
        g_capr2  = 0.e0_wp
        g_coschi = 0.e0_wp
        g_costh  = 0.e0_wp
        g_coul   = 0.e0_wp
        g_d1_w   = 0.e0_wp
        g_dist   = 0.e0_wp
        g_do4z   = 0.e0_wp
        g_dv3c   = 0.e0_wp
        g_dz     = 0.e0_wp
        g_exch   = 0.e0_wp
        g_fcn4   = 0.e0_wp
        g_op3z   = 0.e0_wp
        g_pe     = 0.e0_wp
        g_r      = 0.e0_wp
        g_r132p2 = 0.e0_wp
        g_r1m3   = 0.e0_wp
        g_r1m3p2 = 0.e0_wp
        g_r1p2   = 0.e0_wp
        g_r1p3   = 0.e0_wp
        g_r1p3m2 = 0.e0_wp
        g_r2p2   = 0.e0_wp
        g_r3p2   = 0.e0_wp
        g_rad    = 0.e0_wp
        g_schip2 = 0.e0_wp
        g_sintp2 = 0.e0_wp
        g_sintp6 = 0.e0_wp
        g_sintp8 = 0.e0_wp
        g_top3z  = 0.e0_wp
        g_tzp3   = 0.e0_wp
        g_v3c    = 0.e0_wp
        g_v3cp1  = 0.e0_wp
        g_v3cp2  = 0.e0_wp
        g_v3cp3  = 0.e0_wp
        g_v3cp4  = 0.e0_wp
        g_v3cp5  = 0.e0_wp
        g_v3cp6  = 0.e0_wp
        g_x      = 0.e0_wp
        g_z      = 0.e0_wp
        g_zchi   = 0.e0_wp
        g_zhbrp1 = 0.e0_wp
        g_zhbrp2 = 0.e0_wp
        g_zhhp1  = 0.e0_wp
        g_zhhp2  = 0.e0_wp
        g_zp3    = 0.e0_wp
        g_zpo    = 0.e0_wp
        op3z     = 0.e0_wp
        pe       = 0.e0_wp
        r        = 0.e0_wp
        r132p2   = 0.e0_wp
        r1m3     = 0.e0_wp
        r1m3p2   = 0.e0_wp
        r1p2     = 0.e0_wp
        r1p3     = 0.e0_wp
        r1p3m2   = 0.e0_wp
        r2       = 0.e0_wp
        r2i      = 0.e0_wp
        r2p2     = 0.e0_wp
        r3p2     = 0.e0_wp
        rad      = 0.e0_wp
        re       = 0.e0_wp
        schip2   = 0.e0_wp
        sintp2   = 0.e0_wp
        sintp6   = 0.e0_wp
        sintp8   = 0.e0_wp
        top3z    = 0.e0_wp
        tzp3     = 0.e0_wp
        v3c      = 0.e0_wp
        v3cp1    = 0.e0_wp
        v3cp2    = 0.e0_wp
        v3cp3    = 0.e0_wp
        v3cp4    = 0.e0_wp
        v3cp5    = 0.e0_wp
        v3cp6    = 0.e0_wp
        x        = 0.e0_wp
        xlamb    = 0.e0_wp
        z        = 0.e0_wp
        zchi     = 0.e0_wp
        zhbrp1   = 0.e0_wp
        zhbrp2   = 0.e0_wp
        zhhp1    = 0.e0_wp
        zhhp2    = 0.e0_wp
        zp3      = 0.e0_wp
        zpo      = 0.e0_wp
        d        = [ 90.4473e0_wp, 109.559e0_wp, 90.4473e0_wp]
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        xlamb  = .4568e0_wp / 27.21161e0_wp
        a1     = v3cp1
        a2     = v3cp2
        a3     = v3cp3
        a4     = v3cp4
        a5     = v3cp5
        a6     = v3cp6
        re     = [ 1.41443e0_wp, 0.74144e0_wp, 1.41443e0_wp]
        betaa  = [ 1.81088e0_wp, 1.94198e0_wp, 1.81088e0_wp]
        r2     = 1.41421356e0_wp
        zhhp1  = 2.91230499e-1_wp
        zhhp2  = 4.99160974e0_wp
        zhbrp1 = 3.84262775e-3_wp
        zhbrp2 = 2.87642885e-3_wp
        v3cp1  = 1.67457884e0_wp
        v3cp2  = 7.09518875e-1_wp
        v3cp3  = 1.87582679e-1_wp
        v3cp4  = 2.43517074e-1_wp
        v3cp5  = 2.04145206e-1_wp
        v3cp6  = 2.57833022e0_wp
        !   Convert the constants to atomic units
        do i = 1, 3
            de(i)   = d(i) / charkc
            re(i)   = re(i) / cangau
            beta(i) = betaa(i) * cangau
        enddo
        do k = 1, n
            do g_i_ = 1, g_p_
                g_r(g_i_, 1) = 0.0e0_wp
            enddo
            r(1) = rr(k,1)
            do g_i_ = 1, g_p_
                g_r(g_i_, 2) = 0.0e0_wp
            enddo
            r(2) = rr(k,2)
            do g_i_ = 1, g_p_
                g_r(g_i_, 3) = 0.0e0_wp
            enddo
            r(3)      = rr(k,3)
            g_r(1, 1) = 1.0e0_wp
            g_r(2, 2) = 1.0e0_wp
            g_r(3, 3) = 1.0e0_wp
            !   Compute bond dependent Sato parameters
            d2_b = r(1) + r(1)
            do g_i_ = 1, g_p_
                g_r1p2(g_i_) = d2_b * g_r(g_i_, 1)
            enddo
            r1p2 = r(1) * r(1)
            d2_b = r(2) + r(2)
            do g_i_ = 1, g_p_
                g_r2p2(g_i_) = d2_b * g_r(g_i_, 2)
            enddo
            r2p2 = r(2) * r(2)
            d2_b = r(3) + r(3)
            do g_i_ = 1, g_p_
                g_r3p2(g_i_) = d2_b * g_r(g_i_, 3)
            enddo
            r3p2 = r(3) * r(3)
            r2i  = 1.0e0_wp / r(2)
            do g_i_ = 1, g_p_
                g_zchi(g_i_) = 0.0e0_wp
            enddo
            zchi = 0.0e0_wp
            !   CAPR IS THE Br TO H2 DISTANCE
            d5_b = (-0.5e0_wp) * 0.5e0_wp
            do g_i_ = 1, g_p_
                g_capr2(g_i_) = d5_b * g_r2p2(g_i_) + 0.5e0_wp * g_r3p2(g_i_) + 0.5e0_wp * g_r1p2(g_i_)
            enddo
            capr2 = 0.5e0_wp * (r1p2 + r3p2 - 0.5e0_wp * r2p2)
            if( capr2 .le. 0.e0_wp )then
                do g_i_ = 1, g_p_
                    g_capr2(g_i_) = 0.0e0_wp
                enddo
                capr2 = 0.0e0_wp
                do g_i_ = 1, g_p_
                    g_capr(g_i_) = 0.0e0_wp
                enddo
                capr = 0.0e0_wp
            else
                d2_v = sqrt(capr2)
                if( capr2 .gt. 0.0e0_wp )then
                    d1_p = 1.0e0_wp / (2.0e0_wp *  d2_v)
                else
                    d1_p = 0.0e0_wp
                endif
                do g_i_ = 1, g_p_
                    g_capr(g_i_) = d1_p * g_capr2(g_i_)
                enddo
                capr = d2_v
                !   Compute Chi, the angle between CAPR and RHH
                d7_v = r(2) * capr
                d8_v = 0.5e0_wp * (r3p2 - r1p2) / d7_v
                d3_b = (-d8_v) / d7_v
                d4_b = d3_b * capr
                d5_b = d3_b * r(2)
                d6_b = 1.0e0_wp / d7_v * 0.5e0_wp
                do g_i_ = 1, g_p_
                    g_coschi(g_i_) = d5_b * g_capr(g_i_) + d4_b * g_r(g_i_, 2) + (-d6_b) * g_r1p2(g_i_) + d6_b * &
                        & g_r3p2(g_i_)
                enddo
                coschi = d8_v
                d3_b   = -coschi + (-coschi)
                do g_i_ = 1, g_p_
                    g_schip2(g_i_) = d3_b * g_coschi(g_i_)
                enddo
                schip2 = 1.0e0_wp - coschi * coschi
                d4_v   = zhbrp2 * schip2
                d6_v   = zhbrp1 * schip2 + d4_v * schip2
                d7_b   = capr2 * d4_v + capr2 * schip2 * zhbrp2 + capr2 * zhbrp1
                do g_i_ = 1, g_p_
                    g_zchi(g_i_) = d7_b * g_schip2(g_i_) + d6_v * g_capr2(g_i_ )
                enddo
                zchi = capr2 * d6_v
            endif
            do g_i_ = 1, g_p_
                g_z(g_i_, 1) = g_zchi(g_i_)
            enddo
            z(1) = 0.1675e0_wp + zchi
            do g_i_ = 1, g_p_
                g_z(g_i_, 3) = g_z(g_i_, 1)
            enddo
            z(3) = z(1)
            d3_v = zhhp2 + r2p2
            d4_v = zhhp1 * r2p2 / d3_v
            d4_b = (-d4_v) / d3_v + 1.0e0_wp / d3_v * zhhp1
            do g_i_ = 1, g_p_
                g_z(g_i_, 2) = d4_b * g_r2p2(g_i_)
            enddo
            z(2) = d4_v
            !   Compute the LEPS parameters
            do i = 1, 3
                do g_i_ = 1, g_p_
                    g_zpo(g_i_, i) = g_z(g_i_, i)
                enddo
                zpo(i) = 1.0e0_wp + z(i)
                do g_i_ = 1, g_p_
                    g_op3z(g_i_, i) = 3.0e0_wp * g_z(g_i_, i)
                enddo
                op3z(i) = 1.0e0_wp + 3.0e0_wp * z(i)
                do g_i_ = 1, g_p_
                    g_top3z(g_i_, i) = 2.0e0_wp * g_op3z(g_i_, i)
                enddo
                top3z(i) = 2.0e0_wp * op3z(i)
                do g_i_ = 1, g_p_
                    g_zp3(g_i_, i) = g_z(g_i_, i)
                enddo
                zp3(i) = z(i) + 3.0e0_wp
                do g_i_ = 1, g_p_
                    g_tzp3(g_i_, i) = 2.0e0_wp * g_zp3(g_i_, i)
                enddo
                tzp3(i) = 2.0e0_wp * zp3(i)
                d2_v    = de(i) / 4.0e0_wp / zpo(i)
                d2_b    = (-d2_v) / zpo(i)
                do g_i_ = 1, g_p_
                    g_do4z(g_i_, i) = d2_b * g_zpo(g_i_, i)
                enddo
                do4z(i) = d2_v
                b(i)    = beta(i) * do4z(i) * 2.0e0_wp
            enddo
            do g_i_ = 1, g_p_
                g_pe(g_i_) = 0.0e0_wp
            enddo
            pe = 0.e0_wp
            do i = 1, 3
                do g_i_ = 1, g_p_
                    g_d1_w(g_i_) = (-beta(i)) * g_r(g_i_, i)
                enddo
                d1_w = (-beta(i)) * (r(i) - re(i))
                d2_v = exp(d1_w)
                d1_p = d2_v
                do g_i_ = 1, g_p_
                    g_x(g_i_, i) = d1_p * g_d1_w(g_i_)
                enddo
                x(i) = d2_v
                d6_v = zp3(i) * x(i) - top3z(i)
                d7_v = do4z(i) * d6_v
                d4_b = x(i) * d6_v
                d5_b = x(i) * do4z(i)
                d8_b = d5_b * x(i)
                d3_b = d7_v + d5_b * zp3(i)
                do g_i_ = 1, g_p_
                    g_coul(g_i_, i) = (-d5_b) * g_top3z(g_i_, i) + d3_b * g_x( g_i_, i) + d8_b * g_zp3(g_i_, i) + &
                        & d4_b * g_do4z(g_i_, i)
                enddo
                coul(i) = d7_v * x(i)
                d6_v    = op3z(i) * x(i) - tzp3(i)
                d7_v    = do4z(i) * d6_v
                d4_b    = x(i) * d6_v
                d5_b    = x(i) * do4z(i)
                d8_b    = d5_b * x(i)
                d3_b    = d7_v + d5_b * op3z(i)
                do g_i_ = 1, g_p_
                    g_exch(g_i_, i) = (-d5_b) * g_tzp3(g_i_, i) + d3_b * g_x(g_i_, i) + d8_b * g_op3z(g_i_, i) + &
                        & d4_b * g_do4z(g_i_, i)
                enddo
                exch(i) = d7_v * x(i)
                do g_i_ = 1, g_p_
                    g_pe(g_i_) = g_coul(g_i_, i) + g_pe(g_i_)
                enddo
                pe = pe + coul(i)
            enddo
            d4_v  = (exch(1) - exch(2)) * (exch(1) - exch(2))
            d3_p  = 2.0e0_wp * (exch(1) - exch(2))
            d7_v  = (exch(2) - exch(3)) * (exch(2) - exch(3))
            d2_p  = 2.0e0_wp * (exch(2) - exch(3))
            d10_v = (exch(3) - exch(1)) * (exch(3) - exch(1))
            d1_p  = 2.0e0_wp * (exch(3) - exch(1))
            d5_b  = d1_p + (-d2_p)
            d6_b  = -d1_p + d3_p
            d10_b = d2_p + (-d3_p)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_exch(g_i_, 3) + d10_b * g_exch(g_i_, 2) + d6_b * g_exch(g_i_, 1)
            enddo
            d1_w = d4_v + d7_v + d10_v
            d2_v = sqrt(d1_w)
            if( d1_w .gt. 0.0e0_wp )then
                d1_p = 1.0e0_wp / (2.0e0_wp *  d2_v)
            else
                d1_p = 0.0e0_wp
            endif
            do g_i_ = 1, g_p_
                g_rad(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            rad  = d2_v
            d5_b = -(1.0e0_wp / r2)
            do g_i_ = 1, g_p_
                g_pe(g_i_) = d5_b * g_rad(g_i_) + g_pe(g_i_)
            enddo
            pe = pe - rad / r2 + de(2)
            !   Compute the three-center term
            do g_i_ = 1, g_p_
                g_r1p3(g_i_) = g_r(g_i_, 3) + g_r(g_i_, 1)
            enddo
            r1p3 = r(1) + r(3)
            a1   = v3cp1
            a2   = v3cp2
            a3   = v3cp3
            a4   = v3cp4
            a5   = v3cp5
            a6   = v3cp6
            do g_i_ = 1, g_p_
                g_r1m3(g_i_) = -g_r(g_i_, 3) + g_r(g_i_, 1)
            enddo
            r1m3 = r(1) - r(3)
            do g_i_ = 1, g_p_
                g_r1p3m2(g_i_) = -g_r(g_i_, 2) + g_r1p3(g_i_)
            enddo
            r1p3m2 = r1p3 - r(2)
            d2_b   = r1m3 + r1m3
            do g_i_ = 1, g_p_
                g_r1m3p2(g_i_) = d2_b * g_r1m3(g_i_)
            enddo
            r1m3p2 = r1m3 * r1m3
            d2_b   = r1p3m2 + r1p3m2
            do g_i_ = 1, g_p_
                g_r132p2(g_i_) = d2_b * g_r1p3m2(g_i_)
            enddo
            r132p2 = r1p3m2 * r1p3m2
            !   Compute the bond angle from the three internal coordinates
            d4_b = -r(2) + (-r(2))
            d7_b = r(3) + r(3)
            d8_b = r(1) + r(1)
            do g_i_ = 1, g_p_
                g_dist(g_i_) = d4_b * g_r(g_i_, 2) + d7_b * g_r(g_i_, 3) + d8_b * g_r(g_i_, 1)
            enddo
            dist = r(1) * r(1) + r(3) * r(3) - r(2) * r(2)
            d3_v = 2.0e0_wp * r(1)
            d5_v = d3_v * r(3)
            d6_v = dist / d5_v
            d2_b = 1.0e0_wp / d5_v
            d3_b = (-d6_v) / d5_v
            d5_b = d3_b * d3_v
            d6_b = d3_b * r(3) * 2.0e0_wp
            do g_i_ = 1, g_p_
                g_costh(g_i_) = d5_b * g_r(g_i_, 3) + d6_b * g_r(g_i_, 1) + d2_b * g_dist(g_i_)
            enddo
            costh = d6_v
            d3_b  = -costh + (-costh)
            do g_i_ = 1, g_p_
                g_sintp2(g_i_) = d3_b * g_costh(g_i_)
            enddo
            sintp2 = 1.0e0_wp - costh * costh
            d2_v   = sintp2 * sintp2
            d3_b   = d2_v + sintp2 * sintp2 + sintp2 * sintp2
            do g_i_ = 1, g_p_
                g_sintp6(g_i_) = d3_b * g_sintp2(g_i_)
            enddo
            sintp6 = d2_v * sintp2
            do g_i_ = 1, g_p_
                g_sintp8(g_i_) = sintp6 * g_sintp2(g_i_) + sintp2 * g_sintp6 (g_i_)
            enddo
            sintp8 = sintp6 * sintp2
            do g_i_ = 1, g_p_
                g_fcn4(g_i_) = a6 * g_sintp8(g_i_) + a5 * g_sintp2(g_i_)
            enddo
            fcn4 = 1.e0_wp + a5 * sintp2 + a6 * sintp8
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-a4) * g_r132p2(g_i_) + (-a3) * g_r1m3p2(g_i_) + (-a2) * g_r1p3(g_i_)
            enddo
            d1_w = -(a2 * r1p3 + a3 * r1m3p2 + a4 * r132p2)
            d2_v = a1 * fcn4
            d4_v = exp(d1_w)
            d1_p = d4_v
            d4_b = d2_v * d1_p
            d5_b = d4_v * a1
            do g_i_ = 1, g_p_
                g_v3c(g_i_) = d4_b * g_d1_w(g_i_) + d5_b * g_fcn4(g_i_)
            enddo
            v3c = d2_v * d4_v
            !   Sum all the terms that make up the energy
            do g_i_ = 1, g_p_
                g_pe(g_i_) = g_v3c(g_i_) + g_pe(g_i_)
            enddo
            pe = pe + v3c
            do g_i_ = 1, g_p_
                g_e(g_i_,k) = g_pe(g_i_)
            enddo
            e(k) = pe + xlamb / 3.e0_wp
        enddo
    end subroutine dispatch
end module chempotpy_brh2_leps_ltbz2_1995_dpem_g_pot11_entry_sta_432f3192

!> @brief Compatibility interface for g_pot11.
subroutine g_pot11( g_p_, rr, e, g_e, ldg_e, n )
    use chempotpy_brh2_leps_ltbz2_1995_dpem_g_pot11_entry_sta_432f3192, only: dispatch
    use chempotpy_kinds,                                                only: wp
    implicit none
    integer  :: g_p_
    real(wp) :: rr(n,3)
    real(wp) :: e(n)
    real(wp) :: g_e(ldg_e,*)
    integer  :: ldg_e
    integer  :: n
    call dispatch(g_p_=g_p_, rr=rr, e=e, g_e=g_e(1:ldg_e,1:n), ldg_e=ldg_e, n=n)
end subroutine g_pot11

!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_brh2_leps_ltbz2_1995_dpem_pot12_legacy_state_432f3192
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( rr, e, de, n )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: rr, e, de, n
        real(wp) :: de(:,:), e(:), rr(:,:), xlamb
        integer  :: j, k, n
        xlamb = 0.e0_wp
        xlamb = 0.4568e0_wp/27.21161e0_wp
        do k=1, n
            e(k) = -sqrt(2.e0_wp)*xlamb/3.e0_wp
            do j =1, 3
                de(j,k) = 0.0e0_wp
            enddo
        enddo
    end subroutine dispatch
end module chempotpy_brh2_leps_ltbz2_1995_dpem_pot12_legacy_state_432f3192

!> @brief Compatibility interface for pot12.
subroutine pot12( rr, e, de, n )
    use chempotpy_brh2_leps_ltbz2_1995_dpem_pot12_legacy_state_432f3192, only: dispatch
    use chempotpy_kinds,                                                 only: wp
    implicit none
    real(wp) :: rr(n,3)
    real(wp) :: e(*)
    real(wp) :: de(3,n)
    integer  :: n
    call dispatch(rr=rr, e=e(1:n), de=de, n=n)
end subroutine pot12

!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_brh2_leps_ltbz2_1995_dpem_pot22_legacy_state_432f3192
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( rr, e, der, n )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: rr, e, der, n
        real(wp) :: b(3), beta(3), betaa(3), ccon, coul(3), d(3), de(3), der(:,:), do4z(3), dz(3,3), e(:), exch(3), &
            & op3z(3), pe, r2, rad, re(3), rr(:,:), top3z(3), tzp3(3), vanti1, vanti3, vmorse, x(3), xlamb, y1, y2, y3, z(3), &
            & zhbrp1, zhbrp2, zhhp1, zhhp2, zp3(3), zpo(3)
        integer             :: i, k, n
        real(wp)            :: r(3)
        real(wp), parameter :: cangau = 0.52917706e0_wp
        real(wp), parameter :: charkc = 627.5095e0_wp
        integer             :: g_myrank, g_nprocs
        b      = 0.e0_wp
        beta   = 0.e0_wp
        betaa  = 0.e0_wp
        ccon   = 0.e0_wp
        coul   = 0.e0_wp
        d      = 0.e0_wp
        de     = 0.e0_wp
        do4z   = 0.e0_wp
        dz     = 0.e0_wp
        exch   = 0.e0_wp
        op3z   = 0.e0_wp
        pe     = 0.e0_wp
        r      = 0.e0_wp
        r2     = 0.e0_wp
        rad    = 0.e0_wp
        re     = 0.e0_wp
        top3z  = 0.e0_wp
        tzp3   = 0.e0_wp
        vanti1 = 0.e0_wp
        vanti3 = 0.e0_wp
        vmorse = 0.e0_wp
        x      = 0.e0_wp
        xlamb  = 0.e0_wp
        y1     = 0.e0_wp
        y2     = 0.e0_wp
        y3     = 0.e0_wp
        z      = 0.e0_wp
        zhbrp1 = 0.e0_wp
        zhbrp2 = 0.e0_wp
        zhhp1  = 0.e0_wp
        zhhp2  = 0.e0_wp
        zp3    = 0.e0_wp
        zpo    = 0.e0_wp
        d      = [ 90.4473e0_wp, 109.559e0_wp, 90.4473e0_wp]
        re     = [ 1.41443e0_wp, 0.74144e0_wp, 1.41443e0_wp]
        betaa  = [ 1.81088e0_wp, 1.94198e0_wp, 1.81088e0_wp]
        r2     = 1.41421356e0_wp
        ccon   = 0.5e0_wp
        xlamb  = 0.4568e0_wp/27.21161e0_wp
        zhhp1  = 2.91230499e-1_wp
        zhhp2  = 4.99160974e0_wp
        zhbrp1 = 3.84262775e-3_wp
        zhbrp2 = 2.87642885e-3_wp
        !   Convert the constants to atomic units
        do i = 1, 3
            de(i)   = d(i) / charkc
            re(i)   = re(i) / cangau
            beta(i) = betaa(i) * cangau
        enddo
        do k=1, n
            r(1)     = rr(k,1)
            r(2)     = rr(k,2)
            r(3)     = rr(k,3)
            y2       = exp(-beta(2)*(r(2)-re(2)))
            y1       = exp(-beta(1)*(r(1)-re(1)))
            y3       = exp(-beta(3)*(r(3)-re(3)))
            vmorse   = de(2)*(y2**2-2.e0_wp*y2)
            vanti1   = de(1)*(y1**2+2.e0_wp*y1)
            vanti3   = de(3)*(y3**2+2.e0_wp*y3)
            pe       = de(2)+vmorse+(vanti1+vanti3)*ccon+2.e0_wp*xlamb/3.e0_wp
            der(1,k) = -2.0e0_wp*beta(1)*de(1)*(y1**2+y1)*ccon
            der(2,k) = -2.0e0_wp*beta(2)*de(2)*(y2**2-y2)
            der(3,k) = -2.0e0_wp*beta(3)*de(3)*(y3**2+y3)*ccon
            e(k)     = pe
        enddo
    end subroutine dispatch
end module chempotpy_brh2_leps_ltbz2_1995_dpem_pot22_legacy_state_432f3192

!> @brief Compatibility interface for pot22.
subroutine pot22( rr, e, der, n )
    use chempotpy_brh2_leps_ltbz2_1995_dpem_pot22_legacy_state_432f3192, only: dispatch
    use chempotpy_kinds,                                                 only: wp
    implicit none
    real(wp) :: rr(n,3)
    real(wp) :: e(*)
    real(wp) :: der(3,n)
    integer  :: n
    call dispatch(rr=rr, e=e(1:n), der=der, n=n)
end subroutine pot22
