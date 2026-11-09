!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_brh2_leps_ltbz2_1995_dpem_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, u, ug )
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_legacy_units,       only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: legacy_triatomic_distances_to_bohr
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 2
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: u(nstates,nstates)
        double precision, intent(out) :: ug(nstates,nstates,natoms,3)
        double precision              :: d, g, p
        integer                       :: iatom, idir, nder
        double precision              :: nt, r(1,3), r2(3), v(1)
        double precision              :: u11(1), u12(1), u22(1)
        double precision              :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
        double precision              :: dudr(nstates,nstates,3)
        double precision              :: dudx(nstates,nstates,9)
        double precision              :: t(nstates,nstates)
        double precision              :: tmpmat(nstates,nstates)
        double precision              :: tx(9), drdx(3,9)
        integer                       :: i, j, k, l
        !initialize
        u  = 0.d0; p = 0.d0; g = 0.d0; d = 0.d0
        nt = 1
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is ClHH
        call legacy_triatomic_distances_to_bohr(x, r(1,:))
        nder = igrad
        call pot(r, u11, de_u11, 1, 1)
        call pot(r, u12, de_u12, 1, 2)
        call pot(r, u22, de_u22, 1, 3)
        u(1,1)      = u11(1)*legacy_hartree_to_electronvolt
        u(1,2)      = u12(1)*legacy_hartree_to_electronvolt
        u(2,1)      = u12(1)*legacy_hartree_to_electronvolt
        u(2,2)      = u22(1)*legacy_hartree_to_electronvolt
        dudr(1,1,:) = de_u11(:,1)*51.422067
        dudr(1,2,:) = de_u12(:,1)*51.422067
        dudr(2,1,:) = de_u12(:,1)*51.422067
        dudr(2,2,:) = de_u22(:,1)*51.422067
        r2          = r(1,:)*legacy_bohr_to_angstrom
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

!> @brief Preserve the historical Python-facing dpem entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_brh2_leps_ltbz2_1995_dpem_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 2
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: u(nstates,nstates)
    double precision, intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem

!> @brief Evaluate EvdRdX for the BrH2_LEPS_LTBZ2_1995_DPEM potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    implicit none
    integer, parameter            :: first_atoms(3)  = [1, 2, 1]
    integer, parameter            :: second_atoms(3) = [2, 3, 3]
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!> @brief Evaluate diagonalize for the BrH2_LEPS_LTBZ2_1995_DPEM potential.
subroutine diagonalize( n, a_ss, ev_s, u_ss )
    implicit none
    integer, intent(in)           :: n
    double precision, intent(in)  :: a_ss(n,n)
    double precision, intent(out) :: u_ss(n,n)
    double precision, intent(out) :: ev_s(n)
    integer                       :: io, i
    double precision, allocatable :: work_d(:)
    integer                       :: lwork_d
    u_ss    = a_ss
    lwork_d = -1
    allocate( work_d(1) )
    call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
    lwork_d = int(work_d(1))
    deallocate( work_d )
    allocate( work_d(lwork_d) )
    call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
end subroutine diagonalize

!

!> @brief Evaluate pot for the BrH2_LEPS_LTBZ2_1995_DPEM potential.
subroutine pot( r, e, de, nt, nsurf )
    implicit none
    double precision :: de, e, g_e, r
    dimension r(nt,3), e(nt), de(3,nt)
    integer            :: i, j, k, ldg_r, nsurf, nt
    integer, parameter :: g_p_  = 3
    integer, parameter :: ldg_e = 3
    dimension g_e(ldg_e,nt)
    if( nsurf .eq. 1 )then
        call g_pot11(g_p_, r, e, de, ldg_e, nt)
    else if( nsurf .eq. 2 )then
        call pot12(r, e, de, nt)
    else if( nsurf .eq. 3 )then
        call pot22(r, e, de, nt)
    endif
end subroutine pot

!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_brh2_leps_ltbz2_1995_dpem_g_pot11_entry_sta_432f3192
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, rr, e, g_e, ldg_e, n )
        !   this is a vectorized version of the BrH2SEC potential with derivativ
        !
        !   System:          BrH2
        !   Functional form: extended-LEPS plus three-center term with
        !                    internal-angle-dependent Sato parameters
        !   Common name: BrH2SECD
        !   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, et. al.
        !                    J. Phys. Chem. 99, 207-225 (1995)
        !   Analytic derivatives added, Yuri Volobuev, July-Aug 1998
        !
        !
        !   The potential parameters are included in DATA statements.
        !   For the this potential energy surface the potential energy in the
        !   BC valley is equal to the energy of the HH diatomic.
        !   The energy in the AC valley is equal to the energy of the HBr diatom
        !   which is also the energy value in the BA valley.
        !
        !   This potential is written such that:
        !                  R(1) = R(Br-H)
        !                  R(2) = R(H-H)
        !                  R(3) = R(H-Br)
        !   The zero of energy is defined at Br "infinitely" far from the H2 dia
        !
        implicit none
        optional :: g_p_, rr, e, g_e, ldg_e, n
        double precision :: a1, a2, a3, a4, a5, a6, b, beta, betaa, capr, capr2, coschi, &
            & costh, coul, d, de, deq, derad, dist, do4z, dv3c, dxdr, dz, e, exch, fcn4, op3z, pe, &
            & r132p2, r1m3, r1m3p2, r1p2, r1p3, r1p3m2, r2, r2i, r2p2, r3p2, rad, re, rr, schip2, sintp2, &
            & sintp6, sintp8, top3z, tzp3, v3c, v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, x, xlamb, z, &
            & zchi, zhbrp1, zhbrp2, zhhp1, zhhp2, zp3, zpo
        integer          :: i, k, n
        double precision :: r(3)
        !
        double precision, parameter :: cangau = 0.52917706d0
        double precision, parameter :: charkc = 627.5095d0
        !
        dimension z(3), dz(3,3)
        dimension dv3c(3)
        dimension x(3), coul(3), exch(3)
        dimension zpo(3), op3z(3), zp3(3), tzp3(3), top3z(3), do4z(3), beta(3)
        dimension dxdr(3,3), deq(3), derad(3)
        dimension de(3)
        !
        dimension d(3), re(3), betaa(3), b(3), rr(:,:), e(:)
        !
        !
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_rr, ldg_e
        double precision   :: d3_p, d2_p, d2_v, d10_v, d2_b, d3_v, d4_v, &
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
        a1       = 0.d0
        a2       = 0.d0
        a3       = 0.d0
        a4       = 0.d0
        a5       = 0.d0
        a6       = 0.d0
        b        = 0.d0
        beta     = 0.d0
        betaa    = 0.d0
        capr     = 0.d0
        capr2    = 0.d0
        coschi   = 0.d0
        costh    = 0.d0
        coul     = 0.d0
        d        = 0.d0
        d10_b    = 0.d0
        d10_v    = 0.d0
        d1_p     = 0.d0
        d1_w     = 0.d0
        d2_b     = 0.d0
        d2_p     = 0.d0
        d2_v     = 0.d0
        d3_b     = 0.d0
        d3_p     = 0.d0
        d3_v     = 0.d0
        d4_b     = 0.d0
        d4_v     = 0.d0
        d5_b     = 0.d0
        d5_v     = 0.d0
        d6_b     = 0.d0
        d6_v     = 0.d0
        d7_b     = 0.d0
        d7_v     = 0.d0
        d8_b     = 0.d0
        d8_v     = 0.d0
        de       = 0.d0
        deq      = 0.d0
        derad    = 0.d0
        dist     = 0.d0
        do4z     = 0.d0
        dv3c     = 0.d0
        dxdr     = 0.d0
        dz       = 0.d0
        exch     = 0.d0
        fcn4     = 0.d0
        g_beta   = 0.d0
        g_capr   = 0.d0
        g_capr2  = 0.d0
        g_coschi = 0.d0
        g_costh  = 0.d0
        g_coul   = 0.d0
        g_d1_w   = 0.d0
        g_dist   = 0.d0
        g_do4z   = 0.d0
        g_dv3c   = 0.d0
        g_dz     = 0.d0
        g_exch   = 0.d0
        g_fcn4   = 0.d0
        g_op3z   = 0.d0
        g_pe     = 0.d0
        g_r      = 0.d0
        g_r132p2 = 0.d0
        g_r1m3   = 0.d0
        g_r1m3p2 = 0.d0
        g_r1p2   = 0.d0
        g_r1p3   = 0.d0
        g_r1p3m2 = 0.d0
        g_r2p2   = 0.d0
        g_r3p2   = 0.d0
        g_rad    = 0.d0
        g_schip2 = 0.d0
        g_sintp2 = 0.d0
        g_sintp6 = 0.d0
        g_sintp8 = 0.d0
        g_top3z  = 0.d0
        g_tzp3   = 0.d0
        g_v3c    = 0.d0
        g_v3cp1  = 0.d0
        g_v3cp2  = 0.d0
        g_v3cp3  = 0.d0
        g_v3cp4  = 0.d0
        g_v3cp5  = 0.d0
        g_v3cp6  = 0.d0
        g_x      = 0.d0
        g_z      = 0.d0
        g_zchi   = 0.d0
        g_zhbrp1 = 0.d0
        g_zhbrp2 = 0.d0
        g_zhhp1  = 0.d0
        g_zhhp2  = 0.d0
        g_zp3    = 0.d0
        g_zpo    = 0.d0
        op3z     = 0.d0
        pe       = 0.d0
        r        = 0.d0
        r132p2   = 0.d0
        r1m3     = 0.d0
        r1m3p2   = 0.d0
        r1p2     = 0.d0
        r1p3     = 0.d0
        r1p3m2   = 0.d0
        r2       = 0.d0
        r2i      = 0.d0
        r2p2     = 0.d0
        r3p2     = 0.d0
        rad      = 0.d0
        re       = 0.d0
        schip2   = 0.d0
        sintp2   = 0.d0
        sintp6   = 0.d0
        sintp8   = 0.d0
        top3z    = 0.d0
        tzp3     = 0.d0
        v3c      = 0.d0
        v3cp1    = 0.d0
        v3cp2    = 0.d0
        v3cp3    = 0.d0
        v3cp4    = 0.d0
        v3cp5    = 0.d0
        v3cp6    = 0.d0
        x        = 0.d0
        xlamb    = 0.d0
        z        = 0.d0
        zchi     = 0.d0
        zhbrp1   = 0.d0
        zhbrp2   = 0.d0
        zhhp1    = 0.d0
        zhhp2    = 0.d0
        zp3      = 0.d0
        zpo      = 0.d0
        d        = [ 90.4473d0, 109.559d0, 90.4473d0]
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        xlamb = .4568d0 / 27.21161d0
        a1    = v3cp1
        a2    = v3cp2
        a3    = v3cp3
        a4    = v3cp4
        a5    = v3cp5
        a6    = v3cp6
        re    = [ 1.41443d0, 0.74144d0, 1.41443d0]
        betaa = [ 1.81088d0, 1.94198d0, 1.81088d0]
        r2    = 1.41421356d0
        !
        zhhp1  = 2.91230499d-1
        zhhp2  = 4.99160974d0
        zhbrp1 = 3.84262775d-3
        zhbrp2 = 2.87642885d-3
        v3cp1  = 1.67457884d0
        v3cp2  = 7.09518875d-1
        v3cp3  = 1.87582679d-1
        v3cp4  = 2.43517074d-1
        v3cp5  = 2.04145206d-1
        v3cp6  = 2.57833022d0
        !
        !
        !   Convert the constants to atomic units
        !
        do i = 1, 3
            de(i)   = d(i) / charkc
            re(i)   = re(i) / cangau
            beta(i) = betaa(i) * cangau
        enddo
        !
        !
        !
        !
        !
        do k = 1, n
            do g_i_ = 1, g_p_
                g_r(g_i_, 1) = 0.0d0
            enddo
            r(1) = rr(k,1)
            !--------
            do g_i_ = 1, g_p_
                g_r(g_i_, 2) = 0.0d0
            enddo
            r(2) = rr(k,2)
            !--------
            do g_i_ = 1, g_p_
                g_r(g_i_, 3) = 0.0d0
            enddo
            r(3)      = rr(k,3)
            g_r(1, 1) = 1.0d0
            g_r(2, 2) = 1.0d0
            g_r(3, 3) = 1.0d0
            !--------
            !
            !   Compute bond dependent Sato parameters
            d2_b = r(1) + r(1)
            do g_i_ = 1, g_p_
                g_r1p2(g_i_) = d2_b * g_r(g_i_, 1)
            enddo
            r1p2 = r(1) * r(1)
            !--------
            d2_b = r(2) + r(2)
            do g_i_ = 1, g_p_
                g_r2p2(g_i_) = d2_b * g_r(g_i_, 2)
            enddo
            r2p2 = r(2) * r(2)
            !--------
            d2_b = r(3) + r(3)
            do g_i_ = 1, g_p_
                g_r3p2(g_i_) = d2_b * g_r(g_i_, 3)
            enddo
            r3p2 = r(3) * r(3)
            !--------
            r2i = 1.0d0 / r(2)
            do g_i_ = 1, g_p_
                g_zchi(g_i_) = 0.0d0
            enddo
            zchi = 0.0d0
            !--------
            !   CAPR IS THE Br TO H2 DISTANCE
            d5_b = (-0.5d0) * 0.5d0
            do g_i_ = 1, g_p_
                g_capr2(g_i_) = d5_b * g_r2p2(g_i_) + 0.5d0 * g_r3p2(g_i_) + 0.5d0 * g_r1p2(g_i_)
            enddo
            capr2 = 0.5d0 * (r1p2 + r3p2 - 0.5d0 * r2p2)
            !--------
            if( capr2 .le. 0.d0 )then
                do g_i_ = 1, g_p_
                    g_capr2(g_i_) = 0.0d0
                enddo
                capr2 = 0.0d0
                !--------
                do g_i_ = 1, g_p_
                    g_capr(g_i_) = 0.0d0
                enddo
                capr = 0.0d0
                !--------
            else
                d2_v = sqrt(capr2)
                if( capr2 .gt. 0.0d0 )then
                    d1_p = 1.0d0 / (2.0d0 *  d2_v)
                else
                    d1_p = 0.0d0
                endif
                do g_i_ = 1, g_p_
                    g_capr(g_i_) = d1_p * g_capr2(g_i_)
                enddo
                capr = d2_v
                !--------
                !   Compute Chi, the angle between CAPR and RHH
                d7_v = r(2) * capr
                d8_v = 0.5d0 * (r3p2 - r1p2) / d7_v
                d3_b = (-d8_v) / d7_v
                d4_b = d3_b * capr
                d5_b = d3_b * r(2)
                d6_b = 1.0d0 / d7_v * 0.5d0
                do g_i_ = 1, g_p_
                    g_coschi(g_i_) = d5_b * g_capr(g_i_) + d4_b * g_r(g_i_, 2) + (-d6_b) * g_r1p2(g_i_) + d6_b * &
                        & g_r3p2(g_i_)
                enddo
                coschi = d8_v
                !--------
                d3_b = -coschi + (-coschi)
                do g_i_ = 1, g_p_
                    g_schip2(g_i_) = d3_b * g_coschi(g_i_)
                enddo
                schip2 = 1.0d0 - coschi * coschi
                !--------
                d4_v = zhbrp2 * schip2
                d6_v = zhbrp1 * schip2 + d4_v * schip2
                d7_b = capr2 * d4_v + capr2 * schip2 * zhbrp2 + capr2 * zhbrp1
                do g_i_ = 1, g_p_
                    g_zchi(g_i_) = d7_b * g_schip2(g_i_) + d6_v * g_capr2(g_i_ )
                enddo
                zchi = capr2 * d6_v
                !--------
            endif
            do g_i_ = 1, g_p_
                g_z(g_i_, 1) = g_zchi(g_i_)
            enddo
            z(1) = 0.1675d0 + zchi
            !--------
            do g_i_ = 1, g_p_
                g_z(g_i_, 3) = g_z(g_i_, 1)
            enddo
            z(3) = z(1)
            !--------
            d3_v = zhhp2 + r2p2
            d4_v = zhhp1 * r2p2 / d3_v
            d4_b = (-d4_v) / d3_v + 1.0d0 / d3_v * zhhp1
            do g_i_ = 1, g_p_
                g_z(g_i_, 2) = d4_b * g_r2p2(g_i_)
            enddo
            z(2) = d4_v
            !--------
            !
            !   Compute the LEPS parameters
            !
            do i = 1, 3
                do g_i_ = 1, g_p_
                    g_zpo(g_i_, i) = g_z(g_i_, i)
                enddo
                zpo(i) = 1.0d0 + z(i)
                !--------
                do g_i_ = 1, g_p_
                    g_op3z(g_i_, i) = 3.0d0 * g_z(g_i_, i)
                enddo
                op3z(i) = 1.0d0 + 3.0d0 * z(i)
                !--------
                do g_i_ = 1, g_p_
                    g_top3z(g_i_, i) = 2.0d0 * g_op3z(g_i_, i)
                enddo
                top3z(i) = 2.0d0 * op3z(i)
                !--------
                do g_i_ = 1, g_p_
                    g_zp3(g_i_, i) = g_z(g_i_, i)
                enddo
                zp3(i) = z(i) + 3.0d0
                !--------
                do g_i_ = 1, g_p_
                    g_tzp3(g_i_, i) = 2.0d0 * g_zp3(g_i_, i)
                enddo
                tzp3(i) = 2.0d0 * zp3(i)
                !--------
                d2_v = de(i) / 4.0d0 / zpo(i)
                d2_b = (-d2_v) / zpo(i)
                do g_i_ = 1, g_p_
                    g_do4z(g_i_, i) = d2_b * g_zpo(g_i_, i)
                enddo
                do4z(i) = d2_v
                !--------
                b(i) = beta(i) * do4z(i) * 2.0d0
            enddo
            do g_i_ = 1, g_p_
                g_pe(g_i_) = 0.0d0
            enddo
            pe = 0.d0
            !--------
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
                !--------
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
                !--------
                d6_v = op3z(i) * x(i) - tzp3(i)
                d7_v = do4z(i) * d6_v
                d4_b = x(i) * d6_v
                d5_b = x(i) * do4z(i)
                d8_b = d5_b * x(i)
                d3_b = d7_v + d5_b * op3z(i)
                do g_i_ = 1, g_p_
                    g_exch(g_i_, i) = (-d5_b) * g_tzp3(g_i_, i) + d3_b * g_x(g_i_, i) + d8_b * g_op3z(g_i_, i) + &
                        & d4_b * g_do4z(g_i_, i)
                enddo
                exch(i) = d7_v * x(i)
                !--------
                do g_i_ = 1, g_p_
                    g_pe(g_i_) = g_coul(g_i_, i) + g_pe(g_i_)
                enddo
                pe = pe + coul(i)
                !--------
            enddo
            d4_v  = (exch(1) - exch(2)) * (exch(1) - exch(2))
            d3_p  = 2.0d0 * (exch(1) - exch(2))
            d7_v  = (exch(2) - exch(3)) * (exch(2) - exch(3))
            d2_p  = 2.0d0 * (exch(2) - exch(3))
            d10_v = (exch(3) - exch(1)) * (exch(3) - exch(1))
            d1_p  = 2.0d0 * (exch(3) - exch(1))
            d5_b  = d1_p + (-d2_p)
            d6_b  = -d1_p + d3_p
            d10_b = d2_p + (-d3_p)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_exch(g_i_, 3) + d10_b * g_exch(g_i_, 2) + d6_b * g_exch(g_i_, 1)
            enddo
            d1_w = d4_v + d7_v + d10_v
            d2_v = sqrt(d1_w)
            if( d1_w .gt. 0.0d0 )then
                d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
                d1_p = 0.0d0
            endif
            do g_i_ = 1, g_p_
                g_rad(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            rad = d2_v
            !--------
            d5_b = -(1.0d0 / r2)
            do g_i_ = 1, g_p_
                g_pe(g_i_) = d5_b * g_rad(g_i_) + g_pe(g_i_)
            enddo
            pe = pe - rad / r2 + de(2)
            !--------
            !   Compute the three-center term
            do g_i_ = 1, g_p_
                g_r1p3(g_i_) = g_r(g_i_, 3) + g_r(g_i_, 1)
            enddo
            r1p3 = r(1) + r(3)
            !--------
            a1 = v3cp1
            a2 = v3cp2
            a3 = v3cp3
            a4 = v3cp4
            a5 = v3cp5
            a6 = v3cp6
            do g_i_ = 1, g_p_
                g_r1m3(g_i_) = -g_r(g_i_, 3) + g_r(g_i_, 1)
            enddo
            r1m3 = r(1) - r(3)
            !--------
            do g_i_ = 1, g_p_
                g_r1p3m2(g_i_) = -g_r(g_i_, 2) + g_r1p3(g_i_)
            enddo
            r1p3m2 = r1p3 - r(2)
            !--------
            d2_b = r1m3 + r1m3
            do g_i_ = 1, g_p_
                g_r1m3p2(g_i_) = d2_b * g_r1m3(g_i_)
            enddo
            r1m3p2 = r1m3 * r1m3
            !--------
            d2_b = r1p3m2 + r1p3m2
            do g_i_ = 1, g_p_
                g_r132p2(g_i_) = d2_b * g_r1p3m2(g_i_)
            enddo
            r132p2 = r1p3m2 * r1p3m2
            !--------
            !   Compute the bond angle from the three internal coordinates
            d4_b = -r(2) + (-r(2))
            d7_b = r(3) + r(3)
            d8_b = r(1) + r(1)
            do g_i_ = 1, g_p_
                g_dist(g_i_) = d4_b * g_r(g_i_, 2) + d7_b * g_r(g_i_, 3) + d8_b * g_r(g_i_, 1)
            enddo
            dist = r(1) * r(1) + r(3) * r(3) - r(2) * r(2)
            !--------
            d3_v = 2.0d0 * r(1)
            d5_v = d3_v * r(3)
            d6_v = dist / d5_v
            d2_b = 1.0d0 / d5_v
            d3_b = (-d6_v) / d5_v
            d5_b = d3_b * d3_v
            d6_b = d3_b * r(3) * 2.0d0
            do g_i_ = 1, g_p_
                g_costh(g_i_) = d5_b * g_r(g_i_, 3) + d6_b * g_r(g_i_, 1) + d2_b * g_dist(g_i_)
            enddo
            costh = d6_v
            !--------
            d3_b = -costh + (-costh)
            do g_i_ = 1, g_p_
                g_sintp2(g_i_) = d3_b * g_costh(g_i_)
            enddo
            sintp2 = 1.0d0 - costh * costh
            !--------
            d2_v = sintp2 * sintp2
            d3_b = d2_v + sintp2 * sintp2 + sintp2 * sintp2
            do g_i_ = 1, g_p_
                g_sintp6(g_i_) = d3_b * g_sintp2(g_i_)
            enddo
            sintp6 = d2_v * sintp2
            !--------
            do g_i_ = 1, g_p_
                g_sintp8(g_i_) = sintp6 * g_sintp2(g_i_) + sintp2 * g_sintp6 (g_i_)
            enddo
            sintp8 = sintp6 * sintp2
            !--------
            do g_i_ = 1, g_p_
                g_fcn4(g_i_) = a6 * g_sintp8(g_i_) + a5 * g_sintp2(g_i_)
            enddo
            fcn4 = 1.d0 + a5 * sintp2 + a6 * sintp8
            !--------
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
            !--------
            !   Sum all the terms that make up the energy
            do g_i_ = 1, g_p_
                g_pe(g_i_) = g_v3c(g_i_) + g_pe(g_i_)
            enddo
            pe = pe + v3c
            !--------
            do g_i_ = 1, g_p_
                g_e(g_i_,k) = g_pe(g_i_)
            enddo
            e(k) = pe + xlamb / 3.d0
            !--------
        enddo
        !
    end subroutine dispatch
end module chempotpy_brh2_leps_ltbz2_1995_dpem_g_pot11_entry_sta_432f3192

!> @brief Compatibility interface for g_pot11.
subroutine g_pot11( g_p_, rr, e, g_e, ldg_e, n )
    use chempotpy_brh2_leps_ltbz2_1995_dpem_g_pot11_entry_sta_432f3192, only: dispatch
    implicit none
    integer          :: g_p_
    double precision :: rr(n,3)
    double precision :: e(n)
    double precision :: g_e(ldg_e,*)
    integer          :: ldg_e
    integer          :: n
    call dispatch(g_p_=g_p_, rr=rr, e=e, g_e=g_e(1:ldg_e,1:n), ldg_e=ldg_e, n=n)
end subroutine g_pot11

!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_brh2_leps_ltbz2_1995_dpem_pot12_legacy_state_432f3192
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( rr, e, de, n )
        implicit none
        optional :: rr, e, de, n
        double precision :: de, e, rr, xlamb
        integer          :: j, k, n
        dimension rr(:,:), e(:), de(:,:)
        xlamb = 0.d0
        xlamb = 0.4568d0/27.21161d0
        !
        do k=1, n
            e(k) = -dsqrt(2.d0)*xlamb/3.d0
            do j =1, 3
                de(j,k) = 0.0d0
            enddo
        enddo
    end subroutine dispatch
end module chempotpy_brh2_leps_ltbz2_1995_dpem_pot12_legacy_state_432f3192

!> @brief Compatibility interface for pot12.
subroutine pot12( rr, e, de, n )
    use chempotpy_brh2_leps_ltbz2_1995_dpem_pot12_legacy_state_432f3192, only: dispatch
    implicit none
    double precision :: rr(n,3)
    double precision :: e(*)
    double precision :: de(3,n)
    integer          :: n
    call dispatch(rr=rr, e=e(1:n), de=de, n=n)
end subroutine pot12

!
!
!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_brh2_leps_ltbz2_1995_dpem_pot22_legacy_state_432f3192
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( rr, e, der, n )
        implicit none
        optional :: rr, e, der, n
        double precision :: b, beta, betaa, ccon, coul, d, de, der, do4z, dz, e, exch, &
            & op3z, pe, r2, rad, re, rr, top3z, tzp3, vanti1, vanti3, vmorse, x, xlamb, y1, y2, y3, z, &
            & zhbrp1, zhbrp2, zhhp1, zhhp2, zp3, zpo
        integer :: i, k, n
        !
        double precision            :: r(3)
        double precision, parameter :: cangau = 0.52917706d0
        double precision, parameter :: charkc = 627.5095d0
        !
        dimension z(3), dz(3,3)
        dimension x(3), coul(3), exch(3)
        dimension zpo(3), op3z(3), zp3(3), tzp3(3), top3z(3), do4z(3), beta(3)
        dimension de(3)
        integer :: g_myrank, g_nprocs
        !
        dimension d(3), re(3), betaa(3), b(3), rr(:,:), e(:), der(:,:)
        !
        !
        !
        b      = 0.d0
        beta   = 0.d0
        betaa  = 0.d0
        ccon   = 0.d0
        coul   = 0.d0
        d      = 0.d0
        de     = 0.d0
        do4z   = 0.d0
        dz     = 0.d0
        exch   = 0.d0
        op3z   = 0.d0
        pe     = 0.d0
        r      = 0.d0
        r2     = 0.d0
        rad    = 0.d0
        re     = 0.d0
        top3z  = 0.d0
        tzp3   = 0.d0
        vanti1 = 0.d0
        vanti3 = 0.d0
        vmorse = 0.d0
        x      = 0.d0
        xlamb  = 0.d0
        y1     = 0.d0
        y2     = 0.d0
        y3     = 0.d0
        z      = 0.d0
        zhbrp1 = 0.d0
        zhbrp2 = 0.d0
        zhhp1  = 0.d0
        zhhp2  = 0.d0
        zp3    = 0.d0
        zpo    = 0.d0
        d      = [ 90.4473d0, 109.559d0, 90.4473d0]
        re     = [ 1.41443d0, 0.74144d0, 1.41443d0]
        betaa  = [ 1.81088d0, 1.94198d0, 1.81088d0]
        r2     = 1.41421356d0
        ccon   = 0.5d0
        xlamb  = 0.4568d0/27.21161d0
        zhhp1  = 2.91230499d-1
        zhhp2  = 4.99160974d0
        zhbrp1 = 3.84262775d-3
        zhbrp2 = 2.87642885d-3
        !
        !
        !   Convert the constants to atomic units
        !
        do i = 1, 3
            de(i)   = d(i) / charkc
            re(i)   = re(i) / cangau
            beta(i) = betaa(i) * cangau
        enddo
        !
        !
        !
        do k=1, n
            r(1)     = rr(k,1)
            r(2)     = rr(k,2)
            r(3)     = rr(k,3)
            y2       = exp(-beta(2)*(r(2)-re(2)))
            y1       = exp(-beta(1)*(r(1)-re(1)))
            y3       = exp(-beta(3)*(r(3)-re(3)))
            vmorse   = de(2)*(y2**2-2.d0*y2)
            vanti1   = de(1)*(y1**2+2.d0*y1)
            vanti3   = de(3)*(y3**2+2.d0*y3)
            pe       = de(2)+vmorse+(vanti1+vanti3)*ccon+2.d0*xlamb/3.d0
            der(1,k) = -2.0d0*beta(1)*de(1)*(y1**2+y1)*ccon
            der(2,k) = -2.0d0*beta(2)*de(2)*(y2**2-y2)
            der(3,k) = -2.0d0*beta(3)*de(3)*(y3**2+y3)*ccon
            e(k)     = pe
        enddo
        !
    end subroutine dispatch
end module chempotpy_brh2_leps_ltbz2_1995_dpem_pot22_legacy_state_432f3192

!> @brief Compatibility interface for pot22.
subroutine pot22( rr, e, der, n )
    use chempotpy_brh2_leps_ltbz2_1995_dpem_pot22_legacy_state_432f3192, only: dispatch
    implicit none
    double precision :: rr(n,3)
    double precision :: e(*)
    double precision :: der(3,n)
    integer          :: n
    call dispatch(rr=rr, e=e(1:n), der=der, n=n)
end subroutine pot22
