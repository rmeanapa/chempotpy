!   System:                        NaH2
!   Functional form:               Modified Extended LEPS (2x2 diabatic fit)
!   Common name:                   NaH2 surface set 6
!   Interface:                     3-2V
!   Number of electronic surfaces: 2
!   Number of derivatives:         1
!   References:      M. D. Hack and D. G. Truhlar, J. Chem. Phys.
!                    112, 9716 (2000)
!   Notes:   This is the most accurate NaH2 fit.
!   Protocol: POT evaluates the stateless diabatic surfaces.
!               of the energy with respect to the coordinates for a give
!               geometric configuration
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      Ground diabatic electronic state, first-excited diabatic electron
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Na
!      infinitely far from the H2 diatomic and R(H2) set equal to the
!      H2 equilibrium diatomic value.
!   Parameters:
!   Coordinates:
!      Internal, Definition: R(1) = R(Na-H)
!                            R(2) = R(H-H')
!                            R(3) = R(Na-H')
!  The final version of surface 6, completed on 9/24/98
!  This version contains derivatives, to be used with trajectory
!  and other codes . . .
!  It is based on version 5g of the NaH2 potential, in turn based on 5f.
!  Reference:  P. Halvick and D. G. Truhlar, J. Chem. Phys. 96,2895 (199
!  The "zero" of energy occurs when the Na atom is "infinitely"
!  far from the H2 diatom and the H2 diatom is at its equilibrium geomet

!> @brief Implement NaH2_LEPS_6_2000, a potential energy surface for NaH2.
module chempotpy_nah2_leps_6_2000_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NaH2_LEPS_6_2000 potential energy surface for NaH2.
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
        integer               :: nt
        real(wp)              :: r(1,3), r2(3), v(1)
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
end module chempotpy_nah2_leps_6_2000_pes

!> @brief Expose NaH2_LEPS_6_2000 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_nah2_leps_6_2000_pes, only: evaluate
    use chempotpy_kinds,                only: wp
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

!> @brief Evaluate EvdRdX for the NaH2_LEPS_6_2000 potential.
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

!> @brief Evaluate diagonalize for the NaH2_LEPS_6_2000 potential.
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

!> @brief Evaluate pot for the NaH2_LEPS_6_2000 potential.
subroutine pot( r, e, de, nt, nsurf )
    use chempotpy_kinds, only: wp
    implicit none
    integer, parameter :: g_p_  = 3
    integer, parameter :: ldg_r = 3
    integer, parameter :: ldg_e = 3
    real(wp)           :: de(3,nt), e(nt), g_e(ldg_e,nt), g_r(ldg_r,nt,3), r(nt,3)
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
    endif
end subroutine pot

!  U11 surface - lower diabatic surface
!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_nah2_leps_6_2000_g_pot11_legacy_state_a1118117
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
        real(wp) :: al1, al2, alb, alhh, alt, balp, brho, c2a, c2b, c2c, c6, c60, c61, c63, c690, c6swt1, c6swt3, cevau, &
            & cind1, cind3, cosa, coscos, coul1, coul2, coul3, cplg2, cstsq, cstsq1, cstsq3, dalp, deh2, dexp1, dexp3, div1, &
            & div3, drho, e(:), elec, exch1, exch2, exch3, f0, f90, frac, h2a1, h2a2, h2a3, h2a4, h2a5, h2a6, h2a7, h2a8, h2a9, &
            & h2b1, h2b2, h2b3, h2b4, h2b5, h2b6, h2b7, h2b8, h2b9, h2i1, h2i2, h2i3, h2i4, h2i5, h2i6, h2ie, hie, hnaa1, hnaa2, &
            & hnaa3, hnaa4, hnaa5, hnab1, hnab2, hnab3, hnab4, hnab5, hnai1, hnai2, hnai3, hnai4, hnai5, hnai6, hnau1, hnau2, &
            & hnau3, hnau4, hnau5, hnau6, hpol, onethd, r(:,:), r1s, r2s, r3s, rac2, rb0, rb90, rd, rexp1, rexp2, rexp3, &
            & rmuimu1, rmuimu3, rnag, rnag1, rnag1c, rnag3, rnag3c, rnagc, rnags, rnahie1, rnahie3, rnahmu1, rnahmu3, rnahpar1, &
            & rnahpar3, rnahperp1, rnahperp3, rnahpol1, rnahpol1b, rnahpol3, rnahpol3b, rnaie, rnapol, rsm1, rsm2, rsm3, s1, &
            & s10, s101i, s102i, s103i, s104i, s105i, s106i, s10alp, s10c, s10c1, s10c2, s10c3, s10c4, s10c5, s10c6, s10rho, &
            & s190, s190alp, s190c, s190c1, s190c2, s190c3, s190c6, s190rho, s1a, s1c1, s1c2, s1c3, s1c4, s1c5, s1c6, s1i, &
            & s1ialp, s1ic, s2, s20, s20alp, s20c, s20c1, s20c2, s20c3, s20c4, s20c5, s20c6, s20rho, s290, s290alp, s290c, &
            & s290c1, s290c2, s290c3, s290c4, s290c5, s290c6, s290rho, s2a, s2c1, s2c2, s2c3, s2c4, s2c5, s2c6, s3, s30, s30c, &
            & s390, s390c, s3a, s3i, s3ic, t1, t10, t10ac1, t10ac2, t10ac3, t10ac4, t10ac5, t10ac6, t10alp, t10alt, t10c1, &
            & t10c2, t10c3, t10c4, t10c5, t10c6, t10rho, t190, t190c1, t190c2, t190c3, t190c4, t190c5, t190c6, t1i, t1ic1, &
            & t1ic2, t1ic3, t1ic4, t1ic5, t1ic6, t2, t20, t20alp, t20c, t20c1, t20c2, t20c3, t20c4, t20c5, t20c6, t20rho, t2a, &
            & t2c1, t2c2, t2c3, t2c4, t2c5, t2c6, t3, t30, t30alt, t390, t3i, twothd, w
        integer :: i, j, k, nsurf, nt
        ! data blocks for the NaH singlets . . .
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d3_p, d24_b, d2_p, d15_v, d20_v, d16_v, d22_v, &
            & d2_v, d16_b, d2_b
        real(wp) :: d1_w, d2_w, d3_v, d4_v, d5_v, d6_v, d3_b, d4_b, &
            & d5_b, d6_b
        real(wp) :: d1_p, d7_v, d8_v, d9_v, d10_v, d7_b, d8_b, d9_b, d10_b, d11_v
        real(wp) :: d12_v, d13_v, d19_v, d11_b, d12_b, d13_b, d14_b, g_r1s(g_pmax_), g_r(:,:,:), g_r2s(g_pmax_)
        real(wp) :: g_r3s(g_pmax_), g_d2_w(g_pmax_), g_d1_w(g_pmax_ &
            & ), g_rnag(g_pmax_), g_rnags(g_pmax_), g_rnagc(g_pmax_), g_rnag1(g_pmax_), g_rnag3(g_pmax_), g_rnag1c(g_pmax_), g_rnag3c(g_pmax_)
        real(wp) :: g_rsm1(g_pmax_), g_rsm2(g_pmax_), g_rsm3(g_pmax_), g_cstsq(g_pmax_), g_cosa(g_pmax_), g_coscos(g_pmax_), g_cstsq1 &
            & (g_pmax_), g_cstsq3(g_pmax_), g_s1a(g_pmax_), g_s10c(g_pmax_)
        real(wp) :: g_s10(g_pmax_), g_s1ic(g_pmax_), g_s1i(g_pmax_), g_s190c(g_pmax_), g_s190(g_pmax_), g_s1(g_pmax_), g_t190(g_pmax_ &
            & ), g_t10(g_pmax_), g_t10alt(g_pmax_), g_t1i(g_pmax_)
        real(wp) :: g_t1(g_pmax_), g_s2a(g_pmax_), g_s290c(g_pmax_), g_s290(g_pmax_), g_s20c(g_pmax_), g_s20(g_pmax_), g_s2(g_pmax_), &
            & g_t2a(g_pmax_), g_t20c(g_pmax_), g_t20(g_pmax_)
        real(wp) :: g_t2(g_pmax_), g_s3a(g_pmax_), g_s30c(g_pmax_), &
            & g_s30(g_pmax_), g_s3ic(g_pmax_), g_s3i(g_pmax_), g_s390c(g_pmax_), g_s390(g_pmax_), g_s3(g_pmax_), g_t390(g_pmax_)
        real(wp) :: g_t30(g_pmax_), g_t30alt(g_pmax_), g_t3(g_pmax_ &
            & ), g_t3i(g_pmax_), g_h2ie(g_pmax_), g_al1(g_pmax_), g_al2(g_pmax_), g_rnahie1(g_pmax_), g_rnahie3(g_pmax_), g_rnahmu1(g_pmax_)
        real(wp) :: g_rnahmu3(g_pmax_), g_rnahpar1(g_pmax_), g_rnahpar3(g_pmax_), g_rnahperp1(g_pmax_), g_rnahperp3(g_pmax_), g_rnahpol1(g_pmax_), g_rnahpol1b(g_pmax_), g_rnahpol3(g_pmax_), g_rnahpol3b(g_pmax_), g_div1(g_pmax_)
        real(wp) :: g_rmuimu1(g_pmax_), g_div3(g_pmax_), g_rmuimu3( &
            & g_pmax_), g_c690(g_pmax_), g_c60(g_pmax_), g_c6(g_pmax_), g_c61(g_pmax_), g_dexp1(g_pmax_), g_c63(g_pmax_), g_dexp3(g_pmax_)
        real(wp) :: g_cind1(g_pmax_), g_cind3(g_pmax_), g_c6swt1(g_pmax_), g_c6swt3(g_pmax_), g_elec(g_pmax_), g_coul1(g_pmax_), g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_)
        real(wp) :: g_exch3(g_pmax_), g_w(g_pmax_), g_cplg2(g_pmax_ &
            & ), g_e(:,:)
        integer :: g_ehfid
        ! data blocks for the NaH triplets
        ! data blocks for the H2 singlets
        ! data blocks for the h2 triplets
        ! data blocks for the long range forces
        ! misc constants
        al1         = 0.e0_wp
        al2         = 0.e0_wp
        alb         = 0.e0_wp
        alhh        = 0.e0_wp
        alt         = 0.e0_wp
        balp        = 0.e0_wp
        brho        = 0.e0_wp
        c2a         = 0.e0_wp
        c2b         = 0.e0_wp
        c2c         = 0.e0_wp
        c6          = 0.e0_wp
        c60         = 0.e0_wp
        c61         = 0.e0_wp
        c63         = 0.e0_wp
        c690        = 0.e0_wp
        c6swt1      = 0.e0_wp
        c6swt3      = 0.e0_wp
        cevau       = 0.e0_wp
        cind1       = 0.e0_wp
        cind3       = 0.e0_wp
        cosa        = 0.e0_wp
        coscos      = 0.e0_wp
        coul1       = 0.e0_wp
        coul2       = 0.e0_wp
        coul3       = 0.e0_wp
        cplg2       = 0.e0_wp
        cstsq       = 0.e0_wp
        cstsq1      = 0.e0_wp
        cstsq3      = 0.e0_wp
        d10_b       = 0.e0_wp
        d10_v       = 0.e0_wp
        d11_b       = 0.e0_wp
        d11_v       = 0.e0_wp
        d12_b       = 0.e0_wp
        d12_v       = 0.e0_wp
        d13_b       = 0.e0_wp
        d13_v       = 0.e0_wp
        d14_b       = 0.e0_wp
        d15_v       = 0.e0_wp
        d16_b       = 0.e0_wp
        d16_v       = 0.e0_wp
        d19_v       = 0.e0_wp
        d1_p        = 0.e0_wp
        d1_w        = 0.e0_wp
        d20_v       = 0.e0_wp
        d22_v       = 0.e0_wp
        d24_b       = 0.e0_wp
        d2_b        = 0.e0_wp
        d2_p        = 0.e0_wp
        d2_v        = 0.e0_wp
        d2_w        = 0.e0_wp
        d3_b        = 0.e0_wp
        d3_p        = 0.e0_wp
        d3_v        = 0.e0_wp
        d4_b        = 0.e0_wp
        d4_v        = 0.e0_wp
        d5_b        = 0.e0_wp
        d5_v        = 0.e0_wp
        d6_b        = 0.e0_wp
        d6_v        = 0.e0_wp
        d7_b        = 0.e0_wp
        d7_v        = 0.e0_wp
        d8_b        = 0.e0_wp
        d8_v        = 0.e0_wp
        d9_b        = 0.e0_wp
        d9_v        = 0.e0_wp
        dalp        = 0.e0_wp
        deh2        = 0.e0_wp
        dexp1       = 0.e0_wp
        dexp3       = 0.e0_wp
        div1        = 0.e0_wp
        div3        = 0.e0_wp
        drho        = 0.e0_wp
        elec        = 0.e0_wp
        exch1       = 0.e0_wp
        exch2       = 0.e0_wp
        exch3       = 0.e0_wp
        f0          = 0.e0_wp
        f90         = 0.e0_wp
        frac        = 0.e0_wp
        g_al1       = 0.e0_wp
        g_al2       = 0.e0_wp
        g_c6        = 0.e0_wp
        g_c60       = 0.e0_wp
        g_c61       = 0.e0_wp
        g_c63       = 0.e0_wp
        g_c690      = 0.e0_wp
        g_c6swt1    = 0.e0_wp
        g_c6swt3    = 0.e0_wp
        g_cind1     = 0.e0_wp
        g_cind3     = 0.e0_wp
        g_cosa      = 0.e0_wp
        g_coscos    = 0.e0_wp
        g_coul1     = 0.e0_wp
        g_coul2     = 0.e0_wp
        g_coul3     = 0.e0_wp
        g_cplg2     = 0.e0_wp
        g_cstsq     = 0.e0_wp
        g_cstsq1    = 0.e0_wp
        g_cstsq3    = 0.e0_wp
        g_d1_w      = 0.e0_wp
        g_d2_w      = 0.e0_wp
        g_dexp1     = 0.e0_wp
        g_dexp3     = 0.e0_wp
        g_div1      = 0.e0_wp
        g_div3      = 0.e0_wp
        g_elec      = 0.e0_wp
        g_exch1     = 0.e0_wp
        g_exch2     = 0.e0_wp
        g_exch3     = 0.e0_wp
        g_h2ie      = 0.e0_wp
        g_r1s       = 0.e0_wp
        g_r2s       = 0.e0_wp
        g_r3s       = 0.e0_wp
        g_rmuimu1   = 0.e0_wp
        g_rmuimu3   = 0.e0_wp
        g_rnag      = 0.e0_wp
        g_rnag1     = 0.e0_wp
        g_rnag1c    = 0.e0_wp
        g_rnag3     = 0.e0_wp
        g_rnag3c    = 0.e0_wp
        g_rnagc     = 0.e0_wp
        g_rnags     = 0.e0_wp
        g_rnahie1   = 0.e0_wp
        g_rnahie3   = 0.e0_wp
        g_rnahmu1   = 0.e0_wp
        g_rnahmu3   = 0.e0_wp
        g_rnahpar1  = 0.e0_wp
        g_rnahpar3  = 0.e0_wp
        g_rnahperp1 = 0.e0_wp
        g_rnahperp3 = 0.e0_wp
        g_rnahpol1  = 0.e0_wp
        g_rnahpol1b = 0.e0_wp
        g_rnahpol3  = 0.e0_wp
        g_rnahpol3b = 0.e0_wp
        g_rsm1      = 0.e0_wp
        g_rsm2      = 0.e0_wp
        g_rsm3      = 0.e0_wp
        g_s1        = 0.e0_wp
        g_s10       = 0.e0_wp
        g_s10c      = 0.e0_wp
        g_s190      = 0.e0_wp
        g_s190c     = 0.e0_wp
        g_s1a       = 0.e0_wp
        g_s1i       = 0.e0_wp
        g_s1ic      = 0.e0_wp
        g_s2        = 0.e0_wp
        g_s20       = 0.e0_wp
        g_s20c      = 0.e0_wp
        g_s290      = 0.e0_wp
        g_s290c     = 0.e0_wp
        g_s2a       = 0.e0_wp
        g_s3        = 0.e0_wp
        g_s30       = 0.e0_wp
        g_s30c      = 0.e0_wp
        g_s390      = 0.e0_wp
        g_s390c     = 0.e0_wp
        g_s3a       = 0.e0_wp
        g_s3i       = 0.e0_wp
        g_s3ic      = 0.e0_wp
        g_t1        = 0.e0_wp
        g_t10       = 0.e0_wp
        g_t10alt    = 0.e0_wp
        g_t190      = 0.e0_wp
        g_t1i       = 0.e0_wp
        g_t2        = 0.e0_wp
        g_t20       = 0.e0_wp
        g_t20c      = 0.e0_wp
        g_t2a       = 0.e0_wp
        g_t3        = 0.e0_wp
        g_t30       = 0.e0_wp
        g_t30alt    = 0.e0_wp
        g_t390      = 0.e0_wp
        g_t3i       = 0.e0_wp
        g_w         = 0.e0_wp
        h2a1        = 0.e0_wp
        h2a2        = 0.e0_wp
        h2a3        = 0.e0_wp
        h2a4        = 0.e0_wp
        h2a5        = 0.e0_wp
        h2a6        = 0.e0_wp
        h2a7        = 0.e0_wp
        h2a8        = 0.e0_wp
        h2a9        = 0.e0_wp
        h2b1        = 0.e0_wp
        h2b2        = 0.e0_wp
        h2b3        = 0.e0_wp
        h2b4        = 0.e0_wp
        h2b5        = 0.e0_wp
        h2b6        = 0.e0_wp
        h2b7        = 0.e0_wp
        h2b8        = 0.e0_wp
        h2b9        = 0.e0_wp
        h2i1        = 0.e0_wp
        h2i2        = 0.e0_wp
        h2i3        = 0.e0_wp
        h2i4        = 0.e0_wp
        h2i5        = 0.e0_wp
        h2i6        = 0.e0_wp
        h2ie        = 0.e0_wp
        hie         = 0.e0_wp
        hnaa1       = 0.e0_wp
        hnaa2       = 0.e0_wp
        hnaa3       = 0.e0_wp
        hnaa4       = 0.e0_wp
        hnaa5       = 0.e0_wp
        hnab1       = 0.e0_wp
        hnab2       = 0.e0_wp
        hnab3       = 0.e0_wp
        hnab4       = 0.e0_wp
        hnab5       = 0.e0_wp
        hnai1       = 0.e0_wp
        hnai2       = 0.e0_wp
        hnai3       = 0.e0_wp
        hnai4       = 0.e0_wp
        hnai5       = 0.e0_wp
        hnai6       = 0.e0_wp
        hnau1       = 0.e0_wp
        hnau2       = 0.e0_wp
        hnau3       = 0.e0_wp
        hnau4       = 0.e0_wp
        hnau5       = 0.e0_wp
        hnau6       = 0.e0_wp
        hpol        = 0.e0_wp
        onethd      = 0.e0_wp
        r1s         = 0.e0_wp
        r2s         = 0.e0_wp
        r3s         = 0.e0_wp
        rac2        = 0.e0_wp
        rb0         = 0.e0_wp
        rb90        = 0.e0_wp
        rd          = 0.e0_wp
        rexp1       = 0.e0_wp
        rexp2       = 0.e0_wp
        rexp3       = 0.e0_wp
        rmuimu1     = 0.e0_wp
        rmuimu3     = 0.e0_wp
        rnag        = 0.e0_wp
        rnag1       = 0.e0_wp
        rnag1c      = 0.e0_wp
        rnag3       = 0.e0_wp
        rnag3c      = 0.e0_wp
        rnagc       = 0.e0_wp
        rnags       = 0.e0_wp
        rnahie1     = 0.e0_wp
        rnahie3     = 0.e0_wp
        rnahmu1     = 0.e0_wp
        rnahmu3     = 0.e0_wp
        rnahpar1    = 0.e0_wp
        rnahpar3    = 0.e0_wp
        rnahperp1   = 0.e0_wp
        rnahperp3   = 0.e0_wp
        rnahpol1    = 0.e0_wp
        rnahpol1b   = 0.e0_wp
        rnahpol3    = 0.e0_wp
        rnahpol3b   = 0.e0_wp
        rnaie       = 0.e0_wp
        rnapol      = 0.e0_wp
        rsm1        = 0.e0_wp
        rsm2        = 0.e0_wp
        rsm3        = 0.e0_wp
        s1          = 0.e0_wp
        s10         = 0.e0_wp
        s101i       = 0.e0_wp
        s102i       = 0.e0_wp
        s103i       = 0.e0_wp
        s104i       = 0.e0_wp
        s105i       = 0.e0_wp
        s106i       = 0.e0_wp
        s10alp      = 0.e0_wp
        s10c        = 0.e0_wp
        s10c1       = 0.e0_wp
        s10c2       = 0.e0_wp
        s10c3       = 0.e0_wp
        s10c4       = 0.e0_wp
        s10c5       = 0.e0_wp
        s10c6       = 0.e0_wp
        s10rho      = 0.e0_wp
        s190        = 0.e0_wp
        s190alp     = 0.e0_wp
        s190c       = 0.e0_wp
        s190c1      = 0.e0_wp
        s190c2      = 0.e0_wp
        s190c3      = 0.e0_wp
        s190c6      = 0.e0_wp
        s190rho     = 0.e0_wp
        s1a         = 0.e0_wp
        s1c1        = 0.e0_wp
        s1c2        = 0.e0_wp
        s1c3        = 0.e0_wp
        s1c4        = 0.e0_wp
        s1c5        = 0.e0_wp
        s1c6        = 0.e0_wp
        s1i         = 0.e0_wp
        s1ialp      = 0.e0_wp
        s1ic        = 0.e0_wp
        s2          = 0.e0_wp
        s20         = 0.e0_wp
        s20alp      = 0.e0_wp
        s20c        = 0.e0_wp
        s20c1       = 0.e0_wp
        s20c2       = 0.e0_wp
        s20c3       = 0.e0_wp
        s20c4       = 0.e0_wp
        s20c5       = 0.e0_wp
        s20c6       = 0.e0_wp
        s20rho      = 0.e0_wp
        s290        = 0.e0_wp
        s290alp     = 0.e0_wp
        s290c       = 0.e0_wp
        s290c1      = 0.e0_wp
        s290c2      = 0.e0_wp
        s290c3      = 0.e0_wp
        s290c4      = 0.e0_wp
        s290c5      = 0.e0_wp
        s290c6      = 0.e0_wp
        s290rho     = 0.e0_wp
        s2a         = 0.e0_wp
        s2c1        = 0.e0_wp
        s2c2        = 0.e0_wp
        s2c3        = 0.e0_wp
        s2c4        = 0.e0_wp
        s2c5        = 0.e0_wp
        s2c6        = 0.e0_wp
        s3          = 0.e0_wp
        s30         = 0.e0_wp
        s30c        = 0.e0_wp
        s390        = 0.e0_wp
        s390c       = 0.e0_wp
        s3a         = 0.e0_wp
        s3i         = 0.e0_wp
        s3ic        = 0.e0_wp
        t1          = 0.e0_wp
        t10         = 0.e0_wp
        t10ac1      = 0.e0_wp
        t10ac2      = 0.e0_wp
        t10ac3      = 0.e0_wp
        t10ac4      = 0.e0_wp
        t10ac5      = 0.e0_wp
        t10ac6      = 0.e0_wp
        t10alp      = 0.e0_wp
        t10alt      = 0.e0_wp
        t10c1       = 0.e0_wp
        t10c2       = 0.e0_wp
        t10c3       = 0.e0_wp
        t10c4       = 0.e0_wp
        t10c5       = 0.e0_wp
        t10c6       = 0.e0_wp
        t10rho      = 0.e0_wp
        t190        = 0.e0_wp
        t190c1      = 0.e0_wp
        t190c2      = 0.e0_wp
        t190c3      = 0.e0_wp
        t190c4      = 0.e0_wp
        t190c5      = 0.e0_wp
        t190c6      = 0.e0_wp
        t1i         = 0.e0_wp
        t1ic1       = 0.e0_wp
        t1ic2       = 0.e0_wp
        t1ic3       = 0.e0_wp
        t1ic4       = 0.e0_wp
        t1ic5       = 0.e0_wp
        t1ic6       = 0.e0_wp
        t2          = 0.e0_wp
        t20         = 0.e0_wp
        t20alp      = 0.e0_wp
        t20c        = 0.e0_wp
        t20c1       = 0.e0_wp
        t20c2       = 0.e0_wp
        t20c3       = 0.e0_wp
        t20c4       = 0.e0_wp
        t20c5       = 0.e0_wp
        t20c6       = 0.e0_wp
        t20rho      = 0.e0_wp
        t2a         = 0.e0_wp
        t2c1        = 0.e0_wp
        t2c2        = 0.e0_wp
        t2c3        = 0.e0_wp
        t2c4        = 0.e0_wp
        t2c5        = 0.e0_wp
        t2c6        = 0.e0_wp
        t3          = 0.e0_wp
        t30         = 0.e0_wp
        t30alt      = 0.e0_wp
        t390        = 0.e0_wp
        t3i         = 0.e0_wp
        twothd      = 0.e0_wp
        w           = 0.e0_wp
        s1c1        = 593.88235e0_wp
        s1c2        = 1.72588e0_wp
        s1c3        = -251.80392e0_wp
        s1c4        = 89.98431e0_wp
        s1c5        = -8.83922e0_wp
        s1c6        = 1.383e0_wp
        s10c1       = 514.17323e0_wp
        s10c2       = 1.76693e0_wp
        s10c3       = -143.89764e0_wp
        s10c4       = 13.76378e0_wp
        s10c5       = -2.11024e0_wp
        s10c6       = 1.383e0_wp
        s10alp      = 0.79843e0_wp
        s10rho      = 2.47244e0_wp
        s190c1      = 274.97638e0_wp
        s190c2      = 2.20787e0_wp
        s190c3      = -28.97638e0_wp
        s190c6      = 1.383e0_wp
        s190alp     = 0.66693e0_wp
        s190rho     = 6.95276e0_wp
        s101i       = 479.85039e0_wp
        s102i       = 2.22598e0_wp
        s103i       = -24.20472e0_wp
        s104i       = 30.19685e0_wp
        s105i       = -9.96063e0_wp
        s106i       = 1.383e0_wp
        s1ialp      = 0.15e0_wp
        t10c1       = 335.59055e0_wp
        t10c2       = 4.02126e0_wp
        t10c3       = 78.81890e0_wp
        t10c4       = -27.51969e0_wp
        t10c5       = 4.22835e0_wp
        t10c6       = 1.383e0_wp
        t10ac1      = 365.90551e0_wp
        t10ac2      = 2.46220e0_wp
        t10ac3      = 49.75591e0_wp
        t10ac4      = -50.11024e0_wp
        t10ac5      = 17.32283e0_wp
        t10ac6      = 1.383e0_wp
        t10alp      = 0.77874e0_wp
        t10rho      = 4.46457e0_wp
        t190c1      = 511.51181e0_wp
        t190c2      = 1.91024e0_wp
        t190c3      = -8.53543e0_wp
        t190c4      = -1.58268e0_wp
        t190c5      = 1.53543e0_wp
        t190c6      = 1.383e0_wp
        t1ic1       = 430.31496e0_wp
        t1ic2       = 2.05197e0_wp
        t1ic3       = 24.40945e0_wp
        t1ic4       = 12.91339e0_wp
        t1ic5       = -6.07087e0_wp
        t1ic6       = 1.383e0_wp
        s2c1        = 139.7160e0_wp
        s2c2        = -123.8978e0_wp
        s2c3        = 3.4031e0_wp
        s2c4        = -6.8725e0_wp
        s2c5        = -23.0440e0_wp
        s2c6        = 2.032e0_wp
        s290c1      = 131.65354e0_wp
        s290c2      = 294.56693e0_wp
        s290c3      = 5.20472e0_wp
        s290c4      = -51.69291e0_wp
        s290c5      = 2.16535e0_wp
        s290c6      = 2.032e0_wp
        s290alp     = 0.09843e0_wp
        s290rho     = 6.8e0_wp
        s20c1       = 189.76378e0_wp
        s20c2       = -158.80315e0_wp
        s20c3       = 3.60945e0_wp
        s20c4       = -16.84252e0_wp
        s20c5       = -12.59843e0_wp
        s20c6       = 2.032e0_wp
        s20alp      = 0.8063e0_wp
        s20rho      = 3.37008e0_wp
        t2c1        = 147.214e0_wp
        t2c2        = 3.915e0_wp
        t2c3        = 41.275e0_wp
        t2c4        = 10.505e0_wp
        t2c5        = 4.408e0_wp
        t2c6        = 2.032e0_wp
        t20c1       = 654.51969e0_wp
        t20c2       = 2.45197e0_wp
        t20c3       = -170.03937e0_wp
        t20c4       = 108.26772e0_wp
        t20c5       = -17.70079e0_wp
        t20c6       = 2.032e0_wp
        t20alp      = 0.05236e0_wp
        t20rho      = 52.47683993e0_wp
        h2a1        = 2.38029900e0_wp
        h2a2        = 1.47400334e0_wp
        h2a3        = -0.226331771e0_wp
        h2a4        = -2.17339562e0_wp
        h2a5        = 9.54597998e0_wp
        h2a6        = -23.3366749e0_wp
        h2a7        = 24.7115748e0_wp
        h2a8        = -12.8786102e0_wp
        h2a9        = 2.75711758e0_wp
        h2b1        = 1.35164094e0_wp
        h2b2        = -0.254404564e0_wp
        h2b3        = -10.6627055e0_wp
        h2b4        = -14.0112406e0_wp
        h2b5        = -7.80660398e0_wp
        h2b6        = 5.47697753e0_wp
        h2b7        = -6.95875664e0_wp
        h2b8        = 3.64334771e0_wp
        h2b9        = -0.492505596e0_wp
        h2i1        = 160.484e0_wp
        h2i2        = -20.397e0_wp
        h2i3        = 3.913e0_wp
        h2i4        = 0.777e0_wp
        h2i5        = -7.583
        h2i6        = 1.501
        rnaie       = 5.139e0_wp
        rnapol      = 159.267e0_wp
        hpol        = 4.5e0_wp
        hie         = 13.598e0_wp
        hnai1       = -193.341e0_wp
        hnai2       = 16.968e0_wp
        hnai3       = -77.120e0_wp
        hnai4       = -1.024e0_wp
        hnai5       = 6.009e0_wp
        hnai6       = -1.241e0_wp
        hnau1       = 119.880e0_wp
        hnau2       = -0.620e0_wp
        hnau3       = 31.335e0_wp
        hnau4       = 70.414e0_wp
        hnau5       = -13.016e0_wp
        hnau6       = 0.8681e0_wp
        hnaa1       = 0.67e0_wp
        hnaa2       = -92.229e0_wp
        hnaa3       = 51.01019e0_wp
        hnaa4       = -130.37826e0_wp
        hnaa5       = 6.40285e0_wp
        hnab1       = 0.34333e0_wp
        hnab2       = -92.229e0_wp
        hnab3       = -48.54867e0_wp
        hnab4       = -17.30523e0_wp
        hnab5       = 0.99821e0_wp
        balp        = 2.07234e0_wp
        brho        = 2.0e0_wp
        rb0         = 5.1998697814229e7_wp
        rb90        = 5.0e8_wp
        f90         = 1.0e0_wp
        f0          = 1.0e0_wp
        rd          = 9.1538750508008e8_wp
        dalp        = 1.43307e0_wp
        drho        = 5.71654e0_wp
        c2a         = 0.67801e0_wp
        c2b         = 0.83900e0_wp
        c2c         = 0.15039e0_wp
        deh2        = 4.74772265651191638e0_wp
        g_ehfid     = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        rac2   = 1.e0_wp / sqrt(2.e0_wp)
        cevau  = 1.e0_wp / 27.2113961e0_wp
        onethd = 1.0e0_wp / 3.0e0_wp
        twothd = 2.0e0_wp / 3.0e0_wp
        do i = 1, nt
            ! useful quantities calculated first, r1^2, R2^2, R3^2, atom to diatom
            ! center of mass distances, the H-Na-H bond angle and the Jacobian
            ! angle chi.
            d2_b = r(i,1) + r(i,1)
            do g_i_ = 1, g_p_
                g_r1s(g_i_) = d2_b * g_r(g_i_,i,1)
            enddo
            r1s  = r(i,1) * r(i,1)
            d2_b = r(i,2) + r(i,2)
            do g_i_ = 1, g_p_
                g_r2s(g_i_) = d2_b * g_r(g_i_,i,2)
            enddo
            r2s  = r(i,2) * r(i,2)
            d2_b = r(i,3) + r(i,3)
            do g_i_ = 1, g_p_
                g_r3s(g_i_) = d2_b * g_r(g_i_,i,3)
            enddo
            r3s  = r(i,3) * r(i,3)
            rnag = sqrt(abs(0.5e0_wp*(r1s+r3s-0.5e0_wp*r2s)))
            if( rnag .gt. 0.0e0_wp )then
                g_rnag(1) = 0.5e0_wp*r(i,1)/rnag
                g_rnag(2) = -0.25e0_wp*r(i,2)/rnag
                g_rnag(3) = 0.5e0_wp*r(i,3)/rnag
            else
                g_rnag(1) = 0.0e0_wp
                g_rnag(2) = 0.0e0_wp
                g_rnag(3) = 0.0e0_wp
            endif
            d2_b = rnag + rnag
            do g_i_ = 1, g_p_
                g_rnags(g_i_) = d2_b * g_rnag(g_i_)
            enddo
            rnags = rnag * rnag
            do g_i_ = 1, g_p_
                g_rnagc(g_i_) = rnags * g_rnag(g_i_) + rnag * g_rnags(g_i_)
            enddo
            rnagc = rnags * rnag
            frac  = 0.958e0_wp
            ! H to NaH center of masses
            rnag1 = sqrt(abs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))
            if( rnag1 .gt. 0.0e0_wp )then
                g_rnag1(1) = r(i,1)*frac*(frac-1.0e0_wp)/rnag1
                g_rnag1(2) = r(i,2)*(1.0e0_wp-frac)/rnag1
                g_rnag1(3) = r(i,3)*frac/rnag1
            else
                g_rnag1(1) = 0.0e0_wp
                g_rnag1(2) = 0.0e0_wp
                g_rnag1(3) = 0.0e0_wp
            endif
            rnag3 = sqrt(abs(r2s+frac*(r1s-r2s-r3s)+r3s*frac*frac))
            if( rnag3 .gt. 0.0e0_wp )then
                g_rnag3(1) = r(i,1)*frac/rnag3
                g_rnag3(2) = r(i,2)*(1.0e0_wp-frac)/rnag3
                g_rnag3(3) = r(i,3)*frac*(frac-1.0e0_wp)/rnag3
            else
                g_rnag3(1) = 0.0e0_wp
                g_rnag3(2) = 0.0e0_wp
                g_rnag3(3) = 0.0e0_wp
            endif
            d2_v = rnag1 * rnag1
            d3_b = d2_v + rnag1 * rnag1 + rnag1 * rnag1
            do g_i_ = 1, g_p_
                g_rnag1c(g_i_) = d3_b * g_rnag1(g_i_)
            enddo
            rnag1c = d2_v * rnag1
            d2_v   = rnag3 * rnag3
            d3_b   = d2_v + rnag3 * rnag3 + rnag3 * rnag3
            do g_i_ = 1, g_p_
                g_rnag3c(g_i_) = d3_b * g_rnag3(g_i_)
            enddo
            rnag3c = d2_v * rnag3
            ! the shortest of any two distances . . .
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-0.5e0_wp) * g_r(g_i_,i,2) + 0.5e0_wp * g_r(g_i_,i,3)
            enddo
            d1_w = 0.5e0_wp * (r(i,3) - r(i,2))
            d3_v = r(i,2) - r(i,3)
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d7_v = 0.5e0_wp + 0.5e0_wp * d5_v
            d8_b = d3_v * 0.5e0_wp * d1_p
            d2_b = 1.0e0_wp + (-d7_v)
            do g_i_ = 1, g_p_
                g_rsm1(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_r(g_i_,i,2) + d2_b * g_r(g_i_,i,3)
            enddo
            rsm1 = r(i,3) + d3_v * d7_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-0.5e0_wp) * g_r(g_i_,i,1) + 0.5e0_wp * g_r(g_i_,i,3)
            enddo
            d1_w = 0.5e0_wp * (r(i,3) - r(i,1))
            d3_v = r(i,1) - r(i,3)
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d7_v = 0.5e0_wp + 0.5e0_wp * d5_v
            d8_b = d3_v * 0.5e0_wp * d1_p
            d2_b = 1.0e0_wp + (-d7_v)
            do g_i_ = 1, g_p_
                g_rsm2(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_r(g_i_,i,1) + d2_b * g_r(g_i_,i,3)
            enddo
            rsm2 = r(i,3) + d3_v * d7_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-0.5e0_wp) * g_r(g_i_,i,2) + 0.5e0_wp * g_r(g_i_,i,1)
            enddo
            d1_w = 0.5e0_wp * (r(i,1) - r(i,2))
            d3_v = r(i,2) - r(i,1)
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d7_v = 0.5e0_wp + 0.5e0_wp * d5_v
            d8_b = d3_v * 0.5e0_wp * d1_p
            d2_b = 1.0e0_wp + (-d7_v)
            do g_i_ = 1, g_p_
                g_rsm3(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_r(g_i_,i,2) + d2_b * g_r(g_i_,i,1)
            enddo
            rsm3  = r(i,1) + d3_v * d7_v
            rexp1 = exp(-r(i,1))
            rexp2 = exp(-r(i,2))
            rexp3 = exp(-r(i,3))
            ! checks on angles . . .
            if( rnag .eq. 0.0e0_wp )then
                do g_i_ = 1, g_p_
                    g_cstsq(g_i_) = 0.0e0_wp
                enddo
                cstsq = 1.0e0_wp
            else
                d7_v = r(i,2) * rnag
                d8_v = 0.5e0_wp * (r3s - r1s) / d7_v
                d9_v = d8_v * d8_v
                d1_p = 2.0e0_wp * d8_v
                d4_b = d1_p * ((-d8_v) / d7_v)
                d5_b = d4_b * rnag
                d6_b = d4_b * r(i,2)
                d7_b = d1_p * (1.0e0_wp / d7_v) * 0.5e0_wp
                do g_i_ = 1, g_p_
                    g_cstsq(g_i_) = d6_b * g_rnag(g_i_) + d5_b * g_r(g_i_,i,2) + (-d7_b) * g_r1s(g_i_) + d7_b * &
                        & g_r3s(g_i_)
                enddo
                cstsq = d9_v
            endif
            d7_v  = real(2, wp) * r(i,1)
            d9_v  = d7_v * r(i,3)
            d10_v = (r1s + r3s - r2s) / d9_v
            d2_b  = 1.0e0_wp / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,3) * real(2, wp)
            do g_i_ = 1, g_p_
                g_cosa(g_i_) = d5_b * g_r(g_i_,i,3) + d6_b * g_r(g_i_,i,1) + (-d2_b) * g_r2s(g_i_) + d2_b * &
                    & g_r3s(g_i_) + d2_b * g_r1s(g_i_)
            enddo
            cosa = d10_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-10.0e0_wp) * g_rnags(g_i_)
            enddo
            d1_w = (-10.0e0_wp) * rnags
            d4_v = 0.5e0_wp + 0.5e0_wp * cosa
            d5_v = cstsq * d4_v
            d7_v = exp(d1_w)
            d1_p = d7_v
            d8_v = 1.0e0_wp - d7_v
            d5_b = (-d5_v) * d1_p
            d6_b = d8_v * d4_v
            d9_b = d8_v * cstsq * 0.5e0_wp
            do g_i_ = 1, g_p_
                g_coscos(g_i_) = d5_b * g_d1_w(g_i_) + d9_b * g_cosa(g_i_) + d6_b * g_cstsq(g_i_)
            enddo
            coscos = d5_v * d8_v
            if( rnag1 .eq. 0.0e0_wp )then
                do g_i_ = 1, g_p_
                    g_cstsq1(g_i_) = 0.0e0_wp
                enddo
                cstsq1 = 1.0e0_wp
            else
                d11_v = r(i,1) * rnag1
                d12_v = (r3s - r2s - r1s + 2.0e0_wp * r3s * frac) / d11_v
                d13_v = d12_v * d12_v
                d1_p  = 2.0e0_wp * d12_v
                d3_b  = 0.25e0_wp * d1_p
                d4_b  = d3_b * (1.0e0_wp / d11_v)
                d5_b  = d3_b * ((-d12_v) / d11_v)
                d6_b  = d5_b * rnag1
                d7_b  = d5_b * r(i,1)
                d11_b = d4_b * frac * 2.0e0_wp + d4_b
                do g_i_ = 1, g_p_
                    g_cstsq1(g_i_) = d7_b * g_rnag1(g_i_) + d6_b * g_r(g_i_,i,1) + (-d4_b) * g_r1s(g_i_) + (- d4_b) * &
                        & g_r2s(g_i_) + d11_b * g_r3s(g_i_)
                enddo
                cstsq1 = 0.25e0_wp * d13_v
            endif
            if( rnag3 .eq. 0.0e0_wp )then
                do g_i_ = 1, g_p_
                    g_cstsq3(g_i_) = 0.0e0_wp
                enddo
                cstsq3 = 1.0e0_wp
            else
                d11_v = r(i,3) * rnag3
                d12_v = (r1s - r2s - r3s + 2.0e0_wp * r1s * frac) / d11_v
                d13_v = d12_v * d12_v
                d1_p  = 2.0e0_wp * d12_v
                d3_b  = 0.25e0_wp * d1_p
                d4_b  = d3_b * (1.0e0_wp / d11_v)
                d5_b  = d3_b * ((-d12_v) / d11_v)
                d6_b  = d5_b * rnag3
                d7_b  = d5_b * r(i,3)
                d11_b = d4_b * frac * 2.0e0_wp + d4_b
                do g_i_ = 1, g_p_
                    g_cstsq3(g_i_) = d7_b * g_rnag3(g_i_) + d6_b * g_r(g_i_,i,3) + (-d4_b) * g_r3s(g_i_) + (- d4_b) * &
                        & g_r2s(g_i_) + d11_b * g_r1s(g_i_)
                enddo
                cstsq3 = 0.25e0_wp * d13_v
            endif
            !==============================================================k
            !This section contains the diatomic curves
            ! The first NaH singlet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s1c2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-s1c2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s1c6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-s1c6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s1c3 + s1c4 * r(i,1) + s1c5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * s1c5
            d7_b  = d12_v * d9_v + d8_b * s1c4
            d14_b = s1c1 * d2_p
            do g_i_ = 1, g_p_
                g_s1a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            s1a = s1c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s10c2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-s10c2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s10c6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-s10c6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s10c3 + s10c4 * r(i,1) + s10c5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * s10c5
            d7_b  = d12_v * d9_v + d8_b * s10c4
            d14_b = s10c1 * d2_p
            do g_i_ = 1, g_p_
                g_s10c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            s10c = s10c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s10alp * g_rsm1(g_i_)
            enddo
            d1_w = s10alp * (rsm1 - s10rho)
            d4_v = (s10c - s1a) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.0e0_wp - d6_v
            d7_b = (-d4_v) * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_s10(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s10c(g_i_) + d2_b * g_s1a(g_i_)
            enddo
            s10 = s1a + d4_v * d7_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s102i) * g_r(g_i_,i,1)
            enddo
            d1_w = (-s102i) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s106i) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-s106i) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s103i + s104i * r(i,1) + s105i * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * s105i
            d7_b  = d12_v * d9_v + d8_b * s104i
            d14_b = s101i * d2_p
            do g_i_ = 1, g_p_
                g_s1ic(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            s1ic = s101i * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s1ialp * g_rnags(g_i_)
            enddo
            d1_w = s1ialp * rnags
            d3_v = s1ic - s1a
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d6_v = 1.0e0_wp - d5_v
            d7_b = (-d3_v) * d1_p
            d2_b = 1.0e0_wp + (-d6_v)
            do g_i_ = 1, g_p_
                g_s1i(g_i_) = d7_b * g_d1_w(g_i_) + d6_v * g_s1ic(g_i_) + d2_b * g_s1a(g_i_)
            enddo
            s1i = s1a + d3_v * d6_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s190c2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-s190c2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s190c6) * g_r(g_i_,i,1)
            enddo
            d2_w = (-s190c6) * r(i,1)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d5_v = r(i,1) * s190c3
            d7_v = exp(d2_w)
            d1_p = d7_v
            d6_b = d5_v * d1_p
            d7_b = d7_v * s190c3
            d9_b = s190c1 * d2_p
            do g_i_ = 1, g_p_
                g_s190c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r(g_i_,i,1) + d9_b * g_d1_w(g_i_)
            enddo
            s190c = s190c1 * d2_v + d5_v * d7_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s190alp * g_r(g_i_,i,2)
            enddo
            d1_w = s190alp * (r(i,2) - s190rho)
            d4_v = (s190c - s1a) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.0e0_wp - d6_v
            d7_b = (-d4_v) * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_s190(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s190c(g_i_) + d2_b * g_s1a(g_i_)
            enddo
            s190  = s1a + d4_v * d7_v
            d3_v  = s10 - s190
            d8_v  = 0.5e0_wp + 0.5e0_wp * cosa
            d13_v = 0.5e0_wp - 0.5e0_wp * cosa
            d7_b  = (-s1i) * 0.5e0_wp + s190 * 0.5e0_wp
            d10_b = d8_v + (-coscos)
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d13_v * g_s1i(g_i_) + d7_b * g_cosa(g_i_) + d3_v * g_coscos(g_i_) + d10_b * &
                    & g_s190(g_i_) + coscos * g_s10(g_i_)
            enddo
            s1 = d3_v * coscos + s190 * d8_v + s1i * d13_v
            ! The first NaH triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t190c2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-t190c2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t190c6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-t190c6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t190c3 + t190c4 * r(i,1) + t190c5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * t190c5
            d7_b  = d12_v * d9_v + d8_b * t190c4
            d14_b = t190c1 * d2_p
            do g_i_ = 1, g_p_
                g_t190(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            t190 = t190c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t10c2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-t10c2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t10c6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-t10c6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t10c3 + t10c4 * r(i,1) + t10c5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * t10c5
            d7_b  = d12_v * d9_v + d8_b * t10c4
            d14_b = t10c1 * d2_p
            do g_i_ = 1, g_p_
                g_t10(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            t10 = t10c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t10ac2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-t10ac2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t10ac6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-t10ac6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t10ac3 + t10ac4 * r(i,1) + t10ac5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * t10ac5
            d7_b  = d12_v * d9_v + d8_b * t10ac4
            d14_b = t10ac1 * d2_p
            do g_i_ = 1, g_p_
                g_t10alt(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            t10alt = t10ac1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t10alp * g_r(g_i_,i,2)
            enddo
            d1_w = t10alp * (r(i,2) - t10rho)
            d3_v = t10alt - t10
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d7_v = 0.5e0_wp + 0.5e0_wp * d5_v
            d8_b = d3_v * 0.5e0_wp * d1_p
            d2_b = 1.0e0_wp + (-d7_v)
            do g_i_ = 1, g_p_
                g_t10(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_t10alt(g_i_) + d2_b * g_t10(g_i_)
            enddo
            t10 = t10 + d3_v * d7_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t1ic2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-t1ic2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t1ic6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-t1ic6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t1ic3 + t1ic4 * r(i,1) + t1ic5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * t1ic5
            d7_b  = d12_v * d9_v + d8_b * t1ic4
            d14_b = t1ic1 * d2_p
            do g_i_ = 1, g_p_
                g_t1i(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            t1i   = t1ic1 * d2_v + d10_v * d12_v
            d3_v  = t10 - t190
            d8_v  = 0.5e0_wp + 0.5e0_wp * cosa
            d13_v = 0.5e0_wp - 0.5e0_wp * cosa
            d7_b  = (-t1i) * 0.5e0_wp + t190 * 0.5e0_wp
            d10_b = d8_v + (-coscos)
            do g_i_ = 1, g_p_
                g_t1(g_i_) = d13_v * g_t1i(g_i_) + d7_b * g_cosa(g_i_) + d3_v * g_coscos(g_i_) + d10_b * &
                    & g_t190(g_i_) + coscos * g_t10(g_i_)
            enddo
            t1 = d3_v * coscos + t190 * d8_v + t1i * d13_v
            ! The H2 singlet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s2c3) * g_r(g_i_,i,2)
            enddo
            d1_w = (-s2c3) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s2c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-s2c6) * r(i,2)
            d3_v  = s2c1 + s2c2 * r(i,2)
            d5_v  = exp(d1_w)
            d2_p  = d5_v
            d9_v  = s2c4 + s2c5 * r(i,2)
            d10_v = r2s * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d7_b  = d12_v * d9_v
            d13_b = d3_v * d2_p
            d10_b = d12_v * r2s * s2c5 + d5_v * s2c2
            do g_i_ = 1, g_p_
                g_s2a(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r2s(g_i_) + d13_b * g_d1_w(g_i_) + d10_b * g_r(g_i_,i,2)
            enddo
            s2a = d3_v * d5_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s290c3) * g_r(g_i_,i,2)
            enddo
            d1_w = (-s290c3) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s290c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-s290c6) * r(i,2)
            d3_v  = s290c1 + s290c2 * r(i,2)
            d5_v  = exp(d1_w)
            d2_p  = d5_v
            d9_v  = s290c4 + s290c5 * r(i,2)
            d10_v = r2s * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d7_b  = d12_v * d9_v
            d13_b = d3_v * d2_p
            d10_b = d12_v * r2s * s290c5 + d5_v * s290c2
            do g_i_ = 1, g_p_
                g_s290c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r2s(g_i_) + d13_b * g_d1_w(g_i_) + d10_b * g_r(g_i_,i,2)
            enddo
            s290c = d3_v * d5_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s290alp * g_rnags(g_i_)
            enddo
            d1_w = s290alp * (rnags - s290rho)
            d4_v = (s290c - s2a) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.e0_wp - d6_v
            d7_b = (-d4_v) * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_s290(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s290c(g_i_) + d2_b * g_s2a(g_i_)
            enddo
            s290 = s2a + d4_v * d7_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s20c3) * g_r(g_i_,i,2)
            enddo
            d1_w = (-s20c3) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s20c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-s20c6) * r(i,2)
            d3_v  = s20c1 + s20c2 * r(i,2)
            d5_v  = exp(d1_w)
            d2_p  = d5_v
            d9_v  = s20c4 + s20c5 * r(i,2)
            d10_v = r2s * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d7_b  = d12_v * d9_v
            d13_b = d3_v * d2_p
            d10_b = d12_v * r2s * s20c5 + d5_v * s20c2
            do g_i_ = 1, g_p_
                g_s20c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r2s(g_i_) + d13_b * g_d1_w(g_i_) + d10_b * g_r(g_i_,i,2)
            enddo
            s20c = d3_v * d5_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s20alp * g_rsm2(g_i_)
            enddo
            d1_w = s20alp * (rsm2 - s20rho)
            d4_v = (s20c - s2a) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.e0_wp - d6_v
            d7_b = (-d4_v) * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_s20(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s20c(g_i_) + d2_b * g_s2a(g_i_)
            enddo
            s20  = s2a + d4_v * d7_v
            d3_v = s20 - s290
            d2_b = 1.0e0_wp + (-coscos)
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d3_v * g_coscos(g_i_) + coscos * g_s20(g_i_) + d2_b * g_s290(g_i_)
            enddo
            s2 = s290 + d3_v * coscos
            ! The H2 triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t2c2) * g_r(g_i_,i,2)
            enddo
            d1_w = (-t2c2) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t2c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-t2c6) * r(i,2)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t2c3 + t2c4 * r(i,2) + t2c5 * r2s
            d10_v = r(i,2) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,2)
            d11_b = d8_b * t2c5
            d7_b  = d12_v * d9_v + d8_b * t2c4
            d14_b = t2c1 * d2_p
            do g_i_ = 1, g_p_
                g_t2a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + d7_b * g_r(g_i_,i,2) + d14_b * g_d1_w(g_i_)
            enddo
            t2a = t2c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t20c2) * g_r(g_i_,i,2)
            enddo
            d1_w = (-t20c2) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t20c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-t20c6) * r(i,2)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t20c3 + t20c4 * r(i,2) + t20c5 * r2s
            d10_v = r(i,2) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,2)
            d11_b = d8_b * t20c5
            d7_b  = d12_v * d9_v + d8_b * t20c4
            d14_b = t20c1 * d2_p
            do g_i_ = 1, g_p_
                g_t20c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + d7_b * g_r(g_i_,i,2) + d14_b * g_d1_w(g_i_)
            enddo
            t20c = t20c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t20alp * g_rnags(g_i_)
            enddo
            d1_w = t20alp * (rnags - t20rho)
            d4_v = (t20c - t2a) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.e0_wp - d6_v
            d7_b = (-d4_v) * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_t20(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_t20c(g_i_) + d2_b * g_t2a(g_i_)
            enddo
            t20  = t2a + d4_v * d7_v
            d3_v = t20 - t2a
            d2_b = 1.0e0_wp + (-coscos)
            do g_i_ = 1, g_p_
                g_t2(g_i_) = d3_v * g_coscos(g_i_) + coscos * g_t20(g_i_) + d2_b * g_t2a(g_i_)
            enddo
            t2 = t2a + d3_v * coscos
            ! The second NaH singlet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s1c2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-s1c2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s1c6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-s1c6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s1c3 + s1c4 * r(i,3) + s1c5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * s1c5
            d7_b  = d12_v * d9_v + d8_b * s1c4
            d14_b = s1c1 * d2_p
            do g_i_ = 1, g_p_
                g_s3a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            s3a = s1c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s10c2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-s10c2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s10c6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-s10c6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s10c3 + s10c4 * r(i,3) + s10c5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * s10c5
            d7_b  = d12_v * d9_v + d8_b * s10c4
            d14_b = s10c1 * d2_p
            do g_i_ = 1, g_p_
                g_s30c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            s30c = s10c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s10alp * g_rsm3(g_i_)
            enddo
            d1_w = s10alp * (rsm3 - s10rho)
            d4_v = (s30c - s3a) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.0e0_wp - d6_v
            d7_b = (-d4_v) * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_s30(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s30c(g_i_) + d2_b * g_s3a(g_i_)
            enddo
            s30 = s3a + d4_v * d7_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s102i) * g_r(g_i_,i,3)
            enddo
            d1_w = (-s102i) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s106i) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-s106i) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s103i + s104i * r(i,3) + s105i * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * s105i
            d7_b  = d12_v * d9_v + d8_b * s104i
            d14_b = s101i * d2_p
            do g_i_ = 1, g_p_
                g_s3ic(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            s3ic = s101i * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s1ialp * g_rnags(g_i_)
            enddo
            d1_w = s1ialp * rnags
            d3_v = s3ic - s3a
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d6_v = 1.0e0_wp - d5_v
            d7_b = (-d3_v) * d1_p
            d2_b = 1.0e0_wp + (-d6_v)
            do g_i_ = 1, g_p_
                g_s3i(g_i_) = d7_b * g_d1_w(g_i_) + d6_v * g_s3ic(g_i_) + d2_b * g_s3a(g_i_)
            enddo
            s3i = s3a + d3_v * d6_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s190c2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-s190c2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s190c6) * g_r(g_i_,i,3)
            enddo
            d2_w = (-s190c6) * r(i,3)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d5_v = r(i,3) * s190c3
            d7_v = exp(d2_w)
            d1_p = d7_v
            d6_b = d5_v * d1_p
            d7_b = d7_v * s190c3
            d9_b = s190c1 * d2_p
            do g_i_ = 1, g_p_
                g_s390c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r(g_i_,i,3) + d9_b * g_d1_w(g_i_)
            enddo
            s390c = s190c1 * d2_v + d5_v * d7_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s190alp * g_r(g_i_,i,2)
            enddo
            d1_w = s190alp * (r(i,2) - s190rho)
            d4_v = (s390c - s3a) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.0e0_wp - d6_v
            d7_b = (-d4_v) * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_s390(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s390c(g_i_) + d2_b * g_s3a(g_i_)
            enddo
            s390  = s3a + d4_v * d7_v
            d3_v  = s30 - s390
            d8_v  = 0.5e0_wp + 0.5e0_wp * cosa
            d13_v = 0.5e0_wp - 0.5e0_wp * cosa
            d7_b  = (-s3i) * 0.5e0_wp + s390 * 0.5e0_wp
            d10_b = d8_v + (-coscos)
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d13_v * g_s3i(g_i_) + d7_b * g_cosa(g_i_) + d3_v * g_coscos(g_i_) + d10_b * &
                    & g_s390(g_i_) + coscos * g_s30(g_i_)
            enddo
            s3 = d3_v * coscos + s390 * d8_v + s3i * d13_v
            ! The second NaH triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t190c2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-t190c2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t190c6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-t190c6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t190c3 + t190c4 * r(i,3) + t190c5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * t190c5
            d7_b  = d12_v * d9_v + d8_b * t190c4
            d14_b = t190c1 * d2_p
            do g_i_ = 1, g_p_
                g_t390(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            t390 = t190c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t10c2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-t10c2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t10c6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-t10c6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t10c3 + t10c4 * r(i,3) + t10c5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * t10c5
            d7_b  = d12_v * d9_v + d8_b * t10c4
            d14_b = t10c1 * d2_p
            do g_i_ = 1, g_p_
                g_t30(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            t30 = t10c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t10ac2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-t10ac2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t10ac6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-t10ac6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t10ac3 + t10ac4 * r(i,3) + t10ac5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * t10ac5
            d7_b  = d12_v * d9_v + d8_b * t10ac4
            d14_b = t10ac1 * d2_p
            do g_i_ = 1, g_p_
                g_t30alt(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            t30alt = t10ac1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t10alp * g_r(g_i_,i,2)
            enddo
            d1_w = t10alp * (r(i,2) - t10rho)
            d3_v = t30alt - t30
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d7_v = 0.5e0_wp + 0.5e0_wp * d5_v
            d8_b = d3_v * 0.5e0_wp * d1_p
            d2_b = 1.0e0_wp + (-d7_v)
            do g_i_ = 1, g_p_
                g_t30(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_t30alt(g_i_) + d2_b * g_t30(g_i_)
            enddo
            t30  = t30 + d3_v * d7_v
            d3_v = t30 - t390
            d2_b = 1.0e0_wp + (-cstsq)
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d3_v * g_cstsq(g_i_) + cstsq * g_t30(g_i_) + d2_b * g_t390(g_i_)
            enddo
            t3 = t390 + d3_v * cstsq
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t1ic2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-t1ic2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t1ic6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-t1ic6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t1ic3 + t1ic4 * r(i,3) + t1ic5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * t1ic5
            d7_b  = d12_v * d9_v + d8_b * t1ic4
            d14_b = t1ic1 * d2_p
            do g_i_ = 1, g_p_
                g_t3i(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            t3i  = t1ic1 * d2_v + d10_v * d12_v
            d4_v = 0.5e0_wp * (1.0e0_wp - cosa)
            d6_v = t3i - t3
            d2_b = 1.0e0_wp + (-d4_v)
            d8_b = -(d6_v * 0.5e0_wp)
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d4_v * g_t3i(g_i_) + d8_b * g_cosa(g_i_) + d2_b * g_t3(g_i_)
            enddo
            t3    = t3 + d4_v * d6_v
            d3_v  = t30 - t390
            d8_v  = 0.5e0_wp + 0.5e0_wp * cosa
            d13_v = 0.5e0_wp - 0.5e0_wp * cosa
            d7_b  = (-t3i) * 0.5e0_wp + t390 * 0.5e0_wp
            d10_b = d8_v + (-coscos)
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d13_v * g_t3i(g_i_) + d7_b * g_cosa(g_i_) + d3_v * g_coscos(g_i_) + d10_b * &
                    & g_t390(g_i_) + coscos * g_t30(g_i_)
            enddo
            t3 = d3_v * coscos + t390 * d8_v + t3i * d13_v
            ! this long section contains data pertinent to long range forces
            !  H2 ionization energy:
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-h2i3) * g_r(g_i_,i,2)
            enddo
            d1_w = (-h2i3) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-h2i6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-h2i6) * r(i,2)
            d3_v  = h2i1 + h2i2 * r(i,2)
            d5_v  = exp(d1_w)
            d2_p  = d5_v
            d9_v  = h2i4 + h2i5 * r(i,2)
            d10_v = r2s * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d9_b  = d10_v * d1_p
            d10_b = d12_v * d9_v
            d16_b = d3_v * d2_p
            d13_b = d12_v * r2s * h2i5 + d5_v * h2i2
            do g_i_ = 1, g_p_
                g_h2ie(g_i_) = -g_s2a(g_i_) + d9_b * g_d2_w(g_i_) + d10_b * g_r2s(g_i_) + d16_b * g_d1_w(g_i_) + &
                    & d13_b * g_r(g_i_,i,2)
            enddo
            h2ie = d3_v * d5_v + d10_v * d12_v - s2a + hie
            !  H2 polarizability:
            !   al1 is the parallel polarizability for h2
            !   al2 is the perpendicaular polarizability for h2
            !   alt is the anisotropic? polarizability
            !   alb is the mean polarizability
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-h2a1) * g_r(g_i_,i,2)
            enddo
            d1_w  = h2a1 * (h2a2 - r(i,2))
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = h2a6 * r(i,2)
            d12_v = h2a7 * r2s
            d15_v = h2a8 * r2s
            d16_v = d15_v * r2s
            d19_v = h2a9 * r2s
            d20_v = d19_v * r2s
            d22_v = h2a3 + h2a4 * r(i,2) + h2a5 * r2s + d9_v * r2s + d12_v * r2s + d16_v * r(i,2) + d20_v * r2s
            d7_b  = d2_v * r2s
            d12_b = d2_v * r(i,2)
            d8_b = d2_v * d20_v + d7_b * d19_v + d7_b * r2s * h2a9 + d12_b * d15_v + d12_b * r2s * h2a8 + d2_v * &
                & d12_v + d2_v * r2s * h2a7 + d2_v * d9_v + d2_v * h2a5
            d13_b = d2_v * d16_v + d2_v * r2s * h2a6 + d2_v * h2a4
            d24_b = d22_v * d1_p
            do g_i_ = 1, g_p_
                g_al1(g_i_) = d8_b * g_r2s(g_i_) + d13_b * g_r(g_i_,i,2) + d24_b * g_d1_w(g_i_)
            enddo
            al1 = d2_v * d22_v + 2.0e0_wp * hpol
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-h2b1) * g_r(g_i_,i,2)
            enddo
            d1_w  = h2b1 * (h2b2 - r(i,2))
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = h2b6 * r(i,2)
            d12_v = h2b7 * r2s
            d15_v = h2b8 * r2s
            d16_v = d15_v * r2s
            d19_v = h2b9 * r2s
            d20_v = d19_v * r2s
            d22_v = h2b3 + h2b4 * r(i,2) + h2b5 * r2s + d9_v * r2s + d12_v * r2s + d16_v * r(i,2) + d20_v * r2s
            d7_b  = d2_v * r2s
            d12_b = d2_v * r(i,2)
            d8_b = d2_v * d20_v + d7_b * d19_v + d7_b * r2s * h2b9 + d12_b * d15_v + d12_b * r2s * h2b8 + d2_v * &
                & d12_v + d2_v * r2s * h2b7 + d2_v * d9_v + d2_v * h2b5
            d13_b = d2_v * d16_v + d2_v * r2s * h2b6 + d2_v * h2b4
            d24_b = d22_v * d1_p
            do g_i_ = 1, g_p_
                g_al2(g_i_) = d8_b * g_r2s(g_i_) + d13_b * g_r(g_i_,i,2) + d24_b * g_d1_w(g_i_)
            enddo
            al2  = d2_v * d22_v + 2.0e0_wp * hpol
            alt  = twothd * (al1 - al2)
            alb  = onethd * al1 + twothd * al2
            alhh = alb + alt * (1.5e0_wp * cstsq - 0.5e0_wp)
            !  The NaH + H dispersion and NaH + H induction terms must go on the low
            !t,
            !  so that the adiabat is fitted correctly in the region of the coupling
            ! NaH properties
            d9_v  = hnai6 * r(i,1)
            d11_v = hnai3 + hnai4 * r(i,1) + hnai5 * r1s + d9_v * r1s
            d12_v = (hnai1 + hnai2 * r(i,1)) / d11_v
            d4_b  = (-d12_v) / d11_v
            d8_b  = d4_b * d9_v + d4_b * hnai5
            d9_b  = d4_b * r1s * hnai6 + d4_b * hnai4 + 1.0e0_wp / d11_v * hnai2
            do g_i_ = 1, g_p_
                g_rnahie1(g_i_) = d8_b * g_r1s(g_i_) + d9_b * g_r(g_i_,i,1 )
            enddo
            rnahie1 = d12_v + rnaie
            d9_v    = hnai6 * r(i,3)
            d11_v   = hnai3 + hnai4 * r(i,3) + hnai5 * r3s + d9_v * r3s
            d12_v   = (hnai1 + hnai2 * r(i,3)) / d11_v
            d4_b    = (-d12_v) / d11_v
            d8_b    = d4_b * d9_v + d4_b * hnai5
            d9_b    = d4_b * r3s * hnai6 + d4_b * hnai4 + 1.0e0_wp / d11_v * hnai2
            do g_i_ = 1, g_p_
                g_rnahie3(g_i_) = d8_b * g_r3s(g_i_) + d9_b * g_r(g_i_,i,3 )
            enddo
            rnahie3 = d12_v + rnaie
            d10_v   = hnau6 * r(i,1)
            d12_v   = hnau3 + hnau4 * r(i,1) + hnau5 * r1s + d10_v * r1s
            d13_v   = (hnau1 * r(i,1) + hnau2 * r1s) / d12_v
            d2_b    = 1.0e0_wp / d12_v
            d3_b    = (-d13_v) / d12_v
            d7_b    = d3_b * d10_v + d3_b * hnau5 + d2_b * hnau2
            d8_b    = d3_b * r1s * hnau6 + d3_b * hnau4 + d2_b * hnau1
            do g_i_ = 1, g_p_
                g_rnahmu1(g_i_) = d7_b * g_r1s(g_i_) + d8_b * g_r(g_i_,i,1 )
            enddo
            rnahmu1 = d13_v
            d10_v   = hnau6 * r(i,3)
            d12_v   = hnau3 + hnau4 * r(i,3) + hnau5 * r3s + d10_v * r3s
            d13_v   = (hnau1 * r(i,3) + hnau2 * r3s) / d12_v
            d2_b    = 1.0e0_wp / d12_v
            d3_b    = (-d13_v) / d12_v
            d7_b    = d3_b * d10_v + d3_b * hnau5 + d2_b * hnau2
            d8_b    = d3_b * r3s * hnau6 + d3_b * hnau4 + d2_b * hnau1
            do g_i_ = 1, g_p_
                g_rnahmu3(g_i_) = d7_b * g_r3s(g_i_) + d8_b * g_r(g_i_,i,3 )
            enddo
            rnahmu3 = d13_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-hnaa1) * g_r(g_i_,i,1)
            enddo
            d1_w  = (-hnaa1) * r(i,1)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = hnaa5 * r1s
            d11_v = hnaa2 + hnaa3 * r(i,1) + hnaa4 * r1s + d9_v * r(i,1)
            d9_b  = d2_v * r(i,1) * hnaa5 + d2_v * hnaa4
            d8_b  = d2_v * d9_v + d2_v * hnaa3
            d13_b = d11_v * d1_p
            do g_i_ = 1, g_p_
                g_rnahpar1(g_i_) = d9_b * g_r1s(g_i_) + d8_b * g_r(g_i_,i,1) + d13_b * g_d1_w(g_i_)
            enddo
            rnahpar1 = rnapol + hpol + d2_v * d11_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-hnaa1) * g_r(g_i_,i,3)
            enddo
            d1_w  = (-hnaa1) * r(i,3)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = hnaa5 * r3s
            d11_v = hnaa2 + hnaa3 * r(i,3) + hnaa4 * r3s + d9_v * r(i,3)
            d9_b  = d2_v * r(i,3) * hnaa5 + d2_v * hnaa4
            d8_b  = d2_v * d9_v + d2_v * hnaa3
            d13_b = d11_v * d1_p
            do g_i_ = 1, g_p_
                g_rnahpar3(g_i_) = d9_b * g_r3s(g_i_) + d8_b * g_r(g_i_,i,3) + d13_b * g_d1_w(g_i_)
            enddo
            rnahpar3 = rnapol + hpol + d2_v * d11_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-hnab1) * g_r(g_i_,i,1)
            enddo
            d1_w  = (-hnab1) * r(i,1)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = hnab5 * r1s
            d11_v = hnab2 + hnab3 * r(i,1) + hnab4 * r1s + d9_v * r(i,1)
            d9_b  = d2_v * r(i,1) * hnab5 + d2_v * hnab4
            d8_b  = d2_v * d9_v + d2_v * hnab3
            d13_b = d11_v * d1_p
            do g_i_ = 1, g_p_
                g_rnahperp1(g_i_) = d9_b * g_r1s(g_i_) + d8_b * g_r(g_i_,i,1) + d13_b * g_d1_w(g_i_)
            enddo
            rnahperp1 = rnapol + hpol + d2_v * d11_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-hnab1) * g_r(g_i_,i,3)
            enddo
            d1_w  = (-hnab1) * r(i,3)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = hnab5 * r3s
            d11_v = hnab2 + hnab3 * r(i,3) + hnab4 * r3s + d9_v * r(i,3)
            d9_b  = d2_v * r(i,3) * hnab5 + d2_v * hnab4
            d8_b  = d2_v * d9_v + d2_v * hnab3
            d13_b = d11_v * d1_p
            do g_i_ = 1, g_p_
                g_rnahperp3(g_i_) = d9_b * g_r3s(g_i_) + d8_b * g_r(g_i_,i,3) + d13_b * g_d1_w(g_i_)
            enddo
            rnahperp3 = rnapol + hpol + d2_v * d11_v
            d3_v      = rnahpar1 - rnahperp1
            d2_b      = 1.0e0_wp + (-cstsq1)
            do g_i_ = 1, g_p_
                g_rnahpol1(g_i_) = d3_v * g_cstsq1(g_i_) + cstsq1 * g_rnahpar1(g_i_) + d2_b * g_rnahperp1(g_i_)
            enddo
            rnahpol1 = rnahperp1 + d3_v * cstsq1
            do g_i_ = 1, g_p_
                g_rnahpol1b(g_i_) = twothd * g_rnahperp1(g_i_) + onethd * g_rnahpar1(g_i_)
            enddo
            rnahpol1b = onethd * rnahpar1 + twothd * rnahperp1
            d3_v      = rnahpar3 - rnahperp3
            d2_b      = 1.0e0_wp + (-cstsq3)
            do g_i_ = 1, g_p_
                g_rnahpol3(g_i_) = d3_v * g_cstsq3(g_i_) + cstsq3 * g_rnahpar3(g_i_) + d2_b * g_rnahperp3(g_i_)
            enddo
            rnahpol3 = rnahperp3 + d3_v * cstsq3
            do g_i_ = 1, g_p_
                g_rnahpol3b(g_i_) = twothd * g_rnahperp3(g_i_) + onethd * g_rnahpar3(g_i_)
            enddo
            rnahpol3b = onethd * rnahpar3 + twothd * rnahperp3
            ! induction forces
            ! dipole induced dipole
            d2_v = rnag1c * rnag1c
            d4_v = d2_v * rnag1c + rd
            d5_v = 1.0e0_wp / d4_v
            d3_b = (-d5_v) / d4_v
            d4_b = d3_b * rnag1c
            d5_b = d3_b * d2_v + d4_b * rnag1c + d4_b * rnag1c
            do g_i_ = 1, g_p_
                g_div1(g_i_) = d5_b * g_rnag1c(g_i_)
            enddo
            div1  = d5_v
            d2_v  = rnahmu1 * rnahmu1
            d1_p  = 2.0e0_wp * rnahmu1
            d4_v  = (-d2_v) * hpol
            d7_v  = 1.5e0_wp * cstsq1 + 0.5e0_wp
            d8_v  = d4_v * d7_v
            d10_v = d8_v * rnag1c
            d4_b  = div1 * rnag1c
            d5_b  = div1 * d8_v
            d9_b  = d4_b * d4_v * 1.5e0_wp
            d12_b = (-(d4_b * d7_v * hpol)) * d1_p
            do g_i_ = 1, g_p_
                g_rmuimu1(g_i_) = d10_v * g_div1(g_i_) + d5_b * g_rnag1c(g_i_) + d9_b * g_cstsq1(g_i_) + d12_b * &
                    & g_rnahmu1(g_i_)
            enddo
            rmuimu1 = d10_v * div1
            d2_v    = rnag3c * rnag3c
            d4_v    = d2_v * rnag3c + rd
            d5_v    = 1.0e0_wp / d4_v
            d3_b    = (-d5_v) / d4_v
            d4_b    = d3_b * rnag3c
            d5_b    = d3_b * d2_v + d4_b * rnag3c + d4_b * rnag3c
            do g_i_ = 1, g_p_
                g_div3(g_i_) = d5_b * g_rnag3c(g_i_)
            enddo
            div3  = d5_v
            d2_v  = rnahmu3 * rnahmu3
            d1_p  = 2.0e0_wp * rnahmu3
            d4_v  = (-d2_v) * hpol
            d7_v  = 1.5e0_wp * cstsq3 + 0.5e0_wp
            d8_v  = d4_v * d7_v
            d10_v = d8_v * rnag3c
            d4_b  = div3 * rnag3c
            d5_b  = div3 * d8_v
            d9_b  = d4_b * d4_v * 1.5e0_wp
            d12_b = (-(d4_b * d7_v * hpol)) * d1_p
            do g_i_ = 1, g_p_
                g_rmuimu3(g_i_) = d10_v * g_div3(g_i_) + d5_b * g_rnag3c(g_i_) + d9_b * g_cstsq3(g_i_) + d12_b * &
                    & g_rnahmu3(g_i_)
            enddo
            rmuimu3 = d10_v * div3
            ! dispersion forces
            d3_v  = (-1.5e0_wp) * rnaie * h2ie * rnapol
            d6_v  = rnaie + h2ie
            d8_v  = rnagc * rnagc
            d10_v = d8_v * rnagc + rb90
            d11_v = d6_v * d10_v
            d12_v = d3_v * al2 / d11_v
            d2_b  = 1.0e0_wp / d11_v
            d3_b  = (-d12_v) / d11_v
            d6_b  = d3_b * d6_v
            d7_b  = d6_b * rnagc
            d8_b  = d6_b * d8_v + d7_b * rnagc + d7_b * rnagc
            d11_b = d2_b * d3_v
            d9_b  = d3_b * d10_v + d2_b * al2 * rnapol * ((-1.5e0_wp) * rnaie)
            do g_i_ = 1, g_p_
                g_c690(g_i_) = d8_b * g_rnagc(g_i_) + d11_b * g_al2(g_i_) + d9_b * g_h2ie(g_i_)
            enddo
            c690 = d12_v
            d2_v = c690 * f90
            d4_b = rnagc * f90
            do g_i_ = 1, g_p_
                g_c690(g_i_) = d2_v * g_rnagc(g_i_) + d4_b * g_c690(g_i_)
            enddo
            c690  = d2_v * rnagc
            d3_v  = (-1.5e0_wp) * rnaie * h2ie * rnapol
            d6_v  = rnaie + h2ie
            d8_v  = rnagc * rnagc
            d10_v = d8_v * rnagc + rb0
            d11_v = d6_v * d10_v
            d12_v = d3_v * al1 / d11_v
            d2_b  = 1.0e0_wp / d11_v
            d3_b  = (-d12_v) / d11_v
            d6_b  = d3_b * d6_v
            d7_b  = d6_b * rnagc
            d8_b  = d6_b * d8_v + d7_b * rnagc + d7_b * rnagc
            d11_b = d2_b * d3_v
            d9_b  = d3_b * d10_v + d2_b * al1 * rnapol * ((-1.5e0_wp) * rnaie)
            do g_i_ = 1, g_p_
                g_c60(g_i_) = d8_b * g_rnagc(g_i_) + d11_b * g_al1(g_i_) + d9_b * g_h2ie(g_i_)
            enddo
            c60  = d12_v
            d2_v = c60 * f0
            d4_b = rnagc * f0
            do g_i_ = 1, g_p_
                g_c60(g_i_) = d2_v * g_rnagc(g_i_) + d4_b * g_c60(g_i_)
            enddo
            c60  = d2_v * rnagc
            d2_v = (-39.0625e0_wp) * rnag
            d4_b = rnagc * (-39.0625e0_wp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d2_v * g_rnagc(g_i_) + d4_b * g_rnag(g_i_)
            enddo
            d1_w = d2_v * rnagc
            d3_v = c60 - c690
            d5_v = d3_v * cstsq
            d7_v = exp(d1_w)
            d1_p = d7_v
            d8_v = 1.0e0_wp - d7_v
            d7_b = (-d5_v) * d1_p
            d8_b = d8_v * cstsq
            d9_b = d8_v * d3_v
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_c6(g_i_) = d7_b * g_d1_w(g_i_) + d9_b * g_cstsq(g_i_) + d8_b * g_c60(g_i_) + d2_b * g_c690(g_i_)
            enddo
            c6    = c690 + d5_v * d8_v
            d3_v  = (-1.5e0_wp) * rnahie1 * hie
            d6_v  = rnahie1 + hie
            d7_v  = d3_v * rnag1c / d6_v
            d10_v = d7_v * rnahpol1 * hpol
            d4_b  = div1 * hpol
            d5_b  = d4_b * rnahpol1
            d6_b  = d4_b * d7_v
            d7_b  = d5_b * (1.0e0_wp / d6_v)
            d11_b = d7_b * d3_v
            d9_b  = d5_b * ((-d7_v) / d6_v) + d7_b * rnag1c * hie * (-1.5e0_wp )
            do g_i_ = 1, g_p_
                g_c61(g_i_) = d10_v * g_div1(g_i_) + d6_b * g_rnahpol1(g_i_) + d11_b * g_rnag1c(g_i_) + d9_b * &
                    & g_rnahie1(g_i_)
            enddo
            c61  = d10_v * div1
            d2_v = (-39.0625e0_wp) * rnag1c
            d4_b = rnag1 * (-39.0625e0_wp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d2_v * g_rnag1(g_i_) + d4_b * g_rnag1c(g_i_)
            enddo
            d1_w = d2_v * rnag1
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_dexp1(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            dexp1 = d2_v
            d4_v  = rnahpol1b / rnahpol1
            d5_v  = d4_v - 1.0e0_wp
            d8_v  = 1.0e0_wp + d5_v * dexp1
            d6_b  = c61 * d5_v
            d7_b  = c61 * dexp1
            d8_b  = d7_b * (1.0e0_wp / rnahpol1)
            d9_b  = d7_b * ((-d4_v) / rnahpol1)
            do g_i_ = 1, g_p_
                g_c61(g_i_) = d6_b * g_dexp1(g_i_) + d9_b * g_rnahpol1(g_i_) + d8_b * g_rnahpol1b(g_i_) + d8_v * &
                    & g_c61(g_i_)
            enddo
            c61   = c61 * d8_v
            d3_v  = (-1.5e0_wp) * rnahie3 * hie
            d6_v  = rnahie3 + hie
            d7_v  = d3_v * rnag3c / d6_v
            d10_v = d7_v * rnahpol3 * hpol
            d4_b  = div3 * hpol
            d5_b  = d4_b * rnahpol3
            d6_b  = d4_b * d7_v
            d7_b  = d5_b * (1.0e0_wp / d6_v)
            d11_b = d7_b * d3_v
            d9_b  = d5_b * ((-d7_v) / d6_v) + d7_b * rnag3c * hie * (-1.5e0_wp )
            do g_i_ = 1, g_p_
                g_c63(g_i_) = d10_v * g_div3(g_i_) + d6_b * g_rnahpol3(g_i_) + d11_b * g_rnag3c(g_i_) + d9_b * &
                    & g_rnahie3(g_i_)
            enddo
            c63  = d10_v * div3
            d2_v = (-39.0625e0_wp) * rnag3c
            d4_b = rnag3 * (-39.0625e0_wp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d2_v * g_rnag3(g_i_) + d4_b * g_rnag3c(g_i_)
            enddo
            d1_w = d2_v * rnag3
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_dexp3(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            dexp3 = d2_v
            d4_v  = rnahpol3b / rnahpol3
            d5_v  = d4_v - 1.0e0_wp
            d8_v  = 1.0e0_wp + d5_v * dexp3
            d6_b  = c63 * d5_v
            d7_b  = c63 * dexp3
            d8_b  = d7_b * (1.0e0_wp / rnahpol3)
            d9_b  = d7_b * ((-d4_v) / rnahpol3)
            do g_i_ = 1, g_p_
                g_c63(g_i_) = d6_b * g_dexp3(g_i_) + d9_b * g_rnahpol3(g_i_) + d8_b * g_rnahpol3b(g_i_) + d8_v * &
                    & g_c63(g_i_)
            enddo
            c63   = c63 * d8_v
            d4_v  = 1.5e0_wp * cstsq1 + 0.5e0_wp
            d5_v  = 1.0e0_wp / d4_v
            d6_v  = d5_v - 1.0e0_wp
            d9_v  = 1.0e0_wp + d6_v * dexp1
            d6_b  = rmuimu1 * d6_v
            d10_b = rmuimu1 * dexp1 * ((-d5_v) / d4_v) * 1.5e0_wp
            do g_i_ = 1, g_p_
                g_cind1(g_i_) = d6_b * g_dexp1(g_i_) + d10_b * g_cstsq1(g_i_ ) + d9_v * g_rmuimu1(g_i_)
            enddo
            cind1 = rmuimu1 * d9_v
            d4_v  = 1.5e0_wp * cstsq3 + 0.5e0_wp
            d5_v  = 1.0e0_wp / d4_v
            d6_v  = d5_v - 1.0e0_wp
            d9_v  = 1.0e0_wp + d6_v * dexp3
            d6_b  = rmuimu3 * d6_v
            d10_b = rmuimu3 * dexp3 * ((-d5_v) / d4_v) * 1.5e0_wp
            do g_i_ = 1, g_p_
                g_cind3(g_i_) = d6_b * g_dexp3(g_i_) + d10_b * g_cstsq3(g_i_ ) + d9_v * g_rmuimu3(g_i_)
            enddo
            cind3 = rmuimu3 * d9_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = balp * g_r(g_i_,i,2)
            enddo
            d1_w = balp * (r(i,2) - brho)
            d3_v = c6 * cevau * 0.5e0_wp
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d6_v = 1.0e0_wp - d5_v
            d5_b = (-d3_v) * d1_p
            d7_b = d6_v * 0.5e0_wp * cevau
            do g_i_ = 1, g_p_
                g_c6(g_i_) = d5_b * g_d1_w(g_i_) + d7_b * g_c6(g_i_)
            enddo
            c6 = d3_v * d6_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = dalp * g_r(g_i_,i,1)
            enddo
            d1_w = dalp * (r(i,1) - drho)
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_c6swt1(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            c6swt1 = 0.5e0_wp * (1.0e0_wp - d2_v)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = dalp * g_r(g_i_,i,3)
            enddo
            d1_w = dalp * (r(i,3) - drho)
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_c6swt3(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            c6swt3 = 0.5e0_wp * (1.0e0_wp - d2_v)
            d3_b   = cevau * c6swt1
            d4_b   = cevau * c61
            do g_i_ = 1, g_p_
                g_c61(g_i_) = d4_b * g_c6swt1(g_i_) + d3_b * g_c61(g_i_)
            enddo
            c61  = c61 * c6swt1 * cevau
            d3_b = cevau * c6swt3
            d4_b = cevau * c63
            do g_i_ = 1, g_p_
                g_c63(g_i_) = d4_b * g_c6swt3(g_i_) + d3_b * g_c63(g_i_)
            enddo
            c63 = c63 * c6swt3 * cevau
            do g_i_ = 1, g_p_
                g_cind1(g_i_) = cind1 * g_c6swt1(g_i_) + c6swt1 * g_cind1(g_i_)
            enddo
            cind1 = cind1 * c6swt1
            do g_i_ = 1, g_p_
                g_cind3(g_i_) = cind3 * g_c6swt3(g_i_) + c6swt3 * g_cind3(g_i_)
            enddo
            cind3 = cind3 * c6swt3
            do g_i_ = 1, g_p_
                g_elec(g_i_) = g_c6(g_i_) + g_cind3(g_i_) + g_cind1(g_i_) + g_c63(g_i_) + g_c61(g_i_)
            enddo
            elec = c61 + c63 + cind1 + cind3 + c6
            ! here is where it all gets put together . . .
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
            d8_v = sqrt(d1_w)
            if( d1_w .gt. 0.0e0_wp )then
                d1_p = 1.0e0_wp / (2.0e0_wp *  d8_v)
            else
                d1_p = 0.0e0_wp
            endif
            d6_b = (-cevau) * rac2 * d1_p
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d6_b * g_d1_w(g_i_) + cevau * g_coul3(g_i_) + cevau * g_coul2(g_i_) + cevau * &
                    & g_coul1(g_i_)
            enddo
            e(i) = (coul1 + coul2 + coul3 + deh2 - rac2 * d8_v) * cevau
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = g_elec(g_i_) + g_e(g_i_,i)
            enddo
            e(i) = e(i) + elec
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_6_2000_g_pot11_legacy_state_a1118117

!> @brief Compatibility interface for g_pot11.
subroutine g_pot11( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_nah2_leps_6_2000_g_pot11_legacy_state_a1118117, only: dispatch
    use chempotpy_kinds,                                          only: wp
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

! U12 coupling surface . . .
!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_nah2_leps_6_2000_g_pot12_legacy_state_a1118117
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
        real(wp) :: r(:,:), e(:), amp, erl, ers, eps, rmx, rmod, ralp, frho, thalp, thrho, phalp, phrho, r1s, r2s, r3s, &
            & cevau, pi, f1, f3, rmodswt1, rmodswt3, f1p, f3p, f1mod, f3mod, thetaswt, phi, cosa, rnag, rnagsq, chi, falp, &
            & sinchi, alpha1, alpha3, dtanh
        integer            :: nt, i, j, k
        logical            :: leven
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d2_p, d12_b, d7_b, d8_b, d9_v, d10_b, d10_v, d2_v, d9_b, d2_b
        real(wp)           :: d1_p, d3_v, d4_v, d5_v, d6_v, d3_b, d4_b, d5_b, &
            & d6_b, d1_w
        real(wp) :: d7_v, d8_v, g_r1s(g_pmax_), g_r(:,:,:), &
            & g_r2s(g_pmax_), g_r3s(g_pmax_), g_alpha1(g_pmax_), g_alpha3(g_pmax_), g_d1_w(g_pmax_), g_f1(g_pmax_)
        real(wp) :: g_f3(g_pmax_), g_rnagsq(g_pmax_), g_rnag(g_pmax_), g_chi(g_pmax_), g_cosa(g_pmax_), g_sinchi(g_pmax_), g_rmodswt1 &
            & (g_pmax_), g_rmodswt3(g_pmax_), g_f1p(g_pmax_), g_f3p(g_pmax_)
        real(wp) :: g_f1mod(g_pmax_), g_f3mod(g_pmax_), g_thetaswt( &
            & g_pmax_), g_phi(g_pmax_), g_e(:,:)
        integer :: g_ehfid
        alpha1     = 0.e0_wp
        alpha3     = 0.e0_wp
        amp        = 0.e0_wp
        cevau      = 0.e0_wp
        chi        = 0.e0_wp
        cosa       = 0.e0_wp
        d10_b      = 0.e0_wp
        d10_v      = 0.e0_wp
        d12_b      = 0.e0_wp
        d1_p       = 0.e0_wp
        d1_w       = 0.e0_wp
        d2_b       = 0.e0_wp
        d2_p       = 0.e0_wp
        d2_v       = 0.e0_wp
        d3_b       = 0.e0_wp
        d3_v       = 0.e0_wp
        d4_b       = 0.e0_wp
        d4_v       = 0.e0_wp
        d5_b       = 0.e0_wp
        d5_v       = 0.e0_wp
        d6_b       = 0.e0_wp
        d6_v       = 0.e0_wp
        d7_b       = 0.e0_wp
        d7_v       = 0.e0_wp
        d8_b       = 0.e0_wp
        d8_v       = 0.e0_wp
        d9_b       = 0.e0_wp
        d9_v       = 0.e0_wp
        dtanh      = 0.e0_wp
        eps        = 0.e0_wp
        erl        = 0.e0_wp
        ers        = 0.e0_wp
        f1         = 0.e0_wp
        f1mod      = 0.e0_wp
        f1p        = 0.e0_wp
        f3         = 0.e0_wp
        f3mod      = 0.e0_wp
        f3p        = 0.e0_wp
        falp       = 0.e0_wp
        frho       = 0.e0_wp
        g_alpha1   = 0.e0_wp
        g_alpha3   = 0.e0_wp
        g_chi      = 0.e0_wp
        g_cosa     = 0.e0_wp
        g_d1_w     = 0.e0_wp
        g_f1       = 0.e0_wp
        g_f1mod    = 0.e0_wp
        g_f1p      = 0.e0_wp
        g_f3       = 0.e0_wp
        g_f3mod    = 0.e0_wp
        g_f3p      = 0.e0_wp
        g_phi      = 0.e0_wp
        g_r1s      = 0.e0_wp
        g_r2s      = 0.e0_wp
        g_r3s      = 0.e0_wp
        g_rmodswt1 = 0.e0_wp
        g_rmodswt3 = 0.e0_wp
        g_rnag     = 0.e0_wp
        g_rnagsq   = 0.e0_wp
        g_sinchi   = 0.e0_wp
        g_thetaswt = 0.e0_wp
        phalp      = 0.e0_wp
        phi        = 0.e0_wp
        phrho      = 0.e0_wp
        pi         = 0.e0_wp
        r1s        = 0.e0_wp
        r2s        = 0.e0_wp
        r3s        = 0.e0_wp
        ralp       = 0.e0_wp
        rmod       = 0.e0_wp
        rmodswt1   = 0.e0_wp
        rmodswt3   = 0.e0_wp
        rmx        = 0.e0_wp
        rnag       = 0.e0_wp
        rnagsq     = 0.e0_wp
        sinchi     = 0.e0_wp
        thalp      = 0.e0_wp
        thetaswt   = 0.e0_wp
        thrho      = 0.e0_wp
        amp        = 1.5150e0_wp
        erl        = 0.08877e0_wp
        ers        = 11.630e0_wp
        eps        = 1.0e-4_wp
        rmx        = 2.622e0_wp
        rmod       = 0.5409e0_wp
        falp       = 0.7583e0_wp
        frho       = 5.1858e0_wp
        thalp      = 2.4882e0_wp
        thrho      = 2.0e0_wp
        phalp      = 0.40787e0_wp
        phrho      = 7.39370e0_wp
        g_ehfid    = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        cevau = 1.e0_wp / 27.2113961e0_wp
        pi    = acos(-1.0e0_wp)
        ! set leven true for an even U12 surface, and false for an odd
        ! surface
        leven = .true.
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
            r3s = d2_v
            ! here we do the function evaluations . . .
            d4_v = eps + r(i,1) * r1s
            d5_v = ers / d4_v
            d4_b = (-d5_v) / d4_v
            d5_b = d4_b * r1s
            d6_b = d4_b * r(i,1)
            do g_i_ = 1, g_p_
                g_alpha1(g_i_) = d6_b * g_r1s(g_i_) + d5_b * g_r(g_i_,i,1)
            enddo
            alpha1 = erl + d5_v
            d4_v   = eps + r(i,3) * r3s
            d5_v   = ers / d4_v
            d4_b   = (-d5_v) / d4_v
            d5_b   = d4_b * r3s
            d6_b   = d4_b * r(i,3)
            do g_i_ = 1, g_p_
                g_alpha3(g_i_) = d6_b * g_r3s(g_i_) + d5_b * g_r(g_i_,i,3)
            enddo
            alpha3 = erl + d5_v
            d5_v   = (r(i,1) - rmx) * (r(i,1) - rmx)
            d1_p   = 2.0e0_wp * (r(i,1) - rmx)
            d5_b   = (-alpha1) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_r(g_i_,i,1) + (-d5_v) * g_alpha1(g_i_)
            enddo
            d1_w = (-alpha1) * d5_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            d3_b = amp * d1_p
            do g_i_ = 1, g_p_
                g_f1(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            f1   = amp * d2_v
            d5_v = (r(i,3) - rmx) * (r(i,3) - rmx)
            d1_p = 2.0e0_wp * (r(i,3) - rmx)
            d5_b = (-alpha3) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_r(g_i_,i,3) + (-d5_v) * g_alpha3(g_i_)
            enddo
            d1_w = (-alpha3) * d5_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            d3_b = amp * d1_p
            do g_i_ = 1, g_p_
                g_f3(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            f3 = amp * d2_v
            ! here follows some checking procedures
            do g_i_ = 1, g_p_
                g_rnagsq(g_i_) = (-0.25e0_wp) * g_r2s(g_i_) + 0.5e0_wp * g_r3s(g_i_) + 0.5e0_wp * g_r1s(g_i_)
            enddo
            rnagsq = 0.5e0_wp * r1s + 0.5e0_wp * r3s - 0.25e0_wp * r2s
            if( rnagsq .gt. 0.0e0_wp )then
                rnag = sqrt(rnagsq)
                if( rnag .gt. 0.0e0_wp )then
                    g_rnag(1) = 0.5e0_wp*r(i,1)/rnag
                    g_rnag(2) = -0.25e0_wp*r(i,2)/rnag
                    g_rnag(3) = 0.5e0_wp*r(i,3)/rnag
                else
                    g_rnag(1) = 0.0e0_wp
                    g_rnag(2) = 0.0e0_wp
                    g_rnag(3) = 0.0e0_wp
                endif
                d5_v = 2.0e0_wp * rnag
                d7_v = d5_v * r(i,2)
                d8_v = (r3s - r1s) / d7_v
                d9_v = d8_v * d8_v
                d1_p = 2.0e0_wp * d8_v
                d3_b = d1_p * (1.0e0_wp / d7_v)
                d4_b = d1_p * ((-d8_v) / d7_v)
                d6_b = d4_b * d5_v
                d7_b = d4_b * r(i,2) * 2.0e0_wp
                do g_i_ = 1, g_p_
                    g_chi(g_i_) = d6_b * g_r(g_i_,i,2) + d7_b * g_rnag(g_i_) + (-d3_b) * g_r1s(g_i_) + d3_b * &
                        & g_r3s(g_i_)
                enddo
                chi = d9_v
            else
                do g_i_ = 1, g_p_
                    g_rnag(g_i_) = 0.0e0_wp
                enddo
                rnag = 0.0e0_wp
                do g_i_ = 1, g_p_
                    g_chi(g_i_) = 0.0e0_wp
                enddo
                chi = 0.0e0_wp
            endif
            d7_v  = 2.0e0_wp * r(i,1)
            d9_v  = d7_v * r(i,3)
            d10_v = (r1s + r3s - r2s) / d9_v
            d2_b  = 1.0e0_wp / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,3) * 2.0e0_wp
            do g_i_ = 1, g_p_
                g_cosa(g_i_) = d5_b * g_r(g_i_,i,3) + d6_b * g_r(g_i_,i,1) + (-d2_b) * g_r2s(g_i_) + d2_b * &
                    & g_r3s(g_i_) + d2_b * g_r1s(g_i_)
            enddo
            cosa = d10_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-10.0e0_wp) * g_rnagsq(g_i_)
            enddo
            d1_w  = (-10.0e0_wp) * rnagsq
            d2_v  = 1.0e0_wp - chi
            d5_v  = 0.5e0_wp * cosa + 0.5e0_wp
            d6_v  = d2_v * d5_v
            d8_v  = exp(d1_w)
            d1_p  = d8_v
            d9_v  = 1.0e0_wp - d8_v
            d5_b  = (-d6_v) * d1_p
            d9_b  = d9_v * d2_v * 0.5e0_wp
            d10_b = -(d9_v * d5_v)
            do g_i_ = 1, g_p_
                g_sinchi(g_i_) = d5_b * g_d1_w(g_i_) + d9_b * g_cosa(g_i_) + d10_b * g_chi(g_i_)
            enddo
            sinchi = d6_v * d9_v
            ! here ends the procedure
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = falp * g_r(g_i_,i,3)
            enddo
            d1_w = falp * (r(i,3) - frho)
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_v = 0.5e0_wp - 0.5e0_wp * d2_v
            d7_v = 0.5e0_wp + 0.5e0_wp * cosa
            d5_b = d4_v * 0.5e0_wp
            d8_b = (-d7_v) * 0.5e0_wp * d1_p
            do g_i_ = 1, g_p_
                g_rmodswt1(g_i_) = d5_b * g_cosa(g_i_) + d8_b * g_d1_w(g_i_)
            enddo
            rmodswt1 = d4_v * d7_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = falp * g_r(g_i_,i,1)
            enddo
            d1_w = falp * (r(i,1) - frho)
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_v = 0.5e0_wp - 0.5e0_wp * d2_v
            d7_v = 0.5e0_wp + 0.5e0_wp * cosa
            d5_b = d4_v * 0.5e0_wp
            d8_b = (-d7_v) * 0.5e0_wp * d1_p
            do g_i_ = 1, g_p_
                g_rmodswt3(g_i_) = d5_b * g_cosa(g_i_) + d8_b * g_d1_w(g_i_)
            enddo
            rmodswt3 = d4_v * d7_v
            d5_v     = (r(i,1) - rmx * rmod) * (r(i,1) - rmx * rmod)
            d1_p     = 2.0e0_wp * (r(i,1) - rmx * rmod)
            d5_b     = (-alpha1) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_r(g_i_,i,1) + (-d5_v) * g_alpha1(g_i_)
            enddo
            d1_w = (-alpha1) * d5_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            d3_b = amp * d1_p
            do g_i_ = 1, g_p_
                g_f1p(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            f1p  = amp * d2_v
            d5_v = (r(i,3) - rmx * rmod) * (r(i,3) - rmx * rmod)
            d1_p = 2.0e0_wp * (r(i,3) - rmx * rmod)
            d5_b = (-alpha3) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_r(g_i_,i,3) + (-d5_v) * g_alpha3(g_i_)
            enddo
            d1_w = (-alpha3) * d5_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            d3_b = amp * d1_p
            do g_i_ = 1, g_p_
                g_f3p(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            f3p  = amp * d2_v
            d4_v = f1p - f1
            d3_b = 1.0e0_wp + (-rmodswt1)
            do g_i_ = 1, g_p_
                g_f1mod(g_i_) = d3_b * g_f1(g_i_) + rmodswt1 * g_f1p(g_i_) + d4_v * g_rmodswt1(g_i_)
            enddo
            f1mod = rmodswt1 * d4_v + f1
            d4_v  = f3p - f3
            d3_b  = 1.0e0_wp + (-rmodswt3)
            do g_i_ = 1, g_p_
                g_f3mod(g_i_) = d3_b * g_f3(g_i_) + rmodswt3 * g_f3p(g_i_) + d4_v * g_rmodswt3(g_i_)
            enddo
            f3mod = rmodswt3 * d4_v + f3
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = thalp * g_r(g_i_,i,2)
            enddo
            d1_w = thalp * (r(i,2) - thrho)
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_thetaswt(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            thetaswt = 0.5e0_wp - 0.5e0_wp * d2_v
            d3_v     = f1mod - f3mod
            d6_v     = sinchi - 1.0e0_wp
            d8_v     = thetaswt * d6_v + 1.0e0_wp
            d5_b     = d3_v * d6_v
            d7_b     = d3_v * thetaswt
            do g_i_ = 1, g_p_
                g_phi(g_i_) = d7_b * g_sinchi(g_i_) + d5_b * g_thetaswt(g_i_ ) + (-d8_v) * g_f3mod(g_i_) + d8_v * &
                    & g_f1mod(g_i_)
            enddo
            phi = d3_v * d8_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = phalp * g_r(g_i_,i,2)
            enddo
            d1_w  = phalp * (r(i,2) - phrho)
            d2_v  = tanh (d1_w)
            d1_p  = 1.0e0_wp - ( d2_v *  d2_v)
            d4_v  = 0.5e0_wp * (1.e0_wp - d2_v)
            d7_v  = 0.5e0_wp - 0.5e0_wp * cosa
            d9_v  = 1.0e0_wp - d4_v * d7_v
            d4_b  = cevau * d9_v
            d5_b  = -(cevau * phi)
            d9_b  = (-(d5_b * d4_v)) * 0.5e0_wp
            d12_b = (-(d5_b * d7_v * 0.5e0_wp)) * d1_p
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d4_b * g_phi(g_i_) + d9_b * g_cosa(g_i_) + d12_b * g_d1_w(g_i_)
            enddo
            e(i) = d9_v * phi * cevau
            ! set leven equal to true for an even surface, false for an odd one
            if( leven )then
                d2_v = e(i) ** ( 4 - 2)
                d2_v = d2_v * e(i)
                d2_p = 4 *  d2_v
                d2_v = d2_v * e(i)
                d3_v = e(i) * e(i)
                d1_p = 2.0e0_wp * e(i)
                d4_v = d3_v + 1.0e-8_wp
                d5_v = d2_v / d4_v
                d5_b = (-d5_v) / d4_v * d1_p + 1.0e0_wp / d4_v * d2_p
                do g_i_ = 1, g_p_
                    g_d1_w(g_i_) = d5_b * g_e(g_i_,i)
                enddo
                d1_w = d5_v
                d2_v = sqrt(d1_w)
                if( d1_w .gt. 0.0e0_wp )then
                    d1_p = 1.0e0_wp / (2.0e0_wp *  d2_v)
                else
                    d1_p = 0.0e0_wp
                endif
                do g_i_ = 1, g_p_
                    g_e(g_i_,i) = d1_p * g_d1_w(g_i_)
                enddo
                e(i) = d2_v
            endif
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_6_2000_g_pot12_legacy_state_a1118117

!> @brief Compatibility interface for g_pot12.
subroutine g_pot12( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_nah2_leps_6_2000_g_pot12_legacy_state_a1118117, only: dispatch
    use chempotpy_kinds,                                          only: wp
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

! the upper diabatic surface U22
!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_nah2_leps_6_2000_g_pot22_legacy_state_a1118117
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
        real(wp) :: al1, al2, alb, alhh, alt, balp, brho, c2a, c2b, c2c, c6, c61, c63, c6swt1, c6swt3, cevau, cind1, &
            & cind3, cosa, coscos, coul1, coul2, coul3, cplg2, cstsq, cstsq1, cstsq3, dalp, deh2, dexp1, dexp2, dexp3, div1, &
            & div3, drho, e(:), elec, exch1, exch2, exch3, frac, h2a1, h2a2, h2a3, h2a4, h2a5, h2a6, h2a7, h2a8, h2a9, h2b1, &
            & h2b2, h2b3, h2b4, h2b5, h2b6, h2b7, h2b8, h2b9, h2hc1, h2hc2, h2hc3, h2hc4, h2hc5, h2hex, h2i1, h2i2, h2i3, h2i4, &
            & h2i5, h2i6, h2ie, h2qc1, h2qc2, h2qc3, h2qc4, h2qc5, h2quad, hie, hnaa1, hnaa2, hnaa3, hnaa4, hnaa5, hnab1, hnab2, &
            & hnab3, hnab4, hnab5, hnai1, hnai2, hnai3, hnai4, hnai5, hnai6, hnau1, hnau2, hnau3, hnau4, hnau5, hnau6, hpol, &
            & onethd, pi, qh, qhint, qq, qqint, r(:,:), r1s, r2s, r3s, rb, rd, rmuimu1, rmuimu3, rnag, rnag1, rnag1c, rnag3, &
            & rnag3c, rnagc, rnags, rnahie1, rnahie3, rnahmu1, rnahmu3, rnahpar1, rnahpar3, rnahperp1, rnahperp3, rnahpol1, &
            & rnahpol1b, rnahpol3, rnahpol3b, rnaie, rnapol1, rnaquad, rqh, rqq, s1, s1a, s1c1, s1c2, s1c3, s1c4, s1c5, s1c6, &
            & s1i, s1ialp, s1ic, s1ic1, s1ic2, s1ic3, s1ic4, s1ic5, s1ic6, s2, s290, s290alp, s290c, s290c1, s290c2, s290c3, &
            & s290c4, s290c5, s290c6, s290rho, s2a, s2c1, s2c2, s2c3, s2c4, s2c5, s2c6, s2mod, s3, s3a, s3i, s3ic, t1, t10, &
            & t10c1, t10c2, t10c3, t10c4, t10c5, t10c6, t190, t190c1, t190c2, t190c3, t190c4, t190c5, t190c6, t2, t20, t20alp, &
            & t20c, t20c1, t20c2, t20c3, t20c4, t20c5, t20c6, t20rho, t290, t290alp, t290c, t290c1, t290c2, t290c3, t290c4, &
            & t290c5, t290c6, t290rho, t2a, t2c1, t2c2, t2c3, t2c4, t2c5, t2c6, t3, t30, t390, twothd, una2p, w
        integer :: i, j, k, nsurf, nt
        ! the NaH singlet data blocks
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d3_p, d24_b, d2_p, d15_v, d20_v, d16_v, d22_v, &
            & d2_v, d16_b, d2_b
        real(wp) :: d1_w, d2_w, d3_v, d4_v, d5_v, d6_v, d3_b, d4_b, &
            & d5_b, d6_b
        real(wp) :: d1_p, d7_v, d8_v, d9_v, d7_b, d8_b, d9_b, d10_v, d10_b, d11_v
        real(wp) :: d12_v, d13_v, d19_v, d11_b, d12_b, d13_b, d14_b, g_r1s(g_pmax_), g_r(:,:,:), g_r2s(g_pmax_)
        real(wp) :: g_r3s(g_pmax_), g_d2_w(g_pmax_), g_d1_w(g_pmax_ &
            & ), g_rnag(g_pmax_), g_rnags(g_pmax_), g_rnagc(g_pmax_), g_cstsq(g_pmax_), g_cosa(g_pmax_), g_coscos(g_pmax_), g_rnag1(g_pmax_)
        real(wp) :: g_rnag3(g_pmax_), g_rnag1c(g_pmax_), g_rnag3c(g_pmax_), g_cstsq1(g_pmax_), g_cstsq3(g_pmax_), g_s1a(g_pmax_), g_s1ic(g_pmax_), g_s1i(g_pmax_), g_s1(g_pmax_), g_t10(g_pmax_)
        real(wp) :: g_t190(g_pmax_), g_t1(g_pmax_), g_t2a(g_pmax_), &
            & g_t290c(g_pmax_), g_t290(g_pmax_), g_t20c(g_pmax_), g_t20(g_pmax_ &
            & ), g_t2(g_pmax_), g_s2a(g_pmax_), g_s2mod(g_pmax_)
        real(wp) :: g_s290c(g_pmax_), g_s290(g_pmax_), g_s2(g_pmax_ &
            & ), g_s3a(g_pmax_), g_s3ic(g_pmax_), g_s3i(g_pmax_), g_s3(g_pmax_), &
            & g_t30(g_pmax_), g_t390(g_pmax_), g_t3(g_pmax_)
        real(wp) :: g_h2ie(g_pmax_), g_al1(g_pmax_), g_al2(g_pmax_), g_alt(g_pmax_), g_alb(g_pmax_), g_alhh(g_pmax_), g_h2hex(g_pmax_ &
            & ), g_h2quad(g_pmax_), g_rnahie1(g_pmax_), g_rnahie3(g_pmax_)
        real(wp) :: g_rnahmu1(g_pmax_), g_rnahmu3(g_pmax_), g_rnahpar1(g_pmax_), g_rnahpar3(g_pmax_), g_rnahperp1(g_pmax_), g_rnahperp3(g_pmax_), g_rnahpol1(g_pmax_), g_rnahpol1b(g_pmax_), g_rnahpol3 &
            & (g_pmax_), g_rnahpol3b(g_pmax_)
        real(wp) :: g_qqint(g_pmax_), g_qq(g_pmax_), g_dexp2(g_pmax_), g_qhint(g_pmax_), g_qh(g_pmax_), g_div1(g_pmax_), g_rmuimu1(g_pmax_), g_div3(g_pmax_), g_rmuimu3(g_pmax_), g_c6(g_pmax_)
        real(wp) :: g_c61(g_pmax_), g_dexp1(g_pmax_), g_c63(g_pmax_ &
            & ), g_dexp3(g_pmax_), g_cind1(g_pmax_), g_cind3(g_pmax_), g_elec(g_pmax_), g_c6swt1(g_pmax_), g_c6swt3(g_pmax_), g_coul1(g_pmax_)
        real(wp) :: g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_), g_exch3(g_pmax_), g_w(g_pmax_), g_cplg2(g_pmax_), g_e(:,:)
        integer  :: g_ehfid
        ! the NaH triplet data blocks
        ! The H2 triplet
        ! The H2 singlet
        ! data blocks for the long range forces
        ! misc constants
        al1         = 0.e0_wp
        al2         = 0.e0_wp
        alb         = 0.e0_wp
        alhh        = 0.e0_wp
        alt         = 0.e0_wp
        balp        = 0.e0_wp
        brho        = 0.e0_wp
        c2a         = 0.e0_wp
        c2b         = 0.e0_wp
        c2c         = 0.e0_wp
        c6          = 0.e0_wp
        c61         = 0.e0_wp
        c63         = 0.e0_wp
        c6swt1      = 0.e0_wp
        c6swt3      = 0.e0_wp
        cevau       = 0.e0_wp
        cind1       = 0.e0_wp
        cind3       = 0.e0_wp
        cosa        = 0.e0_wp
        coscos      = 0.e0_wp
        coul1       = 0.e0_wp
        coul2       = 0.e0_wp
        coul3       = 0.e0_wp
        cplg2       = 0.e0_wp
        cstsq       = 0.e0_wp
        cstsq1      = 0.e0_wp
        cstsq3      = 0.e0_wp
        d10_b       = 0.e0_wp
        d10_v       = 0.e0_wp
        d11_b       = 0.e0_wp
        d11_v       = 0.e0_wp
        d12_b       = 0.e0_wp
        d12_v       = 0.e0_wp
        d13_b       = 0.e0_wp
        d13_v       = 0.e0_wp
        d14_b       = 0.e0_wp
        d15_v       = 0.e0_wp
        d16_b       = 0.e0_wp
        d16_v       = 0.e0_wp
        d19_v       = 0.e0_wp
        d1_p        = 0.e0_wp
        d1_w        = 0.e0_wp
        d20_v       = 0.e0_wp
        d22_v       = 0.e0_wp
        d24_b       = 0.e0_wp
        d2_b        = 0.e0_wp
        d2_p        = 0.e0_wp
        d2_v        = 0.e0_wp
        d2_w        = 0.e0_wp
        d3_b        = 0.e0_wp
        d3_p        = 0.e0_wp
        d3_v        = 0.e0_wp
        d4_b        = 0.e0_wp
        d4_v        = 0.e0_wp
        d5_b        = 0.e0_wp
        d5_v        = 0.e0_wp
        d6_b        = 0.e0_wp
        d6_v        = 0.e0_wp
        d7_b        = 0.e0_wp
        d7_v        = 0.e0_wp
        d8_b        = 0.e0_wp
        d8_v        = 0.e0_wp
        d9_b        = 0.e0_wp
        d9_v        = 0.e0_wp
        dalp        = 0.e0_wp
        deh2        = 0.e0_wp
        dexp1       = 0.e0_wp
        dexp2       = 0.e0_wp
        dexp3       = 0.e0_wp
        div1        = 0.e0_wp
        div3        = 0.e0_wp
        drho        = 0.e0_wp
        elec        = 0.e0_wp
        exch1       = 0.e0_wp
        exch2       = 0.e0_wp
        exch3       = 0.e0_wp
        frac        = 0.e0_wp
        g_al1       = 0.e0_wp
        g_al2       = 0.e0_wp
        g_alb       = 0.e0_wp
        g_alhh      = 0.e0_wp
        g_alt       = 0.e0_wp
        g_c6        = 0.e0_wp
        g_c61       = 0.e0_wp
        g_c63       = 0.e0_wp
        g_c6swt1    = 0.e0_wp
        g_c6swt3    = 0.e0_wp
        g_cind1     = 0.e0_wp
        g_cind3     = 0.e0_wp
        g_cosa      = 0.e0_wp
        g_coscos    = 0.e0_wp
        g_coul1     = 0.e0_wp
        g_coul2     = 0.e0_wp
        g_coul3     = 0.e0_wp
        g_cplg2     = 0.e0_wp
        g_cstsq     = 0.e0_wp
        g_cstsq1    = 0.e0_wp
        g_cstsq3    = 0.e0_wp
        g_d1_w      = 0.e0_wp
        g_d2_w      = 0.e0_wp
        g_dexp1     = 0.e0_wp
        g_dexp2     = 0.e0_wp
        g_dexp3     = 0.e0_wp
        g_div1      = 0.e0_wp
        g_div3      = 0.e0_wp
        g_elec      = 0.e0_wp
        g_exch1     = 0.e0_wp
        g_exch2     = 0.e0_wp
        g_exch3     = 0.e0_wp
        g_h2hex     = 0.e0_wp
        g_h2ie      = 0.e0_wp
        g_h2quad    = 0.e0_wp
        g_qh        = 0.e0_wp
        g_qhint     = 0.e0_wp
        g_qq        = 0.e0_wp
        g_qqint     = 0.e0_wp
        g_r1s       = 0.e0_wp
        g_r2s       = 0.e0_wp
        g_r3s       = 0.e0_wp
        g_rmuimu1   = 0.e0_wp
        g_rmuimu3   = 0.e0_wp
        g_rnag      = 0.e0_wp
        g_rnag1     = 0.e0_wp
        g_rnag1c    = 0.e0_wp
        g_rnag3     = 0.e0_wp
        g_rnag3c    = 0.e0_wp
        g_rnagc     = 0.e0_wp
        g_rnags     = 0.e0_wp
        g_rnahie1   = 0.e0_wp
        g_rnahie3   = 0.e0_wp
        g_rnahmu1   = 0.e0_wp
        g_rnahmu3   = 0.e0_wp
        g_rnahpar1  = 0.e0_wp
        g_rnahpar3  = 0.e0_wp
        g_rnahperp1 = 0.e0_wp
        g_rnahperp3 = 0.e0_wp
        g_rnahpol1  = 0.e0_wp
        g_rnahpol1b = 0.e0_wp
        g_rnahpol3  = 0.e0_wp
        g_rnahpol3b = 0.e0_wp
        g_s1        = 0.e0_wp
        g_s1a       = 0.e0_wp
        g_s1i       = 0.e0_wp
        g_s1ic      = 0.e0_wp
        g_s2        = 0.e0_wp
        g_s290      = 0.e0_wp
        g_s290c     = 0.e0_wp
        g_s2a       = 0.e0_wp
        g_s2mod     = 0.e0_wp
        g_s3        = 0.e0_wp
        g_s3a       = 0.e0_wp
        g_s3i       = 0.e0_wp
        g_s3ic      = 0.e0_wp
        g_t1        = 0.e0_wp
        g_t10       = 0.e0_wp
        g_t190      = 0.e0_wp
        g_t2        = 0.e0_wp
        g_t20       = 0.e0_wp
        g_t20c      = 0.e0_wp
        g_t290      = 0.e0_wp
        g_t290c     = 0.e0_wp
        g_t2a       = 0.e0_wp
        g_t3        = 0.e0_wp
        g_t30       = 0.e0_wp
        g_t390      = 0.e0_wp
        g_w         = 0.e0_wp
        h2a1        = 0.e0_wp
        h2a2        = 0.e0_wp
        h2a3        = 0.e0_wp
        h2a4        = 0.e0_wp
        h2a5        = 0.e0_wp
        h2a6        = 0.e0_wp
        h2a7        = 0.e0_wp
        h2a8        = 0.e0_wp
        h2a9        = 0.e0_wp
        h2b1        = 0.e0_wp
        h2b2        = 0.e0_wp
        h2b3        = 0.e0_wp
        h2b4        = 0.e0_wp
        h2b5        = 0.e0_wp
        h2b6        = 0.e0_wp
        h2b7        = 0.e0_wp
        h2b8        = 0.e0_wp
        h2b9        = 0.e0_wp
        h2hc1       = 0.e0_wp
        h2hc2       = 0.e0_wp
        h2hc3       = 0.e0_wp
        h2hc4       = 0.e0_wp
        h2hc5       = 0.e0_wp
        h2hex       = 0.e0_wp
        h2i1        = 0.e0_wp
        h2i2        = 0.e0_wp
        h2i3        = 0.e0_wp
        h2i4        = 0.e0_wp
        h2i5        = 0.e0_wp
        h2i6        = 0.e0_wp
        h2ie        = 0.e0_wp
        h2qc1       = 0.e0_wp
        h2qc2       = 0.e0_wp
        h2qc3       = 0.e0_wp
        h2qc4       = 0.e0_wp
        h2qc5       = 0.e0_wp
        h2quad      = 0.e0_wp
        hie         = 0.e0_wp
        hnaa1       = 0.e0_wp
        hnaa2       = 0.e0_wp
        hnaa3       = 0.e0_wp
        hnaa4       = 0.e0_wp
        hnaa5       = 0.e0_wp
        hnab1       = 0.e0_wp
        hnab2       = 0.e0_wp
        hnab3       = 0.e0_wp
        hnab4       = 0.e0_wp
        hnab5       = 0.e0_wp
        hnai1       = 0.e0_wp
        hnai2       = 0.e0_wp
        hnai3       = 0.e0_wp
        hnai4       = 0.e0_wp
        hnai5       = 0.e0_wp
        hnai6       = 0.e0_wp
        hnau1       = 0.e0_wp
        hnau2       = 0.e0_wp
        hnau3       = 0.e0_wp
        hnau4       = 0.e0_wp
        hnau5       = 0.e0_wp
        hnau6       = 0.e0_wp
        hpol        = 0.e0_wp
        onethd      = 0.e0_wp
        pi          = 0.e0_wp
        qh          = 0.e0_wp
        qhint       = 0.e0_wp
        qq          = 0.e0_wp
        qqint       = 0.e0_wp
        r1s         = 0.e0_wp
        r2s         = 0.e0_wp
        r3s         = 0.e0_wp
        rb          = 0.e0_wp
        rd          = 0.e0_wp
        rmuimu1     = 0.e0_wp
        rmuimu3     = 0.e0_wp
        rnag        = 0.e0_wp
        rnag1       = 0.e0_wp
        rnag1c      = 0.e0_wp
        rnag3       = 0.e0_wp
        rnag3c      = 0.e0_wp
        rnagc       = 0.e0_wp
        rnags       = 0.e0_wp
        rnahie1     = 0.e0_wp
        rnahie3     = 0.e0_wp
        rnahmu1     = 0.e0_wp
        rnahmu3     = 0.e0_wp
        rnahpar1    = 0.e0_wp
        rnahpar3    = 0.e0_wp
        rnahperp1   = 0.e0_wp
        rnahperp3   = 0.e0_wp
        rnahpol1    = 0.e0_wp
        rnahpol1b   = 0.e0_wp
        rnahpol3    = 0.e0_wp
        rnahpol3b   = 0.e0_wp
        rnaie       = 0.e0_wp
        rnapol1     = 0.e0_wp
        rnaquad     = 0.e0_wp
        rqh         = 0.e0_wp
        rqq         = 0.e0_wp
        s1          = 0.e0_wp
        s1a         = 0.e0_wp
        s1c1        = 0.e0_wp
        s1c2        = 0.e0_wp
        s1c3        = 0.e0_wp
        s1c4        = 0.e0_wp
        s1c5        = 0.e0_wp
        s1c6        = 0.e0_wp
        s1i         = 0.e0_wp
        s1ialp      = 0.e0_wp
        s1ic        = 0.e0_wp
        s1ic1       = 0.e0_wp
        s1ic2       = 0.e0_wp
        s1ic3       = 0.e0_wp
        s1ic4       = 0.e0_wp
        s1ic5       = 0.e0_wp
        s1ic6       = 0.e0_wp
        s2          = 0.e0_wp
        s290        = 0.e0_wp
        s290alp     = 0.e0_wp
        s290c       = 0.e0_wp
        s290c1      = 0.e0_wp
        s290c2      = 0.e0_wp
        s290c3      = 0.e0_wp
        s290c4      = 0.e0_wp
        s290c5      = 0.e0_wp
        s290c6      = 0.e0_wp
        s290rho     = 0.e0_wp
        s2a         = 0.e0_wp
        s2c1        = 0.e0_wp
        s2c2        = 0.e0_wp
        s2c3        = 0.e0_wp
        s2c4        = 0.e0_wp
        s2c5        = 0.e0_wp
        s2c6        = 0.e0_wp
        s2mod       = 0.e0_wp
        s3          = 0.e0_wp
        s3a         = 0.e0_wp
        s3i         = 0.e0_wp
        s3ic        = 0.e0_wp
        t1          = 0.e0_wp
        t10         = 0.e0_wp
        t10c1       = 0.e0_wp
        t10c2       = 0.e0_wp
        t10c3       = 0.e0_wp
        t10c4       = 0.e0_wp
        t10c5       = 0.e0_wp
        t10c6       = 0.e0_wp
        t190        = 0.e0_wp
        t190c1      = 0.e0_wp
        t190c2      = 0.e0_wp
        t190c3      = 0.e0_wp
        t190c4      = 0.e0_wp
        t190c5      = 0.e0_wp
        t190c6      = 0.e0_wp
        t2          = 0.e0_wp
        t20         = 0.e0_wp
        t20alp      = 0.e0_wp
        t20c        = 0.e0_wp
        t20c1       = 0.e0_wp
        t20c2       = 0.e0_wp
        t20c3       = 0.e0_wp
        t20c4       = 0.e0_wp
        t20c5       = 0.e0_wp
        t20c6       = 0.e0_wp
        t20rho      = 0.e0_wp
        t290        = 0.e0_wp
        t290alp     = 0.e0_wp
        t290c       = 0.e0_wp
        t290c1      = 0.e0_wp
        t290c2      = 0.e0_wp
        t290c3      = 0.e0_wp
        t290c4      = 0.e0_wp
        t290c5      = 0.e0_wp
        t290c6      = 0.e0_wp
        t290rho     = 0.e0_wp
        t2a         = 0.e0_wp
        t2c1        = 0.e0_wp
        t2c2        = 0.e0_wp
        t2c3        = 0.e0_wp
        t2c4        = 0.e0_wp
        t2c5        = 0.e0_wp
        t2c6        = 0.e0_wp
        t3          = 0.e0_wp
        t30         = 0.e0_wp
        t390        = 0.e0_wp
        twothd      = 0.e0_wp
        una2p       = 0.e0_wp
        w           = 0.e0_wp
        s1c1        = 170.769e0_wp
        s1c2        = 1.514e0_wp
        s1c3        = -87.110e0_wp
        s1c4        = 57.480e0_wp
        s1c5        = -14.228e0_wp
        s1c6        = 1.303e0_wp
        s1ic1       = 460.09449e0_wp
        s1ic2       = 2.45039e0_wp
        s1ic3       = -29.81890e0_wp
        s1ic4       = 42.67717e0_wp
        s1ic5       = -13.11811e0_wp
        s1ic6       = 1.303e0_wp
        s1ialp      = 0.28780e0_wp
        t10c1       = 385.59843e0_wp
        t10c2       = 2.83307e0_wp
        t10c3       = 67.06299e0_wp
        t10c4       = -37.70079e0_wp
        t10c5       = 6.59055e0_wp
        t10c6       = 1.303e0_wp
        t190c1      = 366.53543e0_wp
        t190c2      = 1.74567e0_wp
        t190c3      = -50.0e0_wp
        t190c4      = 4.20472e0_wp
        t190c5      = 1.22047e0_wp
        t190c6      = 1.303e0_wp
        t2c1        = 147.214e0_wp
        t2c2        = 3.915e0_wp
        t2c3        = 41.275e0_wp
        t2c4        = 10.505e0_wp
        t2c5        = 4.408e0_wp
        t2c6        = 2.032e0_wp
        t290c1      = 82.88976e0_wp
        t290c2      = 3.38898e0_wp
        t290c3      = 10.87402e0_wp
        t290c4      = -7.88189e0_wp
        t290c5      = -2.35433e0_wp
        t290c6      = 2.032e0_wp
        t290alp     = 0.09543e0_wp
        t290rho     = 40.81102e0_wp
        t20c1       = 159.34646e0_wp
        t20c2       = 3.51260e0_wp
        t20c3       = 121.65354e0_wp
        t20c4       = 9.41732e0_wp
        t20c5       = -11.03937e0_wp
        t20c6       = 2.032
        t20alp      = 0.11260e0_wp
        t20rho      = 43.22799504e0_wp
        s2c1        = 139.7160e0_wp
        s2c2        = -123.8978e0_wp
        s2c3        = 3.4031e0_wp
        s2c4        = -6.8725e0_wp
        s2c5        = -23.0440e0_wp
        s2c6        = 2.032e0_wp
        s290c1      = 133.63780e0_wp
        s290c2      = 422.23622e0_wp
        s290c3      = 5.56890e0_wp
        s290c4      = -4.22441e0_wp
        s290c5      = -9.79213e0_wp
        s290c6      = 2.032e0_wp
        s290alp     = 0.05906e0_wp
        s290rho     = 12.72449e0_wp
        h2a1        = 2.38029900e0_wp
        h2a2        = 1.47400334e0_wp
        h2a3        = -0.226331771e0_wp
        h2a4        = -2.17339562e0_wp
        h2a5        = 9.54597998e0_wp
        h2a6        = -23.3366749e0_wp
        h2a7        = 24.7115748e0_wp
        h2a8        = -12.8786102e0_wp
        h2a9        = 2.75711758e0_wp
        h2b1        = 1.35164094e0_wp
        h2b2        = -0.254404564e0_wp
        h2b3        = -10.6627055e0_wp
        h2b4        = -14.0112406e0_wp
        h2b5        = -7.80660398e0_wp
        h2b6        = 5.47697753e0_wp
        h2b7        = -6.95875664e0_wp
        h2b8        = 3.64334771e0_wp
        h2b9        = -0.492505596e0_wp
        h2i1        = 160.484e0_wp
        h2i2        = -20.397e0_wp
        h2i3        = 3.913e0_wp
        h2i4        = 0.777e0_wp
        h2i5        = -7.583
        h2i6        = 1.501
        rnaie       = 3.0353
        rnapol1     = 290.4e0_wp
        hpol        = 4.5e0_wp
        hie         = 13.598e0_wp
        rnaquad     = -16.8992e0_wp
        hnai1       = -193.341e0_wp
        hnai2       = 16.968e0_wp
        hnai3       = -77.120e0_wp
        hnai4       = -1.024e0_wp
        hnai5       = 6.009e0_wp
        hnai6       = -1.241e0_wp
        hnau1       = 119.880e0_wp
        hnau2       = -0.620e0_wp
        hnau3       = 31.335e0_wp
        hnau4       = 70.414e0_wp
        hnau5       = -13.016e0_wp
        hnau6       = 0.8681e0_wp
        hnaa1       = 0.67e0_wp
        hnaa2       = -92.229e0_wp
        hnaa3       = 51.01019e0_wp
        hnaa4       = -130.37826e0_wp
        hnaa5       = 6.40285e0_wp
        hnab1       = 0.34333e0_wp
        hnab2       = -92.229e0_wp
        hnab3       = -48.54867e0_wp
        hnab4       = -17.30523e0_wp
        hnab5       = 0.99821e0_wp
        rqq         = 9.766021100e4_wp
        rqh         = 2.355466062927e6_wp
        rb          = 9.3119383371714e5_wp
        balp        = 1.93701e0_wp
        brho        = 1.4e0_wp
        rd          = 9.1351724748364e8_wp
        dalp        = 1.43307e0_wp
        drho        = 5.71654e0_wp
        h2hc1       = 0.677e0_wp
        h2hc2       = 0.378e0_wp
        h2hc3       = 0.281e0_wp
        h2hc4       = -1.022e0_wp
        h2hc5       = 0.665e0_wp
        h2qc1       = 1.532e0_wp
        h2qc2       = -.04635e0_wp
        h2qc3       = 1.826e0_wp
        h2qc4       = -2.508e0_wp
        h2qc5       = 1.856e0_wp
        c2a         = 2.18504e0_wp
        c2b         = 0.10945e0_wp
        c2c         = 0.15039e0_wp
        deh2        = 4.74772265651191638e0_wp
        una2p       = 2.1037
        g_ehfid     = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        cevau  = 1.e0_wp / 27.2113961e0_wp
        onethd = 1.0e0_wp / 3.0e0_wp
        twothd = 2.0e0_wp / 3.0e0_wp
        do i = 1, nt
            ! this section calcluates useful quantities - the A to BC center
            ! of mass distance, the angles between A-BC and BC, the shorest
            ! of any of the three diatomic distances, and the bond angle
            ! H-Na-H
            d2_b = r(i,1) + r(i,1)
            do g_i_ = 1, g_p_
                g_r1s(g_i_) = d2_b * g_r(g_i_,i,1)
            enddo
            r1s  = r(i,1) * r(i,1)
            d2_b = r(i,2) + r(i,2)
            do g_i_ = 1, g_p_
                g_r2s(g_i_) = d2_b * g_r(g_i_,i,2)
            enddo
            r2s  = r(i,2) * r(i,2)
            d2_b = r(i,3) + r(i,3)
            do g_i_ = 1, g_p_
                g_r3s(g_i_) = d2_b * g_r(g_i_,i,3)
            enddo
            r3s  = r(i,3) * r(i,3)
            rnag = sqrt(abs(0.5e0_wp*(r1s+r3s-0.5e0_wp*r2s)))
            if( rnag .gt. 0.0e0_wp )then
                g_rnag(1) = 0.5e0_wp*r(i,1)/rnag
                g_rnag(2) = -0.25e0_wp*r(i,2)/rnag
                g_rnag(3) = 0.5e0_wp*r(i,3)/rnag
            else
                g_rnag(1) = 0.0e0_wp
                g_rnag(2) = 0.0e0_wp
                g_rnag(3) = 0.0e0_wp
            endif
            d2_b = rnag + rnag
            do g_i_ = 1, g_p_
                g_rnags(g_i_) = d2_b * g_rnag(g_i_)
            enddo
            rnags = rnag * rnag
            do g_i_ = 1, g_p_
                g_rnagc(g_i_) = rnags * g_rnag(g_i_) + rnag * g_rnags(g_i_)
            enddo
            rnagc = rnags * rnag
            if( rnag .eq. 0.0e0_wp )then
                do g_i_ = 1, g_p_
                    g_cstsq(g_i_) = 0.0e0_wp
                enddo
                cstsq = 0.0e0_wp
            else
                d7_v = r(i,2) * rnag
                d8_v = 0.5e0_wp * (r3s - r1s) / d7_v
                d9_v = d8_v * d8_v
                d1_p = 2.0e0_wp * d8_v
                d4_b = d1_p * ((-d8_v) / d7_v)
                d5_b = d4_b * rnag
                d6_b = d4_b * r(i,2)
                d7_b = d1_p * (1.0e0_wp / d7_v) * 0.5e0_wp
                do g_i_ = 1, g_p_
                    g_cstsq(g_i_) = d6_b * g_rnag(g_i_) + d5_b * g_r(g_i_,i,2) + (-d7_b) * g_r1s(g_i_) + d7_b * &
                        & g_r3s(g_i_)
                enddo
                cstsq = d9_v
            endif
            d7_v  = real(2, wp) * r(i,1)
            d9_v  = d7_v * r(i,3)
            d10_v = (r1s + r3s - r2s) / d9_v
            d2_b  = 1.0e0_wp / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,3) * real(2, wp)
            do g_i_ = 1, g_p_
                g_cosa(g_i_) = d5_b * g_r(g_i_,i,3) + d6_b * g_r(g_i_,i,1) + (-d2_b) * g_r2s(g_i_) + d2_b * &
                    & g_r3s(g_i_) + d2_b * g_r1s(g_i_)
            enddo
            cosa = d10_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-10.0e0_wp) * g_rnags(g_i_)
            enddo
            d1_w = (-10.0e0_wp) * rnags
            d4_v = 0.5e0_wp + 0.5e0_wp * cosa
            d5_v = cstsq * d4_v
            d7_v = exp(d1_w)
            d1_p = d7_v
            d8_v = 1.0e0_wp - d7_v
            d5_b = (-d5_v) * d1_p
            d6_b = d8_v * d4_v
            d9_b = d8_v * cstsq * 0.5e0_wp
            do g_i_ = 1, g_p_
                g_coscos(g_i_) = d5_b * g_d1_w(g_i_) + d9_b * g_cosa(g_i_) + d6_b * g_cstsq(g_i_)
            enddo
            coscos = d5_v * d8_v
            frac   = 0.958e0_wp
            rnag1  = sqrt(abs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))
            if( rnag1 .gt. 0.0e0_wp )then
                g_rnag1(1) = r(i,1)*frac*(frac-1.0e0_wp)/rnag1
                g_rnag1(2) = r(i,2)*(1.0e0_wp-frac)/rnag1
                g_rnag1(3) = r(i,3)*frac/rnag1
            else
                g_rnag1(1) = 0.0e0_wp
                g_rnag1(2) = 0.0e0_wp
                g_rnag1(3) = 0.0e0_wp
            endif
            if( rnag3 .gt. 0.0e0_wp )then
                g_rnag3(1) = r(i,1)*frac/rnag3
                g_rnag3(2) = r(i,2)*(1.0e0_wp-frac)/rnag3
                g_rnag3(3) = r(i,3)*frac*(frac-1.0e0_wp)/rnag3
            else
                g_rnag3(1) = 0.0e0_wp
                g_rnag3(2) = 0.0e0_wp
                g_rnag3(3) = 0.0e0_wp
            endif
            d2_v = rnag1 * rnag1
            d3_b = d2_v + rnag1 * rnag1 + rnag1 * rnag1
            do g_i_ = 1, g_p_
                g_rnag1c(g_i_) = d3_b * g_rnag1(g_i_)
            enddo
            rnag1c = d2_v * rnag1
            d2_v   = rnag3 * rnag3
            d3_b   = d2_v + rnag3 * rnag3 + rnag3 * rnag3
            do g_i_ = 1, g_p_
                g_rnag3c(g_i_) = d3_b * g_rnag3(g_i_)
            enddo
            rnag3c = d2_v * rnag3
            if( rnag1 .eq. 0.0e0_wp )then
                do g_i_ = 1, g_p_
                    g_cstsq1(g_i_) = 0.0e0_wp
                enddo
                cstsq1 = 1.0e0_wp
            else
                d11_v = r(i,1) * rnag1
                d12_v = (r3s - r2s - r1s + 2.0e0_wp * r3s * frac) / d11_v
                d13_v = d12_v * d12_v
                d1_p  = 2.0e0_wp * d12_v
                d3_b  = 0.25e0_wp * d1_p
                d4_b  = d3_b * (1.0e0_wp / d11_v)
                d5_b  = d3_b * ((-d12_v) / d11_v)
                d6_b  = d5_b * rnag1
                d7_b  = d5_b * r(i,1)
                d11_b = d4_b * frac * 2.0e0_wp + d4_b
                do g_i_ = 1, g_p_
                    g_cstsq1(g_i_) = d7_b * g_rnag1(g_i_) + d6_b * g_r(g_i_,i,1) + (-d4_b) * g_r1s(g_i_) + (- d4_b) * &
                        & g_r2s(g_i_) + d11_b * g_r3s(g_i_)
                enddo
                cstsq1 = 0.25e0_wp * d13_v
            endif
            if( rnag3 .eq. 0.0e0_wp )then
                do g_i_ = 1, g_p_
                    g_cstsq3(g_i_) = 0.0e0_wp
                enddo
                cstsq3 = 1.0e0_wp
            else
                d11_v = r(i,3) * rnag3
                d12_v = (r1s - r2s - r3s + 2.0e0_wp * r1s * frac) / d11_v
                d13_v = d12_v * d12_v
                d1_p  = 2.0e0_wp * d12_v
                d3_b  = 0.25e0_wp * d1_p
                d4_b  = d3_b * (1.0e0_wp / d11_v)
                d5_b  = d3_b * ((-d12_v) / d11_v)
                d6_b  = d5_b * rnag3
                d7_b  = d5_b * r(i,3)
                d11_b = d4_b * frac * 2.0e0_wp + d4_b
                do g_i_ = 1, g_p_
                    g_cstsq3(g_i_) = d7_b * g_rnag3(g_i_) + d6_b * g_r(g_i_,i,3) + (-d4_b) * g_r3s(g_i_) + (- d4_b) * &
                        & g_r2s(g_i_) + d11_b * g_r1s(g_i_)
                enddo
                cstsq3 = 0.25e0_wp * d13_v
            endif
            ! The first NaH singlet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s1c2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-s1c2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s1c6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-s1c6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s1c3 + s1c4 * r(i,1) + s1c5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * s1c5
            d7_b  = d12_v * d9_v + d8_b * s1c4
            d14_b = s1c1 * d2_p
            do g_i_ = 1, g_p_
                g_s1a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            s1a = s1c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s1ic2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-s1ic2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s1ic6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-s1ic6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s1ic3 + s1ic4 * r(i,1) + s1ic5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * s1ic5
            d7_b  = d12_v * d9_v + d8_b * s1ic4
            d14_b = s1ic1 * d2_p
            do g_i_ = 1, g_p_
                g_s1ic(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            s1ic = s1ic1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s1ialp * g_rnags(g_i_)
            enddo
            d1_w = s1ialp * rnags
            d3_v = s1ic - s1a
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d6_v = 1.0e0_wp - d5_v
            d7_b = (-d3_v) * d1_p
            d2_b = 1.0e0_wp + (-d6_v)
            do g_i_ = 1, g_p_
                g_s1i(g_i_) = d7_b * g_d1_w(g_i_) + d6_v * g_s1ic(g_i_) + d2_b * g_s1a(g_i_)
            enddo
            s1i  = s1a + d3_v * d6_v
            d4_v = 0.5e0_wp * (1.0e0_wp - cosa)
            d6_v = s1i - s1a
            d2_b = 1.0e0_wp + (-d4_v)
            d8_b = -(d6_v * 0.5e0_wp)
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d4_v * g_s1i(g_i_) + d8_b * g_cosa(g_i_) + d2_b * g_s1a(g_i_)
            enddo
            s1 = s1a + d4_v * d6_v
            ! The first NaH triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t10c2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-t10c2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t10c6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-t10c6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t10c3 + t10c4 * r(i,1) + t10c5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * t10c5
            d7_b  = d12_v * d9_v + d8_b * t10c4
            d14_b = t10c1 * d2_p
            do g_i_ = 1, g_p_
                g_t10(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            t10 = t10c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t190c2) * g_r(g_i_,i,1)
            enddo
            d1_w = (-t190c2) * r(i,1)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t190c6) * g_r(g_i_,i,1)
            enddo
            d2_w  = (-t190c6) * r(i,1)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t190c3 + t190c4 * r(i,1) + t190c5 * r1s
            d10_v = r(i,1) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,1)
            d11_b = d8_b * t190c5
            d7_b  = d12_v * d9_v + d8_b * t190c4
            d14_b = t190c1 * d2_p
            do g_i_ = 1, g_p_
                g_t190(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7_b * g_r(g_i_,i,1) + d14_b * g_d1_w(g_i_)
            enddo
            t190  = t190c1 * d2_v + d10_v * d12_v
            d4_v  = 0.5e0_wp - 0.5e0_wp * cosa
            d6_v  = t10 - t190
            d9_v  = t10 - t190
            d6_b  = coscos + d4_v
            d7_b  = -coscos + 1.0e0_wp + (-d4_v)
            d12_b = (-d6_v) * 0.5e0_wp
            do g_i_ = 1, g_p_
                g_t1(g_i_) = d9_v * g_coscos(g_i_) + d6_b * g_t10(g_i_) + d12_b * g_cosa(g_i_) + d7_b * g_t190(g_i_)
            enddo
            t1 = t190 + d4_v * d6_v + d9_v * coscos
            ! The H2 triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t2c2) * g_r(g_i_,i,2)
            enddo
            d1_w = (-t2c2) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t2c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-t2c6) * r(i,2)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t2c3 + t2c4 * r(i,2) + t2c5 * r2s
            d10_v = r(i,2) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,2)
            d11_b = d8_b * t2c5
            d7_b  = d12_v * d9_v + d8_b * t2c4
            d14_b = t2c1 * d2_p
            do g_i_ = 1, g_p_
                g_t2a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + d7_b * g_r(g_i_,i,2) + d14_b * g_d1_w(g_i_)
            enddo
            t2a = t2c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t290c2) * g_r(g_i_,i,2)
            enddo
            d1_w = (-t290c2) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t290c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-t290c6) * r(i,2)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t290c3 + r(i,2) * t290c4 + t290c5 * r2s
            d10_v = r(i,2) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,2)
            d11_b = d8_b * t290c5
            d7_b  = d12_v * d9_v + d8_b * t290c4
            d14_b = t290c1 * d2_p
            do g_i_ = 1, g_p_
                g_t290c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + d7_b * g_r(g_i_,i,2) + d14_b * g_d1_w(g_i_)
            enddo
            t290c = t290c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t290alp * g_rnags(g_i_)
            enddo
            d1_w = t290alp * (rnags - t290rho)
            d3_v = t290c - t2a
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d6_v = 1.0e0_wp - d5_v
            d5_b = 0.5e0_wp * d6_v
            d8_b = (-(0.5e0_wp * d3_v)) * d1_p
            d2_b = 1.0e0_wp + (-d5_b)
            do g_i_ = 1, g_p_
                g_t290(g_i_) = d8_b * g_d1_w(g_i_) + d5_b * g_t290c(g_i_) + d2_b * g_t2a(g_i_)
            enddo
            t290 = t2a + d3_v * d6_v * 0.5e0_wp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t20c2) * g_r(g_i_,i,2)
            enddo
            d1_w = (-t20c2) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t20c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-t20c6) * r(i,2)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t20c3 + t20c4 * r(i,2) + t20c5 * r2s
            d10_v = r(i,2) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,2)
            d11_b = d8_b * t20c5
            d7_b  = d12_v * d9_v + d8_b * t20c4
            d14_b = t20c1 * d2_p
            do g_i_ = 1, g_p_
                g_t20c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + d7_b * g_r(g_i_,i,2) + d14_b * g_d1_w(g_i_)
            enddo
            t20c = t20c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t20alp * g_rnags(g_i_)
            enddo
            d1_w = t20alp * (rnags - t20rho)
            d3_v = t20c - t2a
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d6_v = 1.0e0_wp - d5_v
            d5_b = 0.5e0_wp * d6_v
            d8_b = (-(0.5e0_wp * d3_v)) * d1_p
            d2_b = 1.0e0_wp + (-d5_b)
            do g_i_ = 1, g_p_
                g_t20(g_i_) = d8_b * g_d1_w(g_i_) + d5_b * g_t20c(g_i_) + d2_b * g_t2a(g_i_)
            enddo
            t20  = t2a + d3_v * d6_v * 0.5e0_wp
            d3_v = t20 - t290
            d2_b = 1.0e0_wp + (-coscos)
            do g_i_ = 1, g_p_
                g_t2(g_i_) = d3_v * g_coscos(g_i_) + coscos * g_t20(g_i_) + d2_b * g_t290(g_i_)
            enddo
            t2 = t290 + d3_v * coscos
            ! The H2 singlet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s2c3) * g_r(g_i_,i,2)
            enddo
            d1_w = (-s2c3) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s2c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-s2c6) * r(i,2)
            d3_v  = s2c1 + s2c2 * r(i,2)
            d5_v  = exp(d1_w)
            d2_p  = d5_v
            d9_v  = s2c4 + s2c5 * r(i,2)
            d10_v = r2s * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d7_b  = d10_v * d1_p
            d8_b  = d12_v * d9_v
            d14_b = d3_v * d2_p
            d11_b = d12_v * r2s * s2c5 + d5_v * s2c2
            do g_i_ = 1, g_p_
                g_s2a(g_i_) = d7_b * g_d2_w(g_i_) + d8_b * g_r2s(g_i_) + d14_b * g_d1_w(g_i_) + d11_b * g_r(g_i_,i,2)
            enddo
            s2a  = d3_v * d5_v + d10_v * d12_v + una2p
            d2_b = 1.0e0_wp / 0.1e0_wp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-d2_b) * g_t2(g_i_) + d2_b * g_s2a(g_i_)
            enddo
            d1_w = (s2a - t2) / 0.1e0_wp
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d6_v = s2a - t2
            d6_b = (-0.5e0_wp) * d5_v
            d9_b = (-0.5e0_wp) * d6_v * d1_p
            d7_b = d6_b + 0.5e0_wp
            d8_b = -d6_b + 0.5e0_wp
            do g_i_ = 1, g_p_
                g_s2mod(g_i_) = d9_b * g_d1_w(g_i_) + d8_b * g_t2(g_i_) + d7_b * g_s2a(g_i_)
            enddo
            s2mod = 0.5e0_wp * (s2a + t2 - d5_v * d6_v)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s290c3) * g_r(g_i_,i,2)
            enddo
            d1_w = (-s290c3) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s290c6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-s290c6) * r(i,2)
            d3_v  = s290c1 + s290c2 * r(i,2)
            d5_v  = exp(d1_w)
            d2_p  = d5_v
            d9_v  = s290c4 + s290c5 * r(i,2)
            d10_v = r2s * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d7_b  = d12_v * d9_v
            d13_b = d3_v * d2_p
            d10_b = d12_v * r2s * s290c5 + d5_v * s290c2
            do g_i_ = 1, g_p_
                g_s290c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r2s(g_i_) + d13_b * g_d1_w(g_i_) + d10_b * g_r(g_i_,i,2)
            enddo
            s290c = d3_v * d5_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s290alp * g_rnags(g_i_)
            enddo
            d1_w = s290alp * (rnags - s290rho)
            d4_v = (s290c - s2mod) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.0e0_wp - d6_v
            d7_b = (-d4_v) * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_s290(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s290c(g_i_) + d2_b * g_s2mod(g_i_)
            enddo
            s290 = s2mod + d4_v * d7_v
            d3_v = s2mod - s290
            d2_b = 1.0e0_wp + (-coscos)
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d3_v * g_coscos(g_i_) + coscos * g_s2mod(g_i_) + d2_b * g_s290(g_i_)
            enddo
            s2 = s290 + d3_v * coscos
            ! The second NaH singlet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s1c2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-s1c2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s1c6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-s1c6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s1c3 + s1c4 * r(i,3) + s1c5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * s1c5
            d7_b  = d12_v * d9_v + d8_b * s1c4
            d14_b = s1c1 * d2_p
            do g_i_ = 1, g_p_
                g_s3a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            s3a = s1c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s1ic2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-s1ic2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-s1ic6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-s1ic6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = s1ic3 + s1ic4 * r(i,3) + s1ic5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * s1ic5
            d7_b  = d12_v * d9_v + d8_b * s1ic4
            d14_b = s1ic1 * d2_p
            do g_i_ = 1, g_p_
                g_s3ic(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            s3ic = s1ic1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s1ialp * g_rnags(g_i_)
            enddo
            d1_w = s1ialp * rnags
            d3_v = s3ic - s3a
            d5_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d5_v *  d5_v)
            d6_v = 1.0e0_wp - d5_v
            d7_b = (-d3_v) * d1_p
            d2_b = 1.0e0_wp + (-d6_v)
            do g_i_ = 1, g_p_
                g_s3i(g_i_) = d7_b * g_d1_w(g_i_) + d6_v * g_s3ic(g_i_) + d2_b * g_s3a(g_i_)
            enddo
            s3i  = s3a + d3_v * d6_v
            d4_v = 0.5e0_wp * (1.0e0_wp - cosa)
            d6_v = s3i - s3a
            d2_b = 1.0e0_wp + (-d4_v)
            d8_b = -(d6_v * 0.5e0_wp)
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d4_v * g_s3i(g_i_) + d8_b * g_cosa(g_i_) + d2_b * g_s3a(g_i_)
            enddo
            s3 = s3a + d4_v * d6_v
            ! The second NaH triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t10c2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-t10c2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t10c6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-t10c6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t10c3 + t10c4 * r(i,3) + t10c5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * t10c5
            d7_b  = d12_v * d9_v + d8_b * t10c4
            d14_b = t10c1 * d2_p
            do g_i_ = 1, g_p_
                g_t30(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            t30 = t10c1 * d2_v + d10_v * d12_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-t190c2) * g_r(g_i_,i,3)
            enddo
            d1_w = (-t190c2) * r(i,3)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-t190c6) * g_r(g_i_,i,3)
            enddo
            d2_w  = (-t190c6) * r(i,3)
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = t190c3 + t190c4 * r(i,3) + t190c5 * r3s
            d10_v = r(i,3) * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d6_b  = d10_v * d1_p
            d8_b  = d12_v * r(i,3)
            d11_b = d8_b * t190c5
            d7_b  = d12_v * d9_v + d8_b * t190c4
            d14_b = t190c1 * d2_p
            do g_i_ = 1, g_p_
                g_t390(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7_b * g_r(g_i_,i,3) + d14_b * g_d1_w(g_i_)
            enddo
            t390  = t190c1 * d2_v + d10_v * d12_v
            d4_v  = 0.5e0_wp - 0.5e0_wp * cosa
            d6_v  = t30 - t390
            d9_v  = t30 - t390
            d6_b  = coscos + d4_v
            d7_b  = -coscos + 1.0e0_wp + (-d4_v)
            d12_b = (-d6_v) * 0.5e0_wp
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d9_v * g_coscos(g_i_) + d6_b * g_t30(g_i_) + d12_b * g_cosa(g_i_) + d7_b * g_t390(g_i_)
            enddo
            t3 = t390 + d4_v * d6_v + d9_v * coscos
            ! this long section contains data pertinent to long range forces
            !  H2 ionization energy:
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-h2i3) * g_r(g_i_,i,2)
            enddo
            d1_w = (-h2i3) * r(i,2)
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-h2i6) * g_r(g_i_,i,2)
            enddo
            d2_w  = (-h2i6) * r(i,2)
            d3_v  = h2i1 + h2i2 * r(i,2)
            d5_v  = exp(d1_w)
            d2_p  = d5_v
            d9_v  = h2i4 + h2i5 * r(i,2)
            d10_v = r2s * d9_v
            d12_v = exp(d2_w)
            d1_p  = d12_v
            d9_b  = d10_v * d1_p
            d10_b = d12_v * d9_v
            d16_b = d3_v * d2_p
            d13_b = d12_v * r2s * h2i5 + d5_v * h2i2
            do g_i_ = 1, g_p_
                g_h2ie(g_i_) = -g_s2a(g_i_) + d9_b * g_d2_w(g_i_) + d10_b * g_r2s(g_i_) + d16_b * g_d1_w(g_i_) + &
                    & d13_b * g_r(g_i_,i,2)
            enddo
            h2ie = d3_v * d5_v + d10_v * d12_v - s2a + hie
            !  H2 polarizability:
            !   al1 is the parallel polarizability for h2
            !   al2 is the perpendicaular polarizability for h2
            !   alt is the anisotropic? polarizability
            !   alb is the mean polarizability
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-h2a1) * g_r(g_i_,i,2)
            enddo
            d1_w  = h2a1 * (h2a2 - r(i,2))
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = h2a6 * r(i,2)
            d12_v = h2a7 * r2s
            d15_v = h2a8 * r2s
            d16_v = d15_v * r2s
            d19_v = h2a9 * r2s
            d20_v = d19_v * r2s
            d22_v = h2a3 + h2a4 * r(i,2) + h2a5 * r2s + d9_v * r2s + d12_v * r2s + d16_v * r(i,2) + d20_v * r2s
            d7_b  = d2_v * r2s
            d12_b = d2_v * r(i,2)
            d8_b = d2_v * d20_v + d7_b * d19_v + d7_b * r2s * h2a9 + d12_b * d15_v + d12_b * r2s * h2a8 + d2_v * &
                & d12_v + d2_v * r2s * h2a7 + d2_v * d9_v + d2_v * h2a5
            d13_b = d2_v * d16_v + d2_v * r2s * h2a6 + d2_v * h2a4
            d24_b = d22_v * d1_p
            do g_i_ = 1, g_p_
                g_al1(g_i_) = d8_b * g_r2s(g_i_) + d13_b * g_r(g_i_,i,2) + d24_b * g_d1_w(g_i_)
            enddo
            al1 = d2_v * d22_v + 2.0e0_wp * hpol
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-h2b1) * g_r(g_i_,i,2)
            enddo
            d1_w  = h2b1 * (h2b2 - r(i,2))
            d2_v  = exp(d1_w)
            d2_p  = d2_v
            d9_v  = r(i,2) ** ( 3 - 2)
            d9_v  = d9_v * r(i,2)
            d1_p  = 3 *  d9_v
            d9_v  = d9_v * r(i,2)
            d12_v = h2b7 * r2s
            d15_v = h2b8 * r2s
            d16_v = d15_v * r2s
            d19_v = h2b9 * r2s
            d20_v = d19_v * r2s
            d22_v = h2b3 + h2b4 * r(i,2) + h2b5 * r2s + h2b6 * d9_v + d12_v * r2s + d16_v * r(i,2) + d20_v * r2s
            d7_b  = d2_v * r2s
            d12_b = d2_v * r(i,2)
            d8_b = d2_v * d20_v + d7_b * d19_v + d7_b * r2s * h2b9 + d12_b * d15_v + d12_b * r2s * h2b8 + d2_v * &
                & d12_v + d2_v * r2s * h2b7 + d2_v * h2b5
            d13_b = d2_v * d16_v + d2_v * h2b6 * d1_p + d2_v * h2b4
            d24_b = d22_v * d2_p
            do g_i_ = 1, g_p_
                g_al2(g_i_) = d8_b * g_r2s(g_i_) + d13_b * g_r(g_i_,i,2) + d24_b * g_d1_w(g_i_)
            enddo
            al2 = d2_v * d22_v + 2.0e0_wp * hpol
            do g_i_ = 1, g_p_
                g_alt(g_i_) = (-twothd) * g_al2(g_i_) + twothd * g_al1(g_i_)
            enddo
            alt = twothd * (al1 - al2)
            do g_i_ = 1, g_p_
                g_alb(g_i_) = twothd * g_al2(g_i_) + onethd * g_al1(g_i_)
            enddo
            alb  = onethd * al1 + twothd * al2
            d5_v = 1.5e0_wp * cstsq - 0.5e0_wp
            d7_b = alt * 1.5e0_wp
            do g_i_ = 1, g_p_
                g_alhh(g_i_) = d7_b * g_cstsq(g_i_) + d5_v * g_alt(g_i_) + g_alb(g_i_)
            enddo
            alhh = alb + alt * d5_v
            ! here follow properties necesary for the long range forces in
            ! the Na + H2 channel:
            ! H2 properties:
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-h2hc1) * g_r(g_i_,i,2)
            enddo
            d1_w  = (-h2hc1) * r(i,2)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d8_v  = h2hc4 * r(i,2)
            d11_v = h2hc5 * r2s
            d13_v = h2hc2 * r(i,2) + h2hc3 * r2s + d8_v * r2s + d11_v * r2s
            d7_b  = d2_v * d11_v + d2_v * r2s * h2hc5 + d2_v * d8_v + d2_v * h2hc3
            d11_b = d2_v * r2s * h2hc4 + d2_v * h2hc2
            d14_b = d13_v * d1_p
            do g_i_ = 1, g_p_
                g_h2hex(g_i_) = d7_b * g_r2s(g_i_) + d11_b * g_r(g_i_,i,2) + d14_b * g_d1_w(g_i_)
            enddo
            h2hex = d2_v * d13_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-h2qc1) * g_r(g_i_,i,2)
            enddo
            d1_w  = (-h2qc1) * r(i,2)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d8_v  = h2qc4 * r(i,2)
            d11_v = h2qc5 * r2s
            d13_v = h2qc2 * r(i,2) + h2qc3 * r2s + d8_v * r2s + d11_v * r2s
            d7_b  = d2_v * d11_v + d2_v * r2s * h2qc5 + d2_v * d8_v + d2_v * h2qc3
            d11_b = d2_v * r2s * h2qc4 + d2_v * h2qc2
            d14_b = d13_v * d1_p
            do g_i_ = 1, g_p_
                g_h2quad(g_i_) = d7_b * g_r2s(g_i_) + d11_b * g_r(g_i_,i,2 ) + d14_b * g_d1_w(g_i_)
            enddo
            h2quad = d2_v * d13_v
            ! NaH properties
            d9_v  = hnai6 * r(i,1)
            d11_v = hnai3 + hnai4 * r(i,1) + hnai5 * r1s + d9_v * r1s
            d12_v = (hnai1 + hnai2 * r(i,1)) / d11_v
            d5_b  = (-d12_v) / d11_v
            d9_b  = d5_b * d9_v + d5_b * hnai5
            d10_b = d5_b * r1s * hnai6 + d5_b * hnai4 + 1.0e0_wp / d11_v * hnai2
            do g_i_ = 1, g_p_
                g_rnahie1(g_i_) = d9_b * g_r1s(g_i_) + d10_b * g_r(g_i_,i,1)
            enddo
            rnahie1 = d12_v + rnaie + una2p
            d9_v    = hnai6 * r(i,3)
            d11_v   = hnai3 + hnai4 * r(i,3) + hnai5 * r3s + d9_v * r3s
            d12_v   = (hnai1 + hnai2 * r(i,3)) / d11_v
            d5_b    = (-d12_v) / d11_v
            d9_b    = d5_b * d9_v + d5_b * hnai5
            d10_b   = d5_b * r3s * hnai6 + d5_b * hnai4 + 1.0e0_wp / d11_v * hnai2
            do g_i_ = 1, g_p_
                g_rnahie3(g_i_) = d9_b * g_r3s(g_i_) + d10_b * g_r(g_i_,i,3)
            enddo
            rnahie3 = d12_v + rnaie + una2p
            d10_v   = hnau6 * r(i,1)
            d12_v   = hnau3 + hnau4 * r(i,1) + hnau5 * r1s + d10_v * r1s
            d13_v   = (hnau1 * r(i,1) + hnau2 * r1s) / d12_v
            d2_b    = 1.0e0_wp / d12_v
            d3_b    = (-d13_v) / d12_v
            d7_b    = d3_b * d10_v + d3_b * hnau5 + d2_b * hnau2
            d8_b    = d3_b * r1s * hnau6 + d3_b * hnau4 + d2_b * hnau1
            do g_i_ = 1, g_p_
                g_rnahmu1(g_i_) = d7_b * g_r1s(g_i_) + d8_b * g_r(g_i_,i,1 )
            enddo
            rnahmu1 = d13_v
            d10_v   = hnau6 * r(i,3)
            d12_v   = hnau3 + hnau4 * r(i,3) + hnau5 * r3s + d10_v * r3s
            d13_v   = (hnau1 * r(i,3) + hnau2 * r3s) / d12_v
            d2_b    = 1.0e0_wp / d12_v
            d3_b    = (-d13_v) / d12_v
            d7_b    = d3_b * d10_v + d3_b * hnau5 + d2_b * hnau2
            d8_b    = d3_b * r3s * hnau6 + d3_b * hnau4 + d2_b * hnau1
            do g_i_ = 1, g_p_
                g_rnahmu3(g_i_) = d7_b * g_r3s(g_i_) + d8_b * g_r(g_i_,i,3 )
            enddo
            rnahmu3 = d13_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-hnaa1) * g_r(g_i_,i,1)
            enddo
            d1_w  = (-hnaa1) * r(i,1)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = hnaa5 * r1s
            d11_v = hnaa2 + hnaa3 * r(i,1) + hnaa4 * r1s + d9_v * r(i,1)
            d9_b  = d2_v * r(i,1) * hnaa5 + d2_v * hnaa4
            d8_b  = d2_v * d9_v + d2_v * hnaa3
            d13_b = d11_v * d1_p
            do g_i_ = 1, g_p_
                g_rnahpar1(g_i_) = d9_b * g_r1s(g_i_) + d8_b * g_r(g_i_,i,1) + d13_b * g_d1_w(g_i_)
            enddo
            rnahpar1 = 159.267e0_wp + hpol + d2_v * d11_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-hnaa1) * g_r(g_i_,i,3)
            enddo
            d1_w  = (-hnaa1) * r(i,3)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = hnaa5 * r3s
            d11_v = hnaa2 + hnaa3 * r(i,3) + hnaa4 * r3s + d9_v * r(i,3)
            d9_b  = d2_v * r(i,3) * hnaa5 + d2_v * hnaa4
            d8_b  = d2_v * d9_v + d2_v * hnaa3
            d13_b = d11_v * d1_p
            do g_i_ = 1, g_p_
                g_rnahpar3(g_i_) = d9_b * g_r3s(g_i_) + d8_b * g_r(g_i_,i,3) + d13_b * g_d1_w(g_i_)
            enddo
            rnahpar3 = 159.267e0_wp + hpol + d2_v * d11_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-hnab1) * g_r(g_i_,i,1)
            enddo
            d1_w  = (-hnab1) * r(i,1)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = hnab5 * r1s
            d11_v = hnab2 + hnab3 * r(i,1) + hnab4 * r1s + d9_v * r(i,1)
            d9_b  = d2_v * r(i,1) * hnab5 + d2_v * hnab4
            d8_b  = d2_v * d9_v + d2_v * hnab3
            d13_b = d11_v * d1_p
            do g_i_ = 1, g_p_
                g_rnahperp1(g_i_) = d9_b * g_r1s(g_i_) + d8_b * g_r(g_i_,i,1) + d13_b * g_d1_w(g_i_)
            enddo
            rnahperp1 = 159.267e0_wp + hpol + d2_v * d11_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-hnab1) * g_r(g_i_,i,3)
            enddo
            d1_w  = (-hnab1) * r(i,3)
            d2_v  = exp(d1_w)
            d1_p  = d2_v
            d9_v  = hnab5 * r3s
            d11_v = hnab2 + hnab3 * r(i,3) + hnab4 * r3s + d9_v * r(i,3)
            d9_b  = d2_v * r(i,3) * hnab5 + d2_v * hnab4
            d8_b  = d2_v * d9_v + d2_v * hnab3
            d13_b = d11_v * d1_p
            do g_i_ = 1, g_p_
                g_rnahperp3(g_i_) = d9_b * g_r3s(g_i_) + d8_b * g_r(g_i_,i,3) + d13_b * g_d1_w(g_i_)
            enddo
            rnahperp3 = 159.267e0_wp + hpol + d2_v * d11_v
            d3_v      = rnahpar1 - rnahperp1
            d2_b      = 1.0e0_wp + (-cstsq1)
            do g_i_ = 1, g_p_
                g_rnahpol1(g_i_) = d3_v * g_cstsq1(g_i_) + cstsq1 * g_rnahpar1(g_i_) + d2_b * g_rnahperp1(g_i_)
            enddo
            rnahpol1 = rnahperp1 + d3_v * cstsq1
            do g_i_ = 1, g_p_
                g_rnahpol1b(g_i_) = twothd * g_rnahperp1(g_i_) + onethd * g_rnahpar1(g_i_)
            enddo
            rnahpol1b = onethd * rnahpar1 + twothd * rnahperp1
            d3_v      = rnahpar3 - rnahperp3
            d2_b      = 1.0e0_wp + (-cstsq3)
            do g_i_ = 1, g_p_
                g_rnahpol3(g_i_) = d3_v * g_cstsq3(g_i_) + cstsq3 * g_rnahpar3(g_i_) + d2_b * g_rnahperp3(g_i_)
            enddo
            rnahpol3 = rnahperp3 + d3_v * cstsq3
            do g_i_ = 1, g_p_
                g_rnahpol3b(g_i_) = twothd * g_rnahperp3(g_i_) + onethd * g_rnahpar3(g_i_)
            enddo
            rnahpol3b = onethd * rnahpar3 + twothd * rnahperp3
            ! electrostatic forces
            ! quadrupole quadrupole
            do g_i_ = 1, g_p_
                g_qqint(g_i_) = (-7.0e0_wp) * g_cstsq(g_i_)
            enddo
            qqint = 3.0e0_wp - 7.0e0_wp * cstsq
            d5_v  = rnagc * rnag
            d7_v  = d5_v * rnag + rqq
            d8_v  = 0.75e0_wp * rnaquad * h2quad / d7_v
            d6_b  = qqint * ((-d8_v) / d7_v)
            d7_b  = d6_b * rnag
            d9_b  = d7_b * rnag
            d8_b  = d6_b * d5_v + d7_b * rnagc
            d10_b = qqint * (1.0e0_wp / d7_v) * (0.75e0_wp * rnaquad)
            do g_i_ = 1, g_p_
                g_qq(g_i_) = d8_v * g_qqint(g_i_) + d8_b * g_rnag(g_i_) + d9_b * g_rnagc(g_i_) + d10_b * g_h2quad(g_i_)
            enddo
            qq   = d8_v * qqint
            d2_v = (-39.0625e0_wp) * rnagc
            d4_b = rnag * (-39.0625e0_wp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d2_v * g_rnag(g_i_) + d4_b * g_rnagc(g_i_)
            enddo
            d1_w = d2_v * rnag
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_dexp2(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            dexp2 = d2_v
            d3_v  = 3.0e0_wp * qqint
            d4_v  = (-4.0e0_wp) * pi / d3_v
            d5_v  = d4_v - 1.0e0_wp
            d8_v  = 1.0e0_wp + d5_v * dexp2
            d6_b  = qq * d5_v
            d9_b  = qq * dexp2 * ((-d4_v) / d3_v) * 3.0e0_wp
            do g_i_ = 1, g_p_
                g_qq(g_i_) = d6_b * g_dexp2(g_i_) + d9_b * g_qqint(g_i_) + d8_v * g_qq(g_i_)
            enddo
            qq = qq * d8_v
            ! quadrupole hexadecapole
            d2_v = 147.0e0_wp * cstsq
            d4_b = 162.0e0_wp + (-d2_v) + (-cstsq) * 147.0e0_wp
            do g_i_ = 1, g_p_
                g_qhint(g_i_) = d4_b * g_cstsq(g_i_)
            enddo
            qhint = -45.0e0_wp - d2_v * cstsq + 162.0e0_wp * cstsq
            d4_v  = rnagc * rnagc
            d7_v  = d4_v * rnag + rqh
            d8_v  = 0.3125e0_wp * rnaquad * h2hex / d7_v
            d6_b  = qhint * ((-d8_v) / d7_v)
            d7_b  = d6_b * rnag
            d8_b  = d6_b * d4_v
            d9_b  = d7_b * rnagc + d7_b * rnagc
            d10_b = qhint * (1.0e0_wp / d7_v) * (0.3125e0_wp * rnaquad)
            do g_i_ = 1, g_p_
                g_qh(g_i_) = d8_v * g_qhint(g_i_) + d8_b * g_rnag(g_i_) + d9_b * g_rnagc(g_i_) + d10_b * g_h2hex(g_i_)
            enddo
            qh   = d8_v * qhint
            d3_v = (-40.8e0_wp) * pi / qhint
            d4_v = d3_v - 1.0e0_wp
            d7_v = 1.0e0_wp + d4_v * dexp2
            d6_b = qh * d4_v
            d8_b = qh * dexp2 * ((-d3_v) / qhint)
            do g_i_ = 1, g_p_
                g_qh(g_i_) = d6_b * g_dexp2(g_i_) + d8_b * g_qhint(g_i_) + d7_v * g_qh(g_i_)
            enddo
            qh = qh * d7_v
            ! induction forces
            ! dipole induced dipole
            d2_v = rnag1c * rnag1c
            d4_v = d2_v * rnag1c + rd
            d5_v = 1.0e0_wp / d4_v
            d3_b = (-d5_v) / d4_v
            d4_b = d3_b * rnag1c
            d5_b = d3_b * d2_v + d4_b * rnag1c + d4_b * rnag1c
            do g_i_ = 1, g_p_
                g_div1(g_i_) = d5_b * g_rnag1c(g_i_)
            enddo
            div1  = d5_v
            d2_v  = rnahmu1 * rnahmu1
            d1_p  = 2.0e0_wp * rnahmu1
            d4_v  = (-d2_v) * hpol
            d7_v  = 1.5e0_wp * cstsq1 + 0.5e0_wp
            d8_v  = d4_v * d7_v
            d10_v = d8_v * rnag1c
            d4_b  = div1 * rnag1c
            d5_b  = div1 * d8_v
            d9_b  = d4_b * d4_v * 1.5e0_wp
            d12_b = (-(d4_b * d7_v * hpol)) * d1_p
            do g_i_ = 1, g_p_
                g_rmuimu1(g_i_) = d10_v * g_div1(g_i_) + d5_b * g_rnag1c(g_i_) + d9_b * g_cstsq1(g_i_) + d12_b * &
                    & g_rnahmu1(g_i_)
            enddo
            rmuimu1 = d10_v * div1
            d2_v    = rnag3c * rnag3c
            d4_v    = d2_v * rnag3c + rd
            d5_v    = 1.0e0_wp / d4_v
            d3_b    = (-d5_v) / d4_v
            d4_b    = d3_b * rnag3c
            d5_b    = d3_b * d2_v + d4_b * rnag3c + d4_b * rnag3c
            do g_i_ = 1, g_p_
                g_div3(g_i_) = d5_b * g_rnag3c(g_i_)
            enddo
            div3  = d5_v
            d2_v  = rnahmu3 * rnahmu3
            d1_p  = 2.0e0_wp * rnahmu3
            d4_v  = (-d2_v) * hpol
            d7_v  = 1.5e0_wp * cstsq3 + 0.5e0_wp
            d8_v  = d4_v * d7_v
            d10_v = d8_v * rnag3c
            d4_b  = div3 * rnag3c
            d5_b  = div3 * d8_v
            d9_b  = d4_b * d4_v * 1.5e0_wp
            d12_b = (-(d4_b * d7_v * hpol)) * d1_p
            do g_i_ = 1, g_p_
                g_rmuimu3(g_i_) = d10_v * g_div3(g_i_) + d5_b * g_rnag3c(g_i_) + d9_b * g_cstsq3(g_i_) + d12_b * &
                    & g_rnahmu3(g_i_)
            enddo
            rmuimu3 = d10_v * div3
            ! dispersion forces
            d3_v  = rnaie + h2ie
            d4_v  = (-1.5e0_wp) * rnaie * h2ie / d3_v
            d5_v  = d4_v * rnapol1
            d10_v = rnagc * rnagc + rb
            d11_v = d5_v * alhh / d10_v
            d2_b  = 1.0e0_wp / d10_v
            d4_b  = (-d11_v) / d10_v
            d5_b  = d4_b * rnagc + d4_b * rnagc
            d7_b  = d2_b * d5_v
            d8_b  = d2_b * alhh * rnapol1
            d11_b = d8_b * ((-d4_v) / d3_v) + d8_b * (1.0e0_wp / d3_v) * ((-1.5e0_wp) * rnaie)
            do g_i_ = 1, g_p_
                g_c6(g_i_) = d5_b * g_rnagc(g_i_) + d7_b * g_alhh(g_i_) + d11_b * g_h2ie(g_i_)
            enddo
            c6   = d11_v
            d4_v = alb / alhh
            d5_v = d4_v - 1.0e0_wp
            d8_v = 1.0e0_wp + d5_v * dexp2
            d6_b = c6 * d5_v
            d7_b = c6 * dexp2
            d8_b = d7_b * (1.0e0_wp / alhh)
            d9_b = d7_b * ((-d4_v) / alhh)
            do g_i_ = 1, g_p_
                g_c6(g_i_) = d6_b * g_dexp2(g_i_) + d9_b * g_alhh(g_i_) + d8_b * g_alb(g_i_) + d8_v * g_c6(g_i_)
            enddo
            c6    = c6 * d8_v
            d4_v  = rnahie1 + hie
            d5_v  = (-1.5e0_wp) * rnahie1 * hie / d4_v
            d8_v  = d5_v * rnahpol1 * hpol
            d10_v = d8_v * rnag1c
            d5_b  = div1 * d8_v
            d6_b  = div1 * rnag1c * hpol
            d7_b  = d6_b * rnahpol1
            d8_b  = d6_b * d5_v
            d11_b = d7_b * ((-d5_v) / d4_v) + d7_b * (1.0e0_wp / d4_v) * hie * (-1.5e0_wp)
            do g_i_ = 1, g_p_
                g_c61(g_i_) = d10_v * g_div1(g_i_) + d5_b * g_rnag1c(g_i_) + d8_b * g_rnahpol1(g_i_) + d11_b * &
                    & g_rnahie1(g_i_)
            enddo
            c61  = d10_v * div1
            d2_v = (-39.0625e0_wp) * rnag1c
            d4_b = rnag1 * (-39.0625e0_wp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d2_v * g_rnag1(g_i_) + d4_b * g_rnag1c(g_i_)
            enddo
            d1_w = d2_v * rnag1
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_dexp1(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            dexp1 = d2_v
            d4_v  = rnahpol1b / rnahpol1
            d5_v  = d4_v - 1.0e0_wp
            d8_v  = 1.0e0_wp + d5_v * dexp1
            d6_b  = c61 * d5_v
            d7_b  = c61 * dexp1
            d8_b  = d7_b * (1.0e0_wp / rnahpol1)
            d9_b  = d7_b * ((-d4_v) / rnahpol1)
            do g_i_ = 1, g_p_
                g_c61(g_i_) = d6_b * g_dexp1(g_i_) + d9_b * g_rnahpol1(g_i_) + d8_b * g_rnahpol1b(g_i_) + d8_v * &
                    & g_c61(g_i_)
            enddo
            c61   = c61 * d8_v
            d4_v  = rnahie3 + hie
            d5_v  = (-1.5e0_wp) * rnahie3 * hie / d4_v
            d8_v  = d5_v * rnahpol3 * hpol
            d10_v = d8_v * rnag3c
            d5_b  = div3 * d8_v
            d6_b  = div3 * rnag3c * hpol
            d7_b  = d6_b * rnahpol3
            d8_b  = d6_b * d5_v
            d11_b = d7_b * ((-d5_v) / d4_v) + d7_b * (1.0e0_wp / d4_v) * hie * (-1.5e0_wp)
            do g_i_ = 1, g_p_
                g_c63(g_i_) = d10_v * g_div3(g_i_) + d5_b * g_rnag3c(g_i_) + d8_b * g_rnahpol3(g_i_) + d11_b * &
                    & g_rnahie3(g_i_)
            enddo
            c63  = d10_v * div3
            d2_v = (-39.0625e0_wp) * rnag3c
            d4_b = rnag3 * (-39.0625e0_wp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d2_v * g_rnag3(g_i_) + d4_b * g_rnag3c(g_i_)
            enddo
            d1_w = d2_v * rnag3
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_dexp3(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            dexp3 = d2_v
            d4_v  = rnahpol3b / rnahpol3
            d5_v  = d4_v - 1.0e0_wp
            d8_v  = 1.0e0_wp + d5_v * dexp3
            d6_b  = c63 * d5_v
            d7_b  = c63 * dexp3
            d8_b  = d7_b * (1.0e0_wp / rnahpol3)
            d9_b  = d7_b * ((-d4_v) / rnahpol3)
            do g_i_ = 1, g_p_
                g_c63(g_i_) = d6_b * g_dexp3(g_i_) + d9_b * g_rnahpol3(g_i_) + d8_b * g_rnahpol3b(g_i_) + d8_v * &
                    & g_c63(g_i_)
            enddo
            c63   = c63 * d8_v
            d4_v  = 1.5e0_wp * cstsq1 + 0.5e0_wp
            d5_v  = 1.0e0_wp / d4_v
            d6_v  = d5_v - 1.0e0_wp
            d9_v  = 1.0e0_wp + d6_v * dexp1
            d6_b  = rmuimu1 * d6_v
            d10_b = rmuimu1 * dexp1 * ((-d5_v) / d4_v) * 1.5e0_wp
            do g_i_ = 1, g_p_
                g_cind1(g_i_) = d6_b * g_dexp1(g_i_) + d10_b * g_cstsq1(g_i_ ) + d9_v * g_rmuimu1(g_i_)
            enddo
            cind1 = rmuimu1 * d9_v
            d4_v  = 1.5e0_wp * cstsq3 + 0.5e0_wp
            d5_v  = 1.0e0_wp / d4_v
            d6_v  = d5_v - 1.0e0_wp
            d9_v  = 1.0e0_wp + d6_v * dexp3
            d6_b  = rmuimu3 * d6_v
            d10_b = rmuimu3 * dexp3 * ((-d5_v) / d4_v) * 1.5e0_wp
            do g_i_ = 1, g_p_
                g_cind3(g_i_) = d6_b * g_dexp3(g_i_) + d10_b * g_cstsq3(g_i_ ) + d9_v * g_rmuimu3(g_i_)
            enddo
            cind3 = rmuimu3 * d9_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = balp * g_r(g_i_,i,2)
            enddo
            d1_w  = balp * (r(i,2) - brho)
            d7_v  = (qq + qh + c6 * cevau) * 0.5e0_wp
            d9_v  = tanh (d1_w)
            d1_p  = 1.0e0_wp - ( d9_v *  d9_v)
            d10_v = 1.0e0_wp - d9_v
            d5_b  = (-d7_v) * d1_p
            d6_b  = d10_v * 0.5e0_wp
            d9_b  = d6_b * cevau
            do g_i_ = 1, g_p_
                g_elec(g_i_) = d5_b * g_d1_w(g_i_) + d9_b * g_c6(g_i_) + d6_b * g_qh(g_i_) + d6_b * g_qq(g_i_)
            enddo
            elec = d7_v * d10_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = dalp * g_r(g_i_,i,1)
            enddo
            d1_w = dalp * (r(i,1) - drho)
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_c6swt1(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            c6swt1 = 0.5e0_wp * (1.0e0_wp - d2_v)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = dalp * g_r(g_i_,i,3)
            enddo
            d1_w = dalp * (r(i,3) - drho)
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_c6swt3(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            c6swt3 = 0.5e0_wp * (1.0e0_wp - d2_v)
            d3_b   = cevau * c6swt1
            d4_b   = cevau * c61
            do g_i_ = 1, g_p_
                g_c61(g_i_) = d4_b * g_c6swt1(g_i_) + d3_b * g_c61(g_i_)
            enddo
            c61  = c61 * c6swt1 * cevau
            d3_b = cevau * c6swt3
            d4_b = cevau * c63
            do g_i_ = 1, g_p_
                g_c63(g_i_) = d4_b * g_c6swt3(g_i_) + d3_b * g_c63(g_i_)
            enddo
            c63 = c63 * c6swt3 * cevau
            do g_i_ = 1, g_p_
                g_cind1(g_i_) = cind1 * g_c6swt1(g_i_) + c6swt1 * g_cind1(g_i_)
            enddo
            cind1 = cind1 * c6swt1
            do g_i_ = 1, g_p_
                g_cind3(g_i_) = cind3 * g_c6swt3(g_i_) + c6swt3 * g_cind3(g_i_)
            enddo
            cind3 = cind3 * c6swt3
            do g_i_ = 1, g_p_
                g_elec(g_i_) = g_cind3(g_i_) + g_cind1(g_i_) + g_c63(g_i_) + g_c61(g_i_) + g_elec(g_i_)
            enddo
            elec = elec + c61 + c63 + cind1 + cind3
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
            d8_v = sqrt(d1_w)
            if( d1_w .gt. 0.0e0_wp )then
                d1_p = 1.0e0_wp / (2.0e0_wp *  d8_v)
            else
                d1_p = 0.0e0_wp
            endif
            d5_b = (-(1.0e0_wp / sqrt(2.0e0_wp))) * d1_p
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d5_b * g_d1_w(g_i_) + g_coul3(g_i_) + g_coul2 (g_i_) + g_coul1(g_i_)
            enddo
            e(i) = coul1 + coul2 + coul3 + deh2 - d8_v / sqrt(2.0e0_wp)
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = g_elec(g_i_) + cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau + elec
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_6_2000_g_pot22_legacy_state_a1118117

!> @brief Compatibility interface for g_pot22.
subroutine g_pot22( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_nah2_leps_6_2000_g_pot22_legacy_state_a1118117, only: dispatch
    use chempotpy_kinds,                                          only: wp
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
