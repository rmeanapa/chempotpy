!   System:                        NaH2
!   Common name:                   NaH2 surface set 7sd
!   Interface:                     3-2V
!   Number of electronic surfaces: 2
!   Number of derivatives:         1
!   Note:               This version does not contain long range forces,
!                       and has derivatives.  It is based on the surface
!                       6 fit.
!   References:         Unpublished, M. D. Hack and D. G. Truhlar (2000)
!   Protocol: POT evaluates the stateless diabatic surfaces.
!               of the energy with respect to the coordinates for a give
!               geometric configuration
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      Ground diabatic electronic state, first-excited diabatic electron
!and the potential coupling.
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Na
!      infinitely far from the H2 diatomic and R(H2) set equal to the
!      H2 equilibrium diatomic value.
!   Parameters:
!   Coordinates:
!      Internal, Definition: R(1) = R(Na-H)
!                            R(2) = R(H-H)
!                            R(3) = R(Na-H)
! This is NaH2 potential matrix 7sd, finalized on 5/17/99
! This version contains analytic derivatives.
!  all surfaces.  The potential is called by calling POT with the parame
!  r, e, nt, nsurf.  r(nt,3) and ei(nt) are double precision arrays of
!  dimension nt which is an integer representing the number of geometrie
!  input.  The integer nsurf selects the desired surface.  nsurf=1 is th
!  lower diabatic surface, nsurf=2 is the coupling surface, nsurf=3 is
!  the upper diabatic surface, nsurf = 4 selects the lower adiabatic
!  surface, and nsurf = 5 selects the upper adiabatic surface

!> @brief Implement NaH2_LEPS_7SD_2000, a diabatic potential matrix for NaH2.
module chempotpy_nah2_leps_7sd_2000_dpem_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NaH2_LEPS_7SD_2000 diabatic potential matrix for NaH2.
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
end module chempotpy_nah2_leps_7sd_2000_dpem_pes

!> @brief Expose NaH2_LEPS_7SD_2000 through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_nah2_leps_7sd_2000_dpem_pes, only: evaluate
    use chempotpy_kinds,                       only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: u(nstates,nstates)
    real(wp), intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem

!> @brief Evaluate EvdRdX for the NaH2_LEPS_7SD_2000_DPEM potential.
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

!> @brief Evaluate diagonalize for the NaH2_LEPS_7SD_2000_DPEM potential.
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

!> @brief Evaluate pot for the NaH2_LEPS_7SD_2000_DPEM potential.
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
    endif
end subroutine pot

!  U11 surface - lower diabatic surface
!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_nah2_leps_7sd_2000_dpem_g_pot11_legacy_state_d8f976ce
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
        real(wp) :: b0, bect, bect2, beh2, bes1, bf, cevau, cos1, cos2, cosa, cosg, dect, deh2, deh2m, des1, des20, &
            & e(:), gamh2a, gamh2b, onethd, r(:,:), r1s, r2s, r3s, rac2, rave, rect, reh2, res1, rmid, s1, s2, s2alp, s2rho, s3, &
            & sswt, twothd, vmid, xect, xeh2, xes1, xes3
        integer :: i, nt
        ! NaH singlet
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d2_p, d7_b, d12_b, d9_b, d13_v, d12_v, d11_v, d2_v, d1_p, d2_b
        real(wp)           :: d3_v, d4_v, d3_b, d4_b, d5_v, d6_v, d7_v, d1_w, &
            & d9_v, d10_v
        real(wp) :: d5_b, d6_b, g_r1s(g_pmax_), g_r(:,:,:), &
            & g_r2s(g_pmax_), g_r3s(g_pmax_), g_rmid(g_pmax_), g_rave(g_pmax_), &
            & g_cos1(g_pmax_), g_cos2(g_pmax_)
        real(wp) :: g_cosg(g_pmax_), g_cosa(g_pmax_), g_d1_w(g_pmax_), g_xes1(g_pmax_), g_s1(g_pmax_), g_beh2(g_pmax_), g_sswt(g_pmax_), g_deh2m(g_pmax_), g_xeh2(g_pmax_), g_s2(g_pmax_)
        real(wp) :: g_xes3(g_pmax_), g_s3(g_pmax_), g_e(:,:), &
            & g_xect(g_pmax_), g_vmid(g_pmax_)
        integer :: g_ehfid
        ! H2 singlet
        ! H-Na-H term
        !  misc constants
        b0      = 0.e0_wp
        bect    = 0.e0_wp
        bect2   = 0.e0_wp
        beh2    = 0.e0_wp
        bes1    = 0.e0_wp
        bf      = 0.e0_wp
        cevau   = 0.e0_wp
        cos1    = 0.e0_wp
        cos2    = 0.e0_wp
        cosa    = 0.e0_wp
        cosg    = 0.e0_wp
        d10_v   = 0.e0_wp
        d11_v   = 0.e0_wp
        d12_b   = 0.e0_wp
        d12_v   = 0.e0_wp
        d13_v   = 0.e0_wp
        d1_p    = 0.e0_wp
        d1_w    = 0.e0_wp
        d2_b    = 0.e0_wp
        d2_p    = 0.e0_wp
        d2_v    = 0.e0_wp
        d3_b    = 0.e0_wp
        d3_v    = 0.e0_wp
        d4_b    = 0.e0_wp
        d4_v    = 0.e0_wp
        d5_b    = 0.e0_wp
        d5_v    = 0.e0_wp
        d6_b    = 0.e0_wp
        d6_v    = 0.e0_wp
        d7_b    = 0.e0_wp
        d7_v    = 0.e0_wp
        d9_b    = 0.e0_wp
        d9_v    = 0.e0_wp
        dect    = 0.e0_wp
        deh2    = 0.e0_wp
        deh2m   = 0.e0_wp
        des1    = 0.e0_wp
        des20   = 0.e0_wp
        g_beh2  = 0.e0_wp
        g_cos1  = 0.e0_wp
        g_cos2  = 0.e0_wp
        g_cosa  = 0.e0_wp
        g_cosg  = 0.e0_wp
        g_d1_w  = 0.e0_wp
        g_deh2m = 0.e0_wp
        g_r1s   = 0.e0_wp
        g_r2s   = 0.e0_wp
        g_r3s   = 0.e0_wp
        g_rave  = 0.e0_wp
        g_rmid  = 0.e0_wp
        g_s1    = 0.e0_wp
        g_s2    = 0.e0_wp
        g_s3    = 0.e0_wp
        g_sswt  = 0.e0_wp
        g_vmid  = 0.e0_wp
        g_xect  = 0.e0_wp
        g_xeh2  = 0.e0_wp
        g_xes1  = 0.e0_wp
        g_xes3  = 0.e0_wp
        gamh2a  = 0.e0_wp
        gamh2b  = 0.e0_wp
        onethd  = 0.e0_wp
        r1s     = 0.e0_wp
        r2s     = 0.e0_wp
        r3s     = 0.e0_wp
        rac2    = 0.e0_wp
        rave    = 0.e0_wp
        rect    = 0.e0_wp
        reh2    = 0.e0_wp
        res1    = 0.e0_wp
        rmid    = 0.e0_wp
        s1      = 0.e0_wp
        s2      = 0.e0_wp
        s2alp   = 0.e0_wp
        s2rho   = 0.e0_wp
        s3      = 0.e0_wp
        sswt    = 0.e0_wp
        twothd  = 0.e0_wp
        vmid    = 0.e0_wp
        xect    = 0.e0_wp
        xeh2    = 0.e0_wp
        xes1    = 0.e0_wp
        xes3    = 0.e0_wp
        des1    = 0.55118e0_wp
        bes1    = 1.05118e0_wp
        res1    = 3.38583e0_wp
        reh2    = 1.401e0_wp
        bf      = 1.627e0_wp
        b0      = 1.194e0_wp
        gamh2a  = 2.323e0_wp
        gamh2b  = 3.126e0_wp
        deh2    = 4.74806e0_wp
        des20   = 0.996e0_wp
        s2alp   = 3.0e0_wp
        s2rho   = 7.0e0_wp
        dect    = 2.24409e0_wp
        bect    = 0.57874e0_wp
        rect    = 3.46457e0_wp
        bect2   = 0.87143e0_wp
        g_ehfid = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        rac2   = 1.e0_wp / sqrt(2.e0_wp)
        cevau  = 1.e0_wp / 27.2113961e0_wp
        onethd = 1.0e0_wp / 3.0e0_wp
        twothd = 2.0e0_wp / 3.0e0_wp
        do i = 1, nt
            ! useful quantities calculated first
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
            r3s = r(i,3) * r(i,3)
            do g_i_ = 1, g_p_
                g_rmid(g_i_) = 0.5e0_wp * g_r(g_i_,i,3) + 0.5e0_wp * g_r(g_i_,i,1)
            enddo
            rmid = 0.5e0_wp * (r(i,1) + r(i,3))
            do g_i_ = 1, g_p_
                g_rave(g_i_) = (-0.5e0_wp) * g_r(g_i_,i,3) + 0.5e0_wp * g_r(g_i_,i,1)
            enddo
            rave = 0.5e0_wp * (r(i,1) - r(i,3))
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            d7_v  = 2.0e0_wp * r(i,1)
            d9_v  = d7_v * r(i,2)
            d10_v = (r1s + r2s - r3s) / d9_v
            d2_b  = 1.0e0_wp / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0e0_wp
            do g_i_ = 1, g_p_
                g_cos1(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,1) + (-d2_b) * g_r3s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r1s(g_i_)
            enddo
            cos1  = d10_v
            d7_v  = 2.0e0_wp * r(i,3)
            d9_v  = d7_v * r(i,2)
            d10_v = (r3s + r2s - r1s) / d9_v
            d2_b  = 1.0e0_wp / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0e0_wp
            do g_i_ = 1, g_p_
                g_cos2(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,3) + (-d2_b) * g_r1s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r3s(g_i_)
            enddo
            cos2 = d10_v
            do g_i_ = 1, g_p_
                g_cosg(g_i_) = (-0.5e0_wp) * g_cos2(g_i_) + 0.5e0_wp * g_cos1(g_i_ )
            enddo
            cosg  = 0.5e0_wp * (cos1 - cos2)
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
            !==============================================================k
            !This section contains the diatomic curves
            ! The first NaH singlet
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
            d2_v = des1 * onethd * xes1
            d3_v = xes1 + 2.0e0_wp
            d4_b = d2_v + d3_v * (des1 * onethd)
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d4_b * g_xes1(g_i_)
            enddo
            s1 = d2_v * d3_v
            ! The H2 singlet
            reh2   = 1.401e0_wp
            bf     = 1.627e0_wp
            b0     = 1.194e0_wp
            gamh2a = 2.323e0_wp
            gamh2b = 3.126e0_wp
            d5_v   = (r(i,2) / gamh2b) * (r(i,2) / gamh2b)
            d2_p   = 2.0e0_wp * (r(i,2) / gamh2b)
            d11_v  = (r(i,2) / gamh2b) * (r(i,2) / gamh2b)
            d1_p   = 2.0e0_wp * (r(i,2) / gamh2b)
            d12_v  = bf - r(i,2) / gamh2a + d11_v
            d13_v  = bf * (b0 - r(i,2) / gamh2a + d5_v) / d12_v
            d3_b   = (-d13_v) / d12_v
            d9_b   = 1.0e0_wp / d12_v * bf
            d7_b  = d3_b * d1_p * (1.0e0_wp / gamh2b) + (-d3_b) * (1.0e0_wp / gamh2a) + d9_b * d2_p * (1.0e0_wp / &
                & gamh2b) + (- d9_b) * (1.0e0_wp / gamh2a)
            do g_i_ = 1, g_p_
                g_beh2(g_i_) = d7_b * g_r(g_i_,i,2)
            enddo
            beh2 = d13_v
            d3_b = 1.0e0_wp / s2alp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_rmid(g_i_)
            enddo
            d1_w = (rmid - s2rho) / s2alp
            d2_v = tanh (d1_w)
            d2_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_v = 0.5e0_wp * (1.0e0_wp - d2_v)
            d6_v = cosg * cosg
            d1_p = 2.0e0_wp * cosg
            d4_b = d4_v * d1_p
            d7_b = (-(d6_v * 0.5e0_wp)) * d2_p
            do g_i_ = 1, g_p_
                g_sswt(g_i_) = d4_b * g_cosg(g_i_) + d7_b * g_d1_w(g_i_)
            enddo
            sswt = d4_v * d6_v
            d3_b = deh2 * (des20 - 1.0e0_wp)
            do g_i_ = 1, g_p_
                g_deh2m(g_i_) = d3_b * g_sswt(g_i_)
            enddo
            deh2m = deh2 + deh2 * (des20 - 1.0e0_wp) * sswt
            d4_v  = r(i,2) - reh2
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beh2) * g_r(g_i_,i,2) + (-d4_v) * g_beh2( g_i_)
            enddo
            d1_w = (-beh2) * d4_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xeh2(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xeh2 = d2_v
            d3_v = deh2m * xeh2
            d4_v = xeh2 - 2.0e0_wp
            d5_b = d4_v * xeh2
            d4_b = d3_v + d4_v * deh2m
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d4_b * g_xeh2(g_i_) + d5_b * g_deh2m(g_i_)
            enddo
            s2 = d3_v * d4_v
            ! The second NaH singlet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bes1) * g_r(g_i_,i,3)
            enddo
            d1_w = (-bes1) * (r(i,3) - res1)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xes3(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xes3 = d2_v
            d2_v = des1 * onethd * xes3
            d3_v = xes3 + 2.0e0_wp
            d4_b = d2_v + d3_v * (des1 * onethd)
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d4_b * g_xes3(g_i_)
            enddo
            s3 = d2_v * d3_v
            ! here is where it all gets put together . . .
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = g_s3(g_i_) + g_s2(g_i_) + g_s1(g_i_)
            enddo
            e(i) = s1 + s2 + s3 + deh2
            ! here we will try a post-LEPS correction term
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bect) * g_rmid(g_i_)
            enddo
            d1_w = (-bect) * (rmid - rect)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xect(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xect = d2_v
            d2_v = rave * rave
            d1_p = 2.0e0_wp * rave
            d3_b = (-bect2) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_rave(g_i_)
            enddo
            d1_w  = (-bect2) * d2_v
            d2_v  = dect * xect
            d3_v  = xect - 2.0e0_wp
            d4_v  = d2_v * d3_v
            d6_v  = exp(d1_w)
            d2_p  = d6_v
            d7_v  = d4_v * d6_v
            d11_v = (0.5e0_wp * (1.0e0_wp - cosa)) * (0.5e0_wp * (1.0e0_wp - cosa))
            d1_p  = 2.0e0_wp * (0.5e0_wp * (1.0e0_wp - cosa))
            d6_b  = -(d7_v * d1_p * 0.5e0_wp)
            d7_b  = d11_v * d6_v
            d9_b  = d11_v * d4_v * d2_p
            d12_b = d7_b * d2_v + d7_b * d3_v * dect
            do g_i_ = 1, g_p_
                g_vmid(g_i_) = d6_b * g_cosa(g_i_) + d9_b * g_d1_w(g_i_) + d12_b * g_xect(g_i_)
            enddo
            vmid = d7_v * d11_v
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = g_vmid(g_i_) + g_e(g_i_,i)
            enddo
            e(i) = e(i) + vmid
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_7sd_2000_dpem_g_pot11_legacy_state_d8f976ce

!> @brief Compatibility interface for g_pot11.
subroutine g_pot11( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_nah2_leps_7sd_2000_dpem_g_pot11_legacy_state_d8f976ce, only: dispatch
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

! U12 coupling surface . . .
!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_nah2_leps_7sd_2000_dpem_g_pot12_legacy_state_d8f976ce
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
        real(wp)           :: ap, alp, bet, rnagsq, rnag, cosg, rd0, r20
        real(wp)           :: r1s, r2s, r3s, cos1, cos2, theta
        real(wp)           :: gc, gc2
        real(wp)           :: cevau
        integer            :: nt, i, j, k
        logical            :: leven
        real(wp)           :: r(:,:), e(:)
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d2_p, d2_w, d1_w, d4_b, d5_b, d10_v, d9_v, d2_v, d8_b, d2_b
        real(wp)           :: d1_p, d3_v, d4_v, d5_v, d7_b, d7_v, d6_b, d3_b, &
            & g_r1s(g_pmax_), g_r(:,:,:)
        real(wp) :: g_r2s(g_pmax_), g_r3s(g_pmax_), g_rnagsq(g_pmax_), g_rnag(g_pmax_), g_cos1(g_pmax_), g_cos2(g_pmax_), g_cosg(g_pmax_), g_gc(g_pmax_), g_gc2(g_pmax_), g_d1_w(g_pmax_)
        real(wp) :: g_d2_w(g_pmax_), g_e(:,:)
        integer  :: g_ehfid
        alp      = 0.e0_wp
        ap       = 0.e0_wp
        bet      = 0.e0_wp
        cevau    = 0.e0_wp
        cos1     = 0.e0_wp
        cos2     = 0.e0_wp
        cosg     = 0.e0_wp
        d10_v    = 0.e0_wp
        d1_p     = 0.e0_wp
        d1_w     = 0.e0_wp
        d2_b     = 0.e0_wp
        d2_p     = 0.e0_wp
        d2_v     = 0.e0_wp
        d2_w     = 0.e0_wp
        d3_b     = 0.e0_wp
        d3_v     = 0.e0_wp
        d4_b     = 0.e0_wp
        d4_v     = 0.e0_wp
        d5_b     = 0.e0_wp
        d5_v     = 0.e0_wp
        d6_b     = 0.e0_wp
        d7_b     = 0.e0_wp
        d7_v     = 0.e0_wp
        d8_b     = 0.e0_wp
        d9_v     = 0.e0_wp
        g_cos1   = 0.e0_wp
        g_cos2   = 0.e0_wp
        g_cosg   = 0.e0_wp
        g_d1_w   = 0.e0_wp
        g_d2_w   = 0.e0_wp
        g_gc     = 0.e0_wp
        g_gc2    = 0.e0_wp
        g_r1s    = 0.e0_wp
        g_r2s    = 0.e0_wp
        g_r3s    = 0.e0_wp
        g_rnag   = 0.e0_wp
        g_rnagsq = 0.e0_wp
        gc       = 0.e0_wp
        gc2      = 0.e0_wp
        r1s      = 0.e0_wp
        r20      = 0.e0_wp
        r2s      = 0.e0_wp
        r3s      = 0.e0_wp
        rd0      = 0.e0_wp
        rnag     = 0.e0_wp
        rnagsq   = 0.e0_wp
        theta    = 0.e0_wp
        ap       = 0.835e0_wp
        alp      = 4.0e0_wp
        rd0      = 4.0e0_wp
        bet      = 0.5e0_wp
        r20      = -1.3e0_wp
        g_ehfid  = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        theta = 15.6e0_wp * acos(-1.0e0_wp) / 180.0e0_wp
        cevau = 1.e0_wp / 27.2113961e0_wp
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
            do g_i_ = 1, g_p_
                g_rnagsq(g_i_) = (-0.25e0_wp) * g_r2s(g_i_) + 0.5e0_wp * g_r3s(g_i_) + 0.5e0_wp * g_r1s(g_i_)
            enddo
            rnagsq = 0.5e0_wp * r1s + 0.5e0_wp * r3s - 0.25e0_wp * r2s
            if( rnagsq .gt. 0.0e0_wp )then
                d2_v = sqrt(rnagsq)
                if( rnagsq .gt. 0.0e0_wp )then
                    d1_p = 1.0e0_wp / (2.0e0_wp *  d2_v)
                else
                    d1_p = 0.0e0_wp
                endif
                do g_i_ = 1, g_p_
                    g_rnag(g_i_) = d1_p * g_rnagsq(g_i_)
                enddo
                rnag = d2_v
            else
                do g_i_ = 1, g_p_
                    g_rnag(g_i_) = 0.0e0_wp
                enddo
                rnag = 0.0e0_wp
            endif
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            d7_v  = 2.0e0_wp * r(i,1)
            d9_v  = d7_v * r(i,2)
            d10_v = (r1s + r2s - r3s) / d9_v
            d2_b  = 1.0e0_wp / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0e0_wp
            do g_i_ = 1, g_p_
                g_cos1(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,1) + (-d2_b) * g_r3s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r1s(g_i_)
            enddo
            cos1  = d10_v
            d7_v  = 2.0e0_wp * r(i,3)
            d9_v  = d7_v * r(i,2)
            d10_v = (r3s + r2s - r1s) / d9_v
            d2_b  = 1.0e0_wp / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0e0_wp
            do g_i_ = 1, g_p_
                g_cos2(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,3) + (-d2_b) * g_r1s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r3s(g_i_)
            enddo
            cos2 = d10_v
            do g_i_ = 1, g_p_
                g_cosg(g_i_) = (-0.5e0_wp) * g_cos2(g_i_) + 0.5e0_wp * g_cos1(g_i_ )
            enddo
            cosg = 0.5e0_wp * (cos1 - cos2)
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
            d3_v = (gc - rd0) ** ( 4 - 2)
            d3_v = d3_v * (gc - rd0)
            d1_p = 4 *  d3_v
            d3_v = d3_v * (gc - rd0)
            d4_b = (-((1.0e0_wp / alp) ** 4)) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_gc(g_i_)
            enddo
            d1_w = (-((1.0e0_wp / alp) ** 4)) * d3_v
            d3_v = (gc2 - r20) * (gc2 - r20)
            d1_p = 2.0e0_wp * (gc2 - r20)
            d4_b = (-((1.0e0_wp / bet) ** 2)) * d1_p
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = d4_b * g_gc2(g_i_)
            enddo
            d2_w = (-((1.0e0_wp / bet) ** 2)) * d3_v
            d2_v = ap * cosg
            d4_v = exp(d1_w)
            d2_p = d4_v
            d5_v = d2_v * d4_v
            d7_v = exp(d2_w)
            d1_p = d7_v
            d4_b = d5_v * d1_p
            d7_b = d7_v * d2_v * d2_p
            d8_b = d7_v * d4_v * ap
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d4_b * g_d2_w(g_i_) + d7_b * g_d1_w(g_i_) + d8_b * g_cosg(g_i_)
            enddo
            e(i) = d5_v * d7_v
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
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
end module chempotpy_nah2_leps_7sd_2000_dpem_g_pot12_legacy_state_d8f976ce

!> @brief Compatibility interface for g_pot12.
subroutine g_pot12( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_nah2_leps_7sd_2000_dpem_g_pot12_legacy_state_d8f976ce, only: dispatch
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
module chempotpy_nah2_leps_7sd_2000_dpem_g_pot22_legacy_state_d8f976ce
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
        real(wp) :: b0, beh2, bet1, bet10, bet10b, bet190, bet1c, bet20, bet290, beta, bf, c2a, c2b, c2c, cevau, cos1, &
            & cos2, cosa, cosg, coul1, coul2, coul3, cplg2, deh2, deh2m, denah, des290, det1, det190, det1c, det20, det290, &
            & e(:), exch1, exch2, exch3, f10, g10, g30, gamh2a, gamh2b, gamnah, onethd, r(:,:), r1s, r2s, r3s, rave, reh2, &
            & renah, ret1, ret10, ret190, ret1c, ret20, ret290, rmid, s1, s2, s20alp, s20b, s20rho, s2a, s2alp, s2rho, s3, sswt, &
            & t1, t10, t190, t1alp, t1c, t1m, t1rho, t2, t20, t290, t2alp, t2rho, t3, t30, t390, t3c, t3m, twothd, una2p, w, &
            & xeh2, xenah, xet190, xet1c, xet20, xet290, xet3, xet390, xet3c
        integer :: i, nt
        ! NaH triplets
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d3_p, d1_w, d12_b, d11_b, d13_v, d12_v, d11_v, &
            & d2_v, d3_v, d4_v
        real(wp) :: d2_w, d2_b, d3_b, d4_b, d5_v, d6_v, d7_v, d8_v, &
            & d9_v, d10_v
        real(wp) :: d5_b, d6_b, d7_b, d8_b, d9_b, d10_b, d1_p, d2_p, g_rmid(g_pmax_), g_r(:,:,:)
        real(wp) :: g_r1s(g_pmax_), g_r2s(g_pmax_), g_r3s(g_pmax_), &
            & g_cos1(g_pmax_), g_cos2(g_pmax_), g_cosg(g_pmax_), g_beta(g_pmax_ &
            & ), g_d1_w(g_pmax_), g_xenah(g_pmax_), g_s1(g_pmax_)
        real(wp) :: g_g10(g_pmax_), g_d2_w(g_pmax_), g_t10(g_pmax_), g_xet190(g_pmax_), g_t190(g_pmax_), g_xet1c(g_pmax_), g_t1c(g_pmax_), g_t1m(g_pmax_), g_t1(g_pmax_), g_xet20(g_pmax_)
        real(wp) :: g_t20(g_pmax_), g_xet290(g_pmax_), g_t290(g_pmax_), g_t2(g_pmax_), g_beh2(g_pmax_), g_sswt(g_pmax_), g_deh2m(g_pmax_), g_xeh2(g_pmax_), g_s2a(g_pmax_), g_s2(g_pmax_)
        real(wp) :: g_s20b(g_pmax_), g_s3(g_pmax_), g_xet3(g_pmax_), g_t30(g_pmax_), g_g30(g_pmax_), g_xet390(g_pmax_), g_t390(g_pmax_), g_xet3c(g_pmax_), g_t3c(g_pmax_), g_t3m(g_pmax_)
        real(wp) :: g_t3(g_pmax_), g_coul1(g_pmax_), g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_), g_exch3( &
            & g_pmax_), g_w(g_pmax_), g_cplg2(g_pmax_), g_e(:,:)
        integer :: g_ehfid
        ! H2 triplets
        ! H2 singlets
        ! misc constants
        b0       = 0.e0_wp
        beh2     = 0.e0_wp
        bet1     = 0.e0_wp
        bet10    = 0.e0_wp
        bet10b   = 0.e0_wp
        bet190   = 0.e0_wp
        bet1c    = 0.e0_wp
        bet20    = 0.e0_wp
        bet290   = 0.e0_wp
        beta     = 0.e0_wp
        bf       = 0.e0_wp
        c2a      = 0.e0_wp
        c2b      = 0.e0_wp
        c2c      = 0.e0_wp
        cevau    = 0.e0_wp
        cos1     = 0.e0_wp
        cos2     = 0.e0_wp
        cosa     = 0.e0_wp
        cosg     = 0.e0_wp
        coul1    = 0.e0_wp
        coul2    = 0.e0_wp
        coul3    = 0.e0_wp
        cplg2    = 0.e0_wp
        d10_b    = 0.e0_wp
        d10_v    = 0.e0_wp
        d11_b    = 0.e0_wp
        d11_v    = 0.e0_wp
        d12_b    = 0.e0_wp
        d12_v    = 0.e0_wp
        d13_v    = 0.e0_wp
        d1_p     = 0.e0_wp
        d1_w     = 0.e0_wp
        d2_b     = 0.e0_wp
        d2_p     = 0.e0_wp
        d2_v     = 0.e0_wp
        d2_w     = 0.e0_wp
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
        d9_b     = 0.e0_wp
        d9_v     = 0.e0_wp
        deh2     = 0.e0_wp
        deh2m    = 0.e0_wp
        denah    = 0.e0_wp
        des290   = 0.e0_wp
        det1     = 0.e0_wp
        det190   = 0.e0_wp
        det1c    = 0.e0_wp
        det20    = 0.e0_wp
        det290   = 0.e0_wp
        exch1    = 0.e0_wp
        exch2    = 0.e0_wp
        exch3    = 0.e0_wp
        f10      = 0.e0_wp
        g10      = 0.e0_wp
        g30      = 0.e0_wp
        g_beh2   = 0.e0_wp
        g_beta   = 0.e0_wp
        g_cos1   = 0.e0_wp
        g_cos2   = 0.e0_wp
        g_cosg   = 0.e0_wp
        g_coul1  = 0.e0_wp
        g_coul2  = 0.e0_wp
        g_coul3  = 0.e0_wp
        g_cplg2  = 0.e0_wp
        g_d1_w   = 0.e0_wp
        g_d2_w   = 0.e0_wp
        g_deh2m  = 0.e0_wp
        g_exch1  = 0.e0_wp
        g_exch2  = 0.e0_wp
        g_exch3  = 0.e0_wp
        g_g10    = 0.e0_wp
        g_g30    = 0.e0_wp
        g_r1s    = 0.e0_wp
        g_r2s    = 0.e0_wp
        g_r3s    = 0.e0_wp
        g_rmid   = 0.e0_wp
        g_s1     = 0.e0_wp
        g_s2     = 0.e0_wp
        g_s20b   = 0.e0_wp
        g_s2a    = 0.e0_wp
        g_s3     = 0.e0_wp
        g_sswt   = 0.e0_wp
        g_t1     = 0.e0_wp
        g_t10    = 0.e0_wp
        g_t190   = 0.e0_wp
        g_t1c    = 0.e0_wp
        g_t1m    = 0.e0_wp
        g_t2     = 0.e0_wp
        g_t20    = 0.e0_wp
        g_t290   = 0.e0_wp
        g_t3     = 0.e0_wp
        g_t30    = 0.e0_wp
        g_t390   = 0.e0_wp
        g_t3c    = 0.e0_wp
        g_t3m    = 0.e0_wp
        g_w      = 0.e0_wp
        g_xeh2   = 0.e0_wp
        g_xenah  = 0.e0_wp
        g_xet190 = 0.e0_wp
        g_xet1c  = 0.e0_wp
        g_xet20  = 0.e0_wp
        g_xet290 = 0.e0_wp
        g_xet3   = 0.e0_wp
        g_xet390 = 0.e0_wp
        g_xet3c  = 0.e0_wp
        gamh2a   = 0.e0_wp
        gamh2b   = 0.e0_wp
        gamnah   = 0.e0_wp
        onethd   = 0.e0_wp
        r1s      = 0.e0_wp
        r2s      = 0.e0_wp
        r3s      = 0.e0_wp
        rave     = 0.e0_wp
        reh2     = 0.e0_wp
        renah    = 0.e0_wp
        ret1     = 0.e0_wp
        ret10    = 0.e0_wp
        ret190   = 0.e0_wp
        ret1c    = 0.e0_wp
        ret20    = 0.e0_wp
        ret290   = 0.e0_wp
        rmid     = 0.e0_wp
        s1       = 0.e0_wp
        s2       = 0.e0_wp
        s20alp   = 0.e0_wp
        s20b     = 0.e0_wp
        s20rho   = 0.e0_wp
        s2a      = 0.e0_wp
        s2alp    = 0.e0_wp
        s2rho    = 0.e0_wp
        s3       = 0.e0_wp
        sswt     = 0.e0_wp
        t1       = 0.e0_wp
        t10      = 0.e0_wp
        t190     = 0.e0_wp
        t1alp    = 0.e0_wp
        t1c      = 0.e0_wp
        t1m      = 0.e0_wp
        t1rho    = 0.e0_wp
        t2       = 0.e0_wp
        t20      = 0.e0_wp
        t290     = 0.e0_wp
        t2alp    = 0.e0_wp
        t2rho    = 0.e0_wp
        t3       = 0.e0_wp
        t30      = 0.e0_wp
        t390     = 0.e0_wp
        t3c      = 0.e0_wp
        t3m      = 0.e0_wp
        twothd   = 0.e0_wp
        una2p    = 0.e0_wp
        w        = 0.e0_wp
        xeh2     = 0.e0_wp
        xenah    = 0.e0_wp
        xet190   = 0.e0_wp
        xet1c    = 0.e0_wp
        xet20    = 0.e0_wp
        xet290   = 0.e0_wp
        xet3     = 0.e0_wp
        xet390   = 0.e0_wp
        xet3c    = 0.e0_wp
        ret10    = 3.22835e0_wp
        bet10    = 0.57874e0_wp
        bet10b   = 2.44094e0_wp
        f10      = 9.44882e0_wp
        det190   = 4.32441e0_wp
        bet190   = 2.43701e0_wp
        ret190   = 1.95480e0_wp
        t1alp    = 0.30778e0_wp
        t1rho    = 3.21850e0_wp
        det1c    = 4.92126e0_wp
        bet1c    = 0.65748e0_wp
        ret1c    = 2.71654e0_wp
        det20    = 3.29921e0_wp
        bet20    = 2.64567e0_wp
        ret20    = 2.44882e0_wp
        det290   = 3.97638e0_wp
        bet290   = 3.73937e0_wp
        ret290   = 0.96024e0_wp
        t2alp    = 1.0e0_wp
        t2rho    = 9.0e0_wp
        des290   = 0.91433e0_wp
        s2alp    = 0.66778e0_wp
        s2rho    = 5.30472e0_wp
        deh2     = 4.74806e0_wp
        s20alp   = 0.5e0_wp
        s20rho   = 6.0e0_wp
        c2a      = 1.38661e0_wp
        c2b      = 0.27362e0_wp
        c2c      = 0.03764e0_wp
        una2p    = 2.1037e0_wp
        g_ehfid  = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        cevau  = 1.e0_wp / 27.2113961e0_wp
        onethd = 1.0e0_wp / 3.0e0_wp
        twothd = 2.0e0_wp / 3.0e0_wp
        do i = 1, nt
            ! useful quantities calculated first
            do g_i_ = 1, g_p_
                g_rmid(g_i_) = 0.5e0_wp * g_r(g_i_,i,3) + 0.5e0_wp * g_r(g_i_,i,1)
            enddo
            rmid = 0.5e0_wp * (r(i,1) + r(i,3))
            rave = 0.5e0_wp * (r(i,1) - r(i,3))
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
            cosa = (r1s + r3s - r2s) / (2 * r(i,1) * r(i,3))
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            d7_v  = 2.0e0_wp * r(i,1)
            d9_v  = d7_v * r(i,2)
            d10_v = (r1s + r2s - r3s) / d9_v
            d2_b  = 1.0e0_wp / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0e0_wp
            do g_i_ = 1, g_p_
                g_cos1(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,1) + (-d2_b) * g_r3s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r1s(g_i_)
            enddo
            cos1  = d10_v
            d7_v  = 2.0e0_wp * r(i,3)
            d9_v  = d7_v * r(i,2)
            d10_v = (r3s + r2s - r1s) / d9_v
            d2_b  = 1.0e0_wp / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0e0_wp
            do g_i_ = 1, g_p_
                g_cos2(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,3) + (-d2_b) * g_r1s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r3s(g_i_)
            enddo
            cos2 = d10_v
            do g_i_ = 1, g_p_
                g_cosg(g_i_) = (-0.5e0_wp) * g_cos2(g_i_) + 0.5e0_wp * g_cos1(g_i_ )
            enddo
            cosg = 0.5e0_wp * (cos1 - cos2)
            ! The diatomic curves . . .
            ! The first NaH singlet
            denah  = 1.97134878e0_wp
            renah  = 3.566044e0_wp
            bf     = 0.864e0_wp
            b0     = 0.594e0_wp
            gamnah = 7.376e0_wp
            d3_v   = (r(i,1) / gamnah) ** ( 8 - 2)
            d3_v   = d3_v * (r(i,1) / gamnah)
            d2_p   = 8 *  d3_v
            d3_v   = d3_v * (r(i,1) / gamnah)
            d7_v   = (r(i,1) / gamnah) ** ( 8 - 2)
            d7_v   = d7_v * (r(i,1) / gamnah)
            d1_p   = 8 *  d7_v
            d7_v   = d7_v * (r(i,1) / gamnah)
            d8_v   = bf + d7_v
            d9_v   = bf * (b0 + d3_v) / d8_v
            d6_b   = (-d9_v) / d8_v * d1_p * (1.0e0_wp / gamnah) + 1.0e0_wp / d8_v * bf * d2_p * (1.0e0_wp / gamnah)
            do g_i_ = 1, g_p_
                g_beta(g_i_) = d6_b * g_r(g_i_,i,1)
            enddo
            beta = d9_v
            d4_v = r(i,1) - renah
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta) * g_r(g_i_,i,1) + (-d4_v) * g_beta( g_i_)
            enddo
            d1_w = (-beta) * d4_v
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
            s1 = d2_v * d3_v
            ! The first NaH triplet
            do g_i_ = 1, g_p_
                g_g10(g_i_) = g_r(g_i_,i,1)
            enddo
            g10 = r(i,1) - ret10
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet10) * g_g10(g_i_)
            enddo
            d1_w = (-bet10) * g10
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-bet10b) * g_g10(g_i_)
            enddo
            d2_w = (-bet10b) * g10
            d2_v = exp(d1_w)
            d2_p = d2_v
            d5_v = exp(d2_w)
            d1_p = d5_v
            d6_b = f10 * d2_p
            do g_i_ = 1, g_p_
                g_t10(g_i_) = d1_p * g_d2_w(g_i_) + d6_b * g_d1_w(g_i_)
            enddo
            t10  = f10 * d2_v + d5_v
            d2_b = 1.0e0_wp / (1.0e0_wp + f10)
            do g_i_ = 1, g_p_
                g_t10(g_i_) = d2_b * g_t10(g_i_)
            enddo
            t10 = t10 / (1.0e0_wp + f10)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet190) * g_r(g_i_,i,1)
            enddo
            d1_w = (-bet190) * (r(i,1) - ret190)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xet190(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xet190 = d2_v
            d2_v   = det190 * onethd * xet190
            d3_v   = xet190 + 2.0e0_wp
            d4_b   = d2_v + d3_v * (det190 * onethd)
            do g_i_ = 1, g_p_
                g_t190(g_i_) = d4_b * g_xet190(g_i_)
            enddo
            t190 = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet1c) * g_r(g_i_,i,1)
            enddo
            d1_w = (-bet1c) * (r(i,1) - ret1c)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xet1c(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xet1c = d2_v
            d2_v  = det1c * onethd * xet1c
            d3_v  = xet1c + 2.0e0_wp
            d4_b  = d2_v + d3_v * (det1c * onethd)
            do g_i_ = 1, g_p_
                g_t1c(g_i_) = d4_b * g_xet1c(g_i_)
            enddo
            t1c = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t1alp * g_r(g_i_,i,2)
            enddo
            d1_w = t1alp * (r(i,2) - t1rho)
            d4_v = (t1c - t190) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.0e0_wp + d6_v
            d7_b = d4_v * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_t1m(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_t1c(g_i_) + d2_b * g_t190(g_i_)
            enddo
            t1m  = t190 + d4_v * d7_v
            d3_v = t10 - t1m
            d5_v = cosg * cosg
            d1_p = 2.0e0_wp * cosg
            d6_b = d3_v * d1_p
            d2_b = 1.0e0_wp + (-d5_v)
            do g_i_ = 1, g_p_
                g_t1(g_i_) = d6_b * g_cosg(g_i_) + d5_v * g_t10(g_i_) + d2_b * g_t1m(g_i_)
            enddo
            t1 = t1m + d3_v * d5_v
            ! The H2 triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet20) * g_r(g_i_,i,2)
            enddo
            d1_w = (-bet20) * (r(i,2) - ret20)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xet20(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xet20 = d2_v
            d2_v  = det20 * onethd * xet20
            d3_v  = xet20 + 2.0e0_wp
            d4_b  = d2_v + d3_v * (det20 * onethd)
            do g_i_ = 1, g_p_
                g_t20(g_i_) = d4_b * g_xet20(g_i_)
            enddo
            t20 = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet290) * g_r(g_i_,i,2)
            enddo
            d1_w = (-bet290) * (r(i,2) - ret290)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xet290(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xet290 = d2_v
            d2_v   = det290 * onethd * xet290
            d3_v   = xet290 + 2.0e0_wp
            d4_b   = d2_v + d3_v * (det290 * onethd)
            do g_i_ = 1, g_p_
                g_t290(g_i_) = d4_b * g_xet290(g_i_)
            enddo
            t290 = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t2alp * g_rmid(g_i_)
            enddo
            d1_w  = t2alp * (rmid - t2rho)
            d3_v  = t290 - t20
            d5_v  = cosg * cosg
            d2_p  = 2.0e0_wp * cosg
            d6_v  = 1.0e0_wp - d5_v
            d8_v  = d3_v * d6_v * 0.5e0_wp
            d10_v = tanh (d1_w)
            d1_p  = 1.0e0_wp - ( d10_v *  d10_v)
            d11_v = 1.0e0_wp - d10_v
            d7_b  = (-d8_v) * d1_p
            d8_b  = d11_v * 0.5e0_wp
            d9_b  = d8_b * d6_v
            d12_b = (-(d8_b * d3_v)) * d2_p
            d2_b  = 1.0e0_wp + (-d9_b)
            do g_i_ = 1, g_p_
                g_t2(g_i_) = d7_b * g_d1_w(g_i_) + d12_b * g_cosg(g_i_) + d9_b * g_t290(g_i_) + d2_b * g_t20(g_i_)
            enddo
            t2 = t20 + d8_v * d11_v
            ! The H2 singlet
            ! new singlet, designed to minic lower adiabatic behavior
            reh2   = 1.401e0_wp
            bf     = 1.627e0_wp
            b0     = 1.194e0_wp
            gamh2a = 2.323e0_wp
            gamh2b = 3.126e0_wp
            d5_v   = (r(i,2) / gamh2b) * (r(i,2) / gamh2b)
            d2_p   = 2.0e0_wp * (r(i,2) / gamh2b)
            d11_v  = (r(i,2) / gamh2b) * (r(i,2) / gamh2b)
            d1_p   = 2.0e0_wp * (r(i,2) / gamh2b)
            d12_v  = bf - r(i,2) / gamh2a + d11_v
            d13_v  = bf * (b0 - r(i,2) / gamh2a + d5_v) / d12_v
            d3_b   = (-d13_v) / d12_v
            d9_b   = 1.0e0_wp / d12_v * bf
            d7_b  = d3_b * d1_p * (1.0e0_wp / gamh2b) + (-d3_b) * (1.0e0_wp / gamh2a) + d9_b * d2_p * (1.0e0_wp / &
                & gamh2b) + (- d9_b) * (1.0e0_wp / gamh2a)
            do g_i_ = 1, g_p_
                g_beh2(g_i_) = d7_b * g_r(g_i_,i,2)
            enddo
            beh2 = d13_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s2alp * g_rmid(g_i_)
            enddo
            d1_w = s2alp * (rmid - s2rho)
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_sswt(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            sswt = 0.5e0_wp * (1.0e0_wp - d2_v)
            d3_b = deh2 * (des290 - 1.0e0_wp)
            do g_i_ = 1, g_p_
                g_deh2m(g_i_) = d3_b * g_sswt(g_i_)
            enddo
            deh2m = deh2 + deh2 * (des290 - 1.0e0_wp) * sswt
            d4_v  = r(i,2) - reh2
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beh2) * g_r(g_i_,i,2) + (-d4_v) * g_beh2( g_i_)
            enddo
            d1_w = (-beh2) * d4_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xeh2(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xeh2 = d2_v
            d3_v = deh2m * xeh2
            d4_v = xeh2 - 2.0e0_wp
            d6_b = d4_v * xeh2
            d5_b = d3_v + d4_v * deh2m
            do g_i_ = 1, g_p_
                g_s2a(g_i_) = d5_b * g_xeh2(g_i_) + d6_b * g_deh2m(g_i_)
            enddo
            s2a  = d3_v * d4_v + una2p
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
                g_s2(g_i_) = d9_b * g_d1_w(g_i_) + d8_b * g_t2(g_i_) + d7_b * g_s2a(g_i_)
            enddo
            s2   = 0.5e0_wp * (s2a + t2 - d5_v * d6_v)
            d2_v = (-1.34e0_wp) * beh2
            d4_v = r(i,2) - reh2
            d5_b = d4_v * (-1.34e0_wp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d2_v * g_r(g_i_,i,2) + d5_b * g_beh2(g_i_)
            enddo
            d1_w = d2_v * d4_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xeh2(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xeh2 = d2_v
            d2_v = (deh2 - una2p) * xeh2
            d3_v = xeh2 - 2.0e0_wp
            d4_b = d2_v + d3_v * (deh2 - una2p)
            do g_i_ = 1, g_p_
                g_s20b(g_i_) = d4_b * g_xeh2(g_i_)
            enddo
            s20b = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s20alp * g_rmid(g_i_)
            enddo
            d1_w  = s20alp * (rmid - s20rho)
            d3_v  = s20b - s2
            d5_v  = cosg * cosg
            d2_p  = 2.0e0_wp * cosg
            d7_v  = d3_v * d5_v * 0.5e0_wp
            d9_v  = tanh (d1_w)
            d1_p  = 1.0e0_wp - ( d9_v *  d9_v)
            d10_v = 1.0e0_wp - d9_v
            d7_b  = (-d7_v) * d1_p
            d8_b  = d10_v * 0.5e0_wp
            d9_b  = d8_b * d5_v
            d11_b = d8_b * d3_v * d2_p
            d2_b  = 1.0e0_wp + (-d9_b)
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d7_b * g_d1_w(g_i_) + d11_b * g_cosg(g_i_) + d9_b * g_s20b(g_i_) + d2_b * g_s2(g_i_)
            enddo
            s2 = s2 + d7_v * d10_v
            ! The second NaH singlet
            denah  = 1.97134878e0_wp
            renah  = 3.566044e0_wp
            bf     = 0.864e0_wp
            b0     = 0.594e0_wp
            gamnah = 7.376e0_wp
            d3_v   = (r(i,3) / gamnah) ** ( 8 - 2)
            d3_v   = d3_v * (r(i,3) / gamnah)
            d2_p   = 8 *  d3_v
            d3_v   = d3_v * (r(i,3) / gamnah)
            d7_v   = (r(i,3) / gamnah) ** ( 8 - 2)
            d7_v   = d7_v * (r(i,3) / gamnah)
            d1_p   = 8 *  d7_v
            d7_v   = d7_v * (r(i,3) / gamnah)
            d8_v   = bf + d7_v
            d9_v   = bf * (b0 + d3_v) / d8_v
            d6_b   = (-d9_v) / d8_v * d1_p * (1.0e0_wp / gamnah) + 1.0e0_wp / d8_v * bf * d2_p * (1.0e0_wp / gamnah)
            do g_i_ = 1, g_p_
                g_beta(g_i_) = d6_b * g_r(g_i_,i,3)
            enddo
            beta = d9_v
            d4_v = r(i,3) - renah
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-beta) * g_r(g_i_,i,3) + (-d4_v) * g_beta( g_i_)
            enddo
            d1_w = (-beta) * d4_v
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
                g_s3(g_i_) = d4_b * g_xenah(g_i_)
            enddo
            s3 = d2_v * d3_v
            ! The second NaH triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet1) * g_r(g_i_,i,3)
            enddo
            d1_w = (-bet1) * (r(i,3) - ret1)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xet3(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xet3 = d2_v
            d2_v = det1 * onethd * xet3
            d3_v = xet3 + 2.0e0_wp
            d4_b = d2_v + d3_v * (det1 * onethd)
            do g_i_ = 1, g_p_
                g_t30(g_i_) = d4_b * g_xet3(g_i_)
            enddo
            t30 = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_g30(g_i_) = g_r(g_i_,i,3)
            enddo
            g30 = r(i,3) - ret10
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet10) * g_g30(g_i_)
            enddo
            d1_w = (-bet10) * g30
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = (-bet10b) * g_g30(g_i_)
            enddo
            d2_w = (-bet10b) * g30
            d2_v = exp(d1_w)
            d2_p = d2_v
            d5_v = exp(d2_w)
            d1_p = d5_v
            d6_b = f10 * d2_p
            do g_i_ = 1, g_p_
                g_t30(g_i_) = d1_p * g_d2_w(g_i_) + d6_b * g_d1_w(g_i_)
            enddo
            t30  = f10 * d2_v + d5_v
            d2_b = 1.0e0_wp / (1.0e0_wp + f10)
            do g_i_ = 1, g_p_
                g_t30(g_i_) = d2_b * g_t30(g_i_)
            enddo
            t30 = t30 / (1.0e0_wp + f10)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet190) * g_r(g_i_,i,3)
            enddo
            d1_w = (-bet190) * (r(i,3) - ret190)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xet390(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xet390 = d2_v
            d2_v   = det190 * onethd * xet390
            d3_v   = xet390 + 2.0e0_wp
            d4_b   = d2_v + d3_v * (det190 * onethd)
            do g_i_ = 1, g_p_
                g_t390(g_i_) = d4_b * g_xet390(g_i_)
            enddo
            t390 = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet1c) * g_r(g_i_,i,3)
            enddo
            d1_w = (-bet1c) * (r(i,3) - ret1c)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xet3c(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xet3c = d2_v
            d2_v  = det1c * onethd * xet3c
            d3_v  = xet3c + 2.0e0_wp
            d4_b  = d2_v + d3_v * (det1c * onethd)
            do g_i_ = 1, g_p_
                g_t3c(g_i_) = d4_b * g_xet3c(g_i_)
            enddo
            t3c = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t1alp * g_r(g_i_,i,2)
            enddo
            d1_w = t1alp * (r(i,2) - t1rho)
            d4_v = (t3c - t390) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.0e0_wp + d6_v
            d7_b = d4_v * d1_p
            d8_b = d7_v * 0.5e0_wp
            d2_b = 1.0e0_wp + (-d8_b)
            do g_i_ = 1, g_p_
                g_t3m(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_t3c(g_i_) + d2_b * g_t390(g_i_)
            enddo
            t3m  = t390 + d4_v * d7_v
            d3_v = t30 - t3m
            d5_v = cosg * cosg
            d1_p = 2.0e0_wp * cosg
            d6_b = d3_v * d1_p
            d2_b = 1.0e0_wp + (-d5_v)
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d6_b * g_cosg(g_i_) + d5_v * g_t30(g_i_) + d2_b * g_t3m(g_i_)
            enddo
            t3 = t3m + d3_v * d5_v
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
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_7sd_2000_dpem_g_pot22_legacy_state_d8f976ce

!> @brief Compatibility interface for g_pot22.
subroutine g_pot22( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_nah2_leps_7sd_2000_dpem_g_pot22_legacy_state_d8f976ce, only: dispatch
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
