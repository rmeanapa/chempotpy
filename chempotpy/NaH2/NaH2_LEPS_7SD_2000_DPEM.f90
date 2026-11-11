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
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NaH2_LEPS_7SD_2000 diabatic potential matrix for NaH2.
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
        double precision              :: d, g
        integer                       :: iatom, idir, nder
        double precision              :: nt, r(1,3), r2(3), v(1)
        double precision              :: p(nstates)
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
        call diagonalize(nstates, u, p, t)
        r2 = r(1,:)*legacy_bohr_to_angstrom
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
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 2
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: u(nstates,nstates)
    double precision, intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem

!> @brief Evaluate EvdRdX for the NaH2_LEPS_7SD_2000_DPEM potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    implicit none
    integer, parameter            :: first_atoms(3)  = [1, 2, 1]
    integer, parameter            :: second_atoms(3) = [2, 3, 3]
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!> @brief Evaluate diagonalize for the NaH2_LEPS_7SD_2000_DPEM potential.
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

!> @brief Evaluate pot for the NaH2_LEPS_7SD_2000_DPEM potential.
subroutine pot( r, e, de, nt, nsurf )
    implicit none
    integer, parameter :: g_p_  = 3
    integer, parameter :: ldg_r = 3
    integer, parameter :: ldg_e = 3
    double precision :: ddis(3), de(3,nt), de1(3,1), de2(3,1), de3(3,1), e(nt), g_e(ldg_e,nt), g_r(ldg_r,nt,3), r(nt,3)
    integer            :: i, j, k, nsurf, nt
    ! initialize seed matrix here
    do i=1, nt
        do j=1, 3
            do k=1, 3
                if( j .eq. k )then
                    g_r(j,i,k) = 1.0d0
                else
                    g_r(j,i,k) = 0.0d0
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
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
        implicit none
        optional :: g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt
        double precision :: b0, bect, bect2, beh2, bes1, bf, cevau, cos1, cos2, cosa, cosg, dect, deh2, deh2m, des1, des20, &
            & e(:), gamh2a, gamh2b, onethd, r(:,:), r1s, r2s, r3s, rac2, rave, rect, reh2, res1, rmid, s1, s2, s2alp, s2rho, s3, &
            & sswt, twothd, vmid, xect, xeh2, xes1, xes3
        integer :: i, nt
        ! NaH singlet
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        double precision   :: d2_p, d7_b, d12_b, d9_b, d13_v, d12_v, d11_v, d2_v, d1_p, d2_b
        double precision   :: d3_v, d4_v, d3_b, d4_b, d5_v, d6_v, d7_v, d1_w, &
            & d9_v, d10_v
        double precision :: d5_b, d6_b, g_r1s(g_pmax_), g_r(:,:,:), &
            & g_r2s(g_pmax_), g_r3s(g_pmax_), g_rmid(g_pmax_), g_rave(g_pmax_), &
            & g_cos1(g_pmax_), g_cos2(g_pmax_)
        double precision :: g_cosg(g_pmax_), g_cosa(g_pmax_), g_d1_w(g_pmax_), g_xes1(g_pmax_), g_s1(g_pmax_), g_beh2(g_pmax_), g_sswt(g_pmax_), g_deh2m(g_pmax_), g_xeh2(g_pmax_), g_s2(g_pmax_)
        double precision :: g_xes3(g_pmax_), g_s3(g_pmax_), g_e(:,:), &
            & g_xect(g_pmax_), g_vmid(g_pmax_)
        integer :: g_ehfid
        ! H2 singlet
        ! H-Na-H term
        !  misc constants
        b0      = 0.d0
        bect    = 0.d0
        bect2   = 0.d0
        beh2    = 0.d0
        bes1    = 0.d0
        bf      = 0.d0
        cevau   = 0.d0
        cos1    = 0.d0
        cos2    = 0.d0
        cosa    = 0.d0
        cosg    = 0.d0
        d10_v   = 0.d0
        d11_v   = 0.d0
        d12_b   = 0.d0
        d12_v   = 0.d0
        d13_v   = 0.d0
        d1_p    = 0.d0
        d1_w    = 0.d0
        d2_b    = 0.d0
        d2_p    = 0.d0
        d2_v    = 0.d0
        d3_b    = 0.d0
        d3_v    = 0.d0
        d4_b    = 0.d0
        d4_v    = 0.d0
        d5_b    = 0.d0
        d5_v    = 0.d0
        d6_b    = 0.d0
        d6_v    = 0.d0
        d7_b    = 0.d0
        d7_v    = 0.d0
        d9_b    = 0.d0
        d9_v    = 0.d0
        dect    = 0.d0
        deh2    = 0.d0
        deh2m   = 0.d0
        des1    = 0.d0
        des20   = 0.d0
        g_beh2  = 0.d0
        g_cos1  = 0.d0
        g_cos2  = 0.d0
        g_cosa  = 0.d0
        g_cosg  = 0.d0
        g_d1_w  = 0.d0
        g_deh2m = 0.d0
        g_r1s   = 0.d0
        g_r2s   = 0.d0
        g_r3s   = 0.d0
        g_rave  = 0.d0
        g_rmid  = 0.d0
        g_s1    = 0.d0
        g_s2    = 0.d0
        g_s3    = 0.d0
        g_sswt  = 0.d0
        g_vmid  = 0.d0
        g_xect  = 0.d0
        g_xeh2  = 0.d0
        g_xes1  = 0.d0
        g_xes3  = 0.d0
        gamh2a  = 0.d0
        gamh2b  = 0.d0
        onethd  = 0.d0
        r1s     = 0.d0
        r2s     = 0.d0
        r3s     = 0.d0
        rac2    = 0.d0
        rave    = 0.d0
        rect    = 0.d0
        reh2    = 0.d0
        res1    = 0.d0
        rmid    = 0.d0
        s1      = 0.d0
        s2      = 0.d0
        s2alp   = 0.d0
        s2rho   = 0.d0
        s3      = 0.d0
        sswt    = 0.d0
        twothd  = 0.d0
        vmid    = 0.d0
        xect    = 0.d0
        xeh2    = 0.d0
        xes1    = 0.d0
        xes3    = 0.d0
        des1    = 0.55118d0
        bes1    = 1.05118d0
        res1    = 3.38583d0
        reh2    = 1.401d0
        bf      = 1.627d0
        b0      = 1.194d0
        gamh2a  = 2.323d0
        gamh2b  = 3.126d0
        deh2    = 4.74806d0
        des20   = 0.996d0
        s2alp   = 3.0d0
        s2rho   = 7.0d0
        dect    = 2.24409d0
        bect    = 0.57874d0
        rect    = 3.46457d0
        bect2   = 0.87143d0
        g_ehfid = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        rac2   = 1.d0 / dsqrt(2.d0)
        cevau  = 1.d0 / 27.2113961d0
        onethd = 1.0d0 / 3.0d0
        twothd = 2.0d0 / 3.0d0
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
                g_rmid(g_i_) = 0.5d0 * g_r(g_i_,i,3) + 0.5d0 * g_r(g_i_,i,1)
            enddo
            rmid = 0.5d0 * (r(i,1) + r(i,3))
            do g_i_ = 1, g_p_
                g_rave(g_i_) = (-0.5d0) * g_r(g_i_,i,3) + 0.5d0 * g_r(g_i_,i,1)
            enddo
            rave = 0.5d0 * (r(i,1) - r(i,3))
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            d7_v  = 2.0d0 * r(i,1)
            d9_v  = d7_v * r(i,2)
            d10_v = (r1s + r2s - r3s) / d9_v
            d2_b  = 1.0d0 / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0d0
            do g_i_ = 1, g_p_
                g_cos1(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,1) + (-d2_b) * g_r3s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r1s(g_i_)
            enddo
            cos1  = d10_v
            d7_v  = 2.0d0 * r(i,3)
            d9_v  = d7_v * r(i,2)
            d10_v = (r3s + r2s - r1s) / d9_v
            d2_b  = 1.0d0 / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0d0
            do g_i_ = 1, g_p_
                g_cos2(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,3) + (-d2_b) * g_r1s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r3s(g_i_)
            enddo
            cos2 = d10_v
            do g_i_ = 1, g_p_
                g_cosg(g_i_) = (-0.5d0) * g_cos2(g_i_) + 0.5d0 * g_cos1(g_i_ )
            enddo
            cosg  = 0.5d0 * (cos1 - cos2)
            d7_v  = dble(2) * r(i,1)
            d9_v  = d7_v * r(i,3)
            d10_v = (r1s + r3s - r2s) / d9_v
            d2_b  = 1.0d0 / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,3) * dble(2)
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
            d3_v = xes1 + 2.0d0
            d4_b = d2_v + d3_v * (des1 * onethd)
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d4_b * g_xes1(g_i_)
            enddo
            s1 = d2_v * d3_v
            ! The H2 singlet
            reh2   = 1.401d0
            bf     = 1.627d0
            b0     = 1.194d0
            gamh2a = 2.323d0
            gamh2b = 3.126d0
            d5_v   = (r(i,2) / gamh2b) * (r(i,2) / gamh2b)
            d2_p   = 2.0d0 * (r(i,2) / gamh2b)
            d11_v  = (r(i,2) / gamh2b) * (r(i,2) / gamh2b)
            d1_p   = 2.0d0 * (r(i,2) / gamh2b)
            d12_v  = bf - r(i,2) / gamh2a + d11_v
            d13_v  = bf * (b0 - r(i,2) / gamh2a + d5_v) / d12_v
            d3_b   = (-d13_v) / d12_v
            d9_b   = 1.0d0 / d12_v * bf
            d7_b  = d3_b * d1_p * (1.0d0 / gamh2b) + (-d3_b) * (1.0d0 / gamh2a) + d9_b * d2_p * (1.0d0 / gamh2b) + (- &
                & d9_b) * (1.0d0 / gamh2a)
            do g_i_ = 1, g_p_
                g_beh2(g_i_) = d7_b * g_r(g_i_,i,2)
            enddo
            beh2 = d13_v
            d3_b = 1.0d0 / s2alp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_rmid(g_i_)
            enddo
            d1_w = (rmid - s2rho) / s2alp
            d2_v = tanh (d1_w)
            d2_p = 1.0d0 - ( d2_v *  d2_v)
            d4_v = 0.5d0 * (1.0d0 - d2_v)
            d6_v = cosg * cosg
            d1_p = 2.0d0 * cosg
            d4_b = d4_v * d1_p
            d7_b = (-(d6_v * 0.5d0)) * d2_p
            do g_i_ = 1, g_p_
                g_sswt(g_i_) = d4_b * g_cosg(g_i_) + d7_b * g_d1_w(g_i_)
            enddo
            sswt = d4_v * d6_v
            d3_b = deh2 * (des20 - 1.0d0)
            do g_i_ = 1, g_p_
                g_deh2m(g_i_) = d3_b * g_sswt(g_i_)
            enddo
            deh2m = deh2 + deh2 * (des20 - 1.0d0) * sswt
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
            d4_v = xeh2 - 2.0d0
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
            d3_v = xes3 + 2.0d0
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
            d1_p = 2.0d0 * rave
            d3_b = (-bect2) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_rave(g_i_)
            enddo
            d1_w  = (-bect2) * d2_v
            d2_v  = dect * xect
            d3_v  = xect - 2.0d0
            d4_v  = d2_v * d3_v
            d6_v  = exp(d1_w)
            d2_p  = d6_v
            d7_v  = d4_v * d6_v
            d11_v = (0.5d0 * (1.0d0 - cosa)) * (0.5d0 * (1.0d0 - cosa))
            d1_p  = 2.0d0 * (0.5d0 * (1.0d0 - cosa))
            d6_b  = -(d7_v * d1_p * 0.5d0)
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
    implicit none
    integer          :: g_p_
    double precision :: r(nt,3)
    double precision :: g_r(ldg_r,nt,3)
    integer          :: ldg_r
    double precision :: e(nt)
    double precision :: g_e(ldg_e,nt)
    integer          :: ldg_e
    integer          :: nt
    call dispatch(g_p_=g_p_, r=r, g_r=g_r, ldg_r=ldg_r, e=e, g_e=g_e, ldg_e=ldg_e, nt=nt)
end subroutine g_pot11

! U12 coupling surface . . .
!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_nah2_leps_7sd_2000_dpem_g_pot12_legacy_state_d8f976ce
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
        implicit none
        optional :: g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt
        double precision   :: ap, alp, bet, rnagsq, rnag, cosg, rd0, r20
        double precision   :: r1s, r2s, r3s, cos1, cos2, theta
        double precision   :: gc, gc2
        double precision   :: cevau
        integer            :: nt, i, j, k
        logical            :: leven
        double precision   :: r(:,:), e(:)
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        double precision   :: d2_p, d2_w, d1_w, d4_b, d5_b, d10_v, d9_v, d2_v, d8_b, d2_b
        double precision   :: d1_p, d3_v, d4_v, d5_v, d7_b, d7_v, d6_b, d3_b, &
            & g_r1s(g_pmax_), g_r(:,:,:)
        double precision :: g_r2s(g_pmax_), g_r3s(g_pmax_), g_rnagsq(g_pmax_), g_rnag(g_pmax_), g_cos1(g_pmax_), g_cos2(g_pmax_), g_cosg(g_pmax_), g_gc(g_pmax_), g_gc2(g_pmax_), g_d1_w(g_pmax_)
        double precision :: g_d2_w(g_pmax_), g_e(:,:)
        integer          :: g_ehfid
        alp      = 0.d0
        ap       = 0.d0
        bet      = 0.d0
        cevau    = 0.d0
        cos1     = 0.d0
        cos2     = 0.d0
        cosg     = 0.d0
        d10_v    = 0.d0
        d1_p     = 0.d0
        d1_w     = 0.d0
        d2_b     = 0.d0
        d2_p     = 0.d0
        d2_v     = 0.d0
        d2_w     = 0.d0
        d3_b     = 0.d0
        d3_v     = 0.d0
        d4_b     = 0.d0
        d4_v     = 0.d0
        d5_b     = 0.d0
        d5_v     = 0.d0
        d6_b     = 0.d0
        d7_b     = 0.d0
        d7_v     = 0.d0
        d8_b     = 0.d0
        d9_v     = 0.d0
        g_cos1   = 0.d0
        g_cos2   = 0.d0
        g_cosg   = 0.d0
        g_d1_w   = 0.d0
        g_d2_w   = 0.d0
        g_gc     = 0.d0
        g_gc2    = 0.d0
        g_r1s    = 0.d0
        g_r2s    = 0.d0
        g_r3s    = 0.d0
        g_rnag   = 0.d0
        g_rnagsq = 0.d0
        gc       = 0.d0
        gc2      = 0.d0
        r1s      = 0.d0
        r20      = 0.d0
        r2s      = 0.d0
        r3s      = 0.d0
        rd0      = 0.d0
        rnag     = 0.d0
        rnagsq   = 0.d0
        theta    = 0.d0
        ap       = 0.835d0
        alp      = 4.0d0
        rd0      = 4.0d0
        bet      = 0.5d0
        r20      = -1.3d0
        g_ehfid  = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        theta = 15.6d0 * dacos(-1.0d0) / 180.0d0
        cevau = 1.d0 / 27.2113961d0
        ! set leven true for an even U12 surface, and false for an odd
        ! surface
        leven = .true.
        do i = 1, nt
            d2_v = r(i,1) * r(i,1)
            d1_p = 2.0d0 * r(i,1)
            do g_i_ = 1, g_p_
                g_r1s(g_i_) = d1_p * g_r(g_i_,i,1)
            enddo
            r1s  = d2_v
            d2_v = r(i,2) * r(i,2)
            d1_p = 2.0d0 * r(i,2)
            do g_i_ = 1, g_p_
                g_r2s(g_i_) = d1_p * g_r(g_i_,i,2)
            enddo
            r2s  = d2_v
            d2_v = r(i,3) * r(i,3)
            d1_p = 2.0d0 * r(i,3)
            do g_i_ = 1, g_p_
                g_r3s(g_i_) = d1_p * g_r(g_i_,i,3)
            enddo
            r3s = d2_v
            ! here we do the function evaluations . . .
            do g_i_ = 1, g_p_
                g_rnagsq(g_i_) = (-0.25d0) * g_r2s(g_i_) + 0.5d0 * g_r3s(g_i_) + 0.5d0 * g_r1s(g_i_)
            enddo
            rnagsq = 0.5d0 * r1s + 0.5d0 * r3s - 0.25d0 * r2s
            if( rnagsq .gt. 0.0d0 )then
                d2_v = sqrt(rnagsq)
                if( rnagsq .gt. 0.0d0 )then
                    d1_p = 1.0d0 / (2.0d0 *  d2_v)
                else
                    d1_p = 0.0d0
                endif
                do g_i_ = 1, g_p_
                    g_rnag(g_i_) = d1_p * g_rnagsq(g_i_)
                enddo
                rnag = d2_v
            else
                do g_i_ = 1, g_p_
                    g_rnag(g_i_) = 0.0d0
                enddo
                rnag = 0.0d0
            endif
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            d7_v  = 2.0d0 * r(i,1)
            d9_v  = d7_v * r(i,2)
            d10_v = (r1s + r2s - r3s) / d9_v
            d2_b  = 1.0d0 / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0d0
            do g_i_ = 1, g_p_
                g_cos1(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,1) + (-d2_b) * g_r3s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r1s(g_i_)
            enddo
            cos1  = d10_v
            d7_v  = 2.0d0 * r(i,3)
            d9_v  = d7_v * r(i,2)
            d10_v = (r3s + r2s - r1s) / d9_v
            d2_b  = 1.0d0 / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0d0
            do g_i_ = 1, g_p_
                g_cos2(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,3) + (-d2_b) * g_r1s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r3s(g_i_)
            enddo
            cos2 = d10_v
            do g_i_ = 1, g_p_
                g_cosg(g_i_) = (-0.5d0) * g_cos2(g_i_) + 0.5d0 * g_cos1(g_i_ )
            enddo
            cosg = 0.5d0 * (cos1 - cos2)
            d4_b = dsin(theta)
            d5_b = dcos(theta)
            do g_i_ = 1, g_p_
                g_gc(g_i_) = d4_b * g_r(g_i_,i,2) + d5_b * g_rnag(g_i_)
            enddo
            gc   = dcos(theta) * rnag + dsin(theta) * r(i,2)
            d4_b = -dcos(theta)
            d5_b = dsin(theta)
            do g_i_ = 1, g_p_
                g_gc2(g_i_) = d4_b * g_r(g_i_,i,2) + d5_b * g_rnag(g_i_)
            enddo
            gc2  = dsin(theta) * rnag - dcos(theta) * r(i,2)
            d3_v = (gc - rd0) ** ( 4 - 2)
            d3_v = d3_v * (gc - rd0)
            d1_p = 4 *  d3_v
            d3_v = d3_v * (gc - rd0)
            d4_b = (-((1.0d0 / alp) ** 4)) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_gc(g_i_)
            enddo
            d1_w = (-((1.0d0 / alp) ** 4)) * d3_v
            d3_v = (gc2 - r20) * (gc2 - r20)
            d1_p = 2.0d0 * (gc2 - r20)
            d4_b = (-((1.0d0 / bet) ** 2)) * d1_p
            do g_i_ = 1, g_p_
                g_d2_w(g_i_) = d4_b * g_gc2(g_i_)
            enddo
            d2_w = (-((1.0d0 / bet) ** 2)) * d3_v
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
                d1_p = 2.0d0 * e(i)
                d4_v = d3_v + 1.0d-8
                d5_v = d2_v / d4_v
                d5_b = (-d5_v) / d4_v * d1_p + 1.0d0 / d4_v * d2_p
                do g_i_ = 1, g_p_
                    g_d1_w(g_i_) = d5_b * g_e(g_i_,i)
                enddo
                d1_w = d5_v
                d2_v = sqrt(d1_w)
                if( d1_w .gt. 0.0d0 )then
                    d1_p = 1.0d0 / (2.0d0 *  d2_v)
                else
                    d1_p = 0.0d0
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
    implicit none
    integer          :: g_p_
    double precision :: r(nt,3)
    double precision :: g_r(ldg_r,nt,3)
    integer          :: ldg_r
    double precision :: e(nt)
    double precision :: g_e(ldg_e,nt)
    integer          :: ldg_e
    integer          :: nt
    call dispatch(g_p_=g_p_, r=r, g_r=g_r, ldg_r=ldg_r, e=e, g_e=g_e, ldg_e=ldg_e, nt=nt)
end subroutine g_pot12

! the upper diabatic surface  U22
!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_nah2_leps_7sd_2000_dpem_g_pot22_legacy_state_d8f976ce
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
        implicit none
        optional :: g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt
        double precision :: b0, beh2, bet1, bet10, bet10b, bet190, bet1c, bet20, bet290, beta, bf, c2a, c2b, c2c, cevau, cos1, &
            & cos2, cosa, cosg, coul1, coul2, coul3, cplg2, deh2, deh2m, denah, des290, det1, det190, det1c, det20, det290, &
            & e(:), exch1, exch2, exch3, f10, g10, g30, gamh2a, gamh2b, gamnah, onethd, r(:,:), r1s, r2s, r3s, rave, reh2, &
            & renah, ret1, ret10, ret190, ret1c, ret20, ret290, rmid, s1, s2, s20alp, s20b, s20rho, s2a, s2alp, s2rho, s3, sswt, &
            & t1, t10, t190, t1alp, t1c, t1m, t1rho, t2, t20, t290, t2alp, t2rho, t3, t30, t390, t3c, t3m, twothd, una2p, w, &
            & xeh2, xenah, xet190, xet1c, xet20, xet290, xet3, xet390, xet3c
        integer :: i, nt
        ! NaH triplets
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        double precision   :: d3_p, d1_w, d12_b, d11_b, d13_v, d12_v, d11_v, &
            & d2_v, d3_v, d4_v
        double precision :: d2_w, d2_b, d3_b, d4_b, d5_v, d6_v, d7_v, d8_v, &
            & d9_v, d10_v
        double precision :: d5_b, d6_b, d7_b, d8_b, d9_b, d10_b, d1_p, d2_p, g_rmid(g_pmax_), g_r(:,:,:)
        double precision :: g_r1s(g_pmax_), g_r2s(g_pmax_), g_r3s(g_pmax_), &
            & g_cos1(g_pmax_), g_cos2(g_pmax_), g_cosg(g_pmax_), g_beta(g_pmax_ &
            & ), g_d1_w(g_pmax_), g_xenah(g_pmax_), g_s1(g_pmax_)
        double precision :: g_g10(g_pmax_), g_d2_w(g_pmax_), g_t10(g_pmax_), g_xet190(g_pmax_), g_t190(g_pmax_), g_xet1c(g_pmax_), g_t1c(g_pmax_), g_t1m(g_pmax_), g_t1(g_pmax_), g_xet20(g_pmax_)
        double precision :: g_t20(g_pmax_), g_xet290(g_pmax_), g_t290(g_pmax_), g_t2(g_pmax_), g_beh2(g_pmax_), g_sswt(g_pmax_), g_deh2m(g_pmax_), g_xeh2(g_pmax_), g_s2a(g_pmax_), g_s2(g_pmax_)
        double precision :: g_s20b(g_pmax_), g_s3(g_pmax_), g_xet3(g_pmax_), g_t30(g_pmax_), g_g30(g_pmax_), g_xet390(g_pmax_), g_t390(g_pmax_), g_xet3c(g_pmax_), g_t3c(g_pmax_), g_t3m(g_pmax_)
        double precision :: g_t3(g_pmax_), g_coul1(g_pmax_), g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_), g_exch3( &
            & g_pmax_), g_w(g_pmax_), g_cplg2(g_pmax_), g_e(:,:)
        integer :: g_ehfid
        ! H2 triplets
        ! H2 singlets
        ! misc constants
        b0       = 0.d0
        beh2     = 0.d0
        bet1     = 0.d0
        bet10    = 0.d0
        bet10b   = 0.d0
        bet190   = 0.d0
        bet1c    = 0.d0
        bet20    = 0.d0
        bet290   = 0.d0
        beta     = 0.d0
        bf       = 0.d0
        c2a      = 0.d0
        c2b      = 0.d0
        c2c      = 0.d0
        cevau    = 0.d0
        cos1     = 0.d0
        cos2     = 0.d0
        cosa     = 0.d0
        cosg     = 0.d0
        coul1    = 0.d0
        coul2    = 0.d0
        coul3    = 0.d0
        cplg2    = 0.d0
        d10_b    = 0.d0
        d10_v    = 0.d0
        d11_b    = 0.d0
        d11_v    = 0.d0
        d12_b    = 0.d0
        d12_v    = 0.d0
        d13_v    = 0.d0
        d1_p     = 0.d0
        d1_w     = 0.d0
        d2_b     = 0.d0
        d2_p     = 0.d0
        d2_v     = 0.d0
        d2_w     = 0.d0
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
        d9_b     = 0.d0
        d9_v     = 0.d0
        deh2     = 0.d0
        deh2m    = 0.d0
        denah    = 0.d0
        des290   = 0.d0
        det1     = 0.d0
        det190   = 0.d0
        det1c    = 0.d0
        det20    = 0.d0
        det290   = 0.d0
        exch1    = 0.d0
        exch2    = 0.d0
        exch3    = 0.d0
        f10      = 0.d0
        g10      = 0.d0
        g30      = 0.d0
        g_beh2   = 0.d0
        g_beta   = 0.d0
        g_cos1   = 0.d0
        g_cos2   = 0.d0
        g_cosg   = 0.d0
        g_coul1  = 0.d0
        g_coul2  = 0.d0
        g_coul3  = 0.d0
        g_cplg2  = 0.d0
        g_d1_w   = 0.d0
        g_d2_w   = 0.d0
        g_deh2m  = 0.d0
        g_exch1  = 0.d0
        g_exch2  = 0.d0
        g_exch3  = 0.d0
        g_g10    = 0.d0
        g_g30    = 0.d0
        g_r1s    = 0.d0
        g_r2s    = 0.d0
        g_r3s    = 0.d0
        g_rmid   = 0.d0
        g_s1     = 0.d0
        g_s2     = 0.d0
        g_s20b   = 0.d0
        g_s2a    = 0.d0
        g_s3     = 0.d0
        g_sswt   = 0.d0
        g_t1     = 0.d0
        g_t10    = 0.d0
        g_t190   = 0.d0
        g_t1c    = 0.d0
        g_t1m    = 0.d0
        g_t2     = 0.d0
        g_t20    = 0.d0
        g_t290   = 0.d0
        g_t3     = 0.d0
        g_t30    = 0.d0
        g_t390   = 0.d0
        g_t3c    = 0.d0
        g_t3m    = 0.d0
        g_w      = 0.d0
        g_xeh2   = 0.d0
        g_xenah  = 0.d0
        g_xet190 = 0.d0
        g_xet1c  = 0.d0
        g_xet20  = 0.d0
        g_xet290 = 0.d0
        g_xet3   = 0.d0
        g_xet390 = 0.d0
        g_xet3c  = 0.d0
        gamh2a   = 0.d0
        gamh2b   = 0.d0
        gamnah   = 0.d0
        onethd   = 0.d0
        r1s      = 0.d0
        r2s      = 0.d0
        r3s      = 0.d0
        rave     = 0.d0
        reh2     = 0.d0
        renah    = 0.d0
        ret1     = 0.d0
        ret10    = 0.d0
        ret190   = 0.d0
        ret1c    = 0.d0
        ret20    = 0.d0
        ret290   = 0.d0
        rmid     = 0.d0
        s1       = 0.d0
        s2       = 0.d0
        s20alp   = 0.d0
        s20b     = 0.d0
        s20rho   = 0.d0
        s2a      = 0.d0
        s2alp    = 0.d0
        s2rho    = 0.d0
        s3       = 0.d0
        sswt     = 0.d0
        t1       = 0.d0
        t10      = 0.d0
        t190     = 0.d0
        t1alp    = 0.d0
        t1c      = 0.d0
        t1m      = 0.d0
        t1rho    = 0.d0
        t2       = 0.d0
        t20      = 0.d0
        t290     = 0.d0
        t2alp    = 0.d0
        t2rho    = 0.d0
        t3       = 0.d0
        t30      = 0.d0
        t390     = 0.d0
        t3c      = 0.d0
        t3m      = 0.d0
        twothd   = 0.d0
        una2p    = 0.d0
        w        = 0.d0
        xeh2     = 0.d0
        xenah    = 0.d0
        xet190   = 0.d0
        xet1c    = 0.d0
        xet20    = 0.d0
        xet290   = 0.d0
        xet3     = 0.d0
        xet390   = 0.d0
        xet3c    = 0.d0
        ret10    = 3.22835d0
        bet10    = 0.57874d0
        bet10b   = 2.44094d0
        f10      = 9.44882d0
        det190   = 4.32441d0
        bet190   = 2.43701d0
        ret190   = 1.95480d0
        t1alp    = 0.30778d0
        t1rho    = 3.21850d0
        det1c    = 4.92126d0
        bet1c    = 0.65748d0
        ret1c    = 2.71654d0
        det20    = 3.29921d0
        bet20    = 2.64567d0
        ret20    = 2.44882d0
        det290   = 3.97638d0
        bet290   = 3.73937d0
        ret290   = 0.96024d0
        t2alp    = 1.0d0
        t2rho    = 9.0d0
        des290   = 0.91433d0
        s2alp    = 0.66778d0
        s2rho    = 5.30472d0
        deh2     = 4.74806d0
        s20alp   = 0.5d0
        s20rho   = 6.0d0
        c2a      = 1.38661d0
        c2b      = 0.27362d0
        c2c      = 0.03764d0
        una2p    = 2.1037d0
        g_ehfid  = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        cevau  = 1.d0 / 27.2113961d0
        onethd = 1.0d0 / 3.0d0
        twothd = 2.0d0 / 3.0d0
        do i = 1, nt
            ! useful quantities calculated first
            do g_i_ = 1, g_p_
                g_rmid(g_i_) = 0.5d0 * g_r(g_i_,i,3) + 0.5d0 * g_r(g_i_,i,1)
            enddo
            rmid = 0.5d0 * (r(i,1) + r(i,3))
            rave = 0.5d0 * (r(i,1) - r(i,3))
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
            d7_v  = 2.0d0 * r(i,1)
            d9_v  = d7_v * r(i,2)
            d10_v = (r1s + r2s - r3s) / d9_v
            d2_b  = 1.0d0 / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0d0
            do g_i_ = 1, g_p_
                g_cos1(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,1) + (-d2_b) * g_r3s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r1s(g_i_)
            enddo
            cos1  = d10_v
            d7_v  = 2.0d0 * r(i,3)
            d9_v  = d7_v * r(i,2)
            d10_v = (r3s + r2s - r1s) / d9_v
            d2_b  = 1.0d0 / d9_v
            d3_b  = (-d10_v) / d9_v
            d5_b  = d3_b * d7_v
            d6_b  = d3_b * r(i,2) * 2.0d0
            do g_i_ = 1, g_p_
                g_cos2(g_i_) = d5_b * g_r(g_i_,i,2) + d6_b * g_r(g_i_,i,3) + (-d2_b) * g_r1s(g_i_) + d2_b * &
                    & g_r2s(g_i_) + d2_b * g_r3s(g_i_)
            enddo
            cos2 = d10_v
            do g_i_ = 1, g_p_
                g_cosg(g_i_) = (-0.5d0) * g_cos2(g_i_) + 0.5d0 * g_cos1(g_i_ )
            enddo
            cosg = 0.5d0 * (cos1 - cos2)
            ! The diatomic curves . . .
            ! The first NaH singlet
            denah  = 1.97134878d0
            renah  = 3.566044d0
            bf     = 0.864d0
            b0     = 0.594d0
            gamnah = 7.376d0
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
            d6_b   = (-d9_v) / d8_v * d1_p * (1.0d0 / gamnah) + 1.0d0 / d8_v * bf * d2_p * (1.0d0 / gamnah)
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
            d3_v  = xenah - 2.0d0
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
            d2_b = 1.0d0 / (1.0d0 + f10)
            do g_i_ = 1, g_p_
                g_t10(g_i_) = d2_b * g_t10(g_i_)
            enddo
            t10 = t10 / (1.0d0 + f10)
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
            d3_v   = xet190 + 2.0d0
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
            d3_v  = xet1c + 2.0d0
            d4_b  = d2_v + d3_v * (det1c * onethd)
            do g_i_ = 1, g_p_
                g_t1c(g_i_) = d4_b * g_xet1c(g_i_)
            enddo
            t1c = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t1alp * g_r(g_i_,i,2)
            enddo
            d1_w = t1alp * (r(i,2) - t1rho)
            d4_v = (t1c - t190) * 0.5d0
            d6_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d6_v *  d6_v)
            d7_v = 1.0d0 + d6_v
            d7_b = d4_v * d1_p
            d8_b = d7_v * 0.5d0
            d2_b = 1.0d0 + (-d8_b)
            do g_i_ = 1, g_p_
                g_t1m(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_t1c(g_i_) + d2_b * g_t190(g_i_)
            enddo
            t1m  = t190 + d4_v * d7_v
            d3_v = t10 - t1m
            d5_v = cosg * cosg
            d1_p = 2.0d0 * cosg
            d6_b = d3_v * d1_p
            d2_b = 1.0d0 + (-d5_v)
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
            d3_v  = xet20 + 2.0d0
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
            d3_v   = xet290 + 2.0d0
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
            d2_p  = 2.0d0 * cosg
            d6_v  = 1.0d0 - d5_v
            d8_v  = d3_v * d6_v * 0.5d0
            d10_v = tanh (d1_w)
            d1_p  = 1.0d0 - ( d10_v *  d10_v)
            d11_v = 1.0d0 - d10_v
            d7_b  = (-d8_v) * d1_p
            d8_b  = d11_v * 0.5d0
            d9_b  = d8_b * d6_v
            d12_b = (-(d8_b * d3_v)) * d2_p
            d2_b  = 1.0d0 + (-d9_b)
            do g_i_ = 1, g_p_
                g_t2(g_i_) = d7_b * g_d1_w(g_i_) + d12_b * g_cosg(g_i_) + d9_b * g_t290(g_i_) + d2_b * g_t20(g_i_)
            enddo
            t2 = t20 + d8_v * d11_v
            ! The H2 singlet
            ! new singlet, designed to minic lower adiabatic behavior
            reh2   = 1.401d0
            bf     = 1.627d0
            b0     = 1.194d0
            gamh2a = 2.323d0
            gamh2b = 3.126d0
            d5_v   = (r(i,2) / gamh2b) * (r(i,2) / gamh2b)
            d2_p   = 2.0d0 * (r(i,2) / gamh2b)
            d11_v  = (r(i,2) / gamh2b) * (r(i,2) / gamh2b)
            d1_p   = 2.0d0 * (r(i,2) / gamh2b)
            d12_v  = bf - r(i,2) / gamh2a + d11_v
            d13_v  = bf * (b0 - r(i,2) / gamh2a + d5_v) / d12_v
            d3_b   = (-d13_v) / d12_v
            d9_b   = 1.0d0 / d12_v * bf
            d7_b  = d3_b * d1_p * (1.0d0 / gamh2b) + (-d3_b) * (1.0d0 / gamh2a) + d9_b * d2_p * (1.0d0 / gamh2b) + (- &
                & d9_b) * (1.0d0 / gamh2a)
            do g_i_ = 1, g_p_
                g_beh2(g_i_) = d7_b * g_r(g_i_,i,2)
            enddo
            beh2 = d13_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s2alp * g_rmid(g_i_)
            enddo
            d1_w = s2alp * (rmid - s2rho)
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_sswt(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            sswt = 0.5d0 * (1.0d0 - d2_v)
            d3_b = deh2 * (des290 - 1.0d0)
            do g_i_ = 1, g_p_
                g_deh2m(g_i_) = d3_b * g_sswt(g_i_)
            enddo
            deh2m = deh2 + deh2 * (des290 - 1.0d0) * sswt
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
            d4_v = xeh2 - 2.0d0
            d6_b = d4_v * xeh2
            d5_b = d3_v + d4_v * deh2m
            do g_i_ = 1, g_p_
                g_s2a(g_i_) = d5_b * g_xeh2(g_i_) + d6_b * g_deh2m(g_i_)
            enddo
            s2a  = d3_v * d4_v + una2p
            d2_b = 1.0d0 / 0.1d0
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-d2_b) * g_t2(g_i_) + d2_b * g_s2a(g_i_)
            enddo
            d1_w = (s2a - t2) / 0.1d0
            d5_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d5_v *  d5_v)
            d6_v = s2a - t2
            d6_b = (-0.5d0) * d5_v
            d9_b = (-0.5d0) * d6_v * d1_p
            d7_b = d6_b + 0.5d0
            d8_b = -d6_b + 0.5d0
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d9_b * g_d1_w(g_i_) + d8_b * g_t2(g_i_) + d7_b * g_s2a(g_i_)
            enddo
            s2   = 0.5d0 * (s2a + t2 - d5_v * d6_v)
            d2_v = (-1.34d0) * beh2
            d4_v = r(i,2) - reh2
            d5_b = d4_v * (-1.34d0)
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
            d3_v = xeh2 - 2.0d0
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
            d2_p  = 2.0d0 * cosg
            d7_v  = d3_v * d5_v * 0.5d0
            d9_v  = tanh (d1_w)
            d1_p  = 1.0d0 - ( d9_v *  d9_v)
            d10_v = 1.0d0 - d9_v
            d7_b  = (-d7_v) * d1_p
            d8_b  = d10_v * 0.5d0
            d9_b  = d8_b * d5_v
            d11_b = d8_b * d3_v * d2_p
            d2_b  = 1.0d0 + (-d9_b)
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d7_b * g_d1_w(g_i_) + d11_b * g_cosg(g_i_) + d9_b * g_s20b(g_i_) + d2_b * g_s2(g_i_)
            enddo
            s2 = s2 + d7_v * d10_v
            ! The second NaH singlet
            denah  = 1.97134878d0
            renah  = 3.566044d0
            bf     = 0.864d0
            b0     = 0.594d0
            gamnah = 7.376d0
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
            d6_b   = (-d9_v) / d8_v * d1_p * (1.0d0 / gamnah) + 1.0d0 / d8_v * bf * d2_p * (1.0d0 / gamnah)
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
            d3_v  = xenah - 2.0d0
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
            d3_v = xet3 + 2.0d0
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
            d2_b = 1.0d0 / (1.0d0 + f10)
            do g_i_ = 1, g_p_
                g_t30(g_i_) = d2_b * g_t30(g_i_)
            enddo
            t30 = t30 / (1.0d0 + f10)
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
            d3_v   = xet390 + 2.0d0
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
            d3_v  = xet3c + 2.0d0
            d4_b  = d2_v + d3_v * (det1c * onethd)
            do g_i_ = 1, g_p_
                g_t3c(g_i_) = d4_b * g_xet3c(g_i_)
            enddo
            t3c = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = t1alp * g_r(g_i_,i,2)
            enddo
            d1_w = t1alp * (r(i,2) - t1rho)
            d4_v = (t3c - t390) * 0.5d0
            d6_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d6_v *  d6_v)
            d7_v = 1.0d0 + d6_v
            d7_b = d4_v * d1_p
            d8_b = d7_v * 0.5d0
            d2_b = 1.0d0 + (-d8_b)
            do g_i_ = 1, g_p_
                g_t3m(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_t3c(g_i_) + d2_b * g_t390(g_i_)
            enddo
            t3m  = t390 + d4_v * d7_v
            d3_v = t30 - t3m
            d5_v = cosg * cosg
            d1_p = 2.0d0 * cosg
            d6_b = d3_v * d1_p
            d2_b = 1.0d0 + (-d5_v)
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d6_b * g_cosg(g_i_) + d5_v * g_t30(g_i_) + d2_b * g_t3m(g_i_)
            enddo
            t3 = t3m + d3_v * d5_v
            ! The final LEPS form
            do g_i_ = 1, g_p_
                g_coul1(g_i_) = 0.5d0 * g_t1(g_i_) + 0.5d0 * g_s1(g_i_)
            enddo
            coul1 = 0.5d0 * (s1 + t1)
            do g_i_ = 1, g_p_
                g_coul2(g_i_) = 0.5d0 * g_t2(g_i_) + 0.5d0 * g_s2(g_i_)
            enddo
            coul2 = 0.5d0 * (s2 + t2)
            do g_i_ = 1, g_p_
                g_coul3(g_i_) = 0.5d0 * g_t3(g_i_) + 0.5d0 * g_s3(g_i_)
            enddo
            coul3 = 0.5d0 * (s3 + t3)
            do g_i_ = 1, g_p_
                g_exch1(g_i_) = (-0.5d0) * g_t1(g_i_) + 0.5d0 * g_s1(g_i_)
            enddo
            exch1 = 0.5d0 * (s1 - t1)
            do g_i_ = 1, g_p_
                g_exch2(g_i_) = (-0.5d0) * g_t2(g_i_) + 0.5d0 * g_s2(g_i_)
            enddo
            exch2 = 0.5d0 * (s2 - t2)
            do g_i_ = 1, g_p_
                g_exch3(g_i_) = (-0.5d0) * g_t3(g_i_) + 0.5d0 * g_s3(g_i_)
            enddo
            exch3 = 0.5d0 * (s3 - t3)
            d4_v  = (exch1 - exch2) * (exch1 - exch2)
            d3_p  = 2.0d0 * (exch1 - exch2)
            d7_v  = (exch2 - exch3) * (exch2 - exch3)
            d2_p  = 2.0d0 * (exch2 - exch3)
            d10_v = (exch3 - exch1) * (exch3 - exch1)
            d1_p  = 2.0d0 * (exch3 - exch1)
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
            if( d1_w .gt. 0.0d0 )then
                d1_p = 1.0d0 / (2.0d0 *  d8_v)
            else
                d1_p = 0.0d0
            endif
            d5_b = (-(1.0d0 / dsqrt(2.0d0))) * d1_p
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d5_b * g_d1_w(g_i_) + g_coul3(g_i_) + g_coul2 (g_i_) + g_coul1(g_i_)
            enddo
            e(i) = coul1 + coul2 + coul3 + deh2 - d8_v / dsqrt(2.0d0)
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
    implicit none
    integer          :: g_p_
    double precision :: r(nt,3)
    double precision :: g_r(ldg_r,nt,3)
    integer          :: ldg_r
    double precision :: e(nt)
    double precision :: g_e(ldg_e,nt)
    integer          :: ldg_e
    integer          :: nt
    call dispatch(g_p_=g_p_, r=r, g_r=g_r, ldg_r=ldg_r, e=e, g_e=g_e, ldg_e=ldg_e, nt=nt)
end subroutine g_pot22
