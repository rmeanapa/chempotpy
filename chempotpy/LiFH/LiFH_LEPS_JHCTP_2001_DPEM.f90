!   System:                        LiFH
!   Functional form:               Modified Extended LEPS (2x2 diabatic fit)
!   Common name:                   LiFH surface fit H
!   Number of derivatives:         1
!   Number of electronic surfaces: 2
!   Interface:                     3-2V
!   References:      A. W. Jasper, M. D. Hack, A. Chakraborty, D. G. Tru
!                    and P. Piecuch, J. Chem. Phys. 115, 7945 (2001).
!   Notes:        This fit is qualitatively correct.
!                 It is quantitatively incorrect near the saddle point.
!   Protocol: POT evaluates the stateless diabatic surfaces.
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      ground and first excited electronic states, diabatic representati
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Li
!      infinitely far from the HF diatomic and R(HF) set equal to the
!      HF equilibrium diatomic value.
!   Parameters:
!   Coordinates:
!      Internal, Definition: R(1) = R(Li-H)
!                            R(2) = R(H-F)
!                            R(3) = R(Li-F)
!  all surfaces.  The potential is called by calling POT with the parame
!  r, e, nt, nsurf.  r(nt,3) and ei(nt) are double precision arrays of
!  dimension nt which is an integer representing the number of geometrie
!  input.  The integer nsurf selects the desired surface.  nsurf=1 is th
!  lower diabatic surface, nsurf=2 is the coupling surface, nsurf=3 is
!  the upper diabatic surface, nsurf = 4 selects the lower adiabatic
!  surface, and nsurf = 5 selects the upper adiabatic surface

!> @brief Implement LiFH_LEPS_JHCTP_2001, a diabatic potential matrix for LiFH.
module chempotpy_lifh_leps_jhctp_2001_dpem_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the LiFH_LEPS_JHCTP_2001 diabatic potential matrix for LiFH.
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
end module chempotpy_lifh_leps_jhctp_2001_dpem_pes

!> @brief Expose LiFH_LEPS_JHCTP_2001 through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_lifh_leps_jhctp_2001_dpem_pes, only: evaluate
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

!> @brief Evaluate EvdRdX for the LiFH_LEPS_JHCTP_2001_DPEM potential.
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

!> @brief Evaluate diagonalize for the LiFH_LEPS_JHCTP_2001_DPEM potential.
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

!> @brief Evaluate pot for the LiFH_LEPS_JHCTP_2001_DPEM potential.
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
module chempotpy_lifh_leps_jhctp_2001_dpem_g_pot11_entry_sta_99e0597c
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
        integer :: i, nt
        ! misc variables
        real(wp) :: r1s, r2s, r3s
        ! LiH variables
        real(wp) :: s1, relih, bes1, s1c2, s1c3, s1a, xelih, bes1b, &
            & s1c
        real(wp) :: rcut2del, rcut2rho, rcut2, r1r3, r1r3al
        ! HF variables
        real(wp) :: s2, behf, yhf, rehf, dehf, xehf
        ! LiF variables
        real(wp) :: s3, res3, bes3, des3, res3c, bes3c, s3c, s3a, xes3, des3c
        real(wp) :: rcutdel, rcutrho, rcut, r1r3al2
        real(wp) :: cevau
        real(wp) :: r(:,:), e(:)
        !! values of variables:
        ! LiH singlet:
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d1_p, d1_w, d4_v, d5_b, d4_b, d3_b, d2_b, d2_v, &
            & d3_v, g_r1r3(g_pmax_)
        real(wp) :: g_r(:,:,:), g_d1_w(g_pmax_), g_rcut2(g_pmax_), g_rcut(g_pmax_), g_xelih(g_pmax_), g_s1a(g_pmax_), g_s1c(g_pmax_), g_s1(g_pmax_), g_yhf(g_pmax_), g_behf(g_pmax_)
        real(wp) :: g_xehf(g_pmax_), g_s2(g_pmax_), g_xes3(g_pmax_), g_s3a(g_pmax_), g_s3c(g_pmax_), g_s3(g_pmax_), g_e(:,:)
        integer  :: g_ehfid
        ! HF singlet
        ! LiF singlet
        behf     = 0.e0_wp
        bes1     = 0.e0_wp
        bes1b    = 0.e0_wp
        bes3     = 0.e0_wp
        bes3c    = 0.e0_wp
        cevau    = 0.e0_wp
        d1_p     = 0.e0_wp
        d1_w     = 0.e0_wp
        d2_b     = 0.e0_wp
        d2_v     = 0.e0_wp
        d3_b     = 0.e0_wp
        d3_v     = 0.e0_wp
        d4_b     = 0.e0_wp
        d4_v     = 0.e0_wp
        d5_b     = 0.e0_wp
        dehf     = 0.e0_wp
        des3     = 0.e0_wp
        des3c    = 0.e0_wp
        g_behf   = 0.e0_wp
        g_d1_w   = 0.e0_wp
        g_r1r3   = 0.e0_wp
        g_rcut   = 0.e0_wp
        g_rcut2  = 0.e0_wp
        g_s1     = 0.e0_wp
        g_s1a    = 0.e0_wp
        g_s1c    = 0.e0_wp
        g_s2     = 0.e0_wp
        g_s3     = 0.e0_wp
        g_s3a    = 0.e0_wp
        g_s3c    = 0.e0_wp
        g_xehf   = 0.e0_wp
        g_xelih  = 0.e0_wp
        g_xes3   = 0.e0_wp
        g_yhf    = 0.e0_wp
        r1r3     = 0.e0_wp
        r1r3al   = 0.e0_wp
        r1r3al2  = 0.e0_wp
        r1s      = 0.e0_wp
        r2s      = 0.e0_wp
        r3s      = 0.e0_wp
        rcut     = 0.e0_wp
        rcut2    = 0.e0_wp
        rcut2del = 0.e0_wp
        rcut2rho = 0.e0_wp
        rcutdel  = 0.e0_wp
        rcutrho  = 0.e0_wp
        rehf     = 0.e0_wp
        relih    = 0.e0_wp
        res3     = 0.e0_wp
        res3c    = 0.e0_wp
        s1       = 0.e0_wp
        s1a      = 0.e0_wp
        s1c      = 0.e0_wp
        s1c2     = 0.e0_wp
        s1c3     = 0.e0_wp
        s2       = 0.e0_wp
        s3       = 0.e0_wp
        s3a      = 0.e0_wp
        s3c      = 0.e0_wp
        xehf     = 0.e0_wp
        xelih    = 0.e0_wp
        xes3     = 0.e0_wp
        yhf      = 0.e0_wp
        relih    = 1.53333e0_wp
        bes1     = 1.73333e0_wp
        s1c2     = 12.90323e0_wp
        s1c3     = 7.09677e0_wp
        bes1b    = 1.36e0_wp
        rcut2rho = 1.07143e0_wp
        rcut2del = 0.6e0_wp
        r1r3al   = 0.4e0_wp
        rehf     = 1.733e0_wp
        dehf     = 6.122e0_wp
        res3     = 3.66667e0_wp
        bes3     = 1.22667e0_wp
        des3     = 0.14286e0_wp
        res3c    = 2.48571e0_wp
        bes3c    = 1.92857e0_wp
        des3c    = 0.25e0_wp
        rcutrho  = 0.4e0_wp
        rcutdel  = 0.89333e0_wp
        r1r3al2  = 0.7333e0_wp
        cevau    = 0.036749309e0_wp
        g_ehfid  = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        do i = 1, nt
            r1s = r(i,1) ** 2
            r2s = r(i,2) ** 2
            r3s = r(i,3) ** 2
            do g_i_ = 1, g_p_
                g_r1r3(g_i_) = r1r3al * g_r(g_i_,i,2) + (-g_r(g_i_,i,3)) + g_r(g_i_,i,1)
            enddo
            r1r3 = r(i,1) - r(i,3) + r1r3al * r(i,2)
            d3_b = 1.0e0_wp / rcut2del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r1r3(g_i_)
            enddo
            d1_w = (r1r3 - rcut2rho) / rcut2del
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_rcut2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            rcut2 = 0.5e0_wp * (1.0e0_wp - d2_v)
            do g_i_ = 1, g_p_
                g_r1r3(g_i_) = r1r3al2 * g_r(g_i_,i,2) + (-g_r(g_i_,i,3) ) + g_r(g_i_,i,1)
            enddo
            r1r3 = r(i,1) - r(i,3) + r1r3al2 * r(i,2)
            d3_v = (r1r3 - rcutrho) * (r1r3 - rcutrho)
            d1_p = 2.0e0_wp * (r1r3 - rcutrho)
            d5_b = (-(1.0e0_wp / rcutdel ** 2)) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_r1r3(g_i_)
            enddo
            d1_w = (-d3_v) / rcutdel ** 2
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_rcut(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            rcut = d2_v
            !==============================================================k
            !This section contains the diatomic curves
            ! The  MH singlet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bes1) * g_r(g_i_,i,1)
            enddo
            d1_w = (-bes1) * (r(i,1) - relih)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xelih(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xelih = d2_v
            do g_i_ = 1, g_p_
                g_s1a(g_i_) = s1c2 * g_xelih(g_i_)
            enddo
            s1a = s1c2 * xelih
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bes1b) * g_r(g_i_,i,1)
            enddo
            d1_w = (-bes1b) * (r(i,1) - relih)
            d2_v = exp(d1_w)
            d1_p = d2_v
            d3_b = s1c3 * d1_p
            do g_i_ = 1, g_p_
                g_s1c(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            s1c  = s1c3 * d2_v
            d3_v = s1c - s1a
            d2_b = 1.0e0_wp + (-rcut2)
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d3_v * g_rcut2(g_i_) + rcut2 * g_s1c(g_i_) + d2_b * g_s1a(g_i_)
            enddo
            s1 = s1a + d3_v * rcut2
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
            d2_v = des3 * xes3
            d3_v = xes3 - 2.0e0_wp
            d4_b = d2_v + d3_v * des3
            do g_i_ = 1, g_p_
                g_s3a(g_i_) = d4_b * g_xes3(g_i_)
            enddo
            s3a = d2_v * d3_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bes3c) * g_r(g_i_,i,3)
            enddo
            d1_w = (-bes3c) * (r(i,3) - res3c)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xes3(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xes3 = d2_v
            d2_v = des3c * xes3
            d3_v = xes3 - 2.0e0_wp
            d4_b = d2_v + d3_v * des3c
            do g_i_ = 1, g_p_
                g_s3c(g_i_) = d4_b * g_xes3(g_i_)
            enddo
            s3c  = d2_v * d3_v
            d3_v = s3c - s3a
            d2_b = 1.0e0_wp + (-rcut)
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d3_v * g_rcut(g_i_) + rcut * g_s3c(g_i_) + d2_b * g_s3a(g_i_)
            enddo
            s3 = s3a + d3_v * rcut
            ! The final form
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
end module chempotpy_lifh_leps_jhctp_2001_dpem_g_pot11_entry_sta_99e0597c

!> @brief Compatibility interface for g_pot11.
subroutine g_pot11( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_lifh_leps_jhctp_2001_dpem_g_pot11_entry_sta_99e0597c, only: dispatch
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

! U12 coupling surface
!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_lifh_leps_jhctp_2001_dpem_g_pot12_entry_sta_99e0597c
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
        integer  :: i, nt
        real(wp) :: r(:,:), e(:)
        ! LiH variables
        integer  :: ng1
        real(wp) :: g1, g1r, g1a, h1, h1del, h1rho
        ! HF variables
        real(wp) :: j2p, j2, j2pdel, j2rho, j2prho, j2del
        ! LiF variables
        integer  :: ng3
        real(wp) :: g3, g3r, g3a, h3, h3del, h3rho
        ! cutoff variables
        real(wp) :: gcrho, theta, cevau, gcdel
        real(wp) :: cut2, gc
        ! LiH data
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d13_b, d12_b, d5_b, d9_v, d7_b, d8_b, d1_w, d1_p, d2_v, d3_v
        real(wp)           :: d4_v, d2_p, d2_b, d3_b, d4_b, d5_v, d6_v, d11_v, g_d1_w(g_pmax_), g_r(:,:,:)
        real(wp)           :: g_g1(g_pmax_), g_g3(g_pmax_), g_j2(g_pmax_), g_j2p(g_pmax_), g_h1(g_pmax_), g_h3(g_pmax_), g_e(:,:), g_gc(g_pmax_), g_cut2(g_pmax_)
        integer            :: g_ehfid
        ! LiF data
        ! HF data
        ! cutoff data
        cevau   = 0.e0_wp
        cut2    = 0.e0_wp
        d11_v   = 0.e0_wp
        d12_b   = 0.e0_wp
        d13_b   = 0.e0_wp
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
        d6_v    = 0.e0_wp
        d7_b    = 0.e0_wp
        d8_b    = 0.e0_wp
        d9_v    = 0.e0_wp
        g1      = 0.e0_wp
        g1a     = 0.e0_wp
        g1r     = 0.e0_wp
        g3      = 0.e0_wp
        g3a     = 0.e0_wp
        g3r     = 0.e0_wp
        g_cut2  = 0.e0_wp
        g_d1_w  = 0.e0_wp
        g_g1    = 0.e0_wp
        g_g3    = 0.e0_wp
        g_gc    = 0.e0_wp
        g_h1    = 0.e0_wp
        g_h3    = 0.e0_wp
        g_j2    = 0.e0_wp
        g_j2p   = 0.e0_wp
        gc      = 0.e0_wp
        gcdel   = 0.e0_wp
        gcrho   = 0.e0_wp
        h1      = 0.e0_wp
        h1del   = 0.e0_wp
        h1rho   = 0.e0_wp
        h3      = 0.e0_wp
        h3del   = 0.e0_wp
        h3rho   = 0.e0_wp
        j2      = 0.e0_wp
        j2del   = 0.e0_wp
        j2p     = 0.e0_wp
        j2pdel  = 0.e0_wp
        j2prho  = 0.e0_wp
        j2rho   = 0.e0_wp
        theta   = 0.e0_wp
        g1a     = 1.27742e0_wp
        g1r     = 2.5873e0_wp
        ng1     = 6
        h1rho   = 4.27097e0_wp
        h1del   = 2.5e0_wp
        g3a     = 0.48e0_wp
        g3r     = 3.47619e0_wp
        ng3     = 8
        h3rho   = 2.17742e0_wp
        h3del   = 0.56452e0_wp
        j2rho   = 1.15484e0_wp
        j2del   = 1.75806e0_wp
        j2prho  = 1.45161e0_wp
        j2pdel  = 0.98387e0_wp
        theta   = 0.0e0_wp
        gcrho   = 3.87097e0_wp
        gcdel   = 0.45806e0_wp
        cevau   = 0.036749309e0_wp
        g_ehfid = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
        do i = 1, nt
            d4_b = (-real(ng1, wp)) * (1.0e0_wp / g1r)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_r(g_i_,i,1)
            enddo
            d1_w = (-real(ng1, wp)) * (r(i,1) / g1r - 1.0e0_wp)
            if( ng1 .ne. 0 )then
                d3_v = (r(i,1) / g1r) ** ( ng1-1)
                d2_p = ng1 *  d3_v
                d3_v = d3_v * (r(i,1) / g1r)
            else
                !            Maybe this should be  d3_v = (r(i, 1) / g1r) **  ng1
                d3_v = 1.0e0_wp
                d2_p = 0.0e0_wp
            endif
            d4_v = g1a * d3_v
            d6_v = exp(d1_w)
            d1_p = d6_v
            d4_b = d4_v * d1_p
            d7_b = d6_v * g1a * d2_p * (1.0e0_wp / g1r)
            do g_i_ = 1, g_p_
                g_g1(g_i_) = d4_b * g_d1_w(g_i_) + d7_b * g_r(g_i_,i,1)
            enddo
            g1   = d4_v * d6_v
            d4_b = (-real(ng3, wp)) * (1.0e0_wp / g3r)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_r(g_i_,i,3)
            enddo
            d1_w = (-real(ng3, wp)) * (r(i,3) / g3r - 1.0e0_wp)
            if( ng3 .ne. 0 )then
                d3_v = (r(i,3) / g3r) ** ( ng3-1)
                d2_p = ng3 *  d3_v
                d3_v = d3_v * (r(i,3) / g3r)
            else
                !            Maybe this should be  d3_v = (r(i, 3) / g3r) **  ng3
                d3_v = 1.0e0_wp
                d2_p = 0.0e0_wp
            endif
            d4_v = g3a * d3_v
            d6_v = exp(d1_w)
            d1_p = d6_v
            d4_b = d4_v * d1_p
            d7_b = d6_v * g3a * d2_p * (1.0e0_wp / g3r)
            do g_i_ = 1, g_p_
                g_g3(g_i_) = d4_b * g_d1_w(g_i_) + d7_b * g_r(g_i_,i,3)
            enddo
            g3   = d4_v * d6_v
            d3_b = 1.0e0_wp / j2del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,2)
            enddo
            d1_w = (r(i,2) - j2rho) / j2del
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = 0.5e0_wp * d1_p
            do g_i_ = 1, g_p_
                g_j2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            j2   = 0.5e0_wp * (1.0e0_wp + d2_v)
            d3_b = 1.0e0_wp / j2pdel
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,2)
            enddo
            d1_w = (r(i,2) - j2prho) / j2pdel
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = 0.5e0_wp * d1_p
            do g_i_ = 1, g_p_
                g_j2p(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            j2p  = 0.5e0_wp * (1.0e0_wp + d2_v)
            d3_b = 1.0e0_wp / h1del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,1)
            enddo
            d1_w = (r(i,1) - h1rho) / h1del
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_h1(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            h1   = 0.5e0_wp * (1.0e0_wp - d2_v)
            d3_b = 1.0e0_wp / h3del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,3)
            enddo
            d1_w = (r(i,3) - h3rho) / h3del
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_h3(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            h3    = 0.5e0_wp * (1.0e0_wp - d2_v)
            d3_v  = g1 * j2
            d5_v  = 1.0e0_wp - h3
            d9_v  = g3 * j2p
            d11_v = 1.0e0_wp - h1
            d7_b  = d11_v * j2p
            d8_b  = d11_v * g3
            d12_b = d5_v * j2
            d13_b = d5_v * g1
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = (-d9_v) * g_h1(g_i_) + d8_b * g_j2p(g_i_) + d7_b * g_g3(g_i_) + (-d3_v) * g_h3(g_i_) + &
                    & d13_b * g_j2(g_i_) + d12_b * g_g1(g_i_)
            enddo
            e(i) = d3_v * d5_v + d9_v * d11_v
            d5_b = -sin(theta)
            d7_b = cos(theta)
            do g_i_ = 1, g_p_
                g_gc(g_i_) = d5_b * g_r(g_i_,i,3) + d7_b * g_r(g_i_,i,2)
            enddo
            gc   = (r(i,2) - gcrho) * cos(theta) - (r(i,3) - 2.9553e0_wp) * sin(theta)
            d2_b = 1.0e0_wp / gcdel
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d2_b * g_gc(g_i_)
            enddo
            d1_w = gc / gcdel
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_cut2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            cut2 = 0.5e0_wp - 0.5e0_wp * d2_v
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = e(i) * g_cut2(g_i_) + cut2 * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cut2
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
        enddo
    end subroutine dispatch
end module chempotpy_lifh_leps_jhctp_2001_dpem_g_pot12_entry_sta_99e0597c

!> @brief Compatibility interface for g_pot12.
subroutine g_pot12( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_lifh_leps_jhctp_2001_dpem_g_pot12_entry_sta_99e0597c, only: dispatch
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

! the upper diabatic surface  U22
!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_lifh_leps_jhctp_2001_dpem_g_pot22_entry_sta_99e0597c
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
        ! misc variables
        integer  :: i, j, k, nt
        real(wp) :: r(:,:), e(:), r1s, r2s, r3s, cs, cevau, uli2p
        real(wp) :: c2a, c2b, c2c, w, cplg2
        ! LiH variables
        real(wp) :: r1r3, xelih, s1c, bes1, rcut2del, rcut2rho, s1c3
        real(wp) :: s1c1, belihc, bet1, relih, s1c2, t1c1, t1c2
        real(wp) :: r1r3al, rcut2
        real(wp) :: coul1, exch1, s1, t1
        ! HF variables
        real(wp) :: behf, xehfa, dehfm, yhf, xehfa_180, t2a3, t2a2, &
            & bet2
        real(wp) :: t2a1, t2a4, dehf, rehf
        real(wp) :: s2alp, rehft, bet2_180, s2rho, dehfc
        real(wp) :: coul2, exch2, t2, s2
        ! LiF variables
        integer  :: nb
        real(wp) :: bet3, bet3_180, t3a2, t3a3, t3a4, belihc1
        real(wp) :: belifff, gam, belifc2, delifc1, delifc2, belifm
        real(wp) :: s3c, xelihf_180, xelif, s3a, delifm, gc2, cut2
        real(wp) :: del2_0, del2_180, rho2_0, rho2_180, theta2_0, grho
        real(wp) :: theta2_180, beliff, belifmod, belif, xelif_180, &
            & relif
        real(wp) :: belifc1, t3a1, gswt, delif, rho2, del2, theta2, &
            & gdel
        real(wp) :: coul3, exch3, s3, t3
        ! LiH singlets and triplets
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, g_p_, ldg_r, ldg_e
        real(wp)           :: d3_p, d10_v, d17_b, d11_v, d15_v, d14_b, d10_b, &
            & d2_v, d16_v, d2_b
        real(wp) :: d1_p, d3_v, d4_v, d5_v, d3_b, d4_b, d5_b, d6_v, &
            & d6_b, d2_p
        real(wp) :: d1_w, d7_v, d8_v, d9_v, d7_b, d8_b, d9_b, g_r1s &
            & (g_pmax_), g_r(:,:,:), g_r2s(g_pmax_)
        real(wp) :: g_r3s(g_pmax_), g_cs(g_pmax_), g_d1_w(g_pmax_), &
            & g_r1r3(g_pmax_), g_rcut2(g_pmax_), g_xelih(g_pmax_), g_s1(g_pmax_ &
            & ), g_t1(g_pmax_), g_s1c(g_pmax_), g_yhf(g_pmax_)
        real(wp) :: g_behf(g_pmax_), g_xehfa(g_pmax_), g_dehfm(g_pmax_), g_s2(g_pmax_), g_xehfa_180(g_pmax_), g_t2(g_pmax_), g_belif( &
            & g_pmax_), g_gswt(g_pmax_), g_beliff(g_pmax_), g_belifmod(g_pmax_)
        real(wp) :: g_xelif(g_pmax_), g_s3a(g_pmax_), g_delifm(g_pmax_), g_belifm(g_pmax_), g_s3c(g_pmax_), g_xelif_180(g_pmax_), g_t3(g_pmax_), g_rho2(g_pmax_), g_theta2(g_pmax_), g_del2(g_pmax_)
        real(wp) :: g_gc2(g_pmax_), g_cut2(g_pmax_), g_s3(g_pmax_), &
            & g_coul1(g_pmax_), g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_), g_exch3(g_pmax_), g_w(g_pmax_)
        real(wp) :: g_cplg2(g_pmax_), g_e(:,:)
        integer  :: g_ehfid
        ! HF singlet & triplet
        ! LiF singlets and triplets
        ! misc constants
        behf        = 0.e0_wp
        belif       = 0.e0_wp
        belifc1     = 0.e0_wp
        belifc2     = 0.e0_wp
        beliff      = 0.e0_wp
        belifff     = 0.e0_wp
        belifm      = 0.e0_wp
        belifmod    = 0.e0_wp
        belihc      = 0.e0_wp
        belihc1     = 0.e0_wp
        bes1        = 0.e0_wp
        bet1        = 0.e0_wp
        bet2        = 0.e0_wp
        bet2_180    = 0.e0_wp
        bet3        = 0.e0_wp
        bet3_180    = 0.e0_wp
        c2a         = 0.e0_wp
        c2b         = 0.e0_wp
        c2c         = 0.e0_wp
        cevau       = 0.e0_wp
        coul1       = 0.e0_wp
        coul2       = 0.e0_wp
        coul3       = 0.e0_wp
        cplg2       = 0.e0_wp
        cs          = 0.e0_wp
        cut2        = 0.e0_wp
        d10_b       = 0.e0_wp
        d10_v       = 0.e0_wp
        d11_v       = 0.e0_wp
        d14_b       = 0.e0_wp
        d15_v       = 0.e0_wp
        d16_v       = 0.e0_wp
        d17_b       = 0.e0_wp
        d1_p        = 0.e0_wp
        d1_w        = 0.e0_wp
        d2_b        = 0.e0_wp
        d2_p        = 0.e0_wp
        d2_v        = 0.e0_wp
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
        dehf        = 0.e0_wp
        dehfc       = 0.e0_wp
        dehfm       = 0.e0_wp
        del2        = 0.e0_wp
        del2_0      = 0.e0_wp
        del2_180    = 0.e0_wp
        delif       = 0.e0_wp
        delifc1     = 0.e0_wp
        delifc2     = 0.e0_wp
        delifm      = 0.e0_wp
        exch1       = 0.e0_wp
        exch2       = 0.e0_wp
        exch3       = 0.e0_wp
        g_behf      = 0.e0_wp
        g_belif     = 0.e0_wp
        g_beliff    = 0.e0_wp
        g_belifm    = 0.e0_wp
        g_belifmod  = 0.e0_wp
        g_coul1     = 0.e0_wp
        g_coul2     = 0.e0_wp
        g_coul3     = 0.e0_wp
        g_cplg2     = 0.e0_wp
        g_cs        = 0.e0_wp
        g_cut2      = 0.e0_wp
        g_d1_w      = 0.e0_wp
        g_dehfm     = 0.e0_wp
        g_del2      = 0.e0_wp
        g_delifm    = 0.e0_wp
        g_exch1     = 0.e0_wp
        g_exch2     = 0.e0_wp
        g_exch3     = 0.e0_wp
        g_gc2       = 0.e0_wp
        g_gswt      = 0.e0_wp
        g_r1r3      = 0.e0_wp
        g_r1s       = 0.e0_wp
        g_r2s       = 0.e0_wp
        g_r3s       = 0.e0_wp
        g_rcut2     = 0.e0_wp
        g_rho2      = 0.e0_wp
        g_s1        = 0.e0_wp
        g_s1c       = 0.e0_wp
        g_s2        = 0.e0_wp
        g_s3        = 0.e0_wp
        g_s3a       = 0.e0_wp
        g_s3c       = 0.e0_wp
        g_t1        = 0.e0_wp
        g_t2        = 0.e0_wp
        g_t3        = 0.e0_wp
        g_theta2    = 0.e0_wp
        g_w         = 0.e0_wp
        g_xehfa     = 0.e0_wp
        g_xehfa_180 = 0.e0_wp
        g_xelif     = 0.e0_wp
        g_xelif_180 = 0.e0_wp
        g_xelih     = 0.e0_wp
        g_yhf       = 0.e0_wp
        gam         = 0.e0_wp
        gc2         = 0.e0_wp
        gdel        = 0.e0_wp
        grho        = 0.e0_wp
        gswt        = 0.e0_wp
        r1r3        = 0.e0_wp
        r1r3al      = 0.e0_wp
        r1s         = 0.e0_wp
        r2s         = 0.e0_wp
        r3s         = 0.e0_wp
        rcut2       = 0.e0_wp
        rcut2del    = 0.e0_wp
        rcut2rho    = 0.e0_wp
        rehf        = 0.e0_wp
        rehft       = 0.e0_wp
        relif       = 0.e0_wp
        relih       = 0.e0_wp
        rho2        = 0.e0_wp
        rho2_0      = 0.e0_wp
        rho2_180    = 0.e0_wp
        s1          = 0.e0_wp
        s1c         = 0.e0_wp
        s1c1        = 0.e0_wp
        s1c2        = 0.e0_wp
        s1c3        = 0.e0_wp
        s2          = 0.e0_wp
        s2alp       = 0.e0_wp
        s2rho       = 0.e0_wp
        s3          = 0.e0_wp
        s3a         = 0.e0_wp
        s3c         = 0.e0_wp
        t1          = 0.e0_wp
        t1c1        = 0.e0_wp
        t1c2        = 0.e0_wp
        t2          = 0.e0_wp
        t2a1        = 0.e0_wp
        t2a2        = 0.e0_wp
        t2a3        = 0.e0_wp
        t2a4        = 0.e0_wp
        t3          = 0.e0_wp
        t3a1        = 0.e0_wp
        t3a2        = 0.e0_wp
        t3a3        = 0.e0_wp
        t3a4        = 0.e0_wp
        theta2      = 0.e0_wp
        theta2_0    = 0.e0_wp
        theta2_180  = 0.e0_wp
        uli2p       = 0.e0_wp
        w           = 0.e0_wp
        xehfa       = 0.e0_wp
        xehfa_180   = 0.e0_wp
        xelif       = 0.e0_wp
        xelif_180   = 0.e0_wp
        xelih       = 0.e0_wp
        xelihf_180  = 0.e0_wp
        yhf         = 0.e0_wp
        s1c1        = 4.32258e0_wp
        s1c2        = 7.06452e0_wp
        t1c1        = 1.64516e0_wp
        t1c2        = 10.38710e0_wp
        belihc      = 1.36e0_wp
        bet1        = 2.10667e0_wp
        relih       = 1.2e0_wp
        r1r3al      = 1.0e0_wp
        rcut2rho    = 0.72e0_wp
        rcut2del    = 0.5e0_wp
        s1c3        = 14.74194e0_wp
        bes1        = 0.90667e0_wp
        rehf        = 1.733e0_wp
        dehf        = 6.122e0_wp
        t2a1        = 1.26667e0_wp
        t2a2        = 16.06667e0_wp
        t2a3        = 11.61290e0_wp
        t2a4        = 15.51613e0_wp
        bet2        = 2.10667e0_wp
        bet2_180    = 1.73333e0_wp
        rehft       = 1.733e0_wp
        dehfc       = 0.26e0_wp
        s2rho       = 2.38095e0_wp
        s2alp       = 0.5e0_wp
        belifc1     = 0.95484e0_wp
        belifc2     = 0.91613e0_wp
        delifc1     = 2.58065e0_wp
        delifc2     = 5.25806e0_wp
        belifff     = 0.25333e0_wp
        nb          = 8
        gam         = 3.87097e0_wp
        grho        = 1.41935e0_wp
        gdel        = 3.77419e0_wp
        t3a1        = 0.38710e0_wp
        t3a2        = 1.80645e0_wp
        t3a3        = 0.51613e0_wp
        t3a4        = 0.51613e0_wp
        bet3        = 0.89333e0_wp
        bet3_180    = 0.80e0_wp
        rho2_0      = 0.51613e0_wp
        rho2_180    = 1.33333e0_wp
        del2_0      = 0.49032e0_wp
        del2_180    = 0.31613e0_wp
        theta2_0    = 0.09333e0_wp
        theta2_180  = 0.18e0_wp
        delif       = 5.909e0_wp
        relif       = 2.9553e0_wp
        c2a         = 3.5e0_wp
        c2b         = 0.27362e0_wp
        c2c         = 0.15e0_wp
        uli2p       = 1.848e0_wp
        cevau       = 0.036749309e0_wp
        g_ehfid     = 0
        if( g_p_ .gt. g_pmax_ )then
            print *, 'Parameter g_p_ is greater than g_pmax_'
            stop
        endif
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
            do g_i_ = 1, g_p_
                g_cs(g_i_) = -g_r1s(g_i_) + g_r2s(g_i_) + g_r3s(g_i_)
            enddo
            cs   = r3s + r2s - r1s
            d3_v = 2.0e0_wp * r(i,2)
            d5_v = d3_v * r(i,3)
            d6_v = cs / d5_v
            d2_b = 1.0e0_wp / d5_v
            d3_b = (-d6_v) / d5_v
            d5_b = d3_b * d3_v
            d6_b = d3_b * r(i,3) * 2.0e0_wp
            do g_i_ = 1, g_p_
                g_cs(g_i_) = d5_b * g_r(g_i_,i,3) + d6_b * g_r(g_i_,i,2) + d2_b * g_cs(g_i_)
            enddo
            cs   = d6_v
            d2_v = min (cs, 1.0e0_wp)
            if( cs .le. 1.0e0_wp )then
                d1_p = 1.0e0_wp
                d2_p = 0.0e0_wp
            elseif( cs .gt. 1.0e0_wp )then
                d1_p = 0.0e0_wp
                d2_p = 1.0e0_wp
            endif
            do g_i_ = 1, g_p_
                g_cs(g_i_) = d1_p * g_cs(g_i_)
            enddo
            cs = d2_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = 0.0e0_wp
            enddo
            d1_w = -1.0e0_wp
            d3_v = max (cs, d1_w)
            if( cs .ge. d1_w )then
                d1_p = 1.0e0_wp
                d2_p = 0.0e0_wp
            elseif( cs .lt. d1_w )then
                d1_p = 0.0e0_wp
                d2_p = 1.0e0_wp
            endif
            do g_i_ = 1, g_p_
                g_cs(g_i_) = d2_p * g_d1_w(g_i_) + d1_p * g_cs(g_i_)
            enddo
            cs = d3_v
            do g_i_ = 1, g_p_
                g_r1r3(g_i_) = r1r3al * g_r(g_i_,i,2) + (-g_r(g_i_,i,3)) + g_r(g_i_,i,1)
            enddo
            r1r3 = r(i,1) - r(i,3) + r1r3al * r(i,2)
            d3_b = 1.0e0_wp / rcut2del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r1r3(g_i_)
            enddo
            d1_w = (r1r3 - rcut2rho) / rcut2del
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_rcut2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            rcut2 = 0.5e0_wp * (1.0e0_wp - d2_v)
            ! The diatomic curves . . .
            ! The MH singlet & triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-belihc) * g_r(g_i_,i,1)
            enddo
            d1_w = (-belihc) * (r(i,1) - relih)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xelih(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xelih = d2_v
            d2_v  = xelih * xelih
            d1_p  = 2.0e0_wp * xelih
            d4_b  = s1c2 + s1c1 * d1_p
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d4_b * g_xelih(g_i_)
            enddo
            s1 = s1c1 * d2_v + s1c2 * xelih
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet1) * g_r(g_i_,i,1)
            enddo
            d1_w = (-bet1) * (r(i,1) - relih)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xelih(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xelih = d2_v
            d2_v  = xelih * xelih
            d1_p  = 2.0e0_wp * xelih
            d4_b  = t1c2 + t1c1 * d1_p
            do g_i_ = 1, g_p_
                g_t1(g_i_) = d4_b * g_xelih(g_i_)
            enddo
            t1 = t1c1 * d2_v + t1c2 * xelih
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bes1) * g_r(g_i_,i,1)
            enddo
            d1_w = (-bes1) * (r(i,1) - relih)
            d2_v = exp(d1_w)
            d1_p = d2_v
            d3_b = s1c3 * d1_p
            do g_i_ = 1, g_p_
                g_s1c(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            s1c = s1c3 * d2_v
            do g_i_ = 1, g_p_
                g_s1(g_i_) = rcut2 * g_s1c(g_i_) + s1c * g_rcut2(g_i_) + g_s1(g_i_)
            enddo
            s1 = s1 + rcut2 * s1c
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
            d3_b  = 1.0e0_wp / s2alp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,3)
            enddo
            d1_w = (r(i,3) - s2rho) / s2alp
            d4_v = dehfc * 0.5e0_wp * (1.0e0_wp - cs) * 0.5e0_wp
            d6_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d6_v *  d6_v)
            d7_v = 1.0e0_wp - d6_v
            d6_b = d4_v * d1_p
            d9_b = -((-d7_v) * 0.5e0_wp * (dehfc * 0.5e0_wp))
            do g_i_ = 1, g_p_
                g_dehfm(g_i_) = d6_b * g_d1_w(g_i_) + d9_b * g_cs(g_i_)
            enddo
            dehfm = dehf - d4_v * d7_v
            d3_v  = dehfm * xehfa
            d4_v  = xehfa - 2.0e0_wp
            d5_b  = d4_v * xehfa
            d4_b  = d3_v + d4_v * dehfm
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d4_b * g_xehfa(g_i_) + d5_b * g_dehfm(g_i_)
            enddo
            s2 = d3_v * d4_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet2) * g_r(g_i_,i,2)
            enddo
            d1_w = (-bet2) * (r(i,2) - rehft)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xehfa(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xehfa = d2_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet2_180) * g_r(g_i_,i,2)
            enddo
            d1_w = (-bet2_180) * (r(i,2) - rehft)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xehfa_180(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xehfa_180 = d2_v
            d2_v      = xehfa * xehfa
            d2_p      = 2.0e0_wp * xehfa
            d6_v      = (t2a1 * d2_v + t2a2 * xehfa) * 0.5e0_wp
            d8_v      = 1.0e0_wp + cs
            d11_v     = xehfa_180 * xehfa_180
            d1_p      = 2.0e0_wp * xehfa_180
            d15_v     = (t2a3 * d11_v + t2a4 * xehfa_180) * 0.5e0_wp
            d16_v     = 1.0e0_wp - cs
            d7_b      = d16_v * 0.5e0_wp
            d10_b     = d7_b * t2a4 + d7_b * t2a3 * d1_p
            d6_b      = -d15_v + d6_v
            d14_b     = d8_v * 0.5e0_wp
            d17_b     = d14_b * t2a2 + d14_b * t2a1 * d2_p
            do g_i_ = 1, g_p_
                g_t2(g_i_) = d10_b * g_xehfa_180(g_i_) + d6_b * g_cs(g_i_) + d17_b * g_xehfa(g_i_)
            enddo
            t2 = d6_v * d8_v + d15_v * d16_v
            ! The Modified HF singlet
            ! The  LiF singlet & triplet
            d5_v = (4.6498e0_wp + 7.0489e0_wp * r(i,3)) * (4.6498e0_wp + 7.0489e0_wp * r(i,3))
            d1_p = 2.0e0_wp * (4.6498e0_wp + 7.0489e0_wp * r(i,3))
            d6_v = r(i,3) * 103.57e0_wp / d5_v
            d7_b = (-d6_v) / d5_v * d1_p * 7.0489e0_wp + 1.0e0_wp / d5_v * 103.57e0_wp
            do g_i_ = 1, g_p_
                g_belif(g_i_) = d7_b * g_r(g_i_,i,3)
            enddo
            belif = d6_v + real(0.064076, wp)
            d3_b  = 1.0e0_wp / gdel
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,2)
            enddo
            d1_w = (r(i,2) - grho) / gdel
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_gswt(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            gswt = 0.5e0_wp - 0.5e0_wp * d2_v
            d2_v = gswt * 0.5e0_wp
            d4_v = 1.0e0_wp + cs
            d5_b = d4_v * 0.5e0_wp
            do g_i_ = 1, g_p_
                g_gswt(g_i_) = d2_v * g_cs(g_i_) + d5_b * g_gswt(g_i_)
            enddo
            gswt = d2_v * d4_v
            d3_v = belifff - belif
            d2_b = 1.0e0_wp + (-gswt)
            do g_i_ = 1, g_p_
                g_beliff(g_i_) = d3_v * g_gswt(g_i_) + d2_b * g_belif(g_i_)
            enddo
            beliff = belif + gswt * d3_v
            if( nb .ne. 0 )then
                d5_v = (r(i,3) / gam) ** ( nb-1)
                d2_p = nb *  d5_v
                d5_v = d5_v * (r(i,3) / gam)
            else
                !            Maybe this should be  d5_v = (r(i, 3) / gam) **  nb
                d5_v = 1.0e0_wp
                d2_p = 0.0e0_wp
            endif
            d6_v = belif + d5_v
            if( nb .ne. 0 )then
                d9_v = (r(i,3) / gam) ** ( nb-1)
                d1_p = nb *  d9_v
                d9_v = d9_v * (r(i,3) / gam)
            else
                !            Maybe this should be  d9_v = (r(i, 3) / gam) **  nb
                d9_v = 1.0e0_wp
                d1_p = 0.0e0_wp
            endif
            d10_v = beliff + d9_v
            d11_v = beliff * d6_v / d10_v
            d2_b  = 1.0e0_wp / d10_v
            d3_b  = (-d11_v) / d10_v
            d4_b  = d3_b + d2_b * d6_v
            d8_b  = d2_b * beliff
            d7_b  = d3_b * d1_p * (1.0e0_wp / gam) + d8_b * d2_p * (1.0e0_wp / gam)
            do g_i_ = 1, g_p_
                g_belifmod(g_i_) = d7_b * g_r(g_i_,i,3) + d8_b * g_belif(g_i_) + d4_b * g_beliff(g_i_)
            enddo
            belifmod = d11_v
            d4_v     = r(i,3) - relif
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-belifmod) * g_r(g_i_,i,3) + (-d4_v) * g_belifmod(g_i_)
            enddo
            d1_w = (-belifmod) * d4_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xelif(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xelif = d2_v
            d2_v  = (delif + uli2p) * xelif
            d3_v  = xelif - 2.0e0_wp
            d4_b  = d2_v + d3_v * (delif + uli2p)
            do g_i_ = 1, g_p_
                g_s3a(g_i_) = d4_b * g_xelif(g_i_)
            enddo
            s3a  = d2_v * d3_v
            d5_b = delifc2 * 0.5e0_wp + (-(delifc1 * 0.5e0_wp))
            do g_i_ = 1, g_p_
                g_delifm(g_i_) = d5_b * g_cs(g_i_)
            enddo
            delifm = delifc1 * 0.5e0_wp * (1.0e0_wp - cs) + delifc2 * 0.5e0_wp * (1.0e0_wp + cs)
            d5_b   = belifc2 * 0.5e0_wp + (-(belifc1 * 0.5e0_wp))
            do g_i_ = 1, g_p_
                g_belifm(g_i_) = d5_b * g_cs(g_i_)
            enddo
            belifm = belifc1 * 0.5e0_wp * (1.0e0_wp - cs) + belifc2 * 0.5e0_wp * (1.0e0_wp + cs)
            d4_v   = r(i,3) - relif
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-belifm) * g_r(g_i_,i,3) + (-d4_v) * g_belifm(g_i_)
            enddo
            d1_w = (-belifm) * d4_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xelif(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xelif = d2_v
            d2_v  = delifm + uli2p
            d4_v  = d2_v * xelif
            d5_v  = xelif - 2.0e0_wp
            d4_b  = d4_v + d5_v * d2_v
            d6_b  = d5_v * xelif
            do g_i_ = 1, g_p_
                g_s3c(g_i_) = d4_b * g_xelif(g_i_) + d6_b * g_delifm(g_i_)
            enddo
            s3c = d4_v * d5_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet3) * g_r(g_i_,i,3)
            enddo
            d1_w = (-bet3) * (r(i,3) - relif)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xelif(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xelif = d2_v
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet3_180) * g_r(g_i_,i,3)
            enddo
            d1_w = (-bet3_180) * (r(i,3) - relif)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xelif_180(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xelif_180 = d2_v
            d2_v      = xelif * xelif
            d2_p      = 2.0e0_wp * xelif
            d6_v      = (t3a1 * d2_v + t3a2 * xelif) * 0.5e0_wp
            d8_v      = 1.0e0_wp + cs
            d11_v     = xelif_180 * xelif_180
            d1_p      = 2.0e0_wp * xelif_180
            d15_v     = (t3a3 * d11_v + t3a4 * xelif_180) * 0.5e0_wp
            d16_v     = 1.0e0_wp - cs
            d7_b      = d16_v * 0.5e0_wp
            d10_b     = d7_b * t3a4 + d7_b * t3a3 * d1_p
            d6_b      = -d15_v + d6_v
            d14_b     = d8_v * 0.5e0_wp
            d17_b     = d14_b * t3a2 + d14_b * t3a1 * d2_p
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d10_b * g_xelif_180(g_i_) + d6_b * g_cs(g_i_) + d17_b * g_xelif(g_i_)
            enddo
            t3   = d6_v * d8_v + d15_v * d16_v
            d4_b = -((rho2_180 - rho2_0) * 0.5e0_wp)
            do g_i_ = 1, g_p_
                g_rho2(g_i_) = d4_b * g_cs(g_i_)
            enddo
            rho2 = rho2_0 + (rho2_180 - rho2_0) * 0.5e0_wp * (1.0e0_wp - cs)
            d4_b = -((theta2_180 - theta2_0) * 0.5e0_wp)
            do g_i_ = 1, g_p_
                g_theta2(g_i_) = d4_b * g_cs(g_i_)
            enddo
            theta2 = theta2_0 + (theta2_180 - theta2_0) * 0.5e0_wp * (1.0e0_wp - cs)
            d4_b   = -((del2_180 - del2_0) * 0.5e0_wp)
            do g_i_ = 1, g_p_
                g_del2(g_i_) = d4_b * g_cs(g_i_)
            enddo
            del2 = del2_0 + (del2_180 - del2_0) * 0.5e0_wp * (1.0e0_wp - cs)
            d2_v = cos(theta2)
            d2_p = -sin(theta2)
            d4_v = rehf - r(i,2)
            d6_v = sin(theta2)
            d1_p = cos(theta2)
            d9_v = rho2 - r(i,3)
            d8_b = (-d9_v) * d1_p + d4_v * d2_p
            do g_i_ = 1, g_p_
                g_gc2(g_i_) = d6_v * g_r(g_i_,i,3) + (-d6_v) * g_rho2(g_i_ ) + (-d2_v) * g_r(g_i_,i,2) + d8_b * &
                    & g_theta2(g_i_)
            enddo
            gc2  = d2_v * d4_v - d6_v * d9_v
            d3_v = gc2 / del2
            d2_b = 1.0e0_wp / del2
            d3_b = (-d3_v) / del2
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_del2(g_i_) + d2_b * g_gc2(g_i_)
            enddo
            d1_w = d3_v
            d2_v = tanh (d1_w)
            d1_p = 1.0e0_wp - ( d2_v *  d2_v)
            d4_b = (-0.5e0_wp) * d1_p
            do g_i_ = 1, g_p_
                g_cut2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            cut2 = 0.5e0_wp - 0.5e0_wp * d2_v
            d3_v = s3a - s3c
            d2_b = 1.0e0_wp + (-cut2)
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d3_v * g_cut2(g_i_) + cut2 * g_s3a(g_i_) + d2_b * g_s3c(g_i_)
            enddo
            s3 = s3c + d3_v * cut2
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
                d1_p = 0.0e0_wp ! should never occur
            endif
            d7_b = (-(1.0e0_wp / sqrt(2.0e0_wp))) * d1_p
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d7_b * g_d1_w(g_i_) + g_coul3(g_i_) + g_coul2 (g_i_) + g_coul1(g_i_)
            enddo
            e(i) = coul1 + coul2 + coul3 - d7_v / sqrt(2.0e0_wp) + uli2p + dehf
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
        enddo
    end subroutine dispatch
end module chempotpy_lifh_leps_jhctp_2001_dpem_g_pot22_entry_sta_99e0597c

!> @brief Compatibility interface for g_pot22.
subroutine g_pot22( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_lifh_leps_jhctp_2001_dpem_g_pot22_entry_sta_99e0597c, only: dispatch
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
