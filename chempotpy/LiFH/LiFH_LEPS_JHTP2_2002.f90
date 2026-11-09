!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_lifh_leps_jhtp2_2002_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_legacy_units,       only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: legacy_triatomic_distances_to_bohr
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 2
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        integer                       :: iatom, idir, nder
        double precision              :: nt, r(1,3), r2(3), v(1)
        double precision              :: u11(1), u12(1), u22(1)
        double precision              :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
        double precision              :: u(nstates,nstates), dudr(nstates,nstates,3)
        double precision              :: t(nstates,nstates)
        double precision              :: tmpmat(nstates,nstates)
        double precision              :: hr(nstates,nstates,3), gr(nstates,3)
        double precision              :: tx(9), drdx(3,9)
        double precision              :: hx(nstates,nstates,9), gx(nstates,9)
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
        do i=1, 3
            tmpmat(:,:) = dudr(:,:,i)
            tmpmat      = matmul(transpose(t), matmul(tmpmat, t))
            do j=1, nstates
                do k=j+1, nstates
                    if( (p(k)-p(j)) > 1.d-8 )then
                        hr(j,k,i) = tmpmat(j,k)/(p(k)-p(j))
                    else
                        hr(j,k,i) = tmpmat(j,k)/1.d-8
                    endif
                    hr(k,j,i) = -hr(j,k,i)
                enddo
            enddo
            do j=1, nstates
                gr(j,i) = tmpmat(j,j)
            enddo
        enddo
        r2 = r(1,:)*legacy_bohr_to_angstrom
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
end module chempotpy_lifh_leps_jhtp2_2002_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_lifh_leps_jhtp2_2002_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 2
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX for the LiFH_LEPS_JHTP2_2002 potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    implicit none
    integer, parameter            :: first_atoms(3)  = [1, 2, 1]
    integer, parameter            :: second_atoms(3) = [2, 3, 3]
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!> @brief Evaluate diagonalize for the LiFH_LEPS_JHTP2_2002 potential.
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
!   System:                        LiFH
!   Functional form:               Modified Extended LEPS (2x2 diabatic fit)
!   Common name:                   LiFH surface fit J
!   Number of derivatives:         1
!   Number of electronic surfaces: 2
!   Interface:                     3-2V
!
!   References:      A. W. Jasper, M. D. Hack, D. G. Truhlar, and P. Pie
!                    J. Chem. Phys., Vol. 116, 8353 (2002).
!
!   Notes:    This fit contains long-range forces.
!
!   Protocol: POT evaluates the stateless diabatic surfaces.
!
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!
!   Surfaces:
!      ground and first excited electronic states, diabatic representati
!
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Li
!      infinitely far from the HF diatomic and R(HF) set equal to the
!      HF equilibrium diatomic value.
!
!   Parameters:
!
!   Coordinates:
!      Internal, Definition: R(1) = R(Li-H)
!                            R(2) = R(H-F)
!                            R(3) = R(Li-F)
!
!
!
!
!                   LiFH-J POTENTIAL ENERGY MATRIX
!
!  Created by A. Jasper on July 27, 2001.
!  Based on the LiFH-H potential energy matrix created by M. Hack.
!  Based on the ab initio calculations of P. Piecuch.
!
!  Reference: A. W. Jasper, M. D. Hack, D. G. Truhlar, and P. Piecuch,
!  J. Chem. Phys., prepared for publication.
!

!      subroutine pot(r,e,nt,nsurf)
!> @brief Evaluate pot for the LiFH_LEPS_JHTP2_2002 potential.
subroutine pot( r, e, de, nt, nsurf )
    implicit none
    double precision :: r, e, de, g_r
    integer          :: i, j, k, nt, nsurf
    dimension r(nt,3), e(nt), g_r(3,nt,3), de(3,nt)
    !       seed matrix
    do i=1, nt
        do j=1, 3
            do k=1, 3
                if( j .eq. k ) g_r(j,i,k) = 1.0d0
                if( j .ne. k ) g_r(j,i,k) = 0.0d0
            enddo
        enddo
    enddo
    if( nsurf .eq. 1 ) call pot11(3, r, g_r, 3, e, de, 3, nt)
    if( nsurf .eq. 2 ) call pot12(3, r, g_r, 3, e, de, 3, nt)
    if( nsurf .eq. 3 ) call pot22(3, r, g_r, 3, e, de, 3, nt)
    if( nsurf .gt. 3.or.nsurf .lt. 1 ) &
        & write(6, *)' WARNING:  nsurf = ', nsurf, &
        & ' in LiFHJ potential routine!'
end subroutine pot

!                            U11 SURFACE
!
!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_lifh_leps_jhtp2_2002_pot11_legacy_state_fb8669b7
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
        implicit none
        optional :: g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt
        !
        ! MISC VARIABLES
        integer          :: i, nt, j
        double precision :: cevau
        double precision :: r(:,:), e(:)
        double precision :: r1s, r2s, r3s
        !
        ! INTERACTION REGION VARIABLES
        !     LiH
        double precision :: s1, relih, bes1, s1c2, s1c3, s1a, xelih, bes1b, &
            & s1c, rcut2del, rcut2rho, rcut2, r1r3, r1r3al
        !
        !     HF
        double precision :: s2, behf, yhf, rehf, dehf, xehf, s2m, behfm, dehfm, rehfm, xehfm, s2rho, s2alp, dehfc, dehfcc, yhf_r, behf_0, behf_1, behf_2
        !
        !     LiF
        double precision :: s3, res3, bes3, des3, res3c, bes3c, s3c, s3a, xes3, des3c, rcutdel, rcutrho, rcut, r1r3al2, bes3i, bes3z, bes3r
        !
        ! LONG RANGE FORCES VARIABLES
        !     HF physical properties
        double precision :: hfdip_y, hfdip_a, hfdip_re, hfdip_m(11), hfdip, &
            & hfpol, hfpolperp, hfpolpar, hfie, hfpolavg
        !
        !     Li physical properties
        double precision :: li2spol, li2sie
        !
        !     misc long-range variables
        double precision :: rhf0did, rhf0disp, hf_did, hf_disp, e_lr_hf
        double precision :: mass(3), cauang, caudeb, costheta, temp, bigr, &
            & coschi, r2lrcut_a, r2lrcut_r
        !
        ! MISC DATA
        integer            :: g_p_
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, ldg_r, ldg_e
        double precision   :: d12_b, d12_v, d1_w, d3_p, d2_p, d1_p, d11_b, d2_v, d3_v, d4_v
        double precision   :: d5_v, d6_v, d7_v, d8_v, d9_v, d10_v, d11_v, d9_b, d2_b, d3_b
        double precision   :: d4_b, d5_b, d6_b, d7_b, d8_b, g_costheta(g_pmax_), g_r(:,:,:), g_temp(g_pmax_), g_d1_w(g_pmax_), g_bigr(g_pmax_)
        double precision   :: g_coschi(g_pmax_), g_xelih(g_pmax_), g_s1a(g_pmax_), g_s1c(g_pmax_), g_r1r3(g_pmax_), g_rcut2(g_pmax_), g_s1(g_pmax_), g_yhf(g_pmax_), g_behf(g_pmax_), g_xehf(g_pmax_)
        double precision   :: g_s2(g_pmax_), g_dehfm(g_pmax_), g_xehfm(g_pmax_), g_s2m(g_pmax_), g_xes3(g_pmax_), g_s3(g_pmax_), g_e(:,:), g_hfdip_y(g_pmax_), g_hfdip(g_pmax_), g_hfpol(g_pmax_)
        double precision   :: g_hf_did(g_pmax_), g_hf_disp(g_pmax_), g_e_lr_hf(g_pmax_)
        integer            :: g_ehfid
        !
        ! INTERACTION REGION DATA
        !     LiH
        !
        !     HF
        !
        !     LiF
        !
        ! LONG RANGE FORCES DATA
        !     HF dipole moment (angstroms vs au)
        !
        !     HF physical properties (polarization in au, ie in eV)
        !
        !     Li(2s) physical properties (polarization in au, ie in eV)
        !
        !     misc long-range parameters (bohr)
        !
        !
        !
        behf       = 0.d0
        behf_0     = 0.d0
        behf_1     = 0.d0
        behf_2     = 0.d0
        behfm      = 0.d0
        bes1       = 0.d0
        bes1b      = 0.d0
        bes3       = 0.d0
        bes3c      = 0.d0
        bes3i      = 0.d0
        bes3r      = 0.d0
        bes3z      = 0.d0
        bigr       = 0.d0
        cauang     = 0.d0
        caudeb     = 0.d0
        cevau      = 0.d0
        coschi     = 0.d0
        costheta   = 0.d0
        d10_v      = 0.d0
        d11_b      = 0.d0
        d11_v      = 0.d0
        d12_b      = 0.d0
        d12_v      = 0.d0
        d1_p       = 0.d0
        d1_w       = 0.d0
        d2_b       = 0.d0
        d2_p       = 0.d0
        d2_v       = 0.d0
        d3_b       = 0.d0
        d3_p       = 0.d0
        d3_v       = 0.d0
        d4_b       = 0.d0
        d4_v       = 0.d0
        d5_b       = 0.d0
        d5_v       = 0.d0
        d6_b       = 0.d0
        d6_v       = 0.d0
        d7_b       = 0.d0
        d7_v       = 0.d0
        d8_b       = 0.d0
        d8_v       = 0.d0
        d9_b       = 0.d0
        d9_v       = 0.d0
        dehf       = 0.d0
        dehfc      = 0.d0
        dehfcc     = 0.d0
        dehfm      = 0.d0
        des3       = 0.d0
        des3c      = 0.d0
        e_lr_hf    = 0.d0
        g_behf     = 0.d0
        g_bigr     = 0.d0
        g_coschi   = 0.d0
        g_costheta = 0.d0
        g_d1_w     = 0.d0
        g_dehfm    = 0.d0
        g_e_lr_hf  = 0.d0
        g_hf_did   = 0.d0
        g_hf_disp  = 0.d0
        g_hfdip    = 0.d0
        g_hfdip_y  = 0.d0
        g_hfpol    = 0.d0
        g_r1r3     = 0.d0
        g_rcut2    = 0.d0
        g_s1       = 0.d0
        g_s1a      = 0.d0
        g_s1c      = 0.d0
        g_s2       = 0.d0
        g_s2m      = 0.d0
        g_s3       = 0.d0
        g_temp     = 0.d0
        g_xehf     = 0.d0
        g_xehfm    = 0.d0
        g_xelih    = 0.d0
        g_xes3     = 0.d0
        g_yhf      = 0.d0
        hf_did     = 0.d0
        hf_disp    = 0.d0
        hfdip      = 0.d0
        hfdip_a    = 0.d0
        hfdip_m    = 0.d0
        hfdip_re   = 0.d0
        hfdip_y    = 0.d0
        hfie       = 0.d0
        hfpol      = 0.d0
        hfpolavg   = 0.d0
        hfpolpar   = 0.d0
        hfpolperp  = 0.d0
        li2sie     = 0.d0
        li2spol    = 0.d0
        mass       = 0.d0
        r1r3       = 0.d0
        r1r3al     = 0.d0
        r1r3al2    = 0.d0
        r1s        = 0.d0
        r2lrcut_a  = 0.d0
        r2lrcut_r  = 0.d0
        r2s        = 0.d0
        r3s        = 0.d0
        rcut       = 0.d0
        rcut2      = 0.d0
        rcut2del   = 0.d0
        rcut2rho   = 0.d0
        rcutdel    = 0.d0
        rcutrho    = 0.d0
        rehf       = 0.d0
        rehfm      = 0.d0
        relih      = 0.d0
        res3       = 0.d0
        res3c      = 0.d0
        rhf0did    = 0.d0
        rhf0disp   = 0.d0
        s1         = 0.d0
        s1a        = 0.d0
        s1c        = 0.d0
        s1c2       = 0.d0
        s1c3       = 0.d0
        s2         = 0.d0
        s2alp      = 0.d0
        s2m        = 0.d0
        s2rho      = 0.d0
        s3         = 0.d0
        s3a        = 0.d0
        s3c        = 0.d0
        temp       = 0.d0
        xehf       = 0.d0
        xehfm      = 0.d0
        xelih      = 0.d0
        xes3       = 0.d0
        yhf        = 0.d0
        yhf_r      = 0.d0
        cevau      = 0.036749309d0
        relih      = 1.3432062d0
        bes1       = 1.8058651d0
        s1c2       = 12.937438d0
        s1c3       = 11.0913d0
        bes1b      = 1.2111436d0
        rcut2rho   = 1.07143d0
        rcut2del   = 0.6d0
        r1r3al     = 0.4d0
        rehf       = 1.733d0
        dehf       = 6.122d0
        rehfm      = 1.6739d0
        dehfc      = 0.344574d0
        dehfcc     = 0.8660801d0
        s2rho      = 1.63d0
        s2alp      = 2.9941d0
        behfm      = 0.74633431d0
        yhf_r      = 2.1042d0
        behf_0     = 1.1622d0
        behf_1     = 0.025647d0
        behf_2     = 0.059062d0
        res3       = 3.6d0
        bes3       = 1.49061584d0
        des3       = 0.17426687d0
        res3c      = 2.48571d0
        bes3c      = 1.92857d0
        des3c      = 0.25d0
        rcutrho    = 0.4d0
        rcutdel    = 0.89333d0
        r1r3al2    = 0.7333d0
        bes3z      = 0.d0
        bes3r      = 7.4799608d0
        hfdip_a    = 1.13d0
        hfdip_re   = 0.924212d0
        hfdip_m = [ 0.703661d0, 0.516815d0, 0.240628d0, -0.430194d0, -2.213088d0, 5.535609d0, 12.872106d0, - &
            & 42.060172d0, -17.398065d0, 84.207629d0, -41.961465d0]
        hfpolperp = 4.59d0
        hfpolpar  = 5.10d0
        hfie      = 16.044d0
        li2spol   = 165.d0
        li2sie    = 5.392d0
        rhf0did   = 6.0d0
        rhf0disp  = 6.0d0
        r2lrcut_a = 2.0d0
        r2lrcut_r = 3.d0
        cauang    = 0.52917706d0
        caudeb    = 2.54177d0
        mass      = [ 7.016003d0, 1.00783d0, 18.9984d0]
        !
        ! ----------------------------------------------------------------------
        !
        !
        do i = 1, nt
            !
            ! PRECOMPUTE THINGS FOR USE LATER...
            !     cosine of HF bond angle (i.e., Li-F-H angle)
            d2_v  = r(i,3) * r(i,3)
            d3_p  = 2.0d0 * r(i,3)
            d4_v  = r(i,2) * r(i,2)
            d2_p  = 2.0d0 * r(i,2)
            d7_v  = r(i,1) * r(i,1)
            d1_p  = 2.0d0 * r(i,1)
            d9_v  = (-2.d0) * r(i,1)
            d10_v = d9_v * r(i,2)
            d11_v = (d2_v - d4_v - d7_v) / d10_v
            d2_b  = 1.0d0 / d10_v
            d3_b  = (-d11_v) / d10_v
            d6_b  = d3_b * r(i,2) * (-2.d0) + (-d2_b) * d1_p
            d5_b  = d3_b * d9_v + (-d2_b) * d2_p
            d11_b = d2_b * d3_p
            do g_i_ = 1, g_p_
                g_costheta(g_i_) = d6_b * g_r(g_i_,i,1) + d5_b * g_r(g_i_,i,2) + d11_b * g_r(g_i_,i,3)
            enddo
            costheta = d11_v
            !--------
            d2_b = mass(3) / (mass(2) + mass(3))
            do g_i_ = 1, g_p_
                g_temp(g_i_) = d2_b * g_r(g_i_,i,2)
            enddo
            temp = mass(3) / (mass(2) + mass(3)) * r(i,2)
            !--------
            !     Li-FH Jacobi distance
            d2_v = r(i,1) * r(i,1)
            d2_p = 2.0d0 * r(i,1)
            d4_v = temp * temp
            d1_p = 2.0d0 * temp
            d6_v = 2.d0 * r(i,1)
            d7_v = d6_v * temp
            d7_b = (-costheta) * d6_v + d1_p
            d8_b = (-costheta) * temp * 2.d0 + d2_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-d7_v) * g_costheta(g_i_) + d7_b * g_temp(g_i_) + d8_b * g_r(g_i_,i,1)
            enddo
            d1_w = d2_v + d4_v - d7_v * costheta
            d2_v = sqrt(d1_w)
            !          if ( d1_w .gt. 0.0d0 ) then
            d1_p = 1.0d0 / (2.0d0 *  d2_v)
            !          else
            !             call ehufDO (9,d1_w, d2_v, d1_p,
            !     +g_ehfid,
            !     +185)
            !          endif
            do g_i_ = 1, g_p_
                g_bigr(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            bigr = d2_v
            !--------
            !     cosine of HF Jacobi angle (i.e., Li--HF angle)
            d2_v  = r(i,1) * r(i,1)
            d3_p  = 2.0d0 * r(i,1)
            d4_v  = bigr * bigr
            d2_p  = 2.0d0 * bigr
            d7_v  = temp * temp
            d1_p  = 2.0d0 * temp
            d10_v = (-2.d0) * bigr
            d11_v = d10_v * temp
            d12_v = (-(d2_v - d4_v - d7_v)) / d11_v
            d3_b  = (-d12_v) / d11_v
            d7_b  = -(1.0d0 / d11_v)
            d5_b  = d3_b * d10_v + (-d7_b) * d1_p
            d6_b  = d3_b * temp * (-2.d0) + (-d7_b) * d2_p
            d12_b = d7_b * d3_p
            do g_i_ = 1, g_p_
                g_coschi(g_i_) = d5_b * g_temp(g_i_) + d6_b * g_bigr(g_i_) + d12_b * g_r(g_i_,i,1)
            enddo
            coschi = d12_v
            !--------
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = 0.0d0
            enddo
            d1_w = -1.d0
            d3_v = max (coschi, d1_w)
            !          if (coschi .gt.  d1_w) then
            d1_p = 1.0d0
            d2_p = 0.0d0
            !          else if (coschi .lt.  d1_w) then
            !             d1_p = 0.0d0
            !             d2_p = 1.0d0
            !          else
            !             call ehbfDO (7,coschi, d1_w, d3_v, d1_p, d2_p,
            !     +g_ehfid,
            !     +228)
            !             d2_p = 1.0d0 -  d1_p
            !          endif
            do g_i_ = 1, g_p_
                g_coschi(g_i_) = d2_p * g_d1_w(g_i_) + d1_p * g_coschi(g_i_)
            enddo
            coschi = d3_v
            !--------
            d2_v = min (coschi, 1.d0)
            !          if (coschi .lt.  1.d0) then
            d1_p = 1.0d0
            d2_p = 0.0d0
            !          else if (coschi .gt.  1.d0) then
            !             d1_p = 0.0d0
            !             d2_p = 1.0d0
            !          else
            !             call ehbfDO (8,coschi, 1.d0, d2_v, d1_p, d2_p,
            !     +g_ehfid,
            !     +247)
            !             d2_p = 1.0d0 -  d1_p
            !          endif
            do g_i_ = 1, g_p_
                g_coschi(g_i_) = d1_p * g_coschi(g_i_)
            enddo
            coschi = d2_v
            !--------
            !     squared!
            r1s = r(i,1) ** 2
            r2s = r(i,2) ** 2
            r3s = r(i,3) ** 2
            !
            ! INTERACTION REGION
            !     LiH
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
            !--------
            do g_i_ = 1, g_p_
                g_s1a(g_i_) = s1c2 * g_xelih(g_i_)
            enddo
            s1a = s1c2 * xelih
            !--------
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
            s1c = s1c3 * d2_v
            !--------
            !
            do g_i_ = 1, g_p_
                g_r1r3(g_i_) = r1r3al * g_r(g_i_,i,2) + (-g_r(g_i_,i,3)) + g_r(g_i_,i,1)
            enddo
            r1r3 = r(i,1) - r(i,3) + r1r3al * r(i,2)
            !--------
            d3_b = 1.0d0 / rcut2del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r1r3(g_i_)
            enddo
            d1_w = (r1r3 - rcut2rho) / rcut2del
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_rcut2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            rcut2 = 0.5d0 * (1.0d0 - d2_v)
            !--------
            !
            d3_v = s1c - s1a
            d2_b = 1.0d0 + (-rcut2)
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d3_v * g_rcut2(g_i_) + rcut2 * g_s1c(g_i_) + d2_b * g_s1a(g_i_)
            enddo
            s1 = s1a + d3_v * rcut2
            !--------
            !
            !     HF
            do g_i_ = 1, g_p_
                g_yhf(g_i_) = g_r(g_i_,i,2)
            enddo
            yhf = r(i,2) - yhf_r
            !--------
            d4_v = yhf * yhf
            d1_p = 2.0d0 * yhf
            d5_b = behf_2 * d1_p + (-behf_1)
            do g_i_ = 1, g_p_
                g_behf(g_i_) = d5_b * g_yhf(g_i_)
            enddo
            behf = behf_0 - behf_1 * yhf + behf_2 * d4_v
            !--------
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
            !--------
            d2_v = dehf * xehf
            d3_v = xehf - 2.0d0
            d4_b = d2_v + d3_v * dehf
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d4_b * g_xehf(g_i_)
            enddo
            s2 = d2_v * d3_v
            !--------
            !
            d5_b = dehfc * 0.5d0
            do g_i_ = 1, g_p_
                g_dehfm(g_i_) = d5_b * g_coschi(g_i_)
            enddo
            dehfm = dehf - (dehfcc + dehfc * 0.5d0 * (1.d0 - coschi))
            !--------
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-behfm) * g_r(g_i_,i,2)
            enddo
            d1_w = (-behfm) * (r(i,2) - rehfm)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xehfm(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xehfm = d2_v
            !--------
            d3_v = dehfm * xehfm
            d4_v = xehfm - 2.0d0
            d5_b = d4_v * xehfm
            d4_b = d3_v + d4_v * dehfm
            do g_i_ = 1, g_p_
                g_s2m(g_i_) = d4_b * g_xehfm(g_i_) + d5_b * g_dehfm(g_i_)
            enddo
            s2m = d3_v * d4_v
            !--------
            d3_b = 1.0d0 / s2alp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,3)
            enddo
            d1_w = (r(i,3) - s2rho) / s2alp
            d4_v = (s2 - s2m) * 0.5d0
            d6_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d6_v *  d6_v)
            d7_v = 1.d0 + d6_v
            d7_b = d4_v * d1_p
            d8_b = d7_v * 0.5d0
            d2_b = 1.0d0 + (-d8_b)
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s2(g_i_) + d2_b * g_s2m(g_i_)
            enddo
            s2 = s2m + d4_v * d7_v
            !--------
            !
            !     LiF
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
            !--------
            d2_v = des3 * xes3
            d3_v = xes3 - 2.0d0
            d4_b = d2_v + d3_v * des3
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d4_b * g_xes3(g_i_)
            enddo
            s3 = d2_v * d3_v
            !--------
            !
            !     Sum of diatomics
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = g_s3(g_i_) + g_s2(g_i_) + g_s1(g_i_)
            enddo
            e(i) = s1 + s2 + s3 + dehf
            !--------
            !
            ! EVALUATE THE LONG-RANGE FORCES
            ! HF DIPOLE MOMENT
            d4_b = (-hfdip_a) * cauang
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_r(g_i_,i,2)
            enddo
            d1_w = (-hfdip_a) * (r(i,2) * cauang - hfdip_re)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_hfdip_y(g_i_) = (-d1_p) * g_d1_w(g_i_)
            enddo
            hfdip_y = 1.d0 - d2_v
            !--------
            do g_i_ = 1, g_p_
                g_hfdip(g_i_) = 0.0d0
            enddo
            hfdip = 0.d0
            !--------
            do j = 1, 11
                d3_v = dble(j - 1)
                !            if ( hfdip_y .ne. 0.0d0 ) then
                d4_v = hfdip_y ** ( d3_v - 2.0d0)
                d4_v = d4_v * hfdip_y
                d1_p = d3_v *  d4_v
                d4_v = d4_v * hfdip_y
                !            else
                !              (hfdip_y = 0)
                !               d4_v = hfdip_y **  d3_v
                !
                !               if (  d3_v .lt. 1.0d0 ) then
                !                  call ehbfDO (10,hfdip_y, d3_v, d4_v, d1_p, 0.0d0,
                !     +g_ehfid,
                !     +463)
                !               else if (  d3_v .lt. 2.0d0 ) then
                !                  d1_p = 0.0d0
                !                  call ehbfDO (10,hfdip_y, d3_v, d4_v, d1_p, 0.0d0,
                !     +g_ehfid,
                !     +468)
                !               else
                !                  d1_p = 0.0d0
                !               endif
                !            endif
                d6_b = 0.0d0
                d5_b = hfdip_m(j) * d1_p
                do g_i_ = 1, g_p_
                    g_hfdip(g_i_) = d5_b * g_hfdip_y(g_i_) + g_hfdip(g_i_)
                enddo
                hfdip = hfdip + hfdip_m(j) * d4_v
                !--------
            enddo
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = 20.d0 * g_r(g_i_,i,2)
            enddo
            d1_w = 20.d0 * (r(i,2) - 1.2d0)
            d2_v = hfdip * 0.5d0
            d4_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d4_v *  d4_v)
            d5_v = 1.d0 + d4_v
            d5_b = d2_v * d1_p
            d6_b = d5_v * 0.5d0
            do g_i_ = 1, g_p_
                g_hfdip(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_hfdip(g_i_)
            enddo
            hfdip = d2_v * d5_v
            !--------
            !
            ! HF POLARIZATION
            d2_v = coschi * coschi
            d2_p = 2.0d0 * coschi
            d5_v = coschi * coschi
            d1_p = 2.0d0 * coschi
            d5_b = hfpolpar * d1_p + (-hfpolperp) * d2_p
            do g_i_ = 1, g_p_
                g_hfpol(g_i_) = d5_b * g_coschi(g_i_)
            enddo
            hfpol = hfpolperp * (1.d0 - d2_v) + hfpolpar * d5_v
            !--------
            hfpolavg = 2.d0 * hfpolperp / 3.d0 + hfpolpar / 3.d0
            !
            ! HF + Li(2s) arrangement
            !     DIPOLE-INDUCED-DIPOLE
            d2_v = coschi * coschi
            d2_p = 2.0d0 * coschi
            d6_v = (-0.5d0) * (3.d0 * d2_v + 1.d0) * li2spol
            d8_v = hfdip * hfdip
            d1_p = 2.0d0 * hfdip
            d4_b = d6_v * d1_p
            d9_b = d8_v * li2spol * (-0.5d0) * 3.d0 * d2_p
            do g_i_ = 1, g_p_
                g_hf_did(g_i_) = d4_b * g_hfdip(g_i_) + d9_b * g_coschi(g_i_ )
            enddo
            hf_did = d6_v * d8_v
            !--------
            d3_v = bigr ** ( 7 - 2)
            d3_v = d3_v * bigr
            d1_p = 7 *  d3_v
            d3_v = d3_v * bigr
            d4_v = d3_v + rhf0did ** 7
            d5_v = hf_did / d4_v
            d6_v = d5_v / cevau
            d4_b = bigr * (1.0d0 / cevau)
            d5_b = d4_b * (1.0d0 / d4_v)
            d3_b = d6_v + d4_b * ((-d5_v) / d4_v) * d1_p
            do g_i_ = 1, g_p_
                g_hf_did(g_i_) = d3_b * g_bigr(g_i_) + d5_b * g_hf_did(g_i_)
            enddo
            hf_did = d6_v * bigr
            !--------
            !     cut off angular dependence in interaction region
            d3_v = (bigr / rhf0did) ** ( 4 - 2)
            d3_v = d3_v * (bigr / rhf0did)
            d1_p = 4 *  d3_v
            d3_v = d3_v * (bigr / rhf0did)
            d4_b = (-d1_p) * (1.0d0 / rhf0did)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_bigr(g_i_)
            enddo
            d1_w  = -d3_v
            d3_v  = coschi * coschi
            d2_p  = 2.0d0 * coschi
            d5_v  = 3.d0 * d3_v + 1.d0
            d6_v  = 2.5d0 / d5_v
            d7_v  = d6_v - 1.d0
            d9_v  = exp(d1_w)
            d1_p  = d9_v
            d11_v = 1.d0 + d7_v * d9_v
            d7_b  = hf_did * d7_v * d1_p
            d12_b = hf_did * d9_v * ((-d6_v) / d5_v) * 3.d0 * d2_p
            do g_i_ = 1, g_p_
                g_hf_did(g_i_) = d7_b * g_d1_w(g_i_) + d12_b * g_coschi(g_i_ ) + d11_v * g_hf_did(g_i_)
            enddo
            hf_did = hf_did * d11_v
            !--------
            !
            !     DISPERSION
            d3_b = 1.0d0 / (li2sie + hfie) * ((-1.5d0) * li2sie * hfie * li2spol)
            do g_i_ = 1, g_p_
                g_hf_disp(g_i_) = d3_b * g_hfpol(g_i_)
            enddo
            hf_disp = (-1.5d0) * li2sie * hfie * li2spol * hfpol / (li2sie + hfie)
            !--------
            d3_v = bigr ** ( 7 - 2)
            d3_v = d3_v * bigr
            d1_p = 7 *  d3_v
            d3_v = d3_v * bigr
            d4_v = d3_v + rhf0disp ** 7
            d5_v = hf_disp / d4_v
            d4_b = bigr * (1.0d0 / d4_v)
            d3_b = d5_v + bigr * ((-d5_v) / d4_v) * d1_p
            do g_i_ = 1, g_p_
                g_hf_disp(g_i_) = d3_b * g_bigr(g_i_) + d4_b * g_hf_disp(g_i_)
            enddo
            hf_disp = d5_v * bigr
            !--------
            !     cut off angular dependence in interaction region
            d3_v = (bigr / rhf0disp) ** ( 4 - 2)
            d3_v = d3_v * (bigr / rhf0disp)
            d1_p = 4 *  d3_v
            d3_v = d3_v * (bigr / rhf0disp)
            d4_b = (-d1_p) * (1.0d0 / rhf0disp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_bigr(g_i_)
            enddo
            d1_w = -d3_v
            d3_v = hfpolavg / hfpol
            d4_v = d3_v - 1.d0
            d6_v = exp(d1_w)
            d1_p = d6_v
            d8_v = 1.d0 + d4_v * d6_v
            d7_b = hf_disp * d4_v * d1_p
            d9_b = hf_disp * d6_v * ((-d3_v) / hfpol)
            do g_i_ = 1, g_p_
                g_hf_disp(g_i_) = d7_b * g_d1_w(g_i_) + d9_b * g_hfpol(g_i_) + d8_v * g_hf_disp(g_i_)
            enddo
            hf_disp = hf_disp * d8_v
            !--------
            !
            !     SUM FORCES
            do g_i_ = 1, g_p_
                g_e_lr_hf(g_i_) = g_hf_did(g_i_) + g_hf_disp(g_i_)
            enddo
            e_lr_hf = hf_disp + hf_did
            !--------
            !     cut off LR forces for large r[hf] i.e., r(i,2)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = r2lrcut_a * g_r(g_i_,i,2)
            enddo
            d1_w = r2lrcut_a * (r(i,2) - r2lrcut_r)
            d2_v = e_lr_hf * 0.5d0
            d4_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d4_v *  d4_v)
            d5_v = 1.d0 - d4_v
            d5_b = (-d2_v) * d1_p
            d6_b = d5_v * 0.5d0
            do g_i_ = 1, g_p_
                g_e_lr_hf(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_e_lr_hf(g_i_)
            enddo
            e_lr_hf = d2_v * d5_v
            !--------
            !
            ! LiF + H
            !     NO LONG-RANGE FORCES
            !
            ! LiH + F
            !     NO LONG-RANGE FORCES
            !
            ! COMBINE INTERACTION REGION ENERGY AND LONG-RANGE FORCES
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = g_e_lr_hf(g_i_) + g_e(g_i_,i)
            enddo
            e(i) = e(i) + e_lr_hf
            !--------
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
            !--------
            !
        enddo
        !
    end subroutine dispatch
end module chempotpy_lifh_leps_jhtp2_2002_pot11_legacy_state_fb8669b7

!> @brief Compatibility interface for pot11.
subroutine pot11( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_lifh_leps_jhtp2_2002_pot11_legacy_state_fb8669b7, only: dispatch
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
end subroutine pot11

!

!                      U12 COUPLING SURFACE
!
!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_lifh_leps_jhtp2_2002_pot12_legacy_state_fb8669b7
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
        !
        implicit none
        optional :: g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt
        !
        ! MISC VARIABLES
        integer          :: i, nt
        double precision :: r(:,:), e(:)
        !
        ! LiH VARIABLES
        integer          :: ng1
        double precision :: g1, g1r, g1a, h1, h1del, h1rho
        !
        ! HF VARIABLES
        double precision :: j2p, j2, j2pdel, j2rho, j2prho, j2del
        !
        ! LiF VARIABLES
        integer          :: ng3
        double precision :: g3, g3r, g3a, h3, h3del, h3rho
        !
        ! CUTOFF VARIABLES
        double precision :: gcrho, cevau, gcdel
        double precision :: cut2, gc
        !
        ! MISC DATA
        integer            :: g_p_
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, ldg_r, ldg_e
        double precision   :: d13_b, d12_b, d11_v, d2_p, d9_v, d8_b, d1_w, d7_b, d2_v, d3_v
        double precision   :: d4_v, d1_p, d2_b, d3_b, d4_b, d5_v, d6_v, g_d1_w(g_pmax_), g_r(:,:,:), g_g1(g_pmax_)
        double precision   :: g_g3(g_pmax_), g_j2(g_pmax_), g_j2p(g_pmax_), g_h1(g_pmax_), g_h3(g_pmax_), g_e(:,:), g_gc(g_pmax_), g_cut2 &
            & (g_pmax_)
        integer :: g_ehfid
        !
        ! LiH DATA
        !
        ! LiF DATA
        !
        ! HF DATA
        !
        ! CUTOFF DATA
        !
        !
        !
        cevau  = 0.d0
        cut2   = 0.d0
        d11_v  = 0.d0
        d12_b  = 0.d0
        d13_b  = 0.d0
        d1_p   = 0.d0
        d1_w   = 0.d0
        d2_b   = 0.d0
        d2_p   = 0.d0
        d2_v   = 0.d0
        d3_b   = 0.d0
        d3_v   = 0.d0
        d4_b   = 0.d0
        d4_v   = 0.d0
        d5_v   = 0.d0
        d6_v   = 0.d0
        d7_b   = 0.d0
        d8_b   = 0.d0
        d9_v   = 0.d0
        g1     = 0.d0
        g1a    = 0.d0
        g1r    = 0.d0
        g3     = 0.d0
        g3a    = 0.d0
        g3r    = 0.d0
        g_cut2 = 0.d0
        g_d1_w = 0.d0
        g_g1   = 0.d0
        g_g3   = 0.d0
        g_gc   = 0.d0
        g_h1   = 0.d0
        g_h3   = 0.d0
        g_j2   = 0.d0
        g_j2p  = 0.d0
        gc     = 0.d0
        gcdel  = 0.d0
        gcrho  = 0.d0
        h1     = 0.d0
        h1del  = 0.d0
        h1rho  = 0.d0
        h3     = 0.d0
        h3del  = 0.d0
        h3rho  = 0.d0
        j2     = 0.d0
        j2del  = 0.d0
        j2p    = 0.d0
        j2pdel = 0.d0
        j2prho = 0.d0
        j2rho  = 0.d0
        cevau  = 0.036749309d0
        g1a    = 1.27742d0
        g1r    = 3.0048876d0
        ng1    = 6
        h1rho  = 4.9843597d0
        h1del  = 2.0899315d0
        g3a    = 0.48d0
        g3r    = 3.4975562d0
        ng3    = 8
        h3rho  = 2.000978d0
        h3del  = 0.90615835d0
        j2rho  = 1.15484d0
        j2del  = 1.75806d0
        j2prho = 1.45161d0
        j2pdel = 0.98387d0
        gcrho  = 4.87097d0
        gcdel  = 2.d0
        !
        ! ----------------------------------------------------------------------
        !
        !
        do i = 1, nt
            !
            d4_b = (-dble(ng1)) * (1.0d0 / g1r)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_r(g_i_,i,1)
            enddo
            d1_w = (-dble(ng1)) * (r(i,1) / g1r - 1.0d0)
            if( ng1 .ne. 0 )then
                d3_v = (r(i,1) / g1r) ** ( ng1-1)
                d2_p = ng1 *  d3_v
                d3_v = d3_v * (r(i,1) / g1r)
            else
                !            Maybe this should be  d3_v = (r(i, 1) / g1r) **  ng1
                d3_v = 1.0d0
                d2_p = 0.0d0
            endif
            d4_v = g1a * d3_v
            d6_v = exp(d1_w)
            d1_p = d6_v
            d4_b = d4_v * d1_p
            d7_b = d6_v * g1a * d2_p * (1.0d0 / g1r)
            do g_i_ = 1, g_p_
                g_g1(g_i_) = d4_b * g_d1_w(g_i_) + d7_b * g_r(g_i_,i,1)
            enddo
            g1 = d4_v * d6_v
            !--------
            !
            d4_b = (-dble(ng3)) * (1.0d0 / g3r)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_r(g_i_,i,3)
            enddo
            d1_w = (-dble(ng3)) * (r(i,3) / g3r - 1.0d0)
            if( ng3 .ne. 0 )then
                d3_v = (r(i,3) / g3r) ** ( ng3-1)
                d2_p = ng3 *  d3_v
                d3_v = d3_v * (r(i,3) / g3r)
            else
                !            Maybe this should be  d3_v = (r(i, 3) / g3r) **  ng3
                d3_v = 1.0d0
                d2_p = 0.0d0
            endif
            d4_v = g3a * d3_v
            d6_v = exp(d1_w)
            d1_p = d6_v
            d4_b = d4_v * d1_p
            d7_b = d6_v * g3a * d2_p * (1.0d0 / g3r)
            do g_i_ = 1, g_p_
                g_g3(g_i_) = d4_b * g_d1_w(g_i_) + d7_b * g_r(g_i_,i,3)
            enddo
            g3 = d4_v * d6_v
            !--------
            !
            d3_b = 1.0d0 / j2del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,2)
            enddo
            d1_w = (r(i,2) - j2rho) / j2del
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = 0.5d0 * d1_p
            do g_i_ = 1, g_p_
                g_j2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            j2 = 0.5d0 * (1.0d0 + d2_v)
            !--------
            d3_b = 1.0d0 / j2pdel
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,2)
            enddo
            d1_w = (r(i,2) - j2prho) / j2pdel
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = 0.5d0 * d1_p
            do g_i_ = 1, g_p_
                g_j2p(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            j2p = 0.5d0 * (1.0d0 + d2_v)
            !--------
            !
            d3_b = 1.0d0 / h1del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,1)
            enddo
            d1_w = (r(i,1) - h1rho) / h1del
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_h1(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            h1 = 0.5d0 * (1.0d0 - d2_v)
            !--------
            d3_b = 1.0d0 / h3del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,3)
            enddo
            d1_w = (r(i,3) - h3rho) / h3del
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_h3(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            h3 = 0.5d0 * (1.0d0 - d2_v)
            !--------
            !
            d3_v  = g1 * j2
            d5_v  = 1.0d0 - h3
            d9_v  = g3 * j2p
            d11_v = 1.0d0 - h1
            d7_b  = d11_v * j2p
            d8_b  = d11_v * g3
            d12_b = d5_v * j2
            d13_b = d5_v * g1
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = (-d9_v) * g_h1(g_i_) + d8_b * g_j2p(g_i_) + d7_b * g_g3(g_i_) + (-d3_v) * g_h3(g_i_) + &
                    & d13_b * g_j2(g_i_) + d12_b * g_g1(g_i_)
            enddo
            e(i) = d3_v * d5_v + d9_v * d11_v
            !--------
            !
            do g_i_ = 1, g_p_
                g_gc(g_i_) = g_r(g_i_,i,2)
            enddo
            gc = r(i,2) - gcrho
            !--------
            d2_b = 1.0d0 / gcdel
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d2_b * g_gc(g_i_)
            enddo
            d1_w = gc / gcdel
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_cut2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            cut2 = 0.5d0 - 0.5d0 * d2_v
            !--------
            !
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = e(i) * g_cut2(g_i_) + cut2 * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cut2
            !--------
            !
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
            !--------
            !
        enddo
    end subroutine dispatch
end module chempotpy_lifh_leps_jhtp2_2002_pot12_legacy_state_fb8669b7

!> @brief Compatibility interface for pot12.
subroutine pot12( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_lifh_leps_jhtp2_2002_pot12_legacy_state_fb8669b7, only: dispatch
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
end subroutine pot12

!

!                            U22 SURFACE
!
!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_lifh_leps_jhtp2_2002_pot22_legacy_state_fb8669b7
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
        implicit none
        optional :: g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt
        !
        ! MISC VARIABLES
        integer          :: i, j, k, nt
        double precision :: r(:,:), e(:), r1s, r2s, r3s, cs, cevau, uli2p
        double precision :: c2a, c2b, c2c, w, cplg2
        !
        ! INTERACTION REGION VARIABLES
        !     LiH
        double precision :: r1r3, xelih, s1c, bes1, rcut2del, rcut2rho, s1c3
        double precision :: s1c1, belihc, bet1, relih, s1c2, t1c1, t1c2
        double precision :: r1r3al, rcut2
        double precision :: coul1, exch1, s1, t1, reliht
        !
        !     HF
        double precision :: behf, xehfa, dehfm, yhf, xehfa_180, t2a3, t2a2, &
            & bet2
        double precision :: t2a1, t2a4, dehf, rehf
        double precision :: s2alp, rehft, bet2_180, s2rho, dehfc
        double precision :: coul2, exch2, t2, s2, dehfcc, s2cut, xehfm, s2m, s2del, s2cc, s2cc_d, s2cc_r, s2cc_b
        !
        !     LiF
        integer          :: nb
        double precision :: bet3, bet3_180, t3a2, t3a3, t3a4, belihc1
        double precision :: belifff, gam, belifc2, delifc1, delifc2, belifm
        double precision :: s3c, xelihf_180, xelif, s3a, delifm, gc2, cut2
        double precision :: del2_0, del2_180, rho2_0, rho2_180, theta2_0, grho
        double precision :: theta2_180, beliff, belifmod, belif, xelif_180, &
            & relif
        double precision :: belifc1, t3a1, gswt, delif, rho2, del2, theta2, &
            & gdel
        double precision :: coul3, exch3, s3, t3, s3cut, swit, swit_r, swit_a, swit_d
        !
        ! LONG RANGE FORCES VARIABLES
        !     HF physical properties
        double precision :: hfdip_y, hfdip_a, hfdip_re, hfdip_m(11), hfdip, &
            & hfpol, hfpolperp, hfpolpar, hfie, hfpolavg
        double precision :: hfqua_y, hfqua_a, hfqua_re, hfqua_m(10), hfqua
        !
        !     Li(2p) physical properties
        double precision :: li2ppol, li2pie, li2pqua, li2ppolpar, li2ppolperp, li2ppolavg
        !
        !     LiF physical properties
        double precision :: lifpol, lifie
        double precision :: lifdip_y, lifdip_re, lifdip_1, lifdip_2, lifdip_3, lifdip, lifdip_4, lifdip_beta, lifdip_d
        !
        !     H physical properties
        double precision :: hie, hpol
        !
        !     misc long-range variables
        double precision :: rhf0did, rhf0disp, hf_disp, e_lr_hf, e_lr_lif, &
            & rlif0did, rlif0disp, lif_disp, lif_did, hf_qq, hf_dq, rhf0qq, rhf0dq, hf_did
        double precision :: mass(3), cauang, caudeb, costhetahf, costhetalif, bigrhf, bigrlif, coschihf, coschilif, temp, r2lrcut_a, r2lrcut_r, r3lrcut_a, r3lrcut_r
        !
        ! MISC DATA
        integer            :: g_p_
        integer, parameter :: g_pmax_ = 3
        integer            :: g_i_, ldg_r, ldg_e
        double precision   :: d1_w, d17_b, d12_v, d12_b, d14_b, d13_v, d16_v, &
            & d2_v, d3_v, d4_v
        double precision :: d5_v, d6_v, d7_v, d8_v, d9_v, d10_v, d11_v, d15_v, d2_b, d3_b
        double precision :: d4_b, d5_b, d6_b, d7_b, d8_b, d9_b, d10_b, d11_b, d1_p, d2_p
        double precision :: d3_p, g_costhetahf(g_pmax_), g_r(:,:,:), &
            & g_temp(g_pmax_), g_d1_w(g_pmax_), g_bigrhf(g_pmax_), g_coschihf(g_pmax_), g_costhetalif(g_pmax_), g_bigrlif(g_pmax_), g_coschilif(g_pmax_)
        double precision :: g_r1s(g_pmax_), g_r2s(g_pmax_), g_r3s(g_pmax_), &
            & g_cs(g_pmax_), g_xelih(g_pmax_), g_s1(g_pmax_), g_s1c(g_pmax_), g_r1r3(g_pmax_), g_rcut2(g_pmax_), g_t1(g_pmax_)
        double precision :: g_yhf(g_pmax_), g_behf(g_pmax_), g_xehfa(g_pmax_), g_s2(g_pmax_), g_s2cut(g_pmax_), g_s2cc(g_pmax_), g_dehfm(g_pmax_), g_xehfm(g_pmax_), g_s2m(g_pmax_), g_xehfa_180(g_pmax_)
        double precision :: g_t2(g_pmax_), g_belif(g_pmax_), g_gswt(g_pmax_ &
            & ), g_beliff(g_pmax_), g_belifmod(g_pmax_), g_xelif(g_pmax_), g_s3a &
            & (g_pmax_), g_s3cut(g_pmax_), g_swit(g_pmax_), g_delifm(g_pmax_)
        double precision :: g_belifm(g_pmax_), g_s3c(g_pmax_), g_rho2(g_pmax_), g_theta2(g_pmax_), g_del2(g_pmax_), g_gc2(g_pmax_), g_cut2(g_pmax_), g_s3(g_pmax_), g_xelif_180(g_pmax_), g_t3(g_pmax_)
        double precision :: g_coul1(g_pmax_), g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_), g_exch3(g_pmax_), g_w(g_pmax_), g_cplg2(g_pmax_), g_e(:,:), g_hfdip_y(g_pmax_)
        double precision :: g_hfdip(g_pmax_), g_hfpol(g_pmax_), g_hfqua_y(g_pmax_), g_hfqua(g_pmax_), g_lifdip_y(g_pmax_), g_lifdip_beta(g_pmax_), g_lifdip(g_pmax_), g_li2ppol(g_pmax_), g_hf_did(g_pmax_), g_hf_disp(g_pmax_)
        double precision :: g_hf_qq(g_pmax_), g_hf_dq(g_pmax_), g_e_lr_hf(g_pmax_), g_lif_did(g_pmax_), g_lif_disp(g_pmax_), g_e_lr_lif(g_pmax_)
        integer          :: g_ehfid
        !
        ! INTERACTION REGION DATA
        !     LiH
        !
        !     HF
        !
        !     LiF
        !
        ! LONG RANGE FORCES DATA
        !     HF dipole moment
        !
        !     HF quadrupole moment
        !
        !     HF pol (au) and ie (eV)
        !
        !     Li(2p) pol (au) and ie (eV) and quadrupole moment (au)
        !
        !     LiF dipole moment (distance and moment in au)
        !
        !     LiF pol (au) and ie (eV)
        !
        !     H pol (au) and ie (eV)
        !
        !     misc long-range
        !      data rlif0did/7.5d0/,rlif0disp/7.5d0/,r2LRcut_a/2.0d0/,
        !
        !        data g_ehfid /0/
        !
        !
        behf          = 0.d0
        belif         = 0.d0
        belifc1       = 0.d0
        belifc2       = 0.d0
        beliff        = 0.d0
        belifff       = 0.d0
        belifm        = 0.d0
        belifmod      = 0.d0
        belihc        = 0.d0
        belihc1       = 0.d0
        bes1          = 0.d0
        bet1          = 0.d0
        bet2          = 0.d0
        bet2_180      = 0.d0
        bet3          = 0.d0
        bet3_180      = 0.d0
        bigrhf        = 0.d0
        bigrlif       = 0.d0
        c2a           = 0.d0
        c2b           = 0.d0
        c2c           = 0.d0
        cauang        = 0.d0
        caudeb        = 0.d0
        cevau         = 0.d0
        coschihf      = 0.d0
        coschilif     = 0.d0
        costhetahf    = 0.d0
        costhetalif   = 0.d0
        coul1         = 0.d0
        coul2         = 0.d0
        coul3         = 0.d0
        cplg2         = 0.d0
        cs            = 0.d0
        cut2          = 0.d0
        d10_b         = 0.d0
        d10_v         = 0.d0
        d11_b         = 0.d0
        d11_v         = 0.d0
        d12_b         = 0.d0
        d12_v         = 0.d0
        d13_v         = 0.d0
        d14_b         = 0.d0
        d15_v         = 0.d0
        d16_v         = 0.d0
        d17_b         = 0.d0
        d1_p          = 0.d0
        d1_w          = 0.d0
        d2_b          = 0.d0
        d2_p          = 0.d0
        d2_v          = 0.d0
        d3_b          = 0.d0
        d3_p          = 0.d0
        d3_v          = 0.d0
        d4_b          = 0.d0
        d4_v          = 0.d0
        d5_b          = 0.d0
        d5_v          = 0.d0
        d6_b          = 0.d0
        d6_v          = 0.d0
        d7_b          = 0.d0
        d7_v          = 0.d0
        d8_b          = 0.d0
        d8_v          = 0.d0
        d9_b          = 0.d0
        d9_v          = 0.d0
        dehf          = 0.d0
        dehfc         = 0.d0
        dehfcc        = 0.d0
        dehfm         = 0.d0
        del2          = 0.d0
        del2_0        = 0.d0
        del2_180      = 0.d0
        delif         = 0.d0
        delifc1       = 0.d0
        delifc2       = 0.d0
        delifm        = 0.d0
        e_lr_hf       = 0.d0
        e_lr_lif      = 0.d0
        exch1         = 0.d0
        exch2         = 0.d0
        exch3         = 0.d0
        g_behf        = 0.d0
        g_belif       = 0.d0
        g_beliff      = 0.d0
        g_belifm      = 0.d0
        g_belifmod    = 0.d0
        g_bigrhf      = 0.d0
        g_bigrlif     = 0.d0
        g_coschihf    = 0.d0
        g_coschilif   = 0.d0
        g_costhetahf  = 0.d0
        g_costhetalif = 0.d0
        g_coul1       = 0.d0
        g_coul2       = 0.d0
        g_coul3       = 0.d0
        g_cplg2       = 0.d0
        g_cs          = 0.d0
        g_cut2        = 0.d0
        g_d1_w        = 0.d0
        g_dehfm       = 0.d0
        g_del2        = 0.d0
        g_delifm      = 0.d0
        g_e_lr_hf     = 0.d0
        g_e_lr_lif    = 0.d0
        g_exch1       = 0.d0
        g_exch2       = 0.d0
        g_exch3       = 0.d0
        g_gc2         = 0.d0
        g_gswt        = 0.d0
        g_hf_did      = 0.d0
        g_hf_disp     = 0.d0
        g_hf_dq       = 0.d0
        g_hf_qq       = 0.d0
        g_hfdip       = 0.d0
        g_hfdip_y     = 0.d0
        g_hfpol       = 0.d0
        g_hfqua       = 0.d0
        g_hfqua_y     = 0.d0
        g_li2ppol     = 0.d0
        g_lif_did     = 0.d0
        g_lif_disp    = 0.d0
        g_lifdip      = 0.d0
        g_lifdip_beta = 0.d0
        g_lifdip_y    = 0.d0
        g_r1r3        = 0.d0
        g_r1s         = 0.d0
        g_r2s         = 0.d0
        g_r3s         = 0.d0
        g_rcut2       = 0.d0
        g_rho2        = 0.d0
        g_s1          = 0.d0
        g_s1c         = 0.d0
        g_s2          = 0.d0
        g_s2cc        = 0.d0
        g_s2cut       = 0.d0
        g_s2m         = 0.d0
        g_s3          = 0.d0
        g_s3a         = 0.d0
        g_s3c         = 0.d0
        g_s3cut       = 0.d0
        g_swit        = 0.d0
        g_t1          = 0.d0
        g_t2          = 0.d0
        g_t3          = 0.d0
        g_temp        = 0.d0
        g_theta2      = 0.d0
        g_w           = 0.d0
        g_xehfa       = 0.d0
        g_xehfa_180   = 0.d0
        g_xehfm       = 0.d0
        g_xelif       = 0.d0
        g_xelif_180   = 0.d0
        g_xelih       = 0.d0
        g_yhf         = 0.d0
        gam           = 0.d0
        gc2           = 0.d0
        gdel          = 0.d0
        grho          = 0.d0
        gswt          = 0.d0
        hf_did        = 0.d0
        hf_disp       = 0.d0
        hf_dq         = 0.d0
        hf_qq         = 0.d0
        hfdip         = 0.d0
        hfdip_a       = 0.d0
        hfdip_m       = 0.d0
        hfdip_re      = 0.d0
        hfdip_y       = 0.d0
        hfie          = 0.d0
        hfpol         = 0.d0
        hfpolavg      = 0.d0
        hfpolpar      = 0.d0
        hfpolperp     = 0.d0
        hfqua         = 0.d0
        hfqua_a       = 0.d0
        hfqua_m       = 0.d0
        hfqua_re      = 0.d0
        hfqua_y       = 0.d0
        hie           = 0.d0
        hpol          = 0.d0
        li2pie        = 0.d0
        li2ppol       = 0.d0
        li2ppolavg    = 0.d0
        li2ppolpar    = 0.d0
        li2ppolperp   = 0.d0
        li2pqua       = 0.d0
        lif_did       = 0.d0
        lif_disp      = 0.d0
        lifdip        = 0.d0
        lifdip_1      = 0.d0
        lifdip_2      = 0.d0
        lifdip_3      = 0.d0
        lifdip_4      = 0.d0
        lifdip_beta   = 0.d0
        lifdip_d      = 0.d0
        lifdip_re     = 0.d0
        lifdip_y      = 0.d0
        lifie         = 0.d0
        lifpol        = 0.d0
        mass          = 0.d0
        r1r3          = 0.d0
        r1r3al        = 0.d0
        r1s           = 0.d0
        r2lrcut_a     = 0.d0
        r2lrcut_r     = 0.d0
        r2s           = 0.d0
        r3lrcut_a     = 0.d0
        r3lrcut_r     = 0.d0
        r3s           = 0.d0
        rcut2         = 0.d0
        rcut2del      = 0.d0
        rcut2rho      = 0.d0
        rehf          = 0.d0
        rehft         = 0.d0
        relif         = 0.d0
        relih         = 0.d0
        reliht        = 0.d0
        rhf0did       = 0.d0
        rhf0disp      = 0.d0
        rhf0dq        = 0.d0
        rhf0qq        = 0.d0
        rho2          = 0.d0
        rho2_0        = 0.d0
        rho2_180      = 0.d0
        rlif0did      = 0.d0
        rlif0disp     = 0.d0
        s1            = 0.d0
        s1c           = 0.d0
        s1c1          = 0.d0
        s1c2          = 0.d0
        s1c3          = 0.d0
        s2            = 0.d0
        s2alp         = 0.d0
        s2cc          = 0.d0
        s2cc_b        = 0.d0
        s2cc_d        = 0.d0
        s2cc_r        = 0.d0
        s2cut         = 0.d0
        s2del         = 0.d0
        s2m           = 0.d0
        s2rho         = 0.d0
        s3            = 0.d0
        s3a           = 0.d0
        s3c           = 0.d0
        s3cut         = 0.d0
        swit          = 0.d0
        swit_a        = 0.d0
        swit_d        = 0.d0
        swit_r        = 0.d0
        t1            = 0.d0
        t1c1          = 0.d0
        t1c2          = 0.d0
        t2            = 0.d0
        t2a1          = 0.d0
        t2a2          = 0.d0
        t2a3          = 0.d0
        t2a4          = 0.d0
        t3            = 0.d0
        t3a1          = 0.d0
        t3a2          = 0.d0
        t3a3          = 0.d0
        t3a4          = 0.d0
        temp          = 0.d0
        theta2        = 0.d0
        theta2_0      = 0.d0
        theta2_180    = 0.d0
        uli2p         = 0.d0
        w             = 0.d0
        xehfa         = 0.d0
        xehfa_180     = 0.d0
        xehfm         = 0.d0
        xelif         = 0.d0
        xelif_180     = 0.d0
        xelih         = 0.d0
        xelihf_180    = 0.d0
        yhf           = 0.d0
        c2a           = 3.5d0
        c2b           = 0.27362d0
        c2c           = 0.15d0
        uli2p         = 1.848d0
        cevau         = 0.036749309d0
        s1c1          = 4.32258d0
        s1c2          = 7.06452d0
        t1c1          = 1.1935483d0
        t1c2          = 13.548387d0
        belihc        = 1.36d0
        bet1          = 2.41319648d0
        relih         = 1.2d0
        r1r3al        = 1.0d0
        rcut2rho      = 0.72d0
        rcut2del      = 0.5d0
        s1c3          = 14.74194d0
        bes1          = 0.90667d0
        reliht        = 1.203323d0
        rehf          = 1.733d0
        dehf          = 6.122d0
        t2a1          = 1.4242424d0
        t2a2          = 14.203323d0
        t2a3          = 6.15835777d0
        t2a4          = 4.56304985d0
        bet2          = 2.4105572d0
        bet2_180      = 1.4046921d0
        rehft         = 1.61329423d0
        s2rho         = 1.00293d0
        s2alp         = 4.0899d0
        dehfc         = 0.0486803519d0
        dehfcc        = 0.72629521d0
        s2cc_b        = 1.799609d0
        s2cc_r        = 1.60215d0
        s2cc_d        = 2.3841642d0
        belifc1       = 0.819648d0
        belifc2       = 0.7795699d0
        delifc1       = 2.352884d0
        delifc2       = 5.5904203d0
        belifff       = 0.13225806d0
        nb            = 8
        gam           = 4.6304985d0
        grho          = 3.382209d0
        gdel          = 0.49472140d0
        t3a1          = -.0552298d0
        t3a2          = 1.8729228d0
        t3a3          = 0.94662756d0
        t3a4          = 0.2994134d0
        bet3          = 0.9519061d0
        bet3_180      = 0.4531769d0
        rho2_0        = 0.2017595d0
        rho2_180      = 2.d0
        del2_0        = 0.51710655d0
        del2_180      = 0.43695014d0
        theta2_0      = 0.12463343d0
        theta2_180    = 0.14907135d0
        delif         = 5.947d0
        relif         = 2.9553d0
        swit_r        = 3.340762d0
        swit_a        = 0.56353861d0
        swit_d        = 1.0243998d0
        hfdip_a       = 1.13d0
        hfdip_re      = 0.924212d0
        hfdip_m = [ 0.703661d0, 0.516815d0, 0.240628d0, -0.430194d0, -2.213088d0, 5.535609d0, 12.872106d0, - &
            & 42.060172d0, -17.398065d0, 84.207629d0, -41.961465d0]
        hfqua_a  = 1.37d0
        hfqua_re = 0.9168d0
        hfqua_m = [ 1.6527d0, 2.072d0, 2.0521d0, 3.545d0, -9.8031d0, -19.3808d0, 76.4423d0, 0.d0, -161.5896d0, &
            & 105.7460d0]
        hfpolperp   = 4.59d0
        hfpolpar    = 5.10d0
        hfie        = 16.044d0
        li2ppolpar  = 131.d0
        li2ppolperp = 129.d0
        li2pie      = 3.544d0
        li2pqua     = 11.1d0
        lifdip_re   = 10.4994d0
        lifdip_1    = 0.02435d0
        lifdip_2    = 9.9355d-8
        lifdip_3    = 0.0015999d0
        lifdip_4    = 4.471959d0
        lifdip_d    = 9.30039d0
        lifpol      = 72.9d0
        lifie       = 11.3d0
        hie         = 13.598d0
        hpol        = 4.4997d0
        rhf0did     = 6.d0
        rhf0disp    = 6.d0
        rhf0qq      = 6.d0
        rhf0dq      = 6.d0
        rlif0did    = 7.d0
        rlif0disp   = 7.d0
        r2lrcut_a   = 2.0d0
        r2lrcut_r   = 3.d0
        r3lrcut_a   = 2.d0
        r3lrcut_r   = 4.5d0
        cauang      = 0.52917706d0
        caudeb      = 2.54177d0
        mass        = [ 7.016003d0, 1.00783d0, 18.9984d0]
        !
        ! ----------------------------------------------------------------------
        !
        !
        do i = 1, nt
            !
            ! PRECOMPUTE THINGS FOR USE LATER
            !     cosine of HF bond angle (i.e., Li-F-H angle)
            d2_v  = r(i,3) * r(i,3)
            d3_p  = 2.0d0 * r(i,3)
            d4_v  = r(i,2) * r(i,2)
            d2_p  = 2.0d0 * r(i,2)
            d7_v  = r(i,1) * r(i,1)
            d1_p  = 2.0d0 * r(i,1)
            d9_v  = (-2.d0) * r(i,1)
            d10_v = d9_v * r(i,2)
            d11_v = (d2_v - d4_v - d7_v) / d10_v
            d2_b  = 1.0d0 / d10_v
            d3_b  = (-d11_v) / d10_v
            d6_b  = d3_b * r(i,2) * (-2.d0) + (-d2_b) * d1_p
            d5_b  = d3_b * d9_v + (-d2_b) * d2_p
            d11_b = d2_b * d3_p
            do g_i_ = 1, g_p_
                g_costhetahf(g_i_) = d6_b * g_r(g_i_,i,1) + d5_b * g_r(g_i_,i,2) + d11_b * g_r(g_i_,i,3)
            enddo
            costhetahf = d11_v
            !--------
            d2_b = mass(3) / (mass(2) + mass(3))
            do g_i_ = 1, g_p_
                g_temp(g_i_) = d2_b * g_r(g_i_,i,2)
            enddo
            temp = mass(3) / (mass(2) + mass(3)) * r(i,2)
            !--------
            !     Li-FH Jacobi distance
            d2_v = r(i,1) * r(i,1)
            d2_p = 2.0d0 * r(i,1)
            d4_v = temp * temp
            d1_p = 2.0d0 * temp
            d6_v = 2.d0 * r(i,1)
            d7_v = d6_v * temp
            d7_b = (-costhetahf) * d6_v + d1_p
            d8_b = (-costhetahf) * temp * 2.d0 + d2_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-d7_v) * g_costhetahf(g_i_) + d7_b * g_temp( g_i_) + d8_b * g_r(g_i_,i,1)
            enddo
            d1_w = d2_v + d4_v - d7_v * costhetahf
            d2_v = sqrt(d1_w)
            !          if ( d1_w .gt. 0.0d0 ) then
            d1_p = 1.0d0 / (2.0d0 *  d2_v)
            !          else
            !             call ehufDO (9,d1_w, d2_v, d1_p,
            !     +g_ehfid,
            !     +278)
            !          endif
            do g_i_ = 1, g_p_
                g_bigrhf(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            bigrhf = d2_v
            !--------
            !     cosine of HF Jacobi angle (i.e., Li--HF angle)
            d2_v  = r(i,1) * r(i,1)
            d3_p  = 2.0d0 * r(i,1)
            d4_v  = bigrhf * bigrhf
            d2_p  = 2.0d0 * bigrhf
            d7_v  = temp * temp
            d1_p  = 2.0d0 * temp
            d10_v = (-2.d0) * bigrhf
            d11_v = d10_v * temp
            d12_v = (-(d2_v - d4_v - d7_v)) / d11_v
            d3_b  = (-d12_v) / d11_v
            d7_b  = -(1.0d0 / d11_v)
            d5_b  = d3_b * d10_v + (-d7_b) * d1_p
            d6_b  = d3_b * temp * (-2.d0) + (-d7_b) * d2_p
            d12_b = d7_b * d3_p
            do g_i_ = 1, g_p_
                g_coschihf(g_i_) = d5_b * g_temp(g_i_) + d6_b * g_bigrhf(g_i_) + d12_b * g_r(g_i_,i,1)
            enddo
            coschihf = d12_v
            !--------
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = 0.0d0
            enddo
            d1_w = -1.d0
            d3_v = max (d1_w, coschihf)
            !          if (d1_w .gt.  coschihf) then
            !             d1_p = 1.0d0
            !             d2_p = 0.0d0
            !          else if (d1_w .lt.  coschihf) then
            d1_p = 0.0d0
            d2_p = 1.0d0
            !          else
            !             call ehbfDO (7,d1_w, coschihf, d3_v, d1_p, d2_p,
            !     +g_ehfid,
            !     +321)
            !             d2_p = 1.0d0 -  d1_p
            !          endif
            do g_i_ = 1, g_p_
                g_coschihf(g_i_) = d2_p * g_coschihf(g_i_) + d1_p * g_d1_w(g_i_)
            enddo
            coschihf = d3_v
            !--------
            d2_v = min (1.d0, coschihf)
            !          if (1.d0 .lt.  coschihf) then
            !             d1_p = 1.0d0
            !             d2_p = 0.0d0
            !          else if (1.d0 .gt.  coschihf) then
            d1_p = 0.0d0
            d2_p = 1.0d0
            !          else
            !             call ehbfDO (8,1.d0, coschihf, d2_v, d1_p, d2_p,
            !     +g_ehfid,
            !     +341)
            !             d2_p = 1.0d0 -  d1_p
            !          endif
            do g_i_ = 1, g_p_
                g_coschihf(g_i_) = d2_p * g_coschihf(g_i_)
            enddo
            coschihf = d2_v
            !--------
            !
            !     cosine of LiF bond angle (i.e., Li-F-H angle)
            d2_v  = r(i,1) * r(i,1)
            d3_p  = 2.0d0 * r(i,1)
            d4_v  = r(i,2) * r(i,2)
            d2_p  = 2.0d0 * r(i,2)
            d7_v  = r(i,3) * r(i,3)
            d1_p  = 2.0d0 * r(i,3)
            d9_v  = (-2.d0) * r(i,3)
            d10_v = d9_v * r(i,2)
            d11_v = (d2_v - d4_v - d7_v) / d10_v
            d2_b  = 1.0d0 / d10_v
            d3_b  = (-d11_v) / d10_v
            d6_b  = d3_b * r(i,2) * (-2.d0) + (-d2_b) * d1_p
            d5_b  = d3_b * d9_v + (-d2_b) * d2_p
            d11_b = d2_b * d3_p
            do g_i_ = 1, g_p_
                g_costhetalif(g_i_) = d6_b * g_r(g_i_,i,3) + d5_b * g_r(g_i_,i,2) + d11_b * g_r(g_i_,i,1)
            enddo
            costhetalif = d11_v
            !--------
            d2_b = mass(1) / (mass(1) + mass(3))
            do g_i_ = 1, g_p_
                g_temp(g_i_) = d2_b * g_r(g_i_,i,3)
            enddo
            temp = mass(1) / (mass(1) + mass(3)) * r(i,3)
            !--------
            !     H-LiF Jacobi distance
            d2_v = r(i,2) * r(i,2)
            d2_p = 2.0d0 * r(i,2)
            d4_v = temp * temp
            d1_p = 2.0d0 * temp
            d6_v = 2.d0 * r(i,2)
            d7_v = d6_v * temp
            d7_b = (-costhetalif) * d6_v + d1_p
            d8_b = (-costhetalif) * temp * 2.d0 + d2_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-d7_v) * g_costhetalif(g_i_) + d7_b * g_temp (g_i_) + d8_b * g_r(g_i_,i,2)
            enddo
            d1_w = d2_v + d4_v - d7_v * costhetalif
            d2_v = sqrt(d1_w)
            !          if ( d1_w .gt. 0.0d0 ) then
            d1_p = 1.0d0 / (2.0d0 *  d2_v)
            !          else
            !             call ehufDO (9,d1_w, d2_v, d1_p,
            !     +g_ehfid,
            !     +398)
            !          endif
            do g_i_ = 1, g_p_
                g_bigrlif(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            bigrlif = d2_v
            !--------
            !     cosine of LiF Jacobi angle (i.e., H--LiF angle)
            d2_v  = r(i,2) * r(i,2)
            d3_p  = 2.0d0 * r(i,2)
            d4_v  = bigrlif * bigrlif
            d2_p  = 2.0d0 * bigrlif
            d7_v  = temp * temp
            d1_p  = 2.0d0 * temp
            d10_v = (-2.d0) * bigrlif
            d11_v = d10_v * temp
            d12_v = (-(d2_v - d4_v - d7_v)) / d11_v
            d3_b  = (-d12_v) / d11_v
            d7_b  = -(1.0d0 / d11_v)
            d5_b  = d3_b * d10_v + (-d7_b) * d1_p
            d6_b  = d3_b * temp * (-2.d0) + (-d7_b) * d2_p
            d12_b = d7_b * d3_p
            do g_i_ = 1, g_p_
                g_coschilif(g_i_) = d5_b * g_temp(g_i_) + d6_b * g_bigrlif(g_i_) + d12_b * g_r(g_i_,i,2)
            enddo
            coschilif = d12_v
            !--------
            d2_v = min (1.d0, coschilif)
            !             if (1.d0 .lt.  coschilif) then
            !                d1_p = 1.0d0
            !                d2_p = 0.0d0
            !             else if (1.d0 .gt.  coschilif) then
            d1_p = 0.0d0
            d2_p = 1.0d0
            !             else
            !                call ehbfDO (8,1.d0, coschilif, d2_v, d1_p, d2_p,
            !     +g_ehfid,
            !     +437)
            !                   d2_p = 1.0d0 -  d1_p
            !                endif
            do g_i_ = 1, g_p_
                g_coschilif(g_i_) = d2_p * g_coschilif(g_i_)
            enddo
            coschilif = d2_v
            !--------
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = 0.0d0
            enddo
            d1_w = -1.d0
            d3_v = max (d1_w, coschilif)
            !             if (d1_w .gt.  coschilif) then
            !                d1_p = 1.0d0
            !                d2_p = 0.0d0
            !             else if (d1_w .lt.  coschilif) then
            d1_p = 0.0d0
            d2_p = 1.0d0
            !             else
            !                call ehbfDO (7,d1_w, coschilif, d3_v, d1_p, d2_p,
            !     +g_ehfid,
            !     +460)
            !                   d2_p = 1.0d0 -  d1_p
            !                endif
            do g_i_ = 1, g_p_
                g_coschilif(g_i_) = d2_p * g_coschilif(g_i_) + d1_p * g_d1_w (g_i_)
            enddo
            coschilif = d3_v
            !--------
            !     squared!
            d2_v = r(i,1) * r(i,1)
            d1_p = 2.0d0 * r(i,1)
            do g_i_ = 1, g_p_
                g_r1s(g_i_) = d1_p * g_r(g_i_,i,1)
            enddo
            r1s = d2_v
            !--------
            d2_v = r(i,2) * r(i,2)
            d1_p = 2.0d0 * r(i,2)
            do g_i_ = 1, g_p_
                g_r2s(g_i_) = d1_p * g_r(g_i_,i,2)
            enddo
            r2s = d2_v
            !--------
            d2_v = r(i,3) * r(i,3)
            d1_p = 2.0d0 * r(i,3)
            do g_i_ = 1, g_p_
                g_r3s(g_i_) = d1_p * g_r(g_i_,i,3)
            enddo
            r3s = d2_v
            !--------
            !     Li-F-H bondangle
            do g_i_ = 1, g_p_
                g_cs(g_i_) = -g_r1s(g_i_) + g_r2s(g_i_) + g_r3s(g_i_)
            enddo
            cs = r3s + r2s - r1s
            !--------
            d3_v = 2.0d0 * r(i,2)
            d5_v = d3_v * r(i,3)
            d6_v = cs / d5_v
            d2_b = 1.0d0 / d5_v
            d3_b = (-d6_v) / d5_v
            d5_b = d3_b * d3_v
            d6_b = d3_b * r(i,3) * 2.0d0
            do g_i_ = 1, g_p_
                g_cs(g_i_) = d5_b * g_r(g_i_,i,3) + d6_b * g_r(g_i_,i,2) + d2_b * g_cs(g_i_)
            enddo
            cs = d6_v
            !--------
            d2_v = min (cs, 1.0d0)
            !          if (cs .lt.  1.0d0) then
            d1_p = 1.0d0
            d2_p = 0.0d0
            !          else if (cs .gt.  1.0d0) then
            !             d1_p = 0.0d0
            !             d2_p = 1.0d0
            !          else
            !             call ehbfDO (8,cs, 1.0d0, d2_v, d1_p, d2_p,
            !     +g_ehfid,
            !     +521)
            !             d2_p = 1.0d0 -  d1_p
            !          endif
            do g_i_ = 1, g_p_
                g_cs(g_i_) = d1_p * g_cs(g_i_)
            enddo
            cs = d2_v
            !--------
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = 0.0d0
            enddo
            d1_w = -1.0d0
            d3_v = max (cs, d1_w)
            !          if (cs .gt.  d1_w) then
            d1_p = 1.0d0
            d2_p = 0.0d0
            !          else if (cs .lt.  d1_w) then
            !             d1_p = 0.0d0
            !             d2_p = 1.0d0
            !          else
            !             call ehbfDO (7,cs, d1_w, d3_v, d1_p, d2_p,
            !     +g_ehfid,
            !     +544)
            !             d2_p = 1.0d0 -  d1_p
            !          endif
            do g_i_ = 1, g_p_
                g_cs(g_i_) = d2_p * g_d1_w(g_i_) + d1_p * g_cs(g_i_)
            enddo
            cs = d3_v
            !--------
            !
            ! INTERACTION REGION
            !     The LiH singlet
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
            !--------
            d2_v = xelih * xelih
            d1_p = 2.0d0 * xelih
            d4_b = s1c2 + s1c1 * d1_p
            do g_i_ = 1, g_p_
                g_s1(g_i_) = d4_b * g_xelih(g_i_)
            enddo
            s1 = s1c1 * d2_v + s1c2 * xelih
            !--------
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
            !--------
            do g_i_ = 1, g_p_
                g_r1r3(g_i_) = r1r3al * g_r(g_i_,i,2) + (-g_r(g_i_,i,3)) + g_r(g_i_,i,1)
            enddo
            r1r3 = r(i,1) - r(i,3) + r1r3al * r(i,2)
            !--------
            d3_b = 1.0d0 / rcut2del
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r1r3(g_i_)
            enddo
            d1_w = (r1r3 - rcut2rho) / rcut2del
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_rcut2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            rcut2 = 0.5d0 * (1.0d0 - d2_v)
            !--------
            do g_i_ = 1, g_p_
                g_s1(g_i_) = rcut2 * g_s1c(g_i_) + s1c * g_rcut2(g_i_) + g_s1(g_i_)
            enddo
            s1 = s1 + rcut2 * s1c
            !--------
            !
            !     The LiH triplet
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-bet1) * g_r(g_i_,i,1)
            enddo
            d1_w = (-bet1) * (r(i,1) - reliht)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xelih(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xelih = d2_v
            !--------
            d2_v = xelih * xelih
            d1_p = 2.0d0 * xelih
            d4_b = t1c2 + t1c1 * d1_p
            do g_i_ = 1, g_p_
                g_t1(g_i_) = d4_b * g_xelih(g_i_)
            enddo
            t1 = t1c1 * d2_v + t1c2 * xelih
            !--------
            !
            !     The HF singlet
            do g_i_ = 1, g_p_
                g_yhf(g_i_) = g_r(g_i_,i,2)
            enddo
            yhf = r(i,2) - 2.1042d0
            !--------
            d4_v = yhf * yhf
            d1_p = 2.0d0 * yhf
            d5_b = 0.059062d0 * d1_p + (-0.025647d0)
            do g_i_ = 1, g_p_
                g_behf(g_i_) = d5_b * g_yhf(g_i_)
            enddo
            behf = 1.1622d0 - 0.025647d0 * yhf + 0.059062d0 * d4_v
            !--------
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
            !--------
            !
            d2_v = dehf * xehfa
            d3_v = xehfa - 2.0d0
            d5_b = d2_v + d3_v * dehf
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d5_b * g_xehfa(g_i_)
            enddo
            s2 = d2_v * d3_v + uli2p
            !--------
            d3_b = 1.0d0 / 0.5d0
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,2)
            enddo
            d1_w = (r(i,2) - 3.d0) / 0.5d0
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_s2cut(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            s2cut = 0.5d0 - 0.5d0 * d2_v
            !--------
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-s2cc_b) * g_r(g_i_,i,2)
            enddo
            d1_w = (-s2cc_b) * (r(i,2) - s2cc_r)
            d2_v = exp(d1_w)
            d2_p = d2_v
            d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
            d1_p = 2.0d0 * (1.d0 - d2_v)
            d6_b = (-(s2cc_d * d1_p)) * d2_p
            do g_i_ = 1, g_p_
                g_s2cc(g_i_) = d6_b * g_d1_w(g_i_)
            enddo
            s2cc = s2cc_d * d4_v - s2cc_d
            !--------
            d3_v = s2 - s2cc
            d2_b = 1.0d0 + (-s2cut)
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d3_v * g_s2cut(g_i_) + s2cut * g_s2(g_i_) + d2_b * g_s2cc(g_i_)
            enddo
            s2 = s2cc + d3_v * s2cut
            !--------
            !
            d5_b = dehfc * 0.5d0
            do g_i_ = 1, g_p_
                g_dehfm(g_i_) = d5_b * g_cs(g_i_)
            enddo
            dehfm = dehf - (dehfcc + dehfc * 0.5d0 * (1.0d0 - cs))
            !--------
            d4_v = r(i,2) - rehf
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = (-behf) * g_r(g_i_,i,2) + (-d4_v) * g_behf( g_i_)
            enddo
            d1_w = (-behf) * d4_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_xehfm(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            xehfm = d2_v
            !--------
            !
            d3_v = dehfm * xehfm
            d4_v = xehfm - 2.0d0
            d6_b = d4_v * xehfm
            d5_b = d3_v + d4_v * dehfm
            do g_i_ = 1, g_p_
                g_s2m(g_i_) = d5_b * g_xehfm(g_i_) + d6_b * g_dehfm(g_i_)
            enddo
            s2m = d3_v * d4_v + uli2p
            !--------
            d3_b = 1.0d0 / 0.5d0
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,2)
            enddo
            d1_w = (r(i,2) - 2.2d0) / 0.5d0
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_s2cut(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            s2cut = 0.5d0 - 0.5d0 * d2_v
            !--------
            do g_i_ = 1, g_p_
                g_s2m(g_i_) = s2m * g_s2cut(g_i_) + s2cut * g_s2m(g_i_)
            enddo
            s2m = s2m * s2cut
            !--------
            !
            d3_b = 1.0d0 / s2alp
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,3)
            enddo
            d1_w = (r(i,3) - s2rho) / s2alp
            d4_v = (s2 - s2m) * 0.5d0
            d6_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d6_v *  d6_v)
            d7_v = 1.d0 + d6_v
            d7_b = d4_v * d1_p
            d8_b = d7_v * 0.5d0
            d2_b = 1.0d0 + (-d8_b)
            do g_i_ = 1, g_p_
                g_s2(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s2(g_i_) + d2_b * g_s2m(g_i_)
            enddo
            s2 = s2m + d4_v * d7_v
            !--------
            !
            !     The HF triplet
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
            !--------
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
            !--------
            !
            d2_v  = xehfa * xehfa
            d2_p  = 2.0d0 * xehfa
            d6_v  = (t2a1 * d2_v + t2a2 * xehfa) * 0.5d0
            d8_v  = 1.0d0 + cs
            d11_v = xehfa_180 * xehfa_180
            d1_p  = 2.0d0 * xehfa_180
            d15_v = (t2a3 * d11_v + t2a4 * xehfa_180) * 0.5d0
            d16_v = 1.0d0 - cs
            d7_b  = d16_v * 0.5d0
            d10_b = d7_b * t2a4 + d7_b * t2a3 * d1_p
            d6_b  = -d15_v + d6_v
            d14_b = d8_v * 0.5d0
            d17_b = d14_b * t2a2 + d14_b * t2a1 * d2_p
            do g_i_ = 1, g_p_
                g_t2(g_i_) = d10_b * g_xehfa_180(g_i_) + d6_b * g_cs(g_i_) + d17_b * g_xehfa(g_i_)
            enddo
            t2 = d6_v * d8_v + d15_v * d16_v
            !--------
            !
            !     The LiF singlet
            d5_v = (4.6498d0 + 7.0489d0 * r(i,3)) * (4.6498d0 + 7.0489d0 * r(i,3))
            d1_p = 2.0d0 * (4.6498d0 + 7.0489d0 * r(i,3))
            d6_v = r(i,3) * 103.57d0 / d5_v
            d7_b = (-d6_v) / d5_v * d1_p * 7.0489d0 + 1.0d0 / d5_v * 103.57d0
            do g_i_ = 1, g_p_
                g_belif(g_i_) = d7_b * g_r(g_i_,i,3)
            enddo
            belif = d6_v + 0.064076d0
            !--------
            !
            d3_b = 1.0d0 / gdel
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,2)
            enddo
            d1_w = (r(i,2) - grho) / gdel
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_gswt(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            gswt = 0.5d0 - 0.5d0 * d2_v
            !--------
            d2_v = gswt * 0.5d0
            d4_v = 1.0d0 + cs
            d5_b = d4_v * 0.5d0
            do g_i_ = 1, g_p_
                g_gswt(g_i_) = d2_v * g_cs(g_i_) + d5_b * g_gswt(g_i_)
            enddo
            gswt = d2_v * d4_v
            !--------
            d3_v = belifff - belif
            d2_b = 1.0d0 + (-gswt)
            do g_i_ = 1, g_p_
                g_beliff(g_i_) = d3_v * g_gswt(g_i_) + d2_b * g_belif(g_i_)
            enddo
            beliff = belif + gswt * d3_v
            !--------
            !
            if( nb .ne. 0 )then
                d5_v = (r(i,3) / gam) ** ( nb-1)
                d2_p = nb *  d5_v
                d5_v = d5_v * (r(i,3) / gam)
            else
                !            Maybe this should be  d5_v = (r(i, 3) / gam) **  nb
                d5_v = 1.0d0
                d2_p = 0.0d0
            endif
            d6_v = belif + d5_v
            if( nb .ne. 0 )then
                d9_v = (r(i,3) / gam) ** ( nb-1)
                d1_p = nb *  d9_v
                d9_v = d9_v * (r(i,3) / gam)
            else
                !            Maybe this should be  d9_v = (r(i, 3) / gam) **  nb
                d9_v = 1.0d0
                d1_p = 0.0d0
            endif
            d10_v = beliff + d9_v
            d11_v = beliff * d6_v / d10_v
            d2_b  = 1.0d0 / d10_v
            d3_b  = (-d11_v) / d10_v
            d4_b  = d3_b + d2_b * d6_v
            d8_b  = d2_b * beliff
            d7_b  = d3_b * d1_p * (1.0d0 / gam) + d8_b * d2_p * (1.0d0 / gam)
            do g_i_ = 1, g_p_
                g_belifmod(g_i_) = d7_b * g_r(g_i_,i,3) + d8_b * g_belif(g_i_) + d4_b * g_beliff(g_i_)
            enddo
            belifmod = d11_v
            !--------
            !
            d4_v = r(i,3) - relif
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
            !--------
            d2_v = (delif + uli2p) * xelif
            d3_v = xelif - 2.0d0
            d5_b = d2_v + d3_v * (delif + uli2p)
            do g_i_ = 1, g_p_
                g_s3a(g_i_) = d5_b * g_xelif(g_i_)
            enddo
            s3a = d2_v * d3_v + uli2p
            !--------
            d3_b = 1.0d0 / 0.5d0
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,3)
            enddo
            d1_w = (r(i,3) - 13.d0) / 0.5d0
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_s3cut(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            s3cut = 0.5d0 - 0.5d0 * d2_v
            !--------
            do g_i_ = 1, g_p_
                g_s3a(g_i_) = s3a * g_s3cut(g_i_) + s3cut * g_s3a(g_i_)
            enddo
            s3a = s3a * s3cut
            !--------
            !
            d3_b = 1.0d0 / swit_a
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_r(g_i_,i,2)
            enddo
            d1_w = (r(i,2) - swit_r) / swit_a
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d5_b = (1.d0 - swit_d) * 0.5d0 * d1_p
            do g_i_ = 1, g_p_
                g_swit(g_i_) = d5_b * g_d1_w(g_i_)
            enddo
            swit = swit_d + (1.d0 - swit_d) * 0.5d0 * (1.d0 + d2_v)
            !--------
            do g_i_ = 1, g_p_
                g_s3a(g_i_) = s3a * g_swit(g_i_) + swit * g_s3a(g_i_)
            enddo
            s3a = s3a * swit
            !--------
            !
            d5_b = delifc2 * 0.5d0 + (-(delifc1 * 0.5d0))
            do g_i_ = 1, g_p_
                g_delifm(g_i_) = d5_b * g_cs(g_i_)
            enddo
            delifm = delifc1 * 0.5d0 * (1.0d0 - cs) + delifc2 * 0.5d0 * (1.0d0 + cs)
            !--------
            !
            d5_b = belifc2 * 0.5d0 + (-(belifc1 * 0.5d0))
            do g_i_ = 1, g_p_
                g_belifm(g_i_) = d5_b * g_cs(g_i_)
            enddo
            belifm = belifc1 * 0.5d0 * (1.0d0 - cs) + belifc2 * 0.5d0 * (1.0d0 + cs)
            !--------
            !
            d4_v = r(i,3) - relif
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
            !--------
            d2_v = delifm + uli2p
            d4_v = d2_v * xelif
            d5_v = xelif - 2.0d0
            d4_b = d4_v + d5_v * d2_v
            d6_b = d5_v * xelif
            do g_i_ = 1, g_p_
                g_s3c(g_i_) = d4_b * g_xelif(g_i_) + d6_b * g_delifm(g_i_)
            enddo
            s3c = d4_v * d5_v
            !--------
            !
            d4_b = -((rho2_180 - rho2_0) * 0.5d0)
            do g_i_ = 1, g_p_
                g_rho2(g_i_) = d4_b * g_cs(g_i_)
            enddo
            rho2 = rho2_0 + (rho2_180 - rho2_0) * 0.5d0 * (1.0d0 - cs)
            !--------
            d4_b = -((theta2_180 - theta2_0) * 0.5d0)
            do g_i_ = 1, g_p_
                g_theta2(g_i_) = d4_b * g_cs(g_i_)
            enddo
            theta2 = theta2_0 + (theta2_180 - theta2_0) * 0.5d0 * (1.0d0 - cs)
            !--------
            d4_b = -((del2_180 - del2_0) * 0.5d0)
            do g_i_ = 1, g_p_
                g_del2(g_i_) = d4_b * g_cs(g_i_)
            enddo
            del2 = del2_0 + (del2_180 - del2_0) * 0.5d0 * (1.0d0 - cs)
            !--------
            !
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
            gc2 = d2_v * d4_v - d6_v * d9_v
            !--------
            d3_v = gc2 / del2
            d2_b = 1.0d0 / del2
            d3_b = (-d3_v) / del2
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d3_b * g_del2(g_i_) + d2_b * g_gc2(g_i_)
            enddo
            d1_w = d3_v
            d2_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            d4_b = (-0.5d0) * d1_p
            do g_i_ = 1, g_p_
                g_cut2(g_i_) = d4_b * g_d1_w(g_i_)
            enddo
            cut2 = 0.5d0 - 0.5d0 * d2_v
            !--------
            !
            d3_v = s3a - s3c
            d2_b = 1.0d0 + (-cut2)
            do g_i_ = 1, g_p_
                g_s3(g_i_) = d3_v * g_cut2(g_i_) + cut2 * g_s3a(g_i_) + d2_b * g_s3c(g_i_)
            enddo
            s3 = s3c + d3_v * cut2
            !--------
            !
            !     The LiF triplet
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
            !--------
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
            !--------
            d2_v  = xelif * xelif
            d2_p  = 2.0d0 * xelif
            d6_v  = (t3a1 * d2_v + t3a2 * xelif) * 0.5d0
            d8_v  = 1.0d0 + cs
            d11_v = xelif_180 * xelif_180
            d1_p  = 2.0d0 * xelif_180
            d15_v = (t3a3 * d11_v + t3a4 * xelif_180) * 0.5d0
            d16_v = 1.0d0 - cs
            d7_b  = d16_v * 0.5d0
            d10_b = d7_b * t3a4 + d7_b * t3a3 * d1_p
            d6_b  = -d15_v + d6_v
            d14_b = d8_v * 0.5d0
            d17_b = d14_b * t3a2 + d14_b * t3a1 * d2_p
            do g_i_ = 1, g_p_
                g_t3(g_i_) = d10_b * g_xelif_180(g_i_) + d6_b * g_cs(g_i_) + d17_b * g_xelif(g_i_)
            enddo
            t3 = d6_v * d8_v + d15_v * d16_v
            !--------
            !
            !     LEPS
            do g_i_ = 1, g_p_
                g_coul1(g_i_) = 0.5d0 * g_t1(g_i_) + 0.5d0 * g_s1(g_i_)
            enddo
            coul1 = 0.5d0 * (s1 + t1)
            !--------
            do g_i_ = 1, g_p_
                g_coul2(g_i_) = 0.5d0 * g_t2(g_i_) + 0.5d0 * g_s2(g_i_)
            enddo
            coul2 = 0.5d0 * (s2 + t2)
            !--------
            do g_i_ = 1, g_p_
                g_coul3(g_i_) = 0.5d0 * g_t3(g_i_) + 0.5d0 * g_s3(g_i_)
            enddo
            coul3 = 0.5d0 * (s3 + t3)
            !--------
            do g_i_ = 1, g_p_
                g_exch1(g_i_) = (-0.5d0) * g_t1(g_i_) + 0.5d0 * g_s1(g_i_)
            enddo
            exch1 = 0.5d0 * (s1 - t1)
            !--------
            do g_i_ = 1, g_p_
                g_exch2(g_i_) = (-0.5d0) * g_t2(g_i_) + 0.5d0 * g_s2(g_i_)
            enddo
            exch2 = 0.5d0 * (s2 - t2)
            !--------
            do g_i_ = 1, g_p_
                g_exch3(g_i_) = (-0.5d0) * g_t3(g_i_) + 0.5d0 * g_s3(g_i_)
            enddo
            exch3 = 0.5d0 * (s3 - t3)
            !--------
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
            !--------
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
            !--------
            d4_b = cplg2 + cplg2
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_cplg2(g_i_) + g_w(g_i_)
            enddo
            d1_w = w + cplg2 * cplg2
            d7_v = sqrt(d1_w)
            !          if ( d1_w .gt. 0.0d0 ) then
            d1_p = 1.0d0 / (2.0d0 *  d7_v)
            !          else
            !             call ehufDO (9,d1_w, d7_v, d1_p,
            !     +g_ehfid,
            !     +1165)
            !          endif
            d6_b = (-(1.0d0 / dsqrt(2.0d0))) * d1_p
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = d6_b * g_d1_w(g_i_) + g_coul3(g_i_) + g_coul2 (g_i_) + g_coul1(g_i_)
            enddo
            e(i) = coul1 + coul2 + coul3 - d7_v / dsqrt(2.0d0) + dehf
            !--------
            !
            ! LONG-RANGE FORCES * * * * * * * * * * * * * * * * * * * * * * * * * *
            ! HF DIPOLE
            d4_b = (-hfdip_a) * cauang
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_r(g_i_,i,2)
            enddo
            d1_w = (-hfdip_a) * (r(i,2) * cauang - hfdip_re)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_hfdip_y(g_i_) = (-d1_p) * g_d1_w(g_i_)
            enddo
            hfdip_y = 1.d0 - d2_v
            !--------
            do g_i_ = 1, g_p_
                g_hfdip(g_i_) = 0.0d0
            enddo
            hfdip = 0.d0
            !--------
            do j = 1, 11
                d3_v = dble(j - 1)
                !            if ( hfdip_y .ne. 0.0d0 ) then
                d4_v = hfdip_y ** ( d3_v - 2.0d0)
                d4_v = d4_v * hfdip_y
                d1_p = d3_v *  d4_v
                d4_v = d4_v * hfdip_y
                !            else
                !              (hfdip_y = 0)
                !               d4_v = hfdip_y **  d3_v
                !
                !               if (  d3_v .lt. 1.0d0 ) then
                !                  call ehbfDO (10,hfdip_y, d3_v, d4_v, d1_p, 0.0d0,
                !     +g_ehfid,
                !     +1209)
                !               else if (  d3_v .lt. 2.0d0 ) then
                !                  d1_p = 0.0d0
                !                  call ehbfDO (10,hfdip_y, d3_v, d4_v, d1_p, 0.0d0,
                !     +g_ehfid,
                !     +1214)
                !               else
                !                  d1_p = 0.0d0
                !               endif
                !            endif
                d6_b = 0.0d0
                d5_b = hfdip_m(j) * d1_p
                do g_i_ = 1, g_p_
                    g_hfdip(g_i_) = d5_b * g_hfdip_y(g_i_) + g_hfdip(g_i_)
                enddo
                hfdip = hfdip + hfdip_m(j) * d4_v
                !--------
            enddo
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = 20.d0 * g_r(g_i_,i,2)
            enddo
            d1_w = 20.d0 * (r(i,2) - 1.2d0)
            d2_v = hfdip * 0.5d0
            d4_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d4_v *  d4_v)
            d5_v = 1.d0 + d4_v
            d5_b = d2_v * d1_p
            d6_b = d5_v * 0.5d0
            do g_i_ = 1, g_p_
                g_hfdip(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_hfdip(g_i_)
            enddo
            hfdip = d2_v * d5_v
            !--------
            !
            ! HF POLARIZATION
            d2_v = coschihf * coschihf
            d2_p = 2.0d0 * coschihf
            d5_v = coschihf * coschihf
            d1_p = 2.0d0 * coschihf
            d5_b = hfpolpar * d1_p + (-hfpolperp) * d2_p
            do g_i_ = 1, g_p_
                g_hfpol(g_i_) = d5_b * g_coschihf(g_i_)
            enddo
            hfpol = hfpolperp * (1.d0 - d2_v) + hfpolpar * d5_v
            !--------
            hfpolavg = 2.d0 * hfpolperp / 3.d0 + hfpolpar / 3.d0
            !
            ! HF QUADRUPOLE
            d4_b = (-hfqua_a) * cauang
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_r(g_i_,i,2)
            enddo
            d1_w = (-hfqua_a) * (r(i,2) * cauang - hfqua_re)
            d2_v = exp(d1_w)
            d1_p = d2_v
            do g_i_ = 1, g_p_
                g_hfqua_y(g_i_) = (-d1_p) * g_d1_w(g_i_)
            enddo
            hfqua_y = 1.d0 - d2_v
            !--------
            do g_i_ = 1, g_p_
                g_hfqua(g_i_) = 0.0d0
            enddo
            hfqua = 0.d0
            !--------
            do j = 1, 10
                d3_v = dble(j - 1)
                !            if ( hfqua_y .ne. 0.0d0 ) then
                d4_v = hfqua_y ** ( d3_v - 2.0d0)
                d4_v = d4_v * hfqua_y
                d1_p = d3_v *  d4_v
                d4_v = d4_v * hfqua_y
                !            else
                !              (hfqua_y = 0)
                !               d4_v = hfqua_y **  d3_v
                !
                !               if (  d3_v .lt. 1.0d0 ) then
                !                  call ehbfDO (10,hfqua_y, d3_v, d4_v, d1_p, 0.0d0,
                !     +g_ehfid,
                !     +1289)
                !               else if (  d3_v .lt. 2.0d0 ) then
                !                  d1_p = 0.0d0
                !                  call ehbfDO (10,hfqua_y, d3_v, d4_v, d1_p, 0.0d0,
                !     +g_ehfid,
                !     +1294)
                !               else
                !                  d1_p = 0.0d0
                !               endif
                !            endif
                d6_b = 0.0d0
                d5_b = hfqua_m(j) * d1_p
                do g_i_ = 1, g_p_
                    g_hfqua(g_i_) = d5_b * g_hfqua_y(g_i_) + g_hfqua(g_i_)
                enddo
                hfqua = hfqua + hfqua_m(j) * d4_v
                !--------
            enddo
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = 20.d0 * g_r(g_i_,i,2)
            enddo
            d1_w = 20.d0 * (r(i,2) - 1.2d0)
            d2_v = hfqua * 0.5d0
            d4_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d4_v *  d4_v)
            d5_v = 1.d0 + d4_v
            d5_b = d2_v * d1_p
            d6_b = d5_v * 0.5d0
            do g_i_ = 1, g_p_
                g_hfqua(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_hfqua(g_i_)
            enddo
            hfqua = d2_v * d5_v
            !--------
            !
            ! LiF DIPOLE
            do g_i_ = 1, g_p_
                g_lifdip_y(g_i_) = g_r(g_i_,i,3)
            enddo
            lifdip_y = r(i,3) - lifdip_4
            !--------
            d2_v = lifdip_y ** ( 8 - 2)
            d2_v = d2_v * lifdip_y
            d1_p = 8 *  d2_v
            d2_v = d2_v * lifdip_y
            d4_b = lifdip_3 + lifdip_2 * d1_p
            do g_i_ = 1, g_p_
                g_lifdip_beta(g_i_) = d4_b * g_lifdip_y(g_i_)
            enddo
            lifdip_beta = lifdip_1 + lifdip_2 * d2_v + lifdip_3 * lifdip_y
            !--------
            d5_v = (r(i,3) - lifdip_re) * (r(i,3) - lifdip_re)
            d1_p = 2.0d0 * (r(i,3) - lifdip_re)
            d5_b = (-lifdip_beta) * d1_p
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d5_b * g_r(g_i_,i,3) + (-d5_v) * g_lifdip_beta(g_i_)
            enddo
            d1_w = (-lifdip_beta) * d5_v
            d2_v = exp(d1_w)
            d1_p = d2_v
            d3_b = lifdip_d * d1_p
            do g_i_ = 1, g_p_
                g_lifdip(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            lifdip = lifdip_d * d2_v
            !--------
            !
            ! Li(2p) POLARIZATION
            d2_v = coschihf * coschihf
            d2_p = 2.0d0 * coschihf
            d5_v = coschihf * coschihf
            d1_p = 2.0d0 * coschihf
            d5_b = li2ppolpar * d1_p + (-li2ppolperp) * d2_p
            do g_i_ = 1, g_p_
                g_li2ppol(g_i_) = d5_b * g_coschihf(g_i_)
            enddo
            li2ppol = li2ppolperp * (1.d0 - d2_v) + li2ppolpar * d5_v
            !--------
            li2ppolavg = 2.d0 * li2ppolperp / 3.d0 + li2ppolpar / 3.d0
            !
            ! HF + Li(2p)
            !     DIPOLE-INDUCED-DIPOLE
            d2_v  = coschihf * coschihf
            d2_p  = 2.0d0 * coschihf
            d5_v  = (-0.5d0) * (3.d0 * d2_v + 1.d0)
            d7_v  = d5_v * li2ppol
            d9_v  = hfdip * hfdip
            d1_p  = 2.0d0 * hfdip
            d4_b  = d7_v * d1_p
            d6_b  = d9_v * d5_v
            d10_b = d9_v * li2ppol * (-0.5d0) * 3.d0 * d2_p
            do g_i_ = 1, g_p_
                g_hf_did(g_i_) = d4_b * g_hfdip(g_i_) + d6_b * g_li2ppol(g_i_) + d10_b * g_coschihf(g_i_)
            enddo
            hf_did = d7_v * d9_v
            !--------
            d3_v = bigrhf ** ( 7 - 2)
            d3_v = d3_v * bigrhf
            d1_p = 7 *  d3_v
            d3_v = d3_v * bigrhf
            d4_v = d3_v + rhf0did ** 7
            d5_v = hf_did / d4_v
            d6_v = d5_v / cevau
            d4_b = bigrhf * (1.0d0 / cevau)
            d5_b = d4_b * (1.0d0 / d4_v)
            d3_b = d6_v + d4_b * ((-d5_v) / d4_v) * d1_p
            do g_i_ = 1, g_p_
                g_hf_did(g_i_) = d3_b * g_bigrhf(g_i_) + d5_b * g_hf_did(g_i_)
            enddo
            hf_did = d6_v * bigrhf
            !--------
            !     cut off angular dependence in interaction region
            d3_v = (bigrhf / rhf0did) ** ( 4 - 2)
            d3_v = d3_v * (bigrhf / rhf0did)
            d1_p = 4 *  d3_v
            d3_v = d3_v * (bigrhf / rhf0did)
            d4_b = (-d1_p) * (1.0d0 / rhf0did)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_bigrhf(g_i_)
            enddo
            d1_w  = -d3_v
            d3_v  = coschihf * coschihf
            d2_p  = 2.0d0 * coschihf
            d5_v  = 3.d0 * d3_v + 1.d0
            d6_v  = 2.5d0 * li2ppolavg / d5_v
            d8_v  = d6_v / li2ppol
            d9_v  = d8_v - 1.d0
            d11_v = exp(d1_w)
            d1_p  = d11_v
            d13_v = 1.d0 + d9_v * d11_v
            d7_b  = hf_did * d9_v * d1_p
            d8_b  = hf_did * d11_v
            d10_b = d8_b * ((-d8_v) / li2ppol)
            d14_b = d8_b * (1.0d0 / li2ppol) * ((-d6_v) / d5_v) * 3.d0 * d2_p
            do g_i_ = 1, g_p_
                g_hf_did(g_i_) = d7_b * g_d1_w(g_i_) + d10_b * g_li2ppol(g_i_) + d14_b * g_coschihf(g_i_) + d13_v * &
                    & g_hf_did(g_i_)
            enddo
            hf_did = hf_did * d13_v
            !--------
            !
            !     DISPERSION
            d2_v = (-1.5d0) * li2pie * hfie * li2ppol
            d2_b = 1.0d0 / (li2pie + hfie)
            d4_b = d2_b * d2_v
            d5_b = d2_b * hfpol * ((-1.5d0) * li2pie * hfie)
            do g_i_ = 1, g_p_
                g_hf_disp(g_i_) = d4_b * g_hfpol(g_i_) + d5_b * g_li2ppol(g_i_)
            enddo
            hf_disp = d2_v * hfpol / (li2pie + hfie)
            !--------
            d3_v = bigrhf ** ( 7 - 2)
            d3_v = d3_v * bigrhf
            d1_p = 7 *  d3_v
            d3_v = d3_v * bigrhf
            d4_v = d3_v + rhf0disp ** 7
            d5_v = hf_disp / d4_v
            d4_b = bigrhf * (1.0d0 / d4_v)
            d3_b = d5_v + bigrhf * ((-d5_v) / d4_v) * d1_p
            do g_i_ = 1, g_p_
                g_hf_disp(g_i_) = d3_b * g_bigrhf(g_i_) + d4_b * g_hf_disp(g_i_)
            enddo
            hf_disp = d5_v * bigrhf
            !--------
            !     cut off angular dependence in interaction region
            d3_v = (bigrhf / rhf0disp) ** ( 4 - 2)
            d3_v = d3_v * (bigrhf / rhf0disp)
            d1_p = 4 *  d3_v
            d3_v = d3_v * (bigrhf / rhf0disp)
            d4_b = (-d1_p) * (1.0d0 / rhf0disp)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_bigrhf(g_i_)
            enddo
            d1_w  = -d3_v
            d3_v  = li2ppolavg * hfpolavg / hfpol
            d5_v  = d3_v / li2ppol
            d6_v  = d5_v - 1.d0
            d8_v  = exp(d1_w)
            d1_p  = d8_v
            d10_v = 1.d0 + d6_v * d8_v
            d7_b  = hf_disp * d6_v * d1_p
            d8_b  = hf_disp * d8_v
            d10_b = d8_b * ((-d5_v) / li2ppol)
            d11_b = d8_b * (1.0d0 / li2ppol) * ((-d3_v) / hfpol)
            do g_i_ = 1, g_p_
                g_hf_disp(g_i_) = d7_b * g_d1_w(g_i_) + d10_b * g_li2ppol(g_i_) + d11_b * g_hfpol(g_i_) + d10_v * &
                    & g_hf_disp(g_i_)
            enddo
            hf_disp = hf_disp * d10_v
            !--------
            !
            !     QUADRUPOLE-QUADRUPOLE
            d3_v = 3.d0 / 4.d0 * hfqua * li2pqua
            d5_v = coschihf * coschihf
            d1_p = 2.0d0 * coschihf
            d7_v = 3.d0 - 7.d0 * d5_v
            d6_b = (-d3_v) * 7.d0 * d1_p
            d8_b = d7_v * li2pqua * (3.d0 / 4.d0)
            do g_i_ = 1, g_p_
                g_hf_qq(g_i_) = d6_b * g_coschihf(g_i_) + d8_b * g_hfqua(g_i_)
            enddo
            hf_qq = d3_v * d7_v
            !--------
            d3_v = bigrhf ** ( 6 - 2)
            d3_v = d3_v * bigrhf
            d1_p = 6 *  d3_v
            d3_v = d3_v * bigrhf
            d4_v = d3_v + rhf0qq ** 6
            d5_v = hf_qq / d4_v
            d6_v = d5_v / cevau
            d4_b = bigrhf * (1.0d0 / cevau)
            d5_b = d4_b * (1.0d0 / d4_v)
            d3_b = d6_v + d4_b * ((-d5_v) / d4_v) * d1_p
            do g_i_ = 1, g_p_
                g_hf_qq(g_i_) = d3_b * g_bigrhf(g_i_) + d5_b * g_hf_qq(g_i_)
            enddo
            hf_qq = d6_v * bigrhf
            !--------
            !     cut off angular dependence in interaction region
            d3_v = (bigrhf / rhf0qq) ** ( 4 - 2)
            d3_v = d3_v * (bigrhf / rhf0qq)
            d1_p = 4 *  d3_v
            d3_v = d3_v * (bigrhf / rhf0qq)
            d4_b = (-d1_p) * (1.0d0 / rhf0qq)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_bigrhf(g_i_)
            enddo
            d1_w  = -d3_v
            d3_v  = coschihf * coschihf
            d2_p  = 2.0d0 * coschihf
            d5_v  = 3.d0 - 7.d0 * d3_v
            d6_v  = (-.5d0) / d5_v
            d7_v  = d6_v - 1.d0
            d9_v  = exp(d1_w)
            d1_p  = d9_v
            d11_v = 1.d0 + d7_v * d9_v
            d7_b  = hf_qq * d7_v * d1_p
            d12_b = (-(hf_qq * d9_v * ((-d6_v) / d5_v))) * 7.d0 * d2_p
            do g_i_ = 1, g_p_
                g_hf_qq(g_i_) = d7_b * g_d1_w(g_i_) + d12_b * g_coschihf(g_i_) + d11_v * g_hf_qq(g_i_)
            enddo
            hf_qq = hf_qq * d11_v
            !--------
            !
            !     DIPOLE-QUADRUPOLE
            d3_v = 3.d0 / 4.d0 * hfdip * li2pqua
            d5_b = coschihf * li2pqua * (3.d0 / 4.d0)
            do g_i_ = 1, g_p_
                g_hf_dq(g_i_) = d3_v * g_coschihf(g_i_) + d5_b * g_hfdip(g_i_)
            enddo
            hf_dq = d3_v * coschihf
            !--------
            d3_v = bigrhf ** ( 5 - 2)
            d3_v = d3_v * bigrhf
            d1_p = 5 *  d3_v
            d3_v = d3_v * bigrhf
            d4_v = d3_v + rhf0dq ** 5
            d5_v = hf_dq / d4_v
            d6_v = d5_v / cevau
            d4_b = bigrhf * (1.0d0 / cevau)
            d5_b = d4_b * (1.0d0 / d4_v)
            d3_b = d6_v + d4_b * ((-d5_v) / d4_v) * d1_p
            do g_i_ = 1, g_p_
                g_hf_dq(g_i_) = d3_b * g_bigrhf(g_i_) + d5_b * g_hf_dq(g_i_)
            enddo
            hf_dq = d6_v * bigrhf
            !--------
            !     cut off angular dependence in interaction region
            d3_v = (bigrhf / rhf0dq) ** ( 4 - 2)
            d3_v = d3_v * (bigrhf / rhf0dq)
            d1_p = 4 *  d3_v
            d3_v = d3_v * (bigrhf / rhf0dq)
            d4_b = (-d1_p) * (1.0d0 / rhf0dq)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_bigrhf(g_i_)
            enddo
            d1_w = -d3_v
            d3_v = exp(d1_w)
            d1_p = d3_v
            d4_b = hf_qq * d1_p
            do g_i_ = 1, g_p_
                g_hf_qq(g_i_) = d4_b * g_d1_w(g_i_) + d3_v * g_hf_qq(g_i_)
            enddo
            hf_qq = hf_qq * d3_v
            !--------
            !
            !     SUM FORCES
            do g_i_ = 1, g_p_
                g_e_lr_hf(g_i_) = g_hf_did(g_i_) + g_hf_disp(g_i_) + g_hf_qq (g_i_) + g_hf_dq(g_i_)
            enddo
            e_lr_hf = hf_dq + hf_qq + hf_disp + hf_did
            !--------
            !     cut off at large HF distances
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = r2lrcut_a * g_r(g_i_,i,2)
            enddo
            d1_w = r2lrcut_a * (r(i,2) - r2lrcut_r)
            d2_v = e_lr_hf * 0.5d0
            d4_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d4_v *  d4_v)
            d5_v = 1.d0 - d4_v
            d5_b = (-d2_v) * d1_p
            d6_b = d5_v * 0.5d0
            do g_i_ = 1, g_p_
                g_e_lr_hf(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_e_lr_hf(g_i_)
            enddo
            e_lr_hf = d2_v * d5_v
            !--------
            !
            ! LiF + H
            !     DIPOLE-INDUCED-DIPOLE
            d2_v = coschilif * coschilif
            d2_p = 2.0d0 * coschilif
            d6_v = (-0.5d0) * (3.d0 * d2_v + 1.d0) * hpol
            d8_v = lifdip * lifdip
            d1_p = 2.0d0 * lifdip
            d4_b = d6_v * d1_p
            d9_b = d8_v * hpol * (-0.5d0) * 3.d0 * d2_p
            do g_i_ = 1, g_p_
                g_lif_did(g_i_) = d4_b * g_lifdip(g_i_) + d9_b * g_coschilif (g_i_)
            enddo
            lif_did = d6_v * d8_v
            !--------
            d3_v = bigrlif ** ( 6 - 2)
            d3_v = d3_v * bigrlif
            d1_p = 6 *  d3_v
            d3_v = d3_v * bigrlif
            d4_v = d3_v + rlif0did ** 6
            d5_v = lif_did / d4_v
            d2_b = 1.0d0 / cevau
            d3_b = d2_b * (1.0d0 / d4_v)
            d6_b = d2_b * ((-d5_v) / d4_v) * d1_p
            do g_i_ = 1, g_p_
                g_lif_did(g_i_) = d6_b * g_bigrlif(g_i_) + d3_b * g_lif_did( g_i_)
            enddo
            lif_did = d5_v / cevau
            !--------
            !     cut off angular dependence in interaction region
            d3_v = (bigrlif / rlif0did) ** ( 4 - 2)
            d3_v = d3_v * (bigrlif / rlif0did)
            d1_p = 4 *  d3_v
            d3_v = d3_v * (bigrlif / rlif0did)
            d4_b = (-d1_p) * (1.0d0 / rlif0did)
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_bigrlif(g_i_)
            enddo
            d1_w  = -d3_v
            d3_v  = coschilif * coschilif
            d2_p  = 2.0d0 * coschilif
            d5_v  = 3.d0 * d3_v + 1.d0
            d6_v  = 2.5d0 / d5_v
            d7_v  = d6_v - 1.d0
            d9_v  = exp(d1_w)
            d1_p  = d9_v
            d11_v = 1.d0 + d7_v * d9_v
            d7_b  = lif_did * d7_v * d1_p
            d12_b = lif_did * d9_v * ((-d6_v) / d5_v) * 3.d0 * d2_p
            do g_i_ = 1, g_p_
                g_lif_did(g_i_) = d7_b * g_d1_w(g_i_) + d12_b * g_coschilif( g_i_) + d11_v * g_lif_did(g_i_)
            enddo
            lif_did = lif_did * d11_v
            !--------
            !
            !     DISPERSION
            do g_i_ = 1, g_p_
                g_lif_disp(g_i_) = 0.0d0
            enddo
            lif_disp = (-1.5d0) * lifie * hie * lifpol * hpol / (lifie + hie)
            !--------
            d3_v = bigrlif ** ( 6 - 2)
            d3_v = d3_v * bigrlif
            d1_p = 6 *  d3_v
            d3_v = d3_v * bigrlif
            d4_v = d3_v + rlif0disp ** 6
            d5_v = lif_disp / d4_v
            d2_b = 1.0d0 / d4_v
            d5_b = (-d5_v) / d4_v * d1_p
            do g_i_ = 1, g_p_
                g_lif_disp(g_i_) = d5_b * g_bigrlif(g_i_) + d2_b * g_lif_disp(g_i_)
            enddo
            lif_disp = d5_v
            !--------
            !     no angular dependence
            !
            !     SUM FORCES
            do g_i_ = 1, g_p_
                g_e_lr_lif(g_i_) = g_lif_did(g_i_) + g_lif_disp(g_i_)
            enddo
            e_lr_lif = lif_disp + lif_did
            !--------
            !     cut off at large LiF distances
            do g_i_ = 1, g_p_
                g_d1_w(g_i_) = r3lrcut_a * g_r(g_i_,i,3)
            enddo
            d1_w = r3lrcut_a * (r(i,3) - r3lrcut_r)
            d2_v = e_lr_lif * 0.5d0
            d4_v = tanh (d1_w)
            d1_p = 1.0d0 - ( d4_v *  d4_v)
            d5_v = 1.d0 - d4_v
            d5_b = (-d2_v) * d1_p
            d6_b = d5_v * 0.5d0
            do g_i_ = 1, g_p_
                g_e_lr_lif(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_e_lr_lif(g_i_)
            enddo
            e_lr_lif = d2_v * d5_v
            !--------
            !
            ! LiH + F
            ! no long-range forces included
            !
            ! COMBINE INTERACTION REGION ENERGY AND LONG-RANGE FORCES
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = g_e_lr_hf(g_i_) + g_e_lr_lif(g_i_) + g_e(g_i_,i)
            enddo
            e(i) = e(i) + e_lr_lif + e_lr_hf
            !--------
            do g_i_ = 1, g_p_
                g_e(g_i_,i) = cevau * g_e(g_i_,i)
            enddo
            e(i) = e(i) * cevau
            !--------
            !
        enddo
    end subroutine dispatch
end module chempotpy_lifh_leps_jhtp2_2002_pot22_legacy_state_fb8669b7

!> @brief Compatibility interface for pot22.
subroutine pot22( g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt )
    use chempotpy_lifh_leps_jhtp2_2002_pot22_legacy_state_fb8669b7, only: dispatch
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
end subroutine pot22

!
