!   System:                        LiFH
!   Functional form:               Modified Extended LEPS (2x2 diabatic fit)
!   Common name:                   LiFH surface fit JS
!   Number of derivatives:         0
!   Number of electronic surfaces: 2
!   Interface:                     3-1V
!   References:      A. W. Jasper, M. D. Hack, D. G. Truhlar, and P. Pie
!                    J. Chem. Phys., Vol. 116, 8353 (2002).
!   Notes:    This fit contains cut-off long-range forces.
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
!                   LiFHJS POTENTIAL ENERGY MATRIX
!  Created by A. Jasper on August 22, 2001.
!  Based on the LiFHJ potential energy matrix, which was
!  based on the LiFHH potential energy matrix created by M. Hack,
!  which was based on the ab initio calculations of P. Piecuch.
!  This surface has the long-range forces cut off.
!  Reference: A. W. Jasper, M. D. Hack, D. G. Truhlar, and P. Piecuch,
!  J. Chem. Phys., in preparation.

!> @brief Implement LiFH_LEPS_JHTP1_2002, a diabatic potential matrix for LiFH.
module chempotpy_lifh_leps_jhtp1_2002_dpem_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the LiFH_LEPS_JHTP1_2002 diabatic potential matrix for LiFH.
    subroutine evaluate( x, igrad, u, ug )
        use chempotpy_legacy_units,       only: legacy_hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: legacy_triatomic_distances_to_bohr
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 2
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: u(nstates,nstates)
        double precision, intent(out) :: ug(nstates,nstates,natoms,3)
        double precision              :: nt, r(1,3), v(1)
        double precision              :: u11(1), u12(1), u22(1)
        integer                       :: i, j, k, l
        !initialize
        u  = 0.d0; ug = 0.d0
        nt = 1
        ! input cartesian is ClHH
        call legacy_triatomic_distances_to_bohr(x, r(1,:))
        if( igrad == 0 )then
            call pot(r, u11, 1, 1)
            call pot(r, u12, 1, 2)
            call pot(r, u22, 1, 3)
            u(1,1) = u11(1)*legacy_hartree_to_electronvolt
            u(1,2) = u12(1)*legacy_hartree_to_electronvolt
            u(2,1) = u12(1)*legacy_hartree_to_electronvolt
            u(2,2) = u22(1)*legacy_hartree_to_electronvolt
        else
            write(*, *) "Only DPEM is available"
        endif
    end subroutine evaluate
end module chempotpy_lifh_leps_jhtp1_2002_dpem_pes

!> @brief Expose LiFH_LEPS_JHTP1_2002 through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_lifh_leps_jhtp1_2002_dpem_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 2
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: u(nstates,nstates)
    double precision, intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem

!> @brief Evaluate EvdRdX for the LiFH_LEPS_JHTP1_2002_DPEM potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    implicit none
    integer, parameter            :: first_atoms(3)  = [1, 2, 1]
    integer, parameter            :: second_atoms(3) = [2, 3, 3]
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!> @brief Evaluate diagonalize for the LiFH_LEPS_JHTP1_2002_DPEM potential.
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

!> @brief Evaluate pot for the LiFH_LEPS_JHTP1_2002_DPEM potential.
subroutine pot( r, e, nt, nsurf )
    implicit none
    double precision :: r(nt,3), e(nt)
    integer          :: i, nt, nsurf
    if( nsurf .eq. 1 ) call pot11(r, e, nt)
    if( nsurf .eq. 2 ) call pot12(r, e, nt)
    if( nsurf .eq. 3 ) call pot22(r, e, nt)
    if( nsurf .gt. 3.or.nsurf .lt. 1 )then
        write(6, *)' WARNING:  nsurf = ', nsurf, &
            & ' in LiFHJS potential routine!'
    endif
end subroutine pot

!                            U11 SURFACE

!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_lifh_leps_jhtp1_2002_dpem_pot11_legacy_state_a97514e1
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        implicit none
        optional :: r, e, nt
        ! MISC VARIABLES
        integer          :: i, nt, j
        double precision :: cevau
        double precision :: r(:,:), e(:)
        double precision :: r1s, r2s, r3s
        ! INTERACTION REGION VARIABLES
        !     LiH
        double precision :: s1, relih, bes1, s1c2, s1c3, s1a, xelih, bes1b, s1c, &
            & rcut2del, rcut2rho, rcut2, r1r3, r1r3al
        !     HF
        double precision :: s2, behf, yhf, rehf, dehf, xehf, s2m, behfm, dehfm, &
            & rehfm, xehfm, s2rho, s2alp, dehfc, dehfcc, yhf_r, behf_0, behf_1, &
            & behf_2
        !     LiF
        double precision :: s3, res3, bes3, des3, res3c, bes3c, s3c, s3a, xes3, des3c, &
            & rcutdel, rcutrho, rcut, r1r3al2, bes3i, bes3z, bes3r
        ! LONG RANGE FORCES VARIABLES
        !     HF physical properties
        double precision :: hfdip_y, hfdip_a, hfdip_re, hfdip_m(11), hfdip, &
            & hfpol, hfpolperp, hfpolpar, hfie, hfpolavg
        !     Li physical properties
        double precision :: li2spol, rhf0did, rhf0disp, li2sie
        !     misc long-range variables
        double precision :: hf_did, hf_disp, e_lr_hf
        double precision :: mass(3), cauang, caudeb, costheta, temp, bigr, &
            & coschi, r2lrcut_a, r2lrcut_r
        ! MISC DATA
        ! INTERACTION REGION DATA
        !     LiH
        !     HF
        !     LiF
        ! LONG RANGE FORCES DATA
        !     HF dipole moment (angstroms vs au)
        !     HF physical properties (polarization in au, ie in eV)
        !     Li(2s) physical properties (polarization in au, ie in eV)
        !     misc long-range parameters (bohr)
        behf      = 0.d0
        behf_0    = 0.d0
        behf_1    = 0.d0
        behf_2    = 0.d0
        behfm     = 0.d0
        bes1      = 0.d0
        bes1b     = 0.d0
        bes3      = 0.d0
        bes3c     = 0.d0
        bes3i     = 0.d0
        bes3r     = 0.d0
        bes3z     = 0.d0
        bigr      = 0.d0
        cauang    = 0.d0
        caudeb    = 0.d0
        cevau     = 0.d0
        coschi    = 0.d0
        costheta  = 0.d0
        dehf      = 0.d0
        dehfc     = 0.d0
        dehfcc    = 0.d0
        dehfm     = 0.d0
        des3      = 0.d0
        des3c     = 0.d0
        e_lr_hf   = 0.d0
        hf_did    = 0.d0
        hf_disp   = 0.d0
        hfdip     = 0.d0
        hfdip_a   = 0.d0
        hfdip_m   = 0.d0
        hfdip_re  = 0.d0
        hfdip_y   = 0.d0
        hfie      = 0.d0
        hfpol     = 0.d0
        hfpolavg  = 0.d0
        hfpolpar  = 0.d0
        hfpolperp = 0.d0
        li2sie    = 0.d0
        li2spol   = 0.d0
        mass      = 0.d0
        r1r3      = 0.d0
        r1r3al    = 0.d0
        r1r3al2   = 0.d0
        r1s       = 0.d0
        r2lrcut_a = 0.d0
        r2lrcut_r = 0.d0
        r2s       = 0.d0
        r3s       = 0.d0
        rcut      = 0.d0
        rcut2     = 0.d0
        rcut2del  = 0.d0
        rcut2rho  = 0.d0
        rcutdel   = 0.d0
        rcutrho   = 0.d0
        rehf      = 0.d0
        rehfm     = 0.d0
        relih     = 0.d0
        res3      = 0.d0
        res3c     = 0.d0
        rhf0did   = 0.d0
        rhf0disp  = 0.d0
        s1        = 0.d0
        s1a       = 0.d0
        s1c       = 0.d0
        s1c2      = 0.d0
        s1c3      = 0.d0
        s2        = 0.d0
        s2alp     = 0.d0
        s2m       = 0.d0
        s2rho     = 0.d0
        s3        = 0.d0
        s3a       = 0.d0
        s3c       = 0.d0
        temp      = 0.d0
        xehf      = 0.d0
        xehfm     = 0.d0
        xelih     = 0.d0
        xes3      = 0.d0
        yhf       = 0.d0
        yhf_r     = 0.d0
        cevau     = 0.036749309d0
        relih     = 1.3432062d0
        bes1      = 1.8058651d0
        s1c2      = 12.937438d0
        s1c3      = 11.0913d0
        bes1b     = 1.2111436d0
        rcut2rho  = 1.07143d0
        rcut2del  = 0.6d0
        r1r3al    = 0.4d0
        rehf      = 1.733d0
        dehf      = 6.122d0
        rehfm     = 1.6739d0
        dehfc     = 0.344574d0
        dehfcc    = 0.8660801d0
        s2rho     = 1.63d0
        s2alp     = 2.9941d0
        behfm     = 0.74633431d0
        yhf_r     = 2.1042d0
        behf_0    = 1.1622d0
        behf_1    = 0.025647d0
        behf_2    = 0.059062d0
        res3      = 3.6d0
        bes3      = 1.49061584d0
        des3      = 0.17426687d0
        res3c     = 2.48571d0
        bes3c     = 1.92857d0
        des3c     = 0.25d0
        rcutrho   = 0.4d0
        rcutdel   = 0.89333d0
        r1r3al2   = 0.7333d0
        bes3z     = 0.d0
        bes3r     = 7.4799608d0
        hfdip_a   = 1.13d0
        hfdip_re  = 0.924212d0
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
        do i=1, nt
            ! PRECOMPUTE THINGS FOR USE LATER...
            !     cosine of HF bond angle (i.e., Li-F-H angle)
            costheta = (r(i,3)**2-r(i,2)**2-r(i,1)**2)/(-2.d0*r(i,1)*r(i,2))
            temp     = mass(3)/(mass(2)+mass(3))*r(i,2)
            !     Li-FH Jacobi distance
            bigr = dsqrt(r(i,1)**2+temp**2-2.d0*r(i,1)*temp*costheta)
            !     cosine of HF Jacobi angle (i.e., Li--HF angle)
            coschi = -(r(i,1)**2-bigr**2-temp**2)/(-2.d0*bigr*temp)
            coschi = max(coschi, -1.d0)
            coschi = min(coschi, 1.d0)
            !     squared!
            r1s = r(i,1)**2
            r2s = r(i,2)**2
            r3s = r(i,3)**2
            ! INTERACTION REGION
            !     LiH
            xelih = dexp(-bes1*(r(i,1)-relih))
            s1a   = s1c2*xelih
            s1c   = s1c3*dexp(-bes1b*(r(i,1)-relih))
            r1r3  = r(i,1) - r(i,3) + r1r3al*r(i,2)
            rcut2 = 0.5d0*(1.0d0 - dtanh((r1r3-rcut2rho)/rcut2del))
            s1    = s1a  + (s1c-s1a)*rcut2
            !     HF
            yhf   = r(i,2) - yhf_r
            behf  = behf_0-behf_1*yhf+behf_2*yhf**2
            xehf  = dexp(-behf*(r(i,2)-rehf))
            s2    = dehf*xehf*(xehf-2.0d0)
            dehfm = dehf - (dehfcc+dehfc*0.5d0*(1.d0-coschi))
            xehfm = dexp(-behfm*(r(i,2)-rehfm))
            s2m   = dehfm*xehfm*(xehfm-2.0d0)
            s2    = s2m + (s2 - s2m)*0.5d0*(1.d0+dtanh((r(i,3)-s2rho)/s2alp))
            !     LiF
            xes3 = dexp(-bes3*(r(i,3)-res3))
            s3   = des3*xes3*(xes3-2.0d0)
            !     Sum of diatomics
            e(i) = s1 + s2 + s3 + dehf
            ! EVALUATE THE LONG-RANGE FORCES
            ! HF DIPOLE MOMENT
            hfdip_y = 1.d0 - dexp(-hfdip_a*(r(i,2)*cauang-hfdip_re))
            hfdip   = 0.d0
            do j=1, 11
                hfdip = hfdip + hfdip_m(j)*hfdip_y**dble(j-1)
            enddo
            hfdip = hfdip*0.5d0*(1.d0 + dtanh(20.d0*(r(i,2)-1.2d0)))
            ! HF POLARIZATION
            hfpol    = hfpolperp*(1.d0 - coschi**2) + hfpolpar*coschi**2
            hfpolavg = 2.d0*hfpolperp/3.d0 + hfpolpar/3.d0
            ! HF + Li(2s) arrangement
            !     DIPOLE-INDUCED-DIPOLE
            hf_did = -0.5d0*(3.d0*coschi**2 + 1.d0)*li2spol*hfdip**2
            hf_did = hf_did/(bigr**7 + rhf0did**7)/cevau*bigr
            !     cut off angular dependence in interaction region
            hf_did = hf_did*(1.d0 + (2.5d0/(3.d0*coschi**2 + 1.d0) - 1.d0) *dexp(-(bigr/rhf0did)**4))
            !     DISPERSION
            hf_disp = -1.5d0*li2sie*hfie*li2spol*hfpol/(li2sie+hfie)
            hf_disp = hf_disp/(bigr**7 + rhf0disp**7)*bigr
            !     cut off angular dependence in interaction region
            hf_disp = hf_disp*(1.d0 + (hfpolavg/hfpol - 1.d0) *dexp(-(bigr/rhf0disp)**4))
            !     SUM FORCES
            e_lr_hf = hf_disp + hf_did
            !     cut off LR forces for large r[hf] i.e., r(i,2)
            e_lr_hf = e_lr_hf*0.5d0*(1.d0-dtanh(r2lrcut_a*(r(i,2)-r2lrcut_r)))
            ! LiF + H
            !     NO LONG-RANGE FORCES
            ! LiH + F
            !     NO LONG-RANGE FORCES
            ! CUT OFF LONG-RANGE FORCES
            if( bigr .gt. 15.d0 )then
                e_lr_hf = 0.d0
            else
                e_lr_hf = e_lr_hf * dexp( .2d0 / (bigr - 15.d0) )
            endif
            ! COMBINE INTERACTION REGION ENERGY AND LONG-RANGE FORCES
            e(i) = e(i) + e_lr_hf
            e(i) = e(i)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_lifh_leps_jhtp1_2002_dpem_pot11_legacy_state_a97514e1

!> @brief Compatibility interface for pot11.
subroutine pot11( r, e, nt )
    use chempotpy_lifh_leps_jhtp1_2002_dpem_pot11_legacy_state_a97514e1, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    integer          :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot11

!                      U12 COUPLING SURFACE

!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_lifh_leps_jhtp1_2002_dpem_pot12_legacy_state_a97514e1
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        implicit none
        optional :: r, e, nt
        ! MISC VARIABLES
        integer          :: i, nt
        double precision :: r(:,:), e(:)
        ! LiH VARIABLES
        integer          :: ng1
        double precision :: g1, g1r, g1a, h1, h1del, h1rho
        ! HF VARIABLES
        double precision :: j2p, j2, j2pdel, j2rho, j2prho, j2del
        ! LiF VARIABLES
        integer          :: ng3
        double precision :: g3, g3r, g3a, h3, h3del, h3rho
        ! CUTOFF VARIABLES
        double precision :: gcrho, cevau, gcdel
        double precision :: cut2, gc, cut3
        ! MISC DATA
        ! LiH DATA
        ! LiF DATA
        ! HF DATA
        ! CUTOFF DATA
        cevau  = 0.d0
        cut2   = 0.d0
        cut3   = 0.d0
        g1     = 0.d0
        g1a    = 0.d0
        g1r    = 0.d0
        g3     = 0.d0
        g3a    = 0.d0
        g3r    = 0.d0
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
        do i=1, nt
            g1   = g1a*(r(i,1)/g1r)**ng1 *dexp(-dble(ng1)*(r(i,1)/g1r-1.0d0))
            g3   = g3a*(r(i,3)/g3r)**ng3 *dexp(-dble(ng3)*(r(i,3)/g3r-1.0d0))
            j2   = 0.5d0*(1.0d0 + dtanh( (r(i,2) - j2rho)/j2del ) )
            j2p  = 0.5d0*(1.0d0 + dtanh( (r(i,2) - j2prho)/j2pdel ) )
            h1   = 0.5d0*(1.0d0 - dtanh( (r(i,1) - h1rho)/h1del ) )
            h3   = 0.5d0*(1.0d0 - dtanh( (r(i,3) - h3rho)/h3del ) )
            e(i) = g1*j2*(1.0d0-h3) + g3*j2p*(1.0d0-h1)
            gc   = r(i,2)-gcrho
            cut2 = 0.5d0 - 0.5d0*dtanh(gc/gcdel)
            cut3 = 0.5d0 - 0.5d0*dtanh((r(i,2) - 5.5d0)/.2d0)
            e(i) = e(i)*cut2*cut3
            e(i) = e(i)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_lifh_leps_jhtp1_2002_dpem_pot12_legacy_state_a97514e1

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, nt )
    use chempotpy_lifh_leps_jhtp1_2002_dpem_pot12_legacy_state_a97514e1, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    integer          :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot12

!                            U22 SURFACE

!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_lifh_leps_jhtp1_2002_dpem_pot22_legacy_state_a97514e1
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        implicit none
        optional :: r, e, nt
        ! MISC VARIABLES
        integer          :: i, j, k, nt
        double precision :: r(:,:), e(:), r1s, r2s, r3s, cs, cevau, uli2p
        double precision :: c2a, c2b, c2c, w, cplg2
        ! INTERACTION REGION VARIABLES
        !     LiH
        double precision :: r1r3, xelih, s1c, bes1, rcut2del, rcut2rho, s1c3
        double precision :: s1c1, belihc, bet1, relih, s1c2, t1c1, t1c2
        double precision :: r1r3al, rcut2
        double precision :: coul1, exch1, s1, t1, reliht
        !     HF
        double precision :: behf, xehfa, dehfm, yhf, xehfa_180, t2a3, t2a2, bet2
        double precision :: t2a1, t2a4, dehf, rehf
        double precision :: s2alp, rehft, bet2_180, s2rho, dehfc
        double precision :: coul2, exch2, t2, s2, dehfcc, s2cut, xehfm, s2m, &
            & s2del, s2cc, s2cc_d, s2cc_r, s2cc_b
        !     LiF
        integer          :: nb
        double precision :: bet3, bet3_180, t3a2, t3a3, t3a4, belihc1
        double precision :: belifff, gam, belifc2, delifc1, delifc2, belifm
        double precision :: s3c, xelihf_180, xelif, s3a, delifm, gc2, cut2
        double precision :: del2_0, del2_180, rho2_0, rho2_180, theta2_0, grho
        double precision :: theta2_180, beliff, belifmod, belif, xelif_180, relif
        double precision :: belifc1, t3a1, gswt, delif, rho2, del2, theta2, gdel
        double precision :: coul3, exch3, s3, t3, s3cut, swit, swit_r, swit_a, swit_d
        ! LONG RANGE FORCES VARIABLES
        !     HF physical properties
        double precision :: hfdip_y, hfdip_a, hfdip_re, hfdip_m(11), hfdip, &
            & hfpol, hfpolperp, hfpolpar, hfie, hfpolavg
        double precision :: hfqua_y, hfqua_a, hfqua_re, hfqua_m(10), hfqua
        !     Li(2p) physical properties
        double precision :: li2ppol, li2pie, li2pqua, li2ppolpar, &
            & li2ppolperp, li2ppolavg
        !     LiF physical properties
        double precision :: lifpol, lifie
        double precision :: lifdip_y, lifdip_re, lifdip_1, lifdip_2, &
            & lifdip_3, lifdip, lifdip_4, lifdip_beta, lifdip_d
        !     H physical properties
        double precision :: hie, hpol
        !     misc long-range variables
        double precision :: rhf0did, rhf0disp, hf_disp, e_lr_hf, e_lr_lif, &
            & rlif0did, rlif0disp, lif_disp, lif_did, hf_qq, hf_dq, rhf0qq, rhf0dq, &
            & hf_did
        double precision :: mass(3), cauang, caudeb, costhetahf, costhetalif, &
            & bigrhf, bigrlif, coschihf, coschilif, temp, r2lrcut_a, r2lrcut_r, &
            & r3lrcut_a, r3lrcut_r
        ! MISC DATA
        ! INTERACTION REGION DATA
        !     LiH
        !     HF
        !     LiF
        ! LONG RANGE FORCES DATA
        !     HF dipole moment
        !     HF quadrupole moment
        !     HF pol (au) and ie (eV)
        !     Li(2p) pol (au) and ie (eV) and quadrupole moment (au)
        !     LiF dipole moment (distance and moment in au)
        !     LiF pol (au) and ie (eV)
        !     H pol (au) and ie (eV)
        !     misc long-range
        !      data rlif0did/7.5d0/,rlif0disp/7.5d0/,r2LRcut_a/2.0d0/,
        behf        = 0.d0
        belif       = 0.d0
        belifc1     = 0.d0
        belifc2     = 0.d0
        beliff      = 0.d0
        belifff     = 0.d0
        belifm      = 0.d0
        belifmod    = 0.d0
        belihc      = 0.d0
        belihc1     = 0.d0
        bes1        = 0.d0
        bet1        = 0.d0
        bet2        = 0.d0
        bet2_180    = 0.d0
        bet3        = 0.d0
        bet3_180    = 0.d0
        bigrhf      = 0.d0
        bigrlif     = 0.d0
        c2a         = 0.d0
        c2b         = 0.d0
        c2c         = 0.d0
        cauang      = 0.d0
        caudeb      = 0.d0
        cevau       = 0.d0
        coschihf    = 0.d0
        coschilif   = 0.d0
        costhetahf  = 0.d0
        costhetalif = 0.d0
        coul1       = 0.d0
        coul2       = 0.d0
        coul3       = 0.d0
        cplg2       = 0.d0
        cs          = 0.d0
        cut2        = 0.d0
        dehf        = 0.d0
        dehfc       = 0.d0
        dehfcc      = 0.d0
        dehfm       = 0.d0
        del2        = 0.d0
        del2_0      = 0.d0
        del2_180    = 0.d0
        delif       = 0.d0
        delifc1     = 0.d0
        delifc2     = 0.d0
        delifm      = 0.d0
        e_lr_hf     = 0.d0
        e_lr_lif    = 0.d0
        exch1       = 0.d0
        exch2       = 0.d0
        exch3       = 0.d0
        gam         = 0.d0
        gc2         = 0.d0
        gdel        = 0.d0
        grho        = 0.d0
        gswt        = 0.d0
        hf_did      = 0.d0
        hf_disp     = 0.d0
        hf_dq       = 0.d0
        hf_qq       = 0.d0
        hfdip       = 0.d0
        hfdip_a     = 0.d0
        hfdip_m     = 0.d0
        hfdip_re    = 0.d0
        hfdip_y     = 0.d0
        hfie        = 0.d0
        hfpol       = 0.d0
        hfpolavg    = 0.d0
        hfpolpar    = 0.d0
        hfpolperp   = 0.d0
        hfqua       = 0.d0
        hfqua_a     = 0.d0
        hfqua_m     = 0.d0
        hfqua_re    = 0.d0
        hfqua_y     = 0.d0
        hie         = 0.d0
        hpol        = 0.d0
        li2pie      = 0.d0
        li2ppol     = 0.d0
        li2ppolavg  = 0.d0
        li2ppolpar  = 0.d0
        li2ppolperp = 0.d0
        li2pqua     = 0.d0
        lif_did     = 0.d0
        lif_disp    = 0.d0
        lifdip      = 0.d0
        lifdip_1    = 0.d0
        lifdip_2    = 0.d0
        lifdip_3    = 0.d0
        lifdip_4    = 0.d0
        lifdip_beta = 0.d0
        lifdip_d    = 0.d0
        lifdip_re   = 0.d0
        lifdip_y    = 0.d0
        lifie       = 0.d0
        lifpol      = 0.d0
        mass        = 0.d0
        r1r3        = 0.d0
        r1r3al      = 0.d0
        r1s         = 0.d0
        r2lrcut_a   = 0.d0
        r2lrcut_r   = 0.d0
        r2s         = 0.d0
        r3lrcut_a   = 0.d0
        r3lrcut_r   = 0.d0
        r3s         = 0.d0
        rcut2       = 0.d0
        rcut2del    = 0.d0
        rcut2rho    = 0.d0
        rehf        = 0.d0
        rehft       = 0.d0
        relif       = 0.d0
        relih       = 0.d0
        reliht      = 0.d0
        rhf0did     = 0.d0
        rhf0disp    = 0.d0
        rhf0dq      = 0.d0
        rhf0qq      = 0.d0
        rho2        = 0.d0
        rho2_0      = 0.d0
        rho2_180    = 0.d0
        rlif0did    = 0.d0
        rlif0disp   = 0.d0
        s1          = 0.d0
        s1c         = 0.d0
        s1c1        = 0.d0
        s1c2        = 0.d0
        s1c3        = 0.d0
        s2          = 0.d0
        s2alp       = 0.d0
        s2cc        = 0.d0
        s2cc_b      = 0.d0
        s2cc_d      = 0.d0
        s2cc_r      = 0.d0
        s2cut       = 0.d0
        s2del       = 0.d0
        s2m         = 0.d0
        s2rho       = 0.d0
        s3          = 0.d0
        s3a         = 0.d0
        s3c         = 0.d0
        s3cut       = 0.d0
        swit        = 0.d0
        swit_a      = 0.d0
        swit_d      = 0.d0
        swit_r      = 0.d0
        t1          = 0.d0
        t1c1        = 0.d0
        t1c2        = 0.d0
        t2          = 0.d0
        t2a1        = 0.d0
        t2a2        = 0.d0
        t2a3        = 0.d0
        t2a4        = 0.d0
        t3          = 0.d0
        t3a1        = 0.d0
        t3a2        = 0.d0
        t3a3        = 0.d0
        t3a4        = 0.d0
        temp        = 0.d0
        theta2      = 0.d0
        theta2_0    = 0.d0
        theta2_180  = 0.d0
        uli2p       = 0.d0
        w           = 0.d0
        xehfa       = 0.d0
        xehfa_180   = 0.d0
        xehfm       = 0.d0
        xelif       = 0.d0
        xelif_180   = 0.d0
        xelih       = 0.d0
        xelihf_180  = 0.d0
        yhf         = 0.d0
        c2a         = 3.5d0
        c2b         = 0.27362d0
        c2c         = 0.15d0
        uli2p       = 1.848d0
        cevau       = 0.036749309d0
        s1c1        = 4.32258d0
        s1c2        = 7.06452d0
        t1c1        = 1.1935483d0
        t1c2        = 13.548387d0
        belihc      = 1.36d0
        bet1        = 2.41319648d0
        relih       = 1.2d0
        r1r3al      = 1.0d0
        rcut2rho    = 0.72d0
        rcut2del    = 0.5d0
        s1c3        = 14.74194d0
        bes1        = 0.90667d0
        reliht      = 1.203323d0
        rehf        = 1.733d0
        dehf        = 6.122d0
        t2a1        = 1.4242424d0
        t2a2        = 14.203323d0
        t2a3        = 6.15835777d0
        t2a4        = 4.56304985d0
        bet2        = 2.4105572d0
        bet2_180    = 1.4046921d0
        rehft       = 1.61329423d0
        s2rho       = 1.00293d0
        s2alp       = 4.0899d0
        dehfc       = 0.0486803519d0
        dehfcc      = 0.72629521d0
        s2cc_b      = 1.799609d0
        s2cc_r      = 1.60215d0
        s2cc_d      = 2.3841642d0
        belifc1     = 0.819648d0
        belifc2     = 0.7795699d0
        delifc1     = 2.352884d0
        delifc2     = 5.5904203d0
        belifff     = 0.13225806d0
        nb          = 8
        gam         = 4.6304985d0
        grho        = 3.382209d0
        gdel        = 0.49472140d0
        t3a1        = -.0552298d0
        t3a2        = 1.8729228d0
        t3a3        = 0.94662756d0
        t3a4        = 0.2994134d0
        bet3        = 0.9519061d0
        bet3_180    = 0.4531769d0
        rho2_0      = 0.2017595d0
        rho2_180    = 2.d0
        del2_0      = 0.51710655d0
        del2_180    = 0.43695014d0
        theta2_0    = 0.12463343d0
        theta2_180  = 0.14907135d0
        delif       = 5.947d0
        relif       = 2.9553d0
        swit_r      = 3.340762d0
        swit_a      = 0.56353861d0
        swit_d      = 1.0243998d0
        hfdip_a     = 1.13d0
        hfdip_re    = 0.924212d0
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
        do i=1, nt
            ! PRECOMPUTE THINGS FOR USE LATER
            !     cosine of HF bond angle (i.e., Li-F-H angle)
            costhetahf = (r(i,3)**2-r(i,2)**2-r(i,1)**2)/(-2.d0*r(i,1)*r(i,2))
            temp       = mass(3)/(mass(2)+mass(3))*r(i,2)
            !     Li-FH Jacobi distance
            bigrhf = dsqrt(r(i,1)**2+temp**2-2.d0*r(i,1)*temp*costhetahf)
            !     cosine of HF Jacobi angle (i.e., Li--HF angle)
            coschihf = -(r(i,1)**2-bigrhf**2-temp**2)/(-2.d0*bigrhf*temp)
            coschihf = max(-1.d0, coschihf)
            coschihf = min(1.d0, coschihf)
            !     cosine of LiF bond angle (i.e., Li-F-H angle)
            costhetalif = (r(i,1)**2-r(i,2)**2-r(i,3)**2)/(-2.d0*r(i,3)*r(i,2))
            temp        = mass(1)/(mass(1)+mass(3))*r(i,3)
            !     H-LiF Jacobi distance
            bigrlif = dsqrt(r(i,2)**2+temp**2-2.d0*r(i,2)*temp*costhetalif)
            !     cosine of LiF Jacobi angle (i.e., H--LiF angle)
            coschilif = -(r(i,2)**2-bigrlif**2-temp**2)/(-2.d0*bigrlif*temp)
            coschilif = min(1.d0, coschilif)
            coschilif = max(-1.d0, coschilif)
            !     squared!
            r1s = r(i,1)**2
            r2s = r(i,2)**2
            r3s = r(i,3)**2
            !     Li-F-H bondangle
            cs = r3s+r2s-r1s
            cs = cs/(2.0d0*r(i,2)*r(i,3))
            cs = min(cs, 1.0d0)
            cs = max(cs, -1.0d0)
            ! INTERACTION REGION
            !     The LiH singlet
            xelih = dexp(-belihc*(r(i,1)-relih))
            s1    = s1c1*xelih**2  + s1c2*xelih
            s1c   = s1c3*dexp(-bes1*(r(i,1)-relih))
            r1r3  = r(i,1) - r(i,3) + r1r3al*r(i,2)
            rcut2 = 0.5d0*(1.0d0 - dtanh((r1r3-rcut2rho)/rcut2del))
            s1    = s1 + rcut2 * s1c
            !     The LiH triplet
            xelih = dexp(-bet1*(r(i,1)-reliht))
            t1    = t1c1*xelih**2  + t1c2*xelih
            !     The HF singlet
            yhf   = r(i,2) - 2.1042d0
            behf  = 1.1622d0 -0.025647d0*yhf+0.059062d0*yhf**2
            xehfa = dexp(-behf*(r(i,2)-rehf))
            s2    = dehf*xehfa*(xehfa-2.0d0) + uli2p
            s2cut = (0.5d0 - 0.5d0*tanh((r(i,2)-3.d0)/0.5d0))
            s2cc  = s2cc_d*(1.d0-dexp(-s2cc_b*(r(i,2)-s2cc_r)))**2 - s2cc_d
            s2    = s2cc + (s2 - s2cc)*s2cut
            dehfm = dehf - (dehfcc + dehfc*0.5d0*(1.0d0-cs))
            xehfm = dexp(-behf*(r(i,2)-rehf))
            s2m   = dehfm*xehfm*(xehfm-2.0d0) + uli2p
            s2cut = (0.5d0 - 0.5d0*tanh((r(i,2)-2.2d0)/0.5d0))
            s2m   = s2m*s2cut
            s2    = s2m + (s2 - s2m)*0.5d0*(1.d0+dtanh((r(i,3)-s2rho)/s2alp))
            !     The HF triplet
            xehfa     = dexp(-bet2*(r(i,2)-rehft))
            xehfa_180 = dexp(-bet2_180*(r(i,2)-rehft))
            t2        = (t2a1*xehfa**2 + t2a2*xehfa)*0.5d0*(1.0d0+cs) + (t2a3*xehfa_180**2 + t2a4*xehfa_180)* 0.5d0* &
                & (1.0d0- cs)
            !     The LiF singlet
            belif    = r(i,3)*103.57d0/(4.6498d0+7.0489d0*r(i,3))**2+0.064076d0
            gswt     = 0.5d0-0.5d0*dtanh((r(i,2)-grho)/gdel)
            gswt     = gswt * 0.5d0*(1.0d0+cs)
            beliff   = belif + gswt*(belifff - belif)
            belifmod = beliff*(belif +(r(i,3)/gam)**nb)/ (beliff+(r(i,3)/gam)**nb)
            xelif    = dexp(-belifmod*(r(i,3)-relif))
            s3a      = (delif+uli2p)*xelif*(xelif-2.0d0)+uli2p
            s3cut    = 0.5d0 - 0.5d0*tanh((r(i,3)-13.d0)/0.5d0)
            s3a      = s3a*s3cut
            swit     = swit_d  + (1.d0 - swit_d)*0.5d0 *(1.d0 + dtanh((r(i,2)-swit_r)/swit_a))
            s3a      = s3a*swit
            delifm   = delifc1 * 0.5d0*(1.0d0-cs) + delifc2 * 0.5d0*(1.0d0+cs)
            belifm   = belifc1 * 0.5d0*(1.0d0-cs) + belifc2 * 0.5d0*(1.0d0+cs)
            xelif    = dexp(-belifm*(r(i,3)-relif))
            s3c      = (delifm+uli2p)*xelif*(xelif-2.0d0)
            rho2     = rho2_0 + (rho2_180-rho2_0)*0.5d0*(1.0d0-cs)
            theta2   = theta2_0 + (theta2_180-theta2_0)*0.5d0*(1.0d0-cs)
            del2     = del2_0 + (del2_180-del2_0)*0.5d0*(1.0d0-cs)
            gc2      = dcos(theta2)*(rehf-r(i,2)) - dsin(theta2)*(rho2-r(i,3))
            cut2     = 0.5d0 - 0.5d0*dtanh(gc2/del2)
            s3       = s3c + (s3a-s3c)*cut2
            !     The LiF triplet
            xelif     = dexp(-bet3*(r(i,3)-relif))
            xelif_180 = dexp(-bet3_180*(r(i,3)-relif))
            t3        = (t3a1*xelif**2 + t3a2*xelif)*0.5d0*(1.0d0+cs) + (t3a3*xelif_180**2 + t3a4*xelif_180)* 0.5d0* &
                & (1.0d0- cs)
            !     LEPS
            coul1 = 0.5d0*(s1+t1)
            coul2 = 0.5d0*(s2+t2)
            coul3 = 0.5d0*(s3+t3)
            exch1 = 0.5d0*(s1-t1)
            exch2 = 0.5d0*(s2-t2)
            exch3 = 0.5d0*(s3-t3)
            w     = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
            cplg2 = c2a*dexp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
            e(i)  = coul1+coul2+coul3-dsqrt(w+cplg2*cplg2)/dsqrt(2.0d0) + dehf
            ! LONG-RANGE FORCES * * * * * * * * * * * * * * * * * * * * * * * * * *
            ! HF DIPOLE
            hfdip_y = 1.d0 - dexp(-hfdip_a*(r(i,2)*cauang-hfdip_re))
            hfdip   = 0.d0
            do j=1, 11
                hfdip = hfdip + hfdip_m(j)*hfdip_y**dble(j-1)
            enddo
            hfdip = hfdip*0.5d0*(1.d0 + dtanh(20.d0*(r(i,2)-1.2d0)))
            ! HF POLARIZATION
            hfpol    = hfpolperp*(1.d0 - coschihf**2) + hfpolpar*coschihf**2
            hfpolavg = 2.d0*hfpolperp/3.d0 + hfpolpar/3.d0
            ! HF QUADRUPOLE
            hfqua_y = 1.d0 - dexp(-hfqua_a*(r(i,2)*cauang-hfqua_re))
            hfqua   = 0.d0
            do j=1, 10
                hfqua = hfqua + hfqua_m(j)*hfqua_y**dble(j-1)
            enddo
            hfqua = hfqua*0.5d0*(1.d0 + dtanh(20.d0*(r(i,2)-1.2d0)))
            ! LiF DIPOLE
            lifdip_y    = r(i,3) - lifdip_4
            lifdip_beta = lifdip_1 + lifdip_2*lifdip_y**8 + lifdip_3*lifdip_y
            lifdip      = lifdip_d*dexp(-lifdip_beta*(r(i,3) - lifdip_re)**2)
            ! Li(2p) POLARIZATION
            li2ppol    = li2ppolperp*(1.d0-coschihf**2)+li2ppolpar*coschihf**2
            li2ppolavg = 2.d0*li2ppolperp/3.d0 + li2ppolpar/3.d0
            ! HF + Li(2p)
            !     DIPOLE-INDUCED-DIPOLE
            hf_did = -0.5d0*(3.d0*coschihf**2 + 1.d0)*li2ppol*hfdip**2
            hf_did = hf_did/(bigrhf**7 + rhf0did**7)/cevau*bigrhf
            !     cut off angular dependence in interaction region
            hf_did = hf_did*(1.d0 + (2.5d0*li2ppolavg/ (3.d0*coschihf**2 + 1.d0)/li2ppol - 1.d0) *dexp(-(bigrhf/ &
                & rhf0did)** 4))
            !     DISPERSION
            hf_disp = -1.5d0*li2pie*hfie*li2ppol*hfpol/(li2pie+hfie)
            hf_disp = hf_disp/(bigrhf**7 + rhf0disp**7)*bigrhf
            !     cut off angular dependence in interaction region
            hf_disp = hf_disp*(1.d0+(li2ppolavg*hfpolavg/hfpol/li2ppol-1.d0) *dexp(-(bigrhf/rhf0disp)**4))
            !     QUADRUPOLE-QUADRUPOLE
            hf_qq = 3.d0/4.d0*hfqua*li2pqua*(3.d0 - 7.d0*coschihf**2)
            hf_qq = hf_qq/(bigrhf**6 + rhf0qq**6)/cevau*bigrhf
            !     cut off angular dependence in interaction region
            hf_qq = hf_qq*(1.d0 + (-.5d0/(3.d0-7.d0*coschihf**2) - 1.d0) *dexp(-(bigrhf/rhf0qq)**4))
            !     DIPOLE-QUADRUPOLE
            hf_dq = 3.d0/4.d0*hfdip*li2pqua*coschihf
            hf_dq = hf_dq/(bigrhf**5 + rhf0dq**5)/cevau*bigrhf
            !     cut off angular dependence in interaction region
            hf_qq = hf_qq*dexp(-(bigrhf/rhf0dq)**4)
            !     SUM FORCES
            e_lr_hf = hf_dq + hf_qq + hf_disp + hf_did
            !     cut off at large HF distances
            e_lr_hf = e_lr_hf*0.5d0*(1.d0-dtanh(r2lrcut_a*(r(i,2)-r2lrcut_r)))
            ! LiF + H
            !     DIPOLE-INDUCED-DIPOLE
            lif_did = -0.5d0*(3.d0*coschilif**2+1.d0)*hpol*lifdip**2
            lif_did = lif_did/(bigrlif**6 + rlif0did**6)/cevau
            !     cut off angular dependence in interaction region
            lif_did = lif_did*(1.d0+(2.5d0/(3.d0*coschilif**2+1.d0)-1.d0) *dexp(-(bigrlif/rlif0did)**4))
            !     DISPERSION
            lif_disp = -1.5d0*lifie*hie*lifpol*hpol/(lifie+hie)
            lif_disp = lif_disp/(bigrlif**6 + rlif0disp**6)
            !     no angular dependence
            !     SUM FORCES
            e_lr_lif = lif_disp + lif_did
            !     cut off at large LiF distances
            e_lr_lif = e_lr_lif*0.5d0*(1.d0-dtanh(r3lrcut_a*(r(i,3)-r3lrcut_r)))
            ! LiH + F
            ! no long-range forces included
            ! CUT OFF LONG-RANGE FORCES
            if( bigrhf .gt. 15.d0 )then
                e_lr_hf = 0.d0
            else
                e_lr_hf = e_lr_hf * dexp( .2d0 / (bigrhf - 15.d0) )
            endif
            if( bigrlif .gt. 15.d0 )then
                e_lr_lif = 0.d0
            else
                e_lr_lif = e_lr_lif * dexp( .2d0 / (bigrlif - 15.d0) )
            endif
            ! COMBINE INTERACTION REGION ENERGY AND LONG-RANGE FORCES
            e(i) = e(i) + e_lr_lif + e_lr_hf
            e(i) = e(i)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_lifh_leps_jhtp1_2002_dpem_pot22_legacy_state_a97514e1

!> @brief Compatibility interface for pot22.
subroutine pot22( r, e, nt )
    use chempotpy_lifh_leps_jhtp1_2002_dpem_pot22_legacy_state_a97514e1, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    integer          :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot22
