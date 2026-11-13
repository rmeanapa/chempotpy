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
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the LiFH_LEPS_JHTP1_2002 diabatic potential matrix for LiFH.
    subroutine evaluate( x, igrad, u, ug )
        use chempotpy_units,              only: hartree_to_electronvolt
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
        real(wp)              :: nt, r(1,3), v(1)
        real(wp)              :: u11(1), u12(1), u22(1)
        integer               :: i, j, k, l
        !initialize
        u  = 0.e0_wp; ug = 0.e0_wp
        nt = 1
        ! input cartesian is ClHH
        call triatomic_distances_to_bohr(x, r(1,:))
        if( igrad == 0 )then
            call pot(r, u11, 1, 1)
            call pot(r, u12, 1, 2)
            call pot(r, u22, 1, 3)
            u(1,1) = u11(1)*hartree_to_electronvolt
            u(1,2) = u12(1)*hartree_to_electronvolt
            u(2,1) = u12(1)*hartree_to_electronvolt
            u(2,2) = u22(1)*hartree_to_electronvolt
        else
            write(*, *) "Only DPEM is available"
        endif
    end subroutine evaluate
end module chempotpy_lifh_leps_jhtp1_2002_dpem_pes

!> @brief Expose LiFH_LEPS_JHTP1_2002 through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_lifh_leps_jhtp1_2002_dpem_pes, only: evaluate
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

!> @brief Evaluate EvdRdX for the LiFH_LEPS_JHTP1_2002_DPEM potential.
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

!> @brief Evaluate diagonalize for the LiFH_LEPS_JHTP1_2002_DPEM potential.
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

!> @brief Evaluate pot for the LiFH_LEPS_JHTP1_2002_DPEM potential.
subroutine pot( r, e, nt, nsurf )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp) :: r(nt,3), e(nt)
    integer  :: i, nt, nsurf
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
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r, e, nt
        ! MISC VARIABLES
        integer  :: i, nt, j
        real(wp) :: cevau
        real(wp) :: r(:,:), e(:)
        real(wp) :: r1s, r2s, r3s
        ! INTERACTION REGION VARIABLES
        !     LiH
        real(wp) :: s1, relih, bes1, s1c2, s1c3, s1a, xelih, bes1b, s1c, &
            & rcut2del, rcut2rho, rcut2, r1r3, r1r3al
        !     HF
        real(wp) :: s2, behf, yhf, rehf, dehf, xehf, s2m, behfm, dehfm, &
            & rehfm, xehfm, s2rho, s2alp, dehfc, dehfcc, yhf_r, behf_0, behf_1, &
            & behf_2
        !     LiF
        real(wp) :: s3, res3, bes3, des3, res3c, bes3c, s3c, s3a, xes3, des3c, &
            & rcutdel, rcutrho, rcut, r1r3al2, bes3i, bes3z, bes3r
        ! LONG RANGE FORCES VARIABLES
        !     HF physical properties
        real(wp) :: hfdip_y, hfdip_a, hfdip_re, hfdip_m(11), hfdip, &
            & hfpol, hfpolperp, hfpolpar, hfie, hfpolavg
        !     Li physical properties
        real(wp) :: li2spol, rhf0did, rhf0disp, li2sie
        !     misc long-range variables
        real(wp) :: hf_did, hf_disp, e_lr_hf
        real(wp) :: mass(3), cauang, caudeb, costheta, temp, bigr, &
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
        behf      = 0.e0_wp
        behf_0    = 0.e0_wp
        behf_1    = 0.e0_wp
        behf_2    = 0.e0_wp
        behfm     = 0.e0_wp
        bes1      = 0.e0_wp
        bes1b     = 0.e0_wp
        bes3      = 0.e0_wp
        bes3c     = 0.e0_wp
        bes3i     = 0.e0_wp
        bes3r     = 0.e0_wp
        bes3z     = 0.e0_wp
        bigr      = 0.e0_wp
        cauang    = 0.e0_wp
        caudeb    = 0.e0_wp
        cevau     = 0.e0_wp
        coschi    = 0.e0_wp
        costheta  = 0.e0_wp
        dehf      = 0.e0_wp
        dehfc     = 0.e0_wp
        dehfcc    = 0.e0_wp
        dehfm     = 0.e0_wp
        des3      = 0.e0_wp
        des3c     = 0.e0_wp
        e_lr_hf   = 0.e0_wp
        hf_did    = 0.e0_wp
        hf_disp   = 0.e0_wp
        hfdip     = 0.e0_wp
        hfdip_a   = 0.e0_wp
        hfdip_m   = 0.e0_wp
        hfdip_re  = 0.e0_wp
        hfdip_y   = 0.e0_wp
        hfie      = 0.e0_wp
        hfpol     = 0.e0_wp
        hfpolavg  = 0.e0_wp
        hfpolpar  = 0.e0_wp
        hfpolperp = 0.e0_wp
        li2sie    = 0.e0_wp
        li2spol   = 0.e0_wp
        mass      = 0.e0_wp
        r1r3      = 0.e0_wp
        r1r3al    = 0.e0_wp
        r1r3al2   = 0.e0_wp
        r1s       = 0.e0_wp
        r2lrcut_a = 0.e0_wp
        r2lrcut_r = 0.e0_wp
        r2s       = 0.e0_wp
        r3s       = 0.e0_wp
        rcut      = 0.e0_wp
        rcut2     = 0.e0_wp
        rcut2del  = 0.e0_wp
        rcut2rho  = 0.e0_wp
        rcutdel   = 0.e0_wp
        rcutrho   = 0.e0_wp
        rehf      = 0.e0_wp
        rehfm     = 0.e0_wp
        relih     = 0.e0_wp
        res3      = 0.e0_wp
        res3c     = 0.e0_wp
        rhf0did   = 0.e0_wp
        rhf0disp  = 0.e0_wp
        s1        = 0.e0_wp
        s1a       = 0.e0_wp
        s1c       = 0.e0_wp
        s1c2      = 0.e0_wp
        s1c3      = 0.e0_wp
        s2        = 0.e0_wp
        s2alp     = 0.e0_wp
        s2m       = 0.e0_wp
        s2rho     = 0.e0_wp
        s3        = 0.e0_wp
        s3a       = 0.e0_wp
        s3c       = 0.e0_wp
        temp      = 0.e0_wp
        xehf      = 0.e0_wp
        xehfm     = 0.e0_wp
        xelih     = 0.e0_wp
        xes3      = 0.e0_wp
        yhf       = 0.e0_wp
        yhf_r     = 0.e0_wp
        cevau     = 0.036749309e0_wp
        relih     = 1.3432062e0_wp
        bes1      = 1.8058651e0_wp
        s1c2      = 12.937438e0_wp
        s1c3      = 11.0913e0_wp
        bes1b     = 1.2111436e0_wp
        rcut2rho  = 1.07143e0_wp
        rcut2del  = 0.6e0_wp
        r1r3al    = 0.4e0_wp
        rehf      = 1.733e0_wp
        dehf      = 6.122e0_wp
        rehfm     = 1.6739e0_wp
        dehfc     = 0.344574e0_wp
        dehfcc    = 0.8660801e0_wp
        s2rho     = 1.63e0_wp
        s2alp     = 2.9941e0_wp
        behfm     = 0.74633431e0_wp
        yhf_r     = 2.1042e0_wp
        behf_0    = 1.1622e0_wp
        behf_1    = 0.025647e0_wp
        behf_2    = 0.059062e0_wp
        res3      = 3.6e0_wp
        bes3      = 1.49061584e0_wp
        des3      = 0.17426687e0_wp
        res3c     = 2.48571e0_wp
        bes3c     = 1.92857e0_wp
        des3c     = 0.25e0_wp
        rcutrho   = 0.4e0_wp
        rcutdel   = 0.89333e0_wp
        r1r3al2   = 0.7333e0_wp
        bes3z     = 0.e0_wp
        bes3r     = 7.4799608e0_wp
        hfdip_a   = 1.13e0_wp
        hfdip_re  = 0.924212e0_wp
        hfdip_m = [ 0.703661e0_wp, 0.516815e0_wp, 0.240628e0_wp, -0.430194e0_wp, -2.213088e0_wp, 5.535609e0_wp, &
            & 12.872106e0_wp, - 42.060172e0_wp, -17.398065e0_wp, 84.207629e0_wp, -41.961465e0_wp]
        hfpolperp = 4.59e0_wp
        hfpolpar  = 5.10e0_wp
        hfie      = 16.044e0_wp
        li2spol   = 165.e0_wp
        li2sie    = 5.392e0_wp
        rhf0did   = 6.0e0_wp
        rhf0disp  = 6.0e0_wp
        r2lrcut_a = 2.0e0_wp
        r2lrcut_r = 3.e0_wp
        cauang    = 0.52917706e0_wp
        caudeb    = 2.54177e0_wp
        mass      = [ 7.016003e0_wp, 1.00783e0_wp, 18.9984e0_wp]
        do i=1, nt
            ! PRECOMPUTE THINGS FOR USE LATER...
            !     cosine of HF bond angle (i.e., Li-F-H angle)
            costheta = (r(i,3)**2-r(i,2)**2-r(i,1)**2)/(-2.e0_wp*r(i,1)*r(i,2))
            temp     = mass(3)/(mass(2)+mass(3))*r(i,2)
            !     Li-FH Jacobi distance
            bigr = sqrt(r(i,1)**2+temp**2-2.e0_wp*r(i,1)*temp*costheta)
            !     cosine of HF Jacobi angle (i.e., Li--HF angle)
            coschi = -(r(i,1)**2-bigr**2-temp**2)/(-2.e0_wp*bigr*temp)
            coschi = max(coschi, -1.e0_wp)
            coschi = min(coschi, 1.e0_wp)
            !     squared!
            r1s = r(i,1)**2
            r2s = r(i,2)**2
            r3s = r(i,3)**2
            ! INTERACTION REGION
            !     LiH
            xelih = exp(-bes1*(r(i,1)-relih))
            s1a   = s1c2*xelih
            s1c   = s1c3*exp(-bes1b*(r(i,1)-relih))
            r1r3  = r(i,1) - r(i,3) + r1r3al*r(i,2)
            rcut2 = 0.5e0_wp*(1.0e0_wp - tanh((r1r3-rcut2rho)/rcut2del))
            s1    = s1a  + (s1c-s1a)*rcut2
            !     HF
            yhf   = r(i,2) - yhf_r
            behf  = behf_0-behf_1*yhf+behf_2*yhf**2
            xehf  = exp(-behf*(r(i,2)-rehf))
            s2    = dehf*xehf*(xehf-2.0e0_wp)
            dehfm = dehf - (dehfcc+dehfc*0.5e0_wp*(1.e0_wp-coschi))
            xehfm = exp(-behfm*(r(i,2)-rehfm))
            s2m   = dehfm*xehfm*(xehfm-2.0e0_wp)
            s2    = s2m + (s2 - s2m)*0.5e0_wp*(1.e0_wp+tanh((r(i,3)-s2rho)/s2alp))
            !     LiF
            xes3 = exp(-bes3*(r(i,3)-res3))
            s3   = des3*xes3*(xes3-2.0e0_wp)
            !     Sum of diatomics
            e(i) = s1 + s2 + s3 + dehf
            ! EVALUATE THE LONG-RANGE FORCES
            ! HF DIPOLE MOMENT
            hfdip_y = 1.e0_wp - exp(-hfdip_a*(r(i,2)*cauang-hfdip_re))
            hfdip   = 0.e0_wp
            do j=1, 11
                hfdip = hfdip + hfdip_m(j)*hfdip_y**real(j-1, wp)
            enddo
            hfdip = hfdip*0.5e0_wp*(1.e0_wp + tanh(20.e0_wp*(r(i,2)-1.2e0_wp)))
            ! HF POLARIZATION
            hfpol    = hfpolperp*(1.e0_wp - coschi**2) + hfpolpar*coschi**2
            hfpolavg = 2.e0_wp*hfpolperp/3.e0_wp + hfpolpar/3.e0_wp
            ! HF + Li(2s) arrangement
            !     DIPOLE-INDUCED-DIPOLE
            hf_did = -0.5e0_wp*(3.e0_wp*coschi**2 + 1.e0_wp)*li2spol*hfdip**2
            hf_did = hf_did/(bigr**7 + rhf0did**7)/cevau*bigr
            !     cut off angular dependence in interaction region
            hf_did = hf_did*(1.e0_wp + (2.5e0_wp/(3.e0_wp*coschi**2 + 1.e0_wp) - 1.e0_wp) *exp(-(bigr/rhf0did)**4))
            !     DISPERSION
            hf_disp = -1.5e0_wp*li2sie*hfie*li2spol*hfpol/(li2sie+hfie)
            hf_disp = hf_disp/(bigr**7 + rhf0disp**7)*bigr
            !     cut off angular dependence in interaction region
            hf_disp = hf_disp*(1.e0_wp + (hfpolavg/hfpol - 1.e0_wp) *exp(-(bigr/rhf0disp)**4))
            !     SUM FORCES
            e_lr_hf = hf_disp + hf_did
            !     cut off LR forces for large r[hf] i.e., r(i,2)
            e_lr_hf = e_lr_hf*0.5e0_wp*(1.e0_wp-tanh(r2lrcut_a*(r(i,2)-r2lrcut_r)))
            ! LiF + H
            !     NO LONG-RANGE FORCES
            ! LiH + F
            !     NO LONG-RANGE FORCES
            ! CUT OFF LONG-RANGE FORCES
            if( bigr .gt. 15.e0_wp )then
                e_lr_hf = 0.e0_wp
            else
                e_lr_hf = e_lr_hf * exp( .2e0_wp / (bigr - 15.e0_wp) )
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
    use chempotpy_kinds,                                                 only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot11

!                      U12 COUPLING SURFACE

!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_lifh_leps_jhtp1_2002_dpem_pot12_legacy_state_a97514e1
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r, e, nt
        ! MISC VARIABLES
        integer  :: i, nt
        real(wp) :: r(:,:), e(:)
        ! LiH VARIABLES
        integer  :: ng1
        real(wp) :: g1, g1r, g1a, h1, h1del, h1rho
        ! HF VARIABLES
        real(wp) :: j2p, j2, j2pdel, j2rho, j2prho, j2del
        ! LiF VARIABLES
        integer  :: ng3
        real(wp) :: g3, g3r, g3a, h3, h3del, h3rho
        ! CUTOFF VARIABLES
        real(wp) :: gcrho, cevau, gcdel
        real(wp) :: cut2, gc, cut3
        ! MISC DATA
        ! LiH DATA
        ! LiF DATA
        ! HF DATA
        ! CUTOFF DATA
        cevau  = 0.e0_wp
        cut2   = 0.e0_wp
        cut3   = 0.e0_wp
        g1     = 0.e0_wp
        g1a    = 0.e0_wp
        g1r    = 0.e0_wp
        g3     = 0.e0_wp
        g3a    = 0.e0_wp
        g3r    = 0.e0_wp
        gc     = 0.e0_wp
        gcdel  = 0.e0_wp
        gcrho  = 0.e0_wp
        h1     = 0.e0_wp
        h1del  = 0.e0_wp
        h1rho  = 0.e0_wp
        h3     = 0.e0_wp
        h3del  = 0.e0_wp
        h3rho  = 0.e0_wp
        j2     = 0.e0_wp
        j2del  = 0.e0_wp
        j2p    = 0.e0_wp
        j2pdel = 0.e0_wp
        j2prho = 0.e0_wp
        j2rho  = 0.e0_wp
        cevau  = 0.036749309e0_wp
        g1a    = 1.27742e0_wp
        g1r    = 3.0048876e0_wp
        ng1    = 6
        h1rho  = 4.9843597e0_wp
        h1del  = 2.0899315e0_wp
        g3a    = 0.48e0_wp
        g3r    = 3.4975562e0_wp
        ng3    = 8
        h3rho  = 2.000978e0_wp
        h3del  = 0.90615835e0_wp
        j2rho  = 1.15484e0_wp
        j2del  = 1.75806e0_wp
        j2prho = 1.45161e0_wp
        j2pdel = 0.98387e0_wp
        gcrho  = 4.87097e0_wp
        gcdel  = 2.e0_wp
        do i=1, nt
            g1   = g1a*(r(i,1)/g1r)**ng1 *exp(-real(ng1, wp)*(r(i,1)/g1r-1.0e0_wp))
            g3   = g3a*(r(i,3)/g3r)**ng3 *exp(-real(ng3, wp)*(r(i,3)/g3r-1.0e0_wp))
            j2   = 0.5e0_wp*(1.0e0_wp + tanh( (r(i,2) - j2rho)/j2del ) )
            j2p  = 0.5e0_wp*(1.0e0_wp + tanh( (r(i,2) - j2prho)/j2pdel ) )
            h1   = 0.5e0_wp*(1.0e0_wp - tanh( (r(i,1) - h1rho)/h1del ) )
            h3   = 0.5e0_wp*(1.0e0_wp - tanh( (r(i,3) - h3rho)/h3del ) )
            e(i) = g1*j2*(1.0e0_wp-h3) + g3*j2p*(1.0e0_wp-h1)
            gc   = r(i,2)-gcrho
            cut2 = 0.5e0_wp - 0.5e0_wp*tanh(gc/gcdel)
            cut3 = 0.5e0_wp - 0.5e0_wp*tanh((r(i,2) - 5.5e0_wp)/.2e0_wp)
            e(i) = e(i)*cut2*cut3
            e(i) = e(i)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_lifh_leps_jhtp1_2002_dpem_pot12_legacy_state_a97514e1

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, nt )
    use chempotpy_lifh_leps_jhtp1_2002_dpem_pot12_legacy_state_a97514e1, only: dispatch
    use chempotpy_kinds,                                                 only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot12

!                            U22 SURFACE

!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_lifh_leps_jhtp1_2002_dpem_pot22_legacy_state_a97514e1
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r, e, nt
        ! MISC VARIABLES
        integer  :: i, j, k, nt
        real(wp) :: r(:,:), e(:), r1s, r2s, r3s, cs, cevau, uli2p
        real(wp) :: c2a, c2b, c2c, w, cplg2
        ! INTERACTION REGION VARIABLES
        !     LiH
        real(wp) :: r1r3, xelih, s1c, bes1, rcut2del, rcut2rho, s1c3
        real(wp) :: s1c1, belihc, bet1, relih, s1c2, t1c1, t1c2
        real(wp) :: r1r3al, rcut2
        real(wp) :: coul1, exch1, s1, t1, reliht
        !     HF
        real(wp) :: behf, xehfa, dehfm, yhf, xehfa_180, t2a3, t2a2, bet2
        real(wp) :: t2a1, t2a4, dehf, rehf
        real(wp) :: s2alp, rehft, bet2_180, s2rho, dehfc
        real(wp) :: coul2, exch2, t2, s2, dehfcc, s2cut, xehfm, s2m, &
            & s2del, s2cc, s2cc_d, s2cc_r, s2cc_b
        !     LiF
        integer  :: nb
        real(wp) :: bet3, bet3_180, t3a2, t3a3, t3a4, belihc1
        real(wp) :: belifff, gam, belifc2, delifc1, delifc2, belifm
        real(wp) :: s3c, xelihf_180, xelif, s3a, delifm, gc2, cut2
        real(wp) :: del2_0, del2_180, rho2_0, rho2_180, theta2_0, grho
        real(wp) :: theta2_180, beliff, belifmod, belif, xelif_180, relif
        real(wp) :: belifc1, t3a1, gswt, delif, rho2, del2, theta2, gdel
        real(wp) :: coul3, exch3, s3, t3, s3cut, swit, swit_r, swit_a, swit_d
        ! LONG RANGE FORCES VARIABLES
        !     HF physical properties
        real(wp) :: hfdip_y, hfdip_a, hfdip_re, hfdip_m(11), hfdip, &
            & hfpol, hfpolperp, hfpolpar, hfie, hfpolavg
        real(wp) :: hfqua_y, hfqua_a, hfqua_re, hfqua_m(10), hfqua
        !     Li(2p) physical properties
        real(wp) :: li2ppol, li2pie, li2pqua, li2ppolpar, &
            & li2ppolperp, li2ppolavg
        !     LiF physical properties
        real(wp) :: lifpol, lifie
        real(wp) :: lifdip_y, lifdip_re, lifdip_1, lifdip_2, &
            & lifdip_3, lifdip, lifdip_4, lifdip_beta, lifdip_d
        !     H physical properties
        real(wp) :: hie, hpol
        !     misc long-range variables
        real(wp) :: rhf0did, rhf0disp, hf_disp, e_lr_hf, e_lr_lif, &
            & rlif0did, rlif0disp, lif_disp, lif_did, hf_qq, hf_dq, rhf0qq, rhf0dq, &
            & hf_did
        real(wp) :: mass(3), cauang, caudeb, costhetahf, costhetalif, &
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
        bigrhf      = 0.e0_wp
        bigrlif     = 0.e0_wp
        c2a         = 0.e0_wp
        c2b         = 0.e0_wp
        c2c         = 0.e0_wp
        cauang      = 0.e0_wp
        caudeb      = 0.e0_wp
        cevau       = 0.e0_wp
        coschihf    = 0.e0_wp
        coschilif   = 0.e0_wp
        costhetahf  = 0.e0_wp
        costhetalif = 0.e0_wp
        coul1       = 0.e0_wp
        coul2       = 0.e0_wp
        coul3       = 0.e0_wp
        cplg2       = 0.e0_wp
        cs          = 0.e0_wp
        cut2        = 0.e0_wp
        dehf        = 0.e0_wp
        dehfc       = 0.e0_wp
        dehfcc      = 0.e0_wp
        dehfm       = 0.e0_wp
        del2        = 0.e0_wp
        del2_0      = 0.e0_wp
        del2_180    = 0.e0_wp
        delif       = 0.e0_wp
        delifc1     = 0.e0_wp
        delifc2     = 0.e0_wp
        delifm      = 0.e0_wp
        e_lr_hf     = 0.e0_wp
        e_lr_lif    = 0.e0_wp
        exch1       = 0.e0_wp
        exch2       = 0.e0_wp
        exch3       = 0.e0_wp
        gam         = 0.e0_wp
        gc2         = 0.e0_wp
        gdel        = 0.e0_wp
        grho        = 0.e0_wp
        gswt        = 0.e0_wp
        hf_did      = 0.e0_wp
        hf_disp     = 0.e0_wp
        hf_dq       = 0.e0_wp
        hf_qq       = 0.e0_wp
        hfdip       = 0.e0_wp
        hfdip_a     = 0.e0_wp
        hfdip_m     = 0.e0_wp
        hfdip_re    = 0.e0_wp
        hfdip_y     = 0.e0_wp
        hfie        = 0.e0_wp
        hfpol       = 0.e0_wp
        hfpolavg    = 0.e0_wp
        hfpolpar    = 0.e0_wp
        hfpolperp   = 0.e0_wp
        hfqua       = 0.e0_wp
        hfqua_a     = 0.e0_wp
        hfqua_m     = 0.e0_wp
        hfqua_re    = 0.e0_wp
        hfqua_y     = 0.e0_wp
        hie         = 0.e0_wp
        hpol        = 0.e0_wp
        li2pie      = 0.e0_wp
        li2ppol     = 0.e0_wp
        li2ppolavg  = 0.e0_wp
        li2ppolpar  = 0.e0_wp
        li2ppolperp = 0.e0_wp
        li2pqua     = 0.e0_wp
        lif_did     = 0.e0_wp
        lif_disp    = 0.e0_wp
        lifdip      = 0.e0_wp
        lifdip_1    = 0.e0_wp
        lifdip_2    = 0.e0_wp
        lifdip_3    = 0.e0_wp
        lifdip_4    = 0.e0_wp
        lifdip_beta = 0.e0_wp
        lifdip_d    = 0.e0_wp
        lifdip_re   = 0.e0_wp
        lifdip_y    = 0.e0_wp
        lifie       = 0.e0_wp
        lifpol      = 0.e0_wp
        mass        = 0.e0_wp
        r1r3        = 0.e0_wp
        r1r3al      = 0.e0_wp
        r1s         = 0.e0_wp
        r2lrcut_a   = 0.e0_wp
        r2lrcut_r   = 0.e0_wp
        r2s         = 0.e0_wp
        r3lrcut_a   = 0.e0_wp
        r3lrcut_r   = 0.e0_wp
        r3s         = 0.e0_wp
        rcut2       = 0.e0_wp
        rcut2del    = 0.e0_wp
        rcut2rho    = 0.e0_wp
        rehf        = 0.e0_wp
        rehft       = 0.e0_wp
        relif       = 0.e0_wp
        relih       = 0.e0_wp
        reliht      = 0.e0_wp
        rhf0did     = 0.e0_wp
        rhf0disp    = 0.e0_wp
        rhf0dq      = 0.e0_wp
        rhf0qq      = 0.e0_wp
        rho2        = 0.e0_wp
        rho2_0      = 0.e0_wp
        rho2_180    = 0.e0_wp
        rlif0did    = 0.e0_wp
        rlif0disp   = 0.e0_wp
        s1          = 0.e0_wp
        s1c         = 0.e0_wp
        s1c1        = 0.e0_wp
        s1c2        = 0.e0_wp
        s1c3        = 0.e0_wp
        s2          = 0.e0_wp
        s2alp       = 0.e0_wp
        s2cc        = 0.e0_wp
        s2cc_b      = 0.e0_wp
        s2cc_d      = 0.e0_wp
        s2cc_r      = 0.e0_wp
        s2cut       = 0.e0_wp
        s2del       = 0.e0_wp
        s2m         = 0.e0_wp
        s2rho       = 0.e0_wp
        s3          = 0.e0_wp
        s3a         = 0.e0_wp
        s3c         = 0.e0_wp
        s3cut       = 0.e0_wp
        swit        = 0.e0_wp
        swit_a      = 0.e0_wp
        swit_d      = 0.e0_wp
        swit_r      = 0.e0_wp
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
        temp        = 0.e0_wp
        theta2      = 0.e0_wp
        theta2_0    = 0.e0_wp
        theta2_180  = 0.e0_wp
        uli2p       = 0.e0_wp
        w           = 0.e0_wp
        xehfa       = 0.e0_wp
        xehfa_180   = 0.e0_wp
        xehfm       = 0.e0_wp
        xelif       = 0.e0_wp
        xelif_180   = 0.e0_wp
        xelih       = 0.e0_wp
        xelihf_180  = 0.e0_wp
        yhf         = 0.e0_wp
        c2a         = 3.5e0_wp
        c2b         = 0.27362e0_wp
        c2c         = 0.15e0_wp
        uli2p       = 1.848e0_wp
        cevau       = 0.036749309e0_wp
        s1c1        = 4.32258e0_wp
        s1c2        = 7.06452e0_wp
        t1c1        = 1.1935483e0_wp
        t1c2        = 13.548387e0_wp
        belihc      = 1.36e0_wp
        bet1        = 2.41319648e0_wp
        relih       = 1.2e0_wp
        r1r3al      = 1.0e0_wp
        rcut2rho    = 0.72e0_wp
        rcut2del    = 0.5e0_wp
        s1c3        = 14.74194e0_wp
        bes1        = 0.90667e0_wp
        reliht      = 1.203323e0_wp
        rehf        = 1.733e0_wp
        dehf        = 6.122e0_wp
        t2a1        = 1.4242424e0_wp
        t2a2        = 14.203323e0_wp
        t2a3        = 6.15835777e0_wp
        t2a4        = 4.56304985e0_wp
        bet2        = 2.4105572e0_wp
        bet2_180    = 1.4046921e0_wp
        rehft       = 1.61329423e0_wp
        s2rho       = 1.00293e0_wp
        s2alp       = 4.0899e0_wp
        dehfc       = 0.0486803519e0_wp
        dehfcc      = 0.72629521e0_wp
        s2cc_b      = 1.799609e0_wp
        s2cc_r      = 1.60215e0_wp
        s2cc_d      = 2.3841642e0_wp
        belifc1     = 0.819648e0_wp
        belifc2     = 0.7795699e0_wp
        delifc1     = 2.352884e0_wp
        delifc2     = 5.5904203e0_wp
        belifff     = 0.13225806e0_wp
        nb          = 8
        gam         = 4.6304985e0_wp
        grho        = 3.382209e0_wp
        gdel        = 0.49472140e0_wp
        t3a1        = -.0552298e0_wp
        t3a2        = 1.8729228e0_wp
        t3a3        = 0.94662756e0_wp
        t3a4        = 0.2994134e0_wp
        bet3        = 0.9519061e0_wp
        bet3_180    = 0.4531769e0_wp
        rho2_0      = 0.2017595e0_wp
        rho2_180    = 2.e0_wp
        del2_0      = 0.51710655e0_wp
        del2_180    = 0.43695014e0_wp
        theta2_0    = 0.12463343e0_wp
        theta2_180  = 0.14907135e0_wp
        delif       = 5.947e0_wp
        relif       = 2.9553e0_wp
        swit_r      = 3.340762e0_wp
        swit_a      = 0.56353861e0_wp
        swit_d      = 1.0243998e0_wp
        hfdip_a     = 1.13e0_wp
        hfdip_re    = 0.924212e0_wp
        hfdip_m = [ 0.703661e0_wp, 0.516815e0_wp, 0.240628e0_wp, -0.430194e0_wp, -2.213088e0_wp, 5.535609e0_wp, &
            & 12.872106e0_wp, - 42.060172e0_wp, -17.398065e0_wp, 84.207629e0_wp, -41.961465e0_wp]
        hfqua_a  = 1.37e0_wp
        hfqua_re = 0.9168e0_wp
        hfqua_m = [ 1.6527e0_wp, 2.072e0_wp, 2.0521e0_wp, 3.545e0_wp, -9.8031e0_wp, -19.3808e0_wp, 76.4423e0_wp, &
            & 0.e0_wp, -161.5896e0_wp, 105.7460e0_wp]
        hfpolperp   = 4.59e0_wp
        hfpolpar    = 5.10e0_wp
        hfie        = 16.044e0_wp
        li2ppolpar  = 131.e0_wp
        li2ppolperp = 129.e0_wp
        li2pie      = 3.544e0_wp
        li2pqua     = 11.1e0_wp
        lifdip_re   = 10.4994e0_wp
        lifdip_1    = 0.02435e0_wp
        lifdip_2    = 9.9355e-8_wp
        lifdip_3    = 0.0015999e0_wp
        lifdip_4    = 4.471959e0_wp
        lifdip_d    = 9.30039e0_wp
        lifpol      = 72.9e0_wp
        lifie       = 11.3e0_wp
        hie         = 13.598e0_wp
        hpol        = 4.4997e0_wp
        rhf0did     = 6.e0_wp
        rhf0disp    = 6.e0_wp
        rhf0qq      = 6.e0_wp
        rhf0dq      = 6.e0_wp
        rlif0did    = 7.e0_wp
        rlif0disp   = 7.e0_wp
        r2lrcut_a   = 2.0e0_wp
        r2lrcut_r   = 3.e0_wp
        r3lrcut_a   = 2.e0_wp
        r3lrcut_r   = 4.5e0_wp
        cauang      = 0.52917706e0_wp
        caudeb      = 2.54177e0_wp
        mass        = [ 7.016003e0_wp, 1.00783e0_wp, 18.9984e0_wp]
        do i=1, nt
            ! PRECOMPUTE THINGS FOR USE LATER
            !     cosine of HF bond angle (i.e., Li-F-H angle)
            costhetahf = (r(i,3)**2-r(i,2)**2-r(i,1)**2)/(-2.e0_wp*r(i,1)*r(i,2))
            temp       = mass(3)/(mass(2)+mass(3))*r(i,2)
            !     Li-FH Jacobi distance
            bigrhf = sqrt(r(i,1)**2+temp**2-2.e0_wp*r(i,1)*temp*costhetahf)
            !     cosine of HF Jacobi angle (i.e., Li--HF angle)
            coschihf = -(r(i,1)**2-bigrhf**2-temp**2)/(-2.e0_wp*bigrhf*temp)
            coschihf = max(-1.e0_wp, coschihf)
            coschihf = min(1.e0_wp, coschihf)
            !     cosine of LiF bond angle (i.e., Li-F-H angle)
            costhetalif = (r(i,1)**2-r(i,2)**2-r(i,3)**2)/(-2.e0_wp*r(i,3)*r(i,2))
            temp        = mass(1)/(mass(1)+mass(3))*r(i,3)
            !     H-LiF Jacobi distance
            bigrlif = sqrt(r(i,2)**2+temp**2-2.e0_wp*r(i,2)*temp*costhetalif)
            !     cosine of LiF Jacobi angle (i.e., H--LiF angle)
            coschilif = -(r(i,2)**2-bigrlif**2-temp**2)/(-2.e0_wp*bigrlif*temp)
            coschilif = min(1.e0_wp, coschilif)
            coschilif = max(-1.e0_wp, coschilif)
            !     squared!
            r1s = r(i,1)**2
            r2s = r(i,2)**2
            r3s = r(i,3)**2
            !     Li-F-H bondangle
            cs = r3s+r2s-r1s
            cs = cs/(2.0e0_wp*r(i,2)*r(i,3))
            cs = min(cs, 1.0e0_wp)
            cs = max(cs, -1.0e0_wp)
            ! INTERACTION REGION
            !     The LiH singlet
            xelih = exp(-belihc*(r(i,1)-relih))
            s1    = s1c1*xelih**2  + s1c2*xelih
            s1c   = s1c3*exp(-bes1*(r(i,1)-relih))
            r1r3  = r(i,1) - r(i,3) + r1r3al*r(i,2)
            rcut2 = 0.5e0_wp*(1.0e0_wp - tanh((r1r3-rcut2rho)/rcut2del))
            s1    = s1 + rcut2 * s1c
            !     The LiH triplet
            xelih = exp(-bet1*(r(i,1)-reliht))
            t1    = t1c1*xelih**2  + t1c2*xelih
            !     The HF singlet
            yhf   = r(i,2) - 2.1042e0_wp
            behf  = 1.1622e0_wp -0.025647e0_wp*yhf+0.059062e0_wp*yhf**2
            xehfa = exp(-behf*(r(i,2)-rehf))
            s2    = dehf*xehfa*(xehfa-2.0e0_wp) + uli2p
            s2cut = (0.5e0_wp - 0.5e0_wp*tanh((r(i,2)-3.e0_wp)/0.5e0_wp))
            s2cc  = s2cc_d*(1.e0_wp-exp(-s2cc_b*(r(i,2)-s2cc_r)))**2 - s2cc_d
            s2    = s2cc + (s2 - s2cc)*s2cut
            dehfm = dehf - (dehfcc + dehfc*0.5e0_wp*(1.0e0_wp-cs))
            xehfm = exp(-behf*(r(i,2)-rehf))
            s2m   = dehfm*xehfm*(xehfm-2.0e0_wp) + uli2p
            s2cut = (0.5e0_wp - 0.5e0_wp*tanh((r(i,2)-2.2e0_wp)/0.5e0_wp))
            s2m   = s2m*s2cut
            s2    = s2m + (s2 - s2m)*0.5e0_wp*(1.e0_wp+tanh((r(i,3)-s2rho)/s2alp))
            !     The HF triplet
            xehfa     = exp(-bet2*(r(i,2)-rehft))
            xehfa_180 = exp(-bet2_180*(r(i,2)-rehft))
            t2        = (t2a1*xehfa**2 + t2a2*xehfa)*0.5e0_wp*(1.0e0_wp+cs) + (t2a3*xehfa_180**2 + t2a4*xehfa_180)* &
                & 0.5e0_wp* (1.0e0_wp- cs)
            !     The LiF singlet
            belif    = r(i,3)*103.57e0_wp/(4.6498e0_wp+7.0489e0_wp*r(i,3))**2+0.064076e0_wp
            gswt     = 0.5e0_wp-0.5e0_wp*tanh((r(i,2)-grho)/gdel)
            gswt     = gswt * 0.5e0_wp*(1.0e0_wp+cs)
            beliff   = belif + gswt*(belifff - belif)
            belifmod = beliff*(belif +(r(i,3)/gam)**nb)/ (beliff+(r(i,3)/gam)**nb)
            xelif    = exp(-belifmod*(r(i,3)-relif))
            s3a      = (delif+uli2p)*xelif*(xelif-2.0e0_wp)+uli2p
            s3cut    = 0.5e0_wp - 0.5e0_wp*tanh((r(i,3)-13.e0_wp)/0.5e0_wp)
            s3a      = s3a*s3cut
            swit     = swit_d  + (1.e0_wp - swit_d)*0.5e0_wp *(1.e0_wp + tanh((r(i,2)-swit_r)/swit_a))
            s3a      = s3a*swit
            delifm   = delifc1 * 0.5e0_wp*(1.0e0_wp-cs) + delifc2 * 0.5e0_wp*(1.0e0_wp+cs)
            belifm   = belifc1 * 0.5e0_wp*(1.0e0_wp-cs) + belifc2 * 0.5e0_wp*(1.0e0_wp+cs)
            xelif    = exp(-belifm*(r(i,3)-relif))
            s3c      = (delifm+uli2p)*xelif*(xelif-2.0e0_wp)
            rho2     = rho2_0 + (rho2_180-rho2_0)*0.5e0_wp*(1.0e0_wp-cs)
            theta2   = theta2_0 + (theta2_180-theta2_0)*0.5e0_wp*(1.0e0_wp-cs)
            del2     = del2_0 + (del2_180-del2_0)*0.5e0_wp*(1.0e0_wp-cs)
            gc2      = cos(theta2)*(rehf-r(i,2)) - sin(theta2)*(rho2-r(i,3))
            cut2     = 0.5e0_wp - 0.5e0_wp*tanh(gc2/del2)
            s3       = s3c + (s3a-s3c)*cut2
            !     The LiF triplet
            xelif     = exp(-bet3*(r(i,3)-relif))
            xelif_180 = exp(-bet3_180*(r(i,3)-relif))
            t3        = (t3a1*xelif**2 + t3a2*xelif)*0.5e0_wp*(1.0e0_wp+cs) + (t3a3*xelif_180**2 + t3a4*xelif_180)* &
                & 0.5e0_wp* (1.0e0_wp- cs)
            !     LEPS
            coul1 = 0.5e0_wp*(s1+t1)
            coul2 = 0.5e0_wp*(s2+t2)
            coul3 = 0.5e0_wp*(s3+t3)
            exch1 = 0.5e0_wp*(s1-t1)
            exch2 = 0.5e0_wp*(s2-t2)
            exch3 = 0.5e0_wp*(s3-t3)
            w     = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
            cplg2 = c2a*exp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
            e(i)  = coul1+coul2+coul3-sqrt(w+cplg2*cplg2)/sqrt(2.0e0_wp) + dehf
            ! LONG-RANGE FORCES * * * * * * * * * * * * * * * * * * * * * * * * * *
            ! HF DIPOLE
            hfdip_y = 1.e0_wp - exp(-hfdip_a*(r(i,2)*cauang-hfdip_re))
            hfdip   = 0.e0_wp
            do j=1, 11
                hfdip = hfdip + hfdip_m(j)*hfdip_y**real(j-1, wp)
            enddo
            hfdip = hfdip*0.5e0_wp*(1.e0_wp + tanh(20.e0_wp*(r(i,2)-1.2e0_wp)))
            ! HF POLARIZATION
            hfpol    = hfpolperp*(1.e0_wp - coschihf**2) + hfpolpar*coschihf**2
            hfpolavg = 2.e0_wp*hfpolperp/3.e0_wp + hfpolpar/3.e0_wp
            ! HF QUADRUPOLE
            hfqua_y = 1.e0_wp - exp(-hfqua_a*(r(i,2)*cauang-hfqua_re))
            hfqua   = 0.e0_wp
            do j=1, 10
                hfqua = hfqua + hfqua_m(j)*hfqua_y**real(j-1, wp)
            enddo
            hfqua = hfqua*0.5e0_wp*(1.e0_wp + tanh(20.e0_wp*(r(i,2)-1.2e0_wp)))
            ! LiF DIPOLE
            lifdip_y    = r(i,3) - lifdip_4
            lifdip_beta = lifdip_1 + lifdip_2*lifdip_y**8 + lifdip_3*lifdip_y
            lifdip      = lifdip_d*exp(-lifdip_beta*(r(i,3) - lifdip_re)**2)
            ! Li(2p) POLARIZATION
            li2ppol    = li2ppolperp*(1.e0_wp-coschihf**2)+li2ppolpar*coschihf**2
            li2ppolavg = 2.e0_wp*li2ppolperp/3.e0_wp + li2ppolpar/3.e0_wp
            ! HF + Li(2p)
            !     DIPOLE-INDUCED-DIPOLE
            hf_did = -0.5e0_wp*(3.e0_wp*coschihf**2 + 1.e0_wp)*li2ppol*hfdip**2
            hf_did = hf_did/(bigrhf**7 + rhf0did**7)/cevau*bigrhf
            !     cut off angular dependence in interaction region
            hf_did = hf_did*(1.e0_wp + (2.5e0_wp*li2ppolavg/ (3.e0_wp*coschihf**2 + 1.e0_wp)/li2ppol - 1.e0_wp) *exp(- &
                & (bigrhf/ rhf0did)** 4))
            !     DISPERSION
            hf_disp = -1.5e0_wp*li2pie*hfie*li2ppol*hfpol/(li2pie+hfie)
            hf_disp = hf_disp/(bigrhf**7 + rhf0disp**7)*bigrhf
            !     cut off angular dependence in interaction region
            hf_disp = hf_disp*(1.e0_wp+(li2ppolavg*hfpolavg/hfpol/li2ppol-1.e0_wp) *exp(-(bigrhf/rhf0disp)**4))
            !     QUADRUPOLE-QUADRUPOLE
            hf_qq = 3.e0_wp/4.e0_wp*hfqua*li2pqua*(3.e0_wp - 7.e0_wp*coschihf**2)
            hf_qq = hf_qq/(bigrhf**6 + rhf0qq**6)/cevau*bigrhf
            !     cut off angular dependence in interaction region
            hf_qq = hf_qq*(1.e0_wp + (-.5e0_wp/(3.e0_wp-7.e0_wp*coschihf**2) - 1.e0_wp) *exp(-(bigrhf/rhf0qq)**4))
            !     DIPOLE-QUADRUPOLE
            hf_dq = 3.e0_wp/4.e0_wp*hfdip*li2pqua*coschihf
            hf_dq = hf_dq/(bigrhf**5 + rhf0dq**5)/cevau*bigrhf
            !     cut off angular dependence in interaction region
            hf_qq = hf_qq*exp(-(bigrhf/rhf0dq)**4)
            !     SUM FORCES
            e_lr_hf = hf_dq + hf_qq + hf_disp + hf_did
            !     cut off at large HF distances
            e_lr_hf = e_lr_hf*0.5e0_wp*(1.e0_wp-tanh(r2lrcut_a*(r(i,2)-r2lrcut_r)))
            ! LiF + H
            !     DIPOLE-INDUCED-DIPOLE
            lif_did = -0.5e0_wp*(3.e0_wp*coschilif**2+1.e0_wp)*hpol*lifdip**2
            lif_did = lif_did/(bigrlif**6 + rlif0did**6)/cevau
            !     cut off angular dependence in interaction region
            lif_did = lif_did*(1.e0_wp+(2.5e0_wp/(3.e0_wp*coschilif**2+1.e0_wp)-1.e0_wp) *exp(-(bigrlif/rlif0did)**4))
            !     DISPERSION
            lif_disp = -1.5e0_wp*lifie*hie*lifpol*hpol/(lifie+hie)
            lif_disp = lif_disp/(bigrlif**6 + rlif0disp**6)
            !     no angular dependence
            !     SUM FORCES
            e_lr_lif = lif_disp + lif_did
            !     cut off at large LiF distances
            e_lr_lif = e_lr_lif*0.5e0_wp*(1.e0_wp-tanh(r3lrcut_a*(r(i,3)-r3lrcut_r)))
            ! LiH + F
            ! no long-range forces included
            ! CUT OFF LONG-RANGE FORCES
            if( bigrhf .gt. 15.e0_wp )then
                e_lr_hf = 0.e0_wp
            else
                e_lr_hf = e_lr_hf * exp( .2e0_wp / (bigrhf - 15.e0_wp) )
            endif
            if( bigrlif .gt. 15.e0_wp )then
                e_lr_lif = 0.e0_wp
            else
                e_lr_lif = e_lr_lif * exp( .2e0_wp / (bigrlif - 15.e0_wp) )
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
    use chempotpy_kinds,                                                 only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot22
