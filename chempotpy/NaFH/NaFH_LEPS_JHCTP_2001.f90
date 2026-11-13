!   System:                        NaFH
!   Common name:                   NaFH surface fit B
!   Functional form:               Modified Extended LEPS with a LEPS correction funct
!   Number of derivatives:         1
!   Number of electronic surfaces: 2
!   Interface:                     3-2V
!   References:  A. W. Jasper, M. D. Hack, A. Chakraborty, D. G. Truhlar
!                and P. Piecuch, J. Chem. Phys. 115, 7945 (2001).
!   Notes:  This is an improved version of NaFH surface fit A.  The diab
!           coupling has been made to vanish asymptotically and the diat
!           curves have been replaced with experimental curves.
!   Protocol:
!               must be called once before any calls to POT
!      POT    - driver for the evaluation of the energy
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      ground and first excited electronic states, diabatic representati
!  This is version 2 of the NaFH potential energy matrix.
!  Details can be found in
!  "Coupled Quasidiabatic Potential Energy Surfaces for LiFH",
!  M. D. Hack, A. W. Jasper, D. G. Truhlar, P. Piecuch, prepared
!  for publication in J. Chem. Phys.
!  This surface is distributed with the NAT code.
!  This surface uses version 1, but the diagonal diabats are corrected w
!  a correction function which is the difference of two LEPS functions.
!  surfaces are corrected such that they have the experimental bond ener
!  The diabatic coupling is cut off in the HF bond direction so that it
!  asymptotically.

!> @brief Implement NaFH_LEPS_JHCTP_2001, a potential energy surface for NaFH.
module chempotpy_nafh_leps_jhctp_2001_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NaFH_LEPS_JHCTP_2001 potential energy surface for NaFH.
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
        integer               :: iatom, idir
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
end module chempotpy_nafh_leps_jhctp_2001_pes

!> @brief Expose NaFH_LEPS_JHCTP_2001 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_nafh_leps_jhctp_2001_pes, only: evaluate
    use chempotpy_kinds,                    only: wp
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

!> @brief Evaluate EvdRdX for the NaFH_LEPS_JHCTP_2001 potential.
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

!> @brief Evaluate diagonalize for the NaFH_LEPS_JHCTP_2001 potential.
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

!> @brief Evaluate pot for the NaFH_LEPS_JHCTP_2001 potential.
subroutine pot( r, e, de, nt, nsurf )
    !      subroutine pot(r,e,nt,nsurf)
    !         R     is an array with dimensions R(NT,3)
    !               this array contains the geometry(ies) to be computed  (i
    !               Geometries are specified in terms of interatomic distanc
    !               that
    !               R(i,1) should be the distance between atoms 1 and 2 (A-B
    !               R(i,2) should be the distance between atoms 2 and 3 (B-C
    !               R(i,3) should be the distance between atoms 1 and 3 (A-C
    !               where i is an index in the range 1 <= i <= NT.
    !         E     is an array with dimensions E(NT)
    !               this array contains the diabatic potential energy(ies)
    !         DE    is an array of the derivatives of the diabatic potential
    !               energy(ies) with dimensions DE(3,NT) where the second
    !               index labels the derivative with respect to the specific
    !               component of R (interatomic distance)  (output)
    !         NT    is the number of geometries to be computed  (input)
    !         NSURF indicates the surface for which the energy should be
    !               computed  (input)
    !               The parameter NSURF may have the values 1, 2, or 3.
    !               These values select a particular diabatic surface:
    !               NSURF = 1  diabatic surface which is the lower one
    !                          in the reactant channel
    !               NSURF = 2  diabatic coupling surface
    !               NSURF = 3  diabatic surface which is the upper one
    !                          in the reactant channel
    use chempotpy_kinds, only: wp
    implicit none
    real(wp) :: de(3,nt), e(nt), r(nt,3)
    integer  :: nsurf, nt
    if( nsurf .eq. 1 )then
        call pot11(r, e, de, nt)
    elseif( nsurf .eq. 2 )then
        call pot12(r, e, de, nt)
    else
        call pot22(r, e, de, nt)
    endif
end subroutine pot

!  U11 surface - lower diabatic surface
module chempotpy_nafh_leps_jhctp_2001_pot11_legacy_state_3ad9bc64
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, de, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r, e, de, nt
        real(wp) :: abs1, ac(4), acut1(4), adip1(4), ads1, af1, af2, al2(4), al2n(4), aln0(4), aln0n1(4), aln0n2(4), &
            & amf, amh, amna, amu, anorm, ars1, asmall, axs1, bc(0:5), c2a(4), c2b(4), c2c(4), cevau, cfu11, chf(6), cnf(6), &
            & cnf1(5), cog, coul1, coul2, coul3, cplg2, csb, danormd1, danormd2, danormd3, dcd1, dcd2, dcd3, dcfu11(3), dcl1d1, &
            & dcl2d2, dcl3d3, dcog, dcogd1, dcogd2, dcogd3, dcsbd1, dcsbd2, dcsbd3, de(:,:), dehf, denf, dex1d1, dex2d2, dex3d3, &
            & dip1(4), dipn(4), dipnn(4), drnag2d1, drnag2d2, drnag2d3, ds1d1, ds2d2, ds3d3, dsurfd1, dsurfd2, dsurfd3, dt1d1, &
            & dt2d2, dt3d3, dwd1, dwd2, dwd3, e(:), exch1, exch2, exch3, g_r(3,1,3), hlp, hlp1, hlp11, hlp2, hlp22, hlp23, hlp3, &
            & r(:,:), r1s, r2s, r3s, rac2, rcut1(4), rdip1(4), rdum(3), rehf, renf, rn0n(4), rnag2, rootw, s1, s2, s3, surf, t1, &
            & t2, t3, una2p, w, y(4, 3, 4), y2, y3
        integer :: i, j, l, l1, l2, nt
        integer :: g_myrank, g_nprocs
        abs1       = 0.e0_wp
        ac         = 0.e0_wp
        acut1      = 0.e0_wp
        adip1      = 0.e0_wp
        ads1       = 0.e0_wp
        af1        = 0.e0_wp
        af2        = 0.e0_wp
        al2        = 0.e0_wp
        al2n       = 0.e0_wp
        aln0       = 0.e0_wp
        aln0n1     = 0.e0_wp
        aln0n2     = 0.e0_wp
        amf        = 0.e0_wp
        amh        = 0.e0_wp
        amna       = 0.e0_wp
        amu        = 0.e0_wp
        anorm      = 0.e0_wp
        ars1       = 0.e0_wp
        asmall     = 0.e0_wp
        axs1       = 0.e0_wp
        bc         = 0.e0_wp
        c2a        = 0.e0_wp
        c2b        = 0.e0_wp
        c2c        = 0.e0_wp
        cevau      = 0.e0_wp
        cfu11      = 0.e0_wp
        chf        = 0.e0_wp
        cnf        = 0.e0_wp
        cnf1       = 0.e0_wp
        cog        = 0.e0_wp
        coul1      = 0.e0_wp
        coul2      = 0.e0_wp
        coul3      = 0.e0_wp
        cplg2      = 0.e0_wp
        csb        = 0.e0_wp
        danormd1   = 0.e0_wp
        danormd2   = 0.e0_wp
        danormd3   = 0.e0_wp
        dcd1       = 0.e0_wp
        dcd2       = 0.e0_wp
        dcd3       = 0.e0_wp
        dcfu11     = 0.e0_wp
        dcl1d1     = 0.e0_wp
        dcl2d2     = 0.e0_wp
        dcl3d3     = 0.e0_wp
        dcog       = 0.e0_wp
        dcogd1     = 0.e0_wp
        dcogd2     = 0.e0_wp
        dcogd3     = 0.e0_wp
        dcsbd1     = 0.e0_wp
        dcsbd2     = 0.e0_wp
        dcsbd3     = 0.e0_wp
        dehf       = 0.e0_wp
        denf       = 0.e0_wp
        dex1d1     = 0.e0_wp
        dex2d2     = 0.e0_wp
        dex3d3     = 0.e0_wp
        dip1       = 0.e0_wp
        dipn       = 0.e0_wp
        dipnn      = 0.e0_wp
        drnag2d1   = 0.e0_wp
        drnag2d2   = 0.e0_wp
        drnag2d3   = 0.e0_wp
        ds1d1      = 0.e0_wp
        ds2d2      = 0.e0_wp
        ds3d3      = 0.e0_wp
        dsurfd1    = 0.e0_wp
        dsurfd2    = 0.e0_wp
        dsurfd3    = 0.e0_wp
        dt1d1      = 0.e0_wp
        dt2d2      = 0.e0_wp
        dt3d3      = 0.e0_wp
        dwd1       = 0.e0_wp
        dwd2       = 0.e0_wp
        dwd3       = 0.e0_wp
        exch1      = 0.e0_wp
        exch2      = 0.e0_wp
        exch3      = 0.e0_wp
        g_r        = 0.e0_wp
        hlp        = 0.e0_wp
        hlp1       = 0.e0_wp
        hlp11      = 0.e0_wp
        hlp2       = 0.e0_wp
        hlp22      = 0.e0_wp
        hlp23      = 0.e0_wp
        hlp3       = 0.e0_wp
        r1s        = 0.e0_wp
        r2s        = 0.e0_wp
        r3s        = 0.e0_wp
        rac2       = 0.e0_wp
        rcut1      = 0.e0_wp
        rdip1      = 0.e0_wp
        rdum       = 0.e0_wp
        rehf       = 0.e0_wp
        renf       = 0.e0_wp
        rn0n       = 0.e0_wp
        rnag2      = 0.e0_wp
        rootw      = 0.e0_wp
        s1         = 0.e0_wp
        s2         = 0.e0_wp
        s3         = 0.e0_wp
        surf       = 0.e0_wp
        t1         = 0.e0_wp
        t2         = 0.e0_wp
        t3         = 0.e0_wp
        una2p      = 0.e0_wp
        w          = 0.e0_wp
        y          = 0.e0_wp
        y2         = 0.e0_wp
        y3         = 0.e0_wp
        y(1, 1, 1) = 0.00000e0_wp
        y(2, 1, 1) = 0.00000e0_wp
        y(1, 1, 2) = 0.00000e0_wp
        y(2, 1, 2) = 0.00000e0_wp
        y(1, 1, 3) = 0.00000e0_wp
        y(2, 1, 3) = 0.00000e0_wp
        y(1, 1, 4) = 0.00000e0_wp
        y(2, 1, 4) = 0.00000e0_wp
        y(3, 1, 1) = 6.6796e0_wp
        y(4, 1, 1) = 2.7839e0_wp
        y(3, 1, 2) = 6.6796e0_wp
        y(4, 1, 2) = 2.7839e0_wp
        y(3, 1, 3) = 7.1360e0_wp
        y(4, 1, 3) = 1.0376e0_wp
        y(3, 1, 4) = 8.8737e0_wp
        y(4, 1, 4) = 1.9627e0_wp
        y(1, 2, 1) = 6.0459e+1_wp
        y(2, 2, 1) = 3.7503e0_wp
        y(1, 2, 2) = 6.0459e+1_wp
        y(2, 2, 2) = 3.7503e0_wp
        y(1, 2, 3) = 6.0459e+1_wp
        y(2, 2, 3) = 3.7503e0_wp
        y(1, 2, 4) = 8.5430e0_wp
        y(2, 2, 4) = 9.00856e-1_wp
        y(3, 2, 1) = 8.4874e0_wp
        y(4, 2, 1) = 1.0916e+1_wp
        y(3, 2, 2) = 8.4874e0_wp
        y(4, 2, 2) = 9.9850e0_wp
        y(3, 2, 3) = 8.4874e0_wp
        y(4, 2, 3) = 5.5898e0_wp
        y(3, 2, 4) = 9.5899e0_wp
        y(4, 2, 4) = 1.1615e+1_wp
        y(1, 3, 1) = 7.3251e0_wp
        y(2, 3, 1) = 2.8774e0_wp
        y(1, 3, 2) = 7.3251e0_wp
        y(2, 3, 2) = 3.6094e0_wp
        y(1, 3, 3) = 7.3251e0_wp
        y(2, 3, 3) = 3.9288e0_wp
        y(1, 3, 4) = 7.5811e0_wp
        y(2, 3, 4) = 4.4272e0_wp
        y(3, 3, 1) = 1.9300e+3_wp
        y(4, 3, 1) = 2.7091e0_wp
        y(3, 3, 2) = 3.8969e+3_wp
        y(4, 3, 2) = 2.7591e0_wp
        y(3, 3, 3) = 5.9130e+3_wp
        y(4, 3, 3) = 2.6881e0_wp
        y(3, 3, 4) = 5.9007e+3_wp
        y(4, 3, 4) = 3.8929e0_wp
        c2a(1)     = 8.5742e0_wp
        c2b(1)     = 1.5953e0_wp
        c2c(1)     = 3.3745e-1_wp
        c2a(2)     = 9.8509e0_wp
        c2b(2)     = 6.6359e0_wp
        c2c(2)     = 3.0396e-1_wp
        c2a(3)     = 9.5937e0_wp
        c2b(3)     = 1.3441e+1_wp
        c2c(3)     = 2.8812e-1_wp
        c2a(4)     = 1.1386e0_wp
        c2b(4)     = 8.1112e0_wp
        c2c(4)     = 2.6061e0_wp
        dip1(1)    = 1.9651e0_wp
        rdip1(1)   = 8.3109e-1_wp
        dip1(2)    = 1.9651e0_wp
        rdip1(2)   = 8.3109e-1_wp
        dip1(3)    = 2.1803e0_wp
        rdip1(3)   = 7.1647e-1_wp
        dip1(4)    = 2.2483e0_wp
        rdip1(4)   = 9.7128e-1_wp
        adip1(1)   = 4.8674e-2_wp
        adip1(2)   = 5.0603e-2_wp
        adip1(3)   = 4.8499e-2_wp
        adip1(4)   = 8.5924e-2_wp
        acut1(1)   = 6.3761e0_wp
        rcut1(1)   = 5.2007e0_wp
        acut1(2)   = 6.3761e0_wp
        rcut1(2)   = 3.6585e0_wp
        acut1(3)   = 6.0812e0_wp
        rcut1(3)   = 3.4684e0_wp
        acut1(4)   = 5.0516e0_wp
        rcut1(4)   = 1.6434e0_wp
        dipn(1)    = 1.5754e-1_wp
        al2(1)     = 8.4809e-2_wp
        dipn(2)    = 1.2274e-1_wp
        al2(2)     = 8.4809e-2_wp
        dipn(3)    = 4.1823e-2_wp
        al2(3)     = 9.4867e-2_wp
        dipn(4)    = 3.4350e-1_wp
        al2(4)     = 6.5542e-2_wp
        aln0(1)    = 3.1971e-2_wp
        aln0(2)    = 3.5817e-2_wp
        aln0(3)    = 1.6686e-1_wp
        aln0(4)    = 5.9961e-2_wp
        dipnn(1)   = 0.0000e0_wp
        al2n(1)    = 0.0000e-2_wp
        dipnn(2)   = 3.7501e-3_wp
        al2n(2)    = 9.2376e0_wp
        dipnn(3)   = 7.6323e-3_wp
        al2n(3)    = 9.3515e-2_wp
        dipnn(4)   = 0.0000e0_wp
        al2n(4)    = .0000e-2_wp
        rn0n(1)    = 0.0e0_wp
        rn0n(2)    = 6.6001e0_wp
        rn0n(3)    = 5.9949e0_wp
        rn0n(4)    = 0.0e0_wp
        aln0n1(1)  = 0.0000e-2_wp
        aln0n2(1)  = 0.0000e-2_wp
        aln0n1(2)  = 3.2545e-2_wp
        aln0n2(2)  = 4.9001e-1_wp
        aln0n1(3)  = 8.5914e-2_wp
        aln0n2(3)  = 2.1382e0_wp
        aln0n1(4)  = 0.0000e-2_wp
        aln0n2(4)  = 0.0000e-2_wp
        ac         = [ 4.6913e0_wp, 5.4684e-1_wp, 7.5932e-1_wp, 8.0000e-1_wp]
        amh        = 1.007825e0_wp
        amf        = 18.998403e0_wp
        amna       = 22.989767e0_wp
        amu        = 1822.8885e0_wp
        denf       = 4.3830525e0_wp
        renf       = 3.6395e0_wp
        cnf        = [ 0.168364e0_wp, 2.29333e0_wp, -0.848620e-1_wp, -0.109707e-1_wp, 0.246582e-2_wp, 0.377012e0_wp]
        cnf1       = [ 0.620115e0_wp, 0.694459e-1_wp, 0.757570e0_wp, 0.283055e0_wp, 2.69998e0_wp]
        dehf       = 5.77096e0_wp
        rehf       = 1.7328e0_wp
        chf        = [ 0.441258e0_wp, 3.88056e0_wp, -1.27110e0_wp, -1.37766e0_wp, 0.168186e0_wp, 2.07230e0_wp]
        bc         = [ -1.e0_wp, 0.e0_wp, 1.e0_wp, 2.e0_wp, 3.0e0_wp, 5.e0_wp]
        ads1       = 2.91087e-3_wp
        ars1       = 7.37707e0_wp
        abs1       = 7.13223e-1_wp
        una2p      = 2.0973375e0_wp
        asmall     = 0.09e0_wp
        cevau      = 0.036749308867692e0_wp
        g_myrank   = 0
        g_nprocs   = 1
        rac2       = 1.e0_wp/sqrt(2.e0_wp)
        af1        = amh/(amf+amh)
        af2        = amf/(amf+amh)
        if( g_myrank .eq. 0 )then
            !        write(6,100)af1,af2
            !     Excited NaF surface:
            !        write(6,102)(cNF1(i),i=1,5)
            !     Ground HF suface:
            !        write(6,103)DeHF,reHF,(cHF(i),i=1,6)
            !        write(6,104)ads1,ars1,abs1
            !        write(6,105)una2p,asmall
            !        write(6,106)(y(1,1,j),j=1,4)
            !        write(6,107)(y(2,1,j),j=1,4)
            !        write(6,108)(y(3,1,j),j=1,4)
            !        write(6,109)(y(4,1,j),j=1,4)
            !        write(6,110)(y(1,2,j),j=1,4)
            !        write(6,111)(y(2,2,j),j=1,4)
            !        write(6,112)(y(3,2,j),j=1,4)
            !        write(6,113)(y(4,2,j),j=1,4)
            !        write(6,114)(y(1,3,j),j=1,4)
            !        write(6,115)(y(2,3,j),j=1,4)
            !        write(6,116)(y(3,3,j),j=1,4)
            !        write(6,117)(y(4,3,j),j=1,4)
            !        write(6,118)(c2a(j),j=1,4)
            !        write(6,119)(c2b(j),j=1,4)
            !        write(6,120)(c2c(j),j=1,4)
            !        write(6,121)(dip1(j),j=1,4)
            !        write(6,122)(adip1(j),j=1,4)
            !        write(6,123)(rdip1(j),j=1,4)
            !        write(6,124)(acut1(j),j=1,4)
            !        write(6,125)(rcut1(j),j=1,4)
            !        write(6,126)(dipn(j),j=1,4)
            !        write(6,127)(al2(j),j=1,4)
            !        write(6,128)(aln0(j),j=1,4)
            !        write(6,129)(dipnn(j),j=1,4)
            !        write(6,130)(al2n(j),j=1,4)
            !        write(6,131)(aln0n1(j),j=1,4)
            !        write(6,132)(aln0n2(j),j=1,4)
            !        write(6,133)(rn0n(j),j=1,4)
        endif
        do i=1, nt
            r1s      = r(i,1)*r(i,1)
            r2s      = r(i,2)*r(i,2)
            r3s      = r(i,3)*r(i,3)
            rnag2    = af1*r1s+af2*r3s-af1*af2*r2s
            drnag2d1 = 2.e0_wp*af1*r(i,1)
            drnag2d2 = -2.e0_wp*af1*af2*r(i,2)
            drnag2d3 = 2.e0_wp*af2*r(i,3)
            !     Cosine of Na-F-H bond angle:
            hlp    = 1.e0_wp/(r(i,2)*r(i,3))
            csb    = 0.5e0_wp*(r2s+r3s-r1s)*hlp
            dcsbd1 = -2.e0_wp*r(i,1)*hlp
            dcsbd2 = 2.e0_wp*(1.e0_wp/r(i,3)-csb/r(i,2))
            dcsbd3 = 2.e0_wp*(1.e0_wp/r(i,2)-csb/r(i,3))
            csb    = 2.e0_wp*(1.e0_wp+csb)
            if( csb .ge. 4.e0_wp ) csb = 3.99999e0_wp
            l1 = idint(csb)+1
            l2 = l1+1
            if( l1 .eq. 4 ) l2=l1
            if( csb .le. bc(l1) ) csb = bc(l1)+1.e-13_wp
            if( csb .ge. bc(l1+1) ) csb = bc(l1+1)-1.e-13_wp
            axs1     = exp(-abs1*(r(i,1)-ars1))
            s1       = ads1*axs1*(axs1-2.e0_wp)
            ds1d1    = 2.e0_wp*ads1*abs1*axs1*(1.e0_wp-axs1)
            y2       = r(i,2)-rehf
            hlp1     = chf(1)*exp(-chf(2)*y2)
            hlp2     = exp(-chf(6)*y2)
            hlp22    = (-1.e0_wp-chf(1)+chf(3)*y2+chf(4)*y2**2+chf(5)*y2**3)*hlp2
            s2       = dehf*(hlp1+hlp22)
            ds2d2    = dehf*(-chf(2)*hlp1 - chf(6)*hlp22 + (chf(3)+2.e0_wp*chf(4)*y2+3.e0_wp*chf(5)*y2**2)*hlp2)
            y3       = r(i,3)-renf
            hlp1     = exp(-cnf1(3)*y3)
            hlp11    = (cnf1(1)+cnf1(2)*y3)*hlp1
            hlp2     = cnf1(4)*exp(-cnf1(5)*y3)
            s3       = hlp11+hlp2
            ds3d3    = -cnf1(3)*hlp11+cnf1(2)*hlp1-cnf1(5)*hlp2
            e(i)     = 0.e0_wp
            de(1,i)  = 0.e0_wp
            de(2,i)  = 0.e0_wp
            de(3,i)  = 0.e0_wp
            anorm    = 0.e0_wp
            danormd1 = 0.e0_wp
            danormd2 = 0.e0_wp
            danormd3 = 0.e0_wp
            do l = l1, l2
                hlp1   = y(1, 1, l)*exp(-y(2, 1, l)*r(i,1))
                hlp2   = y(3, 1, l)*exp(-y(4, 1, l)*r(i,1))
                t1     = hlp1+hlp2+s1
                dt1d1  = -y(2, 1, l)*hlp1-y(4, 1, l)*hlp2+ds1d1
                hlp1   = y(1, 2, l)*exp(-y(2, 2, l)*r(i,2))
                hlp2   = y(3, 2, l)*exp(-y(4, 2, l)*r(i,2))
                t2     = hlp1+hlp2+s2
                dt2d2  = -y(2, 2, l)*hlp1-y(4, 2, l)*hlp2+ds2d2
                hlp1   = y(1, 3, l)*exp(-y(2, 3, l)*r(i,3))
                hlp2   = y(3, 3, l)*exp(-y(4, 3, l)*r(i,3))
                t3     = hlp1+hlp2+s3
                dt3d3  = -y(2, 3, l)*hlp1-y(4, 3, l)*hlp2+ds3d3
                coul1  = 0.5e0_wp*(s1+t1)
                dcl1d1 = 0.5e0_wp*(ds1d1+dt1d1)
                coul2  = 0.5e0_wp*(s2+t2)
                dcl2d2 = 0.5e0_wp*(ds2d2+dt2d2)
                coul3  = 0.5e0_wp*(s3+t3)
                dcl3d3 = 0.5e0_wp*(ds3d3+dt3d3)
                exch1  = 0.5e0_wp*(s1-t1)
                dex1d1 = 0.5e0_wp*(ds1d1-dt1d1)
                exch2  = 0.5e0_wp*(s2-t2)
                dex2d2 = 0.5e0_wp*(ds2d2-dt2d2)
                exch3  = 0.5e0_wp*(s3-t3)
                dex3d3 = 0.5e0_wp*(ds3d3-dt3d3)
                w      = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
                dwd1   = 2.e0_wp*(2.e0_wp*exch1-exch2-exch3)*dex1d1
                dwd2   = 2.e0_wp*(2.e0_wp*exch2-exch1-exch3)*dex2d2
                dwd3   = 2.e0_wp*(2.e0_wp*exch3-exch1-exch2)*dex3d3
                cplg2  = c2a(l)*exp(-c2b(l)*w-c2c(l)*(r(i,1)+r(i,2)+r(i,3)))
                dcd1   = (-c2b(l)*dwd1 - c2c(l))*cplg2
                dcd2   = (-c2b(l)*dwd2 - c2c(l))*cplg2
                dcd3   = (-c2b(l)*dwd3 - c2c(l))*cplg2
                rootw  = rac2*sqrt(w+cplg2*cplg2)
                surf   = (coul1+coul2+coul3+dehf-rootw)
                ! check for zero by zero division (a case when all three atoms are asymp
                ! separated)
                if( (dwd1+2.e0_wp*cplg2*dcd1) .eq. 0.0e0_wp )then
                    dsurfd1 = dcl1d1
                else
                    dsurfd1 = (dcl1d1-0.25e0_wp*(dwd1+2.e0_wp*cplg2*dcd1)/rootw)
                endif
                if( (dwd2+2.e0_wp*cplg2*dcd2) .eq. 0.0e0_wp )then
                    dsurfd2 = dcl2d2
                else
                    dsurfd2 = (dcl2d2-0.25e0_wp*(dwd2+2.e0_wp*cplg2*dcd2)/rootw)
                endif
                if( (dwd3+2.e0_wp*cplg2*dcd3) .eq. 0.0e0_wp )then
                    dsurfd3 = dcl3d3
                else
                    dsurfd3 = (dcl3d3-0.25e0_wp*(dwd3+2.e0_wp*cplg2*dcd3)/rootw)
                endif
                hlp1    = adip1(l)*(rnag2-rdip1(l))
                hlp2    = acut1(l)*(r(i,2)-rcut1(l))
                hlp11   = 0.5e0_wp/(cosh(hlp1)**2)
                hlp22   = 0.5e0_wp/(cosh(hlp2)**2)
                hlp1    = 0.5e0_wp*(tanh(hlp1)-1.e0_wp)
                hlp2    = 0.5e0_wp*(1.e0_wp-tanh(hlp2))
                surf    = surf+dip1(l)*hlp1*hlp2
                dsurfd1 = dsurfd1+dip1(l)*adip1(l)*drnag2d1*hlp11*hlp2
                dsurfd2 = dsurfd2+dip1(l)*adip1(l)*drnag2d2*hlp11*hlp2 -dip1(l)*hlp1*acut1(l)*hlp22
                dsurfd3 = dsurfd3+dip1(l)*adip1(l)*drnag2d3*hlp11*hlp2
                hlp1    = al2(l)*(r(i,2)-rehf)
                hlp2    = aln0(l)*rnag2
                hlp11   = tanh(hlp1)
                hlp22   = tanh(hlp2)
                hlp1    = 1.e0_wp/cosh(hlp1)
                hlp2    = 1.e0_wp/cosh(hlp2)
                surf    = surf - dipn(l)*hlp1*hlp2
                dsurfd1 = dsurfd1 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*drnag2d1
                dsurfd2 = dsurfd2 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*drnag2d2 + dipn(l)*hlp2*hlp1*hlp11*al2(l)
                dsurfd3 = dsurfd3 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*drnag2d3
                hlp1    = al2n(l)*(r(i,2)-rehf)
                hlp2    = aln0n1(l)*(rnag2-rn0n(l)**2)
                ! try to prevent an overflow
                if( hlp2 .lt. 709.0e0_wp )then
                    hlp2 = exp(hlp2)
                else
                    hlp2 = 1.0e308_wp
                endif
                hlp3    = exp(-aln0n2(l)*(rnag2-rn0n(l)**2))
                hlp23   = 1.e0_wp/(hlp2+hlp3)**2
                hlp11   = tanh(hlp1)
                hlp1    = 1.e0_wp/cosh(hlp1)
                surf    = surf - 2.e0_wp*dipnn(l)*hlp1/(hlp2+hlp3)
                dsurfd1 = dsurfd1 + 2.e0_wp*dipnn(l)*hlp1 *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*drnag2d1
                dsurfd2=dsurfd2 + 2.e0_wp*dipnn(l)*hlp1 *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*drnag2d2 + 2.e0_wp* &
                    & dipnn(l)* hlp1* hlp11*al2n(l)/(hlp2+hlp3)
                dsurfd3  = dsurfd3 + 2.e0_wp*dipnn(l)*hlp1 *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*drnag2d3
                hlp      = 1.e0_wp/( (1.e0_wp-(csb-bc(l))/(bc(l-1)-bc(l))) * (1.e0_wp-(csb-bc(l))/(bc(l+1)-bc(l))) )
                hlp1     = (csb-bc(l))**2*hlp
                cog      = exp(-ac(l)*hlp1)
                dcog     = -cog*ac(l)*(2.e0_wp*(csb-bc(l))*hlp +hlp1*(1.e0_wp/(bc(l-1)-csb)+1.e0_wp/(bc(l+1)-csb)))
                dcogd1   = dcog*dcsbd1
                dcogd2   = dcog*dcsbd2
                dcogd3   = dcog*dcsbd3
                e(i)     = e(i) + surf*cog
                de(1,i)  = de(1,i) + dsurfd1*cog+surf*dcogd1
                de(2,i)  = de(2,i) + dsurfd2*cog+surf*dcogd2
                de(3,i)  = de(3,i) + dsurfd3*cog+surf*dcogd3
                anorm    = anorm + cog
                danormd1 = danormd1 + dcogd1
                danormd2 = danormd2 + dcogd2
                danormd3 = danormd3 + dcogd3
            enddo
            !     if (csb.gt.3.99d0) write(6,77)l1,l2,csb,r(i,3),
            !    &                         r(i,2),e(i),anorm
            hlp     = 1.e0_wp/anorm
            e(i)    = e(i)*hlp
            de(1,i) = (de(1,i) - e(i)*danormd1)*hlp
            de(2,i) = (de(2,i) - e(i)*danormd2)*hlp
            de(3,i) = (de(3,i) - e(i)*danormd3)*hlp
            !  77 format(2(2x,i3),5(2x,e16.8))
            do j=1, 3
                rdum(j) = r(i,j)
            enddo
            call lepscfu11(3, rdum, g_r, 3, cfu11, dcfu11, 3)
            e(i)    = e(i) + cfu11
            de(1,i) = de(1,i) + dcfu11(1)
            de(2,i) = de(2,i) + dcfu11(2)
            de(3,i) = de(3,i) + dcfu11(3)
            e(i)    = e(i)*cevau
            de(1,i) = de(1,i)*cevau
            de(2,i) = de(2,i)*cevau
            de(3,i) = de(3,i)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_nafh_leps_jhctp_2001_pot11_legacy_state_3ad9bc64

!> @brief Compatibility interface for pot11.
subroutine pot11( r, e, de, nt )
    use chempotpy_nafh_leps_jhctp_2001_pot11_legacy_state_3ad9bc64, only: dispatch
    use chempotpy_kinds,                                            only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    real(wp) :: de(3,nt)
    integer  :: nt
    call dispatch(r=r, e=e, de=de, nt=nt)
end subroutine pot11

!  U12 surface - coupling surface
module chempotpy_nafh_leps_jhctp_2001_pot12_legacy_state_3ad9bc64
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, de, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r, e, de, nt
        real(wp) :: abs1, ads1, af1, af2, ahf(3), ahft, amf, amh, amna, amu, anaf(3), anaft, anah(3), anaht, ars1, &
            & asmall, bc(0:5), ccnf(4), ccnh(4), ccnha(4), ccnhi(4), ccnhi1(4), cevau, cf_p, cf_r0, chf(6), cnf(6), cnf1(5), &
            & cstsq, dahftd1, dahftd2, dahftd3, danaftd1, danaftd2, danaftd3, danahtd1, danahtd2, danahtd3, dccnhd2(4), &
            & dccnhd3(4), dcstsqd1, dcstsqd2, dcstsqd3, de(:,:), dehf, denf, dr13td1, dr13td2, dr13td3, dr21td1, dr21td2, &
            & dr21td3, dr23td1, dr23td2, dr23td3, drnag2d1, drnag2d2, drnag2d3, dsw1d1, dsw1d2, dsw1d3, dsw21d1, dsw21d2, &
            & dsw21d3, dsw23d1, dsw23d2, dsw23d3, dswit1, dswit2, dswit3, dunafd3, dunahd1, dunahd2, dunahd3, e(:), hlp, hlp1, &
            & hlp2, r(:,:), r13(3), r13t, r1s, r21(3), r21t, r23(3), r23t, r2s, r3s, rehf, renf, rnag2, sw1, sw21, sw23, swit, &
            & una2p, unaf, unah, xpnh, xpnh1
        integer :: i, nt
        integer :: g_myrank, g_nprocs
        abs1     = 0.e0_wp
        ads1     = 0.e0_wp
        af1      = 0.e0_wp
        af2      = 0.e0_wp
        ahf      = 0.e0_wp
        ahft     = 0.e0_wp
        amf      = 0.e0_wp
        amh      = 0.e0_wp
        amna     = 0.e0_wp
        amu      = 0.e0_wp
        anaf     = 0.e0_wp
        anaft    = 0.e0_wp
        anah     = 0.e0_wp
        anaht    = 0.e0_wp
        ars1     = 0.e0_wp
        asmall   = 0.e0_wp
        bc       = 0.e0_wp
        ccnf     = 0.e0_wp
        ccnh     = 0.e0_wp
        ccnha    = 0.e0_wp
        ccnhi    = 0.e0_wp
        ccnhi1   = 0.e0_wp
        cevau    = 0.e0_wp
        cf_p     = 0.e0_wp
        cf_r0    = 0.e0_wp
        chf      = 0.e0_wp
        cnf      = 0.e0_wp
        cnf1     = 0.e0_wp
        cstsq    = 0.e0_wp
        dahftd1  = 0.e0_wp
        dahftd2  = 0.e0_wp
        dahftd3  = 0.e0_wp
        danaftd1 = 0.e0_wp
        danaftd2 = 0.e0_wp
        danaftd3 = 0.e0_wp
        danahtd1 = 0.e0_wp
        danahtd2 = 0.e0_wp
        danahtd3 = 0.e0_wp
        dccnhd2  = 0.e0_wp
        dccnhd3  = 0.e0_wp
        dcstsqd1 = 0.e0_wp
        dcstsqd2 = 0.e0_wp
        dcstsqd3 = 0.e0_wp
        dehf     = 0.e0_wp
        denf     = 0.e0_wp
        dr13td1  = 0.e0_wp
        dr13td2  = 0.e0_wp
        dr13td3  = 0.e0_wp
        dr21td1  = 0.e0_wp
        dr21td2  = 0.e0_wp
        dr21td3  = 0.e0_wp
        dr23td1  = 0.e0_wp
        dr23td2  = 0.e0_wp
        dr23td3  = 0.e0_wp
        drnag2d1 = 0.e0_wp
        drnag2d2 = 0.e0_wp
        drnag2d3 = 0.e0_wp
        dsw1d1   = 0.e0_wp
        dsw1d2   = 0.e0_wp
        dsw1d3   = 0.e0_wp
        dsw21d1  = 0.e0_wp
        dsw21d2  = 0.e0_wp
        dsw21d3  = 0.e0_wp
        dsw23d1  = 0.e0_wp
        dsw23d2  = 0.e0_wp
        dsw23d3  = 0.e0_wp
        dswit1   = 0.e0_wp
        dswit2   = 0.e0_wp
        dswit3   = 0.e0_wp
        dunafd3  = 0.e0_wp
        dunahd1  = 0.e0_wp
        dunahd2  = 0.e0_wp
        dunahd3  = 0.e0_wp
        hlp      = 0.e0_wp
        hlp1     = 0.e0_wp
        hlp2     = 0.e0_wp
        r13      = 0.e0_wp
        r13t     = 0.e0_wp
        r1s      = 0.e0_wp
        r21      = 0.e0_wp
        r21t     = 0.e0_wp
        r23      = 0.e0_wp
        r23t     = 0.e0_wp
        r2s      = 0.e0_wp
        r3s      = 0.e0_wp
        rehf     = 0.e0_wp
        renf     = 0.e0_wp
        rnag2    = 0.e0_wp
        sw1      = 0.e0_wp
        sw21     = 0.e0_wp
        sw23     = 0.e0_wp
        swit     = 0.e0_wp
        una2p    = 0.e0_wp
        unaf     = 0.e0_wp
        unah     = 0.e0_wp
        xpnh     = 0.e0_wp
        xpnh1    = 0.e0_wp
        amh      = 1.007825e0_wp
        amf      = 18.998403e0_wp
        amna     = 22.989767e0_wp
        amu      = 1822.8885e0_wp
        denf     = 4.3830525e0_wp
        renf     = 3.6395e0_wp
        cnf      = [ 0.168364e0_wp, 2.29333e0_wp, -0.848620e-1_wp, -0.109707e-1_wp, 0.246582e-2_wp, 0.377012e0_wp]
        cnf1     = [ 0.620115e0_wp, 0.694459e-1_wp, 0.757570e0_wp, 0.283055e0_wp, 2.69998e0_wp]
        dehf     = 5.77096e0_wp
        rehf     = 1.7328e0_wp
        chf      = [ 0.441258e0_wp, 3.88056e0_wp, -1.27110e0_wp, -1.37766e0_wp, 0.168186e0_wp, 2.07230e0_wp]
        bc       = [ -1.e0_wp, 0.e0_wp, 1.e0_wp, 2.e0_wp, 3.0e0_wp, 5.e0_wp]
        ads1     = 2.91087e-3_wp
        ars1     = 7.37707e0_wp
        abs1     = 7.13223e-1_wp
        una2p    = 2.0973375e0_wp
        asmall   = 0.09e0_wp
        cevau    = 0.036749308867692e0_wp
        xpnh     = 5.6906e-2_wp
        xpnh1    = 2.3185e-1_wp
        ahf      = [ 1.7177e0_wp, 2.2707e0_wp, 1.2248e-1_wp]
        anaf     = [ 6.7786e-1_wp, 9.8809e-1_wp, -1.0541e-1_wp]
        anah     = [ 6.5680e-1_wp, 1.0141e0_wp, -5.1056e-2_wp]
        r21      = [ 3.5114e0_wp, 1.3540e0_wp, 3.2267e-1_wp]
        r23      = [ -3.3006e-1_wp, 2.0804e-1_wp, 2.4499e-1_wp]
        r13      = [ 1.2431e0_wp, 9.9188e-1_wp, 1.1051e-1_wp]
        ccnf     = [ 3.9960e0_wp, 5.4730e-01_wp, -1.1428e0_wp, 2.0440e-01_wp]
        ccnha    = [ 1.00e0_wp, 0.80e0_wp, -2.67e0_wp, 0.456e0_wp]
        ccnhi    = [ 1.0897e0_wp, -6.5205e-1_wp, 7.7968e-1_wp, -1.2933e-2_wp]
        ccnhi1   = [ -6.1073e-1_wp, -1.7425e-1_wp, -8.3877e-1_wp, -3.5361e-1_wp]
        g_myrank = 0
        g_nprocs = 1
        af1      = amh/(amf+amh)
        af2      = amf/(amf+amh)
        if( g_myrank .eq. 0 )then
            !         write(6,100)af1,af2
            !     NaF coupling:
            !         write(6,101)(ccNF(i),i=1,4)
            !         write(6,102)(ccNHa(i),i=1,4)
            !         write(6,103)(ccNHi(i),i=1,4)
            !         write(6,104)(ccNHi1(i),i=1,4)
            !         write(6,105)xpNH,xpNH1
            !         write(6,106)(ahf(i),i=1,3)
            !         write(6,107)(anaf(i),i=1,3)
            !         write(6,108)(anah(i),i=1,3)
            !j         write(6,109)(r21(i),i=1,3)
            !         write(6,110)(r23(i),i=1,3)
            !         write(6,111)(r13(i),i=1,3)
        endif
        do i=1, nt
            hlp1       = exp(-xpnh*r(i,2)**2)
            hlp2       = exp(-xpnh1*r(i,3)**2)
            ccnh(1)    = ccnha(1)+ccnhi(1)*hlp1+ccnhi1(1)*hlp2
            ccnh(2)    = ccnha(2)+ccnhi(2)*hlp1+ccnhi1(2)*hlp2
            ccnh(3)    = ccnha(3)+ccnhi(3)*hlp1+ccnhi1(3)*hlp2
            ccnh(4)    = ccnha(4)+ccnhi(4)*hlp1+ccnhi1(4)*hlp2
            hlp1       = - 2.e0_wp*xpnh*r(i,2)*hlp1
            hlp2       = - 2.e0_wp*xpnh1*r(i,3)*hlp2
            dccnhd2(1) = ccnhi(1)*hlp1
            dccnhd3(1) = ccnhi1(1)*hlp2
            dccnhd2(2) = ccnhi(2)*hlp1
            dccnhd3(2) = ccnhi1(2)*hlp2
            dccnhd2(3) = ccnhi(3)*hlp1
            dccnhd3(3) = ccnhi1(3)*hlp2
            dccnhd2(4) = ccnhi(4)*hlp1
            dccnhd3(4) = ccnhi1(4)*hlp2
            hlp        = ccnf(3)/(1.e0_wp+ccnf(4)*r(i,3))
            unaf       = exp(-ccnf(1)-ccnf(2)*r(i,3)- hlp*r(i,3))
            unaf       = unaf*27.2113961e0_wp
            dunafd3    = -unaf*(ccnf(2)+hlp-ccnf(4)*hlp*hlp*r(i,3)/ccnf(3))
            hlp        = 1.e0_wp/(1.e0_wp+ccnh(4)*r(i,1))
            unah       = exp(-ccnh(1)-ccnh(2)*r(i,1)- ccnh(3)*hlp*r(i,1))
            dunahd1    = -unah*(ccnh(2)+ccnh(3)*hlp -ccnh(4)*hlp*hlp*ccnh(3)*r(i,1))
            dunahd2    = -unah*(dccnhd2(1)+dccnhd2(2)*r(i,1) + dccnhd2(3)*r(i,1)*hlp-ccnh(3)*r(i,1)*r(i,1)* &
                & dccnhd2(4)* hlp*hlp)
            dunahd3    = -unah*(dccnhd3(1)+dccnhd3(2)*r(i,1) + dccnhd3(3)*r(i,1)*hlp-ccnh(3)*r(i,1)*r(i,1)* &
                & dccnhd3(4)* hlp*hlp)
            r1s      = r(i,1)*r(i,1)
            r2s      = r(i,2)*r(i,2)
            r3s      = r(i,3)*r(i,3)
            rnag2    = af1*r1s+af2*r3s-af1*af2*r2s
            drnag2d1 = 2.e0_wp*af1*r(i,1)
            drnag2d2 = -2.e0_wp*af1*af2*r(i,2)
            drnag2d3 = 2.e0_wp*af2*r(i,3)
            hlp2     = rnag2+asmall
            hlp      = sqrt(hlp2)
            cstsq    = 0.5e0_wp*(r3s-r1s+(af2-af1)*r2s) /(r(i,2)*hlp)
            dcstsqd1 = -2.e0_wp*r(i,1)/(r(i,2)*hlp) - cstsq*drnag2d1/hlp2
            dcstsqd2 = 2.e0_wp*((af2-af1)/hlp-cstsq/r(i,2)) - cstsq*drnag2d2/hlp2
            dcstsqd3 = 2.e0_wp*r(i,3)/(r(i,2)*hlp) - cstsq*drnag2d3/hlp2
            cstsq    = 2.0e0_wp*(1.e0_wp+cstsq)
            hlp      = r21(2)-r21(1)-r21(3)
            r21t     = r21(1)+hlp*cstsq +r21(3)*cstsq*cstsq
            hlp1     = hlp+2.e0_wp*r21(3)*cstsq
            dr21td1  = hlp1*dcstsqd1
            dr21td2  = hlp1*dcstsqd2
            dr21td3  = hlp1*dcstsqd3
            hlp      = r23(2)-r23(1)-r23(3)
            r23t     = r23(1)+hlp*cstsq +r23(3)*cstsq*cstsq
            hlp1     = hlp+2.e0_wp*r23(3)*cstsq
            dr23td1  = hlp1*dcstsqd1
            dr23td2  = hlp1*dcstsqd2
            dr23td3  = hlp1*dcstsqd3
            hlp      = r13(2)-r13(1)-r13(3)
            r13t     = r13(1)+hlp*cstsq +r13(3)*cstsq*cstsq
            hlp1     = hlp+2.e0_wp*r13(3)*cstsq
            dr13td1  = hlp1*dcstsqd1
            dr13td2  = hlp1*dcstsqd2
            dr13td3  = hlp1*dcstsqd3
            hlp      = ahf(2)-ahf(1)-ahf(3)
            ahft     = ahf(1)+hlp*cstsq +ahf(3)*cstsq*cstsq
            hlp1     = hlp+2.e0_wp*ahf(3)*cstsq
            dahftd1  = hlp1*dcstsqd1
            dahftd2  = hlp1*dcstsqd2
            dahftd3  = hlp1*dcstsqd3
            hlp      = anaf(2)-anaf(1)-anaf(3)
            anaft    = anaf(1)+hlp*cstsq +anaf(3)*cstsq*cstsq
            hlp1     = hlp+2.e0_wp*anaf(3)*cstsq
            danaftd1 = hlp1*dcstsqd1
            danaftd2 = hlp1*dcstsqd2
            danaftd3 = hlp1*dcstsqd3
            hlp      = anah(2)-anah(1)-anah(3)
            anaht    = anah(1)+hlp*cstsq +anah(3)*cstsq*cstsq
            hlp1     = hlp+2.e0_wp*anah(3)*cstsq
            danahtd1 = hlp1*dcstsqd1
            danahtd2 = hlp1*dcstsqd2
            danahtd3 = hlp1*dcstsqd3
            hlp      = anaht*r(i,1)-anaft*r(i,3)-r13t
            sw1      = 0.5e0_wp*(1.e0_wp+tanh(hlp))
            hlp1     = 0.5e0_wp/(cosh(hlp)**2)
            dsw1d1   = hlp1*(anaht+danahtd1*r(i,1)-danaftd1*r(i,3)-dr13td1)
            dsw1d2   = hlp1*(danahtd2*r(i,1)-danaftd2*r(i,3)-dr13td2)
            dsw1d3   = hlp1*(danahtd3*r(i,1)-anaft-danaftd3*r(i,3)-dr13td3)
            hlp      = ahft*r(i,2)-anaht*r(i,1)-r21t
            sw21     = 0.5e0_wp*(1.e0_wp+tanh(hlp))
            hlp1     = 0.5e0_wp/(cosh(hlp)**2)
            dsw21d1  = hlp1*(dahftd1*r(i,2)-anaht-danahtd1*r(i,1)-dr21td1)
            dsw21d2  = hlp1*(ahft+dahftd2*r(i,2)-danahtd2*r(i,1)-dr21td2)
            dsw21d3  = hlp1*(dahftd3*r(i,2)-danahtd3*r(i,1)-dr21td3)
            hlp      = ahft*r(i,2)-anaft*r(i,3)-r23t
            sw23     = 0.5e0_wp*(1.e0_wp+tanh(hlp))
            hlp1     = 0.5e0_wp/(cosh(hlp)**2)
            dsw23d1  = hlp1*(dahftd1*r(i,2)-danaftd1*r(i,3)-dr23td1)
            dsw23d2  = hlp1*(ahft+dahftd2*r(i,2)-danaftd2*r(i,3)-dr23td2)
            dsw23d3  = hlp1*(dahftd3*r(i,2)-anaft-danaftd3*r(i,3)-dr23td3)
            e(i)     = unaf*sw1*sw23+unah*(1.e0_wp-sw1)*sw21
            de(1,i) = unaf*dsw1d1*sw23+unaf*sw1*dsw23d1 + dunahd1*(1.e0_wp-sw1)*sw21-unah*dsw1d1*sw21 + unah*(1.e0_wp- &
                & sw1)* dsw21d1
            de(2,i) = unaf*dsw1d2*sw23+unaf*sw1*dsw23d2 + dunahd2*(1.e0_wp-sw1)*sw21-unah*dsw1d2*sw21 + unah*(1.e0_wp- &
                & sw1)* dsw21d2
            de(3,i) = dunafd3*sw1*sw23+unaf*dsw1d3*sw23+unaf*sw1*dsw23d3 + dunahd3*(1.e0_wp-sw1)*sw21-unah*dsw1d3* &
                & sw21 + unah* (1.e0_wp-sw1)*dsw21d3
            !     correction function
            cf_p    = 2.e0_wp
            cf_r0   = 3.5e0_wp
            swit    = 0.5e0_wp*(1.e0_wp+tanh(cf_p*(cf_r0-r(i,2))))
            dswit1  = 0.e0_wp
            dswit2  = -0.5e0_wp*cf_p*1.e0_wp/(cosh(cf_p*(-r(i,2)+cf_r0))**2)
            dswit3  = 0.e0_wp
            de(1,i) = e(i)*dswit1+swit*de(1,i)
            de(2,i) = e(i)*dswit2+swit*de(2,i)
            de(3,i) = e(i)*dswit3+swit*de(3,i)
            e(i)    = e(i)*swit
            e(i)    = e(i)*cevau
            de(1,i) = de(1,i)*cevau
            de(2,i) = de(2,i)*cevau
            de(3,i) = de(3,i)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_nafh_leps_jhctp_2001_pot12_legacy_state_3ad9bc64

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, de, nt )
    use chempotpy_nafh_leps_jhctp_2001_pot12_legacy_state_3ad9bc64, only: dispatch
    use chempotpy_kinds,                                            only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    real(wp) :: de(3,nt)
    integer  :: nt
    call dispatch(r=r, e=e, de=de, nt=nt)
end subroutine pot12

!  U22 surface - upper diabatic surface
module chempotpy_nafh_leps_jhctp_2001_pot22_legacy_state_3ad9bc64
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, de, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r, e, de, nt
        real(wp) :: abs1, ac(4), acut2(4), acut2n(4), adip2(4), adip2n(4), ads1, af1, af2, al20, al20a, al20i(4), &
            & al20i1(4), al2a(4), aly(4), amf, amh, amna, amu, anorm, ars1, asmall, bc(0:5), c2a(4), c2b(4), c2c(4), cevau, &
            & cfu22, chf(6), cnf(6), cnf1(5), cog, coul1, coul2, coul3, cplg2, csb, ctanh, cy, dal20, danormd1, danormd2, &
            & danormd3, dcd1, dcd2, dcd3, dcfu22(3), dcl1d1, dcl2d1, dcl2d2, dcl2d3, dcl3d3, dcog, dcogd1, dcogd2, dcogd3, &
            & dcsbd1, dcsbd2, dcsbd3, dctanh, dcyd3, de(:,:), dehf, denf, dex1d1, dex2d1, dex2d2, dex2d3, dex3d3, dip2(4), &
            & dip2n(4), dr20, drnag2d1, drnag2d2, drnag2d3, ds1d1, ds2ad2, ds2d1, ds2d2, ds2d3, ds3d3, dsurfd1, dsurfd2, &
            & dsurfd3, dt1d1, dt2d1, dt2d2, dt2d3, dt3d3, dwd1, dwd2, dwd3, dy(4), dy1(4), dy2(4), e(:), exch1, exch2, exch3, &
            & g_r(3,1,3), hlp, hlp1, hlp11, hlp2, hlp22, hlp3, psi, r(:,:), r1s, r20, r20a, r20i(4), r2s, r3s, rac2, rcut2(4), &
            & rcut2n(4), rdip2(4), rdip2n(4), rdum(3), rehf, renf, rnag2, rootw, rr0a(4), rsi, ry(4), s1, s2, s2a, s3, surf, t1, &
            & t2, t3, una2p, w, y(4, 3, 4), y2, y3
        integer :: i, j, l, l1, l2, nt
        integer :: g_myrank, g_nprocs
        abs1       = 0.e0_wp
        ac         = 0.e0_wp
        acut2      = 0.e0_wp
        acut2n     = 0.e0_wp
        adip2      = 0.e0_wp
        adip2n     = 0.e0_wp
        ads1       = 0.e0_wp
        af1        = 0.e0_wp
        af2        = 0.e0_wp
        al20       = 0.e0_wp
        al20a      = 0.e0_wp
        al20i      = 0.e0_wp
        al20i1     = 0.e0_wp
        al2a       = 0.e0_wp
        aly        = 0.e0_wp
        amf        = 0.e0_wp
        amh        = 0.e0_wp
        amna       = 0.e0_wp
        amu        = 0.e0_wp
        anorm      = 0.e0_wp
        ars1       = 0.e0_wp
        asmall     = 0.e0_wp
        bc         = 0.e0_wp
        c2a        = 0.e0_wp
        c2b        = 0.e0_wp
        c2c        = 0.e0_wp
        cevau      = 0.e0_wp
        cfu22      = 0.e0_wp
        chf        = 0.e0_wp
        cnf        = 0.e0_wp
        cnf1       = 0.e0_wp
        cog        = 0.e0_wp
        coul1      = 0.e0_wp
        coul2      = 0.e0_wp
        coul3      = 0.e0_wp
        cplg2      = 0.e0_wp
        csb        = 0.e0_wp
        ctanh      = 0.e0_wp
        cy         = 0.e0_wp
        dal20      = 0.e0_wp
        danormd1   = 0.e0_wp
        danormd2   = 0.e0_wp
        danormd3   = 0.e0_wp
        dcd1       = 0.e0_wp
        dcd2       = 0.e0_wp
        dcd3       = 0.e0_wp
        dcfu22     = 0.e0_wp
        dcl1d1     = 0.e0_wp
        dcl2d1     = 0.e0_wp
        dcl2d2     = 0.e0_wp
        dcl2d3     = 0.e0_wp
        dcl3d3     = 0.e0_wp
        dcog       = 0.e0_wp
        dcogd1     = 0.e0_wp
        dcogd2     = 0.e0_wp
        dcogd3     = 0.e0_wp
        dcsbd1     = 0.e0_wp
        dcsbd2     = 0.e0_wp
        dcsbd3     = 0.e0_wp
        dctanh     = 0.e0_wp
        dcyd3      = 0.e0_wp
        dehf       = 0.e0_wp
        denf       = 0.e0_wp
        dex1d1     = 0.e0_wp
        dex2d1     = 0.e0_wp
        dex2d2     = 0.e0_wp
        dex2d3     = 0.e0_wp
        dex3d3     = 0.e0_wp
        dip2       = 0.e0_wp
        dip2n      = 0.e0_wp
        dr20       = 0.e0_wp
        drnag2d1   = 0.e0_wp
        drnag2d2   = 0.e0_wp
        drnag2d3   = 0.e0_wp
        ds1d1      = 0.e0_wp
        ds2ad2     = 0.e0_wp
        ds2d1      = 0.e0_wp
        ds2d2      = 0.e0_wp
        ds2d3      = 0.e0_wp
        ds3d3      = 0.e0_wp
        dsurfd1    = 0.e0_wp
        dsurfd2    = 0.e0_wp
        dsurfd3    = 0.e0_wp
        dt1d1      = 0.e0_wp
        dt2d1      = 0.e0_wp
        dt2d2      = 0.e0_wp
        dt2d3      = 0.e0_wp
        dt3d3      = 0.e0_wp
        dwd1       = 0.e0_wp
        dwd2       = 0.e0_wp
        dwd3       = 0.e0_wp
        dy         = 0.e0_wp
        dy1        = 0.e0_wp
        dy2        = 0.e0_wp
        exch1      = 0.e0_wp
        exch2      = 0.e0_wp
        exch3      = 0.e0_wp
        g_r        = 0.e0_wp
        hlp        = 0.e0_wp
        hlp1       = 0.e0_wp
        hlp11      = 0.e0_wp
        hlp2       = 0.e0_wp
        hlp22      = 0.e0_wp
        hlp3       = 0.e0_wp
        psi        = 0.e0_wp
        r1s        = 0.e0_wp
        r20        = 0.e0_wp
        r20a       = 0.e0_wp
        r20i       = 0.e0_wp
        r2s        = 0.e0_wp
        r3s        = 0.e0_wp
        rac2       = 0.e0_wp
        rcut2      = 0.e0_wp
        rcut2n     = 0.e0_wp
        rdip2      = 0.e0_wp
        rdip2n     = 0.e0_wp
        rdum       = 0.e0_wp
        rehf       = 0.e0_wp
        renf       = 0.e0_wp
        rnag2      = 0.e0_wp
        rootw      = 0.e0_wp
        rr0a       = 0.e0_wp
        rsi        = 0.e0_wp
        ry         = 0.e0_wp
        s1         = 0.e0_wp
        s2         = 0.e0_wp
        s2a        = 0.e0_wp
        s3         = 0.e0_wp
        surf       = 0.e0_wp
        t1         = 0.e0_wp
        t2         = 0.e0_wp
        t3         = 0.e0_wp
        una2p      = 0.e0_wp
        w          = 0.e0_wp
        y          = 0.e0_wp
        y2         = 0.e0_wp
        y3         = 0.e0_wp
        ac         = [ 2.0018e0_wp, 4.1908e-1_wp, 5.4988e0_wp, 3.0178e-1_wp]
        y(1, 1, 1) = 1.0407e+3_wp
        y(2, 1, 1) = 4.5911e0_wp
        y(1, 1, 2) = 1.1608e+3_wp
        y(2, 1, 2) = 3.8886e0_wp
        y(1, 1, 3) = 4.9492e+3_wp
        y(2, 1, 3) = 2.6147e0_wp
        y(1, 1, 4) = 8.3093e+3_wp
        y(2, 1, 4) = 3.2431e0_wp
        y(3, 1, 1) = 1.4941e0_wp
        y(4, 1, 1) = 4.5525e0_wp
        y(3, 1, 2) = 1.4941e0_wp
        y(4, 1, 2) = 4.5525e0_wp
        y(3, 1, 3) = 1.4941e0_wp
        y(4, 1, 3) = 4.5525e0_wp
        y(3, 1, 4) = 2.8361e0_wp
        y(4, 1, 4) = 4.4686e0_wp
        y(1, 2, 1) = 0.0000e0_wp
        y(2, 2, 1) = 0.0000e0_wp
        y(1, 2, 2) = 0.0000e0_wp
        y(2, 2, 2) = 0.0000e0_wp
        y(1, 2, 3) = 0.0000e0_wp
        y(2, 2, 3) = 0.0000e0_wp
        y(1, 2, 4) = 0.0000e0_wp
        y(2, 2, 4) = 0.0000e0_wp
        y(3, 2, 1) = 1.8076e+3_wp
        y(4, 2, 1) = 2.8815e0_wp
        y(3, 2, 2) = 1.8076e+3_wp
        y(4, 2, 2) = 2.8815e0_wp
        y(3, 2, 3) = 1.8076e+3_wp
        y(4, 2, 3) = 2.8815e0_wp
        y(3, 2, 4) = 2.3109e+3_wp
        y(4, 2, 4) = 3.1315e0_wp
        y(1, 3, 1) = 5.7330e0_wp
        y(2, 3, 1) = 3.2053e-1_wp
        y(1, 3, 2) = 5.7330e0_wp
        y(2, 3, 2) = 3.2053e-1_wp
        y(1, 3, 3) = 5.7330e0_wp
        y(2, 3, 3) = 3.2053e-1_wp
        y(1, 3, 4) = 5.0393e0_wp
        y(2, 3, 4) = 2.7567e-1_wp
        y(3, 3, 1) = 2.5774e+5_wp
        y(4, 3, 1) = 4.2267e0_wp
        y(3, 3, 2) = 2.5774e+5_wp
        y(4, 3, 2) = 4.2267e0_wp
        y(3, 3, 3) = 2.5774e+5_wp
        y(4, 3, 3) = 4.2267e0_wp
        y(3, 3, 4) = 1.6506e+5_wp
        y(4, 3, 4) = 3.9246e0_wp
        c2a(1)     = 5.6406e0_wp
        c2b(1)     = 7.7152e-2_wp
        c2a(2)     = 8.2390e0_wp
        c2b(2)     = 1.7133e-1_wp
        c2a(3)     = 7.1504e0_wp
        c2b(3)     = 1.1011e-1_wp
        c2a(4)     = 8.0196e0_wp
        c2b(4)     = 1.0017e-1_wp
        c2c(1)     = 8.4595e-2_wp
        c2c(2)     = 8.4595e-2_wp
        c2c(3)     = 8.4595e-2_wp
        c2c(4)     = 8.4411e-2_wp
        r20i(1)    = 2.9013e0_wp
        al20i(1)   = 1.7119e0_wp
        r20i(2)    = 3.1196e0_wp
        al20i(2)   = 9.4849e-1_wp
        r20i(3)    = 3.0995e0_wp
        al20i(3)   = 0.0e0_wp
        r20i(4)    = 2.6572e0_wp
        al20i(4)   = 2.3214e0_wp
        r20a       = 3.1305e0_wp
        al20i1(1)  = 3.5428e0_wp
        al20a      = 3.024e+1_wp
        al20i1(2)  = 6.7257e0_wp
        psi        = 0.7e0_wp
        al20i1(3)  = 3.0240e+1_wp
        rsi        = 22.2e0_wp
        al20i1(4)  = 9.6978e0_wp
        al2a(1)    = 5.8088e-1_wp
        rr0a(1)    = 9.5543e0_wp
        al2a(2)    = 5.8088e-1_wp
        rr0a(2)    = 9.5543e0_wp
        al2a(3)    = 5.8088e-1_wp
        rr0a(3)    = 9.5543e0_wp
        al2a(4)    = 8.1525e-1_wp
        rr0a(4)    = 8.2500e0_wp
        dip2(1)    = 2.6475e-1_wp
        dip2(2)    = 2.6475e-1_wp
        dip2(3)    = 2.6475e-1_wp
        dip2(4)    = 5.0336e-1_wp
        adip2(1)   = 3.2317e-1_wp
        rdip2(1)   = 3.7260e0_wp
        adip2(2)   = 1.8513e-1_wp
        rdip2(2)   = 4.5515e0_wp
        adip2(3)   = 4.3999e-2_wp
        rdip2(3)   = 2.1191e-1_wp
        adip2(4)   = 2.3003e-1_wp
        rdip2(4)   = 5.0160e0_wp
        acut2(1)   = 2.0799e0_wp
        rcut2(1)   = 1.2496e0_wp
        acut2(2)   = 1.4204e0_wp
        rcut2(2)   = 1.0817e0_wp
        acut2(3)   = 7.9086e-1_wp
        rcut2(3)   = 1.0575e0_wp
        acut2(4)   = 2.8447e0_wp
        rcut2(4)   = 3.0015e0_wp
        dip2n(1)   = 2.2091e-1_wp
        adip2n(1)  = 1.0288e-1_wp
        dip2n(2)   = 2.2091e-1_wp
        adip2n(2)  = 1.0288e-1_wp
        dip2n(3)   = 2.2091e-1_wp
        adip2n(3)  = 1.0288e-1_wp
        dip2n(4)   = 4.5301e-1_wp
        adip2n(4)  = 1.1009e-1_wp
        rdip2n(1)  = 5.1554e0_wp
        rdip2n(2)  = 5.1554e0_wp
        rdip2n(3)  = 5.1554e0_wp
        rdip2n(4)  = 6.7003e0_wp
        acut2n(1)  = 3.5607e0_wp
        rcut2n(1)  = 3.0661e0_wp
        acut2n(2)  = 4.3956e0_wp
        rcut2n(2)  = 3.0661e0_wp
        acut2n(3)  = 4.1924e0_wp
        rcut2n(3)  = 3.0661e0_wp
        acut2n(4)  = 3.0765e0_wp
        rcut2n(4)  = 2.8845e0_wp
        ry(1)      = 1.4328e+1_wp
        aly(1)     = 4.9044e-1_wp
        ry(2)      = 1.2504e+1_wp
        aly(2)     = 6.4718e-1_wp
        ry(3)      = 1.2102e+1_wp
        aly(3)     = 2.5833e-1_wp
        ry(4)      = 7.9381e0_wp
        aly(4)     = 6.4853e-1_wp
        dy(1)      = 7.5126e-1_wp
        dy1(1)     = 0.e0_wp
        dy2(1)     = 0.e0_wp
        dy(2)      = 2.0210e-1_wp
        dy1(2)     = 0.e0_wp
        dy2(2)     = 0.e0_wp
        dy(3)      = 1.1749e-1_wp
        dy1(3)     = 0.e0_wp
        dy2(3)     = 0.e0_wp
        dy(4)      = 8.1959e-2_wp
        dy1(4)     = -5.03e+2_wp
        dy2(4)     = -0.25e0_wp
        amh        = 1.007825e0_wp
        amf        = 18.998403e0_wp
        amna       = 22.989767e0_wp
        amu        = 1822.8885e0_wp
        denf       = 4.3830525e0_wp
        renf       = 3.6395e0_wp
        cnf        = [ 0.168364e0_wp, 2.29333e0_wp, -0.848620e-1_wp, -0.109707e-1_wp, 0.246582e-2_wp, 0.377012e0_wp]
        cnf1       = [ 0.620115e0_wp, 0.694459e-1_wp, 0.757570e0_wp, 0.283055e0_wp, 2.69998e0_wp]
        dehf       = 5.77096e0_wp
        rehf       = 1.7328e0_wp
        chf        = [ 0.441258e0_wp, 3.88056e0_wp, -1.27110e0_wp, -1.37766e0_wp, 0.168186e0_wp, 2.07230e0_wp]
        bc         = [ -1.e0_wp, 0.e0_wp, 1.e0_wp, 2.e0_wp, 3.0e0_wp, 5.e0_wp]
        ads1       = 2.91087e-3_wp
        ars1       = 7.37707e0_wp
        abs1       = 7.13223e-1_wp
        una2p      = 2.0973375e0_wp
        asmall     = 0.09e0_wp
        cevau      = 0.036749308867692e0_wp
        g_myrank   = 0
        g_nprocs   = 1
        rac2       = 1.e0_wp/sqrt(2.e0_wp)
        af1        = amh/(amf+amh)
        af2        = amf/(amf+amh)
        if( g_myrank .eq. 0 )then
            !         write(6,100)af1,af2
            !     Ground NaF surface
            !         write(6,101)DeNF,reNF,(cNF(i),i=1,6)
            !     Ground HF suface:
            !         write(6,103)DeHF,reHF,(cHF(i),i=1,6)
            !         write(6,105)una2p,asmall
            !         write(6,106)(y(1,1,j),j=1,4)
            !         write(6,107)(y(2,1,j),j=1,4)
            !         write(6,108)(y(3,1,j),j=1,4)
            !         write(6,109)(y(4,1,j),j=1,4)
            !         write(6,110)(y(1,2,j),j=1,4)
            !         write(6,111)(y(2,2,j),j=1,4)
            !         write(6,112)(y(3,2,j),j=1,4)
            !         write(6,113)(y(4,2,j),j=1,4)
            !         write(6,114)(y(1,3,j),j=1,4)
            !         write(6,115)(y(2,3,j),j=1,4)
            !         write(6,116)(y(3,3,j),j=1,4)
            !         write(6,117)(y(4,3,j),j=1,4)
            !         write(6,118)(c2a(j),j=1,4)
            !         write(6,119)(c2b(j),j=1,4)
            !         write(6,120)(c2c(j),j=1,4)
            !         write(6,121)(r20i(j),j=1,4)
            !         write(6,122)r20a
            !         write(6,123)(al20i(j),j=1,4)
            !         write(6,124)(al20i1(j),j=1,4)
            !         write(6,125)al20a
            !         write(6,126)(al2a(j),j=1,4)
            !         write(6,127)(rr0a(j),j=1,4)
            !         write(6,128)psi
            !         write(6,129)rsi
            !         write(6,130)(dip2(j),j=1,4)
            !         write(6,131)(adip2(j),j=1,4)
            !         write(6,132)(rdip2(j),j=1,4)
            !         write(6,133)(acut2(j),j=1,4)
            !         write(6,134)(rcut2(j),j=1,4)
            !         write(6,135)(ry(j),j=1,4)
            !         write(6,136)(aly(j),j=1,4)
            !         write(6,137)(dy(j),j=1,4)
            !         write(6,138)(dip2n(j),j=1,4)
            !         write(6,139)(adip2n(j),j=1,4)
            !         write(6,140)(rdip2n(j),j=1,4)
            !         write(6,141)(acut2n(j),j=1,4)
            !         write(6,142)(rcut2n(j),j=1,4)
        endif
        do i=1, nt
            r1s      = r(i,1)*r(i,1)
            r2s      = r(i,2)*r(i,2)
            r3s      = r(i,3)*r(i,3)
            rnag2    = af1*r1s+af2*r3s-af1*af2*r2s
            drnag2d1 = 2.e0_wp*af1*r(i,1)
            drnag2d2 = -2.e0_wp*af1*af2*r(i,2)
            drnag2d3 = 2.e0_wp*af2*r(i,3)
            !     Cosine of Na-F-H bond angle:
            hlp    = 1.e0_wp/(r(i,2)*r(i,3))
            csb    = 0.5e0_wp*(r2s+r3s-r1s)*hlp
            dcsbd1 = -2.e0_wp*r(i,1)*hlp
            dcsbd2 = 2.e0_wp*(1.e0_wp/r(i,3)-csb/r(i,2))
            dcsbd3 = 2.e0_wp*(1.e0_wp/r(i,2)-csb/r(i,3))
            csb    = 2.e0_wp*(1.e0_wp+csb)
            if( csb .ge. 4.e0_wp ) csb = 3.99999e0_wp
            l1 = idint(csb)+1
            l2 = l1+1
            if( l1 .eq. 4 ) l2=l1
            if( csb .le. bc(l1) ) csb = bc(l1)+1.e-13_wp
            if( csb .ge. bc(l1+1) ) csb = bc(l1+1)-1.e-13_wp
            s1       = 0.e0_wp
            ds1d1    = 0.e0_wp
            y2       = r(i,2)-rehf
            hlp1     = chf(1)*exp(-chf(2)*y2)
            hlp2     = exp(-chf(6)*y2)
            hlp22    = (-1.e0_wp-chf(1)+chf(3)*y2+chf(4)*y2**2+chf(5)*y2**3)*hlp2
            s2a      = dehf*(hlp1+hlp22) + una2p
            ds2ad2   = dehf*(-chf(2)*hlp1-chf(6)*hlp22 + (chf(3)+2.e0_wp*chf(4)*y2+3.e0_wp*chf(5)*y2**2)*hlp2)
            y3       = r(i,3)-renf
            hlp1     = cnf(1)*exp(-cnf(2)*y3)
            hlp2     = exp(-cnf(6)*y3)
            hlp22    = (-1.e0_wp-cnf(1)+cnf(3)*y3+cnf(4)*y3*y3+cnf(5)*y3**3)*hlp2
            s3       = denf*(hlp1+hlp22)
            ds3d3    = denf*(-cnf(2)*hlp1-cnf(6)*hlp22 + (cnf(3)+2.e0_wp*cnf(4)*y3+3.e0_wp*cnf(5)*y3**2)*hlp2)
            e(i)     = 0.e0_wp
            de(1,i)  = 0.e0_wp
            de(2,i)  = 0.e0_wp
            de(3,i)  = 0.e0_wp
            anorm    = 0.e0_wp
            danormd1 = 0.e0_wp
            danormd2 = 0.e0_wp
            danormd3 = 0.e0_wp
            do l = l1, l2
                hlp1  = y(1, 1, l)*exp(-y(2, 1, l)*r(i,1))
                hlp2  = y(3, 1, l)*exp(-y(4, 1, l)*r(i,1))
                t1    = hlp1+hlp2+s1
                dt1d1 = -y(2, 1, l)*hlp1-y(4, 1, l)*hlp2+ds1d1
                hlp   = al2a(l)*(r(i,3)+r(i,1)-rr0a(l))
                hlp1  = al2a(l)/(cosh(hlp)**2)
                hlp   = 0.5e0_wp*(1.e0_wp+tanh(hlp))
                hlp2  = r(i,3)+r(i,1)-rsi
                if( hlp2 .gt. 0.e0_wp )then
                    !     write(6,*)'R1+R3=',(r(i,1)+r(i,3))
                    ctanh  = exp(-psi/hlp2)
                    dctanh = psi*ctanh/(r(i,3)+r(i,1)-rsi)**2
                else
                    ctanh  = 0.e0_wp
                    dctanh = 0.e0_wp
                endif
                r20    = r20i(l)+(r20a-r20i(l))*hlp
                dr20   = 0.5e0_wp*(r20a-r20i(l))*hlp1
                al20   = al20i(l)+(al20i1(l)-al20i(l))*hlp + (al20a-al20i1(l))*ctanh
                dal20  = 0.5e0_wp*(al20i1(l)-al20i(l))*hlp1 + (al20a-al20i1(l))*dctanh
                hlp1   = y(1, 2, l)*exp(-y(2, 2, l)*r(i,2))
                hlp2   = exp(-(y(4, 2, l)+dy2(l)*ctanh)*r(i,2))
                hlp3   = (y(3, 2, l)+dy1(l)*ctanh)*hlp2
                t2     = hlp1+hlp3
                dt2d1  = dy1(l)*dctanh*hlp2 - dy2(l)*dctanh*r(i,2)*hlp3
                dt2d2  = -y(2, 2, l)*hlp1-(y(4, 2, l)+dy2(l)*ctanh)*hlp3
                dt2d3  = dt2d1
                hlp    = 0.5e0_wp*(1.e0_wp+tanh(al20*(r(i,2)-r20)))
                s2     = s2a+(t2-s2a)*hlp
                hlp1   = 0.5e0_wp/(cosh(al20*(r(i,2)-r20))**2)
                ds2d1  = (t2-s2a)*hlp1*(dal20*(r(i,2)-r20)-al20*dr20)
                ds2d2  = ds2ad2+(dt2d2-ds2ad2)*hlp + al20*hlp1*(t2-s2a)
                ds2d3  = ds2d1
                hlp    = aly(l)*(r(i,3)-ry(l))
                cy     = y(2, 3, l)+0.5e0_wp*dy(l)*(1.e0_wp+tanh(hlp))
                dcyd3  = 0.5e0_wp*dy(l)*aly(l)/(cosh(hlp)**2)
                hlp1   = y(1, 3, l)*exp(-cy*r(i,3))
                hlp2   = y(3, 3, l)*exp(-y(4, 3, l)*r(i,3))
                t3     = hlp1+hlp2
                dt3d3  = -(cy+dcyd3*r(i,3))*hlp1-y(4, 3, l)*hlp2
                coul1  = 0.5e0_wp*(s1+t1)
                dcl1d1 = 0.5e0_wp*(ds1d1+dt1d1)
                coul2  = 0.5e0_wp*(s2+t2)
                dcl2d1 = 0.5e0_wp*(ds2d1+dt2d1)
                dcl2d2 = 0.5e0_wp*(ds2d2+dt2d2)
                dcl2d3 = 0.5e0_wp*(ds2d3+dt2d3)
                coul3  = 0.5e0_wp*(s3+t3)
                dcl3d3 = 0.5e0_wp*(ds3d3+dt3d3)
                exch1  = 0.5e0_wp*(s1-t1)
                dex1d1 = 0.5e0_wp*(ds1d1-dt1d1)
                exch2  = 0.5e0_wp*(s2-t2)
                dex2d1 = 0.5e0_wp*(ds2d1-dt2d1)
                dex2d2 = 0.5e0_wp*(ds2d2-dt2d2)
                dex2d3 = 0.5e0_wp*(ds2d3-dt2d3)
                exch3  = 0.5e0_wp*(s3-t3)
                dex3d3 = 0.5e0_wp*(ds3d3-dt3d3)
                w      = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
                dwd1     = 2.e0_wp*(exch1-exch2)*(dex1d1-dex2d1) + 2.e0_wp*(exch2-exch3)*dex2d1 - 2.e0_wp*(exch3- &
                    & exch1)* dex1d1
                dwd2 = - 2.e0_wp*(exch1-exch2)*dex2d2 + 2.e0_wp*(exch2-exch3)*dex2d2
                dwd3     = - 2.e0_wp*(exch1-exch2)*dex2d3 + 2.e0_wp*(exch2-exch3)*(dex2d3-dex3d3) + 2.e0_wp*(exch3- &
                    & exch1)* dex3d3
                cplg2    = c2a(l)*exp(-c2b(l)*w-c2c(l)*(r(i,1)+r(i,2)+r(i,3)))
                dcd1     = (-c2b(l)*dwd1 - c2c(l))*cplg2
                dcd2     = (-c2b(l)*dwd2 - c2c(l))*cplg2
                dcd3     = (-c2b(l)*dwd3 - c2c(l))*cplg2
                rootw    = rac2*sqrt(w+cplg2*cplg2)
                surf     = (coul1+coul2+coul3+dehf-rootw)
                dsurfd1  = (dcl1d1+dcl2d1-0.25e0_wp*(dwd1+2.e0_wp*cplg2*dcd1)/rootw)
                dsurfd2  = (dcl2d2-0.25e0_wp*(dwd2+2.e0_wp*cplg2*dcd2)/rootw)
                dsurfd3  = (dcl2d3+dcl3d3-0.25e0_wp*(dwd3+2.e0_wp*cplg2*dcd3)/rootw)
                hlp1     = acut2(l)*(r(i,2)-rcut2(l))
                hlp2     = adip2(l)*(rnag2-rdip2(l)**2)
                hlp11    = tanh(hlp1)
                hlp22    = tanh(hlp2)
                hlp1     = 1.e0_wp/cosh(hlp1)
                hlp2     = 1.e0_wp/cosh(hlp2)
                surf     = surf - dip2(l)*hlp1*hlp2
                dsurfd1  = dsurfd1 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*drnag2d1
                dsurfd2  = dsurfd2 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*drnag2d2 + dip2(l)*hlp2*hlp1*hlp11*acut2(l)
                dsurfd3  = dsurfd3 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*drnag2d3
                hlp1     = acut2n(l)*(r(i,2)-rcut2n(l))
                hlp2     = adip2n(l)*(rnag2-rdip2n(l)**2)
                hlp11    = tanh(hlp1)
                hlp22    = tanh(hlp2)
                hlp1     = 1.e0_wp/cosh(hlp1)
                hlp2     = 1.e0_wp/cosh(hlp2)
                surf     = surf - dip2n(l)*hlp1*hlp2
                dsurfd1  = dsurfd1 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*drnag2d1
                dsurfd2  = dsurfd2 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*drnag2d2 + dip2n(l)*hlp2*hlp1*hlp11* acut2n(l)
                dsurfd3  = dsurfd3 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*drnag2d3
                hlp      = 1.e0_wp/( (1.e0_wp-(csb-bc(l))/(bc(l-1)-bc(l))) * (1.e0_wp-(csb-bc(l))/(bc(l+1)-bc(l))) )
                hlp1     = (csb-bc(l))**2*hlp
                cog      = exp(-ac(l)*hlp1)
                dcog     = -cog*ac(l)*(2.e0_wp*(csb-bc(l))*hlp +hlp1*(1.e0_wp/(bc(l-1)-csb)+1.e0_wp/(bc(l+1)-csb)))
                dcogd1   = dcog*dcsbd1
                dcogd2   = dcog*dcsbd2
                dcogd3   = dcog*dcsbd3
                e(i)     = e(i) + surf*cog
                de(1,i)  = de(1,i) + dsurfd1*cog+surf*dcogd1
                de(2,i)  = de(2,i) + dsurfd2*cog+surf*dcogd2
                de(3,i)  = de(3,i) + dsurfd3*cog+surf*dcogd3
                anorm    = anorm + cog
                danormd1 = danormd1 + dcogd1
                danormd2 = danormd2 + dcogd2
                danormd3 = danormd3 + dcogd3
            enddo
            hlp     = 1.e0_wp/anorm
            e(i)    = e(i)*hlp
            de(1,i) = (de(1,i) - e(i)*danormd1)*hlp
            de(2,i) = (de(2,i) - e(i)*danormd2)*hlp
            de(3,i) = (de(3,i) - e(i)*danormd3)*hlp
            do j=1, 3
                rdum(j) = r(i,j)
            enddo
            call lepscfu22(3, rdum, g_r, 3, cfu22, dcfu22, 3)
            e(i)    = e(i) + cfu22
            de(1,i) = de(1,i) + dcfu22(1)
            de(2,i) = de(2,i) + dcfu22(2)
            de(3,i) = de(3,i) + dcfu22(3)
            e(i)    = e(i)*cevau
            de(1,i) = de(1,i)*cevau
            de(2,i) = de(2,i)*cevau
            de(3,i) = de(3,i)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_nafh_leps_jhctp_2001_pot22_legacy_state_3ad9bc64

!> @brief Compatibility interface for pot22.
subroutine pot22( r, e, de, nt )
    use chempotpy_nafh_leps_jhctp_2001_pot22_legacy_state_3ad9bc64, only: dispatch
    use chempotpy_kinds,                                            only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    real(wp) :: de(3,nt)
    integer  :: nt
    call dispatch(r=r, e=e, de=de, nt=nt)
end subroutine pot22

!> @brief Evaluate lepscfu11 for the NaFH_LEPS_JHCTP_2001 potential.
subroutine lepscfu11( g_p_, r, g_r, ldg_r, cfu11, g_cfu11, ldg_cfu11 )
    !     correction function is the difference of two LEPS surfaces
    !     each triplet is weighted by a-f, which was optimezed with a GA
    use chempotpy_kinds, only: wp
    implicit none
    real(wp) :: a, b, c, d, e, f
    real(wp) :: dehf, j1_ap, rehf, shf_ac, yhf, j1_ac, snah, chi, tnah, q1_ap, q1_ac, j2_ap, j2_ac, w, leps_ac, &
        & leps_ap, q2_ap, c7, c6, shf_ap, q2_ac, thf, r2, r1, r3, snaf, ronaf, j3_ac, q3_ap, j3_ap, denah, renah, q3_ac, t2, &
        & t1, t3, tnaf, t4, y, thf_ap, thf_ac, tnaf_ap, tnaf_ac, tnah_ap, tnah_ac, cfu11, de, re, c5, beta, c1, c2, c3, c4, &
        & r(1,3)
    integer            :: ia, ib, ic
    integer, parameter :: g_pmax_ = 3
    integer            :: g_i_, g_p_, ldg_r, ldg_cfu11
    real(wp)           :: g_r(ldg_r,1,3)
    real(wp)           :: d4_p, d3_p, d17_b, d5_b, d6_b, d7_b, d8_b, d12_v, d2_v, d11_b
    real(wp)           :: d10_b, d1_w, d2_w, d3_v, d4_v, d5_v, d1_p, d7_v, d8_v, d2_p
    real(wp)           :: d10_v, d14_v, d4_b, g_r1(g_pmax_), g_r2(g_pmax_ &
        & ), g_r3(g_pmax_), g_y(g_pmax_), g_d1_w(g_pmax_), g_d2_w(g_pmax_), &
        & g_snaf(g_pmax_)
    real(wp) :: g_tnaf_ac(g_pmax_), g_tnaf_ap(g_pmax_), g_q3_ac &
        & (g_pmax_), g_q3_ap(g_pmax_), g_j3_ac(g_pmax_), g_j3_ap(g_pmax_), g_beta(g_pmax_), g_chi(g_pmax_), g_snah(g_pmax_), g_tnah_ac(g_pmax_ &
        & )
    real(wp) :: g_tnah_ap(g_pmax_), g_q1_ac(g_pmax_), g_q1_ap(g_pmax_), g_j1_ac(g_pmax_), g_j1_ap(g_pmax_), g_yhf(g_pmax_), g_shf_ac(g_pmax_), g_shf_ap(g_pmax_), g_thf_ac(g_pmax_), g_thf_ap(g_pmax_)
    real(wp) :: g_q2_ac(g_pmax_), g_q2_ap(g_pmax_), g_j2_ac(g_pmax_), g_j2_ap(g_pmax_), g_w(g_pmax_), g_leps_ac(g_pmax_), g_leps_ap(g_pmax_), g_cfu11(ldg_cfu11)
    integer  :: g_ehfid
    do ia = 1, 3
        do ib = 1, 1
            do ic = 1, 3
                if( ia .eq. ic ) g_r(ia,ib,ic) = 1.e0_wp
                if( ia .ne. ic ) g_r(ia,ib,ic) = 0.e0_wp
            enddo
        enddo
    enddo
    do g_i_ = 1, g_p_
        g_r1(g_i_) = real(g_r(g_i_,1,1), wp)
    enddo
    r1 = real(r(1,1), wp)
    do g_i_ = 1, g_p_
        g_r2(g_i_) = real(g_r(g_i_,1,2), wp)
    enddo
    r2 = real(r(1,2), wp)
    do g_i_ = 1, g_p_
        g_r3(g_i_) = real(g_r(g_i_,1,3), wp)
    enddo
    r3 = real(r(1,3), wp)
    a  = 3.58749e0_wp
    b  = 2.91300e0_wp
    c  = 2.06256e0_wp
    d  = 3.58749e0_wp
    e  = 3.11828e0_wp
    f  = 2.05279e0_wp
    !     NaF stuff
    !     singlets:  approx and acc the same (from paper)
    c1    = 0.620115e0_wp
    c2    = 0.694459e-1_wp
    c3    = 0.757570e0_wp
    c4    = 0.283055e0_wp
    c5    = 2.69998e0_wp
    ronaf = 3.6395e0_wp
    do g_i_ = 1, g_p_
        g_y(g_i_) = g_r3(g_i_)
    enddo
    y = r3 - ronaf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-c3) * g_y(g_i_)
    enddo
    d1_w = (-c3) * y
    do g_i_ = 1, g_p_
        g_d2_w(g_i_) = (-c5) * g_y(g_i_)
    enddo
    d2_w  = (-c5) * y
    d3_v  = c1 + c2 * y
    d5_v  = exp(d1_w)
    d2_p  = d5_v
    d8_v  = exp(d2_w)
    d1_p  = d8_v
    d5_b  = c4 * d1_p
    d8_b  = d3_v * d2_p
    d10_b = d5_v * c2
    do g_i_ = 1, g_p_
        g_snaf(g_i_) = d5_b * g_d2_w(g_i_) + d8_b * g_d1_w(g_i_) + d10_b * g_y(g_i_)
    enddo
    snaf = d3_v * d5_v + c4 * d8_v
    !     triplets:  approx and acc the same
    !     used 180 degree triplet from the paper
    t1 = 7.3251e0_wp
    t2 = 2.8774e0_wp
    t3 = 1930.0e0_wp
    t4 = 2.7091e0_wp
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-t2) * g_r3(g_i_)
    enddo
    d1_w = (-t2) * r3
    do g_i_ = 1, g_p_
        g_d2_w(g_i_) = (-t4) * g_r3(g_i_)
    enddo
    d2_w  = (-t4) * r3
    d3_v  = exp(d1_w)
    d2_p  = d3_v
    d7_v  = exp(d2_w)
    d1_p  = d7_v
    d6_b  = c * t3 * d1_p
    d10_b = c * t1 * d2_p
    do g_i_ = 1, g_p_
        g_tnaf_ac(g_i_) = d6_b * g_d2_w(g_i_) + d10_b * g_d1_w(g_i_) + c * g_snaf(g_i_)
    enddo
    tnaf_ac = (snaf + t1 * d3_v + t3 * d7_v) * c
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-t2) * g_r3(g_i_)
    enddo
    d1_w = (-t2) * r3
    do g_i_ = 1, g_p_
        g_d2_w(g_i_) = (-t4) * g_r3(g_i_)
    enddo
    d2_w  = (-t4) * r3
    d3_v  = exp(d1_w)
    d2_p  = d3_v
    d7_v  = exp(d2_w)
    d1_p  = d7_v
    d6_b  = f * t3 * d1_p
    d10_b = f * t1 * d2_p
    do g_i_ = 1, g_p_
        g_tnaf_ap(g_i_) = d6_b * g_d2_w(g_i_) + d10_b * g_d1_w(g_i_) + f * g_snaf(g_i_)
    enddo
    tnaf_ap = (snaf + t1 * d3_v + t3 * d7_v) * f
    do g_i_ = 1, g_p_
        g_q3_ac(g_i_) = 0.5e0_wp * g_tnaf_ac(g_i_) + 0.5e0_wp * g_snaf(g_i_)
    enddo
    q3_ac = 0.5e0_wp * (snaf + tnaf_ac)
    do g_i_ = 1, g_p_
        g_q3_ap(g_i_) = 0.5e0_wp * g_tnaf_ap(g_i_) + 0.5e0_wp * g_snaf(g_i_)
    enddo
    q3_ap = 0.5e0_wp * (snaf + tnaf_ap)
    do g_i_ = 1, g_p_
        g_j3_ac(g_i_) = (-0.5e0_wp) * g_tnaf_ac(g_i_) + 0.5e0_wp * g_snaf(g_i_)
    enddo
    j3_ac = 0.5e0_wp * (snaf - tnaf_ac)
    do g_i_ = 1, g_p_
        g_j3_ap(g_i_) = (-0.5e0_wp) * g_tnaf_ap(g_i_) + 0.5e0_wp * g_snaf(g_i_)
    enddo
    j3_ap = 0.5e0_wp * (snaf - tnaf_ap)
    !     NaH asymptote
    !     singlets:  approx and acc the same (from paper)
    renah = 7.37707e0_wp
    denah = 0.00291087e0_wp
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0e0_wp
    enddo
    beta = 0.713223e0_wp
    d4_v = r1 - renah
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r1(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d1_p = d2_v
    do g_i_ = 1, g_p_
        g_chi(g_i_) = d1_p * g_d1_w(g_i_)
    enddo
    chi  = d2_v
    d2_v = denah * chi
    d3_v = chi - 2.e0_wp
    d4_b = d2_v + d3_v * denah
    do g_i_ = 1, g_p_
        g_snah(g_i_) = d4_b * g_chi(g_i_)
    enddo
    snah = d2_v * d3_v
    !     triplets:  same for approx and acc
    !     used the 180 degree triplet
    t1 = 6.6796e0_wp
    t2 = 2.7839e0_wp
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-t2) * g_r1(g_i_)
    enddo
    d1_w = (-t2) * r1
    d3_v = exp(d1_w)
    d1_p = d3_v
    d6_b = a * t1 * d1_p
    do g_i_ = 1, g_p_
        g_tnah_ac(g_i_) = d6_b * g_d1_w(g_i_) + a * g_snah(g_i_)
    enddo
    tnah_ac = (snah + t1 * d3_v) * a
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-t2) * g_r1(g_i_)
    enddo
    d1_w = (-t2) * r1
    d3_v = exp(d1_w)
    d1_p = d3_v
    d6_b = d * t1 * d1_p
    do g_i_ = 1, g_p_
        g_tnah_ap(g_i_) = d6_b * g_d1_w(g_i_) + d * g_snah(g_i_)
    enddo
    tnah_ap = (snah + t1 * d3_v) * d
    do g_i_ = 1, g_p_
        g_q1_ac(g_i_) = 0.5e0_wp * g_tnah_ac(g_i_) + 0.5e0_wp * g_snah(g_i_)
    enddo
    q1_ac = 0.5e0_wp * (snah + tnah_ac)
    do g_i_ = 1, g_p_
        g_q1_ap(g_i_) = 0.5e0_wp * g_tnah_ap(g_i_) + 0.5e0_wp * g_snah(g_i_)
    enddo
    q1_ap = 0.5e0_wp * (snah + tnah_ap)
    do g_i_ = 1, g_p_
        g_j1_ac(g_i_) = (-0.5e0_wp) * g_tnah_ac(g_i_) + 0.5e0_wp * g_snah(g_i_)
    enddo
    j1_ac = 0.5e0_wp * (snah - tnah_ac)
    do g_i_ = 1, g_p_
        g_j1_ap(g_i_) = (-0.5e0_wp) * g_tnah_ap(g_i_) + 0.5e0_wp * g_snah(g_i_)
    enddo
    j1_ap = 0.5e0_wp * (snah - tnah_ap)
    !     HF asymptote
    !     accurate singlet
    dehf = 6.122e0_wp
    rehf = 1.733e0_wp
    c1   = 1.1622e0_wp
    c2   = -0.025647e0_wp
    c3   = 0.059062e0_wp
    c4   = 2.1042e0_wp
    do g_i_ = 1, g_p_
        g_yhf(g_i_) = g_r2(g_i_)
    enddo
    yhf  = r2 - c4
    d4_v = yhf * yhf
    d1_p = 2.0e0_wp * yhf
    d5_b = c3 * d1_p + c2
    do g_i_ = 1, g_p_
        g_beta(g_i_) = d5_b * g_yhf(g_i_)
    enddo
    beta = c1 + c2 * yhf + c3 * d4_v
    d4_v = r2 - rehf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d2_p = d2_v
    d4_v = (1.e0_wp - d2_v) * (1.e0_wp - d2_v)
    d1_p = 2.0e0_wp * (1.e0_wp - d2_v)
    d6_b = (-(dehf * d1_p)) * d2_p
    do g_i_ = 1, g_p_
        g_shf_ac(g_i_) = d6_b * g_d1_w(g_i_)
    enddo
    shf_ac = dehf * d4_v - dehf
    !     approx singlet (from paper)
    dehf = 5.77096e0_wp
    c1   = 0.441258e0_wp
    c2   = 3.88056e0_wp
    c3   = -1.441258e0_wp
    c4   = -1.2711e0_wp
    c5   = -1.37766e0_wp
    c6   = 0.168186e0_wp
    c7   = 2.07230e0_wp
    rehf = 1.7328e0_wp
    do g_i_ = 1, g_p_
        g_yhf(g_i_) = g_r2(g_i_)
    enddo
    yhf = r2 - rehf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-c2) * g_yhf(g_i_)
    enddo
    d1_w = (-c2) * yhf
    do g_i_ = 1, g_p_
        g_d2_w(g_i_) = (-c7) * g_yhf(g_i_)
    enddo
    d2_w  = (-c7) * yhf
    d2_v  = exp(d1_w)
    d4_p  = d2_v
    d7_v  = yhf * yhf
    d3_p  = 2.0e0_wp * yhf
    d10_v = yhf ** ( 3 - 2)
    d10_v = d10_v * yhf
    d2_p  = 3 *  d10_v
    d10_v = d10_v * yhf
    d12_v = c3 + c4 * yhf + c5 * d7_v + c6 * d10_v
    d14_v = exp(d2_w)
    d1_p  = d14_v
    d5_b  = dehf * d14_v
    d7_b  = dehf * d12_v * d1_p
    d11_b = d5_b * c6 * d2_p + d5_b * c5 * d3_p + d5_b * c4
    d17_b = dehf * c1 * d4_p
    do g_i_ = 1, g_p_
        g_shf_ap(g_i_) = d7_b * g_d2_w(g_i_) + d11_b * g_yhf(g_i_) + d17_b * g_d1_w(g_i_)
    enddo
    shf_ap = dehf * (c1 * d2_v + d12_v * d14_v)
    !     triplets:  same for acc and approx, taken from constant beta fit t
    !data
    dehf = 5.77096e0_wp
    rehf = 1.7328e0_wp
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0e0_wp
    enddo
    beta = 1.26686217e0_wp
    d4_v = r2 - rehf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d2_p = d2_v
    d4_v = (1.e0_wp + d2_v) * (1.e0_wp + d2_v)
    d1_p = 2.0e0_wp * (1.e0_wp + d2_v)
    d7_b = b * dehf * d1_p * d2_p
    do g_i_ = 1, g_p_
        g_thf_ac(g_i_) = d7_b * g_d1_w(g_i_)
    enddo
    thf_ac = (dehf * d4_v - dehf) * b
    d4_v   = r2 - rehf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d2_p = d2_v
    d4_v = (1.e0_wp + d2_v) * (1.e0_wp + d2_v)
    d1_p = 2.0e0_wp * (1.e0_wp + d2_v)
    d7_b = e * dehf * d1_p * d2_p
    do g_i_ = 1, g_p_
        g_thf_ap(g_i_) = d7_b * g_d1_w(g_i_)
    enddo
    thf_ap = (dehf * d4_v - dehf) * e
    do g_i_ = 1, g_p_
        g_q2_ac(g_i_) = 0.5e0_wp * g_thf_ac(g_i_) + 0.5e0_wp * g_shf_ac(g_i_ )
    enddo
    q2_ac = 0.5e0_wp * (shf_ac + thf_ac)
    do g_i_ = 1, g_p_
        g_q2_ap(g_i_) = 0.5e0_wp * g_thf_ap(g_i_) + 0.5e0_wp * g_shf_ap(g_i_ )
    enddo
    q2_ap = 0.5e0_wp * (shf_ap + thf_ap)
    do g_i_ = 1, g_p_
        g_j2_ac(g_i_) = (-0.5e0_wp) * g_thf_ac(g_i_) + 0.5e0_wp * g_shf_ac(g_i_)
    enddo
    j2_ac = 0.5e0_wp * (shf_ac - thf_ac)
    do g_i_ = 1, g_p_
        g_j2_ap(g_i_) = (-0.5e0_wp) * g_thf_ap(g_i_) + 0.5e0_wp * g_shf_ap(g_i_)
    enddo
    j2_ap = 0.5e0_wp * (shf_ap - thf_ap)
    !     compute LEPS functions
    do g_i_ = 1, g_p_
        g_j1_ac(g_i_) = g_j1_ac(g_i_)
    enddo
    j1_ac = j1_ac
    do g_i_ = 1, g_p_
        g_j2_ac(g_i_) = g_j2_ac(g_i_)
    enddo
    j2_ac = j2_ac
    do g_i_ = 1, g_p_
        g_j3_ac(g_i_) = g_j3_ac(g_i_)
    enddo
    j3_ac = j3_ac
    d4_b  = -j1_ac + (-j2_ac) + j3_ac + j3_ac
    d8_b  = -j3_ac + (-j1_ac) + j2_ac + j2_ac
    d5_b  = -j3_ac + (-j2_ac) + j1_ac + j1_ac
    do g_i_ = 1, g_p_
        g_w(g_i_) = d4_b * g_j3_ac(g_i_) + d8_b * g_j2_ac(g_i_) + d5_b * g_j1_ac(g_i_)
    enddo
    w    = j1_ac * j1_ac + j2_ac * j2_ac + j3_ac * j3_ac - j1_ac * j2_ac - j2_ac * j3_ac - j3_ac * j1_ac
    d7_v = sqrt(w)
    !        if ( w .gt. 0.0d0 ) then
    d1_p = 1.0e0_wp / (2.0e0_wp *  d7_v)
    !        else
    !           call ehufDO (9,w, d7_v, d1_p,
    !     +g_ehfid,
    !     +482)
    !        endif
    do g_i_ = 1, g_p_
        g_leps_ac(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ac(g_i_) + g_q2_ac(g_i_) + g_q1_ac(g_i_)
    enddo
    leps_ac = q1_ac + q2_ac + q3_ac - d7_v
    do g_i_ = 1, g_p_
        g_j1_ap(g_i_) = g_j1_ap(g_i_)
    enddo
    j1_ap = j1_ap
    do g_i_ = 1, g_p_
        g_j2_ap(g_i_) = g_j2_ap(g_i_)
    enddo
    j2_ap = j2_ap
    do g_i_ = 1, g_p_
        g_j3_ap(g_i_) = g_j3_ap(g_i_)
    enddo
    j3_ap = j3_ap
    d4_b  = -j1_ap + (-j2_ap) + j3_ap + j3_ap
    d8_b  = -j3_ap + (-j1_ap) + j2_ap + j2_ap
    d5_b  = -j3_ap + (-j2_ap) + j1_ap + j1_ap
    do g_i_ = 1, g_p_
        g_w(g_i_) = d4_b * g_j3_ap(g_i_) + d8_b * g_j2_ap(g_i_) + d5_b * g_j1_ap(g_i_)
    enddo
    w    = j1_ap * j1_ap + j2_ap * j2_ap + j3_ap * j3_ap - j1_ap * j2_ap - j2_ap * j3_ap - j3_ap * j1_ap
    d7_v = sqrt(w)
    !        if ( w .gt. 0.0d0 ) then
    d1_p = 1.0e0_wp / (2.0e0_wp *  d7_v)
    !        else
    !           call ehufDO (9,w, d7_v, d1_p,
    !     +g_ehfid,
    !     +524)
    !        endif
    do g_i_ = 1, g_p_
        g_leps_ap(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ap(g_i_) + g_q2_ap(g_i_) + g_q1_ap(g_i_)
    enddo
    leps_ap = q1_ap + q2_ap + q3_ap - d7_v
    do g_i_ = 1, g_p_
        g_cfu11(g_i_) = -g_leps_ap(g_i_) + g_leps_ac(g_i_)
    enddo
    cfu11 = leps_ac - leps_ap + 6.122e0_wp - 5.77096e0_wp
end subroutine lepscfu11

!> @brief Evaluate lepscfu22 for the NaFH_LEPS_JHCTP_2001 potential.
subroutine lepscfu22( g_p_, r, g_r, ldg_r, cfu22, g_cfu22, ldg_cfu22 )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp) :: dehf, j1_ap, rehf, shf_ac, yhf, j1_ac, snah, chi, tnah, q1_ap, q1_ac, j2_ap, j2_ac, w, leps_ac, &
        & leps_ap, q2_ap, c7, c6, shf_ap, q2_ac, thf, r2, r1, r3, snaf, ronaf, a, b, c, d, e, f, j3_ac, q3_ap, j3_ap, denah, &
        & renah, q3_ac, t2, t1, t3, tnaf, t4, y, thf_ap, thf_ac, tnaf_ap, tnaf_ac, tnah_ap, tnah_ac, r(1,3)
    real(wp) :: shf_u11, nap, t3nah, t4nah, r0, swit, shfa, alph0, ynaf, snaf_ap, denaf, renaf, t1nah, t2nah, c2naf, snaf_ac
    real(wp) :: c1nf, c2nf, c3nf, c4nf, c5nf, c6nf, c7nf, hlp1, &
        & hlp2, hlp22, thf_ap_as, t1hf, dt1hf, t2hf, dt2hf, t, b1, b2, e1, &
        & e2, hlp, a3, b3, cc3, d3, e3, f3, cfu22, c5, c4, c3, c2, c1, cc1, &
        & cc2, f1, f2, beta
    integer            :: ia, ib, ic
    integer, parameter :: g_pmax_ = 3
    integer            :: g_i_, g_p_, ldg_r, ldg_cfu22
    real(wp)           :: d4_p, d3_p, d17_b, d8_b, d12_v, d2_p, d2_w, d10_v, d2_v, d11_b
    real(wp)           :: d2_b, d1_w, d3_v, d3_b, d1_p, d4_v, d5_v, d6_v, &
        & d7_v, d8_v
    real(wp) :: d14_v, d4_b, d5_b, d6_b, d7_b, g_r1(g_pmax_), g_r(ldg_r,1,3), g_r2(g_pmax_), g_r3(g_pmax_), g_ynaf(g_pmax_)
    real(wp) :: g_d1_w(g_pmax_), g_hlp1(g_pmax_), g_hlp2(g_pmax_), g_hlp22(g_pmax_), g_snaf_ap(g_pmax_), g_beta(g_pmax_), g_snaf_ac(g_pmax_), g_hlp(g_pmax_), g_tnaf_ac(g_pmax_), g_tnaf_ap(g_pmax_ &
        & )
    real(wp) :: g_q3_ac(g_pmax_), g_q3_ap(g_pmax_), g_j3_ac(g_pmax_), g_j3_ap(g_pmax_), g_d2_w(g_pmax_), g_tnah_ac(g_pmax_), g_tnah_ap(g_pmax_), g_q1_ac(g_pmax_), g_q1_ap(g_pmax_), g_j1_ac(g_pmax_)
    real(wp) :: g_j1_ap(g_pmax_), g_thf_ac(g_pmax_), g_thf_ap(g_pmax_), g_thf(g_pmax_), g_yhf(g_pmax_), g_shf_u11(g_pmax_), g_shfa(g_pmax_), g_swit(g_pmax_), g_shf_ac(g_pmax_), g_thf_ap_as(g_pmax_)
    real(wp) :: g_shf_ap(g_pmax_), g_q2_ac(g_pmax_), g_q2_ap(g_pmax_), g_j2_ac(g_pmax_), g_j2_ap(g_pmax_), g_w(g_pmax_), g_leps_ac(g_pmax_), g_leps_ap(g_pmax_), g_cfu22(ldg_cfu22)
    integer  :: g_ehfid
    do ia = 1, 3
        do ib = 1, 1
            do ic = 1, 3
                if( ia .eq. ic ) g_r(ia,ib,ic) = 1.e0_wp
                if( ia .ne. ic ) g_r(ia,ib,ic) = 0.e0_wp
            enddo
        enddo
    enddo
    do g_i_ = 1, g_p_
        g_r1(g_i_) = g_r(g_i_,1,1)
    enddo
    r1 = r(1,1)
    do g_i_ = 1, g_p_
        g_r2(g_i_) = g_r(g_i_,1,2)
    enddo
    r2 = r(1,2)
    do g_i_ = 1, g_p_
        g_r3(g_i_) = g_r(g_i_,1,3)
    enddo
    r3  = r(1,3)
    b1  = 5.125e0_wp
    b2  = 5.125e0_wp
    cc1 = 1.92e0_wp
    cc2 = 9.01e0_wp
    e1  = 5.54e0_wp
    e2  = 4.15e0_wp
    f1  = 1.62e0_wp
    f2  = 8.80e0_wp
    !     NaF stuff
    !     singlets:  approx (from paper)
    c1nf  = 0.168364e0_wp
    c2nf  = 2.29333e0_wp
    c3nf  = -1.168364e0_wp
    c4nf  = -0.084862e0_wp
    c5nf  = -0.0109707e0_wp
    c6nf  = 0.00246582e0_wp
    c7nf  = 0.377012e0_wp
    denaf = 4.3830525e0_wp
    renaf = 3.6359e0_wp
    do g_i_ = 1, g_p_
        g_ynaf(g_i_) = g_r3(g_i_)
    enddo
    ynaf = r3 - renaf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-c2nf) * g_ynaf(g_i_)
    enddo
    d1_w = (-c2nf) * ynaf
    d2_v = exp(d1_w)
    d1_p = d2_v
    d3_b = c1nf * d1_p
    do g_i_ = 1, g_p_
        g_hlp1(g_i_) = d3_b * g_d1_w(g_i_)
    enddo
    hlp1 = c1nf * d2_v
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-c7nf) * g_ynaf(g_i_)
    enddo
    d1_w = (-c7nf) * ynaf
    d2_v = exp(d1_w)
    d1_p = d2_v
    do g_i_ = 1, g_p_
        g_hlp2(g_i_) = d1_p * g_d1_w(g_i_)
    enddo
    hlp2 = d2_v
    d4_v = ynaf * ynaf
    d2_p = 2.0e0_wp * ynaf
    d7_v = ynaf ** ( 3 - 2)
    d7_v = d7_v * ynaf
    d1_p = 3 *  d7_v
    d7_v = d7_v * ynaf
    d5_b = c6nf * d1_p + c5nf * d2_p + c4nf
    do g_i_ = 1, g_p_
        g_hlp22(g_i_) = d5_b * g_ynaf(g_i_)
    enddo
    hlp22 = c3nf + c4nf * ynaf + c5nf * d4_v + c6nf * d7_v
    d5_b  = denaf * hlp2
    d6_b  = denaf * hlp22
    do g_i_ = 1, g_p_
        g_snaf_ap(g_i_) = d6_b * g_hlp2(g_i_) + d5_b * g_hlp22(g_i_) + denaf * g_hlp1(g_i_)
    enddo
    snaf_ap = denaf * (hlp1 + hlp22 * hlp2)
    !     singlet:  accurate from fit to experimental
    c1    = 0.32453e0_wp
    c2    = 1.5102e0_wp
    c3    = 3.0938e0_wp
    c4    = 1.7107e0_wp
    denaf = 4.94e0_wp
    renaf = 3.6395e0_wp
    if( (r3 / c3) .ne. 0.0e0_wp )then
        d3_v = (r3 / c3) ** ( c4 - 2.0e0_wp)
        d3_v = d3_v * (r3 / c3)
        d2_p = c4 *  d3_v
        d3_v = d3_v * (r3 / c3)
    else
        !          ((r3 / c3) = 0)
        d3_v = (r3 / c3) **  c4
        d2_p = 0.0e0_wp
    endif
    if( (r3 / c3) .ne. 0.0e0_wp )then
        d6_v = (r3 / c3) ** ( c4 - 2.0e0_wp)
        d6_v = d6_v * (r3 / c3)
        d1_p = c4 *  d6_v
        d6_v = d6_v * (r3 / c3)
    else
        !          ((r3 / c3) = 0)
        d6_v = (r3 / c3) **  c4
        d1_p = 0.0e0_wp
    endif
    d7_v = c1 + d6_v
    d8_v = (c2 + d3_v) / d7_v
    d7_b = c1 * ((-d8_v) / d7_v) * d1_p * (1.0e0_wp / c3) + c1 * (1.0e0_wp / d7_v) * d2_p * (1.0e0_wp / c3)
    do g_i_ = 1, g_p_
        g_beta(g_i_) = d7_b * g_r3(g_i_)
    enddo
    beta = c1 * d8_v
    d4_v = r3 - renaf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r3(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d2_p = d2_v
    d4_v = (1.e0_wp - d2_v) * (1.e0_wp - d2_v)
    d1_p = 2.0e0_wp * (1.e0_wp - d2_v)
    d6_b = (-(denaf * d1_p)) * d2_p
    do g_i_ = 1, g_p_
        g_snaf_ac(g_i_) = d6_b * g_d1_w(g_i_)
    enddo
    snaf_ac = denaf * d4_v - denaf
    !     triplets:  approx and acc the same
    !     anti-Morse from a constant beta fit to the approx curve
    denaf = 4.49e0_wp
    renaf = 3.6395e0_wp
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0e0_wp
    enddo
    beta = 0.69614e0_wp
    d4_v = r3 - renaf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r3(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d1_p = d2_v
    do g_i_ = 1, g_p_
        g_hlp(g_i_) = d1_p * g_d1_w(g_i_)
    enddo
    hlp  = d2_v
    d2_v = hlp * hlp
    d1_p = 2.0e0_wp * hlp
    d2_b = 0.5e0_wp * denaf
    d5_b = d2_b * (cc2 * 2.e0_wp) + d2_b * cc1 * d1_p
    do g_i_ = 1, g_p_
        g_tnaf_ac(g_i_) = d5_b * g_hlp(g_i_)
    enddo
    tnaf_ac = 0.5e0_wp * denaf * (cc1 * d2_v + cc2 * 2.e0_wp * hlp)
    d2_v    = hlp * hlp
    d1_p    = 2.0e0_wp * hlp
    d2_b    = 0.5e0_wp * denaf
    d5_b    = d2_b * (f2 * 2.e0_wp) + d2_b * f1 * d1_p
    do g_i_ = 1, g_p_
        g_tnaf_ap(g_i_) = d5_b * g_hlp(g_i_)
    enddo
    tnaf_ap = 0.5e0_wp * denaf * (f1 * d2_v + f2 * 2.e0_wp * hlp)
    do g_i_ = 1, g_p_
        g_q3_ac(g_i_) = 0.5e0_wp * g_tnaf_ac(g_i_) + 0.5e0_wp * g_snaf_ac(g_i_)
    enddo
    q3_ac = 0.5e0_wp * (snaf_ac + tnaf_ac)
    do g_i_ = 1, g_p_
        g_q3_ap(g_i_) = 0.5e0_wp * g_tnaf_ap(g_i_) + 0.5e0_wp * g_snaf_ap(g_i_)
    enddo
    q3_ap = 0.5e0_wp * (snaf_ap + tnaf_ap)
    do g_i_ = 1, g_p_
        g_j3_ac(g_i_) = (-0.5e0_wp) * g_tnaf_ac(g_i_) + 0.5e0_wp * g_snaf_ac (g_i_)
    enddo
    j3_ac = 0.5e0_wp * (snaf_ac - tnaf_ac)
    do g_i_ = 1, g_p_
        g_j3_ap(g_i_) = (-0.5e0_wp) * g_tnaf_ap(g_i_) + 0.5e0_wp * g_snaf_ap (g_i_)
    enddo
    j3_ap = 0.5e0_wp * (snaf_ap - tnaf_ap)
    !     NaH asymptote
    !     singlets:  approx and acc the same (from paper)
    snah = 0.e0_wp
    !     triplets:  same for approx and acc
    !     used the 180 degree triplet
    t1nah = 1040.7e0_wp
    t2nah = 4.5911e0_wp
    t3nah = 1.4941e0_wp
    t4nah = 4.5525e0_wp
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-t2nah) * g_r1(g_i_)
    enddo
    d1_w = (-t2nah) * r1
    do g_i_ = 1, g_p_
        g_d2_w(g_i_) = (-t4nah) * g_r1(g_i_)
    enddo
    d2_w = (-t4nah) * r1
    d2_v = exp(d1_w)
    d2_p = d2_v
    d5_v = exp(d2_w)
    d1_p = d5_v
    d5_b = t3nah * d1_p
    d7_b = t1nah * d2_p
    do g_i_ = 1, g_p_
        g_tnah_ac(g_i_) = d5_b * g_d2_w(g_i_) + d7_b * g_d1_w(g_i_)
    enddo
    tnah_ac = t1nah * d2_v + t3nah * d5_v
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-t2nah) * g_r1(g_i_)
    enddo
    d1_w = (-t2nah) * r1
    do g_i_ = 1, g_p_
        g_d2_w(g_i_) = (-t4nah) * g_r1(g_i_)
    enddo
    d2_w = (-t4nah) * r1
    d2_v = exp(d1_w)
    d2_p = d2_v
    d5_v = exp(d2_w)
    d1_p = d5_v
    d5_b = t3nah * d1_p
    d7_b = t1nah * d2_p
    do g_i_ = 1, g_p_
        g_tnah_ap(g_i_) = d5_b * g_d2_w(g_i_) + d7_b * g_d1_w(g_i_)
    enddo
    tnah_ap = t1nah * d2_v + t3nah * d5_v
    do g_i_ = 1, g_p_
        g_q1_ac(g_i_) = 0.5e0_wp * g_tnah_ac(g_i_)
    enddo
    q1_ac = 0.5e0_wp * (snah + tnah_ac)
    do g_i_ = 1, g_p_
        g_q1_ap(g_i_) = 0.5e0_wp * g_tnah_ap(g_i_)
    enddo
    q1_ap = 0.5e0_wp * (snah + tnah_ap)
    do g_i_ = 1, g_p_
        g_j1_ac(g_i_) = (-0.5e0_wp) * g_tnah_ac(g_i_)
    enddo
    j1_ac = 0.5e0_wp * (snah - tnah_ac)
    do g_i_ = 1, g_p_
        g_j1_ap(g_i_) = (-0.5e0_wp) * g_tnah_ap(g_i_)
    enddo
    j1_ap = 0.5e0_wp * (snah - tnah_ap)
    !     HF asymptote
    !     triplets:
    !     anti-morse from a
    !     constant beta fit to MRDCI data
    !     same for approx and accurate
    dehf = 5.77096e0_wp
    rehf = 1.7328e0_wp
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0e0_wp
    enddo
    beta = 1.2669e0_wp
    d4_v = r2 - rehf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d1_p = d2_v
    do g_i_ = 1, g_p_
        g_hlp(g_i_) = d1_p * g_d1_w(g_i_)
    enddo
    hlp  = d2_v
    d3_v = hlp * hlp
    d1_p = 2.0e0_wp * hlp
    d2_b = 0.5e0_wp * dehf
    d6_b = d2_b * e2 * d1_p + d2_b * (e1 * 2.e0_wp)
    do g_i_ = 1, g_p_
        g_thf_ac(g_i_) = d6_b * g_hlp(g_i_)
    enddo
    thf_ac = 0.5e0_wp * dehf * (e1 * 2.e0_wp * hlp + e2 * d3_v)
    dehf   = 5.77096e0_wp
    rehf   = 1.7328e0_wp
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0e0_wp
    enddo
    beta = 1.2669e0_wp
    d4_v = r2 - rehf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d1_p = d2_v
    do g_i_ = 1, g_p_
        g_hlp(g_i_) = d1_p * g_d1_w(g_i_)
    enddo
    hlp  = d2_v
    d3_v = hlp * hlp
    d1_p = 2.0e0_wp * hlp
    d2_b = 0.5e0_wp * dehf
    d6_b = d2_b * b2 * d1_p + d2_b * (b1 * 2.e0_wp)
    do g_i_ = 1, g_p_
        g_thf_ap(g_i_) = d6_b * g_hlp(g_i_)
    enddo
    thf_ap = 0.5e0_wp * dehf * (b1 * 2.e0_wp * hlp + b2 * d3_v)
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-15.e0_wp) * g_r2(g_i_)
    enddo
    d1_w = (-15.e0_wp) * (r2 - 2.e0_wp)
    d2_v = exp(d1_w)
    d1_p = d2_v
    do g_i_ = 1, g_p_
        g_thf(g_i_) = d1_p * g_d1_w(g_i_)
    enddo
    thf = d2_v
    !     accurate singlet
    !     uses u11 singlet which is
    dehf = 6.122e0_wp
    rehf = 1.733e0_wp
    c1   = 1.1622e0_wp
    c2   = -0.025647e0_wp
    c3   = 0.059062e0_wp
    c4   = 2.1042e0_wp
    do g_i_ = 1, g_p_
        g_yhf(g_i_) = g_r2(g_i_)
    enddo
    yhf  = r2 - c4
    d4_v = yhf * yhf
    d1_p = 2.0e0_wp * yhf
    d5_b = c3 * d1_p + c2
    do g_i_ = 1, g_p_
        g_beta(g_i_) = d5_b * g_yhf(g_i_)
    enddo
    beta = c1 + c2 * yhf + c3 * d4_v
    d4_v = r2 - rehf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d2_p = d2_v
    d4_v = (1.e0_wp - d2_v) * (1.e0_wp - d2_v)
    d1_p = 2.0e0_wp * (1.e0_wp - d2_v)
    d6_b = (-(dehf * d1_p)) * d2_p
    do g_i_ = 1, g_p_
        g_shf_u11(g_i_) = d6_b * g_d1_w(g_i_)
    enddo
    shf_u11 = dehf * d4_v - dehf
    nap     = 2.0973375e0_wp
    do g_i_ = 1, g_p_
        g_shfa(g_i_) = g_shf_u11(g_i_)
    enddo
    shfa  = shf_u11 + nap
    alph0 = 15.e0_wp
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-alph0) * g_thf(g_i_) + alph0 * g_shfa(g_i_)
    enddo
    d1_w = alph0 * (shfa - thf)
    d2_v = tanh (d1_w)
    d1_p = 1.0e0_wp - ( d2_v *  d2_v)
    do g_i_ = 1, g_p_
        g_swit(g_i_) = d1_p * g_d1_w(g_i_)
    enddo
    swit = d2_v
    d6_v = 0.5e0_wp * (shfa - thf)
    d6_b = (-swit) * 0.5e0_wp
    d7_b = d6_b + 0.5e0_wp
    d8_b = -d6_b + 0.5e0_wp
    do g_i_ = 1, g_p_
        g_shf_ac(g_i_) = (-d6_v) * g_swit(g_i_) + d8_b * g_thf(g_i_) + d7_b * g_shfa(g_i_)
    enddo
    shf_ac = 0.5e0_wp * (shfa + thf) - d6_v * swit
    !     approx singlet (from paper)
    !     uses u11 singlet, which is
    dehf = 5.77096e0_wp
    c1   = 0.441258e0_wp
    c2   = 3.88056e0_wp
    c3   = -1.441258e0_wp
    c4   = -1.2711e0_wp
    c5   = -1.37766e0_wp
    c6   = 0.168186e0_wp
    c7   = 2.07230e0_wp
    rehf = 1.7328e0_wp
    do g_i_ = 1, g_p_
        g_yhf(g_i_) = g_r2(g_i_)
    enddo
    yhf = r2 - rehf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-c2) * g_yhf(g_i_)
    enddo
    d1_w = (-c2) * yhf
    do g_i_ = 1, g_p_
        g_d2_w(g_i_) = (-c7) * g_yhf(g_i_)
    enddo
    d2_w  = (-c7) * yhf
    d2_v  = exp(d1_w)
    d4_p  = d2_v
    d7_v  = yhf * yhf
    d3_p  = 2.0e0_wp * yhf
    d10_v = yhf ** ( 3 - 2)
    d10_v = d10_v * yhf
    d2_p  = 3 *  d10_v
    d10_v = d10_v * yhf
    d12_v = c3 + c4 * yhf + c5 * d7_v + c6 * d10_v
    d14_v = exp(d2_w)
    d1_p  = d14_v
    d5_b  = dehf * d14_v
    d7_b  = dehf * d12_v * d1_p
    d11_b = d5_b * c6 * d2_p + d5_b * c5 * d3_p + d5_b * c4
    d17_b = dehf * c1 * d4_p
    do g_i_ = 1, g_p_
        g_shf_u11(g_i_) = d7_b * g_d2_w(g_i_) + d11_b * g_yhf(g_i_) + d17_b * g_d1_w(g_i_)
    enddo
    shf_u11 = dehf * (c1 * d2_v + d12_v * d14_v)
    nap     = 2.0973375e0_wp
    do g_i_ = 1, g_p_
        g_shfa(g_i_) = g_shf_u11(g_i_)
    enddo
    shfa  = shf_u11 + nap
    alph0 = 15.e0_wp
    r0    = 3.1305e0_wp
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = alph0 * g_r2(g_i_)
    enddo
    d1_w = alph0 * (r2 - r0)
    d2_v = tanh (d1_w)
    d1_p = 1.0e0_wp - ( d2_v *  d2_v)
    d4_b = 0.5e0_wp * d1_p
    do g_i_ = 1, g_p_
        g_swit(g_i_) = d4_b * g_d1_w(g_i_)
    enddo
    swit  = 0.5e0_wp * (1.e0_wp + d2_v)
    t1hf  = 1807.6e0_wp
    dt1hf = 0.e0_wp
    t2hf  = 2.8815e0_wp
    dt2hf = 0.e0_wp
    d2_b  = -(t2hf + dt2hf)
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = d2_b * g_r2(g_i_)
    enddo
    d1_w = (-(t2hf + dt2hf)) * r2
    d2_v = exp(d1_w)
    d1_p = d2_v
    d3_b = (t1hf + dt1hf) * d1_p
    do g_i_ = 1, g_p_
        g_thf_ap_as(g_i_) = d3_b * g_d1_w(g_i_)
    enddo
    thf_ap_as = (t1hf + dt1hf) * d2_v
    d3_v      = thf_ap_as - shfa
    d2_b      = 1.0e0_wp + (-swit)
    do g_i_ = 1, g_p_
        g_shf_ap(g_i_) = d3_v * g_swit(g_i_) + swit * g_thf_ap_as(g_i_ ) + d2_b * g_shfa(g_i_)
    enddo
    shf_ap = shfa + d3_v * swit
    !     compute coulomb and exchange integrals
    do g_i_ = 1, g_p_
        g_q2_ac(g_i_) = 0.5e0_wp * g_thf_ac(g_i_) + 0.5e0_wp * g_shf_ac(g_i_ )
    enddo
    q2_ac = 0.5e0_wp * (shf_ac + thf_ac)
    do g_i_ = 1, g_p_
        g_q2_ap(g_i_) = 0.5e0_wp * g_thf_ap(g_i_) + 0.5e0_wp * g_shf_ap(g_i_ )
    enddo
    q2_ap = 0.5e0_wp * (shf_ap + thf_ap)
    do g_i_ = 1, g_p_
        g_j2_ac(g_i_) = (-0.5e0_wp) * g_thf_ac(g_i_) + 0.5e0_wp * g_shf_ac(g_i_)
    enddo
    j2_ac = 0.5e0_wp * (shf_ac - thf_ac)
    do g_i_ = 1, g_p_
        g_j2_ap(g_i_) = (-0.5e0_wp) * g_thf_ap(g_i_) + 0.5e0_wp * g_shf_ap(g_i_)
    enddo
    j2_ap = 0.5e0_wp * (shf_ap - thf_ap)
    !     compute LEPS functions
    do g_i_ = 1, g_p_
        g_j1_ac(g_i_) = g_j1_ac(g_i_)
    enddo
    j1_ac = j1_ac
    do g_i_ = 1, g_p_
        g_j2_ac(g_i_) = g_j2_ac(g_i_)
    enddo
    j2_ac = j2_ac
    do g_i_ = 1, g_p_
        g_j3_ac(g_i_) = g_j3_ac(g_i_)
    enddo
    j3_ac = j3_ac
    d4_b  = -j1_ac + (-j2_ac) + j3_ac + j3_ac
    d8_b  = -j3_ac + (-j1_ac) + j2_ac + j2_ac
    d5_b  = -j3_ac + (-j2_ac) + j1_ac + j1_ac
    do g_i_ = 1, g_p_
        g_w(g_i_) = d4_b * g_j3_ac(g_i_) + d8_b * g_j2_ac(g_i_) + d5_b * g_j1_ac(g_i_)
    enddo
    w    = j1_ac * j1_ac + j2_ac * j2_ac + j3_ac * j3_ac - j1_ac * j2_ac - j2_ac * j3_ac - j3_ac * j1_ac
    d7_v = sqrt(w)
    !        if ( w .gt. 0.0d0 ) then
    d1_p = 1.0e0_wp / (2.0e0_wp *  d7_v)
    !        else
    !           call ehufDO (9,w, d7_v, d1_p,
    !     +g_ehfid,
    !     +704)
    !        endif
    do g_i_ = 1, g_p_
        g_leps_ac(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ac(g_i_) + g_q2_ac(g_i_) + g_q1_ac(g_i_)
    enddo
    leps_ac = q1_ac + q2_ac + q3_ac - d7_v + 6.122e0_wp
    do g_i_ = 1, g_p_
        g_j1_ap(g_i_) = g_j1_ap(g_i_)
    enddo
    j1_ap = j1_ap
    do g_i_ = 1, g_p_
        g_j2_ap(g_i_) = g_j2_ap(g_i_)
    enddo
    j2_ap = j2_ap
    do g_i_ = 1, g_p_
        g_j3_ap(g_i_) = g_j3_ap(g_i_)
    enddo
    j3_ap = j3_ap
    d4_b  = -j1_ap + (-j2_ap) + j3_ap + j3_ap
    d8_b  = -j3_ap + (-j1_ap) + j2_ap + j2_ap
    d5_b  = -j3_ap + (-j2_ap) + j1_ap + j1_ap
    do g_i_ = 1, g_p_
        g_w(g_i_) = d4_b * g_j3_ap(g_i_) + d8_b * g_j2_ap(g_i_) + d5_b * g_j1_ap(g_i_)
    enddo
    w    = j1_ap * j1_ap + j2_ap * j2_ap + j3_ap * j3_ap - j1_ap * j2_ap - j2_ap * j3_ap - j3_ap * j1_ap
    d7_v = sqrt(w)
    !        if ( w .gt. 0.0d0 ) then
    d1_p = 1.0e0_wp / (2.0e0_wp *  d7_v)
    !        else
    !           call ehufDO (9,w, d7_v, d1_p,
    !     +g_ehfid,
    !     +746)
    !        endif
    do g_i_ = 1, g_p_
        g_leps_ap(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ap(g_i_) + g_q2_ap(g_i_) + g_q1_ap(g_i_)
    enddo
    leps_ap = q1_ap + q2_ap + q3_ap - d7_v + 5.77096e0_wp
    do g_i_ = 1, g_p_
        g_cfu22(g_i_) = -g_leps_ap(g_i_) + g_leps_ac(g_i_)
    enddo
    cfu22 = leps_ac - leps_ap
end subroutine lepscfu22
