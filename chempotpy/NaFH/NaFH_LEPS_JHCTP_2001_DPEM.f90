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

!> @brief Implement NaFH_LEPS_JHCTP_2001, a diabatic potential matrix for NaFH.
module chempotpy_nafh_leps_jhctp_2001_dpem_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NaFH_LEPS_JHCTP_2001 diabatic potential matrix for NaFH.
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
end module chempotpy_nafh_leps_jhctp_2001_dpem_pes

!> @brief Expose NaFH_LEPS_JHCTP_2001 through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_nafh_leps_jhctp_2001_dpem_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 2
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: u(nstates,nstates)
    double precision, intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem

!> @brief Evaluate EvdRdX for the NaFH_LEPS_JHCTP_2001_DPEM potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    implicit none
    integer, parameter            :: first_atoms(3)  = [1, 2, 1]
    integer, parameter            :: second_atoms(3) = [2, 3, 3]
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!> @brief Evaluate diagonalize for the NaFH_LEPS_JHCTP_2001_DPEM potential.
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

!> @brief Evaluate pot for the NaFH_LEPS_JHCTP_2001_DPEM potential.
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
    implicit none
    double precision :: de(3,nt), e(nt), r(nt,3)
    integer          :: nsurf, nt
    if( nsurf .eq. 1 )then
        call pot11(r, e, de, nt)
    elseif( nsurf .eq. 2 )then
        call pot12(r, e, de, nt)
    else
        call pot22(r, e, de, nt)
    endif
end subroutine pot

!  U11 surface - lower diabatic surface
module chempotpy_nafh_leps_jhctp_2001_dpem_pot11_legacy_state_0db22763
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, de, nt )
        implicit none
        optional :: r, e, de, nt
        double precision :: abs1, ac(4), acut1(4), adip1(4), ads1, af1, af2, al2(4), al2n(4), aln0(4), aln0n1(4), aln0n2(4), &
            & amf, amh, amna, amu, anorm, ars1, asmall, axs1, bc(0:5), c2a(4), c2b(4), c2c(4), cevau, cfu11, chf(6), cnf(6), &
            & cnf1(5), cog, coul1, coul2, coul3, cplg2, csb, danormd1, danormd2, danormd3, dcd1, dcd2, dcd3, dcfu11(3), dcl1d1, &
            & dcl2d2, dcl3d3, dcog, dcogd1, dcogd2, dcogd3, dcsbd1, dcsbd2, dcsbd3, de(:,:), dehf, denf, dex1d1, dex2d2, dex3d3, &
            & dip1(4), dipn(4), dipnn(4), drnag2d1, drnag2d2, drnag2d3, ds1d1, ds2d2, ds3d3, dsurfd1, dsurfd2, dsurfd3, dt1d1, &
            & dt2d2, dt3d3, dwd1, dwd2, dwd3, e(:), exch1, exch2, exch3, g_r(3,1,3), hlp, hlp1, hlp11, hlp2, hlp22, hlp23, hlp3, &
            & r(:,:), r1s, r2s, r3s, rac2, rcut1(4), rdip1(4), rdum(3), rehf, renf, rn0n(4), rnag2, rootw, s1, s2, s3, surf, t1, &
            & t2, t3, una2p, w, y(4,3,4), y2, y3
        integer :: i, j, l, l1, l2, nt
        integer :: g_myrank, g_nprocs
        abs1      = 0.d0
        ac        = 0.d0
        acut1     = 0.d0
        adip1     = 0.d0
        ads1      = 0.d0
        af1       = 0.d0
        af2       = 0.d0
        al2       = 0.d0
        al2n      = 0.d0
        aln0      = 0.d0
        aln0n1    = 0.d0
        aln0n2    = 0.d0
        amf       = 0.d0
        amh       = 0.d0
        amna      = 0.d0
        amu       = 0.d0
        anorm     = 0.d0
        ars1      = 0.d0
        asmall    = 0.d0
        axs1      = 0.d0
        bc        = 0.d0
        c2a       = 0.d0
        c2b       = 0.d0
        c2c       = 0.d0
        cevau     = 0.d0
        cfu11     = 0.d0
        chf       = 0.d0
        cnf       = 0.d0
        cnf1      = 0.d0
        cog       = 0.d0
        coul1     = 0.d0
        coul2     = 0.d0
        coul3     = 0.d0
        cplg2     = 0.d0
        csb       = 0.d0
        danormd1  = 0.d0
        danormd2  = 0.d0
        danormd3  = 0.d0
        dcd1      = 0.d0
        dcd2      = 0.d0
        dcd3      = 0.d0
        dcfu11    = 0.d0
        dcl1d1    = 0.d0
        dcl2d2    = 0.d0
        dcl3d3    = 0.d0
        dcog      = 0.d0
        dcogd1    = 0.d0
        dcogd2    = 0.d0
        dcogd3    = 0.d0
        dcsbd1    = 0.d0
        dcsbd2    = 0.d0
        dcsbd3    = 0.d0
        dehf      = 0.d0
        denf      = 0.d0
        dex1d1    = 0.d0
        dex2d2    = 0.d0
        dex3d3    = 0.d0
        dip1      = 0.d0
        dipn      = 0.d0
        dipnn     = 0.d0
        drnag2d1  = 0.d0
        drnag2d2  = 0.d0
        drnag2d3  = 0.d0
        ds1d1     = 0.d0
        ds2d2     = 0.d0
        ds3d3     = 0.d0
        dsurfd1   = 0.d0
        dsurfd2   = 0.d0
        dsurfd3   = 0.d0
        dt1d1     = 0.d0
        dt2d2     = 0.d0
        dt3d3     = 0.d0
        dwd1      = 0.d0
        dwd2      = 0.d0
        dwd3      = 0.d0
        exch1     = 0.d0
        exch2     = 0.d0
        exch3     = 0.d0
        g_r       = 0.d0
        hlp       = 0.d0
        hlp1      = 0.d0
        hlp11     = 0.d0
        hlp2      = 0.d0
        hlp22     = 0.d0
        hlp23     = 0.d0
        hlp3      = 0.d0
        r1s       = 0.d0
        r2s       = 0.d0
        r3s       = 0.d0
        rac2      = 0.d0
        rcut1     = 0.d0
        rdip1     = 0.d0
        rdum      = 0.d0
        rehf      = 0.d0
        renf      = 0.d0
        rn0n      = 0.d0
        rnag2     = 0.d0
        rootw     = 0.d0
        s1        = 0.d0
        s2        = 0.d0
        s3        = 0.d0
        surf      = 0.d0
        t1        = 0.d0
        t2        = 0.d0
        t3        = 0.d0
        una2p     = 0.d0
        w         = 0.d0
        y         = 0.d0
        y2        = 0.d0
        y3        = 0.d0
        y(1,1,1)  = 0.00000d0
        y(2,1,1)  = 0.00000d0
        y(1,1,2)  = 0.00000d0
        y(2,1,2)  = 0.00000d0
        y(1,1,3)  = 0.00000d0
        y(2,1,3)  = 0.00000d0
        y(1,1,4)  = 0.00000d0
        y(2,1,4)  = 0.00000d0
        y(3,1,1)  = 6.6796d0
        y(4,1,1)  = 2.7839d0
        y(3,1,2)  = 6.6796d0
        y(4,1,2)  = 2.7839d0
        y(3,1,3)  = 7.1360d0
        y(4,1,3)  = 1.0376d0
        y(3,1,4)  = 8.8737d0
        y(4,1,4)  = 1.9627d0
        y(1,2,1)  = 6.0459d+1
        y(2,2,1)  = 3.7503d0
        y(1,2,2)  = 6.0459d+1
        y(2,2,2)  = 3.7503d0
        y(1,2,3)  = 6.0459d+1
        y(2,2,3)  = 3.7503d0
        y(1,2,4)  = 8.5430d0
        y(2,2,4)  = 9.00856d-1
        y(3,2,1)  = 8.4874d0
        y(4,2,1)  = 1.0916d+1
        y(3,2,2)  = 8.4874d0
        y(4,2,2)  = 9.9850d0
        y(3,2,3)  = 8.4874d0
        y(4,2,3)  = 5.5898d0
        y(3,2,4)  = 9.5899d0
        y(4,2,4)  = 1.1615d+1
        y(1,3,1)  = 7.3251d0
        y(2,3,1)  = 2.8774d0
        y(1,3,2)  = 7.3251d0
        y(2,3,2)  = 3.6094d0
        y(1,3,3)  = 7.3251d0
        y(2,3,3)  = 3.9288d0
        y(1,3,4)  = 7.5811d0
        y(2,3,4)  = 4.4272d0
        y(3,3,1)  = 1.9300d+3
        y(4,3,1)  = 2.7091d0
        y(3,3,2)  = 3.8969d+3
        y(4,3,2)  = 2.7591d0
        y(3,3,3)  = 5.9130d+3
        y(4,3,3)  = 2.6881d0
        y(3,3,4)  = 5.9007d+3
        y(4,3,4)  = 3.8929d0
        c2a(1)    = 8.5742d0
        c2b(1)    = 1.5953d0
        c2c(1)    = 3.3745d-1
        c2a(2)    = 9.8509d0
        c2b(2)    = 6.6359d0
        c2c(2)    = 3.0396d-1
        c2a(3)    = 9.5937d0
        c2b(3)    = 1.3441d+1
        c2c(3)    = 2.8812d-1
        c2a(4)    = 1.1386d0
        c2b(4)    = 8.1112d0
        c2c(4)    = 2.6061d0
        dip1(1)   = 1.9651d0
        rdip1(1)  = 8.3109d-1
        dip1(2)   = 1.9651d0
        rdip1(2)  = 8.3109d-1
        dip1(3)   = 2.1803d0
        rdip1(3)  = 7.1647d-1
        dip1(4)   = 2.2483d0
        rdip1(4)  = 9.7128d-1
        adip1(1)  = 4.8674d-2
        adip1(2)  = 5.0603d-2
        adip1(3)  = 4.8499d-2
        adip1(4)  = 8.5924d-2
        acut1(1)  = 6.3761d0
        rcut1(1)  = 5.2007d0
        acut1(2)  = 6.3761d0
        rcut1(2)  = 3.6585d0
        acut1(3)  = 6.0812d0
        rcut1(3)  = 3.4684d0
        acut1(4)  = 5.0516d0
        rcut1(4)  = 1.6434d0
        dipn(1)   = 1.5754d-1
        al2(1)    = 8.4809d-2
        dipn(2)   = 1.2274d-1
        al2(2)    = 8.4809d-2
        dipn(3)   = 4.1823d-2
        al2(3)    = 9.4867d-2
        dipn(4)   = 3.4350d-1
        al2(4)    = 6.5542d-2
        aln0(1)   = 3.1971d-2
        aln0(2)   = 3.5817d-2
        aln0(3)   = 1.6686d-1
        aln0(4)   = 5.9961d-2
        dipnn(1)  = 0.0000d0
        al2n(1)   = 0.0000d-2
        dipnn(2)  = 3.7501d-3
        al2n(2)   = 9.2376d0
        dipnn(3)  = 7.6323d-3
        al2n(3)   = 9.3515d-2
        dipnn(4)  = 0.0000d0
        al2n(4)   = .0000d-2
        rn0n(1)   = 0.0d0
        rn0n(2)   = 6.6001d0
        rn0n(3)   = 5.9949d0
        rn0n(4)   = 0.0d0
        aln0n1(1) = 0.0000d-2
        aln0n2(1) = 0.0000d-2
        aln0n1(2) = 3.2545d-2
        aln0n2(2) = 4.9001d-1
        aln0n1(3) = 8.5914d-2
        aln0n2(3) = 2.1382d0
        aln0n1(4) = 0.0000d-2
        aln0n2(4) = 0.0000d-2
        ac        = [ 4.6913d0, 5.4684d-1, 7.5932d-1, 8.0000d-1]
        amh       = 1.007825d0
        amf       = 18.998403d0
        amna      = 22.989767d0
        amu       = 1822.8885d0
        denf      = 4.3830525d0
        renf      = 3.6395d0
        cnf       = [ 0.168364d0, 2.29333d0, -0.848620d-1, -0.109707d-1, 0.246582d-2, 0.377012d0]
        cnf1      = [ 0.620115d0, 0.694459d-1, 0.757570d0, 0.283055d0, 2.69998d0]
        dehf      = 5.77096d0
        rehf      = 1.7328d0
        chf       = [ 0.441258d0, 3.88056d0, -1.27110d0, -1.37766d0, 0.168186d0, 2.07230d0]
        bc        = [ -1.d0, 0.d0, 1.d0, 2.d0, 3.0d0, 5.d0]
        ads1      = 2.91087d-3
        ars1      = 7.37707d0
        abs1      = 7.13223d-1
        una2p     = 2.0973375d0
        asmall    = 0.09d0
        cevau     = 0.036749308867692d0
        g_myrank  = 0
        g_nprocs  = 1
        rac2      = 1.d0/dsqrt(2.d0)
        af1       = amh/(amf+amh)
        af2       = amf/(amf+amh)
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
            drnag2d1 = 2.d0*af1*r(i,1)
            drnag2d2 = -2.d0*af1*af2*r(i,2)
            drnag2d3 = 2.d0*af2*r(i,3)
            !     Cosine of Na-F-H bond angle:
            hlp    = 1.d0/(r(i,2)*r(i,3))
            csb    = 0.5d0*(r2s+r3s-r1s)*hlp
            dcsbd1 = -2.d0*r(i,1)*hlp
            dcsbd2 = 2.d0*(1.d0/r(i,3)-csb/r(i,2))
            dcsbd3 = 2.d0*(1.d0/r(i,2)-csb/r(i,3))
            csb    = 2.d0*(1.d0+csb)
            if( csb .ge. 4.d0 ) csb = 3.99999d0
            l1 = idint(csb)+1
            l2 = l1+1
            if( l1 .eq. 4 ) l2=l1
            if( csb .le. bc(l1) ) csb = bc(l1)+1.d-13
            if( csb .ge. bc(l1+1) ) csb = bc(l1+1)-1.d-13
            axs1     = dexp(-abs1*(r(i,1)-ars1))
            s1       = ads1*axs1*(axs1-2.d0)
            ds1d1    = 2.d0*ads1*abs1*axs1*(1.d0-axs1)
            y2       = r(i,2)-rehf
            hlp1     = chf(1)*dexp(-chf(2)*y2)
            hlp2     = dexp(-chf(6)*y2)
            hlp22    = (-1.d0-chf(1)+chf(3)*y2+chf(4)*y2**2+chf(5)*y2**3)*hlp2
            s2       = dehf*(hlp1+hlp22)
            ds2d2    = dehf*(-chf(2)*hlp1 - chf(6)*hlp22 + (chf(3)+2.d0*chf(4)*y2+3.d0*chf(5)*y2**2)*hlp2)
            y3       = r(i,3)-renf
            hlp1     = dexp(-cnf1(3)*y3)
            hlp11    = (cnf1(1)+cnf1(2)*y3)*hlp1
            hlp2     = cnf1(4)*dexp(-cnf1(5)*y3)
            s3       = hlp11+hlp2
            ds3d3    = -cnf1(3)*hlp11+cnf1(2)*hlp1-cnf1(5)*hlp2
            e(i)     = 0.d0
            de(1,i)  = 0.d0
            de(2,i)  = 0.d0
            de(3,i)  = 0.d0
            anorm    = 0.d0
            danormd1 = 0.d0
            danormd2 = 0.d0
            danormd3 = 0.d0
            do l = l1, l2
                hlp1   = y(1,1,l)*dexp(-y(2,1,l)*r(i,1))
                hlp2   = y(3,1,l)*dexp(-y(4,1,l)*r(i,1))
                t1     = hlp1+hlp2+s1
                dt1d1  = -y(2,1,l)*hlp1-y(4,1,l)*hlp2+ds1d1
                hlp1   = y(1,2,l)*dexp(-y(2,2,l)*r(i,2))
                hlp2   = y(3,2,l)*dexp(-y(4,2,l)*r(i,2))
                t2     = hlp1+hlp2+s2
                dt2d2  = -y(2,2,l)*hlp1-y(4,2,l)*hlp2+ds2d2
                hlp1   = y(1,3,l)*dexp(-y(2,3,l)*r(i,3))
                hlp2   = y(3,3,l)*dexp(-y(4,3,l)*r(i,3))
                t3     = hlp1+hlp2+s3
                dt3d3  = -y(2,3,l)*hlp1-y(4,3,l)*hlp2+ds3d3
                coul1  = 0.5d0*(s1+t1)
                dcl1d1 = 0.5d0*(ds1d1+dt1d1)
                coul2  = 0.5d0*(s2+t2)
                dcl2d2 = 0.5d0*(ds2d2+dt2d2)
                coul3  = 0.5d0*(s3+t3)
                dcl3d3 = 0.5d0*(ds3d3+dt3d3)
                exch1  = 0.5d0*(s1-t1)
                dex1d1 = 0.5d0*(ds1d1-dt1d1)
                exch2  = 0.5d0*(s2-t2)
                dex2d2 = 0.5d0*(ds2d2-dt2d2)
                exch3  = 0.5d0*(s3-t3)
                dex3d3 = 0.5d0*(ds3d3-dt3d3)
                w      = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
                dwd1   = 2.d0*(2.d0*exch1-exch2-exch3)*dex1d1
                dwd2   = 2.d0*(2.d0*exch2-exch1-exch3)*dex2d2
                dwd3   = 2.d0*(2.d0*exch3-exch1-exch2)*dex3d3
                cplg2  = c2a(l)*dexp(-c2b(l)*w-c2c(l)*(r(i,1)+r(i,2)+r(i,3)))
                dcd1   = (-c2b(l)*dwd1 - c2c(l))*cplg2
                dcd2   = (-c2b(l)*dwd2 - c2c(l))*cplg2
                dcd3   = (-c2b(l)*dwd3 - c2c(l))*cplg2
                rootw  = rac2*dsqrt(w+cplg2*cplg2)
                surf   = (coul1+coul2+coul3+dehf-rootw)
                ! check for zero by zero division (a case when all three atoms are asymp
                ! separated)
                if( (dwd1+2.d0*cplg2*dcd1) .eq. 0.0d0 )then
                    dsurfd1 = dcl1d1
                else
                    dsurfd1 = (dcl1d1-0.25d0*(dwd1+2.d0*cplg2*dcd1)/rootw)
                endif
                if( (dwd2+2.d0*cplg2*dcd2) .eq. 0.0d0 )then
                    dsurfd2 = dcl2d2
                else
                    dsurfd2 = (dcl2d2-0.25d0*(dwd2+2.d0*cplg2*dcd2)/rootw)
                endif
                if( (dwd3+2.d0*cplg2*dcd3) .eq. 0.0d0 )then
                    dsurfd3 = dcl3d3
                else
                    dsurfd3 = (dcl3d3-0.25d0*(dwd3+2.d0*cplg2*dcd3)/rootw)
                endif
                hlp1    = adip1(l)*(rnag2-rdip1(l))
                hlp2    = acut1(l)*(r(i,2)-rcut1(l))
                hlp11   = 0.5d0/(dcosh(hlp1)**2)
                hlp22   = 0.5d0/(dcosh(hlp2)**2)
                hlp1    = 0.5d0*(dtanh(hlp1)-1.d0)
                hlp2    = 0.5d0*(1.d0-dtanh(hlp2))
                surf    = surf+dip1(l)*hlp1*hlp2
                dsurfd1 = dsurfd1+dip1(l)*adip1(l)*drnag2d1*hlp11*hlp2
                dsurfd2 = dsurfd2+dip1(l)*adip1(l)*drnag2d2*hlp11*hlp2 -dip1(l)*hlp1*acut1(l)*hlp22
                dsurfd3 = dsurfd3+dip1(l)*adip1(l)*drnag2d3*hlp11*hlp2
                hlp1    = al2(l)*(r(i,2)-rehf)
                hlp2    = aln0(l)*rnag2
                hlp11   = dtanh(hlp1)
                hlp22   = dtanh(hlp2)
                hlp1    = 1.d0/dcosh(hlp1)
                hlp2    = 1.d0/dcosh(hlp2)
                surf    = surf - dipn(l)*hlp1*hlp2
                dsurfd1 = dsurfd1 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*drnag2d1
                dsurfd2 = dsurfd2 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*drnag2d2 + dipn(l)*hlp2*hlp1*hlp11*al2(l)
                dsurfd3 = dsurfd3 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*drnag2d3
                hlp1    = al2n(l)*(r(i,2)-rehf)
                hlp2    = aln0n1(l)*(rnag2-rn0n(l)**2)
                ! try to prevent an overflow
                if( hlp2 .lt. 709.0d0 )then
                    hlp2 = dexp(hlp2)
                else
                    hlp2 = 1.0d308
                endif
                hlp3    = dexp(-aln0n2(l)*(rnag2-rn0n(l)**2))
                hlp23   = 1.d0/(hlp2+hlp3)**2
                hlp11   = dtanh(hlp1)
                hlp1    = 1.d0/dcosh(hlp1)
                surf    = surf - 2.d0*dipnn(l)*hlp1/(hlp2+hlp3)
                dsurfd1 = dsurfd1 + 2.d0*dipnn(l)*hlp1 *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*drnag2d1
                dsurfd2=dsurfd2 + 2.d0*dipnn(l)*hlp1 *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*drnag2d2 + 2.d0* dipnn(l)* &
                    & hlp1* hlp11*al2n(l)/(hlp2+hlp3)
                dsurfd3  = dsurfd3 + 2.d0*dipnn(l)*hlp1 *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*drnag2d3
                hlp      = 1.d0/( (1.d0-(csb-bc(l))/(bc(l-1)-bc(l))) * (1.d0-(csb-bc(l))/(bc(l+1)-bc(l))) )
                hlp1     = (csb-bc(l))**2*hlp
                cog      = dexp(-ac(l)*hlp1)
                dcog     = -cog*ac(l)*(2.d0*(csb-bc(l))*hlp +hlp1*(1.d0/(bc(l-1)-csb)+1.d0/(bc(l+1)-csb)))
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
            hlp     = 1.d0/anorm
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
end module chempotpy_nafh_leps_jhctp_2001_dpem_pot11_legacy_state_0db22763

!> @brief Compatibility interface for pot11.
subroutine pot11( r, e, de, nt )
    use chempotpy_nafh_leps_jhctp_2001_dpem_pot11_legacy_state_0db22763, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    double precision :: de(3,nt)
    integer          :: nt
    call dispatch(r=r, e=e, de=de, nt=nt)
end subroutine pot11

!  U12 surface - coupling surface
module chempotpy_nafh_leps_jhctp_2001_dpem_pot12_legacy_state_0db22763
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, de, nt )
        implicit none
        optional :: r, e, de, nt
        double precision :: abs1, ads1, af1, af2, ahf(3), ahft, amf, amh, amna, amu, anaf(3), anaft, anah(3), anaht, ars1, &
            & asmall, bc(0:5), ccnf(4), ccnh(4), ccnha(4), ccnhi(4), ccnhi1(4), cevau, cf_p, cf_r0, chf(6), cnf(6), cnf1(5), &
            & cstsq, dahftd1, dahftd2, dahftd3, danaftd1, danaftd2, danaftd3, danahtd1, danahtd2, danahtd3, dccnhd2(4), &
            & dccnhd3(4), dcstsqd1, dcstsqd2, dcstsqd3, de(:,:), dehf, denf, dr13td1, dr13td2, dr13td3, dr21td1, dr21td2, &
            & dr21td3, dr23td1, dr23td2, dr23td3, drnag2d1, drnag2d2, drnag2d3, dsw1d1, dsw1d2, dsw1d3, dsw21d1, dsw21d2, &
            & dsw21d3, dsw23d1, dsw23d2, dsw23d3, dswit1, dswit2, dswit3, dunafd3, dunahd1, dunahd2, dunahd3, e(:), hlp, hlp1, &
            & hlp2, r(:,:), r13(3), r13t, r1s, r21(3), r21t, r23(3), r23t, r2s, r3s, rehf, renf, rnag2, sw1, sw21, sw23, swit, &
            & una2p, unaf, unah, xpnh, xpnh1
        integer :: i, nt
        integer :: g_myrank, g_nprocs
        abs1     = 0.d0
        ads1     = 0.d0
        af1      = 0.d0
        af2      = 0.d0
        ahf      = 0.d0
        ahft     = 0.d0
        amf      = 0.d0
        amh      = 0.d0
        amna     = 0.d0
        amu      = 0.d0
        anaf     = 0.d0
        anaft    = 0.d0
        anah     = 0.d0
        anaht    = 0.d0
        ars1     = 0.d0
        asmall   = 0.d0
        bc       = 0.d0
        ccnf     = 0.d0
        ccnh     = 0.d0
        ccnha    = 0.d0
        ccnhi    = 0.d0
        ccnhi1   = 0.d0
        cevau    = 0.d0
        cf_p     = 0.d0
        cf_r0    = 0.d0
        chf      = 0.d0
        cnf      = 0.d0
        cnf1     = 0.d0
        cstsq    = 0.d0
        dahftd1  = 0.d0
        dahftd2  = 0.d0
        dahftd3  = 0.d0
        danaftd1 = 0.d0
        danaftd2 = 0.d0
        danaftd3 = 0.d0
        danahtd1 = 0.d0
        danahtd2 = 0.d0
        danahtd3 = 0.d0
        dccnhd2  = 0.d0
        dccnhd3  = 0.d0
        dcstsqd1 = 0.d0
        dcstsqd2 = 0.d0
        dcstsqd3 = 0.d0
        dehf     = 0.d0
        denf     = 0.d0
        dr13td1  = 0.d0
        dr13td2  = 0.d0
        dr13td3  = 0.d0
        dr21td1  = 0.d0
        dr21td2  = 0.d0
        dr21td3  = 0.d0
        dr23td1  = 0.d0
        dr23td2  = 0.d0
        dr23td3  = 0.d0
        drnag2d1 = 0.d0
        drnag2d2 = 0.d0
        drnag2d3 = 0.d0
        dsw1d1   = 0.d0
        dsw1d2   = 0.d0
        dsw1d3   = 0.d0
        dsw21d1  = 0.d0
        dsw21d2  = 0.d0
        dsw21d3  = 0.d0
        dsw23d1  = 0.d0
        dsw23d2  = 0.d0
        dsw23d3  = 0.d0
        dswit1   = 0.d0
        dswit2   = 0.d0
        dswit3   = 0.d0
        dunafd3  = 0.d0
        dunahd1  = 0.d0
        dunahd2  = 0.d0
        dunahd3  = 0.d0
        hlp      = 0.d0
        hlp1     = 0.d0
        hlp2     = 0.d0
        r13      = 0.d0
        r13t     = 0.d0
        r1s      = 0.d0
        r21      = 0.d0
        r21t     = 0.d0
        r23      = 0.d0
        r23t     = 0.d0
        r2s      = 0.d0
        r3s      = 0.d0
        rehf     = 0.d0
        renf     = 0.d0
        rnag2    = 0.d0
        sw1      = 0.d0
        sw21     = 0.d0
        sw23     = 0.d0
        swit     = 0.d0
        una2p    = 0.d0
        unaf     = 0.d0
        unah     = 0.d0
        xpnh     = 0.d0
        xpnh1    = 0.d0
        amh      = 1.007825d0
        amf      = 18.998403d0
        amna     = 22.989767d0
        amu      = 1822.8885d0
        denf     = 4.3830525d0
        renf     = 3.6395d0
        cnf      = [ 0.168364d0, 2.29333d0, -0.848620d-1, -0.109707d-1, 0.246582d-2, 0.377012d0]
        cnf1     = [ 0.620115d0, 0.694459d-1, 0.757570d0, 0.283055d0, 2.69998d0]
        dehf     = 5.77096d0
        rehf     = 1.7328d0
        chf      = [ 0.441258d0, 3.88056d0, -1.27110d0, -1.37766d0, 0.168186d0, 2.07230d0]
        bc       = [ -1.d0, 0.d0, 1.d0, 2.d0, 3.0d0, 5.d0]
        ads1     = 2.91087d-3
        ars1     = 7.37707d0
        abs1     = 7.13223d-1
        una2p    = 2.0973375d0
        asmall   = 0.09d0
        cevau    = 0.036749308867692d0
        xpnh     = 5.6906d-2
        xpnh1    = 2.3185d-1
        ahf      = [ 1.7177d0, 2.2707d0, 1.2248d-1]
        anaf     = [ 6.7786d-1, 9.8809d-1, -1.0541d-1]
        anah     = [ 6.5680d-1, 1.0141d0, -5.1056d-2]
        r21      = [ 3.5114d0, 1.3540d0, 3.2267d-1]
        r23      = [ -3.3006d-1, 2.0804d-1, 2.4499d-1]
        r13      = [ 1.2431d0, 9.9188d-1, 1.1051d-1]
        ccnf     = [ 3.9960d0, 5.4730d-01, -1.1428d0, 2.0440d-01]
        ccnha    = [ 1.00d0, 0.80d0, -2.67d0, 0.456d0]
        ccnhi    = [ 1.0897d0, -6.5205d-1, 7.7968d-1, -1.2933d-2]
        ccnhi1   = [ -6.1073d-1, -1.7425d-1, -8.3877d-1, -3.5361d-1]
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
            hlp1       = dexp(-xpnh*r(i,2)**2)
            hlp2       = dexp(-xpnh1*r(i,3)**2)
            ccnh(1)    = ccnha(1)+ccnhi(1)*hlp1+ccnhi1(1)*hlp2
            ccnh(2)    = ccnha(2)+ccnhi(2)*hlp1+ccnhi1(2)*hlp2
            ccnh(3)    = ccnha(3)+ccnhi(3)*hlp1+ccnhi1(3)*hlp2
            ccnh(4)    = ccnha(4)+ccnhi(4)*hlp1+ccnhi1(4)*hlp2
            hlp1       = - 2.d0*xpnh*r(i,2)*hlp1
            hlp2       = - 2.d0*xpnh1*r(i,3)*hlp2
            dccnhd2(1) = ccnhi(1)*hlp1
            dccnhd3(1) = ccnhi1(1)*hlp2
            dccnhd2(2) = ccnhi(2)*hlp1
            dccnhd3(2) = ccnhi1(2)*hlp2
            dccnhd2(3) = ccnhi(3)*hlp1
            dccnhd3(3) = ccnhi1(3)*hlp2
            dccnhd2(4) = ccnhi(4)*hlp1
            dccnhd3(4) = ccnhi1(4)*hlp2
            hlp        = ccnf(3)/(1.d0+ccnf(4)*r(i,3))
            unaf       = dexp(-ccnf(1)-ccnf(2)*r(i,3)- hlp*r(i,3))
            unaf       = unaf*27.2113961d0
            dunafd3    = -unaf*(ccnf(2)+hlp-ccnf(4)*hlp*hlp*r(i,3)/ccnf(3))
            hlp        = 1.d0/(1.d0+ccnh(4)*r(i,1))
            unah       = dexp(-ccnh(1)-ccnh(2)*r(i,1)- ccnh(3)*hlp*r(i,1))
            dunahd1    = -unah*(ccnh(2)+ccnh(3)*hlp -ccnh(4)*hlp*hlp*ccnh(3)*r(i,1))
            dunahd2    = -unah*(dccnhd2(1)+dccnhd2(2)*r(i,1) + dccnhd2(3)*r(i,1)*hlp-ccnh(3)*r(i,1)*r(i,1)* &
                & dccnhd2(4)* hlp*hlp)
            dunahd3    = -unah*(dccnhd3(1)+dccnhd3(2)*r(i,1) + dccnhd3(3)*r(i,1)*hlp-ccnh(3)*r(i,1)*r(i,1)* &
                & dccnhd3(4)* hlp*hlp)
            r1s      = r(i,1)*r(i,1)
            r2s      = r(i,2)*r(i,2)
            r3s      = r(i,3)*r(i,3)
            rnag2    = af1*r1s+af2*r3s-af1*af2*r2s
            drnag2d1 = 2.d0*af1*r(i,1)
            drnag2d2 = -2.d0*af1*af2*r(i,2)
            drnag2d3 = 2.d0*af2*r(i,3)
            hlp2     = rnag2+asmall
            hlp      = dsqrt(hlp2)
            cstsq    = 0.5d0*(r3s-r1s+(af2-af1)*r2s) /(r(i,2)*hlp)
            dcstsqd1 = -2.d0*r(i,1)/(r(i,2)*hlp) - cstsq*drnag2d1/hlp2
            dcstsqd2 = 2.d0*((af2-af1)/hlp-cstsq/r(i,2)) - cstsq*drnag2d2/hlp2
            dcstsqd3 = 2.d0*r(i,3)/(r(i,2)*hlp) - cstsq*drnag2d3/hlp2
            cstsq    = 2.0d0*(1.d0+cstsq)
            hlp      = r21(2)-r21(1)-r21(3)
            r21t     = r21(1)+hlp*cstsq +r21(3)*cstsq*cstsq
            hlp1     = hlp+2.d0*r21(3)*cstsq
            dr21td1  = hlp1*dcstsqd1
            dr21td2  = hlp1*dcstsqd2
            dr21td3  = hlp1*dcstsqd3
            hlp      = r23(2)-r23(1)-r23(3)
            r23t     = r23(1)+hlp*cstsq +r23(3)*cstsq*cstsq
            hlp1     = hlp+2.d0*r23(3)*cstsq
            dr23td1  = hlp1*dcstsqd1
            dr23td2  = hlp1*dcstsqd2
            dr23td3  = hlp1*dcstsqd3
            hlp      = r13(2)-r13(1)-r13(3)
            r13t     = r13(1)+hlp*cstsq +r13(3)*cstsq*cstsq
            hlp1     = hlp+2.d0*r13(3)*cstsq
            dr13td1  = hlp1*dcstsqd1
            dr13td2  = hlp1*dcstsqd2
            dr13td3  = hlp1*dcstsqd3
            hlp      = ahf(2)-ahf(1)-ahf(3)
            ahft     = ahf(1)+hlp*cstsq +ahf(3)*cstsq*cstsq
            hlp1     = hlp+2.d0*ahf(3)*cstsq
            dahftd1  = hlp1*dcstsqd1
            dahftd2  = hlp1*dcstsqd2
            dahftd3  = hlp1*dcstsqd3
            hlp      = anaf(2)-anaf(1)-anaf(3)
            anaft    = anaf(1)+hlp*cstsq +anaf(3)*cstsq*cstsq
            hlp1     = hlp+2.d0*anaf(3)*cstsq
            danaftd1 = hlp1*dcstsqd1
            danaftd2 = hlp1*dcstsqd2
            danaftd3 = hlp1*dcstsqd3
            hlp      = anah(2)-anah(1)-anah(3)
            anaht    = anah(1)+hlp*cstsq +anah(3)*cstsq*cstsq
            hlp1     = hlp+2.d0*anah(3)*cstsq
            danahtd1 = hlp1*dcstsqd1
            danahtd2 = hlp1*dcstsqd2
            danahtd3 = hlp1*dcstsqd3
            hlp      = anaht*r(i,1)-anaft*r(i,3)-r13t
            sw1      = 0.5d0*(1.d0+dtanh(hlp))
            hlp1     = 0.5d0/(dcosh(hlp)**2)
            dsw1d1   = hlp1*(anaht+danahtd1*r(i,1)-danaftd1*r(i,3)-dr13td1)
            dsw1d2   = hlp1*(danahtd2*r(i,1)-danaftd2*r(i,3)-dr13td2)
            dsw1d3   = hlp1*(danahtd3*r(i,1)-anaft-danaftd3*r(i,3)-dr13td3)
            hlp      = ahft*r(i,2)-anaht*r(i,1)-r21t
            sw21     = 0.5d0*(1.d0+dtanh(hlp))
            hlp1     = 0.5d0/(dcosh(hlp)**2)
            dsw21d1  = hlp1*(dahftd1*r(i,2)-anaht-danahtd1*r(i,1)-dr21td1)
            dsw21d2  = hlp1*(ahft+dahftd2*r(i,2)-danahtd2*r(i,1)-dr21td2)
            dsw21d3  = hlp1*(dahftd3*r(i,2)-danahtd3*r(i,1)-dr21td3)
            hlp      = ahft*r(i,2)-anaft*r(i,3)-r23t
            sw23     = 0.5d0*(1.d0+dtanh(hlp))
            hlp1     = 0.5d0/(dcosh(hlp)**2)
            dsw23d1  = hlp1*(dahftd1*r(i,2)-danaftd1*r(i,3)-dr23td1)
            dsw23d2  = hlp1*(ahft+dahftd2*r(i,2)-danaftd2*r(i,3)-dr23td2)
            dsw23d3  = hlp1*(dahftd3*r(i,2)-anaft-danaftd3*r(i,3)-dr23td3)
            e(i)     = unaf*sw1*sw23+unah*(1.d0-sw1)*sw21
            de(1,i) = unaf*dsw1d1*sw23+unaf*sw1*dsw23d1 + dunahd1*(1.d0-sw1)*sw21-unah*dsw1d1*sw21 + unah*(1.d0- sw1)* &
                & dsw21d1
            de(2,i) = unaf*dsw1d2*sw23+unaf*sw1*dsw23d2 + dunahd2*(1.d0-sw1)*sw21-unah*dsw1d2*sw21 + unah*(1.d0- sw1)* &
                & dsw21d2
            de(3,i) = dunafd3*sw1*sw23+unaf*dsw1d3*sw23+unaf*sw1*dsw23d3 + dunahd3*(1.d0-sw1)*sw21-unah*dsw1d3* sw21 + &
                & unah* (1.d0-sw1)*dsw21d3
            !     correction function
            cf_p    = 2.d0
            cf_r0   = 3.5d0
            swit    = 0.5d0*(1.d0+dtanh(cf_p*(cf_r0-r(i,2))))
            dswit1  = 0.d0
            dswit2  = -0.5d0*cf_p*1.d0/(dcosh(cf_p*(-r(i,2)+cf_r0))**2)
            dswit3  = 0.d0
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
end module chempotpy_nafh_leps_jhctp_2001_dpem_pot12_legacy_state_0db22763

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, de, nt )
    use chempotpy_nafh_leps_jhctp_2001_dpem_pot12_legacy_state_0db22763, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    double precision :: de(3,nt)
    integer          :: nt
    call dispatch(r=r, e=e, de=de, nt=nt)
end subroutine pot12

!  U22 surface - upper diabatic surface
module chempotpy_nafh_leps_jhctp_2001_dpem_pot22_legacy_state_0db22763
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, de, nt )
        implicit none
        optional :: r, e, de, nt
        double precision :: abs1, ac(4), acut2(4), acut2n(4), adip2(4), adip2n(4), ads1, af1, af2, al20, al20a, al20i(4), &
            & al20i1(4), al2a(4), aly(4), amf, amh, amna, amu, anorm, ars1, asmall, bc(0:5), c2a(4), c2b(4), c2c(4), cevau, &
            & cfu22, chf(6), cnf(6), cnf1(5), cog, coul1, coul2, coul3, cplg2, csb, ctanh, cy, dal20, danormd1, danormd2, &
            & danormd3, dcd1, dcd2, dcd3, dcfu22(3), dcl1d1, dcl2d1, dcl2d2, dcl2d3, dcl3d3, dcog, dcogd1, dcogd2, dcogd3, &
            & dcsbd1, dcsbd2, dcsbd3, dctanh, dcyd3, de(:,:), dehf, denf, dex1d1, dex2d1, dex2d2, dex2d3, dex3d3, dip2(4), &
            & dip2n(4), dr20, drnag2d1, drnag2d2, drnag2d3, ds1d1, ds2ad2, ds2d1, ds2d2, ds2d3, ds3d3, dsurfd1, dsurfd2, &
            & dsurfd3, dt1d1, dt2d1, dt2d2, dt2d3, dt3d3, dwd1, dwd2, dwd3, dy(4), dy1(4), dy2(4), e(:), exch1, exch2, exch3, &
            & g_r(3,1,3), hlp, hlp1, hlp11, hlp2, hlp22, hlp3, psi, r(:,:), r1s, r20, r20a, r20i(4), r2s, r3s, rac2, rcut2(4), &
            & rcut2n(4), rdip2(4), rdip2n(4), rdum(3), rehf, renf, rnag2, rootw, rr0a(4), rsi, ry(4), s1, s2, s2a, s3, surf, t1, &
            & t2, t3, una2p, w, y(4,3,4), y2, y3
        integer :: i, j, l, l1, l2, nt
        integer :: g_myrank, g_nprocs
        abs1      = 0.d0
        ac        = 0.d0
        acut2     = 0.d0
        acut2n    = 0.d0
        adip2     = 0.d0
        adip2n    = 0.d0
        ads1      = 0.d0
        af1       = 0.d0
        af2       = 0.d0
        al20      = 0.d0
        al20a     = 0.d0
        al20i     = 0.d0
        al20i1    = 0.d0
        al2a      = 0.d0
        aly       = 0.d0
        amf       = 0.d0
        amh       = 0.d0
        amna      = 0.d0
        amu       = 0.d0
        anorm     = 0.d0
        ars1      = 0.d0
        asmall    = 0.d0
        bc        = 0.d0
        c2a       = 0.d0
        c2b       = 0.d0
        c2c       = 0.d0
        cevau     = 0.d0
        cfu22     = 0.d0
        chf       = 0.d0
        cnf       = 0.d0
        cnf1      = 0.d0
        cog       = 0.d0
        coul1     = 0.d0
        coul2     = 0.d0
        coul3     = 0.d0
        cplg2     = 0.d0
        csb       = 0.d0
        ctanh     = 0.d0
        cy        = 0.d0
        dal20     = 0.d0
        danormd1  = 0.d0
        danormd2  = 0.d0
        danormd3  = 0.d0
        dcd1      = 0.d0
        dcd2      = 0.d0
        dcd3      = 0.d0
        dcfu22    = 0.d0
        dcl1d1    = 0.d0
        dcl2d1    = 0.d0
        dcl2d2    = 0.d0
        dcl2d3    = 0.d0
        dcl3d3    = 0.d0
        dcog      = 0.d0
        dcogd1    = 0.d0
        dcogd2    = 0.d0
        dcogd3    = 0.d0
        dcsbd1    = 0.d0
        dcsbd2    = 0.d0
        dcsbd3    = 0.d0
        dctanh    = 0.d0
        dcyd3     = 0.d0
        dehf      = 0.d0
        denf      = 0.d0
        dex1d1    = 0.d0
        dex2d1    = 0.d0
        dex2d2    = 0.d0
        dex2d3    = 0.d0
        dex3d3    = 0.d0
        dip2      = 0.d0
        dip2n     = 0.d0
        dr20      = 0.d0
        drnag2d1  = 0.d0
        drnag2d2  = 0.d0
        drnag2d3  = 0.d0
        ds1d1     = 0.d0
        ds2ad2    = 0.d0
        ds2d1     = 0.d0
        ds2d2     = 0.d0
        ds2d3     = 0.d0
        ds3d3     = 0.d0
        dsurfd1   = 0.d0
        dsurfd2   = 0.d0
        dsurfd3   = 0.d0
        dt1d1     = 0.d0
        dt2d1     = 0.d0
        dt2d2     = 0.d0
        dt2d3     = 0.d0
        dt3d3     = 0.d0
        dwd1      = 0.d0
        dwd2      = 0.d0
        dwd3      = 0.d0
        dy        = 0.d0
        dy1       = 0.d0
        dy2       = 0.d0
        exch1     = 0.d0
        exch2     = 0.d0
        exch3     = 0.d0
        g_r       = 0.d0
        hlp       = 0.d0
        hlp1      = 0.d0
        hlp11     = 0.d0
        hlp2      = 0.d0
        hlp22     = 0.d0
        hlp3      = 0.d0
        psi       = 0.d0
        r1s       = 0.d0
        r20       = 0.d0
        r20a      = 0.d0
        r20i      = 0.d0
        r2s       = 0.d0
        r3s       = 0.d0
        rac2      = 0.d0
        rcut2     = 0.d0
        rcut2n    = 0.d0
        rdip2     = 0.d0
        rdip2n    = 0.d0
        rdum      = 0.d0
        rehf      = 0.d0
        renf      = 0.d0
        rnag2     = 0.d0
        rootw     = 0.d0
        rr0a      = 0.d0
        rsi       = 0.d0
        ry        = 0.d0
        s1        = 0.d0
        s2        = 0.d0
        s2a       = 0.d0
        s3        = 0.d0
        surf      = 0.d0
        t1        = 0.d0
        t2        = 0.d0
        t3        = 0.d0
        una2p     = 0.d0
        w         = 0.d0
        y         = 0.d0
        y2        = 0.d0
        y3        = 0.d0
        ac        = [ 2.0018d0, 4.1908d-1, 5.4988d0, 3.0178d-1]
        y(1,1,1)  = 1.0407d+3
        y(2,1,1)  = 4.5911d0
        y(1,1,2)  = 1.1608d+3
        y(2,1,2)  = 3.8886d0
        y(1,1,3)  = 4.9492d+3
        y(2,1,3)  = 2.6147d0
        y(1,1,4)  = 8.3093d+3
        y(2,1,4)  = 3.2431d0
        y(3,1,1)  = 1.4941d0
        y(4,1,1)  = 4.5525d0
        y(3,1,2)  = 1.4941d0
        y(4,1,2)  = 4.5525d0
        y(3,1,3)  = 1.4941d0
        y(4,1,3)  = 4.5525d0
        y(3,1,4)  = 2.8361d0
        y(4,1,4)  = 4.4686d0
        y(1,2,1)  = 0.0000d0
        y(2,2,1)  = 0.0000d0
        y(1,2,2)  = 0.0000d0
        y(2,2,2)  = 0.0000d0
        y(1,2,3)  = 0.0000d0
        y(2,2,3)  = 0.0000d0
        y(1,2,4)  = 0.0000d0
        y(2,2,4)  = 0.0000d0
        y(3,2,1)  = 1.8076d+3
        y(4,2,1)  = 2.8815d0
        y(3,2,2)  = 1.8076d+3
        y(4,2,2)  = 2.8815d0
        y(3,2,3)  = 1.8076d+3
        y(4,2,3)  = 2.8815d0
        y(3,2,4)  = 2.3109d+3
        y(4,2,4)  = 3.1315d0
        y(1,3,1)  = 5.7330d0
        y(2,3,1)  = 3.2053d-1
        y(1,3,2)  = 5.7330d0
        y(2,3,2)  = 3.2053d-1
        y(1,3,3)  = 5.7330d0
        y(2,3,3)  = 3.2053d-1
        y(1,3,4)  = 5.0393d0
        y(2,3,4)  = 2.7567d-1
        y(3,3,1)  = 2.5774d+5
        y(4,3,1)  = 4.2267d0
        y(3,3,2)  = 2.5774d+5
        y(4,3,2)  = 4.2267d0
        y(3,3,3)  = 2.5774d+5
        y(4,3,3)  = 4.2267d0
        y(3,3,4)  = 1.6506d+5
        y(4,3,4)  = 3.9246d0
        c2a(1)    = 5.6406d0
        c2b(1)    = 7.7152d-2
        c2a(2)    = 8.2390d0
        c2b(2)    = 1.7133d-1
        c2a(3)    = 7.1504d0
        c2b(3)    = 1.1011d-1
        c2a(4)    = 8.0196d0
        c2b(4)    = 1.0017d-1
        c2c(1)    = 8.4595d-2
        c2c(2)    = 8.4595d-2
        c2c(3)    = 8.4595d-2
        c2c(4)    = 8.4411d-2
        r20i(1)   = 2.9013d0
        al20i(1)  = 1.7119d0
        r20i(2)   = 3.1196d0
        al20i(2)  = 9.4849d-1
        r20i(3)   = 3.0995d0
        al20i(3)  = 0.0d0
        r20i(4)   = 2.6572d0
        al20i(4)  = 2.3214d0
        r20a      = 3.1305d0
        al20i1(1) = 3.5428d0
        al20a     = 3.024d+1
        al20i1(2) = 6.7257d0
        psi       = 0.7d0
        al20i1(3) = 3.0240d+1
        rsi       = 22.2d0
        al20i1(4) = 9.6978d0
        al2a(1)   = 5.8088d-1
        rr0a(1)   = 9.5543d0
        al2a(2)   = 5.8088d-1
        rr0a(2)   = 9.5543d0
        al2a(3)   = 5.8088d-1
        rr0a(3)   = 9.5543d0
        al2a(4)   = 8.1525d-1
        rr0a(4)   = 8.2500d0
        dip2(1)   = 2.6475d-1
        dip2(2)   = 2.6475d-1
        dip2(3)   = 2.6475d-1
        dip2(4)   = 5.0336d-1
        adip2(1)  = 3.2317d-1
        rdip2(1)  = 3.7260d0
        adip2(2)  = 1.8513d-1
        rdip2(2)  = 4.5515d0
        adip2(3)  = 4.3999d-2
        rdip2(3)  = 2.1191d-1
        adip2(4)  = 2.3003d-1
        rdip2(4)  = 5.0160d0
        acut2(1)  = 2.0799d0
        rcut2(1)  = 1.2496d0
        acut2(2)  = 1.4204d0
        rcut2(2)  = 1.0817d0
        acut2(3)  = 7.9086d-1
        rcut2(3)  = 1.0575d0
        acut2(4)  = 2.8447d0
        rcut2(4)  = 3.0015d0
        dip2n(1)  = 2.2091d-1
        adip2n(1) = 1.0288d-1
        dip2n(2)  = 2.2091d-1
        adip2n(2) = 1.0288d-1
        dip2n(3)  = 2.2091d-1
        adip2n(3) = 1.0288d-1
        dip2n(4)  = 4.5301d-1
        adip2n(4) = 1.1009d-1
        rdip2n(1) = 5.1554d0
        rdip2n(2) = 5.1554d0
        rdip2n(3) = 5.1554d0
        rdip2n(4) = 6.7003d0
        acut2n(1) = 3.5607d0
        rcut2n(1) = 3.0661d0
        acut2n(2) = 4.3956d0
        rcut2n(2) = 3.0661d0
        acut2n(3) = 4.1924d0
        rcut2n(3) = 3.0661d0
        acut2n(4) = 3.0765d0
        rcut2n(4) = 2.8845d0
        ry(1)     = 1.4328d+1
        aly(1)    = 4.9044d-1
        ry(2)     = 1.2504d+1
        aly(2)    = 6.4718d-1
        ry(3)     = 1.2102d+1
        aly(3)    = 2.5833d-1
        ry(4)     = 7.9381d0
        aly(4)    = 6.4853d-1
        dy(1)     = 7.5126d-1
        dy1(1)    = 0.d0
        dy2(1)    = 0.d0
        dy(2)     = 2.0210d-1
        dy1(2)    = 0.d0
        dy2(2)    = 0.d0
        dy(3)     = 1.1749d-1
        dy1(3)    = 0.d0
        dy2(3)    = 0.d0
        dy(4)     = 8.1959d-2
        dy1(4)    = -5.03d+2
        dy2(4)    = -0.25d0
        amh       = 1.007825d0
        amf       = 18.998403d0
        amna      = 22.989767d0
        amu       = 1822.8885d0
        denf      = 4.3830525d0
        renf      = 3.6395d0
        cnf       = [ 0.168364d0, 2.29333d0, -0.848620d-1, -0.109707d-1, 0.246582d-2, 0.377012d0]
        cnf1      = [ 0.620115d0, 0.694459d-1, 0.757570d0, 0.283055d0, 2.69998d0]
        dehf      = 5.77096d0
        rehf      = 1.7328d0
        chf       = [ 0.441258d0, 3.88056d0, -1.27110d0, -1.37766d0, 0.168186d0, 2.07230d0]
        bc        = [ -1.d0, 0.d0, 1.d0, 2.d0, 3.0d0, 5.d0]
        ads1      = 2.91087d-3
        ars1      = 7.37707d0
        abs1      = 7.13223d-1
        una2p     = 2.0973375d0
        asmall    = 0.09d0
        cevau     = 0.036749308867692d0
        g_myrank  = 0
        g_nprocs  = 1
        rac2      = 1.d0/dsqrt(2.d0)
        af1       = amh/(amf+amh)
        af2       = amf/(amf+amh)
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
            drnag2d1 = 2.d0*af1*r(i,1)
            drnag2d2 = -2.d0*af1*af2*r(i,2)
            drnag2d3 = 2.d0*af2*r(i,3)
            !     Cosine of Na-F-H bond angle:
            hlp    = 1.d0/(r(i,2)*r(i,3))
            csb    = 0.5d0*(r2s+r3s-r1s)*hlp
            dcsbd1 = -2.d0*r(i,1)*hlp
            dcsbd2 = 2.d0*(1.d0/r(i,3)-csb/r(i,2))
            dcsbd3 = 2.d0*(1.d0/r(i,2)-csb/r(i,3))
            csb    = 2.d0*(1.d0+csb)
            if( csb .ge. 4.d0 ) csb = 3.99999d0
            l1 = idint(csb)+1
            l2 = l1+1
            if( l1 .eq. 4 ) l2=l1
            if( csb .le. bc(l1) ) csb = bc(l1)+1.d-13
            if( csb .ge. bc(l1+1) ) csb = bc(l1+1)-1.d-13
            s1       = 0.d0
            ds1d1    = 0.d0
            y2       = r(i,2)-rehf
            hlp1     = chf(1)*dexp(-chf(2)*y2)
            hlp2     = dexp(-chf(6)*y2)
            hlp22    = (-1.d0-chf(1)+chf(3)*y2+chf(4)*y2**2+chf(5)*y2**3)*hlp2
            s2a      = dehf*(hlp1+hlp22) + una2p
            ds2ad2   = dehf*(-chf(2)*hlp1-chf(6)*hlp22 + (chf(3)+2.d0*chf(4)*y2+3.d0*chf(5)*y2**2)*hlp2)
            y3       = r(i,3)-renf
            hlp1     = cnf(1)*dexp(-cnf(2)*y3)
            hlp2     = dexp(-cnf(6)*y3)
            hlp22    = (-1.d0-cnf(1)+cnf(3)*y3+cnf(4)*y3*y3+cnf(5)*y3**3)*hlp2
            s3       = denf*(hlp1+hlp22)
            ds3d3    = denf*(-cnf(2)*hlp1-cnf(6)*hlp22 + (cnf(3)+2.d0*cnf(4)*y3+3.d0*cnf(5)*y3**2)*hlp2)
            e(i)     = 0.d0
            de(1,i)  = 0.d0
            de(2,i)  = 0.d0
            de(3,i)  = 0.d0
            anorm    = 0.d0
            danormd1 = 0.d0
            danormd2 = 0.d0
            danormd3 = 0.d0
            do l = l1, l2
                hlp1  = y(1,1,l)*dexp(-y(2,1,l)*r(i,1))
                hlp2  = y(3,1,l)*dexp(-y(4,1,l)*r(i,1))
                t1    = hlp1+hlp2+s1
                dt1d1 = -y(2,1,l)*hlp1-y(4,1,l)*hlp2+ds1d1
                hlp   = al2a(l)*(r(i,3)+r(i,1)-rr0a(l))
                hlp1  = al2a(l)/(dcosh(hlp)**2)
                hlp   = 0.5d0*(1.d0+dtanh(hlp))
                hlp2  = r(i,3)+r(i,1)-rsi
                if( hlp2 .gt. 0.d0 )then
                    !     write(6,*)'R1+R3=',(r(i,1)+r(i,3))
                    ctanh  = dexp(-psi/hlp2)
                    dctanh = psi*ctanh/(r(i,3)+r(i,1)-rsi)**2
                else
                    ctanh  = 0.d0
                    dctanh = 0.d0
                endif
                r20      = r20i(l)+(r20a-r20i(l))*hlp
                dr20     = 0.5d0*(r20a-r20i(l))*hlp1
                al20     = al20i(l)+(al20i1(l)-al20i(l))*hlp + (al20a-al20i1(l))*ctanh
                dal20    = 0.5d0*(al20i1(l)-al20i(l))*hlp1 + (al20a-al20i1(l))*dctanh
                hlp1     = y(1,2,l)*dexp(-y(2,2,l)*r(i,2))
                hlp2     = dexp(-(y(4,2,l)+dy2(l)*ctanh)*r(i,2))
                hlp3     = (y(3,2,l)+dy1(l)*ctanh)*hlp2
                t2       = hlp1+hlp3
                dt2d1    = dy1(l)*dctanh*hlp2 - dy2(l)*dctanh*r(i,2)*hlp3
                dt2d2    = -y(2,2,l)*hlp1-(y(4,2,l)+dy2(l)*ctanh)*hlp3
                dt2d3    = dt2d1
                hlp      = 0.5d0*(1.d0+dtanh(al20*(r(i,2)-r20)))
                s2       = s2a+(t2-s2a)*hlp
                hlp1     = 0.5d0/(dcosh(al20*(r(i,2)-r20))**2)
                ds2d1    = (t2-s2a)*hlp1*(dal20*(r(i,2)-r20)-al20*dr20)
                ds2d2    = ds2ad2+(dt2d2-ds2ad2)*hlp + al20*hlp1*(t2-s2a)
                ds2d3    = ds2d1
                hlp      = aly(l)*(r(i,3)-ry(l))
                cy       = y(2,3,l)+0.5d0*dy(l)*(1.d0+dtanh(hlp))
                dcyd3    = 0.5d0*dy(l)*aly(l)/(dcosh(hlp)**2)
                hlp1     = y(1,3,l)*dexp(-cy*r(i,3))
                hlp2     = y(3,3,l)*dexp(-y(4,3,l)*r(i,3))
                t3       = hlp1+hlp2
                dt3d3    = -(cy+dcyd3*r(i,3))*hlp1-y(4,3,l)*hlp2
                coul1    = 0.5d0*(s1+t1)
                dcl1d1   = 0.5d0*(ds1d1+dt1d1)
                coul2    = 0.5d0*(s2+t2)
                dcl2d1   = 0.5d0*(ds2d1+dt2d1)
                dcl2d2   = 0.5d0*(ds2d2+dt2d2)
                dcl2d3   = 0.5d0*(ds2d3+dt2d3)
                coul3    = 0.5d0*(s3+t3)
                dcl3d3   = 0.5d0*(ds3d3+dt3d3)
                exch1    = 0.5d0*(s1-t1)
                dex1d1   = 0.5d0*(ds1d1-dt1d1)
                exch2    = 0.5d0*(s2-t2)
                dex2d1   = 0.5d0*(ds2d1-dt2d1)
                dex2d2   = 0.5d0*(ds2d2-dt2d2)
                dex2d3   = 0.5d0*(ds2d3-dt2d3)
                exch3    = 0.5d0*(s3-t3)
                dex3d3   = 0.5d0*(ds3d3-dt3d3)
                w        = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
                dwd1     = 2.d0*(exch1-exch2)*(dex1d1-dex2d1) + 2.d0*(exch2-exch3)*dex2d1 - 2.d0*(exch3-exch1)* dex1d1
                dwd2     = - 2.d0*(exch1-exch2)*dex2d2 + 2.d0*(exch2-exch3)*dex2d2
                dwd3     = - 2.d0*(exch1-exch2)*dex2d3 + 2.d0*(exch2-exch3)*(dex2d3-dex3d3) + 2.d0*(exch3-exch1)* dex3d3
                cplg2    = c2a(l)*dexp(-c2b(l)*w-c2c(l)*(r(i,1)+r(i,2)+r(i,3)))
                dcd1     = (-c2b(l)*dwd1 - c2c(l))*cplg2
                dcd2     = (-c2b(l)*dwd2 - c2c(l))*cplg2
                dcd3     = (-c2b(l)*dwd3 - c2c(l))*cplg2
                rootw    = rac2*dsqrt(w+cplg2*cplg2)
                surf     = (coul1+coul2+coul3+dehf-rootw)
                dsurfd1  = (dcl1d1+dcl2d1-0.25d0*(dwd1+2.d0*cplg2*dcd1)/rootw)
                dsurfd2  = (dcl2d2-0.25d0*(dwd2+2.d0*cplg2*dcd2)/rootw)
                dsurfd3  = (dcl2d3+dcl3d3-0.25d0*(dwd3+2.d0*cplg2*dcd3)/rootw)
                hlp1     = acut2(l)*(r(i,2)-rcut2(l))
                hlp2     = adip2(l)*(rnag2-rdip2(l)**2)
                hlp11    = dtanh(hlp1)
                hlp22    = dtanh(hlp2)
                hlp1     = 1.d0/dcosh(hlp1)
                hlp2     = 1.d0/dcosh(hlp2)
                surf     = surf - dip2(l)*hlp1*hlp2
                dsurfd1  = dsurfd1 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*drnag2d1
                dsurfd2  = dsurfd2 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*drnag2d2 + dip2(l)*hlp2*hlp1*hlp11*acut2(l)
                dsurfd3  = dsurfd3 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*drnag2d3
                hlp1     = acut2n(l)*(r(i,2)-rcut2n(l))
                hlp2     = adip2n(l)*(rnag2-rdip2n(l)**2)
                hlp11    = dtanh(hlp1)
                hlp22    = dtanh(hlp2)
                hlp1     = 1.d0/dcosh(hlp1)
                hlp2     = 1.d0/dcosh(hlp2)
                surf     = surf - dip2n(l)*hlp1*hlp2
                dsurfd1  = dsurfd1 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*drnag2d1
                dsurfd2  = dsurfd2 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*drnag2d2 + dip2n(l)*hlp2*hlp1*hlp11* acut2n(l)
                dsurfd3  = dsurfd3 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*drnag2d3
                hlp      = 1.d0/( (1.d0-(csb-bc(l))/(bc(l-1)-bc(l))) * (1.d0-(csb-bc(l))/(bc(l+1)-bc(l))) )
                hlp1     = (csb-bc(l))**2*hlp
                cog      = dexp(-ac(l)*hlp1)
                dcog     = -cog*ac(l)*(2.d0*(csb-bc(l))*hlp +hlp1*(1.d0/(bc(l-1)-csb)+1.d0/(bc(l+1)-csb)))
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
            hlp     = 1.d0/anorm
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
end module chempotpy_nafh_leps_jhctp_2001_dpem_pot22_legacy_state_0db22763

!> @brief Compatibility interface for pot22.
subroutine pot22( r, e, de, nt )
    use chempotpy_nafh_leps_jhctp_2001_dpem_pot22_legacy_state_0db22763, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    double precision :: de(3,nt)
    integer          :: nt
    call dispatch(r=r, e=e, de=de, nt=nt)
end subroutine pot22

!> @brief Evaluate lepscfu11 for the NaFH_LEPS_JHCTP_2001_DPEM potential.
subroutine lepscfu11( g_p_, r, g_r, ldg_r, cfu11, g_cfu11, ldg_cfu11 )
    !     correction function is the difference of two LEPS surfaces
    !     each triplet is weighted by a-f, which was optimezed with a GA
    implicit none
    double precision :: a, b, c, d, e, f
    double precision :: dehf, j1_ap, rehf, shf_ac, yhf, j1_ac, snah, chi, tnah, q1_ap, q1_ac, j2_ap, j2_ac, w, leps_ac, &
        & leps_ap, q2_ap, c7, c6, shf_ap, q2_ac, thf, r2, r1, r3, snaf, ronaf, j3_ac, q3_ap, j3_ap, denah, renah, q3_ac, t2, &
        & t1, t3, tnaf, t4, y, thf_ap, thf_ac, tnaf_ap, tnaf_ac, tnah_ap, tnah_ac, cfu11, de, re, c5, beta, c1, c2, c3, c4, &
        & r(1,3)
    integer            :: ia, ib, ic
    integer, parameter :: g_pmax_ = 3
    integer            :: g_i_, g_p_, ldg_r, ldg_cfu11
    double precision   :: g_r(ldg_r,1,3)
    double precision   :: d4_p, d3_p, d17_b, d5_b, d6_b, d7_b, d8_b, d12_v, d2_v, d11_b
    double precision   :: d10_b, d1_w, d2_w, d3_v, d4_v, d5_v, d1_p, d7_v, d8_v, d2_p
    double precision   :: d10_v, d14_v, d4_b, g_r1(g_pmax_), g_r2(g_pmax_ &
        & ), g_r3(g_pmax_), g_y(g_pmax_), g_d1_w(g_pmax_), g_d2_w(g_pmax_), &
        & g_snaf(g_pmax_)
    double precision :: g_tnaf_ac(g_pmax_), g_tnaf_ap(g_pmax_), g_q3_ac &
        & (g_pmax_), g_q3_ap(g_pmax_), g_j3_ac(g_pmax_), g_j3_ap(g_pmax_), g_beta(g_pmax_), g_chi(g_pmax_), g_snah(g_pmax_), g_tnah_ac(g_pmax_ &
        & )
    double precision :: g_tnah_ap(g_pmax_), g_q1_ac(g_pmax_), g_q1_ap(g_pmax_), g_j1_ac(g_pmax_), g_j1_ap(g_pmax_), g_yhf(g_pmax_), g_shf_ac(g_pmax_), g_shf_ap(g_pmax_), g_thf_ac(g_pmax_), g_thf_ap(g_pmax_)
    double precision :: g_q2_ac(g_pmax_), g_q2_ap(g_pmax_), g_j2_ac(g_pmax_), g_j2_ap(g_pmax_), g_w(g_pmax_), g_leps_ac(g_pmax_), g_leps_ap(g_pmax_), g_cfu11(ldg_cfu11)
    integer          :: g_ehfid
    do ia = 1, 3
        do ib = 1, 1
            do ic = 1, 3
                if( ia .eq. ic ) g_r(ia,ib,ic) = 1.d0
                if( ia .ne. ic ) g_r(ia,ib,ic) = 0.d0
            enddo
        enddo
    enddo
    do g_i_ = 1, g_p_
        g_r1(g_i_) = dble(g_r(g_i_,1,1))
    enddo
    r1 = dble(r(1,1))
    do g_i_ = 1, g_p_
        g_r2(g_i_) = dble(g_r(g_i_,1,2))
    enddo
    r2 = dble(r(1,2))
    do g_i_ = 1, g_p_
        g_r3(g_i_) = dble(g_r(g_i_,1,3))
    enddo
    r3 = dble(r(1,3))
    a  = 3.58749d0
    b  = 2.91300d0
    c  = 2.06256d0
    d  = 3.58749d0
    e  = 3.11828d0
    f  = 2.05279d0
    !     NaF stuff
    !     singlets:  approx and acc the same (from paper)
    c1    = 0.620115d0
    c2    = 0.694459d-1
    c3    = 0.757570d0
    c4    = 0.283055d0
    c5    = 2.69998d0
    ronaf = 3.6395d0
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
    t1 = 7.3251d0
    t2 = 2.8774d0
    t3 = 1930.0d0
    t4 = 2.7091d0
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
        g_q3_ac(g_i_) = 0.5d0 * g_tnaf_ac(g_i_) + 0.5d0 * g_snaf(g_i_)
    enddo
    q3_ac = 0.5d0 * (snaf + tnaf_ac)
    do g_i_ = 1, g_p_
        g_q3_ap(g_i_) = 0.5d0 * g_tnaf_ap(g_i_) + 0.5d0 * g_snaf(g_i_)
    enddo
    q3_ap = 0.5d0 * (snaf + tnaf_ap)
    do g_i_ = 1, g_p_
        g_j3_ac(g_i_) = (-0.5d0) * g_tnaf_ac(g_i_) + 0.5d0 * g_snaf(g_i_)
    enddo
    j3_ac = 0.5d0 * (snaf - tnaf_ac)
    do g_i_ = 1, g_p_
        g_j3_ap(g_i_) = (-0.5d0) * g_tnaf_ap(g_i_) + 0.5d0 * g_snaf(g_i_)
    enddo
    j3_ap = 0.5d0 * (snaf - tnaf_ap)
    !     NaH asymptote
    !     singlets:  approx and acc the same (from paper)
    renah = 7.37707d0
    denah = 0.00291087d0
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0d0
    enddo
    beta = 0.713223d0
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
    d3_v = chi - 2.d0
    d4_b = d2_v + d3_v * denah
    do g_i_ = 1, g_p_
        g_snah(g_i_) = d4_b * g_chi(g_i_)
    enddo
    snah = d2_v * d3_v
    !     triplets:  same for approx and acc
    !     used the 180 degree triplet
    t1 = 6.6796d0
    t2 = 2.7839d0
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
        g_q1_ac(g_i_) = 0.5d0 * g_tnah_ac(g_i_) + 0.5d0 * g_snah(g_i_)
    enddo
    q1_ac = 0.5d0 * (snah + tnah_ac)
    do g_i_ = 1, g_p_
        g_q1_ap(g_i_) = 0.5d0 * g_tnah_ap(g_i_) + 0.5d0 * g_snah(g_i_)
    enddo
    q1_ap = 0.5d0 * (snah + tnah_ap)
    do g_i_ = 1, g_p_
        g_j1_ac(g_i_) = (-0.5d0) * g_tnah_ac(g_i_) + 0.5d0 * g_snah(g_i_)
    enddo
    j1_ac = 0.5d0 * (snah - tnah_ac)
    do g_i_ = 1, g_p_
        g_j1_ap(g_i_) = (-0.5d0) * g_tnah_ap(g_i_) + 0.5d0 * g_snah(g_i_)
    enddo
    j1_ap = 0.5d0 * (snah - tnah_ap)
    !     HF asymptote
    !     accurate singlet
    dehf = 6.122d0
    rehf = 1.733d0
    c1   = 1.1622d0
    c2   = -0.025647d0
    c3   = 0.059062d0
    c4   = 2.1042d0
    do g_i_ = 1, g_p_
        g_yhf(g_i_) = g_r2(g_i_)
    enddo
    yhf  = r2 - c4
    d4_v = yhf * yhf
    d1_p = 2.0d0 * yhf
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
    d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
    d1_p = 2.0d0 * (1.d0 - d2_v)
    d6_b = (-(dehf * d1_p)) * d2_p
    do g_i_ = 1, g_p_
        g_shf_ac(g_i_) = d6_b * g_d1_w(g_i_)
    enddo
    shf_ac = dehf * d4_v - dehf
    !     approx singlet (from paper)
    dehf = 5.77096d0
    c1   = 0.441258d0
    c2   = 3.88056d0
    c3   = -1.441258d0
    c4   = -1.2711d0
    c5   = -1.37766d0
    c6   = 0.168186d0
    c7   = 2.07230d0
    rehf = 1.7328d0
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
    d3_p  = 2.0d0 * yhf
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
    dehf = 5.77096d0
    rehf = 1.7328d0
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0d0
    enddo
    beta = 1.26686217d0
    d4_v = r2 - rehf
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
    enddo
    d1_w = (-beta) * d4_v
    d2_v = exp(d1_w)
    d2_p = d2_v
    d4_v = (1.d0 + d2_v) * (1.d0 + d2_v)
    d1_p = 2.0d0 * (1.d0 + d2_v)
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
    d4_v = (1.d0 + d2_v) * (1.d0 + d2_v)
    d1_p = 2.0d0 * (1.d0 + d2_v)
    d7_b = e * dehf * d1_p * d2_p
    do g_i_ = 1, g_p_
        g_thf_ap(g_i_) = d7_b * g_d1_w(g_i_)
    enddo
    thf_ap = (dehf * d4_v - dehf) * e
    do g_i_ = 1, g_p_
        g_q2_ac(g_i_) = 0.5d0 * g_thf_ac(g_i_) + 0.5d0 * g_shf_ac(g_i_ )
    enddo
    q2_ac = 0.5d0 * (shf_ac + thf_ac)
    do g_i_ = 1, g_p_
        g_q2_ap(g_i_) = 0.5d0 * g_thf_ap(g_i_) + 0.5d0 * g_shf_ap(g_i_ )
    enddo
    q2_ap = 0.5d0 * (shf_ap + thf_ap)
    do g_i_ = 1, g_p_
        g_j2_ac(g_i_) = (-0.5d0) * g_thf_ac(g_i_) + 0.5d0 * g_shf_ac(g_i_)
    enddo
    j2_ac = 0.5d0 * (shf_ac - thf_ac)
    do g_i_ = 1, g_p_
        g_j2_ap(g_i_) = (-0.5d0) * g_thf_ap(g_i_) + 0.5d0 * g_shf_ap(g_i_)
    enddo
    j2_ap = 0.5d0 * (shf_ap - thf_ap)
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
    d1_p = 1.0d0 / (2.0d0 *  d7_v)
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
    d1_p = 1.0d0 / (2.0d0 *  d7_v)
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
    cfu11 = leps_ac - leps_ap + 6.122d0 - 5.77096d0
end subroutine lepscfu11

!> @brief Evaluate lepscfu22 for the NaFH_LEPS_JHCTP_2001_DPEM potential.
subroutine lepscfu22( g_p_, r, g_r, ldg_r, cfu22, g_cfu22, ldg_cfu22 )
    implicit none
    double precision :: dehf, j1_ap, rehf, shf_ac, yhf, j1_ac, snah, chi, tnah, q1_ap, q1_ac, j2_ap, j2_ac, w, leps_ac, &
        & leps_ap, q2_ap, c7, c6, shf_ap, q2_ac, thf, r2, r1, r3, snaf, ronaf, a, b, c, d, e, f, j3_ac, q3_ap, j3_ap, denah, &
        & renah, q3_ac, t2, t1, t3, tnaf, t4, y, thf_ap, thf_ac, tnaf_ap, tnaf_ac, tnah_ap, tnah_ac, r(1,3)
    double precision :: shf_u11, nap, t3nah, t4nah, r0, swit, shfa, alph0, ynaf, snaf_ap, denaf, renaf, t1nah, t2nah, c2naf, snaf_ac
    double precision :: c1nf, c2nf, c3nf, c4nf, c5nf, c6nf, c7nf, hlp1, &
        & hlp2, hlp22, thf_ap_as, t1hf, dt1hf, t2hf, dt2hf, t, b1, b2, e1, &
        & e2, hlp, a3, b3, cc3, d3, e3, f3, cfu22, c5, c4, c3, c2, c1, cc1, &
        & cc2, f1, f2, beta
    integer            :: ia, ib, ic
    integer, parameter :: g_pmax_ = 3
    integer            :: g_i_, g_p_, ldg_r, ldg_cfu22
    double precision   :: d4_p, d3_p, d17_b, d8_b, d12_v, d2_p, d2_w, d10_v, d2_v, d11_b
    double precision   :: d2_b, d1_w, d3_v, d3_b, d1_p, d4_v, d5_v, d6_v, &
        & d7_v, d8_v
    double precision :: d14_v, d4_b, d5_b, d6_b, d7_b, g_r1(g_pmax_), g_r(ldg_r,1,3), g_r2(g_pmax_), g_r3(g_pmax_), g_ynaf(g_pmax_)
    double precision :: g_d1_w(g_pmax_), g_hlp1(g_pmax_), g_hlp2(g_pmax_), g_hlp22(g_pmax_), g_snaf_ap(g_pmax_), g_beta(g_pmax_), g_snaf_ac(g_pmax_), g_hlp(g_pmax_), g_tnaf_ac(g_pmax_), g_tnaf_ap(g_pmax_ &
        & )
    double precision :: g_q3_ac(g_pmax_), g_q3_ap(g_pmax_), g_j3_ac(g_pmax_), g_j3_ap(g_pmax_), g_d2_w(g_pmax_), g_tnah_ac(g_pmax_), g_tnah_ap(g_pmax_), g_q1_ac(g_pmax_), g_q1_ap(g_pmax_), g_j1_ac(g_pmax_)
    double precision :: g_j1_ap(g_pmax_), g_thf_ac(g_pmax_), g_thf_ap(g_pmax_), g_thf(g_pmax_), g_yhf(g_pmax_), g_shf_u11(g_pmax_), g_shfa(g_pmax_), g_swit(g_pmax_), g_shf_ac(g_pmax_), g_thf_ap_as(g_pmax_)
    double precision :: g_shf_ap(g_pmax_), g_q2_ac(g_pmax_), g_q2_ap(g_pmax_), g_j2_ac(g_pmax_), g_j2_ap(g_pmax_), g_w(g_pmax_), g_leps_ac(g_pmax_), g_leps_ap(g_pmax_), g_cfu22(ldg_cfu22)
    integer          :: g_ehfid
    do ia = 1, 3
        do ib = 1, 1
            do ic = 1, 3
                if( ia .eq. ic ) g_r(ia,ib,ic) = 1.d0
                if( ia .ne. ic ) g_r(ia,ib,ic) = 0.d0
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
    b1  = 5.125d0
    b2  = 5.125d0
    cc1 = 1.92d0
    cc2 = 9.01d0
    e1  = 5.54d0
    e2  = 4.15d0
    f1  = 1.62d0
    f2  = 8.80d0
    !     NaF stuff
    !     singlets:  approx (from paper)
    c1nf  = 0.168364d0
    c2nf  = 2.29333d0
    c3nf  = -1.168364d0
    c4nf  = -0.084862d0
    c5nf  = -0.0109707d0
    c6nf  = 0.00246582d0
    c7nf  = 0.377012d0
    denaf = 4.3830525d0
    renaf = 3.6359d0
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
    d2_p = 2.0d0 * ynaf
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
    c1    = 0.32453d0
    c2    = 1.5102d0
    c3    = 3.0938d0
    c4    = 1.7107d0
    denaf = 4.94d0
    renaf = 3.6395d0
    if( (r3 / c3) .ne. 0.0d0 )then
        d3_v = (r3 / c3) ** ( c4 - 2.0d0)
        d3_v = d3_v * (r3 / c3)
        d2_p = c4 *  d3_v
        d3_v = d3_v * (r3 / c3)
    else
        !          ((r3 / c3) = 0)
        d3_v = (r3 / c3) **  c4
        d2_p = 0.0d0
    endif
    if( (r3 / c3) .ne. 0.0d0 )then
        d6_v = (r3 / c3) ** ( c4 - 2.0d0)
        d6_v = d6_v * (r3 / c3)
        d1_p = c4 *  d6_v
        d6_v = d6_v * (r3 / c3)
    else
        !          ((r3 / c3) = 0)
        d6_v = (r3 / c3) **  c4
        d1_p = 0.0d0
    endif
    d7_v = c1 + d6_v
    d8_v = (c2 + d3_v) / d7_v
    d7_b = c1 * ((-d8_v) / d7_v) * d1_p * (1.0d0 / c3) + c1 * (1.0d0 / d7_v) * d2_p * (1.0d0 / c3)
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
    d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
    d1_p = 2.0d0 * (1.d0 - d2_v)
    d6_b = (-(denaf * d1_p)) * d2_p
    do g_i_ = 1, g_p_
        g_snaf_ac(g_i_) = d6_b * g_d1_w(g_i_)
    enddo
    snaf_ac = denaf * d4_v - denaf
    !     triplets:  approx and acc the same
    !     anti-Morse from a constant beta fit to the approx curve
    denaf = 4.49d0
    renaf = 3.6395d0
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0d0
    enddo
    beta = 0.69614d0
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
    d1_p = 2.0d0 * hlp
    d2_b = 0.5d0 * denaf
    d5_b = d2_b * (cc2 * 2.d0) + d2_b * cc1 * d1_p
    do g_i_ = 1, g_p_
        g_tnaf_ac(g_i_) = d5_b * g_hlp(g_i_)
    enddo
    tnaf_ac = 0.5d0 * denaf * (cc1 * d2_v + cc2 * 2.d0 * hlp)
    d2_v    = hlp * hlp
    d1_p    = 2.0d0 * hlp
    d2_b    = 0.5d0 * denaf
    d5_b    = d2_b * (f2 * 2.d0) + d2_b * f1 * d1_p
    do g_i_ = 1, g_p_
        g_tnaf_ap(g_i_) = d5_b * g_hlp(g_i_)
    enddo
    tnaf_ap = 0.5d0 * denaf * (f1 * d2_v + f2 * 2.d0 * hlp)
    do g_i_ = 1, g_p_
        g_q3_ac(g_i_) = 0.5d0 * g_tnaf_ac(g_i_) + 0.5d0 * g_snaf_ac(g_i_)
    enddo
    q3_ac = 0.5d0 * (snaf_ac + tnaf_ac)
    do g_i_ = 1, g_p_
        g_q3_ap(g_i_) = 0.5d0 * g_tnaf_ap(g_i_) + 0.5d0 * g_snaf_ap(g_i_)
    enddo
    q3_ap = 0.5d0 * (snaf_ap + tnaf_ap)
    do g_i_ = 1, g_p_
        g_j3_ac(g_i_) = (-0.5d0) * g_tnaf_ac(g_i_) + 0.5d0 * g_snaf_ac (g_i_)
    enddo
    j3_ac = 0.5d0 * (snaf_ac - tnaf_ac)
    do g_i_ = 1, g_p_
        g_j3_ap(g_i_) = (-0.5d0) * g_tnaf_ap(g_i_) + 0.5d0 * g_snaf_ap (g_i_)
    enddo
    j3_ap = 0.5d0 * (snaf_ap - tnaf_ap)
    !     NaH asymptote
    !     singlets:  approx and acc the same (from paper)
    snah = 0.d0
    !     triplets:  same for approx and acc
    !     used the 180 degree triplet
    t1nah = 1040.7d0
    t2nah = 4.5911d0
    t3nah = 1.4941d0
    t4nah = 4.5525d0
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
        g_q1_ac(g_i_) = 0.5d0 * g_tnah_ac(g_i_)
    enddo
    q1_ac = 0.5d0 * (snah + tnah_ac)
    do g_i_ = 1, g_p_
        g_q1_ap(g_i_) = 0.5d0 * g_tnah_ap(g_i_)
    enddo
    q1_ap = 0.5d0 * (snah + tnah_ap)
    do g_i_ = 1, g_p_
        g_j1_ac(g_i_) = (-0.5d0) * g_tnah_ac(g_i_)
    enddo
    j1_ac = 0.5d0 * (snah - tnah_ac)
    do g_i_ = 1, g_p_
        g_j1_ap(g_i_) = (-0.5d0) * g_tnah_ap(g_i_)
    enddo
    j1_ap = 0.5d0 * (snah - tnah_ap)
    !     HF asymptote
    !     triplets:
    !     anti-morse from a
    !     constant beta fit to MRDCI data
    !     same for approx and accurate
    dehf = 5.77096d0
    rehf = 1.7328d0
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0d0
    enddo
    beta = 1.2669d0
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
    d1_p = 2.0d0 * hlp
    d2_b = 0.5d0 * dehf
    d6_b = d2_b * e2 * d1_p + d2_b * (e1 * 2.d0)
    do g_i_ = 1, g_p_
        g_thf_ac(g_i_) = d6_b * g_hlp(g_i_)
    enddo
    thf_ac = 0.5d0 * dehf * (e1 * 2.d0 * hlp + e2 * d3_v)
    dehf   = 5.77096d0
    rehf   = 1.7328d0
    do g_i_ = 1, g_p_
        g_beta(g_i_) = 0.0d0
    enddo
    beta = 1.2669d0
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
    d1_p = 2.0d0 * hlp
    d2_b = 0.5d0 * dehf
    d6_b = d2_b * b2 * d1_p + d2_b * (b1 * 2.d0)
    do g_i_ = 1, g_p_
        g_thf_ap(g_i_) = d6_b * g_hlp(g_i_)
    enddo
    thf_ap = 0.5d0 * dehf * (b1 * 2.d0 * hlp + b2 * d3_v)
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-15.d0) * g_r2(g_i_)
    enddo
    d1_w = (-15.d0) * (r2 - 2.d0)
    d2_v = exp(d1_w)
    d1_p = d2_v
    do g_i_ = 1, g_p_
        g_thf(g_i_) = d1_p * g_d1_w(g_i_)
    enddo
    thf = d2_v
    !     accurate singlet
    !     uses u11 singlet which is
    dehf = 6.122d0
    rehf = 1.733d0
    c1   = 1.1622d0
    c2   = -0.025647d0
    c3   = 0.059062d0
    c4   = 2.1042d0
    do g_i_ = 1, g_p_
        g_yhf(g_i_) = g_r2(g_i_)
    enddo
    yhf  = r2 - c4
    d4_v = yhf * yhf
    d1_p = 2.0d0 * yhf
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
    d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
    d1_p = 2.0d0 * (1.d0 - d2_v)
    d6_b = (-(dehf * d1_p)) * d2_p
    do g_i_ = 1, g_p_
        g_shf_u11(g_i_) = d6_b * g_d1_w(g_i_)
    enddo
    shf_u11 = dehf * d4_v - dehf
    nap     = 2.0973375d0
    do g_i_ = 1, g_p_
        g_shfa(g_i_) = g_shf_u11(g_i_)
    enddo
    shfa  = shf_u11 + nap
    alph0 = 15.d0
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = (-alph0) * g_thf(g_i_) + alph0 * g_shfa(g_i_)
    enddo
    d1_w = alph0 * (shfa - thf)
    d2_v = tanh (d1_w)
    d1_p = 1.0d0 - ( d2_v *  d2_v)
    do g_i_ = 1, g_p_
        g_swit(g_i_) = d1_p * g_d1_w(g_i_)
    enddo
    swit = d2_v
    d6_v = 0.5d0 * (shfa - thf)
    d6_b = (-swit) * 0.5d0
    d7_b = d6_b + 0.5d0
    d8_b = -d6_b + 0.5d0
    do g_i_ = 1, g_p_
        g_shf_ac(g_i_) = (-d6_v) * g_swit(g_i_) + d8_b * g_thf(g_i_) + d7_b * g_shfa(g_i_)
    enddo
    shf_ac = 0.5d0 * (shfa + thf) - d6_v * swit
    !     approx singlet (from paper)
    !     uses u11 singlet, which is
    dehf = 5.77096d0
    c1   = 0.441258d0
    c2   = 3.88056d0
    c3   = -1.441258d0
    c4   = -1.2711d0
    c5   = -1.37766d0
    c6   = 0.168186d0
    c7   = 2.07230d0
    rehf = 1.7328d0
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
    d3_p  = 2.0d0 * yhf
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
    nap     = 2.0973375d0
    do g_i_ = 1, g_p_
        g_shfa(g_i_) = g_shf_u11(g_i_)
    enddo
    shfa  = shf_u11 + nap
    alph0 = 15.d0
    r0    = 3.1305d0
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = alph0 * g_r2(g_i_)
    enddo
    d1_w = alph0 * (r2 - r0)
    d2_v = tanh (d1_w)
    d1_p = 1.0d0 - ( d2_v *  d2_v)
    d4_b = 0.5d0 * d1_p
    do g_i_ = 1, g_p_
        g_swit(g_i_) = d4_b * g_d1_w(g_i_)
    enddo
    swit  = 0.5d0 * (1.d0 + d2_v)
    t1hf  = 1807.6d0
    dt1hf = 0.d0
    t2hf  = 2.8815d0
    dt2hf = 0.d0
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
    d2_b      = 1.0d0 + (-swit)
    do g_i_ = 1, g_p_
        g_shf_ap(g_i_) = d3_v * g_swit(g_i_) + swit * g_thf_ap_as(g_i_ ) + d2_b * g_shfa(g_i_)
    enddo
    shf_ap = shfa + d3_v * swit
    !     compute coulomb and exchange integrals
    do g_i_ = 1, g_p_
        g_q2_ac(g_i_) = 0.5d0 * g_thf_ac(g_i_) + 0.5d0 * g_shf_ac(g_i_ )
    enddo
    q2_ac = 0.5d0 * (shf_ac + thf_ac)
    do g_i_ = 1, g_p_
        g_q2_ap(g_i_) = 0.5d0 * g_thf_ap(g_i_) + 0.5d0 * g_shf_ap(g_i_ )
    enddo
    q2_ap = 0.5d0 * (shf_ap + thf_ap)
    do g_i_ = 1, g_p_
        g_j2_ac(g_i_) = (-0.5d0) * g_thf_ac(g_i_) + 0.5d0 * g_shf_ac(g_i_)
    enddo
    j2_ac = 0.5d0 * (shf_ac - thf_ac)
    do g_i_ = 1, g_p_
        g_j2_ap(g_i_) = (-0.5d0) * g_thf_ap(g_i_) + 0.5d0 * g_shf_ap(g_i_)
    enddo
    j2_ap = 0.5d0 * (shf_ap - thf_ap)
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
    d1_p = 1.0d0 / (2.0d0 *  d7_v)
    !        else
    !           call ehufDO (9,w, d7_v, d1_p,
    !     +g_ehfid,
    !     +704)
    !        endif
    do g_i_ = 1, g_p_
        g_leps_ac(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ac(g_i_) + g_q2_ac(g_i_) + g_q1_ac(g_i_)
    enddo
    leps_ac = q1_ac + q2_ac + q3_ac - d7_v + 6.122d0
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
    d1_p = 1.0d0 / (2.0d0 *  d7_v)
    !        else
    !           call ehufDO (9,w, d7_v, d1_p,
    !     +g_ehfid,
    !     +746)
    !        endif
    do g_i_ = 1, g_p_
        g_leps_ap(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ap(g_i_) + g_q2_ap(g_i_) + g_q1_ap(g_i_)
    enddo
    leps_ap = q1_ap + q2_ap + q3_ap - d7_v + 5.77096d0
    do g_i_ = 1, g_p_
        g_cfu22(g_i_) = -g_leps_ap(g_i_) + g_leps_ac(g_i_)
    enddo
    cfu22 = leps_ac - leps_ap
end subroutine lepscfu22
