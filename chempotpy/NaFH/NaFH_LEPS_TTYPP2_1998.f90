!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_nafh_leps_ttypp2_1998_pes
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
end module chempotpy_nafh_leps_ttypp2_1998_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_nafh_leps_ttypp2_1998_pes, only: evaluate
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

!> @brief Evaluate EvdRdX for the NaFH_LEPS_TTYPP2_1998 potential.
subroutine evdrdx( x, r, drdx )
    implicit none
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    integer                       :: i, j
    ! Initialize dRdX(3,9)
    do i=1, 3
        do j=1, 9
            drdx(i,j) = 0.0d0
        enddo
    enddo
    drdx(1,1) = (x(1)-x(4))/r(1)
    drdx(1,2) = (x(2)-x(5))/r(1)
    drdx(1,3) = (x(3)-x(6))/r(1)
    drdx(1,4) = -drdx(1,1)
    drdx(1,5) = -drdx(1,2)
    drdx(1,6) = -drdx(1,3)
    drdx(2,4) = (x(4)-x(7))/r(2)
    drdx(2,5) = (x(5)-x(8))/r(2)
    drdx(2,6) = (x(6)-x(9))/r(2)
    drdx(2,7) = -drdx(2,4)
    drdx(2,8) = -drdx(2,5)
    drdx(2,9) = -drdx(2,6)
    drdx(3,1) = (x(1)-x(7))/r(3)
    drdx(3,2) = (x(2)-x(8))/r(3)
    drdx(3,3) = (x(3)-x(9))/r(3)
    drdx(3,7) = -drdx(3,1)
    drdx(3,8) = -drdx(3,2)
    drdx(3,9) = -drdx(3,3)
end subroutine evdrdx

!> @brief Evaluate diagonalize for the NaFH_LEPS_TTYPP2_1998 potential.
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

!
!> @brief Evaluate pot for the NaFH_LEPS_TTYPP2_1998 potential.
subroutine pot( r, e, de, nt, nsurf )
    !
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
    !
    implicit none
    double precision :: de, e, r
    integer          :: nsurf, nt
    dimension r(nt,3), e(nt)
    dimension de(3,nt)
    if( nsurf .eq. 1 )then
        call pot11(r, e, de, nt)
    else if( nsurf .eq. 2 )then
        call pot12(r, e, de, nt)
    else
        call pot22(r, e, de, nt)
    endif
end subroutine pot

!
!
!  U11 surface - lower diabatic surface
!
module chempotpy_nafh_leps_ttypp2_1998_pot11_legacy_state_1329b7c2
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, de, nt )
        implicit none
        optional :: r, e, de, nt
        double precision :: abs1, ac, acut1, adip1, ads1, af1, af2, al2, al2n, aln0, aln0n1, aln0n2, amf, &
            & amh, amna, amu, anorm, ars1, asmall, axs1, bc, c2a, c2b, c2c, cevau, chf, cnf, cnf1, cog, &
            & coul1, coul2, coul3, cplg2, csb, danormd1, danormd2, danormd3, dcd1, dcd2, dcd3, dcl1d1, &
            & dcl2d2, dcl3d3, dcog, dcogd1, dcogd2, dcogd3, dcsbd1, dcsbd2, dcsbd3, de, dehf, denf, &
            & dex1d1, dex2d2, dex3d3, dip1, dipn, dipnn, drnag2d1, drnag2d2, drnag2d3, ds1d1, ds2d2, &
            & ds3d3, dsurfd1, dsurfd2, dsurfd3, dt1d1, dt2d2, dt3d3, dwd1, dwd2, dwd3, e, exch1, exch2, &
            & exch3, hlp, hlp1, hlp11, hlp2, hlp22, hlp23, hlp3, r, r1s, r2s, r3s, rac2, rcut1, rdip1, &
            & rehf, renf, rn0n, rnag2, rootw, s1, s2, s3, surf, t1, t2, t3, una2p, w, y, y2, y3
        integer :: i, l, l1, l2, nt
        dimension r(:,:), e(:)
        dimension de(:,:)
        integer :: g_myrank, g_nprocs
        dimension cnf(6), cnf1(5), chf(6), bc(0:5)
        dimension y(4,3,4), c2a(4), c2b(4), c2c(4), dip1(4), adip1(4), rdip1(4), acut1(4), rcut1(4), dipn(4), al2(4), aln0(4), dipnn(4), al2n(4), aln0n1(4), aln0n2(4), rn0n(4)
        dimension ac(4)
        !-----------------------------------------------------------------------
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
        !
        do 10 i=1, nt
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
            axs1  = dexp(-abs1*(r(i,1)-ars1))
            s1    = ads1*axs1*(axs1-2.d0)
            ds1d1 = 2.d0*ads1*abs1*axs1*(1.d0-axs1)
            !-----------------------------------------------------------------------
            y2    = r(i,2)-rehf
            hlp1  = chf(1)*dexp(-chf(2)*y2)
            hlp2  = dexp(-chf(6)*y2)
            hlp22 = (-1.d0-chf(1)+chf(3)*y2+chf(4)*y2**2+chf(5)*y2**3)*hlp2
            s2    = dehf*(hlp1+hlp22)
            ds2d2 = dehf*(-chf(2)*hlp1 - chf(6)*hlp22 + (chf(3)+2.d0*chf(4)*y2+3.d0*chf(5)*y2**2)*hlp2)
            !-----------------------------------------------------------------------
            y3    = r(i,3)-renf
            hlp1  = dexp(-cnf1(3)*y3)
            hlp11 = (cnf1(1)+cnf1(2)*y3)*hlp1
            hlp2  = cnf1(4)*dexp(-cnf1(5)*y3)
            s3    = hlp11+hlp2
            ds3d3 = -cnf1(3)*hlp11+cnf1(2)*hlp1-cnf1(5)*hlp2
            !-----------------------------------------------------------------------
            e(i)     = 0.d0
            de(1,i)  = 0.d0
            de(2,i)  = 0.d0
            de(3,i)  = 0.d0
            anorm    = 0.d0
            danormd1 = 0.d0
            danormd2 = 0.d0
            danormd3 = 0.d0
            do 20 l = l1, l2
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
                !
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
                dsurfd3 = dsurfd3 + 2.d0*dipnn(l)*hlp1 *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*drnag2d3
                hlp     = 1.d0/( (1.d0-(csb-bc(l))/(bc(l-1)-bc(l))) * (1.d0-(csb-bc(l))/(bc(l+1)-bc(l))) )
                hlp1    = (csb-bc(l))**2*hlp
                cog     = dexp(-ac(l)*hlp1)
                dcog    = -cog*ac(l)*(2.d0*(csb-bc(l))*hlp +hlp1*(1.d0/(bc(l-1)-csb)+1.d0/(bc(l+1)-csb)))
                !-----------------------------------------------------------------------
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
            20 continue
            !     if (csb.gt.3.99d0) write(6,77)l1,l2,csb,r(i,3),
            !    &                         r(i,2),e(i),anorm
            hlp     = 1.d0/anorm
            e(i)    = e(i)*hlp
            de(1,i) = (de(1,i) - e(i)*danormd1)*hlp
            de(2,i) = (de(2,i) - e(i)*danormd2)*hlp
            de(3,i) = (de(3,i) - e(i)*danormd3)*hlp
            !  77 format(2(2x,i3),5(2x,e16.8))
            e(i)    = e(i)*cevau
            de(1,i) = de(1,i)*cevau
            de(2,i) = de(2,i)*cevau
            de(3,i) = de(3,i)*cevau
        10 continue
        !-----------------------------------------------------------------------
    end subroutine dispatch
end module chempotpy_nafh_leps_ttypp2_1998_pot11_legacy_state_1329b7c2

!> @brief Compatibility interface for pot11.
subroutine pot11( r, e, de, nt )
    use chempotpy_nafh_leps_ttypp2_1998_pot11_legacy_state_1329b7c2, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    double precision :: de(3,nt)
    integer          :: nt
    call dispatch(r=r, e=e, de=de, nt=nt)
end subroutine pot11

!
!  U12 surface - coupling surface
!
module chempotpy_nafh_leps_ttypp2_1998_pot12_legacy_state_1329b7c2
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, de, nt )
        implicit none
        optional :: r, e, de, nt
        double precision :: abs1, ads1, af1, af2, ahf, ahft, amf, amh, amna, amu, anaf, anaft, anah, &
            & anaht, ars1, asmall, bc, ccnf, ccnh, ccnha, ccnhi, ccnhi1, cevau, chf, cnf, cnf1, cstsq, &
            & dahftd1, dahftd2, dahftd3, danaftd1, danaftd2, danaftd3, danahtd1, danahtd2, danahtd3, &
            & dccnhd2, dccnhd3, dcstsqd1, dcstsqd2, dcstsqd3, de, dehf, denf, dr13td1, dr13td2, dr13td3, &
            & dr21td1, dr21td2, dr21td3, dr23td1, dr23td2, dr23td3, drnag2d1, drnag2d2, drnag2d3, dsw1d1, &
            & dsw1d2, dsw1d3, dsw21d1, dsw21d2, dsw21d3, dsw23d1, dsw23d2, dsw23d3, dunafd3, dunahd1, &
            & dunahd2, dunahd3, e, hlp, hlp1, hlp2, r, r13, r13t, r1s, r21, r21t, r23, r23t, r2s, r3s, &
            & rehf, renf, rnag2, sw1, sw21, sw23, una2p, unaf, unah, xpnh, xpnh1
        integer :: i, nt
        dimension r(:,:), e(:)
        dimension de(:,:)
        integer :: g_myrank, g_nprocs
        dimension cnf(6), cnf1(5), chf(6), bc(0:5)
        dimension ahf(3), anaf(3), anah(3), r21(3), r23(3), r13(3)
        dimension ccnf(4), ccnha(4), ccnhi(4), ccnhi1(4)
        dimension ccnh(4)
        dimension dccnhd2(4), dccnhd3(4)
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
        !
        do 10 i=1, nt
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
            !-----------------------------------------------------------------------
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
            e(i)    = e(i)*cevau
            de(1,i) = de(1,i)*cevau
            de(2,i) = de(2,i)*cevau
            de(3,i) = de(3,i)*cevau
        10 continue
    end subroutine dispatch
end module chempotpy_nafh_leps_ttypp2_1998_pot12_legacy_state_1329b7c2

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, de, nt )
    use chempotpy_nafh_leps_ttypp2_1998_pot12_legacy_state_1329b7c2, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    double precision :: de(3,nt)
    integer          :: nt
    call dispatch(r=r, e=e, de=de, nt=nt)
end subroutine pot12

!
!  U22 surface - upper diabatic surface
!
module chempotpy_nafh_leps_ttypp2_1998_pot22_legacy_state_1329b7c2
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, de, nt )
        implicit none
        optional :: r, e, de, nt
        double precision :: abs1, ac, acut2, acut2n, adip2, adip2n, ads1, af1, af2, al20, al20a, al20i, &
            & al20i1, al2a, aly, amf, amh, amna, amu, anorm, ars1, asmall, bc, c2a, c2b, c2c, cevau, chf, &
            & cnf, cnf1, cog, coul1, coul2, coul3, cplg2, csb, ctanh, cy, dal20, danormd1, danormd2, &
            & danormd3, dcd1, dcd2, dcd3, dcl1d1, dcl2d1, dcl2d2, dcl2d3, dcl3d3, dcog, dcogd1, dcogd2, &
            & dcogd3, dcsbd1, dcsbd2, dcsbd3, dctanh, dcyd3, de, dehf, denf, dex1d1, dex2d1, dex2d2, &
            & dex2d3, dex3d3, dip2, dip2n, dr20, drnag2d1, drnag2d2, drnag2d3, ds1d1, ds2ad2, ds2d1, &
            & ds2d2, ds2d3, ds3d3, dsurfd1, dsurfd2, dsurfd3, dt1d1, dt2d1, dt2d2, dt2d3, dt3d3, dwd1, &
            & dwd2, dwd3, dy, dy1, dy2, e, exch1, exch2, exch3, hlp, hlp1, hlp11, hlp2, hlp22, hlp3, psi, &
            & r, r1s, r20, r20a, r20i, r2s, r3s, rac2, rcut2, rcut2n, rdip2, rdip2n, rehf, renf, rnag2, &
            & rootw, rr0a, rsi, ry, s1, s2, s2a, s3, surf, t1, t2, t3, una2p, w, y, y2, y3
        integer :: i, l, l1, l2, nt
        dimension r(:,:), e(:)
        dimension de(:,:)
        integer :: g_myrank, g_nprocs
        dimension cnf(6), cnf1(5), chf(6), bc(0:5)
        dimension y(4,3,4), c2a(4), c2b(4), c2c(4), r20i(4), al20i(4), al20i1(4), al2a(4), rr0a(4), dip2(4), adip2(4), rdip2(4), acut2(4), rcut2(4), ry(4), aly(4), dy(4), dy1(4), dy2(4), dip2n(4), adip2n(4), rdip2n(4), acut2n(4), rcut2n(4)
        dimension ac(4)
        !-----------------------------------------------------------------------
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
        !
        !
        !
        do 10 i=1, nt
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
            !-----------------------------------------------------------------------
            s1    = 0.d0
            ds1d1 = 0.d0
            !-----------------------------------------------------------------------
            y2     = r(i,2)-rehf
            hlp1   = chf(1)*dexp(-chf(2)*y2)
            hlp2   = dexp(-chf(6)*y2)
            hlp22  = (-1.d0-chf(1)+chf(3)*y2+chf(4)*y2**2+chf(5)*y2**3)*hlp2
            s2a    = dehf*(hlp1+hlp22) + una2p
            ds2ad2 = dehf*(-chf(2)*hlp1-chf(6)*hlp22 + (chf(3)+2.d0*chf(4)*y2+3.d0*chf(5)*y2**2)*hlp2)
            !-----------------------------------------------------------------------
            y3    = r(i,3)-renf
            hlp1  = cnf(1)*dexp(-cnf(2)*y3)
            hlp2  = dexp(-cnf(6)*y3)
            hlp22 = (-1.d0-cnf(1)+cnf(3)*y3+cnf(4)*y3*y3+cnf(5)*y3**3)*hlp2
            s3    = denf*(hlp1+hlp22)
            ds3d3 = denf*(-cnf(2)*hlp1-cnf(6)*hlp22 + (cnf(3)+2.d0*cnf(4)*y3+3.d0*cnf(5)*y3**2)*hlp2)
            !-----------------------------------------------------------------------
            e(i)     = 0.d0
            de(1,i)  = 0.d0
            de(2,i)  = 0.d0
            de(3,i)  = 0.d0
            anorm    = 0.d0
            danormd1 = 0.d0
            danormd2 = 0.d0
            danormd3 = 0.d0
            do 20 l = l1, l2
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
            20 continue
            hlp     = 1.d0/anorm
            e(i)    = e(i)*hlp
            de(1,i) = (de(1,i) - e(i)*danormd1)*hlp
            de(2,i) = (de(2,i) - e(i)*danormd2)*hlp
            de(3,i) = (de(3,i) - e(i)*danormd3)*hlp
            e(i)    = e(i)*cevau
            de(1,i) = de(1,i)*cevau
            de(2,i) = de(2,i)*cevau
            de(3,i) = de(3,i)*cevau
        10 continue
    end subroutine dispatch
end module chempotpy_nafh_leps_ttypp2_1998_pot22_legacy_state_1329b7c2

!> @brief Compatibility interface for pot22.
subroutine pot22( r, e, de, nt )
    use chempotpy_nafh_leps_ttypp2_1998_pot22_legacy_state_1329b7c2, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    double precision :: de(3,nt)
    integer          :: nt
    call dispatch(r=r, e=e, de=de, nt=nt)
end subroutine pot22
