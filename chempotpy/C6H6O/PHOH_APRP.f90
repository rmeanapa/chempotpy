!   System:          C6H6O
!   Functional form: Anchor-Points Reactive Potential (APRP)
!                               (3x3 diabatic fit)
!   Common name:                   PhOH(Coupled PESs)
!   Number of derivatives:         1
!   Number of bodies:              13
!   Number of electronic surfaces: 3
!   Interface:                     Section-2
!   Reference: Ke R. Yang, Xuefei Xu, Jingjing Zheng,
!   and Donald G. Truhlar Chem. Sci. (to be published).
!   Notes:
!    1. This routine calculates three coupled singlet PESs for
!       nonadiabatic photodissociation of PhOH to produce
!       PhO + H.
!    2. The final version of this potential was added to POTLIB
!       on Aug. 3, 2014. A preliminary version in the library before
!       that is not correct.
!   Numbering of atoms:
!        H11      H12
!         \      /
!          C5---C6
!         /      \
!   H10--C4       C1---O7
!         \      /      \
!          C3---C2       H13
!         /      \
!        H9       H8
!   Units:      bohr for length
!               hartree for energy
!   Input:
!   igrad = 0           Energy only calculation
!         = 1           Energy + Analytic gradient
!   repflag             flag to indicate wheter to use diabatic or
!                       adiabatic representation (used by ANT program)
!   xx(1,i)             X coordinate of atom i
!   xx(2,i)             Y coordinate of atom i
!   xx(3,i)             Z coordinate of atom i
!   Output:
!   uu(j,j)             diabatic potential of diabatic state j
!   uu(j,k)             diabatic coupling between state j and k
!   guu(1,i,j,j)        X component of the gradients of diabatical
!                       potential of state j at atom i
!   guu(2,i,j,j)        Y component of the gradients of diabatical
!                       potential of state j at atom i
!   guu(3,i,j,j)        Z component of the gradients of diabatical
!                       potential of state j at atom i
!   guu(1,i,j,k)        X component of the gradients of diabatical
!                       coupling between state j and k at atom i
!   guu(2,i,j,k)        Y component of the gradients of diabatical
!                       coupling between state j and k at atom i
!   guu(3,i,j,k)        Z component of the gradients of diabatical
!                       coupling between state j and k at atom i
!   vv(j)               adiabatic potential of adiabatic state j
!   gvv(1,i,j)          X component of the gradients of adiabatic
!                       potential of state j at atom i
!   gvv(2,i,j)          Y component of the gradients of adiabatic
!                       potential of state j at atom i
!   gvv(3,i,j)          Z component of the gradients of adiabatic
!                       potential of state j at atom i
!   dvec(1,i,j,k)       X component of nonadiabatic coupling between
!                       state j and k at atom i
!   dvec(2,i,j,k)       Y component of nonadiabatic coupling between
!                       state j and k at atom i
!   dvec(3,i,j,k)       Z component of nonadiabatic coupling between
!                       state j and k at atom i
!   cc(3,3)             3*3 orghtonal matrix diagonalizing UU matrix
!                       (CC^T*UU*CC yield diagonal matrix with adiabatic
!                       energies as diagonal elements)
!   Note:       LAPACK library is needed to diagonize diabatic potential
!               matrix (UU) to yield adiabatic potenitials VV(3)

!> @brief Implement PHOH_APRP, a potential energy surface for C6H6O.
module chempotpy_phoh_aprp_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the PHOH_APRP potential energy surface for C6H6O.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units, only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 3
        integer, parameter    :: natoms  = 13
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: xx(3,natoms), uu(nstates,nstates)
        real(wp)              :: guu(3,natoms,nstates,nstates)
        real(wp)              :: vv(nstates), gvv(3,natoms,nstates)
        real(wp)              :: dvec(3,natoms,nstates,nstates)
        real(wp)              :: cc(nstates,nstates)
        integer               :: iatom, idir, j, istate, jstate
        !initialize
        p = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                xx(idir,iatom) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        call pot(igrad, xx, uu, guu, vv, gvv, dvec, cc, 0)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = vv(istate)*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = vv(istate)*hartree_to_electronvolt
            enddo
            do istate=1, nstates
                do iatom=1, natoms
                    do idir=1, 3
                        g(istate,iatom,idir) = gvv(idir,iatom,istate)*hartree_per_bohr_to_ev_per_angstrom
                    enddo
                enddo
            enddo
        elseif( igrad == 2 )then
            do istate=1, nstates
                p(istate) = vv(istate)*hartree_to_electronvolt
            enddo
            do istate=1, nstates
                do iatom=1, natoms
                    do idir=1, 3
                        g(istate,iatom,idir) = gvv(idir,iatom,istate)*hartree_per_bohr_to_ev_per_angstrom
                    enddo
                enddo
            enddo
            do istate=1, nstates
                do jstate=1, nstates
                    do iatom=1, natoms
                        do idir=1, 3
                            d(istate,jstate,iatom,idir) = dvec(idir,iatom,istate,jstate)/0.5
                        enddo
                    enddo
                enddo
            enddo
        endif
    end subroutine evaluate

    !> @brief Evaluate pot for the PHOH_APRP potential.
    subroutine pot( igrad, xx, uu, guu, vv, gvv, dvec, cc, repflag )
        use chempotpy_kinds, only: wp
        implicit none
        integer            :: i, j, k, l, m
        integer            :: igrad, repflag
        integer, parameter :: natm   = 13
        integer, parameter :: nsurft = 3
        real(wp)           :: uu(3,3), vv(3)
        real(wp)           :: guu(3,natm,3,3), uup(3,3), vvp(3), vvm(3)
        real(wp)           :: guup(3,natm,3,3), gvv(3,natm,nsurft)
        real(wp)           :: guum(3,natm,3,3)
        real(wp)           :: xx(3,natm), x(natm), y(natm), z(natm)
        real(wp)           :: dvec(3,natm,nsurft,nsurft)
        !     double precision cc(nsurft,nsurft),ccp(nsurft,nsurft)
        real(wp)            :: cc(3,3), ccp(3,3)
        real(wp)            :: tmp
        real(wp), parameter :: stp = 1.e-4_wp
        ! calculate adiabatic and diabatic energies in hartree
        !     write(6,*) 'input geometry'
        do i = 1, natm
            x(i) = xx(1,i)
            y(i) = xx(2,i)
            z(i) = xx(3,i)
            !        write(6,'(3F10.6)') xx(:,i)*0.5291772192d0
        enddo
        !     vv(:) = 0d0
        !     uu(:,:) = 0d0
        !     guu(:,:,:,:) = 0d0
        !     cc(:,:) = 0d0
        call poten(1, x, y, z, uu, guu, vv, gvv, cc)
        ! Compute the nonadiabatic coupling vectors DVEC.
        call getdvec3(natm, uu, guu, vv, dvec, cc)
        !     write(6,*) 'diabatic energy '
        !     write(6,'(5F10.6,a)') uu(1,1),uu(2,2),uu(3,3),uu(1,3),uu(2,3)
        !     write(6,'(a,5F10.6,a)') 'adiabatic energy ', vv(:)
        !     write(6,'(5F10.6,a)') vv(:)
        !     write(6,*) 'diabatic gradients'
        !     do i = 1,natm
        !        write(6,'(5F12.8)') guu(1,i,1,1),guu(1,i,2,2),guu(1,i,3,3),
        !    $                       guu(1,i,1,3),guu(1,i,2,3)
        !        write(6,'(5F12.8)') guu(2,i,1,1),guu(2,i,2,2),guu(2,i,3,3),
        !    $                       guu(2,i,1,3),guu(2,i,2,3)
        !        write(6,'(5F12.8)') guu(3,i,1,1),guu(3,i,2,2),guu(3,i,3,3),
        !    $                       guu(3,i,1,3),guu(3,i,2,3)
        !     enddo
        !     write(6,*) 'adiabatic gradients'
        !     do i = 1,natm
        !       write(6,'(5F10.6,a)') gvv(1,i,1),gvv(1,i,2),gvv(1,i,3)
        !       write(6,'(5F10.6,a)') gvv(2,i,1),gvv(2,i,2),gvv(2,i,3)
        !       write(6,'(5F10.6,a)') gvv(3,i,1),gvv(3,i,2),gvv(3,i,3)
        !     enddo
        !     write(6,*) 'adiabatic coupling dvec'
        !     do i = 1, natm
        !       write(6,'(5F10.6,a)') dvec(1,i,1,3),dvec(1,i,2,3)
        !       write(6,'(5F10.6,a)') dvec(2,i,1,3),dvec(2,i,2,3)
        !       write(6,'(5F10.6,a)') dvec(3,i,1,3),dvec(3,i,2,3)
        !     enddo
        !     stop
    end subroutine pot

    !> @brief Evaluate getdvec3 for the PHOH_APRP potential.
    subroutine getdvec3( nclu, pemd, gpemd, pema, dvec, cc )
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter :: nsurft = 3
        real(wp)           :: pema(nsurft), pemd(nsurft,nsurft)
        real(wp)           :: gpemd(3,nclu,nsurft,nsurft)
        real(wp)           :: dvec(3,nclu,nsurft,nsurft), cc(nsurft,nsurft)
        integer            :: i1, i2, i, j, k, l, nclu
        !     Program to compute nonadiabatic coupling vectors from the adiabati
        !     the gradients of diabatic potential matrix elements for a multisur
        !     The formula is taken from ANT07 and is described in M.D. Hack et a
        ! compute d
        dvec(:,:,:,:) = 0e0_wp
        do i1 = 1, 3
            do i2 = 1, nclu
                do i = 2, nsurft
                    do j = 1, i-1
                        !            dvec(i1,i2,i,j) = 0.d0
                        do k = 1, nsurft
                            do l = 1, nsurft
                                dvec(i1,i2,i,j) = dvec(i1,i2,i,j) + cc(k,i)*cc(l,j)*gpemd(i1,i2,k,l)
                            enddo
                        enddo
                        if( (pema(j) - pema(i)) .ne. 0.0e0_wp )then
                            dvec(i1,i2,i,j) = dvec(i1,i2,i,j)/(pema(j)-pema(i))
                        else
                            dvec(i1,i2,i,j) = 0.0e0_wp
                        endif
                        !              write(6,*) "dvec ",i1,i2,i,j,dvec(i1,i1,i,j)
                        dvec(i1,i2,j,i) = -dvec(i1,i2,i,j)
                    enddo
                enddo
            enddo
        enddo
        !       stop
    end subroutine getdvec3

    !> @brief Evaluate poten for the PHOH_APRP potential.
    subroutine poten( igrad, x, y, z, uu, guu, vv, gvv, cc )
        ! Little subroutine to call evuii and evuij
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter    :: natm = 13
        integer, parameter    :: nap  = 4
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natm), y(natm), z(natm)
        real(wp), intent(out) :: uu(3,3), vv(3), guu(3,natm,3,3)
        real(wp)              :: u12
        integer               :: i, i1, i2, j, k
        real(wp)              :: rap(nap)
        integer               :: ir
        real(wp)              :: roh
        real(wp)              :: t(nap), dtdr(nap)
        real(wp)              :: u1, u2, u3
        ! Add by KRY
        real(wp) :: ubm
        real(wp) :: dubmdx(3*natm)
        real(wp) :: du1dx(3*natm), du2dx(3*natm), du3dx(3*natm)
        real(wp) :: u13, u23
        !      double precision du13dX(3*natm),du23dX(3*natm)
        real(wp) :: du12dx(3*natm), du13dx(3*natm), du23dx(3*natm)
        ! add by J. Zheng
        integer               :: lwork, info
        real(wp)              :: cc(3,3)
        real(wp), allocatable :: work(:)
        !C add by K. Yang
        real(wp) :: agvv(3*natm,3)
        real(wp) :: gvv(3,natm,3)
        real(wp) :: fad(3*natm,3,3)
        ! Determine rOH and determine ir, t, and dtdr
        rap = reshape([ 1.93206106e0_wp, 2.49443845e0_wp, 4.27078097e0_wp, 9.44863047e0_wp ], shape(rap))
        roh = bndlen(7, 13, x, y, z)
        call tent(nap, rap, roh, ir, t, dtdr)
        do i=1, 3*natm
            du1dx(i) = 0.0e0_wp
            du2dx(i) = 0.0e0_wp
            du3dx(i) = 0.0e0_wp
        enddo
        ! Determine u13 and u13 and their Cartesian gradients
        !      call evuij(nap,natm,igrad,ir,t,dtdr,X,Y,Z,u13,u23,
        !     $           du13dX,du23dX)
        ! Determin u12, u13, and u23 and their Cartesian gradients  !KRY 2014071
        call evuij(nap, natm, igrad, ir, t, dtdr, x, y, z, u12, u13, u23, &
            & du12dx, du13dx, du23dx)
        ! Determine u11, u22, u33, and their Cartesian gradients
        call evuii(igrad, natm, x, y, z, u1, du1dx, 11)
        call evuii(igrad, natm, x, y, z, u2, du2dx, 22)
        call evuii(igrad, natm, x, y, z, u3, du3dx, 33)
        ! KRY Born-Mayer repulsion
        call evubm(igrad, natm, x, y, z, ubm, dubmdx)
        ! add by J. Zheng
        guu(:,:,:,:) = 0e0_wp
        do i = 1, natm
            do j = 1, 3
                !        guu(j,i,1,1) = du1dx(3*i-3+j)
                !        guu(j,i,2,2) = du2dx(3*i-3+j)
                !        guu(j,i,3,3) = du3dx(3*i-3+j)
                guu(j,i,1,1) = du1dx(3*i-3+j)+dubmdx(3*i-3+j)
                guu(j,i,2,2) = du2dx(3*i-3+j)+dubmdx(3*i-3+j)
                guu(j,i,3,3) = du3dx(3*i-3+j)+dubmdx(3*i-3+j)
                ! add by K. R. Yang 20140718
                guu(j,i,1,2) = du12dx(3*i-3+j)
                guu(j,i,2,1) = du12dx(3*i-3+j)   ! add du12dx
                guu(j,i,1,3) = du13dx(3*i-3+j)
                guu(j,i,3,1) = du13dx(3*i-3+j)
                guu(j,i,2,3) = du23dx(3*i-3+j)
                guu(j,i,3,2) = du23dx(3*i-3+j)
            enddo
        enddo
        uu(:,:) = 0e0_wp
        !      uu(1,1) = u1
        !      uu(2,2) = u2
        !      uu(3,3) = u3
        uu(1,1) = u1+ubm
        uu(2,2) = u2+ubm
        uu(3,3) = u3+ubm
        ! add by K. R. Yang 20140718
        uu(1,2) = u12
        uu(2,1) = u12       ! add u12
        uu(1,3) = u13
        uu(3,1) = u13
        uu(2,3) = u23
        uu(3,2) = u23
        ! diagonalization to get adiabatic energies
        ! Diagonalize the coupling matrix (from LAPACK)
        !Calls the LAPACK diagonalization subroutine DSYEV        !
        !input:  u(n,n) = real symmetric matrix to be diagonalized!
        !            n  = size of u                               !
        !output: cc(n,n) = orthonormal eigenvectors of u          !
        !        vv(n) = eigenvalues of u in ascending order      !
        cc(:,:) = uu(:,:)
        lwork   = -1
        allocate (work (1) )
        call dsyev('V', 'L', 3, cc, 3, vv, work, lwork, info)
        lwork = int(work(1))
        deallocate (work)
        allocate( work(lwork) )
        call dsyev('V', 'L', 3, cc, 3, vv, work, lwork, info)
        deallocate (work)
        if( info .ne. 0 )then
            write(6, *) 'Error in diagolization of U matrix'
            stop
        endif
        ! Test
        if( igrad .eq. 1 )then
            !     write(*,*) "C matrix:"
            !     do i=1,3
            !       write(*,'(3F12.8)') (cc(i,j),j=1,3)
            !     enddo
            gvv(:,:,:) = 0e0_wp
            !      do i1=1,natm
            !      do i2=1,3
            do i1=1, 3
                do i2=1, natm
                    do j=1, 3
                        do k=1, 3
                            gvv(i1,i2,1) = gvv(i1,i2,1)+cc(k,1)*cc(j,1)*guu(i1,i2,j,k)
                            gvv(i1,i2,2) = gvv(i1,i2,2)+cc(k,2)*cc(j,2)*guu(i1,i2,j,k)
                            gvv(i1,i2,3) = gvv(i1,i2,3)+cc(k,3)*cc(j,3)*guu(i1,i2,j,k)
                        enddo
                    enddo
                enddo
            enddo
        endif
        ! END
    end subroutine poten

    !> @brief Evaluate evuii for the PHOH_APRP potential.
    subroutine evuii( igrad, natm, x, y, z, u, dudx, ii )
        ! Subroutine to evaluate the diabatic potentials and gradients
        ! The bohr-hartree-radian unit is used throughout the program
        ! qtc:  redundant internal coordinates
        ! bmqt: B-matrix of RICs
        ! ii:   flag to indicate the diabatic states
        !       = 11, U11
        !       = 22, U22
        !       = 33, U33
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: dudph, dudr, dudth, phi, r, theta, ub, ur, ut, uto
        integer             :: i, ii, j
        integer, parameter  :: maxdim = 100
        integer, parameter  :: ldbmt  = maxdim
        integer, parameter  :: ldbms  = maxdim
        integer, parameter  :: ldbnrt = maxdim
        integer, parameter  :: ldct   = maxdim
        integer, parameter  :: nstc   = 9
        real(wp), parameter :: u1min  = 0.00905097e0_wp
        ! Global variables
        integer  :: igrad
        integer  :: natm
        real(wp) :: u
        real(wp) :: x(*), y(*), z(*)
        real(wp) :: dudx(3*natm)
        ! Local variables
        integer  :: nqtc
        integer  :: nblt, nbat, ntot, noopt
        integer  :: iblt(2,12), ibat(3,18), itot(4,24), ioopt(4,0)
        integer  :: nbls, nbas, ntos, noops
        integer  :: ibls(2,1), ibas(3,1), itos(4,1), ioops(4,0)
        real(wp) :: qtc(maxdim), qsc(maxdim)
        real(wp) :: dutdr, dubdr, dutodr, durdr
        real(wp) :: dtodph, dubdth
        real(wp) :: dut(maxdim)
        real(wp) :: bmt(ldbmt,3*natm), bms(ldbms,3*natm)
        ! Evaluate internal coordinates
        nblt  = 12
        nbat  = 18
        ntot  = 24
        noopt = 0
        iblt  = reshape([ 1, 2, 1, 6, 1, 7, 2, 3, 2, 8, 3, 4, 3, 9, 4, 5, 4, 10, 5, 6, 5, 11, 6, 12 ], shape(iblt))
        ibat = reshape([ 2, 1, 6, 2, 1, 7, 6, 1, 7, 1, 2, 3, 1, 2, 8, 3, 2, 8, 2, 3, 4, 2, 3, 9, 4, 3, 9, 3, 4, 5, 3, &
            & 4, 10, 5, 4, 10, 4, 5, 6, 4, 5, 11, 6, 5, 11, 1, 6, 5, 1, 6, 12, 5, 6, 12 ], shape(ibat))
        itot = reshape([ 6, 1, 2, 3, 6, 1, 2, 8, 7, 1, 2, 3, 7, 1, 2, 8, 2, 1, 6, 5, 2, 1, 6, 12, 7, 1, 6, 5, 7, 1, 6, &
            & 12, 1, 2, 3, 4, 1, 2, 3, 9, 8, 2, 3, 4, 8, 2, 3, 9, 2, 3, 4, 5, 2, 3, 4, 10, 9, 3, 4, 5, 9, 3, 4, 10, 3, &
            & 4, 5, 6, 3, 4, 5, 11, 10, 4, 5, 6, 10, 4, 5, 11, 4, 5, 6, 1, 4, 5, 6, 12, 11, 5, 6, 1, 11, 5, 6, 12 ], &
            & shape(itot))
        nbls  = 1
        nbas  = 1
        ntos  = 1
        noops = 0
        ibls  = reshape([ 7, 13 ], shape(ibls))
        ibas  = reshape([ 1, 7, 13 ], shape(ibas))
        itos  = reshape([ 2, 1, 7, 13 ], shape(itos))
        call evintcoord(nblt, iblt, nbat, ibat, ntot, itot, noopt, ioopt, &
            & nqtc, qtc, x, y, z)
        r     = bndlen(7, 13, x, y, z)
        theta = bndang(1, 7, 13, x, y, z)
        phi   = dihedr(2, 1, 7, 13, x, y, z)
        ! Evaluate Wilson's B-matrix if igrad = 1
        if( igrad .eq. 1 )then
            ! Added by J. Zheng
            !        bmt(:,:) = 0d0
            call bmat3(nblt, iblt, nbat, ibat, ntot, itot, noopt, ioopt, natm, &
                & x, y, z, bmt, ldbmt)
            ! Added by J. Zheng
            !        bms(:,:) = 0d0
            call bmat3(nbls, ibls, nbas, ibas, ntos, itos, noops, ioops, natm, &
                & x, y, z, bms, ldbms)
        endif
        ! Evaluate uii from tertiary coordinates
        call evuiit(igrad, r, nqtc, qtc, ut, dut, dutdr, ii)
        ! Evaluate uii from reactive (rOH) and secondary coordinates
        ! (bending C-O-H and torsion C-C-O-H)
        if( ii .eq. 11 )then
            call evu1r(igrad, r, ur, durdr)
            call evu1b(igrad, r, theta, ub, dubdr, dubdth)
            call evu1to(igrad, r, phi, uto, dutodr, dtodph)
        elseif( ii .eq. 22 )then
            call evu2r(igrad, r, ur, durdr)
            call evu2b(igrad, r, theta, ub, dubdr, dubdth)
            call evu2to(igrad, r, phi, uto, dutodr, dtodph)
        elseif( ii .eq. 33 )then
            call evu3r(igrad, r, ur, durdr)
            call evu3b(igrad, r, theta, ub, dubdr, dubdth)
            call evu3to(igrad, r, phi, uto, dutodr, dtodph)
        endif
        ! Readjust zero of energy to have U1 min to be zero
        u     = ur+ub+uto+ut-u1min
        dudr  = durdr+dubdr+dutodr+dutdr
        dudth = dubdth
        dudph = dtodph
        ! Calculated analytical Cartesian gradient if igrad = 1
        if( igrad .eq. 1 )then
            do i=1, 3*natm
                dudx(i) = 0.0e0_wp
            enddo
            ! Contribution of dudX from dudr
            do i=1, 3*natm
                dudx(i) = dudx(i)+dudr*bms(1,i)
            enddo
            ! Contribution of dudX from dudth
            do i=1, 3*natm
                dudx(i) = dudx(i)+dudth*bms(2,i)
            enddo
            ! Contribution of dudX from dudph
            do i=1, 3*natm
                dudx(i) = dudx(i)+dudph*bms(3,i)
            enddo
            ! Contribution of dudX from dut
            do i=1, 3*natm
                do j=1, nqtc
                    dudx(i) = dudx(i)+dut(j)*bmt(j,i)
                enddo
            enddo
        endif
    end subroutine evuii

    !      subroutine evuiit(igrad,r,ntot,nqtc,qtc,ut,dut,dutdr,ii)
    !> @brief Evaluate evuiit for the PHOH_APRP potential.
    subroutine evuiit( igrad, r, nqtc, qtc, ut, dut, dutdr, ii )
        ! Subroutine to evaluate the diabatic potentials and gradients
        ! The bohr-hartree-radian unit is used throughout the program
        ! qtc:  redundant internal coordinates
        ! ii:   flag to indicate the diabatic states
        !       = 11, U11
        !       = 22, U22
        !       = 33, U33
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter :: nap = 4
        !      integer nqa,nqb
        !      parameter (nqa=30,nqb=24)
        !      parameter (nqs=12,nqb=18,nqt=24)
        integer, parameter :: nqs = 12
        integer, parameter :: nqb = 18
        integer, parameter :: nqt = 24
        ! Global variables
        integer  :: igrad, ii
        integer  :: ntot, nqtc
        real(wp) :: r, ut, dutdr
        real(wp) :: qtc(*), dut(*)
        ! Local variables
        integer  :: i, j, k
        integer  :: ir
        real(wp) :: u0, u01, u02, u1, u2
        real(wp) :: vs, vs1, vs2
        real(wp) :: vb, vb1, vb2
        real(wp) :: vt, vt1, vt2
        real(wp) :: vbs, vbs1, vbs2
        real(wp) :: dvsdr, dvbdr, dvtdr, dvbsdr
        real(wp) :: t(nap), dtdr(nap)
        real(wp) :: rap(nap)
        real(wp) :: v0(nap)
        real(wp) :: qs(nqs), qs1(nqs), qs2(nqs)
        real(wp) :: qb(nqb), qb1(nqb), qb2(nqb)
        real(wp) :: sqb(nqb)
        real(wp) :: qt(nqt), qt1(nqt), qt2(nqt)
        real(wp) :: sqt1(nqt), cqt1(nqt)
        real(wp) :: sqt2(nqt), cqt2(nqt)
        real(wp) :: dvs(nqs), dvs1(nqs), dvs2(nqs)
        real(wp) :: dvb(nqb), dvb1(nqb), dvb2(nqb)
        real(wp) :: dvt(nqt), dvt1(nqt), dvt2(nqt)
        real(wp) :: dvbsdb(nqb), dvbsdb1(nqb), dvbsdb2(nqb)
        real(wp) :: dvbsds(nqs), dvbsds1(nqs), dvbsds2(nqs)
        real(wp) :: qse(nqs,nap)
        real(wp) :: qbe(nqb,nap)
        real(wp) :: qte(nqt,nap)
        real(wp) :: qse1(nqs), qse2(nqs)
        real(wp) :: qbe1(nqb), qbe2(nqb)
        real(wp) :: fs1(nqs,nqs), fs2(nqs,nqs)
        real(wp) :: fb1(nqb,nqb), fb2(nqb,nqb)
        real(wp) :: ft1(nqt,nqt), ft2(nqt,nqt)
        real(wp) :: fbs1(nqb,nqs), fbs2(nqb,nqs)
        real(wp) :: hs(nqs,nqs,nap)
        real(wp) :: hb(nqb,nqb,nap)
        real(wp) :: ht(nqt,nqt,nap)
        real(wp) :: hbs(nqb,nqs,nap)
        ! Initialize output variables
        ut = 0.0e0_wp; dutdr = 0.0e0_wp
        do i=1, nqtc
            dut(i) = 0.0e0_wp
        enddo
        ! Get HS for stretchings, HB for bendings, and HT for torsions
        call gethess(nap, hs, nqs, hb, nqb, ht, nqt, hbs, ii)
        ! Get equilibrium parameters
        call getgeom(nap, qse, nqs, qbe, nqb, qte, nqt, ii)
        if( ii .eq. 11 )then
            rap(1) = 1.93206106e0_wp ! rOH=1.0224 using rOH=0.964A
            rap(2) = 2.49443845e0_wp ! rOH=1.32 A
            rap(3) = 3.77945219e0_wp ! rOH=2.00 A
            rap(4) = 9.44863047e0_wp ! rOH=5.00 A
            ! Normalized relaxation energy
            ! dE = E(DFT_opt) - E(CAS_geom)
            v0(1) = 0.0e0_wp
            v0(2) = -0.000346019e0_wp
            v0(3) = 0.000672861e0_wp
            v0(4) = -0.001758495e0_wp
        elseif( ii .eq. 22 )then
            rap(1) = 1.93206106e0_wp ! rOH=1.0224 using rOH=0.964A
            rap(2) = 2.49443845e0_wp ! rOH=1.32 A
            rap(3) = 4.27078097e0_wp ! rOH=2.26 A
            rap(4) = 9.44863047e0_wp ! rOH=5.00 A
            v0(1)  = -0.008940821e0_wp
            v0(2)  = -0.00976278e0_wp
            v0(3)  = -0.00976278e0_wp
            v0(4)  = -0.00976278e0_wp
        elseif( ii .eq. 33 )then
            rap(1) = 1.93206106e0_wp ! rOH=1.0224 using rOH=0.964A
            rap(2) = 2.49443845e0_wp ! rOH=1.32 A
            rap(3) = 4.27078097e0_wp ! rOH=2.26 A
            rap(4) = 9.44863047e0_wp ! rOH=5.00 A
            v0(1)  = -0.012423316e0_wp
            v0(2)  = -0.015044491e0_wp
            v0(3)  = -0.014671263e0_wp
            v0(4)  = -0.01433271e0_wp
        endif
        call tent(nap, rap, r, ir, t, dtdr)
        do i=1, nqs
            qs(i) = qtc(i)
        enddo
        do i=1, nqb
            qb(i)  = qtc(i+nqs)
            sqb(i) = sin(qb(i))
        enddo
        do i=1, nqt
            qt(i) = qtc(i+nqs+nqb)
        enddo
        ! Initialize hessian elements of stretchings, bendings, and torsions
        do i=1, nqs
            do j=1, i
                fs1(i,j) = 0.0e0_wp
                fs2(i,j) = 0.0e0_wp
            enddo
        enddo
        do i=1, nqb
            do j=1, i
                fb1(i,j) = 0.0e0_wp
                fb2(i,j) = 0.0e0_wp
            enddo
        enddo
        do i=1, nqt
            do j=1, i
                ft1(i,j) = 0.0e0_wp
                ft2(i,j) = 0.0e0_wp
            enddo
        enddo
        do i=1, nqb
            do j=1, nqs
                fbs1(i,j) = 0.0e0_wp
                fbs2(i,j) = 0.0e0_wp
            enddo
        enddo
        ! Get force constant matrices for strecthings (FS) and bending (FB)
        ! and torsions (FT)
        if( ir .eq. 1 )then
            u0 = v0(1)
            do i=1, nqs
                !          qs1(i)=qs(i)-qse(i,1)
                qse1(i) = qse(i,1)
                do j=1, i
                    fs1(i,j) = hs(i,j,1)*qse1(i)*qse1(j)
                enddo
            enddo
            do i=1, nqb
                !          qb1(i)=qb(i)-qbe(i,1)
                ! cos(theta0)-cos(theta) as variables
                qbe1(i) = qbe(i,1)
                qb1(i)  = cos(qbe1(i))-cos(qb(i))
                do j=1, i
                    fb1(i,j) = hb(i,j,1)/(sin(qbe1(i))*sin(qbe1(j)))
                enddo
            enddo
            do i=1, nqt
                qt1(i)  = qt(i)-qte(i,1)
                sqt1(i) = sin(qt1(i))
                cqt1(i) = cos(qt1(i))
                do j=1, i
                    ft1(i,j) = ht(i,j,1)
                enddo
            enddo
            ! bend-stretch couplings
            do i=1, nqb
                do j=1, nqs
                    fbs1(i,j) = hbs(i,j,1)*qse1(j)/sin(qbe1(i))
                enddo
            enddo
        elseif( ir .eq. nap+1 )then
            u0 = v0(nap)
            do i=1, nqs
                qse1(i) = qse(i,nap)
                do j=1, i
                    fs1(i,j) = hs(i,j,nap)*qse1(i)*qse1(j)
                enddo
            enddo
            do i=1, nqb
                ! cos(theta0)-cos(theta) as variables
                qbe1(i) = qbe(i,nap)
                qb1(i)  = cos(qbe1(i))-cos(qb(i))
                do j=1, i
                    fb1(i,j) = hb(i,j,nap)/(sin(qbe1(i))*sin(qbe1(j)))
                enddo
            enddo
            do i=1, nqt
                qt1(i)  = qt(i)-qte(i,nap)
                sqt1(i) = sin(qt1(i))
                cqt1(i) = cos(qt1(i))
                do j=1, i
                    ft1(i,j) = ht(i,j,nap)
                enddo
            enddo
            ! bend-stretch couplings
            do i=1, nqb
                do j=1, nqs
                    fbs1(i,j) = hbs(i,j,nap)*qse1(j)/sin(qbe1(i))
                enddo
            enddo
        elseif( (ir .gt. 1) .and. (ir .lt. nap+1) )then
            u01 = v0(ir-1)
            u02 = v0(ir)
            do i=1, nqs
                qse1(i) = qse(i,ir-1)
                qse2(i) = qse(i,ir)
                do j=1, i
                    fs1(i,j) = hs(i,j,ir-1)*qse1(i)*qse1(j)
                    fs2(i,j) = hs(i,j,ir)*qse2(i)*qse2(j)
                enddo
            enddo
            do i=1, nqb
                qbe1(i) = qbe(i,ir-1)
                qbe2(i) = qbe(i,ir)
                qb1(i)  = cos(qbe1(i))-cos(qb(i))
                qb2(i)  = cos(qbe2(i))-cos(qb(i))
                do j=1, i
                    fb1(i,j) = hb(i,j,ir-1)/(sin(qbe1(i))*sin(qbe1(j)))
                    fb2(i,j) = hb(i,j,ir)/(sin(qbe2(i))*sin(qbe2(j)))
                enddo
            enddo
            do i=1, nqt
                qt1(i)  = qt(i)-qte(i,ir-1)
                qt2(i)  = qt(i)-qte(i,ir)
                sqt1(i) = sin(qt1(i))
                sqt2(i) = sin(qt2(i))
                cqt1(i) = cos(qt1(i))
                cqt2(i) = cos(qt2(i))
                do j=1, i
                    ft1(i,j) = ht(i,j,ir-1)
                    ft2(i,j) = ht(i,j,ir)
                enddo
            enddo
            ! bend-stretch couplings
            do i=1, nqb
                do j=1, nqs
                    fbs1(i,j) = hbs(i,j,ir-1)*qse1(j)/sin(qbe1(i))
                    fbs2(i,j) = hbs(i,j,ir)*qse2(j)/sin(qbe2(i))
                enddo
            enddo
        endif
        ! Get the upper parts of force matrix by symmetry
        do i=1, nqs
            do j=i+1, nqs
                fs1(i,j) = fs1(j,i)
                fs2(i,j) = fs2(j,i)
            enddo
        enddo
        do i=1, nqb
            do j=i+1, nqb
                fb1(i,j) = fb1(j,i)
                fb2(i,j) = fb2(j,i)
            enddo
        enddo
        do i=1, nqt
            do j=i+1, nqt
                ft1(i,j) = ft1(j,i)
                ft2(i,j) = ft2(j,i)
            enddo
        enddo
        !      write(*,'(A,I2)') "ir = ",ir
        ! Calculate potential energy
        if( (ir .eq. 1) .or. (ir .eq. nap+1) )then
            vs  = 0.0e0_wp
            vb  = 0.0e0_wp
            vt  = 0.0e0_wp
            vbs = 0.0e0_wp
            do i=1, nqs
                do j=1, nqs
                    !            vs=vs+0.5d0*FS1(i,j)*(qs(i)-qse1(i))*(qs(j)-qse1(j))
                    vs = vs+0.5e0_wp*fs1(i,j) *(qs(i)-qse1(i))*(qs(j)-qse1(j))/(qs(i)*qs(j))
                enddo
            enddo
            do i=1, nqb
                do j=1, nqb
                    vb = vb+0.5e0_wp*qb1(i)*fb1(i,j)*qb1(j)
                enddo
            enddo
            do i=1, nqt
                do j=1, nqt
                    vt = vt+0.5e0_wp*ft1(i,j)*sqt1(i)*sqt1(j)
                enddo
                vt = vt-0.5e0_wp*ft1(i,i)*sqt1(i)*sqt1(i) +2.0e0_wp*ft1(i,i)*sin(qt1(i)/2.0e0_wp)*sin(qt1(i)/2.0e0_wp)
            enddo
            ! bend-stretch couplings
            do i=1, nqb
                do j=1, nqs
                    !            vbs=vbs+FBS1(i,j)*qb1(i)*(qs(j)-qse1(j))
                    vbs = vbs+fbs1(i,j)*qb1(i)*(qs(j)-qse1(j))/qs(j)
                enddo
            enddo
            ut = u0+vs+vb+vt+vbs
            if( igrad .eq. 1 )then
                !          dvsdr=0.0d0
                !          dvbdr=0.0d0
                !          dvtdr=0.0d0
                !          dvbsdr=0.0d0
                do i=1, nqs
                    dvs(i) = 0.0e0_wp
                    do j=1, nqs
                        !              dvs(i)=dvs(i)+FS1(i,j)*(qs(j)-qse1(j))
                        dvs(i) = dvs(i)+fs1(i,j)*qse1(i) *(qs(j)-qse1(j))/(qs(i)*qs(i)*qs(j))
                    enddo
                enddo
                do i=1, nqb
                    dvb(i) = 0.0e0_wp
                    do j=1, nqb
                        dvb(i) = dvb(i)+sqb(i)*fb1(i,j)*qb1(j)
                    enddo
                enddo
                do i=1, nqt
                    dvt(i) = 0.0e0_wp
                    do j=1, nqt
                        dvt(i) = dvt(i)+ft1(i,j)*cqt1(i)*sqt1(j)
                    enddo
                    dvt(i) = dvt(i)-ft1(i,i)*cqt1(i)*sqt1(i) +2.0e0_wp*ft1(i,i)*cos(qt1(i)/2.0e0_wp)*sin(qt1(i)/ &
                        & 2.0e0_wp)
                enddo
                ! contribution from vbs (dvs and dvb)
                do i=1, nqb
                    dvbsdb(i) = 0.0e0_wp
                    do j=1, nqs
                        !              dvbsdb(i)=dvbsdb(i)+FBS1(i,j)*sqb(i)*(qs(j)-qse1(j))
                        dvbsdb(i) = dvbsdb(i)+fbs1(i,j)*sqb(i) *(qs(j)-qse1(j))/qs(j)
                    enddo
                enddo
                do i=1, nqs
                    dvbsds(i) = 0.0e0_wp
                    do j=1, nqb
                        !              dvbsds(i)=dvbsds(i)+FBS1(j,i)*qb1(j)
                        dvbsds(i) = dvbsds(i)+fbs1(j,i)*qse1(i) *qb1(j)/(qs(i)*qs(i))
                    enddo
                enddo
                !          dutdr=dvsdr+dvbdr+dvtdr+dvbsdr
                dutdr = 0.0e0_wp
                do i=1, nqs
                    dut(i) = dvs(i)+dvbsds(i)
                enddo
                do i=1, nqb
                    dut(nqs+i) = dvb(i)+dvbsdb(i)
                enddo
                do i=1, nqt
                    dut(nqs+nqb+i) = dvt(i)
                enddo
            endif
        else
            vs1  = 0.0e0_wp
            vs2  = 0.0e0_wp
            vb1  = 0.0e0_wp
            vb2  = 0.0e0_wp
            vt1  = 0.0e0_wp
            vt2  = 0.0e0_wp
            vbs1 = 0.0e0_wp
            vbs2 = 0.0e0_wp
            do i=1, nqs
                do j=1, nqs
                    !            vs1=vs1+0.5d0*FS1(i,j)*(qs(i)-qse1(i))*(qs(j)-qse1(j))
                    vs1 = vs1+0.5e0_wp*fs1(i,j) *(qs(i)-qse1(i))*(qs(j)-qse1(j))/(qs(i)*qs(j))
                    !            vs2=vs2+0.5d0*FS2(i,j)*(qs(i)-qse2(i))*(qs(j)-qse2(j))
                    vs2 = vs2+0.5e0_wp*fs2(i,j) *(qs(i)-qse2(i))*(qs(j)-qse2(j))/(qs(i)*qs(j))
                enddo
            enddo
            do i=1, nqb
                do j=1, nqb
                    vb1 = vb1+0.5e0_wp*qb1(i)*fb1(i,j)*qb1(j)
                    vb2 = vb2+0.5e0_wp*qb2(i)*fb2(i,j)*qb2(j)
                enddo
            enddo
            do i=1, nqt
                do j=1, nqt
                    vt1 = vt1+0.5e0_wp*sqt1(i)*ft1(i,j)*sqt1(j)
                    vt2 = vt2+0.5e0_wp*sqt2(i)*ft2(i,j)*sqt2(j)
                enddo
                vt1 = vt1-0.5e0_wp*ft1(i,i)*sqt1(i)*sqt1(i) +2.0e0_wp*ft1(i,i)*sin(qt1(i)/2.0e0_wp)*sin(qt1(i)/2.0e0_wp)
                vt2 = vt2-0.5e0_wp*ft2(i,i)*sqt2(i)*sqt2(i) +2.0e0_wp*ft2(i,i)*sin(qt2(i)/2.0e0_wp)*sin(qt2(i)/2.0e0_wp)
            enddo
            ! bend-stretch couplings
            do i=1, nqb
                do j=1, nqs
                    !            vbs1=vbs1+qb1(i)*FBS1(i,j)*(qs(j)-qse1(j))
                    !            vbs2=vbs2+qb2(i)*FBS2(i,j)*(qs(j)-qse2(j))
                    vbs1 = vbs1+fbs1(i,j)*qb1(i)*(qs(j)-qse1(j))/qs(j)
                    vbs2 = vbs2+fbs2(i,j)*qb2(i)*(qs(j)-qse2(j))/qs(j)
                enddo
            enddo
            u1 = u01+vs1+vb1+vt1+vbs1
            u2 = u02+vs2+vb2+vt2+vbs2
            ut = u1*t(ir-1)+u2*t(ir)
            if( igrad .eq. 1 )then
                dutdr = u1*dtdr(ir-1)+u2*dtdr(ir)
                do i=1, nqs
                    dvs1(i) = 0.0e0_wp
                    dvs2(i) = 0.0e0_wp
                    do j=1, nqs
                        !              dvs1(i)=dvs1(i)+FS1(i,j)*(qs(j)-qse1(j))
                        !              dvs2(i)=dvs2(i)+FS2(i,j)*(qs(j)-qse2(j))
                        dvs1(i) = dvs1(i)+fs1(i,j)*qse1(i) *(qs(j)-qse1(j))/(qs(i)*qs(i)*qs(j))
                        dvs2(i) = dvs2(i)+fs2(i,j)*qse2(i) *(qs(j)-qse2(j))/(qs(i)*qs(i)*qs(j))
                    enddo
                    dvs(i) = dvs1(i)*t(ir-1)+dvs2(i)*t(ir)
                enddo
                do i=1, nqb
                    dvb1(i) = 0.0e0_wp
                    dvb2(i) = 0.0e0_wp
                    do j=1, nqb
                        dvb1(i) = dvb1(i)+sqb(i)*fb1(i,j)*qb1(j)
                        dvb2(i) = dvb2(i)+sqb(i)*fb2(i,j)*qb2(j)
                    enddo
                    dvb(i) = dvb1(i)*t(ir-1)+dvb2(i)*t(ir)
                enddo
                do i=1, nqt
                    dvt1(i) = 0.0e0_wp
                    dvt2(i) = 0.0e0_wp
                    do j=1, nqt
                        dvt1(i) = dvt1(i)+ft1(i,j)*cqt1(i)*sqt1(j)
                        dvt2(i) = dvt2(i)+ft2(i,j)*cqt2(i)*sqt2(j)
                    enddo
                    dvt1(i) = dvt1(i)-ft1(i,i)*cqt1(i)*sqt1(i) +2.0e0_wp*ft1(i,i)*sin(qt1(i)/2.0e0_wp)*cos(qt1(i)/ &
                        & 2.0e0_wp)
                    dvt2(i) = dvt2(i)-ft2(i,i)*cqt2(i)*sqt2(i) +2.0e0_wp*ft2(i,i)*sin(qt2(i)/2.0e0_wp)*cos(qt2(i)/ &
                        & 2.0e0_wp)
                    dvt(i) = dvt1(i)*t(ir-1)+dvt2(i)*t(ir)
                enddo
                ! contribution from vbs (dvbsdb and dvbsds)
                do i=1, nqb
                    dvbsdb1(i) = 0.0e0_wp
                    dvbsdb2(i) = 0.0e0_wp
                    do j=1, nqs
                        !              dvbsdb1(i)=dvbsdb1(i)+FBS1(i,j)*sqb(i)*(qs(j)-qse1(j))
                        !              dvbsdb2(i)=dvbsdb2(i)+FBS2(i,j)*sqb(i)*(qs(j)-qse2(j))
                        dvbsdb1(i) = dvbsdb1(i)+fbs1(i,j)*sqb(i) *(qs(j)-qse1(j))/qs(j)
                        dvbsdb2(i) = dvbsdb2(i)+fbs2(i,j)*sqb(i) *(qs(j)-qse2(j))/qs(j)
                    enddo
                    dvbsdb(i) = dvbsdb1(i)*t(ir-1)+dvbsdb2(i)*t(ir)
                enddo
                do i=1, nqs
                    dvbsds1(i) = 0.0e0_wp
                    dvbsds2(i) = 0.0e0_wp
                    do j=1, nqb
                        !              dvbsds1(i)=dvbsds1(i)+FBS1(j,i)*qb1(j)
                        !              dvbsds2(i)=dvbsds2(i)+FBS2(j,i)*qb2(j)
                        dvbsds1(i) = dvbsds1(i)+fbs1(j,i)*qse1(i) *qb1(j)/(qs(i)*qs(i))
                        dvbsds2(i) = dvbsds2(i)+fbs2(j,i)*qse2(i) *qb2(j)/(qs(i)*qs(i))
                    enddo
                    dvbsds(i) = dvbsds1(i)*t(ir-1)+dvbsds2(i)*t(ir)
                enddo
                do i=1, nqs
                    dut(i) = dvs(i)+dvbsds(i)
                enddo
                do i=1, nqb
                    dut(nqs+i) = dvb(i)+dvbsdb(i)
                enddo
                do i=1, nqt
                    dut(nqs+nqb+i) = dvt(i)
                enddo
            endif
        endif
    end subroutine evuiit

    !> @brief Evaluate gethess for the PHOH_APRP potential.
    subroutine gethess( nap, hss, nqs, hbb, nqb, htt, nqt, hbs, ii )
        ! Global variables
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: nap
        integer  :: ii
        integer  :: nqs, nqb, nqt
        real(wp) :: hss(nqs,nqs,nap)
        real(wp) :: hbb(nqb,nqb,nap)
        real(wp) :: htt(nqt,nqt,nap)
        real(wp) :: hbs(nqb,nqs,nap)
        if( ii .eq. 11 )then
            ! Hessian elements for anchor point 1
            hss( 1,1,1)  = 0.40759e+00_wp
            hss( 2,1,1)  = 0.17500e-01_wp
            hss( 2,2,1)  = 0.40901e+00_wp
            hss( 3,1,1)  = 0.29920e-01_wp
            hss( 3,2,1)  = 0.35360e-01_wp
            hss( 3,3,1)  = 0.42477e+00_wp
            hss( 4,1,1)  = 0.18000e-01_wp
            hss( 4,2,1)  = -0.26000e-02_wp
            hss( 4,3,1)  = 0.38000e-03_wp
            hss( 4,4,1)  = 0.41176e+00_wp
            hss( 5,1,1)  = 0.45000e-02_wp
            hss( 5,2,1)  = -0.29000e-03_wp
            hss( 5,3,1)  = 0.30900e-02_wp
            hss( 5,4,1)  = 0.51800e-02_wp
            hss( 5,5,1)  = 0.35919e+00_wp
            hss( 6,1,1)  = -0.34400e-02_wp
            hss( 6,2,1)  = 0.82770e-01_wp
            hss( 6,3,1)  = -0.25600e-02_wp
            hss( 6,4,1)  = 0.21010e-01_wp
            hss( 6,5,1)  = -0.76000e-03_wp
            hss( 6,6,1)  = 0.41463e+00_wp
            hss( 7,1,1)  = -0.55000e-03_wp
            hss( 7,2,1)  = -0.34000e-03_wp
            hss( 7,3,1)  = 0.68000e-03_wp
            hss( 7,4,1)  = 0.44000e-02_wp
            hss( 7,5,1)  = 0.48000e-03_wp
            hss( 7,6,1)  = 0.45900e-02_wp
            hss( 7,7,1)  = 0.36428e+00_wp
            hss( 8,1,1)  = 0.82220e-01_wp
            hss( 8,2,1)  = -0.40500e-02_wp
            hss( 8,3,1)  = -0.11000e-03_wp
            hss( 8,4,1)  = -0.33000e-02_wp
            hss( 8,5,1)  = -0.30000e-03_wp
            hss( 8,6,1)  = 0.20930e-01_wp
            hss( 8,7,1)  = -0.31000e-03_wp
            hss( 8,8,1)  = 0.40778e+00_wp
            hss( 9,1,1)  = -0.48000e-03_wp
            hss( 9,2,1)  = -0.49000e-03_wp
            hss( 9,3,1)  = -0.23000e-03_wp
            hss( 9,4,1)  = -0.57000e-03_wp
            hss( 9,5,1)  = 0.20000e-03_wp
            hss( 9,6,1)  = 0.44500e-02_wp
            hss( 9,7,1)  = 0.47000e-03_wp
            hss( 9,8,1)  = 0.44300e-02_wp
            hss( 9,9,1)  = 0.36674e+00_wp
            hss(10,1,1)  = -0.18900e-02_wp
            hss(10,2,1)  = 0.18160e-01_wp
            hss(10,3,1)  = -0.10000e-02_wp
            hss(10,4,1)  = 0.86110e-01_wp
            hss(10,5,1)  = -0.74000e-03_wp
            hss(10,6,1)  = -0.29100e-02_wp
            hss(10,7,1)  = -0.54000e-03_wp
            hss(10,8,1)  = 0.20960e-01_wp
            hss(10,9,1)  = -0.51000e-03_wp
            hss(10,10,1) = 0.41852e+00_wp
            hss(11,1,1)  = -0.43000e-03_wp
            hss(11,2,1)  = -0.50000e-03_wp
            hss(11,3,1)  = 0.78000e-03_wp
            hss(11,4,1)  = -0.50000e-03_wp
            hss(11,5,1)  = 0.11000e-03_wp
            hss(11,6,1)  = -0.28000e-03_wp
            hss(11,7,1)  = 0.20000e-03_wp
            hss(11,8,1)  = 0.45700e-02_wp
            hss(11,9,1)  = 0.47000e-03_wp
            hss(11,10,1) = 0.45000e-02_wp
            hss(11,11,1) = 0.36422e+00_wp
            hss(12,1,1)  = -0.40000e-03_wp
            hss(12,2,1)  = 0.35500e-02_wp
            hss(12,3,1)  = 0.10000e-04_wp
            hss(12,4,1)  = -0.80000e-03_wp
            hss(12,5,1)  = 0.23000e-03_wp
            hss(12,6,1)  = -0.14000e-03_wp
            hss(12,7,1)  = 0.10000e-03_wp
            hss(12,8,1)  = -0.68000e-03_wp
            hss(12,9,1)  = 0.14000e-03_wp
            hss(12,10,1) = 0.42300e-02_wp
            hss(12,11,1) = 0.44000e-03_wp
            hss(12,12,1) = 0.36727e+00_wp
            ! Hessian elements for anchor point 1
            hbb( 1,1,1)  = 0.86570e-01_wp
            hbb( 2,1,1)  = -0.44420e-01_wp
            hbb( 2,2,1)  = 0.13902e+00_wp
            hbb( 3,1,1)  = -0.42150e-01_wp
            hbb( 3,2,1)  = -0.94600e-01_wp
            hbb( 3,3,1)  = 0.13674e+00_wp
            hbb( 4,1,1)  = -0.36220e-01_wp
            hbb( 4,2,1)  = 0.23510e-01_wp
            hbb( 4,3,1)  = 0.12710e-01_wp
            hbb( 4,4,1)  = 0.79210e-01_wp
            hbb( 5,1,1)  = 0.20780e-01_wp
            hbb( 5,2,1)  = -0.13590e-01_wp
            hbb( 5,3,1)  = -0.71800e-02_wp
            hbb( 5,4,1)  = -0.38760e-01_wp
            hbb( 5,5,1)  = 0.77300e-01_wp
            hbb( 6,1,1)  = 0.15440e-01_wp
            hbb( 6,2,1)  = -0.99200e-02_wp
            hbb( 6,3,1)  = -0.55200e-02_wp
            hbb( 6,4,1)  = -0.40440e-01_wp
            hbb( 6,5,1)  = -0.38530e-01_wp
            hbb( 6,6,1)  = 0.78980e-01_wp
            hbb( 7,1,1)  = -0.10540e-01_wp
            hbb( 7,2,1)  = -0.24300e-02_wp
            hbb( 7,3,1)  = 0.12970e-01_wp
            hbb( 7,4,1)  = -0.31280e-01_wp
            hbb( 7,5,1)  = 0.12070e-01_wp
            hbb( 7,6,1)  = 0.19210e-01_wp
            hbb( 7,7,1)  = 0.79110e-01_wp
            hbb( 8,1,1)  = -0.42000e-03_wp
            hbb( 8,2,1)  = 0.59300e-02_wp
            hbb( 8,3,1)  = -0.55100e-02_wp
            hbb( 8,4,1)  = 0.18520e-01_wp
            hbb( 8,5,1)  = -0.65400e-02_wp
            hbb( 8,6,1)  = -0.11980e-01_wp
            hbb( 8,7,1)  = -0.39530e-01_wp
            hbb( 8,8,1)  = 0.79370e-01_wp
            hbb( 9,1,1)  = 0.10960e-01_wp
            hbb( 9,2,1)  = -0.35000e-02_wp
            hbb( 9,3,1)  = -0.74600e-02_wp
            hbb( 9,4,1)  = 0.12760e-01_wp
            hbb( 9,5,1)  = -0.55300e-02_wp
            hbb( 9,6,1)  = -0.72300e-02_wp
            hbb( 9,7,1)  = -0.39590e-01_wp
            hbb( 9,8,1)  = -0.39850e-01_wp
            hbb( 9,9,1)  = 0.79440e-01_wp
            hbb(10,1,1)  = 0.69200e-02_wp
            hbb(10,2,1)  = -0.37400e-02_wp
            hbb(10,3,1)  = -0.31900e-02_wp
            hbb(10,4,1)  = -0.97700e-02_wp
            hbb(10,5,1)  = 0.10500e-01_wp
            hbb(10,6,1)  = -0.72000e-03_wp
            hbb(10,7,1)  = -0.33920e-01_wp
            hbb(10,8,1)  = 0.13770e-01_wp
            hbb(10,9,1)  = 0.20160e-01_wp
            hbb(10,10,1) = 0.80580e-01_wp
            hbb(11,1,1)  = -0.33900e-02_wp
            hbb(11,2,1)  = 0.23200e-02_wp
            hbb(11,3,1)  = 0.10700e-02_wp
            hbb(11,4,1)  = -0.78000e-03_wp
            hbb(11,5,1)  = -0.37400e-02_wp
            hbb(11,6,1)  = 0.45200e-02_wp
            hbb(11,7,1)  = 0.19900e-01_wp
            hbb(11,8,1)  = -0.73600e-02_wp
            hbb(11,9,1)  = -0.12550e-01_wp
            hbb(11,10,1) = -0.40290e-01_wp
            hbb(11,11,1) = 0.79280e-01_wp
            hbb(12,1,1)  = -0.35300e-02_wp
            hbb(12,2,1)  = 0.14200e-02_wp
            hbb(12,3,1)  = 0.21200e-02_wp
            hbb(12,4,1)  = 0.10550e-01_wp
            hbb(12,5,1)  = -0.67600e-02_wp
            hbb(12,6,1)  = -0.37900e-02_wp
            hbb(12,7,1)  = 0.14020e-01_wp
            hbb(12,8,1)  = -0.64100e-02_wp
            hbb(12,9,1)  = -0.76100e-02_wp
            hbb(12,10,1) = -0.40290e-01_wp
            hbb(12,11,1) = -0.38990e-01_wp
            hbb(12,12,1) = 0.79280e-01_wp
            hbb(13,1,1)  = -0.10320e-01_wp
            hbb(13,2,1)  = 0.13110e-01_wp
            hbb(13,3,1)  = -0.27900e-02_wp
            hbb(13,4,1)  = 0.52200e-02_wp
            hbb(13,5,1)  = -0.23700e-02_wp
            hbb(13,6,1)  = -0.28500e-02_wp
            hbb(13,7,1)  = -0.86400e-02_wp
            hbb(13,8,1)  = 0.10140e-01_wp
            hbb(13,9,1)  = -0.15000e-02_wp
            hbb(13,10,1) = -0.34100e-01_wp
            hbb(13,11,1) = 0.14060e-01_wp
            hbb(13,12,1) = 0.20050e-01_wp
            hbb(13,13,1) = 0.79730e-01_wp
            hbb(14,1,1)  = 0.10860e-01_wp
            hbb(14,2,1)  = -0.76500e-02_wp
            hbb(14,3,1)  = -0.32000e-02_wp
            hbb(14,4,1)  = -0.26900e-02_wp
            hbb(14,5,1)  = 0.11000e-02_wp
            hbb(14,6,1)  = 0.16000e-02_wp
            hbb(14,7,1)  = -0.15200e-02_wp
            hbb(14,8,1)  = -0.35000e-02_wp
            hbb(14,9,1)  = 0.50200e-02_wp
            hbb(14,10,1) = 0.20210e-01_wp
            hbb(14,11,1) = -0.76100e-02_wp
            hbb(14,12,1) = -0.12600e-01_wp
            hbb(14,13,1) = -0.39850e-01_wp
            hbb(14,14,1) = 0.79680e-01_wp
            hbb(15,1,1)  = -0.54000e-03_wp
            hbb(15,2,1)  = -0.54500e-02_wp
            hbb(15,3,1)  = 0.59900e-02_wp
            hbb(15,4,1)  = -0.25300e-02_wp
            hbb(15,5,1)  = 0.12700e-02_wp
            hbb(15,6,1)  = 0.12600e-02_wp
            hbb(15,7,1)  = 0.10150e-01_wp
            hbb(15,8,1)  = -0.66300e-02_wp
            hbb(15,9,1)  = -0.35200e-02_wp
            hbb(15,10,1) = 0.13890e-01_wp
            hbb(15,11,1) = -0.64500e-02_wp
            hbb(15,12,1) = -0.74400e-02_wp
            hbb(15,13,1) = -0.39880e-01_wp
            hbb(15,14,1) = -0.39830e-01_wp
            hbb(15,15,1) = 0.79710e-01_wp
            hbb(16,1,1)  = -0.36420e-01_wp
            hbb(16,2,1)  = 0.13970e-01_wp
            hbb(16,3,1)  = 0.22450e-01_wp
            hbb(16,4,1)  = -0.71600e-02_wp
            hbb(16,5,1)  = -0.22100e-02_wp
            hbb(16,6,1)  = 0.93700e-02_wp
            hbb(16,7,1)  = 0.52700e-02_wp
            hbb(16,8,1)  = -0.24800e-02_wp
            hbb(16,9,1)  = -0.27900e-02_wp
            hbb(16,10,1) = -0.97100e-02_wp
            hbb(16,11,1) = 0.10500e-01_wp
            hbb(16,12,1) = -0.79000e-03_wp
            hbb(16,13,1) = -0.31890e-01_wp
            hbb(16,14,1) = 0.12990e-01_wp
            hbb(16,15,1) = 0.18900e-01_wp
            hbb(16,16,1) = 0.79900e-01_wp
            hbb(17,1,1)  = 0.22220e-01_wp
            hbb(17,2,1)  = -0.60200e-02_wp
            hbb(17,3,1)  = -0.16200e-01_wp
            hbb(17,4,1)  = -0.27800e-02_wp
            hbb(17,5,1)  = 0.56900e-02_wp
            hbb(17,6,1)  = -0.29000e-02_wp
            hbb(17,7,1)  = -0.25500e-02_wp
            hbb(17,8,1)  = 0.13300e-02_wp
            hbb(17,9,1)  = 0.12200e-02_wp
            hbb(17,10,1) = 0.10510e-01_wp
            hbb(17,11,1) = -0.68200e-02_wp
            hbb(17,12,1) = -0.37000e-02_wp
            hbb(17,13,1) = 0.12370e-01_wp
            hbb(17,14,1) = -0.57100e-02_wp
            hbb(17,15,1) = -0.66600e-02_wp
            hbb(17,16,1) = -0.39780e-01_wp
            hbb(17,17,1) = 0.76100e-01_wp
            hbb(18,1,1)  = 0.14200e-01_wp
            hbb(18,2,1)  = -0.79500e-02_wp
            hbb(18,3,1)  = -0.62400e-02_wp
            hbb(18,4,1)  = 0.99400e-02_wp
            hbb(18,5,1)  = -0.34800e-02_wp
            hbb(18,6,1)  = -0.64600e-02_wp
            hbb(18,7,1)  = -0.27200e-02_wp
            hbb(18,8,1)  = 0.11500e-02_wp
            hbb(18,9,1)  = 0.15700e-02_wp
            hbb(18,10,1) = -0.80000e-03_wp
            hbb(18,11,1) = -0.36800e-02_wp
            hbb(18,12,1) = 0.44900e-02_wp
            hbb(18,13,1) = 0.19510e-01_wp
            hbb(18,14,1) = -0.72700e-02_wp
            hbb(18,15,1) = -0.12240e-01_wp
            hbb(18,16,1) = -0.40120e-01_wp
            hbb(18,17,1) = -0.36330e-01_wp
            hbb(18,18,1) = 0.76450e-01_wp
            ! Hessian elements for anchor point 1
            htt( 1,1,1)  = 0.92400e-02_wp
            htt( 2,1,1)  = 0.38700e-02_wp
            htt( 2,2,1)  = 0.84000e-02_wp
            htt( 3,1,1)  = 0.67000e-03_wp
            htt( 3,2,1)  = -0.40800e-02_wp
            htt( 3,3,1)  = 0.11550e-01_wp
            htt( 4,1,1)  = -0.47000e-02_wp
            htt( 4,2,1)  = 0.44000e-03_wp
            htt( 4,3,1)  = 0.68000e-02_wp
            htt( 4,4,1)  = 0.11950e-01_wp
            htt( 5,1,1)  = -0.62200e-02_wp
            htt( 5,2,1)  = -0.44300e-02_wp
            htt( 5,3,1)  = 0.20700e-02_wp
            htt( 5,4,1)  = 0.38700e-02_wp
            htt( 5,5,1)  = 0.89900e-02_wp
            htt( 6,1,1)  = -0.43000e-02_wp
            htt( 6,2,1)  = -0.29300e-02_wp
            htt( 6,3,1)  = 0.28800e-02_wp
            htt( 6,4,1)  = 0.42500e-02_wp
            htt( 6,5,1)  = 0.35700e-02_wp
            htt( 6,6,1)  = 0.80900e-02_wp
            htt( 7,1,1)  = 0.19000e-02_wp
            htt( 7,2,1)  = 0.31100e-02_wp
            htt( 7,3,1)  = -0.82500e-02_wp
            htt( 7,4,1)  = -0.70400e-02_wp
            htt( 7,5,1)  = 0.11300e-02_wp
            htt( 7,6,1)  = -0.32500e-02_wp
            htt( 7,7,1)  = 0.10750e-01_wp
            htt( 8,1,1)  = 0.38300e-02_wp
            htt( 8,2,1)  = 0.46100e-02_wp
            htt( 8,3,1)  = -0.74400e-02_wp
            htt( 8,4,1)  = -0.66600e-02_wp
            htt( 8,5,1)  = -0.43000e-02_wp
            htt( 8,6,1)  = 0.12800e-02_wp
            htt( 8,7,1)  = 0.63800e-02_wp
            htt( 8,8,1)  = 0.11960e-01_wp
            htt( 9,1,1)  = -0.63700e-02_wp
            htt( 9,2,1)  = -0.83000e-03_wp
            htt( 9,3,1)  = -0.30300e-02_wp
            htt( 9,4,1)  = 0.25100e-02_wp
            htt( 9,5,1)  = 0.48000e-03_wp
            htt( 9,6,1)  = 0.20600e-02_wp
            htt( 9,7,1)  = -0.26900e-02_wp
            htt( 9,8,1)  = -0.11100e-02_wp
            htt( 9,9,1)  = 0.94700e-02_wp
            htt(10,1,1)  = -0.40200e-02_wp
            htt(10,2,1)  = 0.76000e-03_wp
            htt(10,3,1)  = -0.16700e-02_wp
            htt(10,4,1)  = 0.31200e-02_wp
            htt(10,5,1)  = 0.20200e-02_wp
            htt(10,6,1)  = 0.15300e-02_wp
            htt(10,7,1)  = -0.22000e-03_wp
            htt(10,8,1)  = -0.70000e-03_wp
            htt(10,9,1)  = 0.31100e-02_wp
            htt(10,10,1) = 0.83700e-02_wp
            htt(11,1,1)  = -0.98000e-03_wp
            htt(11,2,1)  = -0.53700e-02_wp
            htt(11,3,1)  = 0.17300e-02_wp
            htt(11,4,1)  = -0.26500e-02_wp
            htt(11,5,1)  = -0.13300e-02_wp
            htt(11,6,1)  = 0.68000e-03_wp
            htt(11,7,1)  = -0.39000e-02_wp
            htt(11,8,1)  = -0.19000e-02_wp
            htt(11,9,1)  = 0.39100e-02_wp
            htt(11,10,1) = -0.16900e-02_wp
            htt(11,11,1) = 0.83100e-02_wp
            htt(12,1,1)  = 0.13700e-02_wp
            htt(12,2,1)  = -0.37800e-02_wp
            htt(12,3,1)  = 0.31000e-02_wp
            htt(12,4,1)  = -0.20400e-02_wp
            htt(12,5,1)  = 0.22000e-03_wp
            htt(12,6,1)  = 0.16000e-03_wp
            htt(12,7,1)  = -0.14300e-02_wp
            htt(12,8,1)  = -0.14900e-02_wp
            htt(12,9,1)  = -0.24600e-02_wp
            htt(12,10,1) = 0.35600e-02_wp
            htt(12,11,1) = 0.27100e-02_wp
            htt(12,12,1) = 0.87300e-02_wp
            htt(13,1,1)  = 0.53000e-03_wp
            htt(13,2,1)  = -0.16100e-02_wp
            htt(13,3,1)  = 0.26100e-02_wp
            htt(13,4,1)  = 0.47000e-03_wp
            htt(13,5,1)  = 0.24000e-02_wp
            htt(13,6,1)  = 0.89000e-03_wp
            htt(13,7,1)  = 0.43000e-03_wp
            htt(13,8,1)  = -0.10800e-02_wp
            htt(13,9,1)  = -0.66000e-02_wp
            htt(13,10,1) = -0.21000e-03_wp
            htt(13,11,1) = -0.44500e-02_wp
            htt(13,12,1) = 0.19400e-02_wp
            htt(13,13,1) = 0.96100e-02_wp
            htt(14,1,1)  = 0.20400e-02_wp
            htt(14,2,1)  = 0.52000e-03_wp
            htt(14,3,1)  = 0.14800e-02_wp
            htt(14,4,1)  = -0.40000e-04_wp
            htt(14,5,1)  = 0.97000e-03_wp
            htt(14,6,1)  = -0.11900e-02_wp
            htt(14,7,1)  = 0.15000e-02_wp
            htt(14,8,1)  = -0.66000e-03_wp
            htt(14,9,1)  = -0.44200e-02_wp
            htt(14,10,1) = 0.12200e-02_wp
            htt(14,11,1) = -0.29000e-02_wp
            htt(14,12,1) = 0.27500e-02_wp
            htt(14,13,1) = 0.37500e-02_wp
            htt(14,14,1) = 0.85800e-02_wp
            htt(15,1,1)  = -0.18400e-02_wp
            htt(15,2,1)  = -0.32100e-02_wp
            htt(15,3,1)  = 0.12300e-02_wp
            htt(15,4,1)  = -0.15000e-03_wp
            htt(15,5,1)  = 0.85000e-03_wp
            htt(15,6,1)  = 0.14200e-02_wp
            htt(15,7,1)  = -0.20600e-02_wp
            htt(15,8,1)  = -0.14900e-02_wp
            htt(15,9,1)  = -0.19000e-03_wp
            htt(15,10,1) = -0.55100e-02_wp
            htt(15,11,1) = 0.11900e-02_wp
            htt(15,12,1) = -0.41300e-02_wp
            htt(15,13,1) = 0.31600e-02_wp
            htt(15,14,1) = -0.19400e-02_wp
            htt(15,15,1) = 0.85200e-02_wp
            htt(16,1,1)  = -0.33000e-03_wp
            htt(16,2,1)  = -0.10800e-02_wp
            htt(16,3,1)  = 0.10000e-03_wp
            htt(16,4,1)  = -0.66000e-03_wp
            htt(16,5,1)  = -0.59000e-03_wp
            htt(16,6,1)  = -0.66000e-03_wp
            htt(16,7,1)  = -0.99000e-03_wp
            htt(16,8,1)  = -0.10700e-02_wp
            htt(16,9,1)  = 0.19900e-02_wp
            htt(16,10,1) = -0.40700e-02_wp
            htt(16,11,1) = 0.27500e-02_wp
            htt(16,12,1) = -0.33200e-02_wp
            htt(16,13,1) = -0.26900e-02_wp
            htt(16,14,1) = 0.28900e-02_wp
            htt(16,15,1) = 0.34200e-02_wp
            htt(16,16,1) = 0.90000e-02_wp
            htt(17,1,1)  = 0.25300e-02_wp
            htt(17,2,1)  = 0.10400e-02_wp
            htt(17,3,1)  = 0.19000e-03_wp
            htt(17,4,1)  = -0.13000e-02_wp
            htt(17,5,1)  = 0.43000e-03_wp
            htt(17,6,1)  = -0.16400e-02_wp
            htt(17,7,1)  = 0.26600e-02_wp
            htt(17,8,1)  = 0.59000e-03_wp
            htt(17,9,1)  = 0.61000e-03_wp
            htt(17,10,1) = -0.18300e-02_wp
            htt(17,11,1) = 0.21000e-02_wp
            htt(17,12,1) = -0.33000e-03_wp
            htt(17,13,1) = -0.66400e-02_wp
            htt(17,14,1) = -0.70000e-03_wp
            htt(17,15,1) = -0.41800e-02_wp
            htt(17,16,1) = 0.17600e-02_wp
            htt(17,17,1) = 0.96800e-02_wp
            htt(18,1,1)  = 0.74000e-03_wp
            htt(18,2,1)  = 0.13600e-02_wp
            htt(18,3,1)  = -0.19200e-02_wp
            htt(18,4,1)  = -0.13000e-02_wp
            htt(18,5,1)  = -0.16100e-02_wp
            htt(18,6,1)  = -0.30500e-02_wp
            htt(18,7,1)  = 0.91000e-03_wp
            htt(18,8,1)  = -0.53000e-03_wp
            htt(18,9,1)  = 0.21300e-02_wp
            htt(18,10,1) = 0.45000e-03_wp
            htt(18,11,1) = 0.15100e-02_wp
            htt(18,12,1) = -0.17000e-03_wp
            htt(18,13,1) = -0.40600e-02_wp
            htt(18,14,1) = 0.10400e-02_wp
            htt(18,15,1) = -0.23700e-02_wp
            htt(18,16,1) = 0.27300e-02_wp
            htt(18,17,1) = 0.31700e-02_wp
            htt(18,18,1) = 0.85000e-02_wp
            htt(19,1,1)  = 0.10300e-02_wp
            htt(19,2,1)  = -0.10900e-02_wp
            htt(19,3,1)  = 0.13200e-02_wp
            htt(19,4,1)  = -0.79000e-03_wp
            htt(19,5,1)  = 0.18700e-02_wp
            htt(19,6,1)  = 0.45000e-03_wp
            htt(19,7,1)  = 0.15900e-02_wp
            htt(19,8,1)  = 0.17000e-03_wp
            htt(19,9,1)  = -0.15700e-02_wp
            htt(19,10,1) = -0.32700e-02_wp
            htt(19,11,1) = 0.55000e-03_wp
            htt(19,12,1) = -0.11400e-02_wp
            htt(19,13,1) = -0.78000e-03_wp
            htt(19,14,1) = -0.55300e-02_wp
            htt(19,15,1) = 0.93000e-03_wp
            htt(19,16,1) = -0.38200e-02_wp
            htt(19,17,1) = 0.37300e-02_wp
            htt(19,18,1) = -0.19300e-02_wp
            htt(19,19,1) = 0.84800e-02_wp
            htt(20,1,1)  = -0.77000e-03_wp
            htt(20,2,1)  = -0.77000e-03_wp
            htt(20,3,1)  = -0.79000e-03_wp
            htt(20,4,1)  = -0.79000e-03_wp
            htt(20,5,1)  = -0.18000e-03_wp
            htt(20,6,1)  = -0.97000e-03_wp
            htt(20,7,1)  = -0.16000e-03_wp
            htt(20,8,1)  = -0.95000e-03_wp
            htt(20,9,1)  = -0.50000e-04_wp
            htt(20,10,1) = -0.98000e-03_wp
            htt(20,11,1) = -0.50000e-04_wp
            htt(20,12,1) = -0.98000e-03_wp
            htt(20,13,1) = 0.18000e-02_wp
            htt(20,14,1) = -0.37900e-02_wp
            htt(20,15,1) = 0.27400e-02_wp
            htt(20,16,1) = -0.28500e-02_wp
            htt(20,17,1) = -0.27700e-02_wp
            htt(20,18,1) = 0.34000e-02_wp
            htt(20,19,1) = 0.28200e-02_wp
            htt(20,20,1) = 0.89900e-02_wp
            htt(21,1,1)  = 0.30000e-03_wp
            htt(21,2,1)  = 0.19600e-02_wp
            htt(21,3,1)  = -0.25100e-02_wp
            htt(21,4,1)  = -0.85000e-03_wp
            htt(21,5,1)  = -0.60900e-02_wp
            htt(21,6,1)  = -0.58000e-03_wp
            htt(21,7,1)  = -0.34300e-02_wp
            htt(21,8,1)  = 0.20800e-02_wp
            htt(21,9,1)  = 0.24300e-02_wp
            htt(21,10,1) = 0.93000e-03_wp
            htt(21,11,1) = 0.76000e-03_wp
            htt(21,12,1) = -0.74000e-03_wp
            htt(21,13,1) = 0.68000e-03_wp
            htt(21,14,1) = -0.16400e-02_wp
            htt(21,15,1) = 0.21900e-02_wp
            htt(21,16,1) = -0.13000e-03_wp
            htt(21,17,1) = -0.66000e-02_wp
            htt(21,18,1) = -0.36000e-03_wp
            htt(21,19,1) = -0.42700e-02_wp
            htt(21,20,1) = 0.19700e-02_wp
            htt(21,21,1) = 0.92800e-02_wp
            htt(22,1,1)  = -0.16800e-02_wp
            htt(22,2,1)  = 0.42000e-03_wp
            htt(22,3,1)  = -0.33400e-02_wp
            htt(22,4,1)  = -0.12400e-02_wp
            htt(22,5,1)  = -0.50000e-03_wp
            htt(22,6,1)  = -0.52300e-02_wp
            htt(22,7,1)  = 0.10700e-02_wp
            htt(22,8,1)  = -0.36600e-02_wp
            htt(22,9,1)  = 0.80000e-03_wp
            htt(22,10,1) = 0.14300e-02_wp
            htt(22,11,1) = -0.13000e-02_wp
            htt(22,12,1) = -0.68000e-03_wp
            htt(22,13,1) = 0.22300e-02_wp
            htt(22,14,1) = 0.58000e-03_wp
            htt(22,15,1) = 0.16000e-02_wp
            htt(22,16,1) = -0.60000e-04_wp
            htt(22,17,1) = -0.44700e-02_wp
            htt(22,18,1) = 0.11200e-02_wp
            htt(22,19,1) = -0.28100e-02_wp
            htt(22,20,1) = 0.27800e-02_wp
            htt(22,21,1) = 0.36100e-02_wp
            htt(22,22,1) = 0.84800e-02_wp
            htt(23,1,1)  = 0.20800e-02_wp
            htt(23,2,1)  = 0.16400e-02_wp
            htt(23,3,1)  = -0.41000e-03_wp
            htt(23,4,1)  = -0.85000e-03_wp
            htt(23,5,1)  = -0.40600e-02_wp
            htt(23,6,1)  = 0.83000e-03_wp
            htt(23,7,1)  = -0.16900e-02_wp
            htt(23,8,1)  = 0.31900e-02_wp
            htt(23,9,1)  = 0.92000e-03_wp
            htt(23,10,1) = -0.13400e-02_wp
            htt(23,11,1) = 0.13500e-02_wp
            htt(23,12,1) = -0.90000e-03_wp
            htt(23,13,1) = -0.18800e-02_wp
            htt(23,14,1) = -0.33700e-02_wp
            htt(23,15,1) = 0.39000e-03_wp
            htt(23,16,1) = -0.11000e-02_wp
            htt(23,17,1) = -0.14000e-03_wp
            htt(23,18,1) = -0.56500e-02_wp
            htt(23,19,1) = 0.13500e-02_wp
            htt(23,20,1) = -0.41600e-02_wp
            htt(23,21,1) = 0.30900e-02_wp
            htt(23,22,1) = -0.19500e-02_wp
            htt(23,23,1) = 0.85600e-02_wp
            htt(24,1,1)  = 0.10000e-03_wp
            htt(24,2,1)  = 0.11000e-03_wp
            htt(24,3,1)  = -0.12500e-02_wp
            htt(24,4,1)  = -0.12400e-02_wp
            htt(24,5,1)  = 0.15300e-02_wp
            htt(24,6,1)  = -0.38300e-02_wp
            htt(24,7,1)  = 0.28100e-02_wp
            htt(24,8,1)  = -0.25500e-02_wp
            htt(24,9,1)  = -0.71000e-03_wp
            htt(24,10,1) = -0.84000e-03_wp
            htt(24,11,1) = -0.71000e-03_wp
            htt(24,12,1) = -0.84000e-03_wp
            htt(24,13,1) = -0.32000e-03_wp
            htt(24,14,1) = -0.11500e-02_wp
            htt(24,15,1) = -0.20000e-03_wp
            htt(24,16,1) = -0.10200e-02_wp
            htt(24,17,1) = 0.19900e-02_wp
            htt(24,18,1) = -0.41700e-02_wp
            htt(24,19,1) = 0.28200e-02_wp
            htt(24,20,1) = -0.33500e-02_wp
            htt(24,21,1) = -0.25900e-02_wp
            htt(24,22,1) = 0.29200e-02_wp
            htt(24,23,1) = 0.35300e-02_wp
            htt(24,24,1) = 0.90500e-02_wp
            ! Hessian elements for anchor point 1
            hbs( 1,1,1)  = -0.36310e-01_wp
            hbs( 1,2,1)  = -0.36690e-01_wp
            hbs( 1,3,1)  = -0.36890e-01_wp
            hbs( 1,4,1)  = -0.11400e-01_wp
            hbs( 1,5,1)  = 0.55900e-02_wp
            hbs( 1,6,1)  = 0.45710e-01_wp
            hbs( 1,7,1)  = -0.44000e-03_wp
            hbs( 1,8,1)  = 0.44410e-01_wp
            hbs( 1,9,1)  = 0.60000e-04_wp
            hbs( 1,10,1) = -0.10640e-01_wp
            hbs( 1,11,1) = -0.46000e-03_wp
            hbs( 1,12,1) = 0.55900e-02_wp
            hbs( 2,1,1)  = 0.42410e-01_wp
            hbs( 2,2,1)  = -0.54500e-02_wp
            hbs( 2,3,1)  = 0.11600e-01_wp
            hbs( 2,4,1)  = -0.17300e-02_wp
            hbs( 2,5,1)  = -0.28000e-02_wp
            hbs( 2,6,1)  = -0.19690e-01_wp
            hbs( 2,7,1)  = 0.11300e-02_wp
            hbs( 2,8,1)  = -0.26090e-01_wp
            hbs( 2,9,1)  = -0.18000e-03_wp
            hbs( 2,10,1) = 0.11810e-01_wp
            hbs( 2,11,1) = -0.69000e-03_wp
            hbs( 2,12,1) = -0.11800e-02_wp
            hbs( 3,1,1)  = -0.61000e-02_wp
            hbs( 3,2,1)  = 0.42140e-01_wp
            hbs( 3,3,1)  = 0.25290e-01_wp
            hbs( 3,4,1)  = 0.13130e-01_wp
            hbs( 3,5,1)  = -0.27900e-02_wp
            hbs( 3,6,1)  = -0.26020e-01_wp
            hbs( 3,7,1)  = -0.68000e-03_wp
            hbs( 3,8,1)  = -0.18310e-01_wp
            hbs( 3,9,1)  = 0.12000e-03_wp
            hbs( 3,10,1) = -0.11700e-02_wp
            hbs( 3,11,1) = 0.11600e-02_wp
            hbs( 3,12,1) = -0.44100e-02_wp
            hbs( 4,1,1)  = -0.29210e-01_wp
            hbs( 4,2,1)  = -0.82700e-02_wp
            hbs( 4,3,1)  = 0.16570e-01_wp
            hbs( 4,4,1)  = -0.31510e-01_wp
            hbs( 4,5,1)  = -0.97200e-02_wp
            hbs( 4,6,1)  = -0.12980e-01_wp
            hbs( 4,7,1)  = 0.48900e-02_wp
            hbs( 4,8,1)  = 0.44650e-01_wp
            hbs( 4,9,1)  = -0.64000e-03_wp
            hbs( 4,10,1) = 0.43490e-01_wp
            hbs( 4,11,1) = 0.19000e-03_wp
            hbs( 4,12,1) = -0.10200e-02_wp
            hbs( 5,1,1)  = 0.27310e-01_wp
            hbs( 5,2,1)  = 0.47900e-02_wp
            hbs( 5,3,1)  = -0.10360e-01_wp
            hbs( 5,4,1)  = 0.23300e-02_wp
            hbs( 5,5,1)  = 0.65800e-02_wp
            hbs( 5,6,1)  = 0.59300e-02_wp
            hbs( 5,7,1)  = -0.19100e-02_wp
            hbs( 5,8,1)  = -0.19800e-01_wp
            hbs( 5,9,1)  = -0.33000e-03_wp
            hbs( 5,10,1) = -0.24210e-01_wp
            hbs( 5,11,1) = -0.10000e-03_wp
            hbs( 5,12,1) = 0.10100e-02_wp
            hbs( 6,1,1)  = 0.19000e-02_wp
            hbs( 6,2,1)  = 0.34900e-02_wp
            hbs( 6,3,1)  = -0.62000e-02_wp
            hbs( 6,4,1)  = 0.29170e-01_wp
            hbs( 6,5,1)  = 0.31400e-02_wp
            hbs( 6,6,1)  = 0.70400e-02_wp
            hbs( 6,7,1)  = -0.29800e-02_wp
            hbs( 6,8,1)  = -0.24840e-01_wp
            hbs( 6,9,1)  = 0.97000e-03_wp
            hbs( 6,10,1) = -0.19290e-01_wp
            hbs( 6,11,1) = -0.90000e-04_wp
            hbs( 6,12,1) = 0.10000e-04_wp
            hbs( 7,1,1)  = -0.12120e-01_wp
            hbs( 7,2,1)  = 0.42330e-01_wp
            hbs( 7,3,1)  = 0.37900e-02_wp
            hbs( 7,4,1)  = -0.32360e-01_wp
            hbs( 7,5,1)  = 0.55000e-02_wp
            hbs( 7,6,1)  = -0.33060e-01_wp
            hbs( 7,7,1)  = -0.89600e-02_wp
            hbs( 7,8,1)  = -0.12200e-01_wp
            hbs( 7,9,1)  = 0.50200e-02_wp
            hbs( 7,10,1) = 0.46020e-01_wp
            hbs( 7,11,1) = -0.96000e-03_wp
            hbs( 7,12,1) = -0.80000e-04_wp
            hbs( 8,1,1)  = 0.67000e-02_wp
            hbs( 8,2,1)  = -0.23490e-01_wp
            hbs( 8,3,1)  = -0.78000e-03_wp
            hbs( 8,4,1)  = 0.30060e-01_wp
            hbs( 8,5,1)  = -0.33500e-02_wp
            hbs( 8,6,1)  = 0.37600e-02_wp
            hbs( 8,7,1)  = 0.46600e-02_wp
            hbs( 8,8,1)  = 0.52800e-02_wp
            hbs( 8,9,1)  = -0.20300e-02_wp
            hbs( 8,10,1) = -0.20610e-01_wp
            hbs( 8,11,1) = -0.19000e-03_wp
            hbs( 8,12,1) = 0.30000e-04_wp
            hbs( 9,1,1)  = 0.54200e-02_wp
            hbs( 9,2,1)  = -0.18840e-01_wp
            hbs( 9,3,1)  = -0.30100e-02_wp
            hbs( 9,4,1)  = 0.22900e-02_wp
            hbs( 9,5,1)  = -0.21500e-02_wp
            hbs( 9,6,1)  = 0.29300e-01_wp
            hbs( 9,7,1)  = 0.43100e-02_wp
            hbs( 9,8,1)  = 0.69200e-02_wp
            hbs( 9,9,1)  = -0.29800e-02_wp
            hbs( 9,10,1) = -0.25420e-01_wp
            hbs( 9,11,1) = 0.11500e-02_wp
            hbs( 9,12,1) = 0.50000e-04_wp
            hbs(10,1,1)  = 0.44050e-01_wp
            hbs(10,2,1)  = 0.43680e-01_wp
            hbs(10,3,1)  = -0.27100e-02_wp
            hbs(10,4,1)  = -0.13440e-01_wp
            hbs(10,5,1)  = -0.64000e-03_wp
            hbs(10,6,1)  = -0.32210e-01_wp
            hbs(10,7,1)  = 0.53500e-02_wp
            hbs(10,8,1)  = -0.31660e-01_wp
            hbs(10,9,1)  = -0.87900e-02_wp
            hbs(10,10,1) = -0.13050e-01_wp
            hbs(10,11,1) = 0.52900e-02_wp
            hbs(10,12,1) = -0.50000e-03_wp
            hbs(11,1,1)  = -0.24480e-01_wp
            hbs(11,2,1)  = -0.19450e-01_wp
            hbs(11,3,1)  = 0.11700e-02_wp
            hbs(11,4,1)  = 0.71100e-02_wp
            hbs(11,5,1)  = 0.10200e-02_wp
            hbs(11,6,1)  = 0.30010e-01_wp
            hbs(11,7,1)  = -0.31900e-02_wp
            hbs(11,8,1)  = 0.20900e-02_wp
            hbs(11,9,1)  = 0.44300e-02_wp
            hbs(11,10,1) = 0.60900e-02_wp
            hbs(11,11,1) = -0.21400e-02_wp
            hbs(11,12,1) = -0.33000e-03_wp
            hbs(12,1,1)  = -0.19570e-01_wp
            hbs(12,2,1)  = -0.24220e-01_wp
            hbs(12,3,1)  = 0.15500e-02_wp
            hbs(12,4,1)  = 0.63300e-02_wp
            hbs(12,5,1)  = -0.38000e-03_wp
            hbs(12,6,1)  = 0.22000e-02_wp
            hbs(12,7,1)  = -0.21500e-02_wp
            hbs(12,8,1)  = 0.29570e-01_wp
            hbs(12,9,1)  = 0.43600e-02_wp
            hbs(12,10,1) = 0.69600e-02_wp
            hbs(12,11,1) = -0.31500e-02_wp
            hbs(12,12,1) = 0.84000e-03_wp
            hbs(13,1,1)  = 0.42020e-01_wp
            hbs(13,2,1)  = -0.11720e-01_wp
            hbs(13,3,1)  = 0.43200e-02_wp
            hbs(13,4,1)  = 0.45370e-01_wp
            hbs(13,5,1)  = -0.30000e-04_wp
            hbs(13,6,1)  = -0.12300e-01_wp
            hbs(13,7,1)  = -0.10000e-02_wp
            hbs(13,8,1)  = -0.32760e-01_wp
            hbs(13,9,1)  = 0.50300e-02_wp
            hbs(13,10,1) = -0.33540e-01_wp
            hbs(13,11,1) = -0.89700e-02_wp
            hbs(13,12,1) = 0.49300e-02_wp
            hbs(14,1,1)  = -0.18670e-01_wp
            hbs(14,2,1)  = 0.51600e-02_wp
            hbs(14,3,1)  = -0.32600e-02_wp
            hbs(14,4,1)  = -0.25130e-01_wp
            hbs(14,5,1)  = 0.30000e-04_wp
            hbs(14,6,1)  = 0.69900e-02_wp
            hbs(14,7,1)  = 0.11900e-02_wp
            hbs(14,8,1)  = 0.29240e-01_wp
            hbs(14,9,1)  = -0.29900e-02_wp
            hbs(14,10,1) = 0.30100e-02_wp
            hbs(14,11,1) = 0.43500e-02_wp
            hbs(14,12,1) = -0.19200e-02_wp
            hbs(15,1,1)  = -0.23350e-01_wp
            hbs(15,2,1)  = 0.65500e-02_wp
            hbs(15,3,1)  = -0.10600e-02_wp
            hbs(15,4,1)  = -0.20240e-01_wp
            hbs(15,5,1)  = 0.00000e+00_wp
            hbs(15,6,1)  = 0.53100e-02_wp
            hbs(15,7,1)  = -0.18000e-03_wp
            hbs(15,8,1)  = 0.35200e-02_wp
            hbs(15,9,1)  = -0.20400e-02_wp
            hbs(15,10,1) = 0.30530e-01_wp
            hbs(15,11,1) = 0.46200e-02_wp
            hbs(15,12,1) = -0.30100e-02_wp
            hbs(16,1,1)  = -0.84300e-02_wp
            hbs(16,2,1)  = -0.29330e-01_wp
            hbs(16,3,1)  = 0.14920e-01_wp
            hbs(16,4,1)  = 0.43340e-01_wp
            hbs(16,5,1)  = -0.70000e-03_wp
            hbs(16,6,1)  = 0.44840e-01_wp
            hbs(16,7,1)  = 0.17000e-03_wp
            hbs(16,8,1)  = -0.12440e-01_wp
            hbs(16,9,1)  = -0.67000e-03_wp
            hbs(16,10,1) = -0.32290e-01_wp
            hbs(16,11,1) = 0.49100e-02_wp
            hbs(16,12,1) = -0.89200e-02_wp
            hbs(17,1,1)  = 0.44500e-02_wp
            hbs(17,2,1)  = 0.25980e-01_wp
            hbs(17,3,1)  = -0.89400e-02_wp
            hbs(17,4,1)  = -0.24010e-01_wp
            hbs(17,5,1)  = 0.10600e-02_wp
            hbs(17,6,1)  = -0.19750e-01_wp
            hbs(17,7,1)  = -0.90000e-04_wp
            hbs(17,8,1)  = 0.59800e-02_wp
            hbs(17,9,1)  = -0.29000e-03_wp
            hbs(17,10,1) = 0.31400e-02_wp
            hbs(17,11,1) = -0.18900e-02_wp
            hbs(17,12,1) = 0.71100e-02_wp
            hbs(18,1,1)  = 0.39800e-02_wp
            hbs(18,2,1)  = 0.33600e-02_wp
            hbs(18,3,1)  = -0.59800e-02_wp
            hbs(18,4,1)  = -0.19330e-01_wp
            hbs(18,5,1)  = -0.36000e-03_wp
            hbs(18,6,1)  = -0.25090e-01_wp
            hbs(18,7,1)  = -0.80000e-04_wp
            hbs(18,8,1)  = 0.64600e-02_wp
            hbs(18,9,1)  = 0.96000e-03_wp
            hbs(18,10,1) = 0.29140e-01_wp
            hbs(18,11,1) = -0.30200e-02_wp
            hbs(18,12,1) = 0.18100e-02_wp
            ! Hessian elements for anchor point 2
            hss( 1,1,2)  = 0.40931e+00_wp
            hss( 2,1,2)  = 0.17750e-01_wp
            hss( 2,2,2)  = 0.39831e+00_wp
            hss( 3,1,2)  = 0.30090e-01_wp
            hss( 3,2,2)  = 0.39670e-01_wp
            hss( 3,3,2)  = 0.41039e+00_wp
            hss( 4,1,2)  = 0.17380e-01_wp
            hss( 4,2,2)  = -0.24100e-02_wp
            hss( 4,3,2)  = -0.47000e-03_wp
            hss( 4,4,2)  = 0.41129e+00_wp
            hss( 5,1,2)  = 0.42800e-02_wp
            hss( 5,2,2)  = -0.46000e-03_wp
            hss( 5,3,2)  = 0.49400e-02_wp
            hss( 5,4,2)  = 0.55000e-02_wp
            hss( 5,5,2)  = 0.35488e+00_wp
            hss( 6,1,2)  = -0.36200e-02_wp
            hss( 6,2,2)  = 0.81920e-01_wp
            hss( 6,3,2)  = -0.25500e-02_wp
            hss( 6,4,2)  = 0.21100e-01_wp
            hss( 6,5,2)  = -0.76000e-03_wp
            hss( 6,6,2)  = 0.41554e+00_wp
            hss( 7,1,2)  = -0.59000e-03_wp
            hss( 7,2,2)  = -0.36000e-03_wp
            hss( 7,3,2)  = 0.73000e-03_wp
            hss( 7,4,2)  = 0.43500e-02_wp
            hss( 7,5,2)  = 0.47000e-03_wp
            hss( 7,6,2)  = 0.46000e-02_wp
            hss( 7,7,2)  = 0.36429e+00_wp
            hss( 8,1,2)  = 0.81370e-01_wp
            hss( 8,2,2)  = -0.43100e-02_wp
            hss( 8,3,2)  = 0.46000e-03_wp
            hss( 8,4,2)  = -0.31800e-02_wp
            hss( 8,5,2)  = -0.28000e-03_wp
            hss( 8,6,2)  = 0.21080e-01_wp
            hss( 8,7,2)  = -0.32000e-03_wp
            hss( 8,8,2)  = 0.40703e+00_wp
            hss( 9,1,2)  = -0.47000e-03_wp
            hss( 9,2,2)  = -0.50000e-03_wp
            hss( 9,3,2)  = -0.34000e-03_wp
            hss( 9,4,2)  = -0.61000e-03_wp
            hss( 9,5,2)  = 0.23000e-03_wp
            hss( 9,6,2)  = 0.44800e-02_wp
            hss( 9,7,2)  = 0.47000e-03_wp
            hss( 9,8,2)  = 0.44100e-02_wp
            hss( 9,9,2)  = 0.36691e+00_wp
            hss(10,1,2)  = -0.11900e-02_wp
            hss(10,2,2)  = 0.18110e-01_wp
            hss(10,3,2)  = -0.13500e-02_wp
            hss(10,4,2)  = 0.85990e-01_wp
            hss(10,5,2)  = -0.79000e-03_wp
            hss(10,6,2)  = -0.28600e-02_wp
            hss(10,7,2)  = -0.55000e-03_wp
            hss(10,8,2)  = 0.21070e-01_wp
            hss(10,9,2)  = -0.50000e-03_wp
            hss(10,10,2) = 0.41886e+00_wp
            hss(11,1,2)  = -0.47000e-03_wp
            hss(11,2,2)  = -0.52000e-03_wp
            hss(11,3,2)  = 0.68000e-03_wp
            hss(11,4,2)  = -0.53000e-03_wp
            hss(11,5,2)  = 0.12000e-03_wp
            hss(11,6,2)  = -0.24000e-03_wp
            hss(11,7,2)  = 0.21000e-03_wp
            hss(11,8,2)  = 0.45400e-02_wp
            hss(11,9,2)  = 0.47000e-03_wp
            hss(11,10,2) = 0.45200e-02_wp
            hss(11,11,2) = 0.36397e+00_wp
            hss(12,1,2)  = -0.53000e-03_wp
            hss(12,2,2)  = 0.41000e-02_wp
            hss(12,3,2)  = -0.53000e-03_wp
            hss(12,4,2)  = -0.83000e-03_wp
            hss(12,5,2)  = 0.34000e-03_wp
            hss(12,6,2)  = -0.14000e-03_wp
            hss(12,7,2)  = 0.10000e-03_wp
            hss(12,8,2)  = -0.69000e-03_wp
            hss(12,9,2)  = 0.14000e-03_wp
            hss(12,10,2) = 0.42400e-02_wp
            hss(12,11,2) = 0.46000e-03_wp
            hss(12,12,2) = 0.36700e+00_wp
            ! Hessian elements for anchor point 2
            hbb( 1,1,2)  = 0.85740e-01_wp
            hbb( 2,1,2)  = -0.47120e-01_wp
            hbb( 2,2,2)  = 0.14272e+00_wp
            hbb( 3,1,2)  = -0.38620e-01_wp
            hbb( 3,2,2)  = -0.95590e-01_wp
            hbb( 3,3,2)  = 0.13421e+00_wp
            hbb( 4,1,2)  = -0.35340e-01_wp
            hbb( 4,2,2)  = 0.24370e-01_wp
            hbb( 4,3,2)  = 0.10970e-01_wp
            hbb( 4,4,2)  = 0.78390e-01_wp
            hbb( 5,1,2)  = 0.19660e-01_wp
            hbb( 5,2,2)  = -0.13050e-01_wp
            hbb( 5,3,2)  = -0.66100e-02_wp
            hbb( 5,4,2)  = -0.38100e-01_wp
            hbb( 5,5,2)  = 0.76190e-01_wp
            hbb( 6,1,2)  = 0.15680e-01_wp
            hbb( 6,2,2)  = -0.11320e-01_wp
            hbb( 6,3,2)  = -0.43700e-02_wp
            hbb( 6,4,2)  = -0.40290e-01_wp
            hbb( 6,5,2)  = -0.38090e-01_wp
            hbb( 6,6,2)  = 0.78390e-01_wp
            hbb( 7,1,2)  = -0.10880e-01_wp
            hbb( 7,2,2)  = -0.19300e-02_wp
            hbb( 7,3,2)  = 0.12810e-01_wp
            hbb( 7,4,2)  = -0.30710e-01_wp
            hbb( 7,5,2)  = 0.11850e-01_wp
            hbb( 7,6,2)  = 0.18860e-01_wp
            hbb( 7,7,2)  = 0.78510e-01_wp
            hbb( 8,1,2)  = -0.26000e-03_wp
            hbb( 8,2,2)  = 0.56200e-02_wp
            hbb( 8,3,2)  = -0.53700e-02_wp
            hbb( 8,4,2)  = 0.18260e-01_wp
            hbb( 8,5,2)  = -0.64700e-02_wp
            hbb( 8,6,2)  = -0.11790e-01_wp
            hbb( 8,7,2)  = -0.39180e-01_wp
            hbb( 8,8,2)  = 0.78940e-01_wp
            hbb( 9,1,2)  = 0.11130e-01_wp
            hbb( 9,2,2)  = -0.36900e-02_wp
            hbb( 9,3,2)  = -0.74400e-02_wp
            hbb( 9,4,2)  = 0.12450e-01_wp
            hbb( 9,5,2)  = -0.53800e-02_wp
            hbb( 9,6,2)  = -0.70800e-02_wp
            hbb( 9,7,2)  = -0.39340e-01_wp
            hbb( 9,8,2)  = -0.39760e-01_wp
            hbb( 9,9,2)  = 0.79100e-01_wp
            hbb(10,1,2)  = 0.73300e-02_wp
            hbb(10,2,2)  = -0.43700e-02_wp
            hbb(10,3,2)  = -0.29600e-02_wp
            hbb(10,4,2)  = -0.10030e-01_wp
            hbb(10,5,2)  = 0.10560e-01_wp
            hbb(10,6,2)  = -0.53000e-03_wp
            hbb(10,7,2)  = -0.33760e-01_wp
            hbb(10,8,2)  = 0.13640e-01_wp
            hbb(10,9,2)  = 0.20120e-01_wp
            hbb(10,10,2) = 0.80240e-01_wp
            hbb(11,1,2)  = -0.36200e-02_wp
            hbb(11,2,2)  = 0.25500e-02_wp
            hbb(11,3,2)  = 0.10700e-02_wp
            hbb(11,4,2)  = -0.61000e-03_wp
            hbb(11,5,2)  = -0.38000e-02_wp
            hbb(11,6,2)  = 0.44100e-02_wp
            hbb(11,7,2)  = 0.19800e-01_wp
            hbb(11,8,2)  = -0.72600e-02_wp
            hbb(11,9,2)  = -0.12540e-01_wp
            hbb(11,10,2) = -0.40120e-01_wp
            hbb(11,11,2) = 0.79050e-01_wp
            hbb(12,1,2)  = -0.37100e-02_wp
            hbb(12,2,2)  = 0.18200e-02_wp
            hbb(12,3,2)  = 0.18900e-02_wp
            hbb(12,4,2)  = 0.10640e-01_wp
            hbb(12,5,2)  = -0.67600e-02_wp
            hbb(12,6,2)  = -0.38800e-02_wp
            hbb(12,7,2)  = 0.13960e-01_wp
            hbb(12,8,2)  = -0.63800e-02_wp
            hbb(12,9,2)  = -0.75800e-02_wp
            hbb(12,10,2) = -0.40110e-01_wp
            hbb(12,11,2) = -0.38920e-01_wp
            hbb(12,12,2) = 0.79040e-01_wp
            hbb(13,1,2)  = -0.10540e-01_wp
            hbb(13,2,2)  = 0.13360e-01_wp
            hbb(13,3,2)  = -0.28100e-02_wp
            hbb(13,4,2)  = 0.52100e-02_wp
            hbb(13,5,2)  = -0.22200e-02_wp
            hbb(13,6,2)  = -0.29900e-02_wp
            hbb(13,7,2)  = -0.87200e-02_wp
            hbb(13,8,2)  = 0.10150e-01_wp
            hbb(13,9,2)  = -0.14300e-02_wp
            hbb(13,10,2) = -0.34040e-01_wp
            hbb(13,11,2) = 0.14040e-01_wp
            hbb(13,12,2) = 0.19990e-01_wp
            hbb(13,13,2) = 0.79490e-01_wp
            hbb(14,1,2)  = 0.10840e-01_wp
            hbb(14,2,2)  = -0.78900e-02_wp
            hbb(14,3,2)  = -0.29500e-02_wp
            hbb(14,4,2)  = -0.26600e-02_wp
            hbb(14,5,2)  = 0.10300e-02_wp
            hbb(14,6,2)  = 0.16200e-02_wp
            hbb(14,7,2)  = -0.14800e-02_wp
            hbb(14,8,2)  = -0.35100e-02_wp
            hbb(14,9,2)  = 0.49900e-02_wp
            hbb(14,10,2) = 0.20150e-01_wp
            hbb(14,11,2) = -0.75900e-02_wp
            hbb(14,12,2) = -0.12560e-01_wp
            hbb(14,13,2) = -0.39740e-01_wp
            hbb(14,14,2) = 0.79580e-01_wp
            hbb(15,1,2)  = -0.30000e-03_wp
            hbb(15,2,2)  = -0.54600e-02_wp
            hbb(15,3,2)  = 0.57600e-02_wp
            hbb(15,4,2)  = -0.25500e-02_wp
            hbb(15,5,2)  = 0.11900e-02_wp
            hbb(15,6,2)  = 0.13700e-02_wp
            hbb(15,7,2)  = 0.10200e-01_wp
            hbb(15,8,2)  = -0.66400e-02_wp
            hbb(15,9,2)  = -0.35600e-02_wp
            hbb(15,10,2) = 0.13890e-01_wp
            hbb(15,11,2) = -0.64600e-02_wp
            hbb(15,12,2) = -0.74300e-02_wp
            hbb(15,13,2) = -0.39750e-01_wp
            hbb(15,14,2) = -0.39830e-01_wp
            hbb(15,15,2) = 0.79580e-01_wp
            hbb(16,1,2)  = -0.36310e-01_wp
            hbb(16,2,2)  = 0.15710e-01_wp
            hbb(16,3,2)  = 0.20610e-01_wp
            hbb(16,4,2)  = -0.75200e-02_wp
            hbb(16,5,2)  = -0.17400e-02_wp
            hbb(16,6,2)  = 0.92600e-02_wp
            hbb(16,7,2)  = 0.55500e-02_wp
            hbb(16,8,2)  = -0.26200e-02_wp
            hbb(16,9,2)  = -0.29400e-02_wp
            hbb(16,10,2) = -0.97400e-02_wp
            hbb(16,11,2) = 0.10510e-01_wp
            hbb(16,12,2) = -0.77000e-03_wp
            hbb(16,13,2) = -0.31400e-01_wp
            hbb(16,14,2) = 0.12890e-01_wp
            hbb(16,15,2) = 0.18510e-01_wp
            hbb(16,16,2) = 0.79420e-01_wp
            hbb(17,1,2)  = 0.22000e-01_wp
            hbb(17,2,2)  = -0.74700e-02_wp
            hbb(17,3,2)  = -0.14540e-01_wp
            hbb(17,4,2)  = -0.24600e-02_wp
            hbb(17,5,2)  = 0.53500e-02_wp
            hbb(17,6,2)  = -0.28900e-02_wp
            hbb(17,7,2)  = -0.26800e-02_wp
            hbb(17,8,2)  = 0.14000e-02_wp
            hbb(17,9,2)  = 0.12800e-02_wp
            hbb(17,10,2) = 0.10490e-01_wp
            hbb(17,11,2) = -0.68100e-02_wp
            hbb(17,12,2) = -0.36800e-02_wp
            hbb(17,13,2) = 0.12230e-01_wp
            hbb(17,14,2) = -0.56800e-02_wp
            hbb(17,15,2) = -0.65500e-02_wp
            hbb(17,16,2) = -0.39570e-01_wp
            hbb(17,17,2) = 0.76420e-01_wp
            hbb(18,1,2)  = 0.14310e-01_wp
            hbb(18,2,2)  = -0.82400e-02_wp
            hbb(18,3,2)  = -0.60700e-02_wp
            hbb(18,4,2)  = 0.99800e-02_wp
            hbb(18,5,2)  = -0.36100e-02_wp
            hbb(18,6,2)  = -0.63700e-02_wp
            hbb(18,7,2)  = -0.28700e-02_wp
            hbb(18,8,2)  = 0.12200e-02_wp
            hbb(18,9,2)  = 0.16500e-02_wp
            hbb(18,10,2) = -0.75000e-03_wp
            hbb(18,11,2) = -0.37000e-02_wp
            hbb(18,12,2) = 0.44500e-02_wp
            hbb(18,13,2) = 0.19170e-01_wp
            hbb(18,14,2) = -0.72100e-02_wp
            hbb(18,15,2) = -0.11960e-01_wp
            hbb(18,16,2) = -0.39850e-01_wp
            hbb(18,17,2) = -0.36850e-01_wp
            hbb(18,18,2) = 0.76700e-01_wp
            ! Hessian elements for anchor point 2
            htt( 1,1,2)  = 0.90500e-02_wp
            htt( 2,1,2)  = 0.38100e-02_wp
            htt( 2,2,2)  = 0.81400e-02_wp
            htt( 3,1,2)  = 0.25000e-03_wp
            htt( 3,2,2)  = -0.43200e-02_wp
            htt( 3,3,2)  = 0.11860e-01_wp
            htt( 4,1,2)  = -0.49900e-02_wp
            htt( 4,2,2)  = 0.10000e-04_wp
            htt( 4,3,2)  = 0.72800e-02_wp
            htt( 4,4,2)  = 0.12280e-01_wp
            htt( 5,1,2)  = -0.60200e-02_wp
            htt( 5,2,2)  = -0.42900e-02_wp
            htt( 5,3,2)  = 0.24900e-02_wp
            htt( 5,4,2)  = 0.42200e-02_wp
            htt( 5,5,2)  = 0.87300e-02_wp
            htt( 6,1,2)  = -0.42000e-02_wp
            htt( 6,2,2)  = -0.28900e-02_wp
            htt( 6,3,2)  = 0.31600e-02_wp
            htt( 6,4,2)  = 0.44800e-02_wp
            htt( 6,5,2)  = 0.34600e-02_wp
            htt( 6,6,2)  = 0.79900e-02_wp
            htt( 7,1,2)  = 0.21800e-02_wp
            htt( 7,2,2)  = 0.32900e-02_wp
            htt( 7,3,2)  = -0.83100e-02_wp
            htt( 7,4,2)  = -0.72000e-02_wp
            htt( 7,5,2)  = 0.81000e-03_wp
            htt( 7,6,2)  = -0.34000e-02_wp
            htt( 7,7,2)  = 0.10570e-01_wp
            htt( 8,1,2)  = 0.39900e-02_wp
            htt( 8,2,2)  = 0.46900e-02_wp
            htt( 8,3,2)  = -0.76400e-02_wp
            htt( 8,4,2)  = -0.69400e-02_wp
            htt( 8,5,2)  = -0.44700e-02_wp
            htt( 8,6,2)  = 0.11300e-02_wp
            htt( 8,7,2)  = 0.63600e-02_wp
            htt( 8,8,2)  = 0.11960e-01_wp
            htt( 9,1,2)  = -0.62900e-02_wp
            htt( 9,2,2)  = -0.85000e-03_wp
            htt( 9,3,2)  = -0.28600e-02_wp
            htt( 9,4,2)  = 0.25800e-02_wp
            htt( 9,5,2)  = 0.40000e-03_wp
            htt( 9,6,2)  = 0.20200e-02_wp
            htt( 9,7,2)  = -0.27900e-02_wp
            htt( 9,8,2)  = -0.11700e-02_wp
            htt( 9,9,2)  = 0.94300e-02_wp
            htt(10,1,2)  = -0.39800e-02_wp
            htt(10,2,2)  = 0.72000e-03_wp
            htt(10,3,2)  = -0.15100e-02_wp
            htt(10,4,2)  = 0.31900e-02_wp
            htt(10,5,2)  = 0.20000e-02_wp
            htt(10,6,2)  = 0.15200e-02_wp
            htt(10,7,2)  = -0.30000e-03_wp
            htt(10,8,2)  = -0.78000e-03_wp
            htt(10,9,2)  = 0.30800e-02_wp
            htt(10,10,2) = 0.83100e-02_wp
            htt(11,1,2)  = -0.99000e-03_wp
            htt(11,2,2)  = -0.52300e-02_wp
            htt(11,3,2)  = 0.17700e-02_wp
            htt(11,4,2)  = -0.24700e-02_wp
            htt(11,5,2)  = -0.13500e-02_wp
            htt(11,6,2)  = 0.69000e-03_wp
            htt(11,7,2)  = -0.39100e-02_wp
            htt(11,8,2)  = -0.18800e-02_wp
            htt(11,9,2)  = 0.39300e-02_wp
            htt(11,10,2) = -0.16800e-02_wp
            htt(11,11,2) = 0.82200e-02_wp
            htt(12,1,2)  = 0.13200e-02_wp
            htt(12,2,2)  = -0.36600e-02_wp
            htt(12,3,2)  = 0.31200e-02_wp
            htt(12,4,2)  = -0.18600e-02_wp
            htt(12,5,2)  = 0.24000e-03_wp
            htt(12,6,2)  = 0.19000e-03_wp
            htt(12,7,2)  = -0.14300e-02_wp
            htt(12,8,2)  = -0.14900e-02_wp
            htt(12,9,2)  = -0.24200e-02_wp
            htt(12,10,2) = 0.35600e-02_wp
            htt(12,11,2) = 0.26100e-02_wp
            htt(12,12,2) = 0.85900e-02_wp
            htt(13,1,2)  = 0.48000e-03_wp
            htt(13,2,2)  = -0.16000e-02_wp
            htt(13,3,2)  = 0.26700e-02_wp
            htt(13,4,2)  = 0.58000e-03_wp
            htt(13,5,2)  = 0.24300e-02_wp
            htt(13,6,2)  = 0.89000e-03_wp
            htt(13,7,2)  = 0.40000e-03_wp
            htt(13,8,2)  = -0.11500e-02_wp
            htt(13,9,2)  = -0.65700e-02_wp
            htt(13,10,2) = -0.18000e-03_wp
            htt(13,11,2) = -0.44600e-02_wp
            htt(13,12,2) = 0.19300e-02_wp
            htt(13,13,2) = 0.95900e-02_wp
            htt(14,1,2)  = 0.19800e-02_wp
            htt(14,2,2)  = 0.50000e-03_wp
            htt(14,3,2)  = 0.14700e-02_wp
            htt(14,4,2)  = -0.10000e-04_wp
            htt(14,5,2)  = 0.10300e-02_wp
            htt(14,6,2)  = -0.11500e-02_wp
            htt(14,7,2)  = 0.15100e-02_wp
            htt(14,8,2)  = -0.68000e-03_wp
            htt(14,9,2)  = -0.44200e-02_wp
            htt(14,10,2) = 0.12200e-02_wp
            htt(14,11,2) = -0.29300e-02_wp
            htt(14,12,2) = 0.27100e-02_wp
            htt(14,13,2) = 0.37900e-02_wp
            htt(14,14,2) = 0.85500e-02_wp
            htt(15,1,2)  = -0.18400e-02_wp
            htt(15,2,2)  = -0.31900e-02_wp
            htt(15,3,2)  = 0.13100e-02_wp
            htt(15,4,2)  = -0.40000e-04_wp
            htt(15,5,2)  = 0.83000e-03_wp
            htt(15,6,2)  = 0.14000e-02_wp
            htt(15,7,2)  = -0.21100e-02_wp
            htt(15,8,2)  = -0.15400e-02_wp
            htt(15,9,2)  = -0.16000e-03_wp
            htt(15,10,2) = -0.54600e-02_wp
            htt(15,11,2) = 0.12000e-02_wp
            htt(15,12,2) = -0.41000e-02_wp
            htt(15,13,2) = 0.31500e-02_wp
            htt(15,14,2) = -0.19100e-02_wp
            htt(15,15,2) = 0.85000e-02_wp
            htt(16,1,2)  = -0.35000e-03_wp
            htt(16,2,2)  = -0.10800e-02_wp
            htt(16,3,2)  = 0.11000e-03_wp
            htt(16,4,2)  = -0.62000e-03_wp
            htt(16,5,2)  = -0.57000e-03_wp
            htt(16,6,2)  = -0.64000e-03_wp
            htt(16,7,2)  = -0.10000e-02_wp
            htt(16,8,2)  = -0.10700e-02_wp
            htt(16,9,2)  = 0.19800e-02_wp
            htt(16,10,2) = -0.40600e-02_wp
            htt(16,11,2) = 0.27200e-02_wp
            htt(16,12,2) = -0.33200e-02_wp
            htt(16,13,2) = -0.26500e-02_wp
            htt(16,14,2) = 0.28600e-02_wp
            htt(16,15,2) = 0.34400e-02_wp
            htt(16,16,2) = 0.89500e-02_wp
            htt(17,1,2)  = 0.25600e-02_wp
            htt(17,2,2)  = 0.11000e-02_wp
            htt(17,3,2)  = 0.12000e-03_wp
            htt(17,4,2)  = -0.13400e-02_wp
            htt(17,5,2)  = 0.36000e-03_wp
            htt(17,6,2)  = -0.16300e-02_wp
            htt(17,7,2)  = 0.26400e-02_wp
            htt(17,8,2)  = 0.64000e-03_wp
            htt(17,9,2)  = 0.63000e-03_wp
            htt(17,10,2) = -0.18200e-02_wp
            htt(17,11,2) = 0.21000e-02_wp
            htt(17,12,2) = -0.35000e-03_wp
            htt(17,13,2) = -0.66700e-02_wp
            htt(17,14,2) = -0.74000e-03_wp
            htt(17,15,2) = -0.42000e-02_wp
            htt(17,16,2) = 0.17300e-02_wp
            htt(17,17,2) = 0.96700e-02_wp
            htt(18,1,2)  = 0.74000e-03_wp
            htt(18,2,2)  = 0.13700e-02_wp
            htt(18,3,2)  = -0.19900e-02_wp
            htt(18,4,2)  = -0.13600e-02_wp
            htt(18,5,2)  = -0.16000e-02_wp
            htt(18,6,2)  = -0.30000e-02_wp
            htt(18,7,2)  = 0.94000e-03_wp
            htt(18,8,2)  = -0.46000e-03_wp
            htt(18,9,2)  = 0.21100e-02_wp
            htt(18,10,2) = 0.44000e-03_wp
            htt(18,11,2) = 0.14700e-02_wp
            htt(18,12,2) = -0.19000e-03_wp
            htt(18,13,2) = -0.40300e-02_wp
            htt(18,14,2) = 0.10300e-02_wp
            htt(18,15,2) = -0.23600e-02_wp
            htt(18,16,2) = 0.27100e-02_wp
            htt(18,17,2) = 0.31700e-02_wp
            htt(18,18,2) = 0.84700e-02_wp
            htt(19,1,2)  = 0.10700e-02_wp
            htt(19,2,2)  = -0.10000e-02_wp
            htt(19,3,2)  = 0.13200e-02_wp
            htt(19,4,2)  = -0.75000e-03_wp
            htt(19,5,2)  = 0.17600e-02_wp
            htt(19,6,2)  = 0.41000e-03_wp
            htt(19,7,2)  = 0.15300e-02_wp
            htt(19,8,2)  = 0.17000e-03_wp
            htt(19,9,2)  = -0.15200e-02_wp
            htt(19,10,2) = -0.32300e-02_wp
            htt(19,11,2) = 0.58000e-03_wp
            htt(19,12,2) = -0.11300e-02_wp
            htt(19,13,2) = -0.86000e-03_wp
            htt(19,14,2) = -0.55100e-02_wp
            htt(19,15,2) = 0.86000e-03_wp
            htt(19,16,2) = -0.37900e-02_wp
            htt(19,17,2) = 0.37400e-02_wp
            htt(19,18,2) = -0.18900e-02_wp
            htt(19,19,2) = 0.84000e-02_wp
            htt(20,1,2)  = -0.75000e-03_wp
            htt(20,2,2)  = -0.74000e-03_wp
            htt(20,3,2)  = -0.79000e-03_wp
            htt(20,4,2)  = -0.78000e-03_wp
            htt(20,5,2)  = -0.20000e-03_wp
            htt(20,6,2)  = -0.96000e-03_wp
            htt(20,7,2)  = -0.17000e-03_wp
            htt(20,8,2)  = -0.93000e-03_wp
            htt(20,9,2)  = -0.40000e-04_wp
            htt(20,10,2) = -0.97000e-03_wp
            htt(20,11,2) = -0.50000e-04_wp
            htt(20,12,2) = -0.98000e-03_wp
            htt(20,13,2) = 0.17700e-02_wp
            htt(20,14,2) = -0.37400e-02_wp
            htt(20,15,2) = 0.27000e-02_wp
            htt(20,16,2) = -0.28100e-02_wp
            htt(20,17,2) = -0.27500e-02_wp
            htt(20,18,2) = 0.34100e-02_wp
            htt(20,19,2) = 0.27600e-02_wp
            htt(20,20,2) = 0.89200e-02_wp
            htt(21,1,2)  = 0.21000e-03_wp
            htt(21,2,2)  = 0.18300e-02_wp
            htt(21,3,2)  = -0.26900e-02_wp
            htt(21,4,2)  = -0.10600e-02_wp
            htt(21,5,2)  = -0.59300e-02_wp
            htt(21,6,2)  = -0.53000e-03_wp
            htt(21,7,2)  = -0.32400e-02_wp
            htt(21,8,2)  = 0.21600e-02_wp
            htt(21,9,2)  = 0.24300e-02_wp
            htt(21,10,2) = 0.91000e-03_wp
            htt(21,11,2) = 0.79000e-03_wp
            htt(21,12,2) = -0.73000e-03_wp
            htt(21,13,2) = 0.69000e-03_wp
            htt(21,14,2) = -0.16600e-02_wp
            htt(21,15,2) = 0.22200e-02_wp
            htt(21,16,2) = -0.12000e-03_wp
            htt(21,17,2) = -0.65400e-02_wp
            htt(21,18,2) = -0.37000e-03_wp
            htt(21,19,2) = -0.42000e-02_wp
            htt(21,20,2) = 0.19700e-02_wp
            htt(21,21,2) = 0.91600e-02_wp
            htt(22,1,2)  = -0.16600e-02_wp
            htt(22,2,2)  = 0.40000e-03_wp
            htt(22,3,2)  = -0.33700e-02_wp
            htt(22,4,2)  = -0.13200e-02_wp
            htt(22,5,2)  = -0.52000e-03_wp
            htt(22,6,2)  = -0.51900e-02_wp
            htt(22,7,2)  = 0.10800e-02_wp
            htt(22,8,2)  = -0.35900e-02_wp
            htt(22,9,2)  = 0.77000e-03_wp
            htt(22,10,2) = 0.14000e-02_wp
            htt(22,11,2) = -0.13100e-02_wp
            htt(22,12,2) = -0.68000e-03_wp
            htt(22,13,2) = 0.22700e-02_wp
            htt(22,14,2) = 0.59000e-03_wp
            htt(22,15,2) = 0.16300e-02_wp
            htt(22,16,2) = -0.50000e-04_wp
            htt(22,17,2) = -0.44900e-02_wp
            htt(22,18,2) = 0.10700e-02_wp
            htt(22,19,2) = -0.28100e-02_wp
            htt(22,20,2) = 0.27500e-02_wp
            htt(22,21,2) = 0.36200e-02_wp
            htt(22,22,2) = 0.84200e-02_wp
            htt(23,1,2)  = 0.20100e-02_wp
            htt(23,2,2)  = 0.15700e-02_wp
            htt(23,3,2)  = -0.60000e-03_wp
            htt(23,4,2)  = -0.10400e-02_wp
            htt(23,5,2)  = -0.39900e-02_wp
            htt(23,6,2)  = 0.82000e-03_wp
            htt(23,7,2)  = -0.15600e-02_wp
            htt(23,8,2)  = 0.32500e-02_wp
            htt(23,9,2)  = 0.97000e-03_wp
            htt(23,10,2) = -0.13300e-02_wp
            htt(23,11,2) = 0.14100e-02_wp
            htt(23,12,2) = -0.89000e-03_wp
            htt(23,13,2) = -0.19300e-02_wp
            htt(23,14,2) = -0.34100e-02_wp
            htt(23,15,2) = 0.39000e-03_wp
            htt(23,16,2) = -0.11000e-02_wp
            htt(23,17,2) = -0.90000e-04_wp
            htt(23,18,2) = -0.56300e-02_wp
            htt(23,19,2) = 0.13900e-02_wp
            htt(23,20,2) = -0.41400e-02_wp
            htt(23,21,2) = 0.30400e-02_wp
            htt(23,22,2) = -0.18900e-02_wp
            htt(23,23,2) = 0.85300e-02_wp
            htt(24,1,2)  = 0.15000e-03_wp
            htt(24,2,2)  = 0.14000e-03_wp
            htt(24,3,2)  = -0.12900e-02_wp
            htt(24,4,2)  = -0.13000e-02_wp
            htt(24,5,2)  = 0.14300e-02_wp
            htt(24,6,2)  = -0.38400e-02_wp
            htt(24,7,2)  = 0.27700e-02_wp
            htt(24,8,2)  = -0.25000e-02_wp
            htt(24,9,2)  = -0.70000e-03_wp
            htt(24,10,2) = -0.84000e-03_wp
            htt(24,11,2) = -0.68000e-03_wp
            htt(24,12,2) = -0.83000e-03_wp
            htt(24,13,2) = -0.34000e-03_wp
            htt(24,14,2) = -0.11700e-02_wp
            htt(24,15,2) = -0.19000e-03_wp
            htt(24,16,2) = -0.10200e-02_wp
            htt(24,17,2) = 0.19500e-02_wp
            htt(24,18,2) = -0.41900e-02_wp
            htt(24,19,2) = 0.27800e-02_wp
            htt(24,20,2) = -0.33600e-02_wp
            htt(24,21,2) = -0.25000e-02_wp
            htt(24,22,2) = 0.29000e-02_wp
            htt(24,23,2) = 0.36000e-02_wp
            htt(24,24,2) = 0.90000e-02_wp
            ! Hessian elements for anchor point 2
            hbs( 1,1,2)  = -0.37150e-01_wp
            hbs( 1,2,2)  = -0.35080e-01_wp
            hbs( 1,3,2)  = -0.37820e-01_wp
            hbs( 1,4,2)  = -0.11490e-01_wp
            hbs( 1,5,2)  = 0.62300e-02_wp
            hbs( 1,6,2)  = 0.46070e-01_wp
            hbs( 1,7,2)  = -0.43000e-03_wp
            hbs( 1,8,2)  = 0.44740e-01_wp
            hbs( 1,9,2)  = 0.80000e-04_wp
            hbs( 1,10,2) = -0.10140e-01_wp
            hbs( 1,11,2) = -0.44000e-03_wp
            hbs( 1,12,2) = 0.53800e-02_wp
            hbs( 2,1,2)  = 0.45070e-01_wp
            hbs( 2,2,2)  = -0.50600e-02_wp
            hbs( 2,3,2)  = 0.16200e-01_wp
            hbs( 2,4,2)  = -0.18400e-02_wp
            hbs( 2,5,2)  = -0.27200e-02_wp
            hbs( 2,6,2)  = -0.20110e-01_wp
            hbs( 2,7,2)  = 0.11400e-02_wp
            hbs( 2,8,2)  = -0.26860e-01_wp
            hbs( 2,9,2)  = -0.17000e-03_wp
            hbs( 2,10,2) = 0.11580e-01_wp
            hbs( 2,11,2) = -0.64000e-03_wp
            hbs( 2,12,2) = -0.16600e-02_wp
            hbs( 3,1,2)  = -0.79200e-02_wp
            hbs( 3,2,2)  = 0.40140e-01_wp
            hbs( 3,3,2)  = 0.21620e-01_wp
            hbs( 3,4,2)  = 0.13330e-01_wp
            hbs( 3,5,2)  = -0.35100e-02_wp
            hbs( 3,6,2)  = -0.25960e-01_wp
            hbs( 3,7,2)  = -0.71000e-03_wp
            hbs( 3,8,2)  = -0.17880e-01_wp
            hbs( 3,9,2)  = 0.10000e-03_wp
            hbs( 3,10,2) = -0.14400e-02_wp
            hbs( 3,11,2) = 0.10800e-02_wp
            hbs( 3,12,2) = -0.37200e-02_wp
            hbs( 4,1,2)  = -0.29410e-01_wp
            hbs( 4,2,2)  = -0.87500e-02_wp
            hbs( 4,3,2)  = 0.18410e-01_wp
            hbs( 4,4,2)  = -0.31300e-01_wp
            hbs( 4,5,2)  = -0.10370e-01_wp
            hbs( 4,6,2)  = -0.13280e-01_wp
            hbs( 4,7,2)  = 0.49300e-02_wp
            hbs( 4,8,2)  = 0.44510e-01_wp
            hbs( 4,9,2)  = -0.65000e-03_wp
            hbs( 4,10,2) = 0.43310e-01_wp
            hbs( 4,11,2) = 0.18000e-03_wp
            hbs( 4,12,2) = -0.83000e-03_wp
            hbs( 5,1,2)  = 0.28150e-01_wp
            hbs( 5,2,2)  = 0.54900e-02_wp
            hbs( 5,3,2)  = -0.13780e-01_wp
            hbs( 5,4,2)  = 0.21700e-02_wp
            hbs( 5,5,2)  = 0.70200e-02_wp
            hbs( 5,6,2)  = 0.60300e-02_wp
            hbs( 5,7,2)  = -0.19500e-02_wp
            hbs( 5,8,2)  = -0.19880e-01_wp
            hbs( 5,9,2)  = -0.33000e-03_wp
            hbs( 5,10,2) = -0.24160e-01_wp
            hbs( 5,11,2) = -0.12000e-03_wp
            hbs( 5,12,2) = 0.77000e-03_wp
            hbs( 6,1,2)  = 0.12600e-02_wp
            hbs( 6,2,2)  = 0.32500e-02_wp
            hbs( 6,3,2)  = -0.46300e-02_wp
            hbs( 6,4,2)  = 0.29130e-01_wp
            hbs( 6,5,2)  = 0.33500e-02_wp
            hbs( 6,6,2)  = 0.72400e-02_wp
            hbs( 6,7,2)  = -0.29800e-02_wp
            hbs( 6,8,2)  = -0.24630e-01_wp
            hbs( 6,9,2)  = 0.99000e-03_wp
            hbs( 6,10,2) = -0.19140e-01_wp
            hbs( 6,11,2) = -0.70000e-04_wp
            hbs( 6,12,2) = 0.60000e-04_wp
            hbs( 7,1,2)  = -0.11820e-01_wp
            hbs( 7,2,2)  = 0.41350e-01_wp
            hbs( 7,3,2)  = 0.31600e-02_wp
            hbs( 7,4,2)  = -0.32360e-01_wp
            hbs( 7,5,2)  = 0.57200e-02_wp
            hbs( 7,6,2)  = -0.33140e-01_wp
            hbs( 7,7,2)  = -0.90100e-02_wp
            hbs( 7,8,2)  = -0.12310e-01_wp
            hbs( 7,9,2)  = 0.50300e-02_wp
            hbs( 7,10,2) = 0.46100e-01_wp
            hbs( 7,11,2) = -0.96000e-03_wp
            hbs( 7,12,2) = -0.70000e-04_wp
            hbs( 8,1,2)  = 0.65000e-02_wp
            hbs( 8,2,2)  = -0.22990e-01_wp
            hbs( 8,3,2)  = -0.46000e-03_wp
            hbs( 8,4,2)  = 0.29900e-01_wp
            hbs( 8,5,2)  = -0.34100e-02_wp
            hbs( 8,6,2)  = 0.38600e-02_wp
            hbs( 8,7,2)  = 0.47600e-02_wp
            hbs( 8,8,2)  = 0.53400e-02_wp
            hbs( 8,9,2)  = -0.20400e-02_wp
            hbs( 8,10,2) = -0.20620e-01_wp
            hbs( 8,11,2) = -0.19000e-03_wp
            hbs( 8,12,2) = 0.30000e-04_wp
            hbs( 9,1,2)  = 0.53200e-02_wp
            hbs( 9,2,2)  = -0.18360e-01_wp
            hbs( 9,3,2)  = -0.27000e-02_wp
            hbs( 9,4,2)  = 0.24600e-02_wp
            hbs( 9,5,2)  = -0.23100e-02_wp
            hbs( 9,6,2)  = 0.29280e-01_wp
            hbs( 9,7,2)  = 0.42500e-02_wp
            hbs( 9,8,2)  = 0.69700e-02_wp
            hbs( 9,9,2)  = -0.29900e-02_wp
            hbs( 9,10,2) = -0.25480e-01_wp
            hbs( 9,11,2) = 0.11500e-02_wp
            hbs( 9,12,2) = 0.40000e-04_wp
            hbs(10,1,2)  = 0.43830e-01_wp
            hbs(10,2,2)  = 0.43000e-01_wp
            hbs(10,3,2)  = -0.29200e-02_wp
            hbs(10,4,2)  = -0.13560e-01_wp
            hbs(10,5,2)  = -0.60000e-03_wp
            hbs(10,6,2)  = -0.32140e-01_wp
            hbs(10,7,2)  = 0.53900e-02_wp
            hbs(10,8,2)  = -0.31440e-01_wp
            hbs(10,9,2)  = -0.88000e-02_wp
            hbs(10,10,2) = -0.13100e-01_wp
            hbs(10,11,2) = 0.53000e-02_wp
            hbs(10,12,2) = -0.54000e-03_wp
            hbs(11,1,2)  = -0.24340e-01_wp
            hbs(11,2,2)  = -0.19110e-01_wp
            hbs(11,3,2)  = 0.11900e-02_wp
            hbs(11,4,2)  = 0.71000e-02_wp
            hbs(11,5,2)  = 0.10500e-02_wp
            hbs(11,6,2)  = 0.30020e-01_wp
            hbs(11,7,2)  = -0.32200e-02_wp
            hbs(11,8,2)  = 0.19900e-02_wp
            hbs(11,9,2)  = 0.44100e-02_wp
            hbs(11,10,2) = 0.61200e-02_wp
            hbs(11,11,2) = -0.21600e-02_wp
            hbs(11,12,2) = -0.33000e-03_wp
            hbs(12,1,2)  = -0.19490e-01_wp
            hbs(12,2,2)  = -0.23890e-01_wp
            hbs(12,3,2)  = 0.17300e-02_wp
            hbs(12,4,2)  = 0.64600e-02_wp
            hbs(12,5,2)  = -0.45000e-03_wp
            hbs(12,6,2)  = 0.21300e-02_wp
            hbs(12,7,2)  = -0.21700e-02_wp
            hbs(12,8,2)  = 0.29450e-01_wp
            hbs(12,9,2)  = 0.43900e-02_wp
            hbs(12,10,2) = 0.69800e-02_wp
            hbs(12,11,2) = -0.31500e-02_wp
            hbs(12,12,2) = 0.87000e-03_wp
            hbs(13,1,2)  = 0.42040e-01_wp
            hbs(13,2,2)  = -0.12440e-01_wp
            hbs(13,3,2)  = 0.58600e-02_wp
            hbs(13,4,2)  = 0.45620e-01_wp
            hbs(13,5,2)  = -0.90000e-04_wp
            hbs(13,6,2)  = -0.12560e-01_wp
            hbs(13,7,2)  = -0.10300e-02_wp
            hbs(13,8,2)  = -0.32810e-01_wp
            hbs(13,9,2)  = 0.50400e-02_wp
            hbs(13,10,2) = -0.33390e-01_wp
            hbs(13,11,2) = -0.90100e-02_wp
            hbs(13,12,2) = 0.50300e-02_wp
            hbs(14,1,2)  = -0.18640e-01_wp
            hbs(14,2,2)  = 0.53900e-02_wp
            hbs(14,3,2)  = -0.35400e-02_wp
            hbs(14,4,2)  = -0.25230e-01_wp
            hbs(14,5,2)  = 0.70000e-04_wp
            hbs(14,6,2)  = 0.71300e-02_wp
            hbs(14,7,2)  = 0.12100e-02_wp
            hbs(14,8,2)  = 0.29250e-01_wp
            hbs(14,9,2)  = -0.29800e-02_wp
            hbs(14,10,2) = 0.28500e-02_wp
            hbs(14,11,2) = 0.44100e-02_wp
            hbs(14,12,2) = -0.19500e-02_wp
            hbs(15,1,2)  = -0.23400e-01_wp
            hbs(15,2,2)  = 0.70500e-02_wp
            hbs(15,3,2)  = -0.23200e-02_wp
            hbs(15,4,2)  = -0.20390e-01_wp
            hbs(15,5,2)  = 0.20000e-04_wp
            hbs(15,6,2)  = 0.54300e-02_wp
            hbs(15,7,2)  = -0.17000e-03_wp
            hbs(15,8,2)  = 0.35600e-02_wp
            hbs(15,9,2)  = -0.20600e-02_wp
            hbs(15,10,2) = 0.30550e-01_wp
            hbs(15,11,2) = 0.46000e-02_wp
            hbs(15,12,2) = -0.30800e-02_wp
            hbs(16,1,2)  = -0.75000e-02_wp
            hbs(16,2,2)  = -0.28090e-01_wp
            hbs(16,3,2)  = 0.13300e-01_wp
            hbs(16,4,2)  = 0.43100e-01_wp
            hbs(16,5,2)  = -0.88000e-03_wp
            hbs(16,6,2)  = 0.45050e-01_wp
            hbs(16,7,2)  = 0.16000e-03_wp
            hbs(16,8,2)  = -0.12690e-01_wp
            hbs(16,9,2)  = -0.70000e-03_wp
            hbs(16,10,2) = -0.32780e-01_wp
            hbs(16,11,2) = 0.49200e-02_wp
            hbs(16,12,2) = -0.89700e-02_wp
            hbs(17,1,2)  = 0.38600e-02_wp
            hbs(17,2,2)  = 0.25790e-01_wp
            hbs(17,3,2)  = -0.82600e-02_wp
            hbs(17,4,2)  = -0.23900e-01_wp
            hbs(17,5,2)  = 0.12400e-02_wp
            hbs(17,6,2)  = -0.19910e-01_wp
            hbs(17,7,2)  = -0.90000e-04_wp
            hbs(17,8,2)  = 0.61400e-02_wp
            hbs(17,9,2)  = -0.27000e-03_wp
            hbs(17,10,2) = 0.32800e-02_wp
            hbs(17,11,2) = -0.18800e-02_wp
            hbs(17,12,2) = 0.70400e-02_wp
            hbs(18,1,2)  = 0.36300e-02_wp
            hbs(18,2,2)  = 0.22900e-02_wp
            hbs(18,3,2)  = -0.50400e-02_wp
            hbs(18,4,2)  = -0.19200e-01_wp
            hbs(18,5,2)  = -0.37000e-03_wp
            hbs(18,6,2)  = -0.25140e-01_wp
            hbs(18,7,2)  = -0.80000e-04_wp
            hbs(18,8,2)  = 0.65600e-02_wp
            hbs(18,9,2)  = 0.97000e-03_wp
            hbs(18,10,2) = 0.29500e-01_wp
            hbs(18,11,2) = -0.30300e-02_wp
            hbs(18,12,2) = 0.19300e-02_wp
            ! Hessian elements for anchor point 3
            hss( 1,1,3)  = 0.39742e+00_wp
            hss( 2,1,3)  = 0.15810e-01_wp
            hss( 2,2,3)  = 0.37782e+00_wp
            hss( 3,1,3)  = 0.44160e-01_wp
            hss( 3,2,3)  = 0.54830e-01_wp
            hss( 3,3,3)  = 0.38872e+00_wp
            hss( 4,1,3)  = 0.19410e-01_wp
            hss( 4,2,3)  = -0.16000e-02_wp
            hss( 4,3,3)  = -0.34300e-02_wp
            hss( 4,4,3)  = 0.41358e+00_wp
            hss( 5,1,3)  = 0.23700e-02_wp
            hss( 5,2,3)  = -0.70000e-03_wp
            hss( 5,3,3)  = 0.14900e-02_wp
            hss( 5,4,3)  = 0.48900e-02_wp
            hss( 5,5,3)  = 0.36213e+00_wp
            hss( 6,1,3)  = -0.43800e-02_wp
            hss( 6,2,3)  = 0.80920e-01_wp
            hss( 6,3,3)  = -0.20700e-02_wp
            hss( 6,4,3)  = 0.21020e-01_wp
            hss( 6,5,3)  = -0.70000e-03_wp
            hss( 6,6,3)  = 0.41580e+00_wp
            hss( 7,1,3)  = -0.56000e-03_wp
            hss( 7,2,3)  = -0.41000e-03_wp
            hss( 7,3,3)  = 0.64000e-03_wp
            hss( 7,4,3)  = 0.44200e-02_wp
            hss( 7,5,3)  = 0.38000e-03_wp
            hss( 7,6,3)  = 0.45800e-02_wp
            hss( 7,7,3)  = 0.36425e+00_wp
            hss( 8,1,3)  = 0.80680e-01_wp
            hss( 8,2,3)  = -0.46300e-02_wp
            hss( 8,3,3)  = -0.10400e-02_wp
            hss( 8,4,3)  = -0.39500e-02_wp
            hss( 8,5,3)  = 0.11000e-03_wp
            hss( 8,6,3)  = 0.21280e-01_wp
            hss( 8,7,3)  = -0.32000e-03_wp
            hss( 8,8,3)  = 0.41002e+00_wp
            hss( 9,1,3)  = -0.40000e-03_wp
            hss( 9,2,3)  = -0.40000e-03_wp
            hss( 9,3,3)  = -0.55000e-03_wp
            hss( 9,4,3)  = -0.66000e-03_wp
            hss( 9,5,3)  = 0.17000e-03_wp
            hss( 9,6,3)  = 0.44100e-02_wp
            hss( 9,7,3)  = 0.44000e-03_wp
            hss( 9,8,3)  = 0.43500e-02_wp
            hss( 9,9,3)  = 0.36777e+00_wp
            hss(10,1,3)  = -0.11100e-02_wp
            hss(10,2,3)  = 0.19640e-01_wp
            hss(10,3,3)  = -0.29200e-02_wp
            hss(10,4,3)  = 0.86160e-01_wp
            hss(10,5,3)  = -0.75000e-03_wp
            hss(10,6,3)  = -0.35500e-02_wp
            hss(10,7,3)  = -0.56000e-03_wp
            hss(10,8,3)  = 0.21260e-01_wp
            hss(10,9,3)  = -0.57000e-03_wp
            hss(10,10,3) = 0.41917e+00_wp
            hss(11,1,3)  = -0.46000e-03_wp
            hss(11,2,3)  = -0.58000e-03_wp
            hss(11,3,3)  = 0.53000e-03_wp
            hss(11,4,3)  = -0.55000e-03_wp
            hss(11,5,3)  = 0.90000e-04_wp
            hss(11,6,3)  = -0.25000e-03_wp
            hss(11,7,3)  = 0.21000e-03_wp
            hss(11,8,3)  = 0.45200e-02_wp
            hss(11,9,3)  = 0.43000e-03_wp
            hss(11,10,3) = 0.45100e-02_wp
            hss(11,11,3) = 0.36393e+00_wp
            hss(12,1,3)  = 0.20000e-04_wp
            hss(12,2,3)  = 0.49500e-02_wp
            hss(12,3,3)  = -0.25500e-02_wp
            hss(12,4,3)  = -0.86000e-03_wp
            hss(12,5,3)  = 0.25000e-03_wp
            hss(12,6,3)  = -0.25000e-03_wp
            hss(12,7,3)  = 0.90000e-04_wp
            hss(12,8,3)  = -0.82000e-03_wp
            hss(12,9,3)  = 0.13000e-03_wp
            hss(12,10,3) = 0.39900e-02_wp
            hss(12,11,3) = 0.47000e-03_wp
            hss(12,12,3) = 0.36874e+00_wp
            ! Hessian elements for anchor point 3
            hbb( 1,1,3)  = 0.82270e-01_wp
            hbb( 2,1,3)  = -0.48250e-01_wp
            hbb( 2,2,3)  = 0.13087e+00_wp
            hbb( 3,1,3)  = -0.34020e-01_wp
            hbb( 3,2,3)  = -0.82620e-01_wp
            hbb( 3,3,3)  = 0.11664e+00_wp
            hbb( 4,1,3)  = -0.33450e-01_wp
            hbb( 4,2,3)  = 0.20400e-01_wp
            hbb( 4,3,3)  = 0.13040e-01_wp
            hbb( 4,4,3)  = 0.77460e-01_wp
            hbb( 5,1,3)  = 0.17380e-01_wp
            hbb( 5,2,3)  = -0.84000e-02_wp
            hbb( 5,3,3)  = -0.89800e-02_wp
            hbb( 5,4,3)  = -0.37880e-01_wp
            hbb( 5,5,3)  = 0.75680e-01_wp
            hbb( 6,1,3)  = 0.16060e-01_wp
            hbb( 6,2,3)  = -0.12000e-01_wp
            hbb( 6,3,3)  = -0.40600e-02_wp
            hbb( 6,4,3)  = -0.39580e-01_wp
            hbb( 6,5,3)  = -0.37810e-01_wp
            hbb( 6,6,3)  = 0.77390e-01_wp
            hbb( 7,1,3)  = -0.10450e-01_wp
            hbb( 7,2,3)  = 0.14900e-02_wp
            hbb( 7,3,3)  = 0.89500e-02_wp
            hbb( 7,4,3)  = -0.31240e-01_wp
            hbb( 7,5,3)  = 0.12780e-01_wp
            hbb( 7,6,3)  = 0.18460e-01_wp
            hbb( 7,7,3)  = 0.78900e-01_wp
            hbb( 8,1,3)  = -0.47000e-03_wp
            hbb( 8,2,3)  = 0.30400e-02_wp
            hbb( 8,3,3)  = -0.25700e-02_wp
            hbb( 8,4,3)  = 0.18370e-01_wp
            hbb( 8,5,3)  = -0.69100e-02_wp
            hbb( 8,6,3)  = -0.11460e-01_wp
            hbb( 8,7,3)  = -0.39290e-01_wp
            hbb( 8,8,3)  = 0.78920e-01_wp
            hbb( 9,1,3)  = 0.10920e-01_wp
            hbb( 9,2,3)  = -0.45400e-02_wp
            hbb( 9,3,3)  = -0.63800e-02_wp
            hbb( 9,4,3)  = 0.12870e-01_wp
            hbb( 9,5,3)  = -0.58800e-02_wp
            hbb( 9,6,3)  = -0.70000e-02_wp
            hbb( 9,7,3)  = -0.39610e-01_wp
            hbb( 9,8,3)  = -0.39630e-01_wp
            hbb( 9,9,3)  = 0.79240e-01_wp
            hbb(10,1,3)  = 0.65500e-02_wp
            hbb(10,2,3)  = -0.46300e-02_wp
            hbb(10,3,3)  = -0.19200e-02_wp
            hbb(10,4,3)  = -0.92300e-02_wp
            hbb(10,5,3)  = 0.10020e-01_wp
            hbb(10,6,3)  = -0.79000e-03_wp
            hbb(10,7,3)  = -0.34090e-01_wp
            hbb(10,8,3)  = 0.13830e-01_wp
            hbb(10,9,3)  = 0.20260e-01_wp
            hbb(10,10,3) = 0.80030e-01_wp
            hbb(11,1,3)  = -0.33400e-02_wp
            hbb(11,2,3)  = 0.27100e-02_wp
            hbb(11,3,3)  = 0.63000e-03_wp
            hbb(11,4,3)  = -0.98000e-03_wp
            hbb(11,5,3)  = -0.35000e-02_wp
            hbb(11,6,3)  = 0.44800e-02_wp
            hbb(11,7,3)  = 0.19980e-01_wp
            hbb(11,8,3)  = -0.73600e-02_wp
            hbb(11,9,3)  = -0.12620e-01_wp
            hbb(11,10,3) = -0.40030e-01_wp
            hbb(11,11,3) = 0.78930e-01_wp
            hbb(12,1,3)  = -0.32100e-02_wp
            hbb(12,2,3)  = 0.19200e-02_wp
            hbb(12,3,3)  = 0.13000e-02_wp
            hbb(12,4,3)  = 0.10210e-01_wp
            hbb(12,5,3)  = -0.65300e-02_wp
            hbb(12,6,3)  = -0.36800e-02_wp
            hbb(12,7,3)  = 0.14110e-01_wp
            hbb(12,8,3)  = -0.64600e-02_wp
            hbb(12,9,3)  = -0.76400e-02_wp
            hbb(12,10,3) = -0.40000e-01_wp
            hbb(12,11,3) = -0.38900e-01_wp
            hbb(12,12,3) = 0.78900e-01_wp
            hbb(13,1,3)  = -0.96900e-02_wp
            hbb(13,2,3)  = 0.10920e-01_wp
            hbb(13,3,3)  = -0.12400e-02_wp
            hbb(13,4,3)  = 0.49700e-02_wp
            hbb(13,5,3)  = -0.18200e-02_wp
            hbb(13,6,3)  = -0.31500e-02_wp
            hbb(13,7,3)  = -0.88100e-02_wp
            hbb(13,8,3)  = 0.10200e-01_wp
            hbb(13,9,3)  = -0.13900e-02_wp
            hbb(13,10,3) = -0.34270e-01_wp
            hbb(13,11,3) = 0.14180e-01_wp
            hbb(13,12,3) = 0.20100e-01_wp
            hbb(13,13,3) = 0.79280e-01_wp
            hbb(14,1,3)  = 0.10380e-01_wp
            hbb(14,2,3)  = -0.74600e-02_wp
            hbb(14,3,3)  = -0.29200e-02_wp
            hbb(14,4,3)  = -0.26400e-02_wp
            hbb(14,5,3)  = 0.92000e-03_wp
            hbb(14,6,3)  = 0.17200e-02_wp
            hbb(14,7,3)  = -0.14500e-02_wp
            hbb(14,8,3)  = -0.35600e-02_wp
            hbb(14,9,3)  = 0.50100e-02_wp
            hbb(14,10,3) = 0.20210e-01_wp
            hbb(14,11,3) = -0.76400e-02_wp
            hbb(14,12,3) = -0.12580e-01_wp
            hbb(14,13,3) = -0.39700e-01_wp
            hbb(14,14,3) = 0.79390e-01_wp
            hbb(15,1,3)  = -0.69000e-03_wp
            hbb(15,2,3)  = -0.34600e-02_wp
            hbb(15,3,3)  = 0.41500e-02_wp
            hbb(15,4,3)  = -0.23300e-02_wp
            hbb(15,5,3)  = 0.91000e-03_wp
            hbb(15,6,3)  = 0.14300e-02_wp
            hbb(15,7,3)  = 0.10260e-01_wp
            hbb(15,8,3)  = -0.66500e-02_wp
            hbb(15,9,3)  = -0.36100e-02_wp
            hbb(15,10,3) = 0.14060e-01_wp
            hbb(15,11,3) = -0.65400e-02_wp
            hbb(15,12,3) = -0.75200e-02_wp
            hbb(15,13,3) = -0.39580e-01_wp
            hbb(15,14,3) = -0.39690e-01_wp
            hbb(15,15,3) = 0.79270e-01_wp
            hbb(16,1,3)  = -0.35240e-01_wp
            hbb(16,2,3)  = 0.20060e-01_wp
            hbb(16,3,3)  = 0.15180e-01_wp
            hbb(16,4,3)  = -0.85100e-02_wp
            hbb(16,5,3)  = -0.49000e-03_wp
            hbb(16,6,3)  = 0.90000e-02_wp
            hbb(16,7,3)  = 0.56800e-02_wp
            hbb(16,8,3)  = -0.26300e-02_wp
            hbb(16,9,3)  = -0.30500e-02_wp
            hbb(16,10,3) = -0.90000e-02_wp
            hbb(16,11,3) = 0.10190e-01_wp
            hbb(16,12,3) = -0.12000e-02_wp
            hbb(16,13,3) = -0.31490e-01_wp
            hbb(16,14,3) = 0.13200e-01_wp
            hbb(16,15,3) = 0.18290e-01_wp
            hbb(16,16,3) = 0.78560e-01_wp
            hbb(17,1,3)  = 0.20310e-01_wp
            hbb(17,2,3)  = -0.11130e-01_wp
            hbb(17,3,3)  = -0.91900e-02_wp
            hbb(17,4,3)  = -0.11900e-02_wp
            hbb(17,5,3)  = 0.41100e-02_wp
            hbb(17,6,3)  = -0.29200e-02_wp
            hbb(17,7,3)  = -0.26100e-02_wp
            hbb(17,8,3)  = 0.13200e-02_wp
            hbb(17,9,3)  = 0.12900e-02_wp
            hbb(17,10,3) = 0.10050e-01_wp
            hbb(17,11,3) = -0.66300e-02_wp
            hbb(17,12,3) = -0.34200e-02_wp
            hbb(17,13,3) = 0.12780e-01_wp
            hbb(17,14,3) = -0.59900e-02_wp
            hbb(17,15,3) = -0.67900e-02_wp
            hbb(17,16,3) = -0.39350e-01_wp
            hbb(17,17,3) = 0.77640e-01_wp
            hbb(18,1,3)  = 0.14930e-01_wp
            hbb(18,2,3)  = -0.89300e-02_wp
            hbb(18,3,3)  = -0.60000e-02_wp
            hbb(18,4,3)  = 0.97000e-02_wp
            hbb(18,5,3)  = -0.36200e-02_wp
            hbb(18,6,3)  = -0.60800e-02_wp
            hbb(18,7,3)  = -0.30800e-02_wp
            hbb(18,8,3)  = 0.13100e-02_wp
            hbb(18,9,3)  = 0.17700e-02_wp
            hbb(18,10,3) = -0.10500e-02_wp
            hbb(18,11,3) = -0.35600e-02_wp
            hbb(18,12,3) = 0.46100e-02_wp
            hbb(18,13,3) = 0.18710e-01_wp
            hbb(18,14,3) = -0.72100e-02_wp
            hbb(18,15,3) = -0.11500e-01_wp
            hbb(18,16,3) = -0.39210e-01_wp
            hbb(18,17,3) = -0.38290e-01_wp
            hbb(18,18,3) = 0.77500e-01_wp
            ! Hessian elements for anchor point 3
            htt( 1,1,3)  = 0.10290e-01_wp
            htt( 2,1,3)  = 0.63900e-02_wp
            htt( 2,2,3)  = 0.13000e-01_wp
            htt( 3,1,3)  = 0.28000e-03_wp
            htt( 3,2,3)  = -0.46200e-02_wp
            htt( 3,3,3)  = 0.11260e-01_wp
            htt( 4,1,3)  = -0.36100e-02_wp
            htt( 4,2,3)  = 0.19900e-02_wp
            htt( 4,3,3)  = 0.63700e-02_wp
            htt( 4,4,3)  = 0.11970e-01_wp
            htt( 5,1,3)  = -0.93800e-02_wp
            htt( 5,2,3)  = -0.10580e-01_wp
            htt( 5,3,3)  = 0.29200e-02_wp
            htt( 5,4,3)  = 0.17200e-02_wp
            htt( 5,5,3)  = 0.16730e-01_wp
            htt( 6,1,3)  = -0.54800e-02_wp
            htt( 6,2,3)  = -0.54200e-02_wp
            htt( 6,3,3)  = 0.33200e-02_wp
            htt( 6,4,3)  = 0.33700e-02_wp
            htt( 6,5,3)  = 0.66800e-02_wp
            htt( 6,6,3)  = 0.91600e-02_wp
            htt( 7,1,3)  = -0.40000e-04_wp
            htt( 7,2,3)  = -0.31000e-03_wp
            htt( 7,3,3)  = -0.73300e-02_wp
            htt( 7,4,3)  = -0.75900e-02_wp
            htt( 7,5,3)  = 0.52600e-02_wp
            htt( 7,6,3)  = -0.15300e-02_wp
            htt( 7,7,3)  = 0.12050e-01_wp
            htt( 8,1,3)  = 0.38600e-02_wp
            htt( 8,2,3)  = 0.48500e-02_wp
            htt( 8,3,3)  = -0.69300e-02_wp
            htt( 8,4,3)  = -0.59400e-02_wp
            htt( 8,5,3)  = -0.48000e-02_wp
            htt( 8,6,3)  = 0.95000e-03_wp
            htt( 8,7,3)  = 0.52700e-02_wp
            htt( 8,8,3)  = 0.11020e-01_wp
            htt( 9,1,3)  = -0.52600e-02_wp
            htt( 9,2,3)  = 0.85000e-03_wp
            htt( 9,3,3)  = -0.32600e-02_wp
            htt( 9,4,3)  = 0.28400e-02_wp
            htt( 9,5,3)  = -0.16300e-02_wp
            htt( 9,6,3)  = 0.11100e-02_wp
            htt( 9,7,3)  = -0.35000e-02_wp
            htt( 9,8,3)  = -0.76000e-03_wp
            htt( 9,9,3)  = 0.98000e-02_wp
            htt(10,1,3)  = -0.32000e-02_wp
            htt(10,2,3)  = 0.20800e-02_wp
            htt(10,3,3)  = -0.19100e-02_wp
            htt(10,4,3)  = 0.33700e-02_wp
            htt(10,5,3)  = 0.37000e-03_wp
            htt(10,6,3)  = 0.79000e-03_wp
            htt(10,7,3)  = -0.83000e-03_wp
            htt(10,8,3)  = -0.41000e-03_wp
            htt(10,9,3)  = 0.34100e-02_wp
            htt(10,10,3) = 0.86300e-02_wp
            htt(11,1,3)  = -0.12800e-02_wp
            htt(11,2,3)  = -0.59200e-02_wp
            htt(11,3,3)  = 0.17500e-02_wp
            htt(11,4,3)  = -0.28900e-02_wp
            htt(11,5,3)  = -0.41000e-03_wp
            htt(11,6,3)  = 0.10500e-02_wp
            htt(11,7,3)  = -0.32300e-02_wp
            htt(11,8,3)  = -0.17700e-02_wp
            htt(11,9,3)  = 0.35500e-02_wp
            htt(11,10,3) = -0.20000e-02_wp
            htt(11,11,3) = 0.83000e-02_wp
            htt(12,1,3)  = 0.78000e-03_wp
            htt(12,2,3)  = -0.46800e-02_wp
            htt(12,3,3)  = 0.31000e-02_wp
            htt(12,4,3)  = -0.23600e-02_wp
            htt(12,5,3)  = 0.16000e-02_wp
            htt(12,6,3)  = 0.73000e-03_wp
            htt(12,7,3)  = -0.56000e-03_wp
            htt(12,8,3)  = -0.14300e-02_wp
            htt(12,9,3)  = -0.28400e-02_wp
            htt(12,10,3) = 0.32200e-02_wp
            htt(12,11,3) = 0.27500e-02_wp
            htt(12,12,3) = 0.88100e-02_wp
            htt(13,1,3)  = -0.66000e-03_wp
            htt(13,2,3)  = -0.37900e-02_wp
            htt(13,3,3)  = 0.29600e-02_wp
            htt(13,4,3)  = -0.16000e-03_wp
            htt(13,5,3)  = 0.51200e-02_wp
            htt(13,6,3)  = 0.20100e-02_wp
            htt(13,7,3)  = 0.17400e-02_wp
            htt(13,8,3)  = -0.13700e-02_wp
            htt(13,9,3)  = -0.72600e-02_wp
            htt(13,10,3) = -0.74000e-03_wp
            htt(13,11,3) = -0.40600e-02_wp
            htt(13,12,3) = 0.24600e-02_wp
            htt(13,13,3) = 0.10500e-01_wp
            htt(14,1,3)  = 0.16400e-02_wp
            htt(14,2,3)  = -0.60000e-04_wp
            htt(14,3,3)  = 0.14800e-02_wp
            htt(14,4,3)  = -0.22000e-03_wp
            htt(14,5,3)  = 0.17400e-02_wp
            htt(14,6,3)  = -0.85000e-03_wp
            htt(14,7,3)  = 0.18900e-02_wp
            htt(14,8,3)  = -0.70000e-03_wp
            htt(14,9,3)  = -0.45800e-02_wp
            htt(14,10,3) = 0.11000e-02_wp
            htt(14,11,3) = -0.28400e-02_wp
            htt(14,12,3) = 0.28400e-02_wp
            htt(14,13,3) = 0.40400e-02_wp
            htt(14,14,3) = 0.85900e-02_wp
            htt(15,1,3)  = -0.27400e-02_wp
            htt(15,2,3)  = -0.50400e-02_wp
            htt(15,3,3)  = 0.16000e-02_wp
            htt(15,4,3)  = -0.70000e-03_wp
            htt(15,5,3)  = 0.30800e-02_wp
            htt(15,6,3)  = 0.23300e-02_wp
            htt(15,7,3)  = -0.97000e-03_wp
            htt(15,8,3)  = -0.17200e-02_wp
            htt(15,9,3)  = -0.78000e-03_wp
            htt(15,10,3) = -0.60300e-02_wp
            htt(15,11,3) = 0.15700e-02_wp
            htt(15,12,3) = -0.36800e-02_wp
            htt(15,13,3) = 0.39000e-02_wp
            htt(15,14,3) = -0.17100e-02_wp
            htt(15,15,3) = 0.92100e-02_wp
            htt(16,1,3)  = -0.44000e-03_wp
            htt(16,2,3)  = -0.13100e-02_wp
            htt(16,3,3)  = 0.11000e-03_wp
            htt(16,4,3)  = -0.76000e-03_wp
            htt(16,5,3)  = -0.29000e-03_wp
            htt(16,6,3)  = -0.53000e-03_wp
            htt(16,7,3)  = -0.81000e-03_wp
            htt(16,8,3)  = -0.10500e-02_wp
            htt(16,9,3)  = 0.18900e-02_wp
            htt(16,10,3) = -0.41800e-02_wp
            htt(16,11,3) = 0.27800e-02_wp
            htt(16,12,3) = -0.33000e-02_wp
            htt(16,13,3) = -0.25500e-02_wp
            htt(16,14,3) = 0.28300e-02_wp
            htt(16,15,3) = 0.36000e-02_wp
            htt(16,16,3) = 0.89900e-02_wp
            htt(17,1,3)  = 0.15600e-02_wp
            htt(17,2,3)  = -0.56000e-03_wp
            htt(17,3,3)  = 0.34000e-03_wp
            htt(17,4,3)  = -0.17800e-02_wp
            htt(17,5,3)  = 0.25100e-02_wp
            htt(17,6,3)  = -0.75000e-03_wp
            htt(17,7,3)  = 0.36600e-02_wp
            htt(17,8,3)  = 0.39000e-03_wp
            htt(17,9,3)  = 0.16000e-03_wp
            htt(17,10,3) = -0.21700e-02_wp
            htt(17,11,3) = 0.23300e-02_wp
            htt(17,12,3) = 0.00000e+00_wp
            htt(17,13,3) = -0.59100e-02_wp
            htt(17,14,3) = -0.56000e-03_wp
            htt(17,15,3) = -0.35500e-02_wp
            htt(17,16,3) = 0.18000e-02_wp
            htt(17,17,3) = 0.10140e-01_wp
            htt(18,1,3)  = 0.13300e-02_wp
            htt(18,2,3)  = 0.25800e-02_wp
            htt(18,3,3)  = -0.22700e-02_wp
            htt(18,4,3)  = -0.10200e-02_wp
            htt(18,5,3)  = -0.30300e-02_wp
            htt(18,6,3)  = -0.35800e-02_wp
            htt(18,7,3)  = 0.33000e-03_wp
            htt(18,8,3)  = -0.23000e-03_wp
            htt(18,9,3)  = 0.24700e-02_wp
            htt(18,10,3) = 0.75000e-03_wp
            htt(18,11,3) = 0.12000e-02_wp
            htt(18,12,3) = -0.52000e-03_wp
            htt(18,13,3) = -0.44900e-02_wp
            htt(18,14,3) = 0.90000e-03_wp
            htt(18,15,3) = -0.27500e-02_wp
            htt(18,16,3) = 0.26400e-02_wp
            htt(18,17,3) = 0.27400e-02_wp
            htt(18,18,3) = 0.87100e-02_wp
            htt(19,1,3)  = -0.74000e-03_wp
            htt(19,2,3)  = -0.42800e-02_wp
            htt(19,3,3)  = 0.18200e-02_wp
            htt(19,4,3)  = -0.17200e-02_wp
            htt(19,5,3)  = 0.58900e-02_wp
            htt(19,6,3)  = 0.21100e-02_wp
            htt(19,7,3)  = 0.35000e-02_wp
            htt(19,8,3)  = -0.28000e-03_wp
            htt(19,9,3)  = -0.25100e-02_wp
            htt(19,10,3) = -0.40200e-02_wp
            htt(19,11,3) = 0.11200e-02_wp
            htt(19,12,3) = -0.38000e-03_wp
            htt(19,13,3) = 0.54000e-03_wp
            htt(19,14,3) = -0.51100e-02_wp
            htt(19,15,3) = 0.20600e-02_wp
            htt(19,16,3) = -0.35900e-02_wp
            htt(19,17,3) = 0.47900e-02_wp
            htt(19,18,3) = -0.26500e-02_wp
            htt(19,19,3) = 0.10440e-01_wp
            htt(20,1,3)  = -0.97000e-03_wp
            htt(20,2,3)  = -0.11500e-02_wp
            htt(20,3,3)  = -0.78000e-03_wp
            htt(20,4,3)  = -0.96000e-03_wp
            htt(20,5,3)  = 0.35000e-03_wp
            htt(20,6,3)  = -0.72000e-03_wp
            htt(20,7,3)  = 0.18000e-03_wp
            htt(20,8,3)  = -0.89000e-03_wp
            htt(20,9,3)  = -0.20000e-03_wp
            htt(20,10,3) = -0.10900e-02_wp
            htt(20,11,3) = -0.20000e-04_wp
            htt(20,12,3) = -0.90000e-03_wp
            htt(20,13,3) = 0.19600e-02_wp
            htt(20,14,3) = -0.36500e-02_wp
            htt(20,15,3) = 0.28600e-02_wp
            htt(20,16,3) = -0.27500e-02_wp
            htt(20,17,3) = -0.26100e-02_wp
            htt(20,18,3) = 0.33200e-02_wp
            htt(20,19,3) = 0.30100e-02_wp
            htt(20,20,3) = 0.89300e-02_wp
            htt(21,1,3)  = 0.34400e-02_wp
            htt(21,2,3)  = 0.76900e-02_wp
            htt(21,3,3)  = -0.32500e-02_wp
            htt(21,4,3)  = 0.10000e-02_wp
            htt(21,5,3)  = -0.13360e-01_wp
            htt(21,6,3)  = -0.35700e-02_wp
            htt(21,7,3)  = -0.71200e-02_wp
            htt(21,8,3)  = 0.26800e-02_wp
            htt(21,9,3)  = 0.42100e-02_wp
            htt(21,10,3) = 0.23400e-02_wp
            htt(21,11,3) = -0.13000e-03_wp
            htt(21,12,3) = -0.20000e-02_wp
            htt(21,13,3) = -0.18100e-02_wp
            htt(21,14,3) = -0.22900e-02_wp
            htt(21,15,3) = 0.90000e-04_wp
            htt(21,16,3) = -0.39000e-03_wp
            htt(21,17,3) = -0.84600e-02_wp
            htt(21,18,3) = 0.98000e-03_wp
            htt(21,19,3) = -0.79800e-02_wp
            htt(21,20,3) = 0.14600e-02_wp
            htt(21,21,3) = 0.15990e-01_wp
            htt(22,1,3)  = -0.55000e-03_wp
            htt(22,2,3)  = 0.24000e-02_wp
            htt(22,3,3)  = -0.36500e-02_wp
            htt(22,4,3)  = -0.70000e-03_wp
            htt(22,5,3)  = -0.30600e-02_wp
            htt(22,6,3)  = -0.61100e-02_wp
            htt(22,7,3)  = -0.17000e-03_wp
            htt(22,8,3)  = -0.32100e-02_wp
            htt(22,9,3)  = 0.14100e-02_wp
            htt(22,10,3) = 0.19200e-02_wp
            htt(22,11,3) = -0.16200e-02_wp
            htt(22,12,3) = -0.11100e-02_wp
            htt(22,13,3) = 0.13700e-02_wp
            htt(22,14,3) = 0.37000e-03_wp
            htt(22,15,3) = 0.86000e-03_wp
            htt(22,16,3) = -0.15000e-03_wp
            htt(22,17,3) = -0.51200e-02_wp
            htt(22,18,3) = 0.15500e-02_wp
            htt(22,19,3) = -0.41100e-02_wp
            htt(22,20,3) = 0.25600e-02_wp
            htt(22,21,3) = 0.59500e-02_wp
            htt(22,22,3) = 0.90700e-02_wp
            htt(23,1,3)  = 0.36700e-02_wp
            htt(23,2,3)  = 0.46000e-02_wp
            htt(23,3,3)  = -0.68000e-03_wp
            htt(23,4,3)  = 0.25000e-03_wp
            htt(23,5,3)  = -0.79000e-02_wp
            htt(23,6,3)  = -0.77000e-03_wp
            htt(23,7,3)  = -0.38400e-02_wp
            htt(23,8,3)  = 0.32800e-02_wp
            htt(23,9,3)  = 0.19400e-02_wp
            htt(23,10,3) = -0.54000e-03_wp
            htt(23,11,3) = 0.99000e-03_wp
            htt(23,12,3) = -0.14900e-02_wp
            htt(23,13,3) = -0.32100e-02_wp
            htt(23,14,3) = -0.37300e-02_wp
            htt(23,15,3) = -0.70000e-03_wp
            htt(23,16,3) = -0.12200e-02_wp
            htt(23,17,3) = -0.11700e-02_wp
            htt(23,18,3) = -0.49000e-02_wp
            htt(23,19,3) = -0.65000e-03_wp
            htt(23,20,3) = -0.43800e-02_wp
            htt(23,21,3) = 0.66800e-02_wp
            htt(23,22,3) = -0.62000e-03_wp
            htt(23,23,3) = 0.10360e-01_wp
            htt(24,1,3)  = -0.33000e-03_wp
            htt(24,2,3)  = -0.69000e-03_wp
            htt(24,3,3)  = -0.10800e-02_wp
            htt(24,4,3)  = -0.14500e-02_wp
            htt(24,5,3)  = 0.24000e-02_wp
            htt(24,6,3)  = -0.33100e-02_wp
            htt(24,7,3)  = 0.31100e-02_wp
            htt(24,8,3)  = -0.26100e-02_wp
            htt(24,9,3)  = -0.87000e-03_wp
            htt(24,10,3) = -0.97000e-03_wp
            htt(24,11,3) = -0.50000e-03_wp
            htt(24,12,3) = -0.60000e-03_wp
            htt(24,13,3) = -0.30000e-04_wp
            htt(24,14,3) = -0.10800e-02_wp
            htt(24,15,3) = 0.70000e-04_wp
            htt(24,16,3) = -0.98000e-03_wp
            htt(24,17,3) = 0.21800e-02_wp
            htt(24,18,3) = -0.43300e-02_wp
            htt(24,19,3) = 0.32200e-02_wp
            htt(24,20,3) = -0.32900e-02_wp
            htt(24,21,3) = -0.33600e-02_wp
            htt(24,22,3) = 0.25000e-02_wp
            htt(24,23,3) = 0.30700e-02_wp
            htt(24,24,3) = 0.89200e-02_wp
            ! Hessian elements for anchor point 3
            hbs( 1,1,3)  = -0.35620e-01_wp
            hbs( 1,2,3)  = -0.31870e-01_wp
            hbs( 1,3,3)  = -0.41940e-01_wp
            hbs( 1,4,3)  = -0.11410e-01_wp
            hbs( 1,5,3)  = 0.68300e-02_wp
            hbs( 1,6,3)  = 0.45490e-01_wp
            hbs( 1,7,3)  = -0.31000e-03_wp
            hbs( 1,8,3)  = 0.44930e-01_wp
            hbs( 1,9,3)  = 0.00000e+00_wp
            hbs( 1,10,3) = -0.10110e-01_wp
            hbs( 1,11,3) = -0.33000e-03_wp
            hbs( 1,12,3) = 0.43300e-02_wp
            hbs( 2,1,3)  = 0.49480e-01_wp
            hbs( 2,2,3)  = -0.69900e-02_wp
            hbs( 2,3,3)  = 0.24890e-01_wp
            hbs( 2,4,3)  = -0.87000e-03_wp
            hbs( 2,5,3)  = -0.66200e-02_wp
            hbs( 2,6,3)  = -0.20260e-01_wp
            hbs( 2,7,3)  = 0.89000e-03_wp
            hbs( 2,8,3)  = -0.28090e-01_wp
            hbs( 2,9,3)  = -0.10000e-04_wp
            hbs( 2,10,3) = 0.10450e-01_wp
            hbs( 2,11,3) = -0.46000e-03_wp
            hbs( 2,12,3) = -0.19500e-02_wp
            hbs( 3,1,3)  = -0.13860e-01_wp
            hbs( 3,2,3)  = 0.38860e-01_wp
            hbs( 3,3,3)  = 0.17050e-01_wp
            hbs( 3,4,3)  = 0.12280e-01_wp
            hbs( 3,5,3)  = -0.21000e-03_wp
            hbs( 3,6,3)  = -0.25230e-01_wp
            hbs( 3,7,3)  = -0.58000e-03_wp
            hbs( 3,8,3)  = -0.16840e-01_wp
            hbs( 3,9,3)  = 0.20000e-04_wp
            hbs( 3,10,3) = -0.34000e-03_wp
            hbs( 3,11,3) = 0.79000e-03_wp
            hbs( 3,12,3) = -0.23700e-02_wp
            hbs( 4,1,3)  = -0.28260e-01_wp
            hbs( 4,2,3)  = -0.99500e-02_wp
            hbs( 4,3,3)  = 0.16700e-01_wp
            hbs( 4,4,3)  = -0.32020e-01_wp
            hbs( 4,5,3)  = -0.96100e-02_wp
            hbs( 4,6,3)  = -0.12920e-01_wp
            hbs( 4,7,3)  = 0.48000e-02_wp
            hbs( 4,8,3)  = 0.45530e-01_wp
            hbs( 4,9,3)  = -0.62000e-03_wp
            hbs( 4,10,3) = 0.44050e-01_wp
            hbs( 4,11,3) = 0.10000e-03_wp
            hbs( 4,12,3) = -0.40000e-03_wp
            hbs( 5,1,3)  = 0.30480e-01_wp
            hbs( 5,2,3)  = 0.69800e-02_wp
            hbs( 5,3,3)  = -0.14860e-01_wp
            hbs( 5,4,3)  = 0.26800e-02_wp
            hbs( 5,5,3)  = 0.62000e-02_wp
            hbs( 5,6,3)  = 0.57700e-02_wp
            hbs( 5,7,3)  = -0.18500e-02_wp
            hbs( 5,8,3)  = -0.20690e-01_wp
            hbs( 5,9,3)  = -0.33000e-03_wp
            hbs( 5,10,3) = -0.24600e-01_wp
            hbs( 5,11,3) = -0.60000e-04_wp
            hbs( 5,12,3) = 0.44000e-03_wp
            hbs( 6,1,3)  = -0.22300e-02_wp
            hbs( 6,2,3)  = 0.29700e-02_wp
            hbs( 6,3,3)  = -0.18400e-02_wp
            hbs( 6,4,3)  = 0.29340e-01_wp
            hbs( 6,5,3)  = 0.34100e-02_wp
            hbs( 6,6,3)  = 0.71600e-02_wp
            hbs( 6,7,3)  = -0.29500e-02_wp
            hbs( 6,8,3)  = -0.24830e-01_wp
            hbs( 6,9,3)  = 0.94000e-03_wp
            hbs( 6,10,3) = -0.19440e-01_wp
            hbs( 6,11,3) = -0.40000e-04_wp
            hbs( 6,12,3) = -0.30000e-04_wp
            hbs( 7,1,3)  = -0.12340e-01_wp
            hbs( 7,2,3)  = 0.39640e-01_wp
            hbs( 7,3,3)  = 0.73200e-02_wp
            hbs( 7,4,3)  = -0.32120e-01_wp
            hbs( 7,5,3)  = 0.49300e-02_wp
            hbs( 7,6,3)  = -0.33370e-01_wp
            hbs( 7,7,3)  = -0.89700e-02_wp
            hbs( 7,8,3)  = -0.12760e-01_wp
            hbs( 7,9,3)  = 0.49400e-02_wp
            hbs( 7,10,3) = 0.45710e-01_wp
            hbs( 7,11,3) = -0.89000e-03_wp
            hbs( 7,12,3) = 0.50000e-04_wp
            hbs( 8,1,3)  = 0.73400e-02_wp
            hbs( 8,2,3)  = -0.22220e-01_wp
            hbs( 8,3,3)  = -0.35200e-02_wp
            hbs( 8,4,3)  = 0.29850e-01_wp
            hbs( 8,5,3)  = -0.30000e-02_wp
            hbs( 8,6,3)  = 0.38900e-02_wp
            hbs( 8,7,3)  = 0.48800e-02_wp
            hbs( 8,8,3)  = 0.56300e-02_wp
            hbs( 8,9,3)  = -0.19800e-02_wp
            hbs( 8,10,3) = -0.20400e-01_wp
            hbs( 8,11,3) = -0.22000e-03_wp
            hbs( 8,12,3) = -0.50000e-04_wp
            hbs( 9,1,3)  = 0.50100e-02_wp
            hbs( 9,2,3)  = -0.17420e-01_wp
            hbs( 9,3,3)  = -0.38000e-02_wp
            hbs( 9,4,3)  = 0.22700e-02_wp
            hbs( 9,5,3)  = -0.19200e-02_wp
            hbs( 9,6,3)  = 0.29480e-01_wp
            hbs( 9,7,3)  = 0.40900e-02_wp
            hbs( 9,8,3)  = 0.71400e-02_wp
            hbs( 9,9,3)  = -0.29600e-02_wp
            hbs( 9,10,3) = -0.25300e-01_wp
            hbs( 9,11,3) = 0.11100e-02_wp
            hbs( 9,12,3) = -0.10000e-04_wp
            hbs(10,1,3)  = 0.42740e-01_wp
            hbs(10,2,3)  = 0.42010e-01_wp
            hbs(10,3,3)  = -0.32900e-02_wp
            hbs(10,4,3)  = -0.13800e-01_wp
            hbs(10,5,3)  = -0.37000e-03_wp
            hbs(10,6,3)  = -0.32450e-01_wp
            hbs(10,7,3)  = 0.53400e-02_wp
            hbs(10,8,3)  = -0.31850e-01_wp
            hbs(10,9,3)  = -0.86200e-02_wp
            hbs(10,10,3) = -0.13500e-01_wp
            hbs(10,11,3) = 0.52600e-02_wp
            hbs(10,12,3) = -0.66000e-03_wp
            hbs(11,1,3)  = -0.23640e-01_wp
            hbs(11,2,3)  = -0.18640e-01_wp
            hbs(11,3,3)  = 0.15000e-02_wp
            hbs(11,4,3)  = 0.71800e-02_wp
            hbs(11,5,3)  = 0.84000e-03_wp
            hbs(11,6,3)  = 0.30210e-01_wp
            hbs(11,7,3)  = -0.32100e-02_wp
            hbs(11,8,3)  = 0.20600e-02_wp
            hbs(11,9,3)  = 0.42700e-02_wp
            hbs(11,10,3) = 0.63000e-02_wp
            hbs(11,11,3) = -0.21200e-02_wp
            hbs(11,12,3) = -0.27000e-03_wp
            hbs(12,1,3)  = -0.19100e-01_wp
            hbs(12,2,3)  = -0.23380e-01_wp
            hbs(12,3,3)  = 0.17900e-02_wp
            hbs(12,4,3)  = 0.66100e-02_wp
            hbs(12,5,3)  = -0.47000e-03_wp
            hbs(12,6,3)  = 0.22400e-02_wp
            hbs(12,7,3)  = -0.21200e-02_wp
            hbs(12,8,3)  = 0.29790e-01_wp
            hbs(12,9,3)  = 0.43500e-02_wp
            hbs(12,10,3) = 0.72000e-02_wp
            hbs(12,11,3) = -0.31400e-02_wp
            hbs(12,12,3) = 0.93000e-03_wp
            hbs(13,1,3)  = 0.40850e-01_wp
            hbs(13,2,3)  = -0.14500e-01_wp
            hbs(13,3,3)  = 0.10130e-01_wp
            hbs(13,4,3)  = 0.45630e-01_wp
            hbs(13,5,3)  = -0.20000e-03_wp
            hbs(13,6,3)  = -0.12730e-01_wp
            hbs(13,7,3)  = -0.96000e-03_wp
            hbs(13,8,3)  = -0.32840e-01_wp
            hbs(13,9,3)  = 0.49400e-02_wp
            hbs(13,10,3) = -0.32520e-01_wp
            hbs(13,11,3) = -0.89800e-02_wp
            hbs(13,12,3) = 0.51000e-02_wp
            hbs(14,1,3)  = -0.17990e-01_wp
            hbs(14,2,3)  = 0.57300e-02_wp
            hbs(14,3,3)  = -0.40700e-02_wp
            hbs(14,4,3)  = -0.25200e-01_wp
            hbs(14,5,3)  = 0.12000e-03_wp
            hbs(14,6,3)  = 0.71400e-02_wp
            hbs(14,7,3)  = 0.11600e-02_wp
            hbs(14,8,3)  = 0.29210e-01_wp
            hbs(14,9,3)  = -0.29300e-02_wp
            hbs(14,10,3) = 0.21300e-02_wp
            hbs(14,11,3) = 0.42700e-02_wp
            hbs(14,12,3) = -0.19500e-02_wp
            hbs(15,1,3)  = -0.22870e-01_wp
            hbs(15,2,3)  = 0.87800e-02_wp
            hbs(15,3,3)  = -0.60600e-02_wp
            hbs(15,4,3)  = -0.20420e-01_wp
            hbs(15,5,3)  = 0.70000e-04_wp
            hbs(15,6,3)  = 0.55900e-02_wp
            hbs(15,7,3)  = -0.20000e-03_wp
            hbs(15,8,3)  = 0.36300e-02_wp
            hbs(15,9,3)  = -0.20100e-02_wp
            hbs(15,10,3) = 0.30390e-01_wp
            hbs(15,11,3) = 0.47100e-02_wp
            hbs(15,12,3) = -0.31500e-02_wp
            hbs(16,1,3)  = -0.73800e-02_wp
            hbs(16,2,3)  = -0.25320e-01_wp
            hbs(16,3,3)  = 0.11080e-01_wp
            hbs(16,4,3)  = 0.43720e-01_wp
            hbs(16,5,3)  = -0.15800e-02_wp
            hbs(16,6,3)  = 0.45980e-01_wp
            hbs(16,7,3)  = 0.11000e-03_wp
            hbs(16,8,3)  = -0.13000e-01_wp
            hbs(16,9,3)  = -0.64000e-03_wp
            hbs(16,10,3) = -0.33620e-01_wp
            hbs(16,11,3) = 0.48400e-02_wp
            hbs(16,12,3) = -0.84200e-02_wp
            hbs(17,1,3)  = 0.40600e-02_wp
            hbs(17,2,3)  = 0.27060e-01_wp
            hbs(17,3,3)  = -0.92100e-02_wp
            hbs(17,4,3)  = -0.24320e-01_wp
            hbs(17,5,3)  = 0.15600e-02_wp
            hbs(17,6,3)  = -0.20640e-01_wp
            hbs(17,7,3)  = -0.50000e-04_wp
            hbs(17,8,3)  = 0.61500e-02_wp
            hbs(17,9,3)  = -0.30000e-03_wp
            hbs(17,10,3) = 0.36700e-02_wp
            hbs(17,11,3) = -0.18400e-02_wp
            hbs(17,12,3) = 0.60700e-02_wp
            hbs(18,1,3)  = 0.33200e-02_wp
            hbs(18,2,3)  = -0.17400e-02_wp
            hbs(18,3,3)  = -0.18800e-02_wp
            hbs(18,4,3)  = -0.19400e-01_wp
            hbs(18,5,3)  = 0.20000e-04_wp
            hbs(18,6,3)  = -0.25340e-01_wp
            hbs(18,7,3)  = -0.50000e-04_wp
            hbs(18,8,3)  = 0.68600e-02_wp
            hbs(18,9,3)  = 0.93000e-03_wp
            hbs(18,10,3) = 0.29960e-01_wp
            hbs(18,11,3) = -0.30000e-02_wp
            hbs(18,12,3) = 0.23500e-02_wp
            ! Hessian elements for anchor point 4
            hss( 1,1,4)  = 0.38571e+00_wp
            hss( 2,1,4)  = 0.12420e-01_wp
            hss( 2,2,4)  = 0.38571e+00_wp
            hss( 3,1,4)  = 0.51260e-01_wp
            hss( 3,2,4)  = 0.51260e-01_wp
            hss( 3,3,4)  = 0.38504e+00_wp
            hss( 4,1,4)  = 0.20570e-01_wp
            hss( 4,2,4)  = -0.17300e-02_wp
            hss( 4,3,4)  = -0.32200e-02_wp
            hss( 4,4,4)  = 0.41668e+00_wp
            hss( 5,1,4)  = 0.50700e-02_wp
            hss( 5,2,4)  = 0.48000e-03_wp
            hss( 5,3,4)  = -0.24900e-02_wp
            hss( 5,4,4)  = 0.39900e-02_wp
            hss( 5,5,4)  = 0.36905e+00_wp
            hss( 6,1,4)  = -0.36300e-02_wp
            hss( 6,2,4)  = 0.82190e-01_wp
            hss( 6,3,4)  = -0.13300e-02_wp
            hss( 6,4,4)  = 0.21260e-01_wp
            hss( 6,5,4)  = -0.91000e-03_wp
            hss( 6,6,4)  = 0.41243e+00_wp
            hss( 7,1,4)  = -0.51000e-03_wp
            hss( 7,2,4)  = -0.39000e-03_wp
            hss( 7,3,4)  = 0.41000e-03_wp
            hss( 7,4,4)  = 0.44700e-02_wp
            hss( 7,5,4)  = 0.43000e-03_wp
            hss( 7,6,4)  = 0.45000e-02_wp
            hss( 7,7,4)  = 0.36452e+00_wp
            hss( 8,1,4)  = 0.82190e-01_wp
            hss( 8,2,4)  = -0.36300e-02_wp
            hss( 8,3,4)  = -0.13300e-02_wp
            hss( 8,4,4)  = -0.41800e-02_wp
            hss( 8,5,4)  = -0.36000e-03_wp
            hss( 8,6,4)  = 0.20820e-01_wp
            hss( 8,7,4)  = -0.28000e-03_wp
            hss( 8,8,4)  = 0.41243e+00_wp
            hss( 9,1,4)  = -0.39000e-03_wp
            hss( 9,2,4)  = -0.39000e-03_wp
            hss( 9,3,4)  = -0.57000e-03_wp
            hss( 9,4,4)  = -0.62000e-03_wp
            hss( 9,5,4)  = 0.12000e-03_wp
            hss( 9,6,4)  = 0.43600e-02_wp
            hss( 9,7,4)  = 0.43000e-03_wp
            hss( 9,8,4)  = 0.43600e-02_wp
            hss( 9,9,4)  = 0.36801e+00_wp
            hss(10,1,4)  = -0.17300e-02_wp
            hss(10,2,4)  = 0.20570e-01_wp
            hss(10,3,4)  = -0.32200e-02_wp
            hss(10,4,4)  = 0.86500e-01_wp
            hss(10,5,4)  = -0.80000e-03_wp
            hss(10,6,4)  = -0.41800e-02_wp
            hss(10,7,4)  = -0.54000e-03_wp
            hss(10,8,4)  = 0.21260e-01_wp
            hss(10,9,4)  = -0.62000e-03_wp
            hss(10,10,4) = 0.41668e+00_wp
            hss(11,1,4)  = -0.39000e-03_wp
            hss(11,2,4)  = -0.51000e-03_wp
            hss(11,3,4)  = 0.41000e-03_wp
            hss(11,4,4)  = -0.54000e-03_wp
            hss(11,5,4)  = 0.90000e-04_wp
            hss(11,6,4)  = -0.28000e-03_wp
            hss(11,7,4)  = 0.20000e-03_wp
            hss(11,8,4)  = 0.45000e-02_wp
            hss(11,9,4)  = 0.43000e-03_wp
            hss(11,10,4) = 0.44700e-02_wp
            hss(11,11,4) = 0.36452e+00_wp
            hss(12,1,4)  = 0.48000e-03_wp
            hss(12,2,4)  = 0.50700e-02_wp
            hss(12,3,4)  = -0.24900e-02_wp
            hss(12,4,4)  = -0.80000e-03_wp
            hss(12,5,4)  = -0.30000e-04_wp
            hss(12,6,4)  = -0.36000e-03_wp
            hss(12,7,4)  = 0.90000e-04_wp
            hss(12,8,4)  = -0.91000e-03_wp
            hss(12,9,4)  = 0.12000e-03_wp
            hss(12,10,4) = 0.39900e-02_wp
            hss(12,11,4) = 0.43000e-03_wp
            hss(12,12,4) = 0.36905e+00_wp
            ! Hessian elements for anchor point 4
            hbb( 1,1,4)  = 0.80030e-01_wp
            hbb( 2,1,4)  = -0.40020e-01_wp
            hbb( 2,2,4)  = 0.92960e-01_wp
            hbb( 3,1,4)  = -0.40020e-01_wp
            hbb( 3,2,4)  = -0.52950e-01_wp
            hbb( 3,3,4)  = 0.92960e-01_wp
            hbb( 4,1,4)  = -0.33740e-01_wp
            hbb( 4,2,4)  = 0.17860e-01_wp
            hbb( 4,3,4)  = 0.15880e-01_wp
            hbb( 4,4,4)  = 0.77910e-01_wp
            hbb( 5,1,4)  = 0.19460e-01_wp
            hbb( 5,2,4)  = -0.10850e-01_wp
            hbb( 5,3,4)  = -0.86100e-02_wp
            hbb( 5,4,4)  = -0.38850e-01_wp
            hbb( 5,5,4)  = 0.77000e-01_wp
            hbb( 6,1,4)  = 0.14280e-01_wp
            hbb( 6,2,4)  = -0.70200e-02_wp
            hbb( 6,3,4)  = -0.72700e-02_wp
            hbb( 6,4,4)  = -0.39060e-01_wp
            hbb( 6,5,4)  = -0.38150e-01_wp
            hbb( 6,6,4)  = 0.77210e-01_wp
            hbb( 7,1,4)  = -0.89600e-02_wp
            hbb( 7,2,4)  = -0.20000e-04_wp
            hbb( 7,3,4)  = 0.89900e-02_wp
            hbb( 7,4,4)  = -0.32200e-01_wp
            hbb( 7,5,4)  = 0.13030e-01_wp
            hbb( 7,6,4)  = 0.19170e-01_wp
            hbb( 7,7,4)  = 0.79490e-01_wp
            hbb( 8,1,4)  = -0.12600e-02_wp
            hbb( 8,2,4)  = 0.34000e-02_wp
            hbb( 8,3,4)  = -0.21400e-02_wp
            hbb( 8,4,4)  = 0.18730e-01_wp
            hbb( 8,5,4)  = -0.69300e-02_wp
            hbb( 8,6,4)  = -0.11800e-01_wp
            hbb( 8,7,4)  = -0.39640e-01_wp
            hbb( 8,8,4)  = 0.79360e-01_wp
            hbb( 9,1,4)  = 0.10220e-01_wp
            hbb( 9,2,4)  = -0.33800e-02_wp
            hbb( 9,3,4)  = -0.68500e-02_wp
            hbb( 9,4,4)  = 0.13460e-01_wp
            hbb( 9,5,4)  = -0.60900e-02_wp
            hbb( 9,6,4)  = -0.73700e-02_wp
            hbb( 9,7,4)  = -0.39850e-01_wp
            hbb( 9,8,4)  = -0.39720e-01_wp
            hbb( 9,9,4)  = 0.79570e-01_wp
            hbb(10,1,4)  = 0.53800e-02_wp
            hbb(10,2,4)  = -0.26900e-02_wp
            hbb(10,3,4)  = -0.26900e-02_wp
            hbb(10,4,4)  = -0.85500e-02_wp
            hbb(10,5,4)  = 0.98800e-02_wp
            hbb(10,6,4)  = -0.13400e-02_wp
            hbb(10,7,4)  = -0.34400e-01_wp
            hbb(10,8,4)  = 0.14090e-01_wp
            hbb(10,9,4)  = 0.20310e-01_wp
            hbb(10,10,4) = 0.80510e-01_wp
            hbb(11,1,4)  = -0.26900e-02_wp
            hbb(11,2,4)  = 0.16300e-02_wp
            hbb(11,3,4)  = 0.10600e-02_wp
            hbb(11,4,4)  = -0.14300e-02_wp
            hbb(11,5,4)  = -0.33500e-02_wp
            hbb(11,6,4)  = 0.47900e-02_wp
            hbb(11,7,4)  = 0.20180e-01_wp
            hbb(11,8,4)  = -0.75700e-02_wp
            hbb(11,9,4)  = -0.12610e-01_wp
            hbb(11,10,4) = -0.40260e-01_wp
            hbb(11,11,4) = 0.79330e-01_wp
            hbb(12,1,4)  = -0.26900e-02_wp
            hbb(12,2,4)  = 0.10600e-02_wp
            hbb(12,3,4)  = 0.16300e-02_wp
            hbb(12,4,4)  = 0.99800e-02_wp
            hbb(12,5,4)  = -0.65300e-02_wp
            hbb(12,6,4)  = -0.34500e-02_wp
            hbb(12,7,4)  = 0.14220e-01_wp
            hbb(12,8,4)  = -0.65200e-02_wp
            hbb(12,9,4)  = -0.77100e-02_wp
            hbb(12,10,4) = -0.40260e-01_wp
            hbb(12,11,4) = -0.39080e-01_wp
            hbb(12,12,4) = 0.79330e-01_wp
            hbb(13,1,4)  = -0.89600e-02_wp
            hbb(13,2,4)  = 0.89900e-02_wp
            hbb(13,3,4)  = -0.20000e-04_wp
            hbb(13,4,4)  = 0.49400e-02_wp
            hbb(13,5,4)  = -0.20900e-02_wp
            hbb(13,6,4)  = -0.28500e-02_wp
            hbb(13,7,4)  = -0.88700e-02_wp
            hbb(13,8,4)  = 0.10300e-01_wp
            hbb(13,9,4)  = -0.14300e-02_wp
            hbb(13,10,4) = -0.34400e-01_wp
            hbb(13,11,4) = 0.14220e-01_wp
            hbb(13,12,4) = 0.20180e-01_wp
            hbb(13,13,4) = 0.79490e-01_wp
            hbb(14,1,4)  = 0.10220e-01_wp
            hbb(14,2,4)  = -0.68500e-02_wp
            hbb(14,3,4)  = -0.33800e-02_wp
            hbb(14,4,4)  = -0.27200e-02_wp
            hbb(14,5,4)  = 0.10700e-02_wp
            hbb(14,6,4)  = 0.16500e-02_wp
            hbb(14,7,4)  = -0.14300e-02_wp
            hbb(14,8,4)  = -0.36200e-02_wp
            hbb(14,9,4)  = 0.50500e-02_wp
            hbb(14,10,4) = 0.20310e-01_wp
            hbb(14,11,4) = -0.77100e-02_wp
            hbb(14,12,4) = -0.12610e-01_wp
            hbb(14,13,4) = -0.39850e-01_wp
            hbb(14,14,4) = 0.79570e-01_wp
            hbb(15,1,4)  = -0.12600e-02_wp
            hbb(15,2,4)  = -0.21400e-02_wp
            hbb(15,3,4)  = 0.34000e-02_wp
            hbb(15,4,4)  = -0.22200e-02_wp
            hbb(15,5,4)  = 0.10200e-02_wp
            hbb(15,6,4)  = 0.12000e-02_wp
            hbb(15,7,4)  = 0.10300e-01_wp
            hbb(15,8,4)  = -0.66700e-02_wp
            hbb(15,9,4)  = -0.36200e-02_wp
            hbb(15,10,4) = 0.14090e-01_wp
            hbb(15,11,4) = -0.65200e-02_wp
            hbb(15,12,4) = -0.75700e-02_wp
            hbb(15,13,4) = -0.39640e-01_wp
            hbb(15,14,4) = -0.39720e-01_wp
            hbb(15,15,4) = 0.79360e-01_wp
            hbb(16,1,4)  = -0.33740e-01_wp
            hbb(16,2,4)  = 0.15880e-01_wp
            hbb(16,3,4)  = 0.17860e-01_wp
            hbb(16,4,4)  = -0.83600e-02_wp
            hbb(16,5,4)  = -0.14300e-02_wp
            hbb(16,6,4)  = 0.97900e-02_wp
            hbb(16,7,4)  = 0.49400e-02_wp
            hbb(16,8,4)  = -0.22200e-02_wp
            hbb(16,9,4)  = -0.27200e-02_wp
            hbb(16,10,4) = -0.85500e-02_wp
            hbb(16,11,4) = 0.99800e-02_wp
            hbb(16,12,4) = -0.14300e-02_wp
            hbb(16,13,4) = -0.32200e-01_wp
            hbb(16,14,4) = 0.13460e-01_wp
            hbb(16,15,4) = 0.18730e-01_wp
            hbb(16,16,4) = 0.77910e-01_wp
            hbb(17,1,4)  = 0.19460e-01_wp
            hbb(17,2,4)  = -0.86100e-02_wp
            hbb(17,3,4)  = -0.10850e-01_wp
            hbb(17,4,4)  = -0.14300e-02_wp
            hbb(17,5,4)  = 0.48200e-02_wp
            hbb(17,6,4)  = -0.33900e-02_wp
            hbb(17,7,4)  = -0.20900e-02_wp
            hbb(17,8,4)  = 0.10200e-02_wp
            hbb(17,9,4)  = 0.10700e-02_wp
            hbb(17,10,4) = 0.98800e-02_wp
            hbb(17,11,4) = -0.65300e-02_wp
            hbb(17,12,4) = -0.33500e-02_wp
            hbb(17,13,4) = 0.13030e-01_wp
            hbb(17,14,4) = -0.60900e-02_wp
            hbb(17,15,4) = -0.69300e-02_wp
            hbb(17,16,4) = -0.38850e-01_wp
            hbb(17,17,4) = 0.77000e-01_wp
            hbb(18,1,4)  = 0.14280e-01_wp
            hbb(18,2,4)  = -0.72700e-02_wp
            hbb(18,3,4)  = -0.70200e-02_wp
            hbb(18,4,4)  = 0.97900e-02_wp
            hbb(18,5,4)  = -0.33900e-02_wp
            hbb(18,6,4)  = -0.63900e-02_wp
            hbb(18,7,4)  = -0.28500e-02_wp
            hbb(18,8,4)  = 0.12000e-02_wp
            hbb(18,9,4)  = 0.16500e-02_wp
            hbb(18,10,4) = -0.13400e-02_wp
            hbb(18,11,4) = -0.34500e-02_wp
            hbb(18,12,4) = 0.47900e-02_wp
            hbb(18,13,4) = 0.19170e-01_wp
            hbb(18,14,4) = -0.73700e-02_wp
            hbb(18,15,4) = -0.11800e-01_wp
            hbb(18,16,4) = -0.39060e-01_wp
            hbb(18,17,4) = -0.38150e-01_wp
            hbb(18,18,4) = 0.77210e-01_wp
            ! Hessian elements for anchor point 4
            htt( 1,1,4)  = 0.84400e-02_wp
            htt( 2,1,4)  = 0.31800e-02_wp
            htt( 2,2,4)  = 0.76400e-02_wp
            htt( 3,1,4)  = 0.17500e-02_wp
            htt( 3,2,4)  = -0.29500e-02_wp
            htt( 3,3,4)  = 0.10070e-01_wp
            htt( 4,1,4)  = -0.35100e-02_wp
            htt( 4,2,4)  = 0.15200e-02_wp
            htt( 4,3,4)  = 0.53700e-02_wp
            htt( 4,4,4)  = 0.10400e-01_wp
            htt( 5,1,4)  = -0.54900e-02_wp
            htt( 5,2,4)  = -0.37900e-02_wp
            htt( 5,3,4)  = 0.12100e-02_wp
            htt( 5,4,4)  = 0.29000e-02_wp
            htt( 5,5,4)  = 0.84400e-02_wp
            htt( 6,1,4)  = -0.37900e-02_wp
            htt( 6,2,4)  = -0.24900e-02_wp
            htt( 6,3,4)  = 0.23300e-02_wp
            htt( 6,4,4)  = 0.36400e-02_wp
            htt( 6,5,4)  = 0.31800e-02_wp
            htt( 6,6,4)  = 0.76400e-02_wp
            htt( 7,1,4)  = 0.12100e-02_wp
            htt( 7,2,4)  = 0.23300e-02_wp
            htt( 7,3,4)  = -0.71100e-02_wp
            htt( 7,4,4)  = -0.59900e-02_wp
            htt( 7,5,4)  = 0.17500e-02_wp
            htt( 7,6,4)  = -0.29500e-02_wp
            htt( 7,7,4)  = 0.10070e-01_wp
            htt( 8,1,4)  = 0.29000e-02_wp
            htt( 8,2,4)  = 0.36400e-02_wp
            htt( 8,3,4)  = -0.59900e-02_wp
            htt( 8,4,4)  = -0.52500e-02_wp
            htt( 8,5,4)  = -0.35100e-02_wp
            htt( 8,6,4)  = 0.15200e-02_wp
            htt( 8,7,4)  = 0.53700e-02_wp
            htt( 8,8,4)  = 0.10400e-01_wp
            htt( 9,1,4)  = -0.60400e-02_wp
            htt( 9,2,4)  = -0.50000e-03_wp
            htt( 9,3,4)  = -0.37000e-02_wp
            htt( 9,4,4)  = 0.18400e-02_wp
            htt( 9,5,4)  = 0.10000e-04_wp
            htt( 9,6,4)  = 0.17500e-02_wp
            htt( 9,7,4)  = -0.23300e-02_wp
            htt( 9,8,4)  = -0.59000e-03_wp
            htt( 9,9,4)  = 0.94900e-02_wp
            htt(10,1,4)  = -0.38600e-02_wp
            htt(10,2,4)  = 0.97000e-03_wp
            htt(10,3,4)  = -0.21900e-02_wp
            htt(10,4,4)  = 0.26400e-02_wp
            htt(10,5,4)  = 0.17900e-02_wp
            htt(10,6,4)  = 0.13600e-02_wp
            htt(10,7,4)  = 0.12000e-03_wp
            htt(10,8,4)  = -0.31000e-03_wp
            htt(10,9,4)  = 0.31000e-02_wp
            htt(10,10,4) = 0.84100e-02_wp
            htt(11,1,4)  = -0.65000e-03_wp
            htt(11,2,4)  = -0.50600e-02_wp
            htt(11,3,4)  = 0.11100e-02_wp
            htt(11,4,4)  = -0.33000e-02_wp
            htt(11,5,4)  = -0.17200e-02_wp
            htt(11,6,4)  = 0.41000e-03_wp
            htt(11,7,4)  = -0.34800e-02_wp
            htt(11,8,4)  = -0.13500e-02_wp
            htt(11,9,4)  = 0.38300e-02_wp
            htt(11,10,4) = -0.18300e-02_wp
            htt(11,11,4) = 0.83400e-02_wp
            htt(12,1,4)  = 0.15300e-02_wp
            htt(12,2,4)  = -0.35900e-02_wp
            htt(12,3,4)  = 0.26200e-02_wp
            htt(12,4,4)  = -0.25000e-02_wp
            htt(12,5,4)  = 0.60000e-04_wp
            htt(12,6,4)  = 0.20000e-04_wp
            htt(12,7,4)  = -0.10300e-02_wp
            htt(12,8,4)  = -0.10700e-02_wp
            htt(12,9,4)  = -0.25600e-02_wp
            htt(12,10,4) = 0.34800e-02_wp
            htt(12,11,4) = 0.26800e-02_wp
            htt(12,12,4) = 0.87100e-02_wp
            htt(13,1,4)  = 0.74000e-03_wp
            htt(13,2,4)  = -0.14900e-02_wp
            htt(13,3,4)  = 0.26300e-02_wp
            htt(13,4,4)  = 0.40000e-03_wp
            htt(13,5,4)  = 0.23400e-02_wp
            htt(13,6,4)  = 0.85000e-03_wp
            htt(13,7,4)  = 0.45000e-03_wp
            htt(13,8,4)  = -0.10400e-02_wp
            htt(13,9,4)  = -0.68300e-02_wp
            htt(13,10,4) = -0.31000e-03_wp
            htt(13,11,4) = -0.45500e-02_wp
            htt(13,12,4) = 0.19800e-02_wp
            htt(13,13,4) = 0.96800e-02_wp
            htt(14,1,4)  = 0.19800e-02_wp
            htt(14,2,4)  = 0.49000e-03_wp
            htt(14,3,4)  = 0.15200e-02_wp
            htt(14,4,4)  = 0.30000e-04_wp
            htt(14,5,4)  = 0.11000e-02_wp
            htt(14,6,4)  = -0.11300e-02_wp
            htt(14,7,4)  = 0.15600e-02_wp
            htt(14,8,4)  = -0.67000e-03_wp
            htt(14,9,4)  = -0.45100e-02_wp
            htt(14,10,4) = 0.12300e-02_wp
            htt(14,11,4) = -0.29800e-02_wp
            htt(14,12,4) = 0.27600e-02_wp
            htt(14,13,4) = 0.38800e-02_wp
            htt(14,14,4) = 0.86200e-02_wp
            htt(15,1,4)  = -0.14700e-02_wp
            htt(15,2,4)  = -0.29800e-02_wp
            htt(15,3,4)  = 0.11000e-02_wp
            htt(15,4,4)  = -0.41000e-03_wp
            htt(15,5,4)  = 0.54000e-03_wp
            htt(15,6,4)  = 0.12400e-02_wp
            htt(15,7,4)  = -0.20300e-02_wp
            htt(15,8,4)  = -0.13200e-02_wp
            htt(15,9,4)  = -0.36000e-03_wp
            htt(15,10,4) = -0.56900e-02_wp
            htt(15,11,4) = 0.11900e-02_wp
            htt(15,12,4) = -0.41400e-02_wp
            htt(15,13,4) = 0.30600e-02_wp
            htt(15,14,4) = -0.19400e-02_wp
            htt(15,15,4) = 0.84700e-02_wp
            htt(16,1,4)  = -0.22000e-03_wp
            htt(16,2,4)  = -0.10000e-02_wp
            htt(16,3,4)  = -0.10000e-04_wp
            htt(16,4,4)  = -0.78000e-03_wp
            htt(16,5,4)  = -0.71000e-03_wp
            htt(16,6,4)  = -0.74000e-03_wp
            htt(16,7,4)  = -0.92000e-03_wp
            htt(16,8,4)  = -0.95000e-03_wp
            htt(16,9,4)  = 0.19700e-02_wp
            htt(16,10,4) = -0.41500e-02_wp
            htt(16,11,4) = 0.27600e-02_wp
            htt(16,12,4) = -0.33600e-02_wp
            htt(16,13,4) = -0.27400e-02_wp
            htt(16,14,4) = 0.28000e-02_wp
            htt(16,15,4) = 0.34700e-02_wp
            htt(16,16,4) = 0.90000e-02_wp
            htt(17,1,4)  = 0.23400e-02_wp
            htt(17,2,4)  = 0.85000e-03_wp
            htt(17,3,4)  = 0.45000e-03_wp
            htt(17,4,4)  = -0.10400e-02_wp
            htt(17,5,4)  = 0.74000e-03_wp
            htt(17,6,4)  = -0.14900e-02_wp
            htt(17,7,4)  = 0.26300e-02_wp
            htt(17,8,4)  = 0.40000e-03_wp
            htt(17,9,4)  = 0.56000e-03_wp
            htt(17,10,4) = -0.18400e-02_wp
            htt(17,11,4) = 0.20800e-02_wp
            htt(17,12,4) = -0.32000e-03_wp
            htt(17,13,4) = -0.64700e-02_wp
            htt(17,14,4) = -0.67000e-03_wp
            htt(17,15,4) = -0.40300e-02_wp
            htt(17,16,4) = 0.17700e-02_wp
            htt(17,17,4) = 0.96800e-02_wp
            htt(18,1,4)  = 0.54000e-03_wp
            htt(18,2,4)  = 0.12400e-02_wp
            htt(18,3,4)  = -0.20300e-02_wp
            htt(18,4,4)  = -0.13200e-02_wp
            htt(18,5,4)  = -0.14700e-02_wp
            htt(18,6,4)  = -0.29800e-02_wp
            htt(18,7,4)  = 0.11000e-02_wp
            htt(18,8,4)  = -0.41000e-03_wp
            htt(18,9,4)  = 0.22500e-02_wp
            htt(18,10,4) = 0.51000e-03_wp
            htt(18,11,4) = 0.15300e-02_wp
            htt(18,12,4) = -0.22000e-03_wp
            htt(18,13,4) = -0.40300e-02_wp
            htt(18,14,4) = 0.97000e-03_wp
            htt(18,15,4) = -0.22600e-02_wp
            htt(18,16,4) = 0.27400e-02_wp
            htt(18,17,4) = 0.30600e-02_wp
            htt(18,18,4) = 0.84700e-02_wp
            htt(19,1,4)  = 0.11000e-02_wp
            htt(19,2,4)  = -0.11300e-02_wp
            htt(19,3,4)  = 0.15600e-02_wp
            htt(19,4,4)  = -0.67000e-03_wp
            htt(19,5,4)  = 0.19800e-02_wp
            htt(19,6,4)  = 0.49000e-03_wp
            htt(19,7,4)  = 0.15200e-02_wp
            htt(19,8,4)  = 0.30000e-04_wp
            htt(19,9,4)  = -0.17700e-02_wp
            htt(19,10,4) = -0.33800e-02_wp
            htt(19,11,4) = 0.51000e-03_wp
            htt(19,12,4) = -0.11100e-02_wp
            htt(19,13,4) = -0.67000e-03_wp
            htt(19,14,4) = -0.54000e-02_wp
            htt(19,15,4) = 0.97000e-03_wp
            htt(19,16,4) = -0.37700e-02_wp
            htt(19,17,4) = 0.38800e-02_wp
            htt(19,18,4) = -0.19400e-02_wp
            htt(19,19,4) = 0.86200e-02_wp
            htt(20,1,4)  = -0.71000e-03_wp
            htt(20,2,4)  = -0.74000e-03_wp
            htt(20,3,4)  = -0.92000e-03_wp
            htt(20,4,4)  = -0.95000e-03_wp
            htt(20,5,4)  = -0.22000e-03_wp
            htt(20,6,4)  = -0.10000e-02_wp
            htt(20,7,4)  = -0.10000e-04_wp
            htt(20,8,4)  = -0.78000e-03_wp
            htt(20,9,4)  = -0.70000e-04_wp
            htt(20,10,4) = -0.10300e-02_wp
            htt(20,11,4) = -0.40000e-04_wp
            htt(20,12,4) = -0.10000e-02_wp
            htt(20,13,4) = 0.17700e-02_wp
            htt(20,14,4) = -0.37700e-02_wp
            htt(20,15,4) = 0.27400e-02_wp
            htt(20,16,4) = -0.28000e-02_wp
            htt(20,17,4) = -0.27400e-02_wp
            htt(20,18,4) = 0.34700e-02_wp
            htt(20,19,4) = 0.28000e-02_wp
            htt(20,20,4) = 0.90000e-02_wp
            htt(21,1,4)  = 0.10000e-04_wp
            htt(21,2,4)  = 0.17500e-02_wp
            htt(21,3,4)  = -0.23300e-02_wp
            htt(21,4,4)  = -0.59000e-03_wp
            htt(21,5,4)  = -0.60400e-02_wp
            htt(21,6,4)  = -0.50000e-03_wp
            htt(21,7,4)  = -0.37000e-02_wp
            htt(21,8,4)  = 0.18400e-02_wp
            htt(21,9,4)  = 0.27700e-02_wp
            htt(21,10,4) = 0.11000e-02_wp
            htt(21,11,4) = 0.10000e-02_wp
            htt(21,12,4) = -0.68000e-03_wp
            htt(21,13,4) = 0.56000e-03_wp
            htt(21,14,4) = -0.17700e-02_wp
            htt(21,15,4) = 0.22500e-02_wp
            htt(21,16,4) = -0.70000e-04_wp
            htt(21,17,4) = -0.68300e-02_wp
            htt(21,18,4) = -0.36000e-03_wp
            htt(21,19,4) = -0.45100e-02_wp
            htt(21,20,4) = 0.19700e-02_wp
            htt(21,21,4) = 0.94900e-02_wp
            htt(22,1,4)  = -0.17200e-02_wp
            htt(22,2,4)  = 0.41000e-03_wp
            htt(22,3,4)  = -0.34800e-02_wp
            htt(22,4,4)  = -0.13500e-02_wp
            htt(22,5,4)  = -0.65000e-03_wp
            htt(22,6,4)  = -0.50600e-02_wp
            htt(22,7,4)  = 0.11100e-02_wp
            htt(22,8,4)  = -0.33000e-02_wp
            htt(22,9,4)  = 0.10000e-02_wp
            htt(22,10,4) = 0.15400e-02_wp
            htt(22,11,4) = -0.11800e-02_wp
            htt(22,12,4) = -0.64000e-03_wp
            htt(22,13,4) = 0.20800e-02_wp
            htt(22,14,4) = 0.51000e-03_wp
            htt(22,15,4) = 0.15300e-02_wp
            htt(22,16,4) = -0.40000e-04_wp
            htt(22,17,4) = -0.45500e-02_wp
            htt(22,18,4) = 0.11900e-02_wp
            htt(22,19,4) = -0.29800e-02_wp
            htt(22,20,4) = 0.27600e-02_wp
            htt(22,21,4) = 0.38300e-02_wp
            htt(22,22,4) = 0.83400e-02_wp
            htt(23,1,4)  = 0.17900e-02_wp
            htt(23,2,4)  = 0.13600e-02_wp
            htt(23,3,4)  = 0.12000e-03_wp
            htt(23,4,4)  = -0.31000e-03_wp
            htt(23,5,4)  = -0.38600e-02_wp
            htt(23,6,4)  = 0.97000e-03_wp
            htt(23,7,4)  = -0.21900e-02_wp
            htt(23,8,4)  = 0.26400e-02_wp
            htt(23,9,4)  = 0.11000e-02_wp
            htt(23,10,4) = -0.12200e-02_wp
            htt(23,11,4) = 0.15400e-02_wp
            htt(23,12,4) = -0.78000e-03_wp
            htt(23,13,4) = -0.18400e-02_wp
            htt(23,14,4) = -0.33800e-02_wp
            htt(23,15,4) = 0.51000e-03_wp
            htt(23,16,4) = -0.10300e-02_wp
            htt(23,17,4) = -0.31000e-03_wp
            htt(23,18,4) = -0.56900e-02_wp
            htt(23,19,4) = 0.12300e-02_wp
            htt(23,20,4) = -0.41500e-02_wp
            htt(23,21,4) = 0.31000e-02_wp
            htt(23,22,4) = -0.18300e-02_wp
            htt(23,23,4) = 0.84100e-02_wp
            htt(24,1,4)  = 0.60000e-04_wp
            htt(24,2,4)  = 0.20000e-04_wp
            htt(24,3,4)  = -0.10300e-02_wp
            htt(24,4,4)  = -0.10700e-02_wp
            htt(24,5,4)  = 0.15300e-02_wp
            htt(24,6,4)  = -0.35900e-02_wp
            htt(24,7,4)  = 0.26200e-02_wp
            htt(24,8,4)  = -0.25000e-02_wp
            htt(24,9,4)  = -0.68000e-03_wp
            htt(24,10,4) = -0.78000e-03_wp
            htt(24,11,4) = -0.64000e-03_wp
            htt(24,12,4) = -0.74000e-03_wp
            htt(24,13,4) = -0.32000e-03_wp
            htt(24,14,4) = -0.11100e-02_wp
            htt(24,15,4) = -0.22000e-03_wp
            htt(24,16,4) = -0.10000e-02_wp
            htt(24,17,4) = 0.19800e-02_wp
            htt(24,18,4) = -0.41400e-02_wp
            htt(24,19,4) = 0.27600e-02_wp
            htt(24,20,4) = -0.33600e-02_wp
            htt(24,21,4) = -0.25600e-02_wp
            htt(24,22,4) = 0.26800e-02_wp
            htt(24,23,4) = 0.34800e-02_wp
            htt(24,24,4) = 0.87100e-02_wp
            ! Hessian elements for anchor point 4
            hbs( 1,1,4)  = -0.31760e-01_wp
            hbs( 1,2,4)  = -0.31760e-01_wp
            hbs( 1,3,4)  = -0.37700e-01_wp
            hbs( 1,4,4)  = -0.11100e-01_wp
            hbs( 1,5,4)  = 0.40600e-02_wp
            hbs( 1,6,4)  = 0.43670e-01_wp
            hbs( 1,7,4)  = -0.35000e-03_wp
            hbs( 1,8,4)  = 0.43670e-01_wp
            hbs( 1,9,4)  = -0.10000e-04_wp
            hbs( 1,10,4) = -0.11100e-01_wp
            hbs( 1,11,4) = -0.35000e-03_wp
            hbs( 1,12,4) = 0.40600e-02_wp
            hbs( 2,1,4)  = 0.39800e-01_wp
            hbs( 2,2,4)  = -0.80400e-02_wp
            hbs( 2,3,4)  = 0.18850e-01_wp
            hbs( 2,4,4)  = -0.12000e-03_wp
            hbs( 2,5,4)  = -0.21800e-02_wp
            hbs( 2,6,4)  = -0.18210e-01_wp
            hbs( 2,7,4)  = 0.75000e-03_wp
            hbs( 2,8,4)  = -0.25460e-01_wp
            hbs( 2,9,4)  = 0.00000e+00_wp
            hbs( 2,10,4) = 0.11220e-01_wp
            hbs( 2,11,4) = -0.39000e-03_wp
            hbs( 2,12,4) = -0.18800e-02_wp
            hbs( 3,1,4)  = -0.80400e-02_wp
            hbs( 3,2,4)  = 0.39800e-01_wp
            hbs( 3,3,4)  = 0.18850e-01_wp
            hbs( 3,4,4)  = 0.11220e-01_wp
            hbs( 3,5,4)  = -0.18800e-02_wp
            hbs( 3,6,4)  = -0.25460e-01_wp
            hbs( 3,7,4)  = -0.39000e-03_wp
            hbs( 3,8,4)  = -0.18210e-01_wp
            hbs( 3,9,4)  = 0.00000e+00_wp
            hbs( 3,10,4) = -0.12000e-03_wp
            hbs( 3,11,4) = 0.75000e-03_wp
            hbs( 3,12,4) = -0.21800e-02_wp
            hbs( 4,1,4)  = -0.26690e-01_wp
            hbs( 4,2,4)  = -0.99100e-02_wp
            hbs( 4,3,4)  = 0.10200e-01_wp
            hbs( 4,4,4)  = -0.32600e-01_wp
            hbs( 4,5,4)  = -0.82100e-02_wp
            hbs( 4,6,4)  = -0.12130e-01_wp
            hbs( 4,7,4)  = 0.47800e-02_wp
            hbs( 4,8,4)  = 0.46360e-01_wp
            hbs( 4,9,4)  = -0.61000e-03_wp
            hbs( 4,10,4) = 0.44260e-01_wp
            hbs( 4,11,4) = 0.10000e-03_wp
            hbs( 4,12,4) = -0.43000e-03_wp
            hbs( 5,1,4)  = 0.27720e-01_wp
            hbs( 5,2,4)  = 0.58200e-02_wp
            hbs( 5,3,4)  = -0.92100e-02_wp
            hbs( 5,4,4)  = 0.33000e-02_wp
            hbs( 5,5,4)  = 0.59500e-02_wp
            hbs( 5,6,4)  = 0.56400e-02_wp
            hbs( 5,7,4)  = -0.18300e-02_wp
            hbs( 5,8,4)  = -0.20790e-01_wp
            hbs( 5,9,4)  = -0.31000e-03_wp
            hbs( 5,10,4) = -0.24570e-01_wp
            hbs( 5,11,4) = -0.60000e-04_wp
            hbs( 5,12,4) = 0.64000e-03_wp
            hbs( 6,1,4)  = -0.10300e-02_wp
            hbs( 6,2,4)  = 0.40900e-02_wp
            hbs( 6,3,4)  = -0.99000e-03_wp
            hbs( 6,4,4)  = 0.29300e-01_wp
            hbs( 6,5,4)  = 0.22600e-02_wp
            hbs( 6,6,4)  = 0.65000e-02_wp
            hbs( 6,7,4)  = -0.29500e-02_wp
            hbs( 6,8,4)  = -0.25570e-01_wp
            hbs( 6,9,4)  = 0.92000e-03_wp
            hbs( 6,10,4) = -0.19690e-01_wp
            hbs( 6,11,4) = -0.40000e-04_wp
            hbs( 6,12,4) = -0.21000e-03_wp
            hbs( 7,1,4)  = -0.14540e-01_wp
            hbs( 7,2,4)  = 0.40000e-01_wp
            hbs( 7,3,4)  = 0.98300e-02_wp
            hbs( 7,4,4)  = -0.32290e-01_wp
            hbs( 7,5,4)  = 0.50700e-02_wp
            hbs( 7,6,4)  = -0.32910e-01_wp
            hbs( 7,7,4)  = -0.88600e-02_wp
            hbs( 7,8,4)  = -0.12410e-01_wp
            hbs( 7,9,4)  = 0.48800e-02_wp
            hbs( 7,10,4) = 0.45460e-01_wp
            hbs( 7,11,4) = -0.87000e-03_wp
            hbs( 7,12,4) = 0.14000e-03_wp
            hbs( 8,1,4)  = 0.89100e-02_wp
            hbs( 8,2,4)  = -0.22340e-01_wp
            hbs( 8,3,4)  = -0.58300e-02_wp
            hbs( 8,4,4)  = 0.30270e-01_wp
            hbs( 8,5,4)  = -0.31400e-02_wp
            hbs( 8,6,4)  = 0.36400e-02_wp
            hbs( 8,7,4)  = 0.47300e-02_wp
            hbs( 8,8,4)  = 0.54700e-02_wp
            hbs( 8,9,4)  = -0.19600e-02_wp
            hbs( 8,10,4) = -0.20340e-01_wp
            hbs( 8,11,4) = -0.23000e-03_wp
            hbs( 8,12,4) = -0.12000e-03_wp
            hbs( 9,1,4)  = 0.56200e-02_wp
            hbs( 9,2,4)  = -0.17660e-01_wp
            hbs( 9,3,4)  = -0.39900e-02_wp
            hbs( 9,4,4)  = 0.20200e-02_wp
            hbs( 9,5,4)  = -0.19200e-02_wp
            hbs( 9,6,4)  = 0.29280e-01_wp
            hbs( 9,7,4)  = 0.41300e-02_wp
            hbs( 9,8,4)  = 0.69400e-02_wp
            hbs( 9,9,4)  = -0.29200e-02_wp
            hbs( 9,10,4) = -0.25120e-01_wp
            hbs( 9,11,4) = 0.11000e-02_wp
            hbs( 9,12,4) = -0.30000e-04_wp
            hbs(10,1,4)  = 0.42900e-01_wp
            hbs(10,2,4)  = 0.42900e-01_wp
            hbs(10,3,4)  = -0.23400e-02_wp
            hbs(10,4,4)  = -0.13730e-01_wp
            hbs(10,5,4)  = -0.64000e-03_wp
            hbs(10,6,4)  = -0.32580e-01_wp
            hbs(10,7,4)  = 0.52100e-02_wp
            hbs(10,8,4)  = -0.32580e-01_wp
            hbs(10,9,4)  = -0.85400e-02_wp
            hbs(10,10,4) = -0.13730e-01_wp
            hbs(10,11,4) = 0.52100e-02_wp
            hbs(10,12,4) = -0.64000e-03_wp
            hbs(11,1,4)  = -0.23760e-01_wp
            hbs(11,2,4)  = -0.19130e-01_wp
            hbs(11,3,4)  = 0.11700e-02_wp
            hbs(11,4,4)  = 0.73400e-02_wp
            hbs(11,5,4)  = 0.91000e-03_wp
            hbs(11,6,4)  = 0.30250e-01_wp
            hbs(11,7,4)  = -0.31400e-02_wp
            hbs(11,8,4)  = 0.23200e-02_wp
            hbs(11,9,4)  = 0.42700e-02_wp
            hbs(11,10,4) = 0.64000e-02_wp
            hbs(11,11,4) = -0.20700e-02_wp
            hbs(11,12,4) = -0.27000e-03_wp
            hbs(12,1,4)  = -0.19130e-01_wp
            hbs(12,2,4)  = -0.23760e-01_wp
            hbs(12,3,4)  = 0.11700e-02_wp
            hbs(12,4,4)  = 0.64000e-02_wp
            hbs(12,5,4)  = -0.27000e-03_wp
            hbs(12,6,4)  = 0.23200e-02_wp
            hbs(12,7,4)  = -0.20700e-02_wp
            hbs(12,8,4)  = 0.30250e-01_wp
            hbs(12,9,4)  = 0.42700e-02_wp
            hbs(12,10,4) = 0.73400e-02_wp
            hbs(12,11,4) = -0.31400e-02_wp
            hbs(12,12,4) = 0.91000e-03_wp
            hbs(13,1,4)  = 0.40000e-01_wp
            hbs(13,2,4)  = -0.14540e-01_wp
            hbs(13,3,4)  = 0.98300e-02_wp
            hbs(13,4,4)  = 0.45460e-01_wp
            hbs(13,5,4)  = 0.14000e-03_wp
            hbs(13,6,4)  = -0.12410e-01_wp
            hbs(13,7,4)  = -0.87000e-03_wp
            hbs(13,8,4)  = -0.32910e-01_wp
            hbs(13,9,4)  = 0.48800e-02_wp
            hbs(13,10,4) = -0.32290e-01_wp
            hbs(13,11,4) = -0.88600e-02_wp
            hbs(13,12,4) = 0.50700e-02_wp
            hbs(14,1,4)  = -0.17660e-01_wp
            hbs(14,2,4)  = 0.56200e-02_wp
            hbs(14,3,4)  = -0.39900e-02_wp
            hbs(14,4,4)  = -0.25120e-01_wp
            hbs(14,5,4)  = -0.30000e-04_wp
            hbs(14,6,4)  = 0.69400e-02_wp
            hbs(14,7,4)  = 0.11000e-02_wp
            hbs(14,8,4)  = 0.29280e-01_wp
            hbs(14,9,4)  = -0.29200e-02_wp
            hbs(14,10,4) = 0.20200e-02_wp
            hbs(14,11,4) = 0.41300e-02_wp
            hbs(14,12,4) = -0.19200e-02_wp
            hbs(15,1,4)  = -0.22340e-01_wp
            hbs(15,2,4)  = 0.89100e-02_wp
            hbs(15,3,4)  = -0.58300e-02_wp
            hbs(15,4,4)  = -0.20340e-01_wp
            hbs(15,5,4)  = -0.12000e-03_wp
            hbs(15,6,4)  = 0.54700e-02_wp
            hbs(15,7,4)  = -0.23000e-03_wp
            hbs(15,8,4)  = 0.36400e-02_wp
            hbs(15,9,4)  = -0.19600e-02_wp
            hbs(15,10,4) = 0.30270e-01_wp
            hbs(15,11,4) = 0.47300e-02_wp
            hbs(15,12,4) = -0.31400e-02_wp
            hbs(16,1,4)  = -0.99100e-02_wp
            hbs(16,2,4)  = -0.26690e-01_wp
            hbs(16,3,4)  = 0.10200e-01_wp
            hbs(16,4,4)  = 0.44260e-01_wp
            hbs(16,5,4)  = -0.43000e-03_wp
            hbs(16,6,4)  = 0.46360e-01_wp
            hbs(16,7,4)  = 0.10000e-03_wp
            hbs(16,8,4)  = -0.12130e-01_wp
            hbs(16,9,4)  = -0.61000e-03_wp
            hbs(16,10,4) = -0.32600e-01_wp
            hbs(16,11,4) = 0.47800e-02_wp
            hbs(16,12,4) = -0.82100e-02_wp
            hbs(17,1,4)  = 0.58200e-02_wp
            hbs(17,2,4)  = 0.27720e-01_wp
            hbs(17,3,4)  = -0.92100e-02_wp
            hbs(17,4,4)  = -0.24570e-01_wp
            hbs(17,5,4)  = 0.64000e-03_wp
            hbs(17,6,4)  = -0.20790e-01_wp
            hbs(17,7,4)  = -0.60000e-04_wp
            hbs(17,8,4)  = 0.56400e-02_wp
            hbs(17,9,4)  = -0.31000e-03_wp
            hbs(17,10,4) = 0.33000e-02_wp
            hbs(17,11,4) = -0.18300e-02_wp
            hbs(17,12,4) = 0.59500e-02_wp
            hbs(18,1,4)  = 0.40900e-02_wp
            hbs(18,2,4)  = -0.10300e-02_wp
            hbs(18,3,4)  = -0.99000e-03_wp
            hbs(18,4,4)  = -0.19690e-01_wp
            hbs(18,5,4)  = -0.21000e-03_wp
            hbs(18,6,4)  = -0.25570e-01_wp
            hbs(18,7,4)  = -0.40000e-04_wp
            hbs(18,8,4)  = 0.65000e-02_wp
            hbs(18,9,4)  = 0.92000e-03_wp
            hbs(18,10,4) = 0.29300e-01_wp
            hbs(18,11,4) = -0.29500e-02_wp
            hbs(18,12,4) = 0.22600e-02_wp
        elseif( ii .eq. 22 )then
            ! Hessian elements for anchor point 1
            hss( 1,1,1)  = 0.35006e+00_wp
            hss( 2,1,1)  = 0.13620e-01_wp
            hss( 2,2,1)  = 0.36874e+00_wp
            hss( 3,1,1)  = 0.29130e-01_wp
            hss( 3,2,1)  = 0.28430e-01_wp
            hss( 3,3,1)  = 0.46705e+00_wp
            hss( 4,1,1)  = -0.29900e-02_wp
            hss( 4,2,1)  = 0.16080e-01_wp
            hss( 4,3,1)  = -0.77700e-02_wp
            hss( 4,4,1)  = 0.38651e+00_wp
            hss( 5,1,1)  = 0.40700e-02_wp
            hss( 5,2,1)  = -0.10000e-02_wp
            hss( 5,3,1)  = 0.38700e-02_wp
            hss( 5,4,1)  = 0.32400e-02_wp
            hss( 5,5,1)  = 0.36515e+00_wp
            hss( 6,1,1)  = 0.33920e-01_wp
            hss( 6,2,1)  = 0.46890e-01_wp
            hss( 6,3,1)  = -0.39100e-02_wp
            hss( 6,4,1)  = -0.52200e-02_wp
            hss( 6,5,1)  = 0.40000e-03_wp
            hss( 6,6,1)  = 0.36166e+00_wp
            hss( 7,1,1)  = -0.40000e-03_wp
            hss( 7,2,1)  = -0.70000e-03_wp
            hss( 7,3,1)  = 0.45000e-03_wp
            hss( 7,4,1)  = 0.39500e-02_wp
            hss( 7,5,1)  = 0.42000e-03_wp
            hss( 7,6,1)  = 0.39100e-02_wp
            hss( 7,7,1)  = 0.37071e+00_wp
            hss( 8,1,1)  = 0.44240e-01_wp
            hss( 8,2,1)  = 0.34840e-01_wp
            hss( 8,3,1)  = -0.58500e-02_wp
            hss( 8,4,1)  = 0.24280e-01_wp
            hss( 8,5,1)  = -0.11400e-02_wp
            hss( 8,6,1)  = 0.14590e-01_wp
            hss( 8,7,1)  = -0.98000e-03_wp
            hss( 8,8,1)  = 0.36092e+00_wp
            hss( 9,1,1)  = -0.57000e-03_wp
            hss( 9,2,1)  = -0.54000e-03_wp
            hss( 9,3,1)  = 0.41000e-03_wp
            hss( 9,4,1)  = -0.31000e-03_wp
            hss( 9,5,1)  = 0.17000e-03_wp
            hss( 9,6,1)  = 0.39200e-02_wp
            hss( 9,7,1)  = 0.50000e-03_wp
            hss( 9,8,1)  = 0.39400e-02_wp
            hss( 9,9,1)  = 0.36389e+00_wp
            hss(10,1,1)  = 0.19730e-01_wp
            hss(10,2,1)  = -0.83300e-02_wp
            hss(10,3,1)  = -0.46800e-02_wp
            hss(10,4,1)  = 0.63810e-01_wp
            hss(10,5,1)  = -0.33000e-03_wp
            hss(10,6,1)  = 0.22020e-01_wp
            hss(10,7,1)  = -0.10900e-02_wp
            hss(10,8,1)  = -0.21000e-02_wp
            hss(10,9,1)  = -0.19000e-03_wp
            hss(10,10,1) = 0.38332e+00_wp
            hss(11,1,1)  = -0.70000e-03_wp
            hss(11,2,1)  = -0.43000e-03_wp
            hss(11,3,1)  = 0.55000e-03_wp
            hss(11,4,1)  = -0.11400e-02_wp
            hss(11,5,1)  = 0.40000e-04_wp
            hss(11,6,1)  = -0.97000e-03_wp
            hss(11,7,1)  = 0.22000e-03_wp
            hss(11,8,1)  = 0.38800e-02_wp
            hss(11,9,1)  = 0.48000e-03_wp
            hss(11,10,1) = 0.39000e-02_wp
            hss(11,11,1) = 0.37178e+00_wp
            hss(12,1,1)  = -0.14000e-02_wp
            hss(12,2,1)  = 0.34100e-02_wp
            hss(12,3,1)  = 0.37000e-03_wp
            hss(12,4,1)  = -0.13000e-03_wp
            hss(12,5,1)  = 0.24000e-03_wp
            hss(12,6,1)  = -0.10500e-02_wp
            hss(12,7,1)  = 0.20000e-04_wp
            hss(12,8,1)  = 0.55000e-03_wp
            hss(12,9,1)  = 0.12000e-03_wp
            hss(12,10,1) = 0.24400e-02_wp
            hss(12,11,1) = 0.39000e-03_wp
            hss(12,12,1) = 0.37251e+00_wp
            ! Hessian elements for anchor point 1
            hbb( 1,1,1)  = 0.82220e-01_wp
            hbb( 2,1,1)  = -0.42040e-01_wp
            hbb( 2,2,1)  = 0.14128e+00_wp
            hbb( 3,1,1)  = -0.40180e-01_wp
            hbb( 3,2,1)  = -0.99240e-01_wp
            hbb( 3,3,1)  = 0.13942e+00_wp
            hbb( 4,1,1)  = -0.36210e-01_wp
            hbb( 4,2,1)  = 0.22250e-01_wp
            hbb( 4,3,1)  = 0.13960e-01_wp
            hbb( 4,4,1)  = 0.66190e-01_wp
            hbb( 5,1,1)  = 0.21220e-01_wp
            hbb( 5,2,1)  = -0.13320e-01_wp
            hbb( 5,3,1)  = -0.79000e-02_wp
            hbb( 5,4,1)  = -0.33260e-01_wp
            hbb( 5,5,1)  = 0.72740e-01_wp
            hbb( 6,1,1)  = 0.14990e-01_wp
            hbb( 6,2,1)  = -0.89300e-02_wp
            hbb( 6,3,1)  = -0.60600e-02_wp
            hbb( 6,4,1)  = -0.32930e-01_wp
            hbb( 6,5,1)  = -0.39480e-01_wp
            hbb( 6,6,1)  = 0.72410e-01_wp
            hbb( 7,1,1)  = -0.67600e-02_wp
            hbb( 7,2,1)  = -0.46000e-02_wp
            hbb( 7,3,1)  = 0.11370e-01_wp
            hbb( 7,4,1)  = -0.26080e-01_wp
            hbb( 7,5,1)  = 0.93900e-02_wp
            hbb( 7,6,1)  = 0.16690e-01_wp
            hbb( 7,7,1)  = 0.64820e-01_wp
            hbb( 8,1,1)  = -0.24700e-02_wp
            hbb( 8,2,1)  = 0.77600e-02_wp
            hbb( 8,3,1)  = -0.52900e-02_wp
            hbb( 8,4,1)  = 0.15900e-01_wp
            hbb( 8,5,1)  = -0.49900e-02_wp
            hbb( 8,6,1)  = -0.10910e-01_wp
            hbb( 8,7,1)  = -0.31790e-01_wp
            hbb( 8,8,1)  = 0.73440e-01_wp
            hbb( 9,1,1)  = 0.92300e-02_wp
            hbb( 9,2,1)  = -0.31500e-02_wp
            hbb( 9,3,1)  = -0.60800e-02_wp
            hbb( 9,4,1)  = 0.10190e-01_wp
            hbb( 9,5,1)  = -0.44000e-02_wp
            hbb( 9,6,1)  = -0.57900e-02_wp
            hbb( 9,7,1)  = -0.33030e-01_wp
            hbb( 9,8,1)  = -0.41650e-01_wp
            hbb( 9,9,1)  = 0.74680e-01_wp
            hbb(10,1,1)  = 0.38400e-02_wp
            hbb(10,2,1)  = -0.14400e-02_wp
            hbb(10,3,1)  = -0.24000e-02_wp
            hbb(10,4,1)  = -0.57300e-02_wp
            hbb(10,5,1)  = 0.88600e-02_wp
            hbb(10,6,1)  = -0.31300e-02_wp
            hbb(10,7,1)  = -0.33500e-01_wp
            hbb(10,8,1)  = 0.13490e-01_wp
            hbb(10,9,1)  = 0.20010e-01_wp
            hbb(10,10,1) = 0.75070e-01_wp
            hbb(11,1,1)  = -0.18400e-02_wp
            hbb(11,2,1)  = 0.11400e-02_wp
            hbb(11,3,1)  = 0.70000e-03_wp
            hbb(11,4,1)  = -0.27200e-02_wp
            hbb(11,5,1)  = -0.30000e-02_wp
            hbb(11,6,1)  = 0.57200e-02_wp
            hbb(11,7,1)  = 0.20110e-01_wp
            hbb(11,8,1)  = -0.74900e-02_wp
            hbb(11,9,1)  = -0.12620e-01_wp
            hbb(11,10,1) = -0.37620e-01_wp
            hbb(11,11,1) = 0.78370e-01_wp
            hbb(12,1,1)  = -0.20000e-02_wp
            hbb(12,2,1)  = 0.30000e-03_wp
            hbb(12,3,1)  = 0.17000e-02_wp
            hbb(12,4,1)  = 0.84500e-02_wp
            hbb(12,5,1)  = -0.58700e-02_wp
            hbb(12,6,1)  = -0.25800e-02_wp
            hbb(12,7,1)  = 0.13390e-01_wp
            hbb(12,8,1)  = -0.60000e-02_wp
            hbb(12,9,1)  = -0.73900e-02_wp
            hbb(12,10,1) = -0.37460e-01_wp
            hbb(12,11,1) = -0.40760e-01_wp
            hbb(12,12,1) = 0.78220e-01_wp
            hbb(13,1,1)  = -0.61900e-02_wp
            hbb(13,2,1)  = 0.10630e-01_wp
            hbb(13,3,1)  = -0.44400e-02_wp
            hbb(13,4,1)  = -0.14000e-02_wp
            hbb(13,5,1)  = 0.11900e-02_wp
            hbb(13,6,1)  = 0.21000e-03_wp
            hbb(13,7,1)  = 0.29200e-02_wp
            hbb(13,8,1)  = 0.41800e-02_wp
            hbb(13,9,1)  = -0.71000e-02_wp
            hbb(13,10,1) = -0.33570e-01_wp
            hbb(13,11,1) = 0.13390e-01_wp
            hbb(13,12,1) = 0.20180e-01_wp
            hbb(13,13,1) = 0.65200e-01_wp
            hbb(14,1,1)  = 0.89900e-02_wp
            hbb(14,2,1)  = -0.58800e-02_wp
            hbb(14,3,1)  = -0.31100e-02_wp
            hbb(14,4,1)  = 0.77000e-03_wp
            hbb(14,5,1)  = -0.83000e-03_wp
            hbb(14,6,1)  = 0.70000e-04_wp
            hbb(14,7,1)  = -0.71700e-02_wp
            hbb(14,8,1)  = -0.51000e-03_wp
            hbb(14,9,1)  = 0.76800e-02_wp
            hbb(14,10,1) = 0.20000e-01_wp
            hbb(14,11,1) = -0.73300e-02_wp
            hbb(14,12,1) = -0.12670e-01_wp
            hbb(14,13,1) = -0.33080e-01_wp
            hbb(14,14,1) = 0.74340e-01_wp
            hbb(15,1,1)  = -0.28000e-02_wp
            hbb(15,2,1)  = -0.47500e-02_wp
            hbb(15,3,1)  = 0.75500e-02_wp
            hbb(15,4,1)  = 0.64000e-03_wp
            hbb(15,5,1)  = -0.36000e-03_wp
            hbb(15,6,1)  = -0.28000e-03_wp
            hbb(15,7,1)  = 0.42500e-02_wp
            hbb(15,8,1)  = -0.36700e-02_wp
            hbb(15,9,1)  = -0.58000e-03_wp
            hbb(15,10,1) = 0.13570e-01_wp
            hbb(15,11,1) = -0.60500e-02_wp
            hbb(15,12,1) = -0.75200e-02_wp
            hbb(15,13,1) = -0.32120e-01_wp
            hbb(15,14,1) = -0.41260e-01_wp
            hbb(15,15,1) = 0.73390e-01_wp
            hbb(16,1,1)  = -0.36890e-01_wp
            hbb(16,2,1)  = 0.15200e-01_wp
            hbb(16,3,1)  = 0.21690e-01_wp
            hbb(16,4,1)  = 0.32300e-02_wp
            hbb(16,5,1)  = -0.74000e-02_wp
            hbb(16,6,1)  = 0.41700e-02_wp
            hbb(16,7,1)  = -0.13900e-02_wp
            hbb(16,8,1)  = 0.69000e-03_wp
            hbb(16,9,1)  = 0.71000e-03_wp
            hbb(16,10,1) = -0.61100e-02_wp
            hbb(16,11,1) = 0.86800e-02_wp
            hbb(16,12,1) = -0.25700e-02_wp
            hbb(16,13,1) = -0.26960e-01_wp
            hbb(16,14,1) = 0.10490e-01_wp
            hbb(16,15,1) = 0.16470e-01_wp
            hbb(16,16,1) = 0.68110e-01_wp
            hbb(17,1,1)  = 0.22720e-01_wp
            hbb(17,2,1)  = -0.69100e-02_wp
            hbb(17,3,1)  = -0.15810e-01_wp
            hbb(17,4,1)  = -0.78500e-02_wp
            hbb(17,5,1)  = 0.81700e-02_wp
            hbb(17,6,1)  = -0.32000e-03_wp
            hbb(17,7,1)  = 0.10300e-02_wp
            hbb(17,8,1)  = -0.32000e-03_wp
            hbb(17,9,1)  = -0.71000e-03_wp
            hbb(17,10,1) = 0.90700e-02_wp
            hbb(17,11,1) = -0.60300e-02_wp
            hbb(17,12,1) = -0.30400e-02_wp
            hbb(17,13,1) = 0.97900e-02_wp
            hbb(17,14,1) = -0.45700e-02_wp
            hbb(17,15,1) = -0.52200e-02_wp
            hbb(17,16,1) = -0.34760e-01_wp
            hbb(17,17,1) = 0.72170e-01_wp
            hbb(18,1,1)  = 0.14170e-01_wp
            hbb(18,2,1)  = -0.82800e-02_wp
            hbb(18,3,1)  = -0.58800e-02_wp
            hbb(18,4,1)  = 0.46200e-02_wp
            hbb(18,5,1)  = -0.77000e-03_wp
            hbb(18,6,1)  = -0.38500e-02_wp
            hbb(18,7,1)  = 0.36000e-03_wp
            hbb(18,8,1)  = -0.37000e-03_wp
            hbb(18,9,1)  = 0.10000e-04_wp
            hbb(18,10,1) = -0.29600e-02_wp
            hbb(18,11,1) = -0.26500e-02_wp
            hbb(18,12,1) = 0.56000e-02_wp
            hbb(18,13,1) = 0.17170e-01_wp
            hbb(18,14,1) = -0.59200e-02_wp
            hbb(18,15,1) = -0.11250e-01_wp
            hbb(18,16,1) = -0.33360e-01_wp
            hbb(18,17,1) = -0.37420e-01_wp
            hbb(18,18,1) = 0.70770e-01_wp
            ! Hessian elements for anchor point 1
            htt( 1,1,1)  = 0.59700e-02_wp
            htt( 2,1,1)  = 0.40300e-02_wp
            htt( 2,2,1)  = 0.53700e-02_wp
            htt( 3,1,1)  = -0.14500e-02_wp
            htt( 3,2,1)  = -0.35400e-02_wp
            htt( 3,3,1)  = 0.68100e-02_wp
            htt( 4,1,1)  = -0.33900e-02_wp
            htt( 4,2,1)  = -0.21900e-02_wp
            htt( 4,3,1)  = 0.47200e-02_wp
            htt( 4,4,1)  = 0.59200e-02_wp
            htt( 5,1,1)  = -0.57900e-02_wp
            htt( 5,2,1)  = -0.45300e-02_wp
            htt( 5,3,1)  = 0.17600e-02_wp
            htt( 5,4,1)  = 0.30200e-02_wp
            htt( 5,5,1)  = 0.64900e-02_wp
            htt( 6,1,1)  = -0.38500e-02_wp
            htt( 6,2,1)  = -0.41000e-02_wp
            htt( 6,3,1)  = 0.21800e-02_wp
            htt( 6,4,1)  = 0.19300e-02_wp
            htt( 6,5,1)  = 0.34500e-02_wp
            htt( 6,6,1)  = 0.52100e-02_wp
            htt( 7,1,1)  = 0.13600e-02_wp
            htt( 7,2,1)  = 0.27600e-02_wp
            htt( 7,3,1)  = -0.62000e-02_wp
            htt( 7,4,1)  = -0.48000e-02_wp
            htt( 7,5,1)  = -0.78000e-03_wp
            htt( 7,6,1)  = -0.23700e-02_wp
            htt( 7,7,1)  = 0.65100e-02_wp
            htt( 8,1,1)  = 0.32900e-02_wp
            htt( 8,2,1)  = 0.31900e-02_wp
            htt( 8,3,1)  = -0.57800e-02_wp
            htt( 8,4,1)  = -0.58800e-02_wp
            htt( 8,5,1)  = -0.38200e-02_wp
            htt( 8,6,1)  = -0.61000e-03_wp
            htt( 8,7,1)  = 0.49200e-02_wp
            htt( 8,8,1)  = 0.81300e-02_wp
            htt( 9,1,1)  = -0.23500e-02_wp
            htt( 9,2,1)  = -0.73000e-03_wp
            htt( 9,3,1)  = 0.17000e-03_wp
            htt( 9,4,1)  = 0.17900e-02_wp
            htt( 9,5,1)  = 0.15900e-02_wp
            htt( 9,6,1)  = 0.12400e-02_wp
            htt( 9,7,1)  = -0.84000e-03_wp
            htt( 9,8,1)  = -0.11900e-02_wp
            htt( 9,9,1)  = 0.25000e-02_wp
            htt(10,1,1)  = -0.13700e-02_wp
            htt(10,2,1)  = 0.15000e-03_wp
            htt(10,3,1)  = -0.60000e-04_wp
            htt(10,4,1)  = 0.14600e-02_wp
            htt(10,5,1)  = 0.15100e-02_wp
            htt(10,6,1)  = 0.46000e-03_wp
            htt(10,7,1)  = 0.25000e-03_wp
            htt(10,8,1)  = -0.80000e-03_wp
            htt(10,9,1)  = 0.46000e-03_wp
            htt(10,10,1) = 0.37200e-02_wp
            htt(11,1,1)  = -0.36000e-03_wp
            htt(11,2,1)  = -0.21100e-02_wp
            htt(11,3,1)  = 0.23100e-02_wp
            htt(11,4,1)  = 0.56000e-03_wp
            htt(11,5,1)  = 0.29000e-03_wp
            htt(11,6,1)  = 0.15000e-02_wp
            htt(11,7,1)  = -0.22800e-02_wp
            htt(11,8,1)  = -0.10800e-02_wp
            htt(11,9,1)  = 0.83000e-03_wp
            htt(11,10,1) = -0.11000e-02_wp
            htt(11,11,1) = 0.26300e-02_wp
            htt(12,1,1)  = 0.62000e-03_wp
            htt(12,2,1)  = -0.12400e-02_wp
            htt(12,3,1)  = 0.20800e-02_wp
            htt(12,4,1)  = 0.22000e-03_wp
            htt(12,5,1)  = 0.21000e-03_wp
            htt(12,6,1)  = 0.72000e-03_wp
            htt(12,7,1)  = -0.11900e-02_wp
            htt(12,8,1)  = -0.69000e-03_wp
            htt(12,9,1)  = -0.12000e-02_wp
            htt(12,10,1) = 0.21600e-02_wp
            htt(12,11,1) = 0.70000e-03_wp
            htt(12,12,1) = 0.40700e-02_wp
            htt(13,1,1)  = -0.11400e-02_wp
            htt(13,2,1)  = -0.19000e-02_wp
            htt(13,3,1)  = 0.76000e-03_wp
            htt(13,4,1)  = -0.10000e-04_wp
            htt(13,5,1)  = 0.16400e-02_wp
            htt(13,6,1)  = 0.16200e-02_wp
            htt(13,7,1)  = -0.19000e-03_wp
            htt(13,8,1)  = -0.20000e-03_wp
            htt(13,9,1)  = -0.21000e-02_wp
            htt(13,10,1) = 0.24000e-03_wp
            htt(13,11,1) = -0.13100e-02_wp
            htt(13,12,1) = 0.10300e-02_wp
            htt(13,13,1) = 0.49700e-02_wp
            htt(14,1,1)  = 0.70000e-03_wp
            htt(14,2,1)  = 0.47000e-03_wp
            htt(14,3,1)  = -0.32000e-03_wp
            htt(14,4,1)  = -0.55000e-03_wp
            htt(14,5,1)  = -0.18000e-03_wp
            htt(14,6,1)  = -0.82000e-03_wp
            htt(14,7,1)  = 0.80000e-03_wp
            htt(14,8,1)  = 0.16000e-03_wp
            htt(14,9,1)  = -0.22000e-03_wp
            htt(14,10,1) = 0.13600e-02_wp
            htt(14,11,1) = 0.10000e-04_wp
            htt(14,12,1) = 0.15900e-02_wp
            htt(14,13,1) = -0.78000e-03_wp
            htt(14,14,1) = 0.38600e-02_wp
            htt(15,1,1)  = -0.21100e-02_wp
            htt(15,2,1)  = -0.27800e-02_wp
            htt(15,3,1)  = 0.98000e-03_wp
            htt(15,4,1)  = 0.32000e-03_wp
            htt(15,5,1)  = 0.17100e-02_wp
            htt(15,6,1)  = 0.23900e-02_wp
            htt(15,7,1)  = -0.12700e-02_wp
            htt(15,8,1)  = -0.59000e-03_wp
            htt(15,9,1)  = -0.80000e-04_wp
            htt(15,10,1) = -0.29900e-02_wp
            htt(15,11,1) = 0.61000e-03_wp
            htt(15,12,1) = -0.23100e-02_wp
            htt(15,13,1) = 0.26500e-02_wp
            htt(15,14,1) = -0.23400e-02_wp
            htt(15,15,1) = 0.55400e-02_wp
            htt(16,1,1)  = -0.27000e-03_wp
            htt(16,2,1)  = -0.40000e-03_wp
            htt(16,3,1)  = -0.90000e-04_wp
            htt(16,4,1)  = -0.22000e-03_wp
            htt(16,5,1)  = -0.10000e-03_wp
            htt(16,6,1)  = -0.50000e-04_wp
            htt(16,7,1)  = -0.28000e-03_wp
            htt(16,8,1)  = -0.22000e-03_wp
            htt(16,9,1)  = 0.18000e-02_wp
            htt(16,10,1) = -0.18700e-02_wp
            htt(16,11,1) = 0.19300e-02_wp
            htt(16,12,1) = -0.17400e-02_wp
            htt(16,13,1) = -0.30900e-02_wp
            htt(16,14,1) = 0.23000e-02_wp
            htt(16,15,1) = 0.55000e-03_wp
            htt(16,16,1) = 0.59400e-02_wp
            htt(17,1,1)  = 0.13400e-02_wp
            htt(17,2,1)  = 0.14100e-02_wp
            htt(17,3,1)  = -0.45000e-03_wp
            htt(17,4,1)  = -0.38000e-03_wp
            htt(17,5,1)  = -0.94000e-03_wp
            htt(17,6,1)  = -0.20500e-02_wp
            htt(17,7,1)  = 0.78000e-03_wp
            htt(17,8,1)  = -0.32000e-03_wp
            htt(17,9,1)  = 0.13000e-02_wp
            htt(17,10,1) = -0.10000e-03_wp
            htt(17,11,1) = 0.12300e-02_wp
            htt(17,12,1) = -0.17000e-03_wp
            htt(17,13,1) = -0.44500e-02_wp
            htt(17,14,1) = 0.13100e-02_wp
            htt(17,15,1) = -0.30600e-02_wp
            htt(17,16,1) = 0.27000e-02_wp
            htt(17,17,1) = 0.48400e-02_wp
            htt(18,1,1)  = 0.14900e-02_wp
            htt(18,2,1)  = 0.21600e-02_wp
            htt(18,3,1)  = -0.11900e-02_wp
            htt(18,4,1)  = -0.51000e-03_wp
            htt(18,5,1)  = -0.18100e-02_wp
            htt(18,6,1)  = -0.28100e-02_wp
            htt(18,7,1)  = 0.76000e-03_wp
            htt(18,8,1)  = -0.23000e-03_wp
            htt(18,9,1)  = 0.98000e-03_wp
            htt(18,10,1) = 0.75000e-03_wp
            htt(18,11,1) = 0.29000e-03_wp
            htt(18,12,1) = 0.60000e-04_wp
            htt(18,13,1) = -0.31600e-02_wp
            htt(18,14,1) = 0.18600e-02_wp
            htt(18,15,1) = -0.29300e-02_wp
            htt(18,16,1) = 0.20900e-02_wp
            htt(18,17,1) = 0.28200e-02_wp
            htt(18,18,1) = 0.55200e-02_wp
            htt(19,1,1)  = -0.50000e-03_wp
            htt(19,2,1)  = -0.97000e-03_wp
            htt(19,3,1)  = 0.63000e-03_wp
            htt(19,4,1)  = 0.15000e-03_wp
            htt(19,5,1)  = 0.88000e-03_wp
            htt(19,6,1)  = 0.40000e-03_wp
            htt(19,7,1)  = -0.21000e-03_wp
            htt(19,8,1)  = -0.69000e-03_wp
            htt(19,9,1)  = -0.59000e-03_wp
            htt(19,10,1) = -0.12200e-02_wp
            htt(19,11,1) = -0.10000e-03_wp
            htt(19,12,1) = -0.74000e-03_wp
            htt(19,13,1) = 0.13200e-02_wp
            htt(19,14,1) = -0.33400e-02_wp
            htt(19,15,1) = 0.19500e-02_wp
            htt(19,16,1) = -0.27100e-02_wp
            htt(19,17,1) = -0.93000e-03_wp
            htt(19,18,1) = -0.22100e-02_wp
            htt(19,19,1) = 0.37400e-02_wp
            htt(20,1,1)  = -0.35000e-03_wp
            htt(20,2,1)  = -0.22000e-03_wp
            htt(20,3,1)  = -0.11000e-03_wp
            htt(20,4,1)  = 0.20000e-04_wp
            htt(20,5,1)  = 0.10000e-04_wp
            htt(20,6,1)  = -0.36000e-03_wp
            htt(20,7,1)  = -0.23000e-03_wp
            htt(20,8,1)  = -0.60000e-03_wp
            htt(20,9,1)  = -0.91000e-03_wp
            htt(20,10,1) = -0.38000e-03_wp
            htt(20,11,1) = -0.10400e-02_wp
            htt(20,12,1) = -0.51000e-03_wp
            htt(20,13,1) = 0.26000e-02_wp
            htt(20,14,1) = -0.27900e-02_wp
            htt(20,15,1) = 0.20800e-02_wp
            htt(20,16,1) = -0.33200e-02_wp
            htt(20,17,1) = -0.29500e-02_wp
            htt(20,18,1) = 0.49000e-03_wp
            htt(20,19,1) = 0.24600e-02_wp
            htt(20,20,1) = 0.59000e-02_wp
            htt(21,1,1)  = 0.19800e-02_wp
            htt(21,2,1)  = 0.17000e-02_wp
            htt(21,3,1)  = -0.76000e-03_wp
            htt(21,4,1)  = -0.10500e-02_wp
            htt(21,5,1)  = -0.29500e-02_wp
            htt(21,6,1)  = -0.43000e-03_wp
            htt(21,7,1)  = -0.30000e-03_wp
            htt(21,8,1)  = 0.22100e-02_wp
            htt(21,9,1)  = -0.98000e-03_wp
            htt(21,10,1) = -0.73000e-03_wp
            htt(21,11,1) = -0.69000e-03_wp
            htt(21,12,1) = -0.44000e-03_wp
            htt(21,13,1) = 0.11100e-02_wp
            htt(21,14,1) = -0.81000e-03_wp
            htt(21,15,1) = 0.87000e-03_wp
            htt(21,16,1) = -0.10600e-02_wp
            htt(21,17,1) = -0.20800e-02_wp
            htt(21,18,1) = -0.34000e-03_wp
            htt(21,19,1) = -0.15000e-03_wp
            htt(21,20,1) = 0.15900e-02_wp
            htt(21,21,1) = 0.28700e-02_wp
            htt(22,1,1)  = -0.50000e-04_wp
            htt(22,2,1)  = 0.12500e-02_wp
            htt(22,3,1)  = -0.12100e-02_wp
            htt(22,4,1)  = 0.10000e-03_wp
            htt(22,5,1)  = 0.25000e-03_wp
            htt(22,6,1)  = -0.22800e-02_wp
            htt(22,7,1)  = 0.13700e-02_wp
            htt(22,8,1)  = -0.11700e-02_wp
            htt(22,9,1)  = -0.62000e-03_wp
            htt(22,10,1) = 0.37000e-03_wp
            htt(22,11,1) = -0.19600e-02_wp
            htt(22,12,1) = -0.97000e-03_wp
            htt(22,13,1) = 0.11300e-02_wp
            htt(22,14,1) = -0.13000e-03_wp
            htt(22,15,1) = 0.15000e-03_wp
            htt(22,16,1) = -0.11100e-02_wp
            htt(22,17,1) = -0.92000e-03_wp
            htt(22,18,1) = 0.70000e-03_wp
            htt(22,19,1) = 0.35000e-03_wp
            htt(22,20,1) = 0.19700e-02_wp
            htt(22,21,1) = 0.22000e-03_wp
            htt(22,22,1) = 0.28800e-02_wp
            htt(23,1,1)  = 0.18300e-02_wp
            htt(23,2,1)  = 0.94000e-03_wp
            htt(23,3,1)  = -0.20000e-04_wp
            htt(23,4,1)  = -0.91000e-03_wp
            htt(23,5,1)  = -0.20700e-02_wp
            htt(23,6,1)  = 0.34000e-03_wp
            htt(23,7,1)  = -0.28000e-03_wp
            htt(23,8,1)  = 0.21200e-02_wp
            htt(23,9,1)  = -0.66000e-03_wp
            htt(23,10,1) = -0.15900e-02_wp
            htt(23,11,1) = 0.26000e-03_wp
            htt(23,12,1) = -0.67000e-03_wp
            htt(23,13,1) = -0.18000e-03_wp
            htt(23,14,1) = -0.13600e-02_wp
            htt(23,15,1) = 0.74000e-03_wp
            htt(23,16,1) = -0.44000e-03_wp
            htt(23,17,1) = -0.50000e-04_wp
            htt(23,18,1) = -0.30600e-02_wp
            htt(23,19,1) = 0.11300e-02_wp
            htt(23,20,1) = -0.18700e-02_wp
            htt(23,21,1) = 0.11100e-02_wp
            htt(23,22,1) = -0.14200e-02_wp
            htt(23,23,1) = 0.41400e-02_wp
            htt(24,1,1)  = -0.20000e-03_wp
            htt(24,2,1)  = 0.49000e-03_wp
            htt(24,3,1)  = -0.47000e-03_wp
            htt(24,4,1)  = 0.23000e-03_wp
            htt(24,5,1)  = 0.11300e-02_wp
            htt(24,6,1)  = -0.15100e-02_wp
            htt(24,7,1)  = 0.13900e-02_wp
            htt(24,8,1)  = -0.12600e-02_wp
            htt(24,9,1)  = -0.29000e-03_wp
            htt(24,10,1) = -0.48000e-03_wp
            htt(24,11,1) = -0.10100e-02_wp
            htt(24,12,1) = -0.12000e-02_wp
            htt(24,13,1) = -0.17000e-03_wp
            htt(24,14,1) = -0.69000e-03_wp
            htt(24,15,1) = 0.20000e-04_wp
            htt(24,16,1) = -0.50000e-03_wp
            htt(24,17,1) = 0.11100e-02_wp
            htt(24,18,1) = -0.20100e-02_wp
            htt(24,19,1) = 0.16300e-02_wp
            htt(24,20,1) = -0.14900e-02_wp
            htt(24,21,1) = -0.15300e-02_wp
            htt(24,22,1) = 0.12400e-02_wp
            htt(24,23,1) = 0.16100e-02_wp
            htt(24,24,1) = 0.43900e-02_wp
            ! Hessian elements for anchor point 1
            hbs( 1,1,1)  = -0.30890e-01_wp
            hbs( 1,2,1)  = -0.30720e-01_wp
            hbs( 1,3,1)  = -0.30470e-01_wp
            hbs( 1,4,1)  = -0.89900e-02_wp
            hbs( 1,5,1)  = 0.47600e-02_wp
            hbs( 1,6,1)  = 0.39030e-01_wp
            hbs( 1,7,1)  = -0.44000e-03_wp
            hbs( 1,8,1)  = 0.39840e-01_wp
            hbs( 1,9,1)  = -0.25000e-03_wp
            hbs( 1,10,1) = -0.97500e-02_wp
            hbs( 1,11,1) = -0.42000e-03_wp
            hbs( 1,12,1) = 0.49400e-02_wp
            hbs( 2,1,1)  = 0.43860e-01_wp
            hbs( 2,2,1)  = -0.11950e-01_wp
            hbs( 2,3,1)  = 0.67900e-02_wp
            hbs( 2,4,1)  = 0.23800e-02_wp
            hbs( 2,5,1)  = -0.22800e-02_wp
            hbs( 2,6,1)  = -0.26060e-01_wp
            hbs( 2,7,1)  = 0.11000e-02_wp
            hbs( 2,8,1)  = -0.14350e-01_wp
            hbs( 2,9,1)  = 0.40000e-04_wp
            hbs( 2,10,1) = 0.41400e-02_wp
            hbs( 2,11,1) = -0.73000e-03_wp
            hbs( 2,12,1) = -0.48000e-03_wp
            hbs( 3,1,1)  = -0.12970e-01_wp
            hbs( 3,2,1)  = 0.42680e-01_wp
            hbs( 3,3,1)  = 0.23690e-01_wp
            hbs( 3,4,1)  = 0.66100e-02_wp
            hbs( 3,5,1)  = -0.24900e-02_wp
            hbs( 3,6,1)  = -0.12970e-01_wp
            hbs( 3,7,1)  = -0.66000e-03_wp
            hbs( 3,8,1)  = -0.25490e-01_wp
            hbs( 3,9,1)  = 0.21000e-03_wp
            hbs( 3,10,1) = 0.56200e-02_wp
            hbs( 3,11,1) = 0.11500e-02_wp
            hbs( 3,12,1) = -0.44600e-02_wp
            hbs( 4,1,1)  = -0.21190e-01_wp
            hbs( 4,2,1)  = -0.10770e-01_wp
            hbs( 4,3,1)  = 0.12660e-01_wp
            hbs( 4,4,1)  = -0.29880e-01_wp
            hbs( 4,5,1)  = -0.90800e-02_wp
            hbs( 4,6,1)  = -0.51700e-02_wp
            hbs( 4,7,1)  = 0.49300e-02_wp
            hbs( 4,8,1)  = 0.30850e-01_wp
            hbs( 4,9,1)  = -0.40000e-03_wp
            hbs( 4,10,1) = 0.40050e-01_wp
            hbs( 4,11,1) = -0.62000e-03_wp
            hbs( 4,12,1) = 0.12000e-03_wp
            hbs( 5,1,1)  = 0.25060e-01_wp
            hbs( 5,2,1)  = 0.46900e-02_wp
            hbs( 5,3,1)  = -0.97000e-02_wp
            hbs( 5,4,1)  = 0.24100e-02_wp
            hbs( 5,5,1)  = 0.62100e-02_wp
            hbs( 5,6,1)  = 0.30800e-02_wp
            hbs( 5,7,1)  = -0.21300e-02_wp
            hbs( 5,8,1)  = -0.14880e-01_wp
            hbs( 5,9,1)  = -0.39000e-03_wp
            hbs( 5,10,1) = -0.21940e-01_wp
            hbs( 5,11,1) = 0.31000e-03_wp
            hbs( 5,12,1) = 0.35000e-03_wp
            hbs( 6,1,1)  = -0.38700e-02_wp
            hbs( 6,2,1)  = 0.60800e-02_wp
            hbs( 6,3,1)  = -0.29500e-02_wp
            hbs( 6,4,1)  = 0.27480e-01_wp
            hbs( 6,5,1)  = 0.28700e-02_wp
            hbs( 6,6,1)  = 0.21000e-02_wp
            hbs( 6,7,1)  = -0.28000e-02_wp
            hbs( 6,8,1)  = -0.15960e-01_wp
            hbs( 6,9,1)  = 0.80000e-03_wp
            hbs( 6,10,1) = -0.18110e-01_wp
            hbs( 6,11,1) = 0.31000e-03_wp
            hbs( 6,12,1) = -0.47000e-03_wp
            hbs( 7,1,1)  = -0.88700e-02_wp
            hbs( 7,2,1)  = 0.35920e-01_wp
            hbs( 7,3,1)  = 0.96000e-03_wp
            hbs( 7,4,1)  = -0.23260e-01_wp
            hbs( 7,5,1)  = 0.57800e-02_wp
            hbs( 7,6,1)  = -0.29090e-01_wp
            hbs( 7,7,1)  = -0.87400e-02_wp
            hbs( 7,8,1)  = -0.88100e-02_wp
            hbs( 7,9,1)  = 0.49600e-02_wp
            hbs( 7,10,1) = 0.32900e-01_wp
            hbs( 7,11,1) = 0.60000e-04_wp
            hbs( 7,12,1) = -0.10800e-02_wp
            hbs( 8,1,1)  = 0.36600e-02_wp
            hbs( 8,2,1)  = -0.18680e-01_wp
            hbs( 8,3,1)  = 0.82000e-03_wp
            hbs( 8,4,1)  = 0.25710e-01_wp
            hbs( 8,5,1)  = -0.33200e-02_wp
            hbs( 8,6,1)  = -0.40000e-03_wp
            hbs( 8,7,1)  = 0.44400e-02_wp
            hbs( 8,8,1)  = 0.48100e-02_wp
            hbs( 8,9,1)  = -0.20200e-02_wp
            hbs( 8,10,1) = -0.14750e-01_wp
            hbs( 8,11,1) = -0.62000e-03_wp
            hbs( 8,12,1) = 0.56000e-03_wp
            hbs( 9,1,1)  = 0.52100e-02_wp
            hbs( 9,2,1)  = -0.17240e-01_wp
            hbs( 9,3,1)  = -0.17800e-02_wp
            hbs( 9,4,1)  = -0.24500e-02_wp
            hbs( 9,5,1)  = -0.24600e-02_wp
            hbs( 9,6,1)  = 0.29490e-01_wp
            hbs( 9,7,1)  = 0.43000e-02_wp
            hbs( 9,8,1)  = 0.40000e-02_wp
            hbs( 9,9,1)  = -0.29400e-02_wp
            hbs( 9,10,1) = -0.18150e-01_wp
            hbs( 9,11,1) = 0.56000e-03_wp
            hbs( 9,12,1) = 0.52000e-03_wp
            hbs(10,1,1)  = 0.36080e-01_wp
            hbs(10,2,1)  = 0.38610e-01_wp
            hbs(10,3,1)  = 0.39500e-02_wp
            hbs(10,4,1)  = -0.11320e-01_wp
            hbs(10,5,1)  = -0.89000e-03_wp
            hbs(10,6,1)  = -0.26740e-01_wp
            hbs(10,7,1)  = 0.48000e-02_wp
            hbs(10,8,1)  = -0.26910e-01_wp
            hbs(10,9,1)  = -0.88100e-02_wp
            hbs(10,10,1) = -0.11650e-01_wp
            hbs(10,11,1) = 0.47600e-02_wp
            hbs(10,12,1) = -0.70000e-03_wp
            hbs(11,1,1)  = -0.20240e-01_wp
            hbs(11,2,1)  = -0.17250e-01_wp
            hbs(11,3,1)  = -0.21300e-02_wp
            hbs(11,4,1)  = 0.55100e-02_wp
            hbs(11,5,1)  = 0.11200e-02_wp
            hbs(11,6,1)  = 0.27050e-01_wp
            hbs(11,7,1)  = -0.28200e-02_wp
            hbs(11,8,1)  = 0.25000e-03_wp
            hbs(11,9,1)  = 0.44500e-02_wp
            hbs(11,10,1) = 0.58600e-02_wp
            hbs(11,11,1) = -0.19800e-02_wp
            hbs(11,12,1) = -0.20000e-03_wp
            hbs(12,1,1)  = -0.15840e-01_wp
            hbs(12,2,1)  = -0.21360e-01_wp
            hbs(12,3,1)  = -0.18200e-02_wp
            hbs(12,4,1)  = 0.58100e-02_wp
            hbs(12,5,1)  = -0.23000e-03_wp
            hbs(12,6,1)  = -0.31000e-03_wp
            hbs(12,7,1)  = -0.19800e-02_wp
            hbs(12,8,1)  = 0.26660e-01_wp
            hbs(12,9,1)  = 0.43600e-02_wp
            hbs(12,10,1) = 0.57900e-02_wp
            hbs(12,11,1) = -0.27800e-02_wp
            hbs(12,12,1) = 0.90000e-03_wp
            hbs(13,1,1)  = 0.35090e-01_wp
            hbs(13,2,1)  = -0.85500e-02_wp
            hbs(13,3,1)  = 0.10200e-02_wp
            hbs(13,4,1)  = 0.32180e-01_wp
            hbs(13,5,1)  = -0.97000e-03_wp
            hbs(13,6,1)  = -0.94600e-02_wp
            hbs(13,7,1)  = 0.10000e-03_wp
            hbs(13,8,1)  = -0.29520e-01_wp
            hbs(13,9,1)  = 0.50200e-02_wp
            hbs(13,10,1) = -0.23910e-01_wp
            hbs(13,11,1) = -0.86300e-02_wp
            hbs(13,12,1) = 0.52600e-02_wp
            hbs(14,1,1)  = -0.16950e-01_wp
            hbs(14,2,1)  = 0.51900e-02_wp
            hbs(14,3,1)  = -0.19500e-02_wp
            hbs(14,4,1)  = -0.17590e-01_wp
            hbs(14,5,1)  = 0.48000e-03_wp
            hbs(14,6,1)  = 0.43300e-02_wp
            hbs(14,7,1)  = 0.54000e-03_wp
            hbs(14,8,1)  = 0.29620e-01_wp
            hbs(14,9,1)  = -0.29700e-02_wp
            hbs(14,10,1) = -0.21100e-02_wp
            hbs(14,11,1) = 0.42400e-02_wp
            hbs(14,12,1) = -0.22400e-02_wp
            hbs(15,1,1)  = -0.18140e-01_wp
            hbs(15,2,1)  = 0.33600e-02_wp
            hbs(15,3,1)  = 0.94000e-03_wp
            hbs(15,4,1)  = -0.14580e-01_wp
            hbs(15,5,1)  = 0.48000e-03_wp
            hbs(15,6,1)  = 0.51200e-02_wp
            hbs(15,7,1)  = -0.64000e-03_wp
            hbs(15,8,1)  = -0.10000e-03_wp
            hbs(15,9,1)  = -0.20500e-02_wp
            hbs(15,10,1) = 0.26020e-01_wp
            hbs(15,11,1) = 0.43900e-02_wp
            hbs(15,12,1) = -0.30200e-02_wp
            hbs(16,1,1)  = -0.10210e-01_wp
            hbs(16,2,1)  = -0.24490e-01_wp
            hbs(16,3,1)  = 0.11890e-01_wp
            hbs(16,4,1)  = 0.41280e-01_wp
            hbs(16,5,1)  = 0.40000e-03_wp
            hbs(16,6,1)  = 0.31430e-01_wp
            hbs(16,7,1)  = -0.65000e-03_wp
            hbs(16,8,1)  = -0.54400e-02_wp
            hbs(16,9,1)  = -0.51000e-03_wp
            hbs(16,10,1) = -0.27630e-01_wp
            hbs(16,11,1) = 0.48400e-02_wp
            hbs(16,12,1) = -0.85500e-02_wp
            hbs(17,1,1)  = 0.43200e-02_wp
            hbs(17,2,1)  = 0.25240e-01_wp
            hbs(17,3,1)  = -0.81200e-02_wp
            hbs(17,4,1)  = -0.22630e-01_wp
            hbs(17,5,1)  = 0.41000e-03_wp
            hbs(17,6,1)  = -0.14630e-01_wp
            hbs(17,7,1)  = 0.33000e-03_wp
            hbs(17,8,1)  = 0.32400e-02_wp
            hbs(17,9,1)  = -0.33000e-03_wp
            hbs(17,10,1) = 0.18800e-02_wp
            hbs(17,11,1) = -0.21000e-02_wp
            hbs(17,12,1) = 0.69700e-02_wp
            hbs(18,1,1)  = 0.58800e-02_wp
            hbs(18,2,1)  = -0.75000e-03_wp
            hbs(18,3,1)  = -0.37700e-02_wp
            hbs(18,4,1)  = -0.18650e-01_wp
            hbs(18,5,1)  = -0.80000e-03_wp
            hbs(18,6,1)  = -0.16800e-01_wp
            hbs(18,7,1)  = 0.32000e-03_wp
            hbs(18,8,1)  = 0.22000e-02_wp
            hbs(18,9,1)  = 0.85000e-03_wp
            hbs(18,10,1) = 0.25750e-01_wp
            hbs(18,11,1) = -0.27400e-02_wp
            hbs(18,12,1) = 0.15700e-02_wp
            ! Hessian elements for anchor point 2
            hss( 1,1,2)  = 0.34759e+00_wp
            hss( 2,1,2)  = 0.17810e-01_wp
            hss( 2,2,2)  = 0.35623e+00_wp
            hss( 3,1,2)  = 0.29010e-01_wp
            hss( 3,2,2)  = 0.32710e-01_wp
            hss( 3,3,2)  = 0.46124e+00_wp
            hss( 4,1,2)  = -0.19600e-02_wp
            hss( 4,2,2)  = 0.14870e-01_wp
            hss( 4,3,2)  = -0.89100e-02_wp
            hss( 4,4,2)  = 0.38581e+00_wp
            hss( 5,1,2)  = 0.38800e-02_wp
            hss( 5,2,2)  = -0.77000e-03_wp
            hss( 5,3,2)  = 0.71700e-02_wp
            hss( 5,4,2)  = 0.36500e-02_wp
            hss( 5,5,2)  = 0.35807e+00_wp
            hss( 6,1,2)  = 0.33890e-01_wp
            hss( 6,2,2)  = 0.44570e-01_wp
            hss( 6,3,2)  = -0.43100e-02_wp
            hss( 6,4,2)  = -0.36100e-02_wp
            hss( 6,5,2)  = 0.55000e-03_wp
            hss( 6,6,2)  = 0.35931e+00_wp
            hss( 7,1,2)  = -0.54000e-03_wp
            hss( 7,2,2)  = -0.63000e-03_wp
            hss( 7,3,2)  = 0.57000e-03_wp
            hss( 7,4,2)  = 0.39700e-02_wp
            hss( 7,5,2)  = 0.39000e-03_wp
            hss( 7,6,2)  = 0.39800e-02_wp
            hss( 7,7,2)  = 0.37091e+00_wp
            hss( 8,1,2)  = 0.42030e-01_wp
            hss( 8,2,2)  = 0.35060e-01_wp
            hss( 8,3,2)  = -0.52800e-02_wp
            hss( 8,4,2)  = 0.24280e-01_wp
            hss( 8,5,2)  = -0.11100e-02_wp
            hss( 8,6,2)  = 0.17880e-01_wp
            hss( 8,7,2)  = -0.10400e-02_wp
            hss( 8,8,2)  = 0.35888e+00_wp
            hss( 9,1,2)  = -0.53000e-03_wp
            hss( 9,2,2)  = -0.57000e-03_wp
            hss( 9,3,2)  = 0.26000e-03_wp
            hss( 9,4,2)  = -0.39000e-03_wp
            hss( 9,5,2)  = 0.19000e-03_wp
            hss( 9,6,2)  = 0.40100e-02_wp
            hss( 9,7,2)  = 0.49000e-03_wp
            hss( 9,8,2)  = 0.39600e-02_wp
            hss( 9,9,2)  = 0.36356e+00_wp
            hss(10,1,2)  = 0.19520e-01_wp
            hss(10,2,2)  = -0.69300e-02_wp
            hss(10,3,2)  = -0.48800e-02_wp
            hss(10,4,2)  = 0.64010e-01_wp
            hss(10,5,2)  = -0.43000e-03_wp
            hss(10,6,2)  = 0.21230e-01_wp
            hss(10,7,2)  = -0.10900e-02_wp
            hss(10,8,2)  = -0.77000e-03_wp
            hss(10,9,2)  = -0.19000e-03_wp
            hss(10,10,2) = 0.38116e+00_wp
            hss(11,1,2)  = -0.69000e-03_wp
            hss(11,2,2)  = -0.50000e-03_wp
            hss(11,3,2)  = 0.52000e-03_wp
            hss(11,4,2)  = -0.11600e-02_wp
            hss(11,5,2)  = 0.80000e-04_wp
            hss(11,6,2)  = -0.10000e-02_wp
            hss(11,7,2)  = 0.22000e-03_wp
            hss(11,8,2)  = 0.39700e-02_wp
            hss(11,9,2)  = 0.48000e-03_wp
            hss(11,10,2) = 0.39800e-02_wp
            hss(11,11,2) = 0.37152e+00_wp
            hss(12,1,2)  = -0.16500e-02_wp
            hss(12,2,2)  = 0.39400e-02_wp
            hss(12,3,2)  = -0.19000e-03_wp
            hss(12,4,2)  = -0.10000e-03_wp
            hss(12,5,2)  = 0.35000e-03_wp
            hss(12,6,2)  = -0.10400e-02_wp
            hss(12,7,2)  = 0.30000e-04_wp
            hss(12,8,2)  = 0.57000e-03_wp
            hss(12,9,2)  = 0.12000e-03_wp
            hss(12,10,2) = 0.23600e-02_wp
            hss(12,11,2) = 0.40000e-03_wp
            hss(12,12,2) = 0.37254e+00_wp
            ! Hessian elements for anchor point 2
            hbb( 1,1,2)  = 0.82450e-01_wp
            hbb( 2,1,2)  = -0.46090e-01_wp
            hbb( 2,2,2)  = 0.14487e+00_wp
            hbb( 3,1,2)  = -0.36360e-01_wp
            hbb( 3,2,2)  = -0.98780e-01_wp
            hbb( 3,3,2)  = 0.13515e+00_wp
            hbb( 4,1,2)  = -0.36000e-01_wp
            hbb( 4,2,2)  = 0.23310e-01_wp
            hbb( 4,3,2)  = 0.12690e-01_wp
            hbb( 4,4,2)  = 0.66170e-01_wp
            hbb( 5,1,2)  = 0.20240e-01_wp
            hbb( 5,2,2)  = -0.12480e-01_wp
            hbb( 5,3,2)  = -0.77700e-02_wp
            hbb( 5,4,2)  = -0.33030e-01_wp
            hbb( 5,5,2)  = 0.71950e-01_wp
            hbb( 6,1,2)  = 0.15760e-01_wp
            hbb( 6,2,2)  = -0.10840e-01_wp
            hbb( 6,3,2)  = -0.49200e-02_wp
            hbb( 6,4,2)  = -0.33140e-01_wp
            hbb( 6,5,2)  = -0.38920e-01_wp
            hbb( 6,6,2)  = 0.72060e-01_wp
            hbb( 7,1,2)  = -0.71500e-02_wp
            hbb( 7,2,2)  = -0.39000e-02_wp
            hbb( 7,3,2)  = 0.11050e-01_wp
            hbb( 7,4,2)  = -0.25740e-01_wp
            hbb( 7,5,2)  = 0.93100e-02_wp
            hbb( 7,6,2)  = 0.16430e-01_wp
            hbb( 7,7,2)  = 0.64420e-01_wp
            hbb( 8,1,2)  = -0.22800e-02_wp
            hbb( 8,2,2)  = 0.74300e-02_wp
            hbb( 8,3,2)  = -0.51500e-02_wp
            hbb( 8,4,2)  = 0.15690e-01_wp
            hbb( 8,5,2)  = -0.49500e-02_wp
            hbb( 8,6,2)  = -0.10740e-01_wp
            hbb( 8,7,2)  = -0.31490e-01_wp
            hbb( 8,8,2)  = 0.72950e-01_wp
            hbb( 9,1,2)  = 0.94300e-02_wp
            hbb( 9,2,2)  = -0.35300e-02_wp
            hbb( 9,3,2)  = -0.59100e-02_wp
            hbb( 9,4,2)  = 0.10050e-01_wp
            hbb( 9,5,2)  = -0.43500e-02_wp
            hbb( 9,6,2)  = -0.56900e-02_wp
            hbb( 9,7,2)  = -0.32930e-01_wp
            hbb( 9,8,2)  = -0.41460e-01_wp
            hbb( 9,9,2)  = 0.74390e-01_wp
            hbb(10,1,2)  = 0.45500e-02_wp
            hbb(10,2,2)  = -0.20200e-02_wp
            hbb(10,3,2)  = -0.25300e-02_wp
            hbb(10,4,2)  = -0.64100e-02_wp
            hbb(10,5,2)  = 0.91900e-02_wp
            hbb(10,6,2)  = -0.27700e-02_wp
            hbb(10,7,2)  = -0.33410e-01_wp
            hbb(10,8,2)  = 0.13410e-01_wp
            hbb(10,9,2)  = 0.20000e-01_wp
            hbb(10,10,2) = 0.75440e-01_wp
            hbb(11,1,2)  = -0.22100e-02_wp
            hbb(11,2,2)  = 0.12300e-02_wp
            hbb(11,3,2)  = 0.98000e-03_wp
            hbb(11,4,2)  = -0.23100e-02_wp
            hbb(11,5,2)  = -0.32200e-02_wp
            hbb(11,6,2)  = 0.55300e-02_wp
            hbb(11,7,2)  = 0.20020e-01_wp
            hbb(11,8,2)  = -0.74100e-02_wp
            hbb(11,9,2)  = -0.12610e-01_wp
            hbb(11,10,2) = -0.37790e-01_wp
            hbb(11,11,2) = 0.78360e-01_wp
            hbb(12,1,2)  = -0.23400e-02_wp
            hbb(12,2,2)  = 0.79000e-03_wp
            hbb(12,3,2)  = 0.15500e-02_wp
            hbb(12,4,2)  = 0.87200e-02_wp
            hbb(12,5,2)  = -0.59700e-02_wp
            hbb(12,6,2)  = -0.27500e-02_wp
            hbb(12,7,2)  = 0.13390e-01_wp
            hbb(12,8,2)  = -0.60000e-02_wp
            hbb(12,9,2)  = -0.73900e-02_wp
            hbb(12,10,2) = -0.37640e-01_wp
            hbb(12,11,2) = -0.40570e-01_wp
            hbb(12,12,2) = 0.78220e-01_wp
            hbb(13,1,2)  = -0.67600e-02_wp
            hbb(13,2,2)  = 0.10710e-01_wp
            hbb(13,3,2)  = -0.39500e-02_wp
            hbb(13,4,2)  = -0.11500e-02_wp
            hbb(13,5,2)  = 0.12900e-02_wp
            hbb(13,6,2)  = -0.14000e-03_wp
            hbb(13,7,2)  = 0.29100e-02_wp
            hbb(13,8,2)  = 0.41300e-02_wp
            hbb(13,9,2)  = -0.70300e-02_wp
            hbb(13,10,2) = -0.33760e-01_wp
            hbb(13,11,2) = 0.13510e-01_wp
            hbb(13,12,2) = 0.20240e-01_wp
            hbb(13,13,2) = 0.65300e-01_wp
            hbb(14,1,2)  = 0.91200e-02_wp
            hbb(14,2,2)  = -0.60500e-02_wp
            hbb(14,3,2)  = -0.30700e-02_wp
            hbb(14,4,2)  = 0.66000e-03_wp
            hbb(14,5,2)  = -0.85000e-03_wp
            hbb(14,6,2)  = 0.19000e-03_wp
            hbb(14,7,2)  = -0.71100e-02_wp
            hbb(14,8,2)  = -0.51000e-03_wp
            hbb(14,9,2)  = 0.76200e-02_wp
            hbb(14,10,2) = 0.20030e-01_wp
            hbb(14,11,2) = -0.73800e-02_wp
            hbb(14,12,2) = -0.12650e-01_wp
            hbb(14,13,2) = -0.33180e-01_wp
            hbb(14,14,2) = 0.74350e-01_wp
            hbb(15,1,2)  = -0.23600e-02_wp
            hbb(15,2,2)  = -0.46500e-02_wp
            hbb(15,3,2)  = 0.70200e-02_wp
            hbb(15,4,2)  = 0.49000e-03_wp
            hbb(15,5,2)  = -0.44000e-03_wp
            hbb(15,6,2)  = -0.50000e-04_wp
            hbb(15,7,2)  = 0.42000e-02_wp
            hbb(15,8,2)  = -0.36200e-02_wp
            hbb(15,9,2)  = -0.59000e-03_wp
            hbb(15,10,2) = 0.13730e-01_wp
            hbb(15,11,2) = -0.61300e-02_wp
            hbb(15,12,2) = -0.76000e-02_wp
            hbb(15,13,2) = -0.32120e-01_wp
            hbb(15,14,2) = -0.41170e-01_wp
            hbb(15,15,2) = 0.73290e-01_wp
            hbb(16,1,2)  = -0.37090e-01_wp
            hbb(16,2,2)  = 0.17990e-01_wp
            hbb(16,3,2)  = 0.19100e-01_wp
            hbb(16,4,2)  = 0.31400e-02_wp
            hbb(16,5,2)  = -0.70000e-02_wp
            hbb(16,6,2)  = 0.38600e-02_wp
            hbb(16,7,2)  = -0.10300e-02_wp
            hbb(16,8,2)  = 0.54000e-03_wp
            hbb(16,9,2)  = 0.48000e-03_wp
            hbb(16,10,2) = -0.64000e-02_wp
            hbb(16,11,2) = 0.87800e-02_wp
            hbb(16,12,2) = -0.23700e-02_wp
            hbb(16,13,2) = -0.26530e-01_wp
            hbb(16,14,2) = 0.10480e-01_wp
            hbb(16,15,2) = 0.16060e-01_wp
            hbb(16,16,2) = 0.67910e-01_wp
            hbb(17,1,2)  = 0.22610e-01_wp
            hbb(17,2,2)  = -0.90200e-02_wp
            hbb(17,3,2)  = -0.13590e-01_wp
            hbb(17,4,2)  = -0.76000e-02_wp
            hbb(17,5,2)  = 0.77900e-02_wp
            hbb(17,6,2)  = -0.19000e-03_wp
            hbb(17,7,2)  = 0.87000e-03_wp
            hbb(17,8,2)  = -0.25000e-03_wp
            hbb(17,9,2)  = -0.62000e-03_wp
            hbb(17,10,2) = 0.91600e-02_wp
            hbb(17,11,2) = -0.60700e-02_wp
            hbb(17,12,2) = -0.31000e-02_wp
            hbb(17,13,2) = 0.96800e-02_wp
            hbb(17,14,2) = -0.46000e-02_wp
            hbb(17,15,2) = -0.50800e-02_wp
            hbb(17,16,2) = -0.34720e-01_wp
            hbb(17,17,2) = 0.72390e-01_wp
            hbb(18,1,2)  = 0.14470e-01_wp
            hbb(18,2,2)  = -0.89700e-02_wp
            hbb(18,3,2)  = -0.55000e-02_wp
            hbb(18,4,2)  = 0.44600e-02_wp
            hbb(18,5,2)  = -0.79000e-03_wp
            hbb(18,6,2)  = -0.36700e-02_wp
            hbb(18,7,2)  = 0.16000e-03_wp
            hbb(18,8,2)  = -0.29000e-03_wp
            hbb(18,9,2)  = 0.13000e-03_wp
            hbb(18,10,2) = -0.27600e-02_wp
            hbb(18,11,2) = -0.27100e-02_wp
            hbb(18,12,2) = 0.54700e-02_wp
            hbb(18,13,2) = 0.16860e-01_wp
            hbb(18,14,2) = -0.58800e-02_wp
            hbb(18,15,2) = -0.10980e-01_wp
            hbb(18,16,2) = -0.33190e-01_wp
            hbb(18,17,2) = -0.37670e-01_wp
            hbb(18,18,2) = 0.70860e-01_wp
            ! Hessian elements for anchor point 2
            htt( 1,1,2)  = 0.60400e-02_wp
            htt( 2,1,2)  = 0.38400e-02_wp
            htt( 2,2,2)  = 0.50900e-02_wp
            htt( 3,1,2)  = -0.18500e-02_wp
            htt( 3,2,2)  = -0.33100e-02_wp
            htt( 3,3,2)  = 0.75100e-02_wp
            htt( 4,1,2)  = -0.40500e-02_wp
            htt( 4,2,2)  = -0.20700e-02_wp
            htt( 4,3,2)  = 0.60500e-02_wp
            htt( 4,4,2)  = 0.80200e-02_wp
            htt( 5,1,2)  = -0.59200e-02_wp
            htt( 5,2,2)  = -0.42600e-02_wp
            htt( 5,3,2)  = 0.21200e-02_wp
            htt( 5,4,2)  = 0.37800e-02_wp
            htt( 5,5,2)  = 0.65600e-02_wp
            htt( 6,1,2)  = -0.41200e-02_wp
            htt( 6,2,2)  = -0.39000e-02_wp
            htt( 6,3,2)  = 0.26900e-02_wp
            htt( 6,4,2)  = 0.29100e-02_wp
            htt( 6,5,2)  = 0.37400e-02_wp
            htt( 6,6,2)  = 0.55900e-02_wp
            htt( 7,1,2)  = 0.15700e-02_wp
            htt( 7,2,2)  = 0.25400e-02_wp
            htt( 7,3,2)  = -0.67700e-02_wp
            htt( 7,4,2)  = -0.58000e-02_wp
            htt( 7,5,2)  = -0.10700e-02_wp
            htt( 7,6,2)  = -0.27200e-02_wp
            htt( 7,7,2)  = 0.68600e-02_wp
            htt( 8,1,2)  = 0.33800e-02_wp
            htt( 8,2,2)  = 0.29000e-02_wp
            htt( 8,3,2)  = -0.62000e-02_wp
            htt( 8,4,2)  = -0.66800e-02_wp
            htt( 8,5,2)  = -0.38900e-02_wp
            htt( 8,6,2)  = -0.88000e-03_wp
            htt( 8,7,2)  = 0.52100e-02_wp
            htt( 8,8,2)  = 0.82200e-02_wp
            htt( 9,1,2)  = -0.23900e-02_wp
            htt( 9,2,2)  = -0.68000e-03_wp
            htt( 9,3,2)  = 0.45000e-03_wp
            htt( 9,4,2)  = 0.21600e-02_wp
            htt( 9,5,2)  = 0.17400e-02_wp
            htt( 9,6,2)  = 0.14100e-02_wp
            htt( 9,7,2)  = -0.96000e-03_wp
            htt( 9,8,2)  = -0.12900e-02_wp
            htt( 9,9,2)  = 0.24600e-02_wp
            htt(10,1,2)  = -0.15600e-02_wp
            htt(10,2,2)  = 0.23000e-03_wp
            htt(10,3,2)  = 0.34000e-03_wp
            htt(10,4,2)  = 0.21300e-02_wp
            htt(10,5,2)  = 0.17600e-02_wp
            htt(10,6,2)  = 0.79000e-03_wp
            htt(10,7,2)  = -0.40000e-04_wp
            htt(10,8,2)  = -0.10100e-02_wp
            htt(10,9,2)  = 0.57000e-03_wp
            htt(10,10,2) = 0.39100e-02_wp
            htt(11,1,2)  = -0.10000e-03_wp
            htt(11,2,2)  = -0.19700e-02_wp
            htt(11,3,2)  = 0.19800e-02_wp
            htt(11,4,2)  = 0.11000e-03_wp
            htt(11,5,2)  = 0.10000e-04_wp
            htt(11,6,2)  = 0.11900e-02_wp
            htt(11,7,2)  = -0.19700e-02_wp
            htt(11,8,2)  = -0.79000e-03_wp
            htt(11,9,2)  = 0.68000e-03_wp
            htt(11,10,2) = -0.12900e-02_wp
            htt(11,11,2) = 0.26200e-02_wp
            htt(12,1,2)  = 0.72000e-03_wp
            htt(12,2,2)  = -0.10500e-02_wp
            htt(12,3,2)  = 0.18600e-02_wp
            htt(12,4,2)  = 0.90000e-04_wp
            htt(12,5,2)  = 0.30000e-04_wp
            htt(12,6,2)  = 0.57000e-03_wp
            htt(12,7,2)  = -0.10500e-02_wp
            htt(12,8,2)  = -0.51000e-03_wp
            htt(12,9,2)  = -0.12000e-02_wp
            htt(12,10,2) = 0.20500e-02_wp
            htt(12,11,2) = 0.64000e-03_wp
            htt(12,12,2) = 0.39000e-02_wp
            htt(13,1,2)  = -0.10900e-02_wp
            htt(13,2,2)  = -0.19400e-02_wp
            htt(13,3,2)  = 0.59000e-03_wp
            htt(13,4,2)  = -0.26000e-03_wp
            htt(13,5,2)  = 0.15200e-02_wp
            htt(13,6,2)  = 0.15100e-02_wp
            htt(13,7,2)  = -0.70000e-04_wp
            htt(13,8,2)  = -0.90000e-04_wp
            htt(13,9,2)  = -0.21100e-02_wp
            htt(13,10,2) = 0.12000e-03_wp
            htt(13,11,2) = -0.12300e-02_wp
            htt(13,12,2) = 0.10000e-02_wp
            htt(13,13,2) = 0.50300e-02_wp
            htt(14,1,2)  = 0.70000e-03_wp
            htt(14,2,2)  = 0.53000e-03_wp
            htt(14,3,2)  = -0.31000e-03_wp
            htt(14,4,2)  = -0.47000e-03_wp
            htt(14,5,2)  = -0.22000e-03_wp
            htt(14,6,2)  = -0.80000e-03_wp
            htt(14,7,2)  = 0.74000e-03_wp
            htt(14,8,2)  = 0.16000e-03_wp
            htt(14,9,2)  = -0.19000e-03_wp
            htt(14,10,2) = 0.13500e-02_wp
            htt(14,11,2) = -0.10000e-04_wp
            htt(14,12,2) = 0.15300e-02_wp
            htt(14,13,2) = -0.81000e-03_wp
            htt(14,14,2) = 0.39200e-02_wp
            htt(15,1,2)  = -0.19100e-02_wp
            htt(15,2,2)  = -0.28400e-02_wp
            htt(15,3,2)  = 0.70000e-03_wp
            htt(15,4,2)  = -0.23000e-03_wp
            htt(15,5,2)  = 0.15000e-02_wp
            htt(15,6,2)  = 0.21100e-02_wp
            htt(15,7,2)  = -0.98000e-03_wp
            htt(15,8,2)  = -0.37000e-03_wp
            htt(15,9,2)  = -0.24000e-03_wp
            htt(15,10,2) = -0.31800e-02_wp
            htt(15,11,2) = 0.73000e-03_wp
            htt(15,12,2) = -0.22200e-02_wp
            htt(15,13,2) = 0.28200e-02_wp
            htt(15,14,2) = -0.23300e-02_wp
            htt(15,15,2) = 0.57300e-02_wp
            htt(16,1,2)  = -0.12000e-03_wp
            htt(16,2,2)  = -0.37000e-03_wp
            htt(16,3,2)  = -0.19000e-03_wp
            htt(16,4,2)  = -0.45000e-03_wp
            htt(16,5,2)  = -0.25000e-03_wp
            htt(16,6,2)  = -0.19000e-03_wp
            htt(16,7,2)  = -0.17000e-03_wp
            htt(16,8,2)  = -0.12000e-03_wp
            htt(16,9,2)  = 0.16800e-02_wp
            htt(16,10,2) = -0.19500e-02_wp
            htt(16,11,2) = 0.19400e-02_wp
            htt(16,12,2) = -0.16900e-02_wp
            htt(16,13,2) = -0.30200e-02_wp
            htt(16,14,2) = 0.24000e-02_wp
            htt(16,15,2) = 0.58000e-03_wp
            htt(16,16,2) = 0.60000e-02_wp
            htt(17,1,2)  = 0.12100e-02_wp
            htt(17,2,2)  = 0.15200e-02_wp
            htt(17,3,2)  = -0.31000e-03_wp
            htt(17,4,2)  = 0.00000e+00_wp
            htt(17,5,2)  = -0.86000e-03_wp
            htt(17,6,2)  = -0.18800e-02_wp
            htt(17,7,2)  = 0.58000e-03_wp
            htt(17,8,2)  = -0.44000e-03_wp
            htt(17,9,2)  = 0.14500e-02_wp
            htt(17,10,2) = 0.90000e-04_wp
            htt(17,11,2) = 0.11300e-02_wp
            htt(17,12,2) = -0.23000e-03_wp
            htt(17,13,2) = -0.46000e-02_wp
            htt(17,14,2) = 0.13000e-02_wp
            htt(17,15,2) = -0.32400e-02_wp
            htt(17,16,2) = 0.26500e-02_wp
            htt(17,17,2) = 0.49600e-02_wp
            htt(18,1,2)  = 0.14500e-02_wp
            htt(18,2,2)  = 0.22200e-02_wp
            htt(18,3,2)  = -0.10600e-02_wp
            htt(18,4,2)  = -0.29000e-03_wp
            htt(18,5,2)  = -0.17100e-02_wp
            htt(18,6,2)  = -0.27100e-02_wp
            htt(18,7,2)  = 0.67000e-03_wp
            htt(18,8,2)  = -0.32000e-03_wp
            htt(18,9,2)  = 0.10000e-02_wp
            htt(18,10,2) = 0.82000e-03_wp
            htt(18,11,2) = 0.19000e-03_wp
            htt(18,12,2) = 0.10000e-04_wp
            htt(18,13,2) = -0.32400e-02_wp
            htt(18,14,2) = 0.19000e-02_wp
            htt(18,15,2) = -0.30600e-02_wp
            htt(18,16,2) = 0.20800e-02_wp
            htt(18,17,2) = 0.29700e-02_wp
            htt(18,18,2) = 0.56300e-02_wp
            htt(19,1,2)  = -0.59000e-03_wp
            htt(19,2,2)  = -0.96000e-03_wp
            htt(19,3,2)  = 0.59000e-03_wp
            htt(19,4,2)  = 0.21000e-03_wp
            htt(19,5,2)  = 0.89000e-03_wp
            htt(19,6,2)  = 0.43000e-03_wp
            htt(19,7,2)  = -0.23000e-03_wp
            htt(19,8,2)  = -0.69000e-03_wp
            htt(19,9,2)  = -0.48000e-03_wp
            htt(19,10,2) = -0.11500e-02_wp
            htt(19,11,2) = -0.90000e-04_wp
            htt(19,12,2) = -0.76000e-03_wp
            htt(19,13,2) = 0.12600e-02_wp
            htt(19,14,2) = -0.34500e-02_wp
            htt(19,15,2) = 0.19200e-02_wp
            htt(19,16,2) = -0.27800e-02_wp
            htt(19,17,2) = -0.95000e-03_wp
            htt(19,18,2) = -0.21800e-02_wp
            htt(19,19,2) = 0.37700e-02_wp
            htt(20,1,2)  = -0.35000e-03_wp
            htt(20,2,2)  = -0.26000e-03_wp
            htt(20,3,2)  = -0.16000e-03_wp
            htt(20,4,2)  = -0.80000e-04_wp
            htt(20,5,2)  = 0.40000e-04_wp
            htt(20,6,2)  = -0.40000e-03_wp
            htt(20,7,2)  = -0.14000e-03_wp
            htt(20,8,2)  = -0.57000e-03_wp
            htt(20,9,2)  = -0.93000e-03_wp
            htt(20,10,2) = -0.42000e-03_wp
            htt(20,11,2) = -0.10200e-02_wp
            htt(20,12,2) = -0.51000e-03_wp
            htt(20,13,2) = 0.26100e-02_wp
            htt(20,14,2) = -0.28400e-02_wp
            htt(20,15,2) = 0.21100e-02_wp
            htt(20,16,2) = -0.33500e-02_wp
            htt(20,17,2) = -0.29400e-02_wp
            htt(20,18,2) = 0.47000e-03_wp
            htt(20,19,2) = 0.25400e-02_wp
            htt(20,20,2) = 0.59400e-02_wp
            htt(21,1,2)  = 0.21700e-02_wp
            htt(21,2,2)  = 0.15000e-02_wp
            htt(21,3,2)  = -0.99000e-03_wp
            htt(21,4,2)  = -0.16600e-02_wp
            htt(21,5,2)  = -0.30200e-02_wp
            htt(21,6,2)  = -0.69000e-03_wp
            htt(21,7,2)  = -0.20000e-04_wp
            htt(21,8,2)  = 0.23100e-02_wp
            htt(21,9,2)  = -0.11900e-02_wp
            htt(21,10,2) = -0.98000e-03_wp
            htt(21,11,2) = -0.50000e-03_wp
            htt(21,12,2) = -0.28000e-03_wp
            htt(21,13,2) = 0.12600e-02_wp
            htt(21,14,2) = -0.76000e-03_wp
            htt(21,15,2) = 0.10500e-02_wp
            htt(21,16,2) = -0.97000e-03_wp
            htt(21,17,2) = -0.21400e-02_wp
            htt(21,18,2) = -0.48000e-03_wp
            htt(21,19,2) = -0.11000e-03_wp
            htt(21,20,2) = 0.15500e-02_wp
            htt(21,21,2) = 0.28900e-02_wp
            htt(22,1,2)  = 0.27000e-03_wp
            htt(22,2,2)  = 0.11300e-02_wp
            htt(22,3,2)  = -0.15900e-02_wp
            htt(22,4,2)  = -0.74000e-03_wp
            htt(22,5,2)  = -0.60000e-04_wp
            htt(22,6,2)  = -0.26200e-02_wp
            htt(22,7,2)  = 0.17100e-02_wp
            htt(22,8,2)  = -0.85000e-03_wp
            htt(22,9,2)  = -0.85000e-03_wp
            htt(22,10,2) = 0.40000e-04_wp
            htt(22,11,2) = -0.17300e-02_wp
            htt(22,12,2) = -0.85000e-03_wp
            htt(22,13,2) = 0.12800e-02_wp
            htt(22,14,2) = -0.15000e-03_wp
            htt(22,15,2) = 0.40000e-03_wp
            htt(22,16,2) = -0.10300e-02_wp
            htt(22,17,2) = -0.10700e-02_wp
            htt(22,18,2) = 0.56000e-03_wp
            htt(22,19,2) = 0.37000e-03_wp
            htt(22,20,2) = 0.20000e-02_wp
            htt(22,21,2) = 0.44000e-03_wp
            htt(22,22,2) = 0.31200e-02_wp
            htt(23,1,2)  = 0.19200e-02_wp
            htt(23,2,2)  = 0.79000e-03_wp
            htt(23,3,2)  = -0.23000e-03_wp
            htt(23,4,2)  = -0.13700e-02_wp
            htt(23,5,2)  = -0.21700e-02_wp
            htt(23,6,2)  = 0.14000e-03_wp
            htt(23,7,2)  = -0.12000e-03_wp
            htt(23,8,2)  = 0.21900e-02_wp
            htt(23,9,2)  = -0.74000e-03_wp
            htt(23,10,2) = -0.17100e-02_wp
            htt(23,11,2) = 0.45000e-03_wp
            htt(23,12,2) = -0.53000e-03_wp
            htt(23,13,2) = -0.10000e-03_wp
            htt(23,14,2) = -0.13700e-02_wp
            htt(23,15,2) = 0.87000e-03_wp
            htt(23,16,2) = -0.40000e-03_wp
            htt(23,17,2) = -0.15000e-03_wp
            htt(23,18,2) = -0.31500e-02_wp
            htt(23,19,2) = 0.11200e-02_wp
            htt(23,20,2) = -0.18800e-02_wp
            htt(23,21,2) = 0.12200e-02_wp
            htt(23,22,2) = -0.12000e-02_wp
            htt(23,23,2) = 0.42400e-02_wp
            htt(24,1,2)  = 0.30000e-04_wp
            htt(24,2,2)  = 0.42000e-03_wp
            htt(24,3,2)  = -0.84000e-03_wp
            htt(24,4,2)  = -0.45000e-03_wp
            htt(24,5,2)  = 0.79000e-03_wp
            htt(24,6,2)  = -0.18000e-02_wp
            htt(24,7,2)  = 0.16200e-02_wp
            htt(24,8,2)  = -0.97000e-03_wp
            htt(24,9,2)  = -0.39000e-03_wp
            htt(24,10,2) = -0.70000e-03_wp
            htt(24,11,2) = -0.79000e-03_wp
            htt(24,12,2) = -0.11000e-02_wp
            htt(24,13,2) = -0.80000e-04_wp
            htt(24,14,2) = -0.76000e-03_wp
            htt(24,15,2) = 0.22000e-03_wp
            htt(24,16,2) = -0.45000e-03_wp
            htt(24,17,2) = 0.92000e-03_wp
            htt(24,18,2) = -0.21000e-02_wp
            htt(24,19,2) = 0.16000e-02_wp
            htt(24,20,2) = -0.14300e-02_wp
            htt(24,21,2) = -0.12300e-02_wp
            htt(24,22,2) = 0.14800e-02_wp
            htt(24,23,2) = 0.18100e-02_wp
            htt(24,24,2) = 0.45300e-02_wp
            ! Hessian elements for anchor point 2
            hbs( 1,1,2)  = -0.31760e-01_wp
            hbs( 1,2,2)  = -0.29360e-01_wp
            hbs( 1,3,2)  = -0.31820e-01_wp
            hbs( 1,4,2)  = -0.95900e-02_wp
            hbs( 1,5,2)  = 0.52300e-02_wp
            hbs( 1,6,2)  = 0.39930e-01_wp
            hbs( 1,7,2)  = -0.45000e-03_wp
            hbs( 1,8,2)  = 0.40450e-01_wp
            hbs( 1,9,2)  = -0.18000e-03_wp
            hbs( 1,10,2) = -0.96900e-02_wp
            hbs( 1,11,2) = -0.42000e-03_wp
            hbs( 1,12,2) = 0.47400e-02_wp
            hbs( 2,1,2)  = 0.49930e-01_wp
            hbs( 2,2,2)  = -0.13930e-01_wp
            hbs( 2,3,2)  = 0.99600e-02_wp
            hbs( 2,4,2)  = 0.18000e-02_wp
            hbs( 2,5,2)  = -0.12200e-02_wp
            hbs( 2,6,2)  = -0.27810e-01_wp
            hbs( 2,7,2)  = 0.11800e-02_wp
            hbs( 2,8,2)  = -0.14650e-01_wp
            hbs( 2,9,2)  = 0.00000e+00_wp
            hbs( 2,10,2) = 0.41300e-02_wp
            hbs( 2,11,2) = -0.69000e-03_wp
            hbs( 2,12,2) = -0.89000e-03_wp
            hbs( 3,1,2)  = -0.18180e-01_wp
            hbs( 3,2,2)  = 0.43280e-01_wp
            hbs( 3,3,2)  = 0.21860e-01_wp
            hbs( 3,4,2)  = 0.77900e-02_wp
            hbs( 3,5,2)  = -0.40100e-02_wp
            hbs( 3,6,2)  = -0.12120e-01_wp
            hbs( 3,7,2)  = -0.73000e-03_wp
            hbs( 3,8,2)  = -0.25800e-01_wp
            hbs( 3,9,2)  = 0.18000e-03_wp
            hbs( 3,10,2) = 0.55600e-02_wp
            hbs( 3,11,2) = 0.11100e-02_wp
            hbs( 3,12,2) = -0.38500e-02_wp
            hbs( 4,1,2)  = -0.20820e-01_wp
            hbs( 4,2,2)  = -0.12100e-01_wp
            hbs( 4,3,2)  = 0.15180e-01_wp
            hbs( 4,4,2)  = -0.29840e-01_wp
            hbs( 4,5,2)  = -0.97300e-02_wp
            hbs( 4,6,2)  = -0.57000e-02_wp
            hbs( 4,7,2)  = 0.49500e-02_wp
            hbs( 4,8,2)  = 0.31150e-01_wp
            hbs( 4,9,2)  = -0.44000e-03_wp
            hbs( 4,10,2) = 0.40270e-01_wp
            hbs( 4,11,2) = -0.59000e-03_wp
            hbs( 4,12,2) = 0.29000e-03_wp
            hbs( 5,1,2)  = 0.26090e-01_wp
            hbs( 5,2,2)  = 0.58900e-02_wp
            hbs( 5,3,2)  = -0.14520e-01_wp
            hbs( 5,4,2)  = 0.25400e-02_wp
            hbs( 5,5,2)  = 0.69400e-02_wp
            hbs( 5,6,2)  = 0.30900e-02_wp
            hbs( 5,7,2)  = -0.21700e-02_wp
            hbs( 5,8,2)  = -0.15190e-01_wp
            hbs( 5,9,2)  = -0.39000e-03_wp
            hbs( 5,10,2) = -0.22140e-01_wp
            hbs( 5,11,2) = 0.28000e-03_wp
            hbs( 5,12,2) = 0.90000e-04_wp
            hbs( 6,1,2)  = -0.52700e-02_wp
            hbs( 6,2,2)  = 0.62100e-02_wp
            hbs( 6,3,2)  = -0.66000e-03_wp
            hbs( 6,4,2)  = 0.27300e-01_wp
            hbs( 6,5,2)  = 0.27900e-02_wp
            hbs( 6,6,2)  = 0.26100e-02_wp
            hbs( 6,7,2)  = -0.27800e-02_wp
            hbs( 6,8,2)  = -0.15960e-01_wp
            hbs( 6,9,2)  = 0.83000e-03_wp
            hbs( 6,10,2) = -0.18130e-01_wp
            hbs( 6,11,2) = 0.31000e-03_wp
            hbs( 6,12,2) = -0.38000e-03_wp
            hbs( 7,1,2)  = -0.90800e-02_wp
            hbs( 7,2,2)  = 0.35870e-01_wp
            hbs( 7,3,2)  = -0.10000e-04_wp
            hbs( 7,4,2)  = -0.22920e-01_wp
            hbs( 7,5,2)  = 0.61900e-02_wp
            hbs( 7,6,2)  = -0.28970e-01_wp
            hbs( 7,7,2)  = -0.87500e-02_wp
            hbs( 7,8,2)  = -0.92900e-02_wp
            hbs( 7,9,2)  = 0.49700e-02_wp
            hbs( 7,10,2) = 0.32780e-01_wp
            hbs( 7,11,2) = 0.60000e-04_wp
            hbs( 7,12,2) = -0.10700e-02_wp
            hbs( 8,1,2)  = 0.36900e-02_wp
            hbs( 8,2,2)  = -0.18590e-01_wp
            hbs( 8,3,2)  = 0.13800e-02_wp
            hbs( 8,4,2)  = 0.25350e-01_wp
            hbs( 8,5,2)  = -0.35000e-02_wp
            hbs( 8,6,2)  = -0.42000e-03_wp
            hbs( 8,7,2)  = 0.45000e-02_wp
            hbs( 8,8,2)  = 0.50400e-02_wp
            hbs( 8,9,2)  = -0.20100e-02_wp
            hbs( 8,10,2) = -0.14680e-01_wp
            hbs( 8,11,2) = -0.62000e-03_wp
            hbs( 8,12,2) = 0.55000e-03_wp
            hbs( 9,1,2)  = 0.53900e-02_wp
            hbs( 9,2,2)  = -0.17270e-01_wp
            hbs( 9,3,2)  = -0.13700e-02_wp
            hbs( 9,4,2)  = -0.24400e-02_wp
            hbs( 9,5,2)  = -0.26800e-02_wp
            hbs( 9,6,2)  = 0.29390e-01_wp
            hbs( 9,7,2)  = 0.42400e-02_wp
            hbs( 9,8,2)  = 0.42500e-02_wp
            hbs( 9,9,2)  = -0.29600e-02_wp
            hbs( 9,10,2) = -0.18100e-01_wp
            hbs( 9,11,2) = 0.55000e-03_wp
            hbs( 9,12,2) = 0.51000e-03_wp
            hbs(10,1,2)  = 0.36130e-01_wp
            hbs(10,2,2)  = 0.38360e-01_wp
            hbs(10,3,2)  = 0.36300e-02_wp
            hbs(10,4,2)  = -0.11890e-01_wp
            hbs(10,5,2)  = -0.10500e-02_wp
            hbs(10,6,2)  = -0.26840e-01_wp
            hbs(10,7,2)  = 0.48100e-02_wp
            hbs(10,8,2)  = -0.27100e-01_wp
            hbs(10,9,2)  = -0.88200e-02_wp
            hbs(10,10,2) = -0.12140e-01_wp
            hbs(10,11,2) = 0.47600e-02_wp
            hbs(10,12,2) = -0.77000e-03_wp
            hbs(11,1,2)  = -0.20190e-01_wp
            hbs(11,2,2)  = -0.17180e-01_wp
            hbs(11,3,2)  = -0.21300e-02_wp
            hbs(11,4,2)  = 0.57400e-02_wp
            hbs(11,5,2)  = 0.12900e-02_wp
            hbs(11,6,2)  = 0.27000e-01_wp
            hbs(11,7,2)  = -0.28400e-02_wp
            hbs(11,8,2)  = 0.47000e-03_wp
            hbs(11,9,2)  = 0.44300e-02_wp
            hbs(11,10,2) = 0.60700e-02_wp
            hbs(11,11,2) = -0.19800e-02_wp
            hbs(11,12,2) = -0.18000e-03_wp
            hbs(12,1,2)  = -0.15940e-01_wp
            hbs(12,2,2)  = -0.21180e-01_wp
            hbs(12,3,2)  = -0.14900e-02_wp
            hbs(12,4,2)  = 0.61500e-02_wp
            hbs(12,5,2)  = -0.24000e-03_wp
            hbs(12,6,2)  = -0.16000e-03_wp
            hbs(12,7,2)  = -0.19600e-02_wp
            hbs(12,8,2)  = 0.26630e-01_wp
            hbs(12,9,2)  = 0.43900e-02_wp
            hbs(12,10,2) = 0.60700e-02_wp
            hbs(12,11,2) = -0.27800e-02_wp
            hbs(12,12,2) = 0.94000e-03_wp
            hbs(13,1,2)  = 0.35720e-01_wp
            hbs(13,2,2)  = -0.96900e-02_wp
            hbs(13,3,2)  = 0.29800e-02_wp
            hbs(13,4,2)  = 0.32310e-01_wp
            hbs(13,5,2)  = -0.94000e-03_wp
            hbs(13,6,2)  = -0.10140e-01_wp
            hbs(13,7,2)  = 0.70000e-04_wp
            hbs(13,8,2)  = -0.29530e-01_wp
            hbs(13,9,2)  = 0.50200e-02_wp
            hbs(13,10,2) = -0.23420e-01_wp
            hbs(13,11,2) = -0.86600e-02_wp
            hbs(13,12,2) = 0.53400e-02_wp
            hbs(14,1,2)  = -0.17210e-01_wp
            hbs(14,2,2)  = 0.56100e-02_wp
            hbs(14,3,2)  = -0.23300e-02_wp
            hbs(14,4,2)  = -0.17680e-01_wp
            hbs(14,5,2)  = 0.48000e-03_wp
            hbs(14,6,2)  = 0.47700e-02_wp
            hbs(14,7,2)  = 0.56000e-03_wp
            hbs(14,8,2)  = 0.29630e-01_wp
            hbs(14,9,2)  = -0.29700e-02_wp
            hbs(14,10,2) = -0.24000e-02_wp
            hbs(14,11,2) = 0.43000e-02_wp
            hbs(14,12,2) = -0.22600e-02_wp
            hbs(15,1,2)  = -0.18510e-01_wp
            hbs(15,2,2)  = 0.40800e-02_wp
            hbs(15,3,2)  = -0.65000e-03_wp
            hbs(15,4,2)  = -0.14630e-01_wp
            hbs(15,5,2)  = 0.46000e-03_wp
            hbs(15,6,2)  = 0.53700e-02_wp
            hbs(15,7,2)  = -0.63000e-03_wp
            hbs(15,8,2)  = -0.10000e-03_wp
            hbs(15,9,2)  = -0.20500e-02_wp
            hbs(15,10,2) = 0.25820e-01_wp
            hbs(15,11,2) = 0.43700e-02_wp
            hbs(15,12,2) = -0.30900e-02_wp
            hbs(16,1,2)  = -0.10200e-01_wp
            hbs(16,2,2)  = -0.23080e-01_wp
            hbs(16,3,2)  = 0.10040e-01_wp
            hbs(16,4,2)  = 0.41930e-01_wp
            hbs(16,5,2)  = 0.31000e-03_wp
            hbs(16,6,2)  = 0.31710e-01_wp
            hbs(16,7,2)  = -0.63000e-03_wp
            hbs(16,8,2)  = -0.56800e-02_wp
            hbs(16,9,2)  = -0.55000e-03_wp
            hbs(16,10,2) = -0.27810e-01_wp
            hbs(16,11,2) = 0.48500e-02_wp
            hbs(16,12,2) = -0.85400e-02_wp
            hbs(17,1,2)  = 0.41100e-02_wp
            hbs(17,2,2)  = 0.25000e-01_wp
            hbs(17,3,2)  = -0.74400e-02_wp
            hbs(17,4,2)  = -0.22980e-01_wp
            hbs(17,5,2)  = 0.60000e-03_wp
            hbs(17,6,2)  = -0.14860e-01_wp
            hbs(17,7,2)  = 0.32000e-03_wp
            hbs(17,8,2)  = 0.33600e-02_wp
            hbs(17,9,2)  = -0.32000e-03_wp
            hbs(17,10,2) = 0.20400e-02_wp
            hbs(17,11,2) = -0.21000e-02_wp
            hbs(17,12,2) = 0.68800e-02_wp
            hbs(18,1,2)  = 0.61000e-02_wp
            hbs(18,2,2)  = -0.19200e-02_wp
            hbs(18,3,2)  = -0.26000e-02_wp
            hbs(18,4,2)  = -0.18950e-01_wp
            hbs(18,5,2)  = -0.92000e-03_wp
            hbs(18,6,2)  = -0.16850e-01_wp
            hbs(18,7,2)  = 0.31000e-03_wp
            hbs(18,8,2)  = 0.23200e-02_wp
            hbs(18,9,2)  = 0.87000e-03_wp
            hbs(18,10,2) = 0.25770e-01_wp
            hbs(18,11,2) = -0.27500e-02_wp
            hbs(18,12,2) = 0.16600e-02_wp
            ! Hessian elements for anchor point 3
            hss( 1,1,3)  = 0.34759e+00_wp
            hss( 2,1,3)  = 0.17810e-01_wp
            hss( 2,2,3)  = 0.35623e+00_wp
            hss( 3,1,3)  = 0.29010e-01_wp
            hss( 3,2,3)  = 0.32710e-01_wp
            hss( 3,3,3)  = 0.46124e+00_wp
            hss( 4,1,3)  = -0.19600e-02_wp
            hss( 4,2,3)  = 0.14870e-01_wp
            hss( 4,3,3)  = -0.89100e-02_wp
            hss( 4,4,3)  = 0.38581e+00_wp
            hss( 5,1,3)  = 0.38800e-02_wp
            hss( 5,2,3)  = -0.77000e-03_wp
            hss( 5,3,3)  = 0.71700e-02_wp
            hss( 5,4,3)  = 0.36500e-02_wp
            hss( 5,5,3)  = 0.35807e+00_wp
            hss( 6,1,3)  = 0.33890e-01_wp
            hss( 6,2,3)  = 0.44570e-01_wp
            hss( 6,3,3)  = -0.43100e-02_wp
            hss( 6,4,3)  = -0.36100e-02_wp
            hss( 6,5,3)  = 0.55000e-03_wp
            hss( 6,6,3)  = 0.35931e+00_wp
            hss( 7,1,3)  = -0.54000e-03_wp
            hss( 7,2,3)  = -0.63000e-03_wp
            hss( 7,3,3)  = 0.57000e-03_wp
            hss( 7,4,3)  = 0.39700e-02_wp
            hss( 7,5,3)  = 0.39000e-03_wp
            hss( 7,6,3)  = 0.39800e-02_wp
            hss( 7,7,3)  = 0.37091e+00_wp
            hss( 8,1,3)  = 0.42030e-01_wp
            hss( 8,2,3)  = 0.35060e-01_wp
            hss( 8,3,3)  = -0.52800e-02_wp
            hss( 8,4,3)  = 0.24280e-01_wp
            hss( 8,5,3)  = -0.11100e-02_wp
            hss( 8,6,3)  = 0.17880e-01_wp
            hss( 8,7,3)  = -0.10400e-02_wp
            hss( 8,8,3)  = 0.35888e+00_wp
            hss( 9,1,3)  = -0.53000e-03_wp
            hss( 9,2,3)  = -0.57000e-03_wp
            hss( 9,3,3)  = 0.26000e-03_wp
            hss( 9,4,3)  = -0.39000e-03_wp
            hss( 9,5,3)  = 0.19000e-03_wp
            hss( 9,6,3)  = 0.40100e-02_wp
            hss( 9,7,3)  = 0.49000e-03_wp
            hss( 9,8,3)  = 0.39600e-02_wp
            hss( 9,9,3)  = 0.36356e+00_wp
            hss(10,1,3)  = 0.19520e-01_wp
            hss(10,2,3)  = -0.69300e-02_wp
            hss(10,3,3)  = -0.48800e-02_wp
            hss(10,4,3)  = 0.64010e-01_wp
            hss(10,5,3)  = -0.43000e-03_wp
            hss(10,6,3)  = 0.21230e-01_wp
            hss(10,7,3)  = -0.10900e-02_wp
            hss(10,8,3)  = -0.77000e-03_wp
            hss(10,9,3)  = -0.19000e-03_wp
            hss(10,10,3) = 0.38116e+00_wp
            hss(11,1,3)  = -0.69000e-03_wp
            hss(11,2,3)  = -0.50000e-03_wp
            hss(11,3,3)  = 0.52000e-03_wp
            hss(11,4,3)  = -0.11600e-02_wp
            hss(11,5,3)  = 0.80000e-04_wp
            hss(11,6,3)  = -0.10000e-02_wp
            hss(11,7,3)  = 0.22000e-03_wp
            hss(11,8,3)  = 0.39700e-02_wp
            hss(11,9,3)  = 0.48000e-03_wp
            hss(11,10,3) = 0.39800e-02_wp
            hss(11,11,3) = 0.37152e+00_wp
            hss(12,1,3)  = -0.16500e-02_wp
            hss(12,2,3)  = 0.39400e-02_wp
            hss(12,3,3)  = -0.19000e-03_wp
            hss(12,4,3)  = -0.10000e-03_wp
            hss(12,5,3)  = 0.35000e-03_wp
            hss(12,6,3)  = -0.10400e-02_wp
            hss(12,7,3)  = 0.30000e-04_wp
            hss(12,8,3)  = 0.57000e-03_wp
            hss(12,9,3)  = 0.12000e-03_wp
            hss(12,10,3) = 0.23600e-02_wp
            hss(12,11,3) = 0.40000e-03_wp
            hss(12,12,3) = 0.37254e+00_wp
            ! Hessian elements for anchor point 3
            hbb( 1,1,3)  = 0.82450e-01_wp
            hbb( 2,1,3)  = -0.46090e-01_wp
            hbb( 2,2,3)  = 0.14487e+00_wp
            hbb( 3,1,3)  = -0.36360e-01_wp
            hbb( 3,2,3)  = -0.98780e-01_wp
            hbb( 3,3,3)  = 0.13515e+00_wp
            hbb( 4,1,3)  = -0.36000e-01_wp
            hbb( 4,2,3)  = 0.23310e-01_wp
            hbb( 4,3,3)  = 0.12690e-01_wp
            hbb( 4,4,3)  = 0.66170e-01_wp
            hbb( 5,1,3)  = 0.20240e-01_wp
            hbb( 5,2,3)  = -0.12480e-01_wp
            hbb( 5,3,3)  = -0.77700e-02_wp
            hbb( 5,4,3)  = -0.33030e-01_wp
            hbb( 5,5,3)  = 0.71950e-01_wp
            hbb( 6,1,3)  = 0.15760e-01_wp
            hbb( 6,2,3)  = -0.10840e-01_wp
            hbb( 6,3,3)  = -0.49200e-02_wp
            hbb( 6,4,3)  = -0.33140e-01_wp
            hbb( 6,5,3)  = -0.38920e-01_wp
            hbb( 6,6,3)  = 0.72060e-01_wp
            hbb( 7,1,3)  = -0.71500e-02_wp
            hbb( 7,2,3)  = -0.39000e-02_wp
            hbb( 7,3,3)  = 0.11050e-01_wp
            hbb( 7,4,3)  = -0.25740e-01_wp
            hbb( 7,5,3)  = 0.93100e-02_wp
            hbb( 7,6,3)  = 0.16430e-01_wp
            hbb( 7,7,3)  = 0.64420e-01_wp
            hbb( 8,1,3)  = -0.22800e-02_wp
            hbb( 8,2,3)  = 0.74300e-02_wp
            hbb( 8,3,3)  = -0.51500e-02_wp
            hbb( 8,4,3)  = 0.15690e-01_wp
            hbb( 8,5,3)  = -0.49500e-02_wp
            hbb( 8,6,3)  = -0.10740e-01_wp
            hbb( 8,7,3)  = -0.31490e-01_wp
            hbb( 8,8,3)  = 0.72950e-01_wp
            hbb( 9,1,3)  = 0.94300e-02_wp
            hbb( 9,2,3)  = -0.35300e-02_wp
            hbb( 9,3,3)  = -0.59100e-02_wp
            hbb( 9,4,3)  = 0.10050e-01_wp
            hbb( 9,5,3)  = -0.43500e-02_wp
            hbb( 9,6,3)  = -0.56900e-02_wp
            hbb( 9,7,3)  = -0.32930e-01_wp
            hbb( 9,8,3)  = -0.41460e-01_wp
            hbb( 9,9,3)  = 0.74390e-01_wp
            hbb(10,1,3)  = 0.45500e-02_wp
            hbb(10,2,3)  = -0.20200e-02_wp
            hbb(10,3,3)  = -0.25300e-02_wp
            hbb(10,4,3)  = -0.64100e-02_wp
            hbb(10,5,3)  = 0.91900e-02_wp
            hbb(10,6,3)  = -0.27700e-02_wp
            hbb(10,7,3)  = -0.33410e-01_wp
            hbb(10,8,3)  = 0.13410e-01_wp
            hbb(10,9,3)  = 0.20000e-01_wp
            hbb(10,10,3) = 0.75440e-01_wp
            hbb(11,1,3)  = -0.22100e-02_wp
            hbb(11,2,3)  = 0.12300e-02_wp
            hbb(11,3,3)  = 0.98000e-03_wp
            hbb(11,4,3)  = -0.23100e-02_wp
            hbb(11,5,3)  = -0.32200e-02_wp
            hbb(11,6,3)  = 0.55300e-02_wp
            hbb(11,7,3)  = 0.20020e-01_wp
            hbb(11,8,3)  = -0.74100e-02_wp
            hbb(11,9,3)  = -0.12610e-01_wp
            hbb(11,10,3) = -0.37790e-01_wp
            hbb(11,11,3) = 0.78360e-01_wp
            hbb(12,1,3)  = -0.23400e-02_wp
            hbb(12,2,3)  = 0.79000e-03_wp
            hbb(12,3,3)  = 0.15500e-02_wp
            hbb(12,4,3)  = 0.87200e-02_wp
            hbb(12,5,3)  = -0.59700e-02_wp
            hbb(12,6,3)  = -0.27500e-02_wp
            hbb(12,7,3)  = 0.13390e-01_wp
            hbb(12,8,3)  = -0.60000e-02_wp
            hbb(12,9,3)  = -0.73900e-02_wp
            hbb(12,10,3) = -0.37640e-01_wp
            hbb(12,11,3) = -0.40570e-01_wp
            hbb(12,12,3) = 0.78220e-01_wp
            hbb(13,1,3)  = -0.67600e-02_wp
            hbb(13,2,3)  = 0.10710e-01_wp
            hbb(13,3,3)  = -0.39500e-02_wp
            hbb(13,4,3)  = -0.11500e-02_wp
            hbb(13,5,3)  = 0.12900e-02_wp
            hbb(13,6,3)  = -0.14000e-03_wp
            hbb(13,7,3)  = 0.29100e-02_wp
            hbb(13,8,3)  = 0.41300e-02_wp
            hbb(13,9,3)  = -0.70300e-02_wp
            hbb(13,10,3) = -0.33760e-01_wp
            hbb(13,11,3) = 0.13510e-01_wp
            hbb(13,12,3) = 0.20240e-01_wp
            hbb(13,13,3) = 0.65300e-01_wp
            hbb(14,1,3)  = 0.91200e-02_wp
            hbb(14,2,3)  = -0.60500e-02_wp
            hbb(14,3,3)  = -0.30700e-02_wp
            hbb(14,4,3)  = 0.66000e-03_wp
            hbb(14,5,3)  = -0.85000e-03_wp
            hbb(14,6,3)  = 0.19000e-03_wp
            hbb(14,7,3)  = -0.71100e-02_wp
            hbb(14,8,3)  = -0.51000e-03_wp
            hbb(14,9,3)  = 0.76200e-02_wp
            hbb(14,10,3) = 0.20030e-01_wp
            hbb(14,11,3) = -0.73800e-02_wp
            hbb(14,12,3) = -0.12650e-01_wp
            hbb(14,13,3) = -0.33180e-01_wp
            hbb(14,14,3) = 0.74350e-01_wp
            hbb(15,1,3)  = -0.23600e-02_wp
            hbb(15,2,3)  = -0.46500e-02_wp
            hbb(15,3,3)  = 0.70200e-02_wp
            hbb(15,4,3)  = 0.49000e-03_wp
            hbb(15,5,3)  = -0.44000e-03_wp
            hbb(15,6,3)  = -0.50000e-04_wp
            hbb(15,7,3)  = 0.42000e-02_wp
            hbb(15,8,3)  = -0.36200e-02_wp
            hbb(15,9,3)  = -0.59000e-03_wp
            hbb(15,10,3) = 0.13730e-01_wp
            hbb(15,11,3) = -0.61300e-02_wp
            hbb(15,12,3) = -0.76000e-02_wp
            hbb(15,13,3) = -0.32120e-01_wp
            hbb(15,14,3) = -0.41170e-01_wp
            hbb(15,15,3) = 0.73290e-01_wp
            hbb(16,1,3)  = -0.37090e-01_wp
            hbb(16,2,3)  = 0.17990e-01_wp
            hbb(16,3,3)  = 0.19100e-01_wp
            hbb(16,4,3)  = 0.31400e-02_wp
            hbb(16,5,3)  = -0.70000e-02_wp
            hbb(16,6,3)  = 0.38600e-02_wp
            hbb(16,7,3)  = -0.10300e-02_wp
            hbb(16,8,3)  = 0.54000e-03_wp
            hbb(16,9,3)  = 0.48000e-03_wp
            hbb(16,10,3) = -0.64000e-02_wp
            hbb(16,11,3) = 0.87800e-02_wp
            hbb(16,12,3) = -0.23700e-02_wp
            hbb(16,13,3) = -0.26530e-01_wp
            hbb(16,14,3) = 0.10480e-01_wp
            hbb(16,15,3) = 0.16060e-01_wp
            hbb(16,16,3) = 0.67910e-01_wp
            hbb(17,1,3)  = 0.22610e-01_wp
            hbb(17,2,3)  = -0.90200e-02_wp
            hbb(17,3,3)  = -0.13590e-01_wp
            hbb(17,4,3)  = -0.76000e-02_wp
            hbb(17,5,3)  = 0.77900e-02_wp
            hbb(17,6,3)  = -0.19000e-03_wp
            hbb(17,7,3)  = 0.87000e-03_wp
            hbb(17,8,3)  = -0.25000e-03_wp
            hbb(17,9,3)  = -0.62000e-03_wp
            hbb(17,10,3) = 0.91600e-02_wp
            hbb(17,11,3) = -0.60700e-02_wp
            hbb(17,12,3) = -0.31000e-02_wp
            hbb(17,13,3) = 0.96800e-02_wp
            hbb(17,14,3) = -0.46000e-02_wp
            hbb(17,15,3) = -0.50800e-02_wp
            hbb(17,16,3) = -0.34720e-01_wp
            hbb(17,17,3) = 0.72390e-01_wp
            hbb(18,1,3)  = 0.14470e-01_wp
            hbb(18,2,3)  = -0.89700e-02_wp
            hbb(18,3,3)  = -0.55000e-02_wp
            hbb(18,4,3)  = 0.44600e-02_wp
            hbb(18,5,3)  = -0.79000e-03_wp
            hbb(18,6,3)  = -0.36700e-02_wp
            hbb(18,7,3)  = 0.16000e-03_wp
            hbb(18,8,3)  = -0.29000e-03_wp
            hbb(18,9,3)  = 0.13000e-03_wp
            hbb(18,10,3) = -0.27600e-02_wp
            hbb(18,11,3) = -0.27100e-02_wp
            hbb(18,12,3) = 0.54700e-02_wp
            hbb(18,13,3) = 0.16860e-01_wp
            hbb(18,14,3) = -0.58800e-02_wp
            hbb(18,15,3) = -0.10980e-01_wp
            hbb(18,16,3) = -0.33190e-01_wp
            hbb(18,17,3) = -0.37670e-01_wp
            hbb(18,18,3) = 0.70860e-01_wp
            ! Hessian elements for anchor point 3
            htt( 1,1,3)  = 0.60400e-02_wp
            htt( 2,1,3)  = 0.38400e-02_wp
            htt( 2,2,3)  = 0.50900e-02_wp
            htt( 3,1,3)  = -0.18500e-02_wp
            htt( 3,2,3)  = -0.33100e-02_wp
            htt( 3,3,3)  = 0.75100e-02_wp
            htt( 4,1,3)  = -0.40500e-02_wp
            htt( 4,2,3)  = -0.20700e-02_wp
            htt( 4,3,3)  = 0.60500e-02_wp
            htt( 4,4,3)  = 0.80200e-02_wp
            htt( 5,1,3)  = -0.59200e-02_wp
            htt( 5,2,3)  = -0.42600e-02_wp
            htt( 5,3,3)  = 0.21200e-02_wp
            htt( 5,4,3)  = 0.37800e-02_wp
            htt( 5,5,3)  = 0.65600e-02_wp
            htt( 6,1,3)  = -0.41200e-02_wp
            htt( 6,2,3)  = -0.39000e-02_wp
            htt( 6,3,3)  = 0.26900e-02_wp
            htt( 6,4,3)  = 0.29100e-02_wp
            htt( 6,5,3)  = 0.37400e-02_wp
            htt( 6,6,3)  = 0.55900e-02_wp
            htt( 7,1,3)  = 0.15700e-02_wp
            htt( 7,2,3)  = 0.25400e-02_wp
            htt( 7,3,3)  = -0.67700e-02_wp
            htt( 7,4,3)  = -0.58000e-02_wp
            htt( 7,5,3)  = -0.10700e-02_wp
            htt( 7,6,3)  = -0.27200e-02_wp
            htt( 7,7,3)  = 0.68600e-02_wp
            htt( 8,1,3)  = 0.33800e-02_wp
            htt( 8,2,3)  = 0.29000e-02_wp
            htt( 8,3,3)  = -0.62000e-02_wp
            htt( 8,4,3)  = -0.66800e-02_wp
            htt( 8,5,3)  = -0.38900e-02_wp
            htt( 8,6,3)  = -0.88000e-03_wp
            htt( 8,7,3)  = 0.52100e-02_wp
            htt( 8,8,3)  = 0.82200e-02_wp
            htt( 9,1,3)  = -0.23900e-02_wp
            htt( 9,2,3)  = -0.68000e-03_wp
            htt( 9,3,3)  = 0.45000e-03_wp
            htt( 9,4,3)  = 0.21600e-02_wp
            htt( 9,5,3)  = 0.17400e-02_wp
            htt( 9,6,3)  = 0.14100e-02_wp
            htt( 9,7,3)  = -0.96000e-03_wp
            htt( 9,8,3)  = -0.12900e-02_wp
            htt( 9,9,3)  = 0.24600e-02_wp
            htt(10,1,3)  = -0.15600e-02_wp
            htt(10,2,3)  = 0.23000e-03_wp
            htt(10,3,3)  = 0.34000e-03_wp
            htt(10,4,3)  = 0.21300e-02_wp
            htt(10,5,3)  = 0.17600e-02_wp
            htt(10,6,3)  = 0.79000e-03_wp
            htt(10,7,3)  = -0.40000e-04_wp
            htt(10,8,3)  = -0.10100e-02_wp
            htt(10,9,3)  = 0.57000e-03_wp
            htt(10,10,3) = 0.39100e-02_wp
            htt(11,1,3)  = -0.10000e-03_wp
            htt(11,2,3)  = -0.19700e-02_wp
            htt(11,3,3)  = 0.19800e-02_wp
            htt(11,4,3)  = 0.11000e-03_wp
            htt(11,5,3)  = 0.10000e-04_wp
            htt(11,6,3)  = 0.11900e-02_wp
            htt(11,7,3)  = -0.19700e-02_wp
            htt(11,8,3)  = -0.79000e-03_wp
            htt(11,9,3)  = 0.68000e-03_wp
            htt(11,10,3) = -0.12900e-02_wp
            htt(11,11,3) = 0.26200e-02_wp
            htt(12,1,3)  = 0.72000e-03_wp
            htt(12,2,3)  = -0.10500e-02_wp
            htt(12,3,3)  = 0.18600e-02_wp
            htt(12,4,3)  = 0.90000e-04_wp
            htt(12,5,3)  = 0.30000e-04_wp
            htt(12,6,3)  = 0.57000e-03_wp
            htt(12,7,3)  = -0.10500e-02_wp
            htt(12,8,3)  = -0.51000e-03_wp
            htt(12,9,3)  = -0.12000e-02_wp
            htt(12,10,3) = 0.20500e-02_wp
            htt(12,11,3) = 0.64000e-03_wp
            htt(12,12,3) = 0.39000e-02_wp
            htt(13,1,3)  = -0.10900e-02_wp
            htt(13,2,3)  = -0.19400e-02_wp
            htt(13,3,3)  = 0.59000e-03_wp
            htt(13,4,3)  = -0.26000e-03_wp
            htt(13,5,3)  = 0.15200e-02_wp
            htt(13,6,3)  = 0.15100e-02_wp
            htt(13,7,3)  = -0.70000e-04_wp
            htt(13,8,3)  = -0.90000e-04_wp
            htt(13,9,3)  = -0.21100e-02_wp
            htt(13,10,3) = 0.12000e-03_wp
            htt(13,11,3) = -0.12300e-02_wp
            htt(13,12,3) = 0.10000e-02_wp
            htt(13,13,3) = 0.50300e-02_wp
            htt(14,1,3)  = 0.70000e-03_wp
            htt(14,2,3)  = 0.53000e-03_wp
            htt(14,3,3)  = -0.31000e-03_wp
            htt(14,4,3)  = -0.47000e-03_wp
            htt(14,5,3)  = -0.22000e-03_wp
            htt(14,6,3)  = -0.80000e-03_wp
            htt(14,7,3)  = 0.74000e-03_wp
            htt(14,8,3)  = 0.16000e-03_wp
            htt(14,9,3)  = -0.19000e-03_wp
            htt(14,10,3) = 0.13500e-02_wp
            htt(14,11,3) = -0.10000e-04_wp
            htt(14,12,3) = 0.15300e-02_wp
            htt(14,13,3) = -0.81000e-03_wp
            htt(14,14,3) = 0.39200e-02_wp
            htt(15,1,3)  = -0.19100e-02_wp
            htt(15,2,3)  = -0.28400e-02_wp
            htt(15,3,3)  = 0.70000e-03_wp
            htt(15,4,3)  = -0.23000e-03_wp
            htt(15,5,3)  = 0.15000e-02_wp
            htt(15,6,3)  = 0.21100e-02_wp
            htt(15,7,3)  = -0.98000e-03_wp
            htt(15,8,3)  = -0.37000e-03_wp
            htt(15,9,3)  = -0.24000e-03_wp
            htt(15,10,3) = -0.31800e-02_wp
            htt(15,11,3) = 0.73000e-03_wp
            htt(15,12,3) = -0.22200e-02_wp
            htt(15,13,3) = 0.28200e-02_wp
            htt(15,14,3) = -0.23300e-02_wp
            htt(15,15,3) = 0.57300e-02_wp
            htt(16,1,3)  = -0.12000e-03_wp
            htt(16,2,3)  = -0.37000e-03_wp
            htt(16,3,3)  = -0.19000e-03_wp
            htt(16,4,3)  = -0.45000e-03_wp
            htt(16,5,3)  = -0.25000e-03_wp
            htt(16,6,3)  = -0.19000e-03_wp
            htt(16,7,3)  = -0.17000e-03_wp
            htt(16,8,3)  = -0.12000e-03_wp
            htt(16,9,3)  = 0.16800e-02_wp
            htt(16,10,3) = -0.19500e-02_wp
            htt(16,11,3) = 0.19400e-02_wp
            htt(16,12,3) = -0.16900e-02_wp
            htt(16,13,3) = -0.30200e-02_wp
            htt(16,14,3) = 0.24000e-02_wp
            htt(16,15,3) = 0.58000e-03_wp
            htt(16,16,3) = 0.60000e-02_wp
            htt(17,1,3)  = 0.12100e-02_wp
            htt(17,2,3)  = 0.15200e-02_wp
            htt(17,3,3)  = -0.31000e-03_wp
            htt(17,4,3)  = 0.00000e+00_wp
            htt(17,5,3)  = -0.86000e-03_wp
            htt(17,6,3)  = -0.18800e-02_wp
            htt(17,7,3)  = 0.58000e-03_wp
            htt(17,8,3)  = -0.44000e-03_wp
            htt(17,9,3)  = 0.14500e-02_wp
            htt(17,10,3) = 0.90000e-04_wp
            htt(17,11,3) = 0.11300e-02_wp
            htt(17,12,3) = -0.23000e-03_wp
            htt(17,13,3) = -0.46000e-02_wp
            htt(17,14,3) = 0.13000e-02_wp
            htt(17,15,3) = -0.32400e-02_wp
            htt(17,16,3) = 0.26500e-02_wp
            htt(17,17,3) = 0.49600e-02_wp
            htt(18,1,3)  = 0.14500e-02_wp
            htt(18,2,3)  = 0.22200e-02_wp
            htt(18,3,3)  = -0.10600e-02_wp
            htt(18,4,3)  = -0.29000e-03_wp
            htt(18,5,3)  = -0.17100e-02_wp
            htt(18,6,3)  = -0.27100e-02_wp
            htt(18,7,3)  = 0.67000e-03_wp
            htt(18,8,3)  = -0.32000e-03_wp
            htt(18,9,3)  = 0.10000e-02_wp
            htt(18,10,3) = 0.82000e-03_wp
            htt(18,11,3) = 0.19000e-03_wp
            htt(18,12,3) = 0.10000e-04_wp
            htt(18,13,3) = -0.32400e-02_wp
            htt(18,14,3) = 0.19000e-02_wp
            htt(18,15,3) = -0.30600e-02_wp
            htt(18,16,3) = 0.20800e-02_wp
            htt(18,17,3) = 0.29700e-02_wp
            htt(18,18,3) = 0.56300e-02_wp
            htt(19,1,3)  = -0.59000e-03_wp
            htt(19,2,3)  = -0.96000e-03_wp
            htt(19,3,3)  = 0.59000e-03_wp
            htt(19,4,3)  = 0.21000e-03_wp
            htt(19,5,3)  = 0.89000e-03_wp
            htt(19,6,3)  = 0.43000e-03_wp
            htt(19,7,3)  = -0.23000e-03_wp
            htt(19,8,3)  = -0.69000e-03_wp
            htt(19,9,3)  = -0.48000e-03_wp
            htt(19,10,3) = -0.11500e-02_wp
            htt(19,11,3) = -0.90000e-04_wp
            htt(19,12,3) = -0.76000e-03_wp
            htt(19,13,3) = 0.12600e-02_wp
            htt(19,14,3) = -0.34500e-02_wp
            htt(19,15,3) = 0.19200e-02_wp
            htt(19,16,3) = -0.27800e-02_wp
            htt(19,17,3) = -0.95000e-03_wp
            htt(19,18,3) = -0.21800e-02_wp
            htt(19,19,3) = 0.37700e-02_wp
            htt(20,1,3)  = -0.35000e-03_wp
            htt(20,2,3)  = -0.26000e-03_wp
            htt(20,3,3)  = -0.16000e-03_wp
            htt(20,4,3)  = -0.80000e-04_wp
            htt(20,5,3)  = 0.40000e-04_wp
            htt(20,6,3)  = -0.40000e-03_wp
            htt(20,7,3)  = -0.14000e-03_wp
            htt(20,8,3)  = -0.57000e-03_wp
            htt(20,9,3)  = -0.93000e-03_wp
            htt(20,10,3) = -0.42000e-03_wp
            htt(20,11,3) = -0.10200e-02_wp
            htt(20,12,3) = -0.51000e-03_wp
            htt(20,13,3) = 0.26100e-02_wp
            htt(20,14,3) = -0.28400e-02_wp
            htt(20,15,3) = 0.21100e-02_wp
            htt(20,16,3) = -0.33500e-02_wp
            htt(20,17,3) = -0.29400e-02_wp
            htt(20,18,3) = 0.47000e-03_wp
            htt(20,19,3) = 0.25400e-02_wp
            htt(20,20,3) = 0.59400e-02_wp
            htt(21,1,3)  = 0.21700e-02_wp
            htt(21,2,3)  = 0.15000e-02_wp
            htt(21,3,3)  = -0.99000e-03_wp
            htt(21,4,3)  = -0.16600e-02_wp
            htt(21,5,3)  = -0.30200e-02_wp
            htt(21,6,3)  = -0.69000e-03_wp
            htt(21,7,3)  = -0.20000e-04_wp
            htt(21,8,3)  = 0.23100e-02_wp
            htt(21,9,3)  = -0.11900e-02_wp
            htt(21,10,3) = -0.98000e-03_wp
            htt(21,11,3) = -0.50000e-03_wp
            htt(21,12,3) = -0.28000e-03_wp
            htt(21,13,3) = 0.12600e-02_wp
            htt(21,14,3) = -0.76000e-03_wp
            htt(21,15,3) = 0.10500e-02_wp
            htt(21,16,3) = -0.97000e-03_wp
            htt(21,17,3) = -0.21400e-02_wp
            htt(21,18,3) = -0.48000e-03_wp
            htt(21,19,3) = -0.11000e-03_wp
            htt(21,20,3) = 0.15500e-02_wp
            htt(21,21,3) = 0.28900e-02_wp
            htt(22,1,3)  = 0.27000e-03_wp
            htt(22,2,3)  = 0.11300e-02_wp
            htt(22,3,3)  = -0.15900e-02_wp
            htt(22,4,3)  = -0.74000e-03_wp
            htt(22,5,3)  = -0.60000e-04_wp
            htt(22,6,3)  = -0.26200e-02_wp
            htt(22,7,3)  = 0.17100e-02_wp
            htt(22,8,3)  = -0.85000e-03_wp
            htt(22,9,3)  = -0.85000e-03_wp
            htt(22,10,3) = 0.40000e-04_wp
            htt(22,11,3) = -0.17300e-02_wp
            htt(22,12,3) = -0.85000e-03_wp
            htt(22,13,3) = 0.12800e-02_wp
            htt(22,14,3) = -0.15000e-03_wp
            htt(22,15,3) = 0.40000e-03_wp
            htt(22,16,3) = -0.10300e-02_wp
            htt(22,17,3) = -0.10700e-02_wp
            htt(22,18,3) = 0.56000e-03_wp
            htt(22,19,3) = 0.37000e-03_wp
            htt(22,20,3) = 0.20000e-02_wp
            htt(22,21,3) = 0.44000e-03_wp
            htt(22,22,3) = 0.31200e-02_wp
            htt(23,1,3)  = 0.19200e-02_wp
            htt(23,2,3)  = 0.79000e-03_wp
            htt(23,3,3)  = -0.23000e-03_wp
            htt(23,4,3)  = -0.13700e-02_wp
            htt(23,5,3)  = -0.21700e-02_wp
            htt(23,6,3)  = 0.14000e-03_wp
            htt(23,7,3)  = -0.12000e-03_wp
            htt(23,8,3)  = 0.21900e-02_wp
            htt(23,9,3)  = -0.74000e-03_wp
            htt(23,10,3) = -0.17100e-02_wp
            htt(23,11,3) = 0.45000e-03_wp
            htt(23,12,3) = -0.53000e-03_wp
            htt(23,13,3) = -0.10000e-03_wp
            htt(23,14,3) = -0.13700e-02_wp
            htt(23,15,3) = 0.87000e-03_wp
            htt(23,16,3) = -0.40000e-03_wp
            htt(23,17,3) = -0.15000e-03_wp
            htt(23,18,3) = -0.31500e-02_wp
            htt(23,19,3) = 0.11200e-02_wp
            htt(23,20,3) = -0.18800e-02_wp
            htt(23,21,3) = 0.12200e-02_wp
            htt(23,22,3) = -0.12000e-02_wp
            htt(23,23,3) = 0.42400e-02_wp
            htt(24,1,3)  = 0.30000e-04_wp
            htt(24,2,3)  = 0.42000e-03_wp
            htt(24,3,3)  = -0.84000e-03_wp
            htt(24,4,3)  = -0.45000e-03_wp
            htt(24,5,3)  = 0.79000e-03_wp
            htt(24,6,3)  = -0.18000e-02_wp
            htt(24,7,3)  = 0.16200e-02_wp
            htt(24,8,3)  = -0.97000e-03_wp
            htt(24,9,3)  = -0.39000e-03_wp
            htt(24,10,3) = -0.70000e-03_wp
            htt(24,11,3) = -0.79000e-03_wp
            htt(24,12,3) = -0.11000e-02_wp
            htt(24,13,3) = -0.80000e-04_wp
            htt(24,14,3) = -0.76000e-03_wp
            htt(24,15,3) = 0.22000e-03_wp
            htt(24,16,3) = -0.45000e-03_wp
            htt(24,17,3) = 0.92000e-03_wp
            htt(24,18,3) = -0.21000e-02_wp
            htt(24,19,3) = 0.16000e-02_wp
            htt(24,20,3) = -0.14300e-02_wp
            htt(24,21,3) = -0.12300e-02_wp
            htt(24,22,3) = 0.14800e-02_wp
            htt(24,23,3) = 0.18100e-02_wp
            htt(24,24,3) = 0.45300e-02_wp
            ! Hessian elements for anchor point 3
            hbs( 1,1,3)  = -0.31760e-01_wp
            hbs( 1,2,3)  = -0.29360e-01_wp
            hbs( 1,3,3)  = -0.31820e-01_wp
            hbs( 1,4,3)  = -0.95900e-02_wp
            hbs( 1,5,3)  = 0.52300e-02_wp
            hbs( 1,6,3)  = 0.39930e-01_wp
            hbs( 1,7,3)  = -0.45000e-03_wp
            hbs( 1,8,3)  = 0.40450e-01_wp
            hbs( 1,9,3)  = -0.18000e-03_wp
            hbs( 1,10,3) = -0.96900e-02_wp
            hbs( 1,11,3) = -0.42000e-03_wp
            hbs( 1,12,3) = 0.47400e-02_wp
            hbs( 2,1,3)  = 0.49930e-01_wp
            hbs( 2,2,3)  = -0.13930e-01_wp
            hbs( 2,3,3)  = 0.99600e-02_wp
            hbs( 2,4,3)  = 0.18000e-02_wp
            hbs( 2,5,3)  = -0.12200e-02_wp
            hbs( 2,6,3)  = -0.27810e-01_wp
            hbs( 2,7,3)  = 0.11800e-02_wp
            hbs( 2,8,3)  = -0.14650e-01_wp
            hbs( 2,9,3)  = 0.00000e+00_wp
            hbs( 2,10,3) = 0.41300e-02_wp
            hbs( 2,11,3) = -0.69000e-03_wp
            hbs( 2,12,3) = -0.89000e-03_wp
            hbs( 3,1,3)  = -0.18180e-01_wp
            hbs( 3,2,3)  = 0.43280e-01_wp
            hbs( 3,3,3)  = 0.21860e-01_wp
            hbs( 3,4,3)  = 0.77900e-02_wp
            hbs( 3,5,3)  = -0.40100e-02_wp
            hbs( 3,6,3)  = -0.12120e-01_wp
            hbs( 3,7,3)  = -0.73000e-03_wp
            hbs( 3,8,3)  = -0.25800e-01_wp
            hbs( 3,9,3)  = 0.18000e-03_wp
            hbs( 3,10,3) = 0.55600e-02_wp
            hbs( 3,11,3) = 0.11100e-02_wp
            hbs( 3,12,3) = -0.38500e-02_wp
            hbs( 4,1,3)  = -0.20820e-01_wp
            hbs( 4,2,3)  = -0.12100e-01_wp
            hbs( 4,3,3)  = 0.15180e-01_wp
            hbs( 4,4,3)  = -0.29840e-01_wp
            hbs( 4,5,3)  = -0.97300e-02_wp
            hbs( 4,6,3)  = -0.57000e-02_wp
            hbs( 4,7,3)  = 0.49500e-02_wp
            hbs( 4,8,3)  = 0.31150e-01_wp
            hbs( 4,9,3)  = -0.44000e-03_wp
            hbs( 4,10,3) = 0.40270e-01_wp
            hbs( 4,11,3) = -0.59000e-03_wp
            hbs( 4,12,3) = 0.29000e-03_wp
            hbs( 5,1,3)  = 0.26090e-01_wp
            hbs( 5,2,3)  = 0.58900e-02_wp
            hbs( 5,3,3)  = -0.14520e-01_wp
            hbs( 5,4,3)  = 0.25400e-02_wp
            hbs( 5,5,3)  = 0.69400e-02_wp
            hbs( 5,6,3)  = 0.30900e-02_wp
            hbs( 5,7,3)  = -0.21700e-02_wp
            hbs( 5,8,3)  = -0.15190e-01_wp
            hbs( 5,9,3)  = -0.39000e-03_wp
            hbs( 5,10,3) = -0.22140e-01_wp
            hbs( 5,11,3) = 0.28000e-03_wp
            hbs( 5,12,3) = 0.90000e-04_wp
            hbs( 6,1,3)  = -0.52700e-02_wp
            hbs( 6,2,3)  = 0.62100e-02_wp
            hbs( 6,3,3)  = -0.66000e-03_wp
            hbs( 6,4,3)  = 0.27300e-01_wp
            hbs( 6,5,3)  = 0.27900e-02_wp
            hbs( 6,6,3)  = 0.26100e-02_wp
            hbs( 6,7,3)  = -0.27800e-02_wp
            hbs( 6,8,3)  = -0.15960e-01_wp
            hbs( 6,9,3)  = 0.83000e-03_wp
            hbs( 6,10,3) = -0.18130e-01_wp
            hbs( 6,11,3) = 0.31000e-03_wp
            hbs( 6,12,3) = -0.38000e-03_wp
            hbs( 7,1,3)  = -0.90800e-02_wp
            hbs( 7,2,3)  = 0.35870e-01_wp
            hbs( 7,3,3)  = -0.10000e-04_wp
            hbs( 7,4,3)  = -0.22920e-01_wp
            hbs( 7,5,3)  = 0.61900e-02_wp
            hbs( 7,6,3)  = -0.28970e-01_wp
            hbs( 7,7,3)  = -0.87500e-02_wp
            hbs( 7,8,3)  = -0.92900e-02_wp
            hbs( 7,9,3)  = 0.49700e-02_wp
            hbs( 7,10,3) = 0.32780e-01_wp
            hbs( 7,11,3) = 0.60000e-04_wp
            hbs( 7,12,3) = -0.10700e-02_wp
            hbs( 8,1,3)  = 0.36900e-02_wp
            hbs( 8,2,3)  = -0.18590e-01_wp
            hbs( 8,3,3)  = 0.13800e-02_wp
            hbs( 8,4,3)  = 0.25350e-01_wp
            hbs( 8,5,3)  = -0.35000e-02_wp
            hbs( 8,6,3)  = -0.42000e-03_wp
            hbs( 8,7,3)  = 0.45000e-02_wp
            hbs( 8,8,3)  = 0.50400e-02_wp
            hbs( 8,9,3)  = -0.20100e-02_wp
            hbs( 8,10,3) = -0.14680e-01_wp
            hbs( 8,11,3) = -0.62000e-03_wp
            hbs( 8,12,3) = 0.55000e-03_wp
            hbs( 9,1,3)  = 0.53900e-02_wp
            hbs( 9,2,3)  = -0.17270e-01_wp
            hbs( 9,3,3)  = -0.13700e-02_wp
            hbs( 9,4,3)  = -0.24400e-02_wp
            hbs( 9,5,3)  = -0.26800e-02_wp
            hbs( 9,6,3)  = 0.29390e-01_wp
            hbs( 9,7,3)  = 0.42400e-02_wp
            hbs( 9,8,3)  = 0.42500e-02_wp
            hbs( 9,9,3)  = -0.29600e-02_wp
            hbs( 9,10,3) = -0.18100e-01_wp
            hbs( 9,11,3) = 0.55000e-03_wp
            hbs( 9,12,3) = 0.51000e-03_wp
            hbs(10,1,3)  = 0.36130e-01_wp
            hbs(10,2,3)  = 0.38360e-01_wp
            hbs(10,3,3)  = 0.36300e-02_wp
            hbs(10,4,3)  = -0.11890e-01_wp
            hbs(10,5,3)  = -0.10500e-02_wp
            hbs(10,6,3)  = -0.26840e-01_wp
            hbs(10,7,3)  = 0.48100e-02_wp
            hbs(10,8,3)  = -0.27100e-01_wp
            hbs(10,9,3)  = -0.88200e-02_wp
            hbs(10,10,3) = -0.12140e-01_wp
            hbs(10,11,3) = 0.47600e-02_wp
            hbs(10,12,3) = -0.77000e-03_wp
            hbs(11,1,3)  = -0.20190e-01_wp
            hbs(11,2,3)  = -0.17180e-01_wp
            hbs(11,3,3)  = -0.21300e-02_wp
            hbs(11,4,3)  = 0.57400e-02_wp
            hbs(11,5,3)  = 0.12900e-02_wp
            hbs(11,6,3)  = 0.27000e-01_wp
            hbs(11,7,3)  = -0.28400e-02_wp
            hbs(11,8,3)  = 0.47000e-03_wp
            hbs(11,9,3)  = 0.44300e-02_wp
            hbs(11,10,3) = 0.60700e-02_wp
            hbs(11,11,3) = -0.19800e-02_wp
            hbs(11,12,3) = -0.18000e-03_wp
            hbs(12,1,3)  = -0.15940e-01_wp
            hbs(12,2,3)  = -0.21180e-01_wp
            hbs(12,3,3)  = -0.14900e-02_wp
            hbs(12,4,3)  = 0.61500e-02_wp
            hbs(12,5,3)  = -0.24000e-03_wp
            hbs(12,6,3)  = -0.16000e-03_wp
            hbs(12,7,3)  = -0.19600e-02_wp
            hbs(12,8,3)  = 0.26630e-01_wp
            hbs(12,9,3)  = 0.43900e-02_wp
            hbs(12,10,3) = 0.60700e-02_wp
            hbs(12,11,3) = -0.27800e-02_wp
            hbs(12,12,3) = 0.94000e-03_wp
            hbs(13,1,3)  = 0.35720e-01_wp
            hbs(13,2,3)  = -0.96900e-02_wp
            hbs(13,3,3)  = 0.29800e-02_wp
            hbs(13,4,3)  = 0.32310e-01_wp
            hbs(13,5,3)  = -0.94000e-03_wp
            hbs(13,6,3)  = -0.10140e-01_wp
            hbs(13,7,3)  = 0.70000e-04_wp
            hbs(13,8,3)  = -0.29530e-01_wp
            hbs(13,9,3)  = 0.50200e-02_wp
            hbs(13,10,3) = -0.23420e-01_wp
            hbs(13,11,3) = -0.86600e-02_wp
            hbs(13,12,3) = 0.53400e-02_wp
            hbs(14,1,3)  = -0.17210e-01_wp
            hbs(14,2,3)  = 0.56100e-02_wp
            hbs(14,3,3)  = -0.23300e-02_wp
            hbs(14,4,3)  = -0.17680e-01_wp
            hbs(14,5,3)  = 0.48000e-03_wp
            hbs(14,6,3)  = 0.47700e-02_wp
            hbs(14,7,3)  = 0.56000e-03_wp
            hbs(14,8,3)  = 0.29630e-01_wp
            hbs(14,9,3)  = -0.29700e-02_wp
            hbs(14,10,3) = -0.24000e-02_wp
            hbs(14,11,3) = 0.43000e-02_wp
            hbs(14,12,3) = -0.22600e-02_wp
            hbs(15,1,3)  = -0.18510e-01_wp
            hbs(15,2,3)  = 0.40800e-02_wp
            hbs(15,3,3)  = -0.65000e-03_wp
            hbs(15,4,3)  = -0.14630e-01_wp
            hbs(15,5,3)  = 0.46000e-03_wp
            hbs(15,6,3)  = 0.53700e-02_wp
            hbs(15,7,3)  = -0.63000e-03_wp
            hbs(15,8,3)  = -0.10000e-03_wp
            hbs(15,9,3)  = -0.20500e-02_wp
            hbs(15,10,3) = 0.25820e-01_wp
            hbs(15,11,3) = 0.43700e-02_wp
            hbs(15,12,3) = -0.30900e-02_wp
            hbs(16,1,3)  = -0.10200e-01_wp
            hbs(16,2,3)  = -0.23080e-01_wp
            hbs(16,3,3)  = 0.10040e-01_wp
            hbs(16,4,3)  = 0.41930e-01_wp
            hbs(16,5,3)  = 0.31000e-03_wp
            hbs(16,6,3)  = 0.31710e-01_wp
            hbs(16,7,3)  = -0.63000e-03_wp
            hbs(16,8,3)  = -0.56800e-02_wp
            hbs(16,9,3)  = -0.55000e-03_wp
            hbs(16,10,3) = -0.27810e-01_wp
            hbs(16,11,3) = 0.48500e-02_wp
            hbs(16,12,3) = -0.85400e-02_wp
            hbs(17,1,3)  = 0.41100e-02_wp
            hbs(17,2,3)  = 0.25000e-01_wp
            hbs(17,3,3)  = -0.74400e-02_wp
            hbs(17,4,3)  = -0.22980e-01_wp
            hbs(17,5,3)  = 0.60000e-03_wp
            hbs(17,6,3)  = -0.14860e-01_wp
            hbs(17,7,3)  = 0.32000e-03_wp
            hbs(17,8,3)  = 0.33600e-02_wp
            hbs(17,9,3)  = -0.32000e-03_wp
            hbs(17,10,3) = 0.20400e-02_wp
            hbs(17,11,3) = -0.21000e-02_wp
            hbs(17,12,3) = 0.68800e-02_wp
            hbs(18,1,3)  = 0.61000e-02_wp
            hbs(18,2,3)  = -0.19200e-02_wp
            hbs(18,3,3)  = -0.26000e-02_wp
            hbs(18,4,3)  = -0.18950e-01_wp
            hbs(18,5,3)  = -0.92000e-03_wp
            hbs(18,6,3)  = -0.16850e-01_wp
            hbs(18,7,3)  = 0.31000e-03_wp
            hbs(18,8,3)  = 0.23200e-02_wp
            hbs(18,9,3)  = 0.87000e-03_wp
            hbs(18,10,3) = 0.25770e-01_wp
            hbs(18,11,3) = -0.27500e-02_wp
            hbs(18,12,3) = 0.16600e-02_wp
            ! Hessian elements for anchor point 4
            hss( 1,1,4)  = 0.34759e+00_wp
            hss( 2,1,4)  = 0.17810e-01_wp
            hss( 2,2,4)  = 0.35623e+00_wp
            hss( 3,1,4)  = 0.29010e-01_wp
            hss( 3,2,4)  = 0.32710e-01_wp
            hss( 3,3,4)  = 0.46124e+00_wp
            hss( 4,1,4)  = -0.19600e-02_wp
            hss( 4,2,4)  = 0.14870e-01_wp
            hss( 4,3,4)  = -0.89100e-02_wp
            hss( 4,4,4)  = 0.38581e+00_wp
            hss( 5,1,4)  = 0.38800e-02_wp
            hss( 5,2,4)  = -0.77000e-03_wp
            hss( 5,3,4)  = 0.71700e-02_wp
            hss( 5,4,4)  = 0.36500e-02_wp
            hss( 5,5,4)  = 0.35807e+00_wp
            hss( 6,1,4)  = 0.33890e-01_wp
            hss( 6,2,4)  = 0.44570e-01_wp
            hss( 6,3,4)  = -0.43100e-02_wp
            hss( 6,4,4)  = -0.36100e-02_wp
            hss( 6,5,4)  = 0.55000e-03_wp
            hss( 6,6,4)  = 0.35931e+00_wp
            hss( 7,1,4)  = -0.54000e-03_wp
            hss( 7,2,4)  = -0.63000e-03_wp
            hss( 7,3,4)  = 0.57000e-03_wp
            hss( 7,4,4)  = 0.39700e-02_wp
            hss( 7,5,4)  = 0.39000e-03_wp
            hss( 7,6,4)  = 0.39800e-02_wp
            hss( 7,7,4)  = 0.37091e+00_wp
            hss( 8,1,4)  = 0.42030e-01_wp
            hss( 8,2,4)  = 0.35060e-01_wp
            hss( 8,3,4)  = -0.52800e-02_wp
            hss( 8,4,4)  = 0.24280e-01_wp
            hss( 8,5,4)  = -0.11100e-02_wp
            hss( 8,6,4)  = 0.17880e-01_wp
            hss( 8,7,4)  = -0.10400e-02_wp
            hss( 8,8,4)  = 0.35888e+00_wp
            hss( 9,1,4)  = -0.53000e-03_wp
            hss( 9,2,4)  = -0.57000e-03_wp
            hss( 9,3,4)  = 0.26000e-03_wp
            hss( 9,4,4)  = -0.39000e-03_wp
            hss( 9,5,4)  = 0.19000e-03_wp
            hss( 9,6,4)  = 0.40100e-02_wp
            hss( 9,7,4)  = 0.49000e-03_wp
            hss( 9,8,4)  = 0.39600e-02_wp
            hss( 9,9,4)  = 0.36356e+00_wp
            hss(10,1,4)  = 0.19520e-01_wp
            hss(10,2,4)  = -0.69300e-02_wp
            hss(10,3,4)  = -0.48800e-02_wp
            hss(10,4,4)  = 0.64010e-01_wp
            hss(10,5,4)  = -0.43000e-03_wp
            hss(10,6,4)  = 0.21230e-01_wp
            hss(10,7,4)  = -0.10900e-02_wp
            hss(10,8,4)  = -0.77000e-03_wp
            hss(10,9,4)  = -0.19000e-03_wp
            hss(10,10,4) = 0.38116e+00_wp
            hss(11,1,4)  = -0.69000e-03_wp
            hss(11,2,4)  = -0.50000e-03_wp
            hss(11,3,4)  = 0.52000e-03_wp
            hss(11,4,4)  = -0.11600e-02_wp
            hss(11,5,4)  = 0.80000e-04_wp
            hss(11,6,4)  = -0.10000e-02_wp
            hss(11,7,4)  = 0.22000e-03_wp
            hss(11,8,4)  = 0.39700e-02_wp
            hss(11,9,4)  = 0.48000e-03_wp
            hss(11,10,4) = 0.39800e-02_wp
            hss(11,11,4) = 0.37152e+00_wp
            hss(12,1,4)  = -0.16500e-02_wp
            hss(12,2,4)  = 0.39400e-02_wp
            hss(12,3,4)  = -0.19000e-03_wp
            hss(12,4,4)  = -0.10000e-03_wp
            hss(12,5,4)  = 0.35000e-03_wp
            hss(12,6,4)  = -0.10400e-02_wp
            hss(12,7,4)  = 0.30000e-04_wp
            hss(12,8,4)  = 0.57000e-03_wp
            hss(12,9,4)  = 0.12000e-03_wp
            hss(12,10,4) = 0.23600e-02_wp
            hss(12,11,4) = 0.40000e-03_wp
            hss(12,12,4) = 0.37254e+00_wp
            ! Hessian elements for anchor point 4
            hbb( 1,1,4)  = 0.82450e-01_wp
            hbb( 2,1,4)  = -0.46090e-01_wp
            hbb( 2,2,4)  = 0.14487e+00_wp
            hbb( 3,1,4)  = -0.36360e-01_wp
            hbb( 3,2,4)  = -0.98780e-01_wp
            hbb( 3,3,4)  = 0.13515e+00_wp
            hbb( 4,1,4)  = -0.36000e-01_wp
            hbb( 4,2,4)  = 0.23310e-01_wp
            hbb( 4,3,4)  = 0.12690e-01_wp
            hbb( 4,4,4)  = 0.66170e-01_wp
            hbb( 5,1,4)  = 0.20240e-01_wp
            hbb( 5,2,4)  = -0.12480e-01_wp
            hbb( 5,3,4)  = -0.77700e-02_wp
            hbb( 5,4,4)  = -0.33030e-01_wp
            hbb( 5,5,4)  = 0.71950e-01_wp
            hbb( 6,1,4)  = 0.15760e-01_wp
            hbb( 6,2,4)  = -0.10840e-01_wp
            hbb( 6,3,4)  = -0.49200e-02_wp
            hbb( 6,4,4)  = -0.33140e-01_wp
            hbb( 6,5,4)  = -0.38920e-01_wp
            hbb( 6,6,4)  = 0.72060e-01_wp
            hbb( 7,1,4)  = -0.71500e-02_wp
            hbb( 7,2,4)  = -0.39000e-02_wp
            hbb( 7,3,4)  = 0.11050e-01_wp
            hbb( 7,4,4)  = -0.25740e-01_wp
            hbb( 7,5,4)  = 0.93100e-02_wp
            hbb( 7,6,4)  = 0.16430e-01_wp
            hbb( 7,7,4)  = 0.64420e-01_wp
            hbb( 8,1,4)  = -0.22800e-02_wp
            hbb( 8,2,4)  = 0.74300e-02_wp
            hbb( 8,3,4)  = -0.51500e-02_wp
            hbb( 8,4,4)  = 0.15690e-01_wp
            hbb( 8,5,4)  = -0.49500e-02_wp
            hbb( 8,6,4)  = -0.10740e-01_wp
            hbb( 8,7,4)  = -0.31490e-01_wp
            hbb( 8,8,4)  = 0.72950e-01_wp
            hbb( 9,1,4)  = 0.94300e-02_wp
            hbb( 9,2,4)  = -0.35300e-02_wp
            hbb( 9,3,4)  = -0.59100e-02_wp
            hbb( 9,4,4)  = 0.10050e-01_wp
            hbb( 9,5,4)  = -0.43500e-02_wp
            hbb( 9,6,4)  = -0.56900e-02_wp
            hbb( 9,7,4)  = -0.32930e-01_wp
            hbb( 9,8,4)  = -0.41460e-01_wp
            hbb( 9,9,4)  = 0.74390e-01_wp
            hbb(10,1,4)  = 0.45500e-02_wp
            hbb(10,2,4)  = -0.20200e-02_wp
            hbb(10,3,4)  = -0.25300e-02_wp
            hbb(10,4,4)  = -0.64100e-02_wp
            hbb(10,5,4)  = 0.91900e-02_wp
            hbb(10,6,4)  = -0.27700e-02_wp
            hbb(10,7,4)  = -0.33410e-01_wp
            hbb(10,8,4)  = 0.13410e-01_wp
            hbb(10,9,4)  = 0.20000e-01_wp
            hbb(10,10,4) = 0.75440e-01_wp
            hbb(11,1,4)  = -0.22100e-02_wp
            hbb(11,2,4)  = 0.12300e-02_wp
            hbb(11,3,4)  = 0.98000e-03_wp
            hbb(11,4,4)  = -0.23100e-02_wp
            hbb(11,5,4)  = -0.32200e-02_wp
            hbb(11,6,4)  = 0.55300e-02_wp
            hbb(11,7,4)  = 0.20020e-01_wp
            hbb(11,8,4)  = -0.74100e-02_wp
            hbb(11,9,4)  = -0.12610e-01_wp
            hbb(11,10,4) = -0.37790e-01_wp
            hbb(11,11,4) = 0.78360e-01_wp
            hbb(12,1,4)  = -0.23400e-02_wp
            hbb(12,2,4)  = 0.79000e-03_wp
            hbb(12,3,4)  = 0.15500e-02_wp
            hbb(12,4,4)  = 0.87200e-02_wp
            hbb(12,5,4)  = -0.59700e-02_wp
            hbb(12,6,4)  = -0.27500e-02_wp
            hbb(12,7,4)  = 0.13390e-01_wp
            hbb(12,8,4)  = -0.60000e-02_wp
            hbb(12,9,4)  = -0.73900e-02_wp
            hbb(12,10,4) = -0.37640e-01_wp
            hbb(12,11,4) = -0.40570e-01_wp
            hbb(12,12,4) = 0.78220e-01_wp
            hbb(13,1,4)  = -0.67600e-02_wp
            hbb(13,2,4)  = 0.10710e-01_wp
            hbb(13,3,4)  = -0.39500e-02_wp
            hbb(13,4,4)  = -0.11500e-02_wp
            hbb(13,5,4)  = 0.12900e-02_wp
            hbb(13,6,4)  = -0.14000e-03_wp
            hbb(13,7,4)  = 0.29100e-02_wp
            hbb(13,8,4)  = 0.41300e-02_wp
            hbb(13,9,4)  = -0.70300e-02_wp
            hbb(13,10,4) = -0.33760e-01_wp
            hbb(13,11,4) = 0.13510e-01_wp
            hbb(13,12,4) = 0.20240e-01_wp
            hbb(13,13,4) = 0.65300e-01_wp
            hbb(14,1,4)  = 0.91200e-02_wp
            hbb(14,2,4)  = -0.60500e-02_wp
            hbb(14,3,4)  = -0.30700e-02_wp
            hbb(14,4,4)  = 0.66000e-03_wp
            hbb(14,5,4)  = -0.85000e-03_wp
            hbb(14,6,4)  = 0.19000e-03_wp
            hbb(14,7,4)  = -0.71100e-02_wp
            hbb(14,8,4)  = -0.51000e-03_wp
            hbb(14,9,4)  = 0.76200e-02_wp
            hbb(14,10,4) = 0.20030e-01_wp
            hbb(14,11,4) = -0.73800e-02_wp
            hbb(14,12,4) = -0.12650e-01_wp
            hbb(14,13,4) = -0.33180e-01_wp
            hbb(14,14,4) = 0.74350e-01_wp
            hbb(15,1,4)  = -0.23600e-02_wp
            hbb(15,2,4)  = -0.46500e-02_wp
            hbb(15,3,4)  = 0.70200e-02_wp
            hbb(15,4,4)  = 0.49000e-03_wp
            hbb(15,5,4)  = -0.44000e-03_wp
            hbb(15,6,4)  = -0.50000e-04_wp
            hbb(15,7,4)  = 0.42000e-02_wp
            hbb(15,8,4)  = -0.36200e-02_wp
            hbb(15,9,4)  = -0.59000e-03_wp
            hbb(15,10,4) = 0.13730e-01_wp
            hbb(15,11,4) = -0.61300e-02_wp
            hbb(15,12,4) = -0.76000e-02_wp
            hbb(15,13,4) = -0.32120e-01_wp
            hbb(15,14,4) = -0.41170e-01_wp
            hbb(15,15,4) = 0.73290e-01_wp
            hbb(16,1,4)  = -0.37090e-01_wp
            hbb(16,2,4)  = 0.17990e-01_wp
            hbb(16,3,4)  = 0.19100e-01_wp
            hbb(16,4,4)  = 0.31400e-02_wp
            hbb(16,5,4)  = -0.70000e-02_wp
            hbb(16,6,4)  = 0.38600e-02_wp
            hbb(16,7,4)  = -0.10300e-02_wp
            hbb(16,8,4)  = 0.54000e-03_wp
            hbb(16,9,4)  = 0.48000e-03_wp
            hbb(16,10,4) = -0.64000e-02_wp
            hbb(16,11,4) = 0.87800e-02_wp
            hbb(16,12,4) = -0.23700e-02_wp
            hbb(16,13,4) = -0.26530e-01_wp
            hbb(16,14,4) = 0.10480e-01_wp
            hbb(16,15,4) = 0.16060e-01_wp
            hbb(16,16,4) = 0.67910e-01_wp
            hbb(17,1,4)  = 0.22610e-01_wp
            hbb(17,2,4)  = -0.90200e-02_wp
            hbb(17,3,4)  = -0.13590e-01_wp
            hbb(17,4,4)  = -0.76000e-02_wp
            hbb(17,5,4)  = 0.77900e-02_wp
            hbb(17,6,4)  = -0.19000e-03_wp
            hbb(17,7,4)  = 0.87000e-03_wp
            hbb(17,8,4)  = -0.25000e-03_wp
            hbb(17,9,4)  = -0.62000e-03_wp
            hbb(17,10,4) = 0.91600e-02_wp
            hbb(17,11,4) = -0.60700e-02_wp
            hbb(17,12,4) = -0.31000e-02_wp
            hbb(17,13,4) = 0.96800e-02_wp
            hbb(17,14,4) = -0.46000e-02_wp
            hbb(17,15,4) = -0.50800e-02_wp
            hbb(17,16,4) = -0.34720e-01_wp
            hbb(17,17,4) = 0.72390e-01_wp
            hbb(18,1,4)  = 0.14470e-01_wp
            hbb(18,2,4)  = -0.89700e-02_wp
            hbb(18,3,4)  = -0.55000e-02_wp
            hbb(18,4,4)  = 0.44600e-02_wp
            hbb(18,5,4)  = -0.79000e-03_wp
            hbb(18,6,4)  = -0.36700e-02_wp
            hbb(18,7,4)  = 0.16000e-03_wp
            hbb(18,8,4)  = -0.29000e-03_wp
            hbb(18,9,4)  = 0.13000e-03_wp
            hbb(18,10,4) = -0.27600e-02_wp
            hbb(18,11,4) = -0.27100e-02_wp
            hbb(18,12,4) = 0.54700e-02_wp
            hbb(18,13,4) = 0.16860e-01_wp
            hbb(18,14,4) = -0.58800e-02_wp
            hbb(18,15,4) = -0.10980e-01_wp
            hbb(18,16,4) = -0.33190e-01_wp
            hbb(18,17,4) = -0.37670e-01_wp
            hbb(18,18,4) = 0.70860e-01_wp
            ! Hessian elements for anchor point 4
            htt( 1,1,4)  = 0.60400e-02_wp
            htt( 2,1,4)  = 0.38400e-02_wp
            htt( 2,2,4)  = 0.50900e-02_wp
            htt( 3,1,4)  = -0.18500e-02_wp
            htt( 3,2,4)  = -0.33100e-02_wp
            htt( 3,3,4)  = 0.75100e-02_wp
            htt( 4,1,4)  = -0.40500e-02_wp
            htt( 4,2,4)  = -0.20700e-02_wp
            htt( 4,3,4)  = 0.60500e-02_wp
            htt( 4,4,4)  = 0.80200e-02_wp
            htt( 5,1,4)  = -0.59200e-02_wp
            htt( 5,2,4)  = -0.42600e-02_wp
            htt( 5,3,4)  = 0.21200e-02_wp
            htt( 5,4,4)  = 0.37800e-02_wp
            htt( 5,5,4)  = 0.65600e-02_wp
            htt( 6,1,4)  = -0.41200e-02_wp
            htt( 6,2,4)  = -0.39000e-02_wp
            htt( 6,3,4)  = 0.26900e-02_wp
            htt( 6,4,4)  = 0.29100e-02_wp
            htt( 6,5,4)  = 0.37400e-02_wp
            htt( 6,6,4)  = 0.55900e-02_wp
            htt( 7,1,4)  = 0.15700e-02_wp
            htt( 7,2,4)  = 0.25400e-02_wp
            htt( 7,3,4)  = -0.67700e-02_wp
            htt( 7,4,4)  = -0.58000e-02_wp
            htt( 7,5,4)  = -0.10700e-02_wp
            htt( 7,6,4)  = -0.27200e-02_wp
            htt( 7,7,4)  = 0.68600e-02_wp
            htt( 8,1,4)  = 0.33800e-02_wp
            htt( 8,2,4)  = 0.29000e-02_wp
            htt( 8,3,4)  = -0.62000e-02_wp
            htt( 8,4,4)  = -0.66800e-02_wp
            htt( 8,5,4)  = -0.38900e-02_wp
            htt( 8,6,4)  = -0.88000e-03_wp
            htt( 8,7,4)  = 0.52100e-02_wp
            htt( 8,8,4)  = 0.82200e-02_wp
            htt( 9,1,4)  = -0.23900e-02_wp
            htt( 9,2,4)  = -0.68000e-03_wp
            htt( 9,3,4)  = 0.45000e-03_wp
            htt( 9,4,4)  = 0.21600e-02_wp
            htt( 9,5,4)  = 0.17400e-02_wp
            htt( 9,6,4)  = 0.14100e-02_wp
            htt( 9,7,4)  = -0.96000e-03_wp
            htt( 9,8,4)  = -0.12900e-02_wp
            htt( 9,9,4)  = 0.24600e-02_wp
            htt(10,1,4)  = -0.15600e-02_wp
            htt(10,2,4)  = 0.23000e-03_wp
            htt(10,3,4)  = 0.34000e-03_wp
            htt(10,4,4)  = 0.21300e-02_wp
            htt(10,5,4)  = 0.17600e-02_wp
            htt(10,6,4)  = 0.79000e-03_wp
            htt(10,7,4)  = -0.40000e-04_wp
            htt(10,8,4)  = -0.10100e-02_wp
            htt(10,9,4)  = 0.57000e-03_wp
            htt(10,10,4) = 0.39100e-02_wp
            htt(11,1,4)  = -0.10000e-03_wp
            htt(11,2,4)  = -0.19700e-02_wp
            htt(11,3,4)  = 0.19800e-02_wp
            htt(11,4,4)  = 0.11000e-03_wp
            htt(11,5,4)  = 0.10000e-04_wp
            htt(11,6,4)  = 0.11900e-02_wp
            htt(11,7,4)  = -0.19700e-02_wp
            htt(11,8,4)  = -0.79000e-03_wp
            htt(11,9,4)  = 0.68000e-03_wp
            htt(11,10,4) = -0.12900e-02_wp
            htt(11,11,4) = 0.26200e-02_wp
            htt(12,1,4)  = 0.72000e-03_wp
            htt(12,2,4)  = -0.10500e-02_wp
            htt(12,3,4)  = 0.18600e-02_wp
            htt(12,4,4)  = 0.90000e-04_wp
            htt(12,5,4)  = 0.30000e-04_wp
            htt(12,6,4)  = 0.57000e-03_wp
            htt(12,7,4)  = -0.10500e-02_wp
            htt(12,8,4)  = -0.51000e-03_wp
            htt(12,9,4)  = -0.12000e-02_wp
            htt(12,10,4) = 0.20500e-02_wp
            htt(12,11,4) = 0.64000e-03_wp
            htt(12,12,4) = 0.39000e-02_wp
            htt(13,1,4)  = -0.10900e-02_wp
            htt(13,2,4)  = -0.19400e-02_wp
            htt(13,3,4)  = 0.59000e-03_wp
            htt(13,4,4)  = -0.26000e-03_wp
            htt(13,5,4)  = 0.15200e-02_wp
            htt(13,6,4)  = 0.15100e-02_wp
            htt(13,7,4)  = -0.70000e-04_wp
            htt(13,8,4)  = -0.90000e-04_wp
            htt(13,9,4)  = -0.21100e-02_wp
            htt(13,10,4) = 0.12000e-03_wp
            htt(13,11,4) = -0.12300e-02_wp
            htt(13,12,4) = 0.10000e-02_wp
            htt(13,13,4) = 0.50300e-02_wp
            htt(14,1,4)  = 0.70000e-03_wp
            htt(14,2,4)  = 0.53000e-03_wp
            htt(14,3,4)  = -0.31000e-03_wp
            htt(14,4,4)  = -0.47000e-03_wp
            htt(14,5,4)  = -0.22000e-03_wp
            htt(14,6,4)  = -0.80000e-03_wp
            htt(14,7,4)  = 0.74000e-03_wp
            htt(14,8,4)  = 0.16000e-03_wp
            htt(14,9,4)  = -0.19000e-03_wp
            htt(14,10,4) = 0.13500e-02_wp
            htt(14,11,4) = -0.10000e-04_wp
            htt(14,12,4) = 0.15300e-02_wp
            htt(14,13,4) = -0.81000e-03_wp
            htt(14,14,4) = 0.39200e-02_wp
            htt(15,1,4)  = -0.19100e-02_wp
            htt(15,2,4)  = -0.28400e-02_wp
            htt(15,3,4)  = 0.70000e-03_wp
            htt(15,4,4)  = -0.23000e-03_wp
            htt(15,5,4)  = 0.15000e-02_wp
            htt(15,6,4)  = 0.21100e-02_wp
            htt(15,7,4)  = -0.98000e-03_wp
            htt(15,8,4)  = -0.37000e-03_wp
            htt(15,9,4)  = -0.24000e-03_wp
            htt(15,10,4) = -0.31800e-02_wp
            htt(15,11,4) = 0.73000e-03_wp
            htt(15,12,4) = -0.22200e-02_wp
            htt(15,13,4) = 0.28200e-02_wp
            htt(15,14,4) = -0.23300e-02_wp
            htt(15,15,4) = 0.57300e-02_wp
            htt(16,1,4)  = -0.12000e-03_wp
            htt(16,2,4)  = -0.37000e-03_wp
            htt(16,3,4)  = -0.19000e-03_wp
            htt(16,4,4)  = -0.45000e-03_wp
            htt(16,5,4)  = -0.25000e-03_wp
            htt(16,6,4)  = -0.19000e-03_wp
            htt(16,7,4)  = -0.17000e-03_wp
            htt(16,8,4)  = -0.12000e-03_wp
            htt(16,9,4)  = 0.16800e-02_wp
            htt(16,10,4) = -0.19500e-02_wp
            htt(16,11,4) = 0.19400e-02_wp
            htt(16,12,4) = -0.16900e-02_wp
            htt(16,13,4) = -0.30200e-02_wp
            htt(16,14,4) = 0.24000e-02_wp
            htt(16,15,4) = 0.58000e-03_wp
            htt(16,16,4) = 0.60000e-02_wp
            htt(17,1,4)  = 0.12100e-02_wp
            htt(17,2,4)  = 0.15200e-02_wp
            htt(17,3,4)  = -0.31000e-03_wp
            htt(17,4,4)  = 0.00000e+00_wp
            htt(17,5,4)  = -0.86000e-03_wp
            htt(17,6,4)  = -0.18800e-02_wp
            htt(17,7,4)  = 0.58000e-03_wp
            htt(17,8,4)  = -0.44000e-03_wp
            htt(17,9,4)  = 0.14500e-02_wp
            htt(17,10,4) = 0.90000e-04_wp
            htt(17,11,4) = 0.11300e-02_wp
            htt(17,12,4) = -0.23000e-03_wp
            htt(17,13,4) = -0.46000e-02_wp
            htt(17,14,4) = 0.13000e-02_wp
            htt(17,15,4) = -0.32400e-02_wp
            htt(17,16,4) = 0.26500e-02_wp
            htt(17,17,4) = 0.49600e-02_wp
            htt(18,1,4)  = 0.14500e-02_wp
            htt(18,2,4)  = 0.22200e-02_wp
            htt(18,3,4)  = -0.10600e-02_wp
            htt(18,4,4)  = -0.29000e-03_wp
            htt(18,5,4)  = -0.17100e-02_wp
            htt(18,6,4)  = -0.27100e-02_wp
            htt(18,7,4)  = 0.67000e-03_wp
            htt(18,8,4)  = -0.32000e-03_wp
            htt(18,9,4)  = 0.10000e-02_wp
            htt(18,10,4) = 0.82000e-03_wp
            htt(18,11,4) = 0.19000e-03_wp
            htt(18,12,4) = 0.10000e-04_wp
            htt(18,13,4) = -0.32400e-02_wp
            htt(18,14,4) = 0.19000e-02_wp
            htt(18,15,4) = -0.30600e-02_wp
            htt(18,16,4) = 0.20800e-02_wp
            htt(18,17,4) = 0.29700e-02_wp
            htt(18,18,4) = 0.56300e-02_wp
            htt(19,1,4)  = -0.59000e-03_wp
            htt(19,2,4)  = -0.96000e-03_wp
            htt(19,3,4)  = 0.59000e-03_wp
            htt(19,4,4)  = 0.21000e-03_wp
            htt(19,5,4)  = 0.89000e-03_wp
            htt(19,6,4)  = 0.43000e-03_wp
            htt(19,7,4)  = -0.23000e-03_wp
            htt(19,8,4)  = -0.69000e-03_wp
            htt(19,9,4)  = -0.48000e-03_wp
            htt(19,10,4) = -0.11500e-02_wp
            htt(19,11,4) = -0.90000e-04_wp
            htt(19,12,4) = -0.76000e-03_wp
            htt(19,13,4) = 0.12600e-02_wp
            htt(19,14,4) = -0.34500e-02_wp
            htt(19,15,4) = 0.19200e-02_wp
            htt(19,16,4) = -0.27800e-02_wp
            htt(19,17,4) = -0.95000e-03_wp
            htt(19,18,4) = -0.21800e-02_wp
            htt(19,19,4) = 0.37700e-02_wp
            htt(20,1,4)  = -0.35000e-03_wp
            htt(20,2,4)  = -0.26000e-03_wp
            htt(20,3,4)  = -0.16000e-03_wp
            htt(20,4,4)  = -0.80000e-04_wp
            htt(20,5,4)  = 0.40000e-04_wp
            htt(20,6,4)  = -0.40000e-03_wp
            htt(20,7,4)  = -0.14000e-03_wp
            htt(20,8,4)  = -0.57000e-03_wp
            htt(20,9,4)  = -0.93000e-03_wp
            htt(20,10,4) = -0.42000e-03_wp
            htt(20,11,4) = -0.10200e-02_wp
            htt(20,12,4) = -0.51000e-03_wp
            htt(20,13,4) = 0.26100e-02_wp
            htt(20,14,4) = -0.28400e-02_wp
            htt(20,15,4) = 0.21100e-02_wp
            htt(20,16,4) = -0.33500e-02_wp
            htt(20,17,4) = -0.29400e-02_wp
            htt(20,18,4) = 0.47000e-03_wp
            htt(20,19,4) = 0.25400e-02_wp
            htt(20,20,4) = 0.59400e-02_wp
            htt(21,1,4)  = 0.21700e-02_wp
            htt(21,2,4)  = 0.15000e-02_wp
            htt(21,3,4)  = -0.99000e-03_wp
            htt(21,4,4)  = -0.16600e-02_wp
            htt(21,5,4)  = -0.30200e-02_wp
            htt(21,6,4)  = -0.69000e-03_wp
            htt(21,7,4)  = -0.20000e-04_wp
            htt(21,8,4)  = 0.23100e-02_wp
            htt(21,9,4)  = -0.11900e-02_wp
            htt(21,10,4) = -0.98000e-03_wp
            htt(21,11,4) = -0.50000e-03_wp
            htt(21,12,4) = -0.28000e-03_wp
            htt(21,13,4) = 0.12600e-02_wp
            htt(21,14,4) = -0.76000e-03_wp
            htt(21,15,4) = 0.10500e-02_wp
            htt(21,16,4) = -0.97000e-03_wp
            htt(21,17,4) = -0.21400e-02_wp
            htt(21,18,4) = -0.48000e-03_wp
            htt(21,19,4) = -0.11000e-03_wp
            htt(21,20,4) = 0.15500e-02_wp
            htt(21,21,4) = 0.28900e-02_wp
            htt(22,1,4)  = 0.27000e-03_wp
            htt(22,2,4)  = 0.11300e-02_wp
            htt(22,3,4)  = -0.15900e-02_wp
            htt(22,4,4)  = -0.74000e-03_wp
            htt(22,5,4)  = -0.60000e-04_wp
            htt(22,6,4)  = -0.26200e-02_wp
            htt(22,7,4)  = 0.17100e-02_wp
            htt(22,8,4)  = -0.85000e-03_wp
            htt(22,9,4)  = -0.85000e-03_wp
            htt(22,10,4) = 0.40000e-04_wp
            htt(22,11,4) = -0.17300e-02_wp
            htt(22,12,4) = -0.85000e-03_wp
            htt(22,13,4) = 0.12800e-02_wp
            htt(22,14,4) = -0.15000e-03_wp
            htt(22,15,4) = 0.40000e-03_wp
            htt(22,16,4) = -0.10300e-02_wp
            htt(22,17,4) = -0.10700e-02_wp
            htt(22,18,4) = 0.56000e-03_wp
            htt(22,19,4) = 0.37000e-03_wp
            htt(22,20,4) = 0.20000e-02_wp
            htt(22,21,4) = 0.44000e-03_wp
            htt(22,22,4) = 0.31200e-02_wp
            htt(23,1,4)  = 0.19200e-02_wp
            htt(23,2,4)  = 0.79000e-03_wp
            htt(23,3,4)  = -0.23000e-03_wp
            htt(23,4,4)  = -0.13700e-02_wp
            htt(23,5,4)  = -0.21700e-02_wp
            htt(23,6,4)  = 0.14000e-03_wp
            htt(23,7,4)  = -0.12000e-03_wp
            htt(23,8,4)  = 0.21900e-02_wp
            htt(23,9,4)  = -0.74000e-03_wp
            htt(23,10,4) = -0.17100e-02_wp
            htt(23,11,4) = 0.45000e-03_wp
            htt(23,12,4) = -0.53000e-03_wp
            htt(23,13,4) = -0.10000e-03_wp
            htt(23,14,4) = -0.13700e-02_wp
            htt(23,15,4) = 0.87000e-03_wp
            htt(23,16,4) = -0.40000e-03_wp
            htt(23,17,4) = -0.15000e-03_wp
            htt(23,18,4) = -0.31500e-02_wp
            htt(23,19,4) = 0.11200e-02_wp
            htt(23,20,4) = -0.18800e-02_wp
            htt(23,21,4) = 0.12200e-02_wp
            htt(23,22,4) = -0.12000e-02_wp
            htt(23,23,4) = 0.42400e-02_wp
            htt(24,1,4)  = 0.30000e-04_wp
            htt(24,2,4)  = 0.42000e-03_wp
            htt(24,3,4)  = -0.84000e-03_wp
            htt(24,4,4)  = -0.45000e-03_wp
            htt(24,5,4)  = 0.79000e-03_wp
            htt(24,6,4)  = -0.18000e-02_wp
            htt(24,7,4)  = 0.16200e-02_wp
            htt(24,8,4)  = -0.97000e-03_wp
            htt(24,9,4)  = -0.39000e-03_wp
            htt(24,10,4) = -0.70000e-03_wp
            htt(24,11,4) = -0.79000e-03_wp
            htt(24,12,4) = -0.11000e-02_wp
            htt(24,13,4) = -0.80000e-04_wp
            htt(24,14,4) = -0.76000e-03_wp
            htt(24,15,4) = 0.22000e-03_wp
            htt(24,16,4) = -0.45000e-03_wp
            htt(24,17,4) = 0.92000e-03_wp
            htt(24,18,4) = -0.21000e-02_wp
            htt(24,19,4) = 0.16000e-02_wp
            htt(24,20,4) = -0.14300e-02_wp
            htt(24,21,4) = -0.12300e-02_wp
            htt(24,22,4) = 0.14800e-02_wp
            htt(24,23,4) = 0.18100e-02_wp
            htt(24,24,4) = 0.45300e-02_wp
            ! Hessian elements for anchor point 4
            hbs( 1,1,4)  = -0.31760e-01_wp
            hbs( 1,2,4)  = -0.29360e-01_wp
            hbs( 1,3,4)  = -0.31820e-01_wp
            hbs( 1,4,4)  = -0.95900e-02_wp
            hbs( 1,5,4)  = 0.52300e-02_wp
            hbs( 1,6,4)  = 0.39930e-01_wp
            hbs( 1,7,4)  = -0.45000e-03_wp
            hbs( 1,8,4)  = 0.40450e-01_wp
            hbs( 1,9,4)  = -0.18000e-03_wp
            hbs( 1,10,4) = -0.96900e-02_wp
            hbs( 1,11,4) = -0.42000e-03_wp
            hbs( 1,12,4) = 0.47400e-02_wp
            hbs( 2,1,4)  = 0.49930e-01_wp
            hbs( 2,2,4)  = -0.13930e-01_wp
            hbs( 2,3,4)  = 0.99600e-02_wp
            hbs( 2,4,4)  = 0.18000e-02_wp
            hbs( 2,5,4)  = -0.12200e-02_wp
            hbs( 2,6,4)  = -0.27810e-01_wp
            hbs( 2,7,4)  = 0.11800e-02_wp
            hbs( 2,8,4)  = -0.14650e-01_wp
            hbs( 2,9,4)  = 0.00000e+00_wp
            hbs( 2,10,4) = 0.41300e-02_wp
            hbs( 2,11,4) = -0.69000e-03_wp
            hbs( 2,12,4) = -0.89000e-03_wp
            hbs( 3,1,4)  = -0.18180e-01_wp
            hbs( 3,2,4)  = 0.43280e-01_wp
            hbs( 3,3,4)  = 0.21860e-01_wp
            hbs( 3,4,4)  = 0.77900e-02_wp
            hbs( 3,5,4)  = -0.40100e-02_wp
            hbs( 3,6,4)  = -0.12120e-01_wp
            hbs( 3,7,4)  = -0.73000e-03_wp
            hbs( 3,8,4)  = -0.25800e-01_wp
            hbs( 3,9,4)  = 0.18000e-03_wp
            hbs( 3,10,4) = 0.55600e-02_wp
            hbs( 3,11,4) = 0.11100e-02_wp
            hbs( 3,12,4) = -0.38500e-02_wp
            hbs( 4,1,4)  = -0.20820e-01_wp
            hbs( 4,2,4)  = -0.12100e-01_wp
            hbs( 4,3,4)  = 0.15180e-01_wp
            hbs( 4,4,4)  = -0.29840e-01_wp
            hbs( 4,5,4)  = -0.97300e-02_wp
            hbs( 4,6,4)  = -0.57000e-02_wp
            hbs( 4,7,4)  = 0.49500e-02_wp
            hbs( 4,8,4)  = 0.31150e-01_wp
            hbs( 4,9,4)  = -0.44000e-03_wp
            hbs( 4,10,4) = 0.40270e-01_wp
            hbs( 4,11,4) = -0.59000e-03_wp
            hbs( 4,12,4) = 0.29000e-03_wp
            hbs( 5,1,4)  = 0.26090e-01_wp
            hbs( 5,2,4)  = 0.58900e-02_wp
            hbs( 5,3,4)  = -0.14520e-01_wp
            hbs( 5,4,4)  = 0.25400e-02_wp
            hbs( 5,5,4)  = 0.69400e-02_wp
            hbs( 5,6,4)  = 0.30900e-02_wp
            hbs( 5,7,4)  = -0.21700e-02_wp
            hbs( 5,8,4)  = -0.15190e-01_wp
            hbs( 5,9,4)  = -0.39000e-03_wp
            hbs( 5,10,4) = -0.22140e-01_wp
            hbs( 5,11,4) = 0.28000e-03_wp
            hbs( 5,12,4) = 0.90000e-04_wp
            hbs( 6,1,4)  = -0.52700e-02_wp
            hbs( 6,2,4)  = 0.62100e-02_wp
            hbs( 6,3,4)  = -0.66000e-03_wp
            hbs( 6,4,4)  = 0.27300e-01_wp
            hbs( 6,5,4)  = 0.27900e-02_wp
            hbs( 6,6,4)  = 0.26100e-02_wp
            hbs( 6,7,4)  = -0.27800e-02_wp
            hbs( 6,8,4)  = -0.15960e-01_wp
            hbs( 6,9,4)  = 0.83000e-03_wp
            hbs( 6,10,4) = -0.18130e-01_wp
            hbs( 6,11,4) = 0.31000e-03_wp
            hbs( 6,12,4) = -0.38000e-03_wp
            hbs( 7,1,4)  = -0.90800e-02_wp
            hbs( 7,2,4)  = 0.35870e-01_wp
            hbs( 7,3,4)  = -0.10000e-04_wp
            hbs( 7,4,4)  = -0.22920e-01_wp
            hbs( 7,5,4)  = 0.61900e-02_wp
            hbs( 7,6,4)  = -0.28970e-01_wp
            hbs( 7,7,4)  = -0.87500e-02_wp
            hbs( 7,8,4)  = -0.92900e-02_wp
            hbs( 7,9,4)  = 0.49700e-02_wp
            hbs( 7,10,4) = 0.32780e-01_wp
            hbs( 7,11,4) = 0.60000e-04_wp
            hbs( 7,12,4) = -0.10700e-02_wp
            hbs( 8,1,4)  = 0.36900e-02_wp
            hbs( 8,2,4)  = -0.18590e-01_wp
            hbs( 8,3,4)  = 0.13800e-02_wp
            hbs( 8,4,4)  = 0.25350e-01_wp
            hbs( 8,5,4)  = -0.35000e-02_wp
            hbs( 8,6,4)  = -0.42000e-03_wp
            hbs( 8,7,4)  = 0.45000e-02_wp
            hbs( 8,8,4)  = 0.50400e-02_wp
            hbs( 8,9,4)  = -0.20100e-02_wp
            hbs( 8,10,4) = -0.14680e-01_wp
            hbs( 8,11,4) = -0.62000e-03_wp
            hbs( 8,12,4) = 0.55000e-03_wp
            hbs( 9,1,4)  = 0.53900e-02_wp
            hbs( 9,2,4)  = -0.17270e-01_wp
            hbs( 9,3,4)  = -0.13700e-02_wp
            hbs( 9,4,4)  = -0.24400e-02_wp
            hbs( 9,5,4)  = -0.26800e-02_wp
            hbs( 9,6,4)  = 0.29390e-01_wp
            hbs( 9,7,4)  = 0.42400e-02_wp
            hbs( 9,8,4)  = 0.42500e-02_wp
            hbs( 9,9,4)  = -0.29600e-02_wp
            hbs( 9,10,4) = -0.18100e-01_wp
            hbs( 9,11,4) = 0.55000e-03_wp
            hbs( 9,12,4) = 0.51000e-03_wp
            hbs(10,1,4)  = 0.36130e-01_wp
            hbs(10,2,4)  = 0.38360e-01_wp
            hbs(10,3,4)  = 0.36300e-02_wp
            hbs(10,4,4)  = -0.11890e-01_wp
            hbs(10,5,4)  = -0.10500e-02_wp
            hbs(10,6,4)  = -0.26840e-01_wp
            hbs(10,7,4)  = 0.48100e-02_wp
            hbs(10,8,4)  = -0.27100e-01_wp
            hbs(10,9,4)  = -0.88200e-02_wp
            hbs(10,10,4) = -0.12140e-01_wp
            hbs(10,11,4) = 0.47600e-02_wp
            hbs(10,12,4) = -0.77000e-03_wp
            hbs(11,1,4)  = -0.20190e-01_wp
            hbs(11,2,4)  = -0.17180e-01_wp
            hbs(11,3,4)  = -0.21300e-02_wp
            hbs(11,4,4)  = 0.57400e-02_wp
            hbs(11,5,4)  = 0.12900e-02_wp
            hbs(11,6,4)  = 0.27000e-01_wp
            hbs(11,7,4)  = -0.28400e-02_wp
            hbs(11,8,4)  = 0.47000e-03_wp
            hbs(11,9,4)  = 0.44300e-02_wp
            hbs(11,10,4) = 0.60700e-02_wp
            hbs(11,11,4) = -0.19800e-02_wp
            hbs(11,12,4) = -0.18000e-03_wp
            hbs(12,1,4)  = -0.15940e-01_wp
            hbs(12,2,4)  = -0.21180e-01_wp
            hbs(12,3,4)  = -0.14900e-02_wp
            hbs(12,4,4)  = 0.61500e-02_wp
            hbs(12,5,4)  = -0.24000e-03_wp
            hbs(12,6,4)  = -0.16000e-03_wp
            hbs(12,7,4)  = -0.19600e-02_wp
            hbs(12,8,4)  = 0.26630e-01_wp
            hbs(12,9,4)  = 0.43900e-02_wp
            hbs(12,10,4) = 0.60700e-02_wp
            hbs(12,11,4) = -0.27800e-02_wp
            hbs(12,12,4) = 0.94000e-03_wp
            hbs(13,1,4)  = 0.35720e-01_wp
            hbs(13,2,4)  = -0.96900e-02_wp
            hbs(13,3,4)  = 0.29800e-02_wp
            hbs(13,4,4)  = 0.32310e-01_wp
            hbs(13,5,4)  = -0.94000e-03_wp
            hbs(13,6,4)  = -0.10140e-01_wp
            hbs(13,7,4)  = 0.70000e-04_wp
            hbs(13,8,4)  = -0.29530e-01_wp
            hbs(13,9,4)  = 0.50200e-02_wp
            hbs(13,10,4) = -0.23420e-01_wp
            hbs(13,11,4) = -0.86600e-02_wp
            hbs(13,12,4) = 0.53400e-02_wp
            hbs(14,1,4)  = -0.17210e-01_wp
            hbs(14,2,4)  = 0.56100e-02_wp
            hbs(14,3,4)  = -0.23300e-02_wp
            hbs(14,4,4)  = -0.17680e-01_wp
            hbs(14,5,4)  = 0.48000e-03_wp
            hbs(14,6,4)  = 0.47700e-02_wp
            hbs(14,7,4)  = 0.56000e-03_wp
            hbs(14,8,4)  = 0.29630e-01_wp
            hbs(14,9,4)  = -0.29700e-02_wp
            hbs(14,10,4) = -0.24000e-02_wp
            hbs(14,11,4) = 0.43000e-02_wp
            hbs(14,12,4) = -0.22600e-02_wp
            hbs(15,1,4)  = -0.18510e-01_wp
            hbs(15,2,4)  = 0.40800e-02_wp
            hbs(15,3,4)  = -0.65000e-03_wp
            hbs(15,4,4)  = -0.14630e-01_wp
            hbs(15,5,4)  = 0.46000e-03_wp
            hbs(15,6,4)  = 0.53700e-02_wp
            hbs(15,7,4)  = -0.63000e-03_wp
            hbs(15,8,4)  = -0.10000e-03_wp
            hbs(15,9,4)  = -0.20500e-02_wp
            hbs(15,10,4) = 0.25820e-01_wp
            hbs(15,11,4) = 0.43700e-02_wp
            hbs(15,12,4) = -0.30900e-02_wp
            hbs(16,1,4)  = -0.10200e-01_wp
            hbs(16,2,4)  = -0.23080e-01_wp
            hbs(16,3,4)  = 0.10040e-01_wp
            hbs(16,4,4)  = 0.41930e-01_wp
            hbs(16,5,4)  = 0.31000e-03_wp
            hbs(16,6,4)  = 0.31710e-01_wp
            hbs(16,7,4)  = -0.63000e-03_wp
            hbs(16,8,4)  = -0.56800e-02_wp
            hbs(16,9,4)  = -0.55000e-03_wp
            hbs(16,10,4) = -0.27810e-01_wp
            hbs(16,11,4) = 0.48500e-02_wp
            hbs(16,12,4) = -0.85400e-02_wp
            hbs(17,1,4)  = 0.41100e-02_wp
            hbs(17,2,4)  = 0.25000e-01_wp
            hbs(17,3,4)  = -0.74400e-02_wp
            hbs(17,4,4)  = -0.22980e-01_wp
            hbs(17,5,4)  = 0.60000e-03_wp
            hbs(17,6,4)  = -0.14860e-01_wp
            hbs(17,7,4)  = 0.32000e-03_wp
            hbs(17,8,4)  = 0.33600e-02_wp
            hbs(17,9,4)  = -0.32000e-03_wp
            hbs(17,10,4) = 0.20400e-02_wp
            hbs(17,11,4) = -0.21000e-02_wp
            hbs(17,12,4) = 0.68800e-02_wp
            hbs(18,1,4)  = 0.61000e-02_wp
            hbs(18,2,4)  = -0.19200e-02_wp
            hbs(18,3,4)  = -0.26000e-02_wp
            hbs(18,4,4)  = -0.18950e-01_wp
            hbs(18,5,4)  = -0.92000e-03_wp
            hbs(18,6,4)  = -0.16850e-01_wp
            hbs(18,7,4)  = 0.31000e-03_wp
            hbs(18,8,4)  = 0.23200e-02_wp
            hbs(18,9,4)  = 0.87000e-03_wp
            hbs(18,10,4) = 0.25770e-01_wp
            hbs(18,11,4) = -0.27500e-02_wp
            hbs(18,12,4) = 0.16600e-02_wp
        elseif( ii .eq. 33 )then
            ! Hessian elements for anchor point 1
            hss( 1,1,1)  = 0.32935e+00_wp
            hss( 2,1,1)  = 0.10630e-01_wp
            hss( 2,2,1)  = 0.33545e+00_wp
            hss( 3,1,1)  = 0.41480e-01_wp
            hss( 3,2,1)  = 0.47860e-01_wp
            hss( 3,3,1)  = 0.60075e+00_wp
            hss( 4,1,1)  = -0.28200e-02_wp
            hss( 4,2,1)  = 0.17340e-01_wp
            hss( 4,3,1)  = -0.51100e-02_wp
            hss( 4,4,1)  = 0.44790e+00_wp
            hss( 5,1,1)  = 0.39700e-02_wp
            hss( 5,2,1)  = -0.41000e-03_wp
            hss( 5,3,1)  = 0.20000e-03_wp
            hss( 5,4,1)  = 0.13100e-02_wp
            hss( 5,5,1)  = 0.36210e+00_wp
            hss( 6,1,1)  = 0.23160e-01_wp
            hss( 6,2,1)  = 0.46150e-01_wp
            hss( 6,3,1)  = 0.31100e-02_wp
            hss( 6,4,1)  = 0.42160e-01_wp
            hss( 6,5,1)  = -0.19000e-03_wp
            hss( 6,6,1)  = 0.33957e+00_wp
            hss( 7,1,1)  = 0.14000e-03_wp
            hss( 7,2,1)  = -0.96000e-03_wp
            hss( 7,3,1)  = 0.11700e-02_wp
            hss( 7,4,1)  = 0.55900e-02_wp
            hss( 7,5,1)  = -0.68000e-03_wp
            hss( 7,6,1)  = 0.28700e-02_wp
            hss( 7,7,1)  = 0.36718e+00_wp
            hss( 8,1,1)  = 0.48190e-01_wp
            hss( 8,2,1)  = 0.21300e-01_wp
            hss( 8,3,1)  = 0.10200e-02_wp
            hss( 8,4,1)  = -0.55600e-02_wp
            hss( 8,5,1)  = -0.21100e-02_wp
            hss( 8,6,1)  = 0.34440e-01_wp
            hss( 8,7,1)  = -0.39000e-03_wp
            hss( 8,8,1)  = 0.33723e+00_wp
            hss( 9,1,1)  = -0.13000e-03_wp
            hss( 9,2,1)  = 0.50000e-04_wp
            hss( 9,3,1)  = -0.53000e-03_wp
            hss( 9,4,1)  = -0.54000e-03_wp
            hss( 9,5,1)  = -0.64000e-03_wp
            hss( 9,6,1)  = 0.37800e-02_wp
            hss( 9,7,1)  = 0.12000e-03_wp
            hss( 9,8,1)  = 0.36400e-02_wp
            hss( 9,9,1)  = 0.36496e+00_wp
            hss(10,1,1)  = 0.14800e-01_wp
            hss(10,2,1)  = 0.14000e-03_wp
            hss(10,3,1)  = -0.42300e-02_wp
            hss(10,4,1)  = 0.81590e-01_wp
            hss(10,5,1)  = -0.48000e-03_wp
            hss(10,6,1)  = -0.32700e-02_wp
            hss(10,7,1)  = -0.71000e-03_wp
            hss(10,8,1)  = 0.36040e-01_wp
            hss(10,9,1)  = -0.87000e-03_wp
            hss(10,10,1) = 0.46416e+00_wp
            hss(11,1,1)  = -0.11800e-02_wp
            hss(11,2,1)  = -0.64000e-03_wp
            hss(11,3,1)  = 0.90000e-03_wp
            hss(11,4,1)  = -0.78000e-03_wp
            hss(11,5,1)  = 0.18000e-03_wp
            hss(11,6,1)  = -0.11000e-03_wp
            hss(11,7,1)  = 0.12000e-03_wp
            hss(11,8,1)  = 0.28600e-02_wp
            hss(11,9,1)  = 0.42000e-03_wp
            hss(11,10,1) = 0.45200e-02_wp
            hss(11,11,1) = 0.36758e+00_wp
            hss(12,1,1)  = 0.23000e-03_wp
            hss(12,2,1)  = 0.14400e-02_wp
            hss(12,3,1)  = 0.13000e-03_wp
            hss(12,4,1)  = -0.34000e-03_wp
            hss(12,5,1)  = -0.37000e-03_wp
            hss(12,6,1)  = -0.80000e-03_wp
            hss(12,7,1)  = 0.40000e-04_wp
            hss(12,8,1)  = 0.22000e-03_wp
            hss(12,9,1)  = -0.10000e-03_wp
            hss(12,10,1) = 0.23100e-02_wp
            hss(12,11,1) = 0.17000e-03_wp
            hss(12,12,1) = 0.37062e+00_wp
            ! Hessian elements for anchor point 1
            hbb( 1,1,1)  = 0.82090e-01_wp
            hbb( 2,1,1)  = -0.40800e-01_wp
            hbb( 2,2,1)  = 0.13448e+00_wp
            hbb( 3,1,1)  = -0.41290e-01_wp
            hbb( 3,2,1)  = -0.93680e-01_wp
            hbb( 3,3,1)  = 0.13497e+00_wp
            hbb( 4,1,1)  = -0.37560e-01_wp
            hbb( 4,2,1)  = 0.25500e-01_wp
            hbb( 4,3,1)  = 0.12070e-01_wp
            hbb( 4,4,1)  = 0.76550e-01_wp
            hbb( 5,1,1)  = 0.22450e-01_wp
            hbb( 5,2,1)  = -0.15660e-01_wp
            hbb( 5,3,1)  = -0.67900e-02_wp
            hbb( 5,4,1)  = -0.39660e-01_wp
            hbb( 5,5,1)  = 0.77760e-01_wp
            hbb( 6,1,1)  = 0.15120e-01_wp
            hbb( 6,2,1)  = -0.98400e-02_wp
            hbb( 6,3,1)  = -0.52700e-02_wp
            hbb( 6,4,1)  = -0.36880e-01_wp
            hbb( 6,5,1)  = -0.38100e-01_wp
            hbb( 6,6,1)  = 0.74990e-01_wp
            hbb( 7,1,1)  = -0.83500e-02_wp
            hbb( 7,2,1)  = -0.60600e-02_wp
            hbb( 7,3,1)  = 0.14410e-01_wp
            hbb( 7,4,1)  = -0.28860e-01_wp
            hbb( 7,5,1)  = 0.12030e-01_wp
            hbb( 7,6,1)  = 0.16840e-01_wp
            hbb( 7,7,1)  = 0.74530e-01_wp
            hbb( 8,1,1)  = -0.15000e-02_wp
            hbb( 8,2,1)  = 0.75600e-02_wp
            hbb( 8,3,1)  = -0.60500e-02_wp
            hbb( 8,4,1)  = 0.17370e-01_wp
            hbb( 8,5,1)  = -0.66700e-02_wp
            hbb( 8,6,1)  = -0.10710e-01_wp
            hbb( 8,7,1)  = -0.36970e-01_wp
            hbb( 8,8,1)  = 0.78530e-01_wp
            hbb( 9,1,1)  = 0.98600e-02_wp
            hbb( 9,2,1)  = -0.15000e-02_wp
            hbb( 9,3,1)  = -0.83500e-02_wp
            hbb( 9,4,1)  = 0.11490e-01_wp
            hbb( 9,5,1)  = -0.53600e-02_wp
            hbb( 9,6,1)  = -0.61300e-02_wp
            hbb( 9,7,1)  = -0.37550e-01_wp
            hbb( 9,8,1)  = -0.41560e-01_wp
            hbb( 9,9,1)  = 0.79110e-01_wp
            hbb(10,1,1)  = 0.88000e-02_wp
            hbb(10,2,1)  = -0.45300e-02_wp
            hbb(10,3,1)  = -0.42700e-02_wp
            hbb(10,4,1)  = -0.86000e-02_wp
            hbb(10,5,1)  = 0.10210e-01_wp
            hbb(10,6,1)  = -0.16100e-02_wp
            hbb(10,7,1)  = -0.33520e-01_wp
            hbb(10,8,1)  = 0.13620e-01_wp
            hbb(10,9,1)  = 0.19900e-01_wp
            hbb(10,10,1) = 0.75620e-01_wp
            hbb(11,1,1)  = -0.44200e-02_wp
            hbb(11,2,1)  = 0.26700e-02_wp
            hbb(11,3,1)  = 0.17500e-02_wp
            hbb(11,4,1)  = -0.13500e-02_wp
            hbb(11,5,1)  = -0.34500e-02_wp
            hbb(11,6,1)  = 0.48000e-02_wp
            hbb(11,7,1)  = 0.19740e-01_wp
            hbb(11,8,1)  = -0.73600e-02_wp
            hbb(11,9,1)  = -0.12380e-01_wp
            hbb(11,10,1) = -0.37740e-01_wp
            hbb(11,11,1) = 0.80040e-01_wp
            hbb(12,1,1)  = -0.43800e-02_wp
            hbb(12,2,1)  = 0.18600e-02_wp
            hbb(12,3,1)  = 0.25200e-02_wp
            hbb(12,4,1)  = 0.99500e-02_wp
            hbb(12,5,1)  = -0.67600e-02_wp
            hbb(12,6,1)  = -0.31900e-02_wp
            hbb(12,7,1)  = 0.13790e-01_wp
            hbb(12,8,1)  = -0.62600e-02_wp
            hbb(12,9,1)  = -0.75300e-02_wp
            hbb(12,10,1) = -0.37880e-01_wp
            hbb(12,11,1) = -0.42300e-01_wp
            hbb(12,12,1) = 0.80170e-01_wp
            hbb(13,1,1)  = -0.76900e-02_wp
            hbb(13,2,1)  = 0.14150e-01_wp
            hbb(13,3,1)  = -0.64500e-02_wp
            hbb(13,4,1)  = 0.32300e-02_wp
            hbb(13,5,1)  = -0.14400e-02_wp
            hbb(13,6,1)  = -0.17900e-02_wp
            hbb(13,7,1)  = -0.67100e-02_wp
            hbb(13,8,1)  = 0.90400e-02_wp
            hbb(13,9,1)  = -0.23400e-02_wp
            hbb(13,10,1) = -0.33330e-01_wp
            hbb(13,11,1) = 0.13660e-01_wp
            hbb(13,12,1) = 0.19670e-01_wp
            hbb(13,13,1) = 0.73530e-01_wp
            hbb(14,1,1)  = 0.96900e-02_wp
            hbb(14,2,1)  = -0.81500e-02_wp
            hbb(14,3,1)  = -0.15500e-02_wp
            hbb(14,4,1)  = -0.17500e-02_wp
            hbb(14,5,1)  = 0.64000e-03_wp
            hbb(14,6,1)  = 0.11000e-02_wp
            hbb(14,7,1)  = -0.22600e-02_wp
            hbb(14,8,1)  = -0.29800e-02_wp
            hbb(14,9,1)  = 0.52500e-02_wp
            hbb(14,10,1) = 0.19960e-01_wp
            hbb(14,11,1) = -0.75600e-02_wp
            hbb(14,12,1) = -0.12400e-01_wp
            hbb(14,13,1) = -0.37280e-01_wp
            hbb(14,14,1) = 0.78950e-01_wp
            hbb(15,1,1)  = -0.20000e-02_wp
            hbb(15,2,1)  = -0.60000e-02_wp
            hbb(15,3,1)  = 0.80000e-02_wp
            hbb(15,4,1)  = -0.14800e-02_wp
            hbb(15,5,1)  = 0.80000e-03_wp
            hbb(15,6,1)  = 0.69000e-03_wp
            hbb(15,7,1)  = 0.89700e-02_wp
            hbb(15,8,1)  = -0.60600e-02_wp
            hbb(15,9,1)  = -0.29100e-02_wp
            hbb(15,10,1) = 0.13370e-01_wp
            hbb(15,11,1) = -0.61000e-02_wp
            hbb(15,12,1) = -0.72700e-02_wp
            hbb(15,13,1) = -0.36250e-01_wp
            hbb(15,14,1) = -0.41670e-01_wp
            hbb(15,15,1) = 0.77920e-01_wp
            hbb(16,1,1)  = -0.37280e-01_wp
            hbb(16,2,1)  = 0.11740e-01_wp
            hbb(16,3,1)  = 0.25540e-01_wp
            hbb(16,4,1)  = -0.47500e-02_wp
            hbb(16,5,1)  = -0.35800e-02_wp
            hbb(16,6,1)  = 0.83300e-02_wp
            hbb(16,7,1)  = 0.29200e-02_wp
            hbb(16,8,1)  = -0.15600e-02_wp
            hbb(16,9,1)  = -0.13600e-02_wp
            hbb(16,10,1) = -0.89600e-02_wp
            hbb(16,11,1) = 0.10110e-01_wp
            hbb(16,12,1) = -0.11500e-02_wp
            hbb(16,13,1) = -0.29030e-01_wp
            hbb(16,14,1) = 0.11640e-01_wp
            hbb(16,15,1) = 0.17390e-01_wp
            hbb(16,16,1) = 0.77110e-01_wp
            hbb(17,1,1)  = 0.22780e-01_wp
            hbb(17,2,1)  = -0.48100e-02_wp
            hbb(17,3,1)  = -0.17970e-01_wp
            hbb(17,4,1)  = -0.38000e-02_wp
            hbb(17,5,1)  = 0.60900e-02_wp
            hbb(17,6,1)  = -0.22900e-02_wp
            hbb(17,7,1)  = -0.14200e-02_wp
            hbb(17,8,1)  = 0.93000e-03_wp
            hbb(17,9,1)  = 0.49000e-03_wp
            hbb(17,10,1) = 0.10430e-01_wp
            hbb(17,11,1) = -0.67600e-02_wp
            hbb(17,12,1) = -0.36800e-02_wp
            hbb(17,13,1) = 0.11310e-01_wp
            hbb(17,14,1) = -0.50800e-02_wp
            hbb(17,15,1) = -0.62300e-02_wp
            hbb(17,16,1) = -0.39310e-01_wp
            hbb(17,17,1) = 0.76690e-01_wp
            hbb(18,1,1)  = 0.14500e-01_wp
            hbb(18,2,1)  = -0.69300e-02_wp
            hbb(18,3,1)  = -0.75700e-02_wp
            hbb(18,4,1)  = 0.85400e-02_wp
            hbb(18,5,1)  = -0.25100e-02_wp
            hbb(18,6,1)  = -0.60300e-02_wp
            hbb(18,7,1)  = -0.15000e-02_wp
            hbb(18,8,1)  = 0.63000e-03_wp
            hbb(18,9,1)  = 0.87000e-03_wp
            hbb(18,10,1) = -0.14700e-02_wp
            hbb(18,11,1) = -0.33600e-02_wp
            hbb(18,12,1) = 0.48300e-02_wp
            hbb(18,13,1) = 0.17720e-01_wp
            hbb(18,14,1) = -0.65600e-02_wp
            hbb(18,15,1) = -0.11160e-01_wp
            hbb(18,16,1) = -0.37800e-01_wp
            hbb(18,17,1) = -0.37380e-01_wp
            hbb(18,18,1) = 0.75170e-01_wp
            ! Hessian elements for anchor point 1
            htt( 1,1,1)  = 0.12410e-01_wp
            htt( 2,1,1)  = 0.66300e-02_wp
            htt( 2,2,1)  = 0.11350e-01_wp
            htt( 3,1,1)  = -0.60000e-02_wp
            htt( 3,2,1)  = -0.92200e-02_wp
            htt( 3,3,1)  = 0.16800e-01_wp
            htt( 4,1,1)  = -0.11780e-01_wp
            htt( 4,2,1)  = -0.45100e-02_wp
            htt( 4,3,1)  = 0.13580e-01_wp
            htt( 4,4,1)  = 0.20850e-01_wp
            htt( 5,1,1)  = -0.10450e-01_wp
            htt( 5,2,1)  = -0.82200e-02_wp
            htt( 5,3,1)  = 0.95200e-02_wp
            htt( 5,4,1)  = 0.11750e-01_wp
            htt( 5,5,1)  = 0.13600e-01_wp
            htt( 6,1,1)  = -0.73200e-02_wp
            htt( 6,2,1)  = -0.56100e-02_wp
            htt( 6,3,1)  = 0.79400e-02_wp
            htt( 6,4,1)  = 0.96500e-02_wp
            htt( 6,5,1)  = 0.70500e-02_wp
            htt( 6,6,1)  = 0.10410e-01_wp
            htt( 7,1,1)  = 0.68700e-02_wp
            htt( 7,2,1)  = 0.66900e-02_wp
            htt( 7,3,1)  = -0.11920e-01_wp
            htt( 7,4,1)  = -0.12100e-01_wp
            htt( 7,5,1)  = -0.51800e-02_wp
            htt( 7,6,1)  = -0.73000e-02_wp
            htt( 7,7,1)  = 0.12490e-01_wp
            htt( 8,1,1)  = 0.10000e-01_wp
            htt( 8,2,1)  = 0.93000e-02_wp
            htt( 8,3,1)  = -0.13500e-01_wp
            htt( 8,4,1)  = -0.14200e-01_wp
            htt( 8,5,1)  = -0.11730e-01_wp
            htt( 8,6,1)  = -0.39400e-02_wp
            htt( 8,7,1)  = 0.10370e-01_wp
            htt( 8,8,1)  = 0.18160e-01_wp
            htt( 9,1,1)  = -0.62000e-02_wp
            htt( 9,2,1)  = -0.56000e-03_wp
            htt( 9,3,1)  = -0.18500e-02_wp
            htt( 9,4,1)  = 0.37900e-02_wp
            htt( 9,5,1)  = 0.12100e-02_wp
            htt( 9,6,1)  = 0.24300e-02_wp
            htt( 9,7,1)  = -0.28800e-02_wp
            htt( 9,8,1)  = -0.16600e-02_wp
            htt( 9,9,1)  = 0.82000e-02_wp
            htt(10,1,1)  = -0.49900e-02_wp
            htt(10,2,1)  = -0.30000e-04_wp
            htt(10,3,1)  = 0.60000e-03_wp
            htt(10,4,1)  = 0.55600e-02_wp
            htt(10,5,1)  = 0.33700e-02_wp
            htt(10,6,1)  = 0.24600e-02_wp
            htt(10,7,1)  = -0.18800e-02_wp
            htt(10,8,1)  = -0.27900e-02_wp
            htt(10,9,1)  = 0.31400e-02_wp
            htt(10,10,1) = 0.89300e-02_wp
            htt(11,1,1)  = -0.90000e-04_wp
            htt(11,2,1)  = -0.55400e-02_wp
            htt(11,3,1)  = 0.15600e-02_wp
            htt(11,4,1)  = -0.39000e-02_wp
            htt(11,5,1)  = -0.11400e-02_wp
            htt(11,6,1)  = 0.63000e-03_wp
            htt(11,7,1)  = -0.26900e-02_wp
            htt(11,8,1)  = -0.91000e-03_wp
            htt(11,9,1)  = 0.22400e-02_wp
            htt(11,10,1) = -0.21000e-02_wp
            htt(11,11,1) = 0.80000e-02_wp
            htt(12,1,1)  = 0.11300e-02_wp
            htt(12,2,1)  = -0.50100e-02_wp
            htt(12,3,1)  = 0.40000e-02_wp
            htt(12,4,1)  = -0.21300e-02_wp
            htt(12,5,1)  = 0.10100e-02_wp
            htt(12,6,1)  = 0.67000e-03_wp
            htt(12,7,1)  = -0.17000e-02_wp
            htt(12,8,1)  = -0.20400e-02_wp
            htt(12,9,1)  = -0.28200e-02_wp
            htt(12,10,1) = 0.36900e-02_wp
            htt(12,11,1) = 0.36700e-02_wp
            htt(12,12,1) = 0.10180e-01_wp
            htt(13,1,1)  = -0.19600e-02_wp
            htt(13,2,1)  = -0.40500e-02_wp
            htt(13,3,1)  = 0.64000e-02_wp
            htt(13,4,1)  = 0.43100e-02_wp
            htt(13,5,1)  = 0.50600e-02_wp
            htt(13,6,1)  = 0.27500e-02_wp
            htt(13,7,1)  = -0.28000e-02_wp
            htt(13,8,1)  = -0.51000e-02_wp
            htt(13,9,1)  = -0.54600e-02_wp
            htt(13,10,1) = 0.27000e-03_wp
            htt(13,11,1) = -0.32400e-02_wp
            htt(13,12,1) = 0.24800e-02_wp
            htt(13,13,1) = 0.10020e-01_wp
            htt(14,1,1)  = 0.73000e-03_wp
            htt(14,2,1)  = -0.74000e-03_wp
            htt(14,3,1)  = 0.28700e-02_wp
            htt(14,4,1)  = 0.14100e-02_wp
            htt(14,5,1)  = 0.19900e-02_wp
            htt(14,6,1)  = -0.24000e-03_wp
            htt(14,7,1)  = -0.30000e-04_wp
            htt(14,8,1)  = -0.22600e-02_wp
            htt(14,9,1)  = -0.28600e-02_wp
            htt(14,10,1) = 0.17000e-02_wp
            htt(14,11,1) = -0.13100e-02_wp
            htt(14,12,1) = 0.32500e-02_wp
            htt(14,13,1) = 0.23900e-02_wp
            htt(14,14,1) = 0.74100e-02_wp
            htt(15,1,1)  = -0.31700e-02_wp
            htt(15,2,1)  = -0.45800e-02_wp
            htt(15,3,1)  = 0.39600e-02_wp
            htt(15,4,1)  = 0.25500e-02_wp
            htt(15,5,1)  = 0.29100e-02_wp
            htt(15,6,1)  = 0.27200e-02_wp
            htt(15,7,1)  = -0.37900e-02_wp
            htt(15,8,1)  = -0.39800e-02_wp
            htt(15,9,1)  = -0.41000e-03_wp
            htt(15,10,1) = -0.55000e-02_wp
            htt(15,11,1) = 0.10800e-02_wp
            htt(15,12,1) = -0.40100e-02_wp
            htt(15,13,1) = 0.43100e-02_wp
            htt(15,14,1) = -0.21600e-02_wp
            htt(15,15,1) = 0.93900e-02_wp
            htt(16,1,1)  = -0.48000e-03_wp
            htt(16,2,1)  = -0.12700e-02_wp
            htt(16,3,1)  = 0.43000e-03_wp
            htt(16,4,1)  = -0.35000e-03_wp
            htt(16,5,1)  = -0.16000e-03_wp
            htt(16,6,1)  = -0.28000e-03_wp
            htt(16,7,1)  = -0.10200e-02_wp
            htt(16,8,1)  = -0.11300e-02_wp
            htt(16,9,1)  = 0.21900e-02_wp
            htt(16,10,1) = -0.40700e-02_wp
            htt(16,11,1) = 0.30200e-02_wp
            htt(16,12,1) = -0.32400e-02_wp
            htt(16,13,1) = -0.33300e-02_wp
            htt(16,14,1) = 0.28700e-02_wp
            htt(16,15,1) = 0.29200e-02_wp
            htt(16,16,1) = 0.91100e-02_wp
            htt(17,1,1)  = 0.40200e-02_wp
            htt(17,2,1)  = 0.24400e-02_wp
            htt(17,3,1)  = -0.28100e-02_wp
            htt(17,4,1)  = -0.43900e-02_wp
            htt(17,5,1)  = -0.18700e-02_wp
            htt(17,6,1)  = -0.30500e-02_wp
            htt(17,7,1)  = 0.45500e-02_wp
            htt(17,8,1)  = 0.33700e-02_wp
            htt(17,9,1)  = 0.30000e-03_wp
            htt(17,10,1) = -0.19500e-02_wp
            htt(17,11,1) = 0.19700e-02_wp
            htt(17,12,1) = -0.28000e-03_wp
            htt(17,13,1) = -0.68200e-02_wp
            htt(17,14,1) = 0.40000e-03_wp
            htt(17,15,1) = -0.45700e-02_wp
            htt(17,16,1) = 0.26500e-02_wp
            htt(17,17,1) = 0.90000e-02_wp
            htt(18,1,1)  = 0.19000e-02_wp
            htt(18,2,1)  = 0.26000e-02_wp
            htt(18,3,1)  = -0.39700e-02_wp
            htt(18,4,1)  = -0.32700e-02_wp
            htt(18,5,1)  = -0.30200e-02_wp
            htt(18,6,1)  = -0.37000e-02_wp
            htt(18,7,1)  = 0.25000e-02_wp
            htt(18,8,1)  = 0.18200e-02_wp
            htt(18,9,1)  = 0.19300e-02_wp
            htt(18,10,1) = 0.43000e-03_wp
            htt(18,11,1) = 0.11900e-02_wp
            htt(18,12,1) = -0.31000e-03_wp
            htt(18,13,1) = -0.47600e-02_wp
            htt(18,14,1) = 0.15300e-02_wp
            htt(18,15,1) = -0.32600e-02_wp
            htt(18,16,1) = 0.30200e-02_wp
            htt(18,17,1) = 0.36100e-02_wp
            htt(18,18,1) = 0.88300e-02_wp
            htt(19,1,1)  = 0.13400e-02_wp
            htt(19,2,1)  = -0.86000e-03_wp
            htt(19,3,1)  = 0.71000e-03_wp
            htt(19,4,1)  = -0.14900e-02_wp
            htt(19,5,1)  = 0.11900e-02_wp
            htt(19,6,1)  = -0.60000e-04_wp
            htt(19,7,1)  = 0.17900e-02_wp
            htt(19,8,1)  = 0.54000e-03_wp
            htt(19,9,1)  = -0.23000e-02_wp
            htt(19,10,1) = -0.33800e-02_wp
            htt(19,11,1) = 0.40000e-04_wp
            htt(19,12,1) = -0.10500e-02_wp
            htt(19,13,1) = 0.81000e-03_wp
            htt(19,14,1) = -0.46200e-02_wp
            htt(19,15,1) = 0.18900e-02_wp
            htt(19,16,1) = -0.35300e-02_wp
            htt(19,17,1) = 0.17900e-02_wp
            htt(19,18,1) = -0.26700e-02_wp
            htt(19,19,1) = 0.72000e-02_wp
            htt(20,1,1)  = -0.78000e-03_wp
            htt(20,2,1)  = -0.70000e-03_wp
            htt(20,3,1)  = -0.45000e-03_wp
            htt(20,4,1)  = -0.37000e-03_wp
            htt(20,5,1)  = 0.40000e-04_wp
            htt(20,6,1)  = -0.71000e-03_wp
            htt(20,7,1)  = -0.27000e-03_wp
            htt(20,8,1)  = -0.10200e-02_wp
            htt(20,9,1)  = -0.66000e-03_wp
            htt(20,10,1) = -0.10000e-02_wp
            htt(20,11,1) = -0.74000e-03_wp
            htt(20,12,1) = -0.10700e-02_wp
            htt(20,13,1) = 0.28600e-02_wp
            htt(20,14,1) = -0.34900e-02_wp
            htt(20,15,1) = 0.32000e-02_wp
            htt(20,16,1) = -0.31600e-02_wp
            htt(20,17,1) = -0.36000e-02_wp
            htt(20,18,1) = 0.25500e-02_wp
            htt(20,19,1) = 0.27400e-02_wp
            htt(20,20,1) = 0.89000e-02_wp
            htt(21,1,1)  = 0.22500e-02_wp
            htt(21,2,1)  = 0.36500e-02_wp
            htt(21,3,1)  = -0.50600e-02_wp
            htt(21,4,1)  = -0.36600e-02_wp
            htt(21,5,1)  = -0.73600e-02_wp
            htt(21,6,1)  = -0.18700e-02_wp
            htt(21,7,1)  = -0.49000e-03_wp
            htt(21,8,1)  = 0.50000e-02_wp
            htt(21,9,1)  = 0.16800e-02_wp
            htt(21,10,1) = 0.80000e-04_wp
            htt(21,11,1) = 0.20000e-03_wp
            htt(21,12,1) = -0.14000e-02_wp
            htt(21,13,1) = -0.68000e-03_wp
            htt(21,14,1) = -0.25100e-02_wp
            htt(21,15,1) = 0.92000e-03_wp
            htt(21,16,1) = -0.91000e-03_wp
            htt(21,17,1) = -0.45300e-02_wp
            htt(21,18,1) = 0.28000e-03_wp
            htt(21,19,1) = -0.27000e-02_wp
            htt(21,20,1) = 0.21000e-02_wp
            htt(21,21,1) = 0.83700e-02_wp
            htt(22,1,1)  = -0.10700e-02_wp
            htt(22,2,1)  = 0.89000e-03_wp
            htt(22,3,1)  = -0.33900e-02_wp
            htt(22,4,1)  = -0.14300e-02_wp
            htt(22,5,1)  = -0.43000e-03_wp
            htt(22,6,1)  = -0.54300e-02_wp
            htt(22,7,1)  = 0.17600e-02_wp
            htt(22,8,1)  = -0.32400e-02_wp
            htt(22,9,1)  = 0.39000e-03_wp
            htt(22,10,1) = 0.10300e-02_wp
            htt(22,11,1) = -0.16800e-02_wp
            htt(22,12,1) = -0.10400e-02_wp
            htt(22,13,1) = 0.17600e-02_wp
            htt(22,14,1) = -0.15000e-03_wp
            htt(22,15,1) = 0.11200e-02_wp
            htt(22,16,1) = -0.78000e-03_wp
            htt(22,17,1) = -0.32800e-02_wp
            htt(22,18,1) = 0.99000e-03_wp
            htt(22,19,1) = -0.13800e-02_wp
            htt(22,20,1) = 0.29000e-02_wp
            htt(22,21,1) = 0.25500e-02_wp
            htt(22,22,1) = 0.78500e-02_wp
            htt(23,1,1)  = 0.43900e-02_wp
            htt(23,2,1)  = 0.34900e-02_wp
            htt(23,3,1)  = -0.38900e-02_wp
            htt(23,4,1)  = -0.47900e-02_wp
            htt(23,5,1)  = -0.62100e-02_wp
            htt(23,6,1)  = -0.12200e-02_wp
            htt(23,7,1)  = 0.15800e-02_wp
            htt(23,8,1)  = 0.65700e-02_wp
            htt(23,9,1)  = 0.40000e-04_wp
            htt(23,10,1) = -0.23200e-02_wp
            htt(23,11,1) = 0.98000e-03_wp
            htt(23,12,1) = -0.13800e-02_wp
            htt(23,13,1) = -0.27500e-02_wp
            htt(23,14,1) = -0.36400e-02_wp
            htt(23,15,1) = -0.39000e-03_wp
            htt(23,16,1) = -0.12800e-02_wp
            htt(23,17,1) = 0.91000e-03_wp
            htt(23,18,1) = -0.49900e-02_wp
            htt(23,19,1) = 0.18000e-02_wp
            htt(23,20,1) = -0.41000e-02_wp
            htt(23,21,1) = 0.35200e-02_wp
            htt(23,22,1) = -0.17600e-02_wp
            htt(23,23,1) = 0.94700e-02_wp
            htt(24,1,1)  = 0.10700e-02_wp
            htt(24,2,1)  = 0.73000e-03_wp
            htt(24,3,1)  = -0.22200e-02_wp
            htt(24,4,1)  = -0.25600e-02_wp
            htt(24,5,1)  = 0.73000e-03_wp
            htt(24,6,1)  = -0.47800e-02_wp
            htt(24,7,1)  = 0.38300e-02_wp
            htt(24,8,1)  = -0.16800e-02_wp
            htt(24,9,1)  = -0.12500e-02_wp
            htt(24,10,1) = -0.13700e-02_wp
            htt(24,11,1) = -0.90000e-03_wp
            htt(24,12,1) = -0.10100e-02_wp
            htt(24,13,1) = -0.31000e-03_wp
            htt(24,14,1) = -0.12800e-02_wp
            htt(24,15,1) = -0.20000e-03_wp
            htt(24,16,1) = -0.11600e-02_wp
            htt(24,17,1) = 0.21600e-02_wp
            htt(24,18,1) = -0.42700e-02_wp
            htt(24,19,1) = 0.31200e-02_wp
            htt(24,20,1) = -0.33100e-02_wp
            htt(24,21,1) = -0.23000e-02_wp
            htt(24,22,1) = 0.35300e-02_wp
            htt(24,23,1) = 0.41900e-02_wp
            htt(24,24,1) = 0.10010e-01_wp
            ! Hessian elements for anchor point 1
            hbs( 1,1,1)  = -0.26190e-01_wp
            hbs( 1,2,1)  = -0.29060e-01_wp
            hbs( 1,3,1)  = -0.39330e-01_wp
            hbs( 1,4,1)  = -0.78300e-02_wp
            hbs( 1,5,1)  = 0.21800e-02_wp
            hbs( 1,6,1)  = 0.40640e-01_wp
            hbs( 1,7,1)  = -0.41000e-03_wp
            hbs( 1,8,1)  = 0.40740e-01_wp
            hbs( 1,9,1)  = -0.20000e-04_wp
            hbs( 1,10,1) = -0.87600e-02_wp
            hbs( 1,11,1) = -0.35000e-03_wp
            hbs( 1,12,1) = 0.45000e-02_wp
            hbs( 2,1,1)  = 0.36740e-01_wp
            hbs( 2,2,1)  = -0.10920e-01_wp
            hbs( 2,3,1)  = 0.13650e-01_wp
            hbs( 2,4,1)  = 0.34200e-02_wp
            hbs( 2,5,1)  = 0.62000e-03_wp
            hbs( 2,6,1)  = -0.29580e-01_wp
            hbs( 2,7,1)  = 0.13300e-02_wp
            hbs( 2,8,1)  = -0.12060e-01_wp
            hbs( 2,9,1)  = 0.19000e-03_wp
            hbs( 2,10,1) = 0.40800e-02_wp
            hbs( 2,11,1) = -0.39000e-03_wp
            hbs( 2,12,1) = -0.22000e-03_wp
            hbs( 3,1,1)  = -0.10550e-01_wp
            hbs( 3,2,1)  = 0.39980e-01_wp
            hbs( 3,3,1)  = 0.25680e-01_wp
            hbs( 3,4,1)  = 0.44100e-02_wp
            hbs( 3,5,1)  = -0.28000e-02_wp
            hbs( 3,6,1)  = -0.11060e-01_wp
            hbs( 3,7,1)  = -0.91000e-03_wp
            hbs( 3,8,1)  = -0.28680e-01_wp
            hbs( 3,9,1)  = -0.17000e-03_wp
            hbs( 3,10,1) = 0.46800e-02_wp
            hbs( 3,11,1) = 0.74000e-03_wp
            hbs( 3,12,1) = -0.42800e-02_wp
            hbs( 4,1,1)  = -0.26180e-01_wp
            hbs( 4,2,1)  = -0.40500e-02_wp
            hbs( 4,3,1)  = 0.14760e-01_wp
            hbs( 4,4,1)  = -0.35070e-01_wp
            hbs( 4,5,1)  = -0.56800e-02_wp
            hbs( 4,6,1)  = -0.74200e-02_wp
            hbs( 4,7,1)  = 0.38500e-02_wp
            hbs( 4,8,1)  = 0.33760e-01_wp
            hbs( 4,9,1)  = -0.72000e-03_wp
            hbs( 4,10,1) = 0.51570e-01_wp
            hbs( 4,11,1) = 0.20000e-04_wp
            hbs( 4,12,1) = -0.81000e-03_wp
            hbs( 5,1,1)  = 0.28180e-01_wp
            hbs( 5,2,1)  = 0.24300e-02_wp
            hbs( 5,3,1)  = -0.87400e-02_wp
            hbs( 5,4,1)  = 0.50900e-02_wp
            hbs( 5,5,1)  = 0.11940e-01_wp
            hbs( 5,6,1)  = 0.32000e-02_wp
            hbs( 5,7,1)  = -0.88000e-03_wp
            hbs( 5,8,1)  = -0.14770e-01_wp
            hbs( 5,9,1)  = 0.30000e-04_wp
            hbs( 5,10,1) = -0.27450e-01_wp
            hbs( 5,11,1) = 0.00000e+00_wp
            hbs( 5,12,1) = 0.11400e-02_wp
            hbs( 6,1,1)  = -0.20100e-02_wp
            hbs( 6,2,1)  = 0.16200e-02_wp
            hbs( 6,3,1)  = -0.60300e-02_wp
            hbs( 6,4,1)  = 0.29980e-01_wp
            hbs( 6,5,1)  = -0.62600e-02_wp
            hbs( 6,6,1)  = 0.42200e-02_wp
            hbs( 6,7,1)  = -0.29700e-02_wp
            hbs( 6,8,1)  = -0.18990e-01_wp
            hbs( 6,9,1)  = 0.69000e-03_wp
            hbs( 6,10,1) = -0.24120e-01_wp
            hbs( 6,11,1) = -0.10000e-04_wp
            hbs( 6,12,1) = -0.33000e-03_wp
            hbs( 7,1,1)  = -0.59200e-02_wp
            hbs( 7,2,1)  = 0.29370e-01_wp
            hbs( 7,3,1)  = 0.59700e-02_wp
            hbs( 7,4,1)  = -0.34310e-01_wp
            hbs( 7,5,1)  = 0.35700e-02_wp
            hbs( 7,6,1)  = -0.38110e-01_wp
            hbs( 7,7,1)  = -0.77800e-02_wp
            hbs( 7,8,1)  = -0.29600e-02_wp
            hbs( 7,9,1)  = 0.44000e-02_wp
            hbs( 7,10,1) = 0.44330e-01_wp
            hbs( 7,11,1) = -0.52000e-03_wp
            hbs( 7,12,1) = -0.41000e-03_wp
            hbs( 8,1,1)  = 0.27500e-02_wp
            hbs( 8,2,1)  = -0.16290e-01_wp
            hbs( 8,3,1)  = -0.20500e-02_wp
            hbs( 8,4,1)  = 0.31190e-01_wp
            hbs( 8,5,1)  = -0.23800e-02_wp
            hbs( 8,6,1)  = 0.45900e-02_wp
            hbs( 8,7,1)  = 0.44200e-02_wp
            hbs( 8,8,1)  = 0.12900e-02_wp
            hbs( 8,9,1)  = -0.16900e-02_wp
            hbs( 8,10,1) = -0.19840e-01_wp
            hbs( 8,11,1) = -0.28000e-03_wp
            hbs( 8,12,1) = 0.25000e-03_wp
            hbs( 9,1,1)  = 0.31600e-02_wp
            hbs( 9,2,1)  = -0.13080e-01_wp
            hbs( 9,3,1)  = -0.39200e-02_wp
            hbs( 9,4,1)  = 0.31100e-02_wp
            hbs( 9,5,1)  = -0.11900e-02_wp
            hbs( 9,6,1)  = 0.33520e-01_wp
            hbs( 9,7,1)  = 0.33700e-02_wp
            hbs( 9,8,1)  = 0.16700e-02_wp
            hbs( 9,9,1)  = -0.27000e-02_wp
            hbs( 9,10,1) = -0.24490e-01_wp
            hbs( 9,11,1) = 0.80000e-03_wp
            hbs( 9,12,1) = 0.17000e-03_wp
            hbs(10,1,1)  = 0.35610e-01_wp
            hbs(10,2,1)  = 0.35190e-01_wp
            hbs(10,3,1)  = -0.80000e-03_wp
            hbs(10,4,1)  = -0.15230e-01_wp
            hbs(10,5,1)  = 0.40000e-04_wp
            hbs(10,6,1)  = -0.28120e-01_wp
            hbs(10,7,1)  = 0.45800e-02_wp
            hbs(10,8,1)  = -0.27160e-01_wp
            hbs(10,9,1)  = -0.74500e-02_wp
            hbs(10,10,1) = -0.14900e-01_wp
            hbs(10,11,1) = 0.44000e-02_wp
            hbs(10,12,1) = -0.51000e-03_wp
            hbs(11,1,1)  = -0.19830e-01_wp
            hbs(11,2,1)  = -0.15670e-01_wp
            hbs(11,3,1)  = 0.46000e-03_wp
            hbs(11,4,1)  = 0.87500e-02_wp
            hbs(11,5,1)  = 0.12000e-03_wp
            hbs(11,6,1)  = 0.26580e-01_wp
            hbs(11,7,1)  = -0.28700e-02_wp
            hbs(11,8,1)  = 0.69000e-03_wp
            hbs(11,9,1)  = 0.35900e-02_wp
            hbs(11,10,1) = 0.63600e-02_wp
            hbs(11,11,1) = -0.16100e-02_wp
            hbs(11,12,1) = -0.22000e-03_wp
            hbs(12,1,1)  = -0.15780e-01_wp
            hbs(12,2,1)  = -0.19520e-01_wp
            hbs(12,3,1)  = 0.34000e-03_wp
            hbs(12,4,1)  = 0.64800e-02_wp
            hbs(12,5,1)  = -0.15000e-03_wp
            hbs(12,6,1)  = 0.15400e-02_wp
            hbs(12,7,1)  = -0.17100e-02_wp
            hbs(12,8,1)  = 0.26480e-01_wp
            hbs(12,9,1)  = 0.38600e-02_wp
            hbs(12,10,1) = 0.85400e-02_wp
            hbs(12,11,1) = -0.27900e-02_wp
            hbs(12,12,1) = 0.73000e-03_wp
            hbs(13,1,1)  = 0.28450e-01_wp
            hbs(13,2,1)  = -0.45800e-02_wp
            hbs(13,3,1)  = 0.32200e-02_wp
            hbs(13,4,1)  = 0.42370e-01_wp
            hbs(13,5,1)  = -0.53000e-03_wp
            hbs(13,6,1)  = -0.11000e-02_wp
            hbs(13,7,1)  = -0.47000e-03_wp
            hbs(13,8,1)  = -0.36660e-01_wp
            hbs(13,9,1)  = 0.43000e-02_wp
            hbs(13,10,1) = -0.34790e-01_wp
            hbs(13,11,1) = -0.78500e-02_wp
            hbs(13,12,1) = 0.43000e-02_wp
            hbs(14,1,1)  = -0.12720e-01_wp
            hbs(14,2,1)  = 0.26100e-02_wp
            hbs(14,3,1)  = -0.28000e-02_wp
            hbs(14,4,1)  = -0.23620e-01_wp
            hbs(14,5,1)  = 0.19000e-03_wp
            hbs(14,6,1)  = 0.10100e-02_wp
            hbs(14,7,1)  = 0.75000e-03_wp
            hbs(14,8,1)  = 0.32490e-01_wp
            hbs(14,9,1)  = -0.26900e-02_wp
            hbs(14,10,1) = 0.34600e-02_wp
            hbs(14,11,1) = 0.37900e-02_wp
            hbs(14,12,1) = -0.16300e-02_wp
            hbs(15,1,1)  = -0.15730e-01_wp
            hbs(15,2,1)  = 0.19700e-02_wp
            hbs(15,3,1)  = -0.42000e-03_wp
            hbs(15,4,1)  = -0.18750e-01_wp
            hbs(15,5,1)  = 0.33000e-03_wp
            hbs(15,6,1)  = 0.90000e-04_wp
            hbs(15,7,1)  = -0.28000e-03_wp
            hbs(15,8,1)  = 0.41800e-02_wp
            hbs(15,9,1)  = -0.16100e-02_wp
            hbs(15,10,1) = 0.31330e-01_wp
            hbs(15,11,1) = 0.40600e-02_wp
            hbs(15,12,1) = -0.26700e-02_wp
            hbs(16,1,1)  = -0.57800e-02_wp
            hbs(16,2,1)  = -0.26870e-01_wp
            hbs(16,3,1)  = 0.16180e-01_wp
            hbs(16,4,1)  = 0.50070e-01_wp
            hbs(16,5,1)  = 0.42000e-03_wp
            hbs(16,6,1)  = 0.34120e-01_wp
            hbs(16,7,1)  = 0.25000e-03_wp
            hbs(16,8,1)  = -0.77100e-02_wp
            hbs(16,9,1)  = -0.51000e-03_wp
            hbs(16,10,1) = -0.37450e-01_wp
            hbs(16,11,1) = 0.43000e-02_wp
            hbs(16,12,1) = -0.70600e-02_wp
            hbs(17,1,1)  = 0.20300e-02_wp
            hbs(17,2,1)  = 0.26580e-01_wp
            hbs(17,3,1)  = -0.98100e-02_wp
            hbs(17,4,1)  = -0.26570e-01_wp
            hbs(17,5,1)  = 0.45000e-03_wp
            hbs(17,6,1)  = -0.15140e-01_wp
            hbs(17,7,1)  = -0.15000e-03_wp
            hbs(17,8,1)  = 0.35400e-02_wp
            hbs(17,9,1)  = -0.25000e-03_wp
            hbs(17,10,1) = 0.60500e-02_wp
            hbs(17,11,1) = -0.15900e-02_wp
            hbs(17,12,1) = 0.61300e-02_wp
            hbs(18,1,1)  = 0.37500e-02_wp
            hbs(18,2,1)  = 0.29000e-03_wp
            hbs(18,3,1)  = -0.63700e-02_wp
            hbs(18,4,1)  = -0.23500e-01_wp
            hbs(18,5,1)  = -0.87000e-03_wp
            hbs(18,6,1)  = -0.18980e-01_wp
            hbs(18,7,1)  = -0.10000e-03_wp
            hbs(18,8,1)  = 0.41700e-02_wp
            hbs(18,9,1)  = 0.76000e-03_wp
            hbs(18,10,1) = 0.31400e-01_wp
            hbs(18,11,1) = -0.27100e-02_wp
            hbs(18,12,1) = 0.93000e-03_wp
            ! Hessian elements for anchor point 2
            hss( 1,1,2)  = 0.30793e+00_wp
            hss( 2,1,2)  = 0.25300e-02_wp
            hss( 2,2,2)  = 0.31227e+00_wp
            hss( 3,1,2)  = 0.72500e-01_wp
            hss( 3,2,2)  = 0.70950e-01_wp
            hss( 3,3,2)  = 0.63382e+00_wp
            hss( 4,1,2)  = 0.56300e-02_wp
            hss( 4,2,2)  = 0.20240e-01_wp
            hss( 4,3,2)  = -0.11340e-01_wp
            hss( 4,4,2)  = 0.43881e+00_wp
            hss( 5,1,2)  = 0.31200e-02_wp
            hss( 5,2,2)  = -0.58000e-03_wp
            hss( 5,3,2)  = -0.14100e-02_wp
            hss( 5,4,2)  = 0.19700e-02_wp
            hss( 5,5,2)  = 0.36903e+00_wp
            hss( 6,1,2)  = 0.93200e-02_wp
            hss( 6,2,2)  = 0.48920e-01_wp
            hss( 6,3,2)  = 0.11310e-01_wp
            hss( 6,4,2)  = 0.43540e-01_wp
            hss( 6,5,2)  = -0.31000e-03_wp
            hss( 6,6,2)  = 0.34629e+00_wp
            hss( 7,1,2)  = -0.79000e-03_wp
            hss( 7,2,2)  = -0.53000e-03_wp
            hss( 7,3,2)  = 0.11700e-02_wp
            hss( 7,4,2)  = 0.55200e-02_wp
            hss( 7,5,2)  = 0.15000e-03_wp
            hss( 7,6,2)  = 0.35500e-02_wp
            hss( 7,7,2)  = 0.36483e+00_wp
            hss( 8,1,2)  = 0.48320e-01_wp
            hss( 8,2,2)  = 0.88800e-02_wp
            hss( 8,3,2)  = 0.10960e-01_wp
            hss( 8,4,2)  = -0.62700e-02_wp
            hss( 8,5,2)  = -0.47000e-03_wp
            hss( 8,6,2)  = 0.30770e-01_wp
            hss( 8,7,2)  = -0.33000e-03_wp
            hss( 8,8,2)  = 0.34600e+00_wp
            hss( 9,1,2)  = -0.27000e-03_wp
            hss( 9,2,2)  = -0.90000e-04_wp
            hss( 9,3,2)  = -0.37000e-03_wp
            hss( 9,4,2)  = -0.49000e-03_wp
            hss( 9,5,2)  = -0.90000e-04_wp
            hss( 9,6,2)  = 0.39700e-02_wp
            hss( 9,7,2)  = 0.32000e-03_wp
            hss( 9,8,2)  = 0.40100e-02_wp
            hss( 9,9,2)  = 0.36524e+00_wp
            hss(10,1,2)  = 0.17870e-01_wp
            hss(10,2,2)  = 0.45600e-02_wp
            hss(10,3,2)  = -0.11880e-01_wp
            hss(10,4,2)  = 0.81470e-01_wp
            hss(10,5,2)  = -0.59000e-03_wp
            hss(10,6,2)  = -0.42200e-02_wp
            hss(10,7,2)  = -0.61000e-03_wp
            hss(10,8,2)  = 0.41780e-01_wp
            hss(10,9,2)  = -0.69000e-03_wp
            hss(10,10,2) = 0.44702e+00_wp
            hss(11,1,2)  = -0.59000e-03_wp
            hss(11,2,2)  = -0.11700e-02_wp
            hss(11,3,2)  = 0.11000e-02_wp
            hss(11,4,2)  = -0.48000e-03_wp
            hss(11,5,2)  = 0.30000e-04_wp
            hss(11,6,2)  = -0.44000e-03_wp
            hss(11,7,2)  = 0.16000e-03_wp
            hss(11,8,2)  = 0.35200e-02_wp
            hss(11,9,2)  = 0.40000e-03_wp
            hss(11,10,2) = 0.50000e-02_wp
            hss(11,11,2) = 0.36451e+00_wp
            hss(12,1,2)  = 0.14000e-03_wp
            hss(12,2,2)  = 0.21900e-02_wp
            hss(12,3,2)  = 0.30000e-04_wp
            hss(12,4,2)  = -0.45000e-03_wp
            hss(12,5,2)  = 0.10000e-04_wp
            hss(12,6,2)  = -0.43000e-03_wp
            hss(12,7,2)  = 0.60000e-04_wp
            hss(12,8,2)  = 0.80000e-04_wp
            hss(12,9,2)  = 0.30000e-04_wp
            hss(12,10,2) = 0.27300e-02_wp
            hss(12,11,2) = 0.33000e-03_wp
            hss(12,12,2) = 0.36751e+00_wp
            ! Hessian elements for anchor point 2
            hbb( 1,1,2)  = 0.81080e-01_wp
            hbb( 2,1,2)  = -0.44810e-01_wp
            hbb( 2,2,2)  = 0.14731e+00_wp
            hbb( 3,1,2)  = -0.36270e-01_wp
            hbb( 3,2,2)  = -0.10250e+00_wp
            hbb( 3,3,2)  = 0.13877e+00_wp
            hbb( 4,1,2)  = -0.36710e-01_wp
            hbb( 4,2,2)  = 0.25500e-01_wp
            hbb( 4,3,2)  = 0.11210e-01_wp
            hbb( 4,4,2)  = 0.75960e-01_wp
            hbb( 5,1,2)  = 0.22750e-01_wp
            hbb( 5,2,2)  = -0.18060e-01_wp
            hbb( 5,3,2)  = -0.46900e-02_wp
            hbb( 5,4,2)  = -0.38310e-01_wp
            hbb( 5,5,2)  = 0.72880e-01_wp
            hbb( 6,1,2)  = 0.13960e-01_wp
            hbb( 6,2,2)  = -0.74300e-02_wp
            hbb( 6,3,2)  = -0.65200e-02_wp
            hbb( 6,4,2)  = -0.37640e-01_wp
            hbb( 6,5,2)  = -0.34570e-01_wp
            hbb( 6,6,2)  = 0.72210e-01_wp
            hbb( 7,1,2)  = -0.88100e-02_wp
            hbb( 7,2,2)  = -0.49400e-02_wp
            hbb( 7,3,2)  = 0.13750e-01_wp
            hbb( 7,4,2)  = -0.28850e-01_wp
            hbb( 7,5,2)  = 0.11260e-01_wp
            hbb( 7,6,2)  = 0.17590e-01_wp
            hbb( 7,7,2)  = 0.74530e-01_wp
            hbb( 8,1,2)  = -0.15100e-02_wp
            hbb( 8,2,2)  = 0.76700e-02_wp
            hbb( 8,3,2)  = -0.61600e-02_wp
            hbb( 8,4,2)  = 0.17480e-01_wp
            hbb( 8,5,2)  = -0.61800e-02_wp
            hbb( 8,6,2)  = -0.11300e-01_wp
            hbb( 8,7,2)  = -0.36990e-01_wp
            hbb( 8,8,2)  = 0.78710e-01_wp
            hbb( 9,1,2)  = 0.10320e-01_wp
            hbb( 9,2,2)  = -0.27300e-02_wp
            hbb( 9,3,2)  = -0.75900e-02_wp
            hbb( 9,4,2)  = 0.11370e-01_wp
            hbb( 9,5,2)  = -0.50800e-02_wp
            hbb( 9,6,2)  = -0.62900e-02_wp
            hbb( 9,7,2)  = -0.37540e-01_wp
            hbb( 9,8,2)  = -0.41720e-01_wp
            hbb( 9,9,2)  = 0.79260e-01_wp
            hbb(10,1,2)  = 0.10400e-01_wp
            hbb(10,2,2)  = -0.60200e-02_wp
            hbb(10,3,2)  = -0.43800e-02_wp
            hbb(10,4,2)  = -0.90700e-02_wp
            hbb(10,5,2)  = 0.10230e-01_wp
            hbb(10,6,2)  = -0.11600e-02_wp
            hbb(10,7,2)  = -0.34090e-01_wp
            hbb(10,8,2)  = 0.13820e-01_wp
            hbb(10,9,2)  = 0.20260e-01_wp
            hbb(10,10,2) = 0.76010e-01_wp
            hbb(11,1,2)  = -0.53400e-02_wp
            hbb(11,2,2)  = 0.37500e-02_wp
            hbb(11,3,2)  = 0.15900e-02_wp
            hbb(11,4,2)  = -0.11000e-02_wp
            hbb(11,5,2)  = -0.35000e-02_wp
            hbb(11,6,2)  = 0.46000e-02_wp
            hbb(11,7,2)  = 0.20040e-01_wp
            hbb(11,8,2)  = -0.74900e-02_wp
            hbb(11,9,2)  = -0.12560e-01_wp
            hbb(11,10,2) = -0.37950e-01_wp
            hbb(11,11,2) = 0.79720e-01_wp
            hbb(12,1,2)  = -0.50600e-02_wp
            hbb(12,2,2)  = 0.22700e-02_wp
            hbb(12,3,2)  = 0.27900e-02_wp
            hbb(12,4,2)  = 0.10170e-01_wp
            hbb(12,5,2)  = -0.67300e-02_wp
            hbb(12,6,2)  = -0.34400e-02_wp
            hbb(12,7,2)  = 0.14040e-01_wp
            hbb(12,8,2)  = -0.63400e-02_wp
            hbb(12,9,2)  = -0.77100e-02_wp
            hbb(12,10,2) = -0.38060e-01_wp
            hbb(12,11,2) = -0.41770e-01_wp
            hbb(12,12,2) = 0.79830e-01_wp
            hbb(13,1,2)  = -0.81300e-02_wp
            hbb(13,2,2)  = 0.15080e-01_wp
            hbb(13,3,2)  = -0.69500e-02_wp
            hbb(13,4,2)  = 0.40600e-02_wp
            hbb(13,5,2)  = -0.19200e-02_wp
            hbb(13,6,2)  = -0.21400e-02_wp
            hbb(13,7,2)  = -0.67600e-02_wp
            hbb(13,8,2)  = 0.90200e-02_wp
            hbb(13,9,2)  = -0.22600e-02_wp
            hbb(13,10,2) = -0.33920e-01_wp
            hbb(13,11,2) = 0.13940e-01_wp
            hbb(13,12,2) = 0.19980e-01_wp
            hbb(13,13,2) = 0.73910e-01_wp
            hbb(14,1,2)  = 0.96700e-02_wp
            hbb(14,2,2)  = -0.84100e-02_wp
            hbb(14,3,2)  = -0.12600e-02_wp
            hbb(14,4,2)  = -0.20700e-02_wp
            hbb(14,5,2)  = 0.84000e-03_wp
            hbb(14,6,2)  = 0.12300e-02_wp
            hbb(14,7,2)  = -0.21800e-02_wp
            hbb(14,8,2)  = -0.30500e-02_wp
            hbb(14,9,2)  = 0.52300e-02_wp
            hbb(14,10,2) = 0.20350e-01_wp
            hbb(14,11,2) = -0.77300e-02_wp
            hbb(14,12,2) = -0.12620e-01_wp
            hbb(14,13,2) = -0.37510e-01_wp
            hbb(14,14,2) = 0.78920e-01_wp
            hbb(15,1,2)  = -0.15300e-02_wp
            hbb(15,2,2)  = -0.66700e-02_wp
            hbb(15,3,2)  = 0.82100e-02_wp
            hbb(15,4,2)  = -0.20000e-02_wp
            hbb(15,5,2)  = 0.10900e-02_wp
            hbb(15,6,2)  = 0.91000e-03_wp
            hbb(15,7,2)  = 0.89400e-02_wp
            hbb(15,8,2)  = -0.59700e-02_wp
            hbb(15,9,2)  = -0.29700e-02_wp
            hbb(15,10,2) = 0.13570e-01_wp
            hbb(15,11,2) = -0.62100e-02_wp
            hbb(15,12,2) = -0.73600e-02_wp
            hbb(15,13,2) = -0.36400e-01_wp
            hbb(15,14,2) = -0.41400e-01_wp
            hbb(15,15,2) = 0.77800e-01_wp
            hbb(16,1,2)  = -0.37830e-01_wp
            hbb(16,2,2)  = 0.15190e-01_wp
            hbb(16,3,2)  = 0.22640e-01_wp
            hbb(16,4,2)  = -0.53900e-02_wp
            hbb(16,5,2)  = -0.40000e-02_wp
            hbb(16,6,2)  = 0.93900e-02_wp
            hbb(16,7,2)  = 0.39700e-02_wp
            hbb(16,8,2)  = -0.18200e-02_wp
            hbb(16,9,2)  = -0.21500e-02_wp
            hbb(16,10,2) = -0.93400e-02_wp
            hbb(16,11,2) = 0.10420e-01_wp
            hbb(16,12,2) = -0.10800e-02_wp
            hbb(16,13,2) = -0.29170e-01_wp
            hbb(16,14,2) = 0.11750e-01_wp
            hbb(16,15,2) = 0.17420e-01_wp
            hbb(16,16,2) = 0.77760e-01_wp
            hbb(17,1,2)  = 0.23450e-01_wp
            hbb(17,2,2)  = -0.61200e-02_wp
            hbb(17,3,2)  = -0.17330e-01_wp
            hbb(17,4,2)  = -0.38800e-02_wp
            hbb(17,5,2)  = 0.65900e-02_wp
            hbb(17,6,2)  = -0.27000e-02_wp
            hbb(17,7,2)  = -0.18400e-02_wp
            hbb(17,8,2)  = 0.10400e-02_wp
            hbb(17,9,2)  = 0.79000e-03_wp
            hbb(17,10,2) = 0.10460e-01_wp
            hbb(17,11,2) = -0.67900e-02_wp
            hbb(17,12,2) = -0.36700e-02_wp
            hbb(17,13,2) = 0.11210e-01_wp
            hbb(17,14,2) = -0.50900e-02_wp
            hbb(17,15,2) = -0.61300e-02_wp
            hbb(17,16,2) = -0.39400e-01_wp
            hbb(17,17,2) = 0.75830e-01_wp
            hbb(18,1,2)  = 0.14380e-01_wp
            hbb(18,2,2)  = -0.90700e-02_wp
            hbb(18,3,2)  = -0.53200e-02_wp
            hbb(18,4,2)  = 0.92800e-02_wp
            hbb(18,5,2)  = -0.25900e-02_wp
            hbb(18,6,2)  = -0.66900e-02_wp
            hbb(18,7,2)  = -0.21300e-02_wp
            hbb(18,8,2)  = 0.77000e-03_wp
            hbb(18,9,2)  = 0.13500e-02_wp
            hbb(18,10,2) = -0.11200e-02_wp
            hbb(18,11,2) = -0.36200e-02_wp
            hbb(18,12,2) = 0.47500e-02_wp
            hbb(18,13,2) = 0.17950e-01_wp
            hbb(18,14,2) = -0.66600e-02_wp
            hbb(18,15,2) = -0.11290e-01_wp
            hbb(18,16,2) = -0.38360e-01_wp
            hbb(18,17,2) = -0.36420e-01_wp
            hbb(18,18,2) = 0.74790e-01_wp
            ! Hessian elements for anchor point 2
            htt( 1,1,2)  = 0.87500e-02_wp
            htt( 2,1,2)  = 0.39300e-02_wp
            htt( 2,2,2)  = 0.86800e-02_wp
            htt( 3,1,2)  = -0.25000e-03_wp
            htt( 3,2,2)  = -0.40300e-02_wp
            htt( 3,3,2)  = 0.69500e-02_wp
            htt( 4,1,2)  = -0.50700e-02_wp
            htt( 4,2,2)  = 0.71000e-03_wp
            htt( 4,3,2)  = 0.31700e-02_wp
            htt( 4,4,2)  = 0.89500e-02_wp
            htt( 5,1,2)  = -0.63800e-02_wp
            htt( 5,2,2)  = -0.48700e-02_wp
            htt( 5,3,2)  = 0.26100e-02_wp
            htt( 5,4,2)  = 0.41100e-02_wp
            htt( 5,5,2)  = 0.87900e-02_wp
            htt( 6,1,2)  = -0.41700e-02_wp
            htt( 6,2,2)  = -0.31900e-02_wp
            htt( 6,3,2)  = 0.27300e-02_wp
            htt( 6,4,2)  = 0.37100e-02_wp
            htt( 6,5,2)  = 0.34900e-02_wp
            htt( 6,6,2)  = 0.76600e-02_wp
            htt( 7,1,2)  = 0.23700e-02_wp
            htt( 7,2,2)  = 0.28700e-02_wp
            htt( 7,3,2)  = -0.43900e-02_wp
            htt( 7,4,2)  = -0.39000e-02_wp
            htt( 7,5,2)  = 0.60000e-04_wp
            htt( 7,6,2)  = -0.32200e-02_wp
            htt( 7,7,2)  = 0.66300e-02_wp
            htt( 8,1,2)  = 0.45800e-02_wp
            htt( 8,2,2)  = 0.45500e-02_wp
            htt( 8,3,2)  = -0.42700e-02_wp
            htt( 8,4,2)  = -0.43000e-02_wp
            htt( 8,5,2)  = -0.52500e-02_wp
            htt( 8,6,2)  = 0.95000e-03_wp
            htt( 8,7,2)  = 0.33500e-02_wp
            htt( 8,8,2)  = 0.95500e-02_wp
            htt( 9,1,2)  = -0.59200e-02_wp
            htt( 9,2,2)  = -0.90000e-03_wp
            htt( 9,3,2)  = -0.18700e-02_wp
            htt( 9,4,2)  = 0.31600e-02_wp
            htt( 9,5,2)  = 0.12000e-02_wp
            htt( 9,6,2)  = 0.22300e-02_wp
            htt( 9,7,2)  = -0.27400e-02_wp
            htt( 9,8,2)  = -0.17100e-02_wp
            htt( 9,9,2)  = 0.79900e-02_wp
            htt(10,1,2)  = -0.37700e-02_wp
            htt(10,2,2)  = 0.64000e-03_wp
            htt(10,3,2)  = -0.10100e-02_wp
            htt(10,4,2)  = 0.33900e-02_wp
            htt(10,5,2)  = 0.20800e-02_wp
            htt(10,6,2)  = 0.14400e-02_wp
            htt(10,7,2)  = -0.59000e-03_wp
            htt(10,8,2)  = -0.12400e-02_wp
            htt(10,9,2)  = 0.28400e-02_wp
            htt(10,10,2) = 0.84900e-02_wp
            htt(11,1,2)  = -0.84000e-03_wp
            htt(11,2,2)  = -0.59000e-02_wp
            htt(11,3,2)  = 0.21200e-02_wp
            htt(11,4,2)  = -0.29400e-02_wp
            htt(11,5,2)  = -0.38000e-03_wp
            htt(11,6,2)  = 0.11900e-02_wp
            htt(11,7,2)  = -0.32600e-02_wp
            htt(11,8,2)  = -0.16800e-02_wp
            htt(11,9,2)  = 0.26900e-02_wp
            htt(11,10,2) = -0.18100e-02_wp
            htt(11,11,2) = 0.80400e-02_wp
            htt(12,1,2)  = 0.13200e-02_wp
            htt(12,2,2)  = -0.43700e-02_wp
            htt(12,3,2)  = 0.29800e-02_wp
            htt(12,4,2)  = -0.27100e-02_wp
            htt(12,5,2)  = 0.50000e-03_wp
            htt(12,6,2)  = 0.41000e-03_wp
            htt(12,7,2)  = -0.11200e-02_wp
            htt(12,8,2)  = -0.12100e-02_wp
            htt(12,9,2)  = -0.24600e-02_wp
            htt(12,10,2) = 0.38400e-02_wp
            htt(12,11,2) = 0.35400e-02_wp
            htt(12,12,2) = 0.98300e-02_wp
            htt(13,1,2)  = 0.53000e-03_wp
            htt(13,2,2)  = -0.13800e-02_wp
            htt(13,3,2)  = 0.17500e-02_wp
            htt(13,4,2)  = -0.16000e-03_wp
            htt(13,5,2)  = 0.18900e-02_wp
            htt(13,6,2)  = 0.51000e-03_wp
            htt(13,7,2)  = 0.70000e-03_wp
            htt(13,8,2)  = -0.68000e-03_wp
            htt(13,9,2)  = -0.54000e-02_wp
            htt(13,10,2) = -0.17000e-03_wp
            htt(13,11,2) = -0.33800e-02_wp
            htt(13,12,2) = 0.18400e-02_wp
            htt(13,13,2) = 0.80400e-02_wp
            htt(14,1,2)  = 0.19700e-02_wp
            htt(14,2,2)  = 0.54000e-03_wp
            htt(14,3,2)  = 0.79000e-03_wp
            htt(14,4,2)  = -0.65000e-03_wp
            htt(14,5,2)  = 0.36000e-03_wp
            htt(14,6,2)  = -0.13400e-02_wp
            htt(14,7,2)  = 0.15100e-02_wp
            htt(14,8,2)  = -0.19000e-03_wp
            htt(14,9,2)  = -0.30400e-02_wp
            htt(14,10,2) = 0.14100e-02_wp
            htt(14,11,2) = -0.15300e-02_wp
            htt(14,12,2) = 0.29300e-02_wp
            htt(14,13,2) = 0.18300e-02_wp
            htt(14,14,2) = 0.71800e-02_wp
            htt(15,1,2)  = -0.16100e-02_wp
            htt(15,2,2)  = -0.29100e-02_wp
            htt(15,3,2)  = 0.90000e-03_wp
            htt(15,4,2)  = -0.39000e-03_wp
            htt(15,5,2)  = 0.10100e-02_wp
            htt(15,6,2)  = 0.13000e-02_wp
            htt(15,7,2)  = -0.14300e-02_wp
            htt(15,8,2)  = -0.11500e-02_wp
            htt(15,9,2)  = -0.28000e-03_wp
            htt(15,10,2) = -0.57800e-02_wp
            htt(15,11,2) = 0.10900e-02_wp
            htt(15,12,2) = -0.44100e-02_wp
            htt(15,13,2) = 0.28500e-02_wp
            htt(15,14,2) = -0.26000e-02_wp
            htt(15,15,2) = 0.83100e-02_wp
            htt(16,1,2)  = -0.16000e-03_wp
            htt(16,2,2)  = -0.99000e-03_wp
            htt(16,3,2)  = -0.50000e-04_wp
            htt(16,4,2)  = -0.88000e-03_wp
            htt(16,5,2)  = -0.52000e-03_wp
            htt(16,6,2)  = -0.55000e-03_wp
            htt(16,7,2)  = -0.62000e-03_wp
            htt(16,8,2)  = -0.66000e-03_wp
            htt(16,9,2)  = 0.20700e-02_wp
            htt(16,10,2) = -0.42000e-02_wp
            htt(16,11,2) = 0.29500e-02_wp
            htt(16,12,2) = -0.33300e-02_wp
            htt(16,13,2) = -0.33700e-02_wp
            htt(16,14,2) = 0.27500e-02_wp
            htt(16,15,2) = 0.28600e-02_wp
            htt(16,16,2) = 0.89800e-02_wp
            htt(17,1,2)  = 0.19500e-02_wp
            htt(17,2,2)  = 0.45000e-03_wp
            htt(17,3,2)  = 0.64000e-03_wp
            htt(17,4,2)  = -0.85000e-03_wp
            htt(17,5,2)  = 0.51000e-03_wp
            htt(17,6,2)  = -0.12300e-02_wp
            htt(17,7,2)  = 0.17800e-02_wp
            htt(17,8,2)  = 0.40000e-04_wp
            htt(17,9,2)  = 0.55000e-03_wp
            htt(17,10,2) = -0.15600e-02_wp
            htt(17,11,2) = 0.21200e-02_wp
            htt(17,12,2) = 0.10000e-04_wp
            htt(17,13,2) = -0.55500e-02_wp
            htt(17,14,2) = 0.56000e-03_wp
            htt(17,15,2) = -0.34500e-02_wp
            htt(17,16,2) = 0.26500e-02_wp
            htt(17,17,2) = 0.80100e-02_wp
            htt(18,1,2)  = 0.68000e-03_wp
            htt(18,2,2)  = 0.11500e-02_wp
            htt(18,3,2)  = -0.14100e-02_wp
            htt(18,4,2)  = -0.93000e-03_wp
            htt(18,5,2)  = -0.14500e-02_wp
            htt(18,6,2)  = -0.25900e-02_wp
            htt(18,7,2)  = 0.58000e-03_wp
            htt(18,8,2)  = -0.57000e-03_wp
            htt(18,9,2)  = 0.17600e-02_wp
            htt(18,10,2) = 0.54000e-03_wp
            htt(18,11,2) = 0.12600e-02_wp
            htt(18,12,2) = 0.40000e-04_wp
            htt(18,13,2) = -0.35300e-02_wp
            htt(18,14,2) = 0.18700e-02_wp
            htt(18,15,2) = -0.23200e-02_wp
            htt(18,16,2) = 0.30800e-02_wp
            htt(18,17,2) = 0.27300e-02_wp
            htt(18,18,2) = 0.81100e-02_wp
            htt(19,1,2)  = 0.50000e-03_wp
            htt(19,2,2)  = -0.14600e-02_wp
            htt(19,3,2)  = 0.16000e-02_wp
            htt(19,4,2)  = -0.37000e-03_wp
            htt(19,5,2)  = 0.20400e-02_wp
            htt(19,6,2)  = 0.62000e-03_wp
            htt(19,7,2)  = 0.97000e-03_wp
            htt(19,8,2)  = -0.45000e-03_wp
            htt(19,9,2)  = -0.18100e-02_wp
            htt(19,10,2) = -0.31500e-02_wp
            htt(19,11,2) = 0.26000e-03_wp
            htt(19,12,2) = -0.10800e-02_wp
            htt(19,13,2) = 0.67000e-03_wp
            htt(19,14,2) = -0.48000e-02_wp
            htt(19,15,2) = 0.20000e-02_wp
            htt(19,16,2) = -0.34700e-02_wp
            htt(19,17,2) = 0.19000e-02_wp
            htt(19,18,2) = -0.26700e-02_wp
            htt(19,19,2) = 0.73700e-02_wp
            htt(20,1,2)  = -0.77000e-03_wp
            htt(20,2,2)  = -0.76000e-03_wp
            htt(20,3,2)  = -0.45000e-03_wp
            htt(20,4,2)  = -0.44000e-03_wp
            htt(20,5,2)  = 0.80000e-04_wp
            htt(20,6,2)  = -0.74000e-03_wp
            htt(20,7,2)  = -0.23000e-03_wp
            htt(20,8,2)  = -0.10600e-02_wp
            htt(20,9,2)  = -0.59000e-03_wp
            htt(20,10,2) = -0.10400e-02_wp
            htt(20,11,2) = -0.60000e-03_wp
            htt(20,12,2) = -0.10500e-02_wp
            htt(20,13,2) = 0.26900e-02_wp
            htt(20,14,2) = -0.34800e-02_wp
            htt(20,15,2) = 0.31300e-02_wp
            htt(20,16,2) = -0.30400e-02_wp
            htt(20,17,2) = -0.33700e-02_wp
            htt(20,18,2) = 0.27100e-02_wp
            htt(20,19,2) = 0.28000e-02_wp
            htt(20,20,2) = 0.88800e-02_wp
            htt(21,1,2)  = 0.11400e-02_wp
            htt(21,2,2)  = 0.27200e-02_wp
            htt(21,3,2)  = -0.27900e-02_wp
            htt(21,4,2)  = -0.12100e-02_wp
            htt(21,5,2)  = -0.59100e-02_wp
            htt(21,6,2)  = -0.84000e-03_wp
            htt(21,7,2)  = -0.20900e-02_wp
            htt(21,8,2)  = 0.29800e-02_wp
            htt(21,9,2)  = 0.14100e-02_wp
            htt(21,10,2) = 0.53000e-03_wp
            htt(21,11,2) = -0.26000e-03_wp
            htt(21,12,2) = -0.11400e-02_wp
            htt(21,13,2) = 0.60000e-03_wp
            htt(21,14,2) = -0.15900e-02_wp
            htt(21,15,2) = 0.14800e-02_wp
            htt(21,16,2) = -0.71000e-03_wp
            htt(21,17,2) = -0.54000e-02_wp
            htt(21,18,2) = -0.23000e-03_wp
            htt(21,19,2) = -0.32100e-02_wp
            htt(21,20,2) = 0.19600e-02_wp
            htt(21,21,2) = 0.79800e-02_wp
            htt(22,1,2)  = -0.11900e-02_wp
            htt(22,2,2)  = 0.94000e-03_wp
            htt(22,3,2)  = -0.29200e-02_wp
            htt(22,4,2)  = -0.79000e-03_wp
            htt(22,5,2)  = -0.30000e-03_wp
            htt(22,6,2)  = -0.52400e-02_wp
            htt(22,7,2)  = 0.13800e-02_wp
            htt(22,8,2)  = -0.35600e-02_wp
            htt(22,9,2)  = 0.33000e-03_wp
            htt(22,10,2) = 0.12100e-02_wp
            htt(22,11,2) = -0.19200e-02_wp
            htt(22,12,2) = -0.10400e-02_wp
            htt(22,13,2) = 0.20500e-02_wp
            htt(22,14,2) = 0.20000e-03_wp
            htt(22,15,2) = 0.11700e-02_wp
            htt(22,16,2) = -0.67000e-03_wp
            htt(22,17,2) = -0.35500e-02_wp
            htt(22,18,2) = 0.98000e-03_wp
            htt(22,19,2) = -0.17100e-02_wp
            htt(22,20,2) = 0.28200e-02_wp
            htt(22,21,2) = 0.26300e-02_wp
            htt(22,22,2) = 0.78500e-02_wp
            htt(23,1,2)  = 0.24200e-02_wp
            htt(23,2,2)  = 0.20100e-02_wp
            htt(23,3,2)  = -0.73000e-03_wp
            htt(23,4,2)  = -0.11400e-02_wp
            htt(23,5,2)  = -0.39300e-02_wp
            htt(23,6,2)  = 0.53000e-03_wp
            htt(23,7,2)  = -0.87000e-03_wp
            htt(23,8,2)  = 0.35900e-02_wp
            htt(23,9,2)  = 0.18000e-03_wp
            htt(23,10,2) = -0.16000e-02_wp
            htt(23,11,2) = 0.61000e-03_wp
            htt(23,12,2) = -0.11700e-02_wp
            htt(23,13,2) = -0.14400e-02_wp
            htt(23,14,2) = -0.29100e-02_wp
            htt(23,15,2) = 0.33000e-03_wp
            htt(23,16,2) = -0.11400e-02_wp
            htt(23,17,2) = -0.80000e-04_wp
            htt(23,18,2) = -0.56500e-02_wp
            htt(23,19,2) = 0.14000e-02_wp
            htt(23,20,2) = -0.41800e-02_wp
            htt(23,21,2) = 0.27800e-02_wp
            htt(23,22,2) = -0.19400e-02_wp
            htt(23,23,2) = 0.84000e-02_wp
            htt(24,1,2)  = 0.90000e-04_wp
            htt(24,2,2)  = 0.23000e-03_wp
            htt(24,3,2)  = -0.86000e-03_wp
            htt(24,4,2)  = -0.71000e-03_wp
            htt(24,5,2)  = 0.16700e-02_wp
            htt(24,6,2)  = -0.38700e-02_wp
            htt(24,7,2)  = 0.25900e-02_wp
            htt(24,8,2)  = -0.29500e-02_wp
            htt(24,9,2)  = -0.90000e-03_wp
            htt(24,10,2) = -0.92000e-03_wp
            htt(24,11,2) = -0.10500e-02_wp
            htt(24,12,2) = -0.10700e-02_wp
            htt(24,13,2) = 0.10000e-04_wp
            htt(24,14,2) = -0.11200e-02_wp
            htt(24,15,2) = 0.30000e-04_wp
            htt(24,16,2) = -0.11000e-02_wp
            htt(24,17,2) = 0.17600e-02_wp
            htt(24,18,2) = -0.44400e-02_wp
            htt(24,19,2) = 0.29000e-02_wp
            htt(24,20,2) = -0.33100e-02_wp
            htt(24,21,2) = -0.25700e-02_wp
            htt(24,22,2) = 0.32800e-02_wp
            htt(24,23,2) = 0.36800e-02_wp
            htt(24,24,2) = 0.95300e-02_wp
            ! Hessian elements for anchor point 2
            hbs( 1,1,2)  = -0.25450e-01_wp
            hbs( 1,2,2)  = -0.24360e-01_wp
            hbs( 1,3,2)  = -0.45340e-01_wp
            hbs( 1,4,2)  = -0.86200e-02_wp
            hbs( 1,5,2)  = 0.45000e-02_wp
            hbs( 1,6,2)  = 0.42960e-01_wp
            hbs( 1,7,2)  = -0.49000e-03_wp
            hbs( 1,8,2)  = 0.43220e-01_wp
            hbs( 1,9,2)  = -0.40000e-04_wp
            hbs( 1,10,2) = -0.84400e-02_wp
            hbs( 1,11,2) = -0.54000e-03_wp
            hbs( 1,12,2) = 0.48600e-02_wp
            hbs( 2,1,2)  = 0.39900e-01_wp
            hbs( 2,2,2)  = -0.84900e-02_wp
            hbs( 2,3,2)  = 0.28370e-01_wp
            hbs( 2,4,2)  = 0.69300e-02_wp
            hbs( 2,5,2)  = -0.38100e-02_wp
            hbs( 2,6,2)  = -0.29340e-01_wp
            hbs( 2,7,2)  = 0.12900e-02_wp
            hbs( 2,8,2)  = -0.15540e-01_wp
            hbs( 2,9,2)  = 0.18000e-03_wp
            hbs( 2,10,2) = 0.28900e-02_wp
            hbs( 2,11,2) = -0.42000e-03_wp
            hbs( 2,12,2) = -0.40000e-03_wp
            hbs( 3,1,2)  = -0.14450e-01_wp
            hbs( 3,2,2)  = 0.32850e-01_wp
            hbs( 3,3,2)  = 0.16980e-01_wp
            hbs( 3,4,2)  = 0.16900e-02_wp
            hbs( 3,5,2)  = -0.70000e-03_wp
            hbs( 3,6,2)  = -0.13620e-01_wp
            hbs( 3,7,2)  = -0.79000e-03_wp
            hbs( 3,8,2)  = -0.27680e-01_wp
            hbs( 3,9,2)  = -0.14000e-03_wp
            hbs( 3,10,2) = 0.55400e-02_wp
            hbs( 3,11,2) = 0.95000e-03_wp
            hbs( 3,12,2) = -0.44600e-02_wp
            hbs( 4,1,2)  = -0.22690e-01_wp
            hbs( 4,2,2)  = -0.61800e-02_wp
            hbs( 4,3,2)  = 0.14530e-01_wp
            hbs( 4,4,2)  = -0.35450e-01_wp
            hbs( 4,5,2)  = -0.71600e-02_wp
            hbs( 4,6,2)  = -0.92200e-02_wp
            hbs( 4,7,2)  = 0.43300e-02_wp
            hbs( 4,8,2)  = 0.36360e-01_wp
            hbs( 4,9,2)  = -0.65000e-03_wp
            hbs( 4,10,2) = 0.49680e-01_wp
            hbs( 4,11,2) = 0.19000e-03_wp
            hbs( 4,12,2) = -0.97000e-03_wp
            hbs( 5,1,2)  = 0.21990e-01_wp
            hbs( 5,2,2)  = 0.37900e-02_wp
            hbs( 5,3,2)  = -0.71900e-02_wp
            hbs( 5,4,2)  = 0.47500e-02_wp
            hbs( 5,5,2)  = 0.95400e-02_wp
            hbs( 5,6,2)  = 0.43800e-02_wp
            hbs( 5,7,2)  = -0.13300e-02_wp
            hbs( 5,8,2)  = -0.15760e-01_wp
            hbs( 5,9,2)  = -0.15000e-03_wp
            hbs( 5,10,2) = -0.26690e-01_wp
            hbs( 5,11,2) = -0.13000e-03_wp
            hbs( 5,12,2) = 0.10900e-02_wp
            hbs( 6,1,2)  = 0.70000e-03_wp
            hbs( 6,2,2)  = 0.24000e-02_wp
            hbs( 6,3,2)  = -0.73400e-02_wp
            hbs( 6,4,2)  = 0.30700e-01_wp
            hbs( 6,5,2)  = -0.23800e-02_wp
            hbs( 6,6,2)  = 0.48500e-02_wp
            hbs( 6,7,2)  = -0.29900e-02_wp
            hbs( 6,8,2)  = -0.20600e-01_wp
            hbs( 6,9,2)  = 0.80000e-03_wp
            hbs( 6,10,2) = -0.22990e-01_wp
            hbs( 6,11,2) = -0.60000e-04_wp
            hbs( 6,12,2) = -0.12000e-03_wp
            hbs( 7,1,2)  = -0.73400e-02_wp
            hbs( 7,2,2)  = 0.28830e-01_wp
            hbs( 7,3,2)  = 0.66600e-02_wp
            hbs( 7,4,2)  = -0.33160e-01_wp
            hbs( 7,5,2)  = 0.40800e-02_wp
            hbs( 7,6,2)  = -0.36980e-01_wp
            hbs( 7,7,2)  = -0.81400e-02_wp
            hbs( 7,8,2)  = -0.48000e-02_wp
            hbs( 7,9,2)  = 0.47000e-02_wp
            hbs( 7,10,2) = 0.43210e-01_wp
            hbs( 7,11,2) = -0.64000e-03_wp
            hbs( 7,12,2) = -0.30000e-03_wp
            hbs( 8,1,2)  = 0.34500e-02_wp
            hbs( 8,2,2)  = -0.15860e-01_wp
            hbs( 8,3,2)  = -0.19800e-02_wp
            hbs( 8,4,2)  = 0.30660e-01_wp
            hbs( 8,5,2)  = -0.26700e-02_wp
            hbs( 8,6,2)  = 0.43400e-02_wp
            hbs( 8,7,2)  = 0.39600e-02_wp
            hbs( 8,8,2)  = 0.18400e-02_wp
            hbs( 8,9,2)  = -0.18900e-02_wp
            hbs( 8,10,2) = -0.19260e-01_wp
            hbs( 8,11,2) = -0.26000e-03_wp
            hbs( 8,12,2) = 0.16000e-03_wp
            hbs( 9,1,2)  = 0.38900e-02_wp
            hbs( 9,2,2)  = -0.12980e-01_wp
            hbs( 9,3,2)  = -0.46900e-02_wp
            hbs( 9,4,2)  = 0.25000e-02_wp
            hbs( 9,5,2)  = -0.14100e-02_wp
            hbs( 9,6,2)  = 0.32640e-01_wp
            hbs( 9,7,2)  = 0.41800e-02_wp
            hbs( 9,8,2)  = 0.29500e-02_wp
            hbs( 9,9,2)  = -0.28100e-02_wp
            hbs( 9,10,2) = -0.23950e-01_wp
            hbs( 9,11,2) = 0.90000e-03_wp
            hbs( 9,12,2) = 0.14000e-03_wp
            hbs(10,1,2)  = 0.32210e-01_wp
            hbs(10,2,2)  = 0.32440e-01_wp
            hbs(10,3,2)  = 0.24000e-02_wp
            hbs(10,4,2)  = -0.14180e-01_wp
            hbs(10,5,2)  = -0.19000e-03_wp
            hbs(10,6,2)  = -0.29000e-01_wp
            hbs(10,7,2)  = 0.46200e-02_wp
            hbs(10,8,2)  = -0.28530e-01_wp
            hbs(10,9,2)  = -0.80800e-02_wp
            hbs(10,10,2) = -0.14370e-01_wp
            hbs(10,11,2) = 0.45400e-02_wp
            hbs(10,12,2) = -0.45000e-03_wp
            hbs(11,1,2)  = -0.18020e-01_wp
            hbs(11,2,2)  = -0.14150e-01_wp
            hbs(11,3,2)  = -0.11000e-02_wp
            hbs(11,4,2)  = 0.83100e-02_wp
            hbs(11,5,2)  = 0.45000e-03_wp
            hbs(11,6,2)  = 0.27070e-01_wp
            hbs(11,7,2)  = -0.28400e-02_wp
            hbs(11,8,2)  = 0.14300e-02_wp
            hbs(11,9,2)  = 0.40300e-02_wp
            hbs(11,10,2) = 0.60300e-02_wp
            hbs(11,11,2) = -0.17200e-02_wp
            hbs(11,12,2) = -0.30000e-03_wp
            hbs(12,1,2)  = -0.14190e-01_wp
            hbs(12,2,2)  = -0.18290e-01_wp
            hbs(12,3,2)  = -0.13000e-02_wp
            hbs(12,4,2)  = 0.58700e-02_wp
            hbs(12,5,2)  = -0.26000e-03_wp
            hbs(12,6,2)  = 0.19300e-02_wp
            hbs(12,7,2)  = -0.17800e-02_wp
            hbs(12,8,2)  = 0.27100e-01_wp
            hbs(12,9,2)  = 0.40500e-02_wp
            hbs(12,10,2) = 0.83400e-02_wp
            hbs(12,11,2) = -0.28100e-02_wp
            hbs(12,12,2) = 0.75000e-03_wp
            hbs(13,1,2)  = 0.28370e-01_wp
            hbs(13,2,2)  = -0.63900e-02_wp
            hbs(13,3,2)  = 0.52100e-02_wp
            hbs(13,4,2)  = 0.42290e-01_wp
            hbs(13,5,2)  = -0.27000e-03_wp
            hbs(13,6,2)  = -0.35800e-02_wp
            hbs(13,7,2)  = -0.64000e-03_wp
            hbs(13,8,2)  = -0.36310e-01_wp
            hbs(13,9,2)  = 0.46600e-02_wp
            hbs(13,10,2) = -0.33260e-01_wp
            hbs(13,11,2) = -0.81200e-02_wp
            hbs(13,12,2) = 0.45000e-02_wp
            hbs(14,1,2)  = -0.12750e-01_wp
            hbs(14,2,2)  = 0.38700e-02_wp
            hbs(14,3,2)  = -0.40000e-02_wp
            hbs(14,4,2)  = -0.23630e-01_wp
            hbs(14,5,2)  = 0.12000e-03_wp
            hbs(14,6,2)  = 0.23100e-02_wp
            hbs(14,7,2)  = 0.89000e-03_wp
            hbs(14,8,2)  = 0.31980e-01_wp
            hbs(14,9,2)  = -0.28200e-02_wp
            hbs(14,10,2) = 0.30000e-02_wp
            hbs(14,11,2) = 0.42200e-02_wp
            hbs(14,12,2) = -0.17600e-02_wp
            hbs(15,1,2)  = -0.15620e-01_wp
            hbs(15,2,2)  = 0.25200e-02_wp
            hbs(15,3,2)  = -0.12100e-02_wp
            hbs(15,4,2)  = -0.18660e-01_wp
            hbs(15,5,2)  = 0.14000e-03_wp
            hbs(15,6,2)  = 0.12700e-02_wp
            hbs(15,7,2)  = -0.25000e-03_wp
            hbs(15,8,2)  = 0.43300e-02_wp
            hbs(15,9,2)  = -0.18400e-02_wp
            hbs(15,10,2) = 0.30260e-01_wp
            hbs(15,11,2) = 0.39000e-02_wp
            hbs(15,12,2) = -0.27400e-02_wp
            hbs(16,1,2)  = -0.51100e-02_wp
            hbs(16,2,2)  = -0.24340e-01_wp
            hbs(16,3,2)  = 0.16540e-01_wp
            hbs(16,4,2)  = 0.49110e-01_wp
            hbs(16,5,2)  = -0.98000e-03_wp
            hbs(16,6,2)  = 0.35820e-01_wp
            hbs(16,7,2)  = 0.33000e-03_wp
            hbs(16,8,2)  = -0.99500e-02_wp
            hbs(16,9,2)  = -0.59000e-03_wp
            hbs(16,10,2) = -0.36820e-01_wp
            hbs(16,11,2) = 0.45700e-02_wp
            hbs(16,12,2) = -0.76400e-02_wp
            hbs(17,1,2)  = 0.19500e-02_wp
            hbs(17,2,2)  = 0.23750e-01_wp
            hbs(17,3,2)  = -0.99800e-02_wp
            hbs(17,4,2)  = -0.26130e-01_wp
            hbs(17,5,2)  = 0.10100e-02_wp
            hbs(17,6,2)  = -0.15770e-01_wp
            hbs(17,7,2)  = -0.17000e-03_wp
            hbs(17,8,2)  = 0.47200e-02_wp
            hbs(17,9,2)  = -0.26000e-03_wp
            hbs(17,10,2) = 0.55900e-02_wp
            hbs(17,11,2) = -0.16700e-02_wp
            hbs(17,12,2) = 0.59100e-02_wp
            hbs(18,1,2)  = 0.31600e-02_wp
            hbs(18,2,2)  = 0.59000e-03_wp
            hbs(18,3,2)  = -0.65700e-02_wp
            hbs(18,4,2)  = -0.22990e-01_wp
            hbs(18,5,2)  = -0.30000e-04_wp
            hbs(18,6,2)  = -0.20050e-01_wp
            hbs(18,7,2)  = -0.16000e-03_wp
            hbs(18,8,2)  = 0.52300e-02_wp
            hbs(18,9,2)  = 0.85000e-03_wp
            hbs(18,10,2) = 0.31230e-01_wp
            hbs(18,11,2) = -0.29000e-02_wp
            hbs(18,12,2) = 0.17300e-02_wp
            ! Hessian elements for anchor point 3
            hss( 1,1,3)  = 0.30450e+00_wp
            hss( 2,1,3)  = -0.38000e-02_wp
            hss( 2,2,3)  = 0.30517e+00_wp
            hss( 3,1,3)  = 0.88850e-01_wp
            hss( 3,2,3)  = 0.88210e-01_wp
            hss( 3,3,3)  = 0.63344e+00_wp
            hss( 4,1,3)  = 0.83900e-02_wp
            hss( 4,2,3)  = 0.19510e-01_wp
            hss( 4,3,3)  = -0.15960e-01_wp
            hss( 4,4,3)  = 0.43559e+00_wp
            hss( 5,1,3)  = 0.15900e-02_wp
            hss( 5,2,3)  = 0.10000e-04_wp
            hss( 5,3,3)  = -0.79000e-03_wp
            hss( 5,4,3)  = 0.28500e-02_wp
            hss( 5,5,3)  = 0.36540e+00_wp
            hss( 6,1,3)  = 0.33900e-02_wp
            hss( 6,2,3)  = 0.49580e-01_wp
            hss( 6,3,3)  = 0.17260e-01_wp
            hss( 6,4,3)  = 0.43530e-01_wp
            hss( 6,5,3)  = -0.70000e-04_wp
            hss( 6,6,3)  = 0.35132e+00_wp
            hss( 7,1,3)  = -0.12200e-02_wp
            hss( 7,2,3)  = -0.41000e-03_wp
            hss( 7,3,3)  = 0.12400e-02_wp
            hss( 7,4,3)  = 0.52300e-02_wp
            hss( 7,5,3)  = 0.38000e-03_wp
            hss( 7,6,3)  = 0.38200e-02_wp
            hss( 7,7,3)  = 0.36408e+00_wp
            hss( 8,1,3)  = 0.49330e-01_wp
            hss( 8,2,3)  = 0.31600e-02_wp
            hss( 8,3,3)  = 0.17100e-01_wp
            hss( 8,4,3)  = -0.51700e-02_wp
            hss( 8,5,3)  = 0.13000e-03_wp
            hss( 8,6,3)  = 0.29020e-01_wp
            hss( 8,7,3)  = -0.29000e-03_wp
            hss( 8,8,3)  = 0.35081e+00_wp
            hss( 9,1,3)  = -0.28000e-03_wp
            hss( 9,2,3)  = -0.28000e-03_wp
            hss( 9,3,3)  = -0.38000e-03_wp
            hss( 9,4,3)  = -0.48000e-03_wp
            hss( 9,5,3)  = 0.60000e-04_wp
            hss( 9,6,3)  = 0.41500e-02_wp
            hss( 9,7,3)  = 0.40000e-03_wp
            hss( 9,8,3)  = 0.41400e-02_wp
            hss( 9,9,3)  = 0.36450e+00_wp
            hss(10,1,3)  = 0.19190e-01_wp
            hss(10,2,3)  = 0.79700e-02_wp
            hss(10,3,3)  = -0.16090e-01_wp
            hss(10,4,3)  = 0.80800e-01_wp
            hss(10,5,3)  = -0.65000e-03_wp
            hss(10,6,3)  = -0.49200e-02_wp
            hss(10,7,3)  = -0.53000e-03_wp
            hss(10,8,3)  = 0.43630e-01_wp
            hss(10,9,3)  = -0.50000e-03_wp
            hss(10,10,3) = 0.43670e+00_wp
            hss(11,1,3)  = -0.38000e-03_wp
            hss(11,2,3)  = -0.12000e-02_wp
            hss(11,3,3)  = 0.12500e-02_wp
            hss(11,4,3)  = -0.56000e-03_wp
            hss(11,5,3)  = 0.60000e-04_wp
            hss(11,6,3)  = -0.36000e-03_wp
            hss(11,7,3)  = 0.19000e-03_wp
            hss(11,8,3)  = 0.38300e-02_wp
            hss(11,9,3)  = 0.40000e-03_wp
            hss(11,10,3) = 0.51900e-02_wp
            hss(11,11,3) = 0.36404e+00_wp
            hss(12,1,3)  = 0.24000e-03_wp
            hss(12,2,3)  = 0.21300e-02_wp
            hss(12,3,3)  = 0.70000e-04_wp
            hss(12,4,3)  = -0.68000e-03_wp
            hss(12,5,3)  = 0.14000e-03_wp
            hss(12,6,3)  = -0.18000e-03_wp
            hss(12,7,3)  = 0.70000e-04_wp
            hss(12,8,3)  = -0.40000e-04_wp
            hss(12,9,3)  = 0.70000e-04_wp
            hss(12,10,3) = 0.28100e-02_wp
            hss(12,11,3) = 0.41000e-03_wp
            hss(12,12,3) = 0.36629e+00_wp
            ! Hessian elements for anchor point 3
            hbb( 1,1,3)  = 0.81840e-01_wp
            hbb( 2,1,3)  = -0.43040e-01_wp
            hbb( 2,2,3)  = 0.14305e+00_wp
            hbb( 3,1,3)  = -0.38800e-01_wp
            hbb( 3,2,3)  = -0.10001e+00_wp
            hbb( 3,3,3)  = 0.13881e+00_wp
            hbb( 4,1,3)  = -0.37990e-01_wp
            hbb( 4,2,3)  = 0.26100e-01_wp
            hbb( 4,3,3)  = 0.11900e-01_wp
            hbb( 4,4,3)  = 0.77580e-01_wp
            hbb( 5,1,3)  = 0.23090e-01_wp
            hbb( 5,2,3)  = -0.18080e-01_wp
            hbb( 5,3,3)  = -0.50100e-02_wp
            hbb( 5,4,3)  = -0.38630e-01_wp
            hbb( 5,5,3)  = 0.74170e-01_wp
            hbb( 6,1,3)  = 0.14900e-01_wp
            hbb( 6,2,3)  = -0.80200e-02_wp
            hbb( 6,3,3)  = -0.68800e-02_wp
            hbb( 6,4,3)  = -0.38950e-01_wp
            hbb( 6,5,3)  = -0.35540e-01_wp
            hbb( 6,6,3)  = 0.74490e-01_wp
            hbb( 7,1,3)  = -0.82000e-02_wp
            hbb( 7,2,3)  = -0.62200e-02_wp
            hbb( 7,3,3)  = 0.14420e-01_wp
            hbb( 7,4,3)  = -0.29360e-01_wp
            hbb( 7,5,3)  = 0.11220e-01_wp
            hbb( 7,6,3)  = 0.18140e-01_wp
            hbb( 7,7,3)  = 0.74660e-01_wp
            hbb( 8,1,3)  = -0.15800e-02_wp
            hbb( 8,2,3)  = 0.80800e-02_wp
            hbb( 8,3,3)  = -0.65000e-02_wp
            hbb( 8,4,3)  = 0.17600e-01_wp
            hbb( 8,5,3)  = -0.61400e-02_wp
            hbb( 8,6,3)  = -0.11460e-01_wp
            hbb( 8,7,3)  = -0.36880e-01_wp
            hbb( 8,8,3)  = 0.78170e-01_wp
            hbb( 9,1,3)  = 0.97800e-02_wp
            hbb( 9,2,3)  = -0.18600e-02_wp
            hbb( 9,3,3)  = -0.79200e-02_wp
            hbb( 9,4,3)  = 0.11760e-01_wp
            hbb( 9,5,3)  = -0.50800e-02_wp
            hbb( 9,6,3)  = -0.66800e-02_wp
            hbb( 9,7,3)  = -0.37780e-01_wp
            hbb( 9,8,3)  = -0.41290e-01_wp
            hbb( 9,9,3)  = 0.79060e-01_wp
            hbb(10,1,3)  = 0.10690e-01_wp
            hbb(10,2,3)  = -0.58300e-02_wp
            hbb(10,3,3)  = -0.48700e-02_wp
            hbb(10,4,3)  = -0.91900e-02_wp
            hbb(10,5,3)  = 0.10240e-01_wp
            hbb(10,6,3)  = -0.10500e-02_wp
            hbb(10,7,3)  = -0.34370e-01_wp
            hbb(10,8,3)  = 0.13850e-01_wp
            hbb(10,9,3)  = 0.20520e-01_wp
            hbb(10,10,3) = 0.76450e-01_wp
            hbb(11,1,3)  = -0.53700e-02_wp
            hbb(11,2,3)  = 0.34700e-02_wp
            hbb(11,3,3)  = 0.19000e-02_wp
            hbb(11,4,3)  = -0.11400e-02_wp
            hbb(11,5,3)  = -0.35400e-02_wp
            hbb(11,6,3)  = 0.46900e-02_wp
            hbb(11,7,3)  = 0.20260e-01_wp
            hbb(11,8,3)  = -0.74900e-02_wp
            hbb(11,9,3)  = -0.12770e-01_wp
            hbb(11,10,3) = -0.38210e-01_wp
            hbb(11,11,3) = 0.79260e-01_wp
            hbb(12,1,3)  = -0.53200e-02_wp
            hbb(12,2,3)  = 0.23500e-02_wp
            hbb(12,3,3)  = 0.29700e-02_wp
            hbb(12,4,3)  = 0.10330e-01_wp
            hbb(12,5,3)  = -0.67000e-02_wp
            hbb(12,6,3)  = -0.36300e-02_wp
            hbb(12,7,3)  = 0.14110e-01_wp
            hbb(12,8,3)  = -0.63600e-02_wp
            hbb(12,9,3)  = -0.77500e-02_wp
            hbb(12,10,3) = -0.38240e-01_wp
            hbb(12,11,3) = -0.41050e-01_wp
            hbb(12,12,3) = 0.79290e-01_wp
            hbb(13,1,3)  = -0.77600e-02_wp
            hbb(13,2,3)  = 0.14880e-01_wp
            hbb(13,3,3)  = -0.71200e-02_wp
            hbb(13,4,3)  = 0.36800e-02_wp
            hbb(13,5,3)  = -0.17700e-02_wp
            hbb(13,6,3)  = -0.19200e-02_wp
            hbb(13,7,3)  = -0.65800e-02_wp
            hbb(13,8,3)  = 0.88500e-02_wp
            hbb(13,9,3)  = -0.22700e-02_wp
            hbb(13,10,3) = -0.34350e-01_wp
            hbb(13,11,3) = 0.14110e-01_wp
            hbb(13,12,3) = 0.20240e-01_wp
            hbb(13,13,3) = 0.74690e-01_wp
            hbb(14,1,3)  = 0.95100e-02_wp
            hbb(14,2,3)  = -0.81200e-02_wp
            hbb(14,3,3)  = -0.13900e-02_wp
            hbb(14,4,3)  = -0.19200e-02_wp
            hbb(14,5,3)  = 0.75000e-03_wp
            hbb(14,6,3)  = 0.11700e-02_wp
            hbb(14,7,3)  = -0.22400e-02_wp
            hbb(14,8,3)  = -0.29900e-02_wp
            hbb(14,9,3)  = 0.52300e-02_wp
            hbb(14,10,3) = 0.20540e-01_wp
            hbb(14,11,3) = -0.77700e-02_wp
            hbb(14,12,3) = -0.12770e-01_wp
            hbb(14,13,3) = -0.37840e-01_wp
            hbb(14,14,3) = 0.79020e-01_wp
            hbb(15,1,3)  = -0.17500e-02_wp
            hbb(15,2,3)  = -0.67600e-02_wp
            hbb(15,3,3)  = 0.85100e-02_wp
            hbb(15,4,3)  = -0.17600e-02_wp
            hbb(15,5,3)  = 0.10100e-02_wp
            hbb(15,6,3)  = 0.75000e-03_wp
            hbb(15,7,3)  = 0.88200e-02_wp
            hbb(15,8,3)  = -0.58500e-02_wp
            hbb(15,9,3)  = -0.29600e-02_wp
            hbb(15,10,3) = 0.13810e-01_wp
            hbb(15,11,3) = -0.63400e-02_wp
            hbb(15,12,3) = -0.74700e-02_wp
            hbb(15,13,3) = -0.36860e-01_wp
            hbb(15,14,3) = -0.41180e-01_wp
            hbb(15,15,3) = 0.78040e-01_wp
            hbb(16,1,3)  = -0.38580e-01_wp
            hbb(16,2,3)  = 0.14110e-01_wp
            hbb(16,3,3)  = 0.24470e-01_wp
            hbb(16,4,3)  = -0.47200e-02_wp
            hbb(16,5,3)  = -0.41600e-02_wp
            hbb(16,6,3)  = 0.88700e-02_wp
            hbb(16,7,3)  = 0.38500e-02_wp
            hbb(16,8,3)  = -0.18400e-02_wp
            hbb(16,9,3)  = -0.20100e-02_wp
            hbb(16,10,3) = -0.92300e-02_wp
            hbb(16,11,3) = 0.10360e-01_wp
            hbb(16,12,3) = -0.11200e-02_wp
            hbb(16,13,3) = -0.29700e-01_wp
            hbb(16,14,3) = 0.11950e-01_wp
            hbb(16,15,3) = 0.17750e-01_wp
            hbb(16,16,3) = 0.78380e-01_wp
            hbb(17,1,3)  = 0.23970e-01_wp
            hbb(17,2,3)  = -0.47800e-02_wp
            hbb(17,3,3)  = -0.19190e-01_wp
            hbb(17,4,3)  = -0.44800e-02_wp
            hbb(17,5,3)  = 0.69900e-02_wp
            hbb(17,6,3)  = -0.25200e-02_wp
            hbb(17,7,3)  = -0.18900e-02_wp
            hbb(17,8,3)  = 0.10700e-02_wp
            hbb(17,9,3)  = 0.82000e-03_wp
            hbb(17,10,3) = 0.10360e-01_wp
            hbb(17,11,3) = -0.67400e-02_wp
            hbb(17,12,3) = -0.36100e-02_wp
            hbb(17,13,3) = 0.11220e-01_wp
            hbb(17,14,3) = -0.51000e-02_wp
            hbb(17,15,3) = -0.61200e-02_wp
            hbb(17,16,3) = -0.39180e-01_wp
            hbb(17,17,3) = 0.74630e-01_wp
            hbb(18,1,3)  = 0.14610e-01_wp
            hbb(18,2,3)  = -0.93300e-02_wp
            hbb(18,3,3)  = -0.52800e-02_wp
            hbb(18,4,3)  = 0.91900e-02_wp
            hbb(18,5,3)  = -0.28400e-02_wp
            hbb(18,6,3)  = -0.63600e-02_wp
            hbb(18,7,3)  = -0.19600e-02_wp
            hbb(18,8,3)  = 0.77000e-03_wp
            hbb(18,9,3)  = 0.11900e-02_wp
            hbb(18,10,3) = -0.11200e-02_wp
            hbb(18,11,3) = -0.36100e-02_wp
            hbb(18,12,3) = 0.47400e-02_wp
            hbb(18,13,3) = 0.18470e-01_wp
            hbb(18,14,3) = -0.68500e-02_wp
            hbb(18,15,3) = -0.11630e-01_wp
            hbb(18,16,3) = -0.39190e-01_wp
            hbb(18,17,3) = -0.35450e-01_wp
            hbb(18,18,3) = 0.74640e-01_wp
            ! Hessian elements for anchor point 3
            htt( 1,1,3)  = 0.86000e-02_wp
            htt( 2,1,3)  = 0.34800e-02_wp
            htt( 2,2,3)  = 0.76600e-02_wp
            htt( 3,1,3)  = -0.57000e-03_wp
            htt( 3,2,3)  = -0.37600e-02_wp
            htt( 3,3,3)  = 0.75500e-02_wp
            htt( 4,1,3)  = -0.56900e-02_wp
            htt( 4,2,3)  = 0.42000e-03_wp
            htt( 4,3,3)  = 0.43600e-02_wp
            htt( 4,4,3)  = 0.10470e-01_wp
            htt( 5,1,3)  = -0.61300e-02_wp
            htt( 5,2,3)  = -0.40700e-02_wp
            htt( 5,3,3)  = 0.29900e-02_wp
            htt( 5,4,3)  = 0.50400e-02_wp
            htt( 5,5,3)  = 0.84500e-02_wp
            htt( 6,1,3)  = -0.42600e-02_wp
            htt( 6,2,3)  = -0.28700e-02_wp
            htt( 6,3,3)  = 0.30200e-02_wp
            htt( 6,4,3)  = 0.44100e-02_wp
            htt( 6,5,3)  = 0.35900e-02_wp
            htt( 6,6,3)  = 0.78600e-02_wp
            htt( 7,1,3)  = 0.30000e-02_wp
            htt( 7,2,3)  = 0.31400e-02_wp
            htt( 7,3,3)  = -0.51000e-02_wp
            htt( 7,4,3)  = -0.49600e-02_wp
            htt( 7,5,3)  = -0.62000e-03_wp
            htt( 7,6,3)  = -0.36600e-02_wp
            htt( 7,7,3)  = 0.74400e-02_wp
            htt( 8,1,3)  = 0.48700e-02_wp
            htt( 8,2,3)  = 0.43400e-02_wp
            htt( 8,3,3)  = -0.50700e-02_wp
            htt( 8,4,3)  = -0.56000e-02_wp
            htt( 8,5,3)  = -0.54900e-02_wp
            htt( 8,6,3)  = 0.61000e-03_wp
            htt( 8,7,3)  = 0.44100e-02_wp
            htt( 8,8,3)  = 0.10510e-01_wp
            htt( 9,1,3)  = -0.58900e-02_wp
            htt( 9,2,3)  = -0.93000e-03_wp
            htt( 9,3,3)  = -0.18100e-02_wp
            htt( 9,4,3)  = 0.31600e-02_wp
            htt( 9,5,3)  = 0.11100e-02_wp
            htt( 9,6,3)  = 0.22400e-02_wp
            htt( 9,7,3)  = -0.29500e-02_wp
            htt( 9,8,3)  = -0.18200e-02_wp
            htt( 9,9,3)  = 0.80600e-02_wp
            htt(10,1,3)  = -0.39400e-02_wp
            htt(10,2,3)  = 0.44000e-03_wp
            htt(10,3,3)  = -0.58000e-03_wp
            htt(10,4,3)  = 0.37900e-02_wp
            htt(10,5,3)  = 0.23700e-02_wp
            htt(10,6,3)  = 0.17100e-02_wp
            htt(10,7,3)  = -0.97000e-03_wp
            htt(10,8,3)  = -0.16300e-02_wp
            htt(10,9,3)  = 0.28200e-02_wp
            htt(10,10,3) = 0.84500e-02_wp
            htt(11,1,3)  = -0.51000e-03_wp
            htt(11,2,3)  = -0.53200e-02_wp
            htt(11,3,3)  = 0.15400e-02_wp
            htt(11,4,3)  = -0.32600e-02_wp
            htt(11,5,3)  = -0.10500e-02_wp
            htt(11,6,3)  = 0.78000e-03_wp
            htt(11,7,3)  = -0.31000e-02_wp
            htt(11,8,3)  = -0.12700e-02_wp
            htt(11,9,3)  = 0.28400e-02_wp
            htt(11,10,3) = -0.17800e-02_wp
            htt(11,11,3) = 0.78900e-02_wp
            htt(12,1,3)  = 0.14400e-02_wp
            htt(12,2,3)  = -0.39600e-02_wp
            htt(12,3,3)  = 0.27700e-02_wp
            htt(12,4,3)  = -0.26300e-02_wp
            htt(12,5,3)  = 0.21000e-03_wp
            htt(12,6,3)  = 0.25000e-03_wp
            htt(12,7,3)  = -0.11200e-02_wp
            htt(12,8,3)  = -0.10700e-02_wp
            htt(12,9,3)  = -0.24000e-02_wp
            htt(12,10,3) = 0.38400e-02_wp
            htt(12,11,3) = 0.32800e-02_wp
            htt(12,12,3) = 0.95300e-02_wp
            htt(13,1,3)  = 0.46000e-03_wp
            htt(13,2,3)  = -0.12300e-02_wp
            htt(13,3,3)  = 0.19200e-02_wp
            htt(13,4,3)  = 0.23000e-03_wp
            htt(13,5,3)  = 0.19500e-02_wp
            htt(13,6,3)  = 0.60000e-03_wp
            htt(13,7,3)  = 0.50000e-03_wp
            htt(13,8,3)  = -0.85000e-03_wp
            htt(13,9,3)  = -0.54600e-02_wp
            htt(13,10,3) = -0.30000e-04_wp
            htt(13,11,3) = -0.36800e-02_wp
            htt(13,12,3) = 0.17400e-02_wp
            htt(13,13,3) = 0.82300e-02_wp
            htt(14,1,3)  = 0.18200e-02_wp
            htt(14,2,3)  = 0.48000e-03_wp
            htt(14,3,3)  = 0.11000e-02_wp
            htt(14,4,3)  = -0.24000e-03_wp
            htt(14,5,3)  = 0.60000e-03_wp
            htt(14,6,3)  = -0.11700e-02_wp
            htt(14,7,3)  = 0.13200e-02_wp
            htt(14,8,3)  = -0.45000e-03_wp
            htt(14,9,3)  = -0.31800e-02_wp
            htt(14,10,3) = 0.14500e-02_wp
            htt(14,11,3) = -0.17700e-02_wp
            htt(14,12,3) = 0.28600e-02_wp
            htt(14,13,3) = 0.21600e-02_wp
            htt(14,14,3) = 0.73500e-02_wp
            htt(15,1,3)  = -0.14800e-02_wp
            htt(15,2,3)  = -0.25800e-02_wp
            htt(15,3,3)  = 0.70000e-03_wp
            htt(15,4,3)  = -0.40000e-03_wp
            htt(15,5,3)  = 0.70000e-03_wp
            htt(15,6,3)  = 0.11300e-02_wp
            htt(15,7,3)  = -0.14600e-02_wp
            htt(15,8,3)  = -0.10400e-02_wp
            htt(15,9,3)  = -0.26000e-03_wp
            htt(15,10,3) = -0.56100e-02_wp
            htt(15,11,3) = 0.90000e-03_wp
            htt(15,12,3) = -0.44500e-02_wp
            htt(15,13,3) = 0.28500e-02_wp
            htt(15,14,3) = -0.24300e-02_wp
            htt(15,15,3) = 0.81500e-02_wp
            htt(16,1,3)  = -0.11000e-03_wp
            htt(16,2,3)  = -0.87000e-03_wp
            htt(16,3,3)  = -0.11000e-03_wp
            htt(16,4,3)  = -0.87000e-03_wp
            htt(16,5,3)  = -0.64000e-03_wp
            htt(16,6,3)  = -0.64000e-03_wp
            htt(16,7,3)  = -0.64000e-03_wp
            htt(16,8,3)  = -0.64000e-03_wp
            htt(16,9,3)  = 0.20100e-02_wp
            htt(16,10,3) = -0.41300e-02_wp
            htt(16,11,3) = 0.28100e-02_wp
            htt(16,12,3) = -0.33300e-02_wp
            htt(16,13,3) = -0.32200e-02_wp
            htt(16,14,3) = 0.27600e-02_wp
            htt(16,15,3) = 0.28700e-02_wp
            htt(16,16,3) = 0.88500e-02_wp
            htt(17,1,3)  = 0.20300e-02_wp
            htt(17,2,3)  = 0.64000e-03_wp
            htt(17,3,3)  = 0.51000e-03_wp
            htt(17,4,3)  = -0.88000e-03_wp
            htt(17,5,3)  = 0.38000e-03_wp
            htt(17,6,3)  = -0.12800e-02_wp
            htt(17,7,3)  = 0.19000e-02_wp
            htt(17,8,3)  = 0.24000e-03_wp
            htt(17,9,3)  = 0.64000e-03_wp
            htt(17,10,3) = -0.15500e-02_wp
            htt(17,11,3) = 0.21000e-02_wp
            htt(17,12,3) = -0.80000e-04_wp
            htt(17,13,3) = -0.58000e-02_wp
            htt(17,14,3) = 0.29000e-03_wp
            htt(17,15,3) = -0.36300e-02_wp
            htt(17,16,3) = 0.24600e-02_wp
            htt(17,17,3) = 0.82200e-02_wp
            htt(18,1,3)  = 0.77000e-03_wp
            htt(18,2,3)  = 0.11200e-02_wp
            htt(18,3,3)  = -0.14100e-02_wp
            htt(18,4,3)  = -0.10600e-02_wp
            htt(18,5,3)  = -0.15200e-02_wp
            htt(18,6,3)  = -0.26700e-02_wp
            htt(18,7,3)  = 0.66000e-03_wp
            htt(18,8,3)  = -0.50000e-03_wp
            htt(18,9,3)  = 0.17500e-02_wp
            htt(18,10,3) = 0.41000e-03_wp
            htt(18,11,3) = 0.13800e-02_wp
            htt(18,12,3) = 0.40000e-04_wp
            htt(18,13,3) = -0.36300e-02_wp
            htt(18,14,3) = 0.16800e-02_wp
            htt(18,15,3) = -0.23000e-02_wp
            htt(18,16,3) = 0.30100e-02_wp
            htt(18,17,3) = 0.28700e-02_wp
            htt(18,18,3) = 0.81700e-02_wp
            htt(19,1,3)  = 0.67000e-03_wp
            htt(19,2,3)  = -0.10700e-02_wp
            htt(19,3,3)  = 0.13300e-02_wp
            htt(19,4,3)  = -0.41000e-03_wp
            htt(19,5,3)  = 0.17300e-02_wp
            htt(19,6,3)  = 0.49000e-03_wp
            htt(19,7,3)  = 0.10800e-02_wp
            htt(19,8,3)  = -0.16000e-03_wp
            htt(19,9,3)  = -0.16400e-02_wp
            htt(19,10,3) = -0.30300e-02_wp
            htt(19,11,3) = 0.19000e-03_wp
            htt(19,12,3) = -0.12000e-02_wp
            htt(19,13,3) = 0.27000e-03_wp
            htt(19,14,3) = -0.49100e-02_wp
            htt(19,15,3) = 0.16600e-02_wp
            htt(19,16,3) = -0.35200e-02_wp
            htt(19,17,3) = 0.21400e-02_wp
            htt(19,18,3) = -0.24400e-02_wp
            htt(19,19,3) = 0.73200e-02_wp
            htt(20,1,3)  = -0.59000e-03_wp
            htt(20,2,3)  = -0.59000e-03_wp
            htt(20,3,3)  = -0.60000e-03_wp
            htt(20,4,3)  = -0.59000e-03_wp
            htt(20,5,3)  = -0.17000e-03_wp
            htt(20,6,3)  = -0.90000e-03_wp
            htt(20,7,3)  = -0.16000e-03_wp
            htt(20,8,3)  = -0.90000e-03_wp
            htt(20,9,3)  = -0.53000e-03_wp
            htt(20,10,3) = -0.10700e-02_wp
            htt(20,11,3) = -0.54000e-03_wp
            htt(20,12,3) = -0.10700e-02_wp
            htt(20,13,3) = 0.24500e-02_wp
            htt(20,14,3) = -0.35100e-02_wp
            htt(20,15,3) = 0.29800e-02_wp
            htt(20,16,3) = -0.29800e-02_wp
            htt(20,17,3) = -0.32100e-02_wp
            htt(20,18,3) = 0.28600e-02_wp
            htt(20,19,3) = 0.27400e-02_wp
            htt(20,20,3) = 0.88200e-02_wp
            htt(21,1,3)  = 0.99000e-03_wp
            htt(21,2,3)  = 0.20900e-02_wp
            htt(21,3,3)  = -0.29700e-02_wp
            htt(21,4,3)  = -0.18700e-02_wp
            htt(21,5,3)  = -0.57000e-02_wp
            htt(21,6,3)  = -0.91000e-03_wp
            htt(21,7,3)  = -0.17600e-02_wp
            htt(21,8,3)  = 0.30300e-02_wp
            htt(21,9,3)  = 0.14100e-02_wp
            htt(21,10,3) = 0.29000e-03_wp
            htt(21,11,3) = 0.26000e-03_wp
            htt(21,12,3) = -0.86000e-03_wp
            htt(21,13,3) = 0.68000e-03_wp
            htt(21,14,3) = -0.16300e-02_wp
            htt(21,15,3) = 0.18000e-02_wp
            htt(21,16,3) = -0.51000e-03_wp
            htt(21,17,3) = -0.54200e-02_wp
            htt(21,18,3) = -0.27000e-03_wp
            htt(21,19,3) = -0.31100e-02_wp
            htt(21,20,3) = 0.20400e-02_wp
            htt(21,21,3) = 0.79000e-02_wp
            htt(22,1,3)  = -0.98000e-03_wp
            htt(22,2,3)  = 0.81000e-03_wp
            htt(22,3,3)  = -0.30000e-02_wp
            htt(22,4,3)  = -0.12000e-02_wp
            htt(22,5,3)  = -0.56000e-03_wp
            htt(22,6,3)  = -0.54200e-02_wp
            htt(22,7,3)  = 0.14400e-02_wp
            htt(22,8,3)  = -0.34100e-02_wp
            htt(22,9,3)  = 0.22000e-03_wp
            htt(22,10,3) = 0.98000e-03_wp
            htt(22,11,3) = -0.16700e-02_wp
            htt(22,12,3) = -0.91000e-03_wp
            htt(22,13,3) = 0.21100e-02_wp
            htt(22,14,3) = 0.24000e-03_wp
            htt(22,15,3) = 0.13500e-02_wp
            htt(22,16,3) = -0.51000e-03_wp
            htt(22,17,3) = -0.36600e-02_wp
            htt(22,18,3) = 0.95000e-03_wp
            htt(22,19,3) = -0.18000e-02_wp
            htt(22,20,3) = 0.28200e-02_wp
            htt(22,21,3) = 0.28500e-02_wp
            htt(22,22,3) = 0.79700e-02_wp
            htt(23,1,3)  = 0.22700e-02_wp
            htt(23,2,3)  = 0.16000e-02_wp
            htt(23,3,3)  = -0.10300e-02_wp
            htt(23,4,3)  = -0.16900e-02_wp
            htt(23,5,3)  = -0.37900e-02_wp
            htt(23,6,3)  = 0.50000e-03_wp
            htt(23,7,3)  = -0.51000e-03_wp
            htt(23,8,3)  = 0.37700e-02_wp
            htt(23,9,3)  = 0.30000e-03_wp
            htt(23,10,3) = -0.16900e-02_wp
            htt(23,11,3) = 0.10000e-02_wp
            htt(23,12,3) = -0.99000e-03_wp
            htt(23,13,3) = -0.15100e-02_wp
            htt(23,14,3) = -0.30400e-02_wp
            htt(23,15,3) = 0.46000e-03_wp
            htt(23,16,3) = -0.10700e-02_wp
            htt(23,17,3) = -0.20000e-04_wp
            htt(23,18,3) = -0.56100e-02_wp
            htt(23,19,3) = 0.15100e-02_wp
            htt(23,20,3) = -0.40900e-02_wp
            htt(23,21,3) = 0.27100e-02_wp
            htt(23,22,3) = -0.18100e-02_wp
            htt(23,23,3) = 0.83500e-02_wp
            htt(24,1,3)  = 0.29000e-03_wp
            htt(24,2,3)  = 0.33000e-03_wp
            htt(24,3,3)  = -0.10600e-02_wp
            htt(24,4,3)  = -0.10200e-02_wp
            htt(24,5,3)  = 0.13500e-02_wp
            htt(24,6,3)  = -0.40100e-02_wp
            htt(24,7,3)  = 0.27000e-02_wp
            htt(24,8,3)  = -0.26700e-02_wp
            htt(24,9,3)  = -0.90000e-03_wp
            htt(24,10,3) = -0.10000e-02_wp
            htt(24,11,3) = -0.94000e-03_wp
            htt(24,12,3) = -0.10400e-02_wp
            htt(24,13,3) = -0.90000e-04_wp
            htt(24,14,3) = -0.11700e-02_wp
            htt(24,15,3) = 0.10000e-04_wp
            htt(24,16,3) = -0.10700e-02_wp
            htt(24,17,3) = 0.17400e-02_wp
            htt(24,18,3) = -0.43900e-02_wp
            htt(24,19,3) = 0.28200e-02_wp
            htt(24,20,3) = -0.33100e-02_wp
            htt(24,21,3) = -0.23500e-02_wp
            htt(24,22,3) = 0.33100e-02_wp
            htt(24,23,3) = 0.38300e-02_wp
            htt(24,24,3) = 0.94900e-02_wp
            ! Hessian elements for anchor point 3
            hbs( 1,1,3)  = -0.23260e-01_wp
            hbs( 1,2,3)  = -0.22950e-01_wp
            hbs( 1,3,3)  = -0.47950e-01_wp
            hbs( 1,4,3)  = -0.88000e-02_wp
            hbs( 1,5,3)  = 0.57100e-02_wp
            hbs( 1,6,3)  = 0.44410e-01_wp
            hbs( 1,7,3)  = -0.55000e-03_wp
            hbs( 1,8,3)  = 0.44440e-01_wp
            hbs( 1,9,3)  = -0.20000e-04_wp
            hbs( 1,10,3) = -0.86200e-02_wp
            hbs( 1,11,3) = -0.54000e-03_wp
            hbs( 1,12,3) = 0.53300e-02_wp
            hbs( 2,1,3)  = 0.35140e-01_wp
            hbs( 2,2,3)  = -0.98400e-02_wp
            hbs( 2,3,3)  = 0.27820e-01_wp
            hbs( 2,4,3)  = 0.60300e-02_wp
            hbs( 2,5,3)  = -0.65700e-02_wp
            hbs( 2,6,3)  = -0.28710e-01_wp
            hbs( 2,7,3)  = 0.12300e-02_wp
            hbs( 2,8,3)  = -0.16910e-01_wp
            hbs( 2,9,3)  = 0.50000e-04_wp
            hbs( 2,10,3) = 0.28400e-02_wp
            hbs( 2,11,3) = -0.71000e-03_wp
            hbs( 2,12,3) = -0.49000e-03_wp
            hbs( 3,1,3)  = -0.11880e-01_wp
            hbs( 3,2,3)  = 0.32790e-01_wp
            hbs( 3,3,3)  = 0.20130e-01_wp
            hbs( 3,4,3)  = 0.27700e-02_wp
            hbs( 3,5,3)  = 0.86000e-03_wp
            hbs( 3,6,3)  = -0.15690e-01_wp
            hbs( 3,7,3)  = -0.68000e-03_wp
            hbs( 3,8,3)  = -0.27530e-01_wp
            hbs( 3,9,3)  = -0.30000e-04_wp
            hbs( 3,10,3) = 0.57800e-02_wp
            hbs( 3,11,3) = 0.12500e-02_wp
            hbs( 3,12,3) = -0.48400e-02_wp
            hbs( 4,1,3)  = -0.23570e-01_wp
            hbs( 4,2,3)  = -0.60500e-02_wp
            hbs( 4,3,3)  = 0.14830e-01_wp
            hbs( 4,4,3)  = -0.35570e-01_wp
            hbs( 4,5,3)  = -0.79300e-02_wp
            hbs( 4,6,3)  = -0.10760e-01_wp
            hbs( 4,7,3)  = 0.47300e-02_wp
            hbs( 4,8,3)  = 0.37200e-01_wp
            hbs( 4,9,3)  = -0.64000e-03_wp
            hbs( 4,10,3) = 0.48250e-01_wp
            hbs( 4,11,3) = 0.25000e-03_wp
            hbs( 4,12,3) = -0.12100e-02_wp
            hbs( 5,1,3)  = 0.21940e-01_wp
            hbs( 5,2,3)  = 0.33300e-02_wp
            hbs( 5,3,3)  = -0.76700e-02_wp
            hbs( 5,4,3)  = 0.47000e-02_wp
            hbs( 5,5,3)  = 0.57400e-02_wp
            hbs( 5,6,3)  = 0.50100e-02_wp
            hbs( 5,7,3)  = -0.16700e-02_wp
            hbs( 5,8,3)  = -0.16330e-01_wp
            hbs( 5,9,3)  = -0.26000e-03_wp
            hbs( 5,10,3) = -0.25900e-01_wp
            hbs( 5,11,3) = -0.12000e-03_wp
            hbs( 5,12,3) = 0.12000e-02_wp
            hbs( 6,1,3)  = 0.16300e-02_wp
            hbs( 6,2,3)  = 0.27200e-02_wp
            hbs( 6,3,3)  = -0.71600e-02_wp
            hbs( 6,4,3)  = 0.30860e-01_wp
            hbs( 6,5,3)  = 0.21900e-02_wp
            hbs( 6,6,3)  = 0.57500e-02_wp
            hbs( 6,7,3)  = -0.30600e-02_wp
            hbs( 6,8,3)  = -0.20870e-01_wp
            hbs( 6,9,3)  = 0.90000e-03_wp
            hbs( 6,10,3) = -0.22350e-01_wp
            hbs( 6,11,3) = -0.13000e-03_wp
            hbs( 6,12,3) = 0.10000e-04_wp
            hbs( 7,1,3)  = -0.76900e-02_wp
            hbs( 7,2,3)  = 0.28950e-01_wp
            hbs( 7,3,3)  = 0.62100e-02_wp
            hbs( 7,4,3)  = -0.33090e-01_wp
            hbs( 7,5,3)  = 0.43500e-02_wp
            hbs( 7,6,3)  = -0.35980e-01_wp
            hbs( 7,7,3)  = -0.84400e-02_wp
            hbs( 7,8,3)  = -0.54300e-02_wp
            hbs( 7,9,3)  = 0.48600e-02_wp
            hbs( 7,10,3) = 0.43110e-01_wp
            hbs( 7,11,3) = -0.75000e-03_wp
            hbs( 7,12,3) = -0.33000e-03_wp
            hbs( 8,1,3)  = 0.32000e-02_wp
            hbs( 8,2,3)  = -0.15810e-01_wp
            hbs( 8,3,3)  = -0.15400e-02_wp
            hbs( 8,4,3)  = 0.30210e-01_wp
            hbs( 8,5,3)  = -0.26700e-02_wp
            hbs( 8,6,3)  = 0.42000e-02_wp
            hbs( 8,7,3)  = 0.39400e-02_wp
            hbs( 8,8,3)  = 0.20700e-02_wp
            hbs( 8,9,3)  = -0.19800e-02_wp
            hbs( 8,10,3) = -0.19200e-01_wp
            hbs( 8,11,3) = -0.22000e-03_wp
            hbs( 8,12,3) = 0.19000e-03_wp
            hbs( 9,1,3)  = 0.44900e-02_wp
            hbs( 9,2,3)  = -0.13140e-01_wp
            hbs( 9,3,3)  = -0.46700e-02_wp
            hbs( 9,4,3)  = 0.28800e-02_wp
            hbs( 9,5,3)  = -0.16800e-02_wp
            hbs( 9,6,3)  = 0.31780e-01_wp
            hbs( 9,7,3)  = 0.45000e-02_wp
            hbs( 9,8,3)  = 0.33600e-02_wp
            hbs( 9,9,3)  = -0.28900e-02_wp
            hbs( 9,10,3) = -0.23910e-01_wp
            hbs( 9,11,3) = 0.97000e-03_wp
            hbs( 9,12,3) = 0.14000e-03_wp
            hbs(10,1,3)  = 0.31310e-01_wp
            hbs(10,2,3)  = 0.31460e-01_wp
            hbs(10,3,3)  = 0.56600e-02_wp
            hbs(10,4,3)  = -0.13750e-01_wp
            hbs(10,5,3)  = -0.28000e-03_wp
            hbs(10,6,3)  = -0.29490e-01_wp
            hbs(10,7,3)  = 0.47400e-02_wp
            hbs(10,8,3)  = -0.29500e-01_wp
            hbs(10,9,3)  = -0.84400e-02_wp
            hbs(10,10,3) = -0.13800e-01_wp
            hbs(10,11,3) = 0.47800e-02_wp
            hbs(10,12,3) = -0.37000e-03_wp
            hbs(11,1,3)  = -0.17710e-01_wp
            hbs(11,2,3)  = -0.13600e-01_wp
            hbs(11,3,3)  = -0.28200e-02_wp
            hbs(11,4,3)  = 0.80200e-02_wp
            hbs(11,5,3)  = 0.65000e-03_wp
            hbs(11,6,3)  = 0.27440e-01_wp
            hbs(11,7,3)  = -0.28800e-02_wp
            hbs(11,8,3)  = 0.19800e-02_wp
            hbs(11,9,3)  = 0.42300e-02_wp
            hbs(11,10,3) = 0.57200e-02_wp
            hbs(11,11,3) = -0.18700e-02_wp
            hbs(11,12,3) = -0.35000e-03_wp
            hbs(12,1,3)  = -0.13600e-01_wp
            hbs(12,2,3)  = -0.17860e-01_wp
            hbs(12,3,3)  = -0.28400e-02_wp
            hbs(12,4,3)  = 0.57400e-02_wp
            hbs(12,5,3)  = -0.37000e-03_wp
            hbs(12,6,3)  = 0.20500e-02_wp
            hbs(12,7,3)  = -0.18500e-02_wp
            hbs(12,8,3)  = 0.27520e-01_wp
            hbs(12,9,3)  = 0.42000e-02_wp
            hbs(12,10,3) = 0.80800e-02_wp
            hbs(12,11,3) = -0.29100e-02_wp
            hbs(12,12,3) = 0.72000e-03_wp
            hbs(13,1,3)  = 0.28850e-01_wp
            hbs(13,2,3)  = -0.76600e-02_wp
            hbs(13,3,3)  = 0.57200e-02_wp
            hbs(13,4,3)  = 0.43010e-01_wp
            hbs(13,5,3)  = -0.30000e-03_wp
            hbs(13,6,3)  = -0.53400e-02_wp
            hbs(13,7,3)  = -0.74000e-03_wp
            hbs(13,8,3)  = -0.35790e-01_wp
            hbs(13,9,3)  = 0.48800e-02_wp
            hbs(13,10,3) = -0.33150e-01_wp
            hbs(13,11,3) = -0.84400e-02_wp
            hbs(13,12,3) = 0.46600e-02_wp
            hbs(14,1,3)  = -0.13040e-01_wp
            hbs(14,2,3)  = 0.45600e-02_wp
            hbs(14,3,3)  = -0.44600e-02_wp
            hbs(14,4,3)  = -0.23920e-01_wp
            hbs(14,5,3)  = 0.13000e-03_wp
            hbs(14,6,3)  = 0.32800e-02_wp
            hbs(14,7,3)  = 0.98000e-03_wp
            hbs(14,8,3)  = 0.31630e-01_wp
            hbs(14,9,3)  = -0.29000e-02_wp
            hbs(14,10,3) = 0.29700e-02_wp
            hbs(14,11,3) = 0.43300e-02_wp
            hbs(14,12,3) = -0.18500e-02_wp
            hbs(15,1,3)  = -0.15810e-01_wp
            hbs(15,2,3)  = 0.31000e-02_wp
            hbs(15,3,3)  = -0.12600e-02_wp
            hbs(15,4,3)  = -0.19090e-01_wp
            hbs(15,5,3)  = 0.17000e-03_wp
            hbs(15,6,3)  = 0.20600e-02_wp
            hbs(15,7,3)  = -0.24000e-03_wp
            hbs(15,8,3)  = 0.41600e-02_wp
            hbs(15,9,3)  = -0.19700e-02_wp
            hbs(15,10,3) = 0.30180e-01_wp
            hbs(15,11,3) = 0.41100e-02_wp
            hbs(15,12,3) = -0.28100e-02_wp
            hbs(16,1,3)  = -0.56300e-02_wp
            hbs(16,2,3)  = -0.23750e-01_wp
            hbs(16,3,3)  = 0.15530e-01_wp
            hbs(16,4,3)  = 0.48210e-01_wp
            hbs(16,5,3)  = -0.15500e-02_wp
            hbs(16,6,3)  = 0.37160e-01_wp
            hbs(16,7,3)  = 0.26000e-03_wp
            hbs(16,8,3)  = -0.10920e-01_wp
            hbs(16,9,3)  = -0.64000e-03_wp
            hbs(16,10,3) = -0.35790e-01_wp
            hbs(16,11,3) = 0.47000e-02_wp
            hbs(16,12,3) = -0.80900e-02_wp
            hbs(17,1,3)  = 0.27600e-02_wp
            hbs(17,2,3)  = 0.21560e-01_wp
            hbs(17,3,3)  = -0.91700e-02_wp
            hbs(17,4,3)  = -0.25780e-01_wp
            hbs(17,5,3)  = 0.13700e-02_wp
            hbs(17,6,3)  = -0.16120e-01_wp
            hbs(17,7,3)  = -0.14000e-03_wp
            hbs(17,8,3)  = 0.52400e-02_wp
            hbs(17,9,3)  = -0.26000e-03_wp
            hbs(17,10,3) = 0.48900e-02_wp
            hbs(17,11,3) = -0.16800e-02_wp
            hbs(17,12,3) = 0.62500e-02_wp
            hbs(18,1,3)  = 0.28800e-02_wp
            hbs(18,2,3)  = 0.21900e-02_wp
            hbs(18,3,3)  = -0.63500e-02_wp
            hbs(18,4,3)  = -0.22430e-01_wp
            hbs(18,5,3)  = 0.18000e-03_wp
            hbs(18,6,3)  = -0.21030e-01_wp
            hbs(18,7,3)  = -0.12000e-03_wp
            hbs(18,8,3)  = 0.56800e-02_wp
            hbs(18,9,3)  = 0.91000e-03_wp
            hbs(18,10,3) = 0.30890e-01_wp
            hbs(18,11,3) = -0.30200e-02_wp
            hbs(18,12,3) = 0.18400e-02_wp
            ! Hessian elements for anchor point 4
            hss( 1,1,4)  = 0.30630e+00_wp
            hss( 2,1,4)  = -0.48000e-02_wp
            hss( 2,2,4)  = 0.30630e+00_wp
            hss( 3,1,4)  = 0.89520e-01_wp
            hss( 3,2,4)  = 0.89520e-01_wp
            hss( 3,3,4)  = 0.62529e+00_wp
            hss( 4,1,4)  = 0.88300e-02_wp
            hss( 4,2,4)  = 0.19200e-01_wp
            hss( 4,3,4)  = -0.16480e-01_wp
            hss( 4,4,4)  = 0.43562e+00_wp
            hss( 5,1,4)  = 0.21600e-02_wp
            hss( 5,2,4)  = 0.18000e-03_wp
            hss( 5,3,4)  = 0.10000e-03_wp
            hss( 5,4,4)  = 0.28600e-02_wp
            hss( 5,5,4)  = 0.36547e+00_wp
            hss( 6,1,4)  = 0.27000e-02_wp
            hss( 6,2,4)  = 0.49850e-01_wp
            hss( 6,3,4)  = 0.18220e-01_wp
            hss( 6,4,4)  = 0.43220e-01_wp
            hss( 6,5,4)  = -0.10000e-04_wp
            hss( 6,6,4)  = 0.35171e+00_wp
            hss( 7,1,4)  = -0.12200e-02_wp
            hss( 7,2,4)  = -0.39000e-03_wp
            hss( 7,3,4)  = 0.12600e-02_wp
            hss( 7,4,4)  = 0.51900e-02_wp
            hss( 7,5,4)  = 0.41000e-03_wp
            hss( 7,6,4)  = 0.38600e-02_wp
            hss( 7,7,4)  = 0.36379e+00_wp
            hss( 8,1,4)  = 0.49850e-01_wp
            hss( 8,2,4)  = 0.27000e-02_wp
            hss( 8,3,4)  = 0.18220e-01_wp
            hss( 8,4,4)  = -0.52300e-02_wp
            hss( 8,5,4)  = -0.11000e-03_wp
            hss( 8,6,4)  = 0.28880e-01_wp
            hss( 8,7,4)  = -0.34000e-03_wp
            hss( 8,8,4)  = 0.35171e+00_wp
            hss( 9,1,4)  = -0.29000e-03_wp
            hss( 9,2,4)  = -0.29000e-03_wp
            hss( 9,3,4)  = -0.38000e-03_wp
            hss( 9,4,4)  = -0.48000e-03_wp
            hss( 9,5,4)  = 0.80000e-04_wp
            hss( 9,6,4)  = 0.41700e-02_wp
            hss( 9,7,4)  = 0.40000e-03_wp
            hss( 9,8,4)  = 0.41700e-02_wp
            hss( 9,9,4)  = 0.36487e+00_wp
            hss(10,1,4)  = 0.19200e-01_wp
            hss(10,2,4)  = 0.88300e-02_wp
            hss(10,3,4)  = -0.16480e-01_wp
            hss(10,4,4)  = 0.80730e-01_wp
            hss(10,5,4)  = -0.72000e-03_wp
            hss(10,6,4)  = -0.52300e-02_wp
            hss(10,7,4)  = -0.55000e-03_wp
            hss(10,8,4)  = 0.43220e-01_wp
            hss(10,9,4)  = -0.48000e-03_wp
            hss(10,10,4) = 0.43562e+00_wp
            hss(11,1,4)  = -0.39000e-03_wp
            hss(11,2,4)  = -0.12200e-02_wp
            hss(11,3,4)  = 0.12600e-02_wp
            hss(11,4,4)  = -0.55000e-03_wp
            hss(11,5,4)  = 0.70000e-04_wp
            hss(11,6,4)  = -0.34000e-03_wp
            hss(11,7,4)  = 0.20000e-03_wp
            hss(11,8,4)  = 0.38600e-02_wp
            hss(11,9,4)  = 0.40000e-03_wp
            hss(11,10,4) = 0.51900e-02_wp
            hss(11,11,4) = 0.36379e+00_wp
            hss(12,1,4)  = 0.18000e-03_wp
            hss(12,2,4)  = 0.21600e-02_wp
            hss(12,3,4)  = 0.10000e-03_wp
            hss(12,4,4)  = -0.72000e-03_wp
            hss(12,5,4)  = 0.17000e-03_wp
            hss(12,6,4)  = -0.11000e-03_wp
            hss(12,7,4)  = 0.70000e-04_wp
            hss(12,8,4)  = -0.10000e-04_wp
            hss(12,9,4)  = 0.80000e-04_wp
            hss(12,10,4) = 0.28600e-02_wp
            hss(12,11,4) = 0.41000e-03_wp
            hss(12,12,4) = 0.36547e+00_wp
            ! Hessian elements for anchor point 4
            hbb( 1,1,4)  = 0.81570e-01_wp
            hbb( 2,1,4)  = -0.40780e-01_wp
            hbb( 2,2,4)  = 0.13357e+00_wp
            hbb( 3,1,4)  = -0.40780e-01_wp
            hbb( 3,2,4)  = -0.92780e-01_wp
            hbb( 3,3,4)  = 0.13357e+00_wp
            hbb( 4,1,4)  = -0.38240e-01_wp
            hbb( 4,2,4)  = 0.25780e-01_wp
            hbb( 4,3,4)  = 0.12470e-01_wp
            hbb( 4,4,4)  = 0.78290e-01_wp
            hbb( 5,1,4)  = 0.23810e-01_wp
            hbb( 5,2,4)  = -0.19980e-01_wp
            hbb( 5,3,4)  = -0.38300e-02_wp
            hbb( 5,4,4)  = -0.39070e-01_wp
            hbb( 5,5,4)  = 0.74370e-01_wp
            hbb( 6,1,4)  = 0.14430e-01_wp
            hbb( 6,2,4)  = -0.58000e-02_wp
            hbb( 6,3,4)  = -0.86300e-02_wp
            hbb( 6,4,4)  = -0.39220e-01_wp
            hbb( 6,5,4)  = -0.35300e-01_wp
            hbb( 6,6,4)  = 0.74520e-01_wp
            hbb( 7,1,4)  = -0.78000e-02_wp
            hbb( 7,2,4)  = -0.70900e-02_wp
            hbb( 7,3,4)  = 0.14890e-01_wp
            hbb( 7,4,4)  = -0.29810e-01_wp
            hbb( 7,5,4)  = 0.11250e-01_wp
            hbb( 7,6,4)  = 0.18570e-01_wp
            hbb( 7,7,4)  = 0.74890e-01_wp
            hbb( 8,1,4)  = -0.17500e-02_wp
            hbb( 8,2,4)  = 0.84900e-02_wp
            hbb( 8,3,4)  = -0.67400e-02_wp
            hbb( 8,4,4)  = 0.17840e-01_wp
            hbb( 8,5,4)  = -0.61200e-02_wp
            hbb( 8,6,4)  = -0.11720e-01_wp
            hbb( 8,7,4)  = -0.36950e-01_wp
            hbb( 8,8,4)  = 0.78020e-01_wp
            hbb( 9,1,4)  = 0.95600e-02_wp
            hbb( 9,2,4)  = -0.14000e-02_wp
            hbb( 9,3,4)  = -0.81500e-02_wp
            hbb( 9,4,4)  = 0.11970e-01_wp
            hbb( 9,5,4)  = -0.51200e-02_wp
            hbb( 9,6,4)  = -0.68500e-02_wp
            hbb( 9,7,4)  = -0.37940e-01_wp
            hbb( 9,8,4)  = -0.41060e-01_wp
            hbb( 9,9,4)  = 0.79000e-01_wp
            hbb(10,1,4)  = 0.10520e-01_wp
            hbb(10,2,4)  = -0.52600e-02_wp
            hbb(10,3,4)  = -0.52600e-02_wp
            hbb(10,4,4)  = -0.91700e-02_wp
            hbb(10,5,4)  = 0.10330e-01_wp
            hbb(10,6,4)  = -0.11700e-02_wp
            hbb(10,7,4)  = -0.34390e-01_wp
            hbb(10,8,4)  = 0.13830e-01_wp
            hbb(10,9,4)  = 0.20560e-01_wp
            hbb(10,10,4) = 0.76580e-01_wp
            hbb(11,1,4)  = -0.52600e-02_wp
            hbb(11,2,4)  = 0.31800e-02_wp
            hbb(11,3,4)  = 0.20800e-02_wp
            hbb(11,4,4)  = -0.11400e-02_wp
            hbb(11,5,4)  = -0.36300e-02_wp
            hbb(11,6,4)  = 0.47700e-02_wp
            hbb(11,7,4)  = 0.20190e-01_wp
            hbb(11,8,4)  = -0.74500e-02_wp
            hbb(11,9,4)  = -0.12740e-01_wp
            hbb(11,10,4) = -0.38290e-01_wp
            hbb(11,11,4) = 0.79780e-01_wp
            hbb(12,1,4)  = -0.52600e-02_wp
            hbb(12,2,4)  = 0.20800e-02_wp
            hbb(12,3,4)  = 0.31800e-02_wp
            hbb(12,4,4)  = 0.10300e-01_wp
            hbb(12,5,4)  = -0.67000e-02_wp
            hbb(12,6,4)  = -0.36000e-02_wp
            hbb(12,7,4)  = 0.14200e-01_wp
            hbb(12,8,4)  = -0.63700e-02_wp
            hbb(12,9,4)  = -0.78200e-02_wp
            hbb(12,10,4) = -0.38290e-01_wp
            hbb(12,11,4) = -0.41490e-01_wp
            hbb(12,12,4) = 0.79780e-01_wp
            hbb(13,1,4)  = -0.78000e-02_wp
            hbb(13,2,4)  = 0.14890e-01_wp
            hbb(13,3,4)  = -0.70900e-02_wp
            hbb(13,4,4)  = 0.37200e-02_wp
            hbb(13,5,4)  = -0.18400e-02_wp
            hbb(13,6,4)  = -0.18800e-02_wp
            hbb(13,7,4)  = -0.66100e-02_wp
            hbb(13,8,4)  = 0.88400e-02_wp
            hbb(13,9,4)  = -0.22300e-02_wp
            hbb(13,10,4) = -0.34390e-01_wp
            hbb(13,11,4) = 0.14200e-01_wp
            hbb(13,12,4) = 0.20190e-01_wp
            hbb(13,13,4) = 0.74890e-01_wp
            hbb(14,1,4)  = 0.95600e-02_wp
            hbb(14,2,4)  = -0.81500e-02_wp
            hbb(14,3,4)  = -0.14000e-02_wp
            hbb(14,4,4)  = -0.19200e-02_wp
            hbb(14,5,4)  = 0.79000e-03_wp
            hbb(14,6,4)  = 0.11300e-02_wp
            hbb(14,7,4)  = -0.22300e-02_wp
            hbb(14,8,4)  = -0.29700e-02_wp
            hbb(14,9,4)  = 0.52000e-02_wp
            hbb(14,10,4) = 0.20560e-01_wp
            hbb(14,11,4) = -0.78200e-02_wp
            hbb(14,12,4) = -0.12740e-01_wp
            hbb(14,13,4) = -0.37940e-01_wp
            hbb(14,14,4) = 0.79000e-01_wp
            hbb(15,1,4)  = -0.17500e-02_wp
            hbb(15,2,4)  = -0.67400e-02_wp
            hbb(15,3,4)  = 0.84900e-02_wp
            hbb(15,4,4)  = -0.18000e-02_wp
            hbb(15,5,4)  = 0.10500e-02_wp
            hbb(15,6,4)  = 0.75000e-03_wp
            hbb(15,7,4)  = 0.88400e-02_wp
            hbb(15,8,4)  = -0.58700e-02_wp
            hbb(15,9,4)  = -0.29700e-02_wp
            hbb(15,10,4) = 0.13830e-01_wp
            hbb(15,11,4) = -0.63700e-02_wp
            hbb(15,12,4) = -0.74500e-02_wp
            hbb(15,13,4) = -0.36950e-01_wp
            hbb(15,14,4) = -0.41060e-01_wp
            hbb(15,15,4) = 0.78020e-01_wp
            hbb(16,1,4)  = -0.38240e-01_wp
            hbb(16,2,4)  = 0.12470e-01_wp
            hbb(16,3,4)  = 0.25780e-01_wp
            hbb(16,4,4)  = -0.47900e-02_wp
            hbb(16,5,4)  = -0.44800e-02_wp
            hbb(16,6,4)  = 0.92700e-02_wp
            hbb(16,7,4)  = 0.37200e-02_wp
            hbb(16,8,4)  = -0.18000e-02_wp
            hbb(16,9,4)  = -0.19200e-02_wp
            hbb(16,10,4) = -0.91700e-02_wp
            hbb(16,11,4) = 0.10300e-01_wp
            hbb(16,12,4) = -0.11400e-02_wp
            hbb(16,13,4) = -0.29810e-01_wp
            hbb(16,14,4) = 0.11970e-01_wp
            hbb(16,15,4) = 0.17840e-01_wp
            hbb(16,16,4) = 0.78290e-01_wp
            hbb(17,1,4)  = 0.23810e-01_wp
            hbb(17,2,4)  = -0.38300e-02_wp
            hbb(17,3,4)  = -0.19980e-01_wp
            hbb(17,4,4)  = -0.44800e-02_wp
            hbb(17,5,4)  = 0.71700e-02_wp
            hbb(17,6,4)  = -0.26900e-02_wp
            hbb(17,7,4)  = -0.18400e-02_wp
            hbb(17,8,4)  = 0.10500e-02_wp
            hbb(17,9,4)  = 0.79000e-03_wp
            hbb(17,10,4) = 0.10330e-01_wp
            hbb(17,11,4) = -0.67000e-02_wp
            hbb(17,12,4) = -0.36300e-02_wp
            hbb(17,13,4) = 0.11250e-01_wp
            hbb(17,14,4) = -0.51200e-02_wp
            hbb(17,15,4) = -0.61200e-02_wp
            hbb(17,16,4) = -0.39070e-01_wp
            hbb(17,17,4) = 0.74370e-01_wp
            hbb(18,1,4)  = 0.14430e-01_wp
            hbb(18,2,4)  = -0.86300e-02_wp
            hbb(18,3,4)  = -0.58000e-02_wp
            hbb(18,4,4)  = 0.92700e-02_wp
            hbb(18,5,4)  = -0.26900e-02_wp
            hbb(18,6,4)  = -0.65700e-02_wp
            hbb(18,7,4)  = -0.18800e-02_wp
            hbb(18,8,4)  = 0.75000e-03_wp
            hbb(18,9,4)  = 0.11300e-02_wp
            hbb(18,10,4) = -0.11700e-02_wp
            hbb(18,11,4) = -0.36000e-02_wp
            hbb(18,12,4) = 0.47700e-02_wp
            hbb(18,13,4) = 0.18570e-01_wp
            hbb(18,14,4) = -0.68500e-02_wp
            hbb(18,15,4) = -0.11720e-01_wp
            hbb(18,16,4) = -0.39220e-01_wp
            hbb(18,17,4) = -0.35300e-01_wp
            hbb(18,18,4) = 0.74520e-01_wp
            ! Hessian elements for anchor point 4
            htt( 1,1,4)  = 0.86600e-02_wp
            htt( 2,1,4)  = 0.36400e-02_wp
            htt( 2,2,4)  = 0.78600e-02_wp
            htt( 3,1,4)  = -0.53000e-03_wp
            htt( 3,2,4)  = -0.36100e-02_wp
            htt( 3,3,4)  = 0.74400e-02_wp
            htt( 4,1,4)  = -0.55500e-02_wp
            htt( 4,2,4)  = 0.61000e-03_wp
            htt( 4,3,4)  = 0.43600e-02_wp
            htt( 4,4,4)  = 0.10520e-01_wp
            htt( 5,1,4)  = -0.62400e-02_wp
            htt( 5,2,4)  = -0.42900e-02_wp
            htt( 5,3,4)  = 0.29500e-02_wp
            htt( 5,4,4)  = 0.49100e-02_wp
            htt( 5,5,4)  = 0.86600e-02_wp
            htt( 6,1,4)  = -0.42900e-02_wp
            htt( 6,2,4)  = -0.29000e-02_wp
            htt( 6,3,4)  = 0.29600e-02_wp
            htt( 6,4,4)  = 0.43500e-02_wp
            htt( 6,5,4)  = 0.36400e-02_wp
            htt( 6,6,4)  = 0.78600e-02_wp
            htt( 7,1,4)  = 0.29500e-02_wp
            htt( 7,2,4)  = 0.29600e-02_wp
            htt( 7,3,4)  = -0.50200e-02_wp
            htt( 7,4,4)  = -0.50000e-02_wp
            htt( 7,5,4)  = -0.53000e-03_wp
            htt( 7,6,4)  = -0.36100e-02_wp
            htt( 7,7,4)  = 0.74400e-02_wp
            htt( 8,1,4)  = 0.49100e-02_wp
            htt( 8,2,4)  = 0.43500e-02_wp
            htt( 8,3,4)  = -0.50000e-02_wp
            htt( 8,4,4)  = -0.55600e-02_wp
            htt( 8,5,4)  = -0.55500e-02_wp
            htt( 8,6,4)  = 0.61000e-03_wp
            htt( 8,7,4)  = 0.43600e-02_wp
            htt( 8,8,4)  = 0.10520e-01_wp
            htt( 9,1,4)  = -0.58800e-02_wp
            htt( 9,2,4)  = -0.96000e-03_wp
            htt( 9,3,4)  = -0.18400e-02_wp
            htt( 9,4,4)  = 0.30900e-02_wp
            htt( 9,5,4)  = 0.10900e-02_wp
            htt( 9,6,4)  = 0.22400e-02_wp
            htt( 9,7,4)  = -0.29500e-02_wp
            htt( 9,8,4)  = -0.18000e-02_wp
            htt( 9,9,4)  = 0.80500e-02_wp
            htt(10,1,4)  = -0.39100e-02_wp
            htt(10,2,4)  = 0.46000e-03_wp
            htt(10,3,4)  = -0.55000e-03_wp
            htt(10,4,4)  = 0.38200e-02_wp
            htt(10,5,4)  = 0.23300e-02_wp
            htt(10,6,4)  = 0.17200e-02_wp
            htt(10,7,4)  = -0.10300e-02_wp
            htt(10,8,4)  = -0.16300e-02_wp
            htt(10,9,4)  = 0.28100e-02_wp
            htt(10,10,4) = 0.84200e-02_wp
            htt(11,1,4)  = -0.58000e-03_wp
            htt(11,2,4)  = -0.54100e-02_wp
            htt(11,3,4)  = 0.14200e-02_wp
            htt(11,4,4)  = -0.34100e-02_wp
            htt(11,5,4)  = -0.97000e-03_wp
            htt(11,6,4)  = 0.77000e-03_wp
            htt(11,7,4)  = -0.29700e-02_wp
            htt(11,8,4)  = -0.12200e-02_wp
            htt(11,9,4)  = 0.28600e-02_wp
            htt(11,10,4) = -0.18000e-02_wp
            htt(11,11,4) = 0.79500e-02_wp
            htt(12,1,4)  = 0.13900e-02_wp
            htt(12,2,4)  = -0.39900e-02_wp
            htt(12,3,4)  = 0.27100e-02_wp
            htt(12,4,4)  = -0.26800e-02_wp
            htt(12,5,4)  = 0.27000e-03_wp
            htt(12,6,4)  = 0.26000e-03_wp
            htt(12,7,4)  = -0.10500e-02_wp
            htt(12,8,4)  = -0.10500e-02_wp
            htt(12,9,4)  = -0.23900e-02_wp
            htt(12,10,4) = 0.38100e-02_wp
            htt(12,11,4) = 0.32900e-02_wp
            htt(12,12,4) = 0.94800e-02_wp
            htt(13,1,4)  = 0.42000e-03_wp
            htt(13,2,4)  = -0.12700e-02_wp
            htt(13,3,4)  = 0.19200e-02_wp
            htt(13,4,4)  = 0.23000e-03_wp
            htt(13,5,4)  = 0.20100e-02_wp
            htt(13,6,4)  = 0.61000e-03_wp
            htt(13,7,4)  = 0.51000e-03_wp
            htt(13,8,4)  = -0.88000e-03_wp
            htt(13,9,4)  = -0.54500e-02_wp
            htt(13,10,4) = -0.30000e-04_wp
            htt(13,11,4) = -0.36700e-02_wp
            htt(13,12,4) = 0.17500e-02_wp
            htt(13,13,4) = 0.82400e-02_wp
            htt(14,1,4)  = 0.18200e-02_wp
            htt(14,2,4)  = 0.50000e-03_wp
            htt(14,3,4)  = 0.11400e-02_wp
            htt(14,4,4)  = -0.18000e-03_wp
            htt(14,5,4)  = 0.61000e-03_wp
            htt(14,6,4)  = -0.11500e-02_wp
            htt(14,7,4)  = 0.12900e-02_wp
            htt(14,8,4)  = -0.47000e-03_wp
            htt(14,9,4)  = -0.31900e-02_wp
            htt(14,10,4) = 0.14500e-02_wp
            htt(14,11,4) = -0.18000e-02_wp
            htt(14,12,4) = 0.28400e-02_wp
            htt(14,13,4) = 0.21800e-02_wp
            htt(14,14,4) = 0.73800e-02_wp
            htt(15,1,4)  = -0.15300e-02_wp
            htt(15,2,4)  = -0.26700e-02_wp
            htt(15,3,4)  = 0.64000e-03_wp
            htt(15,4,4)  = -0.50000e-03_wp
            htt(15,5,4)  = 0.78000e-03_wp
            htt(15,6,4)  = 0.11300e-02_wp
            htt(15,7,4)  = -0.13900e-02_wp
            htt(15,8,4)  = -0.10500e-02_wp
            htt(15,9,4)  = -0.26000e-03_wp
            htt(15,10,4) = -0.55900e-02_wp
            htt(15,11,4) = 0.94000e-03_wp
            htt(15,12,4) = -0.43900e-02_wp
            htt(15,13,4) = 0.28700e-02_wp
            htt(15,14,4) = -0.24200e-02_wp
            htt(15,15,4) = 0.81600e-02_wp
            htt(16,1,4)  = -0.14000e-03_wp
            htt(16,2,4)  = -0.90000e-03_wp
            htt(16,3,4)  = -0.13000e-03_wp
            htt(16,4,4)  = -0.90000e-03_wp
            htt(16,5,4)  = -0.62000e-03_wp
            htt(16,6,4)  = -0.64000e-03_wp
            htt(16,7,4)  = -0.62000e-03_wp
            htt(16,8,4)  = -0.64000e-03_wp
            htt(16,9,4)  = 0.20000e-02_wp
            htt(16,10,4) = -0.41100e-02_wp
            htt(16,11,4) = 0.28100e-02_wp
            htt(16,12,4) = -0.32900e-02_wp
            htt(16,13,4) = -0.31800e-02_wp
            htt(16,14,4) = 0.27800e-02_wp
            htt(16,15,4) = 0.28700e-02_wp
            htt(16,16,4) = 0.88300e-02_wp
            htt(17,1,4)  = 0.20100e-02_wp
            htt(17,2,4)  = 0.61000e-03_wp
            htt(17,3,4)  = 0.51000e-03_wp
            htt(17,4,4)  = -0.88000e-03_wp
            htt(17,5,4)  = 0.42000e-03_wp
            htt(17,6,4)  = -0.12700e-02_wp
            htt(17,7,4)  = 0.19200e-02_wp
            htt(17,8,4)  = 0.23000e-03_wp
            htt(17,9,4)  = 0.64000e-03_wp
            htt(17,10,4) = -0.15500e-02_wp
            htt(17,11,4) = 0.21100e-02_wp
            htt(17,12,4) = -0.80000e-04_wp
            htt(17,13,4) = -0.58000e-02_wp
            htt(17,14,4) = 0.26000e-03_wp
            htt(17,15,4) = -0.36300e-02_wp
            htt(17,16,4) = 0.24300e-02_wp
            htt(17,17,4) = 0.82400e-02_wp
            htt(18,1,4)  = 0.78000e-03_wp
            htt(18,2,4)  = 0.11300e-02_wp
            htt(18,3,4)  = -0.13900e-02_wp
            htt(18,4,4)  = -0.10500e-02_wp
            htt(18,5,4)  = -0.15300e-02_wp
            htt(18,6,4)  = -0.26700e-02_wp
            htt(18,7,4)  = 0.64000e-03_wp
            htt(18,8,4)  = -0.50000e-03_wp
            htt(18,9,4)  = 0.17400e-02_wp
            htt(18,10,4) = 0.40000e-03_wp
            htt(18,11,4) = 0.13800e-02_wp
            htt(18,12,4) = 0.40000e-04_wp
            htt(18,13,4) = -0.36300e-02_wp
            htt(18,14,4) = 0.16600e-02_wp
            htt(18,15,4) = -0.23000e-02_wp
            htt(18,16,4) = 0.29900e-02_wp
            htt(18,17,4) = 0.28700e-02_wp
            htt(18,18,4) = 0.81600e-02_wp
            htt(19,1,4)  = 0.61000e-03_wp
            htt(19,2,4)  = -0.11500e-02_wp
            htt(19,3,4)  = 0.12900e-02_wp
            htt(19,4,4)  = -0.47000e-03_wp
            htt(19,5,4)  = 0.18200e-02_wp
            htt(19,6,4)  = 0.50000e-03_wp
            htt(19,7,4)  = 0.11400e-02_wp
            htt(19,8,4)  = -0.18000e-03_wp
            htt(19,9,4)  = -0.16200e-02_wp
            htt(19,10,4) = -0.30300e-02_wp
            htt(19,11,4) = 0.24000e-03_wp
            htt(19,12,4) = -0.11700e-02_wp
            htt(19,13,4) = 0.26000e-03_wp
            htt(19,14,4) = -0.49400e-02_wp
            htt(19,15,4) = 0.16600e-02_wp
            htt(19,16,4) = -0.35300e-02_wp
            htt(19,17,4) = 0.21800e-02_wp
            htt(19,18,4) = -0.24200e-02_wp
            htt(19,19,4) = 0.73800e-02_wp
            htt(20,1,4)  = -0.62000e-03_wp
            htt(20,2,4)  = -0.64000e-03_wp
            htt(20,3,4)  = -0.62000e-03_wp
            htt(20,4,4)  = -0.64000e-03_wp
            htt(20,5,4)  = -0.14000e-03_wp
            htt(20,6,4)  = -0.90000e-03_wp
            htt(20,7,4)  = -0.13000e-03_wp
            htt(20,8,4)  = -0.90000e-03_wp
            htt(20,9,4)  = -0.52000e-03_wp
            htt(20,10,4) = -0.10800e-02_wp
            htt(20,11,4) = -0.49000e-03_wp
            htt(20,12,4) = -0.10500e-02_wp
            htt(20,13,4) = 0.24300e-02_wp
            htt(20,14,4) = -0.35300e-02_wp
            htt(20,15,4) = 0.29900e-02_wp
            htt(20,16,4) = -0.29700e-02_wp
            htt(20,17,4) = -0.31800e-02_wp
            htt(20,18,4) = 0.28700e-02_wp
            htt(20,19,4) = 0.27800e-02_wp
            htt(20,20,4) = 0.88300e-02_wp
            htt(21,1,4)  = 0.10900e-02_wp
            htt(21,2,4)  = 0.22400e-02_wp
            htt(21,3,4)  = -0.29500e-02_wp
            htt(21,4,4)  = -0.18000e-02_wp
            htt(21,5,4)  = -0.58800e-02_wp
            htt(21,6,4)  = -0.96000e-03_wp
            htt(21,7,4)  = -0.18400e-02_wp
            htt(21,8,4)  = 0.30900e-02_wp
            htt(21,9,4)  = 0.14200e-02_wp
            htt(21,10,4) = 0.31000e-03_wp
            htt(21,11,4) = 0.21000e-03_wp
            htt(21,12,4) = -0.90000e-03_wp
            htt(21,13,4) = 0.64000e-03_wp
            htt(21,14,4) = -0.16200e-02_wp
            htt(21,15,4) = 0.17400e-02_wp
            htt(21,16,4) = -0.52000e-03_wp
            htt(21,17,4) = -0.54500e-02_wp
            htt(21,18,4) = -0.26000e-03_wp
            htt(21,19,4) = -0.31900e-02_wp
            htt(21,20,4) = 0.20000e-02_wp
            htt(21,21,4) = 0.80500e-02_wp
            htt(22,1,4)  = -0.97000e-03_wp
            htt(22,2,4)  = 0.77000e-03_wp
            htt(22,3,4)  = -0.29700e-02_wp
            htt(22,4,4)  = -0.12200e-02_wp
            htt(22,5,4)  = -0.58000e-03_wp
            htt(22,6,4)  = -0.54100e-02_wp
            htt(22,7,4)  = 0.14200e-02_wp
            htt(22,8,4)  = -0.34100e-02_wp
            htt(22,9,4)  = 0.21000e-03_wp
            htt(22,10,4) = 0.95000e-03_wp
            htt(22,11,4) = -0.16300e-02_wp
            htt(22,12,4) = -0.89000e-03_wp
            htt(22,13,4) = 0.21100e-02_wp
            htt(22,14,4) = 0.24000e-03_wp
            htt(22,15,4) = 0.13800e-02_wp
            htt(22,16,4) = -0.49000e-03_wp
            htt(22,17,4) = -0.36700e-02_wp
            htt(22,18,4) = 0.94000e-03_wp
            htt(22,19,4) = -0.18000e-02_wp
            htt(22,20,4) = 0.28100e-02_wp
            htt(22,21,4) = 0.28600e-02_wp
            htt(22,22,4) = 0.79500e-02_wp
            htt(23,1,4)  = 0.23300e-02_wp
            htt(23,2,4)  = 0.17200e-02_wp
            htt(23,3,4)  = -0.10300e-02_wp
            htt(23,4,4)  = -0.16300e-02_wp
            htt(23,5,4)  = -0.39100e-02_wp
            htt(23,6,4)  = 0.46000e-03_wp
            htt(23,7,4)  = -0.55000e-03_wp
            htt(23,8,4)  = 0.38200e-02_wp
            htt(23,9,4)  = 0.31000e-03_wp
            htt(23,10,4) = -0.16600e-02_wp
            htt(23,11,4) = 0.95000e-03_wp
            htt(23,12,4) = -0.10200e-02_wp
            htt(23,13,4) = -0.15500e-02_wp
            htt(23,14,4) = -0.30300e-02_wp
            htt(23,15,4) = 0.40000e-03_wp
            htt(23,16,4) = -0.10800e-02_wp
            htt(23,17,4) = -0.30000e-04_wp
            htt(23,18,4) = -0.55900e-02_wp
            htt(23,19,4) = 0.14500e-02_wp
            htt(23,20,4) = -0.41100e-02_wp
            htt(23,21,4) = 0.28100e-02_wp
            htt(23,22,4) = -0.18000e-02_wp
            htt(23,23,4) = 0.84200e-02_wp
            htt(24,1,4)  = 0.27000e-03_wp
            htt(24,2,4)  = 0.26000e-03_wp
            htt(24,3,4)  = -0.10500e-02_wp
            htt(24,4,4)  = -0.10500e-02_wp
            htt(24,5,4)  = 0.13900e-02_wp
            htt(24,6,4)  = -0.39900e-02_wp
            htt(24,7,4)  = 0.27100e-02_wp
            htt(24,8,4)  = -0.26800e-02_wp
            htt(24,9,4)  = -0.90000e-03_wp
            htt(24,10,4) = -0.10200e-02_wp
            htt(24,11,4) = -0.89000e-03_wp
            htt(24,12,4) = -0.10100e-02_wp
            htt(24,13,4) = -0.80000e-04_wp
            htt(24,14,4) = -0.11700e-02_wp
            htt(24,15,4) = 0.40000e-04_wp
            htt(24,16,4) = -0.10500e-02_wp
            htt(24,17,4) = 0.17500e-02_wp
            htt(24,18,4) = -0.43900e-02_wp
            htt(24,19,4) = 0.28400e-02_wp
            htt(24,20,4) = -0.32900e-02_wp
            htt(24,21,4) = -0.23900e-02_wp
            htt(24,22,4) = 0.32900e-02_wp
            htt(24,23,4) = 0.38100e-02_wp
            htt(24,24,4) = 0.94800e-02_wp
            ! Hessian elements for anchor point 4
            hbs( 1,1,4)  = -0.22700e-01_wp
            hbs( 1,2,4)  = -0.22700e-01_wp
            hbs( 1,3,4)  = -0.46990e-01_wp
            hbs( 1,4,4)  = -0.88700e-02_wp
            hbs( 1,5,4)  = 0.54300e-02_wp
            hbs( 1,6,4)  = 0.44400e-01_wp
            hbs( 1,7,4)  = -0.55000e-03_wp
            hbs( 1,8,4)  = 0.44400e-01_wp
            hbs( 1,9,4)  = -0.10000e-04_wp
            hbs( 1,10,4) = -0.88700e-02_wp
            hbs( 1,11,4) = -0.55000e-03_wp
            hbs( 1,12,4) = 0.54300e-02_wp
            hbs( 2,1,4)  = 0.32950e-01_wp
            hbs( 2,2,4)  = -0.10250e-01_wp
            hbs( 2,3,4)  = 0.23500e-01_wp
            hbs( 2,4,4)  = 0.60000e-02_wp
            hbs( 2,5,4)  = -0.49200e-02_wp
            hbs( 2,6,4)  = -0.28130e-01_wp
            hbs( 2,7,4)  = 0.12500e-02_wp
            hbs( 2,8,4)  = -0.16270e-01_wp
            hbs( 2,9,4)  = 0.00000e+00_wp
            hbs( 2,10,4) = 0.28700e-02_wp
            hbs( 2,11,4) = -0.70000e-03_wp
            hbs( 2,12,4) = -0.51000e-03_wp
            hbs( 3,1,4)  = -0.10250e-01_wp
            hbs( 3,2,4)  = 0.32950e-01_wp
            hbs( 3,3,4)  = 0.23500e-01_wp
            hbs( 3,4,4)  = 0.28700e-02_wp
            hbs( 3,5,4)  = -0.51000e-03_wp
            hbs( 3,6,4)  = -0.16270e-01_wp
            hbs( 3,7,4)  = -0.70000e-03_wp
            hbs( 3,8,4)  = -0.28130e-01_wp
            hbs( 3,9,4)  = 0.00000e+00_wp
            hbs( 3,10,4) = 0.60000e-02_wp
            hbs( 3,11,4) = 0.12500e-02_wp
            hbs( 3,12,4) = -0.49200e-02_wp
            hbs( 4,1,4)  = -0.24180e-01_wp
            hbs( 4,2,4)  = -0.60700e-02_wp
            hbs( 4,3,4)  = 0.14620e-01_wp
            hbs( 4,4,4)  = -0.35760e-01_wp
            hbs( 4,5,4)  = -0.81300e-02_wp
            hbs( 4,6,4)  = -0.10830e-01_wp
            hbs( 4,7,4)  = 0.47400e-02_wp
            hbs( 4,8,4)  = 0.37380e-01_wp
            hbs( 4,9,4)  = -0.65000e-03_wp
            hbs( 4,10,4) = 0.48130e-01_wp
            hbs( 4,11,4) = 0.25000e-03_wp
            hbs( 4,12,4) = -0.12600e-02_wp
            hbs( 5,1,4)  = 0.21630e-01_wp
            hbs( 5,2,4)  = 0.30300e-02_wp
            hbs( 5,3,4)  = -0.86900e-02_wp
            hbs( 5,4,4)  = 0.48300e-02_wp
            hbs( 5,5,4)  = 0.60700e-02_wp
            hbs( 5,6,4)  = 0.52000e-02_wp
            hbs( 5,7,4)  = -0.17000e-02_wp
            hbs( 5,8,4)  = -0.16220e-01_wp
            hbs( 5,9,4)  = -0.27000e-03_wp
            hbs( 5,10,4) = -0.25720e-01_wp
            hbs( 5,11,4) = -0.13000e-03_wp
            hbs( 5,12,4) = 0.12100e-02_wp
            hbs( 6,1,4)  = 0.25500e-02_wp
            hbs( 6,2,4)  = 0.30300e-02_wp
            hbs( 6,3,4)  = -0.59400e-02_wp
            hbs( 6,4,4)  = 0.30930e-01_wp
            hbs( 6,5,4)  = 0.20600e-02_wp
            hbs( 6,6,4)  = 0.56300e-02_wp
            hbs( 6,7,4)  = -0.30400e-02_wp
            hbs( 6,8,4)  = -0.21160e-01_wp
            hbs( 6,9,4)  = 0.91000e-03_wp
            hbs( 6,10,4) = -0.22410e-01_wp
            hbs( 6,11,4) = -0.12000e-03_wp
            hbs( 6,12,4) = 0.50000e-04_wp
            hbs( 7,1,4)  = -0.78200e-02_wp
            hbs( 7,2,4)  = 0.29170e-01_wp
            hbs( 7,3,4)  = 0.57100e-02_wp
            hbs( 7,4,4)  = -0.32990e-01_wp
            hbs( 7,5,4)  = 0.46200e-02_wp
            hbs( 7,6,4)  = -0.35760e-01_wp
            hbs( 7,7,4)  = -0.84700e-02_wp
            hbs( 7,8,4)  = -0.55800e-02_wp
            hbs( 7,9,4)  = 0.48900e-02_wp
            hbs( 7,10,4) = 0.43210e-01_wp
            hbs( 7,11,4) = -0.75000e-03_wp
            hbs( 7,12,4) = -0.32000e-03_wp
            hbs( 8,1,4)  = 0.31800e-02_wp
            hbs( 8,2,4)  = -0.15950e-01_wp
            hbs( 8,3,4)  = -0.12300e-02_wp
            hbs( 8,4,4)  = 0.30010e-01_wp
            hbs( 8,5,4)  = -0.27700e-02_wp
            hbs( 8,6,4)  = 0.41300e-02_wp
            hbs( 8,7,4)  = 0.40300e-02_wp
            hbs( 8,8,4)  = 0.21600e-02_wp
            hbs( 8,9,4)  = -0.19900e-02_wp
            hbs( 8,10,4) = -0.19160e-01_wp
            hbs( 8,11,4) = -0.23000e-03_wp
            hbs( 8,12,4) = 0.17000e-03_wp
            hbs( 9,1,4)  = 0.46400e-02_wp
            hbs( 9,2,4)  = -0.13220e-01_wp
            hbs( 9,3,4)  = -0.44800e-02_wp
            hbs( 9,4,4)  = 0.29800e-02_wp
            hbs( 9,5,4)  = -0.18500e-02_wp
            hbs( 9,6,4)  = 0.31630e-01_wp
            hbs( 9,7,4)  = 0.44300e-02_wp
            hbs( 9,8,4)  = 0.34300e-02_wp
            hbs( 9,9,4)  = -0.29100e-02_wp
            hbs( 9,10,4) = -0.24050e-01_wp
            hbs( 9,11,4) = 0.98000e-03_wp
            hbs( 9,12,4) = 0.15000e-03_wp
            hbs(10,1,4)  = 0.31600e-01_wp
            hbs(10,2,4)  = 0.31600e-01_wp
            hbs(10,3,4)  = 0.63400e-02_wp
            hbs(10,4,4)  = -0.13710e-01_wp
            hbs(10,5,4)  = -0.34000e-03_wp
            hbs(10,6,4)  = -0.29610e-01_wp
            hbs(10,7,4)  = 0.47800e-02_wp
            hbs(10,8,4)  = -0.29610e-01_wp
            hbs(10,9,4)  = -0.84900e-02_wp
            hbs(10,10,4) = -0.13710e-01_wp
            hbs(10,11,4) = 0.47800e-02_wp
            hbs(10,12,4) = -0.34000e-03_wp
            hbs(11,1,4)  = -0.17850e-01_wp
            hbs(11,2,4)  = -0.13750e-01_wp
            hbs(11,3,4)  = -0.31700e-02_wp
            hbs(11,4,4)  = 0.80700e-02_wp
            hbs(11,5,4)  = 0.73000e-03_wp
            hbs(11,6,4)  = 0.27610e-01_wp
            hbs(11,7,4)  = -0.29000e-02_wp
            hbs(11,8,4)  = 0.20000e-02_wp
            hbs(11,9,4)  = 0.42400e-02_wp
            hbs(11,10,4) = 0.56400e-02_wp
            hbs(11,11,4) = -0.18700e-02_wp
            hbs(11,12,4) = -0.38000e-03_wp
            hbs(12,1,4)  = -0.13750e-01_wp
            hbs(12,2,4)  = -0.17850e-01_wp
            hbs(12,3,4)  = -0.31700e-02_wp
            hbs(12,4,4)  = 0.56400e-02_wp
            hbs(12,5,4)  = -0.38000e-03_wp
            hbs(12,6,4)  = 0.20000e-02_wp
            hbs(12,7,4)  = -0.18700e-02_wp
            hbs(12,8,4)  = 0.27610e-01_wp
            hbs(12,9,4)  = 0.42400e-02_wp
            hbs(12,10,4) = 0.80700e-02_wp
            hbs(12,11,4) = -0.29000e-02_wp
            hbs(12,12,4) = 0.73000e-03_wp
            hbs(13,1,4)  = 0.29170e-01_wp
            hbs(13,2,4)  = -0.78200e-02_wp
            hbs(13,3,4)  = 0.57100e-02_wp
            hbs(13,4,4)  = 0.43210e-01_wp
            hbs(13,5,4)  = -0.32000e-03_wp
            hbs(13,6,4)  = -0.55800e-02_wp
            hbs(13,7,4)  = -0.75000e-03_wp
            hbs(13,8,4)  = -0.35760e-01_wp
            hbs(13,9,4)  = 0.48900e-02_wp
            hbs(13,10,4) = -0.32990e-01_wp
            hbs(13,11,4) = -0.84700e-02_wp
            hbs(13,12,4) = 0.46200e-02_wp
            hbs(14,1,4)  = -0.13220e-01_wp
            hbs(14,2,4)  = 0.46400e-02_wp
            hbs(14,3,4)  = -0.44800e-02_wp
            hbs(14,4,4)  = -0.24050e-01_wp
            hbs(14,5,4)  = 0.15000e-03_wp
            hbs(14,6,4)  = 0.34300e-02_wp
            hbs(14,7,4)  = 0.98000e-03_wp
            hbs(14,8,4)  = 0.31630e-01_wp
            hbs(14,9,4)  = -0.29100e-02_wp
            hbs(14,10,4) = 0.29800e-02_wp
            hbs(14,11,4) = 0.44300e-02_wp
            hbs(14,12,4) = -0.18500e-02_wp
            hbs(15,1,4)  = -0.15950e-01_wp
            hbs(15,2,4)  = 0.31800e-02_wp
            hbs(15,3,4)  = -0.12300e-02_wp
            hbs(15,4,4)  = -0.19160e-01_wp
            hbs(15,5,4)  = 0.17000e-03_wp
            hbs(15,6,4)  = 0.21600e-02_wp
            hbs(15,7,4)  = -0.23000e-03_wp
            hbs(15,8,4)  = 0.41300e-02_wp
            hbs(15,9,4)  = -0.19900e-02_wp
            hbs(15,10,4) = 0.30010e-01_wp
            hbs(15,11,4) = 0.40300e-02_wp
            hbs(15,12,4) = -0.27700e-02_wp
            hbs(16,1,4)  = -0.60700e-02_wp
            hbs(16,2,4)  = -0.24180e-01_wp
            hbs(16,3,4)  = 0.14620e-01_wp
            hbs(16,4,4)  = 0.48130e-01_wp
            hbs(16,5,4)  = -0.12600e-02_wp
            hbs(16,6,4)  = 0.37380e-01_wp
            hbs(16,7,4)  = 0.25000e-03_wp
            hbs(16,8,4)  = -0.10830e-01_wp
            hbs(16,9,4)  = -0.65000e-03_wp
            hbs(16,10,4) = -0.35760e-01_wp
            hbs(16,11,4) = 0.47400e-02_wp
            hbs(16,12,4) = -0.81300e-02_wp
            hbs(17,1,4)  = 0.30300e-02_wp
            hbs(17,2,4)  = 0.21630e-01_wp
            hbs(17,3,4)  = -0.86900e-02_wp
            hbs(17,4,4)  = -0.25720e-01_wp
            hbs(17,5,4)  = 0.12100e-02_wp
            hbs(17,6,4)  = -0.16220e-01_wp
            hbs(17,7,4)  = -0.13000e-03_wp
            hbs(17,8,4)  = 0.52000e-02_wp
            hbs(17,9,4)  = -0.27000e-03_wp
            hbs(17,10,4) = 0.48300e-02_wp
            hbs(17,11,4) = -0.17000e-02_wp
            hbs(17,12,4) = 0.60700e-02_wp
            hbs(18,1,4)  = 0.30300e-02_wp
            hbs(18,2,4)  = 0.25500e-02_wp
            hbs(18,3,4)  = -0.59400e-02_wp
            hbs(18,4,4)  = -0.22410e-01_wp
            hbs(18,5,4)  = 0.50000e-04_wp
            hbs(18,6,4)  = -0.21160e-01_wp
            hbs(18,7,4)  = -0.12000e-03_wp
            hbs(18,8,4)  = 0.56300e-02_wp
            hbs(18,9,4)  = 0.91000e-03_wp
            hbs(18,10,4) = 0.30930e-01_wp
            hbs(18,11,4) = -0.30400e-02_wp
            hbs(18,12,4) = 0.20600e-02_wp
        endif
    end subroutine gethess

    !> @brief Evaluate getgeom for the PHOH_APRP potential.
    subroutine getgeom( nap, qse, nqs, qbe, nqb, qte, nqt, ii )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), parameter :: zero = 0.0e0_wp
        real(wp), parameter :: pi   = 3.1415926535898e0_wp
        ! Global variables
        integer  :: ii
        integer  :: nap
        integer  :: nqs, nqb, nqt
        real(wp) :: qse(nqs,nap), qbe(nqb,nap), qte(nqt,nap)
        ! Local variables
        integer :: i, j
        qte(1,1)  = zero
        qte(2,1)  = pi
        qte(3,1)  = pi
        qte(4,1)  = zero
        qte(5,1)  = zero
        qte(6,1)  = pi
        qte(7,1)  = pi
        qte(8,1)  = zero
        qte(9,1)  = zero
        qte(10,1) = pi
        qte(11,1) = pi
        qte(12,1) = zero
        qte(13,1) = zero
        qte(14,1) = pi
        qte(15,1) = pi
        qte(16,1) = zero
        qte(17,1) = zero
        qte(18,1) = pi
        qte(19,1) = pi
        qte(20,1) = zero
        qte(21,1) = zero
        qte(22,1) = pi
        qte(23,1) = pi
        qte(24,1) = zero
        do i=2, nap
            do j=1, nqt
                qte(j,i) = qte(j,1)
            enddo
        enddo
        ! Now get equilibrium geometrical parameters for qae
        if( ii .eq. 11 )then
            ! equilibrium geometry at anchor point 1
            qse( 1,1) = 0.263789e+01_wp
            qse( 2,1) = 0.263870e+01_wp
            qse( 3,1) = 0.258038e+01_wp
            qse( 4,1) = 0.263573e+01_wp
            qse( 5,1) = 0.206310e+01_wp
            qse( 6,1) = 0.263336e+01_wp
            qse( 7,1) = 0.205963e+01_wp
            qse( 8,1) = 0.263920e+01_wp
            qse( 9,1) = 0.205742e+01_wp
            qse(10,1) = 0.262993e+01_wp
            qse(11,1) = 0.205975e+01_wp
            qse(12,1) = 0.205790e+01_wp
            ! equilibrium geometry at anchor point 1
            qbe( 1,1) = 0.210140e+01_wp
            qbe( 2,1) = 0.213715e+01_wp
            qbe( 3,1) = 0.204464e+01_wp
            qbe( 4,1) = 0.208658e+01_wp
            qbe( 5,1) = 0.209471e+01_wp
            qbe( 6,1) = 0.210190e+01_wp
            qbe( 7,1) = 0.210459e+01_wp
            qbe( 8,1) = 0.208211e+01_wp
            qbe( 9,1) = 0.209649e+01_wp
            qbe(10,1) = 0.208155e+01_wp
            qbe(11,1) = 0.210036e+01_wp
            qbe(12,1) = 0.210127e+01_wp
            qbe(13,1) = 0.210834e+01_wp
            qbe(14,1) = 0.209320e+01_wp
            qbe(15,1) = 0.208164e+01_wp
            qbe(16,1) = 0.208392e+01_wp
            qbe(17,1) = 0.207495e+01_wp
            qbe(18,1) = 0.212431e+01_wp
            ! equilibrium geometry at anchor point 2
            qse( 1,2) = 0.263724e+01_wp
            qse( 2,2) = 0.264694e+01_wp
            qse( 3,2) = 0.257138e+01_wp
            qse( 4,2) = 0.263516e+01_wp
            qse( 5,2) = 0.206589e+01_wp
            qse( 6,2) = 0.263264e+01_wp
            qse( 7,2) = 0.205963e+01_wp
            qse( 8,2) = 0.263972e+01_wp
            qse( 9,2) = 0.205726e+01_wp
            qse(10,2) = 0.262912e+01_wp
            qse(11,2) = 0.205993e+01_wp
            qse(12,2) = 0.205801e+01_wp
            ! equilibrium geometry at anchor point 2
            qbe( 1,2) = 0.209269e+01_wp
            qbe( 2,2) = 0.215700e+01_wp
            qbe( 3,2) = 0.203350e+01_wp
            qbe( 4,2) = 0.209468e+01_wp
            qbe( 5,2) = 0.208443e+01_wp
            qbe( 6,2) = 0.210408e+01_wp
            qbe( 7,2) = 0.210277e+01_wp
            qbe( 8,2) = 0.208276e+01_wp
            qbe( 9,2) = 0.209766e+01_wp
            qbe(10,2) = 0.207979e+01_wp
            qbe(11,2) = 0.210131e+01_wp
            qbe(12,2) = 0.210209e+01_wp
            qbe(13,2) = 0.211103e+01_wp
            qbe(14,2) = 0.209320e+01_wp
            qbe(15,2) = 0.207895e+01_wp
            qbe(16,2) = 0.208541e+01_wp
            qbe(17,2) = 0.207654e+01_wp
            qbe(18,2) = 0.212123e+01_wp
            ! equilibrium geometry at anchor point 3
            qse( 1,3) = 0.264670e+01_wp
            qse( 2,3) = 0.266084e+01_wp
            qse( 3,3) = 0.252665e+01_wp
            qse( 4,3) = 0.263224e+01_wp
            qse( 5,3) = 0.206010e+01_wp
            qse( 6,3) = 0.263214e+01_wp
            qse( 7,3) = 0.205982e+01_wp
            qse( 8,3) = 0.263685e+01_wp
            qse( 9,3) = 0.205662e+01_wp
            qse(10,3) = 0.262787e+01_wp
            qse(11,3) = 0.206007e+01_wp
            qse(12,3) = 0.205634e+01_wp
            ! equilibrium geometry at anchor point 3
            qbe( 1,3) = 0.209445e+01_wp
            qbe( 2,3) = 0.215439e+01_wp
            qbe( 3,3) = 0.203434e+01_wp
            qbe( 4,3) = 0.208659e+01_wp
            qbe( 5,3) = 0.207835e+01_wp
            qbe( 6,3) = 0.211824e+01_wp
            qbe( 7,3) = 0.211238e+01_wp
            qbe( 8,3) = 0.207426e+01_wp
            qbe( 9,3) = 0.209655e+01_wp
            qbe(10,3) = 0.207613e+01_wp
            qbe(11,3) = 0.210353e+01_wp
            qbe(12,3) = 0.210352e+01_wp
            qbe(13,3) = 0.211662e+01_wp
            qbe(14,3) = 0.209569e+01_wp
            qbe(15,3) = 0.207088e+01_wp
            qbe(16,3) = 0.208020e+01_wp
            qbe(17,3) = 0.208110e+01_wp
            qbe(18,3) = 0.212189e+01_wp
            ! equilibrium geometry at anchor point 4
            qse( 1,4) = 0.265159e+01_wp
            qse( 2,4) = 0.265159e+01_wp
            qse( 3,4) = 0.251838e+01_wp
            qse( 4,4) = 0.263062e+01_wp
            qse( 5,4) = 0.205598e+01_wp
            qse( 6,4) = 0.263458e+01_wp
            qse( 7,4) = 0.205967e+01_wp
            qse( 8,4) = 0.263458e+01_wp
            qse( 9,4) = 0.205648e+01_wp
            qse(10,4) = 0.263062e+01_wp
            qse(11,4) = 0.205967e+01_wp
            qse(12,4) = 0.205598e+01_wp
            ! equilibrium geometry at anchor point 4
            qbe( 1,4) = 0.211210e+01_wp
            qbe( 2,4) = 0.208554e+01_wp
            qbe( 3,4) = 0.208554e+01_wp
            qbe( 4,4) = 0.207166e+01_wp
            qbe( 5,4) = 0.208711e+01_wp
            qbe( 6,4) = 0.212441e+01_wp
            qbe( 7,4) = 0.211593e+01_wp
            qbe( 8,4) = 0.207134e+01_wp
            qbe( 9,4) = 0.209592e+01_wp
            qbe(10,4) = 0.207909e+01_wp
            qbe(11,4) = 0.210205e+01_wp
            qbe(12,4) = 0.210205e+01_wp
            qbe(13,4) = 0.211593e+01_wp
            qbe(14,4) = 0.209592e+01_wp
            qbe(15,4) = 0.207134e+01_wp
            qbe(16,4) = 0.207166e+01_wp
            qbe(17,4) = 0.208711e+01_wp
            qbe(18,4) = 0.212441e+01_wp
        elseif( ii .eq. 22 )then
            ! equilibrium geometry at anchor point 1
            qse( 1,1) = 0.269941e+01_wp
            qse( 2,1) = 0.268363e+01_wp
            qse( 3,1) = 0.252916e+01_wp
            qse( 4,1) = 0.268093e+01_wp
            qse( 5,1) = 0.205887e+01_wp
            qse( 6,1) = 0.268592e+01_wp
            qse( 7,1) = 0.205454e+01_wp
            qse( 8,1) = 0.268227e+01_wp
            qse( 9,1) = 0.206068e+01_wp
            qse(10,1) = 0.268469e+01_wp
            qse(11,1) = 0.205375e+01_wp
            qse(12,1) = 0.205424e+01_wp
            ! equilibrium geometry at anchor point 1
            qbe( 1,1) = 0.216375e+01_wp
            qbe( 2,1) = 0.209466e+01_wp
            qbe( 3,1) = 0.202477e+01_wp
            qbe( 4,1) = 0.204773e+01_wp
            qbe( 5,1) = 0.209602e+01_wp
            qbe( 6,1) = 0.213944e+01_wp
            qbe( 7,1) = 0.208156e+01_wp
            qbe( 8,1) = 0.210850e+01_wp
            qbe( 9,1) = 0.209313e+01_wp
            qbe(10,1) = 0.214184e+01_wp
            qbe(11,1) = 0.206797e+01_wp
            qbe(12,1) = 0.207337e+01_wp
            qbe(13,1) = 0.207509e+01_wp
            qbe(14,1) = 0.209807e+01_wp
            qbe(15,1) = 0.211003e+01_wp
            qbe(16,1) = 0.205641e+01_wp
            qbe(17,1) = 0.207231e+01_wp
            qbe(18,1) = 0.215446e+01_wp
            ! equilibrium geometry at anchor point 2
            qse( 1,2) = 0.269988e+01_wp
            qse( 2,2) = 0.269065e+01_wp
            qse( 3,2) = 0.251103e+01_wp
            qse( 4,2) = 0.267864e+01_wp
            qse( 5,2) = 0.206339e+01_wp
            qse( 6,2) = 0.268396e+01_wp
            qse( 7,2) = 0.205439e+01_wp
            qse( 8,2) = 0.268024e+01_wp
            qse( 9,2) = 0.206089e+01_wp
            qse(10,2) = 0.268483e+01_wp
            qse(11,2) = 0.205390e+01_wp
            qse(12,2) = 0.205412e+01_wp
            ! equilibrium geometry at anchor point 2
            qbe( 1,2) = 0.215707e+01_wp
            qbe( 2,2) = 0.211094e+01_wp
            qbe( 3,2) = 0.201518e+01_wp
            qbe( 4,2) = 0.205661e+01_wp
            qbe( 5,2) = 0.208168e+01_wp
            qbe( 6,2) = 0.214490e+01_wp
            qbe( 7,2) = 0.207623e+01_wp
            qbe( 8,2) = 0.211213e+01_wp
            qbe( 9,2) = 0.209483e+01_wp
            qbe(10,2) = 0.214363e+01_wp
            qbe(11,2) = 0.206707e+01_wp
            qbe(12,2) = 0.207249e+01_wp
            qbe(13,2) = 0.207837e+01_wp
            qbe(14,2) = 0.209780e+01_wp
            qbe(15,2) = 0.210701e+01_wp
            qbe(16,2) = 0.205447e+01_wp
            qbe(17,2) = 0.207468e+01_wp
            qbe(18,2) = 0.215403e+01_wp
            ! equilibrium geometry at anchor point 3
            qse( 1,3) = 0.269988e+01_wp
            qse( 2,3) = 0.269065e+01_wp
            qse( 3,3) = 0.251103e+01_wp
            qse( 4,3) = 0.267864e+01_wp
            qse( 5,3) = 0.206339e+01_wp
            qse( 6,3) = 0.268396e+01_wp
            qse( 7,3) = 0.205439e+01_wp
            qse( 8,3) = 0.268024e+01_wp
            qse( 9,3) = 0.206089e+01_wp
            qse(10,3) = 0.268483e+01_wp
            qse(11,3) = 0.205390e+01_wp
            qse(12,3) = 0.205412e+01_wp
            ! equilibrium geometry at anchor point 3
            qbe( 1,3) = 0.215707e+01_wp
            qbe( 2,3) = 0.211094e+01_wp
            qbe( 3,3) = 0.201518e+01_wp
            qbe( 4,3) = 0.205661e+01_wp
            qbe( 5,3) = 0.208168e+01_wp
            qbe( 6,3) = 0.214490e+01_wp
            qbe( 7,3) = 0.207623e+01_wp
            qbe( 8,3) = 0.211213e+01_wp
            qbe( 9,3) = 0.209483e+01_wp
            qbe(10,3) = 0.214363e+01_wp
            qbe(11,3) = 0.206707e+01_wp
            qbe(12,3) = 0.207249e+01_wp
            qbe(13,3) = 0.207837e+01_wp
            qbe(14,3) = 0.209780e+01_wp
            qbe(15,3) = 0.210701e+01_wp
            qbe(16,3) = 0.205447e+01_wp
            qbe(17,3) = 0.207468e+01_wp
            qbe(18,3) = 0.215403e+01_wp
            ! equilibrium geometry at anchor point 4
            qse( 1,4) = 0.269988e+01_wp
            qse( 2,4) = 0.269065e+01_wp
            qse( 3,4) = 0.251103e+01_wp
            qse( 4,4) = 0.267864e+01_wp
            qse( 5,4) = 0.206339e+01_wp
            qse( 6,4) = 0.268396e+01_wp
            qse( 7,4) = 0.205439e+01_wp
            qse( 8,4) = 0.268024e+01_wp
            qse( 9,4) = 0.206089e+01_wp
            qse(10,4) = 0.268483e+01_wp
            qse(11,4) = 0.205390e+01_wp
            qse(12,4) = 0.205412e+01_wp
            ! equilibrium geometry at anchor point 4
            qbe( 1,4) = 0.215707e+01_wp
            qbe( 2,4) = 0.211094e+01_wp
            qbe( 3,4) = 0.201518e+01_wp
            qbe( 4,4) = 0.205661e+01_wp
            qbe( 5,4) = 0.208168e+01_wp
            qbe( 6,4) = 0.214490e+01_wp
            qbe( 7,4) = 0.207623e+01_wp
            qbe( 8,4) = 0.211213e+01_wp
            qbe( 9,4) = 0.209483e+01_wp
            qbe(10,4) = 0.214363e+01_wp
            qbe(11,4) = 0.206707e+01_wp
            qbe(12,4) = 0.207249e+01_wp
            qbe(13,4) = 0.207837e+01_wp
            qbe(14,4) = 0.209780e+01_wp
            qbe(15,4) = 0.210701e+01_wp
            qbe(16,4) = 0.205447e+01_wp
            qbe(17,4) = 0.207468e+01_wp
            qbe(18,4) = 0.215403e+01_wp
        elseif( ii .eq. 33 )then
            ! equilibrium geometry at anchor point 1
            qse( 1,1) = 0.272760e+01_wp
            qse( 2,1) = 0.272219e+01_wp
            qse( 3,1) = 0.243093e+01_wp
            qse( 4,1) = 0.259142e+01_wp
            qse( 5,1) = 0.205558e+01_wp
            qse( 6,1) = 0.266713e+01_wp
            qse( 7,1) = 0.205707e+01_wp
            qse( 8,1) = 0.267731e+01_wp
            qse( 9,1) = 0.205961e+01_wp
            qse(10,1) = 0.258566e+01_wp
            qse(11,1) = 0.205666e+01_wp
            qse(12,1) = 0.205567e+01_wp
            ! equilibrium geometry at anchor point 1
            qbe( 1,1) = 0.209913e+01_wp
            qbe( 2,1) = 0.214547e+01_wp
            qbe( 3,1) = 0.203859e+01_wp
            qbe( 4,1) = 0.207342e+01_wp
            qbe( 5,1) = 0.205792e+01_wp
            qbe( 6,1) = 0.215185e+01_wp
            qbe( 7,1) = 0.210210e+01_wp
            qbe( 8,1) = 0.209263e+01_wp
            qbe( 9,1) = 0.208845e+01_wp
            qbe(10,1) = 0.211486e+01_wp
            qbe(11,1) = 0.208558e+01_wp
            qbe(12,1) = 0.208274e+01_wp
            qbe(13,1) = 0.209348e+01_wp
            qbe(14,1) = 0.208777e+01_wp
            qbe(15,1) = 0.210193e+01_wp
            qbe(16,1) = 0.208337e+01_wp
            qbe(17,1) = 0.205127e+01_wp
            qbe(18,1) = 0.214854e+01_wp
            ! equilibrium geometry at anchor point 2
            qse( 1,2) = 0.274917e+01_wp
            qse( 2,2) = 0.274193e+01_wp
            qse( 3,2) = 0.237981e+01_wp
            qse( 4,2) = 0.259572e+01_wp
            qse( 5,2) = 0.205541e+01_wp
            qse( 6,2) = 0.266481e+01_wp
            qse( 7,2) = 0.205872e+01_wp
            qse( 8,2) = 0.266882e+01_wp
            qse( 9,2) = 0.205898e+01_wp
            qse(10,2) = 0.259216e+01_wp
            qse(11,2) = 0.205873e+01_wp
            qse(12,2) = 0.205800e+01_wp
            ! equilibrium geometry at anchor point 2
            qbe( 1,2) = 0.206394e+01_wp
            qbe( 2,2) = 0.213351e+01_wp
            qbe( 3,2) = 0.208574e+01_wp
            qbe( 4,2) = 0.208869e+01_wp
            qbe( 5,2) = 0.205086e+01_wp
            qbe( 6,2) = 0.214363e+01_wp
            qbe( 7,2) = 0.210572e+01_wp
            qbe( 8,2) = 0.209443e+01_wp
            qbe( 9,2) = 0.208304e+01_wp
            qbe(10,2) = 0.211155e+01_wp
            qbe(11,2) = 0.208615e+01_wp
            qbe(12,2) = 0.208549e+01_wp
            qbe(13,2) = 0.209340e+01_wp
            qbe(14,2) = 0.208795e+01_wp
            qbe(15,2) = 0.210183e+01_wp
            qbe(16,2) = 0.210307e+01_wp
            qbe(17,2) = 0.204251e+01_wp
            qbe(18,2) = 0.213761e+01_wp
            ! equilibrium geometry at anchor point 3
            qse( 1,3) = 0.274751e+01_wp
            qse( 2,3) = 0.274706e+01_wp
            qse( 3,3) = 0.235956e+01_wp
            qse( 4,3) = 0.259829e+01_wp
            qse( 5,3) = 0.205917e+01_wp
            qse( 6,3) = 0.266445e+01_wp
            qse( 7,3) = 0.205922e+01_wp
            qse( 8,3) = 0.266458e+01_wp
            qse( 9,3) = 0.205897e+01_wp
            qse(10,3) = 0.259790e+01_wp
            qse(11,3) = 0.205917e+01_wp
            qse(12,3) = 0.205900e+01_wp
            ! equilibrium geometry at anchor point 3
            qbe( 1,3) = 0.204756e+01_wp
            qbe( 2,3) = 0.212154e+01_wp
            qbe( 3,3) = 0.211409e+01_wp
            qbe( 4,3) = 0.210600e+01_wp
            qbe( 5,3) = 0.204389e+01_wp
            qbe( 6,3) = 0.213330e+01_wp
            qbe( 7,3) = 0.209805e+01_wp
            qbe( 8,3) = 0.210013e+01_wp
            qbe( 9,3) = 0.208501e+01_wp
            qbe(10,3) = 0.211048e+01_wp
            qbe(11,3) = 0.208627e+01_wp
            qbe(12,3) = 0.208644e+01_wp
            qbe(13,3) = 0.209637e+01_wp
            qbe(14,3) = 0.208569e+01_wp
            qbe(15,3) = 0.210112e+01_wp
            qbe(16,3) = 0.210790e+01_wp
            qbe(17,3) = 0.204012e+01_wp
            qbe(18,3) = 0.213517e+01_wp
            ! equilibrium geometry at anchor point 4
            qse( 1,4) = 0.274535e+01_wp
            qse( 2,4) = 0.274535e+01_wp
            qse( 3,4) = 0.236001e+01_wp
            qse( 4,4) = 0.259897e+01_wp
            qse( 5,4) = 0.205919e+01_wp
            qse( 6,4) = 0.266411e+01_wp
            qse( 7,4) = 0.205914e+01_wp
            qse( 8,4) = 0.266411e+01_wp
            qse( 9,4) = 0.205909e+01_wp
            qse(10,4) = 0.259897e+01_wp
            qse(11,4) = 0.205914e+01_wp
            qse(12,4) = 0.205919e+01_wp
            ! equilibrium geometry at anchor point 4
            qbe( 1,4) = 0.204757e+01_wp
            qbe( 2,4) = 0.211781e+01_wp
            qbe( 3,4) = 0.211781e+01_wp
            qbe( 4,4) = 0.210758e+01_wp
            qbe( 5,4) = 0.204104e+01_wp
            qbe( 6,4) = 0.213456e+01_wp
            qbe( 7,4) = 0.209631e+01_wp
            qbe( 8,4) = 0.210152e+01_wp
            qbe( 9,4) = 0.208536e+01_wp
            qbe(10,4) = 0.211101e+01_wp
            qbe(11,4) = 0.208609e+01_wp
            qbe(12,4) = 0.208609e+01_wp
            qbe(13,4) = 0.209631e+01_wp
            qbe(14,4) = 0.208536e+01_wp
            qbe(15,4) = 0.210152e+01_wp
            qbe(16,4) = 0.210758e+01_wp
            qbe(17,4) = 0.204104e+01_wp
            qbe(18,4) = 0.213456e+01_wp
        endif
    end subroutine getgeom

    !> @brief Evaluate tent for the PHOH_APRP potential.
    subroutine tent( nap, rap, r, ir, t, dtdr )
        ! Subroutine to evaluate tent functions
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: dudr, dvdr
        integer  :: ir
        integer  :: nap
        real(wp) :: r
        real(wp) :: rap(*), t(*), dtdr(*)
        ! Local variables
        integer  :: i
        real(wp) :: u, v
        do i=1, nap
            t(i)    = 0.0e0_wp
            dtdr(i) = 0.0e0_wp
        enddo
        ! Check whether the point are end points
        if( r .lt. rap(1) )then
            ir      = 1
            t(1)    = 1.0e0_wp
            dtdr(1) = 0.0e0_wp
        elseif( r .gt. rap(nap) )then
            ir        = nap+1
            t(nap)    = 1.0e0_wp
            dtdr(nap) = 0.0e0_wp
        endif
        ! Determine i for which r(i-1) <= r < r(i)
        do i=2, nap
            if( (rap(i-1) .le. r) .and. (r .lt. rap(i)) )then
                ir        = i
                u         = (r-rap(i-1))*(r-rap(i-1))*(r-rap(i-1))*(r-rap(i-1))
                v         = (r-rap(i))*(r-rap(i))*(r-rap(i))*(r-rap(i))
                dudr      = 4*(r-rap(i-1))*(r-rap(i-1))*(r-rap(i-1))
                dvdr      = 4*(r-rap(i))*(r-rap(i))*(r-rap(i))
                t(i-1)    = v/(u+v)
                t(i)      = u/(u+v)
                dtdr(i-1) = (u*dvdr-v*dudr)/((u+v)*(u+v))
                dtdr(i)   = (v*dudr-u*dvdr)/((u+v)*(u+v))
                exit
            endif
        enddo
    end subroutine tent

    !> @brief Evaluate evu1r for the PHOH_APRP potential.
    subroutine evu1r( igrad, r, v, dvdr )
        ! Subroutine to evaluate Varshini potential and derivative for given r
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), parameter :: de   = 0.17596055e0_wp
        real(wp), parameter :: re   = 1.93206106e0_wp
        real(wp), parameter :: beta = 0.63822703e0_wp
        integer             :: igrad
        real(wp)            :: r, v, dvdr
        real(wp)            :: x, y, dvdy, dydx, dxdr
        x = r/re
        y = 1.0e0_wp-exp(-beta*(x*x-1.0e0_wp))/x
        v = de*y*y
        if( igrad .eq. 1 )then
            dvdy = 2.0e0_wp*de*y
            dydx = (2.0e0_wp*beta*x*x+1.0e0_wp)*exp(-beta*(x*x-1.0e0_wp))/(x*x)
            dxdr = 1.0e0_wp/re
            dvdr = dvdy*dydx*dxdr
        endif
    end subroutine evu1r

    !> @brief Evaluate evu1b for the PHOH_APRP potential.
    subroutine evu1b( igrad, r, theta, vbnd, dvdr, dvdth )
        ! Subroutine to evaluate potential energy from H-O-C bending
        ! and derivatives for given rOH and aHOC
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)           :: cthta, cthta0, dct0dr, dk2dr, dthta, dtxdr, dtxdx, dvdr, dvdth, dxdr, gr, sthta, theta, tx, vbnd, x
        integer            :: igrad
        integer, parameter :: ngk2 = 1
        real(wp)           :: a2(ngk2)
        real(wp)           :: alpha2(ngk2)
        real(wp)           :: r2(ngk2)
        real(wp)           :: cthta1, cthta2, a, r0
        integer            :: i, j, k
        real(wp)           :: r
        real(wp)           :: k2
        a2     = 0.10449356e0_wp
        alpha2 = 0.27051185e0_wp
        r2     = 2.28431658e0_wp
        cthta1 = -0.39964587e0_wp
        cthta2 = -0.24769126e0_wp
        a      = 4.01097890e0_wp
        r0     = 1.83449325e0_wp
        cthta  = cos(theta)
        sthta  = sin(theta)
        x      = a*(r-r0)
        !      tx=(exp(x)-exp(-x))/(exp(x)+exp(-x))
        tx     = (1.0e0_wp-exp(-2.0e0_wp*x))/(1.0e0_wp+exp(-2.0e0_wp*x))
        cthta0 = cthta1+0.5e0_wp*(1+tx)*(cthta2-cthta1)
        dxdr   = a
        dtxdx  = 4.0e0_wp/(2.0e0_wp+exp(2.0e0_wp*x)+exp(-2.0e0_wp*x))
        dtxdr  = dtxdx*dxdr
        dct0dr = 0.5e0_wp*(cthta2-cthta1)*dtxdr
        k2     = 0.0e0_wp
        dk2dr  = 0.0e0_wp
        do i=1, ngk2
            gr = a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
            k2 = k2+gr
            if( igrad .eq. 1 )then
                dk2dr = dk2dr-2.0e0_wp*alpha2(i)*(r-r2(i))*gr
            endif
        enddo
        dthta = cthta-cthta0
        vbnd  = k2*dthta*dthta
        if( igrad .eq. 1 )then
            dvdr  = dk2dr*dthta*dthta - 2.0e0_wp*k2*dthta*dct0dr
            dvdth = -2.0e0_wp*k2*dthta*sthta
        endif
    end subroutine evu1b

    !> @brief Evaluate evu1to for the PHOH_APRP potential.
    subroutine evu1to( igrad, r, phi, v, dvdr, dvdph )
        ! Subroutine to evaluate torsion potentials and derivatives
        ! for given r and phi
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: dx1dr, dx2dr, dydph, x1, x2
        real(wp), parameter :: a1     = 0.02371114e0_wp
        real(wp), parameter :: alpha1 = 0.96558083e0_wp
        real(wp), parameter :: r1     = 3.25505778e0_wp
        real(wp), parameter :: a2     = 0.00349299e0_wp
        real(wp), parameter :: alpha2 = 1.092945551e1_wp
        real(wp), parameter :: r2     = 2.10516007e0_wp
        integer             :: igrad
        real(wp)            :: phi, r, v, dvdph, dvdr
        real(wp)            :: x, y, dvdy, dydx, dxdr
        x1 = a1*exp(-alpha1*(r-r1)*(r-r1))
        x2 = a2*exp(-alpha2*(r-r2)*(r-r2))
        y  = 1.0e0_wp-cos(2.0e0_wp*phi)
        v  = x1*y+x2*y*y
        if( igrad .eq. 1 )then
            dx1dr = -2.0e0_wp*alpha1*(r-r1)*x1
            dx2dr = -2.0e0_wp*alpha2*(r-r2)*x2
            dydph = 2.0e0_wp*sin(2.0e0_wp*phi)
            dvdr  = dx1dr*y+dx2dr*y*y
            dvdph = x1*dydph+2.0e0_wp*x2*y*dydph
        endif
    end subroutine evu1to

    !> @brief Evaluate evu2r for the PHOH_APRP potential.
    subroutine evu2r( igrad, r, v, dvdr )
        ! Subroutine to evaluate Morse potential and derivative for given r
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), parameter :: alpha = 0.67269344e0_wp
        real(wp), parameter :: de    = 0.41866243e0_wp
        real(wp), parameter :: re    = 1.95347670e0_wp
        real(wp), parameter :: v0    = 0.17280158e0_wp
        integer             :: igrad
        real(wp)            :: r, v, dvdr
        real(wp)            :: x, y, dvdy, dydx, dxdr
        x = r-re
        y = 1.0e0_wp-exp(-alpha*x)
        v = de*y*y+v0
        if( igrad .eq. 1 )then
            dvdy = 2.0e0_wp*de*y
            dydx = alpha*exp(-alpha*x)
            dxdr = 1.0e0_wp
            dvdr = dvdy*dydx*dxdr
        endif
    end subroutine evu2r

    !> @brief Evaluate evu2b for the PHOH_APRP potential.
    subroutine evu2b( igrad, r, theta, vbnd, dvdr, dvdth )
        ! Subroutine to evaluate potential energy from H-O-C bending
        ! and derivatives for given rOH and aHOC
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)           :: cthta, cthta0, dct0dr, dk2dr, dthta, dtxdr, dtxdx, dvdr, dvdth, dxdr, gr, sthta, theta, tx, vbnd, x
        integer            :: igrad
        integer, parameter :: ngk2 = 1
        real(wp)           :: a2(ngk2)
        real(wp)           :: alpha2(ngk2)
        real(wp)           :: r2(ngk2)
        real(wp)           :: cthta1, cthta2, a, r0
        integer            :: i, j, k
        real(wp)           :: r
        real(wp)           :: k2
        a2     = 0.10090668e0_wp
        alpha2 = 0.21154184e0_wp
        r2     = 1.80972781e0_wp
        cthta1 = -0.31835353e0_wp
        cthta2 = -0.22230161e0_wp
        a      = 4.73544047e0_wp
        r0     = 2.01933829e0_wp
        cthta  = cos(theta)
        sthta  = sin(theta)
        x      = a*(r-r0)
        !      tx=(exp(x)-exp(-x))/(exp(x)+exp(-x))
        tx     = (1.0e0_wp-exp(-2.0e0_wp*x))/(1.0e0_wp+exp(-2.0e0_wp*x))
        cthta0 = cthta1+0.5e0_wp*(1+tx)*(cthta2-cthta1)
        dxdr   = a
        dtxdx  = 4.0e0_wp/(2.0e0_wp+exp(2.0e0_wp*x)+exp(-2.0e0_wp*x))
        dtxdr  = dtxdx*dxdr
        dct0dr = 0.5e0_wp*(cthta2-cthta1)*dtxdr
        k2     = 0.0e0_wp
        dk2dr  = 0.0e0_wp
        do i=1, ngk2
            gr = a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
            k2 = k2+gr
            if( igrad .eq. 1 )then
                dk2dr = dk2dr-2.0e0_wp*alpha2(i)*(r-r2(i))*gr
            endif
        enddo
        dthta = cthta-cthta0
        vbnd  = k2*dthta*dthta
        if( igrad .eq. 1 )then
            dvdr  = dk2dr*dthta*dthta - 2.0e0_wp*k2*dthta*dct0dr
            dvdth = -2.0e0_wp*k2*dthta*sthta
        endif
    end subroutine evu2b

    !> @brief Evaluate evu2to for the PHOH_APRP potential.
    subroutine evu2to( igrad, r, phi, v, dvdr, dvdph )
        ! Subroutine to evaluate torsion potentials and derivatives
        ! for given r and phi
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: dydph
        real(wp), parameter :: a1     = 0.01925782e0_wp
        real(wp), parameter :: alpha1 = 5.00873578e0_wp
        real(wp), parameter :: r1     = 2.22678700e0_wp
        integer             :: igrad
        real(wp)            :: phi, r, v, dvdph, dvdr
        real(wp)            :: x, y, dvdy, dydx, dxdr
        x = a1*exp(-alpha1*(r-r1)*(r-r1))
        y = 1.0e0_wp-cos(2.0e0_wp*phi)
        v = x*y
        if( igrad .eq. 1 )then
            dxdr  = -2.0e0_wp*alpha1*(r-r1)*x
            dydph = 2.0e0_wp*sin(2.0e0_wp*phi)
            dvdr  = y*dxdr
            dvdph = x*dydph
        endif
    end subroutine evu2to

    !> @brief Evaluate evu3r for the PHOH_APRP potential.
    subroutine evu3r( igrad, r, v, dvdr )
        ! Subroutine to evaluate Morse potential and derivative for given r
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: dx1dr, dx2dr, dx3dr, dy1dx1, dy2dx2, dy3dx3, x1, x2, x3, y1, y2, y3
        real(wp), parameter :: v0     = 0.15976218e0_wp
        real(wp), parameter :: a1     = 0.06496177e0_wp
        real(wp), parameter :: alpha1 = 1.81347726e0_wp
        real(wp), parameter :: r1     = 2.63168269e0_wp
        real(wp), parameter :: a2     = -0.00077880e0_wp
        real(wp), parameter :: alpha2 = 3.24769899e0_wp
        real(wp), parameter :: r2     = 3.71412291e0_wp
        real(wp), parameter :: a3     = 9.3151142407e2_wp
        real(wp), parameter :: alpha3 = 4.73016534e0_wp
        real(wp), parameter :: r3     = -0.04117486e0_wp
        integer             :: igrad
        real(wp)            :: r, v, dvdr
        real(wp)            :: x, y, dvdy, dydx, dxdr
        x1 = r-r1
        y1 = a1*exp(-alpha1*x1)
        x2 = r-r2
        y2 = a2*exp(-alpha2*x2)
        x3 = r-r3
        y3 = a3*exp(-alpha3*x3)
        v  = v0+y1+y2+y3
        if( igrad .eq. 1 )then
            dx1dr  = 1.0e0_wp
            dx2dr  = 1.0e0_wp
            dx3dr  = 1.0e0_wp
            dy1dx1 = -alpha1*y1
            dy2dx2 = -alpha2*y2
            dy3dx3 = -alpha3*y3
            dvdr   = dy1dx1*dx1dr+dy2dx2*dx2dr+dy3dx3*dx3dr
        endif
    end subroutine evu3r

    !> @brief Evaluate evu3b for the PHOH_APRP potential.
    subroutine evu3b( igrad, r, theta, vbnd, dvdr, dvdth )
        ! Subroutine to evaluate potential energy from H-O-C bending
        ! and derivatives for given rOH and aHOC
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)           :: cthta, cthta0, dct0dr, dk2dr, dk3dr, dthta, dtxdr, dtxdx, dvdr, dvdth, dxdr, gr, sthta, theta, tx, vbnd, x
        integer            :: igrad
        integer, parameter :: ngk2 = 1
        integer, parameter :: ngk3 = 1
        real(wp)           :: a2(ngk2), a3(ngk3)
        real(wp)           :: alpha2(ngk2), alpha3(ngk3)
        real(wp)           :: r2(ngk2), r3(ngk3)
        real(wp)           :: cthta1, cthta2, a, r0
        integer            :: i, j, k
        real(wp)           :: r
        real(wp)           :: k2, k3
        a2     = 3.00968695e0_wp
        a3     = 0.06538039e0_wp
        alpha2 = 0.04954325e0_wp
        alpha3 = 14.94793100e0_wp
        r2     = -6.33258333e0_wp
        r3     = 2.09315608e0_wp
        cthta1 = -0.30362526e0_wp
        cthta2 = -0.65443171e0_wp
        a      = 0.90134535e0_wp
        r0     = 2.95896455e0_wp
        cthta  = cos(theta)
        sthta  = sin(theta)
        x      = a*(r-r0)
        !      tx=(exp(x)-exp(-x))/(exp(x)+exp(-x))
        tx     = (1.0e0_wp-exp(-2.0e0_wp*x))/(1.0e0_wp+exp(-2.0e0_wp*x))
        cthta0 = cthta1+0.5e0_wp*(1+tx)*(cthta2-cthta1)
        dxdr   = a
        dtxdx  = 4.0e0_wp/(2.0e0_wp+exp(2.0e0_wp*x)+exp(-2.0e0_wp*x))
        dtxdr  = dtxdx*dxdr
        dct0dr = 0.5e0_wp*(cthta2-cthta1)*dtxdr
        k2     = 0.0e0_wp
        dk2dr  = 0.0e0_wp
        do i=1, ngk2
            gr = a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
            k2 = k2+gr
            if( igrad .eq. 1 )then
                dk2dr = dk2dr-2.0e0_wp*alpha2(i)*(r-r2(i))*gr
            endif
        enddo
        k3    = 0.0e0_wp
        dk3dr = 0.0e0_wp
        do i=1, ngk3
            gr = a3(i)*exp(-alpha3(i)*(r-r3(i))*(r-r3(i)))
            k3 = k3+gr
            if( igrad .eq. 1 )then
                dk3dr = dk3dr-2.0e0_wp*alpha3(i)*(r-r3(i))*gr
            endif
        enddo
        dthta = cthta-cthta0
        vbnd  = k2*dthta*dthta + k3*dthta*dthta*dthta
        if( igrad .eq. 1 )then
            dvdr  = dk2dr*dthta*dthta + dk3dr*dthta*dthta*dthta - 2.0e0_wp*k2*dthta*dct0dr - 3.0e0_wp*k3*dthta*dthta* &
                & dct0dr
            dvdth = -2.0e0_wp*k2*dthta*sthta - 3.0e0_wp*k3*dthta*dthta*sthta
        endif
    end subroutine evu3b

    !> @brief Evaluate evu3to for the PHOH_APRP potential.
    subroutine evu3to( igrad, r, phi, v, dvdr, dvdph )
        ! Subroutine to evaluate torsion potentials and derivatives
        ! for given r and phi
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: dx1dr, dx2dr, dydph, x1, x2
        real(wp), parameter :: a1     = -0.02298368e0_wp
        real(wp), parameter :: alpha1 = 1.06740430e0_wp
        real(wp), parameter :: r1     = 3.28664483e0_wp
        real(wp), parameter :: a2     = 0.01138152e0_wp
        real(wp), parameter :: alpha2 = 0.65339172e0_wp
        real(wp), parameter :: r2     = 2.58341105e0_wp
        integer             :: igrad
        real(wp)            :: phi, r, v, dvdph, dvdr
        real(wp)            :: x, y, dvdy, dydx, dxdr
        x1 = a1*exp(-alpha1*(r-r1)*(r-r1))
        x2 = a2*exp(-alpha2*(r-r2)*(r-r2))
        y  = 1.0e0_wp-cos(2.0e0_wp*phi)
        v  = x1*y+x2*y*y
        if( igrad .eq. 1 )then
            dx1dr = -2.0e0_wp*alpha1*(r-r1)*x1
            dx2dr = -2.0e0_wp*alpha2*(r-r2)*x2
            dydph = 2.0e0_wp*sin(2.0e0_wp*phi)
            dvdr  = dx1dr*y+dx2dr*y*y
            dvdph = x1*dydph+2.0e0_wp*x2*y*dydph
        endif
    end subroutine evu3to

    !> @brief Evaluate evubm for the PHOH_APRP potential.
    subroutine evubm( igrad, natm, x, y, z, ubm, dubmdx )
        ! Subroutine to evulate the Born-Mayer repulsions between non-bonded
        ! C atoms
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), parameter :: b     = 6.693126e1_wp
        real(wp), parameter :: alpha = 1.894454e0_wp
        integer             :: i, j
        integer             :: igrad, natm
        real(wp)            :: ubm
        real(wp)            :: x(natm), y(natm), z(natm)
        real(wp)            :: dubmdx(3*natm)
        ! Local variables
        integer  :: nbl, nba, nto, noop
        integer  :: ibl(2,3), iba(3,0), ito(4,0), ioop(4,0)
        real(wp) :: r14, r25, r36
        real(wp) :: dudr14, dudr25, dudr36
        real(wp) :: bm(3,3*natm)
        nbl  = 3
        nba  = 0
        nto  = 0
        noop = 0
        ibl  = reshape([ 1, 4, 2, 5, 3, 6 ], shape(ibl))
        r14  = bndlen(1, 4, x, y, z)
        r25  = bndlen(2, 5, x, y, z)
        r36  = bndlen(3, 6, x, y, z)
        ubm  = b*(exp(-alpha*r14)+exp(-alpha*r25)+exp(-alpha*r36))
        !      write(*,*) "ubm =",ubm
        if( igrad .eq. 1 )then
            dudr14 = -alpha*b*exp(-alpha*r14)
            dudr25 = -alpha*b*exp(-alpha*r25)
            dudr36 = -alpha*b*exp(-alpha*r36)
            call bmat3(nbl, ibl, nba, iba, nto, ito, noop, ioop, natm, x, y, z, bm, 3)
            do i=1, 3*natm
                dubmdx(i) = 0.e0_wp
            enddo
            do i=1, 3*natm
                dubmdx(i) = dubmdx(i)+dudr14*bm(1,i)+dudr25*bm(2,i) +dudr36*bm(3,i)
            enddo
            !        do i=1,3*natm
            !          write(*,*) dubmdX(i)
            !        enddo
        endif
    end subroutine evubm

    !      subroutine evuij(nap,natm,igrad,ir,t,dtdr,X,Y,Z,u13,u23,
    !     $                 du13dX,du23dX)
    ! Add u12
    !> @brief Evaluate evuij for the PHOH_APRP potential.
    subroutine evuij( nap, natm, igrad, ir, t, dtdr, x, y, z, u12, u13, u23, &
            & du12dx, du13dx, du23dx )
        ! Subroutine to evaluate the diabatic couplings and gradients
        ! The bohr-hartree-radian unit is used throughout the program
        ! qtc:  redundant internal coordinates
        ! stc:  non-redundant internal coordinates
        ! Ct:   transform matrix between RICs and NRICs
        !       stc=Ct*qtc
        ! bmqt: B-matrix of RICs
        ! bmst: B-matrix of NRICs
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: d13dph, d23dph, du12dr, phi, r, u12t, u13s, u13t, u23s, u23t
        integer             :: i, j, nap
        integer, parameter  :: maxdim = 100
        integer, parameter  :: ldbmt  = maxdim
        integer, parameter  :: ldbms  = maxdim
        integer, parameter  :: ldbnrt = maxdim
        integer, parameter  :: ldct   = maxdim
        integer, parameter  :: nstc   = 9
        real(wp), parameter :: pi     = 3.1415926535898e0_wp
        ! Global variables
        integer :: igrad, ir
        integer :: natm
        !      double precision u13,u23
        real(wp) :: u12, u13, u23
        real(wp) :: t(*), dtdr(*)
        real(wp) :: x(*), y(*), z(*)
        !      double precision du13dX(3*natm),du23dX(3*natm)
        real(wp) :: du12dx(3*natm), du13dx(3*natm), du23dx(3*natm)
        ! Local variables
        integer :: nqtc
        integer :: nblt, nbat, ntot, noopt
        !     data(nblt=0,nbat=0,ntot=6,noopt=6)
        integer :: iblt(2,0), ibat(3,0), itot(4,6), ioopt(4,6)
        integer :: nbls, nbas, ntos, noops
        !      data(nbls=1,nbas=0,ntos=1,noops=0)
        integer             :: ibls(2,1), ibas(3,0), itos(4,1), ioops(4,0)
        real(wp), parameter :: half   = 0.5e0_wp
        real(wp), parameter :: oossix = 0.408248290463863e0_wp
        real(wp), parameter :: oostlv = 0.288675134594813e0_wp
        real(wp)            :: qtc(maxdim), qsc(maxdim)
        real(wp)            :: stc(maxdim)
        ! KRY 20140718
        real(wp) :: d12tdr, d12sdr
        real(wp) :: d13tdr, d13sdr, du13dr, du13dph
        real(wp) :: d23tdr, d23sdr, du23dr, du23dph
        !      double precision du13t(maxdim),du23t(maxdim)
        real(wp) :: du12t(maxdim), du13t(maxdim), du23t(maxdim)
        real(wp) :: bmt(ldbmt,3*natm), bms(ldbms,3*natm)
        real(wp) :: bnrt(ldbnrt,3*natm)
        real(wp) :: ct(ldct,maxdim)
        ! Initialize Ct matrix
        nblt  = 0; nbat = 0; ntot = 6; noopt = 6; nblt = 0; nbat = 0; ntot = 6; noopt = 6
        itot  = reshape([ 2, 1, 6, 5, 1, 6, 5, 4, 6, 5, 4, 3, 5, 4, 3, 2, 4, 3, 2, 1, 3, 2, 1, 6 ], shape(itot))
        ioopt = reshape([ 7, 1, 2, 6, 12, 6, 1, 5, 11, 5, 6, 4, 10, 4, 5, 3, 9, 3, 4, 2, 8, 2, 3, 1 ], shape(ioopt))
        nbls  = 1
        nbas  = 0
        ntos  = 1
        noops = 0
        ibls  = reshape([ 7, 13 ], shape(ibls))
        itos  = reshape([ 2, 1, 7, 13 ], shape(itos))
        do i=1, 9
            do j=1, 12
                ct(i,j) = 0.0e0_wp
            enddo
        enddo
        ! Initialize Qtc vectors
        do i=1, nblt+nbat+ntot+noopt
            qtc(i) = 0.0e0_wp
        enddo
        ! Initialize Stc vectors
        do i=1, nstc
            stc(i) = 0.0e0_wp
        enddo
        ct(1,1)  = oossix
        ct(1,2)  = -oossix
        ct(1,3)  = oossix
        ct(1,4)  = -oossix
        ct(1,5)  = oossix
        ct(1,6)  = -oossix
        ct(2,1)  = -oostlv
        ct(2,2)  = 2.0e0_wp*oostlv
        ct(2,3)  = -oostlv
        ct(2,4)  = -oostlv
        ct(2,5)  = 2.0e0_wp*oostlv
        ct(2,6)  = -oostlv
        ct(3,1)  = -half
        ct(3,3)  = half
        ct(3,4)  = -half
        ct(3,6)  = half
        ct(4,7)  = 1.0e0_wp
        ct(5,8)  = 1.0e0_wp
        ct(6,9)  = 1.0e0_wp
        ct(7,10) = 1.0e0_wp
        ct(8,11) = 1.0e0_wp
        ct(9,12) = 1.0e0_wp
        ! Evaluate internal coordinates
        !Add by J. Zheng
        !      qtc(:) = 0d0
        call evintcoord(nblt, iblt, nbat, ibat, ntot, itot, noopt, ioopt, &
            & nqtc, qtc, x, y, z)
        ! Test
        !      do i=1,nqtc
        !        write(*,'(A,I2,A,F8.4)') 'qtc(',i,') = ',qtc(i)*180.d0/pi
        !      enddo
        ! End
        !Add by J. Zheng
        !      stc(:) = 0d0
        call matvect(stc, nstc, ct, ldct, qtc, nqtc)
        r   = bndlen(7, 13, x, y, z)
        phi = dihedr(2, 1, 7, 13, x, y, z)
        ! Evaluate Wilson's B-matrix if igrad = 1
        if( igrad .eq. 1 )then
            !Added by J. Zheng
            !        bmt(:,:) = 0d0
            call bmat3(nblt, iblt, nbat, ibat, ntot, itot, noopt, ioopt, natm, &
                & x, y, z, bmt, ldbmt)
            ! Evalate the Wilson's B-matrix in nonredundant internal coordinates
            call matmult(bnrt, ldbnrt, ct, ldct, nstc, bmt, ldbmt, nqtc, 3*natm)
            !Added by J. Zheng
            !        bms(:,:) = 0d0
            call bmat3(nbls, ibls, nbas, ibas, ntos, itos, noops, ioops, natm, &
                & x, y, z, bms, ldbms)
        endif
        ! Evaluate u12 KRY 20140718
        call evuijt(nap, igrad, ir, t, dtdr, nstc, stc, u12t, du12t, d12tdr, 12)
        ! Evaluate u13 and u23 from tertiary out-of-plane coordinates
        call evuijt(nap, igrad, ir, t, dtdr, nstc, stc, u13t, du13t, d13tdr, 13)
        call evuijt(nap, igrad, ir, t, dtdr, nstc, stc, u23t, du23t, d23tdr, 23)
        ! Evaluate u13 and u23 from secondary out-of-plane coordinates
        call evu13to(r, phi, u13s, d13dph, d13sdr)
        call evu23to(r, phi, u23s, d23dph, d23sdr)
        u12 = u12t                  ! KRY 20140718
        u13 = u13t+u13s
        u23 = u23t+u23s
        ! Calculated analytical Cartesian gradient if igrad = 1
        if( igrad .eq. 1 )then
            do i=1, 3*natm
                du12dx(i) = 0.0e0_wp
                du13dx(i) = 0.0e0_wp
                du23dx(i) = 0.0e0_wp
            enddo
            du12dr = d12tdr
            du13dr = d13tdr+d13sdr
            du23dr = d23tdr+d23sdr
            ! Contribution of duijdX from duijdr
            do i=1, 3*natm
                du12dx(i) = du12dx(i)+du12dr*bms(1,i)
                du13dx(i) = du13dx(i)+du13dr*bms(1,i)
                du23dx(i) = du23dx(i)+du23dr*bms(1,i)
            enddo
            ! Contribution of duijdX from dijdph
            do i=1, 3*natm
                du13dx(i) = du13dx(i)+d13dph*bms(2,i)
                du23dx(i) = du23dx(i)+d23dph*bms(2,i)
            enddo
            ! Contribution of duijdX from duijt
            do i=1, 3*natm
                do j=1, nstc
                    du12dx(i) = du12dx(i)+du12t(j)*bnrt(j,i)
                    du13dx(i) = du13dx(i)+du13t(j)*bnrt(j,i)
                    du23dx(i) = du23dx(i)+du23t(j)*bnrt(j,i)
                enddo
            enddo
        endif
    end subroutine evuij

    !> @brief Evaluate evuijt for the PHOH_APRP potential.
    subroutine evuijt( nap, igrad, ir, t, dtdr, ntc, qtc, uij, &
            & duij, duijdr, ij )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: nap
        integer  :: ij
        integer  :: igrad, ir
        integer  :: ntc
        real(wp) :: t(*), dtdr(*)
        real(wp) :: qtc(*)
        real(wp) :: uij
        real(wp) :: duijdr
        real(wp) :: duij(*)
        ! Local variables
        integer  :: i, j
        real(wp) :: u(nap)
        real(wp) :: a(nap,ntc)
        real(wp) :: uij0(nap)
        real(wp) :: uij1, uij2
        real(wp) :: duij1(ntc), duij2(ntc)
        ! Initialize uij,uij1,uij2
        uij = 0.0e0_wp; uij1 = 0.0e0_wp; uij2 = 0.0e0_wp
        ! Initialize duij,duij1,duij2
        do i=1, ntc
            duij(i)  = 0.0e0_wp
            duij1(i) = 0.0e0_wp
            duij2(i) = 0.0e0_wp
        enddo
        ! Initialize duijdr
        duijdr = 0.0e0_wp
        ! Initialize numerical gradients of diabatic couplings
        do i=1, nap
            uij0(i) = 0.0e0_wp
            do j=1, 9
                a(i,j) = 0.0e0_wp
            enddo
        enddo
        ! Assign u120 for anchor point 1 and 2
        if( ij .eq. 12 )then
            uij0(1) = -0.000616170e0_wp
            uij0(2) = -0.001064286e0_wp
        endif
        ! Assign parameters based ij
        if( ij .eq. 12 )then
            a(3,1) = 0.0e0_wp
            a(3,2) = 0.0e0_wp
            a(3,3) = 0.0e0_wp
            a(3,4) = 0.0e0_wp
            a(3,5) = 0.0e0_wp
            a(3,6) = 0.0e0_wp
            a(3,7) = 0.0e0_wp
            a(3,8) = 0.0e0_wp
            a(3,9) = 0.0e0_wp
        elseif( ij .eq. 13 )then
            a(3,1) = -0.00070e0_wp
            a(3,2) = -0.00966e0_wp
            a(3,3) = -0.00053e0_wp
            a(3,4) = 0.00193e0_wp
            a(3,5) = 0.00056e0_wp
            a(3,6) = -0.00362e0_wp
            a(3,7) = 0.00030e0_wp
            a(3,8) = 0.00377e0_wp
            a(3,9) = 0.00101e0_wp
        elseif( ij .eq. 23 )then
            a(2,1) = 0.00452e0_wp
            a(2,2) = -0.01212e0_wp
            a(2,3) = -0.00088e0_wp
            a(2,4) = 0.00532e0_wp
            a(2,5) = 0.00757e0_wp
            a(2,6) = -0.00348e0_wp
            a(2,7) = 0.00023e0_wp
            a(2,8) = -0.00209e0_wp
            a(2,9) = -0.01501e0_wp
        else
            write(*, *) "Only parameters for U12, U13 and U23 are allowed!"
            stop
        endif
        ! Check for consistence
        if( ntc .ne. 9 )then
            write(*, *) "Inconsistence in the number of out-of-plane tertiary&
                & coordinates"
            stop
        endif
        if( nap .ne. 4 )then
            write(*, *) "Inconsistence in the number of anchor structures"
            stop
        endif
        ! Now calculate coupling contributions from tertiary coordinates
        if( ir .eq. 1 )then
            uij = uij0(1)
            do i=1, ntc
                uij = uij+a(1,i)*qtc(i)
            enddo
        elseif( ir .eq. nap+1 )then
            uij = uij0(nap)
            do i=1, ntc
                uij = uij+a(nap,i)*qtc(i)
            enddo
        elseif( (ir .gt. 1) .and. (ir .lt. nap+1) )then
            uij1 = uij0(ir-1)
            uij2 = uij0(ir)
            do i=1, ntc
                uij1 = uij1+a(ir-1,i)*qtc(i)
                uij2 = uij2+a(ir,i)*qtc(i)
            enddo
            uij = uij1*t(ir-1)+uij2*t(ir)
        endif
        ! Now calculate the derivatives of couplings if igrad =1
        if( igrad .eq. 1 )then
            if( ir .eq. 1 )then
                do i=1, ntc
                    duij(i) = a(1,i)
                enddo
                duijdr = uij*dtdr(1)
            elseif( ir .eq. nap+1 )then
                do i=1, ntc
                    duij(i) = a(nap,i)
                enddo
                duijdr = uij*dtdr(nap)
            elseif( (ir .gt. 1) .and. (ir .lt. nap+1) )then
                do i=1, ntc
                    duij1(i) = a(ir-1,i)
                    duij2(i) = a(ir,i)
                    duij(i)  = duij1(i)*t(ir-1)+duij2(i)*t(ir)
                enddo
                duijdr = uij1*dtdr(ir-1)+uij2*dtdr(ir)
            endif
        endif
    end subroutine evuijt

    !> @brief Evaluate evu13to for the PHOH_APRP potential.
    subroutine evu13to( r, phi, u13to, d13dph, d13dr )
        ! Little subroutine to evaluate U13 from C-C-O-H torsion
        ! U13(r,phi)=f1(r)*sin(phi)+f3(r)*sin(phi)**3+f5(r)*sin(phi)**5
        ! where fi(r) are expanded with a linear combination of Gaussians
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)           :: df1rdr, df3rdr, df5rdr
        integer, parameter :: ngf = 3
        real(wp)           :: a1(ngf), a3(ngf), a5(ngf)
        real(wp)           :: alpha1(ngf), alpha3(ngf), alpha5(ngf)
        real(wp)           :: r1(ngf), r3(ngf), r5(ngf)
        real(wp)           :: r, phi, u13to, d13dph, d13dr
        ! Local variables
        integer  :: i
        real(wp) :: sphi, cphi
        real(wp) :: sphi1, sphi2, sphi3, sphi4, sphi5
        real(wp) :: f1r, f3r, f5r
        real(wp) :: df1dr, df3dr, df5dr
        real(wp) :: g1r, g3r, g5r
        a1     = reshape([ 0.068885e0_wp, 0.034221e0_wp, 0.001544e0_wp ], shape(a1))
        a3     = reshape([ -0.082676e0_wp, -0.021770e0_wp, 0.0e0_wp ], shape(a3))
        a5     = reshape([ 0.053915e0_wp, -0.008539e0_wp, .0e0_wp ], shape(a5))
        alpha1 = reshape([ 0.2956970e1_wp, 0.4695508e1_wp, 0.1167016e1_wp ], shape(alpha1))
        alpha3 = reshape([ 0.5029657e1_wp, 0.8744002e1_wp, 0.0e0_wp ], shape(alpha3))
        alpha5 = reshape([ 0.7263172e1_wp, 0.27942459e2_wp, 0.0e0_wp ], shape(alpha5))
        r1     = reshape([ 0.3162075e1_wp, 0.2255351e1_wp, 0.4480427e1_wp ], shape(r1))
        r3     = reshape([ 0.3286539e1_wp, 0.2060405e1_wp, 0.0e0_wp ], shape(r3))
        r5     = reshape([ 0.3275937e1_wp, 0.2498052e1_wp, 0.0e0_wp ], shape(r5))
        u13to  = 0.0e0_wp
        d13dph = 0.0e0_wp
        d13dr  = 0.0e0_wp
        f1r    = 0.0e0_wp
        f3r    = 0.0e0_wp
        f5r    = 0.0e0_wp
        df1rdr = 0.0e0_wp
        df3rdr = 0.0e0_wp
        df5rdr = 0.0e0_wp
        do i=1, ngf
            g1r    = a1(i)*exp(-alpha1(i)*(r-r1(i))*(r-r1(i)))
            f1r    = f1r+g1r
            df1rdr = df1rdr-2.0e0_wp*alpha1(i)*(r-r1(i))*g1r
            g3r    = a3(i)*exp(-alpha3(i)*(r-r3(i))*(r-r3(i)))
            f3r    = f3r+g3r
            df3rdr = df3rdr-2.0e0_wp*alpha3(i)*(r-r3(i))*g3r
            g5r    = a5(i)*exp(-alpha5(i)*(r-r5(i))*(r-r5(i)))
            f5r    = f5r+g5r
            df5rdr = df5rdr-2.0e0_wp*alpha5(i)*(r-r5(i))*g5r
        enddo
        sphi   = sin(phi)
        cphi   = cos(phi)
        sphi1  = sphi
        sphi2  = sphi*sphi1
        sphi3  = sphi*sphi2
        sphi4  = sphi*sphi3
        sphi5  = sphi*sphi4
        u13to  = f1r*sphi1+f3r*sphi3+f5r*sphi5
        d13dph = f1r*cphi+3.0e0_wp*f3r*sphi2*cphi+5.0e0_wp*f5r*sphi4*cphi
        d13dr  = df1rdr*sphi1+df3rdr*sphi3+df5rdr*sphi5
    end subroutine evu13to

    !> @brief Evaluate evu23to for the PHOH_APRP potential.
    subroutine evu23to( r, phi, u23to, d23dph, d23dr )
        ! Little subroutine to evaluate U23 from C-C-O-H torsion
        ! U23(r,phi)=f2(r)*sin(2*phi)+f4(r)*sin(4*phi)+f6(r)*sin(6*phi)
        ! where fi(r) are expanded with a linear combination of Gaussians
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: c2phi, c4phi, c6phi, df2rdr, df4rdr, df6rdr, s2phi, s4phi, s6phi
        integer, parameter  :: ngf   = 1
        real(wp), parameter :: scale = 1.0e0_wp
        real(wp)            :: a2(ngf), a4(ngf), a6(ngf)
        real(wp)            :: alpha2(ngf), alpha4(ngf), alpha6(ngf)
        real(wp)            :: r2(ngf), r4(ngf), r6(ngf)
        real(wp)            :: r, phi, u23to, d23dph, d23dr
        ! Local variables
        integer  :: i
        real(wp) :: sphi, cphi
        real(wp) :: sphi1, sphi2, sphi3, sphi4, sphi5
        real(wp) :: f2r, f4r, f6r
        real(wp) :: df2dr, df4dr, df6dr
        real(wp) :: g2r, g4r, g6r
        a2     = 0.016924e0_wp
        a4     = 0.007145e0_wp
        a6     = 0.003113e0_wp
        alpha2 = 0.10841076e2_wp
        alpha4 = 0.12574937e2_wp
        alpha6 = 0.25329325e2_wp
        r2     = 0.2151852e1_wp
        r4     = 0.2197436e1_wp
        r6     = 0.2285492e1_wp
        u23to  = 0.0e0_wp
        d23dph = 0.0e0_wp
        d23dr  = 0.0e0_wp
        f2r    = 0.0e0_wp
        f4r    = 0.0e0_wp
        f6r    = 0.0e0_wp
        df2rdr = 0.0e0_wp
        df4rdr = 0.0e0_wp
        df6rdr = 0.0e0_wp
        do i=1, ngf
            g2r    = a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
            f2r    = f2r+g2r
            df2rdr = df2rdr-2.0e0_wp*alpha2(i)*(r-r2(i))*g2r
            g4r    = a4(i)*exp(-alpha4(i)*(r-r4(i))*(r-r4(i)))
            f4r    = f4r+g4r
            df4rdr = df4rdr-2.0e0_wp*alpha4(i)*(r-r4(i))*g4r
            g6r    = a6(i)*exp(-alpha6(i)*(r-r6(i))*(r-r6(i)))
            f6r    = f6r+g6r
            df6rdr = df6rdr-2.0e0_wp*alpha6(i)*(r-r6(i))*g6r
        enddo
        s2phi = sin(2.0e0_wp*phi)
        c2phi = cos(2.0e0_wp*phi)
        s4phi = sin(4.0e0_wp*phi)
        c4phi = cos(4.0e0_wp*phi)
        s6phi = sin(6.0e0_wp*phi)
        c6phi = cos(6.0e0_wp*phi)
        !      u23to=f2r*s2phi+f4r*s4phi+f6r*s6phi
        !      d23dph=2.0d0*f2r*c2phi+4.0d0*f4r*c4phi+6.0d0*f6r*c6phi
        !      d23dr=df2rdr*s2phi+df4rdr*s4phi+df6rdr*s6phi
        u23to  = scale*(f2r*s2phi+f4r*s4phi+f6r*s6phi)
        d23dph = scale*(2.0e0_wp*f2r*c2phi+4.0e0_wp*f4r*c4phi+6.0e0_wp*f6r*c6phi)
        d23dr  = scale*(df2rdr*s2phi+df4rdr*s4phi+df6rdr*s6phi)
    end subroutine evu23to

    !> @brief Evaluate evintcoord for the PHOH_APRP potential.
    subroutine evintcoord( nbl, ibl, nba, iba, nto, ito, noop, ioop, &
            & nint, ric, x, y, z )
        ! Evaluate the the interal coordinates for given connections
        ! and Cartesian coordinates
        ! nbl:          # of bond lengths
        ! nba:          # of bond angles
        ! nto:          # of torsions
        ! noop:         # of out-of-plane bendings
        ! nint:         # of internal coordinates
        ! ibl:          index of bond stretchings
        ! iba:          index of bond angles
        ! ito:          index of torsions
        ! ioop:         index of out-of-plane bendings
        ! RIC:          redundant internal coordinates
        ! X,Y,Z:        Cartesian coordinates
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter  :: maxint = 1000
        real(wp), parameter :: btoa   = 0.5291772192e0_wp
        real(wp), parameter :: pi     = 3.1415926535898e0_wp
        integer             :: i
        integer             :: nbl, nba, noop, nto, nint
        integer             :: ibl(2,*), iba(3,*), ito(4,*), ioop(4,*)
        real(wp)            :: x(*), y(*), z(*)
        real(wp)            :: rab(maxint), theta(maxint), tors(maxint), oopb(maxint)
        real(wp)            :: ric(*)
        ! Initialize RIC
        do i=1, nbl+nba+nto+noop
            ric(i) = 0.0e0_wp
        enddo
        ! Initialize rab,theta,tors, and oopb
        do i=1, maxint
            rab(i)   = 0.0e0_wp
            theta(i) = 0.0e0_wp
            tors(i)  = 0.0e0_wp
            oopb(i)  = 0.0e0_wp
        enddo
        ! Evaluate bond lengths
        do i=1, nbl
            rab(i) = bndlen(ibl(1,i), ibl(2,i), x, y, z)
        enddo
        ! Evaluate bend angles
        do i=1, nba
            theta(i) = bndang(iba(1,i), iba(2,i), iba(3,i), x, y, z)
        enddo
        ! Evaluate torsions
        do i=1, nto
            tors(i) = dihedr(ito(1,i), ito(2,i), ito(3,i), ito(4,i), x, y, z)
        enddo
        ! Evaluate out-of-plane bending angles
        do i=1, noop
            oopb(i) = oopba1(ioop(1,i), ioop(2,i), ioop(3,i), ioop(4,i), x, y, z)
        enddo
        nint = 0
        do i=1, nbl
            nint      = nint+1
            ric(nint) = rab(i)
        enddo
        do i=1, nba
            nint      = nint+1
            ric(nint) = theta(i)
        enddo
        do i=1, nto
            nint      = nint+1
            ric(nint) = tors(i)
        enddo
        do i=1, noop
            nint      = nint+1
            ric(nint) = oopb(i)
        enddo
    end subroutine evintcoord

    !> @brief Evaluate matvect for the PHOH_APRP potential.
    subroutine matvect( y, m, a, lda, x, n )
        ! Subroutien to calculate matrix A multiply vector X
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j
        integer  :: lda, m, n
        real(wp) :: x(*), y(*)
        real(wp) :: a(lda,*)
        ! Initialize Y
        do i=1, m
            y(i) = 0.0e0_wp
        enddo
        do i=1, m
            do j=1, n
                y(i) = y(i)+a(i,j)*x(j)
            enddo
        enddo
    end subroutine matvect

    !> @brief Evaluate bmat3 for the PHOH_APRP potential.
    subroutine bmat3( nbl, ibl, nba, iba, nto, ito, noop, ioop, natm, &
            & x, y, z, bm, ldbm )
        ! Subroutine to evaluate the Wilson B-matrix
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: ai, aijk, ajkl, ak, al, ci, cijk, cjkl, ck, cl, ctheta, rij, rji, rjk, rjl, rkl, si, sijk, sjkl, sk, sl, theta, ttheta
        ! Global variables
        integer  :: ldbm
        integer  :: nbl, nba, noop, nto, natm
        integer  :: ibl(2,nbl), iba(3,nba), ito(4,nto), ioop(4,noop)
        real(wp) :: x(*), y(*), z(*)
        real(wp) :: bm(ldbm,*)
        ! Local variables
        integer  :: i, j, k, l, ij, ijk, ijkl, irow
        real(wp) :: eij(3), ejk(3), ekl(3), eji(3), ejl(3), u(3), v(3)
        real(wp) :: w(3)
        ! Initialize B matrix
        do i=1, ldbm
            do j=1, 3*natm
                bm(i,j) = 0.0e0_wp
            enddo
        enddo
        irow = 0
        ! Evaluate the B matrix contributed from bond stretchings
        do ij=1, nbl
            irow           = irow+1
            i              = ibl(1,ij)
            j              = ibl(2,ij)
            rij            = bndlen(i, j, x, y, z)
            eij(1)         = (x(j)-x(i))/rij
            eij(2)         = (y(j)-y(i))/rij
            eij(3)         = (z(j)-z(i))/rij
            bm(irow,3*i-2) = -eij(1)
            bm(irow,3*i-1) = -eij(2)
            bm(irow,3*i)   = -eij(3)
            bm(irow,3*j-2) = eij(1)
            bm(irow,3*j-1) = eij(2)
            bm(irow,3*j)   = eij(3)
        enddo
        ! Evaluate the B matrix contributed from bendings
        do ijk=1, nba
            irow           = irow+1
            i              = iba(1,ijk)
            j              = iba(2,ijk)
            k              = iba(3,ijk)
            rji            = bndlen(j, i, x, y, z)
            rjk            = bndlen(j, k, x, y, z)
            aijk           = bndang(i, j, k, x, y, z)
            cijk           = cos(aijk)
            sijk           = sin(aijk)
            eji(1)         = (x(i)-x(j))/rji
            eji(2)         = (y(i)-y(j))/rji
            eji(3)         = (z(i)-z(j))/rji
            ejk(1)         = (x(k)-x(j))/rjk
            ejk(2)         = (y(k)-y(j))/rjk
            ejk(3)         = (z(k)-z(j))/rjk
            bm(irow,3*i-2) = (cijk*eji(1)-ejk(1))/(rji*sijk)
            bm(irow,3*i-1) = (cijk*eji(2)-ejk(2))/(rji*sijk)
            bm(irow,3*i)   = (cijk*eji(3)-ejk(3))/(rji*sijk)
            bm(irow,3*k-2) = (cijk*ejk(1)-eji(1))/(rjk*sijk)
            bm(irow,3*k-1) = (cijk*ejk(2)-eji(2))/(rjk*sijk)
            bm(irow,3*k)   = (cijk*ejk(3)-eji(3))/(rjk*sijk)
            bm(irow,3*j-2) = -bm(irow,3*i-2)-bm(irow,3*k-2)
            bm(irow,3*j-1) = -bm(irow,3*i-1)-bm(irow,3*k-1)
            bm(irow,3*j)   = -bm(irow,3*i)  -bm(irow,3*k)
        enddo
        ! Evaluate the B matrix contributed from torsions
        ! The calculations follow I. H. Williams, J. Mol. Spectros. 66, 288 (197
        do ijkl=1, nto
            irow   = irow+1
            i      = ito(1,ijkl)
            j      = ito(2,ijkl)
            k      = ito(3,ijkl)
            l      = ito(4,ijkl)
            rij    = bndlen(i, j, x, y, z)
            rjk    = bndlen(j, k, x, y, z)
            rkl    = bndlen(k, l, x, y, z)
            aijk   = bndang(i, j, k, x, y, z)
            ajkl   = bndang(j, k, l, x, y, z)
            cijk   = cos(aijk)
            sijk   = sin(aijk)
            cjkl   = cos(ajkl)
            sjkl   = sin(ajkl)
            eij(1) = (x(j)-x(i))/rij
            eij(2) = (y(j)-y(i))/rij
            eij(3) = (z(j)-z(i))/rij
            ejk(1) = (x(k)-x(j))/rjk
            ejk(2) = (y(k)-y(j))/rjk
            ejk(3) = (z(k)-z(j))/rjk
            ekl(1) = (x(l)-x(k))/rkl
            ekl(2) = (y(l)-y(k))/rkl
            ekl(3) = (z(l)-z(k))/rkl
            ! Caculate u=eij*ejk, v=ejk*ekl
            call xprod1(eij, ejk, u)
            call xprod1(ejk, ekl, v)
            bm(irow,3*i-2) = -u(1)/(rij*sijk*sijk)
            bm(irow,3*i-1) = -u(2)/(rij*sijk*sijk)
            bm(irow,3*i)   = -u(3)/(rij*sijk*sijk)
            bm(irow,3*j-2) = u(1)*(rjk-rij*cijk)/(rij*rjk*sijk*sijk) -v(1)*cjkl/(rjk*sjkl*sjkl)
            bm(irow,3*j-1) = u(2)*(rjk-rij*cijk)/(rij*rjk*sijk*sijk) -v(2)*cjkl/(rjk*sjkl*sjkl)
            bm(irow,3*j)   = u(3)*(rjk-rij*cijk)/(rij*rjk*sijk*sijk) -v(3)*cjkl/(rjk*sjkl*sjkl)
            bm(irow,3*k-2) = -v(1)*(rjk-rkl*cjkl)/(rjk*rkl*sjkl*sjkl) +u(1)*cijk/(rjk*sijk*sijk)
            bm(irow,3*k-1) = -v(2)*(rjk-rkl*cjkl)/(rjk*rkl*sjkl*sjkl) +u(2)*cijk/(rjk*sijk*sijk)
            bm(irow,3*k)   = -v(3)*(rjk-rkl*cjkl)/(rjk*rkl*sjkl*sjkl) +u(3)*cijk/(rjk*sijk*sijk)
            bm(irow,3*l-2) = v(1)/(rkl*sjkl*sjkl)
            bm(irow,3*l-1) = v(2)/(rkl*sjkl*sjkl)
            bm(irow,3*l)   = v(3)/(rkl*sjkl*sjkl)
        enddo
        ! Evaluate the B matrix contributed from out-of-plane bending
        ! The calculation follows D. F. Mclntosh and K. H. Michaelian,
        ! Can. J. Spectros. 24, 35-40 (1979).
        do ijkl=1, noop
            irow   = irow+1
            i      = ioop(1,ijkl)
            j      = ioop(2,ijkl)
            k      = ioop(3,ijkl)
            l      = ioop(4,ijkl)
            rji    = bndlen(j, i, x, y, z)
            rjk    = bndlen(j, k, x, y, z)
            rjl    = bndlen(j, l, x, y, z)
            ai     = bndang(k, j, l, x, y, z)
            ak     = bndang(i, j, l, x, y, z)
            al     = bndang(i, j, k, x, y, z)
            theta  = oopba1(i, j, k, l, x, y, z)
            ci     = cos(ai)
            si     = sin(ai)
            ck     = cos(ak)
            sk     = sin(ak)
            cl     = cos(al)
            sl     = sin(al)
            ctheta = cos(theta)
            ttheta = tan(theta)
            eji(1) = (x(i)-x(j))/rji
            eji(2) = (y(i)-y(j))/rji
            eji(3) = (z(i)-z(j))/rji
            ejk(1) = (x(k)-x(j))/rjk
            ejk(2) = (y(k)-y(j))/rjk
            ejk(3) = (z(k)-z(j))/rjk
            ejl(1) = (x(l)-x(j))/rjl
            ejl(2) = (y(l)-y(j))/rjl
            ejl(3) = (z(l)-z(j))/rjl
            ! Caculate u=ejk*ejl
            call xprod1(ejk, ejl, u)
            bm(irow,3*i-2) = (1.0e0_wp/rji)*(u(1)/(ctheta*si)-ttheta*eji(1))
            bm(irow,3*i-1) = (1.0e0_wp/rji)*(u(2)/(ctheta*si)-ttheta*eji(2))
            bm(irow,3*i)   = (1.0e0_wp/rji)*(u(3)/(ctheta*si)-ttheta*eji(3))
            bm(irow,3*k-2) = (1.0e0_wp/rjk)*u(1)*(ci*ck-cl)/(ctheta*si*si*si)
            bm(irow,3*k-1) = (1.0e0_wp/rjk)*u(2)*(ci*ck-cl)/(ctheta*si*si*si)
            bm(irow,3*k)   = (1.0e0_wp/rjk)*u(3)*(ci*ck-cl)/(ctheta*si*si*si)
            bm(irow,3*l-2) = (1.0e0_wp/rjl)*u(1)*(ci*cl-ck)/(ctheta*si*si*si)
            bm(irow,3*l-1) = (1.0e0_wp/rjl)*u(2)*(ci*cl-ck)/(ctheta*si*si*si)
            bm(irow,3*l)   = (1.0e0_wp/rjl)*u(3)*(ci*cl-ck)/(ctheta*si*si*si)
            bm(irow,3*j-2) = -bm(irow,3*i-2)-bm(irow,3*k-2)-bm(irow,3*l-2)
            bm(irow,3*j-1) = -bm(irow,3*i-1)-bm(irow,3*k-1)-bm(irow,3*l-1)
            bm(irow,3*j)   = -bm(irow,3*i)  -bm(irow,3*k)  -bm(irow,3*l)
        enddo
    end subroutine bmat3

    !> @brief Evaluate matmult for the PHOH_APRP potential.
    subroutine matmult( c, ldc, a, lda, m, b, ldb, n, l )
        ! Global variables
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: lda, ldb, ldc
        integer  :: m, n, l
        real(wp) :: a(lda,*), b(ldb,*), c(ldc,*)
        ! Local variables
        integer :: i, j, k
        do i=1, m
            do j=1, l
                c(i,j) = 0.0e0_wp
            enddo
        enddo
        do i=1, m
            do j=1, l
                do k=1, n
                    c(i,j) = c(i,j)+a(i,k)*b(k,j)
                enddo
            enddo
        enddo
    end subroutine matmult

    !> @brief Evaluate bndlen for the PHOH_APRP potential.
    real(wp) function bndlen( i, j, x, y, z )
        ! Function to evaluate the bond length between atom i and j
        ! Input: i,j (indices of 2 atoms)
        use chempotpy_kinds, only: wp
        implicit none
        integer             :: i, j
        real(wp), parameter :: tiny = 1.0e-13_wp
        real(wp)            :: x(*), y(*), z(*)
        bndlen = sqrt( (x(i)-x(j))*(x(i)-x(j)) + (y(i)-y(j))*(y(i)-y(j)) + (z(i)-z(j))*(z(i)-z(j)))
        if( bndlen .le. tiny )then
            write(*, '('' Atom'',I3,'' and'',I3,'' are too close to each other!'')') i, j
            !DB  KRY 20140718
            write(*, '(I3,3F18.12)') i, x(i), y(i), z(i)
            write(*, '(I3,3F18.12)') j, x(j), y(j), z(j)
            stop
        endif
    end function bndlen

    !> @brief Evaluate bndang for the PHOH_APRP potential.
    real(wp) function bndang( i, j, k, x, y, z )
        ! Evaluate the bond angle between atoms i, j, and k
        ! Input: i,j,k (indices of 3 atoms, i-j and j-k bonds)
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a
        integer  :: i, j, k
        real(wp) :: x(*), y(*), z(*)
        real(wp) :: eji(3), ejk(3), rji, rjk
        rji    = bndlen(i, j, x, y, z)
        rjk    = bndlen(j, k, x, y, z)
        eji(1) = (x(i)-x(j))/rji
        eji(2) = (y(i)-y(j))/rji
        eji(3) = (z(i)-z(j))/rji
        ejk(1) = (x(k)-x(j))/rjk
        ejk(2) = (y(k)-y(j))/rjk
        ejk(3) = (z(k)-z(j))/rjk
        call sprod(3, eji, ejk, a)
        bndang = acos(a)
    end function bndang

    !> @brief Evaluate dihedr for the PHOH_APRP potential.
    real(wp) function dihedr( i, j, k, l, x, y, z )
        ! Evaluate the dihedral angle between atoms i, j, k, and l
        ! Input: i,j,k,l (indices of 4 atoms)
        !        i--j--k--l
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: a
        integer             :: i, j, k, l
        real(wp), parameter :: tiny = 1.0e-13_wp
        real(wp), parameter :: pi   = 3.1415926535898e0_wp
        real(wp)            :: x(*), y(*), z(*)
        real(wp)            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
        real(wp)            :: rij, rjk, rkl, aijk, ajkl, wejk
        !   eij
        ! i-----j
        !   aijk \ ejk
        !         \
        !          \ ajkl
        !           k----l
        !             kl
        rij  = bndlen(i, j, x, y, z)
        rjk  = bndlen(j, k, x, y, z)
        rkl  = bndlen(k, l, x, y, z)
        aijk = bndang(i, j, k, x, y, z)
        ajkl = bndang(j, k, l, x, y, z)
        ! Check whether i-j-k or j-k-l are collinear
        if( (abs(aijk) .le. tiny) .or. (abs(aijk-pi) .lt. tiny) )then
            write(*, '(''Torsion angle not defined since atoms'',3I3,'' collinear'')') i, j, k
        endif
        if( (abs(ajkl) .le. tiny) .or. (abs(ajkl-pi) .lt. tiny) )then
            write(*, '(''Torsion angle not defined since atoms'',3I3,'' collinear'')') j, k, l
        endif
        ! Calculate eij
        eij(1) = (x(j)-x(i))/rij
        eij(2) = (y(j)-y(i))/rij
        eij(3) = (z(j)-z(i))/rij
        ! Calculate ejk
        ejk(1) = (x(k)-x(j))/rjk
        ejk(2) = (y(k)-y(j))/rjk
        ejk(3) = (z(k)-z(j))/rjk
        ! Calculate ekl
        ekl(1) = (x(l)-x(k))/rkl
        ekl(2) = (y(l)-y(k))/rkl
        ekl(3) = (z(l)-z(k))/rkl
        ! u = eij*ejk normal vector in plane i-j-k
        ! v = ejk*ekl Normal vector in plane j-k-l
        ! w = u*v used to determine the sign of dihedral angle
        call xprod1(eij, ejk, u)
        call xprod1(ejk, ekl, v)
        call xprod1(u, v, w)
        call sprod(3, u, v, a)
        a      = a/(sin(aijk)*sin(ajkl))
        a      = max(a, -1.0e0_wp)
        a      = min(a, 1.0e0_wp)
        dihedr = acos(a)
        ! Determine the sign of dihedr based on u,v,ejk
        ! dihedr >= 0, u,v,ejk right-hand (w*ejk >=0)
        ! dihedr < 0,  u,v,ejk left-hand (w*ejk < 0)
        wejk = w(1)*ejk(1)+w(2)*ejk(2)+w(3)*ejk(3)
        if( wejk .lt. 0.0e0_wp )then
            dihedr = -dihedr
        endif
        ! Setting the range of dihedral angle between -90.0 degree to 270 degree
        ! to avoid discontinuity around 180 degrees
        !      if (dihedr .lt. -0.5d0*pi) then
        !        dihedr=2.0d0*pi+dihedr
        !      endif
    end function dihedr

    !> @brief Evaluate oopba1 for the PHOH_APRP potential.
    real(wp) function oopba1( i, j, k, l, x, y, z )
        ! Evaluate the out-of-plane bending angle between atoms i, j, k, and l
        ! Input: i,j,k,l (indices of 4 atoms)
        !             l
        !            /
        !     i-----j
        !            \
        !             k
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j, k, l
        real(wp) :: tiny
        real(wp) :: x(*), y(*), z(*)
        real(wp) :: eji(3), ejk(3), ejl(3), u(3)
        real(wp) :: v(3), w(3)
        real(wp) :: a1, a2
        real(wp) :: a, ai, si, rji, rjk, rjl
        !               l
        !              /
        !          ak /
        !     i------j ai != 180 deg
        !          al \
        !              \
        !               k
        tiny   = 1.0e-30_wp
        rji    = bndlen(j, i, x, y, z)
        eji(1) = (x(i)-x(j))/rji
        eji(2) = (y(i)-y(j))/rji
        eji(3) = (z(i)-z(j))/rji
        rjk    = bndlen(j, k, x, y, z)
        ejk(1) = (x(k)-x(j))/rjk
        ejk(2) = (y(k)-y(j))/rjk
        ejk(3) = (z(k)-z(j))/rjk
        rjl    = bndlen(j, l, x, y, z)
        ejl(1) = (x(l)-x(j))/rjl
        ejl(2) = (y(l)-y(j))/rjl
        ejl(3) = (z(l)-z(j))/rjl
        ai     = bndang(k, j, l, x, y, z)
        si     = sin(ai)
        call xprod1(ejk, ejl, u)
        ! Check if whether ai = 180
        if( u(1)**2+u(2)**2+u(3)**2 .le. tiny )then
            write(*, *) "Looks like k-j-l collinear."
            write(*, *) "Out-of-plane bending angle is not defined."
            stop
        endif
        call sprod(3, u, eji, a)
        a      = a/si
        a      = max(a, -1.0e0_wp)
        a      = min(a, 1.0e0_wp)
        oopba1 = asin(a)
        ! Check whether needs to go beyond -pi/2 to pi/2
        !      write(*,*) "before correction, oopba1 = ",oopba1
        !      do n=1,3
        !        v(n)=ejk(n)+ejl(n)
        !      enddo
        !      call sprod(3,eji,v,b)
        !      if ((oopba1 .lt. 0.0d0) .and. (b .gt. 0.0d0)) then
        !        oopba1=2.0d0*asin(-1.0d0)-oopba1
        !      endif
        !      if ((oopba1 .ge. 0.0d0) .and. (b .gt. 0.0d0)) then
        !        oopba1=2.0d0*asin(1.0d0)-oopba1
        !      endif
    end function oopba1

    !> @brief Evaluate xprod1 for the PHOH_APRP potential.
    subroutine xprod1( x1, x2, x3 )
        ! Calculate cross-product x1 and x2 and return x3
        ! Input: X1,X2
        ! Output: X3
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x1(3), x2(3), x3(3)
        x3(1) = x1(2)*x2(3)-x1(3)*x2(2)
        x3(2) = x1(3)*x2(1)-x1(1)*x2(3)
        x3(3) = x1(1)*x2(2)-x1(2)*x2(1)
    end subroutine xprod1

    !> @brief Evaluate sprod for the PHOH_APRP potential.
    subroutine sprod( n, a, b, c )
        ! Calculate scalar-product a and b and return c
        ! Input: a(n),b(n)
        ! Output: c
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i
        integer  :: n
        real(wp) :: a(*), b(*), c
        c = 0.0e0_wp
        do i=1, n
            c = c+a(i)*b(i)
        enddo
    end subroutine sprod
end module chempotpy_phoh_aprp_pes

!> @brief Expose PHOH_APRP through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_phoh_aprp_pes, only: evaluate
    use chempotpy_kinds,         only: wp
    implicit none
    integer, parameter    :: natoms  = 13
    integer, parameter    :: nstates = 3
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
