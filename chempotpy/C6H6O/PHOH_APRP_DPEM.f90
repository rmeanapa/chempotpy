!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!   System:                     C6H6O
!   Functional form:            Anchor-Points Reactive Potential (APRP)
!                               (3x3 diabatic fit)
!   Common name:                PhOH(Coupled PESs)
!   Number of derivatives:      1
!   Number of bodies:           13
!   Number of electronic surfaces: 3
!   Interface: Section-2
!
!   Reference: Ke R. Yang, Xuefei Xu, Jingjing Zheng,
!   and Donald G. Truhlar Chem. Sci. (to be published).
!
!   Notes:
!    1. This routine calculates three coupled singlet PESs for
!       nonadiabatic photodissociation of PhOH to produce
!       PhO + H.
!    2. The final version of this potential was added to POTLIB
!       on Aug. 3, 2014. A preliminary version in the library before
!       that is not correct.
!
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
!
!
!   Units:      bohr for length
!               hartree for energy
!
!   Input:
!   igrad = 0           Energy only calculation
!         = 1           Energy + Analytic gradient
!   repflag             flag to indicate wheter to use diabatic or
!                       adiabatic representation (used by ANT program)
!   xx(1,i)             X coordinate of atom i
!   xx(2,i)             Y coordinate of atom i
!   xx(3,i)             Z coordinate of atom i
!
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
!
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
!
!   Note:       LAPACK library is needed to diagonize diabatic potential
!               matrix (UU) to yield adiabatic potenitials VV(3)
!
!***********************************************************************

!> @brief Evaluate dpem for the PHOH_APRP_DPEM potential.
subroutine dpem( x, igrad, u, ug )

    implicit none

    ! number of electronic state
    integer, parameter            :: nstates = 3
    integer, parameter            :: natoms  = 13
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: u(nstates,nstates)
    double precision, intent(out) :: ug(nstates,nstates,natoms,3)

    double precision :: xx(3,natoms), uu(nstates,nstates)
    double precision :: guu(3,natoms,nstates,nstates)
    double precision :: vv(nstates), gvv(3,natoms,nstates)
    double precision :: dvec(3,natoms,nstates,nstates)
    double precision :: cc(nstates,nstates)
    integer          :: iatom, idir, j, istate, jstate
    !initialize
    u  = 0.d0
    ug = 0.d0

    do iatom=1,natoms
        do idir=1,3
            xx(idir,iatom) = x(iatom,idir)/0.529177211
        enddo
    enddo

    call pot(igrad,xx,uu,guu,vv,gvv,dvec,cc,0)

    if (igrad==0) then
        do istate=1,nstates
            do jstate=1,nstates
                u(istate,jstate) = uu(istate,jstate)*27.211386
            enddo
        enddo
    elseif (igrad==1) then
        do istate=1,nstates
            do jstate=1,nstates
                u(istate,jstate) = uu(istate,jstate)*27.211386
            enddo
        enddo
        do istate=1,nstates
            do jstate=1,nstates
                do iatom=1,natoms
                    do idir=1,3
                        ug(istate,jstate,iatom,idir) = guu(idir,iatom,istate,jstate)*51.
                    enddo
                enddo
            enddo
        enddo
    endif

end subroutine dpem

!***********************************************************************
!> @brief Evaluate pot for the PHOH_APRP_DPEM potential.
subroutine pot( igrad, xx, uu, guu, vv, gvv, dvec, cc, repflag )
    !***********************************************************************
    !
    !***********************************************************************
    implicit none

    integer            :: i,j,k,l,m
    integer            :: igrad, repflag
    integer, parameter :: natm   = 13
    integer, parameter :: nsurft = 3

    double precision :: uu(3,3),vv(3)
    double precision :: guu(3,natm,3,3),uup(3,3),vvp(3),vvm(3)
    double precision :: guup(3,natm,3,3),gvv(3,natm,nsurft)
    double precision :: guum(3,natm,3,3)
    double precision :: xx(3,natm),x(natm),y(natm),z(natm)
    double precision :: dvec(3,natm,nsurft,nsurft)
    !     double precision cc(nsurft,nsurft),ccp(nsurft,nsurft)
    double precision            :: cc(3,3),ccp(3,3)
    double precision            :: tmp
    double precision, parameter :: stp = 1.d-4

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
    call poten(1,x,y,z,uu,guu,vv,gvv,cc)
    ! Compute the nonadiabatic coupling vectors DVEC.
    !      call getdvec3(natm,uu,guu,vv,dvec,cc)

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
!> @brief Evaluate getdvec3 for the PHOH_APRP_DPEM potential.
subroutine getdvec3( nclu, pemd, gpemd, pema, dvec, cc )
    implicit none

    integer, parameter :: nsurft = 3
    double precision   :: pema(nsurft),pemd(nsurft,nsurft)
    double precision   :: gpemd(3,nclu,nsurft,nsurft)
    double precision   :: dvec(3,nclu,nsurft,nsurft),cc(nsurft,nsurft)
    integer            :: i1,i2,i,j,k,l,nclu

    !     Program to compute nonadiabatic coupling vectors from the adiabati
    !     the gradients of diabatic potential matrix elements for a multisur
    !     The formula is taken from ANT07 and is described in M.D. Hack et a
    !
    ! compute d
    dvec(:,:,:,:) = 0d0
    do i1 = 1,3
        do i2 = 1,nclu
            do i = 2,nsurft
                do j = 1, i-1
                    !            dvec(i1,i2,i,j) = 0.d0
                    do k = 1, nsurft
                        do l = 1, nsurft
                            dvec(i1,i2,i,j) = dvec(i1,i2,i,j) + cc(k,i)*cc(l,j)*gpemd(i1,i2,k,l)
                        enddo
                    enddo
                    if ((pema(j) - pema(i)) .ne. 0.0d0) then
                        dvec(i1,i2,i,j) = dvec(i1,i2,i,j)/(pema(j)-pema(i))
                    else
                        dvec(i1,i2,i,j) = 0.0d0
                    endif
                    !              write(6,*) "dvec ",i1,i2,i,j,dvec(i1,i1,i,j)
                    dvec(i1,i2,j,i) = -dvec(i1,i2,i,j)
                enddo
            enddo
        enddo
    enddo
    !       stop
end subroutine getdvec3

!> @brief Evaluate poten for the PHOH_APRP_DPEM potential.
subroutine poten( igrad, x, y, z, uu, guu, vv, gvv, cc )
    !***********************************************************************
    ! Little subroutine to call evuii and evuij
    !***********************************************************************

    implicit none

    double precision :: bndlen, u12
    integer          :: i, j

    integer, parameter :: natm = 13
    integer, parameter :: nap  = 4
    double precision   :: rap(nap)
    data rap/1.93206106d0,2.49443845d0,4.27078097d0,9.44863047d0/

    integer, intent(in)          :: igrad
    double precision, intent(in) :: x(natm),y(natm),z(natm)
    integer                      :: ir
    double precision             :: roh
    double precision             :: t(nap),dtdr(nap)
    double precision             :: u1,u2,u3
    ! Add by KRY
    double precision :: ubm
    double precision :: dubmdx(3*natm)
    !
    double precision :: du1dx(3*natm),du2dx(3*natm),du3dx(3*natm)
    double precision :: u13,u23
    !      double precision du13dX(3*natm),du23dX(3*natm)
    double precision :: du12dx(3*natm),du13dx(3*natm),du23dx(3*natm)
    ! add by J. Zheng
    integer                       :: lwork,info
    double precision, intent(out) :: uu(3,3),vv(3), &
        & guu(3,natm,3,3)
    double precision              :: cc(3,3)
    double precision, allocatable :: work(:)
    !C add by K. Yang
    double precision :: agvv(3*natm,3)
    double precision :: gvv(3, natm, 3)
    double precision :: fad(3*natm,3,3)

    ! Determine rOH and determine ir, t, and dtdr
    roh = bndlen(7,13,x,y,z)
    call tent(nap,rap,roh,ir,t,dtdr)

    do i=1,3*natm
        du1dx(i) = 0.0d0
        du2dx(i) = 0.0d0
        du3dx(i) = 0.0d0
    enddo

    ! Determine u13 and u13 and their Cartesian gradients
    !      call evuij(nap,natm,igrad,ir,t,dtdr,X,Y,Z,u13,u23,
    !     $           du13dX,du23dX)
    ! Determin u12, u13, and u23 and their Cartesian gradients  !KRY 2014071
    call evuij(nap,natm,igrad,ir,t,dtdr,x,y,z,u12,u13,u23, &
        & du12dx,du13dx,du23dx)

    ! Determine u11, u22, u33, and their Cartesian gradients
    call evuii(igrad,natm,x,y,z,u1,du1dx,11)
    call evuii(igrad,natm,x,y,z,u2,du2dx,22)
    call evuii(igrad,natm,x,y,z,u3,du3dx,33)

    ! KRY Born-Mayer repulsion
    call evubm(igrad,natm,x,y,z,ubm,dubmdx)
    !
    ! add by J. Zheng
    guu(:,:,:,:) = 0d0
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

    uu(:,:) = 0d0
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

end subroutine poten

!> @brief Evaluate evuii for the PHOH_APRP_DPEM potential.
subroutine evuii( igrad, natm, x, y, z, u, dudx, ii )
    !***********************************************************************
    ! Subroutine to evaluate the diabatic potentials and gradients
    ! The bohr-hartree-radian unit is used throughout the program
    ! qtc:  redundant internal coordinates
    ! bmqt: B-matrix of RICs
    ! ii:   flag to indicate the diabatic states
    !       = 11, U11
    !       = 22, U22
    !       = 33, U33
    !***********************************************************************

    implicit none

    double precision :: bndang, bndlen, dihedr, dudph, dudr, dudth, phi, r, theta, ub, ur, ut, uto
    integer          :: i, ii, j

    integer, parameter :: maxdim = 100

    integer, parameter :: ldbmt  = maxdim
    integer, parameter :: ldbms  = maxdim
    integer, parameter :: ldbnrt = maxdim
    integer, parameter :: ldct   = maxdim

    integer, parameter :: nstc = 9

    double precision, parameter :: u1min = 0.00905097d0

    ! Global variables
    integer          :: igrad
    integer          :: natm
    double precision :: u
    double precision :: x(*),y(*),z(*)
    double precision :: dudx(3*natm)

    ! Local variables
    integer :: nqtc
    integer :: nblt,nbat,ntot,noopt
    data nblt/12/,nbat/18/ntot/24/,noopt/0/
    integer :: iblt(2,12),ibat(3,18),itot(4,24),ioopt(4,0)
    data iblt  /1,2, &
        & 1,6, &
        & 1,7, &
        & 2,3, &
        & 2,8, &
        & 3,4, &
        & 3,9, &
        & 4,5, &
        & 4,10, &
        & 5,6, &
        & 5,11, &
        & 6,12/
    data ibat  /2,1,6, &
        & 2,1,7, &
        & 6,1,7, &
        & 1,2,3, &
        & 1,2,8, &
        & 3,2,8, &
        & 2,3,4, &
        & 2,3,9, &
        & 4,3,9, &
        & 3,4,5, &
        & 3,4,10, &
        & 5,4,10, &
        & 4,5,6, &
        & 4,5,11, &
        & 6,5,11, &
        & 1,6,5, &
        & 1,6,12, &
        & 5,6,12/
    data itot  /6,1,2,3, &
        & 6,1,2,8, &
        & 7,1,2,3, &
        & 7,1,2,8, &
        & 2,1,6,5, &
        & 2,1,6,12, &
        & 7,1,6,5, &
        & 7,1,6,12, &
        & 1,2,3,4, &
        & 1,2,3,9, &
        & 8,2,3,4, &
        & 8,2,3,9, &
        & 2,3,4,5, &
        & 2,3,4,10, &
        & 9,3,4,5, &
        & 9,3,4,10, &
        & 3,4,5,6, &
        & 3,4,5,11, &
        & 10,4,5,6, &
        & 10,4,5,11, &
        & 4,5,6,1, &
        & 4,5,6,12, &
        & 11,5,6,1, &
        & 11,5,6,12/

    integer :: nbls,nbas,ntos,noops
    data nbls/1/,nbas/1/,ntos/1/,noops/0/

    integer :: ibls(2,1),ibas(3,1),itos(4,1),ioops(4,0)
    data ibls /7,13/
    data ibas /1,7,13/
    data itos /2,1,7,13/

    double precision :: qtc(maxdim),qsc(maxdim)
    double precision :: dutdr,dubdr,dutodr,durdr
    double precision :: dtodph,dubdth
    double precision :: dut(maxdim)
    double precision :: bmt(ldbmt,3*natm),bms(ldbms,3*natm)

    ! Evaluate internal coordinates
    call evintcoord(nblt,iblt,nbat,ibat,ntot,itot,noopt,ioopt, &
        & nqtc,qtc,x,y,z)

    r     = bndlen(7,13,x,y,z)
    theta = bndang(1,7,13,x,y,z)
    phi   = dihedr(2,1,7,13,x,y,z)

    ! Evaluate Wilson's B-matrix if igrad = 1
    if (igrad .eq. 1) then
        ! Added by J. Zheng
        !        bmt(:,:) = 0d0
        call bmat3(nblt,iblt,nbat,ibat,ntot,itot,noopt,ioopt,natm, &
            & x,y,z,bmt,ldbmt)
        ! Added by J. Zheng
        !        bms(:,:) = 0d0
        call bmat3(nbls,ibls,nbas,ibas,ntos,itos,noops,ioops,natm, &
            & x,y,z,bms,ldbms)
    endif

    ! Evaluate uii from tertiary coordinates
    call evuiit(igrad,r,nqtc,qtc,ut,dut,dutdr,ii)

    ! Evaluate uii from reactive (rOH) and secondary coordinates
    ! (bending C-O-H and torsion C-C-O-H)
    if (ii .eq. 11) then
        call evu1r(igrad,r,ur,durdr)
        call evu1b(igrad,r,theta,ub,dubdr,dubdth)
        call evu1to(igrad,r,phi,uto,dutodr,dtodph)
    else if (ii .eq. 22) then
        call evu2r(igrad,r,ur,durdr)
        call evu2b(igrad,r,theta,ub,dubdr,dubdth)
        call evu2to(igrad,r,phi,uto,dutodr,dtodph)
    else if (ii .eq. 33) then
        call evu3r(igrad,r,ur,durdr)
        call evu3b(igrad,r,theta,ub,dubdr,dubdth)
        call evu3to(igrad,r,phi,uto,dutodr,dtodph)
    endif

    ! Readujst zero of energy to have U1 min to be zero
    u     = ur+ub+uto+ut-u1min
    dudr  = durdr+dubdr+dutodr+dutdr
    dudth = dubdth
    dudph = dtodph

    ! Calculated analytical Cartesian gradient if igrad = 1
    if (igrad .eq. 1) then
        do i=1,3*natm
            dudx(i) = 0.0d0
        enddo
        ! Contribution of dudX from dudr
        do i=1,3*natm
            dudx(i) = dudx(i)+dudr*bms(1,i)
        enddo
        ! Contribution of dudX from dudth
        do i=1,3*natm
            dudx(i) = dudx(i)+dudth*bms(2,i)
        enddo
        ! Contribution of dudX from dudph
        do i=1,3*natm
            dudx(i) = dudx(i)+dudph*bms(3,i)
        enddo
        ! Contribution of dudX from dut
        do i=1,3*natm
            do j=1,nqtc
                dudx(i) = dudx(i)+dut(j)*bmt(j,i)
            enddo
        enddo
    endif

end subroutine evuii
!      subroutine evuiit(igrad,r,ntot,nqtc,qtc,ut,dut,dutdr,ii)
!> @brief Evaluate evuiit for the PHOH_APRP_DPEM potential.
subroutine evuiit( igrad, r, nqtc, qtc, ut, dut, dutdr, ii )
    !***********************************************************************
    ! Subroutine to evaluate the diabatic potentials and gradients
    ! The bohr-hartree-radian unit is used throughout the program
    ! qtc:  redundant internal coordinates
    ! ii:   flag to indicate the diabatic states
    !       = 11, U11
    !       = 22, U22
    !       = 33, U33
    !***********************************************************************

    implicit none

    integer, parameter :: nap = 4
    !      integer nqa,nqb
    !      parameter (nqa=30,nqb=24)
    !      parameter (nqs=12,nqb=18,nqt=24)
    integer, parameter :: nqs = 12
    integer, parameter :: nqb = 18
    integer, parameter :: nqt = 24

    ! Global variables
    integer          :: igrad,ii
    integer          :: ntot,nqtc
    double precision :: r,ut,dutdr
    double precision :: qtc(*),dut(*)

    ! Local variables
    integer          :: i,j,k
    integer          :: ir
    double precision :: u0,u01,u02,u1,u2
    double precision :: vs,vs1,vs2
    double precision :: vb,vb1,vb2
    double precision :: vt,vt1,vt2
    double precision :: vbs,vbs1,vbs2
    double precision :: dvsdr,dvbdr,dvtdr,dvbsdr
    double precision :: t(nap),dtdr(nap)
    double precision :: rap(nap)
    double precision :: v0(nap)
    double precision :: qs(nqs),qs1(nqs),qs2(nqs)
    double precision :: qb(nqb),qb1(nqb),qb2(nqb)
    double precision :: sqb(nqb)
    double precision :: qt(nqt),qt1(nqt),qt2(nqt)
    double precision :: sqt1(nqt),cqt1(nqt)
    double precision :: sqt2(nqt),cqt2(nqt)
    double precision :: dvs(nqs),dvs1(nqs),dvs2(nqs)
    double precision :: dvb(nqb),dvb1(nqb),dvb2(nqb)
    double precision :: dvt(nqt),dvt1(nqt),dvt2(nqt)
    double precision :: dvbsdb(nqb),dvbsdb1(nqb),dvbsdb2(nqb)
    double precision :: dvbsds(nqs),dvbsds1(nqs),dvbsds2(nqs)
    double precision :: qse(nqs,nap)
    double precision :: qbe(nqb,nap)
    double precision :: qte(nqt,nap)
    double precision :: qse1(nqs),qse2(nqs)
    double precision :: qbe1(nqb),qbe2(nqb)
    double precision :: fs1(nqs,nqs),fs2(nqs,nqs)
    double precision :: fb1(nqb,nqb),fb2(nqb,nqb)
    double precision :: ft1(nqt,nqt),ft2(nqt,nqt)
    double precision :: fbs1(nqb,nqs),fbs2(nqb,nqs)
    double precision :: hs(nqs,nqs,nap)
    double precision :: hb(nqb,nqb,nap)
    double precision :: ht(nqt,nqt,nap)
    double precision :: hbs(nqb,nqs,nap)

    ! Initialize output variables
    ut    = 0.0d0
    dutdr = 0.0d0
    do i=1,nqtc
        dut(i) = 0.0d0
    enddo

    ! Get HS for stretchings, HB for bendings, and HT for torsions
    call gethess(nap,hs,nqs,hb,nqb,ht,nqt,hbs,ii)
    ! Get equilibrium parameters
    call getgeom(nap,qse,nqs,qbe,nqb,qte,nqt,ii)

    if (ii .eq. 11) then
        rap(1) = 1.93206106d0 ! rOH=1.0224 using rOH=0.964A
        rap(2) = 2.49443845d0 ! rOH=1.32 A
        rap(3) = 3.77945219d0 ! rOH=2.00 A
        rap(4) = 9.44863047d0 ! rOH=5.00 A
        ! Normalized relaxation energy
        ! dE = E(DFT_opt) - E(CAS_geom)
        v0(1) = 0.0d0
        v0(2) = -0.000346019d0
        v0(3) = 0.000672861d0
        v0(4) = -0.001758495d0
    else if (ii .eq. 22) then
        rap(1) = 1.93206106d0 ! rOH=1.0224 using rOH=0.964A
        rap(2) = 2.49443845d0 ! rOH=1.32 A
        rap(3) = 4.27078097d0 ! rOH=2.26 A
        rap(4) = 9.44863047d0 ! rOH=5.00 A
        v0(1)  = -0.008940821d0
        v0(2)  = -0.00976278d0
        v0(3)  = -0.00976278d0
        v0(4)  = -0.00976278d0
    else if (ii .eq. 33) then
        rap(1) = 1.93206106d0 ! rOH=1.0224 using rOH=0.964A
        rap(2) = 2.49443845d0 ! rOH=1.32 A
        rap(3) = 4.27078097d0 ! rOH=2.26 A
        rap(4) = 9.44863047d0 ! rOH=5.00 A
        v0(1)  = -0.012423316d0
        v0(2)  = -0.015044491d0
        v0(3)  = -0.014671263d0
        v0(4)  = -0.01433271d0
    endif

    call tent(nap,rap,r,ir,t,dtdr)

    do i=1,nqs
        qs(i) = qtc(i)
    enddo
    do i=1,nqb
        qb(i)  = qtc(i+nqs)
        sqb(i) = sin(qb(i))
    enddo
    do i=1,nqt
        qt(i) = qtc(i+nqs+nqb)
    enddo

    ! Initialize hessian elements of stretchings, bendings, and torsions
    do i=1,nqs
        do j=1,i
            fs1(i,j) = 0.0d0
            fs2(i,j) = 0.0d0
        enddo
    enddo
    do i=1,nqb
        do j=1,i
            fb1(i,j) = 0.0d0
            fb2(i,j) = 0.0d0
        enddo
    enddo
    do i=1,nqt
        do j=1,i
            ft1(i,j) = 0.0d0
            ft2(i,j) = 0.0d0
        enddo
    enddo
    do i=1,nqb
        do j=1,nqs
            fbs1(i,j) = 0.0d0
            fbs2(i,j) = 0.0d0
        enddo
    enddo

    ! Get force constant matrices for strecthings (FS) and bending (FB)
    ! and torsions (FT)
    if (ir .eq. 1) then
        u0 = v0(1)
        do i=1,nqs
            !          qs1(i)=qs(i)-qse(i,1)
            qse1(i) = qse(i,1)
            do j=1,i
                fs1(i,j) = hs(i,j,1)*qse1(i)*qse1(j)
            enddo
        enddo
        do i=1,nqb
            !          qb1(i)=qb(i)-qbe(i,1)
            ! cos(theta0)-cos(theta) as variables
            qbe1(i) = qbe(i,1)
            qb1(i)  = cos(qbe1(i))-cos(qb(i))
            do j=1,i
                fb1(i,j) = hb(i,j,1)/(sin(qbe1(i))*sin(qbe1(j)))
            enddo
        enddo
        do i=1,nqt
            qt1(i)  = qt(i)-qte(i,1)
            sqt1(i) = sin(qt1(i))
            cqt1(i) = cos(qt1(i))
            do j=1,i
                ft1(i,j) = ht(i,j,1)
            enddo
        enddo
        ! bend-stretch couplings
        do i=1,nqb
            do j=1,nqs
                fbs1(i,j) = hbs(i,j,1)*qse1(j)/sin(qbe1(i))
            enddo
        enddo
    else if (ir .eq. nap+1) then
        u0 = v0(nap)
        do i=1,nqs
            qse1(i) = qse(i,nap)
            do j=1,i
                fs1(i,j) = hs(i,j,nap)*qse1(i)*qse1(j)
            enddo
        enddo
        do i=1,nqb
            ! cos(theta0)-cos(theta) as variables
            qbe1(i) = qbe(i,nap)
            qb1(i)  = cos(qbe1(i))-cos(qb(i))
            do j=1,i
                fb1(i,j) = hb(i,j,nap)/(sin(qbe1(i))*sin(qbe1(j)))
            enddo
        enddo
        do i=1,nqt
            qt1(i)  = qt(i)-qte(i,nap)
            sqt1(i) = sin(qt1(i))
            cqt1(i) = cos(qt1(i))
            do j=1,i
                ft1(i,j) = ht(i,j,nap)
            enddo
        enddo
        ! bend-stretch couplings
        do i=1,nqb
            do j=1,nqs
                fbs1(i,j) = hbs(i,j,nap)*qse1(j)/sin(qbe1(i))
            enddo
        enddo
    else if ((ir .gt. 1) .and. (ir .lt. nap+1)) then
        u01 = v0(ir-1)
        u02 = v0(ir)
        do i=1,nqs
            qse1(i) = qse(i,ir-1)
            qse2(i) = qse(i,ir)
            do j=1,i
                fs1(i,j) = hs(i,j,ir-1)*qse1(i)*qse1(j)
                fs2(i,j) = hs(i,j,ir)*qse2(i)*qse2(j)
            enddo
        enddo
        do i=1,nqb
            qbe1(i) = qbe(i,ir-1)
            qbe2(i) = qbe(i,ir)
            qb1(i)  = cos(qbe1(i))-cos(qb(i))
            qb2(i)  = cos(qbe2(i))-cos(qb(i))
            do j=1,i
                fb1(i,j) = hb(i,j,ir-1)/(sin(qbe1(i))*sin(qbe1(j)))
                fb2(i,j) = hb(i,j,ir)/(sin(qbe2(i))*sin(qbe2(j)))
            enddo
        enddo
        do i=1,nqt
            qt1(i)  = qt(i)-qte(i,ir-1)
            qt2(i)  = qt(i)-qte(i,ir)
            sqt1(i) = sin(qt1(i))
            sqt2(i) = sin(qt2(i))
            cqt1(i) = cos(qt1(i))
            cqt2(i) = cos(qt2(i))
            do j=1,i
                ft1(i,j) = ht(i,j,ir-1)
                ft2(i,j) = ht(i,j,ir)
            enddo
        enddo
        ! bend-stretch couplings
        do i=1,nqb
            do j=1,nqs
                fbs1(i,j) = hbs(i,j,ir-1)*qse1(j)/sin(qbe1(i))
                fbs2(i,j) = hbs(i,j,ir)*qse2(j)/sin(qbe2(i))
            enddo
        enddo
    endif

    ! Get the upper parts of force matrix by symmetry
    do i=1,nqs
        do j=i+1,nqs
            fs1(i,j) = fs1(j,i)
            fs2(i,j) = fs2(j,i)
        enddo
    enddo
    do i=1,nqb
        do j=i+1,nqb
            fb1(i,j) = fb1(j,i)
            fb2(i,j) = fb2(j,i)
        enddo
    enddo
    do i=1,nqt
        do j=i+1,nqt
            ft1(i,j) = ft1(j,i)
            ft2(i,j) = ft2(j,i)
        enddo
    enddo

    !      write(*,'(A,I2)') "ir = ",ir
    ! Calculate potential energy
    if ((ir .eq. 1) .or. (ir .eq. nap+1)) then
        vs  = 0.0d0
        vb  = 0.0d0
        vt  = 0.0d0
        vbs = 0.0d0
        do i=1,nqs
            do j=1,nqs
                !            vs=vs+0.5d0*FS1(i,j)*(qs(i)-qse1(i))*(qs(j)-qse1(j))
                vs = vs+0.5d0*fs1(i,j) *(qs(i)-qse1(i))*(qs(j)-qse1(j))/(qs(i)*qs(j))
            enddo
        enddo
        do i=1,nqb
            do j=1,nqb
                vb = vb+0.5d0*qb1(i)*fb1(i,j)*qb1(j)
            enddo
        enddo
        do i=1,nqt
            do j=1,nqt
                vt = vt+0.5d0*ft1(i,j)*sqt1(i)*sqt1(j)
            enddo
            vt = vt-0.5d0*ft1(i,i)*sqt1(i)*sqt1(i) +2.0d0*ft1(i,i)*sin(qt1(i)/2.0d0)*sin(qt1(i)/2.0d0)
        enddo
        ! bend-stretch couplings
        do i=1,nqb
            do j=1,nqs
                !            vbs=vbs+FBS1(i,j)*qb1(i)*(qs(j)-qse1(j))
                vbs = vbs+fbs1(i,j)*qb1(i)*(qs(j)-qse1(j))/qs(j)
            enddo
        enddo

        ut = u0+vs+vb+vt+vbs
        if (igrad .eq. 1) then
            !          dvsdr=0.0d0
            !          dvbdr=0.0d0
            !          dvtdr=0.0d0
            !          dvbsdr=0.0d0
            do i=1,nqs
                dvs(i) = 0.0d0
                do j=1,nqs
                    !              dvs(i)=dvs(i)+FS1(i,j)*(qs(j)-qse1(j))
                    dvs(i) = dvs(i)+fs1(i,j)*qse1(i) *(qs(j)-qse1(j))/(qs(i)*qs(i)*qs(j))
                enddo
            enddo
            do i=1,nqb
                dvb(i) = 0.0d0
                do j=1,nqb
                    dvb(i) = dvb(i)+sqb(i)*fb1(i,j)*qb1(j)
                enddo
            enddo
            do i=1,nqt
                dvt(i) = 0.0d0
                do j=1,nqt
                    dvt(i) = dvt(i)+ft1(i,j)*cqt1(i)*sqt1(j)
                enddo
                dvt(i) = dvt(i)-ft1(i,i)*cqt1(i)*sqt1(i) +2.0d0*ft1(i,i)*cos(qt1(i)/2.0d0)*sin(qt1(i)/2.0d0)
            enddo
            ! contribution from vbs (dvs and dvb)
            do i=1,nqb
                dvbsdb(i) = 0.0d0
                do j=1,nqs
                    !              dvbsdb(i)=dvbsdb(i)+FBS1(i,j)*sqb(i)*(qs(j)-qse1(j))
                    dvbsdb(i) = dvbsdb(i)+fbs1(i,j)*sqb(i) *(qs(j)-qse1(j))/qs(j)
                enddo
            enddo
            do i=1,nqs
                dvbsds(i) = 0.0d0
                do j=1,nqb
                    !              dvbsds(i)=dvbsds(i)+FBS1(j,i)*qb1(j)
                    dvbsds(i) = dvbsds(i)+fbs1(j,i)*qse1(i) *qb1(j)/(qs(i)*qs(i))
                enddo
            enddo

            !          dutdr=dvsdr+dvbdr+dvtdr+dvbsdr
            dutdr = 0.0d0
            do i=1,nqs
                dut(i) = dvs(i)+dvbsds(i)
            enddo
            do i=1,nqb
                dut(nqs+i) = dvb(i)+dvbsdb(i)
            enddo
            do i=1,nqt
                dut(nqs+nqb+i) = dvt(i)
            enddo
        endif
    else
        vs1  = 0.0d0
        vs2  = 0.0d0
        vb1  = 0.0d0
        vb2  = 0.0d0
        vt1  = 0.0d0
        vt2  = 0.0d0
        vbs1 = 0.0d0
        vbs2 = 0.0d0
        do i=1,nqs
            do j=1,nqs
                !            vs1=vs1+0.5d0*FS1(i,j)*(qs(i)-qse1(i))*(qs(j)-qse1(j))
                vs1 = vs1+0.5d0*fs1(i,j) *(qs(i)-qse1(i))*(qs(j)-qse1(j))/(qs(i)*qs(j))
                !            vs2=vs2+0.5d0*FS2(i,j)*(qs(i)-qse2(i))*(qs(j)-qse2(j))
                vs2 = vs2+0.5d0*fs2(i,j) *(qs(i)-qse2(i))*(qs(j)-qse2(j))/(qs(i)*qs(j))
            enddo
        enddo
        do i=1,nqb
            do j=1,nqb
                vb1 = vb1+0.5d0*qb1(i)*fb1(i,j)*qb1(j)
                vb2 = vb2+0.5d0*qb2(i)*fb2(i,j)*qb2(j)
            enddo
        enddo
        do i=1,nqt
            do j=1,nqt
                vt1 = vt1+0.5d0*sqt1(i)*ft1(i,j)*sqt1(j)
                vt2 = vt2+0.5d0*sqt2(i)*ft2(i,j)*sqt2(j)
            enddo
            vt1 = vt1-0.5d0*ft1(i,i)*sqt1(i)*sqt1(i) +2.0d0*ft1(i,i)*sin(qt1(i)/2.0d0)*sin(qt1(i)/2.0d0)
            vt2 = vt2-0.5d0*ft2(i,i)*sqt2(i)*sqt2(i) +2.0d0*ft2(i,i)*sin(qt2(i)/2.0d0)*sin(qt2(i)/2.0d0)
        enddo
        ! bend-stretch couplings
        do i=1,nqb
            do j=1,nqs
                !            vbs1=vbs1+qb1(i)*FBS1(i,j)*(qs(j)-qse1(j))
                !            vbs2=vbs2+qb2(i)*FBS2(i,j)*(qs(j)-qse2(j))
                vbs1 = vbs1+fbs1(i,j)*qb1(i)*(qs(j)-qse1(j))/qs(j)
                vbs2 = vbs2+fbs2(i,j)*qb2(i)*(qs(j)-qse2(j))/qs(j)
            enddo
        enddo

        u1 = u01+vs1+vb1+vt1+vbs1
        u2 = u02+vs2+vb2+vt2+vbs2
        ut = u1*t(ir-1)+u2*t(ir)
        if (igrad .eq. 1) then
            dutdr = u1*dtdr(ir-1)+u2*dtdr(ir)
            do i=1,nqs
                dvs1(i) = 0.0d0
                dvs2(i) = 0.0d0
                do j=1,nqs
                    !              dvs1(i)=dvs1(i)+FS1(i,j)*(qs(j)-qse1(j))
                    !              dvs2(i)=dvs2(i)+FS2(i,j)*(qs(j)-qse2(j))
                    dvs1(i) = dvs1(i)+fs1(i,j)*qse1(i) *(qs(j)-qse1(j))/(qs(i)*qs(i)*qs(j))
                    dvs2(i) = dvs2(i)+fs2(i,j)*qse2(i) *(qs(j)-qse2(j))/(qs(i)*qs(i)*qs(j))
                enddo
                dvs(i) = dvs1(i)*t(ir-1)+dvs2(i)*t(ir)
            enddo
            do i=1,nqb
                dvb1(i) = 0.0d0
                dvb2(i) = 0.0d0
                do j=1,nqb
                    dvb1(i) = dvb1(i)+sqb(i)*fb1(i,j)*qb1(j)
                    dvb2(i) = dvb2(i)+sqb(i)*fb2(i,j)*qb2(j)
                enddo
                dvb(i) = dvb1(i)*t(ir-1)+dvb2(i)*t(ir)
            enddo
            do i=1,nqt
                dvt1(i) = 0.0d0
                dvt2(i) = 0.0d0
                do j=1,nqt
                    dvt1(i) = dvt1(i)+ft1(i,j)*cqt1(i)*sqt1(j)
                    dvt2(i) = dvt2(i)+ft2(i,j)*cqt2(i)*sqt2(j)
                enddo
                dvt1(i) = dvt1(i)-ft1(i,i)*cqt1(i)*sqt1(i) +2.0d0*ft1(i,i)*sin(qt1(i)/2.0d0)*cos(qt1(i)/2.0d0)
                dvt2(i) = dvt2(i)-ft2(i,i)*cqt2(i)*sqt2(i) +2.0d0*ft2(i,i)*sin(qt2(i)/2.0d0)*cos(qt2(i)/2.0d0)
                dvt(i)  = dvt1(i)*t(ir-1)+dvt2(i)*t(ir)
            enddo
            ! contribution from vbs (dvbsdb and dvbsds)
            do i=1,nqb
                dvbsdb1(i) = 0.0d0
                dvbsdb2(i) = 0.0d0
                do j=1,nqs
                    !              dvbsdb1(i)=dvbsdb1(i)+FBS1(i,j)*sqb(i)*(qs(j)-qse1(j))
                    !              dvbsdb2(i)=dvbsdb2(i)+FBS2(i,j)*sqb(i)*(qs(j)-qse2(j))
                    dvbsdb1(i) = dvbsdb1(i)+fbs1(i,j)*sqb(i) *(qs(j)-qse1(j))/qs(j)
                    dvbsdb2(i) = dvbsdb2(i)+fbs2(i,j)*sqb(i) *(qs(j)-qse2(j))/qs(j)
                enddo
                dvbsdb(i) = dvbsdb1(i)*t(ir-1)+dvbsdb2(i)*t(ir)
            enddo
            do i=1,nqs
                dvbsds1(i) = 0.0d0
                dvbsds2(i) = 0.0d0
                do j=1,nqb
                    !              dvbsds1(i)=dvbsds1(i)+FBS1(j,i)*qb1(j)
                    !              dvbsds2(i)=dvbsds2(i)+FBS2(j,i)*qb2(j)
                    dvbsds1(i) = dvbsds1(i)+fbs1(j,i)*qse1(i) *qb1(j)/(qs(i)*qs(i))
                    dvbsds2(i) = dvbsds2(i)+fbs2(j,i)*qse2(i) *qb2(j)/(qs(i)*qs(i))
                enddo
                dvbsds(i) = dvbsds1(i)*t(ir-1)+dvbsds2(i)*t(ir)
            enddo

            do i=1,nqs
                dut(i) = dvs(i)+dvbsds(i)
            enddo
            do i=1,nqb
                dut(nqs+i) = dvb(i)+dvbsdb(i)
            enddo
            do i=1,nqt
                dut(nqs+nqb+i) = dvt(i)
            enddo
        endif
    endif

end subroutine evuiit
!> @brief Evaluate gethess for the PHOH_APRP_DPEM potential.
subroutine gethess( nap, hss, nqs, hbb, nqb, htt, nqt, hbs, ii )

    ! Global variables
    implicit none

    integer          :: nap
    integer          :: ii
    integer          :: nqs,nqb,nqt
    double precision :: hss(nqs,nqs,nap)
    double precision :: hbb(nqb,nqb,nap)
    double precision :: htt(nqt,nqt,nap)
    double precision :: hbs(nqb,nqs,nap)

    if (ii .eq. 11) then
        ! Hessian elements for anchor point 1
        hss( 1, 1, 1) = 0.40759d+00
        hss( 2, 1, 1) = 0.17500d-01
        hss( 2, 2, 1) = 0.40901d+00
        hss( 3, 1, 1) = 0.29920d-01
        hss( 3, 2, 1) = 0.35360d-01
        hss( 3, 3, 1) = 0.42477d+00
        hss( 4, 1, 1) = 0.18000d-01
        hss( 4, 2, 1) = -0.26000d-02
        hss( 4, 3, 1) = 0.38000d-03
        hss( 4, 4, 1) = 0.41176d+00
        hss( 5, 1, 1) = 0.45000d-02
        hss( 5, 2, 1) = -0.29000d-03
        hss( 5, 3, 1) = 0.30900d-02
        hss( 5, 4, 1) = 0.51800d-02
        hss( 5, 5, 1) = 0.35919d+00
        hss( 6, 1, 1) = -0.34400d-02
        hss( 6, 2, 1) = 0.82770d-01
        hss( 6, 3, 1) = -0.25600d-02
        hss( 6, 4, 1) = 0.21010d-01
        hss( 6, 5, 1) = -0.76000d-03
        hss( 6, 6, 1) = 0.41463d+00
        hss( 7, 1, 1) = -0.55000d-03
        hss( 7, 2, 1) = -0.34000d-03
        hss( 7, 3, 1) = 0.68000d-03
        hss( 7, 4, 1) = 0.44000d-02
        hss( 7, 5, 1) = 0.48000d-03
        hss( 7, 6, 1) = 0.45900d-02
        hss( 7, 7, 1) = 0.36428d+00
        hss( 8, 1, 1) = 0.82220d-01
        hss( 8, 2, 1) = -0.40500d-02
        hss( 8, 3, 1) = -0.11000d-03
        hss( 8, 4, 1) = -0.33000d-02
        hss( 8, 5, 1) = -0.30000d-03
        hss( 8, 6, 1) = 0.20930d-01
        hss( 8, 7, 1) = -0.31000d-03
        hss( 8, 8, 1) = 0.40778d+00
        hss( 9, 1, 1) = -0.48000d-03
        hss( 9, 2, 1) = -0.49000d-03
        hss( 9, 3, 1) = -0.23000d-03
        hss( 9, 4, 1) = -0.57000d-03
        hss( 9, 5, 1) = 0.20000d-03
        hss( 9, 6, 1) = 0.44500d-02
        hss( 9, 7, 1) = 0.47000d-03
        hss( 9, 8, 1) = 0.44300d-02
        hss( 9, 9, 1) = 0.36674d+00
        hss(10, 1, 1) = -0.18900d-02
        hss(10, 2, 1) = 0.18160d-01
        hss(10, 3, 1) = -0.10000d-02
        hss(10, 4, 1) = 0.86110d-01
        hss(10, 5, 1) = -0.74000d-03
        hss(10, 6, 1) = -0.29100d-02
        hss(10, 7, 1) = -0.54000d-03
        hss(10, 8, 1) = 0.20960d-01
        hss(10, 9, 1) = -0.51000d-03
        hss(10,10, 1) = 0.41852d+00
        hss(11, 1, 1) = -0.43000d-03
        hss(11, 2, 1) = -0.50000d-03
        hss(11, 3, 1) = 0.78000d-03
        hss(11, 4, 1) = -0.50000d-03
        hss(11, 5, 1) = 0.11000d-03
        hss(11, 6, 1) = -0.28000d-03
        hss(11, 7, 1) = 0.20000d-03
        hss(11, 8, 1) = 0.45700d-02
        hss(11, 9, 1) = 0.47000d-03
        hss(11,10, 1) = 0.45000d-02
        hss(11,11, 1) = 0.36422d+00
        hss(12, 1, 1) = -0.40000d-03
        hss(12, 2, 1) = 0.35500d-02
        hss(12, 3, 1) = 0.10000d-04
        hss(12, 4, 1) = -0.80000d-03
        hss(12, 5, 1) = 0.23000d-03
        hss(12, 6, 1) = -0.14000d-03
        hss(12, 7, 1) = 0.10000d-03
        hss(12, 8, 1) = -0.68000d-03
        hss(12, 9, 1) = 0.14000d-03
        hss(12,10, 1) = 0.42300d-02
        hss(12,11, 1) = 0.44000d-03
        hss(12,12, 1) = 0.36727d+00
        ! Hessian elements for anchor point 1
        hbb( 1, 1, 1) = 0.86570d-01
        hbb( 2, 1, 1) = -0.44420d-01
        hbb( 2, 2, 1) = 0.13902d+00
        hbb( 3, 1, 1) = -0.42150d-01
        hbb( 3, 2, 1) = -0.94600d-01
        hbb( 3, 3, 1) = 0.13674d+00
        hbb( 4, 1, 1) = -0.36220d-01
        hbb( 4, 2, 1) = 0.23510d-01
        hbb( 4, 3, 1) = 0.12710d-01
        hbb( 4, 4, 1) = 0.79210d-01
        hbb( 5, 1, 1) = 0.20780d-01
        hbb( 5, 2, 1) = -0.13590d-01
        hbb( 5, 3, 1) = -0.71800d-02
        hbb( 5, 4, 1) = -0.38760d-01
        hbb( 5, 5, 1) = 0.77300d-01
        hbb( 6, 1, 1) = 0.15440d-01
        hbb( 6, 2, 1) = -0.99200d-02
        hbb( 6, 3, 1) = -0.55200d-02
        hbb( 6, 4, 1) = -0.40440d-01
        hbb( 6, 5, 1) = -0.38530d-01
        hbb( 6, 6, 1) = 0.78980d-01
        hbb( 7, 1, 1) = -0.10540d-01
        hbb( 7, 2, 1) = -0.24300d-02
        hbb( 7, 3, 1) = 0.12970d-01
        hbb( 7, 4, 1) = -0.31280d-01
        hbb( 7, 5, 1) = 0.12070d-01
        hbb( 7, 6, 1) = 0.19210d-01
        hbb( 7, 7, 1) = 0.79110d-01
        hbb( 8, 1, 1) = -0.42000d-03
        hbb( 8, 2, 1) = 0.59300d-02
        hbb( 8, 3, 1) = -0.55100d-02
        hbb( 8, 4, 1) = 0.18520d-01
        hbb( 8, 5, 1) = -0.65400d-02
        hbb( 8, 6, 1) = -0.11980d-01
        hbb( 8, 7, 1) = -0.39530d-01
        hbb( 8, 8, 1) = 0.79370d-01
        hbb( 9, 1, 1) = 0.10960d-01
        hbb( 9, 2, 1) = -0.35000d-02
        hbb( 9, 3, 1) = -0.74600d-02
        hbb( 9, 4, 1) = 0.12760d-01
        hbb( 9, 5, 1) = -0.55300d-02
        hbb( 9, 6, 1) = -0.72300d-02
        hbb( 9, 7, 1) = -0.39590d-01
        hbb( 9, 8, 1) = -0.39850d-01
        hbb( 9, 9, 1) = 0.79440d-01
        hbb(10, 1, 1) = 0.69200d-02
        hbb(10, 2, 1) = -0.37400d-02
        hbb(10, 3, 1) = -0.31900d-02
        hbb(10, 4, 1) = -0.97700d-02
        hbb(10, 5, 1) = 0.10500d-01
        hbb(10, 6, 1) = -0.72000d-03
        hbb(10, 7, 1) = -0.33920d-01
        hbb(10, 8, 1) = 0.13770d-01
        hbb(10, 9, 1) = 0.20160d-01
        hbb(10,10, 1) = 0.80580d-01
        hbb(11, 1, 1) = -0.33900d-02
        hbb(11, 2, 1) = 0.23200d-02
        hbb(11, 3, 1) = 0.10700d-02
        hbb(11, 4, 1) = -0.78000d-03
        hbb(11, 5, 1) = -0.37400d-02
        hbb(11, 6, 1) = 0.45200d-02
        hbb(11, 7, 1) = 0.19900d-01
        hbb(11, 8, 1) = -0.73600d-02
        hbb(11, 9, 1) = -0.12550d-01
        hbb(11,10, 1) = -0.40290d-01
        hbb(11,11, 1) = 0.79280d-01
        hbb(12, 1, 1) = -0.35300d-02
        hbb(12, 2, 1) = 0.14200d-02
        hbb(12, 3, 1) = 0.21200d-02
        hbb(12, 4, 1) = 0.10550d-01
        hbb(12, 5, 1) = -0.67600d-02
        hbb(12, 6, 1) = -0.37900d-02
        hbb(12, 7, 1) = 0.14020d-01
        hbb(12, 8, 1) = -0.64100d-02
        hbb(12, 9, 1) = -0.76100d-02
        hbb(12,10, 1) = -0.40290d-01
        hbb(12,11, 1) = -0.38990d-01
        hbb(12,12, 1) = 0.79280d-01
        hbb(13, 1, 1) = -0.10320d-01
        hbb(13, 2, 1) = 0.13110d-01
        hbb(13, 3, 1) = -0.27900d-02
        hbb(13, 4, 1) = 0.52200d-02
        hbb(13, 5, 1) = -0.23700d-02
        hbb(13, 6, 1) = -0.28500d-02
        hbb(13, 7, 1) = -0.86400d-02
        hbb(13, 8, 1) = 0.10140d-01
        hbb(13, 9, 1) = -0.15000d-02
        hbb(13,10, 1) = -0.34100d-01
        hbb(13,11, 1) = 0.14060d-01
        hbb(13,12, 1) = 0.20050d-01
        hbb(13,13, 1) = 0.79730d-01
        hbb(14, 1, 1) = 0.10860d-01
        hbb(14, 2, 1) = -0.76500d-02
        hbb(14, 3, 1) = -0.32000d-02
        hbb(14, 4, 1) = -0.26900d-02
        hbb(14, 5, 1) = 0.11000d-02
        hbb(14, 6, 1) = 0.16000d-02
        hbb(14, 7, 1) = -0.15200d-02
        hbb(14, 8, 1) = -0.35000d-02
        hbb(14, 9, 1) = 0.50200d-02
        hbb(14,10, 1) = 0.20210d-01
        hbb(14,11, 1) = -0.76100d-02
        hbb(14,12, 1) = -0.12600d-01
        hbb(14,13, 1) = -0.39850d-01
        hbb(14,14, 1) = 0.79680d-01
        hbb(15, 1, 1) = -0.54000d-03
        hbb(15, 2, 1) = -0.54500d-02
        hbb(15, 3, 1) = 0.59900d-02
        hbb(15, 4, 1) = -0.25300d-02
        hbb(15, 5, 1) = 0.12700d-02
        hbb(15, 6, 1) = 0.12600d-02
        hbb(15, 7, 1) = 0.10150d-01
        hbb(15, 8, 1) = -0.66300d-02
        hbb(15, 9, 1) = -0.35200d-02
        hbb(15,10, 1) = 0.13890d-01
        hbb(15,11, 1) = -0.64500d-02
        hbb(15,12, 1) = -0.74400d-02
        hbb(15,13, 1) = -0.39880d-01
        hbb(15,14, 1) = -0.39830d-01
        hbb(15,15, 1) = 0.79710d-01
        hbb(16, 1, 1) = -0.36420d-01
        hbb(16, 2, 1) = 0.13970d-01
        hbb(16, 3, 1) = 0.22450d-01
        hbb(16, 4, 1) = -0.71600d-02
        hbb(16, 5, 1) = -0.22100d-02
        hbb(16, 6, 1) = 0.93700d-02
        hbb(16, 7, 1) = 0.52700d-02
        hbb(16, 8, 1) = -0.24800d-02
        hbb(16, 9, 1) = -0.27900d-02
        hbb(16,10, 1) = -0.97100d-02
        hbb(16,11, 1) = 0.10500d-01
        hbb(16,12, 1) = -0.79000d-03
        hbb(16,13, 1) = -0.31890d-01
        hbb(16,14, 1) = 0.12990d-01
        hbb(16,15, 1) = 0.18900d-01
        hbb(16,16, 1) = 0.79900d-01
        hbb(17, 1, 1) = 0.22220d-01
        hbb(17, 2, 1) = -0.60200d-02
        hbb(17, 3, 1) = -0.16200d-01
        hbb(17, 4, 1) = -0.27800d-02
        hbb(17, 5, 1) = 0.56900d-02
        hbb(17, 6, 1) = -0.29000d-02
        hbb(17, 7, 1) = -0.25500d-02
        hbb(17, 8, 1) = 0.13300d-02
        hbb(17, 9, 1) = 0.12200d-02
        hbb(17,10, 1) = 0.10510d-01
        hbb(17,11, 1) = -0.68200d-02
        hbb(17,12, 1) = -0.37000d-02
        hbb(17,13, 1) = 0.12370d-01
        hbb(17,14, 1) = -0.57100d-02
        hbb(17,15, 1) = -0.66600d-02
        hbb(17,16, 1) = -0.39780d-01
        hbb(17,17, 1) = 0.76100d-01
        hbb(18, 1, 1) = 0.14200d-01
        hbb(18, 2, 1) = -0.79500d-02
        hbb(18, 3, 1) = -0.62400d-02
        hbb(18, 4, 1) = 0.99400d-02
        hbb(18, 5, 1) = -0.34800d-02
        hbb(18, 6, 1) = -0.64600d-02
        hbb(18, 7, 1) = -0.27200d-02
        hbb(18, 8, 1) = 0.11500d-02
        hbb(18, 9, 1) = 0.15700d-02
        hbb(18,10, 1) = -0.80000d-03
        hbb(18,11, 1) = -0.36800d-02
        hbb(18,12, 1) = 0.44900d-02
        hbb(18,13, 1) = 0.19510d-01
        hbb(18,14, 1) = -0.72700d-02
        hbb(18,15, 1) = -0.12240d-01
        hbb(18,16, 1) = -0.40120d-01
        hbb(18,17, 1) = -0.36330d-01
        hbb(18,18, 1) = 0.76450d-01
        ! Hessian elements for anchor point 1
        htt( 1, 1, 1) = 0.92400d-02
        htt( 2, 1, 1) = 0.38700d-02
        htt( 2, 2, 1) = 0.84000d-02
        htt( 3, 1, 1) = 0.67000d-03
        htt( 3, 2, 1) = -0.40800d-02
        htt( 3, 3, 1) = 0.11550d-01
        htt( 4, 1, 1) = -0.47000d-02
        htt( 4, 2, 1) = 0.44000d-03
        htt( 4, 3, 1) = 0.68000d-02
        htt( 4, 4, 1) = 0.11950d-01
        htt( 5, 1, 1) = -0.62200d-02
        htt( 5, 2, 1) = -0.44300d-02
        htt( 5, 3, 1) = 0.20700d-02
        htt( 5, 4, 1) = 0.38700d-02
        htt( 5, 5, 1) = 0.89900d-02
        htt( 6, 1, 1) = -0.43000d-02
        htt( 6, 2, 1) = -0.29300d-02
        htt( 6, 3, 1) = 0.28800d-02
        htt( 6, 4, 1) = 0.42500d-02
        htt( 6, 5, 1) = 0.35700d-02
        htt( 6, 6, 1) = 0.80900d-02
        htt( 7, 1, 1) = 0.19000d-02
        htt( 7, 2, 1) = 0.31100d-02
        htt( 7, 3, 1) = -0.82500d-02
        htt( 7, 4, 1) = -0.70400d-02
        htt( 7, 5, 1) = 0.11300d-02
        htt( 7, 6, 1) = -0.32500d-02
        htt( 7, 7, 1) = 0.10750d-01
        htt( 8, 1, 1) = 0.38300d-02
        htt( 8, 2, 1) = 0.46100d-02
        htt( 8, 3, 1) = -0.74400d-02
        htt( 8, 4, 1) = -0.66600d-02
        htt( 8, 5, 1) = -0.43000d-02
        htt( 8, 6, 1) = 0.12800d-02
        htt( 8, 7, 1) = 0.63800d-02
        htt( 8, 8, 1) = 0.11960d-01
        htt( 9, 1, 1) = -0.63700d-02
        htt( 9, 2, 1) = -0.83000d-03
        htt( 9, 3, 1) = -0.30300d-02
        htt( 9, 4, 1) = 0.25100d-02
        htt( 9, 5, 1) = 0.48000d-03
        htt( 9, 6, 1) = 0.20600d-02
        htt( 9, 7, 1) = -0.26900d-02
        htt( 9, 8, 1) = -0.11100d-02
        htt( 9, 9, 1) = 0.94700d-02
        htt(10, 1, 1) = -0.40200d-02
        htt(10, 2, 1) = 0.76000d-03
        htt(10, 3, 1) = -0.16700d-02
        htt(10, 4, 1) = 0.31200d-02
        htt(10, 5, 1) = 0.20200d-02
        htt(10, 6, 1) = 0.15300d-02
        htt(10, 7, 1) = -0.22000d-03
        htt(10, 8, 1) = -0.70000d-03
        htt(10, 9, 1) = 0.31100d-02
        htt(10,10, 1) = 0.83700d-02
        htt(11, 1, 1) = -0.98000d-03
        htt(11, 2, 1) = -0.53700d-02
        htt(11, 3, 1) = 0.17300d-02
        htt(11, 4, 1) = -0.26500d-02
        htt(11, 5, 1) = -0.13300d-02
        htt(11, 6, 1) = 0.68000d-03
        htt(11, 7, 1) = -0.39000d-02
        htt(11, 8, 1) = -0.19000d-02
        htt(11, 9, 1) = 0.39100d-02
        htt(11,10, 1) = -0.16900d-02
        htt(11,11, 1) = 0.83100d-02
        htt(12, 1, 1) = 0.13700d-02
        htt(12, 2, 1) = -0.37800d-02
        htt(12, 3, 1) = 0.31000d-02
        htt(12, 4, 1) = -0.20400d-02
        htt(12, 5, 1) = 0.22000d-03
        htt(12, 6, 1) = 0.16000d-03
        htt(12, 7, 1) = -0.14300d-02
        htt(12, 8, 1) = -0.14900d-02
        htt(12, 9, 1) = -0.24600d-02
        htt(12,10, 1) = 0.35600d-02
        htt(12,11, 1) = 0.27100d-02
        htt(12,12, 1) = 0.87300d-02
        htt(13, 1, 1) = 0.53000d-03
        htt(13, 2, 1) = -0.16100d-02
        htt(13, 3, 1) = 0.26100d-02
        htt(13, 4, 1) = 0.47000d-03
        htt(13, 5, 1) = 0.24000d-02
        htt(13, 6, 1) = 0.89000d-03
        htt(13, 7, 1) = 0.43000d-03
        htt(13, 8, 1) = -0.10800d-02
        htt(13, 9, 1) = -0.66000d-02
        htt(13,10, 1) = -0.21000d-03
        htt(13,11, 1) = -0.44500d-02
        htt(13,12, 1) = 0.19400d-02
        htt(13,13, 1) = 0.96100d-02
        htt(14, 1, 1) = 0.20400d-02
        htt(14, 2, 1) = 0.52000d-03
        htt(14, 3, 1) = 0.14800d-02
        htt(14, 4, 1) = -0.40000d-04
        htt(14, 5, 1) = 0.97000d-03
        htt(14, 6, 1) = -0.11900d-02
        htt(14, 7, 1) = 0.15000d-02
        htt(14, 8, 1) = -0.66000d-03
        htt(14, 9, 1) = -0.44200d-02
        htt(14,10, 1) = 0.12200d-02
        htt(14,11, 1) = -0.29000d-02
        htt(14,12, 1) = 0.27500d-02
        htt(14,13, 1) = 0.37500d-02
        htt(14,14, 1) = 0.85800d-02
        htt(15, 1, 1) = -0.18400d-02
        htt(15, 2, 1) = -0.32100d-02
        htt(15, 3, 1) = 0.12300d-02
        htt(15, 4, 1) = -0.15000d-03
        htt(15, 5, 1) = 0.85000d-03
        htt(15, 6, 1) = 0.14200d-02
        htt(15, 7, 1) = -0.20600d-02
        htt(15, 8, 1) = -0.14900d-02
        htt(15, 9, 1) = -0.19000d-03
        htt(15,10, 1) = -0.55100d-02
        htt(15,11, 1) = 0.11900d-02
        htt(15,12, 1) = -0.41300d-02
        htt(15,13, 1) = 0.31600d-02
        htt(15,14, 1) = -0.19400d-02
        htt(15,15, 1) = 0.85200d-02
        htt(16, 1, 1) = -0.33000d-03
        htt(16, 2, 1) = -0.10800d-02
        htt(16, 3, 1) = 0.10000d-03
        htt(16, 4, 1) = -0.66000d-03
        htt(16, 5, 1) = -0.59000d-03
        htt(16, 6, 1) = -0.66000d-03
        htt(16, 7, 1) = -0.99000d-03
        htt(16, 8, 1) = -0.10700d-02
        htt(16, 9, 1) = 0.19900d-02
        htt(16,10, 1) = -0.40700d-02
        htt(16,11, 1) = 0.27500d-02
        htt(16,12, 1) = -0.33200d-02
        htt(16,13, 1) = -0.26900d-02
        htt(16,14, 1) = 0.28900d-02
        htt(16,15, 1) = 0.34200d-02
        htt(16,16, 1) = 0.90000d-02
        htt(17, 1, 1) = 0.25300d-02
        htt(17, 2, 1) = 0.10400d-02
        htt(17, 3, 1) = 0.19000d-03
        htt(17, 4, 1) = -0.13000d-02
        htt(17, 5, 1) = 0.43000d-03
        htt(17, 6, 1) = -0.16400d-02
        htt(17, 7, 1) = 0.26600d-02
        htt(17, 8, 1) = 0.59000d-03
        htt(17, 9, 1) = 0.61000d-03
        htt(17,10, 1) = -0.18300d-02
        htt(17,11, 1) = 0.21000d-02
        htt(17,12, 1) = -0.33000d-03
        htt(17,13, 1) = -0.66400d-02
        htt(17,14, 1) = -0.70000d-03
        htt(17,15, 1) = -0.41800d-02
        htt(17,16, 1) = 0.17600d-02
        htt(17,17, 1) = 0.96800d-02
        htt(18, 1, 1) = 0.74000d-03
        htt(18, 2, 1) = 0.13600d-02
        htt(18, 3, 1) = -0.19200d-02
        htt(18, 4, 1) = -0.13000d-02
        htt(18, 5, 1) = -0.16100d-02
        htt(18, 6, 1) = -0.30500d-02
        htt(18, 7, 1) = 0.91000d-03
        htt(18, 8, 1) = -0.53000d-03
        htt(18, 9, 1) = 0.21300d-02
        htt(18,10, 1) = 0.45000d-03
        htt(18,11, 1) = 0.15100d-02
        htt(18,12, 1) = -0.17000d-03
        htt(18,13, 1) = -0.40600d-02
        htt(18,14, 1) = 0.10400d-02
        htt(18,15, 1) = -0.23700d-02
        htt(18,16, 1) = 0.27300d-02
        htt(18,17, 1) = 0.31700d-02
        htt(18,18, 1) = 0.85000d-02
        htt(19, 1, 1) = 0.10300d-02
        htt(19, 2, 1) = -0.10900d-02
        htt(19, 3, 1) = 0.13200d-02
        htt(19, 4, 1) = -0.79000d-03
        htt(19, 5, 1) = 0.18700d-02
        htt(19, 6, 1) = 0.45000d-03
        htt(19, 7, 1) = 0.15900d-02
        htt(19, 8, 1) = 0.17000d-03
        htt(19, 9, 1) = -0.15700d-02
        htt(19,10, 1) = -0.32700d-02
        htt(19,11, 1) = 0.55000d-03
        htt(19,12, 1) = -0.11400d-02
        htt(19,13, 1) = -0.78000d-03
        htt(19,14, 1) = -0.55300d-02
        htt(19,15, 1) = 0.93000d-03
        htt(19,16, 1) = -0.38200d-02
        htt(19,17, 1) = 0.37300d-02
        htt(19,18, 1) = -0.19300d-02
        htt(19,19, 1) = 0.84800d-02
        htt(20, 1, 1) = -0.77000d-03
        htt(20, 2, 1) = -0.77000d-03
        htt(20, 3, 1) = -0.79000d-03
        htt(20, 4, 1) = -0.79000d-03
        htt(20, 5, 1) = -0.18000d-03
        htt(20, 6, 1) = -0.97000d-03
        htt(20, 7, 1) = -0.16000d-03
        htt(20, 8, 1) = -0.95000d-03
        htt(20, 9, 1) = -0.50000d-04
        htt(20,10, 1) = -0.98000d-03
        htt(20,11, 1) = -0.50000d-04
        htt(20,12, 1) = -0.98000d-03
        htt(20,13, 1) = 0.18000d-02
        htt(20,14, 1) = -0.37900d-02
        htt(20,15, 1) = 0.27400d-02
        htt(20,16, 1) = -0.28500d-02
        htt(20,17, 1) = -0.27700d-02
        htt(20,18, 1) = 0.34000d-02
        htt(20,19, 1) = 0.28200d-02
        htt(20,20, 1) = 0.89900d-02
        htt(21, 1, 1) = 0.30000d-03
        htt(21, 2, 1) = 0.19600d-02
        htt(21, 3, 1) = -0.25100d-02
        htt(21, 4, 1) = -0.85000d-03
        htt(21, 5, 1) = -0.60900d-02
        htt(21, 6, 1) = -0.58000d-03
        htt(21, 7, 1) = -0.34300d-02
        htt(21, 8, 1) = 0.20800d-02
        htt(21, 9, 1) = 0.24300d-02
        htt(21,10, 1) = 0.93000d-03
        htt(21,11, 1) = 0.76000d-03
        htt(21,12, 1) = -0.74000d-03
        htt(21,13, 1) = 0.68000d-03
        htt(21,14, 1) = -0.16400d-02
        htt(21,15, 1) = 0.21900d-02
        htt(21,16, 1) = -0.13000d-03
        htt(21,17, 1) = -0.66000d-02
        htt(21,18, 1) = -0.36000d-03
        htt(21,19, 1) = -0.42700d-02
        htt(21,20, 1) = 0.19700d-02
        htt(21,21, 1) = 0.92800d-02
        htt(22, 1, 1) = -0.16800d-02
        htt(22, 2, 1) = 0.42000d-03
        htt(22, 3, 1) = -0.33400d-02
        htt(22, 4, 1) = -0.12400d-02
        htt(22, 5, 1) = -0.50000d-03
        htt(22, 6, 1) = -0.52300d-02
        htt(22, 7, 1) = 0.10700d-02
        htt(22, 8, 1) = -0.36600d-02
        htt(22, 9, 1) = 0.80000d-03
        htt(22,10, 1) = 0.14300d-02
        htt(22,11, 1) = -0.13000d-02
        htt(22,12, 1) = -0.68000d-03
        htt(22,13, 1) = 0.22300d-02
        htt(22,14, 1) = 0.58000d-03
        htt(22,15, 1) = 0.16000d-02
        htt(22,16, 1) = -0.60000d-04
        htt(22,17, 1) = -0.44700d-02
        htt(22,18, 1) = 0.11200d-02
        htt(22,19, 1) = -0.28100d-02
        htt(22,20, 1) = 0.27800d-02
        htt(22,21, 1) = 0.36100d-02
        htt(22,22, 1) = 0.84800d-02
        htt(23, 1, 1) = 0.20800d-02
        htt(23, 2, 1) = 0.16400d-02
        htt(23, 3, 1) = -0.41000d-03
        htt(23, 4, 1) = -0.85000d-03
        htt(23, 5, 1) = -0.40600d-02
        htt(23, 6, 1) = 0.83000d-03
        htt(23, 7, 1) = -0.16900d-02
        htt(23, 8, 1) = 0.31900d-02
        htt(23, 9, 1) = 0.92000d-03
        htt(23,10, 1) = -0.13400d-02
        htt(23,11, 1) = 0.13500d-02
        htt(23,12, 1) = -0.90000d-03
        htt(23,13, 1) = -0.18800d-02
        htt(23,14, 1) = -0.33700d-02
        htt(23,15, 1) = 0.39000d-03
        htt(23,16, 1) = -0.11000d-02
        htt(23,17, 1) = -0.14000d-03
        htt(23,18, 1) = -0.56500d-02
        htt(23,19, 1) = 0.13500d-02
        htt(23,20, 1) = -0.41600d-02
        htt(23,21, 1) = 0.30900d-02
        htt(23,22, 1) = -0.19500d-02
        htt(23,23, 1) = 0.85600d-02
        htt(24, 1, 1) = 0.10000d-03
        htt(24, 2, 1) = 0.11000d-03
        htt(24, 3, 1) = -0.12500d-02
        htt(24, 4, 1) = -0.12400d-02
        htt(24, 5, 1) = 0.15300d-02
        htt(24, 6, 1) = -0.38300d-02
        htt(24, 7, 1) = 0.28100d-02
        htt(24, 8, 1) = -0.25500d-02
        htt(24, 9, 1) = -0.71000d-03
        htt(24,10, 1) = -0.84000d-03
        htt(24,11, 1) = -0.71000d-03
        htt(24,12, 1) = -0.84000d-03
        htt(24,13, 1) = -0.32000d-03
        htt(24,14, 1) = -0.11500d-02
        htt(24,15, 1) = -0.20000d-03
        htt(24,16, 1) = -0.10200d-02
        htt(24,17, 1) = 0.19900d-02
        htt(24,18, 1) = -0.41700d-02
        htt(24,19, 1) = 0.28200d-02
        htt(24,20, 1) = -0.33500d-02
        htt(24,21, 1) = -0.25900d-02
        htt(24,22, 1) = 0.29200d-02
        htt(24,23, 1) = 0.35300d-02
        htt(24,24, 1) = 0.90500d-02
        ! Hessian elements for anchor point 1
        hbs( 1, 1, 1) = -0.36310d-01
        hbs( 1, 2, 1) = -0.36690d-01
        hbs( 1, 3, 1) = -0.36890d-01
        hbs( 1, 4, 1) = -0.11400d-01
        hbs( 1, 5, 1) = 0.55900d-02
        hbs( 1, 6, 1) = 0.45710d-01
        hbs( 1, 7, 1) = -0.44000d-03
        hbs( 1, 8, 1) = 0.44410d-01
        hbs( 1, 9, 1) = 0.60000d-04
        hbs( 1,10, 1) = -0.10640d-01
        hbs( 1,11, 1) = -0.46000d-03
        hbs( 1,12, 1) = 0.55900d-02
        hbs( 2, 1, 1) = 0.42410d-01
        hbs( 2, 2, 1) = -0.54500d-02
        hbs( 2, 3, 1) = 0.11600d-01
        hbs( 2, 4, 1) = -0.17300d-02
        hbs( 2, 5, 1) = -0.28000d-02
        hbs( 2, 6, 1) = -0.19690d-01
        hbs( 2, 7, 1) = 0.11300d-02
        hbs( 2, 8, 1) = -0.26090d-01
        hbs( 2, 9, 1) = -0.18000d-03
        hbs( 2,10, 1) = 0.11810d-01
        hbs( 2,11, 1) = -0.69000d-03
        hbs( 2,12, 1) = -0.11800d-02
        hbs( 3, 1, 1) = -0.61000d-02
        hbs( 3, 2, 1) = 0.42140d-01
        hbs( 3, 3, 1) = 0.25290d-01
        hbs( 3, 4, 1) = 0.13130d-01
        hbs( 3, 5, 1) = -0.27900d-02
        hbs( 3, 6, 1) = -0.26020d-01
        hbs( 3, 7, 1) = -0.68000d-03
        hbs( 3, 8, 1) = -0.18310d-01
        hbs( 3, 9, 1) = 0.12000d-03
        hbs( 3,10, 1) = -0.11700d-02
        hbs( 3,11, 1) = 0.11600d-02
        hbs( 3,12, 1) = -0.44100d-02
        hbs( 4, 1, 1) = -0.29210d-01
        hbs( 4, 2, 1) = -0.82700d-02
        hbs( 4, 3, 1) = 0.16570d-01
        hbs( 4, 4, 1) = -0.31510d-01
        hbs( 4, 5, 1) = -0.97200d-02
        hbs( 4, 6, 1) = -0.12980d-01
        hbs( 4, 7, 1) = 0.48900d-02
        hbs( 4, 8, 1) = 0.44650d-01
        hbs( 4, 9, 1) = -0.64000d-03
        hbs( 4,10, 1) = 0.43490d-01
        hbs( 4,11, 1) = 0.19000d-03
        hbs( 4,12, 1) = -0.10200d-02
        hbs( 5, 1, 1) = 0.27310d-01
        hbs( 5, 2, 1) = 0.47900d-02
        hbs( 5, 3, 1) = -0.10360d-01
        hbs( 5, 4, 1) = 0.23300d-02
        hbs( 5, 5, 1) = 0.65800d-02
        hbs( 5, 6, 1) = 0.59300d-02
        hbs( 5, 7, 1) = -0.19100d-02
        hbs( 5, 8, 1) = -0.19800d-01
        hbs( 5, 9, 1) = -0.33000d-03
        hbs( 5,10, 1) = -0.24210d-01
        hbs( 5,11, 1) = -0.10000d-03
        hbs( 5,12, 1) = 0.10100d-02
        hbs( 6, 1, 1) = 0.19000d-02
        hbs( 6, 2, 1) = 0.34900d-02
        hbs( 6, 3, 1) = -0.62000d-02
        hbs( 6, 4, 1) = 0.29170d-01
        hbs( 6, 5, 1) = 0.31400d-02
        hbs( 6, 6, 1) = 0.70400d-02
        hbs( 6, 7, 1) = -0.29800d-02
        hbs( 6, 8, 1) = -0.24840d-01
        hbs( 6, 9, 1) = 0.97000d-03
        hbs( 6,10, 1) = -0.19290d-01
        hbs( 6,11, 1) = -0.90000d-04
        hbs( 6,12, 1) = 0.10000d-04
        hbs( 7, 1, 1) = -0.12120d-01
        hbs( 7, 2, 1) = 0.42330d-01
        hbs( 7, 3, 1) = 0.37900d-02
        hbs( 7, 4, 1) = -0.32360d-01
        hbs( 7, 5, 1) = 0.55000d-02
        hbs( 7, 6, 1) = -0.33060d-01
        hbs( 7, 7, 1) = -0.89600d-02
        hbs( 7, 8, 1) = -0.12200d-01
        hbs( 7, 9, 1) = 0.50200d-02
        hbs( 7,10, 1) = 0.46020d-01
        hbs( 7,11, 1) = -0.96000d-03
        hbs( 7,12, 1) = -0.80000d-04
        hbs( 8, 1, 1) = 0.67000d-02
        hbs( 8, 2, 1) = -0.23490d-01
        hbs( 8, 3, 1) = -0.78000d-03
        hbs( 8, 4, 1) = 0.30060d-01
        hbs( 8, 5, 1) = -0.33500d-02
        hbs( 8, 6, 1) = 0.37600d-02
        hbs( 8, 7, 1) = 0.46600d-02
        hbs( 8, 8, 1) = 0.52800d-02
        hbs( 8, 9, 1) = -0.20300d-02
        hbs( 8,10, 1) = -0.20610d-01
        hbs( 8,11, 1) = -0.19000d-03
        hbs( 8,12, 1) = 0.30000d-04
        hbs( 9, 1, 1) = 0.54200d-02
        hbs( 9, 2, 1) = -0.18840d-01
        hbs( 9, 3, 1) = -0.30100d-02
        hbs( 9, 4, 1) = 0.22900d-02
        hbs( 9, 5, 1) = -0.21500d-02
        hbs( 9, 6, 1) = 0.29300d-01
        hbs( 9, 7, 1) = 0.43100d-02
        hbs( 9, 8, 1) = 0.69200d-02
        hbs( 9, 9, 1) = -0.29800d-02
        hbs( 9,10, 1) = -0.25420d-01
        hbs( 9,11, 1) = 0.11500d-02
        hbs( 9,12, 1) = 0.50000d-04
        hbs(10, 1, 1) = 0.44050d-01
        hbs(10, 2, 1) = 0.43680d-01
        hbs(10, 3, 1) = -0.27100d-02
        hbs(10, 4, 1) = -0.13440d-01
        hbs(10, 5, 1) = -0.64000d-03
        hbs(10, 6, 1) = -0.32210d-01
        hbs(10, 7, 1) = 0.53500d-02
        hbs(10, 8, 1) = -0.31660d-01
        hbs(10, 9, 1) = -0.87900d-02
        hbs(10,10, 1) = -0.13050d-01
        hbs(10,11, 1) = 0.52900d-02
        hbs(10,12, 1) = -0.50000d-03
        hbs(11, 1, 1) = -0.24480d-01
        hbs(11, 2, 1) = -0.19450d-01
        hbs(11, 3, 1) = 0.11700d-02
        hbs(11, 4, 1) = 0.71100d-02
        hbs(11, 5, 1) = 0.10200d-02
        hbs(11, 6, 1) = 0.30010d-01
        hbs(11, 7, 1) = -0.31900d-02
        hbs(11, 8, 1) = 0.20900d-02
        hbs(11, 9, 1) = 0.44300d-02
        hbs(11,10, 1) = 0.60900d-02
        hbs(11,11, 1) = -0.21400d-02
        hbs(11,12, 1) = -0.33000d-03
        hbs(12, 1, 1) = -0.19570d-01
        hbs(12, 2, 1) = -0.24220d-01
        hbs(12, 3, 1) = 0.15500d-02
        hbs(12, 4, 1) = 0.63300d-02
        hbs(12, 5, 1) = -0.38000d-03
        hbs(12, 6, 1) = 0.22000d-02
        hbs(12, 7, 1) = -0.21500d-02
        hbs(12, 8, 1) = 0.29570d-01
        hbs(12, 9, 1) = 0.43600d-02
        hbs(12,10, 1) = 0.69600d-02
        hbs(12,11, 1) = -0.31500d-02
        hbs(12,12, 1) = 0.84000d-03
        hbs(13, 1, 1) = 0.42020d-01
        hbs(13, 2, 1) = -0.11720d-01
        hbs(13, 3, 1) = 0.43200d-02
        hbs(13, 4, 1) = 0.45370d-01
        hbs(13, 5, 1) = -0.30000d-04
        hbs(13, 6, 1) = -0.12300d-01
        hbs(13, 7, 1) = -0.10000d-02
        hbs(13, 8, 1) = -0.32760d-01
        hbs(13, 9, 1) = 0.50300d-02
        hbs(13,10, 1) = -0.33540d-01
        hbs(13,11, 1) = -0.89700d-02
        hbs(13,12, 1) = 0.49300d-02
        hbs(14, 1, 1) = -0.18670d-01
        hbs(14, 2, 1) = 0.51600d-02
        hbs(14, 3, 1) = -0.32600d-02
        hbs(14, 4, 1) = -0.25130d-01
        hbs(14, 5, 1) = 0.30000d-04
        hbs(14, 6, 1) = 0.69900d-02
        hbs(14, 7, 1) = 0.11900d-02
        hbs(14, 8, 1) = 0.29240d-01
        hbs(14, 9, 1) = -0.29900d-02
        hbs(14,10, 1) = 0.30100d-02
        hbs(14,11, 1) = 0.43500d-02
        hbs(14,12, 1) = -0.19200d-02
        hbs(15, 1, 1) = -0.23350d-01
        hbs(15, 2, 1) = 0.65500d-02
        hbs(15, 3, 1) = -0.10600d-02
        hbs(15, 4, 1) = -0.20240d-01
        hbs(15, 5, 1) = 0.00000d+00
        hbs(15, 6, 1) = 0.53100d-02
        hbs(15, 7, 1) = -0.18000d-03
        hbs(15, 8, 1) = 0.35200d-02
        hbs(15, 9, 1) = -0.20400d-02
        hbs(15,10, 1) = 0.30530d-01
        hbs(15,11, 1) = 0.46200d-02
        hbs(15,12, 1) = -0.30100d-02
        hbs(16, 1, 1) = -0.84300d-02
        hbs(16, 2, 1) = -0.29330d-01
        hbs(16, 3, 1) = 0.14920d-01
        hbs(16, 4, 1) = 0.43340d-01
        hbs(16, 5, 1) = -0.70000d-03
        hbs(16, 6, 1) = 0.44840d-01
        hbs(16, 7, 1) = 0.17000d-03
        hbs(16, 8, 1) = -0.12440d-01
        hbs(16, 9, 1) = -0.67000d-03
        hbs(16,10, 1) = -0.32290d-01
        hbs(16,11, 1) = 0.49100d-02
        hbs(16,12, 1) = -0.89200d-02
        hbs(17, 1, 1) = 0.44500d-02
        hbs(17, 2, 1) = 0.25980d-01
        hbs(17, 3, 1) = -0.89400d-02
        hbs(17, 4, 1) = -0.24010d-01
        hbs(17, 5, 1) = 0.10600d-02
        hbs(17, 6, 1) = -0.19750d-01
        hbs(17, 7, 1) = -0.90000d-04
        hbs(17, 8, 1) = 0.59800d-02
        hbs(17, 9, 1) = -0.29000d-03
        hbs(17,10, 1) = 0.31400d-02
        hbs(17,11, 1) = -0.18900d-02
        hbs(17,12, 1) = 0.71100d-02
        hbs(18, 1, 1) = 0.39800d-02
        hbs(18, 2, 1) = 0.33600d-02
        hbs(18, 3, 1) = -0.59800d-02
        hbs(18, 4, 1) = -0.19330d-01
        hbs(18, 5, 1) = -0.36000d-03
        hbs(18, 6, 1) = -0.25090d-01
        hbs(18, 7, 1) = -0.80000d-04
        hbs(18, 8, 1) = 0.64600d-02
        hbs(18, 9, 1) = 0.96000d-03
        hbs(18,10, 1) = 0.29140d-01
        hbs(18,11, 1) = -0.30200d-02
        hbs(18,12, 1) = 0.18100d-02
        ! Hessian elements for anchor point 2
        hss( 1, 1, 2) = 0.40931d+00
        hss( 2, 1, 2) = 0.17750d-01
        hss( 2, 2, 2) = 0.39831d+00
        hss( 3, 1, 2) = 0.30090d-01
        hss( 3, 2, 2) = 0.39670d-01
        hss( 3, 3, 2) = 0.41039d+00
        hss( 4, 1, 2) = 0.17380d-01
        hss( 4, 2, 2) = -0.24100d-02
        hss( 4, 3, 2) = -0.47000d-03
        hss( 4, 4, 2) = 0.41129d+00
        hss( 5, 1, 2) = 0.42800d-02
        hss( 5, 2, 2) = -0.46000d-03
        hss( 5, 3, 2) = 0.49400d-02
        hss( 5, 4, 2) = 0.55000d-02
        hss( 5, 5, 2) = 0.35488d+00
        hss( 6, 1, 2) = -0.36200d-02
        hss( 6, 2, 2) = 0.81920d-01
        hss( 6, 3, 2) = -0.25500d-02
        hss( 6, 4, 2) = 0.21100d-01
        hss( 6, 5, 2) = -0.76000d-03
        hss( 6, 6, 2) = 0.41554d+00
        hss( 7, 1, 2) = -0.59000d-03
        hss( 7, 2, 2) = -0.36000d-03
        hss( 7, 3, 2) = 0.73000d-03
        hss( 7, 4, 2) = 0.43500d-02
        hss( 7, 5, 2) = 0.47000d-03
        hss( 7, 6, 2) = 0.46000d-02
        hss( 7, 7, 2) = 0.36429d+00
        hss( 8, 1, 2) = 0.81370d-01
        hss( 8, 2, 2) = -0.43100d-02
        hss( 8, 3, 2) = 0.46000d-03
        hss( 8, 4, 2) = -0.31800d-02
        hss( 8, 5, 2) = -0.28000d-03
        hss( 8, 6, 2) = 0.21080d-01
        hss( 8, 7, 2) = -0.32000d-03
        hss( 8, 8, 2) = 0.40703d+00
        hss( 9, 1, 2) = -0.47000d-03
        hss( 9, 2, 2) = -0.50000d-03
        hss( 9, 3, 2) = -0.34000d-03
        hss( 9, 4, 2) = -0.61000d-03
        hss( 9, 5, 2) = 0.23000d-03
        hss( 9, 6, 2) = 0.44800d-02
        hss( 9, 7, 2) = 0.47000d-03
        hss( 9, 8, 2) = 0.44100d-02
        hss( 9, 9, 2) = 0.36691d+00
        hss(10, 1, 2) = -0.11900d-02
        hss(10, 2, 2) = 0.18110d-01
        hss(10, 3, 2) = -0.13500d-02
        hss(10, 4, 2) = 0.85990d-01
        hss(10, 5, 2) = -0.79000d-03
        hss(10, 6, 2) = -0.28600d-02
        hss(10, 7, 2) = -0.55000d-03
        hss(10, 8, 2) = 0.21070d-01
        hss(10, 9, 2) = -0.50000d-03
        hss(10,10, 2) = 0.41886d+00
        hss(11, 1, 2) = -0.47000d-03
        hss(11, 2, 2) = -0.52000d-03
        hss(11, 3, 2) = 0.68000d-03
        hss(11, 4, 2) = -0.53000d-03
        hss(11, 5, 2) = 0.12000d-03
        hss(11, 6, 2) = -0.24000d-03
        hss(11, 7, 2) = 0.21000d-03
        hss(11, 8, 2) = 0.45400d-02
        hss(11, 9, 2) = 0.47000d-03
        hss(11,10, 2) = 0.45200d-02
        hss(11,11, 2) = 0.36397d+00
        hss(12, 1, 2) = -0.53000d-03
        hss(12, 2, 2) = 0.41000d-02
        hss(12, 3, 2) = -0.53000d-03
        hss(12, 4, 2) = -0.83000d-03
        hss(12, 5, 2) = 0.34000d-03
        hss(12, 6, 2) = -0.14000d-03
        hss(12, 7, 2) = 0.10000d-03
        hss(12, 8, 2) = -0.69000d-03
        hss(12, 9, 2) = 0.14000d-03
        hss(12,10, 2) = 0.42400d-02
        hss(12,11, 2) = 0.46000d-03
        hss(12,12, 2) = 0.36700d+00
        ! Hessian elements for anchor point 2
        hbb( 1, 1, 2) = 0.85740d-01
        hbb( 2, 1, 2) = -0.47120d-01
        hbb( 2, 2, 2) = 0.14272d+00
        hbb( 3, 1, 2) = -0.38620d-01
        hbb( 3, 2, 2) = -0.95590d-01
        hbb( 3, 3, 2) = 0.13421d+00
        hbb( 4, 1, 2) = -0.35340d-01
        hbb( 4, 2, 2) = 0.24370d-01
        hbb( 4, 3, 2) = 0.10970d-01
        hbb( 4, 4, 2) = 0.78390d-01
        hbb( 5, 1, 2) = 0.19660d-01
        hbb( 5, 2, 2) = -0.13050d-01
        hbb( 5, 3, 2) = -0.66100d-02
        hbb( 5, 4, 2) = -0.38100d-01
        hbb( 5, 5, 2) = 0.76190d-01
        hbb( 6, 1, 2) = 0.15680d-01
        hbb( 6, 2, 2) = -0.11320d-01
        hbb( 6, 3, 2) = -0.43700d-02
        hbb( 6, 4, 2) = -0.40290d-01
        hbb( 6, 5, 2) = -0.38090d-01
        hbb( 6, 6, 2) = 0.78390d-01
        hbb( 7, 1, 2) = -0.10880d-01
        hbb( 7, 2, 2) = -0.19300d-02
        hbb( 7, 3, 2) = 0.12810d-01
        hbb( 7, 4, 2) = -0.30710d-01
        hbb( 7, 5, 2) = 0.11850d-01
        hbb( 7, 6, 2) = 0.18860d-01
        hbb( 7, 7, 2) = 0.78510d-01
        hbb( 8, 1, 2) = -0.26000d-03
        hbb( 8, 2, 2) = 0.56200d-02
        hbb( 8, 3, 2) = -0.53700d-02
        hbb( 8, 4, 2) = 0.18260d-01
        hbb( 8, 5, 2) = -0.64700d-02
        hbb( 8, 6, 2) = -0.11790d-01
        hbb( 8, 7, 2) = -0.39180d-01
        hbb( 8, 8, 2) = 0.78940d-01
        hbb( 9, 1, 2) = 0.11130d-01
        hbb( 9, 2, 2) = -0.36900d-02
        hbb( 9, 3, 2) = -0.74400d-02
        hbb( 9, 4, 2) = 0.12450d-01
        hbb( 9, 5, 2) = -0.53800d-02
        hbb( 9, 6, 2) = -0.70800d-02
        hbb( 9, 7, 2) = -0.39340d-01
        hbb( 9, 8, 2) = -0.39760d-01
        hbb( 9, 9, 2) = 0.79100d-01
        hbb(10, 1, 2) = 0.73300d-02
        hbb(10, 2, 2) = -0.43700d-02
        hbb(10, 3, 2) = -0.29600d-02
        hbb(10, 4, 2) = -0.10030d-01
        hbb(10, 5, 2) = 0.10560d-01
        hbb(10, 6, 2) = -0.53000d-03
        hbb(10, 7, 2) = -0.33760d-01
        hbb(10, 8, 2) = 0.13640d-01
        hbb(10, 9, 2) = 0.20120d-01
        hbb(10,10, 2) = 0.80240d-01
        hbb(11, 1, 2) = -0.36200d-02
        hbb(11, 2, 2) = 0.25500d-02
        hbb(11, 3, 2) = 0.10700d-02
        hbb(11, 4, 2) = -0.61000d-03
        hbb(11, 5, 2) = -0.38000d-02
        hbb(11, 6, 2) = 0.44100d-02
        hbb(11, 7, 2) = 0.19800d-01
        hbb(11, 8, 2) = -0.72600d-02
        hbb(11, 9, 2) = -0.12540d-01
        hbb(11,10, 2) = -0.40120d-01
        hbb(11,11, 2) = 0.79050d-01
        hbb(12, 1, 2) = -0.37100d-02
        hbb(12, 2, 2) = 0.18200d-02
        hbb(12, 3, 2) = 0.18900d-02
        hbb(12, 4, 2) = 0.10640d-01
        hbb(12, 5, 2) = -0.67600d-02
        hbb(12, 6, 2) = -0.38800d-02
        hbb(12, 7, 2) = 0.13960d-01
        hbb(12, 8, 2) = -0.63800d-02
        hbb(12, 9, 2) = -0.75800d-02
        hbb(12,10, 2) = -0.40110d-01
        hbb(12,11, 2) = -0.38920d-01
        hbb(12,12, 2) = 0.79040d-01
        hbb(13, 1, 2) = -0.10540d-01
        hbb(13, 2, 2) = 0.13360d-01
        hbb(13, 3, 2) = -0.28100d-02
        hbb(13, 4, 2) = 0.52100d-02
        hbb(13, 5, 2) = -0.22200d-02
        hbb(13, 6, 2) = -0.29900d-02
        hbb(13, 7, 2) = -0.87200d-02
        hbb(13, 8, 2) = 0.10150d-01
        hbb(13, 9, 2) = -0.14300d-02
        hbb(13,10, 2) = -0.34040d-01
        hbb(13,11, 2) = 0.14040d-01
        hbb(13,12, 2) = 0.19990d-01
        hbb(13,13, 2) = 0.79490d-01
        hbb(14, 1, 2) = 0.10840d-01
        hbb(14, 2, 2) = -0.78900d-02
        hbb(14, 3, 2) = -0.29500d-02
        hbb(14, 4, 2) = -0.26600d-02
        hbb(14, 5, 2) = 0.10300d-02
        hbb(14, 6, 2) = 0.16200d-02
        hbb(14, 7, 2) = -0.14800d-02
        hbb(14, 8, 2) = -0.35100d-02
        hbb(14, 9, 2) = 0.49900d-02
        hbb(14,10, 2) = 0.20150d-01
        hbb(14,11, 2) = -0.75900d-02
        hbb(14,12, 2) = -0.12560d-01
        hbb(14,13, 2) = -0.39740d-01
        hbb(14,14, 2) = 0.79580d-01
        hbb(15, 1, 2) = -0.30000d-03
        hbb(15, 2, 2) = -0.54600d-02
        hbb(15, 3, 2) = 0.57600d-02
        hbb(15, 4, 2) = -0.25500d-02
        hbb(15, 5, 2) = 0.11900d-02
        hbb(15, 6, 2) = 0.13700d-02
        hbb(15, 7, 2) = 0.10200d-01
        hbb(15, 8, 2) = -0.66400d-02
        hbb(15, 9, 2) = -0.35600d-02
        hbb(15,10, 2) = 0.13890d-01
        hbb(15,11, 2) = -0.64600d-02
        hbb(15,12, 2) = -0.74300d-02
        hbb(15,13, 2) = -0.39750d-01
        hbb(15,14, 2) = -0.39830d-01
        hbb(15,15, 2) = 0.79580d-01
        hbb(16, 1, 2) = -0.36310d-01
        hbb(16, 2, 2) = 0.15710d-01
        hbb(16, 3, 2) = 0.20610d-01
        hbb(16, 4, 2) = -0.75200d-02
        hbb(16, 5, 2) = -0.17400d-02
        hbb(16, 6, 2) = 0.92600d-02
        hbb(16, 7, 2) = 0.55500d-02
        hbb(16, 8, 2) = -0.26200d-02
        hbb(16, 9, 2) = -0.29400d-02
        hbb(16,10, 2) = -0.97400d-02
        hbb(16,11, 2) = 0.10510d-01
        hbb(16,12, 2) = -0.77000d-03
        hbb(16,13, 2) = -0.31400d-01
        hbb(16,14, 2) = 0.12890d-01
        hbb(16,15, 2) = 0.18510d-01
        hbb(16,16, 2) = 0.79420d-01
        hbb(17, 1, 2) = 0.22000d-01
        hbb(17, 2, 2) = -0.74700d-02
        hbb(17, 3, 2) = -0.14540d-01
        hbb(17, 4, 2) = -0.24600d-02
        hbb(17, 5, 2) = 0.53500d-02
        hbb(17, 6, 2) = -0.28900d-02
        hbb(17, 7, 2) = -0.26800d-02
        hbb(17, 8, 2) = 0.14000d-02
        hbb(17, 9, 2) = 0.12800d-02
        hbb(17,10, 2) = 0.10490d-01
        hbb(17,11, 2) = -0.68100d-02
        hbb(17,12, 2) = -0.36800d-02
        hbb(17,13, 2) = 0.12230d-01
        hbb(17,14, 2) = -0.56800d-02
        hbb(17,15, 2) = -0.65500d-02
        hbb(17,16, 2) = -0.39570d-01
        hbb(17,17, 2) = 0.76420d-01
        hbb(18, 1, 2) = 0.14310d-01
        hbb(18, 2, 2) = -0.82400d-02
        hbb(18, 3, 2) = -0.60700d-02
        hbb(18, 4, 2) = 0.99800d-02
        hbb(18, 5, 2) = -0.36100d-02
        hbb(18, 6, 2) = -0.63700d-02
        hbb(18, 7, 2) = -0.28700d-02
        hbb(18, 8, 2) = 0.12200d-02
        hbb(18, 9, 2) = 0.16500d-02
        hbb(18,10, 2) = -0.75000d-03
        hbb(18,11, 2) = -0.37000d-02
        hbb(18,12, 2) = 0.44500d-02
        hbb(18,13, 2) = 0.19170d-01
        hbb(18,14, 2) = -0.72100d-02
        hbb(18,15, 2) = -0.11960d-01
        hbb(18,16, 2) = -0.39850d-01
        hbb(18,17, 2) = -0.36850d-01
        hbb(18,18, 2) = 0.76700d-01
        ! Hessian elements for anchor point 2
        htt( 1, 1, 2) = 0.90500d-02
        htt( 2, 1, 2) = 0.38100d-02
        htt( 2, 2, 2) = 0.81400d-02
        htt( 3, 1, 2) = 0.25000d-03
        htt( 3, 2, 2) = -0.43200d-02
        htt( 3, 3, 2) = 0.11860d-01
        htt( 4, 1, 2) = -0.49900d-02
        htt( 4, 2, 2) = 0.10000d-04
        htt( 4, 3, 2) = 0.72800d-02
        htt( 4, 4, 2) = 0.12280d-01
        htt( 5, 1, 2) = -0.60200d-02
        htt( 5, 2, 2) = -0.42900d-02
        htt( 5, 3, 2) = 0.24900d-02
        htt( 5, 4, 2) = 0.42200d-02
        htt( 5, 5, 2) = 0.87300d-02
        htt( 6, 1, 2) = -0.42000d-02
        htt( 6, 2, 2) = -0.28900d-02
        htt( 6, 3, 2) = 0.31600d-02
        htt( 6, 4, 2) = 0.44800d-02
        htt( 6, 5, 2) = 0.34600d-02
        htt( 6, 6, 2) = 0.79900d-02
        htt( 7, 1, 2) = 0.21800d-02
        htt( 7, 2, 2) = 0.32900d-02
        htt( 7, 3, 2) = -0.83100d-02
        htt( 7, 4, 2) = -0.72000d-02
        htt( 7, 5, 2) = 0.81000d-03
        htt( 7, 6, 2) = -0.34000d-02
        htt( 7, 7, 2) = 0.10570d-01
        htt( 8, 1, 2) = 0.39900d-02
        htt( 8, 2, 2) = 0.46900d-02
        htt( 8, 3, 2) = -0.76400d-02
        htt( 8, 4, 2) = -0.69400d-02
        htt( 8, 5, 2) = -0.44700d-02
        htt( 8, 6, 2) = 0.11300d-02
        htt( 8, 7, 2) = 0.63600d-02
        htt( 8, 8, 2) = 0.11960d-01
        htt( 9, 1, 2) = -0.62900d-02
        htt( 9, 2, 2) = -0.85000d-03
        htt( 9, 3, 2) = -0.28600d-02
        htt( 9, 4, 2) = 0.25800d-02
        htt( 9, 5, 2) = 0.40000d-03
        htt( 9, 6, 2) = 0.20200d-02
        htt( 9, 7, 2) = -0.27900d-02
        htt( 9, 8, 2) = -0.11700d-02
        htt( 9, 9, 2) = 0.94300d-02
        htt(10, 1, 2) = -0.39800d-02
        htt(10, 2, 2) = 0.72000d-03
        htt(10, 3, 2) = -0.15100d-02
        htt(10, 4, 2) = 0.31900d-02
        htt(10, 5, 2) = 0.20000d-02
        htt(10, 6, 2) = 0.15200d-02
        htt(10, 7, 2) = -0.30000d-03
        htt(10, 8, 2) = -0.78000d-03
        htt(10, 9, 2) = 0.30800d-02
        htt(10,10, 2) = 0.83100d-02
        htt(11, 1, 2) = -0.99000d-03
        htt(11, 2, 2) = -0.52300d-02
        htt(11, 3, 2) = 0.17700d-02
        htt(11, 4, 2) = -0.24700d-02
        htt(11, 5, 2) = -0.13500d-02
        htt(11, 6, 2) = 0.69000d-03
        htt(11, 7, 2) = -0.39100d-02
        htt(11, 8, 2) = -0.18800d-02
        htt(11, 9, 2) = 0.39300d-02
        htt(11,10, 2) = -0.16800d-02
        htt(11,11, 2) = 0.82200d-02
        htt(12, 1, 2) = 0.13200d-02
        htt(12, 2, 2) = -0.36600d-02
        htt(12, 3, 2) = 0.31200d-02
        htt(12, 4, 2) = -0.18600d-02
        htt(12, 5, 2) = 0.24000d-03
        htt(12, 6, 2) = 0.19000d-03
        htt(12, 7, 2) = -0.14300d-02
        htt(12, 8, 2) = -0.14900d-02
        htt(12, 9, 2) = -0.24200d-02
        htt(12,10, 2) = 0.35600d-02
        htt(12,11, 2) = 0.26100d-02
        htt(12,12, 2) = 0.85900d-02
        htt(13, 1, 2) = 0.48000d-03
        htt(13, 2, 2) = -0.16000d-02
        htt(13, 3, 2) = 0.26700d-02
        htt(13, 4, 2) = 0.58000d-03
        htt(13, 5, 2) = 0.24300d-02
        htt(13, 6, 2) = 0.89000d-03
        htt(13, 7, 2) = 0.40000d-03
        htt(13, 8, 2) = -0.11500d-02
        htt(13, 9, 2) = -0.65700d-02
        htt(13,10, 2) = -0.18000d-03
        htt(13,11, 2) = -0.44600d-02
        htt(13,12, 2) = 0.19300d-02
        htt(13,13, 2) = 0.95900d-02
        htt(14, 1, 2) = 0.19800d-02
        htt(14, 2, 2) = 0.50000d-03
        htt(14, 3, 2) = 0.14700d-02
        htt(14, 4, 2) = -0.10000d-04
        htt(14, 5, 2) = 0.10300d-02
        htt(14, 6, 2) = -0.11500d-02
        htt(14, 7, 2) = 0.15100d-02
        htt(14, 8, 2) = -0.68000d-03
        htt(14, 9, 2) = -0.44200d-02
        htt(14,10, 2) = 0.12200d-02
        htt(14,11, 2) = -0.29300d-02
        htt(14,12, 2) = 0.27100d-02
        htt(14,13, 2) = 0.37900d-02
        htt(14,14, 2) = 0.85500d-02
        htt(15, 1, 2) = -0.18400d-02
        htt(15, 2, 2) = -0.31900d-02
        htt(15, 3, 2) = 0.13100d-02
        htt(15, 4, 2) = -0.40000d-04
        htt(15, 5, 2) = 0.83000d-03
        htt(15, 6, 2) = 0.14000d-02
        htt(15, 7, 2) = -0.21100d-02
        htt(15, 8, 2) = -0.15400d-02
        htt(15, 9, 2) = -0.16000d-03
        htt(15,10, 2) = -0.54600d-02
        htt(15,11, 2) = 0.12000d-02
        htt(15,12, 2) = -0.41000d-02
        htt(15,13, 2) = 0.31500d-02
        htt(15,14, 2) = -0.19100d-02
        htt(15,15, 2) = 0.85000d-02
        htt(16, 1, 2) = -0.35000d-03
        htt(16, 2, 2) = -0.10800d-02
        htt(16, 3, 2) = 0.11000d-03
        htt(16, 4, 2) = -0.62000d-03
        htt(16, 5, 2) = -0.57000d-03
        htt(16, 6, 2) = -0.64000d-03
        htt(16, 7, 2) = -0.10000d-02
        htt(16, 8, 2) = -0.10700d-02
        htt(16, 9, 2) = 0.19800d-02
        htt(16,10, 2) = -0.40600d-02
        htt(16,11, 2) = 0.27200d-02
        htt(16,12, 2) = -0.33200d-02
        htt(16,13, 2) = -0.26500d-02
        htt(16,14, 2) = 0.28600d-02
        htt(16,15, 2) = 0.34400d-02
        htt(16,16, 2) = 0.89500d-02
        htt(17, 1, 2) = 0.25600d-02
        htt(17, 2, 2) = 0.11000d-02
        htt(17, 3, 2) = 0.12000d-03
        htt(17, 4, 2) = -0.13400d-02
        htt(17, 5, 2) = 0.36000d-03
        htt(17, 6, 2) = -0.16300d-02
        htt(17, 7, 2) = 0.26400d-02
        htt(17, 8, 2) = 0.64000d-03
        htt(17, 9, 2) = 0.63000d-03
        htt(17,10, 2) = -0.18200d-02
        htt(17,11, 2) = 0.21000d-02
        htt(17,12, 2) = -0.35000d-03
        htt(17,13, 2) = -0.66700d-02
        htt(17,14, 2) = -0.74000d-03
        htt(17,15, 2) = -0.42000d-02
        htt(17,16, 2) = 0.17300d-02
        htt(17,17, 2) = 0.96700d-02
        htt(18, 1, 2) = 0.74000d-03
        htt(18, 2, 2) = 0.13700d-02
        htt(18, 3, 2) = -0.19900d-02
        htt(18, 4, 2) = -0.13600d-02
        htt(18, 5, 2) = -0.16000d-02
        htt(18, 6, 2) = -0.30000d-02
        htt(18, 7, 2) = 0.94000d-03
        htt(18, 8, 2) = -0.46000d-03
        htt(18, 9, 2) = 0.21100d-02
        htt(18,10, 2) = 0.44000d-03
        htt(18,11, 2) = 0.14700d-02
        htt(18,12, 2) = -0.19000d-03
        htt(18,13, 2) = -0.40300d-02
        htt(18,14, 2) = 0.10300d-02
        htt(18,15, 2) = -0.23600d-02
        htt(18,16, 2) = 0.27100d-02
        htt(18,17, 2) = 0.31700d-02
        htt(18,18, 2) = 0.84700d-02
        htt(19, 1, 2) = 0.10700d-02
        htt(19, 2, 2) = -0.10000d-02
        htt(19, 3, 2) = 0.13200d-02
        htt(19, 4, 2) = -0.75000d-03
        htt(19, 5, 2) = 0.17600d-02
        htt(19, 6, 2) = 0.41000d-03
        htt(19, 7, 2) = 0.15300d-02
        htt(19, 8, 2) = 0.17000d-03
        htt(19, 9, 2) = -0.15200d-02
        htt(19,10, 2) = -0.32300d-02
        htt(19,11, 2) = 0.58000d-03
        htt(19,12, 2) = -0.11300d-02
        htt(19,13, 2) = -0.86000d-03
        htt(19,14, 2) = -0.55100d-02
        htt(19,15, 2) = 0.86000d-03
        htt(19,16, 2) = -0.37900d-02
        htt(19,17, 2) = 0.37400d-02
        htt(19,18, 2) = -0.18900d-02
        htt(19,19, 2) = 0.84000d-02
        htt(20, 1, 2) = -0.75000d-03
        htt(20, 2, 2) = -0.74000d-03
        htt(20, 3, 2) = -0.79000d-03
        htt(20, 4, 2) = -0.78000d-03
        htt(20, 5, 2) = -0.20000d-03
        htt(20, 6, 2) = -0.96000d-03
        htt(20, 7, 2) = -0.17000d-03
        htt(20, 8, 2) = -0.93000d-03
        htt(20, 9, 2) = -0.40000d-04
        htt(20,10, 2) = -0.97000d-03
        htt(20,11, 2) = -0.50000d-04
        htt(20,12, 2) = -0.98000d-03
        htt(20,13, 2) = 0.17700d-02
        htt(20,14, 2) = -0.37400d-02
        htt(20,15, 2) = 0.27000d-02
        htt(20,16, 2) = -0.28100d-02
        htt(20,17, 2) = -0.27500d-02
        htt(20,18, 2) = 0.34100d-02
        htt(20,19, 2) = 0.27600d-02
        htt(20,20, 2) = 0.89200d-02
        htt(21, 1, 2) = 0.21000d-03
        htt(21, 2, 2) = 0.18300d-02
        htt(21, 3, 2) = -0.26900d-02
        htt(21, 4, 2) = -0.10600d-02
        htt(21, 5, 2) = -0.59300d-02
        htt(21, 6, 2) = -0.53000d-03
        htt(21, 7, 2) = -0.32400d-02
        htt(21, 8, 2) = 0.21600d-02
        htt(21, 9, 2) = 0.24300d-02
        htt(21,10, 2) = 0.91000d-03
        htt(21,11, 2) = 0.79000d-03
        htt(21,12, 2) = -0.73000d-03
        htt(21,13, 2) = 0.69000d-03
        htt(21,14, 2) = -0.16600d-02
        htt(21,15, 2) = 0.22200d-02
        htt(21,16, 2) = -0.12000d-03
        htt(21,17, 2) = -0.65400d-02
        htt(21,18, 2) = -0.37000d-03
        htt(21,19, 2) = -0.42000d-02
        htt(21,20, 2) = 0.19700d-02
        htt(21,21, 2) = 0.91600d-02
        htt(22, 1, 2) = -0.16600d-02
        htt(22, 2, 2) = 0.40000d-03
        htt(22, 3, 2) = -0.33700d-02
        htt(22, 4, 2) = -0.13200d-02
        htt(22, 5, 2) = -0.52000d-03
        htt(22, 6, 2) = -0.51900d-02
        htt(22, 7, 2) = 0.10800d-02
        htt(22, 8, 2) = -0.35900d-02
        htt(22, 9, 2) = 0.77000d-03
        htt(22,10, 2) = 0.14000d-02
        htt(22,11, 2) = -0.13100d-02
        htt(22,12, 2) = -0.68000d-03
        htt(22,13, 2) = 0.22700d-02
        htt(22,14, 2) = 0.59000d-03
        htt(22,15, 2) = 0.16300d-02
        htt(22,16, 2) = -0.50000d-04
        htt(22,17, 2) = -0.44900d-02
        htt(22,18, 2) = 0.10700d-02
        htt(22,19, 2) = -0.28100d-02
        htt(22,20, 2) = 0.27500d-02
        htt(22,21, 2) = 0.36200d-02
        htt(22,22, 2) = 0.84200d-02
        htt(23, 1, 2) = 0.20100d-02
        htt(23, 2, 2) = 0.15700d-02
        htt(23, 3, 2) = -0.60000d-03
        htt(23, 4, 2) = -0.10400d-02
        htt(23, 5, 2) = -0.39900d-02
        htt(23, 6, 2) = 0.82000d-03
        htt(23, 7, 2) = -0.15600d-02
        htt(23, 8, 2) = 0.32500d-02
        htt(23, 9, 2) = 0.97000d-03
        htt(23,10, 2) = -0.13300d-02
        htt(23,11, 2) = 0.14100d-02
        htt(23,12, 2) = -0.89000d-03
        htt(23,13, 2) = -0.19300d-02
        htt(23,14, 2) = -0.34100d-02
        htt(23,15, 2) = 0.39000d-03
        htt(23,16, 2) = -0.11000d-02
        htt(23,17, 2) = -0.90000d-04
        htt(23,18, 2) = -0.56300d-02
        htt(23,19, 2) = 0.13900d-02
        htt(23,20, 2) = -0.41400d-02
        htt(23,21, 2) = 0.30400d-02
        htt(23,22, 2) = -0.18900d-02
        htt(23,23, 2) = 0.85300d-02
        htt(24, 1, 2) = 0.15000d-03
        htt(24, 2, 2) = 0.14000d-03
        htt(24, 3, 2) = -0.12900d-02
        htt(24, 4, 2) = -0.13000d-02
        htt(24, 5, 2) = 0.14300d-02
        htt(24, 6, 2) = -0.38400d-02
        htt(24, 7, 2) = 0.27700d-02
        htt(24, 8, 2) = -0.25000d-02
        htt(24, 9, 2) = -0.70000d-03
        htt(24,10, 2) = -0.84000d-03
        htt(24,11, 2) = -0.68000d-03
        htt(24,12, 2) = -0.83000d-03
        htt(24,13, 2) = -0.34000d-03
        htt(24,14, 2) = -0.11700d-02
        htt(24,15, 2) = -0.19000d-03
        htt(24,16, 2) = -0.10200d-02
        htt(24,17, 2) = 0.19500d-02
        htt(24,18, 2) = -0.41900d-02
        htt(24,19, 2) = 0.27800d-02
        htt(24,20, 2) = -0.33600d-02
        htt(24,21, 2) = -0.25000d-02
        htt(24,22, 2) = 0.29000d-02
        htt(24,23, 2) = 0.36000d-02
        htt(24,24, 2) = 0.90000d-02
        ! Hessian elements for anchor point 2
        hbs( 1, 1, 2) = -0.37150d-01
        hbs( 1, 2, 2) = -0.35080d-01
        hbs( 1, 3, 2) = -0.37820d-01
        hbs( 1, 4, 2) = -0.11490d-01
        hbs( 1, 5, 2) = 0.62300d-02
        hbs( 1, 6, 2) = 0.46070d-01
        hbs( 1, 7, 2) = -0.43000d-03
        hbs( 1, 8, 2) = 0.44740d-01
        hbs( 1, 9, 2) = 0.80000d-04
        hbs( 1,10, 2) = -0.10140d-01
        hbs( 1,11, 2) = -0.44000d-03
        hbs( 1,12, 2) = 0.53800d-02
        hbs( 2, 1, 2) = 0.45070d-01
        hbs( 2, 2, 2) = -0.50600d-02
        hbs( 2, 3, 2) = 0.16200d-01
        hbs( 2, 4, 2) = -0.18400d-02
        hbs( 2, 5, 2) = -0.27200d-02
        hbs( 2, 6, 2) = -0.20110d-01
        hbs( 2, 7, 2) = 0.11400d-02
        hbs( 2, 8, 2) = -0.26860d-01
        hbs( 2, 9, 2) = -0.17000d-03
        hbs( 2,10, 2) = 0.11580d-01
        hbs( 2,11, 2) = -0.64000d-03
        hbs( 2,12, 2) = -0.16600d-02
        hbs( 3, 1, 2) = -0.79200d-02
        hbs( 3, 2, 2) = 0.40140d-01
        hbs( 3, 3, 2) = 0.21620d-01
        hbs( 3, 4, 2) = 0.13330d-01
        hbs( 3, 5, 2) = -0.35100d-02
        hbs( 3, 6, 2) = -0.25960d-01
        hbs( 3, 7, 2) = -0.71000d-03
        hbs( 3, 8, 2) = -0.17880d-01
        hbs( 3, 9, 2) = 0.10000d-03
        hbs( 3,10, 2) = -0.14400d-02
        hbs( 3,11, 2) = 0.10800d-02
        hbs( 3,12, 2) = -0.37200d-02
        hbs( 4, 1, 2) = -0.29410d-01
        hbs( 4, 2, 2) = -0.87500d-02
        hbs( 4, 3, 2) = 0.18410d-01
        hbs( 4, 4, 2) = -0.31300d-01
        hbs( 4, 5, 2) = -0.10370d-01
        hbs( 4, 6, 2) = -0.13280d-01
        hbs( 4, 7, 2) = 0.49300d-02
        hbs( 4, 8, 2) = 0.44510d-01
        hbs( 4, 9, 2) = -0.65000d-03
        hbs( 4,10, 2) = 0.43310d-01
        hbs( 4,11, 2) = 0.18000d-03
        hbs( 4,12, 2) = -0.83000d-03
        hbs( 5, 1, 2) = 0.28150d-01
        hbs( 5, 2, 2) = 0.54900d-02
        hbs( 5, 3, 2) = -0.13780d-01
        hbs( 5, 4, 2) = 0.21700d-02
        hbs( 5, 5, 2) = 0.70200d-02
        hbs( 5, 6, 2) = 0.60300d-02
        hbs( 5, 7, 2) = -0.19500d-02
        hbs( 5, 8, 2) = -0.19880d-01
        hbs( 5, 9, 2) = -0.33000d-03
        hbs( 5,10, 2) = -0.24160d-01
        hbs( 5,11, 2) = -0.12000d-03
        hbs( 5,12, 2) = 0.77000d-03
        hbs( 6, 1, 2) = 0.12600d-02
        hbs( 6, 2, 2) = 0.32500d-02
        hbs( 6, 3, 2) = -0.46300d-02
        hbs( 6, 4, 2) = 0.29130d-01
        hbs( 6, 5, 2) = 0.33500d-02
        hbs( 6, 6, 2) = 0.72400d-02
        hbs( 6, 7, 2) = -0.29800d-02
        hbs( 6, 8, 2) = -0.24630d-01
        hbs( 6, 9, 2) = 0.99000d-03
        hbs( 6,10, 2) = -0.19140d-01
        hbs( 6,11, 2) = -0.70000d-04
        hbs( 6,12, 2) = 0.60000d-04
        hbs( 7, 1, 2) = -0.11820d-01
        hbs( 7, 2, 2) = 0.41350d-01
        hbs( 7, 3, 2) = 0.31600d-02
        hbs( 7, 4, 2) = -0.32360d-01
        hbs( 7, 5, 2) = 0.57200d-02
        hbs( 7, 6, 2) = -0.33140d-01
        hbs( 7, 7, 2) = -0.90100d-02
        hbs( 7, 8, 2) = -0.12310d-01
        hbs( 7, 9, 2) = 0.50300d-02
        hbs( 7,10, 2) = 0.46100d-01
        hbs( 7,11, 2) = -0.96000d-03
        hbs( 7,12, 2) = -0.70000d-04
        hbs( 8, 1, 2) = 0.65000d-02
        hbs( 8, 2, 2) = -0.22990d-01
        hbs( 8, 3, 2) = -0.46000d-03
        hbs( 8, 4, 2) = 0.29900d-01
        hbs( 8, 5, 2) = -0.34100d-02
        hbs( 8, 6, 2) = 0.38600d-02
        hbs( 8, 7, 2) = 0.47600d-02
        hbs( 8, 8, 2) = 0.53400d-02
        hbs( 8, 9, 2) = -0.20400d-02
        hbs( 8,10, 2) = -0.20620d-01
        hbs( 8,11, 2) = -0.19000d-03
        hbs( 8,12, 2) = 0.30000d-04
        hbs( 9, 1, 2) = 0.53200d-02
        hbs( 9, 2, 2) = -0.18360d-01
        hbs( 9, 3, 2) = -0.27000d-02
        hbs( 9, 4, 2) = 0.24600d-02
        hbs( 9, 5, 2) = -0.23100d-02
        hbs( 9, 6, 2) = 0.29280d-01
        hbs( 9, 7, 2) = 0.42500d-02
        hbs( 9, 8, 2) = 0.69700d-02
        hbs( 9, 9, 2) = -0.29900d-02
        hbs( 9,10, 2) = -0.25480d-01
        hbs( 9,11, 2) = 0.11500d-02
        hbs( 9,12, 2) = 0.40000d-04
        hbs(10, 1, 2) = 0.43830d-01
        hbs(10, 2, 2) = 0.43000d-01
        hbs(10, 3, 2) = -0.29200d-02
        hbs(10, 4, 2) = -0.13560d-01
        hbs(10, 5, 2) = -0.60000d-03
        hbs(10, 6, 2) = -0.32140d-01
        hbs(10, 7, 2) = 0.53900d-02
        hbs(10, 8, 2) = -0.31440d-01
        hbs(10, 9, 2) = -0.88000d-02
        hbs(10,10, 2) = -0.13100d-01
        hbs(10,11, 2) = 0.53000d-02
        hbs(10,12, 2) = -0.54000d-03
        hbs(11, 1, 2) = -0.24340d-01
        hbs(11, 2, 2) = -0.19110d-01
        hbs(11, 3, 2) = 0.11900d-02
        hbs(11, 4, 2) = 0.71000d-02
        hbs(11, 5, 2) = 0.10500d-02
        hbs(11, 6, 2) = 0.30020d-01
        hbs(11, 7, 2) = -0.32200d-02
        hbs(11, 8, 2) = 0.19900d-02
        hbs(11, 9, 2) = 0.44100d-02
        hbs(11,10, 2) = 0.61200d-02
        hbs(11,11, 2) = -0.21600d-02
        hbs(11,12, 2) = -0.33000d-03
        hbs(12, 1, 2) = -0.19490d-01
        hbs(12, 2, 2) = -0.23890d-01
        hbs(12, 3, 2) = 0.17300d-02
        hbs(12, 4, 2) = 0.64600d-02
        hbs(12, 5, 2) = -0.45000d-03
        hbs(12, 6, 2) = 0.21300d-02
        hbs(12, 7, 2) = -0.21700d-02
        hbs(12, 8, 2) = 0.29450d-01
        hbs(12, 9, 2) = 0.43900d-02
        hbs(12,10, 2) = 0.69800d-02
        hbs(12,11, 2) = -0.31500d-02
        hbs(12,12, 2) = 0.87000d-03
        hbs(13, 1, 2) = 0.42040d-01
        hbs(13, 2, 2) = -0.12440d-01
        hbs(13, 3, 2) = 0.58600d-02
        hbs(13, 4, 2) = 0.45620d-01
        hbs(13, 5, 2) = -0.90000d-04
        hbs(13, 6, 2) = -0.12560d-01
        hbs(13, 7, 2) = -0.10300d-02
        hbs(13, 8, 2) = -0.32810d-01
        hbs(13, 9, 2) = 0.50400d-02
        hbs(13,10, 2) = -0.33390d-01
        hbs(13,11, 2) = -0.90100d-02
        hbs(13,12, 2) = 0.50300d-02
        hbs(14, 1, 2) = -0.18640d-01
        hbs(14, 2, 2) = 0.53900d-02
        hbs(14, 3, 2) = -0.35400d-02
        hbs(14, 4, 2) = -0.25230d-01
        hbs(14, 5, 2) = 0.70000d-04
        hbs(14, 6, 2) = 0.71300d-02
        hbs(14, 7, 2) = 0.12100d-02
        hbs(14, 8, 2) = 0.29250d-01
        hbs(14, 9, 2) = -0.29800d-02
        hbs(14,10, 2) = 0.28500d-02
        hbs(14,11, 2) = 0.44100d-02
        hbs(14,12, 2) = -0.19500d-02
        hbs(15, 1, 2) = -0.23400d-01
        hbs(15, 2, 2) = 0.70500d-02
        hbs(15, 3, 2) = -0.23200d-02
        hbs(15, 4, 2) = -0.20390d-01
        hbs(15, 5, 2) = 0.20000d-04
        hbs(15, 6, 2) = 0.54300d-02
        hbs(15, 7, 2) = -0.17000d-03
        hbs(15, 8, 2) = 0.35600d-02
        hbs(15, 9, 2) = -0.20600d-02
        hbs(15,10, 2) = 0.30550d-01
        hbs(15,11, 2) = 0.46000d-02
        hbs(15,12, 2) = -0.30800d-02
        hbs(16, 1, 2) = -0.75000d-02
        hbs(16, 2, 2) = -0.28090d-01
        hbs(16, 3, 2) = 0.13300d-01
        hbs(16, 4, 2) = 0.43100d-01
        hbs(16, 5, 2) = -0.88000d-03
        hbs(16, 6, 2) = 0.45050d-01
        hbs(16, 7, 2) = 0.16000d-03
        hbs(16, 8, 2) = -0.12690d-01
        hbs(16, 9, 2) = -0.70000d-03
        hbs(16,10, 2) = -0.32780d-01
        hbs(16,11, 2) = 0.49200d-02
        hbs(16,12, 2) = -0.89700d-02
        hbs(17, 1, 2) = 0.38600d-02
        hbs(17, 2, 2) = 0.25790d-01
        hbs(17, 3, 2) = -0.82600d-02
        hbs(17, 4, 2) = -0.23900d-01
        hbs(17, 5, 2) = 0.12400d-02
        hbs(17, 6, 2) = -0.19910d-01
        hbs(17, 7, 2) = -0.90000d-04
        hbs(17, 8, 2) = 0.61400d-02
        hbs(17, 9, 2) = -0.27000d-03
        hbs(17,10, 2) = 0.32800d-02
        hbs(17,11, 2) = -0.18800d-02
        hbs(17,12, 2) = 0.70400d-02
        hbs(18, 1, 2) = 0.36300d-02
        hbs(18, 2, 2) = 0.22900d-02
        hbs(18, 3, 2) = -0.50400d-02
        hbs(18, 4, 2) = -0.19200d-01
        hbs(18, 5, 2) = -0.37000d-03
        hbs(18, 6, 2) = -0.25140d-01
        hbs(18, 7, 2) = -0.80000d-04
        hbs(18, 8, 2) = 0.65600d-02
        hbs(18, 9, 2) = 0.97000d-03
        hbs(18,10, 2) = 0.29500d-01
        hbs(18,11, 2) = -0.30300d-02
        hbs(18,12, 2) = 0.19300d-02
        ! Hessian elements for anchor point 3
        hss( 1, 1, 3) = 0.39742d+00
        hss( 2, 1, 3) = 0.15810d-01
        hss( 2, 2, 3) = 0.37782d+00
        hss( 3, 1, 3) = 0.44160d-01
        hss( 3, 2, 3) = 0.54830d-01
        hss( 3, 3, 3) = 0.38872d+00
        hss( 4, 1, 3) = 0.19410d-01
        hss( 4, 2, 3) = -0.16000d-02
        hss( 4, 3, 3) = -0.34300d-02
        hss( 4, 4, 3) = 0.41358d+00
        hss( 5, 1, 3) = 0.23700d-02
        hss( 5, 2, 3) = -0.70000d-03
        hss( 5, 3, 3) = 0.14900d-02
        hss( 5, 4, 3) = 0.48900d-02
        hss( 5, 5, 3) = 0.36213d+00
        hss( 6, 1, 3) = -0.43800d-02
        hss( 6, 2, 3) = 0.80920d-01
        hss( 6, 3, 3) = -0.20700d-02
        hss( 6, 4, 3) = 0.21020d-01
        hss( 6, 5, 3) = -0.70000d-03
        hss( 6, 6, 3) = 0.41580d+00
        hss( 7, 1, 3) = -0.56000d-03
        hss( 7, 2, 3) = -0.41000d-03
        hss( 7, 3, 3) = 0.64000d-03
        hss( 7, 4, 3) = 0.44200d-02
        hss( 7, 5, 3) = 0.38000d-03
        hss( 7, 6, 3) = 0.45800d-02
        hss( 7, 7, 3) = 0.36425d+00
        hss( 8, 1, 3) = 0.80680d-01
        hss( 8, 2, 3) = -0.46300d-02
        hss( 8, 3, 3) = -0.10400d-02
        hss( 8, 4, 3) = -0.39500d-02
        hss( 8, 5, 3) = 0.11000d-03
        hss( 8, 6, 3) = 0.21280d-01
        hss( 8, 7, 3) = -0.32000d-03
        hss( 8, 8, 3) = 0.41002d+00
        hss( 9, 1, 3) = -0.40000d-03
        hss( 9, 2, 3) = -0.40000d-03
        hss( 9, 3, 3) = -0.55000d-03
        hss( 9, 4, 3) = -0.66000d-03
        hss( 9, 5, 3) = 0.17000d-03
        hss( 9, 6, 3) = 0.44100d-02
        hss( 9, 7, 3) = 0.44000d-03
        hss( 9, 8, 3) = 0.43500d-02
        hss( 9, 9, 3) = 0.36777d+00
        hss(10, 1, 3) = -0.11100d-02
        hss(10, 2, 3) = 0.19640d-01
        hss(10, 3, 3) = -0.29200d-02
        hss(10, 4, 3) = 0.86160d-01
        hss(10, 5, 3) = -0.75000d-03
        hss(10, 6, 3) = -0.35500d-02
        hss(10, 7, 3) = -0.56000d-03
        hss(10, 8, 3) = 0.21260d-01
        hss(10, 9, 3) = -0.57000d-03
        hss(10,10, 3) = 0.41917d+00
        hss(11, 1, 3) = -0.46000d-03
        hss(11, 2, 3) = -0.58000d-03
        hss(11, 3, 3) = 0.53000d-03
        hss(11, 4, 3) = -0.55000d-03
        hss(11, 5, 3) = 0.90000d-04
        hss(11, 6, 3) = -0.25000d-03
        hss(11, 7, 3) = 0.21000d-03
        hss(11, 8, 3) = 0.45200d-02
        hss(11, 9, 3) = 0.43000d-03
        hss(11,10, 3) = 0.45100d-02
        hss(11,11, 3) = 0.36393d+00
        hss(12, 1, 3) = 0.20000d-04
        hss(12, 2, 3) = 0.49500d-02
        hss(12, 3, 3) = -0.25500d-02
        hss(12, 4, 3) = -0.86000d-03
        hss(12, 5, 3) = 0.25000d-03
        hss(12, 6, 3) = -0.25000d-03
        hss(12, 7, 3) = 0.90000d-04
        hss(12, 8, 3) = -0.82000d-03
        hss(12, 9, 3) = 0.13000d-03
        hss(12,10, 3) = 0.39900d-02
        hss(12,11, 3) = 0.47000d-03
        hss(12,12, 3) = 0.36874d+00
        ! Hessian elements for anchor point 3
        hbb( 1, 1, 3) = 0.82270d-01
        hbb( 2, 1, 3) = -0.48250d-01
        hbb( 2, 2, 3) = 0.13087d+00
        hbb( 3, 1, 3) = -0.34020d-01
        hbb( 3, 2, 3) = -0.82620d-01
        hbb( 3, 3, 3) = 0.11664d+00
        hbb( 4, 1, 3) = -0.33450d-01
        hbb( 4, 2, 3) = 0.20400d-01
        hbb( 4, 3, 3) = 0.13040d-01
        hbb( 4, 4, 3) = 0.77460d-01
        hbb( 5, 1, 3) = 0.17380d-01
        hbb( 5, 2, 3) = -0.84000d-02
        hbb( 5, 3, 3) = -0.89800d-02
        hbb( 5, 4, 3) = -0.37880d-01
        hbb( 5, 5, 3) = 0.75680d-01
        hbb( 6, 1, 3) = 0.16060d-01
        hbb( 6, 2, 3) = -0.12000d-01
        hbb( 6, 3, 3) = -0.40600d-02
        hbb( 6, 4, 3) = -0.39580d-01
        hbb( 6, 5, 3) = -0.37810d-01
        hbb( 6, 6, 3) = 0.77390d-01
        hbb( 7, 1, 3) = -0.10450d-01
        hbb( 7, 2, 3) = 0.14900d-02
        hbb( 7, 3, 3) = 0.89500d-02
        hbb( 7, 4, 3) = -0.31240d-01
        hbb( 7, 5, 3) = 0.12780d-01
        hbb( 7, 6, 3) = 0.18460d-01
        hbb( 7, 7, 3) = 0.78900d-01
        hbb( 8, 1, 3) = -0.47000d-03
        hbb( 8, 2, 3) = 0.30400d-02
        hbb( 8, 3, 3) = -0.25700d-02
        hbb( 8, 4, 3) = 0.18370d-01
        hbb( 8, 5, 3) = -0.69100d-02
        hbb( 8, 6, 3) = -0.11460d-01
        hbb( 8, 7, 3) = -0.39290d-01
        hbb( 8, 8, 3) = 0.78920d-01
        hbb( 9, 1, 3) = 0.10920d-01
        hbb( 9, 2, 3) = -0.45400d-02
        hbb( 9, 3, 3) = -0.63800d-02
        hbb( 9, 4, 3) = 0.12870d-01
        hbb( 9, 5, 3) = -0.58800d-02
        hbb( 9, 6, 3) = -0.70000d-02
        hbb( 9, 7, 3) = -0.39610d-01
        hbb( 9, 8, 3) = -0.39630d-01
        hbb( 9, 9, 3) = 0.79240d-01
        hbb(10, 1, 3) = 0.65500d-02
        hbb(10, 2, 3) = -0.46300d-02
        hbb(10, 3, 3) = -0.19200d-02
        hbb(10, 4, 3) = -0.92300d-02
        hbb(10, 5, 3) = 0.10020d-01
        hbb(10, 6, 3) = -0.79000d-03
        hbb(10, 7, 3) = -0.34090d-01
        hbb(10, 8, 3) = 0.13830d-01
        hbb(10, 9, 3) = 0.20260d-01
        hbb(10,10, 3) = 0.80030d-01
        hbb(11, 1, 3) = -0.33400d-02
        hbb(11, 2, 3) = 0.27100d-02
        hbb(11, 3, 3) = 0.63000d-03
        hbb(11, 4, 3) = -0.98000d-03
        hbb(11, 5, 3) = -0.35000d-02
        hbb(11, 6, 3) = 0.44800d-02
        hbb(11, 7, 3) = 0.19980d-01
        hbb(11, 8, 3) = -0.73600d-02
        hbb(11, 9, 3) = -0.12620d-01
        hbb(11,10, 3) = -0.40030d-01
        hbb(11,11, 3) = 0.78930d-01
        hbb(12, 1, 3) = -0.32100d-02
        hbb(12, 2, 3) = 0.19200d-02
        hbb(12, 3, 3) = 0.13000d-02
        hbb(12, 4, 3) = 0.10210d-01
        hbb(12, 5, 3) = -0.65300d-02
        hbb(12, 6, 3) = -0.36800d-02
        hbb(12, 7, 3) = 0.14110d-01
        hbb(12, 8, 3) = -0.64600d-02
        hbb(12, 9, 3) = -0.76400d-02
        hbb(12,10, 3) = -0.40000d-01
        hbb(12,11, 3) = -0.38900d-01
        hbb(12,12, 3) = 0.78900d-01
        hbb(13, 1, 3) = -0.96900d-02
        hbb(13, 2, 3) = 0.10920d-01
        hbb(13, 3, 3) = -0.12400d-02
        hbb(13, 4, 3) = 0.49700d-02
        hbb(13, 5, 3) = -0.18200d-02
        hbb(13, 6, 3) = -0.31500d-02
        hbb(13, 7, 3) = -0.88100d-02
        hbb(13, 8, 3) = 0.10200d-01
        hbb(13, 9, 3) = -0.13900d-02
        hbb(13,10, 3) = -0.34270d-01
        hbb(13,11, 3) = 0.14180d-01
        hbb(13,12, 3) = 0.20100d-01
        hbb(13,13, 3) = 0.79280d-01
        hbb(14, 1, 3) = 0.10380d-01
        hbb(14, 2, 3) = -0.74600d-02
        hbb(14, 3, 3) = -0.29200d-02
        hbb(14, 4, 3) = -0.26400d-02
        hbb(14, 5, 3) = 0.92000d-03
        hbb(14, 6, 3) = 0.17200d-02
        hbb(14, 7, 3) = -0.14500d-02
        hbb(14, 8, 3) = -0.35600d-02
        hbb(14, 9, 3) = 0.50100d-02
        hbb(14,10, 3) = 0.20210d-01
        hbb(14,11, 3) = -0.76400d-02
        hbb(14,12, 3) = -0.12580d-01
        hbb(14,13, 3) = -0.39700d-01
        hbb(14,14, 3) = 0.79390d-01
        hbb(15, 1, 3) = -0.69000d-03
        hbb(15, 2, 3) = -0.34600d-02
        hbb(15, 3, 3) = 0.41500d-02
        hbb(15, 4, 3) = -0.23300d-02
        hbb(15, 5, 3) = 0.91000d-03
        hbb(15, 6, 3) = 0.14300d-02
        hbb(15, 7, 3) = 0.10260d-01
        hbb(15, 8, 3) = -0.66500d-02
        hbb(15, 9, 3) = -0.36100d-02
        hbb(15,10, 3) = 0.14060d-01
        hbb(15,11, 3) = -0.65400d-02
        hbb(15,12, 3) = -0.75200d-02
        hbb(15,13, 3) = -0.39580d-01
        hbb(15,14, 3) = -0.39690d-01
        hbb(15,15, 3) = 0.79270d-01
        hbb(16, 1, 3) = -0.35240d-01
        hbb(16, 2, 3) = 0.20060d-01
        hbb(16, 3, 3) = 0.15180d-01
        hbb(16, 4, 3) = -0.85100d-02
        hbb(16, 5, 3) = -0.49000d-03
        hbb(16, 6, 3) = 0.90000d-02
        hbb(16, 7, 3) = 0.56800d-02
        hbb(16, 8, 3) = -0.26300d-02
        hbb(16, 9, 3) = -0.30500d-02
        hbb(16,10, 3) = -0.90000d-02
        hbb(16,11, 3) = 0.10190d-01
        hbb(16,12, 3) = -0.12000d-02
        hbb(16,13, 3) = -0.31490d-01
        hbb(16,14, 3) = 0.13200d-01
        hbb(16,15, 3) = 0.18290d-01
        hbb(16,16, 3) = 0.78560d-01
        hbb(17, 1, 3) = 0.20310d-01
        hbb(17, 2, 3) = -0.11130d-01
        hbb(17, 3, 3) = -0.91900d-02
        hbb(17, 4, 3) = -0.11900d-02
        hbb(17, 5, 3) = 0.41100d-02
        hbb(17, 6, 3) = -0.29200d-02
        hbb(17, 7, 3) = -0.26100d-02
        hbb(17, 8, 3) = 0.13200d-02
        hbb(17, 9, 3) = 0.12900d-02
        hbb(17,10, 3) = 0.10050d-01
        hbb(17,11, 3) = -0.66300d-02
        hbb(17,12, 3) = -0.34200d-02
        hbb(17,13, 3) = 0.12780d-01
        hbb(17,14, 3) = -0.59900d-02
        hbb(17,15, 3) = -0.67900d-02
        hbb(17,16, 3) = -0.39350d-01
        hbb(17,17, 3) = 0.77640d-01
        hbb(18, 1, 3) = 0.14930d-01
        hbb(18, 2, 3) = -0.89300d-02
        hbb(18, 3, 3) = -0.60000d-02
        hbb(18, 4, 3) = 0.97000d-02
        hbb(18, 5, 3) = -0.36200d-02
        hbb(18, 6, 3) = -0.60800d-02
        hbb(18, 7, 3) = -0.30800d-02
        hbb(18, 8, 3) = 0.13100d-02
        hbb(18, 9, 3) = 0.17700d-02
        hbb(18,10, 3) = -0.10500d-02
        hbb(18,11, 3) = -0.35600d-02
        hbb(18,12, 3) = 0.46100d-02
        hbb(18,13, 3) = 0.18710d-01
        hbb(18,14, 3) = -0.72100d-02
        hbb(18,15, 3) = -0.11500d-01
        hbb(18,16, 3) = -0.39210d-01
        hbb(18,17, 3) = -0.38290d-01
        hbb(18,18, 3) = 0.77500d-01
        ! Hessian elements for anchor point 3
        htt( 1, 1, 3) = 0.10290d-01
        htt( 2, 1, 3) = 0.63900d-02
        htt( 2, 2, 3) = 0.13000d-01
        htt( 3, 1, 3) = 0.28000d-03
        htt( 3, 2, 3) = -0.46200d-02
        htt( 3, 3, 3) = 0.11260d-01
        htt( 4, 1, 3) = -0.36100d-02
        htt( 4, 2, 3) = 0.19900d-02
        htt( 4, 3, 3) = 0.63700d-02
        htt( 4, 4, 3) = 0.11970d-01
        htt( 5, 1, 3) = -0.93800d-02
        htt( 5, 2, 3) = -0.10580d-01
        htt( 5, 3, 3) = 0.29200d-02
        htt( 5, 4, 3) = 0.17200d-02
        htt( 5, 5, 3) = 0.16730d-01
        htt( 6, 1, 3) = -0.54800d-02
        htt( 6, 2, 3) = -0.54200d-02
        htt( 6, 3, 3) = 0.33200d-02
        htt( 6, 4, 3) = 0.33700d-02
        htt( 6, 5, 3) = 0.66800d-02
        htt( 6, 6, 3) = 0.91600d-02
        htt( 7, 1, 3) = -0.40000d-04
        htt( 7, 2, 3) = -0.31000d-03
        htt( 7, 3, 3) = -0.73300d-02
        htt( 7, 4, 3) = -0.75900d-02
        htt( 7, 5, 3) = 0.52600d-02
        htt( 7, 6, 3) = -0.15300d-02
        htt( 7, 7, 3) = 0.12050d-01
        htt( 8, 1, 3) = 0.38600d-02
        htt( 8, 2, 3) = 0.48500d-02
        htt( 8, 3, 3) = -0.69300d-02
        htt( 8, 4, 3) = -0.59400d-02
        htt( 8, 5, 3) = -0.48000d-02
        htt( 8, 6, 3) = 0.95000d-03
        htt( 8, 7, 3) = 0.52700d-02
        htt( 8, 8, 3) = 0.11020d-01
        htt( 9, 1, 3) = -0.52600d-02
        htt( 9, 2, 3) = 0.85000d-03
        htt( 9, 3, 3) = -0.32600d-02
        htt( 9, 4, 3) = 0.28400d-02
        htt( 9, 5, 3) = -0.16300d-02
        htt( 9, 6, 3) = 0.11100d-02
        htt( 9, 7, 3) = -0.35000d-02
        htt( 9, 8, 3) = -0.76000d-03
        htt( 9, 9, 3) = 0.98000d-02
        htt(10, 1, 3) = -0.32000d-02
        htt(10, 2, 3) = 0.20800d-02
        htt(10, 3, 3) = -0.19100d-02
        htt(10, 4, 3) = 0.33700d-02
        htt(10, 5, 3) = 0.37000d-03
        htt(10, 6, 3) = 0.79000d-03
        htt(10, 7, 3) = -0.83000d-03
        htt(10, 8, 3) = -0.41000d-03
        htt(10, 9, 3) = 0.34100d-02
        htt(10,10, 3) = 0.86300d-02
        htt(11, 1, 3) = -0.12800d-02
        htt(11, 2, 3) = -0.59200d-02
        htt(11, 3, 3) = 0.17500d-02
        htt(11, 4, 3) = -0.28900d-02
        htt(11, 5, 3) = -0.41000d-03
        htt(11, 6, 3) = 0.10500d-02
        htt(11, 7, 3) = -0.32300d-02
        htt(11, 8, 3) = -0.17700d-02
        htt(11, 9, 3) = 0.35500d-02
        htt(11,10, 3) = -0.20000d-02
        htt(11,11, 3) = 0.83000d-02
        htt(12, 1, 3) = 0.78000d-03
        htt(12, 2, 3) = -0.46800d-02
        htt(12, 3, 3) = 0.31000d-02
        htt(12, 4, 3) = -0.23600d-02
        htt(12, 5, 3) = 0.16000d-02
        htt(12, 6, 3) = 0.73000d-03
        htt(12, 7, 3) = -0.56000d-03
        htt(12, 8, 3) = -0.14300d-02
        htt(12, 9, 3) = -0.28400d-02
        htt(12,10, 3) = 0.32200d-02
        htt(12,11, 3) = 0.27500d-02
        htt(12,12, 3) = 0.88100d-02
        htt(13, 1, 3) = -0.66000d-03
        htt(13, 2, 3) = -0.37900d-02
        htt(13, 3, 3) = 0.29600d-02
        htt(13, 4, 3) = -0.16000d-03
        htt(13, 5, 3) = 0.51200d-02
        htt(13, 6, 3) = 0.20100d-02
        htt(13, 7, 3) = 0.17400d-02
        htt(13, 8, 3) = -0.13700d-02
        htt(13, 9, 3) = -0.72600d-02
        htt(13,10, 3) = -0.74000d-03
        htt(13,11, 3) = -0.40600d-02
        htt(13,12, 3) = 0.24600d-02
        htt(13,13, 3) = 0.10500d-01
        htt(14, 1, 3) = 0.16400d-02
        htt(14, 2, 3) = -0.60000d-04
        htt(14, 3, 3) = 0.14800d-02
        htt(14, 4, 3) = -0.22000d-03
        htt(14, 5, 3) = 0.17400d-02
        htt(14, 6, 3) = -0.85000d-03
        htt(14, 7, 3) = 0.18900d-02
        htt(14, 8, 3) = -0.70000d-03
        htt(14, 9, 3) = -0.45800d-02
        htt(14,10, 3) = 0.11000d-02
        htt(14,11, 3) = -0.28400d-02
        htt(14,12, 3) = 0.28400d-02
        htt(14,13, 3) = 0.40400d-02
        htt(14,14, 3) = 0.85900d-02
        htt(15, 1, 3) = -0.27400d-02
        htt(15, 2, 3) = -0.50400d-02
        htt(15, 3, 3) = 0.16000d-02
        htt(15, 4, 3) = -0.70000d-03
        htt(15, 5, 3) = 0.30800d-02
        htt(15, 6, 3) = 0.23300d-02
        htt(15, 7, 3) = -0.97000d-03
        htt(15, 8, 3) = -0.17200d-02
        htt(15, 9, 3) = -0.78000d-03
        htt(15,10, 3) = -0.60300d-02
        htt(15,11, 3) = 0.15700d-02
        htt(15,12, 3) = -0.36800d-02
        htt(15,13, 3) = 0.39000d-02
        htt(15,14, 3) = -0.17100d-02
        htt(15,15, 3) = 0.92100d-02
        htt(16, 1, 3) = -0.44000d-03
        htt(16, 2, 3) = -0.13100d-02
        htt(16, 3, 3) = 0.11000d-03
        htt(16, 4, 3) = -0.76000d-03
        htt(16, 5, 3) = -0.29000d-03
        htt(16, 6, 3) = -0.53000d-03
        htt(16, 7, 3) = -0.81000d-03
        htt(16, 8, 3) = -0.10500d-02
        htt(16, 9, 3) = 0.18900d-02
        htt(16,10, 3) = -0.41800d-02
        htt(16,11, 3) = 0.27800d-02
        htt(16,12, 3) = -0.33000d-02
        htt(16,13, 3) = -0.25500d-02
        htt(16,14, 3) = 0.28300d-02
        htt(16,15, 3) = 0.36000d-02
        htt(16,16, 3) = 0.89900d-02
        htt(17, 1, 3) = 0.15600d-02
        htt(17, 2, 3) = -0.56000d-03
        htt(17, 3, 3) = 0.34000d-03
        htt(17, 4, 3) = -0.17800d-02
        htt(17, 5, 3) = 0.25100d-02
        htt(17, 6, 3) = -0.75000d-03
        htt(17, 7, 3) = 0.36600d-02
        htt(17, 8, 3) = 0.39000d-03
        htt(17, 9, 3) = 0.16000d-03
        htt(17,10, 3) = -0.21700d-02
        htt(17,11, 3) = 0.23300d-02
        htt(17,12, 3) = 0.00000d+00
        htt(17,13, 3) = -0.59100d-02
        htt(17,14, 3) = -0.56000d-03
        htt(17,15, 3) = -0.35500d-02
        htt(17,16, 3) = 0.18000d-02
        htt(17,17, 3) = 0.10140d-01
        htt(18, 1, 3) = 0.13300d-02
        htt(18, 2, 3) = 0.25800d-02
        htt(18, 3, 3) = -0.22700d-02
        htt(18, 4, 3) = -0.10200d-02
        htt(18, 5, 3) = -0.30300d-02
        htt(18, 6, 3) = -0.35800d-02
        htt(18, 7, 3) = 0.33000d-03
        htt(18, 8, 3) = -0.23000d-03
        htt(18, 9, 3) = 0.24700d-02
        htt(18,10, 3) = 0.75000d-03
        htt(18,11, 3) = 0.12000d-02
        htt(18,12, 3) = -0.52000d-03
        htt(18,13, 3) = -0.44900d-02
        htt(18,14, 3) = 0.90000d-03
        htt(18,15, 3) = -0.27500d-02
        htt(18,16, 3) = 0.26400d-02
        htt(18,17, 3) = 0.27400d-02
        htt(18,18, 3) = 0.87100d-02
        htt(19, 1, 3) = -0.74000d-03
        htt(19, 2, 3) = -0.42800d-02
        htt(19, 3, 3) = 0.18200d-02
        htt(19, 4, 3) = -0.17200d-02
        htt(19, 5, 3) = 0.58900d-02
        htt(19, 6, 3) = 0.21100d-02
        htt(19, 7, 3) = 0.35000d-02
        htt(19, 8, 3) = -0.28000d-03
        htt(19, 9, 3) = -0.25100d-02
        htt(19,10, 3) = -0.40200d-02
        htt(19,11, 3) = 0.11200d-02
        htt(19,12, 3) = -0.38000d-03
        htt(19,13, 3) = 0.54000d-03
        htt(19,14, 3) = -0.51100d-02
        htt(19,15, 3) = 0.20600d-02
        htt(19,16, 3) = -0.35900d-02
        htt(19,17, 3) = 0.47900d-02
        htt(19,18, 3) = -0.26500d-02
        htt(19,19, 3) = 0.10440d-01
        htt(20, 1, 3) = -0.97000d-03
        htt(20, 2, 3) = -0.11500d-02
        htt(20, 3, 3) = -0.78000d-03
        htt(20, 4, 3) = -0.96000d-03
        htt(20, 5, 3) = 0.35000d-03
        htt(20, 6, 3) = -0.72000d-03
        htt(20, 7, 3) = 0.18000d-03
        htt(20, 8, 3) = -0.89000d-03
        htt(20, 9, 3) = -0.20000d-03
        htt(20,10, 3) = -0.10900d-02
        htt(20,11, 3) = -0.20000d-04
        htt(20,12, 3) = -0.90000d-03
        htt(20,13, 3) = 0.19600d-02
        htt(20,14, 3) = -0.36500d-02
        htt(20,15, 3) = 0.28600d-02
        htt(20,16, 3) = -0.27500d-02
        htt(20,17, 3) = -0.26100d-02
        htt(20,18, 3) = 0.33200d-02
        htt(20,19, 3) = 0.30100d-02
        htt(20,20, 3) = 0.89300d-02
        htt(21, 1, 3) = 0.34400d-02
        htt(21, 2, 3) = 0.76900d-02
        htt(21, 3, 3) = -0.32500d-02
        htt(21, 4, 3) = 0.10000d-02
        htt(21, 5, 3) = -0.13360d-01
        htt(21, 6, 3) = -0.35700d-02
        htt(21, 7, 3) = -0.71200d-02
        htt(21, 8, 3) = 0.26800d-02
        htt(21, 9, 3) = 0.42100d-02
        htt(21,10, 3) = 0.23400d-02
        htt(21,11, 3) = -0.13000d-03
        htt(21,12, 3) = -0.20000d-02
        htt(21,13, 3) = -0.18100d-02
        htt(21,14, 3) = -0.22900d-02
        htt(21,15, 3) = 0.90000d-04
        htt(21,16, 3) = -0.39000d-03
        htt(21,17, 3) = -0.84600d-02
        htt(21,18, 3) = 0.98000d-03
        htt(21,19, 3) = -0.79800d-02
        htt(21,20, 3) = 0.14600d-02
        htt(21,21, 3) = 0.15990d-01
        htt(22, 1, 3) = -0.55000d-03
        htt(22, 2, 3) = 0.24000d-02
        htt(22, 3, 3) = -0.36500d-02
        htt(22, 4, 3) = -0.70000d-03
        htt(22, 5, 3) = -0.30600d-02
        htt(22, 6, 3) = -0.61100d-02
        htt(22, 7, 3) = -0.17000d-03
        htt(22, 8, 3) = -0.32100d-02
        htt(22, 9, 3) = 0.14100d-02
        htt(22,10, 3) = 0.19200d-02
        htt(22,11, 3) = -0.16200d-02
        htt(22,12, 3) = -0.11100d-02
        htt(22,13, 3) = 0.13700d-02
        htt(22,14, 3) = 0.37000d-03
        htt(22,15, 3) = 0.86000d-03
        htt(22,16, 3) = -0.15000d-03
        htt(22,17, 3) = -0.51200d-02
        htt(22,18, 3) = 0.15500d-02
        htt(22,19, 3) = -0.41100d-02
        htt(22,20, 3) = 0.25600d-02
        htt(22,21, 3) = 0.59500d-02
        htt(22,22, 3) = 0.90700d-02
        htt(23, 1, 3) = 0.36700d-02
        htt(23, 2, 3) = 0.46000d-02
        htt(23, 3, 3) = -0.68000d-03
        htt(23, 4, 3) = 0.25000d-03
        htt(23, 5, 3) = -0.79000d-02
        htt(23, 6, 3) = -0.77000d-03
        htt(23, 7, 3) = -0.38400d-02
        htt(23, 8, 3) = 0.32800d-02
        htt(23, 9, 3) = 0.19400d-02
        htt(23,10, 3) = -0.54000d-03
        htt(23,11, 3) = 0.99000d-03
        htt(23,12, 3) = -0.14900d-02
        htt(23,13, 3) = -0.32100d-02
        htt(23,14, 3) = -0.37300d-02
        htt(23,15, 3) = -0.70000d-03
        htt(23,16, 3) = -0.12200d-02
        htt(23,17, 3) = -0.11700d-02
        htt(23,18, 3) = -0.49000d-02
        htt(23,19, 3) = -0.65000d-03
        htt(23,20, 3) = -0.43800d-02
        htt(23,21, 3) = 0.66800d-02
        htt(23,22, 3) = -0.62000d-03
        htt(23,23, 3) = 0.10360d-01
        htt(24, 1, 3) = -0.33000d-03
        htt(24, 2, 3) = -0.69000d-03
        htt(24, 3, 3) = -0.10800d-02
        htt(24, 4, 3) = -0.14500d-02
        htt(24, 5, 3) = 0.24000d-02
        htt(24, 6, 3) = -0.33100d-02
        htt(24, 7, 3) = 0.31100d-02
        htt(24, 8, 3) = -0.26100d-02
        htt(24, 9, 3) = -0.87000d-03
        htt(24,10, 3) = -0.97000d-03
        htt(24,11, 3) = -0.50000d-03
        htt(24,12, 3) = -0.60000d-03
        htt(24,13, 3) = -0.30000d-04
        htt(24,14, 3) = -0.10800d-02
        htt(24,15, 3) = 0.70000d-04
        htt(24,16, 3) = -0.98000d-03
        htt(24,17, 3) = 0.21800d-02
        htt(24,18, 3) = -0.43300d-02
        htt(24,19, 3) = 0.32200d-02
        htt(24,20, 3) = -0.32900d-02
        htt(24,21, 3) = -0.33600d-02
        htt(24,22, 3) = 0.25000d-02
        htt(24,23, 3) = 0.30700d-02
        htt(24,24, 3) = 0.89200d-02
        ! Hessian elements for anchor point 3
        hbs( 1, 1, 3) = -0.35620d-01
        hbs( 1, 2, 3) = -0.31870d-01
        hbs( 1, 3, 3) = -0.41940d-01
        hbs( 1, 4, 3) = -0.11410d-01
        hbs( 1, 5, 3) = 0.68300d-02
        hbs( 1, 6, 3) = 0.45490d-01
        hbs( 1, 7, 3) = -0.31000d-03
        hbs( 1, 8, 3) = 0.44930d-01
        hbs( 1, 9, 3) = 0.00000d+00
        hbs( 1,10, 3) = -0.10110d-01
        hbs( 1,11, 3) = -0.33000d-03
        hbs( 1,12, 3) = 0.43300d-02
        hbs( 2, 1, 3) = 0.49480d-01
        hbs( 2, 2, 3) = -0.69900d-02
        hbs( 2, 3, 3) = 0.24890d-01
        hbs( 2, 4, 3) = -0.87000d-03
        hbs( 2, 5, 3) = -0.66200d-02
        hbs( 2, 6, 3) = -0.20260d-01
        hbs( 2, 7, 3) = 0.89000d-03
        hbs( 2, 8, 3) = -0.28090d-01
        hbs( 2, 9, 3) = -0.10000d-04
        hbs( 2,10, 3) = 0.10450d-01
        hbs( 2,11, 3) = -0.46000d-03
        hbs( 2,12, 3) = -0.19500d-02
        hbs( 3, 1, 3) = -0.13860d-01
        hbs( 3, 2, 3) = 0.38860d-01
        hbs( 3, 3, 3) = 0.17050d-01
        hbs( 3, 4, 3) = 0.12280d-01
        hbs( 3, 5, 3) = -0.21000d-03
        hbs( 3, 6, 3) = -0.25230d-01
        hbs( 3, 7, 3) = -0.58000d-03
        hbs( 3, 8, 3) = -0.16840d-01
        hbs( 3, 9, 3) = 0.20000d-04
        hbs( 3,10, 3) = -0.34000d-03
        hbs( 3,11, 3) = 0.79000d-03
        hbs( 3,12, 3) = -0.23700d-02
        hbs( 4, 1, 3) = -0.28260d-01
        hbs( 4, 2, 3) = -0.99500d-02
        hbs( 4, 3, 3) = 0.16700d-01
        hbs( 4, 4, 3) = -0.32020d-01
        hbs( 4, 5, 3) = -0.96100d-02
        hbs( 4, 6, 3) = -0.12920d-01
        hbs( 4, 7, 3) = 0.48000d-02
        hbs( 4, 8, 3) = 0.45530d-01
        hbs( 4, 9, 3) = -0.62000d-03
        hbs( 4,10, 3) = 0.44050d-01
        hbs( 4,11, 3) = 0.10000d-03
        hbs( 4,12, 3) = -0.40000d-03
        hbs( 5, 1, 3) = 0.30480d-01
        hbs( 5, 2, 3) = 0.69800d-02
        hbs( 5, 3, 3) = -0.14860d-01
        hbs( 5, 4, 3) = 0.26800d-02
        hbs( 5, 5, 3) = 0.62000d-02
        hbs( 5, 6, 3) = 0.57700d-02
        hbs( 5, 7, 3) = -0.18500d-02
        hbs( 5, 8, 3) = -0.20690d-01
        hbs( 5, 9, 3) = -0.33000d-03
        hbs( 5,10, 3) = -0.24600d-01
        hbs( 5,11, 3) = -0.60000d-04
        hbs( 5,12, 3) = 0.44000d-03
        hbs( 6, 1, 3) = -0.22300d-02
        hbs( 6, 2, 3) = 0.29700d-02
        hbs( 6, 3, 3) = -0.18400d-02
        hbs( 6, 4, 3) = 0.29340d-01
        hbs( 6, 5, 3) = 0.34100d-02
        hbs( 6, 6, 3) = 0.71600d-02
        hbs( 6, 7, 3) = -0.29500d-02
        hbs( 6, 8, 3) = -0.24830d-01
        hbs( 6, 9, 3) = 0.94000d-03
        hbs( 6,10, 3) = -0.19440d-01
        hbs( 6,11, 3) = -0.40000d-04
        hbs( 6,12, 3) = -0.30000d-04
        hbs( 7, 1, 3) = -0.12340d-01
        hbs( 7, 2, 3) = 0.39640d-01
        hbs( 7, 3, 3) = 0.73200d-02
        hbs( 7, 4, 3) = -0.32120d-01
        hbs( 7, 5, 3) = 0.49300d-02
        hbs( 7, 6, 3) = -0.33370d-01
        hbs( 7, 7, 3) = -0.89700d-02
        hbs( 7, 8, 3) = -0.12760d-01
        hbs( 7, 9, 3) = 0.49400d-02
        hbs( 7,10, 3) = 0.45710d-01
        hbs( 7,11, 3) = -0.89000d-03
        hbs( 7,12, 3) = 0.50000d-04
        hbs( 8, 1, 3) = 0.73400d-02
        hbs( 8, 2, 3) = -0.22220d-01
        hbs( 8, 3, 3) = -0.35200d-02
        hbs( 8, 4, 3) = 0.29850d-01
        hbs( 8, 5, 3) = -0.30000d-02
        hbs( 8, 6, 3) = 0.38900d-02
        hbs( 8, 7, 3) = 0.48800d-02
        hbs( 8, 8, 3) = 0.56300d-02
        hbs( 8, 9, 3) = -0.19800d-02
        hbs( 8,10, 3) = -0.20400d-01
        hbs( 8,11, 3) = -0.22000d-03
        hbs( 8,12, 3) = -0.50000d-04
        hbs( 9, 1, 3) = 0.50100d-02
        hbs( 9, 2, 3) = -0.17420d-01
        hbs( 9, 3, 3) = -0.38000d-02
        hbs( 9, 4, 3) = 0.22700d-02
        hbs( 9, 5, 3) = -0.19200d-02
        hbs( 9, 6, 3) = 0.29480d-01
        hbs( 9, 7, 3) = 0.40900d-02
        hbs( 9, 8, 3) = 0.71400d-02
        hbs( 9, 9, 3) = -0.29600d-02
        hbs( 9,10, 3) = -0.25300d-01
        hbs( 9,11, 3) = 0.11100d-02
        hbs( 9,12, 3) = -0.10000d-04
        hbs(10, 1, 3) = 0.42740d-01
        hbs(10, 2, 3) = 0.42010d-01
        hbs(10, 3, 3) = -0.32900d-02
        hbs(10, 4, 3) = -0.13800d-01
        hbs(10, 5, 3) = -0.37000d-03
        hbs(10, 6, 3) = -0.32450d-01
        hbs(10, 7, 3) = 0.53400d-02
        hbs(10, 8, 3) = -0.31850d-01
        hbs(10, 9, 3) = -0.86200d-02
        hbs(10,10, 3) = -0.13500d-01
        hbs(10,11, 3) = 0.52600d-02
        hbs(10,12, 3) = -0.66000d-03
        hbs(11, 1, 3) = -0.23640d-01
        hbs(11, 2, 3) = -0.18640d-01
        hbs(11, 3, 3) = 0.15000d-02
        hbs(11, 4, 3) = 0.71800d-02
        hbs(11, 5, 3) = 0.84000d-03
        hbs(11, 6, 3) = 0.30210d-01
        hbs(11, 7, 3) = -0.32100d-02
        hbs(11, 8, 3) = 0.20600d-02
        hbs(11, 9, 3) = 0.42700d-02
        hbs(11,10, 3) = 0.63000d-02
        hbs(11,11, 3) = -0.21200d-02
        hbs(11,12, 3) = -0.27000d-03
        hbs(12, 1, 3) = -0.19100d-01
        hbs(12, 2, 3) = -0.23380d-01
        hbs(12, 3, 3) = 0.17900d-02
        hbs(12, 4, 3) = 0.66100d-02
        hbs(12, 5, 3) = -0.47000d-03
        hbs(12, 6, 3) = 0.22400d-02
        hbs(12, 7, 3) = -0.21200d-02
        hbs(12, 8, 3) = 0.29790d-01
        hbs(12, 9, 3) = 0.43500d-02
        hbs(12,10, 3) = 0.72000d-02
        hbs(12,11, 3) = -0.31400d-02
        hbs(12,12, 3) = 0.93000d-03
        hbs(13, 1, 3) = 0.40850d-01
        hbs(13, 2, 3) = -0.14500d-01
        hbs(13, 3, 3) = 0.10130d-01
        hbs(13, 4, 3) = 0.45630d-01
        hbs(13, 5, 3) = -0.20000d-03
        hbs(13, 6, 3) = -0.12730d-01
        hbs(13, 7, 3) = -0.96000d-03
        hbs(13, 8, 3) = -0.32840d-01
        hbs(13, 9, 3) = 0.49400d-02
        hbs(13,10, 3) = -0.32520d-01
        hbs(13,11, 3) = -0.89800d-02
        hbs(13,12, 3) = 0.51000d-02
        hbs(14, 1, 3) = -0.17990d-01
        hbs(14, 2, 3) = 0.57300d-02
        hbs(14, 3, 3) = -0.40700d-02
        hbs(14, 4, 3) = -0.25200d-01
        hbs(14, 5, 3) = 0.12000d-03
        hbs(14, 6, 3) = 0.71400d-02
        hbs(14, 7, 3) = 0.11600d-02
        hbs(14, 8, 3) = 0.29210d-01
        hbs(14, 9, 3) = -0.29300d-02
        hbs(14,10, 3) = 0.21300d-02
        hbs(14,11, 3) = 0.42700d-02
        hbs(14,12, 3) = -0.19500d-02
        hbs(15, 1, 3) = -0.22870d-01
        hbs(15, 2, 3) = 0.87800d-02
        hbs(15, 3, 3) = -0.60600d-02
        hbs(15, 4, 3) = -0.20420d-01
        hbs(15, 5, 3) = 0.70000d-04
        hbs(15, 6, 3) = 0.55900d-02
        hbs(15, 7, 3) = -0.20000d-03
        hbs(15, 8, 3) = 0.36300d-02
        hbs(15, 9, 3) = -0.20100d-02
        hbs(15,10, 3) = 0.30390d-01
        hbs(15,11, 3) = 0.47100d-02
        hbs(15,12, 3) = -0.31500d-02
        hbs(16, 1, 3) = -0.73800d-02
        hbs(16, 2, 3) = -0.25320d-01
        hbs(16, 3, 3) = 0.11080d-01
        hbs(16, 4, 3) = 0.43720d-01
        hbs(16, 5, 3) = -0.15800d-02
        hbs(16, 6, 3) = 0.45980d-01
        hbs(16, 7, 3) = 0.11000d-03
        hbs(16, 8, 3) = -0.13000d-01
        hbs(16, 9, 3) = -0.64000d-03
        hbs(16,10, 3) = -0.33620d-01
        hbs(16,11, 3) = 0.48400d-02
        hbs(16,12, 3) = -0.84200d-02
        hbs(17, 1, 3) = 0.40600d-02
        hbs(17, 2, 3) = 0.27060d-01
        hbs(17, 3, 3) = -0.92100d-02
        hbs(17, 4, 3) = -0.24320d-01
        hbs(17, 5, 3) = 0.15600d-02
        hbs(17, 6, 3) = -0.20640d-01
        hbs(17, 7, 3) = -0.50000d-04
        hbs(17, 8, 3) = 0.61500d-02
        hbs(17, 9, 3) = -0.30000d-03
        hbs(17,10, 3) = 0.36700d-02
        hbs(17,11, 3) = -0.18400d-02
        hbs(17,12, 3) = 0.60700d-02
        hbs(18, 1, 3) = 0.33200d-02
        hbs(18, 2, 3) = -0.17400d-02
        hbs(18, 3, 3) = -0.18800d-02
        hbs(18, 4, 3) = -0.19400d-01
        hbs(18, 5, 3) = 0.20000d-04
        hbs(18, 6, 3) = -0.25340d-01
        hbs(18, 7, 3) = -0.50000d-04
        hbs(18, 8, 3) = 0.68600d-02
        hbs(18, 9, 3) = 0.93000d-03
        hbs(18,10, 3) = 0.29960d-01
        hbs(18,11, 3) = -0.30000d-02
        hbs(18,12, 3) = 0.23500d-02
        ! Hessian elements for anchor point 4
        hss( 1, 1, 4) = 0.38571d+00
        hss( 2, 1, 4) = 0.12420d-01
        hss( 2, 2, 4) = 0.38571d+00
        hss( 3, 1, 4) = 0.51260d-01
        hss( 3, 2, 4) = 0.51260d-01
        hss( 3, 3, 4) = 0.38504d+00
        hss( 4, 1, 4) = 0.20570d-01
        hss( 4, 2, 4) = -0.17300d-02
        hss( 4, 3, 4) = -0.32200d-02
        hss( 4, 4, 4) = 0.41668d+00
        hss( 5, 1, 4) = 0.50700d-02
        hss( 5, 2, 4) = 0.48000d-03
        hss( 5, 3, 4) = -0.24900d-02
        hss( 5, 4, 4) = 0.39900d-02
        hss( 5, 5, 4) = 0.36905d+00
        hss( 6, 1, 4) = -0.36300d-02
        hss( 6, 2, 4) = 0.82190d-01
        hss( 6, 3, 4) = -0.13300d-02
        hss( 6, 4, 4) = 0.21260d-01
        hss( 6, 5, 4) = -0.91000d-03
        hss( 6, 6, 4) = 0.41243d+00
        hss( 7, 1, 4) = -0.51000d-03
        hss( 7, 2, 4) = -0.39000d-03
        hss( 7, 3, 4) = 0.41000d-03
        hss( 7, 4, 4) = 0.44700d-02
        hss( 7, 5, 4) = 0.43000d-03
        hss( 7, 6, 4) = 0.45000d-02
        hss( 7, 7, 4) = 0.36452d+00
        hss( 8, 1, 4) = 0.82190d-01
        hss( 8, 2, 4) = -0.36300d-02
        hss( 8, 3, 4) = -0.13300d-02
        hss( 8, 4, 4) = -0.41800d-02
        hss( 8, 5, 4) = -0.36000d-03
        hss( 8, 6, 4) = 0.20820d-01
        hss( 8, 7, 4) = -0.28000d-03
        hss( 8, 8, 4) = 0.41243d+00
        hss( 9, 1, 4) = -0.39000d-03
        hss( 9, 2, 4) = -0.39000d-03
        hss( 9, 3, 4) = -0.57000d-03
        hss( 9, 4, 4) = -0.62000d-03
        hss( 9, 5, 4) = 0.12000d-03
        hss( 9, 6, 4) = 0.43600d-02
        hss( 9, 7, 4) = 0.43000d-03
        hss( 9, 8, 4) = 0.43600d-02
        hss( 9, 9, 4) = 0.36801d+00
        hss(10, 1, 4) = -0.17300d-02
        hss(10, 2, 4) = 0.20570d-01
        hss(10, 3, 4) = -0.32200d-02
        hss(10, 4, 4) = 0.86500d-01
        hss(10, 5, 4) = -0.80000d-03
        hss(10, 6, 4) = -0.41800d-02
        hss(10, 7, 4) = -0.54000d-03
        hss(10, 8, 4) = 0.21260d-01
        hss(10, 9, 4) = -0.62000d-03
        hss(10,10, 4) = 0.41668d+00
        hss(11, 1, 4) = -0.39000d-03
        hss(11, 2, 4) = -0.51000d-03
        hss(11, 3, 4) = 0.41000d-03
        hss(11, 4, 4) = -0.54000d-03
        hss(11, 5, 4) = 0.90000d-04
        hss(11, 6, 4) = -0.28000d-03
        hss(11, 7, 4) = 0.20000d-03
        hss(11, 8, 4) = 0.45000d-02
        hss(11, 9, 4) = 0.43000d-03
        hss(11,10, 4) = 0.44700d-02
        hss(11,11, 4) = 0.36452d+00
        hss(12, 1, 4) = 0.48000d-03
        hss(12, 2, 4) = 0.50700d-02
        hss(12, 3, 4) = -0.24900d-02
        hss(12, 4, 4) = -0.80000d-03
        hss(12, 5, 4) = -0.30000d-04
        hss(12, 6, 4) = -0.36000d-03
        hss(12, 7, 4) = 0.90000d-04
        hss(12, 8, 4) = -0.91000d-03
        hss(12, 9, 4) = 0.12000d-03
        hss(12,10, 4) = 0.39900d-02
        hss(12,11, 4) = 0.43000d-03
        hss(12,12, 4) = 0.36905d+00
        ! Hessian elements for anchor point 4
        hbb( 1, 1, 4) = 0.80030d-01
        hbb( 2, 1, 4) = -0.40020d-01
        hbb( 2, 2, 4) = 0.92960d-01
        hbb( 3, 1, 4) = -0.40020d-01
        hbb( 3, 2, 4) = -0.52950d-01
        hbb( 3, 3, 4) = 0.92960d-01
        hbb( 4, 1, 4) = -0.33740d-01
        hbb( 4, 2, 4) = 0.17860d-01
        hbb( 4, 3, 4) = 0.15880d-01
        hbb( 4, 4, 4) = 0.77910d-01
        hbb( 5, 1, 4) = 0.19460d-01
        hbb( 5, 2, 4) = -0.10850d-01
        hbb( 5, 3, 4) = -0.86100d-02
        hbb( 5, 4, 4) = -0.38850d-01
        hbb( 5, 5, 4) = 0.77000d-01
        hbb( 6, 1, 4) = 0.14280d-01
        hbb( 6, 2, 4) = -0.70200d-02
        hbb( 6, 3, 4) = -0.72700d-02
        hbb( 6, 4, 4) = -0.39060d-01
        hbb( 6, 5, 4) = -0.38150d-01
        hbb( 6, 6, 4) = 0.77210d-01
        hbb( 7, 1, 4) = -0.89600d-02
        hbb( 7, 2, 4) = -0.20000d-04
        hbb( 7, 3, 4) = 0.89900d-02
        hbb( 7, 4, 4) = -0.32200d-01
        hbb( 7, 5, 4) = 0.13030d-01
        hbb( 7, 6, 4) = 0.19170d-01
        hbb( 7, 7, 4) = 0.79490d-01
        hbb( 8, 1, 4) = -0.12600d-02
        hbb( 8, 2, 4) = 0.34000d-02
        hbb( 8, 3, 4) = -0.21400d-02
        hbb( 8, 4, 4) = 0.18730d-01
        hbb( 8, 5, 4) = -0.69300d-02
        hbb( 8, 6, 4) = -0.11800d-01
        hbb( 8, 7, 4) = -0.39640d-01
        hbb( 8, 8, 4) = 0.79360d-01
        hbb( 9, 1, 4) = 0.10220d-01
        hbb( 9, 2, 4) = -0.33800d-02
        hbb( 9, 3, 4) = -0.68500d-02
        hbb( 9, 4, 4) = 0.13460d-01
        hbb( 9, 5, 4) = -0.60900d-02
        hbb( 9, 6, 4) = -0.73700d-02
        hbb( 9, 7, 4) = -0.39850d-01
        hbb( 9, 8, 4) = -0.39720d-01
        hbb( 9, 9, 4) = 0.79570d-01
        hbb(10, 1, 4) = 0.53800d-02
        hbb(10, 2, 4) = -0.26900d-02
        hbb(10, 3, 4) = -0.26900d-02
        hbb(10, 4, 4) = -0.85500d-02
        hbb(10, 5, 4) = 0.98800d-02
        hbb(10, 6, 4) = -0.13400d-02
        hbb(10, 7, 4) = -0.34400d-01
        hbb(10, 8, 4) = 0.14090d-01
        hbb(10, 9, 4) = 0.20310d-01
        hbb(10,10, 4) = 0.80510d-01
        hbb(11, 1, 4) = -0.26900d-02
        hbb(11, 2, 4) = 0.16300d-02
        hbb(11, 3, 4) = 0.10600d-02
        hbb(11, 4, 4) = -0.14300d-02
        hbb(11, 5, 4) = -0.33500d-02
        hbb(11, 6, 4) = 0.47900d-02
        hbb(11, 7, 4) = 0.20180d-01
        hbb(11, 8, 4) = -0.75700d-02
        hbb(11, 9, 4) = -0.12610d-01
        hbb(11,10, 4) = -0.40260d-01
        hbb(11,11, 4) = 0.79330d-01
        hbb(12, 1, 4) = -0.26900d-02
        hbb(12, 2, 4) = 0.10600d-02
        hbb(12, 3, 4) = 0.16300d-02
        hbb(12, 4, 4) = 0.99800d-02
        hbb(12, 5, 4) = -0.65300d-02
        hbb(12, 6, 4) = -0.34500d-02
        hbb(12, 7, 4) = 0.14220d-01
        hbb(12, 8, 4) = -0.65200d-02
        hbb(12, 9, 4) = -0.77100d-02
        hbb(12,10, 4) = -0.40260d-01
        hbb(12,11, 4) = -0.39080d-01
        hbb(12,12, 4) = 0.79330d-01
        hbb(13, 1, 4) = -0.89600d-02
        hbb(13, 2, 4) = 0.89900d-02
        hbb(13, 3, 4) = -0.20000d-04
        hbb(13, 4, 4) = 0.49400d-02
        hbb(13, 5, 4) = -0.20900d-02
        hbb(13, 6, 4) = -0.28500d-02
        hbb(13, 7, 4) = -0.88700d-02
        hbb(13, 8, 4) = 0.10300d-01
        hbb(13, 9, 4) = -0.14300d-02
        hbb(13,10, 4) = -0.34400d-01
        hbb(13,11, 4) = 0.14220d-01
        hbb(13,12, 4) = 0.20180d-01
        hbb(13,13, 4) = 0.79490d-01
        hbb(14, 1, 4) = 0.10220d-01
        hbb(14, 2, 4) = -0.68500d-02
        hbb(14, 3, 4) = -0.33800d-02
        hbb(14, 4, 4) = -0.27200d-02
        hbb(14, 5, 4) = 0.10700d-02
        hbb(14, 6, 4) = 0.16500d-02
        hbb(14, 7, 4) = -0.14300d-02
        hbb(14, 8, 4) = -0.36200d-02
        hbb(14, 9, 4) = 0.50500d-02
        hbb(14,10, 4) = 0.20310d-01
        hbb(14,11, 4) = -0.77100d-02
        hbb(14,12, 4) = -0.12610d-01
        hbb(14,13, 4) = -0.39850d-01
        hbb(14,14, 4) = 0.79570d-01
        hbb(15, 1, 4) = -0.12600d-02
        hbb(15, 2, 4) = -0.21400d-02
        hbb(15, 3, 4) = 0.34000d-02
        hbb(15, 4, 4) = -0.22200d-02
        hbb(15, 5, 4) = 0.10200d-02
        hbb(15, 6, 4) = 0.12000d-02
        hbb(15, 7, 4) = 0.10300d-01
        hbb(15, 8, 4) = -0.66700d-02
        hbb(15, 9, 4) = -0.36200d-02
        hbb(15,10, 4) = 0.14090d-01
        hbb(15,11, 4) = -0.65200d-02
        hbb(15,12, 4) = -0.75700d-02
        hbb(15,13, 4) = -0.39640d-01
        hbb(15,14, 4) = -0.39720d-01
        hbb(15,15, 4) = 0.79360d-01
        hbb(16, 1, 4) = -0.33740d-01
        hbb(16, 2, 4) = 0.15880d-01
        hbb(16, 3, 4) = 0.17860d-01
        hbb(16, 4, 4) = -0.83600d-02
        hbb(16, 5, 4) = -0.14300d-02
        hbb(16, 6, 4) = 0.97900d-02
        hbb(16, 7, 4) = 0.49400d-02
        hbb(16, 8, 4) = -0.22200d-02
        hbb(16, 9, 4) = -0.27200d-02
        hbb(16,10, 4) = -0.85500d-02
        hbb(16,11, 4) = 0.99800d-02
        hbb(16,12, 4) = -0.14300d-02
        hbb(16,13, 4) = -0.32200d-01
        hbb(16,14, 4) = 0.13460d-01
        hbb(16,15, 4) = 0.18730d-01
        hbb(16,16, 4) = 0.77910d-01
        hbb(17, 1, 4) = 0.19460d-01
        hbb(17, 2, 4) = -0.86100d-02
        hbb(17, 3, 4) = -0.10850d-01
        hbb(17, 4, 4) = -0.14300d-02
        hbb(17, 5, 4) = 0.48200d-02
        hbb(17, 6, 4) = -0.33900d-02
        hbb(17, 7, 4) = -0.20900d-02
        hbb(17, 8, 4) = 0.10200d-02
        hbb(17, 9, 4) = 0.10700d-02
        hbb(17,10, 4) = 0.98800d-02
        hbb(17,11, 4) = -0.65300d-02
        hbb(17,12, 4) = -0.33500d-02
        hbb(17,13, 4) = 0.13030d-01
        hbb(17,14, 4) = -0.60900d-02
        hbb(17,15, 4) = -0.69300d-02
        hbb(17,16, 4) = -0.38850d-01
        hbb(17,17, 4) = 0.77000d-01
        hbb(18, 1, 4) = 0.14280d-01
        hbb(18, 2, 4) = -0.72700d-02
        hbb(18, 3, 4) = -0.70200d-02
        hbb(18, 4, 4) = 0.97900d-02
        hbb(18, 5, 4) = -0.33900d-02
        hbb(18, 6, 4) = -0.63900d-02
        hbb(18, 7, 4) = -0.28500d-02
        hbb(18, 8, 4) = 0.12000d-02
        hbb(18, 9, 4) = 0.16500d-02
        hbb(18,10, 4) = -0.13400d-02
        hbb(18,11, 4) = -0.34500d-02
        hbb(18,12, 4) = 0.47900d-02
        hbb(18,13, 4) = 0.19170d-01
        hbb(18,14, 4) = -0.73700d-02
        hbb(18,15, 4) = -0.11800d-01
        hbb(18,16, 4) = -0.39060d-01
        hbb(18,17, 4) = -0.38150d-01
        hbb(18,18, 4) = 0.77210d-01
        ! Hessian elements for anchor point 4
        htt( 1, 1, 4) = 0.84400d-02
        htt( 2, 1, 4) = 0.31800d-02
        htt( 2, 2, 4) = 0.76400d-02
        htt( 3, 1, 4) = 0.17500d-02
        htt( 3, 2, 4) = -0.29500d-02
        htt( 3, 3, 4) = 0.10070d-01
        htt( 4, 1, 4) = -0.35100d-02
        htt( 4, 2, 4) = 0.15200d-02
        htt( 4, 3, 4) = 0.53700d-02
        htt( 4, 4, 4) = 0.10400d-01
        htt( 5, 1, 4) = -0.54900d-02
        htt( 5, 2, 4) = -0.37900d-02
        htt( 5, 3, 4) = 0.12100d-02
        htt( 5, 4, 4) = 0.29000d-02
        htt( 5, 5, 4) = 0.84400d-02
        htt( 6, 1, 4) = -0.37900d-02
        htt( 6, 2, 4) = -0.24900d-02
        htt( 6, 3, 4) = 0.23300d-02
        htt( 6, 4, 4) = 0.36400d-02
        htt( 6, 5, 4) = 0.31800d-02
        htt( 6, 6, 4) = 0.76400d-02
        htt( 7, 1, 4) = 0.12100d-02
        htt( 7, 2, 4) = 0.23300d-02
        htt( 7, 3, 4) = -0.71100d-02
        htt( 7, 4, 4) = -0.59900d-02
        htt( 7, 5, 4) = 0.17500d-02
        htt( 7, 6, 4) = -0.29500d-02
        htt( 7, 7, 4) = 0.10070d-01
        htt( 8, 1, 4) = 0.29000d-02
        htt( 8, 2, 4) = 0.36400d-02
        htt( 8, 3, 4) = -0.59900d-02
        htt( 8, 4, 4) = -0.52500d-02
        htt( 8, 5, 4) = -0.35100d-02
        htt( 8, 6, 4) = 0.15200d-02
        htt( 8, 7, 4) = 0.53700d-02
        htt( 8, 8, 4) = 0.10400d-01
        htt( 9, 1, 4) = -0.60400d-02
        htt( 9, 2, 4) = -0.50000d-03
        htt( 9, 3, 4) = -0.37000d-02
        htt( 9, 4, 4) = 0.18400d-02
        htt( 9, 5, 4) = 0.10000d-04
        htt( 9, 6, 4) = 0.17500d-02
        htt( 9, 7, 4) = -0.23300d-02
        htt( 9, 8, 4) = -0.59000d-03
        htt( 9, 9, 4) = 0.94900d-02
        htt(10, 1, 4) = -0.38600d-02
        htt(10, 2, 4) = 0.97000d-03
        htt(10, 3, 4) = -0.21900d-02
        htt(10, 4, 4) = 0.26400d-02
        htt(10, 5, 4) = 0.17900d-02
        htt(10, 6, 4) = 0.13600d-02
        htt(10, 7, 4) = 0.12000d-03
        htt(10, 8, 4) = -0.31000d-03
        htt(10, 9, 4) = 0.31000d-02
        htt(10,10, 4) = 0.84100d-02
        htt(11, 1, 4) = -0.65000d-03
        htt(11, 2, 4) = -0.50600d-02
        htt(11, 3, 4) = 0.11100d-02
        htt(11, 4, 4) = -0.33000d-02
        htt(11, 5, 4) = -0.17200d-02
        htt(11, 6, 4) = 0.41000d-03
        htt(11, 7, 4) = -0.34800d-02
        htt(11, 8, 4) = -0.13500d-02
        htt(11, 9, 4) = 0.38300d-02
        htt(11,10, 4) = -0.18300d-02
        htt(11,11, 4) = 0.83400d-02
        htt(12, 1, 4) = 0.15300d-02
        htt(12, 2, 4) = -0.35900d-02
        htt(12, 3, 4) = 0.26200d-02
        htt(12, 4, 4) = -0.25000d-02
        htt(12, 5, 4) = 0.60000d-04
        htt(12, 6, 4) = 0.20000d-04
        htt(12, 7, 4) = -0.10300d-02
        htt(12, 8, 4) = -0.10700d-02
        htt(12, 9, 4) = -0.25600d-02
        htt(12,10, 4) = 0.34800d-02
        htt(12,11, 4) = 0.26800d-02
        htt(12,12, 4) = 0.87100d-02
        htt(13, 1, 4) = 0.74000d-03
        htt(13, 2, 4) = -0.14900d-02
        htt(13, 3, 4) = 0.26300d-02
        htt(13, 4, 4) = 0.40000d-03
        htt(13, 5, 4) = 0.23400d-02
        htt(13, 6, 4) = 0.85000d-03
        htt(13, 7, 4) = 0.45000d-03
        htt(13, 8, 4) = -0.10400d-02
        htt(13, 9, 4) = -0.68300d-02
        htt(13,10, 4) = -0.31000d-03
        htt(13,11, 4) = -0.45500d-02
        htt(13,12, 4) = 0.19800d-02
        htt(13,13, 4) = 0.96800d-02
        htt(14, 1, 4) = 0.19800d-02
        htt(14, 2, 4) = 0.49000d-03
        htt(14, 3, 4) = 0.15200d-02
        htt(14, 4, 4) = 0.30000d-04
        htt(14, 5, 4) = 0.11000d-02
        htt(14, 6, 4) = -0.11300d-02
        htt(14, 7, 4) = 0.15600d-02
        htt(14, 8, 4) = -0.67000d-03
        htt(14, 9, 4) = -0.45100d-02
        htt(14,10, 4) = 0.12300d-02
        htt(14,11, 4) = -0.29800d-02
        htt(14,12, 4) = 0.27600d-02
        htt(14,13, 4) = 0.38800d-02
        htt(14,14, 4) = 0.86200d-02
        htt(15, 1, 4) = -0.14700d-02
        htt(15, 2, 4) = -0.29800d-02
        htt(15, 3, 4) = 0.11000d-02
        htt(15, 4, 4) = -0.41000d-03
        htt(15, 5, 4) = 0.54000d-03
        htt(15, 6, 4) = 0.12400d-02
        htt(15, 7, 4) = -0.20300d-02
        htt(15, 8, 4) = -0.13200d-02
        htt(15, 9, 4) = -0.36000d-03
        htt(15,10, 4) = -0.56900d-02
        htt(15,11, 4) = 0.11900d-02
        htt(15,12, 4) = -0.41400d-02
        htt(15,13, 4) = 0.30600d-02
        htt(15,14, 4) = -0.19400d-02
        htt(15,15, 4) = 0.84700d-02
        htt(16, 1, 4) = -0.22000d-03
        htt(16, 2, 4) = -0.10000d-02
        htt(16, 3, 4) = -0.10000d-04
        htt(16, 4, 4) = -0.78000d-03
        htt(16, 5, 4) = -0.71000d-03
        htt(16, 6, 4) = -0.74000d-03
        htt(16, 7, 4) = -0.92000d-03
        htt(16, 8, 4) = -0.95000d-03
        htt(16, 9, 4) = 0.19700d-02
        htt(16,10, 4) = -0.41500d-02
        htt(16,11, 4) = 0.27600d-02
        htt(16,12, 4) = -0.33600d-02
        htt(16,13, 4) = -0.27400d-02
        htt(16,14, 4) = 0.28000d-02
        htt(16,15, 4) = 0.34700d-02
        htt(16,16, 4) = 0.90000d-02
        htt(17, 1, 4) = 0.23400d-02
        htt(17, 2, 4) = 0.85000d-03
        htt(17, 3, 4) = 0.45000d-03
        htt(17, 4, 4) = -0.10400d-02
        htt(17, 5, 4) = 0.74000d-03
        htt(17, 6, 4) = -0.14900d-02
        htt(17, 7, 4) = 0.26300d-02
        htt(17, 8, 4) = 0.40000d-03
        htt(17, 9, 4) = 0.56000d-03
        htt(17,10, 4) = -0.18400d-02
        htt(17,11, 4) = 0.20800d-02
        htt(17,12, 4) = -0.32000d-03
        htt(17,13, 4) = -0.64700d-02
        htt(17,14, 4) = -0.67000d-03
        htt(17,15, 4) = -0.40300d-02
        htt(17,16, 4) = 0.17700d-02
        htt(17,17, 4) = 0.96800d-02
        htt(18, 1, 4) = 0.54000d-03
        htt(18, 2, 4) = 0.12400d-02
        htt(18, 3, 4) = -0.20300d-02
        htt(18, 4, 4) = -0.13200d-02
        htt(18, 5, 4) = -0.14700d-02
        htt(18, 6, 4) = -0.29800d-02
        htt(18, 7, 4) = 0.11000d-02
        htt(18, 8, 4) = -0.41000d-03
        htt(18, 9, 4) = 0.22500d-02
        htt(18,10, 4) = 0.51000d-03
        htt(18,11, 4) = 0.15300d-02
        htt(18,12, 4) = -0.22000d-03
        htt(18,13, 4) = -0.40300d-02
        htt(18,14, 4) = 0.97000d-03
        htt(18,15, 4) = -0.22600d-02
        htt(18,16, 4) = 0.27400d-02
        htt(18,17, 4) = 0.30600d-02
        htt(18,18, 4) = 0.84700d-02
        htt(19, 1, 4) = 0.11000d-02
        htt(19, 2, 4) = -0.11300d-02
        htt(19, 3, 4) = 0.15600d-02
        htt(19, 4, 4) = -0.67000d-03
        htt(19, 5, 4) = 0.19800d-02
        htt(19, 6, 4) = 0.49000d-03
        htt(19, 7, 4) = 0.15200d-02
        htt(19, 8, 4) = 0.30000d-04
        htt(19, 9, 4) = -0.17700d-02
        htt(19,10, 4) = -0.33800d-02
        htt(19,11, 4) = 0.51000d-03
        htt(19,12, 4) = -0.11100d-02
        htt(19,13, 4) = -0.67000d-03
        htt(19,14, 4) = -0.54000d-02
        htt(19,15, 4) = 0.97000d-03
        htt(19,16, 4) = -0.37700d-02
        htt(19,17, 4) = 0.38800d-02
        htt(19,18, 4) = -0.19400d-02
        htt(19,19, 4) = 0.86200d-02
        htt(20, 1, 4) = -0.71000d-03
        htt(20, 2, 4) = -0.74000d-03
        htt(20, 3, 4) = -0.92000d-03
        htt(20, 4, 4) = -0.95000d-03
        htt(20, 5, 4) = -0.22000d-03
        htt(20, 6, 4) = -0.10000d-02
        htt(20, 7, 4) = -0.10000d-04
        htt(20, 8, 4) = -0.78000d-03
        htt(20, 9, 4) = -0.70000d-04
        htt(20,10, 4) = -0.10300d-02
        htt(20,11, 4) = -0.40000d-04
        htt(20,12, 4) = -0.10000d-02
        htt(20,13, 4) = 0.17700d-02
        htt(20,14, 4) = -0.37700d-02
        htt(20,15, 4) = 0.27400d-02
        htt(20,16, 4) = -0.28000d-02
        htt(20,17, 4) = -0.27400d-02
        htt(20,18, 4) = 0.34700d-02
        htt(20,19, 4) = 0.28000d-02
        htt(20,20, 4) = 0.90000d-02
        htt(21, 1, 4) = 0.10000d-04
        htt(21, 2, 4) = 0.17500d-02
        htt(21, 3, 4) = -0.23300d-02
        htt(21, 4, 4) = -0.59000d-03
        htt(21, 5, 4) = -0.60400d-02
        htt(21, 6, 4) = -0.50000d-03
        htt(21, 7, 4) = -0.37000d-02
        htt(21, 8, 4) = 0.18400d-02
        htt(21, 9, 4) = 0.27700d-02
        htt(21,10, 4) = 0.11000d-02
        htt(21,11, 4) = 0.10000d-02
        htt(21,12, 4) = -0.68000d-03
        htt(21,13, 4) = 0.56000d-03
        htt(21,14, 4) = -0.17700d-02
        htt(21,15, 4) = 0.22500d-02
        htt(21,16, 4) = -0.70000d-04
        htt(21,17, 4) = -0.68300d-02
        htt(21,18, 4) = -0.36000d-03
        htt(21,19, 4) = -0.45100d-02
        htt(21,20, 4) = 0.19700d-02
        htt(21,21, 4) = 0.94900d-02
        htt(22, 1, 4) = -0.17200d-02
        htt(22, 2, 4) = 0.41000d-03
        htt(22, 3, 4) = -0.34800d-02
        htt(22, 4, 4) = -0.13500d-02
        htt(22, 5, 4) = -0.65000d-03
        htt(22, 6, 4) = -0.50600d-02
        htt(22, 7, 4) = 0.11100d-02
        htt(22, 8, 4) = -0.33000d-02
        htt(22, 9, 4) = 0.10000d-02
        htt(22,10, 4) = 0.15400d-02
        htt(22,11, 4) = -0.11800d-02
        htt(22,12, 4) = -0.64000d-03
        htt(22,13, 4) = 0.20800d-02
        htt(22,14, 4) = 0.51000d-03
        htt(22,15, 4) = 0.15300d-02
        htt(22,16, 4) = -0.40000d-04
        htt(22,17, 4) = -0.45500d-02
        htt(22,18, 4) = 0.11900d-02
        htt(22,19, 4) = -0.29800d-02
        htt(22,20, 4) = 0.27600d-02
        htt(22,21, 4) = 0.38300d-02
        htt(22,22, 4) = 0.83400d-02
        htt(23, 1, 4) = 0.17900d-02
        htt(23, 2, 4) = 0.13600d-02
        htt(23, 3, 4) = 0.12000d-03
        htt(23, 4, 4) = -0.31000d-03
        htt(23, 5, 4) = -0.38600d-02
        htt(23, 6, 4) = 0.97000d-03
        htt(23, 7, 4) = -0.21900d-02
        htt(23, 8, 4) = 0.26400d-02
        htt(23, 9, 4) = 0.11000d-02
        htt(23,10, 4) = -0.12200d-02
        htt(23,11, 4) = 0.15400d-02
        htt(23,12, 4) = -0.78000d-03
        htt(23,13, 4) = -0.18400d-02
        htt(23,14, 4) = -0.33800d-02
        htt(23,15, 4) = 0.51000d-03
        htt(23,16, 4) = -0.10300d-02
        htt(23,17, 4) = -0.31000d-03
        htt(23,18, 4) = -0.56900d-02
        htt(23,19, 4) = 0.12300d-02
        htt(23,20, 4) = -0.41500d-02
        htt(23,21, 4) = 0.31000d-02
        htt(23,22, 4) = -0.18300d-02
        htt(23,23, 4) = 0.84100d-02
        htt(24, 1, 4) = 0.60000d-04
        htt(24, 2, 4) = 0.20000d-04
        htt(24, 3, 4) = -0.10300d-02
        htt(24, 4, 4) = -0.10700d-02
        htt(24, 5, 4) = 0.15300d-02
        htt(24, 6, 4) = -0.35900d-02
        htt(24, 7, 4) = 0.26200d-02
        htt(24, 8, 4) = -0.25000d-02
        htt(24, 9, 4) = -0.68000d-03
        htt(24,10, 4) = -0.78000d-03
        htt(24,11, 4) = -0.64000d-03
        htt(24,12, 4) = -0.74000d-03
        htt(24,13, 4) = -0.32000d-03
        htt(24,14, 4) = -0.11100d-02
        htt(24,15, 4) = -0.22000d-03
        htt(24,16, 4) = -0.10000d-02
        htt(24,17, 4) = 0.19800d-02
        htt(24,18, 4) = -0.41400d-02
        htt(24,19, 4) = 0.27600d-02
        htt(24,20, 4) = -0.33600d-02
        htt(24,21, 4) = -0.25600d-02
        htt(24,22, 4) = 0.26800d-02
        htt(24,23, 4) = 0.34800d-02
        htt(24,24, 4) = 0.87100d-02
        ! Hessian elements for anchor point 4
        hbs( 1, 1, 4) = -0.31760d-01
        hbs( 1, 2, 4) = -0.31760d-01
        hbs( 1, 3, 4) = -0.37700d-01
        hbs( 1, 4, 4) = -0.11100d-01
        hbs( 1, 5, 4) = 0.40600d-02
        hbs( 1, 6, 4) = 0.43670d-01
        hbs( 1, 7, 4) = -0.35000d-03
        hbs( 1, 8, 4) = 0.43670d-01
        hbs( 1, 9, 4) = -0.10000d-04
        hbs( 1,10, 4) = -0.11100d-01
        hbs( 1,11, 4) = -0.35000d-03
        hbs( 1,12, 4) = 0.40600d-02
        hbs( 2, 1, 4) = 0.39800d-01
        hbs( 2, 2, 4) = -0.80400d-02
        hbs( 2, 3, 4) = 0.18850d-01
        hbs( 2, 4, 4) = -0.12000d-03
        hbs( 2, 5, 4) = -0.21800d-02
        hbs( 2, 6, 4) = -0.18210d-01
        hbs( 2, 7, 4) = 0.75000d-03
        hbs( 2, 8, 4) = -0.25460d-01
        hbs( 2, 9, 4) = 0.00000d+00
        hbs( 2,10, 4) = 0.11220d-01
        hbs( 2,11, 4) = -0.39000d-03
        hbs( 2,12, 4) = -0.18800d-02
        hbs( 3, 1, 4) = -0.80400d-02
        hbs( 3, 2, 4) = 0.39800d-01
        hbs( 3, 3, 4) = 0.18850d-01
        hbs( 3, 4, 4) = 0.11220d-01
        hbs( 3, 5, 4) = -0.18800d-02
        hbs( 3, 6, 4) = -0.25460d-01
        hbs( 3, 7, 4) = -0.39000d-03
        hbs( 3, 8, 4) = -0.18210d-01
        hbs( 3, 9, 4) = 0.00000d+00
        hbs( 3,10, 4) = -0.12000d-03
        hbs( 3,11, 4) = 0.75000d-03
        hbs( 3,12, 4) = -0.21800d-02
        hbs( 4, 1, 4) = -0.26690d-01
        hbs( 4, 2, 4) = -0.99100d-02
        hbs( 4, 3, 4) = 0.10200d-01
        hbs( 4, 4, 4) = -0.32600d-01
        hbs( 4, 5, 4) = -0.82100d-02
        hbs( 4, 6, 4) = -0.12130d-01
        hbs( 4, 7, 4) = 0.47800d-02
        hbs( 4, 8, 4) = 0.46360d-01
        hbs( 4, 9, 4) = -0.61000d-03
        hbs( 4,10, 4) = 0.44260d-01
        hbs( 4,11, 4) = 0.10000d-03
        hbs( 4,12, 4) = -0.43000d-03
        hbs( 5, 1, 4) = 0.27720d-01
        hbs( 5, 2, 4) = 0.58200d-02
        hbs( 5, 3, 4) = -0.92100d-02
        hbs( 5, 4, 4) = 0.33000d-02
        hbs( 5, 5, 4) = 0.59500d-02
        hbs( 5, 6, 4) = 0.56400d-02
        hbs( 5, 7, 4) = -0.18300d-02
        hbs( 5, 8, 4) = -0.20790d-01
        hbs( 5, 9, 4) = -0.31000d-03
        hbs( 5,10, 4) = -0.24570d-01
        hbs( 5,11, 4) = -0.60000d-04
        hbs( 5,12, 4) = 0.64000d-03
        hbs( 6, 1, 4) = -0.10300d-02
        hbs( 6, 2, 4) = 0.40900d-02
        hbs( 6, 3, 4) = -0.99000d-03
        hbs( 6, 4, 4) = 0.29300d-01
        hbs( 6, 5, 4) = 0.22600d-02
        hbs( 6, 6, 4) = 0.65000d-02
        hbs( 6, 7, 4) = -0.29500d-02
        hbs( 6, 8, 4) = -0.25570d-01
        hbs( 6, 9, 4) = 0.92000d-03
        hbs( 6,10, 4) = -0.19690d-01
        hbs( 6,11, 4) = -0.40000d-04
        hbs( 6,12, 4) = -0.21000d-03
        hbs( 7, 1, 4) = -0.14540d-01
        hbs( 7, 2, 4) = 0.40000d-01
        hbs( 7, 3, 4) = 0.98300d-02
        hbs( 7, 4, 4) = -0.32290d-01
        hbs( 7, 5, 4) = 0.50700d-02
        hbs( 7, 6, 4) = -0.32910d-01
        hbs( 7, 7, 4) = -0.88600d-02
        hbs( 7, 8, 4) = -0.12410d-01
        hbs( 7, 9, 4) = 0.48800d-02
        hbs( 7,10, 4) = 0.45460d-01
        hbs( 7,11, 4) = -0.87000d-03
        hbs( 7,12, 4) = 0.14000d-03
        hbs( 8, 1, 4) = 0.89100d-02
        hbs( 8, 2, 4) = -0.22340d-01
        hbs( 8, 3, 4) = -0.58300d-02
        hbs( 8, 4, 4) = 0.30270d-01
        hbs( 8, 5, 4) = -0.31400d-02
        hbs( 8, 6, 4) = 0.36400d-02
        hbs( 8, 7, 4) = 0.47300d-02
        hbs( 8, 8, 4) = 0.54700d-02
        hbs( 8, 9, 4) = -0.19600d-02
        hbs( 8,10, 4) = -0.20340d-01
        hbs( 8,11, 4) = -0.23000d-03
        hbs( 8,12, 4) = -0.12000d-03
        hbs( 9, 1, 4) = 0.56200d-02
        hbs( 9, 2, 4) = -0.17660d-01
        hbs( 9, 3, 4) = -0.39900d-02
        hbs( 9, 4, 4) = 0.20200d-02
        hbs( 9, 5, 4) = -0.19200d-02
        hbs( 9, 6, 4) = 0.29280d-01
        hbs( 9, 7, 4) = 0.41300d-02
        hbs( 9, 8, 4) = 0.69400d-02
        hbs( 9, 9, 4) = -0.29200d-02
        hbs( 9,10, 4) = -0.25120d-01
        hbs( 9,11, 4) = 0.11000d-02
        hbs( 9,12, 4) = -0.30000d-04
        hbs(10, 1, 4) = 0.42900d-01
        hbs(10, 2, 4) = 0.42900d-01
        hbs(10, 3, 4) = -0.23400d-02
        hbs(10, 4, 4) = -0.13730d-01
        hbs(10, 5, 4) = -0.64000d-03
        hbs(10, 6, 4) = -0.32580d-01
        hbs(10, 7, 4) = 0.52100d-02
        hbs(10, 8, 4) = -0.32580d-01
        hbs(10, 9, 4) = -0.85400d-02
        hbs(10,10, 4) = -0.13730d-01
        hbs(10,11, 4) = 0.52100d-02
        hbs(10,12, 4) = -0.64000d-03
        hbs(11, 1, 4) = -0.23760d-01
        hbs(11, 2, 4) = -0.19130d-01
        hbs(11, 3, 4) = 0.11700d-02
        hbs(11, 4, 4) = 0.73400d-02
        hbs(11, 5, 4) = 0.91000d-03
        hbs(11, 6, 4) = 0.30250d-01
        hbs(11, 7, 4) = -0.31400d-02
        hbs(11, 8, 4) = 0.23200d-02
        hbs(11, 9, 4) = 0.42700d-02
        hbs(11,10, 4) = 0.64000d-02
        hbs(11,11, 4) = -0.20700d-02
        hbs(11,12, 4) = -0.27000d-03
        hbs(12, 1, 4) = -0.19130d-01
        hbs(12, 2, 4) = -0.23760d-01
        hbs(12, 3, 4) = 0.11700d-02
        hbs(12, 4, 4) = 0.64000d-02
        hbs(12, 5, 4) = -0.27000d-03
        hbs(12, 6, 4) = 0.23200d-02
        hbs(12, 7, 4) = -0.20700d-02
        hbs(12, 8, 4) = 0.30250d-01
        hbs(12, 9, 4) = 0.42700d-02
        hbs(12,10, 4) = 0.73400d-02
        hbs(12,11, 4) = -0.31400d-02
        hbs(12,12, 4) = 0.91000d-03
        hbs(13, 1, 4) = 0.40000d-01
        hbs(13, 2, 4) = -0.14540d-01
        hbs(13, 3, 4) = 0.98300d-02
        hbs(13, 4, 4) = 0.45460d-01
        hbs(13, 5, 4) = 0.14000d-03
        hbs(13, 6, 4) = -0.12410d-01
        hbs(13, 7, 4) = -0.87000d-03
        hbs(13, 8, 4) = -0.32910d-01
        hbs(13, 9, 4) = 0.48800d-02
        hbs(13,10, 4) = -0.32290d-01
        hbs(13,11, 4) = -0.88600d-02
        hbs(13,12, 4) = 0.50700d-02
        hbs(14, 1, 4) = -0.17660d-01
        hbs(14, 2, 4) = 0.56200d-02
        hbs(14, 3, 4) = -0.39900d-02
        hbs(14, 4, 4) = -0.25120d-01
        hbs(14, 5, 4) = -0.30000d-04
        hbs(14, 6, 4) = 0.69400d-02
        hbs(14, 7, 4) = 0.11000d-02
        hbs(14, 8, 4) = 0.29280d-01
        hbs(14, 9, 4) = -0.29200d-02
        hbs(14,10, 4) = 0.20200d-02
        hbs(14,11, 4) = 0.41300d-02
        hbs(14,12, 4) = -0.19200d-02
        hbs(15, 1, 4) = -0.22340d-01
        hbs(15, 2, 4) = 0.89100d-02
        hbs(15, 3, 4) = -0.58300d-02
        hbs(15, 4, 4) = -0.20340d-01
        hbs(15, 5, 4) = -0.12000d-03
        hbs(15, 6, 4) = 0.54700d-02
        hbs(15, 7, 4) = -0.23000d-03
        hbs(15, 8, 4) = 0.36400d-02
        hbs(15, 9, 4) = -0.19600d-02
        hbs(15,10, 4) = 0.30270d-01
        hbs(15,11, 4) = 0.47300d-02
        hbs(15,12, 4) = -0.31400d-02
        hbs(16, 1, 4) = -0.99100d-02
        hbs(16, 2, 4) = -0.26690d-01
        hbs(16, 3, 4) = 0.10200d-01
        hbs(16, 4, 4) = 0.44260d-01
        hbs(16, 5, 4) = -0.43000d-03
        hbs(16, 6, 4) = 0.46360d-01
        hbs(16, 7, 4) = 0.10000d-03
        hbs(16, 8, 4) = -0.12130d-01
        hbs(16, 9, 4) = -0.61000d-03
        hbs(16,10, 4) = -0.32600d-01
        hbs(16,11, 4) = 0.47800d-02
        hbs(16,12, 4) = -0.82100d-02
        hbs(17, 1, 4) = 0.58200d-02
        hbs(17, 2, 4) = 0.27720d-01
        hbs(17, 3, 4) = -0.92100d-02
        hbs(17, 4, 4) = -0.24570d-01
        hbs(17, 5, 4) = 0.64000d-03
        hbs(17, 6, 4) = -0.20790d-01
        hbs(17, 7, 4) = -0.60000d-04
        hbs(17, 8, 4) = 0.56400d-02
        hbs(17, 9, 4) = -0.31000d-03
        hbs(17,10, 4) = 0.33000d-02
        hbs(17,11, 4) = -0.18300d-02
        hbs(17,12, 4) = 0.59500d-02
        hbs(18, 1, 4) = 0.40900d-02
        hbs(18, 2, 4) = -0.10300d-02
        hbs(18, 3, 4) = -0.99000d-03
        hbs(18, 4, 4) = -0.19690d-01
        hbs(18, 5, 4) = -0.21000d-03
        hbs(18, 6, 4) = -0.25570d-01
        hbs(18, 7, 4) = -0.40000d-04
        hbs(18, 8, 4) = 0.65000d-02
        hbs(18, 9, 4) = 0.92000d-03
        hbs(18,10, 4) = 0.29300d-01
        hbs(18,11, 4) = -0.29500d-02
        hbs(18,12, 4) = 0.22600d-02
    else if(ii .eq. 22) then
        ! Hessian elements for anchor point 1
        hss( 1, 1, 1) = 0.35006d+00
        hss( 2, 1, 1) = 0.13620d-01
        hss( 2, 2, 1) = 0.36874d+00
        hss( 3, 1, 1) = 0.29130d-01
        hss( 3, 2, 1) = 0.28430d-01
        hss( 3, 3, 1) = 0.46705d+00
        hss( 4, 1, 1) = -0.29900d-02
        hss( 4, 2, 1) = 0.16080d-01
        hss( 4, 3, 1) = -0.77700d-02
        hss( 4, 4, 1) = 0.38651d+00
        hss( 5, 1, 1) = 0.40700d-02
        hss( 5, 2, 1) = -0.10000d-02
        hss( 5, 3, 1) = 0.38700d-02
        hss( 5, 4, 1) = 0.32400d-02
        hss( 5, 5, 1) = 0.36515d+00
        hss( 6, 1, 1) = 0.33920d-01
        hss( 6, 2, 1) = 0.46890d-01
        hss( 6, 3, 1) = -0.39100d-02
        hss( 6, 4, 1) = -0.52200d-02
        hss( 6, 5, 1) = 0.40000d-03
        hss( 6, 6, 1) = 0.36166d+00
        hss( 7, 1, 1) = -0.40000d-03
        hss( 7, 2, 1) = -0.70000d-03
        hss( 7, 3, 1) = 0.45000d-03
        hss( 7, 4, 1) = 0.39500d-02
        hss( 7, 5, 1) = 0.42000d-03
        hss( 7, 6, 1) = 0.39100d-02
        hss( 7, 7, 1) = 0.37071d+00
        hss( 8, 1, 1) = 0.44240d-01
        hss( 8, 2, 1) = 0.34840d-01
        hss( 8, 3, 1) = -0.58500d-02
        hss( 8, 4, 1) = 0.24280d-01
        hss( 8, 5, 1) = -0.11400d-02
        hss( 8, 6, 1) = 0.14590d-01
        hss( 8, 7, 1) = -0.98000d-03
        hss( 8, 8, 1) = 0.36092d+00
        hss( 9, 1, 1) = -0.57000d-03
        hss( 9, 2, 1) = -0.54000d-03
        hss( 9, 3, 1) = 0.41000d-03
        hss( 9, 4, 1) = -0.31000d-03
        hss( 9, 5, 1) = 0.17000d-03
        hss( 9, 6, 1) = 0.39200d-02
        hss( 9, 7, 1) = 0.50000d-03
        hss( 9, 8, 1) = 0.39400d-02
        hss( 9, 9, 1) = 0.36389d+00
        hss(10, 1, 1) = 0.19730d-01
        hss(10, 2, 1) = -0.83300d-02
        hss(10, 3, 1) = -0.46800d-02
        hss(10, 4, 1) = 0.63810d-01
        hss(10, 5, 1) = -0.33000d-03
        hss(10, 6, 1) = 0.22020d-01
        hss(10, 7, 1) = -0.10900d-02
        hss(10, 8, 1) = -0.21000d-02
        hss(10, 9, 1) = -0.19000d-03
        hss(10,10, 1) = 0.38332d+00
        hss(11, 1, 1) = -0.70000d-03
        hss(11, 2, 1) = -0.43000d-03
        hss(11, 3, 1) = 0.55000d-03
        hss(11, 4, 1) = -0.11400d-02
        hss(11, 5, 1) = 0.40000d-04
        hss(11, 6, 1) = -0.97000d-03
        hss(11, 7, 1) = 0.22000d-03
        hss(11, 8, 1) = 0.38800d-02
        hss(11, 9, 1) = 0.48000d-03
        hss(11,10, 1) = 0.39000d-02
        hss(11,11, 1) = 0.37178d+00
        hss(12, 1, 1) = -0.14000d-02
        hss(12, 2, 1) = 0.34100d-02
        hss(12, 3, 1) = 0.37000d-03
        hss(12, 4, 1) = -0.13000d-03
        hss(12, 5, 1) = 0.24000d-03
        hss(12, 6, 1) = -0.10500d-02
        hss(12, 7, 1) = 0.20000d-04
        hss(12, 8, 1) = 0.55000d-03
        hss(12, 9, 1) = 0.12000d-03
        hss(12,10, 1) = 0.24400d-02
        hss(12,11, 1) = 0.39000d-03
        hss(12,12, 1) = 0.37251d+00
        ! Hessian elements for anchor point 1
        hbb( 1, 1, 1) = 0.82220d-01
        hbb( 2, 1, 1) = -0.42040d-01
        hbb( 2, 2, 1) = 0.14128d+00
        hbb( 3, 1, 1) = -0.40180d-01
        hbb( 3, 2, 1) = -0.99240d-01
        hbb( 3, 3, 1) = 0.13942d+00
        hbb( 4, 1, 1) = -0.36210d-01
        hbb( 4, 2, 1) = 0.22250d-01
        hbb( 4, 3, 1) = 0.13960d-01
        hbb( 4, 4, 1) = 0.66190d-01
        hbb( 5, 1, 1) = 0.21220d-01
        hbb( 5, 2, 1) = -0.13320d-01
        hbb( 5, 3, 1) = -0.79000d-02
        hbb( 5, 4, 1) = -0.33260d-01
        hbb( 5, 5, 1) = 0.72740d-01
        hbb( 6, 1, 1) = 0.14990d-01
        hbb( 6, 2, 1) = -0.89300d-02
        hbb( 6, 3, 1) = -0.60600d-02
        hbb( 6, 4, 1) = -0.32930d-01
        hbb( 6, 5, 1) = -0.39480d-01
        hbb( 6, 6, 1) = 0.72410d-01
        hbb( 7, 1, 1) = -0.67600d-02
        hbb( 7, 2, 1) = -0.46000d-02
        hbb( 7, 3, 1) = 0.11370d-01
        hbb( 7, 4, 1) = -0.26080d-01
        hbb( 7, 5, 1) = 0.93900d-02
        hbb( 7, 6, 1) = 0.16690d-01
        hbb( 7, 7, 1) = 0.64820d-01
        hbb( 8, 1, 1) = -0.24700d-02
        hbb( 8, 2, 1) = 0.77600d-02
        hbb( 8, 3, 1) = -0.52900d-02
        hbb( 8, 4, 1) = 0.15900d-01
        hbb( 8, 5, 1) = -0.49900d-02
        hbb( 8, 6, 1) = -0.10910d-01
        hbb( 8, 7, 1) = -0.31790d-01
        hbb( 8, 8, 1) = 0.73440d-01
        hbb( 9, 1, 1) = 0.92300d-02
        hbb( 9, 2, 1) = -0.31500d-02
        hbb( 9, 3, 1) = -0.60800d-02
        hbb( 9, 4, 1) = 0.10190d-01
        hbb( 9, 5, 1) = -0.44000d-02
        hbb( 9, 6, 1) = -0.57900d-02
        hbb( 9, 7, 1) = -0.33030d-01
        hbb( 9, 8, 1) = -0.41650d-01
        hbb( 9, 9, 1) = 0.74680d-01
        hbb(10, 1, 1) = 0.38400d-02
        hbb(10, 2, 1) = -0.14400d-02
        hbb(10, 3, 1) = -0.24000d-02
        hbb(10, 4, 1) = -0.57300d-02
        hbb(10, 5, 1) = 0.88600d-02
        hbb(10, 6, 1) = -0.31300d-02
        hbb(10, 7, 1) = -0.33500d-01
        hbb(10, 8, 1) = 0.13490d-01
        hbb(10, 9, 1) = 0.20010d-01
        hbb(10,10, 1) = 0.75070d-01
        hbb(11, 1, 1) = -0.18400d-02
        hbb(11, 2, 1) = 0.11400d-02
        hbb(11, 3, 1) = 0.70000d-03
        hbb(11, 4, 1) = -0.27200d-02
        hbb(11, 5, 1) = -0.30000d-02
        hbb(11, 6, 1) = 0.57200d-02
        hbb(11, 7, 1) = 0.20110d-01
        hbb(11, 8, 1) = -0.74900d-02
        hbb(11, 9, 1) = -0.12620d-01
        hbb(11,10, 1) = -0.37620d-01
        hbb(11,11, 1) = 0.78370d-01
        hbb(12, 1, 1) = -0.20000d-02
        hbb(12, 2, 1) = 0.30000d-03
        hbb(12, 3, 1) = 0.17000d-02
        hbb(12, 4, 1) = 0.84500d-02
        hbb(12, 5, 1) = -0.58700d-02
        hbb(12, 6, 1) = -0.25800d-02
        hbb(12, 7, 1) = 0.13390d-01
        hbb(12, 8, 1) = -0.60000d-02
        hbb(12, 9, 1) = -0.73900d-02
        hbb(12,10, 1) = -0.37460d-01
        hbb(12,11, 1) = -0.40760d-01
        hbb(12,12, 1) = 0.78220d-01
        hbb(13, 1, 1) = -0.61900d-02
        hbb(13, 2, 1) = 0.10630d-01
        hbb(13, 3, 1) = -0.44400d-02
        hbb(13, 4, 1) = -0.14000d-02
        hbb(13, 5, 1) = 0.11900d-02
        hbb(13, 6, 1) = 0.21000d-03
        hbb(13, 7, 1) = 0.29200d-02
        hbb(13, 8, 1) = 0.41800d-02
        hbb(13, 9, 1) = -0.71000d-02
        hbb(13,10, 1) = -0.33570d-01
        hbb(13,11, 1) = 0.13390d-01
        hbb(13,12, 1) = 0.20180d-01
        hbb(13,13, 1) = 0.65200d-01
        hbb(14, 1, 1) = 0.89900d-02
        hbb(14, 2, 1) = -0.58800d-02
        hbb(14, 3, 1) = -0.31100d-02
        hbb(14, 4, 1) = 0.77000d-03
        hbb(14, 5, 1) = -0.83000d-03
        hbb(14, 6, 1) = 0.70000d-04
        hbb(14, 7, 1) = -0.71700d-02
        hbb(14, 8, 1) = -0.51000d-03
        hbb(14, 9, 1) = 0.76800d-02
        hbb(14,10, 1) = 0.20000d-01
        hbb(14,11, 1) = -0.73300d-02
        hbb(14,12, 1) = -0.12670d-01
        hbb(14,13, 1) = -0.33080d-01
        hbb(14,14, 1) = 0.74340d-01
        hbb(15, 1, 1) = -0.28000d-02
        hbb(15, 2, 1) = -0.47500d-02
        hbb(15, 3, 1) = 0.75500d-02
        hbb(15, 4, 1) = 0.64000d-03
        hbb(15, 5, 1) = -0.36000d-03
        hbb(15, 6, 1) = -0.28000d-03
        hbb(15, 7, 1) = 0.42500d-02
        hbb(15, 8, 1) = -0.36700d-02
        hbb(15, 9, 1) = -0.58000d-03
        hbb(15,10, 1) = 0.13570d-01
        hbb(15,11, 1) = -0.60500d-02
        hbb(15,12, 1) = -0.75200d-02
        hbb(15,13, 1) = -0.32120d-01
        hbb(15,14, 1) = -0.41260d-01
        hbb(15,15, 1) = 0.73390d-01
        hbb(16, 1, 1) = -0.36890d-01
        hbb(16, 2, 1) = 0.15200d-01
        hbb(16, 3, 1) = 0.21690d-01
        hbb(16, 4, 1) = 0.32300d-02
        hbb(16, 5, 1) = -0.74000d-02
        hbb(16, 6, 1) = 0.41700d-02
        hbb(16, 7, 1) = -0.13900d-02
        hbb(16, 8, 1) = 0.69000d-03
        hbb(16, 9, 1) = 0.71000d-03
        hbb(16,10, 1) = -0.61100d-02
        hbb(16,11, 1) = 0.86800d-02
        hbb(16,12, 1) = -0.25700d-02
        hbb(16,13, 1) = -0.26960d-01
        hbb(16,14, 1) = 0.10490d-01
        hbb(16,15, 1) = 0.16470d-01
        hbb(16,16, 1) = 0.68110d-01
        hbb(17, 1, 1) = 0.22720d-01
        hbb(17, 2, 1) = -0.69100d-02
        hbb(17, 3, 1) = -0.15810d-01
        hbb(17, 4, 1) = -0.78500d-02
        hbb(17, 5, 1) = 0.81700d-02
        hbb(17, 6, 1) = -0.32000d-03
        hbb(17, 7, 1) = 0.10300d-02
        hbb(17, 8, 1) = -0.32000d-03
        hbb(17, 9, 1) = -0.71000d-03
        hbb(17,10, 1) = 0.90700d-02
        hbb(17,11, 1) = -0.60300d-02
        hbb(17,12, 1) = -0.30400d-02
        hbb(17,13, 1) = 0.97900d-02
        hbb(17,14, 1) = -0.45700d-02
        hbb(17,15, 1) = -0.52200d-02
        hbb(17,16, 1) = -0.34760d-01
        hbb(17,17, 1) = 0.72170d-01
        hbb(18, 1, 1) = 0.14170d-01
        hbb(18, 2, 1) = -0.82800d-02
        hbb(18, 3, 1) = -0.58800d-02
        hbb(18, 4, 1) = 0.46200d-02
        hbb(18, 5, 1) = -0.77000d-03
        hbb(18, 6, 1) = -0.38500d-02
        hbb(18, 7, 1) = 0.36000d-03
        hbb(18, 8, 1) = -0.37000d-03
        hbb(18, 9, 1) = 0.10000d-04
        hbb(18,10, 1) = -0.29600d-02
        hbb(18,11, 1) = -0.26500d-02
        hbb(18,12, 1) = 0.56000d-02
        hbb(18,13, 1) = 0.17170d-01
        hbb(18,14, 1) = -0.59200d-02
        hbb(18,15, 1) = -0.11250d-01
        hbb(18,16, 1) = -0.33360d-01
        hbb(18,17, 1) = -0.37420d-01
        hbb(18,18, 1) = 0.70770d-01
        ! Hessian elements for anchor point 1
        htt( 1, 1, 1) = 0.59700d-02
        htt( 2, 1, 1) = 0.40300d-02
        htt( 2, 2, 1) = 0.53700d-02
        htt( 3, 1, 1) = -0.14500d-02
        htt( 3, 2, 1) = -0.35400d-02
        htt( 3, 3, 1) = 0.68100d-02
        htt( 4, 1, 1) = -0.33900d-02
        htt( 4, 2, 1) = -0.21900d-02
        htt( 4, 3, 1) = 0.47200d-02
        htt( 4, 4, 1) = 0.59200d-02
        htt( 5, 1, 1) = -0.57900d-02
        htt( 5, 2, 1) = -0.45300d-02
        htt( 5, 3, 1) = 0.17600d-02
        htt( 5, 4, 1) = 0.30200d-02
        htt( 5, 5, 1) = 0.64900d-02
        htt( 6, 1, 1) = -0.38500d-02
        htt( 6, 2, 1) = -0.41000d-02
        htt( 6, 3, 1) = 0.21800d-02
        htt( 6, 4, 1) = 0.19300d-02
        htt( 6, 5, 1) = 0.34500d-02
        htt( 6, 6, 1) = 0.52100d-02
        htt( 7, 1, 1) = 0.13600d-02
        htt( 7, 2, 1) = 0.27600d-02
        htt( 7, 3, 1) = -0.62000d-02
        htt( 7, 4, 1) = -0.48000d-02
        htt( 7, 5, 1) = -0.78000d-03
        htt( 7, 6, 1) = -0.23700d-02
        htt( 7, 7, 1) = 0.65100d-02
        htt( 8, 1, 1) = 0.32900d-02
        htt( 8, 2, 1) = 0.31900d-02
        htt( 8, 3, 1) = -0.57800d-02
        htt( 8, 4, 1) = -0.58800d-02
        htt( 8, 5, 1) = -0.38200d-02
        htt( 8, 6, 1) = -0.61000d-03
        htt( 8, 7, 1) = 0.49200d-02
        htt( 8, 8, 1) = 0.81300d-02
        htt( 9, 1, 1) = -0.23500d-02
        htt( 9, 2, 1) = -0.73000d-03
        htt( 9, 3, 1) = 0.17000d-03
        htt( 9, 4, 1) = 0.17900d-02
        htt( 9, 5, 1) = 0.15900d-02
        htt( 9, 6, 1) = 0.12400d-02
        htt( 9, 7, 1) = -0.84000d-03
        htt( 9, 8, 1) = -0.11900d-02
        htt( 9, 9, 1) = 0.25000d-02
        htt(10, 1, 1) = -0.13700d-02
        htt(10, 2, 1) = 0.15000d-03
        htt(10, 3, 1) = -0.60000d-04
        htt(10, 4, 1) = 0.14600d-02
        htt(10, 5, 1) = 0.15100d-02
        htt(10, 6, 1) = 0.46000d-03
        htt(10, 7, 1) = 0.25000d-03
        htt(10, 8, 1) = -0.80000d-03
        htt(10, 9, 1) = 0.46000d-03
        htt(10,10, 1) = 0.37200d-02
        htt(11, 1, 1) = -0.36000d-03
        htt(11, 2, 1) = -0.21100d-02
        htt(11, 3, 1) = 0.23100d-02
        htt(11, 4, 1) = 0.56000d-03
        htt(11, 5, 1) = 0.29000d-03
        htt(11, 6, 1) = 0.15000d-02
        htt(11, 7, 1) = -0.22800d-02
        htt(11, 8, 1) = -0.10800d-02
        htt(11, 9, 1) = 0.83000d-03
        htt(11,10, 1) = -0.11000d-02
        htt(11,11, 1) = 0.26300d-02
        htt(12, 1, 1) = 0.62000d-03
        htt(12, 2, 1) = -0.12400d-02
        htt(12, 3, 1) = 0.20800d-02
        htt(12, 4, 1) = 0.22000d-03
        htt(12, 5, 1) = 0.21000d-03
        htt(12, 6, 1) = 0.72000d-03
        htt(12, 7, 1) = -0.11900d-02
        htt(12, 8, 1) = -0.69000d-03
        htt(12, 9, 1) = -0.12000d-02
        htt(12,10, 1) = 0.21600d-02
        htt(12,11, 1) = 0.70000d-03
        htt(12,12, 1) = 0.40700d-02
        htt(13, 1, 1) = -0.11400d-02
        htt(13, 2, 1) = -0.19000d-02
        htt(13, 3, 1) = 0.76000d-03
        htt(13, 4, 1) = -0.10000d-04
        htt(13, 5, 1) = 0.16400d-02
        htt(13, 6, 1) = 0.16200d-02
        htt(13, 7, 1) = -0.19000d-03
        htt(13, 8, 1) = -0.20000d-03
        htt(13, 9, 1) = -0.21000d-02
        htt(13,10, 1) = 0.24000d-03
        htt(13,11, 1) = -0.13100d-02
        htt(13,12, 1) = 0.10300d-02
        htt(13,13, 1) = 0.49700d-02
        htt(14, 1, 1) = 0.70000d-03
        htt(14, 2, 1) = 0.47000d-03
        htt(14, 3, 1) = -0.32000d-03
        htt(14, 4, 1) = -0.55000d-03
        htt(14, 5, 1) = -0.18000d-03
        htt(14, 6, 1) = -0.82000d-03
        htt(14, 7, 1) = 0.80000d-03
        htt(14, 8, 1) = 0.16000d-03
        htt(14, 9, 1) = -0.22000d-03
        htt(14,10, 1) = 0.13600d-02
        htt(14,11, 1) = 0.10000d-04
        htt(14,12, 1) = 0.15900d-02
        htt(14,13, 1) = -0.78000d-03
        htt(14,14, 1) = 0.38600d-02
        htt(15, 1, 1) = -0.21100d-02
        htt(15, 2, 1) = -0.27800d-02
        htt(15, 3, 1) = 0.98000d-03
        htt(15, 4, 1) = 0.32000d-03
        htt(15, 5, 1) = 0.17100d-02
        htt(15, 6, 1) = 0.23900d-02
        htt(15, 7, 1) = -0.12700d-02
        htt(15, 8, 1) = -0.59000d-03
        htt(15, 9, 1) = -0.80000d-04
        htt(15,10, 1) = -0.29900d-02
        htt(15,11, 1) = 0.61000d-03
        htt(15,12, 1) = -0.23100d-02
        htt(15,13, 1) = 0.26500d-02
        htt(15,14, 1) = -0.23400d-02
        htt(15,15, 1) = 0.55400d-02
        htt(16, 1, 1) = -0.27000d-03
        htt(16, 2, 1) = -0.40000d-03
        htt(16, 3, 1) = -0.90000d-04
        htt(16, 4, 1) = -0.22000d-03
        htt(16, 5, 1) = -0.10000d-03
        htt(16, 6, 1) = -0.50000d-04
        htt(16, 7, 1) = -0.28000d-03
        htt(16, 8, 1) = -0.22000d-03
        htt(16, 9, 1) = 0.18000d-02
        htt(16,10, 1) = -0.18700d-02
        htt(16,11, 1) = 0.19300d-02
        htt(16,12, 1) = -0.17400d-02
        htt(16,13, 1) = -0.30900d-02
        htt(16,14, 1) = 0.23000d-02
        htt(16,15, 1) = 0.55000d-03
        htt(16,16, 1) = 0.59400d-02
        htt(17, 1, 1) = 0.13400d-02
        htt(17, 2, 1) = 0.14100d-02
        htt(17, 3, 1) = -0.45000d-03
        htt(17, 4, 1) = -0.38000d-03
        htt(17, 5, 1) = -0.94000d-03
        htt(17, 6, 1) = -0.20500d-02
        htt(17, 7, 1) = 0.78000d-03
        htt(17, 8, 1) = -0.32000d-03
        htt(17, 9, 1) = 0.13000d-02
        htt(17,10, 1) = -0.10000d-03
        htt(17,11, 1) = 0.12300d-02
        htt(17,12, 1) = -0.17000d-03
        htt(17,13, 1) = -0.44500d-02
        htt(17,14, 1) = 0.13100d-02
        htt(17,15, 1) = -0.30600d-02
        htt(17,16, 1) = 0.27000d-02
        htt(17,17, 1) = 0.48400d-02
        htt(18, 1, 1) = 0.14900d-02
        htt(18, 2, 1) = 0.21600d-02
        htt(18, 3, 1) = -0.11900d-02
        htt(18, 4, 1) = -0.51000d-03
        htt(18, 5, 1) = -0.18100d-02
        htt(18, 6, 1) = -0.28100d-02
        htt(18, 7, 1) = 0.76000d-03
        htt(18, 8, 1) = -0.23000d-03
        htt(18, 9, 1) = 0.98000d-03
        htt(18,10, 1) = 0.75000d-03
        htt(18,11, 1) = 0.29000d-03
        htt(18,12, 1) = 0.60000d-04
        htt(18,13, 1) = -0.31600d-02
        htt(18,14, 1) = 0.18600d-02
        htt(18,15, 1) = -0.29300d-02
        htt(18,16, 1) = 0.20900d-02
        htt(18,17, 1) = 0.28200d-02
        htt(18,18, 1) = 0.55200d-02
        htt(19, 1, 1) = -0.50000d-03
        htt(19, 2, 1) = -0.97000d-03
        htt(19, 3, 1) = 0.63000d-03
        htt(19, 4, 1) = 0.15000d-03
        htt(19, 5, 1) = 0.88000d-03
        htt(19, 6, 1) = 0.40000d-03
        htt(19, 7, 1) = -0.21000d-03
        htt(19, 8, 1) = -0.69000d-03
        htt(19, 9, 1) = -0.59000d-03
        htt(19,10, 1) = -0.12200d-02
        htt(19,11, 1) = -0.10000d-03
        htt(19,12, 1) = -0.74000d-03
        htt(19,13, 1) = 0.13200d-02
        htt(19,14, 1) = -0.33400d-02
        htt(19,15, 1) = 0.19500d-02
        htt(19,16, 1) = -0.27100d-02
        htt(19,17, 1) = -0.93000d-03
        htt(19,18, 1) = -0.22100d-02
        htt(19,19, 1) = 0.37400d-02
        htt(20, 1, 1) = -0.35000d-03
        htt(20, 2, 1) = -0.22000d-03
        htt(20, 3, 1) = -0.11000d-03
        htt(20, 4, 1) = 0.20000d-04
        htt(20, 5, 1) = 0.10000d-04
        htt(20, 6, 1) = -0.36000d-03
        htt(20, 7, 1) = -0.23000d-03
        htt(20, 8, 1) = -0.60000d-03
        htt(20, 9, 1) = -0.91000d-03
        htt(20,10, 1) = -0.38000d-03
        htt(20,11, 1) = -0.10400d-02
        htt(20,12, 1) = -0.51000d-03
        htt(20,13, 1) = 0.26000d-02
        htt(20,14, 1) = -0.27900d-02
        htt(20,15, 1) = 0.20800d-02
        htt(20,16, 1) = -0.33200d-02
        htt(20,17, 1) = -0.29500d-02
        htt(20,18, 1) = 0.49000d-03
        htt(20,19, 1) = 0.24600d-02
        htt(20,20, 1) = 0.59000d-02
        htt(21, 1, 1) = 0.19800d-02
        htt(21, 2, 1) = 0.17000d-02
        htt(21, 3, 1) = -0.76000d-03
        htt(21, 4, 1) = -0.10500d-02
        htt(21, 5, 1) = -0.29500d-02
        htt(21, 6, 1) = -0.43000d-03
        htt(21, 7, 1) = -0.30000d-03
        htt(21, 8, 1) = 0.22100d-02
        htt(21, 9, 1) = -0.98000d-03
        htt(21,10, 1) = -0.73000d-03
        htt(21,11, 1) = -0.69000d-03
        htt(21,12, 1) = -0.44000d-03
        htt(21,13, 1) = 0.11100d-02
        htt(21,14, 1) = -0.81000d-03
        htt(21,15, 1) = 0.87000d-03
        htt(21,16, 1) = -0.10600d-02
        htt(21,17, 1) = -0.20800d-02
        htt(21,18, 1) = -0.34000d-03
        htt(21,19, 1) = -0.15000d-03
        htt(21,20, 1) = 0.15900d-02
        htt(21,21, 1) = 0.28700d-02
        htt(22, 1, 1) = -0.50000d-04
        htt(22, 2, 1) = 0.12500d-02
        htt(22, 3, 1) = -0.12100d-02
        htt(22, 4, 1) = 0.10000d-03
        htt(22, 5, 1) = 0.25000d-03
        htt(22, 6, 1) = -0.22800d-02
        htt(22, 7, 1) = 0.13700d-02
        htt(22, 8, 1) = -0.11700d-02
        htt(22, 9, 1) = -0.62000d-03
        htt(22,10, 1) = 0.37000d-03
        htt(22,11, 1) = -0.19600d-02
        htt(22,12, 1) = -0.97000d-03
        htt(22,13, 1) = 0.11300d-02
        htt(22,14, 1) = -0.13000d-03
        htt(22,15, 1) = 0.15000d-03
        htt(22,16, 1) = -0.11100d-02
        htt(22,17, 1) = -0.92000d-03
        htt(22,18, 1) = 0.70000d-03
        htt(22,19, 1) = 0.35000d-03
        htt(22,20, 1) = 0.19700d-02
        htt(22,21, 1) = 0.22000d-03
        htt(22,22, 1) = 0.28800d-02
        htt(23, 1, 1) = 0.18300d-02
        htt(23, 2, 1) = 0.94000d-03
        htt(23, 3, 1) = -0.20000d-04
        htt(23, 4, 1) = -0.91000d-03
        htt(23, 5, 1) = -0.20700d-02
        htt(23, 6, 1) = 0.34000d-03
        htt(23, 7, 1) = -0.28000d-03
        htt(23, 8, 1) = 0.21200d-02
        htt(23, 9, 1) = -0.66000d-03
        htt(23,10, 1) = -0.15900d-02
        htt(23,11, 1) = 0.26000d-03
        htt(23,12, 1) = -0.67000d-03
        htt(23,13, 1) = -0.18000d-03
        htt(23,14, 1) = -0.13600d-02
        htt(23,15, 1) = 0.74000d-03
        htt(23,16, 1) = -0.44000d-03
        htt(23,17, 1) = -0.50000d-04
        htt(23,18, 1) = -0.30600d-02
        htt(23,19, 1) = 0.11300d-02
        htt(23,20, 1) = -0.18700d-02
        htt(23,21, 1) = 0.11100d-02
        htt(23,22, 1) = -0.14200d-02
        htt(23,23, 1) = 0.41400d-02
        htt(24, 1, 1) = -0.20000d-03
        htt(24, 2, 1) = 0.49000d-03
        htt(24, 3, 1) = -0.47000d-03
        htt(24, 4, 1) = 0.23000d-03
        htt(24, 5, 1) = 0.11300d-02
        htt(24, 6, 1) = -0.15100d-02
        htt(24, 7, 1) = 0.13900d-02
        htt(24, 8, 1) = -0.12600d-02
        htt(24, 9, 1) = -0.29000d-03
        htt(24,10, 1) = -0.48000d-03
        htt(24,11, 1) = -0.10100d-02
        htt(24,12, 1) = -0.12000d-02
        htt(24,13, 1) = -0.17000d-03
        htt(24,14, 1) = -0.69000d-03
        htt(24,15, 1) = 0.20000d-04
        htt(24,16, 1) = -0.50000d-03
        htt(24,17, 1) = 0.11100d-02
        htt(24,18, 1) = -0.20100d-02
        htt(24,19, 1) = 0.16300d-02
        htt(24,20, 1) = -0.14900d-02
        htt(24,21, 1) = -0.15300d-02
        htt(24,22, 1) = 0.12400d-02
        htt(24,23, 1) = 0.16100d-02
        htt(24,24, 1) = 0.43900d-02
        ! Hessian elements for anchor point 1
        hbs( 1, 1, 1) = -0.30890d-01
        hbs( 1, 2, 1) = -0.30720d-01
        hbs( 1, 3, 1) = -0.30470d-01
        hbs( 1, 4, 1) = -0.89900d-02
        hbs( 1, 5, 1) = 0.47600d-02
        hbs( 1, 6, 1) = 0.39030d-01
        hbs( 1, 7, 1) = -0.44000d-03
        hbs( 1, 8, 1) = 0.39840d-01
        hbs( 1, 9, 1) = -0.25000d-03
        hbs( 1,10, 1) = -0.97500d-02
        hbs( 1,11, 1) = -0.42000d-03
        hbs( 1,12, 1) = 0.49400d-02
        hbs( 2, 1, 1) = 0.43860d-01
        hbs( 2, 2, 1) = -0.11950d-01
        hbs( 2, 3, 1) = 0.67900d-02
        hbs( 2, 4, 1) = 0.23800d-02
        hbs( 2, 5, 1) = -0.22800d-02
        hbs( 2, 6, 1) = -0.26060d-01
        hbs( 2, 7, 1) = 0.11000d-02
        hbs( 2, 8, 1) = -0.14350d-01
        hbs( 2, 9, 1) = 0.40000d-04
        hbs( 2,10, 1) = 0.41400d-02
        hbs( 2,11, 1) = -0.73000d-03
        hbs( 2,12, 1) = -0.48000d-03
        hbs( 3, 1, 1) = -0.12970d-01
        hbs( 3, 2, 1) = 0.42680d-01
        hbs( 3, 3, 1) = 0.23690d-01
        hbs( 3, 4, 1) = 0.66100d-02
        hbs( 3, 5, 1) = -0.24900d-02
        hbs( 3, 6, 1) = -0.12970d-01
        hbs( 3, 7, 1) = -0.66000d-03
        hbs( 3, 8, 1) = -0.25490d-01
        hbs( 3, 9, 1) = 0.21000d-03
        hbs( 3,10, 1) = 0.56200d-02
        hbs( 3,11, 1) = 0.11500d-02
        hbs( 3,12, 1) = -0.44600d-02
        hbs( 4, 1, 1) = -0.21190d-01
        hbs( 4, 2, 1) = -0.10770d-01
        hbs( 4, 3, 1) = 0.12660d-01
        hbs( 4, 4, 1) = -0.29880d-01
        hbs( 4, 5, 1) = -0.90800d-02
        hbs( 4, 6, 1) = -0.51700d-02
        hbs( 4, 7, 1) = 0.49300d-02
        hbs( 4, 8, 1) = 0.30850d-01
        hbs( 4, 9, 1) = -0.40000d-03
        hbs( 4,10, 1) = 0.40050d-01
        hbs( 4,11, 1) = -0.62000d-03
        hbs( 4,12, 1) = 0.12000d-03
        hbs( 5, 1, 1) = 0.25060d-01
        hbs( 5, 2, 1) = 0.46900d-02
        hbs( 5, 3, 1) = -0.97000d-02
        hbs( 5, 4, 1) = 0.24100d-02
        hbs( 5, 5, 1) = 0.62100d-02
        hbs( 5, 6, 1) = 0.30800d-02
        hbs( 5, 7, 1) = -0.21300d-02
        hbs( 5, 8, 1) = -0.14880d-01
        hbs( 5, 9, 1) = -0.39000d-03
        hbs( 5,10, 1) = -0.21940d-01
        hbs( 5,11, 1) = 0.31000d-03
        hbs( 5,12, 1) = 0.35000d-03
        hbs( 6, 1, 1) = -0.38700d-02
        hbs( 6, 2, 1) = 0.60800d-02
        hbs( 6, 3, 1) = -0.29500d-02
        hbs( 6, 4, 1) = 0.27480d-01
        hbs( 6, 5, 1) = 0.28700d-02
        hbs( 6, 6, 1) = 0.21000d-02
        hbs( 6, 7, 1) = -0.28000d-02
        hbs( 6, 8, 1) = -0.15960d-01
        hbs( 6, 9, 1) = 0.80000d-03
        hbs( 6,10, 1) = -0.18110d-01
        hbs( 6,11, 1) = 0.31000d-03
        hbs( 6,12, 1) = -0.47000d-03
        hbs( 7, 1, 1) = -0.88700d-02
        hbs( 7, 2, 1) = 0.35920d-01
        hbs( 7, 3, 1) = 0.96000d-03
        hbs( 7, 4, 1) = -0.23260d-01
        hbs( 7, 5, 1) = 0.57800d-02
        hbs( 7, 6, 1) = -0.29090d-01
        hbs( 7, 7, 1) = -0.87400d-02
        hbs( 7, 8, 1) = -0.88100d-02
        hbs( 7, 9, 1) = 0.49600d-02
        hbs( 7,10, 1) = 0.32900d-01
        hbs( 7,11, 1) = 0.60000d-04
        hbs( 7,12, 1) = -0.10800d-02
        hbs( 8, 1, 1) = 0.36600d-02
        hbs( 8, 2, 1) = -0.18680d-01
        hbs( 8, 3, 1) = 0.82000d-03
        hbs( 8, 4, 1) = 0.25710d-01
        hbs( 8, 5, 1) = -0.33200d-02
        hbs( 8, 6, 1) = -0.40000d-03
        hbs( 8, 7, 1) = 0.44400d-02
        hbs( 8, 8, 1) = 0.48100d-02
        hbs( 8, 9, 1) = -0.20200d-02
        hbs( 8,10, 1) = -0.14750d-01
        hbs( 8,11, 1) = -0.62000d-03
        hbs( 8,12, 1) = 0.56000d-03
        hbs( 9, 1, 1) = 0.52100d-02
        hbs( 9, 2, 1) = -0.17240d-01
        hbs( 9, 3, 1) = -0.17800d-02
        hbs( 9, 4, 1) = -0.24500d-02
        hbs( 9, 5, 1) = -0.24600d-02
        hbs( 9, 6, 1) = 0.29490d-01
        hbs( 9, 7, 1) = 0.43000d-02
        hbs( 9, 8, 1) = 0.40000d-02
        hbs( 9, 9, 1) = -0.29400d-02
        hbs( 9,10, 1) = -0.18150d-01
        hbs( 9,11, 1) = 0.56000d-03
        hbs( 9,12, 1) = 0.52000d-03
        hbs(10, 1, 1) = 0.36080d-01
        hbs(10, 2, 1) = 0.38610d-01
        hbs(10, 3, 1) = 0.39500d-02
        hbs(10, 4, 1) = -0.11320d-01
        hbs(10, 5, 1) = -0.89000d-03
        hbs(10, 6, 1) = -0.26740d-01
        hbs(10, 7, 1) = 0.48000d-02
        hbs(10, 8, 1) = -0.26910d-01
        hbs(10, 9, 1) = -0.88100d-02
        hbs(10,10, 1) = -0.11650d-01
        hbs(10,11, 1) = 0.47600d-02
        hbs(10,12, 1) = -0.70000d-03
        hbs(11, 1, 1) = -0.20240d-01
        hbs(11, 2, 1) = -0.17250d-01
        hbs(11, 3, 1) = -0.21300d-02
        hbs(11, 4, 1) = 0.55100d-02
        hbs(11, 5, 1) = 0.11200d-02
        hbs(11, 6, 1) = 0.27050d-01
        hbs(11, 7, 1) = -0.28200d-02
        hbs(11, 8, 1) = 0.25000d-03
        hbs(11, 9, 1) = 0.44500d-02
        hbs(11,10, 1) = 0.58600d-02
        hbs(11,11, 1) = -0.19800d-02
        hbs(11,12, 1) = -0.20000d-03
        hbs(12, 1, 1) = -0.15840d-01
        hbs(12, 2, 1) = -0.21360d-01
        hbs(12, 3, 1) = -0.18200d-02
        hbs(12, 4, 1) = 0.58100d-02
        hbs(12, 5, 1) = -0.23000d-03
        hbs(12, 6, 1) = -0.31000d-03
        hbs(12, 7, 1) = -0.19800d-02
        hbs(12, 8, 1) = 0.26660d-01
        hbs(12, 9, 1) = 0.43600d-02
        hbs(12,10, 1) = 0.57900d-02
        hbs(12,11, 1) = -0.27800d-02
        hbs(12,12, 1) = 0.90000d-03
        hbs(13, 1, 1) = 0.35090d-01
        hbs(13, 2, 1) = -0.85500d-02
        hbs(13, 3, 1) = 0.10200d-02
        hbs(13, 4, 1) = 0.32180d-01
        hbs(13, 5, 1) = -0.97000d-03
        hbs(13, 6, 1) = -0.94600d-02
        hbs(13, 7, 1) = 0.10000d-03
        hbs(13, 8, 1) = -0.29520d-01
        hbs(13, 9, 1) = 0.50200d-02
        hbs(13,10, 1) = -0.23910d-01
        hbs(13,11, 1) = -0.86300d-02
        hbs(13,12, 1) = 0.52600d-02
        hbs(14, 1, 1) = -0.16950d-01
        hbs(14, 2, 1) = 0.51900d-02
        hbs(14, 3, 1) = -0.19500d-02
        hbs(14, 4, 1) = -0.17590d-01
        hbs(14, 5, 1) = 0.48000d-03
        hbs(14, 6, 1) = 0.43300d-02
        hbs(14, 7, 1) = 0.54000d-03
        hbs(14, 8, 1) = 0.29620d-01
        hbs(14, 9, 1) = -0.29700d-02
        hbs(14,10, 1) = -0.21100d-02
        hbs(14,11, 1) = 0.42400d-02
        hbs(14,12, 1) = -0.22400d-02
        hbs(15, 1, 1) = -0.18140d-01
        hbs(15, 2, 1) = 0.33600d-02
        hbs(15, 3, 1) = 0.94000d-03
        hbs(15, 4, 1) = -0.14580d-01
        hbs(15, 5, 1) = 0.48000d-03
        hbs(15, 6, 1) = 0.51200d-02
        hbs(15, 7, 1) = -0.64000d-03
        hbs(15, 8, 1) = -0.10000d-03
        hbs(15, 9, 1) = -0.20500d-02
        hbs(15,10, 1) = 0.26020d-01
        hbs(15,11, 1) = 0.43900d-02
        hbs(15,12, 1) = -0.30200d-02
        hbs(16, 1, 1) = -0.10210d-01
        hbs(16, 2, 1) = -0.24490d-01
        hbs(16, 3, 1) = 0.11890d-01
        hbs(16, 4, 1) = 0.41280d-01
        hbs(16, 5, 1) = 0.40000d-03
        hbs(16, 6, 1) = 0.31430d-01
        hbs(16, 7, 1) = -0.65000d-03
        hbs(16, 8, 1) = -0.54400d-02
        hbs(16, 9, 1) = -0.51000d-03
        hbs(16,10, 1) = -0.27630d-01
        hbs(16,11, 1) = 0.48400d-02
        hbs(16,12, 1) = -0.85500d-02
        hbs(17, 1, 1) = 0.43200d-02
        hbs(17, 2, 1) = 0.25240d-01
        hbs(17, 3, 1) = -0.81200d-02
        hbs(17, 4, 1) = -0.22630d-01
        hbs(17, 5, 1) = 0.41000d-03
        hbs(17, 6, 1) = -0.14630d-01
        hbs(17, 7, 1) = 0.33000d-03
        hbs(17, 8, 1) = 0.32400d-02
        hbs(17, 9, 1) = -0.33000d-03
        hbs(17,10, 1) = 0.18800d-02
        hbs(17,11, 1) = -0.21000d-02
        hbs(17,12, 1) = 0.69700d-02
        hbs(18, 1, 1) = 0.58800d-02
        hbs(18, 2, 1) = -0.75000d-03
        hbs(18, 3, 1) = -0.37700d-02
        hbs(18, 4, 1) = -0.18650d-01
        hbs(18, 5, 1) = -0.80000d-03
        hbs(18, 6, 1) = -0.16800d-01
        hbs(18, 7, 1) = 0.32000d-03
        hbs(18, 8, 1) = 0.22000d-02
        hbs(18, 9, 1) = 0.85000d-03
        hbs(18,10, 1) = 0.25750d-01
        hbs(18,11, 1) = -0.27400d-02
        hbs(18,12, 1) = 0.15700d-02
        ! Hessian elements for anchor point 2
        hss( 1, 1, 2) = 0.34759d+00
        hss( 2, 1, 2) = 0.17810d-01
        hss( 2, 2, 2) = 0.35623d+00
        hss( 3, 1, 2) = 0.29010d-01
        hss( 3, 2, 2) = 0.32710d-01
        hss( 3, 3, 2) = 0.46124d+00
        hss( 4, 1, 2) = -0.19600d-02
        hss( 4, 2, 2) = 0.14870d-01
        hss( 4, 3, 2) = -0.89100d-02
        hss( 4, 4, 2) = 0.38581d+00
        hss( 5, 1, 2) = 0.38800d-02
        hss( 5, 2, 2) = -0.77000d-03
        hss( 5, 3, 2) = 0.71700d-02
        hss( 5, 4, 2) = 0.36500d-02
        hss( 5, 5, 2) = 0.35807d+00
        hss( 6, 1, 2) = 0.33890d-01
        hss( 6, 2, 2) = 0.44570d-01
        hss( 6, 3, 2) = -0.43100d-02
        hss( 6, 4, 2) = -0.36100d-02
        hss( 6, 5, 2) = 0.55000d-03
        hss( 6, 6, 2) = 0.35931d+00
        hss( 7, 1, 2) = -0.54000d-03
        hss( 7, 2, 2) = -0.63000d-03
        hss( 7, 3, 2) = 0.57000d-03
        hss( 7, 4, 2) = 0.39700d-02
        hss( 7, 5, 2) = 0.39000d-03
        hss( 7, 6, 2) = 0.39800d-02
        hss( 7, 7, 2) = 0.37091d+00
        hss( 8, 1, 2) = 0.42030d-01
        hss( 8, 2, 2) = 0.35060d-01
        hss( 8, 3, 2) = -0.52800d-02
        hss( 8, 4, 2) = 0.24280d-01
        hss( 8, 5, 2) = -0.11100d-02
        hss( 8, 6, 2) = 0.17880d-01
        hss( 8, 7, 2) = -0.10400d-02
        hss( 8, 8, 2) = 0.35888d+00
        hss( 9, 1, 2) = -0.53000d-03
        hss( 9, 2, 2) = -0.57000d-03
        hss( 9, 3, 2) = 0.26000d-03
        hss( 9, 4, 2) = -0.39000d-03
        hss( 9, 5, 2) = 0.19000d-03
        hss( 9, 6, 2) = 0.40100d-02
        hss( 9, 7, 2) = 0.49000d-03
        hss( 9, 8, 2) = 0.39600d-02
        hss( 9, 9, 2) = 0.36356d+00
        hss(10, 1, 2) = 0.19520d-01
        hss(10, 2, 2) = -0.69300d-02
        hss(10, 3, 2) = -0.48800d-02
        hss(10, 4, 2) = 0.64010d-01
        hss(10, 5, 2) = -0.43000d-03
        hss(10, 6, 2) = 0.21230d-01
        hss(10, 7, 2) = -0.10900d-02
        hss(10, 8, 2) = -0.77000d-03
        hss(10, 9, 2) = -0.19000d-03
        hss(10,10, 2) = 0.38116d+00
        hss(11, 1, 2) = -0.69000d-03
        hss(11, 2, 2) = -0.50000d-03
        hss(11, 3, 2) = 0.52000d-03
        hss(11, 4, 2) = -0.11600d-02
        hss(11, 5, 2) = 0.80000d-04
        hss(11, 6, 2) = -0.10000d-02
        hss(11, 7, 2) = 0.22000d-03
        hss(11, 8, 2) = 0.39700d-02
        hss(11, 9, 2) = 0.48000d-03
        hss(11,10, 2) = 0.39800d-02
        hss(11,11, 2) = 0.37152d+00
        hss(12, 1, 2) = -0.16500d-02
        hss(12, 2, 2) = 0.39400d-02
        hss(12, 3, 2) = -0.19000d-03
        hss(12, 4, 2) = -0.10000d-03
        hss(12, 5, 2) = 0.35000d-03
        hss(12, 6, 2) = -0.10400d-02
        hss(12, 7, 2) = 0.30000d-04
        hss(12, 8, 2) = 0.57000d-03
        hss(12, 9, 2) = 0.12000d-03
        hss(12,10, 2) = 0.23600d-02
        hss(12,11, 2) = 0.40000d-03
        hss(12,12, 2) = 0.37254d+00
        ! Hessian elements for anchor point 2
        hbb( 1, 1, 2) = 0.82450d-01
        hbb( 2, 1, 2) = -0.46090d-01
        hbb( 2, 2, 2) = 0.14487d+00
        hbb( 3, 1, 2) = -0.36360d-01
        hbb( 3, 2, 2) = -0.98780d-01
        hbb( 3, 3, 2) = 0.13515d+00
        hbb( 4, 1, 2) = -0.36000d-01
        hbb( 4, 2, 2) = 0.23310d-01
        hbb( 4, 3, 2) = 0.12690d-01
        hbb( 4, 4, 2) = 0.66170d-01
        hbb( 5, 1, 2) = 0.20240d-01
        hbb( 5, 2, 2) = -0.12480d-01
        hbb( 5, 3, 2) = -0.77700d-02
        hbb( 5, 4, 2) = -0.33030d-01
        hbb( 5, 5, 2) = 0.71950d-01
        hbb( 6, 1, 2) = 0.15760d-01
        hbb( 6, 2, 2) = -0.10840d-01
        hbb( 6, 3, 2) = -0.49200d-02
        hbb( 6, 4, 2) = -0.33140d-01
        hbb( 6, 5, 2) = -0.38920d-01
        hbb( 6, 6, 2) = 0.72060d-01
        hbb( 7, 1, 2) = -0.71500d-02
        hbb( 7, 2, 2) = -0.39000d-02
        hbb( 7, 3, 2) = 0.11050d-01
        hbb( 7, 4, 2) = -0.25740d-01
        hbb( 7, 5, 2) = 0.93100d-02
        hbb( 7, 6, 2) = 0.16430d-01
        hbb( 7, 7, 2) = 0.64420d-01
        hbb( 8, 1, 2) = -0.22800d-02
        hbb( 8, 2, 2) = 0.74300d-02
        hbb( 8, 3, 2) = -0.51500d-02
        hbb( 8, 4, 2) = 0.15690d-01
        hbb( 8, 5, 2) = -0.49500d-02
        hbb( 8, 6, 2) = -0.10740d-01
        hbb( 8, 7, 2) = -0.31490d-01
        hbb( 8, 8, 2) = 0.72950d-01
        hbb( 9, 1, 2) = 0.94300d-02
        hbb( 9, 2, 2) = -0.35300d-02
        hbb( 9, 3, 2) = -0.59100d-02
        hbb( 9, 4, 2) = 0.10050d-01
        hbb( 9, 5, 2) = -0.43500d-02
        hbb( 9, 6, 2) = -0.56900d-02
        hbb( 9, 7, 2) = -0.32930d-01
        hbb( 9, 8, 2) = -0.41460d-01
        hbb( 9, 9, 2) = 0.74390d-01
        hbb(10, 1, 2) = 0.45500d-02
        hbb(10, 2, 2) = -0.20200d-02
        hbb(10, 3, 2) = -0.25300d-02
        hbb(10, 4, 2) = -0.64100d-02
        hbb(10, 5, 2) = 0.91900d-02
        hbb(10, 6, 2) = -0.27700d-02
        hbb(10, 7, 2) = -0.33410d-01
        hbb(10, 8, 2) = 0.13410d-01
        hbb(10, 9, 2) = 0.20000d-01
        hbb(10,10, 2) = 0.75440d-01
        hbb(11, 1, 2) = -0.22100d-02
        hbb(11, 2, 2) = 0.12300d-02
        hbb(11, 3, 2) = 0.98000d-03
        hbb(11, 4, 2) = -0.23100d-02
        hbb(11, 5, 2) = -0.32200d-02
        hbb(11, 6, 2) = 0.55300d-02
        hbb(11, 7, 2) = 0.20020d-01
        hbb(11, 8, 2) = -0.74100d-02
        hbb(11, 9, 2) = -0.12610d-01
        hbb(11,10, 2) = -0.37790d-01
        hbb(11,11, 2) = 0.78360d-01
        hbb(12, 1, 2) = -0.23400d-02
        hbb(12, 2, 2) = 0.79000d-03
        hbb(12, 3, 2) = 0.15500d-02
        hbb(12, 4, 2) = 0.87200d-02
        hbb(12, 5, 2) = -0.59700d-02
        hbb(12, 6, 2) = -0.27500d-02
        hbb(12, 7, 2) = 0.13390d-01
        hbb(12, 8, 2) = -0.60000d-02
        hbb(12, 9, 2) = -0.73900d-02
        hbb(12,10, 2) = -0.37640d-01
        hbb(12,11, 2) = -0.40570d-01
        hbb(12,12, 2) = 0.78220d-01
        hbb(13, 1, 2) = -0.67600d-02
        hbb(13, 2, 2) = 0.10710d-01
        hbb(13, 3, 2) = -0.39500d-02
        hbb(13, 4, 2) = -0.11500d-02
        hbb(13, 5, 2) = 0.12900d-02
        hbb(13, 6, 2) = -0.14000d-03
        hbb(13, 7, 2) = 0.29100d-02
        hbb(13, 8, 2) = 0.41300d-02
        hbb(13, 9, 2) = -0.70300d-02
        hbb(13,10, 2) = -0.33760d-01
        hbb(13,11, 2) = 0.13510d-01
        hbb(13,12, 2) = 0.20240d-01
        hbb(13,13, 2) = 0.65300d-01
        hbb(14, 1, 2) = 0.91200d-02
        hbb(14, 2, 2) = -0.60500d-02
        hbb(14, 3, 2) = -0.30700d-02
        hbb(14, 4, 2) = 0.66000d-03
        hbb(14, 5, 2) = -0.85000d-03
        hbb(14, 6, 2) = 0.19000d-03
        hbb(14, 7, 2) = -0.71100d-02
        hbb(14, 8, 2) = -0.51000d-03
        hbb(14, 9, 2) = 0.76200d-02
        hbb(14,10, 2) = 0.20030d-01
        hbb(14,11, 2) = -0.73800d-02
        hbb(14,12, 2) = -0.12650d-01
        hbb(14,13, 2) = -0.33180d-01
        hbb(14,14, 2) = 0.74350d-01
        hbb(15, 1, 2) = -0.23600d-02
        hbb(15, 2, 2) = -0.46500d-02
        hbb(15, 3, 2) = 0.70200d-02
        hbb(15, 4, 2) = 0.49000d-03
        hbb(15, 5, 2) = -0.44000d-03
        hbb(15, 6, 2) = -0.50000d-04
        hbb(15, 7, 2) = 0.42000d-02
        hbb(15, 8, 2) = -0.36200d-02
        hbb(15, 9, 2) = -0.59000d-03
        hbb(15,10, 2) = 0.13730d-01
        hbb(15,11, 2) = -0.61300d-02
        hbb(15,12, 2) = -0.76000d-02
        hbb(15,13, 2) = -0.32120d-01
        hbb(15,14, 2) = -0.41170d-01
        hbb(15,15, 2) = 0.73290d-01
        hbb(16, 1, 2) = -0.37090d-01
        hbb(16, 2, 2) = 0.17990d-01
        hbb(16, 3, 2) = 0.19100d-01
        hbb(16, 4, 2) = 0.31400d-02
        hbb(16, 5, 2) = -0.70000d-02
        hbb(16, 6, 2) = 0.38600d-02
        hbb(16, 7, 2) = -0.10300d-02
        hbb(16, 8, 2) = 0.54000d-03
        hbb(16, 9, 2) = 0.48000d-03
        hbb(16,10, 2) = -0.64000d-02
        hbb(16,11, 2) = 0.87800d-02
        hbb(16,12, 2) = -0.23700d-02
        hbb(16,13, 2) = -0.26530d-01
        hbb(16,14, 2) = 0.10480d-01
        hbb(16,15, 2) = 0.16060d-01
        hbb(16,16, 2) = 0.67910d-01
        hbb(17, 1, 2) = 0.22610d-01
        hbb(17, 2, 2) = -0.90200d-02
        hbb(17, 3, 2) = -0.13590d-01
        hbb(17, 4, 2) = -0.76000d-02
        hbb(17, 5, 2) = 0.77900d-02
        hbb(17, 6, 2) = -0.19000d-03
        hbb(17, 7, 2) = 0.87000d-03
        hbb(17, 8, 2) = -0.25000d-03
        hbb(17, 9, 2) = -0.62000d-03
        hbb(17,10, 2) = 0.91600d-02
        hbb(17,11, 2) = -0.60700d-02
        hbb(17,12, 2) = -0.31000d-02
        hbb(17,13, 2) = 0.96800d-02
        hbb(17,14, 2) = -0.46000d-02
        hbb(17,15, 2) = -0.50800d-02
        hbb(17,16, 2) = -0.34720d-01
        hbb(17,17, 2) = 0.72390d-01
        hbb(18, 1, 2) = 0.14470d-01
        hbb(18, 2, 2) = -0.89700d-02
        hbb(18, 3, 2) = -0.55000d-02
        hbb(18, 4, 2) = 0.44600d-02
        hbb(18, 5, 2) = -0.79000d-03
        hbb(18, 6, 2) = -0.36700d-02
        hbb(18, 7, 2) = 0.16000d-03
        hbb(18, 8, 2) = -0.29000d-03
        hbb(18, 9, 2) = 0.13000d-03
        hbb(18,10, 2) = -0.27600d-02
        hbb(18,11, 2) = -0.27100d-02
        hbb(18,12, 2) = 0.54700d-02
        hbb(18,13, 2) = 0.16860d-01
        hbb(18,14, 2) = -0.58800d-02
        hbb(18,15, 2) = -0.10980d-01
        hbb(18,16, 2) = -0.33190d-01
        hbb(18,17, 2) = -0.37670d-01
        hbb(18,18, 2) = 0.70860d-01
        ! Hessian elements for anchor point 2
        htt( 1, 1, 2) = 0.60400d-02
        htt( 2, 1, 2) = 0.38400d-02
        htt( 2, 2, 2) = 0.50900d-02
        htt( 3, 1, 2) = -0.18500d-02
        htt( 3, 2, 2) = -0.33100d-02
        htt( 3, 3, 2) = 0.75100d-02
        htt( 4, 1, 2) = -0.40500d-02
        htt( 4, 2, 2) = -0.20700d-02
        htt( 4, 3, 2) = 0.60500d-02
        htt( 4, 4, 2) = 0.80200d-02
        htt( 5, 1, 2) = -0.59200d-02
        htt( 5, 2, 2) = -0.42600d-02
        htt( 5, 3, 2) = 0.21200d-02
        htt( 5, 4, 2) = 0.37800d-02
        htt( 5, 5, 2) = 0.65600d-02
        htt( 6, 1, 2) = -0.41200d-02
        htt( 6, 2, 2) = -0.39000d-02
        htt( 6, 3, 2) = 0.26900d-02
        htt( 6, 4, 2) = 0.29100d-02
        htt( 6, 5, 2) = 0.37400d-02
        htt( 6, 6, 2) = 0.55900d-02
        htt( 7, 1, 2) = 0.15700d-02
        htt( 7, 2, 2) = 0.25400d-02
        htt( 7, 3, 2) = -0.67700d-02
        htt( 7, 4, 2) = -0.58000d-02
        htt( 7, 5, 2) = -0.10700d-02
        htt( 7, 6, 2) = -0.27200d-02
        htt( 7, 7, 2) = 0.68600d-02
        htt( 8, 1, 2) = 0.33800d-02
        htt( 8, 2, 2) = 0.29000d-02
        htt( 8, 3, 2) = -0.62000d-02
        htt( 8, 4, 2) = -0.66800d-02
        htt( 8, 5, 2) = -0.38900d-02
        htt( 8, 6, 2) = -0.88000d-03
        htt( 8, 7, 2) = 0.52100d-02
        htt( 8, 8, 2) = 0.82200d-02
        htt( 9, 1, 2) = -0.23900d-02
        htt( 9, 2, 2) = -0.68000d-03
        htt( 9, 3, 2) = 0.45000d-03
        htt( 9, 4, 2) = 0.21600d-02
        htt( 9, 5, 2) = 0.17400d-02
        htt( 9, 6, 2) = 0.14100d-02
        htt( 9, 7, 2) = -0.96000d-03
        htt( 9, 8, 2) = -0.12900d-02
        htt( 9, 9, 2) = 0.24600d-02
        htt(10, 1, 2) = -0.15600d-02
        htt(10, 2, 2) = 0.23000d-03
        htt(10, 3, 2) = 0.34000d-03
        htt(10, 4, 2) = 0.21300d-02
        htt(10, 5, 2) = 0.17600d-02
        htt(10, 6, 2) = 0.79000d-03
        htt(10, 7, 2) = -0.40000d-04
        htt(10, 8, 2) = -0.10100d-02
        htt(10, 9, 2) = 0.57000d-03
        htt(10,10, 2) = 0.39100d-02
        htt(11, 1, 2) = -0.10000d-03
        htt(11, 2, 2) = -0.19700d-02
        htt(11, 3, 2) = 0.19800d-02
        htt(11, 4, 2) = 0.11000d-03
        htt(11, 5, 2) = 0.10000d-04
        htt(11, 6, 2) = 0.11900d-02
        htt(11, 7, 2) = -0.19700d-02
        htt(11, 8, 2) = -0.79000d-03
        htt(11, 9, 2) = 0.68000d-03
        htt(11,10, 2) = -0.12900d-02
        htt(11,11, 2) = 0.26200d-02
        htt(12, 1, 2) = 0.72000d-03
        htt(12, 2, 2) = -0.10500d-02
        htt(12, 3, 2) = 0.18600d-02
        htt(12, 4, 2) = 0.90000d-04
        htt(12, 5, 2) = 0.30000d-04
        htt(12, 6, 2) = 0.57000d-03
        htt(12, 7, 2) = -0.10500d-02
        htt(12, 8, 2) = -0.51000d-03
        htt(12, 9, 2) = -0.12000d-02
        htt(12,10, 2) = 0.20500d-02
        htt(12,11, 2) = 0.64000d-03
        htt(12,12, 2) = 0.39000d-02
        htt(13, 1, 2) = -0.10900d-02
        htt(13, 2, 2) = -0.19400d-02
        htt(13, 3, 2) = 0.59000d-03
        htt(13, 4, 2) = -0.26000d-03
        htt(13, 5, 2) = 0.15200d-02
        htt(13, 6, 2) = 0.15100d-02
        htt(13, 7, 2) = -0.70000d-04
        htt(13, 8, 2) = -0.90000d-04
        htt(13, 9, 2) = -0.21100d-02
        htt(13,10, 2) = 0.12000d-03
        htt(13,11, 2) = -0.12300d-02
        htt(13,12, 2) = 0.10000d-02
        htt(13,13, 2) = 0.50300d-02
        htt(14, 1, 2) = 0.70000d-03
        htt(14, 2, 2) = 0.53000d-03
        htt(14, 3, 2) = -0.31000d-03
        htt(14, 4, 2) = -0.47000d-03
        htt(14, 5, 2) = -0.22000d-03
        htt(14, 6, 2) = -0.80000d-03
        htt(14, 7, 2) = 0.74000d-03
        htt(14, 8, 2) = 0.16000d-03
        htt(14, 9, 2) = -0.19000d-03
        htt(14,10, 2) = 0.13500d-02
        htt(14,11, 2) = -0.10000d-04
        htt(14,12, 2) = 0.15300d-02
        htt(14,13, 2) = -0.81000d-03
        htt(14,14, 2) = 0.39200d-02
        htt(15, 1, 2) = -0.19100d-02
        htt(15, 2, 2) = -0.28400d-02
        htt(15, 3, 2) = 0.70000d-03
        htt(15, 4, 2) = -0.23000d-03
        htt(15, 5, 2) = 0.15000d-02
        htt(15, 6, 2) = 0.21100d-02
        htt(15, 7, 2) = -0.98000d-03
        htt(15, 8, 2) = -0.37000d-03
        htt(15, 9, 2) = -0.24000d-03
        htt(15,10, 2) = -0.31800d-02
        htt(15,11, 2) = 0.73000d-03
        htt(15,12, 2) = -0.22200d-02
        htt(15,13, 2) = 0.28200d-02
        htt(15,14, 2) = -0.23300d-02
        htt(15,15, 2) = 0.57300d-02
        htt(16, 1, 2) = -0.12000d-03
        htt(16, 2, 2) = -0.37000d-03
        htt(16, 3, 2) = -0.19000d-03
        htt(16, 4, 2) = -0.45000d-03
        htt(16, 5, 2) = -0.25000d-03
        htt(16, 6, 2) = -0.19000d-03
        htt(16, 7, 2) = -0.17000d-03
        htt(16, 8, 2) = -0.12000d-03
        htt(16, 9, 2) = 0.16800d-02
        htt(16,10, 2) = -0.19500d-02
        htt(16,11, 2) = 0.19400d-02
        htt(16,12, 2) = -0.16900d-02
        htt(16,13, 2) = -0.30200d-02
        htt(16,14, 2) = 0.24000d-02
        htt(16,15, 2) = 0.58000d-03
        htt(16,16, 2) = 0.60000d-02
        htt(17, 1, 2) = 0.12100d-02
        htt(17, 2, 2) = 0.15200d-02
        htt(17, 3, 2) = -0.31000d-03
        htt(17, 4, 2) = 0.00000d+00
        htt(17, 5, 2) = -0.86000d-03
        htt(17, 6, 2) = -0.18800d-02
        htt(17, 7, 2) = 0.58000d-03
        htt(17, 8, 2) = -0.44000d-03
        htt(17, 9, 2) = 0.14500d-02
        htt(17,10, 2) = 0.90000d-04
        htt(17,11, 2) = 0.11300d-02
        htt(17,12, 2) = -0.23000d-03
        htt(17,13, 2) = -0.46000d-02
        htt(17,14, 2) = 0.13000d-02
        htt(17,15, 2) = -0.32400d-02
        htt(17,16, 2) = 0.26500d-02
        htt(17,17, 2) = 0.49600d-02
        htt(18, 1, 2) = 0.14500d-02
        htt(18, 2, 2) = 0.22200d-02
        htt(18, 3, 2) = -0.10600d-02
        htt(18, 4, 2) = -0.29000d-03
        htt(18, 5, 2) = -0.17100d-02
        htt(18, 6, 2) = -0.27100d-02
        htt(18, 7, 2) = 0.67000d-03
        htt(18, 8, 2) = -0.32000d-03
        htt(18, 9, 2) = 0.10000d-02
        htt(18,10, 2) = 0.82000d-03
        htt(18,11, 2) = 0.19000d-03
        htt(18,12, 2) = 0.10000d-04
        htt(18,13, 2) = -0.32400d-02
        htt(18,14, 2) = 0.19000d-02
        htt(18,15, 2) = -0.30600d-02
        htt(18,16, 2) = 0.20800d-02
        htt(18,17, 2) = 0.29700d-02
        htt(18,18, 2) = 0.56300d-02
        htt(19, 1, 2) = -0.59000d-03
        htt(19, 2, 2) = -0.96000d-03
        htt(19, 3, 2) = 0.59000d-03
        htt(19, 4, 2) = 0.21000d-03
        htt(19, 5, 2) = 0.89000d-03
        htt(19, 6, 2) = 0.43000d-03
        htt(19, 7, 2) = -0.23000d-03
        htt(19, 8, 2) = -0.69000d-03
        htt(19, 9, 2) = -0.48000d-03
        htt(19,10, 2) = -0.11500d-02
        htt(19,11, 2) = -0.90000d-04
        htt(19,12, 2) = -0.76000d-03
        htt(19,13, 2) = 0.12600d-02
        htt(19,14, 2) = -0.34500d-02
        htt(19,15, 2) = 0.19200d-02
        htt(19,16, 2) = -0.27800d-02
        htt(19,17, 2) = -0.95000d-03
        htt(19,18, 2) = -0.21800d-02
        htt(19,19, 2) = 0.37700d-02
        htt(20, 1, 2) = -0.35000d-03
        htt(20, 2, 2) = -0.26000d-03
        htt(20, 3, 2) = -0.16000d-03
        htt(20, 4, 2) = -0.80000d-04
        htt(20, 5, 2) = 0.40000d-04
        htt(20, 6, 2) = -0.40000d-03
        htt(20, 7, 2) = -0.14000d-03
        htt(20, 8, 2) = -0.57000d-03
        htt(20, 9, 2) = -0.93000d-03
        htt(20,10, 2) = -0.42000d-03
        htt(20,11, 2) = -0.10200d-02
        htt(20,12, 2) = -0.51000d-03
        htt(20,13, 2) = 0.26100d-02
        htt(20,14, 2) = -0.28400d-02
        htt(20,15, 2) = 0.21100d-02
        htt(20,16, 2) = -0.33500d-02
        htt(20,17, 2) = -0.29400d-02
        htt(20,18, 2) = 0.47000d-03
        htt(20,19, 2) = 0.25400d-02
        htt(20,20, 2) = 0.59400d-02
        htt(21, 1, 2) = 0.21700d-02
        htt(21, 2, 2) = 0.15000d-02
        htt(21, 3, 2) = -0.99000d-03
        htt(21, 4, 2) = -0.16600d-02
        htt(21, 5, 2) = -0.30200d-02
        htt(21, 6, 2) = -0.69000d-03
        htt(21, 7, 2) = -0.20000d-04
        htt(21, 8, 2) = 0.23100d-02
        htt(21, 9, 2) = -0.11900d-02
        htt(21,10, 2) = -0.98000d-03
        htt(21,11, 2) = -0.50000d-03
        htt(21,12, 2) = -0.28000d-03
        htt(21,13, 2) = 0.12600d-02
        htt(21,14, 2) = -0.76000d-03
        htt(21,15, 2) = 0.10500d-02
        htt(21,16, 2) = -0.97000d-03
        htt(21,17, 2) = -0.21400d-02
        htt(21,18, 2) = -0.48000d-03
        htt(21,19, 2) = -0.11000d-03
        htt(21,20, 2) = 0.15500d-02
        htt(21,21, 2) = 0.28900d-02
        htt(22, 1, 2) = 0.27000d-03
        htt(22, 2, 2) = 0.11300d-02
        htt(22, 3, 2) = -0.15900d-02
        htt(22, 4, 2) = -0.74000d-03
        htt(22, 5, 2) = -0.60000d-04
        htt(22, 6, 2) = -0.26200d-02
        htt(22, 7, 2) = 0.17100d-02
        htt(22, 8, 2) = -0.85000d-03
        htt(22, 9, 2) = -0.85000d-03
        htt(22,10, 2) = 0.40000d-04
        htt(22,11, 2) = -0.17300d-02
        htt(22,12, 2) = -0.85000d-03
        htt(22,13, 2) = 0.12800d-02
        htt(22,14, 2) = -0.15000d-03
        htt(22,15, 2) = 0.40000d-03
        htt(22,16, 2) = -0.10300d-02
        htt(22,17, 2) = -0.10700d-02
        htt(22,18, 2) = 0.56000d-03
        htt(22,19, 2) = 0.37000d-03
        htt(22,20, 2) = 0.20000d-02
        htt(22,21, 2) = 0.44000d-03
        htt(22,22, 2) = 0.31200d-02
        htt(23, 1, 2) = 0.19200d-02
        htt(23, 2, 2) = 0.79000d-03
        htt(23, 3, 2) = -0.23000d-03
        htt(23, 4, 2) = -0.13700d-02
        htt(23, 5, 2) = -0.21700d-02
        htt(23, 6, 2) = 0.14000d-03
        htt(23, 7, 2) = -0.12000d-03
        htt(23, 8, 2) = 0.21900d-02
        htt(23, 9, 2) = -0.74000d-03
        htt(23,10, 2) = -0.17100d-02
        htt(23,11, 2) = 0.45000d-03
        htt(23,12, 2) = -0.53000d-03
        htt(23,13, 2) = -0.10000d-03
        htt(23,14, 2) = -0.13700d-02
        htt(23,15, 2) = 0.87000d-03
        htt(23,16, 2) = -0.40000d-03
        htt(23,17, 2) = -0.15000d-03
        htt(23,18, 2) = -0.31500d-02
        htt(23,19, 2) = 0.11200d-02
        htt(23,20, 2) = -0.18800d-02
        htt(23,21, 2) = 0.12200d-02
        htt(23,22, 2) = -0.12000d-02
        htt(23,23, 2) = 0.42400d-02
        htt(24, 1, 2) = 0.30000d-04
        htt(24, 2, 2) = 0.42000d-03
        htt(24, 3, 2) = -0.84000d-03
        htt(24, 4, 2) = -0.45000d-03
        htt(24, 5, 2) = 0.79000d-03
        htt(24, 6, 2) = -0.18000d-02
        htt(24, 7, 2) = 0.16200d-02
        htt(24, 8, 2) = -0.97000d-03
        htt(24, 9, 2) = -0.39000d-03
        htt(24,10, 2) = -0.70000d-03
        htt(24,11, 2) = -0.79000d-03
        htt(24,12, 2) = -0.11000d-02
        htt(24,13, 2) = -0.80000d-04
        htt(24,14, 2) = -0.76000d-03
        htt(24,15, 2) = 0.22000d-03
        htt(24,16, 2) = -0.45000d-03
        htt(24,17, 2) = 0.92000d-03
        htt(24,18, 2) = -0.21000d-02
        htt(24,19, 2) = 0.16000d-02
        htt(24,20, 2) = -0.14300d-02
        htt(24,21, 2) = -0.12300d-02
        htt(24,22, 2) = 0.14800d-02
        htt(24,23, 2) = 0.18100d-02
        htt(24,24, 2) = 0.45300d-02
        ! Hessian elements for anchor point 2
        hbs( 1, 1, 2) = -0.31760d-01
        hbs( 1, 2, 2) = -0.29360d-01
        hbs( 1, 3, 2) = -0.31820d-01
        hbs( 1, 4, 2) = -0.95900d-02
        hbs( 1, 5, 2) = 0.52300d-02
        hbs( 1, 6, 2) = 0.39930d-01
        hbs( 1, 7, 2) = -0.45000d-03
        hbs( 1, 8, 2) = 0.40450d-01
        hbs( 1, 9, 2) = -0.18000d-03
        hbs( 1,10, 2) = -0.96900d-02
        hbs( 1,11, 2) = -0.42000d-03
        hbs( 1,12, 2) = 0.47400d-02
        hbs( 2, 1, 2) = 0.49930d-01
        hbs( 2, 2, 2) = -0.13930d-01
        hbs( 2, 3, 2) = 0.99600d-02
        hbs( 2, 4, 2) = 0.18000d-02
        hbs( 2, 5, 2) = -0.12200d-02
        hbs( 2, 6, 2) = -0.27810d-01
        hbs( 2, 7, 2) = 0.11800d-02
        hbs( 2, 8, 2) = -0.14650d-01
        hbs( 2, 9, 2) = 0.00000d+00
        hbs( 2,10, 2) = 0.41300d-02
        hbs( 2,11, 2) = -0.69000d-03
        hbs( 2,12, 2) = -0.89000d-03
        hbs( 3, 1, 2) = -0.18180d-01
        hbs( 3, 2, 2) = 0.43280d-01
        hbs( 3, 3, 2) = 0.21860d-01
        hbs( 3, 4, 2) = 0.77900d-02
        hbs( 3, 5, 2) = -0.40100d-02
        hbs( 3, 6, 2) = -0.12120d-01
        hbs( 3, 7, 2) = -0.73000d-03
        hbs( 3, 8, 2) = -0.25800d-01
        hbs( 3, 9, 2) = 0.18000d-03
        hbs( 3,10, 2) = 0.55600d-02
        hbs( 3,11, 2) = 0.11100d-02
        hbs( 3,12, 2) = -0.38500d-02
        hbs( 4, 1, 2) = -0.20820d-01
        hbs( 4, 2, 2) = -0.12100d-01
        hbs( 4, 3, 2) = 0.15180d-01
        hbs( 4, 4, 2) = -0.29840d-01
        hbs( 4, 5, 2) = -0.97300d-02
        hbs( 4, 6, 2) = -0.57000d-02
        hbs( 4, 7, 2) = 0.49500d-02
        hbs( 4, 8, 2) = 0.31150d-01
        hbs( 4, 9, 2) = -0.44000d-03
        hbs( 4,10, 2) = 0.40270d-01
        hbs( 4,11, 2) = -0.59000d-03
        hbs( 4,12, 2) = 0.29000d-03
        hbs( 5, 1, 2) = 0.26090d-01
        hbs( 5, 2, 2) = 0.58900d-02
        hbs( 5, 3, 2) = -0.14520d-01
        hbs( 5, 4, 2) = 0.25400d-02
        hbs( 5, 5, 2) = 0.69400d-02
        hbs( 5, 6, 2) = 0.30900d-02
        hbs( 5, 7, 2) = -0.21700d-02
        hbs( 5, 8, 2) = -0.15190d-01
        hbs( 5, 9, 2) = -0.39000d-03
        hbs( 5,10, 2) = -0.22140d-01
        hbs( 5,11, 2) = 0.28000d-03
        hbs( 5,12, 2) = 0.90000d-04
        hbs( 6, 1, 2) = -0.52700d-02
        hbs( 6, 2, 2) = 0.62100d-02
        hbs( 6, 3, 2) = -0.66000d-03
        hbs( 6, 4, 2) = 0.27300d-01
        hbs( 6, 5, 2) = 0.27900d-02
        hbs( 6, 6, 2) = 0.26100d-02
        hbs( 6, 7, 2) = -0.27800d-02
        hbs( 6, 8, 2) = -0.15960d-01
        hbs( 6, 9, 2) = 0.83000d-03
        hbs( 6,10, 2) = -0.18130d-01
        hbs( 6,11, 2) = 0.31000d-03
        hbs( 6,12, 2) = -0.38000d-03
        hbs( 7, 1, 2) = -0.90800d-02
        hbs( 7, 2, 2) = 0.35870d-01
        hbs( 7, 3, 2) = -0.10000d-04
        hbs( 7, 4, 2) = -0.22920d-01
        hbs( 7, 5, 2) = 0.61900d-02
        hbs( 7, 6, 2) = -0.28970d-01
        hbs( 7, 7, 2) = -0.87500d-02
        hbs( 7, 8, 2) = -0.92900d-02
        hbs( 7, 9, 2) = 0.49700d-02
        hbs( 7,10, 2) = 0.32780d-01
        hbs( 7,11, 2) = 0.60000d-04
        hbs( 7,12, 2) = -0.10700d-02
        hbs( 8, 1, 2) = 0.36900d-02
        hbs( 8, 2, 2) = -0.18590d-01
        hbs( 8, 3, 2) = 0.13800d-02
        hbs( 8, 4, 2) = 0.25350d-01
        hbs( 8, 5, 2) = -0.35000d-02
        hbs( 8, 6, 2) = -0.42000d-03
        hbs( 8, 7, 2) = 0.45000d-02
        hbs( 8, 8, 2) = 0.50400d-02
        hbs( 8, 9, 2) = -0.20100d-02
        hbs( 8,10, 2) = -0.14680d-01
        hbs( 8,11, 2) = -0.62000d-03
        hbs( 8,12, 2) = 0.55000d-03
        hbs( 9, 1, 2) = 0.53900d-02
        hbs( 9, 2, 2) = -0.17270d-01
        hbs( 9, 3, 2) = -0.13700d-02
        hbs( 9, 4, 2) = -0.24400d-02
        hbs( 9, 5, 2) = -0.26800d-02
        hbs( 9, 6, 2) = 0.29390d-01
        hbs( 9, 7, 2) = 0.42400d-02
        hbs( 9, 8, 2) = 0.42500d-02
        hbs( 9, 9, 2) = -0.29600d-02
        hbs( 9,10, 2) = -0.18100d-01
        hbs( 9,11, 2) = 0.55000d-03
        hbs( 9,12, 2) = 0.51000d-03
        hbs(10, 1, 2) = 0.36130d-01
        hbs(10, 2, 2) = 0.38360d-01
        hbs(10, 3, 2) = 0.36300d-02
        hbs(10, 4, 2) = -0.11890d-01
        hbs(10, 5, 2) = -0.10500d-02
        hbs(10, 6, 2) = -0.26840d-01
        hbs(10, 7, 2) = 0.48100d-02
        hbs(10, 8, 2) = -0.27100d-01
        hbs(10, 9, 2) = -0.88200d-02
        hbs(10,10, 2) = -0.12140d-01
        hbs(10,11, 2) = 0.47600d-02
        hbs(10,12, 2) = -0.77000d-03
        hbs(11, 1, 2) = -0.20190d-01
        hbs(11, 2, 2) = -0.17180d-01
        hbs(11, 3, 2) = -0.21300d-02
        hbs(11, 4, 2) = 0.57400d-02
        hbs(11, 5, 2) = 0.12900d-02
        hbs(11, 6, 2) = 0.27000d-01
        hbs(11, 7, 2) = -0.28400d-02
        hbs(11, 8, 2) = 0.47000d-03
        hbs(11, 9, 2) = 0.44300d-02
        hbs(11,10, 2) = 0.60700d-02
        hbs(11,11, 2) = -0.19800d-02
        hbs(11,12, 2) = -0.18000d-03
        hbs(12, 1, 2) = -0.15940d-01
        hbs(12, 2, 2) = -0.21180d-01
        hbs(12, 3, 2) = -0.14900d-02
        hbs(12, 4, 2) = 0.61500d-02
        hbs(12, 5, 2) = -0.24000d-03
        hbs(12, 6, 2) = -0.16000d-03
        hbs(12, 7, 2) = -0.19600d-02
        hbs(12, 8, 2) = 0.26630d-01
        hbs(12, 9, 2) = 0.43900d-02
        hbs(12,10, 2) = 0.60700d-02
        hbs(12,11, 2) = -0.27800d-02
        hbs(12,12, 2) = 0.94000d-03
        hbs(13, 1, 2) = 0.35720d-01
        hbs(13, 2, 2) = -0.96900d-02
        hbs(13, 3, 2) = 0.29800d-02
        hbs(13, 4, 2) = 0.32310d-01
        hbs(13, 5, 2) = -0.94000d-03
        hbs(13, 6, 2) = -0.10140d-01
        hbs(13, 7, 2) = 0.70000d-04
        hbs(13, 8, 2) = -0.29530d-01
        hbs(13, 9, 2) = 0.50200d-02
        hbs(13,10, 2) = -0.23420d-01
        hbs(13,11, 2) = -0.86600d-02
        hbs(13,12, 2) = 0.53400d-02
        hbs(14, 1, 2) = -0.17210d-01
        hbs(14, 2, 2) = 0.56100d-02
        hbs(14, 3, 2) = -0.23300d-02
        hbs(14, 4, 2) = -0.17680d-01
        hbs(14, 5, 2) = 0.48000d-03
        hbs(14, 6, 2) = 0.47700d-02
        hbs(14, 7, 2) = 0.56000d-03
        hbs(14, 8, 2) = 0.29630d-01
        hbs(14, 9, 2) = -0.29700d-02
        hbs(14,10, 2) = -0.24000d-02
        hbs(14,11, 2) = 0.43000d-02
        hbs(14,12, 2) = -0.22600d-02
        hbs(15, 1, 2) = -0.18510d-01
        hbs(15, 2, 2) = 0.40800d-02
        hbs(15, 3, 2) = -0.65000d-03
        hbs(15, 4, 2) = -0.14630d-01
        hbs(15, 5, 2) = 0.46000d-03
        hbs(15, 6, 2) = 0.53700d-02
        hbs(15, 7, 2) = -0.63000d-03
        hbs(15, 8, 2) = -0.10000d-03
        hbs(15, 9, 2) = -0.20500d-02
        hbs(15,10, 2) = 0.25820d-01
        hbs(15,11, 2) = 0.43700d-02
        hbs(15,12, 2) = -0.30900d-02
        hbs(16, 1, 2) = -0.10200d-01
        hbs(16, 2, 2) = -0.23080d-01
        hbs(16, 3, 2) = 0.10040d-01
        hbs(16, 4, 2) = 0.41930d-01
        hbs(16, 5, 2) = 0.31000d-03
        hbs(16, 6, 2) = 0.31710d-01
        hbs(16, 7, 2) = -0.63000d-03
        hbs(16, 8, 2) = -0.56800d-02
        hbs(16, 9, 2) = -0.55000d-03
        hbs(16,10, 2) = -0.27810d-01
        hbs(16,11, 2) = 0.48500d-02
        hbs(16,12, 2) = -0.85400d-02
        hbs(17, 1, 2) = 0.41100d-02
        hbs(17, 2, 2) = 0.25000d-01
        hbs(17, 3, 2) = -0.74400d-02
        hbs(17, 4, 2) = -0.22980d-01
        hbs(17, 5, 2) = 0.60000d-03
        hbs(17, 6, 2) = -0.14860d-01
        hbs(17, 7, 2) = 0.32000d-03
        hbs(17, 8, 2) = 0.33600d-02
        hbs(17, 9, 2) = -0.32000d-03
        hbs(17,10, 2) = 0.20400d-02
        hbs(17,11, 2) = -0.21000d-02
        hbs(17,12, 2) = 0.68800d-02
        hbs(18, 1, 2) = 0.61000d-02
        hbs(18, 2, 2) = -0.19200d-02
        hbs(18, 3, 2) = -0.26000d-02
        hbs(18, 4, 2) = -0.18950d-01
        hbs(18, 5, 2) = -0.92000d-03
        hbs(18, 6, 2) = -0.16850d-01
        hbs(18, 7, 2) = 0.31000d-03
        hbs(18, 8, 2) = 0.23200d-02
        hbs(18, 9, 2) = 0.87000d-03
        hbs(18,10, 2) = 0.25770d-01
        hbs(18,11, 2) = -0.27500d-02
        hbs(18,12, 2) = 0.16600d-02
        ! Hessian elements for anchor point 3
        hss( 1, 1, 3) = 0.34759d+00
        hss( 2, 1, 3) = 0.17810d-01
        hss( 2, 2, 3) = 0.35623d+00
        hss( 3, 1, 3) = 0.29010d-01
        hss( 3, 2, 3) = 0.32710d-01
        hss( 3, 3, 3) = 0.46124d+00
        hss( 4, 1, 3) = -0.19600d-02
        hss( 4, 2, 3) = 0.14870d-01
        hss( 4, 3, 3) = -0.89100d-02
        hss( 4, 4, 3) = 0.38581d+00
        hss( 5, 1, 3) = 0.38800d-02
        hss( 5, 2, 3) = -0.77000d-03
        hss( 5, 3, 3) = 0.71700d-02
        hss( 5, 4, 3) = 0.36500d-02
        hss( 5, 5, 3) = 0.35807d+00
        hss( 6, 1, 3) = 0.33890d-01
        hss( 6, 2, 3) = 0.44570d-01
        hss( 6, 3, 3) = -0.43100d-02
        hss( 6, 4, 3) = -0.36100d-02
        hss( 6, 5, 3) = 0.55000d-03
        hss( 6, 6, 3) = 0.35931d+00
        hss( 7, 1, 3) = -0.54000d-03
        hss( 7, 2, 3) = -0.63000d-03
        hss( 7, 3, 3) = 0.57000d-03
        hss( 7, 4, 3) = 0.39700d-02
        hss( 7, 5, 3) = 0.39000d-03
        hss( 7, 6, 3) = 0.39800d-02
        hss( 7, 7, 3) = 0.37091d+00
        hss( 8, 1, 3) = 0.42030d-01
        hss( 8, 2, 3) = 0.35060d-01
        hss( 8, 3, 3) = -0.52800d-02
        hss( 8, 4, 3) = 0.24280d-01
        hss( 8, 5, 3) = -0.11100d-02
        hss( 8, 6, 3) = 0.17880d-01
        hss( 8, 7, 3) = -0.10400d-02
        hss( 8, 8, 3) = 0.35888d+00
        hss( 9, 1, 3) = -0.53000d-03
        hss( 9, 2, 3) = -0.57000d-03
        hss( 9, 3, 3) = 0.26000d-03
        hss( 9, 4, 3) = -0.39000d-03
        hss( 9, 5, 3) = 0.19000d-03
        hss( 9, 6, 3) = 0.40100d-02
        hss( 9, 7, 3) = 0.49000d-03
        hss( 9, 8, 3) = 0.39600d-02
        hss( 9, 9, 3) = 0.36356d+00
        hss(10, 1, 3) = 0.19520d-01
        hss(10, 2, 3) = -0.69300d-02
        hss(10, 3, 3) = -0.48800d-02
        hss(10, 4, 3) = 0.64010d-01
        hss(10, 5, 3) = -0.43000d-03
        hss(10, 6, 3) = 0.21230d-01
        hss(10, 7, 3) = -0.10900d-02
        hss(10, 8, 3) = -0.77000d-03
        hss(10, 9, 3) = -0.19000d-03
        hss(10,10, 3) = 0.38116d+00
        hss(11, 1, 3) = -0.69000d-03
        hss(11, 2, 3) = -0.50000d-03
        hss(11, 3, 3) = 0.52000d-03
        hss(11, 4, 3) = -0.11600d-02
        hss(11, 5, 3) = 0.80000d-04
        hss(11, 6, 3) = -0.10000d-02
        hss(11, 7, 3) = 0.22000d-03
        hss(11, 8, 3) = 0.39700d-02
        hss(11, 9, 3) = 0.48000d-03
        hss(11,10, 3) = 0.39800d-02
        hss(11,11, 3) = 0.37152d+00
        hss(12, 1, 3) = -0.16500d-02
        hss(12, 2, 3) = 0.39400d-02
        hss(12, 3, 3) = -0.19000d-03
        hss(12, 4, 3) = -0.10000d-03
        hss(12, 5, 3) = 0.35000d-03
        hss(12, 6, 3) = -0.10400d-02
        hss(12, 7, 3) = 0.30000d-04
        hss(12, 8, 3) = 0.57000d-03
        hss(12, 9, 3) = 0.12000d-03
        hss(12,10, 3) = 0.23600d-02
        hss(12,11, 3) = 0.40000d-03
        hss(12,12, 3) = 0.37254d+00
        ! Hessian elements for anchor point 3
        hbb( 1, 1, 3) = 0.82450d-01
        hbb( 2, 1, 3) = -0.46090d-01
        hbb( 2, 2, 3) = 0.14487d+00
        hbb( 3, 1, 3) = -0.36360d-01
        hbb( 3, 2, 3) = -0.98780d-01
        hbb( 3, 3, 3) = 0.13515d+00
        hbb( 4, 1, 3) = -0.36000d-01
        hbb( 4, 2, 3) = 0.23310d-01
        hbb( 4, 3, 3) = 0.12690d-01
        hbb( 4, 4, 3) = 0.66170d-01
        hbb( 5, 1, 3) = 0.20240d-01
        hbb( 5, 2, 3) = -0.12480d-01
        hbb( 5, 3, 3) = -0.77700d-02
        hbb( 5, 4, 3) = -0.33030d-01
        hbb( 5, 5, 3) = 0.71950d-01
        hbb( 6, 1, 3) = 0.15760d-01
        hbb( 6, 2, 3) = -0.10840d-01
        hbb( 6, 3, 3) = -0.49200d-02
        hbb( 6, 4, 3) = -0.33140d-01
        hbb( 6, 5, 3) = -0.38920d-01
        hbb( 6, 6, 3) = 0.72060d-01
        hbb( 7, 1, 3) = -0.71500d-02
        hbb( 7, 2, 3) = -0.39000d-02
        hbb( 7, 3, 3) = 0.11050d-01
        hbb( 7, 4, 3) = -0.25740d-01
        hbb( 7, 5, 3) = 0.93100d-02
        hbb( 7, 6, 3) = 0.16430d-01
        hbb( 7, 7, 3) = 0.64420d-01
        hbb( 8, 1, 3) = -0.22800d-02
        hbb( 8, 2, 3) = 0.74300d-02
        hbb( 8, 3, 3) = -0.51500d-02
        hbb( 8, 4, 3) = 0.15690d-01
        hbb( 8, 5, 3) = -0.49500d-02
        hbb( 8, 6, 3) = -0.10740d-01
        hbb( 8, 7, 3) = -0.31490d-01
        hbb( 8, 8, 3) = 0.72950d-01
        hbb( 9, 1, 3) = 0.94300d-02
        hbb( 9, 2, 3) = -0.35300d-02
        hbb( 9, 3, 3) = -0.59100d-02
        hbb( 9, 4, 3) = 0.10050d-01
        hbb( 9, 5, 3) = -0.43500d-02
        hbb( 9, 6, 3) = -0.56900d-02
        hbb( 9, 7, 3) = -0.32930d-01
        hbb( 9, 8, 3) = -0.41460d-01
        hbb( 9, 9, 3) = 0.74390d-01
        hbb(10, 1, 3) = 0.45500d-02
        hbb(10, 2, 3) = -0.20200d-02
        hbb(10, 3, 3) = -0.25300d-02
        hbb(10, 4, 3) = -0.64100d-02
        hbb(10, 5, 3) = 0.91900d-02
        hbb(10, 6, 3) = -0.27700d-02
        hbb(10, 7, 3) = -0.33410d-01
        hbb(10, 8, 3) = 0.13410d-01
        hbb(10, 9, 3) = 0.20000d-01
        hbb(10,10, 3) = 0.75440d-01
        hbb(11, 1, 3) = -0.22100d-02
        hbb(11, 2, 3) = 0.12300d-02
        hbb(11, 3, 3) = 0.98000d-03
        hbb(11, 4, 3) = -0.23100d-02
        hbb(11, 5, 3) = -0.32200d-02
        hbb(11, 6, 3) = 0.55300d-02
        hbb(11, 7, 3) = 0.20020d-01
        hbb(11, 8, 3) = -0.74100d-02
        hbb(11, 9, 3) = -0.12610d-01
        hbb(11,10, 3) = -0.37790d-01
        hbb(11,11, 3) = 0.78360d-01
        hbb(12, 1, 3) = -0.23400d-02
        hbb(12, 2, 3) = 0.79000d-03
        hbb(12, 3, 3) = 0.15500d-02
        hbb(12, 4, 3) = 0.87200d-02
        hbb(12, 5, 3) = -0.59700d-02
        hbb(12, 6, 3) = -0.27500d-02
        hbb(12, 7, 3) = 0.13390d-01
        hbb(12, 8, 3) = -0.60000d-02
        hbb(12, 9, 3) = -0.73900d-02
        hbb(12,10, 3) = -0.37640d-01
        hbb(12,11, 3) = -0.40570d-01
        hbb(12,12, 3) = 0.78220d-01
        hbb(13, 1, 3) = -0.67600d-02
        hbb(13, 2, 3) = 0.10710d-01
        hbb(13, 3, 3) = -0.39500d-02
        hbb(13, 4, 3) = -0.11500d-02
        hbb(13, 5, 3) = 0.12900d-02
        hbb(13, 6, 3) = -0.14000d-03
        hbb(13, 7, 3) = 0.29100d-02
        hbb(13, 8, 3) = 0.41300d-02
        hbb(13, 9, 3) = -0.70300d-02
        hbb(13,10, 3) = -0.33760d-01
        hbb(13,11, 3) = 0.13510d-01
        hbb(13,12, 3) = 0.20240d-01
        hbb(13,13, 3) = 0.65300d-01
        hbb(14, 1, 3) = 0.91200d-02
        hbb(14, 2, 3) = -0.60500d-02
        hbb(14, 3, 3) = -0.30700d-02
        hbb(14, 4, 3) = 0.66000d-03
        hbb(14, 5, 3) = -0.85000d-03
        hbb(14, 6, 3) = 0.19000d-03
        hbb(14, 7, 3) = -0.71100d-02
        hbb(14, 8, 3) = -0.51000d-03
        hbb(14, 9, 3) = 0.76200d-02
        hbb(14,10, 3) = 0.20030d-01
        hbb(14,11, 3) = -0.73800d-02
        hbb(14,12, 3) = -0.12650d-01
        hbb(14,13, 3) = -0.33180d-01
        hbb(14,14, 3) = 0.74350d-01
        hbb(15, 1, 3) = -0.23600d-02
        hbb(15, 2, 3) = -0.46500d-02
        hbb(15, 3, 3) = 0.70200d-02
        hbb(15, 4, 3) = 0.49000d-03
        hbb(15, 5, 3) = -0.44000d-03
        hbb(15, 6, 3) = -0.50000d-04
        hbb(15, 7, 3) = 0.42000d-02
        hbb(15, 8, 3) = -0.36200d-02
        hbb(15, 9, 3) = -0.59000d-03
        hbb(15,10, 3) = 0.13730d-01
        hbb(15,11, 3) = -0.61300d-02
        hbb(15,12, 3) = -0.76000d-02
        hbb(15,13, 3) = -0.32120d-01
        hbb(15,14, 3) = -0.41170d-01
        hbb(15,15, 3) = 0.73290d-01
        hbb(16, 1, 3) = -0.37090d-01
        hbb(16, 2, 3) = 0.17990d-01
        hbb(16, 3, 3) = 0.19100d-01
        hbb(16, 4, 3) = 0.31400d-02
        hbb(16, 5, 3) = -0.70000d-02
        hbb(16, 6, 3) = 0.38600d-02
        hbb(16, 7, 3) = -0.10300d-02
        hbb(16, 8, 3) = 0.54000d-03
        hbb(16, 9, 3) = 0.48000d-03
        hbb(16,10, 3) = -0.64000d-02
        hbb(16,11, 3) = 0.87800d-02
        hbb(16,12, 3) = -0.23700d-02
        hbb(16,13, 3) = -0.26530d-01
        hbb(16,14, 3) = 0.10480d-01
        hbb(16,15, 3) = 0.16060d-01
        hbb(16,16, 3) = 0.67910d-01
        hbb(17, 1, 3) = 0.22610d-01
        hbb(17, 2, 3) = -0.90200d-02
        hbb(17, 3, 3) = -0.13590d-01
        hbb(17, 4, 3) = -0.76000d-02
        hbb(17, 5, 3) = 0.77900d-02
        hbb(17, 6, 3) = -0.19000d-03
        hbb(17, 7, 3) = 0.87000d-03
        hbb(17, 8, 3) = -0.25000d-03
        hbb(17, 9, 3) = -0.62000d-03
        hbb(17,10, 3) = 0.91600d-02
        hbb(17,11, 3) = -0.60700d-02
        hbb(17,12, 3) = -0.31000d-02
        hbb(17,13, 3) = 0.96800d-02
        hbb(17,14, 3) = -0.46000d-02
        hbb(17,15, 3) = -0.50800d-02
        hbb(17,16, 3) = -0.34720d-01
        hbb(17,17, 3) = 0.72390d-01
        hbb(18, 1, 3) = 0.14470d-01
        hbb(18, 2, 3) = -0.89700d-02
        hbb(18, 3, 3) = -0.55000d-02
        hbb(18, 4, 3) = 0.44600d-02
        hbb(18, 5, 3) = -0.79000d-03
        hbb(18, 6, 3) = -0.36700d-02
        hbb(18, 7, 3) = 0.16000d-03
        hbb(18, 8, 3) = -0.29000d-03
        hbb(18, 9, 3) = 0.13000d-03
        hbb(18,10, 3) = -0.27600d-02
        hbb(18,11, 3) = -0.27100d-02
        hbb(18,12, 3) = 0.54700d-02
        hbb(18,13, 3) = 0.16860d-01
        hbb(18,14, 3) = -0.58800d-02
        hbb(18,15, 3) = -0.10980d-01
        hbb(18,16, 3) = -0.33190d-01
        hbb(18,17, 3) = -0.37670d-01
        hbb(18,18, 3) = 0.70860d-01
        ! Hessian elements for anchor point 3
        htt( 1, 1, 3) = 0.60400d-02
        htt( 2, 1, 3) = 0.38400d-02
        htt( 2, 2, 3) = 0.50900d-02
        htt( 3, 1, 3) = -0.18500d-02
        htt( 3, 2, 3) = -0.33100d-02
        htt( 3, 3, 3) = 0.75100d-02
        htt( 4, 1, 3) = -0.40500d-02
        htt( 4, 2, 3) = -0.20700d-02
        htt( 4, 3, 3) = 0.60500d-02
        htt( 4, 4, 3) = 0.80200d-02
        htt( 5, 1, 3) = -0.59200d-02
        htt( 5, 2, 3) = -0.42600d-02
        htt( 5, 3, 3) = 0.21200d-02
        htt( 5, 4, 3) = 0.37800d-02
        htt( 5, 5, 3) = 0.65600d-02
        htt( 6, 1, 3) = -0.41200d-02
        htt( 6, 2, 3) = -0.39000d-02
        htt( 6, 3, 3) = 0.26900d-02
        htt( 6, 4, 3) = 0.29100d-02
        htt( 6, 5, 3) = 0.37400d-02
        htt( 6, 6, 3) = 0.55900d-02
        htt( 7, 1, 3) = 0.15700d-02
        htt( 7, 2, 3) = 0.25400d-02
        htt( 7, 3, 3) = -0.67700d-02
        htt( 7, 4, 3) = -0.58000d-02
        htt( 7, 5, 3) = -0.10700d-02
        htt( 7, 6, 3) = -0.27200d-02
        htt( 7, 7, 3) = 0.68600d-02
        htt( 8, 1, 3) = 0.33800d-02
        htt( 8, 2, 3) = 0.29000d-02
        htt( 8, 3, 3) = -0.62000d-02
        htt( 8, 4, 3) = -0.66800d-02
        htt( 8, 5, 3) = -0.38900d-02
        htt( 8, 6, 3) = -0.88000d-03
        htt( 8, 7, 3) = 0.52100d-02
        htt( 8, 8, 3) = 0.82200d-02
        htt( 9, 1, 3) = -0.23900d-02
        htt( 9, 2, 3) = -0.68000d-03
        htt( 9, 3, 3) = 0.45000d-03
        htt( 9, 4, 3) = 0.21600d-02
        htt( 9, 5, 3) = 0.17400d-02
        htt( 9, 6, 3) = 0.14100d-02
        htt( 9, 7, 3) = -0.96000d-03
        htt( 9, 8, 3) = -0.12900d-02
        htt( 9, 9, 3) = 0.24600d-02
        htt(10, 1, 3) = -0.15600d-02
        htt(10, 2, 3) = 0.23000d-03
        htt(10, 3, 3) = 0.34000d-03
        htt(10, 4, 3) = 0.21300d-02
        htt(10, 5, 3) = 0.17600d-02
        htt(10, 6, 3) = 0.79000d-03
        htt(10, 7, 3) = -0.40000d-04
        htt(10, 8, 3) = -0.10100d-02
        htt(10, 9, 3) = 0.57000d-03
        htt(10,10, 3) = 0.39100d-02
        htt(11, 1, 3) = -0.10000d-03
        htt(11, 2, 3) = -0.19700d-02
        htt(11, 3, 3) = 0.19800d-02
        htt(11, 4, 3) = 0.11000d-03
        htt(11, 5, 3) = 0.10000d-04
        htt(11, 6, 3) = 0.11900d-02
        htt(11, 7, 3) = -0.19700d-02
        htt(11, 8, 3) = -0.79000d-03
        htt(11, 9, 3) = 0.68000d-03
        htt(11,10, 3) = -0.12900d-02
        htt(11,11, 3) = 0.26200d-02
        htt(12, 1, 3) = 0.72000d-03
        htt(12, 2, 3) = -0.10500d-02
        htt(12, 3, 3) = 0.18600d-02
        htt(12, 4, 3) = 0.90000d-04
        htt(12, 5, 3) = 0.30000d-04
        htt(12, 6, 3) = 0.57000d-03
        htt(12, 7, 3) = -0.10500d-02
        htt(12, 8, 3) = -0.51000d-03
        htt(12, 9, 3) = -0.12000d-02
        htt(12,10, 3) = 0.20500d-02
        htt(12,11, 3) = 0.64000d-03
        htt(12,12, 3) = 0.39000d-02
        htt(13, 1, 3) = -0.10900d-02
        htt(13, 2, 3) = -0.19400d-02
        htt(13, 3, 3) = 0.59000d-03
        htt(13, 4, 3) = -0.26000d-03
        htt(13, 5, 3) = 0.15200d-02
        htt(13, 6, 3) = 0.15100d-02
        htt(13, 7, 3) = -0.70000d-04
        htt(13, 8, 3) = -0.90000d-04
        htt(13, 9, 3) = -0.21100d-02
        htt(13,10, 3) = 0.12000d-03
        htt(13,11, 3) = -0.12300d-02
        htt(13,12, 3) = 0.10000d-02
        htt(13,13, 3) = 0.50300d-02
        htt(14, 1, 3) = 0.70000d-03
        htt(14, 2, 3) = 0.53000d-03
        htt(14, 3, 3) = -0.31000d-03
        htt(14, 4, 3) = -0.47000d-03
        htt(14, 5, 3) = -0.22000d-03
        htt(14, 6, 3) = -0.80000d-03
        htt(14, 7, 3) = 0.74000d-03
        htt(14, 8, 3) = 0.16000d-03
        htt(14, 9, 3) = -0.19000d-03
        htt(14,10, 3) = 0.13500d-02
        htt(14,11, 3) = -0.10000d-04
        htt(14,12, 3) = 0.15300d-02
        htt(14,13, 3) = -0.81000d-03
        htt(14,14, 3) = 0.39200d-02
        htt(15, 1, 3) = -0.19100d-02
        htt(15, 2, 3) = -0.28400d-02
        htt(15, 3, 3) = 0.70000d-03
        htt(15, 4, 3) = -0.23000d-03
        htt(15, 5, 3) = 0.15000d-02
        htt(15, 6, 3) = 0.21100d-02
        htt(15, 7, 3) = -0.98000d-03
        htt(15, 8, 3) = -0.37000d-03
        htt(15, 9, 3) = -0.24000d-03
        htt(15,10, 3) = -0.31800d-02
        htt(15,11, 3) = 0.73000d-03
        htt(15,12, 3) = -0.22200d-02
        htt(15,13, 3) = 0.28200d-02
        htt(15,14, 3) = -0.23300d-02
        htt(15,15, 3) = 0.57300d-02
        htt(16, 1, 3) = -0.12000d-03
        htt(16, 2, 3) = -0.37000d-03
        htt(16, 3, 3) = -0.19000d-03
        htt(16, 4, 3) = -0.45000d-03
        htt(16, 5, 3) = -0.25000d-03
        htt(16, 6, 3) = -0.19000d-03
        htt(16, 7, 3) = -0.17000d-03
        htt(16, 8, 3) = -0.12000d-03
        htt(16, 9, 3) = 0.16800d-02
        htt(16,10, 3) = -0.19500d-02
        htt(16,11, 3) = 0.19400d-02
        htt(16,12, 3) = -0.16900d-02
        htt(16,13, 3) = -0.30200d-02
        htt(16,14, 3) = 0.24000d-02
        htt(16,15, 3) = 0.58000d-03
        htt(16,16, 3) = 0.60000d-02
        htt(17, 1, 3) = 0.12100d-02
        htt(17, 2, 3) = 0.15200d-02
        htt(17, 3, 3) = -0.31000d-03
        htt(17, 4, 3) = 0.00000d+00
        htt(17, 5, 3) = -0.86000d-03
        htt(17, 6, 3) = -0.18800d-02
        htt(17, 7, 3) = 0.58000d-03
        htt(17, 8, 3) = -0.44000d-03
        htt(17, 9, 3) = 0.14500d-02
        htt(17,10, 3) = 0.90000d-04
        htt(17,11, 3) = 0.11300d-02
        htt(17,12, 3) = -0.23000d-03
        htt(17,13, 3) = -0.46000d-02
        htt(17,14, 3) = 0.13000d-02
        htt(17,15, 3) = -0.32400d-02
        htt(17,16, 3) = 0.26500d-02
        htt(17,17, 3) = 0.49600d-02
        htt(18, 1, 3) = 0.14500d-02
        htt(18, 2, 3) = 0.22200d-02
        htt(18, 3, 3) = -0.10600d-02
        htt(18, 4, 3) = -0.29000d-03
        htt(18, 5, 3) = -0.17100d-02
        htt(18, 6, 3) = -0.27100d-02
        htt(18, 7, 3) = 0.67000d-03
        htt(18, 8, 3) = -0.32000d-03
        htt(18, 9, 3) = 0.10000d-02
        htt(18,10, 3) = 0.82000d-03
        htt(18,11, 3) = 0.19000d-03
        htt(18,12, 3) = 0.10000d-04
        htt(18,13, 3) = -0.32400d-02
        htt(18,14, 3) = 0.19000d-02
        htt(18,15, 3) = -0.30600d-02
        htt(18,16, 3) = 0.20800d-02
        htt(18,17, 3) = 0.29700d-02
        htt(18,18, 3) = 0.56300d-02
        htt(19, 1, 3) = -0.59000d-03
        htt(19, 2, 3) = -0.96000d-03
        htt(19, 3, 3) = 0.59000d-03
        htt(19, 4, 3) = 0.21000d-03
        htt(19, 5, 3) = 0.89000d-03
        htt(19, 6, 3) = 0.43000d-03
        htt(19, 7, 3) = -0.23000d-03
        htt(19, 8, 3) = -0.69000d-03
        htt(19, 9, 3) = -0.48000d-03
        htt(19,10, 3) = -0.11500d-02
        htt(19,11, 3) = -0.90000d-04
        htt(19,12, 3) = -0.76000d-03
        htt(19,13, 3) = 0.12600d-02
        htt(19,14, 3) = -0.34500d-02
        htt(19,15, 3) = 0.19200d-02
        htt(19,16, 3) = -0.27800d-02
        htt(19,17, 3) = -0.95000d-03
        htt(19,18, 3) = -0.21800d-02
        htt(19,19, 3) = 0.37700d-02
        htt(20, 1, 3) = -0.35000d-03
        htt(20, 2, 3) = -0.26000d-03
        htt(20, 3, 3) = -0.16000d-03
        htt(20, 4, 3) = -0.80000d-04
        htt(20, 5, 3) = 0.40000d-04
        htt(20, 6, 3) = -0.40000d-03
        htt(20, 7, 3) = -0.14000d-03
        htt(20, 8, 3) = -0.57000d-03
        htt(20, 9, 3) = -0.93000d-03
        htt(20,10, 3) = -0.42000d-03
        htt(20,11, 3) = -0.10200d-02
        htt(20,12, 3) = -0.51000d-03
        htt(20,13, 3) = 0.26100d-02
        htt(20,14, 3) = -0.28400d-02
        htt(20,15, 3) = 0.21100d-02
        htt(20,16, 3) = -0.33500d-02
        htt(20,17, 3) = -0.29400d-02
        htt(20,18, 3) = 0.47000d-03
        htt(20,19, 3) = 0.25400d-02
        htt(20,20, 3) = 0.59400d-02
        htt(21, 1, 3) = 0.21700d-02
        htt(21, 2, 3) = 0.15000d-02
        htt(21, 3, 3) = -0.99000d-03
        htt(21, 4, 3) = -0.16600d-02
        htt(21, 5, 3) = -0.30200d-02
        htt(21, 6, 3) = -0.69000d-03
        htt(21, 7, 3) = -0.20000d-04
        htt(21, 8, 3) = 0.23100d-02
        htt(21, 9, 3) = -0.11900d-02
        htt(21,10, 3) = -0.98000d-03
        htt(21,11, 3) = -0.50000d-03
        htt(21,12, 3) = -0.28000d-03
        htt(21,13, 3) = 0.12600d-02
        htt(21,14, 3) = -0.76000d-03
        htt(21,15, 3) = 0.10500d-02
        htt(21,16, 3) = -0.97000d-03
        htt(21,17, 3) = -0.21400d-02
        htt(21,18, 3) = -0.48000d-03
        htt(21,19, 3) = -0.11000d-03
        htt(21,20, 3) = 0.15500d-02
        htt(21,21, 3) = 0.28900d-02
        htt(22, 1, 3) = 0.27000d-03
        htt(22, 2, 3) = 0.11300d-02
        htt(22, 3, 3) = -0.15900d-02
        htt(22, 4, 3) = -0.74000d-03
        htt(22, 5, 3) = -0.60000d-04
        htt(22, 6, 3) = -0.26200d-02
        htt(22, 7, 3) = 0.17100d-02
        htt(22, 8, 3) = -0.85000d-03
        htt(22, 9, 3) = -0.85000d-03
        htt(22,10, 3) = 0.40000d-04
        htt(22,11, 3) = -0.17300d-02
        htt(22,12, 3) = -0.85000d-03
        htt(22,13, 3) = 0.12800d-02
        htt(22,14, 3) = -0.15000d-03
        htt(22,15, 3) = 0.40000d-03
        htt(22,16, 3) = -0.10300d-02
        htt(22,17, 3) = -0.10700d-02
        htt(22,18, 3) = 0.56000d-03
        htt(22,19, 3) = 0.37000d-03
        htt(22,20, 3) = 0.20000d-02
        htt(22,21, 3) = 0.44000d-03
        htt(22,22, 3) = 0.31200d-02
        htt(23, 1, 3) = 0.19200d-02
        htt(23, 2, 3) = 0.79000d-03
        htt(23, 3, 3) = -0.23000d-03
        htt(23, 4, 3) = -0.13700d-02
        htt(23, 5, 3) = -0.21700d-02
        htt(23, 6, 3) = 0.14000d-03
        htt(23, 7, 3) = -0.12000d-03
        htt(23, 8, 3) = 0.21900d-02
        htt(23, 9, 3) = -0.74000d-03
        htt(23,10, 3) = -0.17100d-02
        htt(23,11, 3) = 0.45000d-03
        htt(23,12, 3) = -0.53000d-03
        htt(23,13, 3) = -0.10000d-03
        htt(23,14, 3) = -0.13700d-02
        htt(23,15, 3) = 0.87000d-03
        htt(23,16, 3) = -0.40000d-03
        htt(23,17, 3) = -0.15000d-03
        htt(23,18, 3) = -0.31500d-02
        htt(23,19, 3) = 0.11200d-02
        htt(23,20, 3) = -0.18800d-02
        htt(23,21, 3) = 0.12200d-02
        htt(23,22, 3) = -0.12000d-02
        htt(23,23, 3) = 0.42400d-02
        htt(24, 1, 3) = 0.30000d-04
        htt(24, 2, 3) = 0.42000d-03
        htt(24, 3, 3) = -0.84000d-03
        htt(24, 4, 3) = -0.45000d-03
        htt(24, 5, 3) = 0.79000d-03
        htt(24, 6, 3) = -0.18000d-02
        htt(24, 7, 3) = 0.16200d-02
        htt(24, 8, 3) = -0.97000d-03
        htt(24, 9, 3) = -0.39000d-03
        htt(24,10, 3) = -0.70000d-03
        htt(24,11, 3) = -0.79000d-03
        htt(24,12, 3) = -0.11000d-02
        htt(24,13, 3) = -0.80000d-04
        htt(24,14, 3) = -0.76000d-03
        htt(24,15, 3) = 0.22000d-03
        htt(24,16, 3) = -0.45000d-03
        htt(24,17, 3) = 0.92000d-03
        htt(24,18, 3) = -0.21000d-02
        htt(24,19, 3) = 0.16000d-02
        htt(24,20, 3) = -0.14300d-02
        htt(24,21, 3) = -0.12300d-02
        htt(24,22, 3) = 0.14800d-02
        htt(24,23, 3) = 0.18100d-02
        htt(24,24, 3) = 0.45300d-02
        ! Hessian elements for anchor point 3
        hbs( 1, 1, 3) = -0.31760d-01
        hbs( 1, 2, 3) = -0.29360d-01
        hbs( 1, 3, 3) = -0.31820d-01
        hbs( 1, 4, 3) = -0.95900d-02
        hbs( 1, 5, 3) = 0.52300d-02
        hbs( 1, 6, 3) = 0.39930d-01
        hbs( 1, 7, 3) = -0.45000d-03
        hbs( 1, 8, 3) = 0.40450d-01
        hbs( 1, 9, 3) = -0.18000d-03
        hbs( 1,10, 3) = -0.96900d-02
        hbs( 1,11, 3) = -0.42000d-03
        hbs( 1,12, 3) = 0.47400d-02
        hbs( 2, 1, 3) = 0.49930d-01
        hbs( 2, 2, 3) = -0.13930d-01
        hbs( 2, 3, 3) = 0.99600d-02
        hbs( 2, 4, 3) = 0.18000d-02
        hbs( 2, 5, 3) = -0.12200d-02
        hbs( 2, 6, 3) = -0.27810d-01
        hbs( 2, 7, 3) = 0.11800d-02
        hbs( 2, 8, 3) = -0.14650d-01
        hbs( 2, 9, 3) = 0.00000d+00
        hbs( 2,10, 3) = 0.41300d-02
        hbs( 2,11, 3) = -0.69000d-03
        hbs( 2,12, 3) = -0.89000d-03
        hbs( 3, 1, 3) = -0.18180d-01
        hbs( 3, 2, 3) = 0.43280d-01
        hbs( 3, 3, 3) = 0.21860d-01
        hbs( 3, 4, 3) = 0.77900d-02
        hbs( 3, 5, 3) = -0.40100d-02
        hbs( 3, 6, 3) = -0.12120d-01
        hbs( 3, 7, 3) = -0.73000d-03
        hbs( 3, 8, 3) = -0.25800d-01
        hbs( 3, 9, 3) = 0.18000d-03
        hbs( 3,10, 3) = 0.55600d-02
        hbs( 3,11, 3) = 0.11100d-02
        hbs( 3,12, 3) = -0.38500d-02
        hbs( 4, 1, 3) = -0.20820d-01
        hbs( 4, 2, 3) = -0.12100d-01
        hbs( 4, 3, 3) = 0.15180d-01
        hbs( 4, 4, 3) = -0.29840d-01
        hbs( 4, 5, 3) = -0.97300d-02
        hbs( 4, 6, 3) = -0.57000d-02
        hbs( 4, 7, 3) = 0.49500d-02
        hbs( 4, 8, 3) = 0.31150d-01
        hbs( 4, 9, 3) = -0.44000d-03
        hbs( 4,10, 3) = 0.40270d-01
        hbs( 4,11, 3) = -0.59000d-03
        hbs( 4,12, 3) = 0.29000d-03
        hbs( 5, 1, 3) = 0.26090d-01
        hbs( 5, 2, 3) = 0.58900d-02
        hbs( 5, 3, 3) = -0.14520d-01
        hbs( 5, 4, 3) = 0.25400d-02
        hbs( 5, 5, 3) = 0.69400d-02
        hbs( 5, 6, 3) = 0.30900d-02
        hbs( 5, 7, 3) = -0.21700d-02
        hbs( 5, 8, 3) = -0.15190d-01
        hbs( 5, 9, 3) = -0.39000d-03
        hbs( 5,10, 3) = -0.22140d-01
        hbs( 5,11, 3) = 0.28000d-03
        hbs( 5,12, 3) = 0.90000d-04
        hbs( 6, 1, 3) = -0.52700d-02
        hbs( 6, 2, 3) = 0.62100d-02
        hbs( 6, 3, 3) = -0.66000d-03
        hbs( 6, 4, 3) = 0.27300d-01
        hbs( 6, 5, 3) = 0.27900d-02
        hbs( 6, 6, 3) = 0.26100d-02
        hbs( 6, 7, 3) = -0.27800d-02
        hbs( 6, 8, 3) = -0.15960d-01
        hbs( 6, 9, 3) = 0.83000d-03
        hbs( 6,10, 3) = -0.18130d-01
        hbs( 6,11, 3) = 0.31000d-03
        hbs( 6,12, 3) = -0.38000d-03
        hbs( 7, 1, 3) = -0.90800d-02
        hbs( 7, 2, 3) = 0.35870d-01
        hbs( 7, 3, 3) = -0.10000d-04
        hbs( 7, 4, 3) = -0.22920d-01
        hbs( 7, 5, 3) = 0.61900d-02
        hbs( 7, 6, 3) = -0.28970d-01
        hbs( 7, 7, 3) = -0.87500d-02
        hbs( 7, 8, 3) = -0.92900d-02
        hbs( 7, 9, 3) = 0.49700d-02
        hbs( 7,10, 3) = 0.32780d-01
        hbs( 7,11, 3) = 0.60000d-04
        hbs( 7,12, 3) = -0.10700d-02
        hbs( 8, 1, 3) = 0.36900d-02
        hbs( 8, 2, 3) = -0.18590d-01
        hbs( 8, 3, 3) = 0.13800d-02
        hbs( 8, 4, 3) = 0.25350d-01
        hbs( 8, 5, 3) = -0.35000d-02
        hbs( 8, 6, 3) = -0.42000d-03
        hbs( 8, 7, 3) = 0.45000d-02
        hbs( 8, 8, 3) = 0.50400d-02
        hbs( 8, 9, 3) = -0.20100d-02
        hbs( 8,10, 3) = -0.14680d-01
        hbs( 8,11, 3) = -0.62000d-03
        hbs( 8,12, 3) = 0.55000d-03
        hbs( 9, 1, 3) = 0.53900d-02
        hbs( 9, 2, 3) = -0.17270d-01
        hbs( 9, 3, 3) = -0.13700d-02
        hbs( 9, 4, 3) = -0.24400d-02
        hbs( 9, 5, 3) = -0.26800d-02
        hbs( 9, 6, 3) = 0.29390d-01
        hbs( 9, 7, 3) = 0.42400d-02
        hbs( 9, 8, 3) = 0.42500d-02
        hbs( 9, 9, 3) = -0.29600d-02
        hbs( 9,10, 3) = -0.18100d-01
        hbs( 9,11, 3) = 0.55000d-03
        hbs( 9,12, 3) = 0.51000d-03
        hbs(10, 1, 3) = 0.36130d-01
        hbs(10, 2, 3) = 0.38360d-01
        hbs(10, 3, 3) = 0.36300d-02
        hbs(10, 4, 3) = -0.11890d-01
        hbs(10, 5, 3) = -0.10500d-02
        hbs(10, 6, 3) = -0.26840d-01
        hbs(10, 7, 3) = 0.48100d-02
        hbs(10, 8, 3) = -0.27100d-01
        hbs(10, 9, 3) = -0.88200d-02
        hbs(10,10, 3) = -0.12140d-01
        hbs(10,11, 3) = 0.47600d-02
        hbs(10,12, 3) = -0.77000d-03
        hbs(11, 1, 3) = -0.20190d-01
        hbs(11, 2, 3) = -0.17180d-01
        hbs(11, 3, 3) = -0.21300d-02
        hbs(11, 4, 3) = 0.57400d-02
        hbs(11, 5, 3) = 0.12900d-02
        hbs(11, 6, 3) = 0.27000d-01
        hbs(11, 7, 3) = -0.28400d-02
        hbs(11, 8, 3) = 0.47000d-03
        hbs(11, 9, 3) = 0.44300d-02
        hbs(11,10, 3) = 0.60700d-02
        hbs(11,11, 3) = -0.19800d-02
        hbs(11,12, 3) = -0.18000d-03
        hbs(12, 1, 3) = -0.15940d-01
        hbs(12, 2, 3) = -0.21180d-01
        hbs(12, 3, 3) = -0.14900d-02
        hbs(12, 4, 3) = 0.61500d-02
        hbs(12, 5, 3) = -0.24000d-03
        hbs(12, 6, 3) = -0.16000d-03
        hbs(12, 7, 3) = -0.19600d-02
        hbs(12, 8, 3) = 0.26630d-01
        hbs(12, 9, 3) = 0.43900d-02
        hbs(12,10, 3) = 0.60700d-02
        hbs(12,11, 3) = -0.27800d-02
        hbs(12,12, 3) = 0.94000d-03
        hbs(13, 1, 3) = 0.35720d-01
        hbs(13, 2, 3) = -0.96900d-02
        hbs(13, 3, 3) = 0.29800d-02
        hbs(13, 4, 3) = 0.32310d-01
        hbs(13, 5, 3) = -0.94000d-03
        hbs(13, 6, 3) = -0.10140d-01
        hbs(13, 7, 3) = 0.70000d-04
        hbs(13, 8, 3) = -0.29530d-01
        hbs(13, 9, 3) = 0.50200d-02
        hbs(13,10, 3) = -0.23420d-01
        hbs(13,11, 3) = -0.86600d-02
        hbs(13,12, 3) = 0.53400d-02
        hbs(14, 1, 3) = -0.17210d-01
        hbs(14, 2, 3) = 0.56100d-02
        hbs(14, 3, 3) = -0.23300d-02
        hbs(14, 4, 3) = -0.17680d-01
        hbs(14, 5, 3) = 0.48000d-03
        hbs(14, 6, 3) = 0.47700d-02
        hbs(14, 7, 3) = 0.56000d-03
        hbs(14, 8, 3) = 0.29630d-01
        hbs(14, 9, 3) = -0.29700d-02
        hbs(14,10, 3) = -0.24000d-02
        hbs(14,11, 3) = 0.43000d-02
        hbs(14,12, 3) = -0.22600d-02
        hbs(15, 1, 3) = -0.18510d-01
        hbs(15, 2, 3) = 0.40800d-02
        hbs(15, 3, 3) = -0.65000d-03
        hbs(15, 4, 3) = -0.14630d-01
        hbs(15, 5, 3) = 0.46000d-03
        hbs(15, 6, 3) = 0.53700d-02
        hbs(15, 7, 3) = -0.63000d-03
        hbs(15, 8, 3) = -0.10000d-03
        hbs(15, 9, 3) = -0.20500d-02
        hbs(15,10, 3) = 0.25820d-01
        hbs(15,11, 3) = 0.43700d-02
        hbs(15,12, 3) = -0.30900d-02
        hbs(16, 1, 3) = -0.10200d-01
        hbs(16, 2, 3) = -0.23080d-01
        hbs(16, 3, 3) = 0.10040d-01
        hbs(16, 4, 3) = 0.41930d-01
        hbs(16, 5, 3) = 0.31000d-03
        hbs(16, 6, 3) = 0.31710d-01
        hbs(16, 7, 3) = -0.63000d-03
        hbs(16, 8, 3) = -0.56800d-02
        hbs(16, 9, 3) = -0.55000d-03
        hbs(16,10, 3) = -0.27810d-01
        hbs(16,11, 3) = 0.48500d-02
        hbs(16,12, 3) = -0.85400d-02
        hbs(17, 1, 3) = 0.41100d-02
        hbs(17, 2, 3) = 0.25000d-01
        hbs(17, 3, 3) = -0.74400d-02
        hbs(17, 4, 3) = -0.22980d-01
        hbs(17, 5, 3) = 0.60000d-03
        hbs(17, 6, 3) = -0.14860d-01
        hbs(17, 7, 3) = 0.32000d-03
        hbs(17, 8, 3) = 0.33600d-02
        hbs(17, 9, 3) = -0.32000d-03
        hbs(17,10, 3) = 0.20400d-02
        hbs(17,11, 3) = -0.21000d-02
        hbs(17,12, 3) = 0.68800d-02
        hbs(18, 1, 3) = 0.61000d-02
        hbs(18, 2, 3) = -0.19200d-02
        hbs(18, 3, 3) = -0.26000d-02
        hbs(18, 4, 3) = -0.18950d-01
        hbs(18, 5, 3) = -0.92000d-03
        hbs(18, 6, 3) = -0.16850d-01
        hbs(18, 7, 3) = 0.31000d-03
        hbs(18, 8, 3) = 0.23200d-02
        hbs(18, 9, 3) = 0.87000d-03
        hbs(18,10, 3) = 0.25770d-01
        hbs(18,11, 3) = -0.27500d-02
        hbs(18,12, 3) = 0.16600d-02
        ! Hessian elements for anchor point 4
        hss( 1, 1, 4) = 0.34759d+00
        hss( 2, 1, 4) = 0.17810d-01
        hss( 2, 2, 4) = 0.35623d+00
        hss( 3, 1, 4) = 0.29010d-01
        hss( 3, 2, 4) = 0.32710d-01
        hss( 3, 3, 4) = 0.46124d+00
        hss( 4, 1, 4) = -0.19600d-02
        hss( 4, 2, 4) = 0.14870d-01
        hss( 4, 3, 4) = -0.89100d-02
        hss( 4, 4, 4) = 0.38581d+00
        hss( 5, 1, 4) = 0.38800d-02
        hss( 5, 2, 4) = -0.77000d-03
        hss( 5, 3, 4) = 0.71700d-02
        hss( 5, 4, 4) = 0.36500d-02
        hss( 5, 5, 4) = 0.35807d+00
        hss( 6, 1, 4) = 0.33890d-01
        hss( 6, 2, 4) = 0.44570d-01
        hss( 6, 3, 4) = -0.43100d-02
        hss( 6, 4, 4) = -0.36100d-02
        hss( 6, 5, 4) = 0.55000d-03
        hss( 6, 6, 4) = 0.35931d+00
        hss( 7, 1, 4) = -0.54000d-03
        hss( 7, 2, 4) = -0.63000d-03
        hss( 7, 3, 4) = 0.57000d-03
        hss( 7, 4, 4) = 0.39700d-02
        hss( 7, 5, 4) = 0.39000d-03
        hss( 7, 6, 4) = 0.39800d-02
        hss( 7, 7, 4) = 0.37091d+00
        hss( 8, 1, 4) = 0.42030d-01
        hss( 8, 2, 4) = 0.35060d-01
        hss( 8, 3, 4) = -0.52800d-02
        hss( 8, 4, 4) = 0.24280d-01
        hss( 8, 5, 4) = -0.11100d-02
        hss( 8, 6, 4) = 0.17880d-01
        hss( 8, 7, 4) = -0.10400d-02
        hss( 8, 8, 4) = 0.35888d+00
        hss( 9, 1, 4) = -0.53000d-03
        hss( 9, 2, 4) = -0.57000d-03
        hss( 9, 3, 4) = 0.26000d-03
        hss( 9, 4, 4) = -0.39000d-03
        hss( 9, 5, 4) = 0.19000d-03
        hss( 9, 6, 4) = 0.40100d-02
        hss( 9, 7, 4) = 0.49000d-03
        hss( 9, 8, 4) = 0.39600d-02
        hss( 9, 9, 4) = 0.36356d+00
        hss(10, 1, 4) = 0.19520d-01
        hss(10, 2, 4) = -0.69300d-02
        hss(10, 3, 4) = -0.48800d-02
        hss(10, 4, 4) = 0.64010d-01
        hss(10, 5, 4) = -0.43000d-03
        hss(10, 6, 4) = 0.21230d-01
        hss(10, 7, 4) = -0.10900d-02
        hss(10, 8, 4) = -0.77000d-03
        hss(10, 9, 4) = -0.19000d-03
        hss(10,10, 4) = 0.38116d+00
        hss(11, 1, 4) = -0.69000d-03
        hss(11, 2, 4) = -0.50000d-03
        hss(11, 3, 4) = 0.52000d-03
        hss(11, 4, 4) = -0.11600d-02
        hss(11, 5, 4) = 0.80000d-04
        hss(11, 6, 4) = -0.10000d-02
        hss(11, 7, 4) = 0.22000d-03
        hss(11, 8, 4) = 0.39700d-02
        hss(11, 9, 4) = 0.48000d-03
        hss(11,10, 4) = 0.39800d-02
        hss(11,11, 4) = 0.37152d+00
        hss(12, 1, 4) = -0.16500d-02
        hss(12, 2, 4) = 0.39400d-02
        hss(12, 3, 4) = -0.19000d-03
        hss(12, 4, 4) = -0.10000d-03
        hss(12, 5, 4) = 0.35000d-03
        hss(12, 6, 4) = -0.10400d-02
        hss(12, 7, 4) = 0.30000d-04
        hss(12, 8, 4) = 0.57000d-03
        hss(12, 9, 4) = 0.12000d-03
        hss(12,10, 4) = 0.23600d-02
        hss(12,11, 4) = 0.40000d-03
        hss(12,12, 4) = 0.37254d+00
        ! Hessian elements for anchor point 4
        hbb( 1, 1, 4) = 0.82450d-01
        hbb( 2, 1, 4) = -0.46090d-01
        hbb( 2, 2, 4) = 0.14487d+00
        hbb( 3, 1, 4) = -0.36360d-01
        hbb( 3, 2, 4) = -0.98780d-01
        hbb( 3, 3, 4) = 0.13515d+00
        hbb( 4, 1, 4) = -0.36000d-01
        hbb( 4, 2, 4) = 0.23310d-01
        hbb( 4, 3, 4) = 0.12690d-01
        hbb( 4, 4, 4) = 0.66170d-01
        hbb( 5, 1, 4) = 0.20240d-01
        hbb( 5, 2, 4) = -0.12480d-01
        hbb( 5, 3, 4) = -0.77700d-02
        hbb( 5, 4, 4) = -0.33030d-01
        hbb( 5, 5, 4) = 0.71950d-01
        hbb( 6, 1, 4) = 0.15760d-01
        hbb( 6, 2, 4) = -0.10840d-01
        hbb( 6, 3, 4) = -0.49200d-02
        hbb( 6, 4, 4) = -0.33140d-01
        hbb( 6, 5, 4) = -0.38920d-01
        hbb( 6, 6, 4) = 0.72060d-01
        hbb( 7, 1, 4) = -0.71500d-02
        hbb( 7, 2, 4) = -0.39000d-02
        hbb( 7, 3, 4) = 0.11050d-01
        hbb( 7, 4, 4) = -0.25740d-01
        hbb( 7, 5, 4) = 0.93100d-02
        hbb( 7, 6, 4) = 0.16430d-01
        hbb( 7, 7, 4) = 0.64420d-01
        hbb( 8, 1, 4) = -0.22800d-02
        hbb( 8, 2, 4) = 0.74300d-02
        hbb( 8, 3, 4) = -0.51500d-02
        hbb( 8, 4, 4) = 0.15690d-01
        hbb( 8, 5, 4) = -0.49500d-02
        hbb( 8, 6, 4) = -0.10740d-01
        hbb( 8, 7, 4) = -0.31490d-01
        hbb( 8, 8, 4) = 0.72950d-01
        hbb( 9, 1, 4) = 0.94300d-02
        hbb( 9, 2, 4) = -0.35300d-02
        hbb( 9, 3, 4) = -0.59100d-02
        hbb( 9, 4, 4) = 0.10050d-01
        hbb( 9, 5, 4) = -0.43500d-02
        hbb( 9, 6, 4) = -0.56900d-02
        hbb( 9, 7, 4) = -0.32930d-01
        hbb( 9, 8, 4) = -0.41460d-01
        hbb( 9, 9, 4) = 0.74390d-01
        hbb(10, 1, 4) = 0.45500d-02
        hbb(10, 2, 4) = -0.20200d-02
        hbb(10, 3, 4) = -0.25300d-02
        hbb(10, 4, 4) = -0.64100d-02
        hbb(10, 5, 4) = 0.91900d-02
        hbb(10, 6, 4) = -0.27700d-02
        hbb(10, 7, 4) = -0.33410d-01
        hbb(10, 8, 4) = 0.13410d-01
        hbb(10, 9, 4) = 0.20000d-01
        hbb(10,10, 4) = 0.75440d-01
        hbb(11, 1, 4) = -0.22100d-02
        hbb(11, 2, 4) = 0.12300d-02
        hbb(11, 3, 4) = 0.98000d-03
        hbb(11, 4, 4) = -0.23100d-02
        hbb(11, 5, 4) = -0.32200d-02
        hbb(11, 6, 4) = 0.55300d-02
        hbb(11, 7, 4) = 0.20020d-01
        hbb(11, 8, 4) = -0.74100d-02
        hbb(11, 9, 4) = -0.12610d-01
        hbb(11,10, 4) = -0.37790d-01
        hbb(11,11, 4) = 0.78360d-01
        hbb(12, 1, 4) = -0.23400d-02
        hbb(12, 2, 4) = 0.79000d-03
        hbb(12, 3, 4) = 0.15500d-02
        hbb(12, 4, 4) = 0.87200d-02
        hbb(12, 5, 4) = -0.59700d-02
        hbb(12, 6, 4) = -0.27500d-02
        hbb(12, 7, 4) = 0.13390d-01
        hbb(12, 8, 4) = -0.60000d-02
        hbb(12, 9, 4) = -0.73900d-02
        hbb(12,10, 4) = -0.37640d-01
        hbb(12,11, 4) = -0.40570d-01
        hbb(12,12, 4) = 0.78220d-01
        hbb(13, 1, 4) = -0.67600d-02
        hbb(13, 2, 4) = 0.10710d-01
        hbb(13, 3, 4) = -0.39500d-02
        hbb(13, 4, 4) = -0.11500d-02
        hbb(13, 5, 4) = 0.12900d-02
        hbb(13, 6, 4) = -0.14000d-03
        hbb(13, 7, 4) = 0.29100d-02
        hbb(13, 8, 4) = 0.41300d-02
        hbb(13, 9, 4) = -0.70300d-02
        hbb(13,10, 4) = -0.33760d-01
        hbb(13,11, 4) = 0.13510d-01
        hbb(13,12, 4) = 0.20240d-01
        hbb(13,13, 4) = 0.65300d-01
        hbb(14, 1, 4) = 0.91200d-02
        hbb(14, 2, 4) = -0.60500d-02
        hbb(14, 3, 4) = -0.30700d-02
        hbb(14, 4, 4) = 0.66000d-03
        hbb(14, 5, 4) = -0.85000d-03
        hbb(14, 6, 4) = 0.19000d-03
        hbb(14, 7, 4) = -0.71100d-02
        hbb(14, 8, 4) = -0.51000d-03
        hbb(14, 9, 4) = 0.76200d-02
        hbb(14,10, 4) = 0.20030d-01
        hbb(14,11, 4) = -0.73800d-02
        hbb(14,12, 4) = -0.12650d-01
        hbb(14,13, 4) = -0.33180d-01
        hbb(14,14, 4) = 0.74350d-01
        hbb(15, 1, 4) = -0.23600d-02
        hbb(15, 2, 4) = -0.46500d-02
        hbb(15, 3, 4) = 0.70200d-02
        hbb(15, 4, 4) = 0.49000d-03
        hbb(15, 5, 4) = -0.44000d-03
        hbb(15, 6, 4) = -0.50000d-04
        hbb(15, 7, 4) = 0.42000d-02
        hbb(15, 8, 4) = -0.36200d-02
        hbb(15, 9, 4) = -0.59000d-03
        hbb(15,10, 4) = 0.13730d-01
        hbb(15,11, 4) = -0.61300d-02
        hbb(15,12, 4) = -0.76000d-02
        hbb(15,13, 4) = -0.32120d-01
        hbb(15,14, 4) = -0.41170d-01
        hbb(15,15, 4) = 0.73290d-01
        hbb(16, 1, 4) = -0.37090d-01
        hbb(16, 2, 4) = 0.17990d-01
        hbb(16, 3, 4) = 0.19100d-01
        hbb(16, 4, 4) = 0.31400d-02
        hbb(16, 5, 4) = -0.70000d-02
        hbb(16, 6, 4) = 0.38600d-02
        hbb(16, 7, 4) = -0.10300d-02
        hbb(16, 8, 4) = 0.54000d-03
        hbb(16, 9, 4) = 0.48000d-03
        hbb(16,10, 4) = -0.64000d-02
        hbb(16,11, 4) = 0.87800d-02
        hbb(16,12, 4) = -0.23700d-02
        hbb(16,13, 4) = -0.26530d-01
        hbb(16,14, 4) = 0.10480d-01
        hbb(16,15, 4) = 0.16060d-01
        hbb(16,16, 4) = 0.67910d-01
        hbb(17, 1, 4) = 0.22610d-01
        hbb(17, 2, 4) = -0.90200d-02
        hbb(17, 3, 4) = -0.13590d-01
        hbb(17, 4, 4) = -0.76000d-02
        hbb(17, 5, 4) = 0.77900d-02
        hbb(17, 6, 4) = -0.19000d-03
        hbb(17, 7, 4) = 0.87000d-03
        hbb(17, 8, 4) = -0.25000d-03
        hbb(17, 9, 4) = -0.62000d-03
        hbb(17,10, 4) = 0.91600d-02
        hbb(17,11, 4) = -0.60700d-02
        hbb(17,12, 4) = -0.31000d-02
        hbb(17,13, 4) = 0.96800d-02
        hbb(17,14, 4) = -0.46000d-02
        hbb(17,15, 4) = -0.50800d-02
        hbb(17,16, 4) = -0.34720d-01
        hbb(17,17, 4) = 0.72390d-01
        hbb(18, 1, 4) = 0.14470d-01
        hbb(18, 2, 4) = -0.89700d-02
        hbb(18, 3, 4) = -0.55000d-02
        hbb(18, 4, 4) = 0.44600d-02
        hbb(18, 5, 4) = -0.79000d-03
        hbb(18, 6, 4) = -0.36700d-02
        hbb(18, 7, 4) = 0.16000d-03
        hbb(18, 8, 4) = -0.29000d-03
        hbb(18, 9, 4) = 0.13000d-03
        hbb(18,10, 4) = -0.27600d-02
        hbb(18,11, 4) = -0.27100d-02
        hbb(18,12, 4) = 0.54700d-02
        hbb(18,13, 4) = 0.16860d-01
        hbb(18,14, 4) = -0.58800d-02
        hbb(18,15, 4) = -0.10980d-01
        hbb(18,16, 4) = -0.33190d-01
        hbb(18,17, 4) = -0.37670d-01
        hbb(18,18, 4) = 0.70860d-01
        ! Hessian elements for anchor point 4
        htt( 1, 1, 4) = 0.60400d-02
        htt( 2, 1, 4) = 0.38400d-02
        htt( 2, 2, 4) = 0.50900d-02
        htt( 3, 1, 4) = -0.18500d-02
        htt( 3, 2, 4) = -0.33100d-02
        htt( 3, 3, 4) = 0.75100d-02
        htt( 4, 1, 4) = -0.40500d-02
        htt( 4, 2, 4) = -0.20700d-02
        htt( 4, 3, 4) = 0.60500d-02
        htt( 4, 4, 4) = 0.80200d-02
        htt( 5, 1, 4) = -0.59200d-02
        htt( 5, 2, 4) = -0.42600d-02
        htt( 5, 3, 4) = 0.21200d-02
        htt( 5, 4, 4) = 0.37800d-02
        htt( 5, 5, 4) = 0.65600d-02
        htt( 6, 1, 4) = -0.41200d-02
        htt( 6, 2, 4) = -0.39000d-02
        htt( 6, 3, 4) = 0.26900d-02
        htt( 6, 4, 4) = 0.29100d-02
        htt( 6, 5, 4) = 0.37400d-02
        htt( 6, 6, 4) = 0.55900d-02
        htt( 7, 1, 4) = 0.15700d-02
        htt( 7, 2, 4) = 0.25400d-02
        htt( 7, 3, 4) = -0.67700d-02
        htt( 7, 4, 4) = -0.58000d-02
        htt( 7, 5, 4) = -0.10700d-02
        htt( 7, 6, 4) = -0.27200d-02
        htt( 7, 7, 4) = 0.68600d-02
        htt( 8, 1, 4) = 0.33800d-02
        htt( 8, 2, 4) = 0.29000d-02
        htt( 8, 3, 4) = -0.62000d-02
        htt( 8, 4, 4) = -0.66800d-02
        htt( 8, 5, 4) = -0.38900d-02
        htt( 8, 6, 4) = -0.88000d-03
        htt( 8, 7, 4) = 0.52100d-02
        htt( 8, 8, 4) = 0.82200d-02
        htt( 9, 1, 4) = -0.23900d-02
        htt( 9, 2, 4) = -0.68000d-03
        htt( 9, 3, 4) = 0.45000d-03
        htt( 9, 4, 4) = 0.21600d-02
        htt( 9, 5, 4) = 0.17400d-02
        htt( 9, 6, 4) = 0.14100d-02
        htt( 9, 7, 4) = -0.96000d-03
        htt( 9, 8, 4) = -0.12900d-02
        htt( 9, 9, 4) = 0.24600d-02
        htt(10, 1, 4) = -0.15600d-02
        htt(10, 2, 4) = 0.23000d-03
        htt(10, 3, 4) = 0.34000d-03
        htt(10, 4, 4) = 0.21300d-02
        htt(10, 5, 4) = 0.17600d-02
        htt(10, 6, 4) = 0.79000d-03
        htt(10, 7, 4) = -0.40000d-04
        htt(10, 8, 4) = -0.10100d-02
        htt(10, 9, 4) = 0.57000d-03
        htt(10,10, 4) = 0.39100d-02
        htt(11, 1, 4) = -0.10000d-03
        htt(11, 2, 4) = -0.19700d-02
        htt(11, 3, 4) = 0.19800d-02
        htt(11, 4, 4) = 0.11000d-03
        htt(11, 5, 4) = 0.10000d-04
        htt(11, 6, 4) = 0.11900d-02
        htt(11, 7, 4) = -0.19700d-02
        htt(11, 8, 4) = -0.79000d-03
        htt(11, 9, 4) = 0.68000d-03
        htt(11,10, 4) = -0.12900d-02
        htt(11,11, 4) = 0.26200d-02
        htt(12, 1, 4) = 0.72000d-03
        htt(12, 2, 4) = -0.10500d-02
        htt(12, 3, 4) = 0.18600d-02
        htt(12, 4, 4) = 0.90000d-04
        htt(12, 5, 4) = 0.30000d-04
        htt(12, 6, 4) = 0.57000d-03
        htt(12, 7, 4) = -0.10500d-02
        htt(12, 8, 4) = -0.51000d-03
        htt(12, 9, 4) = -0.12000d-02
        htt(12,10, 4) = 0.20500d-02
        htt(12,11, 4) = 0.64000d-03
        htt(12,12, 4) = 0.39000d-02
        htt(13, 1, 4) = -0.10900d-02
        htt(13, 2, 4) = -0.19400d-02
        htt(13, 3, 4) = 0.59000d-03
        htt(13, 4, 4) = -0.26000d-03
        htt(13, 5, 4) = 0.15200d-02
        htt(13, 6, 4) = 0.15100d-02
        htt(13, 7, 4) = -0.70000d-04
        htt(13, 8, 4) = -0.90000d-04
        htt(13, 9, 4) = -0.21100d-02
        htt(13,10, 4) = 0.12000d-03
        htt(13,11, 4) = -0.12300d-02
        htt(13,12, 4) = 0.10000d-02
        htt(13,13, 4) = 0.50300d-02
        htt(14, 1, 4) = 0.70000d-03
        htt(14, 2, 4) = 0.53000d-03
        htt(14, 3, 4) = -0.31000d-03
        htt(14, 4, 4) = -0.47000d-03
        htt(14, 5, 4) = -0.22000d-03
        htt(14, 6, 4) = -0.80000d-03
        htt(14, 7, 4) = 0.74000d-03
        htt(14, 8, 4) = 0.16000d-03
        htt(14, 9, 4) = -0.19000d-03
        htt(14,10, 4) = 0.13500d-02
        htt(14,11, 4) = -0.10000d-04
        htt(14,12, 4) = 0.15300d-02
        htt(14,13, 4) = -0.81000d-03
        htt(14,14, 4) = 0.39200d-02
        htt(15, 1, 4) = -0.19100d-02
        htt(15, 2, 4) = -0.28400d-02
        htt(15, 3, 4) = 0.70000d-03
        htt(15, 4, 4) = -0.23000d-03
        htt(15, 5, 4) = 0.15000d-02
        htt(15, 6, 4) = 0.21100d-02
        htt(15, 7, 4) = -0.98000d-03
        htt(15, 8, 4) = -0.37000d-03
        htt(15, 9, 4) = -0.24000d-03
        htt(15,10, 4) = -0.31800d-02
        htt(15,11, 4) = 0.73000d-03
        htt(15,12, 4) = -0.22200d-02
        htt(15,13, 4) = 0.28200d-02
        htt(15,14, 4) = -0.23300d-02
        htt(15,15, 4) = 0.57300d-02
        htt(16, 1, 4) = -0.12000d-03
        htt(16, 2, 4) = -0.37000d-03
        htt(16, 3, 4) = -0.19000d-03
        htt(16, 4, 4) = -0.45000d-03
        htt(16, 5, 4) = -0.25000d-03
        htt(16, 6, 4) = -0.19000d-03
        htt(16, 7, 4) = -0.17000d-03
        htt(16, 8, 4) = -0.12000d-03
        htt(16, 9, 4) = 0.16800d-02
        htt(16,10, 4) = -0.19500d-02
        htt(16,11, 4) = 0.19400d-02
        htt(16,12, 4) = -0.16900d-02
        htt(16,13, 4) = -0.30200d-02
        htt(16,14, 4) = 0.24000d-02
        htt(16,15, 4) = 0.58000d-03
        htt(16,16, 4) = 0.60000d-02
        htt(17, 1, 4) = 0.12100d-02
        htt(17, 2, 4) = 0.15200d-02
        htt(17, 3, 4) = -0.31000d-03
        htt(17, 4, 4) = 0.00000d+00
        htt(17, 5, 4) = -0.86000d-03
        htt(17, 6, 4) = -0.18800d-02
        htt(17, 7, 4) = 0.58000d-03
        htt(17, 8, 4) = -0.44000d-03
        htt(17, 9, 4) = 0.14500d-02
        htt(17,10, 4) = 0.90000d-04
        htt(17,11, 4) = 0.11300d-02
        htt(17,12, 4) = -0.23000d-03
        htt(17,13, 4) = -0.46000d-02
        htt(17,14, 4) = 0.13000d-02
        htt(17,15, 4) = -0.32400d-02
        htt(17,16, 4) = 0.26500d-02
        htt(17,17, 4) = 0.49600d-02
        htt(18, 1, 4) = 0.14500d-02
        htt(18, 2, 4) = 0.22200d-02
        htt(18, 3, 4) = -0.10600d-02
        htt(18, 4, 4) = -0.29000d-03
        htt(18, 5, 4) = -0.17100d-02
        htt(18, 6, 4) = -0.27100d-02
        htt(18, 7, 4) = 0.67000d-03
        htt(18, 8, 4) = -0.32000d-03
        htt(18, 9, 4) = 0.10000d-02
        htt(18,10, 4) = 0.82000d-03
        htt(18,11, 4) = 0.19000d-03
        htt(18,12, 4) = 0.10000d-04
        htt(18,13, 4) = -0.32400d-02
        htt(18,14, 4) = 0.19000d-02
        htt(18,15, 4) = -0.30600d-02
        htt(18,16, 4) = 0.20800d-02
        htt(18,17, 4) = 0.29700d-02
        htt(18,18, 4) = 0.56300d-02
        htt(19, 1, 4) = -0.59000d-03
        htt(19, 2, 4) = -0.96000d-03
        htt(19, 3, 4) = 0.59000d-03
        htt(19, 4, 4) = 0.21000d-03
        htt(19, 5, 4) = 0.89000d-03
        htt(19, 6, 4) = 0.43000d-03
        htt(19, 7, 4) = -0.23000d-03
        htt(19, 8, 4) = -0.69000d-03
        htt(19, 9, 4) = -0.48000d-03
        htt(19,10, 4) = -0.11500d-02
        htt(19,11, 4) = -0.90000d-04
        htt(19,12, 4) = -0.76000d-03
        htt(19,13, 4) = 0.12600d-02
        htt(19,14, 4) = -0.34500d-02
        htt(19,15, 4) = 0.19200d-02
        htt(19,16, 4) = -0.27800d-02
        htt(19,17, 4) = -0.95000d-03
        htt(19,18, 4) = -0.21800d-02
        htt(19,19, 4) = 0.37700d-02
        htt(20, 1, 4) = -0.35000d-03
        htt(20, 2, 4) = -0.26000d-03
        htt(20, 3, 4) = -0.16000d-03
        htt(20, 4, 4) = -0.80000d-04
        htt(20, 5, 4) = 0.40000d-04
        htt(20, 6, 4) = -0.40000d-03
        htt(20, 7, 4) = -0.14000d-03
        htt(20, 8, 4) = -0.57000d-03
        htt(20, 9, 4) = -0.93000d-03
        htt(20,10, 4) = -0.42000d-03
        htt(20,11, 4) = -0.10200d-02
        htt(20,12, 4) = -0.51000d-03
        htt(20,13, 4) = 0.26100d-02
        htt(20,14, 4) = -0.28400d-02
        htt(20,15, 4) = 0.21100d-02
        htt(20,16, 4) = -0.33500d-02
        htt(20,17, 4) = -0.29400d-02
        htt(20,18, 4) = 0.47000d-03
        htt(20,19, 4) = 0.25400d-02
        htt(20,20, 4) = 0.59400d-02
        htt(21, 1, 4) = 0.21700d-02
        htt(21, 2, 4) = 0.15000d-02
        htt(21, 3, 4) = -0.99000d-03
        htt(21, 4, 4) = -0.16600d-02
        htt(21, 5, 4) = -0.30200d-02
        htt(21, 6, 4) = -0.69000d-03
        htt(21, 7, 4) = -0.20000d-04
        htt(21, 8, 4) = 0.23100d-02
        htt(21, 9, 4) = -0.11900d-02
        htt(21,10, 4) = -0.98000d-03
        htt(21,11, 4) = -0.50000d-03
        htt(21,12, 4) = -0.28000d-03
        htt(21,13, 4) = 0.12600d-02
        htt(21,14, 4) = -0.76000d-03
        htt(21,15, 4) = 0.10500d-02
        htt(21,16, 4) = -0.97000d-03
        htt(21,17, 4) = -0.21400d-02
        htt(21,18, 4) = -0.48000d-03
        htt(21,19, 4) = -0.11000d-03
        htt(21,20, 4) = 0.15500d-02
        htt(21,21, 4) = 0.28900d-02
        htt(22, 1, 4) = 0.27000d-03
        htt(22, 2, 4) = 0.11300d-02
        htt(22, 3, 4) = -0.15900d-02
        htt(22, 4, 4) = -0.74000d-03
        htt(22, 5, 4) = -0.60000d-04
        htt(22, 6, 4) = -0.26200d-02
        htt(22, 7, 4) = 0.17100d-02
        htt(22, 8, 4) = -0.85000d-03
        htt(22, 9, 4) = -0.85000d-03
        htt(22,10, 4) = 0.40000d-04
        htt(22,11, 4) = -0.17300d-02
        htt(22,12, 4) = -0.85000d-03
        htt(22,13, 4) = 0.12800d-02
        htt(22,14, 4) = -0.15000d-03
        htt(22,15, 4) = 0.40000d-03
        htt(22,16, 4) = -0.10300d-02
        htt(22,17, 4) = -0.10700d-02
        htt(22,18, 4) = 0.56000d-03
        htt(22,19, 4) = 0.37000d-03
        htt(22,20, 4) = 0.20000d-02
        htt(22,21, 4) = 0.44000d-03
        htt(22,22, 4) = 0.31200d-02
        htt(23, 1, 4) = 0.19200d-02
        htt(23, 2, 4) = 0.79000d-03
        htt(23, 3, 4) = -0.23000d-03
        htt(23, 4, 4) = -0.13700d-02
        htt(23, 5, 4) = -0.21700d-02
        htt(23, 6, 4) = 0.14000d-03
        htt(23, 7, 4) = -0.12000d-03
        htt(23, 8, 4) = 0.21900d-02
        htt(23, 9, 4) = -0.74000d-03
        htt(23,10, 4) = -0.17100d-02
        htt(23,11, 4) = 0.45000d-03
        htt(23,12, 4) = -0.53000d-03
        htt(23,13, 4) = -0.10000d-03
        htt(23,14, 4) = -0.13700d-02
        htt(23,15, 4) = 0.87000d-03
        htt(23,16, 4) = -0.40000d-03
        htt(23,17, 4) = -0.15000d-03
        htt(23,18, 4) = -0.31500d-02
        htt(23,19, 4) = 0.11200d-02
        htt(23,20, 4) = -0.18800d-02
        htt(23,21, 4) = 0.12200d-02
        htt(23,22, 4) = -0.12000d-02
        htt(23,23, 4) = 0.42400d-02
        htt(24, 1, 4) = 0.30000d-04
        htt(24, 2, 4) = 0.42000d-03
        htt(24, 3, 4) = -0.84000d-03
        htt(24, 4, 4) = -0.45000d-03
        htt(24, 5, 4) = 0.79000d-03
        htt(24, 6, 4) = -0.18000d-02
        htt(24, 7, 4) = 0.16200d-02
        htt(24, 8, 4) = -0.97000d-03
        htt(24, 9, 4) = -0.39000d-03
        htt(24,10, 4) = -0.70000d-03
        htt(24,11, 4) = -0.79000d-03
        htt(24,12, 4) = -0.11000d-02
        htt(24,13, 4) = -0.80000d-04
        htt(24,14, 4) = -0.76000d-03
        htt(24,15, 4) = 0.22000d-03
        htt(24,16, 4) = -0.45000d-03
        htt(24,17, 4) = 0.92000d-03
        htt(24,18, 4) = -0.21000d-02
        htt(24,19, 4) = 0.16000d-02
        htt(24,20, 4) = -0.14300d-02
        htt(24,21, 4) = -0.12300d-02
        htt(24,22, 4) = 0.14800d-02
        htt(24,23, 4) = 0.18100d-02
        htt(24,24, 4) = 0.45300d-02
        ! Hessian elements for anchor point 4
        hbs( 1, 1, 4) = -0.31760d-01
        hbs( 1, 2, 4) = -0.29360d-01
        hbs( 1, 3, 4) = -0.31820d-01
        hbs( 1, 4, 4) = -0.95900d-02
        hbs( 1, 5, 4) = 0.52300d-02
        hbs( 1, 6, 4) = 0.39930d-01
        hbs( 1, 7, 4) = -0.45000d-03
        hbs( 1, 8, 4) = 0.40450d-01
        hbs( 1, 9, 4) = -0.18000d-03
        hbs( 1,10, 4) = -0.96900d-02
        hbs( 1,11, 4) = -0.42000d-03
        hbs( 1,12, 4) = 0.47400d-02
        hbs( 2, 1, 4) = 0.49930d-01
        hbs( 2, 2, 4) = -0.13930d-01
        hbs( 2, 3, 4) = 0.99600d-02
        hbs( 2, 4, 4) = 0.18000d-02
        hbs( 2, 5, 4) = -0.12200d-02
        hbs( 2, 6, 4) = -0.27810d-01
        hbs( 2, 7, 4) = 0.11800d-02
        hbs( 2, 8, 4) = -0.14650d-01
        hbs( 2, 9, 4) = 0.00000d+00
        hbs( 2,10, 4) = 0.41300d-02
        hbs( 2,11, 4) = -0.69000d-03
        hbs( 2,12, 4) = -0.89000d-03
        hbs( 3, 1, 4) = -0.18180d-01
        hbs( 3, 2, 4) = 0.43280d-01
        hbs( 3, 3, 4) = 0.21860d-01
        hbs( 3, 4, 4) = 0.77900d-02
        hbs( 3, 5, 4) = -0.40100d-02
        hbs( 3, 6, 4) = -0.12120d-01
        hbs( 3, 7, 4) = -0.73000d-03
        hbs( 3, 8, 4) = -0.25800d-01
        hbs( 3, 9, 4) = 0.18000d-03
        hbs( 3,10, 4) = 0.55600d-02
        hbs( 3,11, 4) = 0.11100d-02
        hbs( 3,12, 4) = -0.38500d-02
        hbs( 4, 1, 4) = -0.20820d-01
        hbs( 4, 2, 4) = -0.12100d-01
        hbs( 4, 3, 4) = 0.15180d-01
        hbs( 4, 4, 4) = -0.29840d-01
        hbs( 4, 5, 4) = -0.97300d-02
        hbs( 4, 6, 4) = -0.57000d-02
        hbs( 4, 7, 4) = 0.49500d-02
        hbs( 4, 8, 4) = 0.31150d-01
        hbs( 4, 9, 4) = -0.44000d-03
        hbs( 4,10, 4) = 0.40270d-01
        hbs( 4,11, 4) = -0.59000d-03
        hbs( 4,12, 4) = 0.29000d-03
        hbs( 5, 1, 4) = 0.26090d-01
        hbs( 5, 2, 4) = 0.58900d-02
        hbs( 5, 3, 4) = -0.14520d-01
        hbs( 5, 4, 4) = 0.25400d-02
        hbs( 5, 5, 4) = 0.69400d-02
        hbs( 5, 6, 4) = 0.30900d-02
        hbs( 5, 7, 4) = -0.21700d-02
        hbs( 5, 8, 4) = -0.15190d-01
        hbs( 5, 9, 4) = -0.39000d-03
        hbs( 5,10, 4) = -0.22140d-01
        hbs( 5,11, 4) = 0.28000d-03
        hbs( 5,12, 4) = 0.90000d-04
        hbs( 6, 1, 4) = -0.52700d-02
        hbs( 6, 2, 4) = 0.62100d-02
        hbs( 6, 3, 4) = -0.66000d-03
        hbs( 6, 4, 4) = 0.27300d-01
        hbs( 6, 5, 4) = 0.27900d-02
        hbs( 6, 6, 4) = 0.26100d-02
        hbs( 6, 7, 4) = -0.27800d-02
        hbs( 6, 8, 4) = -0.15960d-01
        hbs( 6, 9, 4) = 0.83000d-03
        hbs( 6,10, 4) = -0.18130d-01
        hbs( 6,11, 4) = 0.31000d-03
        hbs( 6,12, 4) = -0.38000d-03
        hbs( 7, 1, 4) = -0.90800d-02
        hbs( 7, 2, 4) = 0.35870d-01
        hbs( 7, 3, 4) = -0.10000d-04
        hbs( 7, 4, 4) = -0.22920d-01
        hbs( 7, 5, 4) = 0.61900d-02
        hbs( 7, 6, 4) = -0.28970d-01
        hbs( 7, 7, 4) = -0.87500d-02
        hbs( 7, 8, 4) = -0.92900d-02
        hbs( 7, 9, 4) = 0.49700d-02
        hbs( 7,10, 4) = 0.32780d-01
        hbs( 7,11, 4) = 0.60000d-04
        hbs( 7,12, 4) = -0.10700d-02
        hbs( 8, 1, 4) = 0.36900d-02
        hbs( 8, 2, 4) = -0.18590d-01
        hbs( 8, 3, 4) = 0.13800d-02
        hbs( 8, 4, 4) = 0.25350d-01
        hbs( 8, 5, 4) = -0.35000d-02
        hbs( 8, 6, 4) = -0.42000d-03
        hbs( 8, 7, 4) = 0.45000d-02
        hbs( 8, 8, 4) = 0.50400d-02
        hbs( 8, 9, 4) = -0.20100d-02
        hbs( 8,10, 4) = -0.14680d-01
        hbs( 8,11, 4) = -0.62000d-03
        hbs( 8,12, 4) = 0.55000d-03
        hbs( 9, 1, 4) = 0.53900d-02
        hbs( 9, 2, 4) = -0.17270d-01
        hbs( 9, 3, 4) = -0.13700d-02
        hbs( 9, 4, 4) = -0.24400d-02
        hbs( 9, 5, 4) = -0.26800d-02
        hbs( 9, 6, 4) = 0.29390d-01
        hbs( 9, 7, 4) = 0.42400d-02
        hbs( 9, 8, 4) = 0.42500d-02
        hbs( 9, 9, 4) = -0.29600d-02
        hbs( 9,10, 4) = -0.18100d-01
        hbs( 9,11, 4) = 0.55000d-03
        hbs( 9,12, 4) = 0.51000d-03
        hbs(10, 1, 4) = 0.36130d-01
        hbs(10, 2, 4) = 0.38360d-01
        hbs(10, 3, 4) = 0.36300d-02
        hbs(10, 4, 4) = -0.11890d-01
        hbs(10, 5, 4) = -0.10500d-02
        hbs(10, 6, 4) = -0.26840d-01
        hbs(10, 7, 4) = 0.48100d-02
        hbs(10, 8, 4) = -0.27100d-01
        hbs(10, 9, 4) = -0.88200d-02
        hbs(10,10, 4) = -0.12140d-01
        hbs(10,11, 4) = 0.47600d-02
        hbs(10,12, 4) = -0.77000d-03
        hbs(11, 1, 4) = -0.20190d-01
        hbs(11, 2, 4) = -0.17180d-01
        hbs(11, 3, 4) = -0.21300d-02
        hbs(11, 4, 4) = 0.57400d-02
        hbs(11, 5, 4) = 0.12900d-02
        hbs(11, 6, 4) = 0.27000d-01
        hbs(11, 7, 4) = -0.28400d-02
        hbs(11, 8, 4) = 0.47000d-03
        hbs(11, 9, 4) = 0.44300d-02
        hbs(11,10, 4) = 0.60700d-02
        hbs(11,11, 4) = -0.19800d-02
        hbs(11,12, 4) = -0.18000d-03
        hbs(12, 1, 4) = -0.15940d-01
        hbs(12, 2, 4) = -0.21180d-01
        hbs(12, 3, 4) = -0.14900d-02
        hbs(12, 4, 4) = 0.61500d-02
        hbs(12, 5, 4) = -0.24000d-03
        hbs(12, 6, 4) = -0.16000d-03
        hbs(12, 7, 4) = -0.19600d-02
        hbs(12, 8, 4) = 0.26630d-01
        hbs(12, 9, 4) = 0.43900d-02
        hbs(12,10, 4) = 0.60700d-02
        hbs(12,11, 4) = -0.27800d-02
        hbs(12,12, 4) = 0.94000d-03
        hbs(13, 1, 4) = 0.35720d-01
        hbs(13, 2, 4) = -0.96900d-02
        hbs(13, 3, 4) = 0.29800d-02
        hbs(13, 4, 4) = 0.32310d-01
        hbs(13, 5, 4) = -0.94000d-03
        hbs(13, 6, 4) = -0.10140d-01
        hbs(13, 7, 4) = 0.70000d-04
        hbs(13, 8, 4) = -0.29530d-01
        hbs(13, 9, 4) = 0.50200d-02
        hbs(13,10, 4) = -0.23420d-01
        hbs(13,11, 4) = -0.86600d-02
        hbs(13,12, 4) = 0.53400d-02
        hbs(14, 1, 4) = -0.17210d-01
        hbs(14, 2, 4) = 0.56100d-02
        hbs(14, 3, 4) = -0.23300d-02
        hbs(14, 4, 4) = -0.17680d-01
        hbs(14, 5, 4) = 0.48000d-03
        hbs(14, 6, 4) = 0.47700d-02
        hbs(14, 7, 4) = 0.56000d-03
        hbs(14, 8, 4) = 0.29630d-01
        hbs(14, 9, 4) = -0.29700d-02
        hbs(14,10, 4) = -0.24000d-02
        hbs(14,11, 4) = 0.43000d-02
        hbs(14,12, 4) = -0.22600d-02
        hbs(15, 1, 4) = -0.18510d-01
        hbs(15, 2, 4) = 0.40800d-02
        hbs(15, 3, 4) = -0.65000d-03
        hbs(15, 4, 4) = -0.14630d-01
        hbs(15, 5, 4) = 0.46000d-03
        hbs(15, 6, 4) = 0.53700d-02
        hbs(15, 7, 4) = -0.63000d-03
        hbs(15, 8, 4) = -0.10000d-03
        hbs(15, 9, 4) = -0.20500d-02
        hbs(15,10, 4) = 0.25820d-01
        hbs(15,11, 4) = 0.43700d-02
        hbs(15,12, 4) = -0.30900d-02
        hbs(16, 1, 4) = -0.10200d-01
        hbs(16, 2, 4) = -0.23080d-01
        hbs(16, 3, 4) = 0.10040d-01
        hbs(16, 4, 4) = 0.41930d-01
        hbs(16, 5, 4) = 0.31000d-03
        hbs(16, 6, 4) = 0.31710d-01
        hbs(16, 7, 4) = -0.63000d-03
        hbs(16, 8, 4) = -0.56800d-02
        hbs(16, 9, 4) = -0.55000d-03
        hbs(16,10, 4) = -0.27810d-01
        hbs(16,11, 4) = 0.48500d-02
        hbs(16,12, 4) = -0.85400d-02
        hbs(17, 1, 4) = 0.41100d-02
        hbs(17, 2, 4) = 0.25000d-01
        hbs(17, 3, 4) = -0.74400d-02
        hbs(17, 4, 4) = -0.22980d-01
        hbs(17, 5, 4) = 0.60000d-03
        hbs(17, 6, 4) = -0.14860d-01
        hbs(17, 7, 4) = 0.32000d-03
        hbs(17, 8, 4) = 0.33600d-02
        hbs(17, 9, 4) = -0.32000d-03
        hbs(17,10, 4) = 0.20400d-02
        hbs(17,11, 4) = -0.21000d-02
        hbs(17,12, 4) = 0.68800d-02
        hbs(18, 1, 4) = 0.61000d-02
        hbs(18, 2, 4) = -0.19200d-02
        hbs(18, 3, 4) = -0.26000d-02
        hbs(18, 4, 4) = -0.18950d-01
        hbs(18, 5, 4) = -0.92000d-03
        hbs(18, 6, 4) = -0.16850d-01
        hbs(18, 7, 4) = 0.31000d-03
        hbs(18, 8, 4) = 0.23200d-02
        hbs(18, 9, 4) = 0.87000d-03
        hbs(18,10, 4) = 0.25770d-01
        hbs(18,11, 4) = -0.27500d-02
        hbs(18,12, 4) = 0.16600d-02
    else if(ii .eq. 33) then
        ! Hessian elements for anchor point 1
        hss( 1, 1, 1) = 0.32935d+00
        hss( 2, 1, 1) = 0.10630d-01
        hss( 2, 2, 1) = 0.33545d+00
        hss( 3, 1, 1) = 0.41480d-01
        hss( 3, 2, 1) = 0.47860d-01
        hss( 3, 3, 1) = 0.60075d+00
        hss( 4, 1, 1) = -0.28200d-02
        hss( 4, 2, 1) = 0.17340d-01
        hss( 4, 3, 1) = -0.51100d-02
        hss( 4, 4, 1) = 0.44790d+00
        hss( 5, 1, 1) = 0.39700d-02
        hss( 5, 2, 1) = -0.41000d-03
        hss( 5, 3, 1) = 0.20000d-03
        hss( 5, 4, 1) = 0.13100d-02
        hss( 5, 5, 1) = 0.36210d+00
        hss( 6, 1, 1) = 0.23160d-01
        hss( 6, 2, 1) = 0.46150d-01
        hss( 6, 3, 1) = 0.31100d-02
        hss( 6, 4, 1) = 0.42160d-01
        hss( 6, 5, 1) = -0.19000d-03
        hss( 6, 6, 1) = 0.33957d+00
        hss( 7, 1, 1) = 0.14000d-03
        hss( 7, 2, 1) = -0.96000d-03
        hss( 7, 3, 1) = 0.11700d-02
        hss( 7, 4, 1) = 0.55900d-02
        hss( 7, 5, 1) = -0.68000d-03
        hss( 7, 6, 1) = 0.28700d-02
        hss( 7, 7, 1) = 0.36718d+00
        hss( 8, 1, 1) = 0.48190d-01
        hss( 8, 2, 1) = 0.21300d-01
        hss( 8, 3, 1) = 0.10200d-02
        hss( 8, 4, 1) = -0.55600d-02
        hss( 8, 5, 1) = -0.21100d-02
        hss( 8, 6, 1) = 0.34440d-01
        hss( 8, 7, 1) = -0.39000d-03
        hss( 8, 8, 1) = 0.33723d+00
        hss( 9, 1, 1) = -0.13000d-03
        hss( 9, 2, 1) = 0.50000d-04
        hss( 9, 3, 1) = -0.53000d-03
        hss( 9, 4, 1) = -0.54000d-03
        hss( 9, 5, 1) = -0.64000d-03
        hss( 9, 6, 1) = 0.37800d-02
        hss( 9, 7, 1) = 0.12000d-03
        hss( 9, 8, 1) = 0.36400d-02
        hss( 9, 9, 1) = 0.36496d+00
        hss(10, 1, 1) = 0.14800d-01
        hss(10, 2, 1) = 0.14000d-03
        hss(10, 3, 1) = -0.42300d-02
        hss(10, 4, 1) = 0.81590d-01
        hss(10, 5, 1) = -0.48000d-03
        hss(10, 6, 1) = -0.32700d-02
        hss(10, 7, 1) = -0.71000d-03
        hss(10, 8, 1) = 0.36040d-01
        hss(10, 9, 1) = -0.87000d-03
        hss(10,10, 1) = 0.46416d+00
        hss(11, 1, 1) = -0.11800d-02
        hss(11, 2, 1) = -0.64000d-03
        hss(11, 3, 1) = 0.90000d-03
        hss(11, 4, 1) = -0.78000d-03
        hss(11, 5, 1) = 0.18000d-03
        hss(11, 6, 1) = -0.11000d-03
        hss(11, 7, 1) = 0.12000d-03
        hss(11, 8, 1) = 0.28600d-02
        hss(11, 9, 1) = 0.42000d-03
        hss(11,10, 1) = 0.45200d-02
        hss(11,11, 1) = 0.36758d+00
        hss(12, 1, 1) = 0.23000d-03
        hss(12, 2, 1) = 0.14400d-02
        hss(12, 3, 1) = 0.13000d-03
        hss(12, 4, 1) = -0.34000d-03
        hss(12, 5, 1) = -0.37000d-03
        hss(12, 6, 1) = -0.80000d-03
        hss(12, 7, 1) = 0.40000d-04
        hss(12, 8, 1) = 0.22000d-03
        hss(12, 9, 1) = -0.10000d-03
        hss(12,10, 1) = 0.23100d-02
        hss(12,11, 1) = 0.17000d-03
        hss(12,12, 1) = 0.37062d+00
        ! Hessian elements for anchor point 1
        hbb( 1, 1, 1) = 0.82090d-01
        hbb( 2, 1, 1) = -0.40800d-01
        hbb( 2, 2, 1) = 0.13448d+00
        hbb( 3, 1, 1) = -0.41290d-01
        hbb( 3, 2, 1) = -0.93680d-01
        hbb( 3, 3, 1) = 0.13497d+00
        hbb( 4, 1, 1) = -0.37560d-01
        hbb( 4, 2, 1) = 0.25500d-01
        hbb( 4, 3, 1) = 0.12070d-01
        hbb( 4, 4, 1) = 0.76550d-01
        hbb( 5, 1, 1) = 0.22450d-01
        hbb( 5, 2, 1) = -0.15660d-01
        hbb( 5, 3, 1) = -0.67900d-02
        hbb( 5, 4, 1) = -0.39660d-01
        hbb( 5, 5, 1) = 0.77760d-01
        hbb( 6, 1, 1) = 0.15120d-01
        hbb( 6, 2, 1) = -0.98400d-02
        hbb( 6, 3, 1) = -0.52700d-02
        hbb( 6, 4, 1) = -0.36880d-01
        hbb( 6, 5, 1) = -0.38100d-01
        hbb( 6, 6, 1) = 0.74990d-01
        hbb( 7, 1, 1) = -0.83500d-02
        hbb( 7, 2, 1) = -0.60600d-02
        hbb( 7, 3, 1) = 0.14410d-01
        hbb( 7, 4, 1) = -0.28860d-01
        hbb( 7, 5, 1) = 0.12030d-01
        hbb( 7, 6, 1) = 0.16840d-01
        hbb( 7, 7, 1) = 0.74530d-01
        hbb( 8, 1, 1) = -0.15000d-02
        hbb( 8, 2, 1) = 0.75600d-02
        hbb( 8, 3, 1) = -0.60500d-02
        hbb( 8, 4, 1) = 0.17370d-01
        hbb( 8, 5, 1) = -0.66700d-02
        hbb( 8, 6, 1) = -0.10710d-01
        hbb( 8, 7, 1) = -0.36970d-01
        hbb( 8, 8, 1) = 0.78530d-01
        hbb( 9, 1, 1) = 0.98600d-02
        hbb( 9, 2, 1) = -0.15000d-02
        hbb( 9, 3, 1) = -0.83500d-02
        hbb( 9, 4, 1) = 0.11490d-01
        hbb( 9, 5, 1) = -0.53600d-02
        hbb( 9, 6, 1) = -0.61300d-02
        hbb( 9, 7, 1) = -0.37550d-01
        hbb( 9, 8, 1) = -0.41560d-01
        hbb( 9, 9, 1) = 0.79110d-01
        hbb(10, 1, 1) = 0.88000d-02
        hbb(10, 2, 1) = -0.45300d-02
        hbb(10, 3, 1) = -0.42700d-02
        hbb(10, 4, 1) = -0.86000d-02
        hbb(10, 5, 1) = 0.10210d-01
        hbb(10, 6, 1) = -0.16100d-02
        hbb(10, 7, 1) = -0.33520d-01
        hbb(10, 8, 1) = 0.13620d-01
        hbb(10, 9, 1) = 0.19900d-01
        hbb(10,10, 1) = 0.75620d-01
        hbb(11, 1, 1) = -0.44200d-02
        hbb(11, 2, 1) = 0.26700d-02
        hbb(11, 3, 1) = 0.17500d-02
        hbb(11, 4, 1) = -0.13500d-02
        hbb(11, 5, 1) = -0.34500d-02
        hbb(11, 6, 1) = 0.48000d-02
        hbb(11, 7, 1) = 0.19740d-01
        hbb(11, 8, 1) = -0.73600d-02
        hbb(11, 9, 1) = -0.12380d-01
        hbb(11,10, 1) = -0.37740d-01
        hbb(11,11, 1) = 0.80040d-01
        hbb(12, 1, 1) = -0.43800d-02
        hbb(12, 2, 1) = 0.18600d-02
        hbb(12, 3, 1) = 0.25200d-02
        hbb(12, 4, 1) = 0.99500d-02
        hbb(12, 5, 1) = -0.67600d-02
        hbb(12, 6, 1) = -0.31900d-02
        hbb(12, 7, 1) = 0.13790d-01
        hbb(12, 8, 1) = -0.62600d-02
        hbb(12, 9, 1) = -0.75300d-02
        hbb(12,10, 1) = -0.37880d-01
        hbb(12,11, 1) = -0.42300d-01
        hbb(12,12, 1) = 0.80170d-01
        hbb(13, 1, 1) = -0.76900d-02
        hbb(13, 2, 1) = 0.14150d-01
        hbb(13, 3, 1) = -0.64500d-02
        hbb(13, 4, 1) = 0.32300d-02
        hbb(13, 5, 1) = -0.14400d-02
        hbb(13, 6, 1) = -0.17900d-02
        hbb(13, 7, 1) = -0.67100d-02
        hbb(13, 8, 1) = 0.90400d-02
        hbb(13, 9, 1) = -0.23400d-02
        hbb(13,10, 1) = -0.33330d-01
        hbb(13,11, 1) = 0.13660d-01
        hbb(13,12, 1) = 0.19670d-01
        hbb(13,13, 1) = 0.73530d-01
        hbb(14, 1, 1) = 0.96900d-02
        hbb(14, 2, 1) = -0.81500d-02
        hbb(14, 3, 1) = -0.15500d-02
        hbb(14, 4, 1) = -0.17500d-02
        hbb(14, 5, 1) = 0.64000d-03
        hbb(14, 6, 1) = 0.11000d-02
        hbb(14, 7, 1) = -0.22600d-02
        hbb(14, 8, 1) = -0.29800d-02
        hbb(14, 9, 1) = 0.52500d-02
        hbb(14,10, 1) = 0.19960d-01
        hbb(14,11, 1) = -0.75600d-02
        hbb(14,12, 1) = -0.12400d-01
        hbb(14,13, 1) = -0.37280d-01
        hbb(14,14, 1) = 0.78950d-01
        hbb(15, 1, 1) = -0.20000d-02
        hbb(15, 2, 1) = -0.60000d-02
        hbb(15, 3, 1) = 0.80000d-02
        hbb(15, 4, 1) = -0.14800d-02
        hbb(15, 5, 1) = 0.80000d-03
        hbb(15, 6, 1) = 0.69000d-03
        hbb(15, 7, 1) = 0.89700d-02
        hbb(15, 8, 1) = -0.60600d-02
        hbb(15, 9, 1) = -0.29100d-02
        hbb(15,10, 1) = 0.13370d-01
        hbb(15,11, 1) = -0.61000d-02
        hbb(15,12, 1) = -0.72700d-02
        hbb(15,13, 1) = -0.36250d-01
        hbb(15,14, 1) = -0.41670d-01
        hbb(15,15, 1) = 0.77920d-01
        hbb(16, 1, 1) = -0.37280d-01
        hbb(16, 2, 1) = 0.11740d-01
        hbb(16, 3, 1) = 0.25540d-01
        hbb(16, 4, 1) = -0.47500d-02
        hbb(16, 5, 1) = -0.35800d-02
        hbb(16, 6, 1) = 0.83300d-02
        hbb(16, 7, 1) = 0.29200d-02
        hbb(16, 8, 1) = -0.15600d-02
        hbb(16, 9, 1) = -0.13600d-02
        hbb(16,10, 1) = -0.89600d-02
        hbb(16,11, 1) = 0.10110d-01
        hbb(16,12, 1) = -0.11500d-02
        hbb(16,13, 1) = -0.29030d-01
        hbb(16,14, 1) = 0.11640d-01
        hbb(16,15, 1) = 0.17390d-01
        hbb(16,16, 1) = 0.77110d-01
        hbb(17, 1, 1) = 0.22780d-01
        hbb(17, 2, 1) = -0.48100d-02
        hbb(17, 3, 1) = -0.17970d-01
        hbb(17, 4, 1) = -0.38000d-02
        hbb(17, 5, 1) = 0.60900d-02
        hbb(17, 6, 1) = -0.22900d-02
        hbb(17, 7, 1) = -0.14200d-02
        hbb(17, 8, 1) = 0.93000d-03
        hbb(17, 9, 1) = 0.49000d-03
        hbb(17,10, 1) = 0.10430d-01
        hbb(17,11, 1) = -0.67600d-02
        hbb(17,12, 1) = -0.36800d-02
        hbb(17,13, 1) = 0.11310d-01
        hbb(17,14, 1) = -0.50800d-02
        hbb(17,15, 1) = -0.62300d-02
        hbb(17,16, 1) = -0.39310d-01
        hbb(17,17, 1) = 0.76690d-01
        hbb(18, 1, 1) = 0.14500d-01
        hbb(18, 2, 1) = -0.69300d-02
        hbb(18, 3, 1) = -0.75700d-02
        hbb(18, 4, 1) = 0.85400d-02
        hbb(18, 5, 1) = -0.25100d-02
        hbb(18, 6, 1) = -0.60300d-02
        hbb(18, 7, 1) = -0.15000d-02
        hbb(18, 8, 1) = 0.63000d-03
        hbb(18, 9, 1) = 0.87000d-03
        hbb(18,10, 1) = -0.14700d-02
        hbb(18,11, 1) = -0.33600d-02
        hbb(18,12, 1) = 0.48300d-02
        hbb(18,13, 1) = 0.17720d-01
        hbb(18,14, 1) = -0.65600d-02
        hbb(18,15, 1) = -0.11160d-01
        hbb(18,16, 1) = -0.37800d-01
        hbb(18,17, 1) = -0.37380d-01
        hbb(18,18, 1) = 0.75170d-01
        ! Hessian elements for anchor point 1
        htt( 1, 1, 1) = 0.12410d-01
        htt( 2, 1, 1) = 0.66300d-02
        htt( 2, 2, 1) = 0.11350d-01
        htt( 3, 1, 1) = -0.60000d-02
        htt( 3, 2, 1) = -0.92200d-02
        htt( 3, 3, 1) = 0.16800d-01
        htt( 4, 1, 1) = -0.11780d-01
        htt( 4, 2, 1) = -0.45100d-02
        htt( 4, 3, 1) = 0.13580d-01
        htt( 4, 4, 1) = 0.20850d-01
        htt( 5, 1, 1) = -0.10450d-01
        htt( 5, 2, 1) = -0.82200d-02
        htt( 5, 3, 1) = 0.95200d-02
        htt( 5, 4, 1) = 0.11750d-01
        htt( 5, 5, 1) = 0.13600d-01
        htt( 6, 1, 1) = -0.73200d-02
        htt( 6, 2, 1) = -0.56100d-02
        htt( 6, 3, 1) = 0.79400d-02
        htt( 6, 4, 1) = 0.96500d-02
        htt( 6, 5, 1) = 0.70500d-02
        htt( 6, 6, 1) = 0.10410d-01
        htt( 7, 1, 1) = 0.68700d-02
        htt( 7, 2, 1) = 0.66900d-02
        htt( 7, 3, 1) = -0.11920d-01
        htt( 7, 4, 1) = -0.12100d-01
        htt( 7, 5, 1) = -0.51800d-02
        htt( 7, 6, 1) = -0.73000d-02
        htt( 7, 7, 1) = 0.12490d-01
        htt( 8, 1, 1) = 0.10000d-01
        htt( 8, 2, 1) = 0.93000d-02
        htt( 8, 3, 1) = -0.13500d-01
        htt( 8, 4, 1) = -0.14200d-01
        htt( 8, 5, 1) = -0.11730d-01
        htt( 8, 6, 1) = -0.39400d-02
        htt( 8, 7, 1) = 0.10370d-01
        htt( 8, 8, 1) = 0.18160d-01
        htt( 9, 1, 1) = -0.62000d-02
        htt( 9, 2, 1) = -0.56000d-03
        htt( 9, 3, 1) = -0.18500d-02
        htt( 9, 4, 1) = 0.37900d-02
        htt( 9, 5, 1) = 0.12100d-02
        htt( 9, 6, 1) = 0.24300d-02
        htt( 9, 7, 1) = -0.28800d-02
        htt( 9, 8, 1) = -0.16600d-02
        htt( 9, 9, 1) = 0.82000d-02
        htt(10, 1, 1) = -0.49900d-02
        htt(10, 2, 1) = -0.30000d-04
        htt(10, 3, 1) = 0.60000d-03
        htt(10, 4, 1) = 0.55600d-02
        htt(10, 5, 1) = 0.33700d-02
        htt(10, 6, 1) = 0.24600d-02
        htt(10, 7, 1) = -0.18800d-02
        htt(10, 8, 1) = -0.27900d-02
        htt(10, 9, 1) = 0.31400d-02
        htt(10,10, 1) = 0.89300d-02
        htt(11, 1, 1) = -0.90000d-04
        htt(11, 2, 1) = -0.55400d-02
        htt(11, 3, 1) = 0.15600d-02
        htt(11, 4, 1) = -0.39000d-02
        htt(11, 5, 1) = -0.11400d-02
        htt(11, 6, 1) = 0.63000d-03
        htt(11, 7, 1) = -0.26900d-02
        htt(11, 8, 1) = -0.91000d-03
        htt(11, 9, 1) = 0.22400d-02
        htt(11,10, 1) = -0.21000d-02
        htt(11,11, 1) = 0.80000d-02
        htt(12, 1, 1) = 0.11300d-02
        htt(12, 2, 1) = -0.50100d-02
        htt(12, 3, 1) = 0.40000d-02
        htt(12, 4, 1) = -0.21300d-02
        htt(12, 5, 1) = 0.10100d-02
        htt(12, 6, 1) = 0.67000d-03
        htt(12, 7, 1) = -0.17000d-02
        htt(12, 8, 1) = -0.20400d-02
        htt(12, 9, 1) = -0.28200d-02
        htt(12,10, 1) = 0.36900d-02
        htt(12,11, 1) = 0.36700d-02
        htt(12,12, 1) = 0.10180d-01
        htt(13, 1, 1) = -0.19600d-02
        htt(13, 2, 1) = -0.40500d-02
        htt(13, 3, 1) = 0.64000d-02
        htt(13, 4, 1) = 0.43100d-02
        htt(13, 5, 1) = 0.50600d-02
        htt(13, 6, 1) = 0.27500d-02
        htt(13, 7, 1) = -0.28000d-02
        htt(13, 8, 1) = -0.51000d-02
        htt(13, 9, 1) = -0.54600d-02
        htt(13,10, 1) = 0.27000d-03
        htt(13,11, 1) = -0.32400d-02
        htt(13,12, 1) = 0.24800d-02
        htt(13,13, 1) = 0.10020d-01
        htt(14, 1, 1) = 0.73000d-03
        htt(14, 2, 1) = -0.74000d-03
        htt(14, 3, 1) = 0.28700d-02
        htt(14, 4, 1) = 0.14100d-02
        htt(14, 5, 1) = 0.19900d-02
        htt(14, 6, 1) = -0.24000d-03
        htt(14, 7, 1) = -0.30000d-04
        htt(14, 8, 1) = -0.22600d-02
        htt(14, 9, 1) = -0.28600d-02
        htt(14,10, 1) = 0.17000d-02
        htt(14,11, 1) = -0.13100d-02
        htt(14,12, 1) = 0.32500d-02
        htt(14,13, 1) = 0.23900d-02
        htt(14,14, 1) = 0.74100d-02
        htt(15, 1, 1) = -0.31700d-02
        htt(15, 2, 1) = -0.45800d-02
        htt(15, 3, 1) = 0.39600d-02
        htt(15, 4, 1) = 0.25500d-02
        htt(15, 5, 1) = 0.29100d-02
        htt(15, 6, 1) = 0.27200d-02
        htt(15, 7, 1) = -0.37900d-02
        htt(15, 8, 1) = -0.39800d-02
        htt(15, 9, 1) = -0.41000d-03
        htt(15,10, 1) = -0.55000d-02
        htt(15,11, 1) = 0.10800d-02
        htt(15,12, 1) = -0.40100d-02
        htt(15,13, 1) = 0.43100d-02
        htt(15,14, 1) = -0.21600d-02
        htt(15,15, 1) = 0.93900d-02
        htt(16, 1, 1) = -0.48000d-03
        htt(16, 2, 1) = -0.12700d-02
        htt(16, 3, 1) = 0.43000d-03
        htt(16, 4, 1) = -0.35000d-03
        htt(16, 5, 1) = -0.16000d-03
        htt(16, 6, 1) = -0.28000d-03
        htt(16, 7, 1) = -0.10200d-02
        htt(16, 8, 1) = -0.11300d-02
        htt(16, 9, 1) = 0.21900d-02
        htt(16,10, 1) = -0.40700d-02
        htt(16,11, 1) = 0.30200d-02
        htt(16,12, 1) = -0.32400d-02
        htt(16,13, 1) = -0.33300d-02
        htt(16,14, 1) = 0.28700d-02
        htt(16,15, 1) = 0.29200d-02
        htt(16,16, 1) = 0.91100d-02
        htt(17, 1, 1) = 0.40200d-02
        htt(17, 2, 1) = 0.24400d-02
        htt(17, 3, 1) = -0.28100d-02
        htt(17, 4, 1) = -0.43900d-02
        htt(17, 5, 1) = -0.18700d-02
        htt(17, 6, 1) = -0.30500d-02
        htt(17, 7, 1) = 0.45500d-02
        htt(17, 8, 1) = 0.33700d-02
        htt(17, 9, 1) = 0.30000d-03
        htt(17,10, 1) = -0.19500d-02
        htt(17,11, 1) = 0.19700d-02
        htt(17,12, 1) = -0.28000d-03
        htt(17,13, 1) = -0.68200d-02
        htt(17,14, 1) = 0.40000d-03
        htt(17,15, 1) = -0.45700d-02
        htt(17,16, 1) = 0.26500d-02
        htt(17,17, 1) = 0.90000d-02
        htt(18, 1, 1) = 0.19000d-02
        htt(18, 2, 1) = 0.26000d-02
        htt(18, 3, 1) = -0.39700d-02
        htt(18, 4, 1) = -0.32700d-02
        htt(18, 5, 1) = -0.30200d-02
        htt(18, 6, 1) = -0.37000d-02
        htt(18, 7, 1) = 0.25000d-02
        htt(18, 8, 1) = 0.18200d-02
        htt(18, 9, 1) = 0.19300d-02
        htt(18,10, 1) = 0.43000d-03
        htt(18,11, 1) = 0.11900d-02
        htt(18,12, 1) = -0.31000d-03
        htt(18,13, 1) = -0.47600d-02
        htt(18,14, 1) = 0.15300d-02
        htt(18,15, 1) = -0.32600d-02
        htt(18,16, 1) = 0.30200d-02
        htt(18,17, 1) = 0.36100d-02
        htt(18,18, 1) = 0.88300d-02
        htt(19, 1, 1) = 0.13400d-02
        htt(19, 2, 1) = -0.86000d-03
        htt(19, 3, 1) = 0.71000d-03
        htt(19, 4, 1) = -0.14900d-02
        htt(19, 5, 1) = 0.11900d-02
        htt(19, 6, 1) = -0.60000d-04
        htt(19, 7, 1) = 0.17900d-02
        htt(19, 8, 1) = 0.54000d-03
        htt(19, 9, 1) = -0.23000d-02
        htt(19,10, 1) = -0.33800d-02
        htt(19,11, 1) = 0.40000d-04
        htt(19,12, 1) = -0.10500d-02
        htt(19,13, 1) = 0.81000d-03
        htt(19,14, 1) = -0.46200d-02
        htt(19,15, 1) = 0.18900d-02
        htt(19,16, 1) = -0.35300d-02
        htt(19,17, 1) = 0.17900d-02
        htt(19,18, 1) = -0.26700d-02
        htt(19,19, 1) = 0.72000d-02
        htt(20, 1, 1) = -0.78000d-03
        htt(20, 2, 1) = -0.70000d-03
        htt(20, 3, 1) = -0.45000d-03
        htt(20, 4, 1) = -0.37000d-03
        htt(20, 5, 1) = 0.40000d-04
        htt(20, 6, 1) = -0.71000d-03
        htt(20, 7, 1) = -0.27000d-03
        htt(20, 8, 1) = -0.10200d-02
        htt(20, 9, 1) = -0.66000d-03
        htt(20,10, 1) = -0.10000d-02
        htt(20,11, 1) = -0.74000d-03
        htt(20,12, 1) = -0.10700d-02
        htt(20,13, 1) = 0.28600d-02
        htt(20,14, 1) = -0.34900d-02
        htt(20,15, 1) = 0.32000d-02
        htt(20,16, 1) = -0.31600d-02
        htt(20,17, 1) = -0.36000d-02
        htt(20,18, 1) = 0.25500d-02
        htt(20,19, 1) = 0.27400d-02
        htt(20,20, 1) = 0.89000d-02
        htt(21, 1, 1) = 0.22500d-02
        htt(21, 2, 1) = 0.36500d-02
        htt(21, 3, 1) = -0.50600d-02
        htt(21, 4, 1) = -0.36600d-02
        htt(21, 5, 1) = -0.73600d-02
        htt(21, 6, 1) = -0.18700d-02
        htt(21, 7, 1) = -0.49000d-03
        htt(21, 8, 1) = 0.50000d-02
        htt(21, 9, 1) = 0.16800d-02
        htt(21,10, 1) = 0.80000d-04
        htt(21,11, 1) = 0.20000d-03
        htt(21,12, 1) = -0.14000d-02
        htt(21,13, 1) = -0.68000d-03
        htt(21,14, 1) = -0.25100d-02
        htt(21,15, 1) = 0.92000d-03
        htt(21,16, 1) = -0.91000d-03
        htt(21,17, 1) = -0.45300d-02
        htt(21,18, 1) = 0.28000d-03
        htt(21,19, 1) = -0.27000d-02
        htt(21,20, 1) = 0.21000d-02
        htt(21,21, 1) = 0.83700d-02
        htt(22, 1, 1) = -0.10700d-02
        htt(22, 2, 1) = 0.89000d-03
        htt(22, 3, 1) = -0.33900d-02
        htt(22, 4, 1) = -0.14300d-02
        htt(22, 5, 1) = -0.43000d-03
        htt(22, 6, 1) = -0.54300d-02
        htt(22, 7, 1) = 0.17600d-02
        htt(22, 8, 1) = -0.32400d-02
        htt(22, 9, 1) = 0.39000d-03
        htt(22,10, 1) = 0.10300d-02
        htt(22,11, 1) = -0.16800d-02
        htt(22,12, 1) = -0.10400d-02
        htt(22,13, 1) = 0.17600d-02
        htt(22,14, 1) = -0.15000d-03
        htt(22,15, 1) = 0.11200d-02
        htt(22,16, 1) = -0.78000d-03
        htt(22,17, 1) = -0.32800d-02
        htt(22,18, 1) = 0.99000d-03
        htt(22,19, 1) = -0.13800d-02
        htt(22,20, 1) = 0.29000d-02
        htt(22,21, 1) = 0.25500d-02
        htt(22,22, 1) = 0.78500d-02
        htt(23, 1, 1) = 0.43900d-02
        htt(23, 2, 1) = 0.34900d-02
        htt(23, 3, 1) = -0.38900d-02
        htt(23, 4, 1) = -0.47900d-02
        htt(23, 5, 1) = -0.62100d-02
        htt(23, 6, 1) = -0.12200d-02
        htt(23, 7, 1) = 0.15800d-02
        htt(23, 8, 1) = 0.65700d-02
        htt(23, 9, 1) = 0.40000d-04
        htt(23,10, 1) = -0.23200d-02
        htt(23,11, 1) = 0.98000d-03
        htt(23,12, 1) = -0.13800d-02
        htt(23,13, 1) = -0.27500d-02
        htt(23,14, 1) = -0.36400d-02
        htt(23,15, 1) = -0.39000d-03
        htt(23,16, 1) = -0.12800d-02
        htt(23,17, 1) = 0.91000d-03
        htt(23,18, 1) = -0.49900d-02
        htt(23,19, 1) = 0.18000d-02
        htt(23,20, 1) = -0.41000d-02
        htt(23,21, 1) = 0.35200d-02
        htt(23,22, 1) = -0.17600d-02
        htt(23,23, 1) = 0.94700d-02
        htt(24, 1, 1) = 0.10700d-02
        htt(24, 2, 1) = 0.73000d-03
        htt(24, 3, 1) = -0.22200d-02
        htt(24, 4, 1) = -0.25600d-02
        htt(24, 5, 1) = 0.73000d-03
        htt(24, 6, 1) = -0.47800d-02
        htt(24, 7, 1) = 0.38300d-02
        htt(24, 8, 1) = -0.16800d-02
        htt(24, 9, 1) = -0.12500d-02
        htt(24,10, 1) = -0.13700d-02
        htt(24,11, 1) = -0.90000d-03
        htt(24,12, 1) = -0.10100d-02
        htt(24,13, 1) = -0.31000d-03
        htt(24,14, 1) = -0.12800d-02
        htt(24,15, 1) = -0.20000d-03
        htt(24,16, 1) = -0.11600d-02
        htt(24,17, 1) = 0.21600d-02
        htt(24,18, 1) = -0.42700d-02
        htt(24,19, 1) = 0.31200d-02
        htt(24,20, 1) = -0.33100d-02
        htt(24,21, 1) = -0.23000d-02
        htt(24,22, 1) = 0.35300d-02
        htt(24,23, 1) = 0.41900d-02
        htt(24,24, 1) = 0.10010d-01
        ! Hessian elements for anchor point 1
        hbs( 1, 1, 1) = -0.26190d-01
        hbs( 1, 2, 1) = -0.29060d-01
        hbs( 1, 3, 1) = -0.39330d-01
        hbs( 1, 4, 1) = -0.78300d-02
        hbs( 1, 5, 1) = 0.21800d-02
        hbs( 1, 6, 1) = 0.40640d-01
        hbs( 1, 7, 1) = -0.41000d-03
        hbs( 1, 8, 1) = 0.40740d-01
        hbs( 1, 9, 1) = -0.20000d-04
        hbs( 1,10, 1) = -0.87600d-02
        hbs( 1,11, 1) = -0.35000d-03
        hbs( 1,12, 1) = 0.45000d-02
        hbs( 2, 1, 1) = 0.36740d-01
        hbs( 2, 2, 1) = -0.10920d-01
        hbs( 2, 3, 1) = 0.13650d-01
        hbs( 2, 4, 1) = 0.34200d-02
        hbs( 2, 5, 1) = 0.62000d-03
        hbs( 2, 6, 1) = -0.29580d-01
        hbs( 2, 7, 1) = 0.13300d-02
        hbs( 2, 8, 1) = -0.12060d-01
        hbs( 2, 9, 1) = 0.19000d-03
        hbs( 2,10, 1) = 0.40800d-02
        hbs( 2,11, 1) = -0.39000d-03
        hbs( 2,12, 1) = -0.22000d-03
        hbs( 3, 1, 1) = -0.10550d-01
        hbs( 3, 2, 1) = 0.39980d-01
        hbs( 3, 3, 1) = 0.25680d-01
        hbs( 3, 4, 1) = 0.44100d-02
        hbs( 3, 5, 1) = -0.28000d-02
        hbs( 3, 6, 1) = -0.11060d-01
        hbs( 3, 7, 1) = -0.91000d-03
        hbs( 3, 8, 1) = -0.28680d-01
        hbs( 3, 9, 1) = -0.17000d-03
        hbs( 3,10, 1) = 0.46800d-02
        hbs( 3,11, 1) = 0.74000d-03
        hbs( 3,12, 1) = -0.42800d-02
        hbs( 4, 1, 1) = -0.26180d-01
        hbs( 4, 2, 1) = -0.40500d-02
        hbs( 4, 3, 1) = 0.14760d-01
        hbs( 4, 4, 1) = -0.35070d-01
        hbs( 4, 5, 1) = -0.56800d-02
        hbs( 4, 6, 1) = -0.74200d-02
        hbs( 4, 7, 1) = 0.38500d-02
        hbs( 4, 8, 1) = 0.33760d-01
        hbs( 4, 9, 1) = -0.72000d-03
        hbs( 4,10, 1) = 0.51570d-01
        hbs( 4,11, 1) = 0.20000d-04
        hbs( 4,12, 1) = -0.81000d-03
        hbs( 5, 1, 1) = 0.28180d-01
        hbs( 5, 2, 1) = 0.24300d-02
        hbs( 5, 3, 1) = -0.87400d-02
        hbs( 5, 4, 1) = 0.50900d-02
        hbs( 5, 5, 1) = 0.11940d-01
        hbs( 5, 6, 1) = 0.32000d-02
        hbs( 5, 7, 1) = -0.88000d-03
        hbs( 5, 8, 1) = -0.14770d-01
        hbs( 5, 9, 1) = 0.30000d-04
        hbs( 5,10, 1) = -0.27450d-01
        hbs( 5,11, 1) = 0.00000d+00
        hbs( 5,12, 1) = 0.11400d-02
        hbs( 6, 1, 1) = -0.20100d-02
        hbs( 6, 2, 1) = 0.16200d-02
        hbs( 6, 3, 1) = -0.60300d-02
        hbs( 6, 4, 1) = 0.29980d-01
        hbs( 6, 5, 1) = -0.62600d-02
        hbs( 6, 6, 1) = 0.42200d-02
        hbs( 6, 7, 1) = -0.29700d-02
        hbs( 6, 8, 1) = -0.18990d-01
        hbs( 6, 9, 1) = 0.69000d-03
        hbs( 6,10, 1) = -0.24120d-01
        hbs( 6,11, 1) = -0.10000d-04
        hbs( 6,12, 1) = -0.33000d-03
        hbs( 7, 1, 1) = -0.59200d-02
        hbs( 7, 2, 1) = 0.29370d-01
        hbs( 7, 3, 1) = 0.59700d-02
        hbs( 7, 4, 1) = -0.34310d-01
        hbs( 7, 5, 1) = 0.35700d-02
        hbs( 7, 6, 1) = -0.38110d-01
        hbs( 7, 7, 1) = -0.77800d-02
        hbs( 7, 8, 1) = -0.29600d-02
        hbs( 7, 9, 1) = 0.44000d-02
        hbs( 7,10, 1) = 0.44330d-01
        hbs( 7,11, 1) = -0.52000d-03
        hbs( 7,12, 1) = -0.41000d-03
        hbs( 8, 1, 1) = 0.27500d-02
        hbs( 8, 2, 1) = -0.16290d-01
        hbs( 8, 3, 1) = -0.20500d-02
        hbs( 8, 4, 1) = 0.31190d-01
        hbs( 8, 5, 1) = -0.23800d-02
        hbs( 8, 6, 1) = 0.45900d-02
        hbs( 8, 7, 1) = 0.44200d-02
        hbs( 8, 8, 1) = 0.12900d-02
        hbs( 8, 9, 1) = -0.16900d-02
        hbs( 8,10, 1) = -0.19840d-01
        hbs( 8,11, 1) = -0.28000d-03
        hbs( 8,12, 1) = 0.25000d-03
        hbs( 9, 1, 1) = 0.31600d-02
        hbs( 9, 2, 1) = -0.13080d-01
        hbs( 9, 3, 1) = -0.39200d-02
        hbs( 9, 4, 1) = 0.31100d-02
        hbs( 9, 5, 1) = -0.11900d-02
        hbs( 9, 6, 1) = 0.33520d-01
        hbs( 9, 7, 1) = 0.33700d-02
        hbs( 9, 8, 1) = 0.16700d-02
        hbs( 9, 9, 1) = -0.27000d-02
        hbs( 9,10, 1) = -0.24490d-01
        hbs( 9,11, 1) = 0.80000d-03
        hbs( 9,12, 1) = 0.17000d-03
        hbs(10, 1, 1) = 0.35610d-01
        hbs(10, 2, 1) = 0.35190d-01
        hbs(10, 3, 1) = -0.80000d-03
        hbs(10, 4, 1) = -0.15230d-01
        hbs(10, 5, 1) = 0.40000d-04
        hbs(10, 6, 1) = -0.28120d-01
        hbs(10, 7, 1) = 0.45800d-02
        hbs(10, 8, 1) = -0.27160d-01
        hbs(10, 9, 1) = -0.74500d-02
        hbs(10,10, 1) = -0.14900d-01
        hbs(10,11, 1) = 0.44000d-02
        hbs(10,12, 1) = -0.51000d-03
        hbs(11, 1, 1) = -0.19830d-01
        hbs(11, 2, 1) = -0.15670d-01
        hbs(11, 3, 1) = 0.46000d-03
        hbs(11, 4, 1) = 0.87500d-02
        hbs(11, 5, 1) = 0.12000d-03
        hbs(11, 6, 1) = 0.26580d-01
        hbs(11, 7, 1) = -0.28700d-02
        hbs(11, 8, 1) = 0.69000d-03
        hbs(11, 9, 1) = 0.35900d-02
        hbs(11,10, 1) = 0.63600d-02
        hbs(11,11, 1) = -0.16100d-02
        hbs(11,12, 1) = -0.22000d-03
        hbs(12, 1, 1) = -0.15780d-01
        hbs(12, 2, 1) = -0.19520d-01
        hbs(12, 3, 1) = 0.34000d-03
        hbs(12, 4, 1) = 0.64800d-02
        hbs(12, 5, 1) = -0.15000d-03
        hbs(12, 6, 1) = 0.15400d-02
        hbs(12, 7, 1) = -0.17100d-02
        hbs(12, 8, 1) = 0.26480d-01
        hbs(12, 9, 1) = 0.38600d-02
        hbs(12,10, 1) = 0.85400d-02
        hbs(12,11, 1) = -0.27900d-02
        hbs(12,12, 1) = 0.73000d-03
        hbs(13, 1, 1) = 0.28450d-01
        hbs(13, 2, 1) = -0.45800d-02
        hbs(13, 3, 1) = 0.32200d-02
        hbs(13, 4, 1) = 0.42370d-01
        hbs(13, 5, 1) = -0.53000d-03
        hbs(13, 6, 1) = -0.11000d-02
        hbs(13, 7, 1) = -0.47000d-03
        hbs(13, 8, 1) = -0.36660d-01
        hbs(13, 9, 1) = 0.43000d-02
        hbs(13,10, 1) = -0.34790d-01
        hbs(13,11, 1) = -0.78500d-02
        hbs(13,12, 1) = 0.43000d-02
        hbs(14, 1, 1) = -0.12720d-01
        hbs(14, 2, 1) = 0.26100d-02
        hbs(14, 3, 1) = -0.28000d-02
        hbs(14, 4, 1) = -0.23620d-01
        hbs(14, 5, 1) = 0.19000d-03
        hbs(14, 6, 1) = 0.10100d-02
        hbs(14, 7, 1) = 0.75000d-03
        hbs(14, 8, 1) = 0.32490d-01
        hbs(14, 9, 1) = -0.26900d-02
        hbs(14,10, 1) = 0.34600d-02
        hbs(14,11, 1) = 0.37900d-02
        hbs(14,12, 1) = -0.16300d-02
        hbs(15, 1, 1) = -0.15730d-01
        hbs(15, 2, 1) = 0.19700d-02
        hbs(15, 3, 1) = -0.42000d-03
        hbs(15, 4, 1) = -0.18750d-01
        hbs(15, 5, 1) = 0.33000d-03
        hbs(15, 6, 1) = 0.90000d-04
        hbs(15, 7, 1) = -0.28000d-03
        hbs(15, 8, 1) = 0.41800d-02
        hbs(15, 9, 1) = -0.16100d-02
        hbs(15,10, 1) = 0.31330d-01
        hbs(15,11, 1) = 0.40600d-02
        hbs(15,12, 1) = -0.26700d-02
        hbs(16, 1, 1) = -0.57800d-02
        hbs(16, 2, 1) = -0.26870d-01
        hbs(16, 3, 1) = 0.16180d-01
        hbs(16, 4, 1) = 0.50070d-01
        hbs(16, 5, 1) = 0.42000d-03
        hbs(16, 6, 1) = 0.34120d-01
        hbs(16, 7, 1) = 0.25000d-03
        hbs(16, 8, 1) = -0.77100d-02
        hbs(16, 9, 1) = -0.51000d-03
        hbs(16,10, 1) = -0.37450d-01
        hbs(16,11, 1) = 0.43000d-02
        hbs(16,12, 1) = -0.70600d-02
        hbs(17, 1, 1) = 0.20300d-02
        hbs(17, 2, 1) = 0.26580d-01
        hbs(17, 3, 1) = -0.98100d-02
        hbs(17, 4, 1) = -0.26570d-01
        hbs(17, 5, 1) = 0.45000d-03
        hbs(17, 6, 1) = -0.15140d-01
        hbs(17, 7, 1) = -0.15000d-03
        hbs(17, 8, 1) = 0.35400d-02
        hbs(17, 9, 1) = -0.25000d-03
        hbs(17,10, 1) = 0.60500d-02
        hbs(17,11, 1) = -0.15900d-02
        hbs(17,12, 1) = 0.61300d-02
        hbs(18, 1, 1) = 0.37500d-02
        hbs(18, 2, 1) = 0.29000d-03
        hbs(18, 3, 1) = -0.63700d-02
        hbs(18, 4, 1) = -0.23500d-01
        hbs(18, 5, 1) = -0.87000d-03
        hbs(18, 6, 1) = -0.18980d-01
        hbs(18, 7, 1) = -0.10000d-03
        hbs(18, 8, 1) = 0.41700d-02
        hbs(18, 9, 1) = 0.76000d-03
        hbs(18,10, 1) = 0.31400d-01
        hbs(18,11, 1) = -0.27100d-02
        hbs(18,12, 1) = 0.93000d-03
        ! Hessian elements for anchor point 2
        hss( 1, 1, 2) = 0.30793d+00
        hss( 2, 1, 2) = 0.25300d-02
        hss( 2, 2, 2) = 0.31227d+00
        hss( 3, 1, 2) = 0.72500d-01
        hss( 3, 2, 2) = 0.70950d-01
        hss( 3, 3, 2) = 0.63382d+00
        hss( 4, 1, 2) = 0.56300d-02
        hss( 4, 2, 2) = 0.20240d-01
        hss( 4, 3, 2) = -0.11340d-01
        hss( 4, 4, 2) = 0.43881d+00
        hss( 5, 1, 2) = 0.31200d-02
        hss( 5, 2, 2) = -0.58000d-03
        hss( 5, 3, 2) = -0.14100d-02
        hss( 5, 4, 2) = 0.19700d-02
        hss( 5, 5, 2) = 0.36903d+00
        hss( 6, 1, 2) = 0.93200d-02
        hss( 6, 2, 2) = 0.48920d-01
        hss( 6, 3, 2) = 0.11310d-01
        hss( 6, 4, 2) = 0.43540d-01
        hss( 6, 5, 2) = -0.31000d-03
        hss( 6, 6, 2) = 0.34629d+00
        hss( 7, 1, 2) = -0.79000d-03
        hss( 7, 2, 2) = -0.53000d-03
        hss( 7, 3, 2) = 0.11700d-02
        hss( 7, 4, 2) = 0.55200d-02
        hss( 7, 5, 2) = 0.15000d-03
        hss( 7, 6, 2) = 0.35500d-02
        hss( 7, 7, 2) = 0.36483d+00
        hss( 8, 1, 2) = 0.48320d-01
        hss( 8, 2, 2) = 0.88800d-02
        hss( 8, 3, 2) = 0.10960d-01
        hss( 8, 4, 2) = -0.62700d-02
        hss( 8, 5, 2) = -0.47000d-03
        hss( 8, 6, 2) = 0.30770d-01
        hss( 8, 7, 2) = -0.33000d-03
        hss( 8, 8, 2) = 0.34600d+00
        hss( 9, 1, 2) = -0.27000d-03
        hss( 9, 2, 2) = -0.90000d-04
        hss( 9, 3, 2) = -0.37000d-03
        hss( 9, 4, 2) = -0.49000d-03
        hss( 9, 5, 2) = -0.90000d-04
        hss( 9, 6, 2) = 0.39700d-02
        hss( 9, 7, 2) = 0.32000d-03
        hss( 9, 8, 2) = 0.40100d-02
        hss( 9, 9, 2) = 0.36524d+00
        hss(10, 1, 2) = 0.17870d-01
        hss(10, 2, 2) = 0.45600d-02
        hss(10, 3, 2) = -0.11880d-01
        hss(10, 4, 2) = 0.81470d-01
        hss(10, 5, 2) = -0.59000d-03
        hss(10, 6, 2) = -0.42200d-02
        hss(10, 7, 2) = -0.61000d-03
        hss(10, 8, 2) = 0.41780d-01
        hss(10, 9, 2) = -0.69000d-03
        hss(10,10, 2) = 0.44702d+00
        hss(11, 1, 2) = -0.59000d-03
        hss(11, 2, 2) = -0.11700d-02
        hss(11, 3, 2) = 0.11000d-02
        hss(11, 4, 2) = -0.48000d-03
        hss(11, 5, 2) = 0.30000d-04
        hss(11, 6, 2) = -0.44000d-03
        hss(11, 7, 2) = 0.16000d-03
        hss(11, 8, 2) = 0.35200d-02
        hss(11, 9, 2) = 0.40000d-03
        hss(11,10, 2) = 0.50000d-02
        hss(11,11, 2) = 0.36451d+00
        hss(12, 1, 2) = 0.14000d-03
        hss(12, 2, 2) = 0.21900d-02
        hss(12, 3, 2) = 0.30000d-04
        hss(12, 4, 2) = -0.45000d-03
        hss(12, 5, 2) = 0.10000d-04
        hss(12, 6, 2) = -0.43000d-03
        hss(12, 7, 2) = 0.60000d-04
        hss(12, 8, 2) = 0.80000d-04
        hss(12, 9, 2) = 0.30000d-04
        hss(12,10, 2) = 0.27300d-02
        hss(12,11, 2) = 0.33000d-03
        hss(12,12, 2) = 0.36751d+00
        ! Hessian elements for anchor point 2
        hbb( 1, 1, 2) = 0.81080d-01
        hbb( 2, 1, 2) = -0.44810d-01
        hbb( 2, 2, 2) = 0.14731d+00
        hbb( 3, 1, 2) = -0.36270d-01
        hbb( 3, 2, 2) = -0.10250d+00
        hbb( 3, 3, 2) = 0.13877d+00
        hbb( 4, 1, 2) = -0.36710d-01
        hbb( 4, 2, 2) = 0.25500d-01
        hbb( 4, 3, 2) = 0.11210d-01
        hbb( 4, 4, 2) = 0.75960d-01
        hbb( 5, 1, 2) = 0.22750d-01
        hbb( 5, 2, 2) = -0.18060d-01
        hbb( 5, 3, 2) = -0.46900d-02
        hbb( 5, 4, 2) = -0.38310d-01
        hbb( 5, 5, 2) = 0.72880d-01
        hbb( 6, 1, 2) = 0.13960d-01
        hbb( 6, 2, 2) = -0.74300d-02
        hbb( 6, 3, 2) = -0.65200d-02
        hbb( 6, 4, 2) = -0.37640d-01
        hbb( 6, 5, 2) = -0.34570d-01
        hbb( 6, 6, 2) = 0.72210d-01
        hbb( 7, 1, 2) = -0.88100d-02
        hbb( 7, 2, 2) = -0.49400d-02
        hbb( 7, 3, 2) = 0.13750d-01
        hbb( 7, 4, 2) = -0.28850d-01
        hbb( 7, 5, 2) = 0.11260d-01
        hbb( 7, 6, 2) = 0.17590d-01
        hbb( 7, 7, 2) = 0.74530d-01
        hbb( 8, 1, 2) = -0.15100d-02
        hbb( 8, 2, 2) = 0.76700d-02
        hbb( 8, 3, 2) = -0.61600d-02
        hbb( 8, 4, 2) = 0.17480d-01
        hbb( 8, 5, 2) = -0.61800d-02
        hbb( 8, 6, 2) = -0.11300d-01
        hbb( 8, 7, 2) = -0.36990d-01
        hbb( 8, 8, 2) = 0.78710d-01
        hbb( 9, 1, 2) = 0.10320d-01
        hbb( 9, 2, 2) = -0.27300d-02
        hbb( 9, 3, 2) = -0.75900d-02
        hbb( 9, 4, 2) = 0.11370d-01
        hbb( 9, 5, 2) = -0.50800d-02
        hbb( 9, 6, 2) = -0.62900d-02
        hbb( 9, 7, 2) = -0.37540d-01
        hbb( 9, 8, 2) = -0.41720d-01
        hbb( 9, 9, 2) = 0.79260d-01
        hbb(10, 1, 2) = 0.10400d-01
        hbb(10, 2, 2) = -0.60200d-02
        hbb(10, 3, 2) = -0.43800d-02
        hbb(10, 4, 2) = -0.90700d-02
        hbb(10, 5, 2) = 0.10230d-01
        hbb(10, 6, 2) = -0.11600d-02
        hbb(10, 7, 2) = -0.34090d-01
        hbb(10, 8, 2) = 0.13820d-01
        hbb(10, 9, 2) = 0.20260d-01
        hbb(10,10, 2) = 0.76010d-01
        hbb(11, 1, 2) = -0.53400d-02
        hbb(11, 2, 2) = 0.37500d-02
        hbb(11, 3, 2) = 0.15900d-02
        hbb(11, 4, 2) = -0.11000d-02
        hbb(11, 5, 2) = -0.35000d-02
        hbb(11, 6, 2) = 0.46000d-02
        hbb(11, 7, 2) = 0.20040d-01
        hbb(11, 8, 2) = -0.74900d-02
        hbb(11, 9, 2) = -0.12560d-01
        hbb(11,10, 2) = -0.37950d-01
        hbb(11,11, 2) = 0.79720d-01
        hbb(12, 1, 2) = -0.50600d-02
        hbb(12, 2, 2) = 0.22700d-02
        hbb(12, 3, 2) = 0.27900d-02
        hbb(12, 4, 2) = 0.10170d-01
        hbb(12, 5, 2) = -0.67300d-02
        hbb(12, 6, 2) = -0.34400d-02
        hbb(12, 7, 2) = 0.14040d-01
        hbb(12, 8, 2) = -0.63400d-02
        hbb(12, 9, 2) = -0.77100d-02
        hbb(12,10, 2) = -0.38060d-01
        hbb(12,11, 2) = -0.41770d-01
        hbb(12,12, 2) = 0.79830d-01
        hbb(13, 1, 2) = -0.81300d-02
        hbb(13, 2, 2) = 0.15080d-01
        hbb(13, 3, 2) = -0.69500d-02
        hbb(13, 4, 2) = 0.40600d-02
        hbb(13, 5, 2) = -0.19200d-02
        hbb(13, 6, 2) = -0.21400d-02
        hbb(13, 7, 2) = -0.67600d-02
        hbb(13, 8, 2) = 0.90200d-02
        hbb(13, 9, 2) = -0.22600d-02
        hbb(13,10, 2) = -0.33920d-01
        hbb(13,11, 2) = 0.13940d-01
        hbb(13,12, 2) = 0.19980d-01
        hbb(13,13, 2) = 0.73910d-01
        hbb(14, 1, 2) = 0.96700d-02
        hbb(14, 2, 2) = -0.84100d-02
        hbb(14, 3, 2) = -0.12600d-02
        hbb(14, 4, 2) = -0.20700d-02
        hbb(14, 5, 2) = 0.84000d-03
        hbb(14, 6, 2) = 0.12300d-02
        hbb(14, 7, 2) = -0.21800d-02
        hbb(14, 8, 2) = -0.30500d-02
        hbb(14, 9, 2) = 0.52300d-02
        hbb(14,10, 2) = 0.20350d-01
        hbb(14,11, 2) = -0.77300d-02
        hbb(14,12, 2) = -0.12620d-01
        hbb(14,13, 2) = -0.37510d-01
        hbb(14,14, 2) = 0.78920d-01
        hbb(15, 1, 2) = -0.15300d-02
        hbb(15, 2, 2) = -0.66700d-02
        hbb(15, 3, 2) = 0.82100d-02
        hbb(15, 4, 2) = -0.20000d-02
        hbb(15, 5, 2) = 0.10900d-02
        hbb(15, 6, 2) = 0.91000d-03
        hbb(15, 7, 2) = 0.89400d-02
        hbb(15, 8, 2) = -0.59700d-02
        hbb(15, 9, 2) = -0.29700d-02
        hbb(15,10, 2) = 0.13570d-01
        hbb(15,11, 2) = -0.62100d-02
        hbb(15,12, 2) = -0.73600d-02
        hbb(15,13, 2) = -0.36400d-01
        hbb(15,14, 2) = -0.41400d-01
        hbb(15,15, 2) = 0.77800d-01
        hbb(16, 1, 2) = -0.37830d-01
        hbb(16, 2, 2) = 0.15190d-01
        hbb(16, 3, 2) = 0.22640d-01
        hbb(16, 4, 2) = -0.53900d-02
        hbb(16, 5, 2) = -0.40000d-02
        hbb(16, 6, 2) = 0.93900d-02
        hbb(16, 7, 2) = 0.39700d-02
        hbb(16, 8, 2) = -0.18200d-02
        hbb(16, 9, 2) = -0.21500d-02
        hbb(16,10, 2) = -0.93400d-02
        hbb(16,11, 2) = 0.10420d-01
        hbb(16,12, 2) = -0.10800d-02
        hbb(16,13, 2) = -0.29170d-01
        hbb(16,14, 2) = 0.11750d-01
        hbb(16,15, 2) = 0.17420d-01
        hbb(16,16, 2) = 0.77760d-01
        hbb(17, 1, 2) = 0.23450d-01
        hbb(17, 2, 2) = -0.61200d-02
        hbb(17, 3, 2) = -0.17330d-01
        hbb(17, 4, 2) = -0.38800d-02
        hbb(17, 5, 2) = 0.65900d-02
        hbb(17, 6, 2) = -0.27000d-02
        hbb(17, 7, 2) = -0.18400d-02
        hbb(17, 8, 2) = 0.10400d-02
        hbb(17, 9, 2) = 0.79000d-03
        hbb(17,10, 2) = 0.10460d-01
        hbb(17,11, 2) = -0.67900d-02
        hbb(17,12, 2) = -0.36700d-02
        hbb(17,13, 2) = 0.11210d-01
        hbb(17,14, 2) = -0.50900d-02
        hbb(17,15, 2) = -0.61300d-02
        hbb(17,16, 2) = -0.39400d-01
        hbb(17,17, 2) = 0.75830d-01
        hbb(18, 1, 2) = 0.14380d-01
        hbb(18, 2, 2) = -0.90700d-02
        hbb(18, 3, 2) = -0.53200d-02
        hbb(18, 4, 2) = 0.92800d-02
        hbb(18, 5, 2) = -0.25900d-02
        hbb(18, 6, 2) = -0.66900d-02
        hbb(18, 7, 2) = -0.21300d-02
        hbb(18, 8, 2) = 0.77000d-03
        hbb(18, 9, 2) = 0.13500d-02
        hbb(18,10, 2) = -0.11200d-02
        hbb(18,11, 2) = -0.36200d-02
        hbb(18,12, 2) = 0.47500d-02
        hbb(18,13, 2) = 0.17950d-01
        hbb(18,14, 2) = -0.66600d-02
        hbb(18,15, 2) = -0.11290d-01
        hbb(18,16, 2) = -0.38360d-01
        hbb(18,17, 2) = -0.36420d-01
        hbb(18,18, 2) = 0.74790d-01
        ! Hessian elements for anchor point 2
        htt( 1, 1, 2) = 0.87500d-02
        htt( 2, 1, 2) = 0.39300d-02
        htt( 2, 2, 2) = 0.86800d-02
        htt( 3, 1, 2) = -0.25000d-03
        htt( 3, 2, 2) = -0.40300d-02
        htt( 3, 3, 2) = 0.69500d-02
        htt( 4, 1, 2) = -0.50700d-02
        htt( 4, 2, 2) = 0.71000d-03
        htt( 4, 3, 2) = 0.31700d-02
        htt( 4, 4, 2) = 0.89500d-02
        htt( 5, 1, 2) = -0.63800d-02
        htt( 5, 2, 2) = -0.48700d-02
        htt( 5, 3, 2) = 0.26100d-02
        htt( 5, 4, 2) = 0.41100d-02
        htt( 5, 5, 2) = 0.87900d-02
        htt( 6, 1, 2) = -0.41700d-02
        htt( 6, 2, 2) = -0.31900d-02
        htt( 6, 3, 2) = 0.27300d-02
        htt( 6, 4, 2) = 0.37100d-02
        htt( 6, 5, 2) = 0.34900d-02
        htt( 6, 6, 2) = 0.76600d-02
        htt( 7, 1, 2) = 0.23700d-02
        htt( 7, 2, 2) = 0.28700d-02
        htt( 7, 3, 2) = -0.43900d-02
        htt( 7, 4, 2) = -0.39000d-02
        htt( 7, 5, 2) = 0.60000d-04
        htt( 7, 6, 2) = -0.32200d-02
        htt( 7, 7, 2) = 0.66300d-02
        htt( 8, 1, 2) = 0.45800d-02
        htt( 8, 2, 2) = 0.45500d-02
        htt( 8, 3, 2) = -0.42700d-02
        htt( 8, 4, 2) = -0.43000d-02
        htt( 8, 5, 2) = -0.52500d-02
        htt( 8, 6, 2) = 0.95000d-03
        htt( 8, 7, 2) = 0.33500d-02
        htt( 8, 8, 2) = 0.95500d-02
        htt( 9, 1, 2) = -0.59200d-02
        htt( 9, 2, 2) = -0.90000d-03
        htt( 9, 3, 2) = -0.18700d-02
        htt( 9, 4, 2) = 0.31600d-02
        htt( 9, 5, 2) = 0.12000d-02
        htt( 9, 6, 2) = 0.22300d-02
        htt( 9, 7, 2) = -0.27400d-02
        htt( 9, 8, 2) = -0.17100d-02
        htt( 9, 9, 2) = 0.79900d-02
        htt(10, 1, 2) = -0.37700d-02
        htt(10, 2, 2) = 0.64000d-03
        htt(10, 3, 2) = -0.10100d-02
        htt(10, 4, 2) = 0.33900d-02
        htt(10, 5, 2) = 0.20800d-02
        htt(10, 6, 2) = 0.14400d-02
        htt(10, 7, 2) = -0.59000d-03
        htt(10, 8, 2) = -0.12400d-02
        htt(10, 9, 2) = 0.28400d-02
        htt(10,10, 2) = 0.84900d-02
        htt(11, 1, 2) = -0.84000d-03
        htt(11, 2, 2) = -0.59000d-02
        htt(11, 3, 2) = 0.21200d-02
        htt(11, 4, 2) = -0.29400d-02
        htt(11, 5, 2) = -0.38000d-03
        htt(11, 6, 2) = 0.11900d-02
        htt(11, 7, 2) = -0.32600d-02
        htt(11, 8, 2) = -0.16800d-02
        htt(11, 9, 2) = 0.26900d-02
        htt(11,10, 2) = -0.18100d-02
        htt(11,11, 2) = 0.80400d-02
        htt(12, 1, 2) = 0.13200d-02
        htt(12, 2, 2) = -0.43700d-02
        htt(12, 3, 2) = 0.29800d-02
        htt(12, 4, 2) = -0.27100d-02
        htt(12, 5, 2) = 0.50000d-03
        htt(12, 6, 2) = 0.41000d-03
        htt(12, 7, 2) = -0.11200d-02
        htt(12, 8, 2) = -0.12100d-02
        htt(12, 9, 2) = -0.24600d-02
        htt(12,10, 2) = 0.38400d-02
        htt(12,11, 2) = 0.35400d-02
        htt(12,12, 2) = 0.98300d-02
        htt(13, 1, 2) = 0.53000d-03
        htt(13, 2, 2) = -0.13800d-02
        htt(13, 3, 2) = 0.17500d-02
        htt(13, 4, 2) = -0.16000d-03
        htt(13, 5, 2) = 0.18900d-02
        htt(13, 6, 2) = 0.51000d-03
        htt(13, 7, 2) = 0.70000d-03
        htt(13, 8, 2) = -0.68000d-03
        htt(13, 9, 2) = -0.54000d-02
        htt(13,10, 2) = -0.17000d-03
        htt(13,11, 2) = -0.33800d-02
        htt(13,12, 2) = 0.18400d-02
        htt(13,13, 2) = 0.80400d-02
        htt(14, 1, 2) = 0.19700d-02
        htt(14, 2, 2) = 0.54000d-03
        htt(14, 3, 2) = 0.79000d-03
        htt(14, 4, 2) = -0.65000d-03
        htt(14, 5, 2) = 0.36000d-03
        htt(14, 6, 2) = -0.13400d-02
        htt(14, 7, 2) = 0.15100d-02
        htt(14, 8, 2) = -0.19000d-03
        htt(14, 9, 2) = -0.30400d-02
        htt(14,10, 2) = 0.14100d-02
        htt(14,11, 2) = -0.15300d-02
        htt(14,12, 2) = 0.29300d-02
        htt(14,13, 2) = 0.18300d-02
        htt(14,14, 2) = 0.71800d-02
        htt(15, 1, 2) = -0.16100d-02
        htt(15, 2, 2) = -0.29100d-02
        htt(15, 3, 2) = 0.90000d-03
        htt(15, 4, 2) = -0.39000d-03
        htt(15, 5, 2) = 0.10100d-02
        htt(15, 6, 2) = 0.13000d-02
        htt(15, 7, 2) = -0.14300d-02
        htt(15, 8, 2) = -0.11500d-02
        htt(15, 9, 2) = -0.28000d-03
        htt(15,10, 2) = -0.57800d-02
        htt(15,11, 2) = 0.10900d-02
        htt(15,12, 2) = -0.44100d-02
        htt(15,13, 2) = 0.28500d-02
        htt(15,14, 2) = -0.26000d-02
        htt(15,15, 2) = 0.83100d-02
        htt(16, 1, 2) = -0.16000d-03
        htt(16, 2, 2) = -0.99000d-03
        htt(16, 3, 2) = -0.50000d-04
        htt(16, 4, 2) = -0.88000d-03
        htt(16, 5, 2) = -0.52000d-03
        htt(16, 6, 2) = -0.55000d-03
        htt(16, 7, 2) = -0.62000d-03
        htt(16, 8, 2) = -0.66000d-03
        htt(16, 9, 2) = 0.20700d-02
        htt(16,10, 2) = -0.42000d-02
        htt(16,11, 2) = 0.29500d-02
        htt(16,12, 2) = -0.33300d-02
        htt(16,13, 2) = -0.33700d-02
        htt(16,14, 2) = 0.27500d-02
        htt(16,15, 2) = 0.28600d-02
        htt(16,16, 2) = 0.89800d-02
        htt(17, 1, 2) = 0.19500d-02
        htt(17, 2, 2) = 0.45000d-03
        htt(17, 3, 2) = 0.64000d-03
        htt(17, 4, 2) = -0.85000d-03
        htt(17, 5, 2) = 0.51000d-03
        htt(17, 6, 2) = -0.12300d-02
        htt(17, 7, 2) = 0.17800d-02
        htt(17, 8, 2) = 0.40000d-04
        htt(17, 9, 2) = 0.55000d-03
        htt(17,10, 2) = -0.15600d-02
        htt(17,11, 2) = 0.21200d-02
        htt(17,12, 2) = 0.10000d-04
        htt(17,13, 2) = -0.55500d-02
        htt(17,14, 2) = 0.56000d-03
        htt(17,15, 2) = -0.34500d-02
        htt(17,16, 2) = 0.26500d-02
        htt(17,17, 2) = 0.80100d-02
        htt(18, 1, 2) = 0.68000d-03
        htt(18, 2, 2) = 0.11500d-02
        htt(18, 3, 2) = -0.14100d-02
        htt(18, 4, 2) = -0.93000d-03
        htt(18, 5, 2) = -0.14500d-02
        htt(18, 6, 2) = -0.25900d-02
        htt(18, 7, 2) = 0.58000d-03
        htt(18, 8, 2) = -0.57000d-03
        htt(18, 9, 2) = 0.17600d-02
        htt(18,10, 2) = 0.54000d-03
        htt(18,11, 2) = 0.12600d-02
        htt(18,12, 2) = 0.40000d-04
        htt(18,13, 2) = -0.35300d-02
        htt(18,14, 2) = 0.18700d-02
        htt(18,15, 2) = -0.23200d-02
        htt(18,16, 2) = 0.30800d-02
        htt(18,17, 2) = 0.27300d-02
        htt(18,18, 2) = 0.81100d-02
        htt(19, 1, 2) = 0.50000d-03
        htt(19, 2, 2) = -0.14600d-02
        htt(19, 3, 2) = 0.16000d-02
        htt(19, 4, 2) = -0.37000d-03
        htt(19, 5, 2) = 0.20400d-02
        htt(19, 6, 2) = 0.62000d-03
        htt(19, 7, 2) = 0.97000d-03
        htt(19, 8, 2) = -0.45000d-03
        htt(19, 9, 2) = -0.18100d-02
        htt(19,10, 2) = -0.31500d-02
        htt(19,11, 2) = 0.26000d-03
        htt(19,12, 2) = -0.10800d-02
        htt(19,13, 2) = 0.67000d-03
        htt(19,14, 2) = -0.48000d-02
        htt(19,15, 2) = 0.20000d-02
        htt(19,16, 2) = -0.34700d-02
        htt(19,17, 2) = 0.19000d-02
        htt(19,18, 2) = -0.26700d-02
        htt(19,19, 2) = 0.73700d-02
        htt(20, 1, 2) = -0.77000d-03
        htt(20, 2, 2) = -0.76000d-03
        htt(20, 3, 2) = -0.45000d-03
        htt(20, 4, 2) = -0.44000d-03
        htt(20, 5, 2) = 0.80000d-04
        htt(20, 6, 2) = -0.74000d-03
        htt(20, 7, 2) = -0.23000d-03
        htt(20, 8, 2) = -0.10600d-02
        htt(20, 9, 2) = -0.59000d-03
        htt(20,10, 2) = -0.10400d-02
        htt(20,11, 2) = -0.60000d-03
        htt(20,12, 2) = -0.10500d-02
        htt(20,13, 2) = 0.26900d-02
        htt(20,14, 2) = -0.34800d-02
        htt(20,15, 2) = 0.31300d-02
        htt(20,16, 2) = -0.30400d-02
        htt(20,17, 2) = -0.33700d-02
        htt(20,18, 2) = 0.27100d-02
        htt(20,19, 2) = 0.28000d-02
        htt(20,20, 2) = 0.88800d-02
        htt(21, 1, 2) = 0.11400d-02
        htt(21, 2, 2) = 0.27200d-02
        htt(21, 3, 2) = -0.27900d-02
        htt(21, 4, 2) = -0.12100d-02
        htt(21, 5, 2) = -0.59100d-02
        htt(21, 6, 2) = -0.84000d-03
        htt(21, 7, 2) = -0.20900d-02
        htt(21, 8, 2) = 0.29800d-02
        htt(21, 9, 2) = 0.14100d-02
        htt(21,10, 2) = 0.53000d-03
        htt(21,11, 2) = -0.26000d-03
        htt(21,12, 2) = -0.11400d-02
        htt(21,13, 2) = 0.60000d-03
        htt(21,14, 2) = -0.15900d-02
        htt(21,15, 2) = 0.14800d-02
        htt(21,16, 2) = -0.71000d-03
        htt(21,17, 2) = -0.54000d-02
        htt(21,18, 2) = -0.23000d-03
        htt(21,19, 2) = -0.32100d-02
        htt(21,20, 2) = 0.19600d-02
        htt(21,21, 2) = 0.79800d-02
        htt(22, 1, 2) = -0.11900d-02
        htt(22, 2, 2) = 0.94000d-03
        htt(22, 3, 2) = -0.29200d-02
        htt(22, 4, 2) = -0.79000d-03
        htt(22, 5, 2) = -0.30000d-03
        htt(22, 6, 2) = -0.52400d-02
        htt(22, 7, 2) = 0.13800d-02
        htt(22, 8, 2) = -0.35600d-02
        htt(22, 9, 2) = 0.33000d-03
        htt(22,10, 2) = 0.12100d-02
        htt(22,11, 2) = -0.19200d-02
        htt(22,12, 2) = -0.10400d-02
        htt(22,13, 2) = 0.20500d-02
        htt(22,14, 2) = 0.20000d-03
        htt(22,15, 2) = 0.11700d-02
        htt(22,16, 2) = -0.67000d-03
        htt(22,17, 2) = -0.35500d-02
        htt(22,18, 2) = 0.98000d-03
        htt(22,19, 2) = -0.17100d-02
        htt(22,20, 2) = 0.28200d-02
        htt(22,21, 2) = 0.26300d-02
        htt(22,22, 2) = 0.78500d-02
        htt(23, 1, 2) = 0.24200d-02
        htt(23, 2, 2) = 0.20100d-02
        htt(23, 3, 2) = -0.73000d-03
        htt(23, 4, 2) = -0.11400d-02
        htt(23, 5, 2) = -0.39300d-02
        htt(23, 6, 2) = 0.53000d-03
        htt(23, 7, 2) = -0.87000d-03
        htt(23, 8, 2) = 0.35900d-02
        htt(23, 9, 2) = 0.18000d-03
        htt(23,10, 2) = -0.16000d-02
        htt(23,11, 2) = 0.61000d-03
        htt(23,12, 2) = -0.11700d-02
        htt(23,13, 2) = -0.14400d-02
        htt(23,14, 2) = -0.29100d-02
        htt(23,15, 2) = 0.33000d-03
        htt(23,16, 2) = -0.11400d-02
        htt(23,17, 2) = -0.80000d-04
        htt(23,18, 2) = -0.56500d-02
        htt(23,19, 2) = 0.14000d-02
        htt(23,20, 2) = -0.41800d-02
        htt(23,21, 2) = 0.27800d-02
        htt(23,22, 2) = -0.19400d-02
        htt(23,23, 2) = 0.84000d-02
        htt(24, 1, 2) = 0.90000d-04
        htt(24, 2, 2) = 0.23000d-03
        htt(24, 3, 2) = -0.86000d-03
        htt(24, 4, 2) = -0.71000d-03
        htt(24, 5, 2) = 0.16700d-02
        htt(24, 6, 2) = -0.38700d-02
        htt(24, 7, 2) = 0.25900d-02
        htt(24, 8, 2) = -0.29500d-02
        htt(24, 9, 2) = -0.90000d-03
        htt(24,10, 2) = -0.92000d-03
        htt(24,11, 2) = -0.10500d-02
        htt(24,12, 2) = -0.10700d-02
        htt(24,13, 2) = 0.10000d-04
        htt(24,14, 2) = -0.11200d-02
        htt(24,15, 2) = 0.30000d-04
        htt(24,16, 2) = -0.11000d-02
        htt(24,17, 2) = 0.17600d-02
        htt(24,18, 2) = -0.44400d-02
        htt(24,19, 2) = 0.29000d-02
        htt(24,20, 2) = -0.33100d-02
        htt(24,21, 2) = -0.25700d-02
        htt(24,22, 2) = 0.32800d-02
        htt(24,23, 2) = 0.36800d-02
        htt(24,24, 2) = 0.95300d-02
        ! Hessian elements for anchor point 2
        hbs( 1, 1, 2) = -0.25450d-01
        hbs( 1, 2, 2) = -0.24360d-01
        hbs( 1, 3, 2) = -0.45340d-01
        hbs( 1, 4, 2) = -0.86200d-02
        hbs( 1, 5, 2) = 0.45000d-02
        hbs( 1, 6, 2) = 0.42960d-01
        hbs( 1, 7, 2) = -0.49000d-03
        hbs( 1, 8, 2) = 0.43220d-01
        hbs( 1, 9, 2) = -0.40000d-04
        hbs( 1,10, 2) = -0.84400d-02
        hbs( 1,11, 2) = -0.54000d-03
        hbs( 1,12, 2) = 0.48600d-02
        hbs( 2, 1, 2) = 0.39900d-01
        hbs( 2, 2, 2) = -0.84900d-02
        hbs( 2, 3, 2) = 0.28370d-01
        hbs( 2, 4, 2) = 0.69300d-02
        hbs( 2, 5, 2) = -0.38100d-02
        hbs( 2, 6, 2) = -0.29340d-01
        hbs( 2, 7, 2) = 0.12900d-02
        hbs( 2, 8, 2) = -0.15540d-01
        hbs( 2, 9, 2) = 0.18000d-03
        hbs( 2,10, 2) = 0.28900d-02
        hbs( 2,11, 2) = -0.42000d-03
        hbs( 2,12, 2) = -0.40000d-03
        hbs( 3, 1, 2) = -0.14450d-01
        hbs( 3, 2, 2) = 0.32850d-01
        hbs( 3, 3, 2) = 0.16980d-01
        hbs( 3, 4, 2) = 0.16900d-02
        hbs( 3, 5, 2) = -0.70000d-03
        hbs( 3, 6, 2) = -0.13620d-01
        hbs( 3, 7, 2) = -0.79000d-03
        hbs( 3, 8, 2) = -0.27680d-01
        hbs( 3, 9, 2) = -0.14000d-03
        hbs( 3,10, 2) = 0.55400d-02
        hbs( 3,11, 2) = 0.95000d-03
        hbs( 3,12, 2) = -0.44600d-02
        hbs( 4, 1, 2) = -0.22690d-01
        hbs( 4, 2, 2) = -0.61800d-02
        hbs( 4, 3, 2) = 0.14530d-01
        hbs( 4, 4, 2) = -0.35450d-01
        hbs( 4, 5, 2) = -0.71600d-02
        hbs( 4, 6, 2) = -0.92200d-02
        hbs( 4, 7, 2) = 0.43300d-02
        hbs( 4, 8, 2) = 0.36360d-01
        hbs( 4, 9, 2) = -0.65000d-03
        hbs( 4,10, 2) = 0.49680d-01
        hbs( 4,11, 2) = 0.19000d-03
        hbs( 4,12, 2) = -0.97000d-03
        hbs( 5, 1, 2) = 0.21990d-01
        hbs( 5, 2, 2) = 0.37900d-02
        hbs( 5, 3, 2) = -0.71900d-02
        hbs( 5, 4, 2) = 0.47500d-02
        hbs( 5, 5, 2) = 0.95400d-02
        hbs( 5, 6, 2) = 0.43800d-02
        hbs( 5, 7, 2) = -0.13300d-02
        hbs( 5, 8, 2) = -0.15760d-01
        hbs( 5, 9, 2) = -0.15000d-03
        hbs( 5,10, 2) = -0.26690d-01
        hbs( 5,11, 2) = -0.13000d-03
        hbs( 5,12, 2) = 0.10900d-02
        hbs( 6, 1, 2) = 0.70000d-03
        hbs( 6, 2, 2) = 0.24000d-02
        hbs( 6, 3, 2) = -0.73400d-02
        hbs( 6, 4, 2) = 0.30700d-01
        hbs( 6, 5, 2) = -0.23800d-02
        hbs( 6, 6, 2) = 0.48500d-02
        hbs( 6, 7, 2) = -0.29900d-02
        hbs( 6, 8, 2) = -0.20600d-01
        hbs( 6, 9, 2) = 0.80000d-03
        hbs( 6,10, 2) = -0.22990d-01
        hbs( 6,11, 2) = -0.60000d-04
        hbs( 6,12, 2) = -0.12000d-03
        hbs( 7, 1, 2) = -0.73400d-02
        hbs( 7, 2, 2) = 0.28830d-01
        hbs( 7, 3, 2) = 0.66600d-02
        hbs( 7, 4, 2) = -0.33160d-01
        hbs( 7, 5, 2) = 0.40800d-02
        hbs( 7, 6, 2) = -0.36980d-01
        hbs( 7, 7, 2) = -0.81400d-02
        hbs( 7, 8, 2) = -0.48000d-02
        hbs( 7, 9, 2) = 0.47000d-02
        hbs( 7,10, 2) = 0.43210d-01
        hbs( 7,11, 2) = -0.64000d-03
        hbs( 7,12, 2) = -0.30000d-03
        hbs( 8, 1, 2) = 0.34500d-02
        hbs( 8, 2, 2) = -0.15860d-01
        hbs( 8, 3, 2) = -0.19800d-02
        hbs( 8, 4, 2) = 0.30660d-01
        hbs( 8, 5, 2) = -0.26700d-02
        hbs( 8, 6, 2) = 0.43400d-02
        hbs( 8, 7, 2) = 0.39600d-02
        hbs( 8, 8, 2) = 0.18400d-02
        hbs( 8, 9, 2) = -0.18900d-02
        hbs( 8,10, 2) = -0.19260d-01
        hbs( 8,11, 2) = -0.26000d-03
        hbs( 8,12, 2) = 0.16000d-03
        hbs( 9, 1, 2) = 0.38900d-02
        hbs( 9, 2, 2) = -0.12980d-01
        hbs( 9, 3, 2) = -0.46900d-02
        hbs( 9, 4, 2) = 0.25000d-02
        hbs( 9, 5, 2) = -0.14100d-02
        hbs( 9, 6, 2) = 0.32640d-01
        hbs( 9, 7, 2) = 0.41800d-02
        hbs( 9, 8, 2) = 0.29500d-02
        hbs( 9, 9, 2) = -0.28100d-02
        hbs( 9,10, 2) = -0.23950d-01
        hbs( 9,11, 2) = 0.90000d-03
        hbs( 9,12, 2) = 0.14000d-03
        hbs(10, 1, 2) = 0.32210d-01
        hbs(10, 2, 2) = 0.32440d-01
        hbs(10, 3, 2) = 0.24000d-02
        hbs(10, 4, 2) = -0.14180d-01
        hbs(10, 5, 2) = -0.19000d-03
        hbs(10, 6, 2) = -0.29000d-01
        hbs(10, 7, 2) = 0.46200d-02
        hbs(10, 8, 2) = -0.28530d-01
        hbs(10, 9, 2) = -0.80800d-02
        hbs(10,10, 2) = -0.14370d-01
        hbs(10,11, 2) = 0.45400d-02
        hbs(10,12, 2) = -0.45000d-03
        hbs(11, 1, 2) = -0.18020d-01
        hbs(11, 2, 2) = -0.14150d-01
        hbs(11, 3, 2) = -0.11000d-02
        hbs(11, 4, 2) = 0.83100d-02
        hbs(11, 5, 2) = 0.45000d-03
        hbs(11, 6, 2) = 0.27070d-01
        hbs(11, 7, 2) = -0.28400d-02
        hbs(11, 8, 2) = 0.14300d-02
        hbs(11, 9, 2) = 0.40300d-02
        hbs(11,10, 2) = 0.60300d-02
        hbs(11,11, 2) = -0.17200d-02
        hbs(11,12, 2) = -0.30000d-03
        hbs(12, 1, 2) = -0.14190d-01
        hbs(12, 2, 2) = -0.18290d-01
        hbs(12, 3, 2) = -0.13000d-02
        hbs(12, 4, 2) = 0.58700d-02
        hbs(12, 5, 2) = -0.26000d-03
        hbs(12, 6, 2) = 0.19300d-02
        hbs(12, 7, 2) = -0.17800d-02
        hbs(12, 8, 2) = 0.27100d-01
        hbs(12, 9, 2) = 0.40500d-02
        hbs(12,10, 2) = 0.83400d-02
        hbs(12,11, 2) = -0.28100d-02
        hbs(12,12, 2) = 0.75000d-03
        hbs(13, 1, 2) = 0.28370d-01
        hbs(13, 2, 2) = -0.63900d-02
        hbs(13, 3, 2) = 0.52100d-02
        hbs(13, 4, 2) = 0.42290d-01
        hbs(13, 5, 2) = -0.27000d-03
        hbs(13, 6, 2) = -0.35800d-02
        hbs(13, 7, 2) = -0.64000d-03
        hbs(13, 8, 2) = -0.36310d-01
        hbs(13, 9, 2) = 0.46600d-02
        hbs(13,10, 2) = -0.33260d-01
        hbs(13,11, 2) = -0.81200d-02
        hbs(13,12, 2) = 0.45000d-02
        hbs(14, 1, 2) = -0.12750d-01
        hbs(14, 2, 2) = 0.38700d-02
        hbs(14, 3, 2) = -0.40000d-02
        hbs(14, 4, 2) = -0.23630d-01
        hbs(14, 5, 2) = 0.12000d-03
        hbs(14, 6, 2) = 0.23100d-02
        hbs(14, 7, 2) = 0.89000d-03
        hbs(14, 8, 2) = 0.31980d-01
        hbs(14, 9, 2) = -0.28200d-02
        hbs(14,10, 2) = 0.30000d-02
        hbs(14,11, 2) = 0.42200d-02
        hbs(14,12, 2) = -0.17600d-02
        hbs(15, 1, 2) = -0.15620d-01
        hbs(15, 2, 2) = 0.25200d-02
        hbs(15, 3, 2) = -0.12100d-02
        hbs(15, 4, 2) = -0.18660d-01
        hbs(15, 5, 2) = 0.14000d-03
        hbs(15, 6, 2) = 0.12700d-02
        hbs(15, 7, 2) = -0.25000d-03
        hbs(15, 8, 2) = 0.43300d-02
        hbs(15, 9, 2) = -0.18400d-02
        hbs(15,10, 2) = 0.30260d-01
        hbs(15,11, 2) = 0.39000d-02
        hbs(15,12, 2) = -0.27400d-02
        hbs(16, 1, 2) = -0.51100d-02
        hbs(16, 2, 2) = -0.24340d-01
        hbs(16, 3, 2) = 0.16540d-01
        hbs(16, 4, 2) = 0.49110d-01
        hbs(16, 5, 2) = -0.98000d-03
        hbs(16, 6, 2) = 0.35820d-01
        hbs(16, 7, 2) = 0.33000d-03
        hbs(16, 8, 2) = -0.99500d-02
        hbs(16, 9, 2) = -0.59000d-03
        hbs(16,10, 2) = -0.36820d-01
        hbs(16,11, 2) = 0.45700d-02
        hbs(16,12, 2) = -0.76400d-02
        hbs(17, 1, 2) = 0.19500d-02
        hbs(17, 2, 2) = 0.23750d-01
        hbs(17, 3, 2) = -0.99800d-02
        hbs(17, 4, 2) = -0.26130d-01
        hbs(17, 5, 2) = 0.10100d-02
        hbs(17, 6, 2) = -0.15770d-01
        hbs(17, 7, 2) = -0.17000d-03
        hbs(17, 8, 2) = 0.47200d-02
        hbs(17, 9, 2) = -0.26000d-03
        hbs(17,10, 2) = 0.55900d-02
        hbs(17,11, 2) = -0.16700d-02
        hbs(17,12, 2) = 0.59100d-02
        hbs(18, 1, 2) = 0.31600d-02
        hbs(18, 2, 2) = 0.59000d-03
        hbs(18, 3, 2) = -0.65700d-02
        hbs(18, 4, 2) = -0.22990d-01
        hbs(18, 5, 2) = -0.30000d-04
        hbs(18, 6, 2) = -0.20050d-01
        hbs(18, 7, 2) = -0.16000d-03
        hbs(18, 8, 2) = 0.52300d-02
        hbs(18, 9, 2) = 0.85000d-03
        hbs(18,10, 2) = 0.31230d-01
        hbs(18,11, 2) = -0.29000d-02
        hbs(18,12, 2) = 0.17300d-02
        ! Hessian elements for anchor point 3
        hss( 1, 1, 3) = 0.30450d+00
        hss( 2, 1, 3) = -0.38000d-02
        hss( 2, 2, 3) = 0.30517d+00
        hss( 3, 1, 3) = 0.88850d-01
        hss( 3, 2, 3) = 0.88210d-01
        hss( 3, 3, 3) = 0.63344d+00
        hss( 4, 1, 3) = 0.83900d-02
        hss( 4, 2, 3) = 0.19510d-01
        hss( 4, 3, 3) = -0.15960d-01
        hss( 4, 4, 3) = 0.43559d+00
        hss( 5, 1, 3) = 0.15900d-02
        hss( 5, 2, 3) = 0.10000d-04
        hss( 5, 3, 3) = -0.79000d-03
        hss( 5, 4, 3) = 0.28500d-02
        hss( 5, 5, 3) = 0.36540d+00
        hss( 6, 1, 3) = 0.33900d-02
        hss( 6, 2, 3) = 0.49580d-01
        hss( 6, 3, 3) = 0.17260d-01
        hss( 6, 4, 3) = 0.43530d-01
        hss( 6, 5, 3) = -0.70000d-04
        hss( 6, 6, 3) = 0.35132d+00
        hss( 7, 1, 3) = -0.12200d-02
        hss( 7, 2, 3) = -0.41000d-03
        hss( 7, 3, 3) = 0.12400d-02
        hss( 7, 4, 3) = 0.52300d-02
        hss( 7, 5, 3) = 0.38000d-03
        hss( 7, 6, 3) = 0.38200d-02
        hss( 7, 7, 3) = 0.36408d+00
        hss( 8, 1, 3) = 0.49330d-01
        hss( 8, 2, 3) = 0.31600d-02
        hss( 8, 3, 3) = 0.17100d-01
        hss( 8, 4, 3) = -0.51700d-02
        hss( 8, 5, 3) = 0.13000d-03
        hss( 8, 6, 3) = 0.29020d-01
        hss( 8, 7, 3) = -0.29000d-03
        hss( 8, 8, 3) = 0.35081d+00
        hss( 9, 1, 3) = -0.28000d-03
        hss( 9, 2, 3) = -0.28000d-03
        hss( 9, 3, 3) = -0.38000d-03
        hss( 9, 4, 3) = -0.48000d-03
        hss( 9, 5, 3) = 0.60000d-04
        hss( 9, 6, 3) = 0.41500d-02
        hss( 9, 7, 3) = 0.40000d-03
        hss( 9, 8, 3) = 0.41400d-02
        hss( 9, 9, 3) = 0.36450d+00
        hss(10, 1, 3) = 0.19190d-01
        hss(10, 2, 3) = 0.79700d-02
        hss(10, 3, 3) = -0.16090d-01
        hss(10, 4, 3) = 0.80800d-01
        hss(10, 5, 3) = -0.65000d-03
        hss(10, 6, 3) = -0.49200d-02
        hss(10, 7, 3) = -0.53000d-03
        hss(10, 8, 3) = 0.43630d-01
        hss(10, 9, 3) = -0.50000d-03
        hss(10,10, 3) = 0.43670d+00
        hss(11, 1, 3) = -0.38000d-03
        hss(11, 2, 3) = -0.12000d-02
        hss(11, 3, 3) = 0.12500d-02
        hss(11, 4, 3) = -0.56000d-03
        hss(11, 5, 3) = 0.60000d-04
        hss(11, 6, 3) = -0.36000d-03
        hss(11, 7, 3) = 0.19000d-03
        hss(11, 8, 3) = 0.38300d-02
        hss(11, 9, 3) = 0.40000d-03
        hss(11,10, 3) = 0.51900d-02
        hss(11,11, 3) = 0.36404d+00
        hss(12, 1, 3) = 0.24000d-03
        hss(12, 2, 3) = 0.21300d-02
        hss(12, 3, 3) = 0.70000d-04
        hss(12, 4, 3) = -0.68000d-03
        hss(12, 5, 3) = 0.14000d-03
        hss(12, 6, 3) = -0.18000d-03
        hss(12, 7, 3) = 0.70000d-04
        hss(12, 8, 3) = -0.40000d-04
        hss(12, 9, 3) = 0.70000d-04
        hss(12,10, 3) = 0.28100d-02
        hss(12,11, 3) = 0.41000d-03
        hss(12,12, 3) = 0.36629d+00
        ! Hessian elements for anchor point 3
        hbb( 1, 1, 3) = 0.81840d-01
        hbb( 2, 1, 3) = -0.43040d-01
        hbb( 2, 2, 3) = 0.14305d+00
        hbb( 3, 1, 3) = -0.38800d-01
        hbb( 3, 2, 3) = -0.10001d+00
        hbb( 3, 3, 3) = 0.13881d+00
        hbb( 4, 1, 3) = -0.37990d-01
        hbb( 4, 2, 3) = 0.26100d-01
        hbb( 4, 3, 3) = 0.11900d-01
        hbb( 4, 4, 3) = 0.77580d-01
        hbb( 5, 1, 3) = 0.23090d-01
        hbb( 5, 2, 3) = -0.18080d-01
        hbb( 5, 3, 3) = -0.50100d-02
        hbb( 5, 4, 3) = -0.38630d-01
        hbb( 5, 5, 3) = 0.74170d-01
        hbb( 6, 1, 3) = 0.14900d-01
        hbb( 6, 2, 3) = -0.80200d-02
        hbb( 6, 3, 3) = -0.68800d-02
        hbb( 6, 4, 3) = -0.38950d-01
        hbb( 6, 5, 3) = -0.35540d-01
        hbb( 6, 6, 3) = 0.74490d-01
        hbb( 7, 1, 3) = -0.82000d-02
        hbb( 7, 2, 3) = -0.62200d-02
        hbb( 7, 3, 3) = 0.14420d-01
        hbb( 7, 4, 3) = -0.29360d-01
        hbb( 7, 5, 3) = 0.11220d-01
        hbb( 7, 6, 3) = 0.18140d-01
        hbb( 7, 7, 3) = 0.74660d-01
        hbb( 8, 1, 3) = -0.15800d-02
        hbb( 8, 2, 3) = 0.80800d-02
        hbb( 8, 3, 3) = -0.65000d-02
        hbb( 8, 4, 3) = 0.17600d-01
        hbb( 8, 5, 3) = -0.61400d-02
        hbb( 8, 6, 3) = -0.11460d-01
        hbb( 8, 7, 3) = -0.36880d-01
        hbb( 8, 8, 3) = 0.78170d-01
        hbb( 9, 1, 3) = 0.97800d-02
        hbb( 9, 2, 3) = -0.18600d-02
        hbb( 9, 3, 3) = -0.79200d-02
        hbb( 9, 4, 3) = 0.11760d-01
        hbb( 9, 5, 3) = -0.50800d-02
        hbb( 9, 6, 3) = -0.66800d-02
        hbb( 9, 7, 3) = -0.37780d-01
        hbb( 9, 8, 3) = -0.41290d-01
        hbb( 9, 9, 3) = 0.79060d-01
        hbb(10, 1, 3) = 0.10690d-01
        hbb(10, 2, 3) = -0.58300d-02
        hbb(10, 3, 3) = -0.48700d-02
        hbb(10, 4, 3) = -0.91900d-02
        hbb(10, 5, 3) = 0.10240d-01
        hbb(10, 6, 3) = -0.10500d-02
        hbb(10, 7, 3) = -0.34370d-01
        hbb(10, 8, 3) = 0.13850d-01
        hbb(10, 9, 3) = 0.20520d-01
        hbb(10,10, 3) = 0.76450d-01
        hbb(11, 1, 3) = -0.53700d-02
        hbb(11, 2, 3) = 0.34700d-02
        hbb(11, 3, 3) = 0.19000d-02
        hbb(11, 4, 3) = -0.11400d-02
        hbb(11, 5, 3) = -0.35400d-02
        hbb(11, 6, 3) = 0.46900d-02
        hbb(11, 7, 3) = 0.20260d-01
        hbb(11, 8, 3) = -0.74900d-02
        hbb(11, 9, 3) = -0.12770d-01
        hbb(11,10, 3) = -0.38210d-01
        hbb(11,11, 3) = 0.79260d-01
        hbb(12, 1, 3) = -0.53200d-02
        hbb(12, 2, 3) = 0.23500d-02
        hbb(12, 3, 3) = 0.29700d-02
        hbb(12, 4, 3) = 0.10330d-01
        hbb(12, 5, 3) = -0.67000d-02
        hbb(12, 6, 3) = -0.36300d-02
        hbb(12, 7, 3) = 0.14110d-01
        hbb(12, 8, 3) = -0.63600d-02
        hbb(12, 9, 3) = -0.77500d-02
        hbb(12,10, 3) = -0.38240d-01
        hbb(12,11, 3) = -0.41050d-01
        hbb(12,12, 3) = 0.79290d-01
        hbb(13, 1, 3) = -0.77600d-02
        hbb(13, 2, 3) = 0.14880d-01
        hbb(13, 3, 3) = -0.71200d-02
        hbb(13, 4, 3) = 0.36800d-02
        hbb(13, 5, 3) = -0.17700d-02
        hbb(13, 6, 3) = -0.19200d-02
        hbb(13, 7, 3) = -0.65800d-02
        hbb(13, 8, 3) = 0.88500d-02
        hbb(13, 9, 3) = -0.22700d-02
        hbb(13,10, 3) = -0.34350d-01
        hbb(13,11, 3) = 0.14110d-01
        hbb(13,12, 3) = 0.20240d-01
        hbb(13,13, 3) = 0.74690d-01
        hbb(14, 1, 3) = 0.95100d-02
        hbb(14, 2, 3) = -0.81200d-02
        hbb(14, 3, 3) = -0.13900d-02
        hbb(14, 4, 3) = -0.19200d-02
        hbb(14, 5, 3) = 0.75000d-03
        hbb(14, 6, 3) = 0.11700d-02
        hbb(14, 7, 3) = -0.22400d-02
        hbb(14, 8, 3) = -0.29900d-02
        hbb(14, 9, 3) = 0.52300d-02
        hbb(14,10, 3) = 0.20540d-01
        hbb(14,11, 3) = -0.77700d-02
        hbb(14,12, 3) = -0.12770d-01
        hbb(14,13, 3) = -0.37840d-01
        hbb(14,14, 3) = 0.79020d-01
        hbb(15, 1, 3) = -0.17500d-02
        hbb(15, 2, 3) = -0.67600d-02
        hbb(15, 3, 3) = 0.85100d-02
        hbb(15, 4, 3) = -0.17600d-02
        hbb(15, 5, 3) = 0.10100d-02
        hbb(15, 6, 3) = 0.75000d-03
        hbb(15, 7, 3) = 0.88200d-02
        hbb(15, 8, 3) = -0.58500d-02
        hbb(15, 9, 3) = -0.29600d-02
        hbb(15,10, 3) = 0.13810d-01
        hbb(15,11, 3) = -0.63400d-02
        hbb(15,12, 3) = -0.74700d-02
        hbb(15,13, 3) = -0.36860d-01
        hbb(15,14, 3) = -0.41180d-01
        hbb(15,15, 3) = 0.78040d-01
        hbb(16, 1, 3) = -0.38580d-01
        hbb(16, 2, 3) = 0.14110d-01
        hbb(16, 3, 3) = 0.24470d-01
        hbb(16, 4, 3) = -0.47200d-02
        hbb(16, 5, 3) = -0.41600d-02
        hbb(16, 6, 3) = 0.88700d-02
        hbb(16, 7, 3) = 0.38500d-02
        hbb(16, 8, 3) = -0.18400d-02
        hbb(16, 9, 3) = -0.20100d-02
        hbb(16,10, 3) = -0.92300d-02
        hbb(16,11, 3) = 0.10360d-01
        hbb(16,12, 3) = -0.11200d-02
        hbb(16,13, 3) = -0.29700d-01
        hbb(16,14, 3) = 0.11950d-01
        hbb(16,15, 3) = 0.17750d-01
        hbb(16,16, 3) = 0.78380d-01
        hbb(17, 1, 3) = 0.23970d-01
        hbb(17, 2, 3) = -0.47800d-02
        hbb(17, 3, 3) = -0.19190d-01
        hbb(17, 4, 3) = -0.44800d-02
        hbb(17, 5, 3) = 0.69900d-02
        hbb(17, 6, 3) = -0.25200d-02
        hbb(17, 7, 3) = -0.18900d-02
        hbb(17, 8, 3) = 0.10700d-02
        hbb(17, 9, 3) = 0.82000d-03
        hbb(17,10, 3) = 0.10360d-01
        hbb(17,11, 3) = -0.67400d-02
        hbb(17,12, 3) = -0.36100d-02
        hbb(17,13, 3) = 0.11220d-01
        hbb(17,14, 3) = -0.51000d-02
        hbb(17,15, 3) = -0.61200d-02
        hbb(17,16, 3) = -0.39180d-01
        hbb(17,17, 3) = 0.74630d-01
        hbb(18, 1, 3) = 0.14610d-01
        hbb(18, 2, 3) = -0.93300d-02
        hbb(18, 3, 3) = -0.52800d-02
        hbb(18, 4, 3) = 0.91900d-02
        hbb(18, 5, 3) = -0.28400d-02
        hbb(18, 6, 3) = -0.63600d-02
        hbb(18, 7, 3) = -0.19600d-02
        hbb(18, 8, 3) = 0.77000d-03
        hbb(18, 9, 3) = 0.11900d-02
        hbb(18,10, 3) = -0.11200d-02
        hbb(18,11, 3) = -0.36100d-02
        hbb(18,12, 3) = 0.47400d-02
        hbb(18,13, 3) = 0.18470d-01
        hbb(18,14, 3) = -0.68500d-02
        hbb(18,15, 3) = -0.11630d-01
        hbb(18,16, 3) = -0.39190d-01
        hbb(18,17, 3) = -0.35450d-01
        hbb(18,18, 3) = 0.74640d-01
        ! Hessian elements for anchor point 3
        htt( 1, 1, 3) = 0.86000d-02
        htt( 2, 1, 3) = 0.34800d-02
        htt( 2, 2, 3) = 0.76600d-02
        htt( 3, 1, 3) = -0.57000d-03
        htt( 3, 2, 3) = -0.37600d-02
        htt( 3, 3, 3) = 0.75500d-02
        htt( 4, 1, 3) = -0.56900d-02
        htt( 4, 2, 3) = 0.42000d-03
        htt( 4, 3, 3) = 0.43600d-02
        htt( 4, 4, 3) = 0.10470d-01
        htt( 5, 1, 3) = -0.61300d-02
        htt( 5, 2, 3) = -0.40700d-02
        htt( 5, 3, 3) = 0.29900d-02
        htt( 5, 4, 3) = 0.50400d-02
        htt( 5, 5, 3) = 0.84500d-02
        htt( 6, 1, 3) = -0.42600d-02
        htt( 6, 2, 3) = -0.28700d-02
        htt( 6, 3, 3) = 0.30200d-02
        htt( 6, 4, 3) = 0.44100d-02
        htt( 6, 5, 3) = 0.35900d-02
        htt( 6, 6, 3) = 0.78600d-02
        htt( 7, 1, 3) = 0.30000d-02
        htt( 7, 2, 3) = 0.31400d-02
        htt( 7, 3, 3) = -0.51000d-02
        htt( 7, 4, 3) = -0.49600d-02
        htt( 7, 5, 3) = -0.62000d-03
        htt( 7, 6, 3) = -0.36600d-02
        htt( 7, 7, 3) = 0.74400d-02
        htt( 8, 1, 3) = 0.48700d-02
        htt( 8, 2, 3) = 0.43400d-02
        htt( 8, 3, 3) = -0.50700d-02
        htt( 8, 4, 3) = -0.56000d-02
        htt( 8, 5, 3) = -0.54900d-02
        htt( 8, 6, 3) = 0.61000d-03
        htt( 8, 7, 3) = 0.44100d-02
        htt( 8, 8, 3) = 0.10510d-01
        htt( 9, 1, 3) = -0.58900d-02
        htt( 9, 2, 3) = -0.93000d-03
        htt( 9, 3, 3) = -0.18100d-02
        htt( 9, 4, 3) = 0.31600d-02
        htt( 9, 5, 3) = 0.11100d-02
        htt( 9, 6, 3) = 0.22400d-02
        htt( 9, 7, 3) = -0.29500d-02
        htt( 9, 8, 3) = -0.18200d-02
        htt( 9, 9, 3) = 0.80600d-02
        htt(10, 1, 3) = -0.39400d-02
        htt(10, 2, 3) = 0.44000d-03
        htt(10, 3, 3) = -0.58000d-03
        htt(10, 4, 3) = 0.37900d-02
        htt(10, 5, 3) = 0.23700d-02
        htt(10, 6, 3) = 0.17100d-02
        htt(10, 7, 3) = -0.97000d-03
        htt(10, 8, 3) = -0.16300d-02
        htt(10, 9, 3) = 0.28200d-02
        htt(10,10, 3) = 0.84500d-02
        htt(11, 1, 3) = -0.51000d-03
        htt(11, 2, 3) = -0.53200d-02
        htt(11, 3, 3) = 0.15400d-02
        htt(11, 4, 3) = -0.32600d-02
        htt(11, 5, 3) = -0.10500d-02
        htt(11, 6, 3) = 0.78000d-03
        htt(11, 7, 3) = -0.31000d-02
        htt(11, 8, 3) = -0.12700d-02
        htt(11, 9, 3) = 0.28400d-02
        htt(11,10, 3) = -0.17800d-02
        htt(11,11, 3) = 0.78900d-02
        htt(12, 1, 3) = 0.14400d-02
        htt(12, 2, 3) = -0.39600d-02
        htt(12, 3, 3) = 0.27700d-02
        htt(12, 4, 3) = -0.26300d-02
        htt(12, 5, 3) = 0.21000d-03
        htt(12, 6, 3) = 0.25000d-03
        htt(12, 7, 3) = -0.11200d-02
        htt(12, 8, 3) = -0.10700d-02
        htt(12, 9, 3) = -0.24000d-02
        htt(12,10, 3) = 0.38400d-02
        htt(12,11, 3) = 0.32800d-02
        htt(12,12, 3) = 0.95300d-02
        htt(13, 1, 3) = 0.46000d-03
        htt(13, 2, 3) = -0.12300d-02
        htt(13, 3, 3) = 0.19200d-02
        htt(13, 4, 3) = 0.23000d-03
        htt(13, 5, 3) = 0.19500d-02
        htt(13, 6, 3) = 0.60000d-03
        htt(13, 7, 3) = 0.50000d-03
        htt(13, 8, 3) = -0.85000d-03
        htt(13, 9, 3) = -0.54600d-02
        htt(13,10, 3) = -0.30000d-04
        htt(13,11, 3) = -0.36800d-02
        htt(13,12, 3) = 0.17400d-02
        htt(13,13, 3) = 0.82300d-02
        htt(14, 1, 3) = 0.18200d-02
        htt(14, 2, 3) = 0.48000d-03
        htt(14, 3, 3) = 0.11000d-02
        htt(14, 4, 3) = -0.24000d-03
        htt(14, 5, 3) = 0.60000d-03
        htt(14, 6, 3) = -0.11700d-02
        htt(14, 7, 3) = 0.13200d-02
        htt(14, 8, 3) = -0.45000d-03
        htt(14, 9, 3) = -0.31800d-02
        htt(14,10, 3) = 0.14500d-02
        htt(14,11, 3) = -0.17700d-02
        htt(14,12, 3) = 0.28600d-02
        htt(14,13, 3) = 0.21600d-02
        htt(14,14, 3) = 0.73500d-02
        htt(15, 1, 3) = -0.14800d-02
        htt(15, 2, 3) = -0.25800d-02
        htt(15, 3, 3) = 0.70000d-03
        htt(15, 4, 3) = -0.40000d-03
        htt(15, 5, 3) = 0.70000d-03
        htt(15, 6, 3) = 0.11300d-02
        htt(15, 7, 3) = -0.14600d-02
        htt(15, 8, 3) = -0.10400d-02
        htt(15, 9, 3) = -0.26000d-03
        htt(15,10, 3) = -0.56100d-02
        htt(15,11, 3) = 0.90000d-03
        htt(15,12, 3) = -0.44500d-02
        htt(15,13, 3) = 0.28500d-02
        htt(15,14, 3) = -0.24300d-02
        htt(15,15, 3) = 0.81500d-02
        htt(16, 1, 3) = -0.11000d-03
        htt(16, 2, 3) = -0.87000d-03
        htt(16, 3, 3) = -0.11000d-03
        htt(16, 4, 3) = -0.87000d-03
        htt(16, 5, 3) = -0.64000d-03
        htt(16, 6, 3) = -0.64000d-03
        htt(16, 7, 3) = -0.64000d-03
        htt(16, 8, 3) = -0.64000d-03
        htt(16, 9, 3) = 0.20100d-02
        htt(16,10, 3) = -0.41300d-02
        htt(16,11, 3) = 0.28100d-02
        htt(16,12, 3) = -0.33300d-02
        htt(16,13, 3) = -0.32200d-02
        htt(16,14, 3) = 0.27600d-02
        htt(16,15, 3) = 0.28700d-02
        htt(16,16, 3) = 0.88500d-02
        htt(17, 1, 3) = 0.20300d-02
        htt(17, 2, 3) = 0.64000d-03
        htt(17, 3, 3) = 0.51000d-03
        htt(17, 4, 3) = -0.88000d-03
        htt(17, 5, 3) = 0.38000d-03
        htt(17, 6, 3) = -0.12800d-02
        htt(17, 7, 3) = 0.19000d-02
        htt(17, 8, 3) = 0.24000d-03
        htt(17, 9, 3) = 0.64000d-03
        htt(17,10, 3) = -0.15500d-02
        htt(17,11, 3) = 0.21000d-02
        htt(17,12, 3) = -0.80000d-04
        htt(17,13, 3) = -0.58000d-02
        htt(17,14, 3) = 0.29000d-03
        htt(17,15, 3) = -0.36300d-02
        htt(17,16, 3) = 0.24600d-02
        htt(17,17, 3) = 0.82200d-02
        htt(18, 1, 3) = 0.77000d-03
        htt(18, 2, 3) = 0.11200d-02
        htt(18, 3, 3) = -0.14100d-02
        htt(18, 4, 3) = -0.10600d-02
        htt(18, 5, 3) = -0.15200d-02
        htt(18, 6, 3) = -0.26700d-02
        htt(18, 7, 3) = 0.66000d-03
        htt(18, 8, 3) = -0.50000d-03
        htt(18, 9, 3) = 0.17500d-02
        htt(18,10, 3) = 0.41000d-03
        htt(18,11, 3) = 0.13800d-02
        htt(18,12, 3) = 0.40000d-04
        htt(18,13, 3) = -0.36300d-02
        htt(18,14, 3) = 0.16800d-02
        htt(18,15, 3) = -0.23000d-02
        htt(18,16, 3) = 0.30100d-02
        htt(18,17, 3) = 0.28700d-02
        htt(18,18, 3) = 0.81700d-02
        htt(19, 1, 3) = 0.67000d-03
        htt(19, 2, 3) = -0.10700d-02
        htt(19, 3, 3) = 0.13300d-02
        htt(19, 4, 3) = -0.41000d-03
        htt(19, 5, 3) = 0.17300d-02
        htt(19, 6, 3) = 0.49000d-03
        htt(19, 7, 3) = 0.10800d-02
        htt(19, 8, 3) = -0.16000d-03
        htt(19, 9, 3) = -0.16400d-02
        htt(19,10, 3) = -0.30300d-02
        htt(19,11, 3) = 0.19000d-03
        htt(19,12, 3) = -0.12000d-02
        htt(19,13, 3) = 0.27000d-03
        htt(19,14, 3) = -0.49100d-02
        htt(19,15, 3) = 0.16600d-02
        htt(19,16, 3) = -0.35200d-02
        htt(19,17, 3) = 0.21400d-02
        htt(19,18, 3) = -0.24400d-02
        htt(19,19, 3) = 0.73200d-02
        htt(20, 1, 3) = -0.59000d-03
        htt(20, 2, 3) = -0.59000d-03
        htt(20, 3, 3) = -0.60000d-03
        htt(20, 4, 3) = -0.59000d-03
        htt(20, 5, 3) = -0.17000d-03
        htt(20, 6, 3) = -0.90000d-03
        htt(20, 7, 3) = -0.16000d-03
        htt(20, 8, 3) = -0.90000d-03
        htt(20, 9, 3) = -0.53000d-03
        htt(20,10, 3) = -0.10700d-02
        htt(20,11, 3) = -0.54000d-03
        htt(20,12, 3) = -0.10700d-02
        htt(20,13, 3) = 0.24500d-02
        htt(20,14, 3) = -0.35100d-02
        htt(20,15, 3) = 0.29800d-02
        htt(20,16, 3) = -0.29800d-02
        htt(20,17, 3) = -0.32100d-02
        htt(20,18, 3) = 0.28600d-02
        htt(20,19, 3) = 0.27400d-02
        htt(20,20, 3) = 0.88200d-02
        htt(21, 1, 3) = 0.99000d-03
        htt(21, 2, 3) = 0.20900d-02
        htt(21, 3, 3) = -0.29700d-02
        htt(21, 4, 3) = -0.18700d-02
        htt(21, 5, 3) = -0.57000d-02
        htt(21, 6, 3) = -0.91000d-03
        htt(21, 7, 3) = -0.17600d-02
        htt(21, 8, 3) = 0.30300d-02
        htt(21, 9, 3) = 0.14100d-02
        htt(21,10, 3) = 0.29000d-03
        htt(21,11, 3) = 0.26000d-03
        htt(21,12, 3) = -0.86000d-03
        htt(21,13, 3) = 0.68000d-03
        htt(21,14, 3) = -0.16300d-02
        htt(21,15, 3) = 0.18000d-02
        htt(21,16, 3) = -0.51000d-03
        htt(21,17, 3) = -0.54200d-02
        htt(21,18, 3) = -0.27000d-03
        htt(21,19, 3) = -0.31100d-02
        htt(21,20, 3) = 0.20400d-02
        htt(21,21, 3) = 0.79000d-02
        htt(22, 1, 3) = -0.98000d-03
        htt(22, 2, 3) = 0.81000d-03
        htt(22, 3, 3) = -0.30000d-02
        htt(22, 4, 3) = -0.12000d-02
        htt(22, 5, 3) = -0.56000d-03
        htt(22, 6, 3) = -0.54200d-02
        htt(22, 7, 3) = 0.14400d-02
        htt(22, 8, 3) = -0.34100d-02
        htt(22, 9, 3) = 0.22000d-03
        htt(22,10, 3) = 0.98000d-03
        htt(22,11, 3) = -0.16700d-02
        htt(22,12, 3) = -0.91000d-03
        htt(22,13, 3) = 0.21100d-02
        htt(22,14, 3) = 0.24000d-03
        htt(22,15, 3) = 0.13500d-02
        htt(22,16, 3) = -0.51000d-03
        htt(22,17, 3) = -0.36600d-02
        htt(22,18, 3) = 0.95000d-03
        htt(22,19, 3) = -0.18000d-02
        htt(22,20, 3) = 0.28200d-02
        htt(22,21, 3) = 0.28500d-02
        htt(22,22, 3) = 0.79700d-02
        htt(23, 1, 3) = 0.22700d-02
        htt(23, 2, 3) = 0.16000d-02
        htt(23, 3, 3) = -0.10300d-02
        htt(23, 4, 3) = -0.16900d-02
        htt(23, 5, 3) = -0.37900d-02
        htt(23, 6, 3) = 0.50000d-03
        htt(23, 7, 3) = -0.51000d-03
        htt(23, 8, 3) = 0.37700d-02
        htt(23, 9, 3) = 0.30000d-03
        htt(23,10, 3) = -0.16900d-02
        htt(23,11, 3) = 0.10000d-02
        htt(23,12, 3) = -0.99000d-03
        htt(23,13, 3) = -0.15100d-02
        htt(23,14, 3) = -0.30400d-02
        htt(23,15, 3) = 0.46000d-03
        htt(23,16, 3) = -0.10700d-02
        htt(23,17, 3) = -0.20000d-04
        htt(23,18, 3) = -0.56100d-02
        htt(23,19, 3) = 0.15100d-02
        htt(23,20, 3) = -0.40900d-02
        htt(23,21, 3) = 0.27100d-02
        htt(23,22, 3) = -0.18100d-02
        htt(23,23, 3) = 0.83500d-02
        htt(24, 1, 3) = 0.29000d-03
        htt(24, 2, 3) = 0.33000d-03
        htt(24, 3, 3) = -0.10600d-02
        htt(24, 4, 3) = -0.10200d-02
        htt(24, 5, 3) = 0.13500d-02
        htt(24, 6, 3) = -0.40100d-02
        htt(24, 7, 3) = 0.27000d-02
        htt(24, 8, 3) = -0.26700d-02
        htt(24, 9, 3) = -0.90000d-03
        htt(24,10, 3) = -0.10000d-02
        htt(24,11, 3) = -0.94000d-03
        htt(24,12, 3) = -0.10400d-02
        htt(24,13, 3) = -0.90000d-04
        htt(24,14, 3) = -0.11700d-02
        htt(24,15, 3) = 0.10000d-04
        htt(24,16, 3) = -0.10700d-02
        htt(24,17, 3) = 0.17400d-02
        htt(24,18, 3) = -0.43900d-02
        htt(24,19, 3) = 0.28200d-02
        htt(24,20, 3) = -0.33100d-02
        htt(24,21, 3) = -0.23500d-02
        htt(24,22, 3) = 0.33100d-02
        htt(24,23, 3) = 0.38300d-02
        htt(24,24, 3) = 0.94900d-02
        ! Hessian elements for anchor point 3
        hbs( 1, 1, 3) = -0.23260d-01
        hbs( 1, 2, 3) = -0.22950d-01
        hbs( 1, 3, 3) = -0.47950d-01
        hbs( 1, 4, 3) = -0.88000d-02
        hbs( 1, 5, 3) = 0.57100d-02
        hbs( 1, 6, 3) = 0.44410d-01
        hbs( 1, 7, 3) = -0.55000d-03
        hbs( 1, 8, 3) = 0.44440d-01
        hbs( 1, 9, 3) = -0.20000d-04
        hbs( 1,10, 3) = -0.86200d-02
        hbs( 1,11, 3) = -0.54000d-03
        hbs( 1,12, 3) = 0.53300d-02
        hbs( 2, 1, 3) = 0.35140d-01
        hbs( 2, 2, 3) = -0.98400d-02
        hbs( 2, 3, 3) = 0.27820d-01
        hbs( 2, 4, 3) = 0.60300d-02
        hbs( 2, 5, 3) = -0.65700d-02
        hbs( 2, 6, 3) = -0.28710d-01
        hbs( 2, 7, 3) = 0.12300d-02
        hbs( 2, 8, 3) = -0.16910d-01
        hbs( 2, 9, 3) = 0.50000d-04
        hbs( 2,10, 3) = 0.28400d-02
        hbs( 2,11, 3) = -0.71000d-03
        hbs( 2,12, 3) = -0.49000d-03
        hbs( 3, 1, 3) = -0.11880d-01
        hbs( 3, 2, 3) = 0.32790d-01
        hbs( 3, 3, 3) = 0.20130d-01
        hbs( 3, 4, 3) = 0.27700d-02
        hbs( 3, 5, 3) = 0.86000d-03
        hbs( 3, 6, 3) = -0.15690d-01
        hbs( 3, 7, 3) = -0.68000d-03
        hbs( 3, 8, 3) = -0.27530d-01
        hbs( 3, 9, 3) = -0.30000d-04
        hbs( 3,10, 3) = 0.57800d-02
        hbs( 3,11, 3) = 0.12500d-02
        hbs( 3,12, 3) = -0.48400d-02
        hbs( 4, 1, 3) = -0.23570d-01
        hbs( 4, 2, 3) = -0.60500d-02
        hbs( 4, 3, 3) = 0.14830d-01
        hbs( 4, 4, 3) = -0.35570d-01
        hbs( 4, 5, 3) = -0.79300d-02
        hbs( 4, 6, 3) = -0.10760d-01
        hbs( 4, 7, 3) = 0.47300d-02
        hbs( 4, 8, 3) = 0.37200d-01
        hbs( 4, 9, 3) = -0.64000d-03
        hbs( 4,10, 3) = 0.48250d-01
        hbs( 4,11, 3) = 0.25000d-03
        hbs( 4,12, 3) = -0.12100d-02
        hbs( 5, 1, 3) = 0.21940d-01
        hbs( 5, 2, 3) = 0.33300d-02
        hbs( 5, 3, 3) = -0.76700d-02
        hbs( 5, 4, 3) = 0.47000d-02
        hbs( 5, 5, 3) = 0.57400d-02
        hbs( 5, 6, 3) = 0.50100d-02
        hbs( 5, 7, 3) = -0.16700d-02
        hbs( 5, 8, 3) = -0.16330d-01
        hbs( 5, 9, 3) = -0.26000d-03
        hbs( 5,10, 3) = -0.25900d-01
        hbs( 5,11, 3) = -0.12000d-03
        hbs( 5,12, 3) = 0.12000d-02
        hbs( 6, 1, 3) = 0.16300d-02
        hbs( 6, 2, 3) = 0.27200d-02
        hbs( 6, 3, 3) = -0.71600d-02
        hbs( 6, 4, 3) = 0.30860d-01
        hbs( 6, 5, 3) = 0.21900d-02
        hbs( 6, 6, 3) = 0.57500d-02
        hbs( 6, 7, 3) = -0.30600d-02
        hbs( 6, 8, 3) = -0.20870d-01
        hbs( 6, 9, 3) = 0.90000d-03
        hbs( 6,10, 3) = -0.22350d-01
        hbs( 6,11, 3) = -0.13000d-03
        hbs( 6,12, 3) = 0.10000d-04
        hbs( 7, 1, 3) = -0.76900d-02
        hbs( 7, 2, 3) = 0.28950d-01
        hbs( 7, 3, 3) = 0.62100d-02
        hbs( 7, 4, 3) = -0.33090d-01
        hbs( 7, 5, 3) = 0.43500d-02
        hbs( 7, 6, 3) = -0.35980d-01
        hbs( 7, 7, 3) = -0.84400d-02
        hbs( 7, 8, 3) = -0.54300d-02
        hbs( 7, 9, 3) = 0.48600d-02
        hbs( 7,10, 3) = 0.43110d-01
        hbs( 7,11, 3) = -0.75000d-03
        hbs( 7,12, 3) = -0.33000d-03
        hbs( 8, 1, 3) = 0.32000d-02
        hbs( 8, 2, 3) = -0.15810d-01
        hbs( 8, 3, 3) = -0.15400d-02
        hbs( 8, 4, 3) = 0.30210d-01
        hbs( 8, 5, 3) = -0.26700d-02
        hbs( 8, 6, 3) = 0.42000d-02
        hbs( 8, 7, 3) = 0.39400d-02
        hbs( 8, 8, 3) = 0.20700d-02
        hbs( 8, 9, 3) = -0.19800d-02
        hbs( 8,10, 3) = -0.19200d-01
        hbs( 8,11, 3) = -0.22000d-03
        hbs( 8,12, 3) = 0.19000d-03
        hbs( 9, 1, 3) = 0.44900d-02
        hbs( 9, 2, 3) = -0.13140d-01
        hbs( 9, 3, 3) = -0.46700d-02
        hbs( 9, 4, 3) = 0.28800d-02
        hbs( 9, 5, 3) = -0.16800d-02
        hbs( 9, 6, 3) = 0.31780d-01
        hbs( 9, 7, 3) = 0.45000d-02
        hbs( 9, 8, 3) = 0.33600d-02
        hbs( 9, 9, 3) = -0.28900d-02
        hbs( 9,10, 3) = -0.23910d-01
        hbs( 9,11, 3) = 0.97000d-03
        hbs( 9,12, 3) = 0.14000d-03
        hbs(10, 1, 3) = 0.31310d-01
        hbs(10, 2, 3) = 0.31460d-01
        hbs(10, 3, 3) = 0.56600d-02
        hbs(10, 4, 3) = -0.13750d-01
        hbs(10, 5, 3) = -0.28000d-03
        hbs(10, 6, 3) = -0.29490d-01
        hbs(10, 7, 3) = 0.47400d-02
        hbs(10, 8, 3) = -0.29500d-01
        hbs(10, 9, 3) = -0.84400d-02
        hbs(10,10, 3) = -0.13800d-01
        hbs(10,11, 3) = 0.47800d-02
        hbs(10,12, 3) = -0.37000d-03
        hbs(11, 1, 3) = -0.17710d-01
        hbs(11, 2, 3) = -0.13600d-01
        hbs(11, 3, 3) = -0.28200d-02
        hbs(11, 4, 3) = 0.80200d-02
        hbs(11, 5, 3) = 0.65000d-03
        hbs(11, 6, 3) = 0.27440d-01
        hbs(11, 7, 3) = -0.28800d-02
        hbs(11, 8, 3) = 0.19800d-02
        hbs(11, 9, 3) = 0.42300d-02
        hbs(11,10, 3) = 0.57200d-02
        hbs(11,11, 3) = -0.18700d-02
        hbs(11,12, 3) = -0.35000d-03
        hbs(12, 1, 3) = -0.13600d-01
        hbs(12, 2, 3) = -0.17860d-01
        hbs(12, 3, 3) = -0.28400d-02
        hbs(12, 4, 3) = 0.57400d-02
        hbs(12, 5, 3) = -0.37000d-03
        hbs(12, 6, 3) = 0.20500d-02
        hbs(12, 7, 3) = -0.18500d-02
        hbs(12, 8, 3) = 0.27520d-01
        hbs(12, 9, 3) = 0.42000d-02
        hbs(12,10, 3) = 0.80800d-02
        hbs(12,11, 3) = -0.29100d-02
        hbs(12,12, 3) = 0.72000d-03
        hbs(13, 1, 3) = 0.28850d-01
        hbs(13, 2, 3) = -0.76600d-02
        hbs(13, 3, 3) = 0.57200d-02
        hbs(13, 4, 3) = 0.43010d-01
        hbs(13, 5, 3) = -0.30000d-03
        hbs(13, 6, 3) = -0.53400d-02
        hbs(13, 7, 3) = -0.74000d-03
        hbs(13, 8, 3) = -0.35790d-01
        hbs(13, 9, 3) = 0.48800d-02
        hbs(13,10, 3) = -0.33150d-01
        hbs(13,11, 3) = -0.84400d-02
        hbs(13,12, 3) = 0.46600d-02
        hbs(14, 1, 3) = -0.13040d-01
        hbs(14, 2, 3) = 0.45600d-02
        hbs(14, 3, 3) = -0.44600d-02
        hbs(14, 4, 3) = -0.23920d-01
        hbs(14, 5, 3) = 0.13000d-03
        hbs(14, 6, 3) = 0.32800d-02
        hbs(14, 7, 3) = 0.98000d-03
        hbs(14, 8, 3) = 0.31630d-01
        hbs(14, 9, 3) = -0.29000d-02
        hbs(14,10, 3) = 0.29700d-02
        hbs(14,11, 3) = 0.43300d-02
        hbs(14,12, 3) = -0.18500d-02
        hbs(15, 1, 3) = -0.15810d-01
        hbs(15, 2, 3) = 0.31000d-02
        hbs(15, 3, 3) = -0.12600d-02
        hbs(15, 4, 3) = -0.19090d-01
        hbs(15, 5, 3) = 0.17000d-03
        hbs(15, 6, 3) = 0.20600d-02
        hbs(15, 7, 3) = -0.24000d-03
        hbs(15, 8, 3) = 0.41600d-02
        hbs(15, 9, 3) = -0.19700d-02
        hbs(15,10, 3) = 0.30180d-01
        hbs(15,11, 3) = 0.41100d-02
        hbs(15,12, 3) = -0.28100d-02
        hbs(16, 1, 3) = -0.56300d-02
        hbs(16, 2, 3) = -0.23750d-01
        hbs(16, 3, 3) = 0.15530d-01
        hbs(16, 4, 3) = 0.48210d-01
        hbs(16, 5, 3) = -0.15500d-02
        hbs(16, 6, 3) = 0.37160d-01
        hbs(16, 7, 3) = 0.26000d-03
        hbs(16, 8, 3) = -0.10920d-01
        hbs(16, 9, 3) = -0.64000d-03
        hbs(16,10, 3) = -0.35790d-01
        hbs(16,11, 3) = 0.47000d-02
        hbs(16,12, 3) = -0.80900d-02
        hbs(17, 1, 3) = 0.27600d-02
        hbs(17, 2, 3) = 0.21560d-01
        hbs(17, 3, 3) = -0.91700d-02
        hbs(17, 4, 3) = -0.25780d-01
        hbs(17, 5, 3) = 0.13700d-02
        hbs(17, 6, 3) = -0.16120d-01
        hbs(17, 7, 3) = -0.14000d-03
        hbs(17, 8, 3) = 0.52400d-02
        hbs(17, 9, 3) = -0.26000d-03
        hbs(17,10, 3) = 0.48900d-02
        hbs(17,11, 3) = -0.16800d-02
        hbs(17,12, 3) = 0.62500d-02
        hbs(18, 1, 3) = 0.28800d-02
        hbs(18, 2, 3) = 0.21900d-02
        hbs(18, 3, 3) = -0.63500d-02
        hbs(18, 4, 3) = -0.22430d-01
        hbs(18, 5, 3) = 0.18000d-03
        hbs(18, 6, 3) = -0.21030d-01
        hbs(18, 7, 3) = -0.12000d-03
        hbs(18, 8, 3) = 0.56800d-02
        hbs(18, 9, 3) = 0.91000d-03
        hbs(18,10, 3) = 0.30890d-01
        hbs(18,11, 3) = -0.30200d-02
        hbs(18,12, 3) = 0.18400d-02
        ! Hessian elements for anchor point 4
        hss( 1, 1, 4) = 0.30630d+00
        hss( 2, 1, 4) = -0.48000d-02
        hss( 2, 2, 4) = 0.30630d+00
        hss( 3, 1, 4) = 0.89520d-01
        hss( 3, 2, 4) = 0.89520d-01
        hss( 3, 3, 4) = 0.62529d+00
        hss( 4, 1, 4) = 0.88300d-02
        hss( 4, 2, 4) = 0.19200d-01
        hss( 4, 3, 4) = -0.16480d-01
        hss( 4, 4, 4) = 0.43562d+00
        hss( 5, 1, 4) = 0.21600d-02
        hss( 5, 2, 4) = 0.18000d-03
        hss( 5, 3, 4) = 0.10000d-03
        hss( 5, 4, 4) = 0.28600d-02
        hss( 5, 5, 4) = 0.36547d+00
        hss( 6, 1, 4) = 0.27000d-02
        hss( 6, 2, 4) = 0.49850d-01
        hss( 6, 3, 4) = 0.18220d-01
        hss( 6, 4, 4) = 0.43220d-01
        hss( 6, 5, 4) = -0.10000d-04
        hss( 6, 6, 4) = 0.35171d+00
        hss( 7, 1, 4) = -0.12200d-02
        hss( 7, 2, 4) = -0.39000d-03
        hss( 7, 3, 4) = 0.12600d-02
        hss( 7, 4, 4) = 0.51900d-02
        hss( 7, 5, 4) = 0.41000d-03
        hss( 7, 6, 4) = 0.38600d-02
        hss( 7, 7, 4) = 0.36379d+00
        hss( 8, 1, 4) = 0.49850d-01
        hss( 8, 2, 4) = 0.27000d-02
        hss( 8, 3, 4) = 0.18220d-01
        hss( 8, 4, 4) = -0.52300d-02
        hss( 8, 5, 4) = -0.11000d-03
        hss( 8, 6, 4) = 0.28880d-01
        hss( 8, 7, 4) = -0.34000d-03
        hss( 8, 8, 4) = 0.35171d+00
        hss( 9, 1, 4) = -0.29000d-03
        hss( 9, 2, 4) = -0.29000d-03
        hss( 9, 3, 4) = -0.38000d-03
        hss( 9, 4, 4) = -0.48000d-03
        hss( 9, 5, 4) = 0.80000d-04
        hss( 9, 6, 4) = 0.41700d-02
        hss( 9, 7, 4) = 0.40000d-03
        hss( 9, 8, 4) = 0.41700d-02
        hss( 9, 9, 4) = 0.36487d+00
        hss(10, 1, 4) = 0.19200d-01
        hss(10, 2, 4) = 0.88300d-02
        hss(10, 3, 4) = -0.16480d-01
        hss(10, 4, 4) = 0.80730d-01
        hss(10, 5, 4) = -0.72000d-03
        hss(10, 6, 4) = -0.52300d-02
        hss(10, 7, 4) = -0.55000d-03
        hss(10, 8, 4) = 0.43220d-01
        hss(10, 9, 4) = -0.48000d-03
        hss(10,10, 4) = 0.43562d+00
        hss(11, 1, 4) = -0.39000d-03
        hss(11, 2, 4) = -0.12200d-02
        hss(11, 3, 4) = 0.12600d-02
        hss(11, 4, 4) = -0.55000d-03
        hss(11, 5, 4) = 0.70000d-04
        hss(11, 6, 4) = -0.34000d-03
        hss(11, 7, 4) = 0.20000d-03
        hss(11, 8, 4) = 0.38600d-02
        hss(11, 9, 4) = 0.40000d-03
        hss(11,10, 4) = 0.51900d-02
        hss(11,11, 4) = 0.36379d+00
        hss(12, 1, 4) = 0.18000d-03
        hss(12, 2, 4) = 0.21600d-02
        hss(12, 3, 4) = 0.10000d-03
        hss(12, 4, 4) = -0.72000d-03
        hss(12, 5, 4) = 0.17000d-03
        hss(12, 6, 4) = -0.11000d-03
        hss(12, 7, 4) = 0.70000d-04
        hss(12, 8, 4) = -0.10000d-04
        hss(12, 9, 4) = 0.80000d-04
        hss(12,10, 4) = 0.28600d-02
        hss(12,11, 4) = 0.41000d-03
        hss(12,12, 4) = 0.36547d+00
        ! Hessian elements for anchor point 4
        hbb( 1, 1, 4) = 0.81570d-01
        hbb( 2, 1, 4) = -0.40780d-01
        hbb( 2, 2, 4) = 0.13357d+00
        hbb( 3, 1, 4) = -0.40780d-01
        hbb( 3, 2, 4) = -0.92780d-01
        hbb( 3, 3, 4) = 0.13357d+00
        hbb( 4, 1, 4) = -0.38240d-01
        hbb( 4, 2, 4) = 0.25780d-01
        hbb( 4, 3, 4) = 0.12470d-01
        hbb( 4, 4, 4) = 0.78290d-01
        hbb( 5, 1, 4) = 0.23810d-01
        hbb( 5, 2, 4) = -0.19980d-01
        hbb( 5, 3, 4) = -0.38300d-02
        hbb( 5, 4, 4) = -0.39070d-01
        hbb( 5, 5, 4) = 0.74370d-01
        hbb( 6, 1, 4) = 0.14430d-01
        hbb( 6, 2, 4) = -0.58000d-02
        hbb( 6, 3, 4) = -0.86300d-02
        hbb( 6, 4, 4) = -0.39220d-01
        hbb( 6, 5, 4) = -0.35300d-01
        hbb( 6, 6, 4) = 0.74520d-01
        hbb( 7, 1, 4) = -0.78000d-02
        hbb( 7, 2, 4) = -0.70900d-02
        hbb( 7, 3, 4) = 0.14890d-01
        hbb( 7, 4, 4) = -0.29810d-01
        hbb( 7, 5, 4) = 0.11250d-01
        hbb( 7, 6, 4) = 0.18570d-01
        hbb( 7, 7, 4) = 0.74890d-01
        hbb( 8, 1, 4) = -0.17500d-02
        hbb( 8, 2, 4) = 0.84900d-02
        hbb( 8, 3, 4) = -0.67400d-02
        hbb( 8, 4, 4) = 0.17840d-01
        hbb( 8, 5, 4) = -0.61200d-02
        hbb( 8, 6, 4) = -0.11720d-01
        hbb( 8, 7, 4) = -0.36950d-01
        hbb( 8, 8, 4) = 0.78020d-01
        hbb( 9, 1, 4) = 0.95600d-02
        hbb( 9, 2, 4) = -0.14000d-02
        hbb( 9, 3, 4) = -0.81500d-02
        hbb( 9, 4, 4) = 0.11970d-01
        hbb( 9, 5, 4) = -0.51200d-02
        hbb( 9, 6, 4) = -0.68500d-02
        hbb( 9, 7, 4) = -0.37940d-01
        hbb( 9, 8, 4) = -0.41060d-01
        hbb( 9, 9, 4) = 0.79000d-01
        hbb(10, 1, 4) = 0.10520d-01
        hbb(10, 2, 4) = -0.52600d-02
        hbb(10, 3, 4) = -0.52600d-02
        hbb(10, 4, 4) = -0.91700d-02
        hbb(10, 5, 4) = 0.10330d-01
        hbb(10, 6, 4) = -0.11700d-02
        hbb(10, 7, 4) = -0.34390d-01
        hbb(10, 8, 4) = 0.13830d-01
        hbb(10, 9, 4) = 0.20560d-01
        hbb(10,10, 4) = 0.76580d-01
        hbb(11, 1, 4) = -0.52600d-02
        hbb(11, 2, 4) = 0.31800d-02
        hbb(11, 3, 4) = 0.20800d-02
        hbb(11, 4, 4) = -0.11400d-02
        hbb(11, 5, 4) = -0.36300d-02
        hbb(11, 6, 4) = 0.47700d-02
        hbb(11, 7, 4) = 0.20190d-01
        hbb(11, 8, 4) = -0.74500d-02
        hbb(11, 9, 4) = -0.12740d-01
        hbb(11,10, 4) = -0.38290d-01
        hbb(11,11, 4) = 0.79780d-01
        hbb(12, 1, 4) = -0.52600d-02
        hbb(12, 2, 4) = 0.20800d-02
        hbb(12, 3, 4) = 0.31800d-02
        hbb(12, 4, 4) = 0.10300d-01
        hbb(12, 5, 4) = -0.67000d-02
        hbb(12, 6, 4) = -0.36000d-02
        hbb(12, 7, 4) = 0.14200d-01
        hbb(12, 8, 4) = -0.63700d-02
        hbb(12, 9, 4) = -0.78200d-02
        hbb(12,10, 4) = -0.38290d-01
        hbb(12,11, 4) = -0.41490d-01
        hbb(12,12, 4) = 0.79780d-01
        hbb(13, 1, 4) = -0.78000d-02
        hbb(13, 2, 4) = 0.14890d-01
        hbb(13, 3, 4) = -0.70900d-02
        hbb(13, 4, 4) = 0.37200d-02
        hbb(13, 5, 4) = -0.18400d-02
        hbb(13, 6, 4) = -0.18800d-02
        hbb(13, 7, 4) = -0.66100d-02
        hbb(13, 8, 4) = 0.88400d-02
        hbb(13, 9, 4) = -0.22300d-02
        hbb(13,10, 4) = -0.34390d-01
        hbb(13,11, 4) = 0.14200d-01
        hbb(13,12, 4) = 0.20190d-01
        hbb(13,13, 4) = 0.74890d-01
        hbb(14, 1, 4) = 0.95600d-02
        hbb(14, 2, 4) = -0.81500d-02
        hbb(14, 3, 4) = -0.14000d-02
        hbb(14, 4, 4) = -0.19200d-02
        hbb(14, 5, 4) = 0.79000d-03
        hbb(14, 6, 4) = 0.11300d-02
        hbb(14, 7, 4) = -0.22300d-02
        hbb(14, 8, 4) = -0.29700d-02
        hbb(14, 9, 4) = 0.52000d-02
        hbb(14,10, 4) = 0.20560d-01
        hbb(14,11, 4) = -0.78200d-02
        hbb(14,12, 4) = -0.12740d-01
        hbb(14,13, 4) = -0.37940d-01
        hbb(14,14, 4) = 0.79000d-01
        hbb(15, 1, 4) = -0.17500d-02
        hbb(15, 2, 4) = -0.67400d-02
        hbb(15, 3, 4) = 0.84900d-02
        hbb(15, 4, 4) = -0.18000d-02
        hbb(15, 5, 4) = 0.10500d-02
        hbb(15, 6, 4) = 0.75000d-03
        hbb(15, 7, 4) = 0.88400d-02
        hbb(15, 8, 4) = -0.58700d-02
        hbb(15, 9, 4) = -0.29700d-02
        hbb(15,10, 4) = 0.13830d-01
        hbb(15,11, 4) = -0.63700d-02
        hbb(15,12, 4) = -0.74500d-02
        hbb(15,13, 4) = -0.36950d-01
        hbb(15,14, 4) = -0.41060d-01
        hbb(15,15, 4) = 0.78020d-01
        hbb(16, 1, 4) = -0.38240d-01
        hbb(16, 2, 4) = 0.12470d-01
        hbb(16, 3, 4) = 0.25780d-01
        hbb(16, 4, 4) = -0.47900d-02
        hbb(16, 5, 4) = -0.44800d-02
        hbb(16, 6, 4) = 0.92700d-02
        hbb(16, 7, 4) = 0.37200d-02
        hbb(16, 8, 4) = -0.18000d-02
        hbb(16, 9, 4) = -0.19200d-02
        hbb(16,10, 4) = -0.91700d-02
        hbb(16,11, 4) = 0.10300d-01
        hbb(16,12, 4) = -0.11400d-02
        hbb(16,13, 4) = -0.29810d-01
        hbb(16,14, 4) = 0.11970d-01
        hbb(16,15, 4) = 0.17840d-01
        hbb(16,16, 4) = 0.78290d-01
        hbb(17, 1, 4) = 0.23810d-01
        hbb(17, 2, 4) = -0.38300d-02
        hbb(17, 3, 4) = -0.19980d-01
        hbb(17, 4, 4) = -0.44800d-02
        hbb(17, 5, 4) = 0.71700d-02
        hbb(17, 6, 4) = -0.26900d-02
        hbb(17, 7, 4) = -0.18400d-02
        hbb(17, 8, 4) = 0.10500d-02
        hbb(17, 9, 4) = 0.79000d-03
        hbb(17,10, 4) = 0.10330d-01
        hbb(17,11, 4) = -0.67000d-02
        hbb(17,12, 4) = -0.36300d-02
        hbb(17,13, 4) = 0.11250d-01
        hbb(17,14, 4) = -0.51200d-02
        hbb(17,15, 4) = -0.61200d-02
        hbb(17,16, 4) = -0.39070d-01
        hbb(17,17, 4) = 0.74370d-01
        hbb(18, 1, 4) = 0.14430d-01
        hbb(18, 2, 4) = -0.86300d-02
        hbb(18, 3, 4) = -0.58000d-02
        hbb(18, 4, 4) = 0.92700d-02
        hbb(18, 5, 4) = -0.26900d-02
        hbb(18, 6, 4) = -0.65700d-02
        hbb(18, 7, 4) = -0.18800d-02
        hbb(18, 8, 4) = 0.75000d-03
        hbb(18, 9, 4) = 0.11300d-02
        hbb(18,10, 4) = -0.11700d-02
        hbb(18,11, 4) = -0.36000d-02
        hbb(18,12, 4) = 0.47700d-02
        hbb(18,13, 4) = 0.18570d-01
        hbb(18,14, 4) = -0.68500d-02
        hbb(18,15, 4) = -0.11720d-01
        hbb(18,16, 4) = -0.39220d-01
        hbb(18,17, 4) = -0.35300d-01
        hbb(18,18, 4) = 0.74520d-01
        ! Hessian elements for anchor point 4
        htt( 1, 1, 4) = 0.86600d-02
        htt( 2, 1, 4) = 0.36400d-02
        htt( 2, 2, 4) = 0.78600d-02
        htt( 3, 1, 4) = -0.53000d-03
        htt( 3, 2, 4) = -0.36100d-02
        htt( 3, 3, 4) = 0.74400d-02
        htt( 4, 1, 4) = -0.55500d-02
        htt( 4, 2, 4) = 0.61000d-03
        htt( 4, 3, 4) = 0.43600d-02
        htt( 4, 4, 4) = 0.10520d-01
        htt( 5, 1, 4) = -0.62400d-02
        htt( 5, 2, 4) = -0.42900d-02
        htt( 5, 3, 4) = 0.29500d-02
        htt( 5, 4, 4) = 0.49100d-02
        htt( 5, 5, 4) = 0.86600d-02
        htt( 6, 1, 4) = -0.42900d-02
        htt( 6, 2, 4) = -0.29000d-02
        htt( 6, 3, 4) = 0.29600d-02
        htt( 6, 4, 4) = 0.43500d-02
        htt( 6, 5, 4) = 0.36400d-02
        htt( 6, 6, 4) = 0.78600d-02
        htt( 7, 1, 4) = 0.29500d-02
        htt( 7, 2, 4) = 0.29600d-02
        htt( 7, 3, 4) = -0.50200d-02
        htt( 7, 4, 4) = -0.50000d-02
        htt( 7, 5, 4) = -0.53000d-03
        htt( 7, 6, 4) = -0.36100d-02
        htt( 7, 7, 4) = 0.74400d-02
        htt( 8, 1, 4) = 0.49100d-02
        htt( 8, 2, 4) = 0.43500d-02
        htt( 8, 3, 4) = -0.50000d-02
        htt( 8, 4, 4) = -0.55600d-02
        htt( 8, 5, 4) = -0.55500d-02
        htt( 8, 6, 4) = 0.61000d-03
        htt( 8, 7, 4) = 0.43600d-02
        htt( 8, 8, 4) = 0.10520d-01
        htt( 9, 1, 4) = -0.58800d-02
        htt( 9, 2, 4) = -0.96000d-03
        htt( 9, 3, 4) = -0.18400d-02
        htt( 9, 4, 4) = 0.30900d-02
        htt( 9, 5, 4) = 0.10900d-02
        htt( 9, 6, 4) = 0.22400d-02
        htt( 9, 7, 4) = -0.29500d-02
        htt( 9, 8, 4) = -0.18000d-02
        htt( 9, 9, 4) = 0.80500d-02
        htt(10, 1, 4) = -0.39100d-02
        htt(10, 2, 4) = 0.46000d-03
        htt(10, 3, 4) = -0.55000d-03
        htt(10, 4, 4) = 0.38200d-02
        htt(10, 5, 4) = 0.23300d-02
        htt(10, 6, 4) = 0.17200d-02
        htt(10, 7, 4) = -0.10300d-02
        htt(10, 8, 4) = -0.16300d-02
        htt(10, 9, 4) = 0.28100d-02
        htt(10,10, 4) = 0.84200d-02
        htt(11, 1, 4) = -0.58000d-03
        htt(11, 2, 4) = -0.54100d-02
        htt(11, 3, 4) = 0.14200d-02
        htt(11, 4, 4) = -0.34100d-02
        htt(11, 5, 4) = -0.97000d-03
        htt(11, 6, 4) = 0.77000d-03
        htt(11, 7, 4) = -0.29700d-02
        htt(11, 8, 4) = -0.12200d-02
        htt(11, 9, 4) = 0.28600d-02
        htt(11,10, 4) = -0.18000d-02
        htt(11,11, 4) = 0.79500d-02
        htt(12, 1, 4) = 0.13900d-02
        htt(12, 2, 4) = -0.39900d-02
        htt(12, 3, 4) = 0.27100d-02
        htt(12, 4, 4) = -0.26800d-02
        htt(12, 5, 4) = 0.27000d-03
        htt(12, 6, 4) = 0.26000d-03
        htt(12, 7, 4) = -0.10500d-02
        htt(12, 8, 4) = -0.10500d-02
        htt(12, 9, 4) = -0.23900d-02
        htt(12,10, 4) = 0.38100d-02
        htt(12,11, 4) = 0.32900d-02
        htt(12,12, 4) = 0.94800d-02
        htt(13, 1, 4) = 0.42000d-03
        htt(13, 2, 4) = -0.12700d-02
        htt(13, 3, 4) = 0.19200d-02
        htt(13, 4, 4) = 0.23000d-03
        htt(13, 5, 4) = 0.20100d-02
        htt(13, 6, 4) = 0.61000d-03
        htt(13, 7, 4) = 0.51000d-03
        htt(13, 8, 4) = -0.88000d-03
        htt(13, 9, 4) = -0.54500d-02
        htt(13,10, 4) = -0.30000d-04
        htt(13,11, 4) = -0.36700d-02
        htt(13,12, 4) = 0.17500d-02
        htt(13,13, 4) = 0.82400d-02
        htt(14, 1, 4) = 0.18200d-02
        htt(14, 2, 4) = 0.50000d-03
        htt(14, 3, 4) = 0.11400d-02
        htt(14, 4, 4) = -0.18000d-03
        htt(14, 5, 4) = 0.61000d-03
        htt(14, 6, 4) = -0.11500d-02
        htt(14, 7, 4) = 0.12900d-02
        htt(14, 8, 4) = -0.47000d-03
        htt(14, 9, 4) = -0.31900d-02
        htt(14,10, 4) = 0.14500d-02
        htt(14,11, 4) = -0.18000d-02
        htt(14,12, 4) = 0.28400d-02
        htt(14,13, 4) = 0.21800d-02
        htt(14,14, 4) = 0.73800d-02
        htt(15, 1, 4) = -0.15300d-02
        htt(15, 2, 4) = -0.26700d-02
        htt(15, 3, 4) = 0.64000d-03
        htt(15, 4, 4) = -0.50000d-03
        htt(15, 5, 4) = 0.78000d-03
        htt(15, 6, 4) = 0.11300d-02
        htt(15, 7, 4) = -0.13900d-02
        htt(15, 8, 4) = -0.10500d-02
        htt(15, 9, 4) = -0.26000d-03
        htt(15,10, 4) = -0.55900d-02
        htt(15,11, 4) = 0.94000d-03
        htt(15,12, 4) = -0.43900d-02
        htt(15,13, 4) = 0.28700d-02
        htt(15,14, 4) = -0.24200d-02
        htt(15,15, 4) = 0.81600d-02
        htt(16, 1, 4) = -0.14000d-03
        htt(16, 2, 4) = -0.90000d-03
        htt(16, 3, 4) = -0.13000d-03
        htt(16, 4, 4) = -0.90000d-03
        htt(16, 5, 4) = -0.62000d-03
        htt(16, 6, 4) = -0.64000d-03
        htt(16, 7, 4) = -0.62000d-03
        htt(16, 8, 4) = -0.64000d-03
        htt(16, 9, 4) = 0.20000d-02
        htt(16,10, 4) = -0.41100d-02
        htt(16,11, 4) = 0.28100d-02
        htt(16,12, 4) = -0.32900d-02
        htt(16,13, 4) = -0.31800d-02
        htt(16,14, 4) = 0.27800d-02
        htt(16,15, 4) = 0.28700d-02
        htt(16,16, 4) = 0.88300d-02
        htt(17, 1, 4) = 0.20100d-02
        htt(17, 2, 4) = 0.61000d-03
        htt(17, 3, 4) = 0.51000d-03
        htt(17, 4, 4) = -0.88000d-03
        htt(17, 5, 4) = 0.42000d-03
        htt(17, 6, 4) = -0.12700d-02
        htt(17, 7, 4) = 0.19200d-02
        htt(17, 8, 4) = 0.23000d-03
        htt(17, 9, 4) = 0.64000d-03
        htt(17,10, 4) = -0.15500d-02
        htt(17,11, 4) = 0.21100d-02
        htt(17,12, 4) = -0.80000d-04
        htt(17,13, 4) = -0.58000d-02
        htt(17,14, 4) = 0.26000d-03
        htt(17,15, 4) = -0.36300d-02
        htt(17,16, 4) = 0.24300d-02
        htt(17,17, 4) = 0.82400d-02
        htt(18, 1, 4) = 0.78000d-03
        htt(18, 2, 4) = 0.11300d-02
        htt(18, 3, 4) = -0.13900d-02
        htt(18, 4, 4) = -0.10500d-02
        htt(18, 5, 4) = -0.15300d-02
        htt(18, 6, 4) = -0.26700d-02
        htt(18, 7, 4) = 0.64000d-03
        htt(18, 8, 4) = -0.50000d-03
        htt(18, 9, 4) = 0.17400d-02
        htt(18,10, 4) = 0.40000d-03
        htt(18,11, 4) = 0.13800d-02
        htt(18,12, 4) = 0.40000d-04
        htt(18,13, 4) = -0.36300d-02
        htt(18,14, 4) = 0.16600d-02
        htt(18,15, 4) = -0.23000d-02
        htt(18,16, 4) = 0.29900d-02
        htt(18,17, 4) = 0.28700d-02
        htt(18,18, 4) = 0.81600d-02
        htt(19, 1, 4) = 0.61000d-03
        htt(19, 2, 4) = -0.11500d-02
        htt(19, 3, 4) = 0.12900d-02
        htt(19, 4, 4) = -0.47000d-03
        htt(19, 5, 4) = 0.18200d-02
        htt(19, 6, 4) = 0.50000d-03
        htt(19, 7, 4) = 0.11400d-02
        htt(19, 8, 4) = -0.18000d-03
        htt(19, 9, 4) = -0.16200d-02
        htt(19,10, 4) = -0.30300d-02
        htt(19,11, 4) = 0.24000d-03
        htt(19,12, 4) = -0.11700d-02
        htt(19,13, 4) = 0.26000d-03
        htt(19,14, 4) = -0.49400d-02
        htt(19,15, 4) = 0.16600d-02
        htt(19,16, 4) = -0.35300d-02
        htt(19,17, 4) = 0.21800d-02
        htt(19,18, 4) = -0.24200d-02
        htt(19,19, 4) = 0.73800d-02
        htt(20, 1, 4) = -0.62000d-03
        htt(20, 2, 4) = -0.64000d-03
        htt(20, 3, 4) = -0.62000d-03
        htt(20, 4, 4) = -0.64000d-03
        htt(20, 5, 4) = -0.14000d-03
        htt(20, 6, 4) = -0.90000d-03
        htt(20, 7, 4) = -0.13000d-03
        htt(20, 8, 4) = -0.90000d-03
        htt(20, 9, 4) = -0.52000d-03
        htt(20,10, 4) = -0.10800d-02
        htt(20,11, 4) = -0.49000d-03
        htt(20,12, 4) = -0.10500d-02
        htt(20,13, 4) = 0.24300d-02
        htt(20,14, 4) = -0.35300d-02
        htt(20,15, 4) = 0.29900d-02
        htt(20,16, 4) = -0.29700d-02
        htt(20,17, 4) = -0.31800d-02
        htt(20,18, 4) = 0.28700d-02
        htt(20,19, 4) = 0.27800d-02
        htt(20,20, 4) = 0.88300d-02
        htt(21, 1, 4) = 0.10900d-02
        htt(21, 2, 4) = 0.22400d-02
        htt(21, 3, 4) = -0.29500d-02
        htt(21, 4, 4) = -0.18000d-02
        htt(21, 5, 4) = -0.58800d-02
        htt(21, 6, 4) = -0.96000d-03
        htt(21, 7, 4) = -0.18400d-02
        htt(21, 8, 4) = 0.30900d-02
        htt(21, 9, 4) = 0.14200d-02
        htt(21,10, 4) = 0.31000d-03
        htt(21,11, 4) = 0.21000d-03
        htt(21,12, 4) = -0.90000d-03
        htt(21,13, 4) = 0.64000d-03
        htt(21,14, 4) = -0.16200d-02
        htt(21,15, 4) = 0.17400d-02
        htt(21,16, 4) = -0.52000d-03
        htt(21,17, 4) = -0.54500d-02
        htt(21,18, 4) = -0.26000d-03
        htt(21,19, 4) = -0.31900d-02
        htt(21,20, 4) = 0.20000d-02
        htt(21,21, 4) = 0.80500d-02
        htt(22, 1, 4) = -0.97000d-03
        htt(22, 2, 4) = 0.77000d-03
        htt(22, 3, 4) = -0.29700d-02
        htt(22, 4, 4) = -0.12200d-02
        htt(22, 5, 4) = -0.58000d-03
        htt(22, 6, 4) = -0.54100d-02
        htt(22, 7, 4) = 0.14200d-02
        htt(22, 8, 4) = -0.34100d-02
        htt(22, 9, 4) = 0.21000d-03
        htt(22,10, 4) = 0.95000d-03
        htt(22,11, 4) = -0.16300d-02
        htt(22,12, 4) = -0.89000d-03
        htt(22,13, 4) = 0.21100d-02
        htt(22,14, 4) = 0.24000d-03
        htt(22,15, 4) = 0.13800d-02
        htt(22,16, 4) = -0.49000d-03
        htt(22,17, 4) = -0.36700d-02
        htt(22,18, 4) = 0.94000d-03
        htt(22,19, 4) = -0.18000d-02
        htt(22,20, 4) = 0.28100d-02
        htt(22,21, 4) = 0.28600d-02
        htt(22,22, 4) = 0.79500d-02
        htt(23, 1, 4) = 0.23300d-02
        htt(23, 2, 4) = 0.17200d-02
        htt(23, 3, 4) = -0.10300d-02
        htt(23, 4, 4) = -0.16300d-02
        htt(23, 5, 4) = -0.39100d-02
        htt(23, 6, 4) = 0.46000d-03
        htt(23, 7, 4) = -0.55000d-03
        htt(23, 8, 4) = 0.38200d-02
        htt(23, 9, 4) = 0.31000d-03
        htt(23,10, 4) = -0.16600d-02
        htt(23,11, 4) = 0.95000d-03
        htt(23,12, 4) = -0.10200d-02
        htt(23,13, 4) = -0.15500d-02
        htt(23,14, 4) = -0.30300d-02
        htt(23,15, 4) = 0.40000d-03
        htt(23,16, 4) = -0.10800d-02
        htt(23,17, 4) = -0.30000d-04
        htt(23,18, 4) = -0.55900d-02
        htt(23,19, 4) = 0.14500d-02
        htt(23,20, 4) = -0.41100d-02
        htt(23,21, 4) = 0.28100d-02
        htt(23,22, 4) = -0.18000d-02
        htt(23,23, 4) = 0.84200d-02
        htt(24, 1, 4) = 0.27000d-03
        htt(24, 2, 4) = 0.26000d-03
        htt(24, 3, 4) = -0.10500d-02
        htt(24, 4, 4) = -0.10500d-02
        htt(24, 5, 4) = 0.13900d-02
        htt(24, 6, 4) = -0.39900d-02
        htt(24, 7, 4) = 0.27100d-02
        htt(24, 8, 4) = -0.26800d-02
        htt(24, 9, 4) = -0.90000d-03
        htt(24,10, 4) = -0.10200d-02
        htt(24,11, 4) = -0.89000d-03
        htt(24,12, 4) = -0.10100d-02
        htt(24,13, 4) = -0.80000d-04
        htt(24,14, 4) = -0.11700d-02
        htt(24,15, 4) = 0.40000d-04
        htt(24,16, 4) = -0.10500d-02
        htt(24,17, 4) = 0.17500d-02
        htt(24,18, 4) = -0.43900d-02
        htt(24,19, 4) = 0.28400d-02
        htt(24,20, 4) = -0.32900d-02
        htt(24,21, 4) = -0.23900d-02
        htt(24,22, 4) = 0.32900d-02
        htt(24,23, 4) = 0.38100d-02
        htt(24,24, 4) = 0.94800d-02
        ! Hessian elements for anchor point 4
        hbs( 1, 1, 4) = -0.22700d-01
        hbs( 1, 2, 4) = -0.22700d-01
        hbs( 1, 3, 4) = -0.46990d-01
        hbs( 1, 4, 4) = -0.88700d-02
        hbs( 1, 5, 4) = 0.54300d-02
        hbs( 1, 6, 4) = 0.44400d-01
        hbs( 1, 7, 4) = -0.55000d-03
        hbs( 1, 8, 4) = 0.44400d-01
        hbs( 1, 9, 4) = -0.10000d-04
        hbs( 1,10, 4) = -0.88700d-02
        hbs( 1,11, 4) = -0.55000d-03
        hbs( 1,12, 4) = 0.54300d-02
        hbs( 2, 1, 4) = 0.32950d-01
        hbs( 2, 2, 4) = -0.10250d-01
        hbs( 2, 3, 4) = 0.23500d-01
        hbs( 2, 4, 4) = 0.60000d-02
        hbs( 2, 5, 4) = -0.49200d-02
        hbs( 2, 6, 4) = -0.28130d-01
        hbs( 2, 7, 4) = 0.12500d-02
        hbs( 2, 8, 4) = -0.16270d-01
        hbs( 2, 9, 4) = 0.00000d+00
        hbs( 2,10, 4) = 0.28700d-02
        hbs( 2,11, 4) = -0.70000d-03
        hbs( 2,12, 4) = -0.51000d-03
        hbs( 3, 1, 4) = -0.10250d-01
        hbs( 3, 2, 4) = 0.32950d-01
        hbs( 3, 3, 4) = 0.23500d-01
        hbs( 3, 4, 4) = 0.28700d-02
        hbs( 3, 5, 4) = -0.51000d-03
        hbs( 3, 6, 4) = -0.16270d-01
        hbs( 3, 7, 4) = -0.70000d-03
        hbs( 3, 8, 4) = -0.28130d-01
        hbs( 3, 9, 4) = 0.00000d+00
        hbs( 3,10, 4) = 0.60000d-02
        hbs( 3,11, 4) = 0.12500d-02
        hbs( 3,12, 4) = -0.49200d-02
        hbs( 4, 1, 4) = -0.24180d-01
        hbs( 4, 2, 4) = -0.60700d-02
        hbs( 4, 3, 4) = 0.14620d-01
        hbs( 4, 4, 4) = -0.35760d-01
        hbs( 4, 5, 4) = -0.81300d-02
        hbs( 4, 6, 4) = -0.10830d-01
        hbs( 4, 7, 4) = 0.47400d-02
        hbs( 4, 8, 4) = 0.37380d-01
        hbs( 4, 9, 4) = -0.65000d-03
        hbs( 4,10, 4) = 0.48130d-01
        hbs( 4,11, 4) = 0.25000d-03
        hbs( 4,12, 4) = -0.12600d-02
        hbs( 5, 1, 4) = 0.21630d-01
        hbs( 5, 2, 4) = 0.30300d-02
        hbs( 5, 3, 4) = -0.86900d-02
        hbs( 5, 4, 4) = 0.48300d-02
        hbs( 5, 5, 4) = 0.60700d-02
        hbs( 5, 6, 4) = 0.52000d-02
        hbs( 5, 7, 4) = -0.17000d-02
        hbs( 5, 8, 4) = -0.16220d-01
        hbs( 5, 9, 4) = -0.27000d-03
        hbs( 5,10, 4) = -0.25720d-01
        hbs( 5,11, 4) = -0.13000d-03
        hbs( 5,12, 4) = 0.12100d-02
        hbs( 6, 1, 4) = 0.25500d-02
        hbs( 6, 2, 4) = 0.30300d-02
        hbs( 6, 3, 4) = -0.59400d-02
        hbs( 6, 4, 4) = 0.30930d-01
        hbs( 6, 5, 4) = 0.20600d-02
        hbs( 6, 6, 4) = 0.56300d-02
        hbs( 6, 7, 4) = -0.30400d-02
        hbs( 6, 8, 4) = -0.21160d-01
        hbs( 6, 9, 4) = 0.91000d-03
        hbs( 6,10, 4) = -0.22410d-01
        hbs( 6,11, 4) = -0.12000d-03
        hbs( 6,12, 4) = 0.50000d-04
        hbs( 7, 1, 4) = -0.78200d-02
        hbs( 7, 2, 4) = 0.29170d-01
        hbs( 7, 3, 4) = 0.57100d-02
        hbs( 7, 4, 4) = -0.32990d-01
        hbs( 7, 5, 4) = 0.46200d-02
        hbs( 7, 6, 4) = -0.35760d-01
        hbs( 7, 7, 4) = -0.84700d-02
        hbs( 7, 8, 4) = -0.55800d-02
        hbs( 7, 9, 4) = 0.48900d-02
        hbs( 7,10, 4) = 0.43210d-01
        hbs( 7,11, 4) = -0.75000d-03
        hbs( 7,12, 4) = -0.32000d-03
        hbs( 8, 1, 4) = 0.31800d-02
        hbs( 8, 2, 4) = -0.15950d-01
        hbs( 8, 3, 4) = -0.12300d-02
        hbs( 8, 4, 4) = 0.30010d-01
        hbs( 8, 5, 4) = -0.27700d-02
        hbs( 8, 6, 4) = 0.41300d-02
        hbs( 8, 7, 4) = 0.40300d-02
        hbs( 8, 8, 4) = 0.21600d-02
        hbs( 8, 9, 4) = -0.19900d-02
        hbs( 8,10, 4) = -0.19160d-01
        hbs( 8,11, 4) = -0.23000d-03
        hbs( 8,12, 4) = 0.17000d-03
        hbs( 9, 1, 4) = 0.46400d-02
        hbs( 9, 2, 4) = -0.13220d-01
        hbs( 9, 3, 4) = -0.44800d-02
        hbs( 9, 4, 4) = 0.29800d-02
        hbs( 9, 5, 4) = -0.18500d-02
        hbs( 9, 6, 4) = 0.31630d-01
        hbs( 9, 7, 4) = 0.44300d-02
        hbs( 9, 8, 4) = 0.34300d-02
        hbs( 9, 9, 4) = -0.29100d-02
        hbs( 9,10, 4) = -0.24050d-01
        hbs( 9,11, 4) = 0.98000d-03
        hbs( 9,12, 4) = 0.15000d-03
        hbs(10, 1, 4) = 0.31600d-01
        hbs(10, 2, 4) = 0.31600d-01
        hbs(10, 3, 4) = 0.63400d-02
        hbs(10, 4, 4) = -0.13710d-01
        hbs(10, 5, 4) = -0.34000d-03
        hbs(10, 6, 4) = -0.29610d-01
        hbs(10, 7, 4) = 0.47800d-02
        hbs(10, 8, 4) = -0.29610d-01
        hbs(10, 9, 4) = -0.84900d-02
        hbs(10,10, 4) = -0.13710d-01
        hbs(10,11, 4) = 0.47800d-02
        hbs(10,12, 4) = -0.34000d-03
        hbs(11, 1, 4) = -0.17850d-01
        hbs(11, 2, 4) = -0.13750d-01
        hbs(11, 3, 4) = -0.31700d-02
        hbs(11, 4, 4) = 0.80700d-02
        hbs(11, 5, 4) = 0.73000d-03
        hbs(11, 6, 4) = 0.27610d-01
        hbs(11, 7, 4) = -0.29000d-02
        hbs(11, 8, 4) = 0.20000d-02
        hbs(11, 9, 4) = 0.42400d-02
        hbs(11,10, 4) = 0.56400d-02
        hbs(11,11, 4) = -0.18700d-02
        hbs(11,12, 4) = -0.38000d-03
        hbs(12, 1, 4) = -0.13750d-01
        hbs(12, 2, 4) = -0.17850d-01
        hbs(12, 3, 4) = -0.31700d-02
        hbs(12, 4, 4) = 0.56400d-02
        hbs(12, 5, 4) = -0.38000d-03
        hbs(12, 6, 4) = 0.20000d-02
        hbs(12, 7, 4) = -0.18700d-02
        hbs(12, 8, 4) = 0.27610d-01
        hbs(12, 9, 4) = 0.42400d-02
        hbs(12,10, 4) = 0.80700d-02
        hbs(12,11, 4) = -0.29000d-02
        hbs(12,12, 4) = 0.73000d-03
        hbs(13, 1, 4) = 0.29170d-01
        hbs(13, 2, 4) = -0.78200d-02
        hbs(13, 3, 4) = 0.57100d-02
        hbs(13, 4, 4) = 0.43210d-01
        hbs(13, 5, 4) = -0.32000d-03
        hbs(13, 6, 4) = -0.55800d-02
        hbs(13, 7, 4) = -0.75000d-03
        hbs(13, 8, 4) = -0.35760d-01
        hbs(13, 9, 4) = 0.48900d-02
        hbs(13,10, 4) = -0.32990d-01
        hbs(13,11, 4) = -0.84700d-02
        hbs(13,12, 4) = 0.46200d-02
        hbs(14, 1, 4) = -0.13220d-01
        hbs(14, 2, 4) = 0.46400d-02
        hbs(14, 3, 4) = -0.44800d-02
        hbs(14, 4, 4) = -0.24050d-01
        hbs(14, 5, 4) = 0.15000d-03
        hbs(14, 6, 4) = 0.34300d-02
        hbs(14, 7, 4) = 0.98000d-03
        hbs(14, 8, 4) = 0.31630d-01
        hbs(14, 9, 4) = -0.29100d-02
        hbs(14,10, 4) = 0.29800d-02
        hbs(14,11, 4) = 0.44300d-02
        hbs(14,12, 4) = -0.18500d-02
        hbs(15, 1, 4) = -0.15950d-01
        hbs(15, 2, 4) = 0.31800d-02
        hbs(15, 3, 4) = -0.12300d-02
        hbs(15, 4, 4) = -0.19160d-01
        hbs(15, 5, 4) = 0.17000d-03
        hbs(15, 6, 4) = 0.21600d-02
        hbs(15, 7, 4) = -0.23000d-03
        hbs(15, 8, 4) = 0.41300d-02
        hbs(15, 9, 4) = -0.19900d-02
        hbs(15,10, 4) = 0.30010d-01
        hbs(15,11, 4) = 0.40300d-02
        hbs(15,12, 4) = -0.27700d-02
        hbs(16, 1, 4) = -0.60700d-02
        hbs(16, 2, 4) = -0.24180d-01
        hbs(16, 3, 4) = 0.14620d-01
        hbs(16, 4, 4) = 0.48130d-01
        hbs(16, 5, 4) = -0.12600d-02
        hbs(16, 6, 4) = 0.37380d-01
        hbs(16, 7, 4) = 0.25000d-03
        hbs(16, 8, 4) = -0.10830d-01
        hbs(16, 9, 4) = -0.65000d-03
        hbs(16,10, 4) = -0.35760d-01
        hbs(16,11, 4) = 0.47400d-02
        hbs(16,12, 4) = -0.81300d-02
        hbs(17, 1, 4) = 0.30300d-02
        hbs(17, 2, 4) = 0.21630d-01
        hbs(17, 3, 4) = -0.86900d-02
        hbs(17, 4, 4) = -0.25720d-01
        hbs(17, 5, 4) = 0.12100d-02
        hbs(17, 6, 4) = -0.16220d-01
        hbs(17, 7, 4) = -0.13000d-03
        hbs(17, 8, 4) = 0.52000d-02
        hbs(17, 9, 4) = -0.27000d-03
        hbs(17,10, 4) = 0.48300d-02
        hbs(17,11, 4) = -0.17000d-02
        hbs(17,12, 4) = 0.60700d-02
        hbs(18, 1, 4) = 0.30300d-02
        hbs(18, 2, 4) = 0.25500d-02
        hbs(18, 3, 4) = -0.59400d-02
        hbs(18, 4, 4) = -0.22410d-01
        hbs(18, 5, 4) = 0.50000d-04
        hbs(18, 6, 4) = -0.21160d-01
        hbs(18, 7, 4) = -0.12000d-03
        hbs(18, 8, 4) = 0.56300d-02
        hbs(18, 9, 4) = 0.91000d-03
        hbs(18,10, 4) = 0.30930d-01
        hbs(18,11, 4) = -0.30400d-02
        hbs(18,12, 4) = 0.20600d-02
    endif

end subroutine gethess
!> @brief Evaluate getgeom for the PHOH_APRP_DPEM potential.
subroutine getgeom( nap, qse, nqs, qbe, nqb, qte, nqt, ii )

    implicit none

    double precision, parameter :: zero = 0.0d0
    double precision, parameter :: pi   = 3.1415926535898d0

    ! Global variables
    integer          :: ii
    integer          :: nap
    integer          :: nqs,nqb,nqt
    double precision :: qse(nqs,nap),qbe(nqb,nap),qte(nqt,nap)

    ! Local variables
    integer :: i,j

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

    do i=2,nap
        do j=1,nqt
            qte(j,i) = qte(j,1)
        enddo
    enddo

    ! Now get equilibrium geometrical parameters for qae
    if (ii .eq. 11) then
        ! equilibrium geometry at anchor point 1
        qse( 1, 1) = 0.263789d+01
        qse( 2, 1) = 0.263870d+01
        qse( 3, 1) = 0.258038d+01
        qse( 4, 1) = 0.263573d+01
        qse( 5, 1) = 0.206310d+01
        qse( 6, 1) = 0.263336d+01
        qse( 7, 1) = 0.205963d+01
        qse( 8, 1) = 0.263920d+01
        qse( 9, 1) = 0.205742d+01
        qse(10, 1) = 0.262993d+01
        qse(11, 1) = 0.205975d+01
        qse(12, 1) = 0.205790d+01
        ! equilibrium geometry at anchor point 1
        qbe( 1, 1) = 0.210140d+01
        qbe( 2, 1) = 0.213715d+01
        qbe( 3, 1) = 0.204464d+01
        qbe( 4, 1) = 0.208658d+01
        qbe( 5, 1) = 0.209471d+01
        qbe( 6, 1) = 0.210190d+01
        qbe( 7, 1) = 0.210459d+01
        qbe( 8, 1) = 0.208211d+01
        qbe( 9, 1) = 0.209649d+01
        qbe(10, 1) = 0.208155d+01
        qbe(11, 1) = 0.210036d+01
        qbe(12, 1) = 0.210127d+01
        qbe(13, 1) = 0.210834d+01
        qbe(14, 1) = 0.209320d+01
        qbe(15, 1) = 0.208164d+01
        qbe(16, 1) = 0.208392d+01
        qbe(17, 1) = 0.207495d+01
        qbe(18, 1) = 0.212431d+01
        ! equilibrium geometry at anchor point 2
        qse( 1, 2) = 0.263724d+01
        qse( 2, 2) = 0.264694d+01
        qse( 3, 2) = 0.257138d+01
        qse( 4, 2) = 0.263516d+01
        qse( 5, 2) = 0.206589d+01
        qse( 6, 2) = 0.263264d+01
        qse( 7, 2) = 0.205963d+01
        qse( 8, 2) = 0.263972d+01
        qse( 9, 2) = 0.205726d+01
        qse(10, 2) = 0.262912d+01
        qse(11, 2) = 0.205993d+01
        qse(12, 2) = 0.205801d+01
        ! equilibrium geometry at anchor point 2
        qbe( 1, 2) = 0.209269d+01
        qbe( 2, 2) = 0.215700d+01
        qbe( 3, 2) = 0.203350d+01
        qbe( 4, 2) = 0.209468d+01
        qbe( 5, 2) = 0.208443d+01
        qbe( 6, 2) = 0.210408d+01
        qbe( 7, 2) = 0.210277d+01
        qbe( 8, 2) = 0.208276d+01
        qbe( 9, 2) = 0.209766d+01
        qbe(10, 2) = 0.207979d+01
        qbe(11, 2) = 0.210131d+01
        qbe(12, 2) = 0.210209d+01
        qbe(13, 2) = 0.211103d+01
        qbe(14, 2) = 0.209320d+01
        qbe(15, 2) = 0.207895d+01
        qbe(16, 2) = 0.208541d+01
        qbe(17, 2) = 0.207654d+01
        qbe(18, 2) = 0.212123d+01
        ! equilibrium geometry at anchor point 3
        qse( 1, 3) = 0.264670d+01
        qse( 2, 3) = 0.266084d+01
        qse( 3, 3) = 0.252665d+01
        qse( 4, 3) = 0.263224d+01
        qse( 5, 3) = 0.206010d+01
        qse( 6, 3) = 0.263214d+01
        qse( 7, 3) = 0.205982d+01
        qse( 8, 3) = 0.263685d+01
        qse( 9, 3) = 0.205662d+01
        qse(10, 3) = 0.262787d+01
        qse(11, 3) = 0.206007d+01
        qse(12, 3) = 0.205634d+01
        ! equilibrium geometry at anchor point 3
        qbe( 1, 3) = 0.209445d+01
        qbe( 2, 3) = 0.215439d+01
        qbe( 3, 3) = 0.203434d+01
        qbe( 4, 3) = 0.208659d+01
        qbe( 5, 3) = 0.207835d+01
        qbe( 6, 3) = 0.211824d+01
        qbe( 7, 3) = 0.211238d+01
        qbe( 8, 3) = 0.207426d+01
        qbe( 9, 3) = 0.209655d+01
        qbe(10, 3) = 0.207613d+01
        qbe(11, 3) = 0.210353d+01
        qbe(12, 3) = 0.210352d+01
        qbe(13, 3) = 0.211662d+01
        qbe(14, 3) = 0.209569d+01
        qbe(15, 3) = 0.207088d+01
        qbe(16, 3) = 0.208020d+01
        qbe(17, 3) = 0.208110d+01
        qbe(18, 3) = 0.212189d+01
        ! equilibrium geometry at anchor point 4
        qse( 1, 4) = 0.265159d+01
        qse( 2, 4) = 0.265159d+01
        qse( 3, 4) = 0.251838d+01
        qse( 4, 4) = 0.263062d+01
        qse( 5, 4) = 0.205598d+01
        qse( 6, 4) = 0.263458d+01
        qse( 7, 4) = 0.205967d+01
        qse( 8, 4) = 0.263458d+01
        qse( 9, 4) = 0.205648d+01
        qse(10, 4) = 0.263062d+01
        qse(11, 4) = 0.205967d+01
        qse(12, 4) = 0.205598d+01
        ! equilibrium geometry at anchor point 4
        qbe( 1, 4) = 0.211210d+01
        qbe( 2, 4) = 0.208554d+01
        qbe( 3, 4) = 0.208554d+01
        qbe( 4, 4) = 0.207166d+01
        qbe( 5, 4) = 0.208711d+01
        qbe( 6, 4) = 0.212441d+01
        qbe( 7, 4) = 0.211593d+01
        qbe( 8, 4) = 0.207134d+01
        qbe( 9, 4) = 0.209592d+01
        qbe(10, 4) = 0.207909d+01
        qbe(11, 4) = 0.210205d+01
        qbe(12, 4) = 0.210205d+01
        qbe(13, 4) = 0.211593d+01
        qbe(14, 4) = 0.209592d+01
        qbe(15, 4) = 0.207134d+01
        qbe(16, 4) = 0.207166d+01
        qbe(17, 4) = 0.208711d+01
        qbe(18, 4) = 0.212441d+01
    else if (ii .eq. 22) then
        ! equilibrium geometry at anchor point 1
        qse( 1, 1) = 0.269941d+01
        qse( 2, 1) = 0.268363d+01
        qse( 3, 1) = 0.252916d+01
        qse( 4, 1) = 0.268093d+01
        qse( 5, 1) = 0.205887d+01
        qse( 6, 1) = 0.268592d+01
        qse( 7, 1) = 0.205454d+01
        qse( 8, 1) = 0.268227d+01
        qse( 9, 1) = 0.206068d+01
        qse(10, 1) = 0.268469d+01
        qse(11, 1) = 0.205375d+01
        qse(12, 1) = 0.205424d+01
        ! equilibrium geometry at anchor point 1
        qbe( 1, 1) = 0.216375d+01
        qbe( 2, 1) = 0.209466d+01
        qbe( 3, 1) = 0.202477d+01
        qbe( 4, 1) = 0.204773d+01
        qbe( 5, 1) = 0.209602d+01
        qbe( 6, 1) = 0.213944d+01
        qbe( 7, 1) = 0.208156d+01
        qbe( 8, 1) = 0.210850d+01
        qbe( 9, 1) = 0.209313d+01
        qbe(10, 1) = 0.214184d+01
        qbe(11, 1) = 0.206797d+01
        qbe(12, 1) = 0.207337d+01
        qbe(13, 1) = 0.207509d+01
        qbe(14, 1) = 0.209807d+01
        qbe(15, 1) = 0.211003d+01
        qbe(16, 1) = 0.205641d+01
        qbe(17, 1) = 0.207231d+01
        qbe(18, 1) = 0.215446d+01
        ! equilibrium geometry at anchor point 2
        qse( 1, 2) = 0.269988d+01
        qse( 2, 2) = 0.269065d+01
        qse( 3, 2) = 0.251103d+01
        qse( 4, 2) = 0.267864d+01
        qse( 5, 2) = 0.206339d+01
        qse( 6, 2) = 0.268396d+01
        qse( 7, 2) = 0.205439d+01
        qse( 8, 2) = 0.268024d+01
        qse( 9, 2) = 0.206089d+01
        qse(10, 2) = 0.268483d+01
        qse(11, 2) = 0.205390d+01
        qse(12, 2) = 0.205412d+01
        ! equilibrium geometry at anchor point 2
        qbe( 1, 2) = 0.215707d+01
        qbe( 2, 2) = 0.211094d+01
        qbe( 3, 2) = 0.201518d+01
        qbe( 4, 2) = 0.205661d+01
        qbe( 5, 2) = 0.208168d+01
        qbe( 6, 2) = 0.214490d+01
        qbe( 7, 2) = 0.207623d+01
        qbe( 8, 2) = 0.211213d+01
        qbe( 9, 2) = 0.209483d+01
        qbe(10, 2) = 0.214363d+01
        qbe(11, 2) = 0.206707d+01
        qbe(12, 2) = 0.207249d+01
        qbe(13, 2) = 0.207837d+01
        qbe(14, 2) = 0.209780d+01
        qbe(15, 2) = 0.210701d+01
        qbe(16, 2) = 0.205447d+01
        qbe(17, 2) = 0.207468d+01
        qbe(18, 2) = 0.215403d+01
        ! equilibrium geometry at anchor point 3
        qse( 1, 3) = 0.269988d+01
        qse( 2, 3) = 0.269065d+01
        qse( 3, 3) = 0.251103d+01
        qse( 4, 3) = 0.267864d+01
        qse( 5, 3) = 0.206339d+01
        qse( 6, 3) = 0.268396d+01
        qse( 7, 3) = 0.205439d+01
        qse( 8, 3) = 0.268024d+01
        qse( 9, 3) = 0.206089d+01
        qse(10, 3) = 0.268483d+01
        qse(11, 3) = 0.205390d+01
        qse(12, 3) = 0.205412d+01
        ! equilibrium geometry at anchor point 3
        qbe( 1, 3) = 0.215707d+01
        qbe( 2, 3) = 0.211094d+01
        qbe( 3, 3) = 0.201518d+01
        qbe( 4, 3) = 0.205661d+01
        qbe( 5, 3) = 0.208168d+01
        qbe( 6, 3) = 0.214490d+01
        qbe( 7, 3) = 0.207623d+01
        qbe( 8, 3) = 0.211213d+01
        qbe( 9, 3) = 0.209483d+01
        qbe(10, 3) = 0.214363d+01
        qbe(11, 3) = 0.206707d+01
        qbe(12, 3) = 0.207249d+01
        qbe(13, 3) = 0.207837d+01
        qbe(14, 3) = 0.209780d+01
        qbe(15, 3) = 0.210701d+01
        qbe(16, 3) = 0.205447d+01
        qbe(17, 3) = 0.207468d+01
        qbe(18, 3) = 0.215403d+01
        ! equilibrium geometry at anchor point 4
        qse( 1, 4) = 0.269988d+01
        qse( 2, 4) = 0.269065d+01
        qse( 3, 4) = 0.251103d+01
        qse( 4, 4) = 0.267864d+01
        qse( 5, 4) = 0.206339d+01
        qse( 6, 4) = 0.268396d+01
        qse( 7, 4) = 0.205439d+01
        qse( 8, 4) = 0.268024d+01
        qse( 9, 4) = 0.206089d+01
        qse(10, 4) = 0.268483d+01
        qse(11, 4) = 0.205390d+01
        qse(12, 4) = 0.205412d+01
        ! equilibrium geometry at anchor point 4
        qbe( 1, 4) = 0.215707d+01
        qbe( 2, 4) = 0.211094d+01
        qbe( 3, 4) = 0.201518d+01
        qbe( 4, 4) = 0.205661d+01
        qbe( 5, 4) = 0.208168d+01
        qbe( 6, 4) = 0.214490d+01
        qbe( 7, 4) = 0.207623d+01
        qbe( 8, 4) = 0.211213d+01
        qbe( 9, 4) = 0.209483d+01
        qbe(10, 4) = 0.214363d+01
        qbe(11, 4) = 0.206707d+01
        qbe(12, 4) = 0.207249d+01
        qbe(13, 4) = 0.207837d+01
        qbe(14, 4) = 0.209780d+01
        qbe(15, 4) = 0.210701d+01
        qbe(16, 4) = 0.205447d+01
        qbe(17, 4) = 0.207468d+01
        qbe(18, 4) = 0.215403d+01
    else if(ii .eq. 33) then
        ! equilibrium geometry at anchor point 1
        qse( 1, 1) = 0.272760d+01
        qse( 2, 1) = 0.272219d+01
        qse( 3, 1) = 0.243093d+01
        qse( 4, 1) = 0.259142d+01
        qse( 5, 1) = 0.205558d+01
        qse( 6, 1) = 0.266713d+01
        qse( 7, 1) = 0.205707d+01
        qse( 8, 1) = 0.267731d+01
        qse( 9, 1) = 0.205961d+01
        qse(10, 1) = 0.258566d+01
        qse(11, 1) = 0.205666d+01
        qse(12, 1) = 0.205567d+01
        ! equilibrium geometry at anchor point 1
        qbe( 1, 1) = 0.209913d+01
        qbe( 2, 1) = 0.214547d+01
        qbe( 3, 1) = 0.203859d+01
        qbe( 4, 1) = 0.207342d+01
        qbe( 5, 1) = 0.205792d+01
        qbe( 6, 1) = 0.215185d+01
        qbe( 7, 1) = 0.210210d+01
        qbe( 8, 1) = 0.209263d+01
        qbe( 9, 1) = 0.208845d+01
        qbe(10, 1) = 0.211486d+01
        qbe(11, 1) = 0.208558d+01
        qbe(12, 1) = 0.208274d+01
        qbe(13, 1) = 0.209348d+01
        qbe(14, 1) = 0.208777d+01
        qbe(15, 1) = 0.210193d+01
        qbe(16, 1) = 0.208337d+01
        qbe(17, 1) = 0.205127d+01
        qbe(18, 1) = 0.214854d+01
        ! equilibrium geometry at anchor point 2
        qse( 1, 2) = 0.274917d+01
        qse( 2, 2) = 0.274193d+01
        qse( 3, 2) = 0.237981d+01
        qse( 4, 2) = 0.259572d+01
        qse( 5, 2) = 0.205541d+01
        qse( 6, 2) = 0.266481d+01
        qse( 7, 2) = 0.205872d+01
        qse( 8, 2) = 0.266882d+01
        qse( 9, 2) = 0.205898d+01
        qse(10, 2) = 0.259216d+01
        qse(11, 2) = 0.205873d+01
        qse(12, 2) = 0.205800d+01
        ! equilibrium geometry at anchor point 2
        qbe( 1, 2) = 0.206394d+01
        qbe( 2, 2) = 0.213351d+01
        qbe( 3, 2) = 0.208574d+01
        qbe( 4, 2) = 0.208869d+01
        qbe( 5, 2) = 0.205086d+01
        qbe( 6, 2) = 0.214363d+01
        qbe( 7, 2) = 0.210572d+01
        qbe( 8, 2) = 0.209443d+01
        qbe( 9, 2) = 0.208304d+01
        qbe(10, 2) = 0.211155d+01
        qbe(11, 2) = 0.208615d+01
        qbe(12, 2) = 0.208549d+01
        qbe(13, 2) = 0.209340d+01
        qbe(14, 2) = 0.208795d+01
        qbe(15, 2) = 0.210183d+01
        qbe(16, 2) = 0.210307d+01
        qbe(17, 2) = 0.204251d+01
        qbe(18, 2) = 0.213761d+01
        ! equilibrium geometry at anchor point 3
        qse( 1, 3) = 0.274751d+01
        qse( 2, 3) = 0.274706d+01
        qse( 3, 3) = 0.235956d+01
        qse( 4, 3) = 0.259829d+01
        qse( 5, 3) = 0.205917d+01
        qse( 6, 3) = 0.266445d+01
        qse( 7, 3) = 0.205922d+01
        qse( 8, 3) = 0.266458d+01
        qse( 9, 3) = 0.205897d+01
        qse(10, 3) = 0.259790d+01
        qse(11, 3) = 0.205917d+01
        qse(12, 3) = 0.205900d+01
        ! equilibrium geometry at anchor point 3
        qbe( 1, 3) = 0.204756d+01
        qbe( 2, 3) = 0.212154d+01
        qbe( 3, 3) = 0.211409d+01
        qbe( 4, 3) = 0.210600d+01
        qbe( 5, 3) = 0.204389d+01
        qbe( 6, 3) = 0.213330d+01
        qbe( 7, 3) = 0.209805d+01
        qbe( 8, 3) = 0.210013d+01
        qbe( 9, 3) = 0.208501d+01
        qbe(10, 3) = 0.211048d+01
        qbe(11, 3) = 0.208627d+01
        qbe(12, 3) = 0.208644d+01
        qbe(13, 3) = 0.209637d+01
        qbe(14, 3) = 0.208569d+01
        qbe(15, 3) = 0.210112d+01
        qbe(16, 3) = 0.210790d+01
        qbe(17, 3) = 0.204012d+01
        qbe(18, 3) = 0.213517d+01
        ! equilibrium geometry at anchor point 4
        qse( 1, 4) = 0.274535d+01
        qse( 2, 4) = 0.274535d+01
        qse( 3, 4) = 0.236001d+01
        qse( 4, 4) = 0.259897d+01
        qse( 5, 4) = 0.205919d+01
        qse( 6, 4) = 0.266411d+01
        qse( 7, 4) = 0.205914d+01
        qse( 8, 4) = 0.266411d+01
        qse( 9, 4) = 0.205909d+01
        qse(10, 4) = 0.259897d+01
        qse(11, 4) = 0.205914d+01
        qse(12, 4) = 0.205919d+01
        ! equilibrium geometry at anchor point 4
        qbe( 1, 4) = 0.204757d+01
        qbe( 2, 4) = 0.211781d+01
        qbe( 3, 4) = 0.211781d+01
        qbe( 4, 4) = 0.210758d+01
        qbe( 5, 4) = 0.204104d+01
        qbe( 6, 4) = 0.213456d+01
        qbe( 7, 4) = 0.209631d+01
        qbe( 8, 4) = 0.210152d+01
        qbe( 9, 4) = 0.208536d+01
        qbe(10, 4) = 0.211101d+01
        qbe(11, 4) = 0.208609d+01
        qbe(12, 4) = 0.208609d+01
        qbe(13, 4) = 0.209631d+01
        qbe(14, 4) = 0.208536d+01
        qbe(15, 4) = 0.210152d+01
        qbe(16, 4) = 0.210758d+01
        qbe(17, 4) = 0.204104d+01
        qbe(18, 4) = 0.213456d+01
    endif

end subroutine getgeom
!> @brief Evaluate tent for the PHOH_APRP_DPEM potential.
subroutine tent( nap, rap, r, ir, t, dtdr )
    !***********************************************************************
    ! Subroutine to evaluate tent functions
    !***********************************************************************

    implicit none

    double precision :: dudr, dvdr

    integer          :: ir
    integer          :: nap
    double precision :: r
    double precision :: rap(*),t(*),dtdr(*)

    ! Local variables
    integer          :: i
    double precision :: u,v

    do i=1,nap
        t(i)    = 0.0d0
        dtdr(i) = 0.0d0
    enddo

    ! Check whether the point are end points
    if (r .lt. rap(1)) then
        ir      = 1
        t(1)    = 1.0d0
        dtdr(1) = 0.0d0
    else if (r .gt. rap(nap)) then
        ir        = nap+1
        t(nap)    = 1.0d0
        dtdr(nap) = 0.0d0
    endif

    ! Determine i for which r(i-1) <= r < r(i)
    do i=2,nap
        if ((rap(i-1) .le. r) .and. (r .lt. rap(i))) then
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
!> @brief Evaluate evu1r for the PHOH_APRP_DPEM potential.
subroutine evu1r( igrad, r, v, dvdr )
    !**********************************************************************
    ! Subroutine to evaluate Varshini potential and derivative for given r
    !**********************************************************************

    implicit none

    double precision, parameter :: de   = 0.17596055d0
    double precision, parameter :: re   = 1.93206106d0
    double precision, parameter :: beta = 0.63822703d0

    integer          :: igrad
    double precision :: r,v,dvdr
    double precision :: x,y,dvdy,dydx,dxdr

    x = r/re
    y = 1.0d0-exp(-beta*(x*x-1.0d0))/x
    v = de*y*y

    if (igrad .eq. 1) then
        dvdy = 2.0d0*de*y
        dydx = (2.0d0*beta*x*x+1.0d0)*exp(-beta*(x*x-1.0d0))/(x*x)
        dxdr = 1.0d0/re
        dvdr = dvdy*dydx*dxdr
    endif

end subroutine evu1r
!> @brief Evaluate evu1b for the PHOH_APRP_DPEM potential.
subroutine evu1b( igrad, r, theta, vbnd, dvdr, dvdth )
    !**********************************************************************
    ! Subroutine to evaluate potential energy from H-O-C bending
    ! and derivatives for given rOH and aHOC
    !**********************************************************************

    implicit none

    double precision :: cthta, cthta0, dct0dr, dk2dr, dthta, dtxdr, dtxdx, dvdr, dvdth, dxdr, gr, &
        & sthta, theta, tx, vbnd, x
    integer :: igrad

    integer, parameter :: ngk2 = 1

    double precision :: a2(ngk2)
    data a2/0.10449356d0/

    double precision :: alpha2(ngk2)
    data alpha2/0.27051185d0/

    double precision :: r2(ngk2)
    data r2/2.28431658d0/

    double precision :: cthta1,cthta2,a,r0
    data cthta1/-0.39964587d0/,cthta2/-0.24769126d0/
    data a/4.01097890d0/,r0/1.83449325d0/

    integer          :: i,j,k
    double precision :: r
    double precision :: k2

    cthta = cos(theta)
    sthta = sin(theta)

    x = a*(r-r0)
    !      tx=(exp(x)-exp(-x))/(exp(x)+exp(-x))
    tx     = (1.0d0-exp(-2.0d0*x))/(1.0d0+exp(-2.0d0*x))
    cthta0 = cthta1+0.5d0*(1+tx)*(cthta2-cthta1)

    dxdr   = a
    dtxdx  = 4.0d0/(2.0d0+exp(2.0d0*x)+exp(-2.0d0*x))
    dtxdr  = dtxdx*dxdr
    dct0dr = 0.5d0*(cthta2-cthta1)*dtxdr

    k2    = 0.0d0
    dk2dr = 0.0d0
    do i=1,ngk2
        gr = a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
        k2 = k2+gr
        if (igrad .eq. 1) then
            dk2dr = dk2dr-2.0d0*alpha2(i)*(r-r2(i))*gr
        endif
    enddo

    dthta = cthta-cthta0

    vbnd = k2*dthta*dthta

    if (igrad .eq. 1) then
        dvdr  = dk2dr*dthta*dthta - 2.0d0*k2*dthta*dct0dr
        dvdth = -2.0d0*k2*dthta*sthta
    endif

end subroutine evu1b

!> @brief Evaluate evu1to for the PHOH_APRP_DPEM potential.
subroutine evu1to( igrad, r, phi, v, dvdr, dvdph )
    !**********************************************************************
    ! Subroutine to evaluate torsion potentials and derivatives
    ! for given r and phi
    !**********************************************************************

    implicit none

    double precision :: dx1dr, dx2dr, dydph, x1, x2

    double precision, parameter :: a1     = 0.02371114d0
    double precision, parameter :: alpha1 = 0.96558083d0
    double precision, parameter :: r1     = 3.25505778d0
    double precision, parameter :: a2     = 0.00349299d0
    double precision, parameter :: alpha2 = 1.092945551d1
    double precision, parameter :: r2     = 2.10516007d0

    integer          :: igrad
    double precision :: phi,r,v,dvdph,dvdr
    double precision :: x,y,dvdy,dydx,dxdr

    x1 = a1*exp(-alpha1*(r-r1)*(r-r1))
    x2 = a2*exp(-alpha2*(r-r2)*(r-r2))
    y  = 1.0d0-cos(2.0d0*phi)
    v  = x1*y+x2*y*y

    if (igrad .eq. 1) then
        dx1dr = -2.0d0*alpha1*(r-r1)*x1
        dx2dr = -2.0d0*alpha2*(r-r2)*x2
        dydph = 2.0d0*sin(2.0d0*phi)
        dvdr  = dx1dr*y+dx2dr*y*y
        dvdph = x1*dydph+2.0d0*x2*y*dydph
    endif

end subroutine evu1to
!> @brief Evaluate evu2r for the PHOH_APRP_DPEM potential.
subroutine evu2r( igrad, r, v, dvdr )
    !**********************************************************************
    ! Subroutine to evaluate Morse potential and derivative for given r
    !**********************************************************************

    implicit none

    double precision, parameter :: alpha = 0.67269344d0
    double precision, parameter :: de    = 0.41866243d0
    double precision, parameter :: re    = 1.95347670d0
    double precision, parameter :: v0    = 0.17280158d0

    integer          :: igrad
    double precision :: r,v,dvdr
    double precision :: x,y,dvdy,dydx,dxdr

    x = r-re
    y = 1.0d0-exp(-alpha*x)
    v = de*y*y+v0

    if (igrad .eq. 1) then
        dvdy = 2.0d0*de*y
        dydx = alpha*exp(-alpha*x)
        dxdr = 1.0d0
        dvdr = dvdy*dydx*dxdr
    endif

end subroutine evu2r
!> @brief Evaluate evu2b for the PHOH_APRP_DPEM potential.
subroutine evu2b( igrad, r, theta, vbnd, dvdr, dvdth )
    !**********************************************************************
    ! Subroutine to evaluate potential energy from H-O-C bending
    ! and derivatives for given rOH and aHOC
    !**********************************************************************

    implicit none

    double precision :: cthta, cthta0, dct0dr, dk2dr, dthta, dtxdr, dtxdx, dvdr, dvdth, dxdr, gr, &
        & sthta, theta, tx, vbnd, x
    integer :: igrad

    integer, parameter :: ngk2 = 1

    double precision :: a2(ngk2)
    data a2/0.10090668d0/

    double precision :: alpha2(ngk2)
    data alpha2/0.21154184d0/

    double precision :: r2(ngk2)
    data r2/1.80972781d0/

    double precision :: cthta1,cthta2,a,r0
    data cthta1/-0.31835353d0/,cthta2/-0.22230161d0/
    data a/4.73544047d0/,r0/2.01933829d0/

    integer          :: i,j,k
    double precision :: r
    double precision :: k2

    cthta = cos(theta)
    sthta = sin(theta)

    x = a*(r-r0)
    !      tx=(exp(x)-exp(-x))/(exp(x)+exp(-x))
    tx     = (1.0d0-exp(-2.0d0*x))/(1.0d0+exp(-2.0d0*x))
    cthta0 = cthta1+0.5d0*(1+tx)*(cthta2-cthta1)

    dxdr   = a
    dtxdx  = 4.0d0/(2.0d0+exp(2.0d0*x)+exp(-2.0d0*x))
    dtxdr  = dtxdx*dxdr
    dct0dr = 0.5d0*(cthta2-cthta1)*dtxdr

    k2    = 0.0d0
    dk2dr = 0.0d0
    do i=1,ngk2
        gr = a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
        k2 = k2+gr
        if (igrad .eq. 1) then
            dk2dr = dk2dr-2.0d0*alpha2(i)*(r-r2(i))*gr
        endif
    enddo

    dthta = cthta-cthta0

    vbnd = k2*dthta*dthta

    if (igrad .eq. 1) then
        dvdr  = dk2dr*dthta*dthta - 2.0d0*k2*dthta*dct0dr
        dvdth = -2.0d0*k2*dthta*sthta
    endif

end subroutine evu2b

!> @brief Evaluate evu2to for the PHOH_APRP_DPEM potential.
subroutine evu2to( igrad, r, phi, v, dvdr, dvdph )
    !**********************************************************************
    ! Subroutine to evaluate torsion potentials and derivatives
    ! for given r and phi
    !**********************************************************************

    implicit none

    double precision :: dydph

    double precision, parameter :: a1     = 0.01925782d0
    double precision, parameter :: alpha1 = 5.00873578d0
    double precision, parameter :: r1     = 2.22678700d0

    integer          :: igrad
    double precision :: phi,r,v,dvdph,dvdr
    double precision :: x,y,dvdy,dydx,dxdr

    x = a1*exp(-alpha1*(r-r1)*(r-r1))
    y = 1.0d0-cos(2.0d0*phi)
    v = x*y

    if (igrad .eq. 1) then
        dxdr  = -2.0d0*alpha1*(r-r1)*x
        dydph = 2.0d0*sin(2.0d0*phi)
        dvdr  = y*dxdr
        dvdph = x*dydph
    endif

end subroutine evu2to
!> @brief Evaluate evu3r for the PHOH_APRP_DPEM potential.
subroutine evu3r( igrad, r, v, dvdr )
    !**********************************************************************
    ! Subroutine to evaluate Morse potential and derivative for given r
    !**********************************************************************

    implicit none

    double precision :: dx1dr, dx2dr, dx3dr, dy1dx1, dy2dx2, dy3dx3, x1, x2, x3, y1, y2, y3

    double precision, parameter :: v0     = 0.15976218d0
    double precision, parameter :: a1     = 0.06496177d0
    double precision, parameter :: alpha1 = 1.81347726d0
    double precision, parameter :: r1     = 2.63168269d0
    double precision, parameter :: a2     = -0.00077880d0
    double precision, parameter :: alpha2 = 3.24769899d0
    double precision, parameter :: r2     = 3.71412291d0
    double precision, parameter :: a3     = 9.3151142407d2
    double precision, parameter :: alpha3 = 4.73016534d0
    double precision, parameter :: r3     = -0.04117486d0

    integer          :: igrad
    double precision :: r,v,dvdr
    double precision :: x,y,dvdy,dydx,dxdr

    x1 = r-r1
    y1 = a1*exp(-alpha1*x1)
    x2 = r-r2
    y2 = a2*exp(-alpha2*x2)
    x3 = r-r3
    y3 = a3*exp(-alpha3*x3)
    v  = v0+y1+y2+y3

    if (igrad .eq. 1) then
        dx1dr  = 1.0d0
        dx2dr  = 1.0d0
        dx3dr  = 1.0d0
        dy1dx1 = -alpha1*y1
        dy2dx2 = -alpha2*y2
        dy3dx3 = -alpha3*y3
        dvdr   = dy1dx1*dx1dr+dy2dx2*dx2dr+dy3dx3*dx3dr
    endif

end subroutine evu3r
!> @brief Evaluate evu3b for the PHOH_APRP_DPEM potential.
subroutine evu3b( igrad, r, theta, vbnd, dvdr, dvdth )
    !**********************************************************************
    ! Subroutine to evaluate potential energy from H-O-C bending
    ! and derivatives for given rOH and aHOC
    !**********************************************************************

    implicit none

    double precision :: cthta, cthta0, dct0dr, dk2dr, dk3dr, dthta, dtxdr, dtxdx, dvdr, dvdth, dxdr, &
        & gr, sthta, theta, tx, vbnd, x
    integer :: igrad

    integer, parameter :: ngk2 = 1
    integer, parameter :: ngk3 = 1

    double precision :: a2(ngk2),a3(ngk3)
    data a2/3.00968695d0/
    data a3/0.06538039d0/

    double precision :: alpha2(ngk2),alpha3(ngk3)
    data alpha2/0.04954325d0/
    data alpha3/14.94793100d0/

    double precision :: r2(ngk2),r3(ngk3)
    data r2/-6.33258333d0/
    data r3/2.09315608d0/

    double precision :: cthta1,cthta2,a,r0
    data cthta1/-0.30362526d0/,cthta2/-0.65443171d0/
    data a/0.90134535d0/,r0/2.95896455d0/

    integer          :: i,j,k
    double precision :: r
    double precision :: k2,k3

    cthta = cos(theta)
    sthta = sin(theta)

    x = a*(r-r0)
    !      tx=(exp(x)-exp(-x))/(exp(x)+exp(-x))
    tx     = (1.0d0-exp(-2.0d0*x))/(1.0d0+exp(-2.0d0*x))
    cthta0 = cthta1+0.5d0*(1+tx)*(cthta2-cthta1)

    dxdr   = a
    dtxdx  = 4.0d0/(2.0d0+exp(2.0d0*x)+exp(-2.0d0*x))
    dtxdr  = dtxdx*dxdr
    dct0dr = 0.5d0*(cthta2-cthta1)*dtxdr

    k2    = 0.0d0
    dk2dr = 0.0d0
    do i=1,ngk2
        gr = a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
        k2 = k2+gr
        if (igrad .eq. 1) then
            dk2dr = dk2dr-2.0d0*alpha2(i)*(r-r2(i))*gr
        endif
    enddo

    k3    = 0.0d0
    dk3dr = 0.0d0
    do i=1,ngk3
        gr = a3(i)*exp(-alpha3(i)*(r-r3(i))*(r-r3(i)))
        k3 = k3+gr
        if (igrad .eq. 1) then
            dk3dr = dk3dr-2.0d0*alpha3(i)*(r-r3(i))*gr
        endif
    enddo

    dthta = cthta-cthta0
    vbnd  = k2*dthta*dthta + k3*dthta*dthta*dthta
    if (igrad .eq. 1) then
        dvdr  = dk2dr*dthta*dthta + dk3dr*dthta*dthta*dthta - 2.0d0*k2*dthta*dct0dr - 3.0d0*k3*dthta*dthta*dct0dr
        dvdth = -2.0d0*k2*dthta*sthta - 3.0d0*k3*dthta*dthta*sthta
    endif

end subroutine evu3b

!> @brief Evaluate evu3to for the PHOH_APRP_DPEM potential.
subroutine evu3to( igrad, r, phi, v, dvdr, dvdph )
    !**********************************************************************
    ! Subroutine to evaluate torsion potentials and derivatives
    ! for given r and phi
    !**********************************************************************

    implicit none

    double precision :: dx1dr, dx2dr, dydph, x1, x2

    double precision, parameter :: a1     = -0.02298368d0
    double precision, parameter :: alpha1 = 1.06740430d0
    double precision, parameter :: r1     = 3.28664483d0
    double precision, parameter :: a2     = 0.01138152d0
    double precision, parameter :: alpha2 = 0.65339172d0
    double precision, parameter :: r2     = 2.58341105d0

    integer          :: igrad
    double precision :: phi,r,v,dvdph,dvdr
    double precision :: x,y,dvdy,dydx,dxdr

    x1 = a1*exp(-alpha1*(r-r1)*(r-r1))
    x2 = a2*exp(-alpha2*(r-r2)*(r-r2))
    y  = 1.0d0-cos(2.0d0*phi)
    v  = x1*y+x2*y*y

    if (igrad .eq. 1) then
        dx1dr = -2.0d0*alpha1*(r-r1)*x1
        dx2dr = -2.0d0*alpha2*(r-r2)*x2
        dydph = 2.0d0*sin(2.0d0*phi)
        dvdr  = dx1dr*y+dx2dr*y*y
        dvdph = x1*dydph+2.0d0*x2*y*dydph
    endif

end subroutine evu3to
!> @brief Evaluate evubm for the PHOH_APRP_DPEM potential.
subroutine evubm( igrad, natm, x, y, z, ubm, dubmdx )
    !***********************************************************************
    ! Subroutine to evulate the Born-Mayer repulsions between non-bonded
    ! C atoms
    !***********************************************************************

    implicit none

    double precision, parameter :: b     = 6.693126d1
    double precision, parameter :: alpha = 1.894454d0

    integer          :: i,j
    integer          :: igrad,natm
    double precision :: ubm
    double precision :: x(natm),y(natm),z(natm)
    double precision :: dubmdx(3*natm)
    double precision :: bndlen
    external bndlen
    ! Local variables
    integer :: nbl,nba,nto,noop
    data nbl/3/,nba/0/,nto/0/,noop/0/
    integer :: ibl(2,3),iba(3,0),ito(4,0),ioop(4,0)
    data ibl /1,4, 2,5, 3,6/
    double precision :: r14,r25,r36
    double precision :: dudr14,dudr25,dudr36
    double precision :: bm(3,3*natm)

    r14 = bndlen(1,4,x,y,z)
    r25 = bndlen(2,5,x,y,z)
    r36 = bndlen(3,6,x,y,z)

    ubm = b*(exp(-alpha*r14)+exp(-alpha*r25)+exp(-alpha*r36))

    !
    !      write(*,*) "ubm =",ubm
    !
    if (igrad .eq. 1) then
        dudr14 = -alpha*b*exp(-alpha*r14)
        dudr25 = -alpha*b*exp(-alpha*r25)
        dudr36 = -alpha*b*exp(-alpha*r36)

        call bmat3(nbl,ibl,nba,iba,nto,ito,noop,ioop,natm,x,y,z,bm,3)
        do i=1,3*natm
            dubmdx(i) = 0.d0
        enddo

        do i=1,3*natm
            dubmdx(i) = dubmdx(i)+dudr14*bm(1,i)+dudr25*bm(2,i) +dudr36*bm(3,i)
        enddo

        !
        !        do i=1,3*natm
        !          write(*,*) dubmdX(i)
        !        enddo
        !
    endif

end subroutine evubm
!      subroutine evuij(nap,natm,igrad,ir,t,dtdr,X,Y,Z,u13,u23,
!     $                 du13dX,du23dX)
! Add u12
!> @brief Evaluate evuij for the PHOH_APRP_DPEM potential.
subroutine evuij( nap, natm, igrad, ir, t, dtdr, x, y, z, u12, u13, u23, &
        & du12dx, du13dx, du23dx )
    !***********************************************************************
    ! Subroutine to evaluate the diabatic couplings and gradients
    ! The bohr-hartree-radian unit is used throughout the program
    ! qtc:  redundant internal coordinates
    ! stc:  non-redundant internal coordinates
    ! Ct:   transform matrix between RICs and NRICs
    !       stc=Ct*qtc
    ! bmqt: B-matrix of RICs
    ! bmst: B-matrix of NRICs
    !***********************************************************************

    implicit none

    double precision :: bndlen, d13dph, d23dph, dihedr, du12dr, phi, r, u12t, u13s, u13t, u23s, u23t
    integer          :: i, j, nap

    integer, parameter :: maxdim = 100

    integer, parameter :: ldbmt  = maxdim
    integer, parameter :: ldbms  = maxdim
    integer, parameter :: ldbnrt = maxdim
    integer, parameter :: ldct   = maxdim

    integer, parameter :: nstc = 9

    double precision, parameter :: pi = 3.1415926535898d0

    ! Global variables
    integer :: igrad,ir
    integer :: natm
    !      double precision u13,u23
    double precision :: u12,u13,u23
    double precision :: t(*),dtdr(*)
    double precision :: x(*),y(*),z(*)
    !      double precision du13dX(3*natm),du23dX(3*natm)
    double precision :: du12dx(3*natm),du13dx(3*natm),du23dx(3*natm)

    ! Local variables
    integer :: nqtc
    integer :: nblt,nbat,ntot,noopt
    !     data(nblt=0,nbat=0,ntot=6,noopt=6)
    data nblt/0/,nbat/0/,ntot/6/,noopt/6/
    data nblt/0/,nbat/0/,ntot/6/,noopt/6/
    integer :: iblt(2,0),ibat(3,0),itot(4,6),ioopt(4,6)
    data itot  /2,1,6,5, &
        & 1,6,5,4, &
        & 6,5,4,3, &
        & 5,4,3,2, &
        & 4,3,2,1, &
        & 3,2,1,6/
    data ioopt /7,1,2,6, &
        & 12,6,1,5, &
        & 11,5,6,4, &
        & 10,4,5,3, &
        & 9,3,4,2, &
        & 8,2,3,1/

    integer :: nbls,nbas,ntos,noops
    !      data(nbls=1,nbas=0,ntos=1,noops=0)
    data nbls/1/,nbas/0/,ntos/1/,noops/0/

    integer :: ibls(2,1),ibas(3,0),itos(4,1),ioops(4,0)
    data ibls /7,13/
    data itos /2,1,7,13/

    double precision, parameter :: half   = 0.5d0
    double precision, parameter :: oossix = 0.408248290463863d0
    double precision, parameter :: oostlv = 0.288675134594813d0

    double precision :: qtc(maxdim),qsc(maxdim)
    double precision :: stc(maxdim)
    ! KRY 20140718
    double precision :: d12tdr,d12sdr
    double precision :: d13tdr,d13sdr,du13dr,du13dph
    double precision :: d23tdr,d23sdr,du23dr,du23dph
    !      double precision du13t(maxdim),du23t(maxdim)
    double precision :: du12t(maxdim),du13t(maxdim),du23t(maxdim)
    double precision :: bmt(ldbmt,3*natm),bms(ldbms,3*natm)
    double precision :: bnrt(ldbnrt,3*natm)
    double precision :: ct(ldct,maxdim)

    ! Initialize Ct matrix
    do i=1,9
        do j=1,12
            ct(i,j) = 0.0d0
        enddo
    enddo
    ! Initialize Qtc vectors
    do i=1,nblt+nbat+ntot+noopt
        qtc(i) = 0.0d0
    enddo
    ! Initialize Stc vectors
    do i=1,nstc
        stc(i) = 0.0d0
    enddo

    ct(1,1)  = oossix
    ct(1,2)  = -oossix
    ct(1,3)  = oossix
    ct(1,4)  = -oossix
    ct(1,5)  = oossix
    ct(1,6)  = -oossix
    ct(2,1)  = -oostlv
    ct(2,2)  = 2.0d0*oostlv
    ct(2,3)  = -oostlv
    ct(2,4)  = -oostlv
    ct(2,5)  = 2.0d0*oostlv
    ct(2,6)  = -oostlv
    ct(3,1)  = -half
    ct(3,3)  = half
    ct(3,4)  = -half
    ct(3,6)  = half
    ct(4,7)  = 1.0d0
    ct(5,8)  = 1.0d0
    ct(6,9)  = 1.0d0
    ct(7,10) = 1.0d0
    ct(8,11) = 1.0d0
    ct(9,12) = 1.0d0

    ! Evaluate internal coordinates
    !Add by J. Zheng
    !      qtc(:) = 0d0
    call evintcoord(nblt,iblt,nbat,ibat,ntot,itot,noopt,ioopt, &
        & nqtc,qtc,x,y,z)
    ! Test
    !      do i=1,nqtc
    !        write(*,'(A,I2,A,F8.4)') 'qtc(',i,') = ',qtc(i)*180.d0/pi
    !      enddo
    ! End

    !Add by J. Zheng
    !      stc(:) = 0d0
    call matvect(stc,nstc,ct,ldct,qtc,nqtc)

    r   = bndlen(7,13,x,y,z)
    phi = dihedr(2,1,7,13,x,y,z)

    ! Evaluate Wilson's B-matrix if igrad = 1
    if (igrad .eq. 1) then
        !Added by J. Zheng
        !        bmt(:,:) = 0d0
        call bmat3(nblt,iblt,nbat,ibat,ntot,itot,noopt,ioopt,natm, &
            & x,y,z,bmt,ldbmt)
        ! Evalate the Wilson's B-matrix in nonredundant internal coordinates
        call matmult(bnrt,ldbnrt,ct,ldct,nstc,bmt,ldbmt,nqtc,3*natm)
        !Added by J. Zheng
        !        bms(:,:) = 0d0
        call bmat3(nbls,ibls,nbas,ibas,ntos,itos,noops,ioops,natm, &
            & x,y,z,bms,ldbms)
    endif

    ! Evaluate u12 KRY 20140718
    call evuijt(nap,igrad,ir,t,dtdr,nstc,stc,u12t,du12t,d12tdr,12)
    ! Evaluate u13 and u23 from tertiary out-of-plane coordinates
    call evuijt(nap,igrad,ir,t,dtdr,nstc,stc,u13t,du13t,d13tdr,13)
    call evuijt(nap,igrad,ir,t,dtdr,nstc,stc,u23t,du23t,d23tdr,23)

    ! Evaluate u13 and u23 from secondary out-of-plane coordinates
    call evu13to(r,phi,u13s,d13dph,d13sdr)
    call evu23to(r,phi,u23s,d23dph,d23sdr)

    u12 = u12t                  ! KRY 20140718
    u13 = u13t+u13s
    u23 = u23t+u23s

    ! Calculated analytical Cartesian gradient if igrad = 1
    if (igrad .eq. 1) then
        do i=1,3*natm
            du12dx(i) = 0.0d0
            du13dx(i) = 0.0d0
            du23dx(i) = 0.0d0
        enddo
        du12dr = d12tdr
        du13dr = d13tdr+d13sdr
        du23dr = d23tdr+d23sdr
        ! Contribution of duijdX from duijdr
        do i=1,3*natm
            du12dx(i) = du12dx(i)+du12dr*bms(1,i)
            du13dx(i) = du13dx(i)+du13dr*bms(1,i)
            du23dx(i) = du23dx(i)+du23dr*bms(1,i)
        enddo
        ! Contribution of duijdX from dijdph
        do i=1,3*natm
            du13dx(i) = du13dx(i)+d13dph*bms(2,i)
            du23dx(i) = du23dx(i)+d23dph*bms(2,i)
        enddo
        ! Contribution of duijdX from duijt
        do i=1,3*natm
            do j=1,nstc
                du12dx(i) = du12dx(i)+du12t(j)*bnrt(j,i)
                du13dx(i) = du13dx(i)+du13t(j)*bnrt(j,i)
                du23dx(i) = du23dx(i)+du23t(j)*bnrt(j,i)
            enddo
        enddo
    endif

end subroutine evuij
!> @brief Evaluate evuijt for the PHOH_APRP_DPEM potential.
subroutine evuijt( nap, igrad, ir, t, dtdr, ntc, qtc, uij, &
        & duij, duijdr, ij )

    implicit none

    integer :: nap

    integer          :: ij
    integer          :: igrad,ir
    integer          :: ntc
    double precision :: t(*),dtdr(*)
    double precision :: qtc(*)
    double precision :: uij
    double precision :: duijdr
    double precision :: duij(*)

    ! Local variables
    integer          :: i,j
    double precision :: u(nap)
    double precision :: a(nap,ntc)
    double precision :: uij0(nap)
    double precision :: uij1,uij2
    double precision :: duij1(ntc),duij2(ntc)

    ! Initialize uij,uij1,uij2
    uij  = 0.0d0
    uij1 = 0.0d0
    uij2 = 0.0d0
    ! Initialize duij,duij1,duij2
    do i=1,ntc
        duij(i)  = 0.0d0
        duij1(i) = 0.0d0
        duij2(i) = 0.0d0
    enddo
    ! Initialize duijdr
    duijdr = 0.0d0

    ! Initialize numerical gradients of diabatic couplings
    do i=1,nap
        uij0(i) = 0.0d0
        do j=1,9
            a(i,j) = 0.0d0
        enddo
    enddo
    ! Assign u120 for anchor point 1 and 2
    if (ij .eq. 12) then
        uij0(1) = -0.000616170d0
        uij0(2) = -0.001064286d0
    endif

    ! Assign parameters based ij
    if (ij .eq. 12) then
        a(3,1) = 0.0d0
        a(3,2) = 0.0d0
        a(3,3) = 0.0d0
        a(3,4) = 0.0d0
        a(3,5) = 0.0d0
        a(3,6) = 0.0d0
        a(3,7) = 0.0d0
        a(3,8) = 0.0d0
        a(3,9) = 0.0d0
    else if (ij .eq. 13) then
        a(3,1) = -0.00070d0
        a(3,2) = -0.00966d0
        a(3,3) = -0.00053d0
        a(3,4) = 0.00193d0
        a(3,5) = 0.00056d0
        a(3,6) = -0.00362d0
        a(3,7) = 0.00030d0
        a(3,8) = 0.00377d0
        a(3,9) = 0.00101d0
    else if (ij .eq. 23) then
        a(2,1) = 0.00452d0
        a(2,2) = -0.01212d0
        a(2,3) = -0.00088d0
        a(2,4) = 0.00532d0
        a(2,5) = 0.00757d0
        a(2,6) = -0.00348d0
        a(2,7) = 0.00023d0
        a(2,8) = -0.00209d0
        a(2,9) = -0.01501d0
    else
        write(*,*) "Only parameters for U12, U13 and U23 are allowed!"
        stop
    endif

    ! Check for consistence
    if (ntc .ne. 9) then
        write(*,*) "Inconsistence in the number of out-of-plane tertiary&
            & coordinates"
        stop
    endif
    if (nap .ne. 4) then
        write(*,*) "Inconsistence in the number of anchor structures"
        stop
    endif

    ! Now calculate coupling contributions from tertiary coordinates
    if (ir .eq. 1) then
        uij = uij0(1)
        do i=1,ntc
            uij = uij+a(1,i)*qtc(i)
        enddo
    else if (ir .eq. nap+1) then
        uij = uij0(nap)
        do i=1,ntc
            uij = uij+a(nap,i)*qtc(i)
        enddo
    else if ((ir .gt. 1) .and. (ir .lt. nap+1)) then
        uij1 = uij0(ir-1)
        uij2 = uij0(ir)
        do i=1,ntc
            uij1 = uij1+a(ir-1,i)*qtc(i)
            uij2 = uij2+a(ir,i)*qtc(i)
        enddo
        uij = uij1*t(ir-1)+uij2*t(ir)
    endif

    ! Now calculate the derivatives of couplings if igrad =1
    if (igrad .eq. 1) then
        if (ir .eq. 1) then
            do i=1,ntc
                duij(i) = a(1,i)
            enddo
            duijdr = uij*dtdr(1)
        else if (ir .eq. nap+1) then
            do i=1,ntc
                duij(i) = a(nap,i)
            enddo
            duijdr = uij*dtdr(nap)
        else if ((ir .gt. 1) .and. (ir .lt. nap+1)) then
            do i=1,ntc
                duij1(i) = a(ir-1,i)
                duij2(i) = a(ir,i)
                duij(i)  = duij1(i)*t(ir-1)+duij2(i)*t(ir)
            enddo
            duijdr = uij1*dtdr(ir-1)+uij2*dtdr(ir)
        endif
    endif

end subroutine evuijt
!> @brief Evaluate evu13to for the PHOH_APRP_DPEM potential.
subroutine evu13to( r, phi, u13to, d13dph, d13dr )
    !**********************************************************************
    ! Little subroutine to evaluate U13 from C-C-O-H torsion
    ! U13(r,phi)=f1(r)*sin(phi)+f3(r)*sin(phi)**3+f5(r)*sin(phi)**5
    ! where fi(r) are expanded with a linear combination of Gaussians
    !**********************************************************************

    implicit none

    double precision :: df1rdr, df3rdr, df5rdr

    integer, parameter :: ngf = 3

    double precision :: a1(ngf),a3(ngf),a5(ngf)
    data a1/0.068885d0,0.034221d0,0.001544d0/
    data a3/-0.082676d0,-0.021770d0,0.0d0/
    data a5/0.053915d0,-0.008539d0,.0d0/

    double precision :: alpha1(ngf),alpha3(ngf),alpha5(ngf)
    data alpha1/0.2956970d1,0.4695508d1,0.1167016d1/
    data alpha3/0.5029657d1,0.8744002d1,0.0d0/
    data alpha5/0.7263172d1,0.27942459d2,0.0d0/

    double precision :: r1(ngf),r3(ngf),r5(ngf)
    data r1/0.3162075d1,0.2255351d1,0.4480427d1/
    data r3/0.3286539d1,0.2060405d1,0.0d0/
    data r5/0.3275937d1,0.2498052d1,0.0d0/

    double precision :: r,phi,u13to,d13dph,d13dr

    ! Local variables
    integer          :: i
    double precision :: sphi,cphi
    double precision :: sphi1,sphi2,sphi3,sphi4,sphi5
    double precision :: f1r,f3r,f5r
    double precision :: df1dr,df3dr,df5dr
    double precision :: g1r,g3r,g5r

    u13to  = 0.0d0
    d13dph = 0.0d0
    d13dr  = 0.0d0

    f1r    = 0.0d0
    f3r    = 0.0d0
    f5r    = 0.0d0
    df1rdr = 0.0d0
    df3rdr = 0.0d0
    df5rdr = 0.0d0
    do i=1,ngf
        g1r    = a1(i)*exp(-alpha1(i)*(r-r1(i))*(r-r1(i)))
        f1r    = f1r+g1r
        df1rdr = df1rdr-2.0d0*alpha1(i)*(r-r1(i))*g1r
        g3r    = a3(i)*exp(-alpha3(i)*(r-r3(i))*(r-r3(i)))
        f3r    = f3r+g3r
        df3rdr = df3rdr-2.0d0*alpha3(i)*(r-r3(i))*g3r
        g5r    = a5(i)*exp(-alpha5(i)*(r-r5(i))*(r-r5(i)))
        f5r    = f5r+g5r
        df5rdr = df5rdr-2.0d0*alpha5(i)*(r-r5(i))*g5r
    enddo

    sphi  = dsin(phi)
    cphi  = dcos(phi)
    sphi1 = sphi
    sphi2 = sphi*sphi1
    sphi3 = sphi*sphi2
    sphi4 = sphi*sphi3
    sphi5 = sphi*sphi4

    u13to  = f1r*sphi1+f3r*sphi3+f5r*sphi5
    d13dph = f1r*cphi+3.0d0*f3r*sphi2*cphi+5.0d0*f5r*sphi4*cphi
    d13dr  = df1rdr*sphi1+df3rdr*sphi3+df5rdr*sphi5

end subroutine evu13to
!> @brief Evaluate evu23to for the PHOH_APRP_DPEM potential.
subroutine evu23to( r, phi, u23to, d23dph, d23dr )
    !**********************************************************************
    ! Little subroutine to evaluate U23 from C-C-O-H torsion
    ! U23(r,phi)=f2(r)*sin(2*phi)+f4(r)*sin(4*phi)+f6(r)*sin(6*phi)
    ! where fi(r) are expanded with a linear combination of Gaussians
    !**********************************************************************

    implicit none

    double precision :: c2phi, c4phi, c6phi, df2rdr, df4rdr, df6rdr, s2phi, s4phi, s6phi

    integer, parameter :: ngf = 1

    double precision, parameter :: scale = 1.0d0

    double precision :: a2(ngf),a4(ngf),a6(ngf)
    data a2/0.016924d0/
    data a4/0.007145d0/
    data a6/0.003113d0/

    double precision :: alpha2(ngf),alpha4(ngf),alpha6(ngf)
    data alpha2/0.10841076d2/
    data alpha4/0.12574937d2/
    data alpha6/0.25329325d2/

    double precision :: r2(ngf),r4(ngf),r6(ngf)
    data r2/0.2151852d1/
    data r4/0.2197436d1/
    data r6/0.2285492d1/

    double precision :: r,phi,u23to,d23dph,d23dr

    ! Local variables
    integer          :: i
    double precision :: sphi,cphi
    double precision :: sphi1,sphi2,sphi3,sphi4,sphi5
    double precision :: f2r,f4r,f6r
    double precision :: df2dr,df4dr,df6dr
    double precision :: g2r,g4r,g6r

    u23to  = 0.0d0
    d23dph = 0.0d0
    d23dr  = 0.0d0

    f2r    = 0.0d0
    f4r    = 0.0d0
    f6r    = 0.0d0
    df2rdr = 0.0d0
    df4rdr = 0.0d0
    df6rdr = 0.0d0
    do i=1,ngf
        g2r    = a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
        f2r    = f2r+g2r
        df2rdr = df2rdr-2.0d0*alpha2(i)*(r-r2(i))*g2r
        g4r    = a4(i)*exp(-alpha4(i)*(r-r4(i))*(r-r4(i)))
        f4r    = f4r+g4r
        df4rdr = df4rdr-2.0d0*alpha4(i)*(r-r4(i))*g4r
        g6r    = a6(i)*exp(-alpha6(i)*(r-r6(i))*(r-r6(i)))
        f6r    = f6r+g6r
        df6rdr = df6rdr-2.0d0*alpha6(i)*(r-r6(i))*g6r
    enddo

    s2phi = dsin(2.0d0*phi)
    c2phi = dcos(2.0d0*phi)
    s4phi = dsin(4.0d0*phi)
    c4phi = dcos(4.0d0*phi)
    s6phi = dsin(6.0d0*phi)
    c6phi = dcos(6.0d0*phi)

    !      u23to=f2r*s2phi+f4r*s4phi+f6r*s6phi
    !      d23dph=2.0d0*f2r*c2phi+4.0d0*f4r*c4phi+6.0d0*f6r*c6phi
    !      d23dr=df2rdr*s2phi+df4rdr*s4phi+df6rdr*s6phi
    u23to  = scale*(f2r*s2phi+f4r*s4phi+f6r*s6phi)
    d23dph = scale*(2.0d0*f2r*c2phi+4.0d0*f4r*c4phi+6.0d0*f6r*c6phi)
    d23dr  = scale*(df2rdr*s2phi+df4rdr*s4phi+df6rdr*s6phi)

end subroutine evu23to
!> @brief Evaluate evintcoord for the PHOH_APRP_DPEM potential.
subroutine evintcoord( nbl, ibl, nba, iba, nto, ito, noop, ioop, &
        & nint, ric, x, y, z )
    !***********************************************************************
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
    !***********************************************************************

    implicit none

    double precision :: bndang, bndlen, dihedr, oopba1

    integer, parameter          :: maxint = 1000
    double precision, parameter :: btoa   = 0.5291772192d0
    double precision, parameter :: pi     = 3.1415926535898d0

    integer :: i
    integer :: nbl,nba,noop,nto,nint
    integer :: ibl(2,*),iba(3,*),ito(4,*),ioop(4,*)

    double precision :: x(*),y(*),z(*)
    double precision :: rab(maxint),theta(maxint),tors(maxint), &
        & oopb(maxint)
    double precision :: ric(*)

    ! Initialize RIC
    do i=1,nbl+nba+nto+noop
        ric(i) = 0.0d0
    enddo
    ! Initialize rab,theta,tors, and oopb
    do i=1,maxint
        rab(i)   = 0.0d0
        theta(i) = 0.0d0
        tors(i)  = 0.0d0
        oopb(i)  = 0.0d0
    enddo

    ! Evaluate bond lengths
    do i=1,nbl
        rab(i) = bndlen(ibl(1,i),ibl(2,i),x,y,z)
    enddo

    ! Evaluate bend angles
    do i=1,nba
        theta(i) = bndang(iba(1,i),iba(2,i),iba(3,i),x,y,z)
    enddo

    ! Evaluate torsions
    do i=1,nto
        tors(i) = dihedr(ito(1,i),ito(2,i),ito(3,i),ito(4,i),x,y,z)
    enddo

    ! Evaluate out-of-plane bending angles
    do i=1,noop
        oopb(i) = oopba1(ioop(1,i),ioop(2,i),ioop(3,i),ioop(4,i),x,y,z)
    enddo

    nint = 0
    do i=1,nbl
        nint      = nint+1
        ric(nint) = rab(i)
    enddo
    do i=1,nba
        nint      = nint+1
        ric(nint) = theta(i)
    enddo
    do i=1,nto
        nint      = nint+1
        ric(nint) = tors(i)
    enddo
    do i=1,noop
        nint      = nint+1
        ric(nint) = oopb(i)
    enddo

end subroutine evintcoord
!> @brief Evaluate matvect for the PHOH_APRP_DPEM potential.
subroutine matvect( y, m, a, lda, x, n )
    !**********************************************************************
    ! Subroutien to calculate matrix A multiply vector X
    !**********************************************************************

    implicit none

    integer          :: i, j
    integer          :: lda,m,n
    double precision :: x(*),y(*)
    double precision :: a(lda,*)

    ! Initialize Y
    do i=1,m
        y(i) = 0.0d0
    enddo

    do i=1,m
        do j=1,n
            y(i) = y(i)+a(i,j)*x(j)
        enddo
    enddo

end subroutine matvect
!> @brief Evaluate bmat3 for the PHOH_APRP_DPEM potential.
subroutine bmat3( nbl, ibl, nba, iba, nto, ito, noop, ioop, natm, &
        & x, y, z, bm, ldbm )
    !***********************************************************************
    ! Subroutine to evaluate the Wilson B-matrix
    !***********************************************************************

    implicit none

    double precision :: ai, aijk, ajkl, ak, al, bndang, bndlen, ci, cijk, cjkl, ck, cl, ctheta, &
        & oopba1, rij, rji, rjk, rjl, rkl, si, sijk, sjkl, sk, sl, theta, ttheta

    ! Global variables
    integer          :: ldbm
    integer          :: nbl,nba,noop,nto,natm
    integer          :: ibl(2,nbl),iba(3,nba),ito(4,nto),ioop(4,noop)
    double precision :: x(*),y(*),z(*)
    double precision :: bm(ldbm,*)

    ! Local variables
    integer          :: i,j,k,l,ij,ijk,ijkl,irow
    double precision :: eij(3),ejk(3),ekl(3),eji(3),ejl(3),u(3),v(3)
    double precision :: w(3)

    ! Initialize B matrix
    do i=1,ldbm
        do j=1,3*natm
            bm(i,j) = 0.0d0
        enddo
    enddo

    irow = 0
    ! Evaluate the B matrix contributed from bond stretchings
    do ij=1,nbl
        irow = irow+1
        i    = ibl(1,ij)
        j    = ibl(2,ij)

        rij = bndlen(i,j,x,y,z)

        eij(1) = (x(j)-x(i))/rij
        eij(2) = (y(j)-y(i))/rij
        eij(3) = (z(j)-z(i))/rij

        bm(irow,3*i-2) = -eij(1)
        bm(irow,3*i-1) = -eij(2)
        bm(irow,3*i)   = -eij(3)
        bm(irow,3*j-2) = eij(1)
        bm(irow,3*j-1) = eij(2)
        bm(irow,3*j)   = eij(3)
    enddo

    ! Evaluate the B matrix contributed from bendings
    do ijk=1,nba
        irow = irow+1
        i    = iba(1,ijk)
        j    = iba(2,ijk)
        k    = iba(3,ijk)

        rji  = bndlen(j,i,x,y,z)
        rjk  = bndlen(j,k,x,y,z)
        aijk = bndang(i,j,k,x,y,z)
        cijk = cos(aijk)
        sijk = sin(aijk)

        eji(1) = (x(i)-x(j))/rji
        eji(2) = (y(i)-y(j))/rji
        eji(3) = (z(i)-z(j))/rji
        ejk(1) = (x(k)-x(j))/rjk
        ejk(2) = (y(k)-y(j))/rjk
        ejk(3) = (z(k)-z(j))/rjk

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

    !
    ! Evaluate the B matrix contributed from torsions
    ! The calculations follow I. H. Williams, J. Mol. Spectros. 66, 288 (197
    !
    do ijkl=1,nto
        irow = irow+1
        i    = ito(1,ijkl)
        j    = ito(2,ijkl)
        k    = ito(3,ijkl)
        l    = ito(4,ijkl)

        rij  = bndlen(i,j,x,y,z)
        rjk  = bndlen(j,k,x,y,z)
        rkl  = bndlen(k,l,x,y,z)
        aijk = bndang(i,j,k,x,y,z)
        ajkl = bndang(j,k,l,x,y,z)
        cijk = cos(aijk)
        sijk = sin(aijk)
        cjkl = cos(ajkl)
        sjkl = sin(ajkl)

        eij(1) = (x(j)-x(i))/rij
        eij(2) = (y(j)-y(i))/rij
        eij(3) = (z(j)-z(i))/rij
        ejk(1) = (x(k)-x(j))/rjk
        ejk(2) = (y(k)-y(j))/rjk
        ejk(3) = (z(k)-z(j))/rjk
        ekl(1) = (x(l)-x(k))/rkl
        ekl(2) = (y(l)-y(k))/rkl
        ekl(3) = (z(l)-z(k))/rkl
        !
        ! Caculate u=eij*ejk, v=ejk*ekl
        !
        call xprod1(eij,ejk,u)
        call xprod1(ejk,ekl,v)
        !
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

    !
    ! Evaluate the B matrix contributed from out-of-plane bending
    ! The calculation follows D. F. Mclntosh and K. H. Michaelian,
    ! Can. J. Spectros. 24, 35-40 (1979).
    !
    do ijkl=1,noop
        irow = irow+1
        i    = ioop(1,ijkl)
        j    = ioop(2,ijkl)
        k    = ioop(3,ijkl)
        l    = ioop(4,ijkl)

        rji    = bndlen(j,i,x,y,z)
        rjk    = bndlen(j,k,x,y,z)
        rjl    = bndlen(j,l,x,y,z)
        ai     = bndang(k,j,l,x,y,z)
        ak     = bndang(i,j,l,x,y,z)
        al     = bndang(i,j,k,x,y,z)
        theta  = oopba1(i,j,k,l,x,y,z)
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

        !
        ! Caculate u=ejk*ejl
        !
        call xprod1(ejk,ejl,u)
        !
        bm(irow,3*i-2) = (1.0d0/rji)*(u(1)/(ctheta*si)-ttheta*eji(1))
        bm(irow,3*i-1) = (1.0d0/rji)*(u(2)/(ctheta*si)-ttheta*eji(2))
        bm(irow,3*i)   = (1.0d0/rji)*(u(3)/(ctheta*si)-ttheta*eji(3))
        bm(irow,3*k-2) = (1.0d0/rjk)*u(1)*(ci*ck-cl)/(ctheta*si*si*si)
        bm(irow,3*k-1) = (1.0d0/rjk)*u(2)*(ci*ck-cl)/(ctheta*si*si*si)
        bm(irow,3*k)   = (1.0d0/rjk)*u(3)*(ci*ck-cl)/(ctheta*si*si*si)
        bm(irow,3*l-2) = (1.0d0/rjl)*u(1)*(ci*cl-ck)/(ctheta*si*si*si)
        bm(irow,3*l-1) = (1.0d0/rjl)*u(2)*(ci*cl-ck)/(ctheta*si*si*si)
        bm(irow,3*l)   = (1.0d0/rjl)*u(3)*(ci*cl-ck)/(ctheta*si*si*si)
        bm(irow,3*j-2) = -bm(irow,3*i-2)-bm(irow,3*k-2)-bm(irow,3*l-2)
        bm(irow,3*j-1) = -bm(irow,3*i-1)-bm(irow,3*k-1)-bm(irow,3*l-1)
        bm(irow,3*j)   = -bm(irow,3*i)  -bm(irow,3*k)  -bm(irow,3*l)
    enddo

end subroutine bmat3
!> @brief Evaluate matmult for the PHOH_APRP_DPEM potential.
subroutine matmult( c, ldc, a, lda, m, b, ldb, n, l )

    ! Global variables
    implicit none

    integer          :: lda,ldb,ldc
    integer          :: m,n,l
    double precision :: a(lda,*),b(ldb,*),c(ldc,*)

    ! Local variables
    integer :: i,j,k

    do i=1,m
        do j=1,l
            c(i,j) = 0.0d0
        enddo
    enddo

    do i=1,m
        do j=1,l
            do k=1,n
                c(i,j) = c(i,j)+a(i,k)*b(k,j)
            enddo
        enddo
    enddo

end subroutine matmult
!> @brief Evaluate bndlen for the PHOH_APRP_DPEM potential.
double precision function bndlen( i, j, x, y, z )
    !***********************************************************************
    ! Function to evaluate the bond length between atom i and j
    ! Input: i,j (indices of 2 atoms)
    !***********************************************************************

    implicit none

    integer :: i, j

    double precision, parameter :: tiny = 1.0d-13

    double precision :: x(*),y(*),z(*)

    bndlen = sqrt( (x(i)-x(j))*(x(i)-x(j)) + (y(i)-y(j))*(y(i)-y(j)) + (z(i)-z(j))*(z(i)-z(j)))

    if (bndlen .le. tiny) then
        write(*,'('' Atom'',I3,'' and'',I3,'' are too close to each other!'')') i,j
        !DB  KRY 20140718
        write(*,'(I3,3F18.12)') i,x(i),y(i),z(i)
        write(*,'(I3,3F18.12)') j,x(j),y(j),z(j)
        stop
    endif

end function bndlen
!> @brief Evaluate bndang for the PHOH_APRP_DPEM potential.
double precision function bndang( i, j, k, x, y, z )
    !***********************************************************************
    ! Evaluate the bond angle between atoms i, j, and k
    ! Input: i,j,k (indices of 3 atoms, i-j and j-k bonds)
    !***********************************************************************

    implicit none

    double precision :: a, bndlen
    integer          :: i, j, k

    double precision :: x(*),y(*),z(*)
    double precision :: eji(3),ejk(3),rji,rjk

    rji = bndlen(i,j,x,y,z)
    rjk = bndlen(j,k,x,y,z)

    eji(1) = (x(i)-x(j))/rji
    eji(2) = (y(i)-y(j))/rji
    eji(3) = (z(i)-z(j))/rji

    ejk(1) = (x(k)-x(j))/rjk
    ejk(2) = (y(k)-y(j))/rjk
    ejk(3) = (z(k)-z(j))/rjk

    call sprod(3,eji,ejk,a)
    bndang = acos(a)

end function bndang
!> @brief Evaluate dihedr for the PHOH_APRP_DPEM potential.
double precision function dihedr( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Evaluate the dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************

    implicit none

    double precision :: a, bndang, bndlen
    integer          :: i, j, k, l

    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0

    double precision :: x(*),y(*),z(*)
    double precision :: eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
    double precision :: rij,rjk,rkl,aijk,ajkl,wejk

    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    rij  = bndlen(i,j,x,y,z)
    rjk  = bndlen(j,k,x,y,z)
    rkl  = bndlen(k,l,x,y,z)
    aijk = bndang(i,j,k,x,y,z)
    ajkl = bndang(j,k,l,x,y,z)

    ! Check whether i-j-k or j-k-l are collinear
    if ((abs(aijk) .le. tiny) .or. (abs(aijk-pi) .lt. tiny)) then
        write(*,'(''Torsion angle not defined since atoms'',3I3,'' collinear'')') i,j,k
    endif
    if ((abs(ajkl) .le. tiny) .or. (abs(ajkl-pi) .lt. tiny)) then
        write(*,'(''Torsion angle not defined since atoms'',3I3,'' collinear'')') j,k,l
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

    !
    ! u = eij*ejk normal vector in plane i-j-k
    ! v = ejk*ekl Normal vector in plane j-k-l
    ! w = u*v used to determine the sign of dihedral angle
    !
    call xprod1(eij,ejk,u)
    call xprod1(ejk,ekl,v)
    call xprod1(u,v,w)

    call sprod(3,u,v,a)

    a      = a/(dsin(aijk)*dsin(ajkl))
    a      = max(a,-1.0d0)
    a      = min(a,1.0d0)
    dihedr = dacos(a)

    ! Determine the sign of dihedr based on u,v,ejk
    ! dihedr >= 0, u,v,ejk right-hand (w*ejk >=0)
    ! dihedr < 0,  u,v,ejk left-hand (w*ejk < 0)

    wejk = w(1)*ejk(1)+w(2)*ejk(2)+w(3)*ejk(3)
    if (wejk .lt. 0.0d0) then
        dihedr = -dihedr
    endif

    ! Setting the range of dihedral angle between -90.0 degree to 270 degree
    ! to avoid discontinuity around 180 degrees
    !      if (dihedr .lt. -0.5d0*pi) then
    !        dihedr=2.0d0*pi+dihedr
    !      endif

end function dihedr
!> @brief Evaluate oopba1 for the PHOH_APRP_DPEM potential.
double precision function oopba1( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Evaluate the out-of-plane bending angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !             l
    !            /
    !     i-----j
    !            \
    !             k
    !***********************************************************************

    implicit none

    double precision :: bndang, bndlen
    integer          :: i, j, k, l

    double precision :: tiny
    data tiny/1.0d-30/

    double precision :: x(*),y(*),z(*)
    double precision :: eji(3),ejk(3),ejl(3),u(3)
    double precision :: v(3),w(3)
    double precision :: a1,a2
    double precision :: a,ai,si,rji,rjk,rjl

    !               l
    !              /
    !          ak /
    !     i------j ai != 180 deg
    !          al \
    !              \
    !               k
    rji    = bndlen(j,i,x,y,z)
    eji(1) = (x(i)-x(j))/rji
    eji(2) = (y(i)-y(j))/rji
    eji(3) = (z(i)-z(j))/rji

    rjk    = bndlen(j,k,x,y,z)
    ejk(1) = (x(k)-x(j))/rjk
    ejk(2) = (y(k)-y(j))/rjk
    ejk(3) = (z(k)-z(j))/rjk

    rjl    = bndlen(j,l,x,y,z)
    ejl(1) = (x(l)-x(j))/rjl
    ejl(2) = (y(l)-y(j))/rjl
    ejl(3) = (z(l)-z(j))/rjl

    ai = bndang(k,j,l,x,y,z)

    si = sin(ai)

    call xprod1(ejk,ejl,u)
    ! Check if whether ai = 180
    if (u(1)**2+u(2)**2+u(3)**2 .le. tiny) then
        write(*,*) "Looks like k-j-l collinear."
        write(*,*) "Out-of-plane bending angle is not defined."
        stop
    endif

    call sprod(3,u,eji,a)
    a      = a/si
    a      = max(a,-1.0d0)
    a      = min(a,1.0d0)
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
!> @brief Evaluate xprod1 for the PHOH_APRP_DPEM potential.
subroutine xprod1( x1, x2, x3 )
    !***********************************************************************
    ! Calculate cross-product x1 and x2 and return x3
    ! Input: X1,X2
    ! Output: X3
    !***********************************************************************

    implicit none

    double precision :: x1(3),x2(3),x3(3)

    x3(1) = x1(2)*x2(3)-x1(3)*x2(2)
    x3(2) = x1(3)*x2(1)-x1(1)*x2(3)
    x3(3) = x1(1)*x2(2)-x1(2)*x2(1)

end subroutine xprod1
!> @brief Evaluate sprod for the PHOH_APRP_DPEM potential.
subroutine sprod( n, a, b, c )
    !***********************************************************************
    ! Calculate scalar-product a and b and return c
    ! Input: a(n),b(n)
    ! Output: c
    !***********************************************************************

    implicit none

    integer :: i

    integer          :: n
    double precision :: a(*),b(*),c

    c = 0.0d0
    do i=1,n
        c = c+a(i)*b(i)
    enddo

end subroutine sprod
