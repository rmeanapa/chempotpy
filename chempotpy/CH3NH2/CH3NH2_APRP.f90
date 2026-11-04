!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!   System:          CH3NH2
!   Functional form: Anchor-Points Reactive Potential (APRP)
!                               (2x2 diabatic potential energy matrix)
!   Common name:                   Methylamine(Coupled PESs)
!   Number of derivatives:         1
!   Number of bodies:              7
!   Number of electronic surfaces: 2
!   Interface:                     Section-2
!
!   Reference: Kelsey A. Parker,and Donald G. Truhlar
!              J. Chem. Phys. 152, 244309 (2020).
!
!   Notes:
!    1. This routine calculates two coupled singlet PESs for
!       nonadiabatic NH photodissociation of CH3NH2 to produce
!       CH3NH + H. It calculates the diabatic potential energy matrix
!       and its gradient and the adiabatic potential energy surfaces
!       and their gradients, and it calculates the space-frame
!       nonadiabatic couling vector.
!    2. Only the first H is allowed to dissociate.
!    3. The final version of this potential was added to POTLIB
!       in March, 2020.
!
!   Numbering of atoms:
!       H5         H1
!         \       /
!       H6--C4---N3
!         /       \
!       H7         H2
!
!   Units:      input:   angstrom for coordinates
!               output:  hartree for energy
!                        hartree/bohr for gradients
!                        1/bohr for nonadiabatic coupling vector
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
!   guu(1,i,j,j)        X component of the gradients of diabatic
!                       potential of state j at atom i
!   guu(2,i,j,j)        Y component of the gradients of diabatic
!                       potential of state j at atom i
!   guu(3,i,j,j)        Z component of the gradients of diabatic
!                       potential of state j at atom i
!   guu(1,i,j,k)        X component of the gradients of diabatic
!                       coupling between state j and k at atom i
!   guu(2,i,j,k)        Y component of the gradients of diabatic
!                       coupling between state j and k at atom i
!   guu(3,i,j,k)        Z component of the gradients of diabatic
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
!   eigvec(2,2)         2*2 orthogonal matrix diagonalizing UU matrix
!                       (CC^T*UU*CC yield diagonal matrix with adiabatic
!                       energies as diagonal elements)
!
!   Note:       LAPACK library is needed to diagonize diabatic potential
!               matrix (UU) to yield adiabatic potenitials VV(2)
!
!***********************************************************************

!> @brief Evaluate pes for the CH3NH2_APRP potential.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
    implicit none
    ! number of electronic state
    integer, parameter            :: nstates = 2
    integer, parameter            :: natoms  = 7
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    double precision              :: xx(3,natoms), uu(nstates,nstates)
    double precision              :: guu(3,natoms,nstates,nstates)
    double precision              :: vv(nstates), gvv(3,natoms,nstates)
    double precision              :: dvec(3,natoms,nstates,nstates)
    double precision              :: cc(nstates,nstates)
    integer                       :: iatom, idir, j, istate, jstate
    integer                       :: jgrad
    !initialize
    p = 0.d0; g = 0.d0; d = 0.d0
    do iatom=1, natoms
        do idir=1, 3
            xx(idir,iatom) = x(iatom,idir)  !/legacy_bohr_to_angstrom
        enddo
    enddo
    if( igrad == 0 )then
        jgrad = 0
    else
        jgrad = 1
    endif
    call pot(jgrad, xx, uu, guu, vv, gvv, dvec, cc, 0)
    if( igrad == 0 )then
        do istate=1, nstates
            p(istate) = vv(istate)*legacy_hartree_to_electronvolt
        enddo
    else if( igrad == 1 )then
        do istate=1, nstates
            p(istate) = vv(istate)*legacy_hartree_to_electronvolt
        enddo
        do istate=1, nstates
            do iatom=1, natoms
                do idir=1, 3
                    g(istate,iatom,idir) = gvv(idir,iatom,istate)*51.422067
                enddo
            enddo
        enddo
    else if( igrad == 2 )then
        do istate=1, nstates
            p(istate) = vv(istate)*legacy_hartree_to_electronvolt
        enddo
        do istate=1, nstates
            do iatom=1, natoms
                do idir=1, 3
                    g(istate,iatom,idir) = gvv(idir,iatom,istate)*51.422067
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
end subroutine pes

!> @brief Evaluate pot for the CH3NH2_APRP potential.
subroutine pot( igrad, xx, uu, guu, vv, gvv, dvec, eigvec, repflag )
    implicit none
    double precision              :: xx(3,7)
    double precision              :: x(7), y(7), z(7)
    double precision              :: scp(3)
    integer                       :: i, j
    double precision              :: e1
    integer                       :: nrap, ntap
    double precision              :: dx1(7), dy1(7), dz1(7)
    double precision              :: dx(7), dy(7), dz(7)
    integer                       :: igrad
    double precision              :: dux(7,3), duy(7,3), duz(7,3)
    double precision              :: uu(2,2), vv(2)
    double precision              :: tt(1)
    double precision              :: guu(3,7,2,2), gtt(3,7,1,1)
    double precision              :: gvv(3,7,2)
    double precision              :: dvec(3,7,2,2)
    double precision              :: eigvec(2,2)
    integer                       :: lwork, info
    integer                       :: repflag
    double precision, allocatable :: work(:)
    double precision              :: tmpmat22(2,2)
    integer                       :: ist, jst
    double precision, parameter   :: evtohartree = 27.21139d0
    double precision, parameter   :: bohrtoang   = 0.52917721067d0
    double precision              :: v2, v4
    scp = 0d0
    dux = 0d0
    duy = 0d0
    duz = 0d0
    do i=1, 7
        x(i) = xx(1,i)
        y(i) = xx(2,i)
        z(i) = xx(3,i)
    enddo
    ! get tertiary then primary contributions
    do j=1, 3
        nrap = 4
        ntap = 3
        call tert(igrad, j, x, y, z, nrap, ntap, e1, dx1, dy1, dz1)
        scp(j) = e1
        do i=1, 7
            dux(i,j) = dx1(i)
            duy(i,j) = dy1(i)
            duz(i,j) = dz1(i)
        enddo
        call kpev2gm2(igrad, x, y, z, v2, j, dx, dy, dz)
        scp(j) = scp(j)+v2
        do i=1, 7
            dux(i,j) = dux(i,j)+dx(i)
            duy(i,j) = duy(i,j)+dy(i)
            duz(i,j) = duz(i,j)+dz(i)
        enddo
    enddo
    ! change units from eV to Hartree and Angstrom to Bohr
    uu(1,1) = scp(1)/evtohartree
    uu(2,2) = scp(2)/evtohartree
    uu(1,2) = scp(3)/evtohartree
    uu(2,1) = scp(3)/evtohartree
    eigvec  = uu
    lwork   = -1
    allocate (work (1) )
    call dsyev('V', 'L', 2, eigvec, 2, vv, work, lwork, info)
    lwork = int(work(1))
    deallocate (work)
    allocate( work(lwork) )
    call dsyev('V', 'L', 2, eigvec, 2, vv, work, lwork, info)
    deallocate (work)
    if( info .ne. 0 )then
        write(6, *) 'Error in diagolization of U matrix'
        stop
    endif
    do i=1, 7
        guu(1,i,1,1) = (dux(i,1)/evtohartree)*bohrtoang
        guu(2,i,1,1) = (duy(i,1)/evtohartree)*bohrtoang
        guu(3,i,1,1) = (duz(i,1)/evtohartree)*bohrtoang
        !
        guu(1,i,2,2) = (dux(i,2)/evtohartree)*bohrtoang
        guu(2,i,2,2) = (duy(i,2)/evtohartree)*bohrtoang
        guu(3,i,2,2) = (duz(i,2)/evtohartree)*bohrtoang
        !
        guu(1,i,1,2) = (dux(i,3)/evtohartree)*bohrtoang
        guu(2,i,1,2) = (duy(i,3)/evtohartree)*bohrtoang
        guu(3,i,1,2) = (duz(i,3)/evtohartree)*bohrtoang
        !
        guu(1,i,2,1) = (dux(i,3)/evtohartree)*bohrtoang
        guu(2,i,2,1) = (duy(i,3)/evtohartree)*bohrtoang
        guu(3,i,2,1) = (duz(i,3)/evtohartree)*bohrtoang
    enddo
    ! convert diabatic potential energy matrix and its gradients to
    ! adiabatic potentials and gradients
    do i=1, 3
        do j=1, 7
            tmpmat22(:,:) = guu(i,j,:,:)
            tmpmat22      = matmul(transpose(eigvec), matmul(tmpmat22, eigvec))
            do ist=1, 2
                gvv(i,j,ist) = tmpmat22(ist,ist)
            enddo
            !
            ist               = 1
            jst               = 2
            dvec(i,j,ist,jst) = tmpmat22(ist,jst)/ (vv(jst)-vv(ist))
            dvec(i,j,jst,ist) = -dvec(i,j,ist,jst)
        enddo
    enddo
end subroutine pot

!> @brief Evaluate dihedr1 for the CH3NH2_APRP potential.
double precision function dihedr1( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Evaluate the dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    ! I created a modified version of this code to avoid issues when the
    ! amine is flat in a plane
    !***********************************************************************
    implicit none
    double precision            :: a, bndang, bndlen
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
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
    !
    ! u = eij*ejk normal vector in plane i-j-k
    ! v = ejk*ekl Normal vector in plane j-k-l
    ! w = u*v used to determine the sign of dihedral angle
    !
    call xprod1(eij, ejk, u)
    call xprod1(ejk, ekl, v)
    call xprod1(u, v, w)
    call sprod(3, u, v, a)
    a       = a/(dsin(aijk)*dsin(ajkl))
    a       = max(a, -1.0d0)
    a       = min(a, 1.0d0)
    dihedr1 = dacos(a)
    ! Determine the sign of dihedr based on u,v,ejk
    ! dihedr >= 0, u,v,ejk right-hand (w*ejk >=0)
    ! dihedr < 0,  u,v,ejk left-hand (w*ejk < 0)
    wejk = w(1)*ejk(1)+w(2)*ejk(2)+w(3)*ejk(3)
    if( wejk .lt. 0.0d0 )then
        dihedr1 = -dihedr1
    endif
    ! Setting the range of dihedral angle between -90.0 degree to 270 degree
    ! to avoid discontinuity around 180 degrees
    !      if (dihedr .lt. -0.5d0*pi) then
    !        dihedr=2.0d0*pi+dihedr
    !      endif
end function dihedr1

!> @brief Evaluate bndlen for the CH3NH2_APRP potential.
double precision function bndlen( i, j, x, y, z )
    !***********************************************************************
    !* Function to evaluate the bond length between atom i and j
    !* Input: i,j (indices of 2 atoms)
    !***********************************************************************
    implicit none
    integer                     :: i, j
    double precision, parameter :: tiny = 1.0d-13
    double precision            :: x(*), y(*), z(*)
    bndlen = sqrt( (x(i)-x(j))*(x(i)-x(j)) +      (y(i)-y(j))*(y(i)-y(j)) +      (z(i)-z(j))*(z(i)-z(j)))
    if( bndlen .le. tiny )then
        write(*, '('' Atom'',I3,'' and'',I3,'' are too close to each other!'')') i, j
        !CDB  KRY 20140718
        write(*, '(I3,3F18.12)') i, x(i), y(i), z(i)
        write(*, '(I3,3F18.12)') j, x(j), y(j), z(j)
        stop
    endif
end function bndlen

!> @brief Evaluate bndang for the CH3NH2_APRP potential.
double precision function bndang( i, j, k, x, y, z )
    !***********************************************************************
    !* Evaluate the bond angle between atoms i, j, and k
    !* Input: i,j,k (indices of 3 atoms, i-j and j-k bonds)
    !***********************************************************************
    implicit none
    double precision :: a, bndlen
    integer          :: i, j, k
    double precision :: x(*), y(*), z(*)
    double precision :: eji(3), ejk(3), rji, rjk
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

!> @brief Evaluate xprod1 for the CH3NH2_APRP potential.
subroutine xprod1( x1, x2, x3 )
    !***********************************************************************
    !* Calculate cross product of x1 and x2 and return x3
    !* Input: X1,X2
    !* Output: X3
    !***********************************************************************
    implicit none
    double precision :: x1(3), x2(3), x3(3)
    x3(1) = x1(2)*x2(3)-x1(3)*x2(2)
    x3(2) = x1(3)*x2(1)-x1(1)*x2(3)
    x3(3) = x1(1)*x2(2)-x1(2)*x2(1)
end subroutine xprod1

!> @brief Evaluate sprod for the CH3NH2_APRP potential.
subroutine sprod( n, a, b, c )
    !***********************************************************************
    !* Calculate scalar product of a and b and return c
    !* Input: a(n),b(n)
    !* Output: c
    !***********************************************************************
    implicit none
    integer          :: i
    integer          :: n
    double precision :: a(*), b(*), c
    c = 0.0d0
    do i=1, n
        c = c+a(i)*b(i)
    enddo
end subroutine sprod

!> @brief Evaluate convertit for the CH3NH2_APRP potential.
subroutine convertit( x, y, z, dk, dx1, dy1, dz1 )
    ! this subroutine deals with the tertiary derivatives and converts these
    ! from the various internal coords to XYZ
    implicit none
    integer          :: kk, nat
    integer          :: i, j, k
    double precision :: xyz(4,3)
    double precision :: beta
    double precision :: x(7), y(7), z(7)
    double precision :: dk(22)
    double precision :: r1, t2
    double precision :: qc(20)
    double precision :: dx1(7), dy1(7), dz1(7)
    double precision :: bndlen
    double precision :: dngx(7), dngy(7), dngz(7)
    double precision :: dihedrix, dihedriy, dihedriz
    double precision :: dihedrjx, dihedrjy, dihedrjz
    double precision :: dihedrkx, dihedrky, dihedrkz
    double precision :: dihedrlx, dihedrly, dihedrlz
    double precision :: dihedr2
    double precision :: torh1, torh2, rh1hx, rh2hx
    double precision :: ang1, ang2, ang3
    double precision :: len1, len2, len3
    integer          :: bonds(2,8), bends(3,8)
    data bonds   /2, 3, &
        & 3, 4, &
        & 4, 5, &
        & 4, 6, &
        & 4, 7, &
        !
        & 5, 6, &
        & 6, 7, &
        & 7, 5/
    data bends   /1, 3, 4, &
        & 2, 3, 4, &
        & 3, 4, 5, &
        & 3, 4, 6, &
        & 3, 4, 7, &
        !
        & 5, 4, 6, &
        & 6, 4, 7, &
        & 7, 4, 5/
    r1 = bndlen(1, 3, x, y, z)
    t2 = dihedr2(1, 3, 4, 2, x, y, z)
    call calcgeom(x, y, z, qc)
    torh1 = qc(11)
    torh2 = qc(12)
    rh1hx = qc(13)
    rh2hx = qc(14)
    ang1  = qc(15)
    ang2  = qc(16)
    ang3  = qc(17)
    len1  = qc(18)
    len2  = qc(19)
    len3  = qc(20)
    dx1   = 0d0
    dy1   = 0d0
    dz1   = 0d0
    do kk=1, 5
        i = bonds(1,kk)
        j = bonds(2,kk)
        call disder(x, y, z, i, j, dngx, dngy, dngz)
        do nat=1, 7
            dx1(nat) = dx1(nat)+dk(kk)*dngx(nat)
            dy1(nat) = dy1(nat)+dk(kk)*dngy(nat)
            dz1(nat) = dz1(nat)+dk(kk)*dngz(nat)
        enddo
    enddo
    do kk=1, 5
        i = bends(1,kk)
        j = bends(2,kk)
        k = bends(3,kk)
        call angder(x, y, z, i, j, k, dngx, dngy, dngz)
        do nat=1, 7
            dx1(nat) = dx1(nat)+dk(kk+5)*dngx(nat)
            dy1(nat) = dy1(nat)+dk(kk+5)*dngy(nat)
            dz1(nat) = dz1(nat)+dk(kk+5)*dngz(nat)
        enddo
    enddo
    ! 11th/12th dk COORD
    dx1(1) = dx1(1)+dk(11)*((x(1)-x(5))/rh1hx)
    dy1(1) = dy1(1)+dk(11)*((y(1)-y(5))/rh1hx)
    dz1(1) = dz1(1)+dk(11)*((z(1)-z(5))/rh1hx)
    dx1(5) = dx1(5)+dk(11)*((x(5)-x(1))/rh1hx)
    dy1(5) = dy1(5)+dk(11)*((y(5)-y(1))/rh1hx)
    dz1(5) = dz1(5)+dk(11)*((z(5)-z(1))/rh1hx)
    !!!!!!!!
    dx1(2) = dx1(2)+dk(12)*((x(2)-x(5))/rh2hx)
    dy1(2) = dy1(2)+dk(12)*((y(2)-y(5))/rh2hx)
    dz1(2) = dz1(2)+dk(12)*((z(2)-z(5))/rh2hx)
    dx1(5) = dx1(5)+dk(12)*((x(5)-x(2))/rh2hx)
    dy1(5) = dy1(5)+dk(12)*((y(5)-y(2))/rh2hx)
    dz1(5) = dz1(5)+dk(12)*((z(5)-z(2))/rh2hx)
    !!!!
    ! 13th/14th dk COORD
    !!!!
    ! H1
    dx1(1) = dx1(1)+dk(13)*dihedrix(1, 3, 4, 5, x, y, z)
    dy1(1) = dy1(1)+dk(13)*dihedriy(1, 3, 4, 5, x, y, z)
    dz1(1) = dz1(1)+dk(13)*dihedriz(1, 3, 4, 5, x, y, z)
    ! N
    dx1(3) = dx1(3)+dk(13)*dihedrjx(1, 3, 4, 5, x, y, z)
    dy1(3) = dy1(3)+dk(13)*dihedrjy(1, 3, 4, 5, x, y, z)
    dz1(3) = dz1(3)+dk(13)*dihedrjz(1, 3, 4, 5, x, y, z)
    ! C
    dx1(4) = dx1(4)+dk(13)*dihedrkx(1, 3, 4, 5, x, y, z)
    dy1(4) = dy1(4)+dk(13)*dihedrky(1, 3, 4, 5, x, y, z)
    dz1(4) = dz1(4)+dk(13)*dihedrkz(1, 3, 4, 5, x, y, z)
    ! HX
    dx1(5) = dx1(5)+dk(13)*dihedrlx(1, 3, 4, 5, x, y, z)
    dy1(5) = dy1(5)+dk(13)*dihedrly(1, 3, 4, 5, x, y, z)
    dz1(5) = dz1(5)+dk(13)*dihedrlz(1, 3, 4, 5, x, y, z)
    !!!!
    ! H2
    dx1(2) = dx1(2)+dk(14)*dihedrix(2, 3, 4, 5, x, y, z)
    dy1(2) = dy1(2)+dk(14)*dihedriy(2, 3, 4, 5, x, y, z)
    dz1(2) = dz1(2)+dk(14)*dihedriz(2, 3, 4, 5, x, y, z)
    ! N
    dx1(3) = dx1(3)+dk(14)*dihedrjx(2, 3, 4, 5, x, y, z)
    dy1(3) = dy1(3)+dk(14)*dihedrjy(2, 3, 4, 5, x, y, z)
    dz1(3) = dz1(3)+dk(14)*dihedrjz(2, 3, 4, 5, x, y, z)
    ! C
    dx1(4) = dx1(4)+dk(14)*dihedrkx(2, 3, 4, 5, x, y, z)
    dy1(4) = dy1(4)+dk(14)*dihedrky(2, 3, 4, 5, x, y, z)
    dz1(4) = dz1(4)+dk(14)*dihedrkz(2, 3, 4, 5, x, y, z)
    ! HX
    dx1(5) = dx1(5)+dk(14)*dihedrlx(2, 3, 4, 5, x, y, z)
    dy1(5) = dy1(5)+dk(14)*dihedrly(2, 3, 4, 5, x, y, z)
    dz1(5) = dz1(5)+dk(14)*dihedrlz(2, 3, 4, 5, x, y, z)
    do kk=6, 8
        i = bends(1,kk)
        j = bends(2,kk)
        k = bends(3,kk)
        call angder(x, y, z, i, j, k, dngx, dngy, dngz)
        do nat=1, 7
            dx1(nat) = dx1(nat)+dk(9+kk)*dngx(nat)
            dy1(nat) = dy1(nat)+dk(9+kk)*dngy(nat)
            dz1(nat) = dz1(nat)+dk(9+kk)*dngz(nat)
        enddo
    enddo
    do kk=6, 8
        i = bonds(1,kk)
        j = bonds(2,kk)
        call disder(x, y, z, i, j, dngx, dngy, dngz)
        do nat=1, 7
            dx1(nat) = dx1(nat)+dk(12+kk)*dngx(nat)
            dy1(nat) = dy1(nat)+dk(12+kk)*dngy(nat)
            dz1(nat) = dz1(nat)+dk(12+kk)*dngz(nat)
        enddo
    enddo
    ! 21th dk COORD
    dx1(1) = dx1(1)+dk(21)*(x(1)-x(3))/r1
    dy1(1) = dy1(1)+dk(21)*(y(1)-y(3))/r1
    dz1(1) = dz1(1)+dk(21)*(z(1)-z(3))/r1
    dx1(3) = dx1(3)+dk(21)*(x(3)-x(1))/r1
    dy1(3) = dy1(3)+dk(21)*(y(3)-y(1))/r1
    dz1(3) = dz1(3)+dk(21)*(z(3)-z(1))/r1
    ! 22th dk COORD
    dx1(1) = dx1(1)+dk(22)*dihedrix(1, 3, 4, 2, x, y, z)
    dy1(1) = dy1(1)+dk(22)*dihedriy(1, 3, 4, 2, x, y, z)
    dz1(1) = dz1(1)+dk(22)*dihedriz(1, 3, 4, 2, x, y, z)
    dx1(2) = dx1(2)+dk(22)*dihedrlx(1, 3, 4, 2, x, y, z)
    dy1(2) = dy1(2)+dk(22)*dihedrly(1, 3, 4, 2, x, y, z)
    dz1(2) = dz1(2)+dk(22)*dihedrlz(1, 3, 4, 2, x, y, z)
    ! N
    dx1(3) = dx1(3)+dk(22)*dihedrjx(1, 3, 4, 2, x, y, z)
    dy1(3) = dy1(3)+dk(22)*dihedrjy(1, 3, 4, 2, x, y, z)
    dz1(3) = dz1(3)+dk(22)*dihedrjz(1, 3, 4, 2, x, y, z)
    ! C
    dx1(4) = dx1(4)+dk(22)*dihedrkx(1, 3, 4, 2, x, y, z)
    dy1(4) = dy1(4)+dk(22)*dihedrky(1, 3, 4, 2, x, y, z)
    dz1(4) = dz1(4)+dk(22)*dihedrkz(1, 3, 4, 2, x, y, z)
end subroutine convertit

!> @brief Evaluate angder for the CH3NH2_APRP potential.
subroutine angder( x, y, z, i, j, k, dngx, dngy, dngz )
    ! This subroutine converts a bond angle derivative to XYZ
    implicit none
    integer          :: i, j, k
    double precision :: x(7), y(7), z(7)
    double precision :: ra, rb
    double precision :: to, be, go, cgo
    double precision :: dngx(7), dngy(7), dngz(7)
    double precision :: bndlen
    dngx = 0d0
    dngy = 0d0
    dngz = 0d0
    ra   = bndlen(j, i, x, y, z)
    rb   = bndlen(j, k, x, y, z)
    to   = (x(j)-x(i))*(x(j)-x(k)) +(y(j)-y(i))*(y(j)-y(k)) +(z(j)-z(i))*(z(j)-z(k))
    be=((((x(j)-x(i))**2)+((y(j)-y(i))**2)+((z(j)-z(i))**2))**0.5d0) *((((x(j)-x(k))**2)+((y(j)-y(k))**2)+((z(j)- &
        & z(k))** 2))**0.5d0)
    go  = to/be
    cgo = (-1d0/((1d0-go**2d0)**(1d0/2d0)))
    !--------- j -----------!
    dngx(j) = cgo*((2d0*x(j)-x(i)-x(k))*be -( (((x(j)-x(i))/ra)*rb) + (((x(j)-x(k))/rb)*ra) )*to)/(be**2d0)
    dngy(j) = cgo*((2d0*y(j)-y(i)-y(k))*be -( (((y(j)-y(i))/ra)*rb) + (((y(j)-y(k))/rb)*ra) )*to)/(be**2d0)
    dngz(j) = cgo*((2d0*z(j)-z(i)-z(k))*be -( (((z(j)-z(i))/ra)*rb) + (((z(j)-z(k))/rb)*ra) )*to)/(be**2d0)
    !--------- k -----------!
    dngx(k) = cgo*( -1d0*(x(j)-x(i))*be -(((x(k)-x(j))/rb)*ra*to))/(be**2d0)
    dngy(k) = cgo*( -1d0*(y(j)-y(i))*be -(((y(k)-y(j))/rb)*ra*to))/(be**2d0)
    dngz(k) = cgo*( -1d0*(z(j)-z(i))*be -(((z(k)-z(j))/rb)*ra*to))/(be**2d0)
    !--------- i -----------!
    dngx(i) = cgo*( -(x(j)-x(k))*be -(((x(i)-x(j))/ra)*rb*to))/(be**2d0)
    dngy(i) = cgo*( -(y(j)-y(k))*be -(((y(i)-y(j))/ra)*rb*to))/(be**2d0)
    dngz(i) = cgo*( -(z(j)-z(k))*be -(((z(i)-z(j))/ra)*rb*to))/(be**2d0)
end subroutine angder

!> @brief Evaluate disder for the CH3NH2_APRP potential.
subroutine disder( x, y, z, i, j, dngx, dngy, dngz )
    ! This subroutine converts a distance derivative to XYZ
    implicit none
    integer          :: i, j, k
    double precision :: x(7), y(7), z(7)
    double precision :: ra, rb
    double precision :: to, be, go, cgo
    double precision :: dngx(7), dngy(7), dngz(7)
    double precision :: bndlen
    dngx    = 0d0
    dngy    = 0d0
    dngz    = 0d0
    dngx(i) = (x(i)-x(j))/bndlen(i, j, x, y, z)
    dngy(i) = (y(i)-y(j))/bndlen(i, j, x, y, z)
    dngz(i) = (z(i)-z(j))/bndlen(i, j, x, y, z)
    dngx(j) = (x(j)-x(i))/bndlen(i, j, x, y, z)
    dngy(j) = (y(j)-y(i))/bndlen(i, j, x, y, z)
    dngz(j) = (z(j)-z(i))/bndlen(i, j, x, y, z)
end subroutine disder

!> @brief Evaluate solver for the CH3NH2_APRP potential.
subroutine solver( igrad, ra1, ra2, qc, ii, nrap, ntap, ut, dut )
    !***********************************************************************
    ! Subroutine to evaluate the tertiary contributions to the
    ! diabatic potentials and gradients at different anchor points
    ! ii:   flag to indicate the diabatic states
    !       = 1, U11
    !       = 2, U22
    !       = 3, U12
    !***********************************************************************
    implicit none
    integer            :: nrap, ntap, ra1, ra2
    integer, parameter :: nqs = 5
    integer, parameter :: nqb = 5
    integer            :: igrad, ii
    integer            :: ntot, nqtc
    double precision   :: r, ut, dutdr, ra
    ! Local variables
    integer          :: i, j, k
    double precision :: u0
    double precision :: vs, vb, vbs
    double precision :: qs(nqs)
    double precision :: qb(nqb)
    double precision :: sqb(nqb), cqb(nqb)
    double precision :: qc(20), dut(20)
    double precision :: v0(nrap,ntap)
    double precision :: pi
    double precision :: dnc1, dnc2
    double precision :: dch1, dch2
    double precision :: knhb1(nrap,ntap), knhb2(nrap,ntap)
    double precision :: knc1(nrap,ntap), knc2(nrap,ntap)
    double precision :: kch1(nrap,ntap), kch2(nrap,ntap)
    double precision :: kch3(nrap,ntap)
    double precision :: kcnha1(nrap,ntap), kcnha2(nrap,ntap)
    double precision :: kcnhb1(nrap,ntap), kcnhb2(nrap,ntap)
    double precision :: knch1(nrap,ntap), knch2(nrap,ntap)
    !
    double precision :: jnhb1(nrap,ntap), jnhb2(nrap,ntap)
    double precision :: jnhb3(nrap,ntap)
    double precision :: jnc1(nrap,ntap), jnc2(nrap,ntap)
    double precision :: jnc3(nrap,ntap)
    double precision :: jcnha1(nrap,ntap), jcnha2(nrap,ntap)
    double precision :: jcnha3(nrap,ntap)
    double precision :: jcnhb1(nrap,ntap), jcnhb2(nrap,ntap)
    double precision :: jcnhb3(nrap,ntap)
    !
    double precision :: gnhb1(nrap,ntap), gnhb2(nrap,ntap)
    double precision :: gnhb3(nrap,ntap)
    !
    double precision :: gnc1(nrap,ntap), gnc2(nrap,ntap)
    double precision :: gnc3(nrap,ntap)
    !
    double precision :: v3, v4
    double precision :: torh1, torh2, rh1hx, rh2hx
    double precision :: ang1, ang2, ang3
    double precision :: len1, len2, len3
    double precision :: dd1, dd2, rr, ft2, fa2
    double precision :: ft1(nrap,ntap)
    double precision :: fa1(nrap,ntap)
    pi  = 4.d0*datan(1.d0)
    ut  = 0d0
    vs  = 0d0
    vb  = 0d0
    dut = 0d0
    do i=1, nqs
        qs(i) = qc(i)
    enddo
    do i=1, nqb
        qb(i)  = qc(i+nqs)
        sqb(i) = sin(qb(i))
        cqb(i) = cos(qb(i))
    enddo
    torh1 = qc(11)
    torh2 = qc(12)
    rh1hx = qc(13)
    rh2hx = qc(14)
    ang1  = qc(15)
    ang2  = qc(16)
    ang3  = qc(17)
    len1  = qc(18)
    len2  = qc(19)
    len3  = qc(20)
    call getpara(ii, nrap, ntap, &
        & knhb1, knhb2, &
        & knc1, knc2, &
        & kch1, kch2, &
        & kch3, &
        & kcnha1, kcnha2, &
        & kcnhb1, kcnhb2, &
        & knch1, knch2, &
        & jnhb1, jnhb2, &
        & jnhb3, &
        & jnc1, jnc2, jnc3, jcnha1, jcnha2, &
        & jcnha3, jcnhb1, jcnhb2, jcnhb3, &
        & gnhb1, gnhb2, gnhb3, &
        & gnc1, gnc2, gnc3, &
        & ft1, fa1)
    call getv0(ii, nrap, ntap, v0)
    ! Calculate potential energy
    ! CONSTANT
    u0 = v0(ra1,ra2)
    if( (ii .eq. 1) .or. (ii .eq. 2) )then
        ! STRETCHES
        ! STRETCH 1 (NHB)
        !       if (ii .eq. 1) then
        vs = gnhb1(ra1,ra2) *(1d0-exp(gnhb2(ra1,ra2)*(qc(1)-gnhb3(ra1,ra2))))**2d0
        vs = vs+knhb1(ra1,ra2)*(qc(1)-knhb2(ra1,ra2))**4d0
        ! STRETCH 2 (NC)
        vs = vs+gnc1(ra1,ra2) *(1d0-exp(gnc2(ra1,ra2)*(qc(2)-gnc3(ra1,ra2))))**2d0
        vs = vs+knc1(ra1,ra2)*(qc(2)-knc2(ra1,ra2))**4d0
        ! STRETCHES 3-5 (CH)
        do i=3, nqs
            vs = vs+kch1(ra1,ra2) *(1d0-exp(kch2(ra1,ra2)*(qc(i)-kch3(ra1,ra2))))**2d0
        enddo
        ! BENDS
        ! BEND 1 (CNHA)
        vb = kcnha1(ra1,ra2)*(cqb(1)-kcnha2(ra1,ra2))**2d0
        ! BEND 2 (CNHB)
        vb = vb+kcnhb1(ra1,ra2)*(cqb(2)-kcnhb2(ra1,ra2))**2d0
        ! BEND 3-5 (NCH)
        do i=3, nqb
            vb = vb+knch1(ra1,ra2)*(qb(i)-knch2(ra1,ra2))**2d0
        enddo
        ft2 = -3.0d0
        dd1 = 2.2d0
        v3  = ft1(ra1,ra2)*(1d0-(4d0*(torh2**3d0)-3d0*torh2)) * exp(((rh2hx-dd1)**2d0)/ft2)
        dd2 = 1.881d0
        rr  = 1.761d0
        fa2 = -2.0d0
        v4=(fa1(ra1,ra2)*(ang1-dd2)**2d0) *exp(((len1-rr)**2d0)/fa2) +(fa1(ra1,ra2)*(ang2-dd2)**2d0) *exp(((len2- &
            & rr)** 2d0)/fa2) +(fa1(ra1,ra2)*(ang3-dd2)**2d0) *exp(((len3-rr)**2d0)/fa2)
    else if( ii .eq. 3 )then
        ! STRETCHES
        ! STRETCH 1 (NHB)
        vs = jnhb1(ra1,ra2)*exp(jnhb2(ra1,ra2) *(qc(1)-jnhb3(ra1,ra2))**2d0)
        ! STRETCH 2 (NC)
        vs = vs+jnc1(ra1,ra2)*exp(jnc2(ra1,ra2) *(qc(2)-jnc3(ra1,ra2))**2d0)
        ! STRETCHES 3-5 (CH)
        ! BENDS
        ! BEND 1 (CNHA)
        vb = jcnha1(ra1,ra2)*exp(jcnha2(ra1,ra2) *(cqb(1)-jcnha3(ra1,ra2))**2d0)
        ! BEND 2 (CNHB)
        vb = vb+jcnhb1(ra1,ra2)*exp(jcnhb2(ra1,ra2) *(cqb(2)-jcnhb3(ra1,ra2))**2d0)
        ! BEND 3-5 (NCH)
    endif
    ! COMBINE
    ut = vs+vb+u0+v3+v4
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! GRADIENT
    if( igrad .eq. 1 )then
        if( (ii .eq. 1) .or. (ii .eq. 2) )then
            ! STRETCHES
            ! STRETCH 1 (NHB)
            !       if (ii .eq. 1) then
            dut(1)=-2d0*gnhb2(ra1,ra2)*gnhb1(ra1,ra2)* exp(gnhb2(ra1,ra2)*(qc(1)-gnhb3(ra1,ra2)))* (1d0- &
                & exp(gnhb2(ra1,ra2)*(qc(1)-gnhb3(ra1,ra2))))
            dut(1) = dut(1)+4d0*knhb1(ra1,ra2)*(qc(1)-knhb2(ra1,ra2))**3d0
            ! STRETCH 2 (NC)
            dut(2)=-2d0*gnc2(ra1,ra2)*gnc1(ra1,ra2)* exp(gnc2(ra1,ra2)*(qc(2)-gnc3(ra1,ra2)))* (1d0-exp(gnc2(ra1, &
                & ra2)* (qc(2)-gnc3(ra1,ra2))))
            dut(2) = dut(2)+4d0*knc1(ra1,ra2)*(qc(2)-knc2(ra1,ra2))**3d0
            ! STRETCHES 3-5 (CH)
            do i=3, nqs
                dut(i)=-2d0*kch2(ra1,ra2)*kch1(ra1,ra2)* exp(kch2(ra1,ra2)*(qs(i)-kch3(ra1,ra2)))* (1d0- exp(kch2(ra1, &
                    & ra2)*(qs(i)-kch3(ra1,ra2))))
            enddo
            ! BENDS
            ! BEND 1 (CNHA)
            dut(6) = kcnha1(ra1,ra2)*(cqb(1)-kcnha2(ra1,ra2)) *2d0*(-sqb(1))
            ! BEND 2 (CNHB)
            dut(7) = kcnhb1(ra1,ra2)*(cqb(2)-kcnhb2(ra1,ra2)) *2d0*(-sqb(2))
            ! BEND 3-5 (NCH)
            do i=3, nqb
                dut(nqs+i) = knch1(ra1,ra2)*(qb(i)-knch2(ra1,ra2))*2d0
            enddo
            dut(11) = 2d0*((rh1hx-dd1)/ft2)*exp(((rh1hx-dd1)**2d0)/ft2) *ft1(ra1,ra2)*(1d0-(4d0*(torh1**3d0)-3d0* &
                & torh1))
            dut(12) = 2d0*((rh2hx-dd1)/ft2)*exp(((rh2hx-dd1)**2d0)/ft2) *ft1(ra1,ra2)*(1d0-(4d0*(torh2**3d0)-3d0* &
                & torh2))
            dut(13) = ft1(ra1,ra2)*(-12d0*(torh1**2d0)+3d0) * exp(((rh1hx-dd1)**2d0)/ft2)
            dut(14) = ft1(ra1,ra2)*(-12d0*(torh2**2d0)+3d0) * exp(((rh2hx-dd1)**2d0)/ft2)
            dut(15) = 2d0*fa1(ra1,ra2)*(ang1-dd2) *exp(((len1-rr)**2d0)/fa2)
            dut(16) = 2d0*fa1(ra1,ra2)*(ang2-dd2) *exp(((len2-rr)**2d0)/fa2)
            dut(17) = 2d0*fa1(ra1,ra2)*(ang3-dd2) *exp(((len3-rr)**2d0)/fa2)
            dut(18) = (fa1(ra1,ra2)*(ang1-dd2)**2d0) *2d0*((len1-rr)/fa2)*exp(((len1-rr)**2d0)/fa2)
            !
            dut(19) = (fa1(ra1,ra2)*(ang2-dd2)**2d0) *2d0*((len2-rr)/fa2)*exp(((len2-rr)**2d0)/fa2)
            !
            dut(20) = (fa1(ra1,ra2)*(ang3-dd2)**2d0) *2d0*((len3-rr)/fa2)*exp(((len3-rr)**2d0)/fa2)
            dut(11) = 0d0
            dut(13) = 0d0
        else if( ii .eq. 3 )then
            ! STRETCH 1 (NHB)
            dut(1)=jnhb2(ra1,ra2)*(qc(1)-jnhb3(ra1,ra2))*2d0* jnhb1(ra1,ra2)*exp(jnhb2(ra1,ra2) *(qc(1)-jnhb3(ra1, &
                & ra2))**2d0)
            ! STRETCH 2 (NC)
            dut(2) = jnc2(ra1,ra2)*(qc(2)-jnc3(ra1,ra2))*2d0* jnc1(ra1,ra2)*exp(jnc2(ra1,ra2) *(qc(2)-jnc3(ra1, &
                & ra2))** 2d0)
            ! BENDS
            ! BEND 1 (CNHA)
            dut(6)=jcnha2(ra1,ra2)*(cqb(1)-jcnha3(ra1,ra2))*2d0* jcnha1(ra1,ra2)*exp(jcnha2(ra1,ra2) *(cqb(1)- &
                & jcnha3(ra1,ra2))**2d0) *(-sqb(1))
            ! BEND 2 (CNHB)
            dut(7)=jcnhb2(ra1,ra2)*(cqb(2)-jcnhb3(ra1,ra2))*2d0 *jcnhb1(ra1,ra2)*exp(jcnhb2(ra1,ra2) *(cqb(2)- &
                & jcnhb3(ra1,ra2))**2d0) *(-sqb(2))
        endif
    endif
end subroutine solver

!> @brief Evaluate getV0 for the CH3NH2_APRP potential.
subroutine getv0( ii, nrap, ntap, v0 )
    ! subrountine gives constants in tertiary eqns
    implicit none
    integer          :: ii
    integer          :: nrap, ntap
    double precision :: v0(nrap,ntap)
    v0 = 0d0
    if( ii .eq. 1 )then
        v0(1,1) = -0.2734d0
        v0(1,2) = -0.0359d0
        v0(1,3) = -0.2495d0
        v0(2,1) = -0.2243d0
        v0(2,2) = -0.0108d0
        v0(2,3) = -0.3548d0
        !!!!!!!!!!!!!
        v0(3,1) = -0.1393d0
        v0(3,2) = -0.1971d0
        v0(3,3) = -0.4208d0
        !!!!!!!!!!!!
        v0(4,1) = -0.5648d0
        v0(4,2) = -0.5648d0
        v0(4,3) = -0.5648d0
    else if( ii .eq. 2 )then
        v0(1,1) = -0.1853d0
        v0(1,2) = -0.1269d0
        v0(1,3) = -0.3474d0
        v0(2,1) = -0.1265d0
        v0(2,2) = -0.1550d0
        v0(2,3) = -0.2869d0
        !!!!!!!!!!!
        v0(3,1) = -0.1910d0
        v0(3,2) = -0.0865d0
        v0(3,3) = -0.0470d0
        !!!!!!!!!!!
        v0(4,1) = -0.0481d0
        v0(4,2) = -0.0481d0
        v0(4,3) = -0.0481d0
    else if( ii .eq. 3 )then
        v0(1,1) = -0.4000d0
        v0(1,2) = -0.2125d0
        v0(1,3) = -0.0093d0
        v0(2,1) = -0.1812d0
        v0(2,2) = -0.0606d0
        v0(2,3) = -0.0090d0
        !!!!!!!!!!!!
        v0(3,1) = -5.0661d0
        v0(3,2) = -3.5468d0
        v0(3,3) = -0.0369d0
        !!!!!!!!!!!
        v0(4,1) = -0.0000d0
        v0(4,2) = -0.0000d0
        v0(4,3) = -0.0000d0
    endif
end subroutine getv0

!> @brief Evaluate getpara for the CH3NH2_APRP potential.
subroutine getpara( ii, nrap, ntap, &
        & knhb1, knhb2, &
        & knc1, knc2, &
        & kch1, kch2, &
        & kch3, &
        & kcnha1, kcnha2, &
        & kcnhb1, kcnhb2, &
        & knch1, knch2, &
        & jnhb1, jnhb2, &
        & jnhb3, &
        & jnc1, jnc2, jnc3, jcnha1, jcnha2, &
        & jcnha3, jcnhb1, jcnhb2, jcnhb3, &
        & gnhb1, gnhb2, gnhb3, &
        & gnc1, gnc2, gnc3, &
        & ft1, fa1 )
    ! subroutine gives parameters for the tertiary contribution
    implicit none
    integer          :: ii
    integer          :: nrap, i, j, ntap
    double precision :: knhb1(nrap,ntap), knhb2(nrap,ntap)
    double precision :: knc1(nrap,ntap), knc2(nrap,ntap)
    double precision :: kch1(nrap,ntap), kch2(nrap,ntap)
    double precision :: kch3(nrap,ntap)
    double precision :: kcnha1(nrap,ntap), kcnha2(nrap,ntap)
    double precision :: kcnhb1(nrap,ntap), kcnhb2(nrap,ntap)
    double precision :: knch1(nrap,ntap), knch2(nrap,ntap)
    !
    double precision :: jnhb1(nrap,ntap), jnhb2(nrap,ntap)
    double precision :: jnhb3(nrap,ntap)
    double precision :: jnc1(nrap,ntap), jnc2(nrap,ntap)
    double precision :: jnc3(nrap,ntap)
    double precision :: jcnha1(nrap,ntap), jcnha2(nrap,ntap)
    double precision :: jcnha3(nrap,ntap)
    double precision :: jcnhb1(nrap,ntap), jcnhb2(nrap,ntap)
    double precision :: jcnhb3(nrap,ntap)
    !
    double precision :: gnhb1(nrap,ntap), gnhb2(nrap,ntap)
    double precision :: gnhb3(nrap,ntap)
    !
    double precision :: gnc1(nrap,ntap), gnc2(nrap,ntap)
    double precision :: gnc3(nrap,ntap)
    double precision :: ft1(nrap,ntap), fa1(nrap,ntap)
    knhb1  = 0d0
    knhb2  = 0d0
    knc1   = 0d0
    knc2   = 0d0
    kch1   = 0d0
    kch2   = 0d0
    kch3   = 0d0
    kcnha1 = 0d0
    kcnha2 = 0d0
    kcnhb1 = 0d0
    kcnhb2 = 0d0
    knch1  = 0d0
    knch2  = 0d0
    jnhb1  = 0d0
    jnhb2  = 0d0
    jnhb3  = 0d0
    jnc1   = 0d0
    jnc2   = 0d0
    jnc3   = 0d0
    jcnha1 = 0d0
    jcnha2 = 0d0
    jcnha3 = 0d0
    jcnhb1 = 0d0
    jcnhb2 = 0d0
    jcnhb3 = 0d0
    if( ii .eq. 1 )then
        ! STRETCH 1 (NHB)
        gnhb1(1,2) = 6.420621d0
        gnhb1(2,2) = 6.420621d0
        gnhb1(3,2) = 6.420621d0
        gnhb1(4,2) = 6.420621d0
        gnhb1(1,3) = 6.420621d0
        gnhb1(2,3) = 6.420621d0
        gnhb1(3,3) = 6.420621d0
        gnhb1(4,3) = gnhb1(4,2)
        gnhb1(1,1) = 4.3d0
        gnhb1(2,1) = 4.4d0
        gnhb1(3,1) = 8d0
        gnhb1(4,1) = gnhb1(4,2)
        gnhb2(1,2) = -1.85d0
        gnhb2(2,2) = -1.85d0
        gnhb2(3,2) = -1.85d0
        gnhb2(4,2) = -1.96d0
        gnhb2(1,3) = -1.85d0
        gnhb2(2,3) = -1.85d0
        gnhb2(3,3) = -1.85d0
        gnhb2(4,3) = gnhb2(4,2)
        gnhb2(1,1) = -1.85d0
        gnhb2(2,1) = -1.85d0
        gnhb2(3,1) = -2.0d0
        gnhb2(4,1) = gnhb2(4,2)
        gnhb3(1,2) = 1.01574847d0
        gnhb3(2,2) = 1.01574847d0
        gnhb3(3,2) = 1.02d0
        gnhb3(4,2) = 1.00368d0
        gnhb3(1,3) = 1.01574847d0
        gnhb3(2,3) = 1.01574847d0
        gnhb3(3,3) = 1.02d0
        gnhb3(4,3) = gnhb3(4,2)
        gnhb3(1,1) = 1.08d0
        gnhb3(2,1) = 1.07d0
        gnhb3(3,1) = 1.0d0
        gnhb3(4,1) = gnhb3(4,2)
        knhb1(1,2) = 1.0d0
        knhb1(2,2) = 1.0d0
        knhb1(3,2) = 0.5d0
        knhb1(4,2) = 0.5d0
        knhb1(1,3) = 6.0d0
        knhb1(2,3) = 3.0d0
        knhb1(3,3) = 1.0d0
        knhb1(4,3) = knhb1(4,2)
        knhb1(1,1) = 1.0d0
        knhb1(2,1) = 1.0d0
        knhb1(3,1) = 0.5d0
        knhb1(4,1) = knhb1(4,2)
        knhb2(1,2) = 1.01574847d0
        knhb2(2,2) = 1.01574847d0
        knhb2(3,2) = 1.02d0
        knhb2(4,2) = 1.00368d0
        knhb2(1,3) = 1.01574847d0
        knhb2(2,3) = 1.01574847d0
        knhb2(3,3) = 1.02d0
        knhb2(4,3) = knhb2(4,2)
        knhb2(1,1) = 1.08d0
        knhb2(2,1) = 1.07d0
        knhb2(3,1) = 1.0d0
        knhb2(4,1) = knhb2(4,2)
        ! STRETCH 2 (NC)
        gnc1(1,2) = 4.584032444d0
        gnc1(2,2) = 4.543215344d0
        gnc1(3,2) = 4.0d0
        gnc1(4,2) = 2.7d0
        gnc1(1,3) = 4.43273706d0
        gnc1(2,3) = 4.412600624d0
        gnc1(3,3) = 3.1d0
        gnc1(4,3) = gnc1(4,2)
        gnc1(1,1) = 4.90621542d0
        gnc1(2,1) = 4.870024258d0
        gnc1(3,1) = 4.0d0
        gnc1(4,1) = gnc1(4,2)
        gnc2(1,2) = -2.0d0
        gnc2(2,2) = -2.0d0
        gnc2(3,2) = -1.9d0
        gnc2(4,2) = -2.6d0
        gnc2(1,3) = -2.0d0
        gnc2(2,3) = -2.0d0
        gnc2(3,3) = -2.3d0
        gnc2(4,3) = gnc2(4,2)
        gnc2(1,1) = -1.9d0
        gnc2(2,1) = -1.9d0
        gnc2(3,1) = -2.3d0
        gnc2(4,1) = gnc2(4,2)
        gnc3(1,2) = 1.4597436d0
        gnc3(2,2) = 1.4597436d0
        gnc3(3,2) = 1.4597436d0
        gnc3(4,2) = 1.40443d0
        gnc3(1,3) = 1.4597436d0
        gnc3(2,3) = 1.4597436d0
        gnc3(3,3) = 1.44d0
        gnc3(4,3) = gnc3(4,2)
        gnc3(1,1) = 1.4597436d0
        gnc3(2,1) = 1.4597436d0
        gnc3(3,1) = 1.45d0
        gnc3(4,1) = gnc3(4,2)
        knc1(1,2) = 3.0d0
        knc1(2,2) = 5.0d0
        knc1(3,2) = 8.0d0
        knc1(4,2) = 1.7d0
        knc1(1,3) = 2.0d0
        knc1(2,3) = 3.0d0
        knc1(3,3) = 4.5d0
        knc1(4,3) = knc1(4,2)
        knc1(1,1) = 4.0d0
        knc1(2,1) = 5.0d0
        knc1(3,1) = 4.0d0
        knc1(4,1) = knc1(4,2)
        knc2(1,2) = 1.4597436d0
        knc2(2,2) = 1.4597436d0
        knc2(3,2) = 1.4597436d0
        knc2(4,2) = 1.40443d0
        knc2(1,3) = 1.4597436d0
        knc2(2,3) = 1.4597436d0
        knc2(3,3) = 1.44d0
        knc2(4,3) = knc2(4,2)
        knc2(1,1) = 1.4597436d0
        knc2(2,1) = 1.4597436d0
        knc2(3,1) = 1.45d0
        knc2(4,1) = knc2(4,2)
        ! STRETCHES 3-5 (CH)
        kch1(1,2) = 4.680905028d0
        kch1(2,2) = 4.660046401d0
        kch1(3,2) = 3.103188056d0
        kch1(4,2) = 2.783998334d0
        kch1(1,3) = 4.778593954d0
        kch1(2,3) = 4.763899798d0
        kch1(3,3) = 3.552720384d0
        kch1(4,3) = kch1(4,2)
        kch1(1,1) = 4.655054198d0
        kch1(2,1) = 4.633285078d0
        kch1(3,1) = 2.861006596d0
        kch1(4,1) = kch1(4,2)
        kch2(1,2) = -2.0d0
        kch2(2,2) = -1.9d0
        kch2(3,2) = -2.3d0
        kch2(4,2) = -2.3d0
        kch2(1,3) = -2.0d0
        kch2(2,3) = -2.0d0
        kch2(3,3) = -2.2d0
        kch2(4,3) = kch2(4,2)
        kch2(1,1) = -2.0d0
        kch2(2,1) = -2.0d0
        kch2(3,1) = -2.4d0
        kch2(4,1) = kch2(4,2)
        kch3(1,2) = 1.0894707d0
        kch3(2,2) = 1.0894707d0
        kch3(3,2) = 1.0894707d0
        kch3(4,2) = 1.094d0
        kch3(1,3) = 1.0894707d0
        kch3(2,3) = 1.0894707d0
        kch3(3,3) = 1.0894707d0
        kch3(4,3) = kch3(4,2)
        kch3(1,1) = 1.0894707d0
        kch3(2,1) = 1.0894707d0
        kch3(3,1) = 1.0894707d0
        kch3(4,1) = kch3(4,2)
        ! BEND 1 (CNHA)
        kcnha1(1,2) = 3.7d0
        kcnha1(2,2) = 3.4d0
        kcnha1(3,2) = 1.0d0
        kcnha1(4,2) = 0.0d0
        kcnha1(1,3) = 3.7d0
        kcnha1(2,3) = 3.8d0
        kcnha1(3,3) = 1.5d0
        kcnha1(4,3) = kcnha1(4,2)
        kcnha1(1,1) = 3.0d0
        kcnha1(2,1) = 3.3d0
        kcnha1(3,1) = 2.0d0
        kcnha1(4,1) = kcnha1(4,2)
        kcnha2(1,2) = -0.354096319d0
        kcnha2(2,2) = -0.354096319d0
        kcnha2(3,2) = -0.354096319d0
        kcnha2(4,2) = 0.0d0
        kcnha2(1,3) = -0.544639035d0
        kcnha2(2,3) = -0.559192903d0
        kcnha2(3,3) = -0.544639035d0
        kcnha2(4,3) = kcnha2(4,2)
        kcnha2(1,1) = -0.207911691d0
        kcnha2(2,1) = -0.207911691d0
        kcnha2(3,1) = -0.258819045d0
        kcnha2(4,1) = kcnha2(4,2)
        ! BEND 2 (CNHB)
        kcnhb1(1,2) = 3.9d0
        kcnhb1(2,2) = 3.4d0
        kcnhb1(3,2) = 3.3d0
        kcnhb1(4,2) = 2.6d0
        kcnhb1(1,3) = 3.7d0
        kcnhb1(2,3) = 3.3d0
        kcnhb1(3,3) = 4.0d0
        kcnhb1(4,3) = kcnhb1(4,2)
        kcnhb1(1,1) = 3.4d0
        kcnhb1(2,1) = 3.5d0
        kcnhb1(3,1) = 4.0d0
        kcnhb1(4,1) = kcnhb1(4,2)
        kcnhb2(1,2) = -0.325568154d0
        kcnhb2(2,2) = -0.355262388d0
        kcnhb2(3,2) = -0.573576436d0
        kcnhb2(4,2) = -0.788976538d0
        kcnhb2(1,3) = -0.5d0
        kcnhb2(2,3) = -0.573576436d0
        kcnhb2(3,3) = -0.64278761d0
        kcnhb2(4,3) = kcnhb2(4,2)
        kcnhb2(1,1) = -0.173648178d0
        kcnhb2(2,1) = -0.207911691d0
        kcnhb2(3,1) = -0.5d0
        kcnhb2(4,1) = kcnhb2(4,2)
        ! BEND 3-5 (NCH)
        knch1(1,2) = 3.7d0
        knch1(2,2) = 3.4d0
        knch1(3,2) = 3.7d0
        knch1(4,2) = 3.0d0
        knch1(1,3) = 3.5d0
        knch1(2,3) = 3.5d0
        knch1(3,3) = 3.6d0
        knch1(4,3) = knch1(4,2)
        knch1(1,1) = 3.5d0
        knch1(2,1) = 3.6d0
        knch1(3,1) = 3.5d0
        knch1(4,1) = knch1(4,2)
        knch2(1,2) = 1.92d0
        knch2(2,2) = 1.92d0
        knch2(3,2) = 1.97d0
        knch2(4,2) = 1.93d0
        knch2(1,3) = 1.97d0
        knch2(2,3) = 1.97d0
        knch2(3,3) = 1.95d0
        knch2(4,3) = knch2(4,2)
        knch2(1,1) = 1.95d0
        knch2(2,1) = 1.93d0
        knch2(3,1) = 1.93d0
        knch2(4,1) = knch2(4,2)
        ! Methylrot and Methylang
        ft1(1,2) = 0.0305d0
        ft1(2,2) = 0.0305d0
        ft1(3,2) = 0.0035d0
        ft1(4,2) = 0.0035d0
        ft1(1,3) = 0.0305d0
        ft1(2,3) = 0.0305d0
        ft1(3,3) = ft1(3,2)
        ft1(4,3) = ft1(4,2)
        ft1(1,1) = 0.0305d0
        ft1(2,1) = 0.0305d0
        ft1(3,1) = ft1(3,2)
        ft1(4,1) = ft1(4,2)
        fa1(1,2) = 0.5d0
        fa1(2,2) = 0.5d0
        fa1(3,2) = 1.5d0
        fa1(4,2) = 1.5d0
        fa1(1,3) = 0.5d0
        fa1(2,3) = 0.5d0
        fa1(3,3) = fa1(3,2)
        fa1(4,3) = fa1(4,2)
        fa1(1,1) = 0.5d0
        fa1(2,1) = 0.5d0
        fa1(3,1) = fa1(3,2)
        fa1(4,1) = fa1(4,2)
    else if( ii .eq. 2 )then
        ! STRETCH 1 (NHB)
        gnhb1(1,2) = 1.0d0
        gnhb1(2,2) = 1.0d0
        gnhb1(3,2) = 1.5d0
        gnhb1(4,2) = 3.2d0
        gnhb1(1,3) = 1.0d0
        gnhb1(2,3) = 1.0d0
        gnhb1(3,3) = 1.8d0
        gnhb1(4,3) = gnhb1(4,2)
        gnhb1(1,1) = 1.0d0
        gnhb1(2,1) = 1.0d0
        gnhb1(3,1) = 0.8d0
        gnhb1(4,1) = gnhb1(4,2)
        gnhb2(1,2) = -3.8d0
        gnhb2(2,2) = -3.8d0
        gnhb2(3,2) = -3.8d0
        gnhb2(4,2) = -2.5d0
        gnhb2(1,3) = -3.5d0
        gnhb2(2,3) = -4.2d0
        gnhb2(3,3) = -3.6d0
        gnhb2(4,3) = gnhb2(4,2)
        gnhb2(1,1) = -3.1d0
        gnhb2(2,1) = -3.3d0
        gnhb2(3,1) = -3.3d0
        gnhb2(4,1) = gnhb2(4,2)
        gnhb3(1,2) = 1.040695848d0
        gnhb3(2,2) = 1.040695848d0
        gnhb3(3,2) = 1.0d0
        gnhb3(4,2) = 1.030034236d0
        gnhb3(1,3) = 1.040695848d0
        gnhb3(2,3) = 1.055d0
        gnhb3(3,3) = 1.01d0
        gnhb3(4,3) = gnhb3(4,2)
        gnhb3(1,1) = 1.11d0
        gnhb3(2,1) = 1.08d0
        gnhb3(3,1) = 1.08d0
        gnhb3(4,1) = gnhb3(4,2)
        knhb1(1,2) = 2.0d0
        knhb1(2,2) = 2.0d0
        knhb1(3,2) = 3.0d0
        knhb1(4,2) = 11.0d0
        knhb1(1,3) = 2.0d0
        knhb1(2,3) = 2.0d0
        knhb1(3,3) = 6.0d0
        knhb1(4,3) = knhb1(4,2)
        knhb1(1,1) = 2.0d0
        knhb1(2,1) = 2.0d0
        knhb1(3,1) = 3.0d0
        knhb1(4,1) = knhb1(4,2)
        knhb2(1,2) = 1.040695848d0
        knhb2(2,2) = 1.040695848d0
        knhb2(3,2) = 1.0d0
        knhb2(4,2) = 1.030034236d0
        knhb2(1,3) = 1.040695848d0
        knhb2(2,3) = 1.055d0
        knhb2(3,3) = 1.01d0
        knhb2(4,3) = knhb2(4,2)
        knhb2(1,1) = 1.11d0
        knhb2(2,1) = 1.08d0
        knhb2(3,1) = 1.08d0
        knhb2(4,1) = knhb2(4,2)
        ! STRETCH 2 (NC)
        gnc1(1,2) = 3.0d0
        gnc1(2,2) = 3.0d0
        gnc1(3,2) = 3.0d0
        gnc1(4,2) = 3.34292049d0
        gnc1(1,3) = 3.0d0
        gnc1(2,3) = 3.0d0
        gnc1(3,3) = 3.8d0
        gnc1(4,3) = gnc1(4,2)
        gnc1(1,1) = 3.0d0
        gnc1(2,1) = 3.0d0
        gnc1(3,1) = 4.0d0
        gnc1(4,1) = gnc1(4,2)
        gnc2(1,2) = -2.3d0
        gnc2(2,2) = -2.3d0
        gnc2(3,2) = -2.3d0
        gnc2(4,2) = -2.1d0
        gnc2(1,3) = -2.3d0
        gnc2(2,3) = -2.4d0
        gnc2(3,3) = -2.0d0
        gnc2(4,3) = gnc2(4,2)
        gnc2(1,1) = -2.3d0
        gnc2(2,1) = -2.3d0
        gnc2(3,1) = -1.9d0
        gnc2(4,1) = gnc2(4,2)
        gnc3(1,2) = 1.43081d0
        gnc3(2,2) = 1.43081d0
        gnc3(3,2) = 1.43081d0
        gnc3(4,2) = 1.44623d0
        gnc3(1,3) = 1.43081d0
        gnc3(2,3) = 1.43081d0
        gnc3(3,3) = 1.445d0
        gnc3(4,3) = gnc3(4,2)
        gnc3(1,1) = 1.43081d0
        gnc3(2,1) = 1.43081d0
        gnc3(3,1) = 1.43081d0
        gnc3(4,1) = gnc3(4,2)
        knc1(1,2) = 2.0d0
        knc1(2,2) = 2.0d0
        knc1(3,2) = 2.0d0
        knc1(4,2) = 2.0d0
        knc1(1,3) = 2.0d0
        knc1(2,3) = 2.0d0
        knc1(3,3) = 2.5d0
        knc1(4,3) = knc1(4,2)
        knc1(1,1) = 2.0d0
        knc1(2,1) = 2.0d0
        knc1(3,1) = 1.0d0
        knc1(4,1) = knc1(4,2)
        knc2(1,2) = 1.43081d0
        knc2(2,2) = 1.43081d0
        knc2(3,2) = 1.43081d0
        knc2(4,2) = 1.44623d0
        knc2(1,3) = 1.43081d0
        knc2(2,3) = 1.43081d0
        knc2(3,3) = 1.445d0
        knc2(4,3) = knc2(4,2)
        knc2(1,1) = 1.43081d0
        knc2(2,1) = 1.43081d0
        knc2(3,1) = 1.43081d0
        knc2(4,1) = knc2(4,2)
        ! STRETCHES 3-5 (CH)
        kch1(1,2) = 2.811209734d0
        kch1(2,2) = 3.037064354d0
        kch1(3,2) = 4.641992726d0
        kch1(4,2) = 4.605257336d0
        kch1(1,3) = 2.164122642d0
        kch1(2,3) = 2.476781628d0
        kch1(3,3) = 4.53750095d0
        kch1(4,3) = kch1(4,2)
        kch1(1,1) = 3.084140076d0
        kch1(2,1) = 3.237068144d0
        kch1(3,1) = 4.646074436d0
        kch1(4,1) = kch1(4,2)
        kch2(1,2) = -2.4d0
        kch2(2,2) = -2.4d0
        kch2(3,2) = -2.0d0
        kch2(4,2) = -1.9d0
        kch2(1,3) = -2.6d0
        kch2(2,3) = -2.5d0
        kch2(3,3) = -2.0d0
        kch2(4,3) = kch2(4,2)
        kch2(1,1) = -2.3d0
        kch2(2,1) = -2.3d0
        kch2(3,1) = -2.0d0
        kch2(4,1) = kch2(4,2)
        kch3(1,2) = 1.0894707d0
        kch3(2,2) = 1.0894707d0
        kch3(3,2) = 1.0894707d0
        kch3(4,2) = 1.0894707d0
        kch3(1,3) = 1.0894707d0
        kch3(2,3) = 1.0894707d0
        kch3(3,3) = 1.0894707d0
        kch3(4,3) = kch3(4,2)
        kch3(1,1) = 1.0894707d0
        kch3(2,1) = 1.0894707d0
        kch3(3,1) = 1.0894707d0
        kch3(4,1) = kch3(4,2)
        ! BEND 1 (CNHA)
        kcnha1(1,2) = 3.2d0
        kcnha1(2,2) = 3.0d0
        kcnha1(3,2) = 1.0d0
        kcnha1(4,2) = 0.0d0
        kcnha1(1,3) = 3.6d0
        kcnha1(2,3) = 4.0d0
        kcnha1(3,3) = 0.5d0
        kcnha1(4,3) = kcnha1(4,2)
        kcnha1(1,1) = 1.9d0
        kcnha1(2,1) = 2.1d0
        kcnha1(3,1) = 1.0d0
        kcnha1(4,1) = kcnha1(4,2)
        kcnha2(1,2) = -0.4539905d0
        kcnha2(2,2) = -0.48480962d0
        kcnha2(3,2) = -0.48480962d0
        kcnha2(4,2) = 0d0
        kcnha2(1,3) = -0.559192903d0
        kcnha2(2,3) = -0.505340076d0
        kcnha2(3,3) = -0.5d0
        kcnha2(4,3) = kcnha2(4,2)
        kcnha2(1,1) = -0.438371147d0
        kcnha2(2,1) = -0.4539905d0
        kcnha2(3,1) = -0.64278761d0
        kcnha2(4,1) = kcnha2(4,2)
        ! BEND 2 (CNHB)
        kcnhb1(1,2) = 3.7d0
        kcnhb1(2,2) = 3.3d0
        kcnhb1(3,2) = 2.8d0
        kcnhb1(4,2) = 3.3d0
        kcnhb1(1,3) = 3.6d0
        kcnhb1(2,3) = 3.3d0
        kcnhb1(3,3) = 2.5d0
        kcnhb1(4,3) = kcnhb1(4,2)
        kcnhb1(1,1) = 2.3d0
        kcnhb1(2,1) = 2.1d0
        kcnhb1(3,1) = 3.0d0
        kcnhb1(4,1) = kcnhb1(4,2)
        kcnhb2(1,2) = -0.406736643d0
        kcnhb2(2,2) = -0.422618262d0
        kcnhb2(3,2) = -0.406736643d0
        kcnhb2(4,2) = -0.286195791d0
        kcnhb2(1,3) = -0.529919264d0
        kcnhb2(2,3) = -0.528778795d0
        kcnhb2(3,3) = -0.342020143d0
        kcnhb2(4,3) = kcnhb2(4,2)
        kcnhb2(1,1) = -0.342020143d0
        kcnhb2(2,1) = -0.374606593d0
        kcnhb2(3,1) = -0.422618262d0
        kcnhb2(4,1) = kcnhb2(4,2)
        ! BEND 3-5 (NCH)
        knch1(1,2) = 3.5d0
        knch1(2,2) = 3.7d0
        knch1(3,2) = 3.0d0
        knch1(4,2) = 2.9d0
        knch1(1,3) = 3.0d0
        knch1(2,3) = 2.8d0
        knch1(3,3) = 2.6d0
        knch1(4,3) = knch1(4,2)
        knch1(1,1) = 3.6d0
        knch1(2,1) = 3.6d0
        knch1(3,1) = 3.1d0
        knch1(4,1) = knch1(4,2)
        knch2(1,2) = 1.88d0
        knch2(2,2) = 1.88d0
        knch2(3,2) = 1.89d0
        knch2(4,2) = 1.92d0
        knch2(1,3) = 1.87d0
        knch2(2,3) = 1.88d0
        knch2(3,3) = 1.9d0
        knch2(4,3) = knch2(4,2)
        knch2(1,1) = 1.90d0
        knch2(2,1) = 1.91d0
        knch2(3,1) = 1.9d0
        knch2(4,1) = knch2(4,2)
        ! Methylrot and Methylang
        ft1(1,2) = 0.012d0
        ft1(2,2) = 0.012d0
        ft1(3,2) = 0.012d0
        ft1(4,2) = 0.012d0
        ft1(1,3) = 0.012d0
        ft1(2,3) = 0.012d0
        ft1(3,3) = ft1(3,2)
        ft1(4,3) = ft1(4,2)
        ft1(1,1) = 0.012d0
        ft1(2,1) = 0.012d0
        ft1(3,1) = ft1(3,2)
        ft1(4,1) = ft1(4,2)
        fa1(1,2) = 0.5d0
        fa1(2,2) = 0.5d0
        fa1(3,2) = 0.8d0
        fa1(4,2) = 0.8d0
        fa1(1,3) = 0.5d0
        fa1(2,3) = 0.5d0
        fa1(3,3) = fa1(3,2)
        fa1(4,3) = fa1(4,2)
        fa1(1,1) = 0.5d0
        fa1(2,1) = 0.5d0
        fa1(3,1) = fa1(3,2)
        fa1(4,1) = fa1(4,2)
    else if( ii .eq. 3 )then
        ! STRETCH 1 (NHB)
        jnhb1(1,2) = 1.4d0
        jnhb1(2,2) = 1.4d0
        jnhb1(3,2) = 0.9d0
        jnhb1(4,2) = 0d0
        jnhb1(1,3) = 0.28d0
        jnhb1(2,3) = 0.28d0
        jnhb1(3,3) = 1.15d0
        jnhb1(4,3) = jnhb1(4,2)
        jnhb1(1,1) = 1.4d0
        jnhb1(2,1) = 1.44d0
        jnhb1(3,1) = 1.28d0
        jnhb1(4,1) = jnhb1(4,2)
        jnhb2(1,2) = -8.5d0
        jnhb2(2,2) = -10.0d0
        jnhb2(3,2) = -0.4d0
        jnhb2(4,2) = -1d0
        jnhb2(1,3) = -22d0
        jnhb2(2,3) = -22d0
        jnhb2(3,3) = -35d0
        jnhb2(4,3) = jnhb2(4,2)
        jnhb2(1,1) = -6.5d0
        jnhb2(2,1) = -6.5d0
        jnhb2(3,1) = -0.5d0
        jnhb2(4,1) = jnhb2(4,2)
        jnhb3(1,2) = 1.8145815d0
        jnhb3(2,2) = 1.8145815d0
        jnhb3(3,2) = 1.2d0
        jnhb3(4,2) = 1d0
        jnhb3(1,3) = 1.7145815d0
        jnhb3(2,3) = 1.7145815d0
        jnhb3(3,3) = 2d0
        jnhb3(4,3) = jnhb3(4,2)
        jnhb3(1,1) = 2.0145815d0
        jnhb3(2,1) = 2.0145815d0
        jnhb3(3,1) = 1.1d0
        jnhb3(4,1) = jnhb3(4,2)
        ! STRETCH 2 (NC)
        jnc1(1,2) = 1.65d0
        jnc1(2,2) = 1.6d0
        jnc1(3,2) = 0.95
        jnc1(4,2) = 0d0
        jnc1(1,3) = 0.32d0
        jnc1(2,3) = 0.4d0
        jnc1(3,3) = 0.34d0
        jnc1(4,3) = jnc1(4,2)
        jnc1(1,1) = 2.4d0
        jnc1(2,1) = 2.5d0
        jnc1(3,1) = 1.3d0
        jnc1(4,1) = jnc1(4,2)
        jnc2(1,2) = -25.0d0
        jnc2(2,2) = -25.0d0
        jnc2(3,2) = -0.7d0
        jnc2(4,2) = -1
        jnc2(1,3) = -30.0d0
        jnc2(2,3) = -30.0d0
        jnc2(3,3) = -30.0d0
        jnc2(4,3) = jnc2(4,2)
        jnc2(1,1) = -30.0d0
        jnc2(2,1) = -30.0d0
        jnc2(3,1) = -1.2d0
        jnc2(4,1) = jnc2(4,2)
        jnc3(1,2) = 1.9597436d0
        jnc3(2,2) = 2.0597436d0
        jnc3(3,2) = 1.7d0
        jnc3(4,2) = 1d0
        jnc3(1,3) = 2.3597436d0
        jnc3(2,3) = 2.3597436d0
        jnc3(3,3) = 2.6d0
        jnc3(4,3) = jnc3(4,2)
        jnc3(1,1) = 1.9597436d0
        jnc3(2,1) = 1.9597436d0
        jnc3(3,1) = 1.55d0
        jnc3(4,1) = jnc3(4,2)
        ! BEND 1 (CNHA)
        jcnha1(1,2) = 0.1d0
        jcnha1(2,2) = 0.06d0
        jcnha1(3,2) = 0.93d0
        jcnha1(4,2) = 0d0
        jcnha1(1,3) = 0.007d0
        jcnha1(2,3) = 0.005d0
        jcnha1(3,3) = 0.048d0
        jcnha1(4,3) = jcnha1(4,2)
        jcnha1(1,1) = 0.2d0
        jcnha1(2,1) = 0.18d0
        jcnha1(3,1) = 1.3d0
        jcnha1(4,1) = jcnha1(4,2)
        jcnha2(1,2) = -2.5d0
        jcnha2(2,2) = -3.0d0
        jcnha2(3,2) = -1.5d0
        jcnha2(4,2) = -1d0
        jcnha2(1,3) = -1.2d0
        jcnha2(2,3) = -0.5d0
        jcnha2(3,3) = -3.0d0
        jcnha2(4,3) = jcnha2(4,2)
        jcnha2(1,1) = -0.7d0
        jcnha2(2,1) = -2.3d0
        jcnha2(3,1) = -2.0d0
        jcnha2(4,1) = jcnha2(4,2)
        jcnha3(1,2) = -0.354096319d0
        jcnha3(2,2) = 0.160948072d0
        jcnha3(3,2) = -0.173648178
        jcnha3(4,2) = 1
        jcnha3(1,3) = 0.160948072d0
        jcnha3(2,3) = -0.354096319d0
        jcnha3(3,3) = 0.329887214
        jcnha3(4,3) = jcnha3(4,2)
        jcnha3(1,1) = -0.186319468d0
        jcnha3(2,1) = 0.160948072d0
        jcnha3(3,1) = -0.258819045
        jcnha3(4,1) = jcnha3(4,2)
        ! BEND 2 (CNHB)
        jcnhb1(1,2) = 0.15d0
        jcnhb1(2,2) = 0.06d0
        jcnhb1(3,2) = 1.4d0
        jcnhb1(4,2) = 0
        jcnhb1(1,3) = 0.006d0
        jcnhb1(2,3) = 0.005d0
        jcnhb1(3,3) = 0.04d0
        jcnhb1(4,3) = jcnhb1(4,2)
        jcnhb1(1,1) = 0.32d0
        jcnhb1(2,1) = 0.14d0
        jcnhb1(3,1) = 1.9d0
        jcnhb1(4,1) = jcnhb1(4,2)
        jcnhb2(1,2) = -0.8d0
        jcnhb2(2,2) = -2.5d0
        jcnhb2(3,2) = -1d0
        jcnhb2(4,2) = -1
        jcnhb2(1,3) = -3.0d0
        jcnhb2(2,3) = -0.5d0
        jcnhb2(3,3) = -1.0d0
        jcnhb2(4,3) = jcnhb2(4,2)
        jcnhb2(1,1) = -1.0d0
        jcnhb2(2,1) = -1.2d0
        jcnhb2(3,1) = -0.9d0
        jcnhb2(4,1) = jcnhb2(4,2)
        jcnhb3(1,2) = 0.328709624d0
        jcnhb3(2,2) = 0.159717057d0
        jcnhb3(3,2) = 0.342020143d0
        jcnhb3(4,2) = 1d0
        jcnhb3(1,3) = -0.014128431d0
        jcnhb3(2,3) = 0.328709624d0
        jcnhb3(3,3) = 0.328709624d0
        jcnhb3(4,3) = jcnhb3(4,2)
        jcnhb3(1,1) = 0.328709624d0
        jcnhb3(2,1) = 0.328709624d0
        jcnhb3(3,1) = 0.342020143d0
        jcnhb3(4,1) = jcnhb3(4,2)
    endif
end subroutine getpara

!> @brief Evaluate tert for the CH3NH2_APRP potential.
subroutine tert( igrad, ii, x, y, z, nrap, ntap, en, dx1, dy1, dz1 )
    !***********************************************************************
    ! Subroutine to calculate the tertiary contribution
    !***********************************************************************
    implicit none
    integer          :: tk
    integer          :: nrap, ntap
    integer          :: igrad
    double precision :: r1, t2
    double precision :: x(7), y(7), z(7)
    double precision :: xyz(4,3)
    double precision :: dx(22)
    double precision :: rap(nrap,3), tap(ntap,3)
    integer          :: i, ii, j
    integer          :: fg(2)
    double precision :: u, v, f(4)
    double precision :: bndlen, bndang
    double precision :: h(4), en
    double precision :: dh(4,20), dhe(20)
    double precision :: df(4,2)
    double precision :: dtr(2)
    double precision :: dx1(7), dy1(7), dz1(7)
    double precision :: beta
    double precision :: dbeta(12)
    double precision :: dihedr2
    ! The anchor points:
    do i=1, 3
        rap(1,i) = 0.7157485d0
        rap(2,i) = 1.0157485d0
        rap(3,i) = 2.2157485d0
        rap(4,i) = 3.5157485d0
    enddo
    do i=1, 3
        tap(1,i) = -0.193287305d0
        tap(2,i) = 0.482599384d0
        tap(3,i) = 1.0d0
    enddo
    r1 = bndlen(1, 3, x, y, z)
    t2 = dihedr2(1, 3, 4, 2, x, y, z)
    ! calls the 2D tent
    call bigtent(r1, t2, f, rap, nrap, tap, ntap, fg, ii, tk, df)
    ! calcs the contribution at the anchor points
    h = 0d0
    if( tk .eq. 1 )then
        call calchess1(igrad, nrap, ntap, fg, x, y, z, ii, h, dh)
    else if( tk .eq. 2 )then
        call calchess2(igrad, nrap, ntap, fg, x, y, z, ii, h, dh)
    else if( tk .eq. 3 )then
        call calchess3(igrad, nrap, ntap, fg, x, y, z, ii, h, dh)
    endif
    ! Combines the tent functions and contribution from the anchor points
    en = 0d0
    do i=1, 4
        en = en+h(i)*f(i)
    enddo
    dx1 = 0d0
    dy1 = 0d0
    dz1 = 0d0
    !!!!!!!!!!! GET GRADIENT !!!!!!!!!
    if( igrad .eq. 1 )then
        dhe = 0d0
        do j=1, 20
            do i=1, 4
                dhe(j) = dhe(j)+dh(i,j)*f(i)
            enddo
        enddo
        dtr = 0d0
        do j=1, 2
            do i=1, 4
                dtr(j) = dtr(j)+df(i,j)*h(i)
            enddo
        enddo
        do i=1, 20
            dx(i) = dhe(i)
        enddo
        dx(21) = dtr(1)
        dx(22) = dtr(2)
        call convertit(x, y, z, dx, dx1, dy1, dz1)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine tert

!-------------------------------------------------!
!-------------------------------------------------!
!-------------------------------------------------!
!-------------------------------------------------!

!> @brief Evaluate calcHess1 for the CH3NH2_APRP potential.
subroutine calchess1( igrad, nrap, ntap, fg, x, y, z, ii, h, dh )
    ! calcs the tertiary contribution for case 1: corners of the 2D tent
    ! i.e. only need 1 anchor point
    implicit none
    integer          :: igrad
    integer          :: ra1, ra2
    integer          :: nrap, ntap
    double precision :: r1, t2
    double precision :: x(7), y(7), z(7)
    double precision :: qc(20)
    double precision :: e
    double precision :: dx
    integer          :: i, ii
    integer          :: fg(2)
    double precision :: u, v, f(4)
    double precision :: h(4), hs
    double precision :: dhs(20), dh(4,20)
    double precision :: pi
    pi  = 4.d0*datan(1.d0)
    hs  = 0d0
    dhs = 0d0
    dh  = 0d0
    call calcgeom(x, y, z, qc)
    !****** CORNERS
    if( (fg(1)-1 .eq. 0) .and. (fg(2)-1 .eq. 0) )then
        ra1 = 1
        ra2 = 1
        call solver(igrad, ra1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(4) = hs
        do i=1, 20
            dh(4,i) = dhs(i)
        enddo
    else if( (fg(1) .eq. nrap+1) .and. (fg(2) .eq. ntap+1) )then
        ra1 = nrap
        ra2 = ntap
        call solver(igrad, ra1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(1) = hs
        do i=1, 20
            dh(1,i) = dhs(i)
        enddo
    else if( (fg(1)-1 .eq. 0) .and. (fg(2) .eq. ntap+1) )then
        ra1 = 1
        ra2 = ntap
        call solver(igrad, ra1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(3) = hs
        do i=1, 20
            dh(3,i) = dhs(i)
        enddo
    else if( (fg(1) .eq. nrap+1) .and. (fg(2)-1 .eq. 0) )then
        ra1 = nrap
        ra2 = 1
        call solver(igrad, ra1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(2) = hs
        do i=1, 20
            dh(2,i) = dhs(i)
        enddo
    endif
end subroutine calchess1

!> @brief Evaluate calcHess2 for the CH3NH2_APRP potential.
subroutine calchess2( igrad, nrap, ntap, fg, x, y, z, ii, h, dh )
    ! calcs the tertiary contribution for case 2: sides of the 2D tent
    ! i.e. between 2 anchor points
    implicit none
    integer          :: igrad
    integer          :: ra1, ra2
    integer          :: nrap, ntap
    double precision :: r1, t2
    double precision :: x(7), y(7), z(7)
    double precision :: qc(20)
    double precision :: e
    double precision :: dx
    integer          :: i, ii
    integer          :: fg(2)
    double precision :: u, v, f(4)
    double precision :: h(4), hs
    double precision :: dhs(20), dh(4,20)
    double precision :: pi
    pi  = 4.d0*datan(1.d0)
    hs  = 0d0
    dhs = 0d0
    dh  = 0d0
    call calcgeom(x, y, z, qc)
    !****** SIDES
    if( (fg(1)-1 .eq. 0) )then
        ra1 = 1
        ra2 = fg(2)
        call solver(igrad, ra1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(3) = hs
        do i=1, 20
            dh(3,i) = dhs(i)
        enddo
        call solver(igrad, ra1, ra2-1, qc, ii, nrap, ntap, hs, dhs)
        h(4) = hs
        do i=1, 20
            dh(4,i) = dhs(i)
        enddo
    else if( (fg(1) .eq. nrap+1) )then
        ra1 = nrap
        ra2 = fg(2)
        call solver(igrad, ra1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(1) = hs
        do i=1, 20
            dh(1,i) = dhs(i)
        enddo
        call solver(igrad, ra1, ra2-1, qc, ii, nrap, ntap, hs, dhs)
        h(2) = hs
        do i=1, 20
            dh(2,i) = dhs(i)
        enddo
    else if( (fg(2)-1 .eq. 0) )then
        ra2 = 1
        ra1 = fg(1)
        call solver(igrad, ra1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(2) = hs
        do i=1, 20
            dh(2,i) = dhs(i)
        enddo
        call solver(igrad, ra1-1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(4) = hs
        do i=1, 20
            dh(4,i) = dhs(i)
        enddo
    else if( (fg(2) .eq. ntap+1) )then
        ra2 = ntap
        ra1 = fg(1)
        call solver(igrad, ra1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(1) = hs
        do i=1, 20
            dh(1,i) = dhs(i)
        enddo
        call solver(igrad, ra1-1, ra2, qc, ii, nrap, ntap, hs, dhs)
        h(3) = hs
        do i=1, 20
            dh(3,i) = dhs(i)
        enddo
    endif
end subroutine calchess2

!> @brief Evaluate calcHess3 for the CH3NH2_APRP potential.
subroutine calchess3( igrad, nrap, ntap, fg, x, y, z, ii, h, dh )
    ! calcs the tertiary contribution for case 3: between 4 anchor points
    implicit none
    integer          :: igrad
    integer          :: ra1, ra2
    integer          :: nrap, ntap
    double precision :: r1, t2
    double precision :: x(7), y(7), z(7)
    double precision :: qc(20)
    double precision :: e
    double precision :: dx
    integer          :: i, ii
    integer          :: fg(2)
    double precision :: u, v, f(4)
    double precision :: h(4), hs
    double precision :: dhs(20), dh(4,20)
    double precision :: pi
    pi  = 4.d0*datan(1.d0)
    hs  = 0d0
    dhs = 0d0
    dh  = 0d0
    ra1 = fg(1)
    ra2 = fg(2)
    call calcgeom(x, y, z, qc)
    call solver(igrad, ra1, ra2, qc, ii, nrap, ntap, hs, dhs)
    h(1) = hs
    do i=1, 20
        dh(1,i) = dhs(i)
    enddo
    call solver(igrad, ra1, ra2-1, qc, ii, nrap, ntap, hs, dhs)
    h(2) = hs
    do i=1, 20
        dh(2,i) = dhs(i)
    enddo
    call solver(igrad, ra1-1, ra2, qc, ii, nrap, ntap, hs, dhs)
    h(3) = hs
    do i=1, 20
        dh(3,i) = dhs(i)
    enddo
    call solver(igrad, ra1-1, ra2-1, qc, ii, nrap, ntap, hs, dhs)
    h(4) = hs
    do i=1, 20
        dh(4,i) = dhs(i)
    enddo
end subroutine calchess3

!> @brief Evaluate calcgeom for the CH3NH2_APRP potential.
subroutine calcgeom( x, y, z, qc )
    ! get the internal coords for tertiary calcs
    implicit none
    integer          :: i
    double precision :: x(7), y(7), z(7)
    double precision :: qc(20)
    double precision :: bndlen, bndang
    double precision :: dihedr2
    qc(1)  = bndlen(2, 3, x, y, z)
    qc(2)  = bndlen(3, 4, x, y, z)
    qc(3)  = bndlen(4, 5, x, y, z)
    qc(4)  = bndlen(4, 6, x, y, z)
    qc(5)  = bndlen(4, 7, x, y, z)
    qc(6)  = bndang(1, 3, 4, x, y, z)
    qc(7)  = bndang(2, 3, 4, x, y, z)
    qc(8)  = bndang(3, 4, 5, x, y, z)
    qc(9)  = bndang(3, 4, 6, x, y, z)
    qc(10) = bndang(3, 4, 7, x, y, z)
    qc(11) = dihedr2(1, 3, 4, 5, x, y, z)
    qc(12) = dihedr2(2, 3, 4, 5, x, y, z)
    qc(13) = bndlen(1, 5, x, y, z)
    qc(14) = bndlen(2, 5, x, y, z)
    qc(15) = bndang(5, 4, 6, x, y, z)
    qc(16) = bndang(6, 4, 7, x, y, z)
    qc(17) = bndang(7, 4, 5, x, y, z)
    qc(18) = bndlen(5, 6, x, y, z)
    qc(19) = bndlen(6, 7, x, y, z)
    qc(20) = bndlen(7, 5, x, y, z)
end subroutine calcgeom

!> @brief Evaluate bigtent for the CH3NH2_APRP potential.
subroutine bigtent( r1, t2, f, rap, nrap, tap, ntap, fg, si, tk, df )
    !***********************************************************************
    ! Subroutine to evaluate 2D tent functions
    !***********************************************************************
    implicit none
    integer          :: i, si, s, tk
    integer          :: nrap !number of anchor points
    integer          :: ntap !number of anchor points
    double precision :: r1, t2
    double precision :: rap(nrap,3)
    double precision :: tap(ntap,3)
    double precision :: f(4)
    double precision :: un, up, vn, vp, dnom
    double precision :: df(4,2)
    double precision :: dun, dup, dvn, dvp, ddnom(2)
    integer          :: fg(2)
    f  = 0d0
    fg = 0
    ! Check whether the point are end points
    if( r1 .lt. rap(1,si) )then
        fg(1) = 1
    else if( r1 .ge. rap(nrap,si) )then
        fg(1) = nrap+1
    endif
    do i=2, nrap
        if( (rap(i-1,si) .le. r1) .and. (r1 .lt. rap(i,si)) )then
            fg(1) = i
        endif
    enddo
    ! Check whether the point are end points
    if( t2 .lt. tap(1,si) )then
        fg(2) = 1
    else if( t2 .ge. tap(ntap,si) )then
        fg(2) = ntap+1
    endif
    do i=2, ntap
        if( (tap(i-1,si) .le. t2) .and. (t2 .lt. tap(i,si)) )then
            fg(2) = i
        endif
    enddo
    un = 0d0
    up = 0d0
    vn = 0d0
    vp = 0d0
    !****** CORNERS
    if( (fg(1)-1 .eq. 0) .and. (fg(2)-1 .eq. 0) )then
        vp = 1d0
        up = 1d0
        tk = 1
    else if( (fg(1) .eq. nrap+1) .and. (fg(2) .eq. ntap+1) )then
        un = 1d0
        vn = 1d0
        tk = 1
    else if( (fg(1)-1 .eq. 0) .and. (fg(2) .eq. ntap+1) )then
        up = 1d0
        vn = 1d0
        tk = 1
    else if( (fg(1) .eq. nrap+1) .and. (fg(2)-1 .eq. 0) )then
        un = 1d0
        vp = 1d0
        tk = 1
        !****** SIDES
    else if( (fg(1)-1 .eq. 0) )then
        up = (r1-rap(fg(1),si))**4d0
        vn = (t2-tap(fg(2)-1,si))**4d0
        vp = (t2-tap(fg(2),si))**4d0
        tk = 2
    else if( (fg(1) .eq. nrap+1) )then
        un = (r1-rap((fg(1)-1),si))**4d0
        vn = (t2-tap(fg(2)-1,si))**4d0
        vp = (t2-tap(fg(2),si))**4d0
        ! Connection for r1=rap
        if( r1 .eq. rap((fg(1)-1),si) )then
            un = 1d0
        endif
        tk = 2
    else if( (fg(2)-1 .eq. 0) )then
        un = (r1-rap((fg(1)-1),si))**4d0
        up = (r1-rap(fg(1),si))**4d0
        vp = (t2-tap(fg(2),si))**4d0
        tk = 2
    else if( (fg(2) .eq. ntap+1) )then
        un = (r1-rap((fg(1)-1),si))**4d0
        up = (r1-rap(fg(1),si))**4d0
        vn = (t2-tap(fg(2)-1,si))**4d0
        ! Connection for t2=rap
        if( t2 .eq. tap(fg(2)-1,si) )then
            vn = 1d0
        endif
        tk = 2
        !****** MIDDLE
    else
        un = (r1-rap((fg(1)-1),si))**4d0
        up = (r1-rap(fg(1),si))**4d0
        vn = (t2-tap(fg(2)-1,si))**4d0
        vp = (t2-tap(fg(2),si))**4d0
        tk = 3
    endif
    dnom = un*vn + un*vp + up*vn + up*vp
    f(1) = un*vn/dnom
    f(2) = un*vp/dnom
    f(3) = up*vn/dnom
    f(4) = up*vp/dnom
    df   = 0d0
    if( (fg(1)-1 .eq. 0) .and. (fg(2)-1 .eq. 0) )then
    else if( (fg(1) .eq. nrap+1) .and. (fg(2) .eq. ntap+1) )then
    else if( (fg(1)-1 .eq. 0) .and. (fg(2) .eq. ntap+1) )then
    else if( (fg(1) .eq. nrap+1) .and. (fg(2)-1 .eq. 0) )then
        ! THEN nothing for corners
        !****** CORNERS
        ! all derivatives equal zero
        !****** SIDES
        !      if( (fg(1)-1.eq.0) ) then
    else if( (fg(1)-1 .eq. 0) )then
        dup = 4d0*(r1-rap(fg(1),si))**3d0
        dvn = 4d0*(t2-tap(fg(2)-1,si))**3d0
        dvp = 4d0*(t2-tap(fg(2),si))**3d0
    else if( (fg(1) .eq. nrap+1) )then
        dun = 4d0*(r1-rap((fg(1)-1),si))**3d0
        dvn = 4d0*(t2-tap(fg(2)-1,si))**3d0
        dvp = 4d0*(t2-tap(fg(2),si))**3d0
    else if( (fg(2)-1 .eq. 0) )then
        dun = 4d0*(r1-rap((fg(1)-1),si))**3d0
        dup = 4d0*(r1-rap(fg(1),si))**3d0
        dvp = 4d0*(t2-tap(fg(2),si))**3d0
    else if( (fg(2) .eq. nrap+1) )then
        dun = 4d0*(r1-rap((fg(1)-1),si))**3d0
        dup = 4d0*(r1-rap(fg(1),si))**3d0
        dvn = 4d0*(t2-tap(fg(2)-1,si))**3d0
        !****** MIDDLE
    else
        dun = 4d0*(r1-rap((fg(1)-1),si))**3d0
        dup = 4d0*(r1-rap(fg(1),si))**3d0
        dvn = 4d0*(t2-tap(fg(2)-1,si))**3d0
        dvp = 4d0*(t2-tap(fg(2),si))**3d0
    endif
    ddnom(1) = dun*vn + dun*vp + dup*vn + dup*vp
    ddnom(2) = un*dvn + un*dvp + up*dvn + up*dvp
    df(1,1)  = ((dun*vn)*dnom-ddnom(1)*un*vn)/(dnom**2d0)
    df(2,1)  = ((dun*vp)*dnom-ddnom(1)*un*vp)/(dnom**2d0)
    df(3,1)  = ((dup*vn)*dnom-ddnom(1)*up*vn)/(dnom**2d0)
    df(4,1)  = ((dup*vp)*dnom-ddnom(1)*up*vp)/(dnom**2d0)
    df(1,2)  = ((un*dvn)*dnom-ddnom(2)*un*vn)/(dnom**2d0)
    df(2,2)  = ((un*dvp)*dnom-ddnom(2)*un*vp)/(dnom**2d0)
    df(3,2)  = ((up*dvn)*dnom-ddnom(2)*up*vn)/(dnom**2d0)
    df(4,2)  = ((up*dvp)*dnom-ddnom(2)*up*vp)/(dnom**2d0)
end subroutine bigtent

!> @brief Evaluate KPev2gm2 for the CH3NH2_APRP potential.
subroutine kpev2gm2( igrad, x, y, z, v3, imol, dx, dy, dz )
    ! gives the primary coordinate contribution
    implicit none
    integer          :: imol
    double precision :: a12(4)
    double precision :: const(4)
    double precision :: ga1(4), ga2(4), ga3(4)
    double precision :: gb1(4), gb2(4), gb3(4)
    double precision :: ka1(4), ka2(4)
    double precision :: c(0:10)
    double precision :: re, de, rmin, pi, rx
    double precision :: fy, dfdr, u, dydr, lj, dljdr
    double precision :: rr, tt, dv2dr, dv2dt, v3
    integer          :: ri, i, igrad
    double precision :: x(7), y(7), z(7)
    double precision :: dx(7), dy(7), dz(7)
    double precision :: ae, e1, e2, e3, f1, f2, f3
    double precision :: g1, g2, g3, g4, i2, i3, i4
    double precision :: h2, h3, h4, j2, j3, j4
    double precision :: n1, n2, n3, n4
    double precision :: bndlen
    double precision :: dihedr2
    double precision :: ak, tt0, gr1, gr2, gr3
    double precision :: fr1(5), fr2(5), fr3(5)
    double precision :: ft2(5), ft3(5)
    double precision :: dihedrix, dihedriy, dihedriz
    double precision :: dihedrjx, dihedrjy, dihedrjz
    double precision :: dihedrkx, dihedrky, dihedrkz
    double precision :: dihedrlx, dihedrly, dihedrlz
    dv2dr = 0d0
    dv2dt = 0d0
    v3    = 0d0
    dx    = 0d0
    dy    = 0d0
    dz    = 0d0
    rr    = bndlen(1, 3, x, y, z)
    tt    = dihedr2(1, 3, 4, 2, x, y, z)
    if( imol .eq. 1 )then! GROUND STATE
        e1 = 6.420621d0
        e2 = -1.85d0
        e3 = 1.01574847d0
        f1 = 0.6d0
        f2 = -4.0d0
        f3 = 2.23d0
        v3 = e1*(1.0d0-exp(e2*((rr-e3))))**2d0 +f1*exp(f2*((rr-f3)**2d0))
        if( igrad .eq. 1 )then
            dv2dr = -2.0d0*e1*e2*exp(e2*(rr-e3)) *(1-exp(e2*(rr-e3))) + 2d0*f1*f2*(rr-f3)*exp(f2*(rr-f3)**2d0)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ak  = 0.945760042d0
        tt0 = 0.482599385d0
        gr1 = 1.472234651d0
        gr2 = -5.656037749d0
        gr3 = 0.647196191d0
        v3  = v3+ ak*(tt-tt0)**2d0* gr1*exp(gr2*(rr-gr3)**2d0)
        if( igrad .eq. 1 )then
            dv2dr = dv2dr+ ak*(tt-tt0)**2d0* 2d0*gr1*gr2*(rr-gr3)*exp(gr2*(rr-gr3)**2d0)
            dv2dt = ak*(tt-tt0)*2d0* gr1*exp(gr2*(rr-gr3)**2d0)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fr1(1) = 3.067752289d0
        fr1(2) = 1076.860912d0
        fr1(3) = 0.792499074d0
        fr1(4) = -1.373902904d0
        fr2(1) = -0.845025076d0
        fr2(2) = -5.474618306d0
        fr2(3) = -5.214819223d0
        fr2(4) = -7.271506597d0
        fr3(1) = 0.98801149d0
        fr3(2) = 0.89839555d0
        fr3(3) = 1.484837694d0
        fr3(4) = 1.89501344d0
        ft2(1) = -3.340905161d0
        ft2(2) = -14.55654461d0
        ft2(3) = -12.70986496d0
        ft2(4) = -11.06393037d0
        ft3(1) = -0.917131468d0
        ft3(2) = -1.489784616d0
        ft3(3) = 1.002868947d0
        ft3(4) = -0.334959847d0
        do i=1, 4
            v3 = v3+ fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))* exp(ft2(i)*((tt-ft3(i))**2d0))
        enddo
        if( igrad .eq. 1 )then
            do i=1, 4
                dv2dr=dv2dr+ 2d0*fr2(i)*(rr-fr3(i))* fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))* exp(ft2(i)*((tt-ft3(i))** &
                    & 2d0))
            enddo
            do i=1, 4
                dv2dt=dv2dt+ 2d0*ft2(i)*(tt-ft3(i))* fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))* exp(ft2(i)*((tt-ft3(i))** &
                    & 2d0))
            enddo
        endif
    else if( imol .eq. 2 )then! EXCITED STATE
        e1 = 4.509566d0
        n4 = 1.650445177d0
        e2 = -1.460739955d0
        e3 = 0.582944229d0
        f1 = 0.234487181d0
        f2 = -2.282322532d0
        f3 = 0.30831821d0
        n1 = 0.0d0
        n2 = 300d0
        n3 = 1.025d0
        v3 = e1+n4*exp(e2*(rr-e3)) +(n1+n2*(rr-n3)**2d0) *f1*exp(f2*((rr-f3)**2d0))
        if( igrad .eq. 1 )then
            dv2dr=e2*n4*exp(e2*(rr-e3)) + (n2*(rr-n3)*2d0) *f1*exp(f2*((rr-f3)**2d0)) +(n1+n2*(rr-n3)**2d0) *2d0*f1* &
                & f2* (rr-f3)*exp(f2*(rr-f3)**2d0)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ak  = 1.233643227d0
        tt0 = 1.0d0
        gr1 = 0.978308804d0
        gr2 = -0.981638325d0
        gr3 = 0.901675465d0
        v3  = v3+ ak*(tt-tt0)**2d0* gr1*exp(gr2*(rr-gr3)**2d0)
        if( igrad .eq. 1 )then
            dv2dr = dv2dr+ ak*(tt-tt0)**2d0* 2d0*gr1*gr2*(rr-gr3)*exp(gr2*(rr-gr3)**2d0)
            dv2dt = ak*(tt-tt0)*2d0* gr1*exp(gr2*(rr-gr3)**2d0)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fr1(1) = 5.402431748d0
        fr1(2) = -5.542582536d0
        fr1(3) = -3.460058164d0
        fr1(4) = -0.508881823d0
        fr1(5) = 3.692631232d0
        fr2(1) = -5.951923443d0
        fr2(2) = -6.031843482d0
        fr2(3) = -9.434511421d0
        fr2(4) = -2.878092095d0
        fr2(5) = -9.180177476d0
        fr3(1) = 1.867955914d0
        fr3(2) = 1.870582014d0
        fr3(3) = 1.787229743d0
        fr3(4) = 2.302428489d0
        fr3(5) = 1.746230571d0
        ft2(1) = -9.559209372d0
        ft2(2) = -9.696891789d0
        ft2(3) = -5.670799519d0
        ft2(4) = -8.901043282d0
        ft2(5) = -5.346117363d0
        ft3(1) = -0.590397385d0
        ft3(2) = -0.675159631d0
        ft3(3) = 0.465695686d0
        ft3(4) = -0.5580416d0
        ft3(5) = 0.427418036d0
        do i=1, 5
            v3 = v3+ fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))* exp(ft2(i)*((tt-ft3(i))**2d0))
        enddo
        if( igrad .eq. 1 )then
            do i=1, 5
                dv2dr=dv2dr+ 2d0*fr2(i)*(rr-fr3(i))* fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))* exp(ft2(i)*((tt-ft3(i))** &
                    & 2d0))
            enddo
            do i=1, 5
                dv2dt=dv2dt+ 2d0*ft2(i)*(tt-ft3(i))* fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))* exp(ft2(i)*((tt-ft3(i))** &
                    & 2d0))
            enddo
        endif
    else if( imol .eq. 3 )then! COUPLING
        a12(1)   = -0.16290134d0
        ga1(1)   = -0.231270633d0
        ga2(1)   = -42.62845833d0
        ga3(1)   = -0.788781464d0
        gb1(1)   = -0.157578396d0
        gb2(1)   = -4.728738508d0
        gb3(1)   = -0.492785959d0
        const(1) = 0.050222306d0
        a12(2)   = 1.514877043d0
        ga1(2)   = 0.38885665d0
        ga2(2)   = -52.43604225d0
        ga3(2)   = -0.596932233d0
        gb1(2)   = 2.340167593d0
        gb2(2)   = -0.158293136d0
        gb3(2)   = 0.03631615d0
        const(2) = -2.497398102d0
        a12(3)   = 0.96486818d0
        ga1(3)   = -0.567995439d0
        ga2(3)   = -14.16424012d0
        ga3(3)   = -0.534406038d0
        gb1(3)   = 1.682033257d0
        gb2(3)   = -0.103900297d0
        gb3(3)   = 0.059650329d0
        const(3) = -1.859104948d0
        a12(4)   = 0.730431712d0
        ga1(4)   = -1.39153824d0
        ga2(4)   = -0.081123542d0
        ga3(4)   = -1.342253144d0
        gb1(4)   = 4.257325153d0
        gb2(4)   = -1.963596235d0
        gb3(4)   = -2.057979708d0
        const(4) = 0.669000566d0
        ka1(1)   = -2.904320563d0
        ka1(2)   = -5.709140977d0
        ka1(3)   = -5.612288257d0
        ka1(4)   = -2.876200866d0
        ka2(1)   = 0.382538122d0
        ka2(2)   = 1.81713229d0
        ka2(3)   = 2.40922499d0
        ka2(4)   = 3.270326934d0
        do i=1, 4
            v3=v3 + (a12(i)*sqrt(1.1d0-(tt**2d0)) +const(i) + ga1(i)*exp(ga2(i)*((tt-ga3(i))**2d0)) + gb1(i)* &
                & exp(gb2(i)* ((tt-gb3(i))**2d0))) *(exp(ka1(i)*(rr-ka2(i))**2d0))
        enddo
        if( igrad .eq. 1 )then
            do i=1, 4
                dv2dt=dv2dt + (-1.0d0*a12(i)*tt *(1.0d0/sqrt(1.1d0-(tt**2d0))) + ga1(i)*exp(ga2(i)*((tt-ga3(i))** &
                    & 2d0)) * ga2(i)*((tt-ga3(i))*2d0) + gb1(i)*exp(gb2(i)*((tt-gb3(i))**2d0)) *gb2(i)*((tt-gb3(i))* &
                    & 2d0)) * (exp(ka1(i)*(rr-ka2(i))**2d0))
            enddo
            do i=1, 4
                dv2dr=dv2dr + (a12(i)*sqrt(1.1d0-(tt**2d0)) +const(i) + ga1(i)*exp(ga2(i)*((tt-ga3(i))**2d0)) + &
                    & gb1(i)* exp(gb2(i)*((tt-gb3(i))**2d0))) *(exp(ka1(i)*(rr-ka2(i))**2d0)) *(ka1(i)*(rr-ka2(i))*2d0)
            enddo
        endif
    endif
    dx(1) = dv2dr*((x(1)-x(3))/rr)
    dy(1) = dv2dr*((y(1)-y(3))/rr)
    dz(1) = dv2dr*((z(1)-z(3))/rr)
    dx(3) = dv2dr*((x(3)-x(1))/rr)
    dy(3) = dv2dr*((y(3)-y(1))/rr)
    dz(3) = dv2dr*((z(3)-z(1))/rr)
    dx(1) = dx(1)+dv2dt*dihedrix(1, 3, 4, 2, x, y, z)
    dy(1) = dy(1)+dv2dt*dihedriy(1, 3, 4, 2, x, y, z)
    dz(1) = dz(1)+dv2dt*dihedriz(1, 3, 4, 2, x, y, z)
    dx(2) = dv2dt*dihedrlx(1, 3, 4, 2, x, y, z)
    dy(2) = dv2dt*dihedrly(1, 3, 4, 2, x, y, z)
    dz(2) = dv2dt*dihedrlz(1, 3, 4, 2, x, y, z)
    dx(3) = dx(3)+dv2dt*dihedrjx(1, 3, 4, 2, x, y, z)
    dy(3) = dy(3)+dv2dt*dihedrjy(1, 3, 4, 2, x, y, z)
    dz(3) = dz(3)+dv2dt*dihedrjz(1, 3, 4, 2, x, y, z)
    dx(4) = dv2dt*dihedrkx(1, 3, 4, 2, x, y, z)
    dy(4) = dv2dt*dihedrky(1, 3, 4, 2, x, y, z)
    dz(4) = dv2dt*dihedrkz(1, 3, 4, 2, x, y, z)
end subroutine kpev2gm2

!> @brief Evaluate dihedr2 for the CH3NH2_APRP potential.
double precision function dihedr2( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Evaluate the dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    ! modified version of dihedr1 to work through flat amine geoms
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: kkx, kky, kkz, yyx, yyy, yyz
    double precision            :: llx, lly, llz, eex, eey, eez
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    kkx    = x(i)-x(j)
    kky    = y(i)-y(j)
    kkz    = z(i)-z(j)
    yyx    = x(k)-x(j)
    yyy    = y(k)-y(j)
    yyz    = z(k)-z(j)
    llx    = x(l)-x(k)
    lly    = y(l)-y(k)
    llz    = z(l)-z(k)
    eex    = x(j)-x(k)
    eey    = y(j)-y(k)
    eez    = z(j)-z(k)
    oo1(1) = kky*yyz-kkz*yyy
    oo1(2) = kkz*yyx-kkx*yyz
    oo1(3) = kkx*yyy-kky*yyx
    oo2(1) = lly*eez-llz*eey
    oo2(2) = llz*eex-llx*eez
    oo2(3) = llx*eey-lly*eex
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       ra=sab
    !       dihedr2=acos(ra)
    dihedr2 = ra
    !--------------------------------------------------!
end function dihedr2

!> @brief Evaluate dihedrIX for the CH3NH2_APRP potential.
double precision function dihedrix( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !       dihedr2=ra
    !------------------Xi---------------------------
    !top (sab deriv for Xi)
    td = -1d0*kjz*lkz*jkx +kjz*lkx*jkz +kjy*lkx*jky -kjy*lky*jkx
    !bottom (p1p2 deriv for Xi)
    ! deriv of sso1 and bd for Xi
    ld = 2d0*(x(i)-x(j))*(kjz**2d0+kjy**2d0) -2*(kjz*ijz*kjx +kjy*ijy*kjx)
    bd = (sso2*ld)/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrIX=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrix = id!/(1d0-ra**2d0)**(0.5d0)
    !------------------Xi---------------------------
end function dihedrix

!> @brief Evaluate dihedrIY for the CH3NH2_APRP potential.
double precision function dihedriy( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Yi---------------------------
    !top (sab deriv for Yi)
    td = -1d0*kjz*lkz*jky +kjz*lky*jkz +kjx*lky*jkx -kjx*lkx*jky
    ! deriv of sso1 and bd for Yi
    ld = 2d0*(y(i)-y(j))*(kjz**2+kjx**2) -2*(kjz*ijz*kjy +kjx*ijx*kjy)
    bd = (sso2*ld)/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrIY=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedriy = id!/(1d0-ra**2d0)**(0.5d0)
    !------------------Yi---------------------------
end function dihedriy

!> @brief Evaluate dihedrIZ for the CH3NH2_APRP potential.
double precision function dihedriz( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Zi---------------------------
    !top (sab deriv for Zi)
    td = -1d0*kjy*lky*jkz +kjy*lkz*jky +kjx*lkz*jkx -kjx*lkx*jkz
    ! deriv of sso1 and bd for Zi
    ld = 2d0*(z(i)-z(j))*(kjx**2+kjy**2) -2*(kjx*ijx*kjz +kjy*ijy*kjz)
    bd = (sso2*ld)/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrIZ=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedriz = id!/(1d0-ra**2d0)**(0.5d0)
    !------------------Zi---------------------------
end function dihedriz

!> @brief Evaluate dihedrLX for the CH3NH2_APRP potential.
double precision function dihedrlx( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Xl---------------------------
    !top (sab deriv for Xl)
    td = -1d0*ijz*kjx*jkz +ijx*kjz*jkz +ijx*kjy*jky -ijy*kjx*jky
    ! deriv of sso2 for Xl
    rd = 2d0*(x(l)-x(k))*(jky**2+jkz**2) -2*(jky*lky*jkx +jkz*lkz*jkx)
    bd = (sso1*rd)/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrLX=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrlx = id!/(1d0-ra**2d0)**(0.5d0)
    !------------------Xl---------------------------
end function dihedrlx

!> @brief Evaluate dihedrLY for the CH3NH2_APRP potential.
double precision function dihedrly( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Yl---------------------------
    !top (sab deriv for Yl)
    td = -1d0*ijz*kjy*jkz +ijy*kjz*jkz +ijy*kjx*jkx -ijx*kjy*jkx
    ! deriv of sso2 for Yl
    rd = 2d0*(y(l)-y(k))*(jkx**2+jkz**2) -2*(jkx*lkx*jky +jkz*lkz*jky)
    bd = (sso1*rd)/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrLY=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrly = id!/(1d0-ra**2d0)**(0.5d0)
    !------------------Yl---------------------------
end function dihedrly

!> @brief Evaluate dihedrLZ for the CH3NH2_APRP potential.
double precision function dihedrlz( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Zl---------------------------
    !top (sab deriv for Zl)
    td = -1d0*ijy*kjz*jky +ijz*kjy*jky +ijz*kjx*jkx -ijx*kjz*jkx
    ! deriv of sso2 and bd for Zl
    rd = 2d0*(z(l)-z(k))*(jky**2+jkx**2) -2*(jky*lky*jkz +jkx*lkx*jkz)
    bd = (sso1*rd)/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrLZ=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrlz = id!/(1d0-ra**2d0)**(0.5d0)
    !------------------Zl---------------------------
end function dihedrlz

!> @brief Evaluate dihedrJX for the CH3NH2_APRP potential.
double precision function dihedrjx( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Xj---------------------------
    !top (sab deriv for Xj)
    td=(-2d0*x(j)+2d0*x(k))*ijz*lkz +ijz*lkx*jkz -(-2d0*x(j)+x(k)+x(i))*kjz*lkz -kjz*lkx*jkz -kjy*lkx*jky -(-2d0*x(j)+ &
        & x(k)+x(i))*kjy*lky +ijy*lkx*jky +(-2d0*x(j)+2d0*x(k))*ijy*lky
    ! deriv of sso1 for Xj
    ld = 2d0*(x(j)-x(k))*(ijz**2+ijy**2) +2d0*(x(j)-x(i))*(kjz**2+kjy**2) -2*(2d0*x(j)-x(k)-x(i))*(ijz*kjz +ijy*kjy)
    ! deriv of sso2 for Xj
    rd = 2d0*(x(j)-x(k))*(lkz**2+lky**2) -2*(lkz*lkx*jkz +lky*lkx*jky)
    bd = ((ld*sso2)+(sso1*rd))/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrJX=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrjx = id!/(1d0-ra**2d0)**(0.5d0)
    !------------------Xj---------------------------
end function dihedrjx

!> @brief Evaluate dihedrJY for the CH3NH2_APRP potential.
double precision function dihedrjy( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Yj---------------------------
    !top (sab deriv for Yj)
    td=-kjz*lky*jkz -(-2d0*y(j)+y(k)+y(i))*kjz*lkz +ijz*lky*jkz +(-2d0*y(j)+2d0*y(k))*ijz*lkz +(-2d0*y(j)+2d0*y(k))* &
        & ijx* lkx +ijx*lky*jkx -(-2d0*y(j)+y(k)+y(i))*kjx*lkx -kjx*lky*jkx
    ! deriv of sso1 for Yj
    ld = 2d0*(y(j)-y(k))*(ijz**2+ijx**2) +2d0*(y(j)-y(i))*(kjz**2+kjx**2) -2*(2d0*y(j)-y(k)-y(i))*(ijz*kjz +ijx*kjx)
    ! deriv of sso2 for Yj
    rd = 2d0*(y(j)-y(k))*(lkz**2+lkx**2) -2*(lkz*lky*jkz +lkx*lky*jkx)
    bd = ((ld*sso2)+(sso1*rd))/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrJY=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrjy = id!/(1d0-ra**2d0)**(0.5d0)
    !------------------Yj---------------------------
end function dihedrjy

!> @brief Evaluate dihedrJZ for the CH3NH2_APRP potential.
double precision function dihedrjz( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Zj---------------------------
    !top (sab deriv for Zj)
    td=(-2d0*z(j)+2d0*z(k))*ijy*lky +ijy*lkz*jky -(-2d0*z(j)+z(k)+z(i))*kjy*lky -kjy*lkz*jky -kjx*lkz*jkx -(-2d0*z(j)+ &
        & z(k)+z(i))*kjx*lkx +ijx*lkz*jkx +(-2d0*z(j)+2d0*z(k))*ijx*lkx
    ! deriv of sso1 for Zj
    ld = 2d0*(z(j)-z(k))*(ijy**2+ijx**2) +2d0*(z(j)-z(i))*(kjy**2+kjx**2) -2*(2d0*z(j)-z(k)-z(i))*(ijy*kjy +ijx*kjx)
    ! deriv of sso2 for Zj
    rd = 2d0*(z(j)-z(k))*(lky**2+lkx**2) -2*(lky*lkz*jky +lkx*lkz*jkx)
    bd = ((ld*sso2)+(sso1*rd))/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrJZ=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrjz = id!/(1d0-ra**2d0)**(0.5d0)
    !------------------Zj---------------------------
end function dihedrjz

!> @brief Evaluate dihedrKX for the CH3NH2_APRP potential.
double precision function dihedrkx( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Xk---------------------------
    !top (sab deriv for Xk)
    td=(-2d0*x(k)+2d0*x(j))*ijz*lkz -(-2d0*x(k)+x(j)+x(l))*ijz*jkz +ijx*kjz*lkz -ijx*kjz*jkz -ijx*kjy*jky +ijx*kjy* &
        & lky - (-2d0*x(k)+x(j)+x(l))*ijy*jky +(-2d0*x(k)+2d0*x(j))*ijy*lky
    ! deriv of sso1 for Xk
    ld = 2d0*(x(k)-x(j))*(ijz**2+ijy**2) -2*(ijx*ijz*kjz +ijx*ijy*kjy)
    ! deriv of sso2 for Xk
    rd = 2d0*(x(k)-x(j))*(lkz**2+lky**2) +2d0*(x(k)-x(l))*(jkz**2+jky**2) -2*(2d0*x(k)-x(j)-x(l))*(lkz*jkz +lky*jky)
    bd = ((ld*sso2)+(sso1*rd))/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrKX=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrkx = id!(-id)/(1d0-ra**2d0)**(0.5d0)
    !------------------Xk---------------------------
end function dihedrkx

!> @brief Evaluate dihedrKY for the CH3NH2_APRP potential.
double precision function dihedrky( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !        i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Yk---------------------------
    !top (sab deriv for Yk)
    td=(-2d0*y(k)+2d0*y(j))*ijz*lkz -(-2d0*y(k)+y(j)+y(l))*ijz*jkz +ijy*kjz*lkz -ijy*kjz*jkz -ijy*kjx*jkx +ijy*kjx* &
        & lkx - (-2d0*y(k)+y(j)+y(l))*ijx*jkx +(-2d0*y(k)+2d0*y(j))*ijx*lkx
    ! deriv of sso1 for Yk
    ld = 2d0*(y(k)-y(j))*(ijz**2+ijx**2) -2*(ijy*ijz*kjz +ijy*ijx*kjx)
    ! deriv of sso2 for Yk
    rd = 2d0*(y(k)-y(j))*(lkz**2+lkx**2) +2d0*(y(k)-y(l))*(jkz**2+jkx**2) -2*(2d0*y(k)-y(j)-y(l))*(lkz*jkz +lkx*jkx)
    bd = ((ld*sso2)+(sso1*rd))/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrKY=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrky = id!(-id)/(1d0-ra**2d0)**(0.5d0)
    !------------------Yk---------------------------
end function dihedrky

!> @brief Evaluate dihedrKZ for the CH3NH2_APRP potential.
double precision function dihedrkz( i, j, k, l, x, y, z )
    !***********************************************************************
    ! Converts derivative from internal to XYZ for
    ! dihedral angle between atoms i, j, k, and l
    ! Input: i,j,k,l (indices of 4 atoms)
    !         i--j--k--l
    !***********************************************************************
    implicit none
    integer                     :: i, j, k, l
    double precision, parameter :: tiny = 1.0d-13
    double precision, parameter :: pi   = 3.1415926535898d0
    double precision            :: x(*), y(*), z(*)
    double precision            :: eij(3), ejk(3), ekl(3), u(3), v(3), w(3)
    double precision            :: rij, rjk, rkl, aijk, ajkl, wejk
    double precision            :: ijx, ijy, ijz, kjx, kjy, kjz
    double precision            :: lkx, lky, lkz, jkx, jky, jkz
    double precision            :: oo1(3), oo2(3), so1(3), so2(3)
    double precision            :: sso1, sso2, rso1, rso2
    double precision            :: ab(3)
    double precision            :: sab, p1p2, ra
    double precision            :: ld, rd, td, bd, id, fd
    !   eij
    ! i-----j
    !   aijk \ ejk
    !         \
    !          \ ajkl
    !           k----l
    !             kl
    ijx    = x(i)-x(j)
    ijy    = y(i)-y(j)
    ijz    = z(i)-z(j)
    kjx    = x(k)-x(j)
    kjy    = y(k)-y(j)
    kjz    = z(k)-z(j)
    lkx    = x(l)-x(k)
    lky    = y(l)-y(k)
    lkz    = z(l)-z(k)
    jkx    = x(j)-x(k)
    jky    = y(j)-y(k)
    jkz    = z(j)-z(k)
    oo1(1) = ijy*kjz-ijz*kjy
    oo1(2) = ijz*kjx-ijx*kjz
    oo1(3) = ijx*kjy-ijy*kjx
    oo2(1) = lky*jkz-lkz*jky
    oo2(2) = lkz*jkx-lkx*jkz
    oo2(3) = lkx*jky-lky*jkx
    so1(1) = oo1(1)**2d0
    so1(2) = oo1(2)**2d0
    so1(3) = oo1(3)**2d0
    sso1   = so1(1)+so1(2)+so1(3)
    rso1   = (sso1)**(0.5d0)
    so2(1) = oo2(1)**2d0
    so2(2) = oo2(2)**2d0
    so2(3) = oo2(3)**2d0
    sso2   = so2(1)+so2(2)+so2(3)
    rso2   = (sso2)**(0.5d0)
    ab(1)  = oo1(1)*oo2(1)
    ab(2)  = oo1(2)*oo2(2)
    ab(3)  = oo1(3)*oo2(3)
    sab    = ab(1)+ab(2)+ab(3)
    p1p2   = rso1*rso2
    ra     = sab/p1p2
    !       dihedr2=acos(ra)
    !------------------Zk---------------------------
    !top (sab deriv for Zk)
    td=(-2d0*z(k)+2d0*z(j))*ijy*lky -(-2d0*z(k)+z(j)+z(l))*ijy*jky +ijz*kjy*lky -ijz*kjy*jky -ijz*kjx*jkx +ijz*kjx* &
        & lkx - (-2d0*z(k)+z(j)+z(l))*ijx*jkx +(-2d0*z(k)+2d0*z(j))*ijx*lkx
    ! deriv of sso1 for Zk
    ld = 2d0*(z(k)-z(j))*(ijy**2+ijx**2) -2*(ijz*ijy*kjy +ijz*ijx*kjx)
    ! deriv of sso2 for Zk
    rd = 2d0*(z(k)-z(j))*(lky**2+lkx**2) +2d0*(z(k)-z(l))*(jky**2+jkx**2) -2*(2d0*z(k)-z(j)-z(l))*(lky*jky +lkx*jkx)
    bd = ((ld*sso2)+(sso1*rd))/(2d0*p1p2)
    !ra=sab/p1p2
    id = ((td*p1p2)-(bd*sab))/(p1p2**2d0)
    !       dihedrKZ=(-id)/(1d0-ra**2d0)**(0.5d0)
    dihedrkz = id!(-id)/(1d0-ra**2d0)**(0.5d0)
    !------------------Zk---------------------------
end function dihedrkz
