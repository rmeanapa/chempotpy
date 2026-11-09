!     Subroutine to generate values of the 6D PES
!     H2 + HCl --> H2 + HCl inelastic scattering reaction
!     (C)2018 Guo Group, UNM
!     2019.1.20 Qian Yao updated

module nnparam
    implicit none
    private
    public :: alpha, biasa, biasb, biasc, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, nodemax, &
        & nodes, noutput, nscale, nwe, pavga, pavgb, pavgc, pdela, pdelb, pdelc, pi, radian, weighta, &
        & weightb, weightc
    double precision, parameter   :: alpha  = 1.5d0, pi=3.1415926d0, radian=pi/180.0d0
    integer, parameter            :: nbasis = 50, ndim=6, natom=4
    integer, parameter            :: ninput = 50
    integer                       :: noutput, nhid, nlayer, ifunc, nwe, nodemax
    integer                       :: nscale
    integer, allocatable          :: nodes(:)
    double precision, allocatable :: weighta(:,:,:), biasa(:,:)
    double precision, allocatable :: pdela(:), pavga(:)
    double precision, allocatable :: weightb(:,:,:), biasb(:,:)
    double precision, allocatable :: pdelb(:), pavgb(:)
    double precision, allocatable :: weightc(:,:,:), biasc(:,:)
    double precision, allocatable :: pdelc(:), pavgc(:)
end module nnparam

!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_h3cl_pipnn_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, path, p, g, d )
        use nnparam, only: alpha, biasa, biasb, biasc, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, &
            & nodemax, nodes, noutput, nscale, nwe, pavga, pavgb, pavgc, pdela, pdelb, pdelc, pi, radian, &
            & weighta, weightb, weightc
        implicit none
        ! number of electronic state
        integer, parameter              :: nstates = 1
        integer, parameter              :: natoms  = 4
        integer, intent(in)             :: igrad
        character(len=1024), intent(in) :: path
        double precision, intent(in)    :: x(natoms,3)
        double precision, intent(out)   :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out)   :: d(nstates,nstates,natoms,3)
        double precision                :: v
        double precision                :: tx(3,natoms)
        integer                         :: iatom, idir, j, istate
        !initialize
        v = 0.d0; g = 0.d0; d = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        if( igrad == 0 )then
            call pes_init(path)
            call pesh3cl(tx, v, path)
            deallocate(nodes)
            deallocate(weighta)
            deallocate(biasa)
            deallocate(pdela)
            deallocate(pavga)
            deallocate(weightb)
            deallocate(biasb)
            deallocate(pdelb)
            deallocate(pavgb)
            deallocate(weightc)
            deallocate(biasc)
            deallocate(pdelc)
            deallocate(pavgc)
        else
            write (*, *) 'Only energy and gradient are available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
        g = 0.d0
        d = 0.d0
    end subroutine evaluate
end module chempotpy_h3cl_pipnn_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_h3cl_pipnn_pes, only: evaluate
    implicit none
    integer, parameter              :: natoms  = 4
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    double precision, intent(in)    :: x(natoms,3)
    double precision, intent(out)   :: p(nstates)
    double precision, intent(out)   :: g(nstates,natoms,3)
    double precision, intent(out)   :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, path, p, g, d)
end subroutine pes

!> @brief Evaluate pesh3cl for the H3Cl_PIPNN potential.
subroutine pesh3cl( x, vtot, path )
    implicit none
    character(len=1024), intent(in) :: path
    double precision                :: phi, r, r1, r2, r3, s, th1, th2, v, v0, v1, vint, vpes, vpesa, vpesb, vpesc, &
        & vpip, vtot, x
    dimension r(6), x(3,4)
    double precision, parameter :: tocm = 219474.63137d0
    integer, parameter          :: n    = 644
    call xyz_to_jaco(x, r1, r2, r3, th1, th2, phi, r)
    call clh3pipnn(r, vpes, vpesa, vpesb, vpesc)
    call pot0(r1, v0, path)
    call pot1(r2, v1, path)
    call clh3vlr(r1, r2, r3, th1, th2, phi, vint, n, path)
    vpip = (vpes*0.03674931d0+-461.54119162d0)-v0-v1
    vpip = vpip*tocm
    s    = (1d0-dtanh(3d0*(r3-5.5d0)))/2d0
    v    = s*vpip+(1d0-s)*vint
    vtot = (v/tocm+v0+v1--461.54020987d0)*27.21138d0
end subroutine pesh3cl

!> @brief Evaluate POT0 for the H3Cl_PIPNN potential.
subroutine pot0( r0, vv0, path )
    implicit none
    integer, parameter              :: m0 = 1425
    character(len=1024), intent(in) :: path
    integer                         :: file_unit_100
    double precision                :: dy1, dyn, r0, vv0
    double precision                :: xt0(m0), yl0(m0), ym0(m0)
    integer                         :: i, init
    character(len=1024)             :: file_path1
    dy1        = 1.0d30; dyn = 1.0d30
    init       = 0
    file_path1 = trim(path)//"/H3Cl/H2.dat"
    if( init .eq. 0 )then
        open(newunit=file_unit_100, file=file_path1, status='old')
        do i=1, m0
            read(file_unit_100, *)xt0(i), ym0(i)
        enddo
        close(file_unit_100)
        call spline(xt0, ym0, m0, dy1, dyn, yl0)
        init = 1
    endif
    call splint(xt0, ym0, yl0, m0, r0, vv0)
end subroutine pot0

!> @brief Evaluate POT1 for the H3Cl_PIPNN potential.
subroutine pot1( r1, vv1, path )
    implicit none
    integer, parameter              :: m1 = 3200
    character(len=1024), intent(in) :: path
    integer                         :: file_unit_200
    double precision                :: dy1, dyn, r1, vv1
    double precision                :: xt1(m1), yl1(m1), ym1(m1)
    integer                         :: i, init
    character(len=1024)             :: file_path1
    dy1        = 1.0d30; dyn = 1.0d30
    init       = 0
    file_path1 = trim(path)//"/H3Cl/HCl.dat"
    if( init .eq. 0 )then
        open(newunit=file_unit_200, file=file_path1, status='old')
        do i=1, m1
            read(file_unit_200, *)xt1(i), ym1(i)
        enddo
        close(file_unit_200)
        call spline(xt1, ym1, m1, dy1, dyn, yl1)
        init = 1
    endif
    call splint(xt1, ym1, yl1, m1, r1, vv1)
end subroutine pot1

!##################################################################
!# SPLINE ROUTINES
!#            Numerical recipes in fortran
!#            Cambrige University Press
!#            York, 2nd edition, 1992.
!##################################################################
!> @brief Evaluate splint for the H3Cl_PIPNN potential.
subroutine splint( xa, ya, y2a, n, x, y )
    implicit none
    double precision :: a, b, h, x, xa, y, y2a, ya
    integer          :: k, khi, klo, n
    dimension xa(n), y2a(n), ya(n)
    klo = 1
    khi = n
    do while( khi-klo .gt. 1 )
        k = (khi+klo)/2
        if( xa(k) .gt. x )then
            khi = k
        else
            klo = k
        endif
    enddo
    h = xa(khi)-xa(klo)
    if( h .eq. 0.0d0 ) write(6, *) 'bad xa input in splint'
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h** 2)/6.0d0
end subroutine splint

!#######################################################################
!> @brief Evaluate spline for the H3Cl_PIPNN potential.
subroutine spline( x, y, n, yp1, ypn, y2 )
    implicit none
    double precision :: p, qn, sig, u, un, x, y, y2, yp1, ypn
    integer          :: i, k, n
    dimension x(n), y(n), y2(n)
    integer, parameter :: nmax = 4200
    dimension u(nmax)
    if( yp1 .gt. .99d30 )then
        y2(1) = 0.0d0
        u(1)  = 0.0d0
    else
        y2(1) = -0.5d0
        u(1)  = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2, n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1)+2.0d0
        y2(i) = (sig-1.0d0)/p
        u(i)  = (6.0d0*((y(i+1)-y(i))/(x(i+ 1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* u(i-1))/p
    enddo
    if( ypn .gt. .99d30 )then
        qn = 0.0d0
        un = 0.0d0
    else
        qn = 0.5d0
        un = (3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
    do k=n-1, 1, -1
        y2(k) = y2(k)*y2(k+1)+u(k)
    enddo
end subroutine spline

!> @brief Evaluate ClH3vlr for the H3Cl_PIPNN potential.
subroutine clh3vlr( r1, r2, r3, theta1, theta2, phi, vint, n, path )
    implicit none
    integer, parameter              :: m = 644
    character(len=1024), intent(in) :: path
    integer                         :: file_unit_300
    double precision                :: phi, r1, r2, r3, theta1, theta2, vint
    double precision                :: a(m)
    integer                         :: i, init, n
    character(len=1024)             :: file_path1
    double precision                :: afunc(n)
    init       = 0
    file_path1 = trim(path)//"/H3Cl/c.dat"
    if( init .eq. 0 )then
        open(newunit=file_unit_300, file=file_path1, status='old')
        read(file_unit_300, *)a
        close(file_unit_300)
        init = 1
    endif
    call funcs(r1, r2, r3, theta1, theta2, phi, afunc, n)
    vint = 0d0
    do i=1, n
        vint = vint+a(i)*afunc(i)
    enddo
end subroutine clh3vlr

!> @brief Evaluate FUNCS for the H3Cl_PIPNN potential.
subroutine funcs( ra, rb, rc, theta1, theta2, phi, afunc, ma )
    implicit none
    integer, parameter :: m = 3, m1=2*m+1, m2=m1**2! m the order of r1 and r2
    double precision   :: z1(m2), z2(m2), z3(m2), z4(m2), z5(m1), z6(m2), z7(m2), &
        & afunc(ma), a(m2), z8(m2), z9(m2), z10(m2), z11(m2), z12(m2), z13(m2) &
        &, z14(m2), a111(m2)
    integer          :: ma, i, l, n, j, j1
    double precision :: ra, rb, rc, theta1, theta2, phi, r1, r2, r3, th1, th2, a1, a2, a3 &
        &, a11, a22
    r2 = ra;r1=rb;r3=rc;th1=theta2;th2=theta1
    n  = 1
    j  = 0
    do i=-m, m
        a1 = r1**i
        do l=-m, m
            a2 = (r2**l)*a1
            if( n .gt. 0 )then
                j    = j+1
                a(j) = a2
            endif
        enddo;enddo
        j1 = 0
        do i=-m, m
            a11 = r1**i
            do l=-m, m
                a22 = 3d0*(r2**l)*a11/(2d0*(r2**l+a11))
                if( n .gt. 0 )then
                    j1       = j1+1
                    a111(j1) = a22
                endif
            enddo;enddo
            do i=1, m2
                z1(i)=a(i)*(1.5d0*(dcos(th1)*(3d0*(dcos(th2)**2)-1d0) +2d0*dsin(th1)*dsin(th2)*dcos(th2)*dcos(phi)))/ &
                    & (r3**4)
                z2(i)=a(i)*(0.75d0*(1d0-5d0*(dcos(th1)**2)- 5d0*(dcos(th2)**2)+17d0*(dcos(th1)**2)*(dcos(th2)**2)+ &
                    & 2d0* (dsin(th1)**2)*(dsin(th2)**2)*(dcos(phi)**2)+ 16d0*dsin(th1)*dcos(th1)*dsin(th2)*dcos(th2)* &
                    & dcos(phi)))/(r3**5)
                z3(i) = a(i)*(-0.5d0*(3d0*(dcos(th1)**2)+1d0))/(r3**6)
                z4(i)=a(i)*((12d0*(dcos(th1)**2)*(dcos(th2)**2)+3d0* (dsin(th1)**2)*(dsin(th2)**2)*(dcos(phi)**2)-3d0* &
                    & (cos(th1)**2)-1d0 +12d0*dsin(th1)*dcos(th1)*dsin(th2)*dcos(th2)*dcos(phi)) /-6d0)/(r3**6)
                z6(i)=a(i)*(-1.5d0*(6d0*(dcos(th1)**2)*(dcos(th2)**3) -2d0*(dcos(th1)**2)*dcos(th2)+7d0*dsin(th1)* &
                    & dcos(th1)*dsin(th2)* (dcos(th2)**2)*dcos(phi)-dsin(th1)*dcos(th1)*dsin(th2)*dcos(phi) +2d0* &
                    & (dsin(th1)**2)*(dsin(th2)**2)*dcos(th2)*(dcos(phi)**2))) /(r3**7)
                z7(i)=a(i)*(2d0*(6d0*(dcos(th1)**2)*(dcos(th2)**3) -5d0*(dcos(th1)**2)*dcos(th2)+7d0*dsin(th1)* &
                    & dcos(th1)*dsin(th2)* (dcos(th2)**2)*dcos(phi)-2d0*dsin(th1)*dcos(th1)*dsin(th2) *dcos(phi)+2d0* &
                    & (dsin(th1)**2)*(dsin(th2)**2)*dcos(th2) *(dcos(phi)**2)-dcos(th2)))/(r3**7)
                z8(i)             = -a111(i)/(r3**6)
                z9(i)             = (-a111(i)*(1.5d0*(dcos(th2)**2)-0.5d0)/3d0)/(r3**6)
                z10(i)            = (-a111(i)*(1.5d0*(dcos(th1)**2)-0.5d0)/3d0)/(r3**6)
                z11(i)            = -a111(i)*2d0*(dcos(th2)**3)/(r3**7)
                z12(i)            = -a111(i)*2d0*(dcos(th1)**3)/(r3**7)
                z13(i)            = -a111(i)*4d0*(3d0*dcos(th2)-2d0*(dcos(th2)**3))/ (3d0*(r3**7))
                z14(i)            = -a111(i)*4d0*(3d0*dcos(th1)-2d0*(dcos(th1)**3))/ (3d0*(r3**7))
                afunc(i)          = z1(i);afunc(m2+i)=z2(i);afunc(2*m2+i)=z3(i)
                afunc(3*m2+i)     = z4(i);afunc(4*m2+m1+i)=z6(i)
                afunc(5*m2+m1+i)  = z7(i);afunc(6*m2+m1+i)=z8(i)
                afunc(7*m2+m1+i)  = z9(i);afunc(8*m2+m1+i)=z10(i)
                afunc(9*m2+m1+i)  = z11(i);afunc(10*m2+m1+i)=z12(i)
                afunc(11*m2+m1+i) = z13(i);afunc(12*m2+m1+i)=z14(i)
            enddo
            j = 4*m2
            do i=-m, m
                a3 = (r1**i)*(-6d0*(dcos(th1)**3))/(r3**7)
                if( n .gt. 0 )then
                    j        = j+1
                    afunc(j) = a3
                endif
            enddo
end subroutine funcs

!> @brief Evaluate ClH3pipNN for the H3Cl_PIPNN potential.
subroutine clh3pipnn( rb, vpes, vpesa, vpesb, vpesc )
    use nnparam, only: alpha, biasa, biasb, biasc, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, &
        & nodemax, nodes, noutput, nscale, nwe, pavga, pavgb, pavgc, pdela, pdelb, pdelc, pi, radian, &
        & weighta, weightb, weightc
    implicit none
    integer          :: i, j, k, lpes
    double precision :: rb(ndim), xbond(ndim), basis(0:nbasis), tmp1, txinput(nbasis)
    double precision :: vpes, vpesa, vpesb, vpesc, vpesh, vnn
    basis    = 0.d0
    xbond(:) = dexp(-rb(:)/alpha)
    call bemsav134(xbond, basis)
    do j=1, nbasis
        txinput(j) = basis(j)
    enddo
    call getpota(txinput, vpesa)
    call getpotb(txinput, vpesb)
    call getpotc(txinput, vpesc)
    vpes = (vpesa+vpesb+vpesc)/3d0
end subroutine clh3pipnn

!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weights.txt'
!-->  biases saved in 'biases.txt'
!-->  one has to call this subroutine one and only one before calling
!the getpot() subroutine
!> @brief Evaluate pes_init for the H3Cl_PIPNN potential.
subroutine pes_init( path )
    use nnparam, only: alpha, biasa, biasb, biasc, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, &
        & nodemax, nodes, noutput, nscale, nwe, pavga, pavgb, pavgc, pdela, pdelb, pdelc, pi, radian, &
        & weighta, weightb, weightc
    implicit none
    character(len=1024), intent(in) :: path
    integer                         :: file_unit_4
    integer                         :: i, ihid, iwe, inode1, inode2, ilay1, ilay2
    integer                         :: ibasis, npd, iterm, ib, nfile
    character(len=80)               :: f1, line
    character(len=1024)             :: file_path1, file_path2
    character(len=1024)             :: file_path3, file_path4
    character(len=1024)             :: file_path5, file_path6
    file_path1 = trim(path)//"/H3Cl/weights.txt-1"
    file_path2 = trim(path)//"/H3Cl/biases.txt-1"
    open(newunit=file_unit_4, file=file_path2, status='old')
    nfile = 7
    open(nfile, file=file_path1)
    read(nfile, *)nhid, noutput
    nscale = ninput+noutput
    nlayer = nhid+2 !additional one for input layer and one for output
    allocate(nodes(nlayer), pdela(nscale), pavga(nscale))
    nodes(1)      = ninput
    nodes(nlayer) = noutput
    read(nfile, *)(nodes(ihid), ihid=2, nhid+1)
    nodemax = 0
    do i=1, nlayer
        nodemax = max(nodemax, nodes(i))
    enddo
    allocate(weighta(nodemax,nodemax,2:nlayer), &
        & biasa(nodemax,2:nlayer))
    read(nfile, *)ifunc, nwe
    !-->....ifunc hence controls the type of transfer function used for
    !hidden layers
    !-->....At this time, only an equivalent transfer function can be used
    !for all hidden layers
    !-->....and the pure linear function is always applid to the output
    !layer.
    !-->....see function tranfun() for details
    read(nfile, *)(pdela(i), i=1, nscale)
    read(nfile, *)(pavga(i), i=1, nscale)
    iwe = 0
    do ilay1=2, nlayer
        ilay2 = ilay1-1
        do inode1=1, nodes(ilay1)
            do inode2=1, nodes(ilay2) !
                read(nfile, *)weighta(inode2,inode1,ilay1)
                iwe = iwe+1
            enddo
            read(file_unit_4, *)biasa(inode1,ilay1)
            iwe = iwe+1
        enddo
    enddo
    if( iwe .ne. nwe )then
        write(*, *)'provided number of parameters ', nwe
        write(*, *)'actual number of parameters ', iwe
        write(*, *)'nwe not equal to iwe, check input files or code'
        stop
    endif
    close(file_unit_4);close(nfile)
    file_path3 = trim(path)//"/H3Cl/weights.txt-2"
    file_path4 = trim(path)//"/H3Cl/biases.txt-2"
    open(newunit=file_unit_4, file=file_path4, status='old')
    nfile = 8
    open(nfile, file=file_path3)
    read(nfile, *)nhid, noutput
    nscale = ninput+noutput
    nlayer = nhid+2 !additional one for input layer and one for output
    allocate(pdelb(nscale), pavgb(nscale))
    nodes(1)      = ninput
    nodes(nlayer) = noutput
    read(nfile, *)(nodes(ihid), ihid=2, nhid+1)
    nodemax = 0
    do i=1, nlayer
        nodemax = max(nodemax, nodes(i))
    enddo
    allocate(weightb(nodemax,nodemax,2:nlayer), &
        & biasb(nodemax,2:nlayer))
    read(nfile, *)ifunc, nwe
    !-->....ifunc hence controls the type of transfer function used for
    !hidden layers
    !-->....At this time, only an equivalent transfer function can be used
    !for all hidden layers
    !-->....and the pure linear function is always applid to the output
    !layer.
    !-->....see function tranfun() for details
    read(nfile, *)(pdelb(i), i=1, nscale)
    read(nfile, *)(pavgb(i), i=1, nscale)
    iwe = 0
    do ilay1=2, nlayer
        ilay2 = ilay1-1
        do inode1=1, nodes(ilay1)
            do inode2=1, nodes(ilay2) !
                read(nfile, *)weightb(inode2,inode1,ilay1)
                iwe = iwe+1
            enddo
            read(file_unit_4, *)biasb(inode1,ilay1)
            iwe = iwe+1
        enddo
    enddo
    if( iwe .ne. nwe )then
        write(*, *)'provided number of parameters ', nwe
        write(*, *)'actual number of parameters ', iwe
        write(*, *)'nwe not equal to iwe, check input files or code'
        stop
    endif
    close(file_unit_4);close(nfile)
    file_path5 = trim(path)//"/H3Cl/weights.txt-3"
    file_path6 = trim(path)//"/H3Cl/biases.txt-3"
    open(newunit=file_unit_4, file=file_path6, status='old')
    nfile = 9
    open(nfile, file=file_path5)
    read(nfile, *)nhid, noutput
    nscale = ninput+noutput
    nlayer = nhid+2 !additional one for input layer and one for output
    allocate(pdelc(nscale), pavgc(nscale))
    nodes(1)      = ninput
    nodes(nlayer) = noutput
    read(nfile, *)(nodes(ihid), ihid=2, nhid+1)
    nodemax = 0
    do i=1, nlayer
        nodemax = max(nodemax, nodes(i))
    enddo
    allocate(weightc(nodemax,nodemax,2:nlayer), &
        & biasc(nodemax,2:nlayer))
    read(nfile, *)ifunc, nwe
    !-->....ifunc hence controls the type of transfer function used for
    !hidden layers
    !-->....At this time, only an equivalent transfer function can be used
    !for all hidden layers
    !-->....and the pure linear function is always applid to the output
    !layer.
    !-->....see function tranfun() for details
    read(nfile, *)(pdelc(i), i=1, nscale)
    read(nfile, *)(pavgc(i), i=1, nscale)
    iwe = 0
    do ilay1=2, nlayer
        ilay2 = ilay1-1
        do inode1=1, nodes(ilay1)
            do inode2=1, nodes(ilay2)
                read(nfile, *)weightc(inode2,inode1,ilay1)
                iwe = iwe+1
            enddo
            read(file_unit_4, *)biasc(inode1,ilay1)
            iwe = iwe+1
        enddo
    enddo
    if( iwe .ne. nwe )then
        write(*, *)'provided number of parameters ', nwe
        write(*, *)'actual number of parameters ', iwe
        write(*, *)'nwe not equal to iwe, check input files or code'
        stop
    endif
    close(nfile)
    close(file_unit_4)
end subroutine pes_init

!> @brief Evaluate tranfun for the H3Cl_PIPNN potential.
function tranfun( x, ifunc )
    implicit none
    integer          :: ifunc
    double precision :: tranfun, x
    !    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
    !    ifunc=2, transfer function is log sigmoid function, 'logsig'
    !    ifunc=3, transfer function is pure linear function, 'purelin'. It
    !is imposed to the output layer by default
    if( ifunc .eq. 1 )then
        tranfun = dtanh(x)
    else if( ifunc .eq. 2 )then
        tranfun = 1d0/(1d0+dexp(-x))
    else if( ifunc .eq. 3 )then
        tranfun = x
    endif
end function tranfun

!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->subroutine to accelerate~~~~~~

!> @brief Evaluate emsav134 for the H3Cl_PIPNN potential.
function emsav134( x, c ) result(v )
    implicit none
    double precision, dimension(1:6)::x
    double precision, dimension(0:50)::c
    double precision :: v
    ! ::::::::::::::::::::
    double precision, dimension(0:50)::p
    call bemsav134(x, p)
    v = dot_product(p, c)
end function emsav134

!> @brief Evaluate bemsav134 for the H3Cl_PIPNN potential.
subroutine bemsav134( x, p )
    implicit none
    double precision, dimension(1:6), intent(in)::x
    double precision, dimension(0:50), intent(out)::p
    ! ::::::::::::::::::::
    double precision, dimension(0:32)::m
    call evmono134(x, m)
    call evpoly134(m, p)
end subroutine bemsav134

!> @brief Evaluate evmono134 for the H3Cl_PIPNN potential.
subroutine evmono134( x, m )
    implicit none
    double precision, dimension(1:6), intent(in)::x
    double precision, dimension(0:32), intent(out)::m
    m(0)  = 1.d0
    m(1)  = x(6)
    m(2)  = x(5)
    m(3)  = x(4)
    m(4)  = x(3)
    m(5)  = x(2)
    m(6)  = x(1)
    m(7)  = m(1)*m(2)
    m(8)  = m(1)*m(3)
    m(9)  = m(2)*m(3)
    m(10) = m(3)*m(4)
    m(11) = m(2)*m(5)
    m(12) = m(1)*m(6)
    m(13) = m(4)*m(5)
    m(14) = m(4)*m(6)
    m(15) = m(5)*m(6)
    m(16) = m(1)*m(9)
    m(17) = m(1)*m(10)
    m(18) = m(2)*m(10)
    m(19) = m(1)*m(11)
    m(20) = m(3)*m(11)
    m(21) = m(2)*m(12)
    m(22) = m(3)*m(12)
    m(23) = m(2)*m(13)
    m(24) = m(3)*m(13)
    m(25) = m(1)*m(14)
    m(26) = m(3)*m(14)
    m(27) = m(1)*m(15)
    m(28) = m(2)*m(15)
    m(29) = m(4)*m(15)
    m(30) = m(2)*m(24)
    m(31) = m(1)*m(26)
    m(32) = m(1)*m(28)
end subroutine evmono134

!> @brief Evaluate evpoly134 for the H3Cl_PIPNN potential.
subroutine evpoly134( m, p )
    implicit none
    double precision, dimension(0:32), intent(in)::m
    double precision, dimension(0:50), intent(out)::p
    p(0)  = m(0)
    p(1)  = m(1)+m(2)+m(3)
    p(2)  = m(4)+m(5)+m(6)
    p(3)  = m(7)+m(8)+m(9)
    p(4)  = m(10)+m(11)+m(12)
    p(5)  = p(1)*p(2)-p(4)
    p(6)  = m(13)+m(14)+m(15)
    p(7)  = p(1)*p(1)-p(3)-p(3)
    p(8)  = p(2)*p(2)-p(6)-p(6)
    p(9)  = m(16)
    p(10) = m(17)+m(18)+m(19)+m(20)+m(21)+m(22)
    p(11) = p(2)*p(3)-p(10)
    p(12) = m(23)+m(24)+m(25)+m(26)+m(27)+m(28)
    p(13) = p(1)*p(6)-p(12)
    p(14) = m(29)
    p(15) = p(1)*p(3)-p(9)-p(9)-p(9)
    p(16) = p(1)*p(4)-p(10)
    p(17) = p(2)*p(7)-p(16)
    p(18) = p(2)*p(4)-p(12)
    p(19) = p(1)*p(8)-p(18)
    p(20) = p(2)*p(6)-p(14)-p(14)-p(14)
    p(21) = p(1)*p(7)-p(15)
    p(22) = p(2)*p(8)-p(20)
    p(23) = p(9)*p(2)
    p(24) = m(30)+m(31)+m(32)
    p(25) = p(3)*p(6)-p(24)
    p(26) = p(14)*p(1)
    p(27) = p(9)*p(1)
    p(28) = p(3)*p(4)-p(23)
    p(29) = p(1)*p(10)-p(23)-p(28)-p(23)
    p(30) = p(1)*p(11)-p(23)
    p(31) = p(1)*p(12)-p(25)-p(24)-p(24)
    p(32) = p(1)*p(13)-p(25)
    p(33) = p(4)*p(5)-p(25)-p(31)
    p(34) = p(2)*p(11)-p(25)
    p(35) = p(4)*p(6)-p(26)
    p(36) = p(2)*p(12)-p(26)-p(35)-p(26)
    p(37) = p(2)*p(13)-p(26)
    p(38) = p(14)*p(2)
    p(39) = p(3)*p(3)-p(27)-p(27)
    p(40) = p(3)*p(7)-p(27)
    p(41) = p(1)*p(16)-p(28)
    p(42) = p(2)*p(21)-p(41)
    p(43) = p(1)*p(18)-p(33)
    p(44) = p(7)*p(8)-p(43)
    p(45) = p(6)*p(6)-p(38)-p(38)
    p(46) = p(2)*p(18)-p(35)
    p(47) = p(1)*p(22)-p(46)
    p(48) = p(6)*p(8)-p(38)
    p(49) = p(1)*p(21)-p(40)
    p(50) = p(2)*p(22)-p(48)
end subroutine evpoly134

!> @brief Evaluate getpota for the H3Cl_PIPNN potential.
subroutine getpota( x, vpot )
    use nnparam, only: alpha, biasa, biasb, biasc, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, &
        & nodemax, nodes, noutput, nscale, nwe, pavga, pavgb, pavgc, pdela, pdelb, pdelc, pi, radian, &
        & weighta, weightb, weightc
    implicit none
    integer                    :: i, ndriv, inode1, inode2, ilay1, ilay2, j, k, neu1, neu2
    double precision           :: x(ninput), y(nodemax,nlayer), vpot
    double precision, external :: tranfun
    !-->....set up the normalized input layer
    do i=1, ninput
        y(i,1) = (x(i)-pavga(i))/pdela(i)
    enddo
    !-->....evaluate the hidden layer
    do ilay1=2, nlayer-1
        ilay2 = ilay1-1
        do inode1=1, nodes(ilay1)
            y(inode1,ilay1) = biasa(inode1,ilay1)
            do inode2=1, nodes(ilay2)
                y(inode1,ilay1) = y(inode1,ilay1)+y(inode2,ilay2) *weighta(inode2,inode1,ilay1)
            enddo
            y(inode1,ilay1) = tranfun(y(inode1,ilay1), ifunc)
        enddo
    enddo
    !-->....now evaluate the output
    ilay1 = nlayer
    ilay2 = ilay1-1
    do inode1=1, nodes(ilay1)
        y(inode1,ilay1) = biasa(inode1,ilay1)
        do inode2=1, nodes(ilay2)
            y(inode1,ilay1) = y(inode1,ilay1)+y(inode2,ilay2) *weighta(inode2,inode1,ilay1)
        enddo
        !-->....the transfer function is linear y=x for output layer
        !-->....so no operation is needed here
    enddo
    !-->....the value of output layer is the fitted potntial
    vpot = y(nodes(nlayer),nlayer)*pdela(nscale)+pavga(nscale)
end subroutine getpota

!> @brief Evaluate getpotb for the H3Cl_PIPNN potential.
subroutine getpotb( x, vpot )
    use nnparam, only: alpha, biasa, biasb, biasc, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, &
        & nodemax, nodes, noutput, nscale, nwe, pavga, pavgb, pavgc, pdela, pdelb, pdelc, pi, radian, &
        & weighta, weightb, weightc
    implicit none
    integer                    :: i, ndriv, inode1, inode2, ilay1, ilay2, j, k, neu1, neu2
    double precision           :: x(ninput), y(nodemax,nlayer), vpot
    double precision, external :: tranfun
    !-->....set up the normalized input layer
    do i=1, ninput
        y(i,1) = (x(i)-pavgb(i))/pdelb(i)
    enddo
    !-->....evaluate the hidden layer
    do ilay1=2, nlayer-1
        ilay2 = ilay1-1
        do inode1=1, nodes(ilay1)
            y(inode1,ilay1) = biasb(inode1,ilay1)
            do inode2=1, nodes(ilay2)
                y(inode1,ilay1) = y(inode1,ilay1)+y(inode2,ilay2) *weightb(inode2,inode1,ilay1)
            enddo
            y(inode1,ilay1) = tranfun(y(inode1,ilay1), ifunc)
        enddo
    enddo
    !-->....now evaluate the output
    ilay1 = nlayer
    ilay2 = ilay1-1
    do inode1=1, nodes(ilay1)
        y(inode1,ilay1) = biasb(inode1,ilay1)
        do inode2=1, nodes(ilay2)
            y(inode1,ilay1) = y(inode1,ilay1)+y(inode2,ilay2) *weightb(inode2,inode1,ilay1)
        enddo
        !-->....the transfer function is linear y=x for output layer
        !-->....so no operation is needed here
    enddo
    !-->....the value of output layer is the fitted potntial
    vpot = y(nodes(nlayer),nlayer)*pdelb(nscale)+pavgb(nscale)
end subroutine getpotb

!> @brief Evaluate getpotc for the H3Cl_PIPNN potential.
subroutine getpotc( x, vpot )
    use nnparam, only: alpha, biasa, biasb, biasc, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, &
        & nodemax, nodes, noutput, nscale, nwe, pavga, pavgb, pavgc, pdela, pdelb, pdelc, pi, radian, &
        & weighta, weightb, weightc
    implicit none
    integer                    :: i, ndriv, inode1, inode2, ilay1, ilay2, j, k, neu1, neu2
    double precision           :: x(ninput), y(nodemax,nlayer), vpot
    double precision, external :: tranfun
    !-->....set up the normalized input layer
    do i=1, ninput
        y(i,1) = (x(i)-pavgc(i))/pdelc(i)
    enddo
    !-->....evaluate the hidden layer
    do ilay1=2, nlayer-1
        ilay2 = ilay1-1
        do inode1=1, nodes(ilay1)
            y(inode1,ilay1) = biasc(inode1,ilay1)
            do inode2=1, nodes(ilay2)
                y(inode1,ilay1) = y(inode1,ilay1)+y(inode2,ilay2) *weightc(inode2,inode1,ilay1)
            enddo
            y(inode1,ilay1) = tranfun(y(inode1,ilay1), ifunc)
        enddo
    enddo
    !-->....now evaluate the output
    ilay1 = nlayer
    ilay2 = ilay1-1
    do inode1=1, nodes(ilay1)
        y(inode1,ilay1) = biasc(inode1,ilay1)
        do inode2=1, nodes(ilay2)
            y(inode1,ilay1) = y(inode1,ilay1)+y(inode2,ilay2) *weightc(inode2,inode1,ilay1)
        enddo
        !-->....the transfer function is linear y=x for output layer
        !-->....so no operation is needed here
    enddo
    !-->....the value of output layer is the fitted potntial
    vpot = y(nodes(nlayer),nlayer)*pdelc(nscale)+pavgc(nscale)
end subroutine getpotc

!> @brief Evaluate distance for the H3Cl_PIPNN potential.
double precision function distance( a, b )
    implicit none
    double precision :: a, b
    integer          :: i
    dimension a(3), b(3)
    distance = 0.d0
    do i=1, 3
        distance = distance+(a(i)-b(i))**2
    enddo
    distance = dsqrt(distance)
end function distance

!> @brief Evaluate xyz_to_jaco for the H3Cl_PIPNN potential.
subroutine xyz_to_jaco( x, r1, r2, r3, theta1, theta2, phi, r )
    implicit none
    double precision            :: distance, ra1, ra4, rb1, rb4, str11, str12
    double precision            :: x(3,4), a(3), b(3), r(6)
    double precision            :: mass(4)
    double precision            :: r1, r2, r3, theta1, theta2, phi
    double precision, parameter :: pi = 3.14159265d0
    integer                     :: k, m, i, j, n, ndata
    character                   :: long
    mass(1) = 1.00782503d0
    mass(2) = mass(1)
    mass(3) = mass(1)
    mass(4) = 35.4532d0
    r1      = distance(x(:,1), x(:,2))
    r2      = distance(x(:,3), x(:,4))
    a(:)    = x(:,1)+(mass(2)/(mass(1)+mass(2)))*(x(:,2)-x(:,1))
    b(:)    = x(:,4)+(mass(3)/(mass(4)+mass(3)))*(x(:,3)-x(:,4))
    r3      = distance(a(:), b(:))
    ra4     = distance(a(:), x(:,4))
    rb4     = distance(b(:), x(:,4))
    str11   = (r3**2+rb4**2-ra4**2)/(2d0*r3*rb4)
    str11   = min(1d0, str11);str11=max(-1d0, str11)
    theta2  = dacos(str11)
    ra1     = distance(a(:), x(:,1))
    rb1     = distance(b(:), x(:,1))
    str12   = (r3**2+ra1**2-rb1**2)/(2d0*r3*ra1)
    str12   = min(1d0, str12);str12=max(-1d0, str12)
    theta1  = dacos(str12)
    call dihedral(x(:,2), a, b, x(:,4), phi)
    phi = phi*pi/180d0
    if( isnan(phi) )then
        phi = 0d0
    endif
    r(1) = distance(x(:,4), x(:,3))
    r(2) = distance(x(:,4), x(:,2))
    r(3) = distance(x(:,4), x(:,1))
    r(4) = distance(x(:,3), x(:,2))
    r(5) = distance(x(:,3), x(:,1))
    r(6) = distance(x(:,2), x(:,1))
end subroutine xyz_to_jaco

!> @brief Evaluate dihedral for the H3Cl_PIPNN potential.
subroutine dihedral( a, b, c, d, th )
    implicit none
    double precision :: a(3), b(3), c(3), d(3), th, x, y
    double precision :: v_ab(3), v_bc(3), v_cd(3), n1(3), n2(3), s1, s2, s3, m1(3)
    integer          :: i
    do i=1, 3
        v_ab(i) = b(i)-a(i)
        v_bc(i) = c(i)-b(i)
        v_cd(i) = d(i)-c(i)
    enddo
    n1(1) = v_ab(2)*v_bc(3)-v_ab(3)*v_bc(2)
    n1(2) = -v_ab(1)*v_bc(3)+v_ab(3)*v_bc(1)
    n1(3) = v_ab(1)*v_bc(2)-v_ab(2)*v_bc(1)
    n2(1) = v_bc(2)*v_cd(3)-v_bc(3)*v_cd(2)
    n2(2) = -v_bc(1)*v_cd(3)+v_bc(3)*v_cd(1)
    n2(3) = v_bc(1)*v_cd(2)-v_bc(2)*v_cd(1)
    m1(1) = n1(2)*v_bc(3)-n1(3)*v_bc(2)
    m1(2) = -n1(1)*v_bc(3)+n1(3)*v_bc(1)
    m1(3) = n1(1)*v_bc(2)-n1(2)*v_bc(1)
    s1    = 0d0;s2=0d0;s3=0d0
    do i=1, 3
        s1 = s1+n1(i)**2
        s2 = s2+n2(i)**2
        s3 = s3+m1(i)**2
    enddo
    n1 = n1/dsqrt(s1);n2=n2/dsqrt(s2);m1=m1/dsqrt(s3)
    x  = 0d0;y=0d0
    do i=1, 3
        x = x+n1(i)*n2(i)
        y = y+m1(i)*n2(i)
    enddo
    th = -datan2(y, x)
end subroutine dihedral
