module nnparam
    implicit none
    private
    public :: alpha, biasa, ifunc, natom, nbasis, ndata, ndim, nhid, ninput, nlayer, nodemax, nodes, &
        & noutput, nscale, nwe, para2, pavga, pdela, pi, radian, vpesmin, weighta
    double precision, parameter   :: alpha   = 1.0d0, pi=3.1415926d0, radian=pi/180.0d0
    double precision, parameter   :: para2   = 27.2116d0
    integer, parameter            :: nbasis  = 23, ndim=6, natom=4, ndata=33435
    double precision, parameter   :: vpesmin = -399.45863953d0
    integer, parameter            :: ninput  = 22
    integer                       :: noutput, nhid, nlayer, ifunc, nwe, nodemax
    integer                       :: nscale
    integer, allocatable          :: nodes(:)
    double precision, allocatable :: weighta(:,:,:), biasa(:,:)
    double precision, allocatable :: pdela(:), pavga(:)
end module nnparam

!> @brief Implement H3S_PIPNN, a potential energy surface for H3S.
module chempotpy_h3s_pipnn_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3S_PIPNN potential energy surface for H3S.
    subroutine evaluate( x, igrad, path, p, g, d )
        use nnparam, only: alpha, biasa, ifunc, natom, nbasis, ndata, ndim, nhid, ninput, nlayer, &
            & nodemax, nodes, noutput, nscale, nwe, para2, pavga, pdela, pi, radian, vpesmin, weighta
        implicit none
        ! number of electronic state
        integer, parameter              :: nstates = 1
        integer, parameter              :: natoms  = 4
        integer, intent(in)             :: igrad
        character(len=1024), intent(in) :: path
        double precision, intent(in)    :: x(natoms,3)
        double precision, intent(out)   :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out)   :: d(nstates,nstates,natoms,3)
        double precision                :: v, va, vb, vc
        double precision                :: tx(3,natoms), tg(3,natoms)
        integer                         :: iatom, idir, j, istate
        !initialize
        v = 0.d0; g = 0.d0; d = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        call pes_init(path)
        call evvdvdx(tx, v, tg, igrad)
        deallocate(nodes)
        deallocate(weighta)
        deallocate(biasa)
        deallocate(pdela)
        deallocate(pavga)
        do istate=1, nstates
            p(istate) = v/23.0609d0
        enddo
        do iatom=1, natoms
            do idir=1, 3
                g(1,iatom,idir) = tg(idir,iatom)/23.0609d0
            enddo
        enddo
        !g=0.d0
        d = 0.d0
    end subroutine evaluate
end module chempotpy_h3s_pipnn_pes

!> @brief Expose H3S_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_h3s_pipnn_pes, only: evaluate
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

!> @brief Evaluate pipNN.
subroutine pipnn( ct, vpes, vpesa, vpesb, vpesc )
    use chempotpy_geometry, only: all_pair_distances
    use nnparam,            only: alpha, biasa, ifunc, natom, nbasis, ndata, ndim, nhid, ninput, nlayer, &
        & nodemax, nodes, noutput, nscale, nwe, para2, pavga, pdela, pi, radian, vpesmin, weighta
    implicit none
    integer          :: j, lpes
    double precision :: rb(ndim), xbond(ndim), basis(nbasis), tmp1, txinput(nbasis-1)
    double precision :: ct(3,natom), cy(3,natom)
    double precision :: m(30)
    double precision :: vpes, vpesa, vpesb, vpesc, vpesh, vnn
    double precision :: dvdg(ninput)
    integer          :: ndriv
    basis = 0.d0
    call all_pair_distances(transpose(ct), rb)
    xbond(:) = dexp(-rb(:)/alpha)
    call bemsav(xbond, m, basis)
    do j=1, nbasis-1
        txinput(j) = basis(j+1)
    enddo
    call getpota(txinput, vpesa, dvdg, ndriv)
    vnn  = vpesc
    vpes = vnn
    !-->      write(*,*)vpes
end subroutine pipnn

!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weights.txt'
!-->  biases saved in 'biases.txt'
!-->  one has to call this subroutine one and only one before calling the getpot() subroutine
!> @brief Evaluate pes_init.
subroutine pes_init( path )
    use nnparam, only: alpha, biasa, ifunc, natom, nbasis, ndata, ndim, nhid, ninput, nlayer, &
        & nodemax, nodes, noutput, nscale, nwe, para2, pavga, pdela, pi, radian, vpesmin, weighta
    implicit none
    character(len=1024), intent(in) :: path
    integer                         :: file_unit_4
    integer                         :: i, ihid, iwe, inode1, inode2, ilay1, ilay2
    integer                         :: ibasis, npd, iterm, ib, nfile
    character(len=80)               :: f1, line
    character(len=1024)             :: file_path1, file_path2
    file_path1 = trim(path)//"/H3S/weights.txt"
    file_path2 = trim(path)//"/H3S/biase.txt"
    open(newunit=file_unit_4, file=file_path2, status='old')
    read(file_unit_4, '(a80)')line
    nfile = 7
    open(nfile, file=file_path1)
    read(nfile, '(a80)')line
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
        &   biasa(nodemax,2:nlayer))
    read(nfile, *)ifunc, nwe
    !-->....ifunc hence controls the type of transfer function used for hidden layers
    !-->....At this time, only an equivalent transfer function can be used for all hidden layers
    !-->....and the pure linear function is always applid to the output layer.
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
    close(nfile)
    close(file_unit_4)
end subroutine pes_init

!> @brief Evaluate getpota.
subroutine getpota( x, vpot, dvdg, ndriv )
    use nnparam, only: alpha, biasa, ifunc, natom, nbasis, ndata, ndim, nhid, ninput, nlayer, &
        & nodemax, nodes, noutput, nscale, nwe, para2, pavga, pdela, pi, radian, vpesmin, weighta
    implicit none
    integer                       :: i, inode1, inode2, ilay1, ilay2
    integer                       :: j, k, neu1, neu2, ndriv
    double precision              :: x(ninput), y(nodemax,nlayer), vpot
    double precision, allocatable :: nxw3(:), nxw2(:,:), nxw1(:,:), ax(:), bx(:)
    double precision              :: dvdg(ninput), dvtmp
    double precision, external    :: tranfun
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
    if( ndriv .eq. 1 )then
        neu1 = nodes(2);neu2=nodes(3)
        allocate(nxw1(1:ninput,1:neu1), nxw2(1:neu1,1:neu2), &
            &ax(1:neu1), bx(1:neu2), nxw3(1:neu2))
        do i=1, ninput
            do j=1, neu1
                nxw1(i,j) = weighta(i,j,2)
            enddo
        enddo
        do j=1, neu1
            do k=1, neu2
                nxw2(j,k) = weighta(j,k,3)
            enddo
        enddo
        do j=1, neu1
            ax(j) = y(j,2)
        enddo
        do k=1, neu2
            bx(k)   = y(k,3)
            nxw3(k) = weighta(k,1,4)
        enddo
        !        stop
        dvdg  = 0d0
        dvtmp = 0d0
        do i=1, ninput
            do k=1, neu2
                dvtmp = 0d0
                do j=1, neu1
                    dvtmp = dvtmp+nxw2(j,k)*nxw1(i,j)*(1d0-ax(j)**2)
                enddo
                dvdg(i) = dvdg(i)+nxw3(k)*dvtmp*(1d0-bx(k)**2)
            enddo
            dvdg(i) = dvdg(i)*pdela(nscale)/pdela(i)
        enddo
        deallocate(nxw1, nxw2, nxw3, ax, bx)
    endif
end subroutine getpota

!> @brief Evaluate tranfun.
function tranfun( x, ifunc )
    implicit none
    integer          :: ifunc
    double precision :: tranfun, x
    if( ifunc .eq. 1 )then
        tranfun = dtanh(x)
    elseif( ifunc .eq. 2 )then
        tranfun = 1d0/(1d0+exp(-x))
    elseif( ifunc .eq. 3 )then
        tranfun = x
    endif
end

!> @brief Evaluate evvdvdx.
subroutine evvdvdx( xcart, v, dvdx, ndriv )
    use chempotpy_geometry, only: all_pair_distance_jacobian, all_pair_distances
    implicit none
    double precision, parameter :: alpha = 1d0
    integer                     :: j, ndriv
    double precision            :: v, dvdx(1:12), p(0:22), r(1:6)
    double precision            :: xcart(1:3,1:4), m(0:29), xmorse(1:6)
    double precision            :: dvdr(1:6), drdx(1:6,1:12)
    double precision            :: txinput(1:22), dvdg(1:22), dpdr(1:6,0:22)
    call all_pair_distances(transpose(xcart), r)
    xmorse = dexp(-r/alpha)
    call bemsav(xmorse, m, p)
    do j=1, 22
        txinput(j) = p(j)
    enddo
    !      write(*,*)txinput(:)
    !      stop
    call getpota(txinput, v, dvdg, ndriv)
    if( ndriv .eq. 1 )then
        call evdvdr(r, m, p, dpdr)
        call all_pair_distance_jacobian(transpose(xcart), r, drdx)
        call evdvdx(dvdg, dpdr, drdx, dvdx)
    endif
    if( ndriv .eq. 2 )then
        write(*, *) 'Only energy and gradient are available'
    endif
end subroutine evvdvdx

!> @brief Evaluate bemsav.
subroutine bemsav( xmorse, m, p )
    implicit none
    double precision, dimension(1:6), intent(in)::xmorse
    double precision, dimension(0:22), intent(out)::p
    double precision, dimension(0:29)::m
    call evmono(xmorse, m)
    call evpoly(m, p)
end subroutine bemsav

!> @brief Evaluate evmono.
subroutine evmono( x, m )
    implicit none
    double precision, dimension(1:6), intent(in)::x
    double precision, dimension(0:29), intent(out)::m
    m(0)  = 1.d0
    m(1)  = x(6)
    m(2)  = x(5)
    m(3)  = x(3)
    m(4)  = x(4)
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
end subroutine evmono

!> @brief Evaluate evpoly.
subroutine evpoly( m, p )
    implicit none
    double precision, dimension(0:29), intent(in)::m
    double precision, dimension(0:22), intent(out)::p
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
    p(13) = m(29)
    p(14) = p(1)*p(6)-p(12)
    p(15) = p(1)*p(3)-p(9)-p(9)-p(9)
    p(16) = p(1)*p(4)-p(10)
    p(17) = p(2)*p(7)-p(16)
    p(18) = p(2)*p(4)-p(12)
    p(19) = p(1)*p(8)-p(18)
    p(20) = p(2)*p(6)-p(13)-p(13)-p(13)
    p(21) = p(1)*p(7)-p(15)
    p(22) = p(2)*p(8)-p(20)
end subroutine evpoly

!> @brief Evaluate evdvdx.
subroutine evdvdx( dvdg, dpdr, drdx, dvdx )
    use chempotpy_geometry, only: transform_internal_gradient_to_cartesian
    implicit none
    double precision, intent(in)  :: dvdg(1:22), dpdr(1:6,0:22), drdx(1:6,1:12)
    double precision, intent(out) :: dvdx(1:12)
    double precision              :: dvdr(1:6)
    dvdr = matmul(dpdr(:,1:22), dvdg)
    call transform_internal_gradient_to_cartesian(dvdr, drdx, dvdx)
end subroutine evdvdx

!> @brief Evaluate evdvdr.
subroutine evdvdr( r, m, p, dpdr )
    implicit none
    integer          :: i, j
    double precision :: r(1:6), dmsdr(1:6,1:6), dmdr(1:6,0:29)
    double precision :: dvdr(1:6), dpdr(1:6,0:22)
    double precision :: m(0:29), p(0:22)
    dvdr(:) = 0d0
    call evdmsdr(r, dmsdr)
    call evdmdr(m, dmsdr, dmdr)
    call evdpdr(p, dmdr, dpdr)
end subroutine evdvdr

!> @brief Evaluate evdrdx.
subroutine evdrdx( r, x, drdx )
    use chempotpy_geometry, only: all_pair_distance_jacobian
    implicit none
    double precision :: r(1:6), drdx(1:6,1:12), x(1:12)
    call all_pair_distance_jacobian(x, r, drdx)
end subroutine evdrdx

!> @brief Evaluate evdmsdr.
subroutine evdmsdr( r, dmsdr )
    implicit none
    double precision, parameter :: alpha = 1d0
    integer                     :: i, j
    double precision            :: dmsdr(1:6,1:6), r(1:6)
    dmsdr(:,:) = 0d0
    !Morse term ms = exp(-r/alpha)
    ! dmsdr(i,j)=0  i!=j
    ! dmsdr(i,i)= -(1/alpha)*Exp(-r(i)/alpha)
    do i=1, 6
        dmsdr(i,i) = -(1d0/alpha)*dexp(-r(i)/alpha)
    enddo
end subroutine evdmsdr

!> @brief Evaluate evdmdr.
subroutine evdmdr( m, dmsdr, dmdr )
    implicit none
    integer          :: i, j
    double precision :: dmsdr(1:6,1:6), dmdr(1:6,0:29), m(0:29)
    do i=1, 6
        dmdr(i,0)  = 0.d0
        dmdr(i,1)  = dmsdr(i,6)
        dmdr(i,2)  = dmsdr(i,5)
        dmdr(i,3)  = dmsdr(i,3)
        dmdr(i,4)  = dmsdr(i,4)
        dmdr(i,5)  = dmsdr(i,2)
        dmdr(i,6)  = dmsdr(i,1)
        dmdr(i,7)  = dmdr(i,1)*m(2)+m(1)*dmdr(i,2)
        dmdr(i,8)  = dmdr(i,1)*m(3)+m(1)*dmdr(i,3)
        dmdr(i,9)  = dmdr(i,2)*m(3)+m(2)*dmdr(i,3)
        dmdr(i,10) = dmdr(i,3)*m(4)+m(3)*dmdr(i,4)
        dmdr(i,11) = dmdr(i,2)*m(5)+m(2)*dmdr(i,5)
        dmdr(i,12) = dmdr(i,1)*m(6)+m(1)*dmdr(i,6)
        dmdr(i,13) = dmdr(i,4)*m(5)+m(4)*dmdr(i,5)
        dmdr(i,14) = dmdr(i,4)*m(6)+m(4)*dmdr(i,6)
        dmdr(i,15) = dmdr(i,5)*m(6)+m(5)*dmdr(i,6)
        dmdr(i,16) = dmdr(i,1)*m(9)+m(1)*dmdr(i,9)
        dmdr(i,17) = dmdr(i,1)*m(10)+m(1)*dmdr(i,10)
        dmdr(i,18) = dmdr(i,2)*m(10)+m(2)*dmdr(i,10)
        dmdr(i,19) = dmdr(i,1)*m(11)+m(1)*dmdr(i,11)
        dmdr(i,20) = dmdr(i,3)*m(11)+m(3)*dmdr(i,11)
        dmdr(i,21) = dmdr(i,2)*m(12)+m(2)*dmdr(i,12)
        dmdr(i,22) = dmdr(i,3)*m(12)+m(3)*dmdr(i,12)
        dmdr(i,23) = dmdr(i,2)*m(13)+m(2)*dmdr(i,13)
        dmdr(i,24) = dmdr(i,3)*m(13)+m(3)*dmdr(i,13)
        dmdr(i,25) = dmdr(i,1)*m(14)+m(1)*dmdr(i,14)
        dmdr(i,26) = dmdr(i,3)*m(14)+m(3)*dmdr(i,14)
        dmdr(i,27) = dmdr(i,1)*m(15)+m(1)*dmdr(i,15)
        dmdr(i,28) = dmdr(i,2)*m(15)+m(2)*dmdr(i,15)
        dmdr(i,29) = dmdr(i,4)*m(15)+m(4)*dmdr(i,15)
    enddo
end subroutine evdmdr

!> @brief Evaluate evdpdr.
subroutine evdpdr( p, dmdr, dpdr )
    implicit none
    integer          :: i, j
    double precision :: dmdr(1:6,0:29), dpdr(1:6,0:22)
    double precision :: p(0:22)
    do i=1, 6
        dpdr(i,0)  = dmdr(i,0)
        dpdr(i,1)  = dmdr(i,1)+dmdr(i,2)+dmdr(i,3)
        dpdr(i,2)  = dmdr(i,4)+dmdr(i,5)+dmdr(i,6)
        dpdr(i,3)  = dmdr(i,7)+dmdr(i,8)+dmdr(i,9)
        dpdr(i,4)  = dmdr(i,10)+dmdr(i,11)+dmdr(i,12)
        dpdr(i,5)  = dpdr(i,1)*p(2)+p(1)*dpdr(i,2)-dpdr(i,4)
        dpdr(i,6)  = dmdr(i,13)+dmdr(i,14)+dmdr(i,15)
        dpdr(i,7)  = dpdr(i,1)*p(1)+p(1)*dpdr(i,1)-dpdr(i,3)-dpdr(i,3)
        dpdr(i,8)  = dpdr(i,2)*p(2)+p(2)*dpdr(i,2)-dpdr(i,6)-dpdr(i,6)
        dpdr(i,9)  = dmdr(i,16)
        dpdr(i,10) = dmdr(i,17)+dmdr(i,18)+dmdr(i,19)+dmdr(i,20)+dmdr(i,21) +dmdr(i,22)
        dpdr(i,11) = dpdr(i,2)*p(3)+p(2)*dpdr(i,3)-dpdr(i,10)
        dpdr(i,12) = dmdr(i,23)+dmdr(i,24)+dmdr(i,25)+dmdr(i,26)+dmdr(i,27) +dmdr(i,28)
        dpdr(i,13) = dmdr(i,29)
        dpdr(i,14) = dpdr(i,1)*p(6)+p(1)*dpdr(i,6)-dpdr(i,12)
        dpdr(i,15) = dpdr(i,1)*p(3)+p(1)*dpdr(i,3)-dpdr(i,9)-dpdr(i,9)-dpdr (i,9)
        dpdr(i,16) = dpdr(i,1)*p(4)+p(1)*dpdr(i,4)-dpdr(i,10)
        dpdr(i,17) = dpdr(i,2)*p(7)+p(2)*dpdr(i,7)-dpdr(i,16)
        dpdr(i,18) = dpdr(i,2)*p(4)+p(2)*dpdr(i,4)-dpdr(i,12)
        dpdr(i,19) = dpdr(i,1)*p(8)+p(1)*dpdr(i,8)-dpdr(i,18)
        dpdr(i,20) = dpdr(i,2)*p(6)+p(2)*dpdr(i,6)-dpdr(i,13)-dpdr(i,13)-dpdr(i,13)
        dpdr(i,21) = dpdr(i,1)*p(7)+p(1)*dpdr(i,7)-dpdr(i,15)
        dpdr(i,22) = dpdr(i,2)*p(8)+p(2)*dpdr(i,8)-dpdr(i,20)
    enddo
end subroutine evdpdr
