!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
module nnparam
    implicit none
    private
    public :: alpha, biasa, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, nodemax, nodes, &
        & noutput, nscale, nwe, pavga, pdela, pi, radian, vpesmin, weighta
    double precision, parameter   :: alpha   = 1.0d0, pi=3.1415926d0, radian=pi/180.0d0
    integer, parameter            :: nbasis  = 18, ndim=6, natom=4
    double precision, parameter   :: vpesmin = -114.158373618668d0
    integer, parameter            :: ninput  = 17
    integer                       :: noutput, nhid, nlayer, ifunc, nwe, nodemax
    integer                       :: nscale
    integer, allocatable          :: nodes(:)
    double precision, allocatable :: weighta(:,:,:), biasa(:,:)
    double precision, allocatable :: pdela(:), pavga(:)
end module nnparam

!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_ch2o_pipnn_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, path, p, g, d )
        use nnparam, only: alpha, biasa, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, nodemax, &
            & nodes, noutput, nscale, nwe, pavga, pdela, pi, radian, vpesmin, weighta
        implicit none
        ! number of electronic state
        integer, parameter              :: nstates = 1
        integer, parameter              :: natoms  = 4
        integer, intent(in)             :: igrad
        character(len=1024), intent(in) :: path
        double precision, intent(in)    :: x(natoms,3)
        double precision, intent(out)   :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out)   :: d(nstates,nstates,natoms,3)
        double precision                :: v, tx(3,natoms)
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
            call ch2opipnn(tx, v)
            deallocate(nodes)
            deallocate(weighta)
            deallocate(biasa)
            deallocate(pdela)
            deallocate(pavga)
        else
            write (*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
    end subroutine evaluate
end module chempotpy_ch2o_pipnn_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_ch2o_pipnn_pes, only: evaluate
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

!> @brief Evaluate ch2opipNN.
subroutine ch2opipnn( ct, vpes )
    use chempotpy_geometry, only: all_pair_distances
    use nnparam,            only: alpha, biasa, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, nodemax, &
        & nodes, noutput, nscale, nwe, pavga, pdela, pi, radian, vpesmin, weighta
    implicit none
    integer          :: j, lpes
    double precision :: rb(ndim), xbond(ndim), basis(nbasis), tmp1, txinput(nbasis-1)
    double precision :: ct(3,natom), cy(3,natom)
    double precision :: vpes, vpesa, vpesb, vpesc, vpesh, vnn, rfx(6)
    basis = 0.d0
    call all_pair_distances(transpose(ct), rb)
    xbond(:) = dexp(-rb(:)/alpha)
    call bemsav(xbond, basis)
    do j=1, nbasis-1
        txinput(j) = basis(j+1)
    enddo
    call getpota(txinput, vpes)
    !      vnn=min(vpes,5d0)
    vnn  = vpes
    vpes = vnn
end subroutine ch2opipnn

!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weights.txt'
!-->  biases saved in 'biases.txt'
!-->  one has to call this subroutine one and only one before calling the getpot() subroutine
!> @brief Evaluate pes_init.
subroutine pes_init( path )
    use nnparam, only: alpha, biasa, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, nodemax, &
        & nodes, noutput, nscale, nwe, pavga, pdela, pi, radian, vpesmin, weighta
    implicit none
    character(len=1024), intent(in) :: path
    integer                         :: file_unit_4
    integer                         :: i, ihid, iwe, inode1, inode2, ilay1, ilay2
    integer                         :: ibasis, npd, iterm, ib, nfile
    character(len=80)               :: f1, line
    character(len=1024)             :: file_path1, file_path2
    file_path1 = trim(path)//"/CH2O/weights.txt"
    file_path2 = trim(path)//"/CH2O/biases.txt"
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
        &biasa(nodemax,2:nlayer))
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
subroutine getpota( x, vpot )
    use nnparam, only: alpha, biasa, ifunc, natom, nbasis, ndim, nhid, ninput, nlayer, nodemax, &
        & nodes, noutput, nscale, nwe, pavga, pdela, pi, radian, vpesmin, weighta
    implicit none
    integer                    :: i, inode1, inode2, ilay1, ilay2
    integer                    :: j, k, neu1, neu2, ndriv
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

!> @brief Evaluate tranfun.
function tranfun( x, ifunc )
    implicit none
    integer          :: ifunc
    double precision :: tranfun, x
    if( ifunc .eq. 1 )then
        tranfun = dtanh(x)
    else if( ifunc .eq. 2 )then
        tranfun = 1d0/(1d0+dexp(-x))
    else if( ifunc .eq. 3 )then
        tranfun = x
    endif
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> @brief Evaluate emsav.
function emsav( x, c ) result(v)
    implicit none
    double precision, dimension(1:6)::x
    double precision, dimension(0:17)::c
    double precision :: v
    ! ::::::::::::::::::::
    double precision, dimension(0:17)::p
    call bemsav(x, p)
    v = dot_product(p, c)
end function emsav

!> @brief Evaluate bemsav.
subroutine bemsav( x, p )
    implicit none
    double precision, dimension(1:6), intent(in)::x
    double precision, dimension(0:17), intent(out)::p
    ! ::::::::::::::::::::
    double precision, dimension(0:10)::m
    call evmono(x, m)
    call evpoly(m, p)
end subroutine bemsav

!> @brief Evaluate evmono.
subroutine evmono( x, m )
    implicit none
    double precision, dimension(1:6), intent(in)::x
    double precision, dimension(0:10), intent(out)::m
    !::::::::::::::::::::
    m(0)  = 1.0d0
    m(1)  = x(6)
    m(2)  = x(5)
    m(3)  = x(3)
    m(4)  = x(4)
    m(5)  = x(2)
    m(6)  = x(1)
    m(7)  = m(2)*m(3)
    m(8)  = m(3)*m(4)
    m(9)  = m(2)*m(5)
    m(10) = m(4)*m(5)
end subroutine evmono

!> @brief Evaluate evpoly.
subroutine evpoly( m, p )
    implicit none
    double precision, dimension(0:10), intent(in)::m
    double precision, dimension(0:17), intent(out)::p
    !::::::::::::::::::::
    p(0)  = m(0)
    p(1)  = m(1)
    p(2)  = m(2) + m(3)
    p(3)  = m(4) + m(5)
    p(4)  = m(6)
    p(5)  = p(1)*p(2)
    p(6)  = m(7)
    p(7)  = p(1)*p(3)
    p(8)  = m(8) + m(9)
    p(9)  = m(10)
    p(10) = p(2)*p(3) - p(8)
    p(11) = p(1)*p(4)
    p(12) = p(4)*p(2)
    p(13) = p(4)*p(3)
    p(14) = p(1)*p(1)
    p(15) = p(2)*p(2) - p(6) - p(6)
    p(16) = p(3)*p(3) - p(9) - p(9)
    p(17) = p(4)*p(4)
end subroutine evpoly
