!     The ^2A prime PES for C2N system at the MC-PDFT/AVTZ level.
!> @brief Implement C2N_PIPNN_Ap, a potential energy surface for C2N.
module chempotpy_c2n_pipnn_ap_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the C2N_PIPNN_Ap potential energy surface for C2N.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_kinds,    only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: v
        real(wp)              :: rb(3)
        integer               :: iatom, idir, j, istate
        !initialize
        v = 0.e0_wp; p = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        ! order of r_CC,r_NC,r_NC
        rb    = 0.e0_wp
        rb(1) = cartesian_distance(x(1,:), x(2,:))
        rb(2) = cartesian_distance(x(1,:), x(3,:))
        rb(3) = cartesian_distance(x(2,:), x(3,:))
        if( igrad == 0 )then
            call c2n_2a_prime(rb, v)
        else
            write (*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
        g = 0.e0_wp
        d = 0.e0_wp
    end subroutine evaluate

    !     PIP-NN PES of total energy for the C2N(^2A') system.
    !     where,
    !       rNC,rCC: angstrom
    !       theta: degree
    !       Vpes: total energy, eV
    !> @brief Evaluate c2n_2a_prime.
    subroutine c2n_2a_prime( rb, vpes )
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter  :: npip    = 12
        real(wp), parameter :: alpha   = 1.5e0_wp ! angstrom
        real(wp), parameter :: vpescut = 5.e5_wp ! eV
        real(wp), parameter :: pi      = acos(-1.e0_wp)
        integer             :: i
        real(wp)            :: vpesa(1), vpesb(1), vpesc(1), vpes
        real(wp)            :: rb(3), xbond(3), basis(0:npip), txinput(npip)
        if( minval(rb) .lt. 0.79999999e0_wp )then
            vpes  = vpescut
            vpesa = vpescut
            vpesb = vpescut
            vpesc = vpescut
        else
            xbond(:) = exp(-rb(:)/alpha)
            basis    = 0.e0_wp
            call bemsav(xbond, basis)
            do i=1, npip
                txinput(i) = basis(i)
            enddo
            call nnfit(npip, 1, txinput, 1, vpesa)
            call nnfit(npip, 1, txinput, 2, vpesb)
            call nnfit(npip, 1, txinput, 3, vpesc)
            vpes = (vpesa(1)+vpesb(1)+vpesc(1))/3.0e0_wp
            if( vpes .gt. vpescut )then
                vpes  = vpescut
                vpesa = vpescut
                vpesb = vpescut
                vpesc = vpescut
            endif
        endif
    end subroutine c2n_2a_prime

    !> @brief Evaluate one member of the immutable three-network ensemble.
    subroutine nnfit( ndim, ntot, r0, network, vx )
        use chempotpy_c2n_pipnn_ap_parameters, only: hidden_1_biases, hidden_1_count, hidden_2_biases, &
            & hidden_2_count, hidden_hidden_weights, hidden_output_weights, input_hidden_weights, &
            & input_ranges, ninput, output_biases, output_ranges
        use chempotpy_kinds, only: wp
        implicit none
        integer, intent(in)   :: ndim, ntot, network
        real(wp), intent(in)  :: r0(ndim,ntot)
        real(wp), intent(out) :: vx(ntot)
        real(wp)              :: dv(ndim,ntot)
        integer               :: i
        if( ntot .ge. 24 )then
            call nsimx(ntot, r0, vx, 0, dv, ninput, hidden_1_count, hidden_2_count, &
                & input_ranges(:, :, network), input_hidden_weights(:, :, network), &
                & hidden_1_biases(:, network), hidden_hidden_weights(:, :, network), &
                & hidden_2_biases(:, network), hidden_output_weights(:, network), &
                & output_biases(network), output_ranges(:, network))
        else
            do i=1, ntot
                call nsim(r0(1,i), vx(i), 0, dv(:,i), ninput, hidden_1_count, hidden_2_count, &
                    & input_ranges(:, :, network), input_hidden_weights(:, :, network), &
                    & hidden_1_biases(:, network), hidden_hidden_weights(:, :, network), &
                    & hidden_2_biases(:, network), hidden_output_weights(:, network), &
                    & output_biases(network), output_ranges(:, network))
            enddo
        endif
    end subroutine nnfit

    ! Simulate a neural network with two hidden layers: n0-n1-n2-1
    ! Blas routines used in this subroutine: dgemv, ddot
    !> @brief Evaluate nsim.
    subroutine nsim( r0, v, idv, dv, n0, n1, n2, rg, w1, b1, w2, b2, w3, b3, vg )
        use chempotpy_kinds, only: wp
        implicit none
        integer, intent(in)   :: n0, n1, n2, idv
        real(wp), intent(in)  :: r0(n0), rg(2,n0), vg(2)
        real(wp), intent(in)  :: w1(n0,n1), b1(n1)
        real(wp), intent(in)  :: w2(n1,n2), b2(n2)
        real(wp), intent(in)  :: w3(n2), b3
        real(wp), intent(out) :: v, dv(n0)
        integer               :: i, j, k
        real(wp)              :: r(n0), rgg(n0), vgg, ax(n1), bx(n2)
        real(wp)              :: dvtm, rtmp, rt1(n1), rt2(n2)
        real(wp), external    :: ddot
        v   = 0.e0_wp
        r   = r0
        vgg = vg(2)-vg(1)
        ! mapminmax [-1,1]
        do i=1, n0
            rgg(i) = rg(2,i)-rg(1,i)
            r(i)   = 2.e0_wp*(r(i)-rg(1,i))/rgg(i)-1.e0_wp
        enddo
        ! 1st layer
        rt1 = b1
        call dgemv('t', n0, n1, 1.e0_wp, w1, n0, r, 1, 1.e0_wp, rt1, 1)
        ax = tanh(rt1)
        ! 2nd layer
        rt2 = b2
        call dgemv('t', n1, n2, 1.e0_wp, w2, n1, ax, 1, 1.e0_wp, rt2, 1)
        bx = tanh(rt2)
        ! output layer
        v = b3+ddot(n2, w3, 1, bx, 1)
        !reverse map
        v = vgg*(v+1.e0_wp)/2+vg(1)
        if( idv .ne. 1 ) return
        ! calculate first derivatives, dv(i)=dv/dr(i)
        dv = 0.e0_wp
        do i=1, n0
            do k=1, n2
                dvtm = 0.e0_wp
                do j=1, n1
                    dvtm = dvtm+w2(j,k)*w1(i,j)*(1-ax(j)**2)
                enddo
                dv(i) = dv(i)+w3(k)*dvtm*(1-bx(k)**2)
            enddo
            dv(i) = dv(i)*vgg/rgg(i)
        enddo
    end subroutine nsim

    ! Simulate a neural network with two hidden layers:
    !    n0-n1(tansig)-n2(tansig)-1(pureline)
    ! Blas routines used in this subroutine: dgemm, dgemv
    ! The syntax keeps the same with nsim, except input nt at first.
    ! Chen Jun, 2013/11/26
    ! Description of parameters:
    !    nt  : total number of points in r,v and dv
    !    r   : input of nn, r(n0,nt), r(1:n0,i) is the i-th input
    !    v   : output of nn, v(nt), v(i) is the i-th input
    !    idv : if idv==0, calc. v only, else calc. v and derivatives dv
    !    dv  : derivatives of dv/dr, dv(n0,nt), dv(j,i)=dv(i)/dr(j,i)
    !    n0  : neurons of the input layer, dimension of r
    !    n1  : neurons of the first hidden layer
    !    n2  : neurons of the second hidden layer
    !    rg  : input ranges of the training set, rg(2,n0)
    !    w1  : weights connect input layer with 1-st hidden layer
    !    b1  : biases of 1-st hidden layer
    !    w2  : weights connect 1-st hidden layer with 2-nd hidden layer
    !    b2  : biases of 2-nd hidden layer
    !    w3  : weights connect 2-nd hidden layer with output layer
    !    b3  : biase of output layer
    !    vg  : output range of the training set, vg(2)
    !> @brief Evaluate nsimx.
    subroutine nsimx( nt, r, v, idv, dv, n0, n1, n2, rg, w1, b1, w2, b2, w3, b3, vg )
        use chempotpy_kinds, only: wp
        implicit none
        integer, intent(in)   :: nt, idv, n0, n1, n2
        real(wp), intent(in)  :: r(n0,nt)
        real(wp), intent(out) :: v(nt), dv(n0,nt)
        real(wp), intent(in)  :: w1(n0,n1), b1(n1), w2(n1,n2), b2(n2), w3(n2), b3
        real(wp), intent(in)  :: rg(2,n0), vg(2)
        real(wp)              :: x0(n0,nt), x1(nt,n1), x2(nt,n2), tmp, rgg(n0), vgg
        integer               :: i, j, k, n
        v = 0.e0_wp
        do i=1, n0
            rgg(i) = rg(2,i)-rg(1,i)
        enddo
        vgg = vg(2)-vg(1)
        x0  = r
        do i=1, n0
            x0(i,1:nt) = (x0(i,1:nt)-rg(1,i))/rgg(i)*2.e0_wp-1.e0_wp
        enddo
        do i=1, n1
            x1(1:nt,i) = b1(i)
        enddo
        ! ---- x1=trans(x0(n0,nt))*w1(n0,n1)+x1(nt,n1)
        call dgemm('t', 'n', nt, n1, n0, 1.e0_wp, x0, n0, w1, n0, 1.e0_wp, x1, nt)
        x1 = tanh(x1)
        do i=1, n2
            x2(1:nt,i) = b2(i)
        enddo
        ! ---- x2=x1(nt,n1)*w2(n1,n2)+x2(nt,n2)
        call dgemm('n', 'n', nt, n2, n1, 1.e0_wp, x1, nt, w2, n1, 1.e0_wp, x2, nt)
        x2      = tanh(x2)
        v(1:nt) = b3
        ! ---- v=x2(nt,n2)*w3(n2)+v(nt)
        call dgemv('n', nt, n2, 1.e0_wp, x2, nt, w3, 1, 1.e0_wp, v, 1)
        v = (v+1.e0_wp)/2.e0_wp*vgg+vg(1)
        if( idv .ne. 1 ) return
        ! ---- calculate first derivatives
        dv = 0.e0_wp
        do n=1, nt
            do i=1, n0
                do k=1, n2
                    tmp = 0.e0_wp
                    do j=1, n1
                        tmp = tmp+w2(j,k)*w1(i,j)*(1-x1(n,j)**2)
                    enddo
                    dv(i,n) = dv(i,n)+w3(k)*tmp*(1-x2(n,k)**2)
                enddo
                dv(i,n) = dv(i,n)*vgg/(rgg(i))
            enddo
        enddo
    end subroutine nsimx

    !------->This is the subroutine for PIPs up to the order of 3 of A2B1
    !------->type molecule.
    !> @brief Evaluate bemsav.
    subroutine bemsav( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:3)
        real(wp), intent(out) :: p(0:12)
        real(wp)              :: m(0:4)
        call evmono(x, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:3)
        real(wp), intent(out) :: m(0:4)
        m(0) = 1.0e0_wp
        m(1) = x(3)
        m(2) = x(2)
        m(3) = x(1)
        m(4) = m(1)*m(2)
    end subroutine evmono

    !> @brief Evaluate evpoly.
    subroutine evpoly( m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: m(0:4)
        real(wp), intent(out) :: p(0:12)
        p(0)  = m(0)
        p(1)  = m(1)+m(2)
        p(2)  = m(3)
        p(3)  = m(4)
        p(4)  = p(2)*p(1)
        p(5)  = p(1)*p(1)-p(3)-p(3)
        p(6)  = p(2)*p(2)
        p(7)  = p(2)*p(3)
        p(8)  = p(3)*p(1)
        p(9)  = p(2)*p(5)
        p(10) = p(2)*p(4)
        p(11) = p(1)*p(5)-p(8)
        p(12) = p(2)*p(6)
    end subroutine evpoly
end module chempotpy_c2n_pipnn_ap_pes

!> @brief Expose C2N_PIPNN_Ap through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_c2n_pipnn_ap_pes, only: evaluate
    use chempotpy_kinds,            only: wp
    implicit none
    integer, parameter              :: natoms  = 3
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
