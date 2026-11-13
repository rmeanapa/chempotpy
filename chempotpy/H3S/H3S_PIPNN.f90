module chempotpy_h3s_pipnn_pes
    use chempotpy_h3s_pipnn_parameters, only: activation_function, alpha, hidden_1_biases, &
        & hidden_1_count, hidden_2_biases, hidden_2_count, hidden_hidden_weights, &
        & hidden_output_weights, input_hidden_weights, natom, nbasis, ndata, ndim, ninput, nscale, &
        & output_biases, para2, scale_averages, scale_deltas, vpesmin
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3S_PIPNN potential energy surface for H3S.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units, only: electronvolt_to_kilocalorie_per_mole
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 4
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: v, va, vb, vc
        real(wp)              :: tx(3,natoms), tg(3,natoms)
        integer               :: iatom, idir, j, istate
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp; tg = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        call evvdvdx(tx, v, tg, igrad)
        do istate=1, nstates
            p(istate) = v/electronvolt_to_kilocalorie_per_mole
        enddo
        do iatom=1, natoms
            do idir=1, 3
                g(1,iatom,idir) = tg(idir,iatom)/electronvolt_to_kilocalorie_per_mole
            enddo
        enddo
        !g=0.d0
        d = 0.e0_wp
    end subroutine evaluate

    !> @brief Evaluate pipNN.
    subroutine pipnn( ct, vpes, vpesa, vpesb, vpesc )
        use chempotpy_geometry, only: all_pair_distances
        use chempotpy_kinds,    only: wp
        implicit none
        integer  :: j, lpes
        real(wp) :: rb(ndim), xbond(ndim), basis(nbasis), tmp1, txinput(nbasis-1)
        real(wp) :: ct(3,natom), cy(3,natom)
        real(wp) :: m(30)
        real(wp) :: vpes, vpesa, vpesb, vpesc, vpesh, vnn
        real(wp) :: dvdg(ninput)
        integer  :: ndriv
        basis = 0.e0_wp
        call all_pair_distances(transpose(ct), rb)
        xbond(:) = exp(-rb(:)/alpha)
        call bemsav(xbond, m, basis)
        do j=1, nbasis-1
            txinput(j) = basis(j+1)
        enddo
        call getpota(txinput, vpesa, dvdg, ndriv)
        vnn  = vpesc
        vpes = vnn
        !-->      write(*,*)vpes
    end subroutine pipnn

    !> @brief Evaluate the first released H3S neural network and its input derivatives.
    subroutine getpota( x, vpot, dvdg, ndriv )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        integer, intent(in)   :: ndriv
        real(wp), intent(out) :: vpot, dvdg(ninput)
        integer               :: i, j, k
        real(wp)              :: input_values(ninput), hidden_1_values(hidden_1_count)
        real(wp)              :: hidden_2_values(hidden_2_count), output_value, derivative_sum
        do i=1, ninput
            input_values(i) = (x(i)-scale_averages(i, 1))/scale_deltas(i, 1)
        enddo
        do j=1, hidden_1_count
            hidden_1_values(j) = hidden_1_biases(j, 1)
            do i=1, ninput
                hidden_1_values(j) = hidden_1_values(j) + input_values(i)*input_hidden_weights(i, j, 1)
            enddo
            hidden_1_values(j) = tranfun(hidden_1_values(j), activation_function)
        enddo
        do k=1, hidden_2_count
            hidden_2_values(k) = hidden_2_biases(k, 1)
            do j=1, hidden_1_count
                hidden_2_values(k) = hidden_2_values(k) + hidden_1_values(j)*hidden_hidden_weights(j, k, 1)
            enddo
            hidden_2_values(k) = tranfun(hidden_2_values(k), activation_function)
        enddo
        output_value = output_biases(1, 1)
        do k=1, hidden_2_count
            output_value = output_value + hidden_2_values(k)*hidden_output_weights(k, 1, 1)
        enddo
        vpot = output_value*scale_deltas(nscale, 1)+scale_averages(nscale, 1)
        dvdg = 0.e0_wp
        if( ndriv .eq. 1 )then
            do i=1, ninput
                do k=1, hidden_2_count
                    derivative_sum = 0.e0_wp
                    do j=1, hidden_1_count
                        derivative_sum = derivative_sum + hidden_hidden_weights(j, k, 1)* input_hidden_weights(i, j, &
                            & 1)*(1.e0_wp-hidden_1_values(j)**2)
                    enddo
                    dvdg(i) = dvdg(i) + hidden_output_weights(k, 1, 1)*derivative_sum* (1.e0_wp-hidden_2_values(k)**2)
                enddo
                dvdg(i) = dvdg(i)*scale_deltas(nscale, 1)/scale_deltas(i, 1)
            enddo
        endif
    end subroutine getpota

    !> @brief Evaluate tranfun.
    function tranfun( x, ifunc )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: ifunc
        real(wp) :: tranfun, x
        if( ifunc .eq. 1 )then
            tranfun = tanh(x)
        elseif( ifunc .eq. 2 )then
            tranfun = 1e0_wp/(1e0_wp+exp(-x))
        elseif( ifunc .eq. 3 )then
            tranfun = x
        endif
    end
    !> @brief Evaluate evvdvdx.
    subroutine evvdvdx( xcart, v, dvdx, ndriv )
        use chempotpy_geometry, only: all_pair_distance_jacobian, all_pair_distances
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), parameter :: alpha = 1e0_wp
        integer             :: j, ndriv
        real(wp)            :: v, dvdx(1:12), p(0:22), r(1:6)
        real(wp)            :: xcart(1:3,1:4), m(0:29), xmorse(1:6)
        real(wp)            :: dvdr(1:6), drdx(1:6,1:12)
        real(wp)            :: txinput(1:22), dvdg(1:22), dpdr(1:6,0:22)
        call all_pair_distances(transpose(xcart), r)
        xmorse = exp(-r/alpha)
        call bemsav(xmorse, m, p)
        do j=1, 22
            txinput(j) = p(j)
        enddo
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
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: xmorse(1:6)
        real(wp), intent(out) :: p(0:22)
        real(wp)              :: m(0:29)
        call evmono(xmorse, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:6)
        real(wp), intent(out) :: m(0:29)
        m(0)  = 1.e0_wp
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
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: m(0:29)
        real(wp), intent(out) :: p(0:22)
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
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), intent(in)  :: dvdg(1:22), dpdr(1:6,0:22), drdx(1:6,1:12)
        real(wp), intent(out) :: dvdx(1:12)
        real(wp)              :: dvdr(1:6)
        dvdr = matmul(dpdr(:,1:22), dvdg)
        call transform_internal_gradient_to_cartesian(dvdr, drdx, dvdx)
    end subroutine evdvdx

    !> @brief Evaluate evdvdr.
    subroutine evdvdr( r, m, p, dpdr )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j
        real(wp) :: r(1:6), dmsdr(1:6,1:6), dmdr(1:6,0:29)
        real(wp) :: dvdr(1:6), dpdr(1:6,0:22)
        real(wp) :: m(0:29), p(0:22)
        dvdr(:) = 0e0_wp
        call evdmsdr(r, dmsdr)
        call evdmdr(m, dmsdr, dmdr)
        call evdpdr(p, dmdr, dpdr)
    end subroutine evdvdr

    !> @brief Evaluate evdrdx.
    subroutine evdrdx( r, x, drdx )
        use chempotpy_geometry, only: all_pair_distance_jacobian
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp) :: r(1:6), drdx(1:6,1:12), x(1:12)
        call all_pair_distance_jacobian(x, r, drdx)
    end subroutine evdrdx

    !> @brief Evaluate evdmsdr.
    subroutine evdmsdr( r, dmsdr )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), parameter :: alpha = 1e0_wp
        integer             :: i, j
        real(wp)            :: dmsdr(1:6,1:6), r(1:6)
        dmsdr(:,:) = 0e0_wp
        !Morse term ms = exp(-r/alpha)
        ! dmsdr(i,j)=0  i!=j
        ! dmsdr(i,i)= -(1/alpha)*Exp(-r(i)/alpha)
        do i=1, 6
            dmsdr(i,i) = -(1e0_wp/alpha)*exp(-r(i)/alpha)
        enddo
    end subroutine evdmsdr

    !> @brief Evaluate evdmdr.
    subroutine evdmdr( m, dmsdr, dmdr )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j
        real(wp) :: dmsdr(1:6,1:6), dmdr(1:6,0:29), m(0:29)
        do i=1, 6
            dmdr(i,0)  = 0.e0_wp
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
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j
        real(wp) :: dmdr(1:6,0:29), dpdr(1:6,0:22)
        real(wp) :: p(0:22)
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
end module chempotpy_h3s_pipnn_pes

!> @brief Expose H3S_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_h3s_pipnn_pes, only: evaluate
    use chempotpy_kinds,         only: wp
    implicit none
    integer, parameter              :: natoms  = 4
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
