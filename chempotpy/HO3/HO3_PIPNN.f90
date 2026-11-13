!> @brief Implement HO3_PIPNN, a potential energy surface for HO3.
module chempotpy_ho3_pipnn_pes
    use chempotpy_ho3_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_1_count, hidden_2_biases, &
        & hidden_2_count, hidden_hidden_weights, hidden_output_weights, input_hidden_weights, &
        & network_count, ninput, nscale, output_biases, scale_averages, scale_deltas, &
        & vpescut
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the HO3_PIPNN potential energy surface for HO3.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 4
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: v, tx(3,natoms)
        integer               :: iatom, idir, j, istate
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        if( igrad == 0 )then
            call ho3_pes_interface(tx, v)
        else
            write (*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
    end subroutine evaluate

    !> @brief Evaluate ho3_pes_interface.
    subroutine ho3_pes_interface( q, vpes )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: lxc
        real(wp) :: q(3,4) ! q in HOOO order and angstrom.
        real(wp) :: vpes, vpesa, vpesb, vpesc ! in eV
        vpes  = 0.e0_wp
        vpesa = 0.e0_wp
        vpesb = 0.e0_wp
        vpesc = 0.e0_wp
        lxc   = 0
        call before_vpes(q, lxc)
        if( lxc .eq. -1 )then
            vpesa = vpescut
            vpesb = vpescut
            vpesc = vpescut
            vpes  = vpescut
            return
        endif
        if( lxc .eq. 0 )then
            call ho3nn(q, vpes, vpesa, vpesb, vpesc)
            if( vpes .lt. -1.5 .or. vpes .gt. vpescut )then
                vpes  = vpescut
                vpesa = vpescut
                vpesb = vpescut
                vpesc = vpescut
            endif
            return
        endif
    end
    !> @brief Evaluate ho3NN.
    subroutine ho3nn( ct, vpes, vpesa, vpesb, vpesc )
        use chempotpy_geometry, only: all_pair_distances
        use chempotpy_kinds,    only: wp
        implicit none
        integer             :: i
        real(wp)            :: rb(6), xbond(6), basis(0:22)
        real(wp)            :: txinput(22)
        real(wp)            :: ct(3,4), pair_distances(6)
        real(wp)            :: vpesa, vpesb, vpesc, vpes
        real(wp), parameter :: alpha = 2.5e0_wp ! bohr
        real(wp), parameter :: btoa  = 0.529177249e0_wp
        real(wp), parameter :: pi    = acos(-1.e0_wp)
        basis = 0.e0_wp
        call all_pair_distances(transpose(ct), pair_distances)
        rb       = pair_distances([6,5,3,4,2,1])
        xbond(:) = exp(-rb(:)/(alpha*btoa))
        call bemsav(xbond, basis)
        do i=1, 22
            txinput(i) = basis(i)
        enddo
        call getpota(txinput, vpesa)
        call getpotb(txinput, vpesb)
        call getpotc(txinput, vpesc)
        vpes = (vpesa+vpesb+vpesc)/3.0e0_wp
    end subroutine ho3nn

    !> @brief Evaluate one immutable neural network from the released three-model ensemble.
    subroutine evaluate_network( x, network_index, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        integer, intent(in)   :: network_index
        real(wp), intent(out) :: vpot
        integer               :: i, inode1, inode2
        real(wp)              :: input_values(ninput), hidden_1_values(hidden_1_count)
        real(wp)              :: hidden_2_values(hidden_2_count), output_value
        do i=1, ninput
            input_values(i) = (x(i)-scale_averages(i, network_index))/scale_deltas(i, network_index)
        enddo
        do inode1=1, hidden_1_count
            hidden_1_values(inode1) = hidden_1_biases(inode1, network_index)
            do inode2=1, ninput
                hidden_1_values(inode1) = hidden_1_values(inode1) + input_values(inode2)* input_hidden_weights(inode2, &
                    & inode1, network_index)
            enddo
            hidden_1_values(inode1) = tranfun(hidden_1_values(inode1), activation_function)
        enddo
        do inode1=1, hidden_2_count
            hidden_2_values(inode1) = hidden_2_biases(inode1, network_index)
            do inode2=1, hidden_1_count
                hidden_2_values(inode1) = hidden_2_values(inode1) + hidden_1_values(inode2)* &
                    & hidden_hidden_weights(inode2, inode1, network_index)
            enddo
            hidden_2_values(inode1) = tranfun(hidden_2_values(inode1), activation_function)
        enddo
        output_value = output_biases(1, network_index)
        do inode2=1, hidden_2_count
            output_value = output_value + hidden_2_values(inode2)* hidden_output_weights(inode2, 1, network_index)
        enddo
        vpot = output_value*scale_deltas(nscale, network_index)+scale_averages(nscale, network_index)
    end subroutine evaluate_network

    !> @brief Evaluate released neural network 1 from the ensemble.
    subroutine getpota( x, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot
        call evaluate_network(x, 1, vpot)
    end subroutine getpota

    !> @brief Evaluate released neural network 2 from the ensemble.
    subroutine getpotb( x, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot
        call evaluate_network(x, 2, vpot)
    end subroutine getpotb

    !> @brief Evaluate released neural network 3 from the ensemble.
    subroutine getpotc( x, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot
        call evaluate_network(x, 3, vpot)
    end subroutine getpotc

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
    !> @brief Evaluate before_vpes.
    subroutine before_vpes( xct, lxc )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: lxc
        real(wp) :: xct(3,4), xc(3,4), xvec(3,6), rb(6)
        real(wp) :: ro3o4, ro2o4, rh1o4, ro2o3, rh1o3, rh1o2
        ! lxc=-1,  unphysical point
        ! lxc=0,   normal
        lxc       = 0
        xc        = xct
        xvec(:,1) = xc(:,3)-xc(:,4) !  O3->O4
        xvec(:,2) = xc(:,2)-xc(:,4) !  O2->O4
        xvec(:,3) = xc(:,1)-xc(:,4) !  H1->O4
        xvec(:,4) = xc(:,2)-xc(:,3) !  O2->O3
        xvec(:,5) = xc(:,1)-xc(:,3) !  H1->O3
        xvec(:,6) = xc(:,1)-xc(:,2) !  H1->O2
        ro3o4     = sqrt(dot_product(xvec(:,1), xvec(:,1)))
        ro2o4     = sqrt(dot_product(xvec(:,2), xvec(:,2)))
        rh1o4     = sqrt(dot_product(xvec(:,3), xvec(:,3)))
        ro2o3     = sqrt(dot_product(xvec(:,4), xvec(:,4)))
        rh1o3     = sqrt(dot_product(xvec(:,5), xvec(:,5)))
        rh1o2     = sqrt(dot_product(xvec(:,6), xvec(:,6)))
        rb(:)     = (/ro3o4, ro2o4, rh1o4, ro2o3, rh1o3, rh1o2/)
        ! remove some unphysical points
        if(min(rb(3), rb(5), rb(6)).le.0.65e0_wp.or.&
            &     min(rb(1), rb(2), rb(4)).le.0.65e0_wp) then
        lxc = -1
        return
        endif
    end subroutine before_vpes

    !------>
    !This program used to generate PIPs for molecule of A3B1 type
    !up to 3rd order.
    !------>
    !> @brief Evaluate bemsav.
    subroutine bemsav( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:6)
        real(wp), intent(out) :: p(0:22)
        real(wp)              :: m(0:29)
        call evmono(x, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:6)
        real(wp), intent(out) :: m(0:29)
        m(0)  = 1.0e0_wp
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
        p(1)  = m(1) + m(2) + m(3)
        p(2)  = m(4) + m(5) + m(6)
        p(3)  = m(7) + m(8) + m(9)
        p(4)  = m(10) + m(11) + m(12)
        p(5)  = p(1)*p(2) - p(4)
        p(6)  = m(13) + m(14) + m(15)
        p(7)  = p(1)*p(1) - p(3) - p(3)
        p(8)  = p(2)*p(2) - p(6) - p(6)
        p(9)  = m(16)
        p(10) = m(17) + m(18) + m(19) + m(20) + m(21) + m(22)
        p(11) = p(2)*p(3) - p(10)
        p(12) = m(23) + m(24) + m(25) + m(26) + m(27) + m(28)
        p(13) = m(29)
        p(14) = p(1)*p(6) - p(12)
        p(15) = p(1)*p(3) - p(9) - p(9) - p(9)
        p(16) = p(1)*p(4) - p(10)
        p(17) = p(2)*p(7) - p(16)
        p(18) = p(2)*p(4) - p(12)
        p(19) = p(1)*p(8) - p(18)
        p(20) = p(2)*p(6) - p(13) - p(13) - p(13)
        p(21) = p(1)*p(7) - p(15)
        p(22) = p(2)*p(8) - p(20)
    end subroutine evpoly
end module chempotpy_ho3_pipnn_pes

!> @brief Expose HO3_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_ho3_pipnn_pes, only: evaluate
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
