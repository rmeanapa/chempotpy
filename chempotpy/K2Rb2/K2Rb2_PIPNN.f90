!> @brief Implement K2Rb2_PIPNN, a potential energy surface for K2Rb2.
module chempotpy_k2rb2_pipnn_pes
    use chempotpy_k2rb2_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_1_count, hidden_2_biases, &
        & hidden_2_count, hidden_hidden_weights, hidden_output_weights, input_hidden_weights, &
        & network_count, ninput, nscale, output_biases, scale_averages, scale_deltas, &
        & vpescut
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the K2Rb2_PIPNN potential energy surface for K2Rb2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units, only: electronvolt_to_inverse_centimeter
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
            call k2rb2_pes_interface(tx, v)
        else
            write (*, *) 'Only energy is available'
        endif
        v = v/electronvolt_to_inverse_centimeter
        do istate=1, nstates
            p(istate) = v
        enddo
    end subroutine evaluate

    ! q in KKRbRb order and angstrom
    !> @brief Evaluate k2rb2_pes_interface.
    subroutine k2rb2_pes_interface( q, vpes )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: lxc
        real(wp) :: q(3,4)
        real(wp) :: vpes, vpesa, vpesb, vpesc
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
            call k2rb2nn(q, vpes, vpesa, vpesb, vpesc)
            if( vpes .lt. -1.e3_wp .or. vpes .gt. vpescut )then
                vpes  = vpescut
                vpesa = vpescut
                vpesb = vpescut
                vpesc = vpescut
            endif
            return
        endif
    end
    !> @brief Evaluate k2rb2NN.
    subroutine k2rb2nn( ct, vpes, vpesa, vpesb, vpesc )
        use chempotpy_geometry, only: all_pair_distances
        use chempotpy_kinds,    only: wp
        implicit none
        integer  :: i
        real(wp) :: rb(6), xbond(6), basis(0:32)
        real(wp) :: txinput(32)
        ! ct in KKRbRb order and angstrom
        real(wp)            :: ct(3,4)
        real(wp)            :: vpesa, vpesb, vpesc, vpes
        real(wp), parameter :: alpha = 10.e0_wp ! bohr
        ! bohr to angstrom
        real(wp), parameter :: btoa = 0.529177249e0_wp
        real(wp), parameter :: pi   = acos(-1.e0_wp)
        basis = 0.e0_wp
        call all_pair_distances(transpose(ct), rb)
        xbond(:) = exp(-rb(:)/(alpha*btoa))
        call bemsav(xbond, basis)
        do i=1, 32
            txinput(i) = basis(i)
        enddo
        call getpota(txinput, vpesa)
        call getpotb(txinput, vpesb)
        call getpotc(txinput, vpesc)
        vpes = (vpesa+vpesb+vpesc)/3.0e0_wp
    end subroutine k2rb2nn

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
        integer  :: i, lxc
        real(wp) :: xct(3,4), xc(3,4), xvec(3,6), rb(6)
        real(wp) :: rmin
        ! lxc=-1,  unphysical point
        ! lxc=0,   normal
        lxc       = 0
        xc        = xct
        xvec(:,1) = xc(:,2)-xc(:,1)  !  K2 ->K1
        xvec(:,2) = xc(:,3)-xc(:,1)  !  Rb3->K1
        xvec(:,3) = xc(:,4)-xc(:,1)  !  Rb4->K1
        xvec(:,4) = xc(:,3)-xc(:,2)  !  Rb3->K2
        xvec(:,5) = xc(:,4)-xc(:,2)  !  Rb4->K2
        xvec(:,6) = xc(:,4)-xc(:,3)  !  Rb4->Rb3
        rb(1)     = sqrt(dot_product(xvec(:,1), xvec(:,1)))
        rb(2)     = sqrt(dot_product(xvec(:,2), xvec(:,2)))
        rb(3)     = sqrt(dot_product(xvec(:,3), xvec(:,3)))
        rb(4)     = sqrt(dot_product(xvec(:,4), xvec(:,4)))
        rb(5)     = sqrt(dot_product(xvec(:,5), xvec(:,5)))
        rb(6)     = sqrt(dot_product(xvec(:,6), xvec(:,6)))
        ! remove some unphysical points
        rmin = minval(rb)
        if( rmin .le. 0.5e0_wp )then
            lxc = -1
            return
        endif
    end subroutine before_vpes

    !------->
    !------->This is the subroutine for PIPs up to the order
    !------->of 3 of A2B2 type molecule.
    !------->
    !        double precision :: x(1:6)
    !        double precision :: p(0:32)
    !        call bemsav(x,p)
    !> @brief Evaluate bemsav.
    subroutine bemsav( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:6)
        real(wp), intent(out) :: p(0:32)
        real(wp)              :: m(0:16)
        call evmono(x, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:6)
        real(wp), intent(out) :: m(0:16)
        m(0)  = 1.0e0_wp
        m(1)  = x(6)
        m(2)  = x(5)
        m(3)  = x(4)
        m(4)  = x(3)
        m(5)  = x(2)
        m(6)  = x(1)
        m(7)  = m(3)*m(4)
        m(8)  = m(2)*m(5)
        m(9)  = m(2)*m(4)
        m(10) = m(3)*m(5)
        m(11) = m(2)*m(3)
        m(12) = m(4)*m(5)
        m(13) = m(2)*m(7)
        m(14) = m(2)*m(10)
        m(15) = m(2)*m(12)
        m(16) = m(3)*m(12)
    end subroutine evmono

    !> @brief Evaluate evpoly.
    subroutine evpoly( m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: m(0:16)
        real(wp), intent(out) :: p(0:32)
        p(0)  = m(0)
        p(1)  = m(1)
        p(2)  = m(2)+m(3)+m(4)+m(5)
        p(3)  = m(6)
        p(4)  = p(1)*p(2)
        p(5)  = m(7)+m(8)
        p(6)  = m(9)+m(10)
        p(7)  = m(11)+m(12)
        p(8)  = p(1)*p(3)
        p(9)  = p(3)*p(2)
        p(10) = p(1)*p(1)
        p(11) = p(2)*p(2)-p(7)-p(6)-p(5)-p(7)-p(6)- p(5)
        p(12) = p(3)*p(3)
        p(13) = p(1)*p(5)
        p(14) = p(1)*p(6)
        p(15) = p(1)*p(7)
        p(16) = m(13)+m(14)+m(15)+m(16)
        p(17) = p(1)*p(9)
        p(18) = p(3)*p(5)
        p(19) = p(3)*p(6)
        p(20) = p(3)*p(7)
        p(21) = p(1)*p(4)
        p(22) = p(1)*p(11)
        p(23) = p(2)*p(5)-p(16)
        p(24) = p(2)*p(6)-p(16)
        p(25) = p(2)*p(7)-p(16)
        p(26) = p(1)*p(8)
        p(27) = p(3)*p(11)
        p(28) = p(1)*p(12)
        p(29) = p(3)*p(9)
        p(30) = p(1)*p(10)
        p(31) = p(2)*p(11)-p(25)-p(24)-p(23)
        p(32) = p(3)*p(12)
    end subroutine evpoly
end module chempotpy_k2rb2_pipnn_pes

!> @brief Expose K2Rb2_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_k2rb2_pipnn_pes, only: evaluate
    use chempotpy_kinds,           only: wp
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
