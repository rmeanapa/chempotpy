!> @brief Implement CH2O_PIPNN, a potential energy surface for CH2O.
module chempotpy_ch2o_pipnn_pes
    use chempotpy_ch2o_pipnn_parameters, only: activation_function, alpha, hidden_1_biases, &
        & hidden_1_count, hidden_2_biases, hidden_2_count, hidden_hidden_weights, &
        & hidden_output_weights, input_hidden_weights, natom, nbasis, ndim, ninput, nscale, &
        & output_biases, scale_averages, scale_deltas, vpesmin
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the CH2O_PIPNN potential energy surface for CH2O.
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
            call ch2opipnn(tx, v)
        else
            write (*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
    end subroutine evaluate

    !> @brief Evaluate ch2opipNN.
    subroutine ch2opipnn( ct, vpes )
        use chempotpy_geometry, only: all_pair_distances
        use chempotpy_kinds,    only: wp
        implicit none
        integer  :: j, lpes
        real(wp) :: rb(ndim), xbond(ndim), basis(nbasis), tmp1, txinput(nbasis-1)
        real(wp) :: ct(3,natom), cy(3,natom)
        real(wp) :: vpes, vpesa, vpesb, vpesc, vpesh, vnn, rfx(6)
        basis = 0.e0_wp
        call all_pair_distances(transpose(ct), rb)
        xbond(:) = exp(-rb(:)/alpha)
        call bemsav(xbond, basis)
        do j=1, nbasis-1
            txinput(j) = basis(j+1)
        enddo
        call getpota(txinput, vpes)
        !      vnn=min(vpes,5d0)
        vnn  = vpes
        vpes = vnn
    end subroutine ch2opipnn

    !> @brief Evaluate the immutable CH2O_PIPNN neural network.
    subroutine getpota( x, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot
        integer               :: i, inode1, inode2
        real(wp)              :: input_values(ninput), hidden_1_values(hidden_1_count)
        real(wp)              :: hidden_2_values(hidden_2_count), output_value
        do i=1, ninput
            input_values(i) = (x(i)-scale_averages(i, 1))/scale_deltas(i, 1)
        enddo
        do inode1=1, hidden_1_count
            hidden_1_values(inode1) = hidden_1_biases(inode1, 1)
            do inode2=1, ninput
                hidden_1_values(inode1) = hidden_1_values(inode1) + input_values(inode2)* input_hidden_weights(inode2, &
                    & inode1, 1)
            enddo
            hidden_1_values(inode1) = tranfun(hidden_1_values(inode1), activation_function)
        enddo
        do inode1=1, hidden_2_count
            hidden_2_values(inode1) = hidden_2_biases(inode1, 1)
            do inode2=1, hidden_1_count
                hidden_2_values(inode1) = hidden_2_values(inode1) + hidden_1_values(inode2)* &
                    & hidden_hidden_weights(inode2, inode1, 1)
            enddo
            hidden_2_values(inode1) = tranfun(hidden_2_values(inode1), activation_function)
        enddo
        output_value = output_biases(1, 1)
        do inode2=1, hidden_2_count
            output_value = output_value + hidden_2_values(inode2)*hidden_output_weights(inode2, 1, 1)
        enddo
        vpot = output_value*scale_deltas(nscale, 1)+scale_averages(nscale, 1)
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
    !> @brief Evaluate emsav.
    function emsav( x, c ) result(v)
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), dimension(1:6)::x
        real(wp), dimension(0:17)::c
        real(wp) :: v
        real(wp), dimension(0:17)::p
        call bemsav(x, p)
        v = dot_product(p, c)
    end function emsav

    !> @brief Evaluate bemsav.
    subroutine bemsav( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), dimension(1:6), intent(in)::x
        real(wp), dimension(0:17), intent(out)::p
        real(wp), dimension(0:10)::m
        call evmono(x, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), dimension(1:6), intent(in)::x
        real(wp), dimension(0:10), intent(out)::m
        m(0)  = 1.0e0_wp
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
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), dimension(0:10), intent(in)::m
        real(wp), dimension(0:17), intent(out)::p
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
end module chempotpy_ch2o_pipnn_pes

!> @brief Expose CH2O_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_ch2o_pipnn_pes, only: evaluate
    use chempotpy_kinds,          only: wp
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
