!> @brief Implement N2HOCp_PIPNN, a potential energy surface for N2HOCp.
module chempotpy_n2hocp_pipnn_pes
    use chempotpy_n2hocp_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_1_count, hidden_2_biases, &
        & hidden_2_count, hidden_hidden_weights, hidden_output_weights, input_hidden_weights, &
        & network_count, ninput, nscale, output_biases, scale_averages, scale_deltas, &
        & alpha, natom, nbasis, ndim
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the N2HOCp_PIPNN potential energy surface for N2HOCp.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 5
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: v, va, vb, vc
        real(wp)              :: coord(3,natoms), r(10)
        integer               :: iatom, jatom, idir, j, istate
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                coord(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        j = 0
        do iatom=1, natoms
            do jatom=iatom+1, natoms
                j = j+1
                r(j)=sqrt((coord(1,iatom)-coord(1,jatom))**2 + + (coord(2,iatom)-coord(2,jatom))**2 + + (coord(3, &
                    & iatom)- coord(3,jatom))**2)
            enddo
        enddo
        if( igrad == 0 )then
            call n2chopipnn(r, v, va, vb, vc)
        else
            write (*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
        g = 0.e0_wp
        d = 0.e0_wp
    end subroutine evaluate

    !> @brief Evaluate distance.
    real(wp) function distance( a, b )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i
        real(wp) :: a(3), b(3)
        distance = 0.e0_wp
        do i=1, 3
            distance = distance+(a(i)-b(i))**2
        enddo
        distance = sqrt(distance)
    end
    !> @brief Evaluate N2CHOpipNN.
    subroutine n2chopipnn( rb, vpes, vpesa, vpesb, vpesc )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j, k, lpes
        real(wp) :: rb(ndim), xbond(ndim), basis(0:nbasis), tmp1, txinput(nbasis)
        real(wp) :: vpes, vpesa, vpesb, vpesc, vpesh, vnn
        basis    = 0.e0_wp
        xbond(:) = exp(-rb(:)/alpha)
        call bemsav21113(xbond, basis)
        do j=1, nbasis
            txinput(j) = basis(j)
        enddo
        call getpota(txinput, vpesa)
        call getpotb(txinput, vpesb)
        call getpotc(txinput, vpesc)
        vpes = (vpesa+vpesb+vpesc)/3e0_wp
    end subroutine n2chopipnn

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
    subroutine getpot( x, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot
        call evaluate_network(x, 1, vpot)
    end subroutine getpot

    !> @brief Evaluate emsav21113.
    function emsav21113( x, c ) result(v)
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(1:10)
        real(wp) :: c(0:167)
        real(wp) :: v
        real(wp) :: p(0:167)
        call bemsav21113(x, p)
        v = dot_product(p, c)
    end function emsav21113

    !> @brief Evaluate bemsav21113.
    subroutine bemsav21113( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:10)
        real(wp), intent(out) :: p(0:167)
        real(wp)              :: m(0:21)
        call evmono21113(x, m)
        call evpoly21113(m, p)
    end subroutine bemsav21113

    !> @brief Evaluate evmono21113.
    subroutine evmono21113( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:10)
        real(wp), intent(out) :: m(0:21)
        m(0)  = 1.e0_wp
        m(1)  = x(10)
        m(2)  = x(9)
        m(3)  = x(8)
        m(4)  = x(7)
        m(5)  = x(4)
        m(6)  = x(6)
        m(7)  = x(3)
        m(8)  = x(5)
        m(9)  = x(2)
        m(10) = x(1)
        m(11) = m(4)*m(5)
        m(12) = m(5)*m(6)
        m(13) = m(4)*m(7)
        m(14) = m(6)*m(7)
        m(15) = m(5)*m(8)
        m(16) = m(4)*m(9)
        m(17) = m(7)*m(8)
        m(18) = m(6)*m(9)
        m(19) = m(8)*m(9)
        m(20) = m(5)*m(17)
        m(21) = m(4)*m(18)
    end subroutine evmono21113

    !> @brief Evaluate evpoly21113.
    subroutine evpoly21113( m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: m(0:21)
        real(wp), intent(out) :: p(0:167)
        p(0)   = m(0)
        p(1)   = m(1)
        p(2)   = m(2)
        p(3)   = m(3)
        p(4)   = m(4)+m(5)
        p(5)   = m(6)+m(7)
        p(6)   = m(8)+m(9)
        p(7)   = m(10)
        p(8)   = p(1)*p(2)
        p(9)   = p(1)*p(3)
        p(10)  = p(2)*p(3)
        p(11)  = p(1)*p(4)
        p(12)  = p(2)*p(4)
        p(13)  = p(3)*p(4)
        p(14)  = m(11)
        p(15)  = p(1)*p(5)
        p(16)  = p(2)*p(5)
        p(17)  = p(3)*p(5)
        p(18)  = m(12)+m(13)
        p(19)  = m(14)
        p(20)  = p(4)*p(5)-p(18)
        p(21)  = p(1)*p(6)
        p(22)  = p(2)*p(6)
        p(23)  = p(3)*p(6)
        p(24)  = m(15)+m(16)
        p(25)  = m(17)+m(18)
        p(26)  = m(19)
        p(27)  = p(4)*p(6)-p(24)
        p(28)  = p(5)*p(6)-p(25)
        p(29)  = p(1)*p(7)
        p(30)  = p(2)*p(7)
        p(31)  = p(3)*p(7)
        p(32)  = p(7)*p(4)
        p(33)  = p(7)*p(5)
        p(34)  = p(7)*p(6)
        p(35)  = p(1)*p(1)
        p(36)  = p(2)*p(2)
        p(37)  = p(3)*p(3)
        p(38)  = p(4)*p(4)-p(14)-p(14)
        p(39)  = p(5)*p(5)-p(19)-p(19)
        p(40)  = p(6)*p(6)-p(26)-p(26)
        p(41)  = p(7)*p(7)
        p(42)  = p(1)*p(10)
        p(43)  = p(1)*p(12)
        p(44)  = p(1)*p(13)
        p(45)  = p(2)*p(13)
        p(46)  = p(1)*p(14)
        p(47)  = p(2)*p(14)
        p(48)  = p(3)*p(14)
        p(49)  = p(1)*p(16)
        p(50)  = p(1)*p(17)
        p(51)  = p(2)*p(17)
        p(52)  = p(1)*p(18)
        p(53)  = p(2)*p(18)
        p(54)  = p(3)*p(18)
        p(55)  = p(1)*p(19)
        p(56)  = p(2)*p(19)
        p(57)  = p(3)*p(19)
        p(58)  = p(1)*p(20)
        p(59)  = p(2)*p(20)
        p(60)  = p(3)*p(20)
        p(61)  = p(14)*p(5)
        p(62)  = p(19)*p(4)
        p(63)  = p(1)*p(22)
        p(64)  = p(1)*p(23)
        p(65)  = p(2)*p(23)
        p(66)  = p(1)*p(24)
        p(67)  = p(2)*p(24)
        p(68)  = p(3)*p(24)
        p(69)  = p(1)*p(25)
        p(70)  = p(2)*p(25)
        p(71)  = p(3)*p(25)
        p(72)  = m(20)+m(21)
        p(73)  = p(1)*p(26)
        p(74)  = p(2)*p(26)
        p(75)  = p(3)*p(26)
        p(76)  = p(1)*p(27)
        p(77)  = p(2)*p(27)
        p(78)  = p(3)*p(27)
        p(79)  = p(14)*p(6)
        p(80)  = p(4)*p(25)-p(72)
        p(81)  = p(26)*p(4)
        p(82)  = p(1)*p(28)
        p(83)  = p(2)*p(28)
        p(84)  = p(3)*p(28)
        p(85)  = p(5)*p(24)-p(72)
        p(86)  = p(19)*p(6)
        p(87)  = p(26)*p(5)
        p(88)  = p(4)*p(28)-p(85)
        p(89)  = p(1)*p(30)
        p(90)  = p(1)*p(31)
        p(91)  = p(2)*p(31)
        p(92)  = p(1)*p(32)
        p(93)  = p(2)*p(32)
        p(94)  = p(3)*p(32)
        p(95)  = p(7)*p(14)
        p(96)  = p(1)*p(33)
        p(97)  = p(2)*p(33)
        p(98)  = p(3)*p(33)
        p(99)  = p(7)*p(18)
        p(100) = p(7)*p(19)
        p(101) = p(7)*p(20)
        p(102) = p(1)*p(34)
        p(103) = p(2)*p(34)
        p(104) = p(3)*p(34)
        p(105) = p(7)*p(24)
        p(106) = p(7)*p(25)
        p(107) = p(7)*p(26)
        p(108) = p(7)*p(27)
        p(109) = p(7)*p(28)
        p(110) = p(1)*p(8)
        p(111) = p(1)*p(36)
        p(112) = p(1)*p(9)
        p(113) = p(2)*p(10)
        p(114) = p(1)*p(37)
        p(115) = p(2)*p(37)
        p(116) = p(1)*p(11)
        p(117) = p(2)*p(12)
        p(118) = p(3)*p(13)
        p(119) = p(1)*p(38)
        p(120) = p(2)*p(38)
        p(121) = p(3)*p(38)
        p(122) = p(14)*p(4)
        p(123) = p(1)*p(15)
        p(124) = p(2)*p(16)
        p(125) = p(3)*p(17)
        p(126) = p(4)*p(18)-p(61)
        p(127) = p(4)*p(20)-p(61)
        p(128) = p(1)*p(39)
        p(129) = p(2)*p(39)
        p(130) = p(3)*p(39)
        p(131) = p(5)*p(18)-p(62)
        p(132) = p(19)*p(5)
        p(133) = p(4)*p(39)-p(131)
        p(134) = p(1)*p(21)
        p(135) = p(2)*p(22)
        p(136) = p(3)*p(23)
        p(137) = p(4)*p(24)-p(79)
        p(138) = p(5)*p(25)-p(86)
        p(139) = p(4)*p(27)-p(79)
        p(140) = p(5)*p(28)-p(86)
        p(141) = p(1)*p(40)
        p(142) = p(2)*p(40)
        p(143) = p(3)*p(40)
        p(144) = p(6)*p(24)-p(81)
        p(145) = p(6)*p(25)-p(87)
        p(146) = p(26)*p(6)
        p(147) = p(4)*p(40)-p(144)
        p(148) = p(5)*p(40)-p(145)
        p(149) = p(1)*p(29)
        p(150) = p(2)*p(30)
        p(151) = p(3)*p(31)
        p(152) = p(7)*p(38)
        p(153) = p(7)*p(39)
        p(154) = p(7)*p(40)
        p(155) = p(1)*p(41)
        p(156) = p(2)*p(41)
        p(157) = p(3)*p(41)
        p(158) = p(7)*p(32)
        p(159) = p(7)*p(33)
        p(160) = p(7)*p(34)
        p(161) = p(1)*p(35)
        p(162) = p(2)*p(36)
        p(163) = p(3)*p(37)
        p(164) = p(4)*p(38)-p(122)
        p(165) = p(5)*p(39)-p(132)
        p(166) = p(6)*p(40)-p(146)
        p(167) = p(7)*p(41)
    end subroutine evpoly21113

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
end module chempotpy_n2hocp_pipnn_pes

!> @brief Expose N2HOCp_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_n2hocp_pipnn_pes, only: evaluate
    use chempotpy_kinds,            only: wp
    implicit none
    integer, parameter              :: natoms  = 5
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
