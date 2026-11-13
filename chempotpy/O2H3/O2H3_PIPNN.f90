!> @brief Implement O2H3_PIPNN, a potential energy surface for O2H3.
module chempotpy_o2h3_pipnn_pes
    use chempotpy_o2h3_pipnn_parameters, only: activation_function, alpha, hidden_1_biases, &
        & hidden_1_count, hidden_2_biases, hidden_2_count, hidden_hidden_weights, &
        & hidden_output_weights, input_hidden_weights, natom, nbasis, ndim, ninput, nscale, &
        & output_biases, scale_averages, scale_deltas, vpesmin
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the O2H3_PIPNN potential energy surface for O2H3.
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
        real(wp)              :: v, tx(3,natoms)
        integer               :: iatom, idir, j, istate
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        if( igrad == 0 )then
            do iatom=1, natoms
                do idir=1, 3
                    tx(idir,iatom) = x(iatom,idir)
                enddo
            enddo
            call h3o2pipnn(tx, v)
            do istate=1, nstates
                p(istate) = v
            enddo
        else
            write(*, *) 'Only energy and gradient are available'
        endif
    end subroutine evaluate

    !> @brief Evaluate h3o2pipNN.
    subroutine h3o2pipnn( ct, vpes )
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
    end subroutine h3o2pipnn

    !> @brief Evaluate the immutable O2H3_PIPNN neural network.
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
        real(wp), dimension(1:10)::x
        real(wp), dimension(0:138)::c
        real(wp) :: v
        real(wp), dimension(0:138)::p
        call bemsav(x, p)
        v = dot_product(p, c)
    end function emsav

    !> @brief Evaluate bemsav.
    subroutine bemsav( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), dimension(1:10), intent(in)::x
        real(wp), dimension(0:138), intent(out)::p
        real(wp), dimension(0:169)::m
        call evmono(x, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), dimension(1:10), intent(in)::x
        real(wp), dimension(0:169), intent(out)::m
        m(0)   = 1.0e0_wp
        m(1)   = x(10)
        m(2)   = x(9)
        m(3)   = x(8)
        m(4)   = x(7)
        m(5)   = x(6)
        m(6)   = x(4)
        m(7)   = x(3)
        m(8)   = x(5)
        m(9)   = x(2)
        m(10)  = x(1)
        m(11)  = m(3)*m(4)
        m(12)  = m(2)*m(5)
        m(13)  = m(3)*m(6)
        m(14)  = m(5)*m(6)
        m(15)  = m(2)*m(7)
        m(16)  = m(4)*m(7)
        m(17)  = m(2)*m(4)
        m(18)  = m(3)*m(5)
        m(19)  = m(2)*m(6)
        m(20)  = m(4)*m(6)
        m(21)  = m(3)*m(7)
        m(22)  = m(5)*m(7)
        m(23)  = m(2)*m(3)
        m(24)  = m(4)*m(5)
        m(25)  = m(6)*m(7)
        m(26)  = m(6)*m(8)
        m(27)  = m(7)*m(8)
        m(28)  = m(4)*m(9)
        m(29)  = m(5)*m(9)
        m(30)  = m(2)*m(10)
        m(31)  = m(3)*m(10)
        m(32)  = m(8)*m(9)
        m(33)  = m(8)*m(10)
        m(34)  = m(9)*m(10)
        m(35)  = m(3)*m(20)
        m(36)  = m(2)*m(14)
        m(37)  = m(3)*m(14)
        m(38)  = m(2)*m(16)
        m(39)  = m(3)*m(16)
        m(40)  = m(2)*m(22)
        m(41)  = m(2)*m(20)
        m(42)  = m(3)*m(22)
        m(43)  = m(2)*m(11)
        m(44)  = m(2)*m(18)
        m(45)  = m(2)*m(24)
        m(46)  = m(3)*m(24)
        m(47)  = m(2)*m(13)
        m(48)  = m(4)*m(14)
        m(49)  = m(2)*m(21)
        m(50)  = m(4)*m(22)
        m(51)  = m(2)*m(25)
        m(52)  = m(3)*m(25)
        m(53)  = m(4)*m(25)
        m(54)  = m(5)*m(25)
        m(55)  = m(6)*m(27)
        m(56)  = m(4)*m(29)
        m(57)  = m(2)*m(31)
        m(58)  = m(3)*m(26)
        m(59)  = m(5)*m(26)
        m(60)  = m(2)*m(27)
        m(61)  = m(4)*m(27)
        m(62)  = m(3)*m(28)
        m(63)  = m(2)*m(29)
        m(64)  = m(6)*m(29)
        m(65)  = m(7)*m(28)
        m(66)  = m(4)*m(31)
        m(67)  = m(5)*m(30)
        m(68)  = m(6)*m(31)
        m(69)  = m(7)*m(30)
        m(70)  = m(2)*m(26)
        m(71)  = m(4)*m(26)
        m(72)  = m(3)*m(27)
        m(73)  = m(5)*m(27)
        m(74)  = m(2)*m(28)
        m(75)  = m(3)*m(29)
        m(76)  = m(6)*m(28)
        m(77)  = m(7)*m(29)
        m(78)  = m(4)*m(30)
        m(79)  = m(5)*m(31)
        m(80)  = m(6)*m(30)
        m(81)  = m(7)*m(31)
        m(82)  = m(4)*m(32)
        m(83)  = m(5)*m(32)
        m(84)  = m(6)*m(32)
        m(85)  = m(7)*m(32)
        m(86)  = m(2)*m(33)
        m(87)  = m(3)*m(33)
        m(88)  = m(6)*m(33)
        m(89)  = m(7)*m(33)
        m(90)  = m(2)*m(34)
        m(91)  = m(3)*m(34)
        m(92)  = m(4)*m(34)
        m(93)  = m(5)*m(34)
        m(94)  = m(8)*m(34)
        m(95)  = m(2)*m(37)
        m(96)  = m(3)*m(48)
        m(97)  = m(2)*m(39)
        m(98)  = m(2)*m(50)
        m(99)  = m(3)*m(53)
        m(100) = m(2)*m(54)
        m(101) = m(2)*m(35)
        m(102) = m(2)*m(48)
        m(103) = m(2)*m(42)
        m(104) = m(3)*m(50)
        m(105) = m(2)*m(53)
        m(106) = m(3)*m(54)
        m(107) = m(2)*m(46)
        m(108) = m(2)*m(52)
        m(109) = m(4)*m(54)
        m(110) = m(2)*m(55)
        m(111) = m(3)*m(55)
        m(112) = m(4)*m(55)
        m(113) = m(5)*m(55)
        m(114) = m(2)*m(56)
        m(115) = m(3)*m(56)
        m(116) = m(4)*m(64)
        m(117) = m(4)*m(77)
        m(118) = m(2)*m(66)
        m(119) = m(2)*m(79)
        m(120) = m(2)*m(68)
        m(121) = m(2)*m(81)
        m(122) = m(3)*m(71)
        m(123) = m(2)*m(59)
        m(124) = m(3)*m(61)
        m(125) = m(2)*m(73)
        m(126) = m(3)*m(76)
        m(127) = m(3)*m(64)
        m(128) = m(2)*m(65)
        m(129) = m(2)*m(77)
        m(130) = m(5)*m(80)
        m(131) = m(5)*m(68)
        m(132) = m(4)*m(69)
        m(133) = m(4)*m(81)
        m(134) = m(2)*m(58)
        m(135) = m(4)*m(59)
        m(136) = m(2)*m(72)
        m(137) = m(4)*m(73)
        m(138) = m(2)*m(62)
        m(139) = m(2)*m(75)
        m(140) = m(6)*m(65)
        m(141) = m(6)*m(77)
        m(142) = m(4)*m(67)
        m(143) = m(4)*m(79)
        m(144) = m(6)*m(69)
        m(145) = m(6)*m(81)
        m(146) = m(5)*m(84)
        m(147) = m(4)*m(85)
        m(148) = m(3)*m(88)
        m(149) = m(2)*m(89)
        m(150) = m(3)*m(92)
        m(151) = m(2)*m(93)
        m(152) = m(4)*m(84)
        m(153) = m(5)*m(85)
        m(154) = m(2)*m(88)
        m(155) = m(3)*m(89)
        m(156) = m(2)*m(92)
        m(157) = m(3)*m(93)
        m(158) = m(4)*m(83)
        m(159) = m(6)*m(85)
        m(160) = m(2)*m(87)
        m(161) = m(6)*m(89)
        m(162) = m(2)*m(91)
        m(163) = m(4)*m(93)
        m(164) = m(3)*m(35)
        m(165) = m(5)*m(36)
        m(166) = m(6)*m(37)
        m(167) = m(4)*m(39)
        m(168) = m(2)*m(40)
        m(169) = m(7)*m(38)
    end subroutine evmono

    !> @brief Evaluate evpoly.
    subroutine evpoly( m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), dimension(0:169), intent(in)::m
        real(wp), dimension(0:138), intent(out)::p
        p(0)  = m(0)
        p(1)  = m(1)
        p(2)  = m(2) + m(3) + m(4) + m(5) + m(6) + m(7)
        p(3)  = m(8) + m(9) + m(10)
        p(4)  = p(1)*p(2)
        p(5)  = m(11) + m(12) + m(13) + m(14) + m(15) + m(16)
        p(6)  = m(17) + m(18) + m(19) + m(20) + m(21) + m(22)
        p(7)  = m(23) + m(24) + m(25)
        p(8)  = p(1)*p(3)
        p(9)  = m(26) + m(27) + m(28) + m(29) + m(30) + m(31)
        p(10) = p(2)*p(3) - p(9)
        p(11) = m(32) + m(33) + m(34)
        p(12) = p(1)*p(1)
        p(13) = p(2)*p(2) - p(7) - p(6) - p(5) - p(7) - p(6) - p(5)
        p(14) = p(3)*p(3) - p(11) - p(11)
        p(15) = p(1)*p(5)
        p(16) = p(1)*p(6)
        p(17) = m(35) + m(36) + m(37) + m(38) + m(39) + m(40)
        p(18) = m(41) + m(42)
        p(19) = p(1)*p(7)
        p(20) = m(43) + m(44) + m(45) + m(46) + m(47) + m(48) + m(49) + m(50) + m(51) + m(52) + m(53) + m(54)
        p(21) = p(1)*p(9)
        p(22) = m(55) + m(56) + m(57)
        p(23) = p(1)*p(10)
        p(24) = m(58) + m(59) + m(60) + m(61) + m(62) + m(63) + m(64) + m(65) + m(66) + m(67) + m(68) + m(69)
        p(25) = m(70) + m(71) + m(72) + m(73) + m(74) + m(75) + m(76) + m(77) + m(78) + m(79) + m(80) + m(81)
        p(26) = p(3)*p(5) - p(24)
        p(27) = p(3)*p(6) - p(25)
        p(28) = p(3)*p(7) - p(22)
        p(29) = p(1)*p(11)
        p(30) = m(82) + m(83) + m(84) + m(85) + m(86) + m(87) + m(88) + m(89) + m(90) + m(91) + m(92) + m(93)
        p(31) = m(94)
        p(32) = p(2)*p(11) - p(30)
        p(33) = p(1)*p(4)
        p(34) = p(1)*p(13)
        p(35) = p(2)*p(5) - p(20) - p(17) - p(17)
        p(36) = p(2)*p(6) - p(20) - p(18) - p(17) - p(18) - p(18)
        p(37) = p(2)*p(7) - p(20)
        p(38) = p(1)*p(8)
        p(39) = p(2)*p(9) - p(25) - p(24) - p(22) - p(22)
        p(40) = p(3)*p(13) - p(39)
        p(41) = p(1)*p(14)
        p(42) = p(3)*p(9) - p(30)
        p(43) = p(2)*p(14) - p(42)
        p(44) = p(3)*p(11) - p(31) - p(31) - p(31)
        p(45) = p(1)*p(12)
        p(46) = p(2)*p(13) - p(37) - p(36) - p(35)
        p(47) = p(3)*p(14) - p(44)
        p(48) = p(1)*p(17)
        p(49) = p(1)*p(18)
        p(50) = p(1)*p(20)
        p(51) = m(95) + m(96) + m(97) + m(98) + m(99) + m(100)
        p(52) = m(101) + m(102) + m(103) + m(104) + m(105) + m(106)
        p(53) = m(107) + m(108) + m(109)
        p(54) = p(1)*p(22)
        p(55) = p(1)*p(24)
        p(56) = p(1)*p(25)
        p(57)  = m(110) + m(111) + m(112) + m(113) + m(114) + m(115) + m(116) + m(117) + m(118) + m(119) + m(120) + &
            & m(121)
        p(58) = p(1)*p(26)
        p(59)  = m(122) + m(123) + m(124) + m(125) + m(126) + m(127) + m(128) + m(129) + m(130) + m(131) + m(132) + &
            & m(133)
        p(60) = p(1)*p(27)
        p(61) = p(3)*p(17) - p(59)
        p(62) = p(3)*p(18)
        p(63) = p(1)*p(28)
        p(64)  = m(134) + m(135) + m(136) + m(137) + m(138) + m(139) + m(140) + m(141) + m(142) + m(143) + m(144) + &
            & m(145)
        p(65)  = p(3)*p(20) - p(64) - p(57)
        p(66)  = p(1)*p(30)
        p(67)  = m(146) + m(147) + m(148) + m(149) + m(150) + m(151)
        p(68)  = m(152) + m(153) + m(154) + m(155) + m(156) + m(157)
        p(69)  = m(158) + m(159) + m(160) + m(161) + m(162) + m(163)
        p(70)  = p(1)*p(31)
        p(71)  = p(1)*p(32)
        p(72)  = p(5)*p(11) - p(67)
        p(73)  = p(6)*p(11) - p(68)
        p(74)  = p(31)*p(2)
        p(75)  = p(7)*p(11) - p(69)
        p(76)  = p(1)*p(15)
        p(77)  = p(1)*p(16)
        p(78)  = p(1)*p(19)
        p(79)  = p(1)*p(35)
        p(80)  = m(164) + m(165) + m(166) + m(167) + m(168) + m(169)
        p(81)  = p(1)*p(36)
        p(82)  = p(17)*p(2) - p(52) - p(51) - p(80) - p(51)
        p(83)  = p(2)*p(18) - p(52)
        p(84)  = p(5)*p(6) - p(52) - p(82) - p(52)
        p(85)  = p(1)*p(37)
        p(86)  = p(5)*p(7) - p(51)
        p(87)  = p(6)*p(7) - p(52)
        p(88)  = p(1)*p(21)
        p(89)  = p(1)*p(39)
        p(90)  = p(2)*p(22) - p(57)
        p(91)  = p(1)*p(23)
        p(92)  = p(5)*p(9) - p(59) - p(57)
        p(93)  = p(6)*p(9) - p(62) - p(61) - p(57)
        p(94)  = p(1)*p(40)
        p(95)  = p(2)*p(24) - p(64) - p(59) - p(61) - p(57) -p(92)-p(61)
        p(96)  = p(2)*p(25) - p(64) - p(62) - p(59) - p(57) -p(93)-p(62)
        p(97)  = p(2)*p(26) - p(65) - p(59)
        p(98)  = p(3)*p(36) - p(96) - p(93)
        p(99)  = p(3)*p(37) - p(90)
        p(100) = p(1)*p(29)
        p(101) = p(2)*p(30) - p(73) - p(72) - p(69) - p(68)-p(67)-p(69) - p(68) - p(67)
        p(102) = p(11)*p(13) - p(101)
        p(103) = p(1)*p(42)
        p(104) = p(3)*p(22) - p(69)
        p(105) = p(1)*p(43)
        p(106) = p(3)*p(24) - p(72) - p(67) - p(67)
        p(107) = p(3)*p(25) - p(73) - p(68) - p(68)
        p(108) = p(3)*p(26) - p(72)
        p(109) = p(3)*p(27) - p(73)
        p(110) = p(7)*p(14) - p(104)
        p(111) = p(1)*p(44)
        p(112) = p(9)*p(11) - p(74)
        p(113) = p(3)*p(30) - p(74) - p(112) - p(74)
        p(114) = p(31)*p(3)
        p(115) = p(3)*p(32) - p(74)
        p(116) = p(1)*p(33)
        p(117) = p(1)*p(34)
        p(118) = p(5)*p(5) - p(53) - p(51) - p(80) - p(53) -p(51)-p(80)
        p(119) = p(6)*p(6) - p(53) - p(51) - p(83) - p(53) -p(51)-p(83)
        p(120) = p(7)*p(7) - p(53) - p(53)
        p(121) = p(1)*p(46)
        p(122) = p(5)*p(13) - p(87) - p(82)
        p(123) = p(6)*p(13) - p(86) - p(83) - p(80)
        p(124) = p(7)*p(13) - p(84)
        p(125) = p(1)*p(38)
        p(126) = p(2)*p(39) - p(93) - p(92) - p(90)
        p(127) = p(3)*p(46) - p(126)
        p(128) = p(1)*p(41)
        p(129) = p(3)*p(39) - p(101)
        p(130) = p(13)*p(14) - p(129)
        p(131) = p(11)*p(11) - p(114) - p(114)
        p(132) = p(1)*p(47)
        p(133) = p(3)*p(42) - p(112)
        p(134) = p(2)*p(47) - p(133)
        p(135) = p(11)*p(14) - p(114)
        p(136) = p(1)*p(45)
        p(137) = p(2)*p(46) - p(124) - p(123) - p(122)
        p(138) = p(3)*p(47) - p(135)
    end subroutine evpoly
end module chempotpy_o2h3_pipnn_pes

!> @brief Expose O2H3_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_o2h3_pipnn_pes, only: evaluate
    use chempotpy_kinds,          only: wp
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
