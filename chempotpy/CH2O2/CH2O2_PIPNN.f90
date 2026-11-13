!> @brief Implement CH2O2_PIPNN, a potential energy surface for CH2O2.
module chempotpy_ch2o2_pipnn_pes
    use chempotpy_ch2o2_pipnn_parameters, only: activation_function, alpha, hidden_1_biases, &
        & hidden_1_count, hidden_2_biases, hidden_2_count, hidden_hidden_weights, &
        & hidden_output_weights, input_hidden_weights, natom, nbasis, nbond, ninput, nscale, &
        & output_biases, para2, scale_averages, scale_deltas, vpesmin
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the CH2O2_PIPNN potential energy surface for CH2O2.
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
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        if( igrad == 0 )then
            call evvdvdx(tx, v)
        else
            write (*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
    end subroutine evaluate

    !> @brief Evaluate the immutable CH2O2_PIPNN neural network.
    subroutine pot3a1( x, vpot3 )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot3
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
            hidden_1_values(inode1) = tranfun1(hidden_1_values(inode1), activation_function)
        enddo
        do inode1=1, hidden_2_count
            hidden_2_values(inode1) = hidden_2_biases(inode1, 1)
            do inode2=1, hidden_1_count
                hidden_2_values(inode1) = hidden_2_values(inode1) + hidden_1_values(inode2)* &
                    & hidden_hidden_weights(inode2, inode1, 1)
            enddo
            hidden_2_values(inode1) = tranfun1(hidden_2_values(inode1), activation_function)
        enddo
        output_value = output_biases(1, 1)
        do inode2=1, hidden_2_count
            output_value = output_value + hidden_2_values(inode2)*hidden_output_weights(inode2, 1, 1)
        enddo
        vpot3 = output_value*scale_deltas(nscale, 1)+scale_averages(nscale, 1)
    end subroutine pot3a1

    !-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->
    !> @brief Evaluate tranfun1.
    function tranfun1( x, ifunc )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: ifunc
        real(wp) :: tranfun1, x
        !    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
        !    ifunc=2, transfer function is log sigmoid function, 'logsig'
        !    ifunc=3, transfer function is pure linear function, 'purelin'. It is imposed to the output layer by default
        if( ifunc .eq. 1 )then
            tranfun1 = tanh(x)
        elseif( ifunc .eq. 2 )then
            tranfun1 = 1e0_wp/(1e0_wp+exp(-x))
        elseif( ifunc .eq. 3 )then
            tranfun1 = x
        endif
    end function tranfun1

    !> @brief Evaluate LOGSIG1.
    function logsig1( x )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x, logsig1
        logsig1 = 1.e0_wp/(1.e0_wp+exp(-x))
    end function logsig1

    !> @brief Evaluate evvdvdx.
    subroutine evvdvdx( xcart, v )
        use chempotpy_geometry, only: all_pair_distances
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), parameter :: alpha = 1.8e0_wp, para=0.5291772e0_wp
        integer, parameter  :: natom = 5
        integer             :: j, ndriv
        real(wp)            :: v, p(0:101), r(1:10)
        real(wp)            :: xcart(1:3,1:5), m(0:38), xmorse(1:10)
        real(wp)            :: txinput(1:101)
        call all_pair_distances(transpose(xcart), r)
        xmorse = exp(-r/(para*alpha))
        call bemsav1(xmorse, p)
        do j=1, 101
            txinput(j) = p(j)
        enddo
        call pot3a1(txinput, v)
    end subroutine evvdvdx

    !> @brief Evaluate bemsav1.
    subroutine bemsav1( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(1:10)
        real(wp) :: p(0:101)
        real(wp) :: m(0:38)
        call evmono1(x, m)
        call evpoly1(m, p)
    end subroutine bemsav1

    !> @brief Evaluate evmono1.
    subroutine evmono1( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(1:10)
        real(wp) :: m(0:38)
        m(0)  = 1.0e0_wp
        m(1)  = x(10)
        m(2)  = x(9)
        m(3)  = x(8)
        m(4)  = x(7)
        m(5)  = x(4)
        m(6)  = x(6)
        m(7)  = x(5)
        m(8)  = x(3)
        m(9)  = x(2)
        m(10) = x(1)
        m(11) = m(1)*m(2)
        m(12) = m(4)*m(5)
        m(13) = m(2)*m(6)
        m(14) = m(1)*m(7)
        m(15) = m(2)*m(8)
        m(16) = m(1)*m(9)
        m(17) = m(5)*m(6)
        m(18) = m(5)*m(7)
        m(19) = m(4)*m(8)
        m(20) = m(4)*m(9)
        m(21) = m(7)*m(8)
        m(22) = m(6)*m(9)
        m(23) = m(6)*m(8)
        m(24) = m(7)*m(9)
        m(25) = m(6)*m(7)
        m(26) = m(8)*m(9)
        m(27) = m(2)*m(17)
        m(28) = m(1)*m(18)
        m(29) = m(2)*m(19)
        m(30) = m(1)*m(20)
        m(31) = m(2)*m(23)
        m(32) = m(1)*m(24)
        m(33) = m(5)*m(25)
        m(34) = m(4)*m(26)
        m(35) = m(6)*m(21)
        m(36) = m(6)*m(24)
        m(37) = m(6)*m(26)
        m(38) = m(7)*m(26)
    end subroutine evmono1

    !> @brief Evaluate evpoly1.
    subroutine evpoly1( m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(1:10)
        real(wp) :: p(0:101)
        real(wp) :: m(0:38)
        p(0)   = m(0)
        p(1)   = m(1) + m(2)
        p(2)   = m(3)
        p(3)   = m(4) + m(5)
        p(4)   = m(6) + m(7) + m(8) + m(9)
        p(5)   = m(10)
        p(6)   = m(11)
        p(7)   = p(2)*p(1)
        p(8)   = p(1)*p(3)
        p(9)   = p(2)*p(3)
        p(10)  = m(12)
        p(11)  = m(13) + m(14) + m(15) + m(16)
        p(12)  = p(1)*p(4) - p(11)
        p(13)  = p(2)*p(4)
        p(14)  = m(17) + m(18) + m(19) + m(20)
        p(15)  = m(21) + m(22)
        p(16)  = m(23) + m(24)
        p(17)  = p(3)*p(4) - p(14)
        p(18)  = m(25) + m(26)
        p(19)  = p(5)*p(1)
        p(20)  = p(2)*p(5)
        p(21)  = p(5)*p(3)
        p(22)  = p(5)*p(4)
        p(23)  = p(1)*p(1) - p(6) - p(6)
        p(24)  = p(2)*p(2)
        p(25)  = p(3)*p(3) - p(10) - p(10)
        p(26)  = p(4)*p(4) - p(18) - p(16) - p(15) - p(18) - p(16) - p(15)
        p(27)  = p(5)*p(5)
        p(28)  = p(2)*p(6)
        p(29)  = p(6)*p(3)
        p(30)  = p(2)*p(8)
        p(31)  = p(10)*p(1)
        p(32)  = p(2)*p(10)
        p(33)  = p(6)*p(4)
        p(34)  = p(2)*p(11)
        p(35)  = p(2)*p(12)
        p(36)  = m(27) + m(28) + m(29) + m(30)
        p(37)  = p(1)*p(14) - p(36)
        p(38)  = p(2)*p(14)
        p(39)  = p(1)*p(15)
        p(40)  = p(2)*p(15)
        p(41)  = m(31) + m(32)
        p(42)  = p(1)*p(16) - p(41)
        p(43)  = p(2)*p(16)
        p(44)  = p(3)*p(11) - p(36)
        p(45)  = p(1)*p(17) - p(44)
        p(46)  = p(2)*p(17)
        p(47)  = p(10)*p(4)
        p(48)  = p(3)*p(15)
        p(49)  = p(3)*p(16)
        p(50)  = p(1)*p(18)
        p(51)  = p(2)*p(18)
        p(52)  = m(33) + m(34)
        p(53)  = m(35) + m(36) + m(37) + m(38)
        p(54)  = p(3)*p(18) - p(52)
        p(55)  = p(5)*p(6)
        p(56)  = p(2)*p(19)
        p(57)  = p(5)*p(8)
        p(58)  = p(2)*p(21)
        p(59)  = p(5)*p(10)
        p(60)  = p(5)*p(11)
        p(61)  = p(5)*p(12)
        p(62)  = p(2)*p(22)
        p(63)  = p(5)*p(14)
        p(64)  = p(5)*p(15)
        p(65)  = p(5)*p(16)
        p(66)  = p(5)*p(17)
        p(67)  = p(5)*p(18)
        p(68)  = p(6)*p(1)
        p(69)  = p(2)*p(23)
        p(70)  = p(2)*p(7)
        p(71)  = p(3)*p(23)
        p(72)  = p(2)*p(9)
        p(73)  = p(1)*p(25)
        p(74)  = p(2)*p(25)
        p(75)  = p(10)*p(3)
        p(76)  = p(1)*p(11) - p(33)
        p(77)  = p(1)*p(12) - p(33)
        p(78)  = p(2)*p(13)
        p(79)  = p(3)*p(14) - p(47)
        p(80)  = p(3)*p(17) - p(47)
        p(81)  = p(4)*p(11) - p(50) - p(41) - p(39) - p(41)
        p(82)  = p(1)*p(26) - p(81)
        p(83)  = p(2)*p(26)
        p(84)  = p(4)*p(14) - p(52) - p(49) - p(48) - p(52)
        p(85)  = p(4)*p(15) - p(53)
        p(86)  = p(4)*p(16) - p(53)
        p(87)  = p(3)*p(26) - p(84)
        p(88)  = p(4)*p(18) - p(53)
        p(89)  = p(5)*p(23)
        p(90)  = p(2)*p(20)
        p(91)  = p(5)*p(25)
        p(92)  = p(5)*p(26)
        p(93)  = p(5)*p(19)
        p(94)  = p(2)*p(27)
        p(95)  = p(5)*p(21)
        p(96)  = p(5)*p(22)
        p(97)  = p(1)*p(23) - p(68)
        p(98)  = p(2)*p(24)
        p(99)  = p(3)*p(25) - p(75)
        p(100) = p(4)*p(26) - p(88) - p(86) - p(85)
        p(101) = p(5)*p(27)
    end subroutine evpoly1
end module chempotpy_ch2o2_pipnn_pes

!> @brief Expose CH2O2_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_ch2o2_pipnn_pes, only: evaluate
    use chempotpy_kinds,           only: wp
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
