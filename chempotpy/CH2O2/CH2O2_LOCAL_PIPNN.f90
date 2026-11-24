! System:                        CH2O2
! Functional form:               permutationally invariant polynomial neural network
! Common name:                   CH2O2 local PIP-NN
! Reference:                     J. Li et al., J. Phys. Chem. Lett. 5, 2364 (2014)
! Number of derivatives:         0
! Number of bodies:              5
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; neural-network parameters are embedded
! Notes:                         Local surface around the CH2OO Criegee-intermediate minimum.

!> @brief Evaluate the local CH2OO PIP-NN surface of Li and co-workers.
module chempotpy_ch2o2_local_pipnn_pes
    use chempotpy_ch2o2_local_pipnn_parameters, only: activation_function, alpha, hidden_1_biases, &
        & hidden_1_count, hidden_2_biases, hidden_2_count, hidden_hidden_weights, &
        & hidden_output_weights, input_hidden_weights, network_count, ninput, nscale, &
        & output_biases, scale_averages, scale_deltas, vpescut
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the CH2O2_LOCAL_PIPNN potential energy surface for CH2O2.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 5
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: v
        integer               :: istate
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        if( igrad == 0 )then
            call evaluate_energy(x, v)
        else
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives are unavailable')
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
    end subroutine evaluate

    !> @brief Evaluate one immutable network from the released three-model ensemble.
    subroutine evaluate_network( x, network_index, value )
        use chempotpy_neural_network, only: evaluate_two_hidden_network
        use chempotpy_kinds,          only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        integer, intent(in)   :: network_index
        real(wp), intent(out) :: value
        call evaluate_two_hidden_network(x, scale_averages(1:ninput, network_index), &
            & scale_deltas(1:ninput, network_index), input_hidden_weights(:, :, network_index), &
            & hidden_1_biases(:, network_index), hidden_hidden_weights(:, :, network_index), &
            & hidden_2_biases(:, network_index), hidden_output_weights(:, 1, network_index), &
            & output_biases(1, network_index), activation_function, scale_averages(nscale, network_index), &
            & scale_deltas(nscale, network_index), value)
    end subroutine evaluate_network

    !> @brief Evaluate the released local CH2OO ensemble in electronvolts.
    subroutine evaluate_energy( xcart, v )
        use chempotpy_geometry, only: all_pair_distances
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp), intent(in)  :: xcart(5,3)
        real(wp), intent(out) :: v
        integer               :: j, network_index
        real(wp)              :: network_values(network_count), p(0:101), r(10), xmorse(10), txinput(101)
        call all_pair_distances(xcart, r)
        xmorse = exp(-r/alpha)
        call bemsav1(xmorse, p)
        do j=1, 101
            txinput(j) = p(j)
        enddo
        do network_index=1, network_count
            call evaluate_network(txinput, network_index, network_values(network_index))
        enddo
        v = min(sum(network_values)/real(network_count, wp), vpescut)
    end subroutine evaluate_energy

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
end module chempotpy_ch2o2_local_pipnn_pes

!> @brief Expose CH2O2_LOCAL_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_ch2o2_local_pipnn_pes, only: evaluate
    use chempotpy_kinds,                 only: wp
    implicit none
    integer, parameter    :: natoms  = 5
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
