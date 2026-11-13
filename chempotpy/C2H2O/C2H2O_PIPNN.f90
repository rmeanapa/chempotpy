!> @brief Implement C2H2O_PIPNN, a potential energy surface for C2H2O.
module chempotpy_c2h2o_pipnn_pes
    use chempotpy_c2h2o_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_1_count, hidden_2_biases, &
        & hidden_2_count, hidden_hidden_weights, hidden_output_weights, input_hidden_weights, &
        & network_count, ninput, nscale, output_biases, scale_averages, scale_deltas, &
        & vpescut
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the C2H2O_PIPNN potential energy surface for C2H2O.
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
            call c2h2o_pes_interface(tx, v)
        else
            write (*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
    end subroutine evaluate

    !       q in OCCHH order and angstrom
    !       vpes in eV
    !> @brief Evaluate c2h2o_pes_interface.
    subroutine c2h2o_pes_interface( q, vpes )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: lxc
        real(wp) :: q(3,5)
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
            call c2h2onn(q, vpes, vpesa, vpesb, vpesc)
            if( vpes .lt. -1.5 .or. vpes .gt. vpescut )then
                vpes  = vpescut
                vpesa = vpescut
                vpesb = vpescut
                vpesc = vpescut
            endif
            return
        endif
    end
    !> @brief Evaluate c2h2oNN.
    subroutine c2h2onn( ct, vpes, vpesa, vpesb, vpesc )
        use chempotpy_geometry, only: all_pair_distances
        use chempotpy_kinds,    only: wp
        implicit none
        integer  :: i
        real(wp) :: rb(10), xbond(10), basis(0:101)
        real(wp) :: txinput(101)
        ! ct in OCCHH order and angstrom
        real(wp)            :: ct(3,5), pair_distances(10)
        real(wp)            :: vpesa, vpesb, vpesc, vpes
        real(wp), parameter :: alpha = 2.e0_wp ! bohr
        ! bohr to angstrom
        real(wp), parameter :: btoa = 0.529177249e0_wp
        real(wp), parameter :: pi   = acos(-1.e0_wp)
        basis = 0.e0_wp
        call all_pair_distances(transpose(ct), pair_distances)
        rb       = pair_distances([5,6,7,1,8,9,2,10,3,4])
        xbond(:) = exp(-rb(:)/(alpha*btoa))
        call bemsav(xbond, basis)
        do i=1, 101
            txinput(i) = basis(i)
        enddo
        call getpota(txinput, vpesa)
        call getpotb(txinput, vpesb)
        call getpotc(txinput, vpesc)
        vpes = (vpesa+vpesb+vpesc)/3.0e0_wp
    end subroutine c2h2onn

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
        real(wp) :: xct(3,5), xc(3,5), xvec(3,10), rb(10)
        real(wp) :: rc2c3, rc2h4, rc2h5, rc2o1, rc3h4
        real(wp) :: rc3h5, rc3o1, rh4h5, rh4o1, rh5o1
        real(wp) :: min_rch, xx, yy, zz, tt
        ! lxc=-1,  unphysical point
        ! lxc=0,   normal
        lxc        = 0
        xc         = xct
        xvec(:,1)  = xc(:,3)-xc(:,2)  !  C3->C2
        xvec(:,2)  = xc(:,4)-xc(:,2)  !  H4->C2
        xvec(:,3)  = xc(:,5)-xc(:,2)  !  H5->C2
        xvec(:,4)  = xc(:,1)-xc(:,2)  !  O1->C2
        xvec(:,5)  = xc(:,4)-xc(:,3)  !  H4->C3
        xvec(:,6)  = xc(:,5)-xc(:,3)  !  H5->C3
        xvec(:,7)  = xc(:,1)-xc(:,3)  !  O1->C3
        xvec(:,8)  = xc(:,5)-xc(:,4)  !  H5->H4
        xvec(:,9)  = xc(:,1)-xc(:,4)  !  O1->H4
        xvec(:,10) = xc(:,1)-xc(:,5) !  O1->H5
        rc2c3      = sqrt(dot_product(xvec(:,1), xvec(:,1)))    !  C3->C2
        rc2h4      = sqrt(dot_product(xvec(:,2), xvec(:,2)))    !  H4->C2
        rc2h5      = sqrt(dot_product(xvec(:,3), xvec(:,3)))    !  H5->C2
        rc2o1      = sqrt(dot_product(xvec(:,4), xvec(:,4)))    !  O1->C2
        rc3h4      = sqrt(dot_product(xvec(:,5), xvec(:,5)))    !  H4->C3
        rc3h5      = sqrt(dot_product(xvec(:,6), xvec(:,6)))    !  H5->C3
        rc3o1      = sqrt(dot_product(xvec(:,7), xvec(:,7)))    !  O1->C3
        rh4h5      = sqrt(dot_product(xvec(:,8), xvec(:,8)))    !  H5->H4
        rh4o1      = sqrt(dot_product(xvec(:,9), xvec(:,9)))    !  O1->H4
        rh5o1      = sqrt(dot_product(xvec(:,10), xvec(:,10)))  !  O1->H5
        rb(:)      = (/rc2c3, rc2h4, rc2h5, rc2o1, rc3h4, rc3h5, rc3o1, rh4h5, rh4o1, rh5o1/)
        min_rch    = min(rc2h4, rc2h5, rc3h4, rc3h5)
        ! remove some unphysical points
        do i=1, 10
            if( rb(i) .le. 0.7e0_wp )then
                lxc = -1
                return
            endif
        enddo
        ! CH single bond equilibrium distance,1.1199 angstrom,
        !    1.1199*150% = 1.6798 angstrom
        xx = 1.68e0_wp
        ! CC single bond equilibrium distance,1.5360 angstrom,
        !    1.5360*150% = 2.3040 angstrom
        yy = 2.30e0_wp
        ! CO single bond equilibrium distance,1.4270 angstrom,
        !    1.4270*150% = 2.1405 angstrom
        zz = 2.14e0_wp
        ! HO single bond equilibrium distance,0.9560 angstrom,
        !    0.9560*150% = 1.4340 angstrom
        tt = 1.43e0_wp
        if( min_rch .gt. xx )then
            lxc = -1
            return
        else
            if( rc2c3 .gt. yy )then
                if((rc2h4.le.xx .and. rc2h5.le.xx .and. rc3o1.le.zz).or. &
                    &          (rc3h4.le.xx .and. rc3h5.le.xx .and. rc2o1.le.zz))then
                lxc = 0
                return
            else
                lxc = -1
                return
            endif
        else
            if( (rc2h4 .le. xx .and.rc3h5 .le. xx) .or. (rc2h5 .le. xx .and.rc3h4 .le. xx) )then
                lxc = 0
                return
            elseif(((rc2h4.le.xx.or.rc2h5.le.xx).and.rc3o1.le.yy ).or.&
                    &          ((rc3h4.le.xx.or.rc3h5.le.xx) .and. rc2o1.le.yy ))then
                lxc = 0
                return
            elseif(((rc2h4.le.xx.or.rc3h4.le.xx).and.rh5o1.le.tt ).or. &
                    &          ((rc2h5.le.xx.or.rc3h5.le.xx) .and. rh4o1.le.tt ))then
                lxc = 0
                return
            elseif((rc2h4.le.xx.and.rc2h5.le.xx.and.rc3o1.le.yy ).or. &
                    &          (rc3h4.le.xx.and.rc3h5.le.xx.and.rc2o1.le.yy ))then
                lxc = 0
                return
            else
                lxc = -1
                !     $                  HCCO+H/HCC+OH/H2CCO channel"
                return
            endif
        endif
        endif
    end subroutine before_vpes

    !------->
    !------->This is the subroutine for PIPs up to 3rd order
    !------->of A2B2C type molecule.
    !------->
    !> @brief Evaluate bemsav.
    subroutine bemsav( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:10)
        real(wp), intent(out) :: p(0:101)
        real(wp)              :: m(0:38)
        call evmono(x, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:10)
        real(wp), intent(out) :: m(0:38)
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
    end subroutine evmono

    !> @brief Evaluate evpoly.
    subroutine evpoly( m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: m(0:38)
        real(wp), intent(out) :: p(0:101)
        p(0)   = m(0)
        p(1)   = m(1)+m(2)
        p(2)   = m(3)
        p(3)   = m(4)+m(5)
        p(4)   = m(6)+m(7)+m(8)+m(9)
        p(5)   = m(10)
        p(6)   = m(11)
        p(7)   = p(2)*p(1)
        p(8)   = p(1)*p(3)
        p(9)   = p(2)*p(3)
        p(10)  = m(12)
        p(11)  = m(13)+m(14)+m(15)+m(16)
        p(12)  = p(1)*p(4)-p(11)
        p(13)  = p(2)*p(4)
        p(14)  = m(17)+m(18)+m(19)+m(20)
        p(15)  = m(21)+m(22)
        p(16)  = m(23)+m(24)
        p(17)  = p(3)*p(4)-p(14)
        p(18)  = m(25)+m(26)
        p(19)  = p(5)*p(1)
        p(20)  = p(2)*p(5)
        p(21)  = p(5)*p(3)
        p(22)  = p(5)*p(4)
        p(23)  = p(1)*p(1)-p(6)-p(6)
        p(24)  = p(2)*p(2)
        p(25)  = p(3)*p(3)-p(10)-p(10)
        p(26)  = p(4)*p(4)-p(18)-p(16)-p(15)-p(18)-p(16)-p(15)
        p(27)  = p(5)*p(5)
        p(28)  = p(2)*p(6)
        p(29)  = p(6)*p(3)
        p(30)  = p(2)*p(8)
        p(31)  = p(10)*p(1)
        p(32)  = p(2)*p(10)
        p(33)  = p(6)*p(4)
        p(34)  = p(2)*p(11)
        p(35)  = p(2)*p(12)
        p(36)  = m(27)+m(28)+m(29)+m(30)
        p(37)  = p(1)*p(14)-p(36)
        p(38)  = p(2)*p(14)
        p(39)  = p(1)*p(15)
        p(40)  = p(2)*p(15)
        p(41)  = m(31)+m(32)
        p(42)  = p(1)*p(16)-p(41)
        p(43)  = p(2)*p(16)
        p(44)  = p(3)*p(11)-p(36)
        p(45)  = p(1)*p(17)-p(44)
        p(46)  = p(2)*p(17)
        p(47)  = p(10)*p(4)
        p(48)  = p(3)*p(15)
        p(49)  = p(3)*p(16)
        p(50)  = p(1)*p(18)
        p(51)  = p(2)*p(18)
        p(52)  = m(33)+m(34)
        p(53)  = m(35)+m(36)+m(37)+m(38)
        p(54)  = p(3)*p(18)-p(52)
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
        p(76)  = p(1)*p(11)-p(33)
        p(77)  = p(1)*p(12)-p(33)
        p(78)  = p(2)*p(13)
        p(79)  = p(3)*p(14)-p(47)
        p(80)  = p(3)*p(17)-p(47)
        p(81)  = p(4)*p(11)-p(50)-p(41)-p(39)-p(41)
        p(82)  = p(1)*p(26)-p(81)
        p(83)  = p(2)*p(26)
        p(84)  = p(4)*p(14)-p(52)-p(49)-p(48)-p(52)
        p(85)  = p(4)*p(15)-p(53)
        p(86)  = p(4)*p(16)-p(53)
        p(87)  = p(3)*p(26)-p(84)
        p(88)  = p(4)*p(18)-p(53)
        p(89)  = p(5)*p(23)
        p(90)  = p(2)*p(20)
        p(91)  = p(5)*p(25)
        p(92)  = p(5)*p(26)
        p(93)  = p(5)*p(19)
        p(94)  = p(2)*p(27)
        p(95)  = p(5)*p(21)
        p(96)  = p(5)*p(22)
        p(97)  = p(1)*p(23)-p(68)
        p(98)  = p(2)*p(24)
        p(99)  = p(3)*p(25)-p(75)
        p(100) = p(4)*p(26)-p(88)-p(86)-p(85)
        p(101) = p(5)*p(27)
    end subroutine evpoly
end module chempotpy_c2h2o_pipnn_pes

!> @brief Expose C2H2O_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_c2h2o_pipnn_pes, only: evaluate
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
