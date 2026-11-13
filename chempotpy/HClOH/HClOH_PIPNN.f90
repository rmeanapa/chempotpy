!> @brief Implement HClOH_PIPNN, a potential energy surface for HClOH.
module chempotpy_hcloh_pipnn_pes
    use chempotpy_hcloh_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_1_count, hidden_2_biases, &
        & hidden_2_count, hidden_hidden_weights, hidden_output_weights, input_hidden_weights, &
        & network_count, ninput, nscale, output_biases, scale_averages, scale_deltas, &
        & vpescut
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the HClOH_PIPNN potential energy surface for HClOH.
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
            call hcloh_pes_interface(tx, v)
        else
            write (*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
    end subroutine evaluate

    !> @brief Evaluate hcloh_pes_interface.
    subroutine hcloh_pes_interface( ct, vpes )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: lxc
        real(wp) :: xct(3,4), ct(3,4), xvec(3,6), xct2(3,4)
        real(wp) :: vpes, vpesa, vpesb, vpesc
        real(wp) :: r13, r23
        vpes      = 0.e0_wp
        vpesa     = 0.e0_wp
        vpesb     = 0.e0_wp
        vpesc     = 0.e0_wp
        xct       = ct
        xvec(:,2) = xct(:,3)-xct(:,1) !  Cl-H1
        xvec(:,4) = xct(:,3)-xct(:,2) !  H2-Cl
        r13       = sqrt(dot_product(xvec(:,2), xvec(:,2)))
        r23       = sqrt(dot_product(xvec(:,4), xvec(:,4)))
        ! check the linked atom order: make sure that
        ! 1st H linked or near Cl
        if( r13 .gt. r23 )then
            xct2(:,1) = xct(:,2)
            xct2(:,2) = xct(:,1)
            xct2(:,3) = xct(:,3)
            xct2(:,4) = xct(:,4)
        else
            xct2 = xct
        endif
        call before_vpes(xct2, lxc)
        if( lxc .eq. -1 )then
            vpes  = vpescut
            vpesa = vpescut
            vpesb = vpescut
            vpesc = vpescut
            return
        endif
        if( lxc .eq. 0 )then
            call hclohnn(xct2, vpes, vpesa, vpesb, vpesc)
            vpes = min(vpes, vpescut)
            return
        endif
    end
    !> @brief Evaluate hclohNN.
    subroutine hclohnn( ct, vpes, vpesa, vpesb, vpesc )
        ! ct in HHClO order and angstrom
        use chempotpy_geometry, only: all_pair_distances
        use chempotpy_kinds,    only: wp
        implicit none
        integer             :: i
        real(wp)            :: rb(6), xbond(6), basis(0:17)
        real(wp)            :: txinput(17)
        real(wp)            :: ct(3,4)
        real(wp)            :: vpesa, vpesb, vpesc, vpes
        real(wp), parameter :: alpha = 2.e0_wp*0.5291772e0_wp
        basis = 0.e0_wp
        call all_pair_distances(transpose(ct), rb)
        xbond(:) = exp(-rb(:)/alpha)
        call bemsav(xbond, basis)
        do i=1, 17
            txinput(i) = basis(i)
        enddo
        call getpota(txinput, vpesa)
        call getpotb(txinput, vpesb)
        call getpotc(txinput, vpesc)
        vpes = (vpesa+vpesb+vpesc)/3.0e0_wp
        if( vpes .lt. -1.5e0_wp ) vpes=5.0e0_wp
        vpes = min(vpes, 5.0e0_wp)
    end subroutine hclohnn

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
    !> @brief Evaluate bemsav.
    subroutine bemsav( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(6), m(0:10), p(0:17)
        call evmono(x, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(6), m(0:10)
        m(0)  = 1.e0_wp
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
        real(wp) :: m(0:10), p(0:17)
        p(0)  = m(0)
        p(1)  = m(1)
        p(2)  = m(2)+m(3)
        p(3)  = m(4)+m(5)
        p(4)  = m(6)
        p(5)  = p(1)*p(2)
        p(6)  = m(7)
        p(7)  = p(1)*p(3)
        p(8)  = m(8)+m(9)
        p(9)  = m(10)
        p(10) = p(2)*p(3)-p(8)
        p(11) = p(1)*p(4)
        p(12) = p(4)*p(2)
        p(13) = p(4)*p(3)
        p(14) = p(1)*p(1)
        p(15) = p(2)*p(2)-p(6)-p(6)
        p(16) = p(3)*p(3)-p(9)-p(9)
        p(17) = p(4)*p(4)
    end subroutine evpoly

    !> @brief Evaluate before_vpes.
    subroutine before_vpes( xct, lxc )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: lxc, i
        real(wp) :: xct(3,4), xct2(3,4), xvec(3,6), bt(6)
        real(wp) :: r12, r13, r14, r23, r24, r34
        ! lxc=-1, unphysical point
        ! lxc=0,  normal
        ! lxc=1,  HO+HCl, using switch function
        ! lxc=2,  HO+HCl, using dipole interaction correction
        lxc       = 0
        xct2      = xct
        xvec(:,1) = xct2(:,2)-xct2(:,1) ! H2-H1
        xvec(:,2) = xct2(:,3)-xct2(:,1) ! Cl-H1
        xvec(:,3) = xct2(:,4)-xct2(:,1) ! O-H1
        xvec(:,4) = xct2(:,3)-xct2(:,2) ! H2-Cl
        xvec(:,5) = xct2(:,4)-xct2(:,2) ! H2-O
        xvec(:,6) = xct2(:,4)-xct2(:,3) ! Cl-O
        r12       = sqrt(dot_product(xvec(:,1), xvec(:,1)))
        r13       = sqrt(dot_product(xvec(:,2), xvec(:,2)))
        r14       = sqrt(dot_product(xvec(:,3), xvec(:,3)))
        r23       = sqrt(dot_product(xvec(:,4), xvec(:,4)))
        r24       = sqrt(dot_product(xvec(:,5), xvec(:,5)))
        r34       = sqrt(dot_product(xvec(:,6), xvec(:,6)))
        bt(:)     = (/r12, r13, r14, r23, r24, r34/)
        ! remove some unphysical points
        do i=1, 6
            if( bt(i) .le. 0.60e0_wp )then
                lxc = -1
                return
            endif
        enddo
        !--> ClO distance, do not consider  H2 + ClO channel
        !--> ClO equilibrium distance,1.596gstrom, x 130% = 2.07 angstrom
        if( r34 .le. 1.30e0_wp .or. r12 .le. 0.60e0_wp )then
            lxc = -1
            return
        endif
        if( minval(bt) .ge. 2.00e0_wp )then
            lxc = -1
            return
        endif
        if( r13 .le. 0.75e0_wp .or. r23 .le. 0.75e0_wp )then
            lxc = -1
            return
        endif
        ! H+O+HCl or HO + H +Cl or H+H+Cl+O
        if( r14 .lt. r24 .and. r24 .gt. 2.0e0_wp )then
            if( r14 .ge. 2.0e0_wp .or. r23 .ge. 2.0e0_wp )then
                lxc = -1
                return
            endif
        endif
        if( r24 .lt. r14 .and. r14 .gt. 2.0e0_wp )then
            if( r24 .ge. 2.0e0_wp .or. r13 .ge. 2.0e0_wp )then
                lxc = -1
                return
            endif
        endif
        ! OH1 and OH2 are both large
        if( r14 .ge. 2.0e0_wp .and. r24 .ge. 2.0e0_wp )then
            lxc = -1
            return
        endif
        ! ClH1 and ClH2 are both small
        if( r13 .le. 1.1e0_wp .and. r23 .le. 1.1e0_wp )then
            lxc = -1
            return
        endif
        if( r23 .lt. r13 .and. r13 .gt. 1.9e0_wp .and. r24 .gt. 1.8e0_wp )then
            if( r14 .ge. 1.8e0_wp .or. r23 .ge. 1.9e0_wp )then
                lxc = -1
                return
            endif
        endif
        if( r13 .lt. r23 .and. r23 .gt. 1.9e0_wp .and. r14 .gt. 1.8e0_wp )then
            if( r13 .ge. 1.8e0_wp .or. r24 .ge. 1.9e0_wp )then
                lxc = -1
                return
            endif
        endif
        lxc = 0
    end subroutine before_vpes
end module chempotpy_hcloh_pipnn_pes

!> @brief Expose HClOH_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_hcloh_pipnn_pes, only: evaluate
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
