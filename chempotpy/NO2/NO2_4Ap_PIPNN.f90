module chempotpy_no2_4ap_pipnn_pes
    use chempotpy_no2_4ap_pipnn_parameters, only: activation_function, alpha, bohr, hidden_1_biases, &
        & hidden_1_count, hidden_2_biases, hidden_2_count, hidden_hidden_weights, &
        & hidden_output_weights, input_hidden_weights, natom, nbond, ninput, npes, nscale, &
        & output_biases, scale_averages, scale_deltas
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NO2_4Ap_PIPNN potential energy surface for NO2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: v, vpes(3)
        real(wp)              :: tx(3,natoms), tg(3*natoms), tgn(3*natoms,3)
        integer               :: iatom, idir, j, istate
        !initialize
        p = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        tg = 0.e0_wp
        call evvdvdx(tx, v, vpes, tg, tgn, igrad)
        do istate=1, nstates
            p(istate) = v
        enddo
        do istate=1, nstates
            do iatom=1, natoms
                do idir=1, 3
                    j                    = 3*(iatom-1)+idir
                    g(istate,iatom,idir) = tg(j)
                enddo
            enddo
        enddo
        d = 0.e0_wp
    end subroutine evaluate

    !> @brief Evaluate evvdvdx.
    subroutine evvdvdx( xcart, v, vpes, dvdxa, dvdx, ndriv )
        use chempotpy_geometry, only: all_pair_distances, flatten_cartesian_coordinates
        !Subroutine to calculate the average potential energy v and gradient dvdxa
        !Xcart     ==>Cartesian coordinates in angstrom
        !             O1: xcart(1,1), xcart(2,1), xcart(3,1)
        !             O2: xcart(1,2), xcart(2,2), xcart(3,2)
        !             N3: xcart(1,3), xcart(2,3), xcart(3,3)
        !v         ==>Average potential energy(in eV), v=sum(vpes)/npes
        !vpes      ==>Potential energy(in eV) for each PES
        !dvdxa     ==>Average gradient(in eV/ang), dvdxa=sum(dvdx)/npes
        !dvdx      ==>Gradient(in eV/ang) for each PES
        !ndriv     ==>ndriv=1 to compute the gradient, otherwise ndriv=0
        !PES1: totrmse=0.0706, errmax=1.1644, trainrmse=0.0685, validrmse=0.0941, testrmse=0.0795eV
        !PES2: totrmse=0.0881, errmax=1.2347, trainrmse=0.0880, validrmse=0.0875, testrmse=0.0903eV
        !PES3:
        !Average potential energy surface: totormse=0.0720eV
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j, ndriv
        real(wp) :: v, vpes(npes), c(0:ninput), dvdx(1:natom*3,npes)
        real(wp) :: xcart(1:3,1:natom), m(0:4), xmorse(1:nbond), r(1:nbond)
        real(wp) :: dvdr(1:nbond), x(1:natom*3), drdx(1:nbond,1:natom*3)
        real(wp) :: txinput(1:ninput), dvdg(1:ninput,npes), p(0:ninput)
        real(wp) :: dpdr(1:nbond,0:ninput), dvdxa(1:natom*3)
        call flatten_cartesian_coordinates(transpose(xcart), x)
        call all_pair_distances(transpose(xcart), r)
        xmorse = exp(-r/alpha)
        call bemsav(xmorse, m, p)
        do j=1, 12
            txinput(j) = p(j)
        enddo
        call pot3a(txinput, vpes, dvdg, ndriv)
        v = 0e0_wp
        do i=1, npes
            v = v+vpes(i)
        enddo
        v = v/npes
        if( ndriv .eq. 1 )then
            call evdvdr(r, c, m, p, dpdr)
            call evdrdx(xcart, r, x, drdx)
            call evdvdx(dvdg, dpdr, drdx, dvdxa, dvdx)
        endif
        if( ndriv .eq. 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
    end subroutine evvdvdx

    !> @brief Evaluate the immutable three-network NO2 ensemble and optional input derivatives.
    subroutine pot3a( x, vpot3, dvdg, ndriv )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        integer, intent(in)   :: ndriv
        real(wp), intent(out) :: vpot3(npes), dvdg(ninput,npes)
        integer               :: i, j, k, network_index
        real(wp)              :: input_values(ninput), hidden_1_values(hidden_1_count)
        real(wp)              :: hidden_2_values(hidden_2_count), output_value, derivative_sum
        dvdg = 0.e0_wp
        do network_index=1, npes
            do i=1, ninput
                input_values(i) = (x(i)-scale_averages(i, network_index))/scale_deltas(i, network_index)
            enddo
            do j=1, hidden_1_count
                hidden_1_values(j) = hidden_1_biases(j, network_index)
                do i=1, ninput
                    hidden_1_values(j) = hidden_1_values(j) + input_values(i)* input_hidden_weights(i, j, network_index)
                enddo
                hidden_1_values(j) = tranfun(hidden_1_values(j), activation_function)
            enddo
            do k=1, hidden_2_count
                hidden_2_values(k) = hidden_2_biases(k, network_index)
                do j=1, hidden_1_count
                    hidden_2_values(k) = hidden_2_values(k) + hidden_1_values(j)* hidden_hidden_weights(j, k, &
                        & network_index)
                enddo
                hidden_2_values(k) = tranfun(hidden_2_values(k), activation_function)
            enddo
            output_value = output_biases(1, network_index)
            do k=1, hidden_2_count
                output_value = output_value + hidden_2_values(k)* hidden_output_weights(k, 1, network_index)
            enddo
            vpot3(network_index) = output_value*scale_deltas(nscale, network_index) + scale_averages(nscale, &
                & network_index)
            if( ndriv .eq. 1 )then
                do i=1, ninput
                    do k=1, hidden_2_count
                        derivative_sum = 0.e0_wp
                        do j=1, hidden_1_count
                            derivative_sum = derivative_sum + hidden_hidden_weights(j, k, network_index)* &
                                & input_hidden_weights(i, j, network_index)*(1.e0_wp-hidden_1_values(j)**2)
                        enddo
                        dvdg(i,network_index) = dvdg(i,network_index) + hidden_output_weights(k, 1, network_index)* &
                            & derivative_sum* (1.e0_wp-hidden_2_values(k)**2)
                    enddo
                    dvdg(i,network_index) = dvdg(i,network_index)*scale_deltas(nscale, network_index)/ scale_deltas(i, &
                        & network_index)
                enddo
            endif
        enddo
    end subroutine pot3a

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
    end function tranfun

    !> @brief Evaluate LOGSIG.
    function logsig( x )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x, logsig
        logsig = 1.e0_wp/(1.e0_wp+exp(-x))
    end function logsig

    !      real*8 v,vpes(npes),c(0:ninput),dvdx(1:natom*3,npes)
    !      real*8 xcart(1:3,1:natom),m(0:4),xmorse(1:nbond),r(1:nbond)
    !      real*8 dvdr(1:nbond),x(1:natom*3),drdx(1:nbond,1:natom*3)
    !      real*8 txinput(1:ninput),dvdg(1:ninput,npes),p(0:ninput)
    !      real*8 dpdr(1:nbond,0:ninput),dvdxa(1:natom*3)
    !       x((3*(i-1)+1):(3*(i-1)+3))=xcart(:,i)
    !      k=0
    !        k=k+1
    !       r(k)=dot_product(xcart(:,i)-xcart(:,j),
    !     $                  xcart(:,i)-xcart(:,j))
    !        r(k)=dsqrt(r(k))
    !        xmorse(k)=dexp(-r(k)/alpha)
    !      txinput(j)=p(j)
    !      v=0d0
    !       v=v+vpes(i)
    !      v=v/npes
    !> @brief Evaluate bemsav.
    subroutine bemsav( xmorse, m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: xmorse(1:3)
        real(wp) :: p(0:12)
        real(wp) :: m(0:4)
        call evmono(xmorse, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(1:3)
        real(wp) :: m(0:4)
        m(0) = 1.0e0_wp
        m(1) = x(3)
        m(2) = x(2)
        m(3) = x(1)
        m(4) = m(1)*m(2)
    end subroutine evmono

    !> @brief Evaluate evpoly.
    subroutine evpoly( m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: p(0:12)
        real(wp) :: m(0:4)
        p(0)  = m(0)
        p(1)  = m(1) + m(2)
        p(2)  = m(3)
        p(3)  = m(4)
        p(4)  = p(2)*p(1)
        p(5)  = p(1)*p(1) - p(3) - p(3)
        p(6)  = p(2)*p(2)
        p(7)  = p(2)*p(3)
        p(8)  = p(3)*p(1)
        p(9)  = p(2)*p(5)
        p(10) = p(2)*p(4)
        p(11) = p(1)*p(5) - p(8)
        p(12) = p(2)*p(6)
    end subroutine evpoly

    !> @brief Evaluate evdvdx.
    subroutine evdvdx( dvdg, dpdr, drdx, dvdxa, dvdx )
        use chempotpy_geometry, only: transform_internal_gradient_to_cartesian
        use chempotpy_kinds,    only: wp
        implicit none
        integer  :: i, j, k
        real(wp) :: dvdx(1:natom*3,npes), xcart(1:3,1:natom)
        real(wp) :: drdx(1:nbond,1:natom*3), x(1:ninput)
        real(wp) :: dvdr(1:nbond), dpdr(1:nbond,0:ninput), dvdp(0:ninput)
        real(wp) :: dgdx(1:natom*3,1:ninput), dgdr(1:nbond,1:ninput)
        real(wp) :: dvdxa(1:natom*3), dvdg(1:ninput,npes), rtmp, basis_cartesian(ninput,natom*3)
        call transform_internal_gradient_to_cartesian(transpose(dpdr(:,1:ninput)), drdx, basis_cartesian)
        dgdx = transpose(basis_cartesian)
        dvdx = 0e0_wp;dvdxa=0
        do k=1, npes
            do i=1, 9
                do j=1, 12
                    dvdx(i,k) = dvdx(i,k) + dvdg(j,k)*dgdx(i,j)
                enddo
            enddo
        enddo
        do i=1, 9
            do k=1, npes
                dvdxa(i) = dvdxa(i)+dvdx(i,k)
            enddo
        enddo
        dvdxa = dvdxa/real(npes, wp)
    end subroutine evdvdx

    !> @brief Evaluate evdvdr.
    subroutine evdvdr( r, c, m, p, dpdr )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j
        real(wp) :: r(1:nbond), dmsdr(1:nbond,1:nbond), dmdr(1:nbond,0:4)
        real(wp) :: dvdr(1:nbond), dpdr(1:nbond,0:ninput), c(0:ninput)
        real(wp) :: m(0:4), p(0:ninput)
        dvdr(:) = 0e0_wp
        call evdmsdr(r, dmsdr)
        call evdmdr(m, dmsdr, dmdr)
        call evdpdr(p, dmdr, dpdr)
    end subroutine evdvdr

    !> @brief Evaluate evdrdx.
    subroutine evdrdx( xcart, r, x, drdx )
        use chempotpy_geometry, only: ordered_pair_distance_jacobian
        use chempotpy_kinds,    only: wp
        implicit none
        integer, parameter    :: first_atoms(3)  = [1, 1, 2]
        integer, parameter    :: second_atoms(3) = [2, 3, 3]
        real(wp), intent(in)  :: xcart(1:3,1:natom), r(1:nbond), x(1:natom*3)
        real(wp), intent(out) :: drdx(1:nbond,1:natom*3)
        call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
    end subroutine evdrdx

    !> @brief Evaluate evdmsdr.
    subroutine evdmsdr( r, dmsdr )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j
        real(wp) :: dmsdr(1:nbond,1:nbond), r(1:nbond)
        dmsdr(:,:) = 0e0_wp
        !Morse term ms = exp(-r/alpha)
        ! dmsdr(i,j)=0  i!=j
        ! dmsdr(i,i)= -(1/alpha)*Exp(-r(i)/alpha)
        do i=1, nbond
            dmsdr(i,i) = -(1e0_wp/alpha)*exp(-r(i)/alpha)
        enddo
    end subroutine evdmsdr

    !> @brief Evaluate evdmdr.
    subroutine evdmdr( m, dmsdr, dmdr )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j
        real(wp) :: dmsdr(1:nbond,1:nbond), dmdr(1:nbond,0:4), m(0:4)
        do i=1, nbond
            dmdr(i,0) = 0.e0_wp
            dmdr(i,1) = dmsdr(i,3)
            dmdr(i,2) = dmsdr(i,2)
            dmdr(i,3) = dmsdr(i,1)
            dmdr(i,4) = dmdr(i,1)*m(2)+m(1)*dmdr(i,2)
        enddo
    end subroutine evdmdr

    !> @brief Evaluate evdpdr.
    subroutine evdpdr( p, dmdr, dpdr )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j
        real(wp) :: dmdr(1:nbond,0:4), dpdr(1:nbond,0:ninput)
        real(wp) :: p(0:ninput)
        do i=1, nbond
            dpdr(i,0)  = dmdr(i,0)
            dpdr(i,1)  = dmdr(i,1)+dmdr(i,2)
            dpdr(i,2)  = dmdr(i,3)
            dpdr(i,3)  = dmdr(i,4)
            dpdr(i,4)  = dpdr(i,2)*p(1)+p(2)*dpdr(i,1)
            dpdr(i,5)  = dpdr(i,1)*p(1)+p(1)*dpdr(i,1)-dpdr(i,3)-dpdr(i,3)
            dpdr(i,6)  = dpdr(i,2)*p(2)+p(2)*dpdr(i,2)
            dpdr(i,7)  = dpdr(i,2)*p(3)+p(2)*dpdr(i,3)
            dpdr(i,8)  = dpdr(i,3)*p(1)+p(3)*dpdr(i,1)
            dpdr(i,9)  = dpdr(i,2)*p(5)+p(2)*dpdr(i,5)
            dpdr(i,10) = dpdr(i,2)*p(4)+p(2)*dpdr(i,4)
            dpdr(i,11) = dpdr(i,1)*p(5)+p(1)*dpdr(i,5)-dpdr(i,8)
            dpdr(i,12) = dpdr(i,2)*p(6)+p(2)*dpdr(i,6)
        enddo
    end subroutine evdpdr
end module chempotpy_no2_4ap_pipnn_pes

!> @brief Expose NO2_4Ap_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_no2_4ap_pipnn_pes, only: evaluate
    use chempotpy_kinds,             only: wp
    implicit none
    integer, parameter              :: natoms  = 3
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
