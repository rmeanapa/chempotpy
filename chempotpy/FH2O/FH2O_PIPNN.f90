!-->  program to get potential energy for a given geometry after NN
!fitting
!-->  global variables are declared in this module
!> @brief Implement FH2O_PIPNN, a potential energy surface for FH2O.
module chempotpy_fh2o_pipnn_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the FH2O_PIPNN potential energy surface for FH2O.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_fh2o_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_2_biases, hidden_1_count, hidden_2_count, hidden_hidden_weights, &
            & hidden_output_weights, input_hidden_weights, network_count, ninput, nlayer, nodemax, &
            & nodes, noutput, nscale, output_biases, scale_averages, scale_deltas, &
            & alpha, bohr, nbasis, npbasissoc, nrbasissoc, p2soccoef, &
            & pi, r1soccoef, radian, vpesmin
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 4
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: v, va, vb, vc, socg, vgsoc
        double precision              :: tx(3,natoms)
        integer                       :: iatom, idir, j, istate
        !initialize
        v = 0.d0; g = 0.d0; d = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        if( igrad == 0 )then
            call fh2onn(tx, v, va, vb, vc)
        else
            write (*, *) 'Only energy and gradient are available'
        endif
        v = v  !/23.0609
        do istate=1, nstates
            p(istate) = v
        enddo
        g = 0.d0
        d = 0.d0
    end subroutine evaluate

    !> @brief Evaluate fh2o_soc_nn.
    subroutine fh2o_soc_nn( xt, vg, socg, vgsoc )
        use chempotpy_geometry,              only: all_pair_distances
        use chempotpy_fh2o_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_2_biases, hidden_1_count, hidden_2_count, hidden_hidden_weights, &
            & hidden_output_weights, input_hidden_weights, network_count, ninput, nlayer, nodemax, &
            & nodes, noutput, nscale, output_biases, scale_averages, scale_deltas, &
            & alpha, bohr, nbasis, npbasissoc, nrbasissoc, p2soccoef, &
            & pi, r1soccoef, radian, vpesmin
        implicit none
        double precision :: xt(3,4), xt2(3,4), xn(3,4), rint(6), rt(6)
        double precision :: tr1basis(nrbasissoc), tr2basis(nrbasissoc)
        double precision :: rcom, the1, the2, phi
        double precision :: tp2basis(npbasissoc)
        integer          :: n, k, l, m, j, nrcom
        double precision :: vg, vgsoc, socg, socs, vpesa, vpesb, vpesc
        double precision :: socrg, socrs, socpg, socps
        double precision :: r1switch, r2switch, p2switch, switch
        call fh2onn(xt, vg, vpesa, vpesb, vpesc) !xt in HHFO and angstrom
        vg = vg*23.0605d0
        vg = min(vg, 100.0d0)
        if( vg .ge. 100.0d0 )then
            socg  = 0.0d0
            vgsoc = vg+socg/349.75d0
            return
        endif
        xt2 = xt
        call orderxnt(xt2, xn) !xn 1st H closer to F
        call all_pair_distances(transpose(xn), rt)
        xt2(:,1) = xt(:,1)
        xt2(:,2) = xt(:,4)
        xt2(:,3) = xt(:,2)
        xt2(:,4) = xt(:,3)
        call xyztojcb3(xt2, rint) !-->xt2 in HOHF order,angstrom, rt in angstrom and radians
        rcom = rint(3)
        the2 = rint(5)
        phi  = rint(6)
        call soc_react_ground(rcom, the2, phi, socrg)
        xt2(:,1) = xn(:,2)
        xt2(:,2) = xn(:,4)
        xt2(:,3) = xn(:,1)
        xt2(:,4) = xn(:,3)
        call cart_dijacobi(xt2, rint) !xt2 in HOHF order,angstrom, rt in angstrom and radians
        rcom = rint(3)
        the1 = rint(4)
        the2 = rint(5)
        phi  = rint(6)
        call soc_product_ground(rcom, the1, the2, phi, socpg)
        switch = (1.d0-dtanh(20.0d0*(rt(3)-1.6d0)))/2.0d0
        socg   = socrg*switch+socpg*(1.0d0-switch)
        vgsoc  = vg+socg/349.75d0
        vgsoc  = min(100.0d0, vgsoc)
    end subroutine fh2o_soc_nn

    !> @brief Evaluate orderxnt.
    subroutine orderxnt( x1, x2 )
        implicit none
        double precision :: x1(3,4), x2(3,4), xvec(3,2), r1, r2, xt(3,4)
        xt        = x1
        xvec(:,1) = xt(:,1)-xt(:,3)
        xvec(:,2) = xt(:,2)-xt(:,3)
        r1        = dsqrt(dot_product(xvec(:,1), xvec(:,1)))
        r2        = dsqrt(dot_product(xvec(:,2), xvec(:,2)))
        if( r1 .gt. r2 )then
            x2(:,1) = xt(:,2)  !-->1st H closer to F
            x2(:,2) = xt(:,1)
            x2(:,3) = xt(:,3)
            x2(:,4) = xt(:,4)
        else
            x2 = xt
        endif
    end subroutine orderxnt

    !> @brief Evaluate r1function.
    double precision function r1function( r, n )
        implicit none
        double precision, parameter :: r0 = 1.8d0, c=0.45d0
        double precision            :: r
        integer                     :: n
        if( n .eq. 0 )then
            r1function = 1.0d0
            return
        else
            r1function = dtanh(dble(n)*c*(r-r0))
            return
        endif
    end function r1function

    !> @brief Evaluate r2function.
    double precision function r2function( r, n )
        implicit none
        double precision, parameter :: r0 = 2.0d0, c=0.42d0
        double precision            :: r
        integer                     :: n
        if( n .eq. 0 )then
            r2function = 1.0d0
            return
        else
            r2function = dtanh(dble(n)*c*(r-r0))
            return
        endif
    end function r2function

    !> @brief Evaluate p2function.
    double precision function p2function( r, n )
        implicit none
        double precision, parameter :: r0 = 2.3d0, c=0.40d0
        double precision            :: r
        integer                     :: n
        if( n .eq. 0 )then
            p2function = 1.0d0
            return
        else
            p2function = dtanh(dble(n)*c*(r-r0))
            return
        endif
    end function p2function

    !> @brief Evaluate soc_react_ground.
    subroutine soc_react_ground( rcom, t2, phi, socrg )
        use chempotpy_fh2o_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_2_biases, hidden_1_count, hidden_2_count, hidden_hidden_weights, &
            & hidden_output_weights, input_hidden_weights, network_count, ninput, nlayer, nodemax, &
            & nodes, noutput, nscale, output_biases, scale_averages, scale_deltas, &
            & alpha, bohr, nbasis, npbasissoc, nrbasissoc, p2soccoef, &
            & pi, r1soccoef, radian, vpesmin
        implicit none
        double precision :: rcom, t2, phi, socrg
        double precision :: tr1basis(nrbasissoc)
        integer          :: j, n, k, l, nrcom
        tr1basis = 0.0d0
        j        = 0
        do n=0, 6
            do k=0, 2
                do l=0, 6
                    j     = j+1
                    nrcom = n
                    if( l .le. 3 )then
                        tr1basis(j) = r1function(rcom, nrcom)*dcos(dble(k)*phi)* dcos(dble(l)*t2)
                    else
                        tr1basis(j) = r1function(rcom, nrcom)*dcos(dble(k)*phi)* dsin(dble(l-3)*t2)
                    endif
                enddo
            enddo
        enddo
        if( j .ne. nrbasissoc )then
            write(*, *)'error for nrbasissoc'
            stop
        endif
        socrg = dot_product(tr1basis(:), r1soccoef(:))
    end subroutine soc_react_ground

    !> @brief Evaluate soc_product_ground.
    subroutine soc_product_ground( rcom, t1, t2, phi, socpg )
        use chempotpy_fh2o_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_2_biases, hidden_1_count, hidden_2_count, hidden_hidden_weights, &
            & hidden_output_weights, input_hidden_weights, network_count, ninput, nlayer, nodemax, &
            & nodes, noutput, nscale, output_biases, scale_averages, scale_deltas, &
            & alpha, bohr, nbasis, npbasissoc, nrbasissoc, p2soccoef, &
            & pi, r1soccoef, radian, vpesmin
        implicit none
        double precision :: rcom, t1, t2, phi, socpg, socps
        double precision :: tp2basis(npbasissoc)
        integer          :: j, n, k, l, m, nrcom
        tp2basis = 0.0d0
        j        = 0
        do n=0, 5
            do k=0, 2
                do l=0, 3
                    do m=0, 3
                        j           = j+1
                        nrcom       = n
                        tp2basis(j) = p2function(rcom, nrcom)*dcos(dble(k)*phi)* dcos(dble(l)*t1)*dcos(dble(m)*t2)
                    enddo
                enddo
            enddo
        enddo
        if( j .ne. npbasissoc )then
            write(*, *)'error for nrbasissoc'
            stop
        endif
        socps = dot_product(tp2basis(:), p2soccoef(:))
        socpg = -1.0d0*socps
    end subroutine soc_product_ground

    !> @brief Evaluate fh2oNN.
    subroutine fh2onn( ct, vpes, vpesa, vpesb, vpesc ) !ct in HHFO order and angstrom
        use chempotpy_geometry,              only: all_pair_distances
        use chempotpy_fh2o_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_2_biases, hidden_1_count, hidden_2_count, hidden_hidden_weights, &
            & hidden_output_weights, input_hidden_weights, network_count, ninput, nlayer, nodemax, &
            & nodes, noutput, nscale, output_biases, scale_averages, scale_deltas, &
            & alpha, bohr, nbasis, npbasissoc, nrbasissoc, p2soccoef, &
            & pi, r1soccoef, radian, vpesmin
        implicit none
        integer, parameter :: ndim = 6
        integer            :: j
        double precision   :: rb(ndim), xbond(ndim), basis(0:nbasis-1), tmp1
        double precision   :: txinput(1:nbasis-1)
        double precision   :: ct(3,4), vpes
        double precision   :: vpesa, vpesb, vpesc, r12, r13, r14, r23, r24, r34
        basis = 0.d0
        call all_pair_distances(transpose(ct), rb)
        r12 = rb(1);r13=rb(2);r14=rb(3);r23=rb(4);r24=rb(5);r34=rb(6)
        if( minval(rb(:)) .le. 0.67d0.or.minval(rb(:)) .gt. 1.45d0 )then
            vpes = 5.0d0
            return
        endif
        if( r34 .le. 1.5d0.or.r12 .le. 0.8d0.or.minval(rb) .eq. r12 )then
            vpes = 5.0d0
            return
        endif
        if( r14 .lt. r24.and.r24 .gt. 1.9d0 )then
            if( r14 .ge. 1.8d0.or.r23 .ge. 1.8d0 )then
                vpes = 5.0d0
                return
            endif
        endif
        if( r24 .lt. r14.and.r14 .gt. 1.9d0 )then
            if( r24 .ge. 1.8d0.or.r13 .ge. 1.8d0 )then
                vpes = 5.0d0
                return
            endif
        endif
        if( r23 .lt. r13.and.r13 .gt. 1.9d0.and.r24 .gt. 1.7d0 )then
            if( r14 .ge. 1.8d0.or.r23 .ge. 1.8d0 )then
                vpes = 5.0d0
                return
            endif
        endif
        if( r13 .lt. r23.and.r23 .gt. 1.9d0.and.r14 .gt. 1.7d0 )then
            if( r13 .ge. 1.8d0.or.r24 .ge. 1.8d0 )then
                vpes = 5.0d0
                return
            endif
        endif
        xbond(:) = dexp(-rb(:)/alpha)
        !      do j=1,nbasis
        !       tmp1=0.0d0
        !       do k=1,nterm(j)
        !        tmp1=tmp1+
        !    $  (xbond(1)**nindex(j,k,1))*(xbond(2)**nindex(j,k,2))
        !    & *(xbond(3)**nindex(j,k,3))*(xbond(4)**nindex(j,k,4))
        !    & *(xbond(5)**nindex(j,k,5))*(xbond(6)**nindex(j,k,6))
        !       enddo
        !       basis(j)=basis(j)+tmp1
        !      enddo
        call bemsav(xbond, basis)
        do j=1, nbasis-1
            txinput(j) = basis(j)
        enddo
        call evaluate_network(txinput, 1, vpesa)
        call evaluate_network(txinput, 2, vpesb)
        call evaluate_network(txinput, 3, vpesc)
        vpes = (vpesa+vpesb+vpesc)/3.0d0
        if( vpes .lt. -1.5d0 ) vpes=5.0d0
        vpes = min(vpes, 5.0d0)
    end subroutine fh2onn

    !> @brief Evaluate one member of the immutable neural-network ensemble.
    subroutine evaluate_network( x, network_index, potential )
        use chempotpy_fh2o_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_2_biases, hidden_1_count, hidden_2_count, hidden_hidden_weights, &
            & hidden_output_weights, input_hidden_weights, network_count, ninput, nlayer, nodemax, &
            & nodes, noutput, nscale, output_biases, scale_averages, scale_deltas, &
            & alpha, bohr, nbasis, npbasissoc, nrbasissoc, p2soccoef, &
            & pi, r1soccoef, radian, vpesmin
        implicit none
        double precision, intent(in)  :: x(ninput)
        integer, intent(in)           :: network_index
        double precision, intent(out) :: potential
        double precision              :: y(nodemax,nlayer)
        integer                       :: input_index, node_index, previous_index
        y = 0.0d0
        do input_index=1, ninput
            y(input_index,1) = (x(input_index)-scale_averages(input_index, network_index))/ scale_deltas(input_index, &
                & network_index)
        enddo
        do node_index=1, hidden_1_count
            y(node_index,2) = hidden_1_biases(node_index, network_index)
            do previous_index=1, ninput
                y(node_index,2) = y(node_index,2)+y(previous_index,1)* input_hidden_weights(previous_index, &
                    & node_index, network_index)
            enddo
            select case( activation_function )
            case( 1 )
                y(node_index,2) = tanh(y(node_index,2))
            case( 2 )
                y(node_index,2) = 1.0d0/(1.0d0+exp(-y(node_index,2)))
            case( 3 )
            end select
        enddo
        do node_index=1, hidden_2_count
            y(node_index,3) = hidden_2_biases(node_index, network_index)
            do previous_index=1, hidden_1_count
                y(node_index,3) = y(node_index,3)+y(previous_index,2)* hidden_hidden_weights(previous_index, &
                    & node_index, network_index)
            enddo
            select case( activation_function )
            case( 1 )
                y(node_index,3) = tanh(y(node_index,3))
            case( 2 )
                y(node_index,3) = 1.0d0/(1.0d0+exp(-y(node_index,3)))
            case( 3 )
            end select
        enddo
        do node_index=1, noutput
            y(node_index,nlayer) = output_biases(node_index, network_index)
            do previous_index=1, hidden_2_count
                y(node_index,nlayer) = y(node_index,nlayer)+y(previous_index,3)* hidden_output_weights(previous_index, &
                    & node_index, network_index)
            enddo
        enddo
        potential = y(nodes(nlayer),nlayer)*scale_deltas(nscale, network_index)+ scale_averages(nscale, network_index)
    end subroutine evaluate_network

    !> @brief Evaluate emsav.
    function emsav( x, c ) result(v)
        implicit none
        double precision, dimension(1:6)::x
        double precision, dimension(0:17)::c
        double precision :: v
        double precision, dimension(0:17)::p
        call bemsav(x, p)
        v = dot_product(p, c)
    end function emsav

    !> @brief Evaluate bemsav.
    subroutine bemsav( x, p )
        implicit none
        double precision, dimension(1:6), intent(in)::x
        double precision, dimension(0:17), intent(out)::p
        double precision, dimension(0:10)::m
        call evmono(x, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        implicit none
        double precision, dimension(1:6), intent(in)::x
        double precision, dimension(0:10), intent(out)::m
        m(0)  = 1.d0
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
        implicit none
        double precision, dimension(0:10), intent(in)::m
        double precision, dimension(0:17), intent(out)::p
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

    !> @brief Evaluate xyzTOjcb3.
    subroutine xyztojcb3( ct, rt )
        implicit none
        !   D  F
        !           /|\         /
        !     A  /\  |         /
        !     H    \ |        /
        !           \|       /
        !  ----------\------------------/ C  H
        !             \
        !              \
        !               \
        !             O  B
        double precision :: ct(3,4), xab(3), xabc(3), rt(6) !ct in HOH + F= HO + HF order
        double precision :: xvec(3,6), crop1(3), crop2(3)
        double precision :: mass(4)
        mass(1)   = 1.00783d0
        mass(2)   = 15.99491d0
        mass(3)   = 1.00783d0
        mass(4)   = 18.99840d0
        xab(:)    = (mass(1)*ct(:,1)+mass(2)*ct(:,2))/(mass(1)+mass(2))
        xabc(:)   = (mass(1)*ct(:,1)+mass(2)*ct(:,2)+mass(3)*ct(:,3))/ (mass(1)+mass(2)+mass(3))
        xvec(:,1) = ct(:,1)-ct(:,2)
        xvec(:,2) = ct(:,3)-xab(:)
        xvec(:,3) = ct(:,4)-xabc(:)
        rt(1)     = dsqrt(dot_product(xvec(:,1), xvec(:,1)))
        rt(2)     = dsqrt(dot_product(xvec(:,2), xvec(:,2)))
        rt(3)     = dsqrt(dot_product(xvec(:,3), xvec(:,3)))
        rt(4)     = dot_product(xvec(:,1), xvec(:,2))/rt(1)/rt(2)
        rt(5)     = dot_product(xvec(:,2), xvec(:,3))/rt(2)/rt(3)
        rt(4)     = dacos(max(-1.0d0, min(1.0d0, rt(4))))
        rt(5)     = dacos(max(-1.0d0, min(1.0d0, rt(5))))
        call crossp(xvec(:,1), xvec(:,2), crop1)
        call crossp(xvec(:,2), xvec(:,3), crop2)
        rt(6) = dot_product(crop1, crop2)/dsqrt(dot_product(crop1, crop1))/ dsqrt(dot_product(crop2, crop2))
        rt(6) = dacos(max(-1.0d0, min(1.0d0, rt(6))))
        rt(6) = dacos(-1.d0)-rt(6)
    end subroutine xyztojcb3

    !> @brief Evaluate jcb3tocart.
    subroutine jcb3tocart( rint, cart )
        implicit none
        double precision :: rint(6), cart(3,4)
        double precision :: r1, r2, r3, th1, th2, phia, sth1, sth2, sth3, cth1, cth2, cth3
        double precision :: wa, wb, wc, wd, abm, cdm, ama, amb, amc, amd
        !-->-->rint(1-6)=r1,r2,r3,th1,th2,phi
        !-->-->r1 is vector A->B, r2 is a vector pointing from the COM of AB to
        !C
        !-->-->r3 is a vector pointing from the COM of ABC to D
        !-->-->th1 is the angle for r1 and r2
        !-->-->th2 is the angle for r2 and r3
        !-->-->phia is the diheral angle for r1*r2 plane and r2*r3 plane.
        !-->-->cart(3,1-4)=cartesian coordinates of four atoms A,B,C,D.
        double precision :: temp1, temp2, temp3
        ama  = 1.00783d0
        amb  = 15.99491d0
        amc  = 1.00783d0
        amd  = 18.99840d0
        r1   = rint(1)
        r2   = rint(2)
        r3   = rint(3)
        th1  = rint(4)
        th2  = rint(5)
        phia = rint(6)
        !-->-->ABC+D Jacobi coordinates
        !-->-->            Schematic diagram
        !-->-->origin is the center of mass of AB.
        !-->-->r2 lies on Z positive axis, AB is in the XZ positive plane
        !-->-->Xa is always positive and Xb is always negative.
        !-->-->Xc is always positive
        abm       = r1/(ama+amb)
        wa        = amb*abm
        wb        = ama*abm
        sth1      = dsin(th1)
        cth1      = dcos(th1)
        sth2      = dsin(th2)
        cth2      = dcos(th2)
        cth3      = dcos(phia)
        sth3      = dsin(phia)
        cart(1,1) = wa*sth1
        cart(2,1) = 0d0
        cart(3,1) = wa*cth1
        cart(1,2) = -wb*sth1
        cart(2,2) = 0d0
        cart(3,2) = -wb*cth1
        cart(1,3) = 0d0
        cart(2,3) = 0d0
        cart(3,3) = r2
        cart(1,4) = r3*sth2*cth3
        cart(2,4) = r3*sth2*sth3
        cart(3,4) = r3*cth2+amc*r2/(ama+amb+amc)
    end subroutine jcb3tocart

    !> @brief Evaluate cart_diJacobi.
    subroutine cart_dijacobi( c1, r )
        ! convert the diatom-diatom Jacobi coordinates to XYZ coordinates
        !  c in angstrom and ABCD, r in angstrom and radians
        implicit none
        double precision :: r(1:6), c(1:3,1:4), mass(1:4), c1(1:3,1:4)
        double precision :: crop1(1:3), crop2(1:3)
        double precision :: vec(1:3,1:9), dotpp(1:9), bohr
        double precision :: dotp, costh1, costh2, amass, bmass, cmass, dmass, costh3
        c        = c1  ! c in HOHF order angstrom
        bohr     = 1.0d0
        amass    = 2.0141d0
        bmass    = 15.99491d0
        cmass    = amass
        dmass    = 18.9984d0
        vec(:,1) = c(:,1)-c(:,2) !  vector B->A OH
        r(1)     = dsqrt( sum(vec(:,1)**2) )/bohr
        vec(:,2) = c(:,2)+vec(:,1)*amass/(amass+bmass) !vector BX (X COM of BA)
        vec(:,3) = c(:,3)-c(:,4) ! vector D->C OC
        r(2)     = dsqrt( sum(vec(:,3)**2) )/bohr
        vec(:,4) = c(:,4)+vec(:,3)*cmass/(cmass+dmass) !vector Y(Y COM of
        vec(:,5) = vec(:,4)-vec(:,2) !vecor R X->Y
        r(3)     = dsqrt( sum(vec(:,5)**2) )/bohr
        !-
        vec(:,6) = c(:,3)-c(:,2) ! vector  O->C
        vec(:,7) = c(:,2)-c(:,3) ! vector  C->O
        vec(:,8) = c(:,4)-c(:,3) ! vector  C->O'
        !     dompp(5)=dsqrt(dot_product(vec(:,5),vec(:,5)))
        !     dotpp(6)=dsqrt(dot_product(vec(:,6),vec(:,6)))
        !     dotpp(7)=dsqrt(dot_product(vec(:,7),vec(:,7)))
        dotpp(1) = dsqrt(dot_product(vec(:,1), vec(:,1)))
        dotpp(3) = dsqrt(dot_product(vec(:,3), vec(:,3)))
        dotpp(5) = dsqrt(dot_product(vec(:,5), vec(:,5)))
        if( dotpp(1)*dotpp(5) .lt. 1.0d-16 )then
            costh1 = 0.0d0
        else
            costh1 = dot_product(vec(:,1), vec(:,5))/(dotpp(1)*dotpp(5))
        endif
        if( dotpp(3)*dotpp(5) .lt. 1.0d-16 )then
            costh2 = 0.0d0
        else
            costh2 = dot_product(vec(:,3), vec(:,5))/(dotpp(3)*dotpp(5))
        endif
        !  note the round off errors
        r(4) = dacos( max(-1.0d0, min(1.0d0, costh1))) ! in radians
        r(5) = dacos( max(-1.0d0, min(1.0d0, costh2))) ! in radians
        if( dabs(dcos(r(4))*dcos(r(5))) .lt. 1.0d-16 )then
            !       return
        else
            ! compute the Jacobi dihedral angle, based on the norm of the plane
            call crossp(vec(:,1), vec(:,5), crop1) ! crossproduct for the norm
            call crossp(vec(:,3), vec(:,5), crop2) !
            if( dsqrt( sum(crop1(:)**2))*dsqrt(sum(crop2(:)**2)) .lt. 1.0d-10 ) &
                &  then
            r(6) = 0.0d0
        else
            costh3 = dot_product(crop1(:), crop2(:))/(dsqrt(sum(crop1(:)**2)) *dsqrt(sum(crop2(:)**2)))
            !    write(99999,'(f20.10)') costh3
            r(6) = dacos( max(-1.0d0, min(1.0d0, costh3)))
        endif
        endif
    end subroutine cart_dijacobi

    !> @brief Evaluate diJacobi_cart.
    subroutine dijacobi_cart( r, c )
        ! convert the diatom-diatom Jacobi coordinates to XYZ coordinates
        ! r in angstrom and radians,  c in bohr and ABCD
        implicit none
        double precision :: r(1:6), c(1:3,1:4), mass(1:4)
        double precision :: cm(3,2), vec(1:3,1:6), pi
        double precision :: xa, xb, yc, yd, cxsin, dxsin, dotp
        double precision :: rx1, rx2, rx3, rx4, rx5, rx6
        double precision :: amass, bmass, cmass, dmass
        pi      = dacos(-1.0d0)
        amass   = 1.00783d0
        bmass   = 15.99491d0
        cmass   = amass
        dmass   = 18.99834d0
        mass(1) = amass
        mass(2) = bmass
        mass(3) = cmass
        mass(4) = dmass
        !     diatomic-diatomic AB+CD Jacobi is
        !     r1, r2, R, theta1, theta2, phi
        !     r1, distance AB
        !     r2, distance CD
        !     R,  distance between the centers of mass AB(X) and CD(Y)
        !     theta1, angle from vector B->A to vector X->Y
        !     theta2, angle from vector D->C to vector X->Y
        !     phi, dihedral angle between the plane AXY and plane XYD
        !   A/- r1                D\
        !     \                     \ r2
        !      \ theta1    R         \
        !      X\-------------------->\ Y--->
        !        \                     \ th2
        !       B \                    -/ C
        !   put ABY on xy plane, X origin point, X->Y (R)  on the x axis
        !   and also put A in the positive y direction
        !             ^y axis
        !        A    |
        !        /\-  |
        !        r1\  |
        !           \ |
        !            \|     R
        !    ---------X--------->Y-->x axis
        !              \
        !               \
        !                B
        c  = 0.0d0
        xa = bmass/(amass+bmass)*r(1)
        xb = amass/(amass+bmass)*r(1)
        !    A
        c(1,1) = xa*dcos(r(4))
        c(2,1) = xa*dsin(r(4))
        c(3,1) = 0.0d0
        !    B
        c(1,2) = xb*dcos(pi-r(4))
        c(2,2) = -xb*dsin(pi-r(4))
        c(3,2) = 0.d0
        !     put ABY on xy plane, X origin point
        !               ^y
        !          A    |
        !          /\-  |
        !          r1\  |               \D
        !             \ |                \ |
        !              \|     R           \|
        !      ---------X----------------->Y-->x
        !                \                 |\
        !                 \                |-\/C
        !                  B
        !      C in +z and -z phase is identical, so put it in +z
        yc    = dmass/(cmass+dmass)*r(2)
        yd    = cmass/(cmass+dmass)*r(2)
        cxsin = yc*dsin(r(5))
        dxsin = -yd*dsin(r(5))
        !     C
        c(1,3) = r(3)+yc*dcos(r(5))
        c(2,3) = cxsin*dcos(r(6))
        c(3,3) = cxsin*dsin(r(6))
        !     D
        c(1,4) = r(3)-yd*dcos(r(5))
        c(2,4) = dxsin*dcos(r(6))
        c(3,4) = dxsin*dsin(r(6))
    end subroutine dijacobi_cart

    !> @brief Evaluate crossp.
    subroutine crossp( v1, v2, v3 )
        implicit none
        double precision :: v1(1:3), v2(1:3), v3(1:3)
        ! cross product of two 3-dimensional vectors
        v3(1) = v1(2)*v2(3)-v1(3)*v2(2)
        v3(2) = v1(3)*v2(1)-v1(1)*v2(3)
        v3(3) = v1(1)*v2(2)-v1(2)*v2(1)
    end subroutine crossp
end module chempotpy_fh2o_pipnn_pes

!> @brief Expose FH2O_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_fh2o_pipnn_pes, only: evaluate
    implicit none
    integer, parameter              :: natoms  = 4
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    double precision, intent(in)    :: x(natoms,3)
    double precision, intent(out)   :: p(nstates)
    double precision, intent(out)   :: g(nstates,natoms,3)
    double precision, intent(out)   :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
