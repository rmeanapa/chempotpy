module chempotpy_h3cl_pipnn_pes
    use chempotpy_h3cl_pipnn_parameters, only: activation_function, alpha, h2_distance_grid, &
        & h2_potential_grid, hcl_distance_grid, hcl_potential_grid, hidden_1_biases, hidden_1_count, &
        & hidden_2_biases, hidden_2_count, hidden_hidden_weights, hidden_output_weights, &
        & input_hidden_weights, long_range_coefficients, natom, nbasis, ndim, ninput, nscale, &
        & output_biases, scale_averages, scale_deltas
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3Cl_PIPNN potential energy surface for H3Cl.
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
        real(wp)              :: v
        real(wp)              :: tx(3,natoms)
        integer               :: iatom, idir, j, istate
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        if( igrad == 0 )then
            call pesh3cl(tx, v)
        else
            write (*, *) 'Only energy and gradient are available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
        g = 0.e0_wp
        d = 0.e0_wp
    end subroutine evaluate

    !> @brief Evaluate pesh3cl for the H3Cl_PIPNN potential.
    subroutine pesh3cl( x, vtot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), parameter :: tocm = 219474.63137e0_wp
        integer, parameter  :: n    = 644
        real(wp)            :: phi, r(6), r1, r2, r3, s, th1, th2, v, v0, v1, vint, vpes, vpesa, vpesb, vpesc, vpip, &
            & vtot, x(3,4)
        call xyz_to_jaco(x, r1, r2, r3, th1, th2, phi, r)
        call clh3pipnn(r, vpes, vpesa, vpesb, vpesc)
        call pot0(r1, v0)
        call pot1(r2, v1)
        call clh3vlr(r1, r2, r3, th1, th2, phi, vint, n)
        vpip = (vpes*0.03674931e0_wp+-461.54119162e0_wp)-v0-v1
        vpip = vpip*tocm
        s    = (1e0_wp-tanh(3e0_wp*(r3-5.5e0_wp)))/2e0_wp
        v    = s*vpip+(1e0_wp-s)*vint
        vtot = (v/tocm+v0+v1--461.54020987e0_wp)*27.21138e0_wp
    end subroutine pesh3cl

    !> @brief Interpolate the immutable H2 diatomic correction table.
    subroutine pot0( r0, vv0 )
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter    :: m0 = 1425
        real(wp), intent(in)  :: r0
        real(wp), intent(out) :: vv0
        real(wp)              :: dy1, dyn, second_derivatives(m0)
        dy1 = 1.0e30_wp; dyn = 1.0e30_wp
        call spline(h2_distance_grid, h2_potential_grid, m0, dy1, dyn, second_derivatives)
        call splint(h2_distance_grid, h2_potential_grid, second_derivatives, m0, r0, vv0)
    end subroutine pot0

    !> @brief Interpolate the immutable HCl diatomic correction table.
    subroutine pot1( r1, vv1 )
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter    :: m1 = 3200
        real(wp), intent(in)  :: r1
        real(wp), intent(out) :: vv1
        real(wp)              :: dy1, dyn, second_derivatives(m1)
        dy1 = 1.0e30_wp; dyn = 1.0e30_wp
        call spline(hcl_distance_grid, hcl_potential_grid, m1, dy1, dyn, second_derivatives)
        call splint(hcl_distance_grid, hcl_potential_grid, second_derivatives, m1, r1, vv1)
    end subroutine pot1

    !# SPLINE ROUTINES
    !#            Numerical recipes in fortran
    !#            Cambrige University Press
    !#            York, 2nd edition, 1992.
    !> @brief Evaluate splint for the H3Cl_PIPNN potential.
    subroutine splint( xa, ya, y2a, n, x, y )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a, b, h, x, xa(n), y, y2a(n), ya(n)
        integer  :: k, khi, klo, n
        klo = 1
        khi = n
        do while( khi-klo .gt. 1 )
            k = (khi+klo)/2
            if( xa(k) .gt. x )then
                khi = k
            else
                klo = k
            endif
        enddo
        h = xa(khi)-xa(klo)
        if( h .eq. 0.0e0_wp ) write(6, *) 'bad xa input in splint'
        a = (xa(khi)-x)/h
        b = (x-xa(klo))/h
        y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h** 2)/6.0e0_wp
    end subroutine splint

    !> @brief Evaluate spline for the H3Cl_PIPNN potential.
    subroutine spline( x, y, n, yp1, ypn, y2 )
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter :: nmax = 4200
        real(wp)           :: p, qn, sig, u(nmax), un, x(n), y(n), y2(n), yp1, ypn
        integer            :: i, k, n
        if( yp1 .gt. .99e30_wp )then
            y2(1) = 0.0e0_wp
            u(1)  = 0.0e0_wp
        else
            y2(1) = -0.5e0_wp
            u(1)  = (3.0e0_wp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif
        do i=2, n-1
            sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
            p     = sig*y2(i-1)+2.0e0_wp
            y2(i) = (sig-1.0e0_wp)/p
            u(i)  = (6.0e0_wp*((y(i+1)-y(i))/(x(i+ 1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* u(i-1))/p
        enddo
        if( ypn .gt. .99e30_wp )then
            qn = 0.0e0_wp
            un = 0.0e0_wp
        else
            qn = 0.5e0_wp
            un = (3.0e0_wp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
        y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0e0_wp)
        do k=n-1, 1, -1
            y2(k) = y2(k)*y2(k+1)+u(k)
        enddo
    end subroutine spline

    !> @brief Evaluate the immutable H3Cl long-range expansion.
    subroutine clh3vlr( r1, r2, r3, theta1, theta2, phi, vint, n )
        use chempotpy_kinds, only: wp
        implicit none
        integer, intent(in)   :: n
        real(wp), intent(in)  :: phi, r1, r2, r3, theta1, theta2
        real(wp), intent(out) :: vint
        real(wp)              :: afunc(n)
        call funcs(r1, r2, r3, theta1, theta2, phi, afunc, n)
        vint = dot_product(long_range_coefficients(1:n), afunc)
    end subroutine clh3vlr

    !> @brief Evaluate FUNCS for the H3Cl_PIPNN potential.
    subroutine funcs( ra, rb, rc, theta1, theta2, phi, afunc, ma )
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter :: m = 3, m1=2*m+1, m2=m1**2! m the order of r1 and r2
        real(wp)           :: z1(m2), z2(m2), z3(m2), z4(m2), z5(m1), z6(m2), z7(m2), &
            & afunc(ma), a(m2), z8(m2), z9(m2), z10(m2), z11(m2), z12(m2), z13(m2) &
            &, z14(m2), a111(m2)
        integer  :: ma, i, l, n, j, j1
        real(wp) :: ra, rb, rc, theta1, theta2, phi, r1, r2, r3, th1, th2, a1, a2, a3 &
            &, a11, a22
        r2 = ra;r1=rb;r3=rc;th1=theta2;th2=theta1
        n  = 1
        j  = 0
        do i=-m, m
            a1 = r1**i
            do l=-m, m
                a2 = (r2**l)*a1
                if( n .gt. 0 )then
                    j    = j+1
                    a(j) = a2
                endif
            enddo;enddo
            j1 = 0
            do i=-m, m
                a11 = r1**i
                do l=-m, m
                    a22 = 3e0_wp*(r2**l)*a11/(2e0_wp*(r2**l+a11))
                    if( n .gt. 0 )then
                        j1       = j1+1
                        a111(j1) = a22
                    endif
                enddo;enddo
                do i=1, m2
                    z1(i)=a(i)*(1.5e0_wp*(cos(th1)*(3e0_wp*(cos(th2)**2)-1e0_wp) +2e0_wp*sin(th1)*sin(th2)*cos(th2)* &
                        & cos(phi)))/ (r3**4)
                    z2(i)=a(i)*(0.75e0_wp*(1e0_wp-5e0_wp*(cos(th1)**2)- 5e0_wp*(cos(th2)**2)+17e0_wp*(cos(th1)**2)* &
                        & (cos(th2)**2)+ 2e0_wp* (sin(th1)**2)*(sin(th2)**2)*(cos(phi)**2)+ 16e0_wp*sin(th1)*cos(th1)* &
                        & sin(th2)*cos(th2)* cos(phi)))/(r3**5)
                    z3(i) = a(i)*(-0.5e0_wp*(3e0_wp*(cos(th1)**2)+1e0_wp))/(r3**6)
                    z4(i)=a(i)*((12e0_wp*(cos(th1)**2)*(cos(th2)**2)+3e0_wp* (sin(th1)**2)*(sin(th2)**2)*(cos(phi)** &
                        & 2)- 3e0_wp* (cos(th1)**2)-1e0_wp +12e0_wp*sin(th1)*cos(th1)*sin(th2)*cos(th2)*cos(phi)) /- &
                        & 6e0_wp)/ (r3**6)
                    z6(i)=a(i)*(-1.5e0_wp*(6e0_wp*(cos(th1)**2)*(cos(th2)**3) -2e0_wp*(cos(th1)**2)*cos(th2)+7e0_wp* &
                        & sin(th1)* cos(th1)*sin(th2)* (cos(th2)**2)*cos(phi)-sin(th1)*cos(th1)*sin(th2)*cos(phi) + &
                        & 2e0_wp* (sin(th1)**2)*(sin(th2)**2)*cos(th2)*(cos(phi)**2))) /(r3**7)
                    z7(i)=a(i)*(2e0_wp*(6e0_wp*(cos(th1)**2)*(cos(th2)**3) -5e0_wp*(cos(th1)**2)*cos(th2)+7e0_wp* &
                        & sin(th1)* cos(th1)*sin(th2)* (cos(th2)**2)*cos(phi)-2e0_wp*sin(th1)*cos(th1)*sin(th2) * &
                        & cos(phi)+2e0_wp* (sin(th1)**2)*(sin(th2)**2)*cos(th2) *(cos(phi)**2)-cos(th2)))/(r3**7)
                    z8(i)             = -a111(i)/(r3**6)
                    z9(i)             = (-a111(i)*(1.5e0_wp*(cos(th2)**2)-0.5e0_wp)/3e0_wp)/(r3**6)
                    z10(i)            = (-a111(i)*(1.5e0_wp*(cos(th1)**2)-0.5e0_wp)/3e0_wp)/(r3**6)
                    z11(i)            = -a111(i)*2e0_wp*(cos(th2)**3)/(r3**7)
                    z12(i)            = -a111(i)*2e0_wp*(cos(th1)**3)/(r3**7)
                    z13(i)            = -a111(i)*4e0_wp*(3e0_wp*cos(th2)-2e0_wp*(cos(th2)**3))/ (3e0_wp*(r3**7))
                    z14(i)            = -a111(i)*4e0_wp*(3e0_wp*cos(th1)-2e0_wp*(cos(th1)**3))/ (3e0_wp*(r3**7))
                    afunc(i)          = z1(i);afunc(m2+i)=z2(i);afunc(2*m2+i)=z3(i)
                    afunc(3*m2+i)     = z4(i);afunc(4*m2+m1+i)=z6(i)
                    afunc(5*m2+m1+i)  = z7(i);afunc(6*m2+m1+i)=z8(i)
                    afunc(7*m2+m1+i)  = z9(i);afunc(8*m2+m1+i)=z10(i)
                    afunc(9*m2+m1+i)  = z11(i);afunc(10*m2+m1+i)=z12(i)
                    afunc(11*m2+m1+i) = z13(i);afunc(12*m2+m1+i)=z14(i)
                enddo
                j = 4*m2
                do i=-m, m
                    a3 = (r1**i)*(-6e0_wp*(cos(th1)**3))/(r3**7)
                    if( n .gt. 0 )then
                        j        = j+1
                        afunc(j) = a3
                    endif
                enddo
    end subroutine funcs

    !> @brief Evaluate ClH3pipNN for the H3Cl_PIPNN potential.
    subroutine clh3pipnn( rb, vpes, vpesa, vpesb, vpesc )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: i, j, k, lpes
        real(wp) :: rb(ndim), xbond(ndim), basis(0:nbasis), tmp1, txinput(nbasis)
        real(wp) :: vpes, vpesa, vpesb, vpesc, vpesh, vnn
        basis    = 0.e0_wp
        xbond(:) = exp(-rb(:)/alpha)
        call bemsav134(xbond, basis)
        do j=1, nbasis
            txinput(j) = basis(j)
        enddo
        call getpota(txinput, vpesa)
        call getpotb(txinput, vpesb)
        call getpotc(txinput, vpesc)
        vpes = (vpesa+vpesb+vpesc)/3e0_wp
    end subroutine clh3pipnn

    !> @brief Evaluate tranfun for the H3Cl_PIPNN potential.
    function tranfun( x, ifunc )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: ifunc
        real(wp) :: tranfun, x
        !    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
        !    ifunc=2, transfer function is log sigmoid function, 'logsig'
        !    ifunc=3, transfer function is pure linear function, 'purelin'. It
        !is imposed to the output layer by default
        if( ifunc .eq. 1 )then
            tranfun = tanh(x)
        elseif( ifunc .eq. 2 )then
            tranfun = 1e0_wp/(1e0_wp+exp(-x))
        elseif( ifunc .eq. 3 )then
            tranfun = x
        endif
    end function tranfun

    !-->!-->!-->!-->!-->!-->!-->!-->!-->!-->subroutine to accelerate~~~~~~
    !> @brief Evaluate emsav134 for the H3Cl_PIPNN potential.
    function emsav134( x, c ) result(v )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(1:6)
        real(wp) :: c(0:50)
        real(wp) :: v
        real(wp) :: p(0:50)
        call bemsav134(x, p)
        v = dot_product(p, c)
    end function emsav134

    !> @brief Evaluate bemsav134 for the H3Cl_PIPNN potential.
    subroutine bemsav134( x, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:6)
        real(wp), intent(out) :: p(0:50)
        real(wp)              :: m(0:32)
        call evmono134(x, m)
        call evpoly134(m, p)
    end subroutine bemsav134

    !> @brief Evaluate evmono134 for the H3Cl_PIPNN potential.
    subroutine evmono134( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(1:6)
        real(wp), intent(out) :: m(0:32)
        m(0)  = 1.e0_wp
        m(1)  = x(6)
        m(2)  = x(5)
        m(3)  = x(4)
        m(4)  = x(3)
        m(5)  = x(2)
        m(6)  = x(1)
        m(7)  = m(1)*m(2)
        m(8)  = m(1)*m(3)
        m(9)  = m(2)*m(3)
        m(10) = m(3)*m(4)
        m(11) = m(2)*m(5)
        m(12) = m(1)*m(6)
        m(13) = m(4)*m(5)
        m(14) = m(4)*m(6)
        m(15) = m(5)*m(6)
        m(16) = m(1)*m(9)
        m(17) = m(1)*m(10)
        m(18) = m(2)*m(10)
        m(19) = m(1)*m(11)
        m(20) = m(3)*m(11)
        m(21) = m(2)*m(12)
        m(22) = m(3)*m(12)
        m(23) = m(2)*m(13)
        m(24) = m(3)*m(13)
        m(25) = m(1)*m(14)
        m(26) = m(3)*m(14)
        m(27) = m(1)*m(15)
        m(28) = m(2)*m(15)
        m(29) = m(4)*m(15)
        m(30) = m(2)*m(24)
        m(31) = m(1)*m(26)
        m(32) = m(1)*m(28)
    end subroutine evmono134

    !> @brief Evaluate evpoly134 for the H3Cl_PIPNN potential.
    subroutine evpoly134( m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: m(0:32)
        real(wp), intent(out) :: p(0:50)
        p(0)  = m(0)
        p(1)  = m(1)+m(2)+m(3)
        p(2)  = m(4)+m(5)+m(6)
        p(3)  = m(7)+m(8)+m(9)
        p(4)  = m(10)+m(11)+m(12)
        p(5)  = p(1)*p(2)-p(4)
        p(6)  = m(13)+m(14)+m(15)
        p(7)  = p(1)*p(1)-p(3)-p(3)
        p(8)  = p(2)*p(2)-p(6)-p(6)
        p(9)  = m(16)
        p(10) = m(17)+m(18)+m(19)+m(20)+m(21)+m(22)
        p(11) = p(2)*p(3)-p(10)
        p(12) = m(23)+m(24)+m(25)+m(26)+m(27)+m(28)
        p(13) = p(1)*p(6)-p(12)
        p(14) = m(29)
        p(15) = p(1)*p(3)-p(9)-p(9)-p(9)
        p(16) = p(1)*p(4)-p(10)
        p(17) = p(2)*p(7)-p(16)
        p(18) = p(2)*p(4)-p(12)
        p(19) = p(1)*p(8)-p(18)
        p(20) = p(2)*p(6)-p(14)-p(14)-p(14)
        p(21) = p(1)*p(7)-p(15)
        p(22) = p(2)*p(8)-p(20)
        p(23) = p(9)*p(2)
        p(24) = m(30)+m(31)+m(32)
        p(25) = p(3)*p(6)-p(24)
        p(26) = p(14)*p(1)
        p(27) = p(9)*p(1)
        p(28) = p(3)*p(4)-p(23)
        p(29) = p(1)*p(10)-p(23)-p(28)-p(23)
        p(30) = p(1)*p(11)-p(23)
        p(31) = p(1)*p(12)-p(25)-p(24)-p(24)
        p(32) = p(1)*p(13)-p(25)
        p(33) = p(4)*p(5)-p(25)-p(31)
        p(34) = p(2)*p(11)-p(25)
        p(35) = p(4)*p(6)-p(26)
        p(36) = p(2)*p(12)-p(26)-p(35)-p(26)
        p(37) = p(2)*p(13)-p(26)
        p(38) = p(14)*p(2)
        p(39) = p(3)*p(3)-p(27)-p(27)
        p(40) = p(3)*p(7)-p(27)
        p(41) = p(1)*p(16)-p(28)
        p(42) = p(2)*p(21)-p(41)
        p(43) = p(1)*p(18)-p(33)
        p(44) = p(7)*p(8)-p(43)
        p(45) = p(6)*p(6)-p(38)-p(38)
        p(46) = p(2)*p(18)-p(35)
        p(47) = p(1)*p(22)-p(46)
        p(48) = p(6)*p(8)-p(38)
        p(49) = p(1)*p(21)-p(40)
        p(50) = p(2)*p(22)-p(48)
    end subroutine evpoly134

    !> @brief Evaluate one immutable neural network from the released three-model ensemble.
    subroutine evaluate_network( x, network_index, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        integer, intent(in)   :: network_index
        real(wp), intent(out) :: vpot
        integer               :: i, j, k
        real(wp)              :: input_values(ninput), hidden_1_values(hidden_1_count)
        real(wp)              :: hidden_2_values(hidden_2_count), output_value
        do i=1, ninput
            input_values(i) = (x(i)-scale_averages(i, network_index))/scale_deltas(i, network_index)
        enddo
        do j=1, hidden_1_count
            hidden_1_values(j) = hidden_1_biases(j, network_index)
            do i=1, ninput
                hidden_1_values(j) = hidden_1_values(j) + input_values(i)*input_hidden_weights(i, j, network_index)
            enddo
            hidden_1_values(j) = tranfun(hidden_1_values(j), activation_function)
        enddo
        do k=1, hidden_2_count
            hidden_2_values(k) = hidden_2_biases(k, network_index)
            do j=1, hidden_1_count
                hidden_2_values(k) = hidden_2_values(k) + hidden_1_values(j)*hidden_hidden_weights(j, k, network_index)
            enddo
            hidden_2_values(k) = tranfun(hidden_2_values(k), activation_function)
        enddo
        output_value = output_biases(1, network_index)
        do k=1, hidden_2_count
            output_value = output_value + hidden_2_values(k)*hidden_output_weights(k, 1, network_index)
        enddo
        vpot = output_value*scale_deltas(nscale, network_index)+scale_averages(nscale, network_index)
    end subroutine evaluate_network

    !> @brief Evaluate released neural network 1 from the H3Cl ensemble.
    subroutine getpota( x, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot
        call evaluate_network(x, 1, vpot)
    end subroutine getpota

    !> @brief Evaluate released neural network 2 from the H3Cl ensemble.
    subroutine getpotb( x, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot
        call evaluate_network(x, 2, vpot)
    end subroutine getpotb

    !> @brief Evaluate released neural network 3 from the H3Cl ensemble.
    subroutine getpotc( x, vpot )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot
        call evaluate_network(x, 3, vpot)
    end subroutine getpotc

    !> @brief Evaluate distance for the H3Cl_PIPNN potential.
    real(wp) function distance( a, b )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a(3), b(3)
        integer  :: i
        distance = 0.e0_wp
        do i=1, 3
            distance = distance+(a(i)-b(i))**2
        enddo
        distance = sqrt(distance)
    end function distance

    !> @brief Evaluate xyz_to_jaco for the H3Cl_PIPNN potential.
    subroutine xyz_to_jaco( x, r1, r2, r3, theta1, theta2, phi, r )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: ra1, ra4, rb1, rb4, str11, str12
        real(wp)            :: x(3,4), a(3), b(3), r(6)
        real(wp)            :: mass(4)
        real(wp)            :: r1, r2, r3, theta1, theta2, phi
        real(wp), parameter :: pi = 3.14159265e0_wp
        integer             :: k, m, i, j, n, ndata
        character           :: long
        mass(1) = 1.00782503e0_wp
        mass(2) = mass(1)
        mass(3) = mass(1)
        mass(4) = 35.4532e0_wp
        r1      = distance(x(:,1), x(:,2))
        r2      = distance(x(:,3), x(:,4))
        a(:)    = x(:,1)+(mass(2)/(mass(1)+mass(2)))*(x(:,2)-x(:,1))
        b(:)    = x(:,4)+(mass(3)/(mass(4)+mass(3)))*(x(:,3)-x(:,4))
        r3      = distance(a(:), b(:))
        ra4     = distance(a(:), x(:,4))
        rb4     = distance(b(:), x(:,4))
        str11   = (r3**2+rb4**2-ra4**2)/(2e0_wp*r3*rb4)
        str11   = min(1e0_wp, str11);str11=max(-1e0_wp, str11)
        theta2  = acos(str11)
        ra1     = distance(a(:), x(:,1))
        rb1     = distance(b(:), x(:,1))
        str12   = (r3**2+ra1**2-rb1**2)/(2e0_wp*r3*ra1)
        str12   = min(1e0_wp, str12);str12=max(-1e0_wp, str12)
        theta1  = acos(str12)
        call dihedral(x(:,2), a, b, x(:,4), phi)
        phi = phi*pi/180e0_wp
        if( isnan(phi) )then
            phi = 0e0_wp
        endif
        r(1) = distance(x(:,4), x(:,3))
        r(2) = distance(x(:,4), x(:,2))
        r(3) = distance(x(:,4), x(:,1))
        r(4) = distance(x(:,3), x(:,2))
        r(5) = distance(x(:,3), x(:,1))
        r(6) = distance(x(:,2), x(:,1))
    end subroutine xyz_to_jaco

    !> @brief Evaluate dihedral for the H3Cl_PIPNN potential.
    subroutine dihedral( a, b, c, d, th )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a(3), b(3), c(3), d(3), th, x, y
        real(wp) :: v_ab(3), v_bc(3), v_cd(3), n1(3), n2(3), s1, s2, s3, m1(3)
        integer  :: i
        do i=1, 3
            v_ab(i) = b(i)-a(i)
            v_bc(i) = c(i)-b(i)
            v_cd(i) = d(i)-c(i)
        enddo
        n1(1) = v_ab(2)*v_bc(3)-v_ab(3)*v_bc(2)
        n1(2) = -v_ab(1)*v_bc(3)+v_ab(3)*v_bc(1)
        n1(3) = v_ab(1)*v_bc(2)-v_ab(2)*v_bc(1)
        n2(1) = v_bc(2)*v_cd(3)-v_bc(3)*v_cd(2)
        n2(2) = -v_bc(1)*v_cd(3)+v_bc(3)*v_cd(1)
        n2(3) = v_bc(1)*v_cd(2)-v_bc(2)*v_cd(1)
        m1(1) = n1(2)*v_bc(3)-n1(3)*v_bc(2)
        m1(2) = -n1(1)*v_bc(3)+n1(3)*v_bc(1)
        m1(3) = n1(1)*v_bc(2)-n1(2)*v_bc(1)
        s1    = 0e0_wp;s2=0e0_wp;s3=0e0_wp
        do i=1, 3
            s1 = s1+n1(i)**2
            s2 = s2+n2(i)**2
            s3 = s3+m1(i)**2
        enddo
        n1 = n1/sqrt(s1);n2=n2/sqrt(s2);m1=m1/sqrt(s3)
        x  = 0e0_wp;y=0e0_wp
        do i=1, 3
            x = x+n1(i)*n2(i)
            y = y+m1(i)*n2(i)
        enddo
        th = -atan2(y, x)
    end subroutine dihedral
end module chempotpy_h3cl_pipnn_pes

!> @brief Expose H3Cl_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_h3cl_pipnn_pes, only: evaluate
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
