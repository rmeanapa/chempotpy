!> @brief Implement ClH2O_PIPNN, a potential energy surface for ClH2O.
module chempotpy_clh2o_pipnn_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the ClH2O_PIPNN potential energy surface for ClH2O.
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
        integer               :: iatom, idir, istate
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        if( igrad == 0 )then
            call hohcl_pes_interface(tx, v)
        else
            write(*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
        g = 0.e0_wp
        d = 0.e0_wp
    end subroutine evaluate

    !> @brief Evaluate clh2oNN.
    subroutine clh2onn( ct, vpes, vpesa, vpesb, vpesc ) !ct in HHClO order and angstrom
        use chempotpy_clh2o_pipnn_parameters, only: alpha, nbasis
        use chempotpy_geometry,               only: all_pair_distances
        use chempotpy_kinds,                  only: wp
        implicit none
        integer, parameter :: ndim = 6
        integer            :: j
        real(wp)           :: rb(ndim), xbond(ndim), basis(0:nbasis-1), tmp1
        real(wp)           :: txinput(1:nbasis-1)
        real(wp)           :: ct(3,4), vpes
        real(wp)           :: vpesa, vpesb, vpesc, r12, r13, r14, r23, r24, r34
        basis = 0.e0_wp
        call all_pair_distances(transpose(ct), rb)
        r12      = rb(1);r13=rb(2);r14=rb(3);r23=rb(4);r24=rb(5);r34=rb(6)
        xbond(:) = exp(-rb(:)/alpha)
        call bemsav(xbond, basis)
        do j=1, nbasis-1
            txinput(j) = basis(j)
        enddo
        call getpot(txinput, 1, vpesa)
        call getpot(txinput, 2, vpesb)
        call getpot(txinput, 3, vpesc)
        vpes = (vpesa+vpesb+vpesc)/3.0e0_wp
        if( vpes .lt. -1.5e0_wp ) vpes=5.0e0_wp
        vpes = min(vpes, 5.0e0_wp)
    end subroutine clh2onn

    !> @brief Evaluate one member of the three-network ClH2O ensemble.
    subroutine getpot( x, network, vpot )
        use chempotpy_clh2o_pipnn_parameters, only: activation_function, hidden_1_biases, &
            & hidden_1_count, hidden_2_biases, hidden_2_count, hidden_hidden_weights, &
            & hidden_output_weights, input_hidden_weights, network_count, ninput, nscale, &
            & output_biases, scale_averages, scale_deltas
        use chempotpy_kinds, only: wp
        implicit none
        integer, intent(in)   :: network
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot
        integer               :: i, inode1, inode2
        real(wp)              :: hidden_1(hidden_1_count), hidden_2(hidden_2_count), output
        real(wp)              :: scaled_input(ninput)
        if( network < 1 .or. network > network_count )then
            error stop "ClH2O_PIPNN: neural-network index is out of range"
        endif
        do i=1, ninput
            scaled_input(i) = (x(i)-scale_averages(i, network))/scale_deltas(i, network)
        enddo
        do inode1=1, hidden_1_count
            hidden_1(inode1) = hidden_1_biases(inode1, network)
            do inode2=1, ninput
                hidden_1(inode1) = hidden_1(inode1) +scaled_input(inode2)*input_hidden_weights(inode2, inode1, network)
            enddo
            hidden_1(inode1) = tranfun(hidden_1(inode1), activation_function)
        enddo
        do inode1=1, hidden_2_count
            hidden_2(inode1) = hidden_2_biases(inode1, network)
            do inode2=1, hidden_1_count
                hidden_2(inode1) = hidden_2(inode1) +hidden_1(inode2)*hidden_hidden_weights(inode2, inode1, network)
            enddo
            hidden_2(inode1) = tranfun(hidden_2(inode1), activation_function)
        enddo
        output = output_biases(1, network)
        do inode2=1, hidden_2_count
            output = output+hidden_2(inode2)*hidden_output_weights(inode2, 1, network)
        enddo
        vpot = output*scale_deltas(nscale, network)+scale_averages(nscale, network)
    end subroutine getpot

    !> @brief Evaluate tranfun.
    function tranfun( x, ifunc ) result(value)
        use chempotpy_kinds, only: wp
        implicit none
        integer, intent(in)  :: ifunc
        real(wp), intent(in) :: x
        real(wp)             :: value
        if( ifunc .eq. 1 )then
            value = tanh(x)
        elseif( ifunc .eq. 2 )then
            value = 1e0_wp/(1e0_wp+exp(-x))
        elseif( ifunc .eq. 3 )then
            value = x
        endif
    end function tranfun

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
        real(wp), dimension(0:10), intent(in)::m
        real(wp), dimension(0:17), intent(out)::p
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

    !> @brief Evaluate HOHCl_pes_interface.
    subroutine hohcl_pes_interface( xc, vpes )
        use chempotpy_clh2o_pipnn_parameters, only: vpescut
        use chempotpy_kinds,                  only: wp
        implicit none
        integer  :: lxc, i, j
        real(wp) :: xc(3,4), rbond(6), xct(3,4), xvec(3,6), xct2(3,4)
        real(wp) :: dpthcl, dptho, r13, r23, rhpo, rhcl, vfit, vdp
        real(wp) :: vt1, vcorr, vpes, tmp1, ybond(6), va, vb, vc
        real(wp) :: xcom1(3), xcom2(3), rcom, mass(4), skey
        mass(1)   = 1.00794e0_wp
        mass(2)   = mass(1)
        mass(3)   = 35.453e0_wp   !HHClO standard weight
        mass(4)   = 15.9994e0_wp
        vt1       = 0.0e0_wp
        vcorr     = 0.0e0_wp
        xct       = xc
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
        call before_vpes(xct2, lxc)  !xct in order for H1-Cl  and H2-O, if possible
        if( lxc .eq. -1 )then
            vpes = vpescut
            return
        endif
        lxc = 0
        if( lxc .eq. 2 )then
            call ho_hcl_correction(xct2, vpes)
            vpes = min(vpes, vpescut)
            return
        endif
        if( lxc .eq. 0 )then
            call clh2onn(xct2, vpes, va, vb, vc)
            vpes = min(vpes, vpescut)
            return
        endif
        if( lxc .eq. 1 )then
            ! check the linked atom order: make sure that
            ! 1st H linked or near Cl
            xcom1(:) = xct2(:,1)+(mass(3)/(mass(1)+mass(3)))* (xct2(:,3)-xct2(:,1))
            xcom2(:) = xct2(:,2)+(mass(4)/(mass(2)+mass(4)))* (xct2(:,4)-xct2(:,2))
            rcom     = sqrt(dot_product(xcom1-xcom2, xcom1-xcom2))
            skey     = (1.0e0_wp-tanh(0.8e0_wp*(rcom-10.0e0_wp)))/2.0e0_wp
            call clh2onn(xct2, vfit, va, vb, vc)
            call ho_hcl_correction(xct2, vdp)
            vpes = skey*vfit+(1.0e0_wp-skey)*vdp
            vpes = min(vpes, vpescut)
            return
        endif
    end subroutine hohcl_pes_interface

    !make sure that 1st H linked or near to Cl
    !> @brief Evaluate HO_HCl_correction.
    subroutine ho_hcl_correction( xc, vdp )
        use chempotpy_clh2o_pipnn_parameters, only: vpesmin
        use chempotpy_kinds,                  only: wp
        implicit none
        real(wp) :: xc(3,4), xct2(3,4), vdp, vcorr, rbond(6), vt1
        real(wp) :: xvec(3,6), rhpo, rhcl, dpthcl, dptho
        real(wp) :: angb, anga, angab, vechcl(3), vecho(3), vechh(3)
        real(wp) :: pointa(3), pointb(3), vecacl(3), vecbo(3), vecab(3)
        real(wp) :: racl, rbo, rab
        xct2      = xc
        vt1       = 0.0e0_wp
        vcorr     = 0.e0_wp
        xvec(:,1) = xct2(:,2)-xct2(:,1) !  H2->H1
        xvec(:,2) = xct2(:,3)-xct2(:,1) !  Cl->H1
        xvec(:,3) = xct2(:,4)-xct2(:,1) !  O->H1
        xvec(:,4) = xct2(:,3)-xct2(:,2) !  H2->Cl
        xvec(:,5) = xct2(:,4)-xct2(:,2) !  H2->O
        xvec(:,6) = xct2(:,4)-xct2(:,3) !  Cl-O
        rbond(1)  = sqrt(dot_product(xvec(:,1), xvec(:,1)))
        rbond(2)  = sqrt(dot_product(xvec(:,2), xvec(:,2)))
        rbond(3)  = sqrt(dot_product(xvec(:,3), xvec(:,3)))
        rbond(4)  = sqrt(dot_product(xvec(:,4), xvec(:,4)))
        rbond(5)  = sqrt(dot_product(xvec(:,5), xvec(:,5)))
        rbond(6)  = sqrt(dot_product(xvec(:,6), xvec(:,6)))
        rhcl      = rbond(2)
        rhpo      = rbond(5)
        !     Cl /\                       /\ O
        !          \                     /
        !           \ thetaA            /thetaB
        !         A  \-----------------/-----/
        !             \               / B
        !              \             /
        !               \           / H
        !             H  \
        call ho_hcl_pes_spline(rhcl, rhpo, vt1)
        rhcl = rhcl/0.5291772e0_wp
        rhpo = rhpo/0.52917720e0_wp
        call hcl_dipole(rhcl, dpthcl)
        call ho_dipole(rhpo, dptho)
        !      dipole should be for the middle point of HCl and HO
        vechcl(:) = xvec(:,2)
        vecho(:)  = xvec(:,5)
        pointa(:) = xct2(:,1)+vechcl(:)/2.0e0_wp
        pointb(:) = xct2(:,2)+vecho(:)/2.0e0_wp
        vecacl(:) = xct2(:,3)-pointa(:)
        vecbo(:)  = xct2(:,4)-pointb(:)
        vecab(:)  = pointb(:)-pointa(:)
        racl      = sqrt(dot_product(vecacl(:), vecacl(:)))
        rbo       = sqrt(dot_product(vecbo(:), vecbo(:)))
        rab       = sqrt(dot_product(vecab(:), vecab(:)))
        anga      = dot_product(vecacl, vecab)/(racl*rab)
        anga      = acos(max(-1.0e0_wp, min(1.0e0_wp, anga)))
        angb      = dot_product(vecbo, vecab)/(rbo*rab)
        angb      = acos(max(-1.0e0_wp, min(1.0e0_wp, angb)))
        angab     = dot_product(vecacl, vecbo)/(racl*rbo)
        angab     = acos(max(-1.0e0_wp, min(1.0e0_wp, angab)))
        vcorr     = dpthcl*dptho/((rab/0.5291772e0_wp)**3)* (cos(angab)-3.0e0_wp*cos(angb)*cos(anga))
        vdp       = ((vt1+vcorr)-vpesmin)*627.509e0_wp - 0.15e0_wp
        vdp       = vdp/23.0605e0_wp
    end subroutine ho_hcl_correction

    ! determine if need to call correction from dipole interaction of HO and HCl
    ! or, switch HHFO in the order of H1-Cl3, H2-O4 if possible
    !> @brief Evaluate before_vpes.
    subroutine before_vpes( xct, lxc )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: xct(3,4), bt(6), xvec(3,6), xct2(3,4)
        real(wp) :: xcom1(3), xcom2(3), rcom, mass(4)
        real(wp) :: r12, r13, r14, r23, r24, r34
        integer  :: lxc, i, lbc, lcd, lde
        mass(1) = 1.00794e0_wp
        mass(2) = mass(1)
        mass(3) = 35.453e0_wp   !HHClO standard weight
        mass(4) = 15.9994e0_wp
        ! lxc=-1, unphysical point
        ! lxc=0,  normal
        ! lxc=1,  HO+HCl, using switch function
        ! lxc=2,  HO+HCl, using dipole interaction correction
        lxc       = 0
        lbc       = 0
        lcd       = 0
        lde       = 0
        xct2      = xct
        xvec(:,1) = xct2(:,2)-xct2(:,1) !  H2-H1
        xvec(:,2) = xct2(:,3)-xct2(:,1) !  Cl-H1
        xvec(:,3) = xct2(:,4)-xct2(:,1) !  O-H1
        xvec(:,4) = xct2(:,3)-xct2(:,2) !  H2-Cl
        xvec(:,5) = xct2(:,4)-xct2(:,2) !  H2-O
        xvec(:,6) = xct2(:,4)-xct2(:,3) !  Cl-O
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
        if( r34 .le. 1.30e0_wp.or.r12 .le. 0.60e0_wp )then
            lxc = -1
            return
        endif
        if( minval(bt) .ge. 2.00e0_wp )then!
            lxc = -1
            return
        endif
        if( r13 .le. 0.75e0_wp.or.r23 .le. 0.75e0_wp )then
            lxc = -1
            return
        endif
        ! H+O+HCl or HO + H +Cl or H+H+Cl+O
        if( r14 .lt. r24.and.r24 .gt. 2.0e0_wp )then
            if( r14 .ge. 2.0e0_wp.or.r23 .ge. 2.0e0_wp )then
                lxc = -1
                return
            endif
        endif
        if( r24 .lt. r14.and.r14 .gt. 2.0e0_wp )then
            if( r24 .ge. 2.0e0_wp.or.r13 .ge. 2.0e0_wp )then
                lxc = -1
                return
            endif
        endif
        ! OH1 and OH2 are both large
        if( r14 .ge. 2.0e0_wp.and.r24 .ge. 2.0e0_wp )then
            lxc = -1
            return
        endif
        ! ClH1 and ClH2 are both small
        if( r13 .le. 1.1e0_wp.and.r23 .le. 1.1e0_wp )then
            lxc = -1
            return
        endif
        !check if HO and HCl were produced
        !    $ (xct2(:,3)-xct2(:,1))
        !    $ (xct2(:,4)-xct2(:,2))
        lxc = 0
    end subroutine before_vpes

    !> @brief Evaluate ho_dipole.
    subroutine ho_dipole( rt, dpt )
        use chempotpy_clh2o_pipnn_parameters, only: dp0ho, dp0ho2, n0ho, r0ho
        use chempotpy_kinds,                  only: wp
        implicit none
        real(wp) :: y2(n0ho), rt, dpt
        y2 = dp0ho2
        if( rt .lt. r0ho(1) ) rt=r0ho(1)
        if( rt .gt. r0ho(n0ho) ) rt=r0ho(n0ho)
        call splint(r0ho, dp0ho, y2, n0ho, rt, dpt)
    end subroutine ho_dipole

    !> @brief Evaluate hcl_dipole.
    subroutine hcl_dipole( rt, dpt )
        use chempotpy_clh2o_pipnn_parameters, only: dp0hcl, dp0hcl2, n0hcl, r0hcl
        use chempotpy_kinds,                  only: wp
        implicit none
        real(wp) :: y2(n0hcl), rt, dpt
        y2 = dp0hcl2
        if( rt .lt. r0hcl(1) ) rt=r0hcl(1)
        if( rt .gt. r0hcl(n0hcl) ) rt=r0hcl(n0hcl)
        call splint(r0hcl, dp0hcl, y2, n0hcl, rt, dpt)
    end subroutine hcl_dipole

    !> @brief Evaluate HO_HCl_pes_spline.
    subroutine ho_hcl_pes_spline( r1, r2, vpes )
        use chempotpy_clh2o_pipnn_parameters, only: exhohcl, exhohcl2, nxhcl, nxho, rxhcl, rxho
        use chempotpy_kinds,                  only: wp
        implicit none
        integer, parameter :: nn = 100
        real(wp)           :: ytmp(nn), y2tmp(nn), y2a(nxhcl,nxho), yytmp(nn)
        real(wp)           :: x1a(nxhcl), x2a(nxho), ya(nxhcl,nxho)
        integer            :: i, j, k, m, n
        real(wp)           :: r1, r2, x1, x2, y, vpes
        x1  = r1
        x2  = r2
        x1a = rxhcl
        x2a = rxho
        ya  = exhohcl
        m   = nxhcl
        n   = nxho
        y2a = exhohcl2
        do j=1, m
            do k=1, n
                ytmp(k)  = ya(j,k)
                y2tmp(k) = y2a(j,k)
            enddo
            call splint(x2a, ytmp, y2tmp, n, x2, yytmp(j))
        enddo
        call spline(x1a, yytmp, m, 1.0e30_wp, 1.0e30_wp, y2tmp)
        call splint(x1a, yytmp, y2tmp, m, x1, y)
        vpes = y
    end subroutine ho_hcl_pes_spline

    !> @brief Evaluate splint.
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
        y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0e0_wp
    end subroutine splint

    !> @brief Evaluate spline.
    subroutine spline( x, y, n, yp1, ypn, y2 )
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter :: nmax = 100
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
end module chempotpy_clh2o_pipnn_pes

!> @brief Expose ClH2O_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_clh2o_pipnn_pes, only: evaluate
    use chempotpy_kinds,           only: wp
    implicit none
    integer, parameter    :: natoms  = 4
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
