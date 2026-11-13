!       1A' and 1A'' PESs of ArNO at the MRCI+Q level

!> @brief Implement ArNO_Ap_PIPNN, a potential energy surface for ArNO.
module chempotpy_arno_ap_pipnn_pes
    use chempotpy_arno_spline_parameters, only: bvcut_vcut, pesa_ind, pesa_roh1, pesa_roh2, pesa_thth, pesa_ve
    use chempotpy_kinds,                  only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the ArNO_Ap_PIPNN potential energy surface for ArNO.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,    only: electronvolt_to_inverse_centimeter
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_units,    only: bohr_to_angstrom
        use chempotpy_kinds,    only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: vtm, vim
        real(wp)              :: r_no, r_arn, r_aro, cthx
        integer               :: iatom, idir, j, istate
        !initialize
        vtm = 0.0; p = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        ! x: Ar, N, O order
        r_no  = cartesian_distance(x(2,:), x(3,:))
        r_no  = r_no/bohr_to_angstrom
        r_arn = cartesian_distance(x(1,:), x(2,:))
        r_arn = r_arn/bohr_to_angstrom
        r_aro = cartesian_distance(x(1,:), x(3,:))
        r_aro = r_aro/bohr_to_angstrom
        cthx  = (r_arn**2+r_no**2-r_aro**2)/2e0_wp/r_arn/r_no
        if( igrad == 0 )then
            call arnopes(r_no, r_arn, cthx, vtm, vim, 1)
        else
            write (*, *) 'Only energy is available'
        endif
        vtm = vtm/electronvolt_to_inverse_centimeter
        do istate=1, nstates
            p(istate) = vtm
        enddo

        contains

        !> @brief Evaluate ArNOpes.
        subroutine arnopes( r_no, r_arn, cthx, vtot, vint, istat )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: cthx, r_arn, r_no, s, v, v1, v2, vint, vtot
            integer  :: istat
            call arnopes_sh(r_no, r_arn, cthx, v1, istat)
            call longpes(r_no, r_arn, cthx, v2, istat)
            s    = 0.5e0_wp*(1e0_wp-tanh(0.8e0_wp*(r_arn-17e0_wp)))
            vtot = s*v1+(1e0_wp-s)*v2
            call pes1d(r_no, v, istat)
            vint = vtot-v
        end
        !> @brief Evaluate ArNOpes_sh.
        subroutine arnopes_sh( r10, r20, cthi, va, istate )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: cth, cthi, r1, r10, r2, r20, th, va
            integer             :: istate
            real(wp), parameter :: rbohr = 0.5291771
            real(wp), parameter :: pi    = 3.141592653589793e0_wp
            ! : r1=r(N-O), r2=r(Ar-N)
            r1  = r10
            r2  = r20
            cth = cthi
            if( cth .gt. 1e0_wp ) cth=1e0_wp
            if( cth .lt. -1e0_wp ) cth=-1e0_wp
            th = acos(cth)*180.0e0_wp/pi
            if( r1 .lt. 1.60e0_wp ) r1=1.60e0_wp
            if( r1 .gt. 4.0e0_wp ) r1=4.0e0_wp
            if( r2 .lt. 3.50e0_wp ) r2=3.50e0_wp
            if( r2 .gt. 40.e0_wp ) r2=40.e0_wp
            if( istate .ne. 1.and.istate .ne. 2 ) istate=1
            call spl3(r2, th, r1, va, istate)
            !       set zero point at minimum of Ar + NO on A' state
            va = va*219474.63067e0_wp/27.2114e0_wp-31.4930895583566e0_wp
        end
        !> @brief Evaluate spl3.
        subroutine spl3( r1, th, r2, v, istate )
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: nroh1 = 17
            integer, parameter :: nth   = 13
            integer, parameter :: nroh2 = 70
            integer, parameter :: ip    = 2
            integer, parameter :: m     = 100
            integer, parameter :: n     = 100
            integer, parameter :: l     = 100
            real(wp)           :: ddhy(m), ddny(n), dds1(l), dds2(l), dds3(l), ddty(l), dhy(m), dny(n), ds1(l), ds2(l), ds3(l), &
                & dty(l), dy1, dyn, h1(l), h2(n), h3(m), r1, r1a, r2, s1(l), s2(l), s3(l), ss(m), sss(m), th, v, xt(m), y(m), y2(l), &
                & y3, yw2
            integer :: i, iop, istate, j, k, nh, nr1, nthth
            !        do iop=istate,istate
            dy1 = 1.0e30_wp; dyn = 1.0e30_wp
            iop = istate
            !$OMP PARALLEL DEFAULT(SHARED)
            !$OMP&  PRIVATE (i,j,nh,k,xt,y,y2,y3,ss,nthth,yw2)
            !$OMP DO
            do i=1, nroh1
                do j=1, nth
                    nh = 0
                    do k=1, pesa_ind(i, j)
                        nh     = nh+1
                        xt(nh) = pesa_roh2(i, j, k)
                        y(nh)  = pesa_ve(iop, i, j, k)
                    enddo
                    r1a = r1
                    call spline(xt, y, nh, dy1, dyn, y2)
                    call splint(xt, y, y2, nh, r1a, y3)
                    if( y3 .gt. bvcut_vcut ) y3=bvcut_vcut
                    ss(j) = y3
                enddo
                nthth = 0
                do j=1, nth
                    nthth     = nthth+1
                    xt(nthth) = pesa_thth(j)
                    y(nthth)  = ss(j)
                enddo
                call spline(xt, y, nthth, 0.0e0_wp, 0.0e0_wp, y2)
                call splint(xt, y, y2, nthth, th, yw2)
                if( yw2 .gt. bvcut_vcut ) yw2=bvcut_vcut
                sss(i) = yw2
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
            nr1 = nroh1
            !       write(*,*)' vrhh=',(sss(i),i=1,noh)
            !$OMP PARALLEL DEFAULT(SHARED)
            !$OMP&  PRIVATE (i)
            !$OMP DO
            do i=1, nroh1
                !       nr1=nr1+1
                xt(i) = pesa_roh1(i)
                y(i)  = sss(i)
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
            call spline(xt, y, nr1, dy1, dyn, y2)
            call splint(xt, y, y2, nr1, r2, yw2)
            if( yw2 .gt. bvcut_vcut ) yw2=bvcut_vcut
            v = yw2
            !       enddo
        end
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
            y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))* (h**2)/6.0e0_wp
        end
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
                u(i)  = (6.0e0_wp*((y(i+1)-y(i))/(x(i+ 1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))- sig* &
                    & u(i- 1))/p
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
        end
        !> @brief Evaluate longpes.
        subroutine longpes( r_no, r_arn, cthx, v, istat )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a1(4,17), a2(4,17), amn, amo, cth, cthx, de(17), dy1, dyn, r, r1, r2, r_arn, r_no, s, t, tocm, v, &
                & vmin1, x(17), y(17), y2(17)
            integer :: i, istat
            !       R in angstrom
            vmin1 = -656.852768990000e0_wp
            amn   = 14.0030740e0_wp
            amo   = 15.9949146e0_wp
            r2    = r_no
            r1    = r_arn
            t     = amo/(amo+amn)*r2
            s     = sqrt(r1**2+t**2-2e0_wp*r1*t*cthx) ! in bohr
            r     = s*0.5291772083e0_wp
            cth   = (t**2+s**2-r1**2)/2e0_wp/t/s
            cth   = min(1e0_wp, cth);cth=max(-1e0_wp, cth)
            tocm  = 219474.63067e0_wp
            dy1   = 1e30_wp
            dyn   = 1e30_wp
            if( r2 .lt. 1.6e0_wp ) r2=1.6e0_wp
            if( r2 .gt. 4e0_wp ) r2=4e0_wp
            x(1)     = 1.60000000000000e0_wp
            x(2)     = 1.80000000000000e0_wp
            x(3)     = 1.90000000000000e0_wp
            x(4)     = 1.95000000000000e0_wp
            x(5)     = 2.00000000000000e0_wp
            x(6)     = 2.05000000000000e0_wp
            x(7)     = 2.10000000000000e0_wp
            x(8)     = 2.15000000000000e0_wp
            x(9)     = 2.20000000000000e0_wp
            x(10)    = 2.2500000000000e0_wp
            x(11)    = 2.3000000000000e0_wp
            x(12)    = 2.4000000000000e0_wp
            x(13)    = 2.6000000000000e0_wp
            x(14)    = 2.8000000000000e0_wp
            x(15)    = 3.0000000000000e0_wp
            x(16)    = 3.5000000000000e0_wp
            x(17)    = 4.0000000000000e0_wp
            a1(1,1)  = -267943.38590546e0_wp
            a1(2,1)  = -38055.810825624e0_wp
            a1(3,1)  = -22967.743334597e0_wp
            a1(4,1)  = 11194.869520824e0_wp
            a1(1,2)  = -251809.17227766e0_wp
            a1(2,2)  = -27650.194725372e0_wp
            a1(3,2)  = -36839.299113681e0_wp
            a1(4,2)  = 8549.8153161124e0_wp
            a1(1,3)  = -251879.04153352e0_wp
            a1(2,3)  = -23032.822374801e0_wp
            a1(3,3)  = -41744.334275868e0_wp
            a1(4,3)  = 6985.8690896270e0_wp
            a1(1,4)  = -250283.76005970e0_wp
            a1(2,4)  = -20018.193520338e0_wp
            a1(3,4)  = -44240.758806033e0_wp
            a1(4,4)  = 5962.8283874165e0_wp
            a1(1,5)  = -254138.28642784e0_wp
            a1(2,5)  = -19723.912012097e0_wp
            a1(3,5)  = -46256.683875458e0_wp
            a1(4,5)  = 6261.8163891098e0_wp
            a1(1,6)  = -255157.48758001e0_wp
            a1(2,6)  = -18048.919914254e0_wp
            a1(3,6)  = -49265.444654343e0_wp
            a1(4,6)  = 5298.6955381574e0_wp
            a1(1,7)  = -255804.28482087e0_wp
            a1(2,7)  = -15646.998402989e0_wp
            a1(3,7)  = -50782.220878825e0_wp
            a1(4,7)  = 5032.7896739693e0_wp
            a1(1,8)  = -257992.33687205e0_wp
            a1(2,8)  = -14507.636299835e0_wp
            a1(3,8)  = -53263.086294840e0_wp
            a1(4,8)  = 5668.8453773915e0_wp
            a1(1,9)  = -262083.75815139e0_wp
            a1(2,9)  = -14037.602789474e0_wp
            a1(3,9)  = -55936.840259583e0_wp
            a1(4,9)  = 5142.3562977947e0_wp
            a1(1,10) = -264567.17397801e0_wp
            a1(2,10) = -13316.118074956e0_wp
            a1(3,10) = -58104.973869090e0_wp
            a1(4,10) = 5459.9310221469e0_wp
            a1(1,11) = -267843.93076864e0_wp
            a1(2,11) = -12807.022293995e0_wp
            a1(3,11) = -60824.874069809e0_wp
            a1(4,11) = 5830.2661402447e0_wp
            a1(1,12) = -276751.23192033e0_wp
            a1(2,12) = -13497.471138696e0_wp
            a1(3,12) = -66418.984554633e0_wp
            a1(4,12) = 5439.6524419885e0_wp
            a1(1,13) = -290143.29086576e0_wp
            a1(2,13) = -11329.693163825e0_wp
            a1(3,13) = -77691.580833384e0_wp
            a1(4,13) = 5266.1638029908e0_wp
            a1(1,14) = -302025.84562829e0_wp
            a1(2,14) = -9660.7649810443e0_wp
            a1(3,14) = -88948.228567055e0_wp
            a1(4,14) = 6062.8633344077e0_wp
            a1(1,15) = -311769.98048147e0_wp
            a1(2,15) = -7392.2131839275e0_wp
            a1(3,15) = -101015.44505332e0_wp
            a1(4,15) = 5729.7740086266e0_wp
            a1(1,16) = -322755.44544985e0_wp
            a1(2,16) = -7292.9293756397e0_wp
            a1(3,16) = -115151.78247626e0_wp
            a1(4,16) = 9316.8888253736e0_wp
            a1(1,17) = -314058.98722717e0_wp
            a1(2,17) = -9845.6172327698e0_wp
            a1(3,17) = -108095.25443221e0_wp
            a1(4,17) = 12298.028717238e0_wp
            a2(1,1)  = -240431.54496463e0_wp
            a2(2,1)  = -24089.327542574e0_wp
            a2(3,1)  = -47953.116074105e0_wp
            a2(4,1)  = -1821.4704857423e0_wp
            a2(1,2)  = -241713.62149381e0_wp
            a2(2,2)  = -19568.897780005e0_wp
            a2(3,2)  = -45838.247416406e0_wp
            a2(4,2)  = 1062.8908146373e0_wp
            a2(1,3)  = -243933.09613888e0_wp
            a2(2,3)  = -17425.641520976e0_wp
            a2(3,3)  = -47462.042744971e0_wp
            a2(4,3)  = 1884.9096687051e0_wp
            a2(1,4)  = -244670.78972127e0_wp
            a2(2,4)  = -15757.181068769e0_wp
            a2(3,4)  = -48668.574025652e0_wp
            a2(4,4)  = 2025.3375037173e0_wp
            a2(1,5)  = -249752.59229822e0_wp
            a2(2,5)  = -16345.228476962e0_wp
            a2(3,5)  = -49790.778143198e0_wp
            a2(4,5)  = 3270.3444683698e0_wp
            a2(1,6)  = -251769.89335023e0_wp
            a2(2,6)  = -15313.617119895e0_wp
            a2(3,6)  = -52458.348815466e0_wp
            a2(4,6)  = 2325.1240929174e0_wp
            a2(1,7)  = -251678.55830509e0_wp
            a2(2,7)  = -12958.861408430e0_wp
            a2(3,7)  = -53659.913259469e0_wp
            a2(4,7)  = 2575.1831916109e0_wp
            a2(1,8)  = -253437.72857935e0_wp
            a2(2,8)  = -12258.648248398e0_wp
            a2(3,8)  = -55795.170116181e0_wp
            a2(4,8)  = 3732.0847356143e0_wp
            a2(1,9)  = -258471.67365255e0_wp
            a2(2,9)  = -12799.769627073e0_wp
            a2(3,9)  = -57508.633013489e0_wp
            a2(4,9)  = 4303.6465924945e0_wp
            a2(1,10) = -262187.31147592e0_wp
            a2(2,10) = -12307.330100049e0_wp
            a2(3,10) = -59805.051173756e0_wp
            a2(4,10) = 4235.9341267133e0_wp
            a2(1,11) = -265266.31321490e0_wp
            a2(2,11) = -11337.914954709e0_wp
            a2(3,11) = -62140.255879242e0_wp
            a2(4,11) = 5747.4771478058e0_wp
            a2(1,12) = -277136.94531031e0_wp
            a2(2,12) = -13950.852000624e0_wp
            a2(3,12) = -67638.158750528e0_wp
            a2(4,12) = 5348.7496751312e0_wp
            a2(1,13) = -287623.13431309e0_wp
            a2(2,13) = -10926.585393021e0_wp
            a2(3,13) = -78513.457663825e0_wp
            a2(4,13) = 5665.0625303400e0_wp
            a2(1,14) = -300470.77056637e0_wp
            a2(2,14) = -10421.152677201e0_wp
            a2(3,14) = -89919.062773954e0_wp
            a2(4,14) = 7178.0270447738e0_wp
            a2(1,15) = -312097.32976855e0_wp
            a2(2,15) = -10237.054620313e0_wp
            a2(3,15) = -102418.68339132e0_wp
            a2(4,15) = 7772.4692422368e0_wp
            a2(1,16) = -322700.43808882e0_wp
            a2(2,16) = -13197.018758343e0_wp
            a2(3,16) = -118420.01108172e0_wp
            a2(4,16) = 13070.129322234e0_wp
            a2(1,17) = -307994.98282431e0_wp
            a2(2,17) = -12733.943399831e0_wp
            a2(3,17) = -111176.93493823e0_wp
            a2(4,17) = 16996.093766663e0_wp
            if( istat .ne. 1.and.istat .ne. 2 ) istat=1
            if( istat .eq. 1 )then
                de(1)  = -656.45779293e0_wp
                de(2)  = -656.72825267e0_wp
                de(3)  = -656.79453152e0_wp
                de(4)  = -656.81626228e0_wp
                de(5)  = -656.83198474e0_wp
                de(6)  = -656.84267992e0_wp
                de(7)  = -656.84918525e0_wp
                de(8)  = -656.85221531e0_wp
                de(9)  = -656.85237944e0_wp
                de(10) = -656.85019678e0_wp
                de(11) = -656.84610901e0_wp
                de(12) = -656.83366196e0_wp
                de(13) = -656.79902012e0_wp
                de(14) = -656.76032202e0_wp
                de(15) = -656.72354209e0_wp
                de(16) = -656.65642703e0_wp
                de(17) = -656.62721422e0_wp
                do i=1, 17
                    call vfuncs(cth, r, a1(:,i), v)
                    y(i) = v+(de(i)-vmin1)*tocm
                enddo
            else
                de(1)  = -656.45778154e0_wp
                de(2)  = -656.72823764e0_wp
                de(3)  = -656.79451533e0_wp
                de(4)  = -656.81624562e0_wp
                de(5)  = -656.83196768e0_wp
                de(6)  = -656.84266252e0_wp
                de(7)  = -656.84916757e0_wp
                de(8)  = -656.85219739e0_wp
                de(9)  = -656.85236132e0_wp
                de(10) = -656.85017849e0_wp
                de(11) = -656.84609058e0_wp
                de(12) = -656.83364329e0_wp
                de(13) = -656.79900106e0_wp
                de(14) = -656.76030267e0_wp
                de(15) = -656.72352380e0_wp
                de(16) = -656.65639961e0_wp
                de(17) = -656.62718628e0_wp
                do i=1, 17
                    call vfuncs(cth, r, a2(:,i), v)
                    y(i) = v+(de(i)-vmin1)*tocm
                enddo
            endif
            call spline(x, y, 17, dy1, dyn, y2)
            call splint(x, y, y2, 17, r2, v)
            !       set zero point at minimum of Ar + NO on A' state
            v = v-31.4930895583566e0_wp
        end
        !> @brief Evaluate VFUNCS.
        subroutine vfuncs( theta, r, c, v )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: afunc(4), c(4), cth, f, r, theta, v
            integer  :: i
            cth = theta
            call fval(r, cth, f)
            afunc(1) = f/r**6
            afunc(2) = cth*f/r**6
            afunc(3) = 0.5e0_wp*(3e0_wp*cth**2-1e0_wp)*f/r**6
            afunc(4) = 0.5e0_wp*(5e0_wp*cth**3-3e0_wp*cth)*f/r**6
            v        = 0e0_wp
            do i=1, 4
                v = v+afunc(i)*c(i)
            enddo
        end
        !> @brief Evaluate Bval.
        subroutine bval( cth, b )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: cth, b, b0, b1, b2
            b0 = -3.0318e0_wp
            b1 = 0.1299e0_wp
            b2 = 0.0564e0_wp
            b  = b0+b1*cth+b2*0.5e0_wp*(3e0_wp*cth**2-1e0_wp)
        end
        !> @brief Evaluate Fval.
        subroutine fval( r, cth, f )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: r, cth, f, b
            call bval(cth, b)
            b = abs(b)*r
            f = 1e0_wp-exp(-b)*(1e0_wp+b+b**2/2e0_wp+b**3/6e0_wp+b**4/24e0_wp+b**5/120e0_wp+ b**6/720e0_wp)
        end
        !> @brief Evaluate pes1D.
        subroutine pes1d( r_no, v, istat )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a1(4,17), a2(4,17), de(17), dy1, dyn, r2, r_no, tocm, v, vmin1, x(17), y(17), y2(17)
            integer  :: i, istat
            vmin1 = -656.852768990000e0_wp
            r2    = r_no
            tocm  = 219474.63067e0_wp
            dy1   = 1e30_wp
            dyn   = 1e30_wp
            if( r2 .lt. 1.6e0_wp ) r2=1.6e0_wp
            if( r2 .gt. 4e0_wp ) r2=4e0_wp
            x(1)  = 1.60000000000000e0_wp
            x(2)  = 1.80000000000000e0_wp
            x(3)  = 1.90000000000000e0_wp
            x(4)  = 1.95000000000000e0_wp
            x(5)  = 2.00000000000000e0_wp
            x(6)  = 2.05000000000000e0_wp
            x(7)  = 2.10000000000000e0_wp
            x(8)  = 2.15000000000000e0_wp
            x(9)  = 2.20000000000000e0_wp
            x(10) = 2.2500000000000e0_wp
            x(11) = 2.3000000000000e0_wp
            x(12) = 2.4000000000000e0_wp
            x(13) = 2.6000000000000e0_wp
            x(14) = 2.8000000000000e0_wp
            x(15) = 3.0000000000000e0_wp
            x(16) = 3.5000000000000e0_wp
            x(17) = 4.0000000000000e0_wp
            if( istat .ne. 1.and.istat .ne. 2 ) istat=1
            if( istat .eq. 1 )then
                de(1)  = -656.45779293e0_wp
                de(2)  = -656.72825267e0_wp
                de(3)  = -656.79453152e0_wp
                de(4)  = -656.81626228e0_wp
                de(5)  = -656.83198474e0_wp
                de(6)  = -656.84267992e0_wp
                de(7)  = -656.84918525e0_wp
                de(8)  = -656.85221531e0_wp
                de(9)  = -656.85237944e0_wp
                de(10) = -656.85019678e0_wp
                de(11) = -656.84610901e0_wp
                de(12) = -656.83366196e0_wp
                de(13) = -656.79902012e0_wp
                de(14) = -656.76032202e0_wp
                de(15) = -656.72354209e0_wp
                de(16) = -656.65642703e0_wp
                de(17) = -656.62721422e0_wp
                do i=1, 17
                    y(i) = (de(i)-vmin1)*tocm
                enddo
            else
                de(1)  = -656.45778154e0_wp
                de(2)  = -656.72823764e0_wp
                de(3)  = -656.79451533e0_wp
                de(4)  = -656.81624562e0_wp
                de(5)  = -656.83196768e0_wp
                de(6)  = -656.84266252e0_wp
                de(7)  = -656.84916757e0_wp
                de(8)  = -656.85219739e0_wp
                de(9)  = -656.85236132e0_wp
                de(10) = -656.85017849e0_wp
                de(11) = -656.84609058e0_wp
                de(12) = -656.83364329e0_wp
                de(13) = -656.79900106e0_wp
                de(14) = -656.76030267e0_wp
                de(15) = -656.72352380e0_wp
                de(16) = -656.65639961e0_wp
                de(17) = -656.62718628e0_wp
                do i=1, 17
                    y(i) = (de(i)-vmin1)*tocm
                enddo
            endif
            call spline(x, y, 17, dy1, dyn, y2)
            call splint(x, y, y2, 17, r2, v)
            !       set zero point at minimum of Ar + NO on A' state
            v = v-31.4930895583566e0_wp
        end
    end subroutine evaluate
end module chempotpy_arno_ap_pipnn_pes

!> @brief Expose ArNO_Ap_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_arno_ap_pipnn_pes, only: evaluate
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
