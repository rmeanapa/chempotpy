!======================================================================
!       1A' and 1A'' PESs of ArNO at the MRCI+Q level
!======================================================================

!> @brief Evaluate pes.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_geometry,     only: cartesian_distance
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom
    implicit none
    ! number of electronic state
    integer, parameter              :: nstates = 1
    integer, parameter              :: natoms  = 3
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    double precision, intent(in)    :: x(natoms,3)
    double precision, intent(out)   :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out)   :: d(nstates,nstates,natoms,3)
    double precision                :: vtm, vim
    double precision                :: r_no, r_arn, r_aro, cthx
    integer                         :: iatom, idir, j, istate
    !initialize
    vtm = 0.0; p = 0.d0; g = 0.d0; d = 0.d0
    ! x: Ar, N, O order
    r_no  = cartesian_distance(x(2,:), x(3,:))
    r_no  = r_no/legacy_bohr_to_angstrom
    r_arn = cartesian_distance(x(1,:), x(2,:))
    r_arn = r_arn/legacy_bohr_to_angstrom
    r_aro = cartesian_distance(x(1,:), x(3,:))
    r_aro = r_aro/legacy_bohr_to_angstrom
    cthx  = (r_arn**2+r_no**2-r_aro**2)/2d0/r_arn/r_no
    if( igrad == 0 )then
        call potread(path)
        call arnopes(r_no, r_arn, cthx, vtm, vim, 1)
    else
        write (*, *) 'Only energy is available'
    endif
    vtm = vtm/8065.7112013
    do istate=1, nstates
        p(istate) = vtm
    enddo
end subroutine

!> @brief Evaluate ArNOpes.
subroutine arnopes( r_no, r_arn, cthx, vtot, vint, istat )
    implicit none
    double precision :: cthx, r_arn, r_no, s, v, v1, v2, vint, vtot
    integer          :: istat
    call arnopes_sh(r_no, r_arn, cthx, v1, istat)
    call longpes(r_no, r_arn, cthx, v2, istat)
    s    = 0.5d0*(1d0-dtanh(0.8d0*(r_arn-17d0)))
    vtot = s*v1+(1d0-s)*v2
    call pes1d(r_no, v, istat)
    vint = vtot-v
end

!> @brief Evaluate ArNOpes_sh.
subroutine arnopes_sh( r10, r20, cthi, va, istate )
    implicit none
    double precision            :: cth, cthi, r1, r10, r2, r20, th, va, vcut
    integer                     :: istate
    double precision, parameter :: rbohr = 0.5291771
    double precision, parameter :: pi    = 3.141592653589793d0
    common /bvcut/vcut
    ! : r1=r(N-O), r2=r(Ar-N)
    r1  = r10
    r2  = r20
    cth = cthi
    if( cth .gt. 1d0 ) cth=1d0
    if( cth .lt. -1d0 ) cth=-1d0
    th = dacos(cth)*180.0d0/pi
    if( r1 .lt. 1.60d0 ) r1=1.60d0
    if( r1 .gt. 4.0d0 ) r1=4.0d0
    if( r2 .lt. 3.50d0 ) r2=3.50d0
    if( r2 .gt. 40.d0 ) r2=40.d0
    if( istate .ne. 1.and.istate .ne. 2 ) istate=1
    call spl3(r2, th, r1, va, istate)
    !       set zero point at minimum of Ar + NO on A' state
    va = va*219474.63067d0/27.2114d0-31.4930895583566d0
end

!> @brief Evaluate potread.
subroutine potread( path )
    implicit none
    double precision, parameter     :: pi    = 3.141592653589793d0
    integer, parameter              :: nroh1 = 17
    integer, parameter              :: nth   = 13
    integer, parameter              :: nroh2 = 70
    integer, parameter              :: ip    = 2
    integer, parameter              :: m     = 100
    integer, parameter              :: n     = 100
    integer, parameter              :: l     = 100
    character(len=1024), intent(in) :: path
    integer                         :: file_unit_98
    double precision                :: dal, dy1, dyn, pesmin, r1m, r2m, ra1, ra2, roh1, roh2, rra1, &
        & rra2, th1, thm, thth, tth1, vcut, ve, ve1, ve2, vmin1, xt, y, y2
    integer             :: i, ii, ind, j, k
    character(len=1024) :: file_path1
    dimension xt(m), y(m), y2(l)
    data vmin1/-656.852768990000d0/
    data dy1, dyn/1.0d30, 1.0d30/
    common /bvcut/vcut
    common /pesa/pesmin, roh2(nroh1, nth, nroh2), roh1(nroh1), &
        &thth(nth), ve(ip, nroh1, nth, nroh2), ind(nroh1, nth)
    file_path1 = trim(path)//"/ArNO/spline.dat"
    open(newunit=file_unit_98, file=file_path1, status='old')
    !        open(66,file='check.txt')
    vcut   = 12d0 ! in eV
    pesmin = 10000
    !        read(98,*)
    do i=1, nroh1
        do j=1, nth
            ind(i, j) = 0
            do k=1, nroh2
                read(file_unit_98, *)ra1, th1, ra2, ve1, ve2
                if( ra1+ra2 .lt. 0.1 ) goto 7
                if( k .eq. 1 )then
                    rra1 = ra1
                    rra2 = ra2
                    tth1 = th1
                else
                    dal = abs(ra1-rra1)+abs(th1-tth1)
                    if( dal .gt. 0.01 )then
                        write(*, *)' error  ', ii, ra1, rra1, ra2, rra2
                        stop
                    endif
                endif
                if( ve1 .lt. pesmin )then
                    pesmin = ve1
                    r1m    = ra1
                    r2m    = ra2
                    thm    = th1
                endif
                ind(i, j)      = ind(i, j)+1
                roh2(i, j, k)  = ra2
                roh1(i)        = ra1
                thth(j)        = th1
                ve(1, i, j, k) = (ve1-vmin1)*27.2114d0
                ve(2, i, j, k) = (ve2-vmin1)*27.2114d0
            enddo
            7 continue
        enddo
    enddo
    !        write(*,*)r1m,r2m,thm,pesmin
    !        write(*,*)' th1=',(thth(i),i=1,nth)
    !        write(*,*)' roh1=',(roh1(i),i=1,nroh1)
    do i=1, nroh1
        do j=1, nth
            do k=1, ind(i, j)
                if( ve(1, i, j, k) .gt. vcut ) ve(1, i, j, k) = vcut
                if( ve(2, i, j, k) .gt. vcut ) ve(2, i, j, k) = vcut
                !          write(66,1001)roh1(i),thth(j),roh2(i,j,k)
                !     &,ve(1,i,j,k),ve(2,i,j,k)
            enddo
        enddo
    enddo
    close(file_unit_98)
end

!> @brief Evaluate spl3.
subroutine spl3( r1, th, r2, v, istate )
    implicit none
    double precision :: ddhy, ddny, dds1, dds2, dds3, ddty, dhy, dny, ds1, ds2, ds3, dty, dy1, &
        & dyn, h1, h2, h3, pesmin, r1, r1a, r2, roh1, roh2, s1, s2, s3, ss, sss, th, thth, v, &
        & vcut, ve, xt, y, y2, y3, yw2
    integer            :: i, ind, iop, istate, j, k, nh, nr1, nthth
    integer, parameter :: nroh1 = 17
    integer, parameter :: nth   = 13
    integer, parameter :: nroh2 = 70
    integer, parameter :: ip    = 2
    integer, parameter :: m     = 100
    integer, parameter :: n     = 100
    integer, parameter :: l     = 100
    dimension xt(m)
    dimension dty(l), ddty(l), s1(l), ds1(l), dds1(l), h1(l)
    dimension dny(n), ddny(n), s2(l), ds2(l), dds2(l), h2(n)
    dimension dhy(m), ddhy(m), s3(l), ds3(l), dds3(l), h3(m)
    dimension y(m), ss(m), sss(m), y2(l)
    data dy1, dyn/1.0d30, 1.0d30/
    common /bvcut/vcut
    common /pesa/pesmin, roh2(nroh1, nth, nroh2), roh1(nroh1), &
        &thth(nth), ve(ip, nroh1, nth, nroh2), ind(nroh1, nth)
    !        do iop=istate,istate
    iop = istate
    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP&  PRIVATE (i,j,nh,k,xt,y,y2,y3,ss,nthth,yw2)
    !$OMP DO
    do 20 i=1, nroh1
        do 10 j=1, nth
            nh = 0
            do 2 k=1, ind(i, j)
                nh     = nh+1
                xt(nh) = roh2(i, j, k)
                y(nh)  = ve(iop, i, j, k)
            2 continue
            r1a = r1
            call spline(xt, y, nh, dy1, dyn, y2)
            call splint(xt, y, y2, nh, r1a, y3)
            22 continue
            if( y3 .gt. vcut ) y3=vcut
            ss(j) = y3
        10 continue
        nthth = 0
        do 5 j=1, nth
            nthth     = nthth+1
            xt(nthth) = thth(j)
            y(nthth)  = ss(j)
        5 continue
        call spline(xt, y, nthth, 0.0d0, 0.0d0, y2)
        call splint(xt, y, y2, nthth, th, yw2)
        33 continue
        if( yw2 .gt. vcut ) yw2=vcut
        sss(i) = yw2
    20 continue
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    nr1 = nroh1
    !       write(*,*)' vrhh=',(sss(i),i=1,noh)
    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP&  PRIVATE (i)
    !$OMP DO
    do i=1, nroh1
        !       nr1=nr1+1
        xt(i) = roh1(i)
        y(i)  = sss(i)
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    call spline(xt, y, nr1, dy1, dyn, y2)
    call splint(xt, y, y2, nr1, r2, yw2)
    44 continue
    if( yw2 .gt. vcut ) yw2=vcut
    v = yw2
    !       enddo
end

!##################################################################
!# SPLINE ROUTINES
!#            Numerical recipes in fortran
!#            Cambrige University Press
!#            York, 2nd edition, 1992.
!##################################################################
!> @brief Evaluate splint.
subroutine splint( xa, ya, y2a, n, x, y )
    implicit none
    double precision :: a, b, h, x, xa, y, y2a, ya
    integer          :: k, khi, klo, n
    dimension xa(n), y2a(n), ya(n)
    klo = 1
    khi = n
    1 if (khi-klo.gt.1) then
        k = (khi+klo)/2
        if( xa(k) .gt. x )then
            khi = k
        else
            klo = k
        endif
        goto 1
    endif
    h = xa(khi)-xa(klo)
    if( h .eq. 0.0d0 ) write(6, *) 'bad xa input in splint'
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))* (h**2)/6.0d0
end

!##############################################################################
!> @brief Evaluate spline.
subroutine spline( x, y, n, yp1, ypn, y2 )
    implicit none
    double precision :: p, qn, sig, u, un, x, y, y2, yp1, ypn
    integer          :: i, k, n
    dimension x(n), y(n), y2(n)
    integer, parameter :: nmax = 100
    dimension u(nmax)
    if( yp1 .gt. .99d30 )then
        y2(1) = 0.0d0
        u(1)  = 0.0d0
    else
        y2(1) = -0.5d0
        u(1)  = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do 11 i=2, n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1)+2.0d0
        y2(i) = (sig-1.0d0)/p
        u(i)  = (6.0d0*((y(i+1)-y(i))/(x(i+ 1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))- sig*u(i-1))/p
    11 continue
    if( ypn .gt. .99d30 )then
        qn = 0.0d0
        un = 0.0d0
    else
        qn = 0.5d0
        un = (3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
    do 12 k=n-1, 1, -1
        y2(k) = y2(k)*y2(k+1)+u(k)
    12 continue
end

!> @brief Evaluate longpes.
subroutine longpes( r_no, r_arn, cthx, v, istat )
    implicit none
    double precision :: a1, a2, amn, amo, cth, cthx, de, dy1, dyn, r, r1, r2, r_arn, r_no, s, &
        & t, tocm, v, vmin1, x, y, y2
    integer :: i, istat
    dimension x(17), a1(4,17), a2(4,17), y(17), y2(17), de(17)
    data vmin1/-656.852768990000d0/
    !       R in angstrom
    amn  = 14.0030740d0
    amo  = 15.9949146d0
    r2   = r_no
    r1   = r_arn
    t    = amo/(amo+amn)*r2
    s    = dsqrt(r1**2+t**2-2d0*r1*t*cthx) ! in bohr
    r    = s*0.5291772083d0
    cth  = (t**2+s**2-r1**2)/2d0/t/s
    cth  = min(1d0, cth);cth=max(-1d0, cth)
    tocm = 219474.63067d0
    dy1  = 1d30
    dyn  = 1d30
    if( r2 .lt. 1.6d0 ) r2=1.6d0
    if( r2 .gt. 4d0 ) r2=4d0
    x(1)     = 1.60000000000000d0
    x(2)     = 1.80000000000000d0
    x(3)     = 1.90000000000000d0
    x(4)     = 1.95000000000000d0
    x(5)     = 2.00000000000000d0
    x(6)     = 2.05000000000000d0
    x(7)     = 2.10000000000000d0
    x(8)     = 2.15000000000000d0
    x(9)     = 2.20000000000000d0
    x(10)    = 2.2500000000000d0
    x(11)    = 2.3000000000000d0
    x(12)    = 2.4000000000000d0
    x(13)    = 2.6000000000000d0
    x(14)    = 2.8000000000000d0
    x(15)    = 3.0000000000000d0
    x(16)    = 3.5000000000000d0
    x(17)    = 4.0000000000000d0
    a1(1,1)  = -267943.38590546d0
    a1(2,1)  = -38055.810825624d0
    a1(3,1)  = -22967.743334597d0
    a1(4,1)  = 11194.869520824d0
    a1(1,2)  = -251809.17227766d0
    a1(2,2)  = -27650.194725372d0
    a1(3,2)  = -36839.299113681d0
    a1(4,2)  = 8549.8153161124d0
    a1(1,3)  = -251879.04153352d0
    a1(2,3)  = -23032.822374801d0
    a1(3,3)  = -41744.334275868d0
    a1(4,3)  = 6985.8690896270d0
    a1(1,4)  = -250283.76005970d0
    a1(2,4)  = -20018.193520338d0
    a1(3,4)  = -44240.758806033d0
    a1(4,4)  = 5962.8283874165d0
    a1(1,5)  = -254138.28642784d0
    a1(2,5)  = -19723.912012097d0
    a1(3,5)  = -46256.683875458d0
    a1(4,5)  = 6261.8163891098d0
    a1(1,6)  = -255157.48758001d0
    a1(2,6)  = -18048.919914254d0
    a1(3,6)  = -49265.444654343d0
    a1(4,6)  = 5298.6955381574d0
    a1(1,7)  = -255804.28482087d0
    a1(2,7)  = -15646.998402989d0
    a1(3,7)  = -50782.220878825d0
    a1(4,7)  = 5032.7896739693d0
    a1(1,8)  = -257992.33687205d0
    a1(2,8)  = -14507.636299835d0
    a1(3,8)  = -53263.086294840d0
    a1(4,8)  = 5668.8453773915d0
    a1(1,9)  = -262083.75815139d0
    a1(2,9)  = -14037.602789474d0
    a1(3,9)  = -55936.840259583d0
    a1(4,9)  = 5142.3562977947d0
    a1(1,10) = -264567.17397801d0
    a1(2,10) = -13316.118074956d0
    a1(3,10) = -58104.973869090d0
    a1(4,10) = 5459.9310221469d0
    a1(1,11) = -267843.93076864d0
    a1(2,11) = -12807.022293995d0
    a1(3,11) = -60824.874069809d0
    a1(4,11) = 5830.2661402447d0
    a1(1,12) = -276751.23192033d0
    a1(2,12) = -13497.471138696d0
    a1(3,12) = -66418.984554633d0
    a1(4,12) = 5439.6524419885d0
    a1(1,13) = -290143.29086576d0
    a1(2,13) = -11329.693163825d0
    a1(3,13) = -77691.580833384d0
    a1(4,13) = 5266.1638029908d0
    a1(1,14) = -302025.84562829d0
    a1(2,14) = -9660.7649810443d0
    a1(3,14) = -88948.228567055d0
    a1(4,14) = 6062.8633344077d0
    a1(1,15) = -311769.98048147d0
    a1(2,15) = -7392.2131839275d0
    a1(3,15) = -101015.44505332d0
    a1(4,15) = 5729.7740086266d0
    a1(1,16) = -322755.44544985d0
    a1(2,16) = -7292.9293756397d0
    a1(3,16) = -115151.78247626d0
    a1(4,16) = 9316.8888253736d0
    a1(1,17) = -314058.98722717d0
    a1(2,17) = -9845.6172327698d0
    a1(3,17) = -108095.25443221d0
    a1(4,17) = 12298.028717238d0
    a2(1,1)  = -240431.54496463d0
    a2(2,1)  = -24089.327542574d0
    a2(3,1)  = -47953.116074105d0
    a2(4,1)  = -1821.4704857423d0
    a2(1,2)  = -241713.62149381d0
    a2(2,2)  = -19568.897780005d0
    a2(3,2)  = -45838.247416406d0
    a2(4,2)  = 1062.8908146373d0
    a2(1,3)  = -243933.09613888d0
    a2(2,3)  = -17425.641520976d0
    a2(3,3)  = -47462.042744971d0
    a2(4,3)  = 1884.9096687051d0
    a2(1,4)  = -244670.78972127d0
    a2(2,4)  = -15757.181068769d0
    a2(3,4)  = -48668.574025652d0
    a2(4,4)  = 2025.3375037173d0
    a2(1,5)  = -249752.59229822d0
    a2(2,5)  = -16345.228476962d0
    a2(3,5)  = -49790.778143198d0
    a2(4,5)  = 3270.3444683698d0
    a2(1,6)  = -251769.89335023d0
    a2(2,6)  = -15313.617119895d0
    a2(3,6)  = -52458.348815466d0
    a2(4,6)  = 2325.1240929174d0
    a2(1,7)  = -251678.55830509d0
    a2(2,7)  = -12958.861408430d0
    a2(3,7)  = -53659.913259469d0
    a2(4,7)  = 2575.1831916109d0
    a2(1,8)  = -253437.72857935d0
    a2(2,8)  = -12258.648248398d0
    a2(3,8)  = -55795.170116181d0
    a2(4,8)  = 3732.0847356143d0
    a2(1,9)  = -258471.67365255d0
    a2(2,9)  = -12799.769627073d0
    a2(3,9)  = -57508.633013489d0
    a2(4,9)  = 4303.6465924945d0
    a2(1,10) = -262187.31147592d0
    a2(2,10) = -12307.330100049d0
    a2(3,10) = -59805.051173756d0
    a2(4,10) = 4235.9341267133d0
    a2(1,11) = -265266.31321490d0
    a2(2,11) = -11337.914954709d0
    a2(3,11) = -62140.255879242d0
    a2(4,11) = 5747.4771478058d0
    a2(1,12) = -277136.94531031d0
    a2(2,12) = -13950.852000624d0
    a2(3,12) = -67638.158750528d0
    a2(4,12) = 5348.7496751312d0
    a2(1,13) = -287623.13431309d0
    a2(2,13) = -10926.585393021d0
    a2(3,13) = -78513.457663825d0
    a2(4,13) = 5665.0625303400d0
    a2(1,14) = -300470.77056637d0
    a2(2,14) = -10421.152677201d0
    a2(3,14) = -89919.062773954d0
    a2(4,14) = 7178.0270447738d0
    a2(1,15) = -312097.32976855d0
    a2(2,15) = -10237.054620313d0
    a2(3,15) = -102418.68339132d0
    a2(4,15) = 7772.4692422368d0
    a2(1,16) = -322700.43808882d0
    a2(2,16) = -13197.018758343d0
    a2(3,16) = -118420.01108172d0
    a2(4,16) = 13070.129322234d0
    a2(1,17) = -307994.98282431d0
    a2(2,17) = -12733.943399831d0
    a2(3,17) = -111176.93493823d0
    a2(4,17) = 16996.093766663d0
    if( istat .ne. 1.and.istat .ne. 2 ) istat=1
    if( istat .eq. 1 )then
        de(1)  = -656.45779293d0
        de(2)  = -656.72825267d0
        de(3)  = -656.79453152d0
        de(4)  = -656.81626228d0
        de(5)  = -656.83198474d0
        de(6)  = -656.84267992d0
        de(7)  = -656.84918525d0
        de(8)  = -656.85221531d0
        de(9)  = -656.85237944d0
        de(10) = -656.85019678d0
        de(11) = -656.84610901d0
        de(12) = -656.83366196d0
        de(13) = -656.79902012d0
        de(14) = -656.76032202d0
        de(15) = -656.72354209d0
        de(16) = -656.65642703d0
        de(17) = -656.62721422d0
        do i=1, 17
            call vfuncs(cth, r, a1(:,i), v)
            y(i) = v+(de(i)-vmin1)*tocm
        enddo
    else
        de(1)  = -656.45778154d0
        de(2)  = -656.72823764d0
        de(3)  = -656.79451533d0
        de(4)  = -656.81624562d0
        de(5)  = -656.83196768d0
        de(6)  = -656.84266252d0
        de(7)  = -656.84916757d0
        de(8)  = -656.85219739d0
        de(9)  = -656.85236132d0
        de(10) = -656.85017849d0
        de(11) = -656.84609058d0
        de(12) = -656.83364329d0
        de(13) = -656.79900106d0
        de(14) = -656.76030267d0
        de(15) = -656.72352380d0
        de(16) = -656.65639961d0
        de(17) = -656.62718628d0
        do i=1, 17
            call vfuncs(cth, r, a2(:,i), v)
            y(i) = v+(de(i)-vmin1)*tocm
        enddo
    endif
    call spline(x, y, 17, dy1, dyn, y2)
    call splint(x, y, y2, 17, r2, v)
    !       set zero point at minimum of Ar + NO on A' state
    v = v-31.4930895583566d0
end

!***********************************************************************
!> @brief Evaluate VFUNCS.
subroutine vfuncs( theta, r, c, v )
    implicit none
    double precision :: afunc, c, cth, f, r, theta, v
    integer          :: i
    dimension afunc(4), c(4)
    cth = theta
    call fval(r, cth, f)
    afunc(1) = f/r**6
    afunc(2) = cth*f/r**6
    afunc(3) = 0.5d0*(3d0*cth**2-1d0)*f/r**6
    afunc(4) = 0.5d0*(5d0*cth**3-3d0*cth)*f/r**6
    v        = 0d0
    do i=1, 4
        v = v+afunc(i)*c(i)
    enddo
end

!> @brief Evaluate Bval.
subroutine bval( cth, b )
    implicit none
    double precision :: cth, b, b0, b1, b2
    b0 = -3.0318d0
    b1 = 0.1299d0
    b2 = 0.0564d0
    b  = b0+b1*cth+b2*0.5d0*(3d0*cth**2-1d0)
end

!> @brief Evaluate Fval.
subroutine fval( r, cth, f )
    implicit none
    double precision :: r, cth, f, b
    call bval(cth, b)
    b = abs(b)*r
    f = 1d0-dexp(-b)*(1d0+b+b**2/2d0+b**3/6d0+b**4/24d0+b**5/120d0+ b**6/720d0)
end

!> @brief Evaluate pes1D.
subroutine pes1d( r_no, v, istat )
    implicit none
    double precision :: a1, a2, de, dy1, dyn, r2, r_no, tocm, v, vmin1, x, y, y2
    integer          :: i, istat
    dimension x(17), a1(4,17), a2(4,17), y(17), y2(17), de(17)
    data vmin1/-656.852768990000d0/
    r2   = r_no
    tocm = 219474.63067d0
    dy1  = 1d30
    dyn  = 1d30
    if( r2 .lt. 1.6d0 ) r2=1.6d0
    if( r2 .gt. 4d0 ) r2=4d0
    x(1)  = 1.60000000000000d0
    x(2)  = 1.80000000000000d0
    x(3)  = 1.90000000000000d0
    x(4)  = 1.95000000000000d0
    x(5)  = 2.00000000000000d0
    x(6)  = 2.05000000000000d0
    x(7)  = 2.10000000000000d0
    x(8)  = 2.15000000000000d0
    x(9)  = 2.20000000000000d0
    x(10) = 2.2500000000000d0
    x(11) = 2.3000000000000d0
    x(12) = 2.4000000000000d0
    x(13) = 2.6000000000000d0
    x(14) = 2.8000000000000d0
    x(15) = 3.0000000000000d0
    x(16) = 3.5000000000000d0
    x(17) = 4.0000000000000d0
    if( istat .ne. 1.and.istat .ne. 2 ) istat=1
    if( istat .eq. 1 )then
        de(1)  = -656.45779293d0
        de(2)  = -656.72825267d0
        de(3)  = -656.79453152d0
        de(4)  = -656.81626228d0
        de(5)  = -656.83198474d0
        de(6)  = -656.84267992d0
        de(7)  = -656.84918525d0
        de(8)  = -656.85221531d0
        de(9)  = -656.85237944d0
        de(10) = -656.85019678d0
        de(11) = -656.84610901d0
        de(12) = -656.83366196d0
        de(13) = -656.79902012d0
        de(14) = -656.76032202d0
        de(15) = -656.72354209d0
        de(16) = -656.65642703d0
        de(17) = -656.62721422d0
        do i=1, 17
            y(i) = (de(i)-vmin1)*tocm
        enddo
    else
        de(1)  = -656.45778154d0
        de(2)  = -656.72823764d0
        de(3)  = -656.79451533d0
        de(4)  = -656.81624562d0
        de(5)  = -656.83196768d0
        de(6)  = -656.84266252d0
        de(7)  = -656.84916757d0
        de(8)  = -656.85219739d0
        de(9)  = -656.85236132d0
        de(10) = -656.85017849d0
        de(11) = -656.84609058d0
        de(12) = -656.83364329d0
        de(13) = -656.79900106d0
        de(14) = -656.76030267d0
        de(15) = -656.72352380d0
        de(16) = -656.65639961d0
        de(17) = -656.62718628d0
        do i=1, 17
            y(i) = (de(i)-vmin1)*tocm
        enddo
    endif
    call spline(x, y, 17, dy1, dyn, y2)
    call splint(x, y, y2, 17, r2, v)
    !       set zero point at minimum of Ar + NO on A' state
    v = v-31.4930895583566d0
end
