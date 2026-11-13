!> @brief Implement N4_BDTC_2014, a potential energy surface for N4.
module chempotpy_n4_bdtc_2014_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the N4_BDTC_2014 potential energy surface for N4.
    subroutine evaluate( x, igrad, path, p, g, d )
        use chempotpy_units, only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter              :: nstates = 1
        integer, parameter              :: natoms  = 4
        character(len=1024), intent(in) :: path
        integer, intent(in)             :: igrad
        real(wp), intent(in)            :: x(natoms,3)
        real(wp), intent(out)           :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
        real(wp)                        :: x1(4), y1(4), z1(4)
        real(wp)                        :: dx1(4), dy1(4), dz1(4)
        real(wp)                        :: v
        integer                         :: i, j, k, l
        integer                         :: istate, iatom, idir
        !initialize
        v     = 0.e0_wp; p = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        x1(:) = x(:,1)/bohr_to_angstrom
        y1(:) = x(:,2)/bohr_to_angstrom
        z1(:) = x(:,3)/bohr_to_angstrom
        call assemble_model_parameters()
        call pot(x1, y1, z1, v, dx1, dy1, dz1)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = v*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = v*hartree_to_electronvolt
            enddo
            do istate=1, nstates
                g(istate,:,1) = dx1(:)*hartree_per_bohr_to_ev_per_angstrom
                g(istate,:,2) = dy1(:)*hartree_per_bohr_to_ev_per_angstrom
                g(istate,:,3) = dz1(:)*hartree_per_bohr_to_ev_per_angstrom
            enddo
        else
            write(*, *) "Only energy and gradients are available"
        endif
    end subroutine evaluate
end module chempotpy_n4_bdtc_2014_pes

!> @brief Expose N4_BDTC_2014 through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_n4_bdtc_2014_pes, only: evaluate
    use chempotpy_kinds,            only: wp
    implicit none
    integer, parameter              :: natoms  = 4
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, path, p, g, d)
end subroutine pes

module falcons_vars_eval
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: aji, qi
    real(wp), dimension(:, :), allocatable :: qi, aji
end module falcons_vars_eval

module falcons_vars_main
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: aa, cgmeth, cmap, crc, crcmult, dedy, e, epsnu, epsu, epw, eta, i, ind, j, k, l, m, &
        & mcrc, moltype, n, nall, nitcrc, nma, nmincrc, nn, np, nplmax, ntar, ntarmax, ntarmin, &
        & pesform, pnu, ps, pu, pwonly, q, qeq, qp, qps, qr, qs, r, rwfa, rwfb, sd, systo, t, x, xr, &
        & y
    integer :: i, j, k, l, nn, ind, n, m, eta, np, moltype, pesform, &
        & nma, nplmax, nall, nmincrc, mcrc, qps, pwonly, cgmeth, &
        & ntarmin, ntarmax, ntar, crcmult, nitcrc
    integer, dimension(6)     :: sd
    integer, dimension(4, 3)  :: cmap
    integer, dimension(24, 6) :: qp
    real(wp)                  :: pnu, epsnu, ps, pu, epsu, e, epw, systo, r, rwfa, rwfb
    real(wp), dimension(6)    :: t, qr, qs, q, xr, x, aa, qeq
    real(wp), dimension(4, 3) :: y, dedy
    real(wp), dimension(:), allocatable :: crc
end module falcons_vars_main

!> @brief Assemble derived N4 BDTC arrays from immutable model parameters.
subroutine assemble_model_parameters( )
    use chempotpy_n4_bdtc_parameters, only: coefficient_value, parameter_crc => crc, &
        & parameter_crcmult => crcmult, parameter_epsnu => epsnu, parameter_eta => eta, &
        & parameter_equilibrium_distances => equilibrium_distances, geometry_value, &
        & parameter_m => m, parameter_mcrc => mcrc, parameter_moltype => moltype, &
        & parameter_morse_widths => morse_widths, parameter_n => n, parameter_nitcrc => nitcrc, &
        & parameter_nma => nma, parameter_nmincrc => nmincrc, parameter_nplmax => nplmax, &
        & parameter_ntar => ntar, parameter_ntarmax => ntarmax, parameter_ntarmin => ntarmin, &
        & parameter_pesform => pesform, parameter_pnu => pnu, parameter_rwfa => rwfa, &
        & parameter_rwfb => rwfb
    use falcons_vars_eval, only: aji, qi
    use falcons_vars_main, only: aa, cgmeth, cmap, crc, crcmult, dedy, e, epsnu, epsu, epw, eta, i, &
        & ind, j, k, l, m, mcrc, moltype, n, nall, nitcrc, nma, nmincrc, nn, np, nplmax, ntar, &
        & ntarmax, ntarmin, pesform, pnu, ps, pu, pwonly, q, qeq, qp, qps, qr, qs, r, rwfa, rwfb, sd, &
        & systo, t, x, xr, y
    use chempotpy_kinds, only: wp
    implicit none
    integer :: oridval, maskval
    if( allocated(qi) ) deallocate(qi)
    if( allocated(aji) ) deallocate(aji)
    if( allocated(crc) ) deallocate(crc)
    cmap(1, :) = (/1, 2, 3/)
    cmap(2, :) = (/1, 4, 5/)
    cmap(3, :) = (/2, 4, 6/)
    cmap(4, :) = (/3, 5, 6/)
    qps        = 24
    qp(1, :)   = (/1, 2, 3, 4, 5, 6/)
    qp(2, :)   = (/1, 3, 2, 5, 4, 6/)
    qp(3, :)   = (/2, 1, 3, 4, 6, 5/)
    qp(4, :)   = (/2, 3, 1, 6, 4, 5/)
    qp(5, :)   = (/3, 1, 2, 5, 6, 4/)
    qp(6, :)   = (/3, 2, 1, 6, 5, 4/)
    qp(7, :)   = (/1, 4, 5, 2, 3, 6/)
    qp(8, :)   = (/1, 5, 4, 3, 2, 6/)
    qp(9, :)   = (/4, 1, 5, 2, 6, 3/)
    qp(10, :)  = (/4, 5, 1, 6, 2, 3/)
    qp(11, :)  = (/5, 1, 4, 3, 6, 2/)
    qp(12, :)  = (/5, 4, 1, 6, 3, 2/)
    qp(13, :)  = (/2, 4, 6, 1, 3, 5/)
    qp(14, :)  = (/2, 6, 4, 3, 1, 5/)
    qp(15, :)  = (/4, 2, 6, 1, 5, 3/)
    qp(16, :)  = (/4, 6, 2, 5, 1, 3/)
    qp(17, :)  = (/6, 2, 4, 3, 5, 1/)
    qp(18, :)  = (/6, 4, 2, 5, 3, 1/)
    qp(19, :)  = (/3, 5, 6, 1, 2, 4/)
    qp(20, :)  = (/3, 6, 5, 2, 1, 4/)
    qp(21, :)  = (/5, 3, 6, 1, 4, 2/)
    qp(22, :)  = (/5, 6, 3, 4, 1, 2/)
    qp(23, :)  = (/6, 3, 5, 2, 4, 1/)
    qp(24, :)  = (/6, 5, 3, 4, 2, 1/)
    ps         = 4e0_wp
    moltype    = parameter_moltype; pesform = parameter_pesform
    n          = parameter_n; m = parameter_m; eta = parameter_eta
    pnu        = parameter_pnu; epsnu = parameter_epsnu
    aa         = parameter_morse_widths; qeq = parameter_equilibrium_distances
    nma        = parameter_nma; nplmax = parameter_nplmax
    nmincrc    = parameter_nmincrc; mcrc = parameter_mcrc
    allocate(crc(mcrc)); crc = parameter_crc
    ntarmin = parameter_ntarmin; ntarmax = parameter_ntarmax; ntar = parameter_ntar
    crcmult = parameter_crcmult; nitcrc = parameter_nitcrc
    rwfa    = parameter_rwfa; rwfb = parameter_rwfb
    nall    = n+nma
    allocate(qi(n, 6), aji(n, m))
    do i=1, n
        do l=1, 6
            t(l) = geometry_value(i, l)
        enddo
        call default2cartesian(t(1), t(2), t(3), t(4), t(5), t(6), y)
        call cartesian2rawdist(y, qr)
        call rawdist2rawmorse(qr, qi(i, :), qeq, aa)
        do j=1, m
            aji(i, j) = coefficient_value(i, j)
        enddo
    enddo
end subroutine assemble_model_parameters

!> @brief Evaluate pot.
subroutine pot( x_, y_, z_, e_, dedx_, dedy_, dedz_ )
    use falcons_vars_eval, only: aji, qi
    use falcons_vars_main, only: aa, cgmeth, cmap, crc, crcmult, dedy, e, epsnu, epsu, epw, eta, i, &
        & ind, j, k, l, m, mcrc, moltype, n, nall, nitcrc, nma, nmincrc, nn, np, nplmax, ntar, &
        & ntarmax, ntarmin, pesform, pnu, ps, pu, pwonly, q, qeq, qp, qps, qr, qs, r, rwfa, rwfb, sd, &
        & systo, t, x, xr, y
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), dimension(4), intent(in)  :: x_, y_, z_
    real(wp), intent(out)               :: e_
    real(wp), dimension(4), intent(out) :: dedx_, dedy_, dedz_
    real(wp)                            :: cconv, econv, gconv, eref
    real(wp), dimension(4)              :: x__, y__, z__
    pu     = 3e0_wp
    epsu   = 1e-1_wp
    pwonly = 0
    cgmeth = 1
    cconv  = 0.52917721092e0_wp
    econv  = 0.159360144e-2_wp
    gconv  = 0.843297564e-3_wp
    eref   = -218.4080132e0_wp
    x__    = x_*cconv
    y__    = y_*cconv
    z__    = z_*cconv
    do l=1, 4
        y(l, 1) = x__(l)
        y(l, 2) = y__(l)
        y(l, 3) = z__(l)
    enddo
    call calccoordinates(y, qr, xr, q, x, sd, aa, qeq, moltype, qp, qps)
    call evaluatesurface(cgmeth, y, qr, xr, q, x, moltype, pwonly, &
        & cmap, qp, qps, mcrc, crc, n, m, eta, ps, pu, epsu, aa, r, np, e, dedy)
    e_    = e
    dedx_ = dedy(:, 1)
    dedy_ = dedy(:, 2)
    dedz_ = dedy(:, 3)
    e_    = e_*econv+ eref
    dedx_ = dedx_*gconv
    dedy_ = dedy_*gconv
    dedz_ = dedz_*gconv
end subroutine pot

!> @brief Evaluate evaluatesurface.
subroutine evaluatesurface( cg, y, qr, xr, q, x, moltype, pwonly, &
        & cmap, qp, qps, mcrc, crc, n, m, eta, ps, pu, epsu, aa, r, np, e, dedy )
    use falcons_vars_eval, only: aji, qi
    use chempotpy_kinds,   only: wp
    implicit none
    integer, intent(in)                    :: cg, moltype, pwonly, qps, mcrc, n, m, eta
    integer, dimension(4, 3), intent(in)   :: cmap
    integer, dimension(qps, 6), intent(in) :: qp
    real(wp), intent(in)                   :: ps, pu, epsu
    real(wp), dimension(6), intent(in)     :: qr, xr, q, x, aa
    real(wp), dimension(mcrc), intent(in)  :: crc
    real(wp), dimension(4, 3), intent(in)  :: y
    integer, intent(out)                   :: np
    real(wp), intent(out)                  :: r, e
    real(wp), dimension(4, 3), intent(out) :: dedy
    integer                                :: i, j, k, l, a, b, c
    real(wp)                               :: rs, drdqch, dis, wsum, wi, pi, dwidz, tmp, epw
    real(wp), dimension(6)                 :: qii, dwidxr, dwdxrsum, dpidx, dpidxr, &
        & dedxrpart, dedxr, dedqr, ddisdq, depwdqr
    real(wp), dimension(m)       :: aj, bj
    real(wp), dimension(6, 6)    :: dxdxr
    real(wp), dimension(m, 6)    :: dbjdx
    real(wp), dimension(6, 4, 3) :: dqrdy
    np   = 0
    e    = 0e0_wp
    wsum = 0e0_wp !sum of the weight functions
    do a=1, 6
        dedxrpart(a) = 0e0_wp
        dwdxrsum(a)  = 0e0_wp !sums of the weight function derivatives
    enddo
    call calcqch(q, tmp)
    call calccr(tmp, mcrc, crc, r, rs, drdqch)
    call calcbj(bj, dbjdx, x, eta, m, 1)
    if( cg /= 0 )then
        call calcinvariantsjacobian(xr, x, dxdxr, moltype)
        call calcdqrdyjacobian(qr, y, dqrdy)
    endif
    call calcsumpairwise(qr, epw, depwdqr, moltype, cg)
    if( pwonly == 0 )then
        do i=1, n
            do a=1, 6
                qii(a) = qi(i, a)
            enddo
            call distance6ds_perm(q, qii, qp, qps, 0, dis, ddisdq)
            if( dis <= rs )then
                np = np+1
                if( cg /= 0 )then
                    call distance6ds_perm(q, qii, qp, qps, cg, dis, ddisdq)
                endif
                do j=1, m
                    aj(j) = aji(i, j)
                enddo
                call weightw(wi, dwidz, dis, rs, ps, pu, epsu)
                pi = 0e0_wp
                do j=1, m
                    pi = pi+ aj(j)*bj(j)
                enddo
                wsum = wsum+ wi
                e    = e+ wi*pi
                if( cg /= 0 )then
                    do a=1, 6
                        dpidx(a) = 0e0_wp
                    enddo
                    do j=1, m
                        do a=1, 6
                            dpidx(a) = dpidx(a)+ aj(j)*dbjdx(j, a)
                        enddo
                    enddo
                    do a=1, 6
                        tmp = 0e0_wp
                        do k=1, 6
                            tmp = tmp+ dpidx(k)*dxdxr(k, a)
                        enddo
                        dpidxr(a) = tmp
                    enddo
                    tmp = dis*4e0_wp*drdqch/r
                    do a=1, 6
                        dwidxr(a)   = dwidz*(ddisdq(a)- tmp*q(a))/rs
                        dwdxrsum(a) = dwdxrsum(a)+ dwidxr(a)
                    enddo
                    do a=1, 6
                        dedxrpart(a) = dedxrpart(a)+ wi*dpidxr(a)+ pi*dwidxr(a)
                    enddo
                endif
            endif
        enddo  !end loop over data points within the cutoff radius
    endif  !end if 'pwonly==0' computation
    if( wsum == 0e0_wp )then
        wsum = 1e0_wp
    endif
    e = e/wsum
    if( cg /= 0 )then
        do a=1, 6
            dedxr(a) = (dedxrpart(a)- e*dwdxrsum(a))/wsum
        enddo
        do a=1, 6
            dedqr(a) = -dedxr(a)*xr(a)/aa(a)+ depwdqr(a)
        enddo
        do b=1, 4
            do c=1, 3
                tmp = 0e0_wp
                do k=1, 3
                    a   = cmap(b, k)
                    tmp = tmp+ dedqr(a)*dqrdy(a, b, c)
                enddo
                dedy(b, c) = tmp
            enddo
        enddo
    else
        do b=1, 4
            do c=1, 3
                dedy(b, c) = 0e0_wp
            enddo
        enddo
    endif
    e = e+ epw
    if( np == 0 )then
        write(*, '(a)') &
            &  ' WARNING: No points found within cutoff radius for'//&
            &  ' evaluation point:'
        write(*, '(a,6(1x,f15.5))') ' qr= ', (qr(l), l=1, 6)
        write(*, '(a,6(1x,f15.5))') ' q = ', (q(l), l=1, 6)
        write(*, '(a,6(1x,f15.5))') ' x = ', (x(l), l=1, 6)
    endif
end subroutine evaluatesurface

!> @brief Evaluate calcqch.
subroutine calcqch( q, qch )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), dimension(6), intent(in) :: q
    real(wp), intent(out)              :: qch
    integer                            :: i
    real(wp)                           :: val
    val = 0e0_wp
    do i=1, 6
        val = val+ q(i)**2e0_wp
    enddo
    qch = val
end subroutine calcqch

!> @brief Evaluate calccr.
subroutine calccr( qch, mcrc, crc, cr, crs, dcrdqch )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)                   :: mcrc
    real(wp), intent(in)                  :: qch
    real(wp), intent(out)                 :: cr, crs, dcrdqch
    real(wp), dimension(mcrc), intent(in) :: crc
    integer                               :: icrc
    real(wp)                              :: pqch
    pqch    = 1e0_wp
    cr      = crc(1)
    dcrdqch = 0e0_wp
    do icrc= 2, mcrc
        dcrdqch = dcrdqch+ (crc(icrc)*real(icrc-1, wp)*pqch)
        pqch    = pqch*qch
        cr      = cr+ (crc(icrc)* pqch)
    enddo
    crs = cr*cr
end subroutine calccr

!> @brief Evaluate calcbj.
subroutine calcbj( b, dbdx, x, eta, m, dflag )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)                    :: eta, m, dflag
    real(wp), dimension(6), intent(in)     :: x
    real(wp), dimension(m), intent(out)    :: b
    real(wp), dimension(m, 6), intent(out) :: dbdx
    integer                                :: ind, i, j, k, l, ll, a
    real(wp)                               :: tmp, tmp2, tmp3
    real(wp), dimension(6)                 :: xx, x2, xxx
    real(wp), dimension(m)                 :: dbdxa
    real(wp), dimension(6, 6)              :: xy, xy2
    real(wp), dimension(6, 6, 6)           :: xyz
    b    = 0e0_wp
    dbdx = 0e0_wp
    xx   = 0e0_wp !xx(j)     = x(j)*x(j)
    x2   = 0e0_wp !x2(j)     = x(j)*2d0
    xxx  = 0e0_wp !xxx(j)    = x(j)*x(j)*x(j)
    xy   = 0e0_wp !xy(j,k)   = x(j)*x(k)
    xy2  = 0e0_wp !xy2(j)    = x(j)*x(k)*2d0
    xyz  = 0e0_wp !xyz(j,k,l)= x(j)*x(k)*x(l)
    tmp  = 0e0_wp
    tmp2 = 0e0_wp
    tmp3 = 0e0_wp
    ind  = 0
    do j=1, 6
        ind    = ind+1
        b(ind) = x(j) !linear terms
    enddo
    if( eta .ge. 2 )then
        do j= 1, 6
            ind    = ind+1
            tmp    = x(j)*x(j)
            xx(j)  = tmp
            b(ind) = tmp !quadratic terms such as x1^2
            do k= j+1, 6
                ind      = ind+1
                tmp      = x(j)*x(k)
                xy(j, k) = tmp
                b(ind)   = tmp
            enddo
        enddo
    endif
    if( eta .ge. 3 )then
        do j= 1, 6
            ind    = ind+1
            tmp    = xx(j)
            tmp2   = tmp*x(j)
            b(ind) = tmp2 !cubic terms such as x1^3
            xxx(j) = tmp2
            do k= j+1, 6
                ind    = ind+1
                b(ind) = tmp*x(k)
                ind    = ind+1
                b(ind) = x(j)*xx(k)
            enddo
        enddo
        do j= 1, 6
            do k= j+1, 6
                tmp = xy(j, k)
                do l= k+1, 6
                    ind          = ind+1
                    tmp2         = tmp*x(l)
                    b(ind)       = tmp2
                    xyz(j, k, l) = tmp2
                enddo
            enddo
        enddo
    endif
    if( eta .ge. 4 )then
        do j= 1, 6
            ind    = ind+1
            tmp    = xxx(j)
            b(ind) = tmp*x(j) !quartic terms such as x1^4
            do k= j+1, 6
                ind    = ind+1
                b(ind) = tmp*x(k)
                ind    = ind+1
                b(ind) = x(j)*xxx(k)
            enddo
        enddo
        do j= 1, 6
            tmp = xx(j)
            do k= j+1, 6
                tmp2   = xx(k)
                ind    = ind+1
                b(ind) = tmp*tmp2
                tmp3   = xy(j, k)
                do l= k+1, 6
                    ind    = ind+1
                    b(ind) = tmp*xy(k, l)
                    ind    = ind+1
                    b(ind) = tmp2*xy(j, l)
                    ind    = ind+1
                    b(ind) = xx(l)*tmp3
                enddo
            enddo
        enddo
        do j=1, 6
            do k= j+1, 6
                tmp = xy(j, k)
                do l= k+1, 6
                    do ll= l+1, 6
                        ind    = ind+1
                        b(ind) = tmp*xy(l, ll)
                    enddo
                enddo
            enddo
        enddo
    endif
    if( dflag /= 0 )then
        do j=1, 6
            x2(j) = 2e0_wp*x(j)
            do k=j+1, 6
                xy2(j, k) = 2e0_wp*xy(j, k)
            enddo
        enddo
        do a=1, 6
            do j=1, m
                dbdxa(j) = 0e0_wp
            enddo
            dbdxa(a) = 1e0_wp !derivative of linear term
            ind      = 6
            if( eta .ge. 2 )then
                do j= 1, 6
                    ind = ind+1
                    if( j == a )then
                        dbdxa(ind) = x2(j)
                    endif
                    do k= j+1, 6
                        ind = ind+1
                        if( j == a )then
                            dbdxa(ind) = x(k)
                        elseif( k == a )then
                            dbdxa(ind) = x(j)
                        endif
                    enddo
                enddo
            endif
            if( eta .ge. 3 )then
                do j= 1, 6
                    ind = ind+1
                    if( j == a )then
                        dbdxa(ind) = 3e0_wp*xx(a)
                        do k= j+1, 6
                            ind        = ind+1
                            dbdxa(ind) = xy2(j, k)
                            ind        = ind+1
                            dbdxa(ind) = xx(k)
                        enddo
                    else
                        do k= j+1, 6
                            if( k == a )then
                                ind        = ind+1
                                dbdxa(ind) = xx(j)
                                ind        = ind+1
                                dbdxa(ind) = xy2(j, k)
                            else
                                ind = ind+2
                            endif
                        enddo
                    endif
                enddo
                do j= 1, 6
                    do k= j+1, 6
                        do l= k+1, 6
                            ind = ind+1
                            if( j == a )then
                                dbdxa(ind) = xy(k, l)
                            elseif( k == a )then
                                dbdxa(ind) = xy(j, l)
                            elseif( l == a )then
                                dbdxa(ind) = xy(j, k)
                            endif
                        enddo
                    enddo
                enddo
            endif
            if( eta .ge. 4 )then
                do j= 1, 6
                    ind = ind+1
                    if( j == a )then
                        dbdxa(ind) = 4e0_wp*xxx(a)
                        do k= j+1, 6
                            ind        = ind+1
                            dbdxa(ind) = 3e0_wp*xx(j)*x(k)
                            ind        = ind+1
                            dbdxa(ind) = xxx(k)
                        enddo
                    else
                        do k= j+1, 6
                            if( k == a )then
                                ind        = ind+1
                                dbdxa(ind) = xxx(j)
                                ind        = ind+1
                                dbdxa(ind) = 3e0_wp*x(j)*xx(k)
                            else
                                ind = ind+2
                            endif
                        enddo
                    endif
                enddo
                do j= 1, 6
                    tmp = xx(j)
                    do k= j+1, 6
                        tmp2 = xx(k)
                        ind  = ind+1
                        if( j == a )then
                            dbdxa(ind) = x2(j)*tmp2
                        elseif( k == a )then
                            dbdxa(ind) = tmp*x2(k)
                        endif
                        do l= k+1, 6
                            ind = ind+1
                            if( j == a )then
                                dbdxa(ind) = x2(j)*xy(k, l)
                            elseif( k == a )then
                                dbdxa(ind) = tmp*x(l)
                            elseif( l == a )then
                                dbdxa(ind) = tmp*x(k)
                            endif
                            ind = ind+1
                            if( j == a )then
                                dbdxa(ind) = tmp2*x(l)
                            elseif( k == a )then
                                dbdxa(ind) = x2(k)*xy(j, l)
                            elseif( l == a )then
                                dbdxa(ind) = tmp2*x(j)
                            endif
                            ind = ind+1
                            if( j == a )then
                                dbdxa(ind) = xx(l)*x(k)
                            elseif( k == a )then
                                dbdxa(ind) = xx(l)*x(j)
                            elseif( l == a )then
                                dbdxa(ind) = x2(l)*xy(j, k)
                            endif
                        enddo
                    enddo
                enddo
                do j= 1, 6
                    do k= j+1, 6
                        do l= k+1, 6
                            do ll= l+1, 6
                                ind = ind+1
                                if( j == a )then
                                    dbdxa(ind) = xyz(k, l, ll)
                                elseif( k == a )then
                                    dbdxa(ind) = xyz(j, l, ll)
                                elseif( l == a )then
                                    dbdxa(ind) = xyz(j, k, ll)
                                elseif( ll == a )then
                                    dbdxa(ind) = xyz(j, k, l)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            endif
            do j=1, m
                dbdx(j, a) = dbdxa(j)
            enddo
        enddo !end 'do a=1,6'
    endif !end 'if (dflag /= 0) then'
end subroutine calcbj

!> @brief Evaluate distance6ds_perm.
subroutine distance6ds_perm( q, qi, qp, qps, cg, ds, ddsdq )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)                    :: cg, qps
    integer, dimension(qps, 6), intent(in) :: qp
    real(wp), intent(out)                  :: ds
    real(wp), dimension(6), intent(in)     :: q, qi
    real(wp), dimension(6), intent(out)    :: ddsdq
    integer                                :: i, a, iszero
    integer, dimension(6)                  :: sd
    real(wp)                               :: dsval, dssum, mult, dsvali, dsvalis, dssumi
    real(wp), dimension(6)                 :: qdif
    iszero = 0
    dssum  = 0e0_wp
    do a=1, 6
        ddsdq(a) = 0e0_wp
    enddo
    do i=1, qps
        do a=1, 6
            sd(a) = qp(i, a)
        enddo
        do a=1, 6
            qdif(a) = q(a)- qi(sd(a))
        enddo
        dsval = 0e0_wp
        do a=1, 6
            dsval = dsval+ qdif(a)**2e0_wp
        enddo
        if( dsval > 0e0_wp )then
            dsvali  = 1e0_wp/dsval
            dsvalis = dsvali*dsvali !dsvalis= dsval^(-2)
            dssum   = dssum+ dsvalis
            if( cg /= 0 )then
                mult = 2e0_wp*dsvali*dsvalis !mult= 2*dsval^(-3)
                do a=1, 6
                    ddsdq(a) = ddsdq(a)+ mult*qdif(a)
                enddo
            endif
        else
            iszero = 1
        endif
    enddo
    if( iszero /= 1 )then
        dssumi = 1e0_wp/dssum
        ds     = sqrt(dssumi) !ds= dssum^(-0.5)
        if( cg /= 0 )then
            mult = dssumi*ds  !mult= dssum^(-1.5)
            do a=1, 6
                ddsdq(a) = mult*ddsdq(a)
            enddo
        endif
    else
        ds = 0e0_wp
        if( cg /= 0 )then
            do a=1, 6
                ddsdq(a) = 0e0_wp
            enddo
        endif
    endif
end subroutine distance6ds_perm

!> @brief Evaluate weightw.
subroutine weightw( w, dwdz, ds, rs, ps, pu, epsu )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), intent(in)  :: ds, ps, pu, epsu, rs
    real(wp), intent(out) :: w, dwdz
    real(wp)              :: z, uval, sval
    z    = ds/rs
    uval = 1e0_wp/((z**pu)+ (epsu**pu))
    sval = (1e0_wp- z**ps)**4e0_wp
    w    = uval*sval
    dwdz = uval*4e0_wp*((1e0_wp-z**ps)**3e0_wp)*(-ps*(z**(ps-1e0_wp)))+ sval*(uval*uval)*(-(pu*(z**(pu-1e0_wp))))
end subroutine weightw

!> @brief Evaluate calccoordinates.
subroutine calccoordinates( &
        & y, qr, xr, q, x, sd, aa, qeq, moltype, qp, qps )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)                 :: moltype, qps
    integer, dimension(6), intent(out)  :: sd
    real(wp), dimension(6), intent(in)  :: aa, qeq
    real(wp), dimension(6), intent(out) :: qr, xr, q, x
    integer, dimension(qps, 6)          :: qp
    real(wp), dimension(4, 3)           :: y
    call cartesian2rawdist(y, qr)
    call rawdist2rawmorse(qr, xr, qeq, aa)
    q = xr
    call calcinvariants(xr, x, moltype)
end subroutine calccoordinates

!> @brief Evaluate default2cartesian.
subroutine default2cartesian( ra, rb, d, thetaa_d, thetab_d, &
        &  phi_d, y )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), intent(in)                   :: ra, rb, d, thetaa_d, thetab_d, phi_d
    real(wp), dimension(4, 3), intent(out) :: y
    real(wp)                               :: thetaa, thetab, phi, pi, y11, y31, y32
    pi      = 3.1415926536e0_wp
    thetaa  = thetaa_d*pi/180e0_wp
    thetab  = thetab_d*pi/180e0_wp
    phi     = phi_d*pi/180e0_wp
    y11     = ra*sin(thetaa)/2e0_wp
    y(1, 1) = y11
    y(1, 2) = 0e0_wp
    y(1, 3) = (d- ra*cos(thetaa))/2e0_wp
    y(2, 1) = -y11
    y(2, 2) = 0e0_wp
    y(2, 3) = (d+ ra*cos(thetaa))/2e0_wp
    y31     = rb*sin(thetab)*cos(phi)/2e0_wp
    y(3, 1) = y31
    y32     = rb*sin(thetab)*sin(phi)/2e0_wp
    y(3, 2) = y32
    y(3, 3) = (-d+ rb*cos(thetab))/2e0_wp
    y(4, 1) = -y31
    y(4, 2) = -y32
    y(4, 3) = -(d+ rb*cos(thetab))/2e0_wp
end subroutine default2cartesian

!> @brief Evaluate cartesian2rawdist.
subroutine cartesian2rawdist( y, qr )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), dimension(6), intent(out)   :: qr
    real(wp), dimension(4, 3), intent(in) :: y
    integer                               :: i
    call distance3d(y(1, :), y(2, :), qr(1))
    call distance3d(y(1, :), y(3, :), qr(2))
    call distance3d(y(1, :), y(4, :), qr(3))
    call distance3d(y(2, :), y(3, :), qr(4))
    call distance3d(y(2, :), y(4, :), qr(5))
    call distance3d(y(3, :), y(4, :), qr(6))
end subroutine cartesian2rawdist

!> @brief Evaluate distance3d.
subroutine distance3d( p1, p2, d )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), intent(out)              :: d
    real(wp), dimension(3), intent(in) :: p1, p2
    d = sqrt((p1(1)-p2(1))**2e0_wp+(p1(2)-p2(2))**2e0_wp+(p1(3)-p2(3))**2e0_wp)
end subroutine distance3d

!> @brief Evaluate rawdist2rawmorse.
subroutine rawdist2rawmorse( qr, xr, qeq, aa )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), dimension(6), intent(in)  :: qr, qeq, aa
    real(wp), dimension(6), intent(out) :: xr
    integer                             :: i
    do i=1, 6
        xr(i) = exp(-(qr(i)-qeq(i))/aa(i))
    enddo
end subroutine rawdist2rawmorse

!> @brief Evaluate calcinvariants.
subroutine calcinvariants( vr, v, moltype )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)                 :: moltype
    real(wp), dimension(6), intent(in)  :: vr
    real(wp), dimension(6), intent(out) :: v
    integer                             :: i
    real(wp)                            :: e0, e1, f0, f1, vrtmp
    real(wp), dimension(4)              :: vv
    v = 0e0_wp
    select case (moltype)
    case (1)
        v(1)= (vr(1)*vr(2)+vr(1)*vr(3)+vr(1)*vr(4)+ vr(1)*vr(5)+vr(2)*vr(3)+vr(2)*vr(4)+ vr(2)*vr(6)+vr(3)*vr(5)+ &
            & vr(3)* vr(6)+ vr(4)*vr(5)+vr(4)*vr(6)+vr(5)*vr(6))/12e0_wp
        v(1) = v(1)       !LEGACY: **0.5d0
        v(2) = (vr(1)*vr(2)*vr(4)+ vr(1)*vr(3)*vr(5)+ vr(2)*vr(3)*vr(6)+ vr(4)*vr(5)*vr(6))/4e0_wp
        v(2) = v(2)       !LEGACY: **(1d0/3d0)
        v(3)= (vr(1)*vr(2)*(vr(1)+vr(2))+ vr(1)*vr(3)*(vr(1)+vr(3))+ vr(1)*vr(4)*(vr(1)+vr(4))+ vr(1)*vr(5)*(vr(1)+ &
            & vr(5))+ vr(2)*vr(3)*(vr(2)+vr(3))+ vr(2)*vr(4)*(vr(2)+vr(4))+ vr(2)*vr(6)*(vr(2)+vr(6))+ vr(3)*vr(5)* &
            & (vr(3)+vr(5))+ vr(3)*vr(6)*(vr(3)+vr(6))+ vr(4)*vr(5)*(vr(4)+vr(5))+ vr(4)*vr(6)*(vr(4)+vr(6))+ vr(5)* &
            & vr(6)*(vr(5)+vr(6)))/24e0_wp
        v(3) = v(3)       !LEGACY: **(1d0/3d0)
        v(4) = (vr(1)*vr(2)*vr(3)+ vr(1)*vr(4)*vr(5)+ vr(2)*vr(4)*vr(6)+ vr(3)*vr(5)*vr(6))/4e0_wp
        v(4) = v(4)       !LEGACY: **(1d0/3d0)
        v(5)= (vr(1)*vr(3)*vr(4)+ vr(2)*vr(3)*vr(4)+ vr(1)*vr(2)*vr(5)+ vr(2)*vr(3)*vr(5)+ vr(2)*vr(4)*vr(5)+ vr(3)* &
            & vr(4)*vr(5)+ vr(1)*vr(2)*vr(6)+ vr(1)*vr(3)*vr(6)+ vr(1)*vr(4)*vr(6)+ vr(3)*vr(4)*vr(6)+ vr(1)*vr(5)* &
            & vr(6)+ vr(2)*vr(5)*vr(6))/12e0_wp
        v(5) = v(5)       !LEGACY: **(1d0/3d0)
        v(6) = (vr(1)*vr(2)*vr(5)*vr(6)+ vr(1)*vr(3)*vr(4)*vr(6)+ vr(2)*vr(3)*vr(4)*vr(5))/3e0_wp
        v(6) = v(6)       !LEGACY: **0.25d0
    case (2)
        v(1) = (vr(2)*vr(4)+ vr(3)*vr(5))/2e0_wp
        v(2) = (vr(2)*vr(3)+ vr(4)*vr(5))/2e0_wp
        v(3) = (vr(1)*(vr(2)+vr(3)+vr(4)+vr(5)))/4e0_wp
        v(4) = (vr(6)*(vr(2)+vr(3)+vr(4)+vr(5)))/4e0_wp
        v(5) = (vr(1)*(vr(2)*vr(3)+vr(4)*vr(5)))/2e0_wp
        v(6) = (vr(6)*(vr(2)*vr(4)+vr(3)*vr(5)))/2e0_wp
    end select
end subroutine calcinvariants

!> @brief Evaluate calcinvariantsjacobian.
subroutine calcinvariantsjacobian( vr, v, dvdvr, moltype )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)                    :: moltype
    real(wp), dimension(6), intent(in)     :: vr, v
    real(wp), dimension(6, 6), intent(out) :: dvdvr
    integer                                :: i
    real(wp)                               :: mult
    real(wp), dimension(6)                 :: tmp, tmp2
    select case (moltype)
    case (1)
        mult   = 1e0_wp/4e0_wp     !LEGACY: *3d0*(v(4)**2d0))
        tmp(1) = vr(2)*vr(3)+ vr(4)*vr(5)
        tmp(2) = vr(1)*vr(3)+ vr(4)*vr(6)
        tmp(3) = vr(1)*vr(2)+ vr(5)*vr(6)
        tmp(4) = vr(1)*vr(5)+ vr(2)*vr(6)
        tmp(5) = vr(1)*vr(4)+ vr(3)*vr(6)
        tmp(6) = vr(2)*vr(4)+ vr(3)*vr(5)
        do i=1, 6
            dvdvr(4, i) = tmp(i)*mult
        enddo
        mult = 1e0_wp/4e0_wp     !LEGACY: *3d0*(v(2)**2d0))
        do i=1, 6
            dvdvr(2, i) = tmp(7-i)*mult
        enddo
        mult   = 1e0_wp/12e0_wp    !LEGACY: *2d0*v(1))
        tmp(1) = vr(2)+vr(3)+vr(4)+vr(5)
        tmp(2) = vr(1)+vr(3)+vr(4)+vr(6)
        tmp(3) = vr(1)+vr(2)+vr(5)+vr(6)
        tmp(4) = tmp(3)
        tmp(5) = tmp(2)
        tmp(6) = tmp(1)
        do i=1, 6
            dvdvr(1, i) = tmp(i)*mult
        enddo
        mult    = 1e0_wp/12e0_wp    !LEGACY: *3d0*(v(5)**2d0))
        tmp2(1) = vr(3)*vr(4)+ vr(2)*vr(5)
        tmp2(2) = vr(3)*vr(4)+ vr(1)*vr(6)
        tmp2(3) = vr(2)*vr(5)+ vr(1)*vr(6)
        tmp2(4) = tmp2(3)
        tmp2(5) = tmp2(2)
        tmp2(6) = tmp2(1)
        do i=1, 6
            dvdvr(5, i) = (tmp(i)*vr(7-i)+ tmp2(i))*mult
        enddo
        mult = 1e0_wp/3e0_wp     !LEGACY: *4d0*(v(6)**3d0))
        do i=1, 6
            dvdvr(6, i) = vr(7-i)*tmp2(i)*mult
        enddo
        mult    = 1e0_wp/24e0_wp    !LEGACY: *3d0*(v(3)**2d0))
        tmp2(1) = vr(2)**2e0_wp+vr(3)**2e0_wp+vr(4)**2e0_wp+vr(5)**2e0_wp
        tmp2(2) = vr(1)**2e0_wp+vr(3)**2e0_wp+vr(4)**2e0_wp+vr(6)**2e0_wp
        tmp2(3) = vr(1)**2e0_wp+vr(2)**2e0_wp+vr(5)**2e0_wp+vr(6)**2e0_wp
        tmp2(4) = tmp2(3)
        tmp2(5) = tmp2(2)
        tmp2(6) = tmp2(1)
        do i=1, 6
            dvdvr(3, i) = (tmp2(i)+ 2e0_wp*vr(i)*tmp(i))*mult
        enddo
    case (2)
        mult   = 1e0_wp/2e0_wp
        tmp(1) = 0e0_wp
        tmp(2) = vr(4)
        tmp(3) = vr(5)
        tmp(4) = vr(2)
        tmp(5) = vr(3)
        tmp(6) = 0e0_wp
        do i=1, 6
            dvdvr(1, i) = tmp(i)*mult
        enddo
        tmp(1) = 0e0_wp
        tmp(2) = vr(3)
        tmp(3) = vr(2)
        tmp(4) = vr(5)
        tmp(5) = vr(4)
        tmp(6) = 0e0_wp
        do i=1, 6
            dvdvr(2, i) = tmp(i)*mult
        enddo
        tmp(1) = 2e0_wp*v(2)
        tmp(2) = vr(1)*vr(3)
        tmp(3) = vr(1)*vr(2)
        tmp(4) = vr(1)*vr(5)
        tmp(5) = vr(1)*vr(4)
        tmp(6) = 0e0_wp
        do i=1, 6
            dvdvr(5, i) = tmp(i)*mult
        enddo
        tmp(1) = 0e0_wp
        tmp(2) = vr(6)*vr(4)
        tmp(3) = vr(6)*vr(5)
        tmp(4) = vr(6)*vr(2)
        tmp(5) = vr(6)*vr(3)
        tmp(6) = 2e0_wp*v(1)
        do i=1, 6
            dvdvr(6, i) = tmp(i)*mult
        enddo
        mult   = 1e0_wp/4e0_wp
        tmp(1) = vr(2)+vr(3)+vr(4)+vr(5)
        tmp(2) = vr(1)
        tmp(3) = vr(1)
        tmp(4) = vr(1)
        tmp(5) = vr(1)
        tmp(6) = 0e0_wp
        do i=1, 6
            dvdvr(3, i) = tmp(i)*mult
        enddo
        tmp(1) = 0e0_wp
        tmp(2) = vr(6)
        tmp(3) = vr(6)
        tmp(4) = vr(6)
        tmp(5) = vr(6)
        tmp(6) = vr(2)+vr(3)+vr(4)+vr(5)
        do i=1, 6
            dvdvr(4, i) = tmp(i)*mult
        enddo
    end select
end subroutine calcinvariantsjacobian

!> @brief Evaluate calcdqrdyjacobian.
subroutine calcdqrdyjacobian( qr, y, dqrdy )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), dimension(6), intent(in)        :: qr
    real(wp), dimension(4, 3), intent(in)     :: y
    real(wp), dimension(6, 4, 3), intent(out) :: dqrdy
    integer                                   :: i
    real(wp)                                  :: val
    dqrdy = 0e0_wp
    do i=1, 3
        val            = y(1, i)-y(2, i)
        dqrdy(1, 1, i) = val
        dqrdy(1, 2, i) = -val
        val            = y(1, i)-y(3, i)
        dqrdy(2, 1, i) = val
        dqrdy(2, 3, i) = -val
        val            = y(1, i)-y(4, i)
        dqrdy(3, 1, i) = val
        dqrdy(3, 4, i) = -val
        val            = y(2, i)-y(3, i)
        dqrdy(4, 2, i) = val
        dqrdy(4, 3, i) = -val
        val            = y(2, i)-y(4, i)
        dqrdy(5, 2, i) = val
        dqrdy(5, 4, i) = -val
        val            = y(3, i)-y(4, i)
        dqrdy(6, 3, i) = val
        dqrdy(6, 4, i) = -val
    enddo
    do i=1, 6
        dqrdy(i, :, :) = dqrdy(i, :, :)/qr(i)
    enddo
end subroutine calcdqrdyjacobian

!> @brief Evaluate calcpairwise.
subroutine calcpairwise( r, v2, dv2dr, imol, igrad )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in) :: igrad, imol
    real(wp), dimension(0:10) :: c
    real(wp), intent(in)  :: r
    real(wp), intent(out) :: v2, dv2dr
    real(wp)              :: re, de, fy, dfdy, dfdr, u, y, dydr, r2, r4, re2, re4, y2, y4, v
    select case (imol)
    case (1)
        re   = 1.098e0_wp
        de   = 228.7e0_wp
        c(0) = 2.70963254293e0_wp
        c(1) = 1.32620177271e-1_wp
        c(2) = 2.96757048793e-1_wp
        c(3) = 1.97112432229e-1_wp
        c(4) = -5.02002309588e-1_wp
        c(5) = 3.80734244606e-1_wp
        c(6) = 1.21001628750e0_wp
    case (2)
        re   = 1.098e0_wp
        de   = 228.4e0_wp
        c(0) = 2.71405774451e0_wp
        c(1) = 1.32757649829e-1_wp
        c(2) = 2.66756890408e-1_wp
        c(3) = 1.95350725241e-1_wp
        c(4) = -4.08663480982e-1_wp
        c(5) = 3.92451705557e-1_wp
        c(6) = 1.13006674877e0_wp
    case (3)
        re   = 1.1508e0_wp
        de   = 149.9e0_wp
        c(0) = 2.81134495569e0_wp
        c(1) = 1.43241169611e-1_wp
        c(2) = 1.35760038863e-2_wp
        c(3) = 3.92892178507e-1_wp
        c(4) = 9.29495534058e-1_wp
        c(5) = 2.66966672332e-1_wp
        c(6) = -3.68118714223e-1_wp
    case (4)
        re   = 1.208e0_wp
        de   = 120.243e0_wp
        c(0) = 2.69132890094e0_wp
        c(1) = 3.39550045614e-1_wp
        c(2) = 3.46401777195e-1_wp
        c(3) = -7.76983671636e-1_wp
        c(4) = -3.29632972405e-1_wp
        c(5) = 2.40883331247e0_wp
        c(6) = 2.09264029009e0_wp
    end select
    r2  = r*r
    r4  = r2*r2
    re2 = re*re
    re4 = re2*re2
    v   = r4+ re4
    y   = (r4 - re4)/v
    y2  = y*y
    y4  = y2*y2
    fy  = c(0) + c(1)*y + c(2)*y2 + c(3)*y2*y + c(4)*y4 + c(5)*y4*y + c(6)*y4*y2
    u   = exp(-fy*(r-re))
    v2  = de*(1.0e0_wp-u)*(1.0e0_wp-u)- de
    if( igrad .eq. 1 )then
        dfdy  = c(1)+ 2e0_wp*c(2)*y+ 3e0_wp*c(3)*y2+ 4e0_wp*c(4)*y2*y+ 5e0_wp*c(5)*y4 + 6e0_wp*c(6)*y4*y
        dydr  = 8.0e0_wp*r2*r*re4/(v*v)
        dfdr  = dfdy*dydr
        dv2dr = 2.0e0_wp*de*(1e0_wp-u)*u*(dfdr*(r-re)+fy)
    endif
end subroutine calcpairwise

!> @brief Evaluate calcsumpairwise.
subroutine calcsumpairwise( rvec, v2, dv2dr, systype, igrad )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)                 :: systype, igrad
    real(wp), dimension(6), intent(in)  :: rvec
    real(wp), dimension(6), intent(out) :: dv2dr
    integer                             :: a
    real(wp)                            :: v2, tmp, de
    v2 = 0e0_wp
    do a=1, 6
        dv2dr(a) = 0e0_wp
    enddo
    select case (systype)
    case (1)
        do a=1, 6
            call calcpairwise(rvec(a), tmp, dv2dr(a), 1, igrad)
            v2 = v2+ tmp
        enddo
        de = 228.7e0_wp
        v2 = v2+ 2e0_wp*de
    case (2)
        call calcpairwise(rvec(1), tmp, dv2dr(1), 4, igrad)
        v2 = v2+ tmp
        do a=2, 5
            call calcpairwise(rvec(a), tmp, dv2dr(a), 3, igrad)
            v2 = v2+ tmp
        enddo
        call calcpairwise(rvec(6), tmp, dv2dr(6), 2, igrad)
        v2 = v2+ tmp
        v2 = v2+ 228.4e0_wp+ 120.243e0_wp
    end select
end subroutine calcsumpairwise
