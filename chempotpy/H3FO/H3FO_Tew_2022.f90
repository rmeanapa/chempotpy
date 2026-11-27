! System:                        H3FO
! Functional form:               Morse, Coulomb, repulsion, and distributed Gaussian functions
! Common name:                   Tew H2O-HF surface
! Reference:                     D. Viglaska, X.-G. Wang, T. Carrington Jr., and D. P. Tew, J. Mol. Spectrosc. 384, 111587 (2022), doi:10.1016/j.jms.2022.111587
! Number of derivatives:         1
! Number of bodies:              5
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable released fit data are embedded
! Notes:                         Public atom order is O, H, H, H(HF), F; energy is relative to equilibrium.

!> @brief Implement the full-dimensional Tew 2022 H2O-HF potential with analytic gradients.
module chempotpy_h3fo_tew_2022_pes
    use chempotpy_h3fo_tew_parameters, only: coul, gauss, icoul, igauss, imorse, ixrep, morse, ncoul, ndim, &
        & ngauss, nmorse, nsym, nxrep, shift, xrep
    use chempotpy_kinds,  only: wp
    use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable
    use chempotpy_units,  only: bohr_to_angstrom, hartree_per_bohr_to_ev_per_angstrom, hartree_to_electronvolt
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate H3FO_Tew_2022 through the CHEMPOTPY Cartesian interface.
    subroutine evaluate( x, igrad, p, g, d, status )
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms = 5, nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
        real(wp)              :: native_coordinates(3,natoms), native_energy
        real(wp)              :: native_gradient(3*natoms), native_hessian(3*natoms*3*natoms)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad < 0 .or. igrad > 1 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'H3FO_Tew_2022 provides energies and first derivatives')
        else
            native_coordinates = transpose(x)/bohr_to_angstrom
            call evaluate_native(native_coordinates, igrad, native_energy, native_gradient, native_hessian)
            p(1) = native_energy*hartree_to_electronvolt
            if( igrad == 1 ) g(1,:,:) = transpose(reshape(native_gradient, [3, natoms]))* &
                hartree_per_bohr_to_ev_per_angstrom
        endif
    end subroutine evaluate

    ! Generated from the released Tew source, SHA-256 43eabf74c99d225aa902b519689d3f61c0653fbd4649f6725ad7a3954ae0334d.
    !> @brief Evaluate evaluate_native for the H3FO_Tew_2022 analytic potential.
    subroutine evaluate_native( x, iopt, e, g, h )
        implicit none
        integer, parameter :: natoms = 5, ndist=10  ! natoms*(natoms-1)/2
        ! Input:  Cartesian coordiates in bohr - connectivity assumed to be as eq structure
        real(wp), intent(in) :: x(3,natoms)
        ! Input:  iopt = 0,1,2
        integer, intent(in) :: iopt
        ! Output: energy above equilibrium in Hartree
        real(wp), intent(out) :: e, g(3*natoms), h(3*natoms*3*natoms)
        ! Local
        integer  :: i, j, k, l, ii, jj, ij, kl, ix, jx, isym, nd
        real(wp) :: r(4), dist(ndist), rij, v, gg(4), hh(10), rrij, d0
        real(wp) :: gint(ndist), hint((ndist+1)*ndist/2)
        real(wp) :: b(ndist,natoms*3), db(ndist,3*natoms,3*natoms)
        !          Transforms xyz, grad and hess from internal coordinates     !
        !           q = B x      B+ g' = g    B+ h' B + dB+ g' = h             !
        !  dV/dx = dV/dq dq/dx    d2V/dx2 = d2V/dq2 (dq/dx)2 + dV/dq d2q/dx2   !
        ij = 0
        do i = 1, natoms
            do j = 1, i-1
                ij       = ij+1
                r(1:3)   = x(1:3,i)-x(1:3,j)
                rij      = r(1)**2+r(2)**2+r(3)**2
                dist(ij) = sqrt(rij)
            enddo
        enddo
        if( iopt .ge. 1 )then
            ! build B matrix  dr/dx  from r = |xi-xj|
            b = 0.0_wp
            if( iopt .ge. 2 ) db = 0.0_wp
            ij = 0
            do i = 1, natoms
                do j = 1, i-1
                    ij              = ij+1
                    ix              = (i-1)*3
                    jx              = (j-1)*3
                    r(1:3)          = x(1:3,i)-x(1:3,j)
                    rrij            = r(1)**2+r(2)**2+r(3)**2
                    rrij            = 1.0_wp/sqrt(rrij)
                    r(1:3)          = rrij*r(1:3)
                    b(ij,ix+1:ix+3) = r(1:3)
                    b(ij,jx+1:jx+3) = -r(1:3)
                    if( iopt .ge. 2 )then
                        do k = 1, 3
                            db(ij,ix+k,ix+k) = rrij
                            db(ij,jx+k,jx+k) = rrij
                            db(ij,ix+k,jx+k) = -rrij
                            db(ij,jx+k,ix+k) = -rrij
                        enddo
                        do k = 1, 3
                            do l = 1, 3
                                d0               = r(k)*r(l)*rrij
                                db(ij,ix+k,ix+l) = db(ij,ix+k,ix+l) - d0
                                db(ij,jx+k,jx+l) = db(ij,jx+k,jx+l) - d0
                                db(ij,ix+k,jx+l) = db(ij,ix+k,jx+l) + d0
                                db(ij,jx+k,ix+l) = db(ij,jx+k,ix+l) + d0
                            enddo
                        enddo
                    endif
                enddo
            enddo
        endif
        if( iopt .ge. 0 )then
            e = shift
            do i = 1, nmorse
                do isym= 1, nsym
                    ii = imorse(isym, i)
                    if( ii == 0 ) exit
                    r(1) = dist(ii)
                    call v_morse(morse(1, i), morse(2, i), morse(3, i), r(1), v)
                    e = e+v
                enddo
            enddo
            do i = 1, ncoul
                do isym= 1, nsym
                    ii = icoul(isym, i)
                    if( ii == 0 ) exit
                    r(1) = dist(ii)
                    call v_coul(coul(1, i), coul(2, i), r(1), v)
                    e = e+v
                enddo
            enddo
            do i = 1, nxrep
                do isym= 1, nsym
                    ii = ixrep(isym, i)
                    if( ii == 0 ) exit
                    r(1) = dist(ii)
                    call v_xrep(xrep(1, i), xrep(2, i), xrep(3, i), r(1), v)
                    e = e+v
                enddo
            enddo
            do i = 1, ngauss
                do isym= 1, nsym
                    nd = igauss(1, isym, i)
                    if( nd == 0 ) exit
                    do j = 1, nd
                        ii   = igauss(j+1, isym, i)
                        r(j) = dist(ii)
                    enddo
                    call v_gauss(nd, r, gauss(1, i), gauss(ndim+1, i), gauss(2*ndim+1, i), gauss(2*ndim+2, i), v)
                    e = e+v
                enddo
            enddo
        endif
        if( iopt .ge. 1 )then
            g    = 0.0_wp
            gint = 0.0_wp
            do i = 1, nmorse
                do isym= 1, nsym
                    ii = imorse(isym, i)
                    if( ii == 0 ) exit
                    r(1) = dist(ii)
                    call f_morse(morse(1, i), morse(2, i), morse(3, i), r(1), gg(1))
                    gint(ii) = gint(ii)-gg(1)
                enddo
            enddo
            do i = 1, ncoul
                do isym= 1, nsym
                    ii = icoul(isym, i)
                    if( ii == 0 ) exit
                    r(1) = dist(ii)
                    call f_coul(coul(1, i), coul(2, i), r(1), gg(1))
                    gint(ii) = gint(ii)-gg(1)
                enddo
            enddo
            do i = 1, nxrep
                do isym= 1, nsym
                    ii = ixrep(isym, i)
                    if( ii == 0 ) exit
                    r(1) = dist(ii)
                    call f_xrep(xrep(1, i), xrep(2, i), xrep(3, i), r(1), gg(1))
                    gint(ii) = gint(ii)-gg(1)
                enddo
            enddo
            do i = 1, ngauss
                do isym= 1, nsym
                    nd = igauss(1, isym, i)
                    if( nd == 0 ) exit
                    do j = 1, nd
                        ii   = igauss(j+1, isym, i)
                        r(j) = dist(ii)
                    enddo
                    call f_gauss(nd, r, gauss(1, i), gauss(ndim+1, i), gauss(2*ndim+1, i), gg)
                    do j = 1, nd
                        ii       = igauss(j+1, isym, i)
                        gint(ii) = gint(ii)+gg(j)
                    enddo
                enddo
            enddo
            do i = 1, 3*natoms
                do j = 1, ndist
                    g(i) = g(i) + b(j,i) * gint(j)
                enddo
            enddo
        endif
        if( iopt .ge. 2 )then
            h    = 0.0_wp
            hint = 0.0_wp
            do i = 1, nmorse
                do isym= 1, nsym
                    ii = imorse(isym, i)
                    if( ii == 0 ) exit
                    r(1) = dist(ii)
                    call h_morse(morse(1, i), morse(2, i), morse(3, i), r(1), hh(1))
                    ij       = (ii-1)*ii/2 + ii
                    hint(ij) = hint(ij)+hh(1)
                enddo
            enddo
            do i = 1, ncoul
                do isym= 1, nsym
                    ii = icoul(isym, i)
                    if( ii == 0 ) exit
                    r(1) = dist(ii)
                    call h_coul(coul(1, i), coul(2, i), r(1), hh(1))
                    ij       = (ii-1)*ii/2 + ii
                    hint(ij) = hint(ij)+hh(1)
                enddo
            enddo
            do i = 1, nxrep
                do isym= 1, nsym
                    ii = ixrep(isym, i)
                    if( ii == 0 ) exit
                    r(1) = dist(ii)
                    call h_xrep(xrep(1, i), xrep(2, i), xrep(3, i), r(1), hh(1))
                    ij       = (ii-1)*ii/2 + ii
                    hint(ij) = hint(ij)+hh(1)
                enddo
            enddo
            do i = 1, ngauss
                do isym= 1, nsym
                    nd = igauss(1, isym, i)
                    if( nd == 0 ) exit
                    do j = 1, nd
                        ii   = igauss(j+1, isym, i)
                        r(j) = dist(ii)
                    enddo
                    call h_gauss(nd, r, gauss(1, i), gauss(ndim+1, i), gauss(2*ndim+1, i), hh)
                    ij = 0
                    do j = 1, nd
                        do k = 1, j
                            ij = ij+1
                            ii = igauss(j+1, isym, i)
                            jj = igauss(k+1, isym, i)
                            call iorder(jj, ii)
                            kl       = (ii-1)*ii/2 + jj
                            hint(kl) = hint(kl)+hh(ij)
                        enddo
                    enddo
                enddo
            enddo
            ij = 0
            do i = 1, 3*natoms
                do j = 1, i
                    ij = ij +1
                    do k = 1, ndist
                        h(ij) = h(ij) + db(k,j,i) * gint(k)
                    enddo
                enddo
            enddo
            ! reuse dB as scratch
            db(:,:,1) = 0.0_wp
            do j = 1, 3*natoms
                kl = 0
                do k = 1, ndist
                    do l = 1, k-1
                        kl        = kl+1
                        db(k,j,1) = db(k,j,1) + hint(kl) * b(l,j)
                        db(l,j,1) = db(l,j,1) + hint(kl) * b(k,j)
                    enddo
                    kl        = kl+1
                    db(k,j,1) = db(k,j,1) + hint(kl) * b(k,j)
                enddo
            enddo
            ij = 0
            do i = 1, 3*natoms
                do j = 1, i
                    ij = ij + 1
                    do k = 1, ndist
                        h(ij) = h(ij) + b(k,i) * db(k,j,1)
                    enddo
                enddo
            enddo
        endif
    end subroutine evaluate_native

    !> @brief Evaluate v_morse for the H3FO_Tew_2022 analytic potential.
    subroutine v_morse( re, al, de, r, v )
        !       Purpose: Evaluate V(r) = De**2*(1-exp(-alpha*(r-re)))**2       !
        real(wp), intent(in)  :: al, re, de, r
        real(wp), intent(out) :: v
        v = al*(re-r)
        v = exp(v) - 1.0_wp
        v = (de*v)**2
    end subroutine v_morse

    !> @brief Evaluate f_morse for the H3FO_Tew_2022 analytic potential.
    subroutine f_morse( re, al, de, r, f )
        !       Purpose: Evaluate -dV(r)/dr  = dV(r)/d re
        real(wp), intent(in)  :: al, re, de, r
        real(wp), intent(out) :: f
        f = al*(re-r)
        f = exp(f)
        f = f*(f-1.0_wp)
        f = f*2.0_wp*al*de**2
    end subroutine f_morse

    !> @brief Evaluate h_morse for the H3FO_Tew_2022 analytic potential.
    subroutine h_morse( re, al, de, r, h )
        !       Purpose: Evaluate d2V(r)/dr2
        real(wp), intent(in)  :: al, re, de, r
        real(wp), intent(out) :: h
        h = al*(re-r)
        h = exp(h)
        h = h*(2.0_wp*h-1.0_wp)
        h = h*2.0_wp*al**2*de**2
    end subroutine h_morse

    !> @brief Evaluate v_coul for the H3FO_Tew_2022 analytic potential.
    subroutine v_coul( q1, q2, r, v )
        !       Purpose: Evaluate V(r)
        real(wp), intent(in)  :: r, q1, q2
        real(wp), intent(out) :: v
        v = q1 * q2 / r
    end subroutine v_coul

    !> @brief Evaluate f_coul for the H3FO_Tew_2022 analytic potential.
    subroutine f_coul( q1, q2, r, f )
        !       Purpose: Evaluate -dV(r)/dr
        real(wp), intent(in)  :: r, q1, q2
        real(wp), intent(out) :: f
        f = q1 * q2 / r**2
    end subroutine f_coul

    !> @brief Evaluate h_coul for the H3FO_Tew_2022 analytic potential.
    subroutine h_coul( q1, q2, r, h )
        !       Purpose: Evaluate d2V(r)/dr2
        real(wp), intent(in)  :: r, q1, q2
        real(wp), intent(out) :: h
        h = 2.0_wp*q1*q2/r**3
    end subroutine h_coul

    !> @brief Evaluate v_xrep for the H3FO_Tew_2022 analytic potential.
    subroutine v_xrep( re, al, de, r, v )
        !       Purpose: Evaluate V(r)
        real(wp), intent(in)  :: r, al, re, de
        real(wp), intent(out) :: v
        v = al*(re-r)
        v = exp(v) * de**2
    end subroutine v_xrep

    !> @brief Evaluate f_xrep for the H3FO_Tew_2022 analytic potential.
    subroutine f_xrep( re, al, de, r, f )
        !       Purpose: Evaluate -dV(r)/dr  = dV(r)/d re
        real(wp), intent(in)  :: r, al, re, de
        real(wp), intent(out) :: f
        f = al*(re-r)
        f = exp(f)
        f = f*al*de**2
    end subroutine f_xrep

    !> @brief Evaluate h_xrep for the H3FO_Tew_2022 analytic potential.
    subroutine h_xrep( re, al, de, r, h )
        !       Purpose: Evaluate d2V(r)/dr2
        real(wp), intent(in)  :: r, al, re, de
        real(wp), intent(out) :: h
        h = al*(re-r)
        h = exp(h)
        h = h*al**2*de**2
    end subroutine h_xrep

    !> @brief Evaluate v_gauss for the H3FO_Tew_2022 analytic potential.
    subroutine v_gauss( ndim, r, x, alpha, d, shift, v )
        !       Purpose: Evaluate V(r) = D ( exp(- 0.5 a rx**2 ) - shift )
        integer, intent(in)   :: ndim
        real(wp), intent(in)  :: r(ndim), alpha(ndim), x(ndim), d, shift
        real(wp), intent(out) :: v
        integer               :: i
        v = 0.0_wp
        do i = 1, ndim
            v = v + (r(i)-x(i))**2 * alpha(i)
        enddo
        v = v*0.5_wp
        v = exp(-v) - shift
        v = v*d
    end subroutine v_gauss

    !> @brief Evaluate f_gauss for the H3FO_Tew_2022 analytic potential.
    subroutine f_gauss( ndim, r, x, alpha, d, v )
        !       Purpose: Evaluate dV/dr
        integer, intent(in)   :: ndim
        real(wp), intent(in)  :: r(ndim), alpha(ndim), x(ndim), d
        real(wp), intent(out) :: v(ndim)
        real(wp)              :: vv
        integer               :: i
        vv = 0.0_wp
        do i = 1, ndim
            vv = vv + (r(i)-x(i))**2 * alpha(i)
        enddo
        vv = vv*0.5_wp
        vv = exp(-vv)
        vv = vv*d
        do i = 1, ndim
            v(i) = - vv * (r(i)-x(i)) * alpha(i)
        enddo
    end subroutine f_gauss

    !> @brief Evaluate h_gauss for the H3FO_Tew_2022 analytic potential.
    subroutine h_gauss( ndim, r, x, alpha, d, h )
        !       Purpose: Evaluate d2V/dr2
        integer, intent(in)   :: ndim
        real(wp), intent(in)  :: r(ndim), alpha(ndim), x(ndim), d
        real(wp), intent(out) :: h(ndim*(ndim+1)/2)
        real(wp)              :: vv, fi
        integer               :: i, j, ij
        vv = 0.0_wp
        do i = 1, ndim
            vv = vv + (r(i)-x(i))**2 * alpha(i)
        enddo
        vv = vv*0.5_wp
        vv = exp(-vv)
        vv = vv*d
        ij = 0
        do i = 1, ndim
            fi = (r(i)-x(i)) * alpha(i)
            do j = 1, i
                ij    = ij+1
                h(ij) = vv * fi * (r(j)-x(j)) * alpha(j)
            enddo
            h(ij) = h(ij) - alpha(i) * vv
        enddo
    end subroutine h_gauss

    !> @brief Evaluate iorder for the H3FO_Tew_2022 analytic potential.
    subroutine iorder( i, j )
        ! permute i,j such that i < j
        implicit none
        integer :: i, j, k
        if( i > j )then
            k = i
            i = j
            j = k
        endif
    end subroutine iorder
end module chempotpy_h3fo_tew_2022_pes

!> @brief Expose H3FO_Tew_2022 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3fo_tew_2022_pes, only: evaluate
    use chempotpy_kinds,             only: wp
    implicit none
    integer, parameter    :: natoms = 5, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
