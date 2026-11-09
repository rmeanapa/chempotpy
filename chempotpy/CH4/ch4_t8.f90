!   System:                        CH4
!   Functional form:
!   Common name:                   Schwenke-Partridge T8
!   Number of derivatives:         0
!   Number of bodies:              5
!   Number of electronic surfaces: 1
!   Interface:                     Section-2
!   Data file:                     ch4T8datafile.tar
!
!   References:      D. W. Schwenke and H. Partridge, Spectrochim. Acta
!
!   Notes:
!      The data file (ch4T8datafile.tar) contains two files:
!      ch4pes is a data file read once by the PES routine
!      getpot_ch4t8.f90 is a calling routine that accepts cartesian coor
!      (the required input for the T8 PES), calls the surface, zeroes it
!      and returns the energy.
!      There are no derivatives.  One should use nuclear rather than ato
!

!> @brief Evaluate vibpot.
subroutine vibpot( rij, v, n, path )
    implicit none
    integer, parameter              :: id = 2000
    integer, intent(in)             :: n
    character(len=1024), intent(in) :: path
    real(8), intent(in)             :: rij(n,10)
    real(8), intent(out)            :: v(n,1)
    character(len=80)               :: title
    character(len=1024)             :: file_path1
    integer                         :: file_unit, i, ig, j
    integer                         :: nx(8), ixp(id,8,8), ixw(id,8,8)
    real(8)                         :: ccc(id,8), cs(6), cb(6), fmat(9,9)
    real(8)                         :: a02, a03, a04, a06, a08
    real(8)                         :: beta, beta0, bohr, bohri, const, damp
    real(8)                         :: d1, d2, d3, d4, d5, d6, ex, ex0, exa
    real(8)                         :: r01, r02, r03, r04, re0, reau
    real(8)                         :: rr1, rr2, rr3, rr4, s1, s2a, s2b
    real(8)                         :: s3x, s3y, s3z, s4x, s4y, s4z, sr12, srh
    real(8)                         :: t1, t2, t3, t4, t5, t6, t12, t22, t32, t42, t52, t62
    real(8)                         :: thrd, vs, y1, y2, y3, y4
    !
    !     ch4 potential based on Radau coordinates
    !
    !     rij(i,1) - ch1
    !     rij(i,2) - ch2
    !     rij(i,3) - ch3
    !     rij(i,4) - ch4
    !     rij(i,5) - h1ch2
    !     rij(i,6) - h1ch3
    !     rij(i,7) - h1ch4
    !     rij(i,8) - h2ch3
    !     rij(i,9) - h2ch4
    !     rij(i,10) - h3ch4
    !                  ALL UNITS ARE IN au,i.e., bohr, radian and hartree. Y
    bohr       = 0.529177249d0
    bohri      = 1d0/bohr
    file_path1 = trim(path)//"/CH4/ch4pes"
    open(newunit=file_unit, file=file_path1, form='formatted', status='old')
    read(file_unit, '(A80)') title
    read(file_unit, *) re0, ex0, beta0
    read(file_unit, *) r01, r02, r03, r04
    read(file_unit, *) a02, a03, a04, a06, a08
    re0 = re0*bohri
    r01 = r01*6d0; r02 = r02*6d0; r03 = r03*6d0; r04 = r04*6d0
    a02 = a02*4d0; a03 = a03*4d0; a04 = a04*4d0
    call setupvp(file_unit, nx, ixp, ixw, ccc, id, reau, ex, thrd, exa, srh, sr12, cs, cb, beta, const)
    close(file_unit)
    do ig=1, n
        t1  = cos(rij(ig,5))+thrd
        t2  = cos(rij(ig,6))+thrd
        t3  = cos(rij(ig,7))+thrd
        t4  = cos(rij(ig,8))+thrd
        t5  = cos(rij(ig,9))+thrd
        t6  = cos(rij(ig,10))+thrd
        rr1 = rij(ig,1)-re0
        rr2 = rij(ig,2)-re0
        rr3 = rij(ig,3)-re0
        rr4 = rij(ig,4)-re0
        d1  = exp(-beta0*(rr1**2+rr2**2))
        d2  = exp(-beta0*(rr1**2+rr3**2))
        d3  = exp(-beta0*(rr1**2+rr4**2))
        d4  = exp(-beta0*(rr2**2+rr3**2))
        d5  = exp(-beta0*(rr2**2+rr4**2))
        d6  = exp(-beta0*(rr3**2+rr4**2))
        y1  = 1d0-exp(-ex0*rr1)
        y2  = 1d0-exp(-ex0*rr2)
        y3  = 1d0-exp(-ex0*rr3)
        y4  = 1d0-exp(-ex0*rr4)
        t12 = t1*t1
        t22 = t2*t2
        t32 = t3*t3
        t42 = t4*t4
        t52 = t5*t5
        t62 = t6*t6
        vs=y1*(r01+y1*(r02+y1*(r03+y1*r04))) +y2*(r01+y2*(r02+y2*(r03+y2*r04))) +y3*(r01+y3*(r02+y3*(r03+y3*r04))) + &
            & y4*(r01+ y4*(r02+y4*(r03+y4*r04))) +d1*t12*(a02+t1*(a03+t1*(a04+t12*(a06+t12*a08)))) +d2*t22*(a02+t2* &
            & (a03+t2*(a04+t22* (a06+t22*a08)))) +d3*t32*(a02+t3*(a03+t3*(a04+t32*(a06+t32*a08)))) +d4*t42*(a02+t4* &
            & (a03+t4*(a04+t42*(a06+t42* a08)))) +d5*t52*(a02+t5*(a03+t5*(a04+t52*(a06+t52*a08)))) +d6*t62*(a02+t6* &
            & (a03+t6*(a04+t62*(a06+t62*a08))))
        d1 = rij(ig,1)-reau
        d2 = rij(ig,2)-reau
        d3 = rij(ig,3)-reau
        d4 = rij(ig,4)-reau
        if( ex .gt. 0d0 )then
            d1 = 1d0-exp(-ex*d1)
            d2 = 1d0-exp(-ex*d2)
            d3 = 1d0-exp(-ex*d3)
            d4 = 1d0-exp(-ex*d4)
        endif
        damp      = exp(-beta*(d1*d1+d2*d2+d3*d3+d4*d4))
        s1        = 0.5d0*(d1+d2+d3+d4)
        s2a       = sr12*(2d0*t1-t2-t3-t4-t5+2d0*t6)
        s2b       = 0.5d0*(t2-t3-t4+t5)
        s3x       = 0.5d0*(d1-d2+d3-d4)
        s3y       = 0.5d0*(d1-d2-d3+d4)
        s3z       = 0.5d0*(d1+d2-d3-d4)
        s4x       = srh*(t5-t2)
        s4y       = srh*(t4-t3)
        s4z       = srh*(t6-t1)
        fmat(1,1) = 1d0
        fmat(1,2) = 1d0
        fmat(1,3) = 1d0
        fmat(1,4) = 1d0
        fmat(1,5) = 1d0
        fmat(1,6) = 1d0
        fmat(1,7) = 1d0
        fmat(1,8) = 1d0
        fmat(1,9) = 1d0
        fmat(2,1) = s1
        fmat(2,2) = s2a
        fmat(2,3) = s2b
        fmat(2,4) = s3x
        fmat(2,5) = s3y
        fmat(2,6) = s3z
        fmat(2,7) = s4x
        fmat(2,8) = s4y
        fmat(2,9) = s4z
        do i=3, 9
            fmat(i,1) = fmat(i-1,1)*fmat(2,1)
            fmat(i,2) = fmat(i-1,2)*fmat(2,2)
            fmat(i,3) = fmat(i-1,3)*fmat(2,3)
            fmat(i,4) = fmat(i-1,4)*fmat(2,4)
            fmat(i,5) = fmat(i-1,5)*fmat(2,5)
            fmat(i,6) = fmat(i-1,6)*fmat(2,6)
            fmat(i,7) = fmat(i-1,7)*fmat(2,7)
            fmat(i,8) = fmat(i-1,8)*fmat(2,8)
            fmat(i,9) = fmat(i-1,9)*fmat(2,9)
        enddo
        do i=1, nx(1)
            vs = vs+ccc(i,1)*fmat(ixp(i,1,1),ixw(i,1,1))*damp
        enddo
        do i=1, nx(2)
            vs = vs+ccc(i,2)*fmat(ixp(i,1,2),ixw(i,1,2))*damp *fmat(ixp(i,2,2),ixw(i,2,2))
        enddo
        do i=1, nx(3)
            vs = vs+ccc(i,3)*fmat(ixp(i,1,3),ixw(i,1,3))*damp *fmat(ixp(i,2,3),ixw(i,2,3)) *fmat(ixp(i,3,3),ixw(i,3,3))
        enddo
        do i=1, nx(4)
            vs=vs+ccc(i,4)*fmat(ixp(i,1,4),ixw(i,1,4))*damp *fmat(ixp(i,2,4),ixw(i,2,4)) *fmat(ixp(i,3,4),ixw(i,3, &
                & 4)) * fmat(ixp(i,4,4),ixw(i,4,4))
        enddo
        do i=1, nx(5)
            vs=vs+ccc(i,5)*fmat(ixp(i,1,5),ixw(i,1,5))*damp *fmat(ixp(i,2,5),ixw(i,2,5)) *fmat(ixp(i,3,5),ixw(i,3, &
                & 5)) * fmat(ixp(i,4,5),ixw(i,4,5)) *fmat(ixp(i,5,5),ixw(i,5,5))
        enddo
        do i=1, nx(6)
            vs=vs+ccc(i,6)*fmat(ixp(i,1,6),ixw(i,1,6))*damp *fmat(ixp(i,2,6),ixw(i,2,6)) *fmat(ixp(i,3,6),ixw(i,3, &
                & 6)) * fmat(ixp(i,4,6),ixw(i,4,6)) *fmat(ixp(i,5,6),ixw(i,5,6)) *fmat(ixp(i,6,6),ixw(i,6,6))
        enddo
        do i=1, nx(7)
            vs=vs+ccc(i,7)*fmat(ixp(i,1,7),ixw(i,1,7))*damp *fmat(ixp(i,2,7),ixw(i,2,7)) *fmat(ixp(i,3,7),ixw(i,3, &
                & 7)) * fmat(ixp(i,4,7),ixw(i,4,7)) *fmat(ixp(i,5,7),ixw(i,5,7)) *fmat(ixp(i,6,7),ixw(i,6,7)) * &
                & fmat(ixp(i,7,7),ixw(i,7,7))
        enddo
        do i=1, nx(8)
            vs=vs+ccc(i,8)*fmat(ixp(i,1,8),ixw(i,1,8))*damp *fmat(ixp(i,2,8),ixw(i,2,8)) *fmat(ixp(i,3,8),ixw(i,3, &
                & 8)) * fmat(ixp(i,4,8),ixw(i,4,8)) *fmat(ixp(i,5,8),ixw(i,5,8)) *fmat(ixp(i,6,8),ixw(i,6,8)) * &
                & fmat(ixp(i,7,8),ixw(i,7,8)) *fmat(ixp(i,8,8),ixw(i,8,8))
        enddo
        vs      = vs+const*damp
        v(ig,1) = vs
    enddo
end subroutine vibpot

!> @brief Evaluate setupvp.
subroutine setupvp( file_unit, nx, ixp, ixw, ccc, idx, reau, ex, thrd, exa, srh, sr12, cs, cb, beta, const )
    implicit none
    integer, parameter   :: id = 6400
    integer, intent(in)  :: file_unit, idx
    integer, intent(out) :: nx(8), ixp(idx,8,8), ixw(idx,8,8)
    real(8), intent(out) :: ccc(idx,8), reau, ex, thrd, exa, srh, sr12
    real(8), intent(out) :: cs(6), cb(6), beta, const
    integer              :: i, i1, i2, i3, i4, i5, i6, i7, i8, i9
    integer              :: ii, in, ios, iradau, ix1, ix2, j, j1, j2, lm
    integer              :: nexp, nnz, nxzero
    integer              :: ix(id,9)
    real(8)              :: cc(id), coef, re
    read(file_unit, *) re, ex, exa, beta, iradau, lm
    sr12 = sqrt(1d0/12d0)
    srh  = sqrt(0.5d0)
    cs   = 0.0d0; cb = 0.0d0; ccc = 0.0d0
    ixp  = 0; ixw = 0
    nexp = 0
    ix1  = 0
    ix2  = 0
    do
        read(file_unit, *, iostat=ios) i1, i2, i3, i4, i5, i6, i7, i8, i9, coef
        if( ios < 0 ) exit
        if( ios > 0 ) error stop 'SETUPVP: invalid CH4 coefficient record'
        in = i3+i4+i5+i6+i7+i8+i9
        j1 = in+i2; j2 = in+i1
        if( j1 == 0 .and. lm /= 0 )then
            cs(i1) = coef*6d0; ix1 = max(ix1, i1); coef = 0d0
        else if( j2 == 0 .and. lm /= 0 )then
            cb(i2) = coef*4d0; ix2 = max(ix2, i2); coef = 0d0
        endif
        nexp = nexp+1
        if( nexp > id ) error stop 'SETUPVP: coefficient capacity exceeded'
        ix(nexp,:) = [i1, i2, i3, i4, i5, i6, i7, i8, i9]
        cc(nexp)   = coef
    enddo
    reau = re*1.889725989d0
    do i=1, nexp
        do j=1, 9
            ix(i,j) = ix(i,j)+1
        enddo
    enddo
    thrd   = 1d0/3d0
    nx     = 0
    nxzero = 0
    do i=1, nexp
        nnz = 0
        do j=1, 9
            if( ix(i,j) .gt. 1 ) nnz=nnz+1
        enddo
        if( nnz .eq. 0 )then
            nxzero = nxzero+1
            const  = cc(i)
        else
            nx(nnz) = nx(nnz)+1
            if( nx(nnz) .gt. idx ) stop 'idx'
            ccc(nx(nnz),nnz) = cc(i)
            ii               = 1
            do j=1, 9
                if( ix(i,j) .gt. 1 )then
                    ixp(nx(nnz),ii,nnz) = ix(i,j)
                    ixw(nx(nnz),ii,nnz) = j
                    ii                  = ii+1
                endif
            enddo
        endif
    enddo
end subroutine setupvp
