!   System:                        CH4
!   Functional form:
!   Common name:                   Schwenke-Partridge T8
!   Number of derivatives:         0
!   Number of bodies:              5
!   Number of electronic surfaces: 1
!   Interface:                     Section-2
!   Data file:                     ch4T8datafile.tar
!   References:      D. W. Schwenke and H. Partridge, Spectrochim. Acta
!   Notes:
!      The data file (ch4T8datafile.tar) contains two files:
!      the released ch4pes data are embedded in the parameter module
!      getpot_ch4t8.f90 is a calling routine that accepts cartesian coor
!      (the required input for the T8 PES), calls the surface, zeroes it
!      and returns the energy.
!      There are no derivatives.  One should use nuclear rather than ato

!> @brief Evaluate vibpot.
subroutine vibpot( rij, v, n, path )
    use chempotpy_ch4_t8_parameters, only: ccc => coefficients, const => constant, &
        & ixp => exponent_indices, ixw => coordinate_indices, nx, term_offsets, parameter_a02 => a02, &
        & parameter_a03 => a03, parameter_a04 => a04, parameter_a06 => a06, &
        & parameter_a08 => a08, parameter_beta => beta, parameter_beta0 => beta0, &
        & parameter_ex => ex, parameter_ex0 => ex0, parameter_r01 => r01, parameter_r02 => r02, &
        & parameter_r03 => r03, parameter_r04 => r04, parameter_re0 => re0, parameter_reau => reau
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)             :: n
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: rij(n,10)
    real(wp), intent(out)           :: v(n,1)
    integer                         :: i, ig, j
    real(wp)                        :: fmat(9,9)
    real(wp)                        :: a02, a03, a04, a06, a08
    real(wp)                        :: beta, beta0, bohr, bohri, damp
    real(wp)                        :: d1, d2, d3, d4, d5, d6, ex, ex0
    real(wp)                        :: r01, r02, r03, r04, re0, reau
    real(wp)                        :: rr1, rr2, rr3, rr4, s1, s2a, s2b
    real(wp)                        :: s3x, s3y, s3z, s4x, s4y, s4z, sr12, srh
    real(wp)                        :: t1, t2, t3, t4, t5, t6, t12, t22, t32, t42, t52, t62
    real(wp)                        :: thrd, vs, y1, y2, y3, y4
    !     ch4 potential based on Radau coordinates
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
    bohr  = 0.529177249e0_wp
    bohri = 1e0_wp/bohr
    re0   = parameter_re0; ex0 = parameter_ex0; beta0 = parameter_beta0
    r01   = parameter_r01; r02 = parameter_r02; r03 = parameter_r03; r04 = parameter_r04
    a02   = parameter_a02; a03 = parameter_a03; a04 = parameter_a04
    a06   = parameter_a06; a08 = parameter_a08
    reau  = parameter_reau; ex = parameter_ex; beta = parameter_beta
    re0   = re0*bohri
    r01   = r01*6e0_wp; r02 = r02*6e0_wp; r03 = r03*6e0_wp; r04 = r04*6e0_wp
    a02   = a02*4e0_wp; a03 = a03*4e0_wp; a04 = a04*4e0_wp
    thrd  = 1e0_wp/3e0_wp; sr12 = sqrt(1e0_wp/12e0_wp); srh = sqrt(0.5e0_wp)
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
        y1  = 1e0_wp-exp(-ex0*rr1)
        y2  = 1e0_wp-exp(-ex0*rr2)
        y3  = 1e0_wp-exp(-ex0*rr3)
        y4  = 1e0_wp-exp(-ex0*rr4)
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
        if( ex .gt. 0e0_wp )then
            d1 = 1e0_wp-exp(-ex*d1)
            d2 = 1e0_wp-exp(-ex*d2)
            d3 = 1e0_wp-exp(-ex*d3)
            d4 = 1e0_wp-exp(-ex*d4)
        endif
        damp      = exp(-beta*(d1*d1+d2*d2+d3*d3+d4*d4))
        s1        = 0.5e0_wp*(d1+d2+d3+d4)
        s2a       = sr12*(2e0_wp*t1-t2-t3-t4-t5+2e0_wp*t6)
        s2b       = 0.5e0_wp*(t2-t3-t4+t5)
        s3x       = 0.5e0_wp*(d1-d2+d3-d4)
        s3y       = 0.5e0_wp*(d1-d2-d3+d4)
        s3z       = 0.5e0_wp*(d1+d2-d3-d4)
        s4x       = srh*(t5-t2)
        s4y       = srh*(t4-t3)
        s4z       = srh*(t6-t1)
        fmat(1,1) = 1e0_wp
        fmat(1,2) = 1e0_wp
        fmat(1,3) = 1e0_wp
        fmat(1,4) = 1e0_wp
        fmat(1,5) = 1e0_wp
        fmat(1,6) = 1e0_wp
        fmat(1,7) = 1e0_wp
        fmat(1,8) = 1e0_wp
        fmat(1,9) = 1e0_wp
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
            vs = vs+ccc(i, 1)*fmat(ixp(term_offsets(1)+(i-1)*1+1),ixw(term_offsets(1)+(i-1)*1+1))*damp
        enddo
        do i=1, nx(2)
            vs = vs+ccc(i, 2)*fmat(ixp(term_offsets(2)+(i-1)*2+1),ixw(term_offsets(2)+(i-1)*2+1))*damp * &
                & fmat(ixp(term_offsets(2)+(i-1)*2+2),ixw(term_offsets(2)+(i-1)*2+2))
        enddo
        do i=1, nx(3)
            vs = vs+ccc(i, 3)*fmat(ixp(term_offsets(3)+(i-1)*3+1),ixw(term_offsets(3)+(i-1)*3+1))*damp * &
                & fmat(ixp(term_offsets(3)+(i-1)*3+2),ixw(term_offsets(3)+(i-1)*3+2)) *fmat(ixp(term_offsets(3)+(i-1)* &
                & 3+3),ixw(term_offsets(3)+(i-1)*3+3))
        enddo
        do i=1, nx(4)
            vs=vs+ccc(i, 4)*fmat(ixp(term_offsets(4)+(i-1)*4+1),ixw(term_offsets(4)+(i-1)*4+1))*damp * &
                & fmat(ixp(term_offsets(4)+(i-1)*4+2),ixw(term_offsets(4)+(i-1)*4+2)) *fmat(ixp(term_offsets(4)+(i-1)* &
                & 4+3),ixw(term_offsets(4)+(i-1)*4+3)) * fmat(ixp(term_offsets(4)+(i-1)*4+4),ixw(term_offsets(4)+(i- &
                & 1)*4+4))
        enddo
        do i=1, nx(5)
            vs=vs+ccc(i, 5)*fmat(ixp(term_offsets(5)+(i-1)*5+1),ixw(term_offsets(5)+(i-1)*5+1))*damp * &
                & fmat(ixp(term_offsets(5)+(i-1)*5+2),ixw(term_offsets(5)+(i-1)*5+2)) *fmat(ixp(term_offsets(5)+(i-1)* &
                & 5+3),ixw(term_offsets(5)+(i-1)*5+3)) * fmat(ixp(term_offsets(5)+(i-1)*5+4),ixw(term_offsets(5)+(i- &
                & 1)*5+4)) *fmat(ixp(term_offsets(5)+(i-1)*5+5),ixw(term_offsets(5)+(i-1)*5+5))
        enddo
        do i=1, nx(6)
            vs=vs+ccc(i, 6)*fmat(ixp(term_offsets(6)+(i-1)*6+1),ixw(term_offsets(6)+(i-1)*6+1))*damp * &
                & fmat(ixp(term_offsets(6)+(i-1)*6+2),ixw(term_offsets(6)+(i-1)*6+2)) *fmat(ixp(term_offsets(6)+(i-1)* &
                & 6+3),ixw(term_offsets(6)+(i-1)*6+3)) * fmat(ixp(term_offsets(6)+(i-1)*6+4),ixw(term_offsets(6)+(i- &
                & 1)*6+4)) *fmat(ixp(term_offsets(6)+(i-1)*6+5),ixw(term_offsets(6)+(i-1)*6+5)) * &
                & fmat(ixp(term_offsets(6)+(i-1)*6+6),ixw(term_offsets(6)+(i-1)*6+6))
        enddo
        do i=1, nx(7)
            vs=vs+ccc(i, 7)*fmat(ixp(term_offsets(7)+(i-1)*7+1),ixw(term_offsets(7)+(i-1)*7+1))*damp * &
                & fmat(ixp(term_offsets(7)+(i-1)*7+2),ixw(term_offsets(7)+(i-1)*7+2)) *fmat(ixp(term_offsets(7)+(i-1)* &
                & 7+3),ixw(term_offsets(7)+(i-1)*7+3)) * fmat(ixp(term_offsets(7)+(i-1)*7+4),ixw(term_offsets(7)+(i- &
                & 1)*7+4)) *fmat(ixp(term_offsets(7)+(i-1)*7+5),ixw(term_offsets(7)+(i-1)*7+5)) * &
                & fmat(ixp(term_offsets(7)+(i-1)*7+6),ixw(term_offsets(7)+(i-1)*7+6)) * fmat(ixp(term_offsets(7)+(i- &
                & 1)*7+7),ixw(term_offsets(7)+(i-1)*7+7))
        enddo
        do i=1, nx(8)
            vs=vs+ccc(i, 8)*fmat(ixp(term_offsets(8)+(i-1)*8+1),ixw(term_offsets(8)+(i-1)*8+1))*damp * &
                & fmat(ixp(term_offsets(8)+(i-1)*8+2),ixw(term_offsets(8)+(i-1)*8+2)) *fmat(ixp(term_offsets(8)+(i-1)* &
                & 8+3),ixw(term_offsets(8)+(i-1)*8+3)) * fmat(ixp(term_offsets(8)+(i-1)*8+4),ixw(term_offsets(8)+(i- &
                & 1)*8+4)) *fmat(ixp(term_offsets(8)+(i-1)*8+5),ixw(term_offsets(8)+(i-1)*8+5)) * &
                & fmat(ixp(term_offsets(8)+(i-1)*8+6),ixw(term_offsets(8)+(i-1)*8+6)) * fmat(ixp(term_offsets(8)+(i- &
                & 1)*8+7),ixw(term_offsets(8)+(i-1)*8+7)) *fmat(ixp(term_offsets(8)+(i-1)*8+8),ixw(term_offsets(8)+(i- &
                & 1)*8+8))
        enddo
        vs      = vs+const*damp
        v(ig,1) = vs
    enddo
end subroutine vibpot
