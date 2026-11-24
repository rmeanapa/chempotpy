!> @brief Evaluate the released degree-four S4 x S4 x S2 invariant basis for the formic-acid dimer PES.
! The generated algebra is reduced to the degree required by the published fit.
! Source: C. Qu and J. M. Bowman, Phys. Chem. Chem. Phys. 18, 24835 (2016).
module chempotpy_c2h4o4_invariants
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_fad_fit
    integer, parameter :: nk              = 10
    integer, parameter :: nr              = nk*(nk-1)/2
    integer, parameter :: nkj(0:2)        = [4, 4, 2]
    integer, parameter :: mg442_dnpb(0:4) = [1, 6, 34, 145, 567]
    integer, parameter :: mg442_npb(0:4)  = [1, 7, 41, 186, 753]
    integer, parameter :: mg442_nsc(0:4)  = [1, 1, 8, 70, 380]
    integer, parameter :: mg442_nb(0:4)   = [1, 7, 48, 297, 1784]

    contains

    !> @brief Flatten a symmetric distance matrix in reverse-lexicographic pair order.
    subroutine mgx_mkrl1d( nk, d, x )
        ! Simple revlex order
        integer, intent (in)   :: nk
        real(wp), intent (in)  :: d(0:,0:)
        real(wp), intent (out) :: x(0:)
        integer                :: i, j, k
        k = 0
        do j = 0, nk-1
            do i = 0, j-1
                x(k) = d(i,j)
                k    = k+1
            enddo
        enddo
    end subroutine mgx_mkrl1d

    !> @brief Flatten block-grouped distances in the ordering required by the invariant basis.
    subroutine mgx_mk1d( nkj, d, x )
        ! Block revlex order; block sizes nkj(0:)
        integer, intent (in)   :: nkj(0:)
        real(wp), intent (in)  :: d(0:,0:)
        real(wp), intent (out) :: x(0:)
        integer                :: i, j, k, l0, l1, j0, j1
        ! For reasons of speed we use the code for simple revlex in cases
        ! where it will give the same result
        if( size(nkj) .le. 1 )then
            call mgx_mkrl1d(sum(nkj), d, x)
        elseif( nkj(1) .lt. 3.and.all(nkj(2:) .lt. 2) )then
            call mgx_mkrl1d(sum(nkj), d, x)
        else
            k  = 0
            j1 = 0
            do l1 = 0, size(nkj)-1
                j0 = 0
                do l0 = 0, l1
                    do j = j1, j1+nkj(l1)-1
                        do i = j0, min(j-1, j0+nkj(l0)-1)
                            x(k) = d(i,j)
                            k    = k+1
                        enddo
                    enddo
                    j0 = j0+nkj(l0)
                enddo
                j1 = j1+nkj(l1)
            enddo
        endif
    end subroutine mgx_mk1d

    !> @brief Symmetrize a distance matrix and form its requested elementwise powers.
    subroutine mgx_setd( r, d, d2, d3, d4, d5, d6, d7, d8, d9 )
        real(wp), intent (in)            :: r(0:,0:)
        real(wp), intent (out)           :: d(0:,0:)
        real(wp), optional, intent (out) :: d2(0:,0:), d3(0:,0:), &
            d4(0:,0:), d5(0:, 0:), d6(0:, 0:), d7(0:, 0:), d8(0:, 0:), d9(0:, 0:)
        integer :: i, j
        do j = 0, size(r, 2)-1
            do i = 0, j-1
                d(i,j) = (r(i,j)+r(j,i))/2
                d(j,i) = d(i,j)
            enddo
            d(j,j) = 0
        enddo
        if( present(d2) )then
            d2 = d*d
            if( present(d3) )then
                d3 = d2*d
                if( present(d4) )then
                    d4 = d3*d
                    if( present(d5) )then
                        d5 = d4*d
                        if( present(d6) )then
                            d6 = d5*d
                            if( present(d7) )then
                                d7 = d6*d
                                if( present(d8) )then
                                    d8 = d7*d
                                    if( present(d9) )then
                                        d9 = d8*d
                                    endif
                                endif
                            endif
                        endif
                    endif
                endif
            endif
        endif
    end subroutine mgx_setd

    !> @brief Build the six primitive invariants for four permutable atoms.
    subroutine cg4_prims( x, u )
        real(wp), intent (in)  :: x(0:)
        real(wp), intent (out) :: u(0:)
        real(wp)               :: t1(0:3)
        ! Recall revlex order,
        ! x(0:5) = (/ d01, d02, d12, d03, d13, d23 /)
        t1(0) = (x(0)+x(1)+x(3))/3
        t1(1) = (x(0)+x(2)+x(4))/3
        t1(2) = (x(1)+x(2)+x(5))/3
        t1(3) = (x(3)+x(4)+x(5))/3
        u(0)  = sum(x)/size(x)
        u(1)  = sum(t1**2)/size(t1)
        u(2)  = sum(x**2)/size(x)
        u(3)  = sum(t1**3)/size(t1)
        u(4)  = sum(x**3)/size(x)
        u(5)  = sum(t1**4)/size(t1)
    end subroutine cg4_prims

    !> @brief Build the sixteen primitive invariants coupling two four-atom groups.
    subroutine bg44_prims( x, u )
        real(wp), intent (in)  :: x(0:)
        real(wp), intent (out) :: u(0:)
        integer, parameter     :: m = 4, n=4
        integer                :: i, j
        real(wp)               :: e(0:m-1), f(0:n-1), t0, t1, t2, t3
        do i = 0, m-1
            e(i) = sum(x(i:i+m*(n-1):m))/n
        enddo
        do j = 0, n-1
            f(j) = sum(x(m*j:m*(j+1)-1))/m
        enddo
        t0 = 0 ; t1 = 0 ; t2 = 0 ; t3 = 0
        do j = 0, n-1
            do i = 0, m-1
                t0 = t0+x(i+m*j)**2*f(j)
                t1 = t1+x(i+m*j)**2*e(i)
                t2 = t2+x(i+m*j)**3*f(j)
                t3 = t3+x(i+m*j)**3*e(i)
            enddo
        enddo
        u(0)  = sum(x)/size(x)
        u(1)  = sum(e**2)/size(e)
        u(2)  = sum(f**2)/size(f)
        u(3)  = sum(x**2)/size(x)
        u(4)  = sum(e**3)/size(e)
        u(5)  = sum(f**3)/size(f)
        u(6)  = sum(x**3)/size(x)
        u(7)  = t0/(m*n)
        u(8)  = t1/(m*n)
        u(9)  = sum(e**4)/size(e)
        u(10) = sum(f**4)/size(f)
        u(11) = sum(x**4)/size(x)
        u(12) = t2/(m*n)
        u(13) = t3/(m*n)
        u(14) = sum(x**6)/size(x)
        u(15) = sum(x**12)/size(x)
    end subroutine bg44_prims

    !> @brief Build the eight primitive invariants coupling four- and two-atom groups.
    subroutine bg42_prims( x, u )
        real(wp), intent (in)  :: x(0:)
        real(wp), intent (out) :: u(0:)
        integer, parameter     :: m = 4, n=2
        integer                :: i, j
        real(wp)               :: e(0:m-1), f(0:n-1)
        do i = 0, m-1
            e(i) = sum(x(i:i+m*(n-1):m))/n
        enddo
        do j = 0, n-1
            f(j) = sum(x(m*j:m*(j+1)-1))/m
        enddo
        u(0) = sum(x)/size(x)
        u(1) = sum(e**2)/size(e)
        u(2) = sum(f**2)/size(f)
        u(3) = sum(x**2)/size(x)
        u(4) = sum(e**3)/size(e)
        u(5) = sum(e**4)/size(e)
        u(6) = sum(x**4)/size(x)
        u(7) = sum(x**6)/size(x)
    end subroutine bg42_prims

    !> @brief Build the single primitive invariant for two permutable atoms.
    subroutine cg2_prims( x, u )
        real(wp), intent (in)  :: x(0:)
        real(wp), intent (out) :: u(0:)
        ! There is just one variable
        u = x
    end subroutine cg2_prims

    !> @brief Assemble the 45 primitive invariants for the four-four-two atom partition.
    subroutine mg442_prims( r, u )
        real(wp), intent (in)  :: r(0:,0:)
        real(wp), intent (out) :: u(0:)
        integer, parameter     :: m = 4, n=4, l=2, m2=m*(m-1)/2, n2=n*(n-1)/2, &
            l2 = l*(l-1)/2
        real(wp) :: x(0:nr-1), u0(0:m2-1), u1(0:m*n-1), u2(0:n2-1), &
            u3(0:m*l-1), u4(0:n*l-1), u5(0:l2-1)
        call mgx_mk1d(nkj, r, x)
        call cg4_prims(x(0:m2-1), u0)
        call bg44_prims(x(m2:m2+m*n-1), u1)
        call cg4_prims(x(m2+m*n:m2+m*n+n2-1), u2)
        call bg42_prims(x(m2+m*n+n2:m2+m*n+n2+m*l-1), u3)
        call bg42_prims(x(m2+m*n+n2+m*l:m2+m*n+n2+m*l+n*l-1), u4)
        call cg2_prims(x(m2+m*n+n2+m*l+n*l:nr-1), u5)
        u = (/ u0(0), u1(0), u2(0), u3(0), u4(0), u5(0), u0(1), u0(2), u1(1), u1(2), u1(3), u2(1), u2(2), u3(1), &
            & u3(2), u3(3), u4(1), u4(2), u4(3), u0(3), u0(4), u1(4), u1(5), u1(6), u1(7), u1(8), u2(3), u2(4), u3(4), &
            & u4(4), u0(5), u1(9), u1(10), u1(11), u1(12), u1(13), u2(5), u3(5), u3(6), u4(5), u4(6), u1(14), u3(7), &
            & u4(7), u1(15) /)
    end subroutine mg442_prims

    !> @brief Build the primary invariant products through polynomial degree four.
    subroutine mg442_prib( mxd, u, w )
        integer, intent (in)   :: mxd
        real(wp), intent (in)  :: u(0:)
        real(wp), intent (out) :: w(0:mg442_npb(mxd)-1)
        ! The following code was obtained using these parameters
        ! prims, dnpr(0:*): 0 6 13 11 11 0 3 0 0 0
        ! prib,  dnpb(0:*): 1 6 34 145 567 1966 6383 19293 55448 151617
        ! constant term
        w(0) = 1
        ! terms of degree 1
        if( 1 .le. mxd )then
            w(1) = u(0)*w(0)
            w(2) = u(1)*w(0)
            w(3) = u(2)*w(0)
            w(4) = u(3)*w(0)
            w(5) = u(4)*w(0)
            w(6) = u(5)*w(0)
        endif
        ! terms of degree 2
        if( 2 .le. mxd )then
            w(7:12)  = u(0)*w(1:6)
            w(13:17) = u(1)*w(2:6)
            w(18:21) = u(2)*w(3:6)
            w(22:24) = u(3)*w(4:6)
            w(25:26) = u(4)*w(5:6)
            w(27)    = u(5)*w(6)
            w(28)    = u(6)*w(0)
            w(29)    = u(7)*w(0)
            w(30)    = u(8)*w(0)
            w(31)    = u(9)*w(0)
            w(32)    = u(10)*w(0)
            w(33)    = u(11)*w(0)
            w(34)    = u(12)*w(0)
            w(35)    = u(13)*w(0)
            w(36)    = u(14)*w(0)
            w(37)    = u(15)*w(0)
            w(38)    = u(16)*w(0)
            w(39)    = u(17)*w(0)
            w(40)    = u(18)*w(0)
        endif
        ! terms of degree 3
        if( 3 .le. mxd )then
            w(41:74)   = u(0)*w(7:40)
            w(75:102)  = u(1)*w(13:40)
            w(103:125) = u(2)*w(18:40)
            w(126:144) = u(3)*w(22:40)
            w(145:160) = u(4)*w(25:40)
            w(161:174) = u(5)*w(27:40)
            w(175)     = u(19)*w(0)
            w(176)     = u(20)*w(0)
            w(177)     = u(21)*w(0)
            w(178)     = u(22)*w(0)
            w(179)     = u(23)*w(0)
            w(180)     = u(24)*w(0)
            w(181)     = u(25)*w(0)
            w(182)     = u(26)*w(0)
            w(183)     = u(27)*w(0)
            w(184)     = u(28)*w(0)
            w(185)     = u(29)*w(0)
        endif
        ! terms of degree 4
        if( 4 .le. mxd )then
            w(186:330) = u(0)*w(41:185)
            w(331:441) = u(1)*w(75:185)
            w(442:524) = u(2)*w(103:185)
            w(525:584) = u(3)*w(126:185)
            w(585:625) = u(4)*w(145:185)
            w(626:650) = u(5)*w(161:185)
            w(651:663) = u(6)*w(28:40)
            w(664:675) = u(7)*w(29:40)
            w(676:686) = u(8)*w(30:40)
            w(687:696) = u(9)*w(31:40)
            w(697:705) = u(10)*w(32:40)
            w(706:713) = u(11)*w(33:40)
            w(714:720) = u(12)*w(34:40)
            w(721:726) = u(13)*w(35:40)
            w(727:731) = u(14)*w(36:40)
            w(732:735) = u(15)*w(37:40)
            w(736:738) = u(16)*w(38:40)
            w(739:740) = u(17)*w(39:40)
            w(741)     = u(18)*w(40)
            w(742)     = u(30)*w(0)
            w(743)     = u(31)*w(0)
            w(744)     = u(32)*w(0)
            w(745)     = u(33)*w(0)
            w(746)     = u(34)*w(0)
            w(747)     = u(35)*w(0)
            w(748)     = u(36)*w(0)
            w(749)     = u(37)*w(0)
            w(750)     = u(38)*w(0)
            w(751)     = u(39)*w(0)
            w(752)     = u(40)*w(0)
        endif
        ! terms of degree 5
        ! terms of degree 6
        ! terms of degree 7
        ! terms of degree 8
        ! terms of degree 9
    end subroutine mg442_prib

    !> @brief Accumulate the irreducible secondary invariants through degree four.
    subroutine mg442_isecs( mxd, r, pv )
        integer, intent (in)   :: mxd
        real(wp), intent (in)  :: r(0:,0:)
        real(wp), intent (out) :: pv(0:)
        integer, parameter     :: m = 4, n=4, l=2, m2=m*(m-1), mm=m2*(m-2), &
            n2 = n*(n-1), nn=n2*(n-2), ll=l
        integer  :: i0, i1, i2, i3, j0, j1, j2, j3, k0, k1
        real(wp) :: d(0:nk-1,0:nk-1), d2(0:nk-1,0:nk-1), d3(0:nk-1,0:nk-1), d4(0:nk-1,0:nk-1)
        call mgx_setd(r, d, d2, d3, d4)
        pv = 0
        do k0 = m+n, m+n+l-1
            do k1 = m+n, m+n+l-1
                if( k1 .ne. k0 )then
                    do j0 = m, m+n-1
                        do i0 = 0, m-1
                            call mg442_i0j0(mxd, m*n*ll)
                            do i1 = 0, m-1
                                if( i1 .ne. i0 )then
                                    call mg442_i1j0(mxd, m2*n*ll)
                                    do i2 = 0, m-1
                                        if( i2 .ne. i0.and.i2 .ne. i1 )then
                                            do i3 = 0, m-1
                                                if( i3 .ne. i0.and.i3 .ne. i1.and.i3 .ne. i2 )then
                                                    call mg442_iij0(mxd, mm*n*ll)
                                                endif
                                            enddo
                                        endif
                                    enddo
                                endif
                            enddo
                        enddo
                        do j1 = m, m+n-1
                            if( j1 .ne. j0 )then
                                do i0 = 0, m-1
                                    call mg442_i0j1(mxd, m*n2*ll)
                                    do i1 = 0, m-1
                                        if( i1 .ne. i0 )then
                                            call mg442_i1j1(mxd, m2*n2*ll)
                                            do i2 = 0, m-1
                                                if( i2 .ne. i0.and.i2 .ne. i1 )then
                                                    do i3 = 0, m-1
                                                        if( i3 .ne. i0.and.i3 .ne. i1.and.i3 .ne. i2 )then
                                                            call mg442_iij1(mxd, mm*n2*ll)
                                                        endif
                                                    enddo
                                                endif
                                            enddo
                                        endif
                                    enddo
                                enddo
                                do j2 = m, m+n-1
                                    if( j2 .ne. j0.and.j2 .ne. j1 )then
                                        do j3 = m, m+n-1
                                            if( j3 .ne. j0.and.j3 .ne. j1.and.j3 .ne. j2 )then
                                                do i0 = 0, m-1
                                                    call mg442_i0jj(mxd, m*nn*ll)
                                                    do i1 = 0, m-1
                                                        if( i1 .ne. i0 )then
                                                            call mg442_i1jj(mxd, m2*nn*ll)
                                                            do i2 = 0, m-1
                                                                if( i2 .ne. i0.and.i2 .ne. i1 )then
                                                                    do i3 = 0, m-1
                                                                        if( i3 .ne. i0.and.i3 .ne. i1.and.i3 .ne. i2 )then
                                                                        endif
                                                                    enddo
                                                                endif
                                                            enddo
                                                        endif
                                                    enddo
                                                enddo
                                            endif
                                        enddo
                                    endif
                                enddo
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo

        contains

        !> @brief Accumulate secondary invariants for the single-H and single-O index class.
        subroutine mg442_i0j0( mxd, mu )
            integer, intent (in) :: mxd, mu
            if( 2 .le. mxd )then
                pv(3) = pv(3)+d(i0,j0)*d(i0,k0)/mu
                pv(4) = pv(4)+d(i0,j0)*d(j0,k0)/mu
                pv(6) = pv(6)+d(i0,k0)*d(j0,k0)/mu
            endif
            if( 3 .le. mxd )then
                pv(28) = pv(28)+d2(i0,j0)*d(i0,k0)/mu
                pv(34) = pv(34)+d(i0,j0)*d2(i0,k0)/mu
                pv(35) = pv(35)+d3(i0,k0)/mu
                pv(41) = pv(41)+d(i0,j0)*d(i0,k0)*d(i0,k1)/mu
                pv(44) = pv(44)+d2(i0,j0)*d(j0,k0)/mu
                pv(52) = pv(52)+d(i0,j0)*d(i0,k0)*d(j0,k0)/mu
                pv(56) = pv(56)+d2(i0,k0)*d(j0,k0)/mu
                pv(57) = pv(57)+d(i0,j0)*d(i0,k1)*d(j0,k0)/mu
                pv(58) = pv(58)+d(i0,j0)*d2(j0,k0)/mu
                pv(60) = pv(60)+d(i0,k0)*d2(j0,k0)/mu
                pv(61) = pv(61)+d3(j0,k0)/mu
                pv(66) = pv(66)+d(i0,j0)*d(j0,k0)*d(j0,k1)/mu
            endif
            if( 4 .le. mxd )then
                pv(147) = pv(147)+d3(i0,j0)*d(i0,k0)/mu
                pv(175) = pv(175)+d2(i0,j0)*d2(i0,k0)/mu
                pv(181) = pv(181)+d(i0,j0)*d3(i0,k0)/mu
                pv(205) = pv(205)+d2(i0,j0)*d(i0,k0)*d(i0,k1)/mu
                pv(211) = pv(211)+d(i0,j0)*d2(i0,k0)*d(i0,k1)/mu
                pv(212) = pv(212)+d3(i0,k0)*d(i0,k1)/mu
                pv(220) = pv(220)+d3(i0,j0)*d(j0,k0)/mu
                pv(251) = pv(251)+d2(i0,j0)*d(i0,k0)*d(j0,k0)/mu
                pv(276) = pv(276)+d(i0,j0)*d2(i0,k0)*d(j0,k0)/mu
                pv(280) = pv(280)+d3(i0,k0)*d(j0,k0)/mu
                pv(284) = pv(284)+d2(i0,j0)*d(i0,k1)*d(j0,k0)/mu
                pv(290) = pv(290)+d(i0,j0)*d(i0,k0)*d(i0,k1)*d(j0,k0)/mu
                pv(291) = pv(291)+d2(i0,k0)*d(i0,k1)*d(j0,k0)/mu
                pv(293) = pv(293)+d(i0,j0)*d2(i0,k1)*d(j0,k0)/mu
                pv(295) = pv(295)+d2(i0,j0)*d2(j0,k0)/mu
                pv(303) = pv(303)+d(i0,j0)*d(i0,k0)*d2(j0,k0)/mu
                pv(307) = pv(307)+d2(i0,k0)*d2(j0,k0)/mu
                pv(308) = pv(308)+d(i0,j0)*d(i0,k1)*d2(j0,k0)/mu
                pv(309) = pv(309)+d(i0,j0)*d3(j0,k0)/mu
                pv(311) = pv(311)+d(i0,k0)*d3(j0,k0)/mu
                pv(336) = pv(336)+d2(i0,j0)*d(j0,k0)*d(j0,k1)/mu
                pv(343) = pv(343)+d(i0,j0)*d(i0,k0)*d(j0,k0)*d(j0,k1)/mu
                pv(344) = pv(344)+d(i0,j0)*d2(j0,k0)*d(j0,k1)/mu
                pv(346) = pv(346)+d(i0,k0)*d2(j0,k0)*d(j0,k1)/mu
                pv(347) = pv(347)+d3(j0,k0)*d(j0,k1)/mu
            endif
        end subroutine mg442_i0j0

        !> @brief Accumulate secondary invariants for the paired-H and single-O index class.
        subroutine mg442_i1j0( mxd, mu )
            integer, intent (in) :: mxd, mu
            if( 2 .le. mxd )then
                pv(0) = pv(0)+d(i0,i1)*d(i0,j0)/mu
                pv(2) = pv(2)+d(i0,i1)*d(i0,k0)/mu
            endif
            if( 3 .le. mxd )then
                pv(8)  = pv(8)+d2(i0,i1)*d(i0,j0)/mu
                pv(10) = pv(10)+d(i0,i1)*d2(i0,j0)/mu
                pv(11) = pv(11)+d(i0,i1)*d(i0,j0)*d(i1,j0)/mu
                pv(25) = pv(25)+d2(i0,i1)*d(i0,k0)/mu
                pv(27) = pv(27)+d(i0,i1)*d(i0,j0)*d(i0,k0)/mu
                pv(29) = pv(29)+d(i0,i1)*d(i1,j0)*d(i0,k0)/mu
                pv(30) = pv(30)+d(i0,j0)*d(i1,j0)*d(i0,k0)/mu
                pv(33) = pv(33)+d(i0,i1)*d2(i0,k0)/mu
                pv(36) = pv(36)+d(i0,i1)*d(i0,k0)*d(i1,k0)/mu
                pv(38) = pv(38)+d(i0,j0)*d(i0,k0)*d(i1,k0)/mu
                pv(39) = pv(39)+d2(i0,k0)*d(i1,k0)/mu
                pv(40) = pv(40)+d(i0,i1)*d(i0,k0)*d(i0,k1)/mu
                pv(42) = pv(42)+d(i0,i1)*d(i1,k0)*d(i0,k1)/mu
                pv(43) = pv(43)+d(i0,i1)*d(i0,j0)*d(j0,k0)/mu
                pv(45) = pv(45)+d(i0,j0)*d(i1,j0)*d(j0,k0)/mu
                pv(51) = pv(51)+d(i0,i1)*d(i0,k0)*d(j0,k0)/mu
                pv(53) = pv(53)+d(i1,j0)*d(i0,k0)*d(j0,k0)/mu
            endif
            if( 4 .le. mxd )then
                pv(69)  = pv(69)+d4(i0,i1)/mu
                pv(70)  = pv(70)+d3(i0,i1)*d(i0,j0)/mu
                pv(73)  = pv(73)+d2(i0,i1)*d2(i0,j0)/mu
                pv(75)  = pv(75)+d(i0,i1)*d3(i0,j0)/mu
                pv(76)  = pv(76)+d2(i0,i1)*d(i0,j0)*d(i1,j0)/mu
                pv(80)  = pv(80)+d(i0,i1)*d2(i0,j0)*d(i1,j0)/mu
                pv(83)  = pv(83)+d2(i0,j0)*d2(i1,j0)/mu
                pv(141) = pv(141)+d3(i0,i1)*d(i0,k0)/mu
                pv(144) = pv(144)+d2(i0,i1)*d(i0,j0)*d(i0,k0)/mu
                pv(146) = pv(146)+d(i0,i1)*d2(i0,j0)*d(i0,k0)/mu
                pv(148) = pv(148)+d2(i0,i1)*d(i1,j0)*d(i0,k0)/mu
                pv(150) = pv(150)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,k0)/mu
                pv(153) = pv(153)+d2(i0,j0)*d(i1,j0)*d(i0,k0)/mu
                pv(154) = pv(154)+d(i0,i1)*d2(i1,j0)*d(i0,k0)/mu
                pv(155) = pv(155)+d(i0,j0)*d2(i1,j0)*d(i0,k0)/mu
                pv(172) = pv(172)+d2(i0,i1)*d2(i0,k0)/mu
                pv(174) = pv(174)+d(i0,i1)*d(i0,j0)*d2(i0,k0)/mu
                pv(176) = pv(176)+d(i0,i1)*d(i1,j0)*d2(i0,k0)/mu
                pv(177) = pv(177)+d(i0,j0)*d(i1,j0)*d2(i0,k0)/mu
                pv(180) = pv(180)+d(i0,i1)*d3(i0,k0)/mu
                pv(182) = pv(182)+d2(i0,i1)*d(i0,k0)*d(i1,k0)/mu
                pv(186) = pv(186)+d(i0,i1)*d(i0,j0)*d(i0,k0)*d(i1,k0)/mu
                pv(189) = pv(189)+d2(i0,j0)*d(i0,k0)*d(i1,k0)/mu
                pv(190) = pv(190)+d(i0,j0)*d(i1,j0)*d(i0,k0)*d(i1,k0)/mu
                pv(195) = pv(195)+d(i0,i1)*d2(i0,k0)*d(i1,k0)/mu
                pv(198) = pv(198)+d(i0,j0)*d2(i0,k0)*d(i1,k0)/mu
                pv(199) = pv(199)+d(i1,j0)*d2(i0,k0)*d(i1,k0)/mu
                pv(200) = pv(200)+d3(i0,k0)*d(i1,k0)/mu
                pv(201) = pv(201)+d2(i0,k0)*d2(i1,k0)/mu
                pv(202) = pv(202)+d2(i0,i1)*d(i0,k0)*d(i0,k1)/mu
                pv(204) = pv(204)+d(i0,i1)*d(i0,j0)*d(i0,k0)*d(i0,k1)/mu
                pv(206) = pv(206)+d(i0,i1)*d(i1,j0)*d(i0,k0)*d(i0,k1)/mu
                pv(207) = pv(207)+d(i0,j0)*d(i1,j0)*d(i0,k0)*d(i0,k1)/mu
                pv(210) = pv(210)+d(i0,i1)*d2(i0,k0)*d(i0,k1)/mu
                pv(213) = pv(213)+d2(i0,i1)*d(i1,k0)*d(i0,k1)/mu
                pv(214) = pv(214)+d(i0,j0)*d(i1,j0)*d(i1,k0)*d(i0,k1)/mu
                pv(215) = pv(215)+d(i0,i1)*d(i0,k0)*d(i1,k0)*d(i0,k1)/mu
                pv(216) = pv(216)+d2(i0,k0)*d(i1,k0)*d(i0,k1)/mu
                pv(217) = pv(217)+d2(i0,i1)*d(i0,j0)*d(j0,k0)/mu
                pv(219) = pv(219)+d(i0,i1)*d2(i0,j0)*d(j0,k0)/mu
                pv(221) = pv(221)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(j0,k0)/mu
                pv(223) = pv(223)+d2(i0,j0)*d(i1,j0)*d(j0,k0)/mu
                pv(248) = pv(248)+d2(i0,i1)*d(i0,k0)*d(j0,k0)/mu
                pv(250) = pv(250)+d(i0,i1)*d(i0,j0)*d(i0,k0)*d(j0,k0)/mu
                pv(252) = pv(252)+d(i0,i1)*d(i1,j0)*d(i0,k0)*d(j0,k0)/mu
                pv(255) = pv(255)+d(i0,j0)*d(i1,j0)*d(i0,k0)*d(j0,k0)/mu
                pv(256) = pv(256)+d2(i1,j0)*d(i0,k0)*d(j0,k0)/mu
                pv(275) = pv(275)+d(i0,i1)*d2(i0,k0)*d(j0,k0)/mu
                pv(277) = pv(277)+d(i1,j0)*d2(i0,k0)*d(j0,k0)/mu
                pv(281) = pv(281)+d(i0,i1)*d(i0,k0)*d(i1,k0)*d(j0,k0)/mu
                pv(282) = pv(282)+d(i0,j0)*d(i0,k0)*d(i1,k0)*d(j0,k0)/mu
                pv(283) = pv(283)+d(i0,i1)*d(i0,j0)*d(i0,k1)*d(j0,k0)/mu
                pv(285) = pv(285)+d(i0,i1)*d(i1,j0)*d(i0,k1)*d(j0,k0)/mu
                pv(286) = pv(286)+d(i0,j0)*d(i1,j0)*d(i0,k1)*d(j0,k0)/mu
                pv(292) = pv(292)+d(i0,j0)*d(i1,k0)*d(i0,k1)*d(j0,k0)/mu
                pv(294) = pv(294)+d(i0,i1)*d(i0,j0)*d2(j0,k0)/mu
                pv(296) = pv(296)+d(i0,j0)*d(i1,j0)*d2(j0,k0)/mu
                pv(302) = pv(302)+d(i0,i1)*d(i0,k0)*d2(j0,k0)/mu
                pv(304) = pv(304)+d(i1,j0)*d(i0,k0)*d2(j0,k0)/mu
                pv(335) = pv(335)+d(i0,i1)*d(i0,j0)*d(j0,k0)*d(j0,k1)/mu
                pv(337) = pv(337)+d(i0,j0)*d(i1,j0)*d(j0,k0)*d(j0,k1)/mu
            endif
        end subroutine mg442_i1j0

        !> @brief Accumulate secondary invariants for the multi-H and single-O index class.
        subroutine mg442_iij0( mxd, mu )
            integer, intent (in) :: mxd, mu
            if( 3 .le. mxd )then
                pv(7)  = pv(7)+d2(i0,i1)*d(i0,i2)/mu
                pv(9)  = pv(9)+d(i0,i1)*d(i0,i2)*d(i0,j0)/mu
                pv(12) = pv(12)+d(i0,i2)*d(i0,j0)*d(i1,j0)/mu
                pv(26) = pv(26)+d(i0,i1)*d(i0,i2)*d(i0,k0)/mu
                pv(37) = pv(37)+d(i0,i2)*d(i0,k0)*d(i1,k0)/mu
            endif
            if( 4 .le. mxd )then
                pv(71)  = pv(71)+d2(i0,i1)*d(i0,i2)*d(i0,j0)/mu
                pv(72)  = pv(72)+d2(i0,i1)*d(i1,i2)*d(i0,j0)/mu
                pv(74)  = pv(74)+d(i0,i1)*d(i0,i2)*d2(i0,j0)/mu
                pv(77)  = pv(77)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i1,j0)/mu
                pv(78)  = pv(78)+d2(i0,i2)*d(i0,j0)*d(i1,j0)/mu
                pv(79)  = pv(79)+d(i0,i2)*d(i1,i2)*d(i0,j0)*d(i1,j0)/mu
                pv(81)  = pv(81)+d(i0,i2)*d2(i0,j0)*d(i1,j0)/mu
                pv(82)  = pv(82)+d(i1,i2)*d2(i0,j0)*d(i1,j0)/mu
                pv(84)  = pv(84)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i2,j0)/mu
                pv(85)  = pv(85)+d2(i0,j0)*d(i1,j0)*d(i2,j0)/mu
                pv(142) = pv(142)+d2(i0,i1)*d(i0,i2)*d(i0,k0)/mu
                pv(143) = pv(143)+d2(i0,i1)*d(i1,i2)*d(i0,k0)/mu
                pv(145) = pv(145)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,k0)/mu
                pv(149) = pv(149)+d(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,k0)/mu
                pv(151) = pv(151)+d(i0,i2)*d(i0,j0)*d(i1,j0)*d(i0,k0)/mu
                pv(152) = pv(152)+d(i1,i2)*d(i0,j0)*d(i1,j0)*d(i0,k0)/mu
                pv(156) = pv(156)+d(i0,j0)*d(i1,j0)*d(i2,j0)*d(i0,k0)/mu
                pv(173) = pv(173)+d(i0,i1)*d(i0,i2)*d2(i0,k0)/mu
                pv(183) = pv(183)+d(i0,i1)*d(i0,i2)*d(i0,k0)*d(i1,k0)/mu
                pv(184) = pv(184)+d2(i0,i2)*d(i0,k0)*d(i1,k0)/mu
                pv(185) = pv(185)+d(i0,i2)*d(i1,i2)*d(i0,k0)*d(i1,k0)/mu
                pv(187) = pv(187)+d(i0,i2)*d(i0,j0)*d(i0,k0)*d(i1,k0)/mu
                pv(188) = pv(188)+d(i1,i2)*d(i0,j0)*d(i0,k0)*d(i1,k0)/mu
                pv(191) = pv(191)+d(i0,j0)*d(i2,j0)*d(i0,k0)*d(i1,k0)/mu
                pv(196) = pv(196)+d(i0,i2)*d2(i0,k0)*d(i1,k0)/mu
                pv(197) = pv(197)+d(i1,i2)*d2(i0,k0)*d(i1,k0)/mu
                pv(203) = pv(203)+d(i0,i1)*d(i0,i2)*d(i0,k0)*d(i0,k1)/mu
                pv(218) = pv(218)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(j0,k0)/mu
                pv(222) = pv(222)+d(i0,i2)*d(i0,j0)*d(i1,j0)*d(j0,k0)/mu
                pv(224) = pv(224)+d(i0,j0)*d(i1,j0)*d(i2,j0)*d(j0,k0)/mu
                pv(249) = pv(249)+d(i0,i1)*d(i0,i2)*d(i0,k0)*d(j0,k0)/mu
                pv(253) = pv(253)+d(i0,i2)*d(i1,j0)*d(i0,k0)*d(j0,k0)/mu
                pv(254) = pv(254)+d(i1,i2)*d(i1,j0)*d(i0,k0)*d(j0,k0)/mu
                pv(257) = pv(257)+d(i1,j0)*d(i2,j0)*d(i0,k0)*d(j0,k0)/mu
            endif
        end subroutine mg442_iij0

        !> @brief Accumulate secondary invariants for the single-H and paired-O index class.
        subroutine mg442_i0j1( mxd, mu )
            integer, intent (in) :: mxd, mu
            if( 2 .le. mxd )then
                pv(1) = pv(1)+d(i0,j0)*d(j0,j1)/mu
                pv(5) = pv(5)+d(j0,j1)*d(j0,k0)/mu
            endif
            if( 3 .le. mxd )then
                pv(17) = pv(17)+d2(i0,j0)*d(j0,j1)/mu
                pv(19) = pv(19)+d(i0,j0)*d(i0,j1)*d(j0,j1)/mu
                pv(22) = pv(22)+d(i0,j0)*d2(j0,j1)/mu
                pv(31) = pv(31)+d(i0,j0)*d(i0,j1)*d(i0,k0)/mu
                pv(32) = pv(32)+d(i0,j0)*d(j0,j1)*d(i0,k0)/mu
                pv(46) = pv(46)+d(i0,j0)*d(i0,j1)*d(j0,k0)/mu
                pv(47) = pv(47)+d(i0,j0)*d(j0,j1)*d(j0,k0)/mu
                pv(48) = pv(48)+d(i0,j1)*d(j0,j1)*d(j0,k0)/mu
                pv(49) = pv(49)+d2(j0,j1)*d(j0,k0)/mu
                pv(54) = pv(54)+d(i0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(55) = pv(55)+d(j0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(59) = pv(59)+d(j0,j1)*d2(j0,k0)/mu
                pv(62) = pv(62)+d(i0,j0)*d(j0,k0)*d(j1,k0)/mu
                pv(63) = pv(63)+d(j0,j1)*d(j0,k0)*d(j1,k0)/mu
                pv(65) = pv(65)+d2(j0,k0)*d(j1,k0)/mu
                pv(67) = pv(67)+d(j0,j1)*d(j0,k0)*d(j0,k1)/mu
                pv(68) = pv(68)+d(j0,j1)*d(j1,k0)*d(j0,k1)/mu
            endif
            if( 4 .le. mxd )then
                pv(97)  = pv(97)+d2(i0,j0)*d2(i0,j1)/mu
                pv(108) = pv(108)+d3(i0,j0)*d(j0,j1)/mu
                pv(114) = pv(114)+d2(i0,j0)*d(i0,j1)*d(j0,j1)/mu
                pv(127) = pv(127)+d2(i0,j0)*d2(j0,j1)/mu
                pv(129) = pv(129)+d(i0,j0)*d(i0,j1)*d2(j0,j1)/mu
                pv(132) = pv(132)+d(i0,j0)*d3(j0,j1)/mu
                pv(133) = pv(133)+d4(j0,j1)/mu
                pv(158) = pv(158)+d2(i0,j0)*d(i0,j1)*d(i0,k0)/mu
                pv(164) = pv(164)+d2(i0,j0)*d(j0,j1)*d(i0,k0)/mu
                pv(167) = pv(167)+d(i0,j0)*d(i0,j1)*d(j0,j1)*d(i0,k0)/mu
                pv(170) = pv(170)+d(i0,j0)*d2(j0,j1)*d(i0,k0)/mu
                pv(178) = pv(178)+d(i0,j0)*d(i0,j1)*d2(i0,k0)/mu
                pv(179) = pv(179)+d(i0,j0)*d(j0,j1)*d2(i0,k0)/mu
                pv(208) = pv(208)+d(i0,j0)*d(i0,j1)*d(i0,k0)*d(i0,k1)/mu
                pv(209) = pv(209)+d(i0,j0)*d(j0,j1)*d(i0,k0)*d(i0,k1)/mu
                pv(226) = pv(226)+d2(i0,j0)*d(i0,j1)*d(j0,k0)/mu
                pv(229) = pv(229)+d(i0,j0)*d2(i0,j1)*d(j0,k0)/mu
                pv(233) = pv(233)+d2(i0,j0)*d(j0,j1)*d(j0,k0)/mu
                pv(236) = pv(236)+d(i0,j0)*d(i0,j1)*d(j0,j1)*d(j0,k0)/mu
                pv(238) = pv(238)+d2(i0,j1)*d(j0,j1)*d(j0,k0)/mu
                pv(241) = pv(241)+d(i0,j0)*d2(j0,j1)*d(j0,k0)/mu
                pv(242) = pv(242)+d(i0,j1)*d2(j0,j1)*d(j0,k0)/mu
                pv(243) = pv(243)+d3(j0,j1)*d(j0,k0)/mu
                pv(259) = pv(259)+d(i0,j0)*d(i0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(261) = pv(261)+d2(i0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(267) = pv(267)+d(i0,j0)*d(j0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(269) = pv(269)+d(i0,j1)*d(j0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(272) = pv(272)+d2(j0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(278) = pv(278)+d(i0,j1)*d2(i0,k0)*d(j0,k0)/mu
                pv(279) = pv(279)+d(j0,j1)*d2(i0,k0)*d(j0,k0)/mu
                pv(287) = pv(287)+d(i0,j0)*d(i0,j1)*d(i0,k1)*d(j0,k0)/mu
                pv(288) = pv(288)+d(i0,j0)*d(j0,j1)*d(i0,k1)*d(j0,k0)/mu
                pv(289) = pv(289)+d(i0,j1)*d(j0,j1)*d(i0,k1)*d(j0,k0)/mu
                pv(297) = pv(297)+d(i0,j0)*d(i0,j1)*d2(j0,k0)/mu
                pv(298) = pv(298)+d(i0,j0)*d(j0,j1)*d2(j0,k0)/mu
                pv(299) = pv(299)+d(i0,j1)*d(j0,j1)*d2(j0,k0)/mu
                pv(300) = pv(300)+d2(j0,j1)*d2(j0,k0)/mu
                pv(305) = pv(305)+d(i0,j1)*d(i0,k0)*d2(j0,k0)/mu
                pv(306) = pv(306)+d(j0,j1)*d(i0,k0)*d2(j0,k0)/mu
                pv(310) = pv(310)+d(j0,j1)*d3(j0,k0)/mu
                pv(313) = pv(313)+d2(i0,j0)*d(j0,k0)*d(j1,k0)/mu
                pv(315) = pv(315)+d(i0,j0)*d(i0,j1)*d(j0,k0)*d(j1,k0)/mu
                pv(318) = pv(318)+d(i0,j0)*d(j0,j1)*d(j0,k0)*d(j1,k0)/mu
                pv(319) = pv(319)+d2(j0,j1)*d(j0,k0)*d(j1,k0)/mu
                pv(325) = pv(325)+d(i0,j0)*d(i0,k0)*d(j0,k0)*d(j1,k0)/mu
                pv(326) = pv(326)+d(j0,j1)*d(i0,k0)*d(j0,k0)*d(j1,k0)/mu
                pv(327) = pv(327)+d(i0,j0)*d(i0,k1)*d(j0,k0)*d(j1,k0)/mu
                pv(328) = pv(328)+d(i0,j0)*d2(j0,k0)*d(j1,k0)/mu
                pv(329) = pv(329)+d(i0,j1)*d2(j0,k0)*d(j1,k0)/mu
                pv(330) = pv(330)+d(j0,j1)*d2(j0,k0)*d(j1,k0)/mu
                pv(333) = pv(333)+d3(j0,k0)*d(j1,k0)/mu
                pv(334) = pv(334)+d2(j0,k0)*d2(j1,k0)/mu
                pv(338) = pv(338)+d(i0,j0)*d(i0,j1)*d(j0,k0)*d(j0,k1)/mu
                pv(339) = pv(339)+d(i0,j0)*d(j0,j1)*d(j0,k0)*d(j0,k1)/mu
                pv(340) = pv(340)+d(i0,j1)*d(j0,j1)*d(j0,k0)*d(j0,k1)/mu
                pv(341) = pv(341)+d2(j0,j1)*d(j0,k0)*d(j0,k1)/mu
                pv(345) = pv(345)+d(j0,j1)*d2(j0,k0)*d(j0,k1)/mu
                pv(348) = pv(348)+d(i0,j0)*d(i0,j1)*d(j1,k0)*d(j0,k1)/mu
                pv(349) = pv(349)+d2(j0,j1)*d(j1,k0)*d(j0,k1)/mu
                pv(350) = pv(350)+d(j0,j1)*d(j0,k0)*d(j1,k0)*d(j0,k1)/mu
                pv(351) = pv(351)+d2(j0,k0)*d(j1,k0)*d(j0,k1)/mu
            endif
        end subroutine mg442_i0j1

        !> @brief Accumulate secondary invariants for the paired-H and paired-O index class.
        subroutine mg442_i1j1( mxd, mu )
            integer, intent (in) :: mxd, mu
            if( 3 .le. mxd )then
                pv(13) = pv(13)+d(i0,i1)*d(i0,j0)*d(i0,j1)/mu
                pv(14) = pv(14)+d(i0,i1)*d(i1,j0)*d(i0,j1)/mu
                pv(15) = pv(15)+d(i0,j0)*d(i1,j0)*d(i0,j1)/mu
                pv(16) = pv(16)+d(i0,i1)*d(i0,j0)*d(j0,j1)/mu
                pv(18) = pv(18)+d(i0,j0)*d(i1,j0)*d(j0,j1)/mu
                pv(20) = pv(20)+d(i1,j0)*d(i0,j1)*d(j0,j1)/mu
            endif
            if( 4 .le. mxd )then
                pv(86)  = pv(86)+d2(i0,i1)*d(i0,j0)*d(i0,j1)/mu
                pv(88)  = pv(88)+d(i0,i1)*d2(i0,j0)*d(i0,j1)/mu
                pv(89)  = pv(89)+d2(i0,i1)*d(i1,j0)*d(i0,j1)/mu
                pv(90)  = pv(90)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,j1)/mu
                pv(93)  = pv(93)+d2(i0,j0)*d(i1,j0)*d(i0,j1)/mu
                pv(94)  = pv(94)+d(i0,i1)*d2(i1,j0)*d(i0,j1)/mu
                pv(95)  = pv(95)+d(i0,j0)*d2(i1,j0)*d(i0,j1)/mu
                pv(98)  = pv(98)+d(i0,j0)*d(i1,j0)*d2(i0,j1)/mu
                pv(99)  = pv(99)+d(i0,j0)*d(i1,j0)*d(i0,j1)*d(i1,j1)/mu
                pv(105) = pv(105)+d2(i0,i1)*d(i0,j0)*d(j0,j1)/mu
                pv(107) = pv(107)+d(i0,i1)*d2(i0,j0)*d(j0,j1)/mu
                pv(109) = pv(109)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(j0,j1)/mu
                pv(111) = pv(111)+d2(i0,j0)*d(i1,j0)*d(j0,j1)/mu
                pv(113) = pv(113)+d(i0,i1)*d(i0,j0)*d(i0,j1)*d(j0,j1)/mu
                pv(115) = pv(115)+d(i0,i1)*d(i1,j0)*d(i0,j1)*d(j0,j1)/mu
                pv(117) = pv(117)+d(i0,j0)*d(i1,j0)*d(i0,j1)*d(j0,j1)/mu
                pv(118) = pv(118)+d2(i1,j0)*d(i0,j1)*d(j0,j1)/mu
                pv(126) = pv(126)+d(i0,i1)*d(i0,j0)*d2(j0,j1)/mu
                pv(128) = pv(128)+d(i0,j0)*d(i1,j0)*d2(j0,j1)/mu
                pv(130) = pv(130)+d(i1,j0)*d(i0,j1)*d2(j0,j1)/mu
                pv(157) = pv(157)+d(i0,i1)*d(i0,j0)*d(i0,j1)*d(i0,k0)/mu
                pv(159) = pv(159)+d(i0,i1)*d(i1,j0)*d(i0,j1)*d(i0,k0)/mu
                pv(160) = pv(160)+d(i0,j0)*d(i1,j0)*d(i0,j1)*d(i0,k0)/mu
                pv(161) = pv(161)+d(i0,j0)*d(i1,j0)*d(i1,j1)*d(i0,k0)/mu
                pv(163) = pv(163)+d(i0,i1)*d(i0,j0)*d(j0,j1)*d(i0,k0)/mu
                pv(165) = pv(165)+d(i0,i1)*d(i1,j0)*d(j0,j1)*d(i0,k0)/mu
                pv(166) = pv(166)+d(i0,j0)*d(i1,j0)*d(j0,j1)*d(i0,k0)/mu
                pv(168) = pv(168)+d(i1,j0)*d(i0,j1)*d(j0,j1)*d(i0,k0)/mu
                pv(192) = pv(192)+d(i0,j0)*d(i0,j1)*d(i0,k0)*d(i1,k0)/mu
                pv(193) = pv(193)+d(i1,j0)*d(i0,j1)*d(i0,k0)*d(i1,k0)/mu
                pv(194) = pv(194)+d(i0,j0)*d(j0,j1)*d(i0,k0)*d(i1,k0)/mu
                pv(225) = pv(225)+d(i0,i1)*d(i0,j0)*d(i0,j1)*d(j0,k0)/mu
                pv(227) = pv(227)+d(i0,i1)*d(i1,j0)*d(i0,j1)*d(j0,k0)/mu
                pv(228) = pv(228)+d(i0,j0)*d(i1,j0)*d(i0,j1)*d(j0,k0)/mu
                pv(230) = pv(230)+d(i0,j0)*d(i0,j1)*d(i1,j1)*d(j0,k0)/mu
                pv(232) = pv(232)+d(i0,i1)*d(i0,j0)*d(j0,j1)*d(j0,k0)/mu
                pv(234) = pv(234)+d(i0,j0)*d(i1,j0)*d(j0,j1)*d(j0,k0)/mu
                pv(235) = pv(235)+d(i0,i1)*d(i0,j1)*d(j0,j1)*d(j0,k0)/mu
                pv(237) = pv(237)+d(i1,j0)*d(i0,j1)*d(j0,j1)*d(j0,k0)/mu
                pv(258) = pv(258)+d(i0,i1)*d(i0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(260) = pv(260)+d(i1,j0)*d(i0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(262) = pv(262)+d(i0,i1)*d(i1,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(263) = pv(263)+d(i1,j0)*d(i1,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(264) = pv(264)+d(i0,j1)*d(i1,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(266) = pv(266)+d(i0,i1)*d(j0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(268) = pv(268)+d(i1,j0)*d(j0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(270) = pv(270)+d(i1,j1)*d(j0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(312) = pv(312)+d(i0,i1)*d(i0,j0)*d(j0,k0)*d(j1,k0)/mu
                pv(314) = pv(314)+d(i0,j0)*d(i1,j0)*d(j0,k0)*d(j1,k0)/mu
                pv(316) = pv(316)+d(i1,j0)*d(i0,j1)*d(j0,k0)*d(j1,k0)/mu
            endif
        end subroutine mg442_i1j1

        !> @brief Accumulate secondary invariants for the multi-H and paired-O index class.
        subroutine mg442_iij1( mxd, mu )
            integer, intent (in) :: mxd, mu
            if( 4 .le. mxd )then
                pv(87)  = pv(87)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,j1)/mu
                pv(91)  = pv(91)+d(i0,i2)*d(i0,j0)*d(i1,j0)*d(i0,j1)/mu
                pv(92)  = pv(92)+d(i1,i2)*d(i0,j0)*d(i1,j0)*d(i0,j1)/mu
                pv(96)  = pv(96)+d(i0,j0)*d(i1,j0)*d(i2,j0)*d(i0,j1)/mu
                pv(100) = pv(100)+d(i0,j0)*d(i2,j0)*d(i0,j1)*d(i1,j1)/mu
                pv(106) = pv(106)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(j0,j1)/mu
                pv(110) = pv(110)+d(i0,i2)*d(i0,j0)*d(i1,j0)*d(j0,j1)/mu
                pv(112) = pv(112)+d(i0,j0)*d(i1,j0)*d(i2,j0)*d(j0,j1)/mu
                pv(116) = pv(116)+d(i0,i2)*d(i1,j0)*d(i0,j1)*d(j0,j1)/mu
            endif
        end subroutine mg442_iij1

        !> @brief Accumulate secondary invariants for the single-H and multi-O index class.
        subroutine mg442_i0jj( mxd, mu )
            integer, intent (in) :: mxd, mu
            if( 3 .le. mxd )then
                pv(21) = pv(21)+d(i0,j0)*d(i0,j2)*d(j0,j1)/mu
                pv(23) = pv(23)+d(i0,j0)*d(j0,j1)*d(j0,j2)/mu
                pv(24) = pv(24)+d2(j0,j1)*d(j0,j2)/mu
                pv(50) = pv(50)+d(j0,j1)*d(j0,j2)*d(j0,k0)/mu
                pv(64) = pv(64)+d(j0,j2)*d(j0,k0)*d(j1,k0)/mu
            endif
            if( 4 .le. mxd )then
                pv(102) = pv(102)+d2(i0,j0)*d(i0,j1)*d(i0,j2)/mu
                pv(120) = pv(120)+d2(i0,j0)*d(i0,j2)*d(j0,j1)/mu
                pv(123) = pv(123)+d(i0,j0)*d(i0,j1)*d(i0,j2)*d(j0,j1)/mu
                pv(125) = pv(125)+d(i0,j0)*d2(i0,j2)*d(j0,j1)/mu
                pv(131) = pv(131)+d(i0,j0)*d(i0,j2)*d2(j0,j1)/mu
                pv(135) = pv(135)+d2(i0,j0)*d(j0,j1)*d(j0,j2)/mu
                pv(137) = pv(137)+d(i0,j0)*d(i0,j1)*d(j0,j1)*d(j0,j2)/mu
                pv(138) = pv(138)+d(i0,j1)*d(i0,j2)*d(j0,j1)*d(j0,j2)/mu
                pv(139) = pv(139)+d(i0,j0)*d2(j0,j1)*d(j0,j2)/mu
                pv(140) = pv(140)+d(i0,j1)*d2(j0,j1)*d(j0,j2)/mu
                pv(162) = pv(162)+d(i0,j0)*d(i0,j1)*d(i0,j2)*d(i0,k0)/mu
                pv(169) = pv(169)+d(i0,j0)*d(i0,j2)*d(j0,j1)*d(i0,k0)/mu
                pv(171) = pv(171)+d(i0,j0)*d(j0,j1)*d(j0,j2)*d(i0,k0)/mu
                pv(231) = pv(231)+d(i0,j0)*d(i0,j1)*d(i0,j2)*d(j0,k0)/mu
                pv(239) = pv(239)+d(i0,j0)*d(i0,j2)*d(j0,j1)*d(j0,k0)/mu
                pv(240) = pv(240)+d(i0,j1)*d(i0,j2)*d(j0,j1)*d(j0,k0)/mu
                pv(244) = pv(244)+d(i0,j0)*d(j0,j1)*d(j0,j2)*d(j0,k0)/mu
                pv(245) = pv(245)+d(i0,j1)*d(j0,j1)*d(j0,j2)*d(j0,k0)/mu
                pv(246) = pv(246)+d2(j0,j1)*d(j0,j2)*d(j0,k0)/mu
                pv(247) = pv(247)+d2(j0,j1)*d(j1,j2)*d(j0,k0)/mu
                pv(265) = pv(265)+d(i0,j1)*d(i0,j2)*d(i0,k0)*d(j0,k0)/mu
                pv(271) = pv(271)+d(i0,j2)*d(j0,j1)*d(i0,k0)*d(j0,k0)/mu
                pv(273) = pv(273)+d(j0,j1)*d(j0,j2)*d(i0,k0)*d(j0,k0)/mu
                pv(274) = pv(274)+d(i0,j1)*d(j1,j2)*d(i0,k0)*d(j0,k0)/mu
                pv(301) = pv(301)+d(j0,j1)*d(j0,j2)*d2(j0,k0)/mu
                pv(317) = pv(317)+d(i0,j0)*d(i0,j2)*d(j0,k0)*d(j1,k0)/mu
                pv(320) = pv(320)+d(i0,j0)*d(j0,j2)*d(j0,k0)*d(j1,k0)/mu
                pv(321) = pv(321)+d(i0,j1)*d(j0,j2)*d(j0,k0)*d(j1,k0)/mu
                pv(322) = pv(322)+d(j0,j1)*d(j0,j2)*d(j0,k0)*d(j1,k0)/mu
                pv(323) = pv(323)+d2(j0,j2)*d(j0,k0)*d(j1,k0)/mu
                pv(324) = pv(324)+d(j0,j2)*d(j1,j2)*d(j0,k0)*d(j1,k0)/mu
                pv(331) = pv(331)+d(j0,j2)*d2(j0,k0)*d(j1,k0)/mu
                pv(332) = pv(332)+d(j1,j2)*d2(j0,k0)*d(j1,k0)/mu
                pv(342) = pv(342)+d(j0,j1)*d(j0,j2)*d(j0,k0)*d(j0,k1)/mu
            endif
        end subroutine mg442_i0jj

        !> @brief Accumulate secondary invariants for the paired-H and multi-O index class.
        subroutine mg442_i1jj( mxd, mu )
            integer, intent (in) :: mxd, mu
            if( 4 .le. mxd )then
                pv(101) = pv(101)+d(i0,i1)*d(i0,j0)*d(i0,j1)*d(i0,j2)/mu
                pv(103) = pv(103)+d(i0,j0)*d(i1,j0)*d(i0,j1)*d(i0,j2)/mu
                pv(104) = pv(104)+d(i0,j0)*d(i1,j0)*d(i1,j1)*d(i0,j2)/mu
                pv(119) = pv(119)+d(i0,i1)*d(i0,j0)*d(i0,j2)*d(j0,j1)/mu
                pv(121) = pv(121)+d(i0,i1)*d(i1,j0)*d(i0,j2)*d(j0,j1)/mu
                pv(122) = pv(122)+d(i0,j0)*d(i1,j0)*d(i0,j2)*d(j0,j1)/mu
                pv(124) = pv(124)+d(i1,j0)*d(i0,j1)*d(i0,j2)*d(j0,j1)/mu
                pv(134) = pv(134)+d(i0,i1)*d(i0,j0)*d(j0,j1)*d(j0,j2)/mu
                pv(136) = pv(136)+d(i0,j0)*d(i1,j0)*d(j0,j1)*d(j0,j2)/mu
            endif
        end subroutine mg442_i1jj
    end subroutine mg442_isecs

    !> @brief Assemble secondary invariant products through polynomial degree four.
    subroutine mg442_secs( mxd, r, v )
        integer, intent (in)   :: mxd
        real(wp), intent (in)  :: r(0:,0:)
        real(wp), intent (out) :: v(0:)
        real(wp)               :: pv(0:size(v)-1) !! an upper bound size
        call mg442_isecs(mxd, r, pv)
        v(0) = 1
        if( 2 .le. mxd )then
            v(1:7) = pv(0:6)
        endif
        if( 3 .le. mxd )then
            v(8:69) = pv(7:68)
        endif
        if( 4 .le. mxd )then
            v(70)     = pv(0)*pv(0)
            v(71)     = pv(0)*pv(1)
            v(72)     = pv(1)*pv(1)
            v(73)     = pv(0)*pv(2)
            v(74)     = pv(1)*pv(2)
            v(75)     = pv(2)*pv(2)
            v(76)     = pv(0)*pv(3)
            v(77)     = pv(1)*pv(3)
            v(78)     = pv(2)*pv(3)
            v(79)     = pv(3)*pv(3)
            v(80)     = pv(0)*pv(4)
            v(81)     = pv(1)*pv(4)
            v(82)     = pv(2)*pv(4)
            v(83)     = pv(3)*pv(4)
            v(84)     = pv(4)*pv(4)
            v(85)     = pv(0)*pv(5)
            v(86)     = pv(1)*pv(5)
            v(87)     = pv(2)*pv(5)
            v(88)     = pv(3)*pv(5)
            v(89)     = pv(4)*pv(5)
            v(90)     = pv(5)*pv(5)
            v(91)     = pv(0)*pv(6)
            v(92)     = pv(1)*pv(6)
            v(93)     = pv(2)*pv(6)
            v(94)     = pv(3)*pv(6)
            v(95)     = pv(4)*pv(6)
            v(96)     = pv(5)*pv(6)
            v(97:379) = pv(69:351)
        endif
    end subroutine mg442_secs

    !> @brief Assemble the complete 1,784-term four-four-two invariant basis.
    subroutine mg442_base( mxd, r, w )
        integer, intent (in)   :: mxd
        real(wp), intent (in)  :: r(0:,0:)
        real(wp), intent (out) :: w(0:mg442_nb(mxd)-1)
        integer                :: k, l0, m0, l, n, ind
        real(wp)               :: u(0:nr-1), v(0:mg442_nsc(mxd)-1), &
            w0(0:mg442_npb(mxd)-1)
        call mg442_prims(r, u)
        call mg442_prib(mxd, u, w0)
        call mg442_secs(mxd, r, v)
        ind = 0
        do n = 0, mxd
            do k = 0, n
                if( k .eq. 0 )then
                    l0 = 0
                else
                    l0 = mg442_nsc(max(0, k-1))
                endif
                if( k .eq. n )then
                    m0 = 0
                else
                    m0 = mg442_npb(max(0, n-k-1))
                endif
                do l = l0, mg442_nsc(k)-1
                    w(ind:ind+mg442_dnpb(n-k)-1) = v(l)*w0(m0:mg442_npb(n-k)-1)
                    ind                          = ind+mg442_dnpb(n-k)
                enddo
            enddo
        enddo
    end subroutine mg442_base

    !> @brief Contract the immutable formic-acid dimer coefficients with the invariant basis.
    function evaluate_fad_fit( distances, coefficients ) result(energy_hartree)
        real(wp), intent(in) :: distances(nk,nk)
        real(wp), intent(in) :: coefficients(1784)
        real(wp)             :: energy_hartree
        real(wp)             :: transformed(0:nk-1,0:nk-1), basis(0:1783)
        integer              :: i, j
        transformed = 0.0_wp
        do j=1, nk
            do i=1, nk
                if( i /= j ) transformed(i-1,j-1) = exp(-distances(i,j)/2.0_wp)
            enddo
        enddo
        call mg442_base(4, transformed, basis)
        energy_hartree = dot_product(coefficients, basis)
    end function evaluate_fad_fit
end module chempotpy_c2h4o4_invariants
