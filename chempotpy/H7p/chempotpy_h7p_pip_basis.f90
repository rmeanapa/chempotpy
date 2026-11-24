!> @brief Evaluate the permutation-invariant polynomial basis used by H7p_PIP.
module chempotpy_h7p_pip_basis
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: h7p_invariant_basis
    integer, parameter :: nk            = 7
    integer, parameter :: nr            = nk*(nk-1)/2
    integer, parameter :: mg7_dnpb(0:9) = [1, 1, 3, 8, 14, 27, 54, 93, 162, 281]
    integer, parameter :: mg7_npb(-1:9) = [0, 1, 2, 5, 13, 27, 54, 108, 201, 363, 644]
    integer, parameter :: mg7_nsc(-1:9) = [0, 1, 1, 1, 1, 9, 34, 104, 273, 723, 1869]
    integer, parameter :: mg7_nb(0:9)   = [1, 2, 5, 13, 35, 95, 268, 739, 2042, 5552]

    contains

    !> @brief Build the 739 degree-seven invariant basis values from H-H distances in bohr.
    subroutine h7p_invariant_basis( distances, basis )
        real(wp), intent(in)  :: distances(nk,nk)
        real(wp), intent(out) :: basis(mg7_nb(7))
        real(wp)              :: transformed(0:nk-1,0:nk-1)
        real(wp)              :: mean_distance, cutoff, squared_distance_sum
        integer               :: i, j
        transformed          = 0.0_wp
        squared_distance_sum = 0.0_wp
        do j=1, nk
            do i=1, j-1
                transformed(i-1,j-1) = exp(-distances(i,j)/2.0_wp)
                transformed(j-1,i-1) = transformed(i-1,j-1)
                squared_distance_sum = squared_distance_sum+distances(i,j)**2
            enddo
        enddo
        mean_distance = sqrt(squared_distance_sum/real(nr, wp))
        if( mean_distance < 1.0e6_wp )then
            cutoff = (1.0_wp-mean_distance/1.0e6_wp)**3
        else
            cutoff = 0.0_wp
        endif
        call mg7_base(7, transformed, basis)
        basis = cutoff*basis
    end subroutine h7p_invariant_basis

    !> @brief Evaluate the released H7+ cg7_prims invariant-basis operation.
    subroutine cg7_prims( x, u )
        real(wp), intent (in)  :: x(0:)
        real(wp), intent (out) :: u(0:)
        integer                :: i0, i1, i2, i3, k
        real(wp)               :: d(0:nk-1,0:nk-1), t1(0:nk-1)
        if( size(x) .ne. nr.or.size(u) .ne. nr )then
            return
        endif
        ! recall revlex order,
        ! x(0:20) = (/ d01, d02, d12, d03, d13, d23, d04, d14, d24, d34, d05, ... /)
        k = 0
        do i1 = 0, nk-1
            do i0 = 0, i1-1
                d(i0,i1) = x(k)
                d(i1,i0) = x(k)
                k        = k+1
            enddo
            d(i1,i1) = 0
        enddo
        t1(0) = (x(0)+x(1)+x(3)+x(6)+x(10)+x(15))/6
        t1(1) = (x(0)+x(2)+x(4)+x(7)+x(11)+x(16))/6
        t1(2) = (x(1)+x(2)+x(5)+x(8)+x(12)+x(17))/6
        t1(3) = (x(3)+x(4)+x(5)+x(9)+x(13)+x(18))/6
        t1(4) = (x(6)+x(7)+x(8)+x(9)+x(14)+x(19))/6
        t1(5) = (x(10)+x(11)+x(12)+x(13)+x(14)+x(20))/6
        t1(6) = (x(15)+x(16)+x(17)+x(18)+x(19)+x(20))/6
        u     = 0
        u(0)  = sum(x)/size(x)
        u(1)  = sum(t1**2)/size(t1)
        u(2)  = sum(x**2)/size(x)
        u(3)  = sum(t1**3)/size(t1)
        u(4)  = sum(x**3)/size(x)
        u(8)  = sum(t1**4)/size(t1)
        u(9)  = sum(x**4)/size(x)
        u(11) = sum(t1**5)/size(t1)
        u(12) = sum(x**5)/size(x)
        u(14) = sum(t1**6)/size(t1)
        u(15) = sum(x**6)/size(x)
        u(16) = sum(t1**7)/size(t1)
        u(17) = sum(x**7)/size(x)
        u(19) = sum(x**10)/size(x)
        u(20) = sum(x**12)/size(x)
        do i0 = 0, nk-1
            do i1 = 0, nk-1
                if( i1 .ne. i0 )then
                    do i2 = 0, nk-1
                        if( i2 .ne. i0.and.i2 .ne. i1 )then
                            u(5)  = u(5)+d(i0,i1)**2*d(i0,i2)/210
                            u(6)  = u(6)+d(i0,i1)*d(i0,i2)*d(i1,i2)/210
                            u(10) = u(10)+d(i0,i1)**3*d(i0,i2)/210
                            u(13) = u(13)+d(i0,i1)**4*d(i0,i2)/210
                            u(18) = u(18)+d(i0,i1)**6*d(i0,i2)/210
                            do i3 = 0, nk-1
                                if( i3 .ne. i0.and.i3 .ne. i1.and.i3 .ne. i2 )then
                                    u(7) = u(7)+d(i0,i1)*d(i1,i2)*d(i0,i3)/840
                                endif
                            enddo
                        endif
                    enddo
                endif
            enddo
        enddo
    end subroutine cg7_prims

    !> @brief Evaluate the released H7+ mg7_prims invariant-basis operation.
    subroutine mg7_prims( r, u )
        real(wp), intent (in)  :: r(0:,0:)
        real(wp), intent (out) :: u(0:)
        real(wp)               :: x(0:nr-1)
        integer                :: i, j, pair
        if( size(r, 1) .ne. nk.or.size(r, 2) .ne. nk.or.size(u) .ne. nr )then
            return
        endif
        pair = 0
        do j=0, nk-1
            do i=0, j-1
                x(pair) = r(i,j)
                pair    = pair+1
            enddo
        enddo
        call cg7_prims(x, u)
    end subroutine mg7_prims

    !> @brief Evaluate the released H7+ mg7_prib invariant-basis operation.
    subroutine mg7_prib( mxd, u, w )
        integer, intent (in)   :: mxd
        real(wp), intent (in)  :: u(0:)
        real(wp), intent (out) :: w(0:mg7_npb(mxd)-1)
        if( size(u) .ne. nr )then
            return
        endif
        ! the following code was obtained using these parameters
        ! prims, dnpr(0:*): 0 1 2 5 3 3 2 3 0 0
        ! prib,  dnpb(0:*): 1 1 3 8 14 27 54 93 162 281
        ! constant term
        w(0) = 1
        ! terms of degree 1
        if( 1 .le. mxd )then
            w(1) = u(0)*w(0)
        endif
        ! terms of degree 2
        if( 2 .le. mxd )then
            w(2) = u(0)*w(1)
            w(3) = u(1)*w(0)
            w(4) = u(2)*w(0)
        endif
        ! terms of degree 3
        if( 3 .le. mxd )then
            w(5:7) = u(0)*w(2:4)
            w(8)   = u(3)*w(0)
            w(9)   = u(4)*w(0)
            w(10)  = u(5)*w(0)
            w(11)  = u(6)*w(0)
            w(12)  = u(7)*w(0)
        endif
        ! terms of degree 4
        if( 4 .le. mxd )then
            w(13:20) = u(0)*w(5:12)
            w(21:22) = u(1)*w(3:4)
            w(23)    = u(2)*w(4)
            w(24)    = u(8)*w(0)
            w(25)    = u(9)*w(0)
            w(26)    = u(10)*w(0)
        endif
        ! terms of degree 5
        if( 5 .le. mxd )then
            w(27:40) = u(0)*w(13:26)
            w(41:45) = u(1)*w(8:12)
            w(46:50) = u(2)*w(8:12)
            w(51)    = u(11)*w(0)
            w(52)    = u(12)*w(0)
            w(53)    = u(13)*w(0)
        endif
        ! terms of degree 6
        if( 6 .le. mxd )then
            w(54:80)   = u(0)*w(27:53)
            w(81:86)   = u(1)*w(21:26)
            w(87:90)   = u(2)*w(23:26)
            w(91:95)   = u(3)*w(8:12)
            w(96:99)   = u(4)*w(9:12)
            w(100:102) = u(5)*w(10:12)
            w(103:104) = u(6)*w(11:12)
            w(105)     = u(7)*w(12)
            w(106)     = u(14)*w(0)
            w(107)     = u(15)*w(0)
        endif
        ! terms of degree 7
        if( 7 .le. mxd )then
            w(108:161) = u(0)*w(54:107)
            w(162:174) = u(1)*w(41:53)
            w(175:182) = u(2)*w(46:53)
            w(183:185) = u(3)*w(24:26)
            w(186:188) = u(4)*w(24:26)
            w(189:191) = u(5)*w(24:26)
            w(192:194) = u(6)*w(24:26)
            w(195:197) = u(7)*w(24:26)
            w(198)     = u(16)*w(0)
            w(199)     = u(17)*w(0)
            w(200)     = u(18)*w(0)
        endif
        ! terms of degree 8
        if( 8 .le. mxd )then
            w(201:293) = u(0)*w(108:200)
            w(294:320) = u(1)*w(81:107)
            w(321:341) = u(2)*w(87:107)
            w(342:344) = u(3)*w(51:53)
            w(345:347) = u(4)*w(51:53)
            w(348:350) = u(5)*w(51:53)
            w(351:353) = u(6)*w(51:53)
            w(354:356) = u(7)*w(51:53)
            w(357:359) = u(8)*w(24:26)
            w(360:361) = u(9)*w(25:26)
            w(362)     = u(10)*w(26)
        endif
        ! terms of degree 9
        if( 9 .le. mxd )then
            w(363:524) = u(0)*w(201:362)
            w(525:563) = u(1)*w(162:200)
            w(564:589) = u(2)*w(175:200)
            w(590:606) = u(3)*w(91:107)
            w(607:618) = u(4)*w(96:107)
            w(619:626) = u(5)*w(100:107)
            w(627:631) = u(6)*w(103:107)
            w(632:634) = u(7)*w(105:107)
            w(635:637) = u(8)*w(51:53)
            w(638:640) = u(9)*w(51:53)
            w(641:643) = u(10)*w(51:53)
        endif
    end subroutine mg7_prib

    !> @brief Evaluate the released H7+ mgx_setd invariant-basis operation.
    subroutine mgx_setd( r, d, d2, d3, d4, d5, d6, d7, d8, d9 )
        real(wp), intent (in)            :: r(0:,0:)
        real(wp), intent (out)           :: d(0:,0:)
        real(wp), optional, intent (out) :: d2(0:,0:), d3(0:,0:), &
            d4(0:,0:), d5(0:,0:), d6(0:,0:), d7(0:,0:), d8(0:, 0:), d9(0:, 0:)
        integer :: i, j
        if (size(d, 1).ne.size(r, 1).or.size(d, 2).ne.size(r, 2).or. &
            size(r, 1).ne.size(r, 2)) then
        return
        endif
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

    !> @brief Evaluate the released H7+ mg7_isecs invariant-basis operation.
    subroutine mg7_isecs( mxd, r, pv )
        integer, intent (in)   :: mxd
        real(wp), intent (in)  :: r(0:,0:)
        real(wp), intent (out) :: pv(0:)
        integer, parameter     :: m = nk, m2=m*(m-1), m3=m2*(m-2), &
            m4 = m3*(m-3), m5=m4*(m-4), mm=m5*(m-5)
        integer            :: i0, i1, i2, i3, i4, i5
        real(wp)           :: d(0:nk-1,0:nk-1), d2(0:nk-1,0:nk-1), d3(0:nk-1,0:nk-1), d4(0:nk-1,0:nk-1), d5(0:nk-1,0:nk-1), d6(0:nk-1,0:nk-1), d7(0:nk-1,0:nk-1)
        integer, parameter :: aa(17) = (/0, 1, 8, 9, 10, 33, 34, 35, 36, 37, 38, 103, 104, 105, 106, 107, 108/)
        integer, parameter :: bb(75) = (/2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 39, 40, 41, &
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 109, 110, 111, 112, 113, 114, &
            115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, &
            135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145/)
        integer, parameter :: cc(125) = (/7, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 61, 62, 63, 64, 65, 66, 67, &
            68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 146, &
            147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, &
            167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, &
            187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, &
            207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, &
            227/)
        integer, parameter :: dd(55) = (/32, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 228, 229, 230, 231, &
            232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, &
            252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271/)
        if( size(r) .ne. nk*nk )then
            return
        endif
        call mgx_setd(r, d, d2, d3, d4, d5, d6, d7)
        pv = 0
        do i0 = 0, m-1
            do i1 = 0, m-1
                if( i1 .ne. i0 )then
                    do i2 = 0, m-1
                        if( i2 .ne. i0.and.i2 .ne. i1 )then
                            if( 4 .le. mxd )then
                                pv(0) = pv(0)+d2(i0,i1)*d2(i0,i2)
                                pv(1) = pv(1)+d2(i0,i1)*d(i0,i2)*d(i1,i2)
                            endif
                            if( 5 .le. mxd )then
                                pv(8)  = pv(8)+d3(i0,i1)*d2(i0,i2)
                                pv(9)  = pv(9)+d3(i0,i1)*d(i0,i2)*d(i1,i2)
                                pv(10) = pv(10)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)
                            endif
                            if( 6 .le. mxd )then
                                pv(33) = pv(33)+d5(i0,i1)*d(i0,i2)
                                pv(34) = pv(34)+d4(i0,i1)*d2(i0,i2)
                                pv(35) = pv(35)+d3(i0,i1)*d3(i0,i2)
                                pv(36) = pv(36)+d4(i0,i1)*d(i0,i2)*d(i1,i2)
                                pv(37) = pv(37)+d3(i0,i1)*d2(i0,i2)*d(i1,i2)
                                pv(38) = pv(38)+d2(i0,i1)*d2(i0,i2)*d2(i1,i2)
                            endif
                            if( 7 .le. mxd )then
                                pv(103) = pv(103)+d5(i0,i1)*d2(i0,i2)
                                pv(104) = pv(104)+d4(i0,i1)*d3(i0,i2)
                                pv(105) = pv(105)+d5(i0,i1)*d(i0,i2)*d(i1,i2)
                                pv(106) = pv(106)+d4(i0,i1)*d2(i0,i2)*d(i1,i2)
                                pv(107) = pv(107)+d3(i0,i1)*d3(i0,i2)*d(i1,i2)
                                pv(108) = pv(108)+d3(i0,i1)*d2(i0,i2)*d2(i1,i2)
                            endif
                            do i3 = 0, m-1
                                if( i3 .ne. i0.and.i3 .ne. i1.and.i3 .ne. i2 )then
                                    if( 4 .le. mxd )then
                                        pv(2) = pv(2)+d2(i0,i1)*d(i0,i2)*d(i0,i3)
                                        pv(3) = pv(3)+d2(i0,i1)*d(i1,i2)*d(i0,i3)
                                        pv(4) = pv(4)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)
                                        pv(5) = pv(5)+d(i0,i1)*d2(i1,i2)*d(i0,i3)
                                        pv(6) = pv(6)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                    endif
                                    if( 5 .le. mxd )then
                                        pv(11) = pv(11)+d3(i0,i1)*d(i0,i2)*d(i0,i3)
                                        pv(12) = pv(12)+d2(i0,i1)*d2(i0,i2)*d(i0,i3)
                                        pv(13) = pv(13)+d3(i0,i1)*d(i1,i2)*d(i0,i3)
                                        pv(14) = pv(14)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)
                                        pv(15) = pv(15)+d2(i0,i1)*d2(i1,i2)*d(i0,i3)
                                        pv(16) = pv(16)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,i3)
                                        pv(17) = pv(17)+d(i0,i1)*d3(i1,i2)*d(i0,i3)
                                        pv(18) = pv(18)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,i3)
                                        pv(19) = pv(19)+d(i0,i1)*d2(i1,i2)*d2(i0,i3)
                                        pv(20) = pv(20)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(21) = pv(21)+d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                    endif
                                    if( 6 .le. mxd )then
                                        pv(39) = pv(39)+d4(i0,i1)*d(i0,i2)*d(i0,i3)
                                        pv(40) = pv(40)+d3(i0,i1)*d2(i0,i2)*d(i0,i3)
                                        pv(41) = pv(41)+d4(i0,i1)*d(i1,i2)*d(i0,i3)
                                        pv(42) = pv(42)+d3(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)
                                        pv(43) = pv(43)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i0,i3)
                                        pv(44) = pv(44)+d3(i0,i1)*d2(i1,i2)*d(i0,i3)
                                        pv(45) = pv(45)+d2(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,i3)
                                        pv(46) = pv(46)+d2(i0,i1)*d3(i1,i2)*d(i0,i3)
                                        pv(47) = pv(47)+d(i0,i1)*d(i0,i2)*d3(i1,i2)*d(i0,i3)
                                        pv(48) = pv(48)+d(i0,i1)*d4(i1,i2)*d(i0,i3)
                                        pv(49) = pv(49)+d2(i0,i1)*d2(i0,i2)*d2(i0,i3)
                                        pv(50) = pv(50)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,i3)
                                        pv(51) = pv(51)+d2(i0,i1)*d2(i1,i2)*d2(i0,i3)
                                        pv(52) = pv(52)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d2(i0,i3)
                                        pv(53) = pv(53)+d(i0,i1)*d3(i1,i2)*d2(i0,i3)
                                        pv(54) = pv(54)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d3(i0,i3)
                                        pv(55) = pv(55)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(56) = pv(56)+d(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(57) = pv(57)+d3(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(58) = pv(58)+d2(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(59) = pv(59)+d(i0,i2)*d2(i1,i2)*d2(i0,i3)*d(i1,i3)
                                        pv(60) = pv(60)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i2,i3)
                                    endif
                                    if( 7 .le. mxd )then
                                        pv(109) = pv(109)+d5(i0,i1)*d(i0,i2)*d(i0,i3)
                                        pv(110) = pv(110)+d4(i0,i1)*d2(i0,i2)*d(i0,i3)
                                        pv(111) = pv(111)+d3(i0,i1)*d3(i0,i2)*d(i0,i3)
                                        pv(112) = pv(112)+d5(i0,i1)*d(i1,i2)*d(i0,i3)
                                        pv(113) = pv(113)+d4(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)
                                        pv(114) = pv(114)+d3(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i0,i3)
                                        pv(115) = pv(115)+d4(i0,i1)*d2(i1,i2)*d(i0,i3)
                                        pv(116) = pv(116)+d3(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,i3)
                                        pv(117) = pv(117)+d2(i0,i1)*d2(i0,i2)*d2(i1,i2)*d(i0,i3)
                                        pv(118) = pv(118)+d3(i0,i1)*d3(i1,i2)*d(i0,i3)
                                        pv(119) = pv(119)+d2(i0,i1)*d(i0,i2)*d3(i1,i2)*d(i0,i3)
                                        pv(120) = pv(120)+d2(i0,i1)*d4(i1,i2)*d(i0,i3)
                                        pv(121) = pv(121)+d(i0,i1)*d(i0,i2)*d4(i1,i2)*d(i0,i3)
                                        pv(122) = pv(122)+d(i0,i1)*d5(i1,i2)*d(i0,i3)
                                        pv(123) = pv(123)+d3(i0,i1)*d2(i0,i2)*d2(i0,i3)
                                        pv(124) = pv(124)+d3(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,i3)
                                        pv(125) = pv(125)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)*d2(i0,i3)
                                        pv(126) = pv(126)+d3(i0,i1)*d2(i1,i2)*d2(i0,i3)
                                        pv(127) = pv(127)+d2(i0,i1)*d(i0,i2)*d2(i1,i2)*d2(i0,i3)
                                        pv(128) = pv(128)+d2(i0,i1)*d3(i1,i2)*d2(i0,i3)
                                        pv(129) = pv(129)+d(i0,i1)*d(i0,i2)*d3(i1,i2)*d2(i0,i3)
                                        pv(130) = pv(130)+d(i0,i1)*d4(i1,i2)*d2(i0,i3)
                                        pv(131) = pv(131)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d3(i0,i3)
                                        pv(132) = pv(132)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d3(i0,i3)
                                        pv(133) = pv(133)+d(i0,i1)*d3(i1,i2)*d3(i0,i3)
                                        pv(134) = pv(134)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d4(i0,i3)
                                        pv(135) = pv(135)+d3(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(136) = pv(136)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(137) = pv(137)+d(i0,i1)*d3(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(138) = pv(138)+d4(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(139) = pv(139)+d(i0,i1)*d2(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(140) = pv(140)+d3(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i1,i3)
                                        pv(141) = pv(141)+d(i0,i1)*d2(i0,i2)*d(i1,i2)*d2(i0,i3)*d(i1,i3)
                                        pv(142) = pv(142)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d2(i0,i3)*d(i1,i3)
                                        pv(143) = pv(143)+d2(i0,i2)*d2(i1,i2)*d2(i0,i3)*d(i1,i3)
                                        pv(144) = pv(144)+d(i0,i2)*d3(i1,i2)*d2(i0,i3)*d(i1,i3)
                                        pv(145) = pv(145)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i2,i3)
                                    endif
                                    do i4 = 0, m-1
                                        if( i4 .ne. i0.and.i4 .ne. i1.and.i4 .ne. i2.and.i4 .ne. i3 )then
                                            if( 4 .le. mxd )then
                                                pv(7) = pv(7)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,i4)
                                            endif
                                            if( 5 .le. mxd )then
                                                pv(22) = pv(22)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,i4)
                                                pv(23) = pv(23)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(24) = pv(24)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(25) = pv(25)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(26) = pv(26)+d(i0,i1)*d(i1,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(27) = pv(27)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(28) = pv(28)+d2(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(29) = pv(29)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(30) = pv(30)+d(i0,i1)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(31) = pv(31)+d(i1,i2)*d(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                            endif
                                            if( 6 .le. mxd )then
                                                pv(61) = pv(61)+d3(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,i4)
                                                pv(62) = pv(62)+d2(i0,i1)*d2(i0,i2)*d(i0,i3)*d(i0,i4)
                                                pv(63) = pv(63)+d3(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(64) = pv(64)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(65) = pv(65)+d2(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(66) = pv(66)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(67) = pv(67)+d(i0,i1)*d3(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(68) = pv(68)+d2(i0,i1)*d(i1,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(69) = pv(69)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(70) = pv(70)+d(i0,i1)*d2(i1,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(71) = pv(71)+d(i0,i1)*d(i1,i2)*d3(i0,i3)*d(i0,i4)
                                                pv(72) = pv(72)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(73) = pv(73)+d(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(74) = pv(74)+d3(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(75) = pv(75)+d2(i0,i2)*d2(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(76) = pv(76)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(77) = pv(77)+d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(78) = pv(78)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(79) = pv(79)+d(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(80) = pv(80)+d2(i1,i2)*d2(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(81) = pv(81)+d(i0,i1)*d2(i1,i2)*d2(i1,i3)*d(i0,i4)
                                                pv(82) = pv(82)+d(i0,i2)*d2(i1,i2)*d2(i1,i3)*d(i0,i4)
                                                pv(83) = pv(83)+d2(i0,i1)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(84) = pv(84)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(85) = pv(85)+d(i0,i1)*d2(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(86) = pv(86)+d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(87) = pv(87)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d2(i0,i4)
                                                pv(88) = pv(88)+d(i0,i1)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d2(i0,i4)
                                                pv(89) = pv(89)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)*d(i1,i4)
                                                pv(90) = pv(90)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(91) = pv(91)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(92) = pv(92)+d2(i1,i2)*d(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                            endif
                                            if( 7 .le. mxd )then
                                                pv(146) = pv(146)+d4(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,i4)
                                                pv(147) = pv(147)+d3(i0,i1)*d2(i0,i2)*d(i0,i3)*d(i0,i4)
                                                pv(148) = pv(148)+d4(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(149) = pv(149)+d3(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(150) = pv(150)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(151) = pv(151)+d3(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(152) = pv(152)+d2(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(153) = pv(153)+d2(i0,i1)*d3(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(154) = pv(154)+d(i0,i1)*d(i0,i2)*d3(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(155) = pv(155)+d(i0,i1)*d4(i1,i2)*d(i0,i3)*d(i0,i4)
                                                pv(156) = pv(156)+d2(i0,i1)*d2(i0,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(157) = pv(157)+d3(i0,i1)*d(i1,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(158) = pv(158)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(159) = pv(159)+d2(i0,i1)*d2(i1,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(160) = pv(160)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(161) = pv(161)+d(i0,i1)*d3(i1,i2)*d2(i0,i3)*d(i0,i4)
                                                pv(162) = pv(162)+d2(i0,i1)*d(i1,i2)*d3(i0,i3)*d(i0,i4)
                                                pv(163) = pv(163)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d3(i0,i3)*d(i0,i4)
                                                pv(164) = pv(164)+d(i0,i1)*d2(i1,i2)*d3(i0,i3)*d(i0,i4)
                                                pv(165) = pv(165)+d(i0,i1)*d(i1,i2)*d4(i0,i3)*d(i0,i4)
                                                pv(166) = pv(166)+d3(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(167) = pv(167)+d2(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(168) = pv(168)+d(i0,i1)*d3(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(169) = pv(169)+d4(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(170) = pv(170)+d(i0,i1)*d2(i0,i2)*d2(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(171) = pv(171)+d3(i0,i2)*d2(i1,i2)*d(i1,i3)*d(i0,i4)
                                                pv(172) = pv(172)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(173) = pv(173)+d(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(174) = pv(174)+d3(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(175) = pv(175)+d2(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(176) = pv(176)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(177) = pv(177)+d2(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(178) = pv(178)+d(i0,i1)*d3(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(179) = pv(179)+d(i0,i2)*d3(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(180) = pv(180)+d2(i0,i2)*d(i1,i2)*d2(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(181) = pv(181)+d(i0,i1)*d2(i1,i2)*d2(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(182) = pv(182)+d(i0,i2)*d2(i1,i2)*d2(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(183) = pv(183)+d3(i1,i2)*d2(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(184) = pv(184)+d2(i1,i2)*d3(i0,i3)*d(i1,i3)*d(i0,i4)
                                                pv(185) = pv(185)+d2(i0,i1)*d2(i1,i2)*d2(i1,i3)*d(i0,i4)
                                                pv(186) = pv(186)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d2(i1,i3)*d(i0,i4)
                                                pv(187) = pv(187)+d2(i0,i2)*d2(i1,i2)*d2(i1,i3)*d(i0,i4)
                                                pv(188) = pv(188)+d(i0,i1)*d3(i1,i2)*d2(i1,i3)*d(i0,i4)
                                                pv(189) = pv(189)+d(i0,i2)*d3(i1,i2)*d2(i1,i3)*d(i0,i4)
                                                pv(190) = pv(190)+d(i0,i2)*d2(i1,i2)*d(i0,i3)*d2(i1,i3)*d(i0,i4)
                                                pv(191) = pv(191)+d3(i1,i2)*d(i0,i3)*d2(i1,i3)*d(i0,i4)
                                                pv(192) = pv(192)+d3(i0,i1)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(193) = pv(193)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(194) = pv(194)+d2(i0,i1)*d2(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(195) = pv(195)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(196) = pv(196)+d(i0,i1)*d3(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(197) = pv(197)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i2, &
                                                    & i3)* d(i0,i4)
                                                pv(198) = pv(198)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(199) = pv(199)+d3(i1,i2)*d(i0,i3)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(200) = pv(200)+d2(i1,i2)*d2(i0,i3)*d(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(201) = pv(201)+d(i0,i1)*d2(i1,i2)*d2(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(202) = pv(202)+d(i0,i2)*d2(i1,i2)*d2(i1,i3)*d(i2,i3)*d(i0,i4)
                                                pv(203) = pv(203)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,i3)*d2(i0,i4)
                                                pv(204) = pv(204)+d(i0,i1)*d2(i1,i2)*d2(i0,i3)*d2(i0,i4)
                                                pv(205) = pv(205)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d2(i0,i4)
                                                pv(206) = pv(206)+d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d2(i0,i4)
                                                pv(207) = pv(207)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d2(i0,i4)
                                                pv(208) = pv(208)+d(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d2(i0,i4)
                                                pv(209) = pv(209)+d2(i0,i1)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d2(i0,i4)
                                                pv(210) = pv(210)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d2(i0,i4)
                                                pv(211) = pv(211)+d(i0,i1)*d2(i1,i2)*d(i1,i3)*d(i2,i3)*d2(i0,i4)
                                                pv(212) = pv(212)+d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i2,i3)*d2(i0,i4)
                                                pv(213) = pv(213)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d3(i0,i4)
                                                pv(214) = pv(214)+d(i0,i1)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d3(i0,i4)
                                                pv(215) = pv(215)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0, &
                                                    & i4)* d(i1,i4)
                                                pv(216) = pv(216)+d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)*d(i1,i4)
                                                pv(217) = pv(217)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(218) = pv(218)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(219) = pv(219)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i2,i3)*d(i0, &
                                                    & i4)* d(i1,i4)
                                                pv(220) = pv(220)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(221) = pv(221)+d3(i1,i2)*d(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(222) = pv(222)+d(i0,i2)*d(i1,i2)*d2(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(223) = pv(223)+d2(i1,i2)*d2(i0,i3)*d(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(224) = pv(224)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i2,i3)*d(i0, &
                                                    & i4)* d(i1,i4)
                                                pv(225) = pv(225)+d(i0,i1)*d(i0,i2)*d(i0,i3)*d2(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(226) = pv(226)+d(i0,i1)*d(i1,i2)*d(i0,i3)*d2(i2,i3)*d(i0,i4)*d(i1,i4)
                                                pv(227) = pv(227)+d2(i1,i2)*d(i0,i3)*d2(i2,i3)*d(i0,i4)*d(i1,i4)
                                            endif
                                            do i5 = 0, m-1
                                                if (i5.ne.i0.and.i5.ne.i1.and.i5.ne.i2.and.i5.ne.i3.and. &
                                                    i5.ne.i4) then
                                                if( 5 .le. mxd )then
                                                    pv(32) = pv(32)+d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                endif
                                                if( 6 .le. mxd )then
                                                    pv(93) = pv(93)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(94) = pv(94)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(95)  = pv(95)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(96) = pv(96)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(97) = pv(97)+d2(i0,i1)*d(i1,i2)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(98)  = pv(98)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(99) = pv(99)+d2(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(100) = pv(100)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(101) = pv(101)+d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(102) = pv(102)+d(i0,i2)*d(i1,i2)*d(i1,i3)*d2(i0,i4)*d(i0,i5)
                                                endif
                                                if( 7 .le. mxd )then
                                                    pv(228) = pv(228)+d3(i0,i1)*d(i0,i2)*d(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(229) = pv(229)+d2(i0,i1)*d2(i0,i2)*d(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(230) = pv(230)+d3(i0,i1)*d(i1,i2)*d(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(231) = pv(231)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(232) = pv(232)+d2(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(233) = pv(233)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(234) = pv(234)+d(i0,i1)*d3(i1,i2)*d(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(235) = pv(235)+d2(i0,i1)*d(i1,i2)*d2(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(236) = pv(236)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d2(i0,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(237) = pv(237)+d(i0,i1)*d2(i1,i2)*d2(i0,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(238) = pv(238)+d3(i0,i1)*d(i1,i2)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(239) = pv(239)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(240) = pv(240)+d(i0,i1)*d2(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(241) = pv(241)+d3(i0,i2)*d(i1,i2)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(242) = pv(242)+d2(i0,i1)*d2(i1,i2)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(243) = pv(243)+d(i0,i1)*d(i0,i2)*d2(i1,i2)*d(i1,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(244) = pv(244)+d2(i0,i2)*d2(i1,i2)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(245) = pv(245)+d(i0,i2)*d3(i1,i2)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(246) = pv(246)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)* &
                                                        & d(i0,i4)* d(i0,i5)
                                                    pv(247) = pv(247)+d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(248) = pv(248)+d(i0,i1)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(249) = pv(249)+d(i0,i2)*d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(250) = pv(250)+d3(i1,i2)*d(i0,i3)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(251) = pv(251)+d2(i1,i2)*d2(i0,i3)*d(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(252) = pv(252)+d(i0,i1)*d2(i1,i2)*d2(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(253) = pv(253)+d(i0,i2)*d2(i1,i2)*d2(i1,i3)*d(i0,i4)*d(i0,i5)
                                                    pv(254) = pv(254)+d2(i0,i1)*d(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(255) = pv(255)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d(i2,i3)* &
                                                        & d(i0,i4)* d(i0,i5)
                                                    pv(256) = pv(256)+d(i0,i1)*d2(i1,i2)*d(i1,i3)*d(i2,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(257) = pv(257)+d2(i1,i2)*d(i0,i3)*d(i1,i3)*d(i2,i3)*d(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(258) = pv(258)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i1,i3)*d2(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(259) = pv(259)+d2(i0,i2)*d(i1,i2)*d(i1,i3)*d2(i0,i4)*d(i0,i5)
                                                    pv(260) = pv(260)+d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d2(i0,i4)* &
                                                        & d(i0,i5)
                                                    pv(261) = pv(261)+d2(i1,i2)*d(i0,i3)*d(i1,i3)*d2(i0,i4)*d(i0,i5)
                                                    pv(262) = pv(262)+d(i0,i2)*d(i1,i2)*d(i1,i3)*d3(i0,i4)*d(i0,i5)
                                                    pv(263) = pv(263)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)* &
                                                        & d(i1,i4)* d(i0,i5)
                                                    pv(264) = pv(264)+d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i1,i3)*d(i1,i4)* &
                                                        & d(i0,i5)
                                                    pv(265) = pv(265)+d2(i0,i1)*d2(i0,i2)*d(i2,i3)*d(i1,i4)*d(i0,i5)
                                                    pv(266) = pv(266)+d2(i0,i1)*d(i0,i2)*d(i1,i2)*d(i2,i3)*d(i1,i4)* &
                                                        & d(i0,i5)
                                                    pv(267) = pv(267)+d2(i0,i1)*d(i0,i2)*d(i0,i3)*d(i2,i3)*d(i1,i4)* &
                                                        & d(i0,i5)
                                                    pv(268) = pv(268)+d(i0,i1)*d2(i0,i2)*d(i0,i3)*d(i2,i3)*d(i1,i4)* &
                                                        & d(i0,i5)
                                                    pv(269) = pv(269)+d2(i0,i1)*d(i1,i2)*d(i0,i3)*d(i2,i3)*d(i1,i4)* &
                                                        & d(i0,i5)
                                                    pv(270) = pv(270)+d(i0,i1)*d(i0,i2)*d(i1,i2)*d(i0,i3)*d(i2,i3)* &
                                                        & d(i1,i4)* d(i0,i5)
                                                    pv(271) = pv(271)+d2(i0,i2)*d(i1,i2)*d(i0,i3)*d(i2,i3)*d(i1,i4)* &
                                                        & d(i0,i5)
                                                endif
                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    enddo
                endif
            enddo
            endif
        enddo
        enddo
        pv(aa) = pv(aa)/m3
        pv(bb) = pv(bb)/m4
        pv(cc) = pv(cc)/m5
        pv(dd) = pv(dd)/mm
    end subroutine mg7_isecs

    !> @brief Evaluate the released H7+ mg7_secs invariant-basis operation.
    subroutine mg7_secs( mxd, r, v )
        integer, intent (in)   :: mxd
        real(wp), intent (in)  :: r(0:,0:)
        real(wp), intent (out) :: v(0:)
        real(wp)               :: pv(0:size(v)-1) !! an upper bound size
        if( size(r) .ne. nk*nk.or.size(v) .ne. mg7_nsc(mxd) )then
            return
        endif
        call mg7_isecs(mxd, r, pv)
        v(0) = 1
        if( 4 .le. mxd )then
            v(1:8) = pv(0:7)
        endif
        if( 5 .le. mxd )then
            v(9:33) = pv(8:32)
        endif
        if( 6 .le. mxd )then
            v(34:103) = pv(33:102)
        endif
        if( 7 .le. mxd )then
            v(104:272) = pv(103:271)
        endif
        if( 8 .le. mxd )then
            return
        endif
    end subroutine mg7_secs

    !> @brief Evaluate the released H7+ mg7_base invariant-basis operation.
    subroutine mg7_base( mxd, r, w )
        integer, intent (in)   :: mxd
        real(wp), intent (in)  :: r(0:,0:)
        real(wp), intent (out) :: w(0:mg7_nb(mxd)-1)
        integer                :: k, l0, m0, l, n, ind
        real(wp)               :: u(0:nr-1), v(0:mg7_nsc(mxd)-1), &
            w0(0:mg7_npb(mxd)-1)
        if( size(r, 1) .ne. nk.or.size(r, 2) .ne. nk )then
            return
        endif
        call mg7_prims(r, u)
        call mg7_prib(mxd, u, w0)
        call mg7_secs(mxd, r, v)
        ind = 0
        do n = 0, mxd
            do k = 0, n
                l0 = mg7_nsc(k-1)
                m0 = mg7_npb(n-k-1)
                do l = l0, mg7_nsc(k)-1
                    w(ind:ind+mg7_dnpb(n-k)-1) = v(l)*w0(m0:mg7_npb(n-k)-1)
                    ind                        = ind+mg7_dnpb(n-k)
                enddo
            enddo
        enddo
        if( ind .ne. mg7_nb(mxd) )then
            return
        endif
    end subroutine mg7_base
end module chempotpy_h7p_pip_basis
