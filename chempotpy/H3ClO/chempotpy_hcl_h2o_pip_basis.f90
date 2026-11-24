!> @brief Evaluate the specialized degree-six X3Y1Z1 invariant basis for HCl-H2O.
module chempotpy_hcl_h2o_pip_basis
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_hcl_h2o_basis
    integer, parameter :: degree                         = 6
    integer, parameter :: primary_degree_count(0:degree) = [1, 4, 13, 35, 83, 179, 361]
    integer, parameter :: primary_count(0:degree)        = [1, 5, 18, 53, 136, 315, 676]
    integer, parameter :: secondary_count(0:degree)      = [1, 1, 4, 11, 17, 23, 33]

    contains

    !> @brief Evaluate all 1,603 normalized X3Y1Z1 invariant polynomials through degree six.
    pure subroutine evaluate_hcl_h2o_basis( transformed_distances, basis )
        real(wp), intent(in)  :: transformed_distances(0:4,0:4)
        real(wp), intent(out) :: basis(0:1602)
        real(wp)              :: primary(0:675), secondary(0:32)
        integer               :: basis_index, degree_index, primary_start, secondary_index, secondary_start
        call evaluate_primary_basis(transformed_distances, primary)
        call evaluate_secondary_basis(transformed_distances, secondary)
        basis_index = 0
        do degree_index=0, degree
            do secondary_index=0, degree_index
                secondary_start = block_start(secondary_count, secondary_index)
                primary_start   = block_start(primary_count, degree_index-secondary_index)
                do while( secondary_start < secondary_count(secondary_index) )
                    basis(basis_index:basis_index+primary_degree_count(degree_index-secondary_index)- &
                        & 1) = secondary(secondary_start)*primary(primary_start:primary_count(degree_index- &
                        & secondary_index)-1)
                    basis_index     = basis_index+primary_degree_count(degree_index-secondary_index)
                    secondary_start = secondary_start+1
                enddo
            enddo
        enddo
    end subroutine evaluate_hcl_h2o_basis

    !> @brief Return the first index of one degree block in a cumulative count array.
    pure integer function block_start( cumulative_count, block_index ) result(first_index)
        integer, intent(in) :: cumulative_count(0:), block_index
        if( block_index == 0 )then
            first_index = 0
        else
            first_index = cumulative_count(block_index-1)
        endif
    end function block_start

    !> @brief Evaluate the primary-invariant monomial basis through degree six.
    pure subroutine evaluate_primary_basis( r, w )
        real(wp), intent(in)  :: r(0:4,0:4)
        real(wp), intent(out) :: w(0:675)
        real(wp)              :: hydrogen_pairs(3), hydrogen_oxygen(3), hydrogen_chlorine(3), averages(3), u(0:9)
        hydrogen_pairs    = [r(0,1), r(0,2), r(1,2)]
        hydrogen_oxygen   = [r(0,3), r(1,3), r(2,3)]
        hydrogen_chlorine = [r(0,4), r(1,4), r(2,4)]
        averages = [(hydrogen_pairs(1)+hydrogen_pairs(2))/2.0_wp, (hydrogen_pairs(1)+hydrogen_pairs(3))/2.0_wp, &
            & (hydrogen_pairs(2)+hydrogen_pairs(3))/2.0_wp]
        u(0)       = sum(hydrogen_pairs)/3.0_wp
        u(1)       = sum(hydrogen_oxygen)/3.0_wp
        u(2)       = sum(hydrogen_chlorine)/3.0_wp
        u(3)       = r(3,4)
        u(4)       = sum(averages**2)/3.0_wp
        u(5)       = sum(hydrogen_oxygen**2)/3.0_wp
        u(6)       = sum(hydrogen_chlorine**2)/3.0_wp
        u(7)       = sum(averages**3)/3.0_wp
        u(8)       = sum(hydrogen_oxygen**3)/3.0_wp
        u(9)       = sum(hydrogen_chlorine**3)/3.0_wp
        w(0)       = 1.0_wp
        w(1:4)     = u(0:3)*w(0)
        w(5:8)     = u(0)*w(1:4)
        w(9:11)    = u(1)*w(2:4)
        w(12:13)   = u(2)*w(3:4)
        w(14)      = u(3)*w(4)
        w(15:17)   = u(4:6)*w(0)
        w(18:30)   = u(0)*w(5:17)
        w(31:39)   = u(1)*w(9:17)
        w(40:45)   = u(2)*w(12:17)
        w(46:49)   = u(3)*w(14:17)
        w(50:52)   = u(7:9)*w(0)
        w(53:87)   = u(0)*w(18:52)
        w(88:109)  = u(1)*w(31:52)
        w(110:122) = u(2)*w(40:52)
        w(123:129) = u(3)*w(46:52)
        w(130:132) = u(4)*w(15:17)
        w(133:134) = u(5)*w(16:17)
        w(135)     = u(6)*w(17)
        w(136:218) = u(0)*w(53:135)
        w(219:266) = u(1)*w(88:135)
        w(267:292) = u(2)*w(110:135)
        w(293:305) = u(3)*w(123:135)
        w(306:308) = u(4)*w(50:52)
        w(309:311) = u(5)*w(50:52)
        w(312:314) = u(6)*w(50:52)
        w(315:493) = u(0)*w(136:314)
        w(494:589) = u(1)*w(219:314)
        w(590:637) = u(2)*w(267:314)
        w(638:659) = u(3)*w(293:314)
        w(660:665) = u(4)*w(130:135)
        w(666:668) = u(5)*w(133:135)
        w(669)     = u(6)*w(135)
        w(670:672) = u(7)*w(50:52)
        w(673:674) = u(8)*w(51:52)
        w(675)     = u(9)*w(52)
    end subroutine evaluate_primary_basis

    !> @brief Evaluate the secondary invariants through degree six.
    pure subroutine evaluate_secondary_basis( r, v )
        real(wp), intent(in)  :: r(0:4,0:4)
        real(wp), intent(out) :: v(0:32)
        real(wp)              :: pv(0:9), r2(0:4,0:4)
        integer               :: first_hydrogen, second_hydrogen
        r2 = r*r
        pv = 0.0_wp
        do first_hydrogen=0, 2
            pv(2) = pv(2)+r(first_hydrogen,3)*r(first_hydrogen,4)/3.0_wp
            pv(7) = pv(7)+r2(first_hydrogen,3)*r(first_hydrogen,4)/3.0_wp
            pv(9) = pv(9)+r(first_hydrogen,3)*r2(first_hydrogen,4)/3.0_wp
            do second_hydrogen=0, 2
                if( second_hydrogen == first_hydrogen ) cycle
                pv(0) = pv(0)+r(first_hydrogen,second_hydrogen)*r(first_hydrogen,3)/6.0_wp
                pv(1) = pv(1)+r(first_hydrogen,second_hydrogen)*r(first_hydrogen,4)/6.0_wp
                pv(3) = pv(3)+r2(first_hydrogen,second_hydrogen)*r(first_hydrogen,3)/6.0_wp
                pv(4) = pv(4)+r(first_hydrogen,second_hydrogen)*r2(first_hydrogen,3)/6.0_wp
                pv(5) = pv(5)+r2(first_hydrogen,second_hydrogen)*r(first_hydrogen,4)/6.0_wp
                pv(6) = pv(6)+r(first_hydrogen,second_hydrogen)*r(first_hydrogen,3)*r(first_hydrogen,4)/6.0_wp
                pv(8) = pv(8)+r(first_hydrogen,second_hydrogen)*r2(first_hydrogen,4)/6.0_wp
            enddo
        enddo
        v(0)    = 1.0_wp
        v(1:3)  = pv(0:2)
        v(4:10) = pv(3:9)
        v(11)   = pv(0)*pv(0); v(12) = pv(0)*pv(1); v(13) = pv(1)*pv(1)
        v(14)   = pv(0)*pv(2); v(15) = pv(1)*pv(2); v(16) = pv(2)*pv(2)
        v(17)   = pv(1)*pv(4); v(18) = pv(0)*pv(5); v(19) = pv(0)*pv(7)
        v(20)   = pv(0)*pv(8); v(21) = pv(0)*pv(9); v(22) = pv(1)*pv(9)
        v(23)   = pv(3)*pv(4); v(24) = pv(3)*pv(6); v(25) = pv(3)*pv(7)
        v(26)   = pv(3)*pv(8); v(27) = pv(5)*pv(8); v(28) = pv(3)*pv(9)
        v(29)   = pv(4)*pv(9); v(30) = pv(5)*pv(9); v(31) = pv(6)*pv(9); v(32) = pv(7)*pv(9)
    end subroutine evaluate_secondary_basis
end module chempotpy_hcl_h2o_pip_basis
