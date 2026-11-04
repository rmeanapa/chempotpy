!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
module chempotpy_h3_basis
    implicit none
    private
    integer, parameter :: maximum_terms = 5000
    integer, save      :: first_powers(maximum_terms)
    integer, save      :: second_powers(maximum_terms)
    integer, save      :: third_powers(maximum_terms)
    integer, save      :: coefficient_indices(maximum_terms)
    integer, save      :: coefficient_count = 0
    integer, save      :: term_count = 0
    public :: get_h3_basis_term
    public :: initialize_h3_basis
    contains
    !> @brief Build the two symmetric polynomial index partitions used by the 2002 H3 fits.
    subroutine initialize_h3_basis( primary_order, secondary_order, primary_term_count, secondary_term_count )
        integer, intent(in)  :: primary_order
        integer, intent(in)  :: secondary_order
        integer, intent(out) :: primary_term_count
        integer, intent(out) :: secondary_term_count
        coefficient_count = 0
        term_count        = 0
        call append_symmetric_basis(primary_order, primary_order)
        primary_term_count = term_count
        coefficient_count  = 0
        call append_symmetric_basis(secondary_order, secondary_order)
        secondary_term_count = term_count-primary_term_count
    end subroutine initialize_h3_basis
    !> @brief Return the powers and coefficient index for one symmetric H3 basis term.
    subroutine get_h3_basis_term( term_index, first_power, second_power, third_power, coefficient_index )
        integer, intent(in)  :: term_index
        integer, intent(out) :: first_power
        integer, intent(out) :: second_power
        integer, intent(out) :: third_power
        integer, intent(out) :: coefficient_index
        first_power       = first_powers(term_index)
        second_power      = second_powers(term_index)
        third_power       = third_powers(term_index)
        coefficient_index = coefficient_indices(term_index)
    end subroutine get_h3_basis_term
    !> @brief Append the permutation-expanded symmetric basis through the requested total order.
    subroutine append_symmetric_basis( maximum_order, top_order )
        integer, intent(in) :: maximum_order
        integer, intent(in) :: top_order
        integer             :: first_power
        integer             :: second_power
        integer             :: third_power
        integer             :: total_order
        do first_power=0, min(top_order, maximum_order)
            do second_power=first_power, min(top_order, maximum_order)
                do third_power=second_power, min(top_order, maximum_order)
                    total_order = first_power+second_power+third_power
                    if( total_order > maximum_order ) cycle
                    if ( &
                        total_order-first_power <= 0 .or. &
                        total_order-second_power <= 0 .or. &
                        total_order-third_power <= 0 &
                        ) cycle
                    coefficient_count = coefficient_count+1
                    call append_term(first_power, second_power, third_power)
                    if( first_power /= second_power )then
                        if( second_power /= third_power )then
                            call append_term(first_power, third_power, second_power)
                            call append_term(second_power, first_power, third_power)
                            call append_term(second_power, third_power, first_power)
                            call append_term(third_power, first_power, second_power)
                            call append_term(third_power, second_power, first_power)
                        else
                            call append_term(second_power, first_power, third_power)
                            call append_term(second_power, third_power, first_power)
                        endif
                    else if( second_power /= third_power )then
                        call append_term(first_power, third_power, second_power)
                        call append_term(third_power, first_power, second_power)
                    endif
                enddo
            enddo
        enddo
    end subroutine append_symmetric_basis
    !> @brief Append one ordered power triple to the private basis table.
    subroutine append_term( first_power, second_power, third_power )
        integer, intent(in) :: first_power
        integer, intent(in) :: second_power
        integer, intent(in) :: third_power
        term_count                      = term_count+1
        first_powers(term_count)        = first_power
        second_powers(term_count)       = second_power
        third_powers(term_count)        = third_power
        coefficient_indices(term_count) = coefficient_count
    end subroutine append_term
end module chempotpy_h3_basis
