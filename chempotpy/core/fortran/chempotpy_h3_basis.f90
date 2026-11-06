!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
module chempotpy_h3_basis
    implicit none
    private
    integer, parameter, public :: primary_term_count   = 418
    integer, parameter, public :: secondary_term_count = 330
    public :: get_h3_basis_term
    contains
    !> @brief Return one term from the fixed order-12/order-11 symmetric H3 basis.
    pure subroutine get_h3_basis_term( term_index, first_power, second_power, third_power, coefficient_index )
        integer, intent(in)  :: term_index
        integer, intent(out) :: first_power
        integer, intent(out) :: second_power
        integer, intent(out) :: third_power
        integer, intent(out) :: coefficient_index
        integer              :: current_coefficient, current_term, maximum_order
        integer              :: first_candidate, second_candidate, third_candidate, total_order
        integer              :: requested_term
        logical              :: found
        first_power = 0; second_power = 0; third_power = 0; coefficient_index = 0
        if( term_index < 1 .or. term_index > primary_term_count+secondary_term_count ) return
        if( term_index <= primary_term_count )then
            maximum_order  = 12
            requested_term = term_index
        else
            maximum_order  = 11
            requested_term = term_index-primary_term_count
        endif
        current_coefficient = 0
        current_term        = 0
        do first_candidate=0, maximum_order
            do second_candidate=first_candidate, maximum_order
                do third_candidate=second_candidate, maximum_order
                    total_order = first_candidate+second_candidate+third_candidate
                    if( total_order > maximum_order ) cycle
                    if( total_order-first_candidate <= 0 .or. &
                        & total_order-second_candidate <= 0 .or. &
                        & total_order-third_candidate <= 0 ) cycle
                    current_coefficient = current_coefficient+1
                    call select_term(first_candidate, second_candidate, third_candidate, current_coefficient, &
                        & requested_term, current_term, first_power, second_power, third_power, coefficient_index, found)
                    if( found ) return
                    if( first_candidate /= second_candidate )then
                        if( second_candidate /= third_candidate )then
                            call select_term(first_candidate, third_candidate, second_candidate, current_coefficient, &
                                & requested_term, current_term, first_power, second_power, third_power, &
                                & coefficient_index, found)
                            if( found ) return
                            call select_term(second_candidate, first_candidate, third_candidate, current_coefficient, &
                                & requested_term, current_term, first_power, second_power, third_power, &
                                & coefficient_index, found)
                            if( found ) return
                            call select_term(second_candidate, third_candidate, first_candidate, current_coefficient, &
                                & requested_term, current_term, first_power, second_power, third_power, &
                                & coefficient_index, found)
                            if( found ) return
                            call select_term(third_candidate, first_candidate, second_candidate, current_coefficient, &
                                & requested_term, current_term, first_power, second_power, third_power, &
                                & coefficient_index, found)
                            if( found ) return
                            call select_term(third_candidate, second_candidate, first_candidate, current_coefficient, &
                                & requested_term, current_term, first_power, second_power, third_power, &
                                & coefficient_index, found)
                            if( found ) return
                        else
                            call select_term(second_candidate, first_candidate, third_candidate, current_coefficient, &
                                & requested_term, current_term, first_power, second_power, third_power, &
                                & coefficient_index, found)
                            if( found ) return
                            call select_term(second_candidate, third_candidate, first_candidate, current_coefficient, &
                                & requested_term, current_term, first_power, second_power, third_power, &
                                & coefficient_index, found)
                            if( found ) return
                        endif
                    else if( second_candidate /= third_candidate )then
                        call select_term(first_candidate, third_candidate, second_candidate, current_coefficient, &
                            & requested_term, current_term, first_power, second_power, third_power, coefficient_index, found)
                        if( found ) return
                        call select_term(third_candidate, first_candidate, second_candidate, current_coefficient, &
                            & requested_term, current_term, first_power, second_power, third_power, coefficient_index, found)
                        if( found ) return
                    endif
                enddo
            enddo
        enddo
    end subroutine get_h3_basis_term
    !> @brief Match one generated basis term to the requested stateless lookup index.
    pure subroutine select_term( first_candidate, second_candidate, third_candidate, current_coefficient, &
            & requested_term, current_term, first_power, second_power, third_power, coefficient_index, found )
        integer, intent(in)    :: first_candidate
        integer, intent(in)    :: second_candidate
        integer, intent(in)    :: third_candidate
        integer, intent(in)    :: current_coefficient
        integer, intent(in)    :: requested_term
        integer, intent(inout) :: current_term
        integer, intent(out)   :: first_power
        integer, intent(out)   :: second_power
        integer, intent(out)   :: third_power
        integer, intent(out)   :: coefficient_index
        logical, intent(out)   :: found
        current_term = current_term+1
        found        = current_term == requested_term
        if( found )then
            first_power       = first_candidate
            second_power      = second_candidate
            third_power       = third_candidate
            coefficient_index = current_coefficient
        else
            first_power = 0; second_power = 0; third_power = 0; coefficient_index = 0
        endif
    end subroutine select_term
end module chempotpy_h3_basis
