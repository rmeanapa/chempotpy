!> @brief Provide reusable mathematical special functions used by PES angular expansions.
module chempotpy_special_functions
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: associated_legendre_polynomial

    contains

    !> @brief Evaluate an associated Legendre polynomial including the Condon--Shortley phase.
    pure function associated_legendre_polynomial( degree, order, argument ) result(value)
        integer, intent(in)  :: degree, order
        real(wp), intent(in) :: argument
        real(wp)             :: value
        real(wp)             :: factor, previous, previous_previous, root
        integer              :: current_degree, index
        value = 0.0_wp
        if( order < 0 .or. order > degree .or. abs(argument) > 1.0_wp ) return
        previous_previous = 1.0_wp
        if( order > 0 )then
            root   = sqrt(max(0.0_wp, (1.0_wp-argument)*(1.0_wp+argument)))
            factor = 1.0_wp
            do index=1, order
                previous_previous = -previous_previous*factor*root
                factor            = factor+2.0_wp
            enddo
        endif
        if( degree == order )then
            value = previous_previous
            return
        endif
        previous = argument*real(2*order+1, wp)*previous_previous
        if( degree == order+1 )then
            value = previous
            return
        endif
        do current_degree=order+2, degree
            value = (argument*real(2*current_degree-1, wp)*previous - real(current_degree+order-1, wp)* &
                & previous_previous)/real(current_degree-order, wp)
            previous_previous = previous
            previous          = value
        enddo
    end function associated_legendre_polynomial
end module chempotpy_special_functions
