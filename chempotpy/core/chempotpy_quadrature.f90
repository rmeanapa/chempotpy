!> @brief Provide reusable numerical quadrature rules for PES kernels.
module chempotpy_quadrature
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: gauss_legendre_rule

    contains

    !> @brief Generate Gauss--Legendre nodes and weights over a finite interval.
    pure subroutine gauss_legendre_rule( lower_bound, upper_bound, nodes, weights )
        real(wp), intent(in)  :: lower_bound, upper_bound
        real(wp), intent(out) :: nodes(:), weights(size(nodes))
        integer               :: index, iteration, order, polynomial_index
        real(wp)              :: derivative, half_range, midpoint, polynomial, previous_polynomial
        real(wp)              :: previous_root, root, temporary
        order      = size(nodes)
        midpoint   = 0.5_wp*(upper_bound+lower_bound)
        half_range = 0.5_wp*(upper_bound-lower_bound)
        do index=1, (order+1)/2
            root = cos(PI*(real(index, wp)-0.25_wp)/(real(order, wp)+0.5_wp))
            do iteration=1, 100
                polynomial          = 1.0_wp
                previous_polynomial = 0.0_wp
                do polynomial_index=1, order
                    temporary           = previous_polynomial
                    previous_polynomial = polynomial
                    polynomial = ((2.0_wp*real(polynomial_index, wp)-1.0_wp)*root*previous_polynomial - &
                        & (real(polynomial_index, wp)-1.0_wp)*temporary)/real(polynomial_index, wp)
                enddo
                derivative    = real(order, wp)*(root*polynomial-previous_polynomial)/(root**2-1.0_wp)
                previous_root = root
                root          = previous_root-polynomial/derivative
                if( abs(root-previous_root) <= 128.0_wp*epsilon(1.0_wp) ) exit
            enddo
            nodes(index)           = midpoint-half_range*root
            nodes(order+1-index)   = midpoint+half_range*root
            weights(index)         = 2.0_wp*half_range/((1.0_wp-root**2)*derivative**2)
            weights(order+1-index) = weights(index)
        enddo
    end subroutine gauss_legendre_rule
end module chempotpy_quadrature
