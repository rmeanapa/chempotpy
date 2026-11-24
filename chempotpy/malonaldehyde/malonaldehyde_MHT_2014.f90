! System:                        malonaldehyde
! Functional form:               Morse terms and one-, two-, and three-dimensional Gaussian functions
! Common name:                   MHT malonaldehyde surface
! Reference:                     W. Mizukami, S. Habershon, and D. P. Tew, J. Chem. Phys. 141, 144310 (2014)
! Number of derivatives:         1
! Number of bodies:              9
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable fit coefficients are embedded
! Notes:                         Atom order is C-C-O-C-O-H-H-H-H; energy is relative to equilibrium.

!> @brief Evaluate the MHT 2014 enol-malonaldehyde potential energy surface.
module chempotpy_malonaldehyde_mht_2014_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    integer, parameter :: natoms                     = 9, ndist = natoms*(natoms-1)/2
    integer, parameter :: released_pair_order(ndist) = [ &
        1, 2, 9, 3, 10, 16, 4, 11, 17, 22, 5, 12, 18, 23, 27, 6, 13, 19, &
        & 24, 28, 31, 7, 14, 20, 25, 29, 32, 34, 8, 15, 21, 26, 30, 33, 35, 36 ]

    contains

    !> @brief Evaluate relative energy and, when requested, the analytic Cartesian gradient.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_geometry, only: all_pair_distances_and_jacobian, transform_internal_gradient_to_cartesian
        use chempotpy_status,   only: chempotpy_status_t, status_derivative_unavailable
        use chempotpy_units,    only: bohr_to_angstrom, hartree_per_bohr_to_ev_per_angstrom, &
            & hartree_to_electronvolt
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: all_distances(ndist), all_jacobian(ndist,3*natoms)
        real(wp)              :: cartesian_gradient(3*natoms), coordinates_bohr(natoms,3)
        real(wp)              :: distance_gradient(ndist), distance_jacobian(ndist,3*natoms)
        real(wp)              :: distances(ndist), energy_hartree
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad == 0 .or. igrad == 1 )then
            coordinates_bohr = x/bohr_to_angstrom
            call all_pair_distances_and_jacobian(coordinates_bohr, all_distances, all_jacobian)
            distances         = all_distances(released_pair_order)
            distance_jacobian = all_jacobian(released_pair_order,:)
            call evaluate_fit(distances, igrad == 1, energy_hartree, distance_gradient)
            p(1) = energy_hartree*hartree_to_electronvolt
            if( igrad == 1 )then
                call transform_internal_gradient_to_cartesian(distance_gradient, distance_jacobian, cartesian_gradient)
                cartesian_gradient = cartesian_gradient*hartree_per_bohr_to_ev_per_angstrom
                g(1,:,:)           = transpose(reshape(cartesian_gradient, [3, natoms]))
            endif
        elseif( present(status) )then
            call status%set(status_derivative_unavailable, &
                'malonaldehyde_MHT_2014 provides energy and first derivatives')
        endif
    end subroutine evaluate

    !> @brief Contract the released Morse and Gaussian basis in its original summation order.
    pure subroutine evaluate_fit( distances, calculate_gradient, energy, gradient )
        use chempotpy_malonaldehyde_mht_parameters, only: g1d, g2d, g3d, ig1d, ig2d, ig3d, imorse, morse, shift
        real(wp), intent(in)  :: distances(ndist)
        logical, intent(in)   :: calculate_gradient
        real(wp), intent(out) :: energy, gradient(ndist)
        real(wp)              :: term_distances(3), term_energy, term_gradient(3)
        integer               :: coordinate, distance_index, term
        energy   = shift
        gradient = 0.0_wp
        do term=1, size(imorse)
            distance_index = imorse(term)
            call evaluate_morse(morse(:, term), distances(distance_index), term_energy, term_gradient(1))
            energy                                            = energy+term_energy
            if( calculate_gradient ) gradient(distance_index) = gradient(distance_index)+term_gradient(1)
        enddo
        do term=1, size(ig1d)
            distance_index    = ig1d(term)
            term_distances(1) = distances(distance_index)
            call evaluate_gaussian(term_distances(1:1), g1d(1:1, term), g1d(2:2, term), g1d(3, term), &
                & g1d(4, term), term_energy, term_gradient(1:1))
            energy                                            = energy+term_energy
            if( calculate_gradient ) gradient(distance_index) = gradient(distance_index)+term_gradient(1)
        enddo
        do term=1, size(ig2d, 2)
            do coordinate=1, 2
                term_distances(coordinate) = distances(ig2d(coordinate, term))
            enddo
            call evaluate_gaussian(term_distances(1:2), g2d(1:2, term), g2d(3:4, term), g2d(5, term), &
                & g2d(6, term), term_energy, term_gradient(1:2))
            energy = energy+term_energy
            if( calculate_gradient )then
                do coordinate=1, 2
                    distance_index           = ig2d(coordinate, term)
                    gradient(distance_index) = gradient(distance_index)+term_gradient(coordinate)
                enddo
            endif
        enddo
        do term=1, size(ig3d, 2)
            do coordinate=1, 3
                term_distances(coordinate) = distances(ig3d(coordinate, term))
            enddo
            call evaluate_gaussian(term_distances, g3d(1:3, term), g3d(4:6, term), g3d(7, term), &
                & g3d(8, term), term_energy, term_gradient)
            energy = energy+term_energy
            if( calculate_gradient )then
                do coordinate=1, 3
                    distance_index           = ig3d(coordinate, term)
                    gradient(distance_index) = gradient(distance_index)+term_gradient(coordinate)
                enddo
            endif
        enddo
    end subroutine evaluate_fit

    !> @brief Evaluate one released Morse contribution and its distance derivative.
    pure subroutine evaluate_morse( parameters, distance, energy, derivative )
        real(wp), intent(in)  :: parameters(3), distance
        real(wp), intent(out) :: energy, derivative
        real(wp)              :: exponential
        exponential = exp(parameters(2)*(parameters(1)-distance))
        energy      = (parameters(3)*(exponential-1.0_wp))**2
        derivative  = -2.0_wp*parameters(2)*parameters(3)**2*exponential*(exponential-1.0_wp)
    end subroutine evaluate_morse

    !> @brief Evaluate one multidimensional Gaussian contribution and its distance derivatives.
    pure subroutine evaluate_gaussian( distances, centers, exponents, depth, offset, energy, derivative )
        real(wp), intent(in)  :: distances(:), centers(:), exponents(:), depth, offset
        real(wp), intent(out) :: energy, derivative(size(distances))
        real(wp)              :: exponent_sum, gaussian, separation(size(distances))
        integer               :: coordinate
        separation   = distances-centers
        exponent_sum = 0.0_wp
        do coordinate=1, size(distances)
            exponent_sum = exponent_sum+separation(coordinate)**2*exponents(coordinate)
        enddo
        gaussian   = exp(-0.5_wp*exponent_sum)
        energy     = depth*(gaussian-offset)
        derivative = -depth*gaussian*separation*exponents
    end subroutine evaluate_gaussian
end module chempotpy_malonaldehyde_mht_2014_pes

!> @brief Expose malonaldehyde_MHT_2014 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_kinds,                      only: wp
    use chempotpy_malonaldehyde_mht_2014_pes, only: evaluate
    implicit none
    integer, parameter    :: natoms = 9, nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
