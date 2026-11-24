! System:                        HeH2p
! Functional form:               diatomics-in-molecules
! Common name:                   He + H2+ DIM
! Reference:                     W. N. Whitton and P. J. Kurtz, J. Chem. Phys. 64, 3624 (1976)
! Number of derivatives:         1
! Number of bodies:              3
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable parameters are embedded
! Notes:                         The atom order is H, H, He. Energy is relative to separated He + H2+.

!> @brief Evaluate the Whitton-Kurtz diatomics-in-molecules surface for He + H2+.
module chempotpy_heh2p_dim_1976_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    real(wp), parameter :: asymptotic_energy        = -2.8143_wp
    real(wp), parameter :: underflow_limit          = 1.0e-60_wp
    real(wp), parameter :: dissociation_energy(2,3) = reshape([ &
        2.8141_wp, 2.9014_wp, 2.0389_wp, 1.6431_wp, 2.0389_wp, 1.6431_wp], [2, 3])
    real(wp), parameter :: equilibrium_distance(2,3) = reshape([ &
        1.0490_wp, 1.0680_wp, 0.7761_wp, 0.7226_wp, 0.7761_wp, 0.7226_wp], [2, 3])
    real(wp), parameter :: exponential_scale(2,3) = reshape([ &
        1.4049_wp, 1.4968_wp, 2.7936_wp, 3.0664_wp, 2.7936_wp, 3.0664_wp], [2, 3])
    real(wp), parameter :: linear_scale(2,3) = reshape([ &
        -0.1374_wp, -0.0573_wp, 0.0656_wp, -0.0227_wp, 0.0656_wp, -0.0227_wp], [2, 3])
    real(wp), parameter :: quadratic_scale(2,3) = reshape([ &
        0.0514_wp, 0.0015_wp, 0.0453_wp, 0.0681_wp, 0.0453_wp, 0.0681_wp], [2, 3])
    real(wp), parameter :: exponential_sign(2) = [-1.0_wp, 1.0_wp]

    contains

    !> @brief Evaluate energy and optional Cartesian gradient in CHEMPOTPY units.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_geometry, only: ordered_pair_distance_jacobian, &
            & transform_internal_gradient_to_cartesian
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms          = 3
        integer, parameter    :: nstates         = 1
        integer, parameter    :: first_atoms(3)  = [1, 2, 1]
        integer, parameter    :: second_atoms(3) = [2, 3, 3]
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: distance_jacobian(3,9), distances(3), internal_gradient(1,3)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad < 0 .or. igrad > 1 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'HeH2p_DIM_1976 supports energy and first derivatives')
            return
        endif
        call pair_distances(x, distances)
        call evaluate_internal(distances, p(1), internal_gradient(1,:))
        if( igrad == 1 )then
            call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, distances, distance_jacobian)
            call transform_internal_gradient_to_cartesian(internal_gradient, distance_jacobian, g)
        endif
    end subroutine evaluate

    !> @brief Return H-H, H-He, and H-He distances in the published ordering.
    pure subroutine pair_distances( coordinates, distances )
        use chempotpy_geometry, only: cartesian_distance
        real(wp), intent(in)  :: coordinates(3,3)
        real(wp), intent(out) :: distances(3)
        distances(1) = cartesian_distance(coordinates(1,:), coordinates(2,:))
        distances(2) = cartesian_distance(coordinates(2,:), coordinates(3,:))
        distances(3) = cartesian_distance(coordinates(1,:), coordinates(3,:))
    end subroutine pair_distances

    !> @brief Evaluate the DIM fit in its native eV and angstrom units.
    pure subroutine evaluate_internal( distances, energy, derivative )
        real(wp), intent(in)  :: distances(3)
        real(wp), intent(out) :: energy, derivative(3)
        real(wp)              :: coupling(3), coupling_derivative(3), coupling_norm
        real(wp)              :: first_derivative, first_energy, second_derivative, second_energy
        integer               :: pair
        energy = 0.0_wp
        do pair = 1, 3
            call diatomic_term(1, pair, distances(pair), first_energy, first_derivative)
            call diatomic_term(2, pair, distances(pair), second_energy, second_derivative)
            energy                    = energy+0.5_wp*(first_energy+second_energy)
            derivative(pair)          = 0.5_wp*(first_derivative+second_derivative)
            coupling(pair)            = 0.5_wp*(first_energy-second_energy)
            coupling_derivative(pair) = 0.5_wp*(first_derivative-second_derivative)
        enddo
        coupling_norm = max(sqrt(coupling(1)**2+(coupling(2)-coupling(3))**2), underflow_limit)
        energy        = energy-coupling_norm-asymptotic_energy
        derivative(1) = derivative(1)-coupling(1)*coupling_derivative(1)/coupling_norm
        derivative(2) = derivative(2)-(coupling(2)-coupling(3))*coupling_derivative(2)/coupling_norm
        derivative(3) = derivative(3)+(coupling(2)-coupling(3))*coupling_derivative(3)/coupling_norm
    end subroutine evaluate_internal

    !> @brief Evaluate one published extended-Morse diatomic term and its radial derivative.
    pure subroutine diatomic_term( state, pair, distance, energy, derivative )
        integer, intent(in)   :: state, pair
        real(wp), intent(in)  :: distance
        real(wp), intent(out) :: energy, derivative
        real(wp)              :: displacement, exponent, exponent_derivative, exponential
        displacement = distance-equilibrium_distance(state,pair)
        exponent            = exponential_scale(state,pair)*(1.0_wp+linear_scale(state,pair)*displacement + &
            & quadratic_scale(state,pair)*displacement**2)
        exponent_derivative = exponential_scale(state,pair)*(linear_scale(state,pair) + 2.0_wp*quadratic_scale(state, &
            & pair)*displacement)
        exponential = exponential_sign(state)*exp(-exponent*displacement)
        energy      = dissociation_energy(state,pair)*(1.0_wp+exponential)**2 - dissociation_energy(state,pair)
        derivative          = -2.0_wp*dissociation_energy(state,pair)*(1.0_wp+exponential)*exponential* (exponent+ &
            & displacement*exponent_derivative)
    end subroutine diatomic_term
end module chempotpy_heh2p_dim_1976_pes

!> @brief Expose HeH2p_DIM_1976 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_heh2p_dim_1976_pes, only: evaluate
    use chempotpy_kinds,              only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
