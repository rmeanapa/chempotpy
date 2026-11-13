!> @brief Adapt internal-coordinate PES kernels to Cartesian CHEMPOTPY interfaces.
module chempotpy_internal_coordinates
    use chempotpy_geometry, only: transform_internal_gradient_to_cartesian, triatomic_distance_jacobian
    use chempotpy_kinds,    only: wp
    use chempotpy_units,    only: bohr_to_angstrom, hartree_per_bohr_to_ev_per_angstrom, &
        & hartree_to_electronvolt
    implicit none
    private
    public :: finalize_legacy_triatomic_surface
    public :: prepare_legacy_triatomic_coordinates

    contains

    !> @brief Build historical triatomic distances in bohr and their Cartesian Jacobian.
    pure subroutine prepare_legacy_triatomic_coordinates( coordinates_angstrom, distances_bohr, distance_jacobian )
        real(wp), intent(in)  :: coordinates_angstrom(3,3)
        real(wp), intent(out) :: distances_bohr(3), distance_jacobian(3,9)
        real(wp)              :: distances_angstrom(3)
        call triatomic_distance_jacobian(coordinates_angstrom, distances_angstrom, distance_jacobian)
        distances_bohr = distances_angstrom/bohr_to_angstrom
    end subroutine prepare_legacy_triatomic_coordinates

    !> @brief Convert and pack one triatomic internal-coordinate result for the historical PES contract.
    subroutine finalize_legacy_triatomic_surface( derivative_order, energy_hartree, internal_gradient, &
            distance_jacobian, energy_electronvolt, gradient_electronvolt_per_angstrom, coupling_per_angstrom )
        integer, intent(in)   :: derivative_order
        real(wp), intent(in)  :: energy_hartree, internal_gradient(3), distance_jacobian(3,9)
        real(wp), intent(out) :: energy_electronvolt(1)
        real(wp), intent(out) :: gradient_electronvolt_per_angstrom(1,3,3)
        real(wp), intent(out) :: coupling_per_angstrom(1,1,3,3)
        real(wp)              :: cartesian_gradient(9)
        gradient_electronvolt_per_angstrom = 0.0_wp
        coupling_per_angstrom              = 0.0_wp
        energy_electronvolt(1)             = energy_hartree*hartree_to_electronvolt
        if( derivative_order == 1 )then
            call transform_internal_gradient_to_cartesian( &
                internal_gradient*hartree_per_bohr_to_ev_per_angstrom, distance_jacobian, cartesian_gradient)
            gradient_electronvolt_per_angstrom(1,:,:) = transpose(reshape(cartesian_gradient, [3, 3]))
        elseif( derivative_order == 2 )then
            write(*, *) "Only energy and gradient are available"
        endif
    end subroutine finalize_legacy_triatomic_surface
end module chempotpy_internal_coordinates
