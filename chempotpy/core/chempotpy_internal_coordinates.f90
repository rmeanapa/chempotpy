!> @brief Adapt internal-coordinate PES kernels to Cartesian CHEMPOTPY interfaces.
module chempotpy_internal_coordinates
    use chempotpy_geometry,     only: transform_internal_gradient_to_cartesian, triatomic_distance_jacobian
    use chempotpy_kinds,        only: real64
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_gradient_to_electronvolt_per_angstrom, &
        & legacy_hartree_to_electronvolt
    implicit none
    private
    public :: finalize_legacy_triatomic_surface
    public :: prepare_legacy_triatomic_coordinates

    contains

    !> @brief Build historical triatomic distances in bohr and their Cartesian Jacobian.
    pure subroutine prepare_legacy_triatomic_coordinates( coordinates_angstrom, distances_bohr, distance_jacobian )
        real(real64), intent(in)  :: coordinates_angstrom(3,3)
        real(real64), intent(out) :: distances_bohr(3), distance_jacobian(3,9)
        real(real64)              :: distances_angstrom(3)
        call triatomic_distance_jacobian(coordinates_angstrom, distances_angstrom, distance_jacobian)
        distances_bohr = distances_angstrom/legacy_bohr_to_angstrom
    end subroutine prepare_legacy_triatomic_coordinates

    !> @brief Convert and pack one triatomic internal-coordinate result for the historical PES contract.
    subroutine finalize_legacy_triatomic_surface( derivative_order, energy_hartree, internal_gradient, &
            distance_jacobian, energy_electronvolt, gradient_electronvolt_per_angstrom, coupling_per_angstrom )
        integer, intent(in)       :: derivative_order
        real(real64), intent(in)  :: energy_hartree, internal_gradient(3), distance_jacobian(3,9)
        real(real64), intent(out) :: energy_electronvolt(1)
        real(real64), intent(out) :: gradient_electronvolt_per_angstrom(1,3,3)
        real(real64), intent(out) :: coupling_per_angstrom(1,1,3,3)
        real(real64)              :: cartesian_gradient(9)
        gradient_electronvolt_per_angstrom = 0.0_real64
        coupling_per_angstrom              = 0.0_real64
        energy_electronvolt(1)             = energy_hartree*legacy_hartree_to_electronvolt
        if( derivative_order == 1 )then
            call transform_internal_gradient_to_cartesian( &
                internal_gradient*legacy_gradient_to_electronvolt_per_angstrom, distance_jacobian, cartesian_gradient)
            gradient_electronvolt_per_angstrom(1,:,:) = transpose(reshape(cartesian_gradient, [3, 3]))
        elseif( derivative_order == 2 )then
            write(*, *) "Only energy and gradient are available"
        endif
    end subroutine finalize_legacy_triatomic_surface
end module chempotpy_internal_coordinates
