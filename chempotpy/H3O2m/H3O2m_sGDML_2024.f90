! System:                        H3O2m
! Functional form:               symmetric gradient-domain machine learning (sGDML)
! Common name:                   H3O2- sGDML
! Reference:                     M. Arandhara and S. G. Ramesh, J. Phys. Chem. A 128, 1600-1610 (2024),
!                                doi:10.1021/acs.jpca.3c08027
! Number of derivatives:         1
! Number of bodies:              5
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; the released 3,000-point model is embedded
! Notes:                         The public atom order is O, O, H, H, H. The fit uses 10 inverse-distance
!                                descriptors and all 12 permutations of like atoms.

!> @brief Implement the 2024 H3O2- symmetric gradient-domain machine-learning surface.
module chempotpy_h3o2m_sgdml_2024_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    integer, parameter :: distance_order(10) = [1, 2, 5, 3, 6, 8, 4, 7, 9, 10]

    contains

    !> @brief Evaluate the H3O2- sGDML energy and optional analytic gradient.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_geometry, only: all_pair_distances, all_pair_distances_and_jacobian, &
            reshape_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_h3o2m_sgdml_parameters, only: descriptor_coefficients, descriptor_permutations, &
            energy_offset, energy_scale, kernel_sigma, training_descriptors
        use chempotpy_sgdml,  only: evaluate_matern52_sgdml
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: coordinates_to_bohr, hartree_per_bohr_to_ev_per_angstrom, &
            hartree_to_electronvolt
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms  = 5
        integer, parameter    :: nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
        real(wp)              :: coordinates_bohr(natoms,3), distance_jacobian(10,15), distances(10)
        real(wp)              :: descriptor_gradient(10), descriptor_jacobian(10,15), descriptors(10)
        real(wp)              :: energy_hartree, gradient_flat(15), gradient_hartree(natoms,3)
        integer               :: descriptor
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( present(status) ) call status%clear()
        if( igrad < 0.or.igrad > 1 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'H3O2m_sGDML_2024 provides energies and first derivatives')
            return
        endif
        call coordinates_to_bohr(x, coordinates_bohr)
        call all_pair_distances(coordinates_bohr, distances)
        if( any(distances <= sqrt(tiny(1.0_wp))) )then
            if( present(status) ) call status%set(status_invalid_input, &
                'H3O2m_sGDML_2024 requires distinct atomic positions')
            return
        endif
        descriptors = 1.0_wp/distances(distance_order)
        call evaluate_matern52_sgdml(descriptors, training_descriptors, descriptor_coefficients, &
            descriptor_permutations, kernel_sigma, energy_scale, energy_offset, igrad == 1, energy_hartree, &
            descriptor_gradient)
        p(1) = energy_hartree*hartree_to_electronvolt
        if( igrad == 1 )then
            call all_pair_distances_and_jacobian(coordinates_bohr, distances, distance_jacobian)
            do descriptor = 1, 10
                descriptor_jacobian(descriptor,:) = -distance_jacobian(distance_order(descriptor),:) / &
                    & distances(distance_order(descriptor))**2
            enddo
            call transform_internal_gradient_to_cartesian(descriptor_gradient, descriptor_jacobian, gradient_flat)
            call reshape_cartesian_coordinates(gradient_flat, gradient_hartree)
            g(1,:,:) = gradient_hartree*hartree_per_bohr_to_ev_per_angstrom
        endif
    end subroutine evaluate
end module chempotpy_h3o2m_sgdml_2024_pes

!> @brief Expose the H3O2- sGDML surface through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3o2m_sgdml_2024_pes, only: evaluate
    use chempotpy_kinds,                only: wp
    implicit none
    integer, parameter    :: natoms  = 5
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
