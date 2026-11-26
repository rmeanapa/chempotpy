!> @brief Evaluate symmetric gradient-domain machine-learning models.
module chempotpy_sgdml
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_matern52_sgdml

    contains

    !> @brief Evaluate an sGDML model using its permuted Matern-5/2 descriptor kernel.
    pure subroutine evaluate_matern52_sgdml( query_descriptor, training_descriptors, descriptor_coefficients, &
            descriptor_permutations, kernel_sigma, energy_scale, energy_offset, calculate_gradient, energy, &
            descriptor_gradient )
        real(wp), intent(in)  :: query_descriptor(:), training_descriptors(:,:), descriptor_coefficients(:,:)
        integer, intent(in)   :: descriptor_permutations(:,:)
        real(wp), intent(in)  :: kernel_sigma, energy_scale, energy_offset
        logical, intent(in)   :: calculate_gradient
        real(wp), intent(out) :: energy, descriptor_gradient(size(query_descriptor))
        real(wp), allocatable :: energy_kernel(:), projection(:), weighted_projection(:)
        real(wp)              :: coefficient(size(query_descriptor)), difference(size(query_descriptor))
        real(wp)              :: first_sum, kernel_value, norm_difference, second_sum
        real(wp)              :: base_factor, diagonal_factor, sqrt_five
        integer               :: descriptor, model_index, permutation, training
        sqrt_five       = sqrt(5.0_wp)
        base_factor     = 5.0_wp/(3.0_wp*kernel_sigma**3)
        diagonal_factor = 5.0_wp/kernel_sigma
        allocate(energy_kernel(size(training_descriptors, 2)*size(descriptor_permutations, 2)))
        allocate(projection(size(energy_kernel)), weighted_projection(size(energy_kernel)))
        model_index = 0
        do training = 1, size(training_descriptors, 2)
            do permutation = 1, size(descriptor_permutations, 2)
                model_index = model_index+1
                do descriptor = 1, size(query_descriptor)
                    difference(descriptor) = query_descriptor(descriptor) - &
                        & training_descriptors( descriptor_permutations(descriptor,permutation),training)
                    coefficient(descriptor) = descriptor_coefficients( descriptor_permutations(descriptor, &
                        & permutation),training)
                enddo
                norm_difference                  = sqrt_five*sqrt(sum(difference*difference))
                kernel_value                     = exp(-norm_difference/kernel_sigma)*base_factor
                projection(model_index)          = dot_product(difference, coefficient)
                energy_kernel(model_index)       = kernel_value*(norm_difference+kernel_sigma)
                weighted_projection(model_index) = projection(model_index)*kernel_value
            enddo
        enddo
        energy              = dot_product(projection, energy_kernel)*energy_scale+energy_offset
        descriptor_gradient = 0.0_wp
        if( calculate_gradient )then
            do descriptor = 1, size(query_descriptor)
                first_sum   = 0.0_wp
                second_sum  = 0.0_wp
                model_index = 0
                do training = 1, size(training_descriptors, 2)
                    do permutation = 1, size(descriptor_permutations, 2)
                        model_index = model_index+1
                        first_sum = first_sum+weighted_projection(model_index)*(query_descriptor(descriptor) - &
                            & training_descriptors(descriptor_permutations(descriptor,permutation),training))
                        second_sum = second_sum+energy_kernel(model_index)* &
                            & descriptor_coefficients( descriptor_permutations(descriptor,permutation),training)
                    enddo
                enddo
                descriptor_gradient(descriptor) = energy_scale*(second_sum-diagonal_factor*first_sum)
            enddo
        endif
    end subroutine evaluate_matern52_sgdml
end module chempotpy_sgdml
