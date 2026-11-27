!> @brief Provide reusable dense neural-network evaluation utilities for PES fits.
module chempotpy_neural_network
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: activation_linear
    public :: activation_logistic
    public :: activation_tanh
    public :: apply_activation
    public :: evaluate_three_hidden_tanh_network
    public :: evaluate_two_hidden_network
    integer, parameter :: activation_tanh     = 1
    integer, parameter :: activation_logistic = 2
    integer, parameter :: activation_linear   = 3

    contains

    !> @brief Apply the selected scalar activation function.
    pure elemental function apply_activation( argument, activation ) result(value)
        real(wp), intent(in) :: argument
        integer, intent(in)  :: activation
        real(wp)             :: value
        select case( activation )
        case( activation_tanh )
            value = tanh(argument)
        case( activation_logistic )
            value = 1.0_wp/(1.0_wp+exp(-argument))
        case default
            value = argument
        end select
    end function apply_activation

    !> @brief Evaluate a three-hidden-layer tanh network and its full input Jacobian.
    subroutine evaluate_three_hidden_tanh_network( inputs, weights_1, biases_1, weights_2, biases_2, &
            & weights_3, biases_3, output_weights, output_biases, outputs, output_jacobian )
        real(wp), intent(in)  :: inputs(:), weights_1(:,:), biases_1(:)
        real(wp), intent(in)  :: weights_2(:,:), biases_2(:), weights_3(:,:), biases_3(:)
        real(wp), intent(in)  :: output_weights(:,:), output_biases(:)
        real(wp), intent(out) :: outputs(size(output_biases))
        real(wp), intent(out) :: output_jacobian(size(output_biases),size(inputs))
        real(wp)              :: hidden_1(size(biases_1)), hidden_2(size(biases_2)), hidden_3(size(biases_3))
        real(wp)              :: jacobian_1(size(biases_1),size(inputs))
        real(wp)              :: jacobian_2(size(biases_2),size(biases_1))
        real(wp)              :: jacobian_3(size(biases_3),size(biases_2))
        real(wp)              :: output_layer_jacobian(size(output_biases),size(biases_3))
        real(wp)              :: combined_2(size(biases_2),size(inputs))
        real(wp)              :: combined_3(size(biases_3),size(inputs))
        integer               :: input_index, node_index
        hidden_1 = tanh(matmul(inputs, weights_1)+biases_1)
        hidden_2 = tanh(matmul(hidden_1, weights_2)+biases_2)
        hidden_3 = tanh(matmul(hidden_2, weights_3)+biases_3)
        outputs  = matmul(hidden_3, output_weights)+output_biases
        do input_index = 1, size(inputs)
            jacobian_1(:,input_index) = (1.0_wp-hidden_1**2)*weights_1(input_index,:)
        enddo
        do input_index = 1, size(hidden_1)
            jacobian_2(:,input_index) = (1.0_wp-hidden_2**2)*weights_2(input_index,:)
        enddo
        do input_index = 1, size(hidden_2)
            jacobian_3(:,input_index) = (1.0_wp-hidden_3**2)*weights_3(input_index,:)
        enddo
        do node_index = 1, size(hidden_3)
            output_layer_jacobian(:,node_index) = output_weights(node_index,:)
        enddo
        combined_2      = matmul(jacobian_2, jacobian_1)
        combined_3      = matmul(jacobian_3, combined_2)
        output_jacobian = matmul(output_layer_jacobian, combined_3)
    end subroutine evaluate_three_hidden_tanh_network

    !> @brief Evaluate a scaled dense network with two hidden layers and one linear output.
    subroutine evaluate_two_hidden_network( inputs, input_averages, input_deltas, input_hidden_weights, &
            & hidden_1_biases, hidden_hidden_weights, hidden_2_biases, hidden_output_weights, output_bias, &
            & activation, output_average, output_delta, output )
        real(wp), intent(in)  :: inputs(:), input_averages(:), input_deltas(:)
        real(wp), intent(in)  :: input_hidden_weights(:,:), hidden_1_biases(:)
        real(wp), intent(in)  :: hidden_hidden_weights(:,:), hidden_2_biases(:)
        real(wp), intent(in)  :: hidden_output_weights(:), output_bias
        integer, intent(in)   :: activation
        real(wp), intent(in)  :: output_average, output_delta
        real(wp), intent(out) :: output
        integer               :: input_index, node_index, previous_index
        real(wp)              :: input_values(size(inputs))
        real(wp)              :: hidden_1_values(size(hidden_1_biases))
        real(wp)              :: hidden_2_values(size(hidden_2_biases)), output_value
        do input_index=1, size(inputs)
            input_values(input_index) = (inputs(input_index)-input_averages(input_index))/input_deltas(input_index)
        enddo
        do node_index=1, size(hidden_1_biases)
            hidden_1_values(node_index) = hidden_1_biases(node_index)
            do previous_index=1, size(inputs)
                hidden_1_values(node_index) = hidden_1_values(node_index) + input_values(previous_index)* &
                    & input_hidden_weights(previous_index,node_index)
            enddo
            hidden_1_values(node_index) = apply_activation(hidden_1_values(node_index), activation)
        enddo
        do node_index=1, size(hidden_2_biases)
            hidden_2_values(node_index) = hidden_2_biases(node_index)
            do previous_index=1, size(hidden_1_biases)
                hidden_2_values(node_index) = hidden_2_values(node_index) + hidden_1_values(previous_index)* &
                    & hidden_hidden_weights(previous_index,node_index)
            enddo
            hidden_2_values(node_index) = apply_activation(hidden_2_values(node_index), activation)
        enddo
        output_value = output_bias
        do previous_index=1, size(hidden_2_biases)
            output_value = output_value+hidden_2_values(previous_index)*hidden_output_weights(previous_index)
        enddo
        output = output_value*output_delta+output_average
    end subroutine evaluate_two_hidden_network
end module chempotpy_neural_network
