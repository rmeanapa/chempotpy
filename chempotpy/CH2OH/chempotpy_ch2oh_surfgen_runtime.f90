!> @brief Evaluate immutable two- and three-state CH2OH Surfgen Hamiltonians.
module chempotpy_ch2oh_surfgen_runtime
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_ch2oh_adiabatic
    public :: evaluate_ch2oh_diabatic
    real(wp), parameter :: finite_difference_step_angstrom = 1.0e-5_wp

    contains

    !> @brief Evaluate adiabatic energies, gradients, and derivative couplings for one contracted model.
    subroutine evaluate_ch2oh_adiabatic( x, derivative_order, coefficient_order, coefficient_row, &
            & coefficient_column, coefficient_term, coefficient_value, energy, gradient, coupling, status )
        use chempotpy_linear_algebra, only: symmetric_eigendecomposition
        use chempotpy_status,         only: chempotpy_status_t, status_derivative_unavailable, status_numerical_error
        use chempotpy_units,          only: hartree_to_electronvolt
        type(chempotpy_status_t), intent(out), optional :: status
        integer, intent(in)   :: derivative_order
        integer, intent(in)   :: coefficient_order(:), coefficient_row(:), coefficient_column(:), coefficient_term(:)
        real(wp), intent(in)  :: x(:,:), coefficient_value(:)
        real(wp), intent(out) :: energy(:), gradient(:,:,:), coupling(:,:,:,:)
        real(wp)              :: diabatic_gradient(size(x, 1),size(x, 2),size(energy),size(energy))
        real(wp)              :: eigenvectors(size(energy),size(energy)), hamiltonian(size(energy),size(energy))
        real(wp)              :: projected(size(energy),size(energy)), gap
        integer               :: atom, column, first_state, info, second_state, state
        if( present(status) ) call status%clear()
        energy   = 0.0_wp
        gradient = 0.0_wp
        coupling = 0.0_wp
        call evaluate_model_hamiltonian(x, coefficient_order, coefficient_row, coefficient_column, &
            & coefficient_term, coefficient_value, hamiltonian, info)
        if( info /= 0 )then
            if( present(status) ) call status%set(status_numerical_error, &
                & 'CH2OH Surfgen coordinate evaluation failed')
            return
        endif
        call symmetric_eigendecomposition(hamiltonian, energy, eigenvectors, info)
        if( info /= 0 )then
            if( present(status) ) call status%set(status_numerical_error, &
                & 'CH2OH Surfgen Hamiltonian diagonalization failed')
            return
        endif
        do state = 1, size(energy)
            if( eigenvectors(1,state) > 0.0_wp ) eigenvectors(:,state) = -eigenvectors(:,state)
        enddo
        if( derivative_order > 1 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                & 'CH2OH Surfgen provides energies and first derivatives')
            return
        endif
        if( derivative_order == 1 )then
            call finite_difference_diabatic_gradient(x, coefficient_order, coefficient_row, coefficient_column, &
                & coefficient_term, coefficient_value, diabatic_gradient, info)
            if( info /= 0 )then
                if( present(status) ) call status%set(status_numerical_error, &
                    & 'CH2OH Surfgen finite-difference evaluation failed')
                return
            endif
            do atom = 1, size(x, 1)
                do column = 1, size(x, 2)
                    projected = matmul(transpose(eigenvectors), matmul(diabatic_gradient(atom,column,:,:), &
                        & eigenvectors))
                    do first_state = 1, size(energy)
                        gradient(first_state,atom,column) = projected(first_state,first_state) *hartree_to_electronvolt
                        do second_state = first_state+1, size(energy)
                            gap = energy(second_state)-energy(first_state)
                            if( abs(gap) > tiny(gap) )then
                                coupling(first_state,second_state,atom,column) = projected(first_state,second_state)/gap
                                coupling(second_state,first_state,atom,column) = coupling(first_state,second_state, &
                                    & atom,column)
                            endif
                        enddo
                    enddo
                enddo
            enddo
        endif
        energy = energy*hartree_to_electronvolt
    end subroutine evaluate_ch2oh_adiabatic

    !> @brief Evaluate the diabatic Hamiltonian and its Cartesian gradient for one contracted model.
    subroutine evaluate_ch2oh_diabatic( x, derivative_order, coefficient_order, coefficient_row, &
            & coefficient_column, coefficient_term, coefficient_value, hamiltonian, hamiltonian_gradient, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_numerical_error
        use chempotpy_units,  only: hartree_to_electronvolt
        type(chempotpy_status_t), intent(out), optional :: status
        integer, intent(in)   :: derivative_order
        integer, intent(in)   :: coefficient_order(:), coefficient_row(:), coefficient_column(:), coefficient_term(:)
        real(wp), intent(in)  :: x(:,:), coefficient_value(:)
        real(wp), intent(out) :: hamiltonian(:,:), hamiltonian_gradient(:,:,:,:)
        real(wp)              :: cartesian_gradient(size(x, 1),size(x, 2),size(hamiltonian, 1),size(hamiltonian, 2))
        integer               :: atom, column, info
        if( present(status) ) call status%clear()
        hamiltonian          = 0.0_wp
        hamiltonian_gradient = 0.0_wp
        call evaluate_model_hamiltonian(x, coefficient_order, coefficient_row, coefficient_column, &
            & coefficient_term, coefficient_value, hamiltonian, info)
        if( info /= 0 )then
            if( present(status) ) call status%set(status_numerical_error, &
                & 'CH2OH Surfgen coordinate evaluation failed')
            return
        endif
        if( derivative_order > 1 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                & 'CH2OH Surfgen provides energies and first derivatives')
            return
        endif
        if( derivative_order == 1 )then
            call finite_difference_diabatic_gradient(x, coefficient_order, coefficient_row, coefficient_column, &
                & coefficient_term, coefficient_value, cartesian_gradient, info)
            if( info /= 0 )then
                if( present(status) ) call status%set(status_numerical_error, &
                    & 'CH2OH Surfgen finite-difference evaluation failed')
                return
            endif
            do atom = 1, size(x, 1)
                do column = 1, size(x, 2)
                    hamiltonian_gradient(:,:,atom,column) = cartesian_gradient(atom,column,:,:)
                enddo
            enddo
        endif
        hamiltonian          = hamiltonian*hartree_to_electronvolt
        hamiltonian_gradient = hamiltonian_gradient*hartree_to_electronvolt
    end subroutine evaluate_ch2oh_diabatic

    !> @brief Evaluate a contracted CH2OH quasi-diabatic Hamiltonian in hartree.
    subroutine evaluate_model_hamiltonian( x, coefficient_order, coefficient_row, coefficient_column, &
            & coefficient_term, coefficient_value, hamiltonian, info )
        use chempotpy_ch2oh_surfgen_basis, only: atom_indices, coordinate_parameters, coordinate_types, &
            & scaling_types, second_order_direction, second_order_divisor, second_order_parent, &
            & third_order_direction, third_order_divisor, third_order_parent
        use chempotpy_surfgen, only: evaluate_cubic_surfgen_hamiltonian, evaluate_surfgen_coordinates
        use chempotpy_units,   only: angstrom_to_bohr
        integer, intent(in)   :: coefficient_order(:), coefficient_row(:), coefficient_column(:), coefficient_term(:)
        real(wp), intent(in)  :: x(:,:), coefficient_value(:)
        real(wp), intent(out) :: hamiltonian(:,:)
        integer, intent(out)  :: info
        real(wp)              :: coordinates(49), geometry_bohr(size(x, 1),size(x, 2))
        geometry_bohr = x*angstrom_to_bohr
        call evaluate_surfgen_coordinates(geometry_bohr, coordinate_types, scaling_types, atom_indices, &
            & coordinate_parameters, coordinates, info)
        if( info /= 0 )then
            hamiltonian = 0.0_wp
            return
        endif
        call evaluate_cubic_surfgen_hamiltonian(coordinates, second_order_parent, second_order_direction, &
            & second_order_divisor, third_order_parent, third_order_direction, third_order_divisor, &
            & coefficient_order, coefficient_row, coefficient_column, coefficient_term, coefficient_value, hamiltonian)
    end subroutine evaluate_model_hamiltonian

    !> @brief Differentiate the diabatic Hamiltonian with centered Cartesian differences.
    subroutine finite_difference_diabatic_gradient( x, coefficient_order, coefficient_row, coefficient_column, &
            & coefficient_term, coefficient_value, gradient, info )
        integer, intent(in)   :: coefficient_order(:), coefficient_row(:), coefficient_column(:), coefficient_term(:)
        real(wp), intent(in)  :: x(:,:), coefficient_value(:)
        real(wp), intent(out) :: gradient(:,:,:,:)
        integer, intent(out)  :: info
        real(wp)              :: backward(size(gradient, 3),size(gradient, 4))
        real(wp)              :: forward(size(gradient, 3),size(gradient, 4)), shifted(size(x, 1),size(x, 2))
        integer               :: atom, column
        info     = 0
        gradient = 0.0_wp
        do atom = 1, size(x, 1)
            do column = 1, size(x, 2)
                shifted              = x
                shifted(atom,column) = shifted(atom,column)+finite_difference_step_angstrom
                call evaluate_model_hamiltonian(shifted, coefficient_order, coefficient_row, coefficient_column, &
                    & coefficient_term, coefficient_value, forward, info)
                if( info /= 0 ) return
                shifted(atom,column) = shifted(atom,column)-2.0_wp*finite_difference_step_angstrom
                call evaluate_model_hamiltonian(shifted, coefficient_order, coefficient_row, coefficient_column, &
                    & coefficient_term, coefficient_value, backward, info)
                if( info /= 0 ) return
                gradient(atom,column,:,:) = (forward-backward)/(2.0_wp*finite_difference_step_angstrom)
            enddo
        enddo
    end subroutine finite_difference_diabatic_gradient
end module chempotpy_ch2oh_surfgen_runtime
