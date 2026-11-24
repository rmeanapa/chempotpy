!> @brief Provide reusable centered finite-difference validation for Cartesian PES evaluators.
module chempotpy_finite_difference
    use chempotpy_kinds,               only: wp
    use chempotpy_status,              only: chempotpy_status_t, status_invalid_input, status_numerical_error
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: central_finite_difference
    public :: central_finite_difference_gradient
    public :: pes_value_evaluator
    abstract interface
        !> @brief Evaluate one Cartesian geometry and return a flat vector of PES values.
        subroutine pes_value_evaluator( coordinates, values, status )
            import :: chempotpy_status_t, wp
            real(wp), intent(in)  :: coordinates(:,:)
            real(wp), intent(out) :: values(:)
            type(chempotpy_status_t), intent(inout) :: status
        end subroutine pes_value_evaluator
    end interface

    contains

    !> @brief Evaluate values and Cartesian gradients using centered finite differences.
    subroutine central_finite_difference_gradient( evaluator, coordinates, step, values, gradient, status )
        procedure(pes_value_evaluator)              :: evaluator
        real(wp), intent(in)  :: coordinates(:,:), step
        real(wp), intent(out) :: values(:), gradient(:,:)
        type(chempotpy_status_t), intent(inout)     :: status
        real(wp) :: backward(size(values)), forward(size(values))
        real(wp) :: displaced(size(coordinates, 1),size(coordinates, 2))
        integer  :: atom, coordinate, coordinate_count, direction
        values = 0.0_wp; gradient = 0.0_wp
        if( status%failed() ) return
        coordinate_count = 3*size(coordinates, 1)
        if( size(coordinates, 2) /= 3 )then
            call status%set(status_invalid_input, 'finite-difference coordinates must have shape (natoms,3)')
        elseif( step <= 0.0_wp .or. .not. ieee_is_finite(step) )then
            call status%set(status_invalid_input, 'finite-difference step must be positive and finite')
        elseif( size(values) < 1 )then
            call status%set(status_invalid_input, 'finite-difference values must not be empty')
        elseif( any(shape(gradient) /= [size(values), coordinate_count]) )then
            call status%set(status_invalid_input, 'finite-difference gradient has the wrong shape')
        endif
        if( status%failed() ) return
        call evaluator(coordinates, values, status)
        call require_finite_values(values, status)
        if( status%failed() ) return
        do coordinate=1, coordinate_count
            atom                      = (coordinate-1)/3+1
            direction                 = mod(coordinate-1, 3)+1
            displaced                 = coordinates
            displaced(atom,direction) = displaced(atom,direction)+step
            call evaluator(displaced, forward, status)
            call require_finite_values(forward, status)
            if( status%failed() ) return
            displaced(atom,direction) = displaced(atom,direction)-2.0_wp*step
            call evaluator(displaced, backward, status)
            call require_finite_values(backward, status)
            if( status%failed() ) return
            gradient(:,coordinate) = (forward-backward)/(2.0_wp*step)
        enddo
    end subroutine central_finite_difference_gradient

    !> @brief Evaluate values, Cartesian gradients, and Hessians using centered finite differences.
    subroutine central_finite_difference( evaluator, coordinates, step, values, gradient, hessian, status )
        procedure(pes_value_evaluator)              :: evaluator
        real(wp), intent(in)  :: coordinates(:,:), step
        real(wp), intent(out) :: values(:), gradient(:,:), hessian(:,:,:)
        type(chempotpy_status_t), intent(inout)     :: status
        real(wp) :: backward(size(values)), forward(size(values))
        real(wp) :: mixed_mm(size(values)), mixed_mp(size(values))
        real(wp) :: mixed_pm(size(values)), mixed_pp(size(values))
        real(wp) :: displaced(size(coordinates, 1),size(coordinates, 2))
        integer  :: atom_i, atom_j, coordinate_count, coordinate_i, coordinate_j
        integer  :: direction_i, direction_j
        values = 0.0_wp; gradient = 0.0_wp; hessian = 0.0_wp
        if( status%failed() ) return
        coordinate_count = 3*size(coordinates, 1)
        if( size(coordinates, 2) /= 3 )then
            call status%set(status_invalid_input, 'finite-difference coordinates must have shape (natoms,3)')
        elseif( step <= 0.0_wp .or. .not. ieee_is_finite(step) )then
            call status%set(status_invalid_input, 'finite-difference step must be positive and finite')
        elseif( size(values) < 1 )then
            call status%set(status_invalid_input, 'finite-difference values must not be empty')
        elseif( any(shape(gradient) /= [size(values), coordinate_count]) )then
            call status%set(status_invalid_input, 'finite-difference gradient has the wrong shape')
        elseif( any(shape(hessian) /= [size(values), coordinate_count, coordinate_count]) )then
            call status%set(status_invalid_input, 'finite-difference Hessian has the wrong shape')
        endif
        if( status%failed() ) return
        call evaluator(coordinates, values, status)
        call require_finite_values(values, status)
        if( status%failed() ) return
        do coordinate_i = 1, coordinate_count
            atom_i                        = (coordinate_i-1)/3+1
            direction_i                   = mod(coordinate_i-1, 3)+1
            displaced                     = coordinates
            displaced(atom_i,direction_i) = displaced(atom_i,direction_i)+step
            call evaluator(displaced, forward, status)
            call require_finite_values(forward, status)
            if( status%failed() ) return
            displaced(atom_i,direction_i) = displaced(atom_i,direction_i)-2.0_wp*step
            call evaluator(displaced, backward, status)
            call require_finite_values(backward, status)
            if( status%failed() ) return
            gradient(:,coordinate_i)             = (forward-backward)/(2.0_wp*step)
            hessian(:,coordinate_i,coordinate_i) = (forward-2.0_wp*values+backward)/(step*step)
            do coordinate_j = coordinate_i+1, coordinate_count
                atom_j                        = (coordinate_j-1)/3+1
                direction_j                   = mod(coordinate_j-1, 3)+1
                displaced                     = coordinates
                displaced(atom_i,direction_i) = displaced(atom_i,direction_i)+step
                displaced(atom_j,direction_j) = displaced(atom_j,direction_j)+step
                call evaluator(displaced, mixed_pp, status)
                call require_finite_values(mixed_pp, status)
                if( status%failed() ) return
                displaced(atom_j,direction_j) = displaced(atom_j,direction_j)-2.0_wp*step
                call evaluator(displaced, mixed_pm, status)
                call require_finite_values(mixed_pm, status)
                if( status%failed() ) return
                displaced(atom_i,direction_i) = displaced(atom_i,direction_i)-2.0_wp*step
                call evaluator(displaced, mixed_mm, status)
                call require_finite_values(mixed_mm, status)
                if( status%failed() ) return
                displaced(atom_j,direction_j) = displaced(atom_j,direction_j)+2.0_wp*step
                call evaluator(displaced, mixed_mp, status)
                call require_finite_values(mixed_mp, status)
                if( status%failed() ) return
                hessian(:,coordinate_i,coordinate_j) = (mixed_pp-mixed_pm-mixed_mp+mixed_mm)/(4.0_wp*step*step)
                hessian(:,coordinate_j,coordinate_i) = hessian(:,coordinate_i,coordinate_j)
            enddo
        enddo
    end subroutine central_finite_difference

    !> @brief Report a numerical error when an evaluator returns a non-finite value.
    pure subroutine require_finite_values( values, status )
        real(wp), intent(in) :: values(:)
        type(chempotpy_status_t), intent(inout) :: status
        if( status%failed() ) return
        if( .not. all(ieee_is_finite(values)) ) &
            call status%set(status_numerical_error, 'finite-difference evaluator returned a non-finite value')
    end subroutine require_finite_values
end module chempotpy_finite_difference
