!> @brief Provide immutable evaluation kernels for contracted Surfgen polynomial models.
module chempotpy_surfgen
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_cubic_surfgen_hamiltonian
    public :: evaluate_surfgen_coordinates

    contains

    !> @brief Evaluate configured Surfgen internal coordinates without mutable initialization state.
    pure subroutine evaluate_surfgen_coordinates( cartesian, coordinate_types, scaling_types, atom_indices, &
            & parameters, coordinates, info )
        real(wp), intent(in)  :: cartesian(:,:)
        integer, intent(in)   :: coordinate_types(:), scaling_types(:), atom_indices(:,:)
        real(wp), intent(in)  :: parameters(:,:)
        real(wp), intent(out) :: coordinates(size(coordinate_types))
        integer, intent(out)  :: info
        integer               :: coordinate
        info        = 0
        coordinates = 0.0_wp
        if( size(cartesian, 2) /= 3 .or. size(scaling_types) /= size(coordinate_types) .or. &
            & size(atom_indices, 1) < 4 .or. size(atom_indices, 2) /= size(coordinate_types) .or. &
            & size(parameters, 1) < 2 .or. size(parameters, 2) /= size(coordinate_types) )then
        info = 1
        return
        endif
        do coordinate = 1, size(coordinate_types)
            select case( coordinate_types(coordinate) )
            case( 0 )
                coordinates(coordinate) = scaled_distance_between_atoms(cartesian, atom_indices(1,coordinate), &
                    & atom_indices(2,coordinate), scaling_types(coordinate), parameters(:,coordinate))
            case( -1 )
                coordinates(coordinate) = scaled_tetrahedral_coordinate(cartesian, atom_indices(:,coordinate), &
                    & scaling_types(coordinate), parameters(:,coordinate))
            case( -2 )
                coordinates(coordinate) = scaled_umbrella_coordinate(cartesian, atom_indices(:,coordinate), &
                    & scaling_types(coordinate), parameters(:,coordinate))
            case( -3 )
                coordinates(coordinate) = scaled_four_center_dot_product(cartesian, atom_indices(:,coordinate), &
                    & scaling_types(coordinate), parameters(:,coordinate))
            case( 1 )
                coordinates(coordinate) = scaled_bending_coordinate(cartesian, atom_indices(1:3,coordinate), &
                    & scaling_types(coordinate), parameters(:,coordinate))
            case default
                info = 2
                return
            end select
        enddo
    end subroutine evaluate_surfgen_coordinates

    !> @brief Contract normalized monomials through immutable cubic Surfgen coefficients.
    pure subroutine evaluate_cubic_surfgen_hamiltonian( coordinates, second_order_parent, second_order_direction, &
            & second_order_divisor, third_order_parent, third_order_direction, third_order_divisor, coefficient_order, &
            & coefficient_row, coefficient_column, coefficient_term, coefficient_value, hamiltonian )
        real(wp), intent(in)  :: coordinates(:), coefficient_value(:)
        integer, intent(in)   :: second_order_parent(:), second_order_direction(:), second_order_divisor(:)
        integer, intent(in)   :: third_order_parent(:), third_order_direction(:), third_order_divisor(:)
        integer, intent(in)   :: coefficient_order(:), coefficient_row(:), coefficient_column(:), coefficient_term(:)
        real(wp), intent(out) :: hamiltonian(:,:)
        real(wp)              :: second_order_terms(size(second_order_parent))
        real(wp)              :: third_order_terms(size(third_order_parent)), term_value
        integer               :: coefficient, row, column, term
        second_order_terms = coordinates(second_order_direction)*coordinates(second_order_parent) / &
            & real(second_order_divisor, wp)
        third_order_terms = coordinates(third_order_direction)*second_order_terms(third_order_parent) / &
            & real(third_order_divisor, wp)
        hamiltonian = 0.0_wp
        do coefficient = 1, size(coefficient_value)
            term = coefficient_term(coefficient)
            select case( coefficient_order(coefficient) )
            case( 0 )
                term_value = 1.0_wp
            case( 1 )
                term_value = coordinates(term)
            case( 2 )
                term_value = second_order_terms(term)
            case( 3 )
                term_value = third_order_terms(term)
            case default
                term_value = 0.0_wp
            end select
            row                     = coefficient_row(coefficient)
            column                  = coefficient_column(coefficient)
            hamiltonian(row,column) = hamiltonian(row,column)+coefficient_value(coefficient)*term_value
        enddo
        do row = 2, size(hamiltonian, 1)
            hamiltonian(row,1:row-1) = hamiltonian(1:row-1,row)
        enddo
    end subroutine evaluate_cubic_surfgen_hamiltonian

    !> @brief Return a configured two-atom distance coordinate.
    pure function scaled_distance_between_atoms( cartesian, first_atom, second_atom, scaling_type, parameters ) &
            & result(value)
        real(wp), intent(in) :: cartesian(:,:), parameters(:)
        integer, intent(in)  :: first_atom, second_atom, scaling_type
        real(wp)             :: value
        value = scaled_distance(norm2(cartesian(first_atom,:)-cartesian(second_atom,:)), scaling_type, parameters)
    end function scaled_distance_between_atoms

    !> @brief Apply one of the Surfgen radial scaling functions to a distance.
    pure function scaled_distance( distance, scaling_type, parameters ) result(value)
        real(wp), intent(in) :: distance, parameters(:)
        integer, intent(in)  :: scaling_type
        real(wp)             :: value
        select case( scaling_type )
        case( 0 )
            value = (distance-parameters(2))*parameters(1)
        case( 1 )
            value = exp(-parameters(1)*(distance-parameters(2)))
        case( 2 )
            value = exp(-parameters(1)*(distance-parameters(2))**2)
        case( 3 )
            value = exp(-parameters(1)*(distance-parameters(2)))/(distance+1.0e-40_wp)
        case( 4 )
            value = exp(parameters(1)*(parameters(2)-distance)) *(distance/parameters(2))**(parameters(1)* &
                & parameters(2))
        case( 5 )
            value = (parameters(2)/distance)**parameters(1)
        case( 6 )
            value = exp(-parameters(1)*(distance-parameters(2)))/(distance+parameters(2))
        case( 7 )
            value = tanh((distance-parameters(2))/parameters(1))/2.0_wp
        case default
            value = 0.0_wp
        end select
        if( abs(value) < 1.0e-30_wp ) value = sign(1.0e-30_wp, value)
    end function scaled_distance

    !> @brief Return the scaled tetrahedral out-of-plane coordinate used by Surfgen.
    pure function scaled_tetrahedral_coordinate( cartesian, atoms, scaling_type, parameters ) result(value)
        real(wp), intent(in) :: cartesian(:,:), parameters(:)
        integer, intent(in)  :: atoms(:), scaling_type
        real(wp)             :: value
        real(wp)             :: geometry(4,3), vectors(6,3), weights(6), determinant_value, denominator
        integer              :: first, second, pair
        geometry = cartesian(atoms(1:4),:)
        call tetrahedral_vectors(geometry, vectors)
        determinant_value = determinant_3x3(vectors(1:3,:))
        pair              = 0
        do first = 1, 3
            do second = first+1, 4
                pair = pair+1
                weights(pair) = scaled_distance(norm2(geometry(first,:)-geometry(second,:)), max(scaling_type, 0), &
                    & parameters)
            enddo
        enddo
        if( scaling_type < 0 )then
            denominator = (parameters(2)/sum(weights))**(parameters(1)+3.0_wp)
            value       = determinant_value*denominator
        else
            value = determinant_value*product(weights)
        endif
    end function scaled_tetrahedral_coordinate

    !> @brief Return the special Surfgen umbrella coordinate about a distinguished atom.
    pure function scaled_umbrella_coordinate( cartesian, atoms, scaling_type, parameters ) result(value)
        real(wp), intent(in) :: cartesian(:,:), parameters(:)
        integer, intent(in)  :: atoms(:), scaling_type
        real(wp)             :: value
        real(wp)             :: center_distances(3), geometry(4,3), scaled_distances(3), vectors(6,3)
        real(wp)             :: determinant_value, norm_squared(6), pair_products(3), unscaled_umbrella
        integer              :: atom
        geometry = cartesian(atoms(1:4),:)
        call tetrahedral_vectors(geometry, vectors)
        determinant_value = determinant_3x3(vectors(1:3,:))
        do atom = 1, 3
            center_distances(atom) = norm2(geometry(1,:)-geometry(atom+1,:))
            scaled_distances(atom) = scaled_distance(center_distances(atom), abs(scaling_type), parameters)
        enddo
        unscaled_umbrella = determinant_value/product(center_distances)
        if( scaling_type < 0 )then
            pair_products = [scaled_distances(2)*scaled_distances(3), scaled_distances(3)*scaled_distances(1), &
                & scaled_distances(1)*scaled_distances(2)]
            value = unscaled_umbrella*product(scaled_distances)/sum(pair_products)
        elseif( scaling_type == 0 )then
            value = unscaled_umbrella*(parameters(2)**3/product(center_distances))**parameters(1)
        else
            do atom = 1, 6
                norm_squared(atom) = dot_product(vectors(atom,:), vectors(atom,:))
            enddo
            unscaled_umbrella = determinant_value/product(norm_squared)**0.25_wp *2.0_wp**(-1.5_wp)*3.0_wp**1.5_wp
            pair_products = [scaled_distances(2)*scaled_distances(3), scaled_distances(3)*scaled_distances(1), &
                & scaled_distances(1)*scaled_distances(2)]
            value = unscaled_umbrella*product(scaled_distances)/sum(pair_products)
        endif
    end function scaled_umbrella_coordinate

    !> @brief Return a normalized four-center dot product with optional radial scaling.
    pure function scaled_four_center_dot_product( cartesian, atoms, scaling_type, parameters ) result(value)
        real(wp), intent(in) :: cartesian(:,:), parameters(:)
        integer, intent(in)  :: atoms(:), scaling_type
        real(wp)             :: value
        real(wp)             :: first_vector(3), second_vector(3), weights(4), scale
        integer              :: first, second, weight
        first_vector  = cartesian(atoms(1),:)-cartesian(atoms(2),:)
        second_vector = cartesian(atoms(3),:)-cartesian(atoms(4),:)
        value         = dot_product(first_vector, second_vector)/(norm2(first_vector)*norm2(second_vector))
        if( scaling_type == 0 ) return
        weight = 0
        do first = 1, 2
            do second = 3, 4
                weight = weight+1
                weights(weight) = scaled_distance(norm2(cartesian(atoms(first),:)-cartesian(atoms(second),:)), &
                    & abs(scaling_type), parameters)
            enddo
        enddo
        if( scaling_type > 0 )then
            scale = product(weights)
        else
            scale = 4.0_wp/sum(1.0_wp/weights)
        endif
        value = value*scale
    end function scaled_four_center_dot_product

    !> @brief Return a Surfgen bond-angle coordinate with its configured scaling.
    pure function scaled_bending_coordinate( cartesian, atoms, scaling_type, parameters ) result(value)
        real(wp), intent(in) :: cartesian(:,:), parameters(:)
        integer, intent(in)  :: atoms(:), scaling_type
        real(wp)             :: value
        real(wp)             :: first_vector(3), second_vector(3), cosine_value, exponent_value
        first_vector  = cartesian(atoms(1),:)-cartesian(atoms(3),:)
        second_vector = cartesian(atoms(2),:)-cartesian(atoms(3),:)
        cosine_value  = dot_product(first_vector, second_vector)/(norm2(first_vector)*norm2(second_vector))
        select case( scaling_type )
        case( 0 )
            value = (acos(cosine_value)-parameters(2))*parameters(1)
        case( 1 )
            value = cosine_value
        case( 2 )
            exponent_value = parameters(1)*(dot_product(first_vector, first_vector) +dot_product(second_vector, &
                & second_vector)-parameters(2)**2)
            if( exponent_value > 70.0_wp )then
                value = 0.0_wp
            else
                value = cosine_value/(1.0_wp+exp(exponent_value))
            endif
        case default
            value = 0.0_wp
        end select
    end function scaled_bending_coordinate

    !> @brief Form the six oriented displacement vectors of four atoms.
    pure subroutine tetrahedral_vectors( geometry, vectors )
        real(wp), intent(in)  :: geometry(4,3)
        real(wp), intent(out) :: vectors(6,3)
        vectors(1,:) = geometry(1,:)-geometry(4,:)
        vectors(2,:) = geometry(2,:)-geometry(4,:)
        vectors(3,:) = geometry(3,:)-geometry(4,:)
        vectors(4,:) = geometry(1,:)-geometry(3,:)
        vectors(5,:) = geometry(2,:)-geometry(3,:)
        vectors(6,:) = geometry(1,:)-geometry(2,:)
    end subroutine tetrahedral_vectors

    !> @brief Return the determinant of a three-by-three matrix.
    pure function determinant_3x3( matrix ) result(value)
        real(wp), intent(in) :: matrix(3,3)
        real(wp)             :: value
        value = (matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))*matrix(3,1) +(matrix(1,3)*matrix(2,1)-matrix(1,1)* &
            & matrix(2,3))*matrix(3,2) +(matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))*matrix(3,3)
    end function determinant_3x3
end module chempotpy_surfgen
