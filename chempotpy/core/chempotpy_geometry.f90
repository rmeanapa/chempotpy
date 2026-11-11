!> @brief Provide small Cartesian geometry operations used by PES kernels.
module chempotpy_geometry
    use chempotpy_kinds, only: real64
    implicit none
    private
    public :: all_pair_distance_jacobian
    public :: all_pair_distances
    public :: all_pair_distances_and_jacobian
    public :: cartesian_distance
    public :: cosine_between_vectors
    public :: dot_product_3d
    public :: fill_distance_matrix
    public :: flatten_cartesian_coordinates
    public :: ordered_pair_distance_jacobian
    public :: triatomic_distance_jacobian
    public :: transform_internal_gradient_to_cartesian
    interface transform_internal_gradient_to_cartesian
        module procedure transform_gradient_vector
        module procedure transform_gradient_matrix
        module procedure transform_gradient_tensor
    end interface transform_internal_gradient_to_cartesian
    interface all_pair_distances_and_jacobian
        module procedure all_pair_distances_from_matrix
        module procedure all_pair_distances_from_vector
    end interface all_pair_distances_and_jacobian
    interface all_pair_distances
        module procedure pair_distances_from_matrix
        module procedure pair_distances_from_vector
    end interface all_pair_distances
    interface all_pair_distance_jacobian
        module procedure all_pair_jacobian_from_matrix
        module procedure all_pair_jacobian_from_vector
    end interface all_pair_distance_jacobian
    interface ordered_pair_distance_jacobian
        module procedure ordered_pair_jacobian_from_matrix
        module procedure ordered_pair_jacobian_from_vector
    end interface ordered_pair_distance_jacobian
    interface fill_distance_matrix
        module procedure fill_distance_matrix_from_components
        module procedure fill_distance_matrix_from_matrix
    end interface fill_distance_matrix

    contains

    !> @brief Return a Cartesian Jacobian for an explicitly ordered list of atom pairs.
    pure subroutine ordered_pair_jacobian_from_matrix( coordinates, first_atoms, second_atoms, distances, jacobian )
        real(real64), intent(in)  :: coordinates(:,:), distances(:)
        integer, intent(in)       :: first_atoms(:), second_atoms(:)
        real(real64), intent(out) :: jacobian(size(distances),3*size(coordinates, 1))
        real(real64)              :: displacement(3)
        integer                   :: first_coordinate, pair, second_coordinate
        jacobian = 0.0_real64
        do pair = 1, size(distances)
            first_coordinate  = 3*first_atoms(pair)-2
            second_coordinate = 3*second_atoms(pair)-2
            displacement                                          = coordinates(first_atoms(pair),1:3) - &
                & coordinates(second_atoms(pair),1:3)
            jacobian(pair,first_coordinate:first_coordinate+2)   = displacement/distances(pair)
            jacobian(pair,second_coordinate:second_coordinate+2) = -displacement/distances(pair)
        enddo
    end subroutine ordered_pair_jacobian_from_matrix

    !> @brief Return a Cartesian Jacobian for ordered pairs and atom-major vector coordinates.
    pure subroutine ordered_pair_jacobian_from_vector( coordinates, first_atoms, second_atoms, distances, jacobian )
        real(real64), intent(in)  :: coordinates(:), distances(:)
        integer, intent(in)       :: first_atoms(:), second_atoms(:)
        real(real64), intent(out) :: jacobian(size(distances),size(coordinates))
        real(real64)              :: displacement(3)
        integer                   :: first_coordinate, pair, second_coordinate
        jacobian = 0.0_real64
        do pair = 1, size(distances)
            first_coordinate  = 3*first_atoms(pair)-2
            second_coordinate = 3*second_atoms(pair)-2
            displacement                                          = coordinates(first_coordinate:first_coordinate+2) - &
                & coordinates(second_coordinate:second_coordinate+2)
            jacobian(pair,first_coordinate:first_coordinate+2)   = displacement/distances(pair)
            jacobian(pair,second_coordinate:second_coordinate+2) = -displacement/distances(pair)
        enddo
    end subroutine ordered_pair_jacobian_from_vector

    !> @brief Return all unique pair distances from matrix coordinates.
    pure subroutine pair_distances_from_matrix( coordinates, distances )
        real(real64), intent(in)  :: coordinates(:,:)
        real(real64), intent(out) :: distances(size(coordinates, 1)*(size(coordinates, 1)-1)/2)
        real(real64)              :: displacement(3)
        integer                   :: first_atom, pair, second_atom
        pair = 0
        do first_atom = 1, size(coordinates, 1)-1
            do second_atom = first_atom+1, size(coordinates, 1)
                pair            = pair+1
                displacement    = coordinates(first_atom,1:3)-coordinates(second_atom,1:3)
                distances(pair) = sqrt(sum(displacement*displacement))
            enddo
        enddo
    end subroutine pair_distances_from_matrix

    !> @brief Return all unique pair distances from atom-major coordinates.
    pure subroutine pair_distances_from_vector( coordinates, distances )
        real(real64), intent(in)  :: coordinates(:)
        real(real64), intent(out) :: distances((size(coordinates)/3)*(size(coordinates)/3-1)/2)
        real(real64)              :: displacement(3)
        integer                   :: atom_count, first_atom, first_coordinate, pair, second_atom, second_coordinate
        atom_count = size(coordinates)/3
        pair       = 0
        do first_atom = 1, atom_count-1
            first_coordinate = 3*first_atom-2
            do second_atom = first_atom+1, atom_count
                second_coordinate = 3*second_atom-2
                pair              = pair+1
                displacement      = coordinates(first_coordinate:first_coordinate+2) - &
                    & coordinates(second_coordinate:second_coordinate+2)
                distances(pair) = sqrt(sum(displacement*displacement))
            enddo
        enddo
    end subroutine pair_distances_from_vector

    !> @brief Return all unique pair distances and their Cartesian Jacobian from matrix coordinates.
    pure subroutine all_pair_distances_from_matrix( coordinates, distances, jacobian )
        real(real64), intent(in)  :: coordinates(:,:)
        real(real64), intent(out) :: distances(size(coordinates, 1)*(size(coordinates, 1)-1)/2)
        real(real64), intent(out) :: jacobian(size(distances),3*size(coordinates, 1))
        real(real64)              :: displacement(3)
        integer                   :: first_atom, first_coordinate, pair, second_atom, second_coordinate
        jacobian = 0.0_real64
        pair     = 0
        do first_atom = 1, size(coordinates, 1)-1
            first_coordinate = 3*first_atom-2
            do second_atom = first_atom+1, size(coordinates, 1)
                second_coordinate = 3*second_atom-2
                pair              = pair+1
                displacement                                         = coordinates(first_atom,1:3)- &
                    & coordinates(second_atom,1:3)
                distances(pair)                                      = sqrt(sum(displacement*displacement))
                jacobian(pair,first_coordinate:first_coordinate+2)   = displacement/distances(pair)
                jacobian(pair,second_coordinate:second_coordinate+2) = -displacement/distances(pair)
            enddo
        enddo
    end subroutine all_pair_distances_from_matrix

    !> @brief Return all unique pair distances and their Cartesian Jacobian from atom-major coordinates.
    pure subroutine all_pair_distances_from_vector( coordinates, distances, jacobian )
        real(real64), intent(in)  :: coordinates(:)
        real(real64), intent(out) :: distances((size(coordinates)/3)*(size(coordinates)/3-1)/2)
        real(real64), intent(out) :: jacobian(size(distances),size(coordinates))
        real(real64)              :: displacement(3)
        integer                   :: atom_count, first_atom, first_coordinate, pair, second_atom, second_coordinate
        atom_count = size(coordinates)/3
        jacobian   = 0.0_real64
        pair       = 0
        do first_atom = 1, atom_count-1
            first_coordinate = 3*first_atom-2
            do second_atom = first_atom+1, atom_count
                second_coordinate = 3*second_atom-2
                pair              = pair+1
                displacement      = coordinates(first_coordinate:first_coordinate+2) - &
                    & coordinates(second_coordinate:second_coordinate+2)
                distances(pair)                                      = sqrt(sum(displacement*displacement))
                jacobian(pair,first_coordinate:first_coordinate+2)   = displacement/distances(pair)
                jacobian(pair,second_coordinate:second_coordinate+2) = -displacement/distances(pair)
            enddo
        enddo
    end subroutine all_pair_distances_from_vector

    !> @brief Return the Cartesian Jacobian for supplied all-pair distances and matrix coordinates.
    pure subroutine all_pair_jacobian_from_matrix( coordinates, distances, jacobian )
        real(real64), intent(in)  :: coordinates(:,:), distances(:)
        real(real64), intent(out) :: jacobian(size(distances),3*size(coordinates, 1))
        real(real64)              :: displacement(3)
        integer                   :: first_atom, first_coordinate, pair, second_atom, second_coordinate
        jacobian = 0.0_real64
        pair     = 0
        do first_atom = 1, size(coordinates, 1)-1
            first_coordinate = 3*first_atom-2
            do second_atom = first_atom+1, size(coordinates, 1)
                second_coordinate = 3*second_atom-2
                pair              = pair+1
                displacement                                         = coordinates(first_atom,1:3)- &
                    & coordinates(second_atom,1:3)
                jacobian(pair,first_coordinate:first_coordinate+2)   = displacement/distances(pair)
                jacobian(pair,second_coordinate:second_coordinate+2) = -displacement/distances(pair)
            enddo
        enddo
    end subroutine all_pair_jacobian_from_matrix

    !> @brief Return the Cartesian Jacobian for supplied all-pair distances and atom-major coordinates.
    pure subroutine all_pair_jacobian_from_vector( coordinates, distances, jacobian )
        real(real64), intent(in)  :: coordinates(:), distances(:)
        real(real64), intent(out) :: jacobian(size(distances),size(coordinates))
        real(real64)              :: displacement(3)
        integer                   :: atom_count, first_atom, first_coordinate, pair, second_atom, second_coordinate
        atom_count = size(coordinates)/3
        jacobian   = 0.0_real64
        pair       = 0
        do first_atom = 1, atom_count-1
            first_coordinate = 3*first_atom-2
            do second_atom = first_atom+1, atom_count
                second_coordinate = 3*second_atom-2
                pair              = pair+1
                displacement      = coordinates(first_coordinate:first_coordinate+2) - &
                    & coordinates(second_coordinate:second_coordinate+2)
                jacobian(pair,first_coordinate:first_coordinate+2)   = displacement/distances(pair)
                jacobian(pair,second_coordinate:second_coordinate+2) = -displacement/distances(pair)
            enddo
        enddo
    end subroutine all_pair_jacobian_from_vector

    !> @brief Flatten atom-major Cartesian coordinates as x1, y1, z1, x2, y2, z2, and so on.
    pure subroutine flatten_cartesian_coordinates( coordinates, flat_coordinates )
        real(real64), intent(in)  :: coordinates(:,:)
        real(real64), intent(out) :: flat_coordinates(3*size(coordinates, 1))
        integer                   :: atom, first_coordinate
        do atom = 1, size(coordinates, 1)
            first_coordinate                                      = 3*atom-2
            flat_coordinates(first_coordinate:first_coordinate+2) = coordinates(atom,1:3)
        enddo
    end subroutine flatten_cartesian_coordinates

    !> @brief Return the 1-2, 2-3, and 1-3 distances and their Cartesian Jacobian.
    pure subroutine triatomic_distance_jacobian( coordinates, distances, jacobian )
        real(real64), intent(in)  :: coordinates(3,3)
        real(real64), intent(out) :: distances(3), jacobian(3,9)
        real(real64)              :: displacement(3)
        integer                   :: coordinate, pair
        integer, parameter        :: first_atom(3)  = [1, 2, 1]
        integer, parameter        :: second_atom(3) = [2, 3, 3]
        jacobian = 0.0_real64
        do pair = 1, 3
            displacement    = coordinates(first_atom(pair),:)-coordinates(second_atom(pair),:)
            distances(pair) = sqrt(sum(displacement*displacement))
            do coordinate = 1, 3
                jacobian(pair,3*first_atom(pair)-3+coordinate)  = displacement(coordinate)/distances(pair)
                jacobian(pair,3*second_atom(pair)-3+coordinate) = -displacement(coordinate)/distances(pair)
            enddo
        enddo
    end subroutine triatomic_distance_jacobian

    !> @brief Apply dq/dx to one internal-coordinate gradient using the chain rule.
    pure subroutine transform_gradient_vector( internal_gradient, coordinate_jacobian, cartesian_gradient )
        real(real64), intent(in)  :: internal_gradient(:), coordinate_jacobian(:,:)
        real(real64), intent(out) :: cartesian_gradient(size(coordinate_jacobian, 2))
        integer                   :: cartesian_index, internal_index
        cartesian_gradient = 0.0_real64
        do cartesian_index = 1, size(coordinate_jacobian, 2)
            do internal_index = 1, size(coordinate_jacobian, 1)
                cartesian_gradient(cartesian_index) = cartesian_gradient(cartesian_index) + &
                    & internal_gradient(internal_index)*coordinate_jacobian(internal_index,cartesian_index)
            enddo
        enddo
    end subroutine transform_gradient_vector

    !> @brief Apply dq/dx to multiple internal-coordinate gradients using the chain rule.
    pure subroutine transform_gradient_matrix( internal_gradient, coordinate_jacobian, cartesian_gradient )
        real(real64), intent(in)  :: internal_gradient(:,:), coordinate_jacobian(:,:)
        real(real64), intent(out) :: cartesian_gradient(size(internal_gradient, 1),size(coordinate_jacobian, 2))
        integer                   :: cartesian_index, internal_index
        cartesian_gradient = 0.0_real64
        do cartesian_index = 1, size(coordinate_jacobian, 2)
            do internal_index = 1, size(coordinate_jacobian, 1)
                cartesian_gradient(:,cartesian_index) = cartesian_gradient(:,cartesian_index) + internal_gradient(:, &
                    & internal_index)*coordinate_jacobian(internal_index,cartesian_index)
            enddo
        enddo
    end subroutine transform_gradient_matrix

    !> @brief Apply dq/dx to a matrix-valued internal gradient using the chain rule.
    pure subroutine transform_gradient_tensor( internal_gradient, coordinate_jacobian, cartesian_gradient )
        real(real64), intent(in)  :: internal_gradient(:,:,:), coordinate_jacobian(:,:)
        real(real64), intent(out) :: cartesian_gradient(size(internal_gradient, 1),size(internal_gradient, 2),&
            & size(coordinate_jacobian, 2))
        integer :: cartesian_index, internal_index
        cartesian_gradient = 0.0_real64
        do cartesian_index = 1, size(coordinate_jacobian, 2)
            do internal_index = 1, size(coordinate_jacobian, 1)
                cartesian_gradient(:,:,cartesian_index) = cartesian_gradient(:,:,cartesian_index) + &
                    & internal_gradient(:,:,internal_index)*coordinate_jacobian(internal_index,cartesian_index)
            enddo
        enddo
    end subroutine transform_gradient_tensor

    !> @brief Return the Euclidean distance between two Cartesian points.
    pure function cartesian_distance( first_point, second_point ) &
            result(distance)
        real(real64), intent(in) :: first_point(3), second_point(3)
        real(real64)             :: distance
        real(real64)             :: delta_x, delta_y, delta_z
        delta_x  = first_point(1)-second_point(1)
        delta_y  = first_point(2)-second_point(2)
        delta_z  = first_point(3)-second_point(3)
        distance = sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z)
    end function cartesian_distance

    !> @brief Return the Cartesian dot product of two three-component vectors.
    pure function dot_product_3d( first_vector, second_vector ) &
            result(product)
        real(real64), intent(in) :: first_vector(3), second_vector(3)
        real(real64)             :: product
        product = first_vector(1)*second_vector(1) +first_vector(2)*second_vector(2) +first_vector(3)*second_vector(3)
    end function dot_product_3d

    !> @brief Return the cosine of the angle between two Cartesian vectors.
    pure function cosine_between_vectors( &
            first_vector, second_vector ) result(cosine)
        real(real64), intent(in) :: first_vector(3), second_vector(3)
        real(real64)             :: cosine
        real(real64)             :: first_norm, second_norm
        first_norm  = sqrt(dot_product_3d(first_vector, first_vector))
        second_norm = sqrt(dot_product_3d(second_vector, second_vector))
        cosine      = dot_product_3d(first_vector, second_vector) /(first_norm*second_norm)
    end function cosine_between_vectors

    !> @brief Fill a symmetric Cartesian distance matrix for a geometry.
    pure subroutine fill_distance_matrix_from_matrix( coordinates, distances )
        real(real64), intent(in)  :: coordinates(:,:)
        real(real64), intent(out) :: distances(size(coordinates, 1),&
            size(coordinates, 1))
        integer :: first_atom, second_atom
        distances = 0.0_real64
        do first_atom = 1, size(coordinates, 1)-1
            do second_atom = first_atom+1, size(coordinates, 1)
                distances(first_atom,second_atom) = cartesian_distance( coordinates(first_atom,:), &
                    & coordinates(second_atom,:) )
                distances(second_atom,first_atom) = distances(first_atom,second_atom)
            enddo
        enddo
    end subroutine fill_distance_matrix_from_matrix

    !> @brief Fill a symmetric distance matrix from separate Cartesian component arrays.
    pure subroutine fill_distance_matrix_from_components( x_coordinates, y_coordinates, z_coordinates, distances )
        real(real64), intent(in)  :: x_coordinates(:), y_coordinates(:), z_coordinates(:)
        real(real64), intent(out) :: distances(size(x_coordinates),size(x_coordinates))
        real(real64)              :: delta_x, delta_y, delta_z
        integer                   :: first_atom, second_atom
        distances = 0.0_real64
        do first_atom = 1, size(x_coordinates)-1
            do second_atom = first_atom+1, size(x_coordinates)
                delta_x                           = x_coordinates(first_atom)-x_coordinates(second_atom)
                delta_y                           = y_coordinates(first_atom)-y_coordinates(second_atom)
                delta_z                           = z_coordinates(first_atom)-z_coordinates(second_atom)
                distances(first_atom,second_atom) = sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z)
                distances(second_atom,first_atom) = distances(first_atom,second_atom)
            enddo
        enddo
    end subroutine fill_distance_matrix_from_components
end module chempotpy_geometry
