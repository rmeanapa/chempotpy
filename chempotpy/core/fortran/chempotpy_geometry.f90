!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Provide small Cartesian geometry operations used by PES kernels.
module chempotpy_geometry
    use chempotpy_kinds, only: real64
    implicit none
    private
    public :: cartesian_distance
    public :: cosine_between_vectors
    public :: dot_product_3d
    public :: fill_distance_matrix
    public :: flatten_cartesian_coordinates
    public :: triatomic_distance_jacobian
    public :: transform_internal_gradient_to_cartesian
    interface transform_internal_gradient_to_cartesian
        module procedure transform_gradient_vector
        module procedure transform_gradient_matrix
        module procedure transform_gradient_tensor
    end interface transform_internal_gradient_to_cartesian
    contains
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
    pure subroutine fill_distance_matrix( coordinates, distances )
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
    end subroutine fill_distance_matrix
end module chempotpy_geometry
