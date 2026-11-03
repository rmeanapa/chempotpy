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

contains

    !> @brief Return the Euclidean distance between two Cartesian points.
    pure function cartesian_distance( first_point, second_point ) &
        result(distance)
        real(real64), intent(in) :: first_point(3), second_point(3)
        real(real64) :: distance
        real(real64) :: delta_x, delta_y, delta_z

        delta_x = first_point(1)-second_point(1)
        delta_y = first_point(2)-second_point(2)
        delta_z = first_point(3)-second_point(3)
        distance = sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z)
    end function cartesian_distance


    !> @brief Return the Cartesian dot product of two three-component vectors.
    pure function dot_product_3d( first_vector, second_vector ) &
        result(product)
        real(real64), intent(in) :: first_vector(3), second_vector(3)
        real(real64) :: product

        product = first_vector(1)*second_vector(1) &
            +first_vector(2)*second_vector(2) &
            +first_vector(3)*second_vector(3)
    end function dot_product_3d


    !> @brief Return the cosine of the angle between two Cartesian vectors.
    pure function cosine_between_vectors( &
        first_vector, second_vector ) result(cosine)
        real(real64), intent(in) :: first_vector(3), second_vector(3)
        real(real64) :: cosine
        real(real64) :: first_norm, second_norm

        first_norm = sqrt(dot_product_3d(first_vector,first_vector))
        second_norm = sqrt(dot_product_3d(second_vector,second_vector))
        cosine = dot_product_3d(first_vector,second_vector) &
            /(first_norm*second_norm)
    end function cosine_between_vectors


    !> @brief Fill a symmetric Cartesian distance matrix for a geometry.
    pure subroutine fill_distance_matrix( coordinates, distances )
        real(real64), intent(in) :: coordinates(:,:)
        real(real64), intent(out) :: distances(size(coordinates,1), &
            size(coordinates,1))
        integer :: first_atom, second_atom

        distances = 0.0_real64
        do first_atom = 1, size(coordinates,1)-1
            do second_atom = first_atom+1, size(coordinates,1)
                distances(first_atom,second_atom) = cartesian_distance( &
                    coordinates(first_atom,:), coordinates(second_atom,:) )
                distances(second_atom,first_atom) = &
                    distances(first_atom,second_atom)
            enddo
        enddo
    end subroutine fill_distance_matrix
end module chempotpy_geometry
