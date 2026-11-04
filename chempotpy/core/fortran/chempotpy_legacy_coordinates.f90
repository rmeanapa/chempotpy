!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Provide explicit-argument adapters for released POTLIB coordinates.
module chempotpy_legacy_coordinates
    use chempotpy_geometry,     only: cartesian_distance
    use chempotpy_kinds,        only: real64
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom
    implicit none
    private
    public :: convert_cartesian_input
    public :: cartesian_to_potlib_coordinates
    public :: legacy_triatomic_distances_to_bohr
    contains
    !> @brief Return the 1-2, 2-3, and 1-3 distances in legacy bohr.
    pure subroutine legacy_triatomic_distances_to_bohr( &
            coordinates_angstrom, distances_bohr )
        real(real64), intent(in)  :: coordinates_angstrom(:,:)
        real(real64), intent(out) :: distances_bohr(3)
        distances_bohr(1) = cartesian_distance( coordinates_angstrom(1,:), coordinates_angstrom(2,:) )/ &
            & legacy_bohr_to_angstrom
        distances_bohr(2) = cartesian_distance( coordinates_angstrom(2,:), coordinates_angstrom(3,:) )/ &
            & legacy_bohr_to_angstrom
        distances_bohr(3) = cartesian_distance( coordinates_angstrom(1,:), coordinates_angstrom(3,:) )/ &
            & legacy_bohr_to_angstrom
    end subroutine legacy_triatomic_distances_to_bohr
    !> @brief Convert user Cartesian coordinates to POTLIB distance units.
    pure subroutine convert_cartesian_input( cartesian_coordinates, &
            atom_labels, atom_count, distance_conversion, reorder_atoms, &
            internal_coordinates )
        real(real64), intent(in)  :: cartesian_coordinates(:,:)
        integer, intent(in)       :: atom_labels(:)
        integer, intent(in)       :: atom_count
        real(real64), intent(in)  :: distance_conversion
        logical, intent(in)       :: reorder_atoms
        real(real64), intent(out) :: internal_coordinates(:,:)
        integer                   :: atom
        if( reorder_atoms )then
            do atom = 1, atom_count
                internal_coordinates(atom_labels(atom),1:3) = cartesian_coordinates(atom,1:3)/distance_conversion
            enddo
        else
            internal_coordinates(1:atom_count,1:3) = cartesian_coordinates(1:atom_count,1:3)/distance_conversion
        endif
    end subroutine convert_cartesian_input
    !> @brief Build a released POTLIB coordinate vector from Cartesian input.
    pure subroutine cartesian_to_potlib_coordinates( &
            cartesian_coordinates, atom_count, coordinate_mode, coordinates, &
            status )
        real(real64), intent(in)    :: cartesian_coordinates(:,:)
        integer, intent(in)         :: atom_count
        integer, intent(in)         :: coordinate_mode
        real(real64), intent(inout) :: coordinates(:)
        integer, intent(out)        :: status
        real(real64), parameter     :: fluorine_mass = 18.99840_real64
        real(real64), parameter     :: hydrogen_mass = 1.007825_real64
        real(real64)                :: center_1(3), center_2(3)
        real(real64)                :: center_delta(3)
        real(real64)                :: relative_1(3), relative_2(3)
        real(real64)                :: center_distance, hydrogen_distance
        real(real64)                :: projected_distance
        real(real64)                :: q1, q2, q3, cosine_phi, phi
        real(real64)                :: theta_1, theta_2
        integer                     :: atom_1, atom_2, coordinate
        status = 0
        select case (coordinate_mode)
        case (0)
            return
        case (1)
            do atom_1 = 1, atom_count
                coordinate                           = 3*atom_1-2
                coordinates(coordinate:coordinate+2) = cartesian_coordinates(atom_1,1:3)
            enddo
        case (2)
            coordinate = 1
            do atom_1 = 1, atom_count-1
                do atom_2 = atom_1+1, atom_count
                    coordinates(coordinate) = sqrt( (cartesian_coordinates(atom_1,1)-cartesian_coordinates(atom_2, &
                        & 1))** 2 +(cartesian_coordinates(atom_1,2)-cartesian_coordinates(atom_2,2))**2 + &
                        & (cartesian_coordinates(atom_1,3)-cartesian_coordinates(atom_2,3))**2 )
                    coordinate = coordinate+1
                enddo
            enddo
        case (3)
            coordinates(1) = sqrt( (cartesian_coordinates(1,1)-cartesian_coordinates(2,1))**2 + &
                & (cartesian_coordinates(1,2)-cartesian_coordinates(2,2))**2 +(cartesian_coordinates(1,3)- &
                & cartesian_coordinates(2,3))**2 )
            coordinates(2) = sqrt( (cartesian_coordinates(2,1)-cartesian_coordinates(3,1))**2 + &
                & (cartesian_coordinates(2,2)-cartesian_coordinates(3,2))**2 +(cartesian_coordinates(2,3)- &
                & cartesian_coordinates(3,3))**2 )
            coordinates(3) = sqrt( (cartesian_coordinates(1,1)-cartesian_coordinates(3,1))**2 + &
                & (cartesian_coordinates(1,2)-cartesian_coordinates(3,2))**2 +(cartesian_coordinates(1,3)- &
                & cartesian_coordinates(3,3))**2 )
        case (4)
            center_1 = (hydrogen_mass*cartesian_coordinates(1,1:3) +fluorine_mass*cartesian_coordinates(2,1:3)) / &
                & (fluorine_mass+hydrogen_mass)
            center_2 = (hydrogen_mass*cartesian_coordinates(3,1:3) +fluorine_mass*cartesian_coordinates(4,1:3)) / &
                & (fluorine_mass+hydrogen_mass)
            center_delta       = center_2-center_1
            relative_1         = cartesian_coordinates(1,1:3)-center_1
            relative_2         = cartesian_coordinates(3,1:3)-center_2
            theta_1            = sum(relative_1*center_delta)
            theta_1            = theta_1/sqrt(sum(relative_1**2))
            theta_1            = theta_1/sqrt(sum(center_delta**2))
            theta_1            = min(1.0_real64, max(-1.0_real64, theta_1))
            theta_1            = acos(theta_1)
            theta_2            = sum(relative_2*(-center_delta))
            theta_2            = theta_2/sqrt(sum(relative_2**2))
            theta_2            = theta_2/sqrt(sum(center_delta**2))
            theta_2            = min(1.0_real64, max(-1.0_real64, theta_2))
            theta_2            = acos(theta_2)
            theta_2            = acos(-1.0_real64)-theta_2
            q1                 = sqrt(sum(relative_1**2)); q2 = sqrt(sum(relative_2**2))
            center_distance    = sqrt(sum(center_delta**2))
            hydrogen_distance  = sqrt(sum( (cartesian_coordinates(1,1:3)-cartesian_coordinates(3,1:3))**2))
            projected_distance = center_distance-q1*cos(theta_1)+q2*cos(theta_2)
            q3                 = sqrt(abs(hydrogen_distance**2-projected_distance**2))
            q1                 = q1*sin(theta_1); q2 = q2*sin(theta_2)
            cosine_phi         = (q1**2+q2**2-q3**2)/(2.0_real64*q1*q2)
            cosine_phi         = min(1.0_real64, max(-1.0_real64, cosine_phi))
            phi                = acos(cosine_phi)
            coordinates(1)     = center_distance
            coordinates(2)     = sqrt(sum(relative_1**2)) *(fluorine_mass+hydrogen_mass)/fluorine_mass
            coordinates(3)     = sqrt(sum(relative_2**2)) *(fluorine_mass+hydrogen_mass)/fluorine_mass
            coordinates(4)     = theta_1
            coordinates(5)     = theta_2
            coordinates(6)     = phi
        case default
            status = 1
        end select
    end subroutine cartesian_to_potlib_coordinates
end module chempotpy_legacy_coordinates
