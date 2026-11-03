!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the Schwenke-Partridge methane potential.
subroutine pes( x, igrad, path, p, g, d )
    implicit none

    integer, parameter :: nstates = 1
    integer, parameter :: natoms = 5
    real(8), parameter :: bohr_to_angstrom = 0.529177211d0
    real(8), parameter :: hartree_to_ev = 27.211386d0
    real(8), intent(in) :: x(natoms,3)
    integer, intent(in) :: igrad
    character(len=1024), intent(in) :: path
    real(8), intent(out) :: p(nstates)
    real(8), intent(out) :: g(nstates,natoms,3)
    real(8), intent(out) :: d(nstates,nstates,natoms,3)
    real(8) :: coordinates_bohr(natoms,3)
    real(8) :: energy_hartree

    p = 0.0d0
    g = 0.0d0
    d = 0.0d0
    coordinates_bohr = x/bohr_to_angstrom
    call getpot( energy_hartree, coordinates_bohr, path )
    if( igrad == 0 ) then
        p(1) = energy_hartree*hartree_to_ev
    else
        write(*,*) 'Only energy is available'
    endif
end subroutine pes


!> @brief Convert methane Cartesian coordinates and call the T8 kernel.
subroutine getpot( potential_energy, coordinates, path )
    implicit none

    character(len=1024), intent(in) :: path
    real(8), intent(in) :: coordinates(5,3)
    real(8), intent(out) :: potential_energy
    real(8), parameter :: zero_point_shift = &
        -28.11193019d0*4.556335d-6
    real(8) :: energy_hartree
    real(8) :: radau_coordinates(10)

    call methaneradau( coordinates, radau_coordinates )
    call vibpot( radau_coordinates, energy_hartree, 1, path )
    potential_energy = energy_hartree-zero_point_shift
end subroutine getpot


!> @brief Convert CH4 Cartesian coordinates to the T8 Radau coordinates.
subroutine methaneradau( coordinates, radau_coordinates )
    implicit none

    real(8), intent(in) :: coordinates(5,3)
    real(8), intent(out) :: radau_coordinates(10)
    real(8), parameter :: unitconvm = 5.48579903d-4
    real(8), save :: alpha0 = 0.0d0
    logical, save :: initialized = .false.
    real(8) :: carbon_mass, hydrogen_mass
    real(8) :: relative(4,3), shifted(4,3), shift(3)
    integer :: atom_index

    if( .not. initialized ) then
        carbon_mass = 12.0d0/unitconvm-6.0d0
        hydrogen_mass = 1.007825d0/unitconvm-1.0d0
        alpha0 = sqrt(carbon_mass/(carbon_mass+4.0d0*hydrogen_mass))
        initialized = .true.
    endif

    do atom_index = 1, 4
        relative(atom_index,:) = &
            coordinates(atom_index+1,:)-coordinates(1,:)
    enddo
    shift = 0.25d0*(1.0d0-alpha0)*sum(relative, dim=1)
    do atom_index = 1, 4
        shifted(atom_index,:) = relative(atom_index,:)-shift
        radau_coordinates(atom_index) = &
            sqrt(sum(shifted(atom_index,:)**2))
    enddo

    call radau_angle( relative(1,:), relative(2,:), &
        radau_coordinates(1), radau_coordinates(2), &
        radau_coordinates(5) )
    call radau_angle( relative(1,:), relative(3,:), &
        radau_coordinates(1), radau_coordinates(3), &
        radau_coordinates(6) )
    call radau_angle( relative(1,:), relative(4,:), &
        radau_coordinates(1), radau_coordinates(4), &
        radau_coordinates(7) )
    call radau_angle( relative(2,:), relative(3,:), &
        radau_coordinates(2), radau_coordinates(3), &
        radau_coordinates(8) )
    call radau_angle( relative(2,:), relative(4,:), &
        radau_coordinates(2), radau_coordinates(4), &
        radau_coordinates(9) )
    call radau_angle( relative(3,:), relative(4,:), &
        radau_coordinates(3), radau_coordinates(4), &
        radau_coordinates(10) )

contains

    !> @brief Calculate one H-C-H angle from relative Cartesian vectors.
    subroutine radau_angle( first, second, first_norm, second_norm, angle )
        real(8), intent(in) :: first(3), second(3)
        real(8), intent(in) :: first_norm, second_norm
        real(8), intent(out) :: angle
        real(8) :: displacement(3), displacement_norm

        displacement = first-second
        displacement_norm = sqrt(sum(displacement**2))
        angle = acos((first_norm**2+second_norm**2-displacement_norm**2) &
            /(2.0d0*first_norm*second_norm))
    end subroutine radau_angle
end subroutine methaneradau
