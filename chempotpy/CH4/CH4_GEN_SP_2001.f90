!> @brief Implement CH4_GEN_SP_2001, a potential energy surface for CH4.
module chempotpy_ch4_gen_sp_2001_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the Schwenke-Partridge methane potential.
    subroutine evaluate( x, igrad, path, p, g, d )
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter              :: nstates = 1
        integer, parameter              :: natoms  = 5
        real(wp), intent(in)            :: x(natoms,3)
        integer, intent(in)             :: igrad
        character(len=1024), intent(in) :: path
        real(wp), intent(out)           :: p(nstates)
        real(wp), intent(out)           :: g(nstates,natoms,3)
        real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
        real(wp)                        :: coordinates_bohr(natoms,3)
        real(wp)                        :: energy_hartree
        p                = 0.0e0_wp
        g                = 0.0e0_wp
        d                = 0.0e0_wp
        coordinates_bohr = x/bohr_to_angstrom
        call getpot(energy_hartree, coordinates_bohr, path)
        if( igrad == 0 )then
            p(1) = energy_hartree*hartree_to_electronvolt
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Convert methane Cartesian coordinates and call the T8 kernel.
    subroutine getpot( potential_energy, coordinates, path )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), parameter :: zero_point_shift = &
            -28.11193019e0_wp*4.556335e-6_wp
        character(len=1024), intent(in) :: path
        real(wp), intent(in)            :: coordinates(5,3)
        real(wp), intent(out)           :: potential_energy
        real(wp)                        :: energy_hartree
        real(wp)                        :: radau_coordinates(10)
        call methaneradau(coordinates, radau_coordinates)
        call vibpot(radau_coordinates, energy_hartree, 1, path)
        potential_energy = energy_hartree-zero_point_shift
    end subroutine getpot

    !> @brief Convert CH4 Cartesian coordinates to the T8 Radau coordinates.
    subroutine methaneradau( coordinates, radau_coordinates )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), parameter   :: unitconvm = 5.48579903e-4_wp
        real(wp), intent(in)  :: coordinates(5,3)
        real(wp), intent(out) :: radau_coordinates(10)
        real(wp)              :: alpha0
        real(wp)              :: carbon_mass, hydrogen_mass
        real(wp)              :: relative(4,3), shifted(4,3), shift(3)
        integer               :: atom_index
        carbon_mass   = 12.0e0_wp/unitconvm-6.0e0_wp
        hydrogen_mass = 1.007825e0_wp/unitconvm-1.0e0_wp
        alpha0        = sqrt(carbon_mass/(carbon_mass+4.0e0_wp*hydrogen_mass))
        do atom_index = 1, 4
            relative(atom_index,:) = coordinates(atom_index+1,:)-coordinates(1,:)
        enddo
        shift = 0.25e0_wp*(1.0e0_wp-alpha0)*sum(relative, dim=1)
        do atom_index = 1, 4
            shifted(atom_index,:)         = relative(atom_index,:)-shift
            radau_coordinates(atom_index) = sqrt(sum(shifted(atom_index,:)**2))
        enddo
        call radau_angle(relative(1,:), relative(2,:), &
            radau_coordinates(1), radau_coordinates(2), &
            radau_coordinates(5))
        call radau_angle(relative(1,:), relative(3,:), &
            radau_coordinates(1), radau_coordinates(3), &
            radau_coordinates(6))
        call radau_angle(relative(1,:), relative(4,:), &
            radau_coordinates(1), radau_coordinates(4), &
            radau_coordinates(7))
        call radau_angle(relative(2,:), relative(3,:), &
            radau_coordinates(2), radau_coordinates(3), &
            radau_coordinates(8))
        call radau_angle(relative(2,:), relative(4,:), &
            radau_coordinates(2), radau_coordinates(4), &
            radau_coordinates(9))
        call radau_angle(relative(3,:), relative(4,:), &
            radau_coordinates(3), radau_coordinates(4), &
            radau_coordinates(10))

        contains

        !> @brief Calculate one H-C-H angle from relative Cartesian vectors.
        subroutine radau_angle( first, second, first_norm, second_norm, angle )
            real(wp), intent(in)  :: first(3), second(3)
            real(wp), intent(in)  :: first_norm, second_norm
            real(wp), intent(out) :: angle
            real(wp)              :: displacement(3), displacement_norm
            displacement      = first-second
            displacement_norm = sqrt(sum(displacement**2))
            angle             = acos((first_norm**2+second_norm**2-displacement_norm**2) /(2.0e0_wp*first_norm* &
                & second_norm))
        end subroutine radau_angle
    end subroutine methaneradau
end module chempotpy_ch4_gen_sp_2001_pes

!> @brief Expose CH4_GEN_SP_2001 through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_ch4_gen_sp_2001_pes, only: evaluate
    use chempotpy_kinds,               only: wp
    implicit none
    integer, parameter              :: natoms  = 5
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, path, p, g, d)
end subroutine pes
