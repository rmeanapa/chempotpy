!> @brief Evaluate the BrHCl_LEPS1 potential energy surface.
module chempotpy_brhcl_leps1_pes
    use chempotpy_internal_coordinates, only: finalize_legacy_triatomic_surface, &
        & prepare_legacy_triatomic_coordinates
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    real(wp), parameter :: energy_per_hartree         = 627.5095_wp
    real(wp), parameter :: bohr_per_angstrom          = 0.529177106_wp
    real(wp), parameter :: root_two                   = 1.41421356_wp
    real(wp), parameter :: dissociation_energy(3)     = [90.24e0_wp, 106.4e0_wp, 52.09e0_wp]/energy_per_hartree
    real(wp), parameter :: equilibrium_distance(3)    = [1.414e0_wp, 1.275e0_wp, 2.136e0_wp]/bohr_per_angstrom
    real(wp), parameter :: morse_exponent(3)          = [1.851e0_wp, 1.867e0_wp, 1.923e0_wp]*bohr_per_angstrom
    real(wp), parameter :: sato_parameter(3)          = [0.154e0_wp, 0.154e0_wp, 0.154e0_wp]
    real(wp), parameter :: one_plus_sato(3)           = 1.0_wp+sato_parameter
    real(wp), parameter :: one_plus_three_sato(3)     = 1.0_wp+3.0_wp*sato_parameter
    real(wp), parameter :: two_one_plus_three_sato(3) = 2.0_wp*one_plus_three_sato
    real(wp), parameter :: sato_plus_three(3)         = sato_parameter+3.0_wp
    real(wp), parameter :: two_sato_plus_three(3)     = 2.0_wp*sato_plus_three
    real(wp), parameter :: quarter_energy(3)          = dissociation_energy/(4.0_wp*one_plus_sato)
    real(wp), parameter :: derivative_weight(3)       = 2.0_wp*morse_exponent*quarter_energy

    contains

    !> @brief Evaluate the BrHCl_LEPS1 potential energy surface for BrHCl.
    subroutine evaluate( x, igrad, p, g, d )
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(3,3)
        real(wp), intent(out) :: p(1), g(1,3,3), d(1,1,3,3)
        real(wp)              :: distance_jacobian(3,9), distances_bohr(3)
        real(wp)              :: energy_hartree, internal_gradient(3)
        call prepare_legacy_triatomic_coordinates(x, distances_bohr, distance_jacobian)
        call evaluate_internal(distances_bohr, igrad, energy_hartree, internal_gradient)
        call finalize_legacy_triatomic_surface(igrad, energy_hartree, internal_gradient, distance_jacobian, p, g, d)
    end subroutine evaluate

    !> @brief Evaluate the extended LEPS energy and optional internal-coordinate gradient.
    pure subroutine evaluate_internal( distances, derivative_order, energy, gradient )
        real(wp), intent(in)  :: distances(3)
        integer, intent(in)   :: derivative_order
        real(wp), intent(out) :: energy, gradient(3)
        real(wp)              :: coulomb(3), exchange(3), exponential(3)
        real(wp)              :: exchange_norm, exchange_sum
        integer               :: bond
        energy   = 0.0_wp
        gradient = 0.0_wp
        do bond = 1, 3
            exponential(bond) = exp(-morse_exponent(bond)*(distances(bond)-equilibrium_distance(bond)))
            coulomb(bond) = quarter_energy(bond)*(sato_plus_three(bond)*exponential(bond)- &
                & two_one_plus_three_sato(bond))*exponential(bond)
            exchange(bond) = quarter_energy(bond)*(one_plus_three_sato(bond)*exponential(bond)- &
                & two_sato_plus_three(bond))*exponential(bond)
            energy = energy+coulomb(bond)
        enddo
        exchange_norm = sqrt((exchange(1)-exchange(2))**2 + (exchange(2)-exchange(3))**2 + (exchange(3)-exchange(1))**2)
        energy        = energy-exchange_norm/root_two+dissociation_energy(2)
        if( derivative_order == 1 )then
            exchange_sum = sum(exchange)
            do bond = 1, 3
                gradient(bond) = derivative_weight(bond)*exponential(bond)*( (3.0_wp*exchange(bond)-exchange_sum)/ &
                    & root_two * (one_plus_three_sato(bond)*exponential(bond)-sato_plus_three(bond))/exchange_norm - &
                    & sato_plus_three(bond)*exponential(bond)+one_plus_three_sato(bond) )
            enddo
        endif
    end subroutine evaluate_internal
end module chempotpy_brhcl_leps1_pes

!> @brief Expose BrHCl_LEPS1 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_brhcl_leps1_pes, only: evaluate
    use chempotpy_kinds,           only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
