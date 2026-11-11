!> @brief Evaluate the OH2_LEPS_JW1977 potential energy surface.
module chempotpy_oh2_leps_jw1977_pes
    use chempotpy_internal_coordinates, only: finalize_legacy_triatomic_surface, &
        & prepare_legacy_triatomic_coordinates
    use chempotpy_kinds, only: real64
    implicit none
    private
    public :: evaluate
    real(real64), parameter :: energy_per_hartree         = 627.5095_real64
    real(real64), parameter :: bohr_per_angstrom          = 0.529177106_real64
    real(real64), parameter :: root_two                   = 1.41421356_real64
    real(real64), parameter :: dissociation_energy(3)     = [106.6d0, 109.4d0, 106.6d0]/energy_per_hartree
    real(real64), parameter :: equilibrium_distance(3)    = [0.9706d0, 0.7417d0, 0.9706d0]/bohr_per_angstrom
    real(real64), parameter :: morse_exponent(3)          = [2.294d0, 1.942d0, 2.294d0]*bohr_per_angstrom
    real(real64), parameter :: sato_parameter(3)          = [0.0885d0, 0.0885d0, 0.0885d0]
    real(real64), parameter :: one_plus_sato(3)           = 1.0_real64+sato_parameter
    real(real64), parameter :: one_plus_three_sato(3)     = 1.0_real64+3.0_real64*sato_parameter
    real(real64), parameter :: two_one_plus_three_sato(3) = 2.0_real64*one_plus_three_sato
    real(real64), parameter :: sato_plus_three(3)         = sato_parameter+3.0_real64
    real(real64), parameter :: two_sato_plus_three(3)     = 2.0_real64*sato_plus_three
    real(real64), parameter :: quarter_energy(3)          = dissociation_energy/(4.0_real64*one_plus_sato)
    real(real64), parameter :: derivative_weight(3)       = 2.0_real64*morse_exponent*quarter_energy

    contains

    !> @brief Evaluate the OH2_LEPS_JW1977 potential energy surface for OH2.
    subroutine evaluate( x, igrad, p, g, d )
        integer, intent(in)       :: igrad
        real(real64), intent(in)  :: x(3,3)
        real(real64), intent(out) :: p(1), g(1,3,3), d(1,1,3,3)
        real(real64)              :: distance_jacobian(3,9), distances_bohr(3)
        real(real64)              :: energy_hartree, internal_gradient(3)
        call prepare_legacy_triatomic_coordinates(x, distances_bohr, distance_jacobian)
        call evaluate_internal(distances_bohr, igrad, energy_hartree, internal_gradient)
        call finalize_legacy_triatomic_surface(igrad, energy_hartree, internal_gradient, distance_jacobian, p, g, d)
    end subroutine evaluate

    !> @brief Evaluate the extended LEPS energy and optional internal-coordinate gradient.
    pure subroutine evaluate_internal( distances, derivative_order, energy, gradient )
        real(real64), intent(in)  :: distances(3)
        integer, intent(in)       :: derivative_order
        real(real64), intent(out) :: energy, gradient(3)
        real(real64)              :: coulomb(3), exchange(3), exponential(3)
        real(real64)              :: exchange_norm, exchange_sum
        integer                   :: bond
        energy   = 0.0_real64
        gradient = 0.0_real64
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
                gradient(bond) = derivative_weight(bond)*exponential(bond)*( (3.0_real64*exchange(bond)-exchange_sum)/ &
                    & root_two * (one_plus_three_sato(bond)*exponential(bond)-sato_plus_three(bond))/exchange_norm - &
                    & sato_plus_three(bond)*exponential(bond)+one_plus_three_sato(bond) )
            enddo
        endif
    end subroutine evaluate_internal
end module chempotpy_oh2_leps_jw1977_pes

!> @brief Expose OH2_LEPS_JW1977 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh2_leps_jw1977_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
