!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate the Cl2H_LEPS_PK2_1987 potential energy surface.
module chempotpy_cl2h_leps_pk2_1987_pes
    use chempotpy_internal_coordinates, only: finalize_legacy_triatomic_surface, &
        & prepare_legacy_triatomic_coordinates
    use chempotpy_kinds, only: real64
    implicit none
    private
    public :: evaluate
    real(real64), parameter :: energy_per_hartree         = 627.5095_real64
    real(real64), parameter :: bohr_per_angstrom          = 0.52917706_real64
    real(real64), parameter :: root_two                   = 1.41421356_real64
    logical, parameter      :: has_jws_correction         = .false.
    real(real64), parameter :: dissociation_energy(3)     = [106.477d0, 106.477d0, 57.983d0]/energy_per_hartree
    real(real64), parameter :: equilibrium_distance(3)    = [1.275d0, 1.275d0, 1.988d0]/bohr_per_angstrom
    real(real64), parameter :: morse_exponent(3)          = [1.868d0, 1.868d0, 2.002d0]*bohr_per_angstrom
    real(real64), parameter :: sato_parameter(3)          = [0.166d0, 0.166d0, -0.098d0]
    real(real64), parameter :: one_plus_sato(3)           = 1.0_real64+sato_parameter
    real(real64), parameter :: one_plus_three_sato(3)     = 1.0_real64+3.0_real64*sato_parameter
    real(real64), parameter :: two_one_plus_three_sato(3) = 2.0_real64*one_plus_three_sato
    real(real64), parameter :: sato_plus_three(3)         = sato_parameter+3.0_real64
    real(real64), parameter :: two_sato_plus_three(3)     = 2.0_real64*sato_plus_three
    real(real64), parameter :: quarter_energy(3)          = dissociation_energy/(4.0_real64*one_plus_sato)
    real(real64), parameter :: derivative_weight(3)       = 2.0_real64*morse_exponent*quarter_energy
    contains
    !> @brief Evaluate the surface through the common CHEMPOTPY Cartesian contract.
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
        real(real64)              :: correction, exchange_norm, exchange_sum, shift_1, shift_2, shift_3
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
        if( has_jws_correction )then
            shift_1    = distances(1)-2.105_real64
            shift_2    = distances(3)-2.105_real64
            shift_3    = distances(2)-4.21_real64
            correction = 0.04382_real64*exp(-shift_1*shift_1-shift_2*shift_2-0.5_real64*shift_3*shift_3)
            energy     = energy+correction
            if( derivative_order == 1 )then
                gradient(1) = gradient(1)-2.0_real64*shift_1*correction
                gradient(3) = gradient(3)-2.0_real64*shift_2*correction
                gradient(2) = gradient(2)-shift_3*correction
            endif
        endif
    end subroutine evaluate_internal
end module chempotpy_cl2h_leps_pk2_1987_pes

!> @brief Preserve the historical Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_cl2h_leps_pk2_1987_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
