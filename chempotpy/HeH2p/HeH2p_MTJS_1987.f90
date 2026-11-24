! System:                        HeH2p
! Functional form:               Carter-Murrell-Sorbie many-body expansion
! Common name:                   He + H2+ MTJS
! Reference:                     T. Joseph and N. Sathyamurthy, J. Chem. Phys. 86, 704 (1987)
! Number of derivatives:         0
! Number of bodies:              3
! Number of electronic surfaces: 1
! Interface:                     PES
! Data file:                     none; immutable parameters are embedded
! Notes:                         The atom order is H, H, He. The released fit is energy-only.

!> @brief Evaluate the Joseph-Sathyamurthy many-body surface for He + H2+.
module chempotpy_heh2p_mtjs_1987_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    real(wp), parameter :: polynomial_coefficients(26) = [ &
        6.34650_wp, 9.34530_wp, -3.63380_wp, 5.17180_wp, -2.96930_wp, -0.08740_wp, &
        6.93820_wp, 0.61410_wp, 4.78700_wp, 5.88260_wp, -9.74530_wp, -3.57910_wp, &
        4.15630_wp, -1.95470_wp, -0.45990_wp, 6.43060_wp, -4.00110_wp, 4.28260_wp, &
        -3.80790_wp, 7.96950_wp, -0.12490_wp, 0.20040_wp, 2.22250_wp, 1.37490_wp, &
        0.95250_wp, 0.95250_wp]
    real(wp), parameter :: dissociation_energy(3)  = [2.78650_wp, 2.02790_wp, 2.02790_wp]
    real(wp), parameter :: equilibrium_distance(3) = [1.05360_wp, 0.77600_wp, 0.77600_wp]
    real(wp), parameter :: diatomic_shape(3,3)     = transpose(reshape([ &
        2.61800_wp, 1.57460_wp, 1.24510_wp, &
        4.60930_wp, 2.94020_wp, 2.61850_wp, &
        4.60930_wp, 2.94020_wp, 2.61850_wp], [3, 3]))

    contains

    !> @brief Evaluate the energy in eV for Cartesian coordinates in angstrom.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_geometry, only: cartesian_distance
        use chempotpy_status,   only: chempotpy_status_t, status_derivative_unavailable
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        integer, parameter    :: natoms  = 3
        integer, parameter    :: nstates = 1
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: distances(3)
        p = 0.0_wp; g = 0.0_wp; d = 0.0_wp
        if( igrad /= 0 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'HeH2p_MTJS_1987 provides energies but not analytic derivatives')
            return
        endif
        distances(1) = cartesian_distance(x(1,:), x(2,:))
        distances(2) = cartesian_distance(x(2,:), x(3,:))
        distances(3) = cartesian_distance(x(1,:), x(3,:))
        p(1)         = evaluate_internal(distances)
    end subroutine evaluate

    !> @brief Evaluate the published many-body expansion in its native angstrom/eV units.
    pure function evaluate_internal( distances ) result(energy)
        real(wp), intent(in) :: distances(3)
        real(wp)             :: energy
        real(wp)             :: attenuation, diatomic_energy(3), polynomial, radial_displacement(3)
        real(wp)             :: rho1, rho2, rho3, rho1_squared, rho2_squared, rho3_squared
        real(wp)             :: rho12, rho13, rho23, term1, term2, term3, term4, term5
        real(wp)             :: term6, term7, term8, term9, term10, term11
        integer              :: pair
        radial_displacement = distances-equilibrium_distance
        do pair = 1, 3
            diatomic_energy(pair) = -dissociation_energy(pair)*(1.0_wp + diatomic_shape(pair,1)* &
                & radial_displacement(pair) + diatomic_shape(pair,2)*radial_displacement(pair)**2 + &
                & diatomic_shape(pair,3)*radial_displacement(pair)**3)* exp(-diatomic_shape(pair,1)* &
                & radial_displacement(pair))
        enddo
        rho1         = distances(1)-polynomial_coefficients(25)
        rho2         = distances(2)-polynomial_coefficients(26)
        rho3         = distances(3)-polynomial_coefficients(26)
        rho1_squared = rho1*rho1; rho2_squared = rho2*rho2; rho3_squared = rho3*rho3
        rho12        = rho1*rho2; rho13 = rho1*rho3; rho23 = rho2*rho3
        term1 = polynomial_coefficients(1) + polynomial_coefficients(2)*rho1 + polynomial_coefficients(3)*(rho2+ &
            & rho3) + polynomial_coefficients(4)*rho1_squared
        term2 = polynomial_coefficients(5)*(rho2_squared+rho3_squared) + polynomial_coefficients(6)*(rho12+rho13) + &
            & polynomial_coefficients(7)*rho23
        term3 = polynomial_coefficients(8)*rho1*rho1_squared + polynomial_coefficients(9)*(rho2*rho2_squared+rho3* &
            & rho3_squared)
        term4 = polynomial_coefficients(10)*rho1_squared*(rho2+rho3) + polynomial_coefficients(11)*rho1*(rho2_squared+ &
            & rho3_squared)
        term5 = polynomial_coefficients(12)*rho2*(rho23+rho3_squared) + polynomial_coefficients(13)*rho12*rho3
        term6 = polynomial_coefficients(14)*rho1_squared**2 + polynomial_coefficients(15)*(rho2_squared**2+ &
            & rho3_squared**2)
        term7 = polynomial_coefficients(16)*rho1_squared*(rho12+rho13) + polynomial_coefficients(17)*rho1_squared* &
            & (rho2_squared+rho3_squared)
        term8 = polynomial_coefficients(18)*(rho12*rho2_squared+rho13*rho3_squared)
        term9 = polynomial_coefficients(19)*rho2*(rho2_squared*rho3+rho3_squared*rho3)
        term10     = polynomial_coefficients(20)*rho2_squared*rho3_squared + polynomial_coefficients(21)*rho1_squared* &
            & rho23
        term11     = polynomial_coefficients(22)*(rho13*rho2_squared+rho12*rho3_squared)
        polynomial = term1+term2+term3+term4+term5+term6+term7+term8+term9+term10+term11
        attenuation = (1.0_wp-tanh(0.5_wp*polynomial_coefficients(23)*rho1))* (1.0_wp-tanh(0.5_wp* &
            & polynomial_coefficients(24)*(rho2+rho3)))
        energy = polynomial*attenuation+sum(diatomic_energy)+dissociation_energy(1)
    end function evaluate_internal
end module chempotpy_heh2p_mtjs_1987_pes

!> @brief Expose HeH2p_MTJS_1987 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_heh2p_mtjs_1987_pes, only: evaluate
    use chempotpy_kinds,               only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
