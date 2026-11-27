! System:                        HeNO
! Functional form:               Rare-gas--linear-molecule analytic fit
! Common name:                   He-NO RCCSD(T)
! Reference:                     J. Klos et al., J. Chem. Phys. 112, 2195 (2000)
! Number of derivatives:         0
! Number of bodies:              3
! Number of electronic surfaces: 2
! Interface:                     PES
! Data file:                     none; immutable parameters are embedded
! Notes:                         Rigid-NO model; atom order He, N, O; states A-double-prime and A-prime

!> @brief Implement the two-component Klos He--NO(X 2Pi) interaction potential.
module chempotpy_heno_klos_2000_pes
    use chempotpy_kinds, only: real32, wp
    implicit none
    private
    public :: evaluate
    public :: evaluate_components
    integer, parameter :: natoms = 3, nstates = 2
    ! The release used default-real E literals for A-double-prime; retain that published rounding deliberately.
    real(wp), parameter :: a_double_prime_parameters(40) = real([ &
        & 0.93729727e+01_real32, -0.34394679_real32, 0.12075059e+01_real32, 0.87512578e-01_real32, &
        & 0.24040713e+00_real32, 0.14739815e+00_real32, -0.20236282e+01_real32, -0.44819367e-01_real32, &
        & -0.18935935e-01_real32, 0.43253788e-04_real32, -0.44905825e-01_real32, -0.42649430e-01_real32, &
        & -0.25182060e+01_real32, 0.21820896e+01_real32, -0.35647710e+01_real32, 0.56402415e+00_real32, &
        & 0.60600803e+01_real32, 0.24719691e+00_real32, 0.38623219e+01_real32, 0.26295718e+00_real32, &
        & 0.23270757e+01_real32, -0.10360823e+01_real32, -0.39459956e+01_real32, -0.67839402e+00_real32, &
        & -0.53872460e+00_real32, -0.11447268e+00_real32, -0.37394241e+00_real32, 0.24835212e+00_real32, &
        & 0.78682216e+00_real32, 0.24061793e+00_real32, 0.90072824e-02_real32, 0.39021151e-02_real32, &
        & 0.18131326e-01_real32, -0.17611395e-01_real32, -0.49279888e-01_real32, -0.19927465e-01_real32, &
        & -0.12539203e+05_real32, -0.47580271e+04_real32, 0.15366833e+05_real32, &
        & 0.51647119e+04_real32 ], wp)
    real(wp), parameter :: a_prime_parameters(40) = [ &
        & 0.89614790e+01_wp, -0.66184380e+00_wp, 0.78007803e+00_wp, -0.60974193e-02_wp, &
        & -0.35808998e-01_wp, -0.11214555e-01_wp, -0.20242288e+01_wp, -0.34243587e-01_wp, &
        & -0.33588331e-01_wp, -0.37086659e-02_wp, 0.77114573e-02_wp, 0.24071004e-02_wp, &
        & -0.24565438e+00_wp, 0.32383547e+01_wp, -0.18898698e+01_wp, 0.13834898e+01_wp, &
        & 0.23578296e+01_wp, 0.14503735e+01_wp, 0.43039404e+01_wp, 0.22754183e+01_wp, &
        & 0.47018602e+01_wp, 0.24381747e+00_wp, -0.25353828e+00_wp, -0.49050007e+00_wp, &
        & -0.44107679e+00_wp, -0.38335519e+00_wp, -0.55140365e+00_wp, 0.42936877e-01_wp, &
        & -0.19769691e+00_wp, 0.19501492e-01_wp, -0.74599920e-02_wp, -0.40805756e-02_wp, &
        & 0.23613677e-07_wp, -0.56473632e-02_wp, 0.11274139e-01_wp, -0.33655160e-02_wp, &
        & -0.13246868e+05_wp, -0.41907209e+04_wp, 0.20463763e+05_wp, 0.58529744e+03_wp ]

    contains

    !> @brief Evaluate the two adiabatic He--NO surfaces from Cartesian coordinates in angstrom.
    subroutine evaluate( coordinates, derivative_order, energy, gradient, coupling )
        use chempotpy_atomic_elements,      only: atomic_mass_u
        use chempotpy_internal_coordinates, only: rigid_diatom_jacobi_coordinates
        use chempotpy_units,                only: hartree_to_electronvolt
        integer, intent(in)   :: derivative_order
        real(wp), intent(in)  :: coordinates(natoms,3)
        real(wp), intent(out) :: energy(nstates), gradient(nstates,natoms,3)
        real(wp), intent(out) :: coupling(nstates,nstates,natoms,3)
        logical               :: valid
        real(wp)              :: cosine_theta, distance_bohr, v_difference_millihartree, v_sum_millihartree
        gradient = 0.0_wp; coupling = 0.0_wp
        call rigid_diatom_jacobi_coordinates(coordinates, [atomic_mass_u('N'), atomic_mass_u('O')], &
            distance_bohr, cosine_theta, valid)
        if( .not. valid )then
            energy = huge(1.0_wp)
            return
        endif
        call evaluate_components_millihartree(distance_bohr, cosine_theta, v_sum_millihartree, &
            v_difference_millihartree)
        energy(1) = (v_sum_millihartree+v_difference_millihartree)*1.0e-3_wp*hartree_to_electronvolt
        energy(2) = (v_sum_millihartree-v_difference_millihartree)*1.0e-3_wp*hartree_to_electronvolt
    end subroutine evaluate

    !> @brief Evaluate the released Vsum and Vdif components in inverse centimeters.
    pure subroutine evaluate_components( distance_bohr, cosine_theta, v_sum, v_difference )
        use chempotpy_units, only: hartree_to_inverse_centimeter
        real(wp), intent(in)  :: distance_bohr, cosine_theta
        real(wp), intent(out) :: v_sum, v_difference
        real(wp)              :: v_difference_millihartree, v_sum_millihartree
        call evaluate_components_millihartree(distance_bohr, cosine_theta, v_sum_millihartree, &
            v_difference_millihartree)
        v_sum        = v_sum_millihartree*1.0e-3_wp*hartree_to_inverse_centimeter
        v_difference = v_difference_millihartree*1.0e-3_wp*hartree_to_inverse_centimeter
    end subroutine evaluate_components

    !> @brief Evaluate the two fitted interaction components in millihartree.
    pure subroutine evaluate_components_millihartree( distance_bohr, cosine_theta, v_sum, v_difference )
        real(wp), intent(in)  :: distance_bohr, cosine_theta
        real(wp), intent(out) :: v_sum, v_difference
        real(wp)              :: a_double_prime, a_prime
        a_double_prime = evaluate_fitted_state(a_double_prime_parameters, distance_bohr, cosine_theta)
        a_prime        = evaluate_fitted_state(a_prime_parameters, distance_bohr, cosine_theta)
        v_sum          = 0.5_wp*(a_double_prime+a_prime)
        v_difference   = 0.5_wp*(a_double_prime-a_prime)
        if( abs(abs(cosine_theta)-1.0_wp) <= 32.0_wp*epsilon(1.0_wp) ) v_difference = 0.0_wp
    end subroutine evaluate_components_millihartree

    !> @brief Evaluate one Klos rare-gas--linear-molecule analytic state fit.
    pure function evaluate_fitted_state( parameters, distance_bohr, cosine_theta ) result(energy_millihartree)
        use chempotpy_dispersion, only: tang_toennies_damping
        real(wp), intent(in) :: parameters(40), distance_bohr, cosine_theta
        real(wp)             :: energy_millihartree
        real(wp)             :: angular(6), damping_6, damping_7, exponent_constant, exponent_slope
        real(wp)             :: long_range, radial_polynomial, short_range, t
        integer              :: index
        t                 = -cosine_theta
        angular(1)        = 1.0_wp
        angular(2)        = t
        angular(3)        = 0.5_wp*(3.0_wp*t**2-1.0_wp)
        angular(4)        = 0.5_wp*(5.0_wp*t**3-3.0_wp*t)
        angular(5)        = 0.125_wp*(35.0_wp*t**4-30.0_wp*t**2+3.0_wp)
        angular(6)        = 0.125_wp*(63.0_wp*t**5-70.0_wp*t**3+15.0_wp*t)
        exponent_constant = dot_product(parameters(1:6), angular)
        exponent_slope    = dot_product(parameters(7:12), angular)
        radial_polynomial = 0.0_wp
        do index=1, 6
            radial_polynomial = radial_polynomial+(parameters(12+index)+parameters(18+index)*distance_bohr+ &
                & parameters(24+index)*distance_bohr**2+parameters(30+index)*distance_bohr**3)*angular(index)
        enddo
        short_range = radial_polynomial*exp(exponent_constant+exponent_slope*distance_bohr)
        damping_6   = tang_toennies_damping(6, abs(exponent_slope*distance_bohr))
        damping_7   = tang_toennies_damping(7, abs(exponent_slope*distance_bohr))
        long_range  = damping_6*(parameters(37)+parameters(38)*angular(3))/distance_bohr**6+ damping_7* &
            & (parameters(39)*angular(2)+parameters(40)*angular(4))/distance_bohr**7
        energy_millihartree = short_range+long_range
    end function evaluate_fitted_state
end module chempotpy_heno_klos_2000_pes

!> @brief Expose the rigid-rotor HeNO_Klos_2000 potential through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_heno_klos_2000_pes, only: evaluate
    use chempotpy_kinds,              only: wp
    implicit none
    integer, parameter    :: natoms = 3, nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
