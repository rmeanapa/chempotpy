! System:                        KrNO
! Functional form:               Rare-gas--diatom analytic radial/angular fit
! Common name:                   Kr-NO(X 2Pi)
! Reference:                     B. Wen et al., J. Phys. Chem. A 113, 7366--7375 (2009); DOI 10.1021/jp811513j
! Number of derivatives:         0
! Number of bodies:              3
! Number of electronic surfaces: 2
! Interface:                     PES
! Data file:                     none; immutable parameters are embedded
! Notes:                         Released at DOI 10.5281/zenodo.14509374; atom order Kr, N, O

!> @brief Implement the two-component Klos Kr--NO(X 2Pi) interaction potential.
module chempotpy_krno_klos_2009_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    public :: evaluate_components
    integer, parameter :: natoms = 3, nstates = 2
    ! Parameters are transcribed from the CC BY 4.0 release at DOI 10.5281/zenodo.14509374.
    real(wp), parameter :: a_prime_parameters(16,7) = reshape( [ &
        & 1.25717367691647652e+00_wp, -9.78788456532057793e+00_wp, 8.09782459815127496e+06_wp, &
        & -6.39172054113550857e+06_wp, 2.26418936169679975e+06_wp, -4.62915169840692368e+05_wp, &
        & 5.93745330221700788e+04_wp, -4.87470347137828958e+03_wp, 2.49887894195171725e+02_wp, &
        & -7.31566968806472495e+00_wp, 9.40572848879226986e-02_wp, 5.01466237891449600e+03_wp, &
        & 7.26629065829550773e+07_wp, 1.07510322343766665e+09_wp, -3.77273093363379761e+11_wp, &
        & 9.75484777759294141e+12_wp, 1.21709560373223069e+00_wp, -6.48744687364050243e+00_wp, &
        & 7.13342089001057483e+06_wp, -5.82025734719443321e+06_wp, 2.15785530873638671e+06_wp, &
        & -4.60420332186074695e+05_wp, 6.09875287599588555e+04_wp, -5.10545749419488038e+03_wp, &
        & 2.63466875881902524e+02_wp, -7.67493607114191878e+00_wp, 9.71001840139579936e-02_wp, &
        & 5.01466237891449600e+03_wp, 6.77881314838012606e+07_wp, 1.02451091368218958e+09_wp, &
        & -6.40449213019956970e+09_wp, 2.62675971067293121e+11_wp, 1.02259705255180311e+00_wp, &
        & -4.22931059819828281e+00_wp, 8.41298547885173745e+06_wp, -6.31332264600359369e+06_wp, &
        & 2.18610436759517575e+06_wp, -4.42547293676054920e+05_wp, 5.61787256064376488e+04_wp, &
        & -4.52345794400489103e+03_wp, 2.24088668974039706e+02_wp, -6.22606462312523501e+00_wp, &
        & 7.43650128608139144e-02_wp, 5.01466237891449600e+03_wp, 7.18144856596150994e+07_wp, &
        & 6.43028416366354823e+08_wp, 1.28439975235925808e+10_wp, 9.33159597162069511e+09_wp, &
        & 1.10047387161036681e+00_wp, -4.17477036260108214e+00_wp, 1.05253742552225571e+07_wp, &
        & -8.45172598987702467e+06_wp, 3.07399328399921861e+06_wp, -6.46167854240398156e+05_wp, &
        & 8.45963621931602393e+04_wp, -6.99818251716040140e+03_wp, 3.55364580400691864e+02_wp, &
        & -1.01045290571688575e+01_wp, 1.23290384491220054e-01_wp, 5.01466237891449600e+03_wp, &
        & 7.64554333821597844e+07_wp, 3.26477245861719072e+08_wp, 2.74662617013192415e+09_wp, &
        & 2.26057768146365738e+10_wp, 1.04859173931842342e+00_wp, -1.99361355561161768e+00_wp, &
        & 1.98095856103767604e+07_wp, -1.04318662799394950e+07_wp, 2.63400242953600083e+06_wp, &
        & -4.35493335130205785e+05_wp, 5.16171332865437726e+04_wp, -4.28008622962904428e+03_wp, &
        & 2.29231045500385875e+02_wp, -6.99723534295297611e+00_wp, 9.21731633760895924e-02_wp, &
        & 5.01466237891449600e+03_wp, 7.88549245680361092e+07_wp, 3.40226481958905697e+08_wp, &
        & 2.88363815375071907e+09_wp, -4.74655194247720432e+09_wp, 1.72323714588573251e+00_wp, &
        & -1.04335137804348985e+01_wp, 2.96556750017641392e+06_wp, -3.21329979630668834e+06_wp, &
        & 1.50522907864027191e+06_wp, -3.97958929976709886e+05_wp, 6.50979126299272102e+04_wp, &
        & -6.75637220520657047e+03_wp, 4.34363677679029820e+02_wp, -1.57816077242744228e+01_wp, &
        & 2.46768575183440336e-01_wp, 5.01466237891449600e+03_wp, 7.94063926321976483e+07_wp, &
        & 5.36145712996094406e+08_wp, 2.84940809407146025e+09_wp, -4.99354056454523773e+10_wp, &
        & 2.08201600110828355e+00_wp, -1.51212984507481583e+01_wp, -1.76033262758345949e+06_wp, &
        & 2.05200832243337925e+06_wp, -1.03799405725804251e+06_wp, 2.97777292732825445e+05_wp, &
        & -5.30399138263440254e+04_wp, 6.01607821808363860e+03_wp, -4.25243450145618453e+02_wp, &
        & 1.71698260554314395e+01_wp, -3.04211320658061113e-01_wp, 5.01466237891449600e+03_wp, &
        & 8.24480374038421363e+07_wp, 5.87850638550512791e+08_wp, 5.10532235086869717e+09_wp, &
        & -9.72015733737625427e+10_wp ], [16, 7])
    real(wp), parameter :: a_double_prime_parameters(16,7) = reshape( [ &
        & 1.25717367691647652e+00_wp, -9.78788456532057793e+00_wp, 8.09782459815127496e+06_wp, &
        & -6.39172054113550857e+06_wp, 2.26418936169679975e+06_wp, -4.62915169840692368e+05_wp, &
        & 5.93745330221700788e+04_wp, -4.87470347137828958e+03_wp, 2.49887894195171725e+02_wp, &
        & -7.31566968806472495e+00_wp, 9.40572848879226986e-02_wp, 5.01466237891449600e+03_wp, &
        & 7.26629065829550773e+07_wp, 1.07510322343766665e+09_wp, -3.77273093363379761e+11_wp, &
        & 9.75484777759294141e+12_wp, 1.21391485484413342e+00_wp, -6.18514112822369011e+00_wp, &
        & 7.17890078289272264e+06_wp, -5.83068988088009879e+06_wp, 2.16438908965821005e+06_wp, &
        & -4.63025207500261371e+05_wp, 6.14099921976343321e+04_wp, -5.13406755209529911e+03_wp, &
        & 2.63825178834662211e+02_wp, -7.63178950393744415e+00_wp, 9.56585762283299890e-02_wp, &
        & 5.01466237891449600e+03_wp, 7.40978586973972172e+07_wp, 8.73275920214649796e+08_wp, &
        & -5.98873426094741821e+09_wp, 2.17132275411089661e+11_wp, 1.30357209222345216e+00_wp, &
        & -6.42528160139274362e+00_wp, -3.40570866048953403e+06_wp, 3.43700423190768762e+06_wp, &
        & -1.49536593043263350e+06_wp, 3.69980133881654358e+05_wp, -5.68690538999575510e+04_wp, &
        & 5.53045932775761958e+03_wp, -3.30440794088726136e+02_wp, 1.10388424491119537e+01_wp, &
        & -1.57473242110801187e-01_wp, 5.01466237891449600e+03_wp, 8.00033478218421191e+07_wp, &
        & 4.41155671826940060e+08_wp, -3.86695984158222040e+05_wp, -2.95907384086639061e+10_wp, &
        & 1.10373175851425009e+00_wp, -3.82993547997363581e+00_wp, 1.27714755976436660e+07_wp, &
        & -1.02072680012995359e+07_wp, 3.70243061476593651e+06_wp, -7.76374550563270692e+05_wp, &
        & 1.01234301709989159e+05_wp, -8.31663728070779689e+03_wp, 4.17695763365515802e+02_wp, &
        & -1.16860477646175518e+01_wp, 1.39342212810863453e-01_wp, 5.01466237891449600e+03_wp, &
        & 8.27257615260146111e+07_wp, 1.89140044924528003e+08_wp, -1.18829920565854505e+08_wp, &
        & 2.83523646792802353e+10_wp, 9.48352654833558129e-01_wp, -1.31206759701059261e+00_wp, &
        & 3.10028020268645100e+07_wp, -1.67780614492693208e+07_wp, 3.99350815849820059e+06_wp, &
        & -5.55371930208534235e+05_wp, 5.06363478304414020e+04_wp, -3.17846487106930954e+03_wp, &
        & 1.35182121272924405e+02_wp, -3.50124218707751611e+00_wp, 4.11234615852565863e-02_wp, &
        & 5.01466237891449600e+03_wp, 8.58183114753690958e+07_wp, 2.01112935912792683e+08_wp, &
        & 1.90067088858048129e+09_wp, -1.60232331911966596e+07_wp, 2.25936440264973415e+00_wp, &
        & -1.59699880768872440e+01_wp, -1.94974339877647883e+06_wp, 2.42250570941398107e+06_wp, &
        & -1.30913050736845192e+06_wp, 4.02031609402810456e+05_wp, -7.67743489451826172e+04_wp, &
        & 9.34230093326465976e+03_wp, -7.08003686438776299e+02_wp, 3.05838069890516699e+01_wp, &
        & -5.77450976022965246e-01_wp, 5.01466237891449600e+03_wp, 8.42653888527839333e+07_wp, &
        & 4.40035742067516506e+08_wp, 2.04014962315189219e+09_wp, -3.82291448870839081e+10_wp, &
        & 2.08201600110828355e+00_wp, -1.51212984507481583e+01_wp, -1.76033262758345949e+06_wp, &
        & 2.05200832243337925e+06_wp, -1.03799405725804251e+06_wp, 2.97777292732825445e+05_wp, &
        & -5.30399138263440254e+04_wp, 6.01607821808363860e+03_wp, -4.25243450145618453e+02_wp, &
        & 1.71698260554314395e+01_wp, -3.04211320658061113e-01_wp, 5.01466237891449600e+03_wp, &
        & 8.24480374038421363e+07_wp, 5.87850638550512791e+08_wp, 5.10532235086869717e+09_wp, &
        & -9.72015733737625427e+10_wp ], [16, 7])
    real(wp), parameter :: sum_angular_inverse(7,7) = reshape( [ &
        & 1.42857142857142627e-02_wp, 4.28571428571427882e-02_wp, 7.93650793650791664e-02_wp, &
        & 1.18518518518518273e-01_wp, 2.21645021645021617e-01_wp, 3.38624338624338939e-01_wp, &
        & 1.84704184704184926e-01_wp, 1.26984126984126811e-01_wp, 3.29914439536929061e-01_wp, &
        & 3.80952380952381153e-01_wp, 2.56600119639833923e-01_wp, -1.38528138528138278e-01_wp, &
        & -5.86514559176762984e-01_wp, -3.69408369408369797e-01_wp, 2.28571428571428537e-01_wp, &
        & 3.42857142857142638e-01_wp, -1.26984126984126977e-01_wp, -6.81481481481481466e-01_wp, &
        & -4.70995670995671134e-01_wp, 3.38624338624338828e-01_wp, 3.69408369408369575e-01_wp, &
        & 2.60317460317460225e-01_wp, 2.90155171092458558e-16_wp, -6.66666666666666519e-01_wp, &
        & -3.55271367880050113e-16_wp, 7.75757575757575646e-01_wp, 6.51161967875915547e-17_wp, &
        & -3.69408369408369408e-01_wp, 2.28571428571428592e-01_wp, -3.42857142857142916e-01_wp, &
        & -1.26984126984127033e-01_wp, 6.81481481481481688e-01_wp, -4.70995670995670912e-01_wp, &
        & -3.38624338624338772e-01_wp, 3.69408369408369353e-01_wp, 1.26984126984127144e-01_wp, &
        & -3.29914439536929172e-01_wp, 3.80952380952380987e-01_wp, -2.56600119639833479e-01_wp, &
        & -1.38528138528138833e-01_wp, 5.86514559176762651e-01_wp, -3.69408369408369297e-01_wp, &
        & 1.42857142857142072e-02_wp, -4.28571428571427604e-02_wp, 7.93650793650792774e-02_wp, &
        & -1.18518518518518523e-01_wp, 2.21645021645021839e-01_wp, -3.38624338624338717e-01_wp, &
        & 1.84704184704184676e-01_wp ], [7, 7])
    real(wp), parameter :: difference_angular_inverse(5,5) = reshape( [ &
        & 1.03682105514631723e-01_wp, 2.81091347570522510e-01_wp, 5.69483769103915738e-01_wp, &
        & 8.49940344415569204e-01_wp, 5.04707916514797184e-01_wp, 5.18410527573159641e-01_wp, &
        & 8.11440825933579579e-01_wp, 3.94257993995019040e-01_wp, -4.90713286643452040e-01_wp, &
        & -5.04707916514797295e-01_wp, 8.03536317738396955e-01_wp, 1.82775953563301890e-16_wp, &
        & -8.76128875544486596e-01_wp, -3.77449335971199863e-16_wp, 5.04707916514797073e-01_wp, &
        & 5.18410527573159197e-01_wp, -8.11440825933579357e-01_wp, 3.94257993995018485e-01_wp, &
        & 4.90713286643451652e-01_wp, -5.04707916514796517e-01_wp, 1.03682105514632195e-01_wp, &
        & -2.81091347570523009e-01_wp, 5.69483769103916182e-01_wp, -8.49940344415568760e-01_wp, &
        & 5.04707916514796517e-01_wp ], [5, 5])

    contains

    !> @brief Evaluate the rigid-NO A-double-prime and A-prime surfaces from Cartesian coordinates.
    subroutine evaluate( coordinates, derivative_order, energy, gradient, coupling )
        use chempotpy_atomic_elements,      only: atomic_mass_u
        use chempotpy_internal_coordinates, only: rigid_diatom_jacobi_coordinates
        use chempotpy_units,                only: electronvolt_to_inverse_centimeter
        integer, intent(in)   :: derivative_order
        real(wp), intent(in)  :: coordinates(natoms,3)
        real(wp), intent(out) :: energy(nstates), gradient(nstates,natoms,3)
        real(wp), intent(out) :: coupling(nstates,nstates,natoms,3)
        logical               :: valid
        real(wp)              :: cosine_theta, distance_bohr, v_difference, v_sum
        gradient = 0.0_wp; coupling = 0.0_wp
        call rigid_diatom_jacobi_coordinates(coordinates, [atomic_mass_u('N'), atomic_mass_u('O')], &
            distance_bohr, cosine_theta, valid)
        if( .not. valid )then
            energy = huge(1.0_wp)
            return
        endif
        call evaluate_components(distance_bohr, cosine_theta, v_sum, v_difference)
        energy = [v_sum+v_difference, v_sum-v_difference]/electronvolt_to_inverse_centimeter
    end subroutine evaluate

    !> @brief Evaluate the released Kr--NO Vsum and Vdif components in inverse centimeters.
    pure subroutine evaluate_components( distance_bohr, cosine_theta, v_sum, v_difference )
        use chempotpy_special_functions, only: associated_legendre_polynomial
        real(wp), intent(in)  :: distance_bohr, cosine_theta
        real(wp), intent(out) :: v_sum, v_difference
        real(wp)              :: basis_difference(5), basis_sum(7), coefficients_difference(5)
        real(wp)              :: coefficients_sum(7), difference_values(5), state_values(2,7), sum_values(7)
        integer               :: angle_index, degree
        do angle_index=1, 7
            state_values(1,angle_index) = evaluate_radial_fit(a_double_prime_parameters(:,angle_index), distance_bohr)
            state_values(2,angle_index) = evaluate_radial_fit(a_prime_parameters(:,angle_index), distance_bohr)
            sum_values(angle_index)     = 0.5_wp*sum(state_values(:,angle_index))
        enddo
        difference_values       = 0.5_wp*(state_values(1,2:6)-state_values(2,2:6))
        coefficients_sum        = matmul(sum_angular_inverse, sum_values)
        coefficients_difference = matmul(difference_angular_inverse, difference_values)
        do degree=0, 6
            basis_sum(degree+1) = associated_legendre_polynomial(degree, 0, cosine_theta)
        enddo
        do degree=2, 6
            basis_difference(degree-1) = associated_legendre_polynomial(degree, 2, cosine_theta)/ &
                & sqrt(real((degree-1)*degree*(degree+1)*(degree+2), wp))
        enddo
        v_sum        = dot_product(coefficients_sum, basis_sum)
        v_difference = dot_product(coefficients_difference, basis_difference)
    end subroutine evaluate_components

    !> @brief Evaluate one published angle-specific Kr--NO radial fit.
    pure function evaluate_radial_fit( parameters, input_distance_bohr ) result(energy_inverse_centimeter)
        real(wp), intent(in) :: parameters(16), input_distance_bohr
        real(wp)             :: energy_inverse_centimeter
        real(wp)             :: distance_bohr, dispersion, polynomial, switching
        integer              :: power
        distance_bohr = max(4.5_wp, input_distance_bohr)
        polynomial    = 0.0_wp
        do power=0, 8
            polynomial = polynomial+parameters(power+3)*distance_bohr**power
        enddo
        switching = 0.5_wp*(1.0_wp+tanh(1.0_wp+parameters(12)*distance_bohr))
        dispersion = parameters(13)/distance_bohr**6+parameters(14)/distance_bohr**7 + parameters(15)/distance_bohr** &
            & 8+parameters(16)/distance_bohr**9
        ! Preserve the release's rounded microhartree-to-wavenumber factor for numerical compatibility.
        energy_inverse_centimeter = (exp(-parameters(1)*distance_bohr-parameters(2))*polynomial - switching* &
            & dispersion)/4.556335_wp
    end function evaluate_radial_fit
end module chempotpy_krno_klos_2009_pes

!> @brief Expose KrNO_Klos_2009 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_kinds,              only: wp
    use chempotpy_krno_klos_2009_pes, only: evaluate
    implicit none
    integer, parameter    :: natoms = 3, nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3), d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
