!> @brief Provide the D3(BJ) pairwise dispersion correction used by PES kernels.
module chempotpy_dispersion
    use chempotpy_kinds, only: real32, wp
    use chempotpy_units, only: hartree_to_kilocalorie_per_mole
    implicit none
    private
    public :: evaluate_d3_bj_dispersion
    public :: evaluate_legacy_d3_diatomic
    public :: evaluate_legacy_d3_o3
    public :: evaluate_legacy_d3_o4
    public :: tang_toennies_damping
    real(wp), parameter :: d3_bohr_to_angstrom = 0.52917726_wp
    real(wp), parameter :: r2r4(94)            = [ &
        & 2.00734898_wp, 1.56637132_wp, 5.01986934_wp, 3.85379032_wp, 3.64446594_wp, 3.10492822_wp, &
        & 2.71175247_wp, 2.59361680_wp, 2.38825250_wp, 2.21522516_wp, 6.58585536_wp, 5.46295967_wp, &
        & 5.65216669_wp, 4.88284902_wp, 4.29727576_wp, 4.04108902_wp, 3.72932356_wp, 3.44677275_wp, &
        & 7.97762753_wp, 7.07623947_wp, 6.60844053_wp, 6.28791364_wp, 6.07728703_wp, 5.54643096_wp, &
        & 5.80491167_wp, 5.58415602_wp, 5.41374528_wp, 5.28497229_wp, 5.22592821_wp, 5.09817141_wp, &
        & 6.12149689_wp, 5.54083734_wp, 5.06696878_wp, 4.87005108_wp, 4.59089647_wp, 4.31176304_wp, &
        & 9.55461698_wp, 8.67396077_wp, 7.97210197_wp, 7.43439917_wp, 6.58711862_wp, 6.19536215_wp, &
        & 6.01517290_wp, 5.81623410_wp, 5.65710424_wp, 5.52640661_wp, 5.44263305_wp, 5.58285373_wp, &
        & 7.02081898_wp, 6.46815523_wp, 5.98089120_wp, 5.81686657_wp, 5.53321815_wp, 5.25477007_wp, &
        & 11.02204549_wp, 10.15679528_wp, 9.35167836_wp, 9.06926079_wp, 8.97241155_wp, 8.90092807_wp, &
        & 8.85984840_wp, 8.81736827_wp, 8.79317710_wp, 7.89969626_wp, 8.80588454_wp, 8.42439218_wp, &
        & 8.54289262_wp, 8.47583370_wp, 8.45090888_wp, 8.47339339_wp, 7.83525634_wp, 8.20702843_wp, &
        & 7.70559063_wp, 7.32755997_wp, 7.03887381_wp, 6.68978720_wp, 6.05450052_wp, 5.88752022_wp, &
        & 5.70661499_wp, 5.78450695_wp, 7.79780729_wp, 7.26443867_wp, 6.78151984_wp, 6.67883169_wp, &
        & 6.39024318_wp, 6.09527958_wp, 11.79156076_wp, 11.10997644_wp, 9.51377795_wp, 8.67197068_wp, &
        & 8.77140725_wp, 8.65402716_wp, 8.53923501_wp, 8.85024712_wp ]

    contains

    !> @brief Evaluate the numerically stable Tang--Toennies damping function of integer order.
    pure function tang_toennies_damping( order, argument ) result(damping)
        integer, intent(in)  :: order
        real(wp), intent(in) :: argument
        real(wp)             :: damping
        real(wp)             :: series_term, series_sum
        integer              :: index
        if( order < 0 .or. argument < 0.0_wp )then
            damping = 0.0_wp
            return
        endif
        if( argument > 1.0_wp )then
            series_term = 1.0_wp; series_sum = 1.0_wp
            do index=1, order
                series_term = series_term*argument/real(index, wp)
                series_sum  = series_sum+series_term
            enddo
            damping = 1.0_wp-exp(-argument)*series_sum
        else
            series_term = 1.0_wp
            do index=1, order+1
                series_term = series_term*argument/real(index, wp)
            enddo
            series_sum = series_term
            do index=order+2, order+1000
                series_term = series_term*argument/real(index, wp)
                series_sum  = series_sum+series_term
                if( abs(series_term) <= epsilon(1.0_wp)*abs(series_sum) ) exit
            enddo
            damping = exp(-argument)*series_sum
        endif
    end function tang_toennies_damping

    !> @brief Evaluate the released N2, NO, or O2 D3(BJ) diatomic correction.
    subroutine evaluate_legacy_d3_diatomic( distance, energy, derivative, derivative_order, molecule_id )
        real(wp), intent(in)  :: distance(1)
        real(wp), intent(out) :: energy, derivative(1)
        integer, intent(in)   :: derivative_order, molecule_id
        integer               :: atomic_numbers(2)
        real(wp)              :: c6_coefficients(1)
        derivative = 0.0_wp
        select case( molecule_id )
        case( 1, 2 )
            atomic_numbers  = [7, 7]
            c6_coefficients = 19.7_wp
        case( 3 )
            atomic_numbers  = [7, 8]
            c6_coefficients = 16.7_wp
        case( 4 )
            atomic_numbers  = [8, 8]
            c6_coefficients = 12.8_wp
        case default
            energy = 0.0_wp
            return
        end select
        if( derivative_order == 1 )then
            call evaluate_d3_bj_dispersion(distance, atomic_numbers, c6_coefficients, 1.0_wp, 2.0_wp, 0.5299_wp, &
                & 2.20_wp, energy, derivative)
        else
            call evaluate_d3_bj_dispersion(distance, atomic_numbers, c6_coefficients, 1.0_wp, 2.0_wp, 0.5299_wp, &
                & 2.20_wp, energy)
        endif
    end subroutine evaluate_legacy_d3_diatomic

    !> @brief Evaluate the released three-pair O3 D3(BJ) correction.
    subroutine evaluate_legacy_d3_o3( distances, energy, derivatives, derivative_order )
        real(wp), intent(in)  :: distances(3)
        real(wp), intent(out) :: energy, derivatives(3)
        integer, intent(in)   :: derivative_order
        integer, parameter    :: atomic_numbers(3)  = [8, 8, 8]
        real(wp), parameter   :: c6_coefficients(3) = 12.8_wp
        derivatives                                 = 0.0_wp
        if( derivative_order == 1 )then
            call evaluate_d3_bj_dispersion(distances, atomic_numbers, c6_coefficients, 1.0_wp, 2.0_wp, 0.5299_wp, &
                & 2.20_wp, energy, derivatives, .true.)
        else
            call evaluate_d3_bj_dispersion(distances, atomic_numbers, c6_coefficients, 1.0_wp, 2.0_wp, 0.5299_wp, &
                & 2.20_wp, energy, legacy_default_real_table=.true.)
        endif
    end subroutine evaluate_legacy_d3_o3

    !> @brief Evaluate the released six-pair O4 D3(BJ) correction.
    subroutine evaluate_legacy_d3_o4( distances, energy, derivatives, derivative_order )
        real(wp), intent(in)  :: distances(6)
        real(wp), intent(out) :: energy, derivatives(6)
        integer, intent(in)   :: derivative_order
        integer, parameter    :: atomic_numbers(4)  = [8, 8, 8, 8]
        real(wp), parameter   :: c6_coefficients(6) = 12.8_wp
        derivatives                                 = 0.0_wp
        if( derivative_order == 1 )then
            call evaluate_d3_bj_dispersion(distances, atomic_numbers, c6_coefficients, 1.0_wp, 2.0_wp, 0.5299_wp, &
                & 2.20_wp, energy, derivatives, .true.)
        else
            call evaluate_d3_bj_dispersion(distances, atomic_numbers, c6_coefficients, 1.0_wp, 2.0_wp, 0.5299_wp, &
                & 2.20_wp, energy, legacy_default_real_table=.true.)
        endif
    end subroutine evaluate_legacy_d3_o4

    !> @brief Evaluate pairwise D3(BJ) energy and distance derivatives in kcal/mol and angstrom.
    subroutine evaluate_d3_bj_dispersion( distances, atomic_numbers, c6_coefficients, s6, s8, a1, a2, energy, &
            & derivatives, legacy_default_real_table )
        real(wp), intent(in)            :: distances(:), c6_coefficients(:), s6, s8, a1, a2
        integer, intent(in)             :: atomic_numbers(:)
        real(wp), intent(out)           :: energy
        real(wp), intent(out), optional :: derivatives(size(distances))
        logical, intent(in), optional   :: legacy_default_real_table
        integer                         :: atom_1, atom_2, pair_index
        real(wp)                        :: c6, c8, damping_radius, distance_bohr, e6, e8, factor_1, factor_2
        pair_index = 0
        energy     = 0.0_wp
        if( present(derivatives) ) derivatives=0.0_wp
        do atom_1=1, size(atomic_numbers)-1
            do atom_2=atom_1+1, size(atomic_numbers)
                pair_index    = pair_index+1
                distance_bohr = distances(pair_index)/d3_bohr_to_angstrom
                c6            = c6_coefficients(pair_index)
                factor_1      = r2r4(atomic_numbers(atom_1))
                factor_2      = r2r4(atomic_numbers(atom_2))
                if( present(legacy_default_real_table) )then
                    if( legacy_default_real_table )then
                        factor_1 = real(real(factor_1, real32), wp)
                        factor_2 = real(real(factor_2, real32), wp)
                    endif
                endif
                c8             = 3.0_wp*c6*factor_1*factor_2
                damping_radius = a1*sqrt(c8/c6)+a2
                e6             = c6/(distance_bohr**6+damping_radius**6)
                e8             = c8/(distance_bohr**8+damping_radius**8)
                energy         = energy+(-s6*e6-s8*e8)*hartree_to_kilocalorie_per_mole
                if( present(derivatives) )then
                    derivatives(pair_index) = (-s6*c6*(-6.0_wp*distance_bohr**5)/ (distance_bohr**6+damping_radius** &
                        & 6)**2 - s8*c8*(-8.0_wp*distance_bohr**7)/ (distance_bohr**8+damping_radius**8)**2)* &
                        & hartree_to_kilocalorie_per_mole/ d3_bohr_to_angstrom
                endif
            enddo
        enddo
    end subroutine evaluate_d3_bj_dispersion
end module chempotpy_dispersion
