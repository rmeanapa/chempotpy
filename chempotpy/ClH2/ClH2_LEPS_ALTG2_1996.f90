! System:          ClH2
! Functional form: eLEPS (extended London-Eyring-Polanyi-Sato)
!                  plus E3C (a three-center term)
!                  and modifications to improve the saddle point
!                  bend potential; based on GQQ surface
! Interface:       3-1V
! Number of derivatives: 0
! Common name:     G3
! Reference:       Allison, T. C.; Lynch, G. C.; Truhlar, D. G.;
!                  Gordon, M. S.,  J. Phys. Chem., Vol. 100, 13575 (1996
! Note:            This is the modified clh2 potential code.
!                  Vector version, includes no derivatives.
! The coordinates are defined as follows:
!           R(1) = R(H-Cl)
!           R(2) = R(H-H')
!           R(3) = R(H'-Cl)
! The zero of energy for this potential occurs when the Cl atom is
! "infinitely" far from the the H2 diatom, and the H2 diatom is at
! its equilibrium seperation.

!> @brief Implement ClH2_LEPS_ALTG2_1996, a potential energy surface for ClH2.
module chempotpy_clh2_leps_altg2_1996_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the ClH2_LEPS_ALTG2_1996 potential energy surface for ClH2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_legacy_units,       only: legacy_hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: legacy_triatomic_distances_to_bohr
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: v
        integer                       :: istate
        double precision              :: r(1,3), e(1)
        !initialize
        v = 0.d0; g = 0.d0; d = 0.d0
        ! input cartesian is ClHH
        call legacy_triatomic_distances_to_bohr(x, r(1,:))
        call pot(r, e, 1, 0)
        e = e*legacy_hartree_to_electronvolt
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = e(1)
            enddo
        else
            write (*, *) 'Only energy is available'
        endif
    end subroutine evaluate
end module chempotpy_clh2_leps_altg2_1996_pes

!> @brief Expose ClH2_LEPS_ALTG2_1996 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_clh2_leps_altg2_1996_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Own the stateless ClH2 ALTG2 numerical kernel.
module chempotpy_clh2_leps_altg2_1996_pot
    implicit none
    private
    public :: dispatch
    double precision, parameter :: energy_per_hartree      = 627.5095d0
    double precision, parameter :: bohr_per_angstrom       = 0.529177249d0
    double precision, parameter :: dissociation_energy(3)  = [106.447d0, 109.458d0, 106.447d0]/energy_per_hartree
    double precision, parameter :: equilibrium_distance(3) = [1.2732d0, 0.74127d0, 1.2732d0]/bohr_per_angstrom
    double precision, parameter :: morse_exponent(3)       = [1.8674d0, 1.9413d0, 1.8674d0]*bohr_per_angstrom
    double precision, parameter :: sato_parameter(3)       = [0.1835d0, 0.167d0, 0.1835d0]
    double precision, parameter :: anti_morse_weight(3)    = 0.5d0*dissociation_energy* &
        & ((1.d0-sato_parameter)/(1.d0+sato_parameter))
    double precision, parameter :: three_center_energy = 0.0758016022d0
    double precision, parameter :: symmetric_decay     = 0.0008627355d0
    double precision, parameter :: asymmetric_decay    = 0.2981969860d0
    double precision, parameter :: bond_decay          = 0.1439106543d0
    double precision, parameter :: bend_quadratic      = 0.6940323070d0
    double precision, parameter :: bend_quartic        = 1.6963092005d0
    double precision, parameter :: first_modifier      = 3.1053877618397071175758280546d0
    double precision, parameter :: second_modifier     = -1.8942608491155350536449828878d0

    contains

    !> @brief Evaluate a batch of ClH2 ALTG2 internal-coordinate geometries.
    subroutine dispatch( rr, e, nt, nsurf )
        implicit none
        integer, intent(in)           :: nt
        integer, intent(in)           :: nsurf
        double precision, intent(in)  :: rr(nt,3)
        double precision, intent(out) :: e(nt)
        double precision              :: cosine_theta, first_factor, first_leps, first_morse, first_y
        double precision              :: fourth_factor, modifier, radical, second_factor, second_leps
        double precision              :: second_morse, second_y, sine_theta, third_factor, third_leps
        double precision              :: third_morse, third_y, three_center, triangle_asymmetry
        double precision              :: distances(3)
        integer                       :: geometry_index
        do geometry_index=1, nt
            distances     = rr(geometry_index,:)
            first_factor  = exp(-symmetric_decay*(distances(1)+distances(3))**4)
            second_factor = exp(-asymmetric_decay*(distances(1)-distances(3))**2)
            third_factor  = exp(-bond_decay*(distances(1)+distances(3)-distances(2))**2)
            cosine_theta  = (distances(1)**2+distances(3)**2-distances(2)**2)/(2.d0*distances(1)*distances(3))
            if( cosine_theta .lt. -1.d0 )then
                cosine_theta = -1.d0
            elseif( cosine_theta .gt. 1.d0 )then
                cosine_theta = 1.d0
            endif
            sine_theta         = sin(acos(cosine_theta))
            fourth_factor      = 1.d0+sine_theta**2*(bend_quadratic+bend_quartic*sine_theta**2)
            three_center       = fourth_factor*three_center_energy*first_factor*second_factor*third_factor
            triangle_asymmetry = (distances(1)-distances(3))/distances(2)
            modifier           = (0.5d0*(1.d0+cosine_theta))**6
            modifier = 1.d0+first_modifier*modifier*(1.d0-triangle_asymmetry**2) + second_modifier*modifier*(1.d0- &
                & triangle_asymmetry**4)
            first_y = exp(-morse_exponent(1)*(distances(1)-equilibrium_distance(1)))
            first_morse = 0.5d0*(dissociation_energy(1)*(first_y-2.d0) + anti_morse_weight(1)*(modifier*first_y+ &
                & 2.d0))*first_y
            first_leps = 0.5d0*(dissociation_energy(1)*(first_y-2.d0) - anti_morse_weight(1)*(modifier*first_y+ &
                & 2.d0))*first_y
            second_y = exp(-morse_exponent(2)*(distances(2)-equilibrium_distance(2)))
            second_morse = 0.5d0*(dissociation_energy(2)*(second_y-2.d0) + anti_morse_weight(2)*(second_y+2.d0))* &
                & second_y
            second_leps = 0.5d0*(dissociation_energy(2)*(second_y-2.d0) - anti_morse_weight(2)*(second_y+2.d0))* &
                & second_y
            third_y = exp(-morse_exponent(3)*(distances(3)-equilibrium_distance(3)))
            third_morse = 0.5d0*(dissociation_energy(3)*(third_y-2.d0) + anti_morse_weight(3)*(modifier*third_y+ &
                & 2.d0))*third_y
            third_leps = 0.5d0*(dissociation_energy(3)*(third_y-2.d0) - anti_morse_weight(3)*(modifier*third_y+ &
                & 2.d0))*third_y
            radical = sqrt(first_leps*(first_leps-second_leps) + second_leps*(second_leps-third_leps)+third_leps* &
                & (third_leps-first_leps))
            e(geometry_index) = dissociation_energy(2)+first_morse+second_morse+third_morse-radical+three_center
        enddo
    end subroutine dispatch
end module chempotpy_clh2_leps_altg2_1996_pot

!> @brief Compatibility interface for pot.
subroutine pot( rr, e, nt, nsurf )
    use chempotpy_clh2_leps_altg2_1996_pot, only: dispatch
    implicit none
    double precision :: rr(nt,3)
    double precision :: e(nt)
    integer          :: nt
    integer          :: nsurf
    call dispatch(rr, e, nt, nsurf)
end subroutine pot

!> @brief Evaluate POTINFO for the ClH2_LEPS_ALTG2_1996 potential.
