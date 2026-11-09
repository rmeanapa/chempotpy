module chempotpy_clh2_leps_stsbltg1_1989_pes
    implicit none
    private
    public :: evaluate
    double precision, parameter :: energy_per_hartree         = 627.5095d0
    double precision, parameter :: bohr_per_angstrom          = 0.52917706d0
    double precision, parameter :: root_two                   = 1.41421356d0
    double precision, parameter :: dissociation_energy(3)     = [106.447d0, 109.458d0, 106.447d0]/energy_per_hartree
    double precision, parameter :: equilibrium_distance(3)    = [1.2732d0, 0.74127d0, 1.2732d0]/bohr_per_angstrom
    double precision, parameter :: morse_exponent(3)          = [1.8674d0, 1.94130d0, 1.8674d0]*bohr_per_angstrom
    double precision, parameter :: sato_parameter(3)          = [0.187d0, 0.167d0, 0.187d0]
    double precision, parameter :: one_plus_sato(3)           = 1.d0+sato_parameter
    double precision, parameter :: one_plus_three_sato(3)     = 1.d0+3.d0*sato_parameter
    double precision, parameter :: two_one_plus_three_sato(3) = 2.d0*one_plus_three_sato
    double precision, parameter :: sato_plus_three(3)         = sato_parameter+3.d0
    double precision, parameter :: two_sato_plus_three(3)     = 2.d0*sato_plus_three
    double precision, parameter :: quarter_energy(3)          = dissociation_energy/(4.d0*one_plus_sato)
    double precision, parameter :: derivative_weight(3)       = 2.d0*morse_exponent*quarter_energy
    double precision, parameter :: three_center_energy        = 0.07672d0
    double precision, parameter :: asymmetric_decay           = 0.3d0
    double precision, parameter :: symmetric_decay            = 0.000863d0
    double precision, parameter :: bond_decay                 = 0.15d0
    double precision, parameter :: bend_quadratic             = 0.0d0
    double precision, parameter :: bend_quartic               = 1.725d0
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_legacy_units,       only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: legacy_triatomic_distances_to_bohr
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: dedr(3), energy, r(3), v
        integer                       :: iatom, idir, istate, j
        double precision              :: tx(9), r2(3), dedx(9), drdx(3,9)
        !initialize
        v = 0.d0; g = 0.d0; d = 0.d0
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is ClHH
        call legacy_triatomic_distances_to_bohr(x, r)
        dedr = 0.d0
        call dispatch(r, igrad, energy, dedr)
        energy = energy*legacy_hartree_to_electronvolt
        dedr   = dedr*51.422067
        r2     = r*legacy_bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(dedr, drdx, dedx)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = energy
            enddo
        else if( igrad == 1 )then
            do istate=1, nstates
                p(istate) = energy
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    j               = (iatom-1)*3+idir
                    g(1,iatom,idir) = dedx(j)
                enddo
            enddo
        else if( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
    end subroutine evaluate
    !> @brief Evaluate EvdRdX for the ClH2_LEPS_STSBLTG1_1989 potential.
    subroutine evdrdx( x, r, drdx )
        use chempotpy_geometry, only: ordered_pair_distance_jacobian
        implicit none
        integer, parameter            :: first_atoms(3)  = [1, 2, 1]
        integer, parameter            :: second_atoms(3) = [2, 3, 3]
        double precision, intent(in)  :: x(9), r(3)
        double precision, intent(out) :: drdx(3,9)
        call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
    end subroutine evdrdx
    !> @brief Evaluate dispatch.
    subroutine dispatch( r, nder, energy, dedr )
        !
        !   System:          ClH2
        !   Functional form: eLEPS (extended-London-Eyring-Polanyi-Sato)
        !                    plus E3C (a three-center term).
        !   Common name: GQ
        !   Reference:       D. W. Schwenke, S. C. Tucker, R. Steckler, F. B. Br
        !                    G. C. Lynch, D. G. Truhlar, and B. C. Garrett
        !                    J. Chem. Phys. 90, 3110-3120 (1989)
        !
        !   Calling Sequence:
        !               must be called once before any calls to POT
        !      POT    - driver for the evaluation of the energy and the derivati
        !               of the energy with respect to the coordinates for a give
        !               geometric configuration
        !
        !   Units:
        !      energies    - hartrees
        !      coordinates - bohr
        !      derivatives - hartrees/bohr
        !
        !   Surfaces:
        !      ground electronic state
        !
        !   Zero of energy:
        !      The classical potential energy is set equal to zero for the Cl
        !      infinitely far from the H2 diatomic and R(H2) set equal to the
        !      H2 equilibrium diatomic value.
        !
        !   Parameters:
        !
        !   Coordinates:
        !      Internal, Definition: R(1) = R(Cl-H)
        !                            R(2) = R(H-H)
        !                            R(3) = R(Cl-H)
        !
        !   Common Blocks (used between the calling program and this potential):
        !      /PT1CM/ R(3), ENERGY, DEDR(3)
        !        passes the coordinates, ground state electronic energy, and
        !        derivatives of the ground electronic state energy with respect
        !        to the coordinates.
        !      /PT2CM/ NSURF, NDER, NDUM(21)
        !        passes the control flags where
        !        NDER  = 0 => no derivatives are computed
        !              = 1 => derivatives of the energy for the ground electroni
        !                     state with respect to the coordinates are computed
        !        NDUM  - not used
        !      /PT4CM/ IPTPRT, IDUM(19)
        !        IPTPRT passes the FORTRAN unit number used for potential output
        !        IDUM   not used
        !      /PT5CM/ EASYAB, EASYBC, EASYAC
        !        passes the energy in the three asymptotic valleys for an A + BC
        !        The energy in the AB valley, EASYAB, is equal to the energy of
        !        C atom "infinitely" far from the AB diatomic and R(AB) set equa
        !        equilibrium_distance(AB), the equilibrium bond length for the AB diatomic.
        !        Similarly, the terms EASYBC and EASYAC represent the energies
        !        in the BC and the AC valleys, respectively.
        !
        !   Default Parameter Values:
        !      Variable      Default value
        !      NSURF            0
        !      NDER             1
        !      IPTPRT           6
        !
        implicit none
        double precision, intent(in)  :: r(3)
        integer, intent(in)           :: nder
        double precision, intent(out) :: energy, dedr(3)
        double precision              :: abcval, ang1, ang11, ang12, ang13, ang14, ang2, ang21, ang22
        double precision              :: ang3, ang31, ang32, ang33, ang34, ctheta, d3cdr1, d3cdr2, d3cdr3
        double precision              :: dist, e3c, exch, f1val, f2val, f3val, f4val
        double precision              :: r1p2, r1p3, r1p4, r1p5, r3p2, r3p3, r3p4, r3p5
        double precision              :: rad, rco, rdif, rsum3, stheta, sval, theta, trm1, trm2, trm3, tval, x
        integer                       :: i
        dimension x(3), exch(3)
        double precision :: coul(3)
        !
        !WRITE (IPTPRT, 600) D, RE, BETA, Z
        !WRITE (IPTPRT, 610) three_center_energy, symmetric_decay, bond_decay, asymmetric_decay, bend_quadratic, bend_quartic
        !
        !   Check the value of NDER
        !
        !   Compute the eLEPS part of the potential
        !
        do i = 1, 3
            x(i)    = exp(-morse_exponent(i)*(r(i)-equilibrium_distance(i)))
            coul(i) = quarter_energy(i)*(sato_plus_three(i)*x(i)-two_one_plus_three_sato(i))*x(i)
            exch(i) = quarter_energy(i)*(one_plus_three_sato(i)*x(i)-two_sato_plus_three(i))*x(i)
        enddo
        rad    = sqrt((exch(1)-exch(2))**2+(exch(2)-exch(3))**2+ (exch(3)-exch(1))**2)
        energy = coul(1) + coul(2) + coul(3) - (rad/root_two) + dissociation_energy(2)
        !
        !   Compute the derivatives of the eLEPS energy w.r.t. the bond lengths
        !
        if( nder .eq. 1 )then
            sval = exch(1) + exch(2) + exch(3)
            do i = 1, 3
                dedr(i) = 0.d0
                if( x(i) .lt. 1.d-30 ) cycle
                tval = (3.0d0*exch(i)-sval)/root_two*(one_plus_three_sato(i)*x(i)-sato_plus_three(i))
                if( abs(rad) .lt. 1.d-32 .and.abs(tval) .gt. 1.d-12 )then
                else
                    tval = tval/rad
                endif
                dedr(i) = derivative_weight(i)*x(i)*(tval-sato_plus_three(i)*x(i)+one_plus_three_sato(i))
            enddo
        endif
        !
        !   Compute the three-center term
        !
        !   F1VAL is equal to f1 in the literature reference,
        !   similarly F2VAL = f2, and F3VAL = f3
        f1val = exp(-asymmetric_decay*(r(1)-r(3))**2)
        f2val = exp(-symmetric_decay*(r(1)+r(3))**4)
        f3val = exp(-bond_decay*(r(1)+r(3)-r(2))**2)
        e3c   = three_center_energy * f1val * f2val * f3val
        !
        !   Compute the angular dependent term, f4, which is called F4VAL
        dist   = r(1)**2 + r(3)**2 - r(2)**2
        ctheta = dist/(2.0d0*r(1)*r(3))
        if( ctheta .lt. -1.0d0 )then
            theta = acos(-1.0d0)
        else if( ctheta .ge. 1.0d0 )then
            theta = 0.0d0
        else
            theta = acos(ctheta)
        endif
        stheta = sin(theta)
        !
        f4val  = (bend_quadratic*(stheta**2))+(bend_quartic * (stheta**4)) + 1.0d0
        e3c    = e3c * f4val
        energy = e3c + energy
        !
        !   Compute the derivatives of the three-center term
        !
        if( nder .eq. 1 )then
            !   Recall DIST = R(1)**2 + R(3)**2 - R(2)**2
            rsum3   = (r(1) + r(3))**3
            rdif    = r(1) - r(3)
            rco     = r(1) + r(3) - r(2)
            abcval  = e3c/f4val
            r1p5    = r(1)**5
            r1p4    = r(1)**4
            r1p3    = r(1)**3
            r1p2    = r(1)**2
            r3p5    = r(3)**5
            r3p4    = r(3)**4
            r3p3    = r(3)**3
            r3p2    = r(3)**2
            trm1    = 4.0d0*symmetric_decay*rsum3
            trm2    = 2.0d0*bond_decay*rco
            trm3    = 2.0d0*asymmetric_decay*rdif
            ang11   = -(dist**4)/(4.0d0*r1p5*r3p4)
            ang12   = (dist**3)/(2.0d0*r1p3*r3p4)
            ang13   = (dist**2)/(r1p3*r3p2)
            ang14   = (-2.0d0*dist)/(r(1)*r3p2)
            ang1    = bend_quartic * (ang11 + ang12 + ang13 + ang14) + bend_quadratic * ((ang13/2.0d0) + (ang14/ 2.0d0))
            ang31   = -(dist**4)/(4.0d0*r1p4*r3p5)
            ang32   = (dist**3)/(2.0d0*r1p4*r3p3)
            ang33   = (dist**2)/(r1p2*r3p3)
            ang34   = (-2.0d0*dist)/(r1p2*r(3))
            ang3    = bend_quartic * (ang31 + ang32 + ang33 + ang34) + bend_quadratic * (ang33/2.0d0 + ang34/2.0d0)
            ang21   = -r(2)*(dist**3)/(2.0d0*r1p4*r3p4)
            ang22   = r(2)*2.0d0*dist/(r1p2*r3p2)
            ang2    = bend_quartic * (ang21 + ang22) + bend_quadratic * (ang22/2.0d0)
            d3cdr1  = e3c * (-trm1 - trm2 - trm3) + abcval * ang1
            d3cdr2  = e3c * trm2 + abcval * ang2
            d3cdr3  = e3c * (-trm1 - trm2 + trm3) + abcval * ang3
            dedr(1) = dedr(1) + d3cdr1
            dedr(2) = dedr(2) + d3cdr2
            dedr(3) = dedr(3) + d3cdr3
        endif
        !
        !
    end subroutine dispatch
end module chempotpy_clh2_leps_stsbltg1_1989_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_clh2_leps_stsbltg1_1989_pes, only: evaluate
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
