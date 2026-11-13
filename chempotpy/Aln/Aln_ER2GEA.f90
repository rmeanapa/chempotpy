!   System:                        Aln
!   Functional form:               Extended Rydberg plus generalized embedded-atom
!   Common name:                   ER2+GEA
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  The parameters
!           for this PEF have been parameterized against a data set of
!           aluminum cluster energies and bulk data.
!           This PEF has a cutoff at 6.5 angstroms.
!  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!              and D. G. Truhlar, "Analytic Potentials Energy Functions
!              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!                    of the system
!  natom   --- number of atoms in the system
!  maxatom --- maximal number of atoms

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_er2gea_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate an aluminum-cluster potential through the Python contract.
    subroutine evaluate( coordinates_angstrom, atom_count, derivative_order, &
            energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
        use chempotpy_units, only: coordinates_to_bohr, &
            energy_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter    :: state_count = 1
        real(wp), intent(in)  :: coordinates_angstrom(10000,3)
        integer, intent(in)   :: atom_count, derivative_order
        real(wp), intent(out) :: energy_ev(state_count)
        real(wp), intent(out) :: gradient_ev_per_angstrom( &
            state_count,atom_count,3 )
        real(wp), intent(out) :: coupling_per_angstrom( &
            state_count,state_count,atom_count,3 )
        real(wp) :: coordinates_bohr(atom_count,3), energy_hartree
        energy_ev                = 0.0e0_wp
        gradient_ev_per_angstrom = 0.0e0_wp
        coupling_per_angstrom    = 0.0e0_wp
        call coordinates_to_bohr( &
            coordinates_angstrom(1:atom_count,:), coordinates_bohr)
        call pot(coordinates_bohr(:,1), coordinates_bohr(:,2), &
            coordinates_bohr(:,3), energy_hartree, atom_count, atom_count)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate pot for the Aln_ER2GEA potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a1, a2, a3, au2, au23, b, bb, bu2, c(4), co, co3, de, dx, dy, dz, ee, en, ff, fu2(maxatom,maxatom), &
            & poly, r(maxatom, maxatom), re, rho, rr, u2, ud, ud2, v, x(maxatom), y(maxatom), z(maxatom)
        integer             :: i, j, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        de   = 1.71013678553981441/autoev
        re   = 5.08182706399609163
        a1   = 1.24074255007327805
        a2   = 0.551880801172447422
        a3   = 0.129970688812896917
        au2  = 0.143243771372740580
        bu2  = 6.50000000000000000/autoang
        bb   = 4.05080605764533441/autoev
        b    = 1.12872496336101613*autoang
        en   = 0.558622374206155348
        au23 = 0.263351245725451877
        ! potential, includes cutoff
        v  = 0.e0_wp
        u2 = 0.e0_wp
        ff = 0.e0_wp
        do i=1, natom
            ud  = 0.e0_wp
            ud2 = 0.e0_wp
            do j=1, natom
                if( i .eq. j ) cycle
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = sqrt(dx*dx+dy*dy+dz*dz)
                if( rr .ge. bu2 ) cycle
                rho  = rr-re
                poly = 1.e0_wp+a1*rho+a2*rho**2+a3*rho**3
                ee   = exp(-a1*rho)
                co   = 0.e0_wp
                co3  = 0.e0_wp
                if( rr .lt. bu2 ) co=exp(au2*(1.e0_wp-1.e0_wp/(1.e0_wp-rr/bu2)))
                if( rr .lt. bu2 ) co3=exp(au23*(1.e0_wp-1.e0_wp/(1.e0_wp-rr/bu2)))
                u2  = u2-0.5e0_wp*de*ee*poly*co
                ud  = ud+exp(-rr*b)*co3
                ud2 = ud2+(co3*exp(-rr*b))**en
            enddo
            ff = ff-bb*ud**en+bb*ud2
        enddo
        v = u2+ff
    end subroutine pot
end module chempotpy_aln_er2gea_pes

!> @brief Expose Aln_ER2GEA through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_er2gea_pes, only: evaluate
    use chempotpy_kinds,          only: wp
    implicit none
    integer, parameter    :: state_count = 1
    real(wp), intent(in)  :: coordinates_angstrom(10000,3)
    integer, intent(in)   :: atom_count, derivative_order
    real(wp), intent(out) :: energy_ev(state_count)
    real(wp), intent(out) :: gradient_ev_per_angstrom(state_count,atom_count,3)
    real(wp), intent(out) :: coupling_per_angstrom(state_count,state_count,atom_count,3)
    call evaluate(coordinates_angstrom, atom_count, derivative_order, energy_ev, gradient_ev_per_angstrom, &
        & coupling_per_angstrom)
end subroutine pes
