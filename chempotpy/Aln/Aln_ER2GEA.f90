!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_er2gea_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate an aluminum-cluster potential through the Python contract.
    subroutine evaluate( coordinates_angstrom, atom_count, derivative_order, &
            energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
        use chempotpy_legacy_units, only: legacy_coordinates_to_bohr, &
            legacy_energy_to_electronvolt
        use chempotpy_kinds, only: real64
        implicit none
        integer, parameter        :: state_count = 1
        real(real64), intent(in)  :: coordinates_angstrom(10000,3)
        integer, intent(in)       :: atom_count, derivative_order
        real(real64), intent(out) :: energy_ev(state_count)
        real(real64), intent(out) :: gradient_ev_per_angstrom( &
            state_count,atom_count,3 )
        real(real64), intent(out) :: coupling_per_angstrom( &
            state_count,state_count,atom_count,3 )
        real(real64) :: coordinates_bohr(atom_count,3), energy_hartree
        energy_ev                = 0.0d0
        gradient_ev_per_angstrom = 0.0d0
        coupling_per_angstrom    = 0.0d0
        call legacy_coordinates_to_bohr( &
            coordinates_angstrom(1:atom_count,:), coordinates_bohr)
        call pot(coordinates_bohr(:,1), coordinates_bohr(:,2), &
            coordinates_bohr(:,3), energy_hartree, atom_count, atom_count)
        if( derivative_order == 0 )then
            energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate
    !> @brief Evaluate pot for the Aln_ER2GEA potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        !   System:                        Aln
        !   Functional form:               Extended Rydberg plus generalized embedded-atom
        !   Common name:                   ER2+GEA
        !   Number of derivatives:         0
        !   Number of bodies:              variable
        !   Number of electronic surfaces: 1
        !   Interface:                     HO-MM-0
        !
        !   Notes:  Many-body aluminum potential energy function.  The parameters
        !           for this PEF have been parameterized against a data set of
        !           aluminum cluster energies and bulk data.
        !           This PEF has a cutoff at 6.5 angstroms.
        !
        !  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
        !              and D. G. Truhlar, "Analytic Potentials Energy Functions
        !              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
        !
        !  Units:
        !       energies    - hartrees
        !       coordinates - bohrs
        !
        !  v       --- energy of the system
        !  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
        !                    of the system
        !  natom   --- number of atoms in the system
        !  maxatom --- maximal number of atoms
        implicit none
        double precision :: a1, a2, a3, au2, au23, b, bb, bu2, c, co, co3, de, dx, &
            & dy, dz, ee, en, ff, fu2, poly, r, re, rho, rr, u2, ud, ud2, v, x, y, z
        integer :: i, j, maxatom, natom
        dimension x(maxatom), y(maxatom), z(maxatom)
        dimension fu2(maxatom,maxatom), r(maxatom,maxatom), c(4)
        double precision, parameter :: autoev  = 27.2113961d0
        double precision, parameter :: autoang = 0.529177249d0
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
        v  = 0.d0
        u2 = 0.d0
        ff = 0.d0
        do i=1, natom
            ud  = 0.d0
            ud2 = 0.d0
            do j=1, natom
                if( i .eq. j ) cycle
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = dsqrt(dx*dx+dy*dy+dz*dz)
                if( rr .ge. bu2 ) cycle
                rho  = rr-re
                poly = 1.d0+a1*rho+a2*rho**2+a3*rho**3
                ee   = dexp(-a1*rho)
                co   = 0.d0
                co3  = 0.d0
                if( rr .lt. bu2 ) co=dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
                if( rr .lt. bu2 ) co3=dexp(au23*(1.d0-1.d0/(1.d0-rr/bu2)))
                u2  = u2-0.5d0*de*ee*poly*co
                ud  = ud+dexp(-rr*b)*co3
                ud2 = ud2+(co3*dexp(-rr*b))**en
            2 enddo
            ff = ff-bb*ud**en+bb*ud2
        1 enddo
        v = u2+ff
    end subroutine pot
end module chempotpy_aln_er2gea_pes

!> @brief Preserve the historical Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_er2gea_pes, only: evaluate
    use chempotpy_kinds,          only: real64
    implicit none
    integer, parameter        :: state_count = 1
    real(real64), intent(in)  :: coordinates_angstrom(10000,3)
    integer, intent(in)       :: atom_count, derivative_order
    real(real64), intent(out) :: energy_ev(state_count)
    real(real64), intent(out) :: gradient_ev_per_angstrom(state_count,atom_count,3)
    real(real64), intent(out) :: coupling_per_angstrom(state_count,state_count,atom_count,3)
    call evaluate(coordinates_angstrom, atom_count, derivative_order, energy_ev, gradient_ev_per_angstrom, &
        & coupling_per_angstrom)
end subroutine pes
