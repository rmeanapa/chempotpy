!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_regea_pes
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
    !> @brief Evaluate pot for the Aln_reGEA potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        !   System:                        Aln
        !   Functional form:               Embedded atom
        !   Common name:                   reGEA
        !   Number of derivatives:         0
        !   Number of bodies:              variable
        !   Number of electronic surfaces: 1
        !   Interface:                     HO-MM-0
        !
        !   Notes:  Many-body aluminum potential energy function.  The parameters
        !           for this PEF have been parameterized against a data set of
        !           aluminum cluster energies and bulk data.
        !
        !  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
        !              and D. G. Truhlar, "Analytic Potentials Energy Functions
        !              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
        !
        !  Units:
        !       energies    - hartrees
        !       coordinates - bohrs
        !
        !  v       --- energy of the system (output)
        !  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
        !              of the system (input)
        !  natom   --- number of atoms in the system (input)
        !  maxatom --- maximal number of atoms (input)
        implicit none
        double precision :: a, b, cc, dx, dy, dz, rr, u, ud, udd, uu, v, v0, v1, x, y, z
        integer          :: i, j, maxatom, natom
        dimension x(maxatom), y(maxatom), z(maxatom)
        !
        !  calculate the potential energy
        ! new parameters
        a  = 1.82725940400586229
        b  = 0.701319003419638487
        cc = 0.260371275036638983
        v0 = 90.0439667806546140
        ! Original
        !      a=1.589647708
        !      b=0.866262749
        !      cc=0.4865568924
        !      v0=44.8337922
        v1 = 0d0
        do i=1, natom
            ud = 0d0
            do j=1, natom
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = dsqrt(dx*dx+dy*dy+dz*dz)
                if( j .eq. i ) cycle
                udd = 0d0
                udd = dexp(-b*rr)
                ud  = ud+udd
                41 continue
            5 enddo
            v1 = v1+dsqrt(ud)
        4 enddo
        u = 0d0
        do i=1, natom-1
            do j=i+1, natom
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = dsqrt(dx*dx+dy*dy+dz*dz)
                uu = 0d0
                uu = dexp(-a*rr)
                u  = u+uu
            3 enddo
        2 enddo
        v = u*v0-v1*cc
    end subroutine pot
    !
end module chempotpy_aln_regea_pes

!> @brief Preserve the historical Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_regea_pes, only: evaluate
    use chempotpy_kinds,         only: real64
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
