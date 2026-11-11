!   System:                        Aln
!   Functional form:               Lennard-Jones plus Axilrod-Teller
!   Common name:                   deSPH/LJAT
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  For a recent
!           discussion of this potential, see:
!           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!           and D. G. Truhlar, "Analytic Potentials Energy Functions
!           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
!  References: P. de Sainte Claire, G. H. Peslherbe, and W. L. Hase
!              J. Phys. Chem., 99, 8147 (1995).
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximal number of atoms (input)

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_desph_ljat_pes
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

    !> @brief Evaluate pot for the Aln_deSPH_LJAT potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        implicit none
        double precision :: a, aa, b, c, co1, co2, co3, dx, dy, dz, rij, rim, rjm, rr, sig, u2, u3, v, x(maxatom), y(maxatom), &
            & z(maxatom)
        integer                     :: i, j, m, maxatom, natom
        double precision, parameter :: autoev  = 27.2113961d0
        double precision, parameter :: autoang = 0.529177249d0
        double precision, parameter :: autokpm = 627.5095d0
        !  calculate the potential energy
        ! in kJ/mol
        !      a = 2975343.77d0 / ((autoang**12) * autokpm)
        !      b = -17765.823d0 / ((autoang**6)  * autokpm)
        !      c = 81286.093d0  / ((autoang**9)  * autokpm)
        ! in eV
        a   = 129023.16d0   / ((autoang**12) * autoev)
        b   = -770.399d0    / ((autoang**6)  * autoev)
        c   = 3524.899741d0 / ((autoang**9)  * autoev)
        sig = 2.3478d0  /   autoang
        aa  = 4.6d0      /                  autoev
        v   = 0.d0
        !  1. calculate energy of two-body interactions L-J
        u2 = 0.d0
        do i=1, natom-1
            do j=i+1, natom
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = dsqrt(dx*dx+dy*dy+dz*dz)
                u2 = u2 + a/rr**12 + b/rr**6
            enddo
        enddo
        !  2. calculate energy of three-body interactions
        u3 = 0.d0
        if( natom /= 2 )then
            do i=1, natom-2
                do j=i+1, natom-1
                    do m=j+1, natom
                        dx  = x(i)-x(j)
                        dy  = y(i)-y(j)
                        dz  = z(i)-z(j)
                        rij = dsqrt(dx*dx+dy*dy+dz*dz)
                        dx  = x(m)-x(j)
                        dy  = y(m)-y(j)
                        dz  = z(m)-z(j)
                        rjm = dsqrt(dx*dx+dy*dy+dz*dz)
                        dx  = x(m)-x(i)
                        dy  = y(m)-y(i)
                        dz  = z(m)-z(i)
                        rim = dsqrt(dx*dx+dy*dy+dz*dz)
                        call scos(i, j, m, x, y, z, rij, rjm, rim, natom, maxatom, co1, co2, co3)
                        u3 = u3+c*(1.d0+3.d0*co1*co2*co3)/((rij*rjm*rim)**3)
                    enddo
                enddo
            enddo
            !  calculate the potential energy
        endif
        v = u2+u3
    end subroutine pot

    !> @brief Evaluate scos for the Aln_deSPH_LJAT potential.
    subroutine scos( i, j, m, x, y, z, rij, rjm, rim, natom, maxatom, &
            & co1, co2, co3 )
        !  The subroutine calculates the cosines of the bound angles
        !  using scalar products of the position vectors of two
        !  atoms originating from the third atom.
        !  Units:
        !  coordinates - bohrs
        !  i,j,m   --- natural indices running from 1 to natom
        !  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
        !              of the system
        !  r       --- two-dimensional array of atomic distances
        !  natom   --- number of atoms in the system
        !  maxatom --- maximal number of atoms
        !  co1     --- cosine of the angle between bonds ij and im
        !  co2     --- cosine of the angle between bonds ji and jm
        !  co3     --- cosine of the angle between bonds mi and mj
        implicit none
        double precision :: co1, co2, co3, rij, rim, rjm, scalar1, scalar2, scalar3, x(maxatom), y(maxatom), z(maxatom)
        integer          :: i, j, m, maxatom, natom
        !  calculate scalar products of the position vectors of two atoms
        !  originating from the third atom
        scalar1 = (x(j)-x(i))*(x(m)-x(i)) +(y(j)-y(i))*(y(m)-y(i)) +(z(j)-z(i))*(z(m)-z(i))
        scalar2 = (x(i)-x(j))*(x(m)-x(j)) +(y(i)-y(j))*(y(m)-y(j)) +(z(i)-z(j))*(z(m)-z(j))
        scalar3 = (x(i)-x(m))*(x(j)-x(m)) +(y(i)-y(m))*(y(j)-y(m)) +(z(i)-z(m))*(z(j)-z(m))
        !  calculate cosine of the angle between bonds ij and im
        co1 = scalar1/(rij*rim)
        !  calculate cosine of the angle between bonds ji and jm
        co2 = scalar2/(rij*rjm)
        !  calculate cosine of the angle between bonds mi and mj
        co3 = scalar3/(rim*rjm)
    end subroutine scos
end module chempotpy_aln_desph_ljat_pes

!> @brief Expose Aln_deSPH_LJAT through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_desph_ljat_pes, only: evaluate
    use chempotpy_kinds,              only: real64
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
