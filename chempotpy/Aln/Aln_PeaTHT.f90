!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_peatht_pes
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
    !> @brief Evaluate pot for the Aln_PeaTHT potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        !   System:                        Aln
        !   Functional form:               Lennard-Jones plus Axilrod-Teller
        !   Common name:                   PeaTHT
        !   Number of derivatives:         0
        !   Number of bodies:              variable
        !   Number of electronic surfaces: 1
        !   Interface:                     HO-MM-0
        !
        !   Notes:  Many-body aluminum potential energy function.  For a recent
        !           discussion of this potential, see:
        !           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
        !           and D. G. Truhlar, "Analytic Potentials Energy Functions
        !           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
        !
        !  References: E. Pearson, T. Takai, T. Halicioglu, and
        !                 W.A. Tiller, J. Cryst. Growth 70, 33 (1984).
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
        double precision :: aa, cote1, cote2, cote3, dx, dy, dz, em, en, r0, r1, r2, r3, &
            & rel, rr, sig, u, um, un, utr, v, x, y, z, zet
        integer :: i, j, m, maxatom, natom
        dimension x(maxatom), y(maxatom), z(maxatom)
        double precision, parameter :: autoev  = 27.2113961d0
        double precision, parameter :: autoang = 0.529177249d0
        !  calculate the potential energy
        aa  = 4.d0*1.216d0/autoev
        r0  = 2.520d0/autoang
        en  = 6.d0
        em  = 12.d0
        sig = (en/em)**(1.d0/(em-en))*r0
        zet = 2241.d0/(autoev*autoang**9)
        v   = 0d0
        !     calculate energy of two-body interactions
        um = 0.d0
        un = 0.d0
        do i=1, natom-1
            do j=i+1, natom
                dx  = x(i)-x(j)
                dy  = y(i)-y(j)
                dz  = z(i)-z(j)
                rr  = dsqrt(dx*dx+dy*dy+dz*dz)
                rel = sig/rr
                um  = um+rel**em
                un  = un+rel**en
            enddo
        enddo
        u = aa*(um-un)
        !  2. calculate energy of three-body interactions
        utr = 0
        if( natom /= 2 )then
            do i=1, natom-2
                do j=i+1, natom-1
                    do m=j+1, natom
                        dx = x(i)-x(j)
                        dy = y(i)-y(j)
                        dz = z(i)-z(j)
                        r1 = dsqrt(dx*dx+dy*dy+dz*dz)
                        dx = x(m)-x(j)
                        dy = y(m)-y(j)
                        dz = z(m)-z(j)
                        r2 = dsqrt(dx*dx+dy*dy+dz*dz)
                        dx = x(m)-x(i)
                        dy = y(m)-y(i)
                        dz = z(m)-z(i)
                        r3 = dsqrt(dx*dx+dy*dy+dz*dz)
                        call scos(i, j, m, x, y, z, r1, r2, r3, natom, maxatom, cote1, cote2, cote3)
                        utr = utr+zet*(1+3d0*cote1*cote2*cote3)/((r1*r2*r3)**3d0)
                    enddo
                enddo
            enddo
            !  calculate the potential energy
        endif
        v = u+utr
    end subroutine pot
    !
    !
    !
    !> @brief Evaluate scos for the Aln_PeaTHT potential.
    subroutine scos( i, j, m, x, y, z, r1, r2, r3, natom, maxatom, &
            & cote1, cote2, cote3 )
        !  The subroutine calculates the cosines of the bound angles
        !  using scalar products of the position vectors of two
        !  atoms originating from the third atom.
        !
        !  Units:
        !  coordinates - bohrs
        !
        !  i,j,m   --- natural indices running from 1 to natom
        !  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
        !              of the system
        !  r       --- two-dimensional array of atomic distances
        !  natom   --- number of atoms in the system
        !  maxatom --- maximal number of atoms
        !  cote1   --- cosine of the angle between bonds ij and im
        !  cote2   --- cosine of the angle between bonds ji and jm
        !  cote3   --- cosine of the angle between bonds mi and mj
        !
        implicit none
        double precision :: cote1, cote2, cote3, r1, r2, r3, scalar1, scalar2, scalar3, x, y, z
        integer          :: i, j, m, maxatom, natom
        dimension x(maxatom), y(maxatom), z(maxatom)
        !  calculate scalar products of the position vectors of two atoms
        !  originating from the third atom
        scalar1 = (x(j)-x(i))*(x(m)-x(i)) +(y(j)-y(i))*(y(m)-y(i)) +(z(j)-z(i))*(z(m)-z(i))
        scalar2 = (x(i)-x(j))*(x(m)-x(j)) +(y(i)-y(j))*(y(m)-y(j)) +(z(i)-z(j))*(z(m)-z(j))
        scalar3 = (x(i)-x(m))*(x(j)-x(m)) +(y(i)-y(m))*(y(j)-y(m)) +(z(i)-z(m))*(z(j)-z(m))
        !  calculate cosine of the angle between bonds ij and im
        cote1 = scalar1/(r1*r3)
        !  calculate cosine of the angle between bonds ji and jm
        cote2 = scalar2/(r1*r2)
        !  calculate cosine of the angle between bonds mi and mj
        cote3 = scalar3/(r3*r2)
    end subroutine scos
end module chempotpy_aln_peatht_pes

!> @brief Preserve the historical Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_peatht_pes, only: evaluate
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
