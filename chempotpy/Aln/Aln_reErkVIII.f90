!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_reerkviii_pes
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
    !> @brief Evaluate pot for the Aln_reErkVIII potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_geometry, only: fill_distance_matrix
        !   System:                        Aln
        !   Functional form:               Lennard-Jones plus Axilrod-Teller
        !   Common name:                   reErkVIII
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
        !              S. Erkoc, Docga(A1), 9, 203 (1985).
        !              S. Erkoc, Empirical Potential Energy Functions
        !                 Used in the Simulations of Materials Properties,
        !                 Annual Reviews of Computational Physics IX,
        !                 D. Stauffer, ed., World Scientific, Singapore 2001,
        !                 pp. 1-103.
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
        double precision :: a1, a2, arc1, arc2, b1, b2, bra, cot12, cot13, cot2, cot3, &
            & cote1, cote2, cote3, dx, dy, dz, em, en, r, r0, rel, rr, rr1, rr2, rr3, u, uem, uen, v, v1, &
            & w1, w2, war, x, y, z
        integer :: i, j, m, maxatom, natom
        dimension x(maxatom), y(maxatom), z(maxatom), r(maxatom,maxatom)
        double precision, parameter :: autoev  = 27.2113961d0
        double precision, parameter :: autoang = 0.529177249d0
        !  calculate and tabulate the distances between all atoms
        call fill_distance_matrix(x(1:natom), y(1:natom), z(1:natom), r(1:natom,1:natom))
        !
        !  calculate the potential energy
        !
        ! original parameters
        r0 = 2.510d0  /  autoang
        a1 = 1.860d0  /  autoev
        a2 = 3.410d0  /  autoev
        b1 = 1643.d0  / (autoev*autoang**9)
        b2 = 1921.d0  / (autoev*autoang**11)
        em = 11.d0
        en = 6.d0
        ! new parameters
        r0 = 2.510d0  /  autoang
        a1 = 1.49145090376160239/autoev
        a2 = 3.11382510991695138/autoev
        b1 = 2522.22765021983378/(autoev*autoang**9)
        b2 = 86.8978993649242710/ (autoev*autoang**11)
        em = 6.90425012212994638
        en = 4.84269662921348321
        u  = 0.d0
        do i=1, natom-1
            do j=i+1, natom
                rr  = r(i,j)
                rel = r0/rr
                uem = a1*rel**em
                uen = a2*rel**en
                u   = u+(uem-uen)
            enddo
        enddo
        v1 = 0.d0
        if( natom /= 2 )then
            do i=1, natom-2
                do j=i+1, natom-1
                    do m=j+1, natom
                        call scos(i, j, m, x, y, z, r, natom, maxatom, cote1, cote2, cote3)
                        rr1   = 1.d0/r(i,j)
                        rr2   = 1.d0/r(i,m)
                        rr3   = 1.d0/r(j,m)
                        cot13 = 1.d0+3.d0*cote1*cote2*cote3
                        war   = dabs(cote1)-1.0d0
                        if( war <= 1.0d-11 )then
                            arc1 = 0.0d0
                        else
                            arc1 = dacos(cote1)
                        endif
                        war = dabs(cote2)-1.0d0
                        if( war <= 1.0d-11 )then
                            arc2 = 0.0d0
                        else
                            arc2 = dacos(cote2)
                        endif
                        cot12 = 6.d0*(dsin(arc1)*dsin(arc2)+cote1*cote2)
                        cot3  = -25.d0*(4.d0*(cote3)**3-3.d0*cote3)
                        cot2  = 3.d0+5.d0*(2.d0*(cote3)**2-1.d0)
                        bra   = 9.d0*cote3+cot3+cot12*cot2
                        w1    = b1*cot13*(rr1*rr2*rr3)**3     ! AT term
                        w2    = b2*bra*(rr1**3)*(rr2*rr3)**4
                        v1    = v1+w1+w2
                    enddo
                enddo
            enddo
        endif
        v = u+v1
    end subroutine pot
    !
    !
    !> @brief Evaluate scos for the Aln_reErkVIII potential.
    subroutine scos( i, j, m, x, y, z, r, natom, maxatom, cote1, cote2, cote3 )
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
        double precision :: cote1, cote2, cote3, r, scalar1, scalar2, scalar3, x, y, z
        integer          :: i, j, m, maxatom, natom
        dimension x(maxatom), y(maxatom), z(maxatom), r(maxatom,maxatom)
        !  calculate scalar products of the position vectors of two atoms
        !  originating from the third atom
        scalar1 = (x(j)-x(i))*(x(m)-x(i)) +(y(j)-y(i))*(y(m)-y(i)) +(z(j)-z(i))*(z(m)-z(i))
        scalar2 = (x(i)-x(j))*(x(m)-x(j)) +(y(i)-y(j))*(y(m)-y(j)) +(z(i)-z(j))*(z(m)-z(j))
        scalar3 = (x(i)-x(m))*(x(j)-x(m)) +(y(i)-y(m))*(y(j)-y(m)) +(z(i)-z(m))*(z(j)-z(m))
        !  calculate cosine of the angle between bonds ij and im
        cote1 = scalar1/(r(i,j)*r(i,m))
        !  calculate cosine of the angle between bonds ji and jm
        cote2 = scalar2/(r(i,j)*r(j,m))
        !  calculate cosine of the angle between bonds mi and mj
        cote3 = scalar3/(r(i,m)*r(j,m))
    end subroutine scos
end module chempotpy_aln_reerkviii_pes

!> @brief Preserve the historical Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_reerkviii_pes, only: evaluate
    use chempotpy_kinds,             only: real64
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
