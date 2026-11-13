!   System:                        Aln
!   Functional form:               Embedded atom
!   Common name:                   CleR
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  For a recent
!           discussion of this potential, see:
!           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!           and D. G. Truhlar, "Analytic Potentials Energy Functions
!           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
!  References: F. Cleri and V. Rosato, Phys. Rev. B 48, 22 (1993).
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximal number of atoms (input)

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_cler_pes
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

    !> @brief Evaluate pot for the Aln_CleR potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a, a0, dx, dy, dz, p, q, r0, rel, rr, si, u, ud, udn, udw, ur, urr, v, x(maxatom), y(maxatom), &
            & z(maxatom)
        integer             :: i, j, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        !  calculate the potential energy
        v  = 0e0_wp
        a  = 0.1221/autoev
        si = 1.316/autoev
        p  = 8.612
        q  = 2.516
        a0 = 4.050/autoang
        r0 = a0/sqrt(2e0_wp)
        u  = 0e0_wp
        do i=1, natom-1
            do j=i+1, natom
                dx  = x(i)-x(j)
                dy  = y(i)-y(j)
                dz  = z(i)-z(j)
                rr  = sqrt(dx*dx+dy*dy+dz*dz)
                rel = rr/r0
                ur  = 0e0_wp
                ur  = a*exp(-p*(rel-1e0_wp))
                u   = u+ur
            enddo
        enddo
        ud = 0e0_wp
        do i=1, natom
            udw = 0e0_wp
            do j=1, natom
                urr = 0e0_wp
                dx  = x(i)-x(j)
                dy  = y(i)-y(j)
                dz  = z(i)-z(j)
                rr  = sqrt(dx*dx+dy*dy+dz*dz)
                rel = rr/r0
                if( j .eq. i ) cycle
                urr = (si**2)*exp(-2e0_wp*q*(rel-1e0_wp))
                udw = udw+urr
            enddo
            udn = sqrt(udw)
            ud  = ud-udn
        enddo
        v = v+u+ud
    end subroutine pot
end module chempotpy_aln_cler_pes

!> @brief Expose Aln_CleR through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_cler_pes, only: evaluate
    use chempotpy_kinds,        only: wp
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
