!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate an aluminum-cluster potential through the Python contract.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
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
end subroutine pes

!> @brief Evaluate pot for the Aln_CleR potential.
subroutine pot( x, y, z, v, natom, maxatom )
    !   System:                        Aln
    !   Functional form:               Embedded atom
    !   Common name:                   CleR
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
    !  References: F. Cleri and V. Rosato, Phys. Rev. B 48, 22 (1993).
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
    double precision :: a, a0, dx, dy, dz, p, q, r0, rel, rr, si, u, ud, udn, udw, &
        & ur, urr, v, x, y, z
    integer :: i, j, maxatom, natom
    dimension x(maxatom), y(maxatom), z(maxatom)
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
    !
    !  calculate the potential energy
    !
    v  = 0d0
    a  = 0.1221/autoev
    si = 1.316/autoev
    p  = 8.612
    q  = 2.516
    a0 = 4.050/autoang
    r0 = a0/dsqrt(2d0)
    u  = 0d0
    do 2 i=1, natom-1
        do 3 j=i+1, natom
            dx  = x(i)-x(j)
            dy  = y(i)-y(j)
            dz  = z(i)-z(j)
            rr  = dsqrt(dx*dx+dy*dy+dz*dz)
            rel = rr/r0
            ur  = 0d0
            ur  = a*dexp(-p*(rel-1d0))
        3 u=u+ur
    2 continue
    ud = 0d0
    do 5 i=1, natom
        udw = 0d0
        do 4 j=1, natom
            urr = 0d0
            dx  = x(i)-x(j)
            dy  = y(i)-y(j)
            dz  = z(i)-z(j)
            rr  = dsqrt(dx*dx+dy*dy+dz*dz)
            rel = rr/r0
            if( j .eq. i ) go to 41
            urr = (si**2)*dexp(-2d0*q*(rel-1d0))
            udw = udw+urr
            41 continue
        4 continue
        udn = dsqrt(udw)
    5 ud=ud-udn
    v = v+u+ud
end subroutine pot
