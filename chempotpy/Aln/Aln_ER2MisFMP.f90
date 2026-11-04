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

!> @brief Evaluate pot for the Aln_ER2MisFMP potential.
subroutine pot( x, y, z, v, natom, maxatom )
    !   System:                        Aln
    !   Functional form:               Extended Rydberg plus Embedded Atom
    !   Common name:                   ER2+MisFMP
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
    double precision :: a, a1, a2, a3, au2, au23, b1, b2, bu2, co, co3, de, dx, dy, &
        & dz, ee, ef, ef2, f0, f2, poly, q1, q2, q3, r3, r4, re, rh, rho, rhobar, rr, v, v2, x, &
        & xksum, xksum2, y, z
    integer                     :: i, j, maxatom, natom
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
    dimension x(maxatom), y(maxatom), z(maxatom)
    !
    !  calculate the potential energy
    !
    de   = 1.71013678553981441/autoev
    re   = 5.08182706399609163
    a1   = 1.24074255007327805
    a2   = 0.551880801172447422
    a3   = 0.129970688812896917
    au2  = 0.143243771372740580
    bu2  = 6.50000000000000000/autoang
    a    = 0.264359550561797718
    r3   = 3.23400097703957012/autoang
    r4   = -9.98827552515876960/autoang
    b1   = 2.00000000000000000*autoang**2
    b2   = 1.48949682462139710*autoang
    f0   = -2.08891060087933589/autoev
    f2   = 7.20566682950659487/autoev
    q1   = -4.45041524181729375/autoev
    q2   = -18.7274059599413789/autoev
    q3   = -15.3981436248168055/autoev
    au23 = 2.62012701514411361
    v    = 0.d0
    ef   = 0.d0
    ef2  = 0.d0
    do 2 i=1, natom
        v2     = 0.d0
        ef2    = 0.d0
        rhobar = 0.d0
        do 3 j=1, natom
            if( i .eq. j ) go to 3
            dx = x(i)-x(j)
            dy = y(i)-y(j)
            dz = z(i)-z(j)
            rr = dsqrt(dx*dx+dy*dy+dz*dz)
            if( rr .ge. bu2 ) go to 3
            ! 2-body part
            ! ERCO2-6.5
            rh   = rr-re
            poly = 1+a1*rh+a2*rh**2+a3*rh**3
            ee   = dexp(-a1*rh)
            co   = 0.d0
            co3  = 0.d0
            if( rr .lt. bu2 ) co=dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
            if( rr .lt. bu2 ) co3=dexp(au23*(1.d0-1.d0/(1.d0-rr/bu2)))
            v2 = -de*ee*poly*co
            ! rhobar
            rho    = (a*dexp(-b1*(rr-r3)**2)+dexp(-b2*(rr-r4)))*co3
            rhobar = rhobar + rho
            ! subtract two-body part
            xksum2 = q1*(rho-1.d0)**3 + q2*(rho-1.d0)**4 + q3*(rho-1.d0)**5
            ef2    = f0 + 0.5d0*f2*(rho-1.d0)**2 + xksum2
        3 continue
        ! embedding term
        xksum = q1*(rhobar-1.d0)**3 + q2*(rhobar-1.d0)**4 + q3*(rhobar-1.d0)**5
        ef    = f0 + 0.5d0*f2*(rhobar-1.d0)**2 + xksum
        v     = v + 0.5d0*v2 + ef - ef2
    2 continue
end subroutine pot
