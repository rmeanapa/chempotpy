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

!> @brief Evaluate pot for the Aln_ER2EAT potential.
subroutine pot( x, y, z, v, natom, maxatom )
    !   System:                        Aln
    !   Functional form:               Extended Rydberg plus the Extended Axilrod-Teller term
    !   Common name:                   ER2+EAT
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
    double precision :: a1, a2, a3, arc1, arc2, au2, au23, b1, b2, bra, bu2, co, &
        & co1, co2, co3, coij, coim, cojm, cot12, cot2, cot3, de, dx, dy, dz, ee, poly, re, rho, rij, &
        & rim, rjm, rr, rr1, rr2, rr3, u3, v, w2, x, y, z
    integer :: i, j, m, maxatom, natom
    dimension x(maxatom), y(maxatom), z(maxatom)
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
    b1   = 569.076697606  / (autoev*autoang**9)
    b2   = 5.0219833903273/ (autoev*autoang**11)
    au23 = 0.164763067904250132e-01
    de   = 1.71013678553981441/autoev
    re   = 5.08182706399609163
    a1   = 1.24074255007327805
    a2   = 0.551880801172447422
    a3   = 0.129970688812896917
    au2  = 0.143243771372740580
    bu2  = 6.50000000000000000/autoang
    !     potential, includes cutoff
    v = 0.d0
    do 1 i=1, natom
        do 1 j=i+1, natom
            dx = x(i)-x(j)
            dy = y(i)-y(j)
            dz = z(i)-z(j)
            rr = dsqrt(dx*dx+dy*dy+dz*dz)
            if( rr .le. bu2 )then
                co   = dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
                rho  = rr-re
                poly = 1.d0+a1*rho+a2*rho**2+a3*rho**3
                ee   = dexp(-a1*rho)
                v    = v-de*ee*poly*co
            endif
    1 continue
    !  2. calculate energy of three-body interactions
    u3 = 0.d0
    if( natom .eq. 2 ) go to 7
    do 6 i=1, natom
        do 4 j=1, natom
            if( i .eq. j ) go to 4
            dx  = x(i)-x(j)
            dy  = y(i)-y(j)
            dz  = z(i)-z(j)
            rij = dsqrt(dx*dx+dy*dy+dz*dz)
            if( rij .ge. bu2 ) go to 4
            do 5 m=1, natom
                if( i .eq. m ) go to 5
                if( j .eq. m ) go to 5
                dx  = x(m)-x(j)
                dy  = y(m)-y(j)
                dz  = z(m)-z(j)
                rjm = dsqrt(dx*dx+dy*dy+dz*dz)
                dx  = x(m)-x(i)
                dy  = y(m)-y(i)
                dz  = z(m)-z(i)
                rim = dsqrt(dx*dx+dy*dy+dz*dz)
                if( rjm .ge. bu2 ) go to 5
                if( rim .ge. bu2 ) go to 5
                coij = 0.d0
                cojm = 0.d0
                coim = 0.d0
                if( rij .lt. bu2 ) coij=dexp(au23*(1.d0-1.d0/(1.d0-rij/bu2)))
                if( rjm .lt. bu2 ) cojm=dexp(au23*(1.d0-1.d0/(1.d0-rjm/bu2)))
                if( rim .lt. bu2 ) coim=dexp(au23*(1.d0-1.d0/(1.d0-rim/bu2)))
                co = coij*cojm*coim
                call scos(i, j, m, x, y, z, rij, rjm, rim, natom, maxatom, co1, co2, co3)
                ! at term
                u3 = u3+b1*(1.d0+3.d0*co1*co2*co3)/((rij*rjm*rim)**3)*co
                ! next term
                if( (dabs(co1)-1.d0) .le. 10d-12 )then
                    arc1 = 0.d0
                else
                    arc1 = dacos(co1)
                endif
                if( (dabs(co2)-1.d0) .le. 10d-12 )then
                    arc2 = 0.d0
                else
                    arc2 = dacos(co2)
                endif
                cot12 = 6.d0*(dsin(arc1)*dsin(arc2)+co1*co2)
                cot3  = -25.d0*(4.d0*(co3)**3-3.d0*co3)
                cot2  = 3.d0+5.d0*(2.d0*(co3)**2-1.d0)
                bra   = 9.d0*co3+cot3+cot12*cot2
                rr1   = 1.d0/rij
                rr2   = 1.d0/rim
                rr3   = 1.d0/rjm
                w2    = co*b2*bra*(rr1**3)*(rr2*rr3)**4
                u3    = u3+w2
            5 continue
        4 continue
    6 continue
    !  calculate the potential energy
    7 v=v+u3
end subroutine pot

!
!
!
!> @brief Evaluate scos for the Aln_ER2EAT potential.
subroutine scos( i, j, m, x, y, z, rij, rjm, rim, natom, maxatom, &
        & co1, co2, co3 )
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
    !  co1     --- cosine of the angle between bonds ij and im
    !  co2     --- cosine of the angle between bonds ji and jm
    !  co3     --- cosine of the angle between bonds mi and mj
    !
    implicit none
    double precision :: co1, co2, co3, rij, rim, rjm, scalar1, scalar2, scalar3, x, y, z
    integer          :: i, j, m, maxatom, natom
    dimension x(maxatom), y(maxatom), z(maxatom)
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
