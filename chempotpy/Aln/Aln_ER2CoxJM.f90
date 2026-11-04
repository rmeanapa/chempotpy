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

!> @brief Evaluate pot for the Aln_ER2CoxJM potential.
subroutine pot( x, y, z, v, natom, maxatom )
    !   System:                        Aln
    !   Functional form:               Extended Rydberg plus a generalized three-body term
    !   Common name:                   ER2+CoxJM
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
    double precision :: a1, a2, a3, au2, au23, bu2, c, c0, c1, c2, c3, c4, c5, c6, &
        & co, coij, cojk, coki, dd, de, dx, dy, dz, ee, fsech, fu2, oneha, onesi, oneth, poly, q1, &
        & q2, q23, q3, q332, r, re, ree, rho, rij, rjk, rki, rr, twoth, v, v3, x, y, z
    integer :: i, j, k, maxatom, natom
    dimension x(maxatom), y(maxatom), z(maxatom), c(15)
    dimension fu2(maxatom,maxatom), r(maxatom,maxatom)
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
    ! two-body
    de  = 1.71013678553981441/autoev
    re  = 5.08182706399609163
    a1  = 1.24074255007327805
    a2  = 0.551880801172447422
    a3  = 0.129970688812896917
    au2 = 0.143243771372740580
    bu2 = 6.50000000000000000/autoang
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
    ! three-body
    c0  = 2.64020517830972157
    c1  = -0.641182217879824101
    c2  = 1.87493893502686859
    c3  = -2.79115779189057145
    c4  = 0.000000000000000000e+00
    c5  = 0.000000000000000000e+00
    c6  = 0.000000000000000000e+00
    dd  = 1.00000000000000000/autoev
    ree = 2.18563751831949205/autoang
    !        xa3= 0.000000000000000000E+00
    au23  = 0.806253053248656593
    oneth = dsqrt(1.d0/3.d0)
    oneha = dsqrt(1.d0/2.d0)
    twoth = dsqrt(2.d0/3.d0)
    onesi = dsqrt(1.d0/6.d0)
    do 111 i=1, natom
        r(i,i) = 0d0
        do 111 j=i+1, natom
            dx     = x(i)-x(j)
            dy     = y(i)-y(j)
            dz     = z(i)-z(j)
            r(i,j) = dsqrt(dx*dx+dy*dy+dz*dz)
    111 r(j,i)=r(i,j)
    !  calculate the potential energy
    v3 = 0.d0
    do 2 i=1, natom
        do 3 j=i+1, natom
            !  two body term
            rr   = r(i,j)
            coij = 0.d0
            if( rr .ge. bu2 ) go to 3
            coij = dexp(au23*(1.d0-1.d0/(1.d0-rr/bu2)))
            rij  = (rr-ree)/ree
            do 4 k=j+1, natom
                ! three body term
                rr   = r(k,j)
                cojk = 0.d0
                if( rr .ge. bu2 ) go to 4
                cojk = dexp(au23*(1.d0-1.d0/(1.d0-rr/bu2)))
                rjk  = (rr-ree)/ree
                rr   = r(k,i)
                coki = 0.d0
                if( rr .ge. bu2 ) go to 4
                coki = dexp(au23*(1.d0-1.d0/(1.d0-rr/bu2)))
                rki  = (rr-ree)/ree
                q1   = oneth*(rij+rjk+rki)
                q2   = oneha*(rjk-rki)
                q3   = twoth*rij-onesi*(rjk+rki)
                q23  = q2**2+q3**2
                q332 = q3**3-3.d0*q3*q2**2
                poly = c0+c1*q1+c2*q1**2+c3*q23+c4*q1**3+c5*q1*q23+c6*q332
                !      poly=c0+c1*q1+c2*q1**2+c3*q23
                !      fsech=1.d0/dcosh(xa3*q1)
                fsech = 1.d0
                co    = coij*cojk*coki
                v3    = v3+dd*poly*co*fsech
            4 continue
        3 continue
    2 continue
    v = v+v3
end subroutine pot
