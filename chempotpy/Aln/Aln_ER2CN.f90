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
        coordinates_angstrom(1:atom_count,:), coordinates_bohr )
    call pot( coordinates_bohr(:,1), coordinates_bohr(:,2), &
        coordinates_bohr(:,3), energy_hartree, atom_count, atom_count )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*,*) 'Only energy is available'
    endif
end subroutine pes

!> @brief Evaluate pot for the Aln_ER2CN potential.
subroutine pot( x, y, z, v, natom, maxatom )

    !   System:                        Aln
    !   Functional form:               Extended Rydberg plus Coordination Number
    !   Common name:                   ER2+CN
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

    double precision :: a1, a2, a3, abscz1, abscz2, au2, bu2, co, cz1, cz2, de, &
        & du2, dx, dy, dz, ec1, ee, es1, fu2, g1, g1e, g2, gcorr, gu2, gzero, poly, r, re, rho, &
        & rij, rr, term, v, x, y, z
    integer :: i, j, maxatom, natom
    dimension x(maxatom),y(maxatom),z(maxatom),gcorr(maxatom)
    dimension fu2(maxatom,maxatom),r(maxatom,maxatom)
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0

    de  = 1.71013678553981441/autoev
    re  = 5.08182706399609163
    a1  = 1.24074255007327805
    a2  = 0.551880801172447422
    a3  = 0.129970688812896917
    au2 = 0.143243771372740580
    bu2 = 6.50000000000000000/autoang

    do 3 i=1,natom
        r(i,i)   = 0.d0
        fu2(i,i) = 1.d0
        do 3 j=i+1,natom
            dx       = x(i)-x(j)
            dy       = y(i)-y(j)
            dz       = z(i)-z(j)
            r(i,j)   = dsqrt(dx*dx+dy*dy+dz*dz)
            r(j,i)   = r(i,j)
            fu2(j,i) = 1.d0
            fu2(i,j) = 1.d0
    3 continue

    g1    = 1.02042012701514406
    g2    = 6.15510503175378609/autoang
    du2   = 1.00048851978505127
    gu2   = 0.779189057156814857
    gzero = 6.49658036150464113

    ! CALCULATE CN FUNCTION
    !     first calculate g_alpha
    g1e = dexp(g1)
    do 20 i=1,natom                         ! alpha loop
        gcorr(i) = 0d0
        do 21 j=1,natom                        ! beta loop
            if(i.eq.j) go to 21                   ! beta != alpha
            rij = r(i,j)
            if (rij.ge.g2) go to 21               ! cutoff term counts nearb
            term     = dexp(-g1/(1d0-(rij/g2)))         ! cutoff term counts nearb
            gcorr(i) = gcorr(i)+term                ! sum over beta
        21 continue
        gcorr(i) = gcorr(i)*g1e                  ! g_alpha
    20 continue

    !     then calculate CN term
    do 30 i=1,natom                         ! alpha loop
        do 31 j=1,natom                         ! beta loop
            if(i.eq.j) go to 31                    ! beat != alpha
            rij = r(i,j)

            term = 0.d0
            if (rij.ge.g2) go to 32                ! g_alphabeta
            term = dexp(g1*(1d0-1d0/(1d0-(rij/g2)))) ! g_alphabeta
            32 continue

            cz1    = (gcorr(i)-term)/gzero              ! ratio of coordination fo
            abscz1 = dabs(cz1)
            if(abscz1.le.10d-12)cz1=abscz1         ! eliminates negative valu
            es1 = 1d0/(1d0+cz1**gu2)                 ! compute alpha part of f

            cz2    = (gcorr(j)-term)/gzero              ! ratio of coordination fo
            abscz2 = dabs(cz2)
            if(abscz2.le.10d-12)cz2=abscz2         ! eliminates negative valu
            ec1 = 1d0/(1d0+cz2**gu2)                 ! compute beta part of f f

            fu2(i,j) = 1d0+du2*(es1*ec1-1d0)         ! put alpha and beta parts
            !       print *,es1,ec1,fu2(i,j)
        31 continue
    30 continue

    ! potential, includes cutoff
    v = 0.d0
    do 1 i=1,natom
        do 1 j=i+1,natom
            dx   = x(i)-x(j)
            dy   = y(i)-y(j)
            dz   = z(i)-z(j)
            rr   = dsqrt(dx*dx+dy*dy+dz*dz)
            rho  = rr-re
            poly = 1.d0+a1*rho+a2*rho**2+a3*rho**3
            ee   = dexp(-a1*rho)
            co   = 0.d0
            if (rr.le.bu2) co=dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
            v = v-de*ee*poly*co*fu2(i,j)
    1 continue

end subroutine pot
