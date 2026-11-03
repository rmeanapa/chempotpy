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

!> @brief Evaluate pot for the Aln_ER2S potential.
subroutine pot( x, y, z, v, natom, maxatom )

    !   System:                        Aln
    !   Functional form:               Extended Rydberg plus Screening
    !   Common name:                   ER2+S
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


    double precision :: a1, a2, a3, ac, au2, au23, bu2, co, de, dx, dy, dz, ee, &
        & fu2, poly, r, re, rho, rij, ril, rjl, rr, s, term, v, x, xk1, xk2, xk3, y, z
    integer :: i, j, l, maxatom, natom
    dimension x(maxatom),y(maxatom),z(maxatom)
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

    xk1  = 2.14338055691255480
    xk2  = 0.456033219345383489
    xk3  = 2.50512945774303875
    au23 = 0.439912066438690741e-01

    ! do BA term, needs xk1,xk2,xk3
    do 10 i=1,natom                          ! alpha loop
        do 15 j=i+1,natom                      ! beta loop

            s = 0.d0
            do 12 l=1,natom                      ! gamma loop
                if(l.eq.j.or.l.eq.i) go to 12        ! gamma != alpha or beta
                rjl = r(l,j)
                rij = r(i,j)
                ril = r(l,i)
                if(rjl.ge.bu2) go to 12              ! cutoff for Rbeta,delta
                if(rij.ge.bu2) go to 12              ! cutoff for Rbeta,delta
                if(ril.ge.bu2) go to 12              ! cutoff for Rbeta,delta
                ac   = au23*(1d0-1d0/(1d0-rjl/bu2))       ! cutoff for Rbeta,delta
                ac   = ac+au23*(1d0-1d0/(1d0-rij/bu2))    ! cutoff for Rbeta,delta
                ac   = ac+au23*(1d0-1d0/(1d0-ril/bu2))    ! cutoff for Rbeta,delta
                term = xk1*dexp(-xk2*((ril+rjl)/rij)**xk3+ac)        ! cutoff ti
                s    = s+term                         ! sum
            12 continue

            !        e2=dexp(-2d0*s)                       ! equivalent to 1-dtanh(s)
            !        fu2(i,j)=2d0*e2/(1d0+e2)              ! equivalent to 1-dtanh(s)
            fu2(i,j) = 1.d0-dtanh(s)
            fu2(j,i) = fu2(i,j)

        15 continue
        fu2(i,i) = 1.d0                           ! shouldn't get used anywh
    10 continue

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
