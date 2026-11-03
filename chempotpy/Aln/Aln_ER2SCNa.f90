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
        write(*, *) 'Only energy is available'
    end if
end subroutine pes

!> @brief Evaluate pot for the Aln_ER2SCNa potential.
subroutine pot( x, y, z, v, natom, maxatom )
    !   System:                        Aln
    !   Functional form:               Extended Rydberg plus Screening and Coordination Number
    !   Common name:                   ER2+SCNa
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
    double precision :: a1, a2, a3, abscz1, abscz2, ac, au2, au23, bu2, c, co, &
        & cz1, cz2, de, du2, dx, dy, dz, ec1, ee, es1, fcn, fs, g1, g1e, g2, gcorr, gu2, gzero, &
        & poly, r, re, rho, rij, ril, rjl, rr, s, term, tmp, v, x, xk1, xk2, xk3, y, z
    integer :: i, j, l, maxatom, natom
    dimension x(maxatom), y(maxatom), z(maxatom), gcorr(maxatom)
    dimension fs(maxatom,maxatom), r(maxatom,maxatom)
    dimension fcn(maxatom,maxatom), c(10)
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
    de  = 1.71013678553981441/autoev
    re  = 5.08182706399609163
    a1  = 1.24074255007327805
    a2  = 0.551880801172447422
    a3  = 0.129970688812896917
    au2 = 0.143243771372740580
    bu2 = 6.50000000000000000/autoang

    xk1   = 2.44335613092330250
    xk2   = 0.153395212506106493
    xk3   = 4.42027357107962882
    au23  = 0.953541768441621973
    g1    = 0.725354176844162146
    g2    = 6.49853444064484620/autoang
    du2   = 0.774548119198827489
    gu2   = 0.799364924279433287
    gzero = 6.48387884709330731

    do 3 i=1, natom
        r(i,i)   = 0.d0
        fs(i,i)  = 1.d0
        fcn(i,i) = 1.d0
        do 3 j=i+1, natom
            dx       = x(i)-x(j)
            dy       = y(i)-y(j)
            dz       = z(i)-z(j)
            r(i,j)   = dsqrt(dx*dx+dy*dy+dz*dz)
            r(j,i)   = r(i,j)
            fs(j,i)  = 1.d0
            fs(i,j)  = 1.d0
            fcn(j,i) = 1.d0
            fcn(i,j) = 1.d0
    3 continue

    ! do S term, needs xk1,xk2,xk3
    do 510 i=1, natom                          ! alpha loop
        do 515 j=i+1, natom                      ! beta loop

            s = 0.d0
            do 512 l=1, natom                      ! gamma loop
                if(l.eq.j.or.l.eq.i) go to 512        ! gamma != alpha or beta
                rjl = r(l,j)
                rij = r(i,j)
                ril = r(l,i)
                if(rjl.ge.bu2) go to 512              ! cutoff for Rbeta,delta
                if(rij.ge.bu2) go to 512              ! cutoff for Rbeta,delta
                if(ril.ge.bu2) go to 512              ! cutoff for Rbeta,delta
                ac   = au23*(1d0-1d0/(1d0-rjl/bu2))       ! cutoff for Rbeta,delta
                ac   = ac+au23*(1d0-1d0/(1d0-rij/bu2))    ! cutoff for Rbeta,delta
                ac   = ac+au23*(1d0-1d0/(1d0-ril/bu2))    ! cutoff for Rbeta,delta
                term = xk1*dexp(-xk2*((ril+rjl)/rij)**xk3+ac)        ! cutoff ti
                s    = s+term                         ! sum
            512 continue

            fs(i,j) = dtanh(s)
            fs(j,i) = fs(i,j)

        515 continue
        fs(i,i) = 1.d0                           ! shouldn't get used anywhe
    510 continue

    ! CALCULATE CN FUNCTION
    !     first calculate g_alpha
    g1e = dexp(g1)
    do 20 i=1, natom                         ! alpha loop
        gcorr(i) = 0d0
        do 21 j=1, natom                        ! beta loop
            if(i.eq.j) go to 21                   ! beta != alpha
            rij = r(i,j)
            if (rij.ge.g2) go to 21               ! cutoff term counts nearb
            term     = dexp(-g1/(1d0-(rij/g2)))         ! cutoff term counts nearb
            gcorr(i) = gcorr(i)+term                ! sum over beta
        21 continue
        gcorr(i) = gcorr(i)*g1e                  ! g_alpha
    20 continue

    !     then calculate CN term
    do 30 i=1, natom                         ! alpha loop
        do 31 j=1, natom                         ! beta loop
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

            fcn(i,j) = du2*(es1*ec1-1d0)         ! put alpha and beta parts tog
        31 continue
    30 continue

    ! potential, includes cutoff
    v = 0.d0
    do 1 i=1, natom
        do 1 j=i+1, natom
            dx   = x(i)-x(j)
            dy   = y(i)-y(j)
            dz   = z(i)-z(j)
            rr   = dsqrt(dx*dx+dy*dy+dz*dz)
            rho  = rr-re
            poly = 1.d0+a1*rho+a2*rho**2+a3*rho**3
            ee   = dexp(-a1*rho)
            co   = 0.d0
            if (rr.le.bu2) co=dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
            tmp = de*ee*poly*co
            v   = v-tmp*(1.d0-fs(i,j)+fcn(i,j))
    1 continue

end subroutine pot
