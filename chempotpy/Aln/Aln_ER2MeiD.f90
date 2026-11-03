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

!> @brief Evaluate pot for the Aln_ER2MeiD potential.
subroutine pot( x, y, z, v, natom, maxatom )
    !   System:                        Aln
    !   Functional form:               Extended Rydberg plus Embedded Atom
    !   Common name:                   ER2+MeiD
    !   Number of derivatives:         0
    !   Number of bodies:              variable
    !   Number of electronic surfaces: 1
    !   Interface:                     HO-MM-0
    !
    !   Notes:  Many-body aluminum potential energy function.  The parameters
    !           for this PEF have been parameterized against a data set of
    !           aluminum cluster energies and bulk data.
    !           This PEF has a cutoff at 6.5 angstroms.
    !           Note that eq 37 of J. Phys. Chem. B 108, 8896 (2004) has two
    !           typos related to this potential. A minus sign is missing in
    !           front of phi-sub-0, and the sign in front of the variable "d"
    !           should be "+".
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
    double precision :: a1, a2, a3, ab, alpha, au2, beta, bigef, bu2, cc, co, de, &
        & delta, dx, dy, dz, ec, ee, ef, gamma, phi0, poly, q, r0, rc, re, rho, rho_pair, rn, rr, s, &
        & sterm, v, v_pair, ver2, x, xm, xx, y, z
    integer                     :: i, j, l, m, maxatom, natom
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
    dimension x(maxatom), y(maxatom), z(maxatom)
    dimension cc(6), s(3)
    ! accurate two-body
    de  = 1.71013678553981441/autoev
    re  = 5.08182706399609163
    a1  = 1.24074255007327805
    a2  = 0.551880801172447422
    a3  = 0.129970688812896917
    au2 = 0.143243771372740580
    bu2 = 6.50000000000000000/autoang
    !     potential, includes cutoff
    ver2 = 0.d0
    do 1 i=1, natom
        do 1 j=i+1, natom
            dx = x(i)-x(j)
            dy = y(i)-y(j)
            dz = z(i)-z(j)
            rr = dsqrt(dx*dx+dy*dy+dz*dz)
            if (rr.le.bu2) then
                co   = dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
                rho  = rr-re
                poly = 1.d0+a1*rho+a2*rho**2+a3*rho**3
                ee   = dexp(-a1*rho)
                ver2 = ver2-de*ee*poly*co
            end if
    1 continue

    !  calculate the potential energy
    ec    = 2.00214948705422557/autoev
    phi0  = 0.110210063507572065/autoev
    r0    = 2.15915974596971205/autoang
    alpha = 7.06277479237909134
    beta  = 6.58251099169516340
    gamma = 6.24914509037616028
    delta = 6.39066927210552027
    cc(1) = 0.407034684904738653
    cc(2) = -6.44787493893502628
    cc(3) = 36.1221299462628238
    cc(4) = -47.8456277479237926
    cc(5) = 25.8260869565217384
    cc(6) = -5.49975574010747437
    s(1)  = 12.6893014167073765
    s(2)  = 6.47234978016609652
    s(3)  = 26.9702002931118692
    rn    = 1.75d0*r0
    rc    = 1.95d0*r0

    v      = 0.d0
    v_pair = 0.d0
    do 2 i=1, natom
        rho = 0.d0
        do 3 j=1, natom
            if (i.eq.j) go to 3
            dx = x(i)-x(j)
            dy = y(i)-y(j)
            dz = z(i)-z(j)
            rr = dsqrt(dx*dx+dy*dy+dz*dz)

            ! cutoff function
            if (rr.le.rn) then
                q = 1.d0
            elseif(rr.ge.rc) then
                q = 0.d0
            else
                xx = (rr-rn)/(rc-rn)
                q  = ((1.d0-xx)**3)*(1.d0+3.d0*xx+6.d0*xx**2)
            end if

            ! embedding potential rho
            rho_pair = 0.d0
            do l=0, 5
                ef       = (cc(l+1)/12.d0)*(r0/rr)**l
                rho      = rho + ef*q
                rho_pair = rho_pair + ef*q
            end do

            sterm = 0.d0
            if (rho_pair.gt.0.d0) then
                do m=1, 3
                    xm = dsqrt(dble(m))
                    sterm=sterm+s(m)*dexp(-gamma*(xm-1.d0)) *(1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho_pair)) * &
                        & (rho_pair)**(xm*gamma/beta)
                end do
            end if

            ab    = alpha/beta
            bigef = 0.d0
            if (rho_pair.gt.0.d0) then
                bigef = -ec*(1.d0-ab*dlog(rho_pair))*rho_pair**ab + 0.5d0*phi0*sterm
            end if

            v_pair = v_pair + bigef

        3 continue

        sterm = 0.d0
        if (rho.gt.0.d0) then
            do m=1, 3
                xm = dsqrt(dble(m))
                sterm=sterm+s(m)*dexp(-gamma*(xm-1.d0)) *(1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho)) *(rho)**(xm* &
                    & gamma/beta)
            end do
        end if

        ab    = alpha/beta
        bigef = 0.d0
        if (rho.gt.0.d0) then
            bigef = -ec*(1.d0-ab*dlog(rho))*rho**ab + 0.5d0*phi0*sterm
        end if
        v = v + bigef
    2 continue

    v = v + ver2 - v_pair

end subroutine pot
