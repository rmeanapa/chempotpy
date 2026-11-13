!   System:                        Aln
!   Functional form:               Extended Rydberg plus Embedded Atom
!   Common name:                   ER2+MeiD
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  The parameters
!           for this PEF have been parameterized against a data set of
!           aluminum cluster energies and bulk data.
!           This PEF has a cutoff at 6.5 angstroms.
!           Note that eq 37 of J. Phys. Chem. B 108, 8896 (2004) has two
!           typos related to this potential. A minus sign is missing in
!           front of phi-sub-0, and the sign in front of the variable "d"
!           should be "+".
!  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!              and D. G. Truhlar, "Analytic Potentials Energy Functions
!              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!                    of the system
!  natom   --- number of atoms in the system
!  maxatom --- maximal number of atoms

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_er2meid_pes
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

    !> @brief Evaluate pot for the Aln_ER2MeiD potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a1, a2, a3, ab, alpha, au2, beta, bigef, bu2, cc(6), co, de, delta, dx, dy, dz, ec, ee, ef, gamma, &
            & phi0, poly, q, r0, rc, re, rho, rho_pair, rn, rr, s(3), sterm, v, v_pair, ver2, x(maxatom), xm, xx, y(maxatom), &
            & z(maxatom)
        integer             :: i, j, l, m, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        ! accurate two-body
        de  = 1.71013678553981441/autoev
        re  = 5.08182706399609163
        a1  = 1.24074255007327805
        a2  = 0.551880801172447422
        a3  = 0.129970688812896917
        au2 = 0.143243771372740580
        bu2 = 6.50000000000000000/autoang
        !     potential, includes cutoff
        ver2 = 0.e0_wp
        do i=1, natom
            do j=i+1, natom
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = sqrt(dx*dx+dy*dy+dz*dz)
                if( rr .le. bu2 )then
                    co   = exp(au2*(1.e0_wp-1.e0_wp/(1.e0_wp-rr/bu2)))
                    rho  = rr-re
                    poly = 1.e0_wp+a1*rho+a2*rho**2+a3*rho**3
                    ee   = exp(-a1*rho)
                    ver2 = ver2-de*ee*poly*co
                endif
            enddo
        enddo
        !  calculate the potential energy
        ec     = 2.00214948705422557/autoev
        phi0   = 0.110210063507572065/autoev
        r0     = 2.15915974596971205/autoang
        alpha  = 7.06277479237909134
        beta   = 6.58251099169516340
        gamma  = 6.24914509037616028
        delta  = 6.39066927210552027
        cc(1)  = 0.407034684904738653
        cc(2)  = -6.44787493893502628
        cc(3)  = 36.1221299462628238
        cc(4)  = -47.8456277479237926
        cc(5)  = 25.8260869565217384
        cc(6)  = -5.49975574010747437
        s(1)   = 12.6893014167073765
        s(2)   = 6.47234978016609652
        s(3)   = 26.9702002931118692
        rn     = 1.75e0_wp*r0
        rc     = 1.95e0_wp*r0
        v      = 0.e0_wp
        v_pair = 0.e0_wp
        do i=1, natom
            rho = 0.e0_wp
            do j=1, natom
                if( i .eq. j ) cycle
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = sqrt(dx*dx+dy*dy+dz*dz)
                ! cutoff function
                if( rr .le. rn )then
                    q = 1.e0_wp
                elseif( rr .ge. rc )then
                    q = 0.e0_wp
                else
                    xx = (rr-rn)/(rc-rn)
                    q  = ((1.e0_wp-xx)**3)*(1.e0_wp+3.e0_wp*xx+6.e0_wp*xx**2)
                endif
                ! embedding potential rho
                rho_pair = 0.e0_wp
                do l=0, 5
                    ef       = (cc(l+1)/12.e0_wp)*(r0/rr)**l
                    rho      = rho + ef*q
                    rho_pair = rho_pair + ef*q
                enddo
                sterm = 0.e0_wp
                if( rho_pair .gt. 0.e0_wp )then
                    do m=1, 3
                        xm = sqrt(real(m, wp))
                        sterm=sterm+s(m)*exp(-gamma*(xm-1.e0_wp)) *(1.e0_wp+delta*(xm-1.e0_wp)-xm*(delta/beta)* &
                            & log(rho_pair)) * (rho_pair)**(xm*gamma/beta)
                    enddo
                endif
                ab    = alpha/beta
                bigef = 0.e0_wp
                if( rho_pair .gt. 0.e0_wp )then
                    bigef = -ec*(1.e0_wp-ab*log(rho_pair))*rho_pair**ab + 0.5e0_wp*phi0*sterm
                endif
                v_pair = v_pair + bigef
            enddo
            sterm = 0.e0_wp
            if( rho .gt. 0.e0_wp )then
                do m=1, 3
                    xm = sqrt(real(m, wp))
                    sterm=sterm+s(m)*exp(-gamma*(xm-1.e0_wp)) *(1.e0_wp+delta*(xm-1.e0_wp)-xm*(delta/beta)*log(rho)) * &
                        & (rho)** (xm* gamma/beta)
                enddo
            endif
            ab    = alpha/beta
            bigef = 0.e0_wp
            if( rho .gt. 0.e0_wp )then
                bigef = -ec*(1.e0_wp-ab*log(rho))*rho**ab + 0.5e0_wp*phi0*sterm
            endif
            v = v + bigef
        enddo
        v = v + ver2 - v_pair
    end subroutine pot
end module chempotpy_aln_er2meid_pes

!> @brief Expose Aln_ER2MeiD through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_er2meid_pes, only: evaluate
    use chempotpy_kinds,           only: wp
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
