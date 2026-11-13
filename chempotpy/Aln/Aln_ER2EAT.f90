!   System:                        Aln
!   Functional form:               Extended Rydberg plus the Extended Axilrod-Teller term
!   Common name:                   ER2+EAT
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  The parameters
!           for this PEF have been parameterized against a data set of
!           aluminum cluster energies and bulk data.
!           This PEF has a cutoff at 6.5 angstroms.
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
module chempotpy_aln_er2eat_pes
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

    !> @brief Evaluate pot for the Aln_ER2EAT potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a1, a2, a3, arc1, arc2, au2, au23, b1, b2, bra, bu2, co, co1, co2, co3, coij, coim, cojm, cot12, &
            & cot2, cot3, de, dx, dy, dz, ee, poly, re, rho, rij, rim, rjm, rr, rr1, rr2, rr3, u3, v, w2, x(maxatom), &
            & y(maxatom), z(maxatom)
        integer             :: i, j, m, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
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
        v = 0.e0_wp
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
                    v    = v-de*ee*poly*co
                endif
            enddo
        enddo
        !  2. calculate energy of three-body interactions
        u3 = 0.e0_wp
        if( natom /= 2 )then
            do i=1, natom
                do j=1, natom
                    if( i .eq. j ) cycle
                    dx  = x(i)-x(j)
                    dy  = y(i)-y(j)
                    dz  = z(i)-z(j)
                    rij = sqrt(dx*dx+dy*dy+dz*dz)
                    if( rij .ge. bu2 ) cycle
                    do m=1, natom
                        if( i .eq. m ) cycle
                        if( j .eq. m ) cycle
                        dx  = x(m)-x(j)
                        dy  = y(m)-y(j)
                        dz  = z(m)-z(j)
                        rjm = sqrt(dx*dx+dy*dy+dz*dz)
                        dx  = x(m)-x(i)
                        dy  = y(m)-y(i)
                        dz  = z(m)-z(i)
                        rim = sqrt(dx*dx+dy*dy+dz*dz)
                        if( rjm .ge. bu2 ) cycle
                        if( rim .ge. bu2 ) cycle
                        coij = 0.e0_wp
                        cojm = 0.e0_wp
                        coim = 0.e0_wp
                        if( rij .lt. bu2 ) coij=exp(au23*(1.e0_wp-1.e0_wp/(1.e0_wp-rij/bu2)))
                        if( rjm .lt. bu2 ) cojm=exp(au23*(1.e0_wp-1.e0_wp/(1.e0_wp-rjm/bu2)))
                        if( rim .lt. bu2 ) coim=exp(au23*(1.e0_wp-1.e0_wp/(1.e0_wp-rim/bu2)))
                        co = coij*cojm*coim
                        call scos(i, j, m, x, y, z, rij, rjm, rim, natom, maxatom, co1, co2, co3)
                        ! at term
                        u3 = u3+b1*(1.e0_wp+3.e0_wp*co1*co2*co3)/((rij*rjm*rim)**3)*co
                        ! next term
                        if( (abs(co1)-1.e0_wp) .le. 10e-12_wp )then
                            arc1 = 0.e0_wp
                        else
                            arc1 = acos(co1)
                        endif
                        if( (abs(co2)-1.e0_wp) .le. 10e-12_wp )then
                            arc2 = 0.e0_wp
                        else
                            arc2 = acos(co2)
                        endif
                        cot12 = 6.e0_wp*(sin(arc1)*sin(arc2)+co1*co2)
                        cot3  = -25.e0_wp*(4.e0_wp*(co3)**3-3.e0_wp*co3)
                        cot2  = 3.e0_wp+5.e0_wp*(2.e0_wp*(co3)**2-1.e0_wp)
                        bra   = 9.e0_wp*co3+cot3+cot12*cot2
                        rr1   = 1.e0_wp/rij
                        rr2   = 1.e0_wp/rim
                        rr3   = 1.e0_wp/rjm
                        w2    = co*b2*bra*(rr1**3)*(rr2*rr3)**4
                        u3    = u3+w2
                    enddo
                enddo
            enddo
            !  calculate the potential energy
        endif
        v = v+u3
    end subroutine pot

    !> @brief Evaluate scos for the Aln_ER2EAT potential.
    subroutine scos( i, j, m, x, y, z, rij, rjm, rim, natom, maxatom, &
            & co1, co2, co3 )
        !  The subroutine calculates the cosines of the bound angles
        !  using scalar products of the position vectors of two
        !  atoms originating from the third atom.
        !  Units:
        !  coordinates - bohrs
        !  i,j,m   --- natural indices running from 1 to natom
        !  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
        !              of the system
        !  r       --- two-dimensional array of atomic distances
        !  natom   --- number of atoms in the system
        !  maxatom --- maximal number of atoms
        !  co1     --- cosine of the angle between bonds ij and im
        !  co2     --- cosine of the angle between bonds ji and jm
        !  co3     --- cosine of the angle between bonds mi and mj
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: co1, co2, co3, rij, rim, rjm, scalar1, scalar2, scalar3, x(maxatom), y(maxatom), z(maxatom)
        integer  :: i, j, m, maxatom, natom
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
end module chempotpy_aln_er2eat_pes

!> @brief Expose Aln_ER2EAT through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_er2eat_pes, only: evaluate
    use chempotpy_kinds,          only: wp
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
