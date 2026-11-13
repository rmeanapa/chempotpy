!   System:                        Aln
!   Functional form:               Extended Rydberg plus a generalized three-body term
!   Common name:                   ER2+CoxJM
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
module chempotpy_aln_er2coxjm_pes
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

    !> @brief Evaluate pot for the Aln_ER2CoxJM potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_geometry, only: fill_distance_matrix
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp) :: a1, a2, a3, au2, au23, bu2, c(15), c0, c1, c2, c3, c4, c5, c6, co, coij, cojk, coki, dd, de, dx, dy, &
            & dz, ee, fsech, fu2(maxatom, maxatom), oneha, onesi, oneth, poly, q1, q2, q23, q3, q332, r(maxatom, maxatom), re, &
            & ree, rho, rij, rjk, rki, rr, twoth, v, v3, x(maxatom), y(maxatom), z(maxatom)
        integer             :: i, j, k, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        ! two-body
        de  = 1.71013678553981441/autoev
        re  = 5.08182706399609163
        a1  = 1.24074255007327805
        a2  = 0.551880801172447422
        a3  = 0.129970688812896917
        au2 = 0.143243771372740580
        bu2 = 6.50000000000000000/autoang
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
        oneth = sqrt(1.e0_wp/3.e0_wp)
        oneha = sqrt(1.e0_wp/2.e0_wp)
        twoth = sqrt(2.e0_wp/3.e0_wp)
        onesi = sqrt(1.e0_wp/6.e0_wp)
        call fill_distance_matrix(x(1:natom), y(1:natom), z(1:natom), r(1:natom, 1:natom))
        !  calculate the potential energy
        v3 = 0.e0_wp
        do i=1, natom
            do j=i+1, natom
                !  two body term
                rr   = r(i, j)
                coij = 0.e0_wp
                if( rr .ge. bu2 ) cycle
                coij = exp(au23*(1.e0_wp-1.e0_wp/(1.e0_wp-rr/bu2)))
                rij  = (rr-ree)/ree
                do k=j+1, natom
                    ! three body term
                    rr   = r(k, j)
                    cojk = 0.e0_wp
                    if( rr .ge. bu2 ) cycle
                    cojk = exp(au23*(1.e0_wp-1.e0_wp/(1.e0_wp-rr/bu2)))
                    rjk  = (rr-ree)/ree
                    rr   = r(k, i)
                    coki = 0.e0_wp
                    if( rr .ge. bu2 ) cycle
                    coki = exp(au23*(1.e0_wp-1.e0_wp/(1.e0_wp-rr/bu2)))
                    rki  = (rr-ree)/ree
                    q1   = oneth*(rij+rjk+rki)
                    q2   = oneha*(rjk-rki)
                    q3   = twoth*rij-onesi*(rjk+rki)
                    q23  = q2**2+q3**2
                    q332 = q3**3-3.e0_wp*q3*q2**2
                    poly = c0+c1*q1+c2*q1**2+c3*q23+c4*q1**3+c5*q1*q23+c6*q332
                    !      poly=c0+c1*q1+c2*q1**2+c3*q23
                    !      fsech=1.d0/dcosh(xa3*q1)
                    fsech = 1.e0_wp
                    co    = coij*cojk*coki
                    v3    = v3+dd*poly*co*fsech
                enddo
            enddo
        enddo
        v = v+v3
    end subroutine pot
end module chempotpy_aln_er2coxjm_pes

!> @brief Expose Aln_ER2CoxJM through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_er2coxjm_pes, only: evaluate
    use chempotpy_kinds,            only: wp
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
