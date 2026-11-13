!   System:                        Aln
!   Functional form:               Embedded atom
!   Common name:                   MeiD
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.
!  References: J. Mei and J. W. Davenport, Phys. Rev. B, 46, 21, (1992).
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximal number of atoms (input)

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_meid_pes
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

    !> @brief Evaluate pot for the Aln_MeiD potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: ab, alpha, beta, bigef, cc(6), delta, dr, dx, dy, dz, ec, ef, gamma, phi, phi0, q, r0, rc, rho, rn, &
            & rr, s(3), sterm, v, v2, x(maxatom), xm, xx, y(maxatom), z(maxatom)
        integer             :: i, j, l, m, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        !  calculate the potential energy
        ec    = 3.39e0_wp/autoev
        phi0  = 0.1318e0_wp/autoev
        r0    = 2.8638e0_wp/autoang
        alpha = 4.6e0_wp
        beta  = 7.10e0_wp
        gamma = 7.34759e0_wp
        delta = 7.35e0_wp
        cc(1) = 0.64085e0_wp
        cc(2) = -6.83764e0_wp
        cc(3) = 26.75616e0_wp
        cc(4) = -47.16495e0_wp
        cc(5) = 36.18925e0_wp
        cc(6) = -8.60834e0_wp
        s(1)  = 12.e0_wp
        s(2)  = 6.e0_wp
        s(3)  = 24.e0_wp
        rn    = 1.75e0_wp*r0
        rc    = 1.95e0_wp*r0
        v     = 0.e0_wp
        do i=1, natom
            rho = 0.e0_wp
            v2  = 0.e0_wp
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
                ! pair potential phi
                dr  = rr/r0-1.e0_wp
                phi = -0.5e0_wp*phi0*(1.e0_wp+delta*dr)*exp(-gamma*dr)
                v2  = v2 + phi*q
                ! embedding potential rho
                do l=0, 5
                    ef  = (cc(l+1)/12.e0_wp)*(r0/rr)**l
                    rho = rho + ef*q
                enddo
            enddo
            sterm = 0.e0_wp
            if( rho .ne. 0.e0_wp )then
                do m=1, 3
                    xm = sqrt(real(m, wp))
                    sterm=sterm+s(m)*exp(-gamma*(xm-1.e0_wp)) *(1.e0_wp+delta*(xm-1.e0_wp)-xm*(delta/beta)*log(rho)) * &
                        & (rho)** (xm* gamma/beta)
                enddo
            endif
            ab    = alpha/beta
            bigef = 0.e0_wp
            if( rho .ne. 0.e0_wp )then
                bigef = -ec*(1.e0_wp-ab*log(rho))*rho**ab + 0.5e0_wp*phi0*sterm
            endif
            v = v + v2 + bigef
        enddo
    end subroutine pot
end module chempotpy_aln_meid_pes

!> @brief Expose Aln_MeiD through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_meid_pes, only: evaluate
    use chempotpy_kinds,        only: wp
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
