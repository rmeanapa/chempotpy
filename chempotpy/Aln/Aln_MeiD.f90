!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_meid_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate an aluminum-cluster potential through the Python contract.
    subroutine evaluate( coordinates_angstrom, atom_count, derivative_order, &
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
    end subroutine evaluate
    !> @brief Evaluate pot for the Aln_MeiD potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        !   System:                        Aln
        !   Functional form:               Embedded atom
        !   Common name:                   MeiD
        !   Number of derivatives:         0
        !   Number of bodies:              variable
        !   Number of electronic surfaces: 1
        !   Interface:                     HO-MM-0
        !
        !   Notes:  Many-body aluminum potential energy function.
        !
        !  References: J. Mei and J. W. Davenport, Phys. Rev. B, 46, 21, (1992).
        !
        !  Units:
        !       energies    - hartrees
        !       coordinates - bohrs
        !
        !  v       --- energy of the system (output)
        !  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
        !              of the system (input)
        !  natom   --- number of atoms in the system (input)
        !  maxatom --- maximal number of atoms (input)
        implicit none
        double precision :: ab, alpha, beta, bigef, cc, delta, dr, dx, dy, dz, ec, ef, &
            & gamma, phi, phi0, q, r0, rc, rho, rn, rr, s, sterm, v, v2, x, xm, xx, y, z
        integer                     :: i, j, l, m, maxatom, natom
        double precision, parameter :: autoev  = 27.2113961d0
        double precision, parameter :: autoang = 0.529177249d0
        dimension x(maxatom), y(maxatom), z(maxatom)
        dimension cc(6), s(3)
        !  calculate the potential energy
        ec    = 3.39d0/autoev
        phi0  = 0.1318d0/autoev
        r0    = 2.8638d0/autoang
        alpha = 4.6d0
        beta  = 7.10d0
        gamma = 7.34759d0
        delta = 7.35d0
        cc(1) = 0.64085d0
        cc(2) = -6.83764d0
        cc(3) = 26.75616d0
        cc(4) = -47.16495d0
        cc(5) = 36.18925d0
        cc(6) = -8.60834d0
        s(1)  = 12.d0
        s(2)  = 6.d0
        s(3)  = 24.d0
        rn    = 1.75d0*r0
        rc    = 1.95d0*r0
        v     = 0.d0
        do i=1, natom
            rho = 0.d0
            v2  = 0.d0
            do j=1, natom
                if( i .eq. j ) cycle
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = dsqrt(dx*dx+dy*dy+dz*dz)
                ! cutoff function
                if( rr .le. rn )then
                    q = 1.d0
                else if( rr .ge. rc )then
                    q = 0.d0
                else
                    xx = (rr-rn)/(rc-rn)
                    q  = ((1.d0-xx)**3)*(1.d0+3.d0*xx+6.d0*xx**2)
                endif
                ! pair potential phi
                dr  = rr/r0-1.d0
                phi = -0.5d0*phi0*(1.d0+delta*dr)*dexp(-gamma*dr)
                v2  = v2 + phi*q
                ! embedding potential rho
                do l=0, 5
                    ef  = (cc(l+1)/12.d0)*(r0/rr)**l
                    rho = rho + ef*q
                enddo
            3 enddo
            sterm = 0.d0
            if( rho .ne. 0.d0 )then
                do m=1, 3
                    xm = dsqrt(dble(m))
                    sterm=sterm+s(m)*dexp(-gamma*(xm-1.d0)) *(1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho)) *(rho)** &
                        & (xm* gamma/beta)
                enddo
            endif
            ab    = alpha/beta
            bigef = 0.d0
            if( rho .ne. 0.d0 )then
                bigef = -ec*(1.d0-ab*dlog(rho))*rho**ab + 0.5d0*phi0*sterm
            endif
            v = v + v2 + bigef
        2 enddo
    end subroutine pot
end module chempotpy_aln_meid_pes

!> @brief Preserve the historical Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_meid_pes, only: evaluate
    use chempotpy_kinds,        only: real64
    implicit none
    integer, parameter        :: state_count = 1
    real(real64), intent(in)  :: coordinates_angstrom(10000,3)
    integer, intent(in)       :: atom_count, derivative_order
    real(real64), intent(out) :: energy_ev(state_count)
    real(real64), intent(out) :: gradient_ev_per_angstrom(state_count,atom_count,3)
    real(real64), intent(out) :: coupling_per_angstrom(state_count,state_count,atom_count,3)
    call evaluate(coordinates_angstrom, atom_count, derivative_order, energy_ev, gradient_ev_per_angstrom, &
        & coupling_per_angstrom)
end subroutine pes
