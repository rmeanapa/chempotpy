!   System:                        Aln
!   Functional form:               Embedded atom/Effective medium theory
!   Common name:                   Jac
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  For a recent
!           discussion of this potential, see:
!           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!           and D. G. Truhlar, "Analytic Potentials Energy Functions
!           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
!  References: K. W. Jacobsen, Comments Cond. Mat. Phys. 14,
!                 129 (1988).
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximal number of atoms (input)

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_jac_pes
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

    !> @brief Evaluate pot for the Aln_Jac potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: alpha, beta, bs, dx, dy, dz, e0, e2, e3, en, ep, eta, eta1, eta2, fra, pi, rr, s, u, ud, udno, uno, &
            & v, vv, x(maxatom), y(maxatom), z(maxatom)
        integer             :: i, j, m, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        !  calculate the potential energy
        eta1  = 0.3
        en    = 0.007
        s     = 3.0
        eta   = 2.0
        en    = 0.04724  *autoang**3
        s     = 1.5875   /autoang
        eta   = 3.7795   *autoang
        e0    = -3.28     /autoev
        e2    = 1.12     /autoev
        e3    = -0.35     /autoev
        alpha = 189.7   /(autoev*autoang**3)
        eta1  = 0.5669   *autoang
        pi    = acos(-1.e0_wp)
        beta  = (2.e0_wp**(-0.5))*((16.e0_wp*pi/3.e0_wp)**(1.e0_wp/3.e0_wp))
        eta2  = (eta+eta1)/beta
        ep    = eta/(eta+eta1)
        bs    = beta*s
        fra   = 1.e0_wp/12.e0_wp
        v     = 0e0_wp
        do i=1, natom
            u = 0e0_wp
            do j=1, natom
                if( j .eq. i ) cycle
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = sqrt(dx*dx+dy*dy+dz*dz)
                u  = u+exp(-eta2*(rr-bs))   ! EAS(1) w/o exponent
            enddo
            ud = 0e0_wp
            do m=1, natom
                if( m .eq. i ) cycle
                dx = x(i)-x(m)
                dy = y(i)-y(m)
                dz = z(i)-z(m)
                rr = sqrt(dx*dx+dy*dy+dz*dz)
                ud = ud+exp(-eta*(rr/beta-s))  ! EAS(2)
            enddo
            uno  = (fra*u)**ep   ! EAS(1)
            udno = fra*ud       ! EAS(2)
            vv   = e0+e2*(uno-1e0_wp)**2+e3*(uno-1e0_wp)**3 +alpha*en*(uno-udno)
            v    = v+vv
        enddo
    end subroutine pot
end module chempotpy_aln_jac_pes

!> @brief Expose Aln_Jac through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_jac_pes, only: evaluate
    use chempotpy_kinds,       only: wp
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
