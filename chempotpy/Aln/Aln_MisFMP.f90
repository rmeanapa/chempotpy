!   System:                        Aln
!   Functional form:               Embedded atom
!   Common name:                   MisFMP
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  For a recent
!           discussion of this potential, see:
!           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!           and D. G. Truhlar, "Analytic Potentials Energy Functions
!           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
!  References: Y. Mishin, D. Farkas, M. J. Mehl, and D. A. Papaconstantopoulos
!              Mat. Res. Soc. Symp. Proc. 583, 535, (1999).
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximal number of atoms (input)

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_misfmp_pes
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

    !> @brief Evaluate pot for the Aln_MisFMP potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a, alp1, alp2, b1, b2, delta, dx, dy, dz, e1, e2, ef, em1, em2, f0, f2, h, psiarg, psicut, q1, q2, &
            & q3, r1, r2, r3, r4, rc, rho, rhobar, rr, v, v2, x(maxatom), xksum, y(maxatom), z(maxatom)
        integer             :: i, j, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        !  calculate the potential energy
        rc    = 6.78040e0_wp   /autoang
        h     = 1.41584e0_wp   /autoang
        e1    = 2.65219e2_wp   /autoev
        e2    = 7.67238e-3_wp  /autoev
        r1    = 1.0688e0_wp    /autoang
        r2    = 6.45771e0_wp   /autoang
        alp1  = 2.09043e0_wp *autoang
        alp2  = 1.03062e0_wp *autoang
        delta = 0.10432e0_wp/autoev
        a     = 5.47491e-2_wp
        r3    = 2.74668e0_wp    /autoang
        r4    = -6.92647e0_wp   /autoang
        b1    = 1.92116e0_wp    *autoang**2
        b2    = 0.42589e0_wp    *autoang
        f0    = -2.81468e0_wp   /autoev
        f2    = 5.57686e0_wp    /autoev
        q1    = -6.24655e0_wp   /autoev
        q2    = -21.52522e0_wp   /autoev
        q3    = -15.30493e0_wp   /autoev
        v     = 0.e0_wp
        do i=1, natom
            v2     = 0.e0_wp
            rhobar = 0.e0_wp
            do j=1, natom
                if( i .eq. j ) cycle
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = sqrt(dx*dx+dy*dy+dz*dz)
                ! 2-body part
                psiarg = (rr-rc)/h
                if( psiarg .ge. 0.e0_wp )then
                    psicut = 0.e0_wp
                else
                    psicut = psiarg**4/(1.e0_wp+psiarg**4)
                endif
                em1 = exp(-2.e0_wp*alp1*(rr-r1))-2.e0_wp*exp(-alp1*(rr-r1))
                em2 = exp(-2.e0_wp*alp2*(rr-r2))-2.e0_wp*exp(-alp2*(rr-r2))
                v2  = v2 + (e1*em1+e2*em2+delta)*psicut
                ! rhobar
                rho    = (a*exp(-b1*(rr-r3)**2)+exp(-b2*(rr-r4)))*psicut
                rhobar = rhobar + rho  ! rho is 2x because we are looping j>i
            enddo
            ! embedding term
            xksum = q1*(rhobar-1.e0_wp)**3 + q2*(rhobar-1.e0_wp)**4 + q3*(rhobar-1.e0_wp)**5
            ef    = f0 + 0.5e0_wp*f2*(rhobar-1)**2 + xksum
            v     = v + 0.5e0_wp*v2 + ef
        enddo
    end subroutine pot
end module chempotpy_aln_misfmp_pes

!> @brief Expose Aln_MisFMP through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_misfmp_pes, only: evaluate
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
