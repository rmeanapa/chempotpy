!   System:                        Aln
!   Functional form:               Heuristic
!   Common name:                   reErkIV
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  The parameters
!           for this PEF have been parameterized against a data set of
!           aluminum cluster energies and bulk data.
!  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!              and D. G. Truhlar, "Analytic Potentials Energy Functions
!              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
!              S. Erkoc, Phys. Stat. Sol. (b) 161 (1990), 211.
!              S. Erkoc, Empirical Potential Energy Functions
!                 Used in the Simulations of Materials Properties,
!                 Annual Reviews of Computational Physics IX,
!                 D. Stauffer, ed., World Scientific, Singapore 2001,
!                 pp. 1-103.
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximal number of atoms (input)

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_reerkiv_pes
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

    !> @brief Evaluate pot for the Aln_reErkIV potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: aa, alpha, c2, c3, dist(maxatom,maxatom), dx, dy, dz, en, en2, fijk, fikj, fjki, r0, r02, r1, r2, &
            & rr, temp1, temp2, u(maxatom, maxatom), v, vmany, vpair, wijk, x(maxatom), y(maxatom), z(maxatom)
        integer             :: i, j, k, maxatom, natom
        real(wp), parameter :: autoang = 0.5291277249e0_wp
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        alpha = log(2.e0_wp)
        ! Original parameters
        !      en = 2.072844d0
        !      r0 = 2.47d0/autoang
        !      aa = 6.20d0/autoev
        !      c2 = 0.2635582d0
        !      c3 = -0.1766508d0
        ! new parameters
        en    = 2.396042
        r0    = 2.409526/autoang
        aa    = 9.999023/autoev
        c2    = 0.6116756
        c3    = -0.9720078
        vpair = 0.0e0_wp
        vmany = 0.0e0_wp
        r02   = r0**2
        en2   = en*2e0_wp
        !  pair potential, energy is denoted as vpair
        do i=1, natom-1
            do j=i+1, natom
                dx        = x(i)-x(j)
                dy        = y(i)-y(j)
                dz        = z(i)-z(j)
                rr        = sqrt(dx*dx+dy*dy+dz*dz)
                dist(i,j) = rr
                dist(j,i) = rr
                r1        = rr/r0
                r2        = r0/rr
                temp1     = exp(-2.e0_wp*alpha*r1**2)
                temp2     = exp(-1.e0_wp*alpha*r1**2)
                u(i, j)   = aa*((r2**en2)*temp1 - (r2**en)*temp2)
                u(j, i)   = u(i, j)
                vpair     = u(i, j) + vpair
            enddo
        enddo
        if( natom .ge. 3 )then
            do i=1, natom-2
                do j=i+1, natom-1
                    do k=j+1, natom
                        fijk  = exp(-(dist(i,k)**2 + dist(j,k)**2)/r02)
                        fikj  = exp(-(dist(i,j)**2 + dist(k,j)**2)/r02)
                        fjki  = exp(-(dist(j,i)**2 + dist(k,i)**2)/r02)
                        wijk  = (u(i, j)*fijk+u(i, k)*fikj+u(j, k)*fjki)
                        vmany = wijk + vmany
                    enddo
                enddo
            enddo
        endif
        v = c2*vpair + c3*vmany
    end subroutine pot
end module chempotpy_aln_reerkiv_pes

!> @brief Expose Aln_reErkIV through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_reerkiv_pes, only: evaluate
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
