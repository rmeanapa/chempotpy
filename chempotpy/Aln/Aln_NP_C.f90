!   System:                        Aln
!   Functional form:               Embedded atom
!   Common name:                   NP-C
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  The functional
!           form is from Ref. 1.  The parameters were re-optimized in Ref. 2
!           against a data set of energies for aluminum clusters and
!           nanoparticles and bulk data.  Reference 3 provides futher
!           background but is not a required reference for this potential.
!  References: (1) F. H. Streitz, and J. W. Mintmire, Phys. Rev. B, 50, 11996 (1994).
!              (2) A. W. Jasper, N. E. Schultz, and D. G. Truhlar, "Analytic
!              Potential Energy Functions for Simulating Aluminum
!              Nanoparticles," in preparation.
!              (3) A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!              and D. G. Truhlar, "Analytic Potential Energy Functions
!              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximal number of atoms (input)

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_np_c_pes
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

    !> @brief Evaluate pot for the Aln_NP_C potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        implicit none
        double precision :: aa, alpha, bb, beta, cc, dx, dy, dz, ff, p1, p2, phi, rho, rr, rrel, rstar, v, x(maxatom), &
            & y(maxatom), z(maxatom)
        integer :: i, j, maxatom, natom
        double precision, parameter :: autoev  = 27.2113961d0
        double precision, parameter :: autoang = 0.529177249d0
        rstar = 2.85655354324007371 / autoang
        alpha = 2.13479173244937614 * autoang
        beta  = 1.16888101913568532 * autoang
        aa    = 0.716877851705074343 / autoev
        bb    = 0.667088456913949440e-01 / autoev
        cc    = 0.333776366433971738 / autoev
        v     = 0.d0
        ! pairwise term
        phi = 0.d0
        ff  = 0.d0
        do i=1, natom
            rho = 0.d0
            do j=1, natom
                if( j .ne. i )then
                    dx   = x(i)-x(j)
                    dy   = y(i)-y(j)
                    dz   = z(i)-z(j)
                    rr   = dsqrt(dx*dx+dy*dy+dz*dz)
                    rrel = rr-rstar
                    rho  = rho+dexp(-beta*rrel)
                endif
                if( i .gt. j )then
                    p1  = 2.d0*bb*dexp(-0.5d0*beta*rrel)
                    p2  = -cc*(1.d0+alpha*rrel)*dexp(-alpha*rrel)
                    phi = phi+p1+p2
                endif
            enddo
            ff = ff-aa*dsqrt(rho)
        enddo
        v = phi+ff
    end subroutine pot
end module chempotpy_aln_np_c_pes

!> @brief Expose Aln_NP_C through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_np_c_pes, only: evaluate
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
