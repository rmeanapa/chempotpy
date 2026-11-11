!   System:                        Aln
!   Functional form:               Embedded atom
!   Common name:                   reMisFMP
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
!              Y. Mishin, D. Farkas, M. J. Mehl, and D. A. Papaconstantopoulos
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
module chempotpy_aln_remisfmp_pes
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

    !> @brief Evaluate pot for the Aln_reMisFMP potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        implicit none
        double precision :: a, alp1, alp2, b1, b2, delta, dx, dy, dz, e1, e2, ef, em1, em2, f0, f2, h, psiarg, psicut, q1, q2, &
            & q3, r1, r2, r3, r4, rc, rho, rhobar, rr, v, v2, x(maxatom), xksum, y(maxatom), z(maxatom)
        integer                     :: i, j, maxatom, natom
        double precision, parameter :: autoev  = 27.2113961d0
        double precision, parameter :: autoang = 0.529177249d0
        !  calculate the potential energy
        rc    = 6.11700048851978462/autoang
        h     = 4.78798241328773777/autoang
        e1    = 1.38861748900830473/autoev
        e2    = 0.176209086468001944e-02/autoev
        r1    = 2.34000977039570124/autoang
        r2    = 5.13873961895456777/autoang
        alp1  = 1.01670737664875421*autoang
        alp2  = 1.31289692232535415*autoang
        delta = 0.109501709819247672/autoev
        a     = 0.998045920859794850e-01
        r3    = 4.41108939912066411/autoang
        r4    = -6.14313629702002917/autoang
        b1    = 0.574206155349291647*autoang**2
        b2    = 0.261455788959452817*autoang
        f0    = -3.13922813873961903/autoev
        f2    = 4.71763556424035180/autoev
        q1    = -3.62530532486565704/autoev
        q2    = -23.1485100146555922/autoev
        q3    = -20.1148021494870548/autoev
        ! original
        !      rc = 6.78040d0   /autoang
        !      h  = 1.41584d0   /autoang
        !      e1 = 2.65219d2   /autoev
        !      e2 = 7.67238d-3  /autoev
        !      r1 = 1.0688d0    /autoang
        !      r2 = 6.45771d0   /autoang
        !      alp1 = 2.09043d0 *autoang
        !      alp2 = 1.03062d0 *autoang
        !      delta = 0.10432d0/autoev
        !      a  = 5.47491d-2
        !      r3 = 2.74668d0    /autoang
        !      r4 = -6.92647d0   /autoang
        !      b1 = 1.92116d0    *autoang**2
        !      b2 = 0.42589d0    *autoang
        !      f0 = -2.81468d0   /autoev
        !      f2 = 5.57686d0    /autoev
        !      q1 = -6.24655d0   /autoev
        !      q2 = -21.52522d0   /autoev
        !      q3 = -15.30493d0   /autoev
        v = 0.d0
        do i=1, natom
            v2     = 0.d0
            rhobar = 0.d0
            do j=1, natom
                if( i .eq. j ) cycle
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = dsqrt(dx*dx+dy*dy+dz*dz)
                ! 2-body part
                psiarg = (rr-rc)/h
                if( psiarg .ge. 0.d0 )then
                    psicut = 0.d0
                else
                    psicut = psiarg**4/(1.d0+psiarg**4)
                endif
                em1 = dexp(-2.d0*alp1*(rr-r1))-2.d0*dexp(-alp1*(rr-r1))
                em2 = dexp(-2.d0*alp2*(rr-r2))-2.d0*dexp(-alp2*(rr-r2))
                v2  = v2 + (e1*em1+e2*em2+delta)*psicut
                ! rhobar
                rho    = (a*dexp(-b1*(rr-r3)**2)+dexp(-b2*(rr-r4)))*psicut
                rhobar = rhobar + rho  ! rho is 2x because we are looping j>i
            enddo
            ! embedding term
            xksum = q1*(rhobar-1.d0)**3 + q2*(rhobar-1.d0)**4 + q3*(rhobar-1.d0)**5
            ef    = f0 + 0.5d0*f2*(rhobar-1)**2 + xksum
            v     = v + 0.5d0*v2 + ef
        enddo
    end subroutine pot
end module chempotpy_aln_remisfmp_pes

!> @brief Expose Aln_reMisFMP through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_remisfmp_pes, only: evaluate
    use chempotpy_kinds,            only: real64
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
