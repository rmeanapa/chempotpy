!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate an aluminum-cluster potential through the Python contract.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
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
        coordinates_angstrom(1:atom_count,:), coordinates_bohr )
    call pot( coordinates_bohr(:,1), coordinates_bohr(:,2), &
        coordinates_bohr(:,3), energy_hartree, atom_count, atom_count )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*, *) 'Only energy is available'
    end if
end subroutine pes

!> @brief Evaluate pot for the Aln_Jac potential.
subroutine pot( x, y, z, v, natom, maxatom )
    !   System:                        Aln
    !   Functional form:               Embedded atom/Effective medium theory
    !   Common name:                   Jac
    !   Number of derivatives:         0
    !   Number of bodies:              variable
    !   Number of electronic surfaces: 1
    !   Interface:                     HO-MM-0
    !
    !   Notes:  Many-body aluminum potential energy function.  For a recent
    !           discussion of this potential, see:
    !           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
    !           and D. G. Truhlar, "Analytic Potentials Energy Functions
    !           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
    !
    !  References: K. W. Jacobsen, Comments Cond. Mat. Phys. 14,
    !                 129 (1988).
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
    double precision :: alpha, beta, bs, dx, dy, dz, e0, e2, e3, en, ep, eta, eta1, &
        & eta2, fra, pi, rr, s, u, ud, udno, uno, v, vv, x, y, z
    integer :: i, j, m, maxatom, natom
    dimension x(maxatom), y(maxatom), z(maxatom)
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
    !  calculate the potential energy
    eta1 = 0.3
    en   = 0.007
    s    = 3.0
    eta  = 2.0

    en    = 0.04724  *autoang**3
    s     = 1.5875   /autoang
    eta   = 3.7795   *autoang
    e0    = -3.28     /autoev
    e2    = 1.12     /autoev
    e3    = -0.35     /autoev
    alpha = 189.7   /(autoev*autoang**3)
    eta1  = 0.5669   *autoang

    pi   = dacos(-1.d0)
    beta = (2.d0**(-0.5))*((16.d0*pi/3.d0)**(1.d0/3.d0))
    !      beta=1.81
    eta2 = (eta+eta1)/beta
    ep   = eta/(eta+eta1)
    bs   = beta*s
    fra  = 1.d0/12.d0

    v = 0d0
    do 2 i=1, natom
        u = 0d0
        do 3 j=1, natom
            if(j.eq.i) go to 3
            dx = x(i)-x(j)
            dy = y(i)-y(j)
            dz = z(i)-z(j)
            rr = dsqrt(dx*dx+dy*dy+dz*dz)
            u  = u+dexp(-eta2*(rr-bs))   ! EAS(1) w/o exponent
        3 continue
        ud = 0d0
        do 4 m=1, natom
            if(m.eq.i) go to 4
            dx = x(i)-x(m)
            dy = y(i)-y(m)
            dz = z(i)-z(m)
            rr = dsqrt(dx*dx+dy*dy+dz*dz)
            ud = ud+dexp(-eta*(rr/beta-s))  ! EAS(2)
        4 continue
        uno  = (fra*u)**ep   ! EAS(1)
        udno = fra*ud       ! EAS(2)
        vv   = e0+e2*(uno-1d0)**2+e3*(uno-1d0)**3 +alpha*en*(uno-udno)
    2 v=v+vv
end subroutine pot
