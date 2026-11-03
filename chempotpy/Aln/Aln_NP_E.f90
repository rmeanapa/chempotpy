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
        write(*,*) 'Only energy is available'
    endif
end subroutine pes

!> @brief Evaluate pot for the Aln_NP_E potential.
subroutine pot( x, y, z, v, natom, maxatom )

    !   System:                        Aln
    !   Functional form:               Embedded atom
    !   Common name:                   NP-E
    !   Number of derivatives:         0
    !   Number of bodies:              variable
    !   Number of electronic surfaces: 1
    !   Interface:                     HO-MM-0
    !
    !   Notes:  Many-body aluminum potential energy function.  The functional
    !           form is from Ref. 1.  The parameters were re-optimized in Ref. 2
    !           against a data set of energies for aluminum clusters and
    !           nanoparticles and bulk data.  Reference 3 provides futher
    !           background but is not a required reference for this potential.
    !
    !  References: (1) Y. Mishin, D. Farkas, M. J. Mehl, and D. A. Papaconstantopoulos
    !              Mat. Res. Soc. Symp. Proc. 583, 535, (1999).
    !              (2) A. W. Jasper, N. E. Schultz, and D. G. Truhlar, "Analytic
    !              Potential Energy Functions for Simulating Aluminum
    !              Nanoparticles," in preparation.
    !              (3) A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
    !              and D. G. Truhlar, "Analytic Potential Energy Functions
    !              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
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

    double precision :: a, alp1, alp2, b1, b2, delta, dx, dy, dz, e1, e2, ef, em1, &
        & em2, f0, f2, h, psiarg, psicut, q1, q2, q3, r1, r2, r3, r4, rc, rho, rhobar, rr, v, v2, x, &
        & xksum, y, z
    integer                     :: i, j, maxatom, natom
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
    dimension x(maxatom),y(maxatom),z(maxatom)

    rc    = 6.21230561520981794/autoang
    h     = 1.23419121131286391/autoang
    e1    = 323.306599351550119/autoev
    e2    = 0.620692841774021006e-02/autoev
    r1    = 0.849576438380276855/autoang
    r2    = 4.63688776461126384/autoang
    alp1  = 2.99843114651358267*autoang
    alp2  = 1.55910921885380915*autoang
    delta = 0.794681217776774529e-01/autoev
    a     = 0.637166953301545880e-01
    r3    = 2.49022329534347442/autoang
    r4    = -6.69553673168315022/autoang
    b1    = 1.42436502858883451*autoang**2
    b2    = 0.457015505777396858*autoang
    f0    = -3.50601555923931452/autoev
    f2    = 7.33764124826237207/autoev
    q1    = -5.98181779608959907/autoev
    q2    = -19.0233990007299987/autoev
    q3    = -12.7919662283307112/autoev

    ! original parameters
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
    do 2 i=1,natom
        v2     = 0.d0
        rhobar = 0.d0
        do 3 j=1,natom
            if (i.eq.j) go to 3
            dx = x(i)-x(j)
            dy = y(i)-y(j)
            dz = z(i)-z(j)
            rr = dsqrt(dx*dx+dy*dy+dz*dz)

            ! 2-body part
            psiarg = (rr-rc)/h
            if (psiarg.ge.0.d0) then
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

        3 continue
        ! embedding term
        xksum = q1*(rhobar-1.d0)**3 + q2*(rhobar-1.d0)**4 + q3*(rhobar-1.d0)**5

        ef = f0 + 0.5d0*f2*(rhobar-1)**2 + xksum
        v  = v + 0.5d0*v2 + ef
    2 continue

end subroutine pot
