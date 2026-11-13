!   System:                        Aln
!   Functional form:               Extended Rydberg plus Extended Screening
!   Common name:                   ER2+ES
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
module chempotpy_aln_er2es_pes
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

    !> @brief Evaluate pot for the Aln_ER2ES potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a1, a1b, a2, a2b, a3, a3b, ac, au2, au23, au2b, bu2, bu2b, c(10), co, cob, de, deb, dx, dy, dz, ee, &
            & eeb, fu2(maxatom, maxatom), poly, polyb, r(maxatom, maxatom), re, reb, rho, rhob, rij, ril, rjl, rr, s, term, v, &
            & x(maxatom), xk1, xk2, xk3, y(maxatom), z(maxatom)
        integer             :: i, j, l, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        de   = 1.71013678553981441/autoev
        re   = 5.08182706399609163
        a1   = 1.24074255007327805
        a2   = 0.551880801172447422
        a3   = 0.129970688812896917
        au2  = 0.143243771372740580
        bu2  = 6.50000000000000000/autoang
        xk1  = 2.02364435759648265
        xk2  = 0.444699560332193433
        xk3  = 2.37791890571568132
        deb  = 1.77156814851001476/autoev
        reb  = 4.71494870542256983
        a1b  = 1.12725940400586233
        a2b  = 0.496384953590620404
        a3b  = 0.501460674157303332
        au2b = 0.467684416218856813
        bu2b = 6.5e0_wp/autoang
        au23 = 0.665852467024914546e-01
        do i=1, natom
            r(i, i)   = 0.e0_wp
            fu2(i, i) = 0.e0_wp
            do j=i+1, natom
                dx        = x(i)-x(j)
                dy        = y(i)-y(j)
                dz        = z(i)-z(j)
                r(i, j)   = sqrt(dx*dx+dy*dy+dz*dz)
                r(j, i)   = r(i, j)
                fu2(j, i) = 0.e0_wp
                fu2(i, j) = 0.e0_wp
            enddo
        enddo
        ! do BA term, needs xk1,xk2,xk3
        do i=1, natom! alpha loop
            do j=i+1, natom! beta loop
                s = 0.e0_wp
                do l=1, natom! gamma loop
                    if( l .eq. j.or.l .eq. i ) cycle! gamma != alpha or beta
                    rjl = r(l, j)
                    rij = r(i, j)
                    ril = r(l, i)
                    if( rjl .ge. bu2 ) cycle! cutoff for Rbeta,delta
                    if( rij .ge. bu2 ) cycle! cutoff for Rbeta,delta
                    if( ril .ge. bu2 ) cycle! cutoff for Rbeta,delta
                    ac   = au23*(1e0_wp-1e0_wp/(1e0_wp-rjl/bu2))       ! cutoff for Rbeta,delta
                    ac   = ac+au23*(1e0_wp-1e0_wp/(1e0_wp-rij/bu2))    ! cutoff for Rbeta,delta
                    ac   = ac+au23*(1e0_wp-1e0_wp/(1e0_wp-ril/bu2))    ! cutoff for Rbeta,delta
                    term = xk1*exp(-xk2*((ril+rjl)/rij)**xk3+ac)        ! cutoff ti
                    !          if (natom.eq.3) print *,i,j,l,rjl,rij,ril,term
                    s = s+term                         ! sum
                enddo
                !        e2=dexp(-2d0*s)                       ! equivalent to 1-dtanh(s)
                !        fu2(i,j)=2d0*e2/(1d0+e2)              ! equivalent to 1-dtanh(s)
                fu2(i, j) = tanh(s)
                fu2(j, i) = fu2(i, j)
                !        if (natom.eq.3) print *,i,j,fu2(i,j)
            enddo
            fu2(i, i) = 1.e0_wp                           ! shouldn't get used anywh
        enddo
        ! potential, includes cutoff
        v = 0.e0_wp
        do i=1, natom
            do j=i+1, natom
                dx    = x(i)-x(j)
                dy    = y(i)-y(j)
                dz    = z(i)-z(j)
                rr    = sqrt(dx*dx+dy*dy+dz*dz)
                rho   = rr-re
                rhob  = rr-reb
                poly  = 1.e0_wp+a1*rho+a2*rho**2+a3*rho**3
                polyb = 1.e0_wp+a1b*rhob+a2b*rhob**2+a3b*rhob**3
                ee    = exp(-a1*rho)
                eeb   = exp(-a1b*rhob)
                co    = 0.e0_wp
                cob   = 0.e0_wp
                if( rr .le. bu2 ) co=exp(au2*(1.e0_wp-1.e0_wp/(1.e0_wp-rr/bu2)))
                if( rr .le. bu2b ) cob=exp(au2b*(1.e0_wp-1.e0_wp/(1.e0_wp-rr/bu2b)))
                v = v-de*ee*poly*co +deb*eeb*polyb*cob*fu2(i, j)
            enddo
        enddo
    end subroutine pot
end module chempotpy_aln_er2es_pes

!> @brief Expose Aln_ER2ES through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_er2es_pes, only: evaluate
    use chempotpy_kinds,         only: wp
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
