!   System:                        Aln
!   Functional form:               Extended Rydberg plus Extended Coordination Number
!   Common name:                   ER2+ECN
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
module chempotpy_aln_er2ecn_pes
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

    !> @brief Evaluate pot for the Aln_ER2ECN potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a1, a1b, a2, a2b, a3, a3b, abscz1, abscz2, au2, au2b, bu2, bu2b, co, cob, cz1, cz2, de, deb, du2, &
            & dx, dy, dz, ec1, ee, eeb, es1, fu2(maxatom, maxatom), g1, g1e, g2, gcorr(maxatom), gu2, gzero, poly, polyb, &
            & r(maxatom, maxatom), re, reb, rho, rhob, rij, rr, term, v, x(maxatom), y(maxatom), z(maxatom)
        integer             :: i, j, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        de    = 1.71013678553981441/autoev
        re    = 5.08182706399609163
        a1    = 1.24074255007327805
        a2    = 0.551880801172447422
        a3    = 0.129970688812896917
        au2   = 0.143243771372740580
        bu2   = 6.50000000000000000/autoang
        deb   = 1.26365412799218380/autoev
        reb   = 4.51734245236932086
        a1b   = 0.766243282852955510
        a2b   = 0.416023448949682462
        a3b   = 0.260429897410845157
        au2b  = 1.08568637029799708
        bu2b  = 6.5e0_wp/autoang
        g1    = 1.80898876404494380
        g2    = 8.09721543722520742/autoang
        du2   = 2.21812408402540306
        gu2   = 0.853639472398632138
        gzero = 4.62188568637029817
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
        ! CALCULATE CN FUNCTION
        !     first calculate g_alpha
        g1e = exp(g1)
        do i=1, natom! alpha loop
            gcorr(i) = 0e0_wp
            do j=1, natom! beta loop
                if( i .eq. j ) cycle! beta != alpha
                rij = r(i, j)
                if( rij .ge. g2 ) cycle! cutoff term counts nearb
                term     = exp(-g1/(1e0_wp-(rij/g2)))         ! cutoff term counts nearb
                gcorr(i) = gcorr(i)+term                ! sum over beta
            enddo
            gcorr(i) = gcorr(i)*g1e                  ! g_alpha
        enddo
        !     then calculate CN term
        do i=1, natom! alpha loop
            do j=1, natom! beta loop
                if( i .eq. j ) cycle! beat != alpha
                rij  = r(i, j)
                term = 0.e0_wp
                goto_block_32_129: block
                    if( rij .ge. g2 ) exit goto_block_32_129! g_alphabeta
                    term = exp(g1*(1e0_wp-1e0_wp/(1e0_wp-(rij/g2)))) ! g_alphabeta
                end block goto_block_32_129
                cz1    = (gcorr(i)-term)/gzero              ! ratio of coordination fo
                abscz1 = abs(cz1)
                if( abscz1 .le. 10e-12_wp ) cz1=abscz1! eliminates negative valu
                es1    = 1e0_wp/(1e0_wp+cz1**gu2)                 ! compute alpha part of f
                cz2    = (gcorr(j)-term)/gzero              ! ratio of coordination fo
                abscz2 = abs(cz2)
                if( abscz2 .le. 10e-12_wp ) cz2=abscz2! eliminates negative valu
                ec1 = 1e0_wp/(1e0_wp+cz2**gu2)                 ! compute beta part of f f
                !       fu2(i,j)=1d0+du2*(es1*ec1-1d0)         ! put alpha and beta parts together:  1 + f for regular formalism
                fu2(i, j) = du2*(es1*ec1-1e0_wp)         ! put alpha and beta parts tog
                !       print *,es1,ec1,fu2(i,j)
            enddo
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
                v = v-de*ee*poly*co -deb*eeb*polyb*cob*fu2(i, j)
            enddo
        enddo
    end subroutine pot
end module chempotpy_aln_er2ecn_pes

!> @brief Expose Aln_ER2ECN through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_er2ecn_pes, only: evaluate
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
