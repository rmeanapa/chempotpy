!   System:                        Aln
!   Functional form:               Extended Rydberg plus Extended Screening and Coordination Number
!   Common name:                   NP-A
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  The functional
!           form is from Ref. 2.  The parameters were re-optimized in Ref. 1
!           against a data set of energies for aluminum clusters and
!           nanoparticles and bulk data.  This PEF has a cutoff at 6.5 A.
!  References: (1) A. W. Jasper, N. E. Schultz, and D. G. Truhlar, "Analytic
!              Potential Energy Functions for Simulating Aluminum
!              Nanoparticles," in preparation.
!              (2) A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!              and D. G. Truhlar, "Analytic Potential Energy Functions
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
module chempotpy_aln_np_a_pes
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

    !> @brief Evaluate pot for the Aln_NP_A potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a1, a1b, a2, a2b, a3, a3b, abscz1, abscz2, ac, au2, au23, au2b, bu2, bu2b, co, cob, cz1, cz2, de, &
            & deb, du2, dx, dy, dz, ec1, ee, eeb, es1, fcn(maxatom, maxatom), fs(maxatom, maxatom), g1, g1e, g2, gcorr(maxatom), &
            & gu2, gzero, poly, polyb, r(maxatom, maxatom), re, reb, rho, rhob, rij, ril, rjl, rr, s, term, tmp, tmpb, v, &
            & x(maxatom), xk1, xk2, xk3, y(maxatom), z(maxatom)
        integer             :: i, j, l, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        de    = 1.71013678553981441/autoev
        re    = 5.08182706399609163
        a1    = 1.24074255007327805
        a2    = 0.551880801172447422
        a3    = 0.129970688812896917
        au2   = 0.143243771372740580
        bu2   = 6.50000000000000000/autoang
        xk1   = 4.24002677622442103
        xk2   = 0.117656503960960862
        xk3   = 4.78063179546451522
        au23  = 1.63973192904916298
        g1    = 0.708483373073205747
        du2   = 1.13286279334603357
        gu2   = 0.663930057862113232
        gzero = 8.54498572971970027
        g2    = 5.39584023677170066/autoang
        deb   = 1.42526928794948882/autoev
        reb   = 4.87735706664722812
        a1b   = 1.20666644170640880
        a2b   = 0.728296669115275908
        a3b   = 0.215461507389864804
        au2b  = 0.138211749991007299
        bu2b  = 6.50000000000000000/autoang
        do i=1, natom
            r(i, i)   = 0.e0_wp
            fs(i, i)  = 1.e0_wp
            fcn(i, i) = 1.e0_wp
            do j=i+1, natom
                dx        = x(i)-x(j)
                dy        = y(i)-y(j)
                dz        = z(i)-z(j)
                r(i, j)   = sqrt(dx*dx+dy*dy+dz*dz)
                r(j, i)   = r(i, j)
                fs(j, i)  = 1.e0_wp
                fs(i, j)  = 1.e0_wp
                fcn(j, i) = 1.e0_wp
                fcn(i, j) = 1.e0_wp
            enddo
        enddo
        ! do S term, needs xk1,xk2,xk3
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
                    s    = s+term                         ! sum
                enddo
                !        fs(i,j)=1.d0-dtanh(s)
                fs(i, j) = tanh(s)
                fs(j, i) = fs(i, j)
            enddo
            fs(i, i) = 1.e0_wp                           ! shouldn't get used anywhe
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
                goto_block_32_163: block
                    if( rij .ge. g2 ) exit goto_block_32_163! g_alphabeta
                    term = exp(g1*(1e0_wp-1e0_wp/(1e0_wp-(rij/g2)))) ! g_alphabeta
                end block goto_block_32_163
                cz1    = (gcorr(i)-term)/gzero              ! ratio of coordination fo
                abscz1 = abs(cz1)
                if( abscz1 .le. 10e-12_wp ) cz1=abscz1! eliminates negative valu
                es1    = 1e0_wp/(1e0_wp+cz1**gu2)                 ! compute alpha part of f
                cz2    = (gcorr(j)-term)/gzero              ! ratio of coordination fo
                abscz2 = abs(cz2)
                if( abscz2 .le. 10e-12_wp ) cz2=abscz2! eliminates negative valu
                ec1 = 1e0_wp/(1e0_wp+cz2**gu2)                 ! compute beta part of f f
                !       fcn(i,j)=1d0+du2*(es1*ec1-1d0)         ! put alpha and beta parts together, mult by BA term for u2 part
                fcn(i, j) = du2*(es1*ec1-1e0_wp)         ! put alpha and beta parts tog
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
                tmp  = de*ee*poly*co
                tmpb = deb*eeb*polyb*cob
                v    = v-tmp+tmpb*(fs(i, j)-fcn(i, j))
            enddo
        enddo
    end subroutine pot
end module chempotpy_aln_np_a_pes

!> @brief Expose Aln_NP_A through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_np_a_pes, only: evaluate
    use chempotpy_kinds,        only: wp
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
