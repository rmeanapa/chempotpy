!   System:                        Aln
!   Functional form:               Extended Rydberg plus Screening and Coordination Number
!   Common name:                   ER2+SCNm
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
module chempotpy_aln_er2scnm_pes
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

    !> @brief Evaluate pot for the Aln_ER2SCNm potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        implicit none
        double precision :: a1, a2, a3, abscz1, abscz2, ac, au2, au23, bu2, co, cz1, cz2, de, du2, dx, dy, dz, ec1, ee, es1, &
            & fcn(maxatom,maxatom), fs(maxatom,maxatom), g1, g1e, g2, gcorr(maxatom), gu2, gzero, poly, r(maxatom,maxatom), re, &
            & rho, rij, ril, rjl, rr, s, term, v, x(maxatom), xk1, xk2, xk3, y(maxatom), z(maxatom)
        integer :: i, j, l, maxatom, natom
        double precision, parameter :: autoev  = 27.2113961d0
        double precision, parameter :: autoang = 0.529177249d0
        de  = 1.71013678553981441/autoev
        re  = 5.08182706399609163
        a1  = 1.24074255007327805
        a2  = 0.551880801172447422
        a3  = 0.129970688812896917
        au2 = 0.143243771372740580
        bu2 = 6.50000000000000000/autoang
        do i=1, natom
            r(i,i)   = 0.d0
            fs(i,i)  = 1.d0
            fcn(i,i) = 1.d0
            do j=i+1, natom
                dx       = x(i)-x(j)
                dy       = y(i)-y(j)
                dz       = z(i)-z(j)
                r(i,j)   = dsqrt(dx*dx+dy*dy+dz*dz)
                r(j,i)   = r(i,j)
                fs(j,i)  = 1.d0
                fs(i,j)  = 1.d0
                fcn(j,i) = 1.d0
                fcn(i,j) = 1.d0
            enddo
        enddo
        xk1   = 3.84941377625793857
        xk2   = 0.242501221299462610
        xk3   = 3.78133854421104054
        au23  = 0.615705911089399094
        g1    = 0.430117244748412286
        du2   = 0.749975574010747437
        gu2   = 0.833229115779189122
        gzero = 10.4987787005373718
        g2    = 6.49990229604298975/autoang
        ! do S term, needs xk1,xk2,xk3
        do i=1, natom! alpha loop
            do j=i+1, natom! beta loop
                s = 0.d0
                do l=1, natom! gamma loop
                    if( l .eq. j.or.l .eq. i ) cycle! gamma != alpha or beta
                    rjl = r(l,j)
                    rij = r(i,j)
                    ril = r(l,i)
                    if( rjl .ge. bu2 ) cycle! cutoff for Rbeta,delta
                    if( rij .ge. bu2 ) cycle! cutoff for Rbeta,delta
                    if( ril .ge. bu2 ) cycle! cutoff for Rbeta,delta
                    ac   = au23*(1d0-1d0/(1d0-rjl/bu2))       ! cutoff for Rbeta,delta
                    ac   = ac+au23*(1d0-1d0/(1d0-rij/bu2))    ! cutoff for Rbeta,delta
                    ac   = ac+au23*(1d0-1d0/(1d0-ril/bu2))    ! cutoff for Rbeta,delta
                    term = xk1*dexp(-xk2*((ril+rjl)/rij)**xk3+ac)        ! cutoff ti
                    s    = s+term                         ! sum
                enddo
                fs(i,j) = 1.d0-dtanh(s)
                fs(j,i) = fs(i,j)
            enddo
            fs(i,i) = 1.d0                           ! shouldn't get used anywhe
        enddo
        ! CALCULATE CN FUNCTION
        !     first calculate g_alpha
        g1e = dexp(g1)
        do i=1, natom! alpha loop
            gcorr(i) = 0d0
            do j=1, natom! beta loop
                if( i .eq. j ) cycle! beta != alpha
                rij = r(i,j)
                if( rij .ge. g2 ) cycle! cutoff term counts nearb
                term     = dexp(-g1/(1d0-(rij/g2)))         ! cutoff term counts nearb
                gcorr(i) = gcorr(i)+term                ! sum over beta
            enddo
            gcorr(i) = gcorr(i)*g1e                  ! g_alpha
        enddo
        !     then calculate CN term
        do i=1, natom! alpha loop
            do j=1, natom! beta loop
                if( i .eq. j ) cycle! beat != alpha
                rij  = r(i,j)
                term = 0.d0
                goto_block_32_152: block
                    if( rij .ge. g2 ) exit goto_block_32_152! g_alphabeta
                    term = dexp(g1*(1d0-1d0/(1d0-(rij/g2)))) ! g_alphabeta
                end block goto_block_32_152
                cz1    = (gcorr(i)-term)/gzero              ! ratio of coordination fo
                abscz1 = dabs(cz1)
                if( abscz1 .le. 10d-12 ) cz1=abscz1! eliminates negative valu
                es1    = 1d0/(1d0+cz1**gu2)                 ! compute alpha part of f
                cz2    = (gcorr(j)-term)/gzero              ! ratio of coordination fo
                abscz2 = dabs(cz2)
                if( abscz2 .le. 10d-12 ) cz2=abscz2! eliminates negative valu
                ec1      = 1d0/(1d0+cz2**gu2)                 ! compute beta part of f f
                fcn(i,j) = 1d0+du2*(es1*ec1-1d0)         ! put alpha and beta parts
            enddo
        enddo
        ! potential, includes cutoff
        v = 0.d0
        do i=1, natom
            do j=i+1, natom
                dx   = x(i)-x(j)
                dy   = y(i)-y(j)
                dz   = z(i)-z(j)
                rr   = dsqrt(dx*dx+dy*dy+dz*dz)
                rho  = rr-re
                poly = 1.d0+a1*rho+a2*rho**2+a3*rho**3
                ee   = dexp(-a1*rho)
                co   = 0.d0
                if( rr .le. bu2 ) co=dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
                v = v-de*ee*poly*co*fs(i,j)*fcn(i,j)
            enddo
        enddo
    end subroutine pot
end module chempotpy_aln_er2scnm_pes

!> @brief Expose Aln_ER2SCNm through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_er2scnm_pes, only: evaluate
    use chempotpy_kinds,           only: real64
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
