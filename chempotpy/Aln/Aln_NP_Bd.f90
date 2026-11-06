!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_np_bd_pes
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
    !> @brief Evaluate pot for the Aln_NP_Bd potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        !   System:                        Aln
        !   Functional form:               Embedded atom
        !   Common name:                   NP-B
        !   Number of derivatives:         1
        !   Number of bodies:              variable
        !   Number of electronic surfaces: 1
        !   Interface:                     HO-MM-1
        !
        !   Notes:  Many-body aluminum potential energy function.  The functional
        !           form is from Ref. 1.  The parameters were re-optimized in Ref. 2
        !           against a data set of energies for aluminum clusters and
        !           nanoparticles and bulk data.  Reference 3 provides futher
        !           background but is not a required reference for this potential.
        !           This coding contains analytic first derivatives.
        !           Note that eq 37 of Reference 3 has two typos related to this
        !           potential. A minus sign is missing in front of phi-sub-0, and
        !           the sign in front of the variable "d" should be "+".
        !
        !  References: (1) J. Mei and J. W. Davenport, Phys. Rev. B, 46, 21, (1992).
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
        !  dvdx,dvdy,dvdz   --- one-dimensional arrays of the gradients with respect
        !                       to the Cartesian coordinates of the system (output)
        !  natom   --- number of atoms in the system (input)
        !  maxatom --- maximum number of atoms (input)
        implicit none
        double precision :: ab, alpha, beta, bigef, cc, d_rho, d_v2, dbigefdrho, delta, &
            & dphidr, dphiqdr, dqdr, dr, drdr, drhodr, dstermdrho, dvdx, dvdy, dvdz, dx, dxxdr, dy, dz, &
            & ec, ef, gamma, phi, phi0, q, r0, rc, rho, rn, rr, s, sterm, v, v2, x, xm, xx, y, z
        integer                     :: i, j, k, l, m, maxatom, natom
        double precision, parameter :: autoev  = 27.2113961d0
        double precision, parameter :: autoang = 0.529177249d0
        dimension x(maxatom), y(maxatom), z(maxatom)
        dimension cc(6), s(3)
        dimension d_v2(3,maxatom), d_rho(3,maxatom), &  ! Derivs
            & dvdx(maxatom), dvdy(maxatom), dvdz(maxatom)
        ! Original Literature Parameters (Ref. 1)
        !        ec    =   3.39d0/autoev
        !        phi0  =   0.1318d0/autoev
        !        r0    =   2.8638d0/autoang
        !        alpha =   4.6d0
        !        beta  =   7.10d0
        !        gamma =   7.34759d0
        !        delta =   7.35d0
        !        cc(1) =   0.64085d0
        !        cc(2) =  -6.83764d0
        !        cc(3) =  26.75616d0
        !        cc(4) = -47.16495d0
        !        cc(5) =  36.18925d0
        !        cc(6) =  -8.60834d0
        !        s(1)  =  12.d0
        !        s(2)  =   6.d0
        !        s(3)  =  24.d0
        !        rn    =   1.75d0*r0
        !        rc    =   1.95d0*r0
        ! Cluster-parameterized Parameters (Ref. 3)
        !        ec    =   3.648900830d0/autoev
        !        phi0  =   0.395759648d0/autoev
        !        r0    =   2.624670249d0/autoang
        !        alpha =   5.194382022d0
        !        beta  =   4.748803127d0
        !        gamma =   5.754323400d0
        !        delta =   7.663800684d0
        !        cc(1) =   0.010571568d0
        !        cc(2) =  -7.999267220d0
        !        cc(3) =  31.689301416d0
        !        cc(4) = -43.634587200d0
        !        cc(5) =  28.121641426d0
        !        cc(6) =  -6.751099169d0
        !        s(1)  =   8.657547631d0
        !        s(2)  =   4.169614069d0
        !        s(3)  =  27.910600879d0
        !        rn    =   1.75d0*r0
        !        rc    =   1.95d0*r0
        ! Nanoparticle-parameterized Parameters (Ref. 2)
        ec    = 2.833661628d0/autoev
        phi0  = 0.209474568d0/autoev
        r0    = 2.759835989d0/autoang
        alpha = 4.953631991d0
        beta  = 5.202672172d0
        gamma = 5.824302949d0
        delta = 8.968682037d0
        cc(1) = 0.433294196d0
        cc(2) = -7.305279256d0
        cc(3) = 29.818956621d0
        cc(4) = -54.437991632d0
        cc(5) = 48.412067298d0
        cc(6) = -15.525225110d0
        s(1)  = 6.927645227d0
        s(2)  = 3.861172975d0
        s(3)  = 15.498062621d0
        rn    = 1.75d0*r0
        rc    = 1.95d0*r0
        ! Initialize
        v = 0.d0
        do i=1, natom
            dvdx(i) = 0.d0
            dvdy(i) = 0.d0
            dvdz(i) = 0.d0
        enddo
        ! Double loop (main loop)
        do i=1, natom
            rho = 0.d0
            v2  = 0.d0
            do k=1, 3
                do j=1, natom
                    d_rho(k,j) = 0.d0
                    d_v2(k,j)  = 0.d0
                enddo
            enddo
            do j=1, natom
                if( i .eq. j ) cycle
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = dsqrt(dx*dx+dy*dy+dz*dz)
                ! Cutoff function
                if( rr .le. rn )then
                    q    = 1.d0
                    dqdr = 0.d0  ! Derivs
                else if( rr .ge. rc )then
                    cycle
                    q    = 0.d0
                    dqdr = 0.d0  ! Derivs
                else
                    xx    = (rr-rn)/(rc-rn)
                    q     = ((1.d0-xx)**3)*(1.d0+3.d0*xx+6.d0*xx**2)
                    dxxdr = 1.d0/(rc-rn)  ! Derivs
                    dqdr  = ((1.d0-xx)**3)*(3.d0+12.d0*xx)*dxxdr &  ! Derivs
                        & -(3.d0*(1.d0-xx)**2)*(1.d0+3.d0*xx+6.d0*xx**2)*dxxdr
                endif
                ! Pair potential phi
                dr   = rr/r0-1.d0
                phi  = -0.5d0*phi0*(1.d0+delta*dr)*dexp(-gamma*dr)
                v2   = v2 + phi*q
                drdr = 1.d0/r0  ! Derivs
                dphidr  = -0.5d0*phi0*(1.d0+delta*dr)*(-gamma*drdr)*dexp(-gamma*dr) -0.5d0*phi0*delta*drdr*dexp(- &
                    & gamma*dr)
                dphiqdr = dqdr*phi+dphidr*q  ! Derivs
                ! Embedding potential rho
                drhodr = 0.d0  ! Derivs
                do l=0, 5
                    ef  = (cc(l+1)/12.d0)*(r0/rr)**l
                    rho = rho + ef*q
                    if( l .ne. 0 ) drhodr = drhodr &! Derivs
                        & - (q/rr)*dble(l)*(cc(l+1)/12.d0)*(r0/rr)**l+ef*dqdr
                    if( l .eq. 0 ) drhodr = drhodr + ef*dqdr! Derivs
                enddo
                ! v2 (pairwise repulsion) derivs
                d_v2(1,i) = d_v2(1,i) + dphiqdr*dx/rr
                d_v2(1,j) = d_v2(1,j) - dphiqdr*dx/rr
                d_v2(2,i) = d_v2(2,i) + dphiqdr*dy/rr
                d_v2(2,j) = d_v2(2,j) - dphiqdr*dy/rr
                d_v2(3,i) = d_v2(3,i) + dphiqdr*dz/rr
                d_v2(3,j) = d_v2(3,j) - dphiqdr*dz/rr
                ! rho (density) derivs
                d_rho(1,i) = d_rho(1,i) + drhodr*dx/rr
                d_rho(1,j) = d_rho(1,j) - drhodr*dx/rr
                d_rho(2,i) = d_rho(2,i) + drhodr*dy/rr
                d_rho(2,j) = d_rho(2,j) - drhodr*dy/rr
                d_rho(3,i) = d_rho(3,i) + drhodr*dz/rr
                d_rho(3,j) = d_rho(3,j) - drhodr*dz/rr
            3 enddo
            sterm      = 0.d0
            dstermdrho = 0.d0  ! Derivs
            if( rho .gt. 0.d0 )then
                do m=1, 3
                    xm = dsqrt(dble(m))
                    sterm=sterm+s(m)*dexp(-gamma*(xm-1.d0)) *(1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho)) *(rho)** &
                        & (xm* gamma/beta)
                    dstermdrho = dstermdrho + s(m)*dexp(-gamma*(xm-1.d0))*( &  ! Der
                        & (1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho)) &
                        & *(xm*gamma/beta)*(rho)**(xm*gamma/beta-1.d0) &
                        & +(-xm*(delta/beta)*(1.d0/rho))*(rho)**(xm*gamma/beta) )
                enddo
            endif
            ab         = alpha/beta
            bigef      = 0.d0
            dbigefdrho = 0.d0
            if( rho .gt. 0.d0 )then
                bigef      = -ec*(1.d0-ab*dlog(rho))*rho**ab + 0.5d0*phi0*sterm
                dbigefdrho = -ec*(1.d0-ab*dlog(rho))*(ab)*rho**(ab-1.d0) &  ! Deri
                    & -ec*(-ab*(1.d0/rho))*rho**ab &
                    & +0.5d0*phi0*dstermdrho
            endif
            ! Combine derivs
            do j=1, natom
                dvdx(j) = dvdx(j) + d_v2(1,j) + dbigefdrho*d_rho(1,j)
                dvdy(j) = dvdy(j) + d_v2(2,j) + dbigefdrho*d_rho(2,j)
                dvdz(j) = dvdz(j) + d_v2(3,j) + dbigefdrho*d_rho(3,j)
            enddo
            v = v + v2 + bigef
        2 enddo
    end subroutine pot
end module chempotpy_aln_np_bd_pes

!> @brief Preserve the historical Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_np_bd_pes, only: evaluate
    use chempotpy_kinds,         only: real64
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
