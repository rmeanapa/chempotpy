!   System:                        Aln
!   Functional form:               Embedded atom
!   Common name:                   NP-B
!   Number of derivatives:         1
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-1
!   Notes:  Many-body aluminum potential energy function.  The functional
!           form is from Ref. 1.  The parameters were re-optimized in Ref. 2
!           against a data set of energies for aluminum clusters and
!           nanoparticles and bulk data.  Reference 3 provides futher
!           background but is not a required reference for this potential.
!           This coding contains analytic first derivatives.
!           Note that eq 37 of Reference 3 has two typos related to this
!           potential. A minus sign is missing in front of phi-sub-0, and
!           the sign in front of the variable "d" should be "+".
!  References: (1) J. Mei and J. W. Davenport, Phys. Rev. B, 46, 21, (1992).
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
!  dvdx,dvdy,dvdz   --- one-dimensional arrays of the gradients with respect
!                       to the Cartesian coordinates of the system (output)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximum number of atoms (input)

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_np_bd_pes
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

    !> @brief Evaluate pot for the Aln_NP_Bd potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: ab, alpha, beta, bigef, cc(6), d_rho(3,maxatom), d_v2(3,maxatom), dbigefdrho, delta, dphidr, &
            & dphiqdr, dqdr, dr, drdr, drhodr, dstermdrho, dvdx(maxatom), dvdy(maxatom), dvdz(maxatom), dx, dxxdr, dy, dz, ec, &
            & ef, gamma, phi, phi0, q, r0, rc, rho, rn, rr, s(3), sterm, v, v2, x(maxatom), xm, xx, y(maxatom), z(maxatom)
        integer             :: i, j, k, l, m, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
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
        ec    = 2.833661628e0_wp/autoev
        phi0  = 0.209474568e0_wp/autoev
        r0    = 2.759835989e0_wp/autoang
        alpha = 4.953631991e0_wp
        beta  = 5.202672172e0_wp
        gamma = 5.824302949e0_wp
        delta = 8.968682037e0_wp
        cc(1) = 0.433294196e0_wp
        cc(2) = -7.305279256e0_wp
        cc(3) = 29.818956621e0_wp
        cc(4) = -54.437991632e0_wp
        cc(5) = 48.412067298e0_wp
        cc(6) = -15.525225110e0_wp
        s(1)  = 6.927645227e0_wp
        s(2)  = 3.861172975e0_wp
        s(3)  = 15.498062621e0_wp
        rn    = 1.75e0_wp*r0
        rc    = 1.95e0_wp*r0
        ! Initialize
        v = 0.e0_wp
        do i=1, natom
            dvdx(i) = 0.e0_wp
            dvdy(i) = 0.e0_wp
            dvdz(i) = 0.e0_wp
        enddo
        ! Double loop (main loop)
        do i=1, natom
            rho = 0.e0_wp
            v2  = 0.e0_wp
            do k=1, 3
                do j=1, natom
                    d_rho(k,j) = 0.e0_wp
                    d_v2(k,j)  = 0.e0_wp
                enddo
            enddo
            do j=1, natom
                if( i .eq. j ) cycle
                dx = x(i)-x(j)
                dy = y(i)-y(j)
                dz = z(i)-z(j)
                rr = sqrt(dx*dx+dy*dy+dz*dz)
                ! Cutoff function
                if( rr .le. rn )then
                    q    = 1.e0_wp
                    dqdr = 0.e0_wp  ! Derivs
                elseif( rr .ge. rc )then
                    cycle
                    q    = 0.e0_wp
                    dqdr = 0.e0_wp  ! Derivs
                else
                    xx    = (rr-rn)/(rc-rn)
                    q     = ((1.e0_wp-xx)**3)*(1.e0_wp+3.e0_wp*xx+6.e0_wp*xx**2)
                    dxxdr = 1.e0_wp/(rc-rn)  ! Derivs
                    dqdr  = ((1.e0_wp-xx)**3)*(3.e0_wp+12.e0_wp*xx)*dxxdr &  ! Derivs
                        & -(3.e0_wp*(1.e0_wp-xx)**2)*(1.e0_wp+3.e0_wp*xx+6.e0_wp*xx**2)*dxxdr
                endif
                ! Pair potential phi
                dr   = rr/r0-1.e0_wp
                phi  = -0.5e0_wp*phi0*(1.e0_wp+delta*dr)*exp(-gamma*dr)
                v2   = v2 + phi*q
                drdr = 1.e0_wp/r0  ! Derivs
                dphidr  = -0.5e0_wp*phi0*(1.e0_wp+delta*dr)*(-gamma*drdr)*exp(-gamma*dr) -0.5e0_wp*phi0*delta*drdr* &
                    & exp(- gamma*dr)
                dphiqdr = dqdr*phi+dphidr*q  ! Derivs
                ! Embedding potential rho
                drhodr = 0.e0_wp  ! Derivs
                do l=0, 5
                    ef  = (cc(l+1)/12.e0_wp)*(r0/rr)**l
                    rho = rho + ef*q
                    if( l .ne. 0 ) drhodr = drhodr &! Derivs
                        & - (q/rr)*real(l, wp)*(cc(l+1)/12.e0_wp)*(r0/rr)**l+ef*dqdr
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
            enddo
            sterm      = 0.e0_wp
            dstermdrho = 0.e0_wp  ! Derivs
            if( rho .gt. 0.e0_wp )then
                do m=1, 3
                    xm = sqrt(real(m, wp))
                    sterm=sterm+s(m)*exp(-gamma*(xm-1.e0_wp)) *(1.e0_wp+delta*(xm-1.e0_wp)-xm*(delta/beta)*log(rho)) * &
                        & (rho)** (xm* gamma/beta)
                    dstermdrho = dstermdrho + s(m)*exp(-gamma*(xm-1.e0_wp))*( &  ! Der
                        & (1.e0_wp+delta*(xm-1.e0_wp)-xm*(delta/beta)*log(rho)) &
                        & *(xm*gamma/beta)*(rho)**(xm*gamma/beta-1.e0_wp) &
                        & +(-xm*(delta/beta)*(1.e0_wp/rho))*(rho)**(xm*gamma/beta) )
                enddo
            endif
            ab         = alpha/beta
            bigef      = 0.e0_wp
            dbigefdrho = 0.e0_wp
            if( rho .gt. 0.e0_wp )then
                bigef      = -ec*(1.e0_wp-ab*log(rho))*rho**ab + 0.5e0_wp*phi0*sterm
                dbigefdrho = -ec*(1.e0_wp-ab*log(rho))*(ab)*rho**(ab-1.e0_wp) &  ! Deri
                    & -ec*(-ab*(1.e0_wp/rho))*rho**ab &
                    & +0.5e0_wp*phi0*dstermdrho
            endif
            ! Combine derivs
            do j=1, natom
                dvdx(j) = dvdx(j) + d_v2(1,j) + dbigefdrho*d_rho(1,j)
                dvdy(j) = dvdy(j) + d_v2(2,j) + dbigefdrho*d_rho(2,j)
                dvdz(j) = dvdz(j) + d_v2(3,j) + dbigefdrho*d_rho(3,j)
            enddo
            v = v + v2 + bigef
        enddo
    end subroutine pot
end module chempotpy_aln_np_bd_pes

!> @brief Expose Aln_NP_Bd through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_np_bd_pes, only: evaluate
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
