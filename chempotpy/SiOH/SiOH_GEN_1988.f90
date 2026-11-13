!   System:                        SiOH
!   Functional form:               Two and three body interaction plus C
!   Common name:                   SiOH
!   Number of derivatives:         1
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:
!   Notes:  Many-body silica  potential energy function.  The functional
!           form is from Ref. 1.
!  References: (1) JApplPhys89(11) 2001 6013
!              (2) JChemPhys 89(9) 1988 5818
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  dvdx,dvdy,dvdz   --- one-dimensional arrays of the gradients with res
!                       to the Cartesian coordinates of the system (outp
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximum number of atoms (input)
!  ien     --- 1: just take coulombic interaction energy
!          --- 0: take all energy

!> @brief Implement SiOH_GEN_1988, a potential energy surface for SiOH.
module chempotpy_sioh_gen_1988_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the SiOH_GEN_1988 potential energy surface for SiOH.
    subroutine evaluate( xyz, nsi, no, nh, natoms, igrad, p, g, d )
        use chempotpy_units, only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, intent(in)   :: nsi, no, nh
        integer, intent(in)   :: natoms
        real(wp), intent(in)  :: xyz(10000,3)
        integer, intent(in)   :: igrad
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: cx(10000), cy(10000), cz(10000)
        real(wp)              :: dcx(10000), dcy(10000), dcz(10000)
        real(wp)              :: q(10000)
        integer               :: indexlist(10000)
        real(wp)              :: v
        integer               :: istate, iatom
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        q = 0.e0_wp
        do iatom=1, natoms
            cx(iatom) = xyz(iatom,1)/bohr_to_angstrom
            cy(iatom) = xyz(iatom,2)/bohr_to_angstrom
            cz(iatom) = xyz(iatom,3)/bohr_to_angstrom
        enddo
        do iatom=1, nsi
            indexlist(iatom) = 14
        enddo
        do iatom=1, no
            indexlist(iatom+nsi) = 8
        enddo
        do iatom=1, nh
            indexlist(iatom+nsi+no) = 1
        enddo
        call siohpot(cx, cy, cz, q, indexlist, v, dcx, dcy, dcz, natoms, 10000, 0)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = v*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = v*hartree_to_electronvolt
            enddo
            do istate=1, nstates
                do iatom=1, natoms
                    g(istate,iatom,1) = dcx(iatom)*hartree_per_bohr_to_ev_per_angstrom
                    g(istate,iatom,2) = dcy(iatom)*hartree_per_bohr_to_ev_per_angstrom
                    g(istate,iatom,3) = dcz(iatom)*hartree_per_bohr_to_ev_per_angstrom
                enddo
            enddo
        else
            write (*, *) 'Only energy and gradient are available'
        endif
    end subroutine evaluate

    !> @brief Evaluate siohpot for the SiOH_GEN_1988 potential.
    subroutine siohpot( x, y, z, q, indmm, v, dvdx, dvdy, dvdz, natom, &
            & maxatom, ien )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: aa(14,14,3), alambda(14,14,14), alpha(14,14), bb(14,14,3), beta(14,14), cc(14,14,3), &
            & dx(maxatom, maxatom), dy(maxatom, maxatom), dz(maxatom, maxatom), e, e1, e2, e3, f, f3, f4, f5, g, gamma(22, 36), h, &
            & p, q(maxatom), r0(22, 36), rho(14, 14), rij(maxatom, maxatom), theta(14, 14, 14), tmp, v, v1, v2, v3
        integer             :: i, ien, indmm(maxatom), j, k, maxatom, natom
        real(wp), parameter :: rn    = 10.393496649e0_wp
        real(wp), parameter :: pi    = 3.1415926536e0_wp
        real(wp), parameter :: pener = 4.3597482e-11_wp
        real(wp), parameter :: pdis  = 0.52917706e0_wp
        real(wp), parameter :: pang  = 57.2957795e0_wp
        real(wp)            :: dvdy(maxatom), dvdz(maxatom), dvdx(maxatom), &
            & x(maxatom), y(maxatom), z(maxatom)
        ! Nanoparticle-parameterized Parameters in two-body term(Ref. 1)
        alpha(14,14) = 0.18770e-8_wp/pener
        alpha(14,8)  = 0.29620e-8_wp/pener
        alpha(8,14)  = alpha(14,8)
        alpha(14,1)  = 0.00690e-8_wp/pener
        alpha(1,14)  = alpha(14,1)
        alpha(8,8)   = 0.07250e-8_wp/pener
        alpha(8,1)   = 0.03984e-8_wp/pener
        alpha(1,8)   = alpha(8,1)
        alpha(1,1)   = 0.00340e-8_wp/pener
        beta(14,14)  = 2.29e0_wp/pdis
        beta(14,8)   = 2.34e0_wp/pdis
        beta(8,14)   = beta(14,8)
        beta(14,1)   = 2.31e0_wp/pdis
        beta(1,14)   = beta(14,1)
        beta(8,8)    = 2.34e0_wp/pdis
        beta(8,1)    = 2.26e0_wp/pdis
        beta(1,8)    = beta(8,1)
        beta(1,1)    = 2.10e0_wp/pdis
        rho(14, 14)  = 0.29e0_wp/pdis
        rho(14, 8)   = 0.29e0_wp/pdis
        rho(8, 14)   = rho(14, 8)
        rho(14, 1)   = 0.29e0_wp/pdis
        rho(1, 14)   = rho(14, 1)
        rho(8, 8)    = 0.29e0_wp/pdis
        rho(8, 1)    = 0.29e0_wp/pdis
        rho(1, 8)    = rho(8, 1)
        rho(1, 1)    = 0.35e0_wp/pdis
        aa(14,1,1)   = -4.6542e-12_wp/pener
        aa(1,14,1)   = aa(14,1,1)
        aa(14,1,2)   = 0.e0_wp
        aa(1,14,2)   = aa(14,1,2)
        aa(14,1,3)   = 0.e0_wp
        aa(1,14,3)   = aa(14,1,3)
        aa(1,1,1)    = -5.2793e-12_wp/pener
        aa(1,1,2)    = 0.3473e-12_wp/pener
        aa(1,1,3)    = 0.e0_wp
        aa(8,1,1)    = -2.0840e-12_wp/pener
        aa(8,1,2)    = 7.6412e-12_wp/pener
        aa(8,1,3)    = -0.8336e-12_wp/pener
        aa(1,8,1)    = aa(8,1,1)
        aa(1,8,2)    = aa(8,1,2)
        aa(1,8,3)    = aa(8,1,3)
        bb(14,1,1)   = 6.0e0_wp*pdis
        bb(14,1,2)   = 0.e0_wp
        bb(14,1,3)   = 0.e0_wp
        bb(1,14,1)   = bb(14,1,1)
        bb(1,14,2)   = bb(14,1,2)
        bb(1,14,3)   = bb(14,1,3)
        bb(1,1,1)    = 6.0e0_wp*pdis
        bb(1,1,2)    = 2.0e0_wp*pdis
        bb(1,1,3)    = 0.e0_wp
        bb(8,1,1)    = 15.0e0_wp*pdis
        bb(8,1,2)    = 3.2e0_wp*pdis
        bb(8,1,3)    = 5.0e0_wp*pdis
        bb(1,8,1)    = bb(8,1,1)
        bb(1,8,2)    = bb(8,1,2)
        bb(1,8,3)    = bb(8,1,3)
        cc(14,1,1)   = 2.20e0_wp/pdis
        cc(14,1,2)   = 0.e0_wp
        cc(14,1,3)   = 0.e0_wp
        cc(1,14,1)   = cc(14,1,1)
        cc(1,14,2)   = cc(14,1,2)
        cc(1,14,3)   = cc(14,1,3)
        cc(1,1,1)    = 1.51e0_wp/pdis
        cc(1,1,2)    = 2.42e0_wp/pdis
        cc(1,1,3)    = 0.e0_wp
        cc(8,1,1)    = 1.05e0_wp/pdis
        cc(8,1,2)    = 1.50e0_wp/pdis
        cc(8,1,3)    = 2.00e0_wp/pdis
        cc(1,8,1)    = cc(8,1,1)
        cc(1,8,2)    = cc(8,1,2)
        cc(1,8,3)    = cc(8,1,3)
        ! parameters in the three-body term
        alambda(8,14,8)  = 19.0e-11_wp/pener
        alambda(14,8,14) = 0.3e-11_wp/pener
        alambda(14,8,1)  = 5.0e-11_wp/pener
        alambda(1,8,14)  = alambda(14,8,1)
        alambda(1,8,1)   = 35.0e-11_wp/pener
        gamma(22, 30)    = 2.8e0_wp/pdis
        gamma(22, 36)    = 2.0e0_wp/pdis
        gamma(22, 23)    = 2.0e0_wp/pdis
        gamma(9, 23)     = 1.2e0_wp/pdis
        gamma(9, 10)     = 1.3e0_wp/pdis
        r0(22, 30)       = 3.0e0_wp/pdis
        r0(22, 36)       = 2.6e0_wp/pdis
        r0(22, 23)       = 2.6e0_wp/pdis
        r0(9, 10)        = 1.6e0_wp/pdis
        r0(9, 23)        = 1.5e0_wp/pdis
        theta(8, 14, 8)  = 109.5e0_wp/pang
        theta(14, 8, 14) = 109.5e0_wp/pang
        theta(14, 8, 1)  = 109.5e0_wp/pang
        theta(1, 8, 14)  = theta(14, 8, 1)
        theta(1, 8, 1)   = 104.5e0_wp/pang
        ! Initialize
        v1 = 0.e0_wp
        v2 = 0.e0_wp
        v3 = 0.e0_wp
        v  = 0.e0_wp
        do i=1, natom
            dvdx(i) = 0.e0_wp
            dvdy(i) = 0.e0_wp
            dvdz(i) = 0.e0_wp
            do j=i+1, natom
                dx(i, j)  = x(i)-x(j)
                dy(i, j)  = y(i)-y(j)
                dz(i, j)  = z(i)-z(j)
                rij(i, j) = sqrt(dx(i, j)*dx(i, j)+dy(i, j)*dy(i, j)+dz(i, j)*dz(i, j))
                rij(j, i) = rij(i, j)
                dx(j, i)  = -dx(i, j)
                dy(j, i)  = -dy(i, j)
                dz(j, i)  = -dz(i, j)
            enddo
        enddo
        ! Double loop (main loop)
        do i=1, natom
            do j=1, natom
                ! Cutoff function
                if( i .eq. j ) cycle
                if( rij(i, j) .le. rn )then
                    goto_block_111_265: block
                        if( i .gt. j ) exit goto_block_111_265
                        v1=v1+alpha(indmm(i),indmm(j))*exp(-rij(i, j)/rho(indmm(i), indmm(j))) +(q(i)*q(j)/rij(i, j))* &
                            & erfc(rij(i, j)/ beta(indmm(i),indmm(j)))
                        e=-alpha(indmm(i),indmm(j))*exp(-rij(i, j)/rho(indmm(i), indmm(j)))/rho(indmm(i), indmm(j)) - &
                            & q(i)* q(j)* erfc(rij(i, j)/beta(indmm(i),indmm(j)))/ rij(i, j)**2 -2*q(i)*q(j)*exp(- &
                            & (rij(i, j)/ beta(indmm(i),indmm(j)))**2)/(sqrt(pi)*rij(i, j)*beta(indmm(i),indmm(j)))
                        tmp     = e/rij(i, j)
                        dvdx(i) = dvdx(i)+tmp*dx(i, j)
                        dvdy(i) = dvdy(i)+tmp*dy(i, j)
                        dvdz(i) = dvdz(i)+tmp*dz(i, j)
                        dvdx(j) = dvdx(j)-tmp*dx(i, j)
                        dvdy(j) = dvdy(j)-tmp*dy(i, j)
                        dvdz(j) = dvdz(j)-tmp*dz(i, j)
                        e1      = exp(bb(indmm(i),indmm(j),1)*(rij(i, j)-cc(indmm(i),indmm(j),1)))
                        e2      = exp(bb(indmm(i),indmm(j),2)*(rij(i, j)-cc(indmm(i),indmm(j),2)))
                        e3      = exp(bb(indmm(i),indmm(j),3)*(rij(i, j)-cc(indmm(i),indmm(j),3)))
                        v2      = v2+aa(indmm(i),indmm(j),1)/(1+e1)+aa(indmm(i),indmm(j),2) /(1+e2)+aa(indmm(i), &
                            & indmm(j),3)/ (1+e3)
                        f=-aa(indmm(i),indmm(j),1)*bb(indmm(i),indmm(j),1)*e1/(1+e1)**2 -aa(indmm(i),indmm(j),2)* &
                            & bb(indmm(i),indmm(j),2)*e2/(1+e2)**2 -aa(indmm(i),indmm(j),3)*bb(indmm(i),indmm(j),3)* &
                            & e3/ (1+e3)**2
                        tmp     = f/rij(i, j)
                        dvdx(i) = dvdx(i)+tmp*dx(i, j)
                        dvdy(i) = dvdy(i)+tmp*dy(i, j)
                        dvdz(i) = dvdz(i)+tmp*dz(i, j)
                        dvdx(j) = dvdx(j)-tmp*dx(i, j)
                        dvdy(j) = dvdy(j)-tmp*dy(i, j)
                        dvdz(j) = dvdz(j)-tmp*dz(i, j)
                    end block goto_block_111_265
                    v3 = 0.e0_wp
                    do k=j+1, natom
                        if( k .eq. i .or. k .eq. j ) cycle
                        if( indmm(i) .eq. 1 ) cycle
                        if( indmm(i) == 8 .and. (indmm(j) == 8 .or. indmm(k) == 8) ) cycle
                        if( indmm(i) == 14 .and. (indmm(j) /= 8 .or. indmm(k) /= 8) ) cycle
                        if (rij(i, j).ge.r0(indmm(i)+indmm(j), indmm(i)+indmm(j)+indmm(k)) &
                            & .or. rij(i, k).ge.r0(indmm(i)+indmm(k), indmm(i)+indmm(k)+&
                            &indmm(j))) then
                        v3 = 0.e0_wp
                    else
                        g=exp(gamma(indmm(i)+indmm(j), indmm(i)+indmm(j)+indmm(k))/ (rij(i, j)-r0(indmm(i)+indmm(j), &
                            & indmm(i)+ indmm(j)+indmm(k))) +gamma(indmm(i)+indmm(k), indmm(i)+indmm(k)+indmm(j))/ &
                            & (rij(i, k)- r0(indmm(i)+indmm(k), indmm(i)+indmm(k)+indmm(j))))
                        h  = (rij(i, j)**2+rij(i, k)**2-rij(j, k)**2)/(2*rij(i, j)*rij(i, k))
                        p  = h-cos(theta(indmm(j), indmm(i), indmm(k)))
                        v3 = alambda(indmm(j),indmm(i),indmm(k))*g*p**2
                        f3 =-alambda(indmm(j),indmm(i),indmm(k)) *g*gamma(indmm(i)+indmm(j), indmm(i)+indmm(j)+ &
                            & indmm(k))* p** 2 /(rij(i, j)-r0(indmm(i)+indmm(j), indmm(i)+indmm(j)+indmm(k)))**2 +2* &
                            & alambda(indmm(j),indmm(i),indmm(k))*g*p*(rij(i, j)**2 -rij(i, k)**2+rij(j, k)**2)/(2* &
                            & rij(i, k)*rij(i, j)**2)
                        f4=-alambda(indmm(j),indmm(i),indmm(k)) *g*gamma(indmm(i)+indmm(k), indmm(i)+indmm(k)+ &
                            & indmm(j))* p** 2 /(rij(i, k)-r0(indmm(i)+indmm(k), indmm(i)+indmm(k)+indmm(j)))**2 +2* &
                            & alambda(indmm(j),indmm(i),indmm(k))*g*p*(rij(i, k)**2 -rij(i, j)**2+rij(j, k)**2)/(2* &
                            & rij(i, j)*rij(i, k)**2)
                        f5      = -2*alambda(indmm(j),indmm(i),indmm(k))*g*p*rij(j, k) /(rij(i, j)*rij(i, k))
                        tmp     = f3/rij(i, j)
                        dvdx(i) = dvdx(i)+tmp*dx(i, j)
                        dvdy(i) = dvdy(i)+tmp*dy(i, j)
                        dvdz(i) = dvdz(i)+tmp*dz(i, j)
                        dvdx(j) = dvdx(j)-tmp*dx(i, j)
                        dvdy(j) = dvdy(j)-tmp*dy(i, j)
                        dvdz(j) = dvdz(j)-tmp*dz(i, j)
                        tmp     = f4/rij(i, k)
                        dvdx(i) = dvdx(i)+tmp*dx(i, k)
                        dvdy(i) = dvdy(i)+tmp*dy(i, k)
                        dvdz(i) = dvdz(i)+tmp*dz(i, k)
                        dvdx(k) = dvdx(k)-tmp*dx(i, k)
                        dvdy(k) = dvdy(k)-tmp*dy(i, k)
                        dvdz(k) = dvdz(k)-tmp*dz(i, k)
                        tmp     = f5/rij(j, k)
                        dvdx(j) = dvdx(j)+tmp*dx(j, k)
                        dvdy(j) = dvdy(j)+tmp*dy(j, k)
                        dvdz(j) = dvdz(j)+tmp*dz(j, k)
                        dvdx(k) = dvdx(k)-tmp*dx(j, k)
                        dvdy(k) = dvdy(k)-tmp*dy(j, k)
                        dvdz(k) = dvdz(k)-tmp*dz(j, k)
                endif
                v = v+v3
            enddo
            endif
        enddo
        enddo
        if( ien .eq. 1 )then
            v = v1
        else
            v = v+v1+v2
        endif
    end subroutine siohpot
end module chempotpy_sioh_gen_1988_pes

!> @brief Expose SiOH_GEN_1988 through the Python-facing PES entry point.
subroutine pes( xyz, nsi, no, nh, natoms, igrad, p, g, d )
    use chempotpy_sioh_gen_1988_pes, only: evaluate
    use chempotpy_kinds,             only: wp
    implicit none
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: nsi, no, nh, natoms, igrad
    real(wp), intent(in)  :: xyz(10000,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(xyz, nsi, no, nh, natoms, igrad, p, g, d)
end subroutine pes
