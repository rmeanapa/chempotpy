! System:          ClH2
! Functional form: eLEPS (extended London-Eyring-Polanyi-Sato)
!                  plus E3C (a three-center term)
!                  and modifications to improve the saddle point
!                  bend potential; based on GQQ surface
! Common name:     G3
! Reference:       Allison, T. C.; Lynch, G. C.; Truhlar, D. G.;
!                  Gordon, M. S. J. Phys. Chem., to be submitted.
! The evaluator is self-contained; all fitted parameters are embedded in
! the owning module.
! The coordinates are defined as follows:
!           R(1) = R(H-Cl)
!           R(2) = R(H-H')
!           R(3) = R(H'-Cl)
! The zero of energy for this potential occurs when the Cl atom is
! "infinitely" far from the the H2 diatom, and the H2 diatom is at
! its equilibrium seperation.
!   Common Blocks (used between the calling program and this potential):
!      /PT1CM/ R(3), ENERGY, DEDR(3)
!        passes the coordinates, ground state electronic energy, and
!        derivatives of the ground electronic state energy with respect
!        to the coordinates.
!      /PT2CM/ NSURF, NDER, NDUM(21)
!        passes the control flags where
!        NSURF = 0    ground electronic surface
!        NDER  = 0 => no derivatives are computed
!              = 1 => derivatives of the energy for the ground electroni
!                     state with respect to the coordinates are computed
!        NDUM  - not used
!      /PT4CM/ IPTPRT, IDUM(19)
!        IPTPRT passes the FORTRAN unit number used for potential output
!        IDUM   not used
!      /PT5CM/ EASYAB, EASYBC, EASYAC
!        passes the energy in the three asymptotic valleys for an A + BC
!        The energy in the AB valley, EASYAB, is equal to the energy of
!        C atom "infinitely" far from the AB diatomic and R(AB) set equa
!        Re(AB), the equilibrium bond length for the AB diatomic.
!        In this potential the AB valley represents H infinitely far fro
!        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the t
!        EASYBC and EASYAC represent the energies in the H2 and the BrH
!        valleys, respectively.
!   Default Parameter Values:
!      Variable      Default value
!      NSURF            0
!      NDER             1
!      IPTPRT           6

!> @brief Implement ClH2_LEPS_ALTG1_1996, a potential energy surface for ClH2.
module chempotpy_clh2_leps_altg1_1996_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the ClH2_LEPS_ALTG1_1996 potential energy surface for ClH2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,              only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_units,              only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: triatomic_distances_to_bohr
        use chempotpy_kinds,              only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: dedr(3), energy, r(3), v
        integer               :: i, iatom, idir, istate, j
        real(wp)              :: tx(9), r2(3), dedx(9), drdx(3,9)
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is ClHH
        call triatomic_distances_to_bohr(x, r)
        dedr = 0.e0_wp
        call pot(r, igrad, energy, dedr)
        energy = energy*hartree_to_electronvolt
        dedr   = dedr*hartree_per_bohr_to_ev_per_angstrom
        r2     = r*bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(dedr, drdx, dedx)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = energy
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = energy
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    j               = (iatom-1)*3+idir
                    g(1,iatom,idir) = dedx(j)
                enddo
            enddo
        elseif( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
    end subroutine evaluate
end module chempotpy_clh2_leps_altg1_1996_pes

!> @brief Expose ClH2_LEPS_ALTG1_1996 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_clh2_leps_altg1_1996_pes, only: evaluate
    use chempotpy_kinds,                    only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX for the ClH2_LEPS_ALTG1_1996 potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    use chempotpy_kinds,    only: wp
    implicit none
    integer, parameter    :: first_atoms(3)  = [1, 2, 1]
    integer, parameter    :: second_atoms(3) = [2, 3, 3]
    real(wp), intent(in)  :: x(9), r(3)
    real(wp), intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!> @brief Own the released fit state for the ClH2_LEPS_ALTG1_1996 potential.
module chempotpy_clh2_leps_altg1_1996_pot_legacy_state_46e6d832
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    ! Private evaluator state preserves the released unit-conversion order.
    real(wp)            :: de(3) = [106.447e0_wp, 109.458e0_wp, 106.447e0_wp]
    real(wp)            :: re(3) = [1.2732e0_wp, 0.74127e0_wp, 1.2732e0_wp]
    real(wp)            :: b(3)  = [1.8674e0_wp, 1.9413e0_wp, 1.8674e0_wp]
    real(wp), parameter :: s(3)   = [0.1835e0_wp, 0.167e0_wp, 0.1835e0_wp]
    real(wp), parameter :: rj     = 0.0758016022e0_wp
    real(wp), parameter :: ap     = 0.0008627355e0_wp
    real(wp), parameter :: at     = 0.2981969860e0_wp
    real(wp), parameter :: bp     = 0.1439106543e0_wp
    real(wp), parameter :: q2     = 0.6940323070e0_wp
    real(wp), parameter :: q4     = 1.6963092005e0_wp
    real(wp), parameter :: c1     = 3.1053877618397071175758280546e+0_wp
    real(wp), parameter :: c2     = -1.8942608491155350536449828878e+0_wp
    real(wp), parameter :: ckcau  = 627.5095e0_wp
    real(wp), parameter :: cangau = 0.52917706e0_wp

    contains

    !> @brief Evaluate the ClH2 ALTG1 potential for one set of bond distances.
    subroutine evaluate( r, nder, energy, dedr )
        ! Last modified on 3 August 1995 by TCA.
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: r(3)
        integer, intent(in)   :: nder
        real(wp), intent(out) :: energy, dedr(3)
        real(wp)              :: ctheta, df1r1, df2r1, df3r1, &
            & df4rr, dfr1, dfr2, dfr3, dfrr, dfrs, dfrt1, dfrt2, dfrt3, dj1r1, dj1r2, dj1r3, dj2r2, &
            & dj3r1, dj3r2, dj3r3, f1, f123, f124, f134, f2, f234, f3, f4, fmod, g, g5, hradi, pf(3), &
            & pfd(3), pfs(3), r1i, r1s, r2i, r2s, r3i, r3s, rad, st2, stheta, wj1, wj12d, &
            & wj2, wj23d, wj3, wj31d, wq1, wq2, wq3, x, x2, y1, y1s, y2, y3, &
            & y3s
        integer :: i
        ! Constants used by the code.  de is the dissociation energy (kcal/mol),
        ! re is the equilibrium bond length (angstroms), b is the Morse beta
        ! the appropriate diatom (1 = H-Cl, 2 = H=H', 3 = H'-Cl).  The other
        ! constants may be found in the literature reference.  ckcau is the
        ! conversion factor between kcal/mol and hartrees and cangau is the
        ! conversion factor between angstroms and bohr.
        ! Calculate first derivatives when requested by nder.
        ! Convert to atomic units and calculate some useful constants.
        do i=1, 3
            de(i)  = de(i)/ckcau
            re(i)  = re(i)/cangau
            b(i)   = b(i)*cangau
            pf(i)  = 0.5e0_wp*de(i)*((1.e0_wp-s(i))/(1.e0_wp+s(i)))
            pfs(i) = de(i)+pf(i)
            pfd(i) = de(i)-pf(i)
        enddo
        r1i = 1.e0_wp/r(1)
        r2i = 1.e0_wp/r(2)
        r3i = 1.e0_wp/r(3)
        r1s = r(1)*r(1)
        r2s = r(2)*r(2)
        r3s = r(3)*r(3)
        ! Compute three-center term.
        f1     = exp(-ap*(r(1)+r(3))**4)
        f2     = exp(-at*(r(1)-r(3))**2)
        f3     = exp(-bp*(r(1)+r(3)-r(2))**2)
        ctheta = (r1s+r3s-r2s)*(0.5e0_wp*r1i*r3i)
        if( ctheta .lt. -1.e0_wp )then
            ctheta = -1.e0_wp
        elseif( ctheta .gt. 1.e0_wp )then
            ctheta = 1.e0_wp
        endif
        stheta = sin(acos(ctheta))
        st2    = stheta*stheta
        f4     = 1.e0_wp+q2*st2+q4*st2*st2
        ! Compute modification to bending potential.
        x    = (r(1)-r(3))*r2i
        x2   = x*x
        g5   = (0.5e0_wp*(1.e0_wp+ctheta))**5
        g    = g5*(0.5e0_wp*(1.e0_wp+ctheta))
        fmod = 1.e0_wp+g*(c1*(1.e0_wp-x2)+c2*(1.e0_wp-x2*x2))
        ! Compute coulomb and exchange terms.
        y1  = exp(b(1)*(re(1)-r(1)))
        wq1 = (0.5e0_wp*(de(1)+pf(1)*fmod)*y1-pfd(1))*y1
        wj1 = (0.5e0_wp*(de(1)-pf(1)*fmod)*y1-pfs(1))*y1
        y2  = exp(b(2)*(re(2)-r(2)))
        wq2 = (0.5e0_wp*pfs(2)*y2-pfd(2))*y2
        wj2 = (0.5e0_wp*pfd(2)*y2-pfs(2))*y2
        y3  = exp(b(3)*(re(3)-r(3)))
        wq3 = (0.5e0_wp*(de(3)+pf(3)*fmod)*y3-pfd(3))*y3
        wj3 = (0.5e0_wp*(de(3)-pf(3)*fmod)*y3-pfs(3))*y3
        ! Compute potential energy.
        rad    = sqrt(wj1*(wj1-wj2)+wj2*(wj2-wj3)+wj3*(wj3-wj1))
        energy = de(2)+wq1+wq2+wq3-rad+f4*(rj*f1*f2*f3)
        if( nder .eq. 1 )then
            ! Compute contributions to the derivatives.
            dfrr  = (1.e0_wp-x2)*(c1+c2*(1.e0_wp+x2))*3.e0_wp*g5
            dfrs  = 2.e0_wp*x*g*(c1+2.e0_wp*c2*x2)*r2i
            dfrt1 = r3i-ctheta*r1i
            dfrt2 = -r(2)*(r1i*r3i)
            dfrt3 = r1i-ctheta*r3i
            dfr1  = dfrr*dfrt1-dfrs
            dfr2  = dfrr*dfrt2+dfrs*x
            dfr3  = dfrr*dfrt3+dfrs
            y1s   = y1*y1
            y3s   = y3*y3
            dj1r1 = y1*(de(1)*b(1)*(1.e0_wp-y1) -pf(1)*(-b(1)*(1.e0_wp+fmod*y1)+0.5e0_wp*dfr1*y1))
            dj1r2 = -0.5e0_wp*pf(1)*dfr2*y1s
            dj1r3 = -0.5e0_wp*pf(1)*dfr3*y1s
            dj2r2 = b(2)*y2*(pfs(2)-pfd(2)*y2)
            dj3r1 = -0.5e0_wp*pf(3)*dfr1*y3s
            dj3r2 = -0.5e0_wp*pf(3)*dfr2*y3s
            dj3r3 = y3*(de(3)*b(3)*(1.e0_wp-y3) -pf(3)*(-b(3)*(1.e0_wp+fmod*y3)+0.5e0_wp*dfr3*y3))
            df1r1 = -2.e0_wp*at*f1*(r(1)-r(3))
            df2r1 = -4.e0_wp*ap*f2*(r(1)+r(3))**3
            df3r1 = -2.e0_wp*bp*f3*(r(1)+r(3)-r(2))
            df4rr = 2.e0_wp*ctheta*(q2+2.e0_wp*q4*st2)
            f123  = f1*f2*f3
            f124  = f1*f2*f4
            f134  = f1*f3*f4
            f234  = f2*f3*f4
            wj12d = wj1-wj2
            wj23d = wj2-wj3
            wj31d = wj3-wj1
            hradi = 0.5e0_wp/rad
            ! Assemble derivatives.
            dedr(1)=y1*(de(1)*b(1)*(1.e0_wp-y1)+pf(1)*(-b(1)*(1.e0_wp+fmod*y1) +0.5e0_wp*dfr1*y1))+0.5e0_wp*pf(3)* &
                & dfr1*y3s - hradi* (wj12d*dj1r1-wj23d*dj3r1+wj31d*(dj3r1-dj1r1)) +rj*(f124*df3r1+f134*df2r1+f234* &
                & df1r1-f123* df4rr*dfrt1)
            dedr(2)=0.5e0_wp*pf(1)*dfr2*y1s+b(2)*y2*(pfd(2)-pfs(2)*y2) +0.5e0_wp*pf(3)*dfr2*y3s -hradi*(wj12d*(dj1r2- &
                & dj2r2)+ wj23d* (dj2r2-dj3r2) +wj31d*(dj3r2-dj1r2))-rj*(f123*df4rr*dfrt2+f124*df3r1)
            dedr(3)=0.5e0_wp*pf(1)*dfr3*y1s+y3*(de(3)*b(3)*(1.e0_wp-y3) +pf(3)*(-b(3)*(1.e0_wp+fmod*y3)+0.5e0_wp*dfr3* &
                & y3)) - hradi* (wj12d*dj1r3-wj23d*dj3r3+wj31d*(dj3r3-dj1r3)) +rj*(f124*df3r1+f134*df2r1-f234*df1r1- &
                & f123* df4rr*dfrt3)
        endif
    end subroutine evaluate
end module chempotpy_clh2_leps_altg1_1996_pot_legacy_state_46e6d832

!> @brief Evaluate the ClH2 ALTG1 potential for one set of bond distances.
subroutine pot( r, nder, energy, dedr )
    use chempotpy_clh2_leps_altg1_1996_pot_legacy_state_46e6d832, only: evaluate
    use chempotpy_kinds,                                          only: wp
    implicit none
    real(wp), intent(in)  :: r(3)
    integer, intent(in)   :: nder
    real(wp), intent(out) :: energy, dedr(3)
    call evaluate(r, nder, energy, dedr)
end subroutine pot

!> @brief Evaluate POTINFO for the ClH2_LEPS_ALTG1_1996 potential.
