!   System:          ALC
!   Functional form: eLEPS + Va
!                    Extended LEPS plus an auxiliary barrier term
!   Common name: CHC
!   Reference:       unpublished
!   Cross Reference: M. M. Kreevoy, D. Ostovic,
!                    D. G. Truhlar, and B. C. Garrett
!                    J. Phys. Chem. 90, 3766-3774 (1986)
!   Note:            In A-L-C, A = acceptor group, L = H (D),
!                    and C = methyl donor
!   Calling Sequence:
!      PREPOT - initializes the potential's variables and
!               must be called once before any calls to POT
!      POT    - driver for the evaluation of the energy and the derivati
!               of the energy with respect to the coordinates for a give
!               geometric configuration
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      ground electronic state
!   Zero of energy:
!      The classical potential energy is set equal to zero for the A
!      infinitely far from the LC diatomic and R(LC) set equal to the
!      LC equilibrium diatomic value.
!   Parameters:
!      Set in the BLOCK DATA subprogram PTPARM
!   Coordinates:
!      Internal, Definition: R(1) = R(A-L)
!                            R(2) = R(L-C)
!                            R(3) = R(A-C)
!   Common Blocks (used between the calling program and this potential):
!      /PT1CM/ R(3), ENERGY, DEDR(3)
!        passes the coordinates, ground state electronic energy, and
!        derivatives of the ground electronic state energy with respect
!        to the coordinates.
!      /PT2CM/ NSURF, NDER, NDUM(21)
!        passes the control flags where
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
!        Similarly, the terms EASYBC and EASYAC represent the energies
!        in the BC and the AC valleys, respectively.
!   Default Parameter Values:
!      Variable      Default value
!      NSURF            0
!      NDER             1
!      IPTPRT           6

!> @brief Implement ALC_LEPS_ModelSurface, a potential energy surface for ALC.
module chempotpy_alc_leps_modelsurface_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the ALC_LEPS_ModelSurface potential energy surface for ALC.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_legacy_units,       only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: legacy_triatomic_distances_to_bohr
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: dedrr(3), e, rr(3), v
        integer                       :: i, iatom, idir, istate, j
        double precision              :: tx(9), r2(3), dedx(9), drdx(3,9)
        !initialize
        v = 0.d0; g = 0.d0; d = 0.d0
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is ClHH
        call legacy_triatomic_distances_to_bohr(x, rr)
        dedrr = 0.d0
        call pot(rr, igrad, e, dedrr)
        e     = e*legacy_hartree_to_electronvolt
        dedrr = dedrr*51.422067
        r2    = rr*legacy_bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(dedrr, drdx, dedx)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = e
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = e
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
end module chempotpy_alc_leps_modelsurface_pes

!> @brief Expose ALC_LEPS_ModelSurface through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_alc_leps_modelsurface_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX for the ALC_LEPS_ModelSurface potential.
subroutine evdrdx( x, r, drdx )
    use chempotpy_geometry, only: ordered_pair_distance_jacobian
    implicit none
    integer, parameter            :: first_atoms(3)  = [1, 2, 1]
    integer, parameter            :: second_atoms(3) = [2, 3, 3]
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    call ordered_pair_distance_jacobian(x, first_atoms, second_atoms, r, drdx)
end subroutine evdrdx

!> @brief Evaluate PREPOT for the ALC_LEPS_ModelSurface potential.
module chempotpy_alc_leps_modelsurface_pot_legacy_state_05ca95c0
    implicit none
    private
    public :: evaluate
    double precision, parameter :: fit_d(3)    = [105.5d0, 105.5d0, 84.4d0]
    double precision, parameter :: fit_re(3)   = [1.09d0, 1.09d0, 1.54d0]
    double precision, parameter :: fit_beta(3) = [1.78d0, 1.78d0, 1.94d0]
    double precision, parameter :: fit_z(3)    = [0.2d0, 0.2d0, 0.1d0]
    double precision, parameter :: fit_an1     = 0.026d0
    double precision, parameter :: fit_fns     = 0.004d0

    contains

    !> @brief Evaluate dispatch.
    subroutine evaluate( rr, nder, e, dedrr )
        implicit none
        double precision, parameter   :: root2 = 1.41421356d0
        double precision, parameter   :: ckcal = 627.5095d0
        double precision, parameter   :: cangs = 0.529177106d0
        double precision, intent(in)  :: rr(3)
        integer, intent(in)           :: nder
        double precision, intent(out) :: e, dedrr(3)
        double precision              :: an1, b(3), beta(3), coul(3), d(3), delta, do4z(3), ens, ensp1, ensp2, ensp3, &
            & enspr, enspwn, exch(3), exns, fns, op3z(3), r, r2, r3, rad, re(3), s, t, top3z(3), tzp3(3), wn, wn2, &
            & wnp1, wnp2, wnp3, wnt, x(3), x21, x22, x23, z(3), zp3(3), zpo(3)
        integer :: i
        !   Echo the potential parameters to the file linked to FORTRAN unit IPT
        d    = fit_d
        re   = fit_re
        beta = fit_beta
        z    = fit_z
        an1  = fit_an1
        fns  = fit_fns
        do i = 1, 3
            !   CONVERT TO ATOMIC UNITS
            d(i)    = d(i)/ckcal
            re(i)   = re(i)/cangs
            beta(i) = beta(i)*cangs
            !   COMPUTE USEFUL CONSTANTS
            zpo(i)   = 1.0d0 + z(i)
            op3z(i)  = 1.0d0 + 3.0d0*z(i)
            top3z(i) = 2.0d0*op3z(i)
            zp3(i)   = z(i) + 3.0d0
            tzp3(i)  = 2.0d0*zp3(i)
            do4z(i)  = d(i)/4.0d0/zpo(i)
            b(i)     = beta(i)*do4z(i)*2.0d0
        enddo
        !    Compute the eLEPS part of the potential energy
        do i = 1, 3
            x(i)    = exp(-beta(i)*(rr(i)-re(i)))
            coul(i) = do4z(i)*(zp3(i)*x(i)-top3z(i))*x(i)
            exch(i) = do4z(i)*(op3z(i)*x(i)-tzp3(i))*x(i)
        enddo
        rad = sqrt((exch(1)-exch(2))**2+(exch(2)-exch(3))**2+ (exch(3)-exch(1))**2)
        e   = coul(1) + coul(2) + coul(3) - (rad/root2) + d(2)
        !   Compute the derivative of the eLEPS part of the potential, if NDER=1
        goto_block_700_249: block
            if( nder .ne. 1 ) exit goto_block_700_249
            s = exch(1) + exch(2) + exch(3)
            do i = 1, 3
                dedrr(i) = 0.d0
                if( x(i) .lt. 1.d-30 ) cycle
                t = (3.0d0*exch(i)-s)/root2*(op3z(i)*x(i)-zp3(i))
                !   If RAD is very small compared to T print a warning to the file linke
                !   to FORTRAN IPTPRT and set T/RAD equal to T.
                if( abs(rad) .lt. 1.d-32.and.abs(t) .gt. 1.d-12 )then
                    write(*, '(/,1X,''In the potential routine: T, RAD = '',1P,2E15.7, /,1X,''T/RAD has been set equal to' // &
                        ' T'')') t, rad
                elseif( abs(rad) .gt. 1.d-32 )then
                    t = t/rad
                endif
                dedrr(i) = b(i)*x(i)*(t -zp3(i)*x(i)+op3z(i))
            enddo
        end block goto_block_700_249
        !   Compute the auxiliary barrier term, Va, where
        !   Va = a{[(R1-R2)*(R2-R3)*(R3-R1)]^2}*exp[-alpha*{(R1+R2+R3)^3}]
        r    = rr(1)+rr(2)+rr(3)
        r2   = r*r
        r3   = r2*r
        exns = exp(-fns*r3)
        ens  = 0.d0
        wnt  = (rr(1)-rr(2))*(rr(2)-rr(3))*(rr(3)-rr(1))
        wn   = abs(wnt)
        wn2  = wn*wn
        ens  = (an1*wn2)*exns
        !    Add the energy terms together
        e = e+ens
        !    Compute the derivative of the Va term with respect to the
        !    coordinates, if NDER = 1
        if( nder .eq. 1 )then
            delta = -1.d0
            if( wn .eq. wnt ) delta =1.d0
            ensp1  = 0.d0
            ensp2  = 0.d0
            ensp3  = 0.d0
            enspwn = (an1*wn*2.d0)*exns
            enspr  = (-3.d0)*fns*r2*ens
            ! DENS/DXI=(DENS/DWN)(DWN/DXI)+(DENS/DR)(DR/DXI)
            x21      = rr(1)*rr(1)
            x22      = rr(2)*rr(2)
            x23      = rr(3)*rr(3)
            wnp1     = (2.d0*rr(1)*(rr(3)-rr(2))+x22-x23)*delta
            wnp2     = (2.d0*rr(2)*(rr(1)-rr(3))+x23-x21)*delta
            wnp3     = (2.d0*rr(3)*(rr(2)-rr(1))+x21-x22)*delta
            ensp1    = enspwn*wnp1+enspr
            ensp2    = enspwn*wnp2+enspr
            ensp3    = enspwn*wnp3+enspr
            dedrr(1) = dedrr(1)+ensp1
            dedrr(2) = dedrr(2)+ensp2
            dedrr(3) = dedrr(3)+ensp3
        endif
    end subroutine evaluate
end module chempotpy_alc_leps_modelsurface_pot_legacy_state_05ca95c0

!> @brief Compatibility interface for pot.
subroutine pot( rr, nder, e, dedrr )
    use chempotpy_alc_leps_modelsurface_pot_legacy_state_05ca95c0, only: evaluate
    implicit none
    double precision, intent(in)  :: rr(3)
    integer, intent(in)           :: nder
    double precision, intent(out) :: e, dedrr(3)
    call evaluate(rr, nder, e, dedrr)
end subroutine pot

!> @brief Evaluate POTINFO for the ALC_LEPS_ModelSurface potential.
