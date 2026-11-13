!   System:                        NaH2
!   Functional form:               Modified Extended LEPS (2x2 diabatic fit)
!   Common name:                   NaH2 surface set 7l
!   Interface:                     3-1V
!   Number of electronic surfaces: 2
!   Number of derivatives:         0
!   Notes:              This version contains long range forces,
!                       no derivatives.  This fit is based on the surfac
!                       set 6 fit.
!   References:         Unpublished,  M. D. Hack and D. G. Truhlar (2000
!   Protocol: POT evaluates the stateless diabatic surfaces.
!               of the energy with respect to the coordinates for a give
!               geometric configuration
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      Ground diabatic electronic state, first-excited diabatic electron
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Na
!      infinitely far from the H2 diatomic and R(H2) set equal to the
!      H2 equilibrium diatomic value.
!   Parameters:
!   Coordinates:
!      Internal, Definition: R(1) = R(Na-H)
!                            R(2) = R(H-H')
!                            R(3) = R(Na-H')
! This is NaH2 potential matrix 7l, finalized on 5/17/99
! This version contains long range forces.
!  all surfaces.  The potential is called by calling POT with the parame
!  r, e, nt, nsurf.  r(nt,3) and ei(nt) are double precision arrays of
!  dimension nt which is an integer representing the number of geometrie
!  input.  The integer nsurf selects the desired surface.  nsurf=1 is th
!  lower diabatic surface, nsurf=2 is the coupling surface, nsurf=3 is
!  the upper diabatic surface, nsurf = 4 selects the lower adiabatic
!  surface, and nsurf = 5 selects the upper adiabatic surface

!> @brief Implement NaH2_LEPS_7L_2000, a potential energy surface for NaH2.
module chempotpy_nah2_leps_7l_2000_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NaH2_LEPS_7L_2000 potential energy surface for NaH2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,              only: hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: triatomic_distances_to_bohr
        use chempotpy_kinds,              only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 2
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        integer               :: nt
        real(wp)              :: r(1,3), v(1)
        real(wp)              :: u11(1), u12(1), u22(1)
        real(wp)              :: u(nstates,nstates)
        real(wp)              :: t(nstates,nstates)
        integer               :: i, j, k, l
        !initialize
        u  = 0.e0_wp; p = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        nt = 1
        ! input cartesian is ClHH
        call triatomic_distances_to_bohr(x, r(1,:))
        if( igrad == 0 )then
            call pot(r, u11, 1, 1)
            call pot(r, u12, 1, 2)
            call pot(r, u22, 1, 3)
            u(1,1) = u11(1)*hartree_to_electronvolt
            u(1,2) = u12(1)*hartree_to_electronvolt
            u(2,1) = u12(1)*hartree_to_electronvolt
            u(2,2) = u22(1)*hartree_to_electronvolt
            call diagonalize(nstates, u, p, t)
        else
            write(*, *) "Only energy is available"
        endif
    end subroutine evaluate
end module chempotpy_nah2_leps_7l_2000_pes

!> @brief Expose NaH2_LEPS_7L_2000 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_nah2_leps_7l_2000_pes, only: evaluate
    use chempotpy_kinds,                 only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX for the NaH2_LEPS_7L_2000 potential.
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

!> @brief Evaluate diagonalize for the NaH2_LEPS_7L_2000 potential.
subroutine diagonalize( n, a_ss, ev_s, u_ss )
    use chempotpy_kinds, only: wp
    implicit none
    integer, intent(in)   :: n
    real(wp), intent(in)  :: a_ss(n,n)
    real(wp), intent(out) :: u_ss(n,n)
    real(wp), intent(out) :: ev_s(n)
    integer               :: io, i
    real(wp), allocatable :: work_d(:)
    integer               :: lwork_d
    u_ss    = a_ss
    lwork_d = -1
    allocate( work_d(1) )
    call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
    lwork_d = int(work_d(1))
    deallocate( work_d )
    allocate( work_d(lwork_d) )
    call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
end subroutine diagonalize

!> @brief Evaluate pot for the NaH2_LEPS_7L_2000 potential.
subroutine pot( r, e, nt, nsurf )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp) :: e(nt), r(nt,3)
    integer  :: nsurf, nt
    if( nsurf .eq. 1 )then
        call pot11(r, e, nt)
    elseif( nsurf .eq. 2 )then
        call pot12(r, e, nt)
    elseif( nsurf .eq. 3 )then
        call pot22(r, e, nt)
    endif
end subroutine pot

!  U11 surface - lower diabatic surface
!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_nah2_leps_7l_2000_pot11_legacy_state_c810f3bd
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r, e, nt
        real(wp) :: al1, al2, alb, alhh, alt, b0, balp, bect, bect2, beh2, bes1, bf, brho, c6, c60, c61, c63, c690, &
            & c6swt1, c6swt3, cevau, cind1, cind3, cos1, cos2, cosa, coscos, cosg, cstsq, cstsq1, cstsq3, dalp, dect, deh2, &
            & deh2m, des1, des20, dexp1, dexp3, div1, div3, drho, e(:), elec, f0, f90, frac, gamh2a, gamh2b, h2a1, h2a2, h2a3, &
            & h2a4, h2a5, h2a6, h2a7, h2a8, h2a9, h2b1, h2b2, h2b3, h2b4, h2b5, h2b6, h2b7, h2b8, h2b9, h2i1, h2i2, h2i3, h2i4, &
            & h2i5, h2i6, h2ie, hie, hnaa1, hnaa2, hnaa3, hnaa4, hnaa5, hnab1, hnab2, hnab3, hnab4, hnab5, hnai1, hnai2, hnai3, &
            & hnai4, hnai5, hnai6, hnau1, hnau2, hnau3, hnau4, hnau5, hnau6, hpol, onethd, r(:,:), r1s, r2s, r3s, rac2, rave, &
            & rb0, rb90, rd, rect, reh2, res1, rexp1, rexp2, rexp3, rmid, rmuimu1, rmuimu3, rnag, rnag1, rnag1c, rnag3, rnag3c, &
            & rnagc, rnags, rnahie1, rnahie3, rnahmu1, rnahmu3, rnahpar1, rnahpar3, rnahperp1, rnahperp3, rnahpol1, rnahpol1b, &
            & rnahpol3, rnahpol3b, rnaie, rnapol, rsm1, rsm2, rsm3, s1, s2, s2a, s2alp, s2rho, s3, sswt, twothd, vmid, xect, &
            & xeh2, xes1, xes3
        integer :: i, nt
        ! NaH singlet
        ! H2 singlet
        ! H-Na-H term
        ! data blocks for the long range forces
        ! misc constants
        al1       = 0.e0_wp
        al2       = 0.e0_wp
        alb       = 0.e0_wp
        alhh      = 0.e0_wp
        alt       = 0.e0_wp
        b0        = 0.e0_wp
        balp      = 0.e0_wp
        bect      = 0.e0_wp
        bect2     = 0.e0_wp
        beh2      = 0.e0_wp
        bes1      = 0.e0_wp
        bf        = 0.e0_wp
        brho      = 0.e0_wp
        c6        = 0.e0_wp
        c60       = 0.e0_wp
        c61       = 0.e0_wp
        c63       = 0.e0_wp
        c690      = 0.e0_wp
        c6swt1    = 0.e0_wp
        c6swt3    = 0.e0_wp
        cevau     = 0.e0_wp
        cind1     = 0.e0_wp
        cind3     = 0.e0_wp
        cos1      = 0.e0_wp
        cos2      = 0.e0_wp
        cosa      = 0.e0_wp
        coscos    = 0.e0_wp
        cosg      = 0.e0_wp
        cstsq     = 0.e0_wp
        cstsq1    = 0.e0_wp
        cstsq3    = 0.e0_wp
        dalp      = 0.e0_wp
        dect      = 0.e0_wp
        deh2      = 0.e0_wp
        deh2m     = 0.e0_wp
        des1      = 0.e0_wp
        des20     = 0.e0_wp
        dexp1     = 0.e0_wp
        dexp3     = 0.e0_wp
        div1      = 0.e0_wp
        div3      = 0.e0_wp
        drho      = 0.e0_wp
        elec      = 0.e0_wp
        f0        = 0.e0_wp
        f90       = 0.e0_wp
        frac      = 0.e0_wp
        gamh2a    = 0.e0_wp
        gamh2b    = 0.e0_wp
        h2a1      = 0.e0_wp
        h2a2      = 0.e0_wp
        h2a3      = 0.e0_wp
        h2a4      = 0.e0_wp
        h2a5      = 0.e0_wp
        h2a6      = 0.e0_wp
        h2a7      = 0.e0_wp
        h2a8      = 0.e0_wp
        h2a9      = 0.e0_wp
        h2b1      = 0.e0_wp
        h2b2      = 0.e0_wp
        h2b3      = 0.e0_wp
        h2b4      = 0.e0_wp
        h2b5      = 0.e0_wp
        h2b6      = 0.e0_wp
        h2b7      = 0.e0_wp
        h2b8      = 0.e0_wp
        h2b9      = 0.e0_wp
        h2i1      = 0.e0_wp
        h2i2      = 0.e0_wp
        h2i3      = 0.e0_wp
        h2i4      = 0.e0_wp
        h2i5      = 0.e0_wp
        h2i6      = 0.e0_wp
        h2ie      = 0.e0_wp
        hie       = 0.e0_wp
        hnaa1     = 0.e0_wp
        hnaa2     = 0.e0_wp
        hnaa3     = 0.e0_wp
        hnaa4     = 0.e0_wp
        hnaa5     = 0.e0_wp
        hnab1     = 0.e0_wp
        hnab2     = 0.e0_wp
        hnab3     = 0.e0_wp
        hnab4     = 0.e0_wp
        hnab5     = 0.e0_wp
        hnai1     = 0.e0_wp
        hnai2     = 0.e0_wp
        hnai3     = 0.e0_wp
        hnai4     = 0.e0_wp
        hnai5     = 0.e0_wp
        hnai6     = 0.e0_wp
        hnau1     = 0.e0_wp
        hnau2     = 0.e0_wp
        hnau3     = 0.e0_wp
        hnau4     = 0.e0_wp
        hnau5     = 0.e0_wp
        hnau6     = 0.e0_wp
        hpol      = 0.e0_wp
        onethd    = 0.e0_wp
        r1s       = 0.e0_wp
        r2s       = 0.e0_wp
        r3s       = 0.e0_wp
        rac2      = 0.e0_wp
        rave      = 0.e0_wp
        rb0       = 0.e0_wp
        rb90      = 0.e0_wp
        rd        = 0.e0_wp
        rect      = 0.e0_wp
        reh2      = 0.e0_wp
        res1      = 0.e0_wp
        rexp1     = 0.e0_wp
        rexp2     = 0.e0_wp
        rexp3     = 0.e0_wp
        rmid      = 0.e0_wp
        rmuimu1   = 0.e0_wp
        rmuimu3   = 0.e0_wp
        rnag      = 0.e0_wp
        rnag1     = 0.e0_wp
        rnag1c    = 0.e0_wp
        rnag3     = 0.e0_wp
        rnag3c    = 0.e0_wp
        rnagc     = 0.e0_wp
        rnags     = 0.e0_wp
        rnahie1   = 0.e0_wp
        rnahie3   = 0.e0_wp
        rnahmu1   = 0.e0_wp
        rnahmu3   = 0.e0_wp
        rnahpar1  = 0.e0_wp
        rnahpar3  = 0.e0_wp
        rnahperp1 = 0.e0_wp
        rnahperp3 = 0.e0_wp
        rnahpol1  = 0.e0_wp
        rnahpol1b = 0.e0_wp
        rnahpol3  = 0.e0_wp
        rnahpol3b = 0.e0_wp
        rnaie     = 0.e0_wp
        rnapol    = 0.e0_wp
        rsm1      = 0.e0_wp
        rsm2      = 0.e0_wp
        rsm3      = 0.e0_wp
        s1        = 0.e0_wp
        s2        = 0.e0_wp
        s2a       = 0.e0_wp
        s2alp     = 0.e0_wp
        s2rho     = 0.e0_wp
        s3        = 0.e0_wp
        sswt      = 0.e0_wp
        twothd    = 0.e0_wp
        vmid      = 0.e0_wp
        xect      = 0.e0_wp
        xeh2      = 0.e0_wp
        xes1      = 0.e0_wp
        xes3      = 0.e0_wp
        des1      = 0.55118e0_wp
        bes1      = 1.05118e0_wp
        res1      = 3.38583e0_wp
        reh2      = 1.401e0_wp
        bf        = 1.627e0_wp
        b0        = 1.194e0_wp
        gamh2a    = 2.323e0_wp
        gamh2b    = 3.126e0_wp
        deh2      = 4.74806e0_wp
        des20     = 0.996e0_wp
        s2alp     = 3.0e0_wp
        s2rho     = 7.0e0_wp
        dect      = 2.24409e0_wp
        bect      = 0.57874e0_wp
        rect      = 3.46457e0_wp
        bect2     = 0.87143e0_wp
        h2a1      = 2.38029900e0_wp
        h2a2      = 1.47400334e0_wp
        h2a3      = -0.226331771e0_wp
        h2a4      = -2.17339562e0_wp
        h2a5      = 9.54597998e0_wp
        h2a6      = -23.3366749e0_wp
        h2a7      = 24.7115748e0_wp
        h2a8      = -12.8786102e0_wp
        h2a9      = 2.75711758e0_wp
        h2b1      = 1.35164094e0_wp
        h2b2      = -0.254404564e0_wp
        h2b3      = -10.6627055e0_wp
        h2b4      = -14.0112406e0_wp
        h2b5      = -7.80660398e0_wp
        h2b6      = 5.47697753e0_wp
        h2b7      = -6.95875664e0_wp
        h2b8      = 3.64334771e0_wp
        h2b9      = -0.492505596e0_wp
        h2i1      = 160.484e0_wp
        h2i2      = -20.397e0_wp
        h2i3      = 3.913e0_wp
        h2i4      = 0.777e0_wp
        h2i5      = -7.583
        h2i6      = 1.501
        rnaie     = 5.139e0_wp
        rnapol    = 159.267e0_wp
        hpol      = 4.5e0_wp
        hie       = 13.598e0_wp
        hnai1     = -193.341e0_wp
        hnai2     = 16.968e0_wp
        hnai3     = -77.120e0_wp
        hnai4     = -1.024e0_wp
        hnai5     = 6.009e0_wp
        hnai6     = -1.241e0_wp
        hnau1     = 119.880e0_wp
        hnau2     = -0.620e0_wp
        hnau3     = 31.335e0_wp
        hnau4     = 70.414e0_wp
        hnau5     = -13.016e0_wp
        hnau6     = 0.8681e0_wp
        hnaa1     = 0.67e0_wp
        hnaa2     = -92.229e0_wp
        hnaa3     = 51.01019e0_wp
        hnaa4     = -130.37826e0_wp
        hnaa5     = 6.40285e0_wp
        hnab1     = 0.34333e0_wp
        hnab2     = -92.229e0_wp
        hnab3     = -48.54867e0_wp
        hnab4     = -17.30523e0_wp
        hnab5     = 0.99821e0_wp
        balp      = 2.07234e0_wp
        brho      = 2.0e0_wp
        rb0       = 5.1998697814229e7_wp
        rb90      = 5.0e8_wp
        f90       = 1.0e0_wp
        f0        = 1.0e0_wp
        rd        = 9.1538750508008e8_wp
        dalp      = 1.43307e0_wp
        drho      = 5.71654e0_wp
        rac2      = 1.e0_wp/sqrt(2.e0_wp)
        cevau     = 1.e0_wp/27.2113961e0_wp
        onethd    = 1.0e0_wp/3.0e0_wp
        twothd    = 2.0e0_wp/3.0e0_wp
        do i=1, nt
            ! useful quantities calculated first, r1^2, R2^2, R3^2, atom to diatom
            ! center of mass distances, the H-Na-H bond angle and the Jacobian
            ! angle chi.
            r1s   = r(i,1)*r(i,1)
            r2s   = r(i,2)*r(i,2)
            r3s   = r(i,3)*r(i,3)
            rnag  = sqrt(0.5e0_wp*abs(r1s+r3s-0.5e0_wp*r2s))
            rnags = rnag*rnag
            rnagc = rnags*rnag
            rmid  = 0.5e0_wp*(r(i,1)+r(i,3))
            rave  = 0.5e0_wp*(r(i,1)-r(i,3))
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            cos1 = (r1s+r2s-r3s)/(2.0e0_wp*r(i,1)*r(i,2))
            cos2 = (r3s+r2s-r1s)/(2.0e0_wp*r(i,3)*r(i,2))
            cosg = 0.5e0_wp*(cos1-cos2)
            frac = 0.958e0_wp
            ! H to NaH center of masses
            rnag1  = sqrt(abs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))
            rnag3  = sqrt(abs(r2s+frac*(r1s-r2s-r3s)+r3s*frac*frac))
            rnag1c = rnag1*rnag1*rnag1
            rnag3c = rnag3*rnag3*rnag3
            ! the shortest of any two distances . . .
            rsm1  = r(i,3)+(r(i,2)-r(i,3)) *(0.5e0_wp+0.5e0_wp*tanh(0.5e0_wp*(r(i,3)-r(i,2))))
            rsm2  = r(i,3)+(r(i,1)-r(i,3)) *(0.5e0_wp+0.5e0_wp*tanh(0.5e0_wp*(r(i,3)-r(i,1))))
            rsm3  = r(i,1)+(r(i,2)-r(i,1)) *(0.5e0_wp+0.5e0_wp*tanh(0.5e0_wp*(r(i,1)-r(i,2))))
            rexp1 = exp(-r(i,1))
            rexp2 = exp(-r(i,2))
            rexp3 = exp(-r(i,3))
            ! checks on angles . . .
            if( rnag .eq. 0.0e0_wp )then
                cstsq = 1.0e0_wp
            else
                cstsq = (0.5e0_wp*(r3s-r1s)/(r(i,2)*rnag))**2
            endif
            cosa   = (r1s+r3s-r2s)/(2*r(i,1)*r(i,3))
            coscos = cstsq*(0.5e0_wp+0.5e0_wp*cosa)*(1.0e0_wp-exp(-10.0e0_wp*rnags))
            if( rnag1 .eq. 0.0e0_wp )then
                cstsq1 = 1.0e0_wp
            else
                cstsq1 = 0.25e0_wp*((r3s-r2s-r1s+2.0e0_wp*r3s*frac)/(r(i,1)*rnag1)) **2
            endif
            if( rnag3 .eq. 0.0e0_wp )then
                cstsq3 = 1.0e0_wp
            else
                cstsq3 = 0.25e0_wp*((r1s-r2s-r3s+2.0e0_wp*r1s*frac)/(r(i,3)*rnag3)) **2
            endif
            !==============================================================k
            !This section contains the diatomic curves
            ! The first NaH singlet
            xes1 = exp(-bes1*(r(i,1)-res1))
            s1   = des1*onethd*xes1*(xes1+2.0e0_wp)
            ! The H2 singlet
            ! new singlet, designed to minic lower adiabatic behavior
            reh2   = 1.401e0_wp
            bf     = 1.627e0_wp
            b0     = 1.194e0_wp
            gamh2a = 2.323e0_wp
            gamh2b = 3.126e0_wp
            beh2   = bf*(b0-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2) /(bf-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2)
            sswt   = 0.5e0_wp*(1.0e0_wp-tanh((rmid-s2rho)/s2alp))*cosg**2
            deh2m  = deh2 + deh2*(des20 - 1.0e0_wp)*sswt
            xeh2   = exp(-beh2*(r(i,2)-reh2))
            s2     = deh2m*xeh2*(xeh2-2.0e0_wp)
            ! The second NaH singlet
            xes3 = exp(-bes1*(r(i,3)-res1))
            s3   = des1*onethd*xes3*(xes3+2.0e0_wp)
            ! this long section contains data pertinent to long range forces
            !  H2 ionization energy:
            h2ie = (h2i1+h2i2*r(i,2))*exp(-h2i3*r(i,2))+ r2s*(h2i4+h2i5*r(i,2))*exp(-h2i6*r(i,2))-s2a+hie
            !  H2 polarizability:
            !   al1 is the parallel polarizability for h2
            !   al2 is the perpendicaular polarizability for h2
            !   alt is the anisotropic? polarizability
            !   alb is the mean polarizability
            al1=exp(h2a1*(h2a2-r(i,2)))*(h2a3+h2a4*r(i,2)+h2a5*r2s +h2a6*r(i,2)*r2s+h2a7*r2s*r2s+h2a8*r2s*r2s* r(i, &
                & 2) + h2a9* r2s*r2s*r2s) +2.0e0_wp*hpol
            al2=exp(h2b1*(h2b2-r(i,2)))*(h2b3+h2b4*r(i,2)+h2b5*r2s +h2b6*r(i,2)*r2s+h2b7*r2s*r2s+h2b8*r2s*r2s* r(i, &
                & 2) + h2b9* r2s*r2s*r2s) +2.0e0_wp*hpol
            alt  = twothd*(al1-al2)
            alb  = onethd*al1+twothd*al2
            alhh = alb+alt*(1.5e0_wp*cstsq-0.5e0_wp)
            !  The NaH + H dispersion and NaH + H induction terms must go on the low
            !  so that the adiabat is fitted correctly in the region of the coupling
            ! NaH properties
            rnahie1   = (hnai1+hnai2*r(i,1)) /(hnai3+hnai4*r(i,1)+hnai5*r1s+hnai6*r(i,1)*r1s)+rnaie
            rnahie3   = (hnai1+hnai2*r(i,3)) /(hnai3+hnai4*r(i,3)+hnai5*r3s+hnai6*r(i,3)*r3s)+rnaie
            rnahmu1   = (hnau1*r(i,1)+hnau2*r1s) /(hnau3+hnau4*r(i,1)+hnau5*r1s+hnau6*r(i,1)*r1s)
            rnahmu3   = (hnau1*r(i,3)+hnau2*r3s) /(hnau3+hnau4*r(i,3)+hnau5*r3s+hnau6*r(i,3)*r3s)
            rnahpar1  = rnapol+hpol+exp(-hnaa1*r(i,1)) *(hnaa2+hnaa3*r(i,1)+hnaa4*r1s+hnaa5*r1s*r(i,1))
            rnahpar3  = rnapol+hpol+exp(-hnaa1*r(i,3)) *(hnaa2+hnaa3*r(i,3)+hnaa4*r3s+hnaa5*r3s*r(i,3))
            rnahperp1 = rnapol+hpol+exp(-hnab1*r(i,1))* (hnab2+hnab3*r(i,1)+hnab4*r1s+hnab5*r1s*r(i,1))
            rnahperp3 = rnapol+hpol+exp(-hnab1*r(i,3))* (hnab2+hnab3*r(i,3)+hnab4*r3s+hnab5*r3s*r(i,3))
            rnahpol1  = rnahperp1+(rnahpar1-rnahperp1)*cstsq1
            rnahpol1b = onethd*rnahpar1+twothd*rnahperp1
            rnahpol3  = rnahperp3+(rnahpar3-rnahperp3)*cstsq3
            rnahpol3b = onethd*rnahpar3+twothd*rnahperp3
            ! induction forces
            ! dipole induced dipole
            div1    = 1.0e0_wp/(rnag1c*rnag1c*rnag1c+rd)
            rmuimu1 = -rnahmu1**2*hpol*(1.5e0_wp*cstsq1+0.5e0_wp)*rnag1c*div1
            div3    = 1.0e0_wp/(rnag3c*rnag3c*rnag3c+rd)
            rmuimu3 = -rnahmu3**2*hpol*(1.5e0_wp*cstsq3+0.5e0_wp)*rnag3c*div3
            ! dispersion forces
            c690   = -1.5e0_wp*rnaie*h2ie*rnapol*al2 /((rnaie+h2ie)*(rnagc*rnagc*rnagc+rb90))
            c690   = c690*f90*rnagc
            c60    = -1.5e0_wp*rnaie*h2ie*rnapol*al1/ ((rnaie+h2ie)*(rnagc*rnagc*rnagc+rb0))
            c60    = c60*f0*rnagc
            c6     = c690+(c60-c690)*cstsq*(1.0e0_wp-exp(-39.0625e0_wp*rnag*rnagc))
            c61    = -1.5e0_wp*rnahie1*hie*rnag1c/(rnahie1+hie)*rnahpol1*hpol*div1
            dexp1  = exp(-39.0625e0_wp*rnag1c*rnag1)
            c61    = c61*(1.0e0_wp+(rnahpol1b/rnahpol1-1.0e0_wp)*dexp1)
            c63    = -1.5e0_wp*rnahie3*hie*rnag3c/(rnahie3+hie)*rnahpol3*hpol*div3
            dexp3  = exp(-39.0625e0_wp*rnag3c*rnag3)
            c63    = c63*(1.0e0_wp+(rnahpol3b/rnahpol3-1.0e0_wp)*dexp3)
            cind1  = rmuimu1*(1.0e0_wp+(1.0e0_wp/(1.5e0_wp*cstsq1+0.5e0_wp)-1.0e0_wp)*dexp1)
            cind3  = rmuimu3*(1.0e0_wp+(1.0e0_wp/(1.5e0_wp*cstsq3+0.5e0_wp)-1.0e0_wp)*dexp3)
            c6     = c6*cevau*0.5e0_wp*(1.0e0_wp-tanh(balp*(r(i,2)-brho)))
            c6swt1 = 0.5e0_wp*(1.0e0_wp-tanh(dalp*(r(i,1)-drho)))
            c6swt3 = 0.5e0_wp*(1.0e0_wp-tanh(dalp*(r(i,3)-drho)))
            c61    = c61*c6swt1*cevau
            c63    = c63*c6swt3*cevau
            cind1  = cind1*c6swt1
            cind3  = cind3*c6swt3
            elec   = c61+c63+cind1+cind3+c6
            ! here is where it all gets put together . . .
            e(i) = s1 + s2 + s3 + deh2
            elec = elec/cevau
            e(i) = e(i)+elec
            ! here we will try a post-LEPS correction term
            xect = exp(-bect*(rmid-rect))
            vmid = dect*xect*(xect-2.0e0_wp) *exp(-bect2*rave**2)*(0.5e0_wp*(1.0e0_wp-cosa))**2
            e(i) = e(i)+vmid
            e(i) = e(i)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_7l_2000_pot11_legacy_state_c810f3bd

!> @brief Compatibility interface for pot11.
subroutine pot11( r, e, nt )
    use chempotpy_nah2_leps_7l_2000_pot11_legacy_state_c810f3bd, only: dispatch
    use chempotpy_kinds,                                         only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot11

! U12 coupling surface . . .

!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_nah2_leps_7l_2000_pot12_legacy_state_c810f3bd
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r, e, nt
        real(wp) :: ap, alp, bet, rnagsq, rnag, cosg, rd0, r20
        real(wp) :: r1s, r2s, r3s, cos1, cos2, theta
        real(wp) :: gc, gc2
        real(wp) :: cevau
        integer  :: nt, i, j, k
        logical  :: leven
        real(wp) :: r(:,:), e(:)
        alp    = 0.e0_wp
        ap     = 0.e0_wp
        bet    = 0.e0_wp
        cevau  = 0.e0_wp
        cos1   = 0.e0_wp
        cos2   = 0.e0_wp
        cosg   = 0.e0_wp
        gc     = 0.e0_wp
        gc2    = 0.e0_wp
        r1s    = 0.e0_wp
        r20    = 0.e0_wp
        r2s    = 0.e0_wp
        r3s    = 0.e0_wp
        rd0    = 0.e0_wp
        rnag   = 0.e0_wp
        rnagsq = 0.e0_wp
        theta  = 0.e0_wp
        ap     = 0.835e0_wp
        alp    = 4.0e0_wp
        rd0    = 4.0e0_wp
        bet    = 0.5e0_wp
        r20    = -1.3e0_wp
        theta  = 15.6e0_wp*acos(-1.0e0_wp)/180.0e0_wp
        cevau  = 1.e0_wp/27.2113961e0_wp
        ! set leven true for an even U12 surface, and false for an odd
        ! surface
        leven = .true.
        do i=1, nt
            r1s = r(i,1)**2
            r2s = r(i,2)**2
            r3s = r(i,3)**2
            ! here we do the function evaluations . . .
            rnagsq = 0.5e0_wp*r1s+0.5e0_wp*r3s-0.25e0_wp*r2s
            if( rnagsq .gt. 0.0e0_wp )then
                rnag = sqrt(rnagsq)
            else
                rnag = 0.0e0_wp
            endif
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            cos1 = (r1s+r2s-r3s)/(2.0e0_wp*r(i,1)*r(i,2))
            cos2 = (r3s+r2s-r1s)/(2.0e0_wp*r(i,3)*r(i,2))
            cosg = 0.5e0_wp*(cos1-cos2)
            gc   = cos(theta)*rnag + sin(theta)*r(i,2)
            gc2  = sin(theta)*rnag - cos(theta)*r(i,2)
            e(i) = ap*cosg *exp(-(1.0e0_wp/alp)**4*(gc - rd0)**4) *exp(-(1.0e0_wp/bet)**2*(gc2 - r20)**2)
            e(i) = e(i)*cevau
            ! set leven equal to true for an even surface, false for an odd one
            if( leven )then
                e(i) = sqrt(e(i)**4/(e(i)**2+1.0e-8_wp))
            endif
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_7l_2000_pot12_legacy_state_c810f3bd

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, nt )
    use chempotpy_nah2_leps_7l_2000_pot12_legacy_state_c810f3bd, only: dispatch
    use chempotpy_kinds,                                         only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot12

! the upper diabatic surface  U22

!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_nah2_leps_7l_2000_pot22_legacy_state_c810f3bd
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: dispatch

    contains

    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        use chempotpy_kinds, only: wp
        implicit none
        optional :: r, e, nt
        real(wp) :: al1, al2, alb, alhh, alt, b0, balp, beh2, bet1, bet10, bet10b, bet190, bet1c, bet20, bet290, beta, &
            & bf, brho, c2a, c2b, c2c, c6, c61, c63, c6swt1, c6swt3, cevau, cind1, cind3, cos1, cos2, cosa, coscos, cosg, coul1, &
            & coul2, coul3, cplg2, cstsq, cstsq1, cstsq3, dalp, deh2, deh2m, denah, des290, det1, det190, det1c, det20, det290, &
            & dexp1, dexp2, dexp3, div1, div3, drho, e(:), elec, exch1, exch2, exch3, f10, frac, g10, g30, gamh2a, gamh2b, &
            & gamnah, h2a1, h2a2, h2a3, h2a4, h2a5, h2a6, h2a7, h2a8, h2a9, h2b1, h2b2, h2b3, h2b4, h2b5, h2b6, h2b7, h2b8, &
            & h2b9, h2hc1, h2hc2, h2hc3, h2hc4, h2hc5, h2hex, h2i1, h2i2, h2i3, h2i4, h2i5, h2i6, h2ie, h2qc1, h2qc2, h2qc3, &
            & h2qc4, h2qc5, h2quad, hie, hnaa1, hnaa2, hnaa3, hnaa4, hnaa5, hnab1, hnab2, hnab3, hnab4, hnab5, hnai1, hnai2, &
            & hnai3, hnai4, hnai5, hnai6, hnau1, hnau2, hnau3, hnau4, hnau5, hnau6, hpol, onethd, pi, qh, qhint, qq, qqint, &
            & r(:,:), r1s, r2s, r3s, rave, rb, rd, reh2, renah, ret1, ret10, ret190, ret1c, ret20, ret290, rmid, rmuimu1, &
            & rmuimu3, rnag, rnag1, rnag1c, rnag3, rnag3c, rnagc, rnags, rnahie1, rnahie3, rnahmu1, rnahmu3, rnahpar1, rnahpar3, &
            & rnahperp1, rnahperp3, rnahpol1, rnahpol1b, rnahpol3, rnahpol3b, rnaie, rnapol1, rnaquad, rqh, rqq, s1, s2, s20alp, &
            & s20b, s20rho, s2a, s2alp, s2rho, s3, sswt, t1, t10, t190, t1alp, t1c, t1m, t1rho, t2, t20, t290, t2alp, t2rho, t3, &
            & t30, t390, t3c, t3m, twothd, una2p, w, xeh2, xenah, xet190, xet1c, xet20, xet290, xet3, xet390, xet3c
        integer :: i, nt
        ! NaH triplets
        ! H2 triplets
        ! H2 singlets
        ! misc constants
        ! data blocks for the long range forces
        al1       = 0.e0_wp
        al2       = 0.e0_wp
        alb       = 0.e0_wp
        alhh      = 0.e0_wp
        alt       = 0.e0_wp
        b0        = 0.e0_wp
        balp      = 0.e0_wp
        beh2      = 0.e0_wp
        bet1      = 0.e0_wp
        bet10     = 0.e0_wp
        bet10b    = 0.e0_wp
        bet190    = 0.e0_wp
        bet1c     = 0.e0_wp
        bet20     = 0.e0_wp
        bet290    = 0.e0_wp
        beta      = 0.e0_wp
        bf        = 0.e0_wp
        brho      = 0.e0_wp
        c2a       = 0.e0_wp
        c2b       = 0.e0_wp
        c2c       = 0.e0_wp
        c6        = 0.e0_wp
        c61       = 0.e0_wp
        c63       = 0.e0_wp
        c6swt1    = 0.e0_wp
        c6swt3    = 0.e0_wp
        cevau     = 0.e0_wp
        cind1     = 0.e0_wp
        cind3     = 0.e0_wp
        cos1      = 0.e0_wp
        cos2      = 0.e0_wp
        cosa      = 0.e0_wp
        coscos    = 0.e0_wp
        cosg      = 0.e0_wp
        coul1     = 0.e0_wp
        coul2     = 0.e0_wp
        coul3     = 0.e0_wp
        cplg2     = 0.e0_wp
        cstsq     = 0.e0_wp
        cstsq1    = 0.e0_wp
        cstsq3    = 0.e0_wp
        dalp      = 0.e0_wp
        deh2      = 0.e0_wp
        deh2m     = 0.e0_wp
        denah     = 0.e0_wp
        des290    = 0.e0_wp
        det1      = 0.e0_wp
        det190    = 0.e0_wp
        det1c     = 0.e0_wp
        det20     = 0.e0_wp
        det290    = 0.e0_wp
        dexp1     = 0.e0_wp
        dexp2     = 0.e0_wp
        dexp3     = 0.e0_wp
        div1      = 0.e0_wp
        div3      = 0.e0_wp
        drho      = 0.e0_wp
        elec      = 0.e0_wp
        exch1     = 0.e0_wp
        exch2     = 0.e0_wp
        exch3     = 0.e0_wp
        f10       = 0.e0_wp
        frac      = 0.e0_wp
        g10       = 0.e0_wp
        g30       = 0.e0_wp
        gamh2a    = 0.e0_wp
        gamh2b    = 0.e0_wp
        gamnah    = 0.e0_wp
        h2a1      = 0.e0_wp
        h2a2      = 0.e0_wp
        h2a3      = 0.e0_wp
        h2a4      = 0.e0_wp
        h2a5      = 0.e0_wp
        h2a6      = 0.e0_wp
        h2a7      = 0.e0_wp
        h2a8      = 0.e0_wp
        h2a9      = 0.e0_wp
        h2b1      = 0.e0_wp
        h2b2      = 0.e0_wp
        h2b3      = 0.e0_wp
        h2b4      = 0.e0_wp
        h2b5      = 0.e0_wp
        h2b6      = 0.e0_wp
        h2b7      = 0.e0_wp
        h2b8      = 0.e0_wp
        h2b9      = 0.e0_wp
        h2hc1     = 0.e0_wp
        h2hc2     = 0.e0_wp
        h2hc3     = 0.e0_wp
        h2hc4     = 0.e0_wp
        h2hc5     = 0.e0_wp
        h2hex     = 0.e0_wp
        h2i1      = 0.e0_wp
        h2i2      = 0.e0_wp
        h2i3      = 0.e0_wp
        h2i4      = 0.e0_wp
        h2i5      = 0.e0_wp
        h2i6      = 0.e0_wp
        h2ie      = 0.e0_wp
        h2qc1     = 0.e0_wp
        h2qc2     = 0.e0_wp
        h2qc3     = 0.e0_wp
        h2qc4     = 0.e0_wp
        h2qc5     = 0.e0_wp
        h2quad    = 0.e0_wp
        hie       = 0.e0_wp
        hnaa1     = 0.e0_wp
        hnaa2     = 0.e0_wp
        hnaa3     = 0.e0_wp
        hnaa4     = 0.e0_wp
        hnaa5     = 0.e0_wp
        hnab1     = 0.e0_wp
        hnab2     = 0.e0_wp
        hnab3     = 0.e0_wp
        hnab4     = 0.e0_wp
        hnab5     = 0.e0_wp
        hnai1     = 0.e0_wp
        hnai2     = 0.e0_wp
        hnai3     = 0.e0_wp
        hnai4     = 0.e0_wp
        hnai5     = 0.e0_wp
        hnai6     = 0.e0_wp
        hnau1     = 0.e0_wp
        hnau2     = 0.e0_wp
        hnau3     = 0.e0_wp
        hnau4     = 0.e0_wp
        hnau5     = 0.e0_wp
        hnau6     = 0.e0_wp
        hpol      = 0.e0_wp
        onethd    = 0.e0_wp
        pi        = 0.e0_wp
        qh        = 0.e0_wp
        qhint     = 0.e0_wp
        qq        = 0.e0_wp
        qqint     = 0.e0_wp
        r1s       = 0.e0_wp
        r2s       = 0.e0_wp
        r3s       = 0.e0_wp
        rave      = 0.e0_wp
        rb        = 0.e0_wp
        rd        = 0.e0_wp
        reh2      = 0.e0_wp
        renah     = 0.e0_wp
        ret1      = 0.e0_wp
        ret10     = 0.e0_wp
        ret190    = 0.e0_wp
        ret1c     = 0.e0_wp
        ret20     = 0.e0_wp
        ret290    = 0.e0_wp
        rmid      = 0.e0_wp
        rmuimu1   = 0.e0_wp
        rmuimu3   = 0.e0_wp
        rnag      = 0.e0_wp
        rnag1     = 0.e0_wp
        rnag1c    = 0.e0_wp
        rnag3     = 0.e0_wp
        rnag3c    = 0.e0_wp
        rnagc     = 0.e0_wp
        rnags     = 0.e0_wp
        rnahie1   = 0.e0_wp
        rnahie3   = 0.e0_wp
        rnahmu1   = 0.e0_wp
        rnahmu3   = 0.e0_wp
        rnahpar1  = 0.e0_wp
        rnahpar3  = 0.e0_wp
        rnahperp1 = 0.e0_wp
        rnahperp3 = 0.e0_wp
        rnahpol1  = 0.e0_wp
        rnahpol1b = 0.e0_wp
        rnahpol3  = 0.e0_wp
        rnahpol3b = 0.e0_wp
        rnaie     = 0.e0_wp
        rnapol1   = 0.e0_wp
        rnaquad   = 0.e0_wp
        rqh       = 0.e0_wp
        rqq       = 0.e0_wp
        s1        = 0.e0_wp
        s2        = 0.e0_wp
        s20alp    = 0.e0_wp
        s20b      = 0.e0_wp
        s20rho    = 0.e0_wp
        s2a       = 0.e0_wp
        s2alp     = 0.e0_wp
        s2rho     = 0.e0_wp
        s3        = 0.e0_wp
        sswt      = 0.e0_wp
        t1        = 0.e0_wp
        t10       = 0.e0_wp
        t190      = 0.e0_wp
        t1alp     = 0.e0_wp
        t1c       = 0.e0_wp
        t1m       = 0.e0_wp
        t1rho     = 0.e0_wp
        t2        = 0.e0_wp
        t20       = 0.e0_wp
        t290      = 0.e0_wp
        t2alp     = 0.e0_wp
        t2rho     = 0.e0_wp
        t3        = 0.e0_wp
        t30       = 0.e0_wp
        t390      = 0.e0_wp
        t3c       = 0.e0_wp
        t3m       = 0.e0_wp
        twothd    = 0.e0_wp
        una2p     = 0.e0_wp
        w         = 0.e0_wp
        xeh2      = 0.e0_wp
        xenah     = 0.e0_wp
        xet190    = 0.e0_wp
        xet1c     = 0.e0_wp
        xet20     = 0.e0_wp
        xet290    = 0.e0_wp
        xet3      = 0.e0_wp
        xet390    = 0.e0_wp
        xet3c     = 0.e0_wp
        ret10     = 3.22835e0_wp
        bet10     = 0.57874e0_wp
        bet10b    = 2.44094e0_wp
        f10       = 9.44882e0_wp
        det190    = 4.32441e0_wp
        bet190    = 2.43701e0_wp
        ret190    = 1.95480e0_wp
        t1alp     = 0.30778e0_wp
        t1rho     = 3.21850e0_wp
        det1c     = 4.92126e0_wp
        bet1c     = 0.65748e0_wp
        ret1c     = 2.71654e0_wp
        det20     = 3.29921e0_wp
        bet20     = 2.64567e0_wp
        ret20     = 2.44882e0_wp
        det290    = 3.97638e0_wp
        bet290    = 3.73937e0_wp
        ret290    = 0.96024e0_wp
        t2alp     = 1.0e0_wp
        t2rho     = 9.0e0_wp
        des290    = 0.91433e0_wp
        s2alp     = 0.66778e0_wp
        s2rho     = 5.30472e0_wp
        deh2      = 4.74806e0_wp
        s20alp    = 0.5e0_wp
        s20rho    = 6.0e0_wp
        c2a       = 1.38661e0_wp
        c2b       = 0.27362e0_wp
        c2c       = 0.03764e0_wp
        una2p     = 2.1037e0_wp
        h2a1      = 2.38029900e0_wp
        h2a2      = 1.47400334e0_wp
        h2a3      = -0.226331771e0_wp
        h2a4      = -2.17339562e0_wp
        h2a5      = 9.54597998e0_wp
        h2a6      = -23.3366749e0_wp
        h2a7      = 24.7115748e0_wp
        h2a8      = -12.8786102e0_wp
        h2a9      = 2.75711758e0_wp
        h2b1      = 1.35164094e0_wp
        h2b2      = -0.254404564e0_wp
        h2b3      = -10.6627055e0_wp
        h2b4      = -14.0112406e0_wp
        h2b5      = -7.80660398e0_wp
        h2b6      = 5.47697753e0_wp
        h2b7      = -6.95875664e0_wp
        h2b8      = 3.64334771e0_wp
        h2b9      = -0.492505596e0_wp
        h2i1      = 160.484e0_wp
        h2i2      = -20.397e0_wp
        h2i3      = 3.913e0_wp
        h2i4      = 0.777e0_wp
        h2i5      = -7.583
        h2i6      = 1.501
        rnaie     = 3.0353
        rnapol1   = 290.4e0_wp
        hpol      = 4.5e0_wp
        hie       = 13.598e0_wp
        rnaquad   = -16.8992e0_wp
        hnai1     = -193.341e0_wp
        hnai2     = 16.968e0_wp
        hnai3     = -77.120e0_wp
        hnai4     = -1.024e0_wp
        hnai5     = 6.009e0_wp
        hnai6     = -1.241e0_wp
        hnau1     = 119.880e0_wp
        hnau2     = -0.620e0_wp
        hnau3     = 31.335e0_wp
        hnau4     = 70.414e0_wp
        hnau5     = -13.016e0_wp
        hnau6     = 0.8681e0_wp
        hnaa1     = 0.67e0_wp
        hnaa2     = -92.229e0_wp
        hnaa3     = 51.01019e0_wp
        hnaa4     = -130.37826e0_wp
        hnaa5     = 6.40285e0_wp
        hnab1     = 0.34333e0_wp
        hnab2     = -92.229e0_wp
        hnab3     = -48.54867e0_wp
        hnab4     = -17.30523e0_wp
        hnab5     = 0.99821e0_wp
        rqq       = 9.766021100e4_wp
        rqh       = 2.355466062927e6_wp
        rb        = 9.3119383371714e5_wp
        balp      = 1.93701e0_wp
        brho      = 1.4e0_wp
        rd        = 9.1351724748364e8_wp
        dalp      = 1.43307e0_wp
        drho      = 5.71654e0_wp
        h2hc1     = 0.677e0_wp
        h2hc2     = 0.378e0_wp
        h2hc3     = 0.281e0_wp
        h2hc4     = -1.022e0_wp
        h2hc5     = 0.665e0_wp
        h2qc1     = 1.532e0_wp
        h2qc2     = -.04635e0_wp
        h2qc3     = 1.826e0_wp
        h2qc4     = -2.508e0_wp
        h2qc5     = 1.856e0_wp
        cevau     = 1.e0_wp/27.2113961e0_wp
        onethd    = 1.0e0_wp/3.0e0_wp
        twothd    = 2.0e0_wp/3.0e0_wp
        do i=1, nt
            ! this section calcluates useful quantities - the A to BC center
            ! of mass distance, the angles between A-BC and BC, the shorest
            ! of any of the three diatomic distances, and the bond angle
            ! H-Na-H
            rmid  = 0.5e0_wp*(r(i,1)+r(i,3))
            rave  = 0.5e0_wp*(r(i,1)-r(i,3))
            r1s   = r(i,1)*r(i,1)
            r2s   = r(i,2)*r(i,2)
            r3s   = r(i,3)*r(i,3)
            rnag  = sqrt(0.5e0_wp*abs(r1s+r3s-0.5e0_wp*r2s))
            rnags = rnag*rnag
            rnagc = rnags*rnag
            if( rnag .eq. 0.0e0_wp )then
                cstsq = 0.0e0_wp
            else
                cstsq = (0.5e0_wp*(r3s-r1s)/(r(i,2)*rnag))**2
            endif
            cosa   = (r1s+r3s-r2s)/(2*r(i,1)*r(i,3))
            coscos = cstsq*(0.5e0_wp+0.5e0_wp*cosa)*(1.0e0_wp-exp(-10.0e0_wp*rnags))
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            cos1   = (r1s+r2s-r3s)/(2.0e0_wp*r(i,1)*r(i,2))
            cos2   = (r3s+r2s-r1s)/(2.0e0_wp*r(i,3)*r(i,2))
            cosg   = 0.5e0_wp*(cos1-cos2)
            frac   = 0.958e0_wp
            rnag1  = sqrt(abs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))
            rnag3  = sqrt(abs(r2s+frac*(r1s-r2s-r3s)+r3s*frac*frac))
            rnag1c = rnag1*rnag1*rnag1
            rnag3c = rnag3*rnag3*rnag3
            if( rnag1 .eq. 0.0e0_wp )then
                cstsq1 = 1.0e0_wp
            else
                cstsq1 = 0.25e0_wp*((r3s-r2s-r1s+2.0e0_wp*r3s*frac)/(r(i,1)*rnag1)) **2
            endif
            if( rnag3 .eq. 0.0e0_wp )then
                cstsq3 = 1.0e0_wp
            else
                cstsq3 = 0.25e0_wp*((r1s-r2s-r3s+2.0e0_wp*r1s*frac)/(r(i,3)*rnag3)) **2
            endif
            ! The first NaH singlet
            ! new singlet, designed to minic lower adiabatic behavior
            denah  = 1.97134878e0_wp
            renah  = 3.566044e0_wp
            bf     = 0.864e0_wp
            b0     = 0.594e0_wp
            gamnah = 7.376e0_wp
            beta   = bf*(b0+(r(i,1)/gamnah)**8)/(bf+(r(i,1)/gamnah)**8)
            xenah  = exp(-beta*(r(i,1)-renah))
            s1     = denah*xenah*(xenah-2.0e0_wp)
            ! The first NaH triplet
            g10    = r(i,1)-ret10
            t10    = f10*exp(-bet10*g10)+exp(-bet10b*g10)
            t10    = t10/(1.0e0_wp+f10)
            xet190 = exp(-bet190*(r(i,1)-ret190))
            t190   = det190*onethd*xet190*(xet190 + 2.0e0_wp)
            xet1c  = exp(-bet1c*(r(i,1)-ret1c))
            t1c    = det1c*onethd*xet1c*(xet1c + 2.0e0_wp)
            t1m    = t190 +(t1c-t190)* 0.5e0_wp*(1.0e0_wp+tanh(t1alp*(r(i,2)-t1rho)))
            t1     = t1m + (t10-t1m)*cosg**2
            ! The H2 triplet
            xet20  = exp(-bet20*(r(i,2)-ret20))
            t20    = det20*onethd*xet20*(xet20 + 2.0e0_wp)
            xet290 = exp(-bet290*(r(i,2)-ret290))
            t290   = det290*onethd*xet290*(xet290 + 2.0e0_wp)
            t2     = t20 + (t290-t20)*(1.0e0_wp-cosg**2)* 0.5e0_wp*(1.0e0_wp-tanh(t2alp*(rmid-t2rho)))
            ! The H2 singlet
            ! new singlet, designed to minic lower adiabatic behavior
            reh2   = 1.401e0_wp
            bf     = 1.627e0_wp
            b0     = 1.194e0_wp
            gamh2a = 2.323e0_wp
            gamh2b = 3.126e0_wp
            beh2   = bf*(b0-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2) /(bf-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2)
            sswt   = 0.5e0_wp*(1.0e0_wp-tanh(s2alp*(rmid-s2rho)))
            deh2m  = deh2 + deh2*(des290 - 1.0e0_wp)*sswt
            xeh2   = exp(-beh2*(r(i,2)-reh2))
            s2a    = deh2m*xeh2*(xeh2-2.0e0_wp) + una2p
            s2     = 0.5e0_wp*(s2a+t2-tanh((s2a-t2)/0.1e0_wp)*(s2a-t2))
            xeh2   = exp(-1.34e0_wp*beh2*(r(i,2)-reh2))
            s20b   = (deh2-una2p)*xeh2*(xeh2-2.0e0_wp)
            s2     = s2 + (s20b - s2)*cosg**2* 0.5e0_wp*(1.0e0_wp-tanh(s20alp*(rmid-s20rho)))
            ! The second NaH singlet
            ! new singlet, designed to minic lower adiabatic behavior
            denah  = 1.97134878e0_wp
            renah  = 3.566044e0_wp
            bf     = 0.864e0_wp
            b0     = 0.594e0_wp
            gamnah = 7.376e0_wp
            beta   = bf*(b0+(r(i,3)/gamnah)**8)/(bf+(r(i,3)/gamnah)**8)
            xenah  = exp(-beta*(r(i,3)-renah))
            s3     = denah*xenah*(xenah-2.0e0_wp)
            ! The second NaH triplet
            xet3   = exp(-bet1*(r(i,3)-ret1))
            t30    = det1*onethd*xet3*(xet3 + 2.0e0_wp)
            g30    = r(i,3)-ret10
            t30    = f10*exp(-bet10*g30)+exp(-bet10b*g30)
            t30    = t30/(1.0e0_wp+f10)
            xet390 = exp(-bet190*(r(i,3)-ret190))
            t390   = det190*onethd*xet390*(xet390 + 2.0e0_wp)
            xet3c  = exp(-bet1c*(r(i,3)-ret1c))
            t3c    = det1c*onethd*xet3c*(xet3c + 2.0e0_wp)
            t3m    = t390 + (t3c - t390)* 0.5e0_wp*(1.0e0_wp+tanh(t1alp*(r(i,2)-t1rho)))
            t3     = t3m + (t30 - t3m)*cosg**2
            ! this long section contains data pertinent to long range forces
            !  H2 ionization energy:
            h2ie = (h2i1+h2i2*r(i,2))*exp(-h2i3*r(i,2))+ r2s*(h2i4+h2i5*r(i,2))*exp(-h2i6*r(i,2))-s2a+hie
            !  H2 polarizability:
            !   al1 is the parallel polarizability for h2
            !   al2 is the perpendicaular polarizability for h2
            !   alt is the anisotropic? polarizability
            !   alb is the mean polarizability
            al1=exp(h2a1*(h2a2-r(i,2)))*(h2a3+h2a4*r(i,2)+h2a5*r2s +h2a6*r(i,2)*r2s+h2a7*r2s*r2s+h2a8*r2s*r2s* r(i, &
                & 2) + h2a9* r2s*r2s*r2s) +2.0e0_wp*hpol
            al2=exp(h2b1*(h2b2-r(i,2)))*(h2b3+h2b4*r(i,2)+h2b5*r2s +h2b6*r(i,2)**3+h2b7*r2s*r2s+h2b8*r2s*r2s*r(i,2) + &
                & h2b9* r2s*r2s*r2s) +2.0e0_wp*hpol
            alt  = twothd*(al1-al2)
            alb  = onethd*al1+twothd*al2
            alhh = alb+alt*(1.5e0_wp*cstsq-0.5e0_wp)
            ! here follow properties necesary for the long range forces in
            ! the Na + H2 channel:
            ! H2 properties:
            h2hex  = exp(-h2hc1*r(i,2)) *(h2hc2*r(i,2)+h2hc3*r2s+h2hc4*r(i,2)*r2s+h2hc5*r2s*r2s)
            h2quad = exp(-h2qc1*r(i,2)) *(h2qc2*r(i,2)+h2qc3*r2s+h2qc4*r(i,2)*r2s+h2qc5*r2s*r2s)
            ! NaH properties
            rnahie1   = (hnai1+hnai2*r(i,1)) /(hnai3+hnai4*r(i,1)+hnai5*r1s+hnai6*r(i,1)*r1s)+rnaie+una2p
            rnahie3   = (hnai1+hnai2*r(i,3)) /(hnai3+hnai4*r(i,3)+hnai5*r3s+hnai6*r(i,3)*r3s)+rnaie+una2p
            rnahmu1   = (hnau1*r(i,1)+hnau2*r1s) /(hnau3+hnau4*r(i,1)+hnau5*r1s+hnau6*r(i,1)*r1s)
            rnahmu3   = (hnau1*r(i,3)+hnau2*r3s) /(hnau3+hnau4*r(i,3)+hnau5*r3s+hnau6*r(i,3)*r3s)
            rnahpar1  = 159.267e0_wp+hpol+exp(-hnaa1*r(i,1)) *(hnaa2+hnaa3*r(i,1)+hnaa4*r1s+hnaa5*r1s*r(i,1))
            rnahpar3  = 159.267e0_wp+hpol+exp(-hnaa1*r(i,3)) *(hnaa2+hnaa3*r(i,3)+hnaa4*r3s+hnaa5*r3s*r(i,3))
            rnahperp1 = 159.267e0_wp+hpol+exp(-hnab1*r(i,1))* (hnab2+hnab3*r(i,1)+hnab4*r1s+hnab5*r1s*r(i,1))
            rnahperp3 = 159.267e0_wp+hpol+exp(-hnab1*r(i,3))* (hnab2+hnab3*r(i,3)+hnab4*r3s+hnab5*r3s*r(i,3))
            rnahpol1  = rnahperp1+(rnahpar1-rnahperp1)*cstsq1
            rnahpol1b = onethd*rnahpar1+twothd*rnahperp1
            rnahpol3  = rnahperp3+(rnahpar3-rnahperp3)*cstsq3
            rnahpol3b = onethd*rnahpar3+twothd*rnahperp3
            ! electrostatic forces
            ! quadrupole quadrupole
            qqint = 3.0e0_wp-7.0e0_wp*cstsq
            qq    = 0.75e0_wp*rnaquad*h2quad/(rnagc*rnag*rnag+rqq)*qqint
            dexp2 = exp(-39.0625e0_wp*rnagc*rnag)
            qq    = qq*(1.0e0_wp+(-4.0e0_wp*pi/(3.0e0_wp*qqint)-1.0e0_wp)*dexp2)
            ! quadrupole hexadecapole
            qhint = -45.0e0_wp-147.0e0_wp*cstsq*cstsq+162.0e0_wp*cstsq
            qh    = 0.3125e0_wp*rnaquad*h2hex/(rnagc*rnagc*rnag+rqh)*(qhint)
            qh    = qh*(1.0e0_wp+(-40.8e0_wp*pi/qhint-1.0e0_wp)*dexp2)
            ! induction forces
            ! dipole induced dipole
            div1    = 1.0e0_wp/(rnag1c*rnag1c*rnag1c+rd)
            rmuimu1 = -rnahmu1**2*hpol*(1.5e0_wp*cstsq1+0.5e0_wp)*rnag1c*div1
            div3    = 1.0e0_wp/(rnag3c*rnag3c*rnag3c+rd)
            rmuimu3 = -rnahmu3**2*hpol*(1.5e0_wp*cstsq3+0.5e0_wp)*rnag3c*div3
            ! dispersion forces
            c6     = -1.5e0_wp*rnaie*h2ie/(rnaie+h2ie)*rnapol1*alhh/(rnagc*rnagc+rb)
            c6     = c6*(1.0e0_wp+(alb/alhh-1.0e0_wp)*dexp2)
            c61    = -1.5e0_wp*rnahie1*hie/(rnahie1+hie)*rnahpol1*hpol*rnag1c*div1
            dexp1  = exp(-39.0625e0_wp*rnag1c*rnag1)
            c61    = c61*(1.0e0_wp+(rnahpol1b/rnahpol1-1.0e0_wp)*dexp1)
            c63    = -1.5e0_wp*rnahie3*hie/(rnahie3+hie)*rnahpol3*hpol*rnag3c*div3
            dexp3  = exp(-39.0625e0_wp*rnag3c*rnag3)
            c63    = c63*(1.0e0_wp+(rnahpol3b/rnahpol3-1.0e0_wp)*dexp3)
            cind1  = rmuimu1*(1.0e0_wp+(1.0e0_wp/(1.5e0_wp*cstsq1+0.5e0_wp)-1.0e0_wp)*dexp1)
            cind3  = rmuimu3*(1.0e0_wp+(1.0e0_wp/(1.5e0_wp*cstsq3+0.5e0_wp)-1.0e0_wp)*dexp3)
            elec   = (qq+qh+c6*cevau)*0.5e0_wp *(1.0e0_wp-tanh(balp*(r(i,2)-brho)))
            c6swt1 = 0.5e0_wp*(1.0e0_wp-tanh(dalp*(r(i,1)-drho)))
            c6swt3 = 0.5e0_wp*(1.0e0_wp-tanh(dalp*(r(i,3)-drho)))
            c61    = c61*c6swt1*cevau
            c63    = c63*c6swt3*cevau
            cind1  = cind1*c6swt1
            cind3  = cind3*c6swt3
            elec   = elec+c61+c63+cind1+cind3
            coul1  = 0.5e0_wp*(s1+t1)
            coul2  = 0.5e0_wp*(s2+t2)
            coul3  = 0.5e0_wp*(s3+t3)
            exch1  = 0.5e0_wp*(s1-t1)
            exch2  = 0.5e0_wp*(s2-t2)
            exch3  = 0.5e0_wp*(s3-t3)
            w      = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
            cplg2  = c2a*exp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
            e(i)   = coul1+coul2+coul3+deh2-sqrt(w+cplg2*cplg2)/sqrt(2.0e0_wp)
            elec   = elec/cevau
            e(i)   = (e(i)+elec)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_7l_2000_pot22_legacy_state_c810f3bd

!> @brief Compatibility interface for pot22.
subroutine pot22( r, e, nt )
    use chempotpy_nah2_leps_7l_2000_pot22_legacy_state_c810f3bd, only: dispatch
    use chempotpy_kinds,                                         only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot22
