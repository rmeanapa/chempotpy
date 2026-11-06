!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_nah2_leps_7l_2000_dpem_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, u, ug )
        use chempotpy_legacy_units,       only: legacy_hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: legacy_triatomic_distances_to_bohr
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 2
        integer, parameter            :: natoms  = 3
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: u(nstates,nstates)
        double precision, intent(out) :: ug(nstates,nstates,natoms,3)
        double precision              :: nt, r(1,3), v(1)
        double precision              :: u11(1), u12(1), u22(1)
        integer                       :: i, j, k, l
        !initialize
        u  = 0.d0; ug = 0.d0
        nt = 1
        ! input cartesian is ClHH
        call legacy_triatomic_distances_to_bohr(x, r(1,:))
        if( igrad == 0 )then
            call pot(r, u11, 1, 1)
            call pot(r, u12, 1, 2)
            call pot(r, u22, 1, 3)
            u(1,1) = u11(1)*legacy_hartree_to_electronvolt
            u(1,2) = u12(1)*legacy_hartree_to_electronvolt
            u(2,1) = u12(1)*legacy_hartree_to_electronvolt
            u(2,2) = u22(1)*legacy_hartree_to_electronvolt
        else
            write(*, *) "Only DPEM is available"
        endif
    end subroutine evaluate
end module chempotpy_nah2_leps_7l_2000_dpem_pes

!> @brief Preserve the historical Python-facing dpem entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_nah2_leps_7l_2000_dpem_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 2
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: u(nstates,nstates)
    double precision, intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem

!> @brief Evaluate EvdRdX for the NaH2_LEPS_7L_2000_DPEM potential.
subroutine evdrdx( x, r, drdx )
    implicit none
    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)
    integer                       :: i, j
    ! Initialize dRdX(3,9)
    do i=1, 3
        do j=1, 9
            drdx(i,j) = 0.0d0
        enddo
    enddo
    drdx(1,1) = (x(1)-x(4))/r(1)
    drdx(1,2) = (x(2)-x(5))/r(1)
    drdx(1,3) = (x(3)-x(6))/r(1)
    drdx(1,4) = -drdx(1,1)
    drdx(1,5) = -drdx(1,2)
    drdx(1,6) = -drdx(1,3)
    drdx(2,4) = (x(4)-x(7))/r(2)
    drdx(2,5) = (x(5)-x(8))/r(2)
    drdx(2,6) = (x(6)-x(9))/r(2)
    drdx(2,7) = -drdx(2,4)
    drdx(2,8) = -drdx(2,5)
    drdx(2,9) = -drdx(2,6)
    drdx(3,1) = (x(1)-x(7))/r(3)
    drdx(3,2) = (x(2)-x(8))/r(3)
    drdx(3,3) = (x(3)-x(9))/r(3)
    drdx(3,7) = -drdx(3,1)
    drdx(3,8) = -drdx(3,2)
    drdx(3,9) = -drdx(3,3)
end subroutine evdrdx

!> @brief Evaluate diagonalize for the NaH2_LEPS_7L_2000_DPEM potential.
subroutine diagonalize( n, a_ss, ev_s, u_ss )
    implicit none
    integer, intent(in)           :: n
    double precision, intent(in)  :: a_ss(n,n)
    double precision, intent(out) :: u_ss(n,n)
    double precision, intent(out) :: ev_s(n)
    integer                       :: io, i
    double precision, allocatable :: work_d(:)
    integer                       :: lwork_d
    u_ss    = a_ss
    lwork_d = -1
    allocate( work_d(1) )
    call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
    lwork_d = int(work_d(1))
    deallocate( work_d )
    allocate( work_d(lwork_d) )
    call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
end subroutine diagonalize

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
!
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!
!   Surfaces:
!      Ground diabatic electronic state, first-excited diabatic electron
!
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Na
!      infinitely far from the H2 diatomic and R(H2) set equal to the
!      H2 equilibrium diatomic value.
!
!   Parameters:
!
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

!

!
!> @brief Evaluate pot for the NaH2_LEPS_7L_2000_DPEM potential.
subroutine pot( r, e, nt, nsurf )
    implicit none
    double precision :: e, r
    integer          :: nsurf, nt
    dimension r(nt,3), e(nt)
    if( nsurf .eq. 1 )then
        call pot11(r, e, nt)
    else if( nsurf .eq. 2 )then
        call pot12(r, e, nt)
    else if( nsurf .eq. 3 )then
        call pot22(r, e, nt)
    endif
end subroutine pot

!  U11 surface - lower diabatic surface
!
!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_nah2_leps_7l_2000_dpem_pot11_legacy_state_f8f175af
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        implicit none
        optional :: r, e, nt
        double precision :: al1, al2, alb, alhh, alt, b0, balp, bect, bect2, beh2, bes1, bf, brho, c6, &
            & c60, c61, c63, c690, c6swt1, c6swt3, cevau, cind1, cind3, cos1, cos2, cosa, coscos, cosg, &
            & cstsq, cstsq1, cstsq3, dalp, dect, deh2, deh2m, des1, des20, dexp1, dexp3, div1, div3, &
            & drho, e, elec, f0, f90, frac, gamh2a, gamh2b, h2a1, h2a2, h2a3, h2a4, h2a5, h2a6, h2a7, &
            & h2a8, h2a9, h2b1, h2b2, h2b3, h2b4, h2b5, h2b6, h2b7, h2b8, h2b9, h2i1, h2i2, h2i3, h2i4, &
            & h2i5, h2i6, h2ie, hie, hnaa1, hnaa2, hnaa3, hnaa4, hnaa5, hnab1, hnab2, hnab3, hnab4, &
            & hnab5, hnai1, hnai2, hnai3, hnai4, hnai5, hnai6, hnau1, hnau2, hnau3, hnau4, hnau5, hnau6, &
            & hpol, onethd, r, r1s, r2s, r3s, rac2, rave, rb0, rb90, rd, rect, reh2, res1, rexp1, rexp2, &
            & rexp3, rmid, rmuimu1, rmuimu3, rnag, rnag1, rnag1c, rnag3, rnag3c, rnagc, rnags, rnahie1, &
            & rnahie3, rnahmu1, rnahmu3, rnahpar1, rnahpar3, rnahperp1, rnahperp3, rnahpol1, rnahpol1b, &
            & rnahpol3, rnahpol3b, rnaie, rnapol, rsm1, rsm2, rsm3, s1, s2, s2a, s2alp, s2rho, s3, sswt, &
            & twothd, vmid, xect, xeh2, xes1, xes3
        integer :: i, nt
        dimension r(:,:), e(:)
        ! NaH singlet
        ! H2 singlet
        ! H-Na-H term
        ! data blocks for the long range forces
        ! misc constants
        al1       = 0.d0
        al2       = 0.d0
        alb       = 0.d0
        alhh      = 0.d0
        alt       = 0.d0
        b0        = 0.d0
        balp      = 0.d0
        bect      = 0.d0
        bect2     = 0.d0
        beh2      = 0.d0
        bes1      = 0.d0
        bf        = 0.d0
        brho      = 0.d0
        c6        = 0.d0
        c60       = 0.d0
        c61       = 0.d0
        c63       = 0.d0
        c690      = 0.d0
        c6swt1    = 0.d0
        c6swt3    = 0.d0
        cevau     = 0.d0
        cind1     = 0.d0
        cind3     = 0.d0
        cos1      = 0.d0
        cos2      = 0.d0
        cosa      = 0.d0
        coscos    = 0.d0
        cosg      = 0.d0
        cstsq     = 0.d0
        cstsq1    = 0.d0
        cstsq3    = 0.d0
        dalp      = 0.d0
        dect      = 0.d0
        deh2      = 0.d0
        deh2m     = 0.d0
        des1      = 0.d0
        des20     = 0.d0
        dexp1     = 0.d0
        dexp3     = 0.d0
        div1      = 0.d0
        div3      = 0.d0
        drho      = 0.d0
        elec      = 0.d0
        f0        = 0.d0
        f90       = 0.d0
        frac      = 0.d0
        gamh2a    = 0.d0
        gamh2b    = 0.d0
        h2a1      = 0.d0
        h2a2      = 0.d0
        h2a3      = 0.d0
        h2a4      = 0.d0
        h2a5      = 0.d0
        h2a6      = 0.d0
        h2a7      = 0.d0
        h2a8      = 0.d0
        h2a9      = 0.d0
        h2b1      = 0.d0
        h2b2      = 0.d0
        h2b3      = 0.d0
        h2b4      = 0.d0
        h2b5      = 0.d0
        h2b6      = 0.d0
        h2b7      = 0.d0
        h2b8      = 0.d0
        h2b9      = 0.d0
        h2i1      = 0.d0
        h2i2      = 0.d0
        h2i3      = 0.d0
        h2i4      = 0.d0
        h2i5      = 0.d0
        h2i6      = 0.d0
        h2ie      = 0.d0
        hie       = 0.d0
        hnaa1     = 0.d0
        hnaa2     = 0.d0
        hnaa3     = 0.d0
        hnaa4     = 0.d0
        hnaa5     = 0.d0
        hnab1     = 0.d0
        hnab2     = 0.d0
        hnab3     = 0.d0
        hnab4     = 0.d0
        hnab5     = 0.d0
        hnai1     = 0.d0
        hnai2     = 0.d0
        hnai3     = 0.d0
        hnai4     = 0.d0
        hnai5     = 0.d0
        hnai6     = 0.d0
        hnau1     = 0.d0
        hnau2     = 0.d0
        hnau3     = 0.d0
        hnau4     = 0.d0
        hnau5     = 0.d0
        hnau6     = 0.d0
        hpol      = 0.d0
        onethd    = 0.d0
        r1s       = 0.d0
        r2s       = 0.d0
        r3s       = 0.d0
        rac2      = 0.d0
        rave      = 0.d0
        rb0       = 0.d0
        rb90      = 0.d0
        rd        = 0.d0
        rect      = 0.d0
        reh2      = 0.d0
        res1      = 0.d0
        rexp1     = 0.d0
        rexp2     = 0.d0
        rexp3     = 0.d0
        rmid      = 0.d0
        rmuimu1   = 0.d0
        rmuimu3   = 0.d0
        rnag      = 0.d0
        rnag1     = 0.d0
        rnag1c    = 0.d0
        rnag3     = 0.d0
        rnag3c    = 0.d0
        rnagc     = 0.d0
        rnags     = 0.d0
        rnahie1   = 0.d0
        rnahie3   = 0.d0
        rnahmu1   = 0.d0
        rnahmu3   = 0.d0
        rnahpar1  = 0.d0
        rnahpar3  = 0.d0
        rnahperp1 = 0.d0
        rnahperp3 = 0.d0
        rnahpol1  = 0.d0
        rnahpol1b = 0.d0
        rnahpol3  = 0.d0
        rnahpol3b = 0.d0
        rnaie     = 0.d0
        rnapol    = 0.d0
        rsm1      = 0.d0
        rsm2      = 0.d0
        rsm3      = 0.d0
        s1        = 0.d0
        s2        = 0.d0
        s2a       = 0.d0
        s2alp     = 0.d0
        s2rho     = 0.d0
        s3        = 0.d0
        sswt      = 0.d0
        twothd    = 0.d0
        vmid      = 0.d0
        xect      = 0.d0
        xeh2      = 0.d0
        xes1      = 0.d0
        xes3      = 0.d0
        des1      = 0.55118d0
        bes1      = 1.05118d0
        res1      = 3.38583d0
        reh2      = 1.401d0
        bf        = 1.627d0
        b0        = 1.194d0
        gamh2a    = 2.323d0
        gamh2b    = 3.126d0
        deh2      = 4.74806d0
        des20     = 0.996d0
        s2alp     = 3.0d0
        s2rho     = 7.0d0
        dect      = 2.24409d0
        bect      = 0.57874d0
        rect      = 3.46457d0
        bect2     = 0.87143d0
        h2a1      = 2.38029900d0
        h2a2      = 1.47400334d0
        h2a3      = -0.226331771d0
        h2a4      = -2.17339562d0
        h2a5      = 9.54597998d0
        h2a6      = -23.3366749d0
        h2a7      = 24.7115748d0
        h2a8      = -12.8786102d0
        h2a9      = 2.75711758d0
        h2b1      = 1.35164094d0
        h2b2      = -0.254404564d0
        h2b3      = -10.6627055d0
        h2b4      = -14.0112406d0
        h2b5      = -7.80660398d0
        h2b6      = 5.47697753d0
        h2b7      = -6.95875664d0
        h2b8      = 3.64334771d0
        h2b9      = -0.492505596d0
        h2i1      = 160.484d0
        h2i2      = -20.397d0
        h2i3      = 3.913d0
        h2i4      = 0.777d0
        h2i5      = -7.583
        h2i6      = 1.501
        rnaie     = 5.139d0
        rnapol    = 159.267d0
        hpol      = 4.5d0
        hie       = 13.598d0
        hnai1     = -193.341d0
        hnai2     = 16.968d0
        hnai3     = -77.120d0
        hnai4     = -1.024d0
        hnai5     = 6.009d0
        hnai6     = -1.241d0
        hnau1     = 119.880d0
        hnau2     = -0.620d0
        hnau3     = 31.335d0
        hnau4     = 70.414d0
        hnau5     = -13.016d0
        hnau6     = 0.8681d0
        hnaa1     = 0.67d0
        hnaa2     = -92.229d0
        hnaa3     = 51.01019d0
        hnaa4     = -130.37826d0
        hnaa5     = 6.40285d0
        hnab1     = 0.34333d0
        hnab2     = -92.229d0
        hnab3     = -48.54867d0
        hnab4     = -17.30523d0
        hnab5     = 0.99821d0
        balp      = 2.07234d0
        brho      = 2.0d0
        rb0       = 5.1998697814229d7
        rb90      = 5.0d8
        f90       = 1.0d0
        f0        = 1.0d0
        rd        = 9.1538750508008d8
        dalp      = 1.43307d0
        drho      = 5.71654d0
        rac2      = 1.d0/dsqrt(2.d0)
        cevau     = 1.d0/27.2113961d0
        onethd    = 1.0d0/3.0d0
        twothd    = 2.0d0/3.0d0
        !
        do 10 i=1, nt
            ! ============================================================
            ! useful quantities calculated first, r1^2, R2^2, R3^2, atom to diatom
            ! center of mass distances, the H-Na-H bond angle and the Jacobian
            ! angle chi.
            r1s   = r(i,1)*r(i,1)
            r2s   = r(i,2)*r(i,2)
            r3s   = r(i,3)*r(i,3)
            rnag  = dsqrt(0.5d0*dabs(r1s+r3s-0.5d0*r2s))
            rnags = rnag*rnag
            rnagc = rnags*rnag
            rmid  = 0.5d0*(r(i,1)+r(i,3))
            rave  = 0.5d0*(r(i,1)-r(i,3))
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            cos1 = (r1s+r2s-r3s)/(2.0d0*r(i,1)*r(i,2))
            cos2 = (r3s+r2s-r1s)/(2.0d0*r(i,3)*r(i,2))
            cosg = 0.5d0*(cos1-cos2)
            frac = 0.958d0
            ! H to NaH center of masses
            rnag1  = dsqrt(dabs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))
            rnag3  = dsqrt(dabs(r2s+frac*(r1s-r2s-r3s)+r3s*frac*frac))
            rnag1c = rnag1*rnag1*rnag1
            rnag3c = rnag3*rnag3*rnag3
            ! the shortest of any two distances . . .
            rsm1  = r(i,3)+(r(i,2)-r(i,3)) *(0.5d0+0.5d0*dtanh(0.5d0*(r(i,3)-r(i,2))))
            rsm2  = r(i,3)+(r(i,1)-r(i,3)) *(0.5d0+0.5d0*dtanh(0.5d0*(r(i,3)-r(i,1))))
            rsm3  = r(i,1)+(r(i,2)-r(i,1)) *(0.5d0+0.5d0*dtanh(0.5d0*(r(i,1)-r(i,2))))
            rexp1 = dexp(-r(i,1))
            rexp2 = dexp(-r(i,2))
            rexp3 = dexp(-r(i,3))
            ! checks on angles . . .
            if( rnag .eq. 0.0d0 )then
                cstsq = 1.0d0
            else
                cstsq = (0.5d0*(r3s-r1s)/(r(i,2)*rnag))**2
            endif
            cosa   = (r1s+r3s-r2s)/(2*r(i,1)*r(i,3))
            coscos = cstsq*(0.5d0+0.5d0*cosa)*(1.0d0-dexp(-10.0d0*rnags))
            if( rnag1 .eq. 0.0d0 )then
                cstsq1 = 1.0d0
            else
                cstsq1 = 0.25d0*((r3s-r2s-r1s+2.0d0*r3s*frac)/(r(i,1)*rnag1)) **2
            endif
            if( rnag3 .eq. 0.0d0 )then
                cstsq3 = 1.0d0
            else
                cstsq3 = 0.25d0*((r1s-r2s-r3s+2.0d0*r1s*frac)/(r(i,3)*rnag3)) **2
            endif
            !==============================================================k
            !This section contains the diatomic curves
            ! The first NaH singlet
            xes1 = dexp(-bes1*(r(i,1)-res1))
            s1   = des1*onethd*xes1*(xes1+2.0d0)
            ! The H2 singlet
            ! new singlet, designed to minic lower adiabatic behavior
            reh2   = 1.401d0
            bf     = 1.627d0
            b0     = 1.194d0
            gamh2a = 2.323d0
            gamh2b = 3.126d0
            beh2   = bf*(b0-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2) /(bf-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2)
            sswt   = 0.5d0*(1.0d0-dtanh((rmid-s2rho)/s2alp))*cosg**2
            deh2m  = deh2 + deh2*(des20 - 1.0d0)*sswt
            xeh2   = dexp(-beh2*(r(i,2)-reh2))
            s2     = deh2m*xeh2*(xeh2-2.0d0)
            ! The second NaH singlet
            xes3 = dexp(-bes1*(r(i,3)-res1))
            s3   = des1*onethd*xes3*(xes3+2.0d0)
            ! ===============================================================
            ! this long section contains data pertinent to long range forces
            !  H2 ionization energy:
            h2ie = (h2i1+h2i2*r(i,2))*dexp(-h2i3*r(i,2))+ r2s*(h2i4+h2i5*r(i,2))*dexp(-h2i6*r(i,2))-s2a+hie
            !  H2 polarizability:
            !   al1 is the parallel polarizability for h2
            !   al2 is the perpendicaular polarizability for h2
            !   alt is the anisotropic? polarizability
            !   alb is the mean polarizability
            al1=dexp(h2a1*(h2a2-r(i,2)))*(h2a3+h2a4*r(i,2)+h2a5*r2s +h2a6*r(i,2)*r2s+h2a7*r2s*r2s+h2a8*r2s*r2s* r(i, &
                & 2) + h2a9* r2s*r2s*r2s) +2.0d0*hpol
            al2=dexp(h2b1*(h2b2-r(i,2)))*(h2b3+h2b4*r(i,2)+h2b5*r2s +h2b6*r(i,2)*r2s+h2b7*r2s*r2s+h2b8*r2s*r2s* r(i, &
                & 2) + h2b9* r2s*r2s*r2s) +2.0d0*hpol
            alt  = twothd*(al1-al2)
            alb  = onethd*al1+twothd*al2
            alhh = alb+alt*(1.5d0*cstsq-0.5d0)
            !  The NaH + H dispersion and NaH + H induction terms must go on the low
            !  so that the adiabat is fitted correctly in the region of the coupling
            ! NaH properties
            rnahie1   = (hnai1+hnai2*r(i,1)) /(hnai3+hnai4*r(i,1)+hnai5*r1s+hnai6*r(i,1)*r1s)+rnaie
            rnahie3   = (hnai1+hnai2*r(i,3)) /(hnai3+hnai4*r(i,3)+hnai5*r3s+hnai6*r(i,3)*r3s)+rnaie
            rnahmu1   = (hnau1*r(i,1)+hnau2*r1s) /(hnau3+hnau4*r(i,1)+hnau5*r1s+hnau6*r(i,1)*r1s)
            rnahmu3   = (hnau1*r(i,3)+hnau2*r3s) /(hnau3+hnau4*r(i,3)+hnau5*r3s+hnau6*r(i,3)*r3s)
            rnahpar1  = rnapol+hpol+dexp(-hnaa1*r(i,1)) *(hnaa2+hnaa3*r(i,1)+hnaa4*r1s+hnaa5*r1s*r(i,1))
            rnahpar3  = rnapol+hpol+dexp(-hnaa1*r(i,3)) *(hnaa2+hnaa3*r(i,3)+hnaa4*r3s+hnaa5*r3s*r(i,3))
            rnahperp1 = rnapol+hpol+dexp(-hnab1*r(i,1))* (hnab2+hnab3*r(i,1)+hnab4*r1s+hnab5*r1s*r(i,1))
            rnahperp3 = rnapol+hpol+dexp(-hnab1*r(i,3))* (hnab2+hnab3*r(i,3)+hnab4*r3s+hnab5*r3s*r(i,3))
            rnahpol1  = rnahperp1+(rnahpar1-rnahperp1)*cstsq1
            rnahpol1b = onethd*rnahpar1+twothd*rnahperp1
            rnahpol3  = rnahperp3+(rnahpar3-rnahperp3)*cstsq3
            rnahpol3b = onethd*rnahpar3+twothd*rnahperp3
            ! induction forces
            ! dipole induced dipole
            div1    = 1.0d0/(rnag1c*rnag1c*rnag1c+rd)
            rmuimu1 = -rnahmu1**2*hpol*(1.5d0*cstsq1+0.5d0)*rnag1c*div1
            div3    = 1.0d0/(rnag3c*rnag3c*rnag3c+rd)
            rmuimu3 = -rnahmu3**2*hpol*(1.5d0*cstsq3+0.5d0)*rnag3c*div3
            ! dispersion forces
            c690   = -1.5d0*rnaie*h2ie*rnapol*al2 /((rnaie+h2ie)*(rnagc*rnagc*rnagc+rb90))
            c690   = c690*f90*rnagc
            c60    = -1.5d0*rnaie*h2ie*rnapol*al1/ ((rnaie+h2ie)*(rnagc*rnagc*rnagc+rb0))
            c60    = c60*f0*rnagc
            c6     = c690+(c60-c690)*cstsq*(1.0d0-dexp(-39.0625d0*rnag*rnagc))
            c61    = -1.5d0*rnahie1*hie*rnag1c/(rnahie1+hie)*rnahpol1*hpol*div1
            dexp1  = dexp(-39.0625d0*rnag1c*rnag1)
            c61    = c61*(1.0d0+(rnahpol1b/rnahpol1-1.0d0)*dexp1)
            c63    = -1.5d0*rnahie3*hie*rnag3c/(rnahie3+hie)*rnahpol3*hpol*div3
            dexp3  = dexp(-39.0625d0*rnag3c*rnag3)
            c63    = c63*(1.0d0+(rnahpol3b/rnahpol3-1.0d0)*dexp3)
            cind1  = rmuimu1*(1.0d0+(1.0d0/(1.5d0*cstsq1+0.5d0)-1.0d0)*dexp1)
            cind3  = rmuimu3*(1.0d0+(1.0d0/(1.5d0*cstsq3+0.5d0)-1.0d0)*dexp3)
            c6     = c6*cevau*0.5d0*(1.0d0-dtanh(balp*(r(i,2)-brho)))
            c6swt1 = 0.5d0*(1.0d0-dtanh(dalp*(r(i,1)-drho)))
            c6swt3 = 0.5d0*(1.0d0-dtanh(dalp*(r(i,3)-drho)))
            c61    = c61*c6swt1*cevau
            c63    = c63*c6swt3*cevau
            cind1  = cind1*c6swt1
            cind3  = cind3*c6swt3
            elec   = c61+c63+cind1+cind3+c6
            ! =========================================================
            ! here is where it all gets put together . . .
            e(i) = s1 + s2 + s3 + deh2
            elec = elec/cevau
            e(i) = e(i)+elec
            ! here we will try a post-LEPS correction term
            xect = dexp(-bect*(rmid-rect))
            vmid = dect*xect*(xect-2.0d0) *dexp(-bect2*rave**2)*(0.5d0*(1.0d0-cosa))**2
            e(i) = e(i)+vmid
            e(i) = e(i)*cevau
        10 continue
    end subroutine dispatch
end module chempotpy_nah2_leps_7l_2000_dpem_pot11_legacy_state_f8f175af

!> @brief Compatibility interface for pot11.
subroutine pot11( r, e, nt )
    use chempotpy_nah2_leps_7l_2000_dpem_pot11_legacy_state_f8f175af, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    integer          :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot11

! U12 coupling surface . . .

!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_nah2_leps_7l_2000_dpem_pot12_legacy_state_f8f175af
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        implicit none
        optional :: r, e, nt
        double precision :: ap, alp, bet, rnagsq, rnag, cosg, rd0, r20
        double precision :: r1s, r2s, r3s, cos1, cos2, theta
        double precision :: gc, gc2
        double precision :: cevau
        integer          :: nt, i, j, k
        logical          :: leven
        double precision :: r(:,:), e(:)
        alp    = 0.d0
        ap     = 0.d0
        bet    = 0.d0
        cevau  = 0.d0
        cos1   = 0.d0
        cos2   = 0.d0
        cosg   = 0.d0
        gc     = 0.d0
        gc2    = 0.d0
        r1s    = 0.d0
        r20    = 0.d0
        r2s    = 0.d0
        r3s    = 0.d0
        rd0    = 0.d0
        rnag   = 0.d0
        rnagsq = 0.d0
        theta  = 0.d0
        ap     = 0.835d0
        alp    = 4.0d0
        rd0    = 4.0d0
        bet    = 0.5d0
        r20    = -1.3d0
        theta  = 15.6d0*dacos(-1.0d0)/180.0d0
        cevau  = 1.d0/27.2113961d0
        ! set leven true for an even U12 surface, and false for an odd
        ! surface
        leven = .true.
        do i=1, nt
            r1s = r(i,1)**2
            r2s = r(i,2)**2
            r3s = r(i,3)**2
            ! here we do the function evaluations . . .
            rnagsq = 0.5d0*r1s+0.5d0*r3s-0.25d0*r2s
            if( rnagsq .gt. 0.0d0 )then
                rnag = dsqrt(rnagsq)
            else
                rnag = 0.0d0
            endif
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            cos1 = (r1s+r2s-r3s)/(2.0d0*r(i,1)*r(i,2))
            cos2 = (r3s+r2s-r1s)/(2.0d0*r(i,3)*r(i,2))
            cosg = 0.5d0*(cos1-cos2)
            gc   = dcos(theta)*rnag + dsin(theta)*r(i,2)
            gc2  = dsin(theta)*rnag - dcos(theta)*r(i,2)
            e(i) = ap*cosg *dexp(-(1.0d0/alp)**4*(gc - rd0)**4) *dexp(-(1.0d0/bet)**2*(gc2 - r20)**2)
            e(i) = e(i)*cevau
            ! set leven equal to true for an even surface, false for an odd one
            if( leven )then
                e(i) = dsqrt(e(i)**4/(e(i)**2+1.0d-8))
            endif
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_7l_2000_dpem_pot12_legacy_state_f8f175af

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, nt )
    use chempotpy_nah2_leps_7l_2000_dpem_pot12_legacy_state_f8f175af, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    integer          :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot12

! the upper diabatic surface  U22

!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_nah2_leps_7l_2000_dpem_pot22_legacy_state_f8f175af
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        implicit none
        optional :: r, e, nt
        double precision :: al1, al2, alb, alhh, alt, b0, balp, beh2, bet1, bet10, bet10b, bet190, bet1c, &
            & bet20, bet290, beta, bf, brho, c2a, c2b, c2c, c6, c61, c63, c6swt1, c6swt3, cevau, cind1, &
            & cind3, cos1, cos2, cosa, coscos, cosg, coul1, coul2, coul3, cplg2, cstsq, cstsq1, cstsq3, &
            & dalp, deh2, deh2m, denah, des290, det1, det190, det1c, det20, det290, dexp1, dexp2, dexp3, &
            & div1, div3, drho, e, elec, exch1, exch2, exch3, f10, frac, g10, g30, gamh2a, gamh2b, &
            & gamnah, h2a1, h2a2, h2a3, h2a4, h2a5, h2a6, h2a7, h2a8, h2a9, h2b1, h2b2, h2b3, h2b4, h2b5, &
            & h2b6, h2b7, h2b8, h2b9, h2hc1, h2hc2, h2hc3, h2hc4, h2hc5, h2hex, h2i1, h2i2, h2i3, h2i4, &
            & h2i5, h2i6, h2ie, h2qc1, h2qc2, h2qc3, h2qc4, h2qc5, h2quad, hie, hnaa1, hnaa2, hnaa3, &
            & hnaa4, hnaa5, hnab1, hnab2, hnab3, hnab4, hnab5, hnai1, hnai2, hnai3, hnai4, hnai5, hnai6, &
            & hnau1, hnau2, hnau3, hnau4, hnau5, hnau6, hpol, onethd, pi, qh, qhint, qq, qqint, r, r1s, &
            & r2s, r3s, rave, rb, rd, reh2, renah, ret1, ret10, ret190, ret1c, ret20, ret290, rmid, &
            & rmuimu1, rmuimu3, rnag, rnag1, rnag1c, rnag3, rnag3c, rnagc, rnags, rnahie1, rnahie3, &
            & rnahmu1, rnahmu3, rnahpar1, rnahpar3, rnahperp1, rnahperp3, rnahpol1, rnahpol1b, rnahpol3, &
            & rnahpol3b, rnaie, rnapol1, rnaquad, rqh, rqq, s1, s2, s20alp, s20b, s20rho, s2a, s2alp, &
            & s2rho, s3, sswt, t1, t10, t190, t1alp, t1c, t1m, t1rho, t2, t20, t290, t2alp, t2rho, t3, &
            & t30, t390, t3c, t3m, twothd, una2p, w, xeh2, xenah, xet190, xet1c, xet20, xet290, xet3, &
            & xet390, xet3c
        integer :: i, nt
        dimension r(:,:), e(:)
        ! NaH triplets
        ! H2 triplets
        ! H2 singlets
        ! misc constants
        ! data blocks for the long range forces
        al1       = 0.d0
        al2       = 0.d0
        alb       = 0.d0
        alhh      = 0.d0
        alt       = 0.d0
        b0        = 0.d0
        balp      = 0.d0
        beh2      = 0.d0
        bet1      = 0.d0
        bet10     = 0.d0
        bet10b    = 0.d0
        bet190    = 0.d0
        bet1c     = 0.d0
        bet20     = 0.d0
        bet290    = 0.d0
        beta      = 0.d0
        bf        = 0.d0
        brho      = 0.d0
        c2a       = 0.d0
        c2b       = 0.d0
        c2c       = 0.d0
        c6        = 0.d0
        c61       = 0.d0
        c63       = 0.d0
        c6swt1    = 0.d0
        c6swt3    = 0.d0
        cevau     = 0.d0
        cind1     = 0.d0
        cind3     = 0.d0
        cos1      = 0.d0
        cos2      = 0.d0
        cosa      = 0.d0
        coscos    = 0.d0
        cosg      = 0.d0
        coul1     = 0.d0
        coul2     = 0.d0
        coul3     = 0.d0
        cplg2     = 0.d0
        cstsq     = 0.d0
        cstsq1    = 0.d0
        cstsq3    = 0.d0
        dalp      = 0.d0
        deh2      = 0.d0
        deh2m     = 0.d0
        denah     = 0.d0
        des290    = 0.d0
        det1      = 0.d0
        det190    = 0.d0
        det1c     = 0.d0
        det20     = 0.d0
        det290    = 0.d0
        dexp1     = 0.d0
        dexp2     = 0.d0
        dexp3     = 0.d0
        div1      = 0.d0
        div3      = 0.d0
        drho      = 0.d0
        elec      = 0.d0
        exch1     = 0.d0
        exch2     = 0.d0
        exch3     = 0.d0
        f10       = 0.d0
        frac      = 0.d0
        g10       = 0.d0
        g30       = 0.d0
        gamh2a    = 0.d0
        gamh2b    = 0.d0
        gamnah    = 0.d0
        h2a1      = 0.d0
        h2a2      = 0.d0
        h2a3      = 0.d0
        h2a4      = 0.d0
        h2a5      = 0.d0
        h2a6      = 0.d0
        h2a7      = 0.d0
        h2a8      = 0.d0
        h2a9      = 0.d0
        h2b1      = 0.d0
        h2b2      = 0.d0
        h2b3      = 0.d0
        h2b4      = 0.d0
        h2b5      = 0.d0
        h2b6      = 0.d0
        h2b7      = 0.d0
        h2b8      = 0.d0
        h2b9      = 0.d0
        h2hc1     = 0.d0
        h2hc2     = 0.d0
        h2hc3     = 0.d0
        h2hc4     = 0.d0
        h2hc5     = 0.d0
        h2hex     = 0.d0
        h2i1      = 0.d0
        h2i2      = 0.d0
        h2i3      = 0.d0
        h2i4      = 0.d0
        h2i5      = 0.d0
        h2i6      = 0.d0
        h2ie      = 0.d0
        h2qc1     = 0.d0
        h2qc2     = 0.d0
        h2qc3     = 0.d0
        h2qc4     = 0.d0
        h2qc5     = 0.d0
        h2quad    = 0.d0
        hie       = 0.d0
        hnaa1     = 0.d0
        hnaa2     = 0.d0
        hnaa3     = 0.d0
        hnaa4     = 0.d0
        hnaa5     = 0.d0
        hnab1     = 0.d0
        hnab2     = 0.d0
        hnab3     = 0.d0
        hnab4     = 0.d0
        hnab5     = 0.d0
        hnai1     = 0.d0
        hnai2     = 0.d0
        hnai3     = 0.d0
        hnai4     = 0.d0
        hnai5     = 0.d0
        hnai6     = 0.d0
        hnau1     = 0.d0
        hnau2     = 0.d0
        hnau3     = 0.d0
        hnau4     = 0.d0
        hnau5     = 0.d0
        hnau6     = 0.d0
        hpol      = 0.d0
        onethd    = 0.d0
        pi        = 0.d0
        qh        = 0.d0
        qhint     = 0.d0
        qq        = 0.d0
        qqint     = 0.d0
        r1s       = 0.d0
        r2s       = 0.d0
        r3s       = 0.d0
        rave      = 0.d0
        rb        = 0.d0
        rd        = 0.d0
        reh2      = 0.d0
        renah     = 0.d0
        ret1      = 0.d0
        ret10     = 0.d0
        ret190    = 0.d0
        ret1c     = 0.d0
        ret20     = 0.d0
        ret290    = 0.d0
        rmid      = 0.d0
        rmuimu1   = 0.d0
        rmuimu3   = 0.d0
        rnag      = 0.d0
        rnag1     = 0.d0
        rnag1c    = 0.d0
        rnag3     = 0.d0
        rnag3c    = 0.d0
        rnagc     = 0.d0
        rnags     = 0.d0
        rnahie1   = 0.d0
        rnahie3   = 0.d0
        rnahmu1   = 0.d0
        rnahmu3   = 0.d0
        rnahpar1  = 0.d0
        rnahpar3  = 0.d0
        rnahperp1 = 0.d0
        rnahperp3 = 0.d0
        rnahpol1  = 0.d0
        rnahpol1b = 0.d0
        rnahpol3  = 0.d0
        rnahpol3b = 0.d0
        rnaie     = 0.d0
        rnapol1   = 0.d0
        rnaquad   = 0.d0
        rqh       = 0.d0
        rqq       = 0.d0
        s1        = 0.d0
        s2        = 0.d0
        s20alp    = 0.d0
        s20b      = 0.d0
        s20rho    = 0.d0
        s2a       = 0.d0
        s2alp     = 0.d0
        s2rho     = 0.d0
        s3        = 0.d0
        sswt      = 0.d0
        t1        = 0.d0
        t10       = 0.d0
        t190      = 0.d0
        t1alp     = 0.d0
        t1c       = 0.d0
        t1m       = 0.d0
        t1rho     = 0.d0
        t2        = 0.d0
        t20       = 0.d0
        t290      = 0.d0
        t2alp     = 0.d0
        t2rho     = 0.d0
        t3        = 0.d0
        t30       = 0.d0
        t390      = 0.d0
        t3c       = 0.d0
        t3m       = 0.d0
        twothd    = 0.d0
        una2p     = 0.d0
        w         = 0.d0
        xeh2      = 0.d0
        xenah     = 0.d0
        xet190    = 0.d0
        xet1c     = 0.d0
        xet20     = 0.d0
        xet290    = 0.d0
        xet3      = 0.d0
        xet390    = 0.d0
        xet3c     = 0.d0
        ret10     = 3.22835d0
        bet10     = 0.57874d0
        bet10b    = 2.44094d0
        f10       = 9.44882d0
        det190    = 4.32441d0
        bet190    = 2.43701d0
        ret190    = 1.95480d0
        t1alp     = 0.30778d0
        t1rho     = 3.21850d0
        det1c     = 4.92126d0
        bet1c     = 0.65748d0
        ret1c     = 2.71654d0
        det20     = 3.29921d0
        bet20     = 2.64567d0
        ret20     = 2.44882d0
        det290    = 3.97638d0
        bet290    = 3.73937d0
        ret290    = 0.96024d0
        t2alp     = 1.0d0
        t2rho     = 9.0d0
        des290    = 0.91433d0
        s2alp     = 0.66778d0
        s2rho     = 5.30472d0
        deh2      = 4.74806d0
        s20alp    = 0.5d0
        s20rho    = 6.0d0
        c2a       = 1.38661d0
        c2b       = 0.27362d0
        c2c       = 0.03764d0
        una2p     = 2.1037d0
        h2a1      = 2.38029900d0
        h2a2      = 1.47400334d0
        h2a3      = -0.226331771d0
        h2a4      = -2.17339562d0
        h2a5      = 9.54597998d0
        h2a6      = -23.3366749d0
        h2a7      = 24.7115748d0
        h2a8      = -12.8786102d0
        h2a9      = 2.75711758d0
        h2b1      = 1.35164094d0
        h2b2      = -0.254404564d0
        h2b3      = -10.6627055d0
        h2b4      = -14.0112406d0
        h2b5      = -7.80660398d0
        h2b6      = 5.47697753d0
        h2b7      = -6.95875664d0
        h2b8      = 3.64334771d0
        h2b9      = -0.492505596d0
        h2i1      = 160.484d0
        h2i2      = -20.397d0
        h2i3      = 3.913d0
        h2i4      = 0.777d0
        h2i5      = -7.583
        h2i6      = 1.501
        rnaie     = 3.0353
        rnapol1   = 290.4d0
        hpol      = 4.5d0
        hie       = 13.598d0
        rnaquad   = -16.8992d0
        hnai1     = -193.341d0
        hnai2     = 16.968d0
        hnai3     = -77.120d0
        hnai4     = -1.024d0
        hnai5     = 6.009d0
        hnai6     = -1.241d0
        hnau1     = 119.880d0
        hnau2     = -0.620d0
        hnau3     = 31.335d0
        hnau4     = 70.414d0
        hnau5     = -13.016d0
        hnau6     = 0.8681d0
        hnaa1     = 0.67d0
        hnaa2     = -92.229d0
        hnaa3     = 51.01019d0
        hnaa4     = -130.37826d0
        hnaa5     = 6.40285d0
        hnab1     = 0.34333d0
        hnab2     = -92.229d0
        hnab3     = -48.54867d0
        hnab4     = -17.30523d0
        hnab5     = 0.99821d0
        rqq       = 9.766021100d4
        rqh       = 2.355466062927d6
        rb        = 9.3119383371714d5
        balp      = 1.93701d0
        brho      = 1.4d0
        rd        = 9.1351724748364d8
        dalp      = 1.43307d0
        drho      = 5.71654d0
        h2hc1     = 0.677d0
        h2hc2     = 0.378d0
        h2hc3     = 0.281d0
        h2hc4     = -1.022d0
        h2hc5     = 0.665d0
        h2qc1     = 1.532d0
        h2qc2     = -.04635d0
        h2qc3     = 1.826d0
        h2qc4     = -2.508d0
        h2qc5     = 1.856d0
        cevau     = 1.d0/27.2113961d0
        onethd    = 1.0d0/3.0d0
        twothd    = 2.0d0/3.0d0
        do 10 i=1, nt
            ! =========================================================
            ! this section calcluates useful quantities - the A to BC center
            ! of mass distance, the angles between A-BC and BC, the shorest
            ! of any of the three diatomic distances, and the bond angle
            ! H-Na-H
            rmid  = 0.5d0*(r(i,1)+r(i,3))
            rave  = 0.5d0*(r(i,1)-r(i,3))
            r1s   = r(i,1)*r(i,1)
            r2s   = r(i,2)*r(i,2)
            r3s   = r(i,3)*r(i,3)
            rnag  = dsqrt(0.5d0*dabs(r1s+r3s-0.5d0*r2s))
            rnags = rnag*rnag
            rnagc = rnags*rnag
            if( rnag .eq. 0.0d0 )then
                cstsq = 0.0d0
            else
                cstsq = (0.5d0*(r3s-r1s)/(r(i,2)*rnag))**2
            endif
            cosa   = (r1s+r3s-r2s)/(2*r(i,1)*r(i,3))
            coscos = cstsq*(0.5d0+0.5d0*cosa)*(1.0d0-dexp(-10.0d0*rnags))
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            cos1   = (r1s+r2s-r3s)/(2.0d0*r(i,1)*r(i,2))
            cos2   = (r3s+r2s-r1s)/(2.0d0*r(i,3)*r(i,2))
            cosg   = 0.5d0*(cos1-cos2)
            frac   = 0.958d0
            rnag1  = dsqrt(dabs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))
            rnag3  = dsqrt(dabs(r2s+frac*(r1s-r2s-r3s)+r3s*frac*frac))
            rnag1c = rnag1*rnag1*rnag1
            rnag3c = rnag3*rnag3*rnag3
            if( rnag1 .eq. 0.0d0 )then
                cstsq1 = 1.0d0
            else
                cstsq1 = 0.25d0*((r3s-r2s-r1s+2.0d0*r3s*frac)/(r(i,1)*rnag1)) **2
            endif
            if( rnag3 .eq. 0.0d0 )then
                cstsq3 = 1.0d0
            else
                cstsq3 = 0.25d0*((r1s-r2s-r3s+2.0d0*r1s*frac)/(r(i,3)*rnag3)) **2
            endif
            ! =========================================================
            ! The first NaH singlet
            ! new singlet, designed to minic lower adiabatic behavior
            denah  = 1.97134878d0
            renah  = 3.566044d0
            bf     = 0.864d0
            b0     = 0.594d0
            gamnah = 7.376d0
            beta   = bf*(b0+(r(i,1)/gamnah)**8)/(bf+(r(i,1)/gamnah)**8)
            xenah  = dexp(-beta*(r(i,1)-renah))
            s1     = denah*xenah*(xenah-2.0d0)
            ! The first NaH triplet
            g10    = r(i,1)-ret10
            t10    = f10*dexp(-bet10*g10)+dexp(-bet10b*g10)
            t10    = t10/(1.0d0+f10)
            xet190 = dexp(-bet190*(r(i,1)-ret190))
            t190   = det190*onethd*xet190*(xet190 + 2.0d0)
            xet1c  = dexp(-bet1c*(r(i,1)-ret1c))
            t1c    = det1c*onethd*xet1c*(xet1c + 2.0d0)
            t1m    = t190 +(t1c-t190)* 0.5d0*(1.0d0+dtanh(t1alp*(r(i,2)-t1rho)))
            t1     = t1m + (t10-t1m)*cosg**2
            ! The H2 triplet
            xet20  = dexp(-bet20*(r(i,2)-ret20))
            t20    = det20*onethd*xet20*(xet20 + 2.0d0)
            xet290 = dexp(-bet290*(r(i,2)-ret290))
            t290   = det290*onethd*xet290*(xet290 + 2.0d0)
            t2     = t20 + (t290-t20)*(1.0d0-cosg**2)* 0.5d0*(1.0d0-dtanh(t2alp*(rmid-t2rho)))
            ! The H2 singlet
            ! new singlet, designed to minic lower adiabatic behavior
            reh2   = 1.401d0
            bf     = 1.627d0
            b0     = 1.194d0
            gamh2a = 2.323d0
            gamh2b = 3.126d0
            beh2   = bf*(b0-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2) /(bf-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2)
            sswt   = 0.5d0*(1.0d0-dtanh(s2alp*(rmid-s2rho)))
            deh2m  = deh2 + deh2*(des290 - 1.0d0)*sswt
            xeh2   = dexp(-beh2*(r(i,2)-reh2))
            s2a    = deh2m*xeh2*(xeh2-2.0d0) + una2p
            s2     = 0.5d0*(s2a+t2-dtanh((s2a-t2)/0.1d0)*(s2a-t2))
            xeh2   = dexp(-1.34d0*beh2*(r(i,2)-reh2))
            s20b   = (deh2-una2p)*xeh2*(xeh2-2.0d0)
            s2     = s2 + (s20b - s2)*cosg**2* 0.5d0*(1.0d0-dtanh(s20alp*(rmid-s20rho)))
            ! The second NaH singlet
            ! new singlet, designed to minic lower adiabatic behavior
            denah  = 1.97134878d0
            renah  = 3.566044d0
            bf     = 0.864d0
            b0     = 0.594d0
            gamnah = 7.376d0
            beta   = bf*(b0+(r(i,3)/gamnah)**8)/(bf+(r(i,3)/gamnah)**8)
            xenah  = dexp(-beta*(r(i,3)-renah))
            s3     = denah*xenah*(xenah-2.0d0)
            ! The second NaH triplet
            xet3   = dexp(-bet1*(r(i,3)-ret1))
            t30    = det1*onethd*xet3*(xet3 + 2.0d0)
            g30    = r(i,3)-ret10
            t30    = f10*dexp(-bet10*g30)+dexp(-bet10b*g30)
            t30    = t30/(1.0d0+f10)
            xet390 = dexp(-bet190*(r(i,3)-ret190))
            t390   = det190*onethd*xet390*(xet390 + 2.0d0)
            xet3c  = dexp(-bet1c*(r(i,3)-ret1c))
            t3c    = det1c*onethd*xet3c*(xet3c + 2.0d0)
            t3m    = t390 + (t3c - t390)* 0.5d0*(1.0d0+dtanh(t1alp*(r(i,2)-t1rho)))
            t3     = t3m + (t30 - t3m)*cosg**2
            ! =========================================================
            ! this long section contains data pertinent to long range forces
            !  H2 ionization energy:
            h2ie = (h2i1+h2i2*r(i,2))*dexp(-h2i3*r(i,2))+ r2s*(h2i4+h2i5*r(i,2))*dexp(-h2i6*r(i,2))-s2a+hie
            !  H2 polarizability:
            !   al1 is the parallel polarizability for h2
            !   al2 is the perpendicaular polarizability for h2
            !   alt is the anisotropic? polarizability
            !   alb is the mean polarizability
            al1=dexp(h2a1*(h2a2-r(i,2)))*(h2a3+h2a4*r(i,2)+h2a5*r2s +h2a6*r(i,2)*r2s+h2a7*r2s*r2s+h2a8*r2s*r2s* r(i, &
                & 2) + h2a9* r2s*r2s*r2s) +2.0d0*hpol
            al2=dexp(h2b1*(h2b2-r(i,2)))*(h2b3+h2b4*r(i,2)+h2b5*r2s +h2b6*r(i,2)**3+h2b7*r2s*r2s+h2b8*r2s*r2s*r(i, &
                & 2) + h2b9* r2s*r2s*r2s) +2.0d0*hpol
            alt  = twothd*(al1-al2)
            alb  = onethd*al1+twothd*al2
            alhh = alb+alt*(1.5d0*cstsq-0.5d0)
            ! here follow properties necesary for the long range forces in
            ! the Na + H2 channel:
            ! H2 properties:
            h2hex  = dexp(-h2hc1*r(i,2)) *(h2hc2*r(i,2)+h2hc3*r2s+h2hc4*r(i,2)*r2s+h2hc5*r2s*r2s)
            h2quad = dexp(-h2qc1*r(i,2)) *(h2qc2*r(i,2)+h2qc3*r2s+h2qc4*r(i,2)*r2s+h2qc5*r2s*r2s)
            ! NaH properties
            rnahie1   = (hnai1+hnai2*r(i,1)) /(hnai3+hnai4*r(i,1)+hnai5*r1s+hnai6*r(i,1)*r1s)+rnaie+una2p
            rnahie3   = (hnai1+hnai2*r(i,3)) /(hnai3+hnai4*r(i,3)+hnai5*r3s+hnai6*r(i,3)*r3s)+rnaie+una2p
            rnahmu1   = (hnau1*r(i,1)+hnau2*r1s) /(hnau3+hnau4*r(i,1)+hnau5*r1s+hnau6*r(i,1)*r1s)
            rnahmu3   = (hnau1*r(i,3)+hnau2*r3s) /(hnau3+hnau4*r(i,3)+hnau5*r3s+hnau6*r(i,3)*r3s)
            rnahpar1  = 159.267d0+hpol+dexp(-hnaa1*r(i,1)) *(hnaa2+hnaa3*r(i,1)+hnaa4*r1s+hnaa5*r1s*r(i,1))
            rnahpar3  = 159.267d0+hpol+dexp(-hnaa1*r(i,3)) *(hnaa2+hnaa3*r(i,3)+hnaa4*r3s+hnaa5*r3s*r(i,3))
            rnahperp1 = 159.267d0+hpol+dexp(-hnab1*r(i,1))* (hnab2+hnab3*r(i,1)+hnab4*r1s+hnab5*r1s*r(i,1))
            rnahperp3 = 159.267d0+hpol+dexp(-hnab1*r(i,3))* (hnab2+hnab3*r(i,3)+hnab4*r3s+hnab5*r3s*r(i,3))
            rnahpol1  = rnahperp1+(rnahpar1-rnahperp1)*cstsq1
            rnahpol1b = onethd*rnahpar1+twothd*rnahperp1
            rnahpol3  = rnahperp3+(rnahpar3-rnahperp3)*cstsq3
            rnahpol3b = onethd*rnahpar3+twothd*rnahperp3
            ! electrostatic forces
            ! quadrupole quadrupole
            qqint = 3.0d0-7.0d0*cstsq
            qq    = 0.75d0*rnaquad*h2quad/(rnagc*rnag*rnag+rqq)*qqint
            dexp2 = dexp(-39.0625d0*rnagc*rnag)
            qq    = qq*(1.0d0+(-4.0d0*pi/(3.0d0*qqint)-1.0d0)*dexp2)
            ! quadrupole hexadecapole
            qhint = -45.0d0-147.0d0*cstsq*cstsq+162.0d0*cstsq
            qh    = 0.3125d0*rnaquad*h2hex/(rnagc*rnagc*rnag+rqh)*(qhint)
            qh    = qh*(1.0d0+(-40.8d0*pi/qhint-1.0d0)*dexp2)
            ! induction forces
            ! dipole induced dipole
            div1    = 1.0d0/(rnag1c*rnag1c*rnag1c+rd)
            rmuimu1 = -rnahmu1**2*hpol*(1.5d0*cstsq1+0.5d0)*rnag1c*div1
            div3    = 1.0d0/(rnag3c*rnag3c*rnag3c+rd)
            rmuimu3 = -rnahmu3**2*hpol*(1.5d0*cstsq3+0.5d0)*rnag3c*div3
            ! dispersion forces
            c6     = -1.5d0*rnaie*h2ie/(rnaie+h2ie)*rnapol1*alhh/(rnagc*rnagc+rb)
            c6     = c6*(1.0d0+(alb/alhh-1.0d0)*dexp2)
            c61    = -1.5d0*rnahie1*hie/(rnahie1+hie)*rnahpol1*hpol*rnag1c*div1
            dexp1  = dexp(-39.0625d0*rnag1c*rnag1)
            c61    = c61*(1.0d0+(rnahpol1b/rnahpol1-1.0d0)*dexp1)
            c63    = -1.5d0*rnahie3*hie/(rnahie3+hie)*rnahpol3*hpol*rnag3c*div3
            dexp3  = dexp(-39.0625d0*rnag3c*rnag3)
            c63    = c63*(1.0d0+(rnahpol3b/rnahpol3-1.0d0)*dexp3)
            cind1  = rmuimu1*(1.0d0+(1.0d0/(1.5d0*cstsq1+0.5d0)-1.0d0)*dexp1)
            cind3  = rmuimu3*(1.0d0+(1.0d0/(1.5d0*cstsq3+0.5d0)-1.0d0)*dexp3)
            elec   = (qq+qh+c6*cevau)*0.5d0 *(1.0d0-dtanh(balp*(r(i,2)-brho)))
            c6swt1 = 0.5d0*(1.0d0-dtanh(dalp*(r(i,1)-drho)))
            c6swt3 = 0.5d0*(1.0d0-dtanh(dalp*(r(i,3)-drho)))
            c61    = c61*c6swt1*cevau
            c63    = c63*c6swt3*cevau
            cind1  = cind1*c6swt1
            cind3  = cind3*c6swt3
            elec   = elec+c61+c63+cind1+cind3
            ! ===============================================================
            coul1 = 0.5d0*(s1+t1)
            coul2 = 0.5d0*(s2+t2)
            coul3 = 0.5d0*(s3+t3)
            exch1 = 0.5d0*(s1-t1)
            exch2 = 0.5d0*(s2-t2)
            exch3 = 0.5d0*(s3-t3)
            w     = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
            cplg2 = c2a*dexp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
            e(i)  = coul1+coul2+coul3+deh2-dsqrt(w+cplg2*cplg2)/dsqrt(2.0d0)
            elec  = elec/cevau
            e(i)  = (e(i)+elec)*cevau
        10 continue
    end subroutine dispatch
end module chempotpy_nah2_leps_7l_2000_dpem_pot22_legacy_state_f8f175af

!> @brief Compatibility interface for pot22.
subroutine pot22( r, e, nt )
    use chempotpy_nah2_leps_7l_2000_dpem_pot22_legacy_state_f8f175af, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    integer          :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot22
