!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_nah2_leps_7s_2000_dpem_pes
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
end module chempotpy_nah2_leps_7s_2000_dpem_pes

!> @brief Preserve the historical Python-facing dpem entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_nah2_leps_7s_2000_dpem_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 3
    integer, parameter            :: nstates = 2
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: u(nstates,nstates)
    double precision, intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem

!> @brief Evaluate EvdRdX for the NaH2_LEPS_7S_2000_DPEM potential.
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

!> @brief Evaluate diagonalize for the NaH2_LEPS_7S_2000_DPEM potential.
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
!   Common name:                   NaH2 surface set 7s
!   Interface:                     3-1V
!   Number of electronic surfaces: 2
!   Number of derivatives:         0

!   Notes:               This version does not contain long range forces
!                        and has no derivatives.  It is based on the sur
!                        set 6 fit.

!   References:         Unpublished, M. D. Hack and D. G. Truhlar (2000)

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
!   and the potential coupling.
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
!                            R(2) = R(H-H)
!                            R(3) = R(Na-H)

!
! This is NaH2 potential matrix 7s, finalized on 5/17/99
! This version does not contain any long range forces.

!
!  all surfaces.  The potential is called by calling POT with the parame
!  r, e, nt, nsurf.  r(nt,3) and ei(nt) are double precision arrays of
!  dimension nt which is an integer representing the number of geometrie
!  input.  The integer nsurf selects the desired surface.  nsurf=1 is th
!  lower diabatic surface, nsurf=2 is the coupling surface, nsurf=3 is
!  the upper diabatic surface, nsurf = 4 selects the lower adiabatic
!  surface, and nsurf = 5 selects the upper adiabatic surface
!

!
!> @brief Evaluate pot for the NaH2_LEPS_7S_2000_DPEM potential.
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
module chempotpy_nah2_leps_7s_2000_dpem_pot11_legacy_state_ac46def9
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        implicit none
        optional :: r, e, nt
        double precision :: b0, bect, bect2, beh2, bes1, bf, cevau, cos1, cos2, cosa, cosg, dect, deh2, &
            & deh2m, des1, des20, e, gamh2a, gamh2b, onethd, r, r1s, r2s, r3s, rac2, rave, rect, reh2, &
            & res1, rmid, s1, s2, s2alp, s2rho, s3, sswt, twothd, vmid, xect, xeh2, xes1, xes3
        integer :: i, nt
        dimension r(:,:), e(:)
        ! NaH singlet
        ! H2 singlet
        ! H-Na-H term
        !  misc constants
        b0     = 0.d0
        bect   = 0.d0
        bect2  = 0.d0
        beh2   = 0.d0
        bes1   = 0.d0
        bf     = 0.d0
        cevau  = 0.d0
        cos1   = 0.d0
        cos2   = 0.d0
        cosa   = 0.d0
        cosg   = 0.d0
        dect   = 0.d0
        deh2   = 0.d0
        deh2m  = 0.d0
        des1   = 0.d0
        des20  = 0.d0
        gamh2a = 0.d0
        gamh2b = 0.d0
        onethd = 0.d0
        r1s    = 0.d0
        r2s    = 0.d0
        r3s    = 0.d0
        rac2   = 0.d0
        rave   = 0.d0
        rect   = 0.d0
        reh2   = 0.d0
        res1   = 0.d0
        rmid   = 0.d0
        s1     = 0.d0
        s2     = 0.d0
        s2alp  = 0.d0
        s2rho  = 0.d0
        s3     = 0.d0
        sswt   = 0.d0
        twothd = 0.d0
        vmid   = 0.d0
        xect   = 0.d0
        xeh2   = 0.d0
        xes1   = 0.d0
        xes3   = 0.d0
        des1   = 0.55118d0
        bes1   = 1.05118d0
        res1   = 3.38583d0
        reh2   = 1.401d0
        bf     = 1.627d0
        b0     = 1.194d0
        gamh2a = 2.323d0
        gamh2b = 3.126d0
        deh2   = 4.74806d0
        des20  = 0.996d0
        s2alp  = 3.0d0
        s2rho  = 7.0d0
        dect   = 2.24409d0
        bect   = 0.57874d0
        rect   = 3.46457d0
        bect2  = 0.87143d0
        rac2   = 1.d0/dsqrt(2.d0)
        cevau  = 1.d0/27.2113961d0
        onethd = 1.0d0/3.0d0
        twothd = 2.0d0/3.0d0
        do 10 i=1, nt
            ! ============================================================
            ! useful quantities calculated first
            r1s  = r(i,1)*r(i,1)
            r2s  = r(i,2)*r(i,2)
            r3s  = r(i,3)*r(i,3)
            rmid = 0.5d0*(r(i,1)+r(i,3))
            rave = 0.5d0*(r(i,1)-r(i,3))
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            cos1 = (r1s+r2s-r3s)/(2.0d0*r(i,1)*r(i,2))
            cos2 = (r3s+r2s-r1s)/(2.0d0*r(i,3)*r(i,2))
            cosg = 0.5d0*(cos1-cos2)
            cosa = (r1s+r3s-r2s)/(2*r(i,1)*r(i,3))
            !==============================================================k
            !This section contains the diatomic curves
            ! The first NaH singlet
            xes1 = dexp(-bes1*(r(i,1)-res1))
            s1   = des1*onethd*xes1*(xes1+2.0d0)
            ! The H2 singlet
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
            ! =========================================================
            ! here is where it all gets put together . . .
            e(i) = s1 + s2 + s3 + deh2
            ! here we will try a post-LEPS correction term
            xect = dexp(-bect*(rmid-rect))
            vmid = dect*xect*(xect-2.0d0) *dexp(-bect2*rave**2)*(0.5d0*(1.0d0-cosa))**2
            e(i) = e(i)+vmid
            e(i) = e(i)*cevau
        10 continue
    end subroutine dispatch
end module chempotpy_nah2_leps_7s_2000_dpem_pot11_legacy_state_ac46def9

!> @brief Compatibility interface for pot11.
subroutine pot11( r, e, nt )
    use chempotpy_nah2_leps_7s_2000_dpem_pot11_legacy_state_ac46def9, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    integer          :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot11

! U12 coupling surface . . .

!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_nah2_leps_7s_2000_dpem_pot12_legacy_state_ac46def9
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
end module chempotpy_nah2_leps_7s_2000_dpem_pot12_legacy_state_ac46def9

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, nt )
    use chempotpy_nah2_leps_7s_2000_dpem_pot12_legacy_state_ac46def9, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    integer          :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot12

! the upper diabatic surface  U22

!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_nah2_leps_7s_2000_dpem_pot22_legacy_state_ac46def9
    implicit none
    private
    public :: dispatch
    contains
    !> @brief Evaluate dispatch.
    subroutine dispatch( r, e, nt )
        implicit none
        optional :: r, e, nt
        double precision :: b0, beh2, bet1, bet10, bet10b, bet190, bet1c, bet20, bet290, beta, bf, c2a, &
            & c2b, c2c, cevau, cos1, cos2, cosa, cosg, coul1, coul2, coul3, cplg2, deh2, deh2m, denah, &
            & des290, det1, det190, det1c, det20, det290, e, exch1, exch2, exch3, f10, g10, g30, gamh2a, &
            & gamh2b, gamnah, onethd, r, r1s, r2s, r3s, rave, reh2, renah, ret1, ret10, ret190, ret1c, &
            & ret20, ret290, rmid, s1, s2, s20alp, s20b, s20rho, s2a, s2alp, s2rho, s3, sswt, t1, t10, &
            & t190, t1alp, t1c, t1m, t1rho, t2, t20, t290, t2alp, t2rho, t3, t30, t390, t3c, t3m, twothd, &
            & una2p, w, xeh2, xenah, xet190, xet1c, xet20, xet290, xet3, xet390, xet3c
        integer :: i, nt
        dimension r(:,:), e(:)
        ! NaH triplets
        ! H2 triplets
        ! H2 singlets
        ! misc constants
        b0     = 0.d0
        beh2   = 0.d0
        bet1   = 0.d0
        bet10  = 0.d0
        bet10b = 0.d0
        bet190 = 0.d0
        bet1c  = 0.d0
        bet20  = 0.d0
        bet290 = 0.d0
        beta   = 0.d0
        bf     = 0.d0
        c2a    = 0.d0
        c2b    = 0.d0
        c2c    = 0.d0
        cevau  = 0.d0
        cos1   = 0.d0
        cos2   = 0.d0
        cosa   = 0.d0
        cosg   = 0.d0
        coul1  = 0.d0
        coul2  = 0.d0
        coul3  = 0.d0
        cplg2  = 0.d0
        deh2   = 0.d0
        deh2m  = 0.d0
        denah  = 0.d0
        des290 = 0.d0
        det1   = 0.d0
        det190 = 0.d0
        det1c  = 0.d0
        det20  = 0.d0
        det290 = 0.d0
        exch1  = 0.d0
        exch2  = 0.d0
        exch3  = 0.d0
        f10    = 0.d0
        g10    = 0.d0
        g30    = 0.d0
        gamh2a = 0.d0
        gamh2b = 0.d0
        gamnah = 0.d0
        onethd = 0.d0
        r1s    = 0.d0
        r2s    = 0.d0
        r3s    = 0.d0
        rave   = 0.d0
        reh2   = 0.d0
        renah  = 0.d0
        ret1   = 0.d0
        ret10  = 0.d0
        ret190 = 0.d0
        ret1c  = 0.d0
        ret20  = 0.d0
        ret290 = 0.d0
        rmid   = 0.d0
        s1     = 0.d0
        s2     = 0.d0
        s20alp = 0.d0
        s20b   = 0.d0
        s20rho = 0.d0
        s2a    = 0.d0
        s2alp  = 0.d0
        s2rho  = 0.d0
        s3     = 0.d0
        sswt   = 0.d0
        t1     = 0.d0
        t10    = 0.d0
        t190   = 0.d0
        t1alp  = 0.d0
        t1c    = 0.d0
        t1m    = 0.d0
        t1rho  = 0.d0
        t2     = 0.d0
        t20    = 0.d0
        t290   = 0.d0
        t2alp  = 0.d0
        t2rho  = 0.d0
        t3     = 0.d0
        t30    = 0.d0
        t390   = 0.d0
        t3c    = 0.d0
        t3m    = 0.d0
        twothd = 0.d0
        una2p  = 0.d0
        w      = 0.d0
        xeh2   = 0.d0
        xenah  = 0.d0
        xet190 = 0.d0
        xet1c  = 0.d0
        xet20  = 0.d0
        xet290 = 0.d0
        xet3   = 0.d0
        xet390 = 0.d0
        xet3c  = 0.d0
        ret10  = 3.22835d0
        bet10  = 0.57874d0
        bet10b = 2.44094d0
        f10    = 9.44882d0
        det190 = 4.32441d0
        bet190 = 2.43701d0
        ret190 = 1.95480d0
        t1alp  = 0.30778d0
        t1rho  = 3.21850d0
        det1c  = 4.92126d0
        bet1c  = 0.65748d0
        ret1c  = 2.71654d0
        det20  = 3.29921d0
        bet20  = 2.64567d0
        ret20  = 2.44882d0
        det290 = 3.97638d0
        bet290 = 3.73937d0
        ret290 = 0.96024d0
        t2alp  = 1.0d0
        t2rho  = 9.0d0
        des290 = 0.91433d0
        s2alp  = 0.66778d0
        s2rho  = 5.30472d0
        deh2   = 4.74806d0
        s20alp = 0.5d0
        s20rho = 6.0d0
        c2a    = 1.38661d0
        c2b    = 0.27362d0
        c2c    = 0.03764d0
        una2p  = 2.1037d0
        cevau  = 1.d0/27.2113961d0
        onethd = 1.0d0/3.0d0
        twothd = 2.0d0/3.0d0
        do 10 i=1, nt
            ! ============================================================
            ! useful quantities calculated first
            rmid = 0.5d0*(r(i,1)+r(i,3))
            rave = 0.5d0*(r(i,1)-r(i,3))
            r1s  = r(i,1)*r(i,1)
            r2s  = r(i,2)*r(i,2)
            r3s  = r(i,3)*r(i,3)
            cosa = (r1s+r3s-r2s)/(2*r(i,1)*r(i,3))
            ! cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
            cos1 = (r1s+r2s-r3s)/(2.0d0*r(i,1)*r(i,2))
            cos2 = (r3s+r2s-r1s)/(2.0d0*r(i,3)*r(i,2))
            cosg = 0.5d0*(cos1-cos2)
            ! =========================================================
            ! The diatomic curves . . .
            ! The first NaH singlet
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
            ! ===============================================================
            ! The final LEPS form
            coul1 = 0.5d0*(s1+t1)
            coul2 = 0.5d0*(s2+t2)
            coul3 = 0.5d0*(s3+t3)
            exch1 = 0.5d0*(s1-t1)
            exch2 = 0.5d0*(s2-t2)
            exch3 = 0.5d0*(s3-t3)
            w     = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
            cplg2 = c2a*dexp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
            e(i)  = coul1+coul2+coul3+deh2-dsqrt(w+cplg2*cplg2)/dsqrt(2.0d0)
            e(i)  = e(i)*cevau
        10 continue
    end subroutine dispatch
end module chempotpy_nah2_leps_7s_2000_dpem_pot22_legacy_state_ac46def9

!> @brief Compatibility interface for pot22.
subroutine pot22( r, e, nt )
    use chempotpy_nah2_leps_7s_2000_dpem_pot22_legacy_state_ac46def9, only: dispatch
    implicit none
    double precision :: r(nt,3)
    double precision :: e(nt)
    integer          :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot22
