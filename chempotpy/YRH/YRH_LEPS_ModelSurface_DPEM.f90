!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate dpem for the YRH_LEPS_ModelSurface_DPEM potential.
subroutine dpem( x, igrad, u, ug )

    implicit none

    ! number of electronic state
    integer, parameter :: nstates = 2
    integer, parameter :: natoms  = 3

    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: u(nstates,nstates)
    double precision, intent(out) :: ug(nstates,nstates,natoms,3)

    double precision :: d, g
    integer          :: iatom, idir, nder
    double precision :: nt, r(1,3), r2(3), v(1)
    double precision :: p(nstates)
    double precision :: u11(1), u12(1), u22(1)
    double precision :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
    double precision :: dudr(nstates,nstates,3)
    double precision :: dudx(nstates,nstates,9)
    double precision :: t(nstates,nstates)
    double precision :: tmpmat(nstates,nstates)
    double precision :: tx(9), drdx(3,9)
    logical, save    :: first_time_data=.true.
    integer          :: i, j, k, l

    !initialize
    u = 0.d0
    p = 0.d0
    g = 0.d0
    d = 0.d0

    nt = 1
    do iatom=1, natoms
        do idir=1,3
            j     = (iatom-1)*3+idir
            tx(j) = x(iatom,idir)
        enddo
    enddo
    ! input cartesian is ClHH
    r(1,1) = sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2 +(x(1,3)-x(2,3))**2)/0.529177211
    r(1,2) = sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2 +(x(2,3)-x(3,3))**2)/0.529177211
    r(1,3) = sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2 +(x(1,3)-x(3,3))**2)/0.529177211

    nder = igrad

    if(first_time_data) then
        call prepot
        first_time_data = .false.
    endif

    call pot(r,u11,de_u11,1,1)
    call pot(r,u12,de_u12,1,2)
    call pot(r,u22,de_u22,1,3)
    u(1,1)      = u11(1)*27.211386
    u(1,2)      = u12(1)*27.211386
    u(2,1)      = u12(1)*27.211386
    u(2,2)      = u22(1)*27.211386
    dudr(1,1,:) = de_u11(:,1)*51.422067
    dudr(1,2,:) = de_u12(:,1)*51.422067
    dudr(2,1,:) = de_u12(:,1)*51.422067
    dudr(2,2,:) = de_u22(:,1)*51.422067
    call diagonalize(nstates,u,p,t)

    r2 = r(1,:)*0.529177211
    call evdrdx(tx, r2, drdx)
    dudx = 0.d0
    do i=1,9
        do j=1,3
            dudx(:,:,i) = dudx(:,:,i)+dudr(:,:,j)*drdx(j,i)
        enddo
    enddo

    do iatom=1, natoms
        do idir=1,3
            j                  = (iatom-1)*3+idir
            ug(:,:,iatom,idir) = dudx(:,:,j)
        enddo
    enddo

end subroutine dpem

!> @brief Evaluate EvdRdX for the YRH_LEPS_ModelSurface_DPEM potential.
subroutine evdrdx( x, r, drdx )

    implicit none

    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)

    integer :: i,j

    ! Initialize dRdX(3,9)
    do i=1,3
        do j=1,9
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

!> @brief Evaluate diagonalize for the YRH_LEPS_ModelSurface_DPEM potential.
subroutine diagonalize( n, a_ss, ev_s, u_ss )
    implicit none

    integer, intent(in)           :: n
    double precision, intent(in)  :: a_ss(n,n)
    double precision, intent(out) :: u_ss(n,n)
    double precision, intent(out) :: ev_s(n)
    integer                       :: io,i
    double precision,allocatable  :: work_d(:)
    integer, save                 :: lwork_d
    logical, save                 :: first_time_diag=.true.

    u_ss = a_ss
    if (first_time_diag) then
        lwork_d = -1
        allocate( work_d(2) )
        call dsyev('V','L',n,u_ss,n,ev_s,work_d,lwork_d,io)
        lwork_d = int(work_d(1))
        deallocate ( work_d )
        first_time_diag = .false.
    endif
    allocate( work_d(lwork_d) )
    call dsyev('V','L',n,u_ss,n,ev_s,work_d,lwork_d,io)
end subroutine diagonalize

!
!   System:          YRH
!   Common name:     YRH
!   Functional form: Modified Extended LEPS (2x2 diabatic fit)
!   Number of derivatives: 1
!   Number of electronic surfaces: 2
!   Interface: 3-2V

!   References:   A. W. Jasper, M. D. Hack, and D. G. Truhlar, J. Chem.
!                 Vol. 115, p. 1804, 2001.

!   Notes:  This is a model surface set where the diabatic surfaces are
!           weakly coupled and do not cross.  The parameter SCALE determ
!           the strength of the diabatic coupling.
!
!   Protocol:
!      PREPOT - initializes the potential's variables
!               must be called once before any calls to POT
!      POT    - driver for the evaluation of the energy
!
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!
!   Surfaces:
!
!   Zero of energy:
!      The classical potential energy is set equal to zero for the Y
!      infinitely far from the RH diatomic and R(RH) set equal to the
!      RH equilibrium diatomic value.
!
!   Parameters:
!
!   Coordinates:
!      Internal, Definition: R(1) = R(Y-R)
!                            R(2) = R(R-H)
!                            R(3) = R(H-Y)

!
!
!  The potential is called by first issuing a call to PREPOT which initi
!  all surfaces.  The potential is called by calling POT with the parame
!  r, e, nt, nsurf.  r(nt,3) and ei(nt) are double precision arrays of
!  dimension nt which is an integer representing the number of geometrie
!  input.  The integer nsurf selects the desired surface.  nsurf=1 is th
!  lower diabatic surface, nsurf=2 is the coupling surface, nsurf=3 is
!  the upper diabatic surface.

!*******************************************

!> @brief Evaluate prepot for the YRH_LEPS_ModelSurface_DPEM potential.
subroutine prepot
    implicit none

    integer, parameter :: g_p_ = 3
    call g_prepot11(g_p_)
    call g_prepot12(g_p_)
    call g_prepot22(g_p_)
end subroutine prepot

!*******************************************

!> @brief Evaluate pot for the YRH_LEPS_ModelSurface_DPEM potential.
subroutine pot( ro, e, de, nt, nsurf )
    !      subroutine pot(ro,nt,nsurf)

    implicit none

    double precision :: de, deo, e, g_e, g_r, r, ro
    integer          :: i, j, k, nsurf, nt
    dimension r(nt,3),de(3,nt),e(nt),ro(nt,3)
    dimension deo(3,nt)

    integer, parameter :: g_p_  = 3
    integer, parameter :: ldg_r = 3
    integer, parameter :: ldg_e = 3
    dimension g_r(ldg_r,nt,3),g_e(ldg_e,nt)

    !     permute the bonds
    !     want r1 = YR
    !          r2 = RH
    !          r3 = YH

    !     so arr #1 = Y + RH
    !            #2 = R + YH
    !            #3 = H + YR

    !     in code:  r1 = RH
    !               r2 = YH
    !               r3 = YR

    !     so...

    do i=1,nt
        r(i,1) = ro(i,2)
        r(i,2) = ro(i,3)
        r(i,3) = ro(i,1)
    enddo

    !     seed matrix

    do i=1,nt
        do j=1,3
            do k=1,3
                if (j .eq. k) g_r(j,i,k)=1.0d0
                if (j .ne. k) g_r(j,i,k)=0.0d0
            enddo
        enddo
    enddo

    if (nsurf.eq.1) then
        call g_pot11(g_p_, r, g_r, ldg_r, e, deo, ldg_e, nt)
    else if (nsurf.eq.2) then
        call g_pot12(g_p_, r, g_r, ldg_r, e, deo, ldg_e, nt)
    else
        call g_pot22(g_p_, r, g_r, ldg_r, e, deo, ldg_e, nt)
    endif

    do i=1,nt
        de(1,i) = deo(3,i)
        de(2,i) = deo(1,i)
        de(3,i) = deo(2,i)
    enddo

end subroutine pot

!******************************************************************
!******************************************************************

!                           DISCLAIMER
!
!   This file was generated on 08/09/99 by the version of
!   ADIFOR compiled on June, 1998.
!
!   ADIFOR was prepared as an account of work sponsored by an
!   agency of the United States Government, Rice University, and
!   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
!   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
!   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
!   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
!   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
!   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
!   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
!
!*******************************************
! U11- lower diabat
!
!> @brief Evaluate g_prepot11 for the YRH_LEPS_ModelSurface_DPEM potential.
subroutine g_prepot11( g_p_ )
    !
    implicit none

    double precision :: beta1_11, beta2_11, beta3_11, c1, ca_11, cb_11, cc_11, cevau, coul1, coul2, &
        & coul3, de1_11, de2_11, de3_11, e, exch1, exch2, exch3, r, r1, r2, r3, re1_11, re2_11, &
        & re3_11, rootw, s1, s2, s3, sp1_11, sp2_11, sp3_11, t1, t2, t3, w
    integer :: i, nt
    dimension r(nt, 3), e(nt)
    !
    !     data is in bohr, and eV units                                    x
    !
    save
    integer, parameter :: g_pmax_ = 3
    integer            :: g_i_, g_p_, ldg_r, ldg_e
    double precision   :: d3_p, d7_b, d10_b, d7_v, d2_p, d5_b, d1_w, d10_&
        &v, d2_v, d6_b
    double precision :: d1_p, d2_b, d3_b, d4_v, g_r1(g_pmax_), g_r(ldg_&
        &r, nt, 3), g_r2(g_pmax_), g_r3(g_pmax_), g_d1_w(g_pmax_), g_s1(g_p&
        &max_)
    double precision :: g_s2(g_pmax_), g_s3(g_pmax_), g_t1(g_pmax_), g_&
        &t2(g_pmax_), g_t3(g_pmax_), g_coul1(g_pmax_), g_coul2(g_pmax_), g_&
        &coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_)
    double precision :: g_exch3(g_pmax_), g_w(g_pmax_), g_c1(g_pmax_), &
        & g_rootw(g_pmax_), g_e(ldg_e, nt)
    integer :: g_ehfid
    save g_coul1, g_coul2, g_coul3, g_exch1, g_exch2, g_exch3, g_w, &
        & g_c1, g_rootw
    save g_r1, g_r2, g_r3, g_d1_w, g_s1, g_s2, g_s3, g_t1, g_t2, g_t&
        &3
    data de1_11 /3.9d0/, de2_11 /4.3d0/, de3_11 /0.4d0/
    data re1_11 /2.1d0/, re2_11 /2.1d0/, re3_11 /2.5d0/
    data beta1_11 /1.0d0/, beta2_11 /1.0d0/, beta3_11 /1.5d0/
    data ca_11 /1.5d0/, cb_11 /1.0d0/, cc_11 /0.15d0/
    data sp1_11 /0.9d0/, sp2_11 /0.9d0/, sp3_11 /0.2d0/
    data g_ehfid /0/
    !
    !
    if (g_p_ .gt. g_pmax_) then
        print *, 'Parameter g_p_ is greater than g_pmax_'
        stop
    endif
    cevau = 1.d0 / 27.2113961d0
    !
    return
    !
    entry g_pot11(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
    !
    do i = 1, nt
        !
        do g_i_ = 1, g_p_
            g_r1(g_i_) = g_r(g_i_, i, 1)
        enddo
        r1 = r(i, 1)
        !--------
        do g_i_ = 1, g_p_
            g_r2(g_i_) = g_r(g_i_, i, 2)
        enddo
        r2 = r(i, 2)
        !--------
        do g_i_ = 1, g_p_
            g_r3(g_i_) = g_r(g_i_, i, 3)
        enddo
        r3 = r(i, 3)
        !--------
        !
        !   SINGLETS (s1, s2, s3)
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta1_11) * g_r1(g_i_)
        enddo
        d1_w = (-beta1_11) * (r1 - re1_11)
        d2_v = exp(d1_w)
        d2_p = d2_v
        d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
        d1_p = 2.0d0 * (1.d0 - d2_v)
        d6_b = (-(de1_11 * d1_p)) * d2_p
        do g_i_ = 1, g_p_
            g_s1(g_i_) = d6_b * g_d1_w(g_i_)
        enddo
        s1 = de1_11 * d4_v - de1_11
        !--------
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta2_11) * g_r2(g_i_)
        enddo
        d1_w = (-beta2_11) * (r2 - re2_11)
        d2_v = exp(d1_w)
        d2_p = d2_v
        d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
        d1_p = 2.0d0 * (1.d0 - d2_v)
        d6_b = (-(de2_11 * d1_p)) * d2_p
        do g_i_ = 1, g_p_
            g_s2(g_i_) = d6_b * g_d1_w(g_i_)
        enddo
        s2 = de2_11 * d4_v - de2_11
        !--------
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta3_11) * g_r3(g_i_)
        enddo
        d1_w = (-beta3_11) * (r3 - re3_11)
        d2_v = exp(d1_w)
        d2_p = d2_v
        d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
        d1_p = 2.0d0 * (1.d0 - d2_v)
        d6_b = (-(de3_11 * d1_p)) * d2_p
        do g_i_ = 1, g_p_
            g_s3(g_i_) = d6_b * g_d1_w(g_i_)
        enddo
        s3 = de3_11 * d4_v - de3_11
        !--------
        !
        !   TRIPLETS (t1, t2, t3)
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta1_11) * g_r1(g_i_)
        enddo
        d1_w = (-beta1_11) * (r1 - re1_11)
        d2_v = exp(d1_w)
        d2_p = d2_v
        d4_v = (1.d0 + d2_v) * (1.d0 + d2_v)
        d1_p = 2.0d0 * (1.d0 + d2_v)
        d7_b = 0.5d0 * de1_11 * d1_p * d2_p
        do g_i_ = 1, g_p_
            g_t1(g_i_) = d7_b * g_d1_w(g_i_)
        enddo
        t1 = 0.5d0 * (de1_11 * d4_v - de1_11)
        !--------
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta2_11) * g_r2(g_i_)
        enddo
        d1_w = (-beta2_11) * (r2 - re2_11)
        d2_v = exp(d1_w)
        d2_p = d2_v
        d4_v = (1.d0 + d2_v) * (1.d0 + d2_v)
        d1_p = 2.0d0 * (1.d0 + d2_v)
        d7_b = 0.5d0 * de2_11 * d1_p * d2_p
        do g_i_ = 1, g_p_
            g_t2(g_i_) = d7_b * g_d1_w(g_i_)
        enddo
        t2 = 0.5d0 * (de2_11 * d4_v - de2_11)
        !--------
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta3_11) * g_r3(g_i_)
        enddo
        d1_w = (-beta3_11) * (r3 - re3_11)
        d2_v = exp(d1_w)
        d2_p = d2_v
        d4_v = (1.d0 + d2_v) * (1.d0 + d2_v)
        d1_p = 2.0d0 * (1.d0 + d2_v)
        d7_b = 0.5d0 * de3_11 * d1_p * d2_p
        do g_i_ = 1, g_p_
            g_t3(g_i_) = d7_b * g_d1_w(g_i_)
        enddo
        t3 = 0.5d0 * (de3_11 * d4_v - de3_11)
        !--------
        !
        do g_i_ = 1, g_p_
            g_t1(g_i_) = sp1_11 * g_t1(g_i_)
        enddo
        t1 = t1 * sp1_11
        !--------
        do g_i_ = 1, g_p_
            g_t2(g_i_) = sp2_11 * g_t2(g_i_)
        enddo
        t2 = t2 * sp2_11
        !--------
        do g_i_ = 1, g_p_
            g_t3(g_i_) = sp3_11 * g_t3(g_i_)
        enddo
        t3 = t3 * sp3_11
        !--------
        !
        !  EXCHANGE AND COULOMB  (exch1, coul1, etc.)
        do g_i_ = 1, g_p_
            g_coul1(g_i_) = 0.5d0 * g_t1(g_i_) + 0.5d0 * g_s1(g_i_)
        enddo
        coul1 = 0.5d0 * (s1 + t1)
        !--------
        do g_i_ = 1, g_p_
            g_coul2(g_i_) = 0.5d0 * g_t2(g_i_) + 0.5d0 * g_s2(g_i_)
        enddo
        coul2 = 0.5d0 * (s2 + t2)
        !--------
        do g_i_ = 1, g_p_
            g_coul3(g_i_) = 0.5d0 * g_t3(g_i_) + 0.5d0 * g_s3(g_i_)
        enddo
        coul3 = 0.5d0 * (s3 + t3)
        !--------
        do g_i_ = 1, g_p_
            g_exch1(g_i_) = (-0.5d0) * g_t1(g_i_) + 0.5d0 * g_s1(g_i_)
        enddo
        exch1 = 0.5d0 * (s1 - t1)
        !--------
        do g_i_ = 1, g_p_
            g_exch2(g_i_) = (-0.5d0) * g_t2(g_i_) + 0.5d0 * g_s2(g_i_)
        enddo
        exch2 = 0.5d0 * (s2 - t2)
        !--------
        do g_i_ = 1, g_p_
            g_exch3(g_i_) = (-0.5d0) * g_t3(g_i_) + 0.5d0 * g_s3(g_i_)
        enddo
        exch3 = 0.5d0 * (s3 - t3)
        !--------
        !
        d4_v  = (exch1 - exch2) * (exch1 - exch2)
        d3_p  = 2.0d0 * (exch1 - exch2)
        d7_v  = (exch2 - exch3) * (exch2 - exch3)
        d2_p  = 2.0d0 * (exch2 - exch3)
        d10_v = (exch3 - exch1) * (exch3 - exch1)
        d1_p  = 2.0d0 * (exch3 - exch1)
        d5_b  = d1_p + (-d2_p)
        d6_b  = -d1_p + d3_p
        d10_b = d2_p + (-d3_p)
        do g_i_ = 1, g_p_
            g_w(g_i_) = d5_b * g_exch3(g_i_) + d10_b * g_exch2(g_i_) + d6_b * g_exch1(g_i_)
        enddo
        w = d4_v + d7_v + d10_v
        !--------
        !
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-cc_11) * g_r3(g_i_) + (-cc_11) * g_r2(g_i_) + (-cc_11) * g_r1(g_i_) + (-cb_11) * g_w(g_i_)
        enddo
        d1_w = (-cb_11) * w - cc_11 * (r1 + r2 + r3)
        d2_v = exp(d1_w)
        d1_p = d2_v
        d3_b = ca_11 * d1_p
        do g_i_ = 1, g_p_
            g_c1(g_i_) = d3_b * g_d1_w(g_i_)
        enddo
        c1 = ca_11 * d2_v
        !--------
        !
        d2_b = 1.0d0 / 2.d0
        d5_b = d2_b * c1 + d2_b * c1
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_c1(g_i_) + d2_b * g_w(g_i_)
        enddo
        d1_w = (w + c1 * c1) / 2.d0
        d2_v = sqrt(d1_w)

        if ( d1_w .gt. 0.0d0 ) then
            d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
            d1_p = 0.d0
        endif
        do g_i_ = 1, g_p_
            g_rootw(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        rootw = d2_v
        !--------
        !
        do g_i_ = 1, g_p_
            g_e(g_i_, i) = -g_rootw(g_i_) + g_coul3(g_i_) + g_coul2(g_i_ ) + g_coul1(g_i_)
        enddo
        e(i) = coul1 + coul2 + coul3 - rootw + de2_11
        !--------
        do g_i_ = 1, g_p_
            g_e(g_i_, i) = cevau * g_e(g_i_, i)
        enddo
        e(i) = e(i) * cevau
        !--------
        !
        !  end the nt loop here:
    enddo
    !
end subroutine g_prepot11
!
!**************************************************

!**************************************************

!  U22- upper diabat
!
!> @brief Evaluate g_prepot22 for the YRH_LEPS_ModelSurface_DPEM potential.
subroutine g_prepot22( g_p_ )
    !
    implicit none

    double precision :: beta1_22, beta2_22, beta3_22, cevau, de1_22, de2_11, de2_22, de3_22, e, &
        & excite, r, r1, r2, r3, re1_22, re2_22, re3_22, s1, t2, t3, wgt
    integer :: i, nt
    dimension r(nt, 3), e(nt)
    !
    !     data is in bohr and eV units
    !
    save
    integer, parameter :: g_pmax_ = 3
    integer            :: g_i_, g_p_, ldg_r, ldg_e
    double precision   :: d6_b, d2_v, d7_b, d4_v, d2_p, d1_p, d1_w, g_r1( &
        & g_pmax_), g_r(ldg_r, nt, 3), g_r2(g_pmax_)
    double precision :: g_r3(g_pmax_), g_d1_w(g_pmax_), g_s1(g_pmax_), &
        & g_t2(g_pmax_), g_t3(g_pmax_), g_e(ldg_e, nt)
    integer :: g_ehfid
    save g_r1, g_r2, g_r3, g_d1_w, g_s1, g_t2, g_t3
    data de2_11 /4.3d0/, excite /0.36d0/, wgt /0.2d0/
    data de1_22 /3.9d0/, de2_22 /4.3d0/, de3_22 /0.4d0/
    data re1_22 /2.1d0/, re2_22 /2.1d0/, re3_22 /2.5d0/
    data beta1_22 /1.d0/, beta2_22 /1.0d0/, beta3_22 /1.5d0/
    data g_ehfid /0/
    !
    !
    if (g_p_ .gt. g_pmax_) then
        print *, 'Parameter g_p_ is greater than g_pmax_'
        stop
    endif
    cevau = 1.d0 / 27.2113961d0
    !
    return
    !
    entry g_pot22(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
    !
    do i = 1, nt
        !
        do g_i_ = 1, g_p_
            g_r1(g_i_) = g_r(g_i_, i, 1)
        enddo
        r1 = r(i, 1)
        !--------
        do g_i_ = 1, g_p_
            g_r2(g_i_) = g_r(g_i_, i, 2)
        enddo
        r2 = r(i, 2)
        !--------
        do g_i_ = 1, g_p_
            g_r3(g_i_) = g_r(g_i_, i, 3)
        enddo
        r3 = r(i, 3)
        !--------
        !
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta1_22) * g_r1(g_i_)
        enddo
        d1_w = (-beta1_22) * (r1 - re1_22)
        d2_v = exp(d1_w)
        d2_p = d2_v
        d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
        d1_p = 2.0d0 * (1.d0 - d2_v)
        d6_b = (-(de1_22 * d1_p)) * d2_p
        do g_i_ = 1, g_p_
            g_s1(g_i_) = d6_b * g_d1_w(g_i_)
        enddo
        s1 = de1_22 * d4_v - de1_22
        !--------
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta2_22) * g_r2(g_i_)
        enddo
        d1_w = (-beta2_22) * (r2 - re2_22)
        d2_v = exp(d1_w)
        d2_p = d2_v
        d4_v = (1.d0 + d2_v) * (1.d0 + d2_v)
        d1_p = 2.0d0 * (1.d0 + d2_v)
        d7_b = 0.5d0 * de2_22 * d1_p * d2_p
        do g_i_ = 1, g_p_
            g_t2(g_i_) = d7_b * g_d1_w(g_i_)
        enddo
        t2 = 0.5d0 * (de2_22 * d4_v - de2_22)
        !--------
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta3_22) * g_r3(g_i_)
        enddo
        d1_w = (-beta3_22) * (r3 - re3_22)
        d2_v = exp(d1_w)
        d2_p = d2_v
        d4_v = (1.d0 + d2_v) * (1.d0 + d2_v)
        d1_p = 2.0d0 * (1.d0 + d2_v)
        d7_b = 0.5d0 * de3_22 * d1_p * d2_p
        do g_i_ = 1, g_p_
            g_t3(g_i_) = d7_b * g_d1_w(g_i_)
        enddo
        t3 = 0.5d0 * (de3_22 * d4_v - de3_22)
        !--------
        !
        do g_i_ = 1, g_p_
            g_e(g_i_, i) = wgt * g_t3(g_i_) + wgt * g_t2(g_i_) + g_s1(g_i_)
        enddo
        e(i) = s1 + (t2 + t3) * wgt + excite + de2_11
        !--------
        do g_i_ = 1, g_p_
            g_e(g_i_, i) = cevau * g_e(g_i_, i)
        enddo
        e(i) = e(i) * cevau
        !--------
        !
        !  end the nt loop here:
    enddo
    !
end subroutine g_prepot22
!
!*****************************************

!  U12- diabatic coupling
!
!> @brief Evaluate g_prepot12 for the YRH_LEPS_ModelSurface_DPEM potential.
subroutine g_prepot12( g_p_ )
    !
    implicit none

    double precision :: a1, a2, cevau, csth, e, r, r1, r2, r3, re1_12, re2_12, s1, s2, scale, sith, &
        & thrad, thshift, tmp1, tmp2
    integer :: i, nt
    !      dimension r(nt,3),De(3,nt),e(nt)
    dimension r(nt, 3), e(nt)
    !
    save
    integer, parameter :: g_pmax_ = 3
    integer            :: g_i_, g_p_, ldg_r, ldg_e
    double precision   :: d1_p, d1_w, d4_b, d3_b, d2_v, g_r1(g_pmax_), g_&
        &r(ldg_r, nt, 3), g_r2(g_pmax_), g_s1(g_pmax_), g_s2(g_pmax_)
    double precision :: g_d1_w(g_pmax_), g_tmp1(g_pmax_), g_tmp2(g_pmax&
        &_), g_e(ldg_e, nt)
    integer :: g_ehfid
    save g_r1, g_r2, g_s1, g_s2, g_d1_w, g_tmp1, g_tmp2
    data a1 /0.8d0/, a2 /1.0d0/
    data re1_12 /2.2d0/, re2_12 /4.2d0/
    data thshift /45.d0/
    data scale /0.2d0/
    data g_ehfid /0/
    !
    !
    if (g_p_ .gt. g_pmax_) then
        print *, 'Parameter g_p_ is greater than g_pmax_'
        stop
    endif
    cevau = 1.d0 / 27.2113961d0
    !
    return
    !
    entry g_pot12(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
    !
    do i = 1, nt
        !
        do g_i_ = 1, g_p_
            g_r1(g_i_) = g_r(g_i_, i, 1)
        enddo
        r1 = r(i, 1)
        !--------
        do g_i_ = 1, g_p_
            g_r2(g_i_) = g_r(g_i_, i, 2)
        enddo
        r2 = r(i, 2)
        !--------
        r3 = r(i, 3)
        !
        thrad = dacos(-1.d0) / 180.d0 * thshift
        !
        csth = dcos(thrad)
        sith = dsin(thrad)
        !
        do g_i_ = 1, g_p_
            g_s1(g_i_) = sith * g_r2(g_i_) + csth * g_r1(g_i_)
        enddo
        s1 = (r1 - re1_12) * csth + (r2 - re2_12) * sith
        !--------
        do g_i_ = 1, g_p_
            g_s2(g_i_) = (-csth) * g_r2(g_i_) + sith * g_r1(g_i_)
        enddo
        s2 = (r1 - re1_12) * sith - (r2 - re2_12) * csth
        !--------
        !
        d2_v = s1 * s1
        d1_p = 2.0d0 * s1
        d3_b = (-a1) * d1_p
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_s1(g_i_)
        enddo
        d1_w = (-a1) * d2_v
        d2_v = exp(d1_w)
        d1_p = d2_v
        do g_i_ = 1, g_p_
            g_tmp1(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        tmp1 = d2_v
        !--------
        d2_v = s2 * s2
        d1_p = 2.0d0 * s2
        d3_b = (-a2) * d1_p
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_s2(g_i_)
        enddo
        d1_w = (-a2) * d2_v
        d2_v = exp(d1_w)
        d1_p = d2_v
        do g_i_ = 1, g_p_
            g_tmp2(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        tmp2 = d2_v
        !--------
        !
        d3_b = scale * tmp2
        d4_b = scale * tmp1
        do g_i_ = 1, g_p_
            g_e(g_i_, i) = d4_b * g_tmp2(g_i_) + d3_b * g_tmp1(g_i_)
        enddo
        e(i) = tmp1 * tmp2 * scale
        !--------
        !
        do g_i_ = 1, g_p_
            g_e(g_i_, i) = cevau * g_e(g_i_, i)
        enddo
        e(i) = e(i) * cevau
        !--------
        !
        !   end the nt loop here:
    enddo
    !
end subroutine g_prepot12
!*****************************************
!
!
