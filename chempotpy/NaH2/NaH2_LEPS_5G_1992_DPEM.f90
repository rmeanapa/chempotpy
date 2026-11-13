!> @brief Implement NaH2_LEPS_5G_1992, a diabatic potential matrix for NaH2.
module chempotpy_nah2_leps_5g_1992_dpem_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NaH2_LEPS_5G_1992 diabatic potential matrix for NaH2.
    subroutine evaluate( x, igrad, u, ug )
        use chempotpy_units,              only: hartree_to_electronvolt
        use chempotpy_legacy_coordinates, only: triatomic_distances_to_bohr
        use chempotpy_kinds,              only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 2
        integer, parameter    :: natoms  = 3
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: u(nstates,nstates)
        real(wp), intent(out) :: ug(nstates,nstates,natoms,3)
        real(wp)              :: nt, r(1,3), v(1)
        real(wp)              :: u11(1), u12(1), u22(1)
        integer               :: i, j, k, l
        !initialize
        u  = 0.e0_wp; ug = 0.e0_wp
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
        else
            write(*, *) "Only DPEM is available"
        endif
    end subroutine evaluate
end module chempotpy_nah2_leps_5g_1992_dpem_pes

!> @brief Expose NaH2_LEPS_5G_1992 through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_nah2_leps_5g_1992_dpem_pes, only: evaluate
    use chempotpy_kinds,                      only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 2
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: u(nstates,nstates)
    real(wp), intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem

!> @brief Evaluate EvdRdX for the NaH2_LEPS_5G_1992_DPEM potential.
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

!> @brief Evaluate diagonalize for the NaH2_LEPS_5G_1992_DPEM potential.
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

!> @brief Evaluate pot for the NaH2_LEPS_5G_1992_DPEM potential.
subroutine pot( r, e, nt, nsurf )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp) :: e(nt), r(nt,3)
    integer  :: nsurf, nt
    if( nsurf .eq. 1 )then
        call pot11(r, e, nt)
    elseif( nsurf .eq. 2 )then
        call pot12(r, e, nt)
    else
        call pot22(r, e, nt)
    endif
end subroutine pot

!  U11 surface - lower diabatic surface
!> @brief Own the stateless g_pot11 evaluator.
module chempotpy_nah2_leps_5g_1992_dpem_pot11_legacy_state_dd307a13
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
        real(wp) :: abs1, ads1, ars1, as1, as3, axs1, axs3, betat2, c1, c2, c2a, c2b, c2c, c3, c4, cbs1, cds1, cevau, &
            & cfm, coul1, coul2, coul3, cplg2, crs1, cs1, cs3, cstsq, cxs1, cxs3, deh2, dt2, e(:), exch1, exch2, exch3, r(:,:), &
            & r1s, r2s, r3s, rac2, ret2, rnag, s1, s2, s3, swa, swa2, swr0, swr20, swt, t1, t2, t3, w, xt2, y0c1, y0c2, y0c3, &
            & y0c4, y90c1, y90c2, y90c3, y90c4
        integer :: i, nt
        abs1   = 0.e0_wp
        ads1   = 0.e0_wp
        ars1   = 0.e0_wp
        as1    = 0.e0_wp
        as3    = 0.e0_wp
        axs1   = 0.e0_wp
        axs3   = 0.e0_wp
        betat2 = 0.e0_wp
        c1     = 0.e0_wp
        c2     = 0.e0_wp
        c2a    = 0.e0_wp
        c2b    = 0.e0_wp
        c2c    = 0.e0_wp
        c3     = 0.e0_wp
        c4     = 0.e0_wp
        cbs1   = 0.e0_wp
        cds1   = 0.e0_wp
        cevau  = 0.e0_wp
        cfm    = 0.e0_wp
        coul1  = 0.e0_wp
        coul2  = 0.e0_wp
        coul3  = 0.e0_wp
        cplg2  = 0.e0_wp
        crs1   = 0.e0_wp
        cs1    = 0.e0_wp
        cs3    = 0.e0_wp
        cstsq  = 0.e0_wp
        cxs1   = 0.e0_wp
        cxs3   = 0.e0_wp
        deh2   = 0.e0_wp
        dt2    = 0.e0_wp
        exch1  = 0.e0_wp
        exch2  = 0.e0_wp
        exch3  = 0.e0_wp
        r1s    = 0.e0_wp
        r2s    = 0.e0_wp
        r3s    = 0.e0_wp
        rac2   = 0.e0_wp
        ret2   = 0.e0_wp
        rnag   = 0.e0_wp
        s1     = 0.e0_wp
        s2     = 0.e0_wp
        s3     = 0.e0_wp
        swa    = 0.e0_wp
        swa2   = 0.e0_wp
        swr0   = 0.e0_wp
        swr20  = 0.e0_wp
        swt    = 0.e0_wp
        t1     = 0.e0_wp
        t2     = 0.e0_wp
        t3     = 0.e0_wp
        w      = 0.e0_wp
        xt2    = 0.e0_wp
        y0c1   = 0.e0_wp
        y0c2   = 0.e0_wp
        y0c3   = 0.e0_wp
        y0c4   = 0.e0_wp
        y90c1  = 0.e0_wp
        y90c2  = 0.e0_wp
        y90c3  = 0.e0_wp
        y90c4  = 0.e0_wp
        dt2    = 2.885656e-2_wp
        ret2   = 1.237624e0_wp
        betat2 = 1.038036e+1_wp
        y90c1  = 6.279847e0_wp
        y90c2  = 3.849003e-3_wp
        y90c3  = 5.546257e0_wp
        y90c4  = 5.573967e-1_wp
        y0c1   = 1.060311e+1_wp
        y0c2   = 6.058605e-3_wp
        y0c3   = 4.573378e0_wp
        y0c4   = 5.182708e-1_wp
        ads1   = 0.08e0_wp
        ars1   = 5.015440e0_wp
        abs1   = 7.536929e-1_wp
        cds1   = 1.323042e0_wp
        crs1   = 4.00e0_wp
        cbs1   = 5.4070e-1_wp
        swa    = 0.7e0_wp
        swr0   = 6.5e0_wp
        swa2   = 1.1e0_wp
        swr20  = 1.9e0_wp
        c2a    = 1.5e0_wp
        c2b    = 1.0e0_wp
        c2c    = 0.15e0_wp
        deh2   = 4.7469e0_wp
        cfm    = 1.8427829108061610e-4_wp
        rac2   = 1.e0_wp/sqrt(2.e0_wp)
        cevau  = 1.e0_wp/27.2113961e0_wp
        write(6, '(/24(''=''),'' NaH2 - U11 potential '',24(''=''), /''   dt2 ='',1pe14.6,2x,''  ret2 ='',e14.6,2x,''betat' // &
            '2 ='',e14.6, /'' y90c1 ='',e14.6,2x,'' y90c2 ='',e14.6,2x,'' y90c3 ='',e14.6, /'' y9' // &
            '0c4 ='',e14.6,2x,''  y0c1 ='',e14.6,2x,''  y0c2 ='',e14.6, /''  y0c3 ='',e14.6,2x,''' // &
            '  y0c4 ='',e14.6,2x,''  ads1 ='',e14.6, /''  ars1 ='',e14.6,2x,''  abs1 ='',e14.6,2x' // &
            ',''  cds1 ='',e14.6, /''  crs1 ='',e14.6,2x,''  cbs1 ='',e14.6,2x,''   swa ='',e14.6' // &
            ', /''  swr0 ='',e14.6,2x,''  swa2 ='',e14.6,2x,'' swr20 ='',e14.6, /''   c2a ='',e14' // &
            '.6,2x,''   c2b ='',e14.6,2x,''   c2c ='',e14.6, /70(''=''))') dt2, ret2, betat2, y90c1, y90c2, y90c3, y90c4, y0c1, y0c2, &
            & y0c3, y0c4, ads1, ars1, abs1, cds1, crs1, cbs1, swa, swr0, swa2, swr20, &
            & c2a, c2b, c2c
        do i=1, nt
            r1s  = r(i,1)*r(i,1)
            r2s  = r(i,2)*r(i,2)
            r3s  = r(i,3)*r(i,3)
            rnag = sqrt(0.5e0_wp*abs(r1s+r3s-0.5e0_wp*r2s))
            if( rnag .lt. 1.e-10_wp ) stop 'RNAG < 1.D-10'
            if( r(i,2) .ne. 0.e0_wp )then
                cstsq = (0.5e0_wp*(r1s-r3s)/(r(i,2)*rnag))**2
            else
                cstsq = 0.e0_wp
            endif
            c1 = y90c1+(y0c1-y90c1)*cstsq
            c2 = y90c2+(y0c2-y90c2)*cstsq
            c3 = y90c3+(y0c3-y90c3)*cstsq
            c4 = y90c4+(y0c4-y90c4)*cstsq
            swt=0.5e0_wp*(1.e0_wp-tanh(swa*(r(i,1)-swr0))) *0.5e0_wp*(1.e0_wp-tanh(swa*(r(i,3)-swr0))) *0.5e0_wp* &
                & (1.e0_wp+ tanh(swa2* (r(i,2)-swr20)))
            axs1 = exp(-abs1*(r(i,1)-ars1))
            as1  = ads1*axs1*(axs1-2.e0_wp)
            cxs1 = exp(-cbs1*(r(i,1)-crs1))
            cs1  = cds1*cxs1*(cxs1-2.e0_wp)
            s1   = as1+(cs1-as1)*swt
            s2=(149.480e0_wp-59.6557e0_wp*r(i,2))*exp(-4.13792e0_wp*r(i,2)) +r2s*(-23.7299e0_wp+3.91747e0_wp*r(i,2))* &
                & exp(- 1.41350e0_wp* r(i,2))
            axs3  = exp(-abs1*(r(i,3)-ars1))
            as3   = ads1*axs3*(axs3-2.e0_wp)
            cxs3  = exp(-cbs1*(r(i,3)-crs1))
            cs3   = cds1*cxs3*(cxs3-2.e0_wp)
            s3    = as3+(cs3-as3)*swt
            t1    = c1*exp(-c2*r(i,1)**6)+c3*exp(-c4*r(i,1))
            xt2   = exp(-betat2*(r(i,2)-ret2))
            t2    = dt2*xt2*(xt2-2.e0_wp)
            t3    = c1*exp(-c2*r(i,3)**6)+c3*exp(-c4*r(i,3))
            coul1 = 0.5e0_wp*(s1+t1)
            coul2 = 0.5e0_wp*(s2+t2)
            coul3 = 0.5e0_wp*(s3+t3)
            exch1 = 0.5e0_wp*(s1-t1)
            exch2 = 0.5e0_wp*(s2-t2)
            exch3 = 0.5e0_wp*(s3-t3)
            w     = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
            cplg2 = c2a*exp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
            e(i)  = (coul1+coul2+coul3+deh2-rac2*sqrt(w+cplg2*cplg2))*cevau-cfm
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_5g_1992_dpem_pot11_legacy_state_dd307a13

!> @brief Compatibility interface for pot11.
subroutine pot11( r, e, nt )
    use chempotpy_nah2_leps_5g_1992_dpem_pot11_legacy_state_dd307a13, only: dispatch
    use chempotpy_kinds,                                              only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot11

!  U12 surface - coupling surface
!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_nah2_leps_5g_1992_dpem_pot12_legacy_state_dd307a13
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
        real(wp) :: amp, ampe, ampn, ampp, ampr, ang, cevau, e(:), eps, erl, erle, erln, erlp, erlr, ers, erse, ersn, &
            & ersp, ersr, f1, f3, ff, fg, fh, hr, r(:,:), rmx, rmxe, rmxn, rmxp, rmxr, rnah16, rnah36
        integer :: i, nt
        amp    = 0.e0_wp
        ampe   = 0.e0_wp
        ampn   = 0.e0_wp
        ampp   = 0.e0_wp
        ampr   = 0.e0_wp
        ang    = 0.e0_wp
        cevau  = 0.e0_wp
        eps    = 0.e0_wp
        erl    = 0.e0_wp
        erle   = 0.e0_wp
        erln   = 0.e0_wp
        erlp   = 0.e0_wp
        erlr   = 0.e0_wp
        ers    = 0.e0_wp
        erse   = 0.e0_wp
        ersn   = 0.e0_wp
        ersp   = 0.e0_wp
        ersr   = 0.e0_wp
        f1     = 0.e0_wp
        f3     = 0.e0_wp
        ff     = 0.e0_wp
        fg     = 0.e0_wp
        fh     = 0.e0_wp
        hr     = 0.e0_wp
        rmx    = 0.e0_wp
        rmxe   = 0.e0_wp
        rmxn   = 0.e0_wp
        rmxp   = 0.e0_wp
        rmxr   = 0.e0_wp
        rnah16 = 0.e0_wp
        rnah36 = 0.e0_wp
        eps    = 1.e-4_wp
        ampn   = 2.35737419e-1_wp
        ampp   = 1.51495000e+0_wp
        ampe   = 5.02447648e-1_wp
        ampr   = 4.05991868e+0_wp
        erln   = 2.30065420e+2_wp
        erlp   = 7.10000000e-2_wp
        erle   = 3.46985415e+0_wp
        erlr   = 8.95373168e-1_wp
        ersn   = 2.00000000e+0_wp
        ersp   = 1.16830000e+1_wp
        erse   = 4.54595148e+1_wp
        ersr   = -2.38192975e+0_wp
        rmxn   = 2.33484458e+0_wp
        rmxp   = 2.62222000e+0_wp
        rmxe   = 2.04523701e+1_wp
        rmxr   = 1.98941949e+0_wp
        cevau  = 1.e0_wp/27.211611e0_wp
        write(6, '(/28(''=''),'' NaH2 - U12 coupling '',29(''=''), /''ampn='',1pe13.5,''  ampp='',e13.5,''  ampe='',e13.5,' // &
            '''  ampr='',e13.5 /''erln='',e13.5,''  erlp='',e13.5,''  erle='',e13.5,''  erlr='',e' // &
            '13.5, /''ersn='',e13.5,''  ersp='',e13.5,''  erse='',e13.5,''  ersr='',e13.5, /''rmx' // &
            'n='',e13.5,''  rmxp='',e13.5,''  rmxe='',e13.5,''  rmxr='',e13.5, /78(''=''))') ampn, ampp, ampe, ampr, erln, erlp, erle, erlr, &
            & ersn, ersp, erse, ersr, rmxn, rmxp, rmxe, rmxr
        do i=1, nt
            amp    = ampn+(ampp-ampn)*0.5e0_wp*(tanh(ampe*(r(i,2)-ampr))+1.e0_wp)
            erl    = erln+(erlp-erln)*0.5e0_wp*(tanh(erle*(r(i,2)-erlr))+1.e0_wp)
            ers    = ersn+(ersp-ersn)*0.5e0_wp*(tanh(erse*(r(i,2)-ersr))+1.e0_wp)
            rmx    = rmxn+(rmxp-rmxn)*0.5e0_wp*(tanh(rmxe*(r(i,2)-rmxr))+1.e0_wp)
            f1     = amp*exp(-(erl+ers/(eps+r(i,1)**3))*(r(i,1)-rmx)**2)
            f3     = amp*exp(-(erl+ers/(eps+r(i,3)**3))*(r(i,3)-rmx)**2)
            rnah16 = r(i,1)**6
            rnah36 = r(i,3)**6
            ff     = ((r(i,1)-r(i,3))/r(i,2))**2
            fg     = ((rnah16-rnah36)/(rnah16+rnah36))**2
            hr     = sqrt(r(i,1)**2+r(i,2)**2+r(i,3)**2)
            fh     = 0.5e0_wp*(1.e0_wp+tanh(0.75e0_wp*(11.e0_wp-hr)))
            ang    = (ff-fg)*fh+fg
            e(i)   = ang*(f1+f3)*cevau
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_5g_1992_dpem_pot12_legacy_state_dd307a13

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, nt )
    use chempotpy_nah2_leps_5g_1992_dpem_pot12_legacy_state_dd307a13, only: dispatch
    use chempotpy_kinds,                                              only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot12

!  U22 surface - upper diabatic surface
!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_nah2_leps_5g_1992_dpem_pot22_legacy_state_dd307a13
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
        real(wp) :: a0g1, a0g2, a0g3, a90g1, a90g2, a90g3, al1, al2, alb, alhh, alt, b0g1, b0g2, b0g3, b90g1, b90g2, &
            & b90g3, c1s1, c2a, c2b, c2c, c2s1, c3s1, c4s1, c5s1, c6, c6s1, cevau, cfm, cof, coul1, coul2, coul3, cplg2, cprime, &
            & cstsq, deh2, e(:), exch1, exch2, exch3, onethd, r(:,:), r0g1, r0g2, r0g3, r1s, r2s, r3s, r90g1, r90g2, r90g3, &
            & rac2, rat, rf11, rf12, rf13, rf21, rf22, rf23, rnag, s1, s2, s2a, s2b, s3, sntsq, swa, swr0, swt, t01, t02a, t03, &
            & t1, t2, t2a, t3, t901, t902a, t903, tr0, twothd, una2p, vdnahh, w, x0t1, x0t2, x0t3, x90t1, x90t2, x90t3, y0b1, &
            & y0b2, y0d1, y0d2, y0r1, y0r2, y0swa, y0swr0, y90b1, y90b2, y90d1, y90d2, y90r1, y90r2, y90swa, y90swr0
        integer :: i, nt
        a0g1    = 0.e0_wp
        a0g2    = 0.e0_wp
        a0g3    = 0.e0_wp
        a90g1   = 0.e0_wp
        a90g2   = 0.e0_wp
        a90g3   = 0.e0_wp
        al1     = 0.e0_wp
        al2     = 0.e0_wp
        alb     = 0.e0_wp
        alhh    = 0.e0_wp
        alt     = 0.e0_wp
        b0g1    = 0.e0_wp
        b0g2    = 0.e0_wp
        b0g3    = 0.e0_wp
        b90g1   = 0.e0_wp
        b90g2   = 0.e0_wp
        b90g3   = 0.e0_wp
        c1s1    = 0.e0_wp
        c2a     = 0.e0_wp
        c2b     = 0.e0_wp
        c2c     = 0.e0_wp
        c2s1    = 0.e0_wp
        c3s1    = 0.e0_wp
        c4s1    = 0.e0_wp
        c5s1    = 0.e0_wp
        c6      = 0.e0_wp
        c6s1    = 0.e0_wp
        cevau   = 0.e0_wp
        cfm     = 0.e0_wp
        cof     = 0.e0_wp
        coul1   = 0.e0_wp
        coul2   = 0.e0_wp
        coul3   = 0.e0_wp
        cplg2   = 0.e0_wp
        cprime  = 0.e0_wp
        cstsq   = 0.e0_wp
        deh2    = 0.e0_wp
        exch1   = 0.e0_wp
        exch2   = 0.e0_wp
        exch3   = 0.e0_wp
        onethd  = 0.e0_wp
        r0g1    = 0.e0_wp
        r0g2    = 0.e0_wp
        r0g3    = 0.e0_wp
        r1s     = 0.e0_wp
        r2s     = 0.e0_wp
        r3s     = 0.e0_wp
        r90g1   = 0.e0_wp
        r90g2   = 0.e0_wp
        r90g3   = 0.e0_wp
        rac2    = 0.e0_wp
        rat     = 0.e0_wp
        rf11    = 0.e0_wp
        rf12    = 0.e0_wp
        rf13    = 0.e0_wp
        rf21    = 0.e0_wp
        rf22    = 0.e0_wp
        rf23    = 0.e0_wp
        rnag    = 0.e0_wp
        s1      = 0.e0_wp
        s2      = 0.e0_wp
        s2a     = 0.e0_wp
        s2b     = 0.e0_wp
        s3      = 0.e0_wp
        sntsq   = 0.e0_wp
        swa     = 0.e0_wp
        swr0    = 0.e0_wp
        swt     = 0.e0_wp
        t01     = 0.e0_wp
        t02a    = 0.e0_wp
        t03     = 0.e0_wp
        t1      = 0.e0_wp
        t2      = 0.e0_wp
        t2a     = 0.e0_wp
        t3      = 0.e0_wp
        t901    = 0.e0_wp
        t902a   = 0.e0_wp
        t903    = 0.e0_wp
        tr0     = 0.e0_wp
        twothd  = 0.e0_wp
        una2p   = 0.e0_wp
        vdnahh  = 0.e0_wp
        w       = 0.e0_wp
        x0t1    = 0.e0_wp
        x0t2    = 0.e0_wp
        x0t3    = 0.e0_wp
        x90t1   = 0.e0_wp
        x90t2   = 0.e0_wp
        x90t3   = 0.e0_wp
        y0b1    = 0.e0_wp
        y0b2    = 0.e0_wp
        y0d1    = 0.e0_wp
        y0d2    = 0.e0_wp
        y0r1    = 0.e0_wp
        y0r2    = 0.e0_wp
        y0swa   = 0.e0_wp
        y0swr0  = 0.e0_wp
        y90b1   = 0.e0_wp
        y90b2   = 0.e0_wp
        y90d1   = 0.e0_wp
        y90d2   = 0.e0_wp
        y90r1   = 0.e0_wp
        y90r2   = 0.e0_wp
        y90swa  = 0.e0_wp
        y90swr0 = 0.e0_wp
        c1s1    = 2.614452e+2_wp
        c2s1    = -3.155166e-4_wp
        c3s1    = -1.785467e+0_wp
        c4s1    = -3.548074e+0_wp
        c5s1    = 3.356766e-1_wp
        c6s1    = -7.621511e-1_wp
        deh2    = 4.7469e0_wp
        una2p   = 2.1037e0_wp
        c2a     = 0.8e0_wp
        c2b     = 0.8e0_wp
        c2c     = 0.15e0_wp
        cprime  = 1.0e-1_wp
        y0d1    = 2.417800e-1_wp
        y0r1    = 3.708055e0_wp
        y0b1    = 6.960127e-1_wp
        y90d1   = 1.372065e-2_wp
        y90r1   = 7.233552e0_wp
        y90b1   = 5.127734e-1_wp
        y0d2    = 8.870172e-1_wp
        y0r2    = 2.300265e0_wp
        y0b2    = 8.497374e-1_wp
        y90d2   = 1.106187e0_wp
        y90r2   = 2.002516e0_wp
        y90b2   = 9.514887e-1_wp
        y0swa   = 0.614393e0_wp
        y0swr0  = 9.006323e0_wp
        y90swa  = 0.872116e0_wp
        y90swr0 = 8.194467e0_wp
        a0g1    = -8.210245e-2_wp
        b0g1    = 4.658698e0_wp
        r0g1    = 4.552860e0_wp
        a0g2    = 3.546427e-2_wp
        b0g2    = 6.456787e-1_wp
        r0g2    = 6.332225e0_wp
        a0g3    = -1.562327e-2_wp
        b0g3    = 1.170537e-1_wp
        r0g3    = 1.097310e+1_wp
        a90g1   = -9.929297e-2_wp
        b90g1   = 4.367220e0_wp
        r90g1   = 4.402447e0_wp
        a90g2   = 4.696181e-2_wp
        b90g2   = 5.470459e-1_wp
        r90g2   = 6.289526e0_wp
        a90g3   = -1.219900e-2_wp
        b90g3   = 1.002291e-2_wp
        r90g3   = 1.118725e+1_wp
        cfm     = 1.8427829108061610e-4_wp
        onethd  = 1.e0_wp/3.e0_wp
        twothd  = 2.e0_wp/3.e0_wp
        rac2    = 1.e0_wp/sqrt(2.e0_wp)
        cevau   = 1.e0_wp/27.2113961e0_wp
        write(6, '(/24(''=''),'' NaH2 - U22 potential '',25(''=''), /''  y90d1 ='',1pe14.6,''    y90r1 ='',e14.6,''  y90b1' // &
            ' ='',e14.6, /''  y90d2 ='',e14.6,  ''    y90r2 ='',e14.6,''  y90b2 ='',e14.6, /'' y9' // &
            '0swa ='',e14.6,  ''  y90swr0 ='',e14.6, /''   y0d1 ='',1pe14.6,''     y0r1 ='',e14.6' // &
            ',''   y0b1 ='',e14.6, /''   y0d2 ='',e14.6,  ''     y0r2 ='',e14.6,''   y0b2 ='',e14' // &
            '.6, /''  y0swa ='',e14.6,  ''   y0swr0 ='',e14.6, /''   c1s1 ='',1pe14.6,''     c2s1' // &
            ' ='',e14.6,''   c3s1 ='',e14.6, /''   c4s1 ='',1pe14.6,''     c5s1 ='',e14.6,''   c6' // &
            's1 ='',e14.6, /''    c2a ='',1pe14.6,''      c2b ='',e14.6,''    c2c ='',e14.6, /'' ' // &
            '  a0g1 ='',1pe14.6,''     b0g1 ='',e14.6,''   r0g1 ='',e14.6, /''   a0g2 ='',1pe14.6' // &
            ',''     b0g2 ='',e14.6,''   r0g2 ='',e14.6, /''   a0g3 ='',1pe14.6,''     b0g3 ='',e' // &
            '14.6,''   r0g3 ='',e14.6, /''  a90g1 ='',1pe14.6,''    b90g1 ='',e14.6,''  r90g1 =''' // &
            ',e14.6, /''  a90g2 ='',1pe14.6,''    b90g2 ='',e14.6,''  r90g2 ='',e14.6, /''  a90g3' // &
            ' ='',1pe14.6,''    b90g3 ='',e14.6,''  r90g3 ='',e14.6, /71(''=''))') y90d1, y90r1, y90b1, y90d2, y90r2, y90b2, y90swa, y90swr0, &
            & y0d1, y0r1, y0b1, y0d2, y0r2, y0b2, y0swa, y0swr0, c1s1, c2s1, c3s1, c4s1, &
            & c5s1, c6s1, c2a, c2b, c2c, a0g1, b0g1, r0g1, a0g2, b0g2, r0g2, a0g3, b0g3, &
            & r0g3, a90g1, b90g1, r90g1, a90g2, b90g2, r90g2, a90g3, b90g3, r90g3
        do i=1, nt
            r1s  = r(i,1)*r(i,1)
            r2s  = r(i,2)*r(i,2)
            r3s  = r(i,3)*r(i,3)
            rnag = sqrt(0.5e0_wp*abs(r1s+r3s-0.5e0_wp*r2s))
            if( rnag .lt. 1.e-10_wp ) stop 'RNAG < 1.D-10'
            cstsq = (0.5e0_wp*(r1s-r3s)/(r(i,2)*rnag))**2
            swa   = y90swa+(y0swa-y90swa)*cstsq
            swr0  = y90swr0+(y0swr0-y90swr0)*cstsq
            s1    = (c1s1+c2s1*r(i,1))*exp(c3s1*r(i,1)) +r1s*(c4s1+c5s1*r(i,1))*exp(c6s1*r(i,1))
            x90t2 = exp(-y90b2*(r(i,2)-y90r2))
            t902a = y90d2*x90t2*(x90t2-2.e0_wp)
            x0t2  = exp(-y0b2*(r(i,2)-y0r2))
            t02a  = y0d2*x0t2*(x0t2-2.e0_wp)
            t2a   = t902a+(t02a-t902a)*cstsq
            swt   = 0.5e0_wp*(1.e0_wp-tanh( swa*(0.5e0_wp*(r(i,1)+r(i,3))-swr0)))
            s2a=(149.480e0_wp-59.6557e0_wp*r(i,2))*exp(-4.13792e0_wp*r(i,2))+r2s *(-23.7299e0_wp+3.91747e0_wp*r(i,2))* &
                & exp(- 1.41350e0_wp* r(i,2))+una2p
            s2b=144.893e0_wp*exp(-3.85716e0_wp*r(i,2))+r(i,2)*(37.5919e0_wp -4.32985e0_wp*r(i,2)-0.003807e0_wp*r2s)* &
                & exp(- 1.52496e0_wp* r(i,2))
            t2    = swt*(t2a-s2b)+s2b
            s2    = 0.5e0_wp*(s2a+t2-tanh((s2a-t2)/cprime)*(s2a-t2))
            s3    = (c1s1+c2s1*r(i,3))*exp(c3s1*r(i,3)) +r3s*(c4s1+c5s1*r(i,3))*exp(c6s1*r(i,3))
            x90t1 = exp(-y90b1*(r(i,1)-y90r1))
            t901  = y90d1*x90t1*(x90t1+2.e0_wp)
            x0t1  = exp(-y0b1*(r(i,1)-y0r1))
            t01   = y0d1*x0t1*(x0t1+2.e0_wp)
            t1    = t901+(t01-t901)*cstsq
            x90t3 = exp(-y90b1*(r(i,3)-y90r1))
            t903  = y90d1*x90t3*(x90t3+2.e0_wp)
            x0t3  = exp(-y0b1*(r(i,3)-y0r1))
            t03   = y0d1*x0t3*(x0t3+2.e0_wp)
            t3    = t903+(t03-t903)*cstsq
            coul1 = 0.5e0_wp*(s1+t1)
            coul2 = 0.5e0_wp*(s2+t2)
            coul3 = 0.5e0_wp*(s3+t3)
            exch1 = 0.5e0_wp*(s1-t1)
            exch2 = 0.5e0_wp*(s2-t2)
            exch3 = 0.5e0_wp*(s3-t3)
            w     = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
            cplg2 = c2a*exp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
            e(i)  = coul1+coul2+coul3+deh2-rac2*sqrt(w+cplg2*cplg2)
            rf11  = 2.59329096e+0_wp*(1.e0_wp+tanh(9.93116245e+0_wp*(r(i,2) +6.85093745e+0_wp)))
            rf12  = 1.00540010e-1_wp*(1.e0_wp+tanh(-3.09412207e+0_wp*(r(i,2) -1.74921175e+0_wp)))
            rf13  = 6.80366735e-2_wp*(1.e0_wp+tanh(4.47330630e+0_wp*(r(i,2) -1.70178393e+0_wp)))
            tr0   = rf11+cstsq*(rf12+rf13*cstsq)
            sntsq = 1.e0_wp-cstsq
            rf21  = 4.52729085e+0_wp*(1.e0_wp+tanh(-7.15957152e-1_wp*(r(i,2) +5.34701010e-1_wp)))
            rf22  = 9.41855078e-2_wp*(1.e0_wp+tanh(-7.86831467e+0_wp*(r(i,2) -1.90926538e+0_wp)))
            rf23  = 2.84176599e-2_wp*(1.e0_wp+tanh(-6.59302513e+0_wp*(r(i,2) -6.09626510e+0_wp)))
            rat   = rf21+sntsq*(rf22+rf23*sntsq)
            al1=9.e0_wp+exp(0.238029900e+1_wp*(0.147400334e+1_wp-r(i,2))) *(-0.226331771e+0_wp-0.217339562e+1_wp*r(i, &
                & 2)+ 0.954597998e+1_wp * r(i,2)**2-0.233366749e+2_wp*r(i,2)**3+0.247115748e+2_wp*r(i,2)**4 - &
                & 0.128786102e+2_wp* r(i,2)**5+ 0.275711758e+1_wp* r(i,2)**6)
            al2=9.e0_wp+exp(0.135164094e+1_wp*(-0.254404564e+0_wp-r(i,2))) *(-0.106627055e+2_wp-0.140112406e+2_wp*r(i, &
                & 2)- 0.780660398e+1_wp * r(i,2)**2+0.547697753e+1_wp*r(i,2)**3-0.695875664e+1_wp*r(i,2)**4 + &
                & 0.364334771e+1_wp* r(i,2)**5- 0.492505596e+0_wp*r(i,2)**6)
            alt    = twothd*(al1-al2)
            alb    = onethd*al1+twothd*al2
            alhh   = alb+alt*(1.5e0_wp*cstsq-0.5e0_wp)
            c6     = 505.5e0_wp*alhh/(18.35755920e0_wp+sqrt(alb*0.5e0_wp))
            vdnahh = (-c6*rnag**(-6))*(1.e0_wp-exp((-rnag/tr0)**7))
            cof    = 0.5e0_wp*(1.e0_wp+tanh(rnag-r(i,2)))
            e(i)   = e(i)*cevau+vdnahh*cof*rat-cfm
        enddo
    end subroutine dispatch
end module chempotpy_nah2_leps_5g_1992_dpem_pot22_legacy_state_dd307a13

!> @brief Compatibility interface for pot22.
subroutine pot22( r, e, nt )
    use chempotpy_nah2_leps_5g_1992_dpem_pot22_legacy_state_dd307a13, only: dispatch
    use chempotpy_kinds,                                              only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot22
