!> @brief Implement NaFH_LEPS_TTYPP1_1998, a potential energy surface for NaFH.
module chempotpy_nafh_leps_ttypp1_1998_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NaFH_LEPS_TTYPP1_1998 potential energy surface for NaFH.
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
        real(wp)              :: nt, r(1,3), v(1)
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
end module chempotpy_nafh_leps_ttypp1_1998_pes

!> @brief Expose NaFH_LEPS_TTYPP1_1998 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_nafh_leps_ttypp1_1998_pes, only: evaluate
    use chempotpy_kinds,                     only: wp
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

!> @brief Evaluate diagonalize for the NaFH_LEPS_TTYPP1_1998 potential.
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

!> @brief Evaluate pot for the NaFH_LEPS_TTYPP1_1998 potential.
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
module chempotpy_nafh_leps_ttypp1_1998_pot11_legacy_state_b37154a3
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
        real(wp) :: a, abs1, ads1, ars1, as1, as3, axs1, axs3, betat2, c1, c2, c2a, c2b, c2c, c3, c4, cbs1, cds1, cevau, &
            & cfm, coul1, coul2, coul3, cplg2, cprime, crs1, cs1, cs3, cxs1, cxs3, deh2, dn, dt2, e(:), exch1, exch2, exch3, &
            & hlp, hlp1, hlp2, hlp3, r(:,:), r2s, rac2, ret2, rootw, rr1, rr1s, rr3, rr3s, s1, s2, s2a, s2h, s2i, s3, scale, &
            & swa, swa2, swr0, swr20, swt, t1, t2, t3, w, xt2
        integer :: i, nt
        a      = 0.e0_wp
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
        cprime = 0.e0_wp
        crs1   = 0.e0_wp
        cs1    = 0.e0_wp
        cs3    = 0.e0_wp
        cxs1   = 0.e0_wp
        cxs3   = 0.e0_wp
        deh2   = 0.e0_wp
        dn     = 0.e0_wp
        dt2    = 0.e0_wp
        exch1  = 0.e0_wp
        exch2  = 0.e0_wp
        exch3  = 0.e0_wp
        hlp    = 0.e0_wp
        hlp1   = 0.e0_wp
        hlp2   = 0.e0_wp
        hlp3   = 0.e0_wp
        r2s    = 0.e0_wp
        rac2   = 0.e0_wp
        ret2   = 0.e0_wp
        rootw  = 0.e0_wp
        rr1    = 0.e0_wp
        rr1s   = 0.e0_wp
        rr3    = 0.e0_wp
        rr3s   = 0.e0_wp
        s1     = 0.e0_wp
        s2     = 0.e0_wp
        s2a    = 0.e0_wp
        s2h    = 0.e0_wp
        s2i    = 0.e0_wp
        s3     = 0.e0_wp
        scale  = 0.e0_wp
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
        scale  = 1.5e0_wp
        cprime = 4.0e-1_wp
        dt2    = 2.885656e-2_wp
        ret2   = 1.237624e0_wp
        betat2 = 1.038036e+1_wp
        c1     = 1255969.4e0_wp
        c2     = 4.7e0_wp
        c3     = 5.546257e0_wp
        c4     = 5.573967e-1_wp
        ads1   = 0.08e0_wp
        ars1   = 4.0e0_wp
        abs1   = 1.55e0_wp
        cds1   = 2.40e0_wp
        crs1   = 4.15e0_wp
        cbs1   = 1.20e0_wp
        swa    = 0.7e0_wp
        swr0   = 6.5e0_wp
        swa2   = 1.1e0_wp
        swr20  = 1.9e0_wp
        c2a    = 1.5e0_wp
        c2b    = 1.0e0_wp
        c2c    = 0.15e0_wp
        deh2   = 3.39687e0_wp
        cfm    = 0.0e0_wp
        rac2   = 1.e0_wp/sqrt(2.e0_wp)
        cevau  = 1.e0_wp/27.2113961e0_wp
        do i=1, nt
            r2s  = r(i,2)*r(i,2)
            rr1  = r(i,1)*scale
            rr3  = r(i,3)*scale
            rr1s = rr1**2
            rr3s = rr3**2
            !*****SWITCHING FUNCTION FOR THE SINGLETS BELOW, AND ITS DERIVATIVES:
            hlp1 = 1.e0_wp-tanh(swa*(rr1-swr0))
            hlp2 = 1.e0_wp+tanh(swa2*(r(i,2)-swr20))
            hlp3 = 1.e0_wp-tanh(swa*(rr3-swr0))
            swt  = 0.125e0_wp*hlp1*hlp2*hlp3
            !*****SINGLETS AND THEIR DERIVATIVES:
            !*****MH SINGLET:
            axs1 = exp(-0.5e0_wp*abs1*(rr1-ars1))
            as1  = 0.5e0_wp*ads1*axs1*(axs1-2.e0_wp)
            cxs1 = exp(-2.80e0_wp*abs1*(rr1-3.500e0_wp))
            cs1  = 0.35e0_wp*cxs1*(cxs1+2.e0_wp)
            s1   = as1+(cs1-as1)*swt
            !*****H2 SINGLET, AND ITS DERIVATIVES:
            !     ORIGINAL
            s2a=(109.384e0_wp+463.033e0_wp*r(i,2))*exp(-6.240e0_wp*r(i,2)) -r2s*(13.654e0_wp+20.717e0_wp*r(i,2))*exp(- &
                & 2.032e0_wp* r(i,2))+ 1.34357e0_wp
            hlp1 = exp(-8.e0_wp*(r(i,2)-1.411655e0_wp))
            hlp2 = 0.2e0_wp*exp(-2.032e0_wp*r(i,2))
            s2h  = hlp1 - hlp2
            hlp  = tanh((s2a-s2h)/cprime)
            s2a  = 0.5e0_wp*(s2a+s2h-hlp*(s2a-s2h))
            a    = 1.00075e0_wp
            dn   = 3.39687e0_wp
            hlp  = exp(-a*(r(i,2)-1.42e0_wp))
            s2i  = dn*hlp*(hlp-2.e0_wp)
            hlp  = 0.5e0_wp*(r(i,1)+r(i,3)-8.e0_wp)
            hlp  = 0.5e0_wp*(1.e0_wp-tanh(hlp))
            s2   = s2a+(s2i-s2a)*hlp
            !*****MF SINGLET:
            axs3 = exp(-1.3e0_wp*abs1*(rr3-ars1))
            as3  = ads1*axs3*(axs3-2.e0_wp)
            cxs3 = exp(-cbs1*(rr3-crs1))
            cs3  = cds1*cxs3*(cxs3-2.e0_wp)
            s3   = as3+(cs3-as3)*swt
            !*****TRIPLETS AND THEIR FIRST DERIVATIVES:
            !     MH TRIPLET
            hlp1 = exp(-1.10e0_wp*c2*(rr1-0.6e0_wp))
            hlp2 = exp(-1.65e0_wp*c4*(rr1-0.6e0_wp))
            t1   = c1*hlp1+c3*hlp2
            !     HH TRIPLET
            xt2 = exp(-betat2*(r(i,2)-ret2))
            t2  = dt2*xt2*(xt2-2.e0_wp)
            !     MF TRIPLET
            hlp1 = exp(-c2*(rr3-0.55e0_wp))
            hlp2 = exp(-2.5e0_wp*c4*(rr3-0.55e0_wp))
            t3   = c1*hlp1+c3*hlp2
            !*****COULOMB INTEGRALS AND THEIR FIRST DERIVATIVES:
            coul1 = 0.5e0_wp*(s1+t1)
            coul2 = 0.5e0_wp*(s2+t2)
            coul3 = 0.5e0_wp*(s3+t3)
            !*****EXCHANGE INTEGRALS AND THEIR FIRST DERIVATIVES:
            exch1 = 0.5e0_wp*(s1-t1)
            exch2 = 0.5e0_wp*(s2-t2)
            exch3 = 0.5e0_wp*(s3-t3)
            w     = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
            cplg2 = c2a*exp(-c2b*w-c2c*(rr1+r(i,2)+rr3))
            rootw = rac2*sqrt(w+cplg2*cplg2)
            e(i)  = ((coul1+coul2+coul3+deh2-rootw)*cevau-cfm)
        enddo
    end subroutine dispatch
end module chempotpy_nafh_leps_ttypp1_1998_pot11_legacy_state_b37154a3

!> @brief Compatibility interface for pot11.
subroutine pot11( r, e, nt )
    use chempotpy_nafh_leps_ttypp1_1998_pot11_legacy_state_b37154a3, only: dispatch
    use chempotpy_kinds,                                             only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot11

!  U12 surface - coupling surface
!> @brief Own the stateless g_pot12 evaluator.
module chempotpy_nafh_leps_ttypp1_1998_pot12_legacy_state_b37154a3
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
        real(wp) :: amc, cevau, e(:), gmc, hlp, r(:,:), rmc, sw1, sw21, sw23, unaf, unah
        integer  :: i, nt
        !      BLOCK DATA UONETWO
        real(wp) :: ahf, anaf, anah, r21, r23, r13
        real(wp) :: ccmf(4), ccmh(4)
        !      END
        !      save /ucoupl/,/ucoupl1/
        ahf   = 0.e0_wp
        amc   = 0.e0_wp
        anaf  = 0.e0_wp
        anah  = 0.e0_wp
        ccmf  = 0.e0_wp
        ccmh  = 0.e0_wp
        cevau = 0.e0_wp
        gmc   = 0.e0_wp
        hlp   = 0.e0_wp
        r13   = 0.e0_wp
        r21   = 0.e0_wp
        r23   = 0.e0_wp
        rmc   = 0.e0_wp
        sw1   = 0.e0_wp
        sw21  = 0.e0_wp
        sw23  = 0.e0_wp
        unaf  = 0.e0_wp
        unah  = 0.e0_wp
        ccmf  = [ 3.9960e0_wp, 5.4730e-01_wp, -1.1428e0_wp, 2.0440e-01_wp]
        ccmh  = [ 1.00e0_wp, 0.80e0_wp, -2.67e0_wp, 0.456e0_wp]
        ahf   = 3.0000e0_wp
        anaf  = 18.0000e-1_wp
        anah  = 16.0000e-1_wp
        r21   = 0.7000e0_wp
        r23   = -3.3006e-1_wp
        r13   = 1.2431e0_wp
        cevau = 1.e0_wp/27.2113961e0_wp
        do i=1, nt
            hlp  = ccmf(3)/(1.e0_wp+ccmf(4)*r(i,3))
            unaf = exp(-ccmf(1)-ccmf(2)*r(i,3)- hlp*r(i,3))
            unaf = unaf*27.2113961e0_wp
            hlp  = ccmh(3)/(1.e0_wp+ccmh(4)*r(i,1))
            unah = exp(-ccmh(1)-ccmh(2)*r(i,1)- hlp*r(i,1))
            hlp  = anah*r(i,1)-anaf*r(i,3)-r13
            sw1  = 0.5e0_wp*(1.e0_wp+tanh(hlp))
            hlp  = ahf*r(i,2)-anah*r(i,1)-r21
            sw21 = 0.5e0_wp*(1.e0_wp+tanh(hlp))
            hlp  = ahf*r(i,2)-anaf*r(i,3)-r23
            sw23 = 0.5e0_wp*(1.e0_wp+tanh(hlp))
            e(i) = unaf*sw1*sw23+unah*(1.e0_wp-sw1)*sw21
            e(i) = e(i)*cevau
            !     MODIFICATION: coupling cut-off at large H-H separations
            gmc  = 0.3e0_wp
            amc  = 0.40e0_wp
            rmc  = 2.7e0_wp
            hlp  = gmc*(1.e0_wp-tanh(amc*(r(i,2)-rmc)))
            e(i) = e(i)*hlp
        enddo
    end subroutine dispatch
end module chempotpy_nafh_leps_ttypp1_1998_pot12_legacy_state_b37154a3

!> @brief Compatibility interface for pot12.
subroutine pot12( r, e, nt )
    use chempotpy_nafh_leps_ttypp1_1998_pot12_legacy_state_b37154a3, only: dispatch
    use chempotpy_kinds,                                             only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot12

!  U22 surface - upper diabatic surface
!> @brief Own the stateless g_pot22 evaluator.
module chempotpy_nafh_leps_ttypp1_1998_pot22_legacy_state_b37154a3
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
        real(wp) :: a, a0g1, a0g2, a0g3, a90g1, a90g2, a90g3, aaa, aaa1, af1, b, b0g1, b0g2, b0g3, b90g1, b90g2, b90g3, &
            & c1s1, c2a, c2b, c2c, c2s1, c3s1, c4s1, c5s1, c6s1, cevau, coul1, coul2, coul3, cplg2, cprime, daa, deh2, dn, &
            & dnaht1, dnaht11, dnaht2, dnaht3, dnaht4, dnaht5, dnaht6, dnaht61, dnaht7, e(:), exch1, exch2, exch3, fff, h2t1, &
            & h2t2, h2t3, h2t4, h2t5, h2t6, h2t7, hlp, hlp1, hlp2, hlpnah1, hlpnah11, hlpnah3, hlpnah33, r(:,:), r00, r0g1, &
            & r0g2, r0g3, r1s, r2s, r3s, r90g1, r90g2, r90g3, raaa, rac2, rf1, rootw, rr1, rr1s, rr3, rr3s, s1, s2, s2a, s2b, &
            & s2i, s3, scale, swa, swr0, swt, t1, t2, t2a, t3, una2p, w, y0b1, y0b2, y0d1, y0d2, y0r1, y0r2, y90b1, y90b2, &
            & y90d1, y90d2, y90r1, y90r2
        integer :: i, nt
        a        = 0.e0_wp
        a0g1     = 0.e0_wp
        a0g2     = 0.e0_wp
        a0g3     = 0.e0_wp
        a90g1    = 0.e0_wp
        a90g2    = 0.e0_wp
        a90g3    = 0.e0_wp
        aaa      = 0.e0_wp
        aaa1     = 0.e0_wp
        af1      = 0.e0_wp
        b        = 0.e0_wp
        b0g1     = 0.e0_wp
        b0g2     = 0.e0_wp
        b0g3     = 0.e0_wp
        b90g1    = 0.e0_wp
        b90g2    = 0.e0_wp
        b90g3    = 0.e0_wp
        c1s1     = 0.e0_wp
        c2a      = 0.e0_wp
        c2b      = 0.e0_wp
        c2c      = 0.e0_wp
        c2s1     = 0.e0_wp
        c3s1     = 0.e0_wp
        c4s1     = 0.e0_wp
        c5s1     = 0.e0_wp
        c6s1     = 0.e0_wp
        cevau    = 0.e0_wp
        coul1    = 0.e0_wp
        coul2    = 0.e0_wp
        coul3    = 0.e0_wp
        cplg2    = 0.e0_wp
        cprime   = 0.e0_wp
        daa      = 0.e0_wp
        deh2     = 0.e0_wp
        dn       = 0.e0_wp
        dnaht1   = 0.e0_wp
        dnaht11  = 0.e0_wp
        dnaht2   = 0.e0_wp
        dnaht3   = 0.e0_wp
        dnaht4   = 0.e0_wp
        dnaht5   = 0.e0_wp
        dnaht6   = 0.e0_wp
        dnaht61  = 0.e0_wp
        dnaht7   = 0.e0_wp
        exch1    = 0.e0_wp
        exch2    = 0.e0_wp
        exch3    = 0.e0_wp
        fff      = 0.e0_wp
        h2t1     = 0.e0_wp
        h2t2     = 0.e0_wp
        h2t3     = 0.e0_wp
        h2t4     = 0.e0_wp
        h2t5     = 0.e0_wp
        h2t6     = 0.e0_wp
        h2t7     = 0.e0_wp
        hlp      = 0.e0_wp
        hlp1     = 0.e0_wp
        hlp2     = 0.e0_wp
        hlpnah1  = 0.e0_wp
        hlpnah11 = 0.e0_wp
        hlpnah3  = 0.e0_wp
        hlpnah33 = 0.e0_wp
        r00      = 0.e0_wp
        r0g1     = 0.e0_wp
        r0g2     = 0.e0_wp
        r0g3     = 0.e0_wp
        r1s      = 0.e0_wp
        r2s      = 0.e0_wp
        r3s      = 0.e0_wp
        r90g1    = 0.e0_wp
        r90g2    = 0.e0_wp
        r90g3    = 0.e0_wp
        raaa     = 0.e0_wp
        rac2     = 0.e0_wp
        rf1      = 0.e0_wp
        rootw    = 0.e0_wp
        rr1      = 0.e0_wp
        rr1s     = 0.e0_wp
        rr3      = 0.e0_wp
        rr3s     = 0.e0_wp
        s1       = 0.e0_wp
        s2       = 0.e0_wp
        s2a      = 0.e0_wp
        s2b      = 0.e0_wp
        s2i      = 0.e0_wp
        s3       = 0.e0_wp
        scale    = 0.e0_wp
        swa      = 0.e0_wp
        swr0     = 0.e0_wp
        swt      = 0.e0_wp
        t1       = 0.e0_wp
        t2       = 0.e0_wp
        t2a      = 0.e0_wp
        t3       = 0.e0_wp
        una2p    = 0.e0_wp
        w        = 0.e0_wp
        y0b1     = 0.e0_wp
        y0b2     = 0.e0_wp
        y0d1     = 0.e0_wp
        y0d2     = 0.e0_wp
        y0r1     = 0.e0_wp
        y0r2     = 0.e0_wp
        y90b1    = 0.e0_wp
        y90b2    = 0.e0_wp
        y90d1    = 0.e0_wp
        y90d2    = 0.e0_wp
        y90r1    = 0.e0_wp
        y90r2    = 0.e0_wp
        scale    = 1.5e0_wp
        c1s1     = 203.128e0_wp
        c2s1     = 1.6167e0_wp
        c3s1     = -100.124e0_wp
        c4s1     = 68.078e0_wp
        c5s1     = -15.832e0_wp
        c6s1     = 1.303e0_wp
        dnaht1   = 3.35000e0_wp
        dnaht2   = 10.000e0_wp
        dnaht3   = 54.6393e0_wp
        dnaht4   = -17.63343e0_wp
        dnaht5   = 3.0085e0_wp
        dnaht6   = 1.450e0_wp
        dnaht7   = -2.87977e0_wp
        dnaht61  = 1.300e0_wp
        dnaht11  = 2.85000e0_wp
        deh2     = 4.7404422e0_wp
        una2p    = 2.1037e0_wp
        c2a      = 0.8e0_wp
        c2b      = 0.8e0_wp
        c2c      = 0.15e0_wp
        cprime   = 8.0e-1_wp
        y0d1     = 2.417800e-1_wp
        y0r1     = 3.708055e0_wp
        y0b1     = 6.960127e-1_wp
        y90d1    = 1.372065e-2_wp
        y90r1    = 7.233552e0_wp
        y90b1    = 5.127734e-1_wp
        y0d2     = 8.870172e-1_wp
        y0r2     = 2.300265e0_wp
        y0b2     = 8.497374e-1_wp
        y90d2    = 1.106187e0_wp
        y90r2    = 2.002516e0_wp
        y90b2    = 9.514887e-1_wp
        swa      = 0.78495e0_wp
        swr0     = 7.47654e0_wp
        a0g1     = -8.210245e-2_wp
        b0g1     = 4.658698e0_wp
        r0g1     = 4.552860e0_wp
        a0g2     = 3.546427e-2_wp
        b0g2     = 6.456787e-1_wp
        r0g2     = 6.332225e0_wp
        a0g3     = -1.562327e-2_wp
        b0g3     = 1.170537e-1_wp
        r0g3     = 1.097310e+1_wp
        a90g1    = -9.929297e-2_wp
        b90g1    = 4.367220e0_wp
        r90g1    = 4.402447e0_wp
        a90g2    = 4.696181e-2_wp
        b90g2    = 5.470459e-1_wp
        r90g2    = 6.289526e0_wp
        a90g3    = -1.219900e-2_wp
        b90g3    = 1.002291e-2_wp
        r90g3    = 1.118725e+1_wp
        h2t1     = 1.5e0_wp
        h2t2     = 4.47009e0_wp
        h2t3     = 4.48534e0_wp
        h2t4     = 1.77908e0_wp
        h2t5     = -6.46276e0_wp
        h2t6     = 3.15396e0_wp
        h2t7     = -1.59912e0_wp
        rac2     = 1.e0_wp/sqrt(2.e0_wp)
        cevau    = 1.e0_wp/27.2113961e0_wp
        do i=1, nt
            r1s  = r(i,1)*r(i,1)
            r2s  = r(i,2)*r(i,2)
            r3s  = r(i,3)*r(i,3)
            rr1  = r(i,1)*scale
            rr3  = r(i,3)*scale
            rr1s = rr1**2
            rr3s = rr3**2
            !*****SWITCHING FUNCTION FOR THE H2 TRIPLET, AND ITS DERIVATIVES:
            swt = 0.5e0_wp*(1.e0_wp-tanh( swa*(0.5e0_wp*(rr1+rr3)-swr0)))
            !*****SINGLETS AND THEIR DERIVATIVES:
            !*****MH SINGLET:
            raaa     = 2.700e0_wp
            aaa      = 2.40e0_wp
            aaa1     = 2.60e0_wp
            daa      = 1.80e0_wp
            hlpnah1  = exp(-aaa1*(r(i,1)-raaa))
            hlpnah11 = exp(-1.65e0_wp*aaa*(r(i,1)-raaa))
            s1       = 0.30e0_wp*daa*(hlpnah11-2.e0_wp*hlpnah1)
            !*****H2 SINGLET
            s2a=(109.384e0_wp+463.033e0_wp*r(i,2))*exp(-6.240e0_wp*r(i,2)) -r2s*(13.654e0_wp+20.717e0_wp*r(i,2))*exp(- &
                & 2.032e0_wp* r(i,2))+ una2p
            b    = 1.4e0_wp
            a    = 2.e0_wp*b
            r00  = 1.411655e0_wp
            dn   = -4.74044e0_wp/(1.e0_wp-a/b)
            hlp1 = exp(-a*(r(i,2)-r00))
            hlp2 = exp(-b*(r(i,2)-r00))
            s2i  = dn*(hlp1-(a/b)*hlp2)+una2p
            hlp  = 0.5e0_wp*(r(i,1)+r(i,3)-7.0e0_wp)
            hlp  = 0.5e0_wp*(1.e0_wp-tanh(hlp))
            s2a  = s2a+(s2i-s2a)*hlp
            !*****MF SINGLET:
            hlpnah3  = exp(-aaa*(r(i,3)-raaa))
            hlpnah33 = exp(-1.65e0_wp*aaa*(r(i,3)-raaa))
            s3       = daa*(hlpnah33-2.e0_wp*hlpnah3)
            !*****TRIPLETS AND THEIR DERIVATIVES:
            !*****H2 TRIPLET (interaction region):
            t2a = 0.5e0_wp*(exp(-h2t2*(r(i,2)-h2t1))+(h2t3+h2t5*r(i,2) +h2t6*r2s+h2t7*r2s*r(i,2))*exp(-h2t4*r(i,2)))
            !*****H2 TRIPLET (asymptotic)
            s2b = 147.214e0_wp*exp(-3.915e0_wp*r(i,2))+r(i,2)*(41.275e0_wp+10.505e0_wp *r(i,2)+4.408e0_wp*r2s)*exp(- &
                & 2.032e0_wp* r(i,2))
            !*****H2 TRIPLET: TOTAL (switched between the asymptotic and the ineract
            t2 = swt*(t2a-s2b)+s2b
            !*****TRIPLETS AND THEIR DERIVATIVES:
            !!!!!!MH TRIPLET:
            t1 = 3.e0_wp*(exp(-dnaht2*(rr1-dnaht11))+ (dnaht7+dnaht3*rr1+dnaht4*rr1s+dnaht5*rr1*rr1s) *exp(- dnaht61* &
                & rr1))
            !!!!!!MF TRIPLET: (should be the same as the first one, t1)
            t3 = 3.e0_wp*(exp(-dnaht2*(rr3-dnaht1))+ (dnaht7+dnaht3*rr3+dnaht4*rr3s+dnaht5*rr3*rr3s) *exp(-dnaht6* rr3))
            !*****TOTAL H2 SINGLET:
            !!!!!!H2 SINGLET: total (switched between the true singlet and a triplet
            !                        in the asymptotic region at large HH distances)
            hlp = tanh((s2a-t2)/cprime)
            s2  = 0.5e0_wp*(s2a+t2-hlp*(s2a-t2))
            !*****COULOMB INTEGRALS AND THEIR FIRST DERIVATIVES:
            coul1 = 0.5e0_wp*(s1+t1)
            coul2 = 0.5e0_wp*(s2+t2)
            coul3 = 0.5e0_wp*(s3+t3)
            !*****EXCHANGE INTEGRALS AND THEIR FIRST DERIVATIVES:
            exch1 = 0.5e0_wp*(s1-t1)
            exch2 = 0.5e0_wp*(s2-t2)
            exch3 = 0.5e0_wp*(s3-t3)
            w     = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
            cplg2 = c2a*exp(-c2b*w-c2c*(rr1+r(i,2)+rr3))
            rootw = rac2*sqrt(w+cplg2*cplg2)
            e(i)  = coul1+coul2+coul3+deh2-rootw
            rf1   = 3.00e0_wp
            af1   = 0.60e0_wp
            fff   = 1.e0_wp - 0.5e0_wp*(1.e0_wp+tanh(af1*(r(i,3)-rf1)))
            e(i)  = ((e(i)-1.35e0_wp)/(1.00e0_wp+1.80e0_wp*fff))*cevau
        enddo
    end subroutine dispatch
end module chempotpy_nafh_leps_ttypp1_1998_pot22_legacy_state_b37154a3

!> @brief Compatibility interface for pot22.
subroutine pot22( r, e, nt )
    use chempotpy_nafh_leps_ttypp1_1998_pot22_legacy_state_b37154a3, only: dispatch
    use chempotpy_kinds,                                             only: wp
    implicit none
    real(wp) :: r(nt,3)
    real(wp) :: e(nt)
    integer  :: nt
    call dispatch(r=r, e=e, nt=nt)
end subroutine pot22
