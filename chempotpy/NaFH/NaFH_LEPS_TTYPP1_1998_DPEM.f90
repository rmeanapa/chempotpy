!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate dpem for the NaFH_LEPS_TTYPP1_1998_DPEM potential.
subroutine dpem( x, igrad, u, ug )
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
    logical, save                 :: first_time_data=.true.
    integer                       :: i, j, k, l
    !initialize
    u = 0.d0; ug = 0.d0

    nt = 1
    ! input cartesian is ClHH
    call legacy_triatomic_distances_to_bohr( x, r(1,:) )

    if(first_time_data) then
        call prepot
        first_time_data = .false.
    end if

    if (igrad==0) then
        call pot(r, u11, 1, 1)
        call pot(r, u12, 1, 2)
        call pot(r, u22, 1, 3)
        u(1,1) = u11(1)*legacy_hartree_to_electronvolt
        u(1,2) = u12(1)*legacy_hartree_to_electronvolt
        u(2,1) = u12(1)*legacy_hartree_to_electronvolt
        u(2,2) = u22(1)*legacy_hartree_to_electronvolt
    else
        write(*, *) "Only DPEM is available"
    end if

end subroutine dpem

!> @brief Evaluate diagonalize for the NaFH_LEPS_TTYPP1_1998_DPEM potential.
subroutine diagonalize( n, a_ss, ev_s, u_ss )
    implicit none
    integer, intent(in)           :: n
    double precision, intent(in)  :: a_ss(n,n)
    double precision, intent(out) :: u_ss(n,n)
    double precision, intent(out) :: ev_s(n)
    integer                       :: io, i
    double precision, allocatable :: work_d(:)
    integer, save                 :: lwork_d
    logical, save                 :: first_time_diag=.true.
    u_ss = a_ss
    if (first_time_diag) then
        lwork_d = -1
        allocate( work_d(2) )
        call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
        lwork_d = int(work_d(1))
        deallocate ( work_d )
        first_time_diag = .false.
    end if
    allocate( work_d(lwork_d) )
    call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
end subroutine diagonalize

!> @brief Evaluate prepot for the NaFH_LEPS_TTYPP1_1998_DPEM potential.
subroutine prepot
    !
    !   System:                        NaFH
    !   Common name:                   N3
    !   Functional form:               Modified extended LEPS
    !   Number of derivatives:         0
    !   Number of electronic surfaces: 2
    !   Interface:                     3-1V

    !   Reference:  Unpublished, M. S. Topaler and D. G. Truhlar

    !   Notes:  Development version of NaFH potential.
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
    !      ground and first excited electronic states, diabatic representati
    !
    !   Zero of energy:
    !      The classical potential energy is set equal to zero for the Na
    !      infinitely far from the HF diatomic and R(HF) set equal to the
    !      HF equilibrium diatomic value.
    !
    !   Parameters:
    !
    !   Coordinates:
    !      Internal, Definition: R(1) = R(Na-H)
    !                            R(2) = R(H-F)
    !                            R(3) = R(Na-F)

    implicit none
    call prepot11
    call prepot12
    call prepot22

end subroutine prepot

!
!> @brief Evaluate pot for the NaFH_LEPS_TTYPP1_1998_DPEM potential.
subroutine pot( r, e, nt, nsurf )
    implicit none
    double precision :: e, r
    integer          :: nsurf, nt
    dimension r(nt,3), e(nt)
    if (nsurf .eq. 1) then
        call pot11(r, e, nt)
    else if (nsurf .eq. 2) then
        call pot12(r, e, nt)
    else
        call pot22(r, e, nt)
    end if

end subroutine pot

!
!  U11 surface - lower diabatic surface
!
!> @brief Evaluate prepot11 for the NaFH_LEPS_TTYPP1_1998_DPEM potential.
subroutine prepot11
    implicit none
    double precision :: a, abs1, ads1, ars1, as1, as3, axs1, axs3, betat2, c1, c2, c2a, c2b, c2c, c3, &
        & c4, cbs1, cds1, cevau, cfm, coul1, coul2, coul3, cplg2, cprime, crs1, cs1, cs3, cxs1, cxs3, &
        & deh2, dn, dt2, e, exch1, exch2, exch3, hlp, hlp1, hlp2, hlp3, r, r2s, rac2, ret2, rootw, &
        & rr1, rr1s, rr3, rr3s, s1, s2, s2a, s2h, s2i, s3, scale, swa, swa2, swr0, swr20, swt, t1, &
        & t2, t3, w, xt2
    integer :: i, nt
    dimension r(nt,3), e(nt)
    save
    data scale/1.5d0/
    data cprime/4.0d-1/
    data dt2/2.885656d-2/, ret2/1.237624d0/, betat2/1.038036d+1/, &
        & c1/1255969.4d0/, c2/4.7d0/, c3/5.546257d0/, &
        & c4/5.573967d-1/
    data ads1/0.08d0/, ars1/4.0d0/, abs1/1.55d0/, &
        & cds1/2.40d0/, crs1/4.15d0/, cbs1/1.20d0/
    data swa/0.7d0/, swr0/6.5d0/, swa2/1.1d0/, swr20/1.9d0/, &
        & c2a/1.5d0/, c2b/1.0d0/, c2c/0.15d0/, deh2/3.39687d0/
    data cfm/0.0d0/
    rac2  = 1.d0/dsqrt(2.d0)
    cevau = 1.d0/27.2113961d0

    return
    !
    entry pot11(r, e, nt)

    do 10 i=1, nt

        !***********************************************************************
        r2s  = r(i,2)*r(i,2)
        rr1  = r(i,1)*scale
        rr3  = r(i,3)*scale
        rr1s = rr1**2
        rr3s = rr3**2
        !***********************************************************************

        !***********************************************************************
        !*****SWITCHING FUNCTION FOR THE SINGLETS BELOW, AND ITS DERIVATIVES:
        hlp1 = 1.d0-dtanh(swa*(rr1-swr0))
        hlp2 = 1.d0+dtanh(swa2*(r(i,2)-swr20))
        hlp3 = 1.d0-dtanh(swa*(rr3-swr0))
        swt  = 0.125d0*hlp1*hlp2*hlp3
        !***********************************************************************

        !***********************************************************************
        !*****SINGLETS AND THEIR DERIVATIVES:
        !-----------------------------------------------------------------------
        !*****MH SINGLET:
        axs1 = dexp(-0.5d0*abs1*(rr1-ars1))
        as1  = 0.5d0*ads1*axs1*(axs1-2.d0)

        cxs1 = dexp(-2.80d0*abs1*(rr1-3.500d0))
        cs1  = 0.35d0*cxs1*(cxs1+2.d0)

        s1 = as1+(cs1-as1)*swt
        !-----------------------------------------------------------------------
        !*****H2 SINGLET, AND ITS DERIVATIVES:
        !     ORIGINAL
        s2a=(109.384d0+463.033d0*r(i,2))*dexp(-6.240d0*r(i,2)) -r2s*(13.654d0+20.717d0*r(i,2))*dexp(-2.032d0*r(i,2))+ &
            & 1.34357d0

        hlp1 = dexp(-8.d0*(r(i,2)-1.411655d0))
        hlp2 = 0.2d0*dexp(-2.032d0*r(i,2))
        s2h  = hlp1 - hlp2

        hlp = dtanh((s2a-s2h)/cprime)
        !      hlp1=(s2a-s2h)/(cprime*dcosh((s2a-s2h)/cprime)**2)
        s2a = 0.5d0*(s2a+s2h-hlp*(s2a-s2h))

        a   = 1.00075d0
        dn  = 3.39687d0
        hlp = dexp(-a*(r(i,2)-1.42d0))
        s2i = dn*hlp*(hlp-2.d0)

        hlp = 0.5d0*(r(i,1)+r(i,3)-8.d0)
        hlp = 0.5d0*(1.d0-dtanh(hlp))
        s2  = s2a+(s2i-s2a)*hlp

        !-----------------------------------------------------------------------
        !*****MF SINGLET:
        axs3 = dexp(-1.3d0*abs1*(rr3-ars1))
        as3  = ads1*axs3*(axs3-2.d0)

        cxs3 = dexp(-cbs1*(rr3-crs1))
        cs3  = cds1*cxs3*(cxs3-2.d0)

        s3 = as3+(cs3-as3)*swt
        !***********************************************************************

        !***********************************************************************
        !*****TRIPLETS AND THEIR FIRST DERIVATIVES:
        !-----------------------------------------------------------------------
        !     MH TRIPLET
        hlp1 = dexp(-1.10d0*c2*(rr1-0.6d0))
        hlp2 = dexp(-1.65d0*c4*(rr1-0.6d0))
        t1   = c1*hlp1+c3*hlp2

        !-----------------------------------------------------------------------
        !     HH TRIPLET
        xt2 = dexp(-betat2*(r(i,2)-ret2))
        t2  = dt2*xt2*(xt2-2.d0)

        !-----------------------------------------------------------------------
        !     MF TRIPLET
        hlp1 = dexp(-c2*(rr3-0.55d0))
        hlp2 = dexp(-2.5d0*c4*(rr3-0.55d0))
        t3   = c1*hlp1+c3*hlp2
        !***********************************************************************

        !***********************************************************************
        !-----------------------------------------------------------------------
        !*****COULOMB INTEGRALS AND THEIR FIRST DERIVATIVES:

        coul1 = 0.5d0*(s1+t1)
        coul2 = 0.5d0*(s2+t2)
        coul3 = 0.5d0*(s3+t3)

        !-----------------------------------------------------------------------
        !*****EXCHANGE INTEGRALS AND THEIR FIRST DERIVATIVES:

        exch1 = 0.5d0*(s1-t1)
        exch2 = 0.5d0*(s2-t2)
        exch3 = 0.5d0*(s3-t3)

        !***********************************************************************

        !***********************************************************************
        w = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2

        cplg2 = c2a*dexp(-c2b*w-c2c*(rr1+r(i,2)+rr3))

        rootw = rac2*dsqrt(w+cplg2*cplg2)

        e(i) = ((coul1+coul2+coul3+deh2-rootw)*cevau-cfm)

        !***********************************************************************

    10 continue

end subroutine prepot11

!
!  U12 surface - coupling surface
!
!> @brief Evaluate prepot12 for the NaFH_LEPS_TTYPP1_1998_DPEM potential.
subroutine prepot12
    implicit none
    double precision :: amc, cevau, e, gmc, hlp, r, rmc, sw1, sw21, sw23, unaf, unah
    integer          :: i, nt
    dimension r(nt,3), e(nt)
    !      BLOCK DATA UONETWO
    double precision :: ahf, anaf, anah, r21, r23, r13
    double precision :: ccmf(4), ccmh(4)
    !      common/ucoupl/ahf,anaf,anah,r21,r23,r13
    !      common/ucoupl1/ccMF(4),ccMH(4)
    save
    data ccmf/3.9960d0, 5.4730d-01, -1.1428d0, 2.0440d-01/
    data ccmh/1.00d0, 0.80d0, -2.67d0, 0.456d0/
    data ahf/3.0000d0/, &
        & anaf/18.0000d-1/, &
        & anah/16.0000d-1/
    data r21/0.7000d0/, &
        & r23/-3.3006d-1/, &
        & r13/1.2431d0/
    !      END
    !      save /ucoupl/,/ucoupl1/
    !      common/ucoupl/ahf,anaf,anah,r21,r23,r13
    !      common/ucoupl1/ccMF(4),ccMH(4)
    cevau = 1.d0/27.2113961d0

    return
    !
    entry pot12(r, e, nt)

    do 10 i=1, nt

        hlp  = ccmf(3)/(1.d0+ccmf(4)*r(i,3))
        unaf = dexp(-ccmf(1)-ccmf(2)*r(i,3)- hlp*r(i,3))
        unaf = unaf*27.2113961d0

        hlp  = ccmh(3)/(1.d0+ccmh(4)*r(i,1))
        unah = dexp(-ccmh(1)-ccmh(2)*r(i,1)- hlp*r(i,1))
        !-----------------------------------------------------------------------

        hlp = anah*r(i,1)-anaf*r(i,3)-r13
        sw1 = 0.5d0*(1.d0+dtanh(hlp))

        hlp  = ahf*r(i,2)-anah*r(i,1)-r21
        sw21 = 0.5d0*(1.d0+dtanh(hlp))

        hlp  = ahf*r(i,2)-anaf*r(i,3)-r23
        sw23 = 0.5d0*(1.d0+dtanh(hlp))

        e(i) = unaf*sw1*sw23+unah*(1.d0-sw1)*sw21

        e(i) = e(i)*cevau

        !***********************************************************************
        !     MODIFICATION: coupling cut-off at large H-H separations
        gmc  = 0.3d0
        amc  = 0.40d0
        rmc  = 2.7d0
        hlp  = gmc*(1.d0-dtanh(amc*(r(i,2)-rmc)))
        e(i) = e(i)*hlp
        !***********************************************************************

    10 continue

end subroutine prepot12

!
!  U22 surface - upper diabatic surface
!
!> @brief Evaluate prepot22 for the NaFH_LEPS_TTYPP1_1998_DPEM potential.
subroutine prepot22
    implicit none
    double precision :: a, a0g1, a0g2, a0g3, a90g1, a90g2, a90g3, aaa, aaa1, af1, b, b0g1, b0g2, &
        & b0g3, b90g1, b90g2, b90g3, c1s1, c2a, c2b, c2c, c2s1, c3s1, c4s1, c5s1, c6s1, cevau, coul1, &
        & coul2, coul3, cplg2, cprime, daa, deh2, dn, dnaht1, dnaht11, dnaht2, dnaht3, dnaht4, &
        & dnaht5, dnaht6, dnaht61, dnaht7, e, exch1, exch2, exch3, fff, h2t1, h2t2, h2t3, h2t4, h2t5, &
        & h2t6, h2t7, hlp, hlp1, hlp2, hlpnah1, hlpnah11, hlpnah3, hlpnah33, r, r00, r0g1, r0g2, &
        & r0g3, r1s, r2s, r3s, r90g1, r90g2, r90g3, raaa, rac2, rf1, rootw, rr1, rr1s, rr3, rr3s, s1, &
        & s2, s2a, s2b, s2i, s3, scale, swa, swr0, swt, t1, t2, t2a, t3, una2p, w, y0b1, y0b2, y0d1, &
        & y0d2, y0r1, y0r2, y90b1, y90b2, y90d1, y90d2, y90r1, y90r2
    integer :: i, nt
    dimension r(nt,3), e(nt)
    save
    data scale/1.5d0/
    data c1s1/203.128d0/, c2s1/1.6167d0/, c3s1/-100.124d0/, &
        & c4s1/68.078d0/, c5s1/-15.832d0/, c6s1/1.303d0/
    data dnaht1/3.35000d0/, dnaht2/10.000d0/, dnaht3/54.6393d0/, &
        & dnaht4/-17.63343d0/, dnaht5/3.0085d0/, dnaht6/1.450d0/, &
        & dnaht7/-2.87977d0/, dnaht61/1.300d0/, dnaht11/2.85000d0/
    data deh2/4.7404422d0/, una2p/2.1037d0/, c2a/0.8d0/, c2b/0.8d0/, &
        & c2c/0.15d0/, cprime/8.0d-1/
    data y0d1/2.417800d-1/, y0r1/3.708055d0/, y0b1/6.960127d-1/, &
        & y90d1/1.372065d-2/, y90r1/7.233552d0/, y90b1/5.127734d-1/, &
        & y0d2/8.870172d-1/, y0r2/2.300265d0/, y0b2/8.497374d-1/, &
        & y90d2/1.106187d0/, y90r2/2.002516d0/, y90b2/9.514887d-1/, &
        & swa/0.78495d0/, swr0/7.47654d0/
    data a0g1/-8.210245d-2/, b0g1/4.658698d0/, r0g1/4.552860d0/, &
        & a0g2/3.546427d-2/, b0g2/6.456787d-1/, r0g2/6.332225d0/, &
        & a0g3/-1.562327d-2/, b0g3/1.170537d-1/, r0g3/1.097310d+1/
    data a90g1/-9.929297d-2/, b90g1/4.367220d0/, r90g1/4.402447d0/, &
        & a90g2/4.696181d-2/, b90g2/5.470459d-1/, r90g2/6.289526d0/, &
        & a90g3/-1.219900d-2/, b90g3/1.002291d-2/, r90g3/1.118725d+1/
    data h2t1/1.5d0/, h2t2/4.47009d0/, h2t3/4.48534d0/, &
        & h2t4/1.77908d0/, h2t5/-6.46276d0/, h2t6/3.15396d0/, &
        & h2t7/-1.59912d0/
    rac2  = 1.d0/dsqrt(2.d0)
    cevau = 1.d0/27.2113961d0

    return
    !
    entry pot22(r, e, nt)

    do 10 i=1, nt

        !***********************************************************************
        r1s  = r(i,1)*r(i,1)
        r2s  = r(i,2)*r(i,2)
        r3s  = r(i,3)*r(i,3)
        rr1  = r(i,1)*scale
        rr3  = r(i,3)*scale
        rr1s = rr1**2
        rr3s = rr3**2
        !***********************************************************************
        !*****SWITCHING FUNCTION FOR THE H2 TRIPLET, AND ITS DERIVATIVES:

        swt = 0.5d0*(1.d0-dtanh( swa*(0.5d0*(rr1+rr3)-swr0)))
        !***********************************************************************

        !***********************************************************************
        !*****SINGLETS AND THEIR DERIVATIVES:

        !-----------------------------------------------------------------------
        !*****MH SINGLET:
        raaa     = 2.700d0
        aaa      = 2.40d0
        aaa1     = 2.60d0
        daa      = 1.80d0
        hlpnah1  = dexp(-aaa1*(r(i,1)-raaa))
        hlpnah11 = dexp(-1.65d0*aaa*(r(i,1)-raaa))
        s1       = 0.30d0*daa*(hlpnah11-2.d0*hlpnah1)
        !-----------------------------------------------------------------------
        !*****H2 SINGLET
        s2a=(109.384d0+463.033d0*r(i,2))*dexp(-6.240d0*r(i,2)) -r2s*(13.654d0+20.717d0*r(i,2))*dexp(-2.032d0*r(i,2))+ &
            & una2p

        b    = 1.4d0
        a    = 2.d0*b
        r00  = 1.411655d0
        dn   = -4.74044d0/(1.d0-a/b)
        hlp1 = dexp(-a*(r(i,2)-r00))
        hlp2 = dexp(-b*(r(i,2)-r00))
        s2i  = dn*(hlp1-(a/b)*hlp2)+una2p

        hlp = 0.5d0*(r(i,1)+r(i,3)-7.0d0)
        hlp = 0.5d0*(1.d0-dtanh(hlp))
        s2a = s2a+(s2i-s2a)*hlp

        !-----------------------------------------------------------------------
        !*****MF SINGLET:
        hlpnah3  = dexp(-aaa*(r(i,3)-raaa))
        hlpnah33 = dexp(-1.65d0*aaa*(r(i,3)-raaa))
        s3       = daa*(hlpnah33-2.d0*hlpnah3)
        !***********************************************************************
        !*****TRIPLETS AND THEIR DERIVATIVES:

        !-----------------------------------------------------------------------
        !*****H2 TRIPLET (interaction region):

        t2a = 0.5d0*(dexp(-h2t2*(r(i,2)-h2t1))+(h2t3+h2t5*r(i,2) +h2t6*r2s+h2t7*r2s*r(i,2))*dexp(-h2t4*r(i,2)))

        !*****H2 TRIPLET (asymptotic)
        s2b = 147.214d0*dexp(-3.915d0*r(i,2))+r(i,2)*(41.275d0+10.505d0 *r(i,2)+4.408d0*r2s)*dexp(-2.032d0*r(i,2))

        !*****H2 TRIPLET: TOTAL (switched between the asymptotic and the ineract
        t2 = swt*(t2a-s2b)+s2b
        !-----------------------------------------------------------------------
        !*****TRIPLETS AND THEIR DERIVATIVES:
        !!!!!!MH TRIPLET:
        t1 = 3.d0*(dexp(-dnaht2*(rr1-dnaht11))+ (dnaht7+dnaht3*rr1+dnaht4*rr1s+dnaht5*rr1*rr1s) *dexp(-dnaht61*rr1))
        !-----------------------------------------------------------------------
        !!!!!!MF TRIPLET: (should be the same as the first one, t1)

        t3 = 3.d0*(dexp(-dnaht2*(rr3-dnaht1))+ (dnaht7+dnaht3*rr3+dnaht4*rr3s+dnaht5*rr3*rr3s) *dexp(-dnaht6*rr3))
        !-----------------------------------------------------------------------
        !***********************************************************************
        !*****TOTAL H2 SINGLET:
        !!!!!!H2 SINGLET: total (switched between the true singlet and a triplet
        !                        in the asymptotic region at large HH distances)
        hlp = dtanh((s2a-t2)/cprime)
        s2  = 0.5d0*(s2a+t2-hlp*(s2a-t2))
        !***********************************************************************
        !***********************************************************************
        !-----------------------------------------------------------------------
        !*****COULOMB INTEGRALS AND THEIR FIRST DERIVATIVES:

        coul1 = 0.5d0*(s1+t1)
        coul2 = 0.5d0*(s2+t2)
        coul3 = 0.5d0*(s3+t3)
        !-----------------------------------------------------------------------
        !*****EXCHANGE INTEGRALS AND THEIR FIRST DERIVATIVES:

        exch1 = 0.5d0*(s1-t1)
        exch2 = 0.5d0*(s2-t2)
        exch3 = 0.5d0*(s3-t3)
        !***********************************************************************

        w     = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
        cplg2 = c2a*dexp(-c2b*w-c2c*(rr1+r(i,2)+rr3))
        rootw = rac2*dsqrt(w+cplg2*cplg2)
        e(i)  = coul1+coul2+coul3+deh2-rootw

        !***********************************************************************
        rf1  = 3.00d0
        af1  = 0.60d0
        fff  = 1.d0 - 0.5d0*(1.d0+dtanh(af1*(r(i,3)-rf1)))
        e(i) = ((e(i)-1.35d0)/(1.00d0+1.80d0*fff))*cevau
        !***********************************************************************
        !***********************************************************************
        !***********************************************************************

    10 continue

end subroutine prepot22
