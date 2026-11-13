!> @brief Implement H3_GEN_BKMP3_1996, a potential energy surface for H3.
module chempotpy_h3_gen_bkmp3_1996_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_BKMP3_1996 potential energy surface for H3.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,              only: hartree_to_electronvolt
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
        integer               :: istate
        real(wp)              :: r(1,3), e(1), v
        integer               :: iatom, idir, i, j
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        ! input cartesian is HHH
        call triatomic_distances_to_bohr(x, r(1,:))
        call pot(r, e, 1, 1)
        v = e(1)*hartree_to_electronvolt
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = v
            enddo
        else
            write (*, *) 'Only energy is available'
        endif
    end subroutine

    !> @brief Evaluate a batch of BKMP3 geometries.
    subroutine pot( rvp, evp, nt, nsurf )
        use chempotpy_kinds, only: wp
        implicit none
        integer, intent(in)   :: nt, nsurf
        real(wp), intent(in)  :: rvp(nt,3)
        real(wp), intent(out) :: evp(nt)
        real(wp)              :: eee, eshift, vp(3), xxx(3)
        integer               :: ivp, ivp1
        eshift = -0.174495770894e0_wp
        do ivp = 1, nt
            do ivp1 = 1, 3
                xxx(ivp1) = rvp(ivp,ivp1)
            enddo
            call bkmp2(xxx, eee, vp, -1)
            evp(ivp) = eee-eshift
        enddo
    end subroutine pot

    !> @brief Evaluate surfid.
    subroutine surfid( h3id )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: h3id
        h3id = 950621.0
    end
    !> @brief Evaluate bkmp2.
    subroutine bkmp2( r, vtot, dvtot, id )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: cal, cas, cb(2,25), cbnda, cbndb, dcal(3), dcas(3), dcbnda(3), dcbndb(3), dt(3), dvas(3), dvbnda(3), &
            & dvbndb(3), dvlon(3), dvtot(3), r(3), t(3), vas, vb(2, 25), vbnda, vbndb, vlon, vtot
        integer :: i, icompc, id, ipr
        ipr   = 0
        vlon  = 0.e0_wp
        vas   = 0.e0_wp
        vbnda = 0.e0_wp
        vbndb = 0.e0_wp
        cal   = 0.e0_wp
        cas   = 0.e0_wp
        cbnda = 0.e0_wp
        cbndb = 0.e0_wp
        do i=1, 3
            dvlon(i)  = 0.e0_wp
            dvas(i)   = 0.e0_wp
            dvbnda(i) = 0.e0_wp
            dvbndb(i) = 0.e0_wp
            dcal(i)   = 0.e0_wp
            dcas(i)   = 0.e0_wp
            dcbnda(i) = 0.e0_wp
            dcbndb(i) = 0.e0_wp
        enddo
        do i=1, 25
            vb(1, i) = 0.e0_wp
            vb(2, i) = 0.e0_wp
            cb(1,i)  = 0.e0_wp
            cb(2,i)  = 0.e0_wp
        enddo
        call h3lond95(r, vlon, dvlon)
        call vascal95(r, vas, dvas)
        call compac95(r, icompc, t, dt)
        if( icompc .ge. 1 )then
            call csym95  (r, cal, dcal)
            call casym95 (r, cas, dcas, 1, t, dt)
        endif
        call vbcb95(r, icompc, t, dt, ipr, &
            &          vbnda, vbndb, dvbnda, dvbndb, &
            &          cbnda, cbndb, dcbnda, dcbndb)
        vtot     = vlon + vas + vbnda + vbndb + cal + cas + cbnda + cbndb
        dvtot(1) = dvlon(1) + dvas(1) + dvbnda(1) + dvbndb(1) +  dcal(1) + dcas(1) + dcbnda(1) + dcbndb(1)
        dvtot(2) = dvlon(2) + dvas(2) + dvbnda(2) + dvbndb(2) +  dcal(2) + dcas(2) + dcbnda(2) + dcbndb(2)
        dvtot(3) = dvlon(3) + dvas(3) + dvbnda(3) + dvbndb(3) +  dcal(3) + dcas(3) + dcbnda(3) + dcbndb(3)
        if( ipr .gt. 0 )then
            write(7, *)
            write(7, *) 'h3tot enter --------------------------------'
            write(7, '(5x,''    r ='',3(1x,f16.10),'' icompac = '',i1,/, 5x,''Vtot  = '',f18.12,/, 5x,''Vlon  = '',f18.12,' // &
                '''       Vas   = '',g18.12,/, 5x,''Cal   = '',g18.12,''       Cas   = '',g18.12,/, 5' // &
                'x,''Vbnda = '',g18.12,''       Vbndb = '',g18.12,/, 5x,''Cbnda = '',g18.12,''       ' // &
                'Cbndb = '',g18.12)') r, icompc, vtot, &
                &             vlon, vas, cal, cas, vbnda, vbndb, cbnda, cbndb
            if( ipr .gt. 1 )then
                write(7, '('' dVtot   = '',3(1x,g18.12),/, '' dVlon   = '',3(1x,g18.12),/, '' dVas    = '',3(1x,g18.12),/,' // &
                    ' '' dVbnda  = '',3(1x,g18.12),/, '' dVbndb  = '',3(1x,g18.12),/, '' dCal    = '',3(1' // &
                    'x,g18.12),/, '' dCas    = '',3(1x,g18.12),/, '' dCbnda  = '',3(1x,g18.12),/, '' dCbn' // &
                    'db  = '',3(1x,g18.12))') dvtot, dvlon, dvas, dvbnda, dvbndb, &
                    &                           dcal, dcas, dcbnda, dcbndb
            endif
            write(7, '(''Vba values: '',5(1x,g12.6),/,''Vbb values: '',5(1x,g12.6))') (vb(1, i), i=1, 5), (vb(2, i), i=1, 5)
            write(7, '(''Cba values: '',3(1x,f16.8),/, ''            '',3(1x,f16.8),/, ''            '',3(1x,f16.8),/, ''C' // &
                'bb values: '',3(1x,f16.8),/, ''            '',3(1x,f16.8),/, ''            '',3(1x,f' // &
                '16.8))') (cb(1,i), i=1, 9), (cb(2,i), i=1, 9)
            write(7, '(''   Vbnda        Vbndb        Cbnda        Cbndb '',/, 4(e12.6,2x),/, ''dVbnda: '',3(e16.10,1x),/,' // &
                ' ''dVbndb: '',3(e16.10,1x),/, ''dCbnda: '',3(e16.10,1x),/, ''dCbndb: '',3(e16.10,1x)' // &
                ')') vbnda, vbndb, cbnda, cbndb, &
                &                  dvbnda, dvbndb, dcbnda, dcbndb
            write(7, *) 'exiting subr.h3tot'
            write(7, *) 'h3tot exit ---------------------------------'
        endif
    end
    !> @brief Evaluate triplet95.
    subroutine triplet95( r, e3, id )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: dr, e1(3), e3(3), exdr, r, ra6, ra61, rsq
        integer             :: id, ipr
        real(wp), parameter :: rl = 0.95e0_wp
        real(wp), parameter :: rr = 1.15e0_wp
        real(wp), parameter :: z1 = 1.e0_wp
        real(wp), parameter :: z2 = 2.e0_wp
        real(wp), parameter :: z4 = 4.e0_wp
        real(wp), parameter :: a1 = -0.0298546962
        real(wp), parameter :: a2 = -23.9604445036
        real(wp), parameter :: a3 = -42.5185569474
        real(wp), parameter :: a4 = 2.0382390988
        real(wp), parameter :: a5 = -11.5214861455
        real(wp), parameter :: a6 = 1.5309487826
        real(wp), parameter :: c1 = -0.4106358351531854
        real(wp), parameter :: c2 = -0.0770355790707090
        real(wp), parameter :: c3 = 0.4303193846943223  !jun21 fit
        ipr   = 0
        e3(3) = 0.e0_wp
        if( r .ge. rr )then
            exdr  = exp( -a4*r )
            rsq   = r*r
            ra6   = r**(-a6)
            e3(1) = a1* ( a2 + r + a3*rsq + a5*ra6 )*exdr
            ra61  = r**(-a6-z1)
            e3(2) = a1*exdr* ( z1 -a2*a4 +(z2*a3-a4)*r -a3*a4*rsq -a5*a6*ra61 -a4*a5*ra6 )
        endif
        if( r .lt. rr )then
            dr = r- rl
            call vh2opt95(r, e1, 2)
            if( r .le. rl )then
                e3(1) = e1(1) + c2*dr + c3
                e3(2) = e1(2) + c2
            else
                e3(1) = e1(1) + c1*dr*dr*dr + c2*dr + c3
                e3(2) = e1(2) + 3.e0_wp*c1*dr*dr + c2
            endif
        endif
        if( ipr .gt. 0 ) write(7, '(''  triplet:  E3 = '',3(1x,f12.8))') e3
    end
    !> @brief Evaluate vH2opt95.
    subroutine vh2opt95( r, e, ideriv )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: alphar, dalphr, ddalph, ddvlr, ddvsr, dvlr, dvsr, e(3), exalph, r, r02, r04, r06, r10, r11, r12, &
            & r13, r14, r15, r2, r3, r4, r5, r6, r7, r8, r9, rr2, rr25, rr26, rr27, rr4, rr43, rr44, rr6, rr62, rr63, vlr, vsr
        integer             :: ideriv
        real(wp), parameter :: a0  = 0.03537359271649620
        real(wp), parameter :: a1  = 2.013977588700072
        real(wp), parameter :: a2  = -2.827452449964767
        real(wp), parameter :: a3  = 2.713257715593500
        real(wp), parameter :: a4  = -2.792039234205731
        real(wp), parameter :: a5  = 2.166542078766724
        real(wp), parameter :: a6  = -1.272679684173909
        real(wp), parameter :: a7  = 0.5630423099212294
        real(wp), parameter :: a8  = -0.1879397372273814
        real(wp), parameter :: a9  = 0.04719891893374140
        real(wp), parameter :: a10 = -0.008851622656489644
        real(wp), parameter :: a11 = 0.001224998776243630
        real(wp), parameter :: a12 = -1.227820520228028e-04_wp
        real(wp), parameter :: a13 = 8.638783190083473e-06_wp
        real(wp), parameter :: a14 = -4.036967926499151e-07_wp
        real(wp), parameter :: a15 = 1.123286608335365e-08_wp
        real(wp), parameter :: a16 = -1.406619156782167e-10_wp
        real(wp), parameter :: r0  = 3.5284882e0_wp
        real(wp), parameter :: dd  = 0.160979391e0_wp
        real(wp), parameter :: c6  = 6.499027e0_wp
        real(wp), parameter :: c8  = 124.3991e0_wp
        real(wp), parameter :: c10 = 3285.828e0_wp
        e(1) = -999.e0_wp
        e(2) = -999.e0_wp
        e(3) = -999.e0_wp
        r2   = r  *r
        r3   = r2 *r
        r4   = r3 *r
        r5   = r4 *r
        r6   = r5 *r
        r7   = r6 *r
        r8   = r7 *r
        r9   = r8 *r
        r10  = r9 *r
        r11  = r10*r
        r12  = r11*r
        r13  = r12*r
        r14  = r13*r
        r15  = r14*r
        r02  = r0*r0
        r04  = r02*r02
        r06  = r04*r02
        rr2  = r2 + r02
        rr4  = r4 + r04
        rr6  = r6 + r06
        rr25 = rr2*rr2*rr2*rr2*rr2
        alphar =  a0/r + a1 + a2 *r   + a3 *r2  + a4 *r3  + a5 *r4  + a6 *r5 + a7 *r6  + a8 *r7  + a9 *r8  + a10*r9  + &
            & a11* r10 + a12*r11 + a13*r12 + a14*r13 + a15*r14 + a16*r15
        exalph = exp(alphar)
        vsr    = dd*(exalph-1.e0_wp)*(exalph-1.e0_wp) - dd
        vlr    = -c6/rr6   -c8/(rr4*rr4)   -c10/rr25
        e(1)   = vsr +  vlr
        if( ideriv .ge. 1 )then
            r3   = r2*r
            r5   = r4*r
            rr26 = rr25*rr2
            rr43 = rr4*rr4*rr4
            dalphr = -a0/r2 + a2 + 2.e0_wp*a3*r +  3.e0_wp *a4 *r2  +  4.e0_wp *a5 *r3  +  5.e0_wp *a6 *r4 + 6.e0_wp * &
                & a7 *r5  + 7.e0_wp * a8 *r6  +  8.e0_wp *a9 *r7 +  9.e0_wp *a10*r8  + 10.e0_wp *a11*r9  + 11.e0_wp * &
                & a12*r10 + 12.e0_wp *a13*r11 + 13.e0_wp * a14* r12 + 14.e0_wp *a15*r13 + 15.e0_wp *a16*r14
            dvsr = 2.e0_wp *dd *(exalph-1.e0_wp) *exalph *dalphr
            dvlr = 6.e0_wp*c6*r5 / (rr6*rr6) +  8.e0_wp*c8*r3 / rr43  + 10.e0_wp*c10*r / rr26
            e(2) = dvsr + dvlr
        endif
        if( ideriv .ge. 2 )then
            r10  = r6*r4
            rr27 = rr26*rr2
            rr44 = rr43*rr4
            rr62 = rr6*rr6
            rr63 = rr62*rr6
            ddalph = 2.e0_wp*a0/r3 +  2.e0_wp*a3      +  6.e0_wp*a4 *r + 12.e0_wp*a5 *r2    + 20.e0_wp*a6 *r3  + &
                & 30.e0_wp*a7 *r4 + 42.e0_wp* a8 * r5    + 56.e0_wp*a9 *r6  + 72.e0_wp*a10*r7 + 90.e0_wp*a11*r8    + &
                & 110.e0_wp*a12*r9  +132.e0_wp*a13*r10 + 156.e0_wp*a14* r11   + 182.e0_wp*a15*r12 +210.e0_wp*a16*r13
            ddvsr = 2.e0_wp *dd *exalph *( (2.e0_wp*exalph-1.e0_wp)*dalphr*dalphr + (exalph-1.e0_wp)*ddalph )
            ddvlr =- 72.e0_wp* c6 *r10 / rr63 -96.e0_wp *c8 *r6 / rr44 -120.e0_wp*c10 *r2  / rr27 +30.e0_wp *c6 *r4 / &
                & rr62 + 24.e0_wp* c8 * r2  / rr43 +10.e0_wp *c10    / rr26
            e(3) = ddvsr + ddvlr
        endif
    end
    !> @brief Evaluate h3lond95.
    subroutine h3lond95( r, vlon, dvlon )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: de1(3), de3(3), dvlon(3), e1(3), e3(3), esing(3), etrip(3), r(3), rootjt, sumj, sumq, vlon
        integer             :: i, ipr
        real(wp)            :: q(3), j(3), jt
        real(wp), parameter :: half = 0.5e0_wp
        real(wp), parameter :: two  = 2.e0_wp
        real(wp), parameter :: eps2 = 1.e-12_wp
        ipr = 0
        do i=1, 3
            call vh2opt95(r(i), esing, 2)
            e1(i)  = esing(1)
            de1(i) = esing(2)
            call triplet95(r(i), etrip, 2)
            if( ipr .gt. 0 ) write(7, '(''from subr.london: '',/, ''  using r'',i1,'':   Esinglet='',3(1x,f12.8),/, ''        ' // &
                ' '',1x,''    Etriplet='',3(1x,f12.8))') i, esing, etrip
            e3(i)  = etrip(1)
            de3(i) = etrip(2)
            q(i)   = half*(e1(i) + e3(i))
            j(i)   = half*(e1(i) - e3(i))
        enddo
        sumq   = q(1) + q(2) + q(3)
        sumj   = abs( j(2)-j(1) )**2 + abs( j(3)-j(2) )**2 + abs( j(3)-j(1) )**2
        jt     = half*sumj + eps2
        rootjt = sqrt(jt)
        vlon   = sumq - rootjt
        if( ipr .gt. 0 )then
            write(7, '(''         sumQ = '',g12.6,''        sumJ = '',g12.6)') sumq, sumj
            write(7, '(''         Vlon = '',f12.8,''      rootJt = '',g12.6)') vlon, rootjt
        endif
        dvlon(1) = half*(de1(1)+de3(1)) - 0.25e0_wp*(two*j(1)-j(2)-j(3))*(de1(1)-de3(1))/rootjt
        dvlon(2) = half*(de1(2)+de3(2)) - 0.25e0_wp*(two*j(2)-j(3)-j(1))*(de1(2)-de3(2))/rootjt
        dvlon(3) = half*(de1(3)+de3(3)) - 0.25e0_wp*(two*j(3)-j(1)-j(2))*(de1(3)-de3(3))/rootjt
        if( ipr .gt. 0 )then
            write(7, '(''             r = '',3(1x,f12.6),/, ''            E1 = '',3(1x,f12.8),/, ''            E3 = '',3(1' // &
                'x,f12.8),/, ''          vlon = '',1x,f12.8)') r, e1, e3, vlon
            write(7, '(13x,''q = '',3(1x,e12.6),/,13x,''j = '',3(1x,e12.6))') q, j
            write(7, '(''         dVlon = '',3(1x,g12.6))') dvlon
        endif
    end
    !> @brief Evaluate Vascal95.
    subroutine vascal95( rpass, vas, dvas )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: a, a2, a3, a4, a5, da(3), ds(3), dvas(3), exp1, exp6, r, r1, r2, r3, rcu, rpass(3), rsq, s, vas
        integer             :: ipr
        real(wp), parameter :: aa1 = 0.3788951192e-02
        real(wp), parameter :: aa2 = 0.1478100901e-02
        real(wp), parameter :: aa3 = -.1848513849e-03
        real(wp), parameter :: aa4 = 0.9230803609e-05
        real(wp), parameter :: aa5 = -.1293180255e-06
        real(wp), parameter :: aa6 = 0.5237179303e+00
        real(wp), parameter :: aa7 = -.1112326215e-02  !jun21 fit
        ipr = 0
        r1  = rpass(1)
        r2  = rpass(2)
        r3  = rpass(3)
        r   = r1 + r2 + r3
        rsq = r*r
        rcu = rsq*r
        call acalc95(r1, r2, r3, a, da)
        a2      = a *a
        a3      = a2*a
        a4      = a3*a
        a5      = a4*a
        exp1    = exp(-aa1*rcu)
        exp6    = exp(-aa6*r)
        s       = aa2*a2 + aa3*a3 + aa4*a4 + aa5*a5
        vas     = s*exp1 +  aa7*a2 *exp6 / r
        ds(1)   = ( 2.e0_wp*aa2*a  +3.e0_wp*aa3*a2 +4.e0_wp*aa4*a3 +5.e0_wp*aa5*a4) * da(1)
        dvas(1) = -3.e0_wp*aa1*rsq*s*exp1 + ds(1)*exp1 -aa7*a2*exp6/rsq  + 2.e0_wp*aa7*a*da(1)*exp6/r -aa6*aa7*a2*exp6/r
        ds(2)   = ( 2.e0_wp*aa2*a  +3.e0_wp*aa3*a2 +4.e0_wp*aa4*a3 +5.e0_wp*aa5*a4) * da(2)
        dvas(2) = -3.e0_wp*aa1*rsq*s*exp1 + ds(2)*exp1 -aa7*a2*exp6/rsq  + 2.e0_wp*aa7*a*da(2)*exp6/r -aa6*aa7*a2*exp6/r
        ds(3)   = ( 2.e0_wp*aa2*a  +3.e0_wp*aa3*a2 +4.e0_wp*aa4*a3 +5.e0_wp*aa5*a4) * da(3)
        dvas(3) = -3.e0_wp*aa1*rsq*s*exp1 + ds(3)*exp1 -aa7*a2*exp6/rsq  + 2.e0_wp*aa7*a*da(3)*exp6/r -aa6*aa7*a2*exp6/r
        if( ipr .gt. 0 ) write(7, '(''    from subr.vascal:   Vas  = '',  1x,g12.6, /, ''                        dS   = '',3(1' // &
            'x,g12.6),/, ''                        dVas = '',3(1x,g12.6))') vas, ds, dvas
    end
    !> @brief Evaluate acalc95.
    subroutine acalc95( r1, r2, r3, a, da )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a, da(3), r1, r2, r3
        integer  :: ipr
        ipr   = 0
        a     = (r1-r2)*(r2-r3)*(r3-r1)
        da(1) = ( -2.e0_wp*r1 + r2 + r3 )*(r2-r3)
        da(2) = ( -2.e0_wp*r2 + r3 + r1 )*(r3-r1)
        da(3) = ( -2.e0_wp*r3 + r1 + r2 )*(r1-r2)
        if( a .lt. 0.e0_wp )then
            a     = -a
            da(1) = -da(1)
            da(2) = -da(2)
            da(3) = -da(3)
        endif
        if( ipr .gt. 0 ) write(7, '(''Acalc enter ---------------------------------------'',/, ''     A = '',f15.10,''   dA = ' // &
            ''',3(1x,g12.6),/, ''Acalc exit ----------------------------------------'')') a, da
    end
    !> @brief Evaluate compac95.
    subroutine compac95( r, icompc, t, dt )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: bot, bot2, dt(3), r(3), t(3), top, top2, top3
        integer             :: i, icompc
        real(wp), parameter :: rr = 1.15e0_wp
        real(wp), parameter :: rp = 1.25e0_wp
        do i=1, 3
            t(i)  = 0.e0_wp
            dt(i) = 0.e0_wp
        enddo
        icompc = 0
        if( r(1) .lt. rr ) icompc = icompc + 1
        if( r(2) .lt. rr ) icompc = icompc + 1
        if( r(3) .lt. rr ) icompc = icompc + 1
        if( icompc .eq. 0 ) return
        do i=1, 3
            if( r(i) .lt. rr )then
                top   = rr-r(i)
                bot   = rp-r(i)
                top2  = top  * top
                top3  = top2 * top
                bot2  = bot  * bot
                t(i)  = top3/bot
                dt(i) = -3.e0_wp*top2/bot + top3/bot2
            endif
        enddo
    end
    !> @brief Evaluate casym95.
    subroutine casym95( r, cas, dcas, id, t, dt )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: a, a2, cas, da(3), dcas(3), determ(3), dpr(3), dseries(3), dsumt(3), dt(3), dterm1(3), eterm, pr, &
            & pr2, pr3, r(3), series, sr, sr2, sr3, sumt, t(3), term1
        integer             :: i, id, ipr
        real(wp), parameter :: u1  = 0.2210243144e+00
        real(wp), parameter :: u2  = 0.4367417579e+00
        real(wp), parameter :: u3  = 0.6994985432e-02
        real(wp), parameter :: u4  = 0.1491096501e+01
        real(wp), parameter :: u5  = 0.1602896673e+01
        real(wp), parameter :: u6  = -.2821747323e+01
        real(wp), parameter :: u7  = 0.4948310833e+00
        real(wp), parameter :: u8  = -.3540394679e-01
        real(wp), parameter :: u9  = -.3305809954e+01
        real(wp), parameter :: u10 = 0.3644382172e+01
        real(wp), parameter :: u11 = -.9997570970e+00
        real(wp), parameter :: u12 = 0.7989919534e-01
        real(wp), parameter :: u13 = -.1075807322e-02  !jun21 fit
        ipr     = 0
        cas     = 0.e0_wp
        dcas(1) = 0.e0_wp
        dcas(2) = 0.e0_wp
        dcas(3) = 0.e0_wp
        call acalc95(r(1), r(2), r(3), a, da)
        a2     = a*a
        sumt   = t(1) + t(2) + t(3)
        sr     = r(1) + r(2) + r(3)
        pr     = r(1) * r(2) * r(3)
        sr2    = sr*sr
        sr3    = sr2*sr
        pr2    = pr*pr
        pr3    = pr2*pr
        series = 1.e0_wp + u4/pr2 +  u5/pr + u6 +   u7*pr +  u8*pr2 + a*(u9/pr2 + u10/pr + u11 + u12*pr + u13*pr2)
        term1  = u1/pr**u2
        eterm  = exp(-u3*sr3)
        cas    = sumt * a2 * term1 * series * eterm
        if( ipr .gt. 0 )then
            write(7, *)
            write(7, *) 'Casym: ------------------'
            write(7, *) ' T(i) = ', t
            write(7, *) 'dT(i) = ', dt
            write(7, *) ' id = ', id
            write(7, '(''subr.Cas:   series = '',g12.6,''     term1 = '',g12.6,/, ''             eterm = '',g12.6,''      ' // &
                'Cas  = '',g12.6)') series, term1, eterm, cas
            write(7, *) ' Cas = ', cas
        endif
        if( id .gt. 0 )then
            dpr(1) = r(2)*r(3)
            dpr(2) = r(3)*r(1)
            dpr(3) = r(1)*r(2)
            do i=1, 3
                dsumt(i)  = dt(i)
                dterm1(i) = -1.e0_wp*u1*u2*pr**(-u2-1.e0_wp)*dpr(i)
                determ(i) = eterm *(-3.e0_wp*u3*sr2)
                dseries(i)= dpr(i)*( -2.e0_wp*u4/pr3  -u5/pr2 +u7  +2.e0_wp*u8*pr    ) +da(i)*( u9/pr2 +u10/pr +u11 + &
                    & u12* pr + u13* pr2    ) +a*dpr(i)*( -2.e0_wp*u9/pr3 -u10/pr2 +u12 + 2.e0_wp*u13*pr )
                dcas(i) = dsumt(i)      * a2    * term1 * series * eterm + 2.e0_wp*a*da(i)  * sumt  * term1 * series * &
                    & eterm + dterm1(i)     * sumt  * a2    * series * eterm + dseries(i)    * sumt  * a2    * &
                    & term1  * eterm + determ(i)     * sumt  * a2    * term1  * series
            enddo
            if( ipr .gt. 0 )then
                write(7, '(''    derivatives: dPr = '',3(1x,g12.6),/, ''                 dT  = '',3(1x,g12.6),/, ''       ' // &
                    '       dterm1 = '',3(1x,g12.6),/, ''              determ = '',3(1x,g12.6),/, ''     ' // &
                    '        dseries = '',3(1x,g12.6),/, ''              dCas   = '',3(1x,g12.6))') dpr, dt, dterm1, determ, dseries, dcas
                write(7, *) ' dCas = ', dcas
            endif
        endif
    end
    !> @brief Evaluate csym95.
    subroutine csym95( r, cal, dcal )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: cal, dcal(3), dexp3, dg(3), exp3, g(3), r(3), ri, rpri, rpri2, rrri, rrri2, rrri3, sr, sr2, sr3, &
            & sumg, sumv(3)
        integer             :: i, ipr
        real(wp), parameter :: rr = 1.15e0_wp
        real(wp), parameter :: rp = 1.25e0_wp
        real(wp), parameter :: v1 = -.2071708868e+00
        real(wp), parameter :: v2 = -.5672350377e+00
        real(wp), parameter :: v3 = 0.9058780367e-02  !jun21 fit
        cal   = 0.e0_wp
        ipr   = 0
        sr    = r(1)+r(2)+r(3)
        sr2   = sr*sr
        sr3   = sr*sr2
        exp3  = exp( -v3*sr3 )
        dexp3 = -3.e0_wp*v3*sr2*exp3
        do i=1, 3
            ri      = r(i)
            rrri    = rr-ri
            rrri2   = rrri*rrri
            rrri3   = rrri*rrri2
            rpri    = rp-ri
            rpri2   = rpri*rpri
            g(i)    = 0.e0_wp
            dg(i)   = 0.e0_wp
            sumv(i) = v1+v1*v2*ri
            if( ri .lt. rr )then
                g(i)  = (rrri3/rpri) *sumv(i)
                dg(i) = (rrri3/rpri2)*sumv(i) -3.e0_wp*(rrri2/rpri)*sumv(i) + (rrri3/rpri)*v1*v2
            endif
        enddo
        sumg    = g(1) + g(2) + g(3)
        cal     = sumg*exp3
        dcal(1) = dg(1)*exp3 + sumg*dexp3
        dcal(2) = dg(2)*exp3 + sumg*dexp3
        dcal(3) = dg(3)*exp3 + sumg*dexp3
        if( ipr .gt. 0 )then
            write(7, *)
            write(7, '(''Csym enter ----------------------------------------'',/, ''    Cal = '',1x,g20.14)') cal
            write(7, '(''      G = '',3(1x,g12.6),/, ''     dG = '',3(1x,g12.6),/, ''   sumv = '',3(1x,g12.6),/, ''   dCal' // &
                ' = '',3(1x,g16.10),/, ''Csym exit -----------------------------------------'')') g, dg, sumv, dcal
            write(7, *) 'exiting subr.csym'
        endif
    end
    !> @brief Evaluate vbcb95.
    subroutine vbcb95( rpass, icompc, t, dt, ipr, &
            &                  vbnda, vbndb, dvbnda, dvbndb, &
            &                  cbnda, cbndb, dcbnda, dcbndb )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: aasum, asum, b1, b12, b13, b14, b15, b1a, b1b, b2, b3, b3b, bsum, c1, c1cube, c2, c2cube, c3, &
            & c3cube, cb(2,25), cbnda, cbndb, cbnds(2), cos3t1, cos3t2, cos3t3, csum, cx1, d1, d2, d3, daasum, dasum, db1a(3), &
            & db1b(3), db2(3), db3(6), dbsum, dc1dr1, dc1dr2, dc1dr3, dc2dr1, dc2dr2, dc2dr3, dc3dr1, dc3dr2, dc3dr3, dcb1(3), &
            & dcb2(3), dcb3(3), dcb4(3), dcb5(3), dcb6(3), dcb7(3), dcb8(3), dcbnda(3), dcbndb(3), dcbnds(2, 3), ddsum1, ddsum2, &
            & dexp1, dexp2, dexp7, dffsum, dfsum, dgsum, dp(3), dsum1, dsum2, dt(3), dvb1(3), dvb2(3), dvb3(3), dvb4(3), &
            & dvb5(3), dvbnda(3), dvbndb(3), dx1, eps2, exp1, exp2, exp7, ffsum, fsum, gsum, p, pcu, psq, r, r1, r2, r3, rcu, &
            & rpass(3), rsq, sum, sumb, sumt, t(3), t1, t2, t3, vb(2, 25), vbnd(2), vbnda, vbndb
        integer             :: icompc, ipr
        real(wp), parameter :: z58   = 0.625e0_wp
        real(wp), parameter :: z38   = 0.375e0_wp
        real(wp), parameter :: beta1 = 0.52e0_wp
        real(wp), parameter :: beta2 = 0.052e0_wp
        real(wp), parameter :: beta3 = 0.79e0_wp
        real(wp), parameter :: a11   = -.1838073394e+03
        real(wp), parameter :: a12   = 0.1334593242e+02
        real(wp), parameter :: a13   = -.2358129537e+00
        real(wp), parameter :: a21   = -.4668193478e+01
        real(wp), parameter :: a22   = 0.7197506670e+01
        real(wp), parameter :: a23   = 0.2162004275e+02
        real(wp), parameter :: a24   = 0.2106294028e+02
        real(wp), parameter :: a31   = 0.4242962586e+01
        real(wp), parameter :: a32   = 0.4453505045e+01
        real(wp), parameter :: a41   = -.1456918088e+00
        real(wp), parameter :: a42   = -.1692657366e-01
        real(wp), parameter :: a43   = 0.1279520698e+01
        real(wp), parameter :: a44   = -.4898940075e+00
        real(wp), parameter :: a51   = 0.1742295219e+03
        real(wp), parameter :: a52   = 0.3142175348e+02
        real(wp), parameter :: a53   = 0.5152903406e+01  !jun21 fit
        real(wp), parameter :: g11   = -.4765732725e+02
        real(wp), parameter :: g12   = 0.3648933563e+01
        real(wp), parameter :: g13   = -.7141145244e-01
        real(wp), parameter :: g21   = 0.1002349176e-01
        real(wp), parameter :: g22   = 0.9989856329e-02
        real(wp), parameter :: g23   = -.4161953634e-02
        real(wp), parameter :: g24   = 0.9075807910e-03
        real(wp), parameter :: g31   = -.2693628729e+00
        real(wp), parameter :: g32   = -.1399065763e-01
        real(wp), parameter :: g41   = -.1417634346e-01
        real(wp), parameter :: g42   = -.4870024792e-03
        real(wp), parameter :: g43   = 0.1312231847e+00
        real(wp), parameter :: g44   = -.4409850519e-01
        real(wp), parameter :: g51   = 0.5382970863e+02
        real(wp), parameter :: g52   = 0.4587102824e+01
        real(wp), parameter :: g53   = 0.1768550515e+01  !jun21 fit
        real(wp), parameter :: c11   = 0.1860299931e+04
        real(wp), parameter :: c12   = -.6134458037e+03
        real(wp), parameter :: c13   = 0.7337207161e+02
        real(wp), parameter :: c14   = -.2676717625e+04
        real(wp), parameter :: c15   = 0.1344099415e+04
        real(wp), parameter :: c21   = 0.1538913137e+03
        real(wp), parameter :: c22   = 0.4348007369e+02
        real(wp), parameter :: c23   = 0.1719720677e+03
        real(wp), parameter :: c24   = 0.2115963042e+03
        real(wp), parameter :: c31   = -.7026089414e+02
        real(wp), parameter :: c32   = -.1300938992e+03
        real(wp), parameter :: c41   = 0.1310273564e+01
        real(wp), parameter :: c42   = -.6175149574e+00
        real(wp), parameter :: c43   = -.2679089358e+02
        real(wp), parameter :: c44   = 0.5577477171e+01
        real(wp), parameter :: c51   = -.3543353539e+04
        real(wp), parameter :: c52   = -.3740709591e+03
        real(wp), parameter :: c53   = 0.7979303144e+02
        real(wp), parameter :: c61   = -.1104230585e+04
        real(wp), parameter :: c62   = 0.4603572025e+04
        real(wp), parameter :: c63   = -.5593496634e+04
        real(wp), parameter :: c71   = -.1069406434e+02
        real(wp), parameter :: c72   = 0.1021807153e+01
        real(wp), parameter :: c73   = 0.6669828341e-01
        real(wp), parameter :: c74   = 0.4168542348e+02
        real(wp), parameter :: c75   = 0.1751608567e+02
        real(wp), parameter :: c81   = 0.9486883238e+02
        real(wp), parameter :: c82   = -.1519334221e+02
        real(wp), parameter :: c83   = 0.4024697252e+04
        real(wp), parameter :: c84   = -.2225159395e+02  !jun21
        real(wp), parameter :: d11   = 0.4203543357e+03
        real(wp), parameter :: d12   = -.4922474096e+02
        real(wp), parameter :: d13   = 0.3362942544e+00
        real(wp), parameter :: d14   = -.3827423082e+03
        real(wp), parameter :: d15   = 0.1746726001e+03
        real(wp), parameter :: d21   = 0.1699995737e-01
        real(wp), parameter :: d22   = 0.1513036778e-01
        real(wp), parameter :: d23   = 0.2659119354e-01
        real(wp), parameter :: d24   = -.5760387483e-02
        real(wp), parameter :: d31   = 0.1020622621e+02
        real(wp), parameter :: d32   = 0.1050536271e-01
        real(wp), parameter :: d41   = 0.6836172780e+00
        real(wp), parameter :: d42   = -.1627858240e+00
        real(wp), parameter :: d43   = -.6925485045e+01
        real(wp), parameter :: d44   = 0.1632567385e+01
        real(wp), parameter :: d51   = 0.1083595009e+04
        real(wp), parameter :: d52   = 0.4641431791e+01
        real(wp), parameter :: d53   = -.8233144461e+00
        real(wp), parameter :: d61   = -.6157225942e+02
        real(wp), parameter :: d62   = 0.3094361471e+03
        real(wp), parameter :: d63   = -.3299631143e+03
        real(wp), parameter :: d71   = 0.8866227120e+01
        real(wp), parameter :: d72   = -.1382126854e+01
        real(wp), parameter :: d73   = 0.7620770145e-01
        real(wp), parameter :: d74   = -.5145757859e+02
        real(wp), parameter :: d75   = 0.2046097265e+01
        real(wp), parameter :: d81   = 0.2540775558e+01
        real(wp), parameter :: d82   = -.4889246569e+00
        real(wp), parameter :: d83   = -.1127439280e+04
        real(wp), parameter :: d84   = -.2269932295e+01  !jun21
        cx1      = c51 + c83
        dx1      = d51 + d83
        r1       = rpass(1)
        r2       = rpass(2)
        r3       = rpass(3)
        t1       = r1*r1 -r2*r2 -r3*r3
        t2       = r2*r2 -r3*r3 -r1*r1
        t3       = r3*r3 -r1*r1 -r2*r2
        c1       = t1 / (-2.e0_wp*r2*r3)
        c2       = t2 / (-2.e0_wp*r3*r1)
        c3       = t3 / (-2.e0_wp*r1*r2)
        sum      = c1 + c2 + c3
        b1a      = 1.e0_wp - sum
        c1cube   = c1*c1*c1
        c2cube   = c2*c2*c2
        c3cube   = c3*c3*c3
        cos3t1   = 4.e0_wp*c1cube - 3.e0_wp*c1
        cos3t2   = 4.e0_wp*c2cube - 3.e0_wp*c2
        cos3t3   = 4.e0_wp*c3cube - 3.e0_wp*c3
        sumb     = cos3t1 + cos3t2 + cos3t3
        b1b      = 1.e0_wp - ( z58*sumb + z38*sum )
        dc1dr1   = -r1/(r2*r3)
        dc2dr2   = -r2/(r1*r3)
        dc3dr3   = -r3/(r1*r2)
        dc1dr2   = ( t1/(r2*r2) + 2.e0_wp )/(2.e0_wp*r3)
        dc1dr3   = ( t1/(r3*r3) + 2.e0_wp )/(2.e0_wp*r2)
        dc2dr1   = ( t2/(r1*r1) + 2.e0_wp )/(2.e0_wp*r3)
        dc2dr3   = ( t2/(r3*r3) + 2.e0_wp )/(2.e0_wp*r1)
        dc3dr1   = ( t3/(r1*r1) + 2.e0_wp )/(2.e0_wp*r2)
        dc3dr2   = ( t3/(r2*r2) + 2.e0_wp )/(2.e0_wp*r1)
        db1a(1)  = -1.e0_wp*( dc1dr1 + dc2dr1 + dc3dr1 )
        db1a(2)  = -1.e0_wp*( dc1dr2 + dc2dr2 + dc3dr2 )
        db1a(3)  = -1.e0_wp*( dc1dr3 + dc2dr3 + dc3dr3 )
        d1       = 12.e0_wp*c1*c1 - 3.e0_wp
        d2       = 12.e0_wp*c2*c2 - 3.e0_wp
        d3       = 12.e0_wp*c3*c3 - 3.e0_wp
        db1b(1)  = -z58*( d1*dc1dr1 + d2*dc2dr1 + d3*dc3dr1 ) -z38*(    dc1dr1 +    dc2dr1 +    dc3dr1 )
        db1b(2)  = -z58*( d1*dc1dr2 + d2*dc2dr2 + d3*dc3dr2 ) -z38*(    dc1dr2 +    dc2dr2 +    dc3dr2 )
        db1b(3)  = -z58*( d1*dc1dr3 + d2*dc2dr3 + d3*dc3dr3 ) -z38*(    dc1dr3 +    dc2dr3 +    dc3dr3 )
        r        = r1 + r2 + r3
        rsq      = r*r
        b2       = 1.e0_wp/r1 + 1.e0_wp/r2 + 1.e0_wp/r3
        b3       = (r2-r1)*(r2-r1) + (r3-r2)*(r3-r2) + (r1-r3)*(r1-r3)
        eps2     = 1.e-12_wp
        b3b      = sqrt( b3 + eps2 )
        db2(1)   = -1.e0_wp/(r1*r1)
        db2(2)   = -1.e0_wp/(r2*r2)
        db2(3)   = -1.e0_wp/(r3*r3)
        db3(1)   = 4.e0_wp*r1 -2.e0_wp*r2 -2.e0_wp*r3
        db3(2)   = 4.e0_wp*r2 -2.e0_wp*r3 -2.e0_wp*r1
        db3(3)   = 4.e0_wp*r3 -2.e0_wp*r1 -2.e0_wp*r2
        db3(4)   = 0.5e0_wp*db3(1)/b3b
        db3(5)   = 0.5e0_wp*db3(2)/b3b
        db3(6)   = 0.5e0_wp*db3(3)/b3b
        exp1     = exp(-beta1*r)
        exp2     = exp(-beta2*rsq)
        exp7     = exp(-beta3*r)
        dexp1    = -beta1*exp1
        dexp2    = -2.e0_wp*beta2*r*exp2
        dexp7    = -beta3*exp7
        b1       = b1a
        b12      = b1*b1
        b13      = b12*b1
        b14      = b13*b1
        b15      = b14*b1
        asum     = a11 +  a12*r +  a13*rsq
        bsum     = a21*b12 +  a22*b13 +  a23*b14 +  a24*b15
        csum     = a31*b1*exp1 +  a32*b12*exp2
        dsum1    = a41*exp1 +  a42*exp2
        dsum2    = a43*exp1 +  a44*exp2
        fsum     = a51 +  a52*r +  a53*rsq
        vb(1, 1) = b1*asum*exp1
        vb(1, 2) = bsum * exp2
        vb(1, 3) = b2*csum
        vb(1, 4) = b1*b3*dsum1 + b1*b3b*dsum2
        vb(1, 5) = b1* fsum  *exp7
        vbnd(1)  = vb(1, 1)+vb(1, 2)+vb(1, 3)+vb(1, 4)+vb(1, 5)
        dasum    = a12+2.e0_wp*a13*r
        dbsum    = 2.e0_wp* a21*b1  + 3.e0_wp* a22*b12 + 4.e0_wp* a23*b13 + 5.e0_wp* a24*b14
        ddsum1   = a41*dexp1 +  a42*dexp2
        ddsum2   = a43*dexp1 +  a44*dexp2
        dfsum    = a52 +2.e0_wp*a53*r
        dvb1(1)  = db1a(1)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
        dvb1(2)  = db1a(2)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
        dvb1(3)  = db1a(3)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
        dvb2(1)  = dbsum*db1a(1)*exp2 + bsum*dexp2
        dvb2(2)  = dbsum*db1a(2)*exp2 + bsum*dexp2
        dvb2(3)  = dbsum*db1a(3)*exp2 + bsum*dexp2
        dvb3(1)  = db2(1)*csum + b2*(  a31*db1a(1)*exp1 +  a31*b1*dexp1 +2.e0_wp* a32*b1*db1a(1)*exp2 +  a32*b12*dexp2)
        dvb3(2)  = db2(2)*csum + b2*(  a31*db1a(2)*exp1 +  a31*b1*dexp1 +2.e0_wp* a32*b1*db1a(2)*exp2 +  a32*b12*dexp2)
        dvb3(3)  = db2(3)*csum + b2*(  a31*db1a(3)*exp1 +  a31*b1*dexp1 +2.e0_wp* a32*b1*db1a(3)*exp2 +  a32*b12*dexp2)
        dvb4(1)   = db1a(1)*b3 *dsum1+  b1*db3(1)*dsum1 + b1*b3 *ddsum1 + db1a(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1* &
            & b3b* ddsum2
        dvb4(2)   = db1a(2)*b3 *dsum1+  b1*db3(2)*dsum1 + b1*b3 *ddsum1 + db1a(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1* &
            & b3b* ddsum2
        dvb4(3)   = db1a(3)*b3 *dsum1+  b1*db3(3)*dsum1 + b1*b3 *ddsum1 + db1a(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1* &
            & b3b* ddsum2
        dvb5(1)   = db1a(1)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
        dvb5(2)   = db1a(2)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
        dvb5(3)   = db1a(3)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
        dvbnda(1) = dvb1(1)+dvb2(1)+dvb3(1)+dvb4(1)+dvb5(1)
        dvbnda(2) = dvb1(2)+dvb2(2)+dvb3(2)+dvb4(2)+dvb5(2)
        dvbnda(3) = dvb1(3)+dvb2(3)+dvb3(3)+dvb4(3)+dvb5(3)
        vbnda     = vbnd(1)
        b1        = b1b
        b12       = b1*b1
        b13       = b12*b1
        b14       = b13*b1
        b15       = b14*b1
        asum      = g11 +  g12*r +  g13*rsq
        bsum      = g21*b12 +  g22*b13 +  g23*b14 +  g24*b15
        csum      = g31*b1*exp1 +  g32*b12*exp2
        dsum1     = g41*exp1 +  g42*exp2
        dsum2     = g43*exp1 +  g44*exp2
        fsum      = g51 +  g52*r +  g53*rsq
        vb(2, 1)  = b1*asum*exp1
        vb(2, 2)  = bsum * exp2
        vb(2, 3)  = b2*csum
        vb(2, 4)  = b1*b3*dsum1 + b1*b3b*dsum2
        vb(2, 5)  = b1* fsum  *exp7
        vbnd(2)   = vb(2, 1)+vb(2, 2)+vb(2, 3)+vb(2, 4)+vb(2, 5)
        dasum     = g12+2.e0_wp*g13*r
        dbsum     = 2.e0_wp* g21*b1  + 3.e0_wp* g22*b12 +4.e0_wp* g23*b13 + 5.e0_wp* g24*b14
        ddsum1    = g41*dexp1 +  g42*dexp2
        ddsum2    = g43*dexp1 +  g44*dexp2
        dfsum     = g52 +2.e0_wp*g53*r
        dvb1(1)   = db1b(1)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
        dvb1(2)   = db1b(2)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
        dvb1(3)   = db1b(3)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
        dvb2(1)   = dbsum*db1b(1)*exp2 + bsum*dexp2
        dvb2(2)   = dbsum*db1b(2)*exp2 + bsum*dexp2
        dvb2(3)   = dbsum*db1b(3)*exp2 + bsum*dexp2
        dvb3(1)   = db2(1)*csum + b2*(  g31*db1b(1)*exp1 +  g31*b1*dexp1 +2.e0_wp* g32*b1*db1b(1)*exp2 +  g32*b12*dexp2)
        dvb3(2)   = db2(2)*csum + b2*(  g31*db1b(2)*exp1 +  g31*b1*dexp1 +2.e0_wp* g32*b1*db1b(2)*exp2 +  g32*b12*dexp2)
        dvb3(3)   = db2(3)*csum + b2*(  g31*db1b(3)*exp1 +  g31*b1*dexp1 +2.e0_wp* g32*b1*db1b(3)*exp2 +  g32*b12*dexp2)
        dvb4(1)   = db1b(1)*b3 *dsum1+  b1*db3(1)*dsum1 + b1*b3 *ddsum1 + db1b(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1* &
            & b3b* ddsum2
        dvb4(2)   = db1b(2)*b3 *dsum1+  b1*db3(2)*dsum1 + b1*b3 *ddsum1 + db1b(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1* &
            & b3b* ddsum2
        dvb4(3)   = db1b(3)*b3 *dsum1+  b1*db3(3)*dsum1 + b1*b3 *ddsum1 + db1b(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1* &
            & b3b* ddsum2
        dvb5(1)   = db1b(1)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
        dvb5(2)   = db1b(2)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
        dvb5(3)   = db1b(3)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
        dvbndb(1) = dvb1(1)+dvb2(1)+dvb3(1)+dvb4(1)+dvb5(1)
        dvbndb(2) = dvb1(2)+dvb2(2)+dvb3(2)+dvb4(2)+dvb5(2)
        dvbndb(3) = dvb1(3)+dvb2(3)+dvb3(3)+dvb4(3)+dvb5(3)
        vbndb     = vbnd(2)
        if( icompc .eq. 0 ) return
        sumt      = t(1) + t(2) + t(3)
        rcu       = rsq*r
        p         = r1*r2*r3
        psq       = p*p
        pcu       = psq*p
        dp(1)     = r2*r3
        dp(2)     = r3*r1
        dp(3)     = r1*r2
        cbnds(1)  = 0.e0_wp
        cbnds(2)  = 0.e0_wp
        dcbnda(1) = 0.e0_wp
        dcbnda(2) = 0.e0_wp
        dcbnda(3) = 0.e0_wp
        dcbndb(1) = 0.e0_wp
        dcbndb(2) = 0.e0_wp
        dcbndb(3) = 0.e0_wp
        exp7      = exp( -beta2*rcu )
        dexp7     = -3.e0_wp*beta2*rsq*exp7
        b1        = b1a
        b12       = b1*b1
        b13       = b12*b1
        b14       = b13*b1
        b15       = b14*b1
        asum      = c11 + c12*r + c13*rsq + c14/r + c15/rsq
        bsum      = c21*b12 + c22*b13 + c23*b14 + c24*b15
        csum      = c31*b1 *exp1  + c32*b12*exp2
        dsum1     = c41*exp1 + c42*exp2
        dsum2     = c43*exp1 + c44*exp2
        fsum      = cx1 + c52*r + c53*rsq
        gsum      = c61 + c62/r + c63/rsq
        aasum     = c71 + c72*p + c73*psq + c74/p + c75/psq
        ffsum     = c81*p + c82*psq + c84/psq
        dasum     = c12 + 2.e0_wp*c13*r -c14/rsq -2.e0_wp*c15/rcu
        dbsum     = 2.e0_wp*c21*b1 +3.e0_wp*c22*b12 +4.e0_wp*c23*b13 +5.e0_wp*c24*b14
        ddsum1    = c41*dexp1 + c42*dexp2
        ddsum2    = c43*dexp1 + c44*dexp2
        dfsum     = c52 + 2.e0_wp*c53*r
        dgsum     = -c62/rsq - 2.e0_wp*c63/rcu
        daasum    = c72 +2.e0_wp*c73*p -c74/psq -2.e0_wp*c75/pcu
        dffsum    = c81 + 2.e0_wp*c82*p -2.e0_wp*c84/pcu
        cb(1,1)   = b1*asum * exp1  /p
        cb(1,2)   = bsum * exp2
        cb(1,3)   = b2*csum
        cb(1,4)   = b1*b3*dsum1 + b1*b3b*dsum2
        cb(1,5)   = b1 * fsum * exp7 /p
        cb(1,6)   = b1 * gsum * exp7
        cb(1,7)   = b1 * aasum * exp2
        cb(1,8)   = b1 * ffsum * exp7
        cbnds(1)  = cb(1,1)+cb(1,2)+cb(1,3)+cb(1,4)+cb(1,5) +cb(1,6)+cb(1,7)+cb(1,8)
        dcb1(1)   = db1a(1)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(1)/psq
        dcb1(2)   = db1a(2)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(2)/psq
        dcb1(3)   = db1a(3)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(3)/psq
        dcb2(1)   = dbsum*db1a(1)*exp2 + bsum*dexp2
        dcb2(2)   = dbsum*db1a(2)*exp2 + bsum*dexp2
        dcb2(3)   = dbsum*db1a(3)*exp2 + bsum*dexp2
        dcb3(1)   = db2(1)*csum + b2*( c31*db1a(1)*exp1 +c31*b1*dexp1 +c32 *2.e0_wp*b1*db1a(1)*exp2 + c32 *b12*dexp2 )
        dcb3(2)   = db2(2)*csum + b2*( c31 *db1a(2)*exp1 +c31 *b1*dexp1 +c32 *2.e0_wp*b1*db1a(2)*exp2 + c32 *b12*dexp2 )
        dcb3(3)   = db2(3)*csum + b2*( c31 *db1a(3)*exp1 +c31 *b1*dexp1 +c32 *2.e0_wp*b1*db1a(3)*exp2 + c32 *b12*dexp2 )
        dcb4(1)     = db1a(1)*b3 *dsum1 + b1*db3(1)*dsum1 + b1*b3 *ddsum1 + db1a(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1* &
            & b3b* ddsum2
        dcb4(2)     = db1a(2)*b3 *dsum1 + b1*db3(2)*dsum1 + b1*b3 *ddsum1 + db1a(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1* &
            & b3b* ddsum2
        dcb4(3)     = db1a(3)*b3 *dsum1 + b1*db3(3)*dsum1 + b1*b3 *ddsum1 + db1a(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1* &
            & b3b* ddsum2
        dcb5(1)      = db1a(1)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(1)/psq
        dcb5(2)      = db1a(2)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(2)/psq
        dcb5(3)      = db1a(3)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(3)/psq
        dcb6(1)      = db1a(1)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
        dcb6(2)      = db1a(2)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
        dcb6(3)      = db1a(3)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
        dcb7(1)      = db1a(1)*aasum*exp2 + b1*daasum*dp(1)*exp2 + b1* aasum*dexp2
        dcb7(2)      = db1a(2)*aasum*exp2 + b1*daasum*dp(2)*exp2 + b1* aasum*dexp2
        dcb7(3)      = db1a(3)*aasum*exp2 + b1*daasum*dp(3)*exp2 + b1* aasum*dexp2
        dcb8(1)      = db1a(1)*ffsum*exp7  +b1*dffsum*dp(1)*exp7 +b1* ffsum*dexp7
        dcb8(2)      = db1a(2)*ffsum*exp7  +b1*dffsum*dp(2)*exp7 +b1* ffsum*dexp7
        dcb8(3)      = db1a(3)*ffsum*exp7  +b1*dffsum*dp(3)*exp7 +b1* ffsum*dexp7
        dcbnds(1, 1) = dcb1(1) +dcb2(1) +dcb3(1) +dcb4(1) +dcb5(1) +dcb6(1) +dcb7(1) +dcb8(1)
        dcbnds(1, 2) = dcb1(2) +dcb2(2) +dcb3(2) +dcb4(2) +dcb5(2) +dcb6(2) +dcb7(2) +dcb8(2)
        dcbnds(1, 3) = dcb1(3) +dcb2(3) +dcb3(3) +dcb4(3) +dcb5(3) +dcb6(3) +dcb7(3) +dcb8(3)
        b1           = b1b
        b12          = b1*b1
        b13          = b12*b1
        b14          = b13*b1
        b15          = b14*b1
        asum         = d11 + d12*r + d13*rsq + d14/r + d15/rsq
        bsum         = d21*b12 + d22*b13 + d23*b14 + d24*b15
        csum         = d31*b1 *exp1  + d32*b12*exp2
        dsum1        = d41*exp1 + d42*exp2
        dsum2        = d43*exp1 + d44*exp2
        fsum         = dx1 + d52*r + d53*rsq
        gsum         = d61 + d62/r + d63/rsq
        aasum        = d71 + d72*p + d73*psq + d74/p + d75/psq
        ffsum        = d81*p + d82*psq + d84/psq
        dasum        = d12 + 2.e0_wp*d13*r -d14/rsq -2.e0_wp*d15/rcu
        dbsum        = 2.e0_wp*d21*b1 +3.e0_wp*d22*b12 +4.e0_wp*d23*b13 +5.e0_wp*d24*b14
        ddsum1       = d41*dexp1 + d42*dexp2
        ddsum2       = d43*dexp1 + d44*dexp2
        dfsum        = d52 + 2.e0_wp*d53*r
        dgsum        = -d62/rsq - 2.e0_wp*d63/rcu
        daasum       = d72 +2.e0_wp*d73*p -d74/psq -2.e0_wp*d75/pcu
        dffsum       = d81 + 2.e0_wp*d82*p -2.e0_wp*d84/pcu
        cb(2,1)      = b1*asum * exp1  /p
        cb(2,2)      = bsum * exp2
        cb(2,3)      = b2*csum
        cb(2,4)      = b1*b3*dsum1 + b1*b3b*dsum2
        cb(2,5)      = b1 * fsum * exp7 /p
        cb(2,6)      = b1 * gsum * exp7
        cb(2,7)      = b1 * aasum * exp2
        cb(2,8)      = b1 * ffsum * exp7
        cbnds(2)     = cb(2,1)+cb(2,2)+cb(2,3)+cb(2,4)+cb(2,5) +cb(2,6)+cb(2,7)+cb(2,8)
        dcb1(1)      = db1b(1)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(1)/psq
        dcb1(2)      = db1b(2)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(2)/psq
        dcb1(3)      = db1b(3)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(3)/psq
        dcb2(1)      = dbsum*db1b(1)*exp2 + bsum*dexp2
        dcb2(2)      = dbsum*db1b(2)*exp2 + bsum*dexp2
        dcb2(3)      = dbsum*db1b(3)*exp2 + bsum*dexp2
        dcb3(1)      = db2(1)*csum + b2*( d31 *db1b(1)*exp1 +d31 *b1*dexp1 +d32 *2.e0_wp*b1*db1b(1)*exp2 + d32 *b12* &
            & dexp2 )
        dcb3(2)      = db2(2)*csum + b2*( d31 *db1b(2)*exp1 +d31 *b1*dexp1 +d32 *2.e0_wp*b1*db1b(2)*exp2 + d32 *b12* &
            & dexp2 )
        dcb3(3)      = db2(3)*csum + b2*( d31 *db1b(3)*exp1 +d31 *b1*dexp1 +d32 *2.e0_wp*b1*db1b(3)*exp2 + d32 *b12* &
            & dexp2 )
        dcb4(1)     = db1b(1)*b3 *dsum1 + b1*db3(1)*dsum1 + b1*b3 *ddsum1 + db1b(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1* &
            & b3b* ddsum2
        dcb4(2)     = db1b(2)*b3 *dsum1 + b1*db3(2)*dsum1 + b1*b3 *ddsum1 + db1b(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1* &
            & b3b* ddsum2
        dcb4(3)     = db1b(3)*b3 *dsum1 + b1*db3(3)*dsum1 + b1*b3 *ddsum1 + db1b(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1* &
            & b3b* ddsum2
        dcb5(1)      = db1b(1)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(1)/psq
        dcb5(2)      = db1b(2)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(2)/psq
        dcb5(3)      = db1b(3)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(3)/psq
        dcb6(1)      = db1b(1)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
        dcb6(2)      = db1b(2)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
        dcb6(3)      = db1b(3)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
        dcb7(1)      = db1b(1)*aasum*exp2 + b1*daasum*dp(1)*exp2 + b1* aasum*dexp2
        dcb7(2)      = db1b(2)*aasum*exp2 + b1*daasum*dp(2)*exp2 + b1* aasum*dexp2
        dcb7(3)      = db1b(3)*aasum*exp2 + b1*daasum*dp(3)*exp2 + b1* aasum*dexp2
        dcb8(1)      = db1b(1)*ffsum*exp7  +b1*dffsum*dp(1)*exp7 +b1* ffsum*dexp7
        dcb8(2)      = db1b(2)*ffsum*exp7  +b1*dffsum*dp(2)*exp7 +b1* ffsum*dexp7
        dcb8(3)      = db1b(3)*ffsum*exp7  +b1*dffsum*dp(3)*exp7 +b1* ffsum*dexp7
        dcbnds(2, 1) = dcb1(1) +dcb2(1) +dcb3(1) +dcb4(1) +dcb5(1) +dcb6(1) +dcb7(1) +dcb8(1)
        dcbnds(2, 2) = dcb1(2) +dcb2(2) +dcb3(2) +dcb4(2) +dcb5(2) +dcb6(2) +dcb7(2) +dcb8(2)
        dcbnds(2, 3) = dcb1(3) +dcb2(3) +dcb3(3) +dcb4(3) +dcb5(3) +dcb6(3) +dcb7(3) +dcb8(3)
        dcbnda(1)    = dt(1) * cbnds(1)    + sumt * dcbnds(1, 1)
        dcbndb(1)    = dt(1) * cbnds(2)    + sumt * dcbnds(2, 1)
        dcbnda(2)    = dt(2) * cbnds(1)    + sumt * dcbnds(1, 2)
        dcbndb(2)    = dt(2) * cbnds(2)    + sumt * dcbnds(2, 2)
        dcbnda(3)    = dt(3) * cbnds(1)    + sumt * dcbnds(1, 3)
        dcbndb(3)    = dt(3) * cbnds(2)    + sumt * dcbnds(2, 3)
        cbnda        = sumt*cbnds(1)
        cbndb        = sumt*cbnds(2)
    end
    !> @brief Evaluate chgeom.
    subroutine chgeom( r, ivalid )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp)            :: r(3), rhi, rlo, rmid
        integer             :: ivalid
        real(wp), parameter :: derr   = 1.e-5_wp
        integer, parameter  :: istop1 = 1
        integer, parameter  :: istop2 = 1
        ivalid = 1
        rhi    = max(r(1), r(2))
        if( r(3) .ge. rhi )then
            rmid = rhi
            rhi  = r(3)
        else
            rmid = max( r(3), min(r(1), r(2)) )
        endif
        rlo = min(r(1), r(2), r(3))
        if( rlo+rmid+derr .lt. rhi )then
            write(6, *) 'Warning: invalid geometry ', r
            if( istop1 .gt. 0 ) stop ' STOP -- geometry error '
        endif
        if( rlo .lt. 0.2e0_wp )then
            write(6, *) 'Warning: invalid geometry ', r
            if( istop2 .gt. 0 ) stop ' STOP -- geometry error '
        endif
    end
end module chempotpy_h3_gen_bkmp3_1996_pes

!> @brief Expose H3_GEN_BKMP3_1996 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_bkmp3_1996_pes, only: evaluate
    use chempotpy_kinds,                 only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
