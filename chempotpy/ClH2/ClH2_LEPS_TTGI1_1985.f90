!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_clh2_leps_ttgi1_1985_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
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
        double precision              :: dq(3), q(3), r2(3), tx(9), dedx(9), drdx(3,9), v, xter
        integer                       :: iatom, idir, istate, j
        !initialize
        v = 0.d0; g = 0.d0; d = 0.d0
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is ClHH
        call legacy_triatomic_distances_to_bohr(x, q)
        dq = 0.d0
        call pot(q, igrad, xter, dq)
        xter = xter*legacy_hartree_to_electronvolt
        dq   = dq*51.422067
        r2   = q*legacy_bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(dq, drdx, dedx)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = xter
            enddo
        else if( igrad == 1 )then
            do istate=1, nstates
                p(istate) = xter
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    j               = (iatom-1)*3+idir
                    g(1,iatom,idir) = dedx(j)
                enddo
            enddo
        else if( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
    end subroutine evaluate
end module chempotpy_clh2_leps_ttgi1_1985_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_clh2_leps_ttgi1_1985_pes, only: evaluate
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

!> @brief Evaluate EvdRdX for the ClH2_LEPS_TTGI1_1985 potential.
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

!
module chempotpy_clh2_leps_ttgi1_1985_pot_legacy_state_bf7e4ade
    implicit none
    private
    public :: evaluate_internal
    contains
    !> @brief Evaluate dispatch.
    subroutine evaluate_internal( q, nder, xter, dq )
        !
        !   System:          ClH2
        !   Functional form: Isaacson and Muckerman's diatomics-in-molecules ver
        !   Common name:     ClH2 DIMA
        !   Reference:       unpublished
        !   Cross Reference: S. C. Tucker, D. G. Truhlar,
        !                    G. C. Garrett, and A. D. Isaacson
        !                    J. Chem. Phys. 82, 4102 (1985)
        !
        !   Calling Sequence:
        !               must be called once before any calls to POT
        !      POT    - driver for the evaluation of the energy and the derivati
        !               of the energy with respect to the coordinates for a give
        !               geometric configuration
        !
        !   Units:
        !      energies    - hartrees
        !      coordinates - bohr
        !      derivatives - hartrees/bohr
        !
        !   Surfaces:
        !      ground electronic state
        !
        !   Zero of energy:
        !      The classical potential energy is set equal to zero for the Cl
        !      infinitely far from the H2 diatomic and R(H2) set equal to the
        !      H2 equilibrium diatomic value.
        !
        !   Parameters:
        !      Set in the BLOCK DATA subprogram PTPARM
        !
        !   Coordinates:
        !      Internal, Definition: R(1) = R(Cl-H)
        !                            R(2) = R(H-H)
        !                            R(3) = R(Cl-H)
        !
        !   Common Blocks (used between the calling program and this potential):
        !      /PT1CM/ R(3), ENERGY, DEDR(3)
        !        passes the coordinates, ground state electronic energy, and
        !        derivatives of the ground electronic state energy with respect
        !        to the coordinates.
        !      /PT2CM/ NDER, NDUM, IDEBUG, IDUM(19)
        !        passes the control flags where
        !        NDER   = 0 => no derivatives are computed
        !               = 1 => derivatives of the energy for the ground electron
        !                      state with respect to the coordinates are compute
        !        NDUM   - not used
        !        IDEBUG = 1 => print intermediate values to the file linked to
        !                      FORTRAN unit IPTPRT
        !        IDUM   - not used
        !      /PT4CM/ IPTPRT
        !        passes the FORTRAN unit number used for potential output
        !      /PT5CM/ EASYAB, EASYBC, EASYAC
        !        passes the energy in the three asymptotic valleys for an A + BC
        !        The energy in the AB valley, EASYAB, is equal to the energy of
        !        C atom "infinitely" far from the AB diatomic and R(AB) set equa
        !        Re(AB), the equilibrium bond length for the AB diatomic.
        !        In this potential the AB valley represents H infinitely far fro
        !        the ClH diatomic and R(ClH) equal to Re(ClH).  Similarly, the t
        !        EASYBC and EASYAC represent the energies in the H2 and the othe
        !        valleys, respectively.
        !
        !   Default Parameter Values:
        !      Variable      Default value
        !      NDER             1
        !      IDEBUG           0
        !      IPTPRT           6
        !
        !
        implicit none
        double precision, intent(in)  :: q(3)
        integer, intent(in)           :: nder
        double precision, intent(out) :: xter, dq(3)
        double precision              :: v1h2(79), v2h2(87), v3h2(126), v4h2(126)
        double precision              :: r1h2(79), r2h2(87), r3h2(126)
        double precision              :: s2ss(79), s3ss(79), delss(79)
        double precision              :: s2ts(87), s3ts(87), delts(87)
        double precision              :: s2dg(126), s3dg(126), deldg(126)
        double precision              :: s2du(126), s3du(126), deldu(126)
        double precision              :: quar, rt3
        dimension c(3)
        double precision :: arg, c, cs, ct, ct2, dd, ddgab, ddsxa, ddsxb, dduab, dgab, dh, dh1, dh11, dh12, dh13, &
            & dh14, dh15, dh16, dh2, dh22, dh24, dh25, dh26, dh33, dh34, dh44, dh45, dh46, dh55, dh56, dh66, dnumb, &
            & dpux, dsgh, dspxa, dspxb, dss11a, dss11b, dss12a, dss12b, dss22a, dss22b, dssab, dsxa, dsxb, dtdr1, &
            & dtdr2, dtdr3, dtpxa, dtpxb, dtsab, dtsxa, dtsxb, duab, easyab, easyac, easybc, eht, en, factr, for, &
            & forth, gvec, h1, h2, half, one, r, rt2i, six, spxa, spxb, sqrt3, ss11a, ss11b, ss12a, ss12b, ss22a, &
            & ss22b, ssab, ssgh, ssgx, st, st2, test, th, tpxa, tpxb, tre, tsab, tsxa, tsxb, two, vec, w, xtra, zero
        integer :: i, idebug, idum, iflag, indx, iptprt, ix, j, k, ndim, ndum, nroot
        dimension h1(2,2), h2(2,2), dh1(2,2), dh2(2,2)
        dimension th(6,6), vec(6,6), w(6), en(6), gvec(6), dh(6,6), r(3), dd(3)
        dimension dh11(3), dh12(3), dh13(3), dh14(3), dh15(3), dh16(3), dh22(3)
        dimension dh24(3), dh25(3), dh26(3), dh33(3), dh34(3), dh44(3), dh45(3)
        dimension dh46(3), dh55(3), dh56(3), dh66(3)
        character(len=80) :: title(2)
        dimension indx(3)
        equivalence (h1(1,1), ss11a), (h1(1,2), ss12a), (h1(2,2), ss22a)
        equivalence (h2(1,1), ss11b), (h2(1,2), ss12b), (h2(2,2), ss22b)
        equivalence (dh1(1,1), dss11a), (dh1(1,2), dss12a), (dh1(2,2), dss22a)
        equivalence (dh2(1,1), dss11b), (dh2(1,2), dss12b), (dh2(2,2), dss22b)
        equivalence (th(1,1), dh(1,1))
        !   Initialize the array INDX which specifies the order of the
        !   input bond distances.  INDX=1,2,3 for Cl + H2 and INDX=132 for H + C
        !
        !   Echo the name of the potential to the file linked to FORTRAN unit IP
        !WRITE (IPTPRT, 600) C
        !WRITE (IPTPRT, 610) TITLE
        !
        !
        !     SET UP SPLINE FITS TO H2 CURVES
        arg    = 0.d0
        c      = 0.d0
        cs     = 0.d0
        ct     = 0.d0
        ct2    = 0.d0
        dd     = 0.d0
        ddgab  = 0.d0
        ddsxa  = 0.d0
        ddsxb  = 0.d0
        dduab  = 0.d0
        deldg  = 0.d0
        deldu  = 0.d0
        delss  = 0.d0
        delts  = 0.d0
        dgab   = 0.d0
        dh     = 0.d0
        dh1    = 0.d0
        dh11   = 0.d0
        dh12   = 0.d0
        dh13   = 0.d0
        dh14   = 0.d0
        dh15   = 0.d0
        dh16   = 0.d0
        dh2    = 0.d0
        dh22   = 0.d0
        dh24   = 0.d0
        dh25   = 0.d0
        dh26   = 0.d0
        dh33   = 0.d0
        dh34   = 0.d0
        dh44   = 0.d0
        dh45   = 0.d0
        dh46   = 0.d0
        dh55   = 0.d0
        dh56   = 0.d0
        dh66   = 0.d0
        dnumb  = 0.d0
        dpux   = 0.d0
        dsgh   = 0.d0
        dspxa  = 0.d0
        dspxb  = 0.d0
        dss11a = 0.d0
        dss11b = 0.d0
        dss12a = 0.d0
        dss12b = 0.d0
        dss22a = 0.d0
        dss22b = 0.d0
        dssab  = 0.d0
        dsxa   = 0.d0
        dsxb   = 0.d0
        dtdr1  = 0.d0
        dtdr2  = 0.d0
        dtdr3  = 0.d0
        dtpxa  = 0.d0
        dtpxb  = 0.d0
        dtsab  = 0.d0
        dtsxa  = 0.d0
        dtsxb  = 0.d0
        duab   = 0.d0
        easyab = 0.d0
        easyac = 0.d0
        easybc = 0.d0
        eht    = 0.d0
        en     = 0.d0
        factr  = 0.d0
        for    = 0.d0
        forth  = 0.d0
        gvec   = 0.d0
        h1     = 0.d0
        h2     = 0.d0
        half   = 0.d0
        one    = 0.d0
        quar   = 0.d0
        r      = 0.d0
        r1h2   = 0.d0
        r2h2   = 0.d0
        r3h2   = 0.d0
        rt2i   = 0.d0
        rt3    = 0.d0
        s2dg   = 0.d0
        s2du   = 0.d0
        s2ss   = 0.d0
        s2ts   = 0.d0
        s3dg   = 0.d0
        s3du   = 0.d0
        s3ss   = 0.d0
        s3ts   = 0.d0
        six    = 0.d0
        spxa   = 0.d0
        spxb   = 0.d0
        sqrt3  = 0.d0
        ss11a  = 0.d0
        ss11b  = 0.d0
        ss12a  = 0.d0
        ss12b  = 0.d0
        ss22a  = 0.d0
        ss22b  = 0.d0
        ssab   = 0.d0
        ssgh   = 0.d0
        ssgx   = 0.d0
        st     = 0.d0
        st2    = 0.d0
        test   = 0.d0
        th     = 0.d0
        tpxa   = 0.d0
        tpxb   = 0.d0
        tre    = 0.d0
        tsab   = 0.d0
        tsxa   = 0.d0
        tsxb   = 0.d0
        two    = 0.d0
        v1h2   = 0.d0
        v2h2   = 0.d0
        v3h2   = 0.d0
        v4h2   = 0.d0
        vec    = 0.d0
        w      = 0.d0
        xtra   = 0.d0
        zero   = 0.d0
        iptprt = 6
        idebug = 0
        zero   = 0.d0
        one    = 1.d0
        two    = 2.d0
        tre    = 3.d0
        for    = 4.d0
        six    = 6.d0
        eht    = 8.d0
        half   = 0.5d0
        quar   = 0.25d0
        rt2i   = 7.07106781186548d-1
        rt3    = 1.73205080756888d0
        dpux   = 0.d0
        ssgx   = -0.132774d0
        dsgh   = 0.d0
        ssgh   = 0.5d0
        c      = [ 0.169552d0, 0.1744734d0, 0.169552d0]
        v1h2 = [ -.1202028d0, -.5266270d0, -.7696341d0, -.9220261d0, -1.0200556d0, -1.0836422d0, -1.1245385d0, - &
            & 1.1500562d0, -1.1649342d0, -1.1723459d0, -1.1744744d0, -1.1728537d0, -1.1685799d0, -1.1624570d0, - &
            & 1.1550670d0, -1.1468496d0, -1.1381312d0, -1.1291562d0, -1.1201233d0, -1.1111725d0, -1.1024127d0, - &
            & 1.0939273d0, -1.0857810d0, -1.0780164d0, -1.0706700d0, -1.0637641d0, -1.0573118d0, -1.0513185d0, - &
            & 1.0457832d0, -1.0407003d0, -1.0360578d0, -1.0318402d0, -1.0280272d0, -1.0245978d0, -1.0215297d0, - &
            & 1.0187967d0, -1.0163689d0, -1.0142247d0, -1.0123371d0, -1.0106810d0, -1.0092303d0, -1.0079682d0, - &
            & 1.0068703d0, -1.0059178d0, -1.0050923d0, -1.0043782d0, -1.0037626d0, -1.0032309d0, -1.0027740d0, - &
            & 1.0023800d0, -1.0020423d0, -1.0017521d0, -1.0015030d0, -1.0012899d0, -1.0011069d0, -1.0009498d0, - &
            & 1.0008150d0, -1.0007002d0, -1.0006030d0, -1.0005162d0, -1.0004466d0, -1.0003864d0, -1.0003328d0, - &
            & 1.0002906d0, -1.0002466d0, -1.0002154d0, -1.0001889d0, -1.0001434d0, -1.0001086d0, -1.0000868d0, - &
            & 1.0000682d0, -1.0000528d0, -1.0000404d0, -1.0000314d0, -1.0000185d0, -1.0000121d0, -1.0000091d0, - &
            & 1.0000045d0, -1.0000025d0]
        v2h2 = [ -.6215227d0, -.6757268d0, -.7189640d0, -.7544033d0, -.7841501d0, -.8096095d0, -.8317238d0, - &
            & .8511347d0, -.8682913d0, -.8835186d0, -.8970636d0, -.9091230d0, -.9198593d0, -.9294123d0, -.9379051d0, - &
            & .9454463d0, -.9521346d0, -.9580585d0, -.9632985d0, -.9679271d0, -.9720104d0, -.9756071d0, -.9787717d0, - &
            & .9815517d0, -.9839910d0, -.9861279d0, -.9879977d0, -.9896312d0, -.9910564d0, -.9922979d0, -.9933781d0, - &
            & .9943164d0, -.9951304d0, -.9958355d0, -.9964456d0, -.9969715d0, -.9974252d0, -.9978159d0, -.9981500d0, - &
            & .9984382d0, -.9986849d0, -.9988960d0, -.9990763d0, -.9992292d0, -.9993597d0, -.9994703d0, -.9995635d0, - &
            & .9996432d0, -.9997094d0, -.9997655d0, -.9998125d0, -.9998515d0, -.9998842d0, -.9999112d0, -.9999337d0, - &
            & .9999521d0, -.9999670d0, -.9999795d0, -.9999887d0, -.9999966d0, -1.0000030d0, -1.0000076d0, - &
            & 1.0000117d0, -1.0000143d0, -1.0000167d0, -1.0000180d0, -1.0000191d0, -1.0000195d0, -1.0000197d0, - &
            & 1.0000194d0, -1.0000196d0, -1.0000192d0, -1.0000191d0, -1.0000180d0, -1.0000173d0, -1.0000164d0, - &
            & 1.0000158d0, -1.0000152d0, -1.0000143d0, -1.0000143d0, -1.0000127d0, -1.0000109d0, -1.0000095d0, - &
            & 1.0000076d0, -1.0000067d0, -1.0000043d0, -1.0000025d0]
        v3h2(1:64) = [ 8.0217579d0, 3.0713797d0, 1.4666292d0, .6992459d0, .2650120d0, -.0048180d0, -.1826248d0, - &
            & .3044801d0, -.3902705d0, -.4517863d0, -.4964119d0, -.5289745d0, -.5527406d0, -.5699835d0, -.5823232d0, - &
            & .5909372d0, -.5966963d0, -.6002536d0, -.6021058d0, -.6026342d0, -.6021349d0, -.6008396d0, -.5989309d0, - &
            & .5965536d0, -.5938235d0, -.5908332d0, -.5876573d0, -.5843560d0, -.5809780d0, -.5775629d0, -.5741424d0, - &
            & .5707426d0, -.5673842d0, -.5640840d0, -.5608555d0, -.5577093d0, -.5546536d0, -.5516947d0, -.5488374d0, - &
            & .5460849d0, -.5434394d0, -.5409022d0, -.5384736d0, -.5361531d0, -.5339400d0, -.5318328d0, -.5298296d0, - &
            & .5279281d0, -.5261259d0, -.5244203d0, -.5228082d0, -.5212866d0, -.5198521d0, -.5185016d0, -.5172315d0, - &
            & .5160385d0, -.5149192d0, -.5138701d0, -.5128878d0, -.5119690d0, -.5111105d0, -.5103089d0, -.5095612d0, - &
            & .5088644d0]
        v3h2(65:126) = [ -.5082155d0, -.5076116d0, -.5070501d0, -.5065283d0, -.5060437d0, -.5055940d0, -.5051769d0, - &
            & .5047902d0, -.5044319d0, -.5041001d0, -.5037930d0, -.5035087d0, -.5032458d0, -.5030027d0, -.5027780d0, - &
            & .5025704d0, -.5023785d0, -.5022014d0, -.5020377d0, -.5018866d0, -.5017472d0, -.5016184d0, -.5014996d0, - &
            & .5013900d0, -.5012888d0, -.5011954d0, -.5011093d0, -.5010298d0, -.5009564d0, -.5008887d0, -.5008262d0, - &
            & .5007685d0, -.5007153d0, -.5006661d0, -.5006207d0, -.5005787d0, -.5004122d0, -.5002992d0, -.5002219d0, - &
            & .5001683d0, -.5001306d0, -.5001035d0, -.5000837d0, -.5000690d0, -.5000577d0, -.5000489d0, -.5000420d0, - &
            & .5000364d0, -.5000318d0, -.5000279d0, -.5000247d0, -.5000220d0, -.5000196d0, -.5000176d0, -.5000158d0, - &
            & .5000142d0, -.5000129d0, -.5000117d0, -.5000106d0, -.5000097d0, -.5000089d0, -.5000081d0]
        v4h2(1:64) = [ 9.4993326d0, 4.4973226d0, 2.8272874d0, 1.9892158d0, 1.4831145d0, 1.1423563d0, .8955959d0, &
            & .7072541d0, .5576760d0, .4351864d0, .3324666d0, .2447306d0, .1687347d0, .1022057d0, .0434985d0, - &
            & .0086173d0, -.0550906d0, -.0966755d0, -.1339863d0, -.1675344d0, -.1977530d0, -.2250137d0, -.2496386d0, - &
            & .2719091d0, -.2920720d0, -.3103453d0, -.3269215d0, -.3419716d0, -.3556477d0, -.3680850d0, -.3794043d0, - &
            & .3897134d0, -.3991089d0, -.4076771d0, -.4154957d0, -.4226343d0, -.4291554d0, -.4351155d0, -.4405653d0, - &
            & .4455506d0, -.4501129d0, -.4542896d0, -.4581145d0, -.4616184d0, -.4648290d0, -.4677716d0, -.4704692d0, - &
            & .4729427d0, -.4752111d0, -.4772916d0, -.4792001d0, -.4809509d0, -.4825572d0, -.4840309d0, -.4853831d0, - &
            & .4866238d0, -.4877622d0, -.4888066d0, -.4897648d0, -.4906439d0, -.4914502d0, -.4921898d0, -.4928681d0, - &
            & .4934901d0]
        v4h2(65:126) = [ -.4940603d0, -.4945830d0, -.4950621d0, -.4955011d0, -.4959033d0, -.4962717d0, -.4966090d0, - &
            & .4969178d0, -.4972005d0, -.4974591d0, -.4976956d0, -.4979119d0, -.4981096d0, -.4982903d0, -.4984553d0, - &
            & .4986060d0, -.4987435d0, -.4988690d0, -.4989834d0, -.4990877d0, -.4991827d0, -.4992692d0, -.4993480d0, - &
            & .4994196d0, -.4994847d0, -.4995438d0, -.4995975d0, -.4996462d0, -.4996904d0, -.4997303d0, -.4997665d0, - &
            & .4997992d0, -.4998288d0, -.4998554d0, -.4998794d0, -.4999011d0, -.4999800d0, -.5000244d0, -.5000475d0, - &
            & .5000579d0, -.5000608d0, -.5000595d0, -.5000560d0, -.5000515d0, -.5000467d0, -.5000421d0, -.5000377d0, - &
            & .5000337d0, -.5000301d0, -.5000269d0, -.5000240d0, -.5000216d0, -.5000194d0, -.5000174d0, -.5000157d0, - &
            & .5000142d0, -.5000129d0, -.5000117d0, -.5000106d0, -.5000097d0, -.5000089d0, -.5000081d0]
        forth    = quar
        sqrt3    = rt3
        ndim     = 6
        nroot    = 1
        test     = 0.999999999d0
        indx     = [ 1, 2, 3]
        title(1) = 'Final BNL DIM surface for CL+H2 -- 5 and7 kcal/mole barriers'
        title(2) = 'Spline fits replaced with analytic forms where possible'
        call setspl
        !
        !   Initialize the energy in the three asymptotic valleys
        easyab = c(1)
        easybc = c(2)
        easyac = c(3)
        !
        !
        !
        !
        !     ENTRY POT GIVE ENERGY AND PARTIALS W.R.T 3 DISTANCES
        !
        !     INTERNUCLEAR DISTANCES STORED IN Q(3) AND R(3)
        !     NROOT-TH STATE ENERGY STORED IN XTER (OR TERM FOR DIATOMIC)
        !     THREE PARTIALS STORED IN DD(3)
        !
        !
        !   Check the value of NDER
        !
        !     IAPP=1
        !
        !     SET FLAG FOR EIGENVECTOR OUTPUT
        iflag = sign(1, nroot)
        nroot = abs(nroot)
        !
        !     FOLLOWING COORDINATE DEFINITIONS PROVIDE FOR PROPER SYSTEM
        do 15 ix=1, 3
        15 r(ix)=q(indx(ix))
        if( idebug .eq. 1 ) write(iptprt, '(1X,''DIMIMA, Q='',1P3E13.5,/,9X,'' R='',3E13.5)') q, r
        !  CHECK FOR NEGATIVE R'S
        xter  = 1.d30
        dq(1) = -xter
        dq(2) = -xter
        dq(3) = -xter
        if( r(1) .lt. 0.d0 .or. r(2) .lt. 0.d0 .or. r(3) .lt. 0.d0 ) return
        !
        !     ZERO OUT UNUSED HAMILTONIAN ELEMENTS
        th(2,3) = zero
        th(3,5) = zero
        th(3,6) = zero
        !
        !     COMPUTE SINE AND COSINE OF ANGLE AXB, AND POWERS
        ct = (r(1)*r(1)+r(3)*r(3)-r(2)*r(2))/(two*r(1)*r(3))
        if( ct .gt. -test .and. ct .lt. test ) go to 10
        ct = sign(test, ct)
        10 continue
        arg = one-ct*ct
        st  = sqrt(arg)
        ct2 = ct*ct
        st2 = st*st
        cs  = ct*st
        !     COMPUTE DIATOMIC ENERGIES FOR R(3) ARRAY
        call hcl(r(1), h1, tsxa, spxa, tpxa, dh1, dtsxa, dspxa, dtpxa)
        call hcl(r(3), h2, tsxb, spxb, tpxb, dh2, dtsxb, dspxb, dtpxb)
        call hclm(r(1), dsxa, ddsxa)
        call hclm(r(3), dsxb, ddsxb)
        call h2curv(r(2), ssab, tsab, dgab, duab, dssab, dtsab, ddgab, dduab)
        if( idebug .eq. 1 ) write(iptprt, '(1X,''SS11A,SS11B,DSXA,DSXB,SSAB='',1P5E13.5)') ss11a, ss11b, dsxa, dsxb, ssab
        !
        !     CONSTRUCT DIM HAMILTONIAN IN ORTHOGONAL BASIS
        !
        th(1,1) = ss11a+forth*((ssab+tre*tsab)+ct2*(ss11b+tre*tsxb) +st2*(spxb+tre*tpxb))-dpux-two*dsgh
        th(1,2) = forth*sqrt3*(ssab-tsab-ct2*(ss11b-tsxb)-st2*(spxb-tpxb))
        th(1,3) = ss12a
        th(1,4) = -half*ct*ss12b
        th(1,5) = forth*cs*(ss11b+tre*tsxb-spxb-tre*tpxb)
        th(1,6) = -forth*sqrt3*cs*(ss11b-tsxb-spxb+tpxb)
        th(2,2) = tsxa+forth*((tre*ssab+tsab)+ct2*(tre*ss11b+tsxb) +st2*(tre*spxb+tpxb))-dpux-two*dsgh
        th(2,4) = half*sqrt3*ct*ss12b
        th(2,5) = th(1,6)
        th(2,6) = forth*cs*(tre*ss11b+tsxb-tre*spxb-tpxb)
        th(3,3) = ss22a+half*(dgab+duab)+dsxb-ssgx-ssgh-dsgh
        th(3,4) = half*(dgab-duab)
        th(4,4) = dsxa+half*(dgab+duab)+ss22b-ssgx-ssgh-dsgh
        th(4,5) = -half*st*ss12b
        th(4,6) = half*sqrt3*st*ss12b
        th(5,5) = spxa+forth*((ssab+tre*tsab)+st2*(ss11b+tre*tsxb) +ct2*(spxb+tre*tpxb))-dpux-two*dsgh
        th(5,6) = forth*sqrt3*(ssab-tsab-st2*(ss11b-tsxb)-ct2*(spxb-tpxb))
        th(6,6) = tpxa+forth*((tre*ssab+tsab)+st2*(tre*ss11b+tsxb) +ct2*(tre*spxb+tpxb))-dpux-two*dsgh
        !
        !     FILL OUT MATRIX AND DIAGONALIZE
        !
        do 30 i=1, ndim
            do 30 j=1, i
                th(i,j) = th(j,i)
        30 continue
        !
        call eign(ndim, th, vec, en, w, ndim)
        !     FOLLOWING STATEMENT ASSUMES ENERGY ZERO AT A + BC(R=RE)
        xter = en(nroot)+c(2)
        if( idebug .eq. 1 ) &
            & write(iptprt, '(1X,''NROOT,EN(1),EN(NROOT),C(2),XTER='',I5,1P4E13.5)') nroot, en(1), en(nroot), c(2), xter
        if( iflag .lt. zero ) write(iptprt, '(1X,''VECTOR='',6D20.8)') (vec(i,nroot), i=1, ndim)
        !
        !************************* END OF FPOT3N *******************************
        !
        !     IF(IAPP.EQ.0) RETURN
        if( nder .ne. 1 ) return
        !
        !     SET  VECTOR FOR STATE OF INTEREST
        do 40 i=1, ndim
            gvec(i) = vec(i,nroot)
        40 continue
        !     COMPUTE PARTIALS OF U=COS(THETA)
        dtdr1 = (r(1)*r(1)+r(2)*r(2)-r(3)*r(3))/(two*r(1)*r(1)*r(3))
        dtdr2 = -r(2)/(r(1)*r(3))
        dtdr3 = (r(2)*r(2)+r(3)*r(3)-r(1)*r(1))/(two*r(3)*r(3)*r(1))
        factr = -one/st
        if( abs(st) .lt. 1.d-6 ) write(iptprt, &
            & '(/,2X,T5,''Warning: Sine(THETA) = '',1PE20.10, '' in the derivative calculatio' // &
            & 'n'', /,2X,T14,''a serious error may result'')') st
        !     COMPUTE PARTIALS OF HAMILTONIAN ELEMENTS
        xtra    = +half*ct*(ss11b+tre*tsxb-spxb-tre*tpxb)
        dh11(1) = dss11a+xtra*dtdr1
        dh11(2) = forth*(dssab+tre*dtsab)+xtra*dtdr2
        dh11(3) = forth*(ct2*(dss11b+tre*dtsxb)+st2*(dspxb+tre*dtpxb))+xtra* dtdr3
        xtra    = -half*sqrt3*ct*(ss11b-tsxb-spxb+tpxb)
        dh12(1) = xtra*dtdr1
        dh12(2) = forth*sqrt3*(dssab-dtsab)+xtra*dtdr2
        dh12(3) = -forth*sqrt3*(ct2*(dss11b-dtsxb)+st2*(dspxb-dtpxb))+xtra* dtdr3
        dh13(1) = dss12a
        dh13(2) = zero
        dh13(3) = zero
        xtra    = -half*ss12b
        dh14(1) = xtra*dtdr1
        dh14(2) = xtra*dtdr2
        dh14(3) = -half*ct*dss12b+xtra*dtdr3
        xtra    = forth*(ct2-st2)*(ss11b+tre*tsxb-spxb-tre*tpxb)*factr
        dh15(1) = xtra*dtdr1
        dh15(2) = xtra*dtdr2
        dh15(3) = forth*cs*(dss11b+tre*dtsxb-dspxb-tre*dtpxb)+xtra*dtdr3
        xtra    = -forth*sqrt3*(ct2-st2)*(ss11b-tsxb-spxb+tpxb)*factr
        dh16(1) = xtra*dtdr1
        dh16(2) = xtra*dtdr2
        dh16(3) = -forth*sqrt3*cs*(dss11b-dtsxb-dspxb+dtpxb)+xtra*dtdr3
        xtra    = +half*ct*(tre*ss11b+tsxb-tre*spxb-tpxb)
        dh22(1) = dtsxa+xtra*dtdr1
        dh22(2) = forth*(tre*dssab+dtsab)+xtra*dtdr2
        dh22(3) = forth*(ct2*(tre*dss11b+dtsxb)+st2*(tre*dspxb+dtpxb))+xtra* dtdr3
        xtra    = +half*sqrt3*ss12b
        dh24(1) = xtra*dtdr1
        dh24(2) = xtra*dtdr2
        dh24(3) = half*sqrt3*ct*dss12b+xtra*dtdr3
        dh25(1) = dh16(1)
        dh25(2) = dh16(2)
        dh25(3) = dh16(3)
        xtra    = forth*(ct2-st2)*(tre*ss11b+tsxb-tre*spxb-tpxb)*factr
        dh26(1) = xtra*dtdr1
        dh26(2) = xtra*dtdr2
        dh26(3) = forth*cs*(tre*dss11b+dtsxb-tre*dspxb-dtpxb)+xtra*dtdr3
        dh33(1) = dss22a
        dh33(2) = half*(ddgab+dduab)
        dh33(3) = ddsxb
        dh34(1) = zero
        dh34(2) = half*(ddgab-dduab)
        dh34(3) = zero
        dh44(1) = ddsxa
        dh44(2) = half*(ddgab+dduab)
        dh44(3) = dss22b
        xtra    = -half*ct*ss12b*factr
        dh45(1) = xtra*dtdr1
        dh45(2) = xtra*dtdr2
        dh45(3) = -half*st*dss12b+xtra*dtdr3
        xtra    = half*sqrt3*ct*ss12b*factr
        dh46(1) = xtra*dtdr1
        dh46(2) = xtra*dtdr2
        dh46(3) = half*sqrt3*st*dss12b+xtra*dtdr3
        xtra    = -half*ct*(ss11b+tre*tsxb-spxb-tre*tpxb)
        dh55(1) = dspxa+xtra*dtdr1
        dh55(2) = forth*(dssab+tre*dtsab)+xtra*dtdr2
        dh55(3) = forth*(st2*(dss11b+tre*dtsxb)+ct2*(dspxb+tre*dtpxb))+xtra* dtdr3
        xtra    = +half*sqrt3*ct*(ss11b-tsxb-spxb+tpxb)
        dh56(1) = xtra*dtdr1
        dh56(2) = forth*sqrt3*(dssab-dtsab)+xtra*dtdr2
        dh56(3) = -forth*sqrt3*(st2*(dss11b-dtsxb)+ct2*(dspxb-dtpxb))+xtra* dtdr3
        xtra    = -half*ct*(tre*ss11b+tsxb-tre*spxb-tpxb)
        dh66(1) = dtpxa+xtra*dtdr1
        dh66(2) = forth*(tre*dssab+dtsab)+xtra*dtdr2
        dh66(3) = forth*(st2*(tre*dss11b+dtsxb)+ct2*(tre*dspxb+dtpxb))+xtra* dtdr3
        !
        !     LOOP OVER R(I) DERIVATIVES
        !     FIRST ZERO OUT UNUSED ELEMENTS
        dh(2,3) = zero
        dh(3,5) = zero
        dh(3,6) = zero
        do 50 i=1, 3
            !
            !     SET UP DERIVATIVE MATRIX
            !
            dh(1,1) = dh11(i)
            dh(1,2) = dh12(i)
            dh(1,3) = dh13(i)
            dh(1,4) = dh14(i)
            dh(1,5) = dh15(i)
            dh(1,6) = dh16(i)
            dh(2,2) = dh22(i)
            dh(2,4) = dh24(i)
            dh(2,5) = dh25(i)
            dh(2,6) = dh26(i)
            dh(3,3) = dh33(i)
            dh(3,4) = dh34(i)
            dh(4,4) = dh44(i)
            dh(4,5) = dh45(i)
            dh(4,6) = dh46(i)
            dh(5,5) = dh55(i)
            dh(5,6) = dh56(i)
            dh(6,6) = dh66(i)
            !
            !     FILL OUT DERIVATIVE MATRIX
            do 70 j=1, ndim
                do 70 k=1, j
                    dh(j,k) = dh(k,j)
            70 continue
            !
            dnumb = zero
            do 80 j=1, ndim
                do 80 k=1, ndim
                    dnumb = dnumb+gvec(j)*dh(j,k)*gvec(k)
            80 continue
            dd(i) = dnumb
            !     END OF DERIVATIVE LOOP
        50 continue
        !     PUT DERIVATIVES INTO CORRESPONDING SPOTS FOR CORRECT COORD CHOICE.
        do 90 ix=1, 3
        90 dq(indx(ix))=dd(ix)
        !
        !************************ END OF DFPOT *********************************
        !
        !
        contains
        !> @brief Evaluate SETSPL for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine setspl
            !     SETS UP SPLINE FITS FOR H2 CURVES
            implicit none
            double precision :: epsln, ss, ss1, ss2, t
            integer          :: i
            !     NOW COMPUTE SPLINE FIT DATA
            !     H2 SINGLET STATE
            do 10 i=1, 67
            10 r1h2(i)=dble(i+3)*0.1d0
            r1h2(68) = 7.2d-00
            r1h2(69) = 7.4d-00
            r1h2(70) = 7.6d-00
            r1h2(71) = 7.8d-00
            r1h2(72) = 8.0d-00
            r1h2(73) = 8.25d-00
            r1h2(74) = 8.5d-00
            r1h2(75) = 9.0d-00
            r1h2(76) = 9.5d-00
            r1h2(77) = 10.d-00
            r1h2(78) = 11.d-00
            r1h2(79) = 12.d-00
            do 30 i=1, 79
            30 v1h2(i)=v1h2(i)+one
            s2ss(1)  = 3.2006875d+1
            s2ss(79) = -7.93478d-7
            call splnin(79, r1h2, v1h2, t, ss, ss1, ss2, epsln, s2ss, s3ss, delss)
            !     H2 TRIPLET STATE
            do 70 i=1, 81
            70 r2h2(i)=dble(i+9)*0.1d-00
            r2h2(82) = 9.25d-00
            r2h2(83) = 9.5d-00
            r2h2(84) = 9.75d-00
            r2h2(85) = 10.d-00
            r2h2(86) = 11.d-00
            r2h2(87) = 12.d-00
            do 90 i=1, 87
            90 v2h2(i)=v2h2(i)+one
            s2ts(1)  = 1.3039390d0
            s2ts(87) = -7.93478d-7
            call splnin(87, r2h2, v2h2, t, ss, ss1, ss2, epsln, s2ts, s3ts, delts)
            !     H2+ G STATE
            do 130 i=1, 100
            130 r3h2(i)=dble(i)*0.1d-00
            do 140 i=1, 26
            140 r3h2(i+100)=dble(i)*half+10.d-00
            do 150 i=1, 126
            150 v3h2(i)=v3h2(i)+one
            s2dg(1)   = 1.0742890d+3
            s2dg(126) = -3.03938d-7
            call splnin(126, r3h2, v3h2, t, ss, ss1, ss2, epsln, s2dg, s3dg, deldg)
            !     H2+ U STATE
            do 170 i=1, 126
            170 v4h2(i)=v4h2(i)+one
            s2du(1)   = 9.0876397d+2
            s2du(126) = -3.0398d-7
            call splnin(126, r3h2, v4h2, t, ss, ss1, ss2, epsln, s2du, s3du, deldu)
        end subroutine setspl
        !> @brief Evaluate HCLM for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine hclm( r, xhm2s, dxhm2s )
            implicit none
            double precision :: a, a1, arg, b, b1, c, c1, dx, dxhm2s, dy, r, rsmall, x, xhm2s, y
            !     DOUBLET SIGMA STATE OF HCL-IN A.U.
            a      = 2.52966d+1
            b      = -2.23d+0
            c      = -1.32774d-1
            a1     = 2.4913064d+4
            b1     = 7.1248778d+0
            c1     = 2.8830121
            rsmall = 0.94486d+0
            if( r .ge. rsmall ) go to 10
            !     SMALL R FIT USED TO BE CONSISTENT WITH BNL SPLINE FIT VERSION
            arg    = a1*exp(-b1*r)
            xhm2s  = arg+c1
            dxhm2s = -b1*arg
            go to 20
            10 dx=one/r**7
            x  = dx*r
            dy = x*r
            y  = dy*r
            !     COMPUTE ENERGY AND DERIVATIVE
            xhm2s  = a*x+b*y+c
            dxhm2s = -six*a*dx-for*b*dy
            20 return
        end subroutine hclm
        !> @brief Evaluate HCL for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine hcl( r, h, xh3s, xh1p, xh3p, dh, dxh3s, dxh1p, dxh3p )
            implicit none
            double precision :: a2, aa, ab, alp, ang, angcon, arg, b0, b2, ba, bb, bet, bt, bv, bx, c, c0, c2, c4, c6, c8, &
                & cc, cs, cv, d0, da, da2, dab, dalp, dang, db2, dbet, dc, dc2, dcs, de1, de12, de2, dfac1, dfac2, dg, dg11, &
                & dg12, dg22, dh, dh11, dh12, dh22, dovlp, ds, ds12, ds2, dt, dt1, dt2, dv, dx, dxh1p, dxh3p, dxh3s, e1, &
                & e12, e2, ex, ex1, ex2, fac1, fac2, g, g11, g12, g22, gb, gmax, h, h11, h12, h22, omex, opex, ovlp, r, r0, &
                & r2, ra, ri, ri2, ri4, ri6, ri8, rmax, root, rops2, rt, rv, rv1, rx, s, s12, s2, t1, t2, ttt, xh1p, xh3p, &
                & xh3s, y, z1, z2, z3, z4, z5, z6, zeta
            !     CALCULATES 2X2 SINGLET SIGMA BLOCK, AND TRIPLET SIGMA, SINGLET
            !     PI, AND TRIPLET PI CURVES FOR HCL IN A.U.
            dimension h(2,2), dh(2,2)
            !
            !     TRIPLET SIGMA DE ALREADY CONTAINS SATO=0.06
            !     ANGCON = +1 FOR #2 PERP TO #1, = -1 FOR #1 PERP TO #2,
            !        AND = 0 FOR SYMMETRIC ORTHOGONALIZATION
            !     SINGLET SIGMA BLOCK
            !     GROUND STATE ADIABAT
            dx     = 1.69552d-1
            rx     = 2.4082d+0
            bx     = 9.89104d-1
            dv     = 1.84605d-1
            rv     = 4.746d+0
            bv     = 2.77985d-1
            cv     = 1.82621d-1
            c8     = -2.018460d+4
            c6     = 2.399677d+3
            c4     = -5.351085d+1
            c0     = 3.67226d-1
            aa     = 25.20626238d-0
            bb     = 1.910495693d-0
            cc     = 1.828090476d-1
            rt     = 2.53d-0
            bt     = 1.047d-0
            dt     = 6.685927358d-2
            ra     = 2.274408d+0
            ba     = 1.069975d+0
            da     = 5.3556d-2
            rmax   = 2.5d+0
            gmax   = 2.0212d-2
            gb     = 4.654d-1
            rv1    = 3.5d-0
            angcon = 1.0d-0
            zeta   = 2.0387d0
            z1     = 7.222580119d+2
            z2     = -2.749954467d+2
            z3     = 8.648352007d+2
            z4     = 2.437196153d+2
            z5     = 3.144371653d+1
            z6     = 1.128252224d-1
            r0     = 2.53d-0
            b0     = 0.9897d-0
            d0     = 0.150789d-0
            ex     = exp(-bx*(r-rx))
            omex   = one-ex
            e1     = dx*omex*omex-dx
            de1    = two*dx*bx*ex*omex
            !     EXCITED STATE ADIABAT
            if( r .gt. rv ) go to 20
            if( r .lt. rv1 ) go to 60
            !     MORSE FIT
            ex   = exp(-bv*(r-rv))
            omex = one-ex
            e2   = dv*omex*omex+cv
            de2  = two*dv*bv*ex*omex
            go to 30
            !     LONG RANGE FIT
            20 ri=one/r
            ri2 = ri*ri
            ri4 = ri2*ri2
            ri6 = ri4*ri2
            ri8 = ri4*ri4
            e2  = c8*ri8+c6*ri6+c4*ri4-ri+c0
            de2 = -ri*(eht*c8*ri8+six*c6*ri6+for*c4*ri4-ri)
            go to 30
            !     SHORT RANGE FIT
            60 arg=aa*exp(-bb*r)
            e2  = arg+cc
            de2 = -bb*arg
            30 continue
            !     LOWER COVALENT DIABAT
            ex      = exp(-b0*(r-r0))
            omex    = one-ex
            h(1,1)  = d0*omex*omex-d0
            dh(1,1) = two*d0*b0*ex*omex
            !     OVERLAP
            r2    = r*r
            t1    = z1*(one+r)+z2*r2
            t2    = z1*(one+zeta*r)+z3*r2+z4*r*r2+z5*r2*r2
            ex1   = exp(-r)
            ex2   = exp(-zeta*r)
            dt1   = z1+two*z2*r
            dt2   = z1*zeta+two*z3*r+tre*z4*r2+for*z5*r2*r
            fac1  = -t1*ex1+t2*ex2
            dfac1 = (-dt1+t1)*ex1+(dt2-zeta*t2)*ex2
            ovlp  = z6*fac1/r2
            dovlp = z6*(dfac1-two*fac1/r)/r2
            rops2 = one/sqrt(one+ovlp*ovlp)
            s12   = ovlp*rops2/rt2i
            if( abs(s12) .gt. 1.0d0 ) write(iptprt, '(1X,''In the DIM potential in the subprogram HCL: '', ''R,S12='',F10.5,1PE1' // &
                '5.7)') r, s12
            ds12 = dovlp*rops2**3/rt2i
            !     H12 AND H22 DIABATS IN NON-ORTHOGONAL REPRESENTATION
            e12   = e1+e2
            de12  = de1+de2
            fac1  = e12*h(1,1)-e1*e2-h(1,1)*h(1,1)
            fac2  = one-s12*s12
            dfac2 = -two*s12*ds12
            ttt   = fac1*fac2
            root  = 0.d0
            if( ttt .gt. 0.d0 ) root=sqrt(ttt)
            h(1,2)  = s12*h(1,1)+root
            dh(1,2) = ds12*h(1,1)+s12*dh(1,1)
            ttt                          = half*(dfac2*fac1+fac2*(de12*h(1,1) +e12*dh(1,1)-e1*de2-de1*e2-two*h(1,1)* &
                & dh(1,1)))
            if( root .ne. 0.d0 ) dh(1,2) = dh(1,2) + ttt/root
            h(2,2)                       = e12*fac2+two*s12*h(1,2)-h(1,1)
            dh(2,2)                      = de12*fac2+e12*dfac2+two*(ds12*h(1,2)+s12*dh(1,2)) -dh(1,1)
            !     TRANSFORM TO SYMMETRIC ORTHOGONAL REPRESENTATION
            t1      = one/sqrt(one+s12)
            t2      = one/sqrt(one-s12)
            dt1     = -half*ds12*t1**3
            dt2     = half*ds12*t2**3
            alp     = half*(t1+t2)
            bet     = half*(t1-t2)
            dalp    = half*(dt1+dt2)
            dbet    = half*(dt1-dt2)
            a2      = alp*alp
            b2      = bet*bet
            ab      = alp*bet
            da2     = two*alp*dalp
            db2     = two*bet*dbet
            dab     = alp*dbet+dalp*bet
            g11     = a2*h(1,1)+two*ab*h(1,2)+b2*h(2,2)
            g12     = (a2+b2)*h(1,2)+ab*(h(1,1)+h(2,2))
            g22     = b2*h(1,1)+two*ab*h(1,2)+a2*h(2,2)
            arg     = two*(ab*dh(1,2)+dab*h(1,2))
            dg11    = a2*dh(1,1)+da2*h(1,1)+arg+b2*dh(2,2)+db2*h(2,2)
            dg12    = (a2+b2)*dh(1,2)+(da2+db2)*h(1,2)+ab*(dh(1,1)+dh(2,2)) +dab*(h(1,1)+h(2,2))
            dg22    = b2*dh(1,1)+db2*h(1,1)+arg+a2*dh(2,2)+da2*h(2,2)
            h(1,1)  = g11
            h(1,2)  = g12
            h(2,2)  = g22
            dh(1,1) = dg11
            dh(1,2) = dg12
            dh(2,2) = dg22
            !     PERFORM UNITARY TRANSFORMATION ON H
            y = half*(one/t1+one/t2)
            if( abs(y) .gt. 1.0d0 ) y = sign(1.0d0, y)
            ang  = angcon*acos(y)
            dang = zero
            if( abs(ang) .lt. 1.d-14 ) go to 90
            dang    = half*angcon*(dt1/(t1*t1)+dt2/(t2*t2))/sqrt(one-y*y)
            c       = cos(ang)
            s       = sin(ang)
            dc      = -s*dang
            ds      = c*dang
            c2      = c*c
            s2      = s*s
            cs      = c*s
            dc2     = two*c*dc
            ds2     = two*s*ds
            dcs     = c*ds+s*dc
            h11     = c2*h(1,1)+two*cs*h(1,2)+s2*h(2,2)
            dh11    = c2*dh(1,1)+two*cs*dh(1,2)+s2*dh(2,2) +dc2*h(1,1)+two*dcs*h(1,2)+ds2*h(2,2)
            h12     = (c2-s2)*h(1,2)+cs*(h(2,2)-h(1,1))
            dh12    = (c2-s2)*dh(1,2)+cs*(dh(2,2)-dh(1,1)) +(dc2-ds2)*h(1,2)+dcs*(h(2,2)-h(1,1))
            h22     = s2*h(1,1)-two*cs*h(1,2)+c2*h(2,2)
            dh22    = s2*dh(1,1)-two*cs*dh(1,2)+c2*dh(2,2) +ds2*h(1,1)-two*dcs*h(1,2)+dc2*h(2,2)
            h(1,1)  = h11
            h(2,2)  = h22
            h(1,2)  = h12
            h(2,1)  = h12
            dh(1,1) = dh11
            dh(2,2) = dh22
            dh(1,2) = dh12
            90 dh(2,1)=dh12
            !     TRIPLET SIGMA CURVE   (T STATE)
            ex    = exp(-bt*(r-rt))
            opex  = one+ex
            xh3s  = dt*opex*opex-dt
            dxh3s = -two*dt*bt*ex*opex
            !     SINGLET PI CURVE   (A STATE)
            ex    = exp(-ba*(r-ra))
            opex  = one+ex
            xh1p  = da*opex*opex-da
            dxh1p = -two*da*ba*ex*opex
            !     TRIPLET PI CURVE   (SMALL A STATE)
            g  = gmax
            dg = zero
            if( r .le. rmax ) go to 100
            g  = gmax*exp(-gb*(r-rmax)**2)
            dg = -two*gb*(r-rmax)*g
            100 xh3p=xh1p-g
            dxh3p = dxh1p-dg
        end subroutine hcl
        !> @brief Evaluate H2CURV for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine h2curv( r, ab1s, ab3s, ab2g, ab2u, dab1s, dab3s, dab2g, dab2u )
            !     CALCULATES ENERGY AND DERIVATIVES FOR H2 SINGLET AND TRIPLET
            !     SIGMA AND H2+ DOUBLET SIGMA G AND U STATES IN A.U.
            implicit none
            double precision :: a1, ab1s, ab2g, ab2u, ab3s, arg, arg2, arg3, b1, be, bfac, c10, c4, c6, c8, dab1s, dab2g, &
                & dab2u, dab3s, de, df1, epsln, ex, f1, f2, opex, r, r10, r2, r4, r6, r8, re, ss2
            !     H2 EXPONENTIALS W.R.T. FINAL LIMITS, H2+ W.R.T. 0.0 LIMIT
            !**   DATA A2,B2,C2/2.41893D00,2.1329115D00,0.0918539D00/
            !     H2 CURVES
            !     SINGLET SIGMA
            b1  = 6.0315697d-00
            a1  = 9.8214033d-00
            c6  = 6.49903d-00
            c8  = 124.399d-00
            c10 = 1135.21d-00
            de  = 3.21533d-2
            be  = 8.3d-1
            re  = 1.60d-0
            c4  = 2.25d0
            if( r .gt. 12.d0 ) go to 20
            if( r .ge. 0.4d0 ) go to 10
            !     SHORT-RANGE EXPONENTIAL FIT
            ab1s  = a1*exp(-b1*r)
            dab1s = -b1*ab1s
            go to 30
            !     SPLINE FIT
            10 call spline(79, r1h2, v1h2, r, ab1s, dab1s, ss2, epsln, s2ss, s3ss, delss)
            go to 30
            !     LONG-RANGE FORM
            20 r2=one/(r*r)
            r6    = r2**3
            r8    = r6*r2
            r10   = r8*r2
            ab1s  = -(c6*r6+c8*r8+c10*r10)
            dab1s = (six*c6*r6+eht*c8*r8+10.d0*c10*r10)/r
            30 continue
            !     TRIPLET SIGMA
            !     LINES STARTING WITH C** ARE FOR CORRECT TRIPLET CURVE
            !**   IF(R.GT.12.D0) GO TO 50
            !**   IF(R.GT.ONE) GO TO 40
            !     SHORT-RANGE EXPONENTIAL FIT
            !**   AB3S=A2*EXP(-B2*R)+C2
            !**   DAB3S=-B2*(AB3S-C2)
            !**   GO TO 60
            !     SPLINE FIT
            !**40 CALL SPLINE(87,R2H2,V2H2,R,AB3S,DAB3S,SS2,EPSLN,S2TS,S3TS,DELTS)
            !**   GO TO 60
            !     LONG-RANGE FORM -- SAME AS FOR SINGLET STATE
            !**   AB3S=AB1S
            !**   DAB3S=DAB1S
            !     REPLACE TRIPLET CURVE WITH MODIFIED ANTI-MORSE
            bfac = 2.6d-2
            if( r .lt. re ) bfac=2.026d0
            arg   = r-re
            ex    = exp(-be*arg)
            opex  = one+ex
            f1    = de*opex*opex-de
            df1   = -two*de*be*ex*opex
            arg2  = arg*arg
            arg3  = arg2*arg
            f2    = exp(-bfac*arg3)
            ab3s  = f1*f2
            dab3s = (df1-f1*tre*bfac*arg2)*f2
            !     H2+ CURVES
            !     DOUBLET G
            60 if(r.gt.23.d0) go to 70
            !     SPLINE FIT
            call spline(126, r3h2, v3h2, r, ab2g, dab2g, ss2, epsln, s2dg, s3dg, deldg)
            go to 80
            !     LONG-RANGE FORM
            70 r4=one/r**4
            ab2g  = -c4*r4+half
            dab2g = 9.d0*r4/r
            !     DOUBLET U
            80 if(r.gt.23.d0) go to 90
            !     SPLINE FIT
            call spline(126, r3h2, v4h2, r, ab2u, dab2u, ss2, epsln, s2du, s3du, deldu)
            go to 100
            !     LONG-RANGE FORM -- SAME AS FOR DOUBLET G STATE
            90 ab2u=ab2g
            dab2u = dab2g
            100 return
        end subroutine h2curv
        !> @brief Evaluate SPLNIN for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine splnin( n, x, y, t, ss, ss1, ss2, epsln, s2, s3, dely )
            !     JIM STINE'S INITIALIZATION ROUTINE FOR SPLINE FITTING
            !     N=NO. POINTS (X,Y); T=ARGUMENT; SS=VALUE; SS1=1ST DERIV;
            !     SS2=2ND DERIV; S2,S3, AND DELY ARE COEFFICIENT ARRAYS
            !     SECOND DERIV'S AT ENDS STORED IN S2(1) AND S2(N) ON INPUT
            !     EPSLN IS OVERALL FIT CRITERION
            implicit none
            double precision :: b, c, delsqy, dely, epsln, eta, ex, h, h2, half, omega, s2, s3, ss, ss1, ss2, t, w, x, y
            integer          :: i, knt, n, n1
            dimension c(126), h(126), h2(126), b(126)
            dimension x(n), y(n), s2(n), s3(n), dely(n)
            ex    = 1.1d0
            omega = 1.0717968d0
            half  = 0.5d0
            epsln = 1.d-12
            n1    = n-1
            knt   = 0
            do 51 i=1, n1
                h(i) = x(i+1)-x(i)
            51 dely(i)=(y(i+1)-y(i))/h(i)
            do 52 i=2, n1
                h2(i)  = h(i-1)+h(i)
                b(i)   = half*h(i-1)/h2(i)
                delsqy = (dely(i)-dely(i-1))/h2(i)
                s2(i)  = delsqy+delsqy
            52 c(i)=delsqy+s2(i)
            5 eta=0.0d0
            do 10 i=2, n1
                w = (c(i)-b(i)*s2(i-1)-(half-b(i))*s2(i+1)-s2(i))*omega
                if( abs(w)-eta ) 10, 10, 9
                9 eta=abs(w)
            10 s2(i)=s2(i)+w
            knt = knt+1
            if( knt .gt. 10 ) epsln=epsln*ex
            if( eta-epsln ) 14, 5, 5
            14 do 53 i=1, n1
            53 s3(i)=(s2(i+1)-s2(i))/h(i)
            if( idebug .eq. 1 ) write(iptprt, '(5X,11HFinal EPS =,D15.6)') epsln
        end subroutine splnin
        !> @brief Evaluate SPLINE for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine spline( n, x, y, t, ss, ss1, ss2, epsln, s2, s3, dely )
            !     JIM STINE'S SPLINE INTERPOLATOR
            implicit none
            double precision :: delsqs, dely, epsln, ht1, ht2, prod, s2, s3, ss, ss1, ss2, t, x, y
            integer          :: i, n
            dimension x(n), y(n), s2(n), s3(n), dely(n)
            i = 1
            if( t .gt. x(1) ) go to 10
            ss  = y(1)
            ss1 = dely(1)
            write(iptprt, '(1X,''In the DIMIMA SPLINE routine ARG is out of range:'', /,1X,'' ARG ='',1PE13.5,'', X='',1PE13.5' // &
                ')') t, x(1)
            return
            10 if(t.lt.x(n)) go to 20
            ss  = y(n)
            ss1 = dely(n)
            write(iptprt, '(1X,''In the DIMIMA SPLINE routine ARG is out of range:'', /,1X,'' ARG ='',1PE13.5,'', X='',1PE13.5' // &
                ')') t, x(1)
            return
            20 continue
            56 if(t-x(i)) 60, 17, 57
            57 i=i+1
            go to 56
            59 i=n
            60 i=i-1
            17 ht1=t-x(i)
            ht2    = t-x(i+1)
            prod   = ht1*ht2
            ss2    = s2(i)+ht1*s3(i)
            delsqs = (s2(i)+s2(i+1)+ss2)/6.d-00
            ss     = y(i)+ht1*dely(i)+prod*delsqs
            ss1    = dely(i)+(ht1+ht2)*delsqs+prod*s3(i)/6.d-00
            61 return
        end subroutine spline
        !> @brief Evaluate EIGN for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine eign( nn, a, vec, eig, w, nd )
            !     MATRIX DIAGONALIZATION ROUTINE FOR REAL SYMMETRIC CASE
            !     HOUSEHOLDER METHOD
            !     RHO=UPPERLIMIT FOR OFF-DIAGONAL ELEMENT
            !     NN=SIZE OF MATRIX
            !     A=MATRIX (ONLY LOWER TRIANGLE IS USED ,THIS IS DESTROYED)
            !     EIG=RETURNED EIGENVALUES IN ALGEBRAIC ASCENDING ORDER
            !     VEC=RETURNED EIGENVECTORS IN COLUMNS
            implicit none
            double precision :: a, a2, b, beta, betasq, cosa, cosap, d, dia, dif, eig, g, gamma, p, pp, ppbr, ppbs, q, qj, &
                & r1, r12, r2, rho, rhosq, s, sgn, shift, sina, sina2, sqrts, sum, temp, u, vec, w, wj, wtaw
            integer :: i, i1, iord, iposv, itemp, ivpos, j, k, l, lv, m, n, n1, n2, nd, nn, np, npas, nr, nrr, nt, nv
            dimension eig(nd), w(nd)
            dimension gamma(6), beta(6), betasq(6)
            dimension p(6), q(6)
            equivalence (p(1), beta(1)), (q(1), beta(1))
            dimension iposv(6), ivpos(6), iord(6)
            dimension a(nd,nd), vec(nd,nd)
            rho   = 1.0d-14
            rhosq = rho*rho
            n     = nn
            if( n .eq. 0 ) go to 560
            1 n1=n-1
            n2       = n-2
            gamma(1) = a(1,1)
            if( n2 ) 280, 270, 120
            120 do 260 nr=1, n2
                b = a(nr+1,nr)
                s = 0.0d0
                do 130 i=nr, n2
                    !     PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION
                130 s=s+a(i+2,nr)**2
                a(nr+1,nr) = 0.0d0
                if( s ) 250, 250, 140
                140 s=s+b*b
                sgn = +1.0d0
                if( b ) 150, 160, 160
                150 sgn=-1.0d0
                160 sqrts=sqrt(s)
                d          = sgn/(sqrts+sqrts)
                temp       = sqrt(.5d0+b*d)
                w(nr)      = temp
                a(nr+1,nr) = temp
                d          = d/temp
                b          = -sgn*sqrts
                !     D IS FACTOR OF PROPORTIONALITY  NOW COMPUTE AND SAVE W VECTOR
                !     EXTRA SINGLY SUBSCRIPTED W VECTOR USED FOR SPEED
                do 170 i=nr, n2
                    temp   = d*a(i+2,nr)
                    w(i+1) = temp
                170 a(i+2,nr)=temp
                !     PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR
                !     SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K)
                wtaw = 0.0d0
                do 220 i=nr, n1
                    sum = 0.0d0
                    do 180 j=nr, i
                    180 sum=sum+a(i+1,j+1)*w(j)
                    i1 = i+1
                    if( n1-i1 ) 210, 190, 190
                    190 do 200 j=i1, n1
                    200 sum=sum+a(j+1,i+1)*w(j)
                    210 p(i)=sum
                220 wtaw=wtaw+sum*w(i)
                !     P VECTOR AND SCALAR K NOW STORED, NEXT COMPUTE Q VECTOR.
                do 230 i=nr, n1
                    !     NOW FORM PAP MATRIX, REQUIRED PART
                230 q(i)=p(i)-wtaw*w(i)
                do 240 j=nr, n1
                    qj = q(j)
                    wj = w(j)
                    do 240 i=j, n1
                240 a(i+1,j+1)=a(i+1,j+1)-2.0d0*(w(i)*qj+wj*q(i))
                250 beta(nr)=b
                betasq(nr) = b*b
            260 gamma(nr+1)=a(nr+1,nr+1)
            270 b=a(n,n-1)
            beta(n-1)   = b
            betasq(n-1) = b*b
            gamma(n)    = a(n,n)
            280 betasq(n)=0.0d0
            !     ADJOIN AN IDENTITY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS
            do 290 j=1, n
                do 290 i=1, n
            290 vec(i,j)=0.0d0
            do 300 i=1, n
            300 vec(i,i)=1.0d0
            m    = n
            sum  = 0.0d0
            npas = 1
            go to 400
            310 sum=sum+shift
            cosa = 1.0d0
            g    = gamma(1)-shift
            pp   = g
            ppbs = pp*pp+betasq(1)
            ppbr = sqrt(ppbs)
            do 370 j=1, m
                cosap = cosa
                if( ppbs .ne. 0.0d0 ) go to 320
                311 sina=0.0d0
                sina2 = 0.0d0
                cosa  = 1.0d0
                go to 350
                320 sina=beta(j)/ppbr
                sina2 = betasq(j)/ppbs
                cosa  = pp/ppbr
                !     POSTMULTIPLY IDAENTITY BY P-TRANSPOSE MATRIX
                nt = j+npas
                if( nt .lt. n ) go to 330
                321 nt=n
                330 do 340 i=1, nt
                    temp       = cosa*vec(i,j)+sina*vec(i,j+1)
                    vec(i,j+1) = -sina*vec(i,j)+cosa*vec(i,j+1)
                340 vec(i,j)=temp
                350 dia=gamma(j+1)-shift
                u        = sina2*(g+dia)
                gamma(j) = g+u
                g        = dia-u
                pp       = dia*cosa-sina*cosap*beta(j)
                if( j .ne. m ) go to 360
                351 beta(j)=sina*pp
                betasq(j) = sina2*pp*pp
                go to 380
                360 ppbs=pp*pp+betasq(j+1)
                ppbr    = sqrt(ppbs)
                beta(j) = sina*ppbr
            370 betasq(j)=sina2*ppbs
            380 gamma(m+1)=g
            !     TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT
            npas = npas+1
            if( betasq(m) .gt. rhosq ) go to 410
            390 eig(m+1)=gamma(m+1)+sum
            400 beta(m)=0.0d0
            betasq(m) = 0.0d0
            m         = m-1
            if( m .eq. 0 ) go to 430
            401 if(betasq(m).le.rhosq) go to 390
            !     TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN VALUE
            !     AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT
            410 a2=gamma(m+1)
            r2   = .5d0*a2
            r1   = .5d0*gamma(m)
            r12  = r1+r2
            dif  = r1-r2
            temp = sqrt(dif*dif+betasq(m))
            r1   = r12+temp
            r2   = r12-temp
            dif  = abs(a2-r1)-abs(a2-r2)
            if( dif .lt. 0.0d0 ) go to 420
            411 shift=r2
            go to 310
            420 shift=r1
            go to 310
            430 eig(1)=gamma(1)+sum
            !     INITIALIZE AUXILARY TABLES REQUIRED FOR REARRANGING THE VECTORS
            do 440 j=1, n
                iposv(j) = j
                ivpos(j) = j
            440 iord(j)=j
            !     USE A TRANSPOSITON SORT TO ORDER THE EIGENVALUES
            m = n
            go to 470
            450 do 460 j=1, m
                if( eig(j) .le. eig(j+1) ) go to 460
                451 temp=eig(j)
                eig(j)    = eig(j+1)
                eig(j+1)  = temp
                itemp     = iord(j)
                iord(j)   = iord(j+1)
                iord(j+1) = itemp
            460 continue
            470 m=m-1
            if( m .ne. 0 ) go to 450
            471 if(n1.eq.0) go to 500
            472 do 490 l=1, n1
                nv = iord(l)
                np = iposv(nv)
                if( np .eq. l ) go to 490
                473 lv=ivpos(l)
                ivpos(np) = lv
                iposv(lv) = np
                do 480 i=1, n
                    temp     = vec(i,l)
                    vec(i,l) = vec(i,np)
                480 vec(i,np)=temp
            490 continue
            !     BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX
            500 do 550 nrr=1, n
                k = n1
                510 k=k-1
                if( k .le. 0 ) go to 550
                511 sum=0.0d0
                do 520 i=k, n1
                520 sum=sum+vec(i+1,nrr)*a(i+1,k)
                sum = sum+sum
                do 530 i=k, n1
                530 vec(i+1,nrr)=vec(i+1,nrr)-sum*a(i+1,k)
                go to 510
            550 continue
            560 return
            !  END OF EIGN
        end subroutine eign
    end subroutine evaluate_internal
end module chempotpy_clh2_leps_ttgi1_1985_pot_legacy_state_bf7e4ade

!> @brief Evaluate TTGI1 for one set of internal distances.
subroutine pot( q, nder, xter, dq )
    use chempotpy_clh2_leps_ttgi1_1985_pot_legacy_state_bf7e4ade, only: evaluate_internal
    implicit none
    double precision, intent(in)  :: q(3)
    integer, intent(in)           :: nder
    double precision, intent(out) :: xter, dq(3)
    call evaluate_internal(q, nder, xter, dq)
end subroutine pot

!> @brief Evaluate POTINFO for the ClH2_LEPS_TTGI1_1985 potential.
