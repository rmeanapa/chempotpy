!   System:          ClH2
!   Functional form: Isaacson and Muckerman's diatomics-in-molecules ver
!   Common name:     ClH2 DIMA
!   Reference:       unpublished
!   Cross Reference: S. C. Tucker, D. G. Truhlar,
!                    G. C. Garrett, and A. D. Isaacson
!                    J. Chem. Phys. 82, 4102 (1985)
!   Calling Sequence:
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
!      The classical potential energy is set equal to zero for the Cl
!      infinitely far from the H2 diatomic and R(H2) set equal to the
!      H2 equilibrium diatomic value.
!   Parameters:
!      Set in the BLOCK DATA subprogram PTPARM
!   Coordinates:
!      Internal, Definition: R(1) = R(Cl-H)
!                            R(2) = R(H-H)
!                            R(3) = R(Cl-H)
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
!   Default Parameter Values:
!      Variable      Default value
!      NDER             1
!      IDEBUG           0
!      IPTPRT           6

!> @brief Implement ClH2_LEPS_TTGI1_1985, a potential energy surface for ClH2.
module chempotpy_clh2_leps_ttgi1_1985_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the ClH2_LEPS_TTGI1_1985 potential energy surface for ClH2.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,              only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_geometry,           only: flatten_cartesian_coordinates, transform_internal_gradient_to_cartesian
        use chempotpy_units,              only: bohr_to_angstrom, hartree_to_electronvolt
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
        real(wp)              :: dq(3), q(3), r2(3), tx(9), dedx(9), drdx(3,9), v, xter
        integer               :: iatom, idir, istate, j
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        call flatten_cartesian_coordinates(x, tx)
        ! input cartesian is ClHH
        call triatomic_distances_to_bohr(x, q)
        dq = 0.e0_wp
        call pot(q, igrad, xter, dq)
        xter = xter*hartree_to_electronvolt
        dq   = dq*hartree_per_bohr_to_ev_per_angstrom
        r2   = q*bohr_to_angstrom
        call evdrdx(tx, r2, drdx)
        call transform_internal_gradient_to_cartesian(dq, drdx, dedx)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = xter
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = xter
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
end module chempotpy_clh2_leps_ttgi1_1985_pes

!> @brief Expose ClH2_LEPS_TTGI1_1985 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_clh2_leps_ttgi1_1985_pes, only: evaluate
    use chempotpy_kinds,                    only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes

!> @brief Evaluate EvdRdX for the ClH2_LEPS_TTGI1_1985 potential.
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

module chempotpy_clh2_leps_ttgi1_1985_pot_legacy_state_bf7e4ade
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_internal

    contains

    !> @brief Evaluate dispatch.
    subroutine evaluate_internal( q, nder, xter, dq )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: q(3)
        integer, intent(in)   :: nder
        real(wp), intent(out) :: xter, dq(3)
        real(wp)              :: v1h2(79), v2h2(87), v3h2(126), v4h2(126)
        real(wp)              :: r1h2(79), r2h2(87), r3h2(126)
        real(wp)              :: s2ss(79), s3ss(79), delss(79)
        real(wp)              :: s2ts(87), s3ts(87), delts(87)
        real(wp)              :: s2dg(126), s3dg(126), deldg(126)
        real(wp)              :: s2du(126), s3du(126), deldu(126)
        real(wp)              :: quar, rt3
        real(wp)              :: arg, c(3), cs, ct, ct2, dd(3), ddgab, ddsxa, ddsxb, dduab, dgab, dh(6,6), dh1(2,2), dh11(3), &
            & dh12(3), dh13(3), dh14(3), dh15(3), dh16(3), dh2(2, 2), dh22(3), dh24(3), dh25(3), dh26(3), dh33(3), dh34(3), &
            & dh44(3), dh45(3), dh46(3), dh55(3), dh56(3), dh66(3), dnumb, dpux, dsgh, dspxa, dspxb, dss11a, dss11b, dss12a, &
            & dss12b, dss22a, dss22b, dssab, dsxa, dsxb, dtdr1, dtdr2, dtdr3, dtpxa, dtpxb, dtsab, dtsxa, dtsxb, duab, easyab, &
            & easyac, easybc, eht, en(6), factr, for, forth, gvec(6), h1(2, 2), h2(2,2), half, one, r(3), rt2i, six, spxa, spxb, &
            & sqrt3, ss11a, ss11b, ss12a, ss12b, ss22a, ss22b, ssab, ssgh, ssgx, st, st2, test, th(6, 6), tpxa, tpxb, tre, tsab, &
            & tsxa, tsxb, two, vec(6,6), w(6), xtra, zero
        integer           :: i, idebug, idum, iflag, indx(3), iptprt, ix, j, k, ndim, ndum, nroot
        character(len=80) :: title(2)
        equivalence (h1(1, 1), ss11a), (h1(1, 2), ss12a), (h1(2, 2), ss22a)
        equivalence (h2(1,1), ss11b), (h2(1,2), ss12b), (h2(2,2), ss22b)
        equivalence (dh1(1,1), dss11a), (dh1(1,2), dss12a), (dh1(2,2), dss22a)
        equivalence (dh2(1, 1), dss11b), (dh2(1, 2), dss12b), (dh2(2, 2), dss22b)
        equivalence (th(1, 1), dh(1,1))
        !   Initialize the array INDX which specifies the order of the
        !   input bond distances.  INDX=1,2,3 for Cl + H2 and INDX=132 for H + C
        !   Echo the name of the potential to the file linked to FORTRAN unit IP
        !     SET UP SPLINE FITS TO H2 CURVES
        arg    = 0.e0_wp
        c      = 0.e0_wp
        cs     = 0.e0_wp
        ct     = 0.e0_wp
        ct2    = 0.e0_wp
        dd     = 0.e0_wp
        ddgab  = 0.e0_wp
        ddsxa  = 0.e0_wp
        ddsxb  = 0.e0_wp
        dduab  = 0.e0_wp
        deldg  = 0.e0_wp
        deldu  = 0.e0_wp
        delss  = 0.e0_wp
        delts  = 0.e0_wp
        dgab   = 0.e0_wp
        dh     = 0.e0_wp
        dh1    = 0.e0_wp
        dh11   = 0.e0_wp
        dh12   = 0.e0_wp
        dh13   = 0.e0_wp
        dh14   = 0.e0_wp
        dh15   = 0.e0_wp
        dh16   = 0.e0_wp
        dh2    = 0.e0_wp
        dh22   = 0.e0_wp
        dh24   = 0.e0_wp
        dh25   = 0.e0_wp
        dh26   = 0.e0_wp
        dh33   = 0.e0_wp
        dh34   = 0.e0_wp
        dh44   = 0.e0_wp
        dh45   = 0.e0_wp
        dh46   = 0.e0_wp
        dh55   = 0.e0_wp
        dh56   = 0.e0_wp
        dh66   = 0.e0_wp
        dnumb  = 0.e0_wp
        dpux   = 0.e0_wp
        dsgh   = 0.e0_wp
        dspxa  = 0.e0_wp
        dspxb  = 0.e0_wp
        dss11a = 0.e0_wp
        dss11b = 0.e0_wp
        dss12a = 0.e0_wp
        dss12b = 0.e0_wp
        dss22a = 0.e0_wp
        dss22b = 0.e0_wp
        dssab  = 0.e0_wp
        dsxa   = 0.e0_wp
        dsxb   = 0.e0_wp
        dtdr1  = 0.e0_wp
        dtdr2  = 0.e0_wp
        dtdr3  = 0.e0_wp
        dtpxa  = 0.e0_wp
        dtpxb  = 0.e0_wp
        dtsab  = 0.e0_wp
        dtsxa  = 0.e0_wp
        dtsxb  = 0.e0_wp
        duab   = 0.e0_wp
        easyab = 0.e0_wp
        easyac = 0.e0_wp
        easybc = 0.e0_wp
        eht    = 0.e0_wp
        en     = 0.e0_wp
        factr  = 0.e0_wp
        for    = 0.e0_wp
        forth  = 0.e0_wp
        gvec   = 0.e0_wp
        h1     = 0.e0_wp
        h2     = 0.e0_wp
        half   = 0.e0_wp
        one    = 0.e0_wp
        quar   = 0.e0_wp
        r      = 0.e0_wp
        r1h2   = 0.e0_wp
        r2h2   = 0.e0_wp
        r3h2   = 0.e0_wp
        rt2i   = 0.e0_wp
        rt3    = 0.e0_wp
        s2dg   = 0.e0_wp
        s2du   = 0.e0_wp
        s2ss   = 0.e0_wp
        s2ts   = 0.e0_wp
        s3dg   = 0.e0_wp
        s3du   = 0.e0_wp
        s3ss   = 0.e0_wp
        s3ts   = 0.e0_wp
        six    = 0.e0_wp
        spxa   = 0.e0_wp
        spxb   = 0.e0_wp
        sqrt3  = 0.e0_wp
        ss11a  = 0.e0_wp
        ss11b  = 0.e0_wp
        ss12a  = 0.e0_wp
        ss12b  = 0.e0_wp
        ss22a  = 0.e0_wp
        ss22b  = 0.e0_wp
        ssab   = 0.e0_wp
        ssgh   = 0.e0_wp
        ssgx   = 0.e0_wp
        st     = 0.e0_wp
        st2    = 0.e0_wp
        test   = 0.e0_wp
        th     = 0.e0_wp
        tpxa   = 0.e0_wp
        tpxb   = 0.e0_wp
        tre    = 0.e0_wp
        tsab   = 0.e0_wp
        tsxa   = 0.e0_wp
        tsxb   = 0.e0_wp
        two    = 0.e0_wp
        v1h2   = 0.e0_wp
        v2h2   = 0.e0_wp
        v3h2   = 0.e0_wp
        v4h2   = 0.e0_wp
        vec    = 0.e0_wp
        w      = 0.e0_wp
        xtra   = 0.e0_wp
        zero   = 0.e0_wp
        iptprt = 6
        idebug = 0
        zero   = 0.e0_wp
        one    = 1.e0_wp
        two    = 2.e0_wp
        tre    = 3.e0_wp
        for    = 4.e0_wp
        six    = 6.e0_wp
        eht    = 8.e0_wp
        half   = 0.5e0_wp
        quar   = 0.25e0_wp
        rt2i   = 7.07106781186548e-1_wp
        rt3    = 1.73205080756888e0_wp
        dpux   = 0.e0_wp
        ssgx   = -0.132774e0_wp
        dsgh   = 0.e0_wp
        ssgh   = 0.5e0_wp
        c      = [ 0.169552e0_wp, 0.1744734e0_wp, 0.169552e0_wp]
        v1h2 = [ -.1202028e0_wp, -.5266270e0_wp, -.7696341e0_wp, -.9220261e0_wp, -1.0200556e0_wp, -1.0836422e0_wp, - &
            & 1.1245385e0_wp, - 1.1500562e0_wp, -1.1649342e0_wp, -1.1723459e0_wp, -1.1744744e0_wp, -1.1728537e0_wp, - &
            & 1.1685799e0_wp, -1.1624570e0_wp, - 1.1550670e0_wp, -1.1468496e0_wp, -1.1381312e0_wp, -1.1291562e0_wp, - &
            & 1.1201233e0_wp, -1.1111725e0_wp, -1.1024127e0_wp, - 1.0939273e0_wp, -1.0857810e0_wp, -1.0780164e0_wp, - &
            & 1.0706700e0_wp, -1.0637641e0_wp, -1.0573118e0_wp, -1.0513185e0_wp, - 1.0457832e0_wp, -1.0407003e0_wp, - &
            & 1.0360578e0_wp, -1.0318402e0_wp, -1.0280272e0_wp, -1.0245978e0_wp, -1.0215297e0_wp, - 1.0187967e0_wp, - &
            & 1.0163689e0_wp, -1.0142247e0_wp, -1.0123371e0_wp, -1.0106810e0_wp, -1.0092303e0_wp, -1.0079682e0_wp, - &
            & 1.0068703e0_wp, -1.0059178e0_wp, -1.0050923e0_wp, -1.0043782e0_wp, -1.0037626e0_wp, -1.0032309e0_wp, - &
            & 1.0027740e0_wp, - 1.0023800e0_wp, -1.0020423e0_wp, -1.0017521e0_wp, -1.0015030e0_wp, -1.0012899e0_wp, - &
            & 1.0011069e0_wp, -1.0009498e0_wp, - 1.0008150e0_wp, -1.0007002e0_wp, -1.0006030e0_wp, -1.0005162e0_wp, - &
            & 1.0004466e0_wp, -1.0003864e0_wp, -1.0003328e0_wp, - 1.0002906e0_wp, -1.0002466e0_wp, -1.0002154e0_wp, - &
            & 1.0001889e0_wp, -1.0001434e0_wp, -1.0001086e0_wp, -1.0000868e0_wp, - 1.0000682e0_wp, -1.0000528e0_wp, - &
            & 1.0000404e0_wp, -1.0000314e0_wp, -1.0000185e0_wp, -1.0000121e0_wp, -1.0000091e0_wp, - 1.0000045e0_wp, - &
            & 1.0000025e0_wp]
        v2h2 = [ -.6215227e0_wp, -.6757268e0_wp, -.7189640e0_wp, -.7544033e0_wp, -.7841501e0_wp, -.8096095e0_wp, - &
            & .8317238e0_wp, - .8511347e0_wp, -.8682913e0_wp, -.8835186e0_wp, -.8970636e0_wp, -.9091230e0_wp, - &
            & .9198593e0_wp, -.9294123e0_wp, -.9379051e0_wp, - .9454463e0_wp, -.9521346e0_wp, -.9580585e0_wp, - &
            & .9632985e0_wp, -.9679271e0_wp, -.9720104e0_wp, -.9756071e0_wp, -.9787717e0_wp, - .9815517e0_wp, - &
            & .9839910e0_wp, -.9861279e0_wp, -.9879977e0_wp, -.9896312e0_wp, -.9910564e0_wp, -.9922979e0_wp, - &
            & .9933781e0_wp, - .9943164e0_wp, -.9951304e0_wp, -.9958355e0_wp, -.9964456e0_wp, -.9969715e0_wp, - &
            & .9974252e0_wp, -.9978159e0_wp, -.9981500e0_wp, - .9984382e0_wp, -.9986849e0_wp, -.9988960e0_wp, - &
            & .9990763e0_wp, -.9992292e0_wp, -.9993597e0_wp, -.9994703e0_wp, -.9995635e0_wp, - .9996432e0_wp, - &
            & .9997094e0_wp, -.9997655e0_wp, -.9998125e0_wp, -.9998515e0_wp, -.9998842e0_wp, -.9999112e0_wp, - &
            & .9999337e0_wp, - .9999521e0_wp, -.9999670e0_wp, -.9999795e0_wp, -.9999887e0_wp, -.9999966e0_wp, - &
            & 1.0000030e0_wp, -1.0000076e0_wp, - 1.0000117e0_wp, -1.0000143e0_wp, -1.0000167e0_wp, -1.0000180e0_wp, - &
            & 1.0000191e0_wp, -1.0000195e0_wp, -1.0000197e0_wp, - 1.0000194e0_wp, -1.0000196e0_wp, -1.0000192e0_wp, - &
            & 1.0000191e0_wp, -1.0000180e0_wp, -1.0000173e0_wp, -1.0000164e0_wp, - 1.0000158e0_wp, -1.0000152e0_wp, - &
            & 1.0000143e0_wp, -1.0000143e0_wp, -1.0000127e0_wp, -1.0000109e0_wp, -1.0000095e0_wp, - 1.0000076e0_wp, - &
            & 1.0000067e0_wp, -1.0000043e0_wp, -1.0000025e0_wp]
        v3h2(1:64) = [ 8.0217579e0_wp, 3.0713797e0_wp, 1.4666292e0_wp, .6992459e0_wp, .2650120e0_wp, -.0048180e0_wp, - &
            & .1826248e0_wp, - .3044801e0_wp, -.3902705e0_wp, -.4517863e0_wp, -.4964119e0_wp, -.5289745e0_wp, - &
            & .5527406e0_wp, -.5699835e0_wp, -.5823232e0_wp, - .5909372e0_wp, -.5966963e0_wp, -.6002536e0_wp, - &
            & .6021058e0_wp, -.6026342e0_wp, -.6021349e0_wp, -.6008396e0_wp, -.5989309e0_wp, - .5965536e0_wp, - &
            & .5938235e0_wp, -.5908332e0_wp, -.5876573e0_wp, -.5843560e0_wp, -.5809780e0_wp, -.5775629e0_wp, - &
            & .5741424e0_wp, - .5707426e0_wp, -.5673842e0_wp, -.5640840e0_wp, -.5608555e0_wp, -.5577093e0_wp, - &
            & .5546536e0_wp, -.5516947e0_wp, -.5488374e0_wp, - .5460849e0_wp, -.5434394e0_wp, -.5409022e0_wp, - &
            & .5384736e0_wp, -.5361531e0_wp, -.5339400e0_wp, -.5318328e0_wp, -.5298296e0_wp, - .5279281e0_wp, - &
            & .5261259e0_wp, -.5244203e0_wp, -.5228082e0_wp, -.5212866e0_wp, -.5198521e0_wp, -.5185016e0_wp, - &
            & .5172315e0_wp, - .5160385e0_wp, -.5149192e0_wp, -.5138701e0_wp, -.5128878e0_wp, -.5119690e0_wp, - &
            & .5111105e0_wp, -.5103089e0_wp, -.5095612e0_wp, - .5088644e0_wp]
        v3h2(65:126) = [ -.5082155e0_wp, -.5076116e0_wp, -.5070501e0_wp, -.5065283e0_wp, -.5060437e0_wp, - &
            & .5055940e0_wp, -.5051769e0_wp, - .5047902e0_wp, -.5044319e0_wp, -.5041001e0_wp, -.5037930e0_wp, - &
            & .5035087e0_wp, -.5032458e0_wp, -.5030027e0_wp, -.5027780e0_wp, - .5025704e0_wp, -.5023785e0_wp, - &
            & .5022014e0_wp, -.5020377e0_wp, -.5018866e0_wp, -.5017472e0_wp, -.5016184e0_wp, -.5014996e0_wp, - &
            & .5013900e0_wp, -.5012888e0_wp, -.5011954e0_wp, -.5011093e0_wp, -.5010298e0_wp, -.5009564e0_wp, - &
            & .5008887e0_wp, -.5008262e0_wp, - .5007685e0_wp, -.5007153e0_wp, -.5006661e0_wp, -.5006207e0_wp, - &
            & .5005787e0_wp, -.5004122e0_wp, -.5002992e0_wp, -.5002219e0_wp, - .5001683e0_wp, -.5001306e0_wp, - &
            & .5001035e0_wp, -.5000837e0_wp, -.5000690e0_wp, -.5000577e0_wp, -.5000489e0_wp, -.5000420e0_wp, - &
            & .5000364e0_wp, -.5000318e0_wp, -.5000279e0_wp, -.5000247e0_wp, -.5000220e0_wp, -.5000196e0_wp, - &
            & .5000176e0_wp, -.5000158e0_wp, - .5000142e0_wp, -.5000129e0_wp, -.5000117e0_wp, -.5000106e0_wp, - &
            & .5000097e0_wp, -.5000089e0_wp, -.5000081e0_wp]
        v4h2(1:64) = [ 9.4993326e0_wp, 4.4973226e0_wp, 2.8272874e0_wp, 1.9892158e0_wp, 1.4831145e0_wp, 1.1423563e0_wp, &
            & .8955959e0_wp, .7072541e0_wp, .5576760e0_wp, .4351864e0_wp, .3324666e0_wp, .2447306e0_wp, .1687347e0_wp, &
            & .1022057e0_wp, .0434985e0_wp, - .0086173e0_wp, -.0550906e0_wp, -.0966755e0_wp, -.1339863e0_wp, - &
            & .1675344e0_wp, -.1977530e0_wp, -.2250137e0_wp, -.2496386e0_wp, - .2719091e0_wp, -.2920720e0_wp, - &
            & .3103453e0_wp, -.3269215e0_wp, -.3419716e0_wp, -.3556477e0_wp, -.3680850e0_wp, -.3794043e0_wp, - &
            & .3897134e0_wp, -.3991089e0_wp, -.4076771e0_wp, -.4154957e0_wp, -.4226343e0_wp, -.4291554e0_wp, - &
            & .4351155e0_wp, -.4405653e0_wp, - .4455506e0_wp, -.4501129e0_wp, -.4542896e0_wp, -.4581145e0_wp, - &
            & .4616184e0_wp, -.4648290e0_wp, -.4677716e0_wp, -.4704692e0_wp, - .4729427e0_wp, -.4752111e0_wp, - &
            & .4772916e0_wp, -.4792001e0_wp, -.4809509e0_wp, -.4825572e0_wp, -.4840309e0_wp, -.4853831e0_wp, - &
            & .4866238e0_wp, -.4877622e0_wp, -.4888066e0_wp, -.4897648e0_wp, -.4906439e0_wp, -.4914502e0_wp, - &
            & .4921898e0_wp, -.4928681e0_wp, - .4934901e0_wp]
        v4h2(65:126) = [ -.4940603e0_wp, -.4945830e0_wp, -.4950621e0_wp, -.4955011e0_wp, -.4959033e0_wp, - &
            & .4962717e0_wp, -.4966090e0_wp, - .4969178e0_wp, -.4972005e0_wp, -.4974591e0_wp, -.4976956e0_wp, - &
            & .4979119e0_wp, -.4981096e0_wp, -.4982903e0_wp, -.4984553e0_wp, - .4986060e0_wp, -.4987435e0_wp, - &
            & .4988690e0_wp, -.4989834e0_wp, -.4990877e0_wp, -.4991827e0_wp, -.4992692e0_wp, -.4993480e0_wp, - &
            & .4994196e0_wp, -.4994847e0_wp, -.4995438e0_wp, -.4995975e0_wp, -.4996462e0_wp, -.4996904e0_wp, - &
            & .4997303e0_wp, -.4997665e0_wp, - .4997992e0_wp, -.4998288e0_wp, -.4998554e0_wp, -.4998794e0_wp, - &
            & .4999011e0_wp, -.4999800e0_wp, -.5000244e0_wp, -.5000475e0_wp, - .5000579e0_wp, -.5000608e0_wp, - &
            & .5000595e0_wp, -.5000560e0_wp, -.5000515e0_wp, -.5000467e0_wp, -.5000421e0_wp, -.5000377e0_wp, - &
            & .5000337e0_wp, -.5000301e0_wp, -.5000269e0_wp, -.5000240e0_wp, -.5000216e0_wp, -.5000194e0_wp, - &
            & .5000174e0_wp, -.5000157e0_wp, - .5000142e0_wp, -.5000129e0_wp, -.5000117e0_wp, -.5000106e0_wp, - &
            & .5000097e0_wp, -.5000089e0_wp, -.5000081e0_wp]
        forth    = quar
        sqrt3    = rt3
        ndim     = 6
        nroot    = 1
        test     = 0.999999999e0_wp
        indx     = [ 1, 2, 3]
        title(1) = 'Final BNL DIM surface for CL+H2 -- 5 and7 kcal/mole barriers'
        title(2) = 'Spline fits replaced with analytic forms where possible'
        call setspl
        !   Initialize the energy in the three asymptotic valleys
        easyab = c(1)
        easybc = c(2)
        easyac = c(3)
        !     ENTRY POT GIVE ENERGY AND PARTIALS W.R.T 3 DISTANCES
        !     INTERNUCLEAR DISTANCES STORED IN Q(3) AND R(3)
        !     NROOT-TH STATE ENERGY STORED IN XTER (OR TERM FOR DIATOMIC)
        !     THREE PARTIALS STORED IN DD(3)
        !   Check the value of NDER
        !     SET FLAG FOR EIGENVECTOR OUTPUT
        iflag = sign(1, nroot)
        nroot = abs(nroot)
        !     FOLLOWING COORDINATE DEFINITIONS PROVIDE FOR PROPER SYSTEM
        do ix=1, 3
            r(ix) = q(indx(ix))
        enddo
        if( idebug .eq. 1 ) write(iptprt, '(1X,''DIMIMA, Q='',1P3E13.5,/,9X,'' R='',3E13.5)') q, r
        !  CHECK FOR NEGATIVE R'S
        xter  = 1.e30_wp
        dq(1) = -xter
        dq(2) = -xter
        dq(3) = -xter
        if( r(1) .lt. 0.e0_wp .or. r(2) .lt. 0.e0_wp .or. r(3) .lt. 0.e0_wp ) return
        !     ZERO OUT UNUSED HAMILTONIAN ELEMENTS
        th(2, 3) = zero
        th(3, 5) = zero
        th(3, 6) = zero
        !     COMPUTE SINE AND COSINE OF ANGLE AXB, AND POWERS
        ct = (r(1)*r(1)+r(3)*r(3)-r(2)*r(2))/(two*r(1)*r(3))
        goto_block_10_466: block
            if( ct .gt. -test .and. ct .lt. test ) exit goto_block_10_466
            ct = sign(test, ct)
        end block goto_block_10_466
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
        !     CONSTRUCT DIM HAMILTONIAN IN ORTHOGONAL BASIS
        th(1, 1) = ss11a+forth*((ssab+tre*tsab)+ct2*(ss11b+tre*tsxb) +st2*(spxb+tre*tpxb))-dpux-two*dsgh
        th(1, 2) = forth*sqrt3*(ssab-tsab-ct2*(ss11b-tsxb)-st2*(spxb-tpxb))
        th(1, 3) = ss12a
        th(1, 4) = -half*ct*ss12b
        th(1, 5) = forth*cs*(ss11b+tre*tsxb-spxb-tre*tpxb)
        th(1, 6) = -forth*sqrt3*cs*(ss11b-tsxb-spxb+tpxb)
        th(2, 2) = tsxa+forth*((tre*ssab+tsab)+ct2*(tre*ss11b+tsxb) +st2*(tre*spxb+tpxb))-dpux-two*dsgh
        th(2, 4) = half*sqrt3*ct*ss12b
        th(2, 5) = th(1, 6)
        th(2, 6) = forth*cs*(tre*ss11b+tsxb-tre*spxb-tpxb)
        th(3, 3) = ss22a+half*(dgab+duab)+dsxb-ssgx-ssgh-dsgh
        th(3, 4) = half*(dgab-duab)
        th(4, 4) = dsxa+half*(dgab+duab)+ss22b-ssgx-ssgh-dsgh
        th(4, 5) = -half*st*ss12b
        th(4, 6) = half*sqrt3*st*ss12b
        th(5, 5) = spxa+forth*((ssab+tre*tsab)+st2*(ss11b+tre*tsxb) +ct2*(spxb+tre*tpxb))-dpux-two*dsgh
        th(5, 6) = forth*sqrt3*(ssab-tsab-st2*(ss11b-tsxb)-ct2*(spxb-tpxb))
        th(6, 6) = tpxa+forth*((tre*ssab+tsab)+st2*(tre*ss11b+tsxb) +ct2*(tre*spxb+tpxb))-dpux-two*dsgh
        !     FILL OUT MATRIX AND DIAGONALIZE
        do i=1, ndim
            do j=1, i
                th(i, j) = th(j, i)
            enddo
        enddo
        call eign(ndim, th, vec, en, w, ndim)
        !     FOLLOWING STATEMENT ASSUMES ENERGY ZERO AT A + BC(R=RE)
        xter = en(nroot)+c(2)
        if( idebug .eq. 1 ) &
            & write(iptprt, '(1X,''NROOT,EN(1),EN(NROOT),C(2),XTER='',I5,1P4E13.5)') nroot, en(1), en(nroot), c(2), xter
        if( iflag .lt. zero ) write(iptprt, '(1X,''VECTOR='',6D20.8)') (vec(i,nroot), i=1, ndim)
        !************************* END OF FPOT3N *******************************
        if( nder .ne. 1 ) return
        !     SET  VECTOR FOR STATE OF INTEREST
        do i=1, ndim
            gvec(i) = vec(i,nroot)
        enddo
        !     COMPUTE PARTIALS OF U=COS(THETA)
        dtdr1 = (r(1)*r(1)+r(2)*r(2)-r(3)*r(3))/(two*r(1)*r(1)*r(3))
        dtdr2 = -r(2)/(r(1)*r(3))
        dtdr3 = (r(2)*r(2)+r(3)*r(3)-r(1)*r(1))/(two*r(3)*r(3)*r(1))
        factr = -one/st
        if( abs(st) .lt. 1.e-6_wp ) write(iptprt, &
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
        !     LOOP OVER R(I) DERIVATIVES
        !     FIRST ZERO OUT UNUSED ELEMENTS
        dh(2,3) = zero
        dh(3,5) = zero
        dh(3,6) = zero
        do i=1, 3
            !     SET UP DERIVATIVE MATRIX
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
            !     FILL OUT DERIVATIVE MATRIX
            do j=1, ndim
                do k=1, j
                    dh(j,k) = dh(k,j)
                enddo
            enddo
            dnumb = zero
            do j=1, ndim
                do k=1, ndim
                    dnumb = dnumb+gvec(j)*dh(j,k)*gvec(k)
                enddo
            enddo
            dd(i) = dnumb
            !     END OF DERIVATIVE LOOP
        enddo
        !     PUT DERIVATIVES INTO CORRESPONDING SPOTS FOR CORRECT COORD CHOICE.
        do ix=1, 3
            dq(indx(ix)) = dd(ix)
        enddo
        !************************ END OF DFPOT *********************************

        contains

        !> @brief Evaluate SETSPL for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine setspl
            !     SETS UP SPLINE FITS FOR H2 CURVES
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: epsln, ss, ss1, ss2, t
            integer  :: i
            !     NOW COMPUTE SPLINE FIT DATA
            !     H2 SINGLET STATE
            do i=1, 67
                r1h2(i) = real(i+3, wp)*0.1e0_wp
            enddo
            r1h2(68) = 7.2e-00_wp
            r1h2(69) = 7.4e-00_wp
            r1h2(70) = 7.6e-00_wp
            r1h2(71) = 7.8e-00_wp
            r1h2(72) = 8.0e-00_wp
            r1h2(73) = 8.25e-00_wp
            r1h2(74) = 8.5e-00_wp
            r1h2(75) = 9.0e-00_wp
            r1h2(76) = 9.5e-00_wp
            r1h2(77) = 10.e-00_wp
            r1h2(78) = 11.e-00_wp
            r1h2(79) = 12.e-00_wp
            do i=1, 79
                v1h2(i) = v1h2(i)+one
            enddo
            s2ss(1)  = 3.2006875e+1_wp
            s2ss(79) = -7.93478e-7_wp
            call splnin(79, r1h2, v1h2, t, ss, ss1, ss2, epsln, s2ss, s3ss, delss)
            !     H2 TRIPLET STATE
            do i=1, 81
                r2h2(i) = real(i+9, wp)*0.1e-00_wp
            enddo
            r2h2(82) = 9.25e-00_wp
            r2h2(83) = 9.5e-00_wp
            r2h2(84) = 9.75e-00_wp
            r2h2(85) = 10.e-00_wp
            r2h2(86) = 11.e-00_wp
            r2h2(87) = 12.e-00_wp
            do i=1, 87
                v2h2(i) = v2h2(i)+one
            enddo
            s2ts(1)  = 1.3039390e0_wp
            s2ts(87) = -7.93478e-7_wp
            call splnin(87, r2h2, v2h2, t, ss, ss1, ss2, epsln, s2ts, s3ts, delts)
            !     H2+ G STATE
            do i=1, 100
                r3h2(i) = real(i, wp)*0.1e-00_wp
            enddo
            do i=1, 26
                r3h2(i+100) = real(i, wp)*half+10.e-00_wp
            enddo
            do i=1, 126
                v3h2(i) = v3h2(i)+one
            enddo
            s2dg(1)   = 1.0742890e+3_wp
            s2dg(126) = -3.03938e-7_wp
            call splnin(126, r3h2, v3h2, t, ss, ss1, ss2, epsln, s2dg, s3dg, deldg)
            !     H2+ U STATE
            do i=1, 126
                v4h2(i) = v4h2(i)+one
            enddo
            s2du(1)   = 9.0876397e+2_wp
            s2du(126) = -3.0398e-7_wp
            call splnin(126, r3h2, v4h2, t, ss, ss1, ss2, epsln, s2du, s3du, deldu)
        end subroutine setspl

        !> @brief Evaluate HCLM for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine hclm( r, xhm2s, dxhm2s )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a, a1, arg, b, b1, c, c1, dx, dxhm2s, dy, r, rsmall, x, xhm2s, y
            !     DOUBLET SIGMA STATE OF HCL-IN A.U.
            a      = 2.52966e+1_wp
            b      = -2.23e+0_wp
            c      = -1.32774e-1_wp
            a1     = 2.4913064e+4_wp
            b1     = 7.1248778e+0_wp
            c1     = 2.8830121
            rsmall = 0.94486e+0_wp
            if( r < rsmall )then
                arg    = a1*exp(-b1*r)
                xhm2s  = arg+c1
                dxhm2s = -b1*arg
            else
                dx     = one/r**7
                x      = dx*r
                dy     = x*r
                y      = dy*r
                xhm2s  = a*x+b*y+c
                dxhm2s = -six*a*dx-for*b*dy
            endif
        end subroutine hclm

        !> @brief Evaluate HCL for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine hcl( r, h, xh3s, xh1p, xh3p, dh, dxh3s, dxh1p, dxh3p )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a2, aa, ab, alp, ang, angcon, arg, b0, b2, ba, bb, bet, bt, bv, bx, c, c0, c2, c4, c6, c8, cc, cs, &
                & cv, d0, da, da2, dab, dalp, dang, db2, dbet, dc, dc2, dcs, de1, de12, de2, dfac1, dfac2, dg, dg11, dg12, dg22, &
                & dh(2,2), dh11, dh12, dh22, dovlp, ds, ds12, ds2, dt, dt1, dt2, dv, dx, dxh1p, dxh3p, dxh3s, e1, e12, e2, ex, ex1, &
                & ex2, fac1, fac2, g, g11, g12, g22, gb, gmax, h(2,2), h11, h12, h22, omex, opex, ovlp, r, r0, r2, ra, ri, ri2, ri4, &
                & ri6, ri8, rmax, root, rops2, rt, rv, rv1, rx, s, s12, s2, t1, t2, ttt, xh1p, xh3p, xh3s, y, z1, z2, z3, z4, z5, &
                & z6, zeta
            !     CALCULATES 2X2 SINGLET SIGMA BLOCK, AND TRIPLET SIGMA, SINGLET
            !     PI, AND TRIPLET PI CURVES FOR HCL IN A.U.
            !     TRIPLET SIGMA DE ALREADY CONTAINS SATO=0.06
            !     ANGCON = +1 FOR #2 PERP TO #1, = -1 FOR #1 PERP TO #2,
            !        AND = 0 FOR SYMMETRIC ORTHOGONALIZATION
            !     SINGLET SIGMA BLOCK
            !     GROUND STATE ADIABAT
            dx     = 1.69552e-1_wp
            rx     = 2.4082e+0_wp
            bx     = 9.89104e-1_wp
            dv     = 1.84605e-1_wp
            rv     = 4.746e+0_wp
            bv     = 2.77985e-1_wp
            cv     = 1.82621e-1_wp
            c8     = -2.018460e+4_wp
            c6     = 2.399677e+3_wp
            c4     = -5.351085e+1_wp
            c0     = 3.67226e-1_wp
            aa     = 25.20626238e-0_wp
            bb     = 1.910495693e-0_wp
            cc     = 1.828090476e-1_wp
            rt     = 2.53e-0_wp
            bt     = 1.047e-0_wp
            dt     = 6.685927358e-2_wp
            ra     = 2.274408e+0_wp
            ba     = 1.069975e+0_wp
            da     = 5.3556e-2_wp
            rmax   = 2.5e+0_wp
            gmax   = 2.0212e-2_wp
            gb     = 4.654e-1_wp
            rv1    = 3.5e-0_wp
            angcon = 1.0e-0_wp
            zeta   = 2.0387e0_wp
            z1     = 7.222580119e+2_wp
            z2     = -2.749954467e+2_wp
            z3     = 8.648352007e+2_wp
            z4     = 2.437196153e+2_wp
            z5     = 3.144371653e+1_wp
            z6     = 1.128252224e-1_wp
            r0     = 2.53e-0_wp
            b0     = 0.9897e-0_wp
            d0     = 0.150789e-0_wp
            ex     = exp(-bx*(r-rx))
            omex   = one-ex
            e1     = dx*omex*omex-dx
            de1    = two*dx*bx*ex*omex
            !     EXCITED STATE ADIABAT
            if( r > rv )then
                ri  = one/r
                ri2 = ri*ri
                ri4 = ri2*ri2
                ri6 = ri4*ri2
                ri8 = ri4*ri4
                e2  = c8*ri8+c6*ri6+c4*ri4-ri+c0
                de2 = -ri*(eht*c8*ri8+six*c6*ri6+for*c4*ri4-ri)
            elseif( r < rv1 )then
                arg = aa*exp(-bb*r)
                e2  = arg+cc
                de2 = -bb*arg
            else
                ex   = exp(-bv*(r-rv))
                omex = one-ex
                e2   = dv*omex*omex+cv
                de2  = two*dv*bv*ex*omex
            endif
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
            if( abs(s12) .gt. 1.0e0_wp ) write(iptprt, '(1X,''In the DIM potential in the subprogram HCL: '', ''R,S12='',F10.5,1PE1' // &
                '5.7)') r, s12
            ds12 = dovlp*rops2**3/rt2i
            !     H12 AND H22 DIABATS IN NON-ORTHOGONAL REPRESENTATION
            e12   = e1+e2
            de12  = de1+de2
            fac1  = e12*h(1,1)-e1*e2-h(1,1)*h(1,1)
            fac2  = one-s12*s12
            dfac2 = -two*s12*ds12
            ttt   = fac1*fac2
            root  = 0.e0_wp
            if( ttt .gt. 0.e0_wp ) root=sqrt(ttt)
            h(1,2)  = s12*h(1,1)+root
            dh(1,2) = ds12*h(1,1)+s12*dh(1,1)
            ttt                          = half*(dfac2*fac1+fac2*(de12*h(1,1) +e12*dh(1,1)-e1*de2-de1*e2-two* h(1,1)* &
                & dh(1,1)))
            if( root .ne. 0.e0_wp ) dh(1,2) = dh(1,2) + ttt/root
            h(2,2)                          = e12*fac2+two*s12*h(1,2)-h(1,1)
            dh(2,2)                         = de12*fac2+e12*dfac2+two*(ds12*h(1,2)+s12*dh(1,2)) -dh(1,1)
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
            if( abs(y) .gt. 1.0e0_wp ) y = sign(1.0e0_wp, y)
            ang  = angcon*acos(y)
            dang = zero
            if( abs(ang) >= 1.0e-14_wp ) &
                dang = half*angcon*(dt1/(t1*t1)+dt2/(t2*t2))/sqrt(one-y*y)
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
            dh(2,1) = dh12
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
            if( r > rmax )then
                g  = gmax*exp(-gb*(r-rmax)**2)
                dg = -two*gb*(r-rmax)*g
            endif
            xh3p  = xh1p-g
            dxh3p = dxh1p-dg
        end subroutine hcl

        !> @brief Evaluate H2CURV for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine h2curv( r, ab1s, ab3s, ab2g, ab2u, dab1s, dab3s, dab2g, dab2u )
            !     CALCULATES ENERGY AND DERIVATIVES FOR H2 SINGLET AND TRIPLET
            !     SIGMA AND H2+ DOUBLET SIGMA G AND U STATES IN A.U.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a1, ab1s, ab2g, ab2u, ab3s, arg, arg2, arg3, b1, be, bfac, c10, c4, c6, c8, dab1s, dab2g, &
                & dab2u, dab3s, de, df1, epsln, ex, f1, f2, opex, r, r10, r2, r4, r6, r8, re, ss2
            !     H2 EXPONENTIALS W.R.T. FINAL LIMITS, H2+ W.R.T. 0.0 LIMIT
            !**   DATA A2,B2,C2/2.41893D00,2.1329115D00,0.0918539D00/
            !     H2 CURVES
            !     SINGLET SIGMA
            b1  = 6.0315697e-00_wp
            a1  = 9.8214033e-00_wp
            c6  = 6.49903e-00_wp
            c8  = 124.399e-00_wp
            c10 = 1135.21e-00_wp
            de  = 3.21533e-2_wp
            be  = 8.3e-1_wp
            re  = 1.60e-0_wp
            c4  = 2.25e0_wp
            if( r > 12.0e0_wp )then
                r2    = one/(r*r)
                r6    = r2**3
                r8    = r6*r2
                r10   = r8*r2
                ab1s  = -(c6*r6+c8*r8+c10*r10)
                dab1s = (six*c6*r6+eht*c8*r8+10.0e0_wp*c10*r10)/r
            elseif( r >= 0.4e0_wp )then
                call spline(79, r1h2, v1h2, r, ab1s, dab1s, ss2, epsln, s2ss, s3ss, delss)
            else
                ab1s  = a1*exp(-b1*r)
                dab1s = -b1*ab1s
            endif
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
            bfac = 2.6e-2_wp
            if( r .lt. re ) bfac=2.026e0_wp
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
            if( r > 23.0e0_wp )then
                r4    = one/r**4
                ab2g  = -c4*r4+half
                dab2g = 9.0e0_wp*r4/r
                ab2u  = ab2g
                dab2u = dab2g
            else
                call spline(126, r3h2, v3h2, r, ab2g, dab2g, ss2, epsln, s2dg, s3dg, deldg)
                call spline(126, r3h2, v4h2, r, ab2u, dab2u, ss2, epsln, s2du, s3du, deldu)
            endif
        end subroutine h2curv

        !> @brief Evaluate SPLNIN for the ClH2_LEPS_TTGI1_1985 potential.
        subroutine splnin( n, x, y, t, ss, ss1, ss2, epsln, s2, s3, dely )
            !     JIM STINE'S INITIALIZATION ROUTINE FOR SPLINE FITTING
            !     N=NO. POINTS (X,Y); T=ARGUMENT; SS=VALUE; SS1=1ST DERIV;
            !     SS2=2ND DERIV; S2,S3, AND DELY ARE COEFFICIENT ARRAYS
            !     SECOND DERIV'S AT ENDS STORED IN S2(1) AND S2(N) ON INPUT
            !     EPSLN IS OVERALL FIT CRITERION
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: b(126), c(126), delsqy, dely(n), epsln, eta, ex, h(126), h2(126), half, omega, s2(n), s3(n), ss, &
                & ss1, ss2, t, w, x(n), y(n)
            integer :: i, knt, n, n1
            ex    = 1.1e0_wp
            omega = 1.0717968e0_wp
            half  = 0.5e0_wp
            epsln = 1.e-12_wp
            n1    = n-1
            knt   = 0
            do i=1, n1
                h(i)    = x(i+1)-x(i)
                dely(i) = (y(i+1)-y(i))/h(i)
            enddo
            do i=2, n1
                h2(i)  = h(i-1)+h(i)
                b(i)   = half*h(i-1)/h2(i)
                delsqy = (dely(i)-dely(i-1))/h2(i)
                s2(i)  = delsqy+delsqy
                c(i)   = delsqy+s2(i)
            enddo
            do
                eta = 0.0e0_wp
                do i=2, n1
                    w     = (c(i)-b(i)*s2(i-1)-(half-b(i))*s2(i+1)-s2(i))*omega
                    eta   = max(eta, abs(w))
                    s2(i) = s2(i)+w
                enddo
                knt = knt+1
                if( knt .gt. 10 ) epsln=epsln*ex
                if( eta < epsln ) exit
            enddo
            do i=1, n1
                s3(i) = (s2(i+1)-s2(i))/h(i)
            enddo
            if( idebug .eq. 1 ) write(iptprt, '(5X,"Final EPS =",D15.6)') epsln
        end subroutine splnin

        !> @brief Interpolate a tabulated curve and its first two derivatives.
        subroutine spline( n, x, y, t, ss, ss1, ss2, epsln, s2, s3, dely )
            use chempotpy_kinds, only: wp
            implicit none
            integer, intent(in)   :: n
            real(wp), intent(in)  :: x(n), y(n), t, epsln, s2(n), s3(n), dely(n)
            real(wp), intent(out) :: ss, ss1, ss2
            integer               :: i
            real(wp)              :: delsqs, ht1, ht2, prod
            if( t <= x(1) )then
                ss  = y(1)
                ss1 = dely(1)
                ss2 = s2(1)
                return
            elseif( t >= x(n) )then
                ss  = y(n)
                ss1 = dely(n)
                ss2 = s2(n)
                return
            endif
            i = 1
            do while( t > x(i) )
                i = i+1
            enddo
            if( t < x(i) ) i = i-1
            ht1    = t-x(i)
            ht2    = t-x(i+1)
            prod   = ht1*ht2
            ss2    = s2(i)+ht1*s3(i)
            delsqs = (s2(i)+s2(i+1)+ss2)/6.0e0_wp
            ss     = y(i)+ht1*dely(i)+prod*delsqs
            ss1    = dely(i)+(ht1+ht2)*delsqs+prod*s3(i)/6.0e0_wp
        end subroutine spline

        !> @brief Evaluate EIGN for the ClH2_LEPS_TTGI1_1985 potential.
        !> @brief Diagonalize a real symmetric matrix with the shared core eigensolver.
        subroutine eign( nn, a, vec, eig, w, nd )
            use chempotpy_kinds, only: wp
            implicit none
            integer, intent(in)     :: nn, nd
            real(wp), intent(inout) :: a(nd,nd)
            real(wp), intent(out)   :: vec(nd,nd), eig(nd), w(nd)
            integer                 :: ierr, lwork
            real(wp), allocatable   :: work(:)
            external                         :: dsyev
            vec = 0.0e0_wp
            eig = 0.0e0_wp
            w   = 0.0e0_wp
            if( nn == 0 ) return
            lwork = max(1, 3*nn-1)
            allocate(work(lwork))
            call dsyev('V', 'L', nn, a, nd, eig, work, lwork, ierr)
            if( ierr /= 0 ) error stop 'EIGN: symmetric eigensolver failed to converge'
            vec(1:nn,1:nn) = a(1:nn,1:nn)
            deallocate(work)
        end subroutine eign
    end subroutine evaluate_internal
end module chempotpy_clh2_leps_ttgi1_1985_pot_legacy_state_bf7e4ade

!> @brief Evaluate TTGI1 for one set of internal distances.
subroutine pot( q, nder, xter, dq )
    use chempotpy_clh2_leps_ttgi1_1985_pot_legacy_state_bf7e4ade, only: evaluate_internal
    use chempotpy_kinds,                                          only: wp
    implicit none
    real(wp), intent(in)  :: q(3)
    integer, intent(in)   :: nder
    real(wp), intent(out) :: xter, dq(3)
    call evaluate_internal(q, nder, xter, dq)
end subroutine pot

!> @brief Evaluate POTINFO for the ClH2_LEPS_TTGI1_1985 potential.
