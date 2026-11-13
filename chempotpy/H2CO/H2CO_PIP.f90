module shell
    use chempotpy_h2co_pip_parameters, only: cof => fitted_coefficients, power => powers
    use chempotpy_kinds,               only: wp
    implicit none
    private
    public :: cof, coff, power, total
    real(wp), parameter :: total = 8.0_wp
    integer, parameter  :: coff  = 1561
end module shell

!> @brief Evaluate getpot.
subroutine getpot( cood, v )
    use shell,           only: cof, coff, power, total
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: i
    real(wp) :: v, f
    real(wp) :: cood(3,4), com(3)
    real(wp) :: y(6), r(6)
    real(wp) :: bas(coff)
    real(wp) :: smallch, smalloh
    real(wp) :: s, x, z, x3, x2, z1, z2
    !       This probram expects cood to be in the order C,O,H,H
    call edis(r(1), y(1), cood(:,3), cood(:,4))
    call edis(r(2), y(2), cood(:,2), cood(:,3))
    call edis(r(3), y(3), cood(:,1), cood(:,3))
    call edis(r(4), y(4), cood(:,1), cood(:,4))
    call edis(r(5), y(5), cood(:,2), cood(:,4))
    call edis(r(6), y(6), cood(:,1), cood(:,2))
    call basis(bas, y)
    if( r(3) < r(4) )then
        smallch = r(3)
        smalloh = r(2)
    else
        smallch = r(4)
        smalloh = r(5)
    endif
    ! Switching function has been used here |
    x  = r(1)
    z  = x-8
    z1 = 10-8
    z2 = z/z1
    x2 = z2*z2
    x3 = x2*z2
    if( x .ge. 10 )then
        s = 1.0
    elseif( x .le. 8 )then
        s = 0.0e0_wp
    else
        s = 10*x3-15*x3*z2+6*x3*x2
    endif
    f = 0
    do i=1, coff
        f = f+cof(i)*bas(i)
    enddo
    f = f+114.332958863-1.059892782251382e-004
    if( x .ge. 8 )then
        call hcopot(v, smallch, r(6), smalloh)
        f = f*(1-s)+v*s
    endif
    v = f
end subroutine getpot

!> @brief Evaluate edis.
subroutine edis( dis, edist, atomi, atomj )
    use shell,           only: cof, coff, power, total
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), intent(out) :: edist
    real(wp)              :: atomi(3), atomj(3)
    real(wp)              :: dis, alfa
    alfa  = 2.0
    dis   = sqrt((atomi(1)-atomj(1))**2+(atomi(2)-atomj(2))**2+(atomi(3)-atomj(3))**2)
    edist = exp(-dis/alfa)
end subroutine edis

!> @brief Evaluate basis.
subroutine basis( bas, y )
    use shell,           only: cof, coff, power, total
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: m, n, p, q, i, j, k
    real(wp) :: bas(coff)
    real(wp) :: y(6)
    do k=1, coff
        bas(k)=(y(1)**power(k, 1))*(y(6)**power(k, 6))*((y(2)**power(k, 2))* (y(3)**power(k, 3))*(y(4)**power(k, 4))* &
            & (y(5)** power(k, 5))+(y(2)**power(k, 5))* (y(3)**power(k, 4))*(y(4)**power(k, 3))*(y(5)**power(k, 2)))
    enddo
end subroutine

!> @brief Evaluate hcopot.
subroutine hcopot( f, rhc, rco, roh )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp), parameter :: au2ev = 27.2116e0_wp
    real(wp), parameter :: a0    = 0.52918e0_wp
    real(wp), parameter :: pi    = 3.1415926536e0_wp
    integer             :: j, loop, iv
    real(wp)            :: rmin(3), rmax(3), z1, z2, z3, m(3), h(3)
    real(wp)            :: pot, potmin, ma, mb, mc
    real(wp)            :: pot3, c3, pothco, f
    real(wp)            :: v, rhc, rco, cgamma, roh
    integer             :: istate
    ma     = 1.0    !  H
    mb     = 12.0    !  C
    mc     = 16.0    !  O
    istate = 1      !  X-Flaeche
    cgamma = (rhc*rhc+rco*rco-roh*roh)/(2*rhc*rco)
    pot    = 0.0e0_wp
    pot    = pothco(rhc, rco, cgamma, istate)
    pot    = pot+.18140000
    f      = pot
end subroutine hcopot

!> @brief Evaluate POTHCO.
function pothco( rhc, rco, rcos, istate )
    use chempotpy_h2co_pip_parameters, only: parameter_aa => aa, parameter_easym => easym, &
        & parameter_gnhc => gnhc, parameter_gnco => gnco, parameter_gnwin => gnwin, &
        & parameter_gwhc => gwhc, parameter_gwco => gwco, parameter_gwwin => gwwin, &
        & parameter_gn2co => gn2co, parameter_gn2hc => gn2hc, parameter_gn2oh => gn2oh, &
        & parameter_gzhc => gzhc, parameter_gzco => gzco, parameter_gzwin => gzwin, &
        & parameter_gz2co => gz2co, parameter_gz2hc => gz2hc, parameter_gz2oh => gz2oh, &
        & parameter_ixco => ixco, parameter_ixhc => ixhc, parameter_ixoh => ixoh, &
        & parameter_ix2co => ix2co, parameter_ix2hc => ix2hc, parameter_ix2oh => ix2oh, &
        & parameter_kk => kk, parameter_nenadd => nenadd, parameter_wuradd => wuradd
    use chempotpy_kinds, only: wp
    implicit none
    integer  :: istate ! Initialiserung falls <>0, INPUT
    real(wp) :: rhc, rco, rcos, pothco
    ! INPUT: R_HC in a.u., R_CO in a.u.,  COS(gamma_HCO),
    ! OUTPUT: Potentialwert in a.u.
    real(wp), parameter :: au2ang = 0.529177e0_wp  ! a_0 / Angstroem
    real(wp), parameter :: au2ev  = 27.2116e0_wp  ! hartree / eV
    real(wp), parameter :: ev2cm  = 8065.44e0_wp  ! eV / cm^-1
    real(wp), parameter :: kc2ev  = 23.061e0_wp  ! {kcal/mol} / eV
    real(wp), parameter :: au2cm  = au2ev*ev2cm  ! hartree / cm^-1
    integer, parameter  :: nhc    = 12, nco=9, nwin=12
    integer             :: gzhc, gzco, gzwin, mlin(2) ! GRAD IM ZAEHLER
    integer             :: gnhc, gnco, gnwin, mtan(2) ! GRAD IM NENNER
    integer             :: ixhc, ixco, ixoh ! FESTE VORFAKTOREN IM NENNER.
    integer             :: gwhc, gwco, gwwin ! POLYNOME VOR WURZEL (KONISCHE D.)
    integer             :: gz2co, gn2co, ix2co ! DIATOM-POTENTIAL v. C-O.
    integer             :: gz2hc, gn2hc, ix2hc ! DIATOM-POTENTIAL v. H-C.
    integer             :: gz2oh, gn2oh, ix2oh ! DIATOM-POTENTIAL v. O-H.
    integer, parameter  :: maxlin = 500
    integer, parameter  :: maxtan = 200  ! Maximalwerte fuer MLIN und MTAN.
    real(wp)            :: wuradd, nenadd, easym
    real(wp), parameter :: nuhc(12) = (/2.3e0_wp, 2.5e0_wp, 2.1e0_wp, 3.0e0_wp, 1.9e0_wp, 3.5e0_wp, 1.7e0_wp, &
        & 4.0e0_wp, 1.5e0_wp, 5.0e0_wp, 1.3e0_wp, 7.0e0_wp/)
    real(wp), parameter :: nuco(9) = (/2.2e0_wp, 1.9e0_wp, 2.5e0_wp, 1.7e0_wp, 3.0e0_wp, 2.05e0_wp, 2.35e0_wp, &
        & 3.5e0_wp, 4.0e0_wp/)
    real(wp), parameter :: nuwin(12) = (/100e0_wp, 120e0_wp, 80e0_wp, 140e0_wp, 60e0_wp, 160e0_wp, &    ! NEUE NULLST.
        &40e0_wp, 170e0_wp, 20e0_wp, 180e0_wp, 10e0_wp, 0e0_wp/)
    real(wp) :: xcos(nwin), nucos(nwin), kk(maxlin,2), aa(maxtan,2), v2ch, v2oh
    external v2ch, v2oh
    integer     :: i, j, k, l, ll, jj, j1, j2, j3, ks, i1, i2
    logical     :: qinit(2) = [.false., .false.] ! zeigt an ob Parameter initialisiert.
    real(wp)    :: pq, a, b, c, nenn, rhcsin, rhccos, konus
    real(wp)    :: phc, pco, pcos, poh, rhci, rcoi, rcosi, rohi
    real(wp)    :: thy, aco, ahc, vplus, pi
    complex(wp) :: qq1, qq2, qq3, qq4, wurz
    gzhc        = parameter_gzhc; gzco = parameter_gzco; gzwin = parameter_gzwin
    gnhc        = parameter_gnhc; gnco = parameter_gnco; gnwin = parameter_gnwin
    ixhc        = parameter_ixhc; ixco = parameter_ixco; ixoh = parameter_ixoh
    gwhc        = parameter_gwhc; gwco = parameter_gwco; gwwin = parameter_gwwin
    gz2hc       = parameter_gz2hc; gn2hc = parameter_gn2hc; ix2hc = parameter_ix2hc
    gz2co       = parameter_gz2co; gn2co = parameter_gn2co; ix2co = parameter_ix2co
    gz2oh       = parameter_gz2oh; gn2oh = parameter_gn2oh; ix2oh = parameter_ix2oh
    wuradd      = parameter_wuradd; nenadd = parameter_nenadd; easym = parameter_easym
    kk          = 0.e0_wp; aa = 0.e0_wp
    kk(1:230,:) = parameter_kk; aa(1:132,:) = parameter_aa
    mlin(1)     = (gzhc+1)*(gzco+1)*(gzwin+1)+(gz2co+1)+(gwhc+1)*(gwco+1)*(gwwin+1)
    mlin(2)     = (gzhc+1)*(gzco+1)*(gzwin+1)+(gz2co+1)
    mtan(1)     = 2*(gnhc+1)*(gnco+1)*(gnwin+1)-2+2*gn2co+4
    mtan(2)     = 2*(gnhc+1)*(gnco+1)*(gnwin+1)-2+2*gn2co
    pi          = acos(-1.e0_wp)
    do j3=1, nwin
        nucos(j3) = cos(pi/180.e0_wp*nuwin(j3))
    enddo
    if( istate .le. 0 )then
        pothco = 0.e0_wp
        return
    else
        rhci  = rhc                                      ! interne Kopien.
        rcoi  = rco
        rcosi = rcos
        vplus = 0e0_wp                                     ! Potential-Offset.
        thy   = tanh(1.1e0_wp*(rhci-4.4e0_wp))                  ! innen -1, aussen 1.
        a     = 0.75e0_wp - 0.25e0_wp * thy                       ! innen 1, aussen 1/2.
        ahc   = 1e0_wp + (-0.002e0_wp + 0.014e0_wp * rcosi ) * a   ! Skalierung in R_HC.
        aco   = 1e0_wp +         0.012e0_wp           * a       ! Skalierung in R_CO.
        rhci  = ahc * (rhci-2.11e0_wp) + 2.11e0_wp
        rcoi  = aco * (rcoi-2.23e0_wp) + 2.23e0_wp
        vplus = ( 0.5e0_wp * thy - 0.5e0_wp ) * 50e-3_wp / 27.2116e0_wp ! innen -50 meV.
        rhci  = max(rhci, 1.3e0_wp)            ! R_HC > 1.3 a.u.
        rcoi  = max(rcoi, 1.5e0_wp)            ! R_CO > 1.5 a.u.
        rohi  = sqrt( rhci**2 + rcoi**2 - 2e0_wp*rhci*rcoi*rcosi )
        rohi  = max(rohi, 0.01e0_wp)          ! zur numerischen Stabilitaet.
        ks    = istate
        jj    = 1                           ! Index in AA(JJ,KS) (Nenner).
        qq2   = (1e0_wp, 0e0_wp)                   ! *** QQ_CO(R_CO)
        pco   = 1e0_wp
        do j=1, gn2co
            pco = pco * ( rcoi - nuco(j) ) ! P_{CO,J}(R_CO)
            qq2 = qq2 + cmplx( aa(jj,ks)*pco, aa(jj+1,ks)*pco, kind=wp)
            jj  = jj + 2
        enddo
        qq4 = (1e0_wp, 0e0_wp)                   ! *** QQ_HCO(R_HC,R_CO,COS())
        phc = 1e0_wp
        do j=0, gnhc
            pco = phc
            do k=0, gnco
                pcos = pco
                if( j .gt. 0.or.k .gt. 0 )then
                    qq4 = qq4 + cmplx( aa(jj,ks)*pco, aa(jj+1,ks)*pco, kind=wp)
                    jj  = jj + 2
                endif
                do l=1, gnwin
                    pcos = pcos * ( rcosi - nucos(l) )
                    qq4  = qq4 + cmplx( aa(jj,ks)*pcos, aa(jj+1,ks)*pcos, kind=wp)
                    jj   = jj + 2
                enddo
                pco = pco * ( rcoi - nuco(k+1) )
            enddo
            phc = phc * ( rhci - nuhc(j+1) )
        enddo
        if( ks .eq. 1 )then
            rhccos = rhci * rcosi              ! R_HC * cos(gamma)
            rhcsin = rhci * sqrt(1e0_wp-rcosi**2) ! R_HC * sin(gamma)
            a      = 1e0_wp + aa(jj+1,ks)*rhccos&                 ! 1 + a2*RHCCOS +
                &+( aa(jj+2,ks) + aa(jj+3,ks)*rcoi ) * rcoi ! + a3*RCOI+a4*RCOI**2
            wurz = sqrt( a**2 + (aa(jj,ks)*rhcsin)**2 + wuradd ) ! WURADD ist
            jj   = jj + 4
        endif
        if( jj .ne. mtan(ks)+1 ) stop 'POTHCO: MTAN FALSCH !'
        ll  = 1
        pq  = easym
        pq  = pq + v2ch(rhci)
        pco = rcoi**ix2co / abs(qq2)**2  ! *** PP_CO(R_CO)
        do j=0, gz2co
            pq  = pq + kk(ll,ks) * pco
            ll  = ll + 1
            pco = pco * ( rcoi - nuco(j+1) ) ! P_{CO,J}(R_CO)
        enddo
        pq = pq + v2oh(rohi)
        nenn= rhci**(ixhc-ixoh) * rcoi**(ixco-ixoh) * rohi**ixoh&
            &/ ( abs(qq4)**2 + nenadd ) ! *** PP_HCO(R_HC,R_CO,COS())
        konus = 0e0_wp                      ! Konische Durchschneidung:
        if( ks .eq. 1 )then! (nur auf 1. Flaeche)
            phc = wurz * nenn
            do j=0, gwhc
                pco = phc
                do k=0, gwco
                    pcos = pco
                    do l=0, gwwin
                        konus = konus + kk(ll,ks) * pcos
                        ll    = ll + 1
                        pcos  = pcos * ( rcosi - nucos(l+1) )
                    enddo
                    pco = pco * ( rcoi - nuco(k+1) )
                enddo
                phc = phc * ( rhci - nuhc(j+1) )
            enddo
            pq = pq + konus
        endif
        phc = nenn
        do j=0, gzhc                     ! normale Polynome:
            pco = phc
            do k=0, gzco
                pcos = pco
                do l=0, gzwin
                    pq   = pq + kk(ll,ks) * pcos
                    ll   = ll + 1
                    pcos = pcos * ( rcosi - nucos(l+1) )
                enddo
                pco = pco * ( rcoi - nuco(k+1) )
            enddo
            phc = phc * ( rhci - nuhc(j+1) )
        enddo
        if( ll .ne. mlin(ks)+1 ) stop ' POTHCO: MLIN FALSCH '
        pq = pq + 0.055e0_wp/au2ev + vplus         ! OFFSET(Dez.1994)
        if (rohi.lt.1.1e0_wp.or.rhci.lt.1.3e0_wp.or.&
            &rcoi.lt.1.5e0_wp.or.rcoi+rohi.lt.2.9e0_wp) then
        pq = max( 6e0_wp/au2ev, pq )      ! >= 6 eV.
    endif
    pothco = min( 20e0_wp/au2ev, pq )   ! nach oben begrenzen (20eV)
    endif
end

!> @brief Evaluate V2OH.
function v2oh( r )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp)            :: v2oh, r
    real(wp), parameter :: au2ang = 0.529177e0_wp  ! a_0 / ANGSTROEM
    real(wp), parameter :: au2ev  = 27.2116e0_wp  ! HARTREE / EV
    real(wp)            :: a1, a2, a3, de, re, rho
    a1   = 4.507                    ! Huxley etc. nach Daten aus Herzberg.
    a2   = 4.884
    a3   = 3.795
    de   = 4.621 / au2ev
    re   = 0.9696
    rho  = ( r*au2ang - re )
    v2oh = -de*(1+a1*rho+a2*rho**2+a3*rho**3)*exp(-a1*rho)
end

!> @brief Evaluate V2CH.
function v2ch( r )
    use chempotpy_kinds, only: wp
    implicit none
    real(wp)            :: v2ch, r
    real(wp), parameter :: au2ang = 0.529177e0_wp  ! a_0 / ANGSTROEM
    real(wp), parameter :: au2ev  = 27.2116e0_wp  ! HARTREE / EV
    real(wp)            :: a1, a2, a3, de, re, rho
    a1   = 3.836
    a2   = 3.511
    a3   = 2.268
    de   = 3.631 / au2ev
    re   = 1.1199
    rho  = ( r*au2ang - re )
    v2ch = -de*(1+a1*rho+a2*rho**2+a3*rho**3 )*exp(-a1*rho)
end

!> @brief Implement H2CO_PIP, a potential energy surface for H2CO.
module chempotpy_h2co_pip_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H2CO_PIP potential energy surface for H2CO.
    subroutine evaluate( x, igrad, p, gradient_out, d )
        use shell,           only: cof, coff, power, total
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 4
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates)
        real(wp), intent(out) :: gradient_out(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: tx(3,natoms), v
        integer               :: iatom, idir, j, istate
        gradient_out = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        if( igrad == 0 )then
            call getpot(tx, v)
            v = v*hartree_to_electronvolt
            do istate=1, nstates
                p(istate) = v
            enddo
        else
            write (*, *) 'Only energy is available'
        endif
    end subroutine evaluate
end module chempotpy_h2co_pip_pes

!> @brief Expose H2CO_PIP through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, gradient_out, d )
    use chempotpy_h2co_pip_pes, only: evaluate
    use chempotpy_kinds,        only: wp
    implicit none
    integer, parameter              :: natoms  = 4
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: gradient_out(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, gradient_out, d)
end subroutine pes
