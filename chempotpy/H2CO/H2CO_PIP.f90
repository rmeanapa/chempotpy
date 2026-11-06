module shell
    implicit none
    private
    public :: cof, coff, power, total
    real, parameter    :: total = 8
    integer, parameter :: coff  = 1561
    !real,parameter::alpha(4)=(/2.0,2.0,2.0,2.0/)
    double precision :: cof(coff), power(coff,6)
end module shell

!> @brief Evaluate pes_init.
subroutine pes_init( path )
    use shell, only: cof, coff, power, total
    implicit none
    character(len=1024), intent(in) :: path
    integer                         :: file_unit_112
    integer                         :: file_unit_113
    character(len=1024)             :: file_path1, file_path2
    integer                         :: i
    file_path1 = trim(path)//"/H2CO/coeff8"
    file_path2 = trim(path)//"/H2CO/index8"
    open(newunit=file_unit_112, status='old', file=file_path1)
    open(newunit=file_unit_113, status='old', file=file_path2)
    read(file_unit_112, *)
    do i=1, coff
        read(file_unit_112, *) cof(i)
    enddo
    do i=1, coff
        read(file_unit_113, *) power(i,1), power(i,2), power(i,3), power(i,4), &
            &power(i,5), power(i,6)
    enddo
    close(file_unit_113)
    close(file_unit_112)
end subroutine pes_init

!> @brief Evaluate getpot.
subroutine getpot( cood, v, path )
    use shell, only: cof, coff, power, total
    implicit none
    character(len=1024), intent(in) :: path
    integer                         :: i
    double precision                :: v, f
    double precision                :: cood(3,4), com(3)
    double precision                :: y(6), r(6)
    double precision                :: bas(coff)
    double precision                :: smallch, smalloh
    double precision                :: s, x, z, x3, x2, z1, z2
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
    else if( x .le. 8 )then
        s = 0.0d0
    else
        s = 10*x3-15*x3*z2+6*x3*x2
    endif
    f = 0
    do i=1, coff
        f = f+cof(i)*bas(i)
    enddo
    f = f+114.332958863-1.059892782251382e-004
    if( x .ge. 8 )then
        call hcopot(v, smallch, r(6), smalloh, path)
        f = f*(1-s)+v*s
    endif
    v = f
end subroutine getpot

!> @brief Evaluate edis.
subroutine edis( dis, edist, atomi, atomj )
    use shell, only: cof, coff, power, total
    implicit none
    double precision, intent(out) :: edist
    double precision              :: atomi(3), atomj(3)
    double precision              :: dis, alfa
    alfa  = 2.0
    dis   = sqrt((atomi(1)-atomj(1))**2+(atomi(2)-atomj(2))**2+(atomi(3)-atomj(3))**2)
    edist = exp(-dis/alfa)
end subroutine edis

!> @brief Evaluate basis.
subroutine basis( bas, y )
    use shell, only: cof, coff, power, total
    implicit none
    integer          :: m, n, p, q, i, j, k
    double precision :: bas(coff)
    double precision :: y(6)
    do k=1, coff
        bas(k)=(y(1)**power(k,1))*(y(6)**power(k,6))*((y(2)**power(k,2))* (y(3)**power(k,3))*(y(4)**power(k,4))* &
            & (y(5)** power(k,5))+(y(2)**power(k,5))* (y(3)**power(k,4))*(y(4)**power(k,3))*(y(5)**power(k,2)))
    enddo
end subroutine

!> @brief Evaluate hcopot.
subroutine hcopot( f, rhc, rco, roh, path )
    implicit none
    double precision, parameter     :: au2ev = 27.2116d0
    double precision, parameter     :: a0    = 0.52918d0
    double precision, parameter     :: pi    = 3.1415926536d0
    character(len=1024), intent(in) :: path
    integer                         :: j, loop, iv
    double precision                :: rmin(3), rmax(3), z1, z2, z3, m(3), h(3)
    double precision                :: pot, potmin, ma, mb, mc
    double precision                :: pot3, c3, pothco, f
    double precision                :: v, rhc, rco, cgamma, roh
    integer                         :: istate
    ma     = 1.0    !  H
    mb     = 12.0    !  C
    mc     = 16.0    !  O
    istate = 1      !  X-Flaeche
    call pot3in(path)
    cgamma = (rhc*rhc+rco*rco-roh*roh)/(2*rhc*rco)
    pot    = 0.0d0
    pot    = pothco(rhc, rco, cgamma, istate)
    pot    = pot+.18140000
    f      = pot
end subroutine hcopot

!> @brief Evaluate POT3IN.
subroutine pot3in( path )
    implicit none
    character(len=1024), intent(in) :: path
    integer                         :: file_unit_10
    character(len=1024)             :: file_path1
    double precision                :: v, pothco
    double precision, save          :: g = 0d0, r = 1d0, rkl = 1d0
    external pothco
    file_path1 = trim(path)//"/H2CO/pothco_hmk.para"
    open(newunit=file_unit_10, file=file_path1, status='old')
    v = pothco(0d0, 0d0, 0d0, 0)
    close(file_unit_10)
end

!> @brief Evaluate POTHCO.
function pothco( rhc, rco, rcos, istate )
    use iso_fortran_env, only: real64
    implicit none
    integer          :: istate ! Initialiserung falls <>0, INPUT
    double precision :: rhc, rco, rcos, pothco
    ! INPUT: R_HC in a.u., R_CO in a.u.,  COS(gamma_HCO),
    ! OUTPUT: Potentialwert in a.u.
    double precision, parameter :: au2ang = 0.529177d0  ! a_0 / Angstroem
    double precision, parameter :: au2ev  = 27.2116d0  ! hartree / eV
    double precision, parameter :: ev2cm  = 8065.44d0  ! eV / cm^-1
    double precision, parameter :: kc2ev  = 23.061d0  ! {kcal/mol} / eV
    double precision, parameter :: au2cm  = au2ev*ev2cm  ! hartree / cm^-1
    integer, parameter          :: nhc    = 12, nco=9, nwin=12
    integer                     :: gzhc, gzco, gzwin, mlin(2) ! GRAD IM ZAEHLER
    integer                     :: gnhc, gnco, gnwin, mtan(2) ! GRAD IM NENNER
    integer                     :: ixhc, ixco, ixoh ! FESTE VORFAKTOREN IM NENNER.
    integer                     :: gwhc, gwco, gwwin ! POLYNOME VOR WURZEL (KONISCHE D.)
    integer                     :: gz2co, gn2co, ix2co ! DIATOM-POTENTIAL v. C-O.
    integer                     :: gz2hc, gn2hc, ix2hc ! DIATOM-POTENTIAL v. H-C.
    integer                     :: gz2oh, gn2oh, ix2oh ! DIATOM-POTENTIAL v. O-H.
    integer, parameter          :: maxlin = 500
    integer, parameter          :: maxtan = 200  ! Maximalwerte fuer MLIN und MTAN.
    double precision            :: wuradd, nenadd, easym
    double precision, parameter :: nuhc(12) = (/2.3d0, 2.5d0, 2.1d0, 3.0d0, 1.9d0, 3.5d0, &
        &1.7d0, 4.0d0, 1.5d0, 5.0d0, 1.3d0, 7.0d0/)
    double precision, parameter :: nuco(9) = (/2.2d0, 1.9d0, 2.5d0, 1.7d0, 3.0d0, 2.05d0, &
        &2.35d0, 3.5d0, 4.0d0/)
    double precision, parameter :: nuwin(12) = (/100d0, 120d0, 80d0, 140d0, 60d0, 160d0, &    ! NEUE NULLST.
        &40d0, 170d0, 20d0, 180d0, 10d0, 0d0/)
    double precision :: xcos(nwin), nucos(nwin), kk(maxlin,2), aa(maxtan,2), v2ch, v2oh
    external v2ch, v2oh
    integer          :: i, j, k, l, ll, jj, j1, j2, j3, ks, i1, i2
    logical, save    :: qinit(2) = [.false., .false.] ! zeigt an ob Parameter initialisiert.
    double precision :: pq, a, b, c, nenn, rhcsin, rhccos, konus
    double precision :: phc, pco, pcos, poh, rhci, rcoi, rcosi, rohi
    double precision :: thy, aco, ahc, vplus, pi
    complex(real64)  :: qq1, qq2, qq3, qq4, wurz
    if( istate .le. 0 )then
        if( istate .eq. 0 )then
            i1 = 1                         ! 1. Flaeche in Datei (Unit 10).
            i2 = 2                         ! letzte Flaeche     "           .
        else if( istate .eq. -1 )then
            i1 = 1                         ! 1. Flaeche in Datei (Unit 10).
            i2 = 1                         ! letzte Flaeche        "        .
        else if( istate .eq. -2 )then
            i1 = 2                         ! 1. Flaeche in Datei (Unit 10).
            i2 = 2                         ! letzte Flaeche        "        .
        endif
        pi = acos(-1d0)
        do j1=1, nhc
        enddo
        do j2=1, nco
        enddo
        do j3=1, nwin
            nucos(j3) = cos(pi/180d0*nuwin(j3)) ! Nullstellen d. 1d-Polynome.
        enddo
        read(10, *)
        read(10, *)
        read(10, *)   gzhc, gzco, gzwin
        read(10, *)
        read(10, *)   gnhc, gnco, gnwin
        read(10, *)
        read(10, *)   ixhc, ixco, ixoh
        read(10, *)
        read(10, *)   gwhc, gwco, gwwin
        read(10, *)
        read(10, *)   gz2hc, gn2hc, ix2hc
        read(10, *)
        read(10, *)   gz2co, gn2co, ix2co
        read(10, *)
        read(10, *)   gz2oh, gn2oh, ix2oh
        read(10, *)
        read(10, *)   wuradd, nenadd, easym
        do ks=i1, i2
            if( qinit(ks) )then! Test ob schon mal initialisiert.
            else
                qinit(ks) = .true.
            endif
            mlin(ks) = (gzhc+1)*(gzco+1)*(gzwin+1)+(gz2co+1)
            mtan(ks) = 2*(gnhc+1)*(gnco+1)*(gnwin+1)-2+2*gn2co
            if( ks .eq. 1 )then! Mit Konus:
                mlin(ks) = mlin(ks) + (gwhc+1)*(gwco+1)*(gwwin+1)
                mtan(ks) = mtan(ks) + 4
            endif
            read(10, *)
            do i=1, mlin(ks)
                read(10, *) j, kk(i,ks)
                if( i .ne. j ) stop ' POTHCO: INDEX FALSCH '
            enddo
            read(10, *)
            do i=1, mtan(ks)
                read(10, *) j, aa(i,ks)
                if( i .ne. j ) stop ' POTHCO: INDEX FALSCH '
            enddo
        enddo
        pothco = 0d0
    else
        rhci  = rhc                                      ! interne Kopien.
        rcoi  = rco
        rcosi = rcos
        vplus = 0d0                                     ! Potential-Offset.
        thy   = tanh(1.1d0*(rhci-4.4d0))                  ! innen -1, aussen 1.
        a     = 0.75d0 - 0.25d0 * thy                       ! innen 1, aussen 1/2.
        ahc   = 1d0 + (-0.002d0 + 0.014d0 * rcosi ) * a   ! Skalierung in R_HC.
        aco   = 1d0 +         0.012d0           * a       ! Skalierung in R_CO.
        rhci  = ahc * (rhci-2.11d0) + 2.11d0
        rcoi  = aco * (rcoi-2.23d0) + 2.23d0
        vplus = ( 0.5d0 * thy - 0.5d0 ) * 50d-3 / 27.2116d0 ! innen -50 meV.
        rhci  = max(rhci, 1.3d0)            ! R_HC > 1.3 a.u.
        rcoi  = max(rcoi, 1.5d0)            ! R_CO > 1.5 a.u.
        rohi  = sqrt( rhci**2 + rcoi**2 - 2d0*rhci*rcoi*rcosi )
        rohi  = max(rohi, 0.01d0)          ! zur numerischen Stabilitaet.
        ks    = istate
        jj    = 1                           ! Index in AA(JJ,KS) (Nenner).
        qq2   = (1d0, 0d0)                   ! *** QQ_CO(R_CO)
        pco   = 1d0
        do j=1, gn2co
            pco = pco * ( rcoi - nuco(j) ) ! P_{CO,J}(R_CO)
            qq2 = qq2 + dcmplx( aa(jj,ks)*pco, aa(jj+1,ks)*pco )
            jj  = jj + 2
        enddo
        qq4 = (1d0, 0d0)                   ! *** QQ_HCO(R_HC,R_CO,COS())
        phc = 1d0
        do j=0, gnhc
            pco = phc
            do k=0, gnco
                pcos = pco
                if( j .gt. 0.or.k .gt. 0 )then
                    qq4 = qq4 + dcmplx( aa(jj,ks)*pco, aa(jj+1,ks)*pco )
                    jj  = jj + 2
                endif
                do l=1, gnwin
                    pcos = pcos * ( rcosi - nucos(l) )
                    qq4  = qq4 + dcmplx( aa(jj,ks)*pcos, aa(jj+1,ks)*pcos )
                    jj   = jj + 2
                enddo
                pco = pco * ( rcoi - nuco(k+1) )
            enddo
            phc = phc * ( rhci - nuhc(j+1) )
        enddo
        if( ks .eq. 1 )then
            rhccos = rhci * rcosi              ! R_HC * cos(gamma)
            rhcsin = rhci * sqrt(1d0-rcosi**2) ! R_HC * sin(gamma)
            a      = 1d0 + aa(jj+1,ks)*rhccos&                 ! 1 + a2*RHCCOS +
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
        konus = 0d0                      ! Konische Durchschneidung:
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
        pq = pq + 0.055d0/au2ev + vplus         ! OFFSET(Dez.1994)
        if (rohi.lt.1.1d0.or.rhci.lt.1.3d0.or.&
            &rcoi.lt.1.5d0.or.rcoi+rohi.lt.2.9d0) then
        pq = max( 6d0/au2ev, pq )      ! >= 6 eV.
    endif
    pothco = min( 20d0/au2ev, pq )   ! nach oben begrenzen (20eV)
    endif
    999 return
end

!> @brief Evaluate V2OH.
function v2oh( r )
    implicit none
    double precision            :: v2oh, r
    double precision, parameter :: au2ang = 0.529177d0  ! a_0 / ANGSTROEM
    double precision, parameter :: au2ev  = 27.2116d0  ! HARTREE / EV
    double precision            :: a1, a2, a3, de, re, rho
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
    implicit none
    double precision            :: v2ch, r
    double precision, parameter :: au2ang = 0.529177d0  ! a_0 / ANGSTROEM
    double precision, parameter :: au2ev  = 27.2116d0  ! HARTREE / EV
    double precision            :: a1, a2, a3, de, re, rho
    a1   = 3.836
    a2   = 3.511
    a3   = 2.268
    de   = 3.631 / au2ev
    re   = 1.1199
    rho  = ( r*au2ang - re )
    v2ch = -de*(1+a1*rho+a2*rho**2+a3*rho**3 )*exp(-a1*rho)
end

!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_h2co_pip_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, path, p, gradient_out, d )
        use shell,                  only: cof, coff, power, total
        use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        implicit none
        ! number of electronic state
        integer, parameter              :: nstates = 1
        integer, parameter              :: natoms  = 4
        integer, intent(in)             :: igrad
        character(len=1024), intent(in) :: path
        double precision, intent(in)    :: x(natoms,3)
        double precision, intent(out)   :: p(nstates)
        double precision, intent(out)   :: gradient_out(nstates,natoms,3)
        double precision, intent(out)   :: d(nstates,nstates,natoms,3)
        double precision                :: tx(3,natoms), v
        integer                         :: iatom, idir, j, istate
        logical, save                   :: first_time_data=.true.
        !initialize
        gradient_out = 0.d0; d = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)/legacy_bohr_to_angstrom
            enddo
        enddo
        if( igrad == 0 )then
            if( first_time_data )then
                call pes_init(path)
                first_time_data = .false.
            endif
            call getpot(tx, v, path)
            v = v*legacy_hartree_to_electronvolt
            do istate=1, nstates
                p(istate) = v
            enddo
        else
            write (*, *) 'Only energy is available'
        endif
    end subroutine evaluate
end module chempotpy_h2co_pip_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, path, p, gradient_out, d )
    use chempotpy_h2co_pip_pes, only: evaluate
    implicit none
    integer, parameter              :: natoms  = 4
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    double precision, intent(in)    :: x(natoms,3)
    double precision, intent(out)   :: p(nstates)
    double precision, intent(out)   :: gradient_out(nstates,natoms,3)
    double precision, intent(out)   :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, path, p, gradient_out, d)
end subroutine pes
