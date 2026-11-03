!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the H3_GEN_BKMP_1991 potential.
subroutine pes( x, igrad, p, g, d )

    implicit none

    ! number of electronic state
    integer, parameter :: nstates = 1
    integer, parameter :: natoms  = 3

    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)

    double precision :: dedr, energy, r1, r2, r3
    double precision :: r(3), v
    double precision :: tx(9), rx(3), drdx(3,9)
    double precision :: gx(nstates,9)
    integer          :: iatom, idir, i, j
    logical, save    :: first_time_data=.true.

    common /potcm/ r1, r2, r3, energy, dedr(3)

    !initialize
    v = 0.d0
    g = 0.d0
    d = 0.d0

    do iatom=1, natoms
        do idir=1,3
            j     = (iatom-1)*3+idir
            tx(j) = x(iatom,idir)
        enddo
    enddo
    ! input cartesian is HHH
    r(1) = sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2 +(x(1,3)-x(2,3))**2)/0.529177211
    r(2) = sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2 +(x(2,3)-x(3,3))**2)/0.529177211
    r(3) = sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2 +(x(1,3)-x(3,3))**2)/0.529177211

    r1 = r(1)
    r2 = r(2)
    r3 = r(3)
    rx = r*0.529177211

    call prepot
    call pot

    p    = energy*27.211386
    dedr = dedr*51.422067

    call evdrdx(tx, rx, drdx)
    gx = 0.d0
    do i=1,9
        do j=1,3
            gx(:,i) = gx(:,i)+dedr(j)*drdx(j,i)
        enddo
    enddo

    do iatom=1, natoms
        do idir=1,3
            j               = (iatom-1)*3+idir
            g(:,iatom,idir) = gx(:,j)
        enddo
    enddo

    if (igrad==2) then
        write (*,*) 'Only energy and gradient are available'
    endif

end subroutine pes

!> @brief Evaluate EvdRdX for the H3_GEN_BKMP_1991 potential.
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

!
!> @brief Evaluate PREPOT for the H3_GEN_BKMP_1991 potential.
subroutine prepot
    !
    !   System:          H3
    !   Common name:     BKMP
    !   Reference:       A. I. Boothroyd, W. J. Keogh, P. G. Martin,
    !                    and M. R. Peterson
    !                    J. Chem. Phys. 95, 4343 (1991).
    !
    !   Interface:       3-1S
    !
    !   PREPOT must be called once before any calls to POT.
    !   The potential parameters are included in DATA statements.
    !   Coordinates, potential energy, and derivatives are passed
    !   through the common block POTCM:
    !                  COMMON /POTCM/ R(3), ENERGY, DEDR(3)

    !   The potential energy in the three asymptotic valleys are
    !   stored in the common block ASYCOM:
    !                  COMMON /ASYCOM/ EASYAB, EASYBC, EASYAC
    !   The potential energy in the AB valley, EASYAB, is equal to the poten
    !   energy of the H "infinitely" far from the H2 diatomic, with the
    !   H2 diatomic at its equilibrium configuration.  For this potential en
    !   surface the potential energy in the three asymptotic valleys are equ
    !   All the information passed through the common blocks POTCM and ASYCO
    !   is in Hartree atomic units.
    !
    !   The potential energy is defined to be equal to zero when one H is
    !   "infinitely" far from the H2 diatomic, and the H2 diatomic bond
    !   length is equal to the H2 equilibrium bond length.
    !
    !   The the flags that indicate what calculations should be carried out
    !   the potential routine are passed through the common block POTCCM:
    !                  /POTCCM/ NSURF, NDER, IDBUG, NDUM(7)
    !   where:
    !        NSURF - which electronic state should be used.
    !                This option is not used for this potential as only the
    !                ground electronic state is available.
    !        NDER  = 0 => no derivatives should be calculated
    !        NDER  = 1 => calculate first derivatives
    !        IDBUG = 0 => do not print extra information
    !        IDBUG > 0 => print details of the energy calculation
    !        IDBUG > 1 => print details of the derivative calculation
    !                     All extra output is written to FORTRAN unit 7
    !        NDUM  - these 7 integer values can be used to flag options
    !                within the potential; in this potential these options
    !                are not used.
    !
    !   The common block IOCOM contains the FORTRAN unit numbers for the
    !   potential output.  In this potential IOCOM contains one variable, IP
    !                      /IOCOM/ IPRT
    !
    !   Potential parameters' default settings
    !                  Variable            Default value
    !                  NSURF               0
    !                  NDER                1
    !                  IDBUG               0
    !                  IPRT                6
    !
    implicit none

    double precision            :: der, easyab, easyac, easybc, energy, rab, rac, rbc
    integer                     :: idbug, iprt, nder, ndum, nsurf
    double precision, parameter :: verson = 706.0d0
    common /potcm/ rab,rbc,rac,energy,der(3)
    common /potccm/ nsurf, nder, idbug, ndum(7)
    common /iocom/ iprt
    common /asycom/ easyab, easybc, easyac
    !
    easyab = 0.174495771d0
    easybc = easyab
    easyac = easyab
    !
    nsurf = 0
    nder  = 1
    idbug = 0
    iprt  = 6
    !
    !
end subroutine prepot
!
!> @brief Evaluate POT for the H3_GEN_BKMP_1991 potential.
subroutine pot
    implicit none

    double precision :: dedr, easyab, easyac, easybc, energy, r1, r2, r3, rval
    integer          :: idbug, iprt, nder, ndum, nsurf
    common /potcm/ r1, r2, r3, energy, dedr(3)
    common /potccm/ nsurf, nder, idbug, ndum(7)
    common /iocom/ iprt
    common /asycom/ easyab, easybc, easyac
    dimension rval(3)
    rval(1) = r1
    rval(2) = r2
    rval(3) = r3
    !   Check the values of NSURF and NDER for validity.
    if (nsurf .ne. 0) then
        write (iprt, '(/,2X,T5,''NSURF has been set equal to '',I5, /,2X,T5,''This value of NSURF is not allowed for ' // &
            'this '', ''potential, '', /,2X,T5,''only the ground electronic surface, NSURF = 0, ' // &
            ''', ''is available'')') nsurf
        stop 'POT 1'
    endif
    if (nder .gt. 1) then
        write (iprt, '(/, 2X,''POT has been called with NDER = '',I5, /, 2X, ''This value of NDER is not allowed in t' // &
            'his '', ''version of the potential.'')') nder
        stop 'POT 2'
    endif
    !
    call h3tot (rval, energy, dedr, nder, idbug)
    !    Put the zero of energy at the BC asymptotic valley, i.e., at
    !    H infinitely far from the H2 diatomic
    energy = energy + easybc
    !
    !
end subroutine pot
!
!     subroutine h3tot(r,Vtot,dVtot,id)
!> @brief Evaluate h3tot for the H3_GEN_BKMP_1991 potential.
subroutine h3tot( r, vtot, dvtot, id, ipr )
    !-----------------------------------------------------------------
    !     calculate total h3 potential from all of its parts
    !     calculate the derivative too (if id.gt.0)
    ! version of july27/91 ... surf706d.out Cbend values put in
    ! For a discussion of this surface, see:
    !     boothroyd,keogh,martin,peterson j.chem.phys. sept1991 ??
    !     Note: the surface parameter values as published lead to an
    !           anomolously deep van der Waals well for a very compact
    !           H2 molecule (say r=0.8).  After that paper was submitted,
    !           this problem was fixed and the corrected Cbend coefficients
    !           are used in this version of the surface (version 706).
    !           The old coefficients are still in subr.vbcb but have been
    !           commented out.
    ! any questions/problems could be addressed to:
    !     wkeogh@alchemy.chem.utoronto.ca
    ! all distances are in bohrs and all energies are in hartrees
    !
    ! vector 'ipass' is just used to pass print flags around
    !  if ipr.gt.0 print out some details of various contributions to Vtot
    !  if ipr.gt.1  also print out some details about derivatives
    !----------------------------------------------------------------------
    implicit none

    double precision :: cal, cas, cb, cbnda, cbndb, cbp, db1a, db1b, db2, db3, dcal, dcas, dcbnda, &
        & dcbndb, dt, dvas, dvbnda, dvbndb, dvlon, dvtot, expass, r, t, vas, vb, vbnda, vbndb, vbp, &
        & vlon, vtot
    integer :: i, icompc, id, ipass, ipr
    dimension r(3),dvtot(3)
    dimension dvlon(3),dvas(3),dvbnda(3),dvbndb(3), &
        & dcal(3), dcas(3),dcbnda(3),dcbndb(3)
    !     dimension vb(2,25),cb(2,25),vbp(2,25),cbp(2,25)    {debugging}
    dimension vb(2,25),cb(2,25),vbp(2,25),cbp(2,25)
    dimension db1a(3),db1b(3),dt(3),t(3)
    dimension db2(3),db3(6),expass(4)
    common/pass/ipass(100)
    !     ipr = ipass(1)
    !  zero everything to avoid any 'funny' values:
    vlon  = 0.d0
    vas   = 0.d0
    vbnda = 0.d0
    vbndb = 0.d0
    cal   = 0.d0
    cas   = 0.d0
    cbnda = 0.d0
    cbndb = 0.d0
    !     do i=1,3
    do 700 i=1,3
        dvlon(i)  = 0.d0
        dvas(i)   = 0.d0
        dvbnda(i) = 0.d0
        dvbndb(i) = 0.d0
        dcal(i)   = 0.d0
        dcas(i)   = 0.d0
        dcbnda(i) = 0.d0
        dcbndb(i) = 0.d0
        !     end do
    700 continue
    !  zero the vb and cb arrays (used only for numerical derivatives)
    !     do i=1,25
    do 701 i=1,25
        vb(1,i) = 0.d0
        vb(2,i) = 0.d0
        cb(1,i) = 0.d0
        cb(2,i) = 0.d0
        !     end do
    701 continue
    call h3lond( r, vlon, dvlon )
    call vascal( r, vas, dvas )
    !    now do any corrections required for compact geometries:
    call compac(r,icompc,t,dt)
    !    oct.3/90 compact routines only called for compact geometries:
    if(icompc.ge.1)then
        call csym  ( r, cal, dcal )
        call casym ( r, cas, dcas,  1 ,t,dt)
    end if
    call vbcb(r,icompc,t,dt,ipr, &
        & vbnda,vbndb,dvbnda,dvbndb, &
        & cbnda,cbndb,dcbnda,dcbndb)
    !    add up the various parts of the potential:
    vtot = vlon + vas + vbnda + vbndb + cal + cas + cbnda + cbndb
    !    add up the various parts of the derivative:
    dvtot(1) = dvlon(1) + dvas(1) + dvbnda(1) + dvbndb(1) +  dcal(1) + dcas(1) + dcbnda(1) + dcbndb(1)
    dvtot(2) = dvlon(2) + dvas(2) + dvbnda(2) + dvbndb(2) +  dcal(2) + dcas(2) + dcbnda(2) + dcbndb(2)
    dvtot(3) = dvlon(3) + dvas(3) + dvbnda(3) + dvbndb(3) +  dcal(3) + dcas(3) + dcbnda(3) + dcbndb(3)

    if(ipr.gt.0) then
        write(7,*)
        write(7,*) 'h3tot enter --------------------------------'
        write(7,'(5x,''    r ='',3(1x,f16.10),'' icompac = '',i1,/, 5x,''Vtot  = '',f18.12,/, 5x,''Vlon  = '',f18.12,' // &
            '''       Vas   = '',g18.12,/, 5x,''Cal   = '',g18.12,''       Cas   = '',g18.12,/, 5' // &
            'x,''Vbnda = '',g18.12,''       Vbndb = '',g18.12,/, 5x,''Cbnda = '',g18.12,''       ' // &
            'Cbndb = '',g18.12)') r,icompc,vtot, &
            & vlon,vas,cal,cas,vbnda,vbndb,cbnda,cbndb
        if(ipr.gt.1)then
            write(7,'('' dVtot   = '',3(1x,g18.12),/, '' dVlon   = '',3(1x,g18.12),/, '' dVas    = '',3(1x,g18.12),/,' // &
                ' '' dVbnda  = '',3(1x,g18.12),/, '' dVbndb  = '',3(1x,g18.12),/, '' dCal    = '',3(1' // &
                'x,g18.12),/, '' dCas    = '',3(1x,g18.12),/, '' dCbnda  = '',3(1x,g18.12),/, '' dCbn' // &
                'db  = '',3(1x,g18.12))') dvtot,dvlon,dvas,dvbnda,dvbndb, &
                & dcal,dcas,dcbnda,dcbndb
        end if
        write(7,'(''Vba values: '',5(1x,g12.6),/,''Vbb values: '',5(1x,g12.6))') (vb(1,i),i=1,5),(vb(2,i),i=1,5)
        write(7,'(''Cba values: '',3(1x,f16.8),/, ''            '',3(1x,f16.8),/, ''            '',3(1x,f16.8),/, ''C' // &
            'bb values: '',3(1x,f16.8),/, ''            '',3(1x,f16.8),/, ''            '',3(1x,f' // &
            '16.8))') (cb(1,i),i=1,9),(cb(2,i),i=1,9)
        write(7,'(''   Vbnda        Vbndb        Cbnda        Cbndb '',/, 4(e12.6,2x),/, ''dVbnda: '',3(e16.10,1x),/,' // &
            ' ''dVbndb: '',3(e16.10,1x),/, ''dCbnda: '',3(e16.10,1x),/, ''dCbndb: '',3(e16.10,1x)' // &
            ')') vbnda,vbndb,cbnda,cbndb, &
            & dvbnda,dvbndb,dcbnda,dcbndb
        write(7,*) 'exiting subr.h3tot'
        write(7,*) 'h3tot exit ---------------------------------'
    end if
end subroutine h3tot
!
!     subroutine triplet(r,E3,id)
!> @brief Evaluate triplt for the H3_GEN_BKMP_1991 potential.
subroutine triplt( r, e3, id )
    !--------------------------------------------------------------------
    !  version of oct4/90  surface626 values
    !  H2 triplet curve and derivatives:
    !  calculates triplet potential and first derivative
    !  uses truhlar horowitz equation with our extension
    !  uses the Johnson correction at short distances (r < rr)
    !     if r .ge. rr         use modified t/h triplet equation
    !     if r .le. rl         use the Johnson correction
    !     in between           use the transition equation
    !--------------------------------------------------------------------
    implicit none

    double precision :: dr, e1, e3, exdr, r, ra6, ra61, rsq
    integer          :: id, ipass, ipr
    dimension e3(3),e1(3)
    common/pass/ipass(100)
    double precision, parameter :: rl = 0.95d0
    double precision, parameter :: rr = 1.15d0
    double precision, parameter :: z1 = 1.d0
    double precision, parameter :: z2 = 2.d0
    double precision, parameter :: z4 = 4.d0
    !  triplet and johnson values from fit621:
    double precision, parameter :: a1 = -0.0253496194d0
    double precision, parameter :: a2 = -29.2302126444d0
    double precision, parameter :: a3 = -50.7225015503d0
    double precision, parameter :: a4 = 2.0452676876d0
    double precision, parameter :: a5 = -12.2408908509d0
    double precision, parameter :: a6 = 1.6733157383d0
    double precision, parameter :: c1 = -0.4170298146519658d0
    double precision, parameter :: c2 = -0.0746027774843370d0
    double precision, parameter :: c3 = 0.4297899952237434d0
    ! surface parameters from fit601.out
    !     parameter(a1=-0.66129429,a2=-1.99434198,a3=-2.37604328,
    !    .          a4= 2.08107802,a5=-0.0313032510,a6=3.76546699,
    !    .          C1=-0.4222590135447196,C2=-0.0731117796738824,
    !    .          C3= 0.4295918082189010 )
    ipr   = ipass(9)
    e3(3) = 0.d0
    if(r.ge.rr )then
        !       modified truhlar/horowitz triplet equation:
        exdr  = exp( -a4*r )
        rsq   = r*r
        ra6   = r**(-a6)
        e3(1) = a1* ( a2 + r + a3*rsq + a5*ra6 )*exdr
        !       first derivative of triplet curve:
        ra61  = r**(-a6-z1)
        e3(2) = a1*exdr* ( z1 -a2*a4 +(z2*a3-a4)*r -a3*a4*rsq -a5*a6*ra61 -a4*a5*ra6 )
    end if
    if(r.lt.rr )then
        dr = r- rl
        call vh2opt(r,e1,2)
        if( r.le.rl ) then
            !          Johnson triplet equation:
            e3(1) = e1(1) + c2*dr + c3
            e3(2) = e1(2) + c2
        else
            !          Transition equation:
            e3(1) = e1(1) + c1*dr*dr*dr + c2*dr + c3
            e3(2) = e1(2) + 3.d0*c1*dr*dr + c2
        end if
    end if
    if(ipr.gt.0) write(7,'(''  triplet:  E3 = '',3(1x,f12.8))') e3
end subroutine triplt
!
!
!> @brief Evaluate vH2opt for the H3_GEN_BKMP_1991 potential.
subroutine vh2opt( r, e, ideriv )
    !-----------------------------------------------------------------
    !  version of jul 09/90 ... super duper speedy version
    ! self-contained version of schwenke's H2 potential
    ! all distances in bohrs and all energies in Hartrees
    ! (1st deriv added on May 2 1989; 2nd deriv on May 28 1989)
    !-----------------------------------------------------------------
    implicit none

    double precision :: alphar, dalphr, ddalph, ddvlr, ddvsr, dvlr, dvsr, e, exalph, r, r02, r04, r06, r10, r11, r12, &
        & r13, r14, r15, r2, r3, r4, r5, r6, r7, r8, r9, rr2, rr25, rr26, rr27, rr4, rr43, rr44, rr6, rr62, rr63, vlr, &
        & vsr
    integer :: ideriv, ipass
    dimension e(3)
    common/pass/ipass(100)
    double precision, parameter :: a0  = 0.03537359271649620d0
    double precision, parameter :: a1  = 2.013977588700072d0
    double precision, parameter :: a2  = -2.827452449964767d0
    double precision, parameter :: a3  = 2.713257715593500d0
    double precision, parameter :: a4  = -2.792039234205731d0
    double precision, parameter :: a5  = 2.166542078766724d0
    double precision, parameter :: a6  = -1.272679684173909d0
    double precision, parameter :: a7  = 0.5630423099212294d0
    double precision, parameter :: a8  = -0.1879397372273814d0
    double precision, parameter :: a9  = 0.04719891893374140d0
    double precision, parameter :: a10 = -0.008851622656489644d0
    double precision, parameter :: a11 = 0.001224998776243630d0
    double precision, parameter :: a12 = -1.227820520228028d-04
    double precision, parameter :: a13 = 8.638783190083473d-06
    double precision, parameter :: a14 = -4.036967926499151d-07
    double precision, parameter :: a15 = 1.123286608335365d-08
    double precision, parameter :: a16 = -1.406619156782167d-10
    double precision, parameter :: r0  = 3.5284882d0
    double precision, parameter :: dd  = 0.160979391d0
    double precision, parameter :: c6  = 6.499027d0
    double precision, parameter :: c8  = 124.3991d0
    double precision, parameter :: c10 = 3285.828d0
    !     Ediss = 0.174445d0
    e(1) = -999.d0
    e(2) = -999.d0
    e(3) = -999.d0
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
    !     general term:  a(i)*r(i-1),  i=0,16
    alphar =  a0/r + a1 + a2 *r   + a3 *r2  + a4 *r3  + a5 *r4  + a6 *r5 + a7 *r6  + a8 *r7  + a9 *r8  + a10*r9  + &
        & a11* r10 + a12*r11 + a13*r12 + a14*r13 + a15*r14 + a16*r15
    exalph = exp(alphar)
    vsr    = dd*(exalph-1.d0)*(exalph-1.d0) - dd
    vlr    = -c6/rr6   -c8/(rr4*rr4)   -c10/rr25
    e(1)   = vsr +  vlr
    !
    !  calculate first derivative if required:
    if(ideriv.ge.1)then
        r3   = r2*r
        r5   = r4*r
        rr26 = rr25*rr2
        rr43 = rr4*rr4*rr4
        !        general term:  (i-1)*a(i)*r**(i-2)  ,  i=0,16
        dalphr = -a0/r2 + a2 + 2.d0*a3*r +  3.d0 *a4 *r2  +  4.d0 *a5 *r3  +  5.d0 *a6 *r4 +  6.d0 *a7 *r5  +  7.d0 * &
            & a8 *r6  +  8.d0 *a9 *r7 +  9.d0 *a10*r8  + 10.d0 *a11*r9  + 11.d0 *a12*r10 + 12.d0 *a13*r11 + 13.d0 * &
            & a14* r12 + 14.d0 *a15*r13 + 15.d0 *a16*r14
        dvsr = 2.d0 *dd *(exalph-1.d0) *exalph *dalphr
        dvlr = 6.d0*c6*r5 / (rr6*rr6) +  8.d0*c8*r3 / rr43  + 10.d0*c10*r / rr26
        e(2) = dvsr + dvlr
    end if
    !
    !  calculate second derivative if required:
    if(ideriv.ge.2)then
        r10  = r6*r4
        rr27 = rr26*rr2
        rr44 = rr43*rr4
        rr62 = rr6*rr6
        rr63 = rr62*rr6
        !        general term: (i-1)*(i-2)*a(i)*r**(i-3),    i=0,16
        ddalph = 2.d0*a0/r3 +  2.d0*a3      +  6.d0*a4 *r + 12.d0*a5 *r2    + 20.d0*a6 *r3  + 30.d0*a7 *r4 + 42.d0* &
            & a8 * r5    + 56.d0*a9 *r6  + 72.d0*a10*r7 + 90.d0*a11*r8    +110.d0*a12*r9  +132.d0*a13*r10 +156.d0*a14* &
            & r11   + 182.d0*a15*r12 +210.d0*a16*r13
        ddvsr = 2.d0 *dd *exalph *( (2.d0*exalph-1.d0)*dalphr*dalphr + (exalph-1.d0)*ddalph )
        ddvlr =- 72.d0* c6 *r10 / rr63 -96.d0 *c8 *r6 / rr44 -120.d0*c10 *r2  / rr27 +30.d0 *c6 *r4 / rr62 + 24.d0* &
            & c8 * r2  / rr43 +10.d0 *c10    / rr26
        e(3) = ddvsr + ddvlr
    end if
end subroutine vh2opt
!
!> @brief Evaluate h3lond for the H3_GEN_BKMP_1991 potential.
subroutine h3lond( r, vlon, dvlon )
    !----------------------------------------------------------------------
    ! version of may 12/90  ... derivatives corrected (0.5 changed to 0.25)
    ! calculates the h3 london terms and derivatives
    ! modified oct 7/89 to include eps**2 term which rounds off the
    ! cusp in the h3 potential which occurs at equilateral triangle
    ! configurations
    !----------------------------------------------------------------------
    implicit none

    double precision :: de1, de3, dj1, dj2, dj3, dvlon, e1, e3, esing, etrip, r, rootjt, sumj, sumq, vlon
    integer          :: i, ipass, ipr
    !     real*8 Q(3),J(3),Jt
    double precision :: q(3),j(3),jt
    dimension r(3),e1(3),e3(3),esing(3),etrip(3),dvlon(3)
    dimension de1(3),de3(3),dj1(3),dj2(3),dj3(3)
    common/pass/ipass(100)
    double precision, parameter :: half = 0.5d0
    double precision, parameter :: two  = 2.d0
    double precision, parameter :: eps2 = 1.d-12
    ipr = ipass(2)
    !     do i=1,3
    do 702 i=1,3
        call vh2opt(r(i),esing,2)
        e1(i)  = esing(1)
        de1(i) = esing(2)
        !        call triplet(r(i),etrip,2)
        call triplt(r(i),etrip,2)
        if(ipr.gt.0) write(7,'(''from subr.london: '',/, ''  using r'',i1,'':   Esinglet='',3(1x,f12.8),/, ''        ' // &
            ' '',1x,''    Etriplet='',3(1x,f12.8))') i,esing,etrip
        e3(i)  = etrip(1)
        de3(i) = etrip(2)
        q(i)   = half*(e1(i) + e3(i))
        j(i)   = half*(e1(i) - e3(i))
        !     end do
    702 continue
    sumq   = q(1) + q(2) + q(3)
    sumj   = abs( j(2)-j(1) )**2 + abs( j(3)-j(2) )**2 + abs( j(3)-j(1) )**2
    jt     = half*sumj + eps2
    rootjt = sqrt(jt)
    vlon   = sumq - rootjt
    if(ipr.gt.0) then
        write(7,'(''         sumQ = '',g12.6,''        sumJ = '',g12.6)') sumq,sumj
        write(7,'(''         Vlon = '',f12.8,''      rootJt = '',g12.6)') vlon,rootjt
    end if
    !  calculate the derivatives with respect to r(i):
    dvlon(1) = half*(de1(1)+de3(1)) - 0.25d0*(two*j(1)-j(2)-j(3))*(de1(1)-de3(1))/rootjt
    dvlon(2) = half*(de1(2)+de3(2)) - 0.25d0*(two*j(2)-j(3)-j(1))*(de1(2)-de3(2))/rootjt
    dvlon(3) = half*(de1(3)+de3(3)) - 0.25d0*(two*j(3)-j(1)-j(2))*(de1(3)-de3(3))/rootjt
    if(ipr.gt.0) then
        write(7,'(''             r = '',3(1x,f12.6),/, ''            E1 = '',3(1x,f12.8),/, ''            E3 = '',3(1' // &
            'x,f12.8),/, ''          vlon = '',1x,f12.8)') r,e1,e3,vlon
        write(7,'(13x,''q = '',3(1x,e12.6),/,13x,''j = '',3(1x,e12.6))') q,j
        write(7,'(''         dVlon = '',3(1x,g12.6))') dvlon
    end if
end subroutine h3lond
!
!> @brief Evaluate Vascal for the H3_GEN_BKMP_1991 potential.
subroutine vascal( rpass, vas, dvas )
    !------------------------------------------------------------------
    !  version of oct11/90  fit632c.out values
    !  version of oct5/90   surf626 values
    !  calculate the asymmetric correction term and its derivatives
    !  see equations [14] to [16] of truhlar/horowitz 1978 paper
    !------------------------------------------------------------------
    implicit none

    double precision :: a, a2, a3, a4, a5, da, ds, dvas, exp1, exp6, r, r1, r2, r3, rcu, rpass, rsq, s, vas
    integer          :: ipass, ipr
    common/pass/ipass(100)
    dimension dvas(3),da(3),ds(3),rpass(3)
    !--- Vasym values from fit632c ----------------------
    double precision, parameter :: aa1 = 0.3438222224d-02
    double precision, parameter :: aa2 = 0.1398145763d-02
    double precision, parameter :: aa3 = -.1923999449d-03
    double precision, parameter :: aa4 = 0.9712737075d-05
    double precision, parameter :: aa5 = -.1263794562d-06
    double precision, parameter :: aa6 = 0.5181432712d+00
    double precision, parameter :: aa7 = -.9487002995d-03
    ipr = ipass(3)
    r1  = rpass(1)
    r2  = rpass(2)
    r3  = rpass(3)
    r   = r1 + r2 + r3
    rsq = r*r
    rcu = rsq*r
    !   calculate the Vas term first (eq.14 of truhlar/horowitz)
    call acalc(r1,r2,r3,a,da)
    a2      = a *a
    a3      = a2*a
    a4      = a3*a
    a5      = a4*a
    exp1    = exp(-aa1*rcu)
    exp6    = exp(-aa6*r)
    s       = aa2*a2 + aa3*a3 + aa4*a4 + aa5*a5
    vas     = s*exp1 +  aa7*a2 *exp6 / r
    ds(1)   = ( 2.d0*aa2*a  +3.d0*aa3*a2 +4.d0*aa4*a3 +5.d0*aa5*a4) * da(1)
    dvas(1) = -3.d0*aa1*rsq*s*exp1 + ds(1)*exp1 -aa7*a2*exp6/rsq  + 2.d0*aa7*a*da(1)*exp6/r -aa6*aa7*a2*exp6/r
    ds(2)   = ( 2.d0*aa2*a  +3.d0*aa3*a2 +4.d0*aa4*a3 +5.d0*aa5*a4) * da(2)
    dvas(2) = -3.d0*aa1*rsq*s*exp1 + ds(2)*exp1 -aa7*a2*exp6/rsq  + 2.d0*aa7*a*da(2)*exp6/r -aa6*aa7*a2*exp6/r
    ds(3)   = ( 2.d0*aa2*a  +3.d0*aa3*a2 +4.d0*aa4*a3 +5.d0*aa5*a4) * da(3)
    dvas(3) = -3.d0*aa1*rsq*s*exp1 + ds(3)*exp1 -aa7*a2*exp6/rsq  + 2.d0*aa7*a*da(3)*exp6/r -aa6*aa7*a2*exp6/r
    if(ipr.gt.0) write(7,'(''    from subr.vascal:   Vas  = '',  1x,g12.6, /, ''                        dS   = '',3(1' // &
        'x,g12.6),/, ''                        dVas = '',3(1x,g12.6))') vas,ds,dvas
end subroutine vascal
!
!> @brief Evaluate acalc for the H3_GEN_BKMP_1991 potential.
subroutine acalc( r1, r2, r3, a, da )
    !---------------------------------------------------------------------
    !  version of may 1, 1990
    !  based on equations from notes date april 4, 1990
    !  calculate the A value and its derivatives
    !           A  = dabs[ (r1-r2)*(r2-r3)*(r3-r1) ]
    !---------------------------------------------------------------------
    implicit none

    double precision :: a, da, r1, r2, r3
    integer          :: ipass, ipr
    dimension da(3)
    common/pass/ipass(100)
    ipr   = ipass(11)
    a     = (r1-r2)*(r2-r3)*(r3-r1)
    da(1) = ( -2.d0*r1 + r2 + r3 )*(r2-r3)
    da(2) = ( -2.d0*r2 + r3 + r1 )*(r3-r1)
    da(3) = ( -2.d0*r3 + r1 + r2 )*(r1-r2)
    if(a.lt.0.d0)then
        a     = -a
        da(1) = -da(1)
        da(2) = -da(2)
        da(3) = -da(3)
    end if
    if(ipr.gt.0) write(7,'(''Acalc enter ---------------------------------------'',/, ''     A = '',f15.10,''   dA = ' // &
        ''',3(1x,g12.6),/, ''Acalc exit ----------------------------------------'')') a,da
end subroutine acalc
!
!> @brief Evaluate compac for the H3_GEN_BKMP_1991 potential.
subroutine compac( r, icompc, t, dt )
    !---------------------------------------------------------------
    !  version of may 14/90 ... calculates T and dT values also
    !  decide whether or not this particular geometry is compact,
    !  that is, are any of the three distances smaller than the
    !  rr value from the johnson correction.
    !-----------------------------------------------------------------
    implicit none

    double precision :: bot, bot2, dt, r, t, top, top2, top3
    integer          :: i, icompc, ipass
    dimension r(3),t(3),dt(3)
    common/pass/ipass(100)
    double precision, parameter :: rr = 1.15d0
    double precision, parameter :: rp = 1.25d0
    !  first see if this is a compact geometry:
    !     do i=1,3
    do 703 i=1,3
        t(i)  = 0.d0
        dt(i) = 0.d0
        !     end do
    703 continue
    icompc = 0
    if(r(1).lt.rr) icompc = icompc + 1
    if(r(2).lt.rr) icompc = icompc + 1
    if(r(3).lt.rr) icompc = icompc + 1
    if(icompc.eq.0) return
    ! calculate the T(i) values:
    !     do i=1,3
    do 704 i=1,3
        if(r(i).lt.rr) then
            top   = rr-r(i)
            bot   = rp-r(i)
            top2  = top  * top
            top3  = top2 * top
            bot2  = bot  * bot
            t(i)  = top3/bot
            dt(i) = -3.d0*top2/bot + top3/bot2
        end if
        !     end do
    704 continue
end subroutine compac
!
!> @brief Evaluate casym for the H3_GEN_BKMP_1991 potential.
subroutine casym ( r, cas, dcas, id, t, dt )
    !---------------------------------------------------------------
    !  version of sept14/90
    !  the compact asymmetric correction term and derivatives
    !---------------------------------------------------------------
    implicit none

    double precision :: a, a2, cas, da, dcas, determ, dpr, dseris, dsumt, dt, dterm1, eterm, pr, pr2, pr3, r, series, sr, &
        & sr2, sr3, sumt, t, term1
    integer :: i, id, ipass, ipr
    dimension r(3),t(3)
    !     dimension dPr(3),dT(3),dsumT(3),dterm1(3),determ(3),dseries(3),
    dimension dpr(3),dt(3),dsumt(3),dterm1(3),determ(3),dseris(3), &
        & dcas(3),da(3)
    common/pass/ipass(100)
    !---- Casym values from fit633a
    double precision, parameter :: u1  = 0.1537481166d+00
    double precision, parameter :: u2  = 0.2745950036d+00
    double precision, parameter :: u3  = 0.7501206780d-02
    double precision, parameter :: u4  = 0.3119136023d+01
    double precision, parameter :: u5  = 0.9969170798d+00
    double precision, parameter :: u6  = -.3373682823d+01
    double precision, parameter :: u7  = 0.6807215913d+00
    double precision, parameter :: u8  = -.4920325491d-01
    double precision, parameter :: u9  = -.3919467989d+01
    double precision, parameter :: u10 = 0.5085532326d+01
    double precision, parameter :: u11 = -.1415264778d+01
    double precision, parameter :: u12 = 0.1138681785d+00
    double precision, parameter :: u13 = -.1525367566d-02
    ipr     = ipass(6)
    cas     = 0.d0
    dcas(1) = 0.d0
    dcas(2) = 0.d0
    dcas(3) = 0.d0
    call acalc(r(1),r(2),r(3),a,da)
    a2 = a*a
    !     if(A.eq.0.d0) return
    sumt = t(1) + t(2) + t(3)
    sr   = r(1) + r(2) + r(3)
    pr   = r(1) * r(2) * r(3)
    sr2  = sr*sr
    sr3  = sr2*sr
    pr2  = pr*pr
    pr3  = pr2*pr
    !    write out the series explicitly:
    series = 1.d0 + u4/pr2 +  u5/pr + u6 +   u7*pr +  u8*pr2 + a*(u9/pr2 + u10/pr + u11 + u12*pr + u13*pr2)
    term1  = u1/pr**u2
    eterm  = exp(-u3*sr3)
    cas    = sumt * a2 * term1 * series * eterm
    if(ipr.gt.0) then
        write(7,*)
        write(7,*) 'Casym: ------------------'
        write(7,*) ' T(i) = ',t
        write(7,*) 'dT(i) = ',dt
        write(7,*) ' id = ',id
        write(7,'(''subr.Cas:   series = '',g12.6,''     term1 = '',g12.6,/, ''             eterm = '',g12.6,''      ' // &
            'Cas  = '',g12.6)') series,term1,eterm,cas
        write(7,*) ' Cas = ',cas
    end if
    !     if(id.gt.0)then {calculate the derivative as well}
    if(id.gt.0)then
        dpr(1) = r(2)*r(3)
        dpr(2) = r(3)*r(1)
        dpr(3) = r(1)*r(2)
        !        do i=1,3
        do 705 i=1,3
            dsumt(i)  = dt(i)
            dterm1(i) = -1.d0*u1*u2*pr**(-u2-1.d0)*dpr(i)
            determ(i) = eterm *(-3.d0*u3*sr2)
            !         dseries(i)=
            dseris(i)= dpr(i)*( -2.d0*u4/pr3  -u5/pr2 +u7  +2.d0*u8*pr    ) +da(i)*( u9/pr2 +u10/pr +u11 +u12*pr +u13* &
                & pr2    ) +a*dpr(i)*( -2.d0*u9/pr3 -u10/pr2 +u12 + 2.d0*u13*pr )
            dcas(i) = &
                & dsumt(i)      * a2    * term1 * series * eterm&
                &+ 2.d0*a*da(i)  * sumt  * term1 * series * eterm&
                &+ dterm1(i)     * sumt  * a2    * series * eterm&
                !    .      + dseries(i)    * sumt  * a2    * term1  * eterm
                &+ dseris(i)     * sumt  * a2    * term1  * eterm&
                &+ determ(i)     * sumt  * a2    * term1  * series
            !        end do
        705 continue
        if(ipr.gt.0) then
            !           write(7,7500) dPr,dt,dterm1,determ,dseries,dCas
            write(7,'(''    derivatives: dPr = '',3(1x,g12.6),/, ''                 dT  = '',3(1x,g12.6),/, ''       ' // &
                '       dterm1 = '',3(1x,g12.6),/, ''              determ = '',3(1x,g12.6),/,  ''    ' // &
                '          dseris = '',3(1x,g12.6),/, ''              dCas   = '',3(1x,g12.6))') dpr,dt,dterm1,determ,dseris,dcas
            write(7,*) ' dCas = ',dcas
        end if
    end if
end subroutine casym
!
!> @brief Evaluate csym for the H3_GEN_BKMP_1991 potential.
subroutine csym( r, cal, dcal )
    !---------------------------------------------------------------
    !  version of oct12/90  surf636 values
    !  calculate the 'compact all' correction term and derivatives
    !  a correction term (added sept 11/89), which adds a small
    !  correction to all compact geometries
    !---------------------------------------------------------------
    implicit none

    double precision :: cal, dcal, dexp3, dg, exp3, g, r, ri, rpri, rpri2, rrri, rrri2, rrri3, sr, sr2, sr3, sumg, sumv, &
        & t
    integer :: i, ipass, ipr
    dimension r(3),t(3),dcal(3),g(3)
    dimension dg(3),sumv(3)
    common/pass/ipass(100)
    double precision, parameter :: rr = 1.15d0
    double precision, parameter :: rp = 1.25d0
    !---- Csym values from fit633a
    double precision, parameter :: v1 = -.2070776049d+00
    double precision, parameter :: v2 = -.5350898737d+00
    double precision, parameter :: v3 = 0.1011861942d-01
    cal   = 0.d0
    ipr   = ipass(5)
    sr    = r(1)+r(2)+r(3)
    sr2   = sr*sr
    sr3   = sr*sr2
    exp3  = exp( -v3*sr3 )
    dexp3 = -3.d0*v3*sr2*exp3
    !     do i=1,3
    do 706 i=1,3
        ri      = r(i)
        rrri    = rr-ri
        rrri2   = rrri*rrri
        rrri3   = rrri*rrri2
        rpri    = rp-ri
        rpri2   = rpri*rpri
        g(i)    = 0.d0
        dg(i)   = 0.d0
        sumv(i) = v1+v1*v2*ri
        if(ri.lt.rr) then
            g(i)  = (rrri3/rpri) *sumv(i)
            dg(i) = (rrri3/rpri2)*sumv(i) -3.d0*(rrri2/rpri)*sumv(i) + (rrri3/rpri)*v1*v2
        end if
        !     end do
    706 continue
    sumg    = g(1) + g(2) + g(3)
    cal     = sumg*exp3
    dcal(1) = dg(1)*exp3 + sumg*dexp3
    dcal(2) = dg(2)*exp3 + sumg*dexp3
    dcal(3) = dg(3)*exp3 + sumg*dexp3
    if(ipr.gt.0)then
        write(7,*)
        write(7,'(''Csym enter ----------------------------------------'',/, ''    Cal = '',1x,g20.14)') cal
        write(7,'(''      G = '',3(1x,g12.6),/, ''     dG = '',3(1x,g12.6),/, ''   sumv = '',3(1x,g12.6),/, ''   dCal' // &
            ' = '',3(1x,g16.10),/, ''Csym exit -----------------------------------------'')') g,dg,sumv,dcal
        write(7,*) 'exiting subr.csym'
    end if
end subroutine csym
!
!> @brief Evaluate VbCb for the H3_GEN_BKMP_1991 potential.
subroutine vbcb( rpass, icompc, t, dt, ipr, &
        & vbnda, vbndb, dvbnda, dvbndb, &
        & cbnda, cbndb, dcbnda, dcbndb )
    !---------------------------------------------------------------
    !  jul24/91 fit705 cbend values added
    !  feb27/91 modified to match equation in h3 paper more closely
    !  nov 4/90 Vbend coefficients now a,g  Cbend coeff's still c,d
    !  in this version, the derivatives are always calculated
    !  subroutines vbend, cbend and B1ab all combined into this one
    !  module in order to improve efficiency by not having to pass
    !  around B1a,B1b,B2,B3a,B3b functions and derivatives
    !   B1a = 1 - sum of [ P1(cos(theta(i))) ]
    !   B1b = 1 - sum of [ P3(cos(theta(i))) ]
    !--------------------------------------------------------------
    implicit none

    double precision :: aasum, asum, b1, b12, b13, b14, b15, b1a, b1b, b2, b3, b3b, bsum, c1, c1cube, c2, c2cube, c3, &
        & c3cube, cb, cbnda, cbndb, cbnds, cos3t1, cos3t2, cos3t3, csum, cx1, d1, d2, d3, daasum, dasum, db1a, db1b, db2, &
        & db3, dbsum, dc1dr1, dc1dr2, dc1dr3, dc2dr1, dc2dr2, dc2dr3, dc3dr1, dc3dr2, dc3dr3, dcb1, dcb2, dcb3, dcb4, &
        & dcb5, dcb6, dcb7, dcb8, dcbnda, dcbndb, dcbnds, ddsum1, ddsum2, dexp1, dexp2, dexp7, dffsum, dfsum, dgsum, dp, &
        & dsum1, dsum2, dt, dvb1, dvb2, dvb3, dvb4, dvb5, dvbnd, dvbnda, dvbndb, dx1, eps2, exp1, exp2, exp7, ffsum, &
        & fsum, gsum, p, pcu, psq, r, r1, r2, r3, rcu, rpass, rsq, sum, sumb, sumt, t, t1, t2, t3, th, vb, vbnd, vbnda, &
        & vbndb
    integer :: icompc, ipr
    !    arrays required for B1,B2,B3 calculations:
    dimension rpass(3),db1a(3),db1b(3),th(3), &
        & db2(3),db3(6)
    !    arrays required for Vbend calculations:
    dimension dvbnda(3),dvbndb(3)
    dimension dvb1(3),dvb2(3),dvb3(3),dvb4(3),dvb5(3),dvbnd(3)
    dimension vbnd(2)
    dimension vb(2,25)
    !     dimension vb(2,25)   {passed out for debugging only}
    !    arrays required for Cbend calculations:
    dimension cbnds(2),dcbnds(2,3)
    dimension dcb1(3),dcb2(3),dcb3(3),dcb4(3),dcb5(3),dcb6(3), &
        & dcb7(3),dcb8(3)
    dimension t(3),dp(3)
    dimension dt(3),dcbnda(3),dcbndb(3)
    !     dimension cb(2,25)   {passed out for debugging only}
    dimension cb(2,25)
    !     parameter( z58=0.625d0, z38=0.375d0 )  {3/8 and 5/8}
    double precision, parameter :: z58 = 0.625d0
    double precision, parameter :: z38 = 0.375d0
    !   parameters required for Vbend calculations:
    double precision, parameter :: beta1 = 0.52d0
    double precision, parameter :: beta2 = 0.052d0
    double precision, parameter :: beta3 = 0.79d0
    !-- Vbenda values from fit634b
    double precision, parameter :: a11 = -.2918252280d+03
    double precision, parameter :: a12 = 0.2164569141d+02
    double precision, parameter :: a13 = -.4005497543d+00
    double precision, parameter :: a21 = -.2890774947d+01
    double precision, parameter :: a22 = 0.1032032542d+02
    double precision, parameter :: a23 = 0.2681748056d+02
    double precision, parameter :: a24 = 0.2633751055d+02
    double precision, parameter :: a31 = 0.6180641351d+01
    double precision, parameter :: a32 = 0.5037667041d+01
    double precision, parameter :: a41 = -.1125570079d+00
    double precision, parameter :: a42 = -.3176529304d-02
    double precision, parameter :: a43 = 0.9068915355d+00
    double precision, parameter :: a44 = -.7228418516d+00
    double precision, parameter :: a51 = 0.2785898232d+03
    double precision, parameter :: a52 = 0.4764442446d+02
    double precision, parameter :: a53 = 0.8621545662d+01
    !-- Vbendb values from fit634b
    double precision, parameter :: g11 = -.4241912309d+02
    double precision, parameter :: g12 = 0.2951365281d+01
    double precision, parameter :: g13 = -.4840201514d-01
    double precision, parameter :: g21 = -.1159168549d-03
    double precision, parameter :: g22 = 0.8688003567d-02
    double precision, parameter :: g23 = -.3486923900d-02
    double precision, parameter :: g24 = 0.8312212839d-03
    double precision, parameter :: g31 = 0.5621810473d-01
    double precision, parameter :: g32 = -.9776930747d-02
    double precision, parameter :: g41 = -.1178456251d-01
    double precision, parameter :: g42 = 0.3491086729d-02
    double precision, parameter :: g43 = 0.7430516993d-01
    double precision, parameter :: g44 = -.9643636957d-01
    double precision, parameter :: g51 = 0.4735533782d+02
    double precision, parameter :: g52 = 0.3001038808d+01
    double precision, parameter :: g53 = 0.1896630453d+01
    ! parameters required for Cbend calculations:
    !-- Cbenda values from fit702c -------------------
    !x    parameter(
    !x   .c11=0.3408013485D+04,c12=-.2569675505D+04,c13=-.4171156109D+03,
    !x   .c14=-.6179297444D+04,c15=0.3235223485D+04,c21=-.3409049347D+03,
    !x   .c22=0.9612365597D+02,c23=0.4104361425D+03,c24=0.6238298637D+03,
    !x   .c31=0.2655770314D+04,c32=0.5314085916D+03,c41=0.2222002477D+02,
    !x   .c42=-.6611541362D+01,c43=-.2541255632D+03,c44=0.7752854954D+02,
    !x   .c51=0.9613217362D+04,c52=-.7010043704D+03,c53=0.1406657509D+03,
    !x   .c61=-.1677846467D+04,c62=0.7411196008D+04,c63=-.1104074171D+05,
    !x   .c71=-.1077399844D+03,c72=-.6598646429D+01,c73=0.2035775472D+00,
    !x   .c74=0.1097670442D+04,c75=0.4873456879D+02,c81=0.1005135589D+03,
    !x   .c82=-.1972894789D+02,c83=-.8495015577D+04,c84=-.6106570657D+02)
    !-- Cbendb values from fit702c -------------------
    !x    parameter(
    !x   .d11=0.8450990097D+03,d12=0.1466587231D+03,d13=0.1239413655D+02,
    !x   .d14=-.5649327722D+03,d15=0.3112693920D+03,d21=0.8029625257D-01,
    !x   .d22=0.2379962852D-01,d23=0.1233566879D+00,d24=-.3582414079D-01,
    !x   .d31=-.1012436386D+03,d32=0.8904152887D-02,d41=0.8740633605D+01,
    !x   .d42=-.2937137124D+01,d43=-.5497174488D+02,d44=0.1656850593D+02,
    !x   .d51=0.3859496208D+05,d52=0.8458336162D+02,d53=-.1431690434D+02,
    !x   .d61=0.5487645288D+02,d62=-.1954058540D+03,d63=-.7226235345D+02,
    !x   .d71=-.1895543609D+02,d72=0.1281231795D+00,d73=0.1291641223D+00,
    !x   .d74=-.3052282568D+03,d75=-.6194562183D+01,d81=-.1746535079D+02,
    !x   .d82=-.8402855300D+00,d83=-.3878264141D+05,d84=0.5655392575D+01)
    !-- Cbenda values from fit705a.out -------------
    !xx   parameter(
    !xx  .c11=0.7138332997D+04,c12=-.3941467854D+04,c13=0.1214839191D+03,
    !xx  .c14=-.9715171923D+04,c15=0.4549112075D+04,c21=-.7153573233D+02,
    !xx  .c22=0.1410801562D+03,c23=0.5310286989D+03,c24=0.6632669305D+03,
    !xx  .c31=0.9648491179D+03,c32=0.1336994724D+03,c41=0.3986106227D+01,
    !xx  .c42=-.2071394753D+01,c43=-.1304296435D+03,c44=0.4433517239D+02,
    !xx  .c51=-.2160974730D+05,c52=-.5040404494D+03,c53=0.1150251071D+03,
    !xx  .c61=-.1900650362D+04,c62=0.8065979973D+04,c63=-.1022333250D+05,
    !xx  .c71=-.4023633922D+02,c72=-.5593869000D+01,c73=0.2901552041D+00,
    !xx  .c74=0.1332299065D+04,c75=0.4830933142D+02,c81=0.1527546398D+03,
    !xx  .c82=-.2277696550D+02,c83=0.2214880233D+05,c84=-.5951714583D+02)
    !-- Cbendb values from fit705a.out -------------
    !xx   parameter(
    !xx  .d11=0.1733880674D+04,d12=-.5105474641D+03,d13=0.5739180116D+02,
    !xx  .d14=-.1549259263D+04,d15=0.6494865366D+03,d21=-.1180847496D+00,
    !xx  .d22=0.3847529877D-01,d23=0.1035799123D+00,d24=-.2783059123D-01,
    !xx  .d31=-.4049790116D+02,d32=0.3623102614D-01,d41=0.2132254354D+01,
    !xx  .d42=-.1337841873D+01,d43=-.2911759497D+02,d44=0.9698950726D+01,
    !xx  .d51=0.7851192036D+05,d52=0.1244750132D+03,d53=-.2157749426D+02,
    !xx  .d61=0.1945719436D+02,d62=-.4336295493D+02,d63=0.8025442515D+02,
    !xx  .d71=0.5232787178D+01,d72=-.2082854432D+01,d73=0.1740967673D+00,
    !xx  .d74=-.3069559809D+02,d75=0.9099007208D+01,d81=-.6088433744D+01,
    !xx  .d82=-.9840250626D+00,d83=-.7880550115D+05,d84=-.9084945347D+01)
    !-- Cbenda values from fit706d.out
    double precision, parameter :: c11 = 0.7107064647d+04
    double precision, parameter :: c12 = -.3728421267d+04
    double precision, parameter :: c13 = 0.1757654624d+03
    double precision, parameter :: c14 = -.9725998132d+04
    double precision, parameter :: c15 = 0.4665086074d+04
    double precision, parameter :: c21 = -.4435165986d+02
    double precision, parameter :: c22 = 0.1604477309d+03
    double precision, parameter :: c23 = 0.5805142925d+03
    double precision, parameter :: c24 = 0.6892349445d+03
    double precision, parameter :: c31 = 0.6581730442d+03
    double precision, parameter :: c32 = 0.6078389739d+02
    double precision, parameter :: c41 = 0.2885182566d+01
    double precision, parameter :: c42 = -.1728169916d+01
    double precision, parameter :: c43 = -.1119535503d+03
    double precision, parameter :: c44 = 0.4052536250d+02
    double precision, parameter :: c51 = 0.2540505673d+03
    double precision, parameter :: c52 = -.5762083627d+03
    double precision, parameter :: c53 = 0.1295901320d+03
    double precision, parameter :: c61 = -.2131706075d+04
    double precision, parameter :: c62 = 0.9084452020d+04
    double precision, parameter :: c63 = -.1138253963d+05
    double precision, parameter :: c71 = -.3964298833d+02
    double precision, parameter :: c72 = -.5019979693d+01
    double precision, parameter :: c73 = 0.2906541488d+00
    double precision, parameter :: c74 = 0.1212943686d+04
    double precision, parameter :: c75 = 0.4140463415d+02
    double precision, parameter :: c81 = 0.1752855549d+03
    double precision, parameter :: c82 = -.2496320107d+02
    double precision, parameter :: c83 = 0.3765413052d+03
    double precision, parameter :: c84 = -.5480488130d+02
    !-- Cbendb values from fit706d.out
    double precision, parameter :: d11 = 0.1917166552d+04
    double precision, parameter :: d12 = -.6542563392d+03
    double precision, parameter :: d13 = 0.6793758367d+02
    double precision, parameter :: d14 = -.1694968847d+04
    double precision, parameter :: d15 = 0.6866649703d+03
    double precision, parameter :: d21 = -.2137567948d+00
    double precision, parameter :: d22 = 0.4975938228d-01
    double precision, parameter :: d23 = 0.9364998295d-01
    double precision, parameter :: d24 = -.2444320779d-01
    double precision, parameter :: d31 = -.2863126914d+02
    double precision, parameter :: d32 = 0.5443219625d-01
    double precision, parameter :: d41 = 0.9673956120d+00
    double precision, parameter :: d42 = -.1160159706d+01
    double precision, parameter :: d43 = -.2424199759d+02
    double precision, parameter :: d44 = 0.8569424490d+01
    double precision, parameter :: d51 = -.6517635862d+04
    double precision, parameter :: d52 = 0.1518098147d+03
    double precision, parameter :: d53 = -.2706514366d+02
    double precision, parameter :: d61 = 0.4308392956d+02
    double precision, parameter :: d62 = -.1234851732d+03
    double precision, parameter :: d63 = 0.2320626055d+03
    double precision, parameter :: d71 = 0.1049541418d+02
    double precision, parameter :: d72 = -.2424169341d+01
    double precision, parameter :: d73 = 0.1745646946d+00
    double precision, parameter :: d74 = 0.2603615561d+02
    double precision, parameter :: d75 = 0.1345799970d+02
    double precision, parameter :: d81 = -.6653710513d+01
    double precision, parameter :: d82 = -.2576854447d+00
    double precision, parameter :: d83 = 0.6172608425d+04
    double precision, parameter :: d84 = -.1328142473d+02
    !   feb27/91 new c51 = c51+c83;   new d51 = d51+d83
    cx1 = c51 + c83
    dx1 = d51 + d83
    r1  = rpass(1)
    r2  = rpass(2)
    r3  = rpass(3)
    t1  = r1*r1 -r2*r2 -r3*r3
    t2  = r2*r2 -r3*r3 -r1*r1
    t3  = r3*r3 -r1*r1 -r2*r2
    !  calculate the cosines of the three internal angles:
    c1     = t1 / (-2.d0*r2*r3)
    c2     = t2 / (-2.d0*r3*r1)
    c3     = t3 / (-2.d0*r1*r2)
    sum    = c1 + c2 + c3
    b1a    = 1.d0 - sum
    c1cube = c1*c1*c1
    c2cube = c2*c2*c2
    c3cube = c3*c3*c3
    cos3t1 = 4.d0*c1cube-3.d0*c1
    cos3t2 = 4.d0*c2cube-3.d0*c2
    cos3t3 = 4.d0*c3cube-3.d0*c3
    sumb   = cos3t1 + cos3t2 + cos3t3
    b1b    = 1.d0 - ( z58*sumb + z38*sum )
    !  calculate derivatives if desired
    dc1dr1  = -r1/(r2*r3)
    dc2dr2  = -r2/(r1*r3)
    dc3dr3  = -r3/(r1*r2)
    dc1dr2  = ( t1/(r2*r2) + 2.d0 )/(2.d0*r3)
    dc1dr3  = ( t1/(r3*r3) + 2.d0 )/(2.d0*r2)
    dc2dr1  = ( t2/(r1*r1) + 2.d0 )/(2.d0*r3)
    dc2dr3  = ( t2/(r3*r3) + 2.d0 )/(2.d0*r1)
    dc3dr1  = ( t3/(r1*r1) + 2.d0 )/(2.d0*r2)
    dc3dr2  = ( t3/(r2*r2) + 2.d0 )/(2.d0*r1)
    db1a(1) = -1.d0*( dc1dr1 + dc2dr1 + dc3dr1 )
    db1a(2) = -1.d0*( dc1dr2 + dc2dr2 + dc3dr2 )
    db1a(3) = -1.d0*( dc1dr3 + dc2dr3 + dc3dr3 )
    d1      = 12.d0*c1*c1 - 3.d0
    d2      = 12.d0*c2*c2 - 3.d0
    d3      = 12.d0*c3*c3 - 3.d0
    db1b(1) = -z58*( d1*dc1dr1 + d2*dc2dr1 + d3*dc3dr1 ) -z38*(    dc1dr1 +    dc2dr1 +    dc3dr1 )
    db1b(2) = -z58*( d1*dc1dr2 + d2*dc2dr2 + d3*dc3dr2 ) -z38*(    dc1dr2 +    dc2dr2 +    dc3dr2 )
    db1b(3) = -z58*( d1*dc1dr3 + d2*dc2dr3 + d3*dc3dr3 ) -z38*(    dc1dr3 +    dc2dr3 +    dc3dr3 )
    !d     if(ipr.gt.1)then
    !d         write(7,*)
    !d         write(7,*) ' from subr.vbcb -----------------------------'
    !d         write(7,6000) rpass
    !d         write(7,6010) b1a,b1b
    !d         write(7,6020) db1a,db1b
    !d         if(ipr.gt.2)then
    !d            write(7,7000) dc1dr1,dc1dr2,dc1dr3,
    !d    .                     dc2dr1,dc2dr2,dc2dr3,
    !d    .                     dc3dr1,dc3dr2,dc3dr3
    !d            write(7,7010) t1,t2,t3,c1,c2,c3,d1,d2,d3
    !d         end if
    !d     end if
    !
    !  calculate the quantities used by both Vbend and Cbend:
    !
    r      = r1 + r2 + r3
    rsq    = r*r
    b2     = 1.d0/r1 + 1.d0/r2 + 1.d0/r3
    b3     = (r2-r1)*(r2-r1) + (r3-r2)*(r3-r2) + (r1-r3)*(r1-r3)
    eps2   = 1.d-12
    b3b    = sqrt( b3 + eps2 )
    db2(1) = -1.d0/(r1*r1)
    db2(2) = -1.d0/(r2*r2)
    db2(3) = -1.d0/(r3*r3)
    db3(1) = 4.d0*r1 -2.d0*r2 -2.d0*r3
    db3(2) = 4.d0*r2 -2.d0*r3 -2.d0*r1
    db3(3) = 4.d0*r3 -2.d0*r1 -2.d0*r2
    db3(4) = 0.5d0*db3(1)/b3b
    db3(5) = 0.5d0*db3(2)/b3b
    db3(6) = 0.5d0*db3(3)/b3b
    exp1   = exp(-beta1*r)
    exp2   = exp(-beta2*rsq)
    exp7   = exp(-beta3*r)
    dexp1  = -beta1*exp1
    dexp2  = -2.d0*beta2*r*exp2
    dexp7  = -beta3*exp7
    !   do the vbnda calculations:
    b1      = b1a
    b12     = b1*b1
    b13     = b12*b1
    b14     = b13*b1
    b15     = b14*b1
    asum    = a11 +  a12*r +  a13*rsq
    bsum    = a21*b12 +  a22*b13 +  a23*b14 +  a24*b15
    csum    = a31*b1*exp1 +  a32*b12*exp2
    dsum1   = a41*exp1 +  a42*exp2
    dsum2   = a43*exp1 +  a44*exp2
    fsum    = a51 +  a52*r +  a53*rsq
    vb(1,1) = b1*asum*exp1
    vb(1,2) = bsum * exp2
    vb(1,3) = b2*csum
    vb(1,4) = b1*b3*dsum1 + b1*b3b*dsum2
    vb(1,5) = b1* fsum  *exp7
    vbnd(1) = vb(1,1)+vb(1,2)+vb(1,3)+vb(1,4)+vb(1,5)
    !
    dasum   = a12+2.d0*a13*r
    dbsum   = 2.d0* a21*b1  + 3.d0* a22*b12 + 4.d0* a23*b13 + 5.d0* a24*b14
    ddsum1  = a41*dexp1 +  a42*dexp2
    ddsum2  = a43*dexp1 +  a44*dexp2
    dfsum   = a52 +2.d0*a53*r
    dvb1(1) = db1a(1)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
    dvb1(2) = db1a(2)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
    dvb1(3) = db1a(3)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
    !
    dvb2(1) = dbsum*db1a(1)*exp2 + bsum*dexp2
    dvb2(2) = dbsum*db1a(2)*exp2 + bsum*dexp2
    dvb2(3) = dbsum*db1a(3)*exp2 + bsum*dexp2
    !
    ! calculate the Vb3 derivatives:
    dvb3(1) = db2(1)*csum + b2*(  a31*db1a(1)*exp1 +  a31*b1*dexp1 +2.d0* a32*b1*db1a(1)*exp2 +  a32*b12*dexp2)
    dvb3(2) = db2(2)*csum + b2*(  a31*db1a(2)*exp1 +  a31*b1*dexp1 +2.d0* a32*b1*db1a(2)*exp2 +  a32*b12*dexp2)
    dvb3(3) = db2(3)*csum + b2*(  a31*db1a(3)*exp1 +  a31*b1*dexp1 +2.d0* a32*b1*db1a(3)*exp2 +  a32*b12*dexp2)
    !
    ! calculate the Vb4 derivatives (may 27/90):
    dvb4(1) = db1a(1)*b3 *dsum1+  b1*db3(1)*dsum1 + b1*b3 *ddsum1 + db1a(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1*b3b*ddsum2
    dvb4(2) = db1a(2)*b3 *dsum1+  b1*db3(2)*dsum1 + b1*b3 *ddsum1 + db1a(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1*b3b*ddsum2
    dvb4(3) = db1a(3)*b3 *dsum1+  b1*db3(3)*dsum1 + b1*b3 *ddsum1 + db1a(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1*b3b*ddsum2
    !
    ! calculate the Vb5 derivatives:
    dvb5(1) = db1a(1)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
    dvb5(2) = db1a(2)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
    dvb5(3) = db1a(3)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
    !
    ! calculate the overall derivatives:
    dvbnda(1) = dvb1(1)+dvb2(1)+dvb3(1)+dvb4(1)+dvb5(1)
    dvbnda(2) = dvb1(2)+dvb2(2)+dvb3(2)+dvb4(2)+dvb5(2)
    dvbnda(3) = dvb1(3)+dvb2(3)+dvb3(3)+dvb4(3)+dvb5(3)
    vbnda     = vbnd(1)
    !
    !---- now do the vbendb calculations --------
    b1      = b1b
    b12     = b1*b1
    b13     = b12*b1
    b14     = b13*b1
    b15     = b14*b1
    asum    = g11 +  g12*r +  g13*rsq
    bsum    = g21*b12 +  g22*b13 +  g23*b14 +  g24*b15
    csum    = g31*b1*exp1 +  g32*b12*exp2
    dsum1   = g41*exp1 +  g42*exp2
    dsum2   = g43*exp1 +  g44*exp2
    fsum    = g51 +  g52*r +  g53*rsq
    vb(2,1) = b1*asum*exp1
    vb(2,2) = bsum * exp2
    vb(2,3) = b2*csum
    vb(2,4) = b1*b3*dsum1 + b1*b3b*dsum2
    vb(2,5) = b1* fsum  *exp7
    vbnd(2) = vb(2,1)+vb(2,2)+vb(2,3)+vb(2,4)+vb(2,5)
    !
    dasum   = g12+2.d0*g13*r
    dbsum   = 2.d0* g21*b1  + 3.d0* g22*b12 +4.d0* g23*b13 + 5.d0* g24*b14
    ddsum1  = g41*dexp1 +  g42*dexp2
    ddsum2  = g43*dexp1 +  g44*dexp2
    dfsum   = g52 +2.d0*g53*r
    dvb1(1) = db1b(1)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
    dvb1(2) = db1b(2)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
    dvb1(3) = db1b(3)*asum*exp1 + b1*dasum*exp1 +b1*asum*dexp1
    !
    dvb2(1) = dbsum*db1b(1)*exp2 + bsum*dexp2
    dvb2(2) = dbsum*db1b(2)*exp2 + bsum*dexp2
    dvb2(3) = dbsum*db1b(3)*exp2 + bsum*dexp2
    !
    ! calculate the Vb3 derivatives:
    dvb3(1) = db2(1)*csum + b2*(  g31*db1b(1)*exp1 +  g31*b1*dexp1 +2.d0* g32*b1*db1b(1)*exp2 +  g32*b12*dexp2)
    dvb3(2) = db2(2)*csum + b2*(  g31*db1b(2)*exp1 +  g31*b1*dexp1 +2.d0* g32*b1*db1b(2)*exp2 +  g32*b12*dexp2)
    dvb3(3) = db2(3)*csum + b2*(  g31*db1b(3)*exp1 +  g31*b1*dexp1 +2.d0* g32*b1*db1b(3)*exp2 +  g32*b12*dexp2)
    !
    ! calculate the Vb4 derivatives (may 27/90):
    dvb4(1) = db1b(1)*b3 *dsum1+  b1*db3(1)*dsum1 + b1*b3 *ddsum1 + db1b(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1*b3b*ddsum2
    dvb4(2) = db1b(2)*b3 *dsum1+  b1*db3(2)*dsum1 + b1*b3 *ddsum1 + db1b(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1*b3b*ddsum2
    dvb4(3) = db1b(3)*b3 *dsum1+  b1*db3(3)*dsum1 + b1*b3 *ddsum1 + db1b(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1*b3b*ddsum2
    !
    ! calculate the Vb5 derivatives:
    dvb5(1) = db1b(1)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
    dvb5(2) = db1b(2)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
    dvb5(3) = db1b(3)*fsum*exp7 + b1*dfsum*exp7 + b1*fsum*dexp7
    !
    ! calculate the overall derivatives:
    dvbndb(1) = dvb1(1)+dvb2(1)+dvb3(1)+dvb4(1)+dvb5(1)
    dvbndb(2) = dvb1(2)+dvb2(2)+dvb3(2)+dvb4(2)+dvb5(2)
    dvbndb(3) = dvb1(3)+dvb2(3)+dvb3(3)+dvb4(3)+dvb5(3)
    vbndb     = vbnd(2)
    !
    !  now calculate the Cbend correction terms:
    !  Cbnda uses the B1a formula and Cbndb uses the B1b formula
    if(icompc.eq.0) return
    sumt      = t(1) + t(2) + t(3)
    rcu       = rsq*r
    p         = r1*r2*r3
    psq       = p*p
    pcu       = psq*p
    dp(1)     = r2*r3
    dp(2)     = r3*r1
    dp(3)     = r1*r2
    cbnds(1)  = 0.d0
    cbnds(2)  = 0.d0
    dcbnda(1) = 0.d0
    dcbnda(2) = 0.d0
    dcbnda(3) = 0.d0
    dcbndb(1) = 0.d0
    dcbndb(2) = 0.d0
    dcbndb(3) = 0.d0
    !  calculate exponentials and derivatives (common to cbnda and cbndb):
    exp7  = exp( -beta2*rcu )
    dexp7 = -3.d0*beta2*rsq*exp7
    ! ----- calculate the Cbnda correction term ----------------
    b1       = b1a
    b12      = b1*b1
    b13      = b12*b1
    b14      = b13*b1
    b15      = b14*b1
    asum     = c11 + c12*r + c13*rsq + c14/r + c15/rsq
    bsum     = c21*b12 + c22*b13 + c23*b14 + c24*b15
    csum     = c31*b1 *exp1  + c32*b12*exp2
    dsum1    = c41*exp1 + c42*exp2
    dsum2    = c43*exp1 + c44*exp2
    fsum     = cx1 + c52*r + c53*rsq
    gsum     = c61 + c62/r + c63/rsq
    aasum    = c71 + c72*p + c73*psq + c74/p + c75/psq
    ffsum    = c81*p + c82*psq + c84/psq
    dasum    = c12 + 2.d0*c13*r -c14/rsq -2.d0*c15/rcu
    dbsum    = 2.d0*c21*b1 +3.d0*c22*b12 +4.d0*c23*b13 +5.d0*c24*b14
    ddsum1   = c41*dexp1 + c42*dexp2
    ddsum2   = c43*dexp1 + c44*dexp2
    dfsum    = c52 + 2.d0*c53*r
    dgsum    = -c62/rsq - 2.d0*c63/rcu
    daasum   = c72 +2.d0*c73*p -c74/psq -2.d0*c75/pcu
    dffsum   = c81 + 2.d0*c82*p -2.d0*c84/pcu
    cb(1,1)  = b1*asum * exp1  /p
    cb(1,2)  = bsum * exp2
    cb(1,3)  = b2*csum
    cb(1,4)  = b1*b3*dsum1 + b1*b3b*dsum2
    cb(1,5)  = b1 * fsum * exp7 /p
    cb(1,6)  = b1 * gsum * exp7
    cb(1,7)  = b1 * aasum * exp2
    cb(1,8)  = b1 * ffsum * exp7
    cbnds(1) = cb(1,1)+cb(1,2)+cb(1,3)+cb(1,4)+cb(1,5) +cb(1,6)+cb(1,7)+cb(1,8)
    !       calculate the derivatives:
    dcb1(1) = db1a(1)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(1)/psq
    dcb1(2) = db1a(2)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(2)/psq
    dcb1(3) = db1a(3)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(3)/psq
    !
    dcb2(1) = dbsum*db1a(1)*exp2 + bsum*dexp2
    dcb2(2) = dbsum*db1a(2)*exp2 + bsum*dexp2
    dcb2(3) = dbsum*db1a(3)*exp2 + bsum*dexp2
    !
    dcb3(1) = db2(1)*csum + b2*( c31*db1a(1)*exp1 +c31*b1*dexp1 +c32 *2.d0*b1*db1a(1)*exp2 + c32 *b12*dexp2 )
    dcb3(2) = db2(2)*csum + b2*( c31 *db1a(2)*exp1 +c31 *b1*dexp1 +c32 *2.d0*b1*db1a(2)*exp2 + c32 *b12*dexp2 )
    dcb3(3) = db2(3)*csum + b2*( c31 *db1a(3)*exp1 +c31 *b1*dexp1 +c32 *2.d0*b1*db1a(3)*exp2 + c32 *b12*dexp2 )
    !          dCb4 equations may 27/90:
    dcb4(1) = db1a(1)*b3 *dsum1 + b1*db3(1)*dsum1 + b1*b3 *ddsum1 + db1a(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1*b3b*ddsum2
    dcb4(2) = db1a(2)*b3 *dsum1 + b1*db3(2)*dsum1 + b1*b3 *ddsum1 + db1a(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1*b3b*ddsum2
    dcb4(3) = db1a(3)*b3 *dsum1 + b1*db3(3)*dsum1 + b1*b3 *ddsum1 + db1a(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1*b3b*ddsum2
    !          dCb5 equations may 27/90: (corrected jun01)
    dcb5(1) = db1a(1)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(1)/psq
    dcb5(2) = db1a(2)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(2)/psq
    dcb5(3) = db1a(3)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(3)/psq
    !
    dcb6(1) = db1a(1)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
    dcb6(2) = db1a(2)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
    dcb6(3) = db1a(3)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
    !
    dcb7(1)     = db1a(1)*aasum*exp2 + b1*daasum*dp(1)*exp2 + b1* aasum*dexp2
    dcb7(2)     = db1a(2)*aasum*exp2 + b1*daasum*dp(2)*exp2 + b1* aasum*dexp2
    dcb7(3)     = db1a(3)*aasum*exp2 + b1*daasum*dp(3)*exp2 + b1* aasum*dexp2
    dcb8(1)     = db1a(1)*ffsum*exp7  +b1*dffsum*dp(1)*exp7 +b1* ffsum*dexp7
    dcb8(2)     = db1a(2)*ffsum*exp7  +b1*dffsum*dp(2)*exp7 +b1* ffsum*dexp7
    dcb8(3)     = db1a(3)*ffsum*exp7  +b1*dffsum*dp(3)*exp7 +b1* ffsum*dexp7
    dcbnds(1,1) = dcb1(1) +dcb2(1) +dcb3(1) +dcb4(1) +dcb5(1) +dcb6(1) +dcb7(1) +dcb8(1)
    dcbnds(1,2) = dcb1(2) +dcb2(2) +dcb3(2) +dcb4(2) +dcb5(2) +dcb6(2) +dcb7(2) +dcb8(2)
    dcbnds(1,3) = dcb1(3) +dcb2(3) +dcb3(3) +dcb4(3) +dcb5(3) +dcb6(3) +dcb7(3) +dcb8(3)
    !
    ! ----- calculate the Cbndb correction term ----------------
    b1       = b1b
    b12      = b1*b1
    b13      = b12*b1
    b14      = b13*b1
    b15      = b14*b1
    asum     = d11 + d12*r + d13*rsq + d14/r + d15/rsq
    bsum     = d21*b12 + d22*b13 + d23*b14 + d24*b15
    csum     = d31*b1 *exp1  + d32*b12*exp2
    dsum1    = d41*exp1 + d42*exp2
    dsum2    = d43*exp1 + d44*exp2
    fsum     = dx1 + d52*r + d53*rsq
    gsum     = d61 + d62/r + d63/rsq
    aasum    = d71 + d72*p + d73*psq + d74/p + d75/psq
    ffsum    = d81*p + d82*psq + d84/psq
    dasum    = d12 + 2.d0*d13*r -d14/rsq -2.d0*d15/rcu
    dbsum    = 2.d0*d21*b1 +3.d0*d22*b12 +4.d0*d23*b13 +5.d0*d24*b14
    ddsum1   = d41*dexp1 + d42*dexp2
    ddsum2   = d43*dexp1 + d44*dexp2
    dfsum    = d52 + 2.d0*d53*r
    dgsum    = -d62/rsq - 2.d0*d63/rcu
    daasum   = d72 +2.d0*d73*p -d74/psq -2.d0*d75/pcu
    dffsum   = d81 + 2.d0*d82*p -2.d0*d84/pcu
    cb(2,1)  = b1*asum * exp1  /p
    cb(2,2)  = bsum * exp2
    cb(2,3)  = b2*csum
    cb(2,4)  = b1*b3*dsum1 + b1*b3b*dsum2
    cb(2,5)  = b1 * fsum * exp7 /p
    cb(2,6)  = b1 * gsum * exp7
    cb(2,7)  = b1 * aasum * exp2
    cb(2,8)  = b1 * ffsum * exp7
    cbnds(2) = cb(2,1)+cb(2,2)+cb(2,3)+cb(2,4)+cb(2,5) +cb(2,6)+cb(2,7)+cb(2,8)
    !       calculate the derivatives:
    dcb1(1) = db1b(1)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(1)/psq
    dcb1(2) = db1b(2)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(2)/psq
    dcb1(3) = db1b(3)*asum*exp1/p + b1*dasum*exp1/p +b1*asum*dexp1/p - b1*asum*exp1*dp(3)/psq
    !
    dcb2(1) = dbsum*db1b(1)*exp2 + bsum*dexp2
    dcb2(2) = dbsum*db1b(2)*exp2 + bsum*dexp2
    dcb2(3) = dbsum*db1b(3)*exp2 + bsum*dexp2
    !
    dcb3(1) = db2(1)*csum + b2*( d31 *db1b(1)*exp1 +d31 *b1*dexp1 +d32 *2.d0*b1*db1b(1)*exp2 + d32 *b12*dexp2 )
    dcb3(2) = db2(2)*csum + b2*( d31 *db1b(2)*exp1 +d31 *b1*dexp1 +d32 *2.d0*b1*db1b(2)*exp2 + d32 *b12*dexp2 )
    dcb3(3) = db2(3)*csum + b2*( d31 *db1b(3)*exp1 +d31 *b1*dexp1 +d32 *2.d0*b1*db1b(3)*exp2 + d32 *b12*dexp2 )
    !          dCb4 equations may 27/90:
    dcb4(1) = db1b(1)*b3 *dsum1 + b1*db3(1)*dsum1 + b1*b3 *ddsum1 + db1b(1)*b3b*dsum2 + b1*db3(4)*dsum2 + b1*b3b*ddsum2
    dcb4(2) = db1b(2)*b3 *dsum1 + b1*db3(2)*dsum1 + b1*b3 *ddsum1 + db1b(2)*b3b*dsum2 + b1*db3(5)*dsum2 + b1*b3b*ddsum2
    dcb4(3) = db1b(3)*b3 *dsum1 + b1*db3(3)*dsum1 + b1*b3 *ddsum1 + db1b(3)*b3b*dsum2 + b1*db3(6)*dsum2 + b1*b3b*ddsum2
    !          dCb5 equations may 27/90:
    dcb5(1) = db1b(1)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(1)/psq
    dcb5(2) = db1b(2)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(2)/psq
    dcb5(3) = db1b(3)*fsum*exp7/p + b1*dfsum*exp7/p + b1*fsum*dexp7/p - b1*fsum*exp7*dp(3)/psq
    !
    dcb6(1) = db1b(1)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
    dcb6(2) = db1b(2)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
    dcb6(3) = db1b(3)*gsum*exp7 + b1*dgsum*exp7 + b1*gsum*dexp7
    !
    dcb7(1)     = db1b(1)*aasum*exp2 + b1*daasum*dp(1)*exp2 + b1* aasum*dexp2
    dcb7(2)     = db1b(2)*aasum*exp2 + b1*daasum*dp(2)*exp2 + b1* aasum*dexp2
    dcb7(3)     = db1b(3)*aasum*exp2 + b1*daasum*dp(3)*exp2 + b1* aasum*dexp2
    dcb8(1)     = db1b(1)*ffsum*exp7  +b1*dffsum*dp(1)*exp7 +b1* ffsum*dexp7
    dcb8(2)     = db1b(2)*ffsum*exp7  +b1*dffsum*dp(2)*exp7 +b1* ffsum*dexp7
    dcb8(3)     = db1b(3)*ffsum*exp7  +b1*dffsum*dp(3)*exp7 +b1* ffsum*dexp7
    dcbnds(2,1) = dcb1(1) +dcb2(1) +dcb3(1) +dcb4(1) +dcb5(1) +dcb6(1) +dcb7(1) +dcb8(1)
    dcbnds(2,2) = dcb1(2) +dcb2(2) +dcb3(2) +dcb4(2) +dcb5(2) +dcb6(2) +dcb7(2) +dcb8(2)
    dcbnds(2,3) = dcb1(3) +dcb2(3) +dcb3(3) +dcb4(3) +dcb5(3) +dcb6(3) +dcb7(3) +dcb8(3)
    !
    !  calculate the total derivative from the pieces:
    dcbnda(1) = dt(1) * cbnds(1)    + sumt * dcbnds(1,1)
    dcbndb(1) = dt(1) * cbnds(2)    + sumt * dcbnds(2,1)
    dcbnda(2) = dt(2) * cbnds(1)    + sumt * dcbnds(1,2)
    dcbndb(2) = dt(2) * cbnds(2)    + sumt * dcbnds(2,2)
    dcbnda(3) = dt(3) * cbnds(1)    + sumt * dcbnds(1,3)
    dcbndb(3) = dt(3) * cbnds(2)    + sumt * dcbnds(2,3)
    cbnda     = sumt*cbnds(1)
    cbndb     = sumt*cbnds(2)
    !  some debugging print statements:
    !d     if(ipr.gt.1)then
    !d        write(7,7650) B2,B3
    !d        write(7,7700) db2,db3
    !d        write(7,*)
    !d        write(7,*) 'Vbend values :'
    !d        k=1
    !d        write(7,7150) k,dvb1,dvb2,dvb3,dvb4,dvb5,dVbnda
    !d        write(7,7155) (vb(1,i),i=1,5)
    !d        k=2
    !d        write(7,7150) k,dvb1,dvb2,dvb3,dvb4,dvb5,dVbndb
    !d        write(7,7155) (vb(2,i),i=1,5)
    !d        write(7,7100) Vbnda,Vbndb
    !d        write(7,*)
    !d        write(7,*) 'Cbend values:'
    !d        write(7,*) 'c81,c81,c82=',c81,c81,c82
    !d        write(7,*) 'c91,c83,c84=',c91,c83,c84
    !d        write(7,*) 'd81,d81,d82=',d81,d81,d82
    !d        write(7,*) 'd91,d83,d84=',d91,d83,d84
    !d        k=1
    !d           write(7,2020) k,(cb(1,j),j=1,9)
    !d           write(7,2033) sumT,Cbnds(1)
    !d        k=2
    !d           write(7,2020) k,(cb(k,j),j=1,9)
    !d           write(7,2033) sumT,Cbnds(2)
    !d        k=1
    !d        write(7,7152) k,dcb1,dcb2,dcb3,dcb4,dcb5,dcb6,
    !d    .               dcb7,dcb8,dcb9,(dCbnds(1,i),i=1,3)
    !d        k=2
    !d        write(7,7152) k,dcb1,dcb2,dcb3,dcb4,dcb5,dcb6,
    !d    .                 dcb7,dcb8,dcb9,(dCbnds(k,i),i=1,3)
    !d        write(7,7200) dT,dCbnda,dCbndb
    !d        write(7,*) 'VbCb  exit -------------------------'
    !d        write(7,*)
    !d     end if
end subroutine vbcb
!
