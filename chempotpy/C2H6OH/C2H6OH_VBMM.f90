!   System: C2H6 + OH --> C2H5 + H2O.
!              Joaquin Espinosa-Garcia Mayo 2019
!   All the information passed to and from the potential energy surface
!   routine is in hartree atomic units.
!        This potential is written such that:
!                       X(1)  - X(3)  : X, Y, Z for C1
!                       X(4)  - X(6)  : X, Y, Z for C2
!                       X(7)  - X(9)  : X, Y, Z for H1
!                       X(10) - X(12) : X, Y, Z for H2
!                       X(13) - X(15) : X, Y, Z for H3
!                       X(16) - X(18) : X, Y, Z for H4
!                       X(19) - X(21) : X, Y, Z for H5
!                       X(22) - X(24) : X, Y, Z for H6
!                       X(25) - X(27) : X, Y, Z for O
!                       X(28) - X(30) : X, Y, Z for HO

!> @brief Implement C2H6OH_VBMM, a potential energy surface for C2H6OH.
module chempotpy_c2h6oh_vbmm_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the C2H6OH_VBMM potential energy surface for C2H6OH.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_legacy_units, only: legacy_hartree_to_electronvolt
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 10
        integer, parameter            :: natom   = 25
        integer, parameter            :: isurf   = 5
        integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter            :: n3atom  = 75
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: a1c
        double precision              :: a1c1h
        double precision              :: a1c2h
        double precision              :: a1cc
        double precision              :: a1s
        double precision              :: a1scc
        double precision              :: a2c
        double precision              :: a2s
        double precision              :: a2scc
        double precision              :: a3s
        double precision              :: a3scc
        double precision              :: aa1
        double precision              :: aa2
        double precision              :: aa3
        double precision              :: aa4
        double precision              :: acb
        double precision              :: acc
        double precision              :: ahh
        double precision              :: ak
        double precision              :: alph2o
        double precision              :: angh2o(6)
        double precision              :: angh2oeq
        double precision              :: anuzero
        double precision              :: aphi
        double precision              :: atheta
        double precision              :: b1c1h
        double precision              :: b1c2h
        double precision              :: b1cc
        double precision              :: b1s
        double precision              :: b1scc
        double precision              :: b2s
        double precision              :: b2scc
        double precision              :: b3s
        double precision              :: b3scc
        double precision              :: bk
        double precision              :: bphi
        double precision              :: btheta
        double precision              :: c1c1h
        double precision              :: c1c2h
        double precision              :: c1cc
        double precision              :: cart(natom,3)
        double precision              :: cartnu(natom,3)
        double precision              :: cnvrtd
        double precision              :: cnvrtde
        double precision              :: cnvrte
        double precision              :: cphi
        double precision              :: ctheta
        double precision              :: d1c1h
        double precision              :: d1c2h
        double precision              :: d1cb
        double precision              :: d1cc
        double precision              :: d1ch
        double precision              :: d1hh
        double precision              :: d3c1h
        double precision              :: d3c2h
        double precision              :: d3cb
        double precision              :: d3ch
        double precision              :: d3hh
        double precision              :: deesdr(n3atom,isurf)
        double precision              :: degsdr(n3atom)
        double precision              :: deijdr(n3atom,jsurf)
        double precision              :: descart(natom,3,isurf)
        double precision              :: descartnu(natom,3,isurf)
        double precision              :: dfdc(8,8,8)
        double precision              :: dfdelta(8,8)
        double precision              :: dfdh(8,8,8)
        double precision              :: dgscart(natom,3)
        double precision              :: dgscartnu(natom,3)
        double precision              :: dhdelta(8,8)
        double precision              :: dijcart(natom,3,jsurf)
        double precision              :: dijcartnu(natom,3,jsurf)
        double precision              :: ds1(8)
        double precision              :: ds2(8)
        double precision              :: ds3(8)
        double precision              :: dsphi(8)
        double precision              :: dstheta(8)
        double precision              :: dtheta0(8,8,8)
        double precision              :: engyes(isurf)
        double precision              :: engygs
        double precision              :: engyij(jsurf)
        double precision              :: ezero(isurf+1)
        double precision              :: f1(8)
        double precision              :: fch3
        double precision              :: fdelta(8)
        double precision              :: fk0(8,8)
        double precision              :: fkh2o(6)
        double precision              :: fkh2oeq
        double precision              :: fkinf
        double precision              :: hch3
        double precision              :: hdelta(8)
        integer                       :: icartr
        integer                       :: indexes(natom)
        integer                       :: irctnt
        integer                       :: ireorder
        integer                       :: kediag
        integer                       :: keoffd
        integer                       :: ksdiag
        integer                       :: ksoffd
        integer                       :: legacy_natoms
        integer                       :: mder
        integer                       :: msurf
        integer                       :: nasurf(isurf+1,isurf+1)
        integer                       :: nc(3)
        integer                       :: nc1(3)
        integer                       :: nc2(3)
        integer                       :: nder
        integer                       :: nflag(20)
        integer                       :: nh(8,3)
        integer                       :: nhb(3)
        integer                       :: nnb
        integer                       :: nnc1
        integer                       :: nnc2
        integer                       :: nnh(6)
        integer                       :: nno
        integer                       :: no(3)
        integer                       :: nulbl(natom)
        double precision              :: pdot(30)
        double precision              :: pengyes(isurf)
        double precision              :: pengygs
        double precision              :: pengyij(jsurf)
        double precision              :: q(30)
        double precision              :: r(n3atom)
        double precision              :: r0c1h
        double precision              :: r0c1hp
        double precision              :: r0c1hr
        double precision              :: r0c2h
        double precision              :: r0c2hp
        double precision              :: r0c2hr
        double precision              :: r0cb
        double precision              :: r0cc
        double precision              :: r0hh
        double precision              :: rbh(6)
        double precision              :: rc1b
        double precision              :: rc1h(3)
        double precision              :: rc2b
        double precision              :: rc2h(3)
        double precision              :: rcc
        character                     :: ref(5)
        double precision              :: rhh(3,6)
        double precision              :: rno
        double precision              :: s1(8)
        double precision              :: s2(8)
        double precision              :: s3(8)
        double precision              :: sphi(8)
        double precision              :: stheta(8)
        double precision              :: tbh(6,3)
        double precision              :: tc1b(3)
        double precision              :: tc1h(3,3)
        double precision              :: tc2b(3)
        double precision              :: tc2h(3,3)
        double precision              :: tcb(3)
        double precision              :: tcc(3)
        double precision              :: theta0(8,8)
        double precision              :: thetahcch(9)
        double precision              :: thh(3,6,3)
        double precision              :: tno(3)
        double precision              :: v3
        double precision              :: w1
        double precision              :: w2
        double precision              :: w3
        double precision              :: w4
        double precision              :: w5
        double precision              :: w6
        double precision              :: v
        integer                       :: iatom, idir, istate
        !initialize
        a1c       = 0.d0; a1c1h = 0.d0; a1c2h = 0.d0; a1cc = 0.d0; a1s = 0.d0; a1scc = 0.d0; a2c = 0.d0; a2s = 0.d0
        a2scc     = 0.d0; a3s = 0.d0; a3scc = 0.d0; aa1 = 0.d0; aa2 = 0.d0; aa3 = 0.d0; aa4 = 0.d0; acb = 0.d0; acc = 0.d0
        ahh       = 0.d0; ak = 0.d0; alph2o = 0.d0; angh2o = 0.d0; angh2oeq = 0.d0; anuzero = 0.d0; aphi = 0.d0; atheta = 0.d0
        b1c1h     = 0.d0; b1c2h = 0.d0; b1cc = 0.d0; b1s = 0.d0; b1scc = 0.d0; b2s = 0.d0; b2scc = 0.d0; b3s = 0.d0
        b3scc     = 0.d0; bk = 0.d0; bphi = 0.d0; btheta = 0.d0; c1c1h = 0.d0; c1c2h = 0.d0; c1cc = 0.d0; cart = 0.d0
        cartnu    = 0.d0; cnvrtd = 0.d0; cnvrtde = 0.d0; cnvrte = 0.d0; cphi = 0.d0; ctheta = 0.d0; d1c1h = 0.d0
        d1c2h     = 0.d0; d1cb = 0.d0; d1cc = 0.d0; d1ch = 0.d0; d1hh = 0.d0; d3c1h = 0.d0; d3c2h = 0.d0; d3cb = 0.d0
        d3ch      = 0.d0; d3hh = 0.d0; deesdr = 0.d0; degsdr = 0.d0; deijdr = 0.d0; descart = 0.d0; descartnu = 0.d0
        dfdc      = 0.d0; dfdelta = 0.d0; dfdh = 0.d0; dgscart = 0.d0; dgscartnu = 0.d0; dhdelta = 0.d0; dijcart = 0.d0
        dijcartnu = 0.d0; ds1 = 0.d0; ds2 = 0.d0; ds3 = 0.d0; dsphi = 0.d0; dstheta = 0.d0; dtheta0 = 0.d0
        engyes    = 0.d0; engygs = 0.d0; engyij = 0.d0; ezero = 0.d0; f1 = 0.d0; fch3 = 0.d0; fdelta = 0.d0; fk0 = 0.d0
        fkh2o     = 0.d0; fkh2oeq = 0.d0; fkinf = 0.d0; hch3 = 0.d0; hdelta = 0.d0; icartr = 0; indexes = 0; irctnt = 0
        ireorder  = 0; kediag = 0; keoffd = 0; ksdiag = 0; ksoffd = 0; legacy_natoms = 0; mder = 0; msurf = 0; nasurf = 0
        nc        = 0; nc1 = 0; nc2 = 0; nder = 0; nflag = 0; nh = 0; nhb = 0; nnb = 0; nnc1 = 0; nnc2 = 0; nnh = 0; nno = 0
        no        = 0; nulbl = 0; pdot = 0.d0; pengyes = 0.d0; pengygs = 0.d0; pengyij = 0.d0; q = 0.d0; r = 0.d0; r0c1h = 0.d0
        r0c1hp    = 0.d0; r0c1hr = 0.d0; r0c2h = 0.d0; r0c2hp = 0.d0; r0c2hr = 0.d0; r0cb = 0.d0; r0cc = 0.d0; r0hh = 0.d0
        rbh       = 0.d0; rc1b = 0.d0; rc1h = 0.d0; rc2b = 0.d0; rc2h = 0.d0; rcc = 0.d0
        ref       = ' '
        rhh       = 0.d0
        rno       = 0.d0
        s1        = 0.d0
        s2        = 0.d0
        s3        = 0.d0
        sphi      = 0.d0
        stheta    = 0.d0
        tbh       = 0.d0
        tc1b      = 0.d0
        tc1h      = 0.d0
        tc2b      = 0.d0
        tc2h      = 0.d0
        tcb       = 0.d0
        tcc       = 0.d0
        theta0    = 0.d0
        thetahcch = 0.d0
        thh       = 0.d0
        tno       = 0.d0
        v3        = 0.d0
        w1        = 0.d0
        w2        = 0.d0
        w3        = 0.d0
        w4        = 0.d0
        w5        = 0.d0
        w6        = 0.d0
        nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0], shape(nasurf))
        nder          = 1
        nflag         = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        anuzero       = 0.0d0
        icartr        = 1
        msurf         = 0
        mder          = 1
        nulbl         = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        legacy_natoms = 10
        nnc1          = 1
        nnc2          = 2
        nnb           = 9
        nnh           = [3, 4, 5, 6, 7, 8]
        nno           = 10
        r0c1hr        = 1.08897d0
        r0c1hp        = 1.07577d0
        w3            = 1.00000d0
        w4            = 1.08897d0
        d1c1h         = 109.23000d0
        d3c1h         = 22.00400d0
        a1c1h         = 1.68000d0
        b1c1h         = 0.14500d0
        c1c1h         = 8.21404d0
        r0c2hr        = 1.08897d0
        r0c2hp        = 1.07577d0
        w5            = 1.00000d0
        w6            = 1.08897d0
        d1c2h         = 109.23000d0
        d3c2h         = 22.00400d0
        a1c2h         = 1.68000d0
        b1c2h         = 0.14500d0
        c1c2h         = 8.21404d0
        r0cc          = 1.50007d0
        d1cc          = 117.94900d0
        acc           = 1.93400d0
        a1c           = 0.90007d0
        a2c           = 0.94900d0
        r0hh          = 0.97141d0
        d1hh          = 125.85800d0
        d3hh          = 29.56400d0
        ahh           = 2.10570d0
        r0cb          = 1.90397d0
        d1cb          = 52.50900d0
        d3cb          = 14.86000d0
        acb           = 0.8310285d0
        a3s           = 1.8419147d0
        b3s           = 1.0889300d0
        aphi          = 1.0087903d0
        bphi          = 1.2006638d0
        cphi          = 1.0889937d0
        atheta        = 1.2078714d0
        btheta        = 0.3548859d0
        ctheta        = 1.8915497d0
        fch3          = 0.1357500d0
        hch3          = 0.2615000d0
        fkinf         = 0.2977000d0
        ak            = 0.3060000d0
        bk            = 0.2032000d0
        aa1           = 0.800952d0
        aa2           = 0.499963d0
        aa3           = 6.465953d0
        aa4           = 1.569977d0
        fkh2oeq       = 0.6500000d0
        alph2o        = 3.2080000d0
        angh2oeq      = 104.7132000d0
        v3            = 13.730000d0
        w1            = 9.3548859d0
        w2            = 0.913970d0
        a1s           = 2.5313681d0
        b1s           = 1.088625d0
        a2s           = 2.0147402d0
        b2s           = 1.088980d0
        a1scc         = 2.5313681d0
        b1scc         = 1.550070d0
        a2scc         = 2.0147402d0
        b2scc         = 1.550070d0
        a3scc         = 0.531368d0
        b3scc         = 1.550070d0
        v             = 0.d0; g = 0.d0; d = 0.d0
        cart          = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                cart(iatom,idir) = x(iatom,idir)
            enddo
        enddo
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = pengygs*legacy_hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            write (*, *) 'Only energy is available'
            !       do istate=1,nstates
            !         p(istate)=PENGYGS*legacy_hartree_to_electronvolt
            !       enddo
            !       do iatom=1,natoms
            !       do idir=1,3
            !         g(1,iatom,idir)=DGSCART(iatom,idir)*51.422067
            !       enddo
            !       enddo
        elseif( igrad == 2 )then
            !       write (*,*) 'Only energy and gradient are available'
            write (*, *) 'Only energy is available'
        endif

        contains

        !> @brief Evaluate POT for the C2H6OH_VBMM potential.
        subroutine pot
            implicit none
            double precision :: coord(30), dx(30), en, vip, vop, vstr, vtor
            integer          :: i, n
            !      DIMENSION COORD(27),DX(27)
            !      DIMENSION dxvstr(27),dxvtor(27),dxvop(27),dxvip(27)
            !C      include 'c2h6cl.inc'
            !     write(97,*) 'pasa'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = int(isurf*(isurf+1)/2)
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            ! Start of changes by EMN.
            ! The lines below do not seem to be needed because the Cartesian
            ! coordinates are in common qpdot.
            ! Besides, BOHR_A in utility.f does not match the value below.
            ! As a consequence, q(i) changes and energy is not conserved
            !c call de venus
            call cartou
            call carttor
            !     PUT COORDINATES IN PROPER ARRAYS
            ! Inicialization de variables
            do i=1, 27
                q(i) = r(i)
                !          q(I) = R(I)*0.52918d0
                pdot(i) = 0.d0
            enddo
            ! Changing units from angstrom to bohr
            !      DO  I = 1, 27
            !c         q(I) = COORD(I)*0.52918d0
            !        q(I) = R(I)*0.52918d0
            !      enddo
            ! End of corrections by EMN
            en = 0.0d0
            call coorden
            !  calculate switching functions
            call switchf
            !  calculate reference angles and their derivatives
            call refangles
            ! CRD 2011 incluimos los terminos de potencial de torsion
            call torsion(vtor)
            !     write(*,*) 'vtor', vtor*0.03812D0
            !  calculate stretching potential
            call stretch(vstr)
            !     write(*,*) 'vstr', vstr*0.03812D0
            !  calculate out of plane bending potential
            call opbend(vop)
            !     write(*,*) 'vop', vop*0.03812D0
            !  calculate in plane bending potential
            call ipbend(vip)
            !     write(*,*) 'vip', vip*0.03812D0
            en = vstr+vop+vip+vtor
            !  convert from 10(5) j/mol to au
            en = en*0.03812d0
            !      V = en
            engygs = en
            call eunitzero
            if( nder .ne. 0 )then
                !  Initialize derivatives
                do n=1, 30
                    !          DX(n)=0.0d0
                    degsdr(i) = pdot(i)*0.0201723d0
                enddo
                ! Start of changes by EMN.
                ! The lines below are not needed because the derivatives are calculated
                ! in venuspotlib_3.f. Anyway NDER=1
                !      PASO= 1.0D-8
                !     PASO= 1.0D-6
                !      do I=1,27
                !        q(I)=q(I) + PASO
                !        call coorden
                !        call switchf
                !        call refangles
                !        call torsion(vtor)
                !        call stretch(vstr)
                !        call opbend(vop)
                !        call ipbend(vip)
                !        en=vstr+vop+vip+vtor
                !C
                !c  convert from 10(5) j/mol to au
                !c
                !        en=en*0.03812D0
                !c
                !C       DX(I)=(en-V)/PASO
                !C       DX(I)=DX(I)*0.52918d0
                !C       pdot(I)=DX(I)
                !        DEGSDR(I)=(en-ENGYGS)/PASO
                !        DEGSDR(I)=DEGSDR(I)*0.52918d0
                !        q(I)=q(I)-PASO
                !      enddo
                !c
                ! 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
                ! pdot sale de la derivada en hartrees/bohr por eso no lo multiplico
                !        do i=1,27
                !           DEGSDR(i)=pdot(i)*0.0201723d0
                !           DEGSDR(i)=pdot(i)
                !        enddo
                ! End of corrections by EMN
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate coorden for the C2H6OH_VBMM potential.
        subroutine coorden
            !  calculates relative coordinates and bond lengths
            implicit none
            double precision :: argmax, argp3, argp4, p3, p4, t3tmp, t4tmp
            integer          :: i, ind
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  calculate relative coordinates
            do ind=1, 3
                tc1b(ind) = q(nc1(ind))-q(nhb(ind))
                tc2b(ind) = q(nc2(ind))-q(nhb(ind))
                tcc(ind)  = q(nc1(ind))-q(nc2(ind))
                ! CRD 2019
                tno(ind) = q(no(ind))-q(nhb(ind))
            enddo
            do ind=1, 3
                do i=1, 3
                    tc1h(i,ind)  = q(nc1(ind))-q(nh(i,ind))
                    tc2h(i,ind)  = q(nc2(ind))-q(nh(i+3,ind))
                    tbh(i,ind)   = q(nhb(ind))-q(nh(i,ind))
                    tbh(i+3,ind) = q(nhb(ind))-q(nh(i+3,ind))
                enddo
            enddo
            ! coordenadas entre los hidrogenos de los dos carbonos
            !      do i=1,3
            !         do j=4,6
            !            do ind=1,3
            !          thh(i,j,ind)=q(nh(j,ind))-q(nh(i,ind))
            !            enddo
            !         enddo
            !      enddo
            !  calculate bond lengths
            rc1b = sqrt(tc1b(1)*tc1b(1)+tc1b(2)*tc1b(2)+tc1b(3)*tc1b(3))
            rc2b = sqrt(tc2b(1)*tc2b(1)+tc2b(2)*tc2b(2)+tc2b(3)*tc2b(3))
            rcc  = sqrt(tcc(1)*tcc(1)+tcc(2)*tcc(2)+tcc(3)*tcc(3))
            ! CRD 2019
            rno = sqrt(tno(1)*tno(1)+tno(2)*tno(2)+tno(3)*tno(3))
            ! fin CRD 2019
            !      write(*,*) rc1b, rc2b
            do i=1, 3
                rc1h(i) = sqrt(tc1h(i,1)*tc1h(i,1)+tc1h(i,2)*tc1h(i,2)+ tc1h(i,3)*tc1h(i,3))
                ! redondeo para corregir las frecuencias
                !       rc1h(i)=(anint(rc1h(i)*100000))/100000
                rc2h(i) = sqrt(tc2h(i,1)*tc2h(i,1)+tc2h(i,2)*tc2h(i,2)+ tc2h(i,3)*tc2h(i,3))
                ! redondeo para corregir las frecuencias
                !       rc2h(i)=(anint(rc2h(i)*100000))/100000
            enddo
            !      write(*,*) 'rc1h(1)',rc1h(1),'rc1h(2)',rc1h(2),'rc1h(3)',rc1h(3)
            !      write(*,*) 'rc2h(1)',rc2h(1),'rc2h(2)',rc2h(2),'rc2h(3)',rc2h(3)
            do i=1, 6
                rbh(i) = sqrt(tbh(i,1)*tbh(i,1)+tbh(i,2)*tbh(i,2)+ tbh(i,3)*tbh(i,3))
            enddo
            ! calculamos las distancias H-H
            !      do i=1,3
            !         do j=4,6
            !        rhh(i,j)=sqrt(thh(i,j,1)*thh(i,j,1)+thh(i,j,2)*thh(i,j,2)+
            !    *                thh(i,j,3)*thh(i,j,3))
            !      write(55,*) 'rhh', i, j, '=', rhh(i,j)
            !         enddo
            !      enddo
            ! crd 2013 modificaciones de las distancias
            argmax = 19.d0
            p3     = 1.d0
            do i=1, 3
                argp3 = (w3*(rc1h(i)-w4))
                if( argp3 .lt. argmax )then
                    t3tmp = 1-tanh(argp3)
                else
                    t3tmp = 0.d0
                endif
                p3 = p3*t3tmp
            enddo
            r0c1h = p3*r0c1hr+(1-p3)*r0c1hp
            !      r0c1h=r0c1hr
            p4 = 1.d0
            do i=1, 3
                argp4 = (w5*(rc2h(i)-w6))
                if( argp4 .lt. argmax )then
                    t4tmp = 1-tanh(argp4)
                else
                    t4tmp = 0.d0
                endif
                p4 = p4*t4tmp
            enddo
            r0c2h = p4*r0c2hr+(1-p4)*r0c2hp
            !      r0c2h=r0c2hr
            !      write(*,*) 'r0c1h', r0c1h, 'r0c2h',r0c2h
        end subroutine coorden

        !> @brief Evaluate torsion for the C2H6OH_VBMM potential.
        subroutine torsion( vtor )
            implicit none
            double precision :: dihed, t1(3), t2(3), vtor
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !      double precision pi
            !      pi=4.0d0*atan(1.0d0)
            ! CRD 2011 incluimos un apartado para calcular los angulos dihedro
            ! para poder calcular despues vtorsion
            call dihedro(3, 1, 2, 6, dihed)
            thetahcch(1) = dihed*pi/180.0d0
            call dihedro(3, 1, 2, 7, dihed)
            thetahcch(2) = dihed*pi/180.0d0
            call dihedro(3, 1, 2, 8, dihed)
            thetahcch(3) = dihed*pi/180.0d0
            call dihedro(4, 1, 2, 6, dihed)
            thetahcch(4) = dihed*pi/180.0d0
            call dihedro(4, 1, 2, 7, dihed)
            thetahcch(5) = dihed*pi/180.0d0
            call dihedro(4, 1, 2, 8, dihed)
            thetahcch(6) = dihed*pi/180.0d0
            call dihedro(5, 1, 2, 6, dihed)
            thetahcch(7) = dihed*pi/180.0d0
            call dihedro(5, 1, 2, 7, dihed)
            thetahcch(8) = dihed*pi/180.0d0
            call dihedro(5, 1, 2, 8, dihed)
            thetahcch(9) = dihed*pi/180.0d0
            ! CRD 2011 calculamos los terminos switching t1 y t2
            do i=1, 3
                t1(i) = 0.5d0*(1-tanh(w1*(rc1h(i)-w2)))
                t2(i) = 0.5d0*(1-tanh(w1*(rc2h(i)-w2)))
                !       write(*,*) 't1',i,t1(i),'t2',i,t2(i)
                !       t1(i)=1-tanh(w1*(rc1h(i)-w2))
                !       t2(i)=1-tanh(w1*(rc2h(i)-w2))
            enddo
            ! CRD 2011 inicializamos vtor
            vtor = 0.0d0
            vtor = vtor+(v3/3.0d0)*(1+cos(3*thetahcch(1)))*t1(1)*t2(1)
            vtor = vtor+(v3/3.0d0)*(1+cos(3*thetahcch(2)))*t1(1)*t2(2)
            vtor = vtor+(v3/3.0d0)*(1+cos(3*thetahcch(3)))*t1(1)*t2(3)
            vtor = vtor+(v3/3.0d0)*(1+cos(3*thetahcch(4)))*t1(2)*t2(1)
            vtor = vtor+(v3/3.0d0)*(1+cos(3*thetahcch(5)))*t1(2)*t2(2)
            vtor = vtor+(v3/3.0d0)*(1+cos(3*thetahcch(6)))*t1(2)*t2(3)
            vtor = vtor+(v3/3.0d0)*(1+cos(3*thetahcch(7)))*t1(3)*t2(1)
            vtor = vtor+(v3/3.0d0)*(1+cos(3*thetahcch(8)))*t1(3)*t2(2)
            vtor = vtor+(v3/3.0d0)*(1+cos(3*thetahcch(9)))*t1(3)*t2(3)
        end subroutine torsion

        !> @brief Evaluate dihedro for the C2H6OH_VBMM potential.
        subroutine dihedro( i, j, k, l, dihed )
            implicit none
            double precision :: dihed
            !        include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            integer :: i, j, k, l
            integer :: m
            !        double precision x(300)
            !        double precision dista
            !        double precision d12,d13,d14,d23,d24,d34
            !        double precision p,q
            !        double precision p
            !        double precision pi
            double precision :: rij(3), rjk(3), rkj(3), rkl(3)
            double precision :: zijjk(3), zijkj(3), zkjkl(3), zjkkl(3)
            double precision :: tmpsign(3), tmp, arg1, arg2
            !        pi=4.0d0*atan(1.0d0)
            do m = 1, 3
                rij(m) = q((i-1)*3+m)-q((j-1)*3+m)
                rjk(m) = q((j-1)*3+m)-q((k-1)*3+m)
                rkl(m) = q((k-1)*3+m)-q((l-1)*3+m)
            enddo
            call crosprod(rjk, rkl, zjkkl)
            call crosprod(rij, rjk, zijjk)
            arg2 = anorm(rjk)*dotprod(rij, zjkkl)
            arg1 = dotprod(zijjk, zjkkl)
            if( abs(arg1) .lt. 1.d-6 )then
                if( arg2 .gt. 0.d0 ) dihed = 90.d0
                if( arg2 .lt. 0.d0 ) dihed = -90.d0
                return
            endif
            if( arg1 .gt. 0.d0 )then
                dihed = atan(arg2/arg1) *180.d0/pi
                return
            endif
            if( arg2 .ge. 0.d0 )then
                dihed = 180.d0 + atan(arg2/arg1) *180.d0/pi
                return
            endif
            if( arg2 .lt. 0.d0 )then
                dihed = -180.d0 + atan(arg2/arg1) *180.d0/pi
                return
            endif
        end subroutine dihedro

        !> @brief Evaluate crosprod for the C2H6OH_VBMM potential.
        subroutine crosprod( x, y, z )
            implicit none
            double precision :: x(3), y(3), z(3)
            z(1) = x(2)*y(3)-x(3)*y(2)
            z(2) = x(3)*y(1)-x(1)*y(3)
            z(3) = x(1)*y(2)-x(2)*y(1)
        end subroutine crosprod

        !> @brief Evaluate dotprod for the C2H6OH_VBMM potential.
        double precision function dotprod( x, y )
            implicit none
            integer          :: k
            double precision :: x(3), y(3), tmp
            tmp = 0.d0
            do k = 1, 3
                tmp = tmp + x(k)*y(k)
            enddo
            dotprod = tmp
        end function dotprod

        !> @brief Evaluate anorm for the C2H6OH_VBMM potential.
        double precision function anorm( x )
            implicit none
            integer          :: k
            double precision :: x(3), tmp
            tmp = 0.d0
            do k = 1, 3
                tmp = tmp + x(k)*x(k)
            enddo
            anorm = dsqrt(tmp)
        end function anorm

        !> @brief Evaluate refangles for the C2H6OH_VBMM potential.
        subroutine refangles
            !  subroutine calculates reference angles for the "in-plane" potential
            implicit none
            double precision :: ddr(4,4), halfpi, sumd2(4), sumd4(4), tau, taucc, twopi
            integer          :: i, j
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            tau = acos(-1.0d0/3.0d0)
            ! CRD 2012 taucc para el angulo H-C-C
            taucc = acos(-0.3616d0)
            !      write(*,*) 'tau', tau, 'acostau', tau*57.2957795
            !      write(*,*) 'taucc', taucc, 'acostaucc', taucc*57.2957795
            !      pi=4.0d0*atan(1.0d0)
            halfpi = 0.5d0*pi
            twopi  = 2.0d0*pi
            !      write(38,*) tau, pi, halfpi, twopi
            ! jcc-2010
            ! valor para ajustar el angulo
            !      anghch=0.649d0*pi
            !      tausih = pi - asin ( sin(anghch/2.d0) / sin(pi/3.d0) )
            !      tausih=acos(-0.4648d0)
            ! jcc-2010
            !  set diagonal elements to zero
            do i=1, 8
                theta0(i,i) = 0.0d0
            enddo
            !  calculate reference angles
            !      theta0(1,2)=tau+(tau-tausih)*(sphi(1)*sphi(2)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
            !      theta0(1,3)=tau+(tau-tausih)*(sphi(1)*sphi(3)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
            !      theta0(1,4)=taucc+(tau-tausih)*(sphi(1)*sphi(4)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
            !      theta0(2,3)=tau+(tau-tausih)*(sphi(2)*sphi(3)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
            !      theta0(2,4)=taucc+(tau-tausih)*(sphi(2)*sphi(4)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
            !      theta0(3,4)=taucc+(tau-tausih)*(sphi(3)*sphi(4)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
            !      theta0(5,6)=tau+(tau-tausih)*(sphi(5)*sphi(6)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(7)*stheta(8)-1.0d0)
            !      theta0(5,7)=tau+(tau-tausih)*(sphi(5)*sphi(7)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(6)*stheta(8)-1.0d0)
            !      theta0(5,8)=taucc+(tau-tausih)*(sphi(5)*sphi(8)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(6)*stheta(7)-1.0d0)
            !      theta0(6,7)=tau+(tau-tausih)*(sphi(6)*sphi(7)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(5)*stheta(8)-1.0d0)
            !      theta0(6,8)=taucc+(tau-tausih)*(sphi(6)*sphi(8)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(5)*stheta(7)-1.0d0)
            !      theta0(7,8)=taucc+(tau-tausih)*(sphi(7)*sphi(8)-1.0d0)
            !    *             +(tau-anghch/3.0d0)*(stheta(5)*stheta(6)-1.0d0)
            ! el atomo 4 es el c2
            theta0(1,2) = tau+(tau-halfpi)*(sphi(1)*sphi(2)-1.0d0) +(tau-twopi/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
            theta0(1,3) = tau+(tau-halfpi)*(sphi(1)*sphi(3)-1.0d0) +(tau-twopi/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
            theta0(1,4) = taucc+(tau-halfpi)*(sphi(1)*sphi(4)-1.0d0) +(tau-twopi/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
            theta0(2,3) = tau+(tau-halfpi)*(sphi(2)*sphi(3)-1.0d0) +(tau-twopi/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
            theta0(2,4) = taucc+(tau-halfpi)*(sphi(2)*sphi(4)-1.0d0) +(tau-twopi/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
            theta0(3,4) = taucc+(tau-halfpi)*(sphi(3)*sphi(4)-1.0d0) +(tau-twopi/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
            !  el atomo 8 es el c1
            theta0(5,6) = tau+(tau-halfpi)*(sphi(5)*sphi(6)-1.0d0) +(tau-twopi/3.0d0)*(stheta(7)*stheta(8)-1.0d0)
            theta0(5,7) = tau+(tau-halfpi)*(sphi(5)*sphi(7)-1.0d0) +(tau-twopi/3.0d0)*(stheta(6)*stheta(8)-1.0d0)
            theta0(5,8) = taucc+(tau-halfpi)*(sphi(5)*sphi(8)-1.0d0) +(tau-twopi/3.0d0)*(stheta(6)*stheta(7)-1.0d0)
            theta0(6,7) = tau+(tau-halfpi)*(sphi(6)*sphi(7)-1.0d0) +(tau-twopi/3.0d0)*(stheta(5)*stheta(8)-1.0d0)
            theta0(6,8) = taucc+(tau-halfpi)*(sphi(6)*sphi(8)-1.0d0) +(tau-twopi/3.0d0)*(stheta(5)*stheta(7)-1.0d0)
            theta0(7,8) = taucc+(tau-halfpi)*(sphi(7)*sphi(8)-1.0d0) +(tau-twopi/3.0d0)*(stheta(5)*stheta(6)-1.0d0)
            !  fill in the other half of the matrix
            do i=1, 3
                do j=i+1, 4
                    theta0(j,i) = theta0(i,j)
                enddo
            enddo
            do i=5, 7
                do j=i+1, 8
                    theta0(j,i) = theta0(i,j)
                enddo
            enddo
            !      do i=1,3
            !         do j=i+1,4
            !          write(*,*) 'theta',i,j,'=', theta0(i,j)*180/pi
            !         enddo
            !      enddo
            !      do i=5,7
            !         do j=i+1,8
            !          write(*,*) 'theta',i,j,'=', theta0(i,j)*180/pi
            !         enddo
            !      enddo
        end subroutine refangles

        !> @brief Evaluate stretch for the C2H6OH_VBMM potential.
        subroutine stretch( vstr )
            !  subroutine to calculate leps-type stretching potential and its
            !  derivatives
            implicit none
            double precision :: ach1, ach2, achdc(3), achdh(4,3), achdx(3), arga1, arga2, dt, e1, e3, expterm, factj(4), rav1, &
                & rav2, vcc, vj(3), vjbh(6), vjc1b, vjc1h(3), vjc2b, vjc2h(3), vno, vq(3), vqbh(6), vqc1b, vqc1h(3), vqc2b, &
                & vqc2h(3), vstr
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !      write(*,*) nnc1,nnc2,nnb
            !      write(*,*) r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h
            !      write(*,*) a1c1h,b1c1h,c1c1h
            !      write(*,*) r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h
            !      write(*,*) a1c2h,b1c2h,c1c2h
            !      write(*,*) r0cc,rcc,d1cc,acc
            !      write(*,*) a1c,a2c
            !      write(*,*) r0hh,d1hh,d3hh,ahh
            !      write(*,*) r0cb,d1cb,d3cb,acb
            !      write(*,*) a3s,b3s,aphi,bphi,cphi
            !      write(*,*) atheta,btheta,ctheta
            !      write(*,*) fch3,hch3
            !      write(*,*) fkinf,ak,bk,aa1,aa2,aa3,aa4
            !      write(*,*) V3,w1,w2
            !      write(*,*) a1s,b1s,a2s,b2s
            !      write(*,*) a1scc,b1scc,a2scc,b2scc,a3scc,b3scc
            !  calculate avergage bond length for the methane moiety
            !C       rav=(rc1h(1)+rc1h(2)+rc1h(3)+rc2h(1)+rc2h(2)+rc2h(3))/6.0d0
            rav1 = (rc1h(1)+rc1h(2)+rc1h(3))/3.0d0
            rav2 = (rc2h(1)+rc2h(2)+rc2h(3))/3.0d0
            !      write(*,*) 'rav1',rav1,'rav2', rav2
            !  initialise:
            vstr = 0.0d0
            !  ach:
            !  nb: in double precision tanh(19.0d0)=1.0d0 and we put the if statemen
            !  in to avoid overflow/underflow errors
            !     write(*,*) 'r0c1h',r0c1h,'r0c2h',r0c2h
            arga1 = c1c1h*(rav1-r0c1h)
            arga2 = c1c2h*(rav2-r0c2h)
            !      write(*,*) 'arga1', arga1,'arga2',arga2
            if( arga1 .lt. 19.0d0 )then
                ach1 = a1c1h+b1c1h*(tanh(arga1)+1.0d0)*0.5d0
            else
                ach1 = a1c1h+b1c1h
            endif
            if( arga2 .lt. 19.0d0 )then
                ach2 = a1c2h+b1c2h*(tanh(arga2)+1.0d0)*0.5d0
            else
                ach2 = a1c2h+b1c2h
            endif
            !      write(6,*) r0c1h, rav1, r0c2h, rav2, ach1, ach2
            !      argac=c1cc*(rcc-r0cc)
            !      if(argac.lt.19.0d0)then
            !        acc=a1cc+b1cc*(tanh(argac)+1.0d0)*0.5d0
            !      else
            !        acc=a1cc+b1cc
            !      endif
            !  calculate singlet: e1, triplet: e3 energies and vq and vj
            !  terms for each bond
            e1    = d1cb*(exp(-2.0d0*acb*(rc1b-r0cb))-2.0d0*exp(-acb*(rc1b-r0cb)))
            e3    = d3cb*(exp(-2.0d0*acb*(rc1b-r0cb))+2.0d0*exp(-acb*(rc1b-r0cb)))
            vqc1b = (e1+e3)*0.5d0
            vjc1b = (e1-e3)*0.5d0
            e1    = d1cb*(exp(-2.0d0*acb*(rc2b-r0cb))-2.0d0*exp(-acb*(rc2b-r0cb)))
            e3    = d3cb*(exp(-2.0d0*acb*(rc2b-r0cb))+2.0d0*exp(-acb*(rc2b-r0cb)))
            vqc2b = (e1+e3)*0.5d0
            vjc2b = (e1-e3)*0.5d0
            !     write(*,*) 'vqc1b',vqc1b,'vqc2b',vqc2b,'vjc1b',vjc1b,'vjc2b',vjc2b
            !      e1=d1cc*(exp(-2.0d0*acc*(rcc-r0cc))-2.0d0*exp(-acc*(rcc-r0cc)))
            !      e3=d3cc*(exp(-2.0d0*acc*(rcc-r0cc))+2.0d0*exp(-acc*(rcc-r0cc)))
            !      vqcc=(e1+e3)*0.5d0
            !      vjcc=(e1-e3)*0.5d0
            do i=1, 3
                e1       = d1c1h*(exp(-2.0d0*ach1*(rc1h(i)-r0c1h)) -2.0d0*exp(-ach1*(rc1h(i)-r0c1h)))
                e3       = d3c1h*(exp(-2.0d0*ach1*(rc1h(i)-r0c1h)) +2.0d0*exp(-ach1*(rc1h(i)-r0c1h)))
                vqc1h(i) = (e1+e3)*0.5d0
                vjc1h(i) = (e1-e3)*0.5d0
                e1       = d1c2h*(exp(-2.0d0*ach2*(rc2h(i)-r0c2h)) -2.0d0*exp(-ach2*(rc2h(i)-r0c2h)))
                e3       = d3c2h*(exp(-2.0d0*ach2*(rc2h(i)-r0c2h)) +2.0d0*exp(-ach2*(rc2h(i)-r0c2h)))
                vqc2h(i) = (e1+e3)*0.5d0
                vjc2h(i) = (e1-e3)*0.5d0
            enddo
            do i=1, 6
                e1      = d1hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh)) -2.0d0*exp(-ahh*(rbh(i)-r0hh)))
                e3      = d3hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh)) +2.0d0*exp(-ahh*(rbh(i)-r0hh)))
                vqbh(i) = (e1+e3)*0.5d0
                vjbh(i) = (e1-e3)*0.5d0
            enddo
            !  calculate 3 body potential
            do i=1, 3
                vq(i) = vqc1h(i)+vqc1b+vqbh(i)
                vj(i) = -sqrt(((vjc1h(i)-vjc1b)**2+(vjc1b-vjbh(i))**2 +(vjbh(i)-vjc1h(i))**2)*0.5d0)
                vstr  = vstr+vq(i)+vj(i)
            enddo
            do i=1, 3
                vq(i) = vqc2h(i)+vqc2b+vqbh(i+3)
                vj(i) = -sqrt(((vjc2h(i)-vjc2b)**2+(vjc2b-vjbh(i+3))**2 +(vjbh(i+3)-vjc2h(i))**2)*0.5d0)
                vstr  = vstr+vq(i)+vj(i)
            enddo
            ! CRD 2012 potencial morse para el enlace C-C
            dt      = (rcc-r0cc)
            expterm = exp(-acc*dt)
            vcc     = d1cc*(1.0d0-expterm)**2.0d0
            vstr    = vstr + vcc
            ! CRD 2019
            ! O-H  new term (simple morse term)
            dt      = (rno-r0hh)
            expterm = exp(-ahh*dt)
            vno     = d1hh*(1.d0-expterm)**2.d0
            vstr    = vstr + vno
            ! fin CRD 2019
        end subroutine stretch

        !> @brief Evaluate opbend for the C2H6OH_VBMM potential.
        subroutine opbend( vop )
            !  subroutine calculates symmetrized vop potential and derivatives
            implicit none
            double precision :: a(3), argd(8), axb(3), b(3), c(8,3), sum2, sum4, sumd2(8), sumd4(8), vop
            integer          :: i, in(3), j, k, l
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            double precision :: norma
            vop = 0.0d0
            !  calculate force constants and their derivatives
            call opforce
            !  calculate out-of-plane angle and derivatives
            ! CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C2
            !      do i=1,4
            i = 4
            !      do i=1,3
            j = i+1
            if( j .gt. 4 ) j=j-4
            k = j+1
            if( k .gt. 4 ) k=k-4
            l = k+1
            if( l .gt. 4 ) l=l-4
            call calcdelta1(i, j, k, l, sum2, sum4)
            sumd2(i) = sum2
            sumd4(i) = sum4
            !      write(*,*) 'sum2',sum2,'sum4',sum4
            !      write(*,*) 'sumd2',sumd2(i),'sumd4',sumd4(i)
            vop = vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
            !      write(*,*) i,j,k,l,sumd2(i),fdelta(i),sumd4(i),hdelta(i)
            !      enddo
            ! CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C1
            !      do i=5,8
            i = 8
            !      do i=4,6
            j = i+1
            if( j .gt. 8 ) j=j-4
            k = j+1
            if( k .gt. 8 ) k=k-4
            l = k+1
            if( l .gt. 8 ) l=l-4
            call calcdelta2(i, j, k, l, sum2, sum4)
            sumd2(i) = sum2
            sumd4(i) = sum4
            !      write(*,*) 'sum2',sum2,'sum4',sum4
            !      write(*,*) 'sumd2',sumd2(i),'sumd4',sumd4(i)
            vop = vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
            !      write(*,*) i,j,k,l,sumd2(i),fdelta(i),sumd4(i),hdelta(i)
            !      enddo
        end subroutine opbend

        !> @brief Evaluate ipbend for the C2H6OH_VBMM potential.
        subroutine ipbend( vip )
            !  subroutine calculates symmetrised in plane bend term
            !  and its derivatives
            implicit none
            double precision :: arga, cosine, costh(8,8), dang, dot, dth(8,8), theta(8,8), vip
            integer          :: i, j
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  initialise
            vip = 0.0d0
            !  calculate force constants: fk0(i,j), f1(i)
            !  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
            call ipforce
            !      do i=1,8
            !      write(*,*) f1(i)
            !      enddo
            !      write(*,*) 'fin f1(i)'
            !  calculate theta(i,j) and in plane bend potential
            ! bucle para el c1
            do i=1, 3
                do j=i+1, 4
                    if( j .lt. 4.0 )then
                        costh(i,j) = tc1h(i,1)*tc1h(j,1)+tc1h(i,2)*tc1h(j,2) +tc1h(i,3)*tc1h(j,3)
                        costh(i,j) = costh(i,j)/rc1h(i)/rc1h(j)
                        theta(i,j) = acos(costh(i,j))
                        dth(i,j)   = theta(i,j)-theta0(i,j)
                        vip        = vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
                    else
                        costh(i,j) = tc1h(i,1)*tcc(1)+tc1h(i,2)*tcc(2) +tc1h(i,3)*tcc(3)
                        !      write(92,*) 'rc1h',i,'=', rc1h(i), 'rcc=', rcc
                        costh(i,j) = costh(i,j)/rc1h(i)/rcc
                        theta(i,j) = acos(costh(i,j))
                        dth(i,j)   = theta(i,j)-theta0(i,j)
                        vip        = vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
                    endif
                    !      write(*,*) 'fk0',i,j,'=',fk0(i,j)
                    !      write(*,*) 'f1',i,'=',f1(i)
                    !      write(*,*) 'f1',j,'=',f1(j)
                    !      write(*,*) 'rcc','=',rcc
                    !      write(*,*) 'costh',i,j,'=',costh(i,j)
                    !      write(*,*) 'theta',i,j,'=',theta(i,j)
                    !      write(*,*) 'theta0',i,j,'=',theta0(i,j)
                    !      write(*,*) 'dth',i,j,'=',dth(i,j)
                    !      write(*,*) 'costh',i,j,'=',costh(i,j)
                enddo
            enddo
            ! bucle para el c2 para los theta 5,6,7,8
            do i=1, 3
                do j=i+1, 4
                    if( j .lt. 4.0 )then
                        costh(i+4,j+4) = tc2h(i,1)*tc2h(j,1)+tc2h(i,2)*tc2h(j,2) +tc2h(i,3)*tc2h(j,3)
                        costh(i+4,j+4) = costh(i+4,j+4)/rc2h(i)/rc2h(j)
                        theta(i+4,j+4) = acos(costh(i+4,j+4))
                        dth(i+4,j+4)   = theta(i+4,j+4)-theta0(i+4,j+4)
                        vip            = vip+0.5d0*fk0(i+4,j+4)*f1(i+4)*f1(j+4)*dth(i+4,j+4)**2
                    else
                        ! modificacion para corregir los angulos del c2
                        !          costh(i,j)=tc2h(i,1)*tcc(1)+tc2h(i,2)*tcc(2)
                        !    *                       +tc2h(i,3)*tcc(3)
                        costh(i+4,j+4) = tc2h(i,1)*(-tcc(1))+tc2h(i,2)*(-tcc(2)) +tc2h(i,3)*(-tcc(3))
                        costh(i+4,j+4) = costh(i+4,j+4)/rc2h(i)/rcc
                        theta(i+4,j+4) = acos(costh(i+4,j+4))
                        dth(i+4,j+4)   = theta(i+4,j+4)-theta0(i+4,j+4)
                        vip            = vip+0.5d0*fk0(i+4,j+4)*f1(i+4)*f1(j+4)*dth(i+4,j+4)**2
                    endif
                enddo
            enddo
            ! CRD 2019
            !  Now, calculate h2o bending energies.
            !  First, calculate the angles:
            do i = 1, 6
                dot = 0.d0
                do j = 1, 3
                    dot = dot - tno(j)*tbh(i,j)
                enddo
                cosine    = dot / (rno*rbh(i))
                cosine    = min(1.d0, max(-1.d0, cosine))
                angh2o(i) = acos(cosine)
            enddo
            !      write (*,*) "angulos:"
            !      do i=1,4
            !      write (*,*) i, angh2o(i)*360.d0/(2.d0*3.1415926d0)
            !      enddo
            ! Now, calculate each force constant fkh2o(i) as function of the O-H(I)
            !  distance, rbh(i)
            do i=1, 6
                arga = alph2o * (rbh(i) - r0hh)
                if( arga .lt. 19.0d0 )then
                    fkh2o(i) = fkh2oeq * (1 - tanh(arga))
                else
                    fkh2o(i) = 0.0d0
                endif
            enddo
            !  Now calculate the contribution to the energy for each H(I)-O-H(O)
            !  harmonic bending term
            do i=1, 6
                dang = (angh2o(i) - angh2oeq)
                vip  = vip+0.5d0* fkh2o(i) * dang * dang
            enddo
            ! fin CRD 2019
        end subroutine ipbend

        !> @brief Evaluate calcdelta1 for the C2H6OH_VBMM potential.
        subroutine calcdelta1( i, j, k, l, sum2, sum4 )
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            implicit none
            double precision :: a(3), argd(8), atemp2(3), axb(3), b(3), c(8,3), cdot(8,3,3), daxb(8,3,3), delta(8), sum2, sum4
            integer          :: i, ii, in(3), ind, j, k, l
            double precision :: norma
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  initialise
            sum2 = 0.0d0; sum4 = 0.0d0
            !  set j,k,l indices
            in(1) = j
            in(2) = k
            in(3) = l
            ! CRD 2005 incluimos bucles if para asegurar que siempre
            ! que una de los valores j, k o l sea 4 las coordenadas
            ! corresponden al atomo nuevo.
            !  vector a is rk-rj, vector b is rl-rj
            if( j .eq. 4 )then
                do ind=1, 3
                    a(ind) = q(nh(k,ind))-q(nc2(ind))
                    b(ind) = q(nh(l,ind))-q(nc2(ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            if( k .eq. 4 )then
                do ind=1, 3
                    a(ind) = q(nc2(ind))-q(nh(j,ind))
                    b(ind) = q(nh(l,ind))-q(nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            if( l .eq. 4 )then
                do ind=1, 3
                    a(ind) = q(nh(k,ind))-q(nh(j,ind))
                    b(ind) = q(nc2(ind))-q(nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                ! modificacion para que hiciese vxu
                !      axb(1)=b(2)*a(3)-b(3)*a(2)
                !      axb(2)=a(1)*b(3)-b(1)*a(3)
                !      axb(3)=b(1)*a(2)-a(1)*b(2)
                norma = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma = sqrt(norma)
            else
            endif
            ! CRD 2012 eliminamos porque seria la salida del C2
            if( k .lt. 4.and.j .lt. 4.and.l .lt. 4 )then
                do ind=1, 3
                    a(ind) = q(nh(k,ind))-q(nh(j,ind))
                    b(ind) = q(nh(l,ind))-q(nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
                !      write(*,*) 'norma=',norma
            else
            endif
            !  c is position vector of h(ii): calculate c(j),c(k),c(l)
            do ii=1, 3
                if( in(ii) .eq. 4 )then
                    do ind=1, 3
                        c(in(ii),ind) = -tcc(ind)/rcc
                    enddo
                else
                    do ind=1, 3
                        c(in(ii),ind) = -tc1h(in(ii),ind)/rc1h(in(ii))
                        !      write(*,*) 'calcdelta1 in(ii)', in(ii),c(in(ii),ind)
                    enddo
                endif
            enddo
            !  argd is the dot product axb dot c
            do ii=1, 3
                argd(in(ii)) = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                argd(in(ii)) = argd(in(ii))/norma
                ! Start of corrections by EMN
                ! The two lines below are added to avoid problems when acos is
                ! calculated
                if( argd(in(ii)) > 1.d0 ) argd(in(ii))  = 1
                if( argd(in(ii)) < -1.d0 ) argd(in(ii)) = -1
                ! End of corrections by EMN
                delta(in(ii)) = acos(argd(in(ii)))-theta0(i,in(ii))
                sum2          = sum2+delta(in(ii))**2
                sum4          = sum4+delta(in(ii))**4
                !      write(*,*) 'calcdelta1 sum2', sum2
                !      write(*,*) 'calcdelta1 sum4', sum4
            enddo
        end subroutine calcdelta1

        !> @brief Evaluate calcdelta2 for the C2H6OH_VBMM potential.
        subroutine calcdelta2( i, j, k, l, sum2, sum4 )
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            implicit none
            double precision :: a(3), argd(8), atemp2(3), axb(3), b(3), c(8,3), cdot(8,3,3), daxb(8,3,3), delta(8), sum2, sum4
            integer          :: i, ii, in(3), ind, j, k, l
            double precision :: norma
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  initialise
            sum2 = 0.0d0; sum4 = 0.0d0
            !  set j,k,l indices
            !      in(1)=j
            !      in(2)=k
            !      in(3)=l
            in(1) = 4
            in(2) = 5
            in(3) = 6
            ! CRD 2005 incluimos bucles if para asegurar que siempre
            ! que una de los valores j, k o l sea 4 las coordenadas
            ! corresponden al atomo nuevo.
            !  vector a is rk-rj, vector b is rl-rj
            if( j .eq. 8 )then
                do ind=1, 3
                    a(ind) = q(nh(k,ind))-q(nc1(ind))
                    b(ind) = q(nh(l,ind))-q(nc1(ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            if( k .eq. 8 )then
                do ind=1, 3
                    a(ind) = q(nc1(ind))-q(nh(j,ind))
                    b(ind) = q(nh(l,ind))-q(nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            if( l .eq. 8 )then
                do ind=1, 3
                    a(ind) = q(nh(k,ind))-q(nh(j,ind))
                    b(ind) = q(nc1(ind))-q(nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            ! CRD 2012 eliminamos porque seria la salida del C1
            if( k .lt. 8.and.j .lt. 8.and.l .lt. 8 )then
                do ind=1, 3
                    a(ind) = q(nh(k,ind))-q(nh(j,ind))
                    b(ind) = q(nh(l,ind))-q(nh(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
                !      write(*,*) 'calcdelta2 norma',norma
            else
            endif
            !  c is position vector of h(ii): calculate c(j),c(k),c(l)
            do ii=1, 3
                if( in(ii) .eq. 4 )then
                    do ind=1, 3
                        c(in(ii),ind) = -tcc(ind)/rcc
                    enddo
                else
                    do ind=1, 3
                        c(in(ii),ind) = -tc2h(in(ii),ind)/rc2h(in(ii))
                    enddo
                endif
            enddo
            !  argd is the dot product axb dot c
            do ii=1, 3
                argd(in(ii)) = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                argd(in(ii)) = argd(in(ii))/norma
                ! Start of corrections by EMN
                ! The two lines below are added to avoid problems when acos is
                ! calculated
                if( argd(in(ii)) > 1.d0 ) argd(in(ii))  = 1
                if( argd(in(ii)) < -1.d0 ) argd(in(ii)) = -1
                ! End of corrections by EMN
                delta(in(ii)) = acos(argd(in(ii)))-theta0(i,in(ii))
                sum2          = sum2+delta(in(ii))**2
                sum4          = sum4+delta(in(ii))**4
            enddo
        end subroutine calcdelta2

        !> @brief Evaluate opforce for the C2H6OH_VBMM potential.
        subroutine opforce
            !  calculates the out-of-plane bending force constants
            !  and their derivatives
            implicit none
            double precision :: switch(8)
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  calculate switching functions:
            switch(1) = (1.0d0-s3(1))*s3(2)*s3(3)*s3(4)
            switch(2) = (1.0d0-s3(2))*s3(3)*s3(4)*s3(1)
            switch(3) = (1.0d0-s3(3))*s3(4)*s3(1)*s3(2)
            switch(4) = (1.0d0-s3(4))*s3(1)*s3(2)*s3(3)
            switch(5) = (1.0d0-s3(5))*s3(6)*s3(7)*s3(8)
            switch(6) = (1.0d0-s3(6))*s3(7)*s3(8)*s3(5)
            switch(7) = (1.0d0-s3(7))*s3(8)*s3(5)*s3(6)
            switch(8) = (1.0d0-s3(8))*s3(5)*s3(6)*s3(7)
            !      do i=1,8
            !       write(*,*) 'switch',i,'=',switch(i)
            !      enddo
            !       write(*,*) 'pasa'
            !  calculate the force constants
            do i=1, 8
                fdelta(i) = switch(i)*fch3
                hdelta(i) = switch(i)*hch3
            enddo
        end subroutine opforce

        !> @brief Evaluate ipforce for the C2H6OH_VBMM potential.
        subroutine ipforce
            !  calculates the symmetrised in plane bend force constants and
            !  all partial derivatives involving them
            implicit none
            double precision :: a1, a2, arga1, arga2, f0, f2
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  set force constant at asymptotes
            f0       = fkinf+ak
            f2       = fkinf
            fk0(1,2) = f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)*s2(4)-1.0d0)
            fk0(1,3) = f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)*s2(4)-1.0d0)
            fk0(1,4) = f0+f0*(s1(1)*s1(4)-1.0d0)+(f0-f2)*(s2(2)*s2(3)-1.0d0)
            fk0(2,3) = f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)*s2(4)-1.0d0)
            fk0(2,4) = f0+f0*(s1(2)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(3)-1.0d0)
            fk0(3,4) = f0+f0*(s1(3)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(2)-1.0d0)
            fk0(5,6) = f0+f0*(s1(5)*s1(6)-1.0d0)+(f0-f2)*(s2(7)*s2(8)-1.0d0)
            fk0(5,7) = f0+f0*(s1(5)*s1(7)-1.0d0)+(f0-f2)*(s2(6)*s2(8)-1.0d0)
            fk0(5,8) = f0+f0*(s1(5)*s1(8)-1.0d0)+(f0-f2)*(s2(6)*s2(7)-1.0d0)
            fk0(6,7) = f0+f0*(s1(6)*s1(7)-1.0d0)+(f0-f2)*(s2(5)*s2(8)-1.0d0)
            fk0(6,8) = f0+f0*(s1(6)*s1(8)-1.0d0)+(f0-f2)*(s2(5)*s2(7)-1.0d0)
            fk0(7,8) = f0+f0*(s1(7)*s1(8)-1.0d0)+(f0-f2)*(s2(5)*s2(6)-1.0d0)
            ! f1 para los h del c1
            do i=1, 3
                !  calc derivatives of fk0 wrt each of the rch(i) bonds
                !  calculate the terms f1(i)
                arga1 = aa1*rbh(i)*rbh(i)
                arga2 = aa4*(rbh(i)-r0hh)*(rbh(i)-r0hh)
                a1    = 1.0d0-exp(-arga1)
                a2    = aa2+aa3*exp(-arga2)
                f1(i) = a1*exp(-a2*(rc1h(i)-r0c1h)**2)
            enddo
            ! f1 para el c2
            ! CRD 2012 eliminamos el calculo de a1c y a2c pq es constante
            ! lo incluimos en el CONST
            !        arga1c=aa1*rbc*rbc
            !        arga2c=aa4*(rbc-r0cb)*(rbc-r0cb)
            !        a1c=1.0d0-exp(-arga1c)
            !        a2c=aa2+aa3*exp(-arga2c)
            f1(4) = a1c*exp(-a2c*(rcc-r0cc)**2)
            f1(8) = f1(4)
            ! f1 para los h del c2
            do i=5, 7
                !  calc derivatives of fk0 wrt each of the rch(i) bonds
                !  calculate the terms f1(i)
                arga1 = aa1*rbh(i-1)*rbh(i-1)
                arga2 = aa4*(rbh(i-1)-r0hh)*(rbh(i-1)-r0hh)
                a1    = 1.0d0-exp(-arga1)
                a2    = aa2+aa3*exp(-arga2)
                f1(i) = a1*exp(-a2*(rc2h(i-4)-r0c2h)**2)
            enddo
            ! f1 para el c2
            ! CRD 2012 eliminamos el calculo de a1c y a2c pq es constante
            ! lo incluimos en el CONST
            !        arga1c=aa1*rcb*rcb
            !        arga2c=aa4*(rcb-r0cb)*(rcb-r0cb)
            !        a1c=1.0d0-exp(-arga1c)
            !        a2c=aa2+aa3*exp(-arga2c)
            !        f1(8)=a1c*exp(-a2c*(rcc-r0cc)**2)
        end subroutine ipforce

        !> @brief Evaluate switchf for the C2H6OH_VBMM potential.
        subroutine switchf
            !  calculates switching functions: s3,sphi,stheta
            !  and their derivatives ds3,dsphi,dstheta
            implicit none
            double precision :: argmax, args1, args2, args3, argsc1, argsc2, argsc3, argsphi, argsphic, argstheta, argsthetac
            integer          :: i
            !      include 'c2h6cl.inc'
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  nb remember that integration units are:
            !  energy in   1.0d+05 j/mol
            !  time in     1.0d-14 s
            ! original     a1s=1.5313681d-7
            !      a1s=1.5313681d-8
            !      b1s=-4.6696246d0
            ! original      a2s=1.0147402d-7
            !      a2s=1.0147402d-8
            !      b2s=-12.363798d0
            !  use double precision criterion:
            !  tanh(19.0d0)=1.0d0
            argmax = 19.0d0
            !  calculate s1 and ds1
            do i=1, 3
                args1 = a1s*(rc1h(i)-r0c1h)*(rc1h(i)-b1s)**8
                if( args1 .lt. argmax )then
                    s1(i) = 1.0d0-tanh(args1)
                else
                    s1(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                args1 = a1s*(rc2h(i)-r0c2h)*(rc2h(i)-b1s)**8
                if( args1 .lt. argmax )then
                    s1(i+4) = 1.0d0-tanh(args1)
                else
                    s1(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            !      s1(4)=1.0d0
            !      s1(8)=1.0d0
            argsc1 = a1scc*(rcc-r0cc)*(rcc-b1scc)**8
            if( argsc1 .lt. argmax )then
                s1(4) = 1.0d0-tanh(argsc1)
                s1(8) = s1(4)
            else
                s1(4) = 0.0d0
                s1(8) = s1(4)
            endif
            !  calculate s2 and ds2
            do i=1, 3
                args2 = a2s*(rc1h(i)-r0c1h)*(rc1h(i)-b2s)**6
                if( args2 .lt. argmax )then
                    s2(i) = 1.0d0-tanh(args2)
                else
                    s2(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                args2 = a2s*(rc2h(i)-r0c2h)*(rc2h(i)-b2s)**6
                if( args2 .lt. argmax )then
                    s2(i+4) = 1.0d0-tanh(args2)
                else
                    s2(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2005 s2(4)for the new atom
            ! CRD 2012 todos los S de los carbonos son 1
            !      s2(4)=1.0d0
            !      s2(8)=1.0d0
            argsc2 = a2scc*(rcc-r0cc)*(rcc-b2scc)**6
            if( argsc2 .lt. argmax )then
                s2(4) = 1.0d0-tanh(argsc2)
                s2(8) = s2(4)
            else
                s2(4) = 0.0d0
                s2(8) = s2(4)
            endif
            !  calculate s3 and ds3
            do i=1, 3
                args3 = a3s*(rc1h(i)-r0c1h)*(rc1h(i)-b3s)**2
                if( args3 .lt. argmax )then
                    s3(i) = 1.0d0-tanh(args3)
                else
                    s3(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                args3 = a3s*(rc2h(i)-r0c2h)*(rc2h(i)-b3s)**2
                if( args3 .lt. argmax )then
                    s3(i+4) = 1.0d0-tanh(args3)
                else
                    s3(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2005 s3(4) for the new atom
            ! CRD 2012 todos los S de los carbonos son 1
            !      s3(4)=1.0d0
            !      s3(8)=1.0d0
            argsc3 = a3scc*(rcc-r0cc)*(rcc-b3scc)**2
            if( argsc3 .lt. argmax )then
                s3(4) = 1.0d0-tanh(argsc3)
                s3(8) = s3(4)
            else
                s3(4) = 0.0d0
                s3(8) = s3(4)
            endif
            ! CRD 2019 fijamos el S3
            s3(4) = 1.0d0
            s3(8) = 1.0d0
            !  calculate sphi and dsphi
            !  condition here is on the bondlength rch(i)
            !  st argsphi is lt approx 19.0d0
            do i=1, 3
                if( rc1h(i) .lt. 3.8d0 )then
                    argsphi = aphi*(rc1h(i)-r0c1h)*exp(bphi*(rc1h(i)-cphi)**3)
                    sphi(i) = 1.0d0-tanh(argsphi)
                else
                    sphi(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                if( rc2h(i) .lt. 3.8d0 )then
                    argsphi   = aphi*(rc2h(i)-r0c2h)*exp(bphi*(rc2h(i)-cphi)**3)
                    sphi(i+4) = 1.0d0-tanh(argsphi)
                else
                    sphi(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            !      sphi(4)=1.0d0
            !      sphi(8)=1.0d0
            if( rcc .lt. 3.8d0 )then
                argsphic = aphi*(rcc-r0cc)*exp(bphi*(rcc-cphi)**3)
                sphi(4)  = 1.0d0-tanh(argsphic)
                sphi(8)  = sphi(4)
            else
                sphi(4) = 0.0d0
                sphi(8) = sphi(4)
            endif
            !  calculate stheta and dstheta
            do i=1, 3
                if( rc1h(i) .lt. 3.8d0 )then
                    argstheta = atheta*(rc1h(i)-r0c1h)*exp(btheta*(rc1h(i)-ctheta)**3)
                    stheta(i) = 1.0d0-tanh(argstheta)
                else
                    stheta(i) = 0.0d0
                endif
            enddo
            do i=1, 3
                if( rc2h(i) .lt. 3.8d0 )then
                    argstheta   = atheta*(rc2h(i)-r0c2h)*exp(btheta*(rc2h(i)-ctheta)**3)
                    stheta(i+4) = 1.0d0-tanh(argstheta)
                else
                    stheta(i+4) = 0.0d0
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            !      stheta(4)=1.0d0
            !      stheta(8)=1.0d0
            if( rcc .lt. 3.8d0 )then
                argsthetac = atheta*(rcc-r0cc)*exp(btheta*(rcc-ctheta)**3)
                stheta(4)  = 1.0d0-tanh(argsthetac)
                stheta(8)  = stheta(4)
            else
                stheta(4) = 0.0d0
                stheta(8) = stheta(4)
            endif
        end subroutine switchf

        !> @brief Evaluate prepare_local_context for the C2H6OH_VBMM potential.
        subroutine prepare_local_context
            !   N3TMMN = 3 * legacy_natoms
            !   legacy_natoms = the number of atoms represented by this potential function
            !   The variable N3TMMN is the minimum value of N3TM allowed to be
            !   passed by the calling routine for the number of cartesian
            !   coordinates needed to represent the full system represented by this
            !   potential energy surface routine.
            !   N3TM must be greater than or equal to N3TMMN.
            implicit none
            double precision :: d1bc, d3bc, fact1, fact2, fact3, vdwc2
            integer          :: i, icount, ind
            !C INCLUDE DEL VENUS
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            integer, parameter          :: n3tmmn = 30
            !c include CONST
            !c  include file for C2H6+H PES
            !      common /bonds/ r0ch
            !c call de venus
            !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
            if( legacy_natoms .gt. 25 )then
                write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            do i=1, 5
                ref(i) = ' '
            enddo
            !      REF(1)='J. Espinosa-Garcia and C. Rangel'
            !      REF(2)='No publicados, 2009'
            indexes(1)  = 6
            indexes(2)  = 6
            indexes(3)  = 1
            indexes(4)  = 1
            indexes(5)  = 1
            indexes(6)  = 1
            indexes(7)  = 1
            indexes(8)  = 1
            indexes(9)  = 8
            indexes(10) = 1
            irctnt      = 9
            !CALL POTINFO
            call ancvrt
            !  calculate indexes for coordinates
            do ind=1, 3
                icount   = ind-3
                nc1(ind) = 3*nnc1+icount
                nc2(ind) = 3*nnc2+icount
                nhb(ind) = 3*nnb+icount
                ! CRD 2019 añadimos el caluculo de las coordenadas de H
                no(ind) = 3*nno+icount
                ! fin CRD 2019
                do i=1, 6
                    nh(i,ind) = 3*nnh(i)+icount
                enddo
            enddo
            !  convert to appropriate units:
            !  From:
            !  kcal/mol
            !  mdyn A-1        -> 1.0d+05 j/mol...
            !  mdyn A rad-1
            !  To
            !  energy   in 1.0d+05 j/mol
            !  distance in 1.0d-10 m
            !  angles   in radians
            !  mass     in amu
            fact1 = 0.041840d0
            fact2 = 6.022045d0
            ! CRD 2019 fact3
            fact3 = 2.d0*3.1415926d0/360.d0
            ! fin CRD 2019
            d1ch  = d1ch*fact1
            d3ch  = d3ch*fact1
            d1c1h = d1c1h*fact1
            d3c1h = d3c1h*fact1
            d1c2h = d1c2h*fact1
            d3c2h = d3c2h*fact1
            d1cc  = d1cc*fact1
            !      d3cc=d3cc*fact1
            d1bc  = d1bc*fact1
            d3bc  = d3bc*fact1
            d1cb  = d1cb*fact1
            d3cb  = d3cb*fact1
            d1hh  = d1hh*fact1
            d3hh  = d3hh*fact1
            fch3  = fch3*fact2
            hch3  = hch3*fact2
            fkinf = fkinf*fact2
            ak    = ak*fact2
            ! CRD 2019
            fkh2oeq  = fkh2oeq*fact2
            angh2oeq = angh2oeq*fact3
            ! fin CRD 2019
            a1s   = a1s*1.d-8
            a2s   = a2s*1.d-8
            a1scc = a1scc*1.d-8
            a2scc = a2scc*1.d-8
            vdwc2 = vdwc2*1.d+4
        end subroutine prepare_local_context

        !> @brief Evaluate POTINFO for the C2H6OH_VBMM potential.
        subroutine potinfo
            implicit none
            integer                     :: i, j, k, kmax
            integer, parameter          :: n3atom = 75
            integer, parameter          :: natom  = 25
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = int(isurf*(isurf+1)/2)
            double precision, parameter :: pi     = 3.141592653589793d0
            write(nflag(18), '(/)')
            do i =1, 5
                write(nflag(18), '(2x,a75)') ref(i)
            enddo
            write(nflag(18), '(/)')
            kmax = 0
            do i = 1, isurf+1
                do j = 1, isurf+1
                    if( nasurf(i,j) .ne. 0.and.kmax .lt. max(i, j) ) kmax = max(i, j)
                enddo
            enddo
            write(nflag(18), '(2x,'' MAX. AND ACTUAL NO. OF EXCITED SURFACES: '',I3,5x,I3)') msurf, kmax-1
            if( kmax-1 .gt. msurf )then
                write(6, *) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
                stop
            endif
            ksdiag = 0
            kediag = 0
            do i = 2, isurf+1
                if( nasurf(i,i) .ne. 0 )then
                    kediag = i-1
                    if( ksdiag .eq. 0 ) ksdiag = i-1
                endif
            enddo
            ksoffd = 0
            keoffd = 0
            k      = 0
            do i = 1, isurf
                do j = i+1, isurf+1
                    k = k+1
                    if( nasurf(i,j)+nasurf(j,i) .ne. 0 )then
                        keoffd = k
                        if( ksoffd .eq. 0 ) ksoffd = k
                    endif
                enddo
            enddo
            write(nflag(18), '(2x,'' MAX. AND ACTUAL ORDER OF DERIVATIVES:    '',I3,5x,I3)') mder, nder
            if( nder .gt. mder )then
                write(6, *) ' WRONG INPUT ON ORDER OF DERIVATIVES'
                stop
            endif
            if( nflag(19) .eq. 1 )then
                write(nflag(18), '(/)')
                write(nflag(18), '(2x,''Cartesian coordinates are supplied by'',/, 2x,''the user in the array CART.'',//)')
                write(nflag(18), '(2x,''Provide cartesian coordinates in the'',/, 2x,''following order using the array CART'',' // &
                    '//, 2x,'' CART(1,1)...CART(1,3)   => ATOM 1'',/, 2x,'' CART(2,1)...CART(2,3)   => AT' // &
                    'OM 2'',/, 2x,'' CART(3,1)...CART(3,3)   => ATOM 3'',/, 2x,'' CART(N,1)...CART(N,3)  ' // &
                    ' => ATOM N'',/, 2x,''CART(25,1)...CART(25,3)  => ATOM 25'',/)')
                write(nflag(18), '(2x,''If the user wishes to relabel the atoms,'',/, 2x,''set the variable IREORDER equal to ' // &
                    '1'',/, 2x,''in the PARAMETER statement.  The user'',/, 2x,''must also specify the ne' // &
                    'w labeling'',/, 2x,''scheme.  This is done using the array'',/, 2x,''NULBL in the fo' // &
                    'llowing manner:'',//, 2x,''NULBL(i) = j'',/, 2x,''where:  i => old index'',/, 2x,'' ' // &
                    '       j => new index'',//)')
                write(nflag(18), '(2x,''Cartesian coordinates can be provided to'',/, 2x,''the potential routine in a variety ' // &
                    'of units.'',/, 2x,''The input units will be converted to Bohr'',/, 2x,''based on the' // &
                    ' following values of the NFLAG'',/, 2x,''variable:'',//, 2x,''NFLAG(1)  =  1  =>  CA' // &
                    'RTESIANS IN BOHR (no'',/, 2x,''                    conversion required)'',/, 2x,''NF' // &
                    'LAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS'',//)')
                write(nflag(18), '(2x,''The value of the energy and derivatives'',/, 2x,''(if computed) can be reported in a v' // &
                    'ariety'',/, 2x,''units.  A units conversion will take place'',/, 2x,''as specified b' // &
                    'y the following values of the'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(2) = 1 =>  ' // &
                    'ENERGIES REPORTED IN HARTEEE'',/, 2x,''NFLAG(2) = 2 =>  ENERGIES REPORTED IN mHARTRE' // &
                    'E'',/, 2x,''NFLAG(2) = 3 =>  ENERGIES REPORTED IN eV'',/, 2x,''NFLAG(2) = 4 =>  ENER' // &
                    'GIES REPORTED IN kcal/mol'',/, 2x,''NFLAG(2) = 5 =>  ENERGIES REPORTED IN cm**-1'',/' // &
                    '/)')
                write(nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
                    'he'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(1)=1 & NFLAG(2)=1 => DERIVATIVES REPOR' // &
                    'TED IN'',/, 2x,''                           HARTEEE/BOHR'',/, 2x,''NFLAG(1)=1 & NFLA' // &
                    'G(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/BOHR' // &
                    ''',/, 2x,''NFLAG(1)=1 & NFLAG(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''            ' // &
                    '               eV/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=4 => DERIVATIVES REPORTED IN' // &
                    ''',/, 2x,''                           kcal/mol/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=' // &
                    '5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm**-1/BOHR'',//)')
                write(nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
                    'he'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(1)=2 & NFLAG(2)=1 => DERIVATIVES REPOR' // &
                    'TED IN'',/, 2x,''                           HARTEEE/ANGSTROM'',/, 2x,''NFLAG(1)=2 & ' // &
                    'NFLAG(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/' // &
                    'ANGSTROM'',/, 2x,''NFLAG(1)=2 & NFLAG(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''    ' // &
                    '                       eV/ANGSTROM'',/, 2x,''NFLAG(1)=2 & NFLAG(2)=4 => DERIVATIVES ' // &
                    'REPORTED IN'',/, 2x,''                           kcal/mol/ANGSTROM'',/, 2x,''NFLAG(1' // &
                    ')=2 & NFLAG(2)=5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm*' // &
                    '*-1/ANGSTROM'',//)')
            endif
        end subroutine potinfo

        !> @brief Evaluate ANCVRT for the C2H6OH_VBMM potential.
        subroutine ancvrt
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi         = 3.141592653589793d0
            double precision, parameter :: clight     = 2.99792458d08
            double precision, parameter :: cmu_0      = 4.0d0*pi*1.0d-07
            double precision, parameter :: cepsilon_0 = 1.0d0/(cmu_0*clight**2)
            double precision, parameter :: ce         = 1.602176462d-19
            double precision, parameter :: cplanck    = 6.62606876d-34
            double precision, parameter :: cm_e       = 9.10938188d-31
            double precision, parameter :: cang       = 1.0d-10
            double precision, parameter :: cavogadro  = 6.02214199d23
            double precision, parameter :: ckcal      = 4.184d10
            double precision, parameter :: htomillh   = 1000.d0
            double precision, parameter :: htoev      = 27.2113834d0
            double precision, parameter :: htokcal    = 627.509470d0
            double precision, parameter :: htowave    = 219474.631d0
            double precision, parameter :: htokj      = 2625.49962d0
            double precision, parameter :: bohr_a     = .5291772083d0
            integer            :: i, ianum(7,32), inc1, inc2, ind, ind2, inu, isave(natom), isum, j, jsave(natom)
            character(len=3)   :: periodic_1(7,32)
            character(len=2)   :: name1(natom)
            character(len=2)   :: name2(natom)
            character(len=1)   :: iblank
            character(len=20)  :: distance
            character(len=20)  :: units
            do i=1, 7
                do j=1, 32
                    ianum(i,j)      = 0
                    periodic_1(i,j) = ' '
                enddo
            enddo
            distance         = 'BOHR                '
            units            = 'HARTREE             '
            ianum(1,1)       = 1
            ianum(1,32)      = 2
            ianum(2,1)       = 3
            ianum(2,2)       = 4
            ianum(2,27)      = 5
            ianum(2,28)      = 6
            ianum(2,29)      = 7
            ianum(2,30)      = 8
            ianum(2,31)      = 9
            ianum(2,32)      = 10
            ianum(3,1)       = 11
            ianum(3,2)       = 12
            ianum(3,27)      = 13
            ianum(3,28)      = 14
            ianum(3,29)      = 15
            ianum(3,30)      = 16
            ianum(3,31)      = 17
            ianum(3,32)      = 18
            ianum(4,1)       = 19
            ianum(4,2)       = 20
            ianum(4,17)      = 21
            ianum(4,18)      = 22
            ianum(4,19)      = 23
            ianum(4,20)      = 24
            ianum(4,21)      = 25
            ianum(4,22)      = 26
            ianum(4,23)      = 27
            ianum(4,24)      = 28
            ianum(4,25)      = 29
            ianum(4,26)      = 30
            ianum(4,27)      = 31
            ianum(4,28)      = 32
            ianum(4,29)      = 33
            ianum(4,30)      = 34
            ianum(4,31)      = 35
            ianum(4,32)      = 36
            ianum(5,1)       = 37
            ianum(5,2)       = 38
            ianum(5,17)      = 39
            ianum(5,18)      = 40
            ianum(5,19)      = 41
            ianum(5,20)      = 42
            ianum(5,21)      = 43
            ianum(5,22)      = 44
            ianum(5,23)      = 45
            ianum(5,24)      = 46
            ianum(5,25)      = 47
            ianum(5,26)      = 48
            ianum(5,27)      = 49
            ianum(5,28)      = 50
            ianum(5,29)      = 51
            ianum(5,30)      = 52
            ianum(5,31)      = 53
            ianum(5,32)      = 54
            ianum(6,1)       = 55
            ianum(6,2)       = 56
            ianum(6,3)       = 57
            ianum(6,4)       = 58
            ianum(6,5)       = 59
            ianum(6,6)       = 60
            ianum(6,7)       = 61
            ianum(6,8)       = 62
            ianum(6,9)       = 63
            ianum(6,10)      = 64
            ianum(6,11)      = 65
            ianum(6,12)      = 66
            ianum(6,13)      = 67
            ianum(6,14)      = 68
            ianum(6,15)      = 69
            ianum(6,16)      = 70
            ianum(6,17)      = 71
            ianum(6,18)      = 72
            ianum(6,19)      = 73
            ianum(6,20)      = 74
            ianum(6,21)      = 75
            ianum(6,22)      = 76
            ianum(6,23)      = 77
            ianum(6,24)      = 78
            ianum(6,25)      = 79
            ianum(6,26)      = 80
            ianum(6,27)      = 81
            ianum(6,28)      = 82
            ianum(6,29)      = 83
            ianum(6,30)      = 84
            ianum(6,31)      = 85
            ianum(6,32)      = 86
            ianum(7,1)       = 87
            ianum(7,2)       = 88
            ianum(7,3)       = 89
            ianum(7,4)       = 90
            ianum(7,5)       = 91
            ianum(7,6)       = 92
            ianum(7,7)       = 93
            ianum(7,8)       = 94
            ianum(7,9)       = 95
            ianum(7,10)      = 96
            ianum(7,11)      = 97
            ianum(7,12)      = 98
            ianum(7,13)      = 99
            ianum(7,14)      = 100
            ianum(7,15)      = 101
            ianum(7,16)      = 102
            ianum(7,17)      = 103
            ianum(7,18)      = 104
            ianum(7,19)      = 105
            ianum(7,20)      = 106
            ianum(7,21)      = 107
            ianum(7,22)      = 108
            ianum(7,23)      = 109
            ianum(7,24)      = 110
            ianum(7,25)      = 111
            ianum(7,26)      = 112
            ianum(7,27)      = 113
            ianum(7,28)      = 114
            ianum(7,29)      = 115
            ianum(7,30)      = 116
            ianum(7,31)      = 117
            ianum(7,32)      = 120
            periodic_1(1,1)  = 'H  '
            periodic_1(1,32) = 'He '
            periodic_1(2,1)  = 'Li '
            periodic_1(2,2)  = 'Be '
            periodic_1(2,27) = 'B  '
            periodic_1(2,28) = 'C  '
            periodic_1(2,29) = 'N  '
            periodic_1(2,30) = 'O  '
            periodic_1(2,31) = 'F  '
            periodic_1(2,32) = 'Ne '
            periodic_1(3,1)  = 'Na '
            periodic_1(3,2)  = 'Mg '
            periodic_1(3,27) = 'Al '
            periodic_1(3,28) = 'Si '
            periodic_1(3,29) = 'P  '
            periodic_1(3,30) = 'S  '
            periodic_1(3,31) = 'Cl '
            periodic_1(3,32) = 'Ar '
            periodic_1(4,1)  = 'K  '
            periodic_1(4,2)  = 'Ca '
            periodic_1(4,17) = 'Sc '
            periodic_1(4,18) = 'Ti '
            periodic_1(4,19) = 'V  '
            periodic_1(4,20) = 'Cr '
            periodic_1(4,21) = 'Mn '
            periodic_1(4,22) = 'Fe '
            periodic_1(4,23) = 'Co '
            periodic_1(4,24) = 'Ni '
            periodic_1(4,25) = 'Cu '
            periodic_1(4,26) = 'Zn '
            periodic_1(4,27) = 'Ga '
            periodic_1(4,28) = 'Ge '
            periodic_1(4,29) = 'As '
            periodic_1(4,30) = 'Se '
            periodic_1(4,31) = 'Br '
            periodic_1(4,32) = 'Kr '
            periodic_1(5,1)  = 'Rb '
            periodic_1(5,2)  = 'Sr '
            periodic_1(5,17) = 'Y  '
            periodic_1(5,18) = 'Zr '
            periodic_1(5,19) = 'Nb '
            periodic_1(5,20) = 'Mo '
            periodic_1(5,21) = 'Tc '
            periodic_1(5,22) = 'Ru '
            periodic_1(5,23) = 'Rh '
            periodic_1(5,24) = 'Pd '
            periodic_1(5,25) = 'Ag '
            periodic_1(5,26) = 'Cd '
            periodic_1(5,27) = 'In '
            periodic_1(5,28) = 'Sn '
            periodic_1(5,29) = 'Sb '
            periodic_1(5,30) = 'Te '
            periodic_1(5,31) = 'I  '
            periodic_1(5,32) = 'Xe '
            periodic_1(5,32) = 'Xe '
            do i=1, legacy_natoms
                isave(i) = 0
                jsave(i) = 0
                name1(i) = '  '
                name2(i) = '  '
            enddo
            iblank = ' '
            do ind=1, legacy_natoms
                do i=1, 7
                    do j=1, 32
                        if( indexes(ind) .eq. ianum(i,j) )then
                            isave(ind) = i
                            jsave(ind) = j
                        endif
                    enddo
                enddo
            enddo
            do ind=1, legacy_natoms
                ind2 = nulbl(ind)
                if( ind2 .eq. 0 ) ind2=ind
            enddo
            inc1 = 0
            do ind=1, irctnt-1
                inc1        = inc1+1
                name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
            enddo
            inc2 = 0
            do ind=irctnt, legacy_natoms
                inc2        = inc2+1
                name2(inc2) = periodic_1(isave(ind),jsave(ind))(:2)
            enddo
            if( nflag(1) .eq. 2 ) distance = 'ANGSTROMS           '
            if( nflag(2) .eq. 2 )then
                units = 'MILLIHARTREE        '
            elseif( nflag(2) .eq. 3 )then
                units = 'EV                  '
            elseif( nflag(2) .eq. 4 )then
                units = 'KCAL PER MOLE       '
            elseif( nflag(2) .eq. 5 )then
                units = 'WAVENUMBERS         '
            elseif( nflag(2) .eq. 6 )then
                units = 'KILOJOULES PER MOLE '
            endif
            cnvrtd  = 1.d0
            cnvrte  = 1.d0
            cnvrtde = 1.d0
            if( nflag(1) .eq. 2 ) cnvrtd = bohr_a
            if( nflag(2) .eq. 2 )then
                cnvrte = cnvrte*htomillh
            elseif( nflag(2) .eq. 3 )then
                cnvrte = cnvrte*htoev
            elseif( nflag(2) .eq. 4 )then
                cnvrte = cnvrte*htokcal
            elseif( nflag(2) .eq. 5 )then
                cnvrte = cnvrte*htowave
            elseif( nflag(2) .eq. 6 )then
                cnvrte = cnvrte*htokj
            endif
            cnvrtde = cnvrte/cnvrtd
            isum    = 0
            do inu=1, 25
                isum = isum + nulbl(inu)
            enddo
            ireorder = 0
            if( isum .ne. 0 ) ireorder = 1
        end subroutine ancvrt

        !> @brief Evaluate CARTOU for the C2H6OH_VBMM potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            implicit none
            integer            :: i, j
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            call convert_cartesian_input(cart, nulbl, legacy_natoms, cnvrtd, ireorder == 1, cartnu)
        end subroutine cartou

        !> @brief Evaluate CARTTOR for the C2H6OH_VBMM potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer            :: coordinate_status
            double precision   :: cmm, cphi, flm, hhd, hym, phi, pi, q, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
            integer            :: i, ind, k, l
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            call cartesian_to_potlib_coordinates(cartnu, legacy_natoms, icartr, r, coordinate_status)
            if( coordinate_status /= 0 )then
                write(nflag(18), '(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') icartr
                stop
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the C2H6OH_VBMM potential.
        subroutine eunitzero
            implicit none
            integer            :: i, j
            character(len=75)  :: ref(5)
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            pengygs = engygs * cnvrte - anuzero
            if( ksdiag .ne. 0 )then
                do i=ksdiag, kediag
                    pengyes(i) = engyes(i) * cnvrte - anuzero
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do j=ksoffd, keoffd
                    pengyij(j) = engyij(j) * cnvrte
                enddo
            endif
        end subroutine eunitzero

        !> @brief Evaluate RTOCART for the C2H6OH_VBMM potential.
        subroutine rtocart
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            double precision   :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, &
                & termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes(n3atom,isurf), &
                & ygs(n3atom), yij(n3atom,jsurf), za, zb, zc, zd
            integer            :: i, ind, j, j1, j2, k, m1
            if( icartr .eq. 1 )then
                do i = 1, legacy_natoms
                    ind            = 3*i-2
                    dgscartnu(i,1) = degsdr(ind)
                    dgscartnu(i,2) = degsdr(ind+1)
                    dgscartnu(i,3) = degsdr(ind+2)
                    if( ksdiag .ne. 0 )then
                        do j = ksdiag, kediag
                            descartnu(i,1,j) = deesdr(ind,j)
                            descartnu(i,2,j) = deesdr(ind+1,j)
                            descartnu(i,3,j) = deesdr(ind+2,j)
                        enddo
                    endif
                    if( keoffd .ne. 0 )then
                        do k = ksoffd, keoffd
                            dijcartnu(i,1,k) = deijdr(ind,k)
                            dijcartnu(i,2,k) = deijdr(ind+1,k)
                            dijcartnu(i,3,k) = deijdr(ind+2,k)
                        enddo
                    endif
                enddo
            elseif( icartr .eq. 2 )then
                do i = 1, legacy_natoms
                    dgscartnu(i,1) = 0.d0
                    dgscartnu(i,2) = 0.d0
                    dgscartnu(i,3) = 0.d0
                    if( ksdiag .ne. 0 )then
                        do j1=ksdiag, kediag
                            descartnu(i,1,j1) = 0.d0
                            descartnu(i,2,j1) = 0.d0
                            descartnu(i,3,j1) = 0.d0
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do j2=ksoffd, keoffd
                            dijcartnu(i,1,j2) = 0.d0
                            dijcartnu(i,2,j2) = 0.d0
                            dijcartnu(i,3,j2) = 0.d0
                        enddo
                    endif
                    do j = 1, legacy_natoms
                        if( j .lt. i )then
                            m1 = legacy_natoms*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
                            m1 = legacy_natoms*(i-1) - (i*(i-1))/2 + j-i
                        else
                            cycle
                        endif
                        y              = degsdr(m1)
                        termx          = (cartnu(i,1)-cartnu(j,1))/r(m1)
                        termy          = (cartnu(i,2)-cartnu(j,2))/r(m1)
                        termz          = (cartnu(i,3)-cartnu(j,3))/r(m1)
                        dgscartnu(i,1) = dgscartnu(i,1) + termx*y
                        dgscartnu(i,2) = dgscartnu(i,2) + termy*y
                        dgscartnu(i,3) = dgscartnu(i,3) + termz*y
                        if( ksdiag .gt. 0 )then
                            y = deesdr(m1,j1)
                            do j1=ksdiag, kediag
                                descartnu(i,1,j1) = descartnu(i,1,j1) + termx*y
                                descartnu(i,2,j1) = descartnu(i,2,j1) + termy*y
                                descartnu(i,3,j1) = descartnu(i,3,j1) + termz*y
                            enddo
                        elseif( ksoffd .gt. 0 )then
                            do j2=ksoffd, keoffd
                                y                 = deijdr(m1,j2)
                                dijcartnu(i,1,j2) = dijcartnu(i,1,j2) + termx*y
                                dijcartnu(i,2,j2) = dijcartnu(i,2,j2) + termy*y
                                dijcartnu(i,3,j2) = dijcartnu(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            elseif( icartr .eq. 3 )then
                do i = 1, legacy_natoms
                    ygs(i) = degsdr(i)/r(i)
                    if( ksdiag .ne. 0 )then
                        do j=ksdiag, kediag
                            yes(i,j) = deesdr(i,j)/r(i)
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do k=ksoffd, keoffd
                            yij(i,k) = deijdr(i,k)/r(i)
                        enddo
                    endif
                enddo
                do k = 1, 3
                    term12         = cartnu(1,k)-cartnu(2,k)
                    term23         = cartnu(2,k)-cartnu(3,k)
                    term13         = cartnu(1,k)-cartnu(3,k)
                    dgscartnu(1,k) = term12*ygs(1) + term13*ygs(3)
                    dgscartnu(2,k) = -term12*ygs(1) + term23*ygs(2)
                    dgscartnu(3,k) = -term13*ygs(3) - term23*ygs(2)
                    if( ksdiag .ne. 0 )then
                        do j1=ksdiag, kediag
                            descartnu(1,k,j1) = term12*yes(1,j1) + term13*yes(3,j1)
                            descartnu(2,k,j1) = -term12*yes(1,j1) + term23*yes(2,j1)
                            descartnu(3,k,j1) = -term13*yes(3,j1) - term23*yes(2,j1)
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do j2=ksoffd, keoffd
                            dijcartnu(1,k,j2) = term12*yij(1,j2) + term13*yij(3,j2)
                            dijcartnu(2,k,j2) = -term12*yij(1,j2) + term23*yij(2,j2)
                            dijcartnu(3,k,j2) = -term13*yij(3,j2) - term23*yij(2,j2)
                        enddo
                    endif
                enddo
            elseif( icartr .eq. 4 )then
                wh   = 1.007825d0
                wf   = 18.99840d0
                sum  = wh+wf
                eps  = wf/sum
                epsp = wh/sum
                u1   = cos(theta1)
                u2   = cos(theta2)
                u3   = cos(phi)
                ss1  = sin(theta1)
                ss2  = sin(theta2)
                ss3  = sin(phi)
                ya   = 0.0d0
                yb   = 0.0d0
                t0   = r1*u1
                za   = -epsp*t0
                zb   = eps*t0
                t0   = r1*ss1
                xa   = -epsp*t0
                xbb  = eps*t0
                t0   = r2*ss2
                t1   = t0*u3
                xc   = -epsp*t1
                xd   = eps*t1
                t1   = t0*ss3
                yc   = -epsp*t1
                yd   = eps*t1
                t0   = r2*u2
                zc   = -epsp*t0+rcm
                zd   = eps*t0+rcm
                rff  = sqrt((xa-xc)**2+yc**2+(za-zc)**2)
            else
                write(nflag(18), '(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') icartr
                stop
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the C2H6OH_VBMM potential.
        subroutine dedcou
            implicit none
            integer            :: i, j, k
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            if( ireorder .eq. 1 )then
                do i = 1, legacy_natoms
                    dgscart(i,1) = dgscartnu(nulbl(i),1) * cnvrtde
                    dgscart(i,2) = dgscartnu(nulbl(i),2) * cnvrtde
                    dgscart(i,3) = dgscartnu(nulbl(i),3) * cnvrtde
                    if( ksdiag .ne. 0 )then
                        do j=ksdiag, kediag
                            descart(i,1,j) = descartnu(nulbl(i),1,j) * cnvrtde
                            descart(i,2,j) = descartnu(nulbl(i),2,j) * cnvrtde
                            descart(i,3,j) = descartnu(nulbl(i),3,j) * cnvrtde
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do k=ksoffd, keoffd
                            dijcart(i,1,k) = dijcartnu(nulbl(i),1,k) * cnvrtde
                            dijcart(i,2,k) = dijcartnu(nulbl(i),2,k) * cnvrtde
                            dijcart(i,3,k) = dijcartnu(nulbl(i),3,k) * cnvrtde
                        enddo
                    endif
                enddo
            else
                do i = 1, legacy_natoms
                    dgscart(i,1) = dgscartnu(i,1) * cnvrtde
                    dgscart(i,2) = dgscartnu(i,2) * cnvrtde
                    dgscart(i,3) = dgscartnu(i,3) * cnvrtde
                    if( ksdiag .ne. 0 )then
                        do j=ksdiag, kediag
                            descart(i,1,j) = descartnu(i,1,j) * cnvrtde
                            descart(i,2,j) = descartnu(i,2,j) * cnvrtde
                            descart(i,3,j) = descartnu(i,3,j) * cnvrtde
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do k=ksoffd, keoffd
                            dijcart(i,1,k) = dijcartnu(i,1,k) * cnvrtde
                            dijcart(i,2,k) = dijcartnu(i,2,k) * cnvrtde
                            dijcart(i,3,k) = dijcartnu(i,3,k) * cnvrtde
                        enddo
                    endif
                enddo
            endif
        end subroutine dedcou
    end subroutine evaluate
end module chempotpy_c2h6oh_vbmm_pes

!> @brief Expose C2H6OH_VBMM through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_c2h6oh_vbmm_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 10
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
