!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_ch5_vbmm_cbg_2009_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_legacy_units, only: legacy_hartree_to_electronvolt
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 6
        integer, parameter            :: natom   = 25
        integer, parameter            :: isurf   = 5
        integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter            :: n3atom  = 75
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: a1ch
        double precision              :: a1s
        double precision              :: a2s
        double precision              :: a3s
        double precision              :: aa1
        double precision              :: aa2
        double precision              :: aa3
        double precision              :: aa4
        double precision              :: acb
        double precision              :: ahh
        double precision              :: ak
        double precision              :: anuzero
        double precision              :: aphi
        double precision              :: atheta
        double precision              :: b1ch
        double precision              :: b1s
        double precision              :: b2s
        double precision              :: b3s
        double precision              :: bk
        double precision              :: bphi
        double precision              :: btheta
        double precision              :: c1ch
        double precision              :: cart(natom,3)
        double precision              :: cartnu(natom,3)
        double precision              :: cnvrtd
        double precision              :: cnvrtde
        double precision              :: cnvrte
        double precision              :: cphi
        double precision              :: ctheta
        double precision              :: d1cb
        double precision              :: d1ch
        double precision              :: d1hh
        double precision              :: d3cb
        double precision              :: d3ch
        double precision              :: d3hh
        double precision              :: deesdr(n3atom,isurf)
        double precision              :: degsdr(n3atom)
        double precision              :: deijdr(n3atom,jsurf)
        double precision              :: descart(natom,3,isurf)
        double precision              :: descartnu(natom,3,isurf)
        double precision              :: dfdc(4,4,4)
        double precision              :: dfdelta(4,4)
        double precision              :: dfdh(4,4,4)
        double precision              :: dgscart(natom,3)
        double precision              :: dgscartnu(natom,3)
        double precision              :: dhdelta(4,4)
        double precision              :: dijcart(natom,3,jsurf)
        double precision              :: dijcartnu(natom,3,jsurf)
        double precision              :: ds1(4)
        double precision              :: ds2(4)
        double precision              :: ds3(4)
        double precision              :: dsphi(4)
        double precision              :: dstheta(4)
        double precision              :: dtheta0(4,4,4)
        double precision              :: engyes(isurf)
        double precision              :: engygs
        double precision              :: engyij(jsurf)
        double precision              :: ezero(isurf+1)
        double precision              :: f1(4)
        double precision              :: fch3
        double precision              :: fdelta(4)
        double precision              :: fk0(4,4)
        double precision              :: fkinf
        double precision              :: hch3
        double precision              :: hdelta(4)
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
        integer                       :: nder
        integer                       :: nflag(20)
        integer                       :: nh(4,3)
        integer                       :: nhb(3)
        integer                       :: nnb
        integer                       :: nnc
        integer                       :: nnh(4)
        integer                       :: nulbl(natom)
        double precision              :: pdot(150)
        double precision              :: pengyes(isurf)
        double precision              :: pengygs
        double precision              :: pengyij(jsurf)
        double precision              :: q(150)
        double precision              :: r(n3atom)
        double precision              :: r0cb
        double precision              :: r0ch
        double precision              :: r0hh
        double precision              :: rbh(4)
        double precision              :: rcb
        double precision              :: rch(4)
        character                     :: ref(5)
        double precision              :: s1(4)
        double precision              :: s2(4)
        double precision              :: s3(4)
        double precision              :: sphi(4)
        double precision              :: stheta(4)
        double precision              :: tbh(4,3)
        double precision              :: tcb(3)
        double precision              :: tch(4,3)
        double precision              :: theta0(4,4)
        double precision              :: v
        integer                       :: iatom, idir, istate
        !initialize
        a1ch      = 0.d0; a1s = 0.d0; a2s = 0.d0; a3s = 0.d0; aa1 = 0.d0; aa2 = 0.d0; aa3 = 0.d0; aa4 = 0.d0; acb = 0.d0
        ahh       = 0.d0; ak = 0.d0; anuzero = 0.d0; aphi = 0.d0; atheta = 0.d0; b1ch = 0.d0; b1s = 0.d0; b2s = 0.d0
        b3s       = 0.d0; bk = 0.d0; bphi = 0.d0; btheta = 0.d0; c1ch = 0.d0; cart = 0.d0; cartnu = 0.d0; cnvrtd = 0.d0
        cnvrtde   = 0.d0; cnvrte = 0.d0; cphi = 0.d0; ctheta = 0.d0; d1cb = 0.d0; d1ch = 0.d0; d1hh = 0.d0; d3cb = 0.d0
        d3ch      = 0.d0; d3hh = 0.d0; deesdr = 0.d0; degsdr = 0.d0; deijdr = 0.d0; descart = 0.d0; descartnu = 0.d0
        dfdc      = 0.d0; dfdelta = 0.d0; dfdh = 0.d0; dgscart = 0.d0; dgscartnu = 0.d0; dhdelta = 0.d0; dijcart = 0.d0
        dijcartnu = 0.d0; ds1 = 0.d0; ds2 = 0.d0; ds3 = 0.d0; dsphi = 0.d0; dstheta = 0.d0; dtheta0 = 0.d0
        engyes    = 0.d0; engygs = 0.d0; engyij = 0.d0; ezero = 0.d0; f1 = 0.d0; fch3 = 0.d0; fdelta = 0.d0; fk0 = 0.d0
        fkinf     = 0.d0; hch3 = 0.d0; hdelta = 0.d0; icartr = 0; indexes = 0; irctnt = 0; ireorder = 0; kediag = 0
        keoffd    = 0; ksdiag = 0; ksoffd = 0; legacy_natoms = 0; mder = 0; msurf = 0; nasurf = 0; nc = 0; nder = 0
        nflag     = 0; nh = 0; nhb = 0; nnb = 0; nnc = 0; nnh = 0; nulbl = 0; pdot = 0.d0; pengyes = 0.d0; pengygs = 0.d0
        pengyij   = 0.d0; q = 0.d0; r = 0.d0; r0cb = 0.d0; r0ch = 0.d0; r0hh = 0.d0; rbh = 0.d0; rcb = 0.d0; rch = 0.d0
        ref       = ' '
        s1        = 0.d0
        s2        = 0.d0
        s3        = 0.d0
        sphi      = 0.d0
        stheta    = 0.d0
        tbh       = 0.d0
        tcb       = 0.d0
        tch       = 0.d0
        theta0    = 0.d0
        nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0], shape(nasurf))
        nder          = 1
        nflag         = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        anuzero       = 0.0d0
        icartr        = 1
        msurf         = 0
        mder          = 1
        nulbl         = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        legacy_natoms = 6
        nnc           = 2
        nnb           = 6
        nnh           = [3, 4, 5, 1]
        r0ch          = 1.08898d0
        d1ch          = 111.266d0
        d3ch          = 48.96226d0
        a1ch          = 1.78374d0
        b1ch          = 0.14201d0
        c1ch          = 2.21773d0
        r0hh          = 0.74239d0
        d1hh          = 108.382d0
        d3hh          = 38.42657d0
        ahh           = 1.9589d0
        r0cb          = 1.08898d0
        d1cb          = 56.505d0
        d3cb          = 19.612d0
        acb           = 1.4173200d0
        a3s           = 0.1475300d0
        b3s           = -2.9926300d0
        aphi          = 0.5307000d0
        bphi          = 0.4012200d0
        cphi          = 1.9235100d0
        atheta        = 0.9119800d0
        btheta        = 0.3537500d0
        ctheta        = 1.8970500d0
        fch3          = 0.0693700d0
        hch3          = 0.1387400d0
        fkinf         = 0.4291400d0
        ak            = 0.1353000d0
        bk            = 10.7132d0
        aa1           = 1.265960d0
        aa2           = 0.000710d0
        aa3           = 0.985920d0
        aa4           = 2.785060d0
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
        else if( igrad == 1 )then
            do istate=1, nstates
                p(istate) = pengygs*legacy_hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = dgscart(iatom,idir)*51.422067
                enddo
            enddo
        else if( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
        contains
        !> @brief Evaluate POT for the CH5_VBMM_CBG_2009 potential.
        subroutine pot
            !
            !        This potential is written such that:
            !
            !                       X(1)  - X(3)  : X, Y, Z for H1
            !                       X(4)  - X(6)  : X, Y, Z for C
            !                       X(7)  - X(9)  : X, Y, Z for H3
            !                       X(10) - X(12) : X, Y, Z for H4
            !                       X(13) - X(15) : X, Y, Z for H2
            !                       X(16) - X(18) : X, Y, Z for H5
            !
            implicit none
            double precision :: en, vip, vop, vstr
            integer          :: i
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !      DIMENSION COORD(N3TM),DX(N3TM)
            !
            !
            call cartou
            call carttor
            !
            ! Changing units from angstrom to bohr and initialize
            !
            do i = 1, 18
                q(i)    = r(i)
                pdot(i) = 0.d0
            enddo
            !
            !  calculate relative coordinates and bond lengths
            !
            en = 0.0d0
            call coorden
            !
            !  calculate switching functions
            !
            call switchf
            !
            !  calculate reference angles and their derivatives
            !
            call refangles
            !
            !  calculate stretching potential
            !
            call stretch(vstr)
            !
            !  calculate out of plane bending potential
            !
            call opbend(vop)
            !
            !  calculate in plane bending potential
            !
            call ipbend(vip)
            !
            !  total potential energy is vstr+vop+vip
            !
            en = vstr+vop+vip
            !
            ! 0.03812 conversion factor from 10(5) j/mol to hartrees
            !
            en     = en*0.03812d0
            engygs = en
            !
            call eunitzero
            if( nder .ne. 0 )then
                !
                ! 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
                !
                do i=1, 18
                    degsdr(i) = pdot(i)*0.0201723d0
                enddo
                call rtocart
                call dedcou
            endif
            !
        end subroutine pot
        !> @brief Evaluate coorden for the CH5_VBMM_CBG_2009 potential.
        subroutine coorden
            !
            !  calculates relative coordinates and bond lengths
            !
            implicit none
            integer :: i, ind
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !
            !  calculate relative coordinates
            !
            do ind=1, 3
                tcb(ind) = q(nc(ind))-q(nhb(ind))
                do i=1, 4
                    tch(i,ind) = q(nc(ind))-q(nh(i,ind))
                    tbh(i,ind) = q(nhb(ind))-q(nh(i,ind))
                enddo
            enddo
            !
            !  calculate bond lengths
            !
            rcb = sqrt(tcb(1)*tcb(1)+tcb(2)*tcb(2)+tcb(3)*tcb(3))
            do i=1, 4
                rch(i) = sqrt(tch(i,1)*tch(i,1)+tch(i,2)*tch(i,2)+ tch(i,3)*tch(i,3))
                rbh(i) = sqrt(tbh(i,1)*tbh(i,1)+tbh(i,2)*tbh(i,2)+ tbh(i,3)*tbh(i,3))
            enddo
        end subroutine coorden
        !> @brief Evaluate refangles for the CH5_VBMM_CBG_2009 potential.
        subroutine refangles
            !
            !  subroutine calculates reference angles for the "in-plane" potential
            !
            implicit none
            double precision :: halfpi, tau, twopi
            integer          :: i, j, k
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !
            tau = acos(-1.0d0/3.0d0)
            !       pi=4.0d0*atan(1.0d0)
            halfpi = 0.5d0*pi
            twopi  = 2.0d0*pi
            !
            !  set diagonal elements to zero
            !
            do i=1, 4
                theta0(i,i) = 0.0d0
                do k=1, 4
                    dtheta0(i,i,k) = 0.0d0
                enddo
            enddo
            !
            !  calculate reference angles
            !
            theta0(1,2) = tau+(tau-halfpi)*(sphi(1)*sphi(2)-1.0d0) +(tau-twopi/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
            theta0(1,3) = tau+(tau-halfpi)*(sphi(1)*sphi(3)-1.0d0) +(tau-twopi/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
            theta0(1,4) = tau+(tau-halfpi)*(sphi(1)*sphi(4)-1.0d0) +(tau-twopi/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
            theta0(2,3) = tau+(tau-halfpi)*(sphi(2)*sphi(3)-1.0d0) +(tau-twopi/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
            theta0(2,4) = tau+(tau-halfpi)*(sphi(2)*sphi(4)-1.0d0) +(tau-twopi/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
            theta0(3,4) = tau+(tau-halfpi)*(sphi(3)*sphi(4)-1.0d0) +(tau-twopi/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
            !
            !  calculate the derivatives of theta0(i,j) in terms of rch(k)
            !  quantity calulated is dtheta0(i,j,k)
            !
            !  derivatives wrt rch(1)
            !
            dtheta0(1,2,1) = (tau-halfpi)*dsphi(1)*sphi(2)
            dtheta0(1,3,1) = (tau-halfpi)*dsphi(1)*sphi(3)
            dtheta0(1,4,1) = (tau-halfpi)*dsphi(1)*sphi(4)
            dtheta0(2,3,1) = (tau-twopi/3.0d0)*dstheta(1)*stheta(4)
            dtheta0(2,4,1) = (tau-twopi/3.0d0)*dstheta(1)*stheta(3)
            dtheta0(3,4,1) = (tau-twopi/3.0d0)*dstheta(1)*stheta(2)
            !
            !  derivatives wrt rch(2)
            !
            dtheta0(1,2,2) = (tau-halfpi)*sphi(1)*dsphi(2)
            dtheta0(1,3,2) = (tau-twopi/3.0d0)*dstheta(2)*stheta(4)
            dtheta0(1,4,2) = (tau-twopi/3.0d0)*dstheta(2)*stheta(3)
            dtheta0(2,3,2) = (tau-halfpi)*dsphi(2)*sphi(3)
            dtheta0(2,4,2) = (tau-halfpi)*dsphi(2)*sphi(4)
            dtheta0(3,4,2) = (tau-twopi/3.0d0)*stheta(1)*dstheta(2)
            !
            !  derivatives wrt rch(3)
            !
            dtheta0(1,2,3) = (tau-twopi/3.0d0)*dstheta(3)*stheta(4)
            dtheta0(1,3,3) = (tau-halfpi)*sphi(1)*dsphi(3)
            dtheta0(1,4,3) = (tau-twopi/3.0d0)*stheta(2)*dstheta(3)
            dtheta0(2,3,3) = (tau-halfpi)*sphi(2)*dsphi(3)
            dtheta0(2,4,3) = (tau-twopi/3.0d0)*stheta(1)*dstheta(3)
            dtheta0(3,4,3) = (tau-halfpi)*dsphi(3)*sphi(4)
            !
            !  derivatives wrt rch(4)
            !
            dtheta0(1,2,4) = (tau-twopi/3.0d0)*stheta(3)*dstheta(4)
            dtheta0(1,3,4) = (tau-twopi/3.0d0)*stheta(2)*dstheta(4)
            dtheta0(1,4,4) = (tau-halfpi)*sphi(1)*dsphi(4)
            dtheta0(2,3,4) = (tau-twopi/3.0d0)*stheta(1)*dstheta(4)
            dtheta0(2,4,4) = (tau-halfpi)*sphi(2)*dsphi(4)
            dtheta0(3,4,4) = (tau-halfpi)*sphi(3)*dsphi(4)
            !
            !  fill in the other half of the matrix
            !
            do i=1, 3
                do j=i+1, 4
                    theta0(j,i) = theta0(i,j)
                    do k=1, 4
                        dtheta0(j,i,k) = dtheta0(i,j,k)
                    enddo
                enddo
            enddo
        end subroutine refangles
        !> @brief Evaluate stretch for the CH5_VBMM_CBG_2009 potential.
        subroutine stretch( vstr )
            !
            !  subroutine to calculate leps-type stretching potential and its
            !  derivatives
            !
            implicit none
            double precision :: ach, achdc, achdh, arga, dumach, dumjbh, dumjcb, dumjch, dumjhi, dumjhj, dumqbh, dumqcb, dumqch, dumqhi, dumqhj, e1, e3, factj, rav, vj, vjbh, vjcb, vjch, vq, vqbh, vqcb, vqch, vstr
            integer          :: i, ind, j, k
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !
            dimension vqch(4), vjch(4), vqbh(4), vjbh(4), vq(4), vj(4), &
                & achdc(3), achdh(4, 3)
            !
            !  calculate avergage bond length for the methane moiety
            !
            rav = (rch(1)+rch(2)+rch(3)+rch(4))/4.0d0
            !
            !  initialise:
            !
            vstr = 0.0d0
            !
            !  ach:
            !
            !  in double precision tanh(19.0d0)=1.0d0 and we put the if statement
            !  in to avoid overflow/underflow errors
            !
            arga = c1ch*(rav-r0ch)
            if( arga .lt. 19.0d0 )then
                ach    = a1ch+b1ch*(tanh(arga)+1.0d0)*0.5d0
                dumach = b1ch*c1ch/(2.0d0*cosh(arga)**2)
            else
                ach    = a1ch+b1ch
                dumach = 0.0d0
            endif
            !
            !  calculate singlet: e1, triplet: e3 energies and vq and vj
            !  terms for each bond
            !
            e1   = d1cb*(exp(-2.0d0*acb*(rcb-r0cb))-2.0d0*exp(-acb*(rcb-r0cb)))
            e3   = d3cb*(exp(-2.0d0*acb*(rcb-r0cb))+2.0d0*exp(-acb*(rcb-r0cb)))
            vqcb = (e1+e3)*0.5d0
            vjcb = (e1-e3)*0.5d0
            do i=1, 4
                e1      = d1ch*(exp(-2.0d0*ach*(rch(i)-r0ch)) -2.0d0*exp(-ach*(rch(i)-r0ch)))
                e3      = d3ch*(exp(-2.0d0*ach*(rch(i)-r0ch)) +2.0d0*exp(-ach*(rch(i)-r0ch)))
                vqch(i) = (e1+e3)*0.5d0
                vjch(i) = (e1-e3)*0.5d0
                e1      = d1hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh)) -2.0d0*exp(-ahh*(rbh(i)-r0hh)))
                e3      = d3hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh)) +2.0d0*exp(-ahh*(rbh(i)-r0hh)))
                vqbh(i) = (e1+e3)*0.5d0
                vjbh(i) = (e1-e3)*0.5d0
                !
                !  calculate 3 body potential
                !
                vq(i) = vqch(i)+vqcb+vqbh(i)
                vj(i) = -sqrt(((vjch(i)-vjcb)**2+(vjcb-vjbh(i))**2 +(vjbh(i)-vjch(i))**2)*0.5d0)
                vstr  = vstr+vq(i)+vj(i)
            enddo
            !
            !  partial derivatives
            !  first we need the derivative of ach:
            !
            do ind=1, 3
                achdc(ind) = dumach*(tch(1,ind)/rch(1)+tch(2,ind)/rch(2) +tch(3,ind)/rch(3)+tch(4,ind)/rch(4))/4.0d0
                do i=1, 4
                    achdh(i, ind) = -dumach*tch(i,ind)/rch(i)/4.0d0
                enddo
            enddo
            dumqcb = -acb*((d1cb+d3cb)*exp(-2.0d0*acb*(rcb-r0cb))- (d1cb-d3cb)*exp(-acb*(rcb-r0cb)))/rcb
            !
            !  calculate cartesian derivatives:
            !  looping over ch(i) and bh(i)
            !
            do i=1, 4
                dumqbh = -ahh*((d1hh+d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh))- (d1hh-d3hh)*exp(-ahh*(rbh(i)-r0hh)))/rbh(i)
                factj  = 0.5d0/vj(i)
                dumjcb = -acb*((d1cb-d3cb)*exp(-2.0d0*acb*(rcb-r0cb)) -(d1cb+d3cb)*exp(-acb*(rcb-r0cb)))*factj/rcb
                dumjbh = -ahh*((d1hh-d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh)) -(d1hh+d3hh)*exp(-ahh*(rbh(i)-r0hh)))*factj/ &
                    & rbh(i)
                do ind=1, 3
                    !
                    !  deriv wrt hb:
                    !
                    pdot(nhb(ind))=pdot(nhb(ind)) -tcb(ind)*dumqcb+tbh(i,ind)*dumqbh +(vjch(i)-vjcb)*(dumjcb* &
                        & tcb(ind)) + (vjcb- vjbh(i))*(-dumjcb*tcb(ind)-dumjbh*tbh(i,ind)) +(vjbh(i)-vjch(i))*dumjbh* &
                        & tbh(i,ind)
                    !
                    !  dvqch(i)/dc
                    !
                    dumqch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch)) *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)- &
                        & r0ch)) - (d1ch- d3ch)*exp(-ach*(rch(i)-r0ch)))
                    pdot(nc(ind)) = pdot(nc(ind))+dumqch+tcb(ind)*dumqcb
                    !
                    !  dvqch(i)/dh(i)
                    !
                    dumqhi=(ach*tch(i,ind)/rch(i)-achdh(i, ind)*(rch(i)-r0ch)) *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)- &
                        & r0ch)) - (d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
                    pdot(nh(i,ind)) = pdot(nh(i,ind))+dumqhi-tbh(i,ind)*dumqbh
                    !
                    !  dvjch(i)/dc
                    !
                    dumjch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch)) *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)- &
                        & r0ch)) - (d1ch+ d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
                    !
                    !  dvj(i)/dnc(ind)
                    !
                    pdot(nc(ind))=pdot(nc(ind)) +(vjch(i)-vjcb)*(dumjch-dumjcb*tcb(ind)) +(vjcb-vjbh(i))*dumjcb* &
                        & tcb(ind) - (vjbh(i)-vjch(i))*dumjch
                    !
                    !  dvjch(i)/dh(i)
                    !
                    dumjhi=(ach*tch(i,ind)/rch(i)-achdh(i, ind)*(rch(i)-r0ch)) *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)- &
                        & r0ch)) - (d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
                    !
                    !  dvj(i)/dnh(i,ind)
                    !
                    pdot(nh(i,ind))=pdot(nh(i,ind)) +(vjch(i)-vjcb)*dumjhi +(vjcb-vjbh(i))*dumjbh*tbh(i,ind) + &
                        & (vjbh(i)- vjch(i))* (-dumjbh*tbh(i,ind)-dumjhi)
                    !
                    !  dv(i)/dh(j)
                    !
                    do k=1, 3
                        j = i+k
                        if( j .gt. 4 ) j=j-4
                        dumqhj=-achdh(j, ind)*(rch(i)-r0ch) *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch)) -(d1ch-d3ch)* &
                            & exp(-ach* (rch(i)-r0ch)))
                        dumjhj=-achdh(j, ind)*(rch(i)-r0ch) *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch)) -(d1ch+d3ch)* &
                            & exp(-ach* (rch(i)-r0ch)))*factj
                        pdot(nh(j,ind)) = pdot(nh(j,ind))+dumqhj +(vjch(i)-vjcb)*dumjhj -(vjbh(i)-vjch(i))*dumjhj
                    enddo
                enddo
            enddo
        end subroutine stretch
        !> @brief Evaluate opbend for the CH5_VBMM_CBG_2009 potential.
        subroutine opbend( vop )
            !
            !  subroutine calculates symmetrized vop potential and derivatives
            !
            implicit none
            double precision :: a, argd, axb, b, c, ddr, sum2, sum4, sumd2, sumd4, vop
            integer          :: i, ii, in, ind, itemp, j, k, l
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !
            double precision :: norma
            dimension sumd2(4), sumd4(4)
            dimension in(3), a(3), b(3), axb(3), c(4,3), argd(4)
            !
            !
            vop = 0.0d0
            !
            !  calculate force constants and their derivatives
            !
            call opforce
            !
            !  calculate out-of-plane angle and derivatives
            !
            do i=1, 4
                j = i+1
                if( j .gt. 4 ) j=j-4
                k = j+1
                if( k .gt. 4 ) k=k-4
                l = k+1
                if( l .gt. 4 ) l=l-4
                !
                !  modification to ensure that the set of methane CH bond vector
                !  (rj,rk,rl) is a right-handed set
                !
                in(1) = j
                in(2) = k
                in(3) = l
                !
                !  vector a is rk-rj, vector b is rl-rj
                !
                do ind=1, 3
                    a(ind) = q(nh(k,ind))-q(nh(j,ind))
                    b(ind) = q(nh(l,ind))-q(nh(j,ind))
                enddo
                !
                !  axb is vector a cross b
                !
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
                !
                !  c is position vector of h(ii): calculate c(j),c(k),c(l)
                !
                do ii=1, 3
                    do ind=1, 3
                        c(in(ii),ind) = -tch(in(ii),ind)/rch(in(ii))
                    enddo
                enddo
                !
                !  argd is the dot product axb dot c
                !
                do ii=1, 3
                    argd(in(ii)) = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                    argd(in(ii)) = argd(in(ii))/norma
                    !
                    !  if argd > 0 we need to switch vectors k and l around
                    !
                    if( argd(in(ii)) .gt. 0.d0 )then
                        itemp = k
                        k     = l
                        l     = itemp
                    endif
                enddo
                !
                !  subroutine performs sum over j, k, l
                !  sum2 = sum delta**2
                !  sum4 = sum delta**4
                !
                call calcdelta(i, j, k, l, sum2, sum4)
                sumd2(i) = sum2
                sumd4(i) = sum4
                vop      = vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
            enddo
            do i=1, 4
                do j=1, 4
                    !
                    !  overall derivatives of force constants i wrt the bond-length rch(j)
                    !
                    ddr = dfdelta(i,j)*sumd2(i)+dhdelta(i,j)*sumd4(i)
                    !
                    !  calculate derivatives in terms of cartesian coordinates:
                    !
                    do ind=1, 3
                        pdot(nh(j,ind)) = pdot(nh(j,ind))-tch(j,ind)*ddr/rch(j)
                        pdot(nc(ind))   = pdot(nc(ind))+tch(j,ind)*ddr/rch(j)
                    enddo
                enddo
            enddo
        end subroutine opbend
        !> @brief Evaluate ipbend for the CH5_VBMM_CBG_2009 potential.
        subroutine ipbend( vip )
            !
            !  subroutine calculates symmetrised in plane bend term
            !  and its derivatives
            !
            implicit none
            double precision :: costh, dth, dth0c, dth0k, dthc, dthi, dthj, termth, theta, vip
            integer          :: i, ind, j, k
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !
            dimension costh(4,4), theta(4,4), dth(4,4)
            !
            !  initialise
            !
            vip = 0.0d0
            !
            !  calculate force constants: fk0(i,j), f1(i)
            !  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
            !
            call ipforce
            !
            !  calculate theta(i,j) and in plane bend potential
            !
            do i=1, 3
                do j=i+1, 4
                    costh(i,j) = tch(i,1)*tch(j,1)+tch(i,2)*tch(j,2) +tch(i,3)*tch(j,3)
                    costh(i,j) = costh(i,j)/rch(i)/rch(j)
                    theta(i,j) = acos(costh(i,j))
                    dth(i,j)   = theta(i,j)-theta0(i,j)
                    vip        = vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
                    !
                    !  calculate partial derivatives wrt cartesian coordinates
                    !
                    !  calculate pdots wrt theta:
                    !
                    termth = -1.0d0/sqrt(1.0d0-costh(i,j)*costh(i,j))
                    do ind=1, 3
                        dthi            = -tch(j,ind)/rch(i)/rch(j) +costh(i,j)*tch(i,ind)/rch(i)/rch(i)
                        dthi            = dthi*termth
                        dthj            = -tch(i,ind)/rch(i)/rch(j) +costh(i,j)*tch(j,ind)/rch(j)/rch(j)
                        dthj            = dthj*termth
                        dthc            = -(dthi+dthj)
                        pdot(nh(i,ind)) = pdot(nh(i,ind)) +fk0(i,j)*f1(i)*f1(j)*dthi*dth(i,j)
                        pdot(nh(j,ind)) = pdot(nh(j,ind)) +fk0(i,j)*f1(i)*f1(j)*dthj*dth(i,j)
                        pdot(nc(ind))   = pdot(nc(ind)) +fk0(i,j)*f1(i)*f1(j)*dthc*dth(i,j)
                        do k=1, 4
                            !
                            !  calculate pdots wrt force constants and wrt theta0
                            !
                            dth0k = -dtheta0(i,j,k)*tch(k,ind)/rch(k)
                            dth0c = -dth0k
                            pdot(nh(k,ind))=pdot(nh(k,ind)) -0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k) -0.5d0* &
                                & tbh(k,ind)* dfdh(i,j,k)*dth(i,j)**2/rbh(k) -fk0(i,j)*f1(i)*f1(j)*dth0k*dth(i,j)
                            pdot(nc(ind))=pdot(nc(ind)) +0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k) -fk0(i,j)* &
                                & f1(i)* f1(j)* dth0c*dth(i,j)
                            pdot(nhb(ind)) = pdot(nhb(ind)) +0.5d0*tbh(k,ind)*dfdh(i,j,k)*dth(i,j)**2/rbh(k)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine ipbend
        !> @brief Evaluate calcdelta for the CH5_VBMM_CBG_2009 potential.
        subroutine calcdelta( i, j, k, l, sum2, sum4 )
            !
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !
            !   also calculates the derivatives wrt delta
            !
            implicit none
            double precision :: a, argd, atemp1, atemp2, atemp3, atemp4, atemp5, axb, b, c, cdot, daxb, deldot, delta, sum2, sum4
            integer          :: i, ii, in, ind, j, jj, k, l
            double precision :: norma
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !
            dimension  delta(4), in(3), a(3), b(3), axb(3), c(4,3), argd(4), &
                & daxb(4, 3, 3), cdot(4, 3, 3), atemp2(3)
            !
            !  initialise
            !
            sum2 = 0.0d0
            sum4 = 0.0d0
            !
            !  set j,k,l indices
            !
            in(1) = j
            in(2) = k
            in(3) = l
            !
            !  vector a is rk-rj, vector b is rl-rj
            !
            do ind=1, 3
                a(ind) = q(nh(k,ind))-q(nh(j,ind))
                b(ind) = q(nh(l,ind))-q(nh(j,ind))
            enddo
            !
            !  axb is vector a cross b
            !
            axb(1) = a(2)*b(3)-a(3)*b(2)
            axb(2) = a(3)*b(1)-a(1)*b(3)
            axb(3) = a(1)*b(2)-a(2)*b(1)
            norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
            norma  = sqrt(norma)
            !
            !  c is position vector of h(ii): calculate c(j),c(k),c(l)
            !
            do ii=1, 3
                do ind=1, 3
                    c(in(ii),ind) = -tch(in(ii),ind)/rch(in(ii))
                enddo
            enddo
            !
            !  argd is the dot product axb dot c
            !
            do ii=1, 3
                argd(in(ii))  = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                argd(in(ii))  = argd(in(ii))/norma
                delta(in(ii)) = acos(argd(in(ii)))-theta0(i,in(ii))
                !        write(*,*) 'theta,delta ',theta0(i,in(ii)),delta(in(ii))
                sum2 = sum2+delta(in(ii))**2
                sum4 = sum4+delta(in(ii))**4
            enddo
            !
            !  derivatives of axb wrt hj:
            !
            daxb(j, 1, 1) = 0.0d0
            daxb(j, 1, 2) = b(3)-a(3)
            daxb(j, 1, 3) = -b(2)+a(2)
            daxb(j, 2, 1) = -b(3)+a(3)
            daxb(j, 2, 2) = 0.0d0
            daxb(j, 2, 3) = b(1)-a(1)
            daxb(j, 3, 1) = b(2)-a(2)
            daxb(j, 3, 2) = -b(1)+a(1)
            daxb(j, 3, 3) = 0.0d0
            !
            !  derivatives of axb wrt hk:
            !
            daxb(k, 1, 1) = 0.0d0
            daxb(k, 1, 2) = -b(3)
            daxb(k, 1, 3) = b(2)
            daxb(k, 2, 1) = b(3)
            daxb(k, 2, 2) = 0.0d0
            daxb(k, 2, 3) = -b(1)
            daxb(k, 3, 1) = -b(2)
            daxb(k, 3, 2) = b(1)
            daxb(k, 3, 3) = 0.0d0
            !
            !  derivatives of axb wrt hl:
            !
            daxb(l, 1, 1) = 0.0d0
            daxb(l, 1, 2) = a(3)
            daxb(l, 1, 3) = -a(2)
            daxb(l, 2, 1) = -a(3)
            daxb(l, 2, 2) = 0.0d0
            daxb(l, 2, 3) = a(1)
            daxb(l, 3, 1) = a(2)
            daxb(l, 3, 2) = -a(1)
            daxb(l, 3, 3) = 0.0d0
            !
            !   loop over cdot(in(ii),ind,jind) where we consider deriv of c(in(ii))
            !   wrt h(in(ii),jind) with components jind
            !
            do ii=1, 3
                !
                !  deriv of cdot(in(ii),x) wrt x, y, z
                !
                cdot(in(ii), 1, 1) = 1.0d0/rch(in(ii)) +tch(in(ii),1)*c(in(ii),1)/rch(in(ii))**2
                cdot(in(ii), 1, 2) = tch(in(ii),2)*c(in(ii),1)/rch(in(ii))**2
                cdot(in(ii), 1, 3) = tch(in(ii),3)*c(in(ii),1)/rch(in(ii))**2
                !
                !  deriv of cdot(in(ii),y) wrt x, y, z
                !
                cdot(in(ii), 2, 1) = tch(in(ii),1)*c(in(ii),2)/rch(in(ii))**2
                cdot(in(ii), 2, 2) = 1.0d0/rch(in(ii)) +tch(in(ii),2)*c(in(ii),2)/rch(in(ii))**2
                cdot(in(ii), 2, 3) = tch(in(ii),3)*c(in(ii),2)/rch(in(ii))**2
                !
                !  deriv of cdot(in(ii),z) wrt x, y, z
                !
                cdot(in(ii), 3, 1) = tch(in(ii),1)*c(in(ii),3)/rch(in(ii))**2
                cdot(in(ii), 3, 2) = tch(in(ii),2)*c(in(ii),3)/rch(in(ii))**2
                cdot(in(ii), 3, 3) = 1.0d0/rch(in(ii)) +tch(in(ii),3)*c(in(ii),3)/rch(in(ii))**2
            enddo
            !
            do ii=1, 3
                do ind=1, 3
                    deldot = -dtheta0(i,in(ii),i)
                    !
                    !  derivative wrt h(i,ind)
                    !  for  rch(i) only terms are from the derivatives of theta0
                    !
                    deldot = -deldot*tch(i,ind)/rch(i)
                    pdot(nh(i,ind))=pdot(nh(i,ind)) +2.0d0*fdelta(i)*delta(in(ii))*deldot +4.0d0*hdelta(i)* &
                        & delta(in(ii))**3* deldot
                    !
                    !  derivative wrt c(ind)
                    !
                    deldot = -deldot
                    pdot(nc(ind)) = pdot(nc(ind)) +2.0d0*fdelta(i)*delta(in(ii))*deldot +4.0d0*hdelta(i)* &
                        & delta(in(ii))**3*deldot
                    do jj=1, 3
                        !
                        !  partial derivatives wrt h(in(jj),ind), loop over delta(i,in(ii))
                        !
                        !   atemp1 is axb dot daxb wrt h(in(jj))
                        !
                        atemp1 = axb(1)*daxb(in(jj), ind, 1) +axb(2)*daxb(in(jj), ind, 2) +axb(3)*daxb(in(jj), ind, 3)
                        atemp1 = atemp1/(norma**3)
                        !
                        !  atemp2 is deriv of normalised axb
                        !
                        atemp2(1) = daxb(in(jj), ind, 1)/norma-atemp1*axb(1)
                        atemp2(2) = daxb(in(jj), ind, 2)/norma-atemp1*axb(2)
                        atemp2(3) = daxb(in(jj), ind, 3)/norma-atemp1*axb(3)
                        !
                        !  atemp3 is daxb dot c(in(ii))
                        atemp3 = atemp2(1)*c(in(ii),1)+atemp2(2)*c(in(ii),2) +atemp2(3)*c(in(ii),3)
                        !
                        !  atemp4 is axb dot cdot
                        !
                        atemp4 = 0.0d0
                        if( ii .eq. jj )then
                            !
                            !  ie deriv of c(in(ii)) wrt h(in(jj)) is non zero only for ii = jj
                            !
                            atemp4 = axb(1)*cdot(in(ii), 1, ind) +axb(2)*cdot(in(ii), 2, ind) +axb(3)*cdot(in(ii), 3, &
                                & ind)
                            atemp4 = atemp4/norma
                        endif
                        !
                        !  atemp5 is deriv of theta0(i,in(ii)) wrt to nh(in(jj),ind)
                        !
                        atemp5 = -dtheta0(i,in(ii),in(jj))
                        !
                        !  deriv wrt h(in(jj)),ind):
                        !
                        atemp5 = -atemp5*tch(in(jj),ind)/rch(in(jj))
                        deldot = atemp3+atemp4
                        deldot = -1.0d0/sqrt(1.0d0-argd(in(ii))**2)*deldot
                        deldot = deldot+atemp5
                        pdot(nh(in(jj),ind))=pdot(nh(in(jj),ind)) +2.0d0*fdelta(i)*delta(in(ii))*deldot +4.0d0* &
                            & hdelta(i)* delta(in(ii))**3*deldot
                        !
                        !  for carbon the only contributions are from axb dot cdot term and
                        !  from theta0 and derivative cdot wrt carbon=-cdot wrt hydrogen
                        !
                        deldot = 1.0d0/sqrt(1.0d0-argd(in(ii))**2)*atemp4
                        deldot = deldot-atemp5
                        pdot(nc(ind))=pdot(nc(ind)) +2.0d0*fdelta(i)*delta(in(ii))*deldot +4.0d0*hdelta(i)* &
                            & delta(in(ii))**3* deldot
                    enddo
                enddo
            enddo
        end subroutine calcdelta
        !> @brief Evaluate opforce for the CH5_VBMM_CBG_2009 potential.
        subroutine opforce
            !
            !  calculates the out-of-plane bending force constants
            !  and their derivatives
            !
            implicit none
            double precision :: dswitch, switch
            integer          :: i, j
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !
            dimension switch(4), dswitch(4,4)
            !
            !  calculate switching functions:
            !
            switch(1) = (1.0d0-s3(1))*s3(2)*s3(3)*s3(4)
            switch(2) = (1.0d0-s3(2))*s3(3)*s3(4)*s3(1)
            switch(3) = (1.0d0-s3(3))*s3(4)*s3(1)*s3(2)
            switch(4) = (1.0d0-s3(4))*s3(1)*s3(2)*s3(3)
            !
            !  calculate derivatives:
            !  derivative of switch(1) wrt the 4 rch bond lengths
            !
            dswitch(1,1) = -ds3(1)*s3(2)*s3(3)*s3(4)
            dswitch(1,2) = (1.0d0-s3(1))*ds3(2)*s3(3)*s3(4)
            dswitch(1,3) = (1.0d0-s3(1))*s3(2)*ds3(3)*s3(4)
            dswitch(1,4) = (1.0d0-s3(1))*s3(2)*s3(3)*ds3(4)
            !
            !  derivative of switch(2) wrt the 4 rch bond lengths
            !
            dswitch(2,1) = (1.0d0-s3(2))*s3(3)*s3(4)*ds3(1)
            dswitch(2,2) = -ds3(2)*s3(3)*s3(4)*s3(1)
            dswitch(2,3) = (1.0d0-s3(2))*ds3(3)*s3(4)*s3(1)
            dswitch(2,4) = (1.0d0-s3(2))*s3(3)*ds3(4)*s3(1)
            !
            !  derivative of switch(3) wrt the 4 rch bond lengths
            !
            dswitch(3,1) = (1.0d0-s3(3))*s3(4)*ds3(1)*s3(2)
            dswitch(3,2) = (1.0d0-s3(3))*s3(4)*s3(1)*ds3(2)
            dswitch(3,3) = -ds3(3)*s3(4)*s3(1)*s3(2)
            dswitch(3,4) = (1.0d0-s3(3))*ds3(4)*s3(1)*s3(2)
            !
            !  derivative of switch(3) wrt the 4 rch bond lengths
            !
            dswitch(4,1) = (1.0d0-s3(4))*ds3(1)*s3(2)*s3(3)
            dswitch(4,2) = (1.0d0-s3(4))*s3(1)*ds3(2)*s3(3)
            dswitch(4,3) = (1.0d0-s3(4))*s3(1)*s3(2)*ds3(3)
            dswitch(4,4) = -ds3(4)*s3(1)*s3(2)*s3(3)
            !
            !  calculate the force constants and their derivatives
            !
            do i=1, 4
                fdelta(i) = switch(i)*fch3
                hdelta(i) = switch(i)*hch3
                do j=1, 4
                    dfdelta(i,j) = dswitch(i,j)*fch3
                    dhdelta(i,j) = dswitch(i,j)*hch3
                enddo
            enddo
        end subroutine opforce
        !> @brief Evaluate ipforce for the CH5_VBMM_CBG_2009 potential.
        subroutine ipforce
            !
            !  calculates the symmetrised in plane bend force constants and
            !  all partial derivatives involving them
            !
            implicit none
            double precision :: a1, a2, arga1, arga2, df1dc, df1dh, dfk0, duma1, duma2, f0, f2
            integer          :: i
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !
            dimension dfk0(4,4,4), df1dc(4), df1dh(4)
            !
            !  set force constant at asymptotes
            !
            f0 = fkinf+ak
            f2 = fkinf
            !
            fk0(1,2) = f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)*s2(4)-1.0d0)
            fk0(1,3) = f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)*s2(4)-1.0d0)
            fk0(1,4) = f0+f0*(s1(1)*s1(4)-1.0d0)+(f0-f2)*(s2(2)*s2(3)-1.0d0)
            fk0(2,3) = f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)*s2(4)-1.0d0)
            fk0(2,4) = f0+f0*(s1(2)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(3)-1.0d0)
            fk0(3,4) = f0+f0*(s1(3)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(2)-1.0d0)
            !
            !  derivative of fk0
            !
            dfk0(1,2,1) = f0*ds1(1)*s1(2)
            dfk0(1,2,2) = f0*s1(1)*ds1(2)
            dfk0(1,2,3) = (f0-f2)*ds2(3)*s2(4)
            dfk0(1,2,4) = (f0-f2)*s2(3)*ds2(4)
            !
            dfk0(1,3,1) = f0*ds1(1)*s1(3)
            dfk0(1,3,2) = (f0-f2)*ds2(2)*s2(4)
            dfk0(1,3,3) = f0*s1(1)*ds1(3)
            dfk0(1,3,4) = (f0-f2)*s2(2)*ds2(4)
            !
            dfk0(1,4,1) = f0*ds1(1)*s1(4)
            dfk0(1,4,2) = (f0-f2)*ds2(2)*s2(3)
            dfk0(1,4,3) = (f0-f2)*s2(2)*ds2(3)
            dfk0(1,4,4) = f0*s1(1)*ds1(4)
            !
            dfk0(2,3,1) = (f0-f2)*ds2(1)*s2(4)
            dfk0(2,3,2) = f0*ds1(2)*s1(3)
            dfk0(2,3,3) = f0*s1(2)*ds1(3)
            dfk0(2,3,4) = (f0-f2)*s2(1)*ds2(4)
            !
            dfk0(2,4,1) = (f0-f2)*ds2(1)*s2(3)
            dfk0(2,4,2) = f0*ds1(2)*s1(4)
            dfk0(2,4,3) = (f0-f2)*s2(1)*ds2(3)
            dfk0(2,4,4) = f0*s1(2)*ds1(4)
            !
            dfk0(3,4,1) = (f0-f2)*ds2(1)*s2(2)
            dfk0(3,4,2) = (f0-f2)*s2(1)*ds2(2)
            dfk0(3,4,3) = f0*ds1(3)*s1(4)
            dfk0(3,4,4) = f0*s1(3)*ds1(4)
            !
            do i=1, 4
                !
                !  calculate the terms f1(i)
                !
                arga1 = aa1*rbh(i)*rbh(i)
                arga2 = aa4*(rbh(i)-r0hh)*(rbh(i)-r0hh)
                a1    = 1.0d0-exp(-arga1)
                a2    = aa2+aa3*exp(-arga2)
                f1(i) = a1*exp(-a2*(rch(i)-r0ch)**2)
                !
                !  and calculate the derivatives wrt rch(i) and rbh(i)
                !
                duma1    = 2.0d0*aa1*rbh(i)*exp(-arga1)
                duma2    = -2.0d0*aa3*aa4*(rbh(i)-r0hh)*exp(-arga2)
                df1dc(i) = -2.0d0*(rch(i)-r0ch)*a1*a2*exp(-a2*(rch(i)-r0ch)**2)
                df1dh(i) = duma1*exp(-a2*(rch(i)-r0ch)**2) -duma2*(rch(i)-r0ch)**2*a1*exp(-a2*(rch(i)-r0ch)**2)
            enddo
            !
            !  derivative of total force constant f(i,j) wrt bond length rch(k)
            !  is given by dfdc(i,j,k)
            !
            dfdc(1,2,1) = dfk0(1,2,1)*f1(1)*f1(2)+fk0(1,2)*df1dc(1)*f1(2)
            dfdc(1,2,2) = dfk0(1,2,2)*f1(1)*f1(2)+fk0(1,2)*f1(1)*df1dc(2)
            dfdc(1,2,3) = dfk0(1,2,3)*f1(1)*f1(2)
            dfdc(1,2,4) = dfk0(1,2,4)*f1(1)*f1(2)
            !
            dfdc(1,3,1) = dfk0(1,3,1)*f1(1)*f1(3)+fk0(1,3)*df1dc(1)*f1(3)
            dfdc(1,3,2) = dfk0(1,3,2)*f1(1)*f1(3)
            dfdc(1,3,3) = dfk0(1,3,3)*f1(1)*f1(3)+fk0(1,3)*f1(1)*df1dc(3)
            dfdc(1,3,4) = dfk0(1,3,4)*f1(1)*f1(3)
            !
            dfdc(1,4,1) = dfk0(1,4,1)*f1(1)*f1(4)+fk0(1,4)*df1dc(1)*f1(4)
            dfdc(1,4,2) = dfk0(1,4,2)*f1(1)*f1(4)
            dfdc(1,4,3) = dfk0(1,4,3)*f1(1)*f1(4)
            dfdc(1,4,4) = dfk0(1,4,4)*f1(1)*f1(4)+fk0(1,4)*f1(1)*df1dc(4)
            !
            dfdc(2,3,1) = dfk0(2,3,1)*f1(2)*f1(3)
            dfdc(2,3,2) = dfk0(2,3,2)*f1(2)*f1(3)+fk0(2,3)*df1dc(2)*f1(3)
            dfdc(2,3,3) = dfk0(2,3,3)*f1(2)*f1(3)+fk0(2,3)*f1(2)*df1dc(3)
            dfdc(2,3,4) = dfk0(2,3,4)*f1(2)*f1(3)
            !
            dfdc(2,4,1) = dfk0(2,4,1)*f1(2)*f1(4)
            dfdc(2,4,2) = dfk0(2,4,2)*f1(2)*f1(4)+fk0(2,4)*df1dc(2)*f1(4)
            dfdc(2,4,3) = dfk0(2,4,3)*f1(2)*f1(4)
            dfdc(2,4,4) = dfk0(2,4,4)*f1(2)*f1(4)+fk0(2,4)*f1(2)*df1dc(4)
            !
            dfdc(3,4,1) = dfk0(3,4,1)*f1(3)*f1(4)
            dfdc(3,4,2) = dfk0(3,4,2)*f1(3)*f1(4)
            dfdc(3,4,3) = dfk0(3,4,3)*f1(3)*f1(4)+fk0(3,4)*df1dc(3)*f1(4)
            dfdc(3,4,4) = dfk0(3,4,4)*f1(3)*f1(4)+fk0(3,4)*f1(3)*df1dc(4)
            !
            !  derivative of total force constant f(i,j) wrt bond length rbh(k)
            !  is given by dfdh(i,j,k)
            !
            !  only non-zero derivatives are those from rbh(i) and rbh(j)
            !
            dfdh(1,2,1) = fk0(1,2)*df1dh(1)*f1(2)
            dfdh(1,2,2) = fk0(1,2)*f1(1)*df1dh(2)
            dfdh(1,2,3) = 0.0d0
            dfdh(1,2,4) = 0.0d0
            !
            dfdh(1,3,1) = fk0(1,3)*df1dh(1)*f1(3)
            dfdh(1,3,2) = 0.0d0
            dfdh(1,3,3) = fk0(1,3)*f1(1)*df1dh(3)
            dfdh(1,3,4) = 0.0d0
            !
            dfdh(1,4,1) = fk0(1,4)*df1dh(1)*f1(4)
            dfdh(1,4,2) = 0.0d0
            dfdh(1,4,3) = 0.0d0
            dfdh(1,4,4) = fk0(1,4)*f1(1)*df1dh(4)
            !
            dfdh(2,3,1) = 0.0d0
            dfdh(2,3,2) = fk0(2,3)*df1dh(2)*f1(3)
            dfdh(2,3,3) = fk0(2,3)*f1(2)*df1dh(3)
            dfdh(2,3,4) = 0.0d0
            !
            dfdh(2,4,1) = 0.0d0
            dfdh(2,4,2) = fk0(2,4)*df1dh(2)*f1(4)
            dfdh(2,4,3) = 0.0d0
            dfdh(2,4,4) = fk0(2,4)*f1(2)*df1dh(4)
            !
            dfdh(3,4,1) = 0.0d0
            dfdh(3,4,2) = 0.0d0
            dfdh(3,4,3) = fk0(3,4)*df1dh(3)*f1(4)
            dfdh(3,4,4) = fk0(3,4)*f1(3)*df1dh(4)
            !
        end subroutine ipforce
        !> @brief Evaluate switchf for the CH5_VBMM_CBG_2009 potential.
        subroutine switchf
            !
            !  calculates switching functions: s3,sphi,stheta
            !  and their derivatives ds3,dsphi,dstheta
            !
            implicit none
            double precision :: argmax, args1, args2, args3, argsphi, argstheta
            integer          :: i
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            !
            !
            !
            !
            !
            !
            !
            ! jcc valores originales
            !      a1s=1.5313681d-7
            !      b1s=-4.6696246d0
            !      a2s=1.0147402d-7
            !      b2s=-12.363798d0
            ! jcc valores ajustados a CCSDT
            a1s = 1.5132681d-7
            b1s = -4.3792246d0
            a2s = 1.9202402d-7
            b2s = -12.323018d0
            !
            !  use double precision criterion:
            !
            !  tanh(19.0d0)=1.0d0
            !
            argmax = 19.0d0
            do i=1, 4
                args1 = a1s*(rch(i)-r0ch)*(rch(i)-b1s)**8
                if( args1 .lt. argmax )then
                    s1(i)  = 1.0d0-tanh(args1)
                    ds1(i) = a1s*((rch(i)-b1s)**8 +8.0d0*(rch(i)-r0ch)*(rch(i)-b1s)**7)
                    ds1(i) = -ds1(i)/cosh(args1)**2
                else
                    s1(i)  = 0.0d0
                    ds1(i) = 0.0d0
                endif
                !
                args2 = a2s*(rch(i)-r0ch)*(rch(i)-b2s)**6
                if( args2 .lt. argmax )then
                    s2(i)  = 1.0d0-tanh(args2)
                    ds2(i) = a2s*((rch(i)-b2s)**6 +6.0d0*(rch(i)-r0ch)*(rch(i)-b2s)**5)
                    ds2(i) = -ds2(i)/cosh(args2)**2
                else
                    s2(i)  = 0.0d0
                    ds2(i) = 0.0d0
                endif
                !
                !  calculate s3 and ds3
                !
                args3 = a3s*(rch(i)-r0ch)*(rch(i)-b3s)**2
                if( args3 .lt. argmax )then
                    s3(i)  = 1.0d0-tanh(args3)
                    ds3(i) = a3s*(3.0d0*rch(i)**2-2.0d0*rch(i)*(r0ch+2.0d0*b3s) +b3s*(b3s+2.0d0*r0ch))
                    ds3(i) = -ds3(i)/cosh(args3)**2
                else
                    s3(i)  = 0.0d0
                    ds3(i) = 0.0d0
                endif
                !
                !  calculate sphi and dsphi
                !
                !  condition here is on the bondlength rch(i)
                !  st argsphi is lt approx 19.0d0
                !
                if( rch(i) .lt. 3.8d0 )then
                    argsphi  = aphi*(rch(i)-r0ch)*exp(bphi*(rch(i)-cphi)**3)
                    sphi(i)  = 1.0d0-tanh(argsphi)
                    dsphi(i) = aphi*(1.0d0+3.0d0*bphi*(rch(i)-r0ch) *(rch(i)-cphi)**2)
                    dsphi(i) = dsphi(i)*exp(bphi*(rch(i)-cphi)**3)
                    dsphi(i) = -dsphi(i)/cosh(argsphi)**2
                else
                    sphi(i)  = 0.0d0
                    dsphi(i) = 0.0d0
                endif
                !
                !  calculate stheta and dstheta
                !
                if( rch(i) .lt. 3.8d0 )then
                    argstheta  = atheta*(rch(i)-r0ch)*exp(btheta*(rch(i)-ctheta)**3)
                    stheta(i)  = 1.0d0-tanh(argstheta)
                    dstheta(i) = atheta*(1.0d0+3.0d0*btheta*(rch(i)-r0ch) *(rch(i)-ctheta)**2)
                    dstheta(i) = dstheta(i)*exp(btheta*(rch(i)-ctheta)**3)
                    dstheta(i) = -dstheta(i)/cosh(argstheta)**2
                else
                    stheta(i)  = 0.0d0
                    dstheta(i) = 0.0d0
                endif
            enddo
        end subroutine switchf
        !> @brief Evaluate prepare_local_context for the CH5_VBMM_CBG_2009 potential.
        subroutine prepare_local_context
            !
            !   N3TMMN = 3 * legacy_natoms
            !   legacy_natoms = the number of atoms represented by this potential function
            !
            !   The variable N3TMMN is the minimum value of N3TM allowed to be
            !   passed by the calling routine for the number of cartesian
            !   coordinates needed to represent the full system represented by this
            !   potential energy surface routine.
            !   N3TM must be greater than or equal to N3TMMN.
            !
            implicit none
            double precision :: fact1, fact2
            integer          :: i, icount, ind
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            integer, parameter          :: n3tmmn = 18
            !
            !
            !
            !
            !
            !
            !
            !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
            !
            if( legacy_natoms .gt. 25 )then
                write(nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !
            do i=1, 5
                ref(i) = ' '
            enddo
            !
            ref(1) = 'J. Espinosa-Garcia and J. C. Corchado'
            ref(2) = 'Unpublished results, 1999.'
            !
            indexes(1) = 1
            indexes(2) = 6
            indexes(3) = 1
            indexes(4) = 1
            indexes(5) = 1
            indexes(6) = 1
            !
            !
            !
            irctnt = 6
            !
            !      CALL POTINFO
            !
            call ancvrt
            !
            !  calculate indexes for coordinates
            !
            do ind=1, 3
                icount   = ind-3
                nc(ind)  = 3*nnc+icount
                nhb(ind) = 3*nnb+icount
                do i=1, 4
                    nh(i,ind) = 3*nnh(i)+icount
                enddo
            enddo
            !
            !  convert to appropriate units:
            !
            !  energy   in 1.0d+05 j/mol
            !  time     in 1.0d-14 s
            !  distance in 1.0d-10 m
            !  angles   in radians
            !  mass     in amu
            !
            fact1 = 0.041840d0
            fact2 = 6.022045d0
            !
            d1ch  = d1ch*fact1
            d3ch  = d3ch*fact1
            d1cb  = d1cb*fact1
            d3cb  = d3cb*fact1
            d1hh  = d1hh*fact1
            d3hh  = d3hh*fact1
            fch3  = fch3*fact2
            hch3  = hch3*fact2
            fkinf = fkinf*fact2
            ak    = ak*fact2
            !
        end subroutine prepare_local_context
        !> @brief Evaluate ANCVRT for the CH5_VBMM_CBG_2009 potential.
        subroutine ancvrt
            implicit none
            integer            :: i, ianum, inc1, inc2, ind, ind2, inu, isave, isum, j, jsave
            character(len=3)   :: periodic_1(7,32)
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            character(len=2)   :: name1(natom)
            character(len=2)   :: name2(natom)
            character(len=1)   :: iblank
            character(len=20)  :: distance
            character(len=20)  :: units
            dimension ianum(7,32)
            dimension isave(natom), jsave(natom)
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
            else if( nflag(2) .eq. 3 )then
                units = 'EV                  '
            else if( nflag(2) .eq. 4 )then
                units = 'KCAL PER MOLE       '
            else if( nflag(2) .eq. 5 )then
                units = 'WAVENUMBERS         '
            else if( nflag(2) .eq. 6 )then
                units = 'KILOJOULES PER MOLE '
            endif
            cnvrtd  = 1.d0
            cnvrte  = 1.d0
            cnvrtde = 1.d0
            if( nflag(1) .eq. 2 ) cnvrtd = bohr_a
            if( nflag(2) .eq. 2 )then
                cnvrte = cnvrte*htomillh
            else if( nflag(2) .eq. 3 )then
                cnvrte = cnvrte*htoev
            else if( nflag(2) .eq. 4 )then
                cnvrte = cnvrte*htokcal
            else if( nflag(2) .eq. 5 )then
                cnvrte = cnvrte*htowave
            else if( nflag(2) .eq. 6 )then
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
        !> @brief Evaluate CARTOU for the CH5_VBMM_CBG_2009 potential.
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
        !> @brief Evaluate CARTTOR for the CH5_VBMM_CBG_2009 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer            :: coordinate_status
            double precision   :: cmm, flm, hhd, hym, phi, pi, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
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
        !> @brief Evaluate EUNITZERO for the CH5_VBMM_CBG_2009 potential.
        subroutine eunitzero
            implicit none
            integer            :: i, j
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
        !> @brief Evaluate RTOCART for the CH5_VBMM_CBG_2009 potential.
        subroutine rtocart
            implicit none
            double precision   :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes, ygs, yij, za, zb, zc, zd
            integer            :: i, ind, j, j1, j2, k, m1
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            dimension ygs(n3atom), yes(n3atom,isurf), yij(n3atom,jsurf)
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
            else if( icartr .eq. 2 )then
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
                        else if( j .gt. i )then
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
                        else if( ksoffd .gt. 0 )then
                            do j2=ksoffd, keoffd
                                y                 = deijdr(m1,j2)
                                dijcartnu(i,1,j2) = dijcartnu(i,1,j2) + termx*y
                                dijcartnu(i,2,j2) = dijcartnu(i,2,j2) + termy*y
                                dijcartnu(i,3,j2) = dijcartnu(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            else if( icartr .eq. 3 )then
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
            else if( icartr .eq. 4 )then
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
        !> @brief Evaluate DEDCOU for the CH5_VBMM_CBG_2009 potential.
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
end module chempotpy_ch5_vbmm_cbg_2009_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_ch5_vbmm_cbg_2009_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 6
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
