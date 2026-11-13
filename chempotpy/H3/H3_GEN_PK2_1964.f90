!> @brief Implement H3_GEN_PK2_1964, a potential energy surface for H3.
module chempotpy_h3_gen_pk2_1964_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3_GEN_PK2_1964 potential energy surface for H3.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units, only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter    :: nstates               = 1
        integer, parameter    :: number_of_atoms       = 3
        integer, parameter    :: natom                 = 25
        integer, parameter    :: n3atom                = 75
        integer, parameter    :: isurf                 = 5
        integer, parameter    :: jsurf                 = isurf*(isurf+1)/2
        real(wp), parameter   :: legacy_ev_per_hartree = 27.21161e0_wp
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(number_of_atoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,number_of_atoms,3)
        real(wp), intent(out) :: d(nstates,nstates,number_of_atoms,3)
        real(wp)              :: alph, anuzero, beta, cnvrtd, cnvrtde, cnvrte, d1, d3, del, easyab
        real(wp)              :: easyac, easybc, engygs, eps, fkap, gam, pengygs, re, v
        real(wp)              :: cart(natom,3), cartnu(natom,3), deesdr(n3atom,isurf)
        real(wp)              :: degsdr(n3atom), deijdr(n3atom,jsurf), descart(natom,3,isurf)
        real(wp)              :: descartnu(natom,3,isurf), dgscart(natom,3), dgscartnu(natom,3)
        real(wp)              :: dijcart(natom,3,jsurf), dijcartnu(natom,3,jsurf), engyes(isurf)
        real(wp)              :: engyij(jsurf), ezero(isurf+1), pengyes(isurf), pengyij(jsurf)
        real(wp)              :: r(n3atom)
        integer               :: iatom, icartr, idir, indexes(natom), irctnt, ireorder, istate
        integer               :: kediag, keoffd, ksdiag, ksoffd, mder, msurf, natoms, nder
        integer               :: nasurf(isurf+1,isurf+1), nflag(20), nulbl(natom)
        character(len=75)     :: ref(5)
        v       = 0.e0_wp; p = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        alph    = 0.e0_wp; anuzero = 0.e0_wp; beta = 0.e0_wp; cnvrtd = 0.e0_wp; cnvrtde = 0.e0_wp; cnvrte = 0.e0_wp
        d1      = 0.e0_wp; d3 = 0.e0_wp; del = 0.e0_wp; easyab = 0.e0_wp; easyac = 0.e0_wp; easybc = 0.e0_wp
        engygs  = 0.e0_wp; eps = 0.e0_wp; fkap = 0.e0_wp; gam = 0.e0_wp; pengygs = 0.e0_wp; re = 0.e0_wp
        cart    = 0.e0_wp; cartnu = 0.e0_wp; deesdr = 0.e0_wp; degsdr = 0.e0_wp; deijdr = 0.e0_wp
        descart = 0.e0_wp; descartnu = 0.e0_wp; dgscart = 0.e0_wp; dgscartnu = 0.e0_wp
        dijcart = 0.e0_wp; dijcartnu = 0.e0_wp; engyes = 0.e0_wp; engyij = 0.e0_wp; ezero = 0.e0_wp
        pengyes = 0.e0_wp; pengyij = 0.e0_wp; r = 0.e0_wp
        indexes = 0; irctnt = 0; ireorder = 0; kediag = 0; keoffd = 0; ksdiag = 0; ksoffd = 0
        mder    = 1; msurf = 0; natoms = number_of_atoms; nder = 0; nasurf = 0; nasurf(1,1) = 1
        nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        nulbl   = 0; icartr = 3; ref = ' '
        d1      = 4.7466e0_wp; d3 = 1.9668e0_wp; re = 1.40083e0_wp; alph = 1.04435e0_wp; beta = 1.000122e0_wp
        eps     = -17.5e0_wp; del = 28.2e0_wp; fkap = 0.60e0_wp; gam = 0.65e0_wp
        do iatom=1, number_of_atoms
            do idir=1, 3
                cart(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        nder         = igrad
        d1           = d1/legacy_ev_per_hartree; d3 = d3/legacy_ev_per_hartree
        eps          = eps/legacy_ev_per_hartree; del = del/legacy_ev_per_hartree
        easyab       = d1; easybc = easyab; easyac = easyab; ezero(1) = d1
        ref(1)       = 'R. N. Porter and M. Karplus'
        ref(2)       = 'J. Chem. Phys. 40, 1105(1964)'
        indexes(1:3) = 1; irctnt = 2
        call ancvrt
        call evaluate_kernel
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = pengygs*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = pengygs*hartree_to_electronvolt
            enddo
            do iatom=1, number_of_atoms
                do idir=1, 3
                    g(1,iatom,idir) = dgscart(iatom,idir)*hartree_per_bohr_to_ev_per_angstrom
                enddo
            enddo
        elseif( igrad == 2 )then
            write(*, *) 'Only energy and gradient are available'
        endif

        contains

        !> @brief Evaluate the local PK2 potential kernel.
        subroutine evaluate_kernel
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the H2 diatomic, with the
            !   H2 diatomic at its equilibrium configuration.  For this potential en
            !   surface the potential energy in the three asymptotic valleys are equ
            !   The potential energy is defined to be equal to zero when one H is
            !   "infinitely" far from the H2 diatomic, and the H2 diatomic bond
            !   length is equal to the H2 equilibrium bond length.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a, a1, a2, a3, b, b1, b2, b3, bac, c, c1, c2, c3, cs123, de1(3), de3(3), delj(3), &
                & dife(3), dj(3), dj112, dj113, dj12, dj123, dj13, dj21, dj212, dj213, dj223, dj23, dj31, dj312, dj313, dj32, dj323, &
                & dq(3), ds(3), e1(3), e3(3), engygs0, ex(3), exa, exb, expg, expr, f, fj(3), fj12, fj13, fj23, gamr, q(3), qt, &
                & rbac, rel, rho, rr(3), s(3), s12, s123, s13, s23, ss(3), sume(3), t, zeta
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !   Conversion constants
            !   CEV (eV) converts from Hartree to kcal/mol
            real(wp), parameter :: cev    = 27.21161e0_wp
            real(wp), parameter :: hartre = 627.5095e0_wp
            call cartou
            call carttor
            !   Check the values of NASURF and NDER for validity.
            do i = 1, 3
                gamr    = gam*r(i)
                expg    = exp(-gamr)
                zeta    = 1.0e0_wp + fkap*expg
                rho     = r(i)*zeta
                rr(i)   = 1.0e0_wp/r(i)
                ex(i)   = del*exp(-r(i)-r(i))
                expr    = exp(-rho)
                s(i)    = (1.0e0_wp+rho+(rho*rho)/3.0e0_wp)*expr
                ss(i)   = s(i)*s(i)
                rel     = re-r(i)
                exa     = exp(alph*rel)
                exb     = exp(beta*rel)
                e1(i)   = d1*(1.0e0_wp-exa)**2-d1
                e3(i)   = d3*(1.0e0_wp+exb)**2-d3
                sume(i) = e1(i)+e3(i)
                dife(i) = e1(i)-e3(i)
                q(i)    = 0.5e0_wp*(sume(i)+ss(i)*dife(i))
                fj(i)   = 0.5e0_wp*(dife(i)+ss(i)*sume(i))
                de1(i)  = 2.0e0_wp*d1*alph*exa*(1.0e0_wp-exa)
                de3(i)  = 2.0e0_wp*d3*beta*exb*(-1.0e0_wp-exb)
                ds(i)   = -rho*(1.0e0_wp+rho)*expr*(zeta+gamr*(1.0e0_wp-zeta))/3.0e0_wp
                dj(i)   = s(i)*sume(i)*ds(i)+0.5e0_wp*((1.0e0_wp+ss(i))*de1(i)-(1.0e0_wp-ss(i ))* de3(i))
                delj(i) = (1.0e0_wp+rr(i))*ex(i)
                dq(i)   = s(i)*dife(i)*ds(i)+0.5e0_wp*((1.0e0_wp+ss(i))*de1(i)+(1.0e0_wp-ss(i ))* de3(i))
            enddo
            dj(1)   = dj(1) + 2.0e0_wp*s(1)*ds(1)*(delj(2)+delj(3))
            dj(2)   = dj(2) + 2.0e0_wp*s(2)*ds(2)*(delj(1)+delj(3))
            dj(3)   = dj(3) + 2.0e0_wp*s(3)*ds(3)*(delj(1)+delj(2))
            fj(1)   = fj(1) + ss(1)*(delj(2)+delj(3))
            fj(2)   = fj(2) + ss(2)*(delj(1)+delj(3))
            fj(3)   = fj(3) + ss(3)*(delj(1)+delj(2))
            s12     = ss(1)-ss(2)
            s23     = ss(2)-ss(3)
            s13     = s12+s23
            fj12    = fj(1)-fj(2)
            fj23    = fj(2)-fj(3)
            fj13    = fj12+fj23
            dj12    = -ss(2)*ex(1)*(2.0e0_wp+2.0e0_wp*rr(1)+rr(1)*rr(1))
            dj13    = -ss(3)*ex(1)*(2.0e0_wp+2.0e0_wp*rr(1)+rr(1)*rr(1))
            dj21    = -ss(1)*ex(2)*(2.0e0_wp+2.0e0_wp*rr(2)+rr(2)*rr(2))
            dj23    = -ss(3)*ex(2)*(2.0e0_wp+2.0e0_wp*rr(2)+rr(2)*rr(2))
            dj31    = -ss(1)*ex(3)*(2.0e0_wp+2.0e0_wp*rr(3)+rr(3)*rr(3))
            dj32    = -ss(2)*ex(3)*(2.0e0_wp+2.0e0_wp*rr(3)+rr(3)*rr(3))
            dj112   = dj(1)-dj12
            dj113   = dj(1)-dj13
            dj212   = dj21-dj(2)
            dj223   = dj(2)-dj23
            dj313   = dj31-dj(3)
            dj323   = dj32-dj(3)
            dj123   = dj12-dj13
            dj213   = dj21-dj23
            dj312   = dj31-dj32
            s123    = s(1)*s(2)*s(3)
            qt      = q(1)+q(2)+q(3)-eps*s123
            cs123   = 1.0e0_wp-s123
            a       = cs123*cs123-0.5e0_wp*(s12*s12+s23*s23+s13*s13)
            b       = -qt*cs123+0.5e0_wp*(fj12*s12+fj23*s23+fj13*s13)
            c       = qt*qt-0.5e0_wp*(fj12*fj12+fj23*fj23+fj13*fj13)
            a1      = -(ds(1)+ds(1))*(cs123*s(2)*s(3)+s(1)*(s12+s13))
            a2      = -(ds(2)+ds(2))*(cs123*s(1)*s(3)+s(2)*(s23-s12))
            a3      = -(ds(3)+ds(3))*(cs123*s(1)*s(2)-s(3)*(s13+s23))
            f       = cs123*eps+qt
            b1      = -cs123*dq(1)+ds(1)*(s(2)*s(3)*f+s(1)*(fj12+fj13)) +0.5e0_wp*(s12*dj112+s23*dj123+s13*dj113)
            b2      = -cs123*dq(2)+ds(2)*(s(1)*s(3)*f+s(2)*(fj23-fj12)) +0.5e0_wp*(s12*dj212+s23*dj223+s13*dj213)
            b3      = -cs123*dq(3)+ds(3)*(s(1)*s(2)*f-s(3)*(fj13+fj23)) +0.5e0_wp*(s12*dj312+s23*dj323+s13*dj313)
            f       = qt+qt
            c1      = f*(dq(1)-eps*s(2)*s(3)*ds(1))-(fj12*dj112+fj23*dj123 +fj13*dj113)
            c2      = f*(dq(2)-eps*s(1)*s(3)*ds(2))-(fj12*dj212+fj23*dj223 +fj13*dj213)
            c3      = f*(dq(3)-eps*s(2)*s(1)*ds(3))-(fj12*dj312+fj23*dj323 +fj13*dj313)
            bac     = b*b-a*c
            rbac    = sqrt(bac)
            engygs0 = -b - rbac
            if( engygs0 .ne. 0.e0_wp ) engygs0 = engygs0/a
            engygs = engygs0 + ezero(1)
            call eunitzero
            if( nder .ne. 0 )then
                t                                      = 2.0e0_wp*rbac
                degsdr(1)                              = engygs0*engygs0*a1+(engygs0+engygs0)*b1+c1
                if( degsdr(1) .ne. 0.e0_wp ) degsdr(1) = degsdr(1)/t
                degsdr(2)                              = engygs0*engygs0*a2+(engygs0+engygs0)*b2+c2
                if( degsdr(2) .ne. 0.e0_wp ) degsdr(2) = degsdr(2)/t
                degsdr(3)                              = engygs0*engygs0*a3+(engygs0+engygs0)*b3+c3
                if( degsdr(3) .ne. 0.e0_wp ) degsdr(3) = degsdr(3)/t
                call rtocart
                call dedcou
            endif
        end subroutine evaluate_kernel

        !> @brief Evaluate ANCVRT for the H3_GEN_PK2_1964 potential.
        subroutine ancvrt
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter  :: n3atom     = 75
            integer, parameter  :: natom      = 25
            integer, parameter  :: isurf      = 5
            integer, parameter  :: jsurf      = isurf*(isurf+1)/2
            real(wp), parameter :: pi         = 3.141592653589793e0_wp
            real(wp), parameter :: clight     = 2.99792458e08_wp
            real(wp), parameter :: cmu_0      = 4.0e0_wp*pi*1.0e-07_wp
            real(wp), parameter :: cepsilon_0 = 1.0e0_wp/(cmu_0*clight**2)
            real(wp), parameter :: ce         = 1.602176462e-19_wp
            real(wp), parameter :: cplanck    = 6.62606876e-34_wp
            real(wp), parameter :: cm_e       = 9.10938188e-31_wp
            real(wp), parameter :: cang       = 1.0e-10_wp
            real(wp), parameter :: cavogadro  = 6.02214199e23_wp
            real(wp), parameter :: ckcal      = 4.184e10_wp
            real(wp), parameter :: htomillh   = 1000.e0_wp
            real(wp), parameter :: htoev      = 27.2113834e0_wp
            real(wp), parameter :: htokcal    = 627.509470e0_wp
            real(wp), parameter :: htowave    = 219474.631e0_wp
            real(wp), parameter :: htokj      = 2625.49962e0_wp
            real(wp), parameter :: bohr_a     = .5291772083e0_wp
            integer             :: i, ianum(7,32), inc1, inc2, ind, ind2, inu, isave(natom), isum, j, jsave(natom)
            character(len=3)    :: periodic_1(7,32)
            character(len=2)    :: name1(natom)
            character(len=2)    :: name2(natom)
            character(len=1)    :: iblank
            character(len=20)   :: distance
            character(len=20)   :: units
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
            do i=1, natoms
                isave(i) = 0
                jsave(i) = 0
                name1(i) = '  '
                name2(i) = '  '
            enddo
            iblank = ' '
            do ind=1, natoms
                do i=1, 7
                    do j=1, 32
                        if( indexes(ind) .eq. ianum(i,j) )then
                            isave(ind) = i
                            jsave(ind) = j
                        endif
                    enddo
                enddo
            enddo
            do ind=1, natoms
                ind2 = nulbl(ind)
                if( ind2 .eq. 0 ) ind2=ind
            enddo
            inc1 = 0
            do ind=1, irctnt-1
                inc1        = inc1+1
                name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
            enddo
            inc2 = 0
            do ind=irctnt, natoms
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
            cnvrtd  = 1.e0_wp
            cnvrte  = 1.e0_wp
            cnvrtde = 1.e0_wp
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

        !> @brief Evaluate CARTOU for the H3_GEN_PK2_1964 potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            use chempotpy_kinds,              only: wp
            implicit none
            integer            :: i, j
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            call convert_cartesian_input(cart, nulbl, natoms, cnvrtd, ireorder == 1, cartnu)
        end subroutine cartou

        !> @brief Evaluate CARTTOR for the H3_GEN_PK2_1964 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            use chempotpy_kinds,              only: wp
            implicit none
            integer            :: coordinate_status
            real(wp)           :: cmm, cphi, flm, hhd, hym, phi, pi, q, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
            integer            :: i, ind, k, l
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            call cartesian_to_potlib_coordinates(cartnu, natoms, icartr, r, coordinate_status)
            if( coordinate_status /= 0 )then
                write(nflag(18), '(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') icartr
                stop
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the H3_GEN_PK2_1964 potential.
        subroutine eunitzero
            use chempotpy_kinds, only: wp
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

        !> @brief Evaluate RTOCART for the H3_GEN_PK2_1964 potential.
        subroutine rtocart
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            real(wp)           :: epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, termy, &
                & termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes(n3atom, isurf), ygs(n3atom), &
                & yij(n3atom, jsurf), za, zb, zc, zd
            integer :: i, ind, j, j1, j2, k, m1
            if( icartr .eq. 1 )then
                do i = 1, natoms
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
                do i = 1, natoms
                    dgscartnu(i,1) = 0.e0_wp
                    dgscartnu(i,2) = 0.e0_wp
                    dgscartnu(i,3) = 0.e0_wp
                    if( ksdiag .ne. 0 )then
                        do j1=ksdiag, kediag
                            descartnu(i,1,j1) = 0.e0_wp
                            descartnu(i,2,j1) = 0.e0_wp
                            descartnu(i,3,j1) = 0.e0_wp
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do j2=ksoffd, keoffd
                            dijcartnu(i,1,j2) = 0.e0_wp
                            dijcartnu(i,2,j2) = 0.e0_wp
                            dijcartnu(i,3,j2) = 0.e0_wp
                        enddo
                    endif
                    do j = 1, natoms
                        if( j .lt. i )then
                            m1 = natoms*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
                            m1 = natoms*(i-1) - (i*(i-1))/2 + j-i
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
                do i = 1, natoms
                    ygs(i) = degsdr(i)/r(i)
                    if( ksdiag .ne. 0 )then
                        do j=ksdiag, kediag
                            yes(i, j) = deesdr(i,j)/r(i)
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do k=ksoffd, keoffd
                            yij(i, k) = deijdr(i,k)/r(i)
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
                            descartnu(1,k,j1) = term12*yes(1, j1) + term13*yes(3, j1)
                            descartnu(2,k,j1) = -term12*yes(1, j1) + term23*yes(2, j1)
                            descartnu(3,k,j1) = -term13*yes(3, j1) - term23*yes(2, j1)
                        enddo
                    endif
                    if( ksoffd .ne. 0 )then
                        do j2=ksoffd, keoffd
                            dijcartnu(1,k,j2) = term12*yij(1, j2) + term13*yij(3, j2)
                            dijcartnu(2,k,j2) = -term12*yij(1, j2) + term23*yij(2, j2)
                            dijcartnu(3,k,j2) = -term13*yij(3, j2) - term23*yij(2, j2)
                        enddo
                    endif
                enddo
            elseif( icartr .eq. 4 )then
                wh   = 1.007825e0_wp
                wf   = 18.99840e0_wp
                sum  = wh+wf
                eps  = wf/sum
                epsp = wh/sum
                u1   = cos(theta1)
                u2   = cos(theta2)
                u3   = cos(phi)
                ss1  = sin(theta1)
                ss2  = sin(theta2)
                ss3  = sin(phi)
                ya   = 0.0e0_wp
                yb   = 0.0e0_wp
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

        !> @brief Evaluate DEDCOU for the H3_GEN_PK2_1964 potential.
        subroutine dedcou
            use chempotpy_kinds, only: wp
            implicit none
            integer            :: i, j, k
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            if( ireorder .eq. 1 )then
                do i = 1, natoms
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
                do i = 1, natoms
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
end module chempotpy_h3_gen_pk2_1964_pes

!> @brief Expose H3_GEN_PK2_1964 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_h3_gen_pk2_1964_pes, only: evaluate
    use chempotpy_kinds,               only: wp
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
