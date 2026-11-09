!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_clh2_leps_lb_1980_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 3
        integer, parameter            :: natom   = 25
        integer, parameter            :: isurf   = 5
        integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter            :: n3atom  = 75
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: asycm_easyab
        double precision              :: asycm_easyac
        double precision              :: asycm_easybc
        double precision              :: infocm_cartnu(natom,3)
        integer                       :: infocm_icartr
        integer                       :: infocm_indexes(natom)
        integer                       :: infocm_irctnt
        integer                       :: infocm_mder
        integer                       :: infocm_msurf
        integer                       :: infocm_natoms
        character                     :: infocm_ref(5)
        double precision              :: potcm_ahh
        double precision              :: potcm_ahx
        double precision              :: potcm_alfw
        double precision              :: potcm_bhh
        double precision              :: potcm_bhx
        double precision              :: potcm_dhh
        double precision              :: potcm_dhx
        double precision              :: potcm_eps
        double precision              :: potcm_eta1p
        double precision              :: potcm_eta3p
        double precision              :: potcm_eta3s
        double precision              :: potcm_etahh
        double precision              :: potcm_g
        double precision              :: potcm_one
        double precision              :: potcm_rhh
        double precision              :: potcm_rhx
        double precision              :: potcm_xih
        double precision              :: potcm_xix
        double precision              :: precm_c3
        double precision              :: precm_denom
        double precision              :: precm_dh1(10)
        double precision              :: precm_dh2(10)
        double precision              :: precm_dh3(10)
        double precision              :: precm_e(4)
        double precision              :: precm_h(10)
        double precision              :: precm_rt34
        double precision              :: precm_rt38
        double precision              :: precm_scr1(4)
        double precision              :: precm_scr2(4)
        double precision              :: precm_u(4,4)
        double precision              :: precm_xib
        double precision              :: pt1cm_degsdr(n3atom)
        double precision              :: pt1cm_engygs
        double precision              :: pt1cm_r(n3atom)
        double precision              :: pt3cm_ezero(isurf+1)
        double precision              :: pt4cm_deesdr(n3atom,isurf)
        double precision              :: pt4cm_engyes(isurf)
        double precision              :: pt5cm_deijdr(n3atom,jsurf)
        double precision              :: pt5cm_engyij(jsurf)
        double precision              :: usricm_anuzero
        double precision              :: usricm_cart(natom,3)
        integer                       :: usricm_nasurf(isurf+1,isurf+1)
        integer                       :: usricm_nder
        integer                       :: usricm_nflag(20)
        integer                       :: usricm_nulbl(natom)
        double precision              :: usrocm_descart(natom,3,isurf)
        double precision              :: usrocm_dgscart(natom,3)
        double precision              :: usrocm_dijcart(natom,3,jsurf)
        double precision              :: usrocm_pengyes(isurf)
        double precision              :: usrocm_pengygs
        double precision              :: usrocm_pengyij(jsurf)
        double precision              :: utilcm_cnvrtd
        double precision              :: utilcm_cnvrtde
        double precision              :: utilcm_cnvrte
        double precision              :: utilcm_descartnu(natom,3,isurf)
        double precision              :: utilcm_dgscartnu(natom,3)
        double precision              :: utilcm_dijcartnu(natom,3,jsurf)
        integer                       :: utilcm_ireorder
        integer                       :: utilcm_kediag
        integer                       :: utilcm_keoffd
        integer                       :: utilcm_ksdiag
        integer                       :: utilcm_ksoffd
        double precision              :: v
        integer                       :: iatom, idir, istate
        !initialize
        asycm_easyab     = 0.d0; asycm_easyac = 0.d0; asycm_easybc = 0.d0; infocm_cartnu = 0.d0; infocm_icartr = 0
        infocm_indexes   = 0; infocm_irctnt = 0; infocm_mder = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref       = ' '
        potcm_ahh        = 0.d0
        potcm_ahx        = 0.d0
        potcm_alfw       = 0.d0
        potcm_bhh        = 0.d0
        potcm_bhx        = 0.d0
        potcm_dhh        = 0.d0
        potcm_dhx        = 0.d0
        potcm_eps        = 0.d0
        potcm_eta1p      = 0.d0
        potcm_eta3p      = 0.d0
        potcm_eta3s      = 0.d0
        potcm_etahh      = 0.d0
        potcm_g          = 0.d0
        potcm_one        = 0.d0
        potcm_rhh        = 0.d0
        potcm_rhx        = 0.d0
        potcm_xih        = 0.d0
        potcm_xix        = 0.d0
        precm_c3         = 0.d0
        precm_denom      = 0.d0
        precm_dh1        = 0.d0
        precm_dh2        = 0.d0
        precm_dh3        = 0.d0
        precm_e          = 0.d0
        precm_h          = 0.d0
        precm_rt34       = 0.d0
        precm_rt38       = 0.d0
        precm_scr1       = 0.d0
        precm_scr2       = 0.d0
        precm_u          = 0.d0
        precm_xib        = 0.d0
        pt1cm_degsdr     = 0.d0
        pt1cm_engygs     = 0.d0
        pt1cm_r          = 0.d0
        pt3cm_ezero      = 0.d0
        pt4cm_deesdr     = 0.d0
        pt4cm_engyes     = 0.d0
        pt5cm_deijdr     = 0.d0
        pt5cm_engyij     = 0.d0
        usricm_anuzero   = 0.d0
        usricm_cart      = 0.d0
        usricm_nasurf    = 0
        usricm_nder      = 0
        usricm_nflag     = 0
        usricm_nulbl     = 0
        usrocm_descart   = 0.d0
        usrocm_dgscart   = 0.d0
        usrocm_dijcart   = 0.d0
        usrocm_pengyes   = 0.d0
        usrocm_pengygs   = 0.d0
        usrocm_pengyij   = 0.d0
        utilcm_cnvrtd    = 0.d0
        utilcm_cnvrtde   = 0.d0
        utilcm_cnvrte    = 0.d0
        utilcm_descartnu = 0.d0
        utilcm_dgscartnu = 0.d0
        utilcm_dijcartnu = 0.d0
        utilcm_ireorder  = 0
        utilcm_kediag    = 0
        utilcm_keoffd    = 0
        utilcm_ksdiag    = 0
        utilcm_ksoffd    = 0
        usricm_nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0], shape(usricm_nasurf))
        usricm_nder    = 0
        usricm_nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero = 0.0d0
        infocm_icartr  = 3
        infocm_msurf   = 0
        infocm_mder    = 1
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 3
        potcm_eps      = 1.d-6
        potcm_one      = 1.0d0
        potcm_rhh      = 1.4016d0
        potcm_ahh      = 1.0291d0
        potcm_dhh      = 0.17445d0
        potcm_bhh      = 0.018d0
        potcm_xih      = 1.d0
        potcm_xix      = 2.033d0
        potcm_g        = 0.33d0
        potcm_alfw     = 1.0d0
        potcm_rhx      = 2.409d0
        potcm_ahx      = 0.986d0
        potcm_dhx      = 0.1696d0
        potcm_bhx      = 0.013d0
        potcm_etahh    = 0.393764d0
        potcm_eta3s    = 0.321d0
        potcm_eta1p    = 0.20286d0
        potcm_eta3p    = 0.1771d0
        v              = 0.d0; g = 0.d0; d = 0.d0
        usricm_cart    = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                usricm_cart(iatom,idir) = x(iatom,idir)/legacy_bohr_to_angstrom
            enddo
        enddo
        usricm_nder = igrad
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*legacy_hartree_to_electronvolt
            enddo
        else if( igrad == 1 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*legacy_hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = usrocm_dgscart(iatom,idir)*51.422067
                enddo
            enddo
        else if( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif
        contains
        !> @brief Evaluate prepare_local_context for the ClH2_LEPS_LB_1980 potential.
        subroutine prepare_local_context
            !
            !   System:          ClH2
            !   Functional form: Diatomics-in-molecules plus three-center term
            !   Common name:     ClH2 DIM-3C
            !   Reference:       I. Last and M. Baer
            !                    Chem. Phys. Lett. 73, 514 (1980)
            !                    J. Chem. Phys. 75, 288 (1981)
            !
            !   Calling Sequence:
            !      prepare_local_context - initializes the potential's variables and
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
            !      Set in data statements.
            !
            !   Coordinates:
            !      Internal, Definition: R(1) = R(Cl-precm_h)
            !                            R(2) = R(precm_h-precm_h)
            !                            R(3) = R(Cl-precm_h)
            !
            !   Common Blocks (used between the calling program and this potential):
            !        passes the coordinates, ground state electronic energy, and
            !        derivatives of the ground electronic state energy with respect
            !        to the coordinates.
            !        passes the control flags where
            !        usricm_nder  = 0 => no derivatives are computed
            !              = 1 => derivatives of the energy for the ground electroni
            !                     state with respect to the coordinates are computed
            !        usricm_nflag  - Control flags
            !      usricm_nflag(18-20)
            !        passes the FORTRAN unit number used for potential output
            !      /ASYCM/ asycm_easyab, asycm_easybc, asycm_easyac
            !        passes the energy in the three asymptotic valleys for an A + BC
            !        The energy in the AB valley, asycm_easyab, is equal to the energy of
            !        C atom "infinitely" far from the AB diatomic and R(AB) set equa
            !        Re(AB), the equilibrium bond length for the AB diatomic.
            !        In this potential the AB valley represents precm_h infinitely far fro
            !        the ClH diatomic and R(ClH) equal to Re(ClH).  Similarly, the t
            !        asycm_easybc and asycm_easyac represent the energies in the H2 and the othe
            !        valleys, respectively.
            !
            !   Default Parameter Values:
            !      Variable      Default value
            !      usricm_nder             1
            !      usricm_nflag(18)        6
            !
            !
            implicit none
            double precision :: t
            integer          :: i
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            !
            !      LOGICAL LCOL
            !      DIMENSION precm_h(10),precm_dh1(10),precm_dh2(10),precm_dh3(10),precm_u(4,4),precm_e(4),precm_scr1(4),
            !     1           precm_scr2(4)
            !      COMMON /POTCOM/ potcm_eps,potcm_one,potcm_rhh,potcm_ahh,potcm_dhh,potcm_bhh,potcm_xih,potcm_xix,potcm_g,potcm_alfw,potcm_rhx,potcm_ahx,
            !     &                potcm_dhx,potcm_bhx,potcm_etahh,potcm_eta3s,potcm_eta1p,potcm_eta3p
            !
            !   Echo potential parameters to UNIT usricm_nflag(18)
            !
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !
            !WRITE (usricm_nflag(18), 600) potcm_rhh,potcm_ahh,potcm_dhh,potcm_bhh,potcm_rhx,potcm_ahx,potcm_dhx,potcm_bhx
            !WRITE (usricm_nflag(18), 610) potcm_etahh,potcm_eta3s,potcm_eta1p,potcm_eta3p,
            ! *                    potcm_xih,potcm_xix,potcm_g,potcm_alfw
            !
            asycm_easyab = potcm_dhx
            asycm_easybc = potcm_dhh
            asycm_easyac = potcm_dhx
            !
            precm_xib   = 0.5d0*(potcm_xih + potcm_xix)
            t           = precm_xib*potcm_rhx
            precm_c3    = 1.d0/3.d0
            precm_denom = 2.d0*potcm_rhx*(1.d0+t*(1.d0+precm_c3*t))*exp(-t)
            potcm_g     = potcm_g/potcm_rhx
            precm_rt34  = 0.25d0*sqrt(3.d0)
            precm_rt38  = 0.5d0*precm_rt34
            !
            !
            !
            pt3cm_ezero(1) = potcm_dhh
            !
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !
            infocm_ref(1) = 'I. Last and M. Baer'
            infocm_ref(2) = 'Chem. Phys. Lett. 73, 514(1980)'
            infocm_ref(3) = 'J. Chem. Phys. 75, 288 (1981)'
            !
            infocm_indexes(1) = 17
            infocm_indexes(2) = 1
            infocm_indexes(3) = 1
            !
            !
            !
            infocm_irctnt = 2
            !
            !CALL POTINFO
            !
            call ancvrt
            !
        end subroutine prepare_local_context
        !> @brief Evaluate POT for the ClH2_LEPS_LB_1980 potential.
        subroutine pot
            !
            !   Zero of energy:
            !      The classical potential energy is set equal to zero for the Cl
            !      infinitely far from the H2 diatomic and pt1cm_r(H2) set equal to the
            !      H2 equilibrium diatomic value.
            !
            !   Coordinates:
            !      Internal, Definition: pt1cm_r(1) = pt1cm_r(Cl-precm_h)
            !                            pt1cm_r(2) = pt1cm_r(precm_h-precm_h)
            !                            pt1cm_r(3) = pt1cm_r(Cl-precm_h)
            !
            !        The energy in the AB valley, asycm_easyab, is equal to the energy of
            !        C atom "infinitely" far from the AB diatomic and pt1cm_r(AB) set equa
            !        Re(AB), the equilibrium bond length for the AB diatomic.
            !        In this potential the AB valley represents precm_h infinitely far fro
            !        the ClH diatomic and pt1cm_r(ClH) equal to Re(ClH).  Similarly, the t
            !        asycm_easybc and asycm_easyac represent the energies in the H2 and the othe
            !        valleys, respectively.
            !
            implicit none
            double precision :: csa, csa2, csb, csb2, csg2, d1, d2, d3, dcsa21, dcsa22, dcsa23, dcsb21, dcsb22, dcsb23, dcsg21, dcsg22, dcsg23, dp11, dp12, dp21, dp22, dp31, dp32, ds11, ds12, ds21, ds22, ds31, ds32, dshh, dshx1, dshx2, dsn2a1, dsn2a2, dsn2a3, dsn2b1, dsn2b2, dsn2b3, dv1hh, dv1p1, dv1p2, dv1s1, dv1s2, dv3hh, dv3p1, dv3p2, dv3s1, dv3s2, ex1, ex2, p11, p12, p21, p22, p31, p32, r1, r12, r13, r1s, r2, r23, r2s, r3, r3s, rdif, rdif2, s11, s12, s21, s22, s31, s32, shh, shx1, shx2, sn2a, sn2b, sna, sna2, snb, snb2, t, t1, t11, t12, t2, t21, t22, t3, t31, t32, v1hh, v1p1, v1p2, v1s1, v1s2, v3hh, v3p1, v3p2, v3s1, v3s2
            integer          :: i, ierr, ind, index, j, jj, jm, k, kk, l, ll
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            logical :: lcol, matrix_is_diagonal
            !
            call cartou
            call carttor
            !
            !   The potential routines define R1=pt1cm_r(Cl-precm_h), R2=pt1cm_r(precm_h-Cl), and R3=pt1cm_r(precm_h-precm_h);
            !   rearrange the input coordinates to match these definitions
            !
            r1 = pt1cm_r(1)
            r2 = pt1cm_r(3)
            r3 = pt1cm_r(2)
            !
            if( r1 .gt. r2+r3 ) r1 = r2+r3
            if( r2 .gt. r1+r3 ) r2 = r1+r3
            if( r3 .gt. r1+r2 ) r3 = r1+r2
            r1s = r1*r1
            r2s = r2*r2
            r3s = r3*r3
            r12 = r1*r2
            r13 = r1*r3
            r23 = r2*r3
            t3  = (r1s+r2s-r3s)*0.5d0
            t2  = (r1s+r3s-r2s)*0.5d0
            t1  = (r2s+r3s-r1s)*0.5d0
            csa = t2/r13
            if( abs(csa) .gt. potcm_one ) csa = sign(potcm_one, csa)
            csa2 = csa*csa
            sna2 = 1.d0-csa2
            lcol = sna2 .lt. potcm_eps
            sna  = sqrt(sna2)
            t22  = 2.d0*csa
            t11  = 0.d0
            if( .not.lcol ) t11 = 2.d0*(sna2-csa2)/sna
            sn2a   = t22*sna
            t      = t3/(r1*r13)
            dcsa21 = t22*t
            dsn2a1 = t11*t
            t      = -r2/r13
            dcsa22 = t22*t
            dsn2a2 = t11*t
            t      = t1/(r13*r3)
            dcsa23 = t22*t
            dsn2a3 = t11*t
            csb    = t1/r23
            if( abs(csb) .gt. potcm_one ) csb = sign(potcm_one, csb)
            csb2 = csb*csb
            snb2 = 1.d0-csb2
            snb  = sqrt(snb2)
            t11  = 0.d0
            if( .not.lcol ) t11 = 2.d0*(snb2-csb2)/snb
            t22    = 2.d0*csb
            sn2b   = t22*snb
            t      = -r1/r23
            dcsb21 = t22*t
            dsn2b1 = t11*t
            t      = t3/(r2*r23)
            dcsb22 = t22*t
            dsn2b2 = t11*t
            t      = t2/(r23*r3)
            dcsb23 = t22*t
            dsn2b3 = t11*t
            t      = t3/r12
            csg2   = t*t
            t      = 2.d0*t
            dcsg21 = t*t2/(r1*r12)
            dcsg22 = t*t1/(r12*r2)
            dcsg23 = -t*r3/r12
            !
            !  DIATOMIC CURVES
            !
            !    HH
            !
            rdif  = r3 - potcm_rhh
            ex1   = exp(-potcm_ahh*rdif)
            rdif2 = rdif*rdif
            ex2   = exp(-potcm_bhh*rdif*rdif2)
            t1    = potcm_dhh*ex1*ex2
            v1hh  = t1*(ex1-2.d0)
            v3hh  = potcm_etahh*t1*(ex1+2.d0)
            t1    = 2.d0*potcm_ahh*t1
            t2    = 3.d0*potcm_bhh*rdif2
            dv1hh = t1*(1.d0-ex1) - t2*v1hh
            dv3hh = -t1*potcm_etahh*(1.d0+ex1) - t2*v3hh
            !
            !    HX
            !
            call vhx(r1, v1s1, v3s1, v1p1, v3p1, dv1s1, dv3s1, dv1p1, dv3p1)
            call vhx(r2, v1s2, v3s2, v1p2, v3p2, dv1s2, dv3s2, dv1p2, dv3p2)
            !
            s11 = v1s1 + 3.d0*v3s1
            s12 = v1s2 + 3.d0*v3s2
            s21 = 3.d0*v1s1 + v3s1
            s22 = 3.d0*v1s2 + v3s2
            s31 = v1s1 - v3s1
            s32 = v1s2 - v3s2
            p11 = v1p1 + 3.d0*v3p1
            p12 = v1p2 + 3.d0*v3p2
            p21 = 3.d0*v1p1 + v3p1
            p22 = 3.d0*v1p2 + v3p2
            p31 = v1p1 - v3p1
            p32 = v1p2 - v3p2
            !
            ds11 = dv1s1 + 3.d0*dv3s1
            ds12 = dv1s2 + 3.d0*dv3s2
            ds21 = 3.d0*dv1s1 + dv3s1
            ds22 = 3.d0*dv1s2 + dv3s2
            ds31 = dv1s1 - dv3s1
            ds32 = dv1s2 - dv3s2
            dp11 = dv1p1 + 3.d0*dv3p1
            dp12 = dv1p2 + 3.d0*dv3p2
            dp21 = 3.d0*dv1p1 + dv3p1
            dp22 = 3.d0*dv1p2 + dv3p2
            dp31 = dv1p1 - dv3p1
            dp32 = dv1p2 - dv3p2
            !
            !    CONSTRUCT 4X4 HAMILTONIAN PUT IN PACKED ARRAY
            !
            precm_h(1)                 = v1hh + 0.25d0*(s11*csa2 + p11*sna2 + s12*csb2 + p12*snb2)
            precm_h(3)                 = v1hh + 0.25d0*(s11*sna2 + p11*csa2 + s12*snb2 + p12*csb2)
            precm_h(6)                 = v3hh + 0.25d0*(s21*csa2 + p21*sna2 + s22*csb2 + p22*snb2)
            precm_h(10)                = v3hh + 0.25d0*(s21*sna2 + p21*csa2 + s22*snb2 + p22*csb2)
            t11                        = s11 - p11
            t12                        = s12 - p12
            precm_h(2)                 = 0.d0
            if( .not.lcol ) precm_h(2) = 0.125d0*(t11*sn2a - t12*sn2b)
            precm_h(4)                 = precm_rt34*(s31*csa2 + p31*sna2 - s32*csb2 - p32*snb2)
            t31                        = s31 - p31
            t32                        = s32 - p32
            precm_h(5)                 = 0.d0
            if( .not.lcol ) precm_h(5) = precm_rt38*(t31*sn2a + t32*sn2b)
            precm_h(7)                 = precm_h(5)
            precm_h(8)                 = precm_rt34*(s31*sna2 + p31*csa2 - s32*snb2 - p32*csb2)
            t21                        = s21 - p21
            t22                        = s22 - p22
            precm_h(9)                 = 0.d0
            if( .not.lcol ) precm_h(9) = 0.125d0*(t21*sn2a - t22*sn2b)
            !
            !   CONSTRUCT DERIVATIVE OF HAMILTONIAN MATRIX
            !
            t             = t11*dcsa21 + t12*dcsb21
            precm_dh1(1)  = 0.25d0*(ds11*csa2 + dp11*sna2 + t)
            precm_dh1(2)  = 0.125d0*((ds11-dp11)*sn2a + t11*dsn2a1 - t12*dsn2b1)
            precm_dh1(3)  = 0.25d0*(ds11*sna2 + dp11*csa2 - t)
            t             = t11*dcsa22 + t12*dcsb22
            precm_dh2(1)  = 0.25d0*(ds12*csb2 + dp12*snb2 + t)
            precm_dh2(2)  = 0.125d0*(-(ds12-dp12)*sn2b + t11*dsn2a2 - t12*dsn2b2)
            precm_dh2(3)  = 0.25d0*(ds12*snb2 + dp12*csb2 - t)
            t             = 0.25d0*(t11*dcsa23 + t12*dcsb23)
            precm_dh3(1)  = dv1hh + t
            precm_dh3(2)  = 0.125d0*(t11*dsn2a3 - t12*dsn2b3)
            precm_dh3(3)  = dv1hh - t
            t             = t21*dcsa21 + t22*dcsb21
            precm_dh1(6)  = 0.25d0*(ds21*csa2 + dp21*sna2 + t)
            precm_dh1(9)  = 0.125d0*((ds21-dp21)*sn2a + t21*dsn2a1 - t22*dsn2b1)
            precm_dh1(10) = 0.25d0*(ds21*sna2 + dp21*csa2 - t)
            t             = t21*dcsa22 + t22*dcsb22
            precm_dh2(6)  = 0.25d0*(ds22*csb2 + dp22*snb2 + t)
            precm_dh2(9)  = 0.125d0*((dp22-ds22)*sn2b + t21*dsn2a2 - t22*dsn2b2)
            precm_dh2(10) = 0.25d0*(ds22*snb2 + dp22*csb2 - t)
            t             = 0.25d0*(t21*dcsa23 + t22*dcsb23)
            precm_dh3(6)  = dv3hh + t
            precm_dh3(9)  = 0.125d0*(t21*dsn2a3 - t22*dsn2b3)
            precm_dh3(10) = dv3hh - t
            t             = t31*dcsa21 - t32*dcsb21
            precm_dh1(4)  = precm_rt34*(ds31*csa2 + dp31*sna2 + t)
            precm_dh1(5)  = precm_rt38*((ds31-dp31)*sn2a + t31*dsn2a1 + t32*dsn2b1)
            precm_dh1(7)  = precm_dh1(5)
            precm_dh1(8)  = precm_rt34*(ds31*sna2 + dp31*csa2 - t)
            t             = t31*dcsa22 - t32*dcsb22
            precm_dh2(4)  = precm_rt34*(-ds32*csb2 - dp32*snb2 + t)
            precm_dh2(5)  = precm_rt38*((ds32-dp32)*sn2b + t31*dsn2a2 + t32*dsn2b2)
            precm_dh2(7)  = precm_dh2(5)
            precm_dh2(8)  = precm_rt34*(-ds32*snb2 - dp32*csb2 - t)
            t             = t31*dcsa23 - t32*dcsb23
            precm_dh3(4)  = precm_rt34*t
            precm_dh3(5)  = precm_rt38*(t31*dsn2a3 + t32*dsn2b3)
            precm_dh3(7)  = precm_dh3(5)
            precm_dh3(8)  = -precm_dh3(4)
            !
            !  CHECK IF precm_h IS ALREADY DIAGONAL
            !
            matrix_is_diagonal = .true.
            do j = 2, 4
                jm = j - 1
                do i = 1, jm
                    index = (j*(j-1))/2 + i
                    if( abs(precm_h(index)) > 1.0d-30 )then
                        matrix_is_diagonal = .false.
                        exit
                    endif
                enddo
                if( .not.matrix_is_diagonal ) exit
            enddo
            if( matrix_is_diagonal )then
                precm_e(1) = precm_h(1)
                do j = 2, 4
                    jj         = (j*(j+1))/2
                    precm_e(1) = min(precm_e(1), precm_h(jj))
                enddo
                d1 = precm_dh1(1)
                d2 = precm_dh2(1)
                d3 = precm_dh3(1)
            else
                call rsp(4, 4, 10, precm_h, precm_e, 1, precm_u, precm_scr1, precm_scr2, ierr)
                d1 = 0.0d0
                d2 = 0.0d0
                d3 = 0.0d0
                do l = 1, 4
                    t1 = precm_u(l,1)
                    ll = (l*(l-1))/2
                    do k = 1, 4
                        t2 = t1*precm_u(k,1)
                        kk = (k*(k-1))/2
                        if( l >= k ) ind = ll+k
                        if( k > l ) ind  = kk+l
                        d1 = d1+t2*precm_dh1(ind)
                        d2 = d2+t2*precm_dh2(ind)
                        d3 = d3+t2*precm_dh3(ind)
                    enddo
                enddo
            endif
            t            = potcm_xih*r3
            ex1          = exp(-t)
            shh          = (1.d0+t*(1.d0+precm_c3*t))*ex1
            dshh         = -potcm_xih*t*(1.d0+t)*precm_c3*ex1
            t            = precm_xib*r1
            ex1          = exp(-t)
            shx1         = r1*(1.d0+t*(1.d0+precm_c3*t))*ex1/precm_denom
            dshx1        = (1.d0+t*(1.d0-precm_c3*t*t))*ex1/precm_denom
            t            = precm_xib*r2
            ex1          = exp(-t)
            shx2         = r2*(1.d0+t*(1.d0+precm_c3*t))*ex1/precm_denom
            dshx2        = (1.d0+t*(1.d0-precm_c3*t*t))*ex1/precm_denom
            rdif         = r1 - r2
            t1           = potcm_g*exp(-potcm_alfw*rdif*rdif)
            t2           = shh*(shx1+shx2) + shx1*shx2
            pt1cm_engygs = precm_e(1) + t2*t1*csg2 + pt3cm_ezero(1)
            !
            call eunitzero
            if( usricm_nder .ne. 0 )then
                t3              = 2.d0*potcm_alfw*rdif*csg2
                pt1cm_degsdr(1) = d1 + (dshx1*(shh+shx2)*csg2 + t2*(dcsg21-t3))*t1
                pt1cm_degsdr(3) = d2 + (dshx2*(shh+shx1)*csg2 + t2*(dcsg22+t3))*t1
                pt1cm_degsdr(2) = d3 + (dshh*(shx1+shx2)*csg2 + t2*dcsg23)*t1
                call rtocart
                call dedcou
            endif
            !
            !
        end subroutine pot
        !> @brief Evaluate VHX for the ClH2_LEPS_LB_1980 potential.
        subroutine vhx( rr, v1s, v3s, v1p, v3p, dv1s, dv3s, dv1p, dv3p )
            implicit none
            double precision :: dv1, dv1p, dv1s, dv3p, dv3s, ex1, ex2, rdif, rdif2, rr, t1, t2, v1, v1p, v1s, v3p, v3s
            !
            !
            integer, parameter          :: n3atom = 75
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            integer, parameter          :: natom  = 25
            !
            !
            !
            !
            !
            !
            rdif  = rr - potcm_rhx
            rdif2 = rdif*rdif
            ex1   = exp(-potcm_ahx*rdif)
            ex2   = exp(-potcm_bhx*rdif2*rdif)
            t1    = potcm_dhx*ex1*ex2
            v1s   = t1*(ex1-2.d0)
            v1    = t1*(ex1+2.d0)
            v3s   = potcm_eta3s*v1
            v1p   = potcm_eta1p*v1
            v3p   = potcm_eta3p*v1
            t1    = 2.d0*potcm_ahx*t1
            t2    = 3.d0*potcm_bhx*rdif2
            dv1s  = t1*(1.d0-ex1) - t2*v1s
            dv1   = -t1*(1.d0+ex1) - t2*v1
            dv3s  = potcm_eta3s*dv1
            dv1p  = potcm_eta1p*dv1
            dv3p  = potcm_eta3p*dv1
        end subroutine vhx
        !> @brief Evaluate RSP for the ClH2_LEPS_LB_1980 potential.
        !> @brief Diagonalize a real symmetric packed matrix with the shared core eigensolver.
        subroutine rsp( nm, n, nv, a, w, matz, z, fv1, fv2, ierr )
            implicit none
            integer, intent(in)           :: nm, n, nv, matz
            double precision, intent(in)  :: a(nv)
            double precision, intent(out) :: w(n), z(nm,n), fv1(n), fv2(n)
            integer, intent(out)          :: ierr
            integer                       :: i, j, packed_index, lwork
            double precision, allocatable :: matrix(:,:), work(:)
            external                      :: dsyev
            ierr = 0
            z    = 0.0d0
            fv1  = 0.0d0
            fv2  = 0.0d0
            if( n > nm )then
                ierr = 10*n
                return
            endif
            if( nv < n*(n+1)/2 )then
                ierr = 20*n
                return
            endif
            lwork = max(1, 3*n-1)
            allocate(matrix(n,n), work(lwork))
            matrix       = 0.0d0
            packed_index = 0
            do i = 1, n
                do j = 1, i
                    packed_index = packed_index+1
                    matrix(i,j)  = a(packed_index)
                    matrix(j,i)  = a(packed_index)
                enddo
            enddo
            if( matz == 0 )then
                call dsyev('N', 'L', n, matrix, n, w, work, lwork, ierr)
            else
                call dsyev('V', 'L', n, matrix, n, w, work, lwork, ierr)
                if( ierr == 0 ) z(1:n,1:n) = matrix
            endif
            deallocate(matrix, work)
        end subroutine rsp
        !> @brief Evaluate POTINFO for the ClH2_LEPS_LB_1980 potential.
        subroutine potinfo
            implicit none
            integer                     :: i, j, k, kmax
            integer, parameter          :: n3atom = 75
            integer, parameter          :: natom  = 25
            integer, parameter          :: isurf  = 5
            integer, parameter          :: jsurf  = isurf*(isurf+1)/2
            double precision, parameter :: pi     = 3.141592653589793d0
            write(usricm_nflag(18), '(/)')
            do i =1, 5
                write(usricm_nflag(18), '(2x,a75)') infocm_ref(i)
            enddo
            write(usricm_nflag(18), '(/)')
            kmax = 0
            do i = 1, isurf+1
                do j = 1, isurf+1
                    if( usricm_nasurf(i,j) .ne. 0.and.kmax .lt. max(i, j) ) kmax = max(i, j)
                enddo
            enddo
            write(usricm_nflag(18), '(2x,'' MAX. AND ACTUAL NO. OF EXCITED SURFACES: '',I3,5x,I3)') infocm_msurf, kmax-1
            if( kmax-1 .gt. infocm_msurf )then
                write(6, *) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
                stop
            endif
            utilcm_ksdiag = 0
            utilcm_kediag = 0
            do i = 2, isurf+1
                if( usricm_nasurf(i,i) .ne. 0 )then
                    utilcm_kediag = i-1
                    if( utilcm_ksdiag .eq. 0 ) utilcm_ksdiag = i-1
                endif
            enddo
            utilcm_ksoffd = 0
            utilcm_keoffd = 0
            k             = 0
            do i = 1, isurf
                do j = i+1, isurf+1
                    k = k+1
                    if( usricm_nasurf(i,j)+usricm_nasurf(j,i) .ne. 0 )then
                        utilcm_keoffd = k
                        if( utilcm_ksoffd .eq. 0 ) utilcm_ksoffd = k
                    endif
                enddo
            enddo
            write(usricm_nflag(18), '(2x,'' MAX. AND ACTUAL ORDER OF DERIVATIVES:    '',I3,5x,I3)') infocm_mder, usricm_nder
            if( usricm_nder .gt. infocm_mder )then
                write(6, *) ' WRONG INPUT ON ORDER OF DERIVATIVES'
                stop
            endif
            if( usricm_nflag(19) .eq. 1 )then
                write(usricm_nflag(18), '(/)')
                write(usricm_nflag(18), '(2x,''Cartesian coordinates are supplied by'',/, 2x,''the user in the array usricm_cart.'',//)')
                write(usricm_nflag(18), '(2x,''Provide cartesian coordinates in the'',/, 2x,''following order using the array usricm_cart'',' // &
                    '//, 2x,'' usricm_cart(1,1)...usricm_cart(1,3)   => ATOM 1'',/, 2x,'' usricm_cart(2,1)...usricm_cart(2,3)   => AT' // &
                    'OM 2'',/, 2x,'' usricm_cart(3,1)...usricm_cart(3,3)   => ATOM 3'',/, 2x,'' usricm_cart(N,1)...usricm_cart(N,3)  ' // &
                    ' => ATOM N'',/, 2x,''usricm_cart(25,1)...usricm_cart(25,3)  => ATOM 25'',/)')
                write(usricm_nflag(18), '(2x,''If the user wishes to relabel the atoms,'',/, 2x,''set the variable utilcm_ireorder equal to ' // &
                    '1'',/, 2x,''in the PARAMETER statement.  The user'',/, 2x,''must also specify the ne' // &
                    'w labeling'',/, 2x,''scheme.  This is done using the array'',/, 2x,''usricm_nulbl in the fo' // &
                    'llowing manner:'',//, 2x,''usricm_nulbl(i) = j'',/, 2x,''where:  i => old index'',/, 2x,'' ' // &
                    '       j => new index'',//)')
                write(usricm_nflag(18), '(2x,''Cartesian coordinates can be provided to'',/, 2x,''the potential routine in a variety ' // &
                    'of units.'',/, 2x,''The input units will be converted to Bohr'',/, 2x,''based on the' // &
                    ' following values of the usricm_nflag'',/, 2x,''variable:'',//, 2x,''usricm_nflag(1)  =  1  =>  CA' // &
                    'RTESIANS IN BOHR (no'',/, 2x,''                    conversion required)'',/, 2x,''NF' // &
                    'LAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS'',//)')
                write(usricm_nflag(18), '(2x,''The value of the energy and derivatives'',/, 2x,''(if computed) can be reported in a v' // &
                    'ariety'',/, 2x,''units.  A units conversion will take place'',/, 2x,''as specified b' // &
                    'y the following values of the'',/, 2x,''usricm_nflag variable:'',//, 2x,''usricm_nflag(2) = 1 =>  ' // &
                    'ENERGIES REPORTED IN HARTEEE'',/, 2x,''usricm_nflag(2) = 2 =>  ENERGIES REPORTED IN mHARTRE' // &
                    'E'',/, 2x,''usricm_nflag(2) = 3 =>  ENERGIES REPORTED IN eV'',/, 2x,''usricm_nflag(2) = 4 =>  ENER' // &
                    'GIES REPORTED IN kcal/mol'',/, 2x,''usricm_nflag(2) = 5 =>  ENERGIES REPORTED IN cm**-1'',/' // &
                    '/)')
                write(usricm_nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
                    'he'',/, 2x,''usricm_nflag variable:'',//, 2x,''usricm_nflag(1)=1 & usricm_nflag(2)=1 => DERIVATIVES REPOR' // &
                    'TED IN'',/, 2x,''                           HARTEEE/BOHR'',/, 2x,''usricm_nflag(1)=1 & NFLA' // &
                    'G(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/BOHR' // &
                    ''',/, 2x,''usricm_nflag(1)=1 & usricm_nflag(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''            ' // &
                    '               eV/BOHR'',/, 2x,''usricm_nflag(1)=1 & usricm_nflag(2)=4 => DERIVATIVES REPORTED IN' // &
                    ''',/, 2x,''                           kcal/mol/BOHR'',/, 2x,''usricm_nflag(1)=1 & usricm_nflag(2)=' // &
                    '5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm**-1/BOHR'',//)')
                write(usricm_nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
                    'he'',/, 2x,''usricm_nflag variable:'',//, 2x,''usricm_nflag(1)=2 & usricm_nflag(2)=1 => DERIVATIVES REPOR' // &
                    'TED IN'',/, 2x,''                           HARTEEE/ANGSTROM'',/, 2x,''usricm_nflag(1)=2 & ' // &
                    'usricm_nflag(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/' // &
                    'ANGSTROM'',/, 2x,''usricm_nflag(1)=2 & usricm_nflag(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''    ' // &
                    '                       eV/ANGSTROM'',/, 2x,''usricm_nflag(1)=2 & usricm_nflag(2)=4 => DERIVATIVES ' // &
                    'REPORTED IN'',/, 2x,''                           kcal/mol/ANGSTROM'',/, 2x,''usricm_nflag(1' // &
                    ')=2 & usricm_nflag(2)=5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm*' // &
                    '*-1/ANGSTROM'',//)')
            endif
        end subroutine potinfo
        !> @brief Evaluate ANCVRT for the ClH2_LEPS_LB_1980 potential.
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
            do i=1, infocm_natoms
                isave(i) = 0
                jsave(i) = 0
                name1(i) = '  '
                name2(i) = '  '
            enddo
            iblank = ' '
            do ind=1, infocm_natoms
                do i=1, 7
                    do j=1, 32
                        if( infocm_indexes(ind) .eq. ianum(i,j) )then
                            isave(ind) = i
                            jsave(ind) = j
                        endif
                    enddo
                enddo
            enddo
            do ind=1, infocm_natoms
                ind2 = usricm_nulbl(ind)
                if( ind2 .eq. 0 ) ind2=ind
            enddo
            inc1 = 0
            do ind=1, infocm_irctnt-1
                inc1        = inc1+1
                name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
            enddo
            inc2 = 0
            do ind=infocm_irctnt, infocm_natoms
                inc2        = inc2+1
                name2(inc2) = periodic_1(isave(ind),jsave(ind))(:2)
            enddo
            if( usricm_nflag(1) .eq. 2 ) distance = 'ANGSTROMS           '
            if( usricm_nflag(2) .eq. 2 )then
                units = 'MILLIHARTREE        '
            else if( usricm_nflag(2) .eq. 3 )then
                units = 'EV                  '
            else if( usricm_nflag(2) .eq. 4 )then
                units = 'KCAL PER MOLE       '
            else if( usricm_nflag(2) .eq. 5 )then
                units = 'WAVENUMBERS         '
            else if( usricm_nflag(2) .eq. 6 )then
                units = 'KILOJOULES PER MOLE '
            endif
            utilcm_cnvrtd  = 1.d0
            utilcm_cnvrte  = 1.d0
            utilcm_cnvrtde = 1.d0
            if( usricm_nflag(1) .eq. 2 ) utilcm_cnvrtd = bohr_a
            if( usricm_nflag(2) .eq. 2 )then
                utilcm_cnvrte = utilcm_cnvrte*htomillh
            else if( usricm_nflag(2) .eq. 3 )then
                utilcm_cnvrte = utilcm_cnvrte*htoev
            else if( usricm_nflag(2) .eq. 4 )then
                utilcm_cnvrte = utilcm_cnvrte*htokcal
            else if( usricm_nflag(2) .eq. 5 )then
                utilcm_cnvrte = utilcm_cnvrte*htowave
            else if( usricm_nflag(2) .eq. 6 )then
                utilcm_cnvrte = utilcm_cnvrte*htokj
            endif
            utilcm_cnvrtde = utilcm_cnvrte/utilcm_cnvrtd
            isum           = 0
            do inu=1, 25
                isum = isum + usricm_nulbl(inu)
            enddo
            utilcm_ireorder = 0
            if( isum .ne. 0 ) utilcm_ireorder = 1
        end subroutine ancvrt
        !> @brief Evaluate CARTOU for the ClH2_LEPS_LB_1980 potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            implicit none
            integer            :: i, j
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            call convert_cartesian_input(usricm_cart, usricm_nulbl, infocm_natoms, utilcm_cnvrtd, utilcm_ireorder == 1, infocm_cartnu)
        end subroutine cartou
        !> @brief Evaluate CARTTOR for the ClH2_LEPS_LB_1980 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer            :: coordinate_status
            double precision   :: cmm, cphi, flm, hhd, hym, phi, pi, q, q1, q2, q3, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, zcm1, zcm2, zcm3, zrm1, zrm2
            integer            :: i, ind, k, l
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            call cartesian_to_potlib_coordinates(infocm_cartnu, infocm_natoms, infocm_icartr, pt1cm_r, coordinate_status)
            if( coordinate_status /= 0 )then
                write(usricm_nflag(18), '(2x,''WRONG infocm_icartr FOR infocm_cartnu; infocm_icartr ='',i5//)') infocm_icartr
                stop
            endif
        end subroutine carttor
        !> @brief Evaluate EUNITZERO for the ClH2_LEPS_LB_1980 potential.
        subroutine eunitzero
            implicit none
            integer            :: i, j
            character(len=75)  :: ref(5)
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            usrocm_pengygs = pt1cm_engygs * utilcm_cnvrte - usricm_anuzero
            if( utilcm_ksdiag .ne. 0 )then
                do i=utilcm_ksdiag, utilcm_kediag
                    usrocm_pengyes(i) = pt4cm_engyes(i) * utilcm_cnvrte - usricm_anuzero
                enddo
            endif
            if( utilcm_ksoffd .ne. 0 )then
                do j=utilcm_ksoffd, utilcm_keoffd
                    usrocm_pengyij(j) = pt5cm_engyij(j) * utilcm_cnvrte
                enddo
            endif
        end subroutine eunitzero
        !> @brief Evaluate RTOCART for the ClH2_LEPS_LB_1980 potential.
        subroutine rtocart
            implicit none
            double precision   :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes, ygs, yij, za, zb, zc, zd
            integer            :: i, ind, j, j1, j2, k, m1
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            dimension ygs(n3atom), yes(n3atom,isurf), yij(n3atom,jsurf)
            if( infocm_icartr .eq. 1 )then
                do i = 1, infocm_natoms
                    ind                   = 3*i-2
                    utilcm_dgscartnu(i,1) = pt1cm_degsdr(ind)
                    utilcm_dgscartnu(i,2) = pt1cm_degsdr(ind+1)
                    utilcm_dgscartnu(i,3) = pt1cm_degsdr(ind+2)
                    if( utilcm_ksdiag .ne. 0 )then
                        do j = utilcm_ksdiag, utilcm_kediag
                            utilcm_descartnu(i,1,j) = pt4cm_deesdr(ind,j)
                            utilcm_descartnu(i,2,j) = pt4cm_deesdr(ind+1,j)
                            utilcm_descartnu(i,3,j) = pt4cm_deesdr(ind+2,j)
                        enddo
                    endif
                    if( utilcm_keoffd .ne. 0 )then
                        do k = utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(i,1,k) = pt5cm_deijdr(ind,k)
                            utilcm_dijcartnu(i,2,k) = pt5cm_deijdr(ind+1,k)
                            utilcm_dijcartnu(i,3,k) = pt5cm_deijdr(ind+2,k)
                        enddo
                    endif
                enddo
            else if( infocm_icartr .eq. 2 )then
                do i = 1, infocm_natoms
                    utilcm_dgscartnu(i,1) = 0.d0
                    utilcm_dgscartnu(i,2) = 0.d0
                    utilcm_dgscartnu(i,3) = 0.d0
                    if( utilcm_ksdiag .ne. 0 )then
                        do j1=utilcm_ksdiag, utilcm_kediag
                            utilcm_descartnu(i,1,j1) = 0.d0
                            utilcm_descartnu(i,2,j1) = 0.d0
                            utilcm_descartnu(i,3,j1) = 0.d0
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do j2=utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(i,1,j2) = 0.d0
                            utilcm_dijcartnu(i,2,j2) = 0.d0
                            utilcm_dijcartnu(i,3,j2) = 0.d0
                        enddo
                    endif
                    do j = 1, infocm_natoms
                        if( j .lt. i )then
                            m1 = infocm_natoms*(j-1) - (j*(j-1))/2 + i-j
                        else if( j .gt. i )then
                            m1 = infocm_natoms*(i-1) - (i*(i-1))/2 + j-i
                        else
                            cycle
                        endif
                        y                     = pt1cm_degsdr(m1)
                        termx                 = (infocm_cartnu(i,1)-infocm_cartnu(j,1))/pt1cm_r(m1)
                        termy                 = (infocm_cartnu(i,2)-infocm_cartnu(j,2))/pt1cm_r(m1)
                        termz                 = (infocm_cartnu(i,3)-infocm_cartnu(j,3))/pt1cm_r(m1)
                        utilcm_dgscartnu(i,1) = utilcm_dgscartnu(i,1) + termx*y
                        utilcm_dgscartnu(i,2) = utilcm_dgscartnu(i,2) + termy*y
                        utilcm_dgscartnu(i,3) = utilcm_dgscartnu(i,3) + termz*y
                        if( utilcm_ksdiag .gt. 0 )then
                            y = pt4cm_deesdr(m1,j1)
                            do j1=utilcm_ksdiag, utilcm_kediag
                                utilcm_descartnu(i,1,j1) = utilcm_descartnu(i,1,j1) + termx*y
                                utilcm_descartnu(i,2,j1) = utilcm_descartnu(i,2,j1) + termy*y
                                utilcm_descartnu(i,3,j1) = utilcm_descartnu(i,3,j1) + termz*y
                            enddo
                        else if( utilcm_ksoffd .gt. 0 )then
                            do j2=utilcm_ksoffd, utilcm_keoffd
                                y                        = pt5cm_deijdr(m1,j2)
                                utilcm_dijcartnu(i,1,j2) = utilcm_dijcartnu(i,1,j2) + termx*y
                                utilcm_dijcartnu(i,2,j2) = utilcm_dijcartnu(i,2,j2) + termy*y
                                utilcm_dijcartnu(i,3,j2) = utilcm_dijcartnu(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            else if( infocm_icartr .eq. 3 )then
                do i = 1, infocm_natoms
                    ygs(i) = pt1cm_degsdr(i)/pt1cm_r(i)
                    if( utilcm_ksdiag .ne. 0 )then
                        do j=utilcm_ksdiag, utilcm_kediag
                            yes(i,j) = pt4cm_deesdr(i,j)/pt1cm_r(i)
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do k=utilcm_ksoffd, utilcm_keoffd
                            yij(i,k) = pt5cm_deijdr(i,k)/pt1cm_r(i)
                        enddo
                    endif
                enddo
                do k = 1, 3
                    term12                = infocm_cartnu(1,k)-infocm_cartnu(2,k)
                    term23                = infocm_cartnu(2,k)-infocm_cartnu(3,k)
                    term13                = infocm_cartnu(1,k)-infocm_cartnu(3,k)
                    utilcm_dgscartnu(1,k) = term12*ygs(1) + term13*ygs(3)
                    utilcm_dgscartnu(2,k) = -term12*ygs(1) + term23*ygs(2)
                    utilcm_dgscartnu(3,k) = -term13*ygs(3) - term23*ygs(2)
                    if( utilcm_ksdiag .ne. 0 )then
                        do j1=utilcm_ksdiag, utilcm_kediag
                            utilcm_descartnu(1,k,j1) = term12*yes(1,j1) + term13*yes(3,j1)
                            utilcm_descartnu(2,k,j1) = -term12*yes(1,j1) + term23*yes(2,j1)
                            utilcm_descartnu(3,k,j1) = -term13*yes(3,j1) - term23*yes(2,j1)
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do j2=utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(1,k,j2) = term12*yij(1,j2) + term13*yij(3,j2)
                            utilcm_dijcartnu(2,k,j2) = -term12*yij(1,j2) + term23*yij(2,j2)
                            utilcm_dijcartnu(3,k,j2) = -term13*yij(3,j2) - term23*yij(2,j2)
                        enddo
                    endif
                enddo
            else if( infocm_icartr .eq. 4 )then
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
                write(usricm_nflag(18), '(2X,'' WRONG infocm_icartr FOR DERIVATIVE; infocm_icartr ='',I5//)') infocm_icartr
                stop
            endif
        end subroutine rtocart
        !> @brief Evaluate DEDCOU for the ClH2_LEPS_LB_1980 potential.
        subroutine dedcou
            implicit none
            integer            :: i, j, k
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            if( utilcm_ireorder .eq. 1 )then
                do i = 1, infocm_natoms
                    usrocm_dgscart(i,1) = utilcm_dgscartnu(usricm_nulbl(i),1) * utilcm_cnvrtde
                    usrocm_dgscart(i,2) = utilcm_dgscartnu(usricm_nulbl(i),2) * utilcm_cnvrtde
                    usrocm_dgscart(i,3) = utilcm_dgscartnu(usricm_nulbl(i),3) * utilcm_cnvrtde
                    if( utilcm_ksdiag .ne. 0 )then
                        do j=utilcm_ksdiag, utilcm_kediag
                            usrocm_descart(i,1,j) = utilcm_descartnu(usricm_nulbl(i),1,j) * utilcm_cnvrtde
                            usrocm_descart(i,2,j) = utilcm_descartnu(usricm_nulbl(i),2,j) * utilcm_cnvrtde
                            usrocm_descart(i,3,j) = utilcm_descartnu(usricm_nulbl(i),3,j) * utilcm_cnvrtde
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do k=utilcm_ksoffd, utilcm_keoffd
                            usrocm_dijcart(i,1,k) = utilcm_dijcartnu(usricm_nulbl(i),1,k) * utilcm_cnvrtde
                            usrocm_dijcart(i,2,k) = utilcm_dijcartnu(usricm_nulbl(i),2,k) * utilcm_cnvrtde
                            usrocm_dijcart(i,3,k) = utilcm_dijcartnu(usricm_nulbl(i),3,k) * utilcm_cnvrtde
                        enddo
                    endif
                enddo
            else
                do i = 1, infocm_natoms
                    usrocm_dgscart(i,1) = utilcm_dgscartnu(i,1) * utilcm_cnvrtde
                    usrocm_dgscart(i,2) = utilcm_dgscartnu(i,2) * utilcm_cnvrtde
                    usrocm_dgscart(i,3) = utilcm_dgscartnu(i,3) * utilcm_cnvrtde
                    if( utilcm_ksdiag .ne. 0 )then
                        do j=utilcm_ksdiag, utilcm_kediag
                            usrocm_descart(i,1,j) = utilcm_descartnu(i,1,j) * utilcm_cnvrtde
                            usrocm_descart(i,2,j) = utilcm_descartnu(i,2,j) * utilcm_cnvrtde
                            usrocm_descart(i,3,j) = utilcm_descartnu(i,3,j) * utilcm_cnvrtde
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do k=utilcm_ksoffd, utilcm_keoffd
                            usrocm_dijcart(i,1,k) = utilcm_dijcartnu(i,1,k) * utilcm_cnvrtde
                            usrocm_dijcart(i,2,k) = utilcm_dijcartnu(i,2,k) * utilcm_cnvrtde
                            usrocm_dijcart(i,3,k) = utilcm_dijcartnu(i,3,k) * utilcm_cnvrtde
                        enddo
                    endif
                enddo
            endif
        end subroutine dedcou
    end subroutine evaluate
end module chempotpy_clh2_leps_lb_1980_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_clh2_leps_lb_1980_pes, only: evaluate
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
