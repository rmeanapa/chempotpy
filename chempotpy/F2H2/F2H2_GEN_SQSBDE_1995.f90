!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_f2h2_gen_sqsbde_1995_pes
    implicit none
    private
    public :: evaluate
    contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_legacy_units, only: legacy_bohr_to_angstrom
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 4
        integer, parameter            :: natom   = 25
        integer, parameter            :: isurf   = 5
        integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter            :: ma      = 135
        integer, parameter            :: n3atom  = 75
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: ang_g(0:4,0:4,0:2*4,0:4)
        double precision              :: ang_pf(0:4,0:4)
        double precision              :: infocm_cartnu(natom,3)
        integer                       :: infocm_icartr
        integer                       :: infocm_indexes(natom)
        integer                       :: infocm_irctnt
        integer                       :: infocm_mder
        integer                       :: infocm_msurf
        integer                       :: infocm_natoms
        character                     :: infocm_ref(5)
        double precision              :: potcm_a(ma)
        double precision              :: potcm_ac(3,0:4,0:4,0:2*4)
        double precision              :: potcm_ar(0:4,0:4,0:2*4)
        double precision              :: potcm_da(ma)
        double precision              :: potcm_efct(401)
        double precision              :: potcm_fac(0:2*4)
        double precision              :: potcm_fct(401)
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
        ang_g            = 0.d0; ang_pf = 0.d0; infocm_cartnu = 0.d0; infocm_icartr = 0; infocm_indexes = 0; infocm_irctnt = 0
        infocm_mder      = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref       = ' '
        potcm_a          = 0.d0
        potcm_ac         = 0.d0
        potcm_ar         = 0.d0
        potcm_da         = 0.d0
        potcm_efct       = 0.d0
        potcm_fac        = 0.d0
        potcm_fct        = 0.d0
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
        infocm_icartr  = 4
        infocm_msurf   = 0
        infocm_mder    = 0
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 4
        potcm_a(1)     = 1558.0d0
        potcm_da(1)    = 0.0d0
        potcm_a(2)     = -202900000.0d0
        potcm_da(2)    = 0.0d0
        potcm_a(3)     = -3703000000.0d0
        potcm_da(3)    = 0.0d0
        potcm_a(4)     = -11350000000.0d0
        potcm_da(4)    = 0.0d0
        potcm_a(5)     = -1736693.09500000d0
        potcm_da(5)    = 0.0d0
        potcm_a(6)     = -4416487.722000003d0
        potcm_da(6)    = 0.0d0
        potcm_a(7)     = 16012389.5100000d0
        potcm_da(7)    = 0.0d0
        potcm_a(8)     = -6389444.85500000d0
        potcm_da(8)    = 0.0d0
        potcm_a(9)     = -866.867375551329d0
        potcm_da(9)    = 117.589789927560d0
        potcm_a(10)    = 129789.363470169d0
        potcm_da(10)   = 1937.56316874928d0
        potcm_a(11)    = 0.0d0
        potcm_da(11)   = 0.0d0
        potcm_a(12)    = 0.0d0
        potcm_da(12)   = 0.0d0
        potcm_a(13)    = 2622184.75665877d0
        potcm_da(13)   = 29726.6301460444d0
        potcm_a(14)    = 0.0d0
        potcm_da(14)   = 0.0d0
        potcm_a(15)    = 0.0d0
        potcm_da(15)   = 0.0d0
        potcm_a(16)    = 0.0d0
        potcm_da(16)   = 0.0d0
        potcm_a(17)    = 184254701843.311d0
        potcm_da(17)   = 17553092934.4463d0
        potcm_a(18)    = -320000.0d0
        potcm_da(18)   = 0.0d0
        potcm_a(19)    = 3.326438943860246d-2
        potcm_da(19)   = 1.353645657747023d-3
        potcm_a(20)    = -4.576887206887037d-3
        potcm_da(20)   = 1.865017706495822d-4
        potcm_a(21)    = 0.409730666039653d0
        potcm_da(21)   = 2.777020547964313d-2
        potcm_a(22)    = 32623551943283.8d0
        potcm_da(22)   = 3630087056915.05d0
        potcm_a(23)    = 0.0d0
        potcm_da(23)   = 0.0d0
        potcm_a(24)    = 71754497894201.8d0
        potcm_da(24)   = 2525041881368.02d0
        potcm_a(25)    = 0.0d0
        potcm_da(25)   = 0.0d0
        potcm_a(26)    = 5.032387872479119d-2
        potcm_da(26)   = 1.860431838521344d-3
        potcm_a(27)    = 1.46832057915281d0
        potcm_da(27)   = 9.713248031677595d-2
        potcm_a(28)    = 3384.63155928653d0
        potcm_da(28)   = 85.5567822945890d0
        potcm_a(29)    = 0.0d0
        potcm_da(29)   = 0.0d0
        potcm_a(30)    = 0.0d0
        potcm_da(30)   = 0.0d0
        potcm_a(31)    = -160267.101272529d0
        potcm_da(31)   = 4055.33877671770d0
        potcm_a(32)    = 0.0d0
        potcm_da(32)   = 0.0d0
        potcm_a(33)    = 0.0d0
        potcm_da(33)   = 0.0d0
        potcm_a(34)    = 0.0d0
        potcm_da(34)   = 0.0d0
        potcm_a(35)    = 0.0d0
        potcm_da(35)   = 0.0d0
        potcm_a(36)    = -16268.8946598104d0
        potcm_da(36)   = 541.139677344219d0
        potcm_a(37)    = 0.0d0
        potcm_da(37)   = 0.0d0
        potcm_a(38)    = 0.0d0
        potcm_da(38)   = 0.0d0
        potcm_a(39)    = 2802.57610197247d0
        potcm_da(39)   = 56.8038283168248d0
        potcm_a(40)    = 0.0d0
        potcm_da(40)   = 0.0d0
        potcm_a(41)    = 0.0d0
        potcm_da(41)   = 0.0d0
        potcm_a(42)    = 0.0d0
        potcm_da(42)   = 0.0d0
        potcm_a(43)    = 0.0d0
        potcm_da(43)   = 0.0d0
        potcm_a(44)    = 0.0d0
        potcm_da(44)   = 0.0d0
        potcm_a(45)    = 0.0d0
        potcm_da(45)   = 0.0d0
        potcm_a(46)    = 2860.06394339599d0
        potcm_da(46)   = 44.3633384808498d0
        potcm_a(47)    = 0.0d0
        potcm_da(47)   = 0.0d0
        potcm_a(48)    = 0.0d0
        potcm_da(48)   = 0.0d0
        potcm_a(49)    = 0.0d0
        potcm_da(49)   = 0.0d0
        potcm_a(50)    = 0.0d0
        potcm_da(50)   = 0.0d0
        potcm_a(51)    = 0.0d0
        potcm_da(51)   = 0.0d0
        potcm_a(52)    = 509.102022118840d0
        potcm_da(52)   = 18.7375880046217d0
        potcm_a(53)    = 0.0d0
        potcm_da(53)   = 0.0d0
        potcm_a(54)    = 0.0d0
        potcm_da(54)   = 0.0d0
        potcm_a(55)    = 0.0d0
        potcm_da(55)   = 0.0d0
        potcm_a(56)    = 0.0d0
        potcm_da(56)   = 0.0d0
        potcm_a(57)    = 0.0d0
        potcm_da(57)   = 0.0d0
        potcm_a(58)    = 0.0d0
        potcm_da(58)   = 0.0d0
        potcm_a(59)    = 0.0d0
        potcm_da(59)   = 0.0d0
        potcm_a(60)    = 0.0d0
        potcm_da(60)   = 0.0d0
        potcm_a(61)    = 0.0d0
        potcm_da(61)   = 0.0d0
        potcm_a(62)    = 0.0d0
        potcm_da(62)   = 0.0d0
        potcm_a(63)    = 0.0d0
        potcm_da(63)   = 0.0d0
        potcm_a(64)    = 0.0d0
        potcm_da(64)   = 0.0d0
        potcm_a(65)    = 0.0d0
        potcm_da(65)   = 0.0d0
        potcm_a(66)    = 0.0d0
        potcm_da(66)   = 0.0d0
        potcm_a(67)    = 0.0d0
        potcm_da(67)   = 0.0d0
        potcm_a(68)    = -3467.23417017878d0
        potcm_da(68)   = 52.9105695619519d0
        potcm_a(69)    = 0.0d0
        potcm_da(69)   = 0.0d0
        potcm_a(70)    = 0.0d0
        potcm_da(70)   = 0.0d0
        potcm_a(71)    = 0.0d0
        potcm_da(71)   = 0.0d0
        potcm_a(72)    = 0.0d0
        potcm_da(72)   = 0.0d0
        potcm_a(73)    = 2304.6328580383d0
        potcm_da(73)   = 44.0232721532835d0
        potcm_a(74)    = 0.0d0
        potcm_da(74)   = 0.0d0
        potcm_a(75)    = 0.0d0
        potcm_da(75)   = 0.0d0
        potcm_a(76)    = 0.0d0
        potcm_da(76)   = 0.0d0
        potcm_a(77)    = 0.0d0
        potcm_da(77)   = 0.0d0
        potcm_a(78)    = 0.0d0
        potcm_da(78)   = 0.0d0
        potcm_a(79)    = 0.0d0
        potcm_da(79)   = 0.0d0
        potcm_a(80)    = 0.0d0
        potcm_da(80)   = 0.0d0
        potcm_a(81)    = 0.0d0
        potcm_da(81)   = 0.0d0
        potcm_a(82)    = 0.0d0
        potcm_da(82)   = 0.0d0
        potcm_a(83)    = 0.0d0
        potcm_da(83)   = 0.0d0
        potcm_a(84)    = -644.815083021031d0
        potcm_da(84)   = 18.1151965381963d0
        potcm_a(85)    = 0.0d0
        potcm_da(85)   = 0.0d0
        potcm_a(86)    = 0.0d0
        potcm_da(86)   = 0.0d0
        potcm_a(87)    = 0.0d0
        potcm_da(87)   = 0.0d0
        potcm_a(88)    = 0.0d0
        potcm_da(88)   = 0.0d0
        potcm_a(89)    = 0.0d0
        potcm_da(89)   = 0.0d0
        potcm_a(90)    = 0.0d0
        potcm_da(90)   = 0.0d0
        potcm_a(91)    = 0.0d0
        potcm_da(91)   = 0.0d0
        potcm_a(92)    = 0.0d0
        potcm_da(92)   = 0.0d0
        potcm_a(93)    = 0.0d0
        potcm_da(93)   = 0.0d0
        potcm_a(94)    = 0.0d0
        potcm_da(94)   = 0.0d0
        potcm_a(95)    = 0.0d0
        potcm_da(95)   = 0.0d0
        potcm_a(96)    = 0.0d0
        potcm_da(96)   = 0.0d0
        potcm_a(97)    = 0.0d0
        potcm_da(97)   = 0.0d0
        potcm_a(98)    = 0.0d0
        potcm_da(98)   = 0.0d0
        potcm_a(99)    = 0.0d0
        potcm_da(99)   = 0.0d0
        potcm_a(100)   = 0.0d0
        potcm_da(100)  = 0.0d0
        potcm_a(101)   = 22530.9324112547d0
        potcm_da(101)  = 989.768118754483d0
        potcm_a(102)   = 1630.40722103571d0
        potcm_da(102)  = 111.781936229583d0
        potcm_a(103)   = 0.0d0
        potcm_da(103)  = 0.0d0
        potcm_a(104)   = 1042.90771467500d0
        potcm_da(104)  = 18.4325305864792d0
        potcm_a(105)   = -366.368726725276d0
        potcm_da(105)  = 12.8434812694663d0
        potcm_a(106)   = 2428.17244795167d0
        potcm_da(106)  = 246.109466713404d0
        potcm_a(107)   = -1882.37787794149d0
        potcm_da(107)  = 70.8301781102542d0
        potcm_a(108)   = 0.0d0
        potcm_da(108)  = 0.0d0
        potcm_a(109)   = 0.0d0
        potcm_da(109)  = 0.0d0
        potcm_a(110)   = -326.002341194829d0
        potcm_da(110)  = 19.3580274207174d0
        potcm_a(111)   = 0.0d0
        potcm_da(111)  = 0.0d0
        potcm_a(112)   = 0.0d0
        potcm_da(112)  = 0.0d0
        potcm_a(113)   = 116.704136523402d0
        potcm_da(113)  = 11.3452311839242d0
        potcm_a(114)   = 0.0d0
        potcm_da(114)  = 0.0d0
        potcm_a(115)   = 572.632667044664d0
        potcm_da(115)  = 33.6908429103859d0
        potcm_a(116)   = 0.0d0
        potcm_da(116)  = 0.0d0
        potcm_a(117)   = 0.0d0
        potcm_da(117)  = 0.0d0
        potcm_a(118)   = 0.0d0
        potcm_da(118)  = 0.0d0
        potcm_a(119)   = 0.0d0
        potcm_da(119)  = 0.0d0
        potcm_a(120)   = 0.0d0
        potcm_da(120)  = 0.0d0
        potcm_a(121)   = 0.0d0
        potcm_da(121)  = 0.0d0
        potcm_a(122)   = 0.0d0
        potcm_da(122)  = 0.0d0
        potcm_a(123)   = 0.0d0
        potcm_da(123)  = 0.0d0
        potcm_a(124)   = 0.0d0
        potcm_da(124)  = 0.0d0
        potcm_a(125)   = 0.0d0
        potcm_da(125)  = 0.0d0
        potcm_a(126)   = 0.0d0
        potcm_da(126)  = 0.0d0
        potcm_a(127)   = 0.0d0
        potcm_da(127)  = 0.0d0
        potcm_a(128)   = 0.0d0
        potcm_da(128)  = 0.0d0
        potcm_a(129)   = 0.0d0
        potcm_da(129)  = 0.0d0
        potcm_a(130)   = 0.0d0
        potcm_da(130)  = 0.0d0
        potcm_a(131)   = 0.0d0
        potcm_da(131)  = 0.0d0
        potcm_a(132)   = 0.0d0
        potcm_da(132)  = 0.0d0
        potcm_a(133)   = 0.0d0
        potcm_da(133)  = 0.0d0
        potcm_a(134)   = 0.0d0
        potcm_da(134)  = 0.0d0
        potcm_a(135)   = 0.0d0
        potcm_da(135)  = 0.0d0
        v              = 0.d0; g = 0.d0; d = 0.d0
        usricm_cart    = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                usricm_cart(iatom,idir) = x(iatom,idir)/legacy_bohr_to_angstrom
            enddo
        enddo
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs/8065.7112013
            enddo
        else
            write (*, *) 'Only energy is available'
        endif
        contains
        !> @brief Evaluate prepare_local_context for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine prepare_local_context
            !
            !   System:                        F2H2
            !   Common name:                   SQSBDE
            !   Number of derivatives:         0
            !   Number of bodies:              4
            !   Number of electronic surfaces: 1
            !   Interface:                     potlib2001
            !   References:      M. Quack, M. Suhm Mol. Phys. 69 (1990) P791
            !                    Kofranek et al Chem. Phys. 121 (1988) 137
            !                    Rijks and Wormer J. Chem. Phys. 90(1989) 6507
            !
            !   Protocol:
            !
            !      prepare_local_context - initializes the potential's variables
            !               must be called once before any calls to POT
            !      POT    - driver for the evaluation of the energy
            !
            !   Units:
            !
            !      energies    - cm-1
            !      coordinates - distances: bohr; angles: degrees
            !
            !   Surfaces:
            !
            !      ground electronic state
            !
            ! M. A. SUHM  SEPTEMBER 1990
            !
            ! INPUT : SIX INTERNAL COORDINATES
            ! --R1,R2,R IN ATOMIC UNITS (BOHR)
            ! --TH1,TH2,TAU IN DEGREES
            ! DEFINITION AS IN QUACK+SUHM MOL.PHYS. 69 (1990) P791
            ! OUTPUT : HF DIMER (SQSBDE) POTENTIAL ENERGY IN CM-1
            !
            ! THE PARAMETERS OF THE POTENTIAL ARE READ IN FROM UNIT 1
            ! (FILENAME GA.IN)
            !
            ! THE 6-D POTENTIAL IS BASED ON 1070 AB INITIO POINTS FROM
            ! KOFRANEK ET AL (CHEM PHYS 121(1988) 137) UP TO 42000 CM-1
            ! AND A PERTURBATION TREATMENT BY RIJKS AND WORMER
            ! (JCP 90(1989) 6507) AS WELL AS MONOMER We,WeXe FROM
            ! HUBER/HERZBERG. 29 FREE, 15 MONOMER OR AB INITIO
            ! BASED AND 18 CONSTRAINED PARAMETERS ARE USED.
            ! SQSBDE IS AN EMPIRICALLY REFINED VERSION WHICH HAS THE
            ! CORRECT EXPERIMENTAL HYDROGEN BOND WAVENUMBER OF
            ! 1062CM-1 (MILLER, ACC. CHEM. RES. 23 (1990) P10) AND
            ! B ROTATIONAL CONSTANT. FOR THIS PURPOSE, THE AB INITIO
            ! RAB COORDINATES WERE SCALED BY 1/1.035.
            ! THE RMS OF THE WEIGHTED FIT (W=1 FOR R1,2=/=1.7374 OR
            ! R1,2=1.7374 AND E<2220CM-1, W=0.1 FOR R1,2=1.7374 AND
            ! 2220<E<3210CM-1, W=0.001 FOR R1,2=1.7374 AND E>3210CM-1)
            ! IS 32.8CM-1.
            !
            ! E(1.7445,1.7404,5.1440,9.0007,64.1361,180)=-1.307CM-1 IS
            ! THE ABSOLUTE MINIMUM OF THE POTENTIAL
            !
            ! M. A. SUHM  SEPTEMBER 1989
            !
            ! preparatory calculations for gssqs-routine
            !
            implicit none
            double precision :: cg3j, e, esum, f, sum, x
            integer          :: i, i1, i2, ia, ib, ic, id, ie, if, j, j1, j2, j3, l, m, m1, m2, m3, nmax
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            integer, parameter          :: ma    = 135
            !
            !
            !
            !
            !
            !
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !
            !   References:      M. Quack, M. Suhm Mol. Phys. 69 (1990) P791
            !                    Kofranek et al Chem. Phys. 121 (1988) 137
            !                    Rijks and Wormer J. Chem. Phys. 90(1989) 6507
            infocm_ref(1) = 'M. Quack, M. Suhm, Mol. Phys. 69, P791(1990)'
            infocm_ref(2) = 'Kofranek et al., Chem. Phys. 121, 137(1988)'
            infocm_ref(3) = 'Rijks and Wormer, J. Chem. Phys. 90, 6507(1989)'
            !
            infocm_indexes(1) = 1
            infocm_indexes(2) = 9
            infocm_indexes(3) = 1
            infocm_indexes(4) = 9
            !
            infocm_irctnt = 3
            !
            !      CALL POTINFO
            !
            call ancvrt
            !
            ! prefactors for polynomials
            !
            potcm_fac(0) = 1.
            do i=1, 2*4
                potcm_fac(i) = potcm_fac(i-1)*i
            101 enddo
            do l=0, 4
                do m=0, l
                    ang_pf(l,m) = (-1)**m*sqrt((2*l+1)*potcm_fac(l-m)/(4*pi*potcm_fac(l+m)))
                103 enddo
            102 enddo
            !
            ! 3j-symbols for potential energy surface
            ! (abc/d-d0)=g(a,b,c,d)
            !
            !
            potcm_fct(1)  = 1.0d0
            potcm_efct(1) = 0.0d0
            do i=1, 400
                x             = i
                j             = i+1
                potcm_efct(j) = potcm_efct(i)
                potcm_fct(j)  = potcm_fct(i)*x
                do while( potcm_fct(j) .ge. 10.0d0 )
                    potcm_fct(j)  = 0.1d0*potcm_fct(j)
                    potcm_efct(j) = potcm_efct(j)+1.0d0
                enddo
            104 enddo
            !
            !
            m2 = 0
            m3 = 0
            do j1=0, 4*2, 2
                do j2=0, 4*2, 2
                    do j3=iabs(j1-j2), iabs(j1+j2), 4
                        do m1=0, min(j1, j2), 2
                            m2 = -m1
                            i1 = max0(0, (j2-j3-m1)/2, (j1-j3+m2)/2)
                            i2 = min0((j1+j2-j3)/2, (j1-m1)/2, (j2+m2)/2)
                            ia = i1
                            ib = (j1+j2-j3)/2-i1
                            ic = (j1-m1)/2-i1
                            id = (j2+m2)/2-i1
                            ie = (j3-j2+m1)/2+i1
                            if = (j3-j1-m2)/2+i1
                            sum  = (-1)**i1/(potcm_fct(ia+1)*potcm_fct(ib+1)*potcm_fct(ic+1)*potcm_fct(id+1)* &
                                & potcm_fct(ie+1)*potcm_fct(if+1))
                            esum = -potcm_efct(ia+1)-potcm_efct(ib+1)-potcm_efct(ic+1)-potcm_efct(id+1)- &
                                & potcm_efct(ie+1)-potcm_efct(if+1)
                            nmax = i2-i1
                            f    = 1.0d0
                            if( nmax .gt. 0 )then
                                do i=nmax, 1, -1
                                    x = (ib-i+1)*(ic-i+1)*(id-i+1)
                                    x = x/((ia+i)*(ie+i)*(if+i))
                                    f = 1.0d0-x*f
                                109 enddo
                            endif
                            sum = sum*f
                            f=potcm_fct((j1+j2-j3)/2+1)*potcm_fct((j1-j2+j3)/2+1)* potcm_fct((j2-j1+j3)/2+1)* &
                                & potcm_fct((j1+m1)/2+1)* potcm_fct((j1-m1)/2+1)* potcm_fct((j2+m2)/2+1)* &
                                & potcm_fct((j2-m2)/2+1)* potcm_fct((j3+m3)/2+1)* potcm_fct((j3-m3)/2+1)/ &
                                & potcm_fct((j1+j2+j3)/2+2)
                            e=potcm_efct((j1+j2-j3)/2+1)+potcm_efct((j1-j2+j3)/2+1)+ potcm_efct((j2-j1+j3)/2+1)+ &
                                & potcm_efct((j1+m1)/2+1)+ potcm_efct((j1-m1)/ 2+ 1)+ potcm_efct((j2+m2)/2+1)+ &
                                & potcm_efct((j2-m2)/2+1)+ potcm_efct((j3+m3)/2+1)+ potcm_efct((j3-m3)/2+1)- &
                                & potcm_efct((j1+ j2+j3)/ 2+2)
                            cg3j                       = sum*dsqrt(f)*10.0d0**(esum+0.5d0*e)* (-1)**((j2-j1-m3)/2)
                            ang_g(j1/2,j2/2,j3/2,m1/2) = cg3j
                        105 enddo
                    106 enddo
                107 enddo
            108 enddo
        end subroutine prepare_local_context
        !> @brief Evaluate gsrsqs for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine gsrsqs
            !
            implicit none
            double precision :: dab1, dab2, dab3
            integer          :: j1, j2, j3
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            integer, parameter :: ndmax  = 1100
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            integer, parameter          :: ma    = 135
            !
            !
            !
            !
            !
            !
            !      dab1=r-3.d0
            !      dab2=r-4.d0
            !      dab3=r-4.9d0
            !
            dab1 = pt1cm_r(1)-3.0d0
            dab2 = pt1cm_r(1)-4.0d0
            dab3 = pt1cm_r(1)-4.9d0
            do j1=0, 4, 1
                do j2=j1, 4, 1
                    do j3=(j2-j1), (j1+j2), 2
                        potcm_ar(j1,j2,j3)=potcm_ac(1,j1,j2,j3)*exp(-0.20d0*dab1**2)+ potcm_ac(2,j1,j2,j3)*exp(- &
                            & 0.35d0*dab2**2)+ potcm_ac(3,j1,j2,j3)* exp(-0.80d0*dab3**2)
                        potcm_ar(j2,j1,j3) = (-1)**(j1+j2)*potcm_ar(j1,j2,j3)
                    303 enddo
                302 enddo
            301 enddo
            !
            !
        end subroutine gsrsqs
        !> @brief Evaluate pot for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine pot
            !
            implicit none
            double precision :: amf, amh, c, c2, c3, c4, c5, cs, e, fac2, p1, p2, phicon, pif, r1, r12, r14, r2, r23, r34, rcm, s, s2, s3, s4, s5, tau, th1, th2, vp, x, y, z
            integer          :: icu, j, j1, j2, j3, m
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            integer, parameter          :: ma    = 135
            !
            !
            !
            !
            !
            !
            dimension x(4), y(4), z(4)
            dimension p1(0:4,0:4), p2(0:4,0:4), cs(4)
            dimension vp(0:4,0:4,0:2*4)
            !
            data amh/1837.152756d0/amf/34631.9678d0/
            !
            !     PUT COORDINATES IN PROPER ARRAYS
            !
            call cartou
            call carttor
            !
            !                                *
            !   Assign entries from R-array  *
            !   to variables used internally *
            !   by the POT subroutine        *
            !                                *
            !
            rcm    = pt1cm_r(1)
            r1     = pt1cm_r(2)
            r2     = pt1cm_r(3)
            th1    = pt1cm_r(4)
            th2    = pt1cm_r(5)
            tau    = pt1cm_r(6)
            phicon = 0.01745329252d0
            pif    = 1.0d0/(2.0d0*sqrt(pi))
            fac2   = amh/(amh+amf)
            x(3)   = 0.0d0-r1*cos(th1)*fac2
            y(3)   = 0.0d0-r1*sin(th1)*fac2
            z(3)   = 0.0d0
            x(1)   = r1*cos(th1)*(1.0d0-fac2)
            y(1)   = r1*sin(th1)*(1.0d0-fac2)
            z(1)   = 0.0d0
            x(4)   = rcm-r2*cos(th2)*fac2
            y(4)   = 0.0d0-r2*sin(th2)*cos(tau)*fac2
            z(4)   = 0.0d0-r2*sin(th2)*sin(tau)*fac2
            x(2)   = rcm+r2*cos(th2)*(1.0d0-fac2)
            y(2)   = 0.0d0+r2*sin(th2)*cos(tau)*(1.0d0-fac2)
            z(2)   = 0.0d0+r2*sin(th2)*sin(tau)*(1.0d0-fac2)
            r12    = sqrt((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
            r23    = sqrt((x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2)
            r14    = sqrt((x(4)-x(1))**2+(y(4)-y(1))**2+(z(4)-z(1))**2)
            r34    = sqrt((x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2)
            do j=1, 4
                cs(j) = cos(j*tau)
            101 enddo
            s       = sin(th1)
            c       = cos(th1)
            s2      = s*s
            s3      = s2*s
            s4      = s3*s
            s5      = s4*s
            c2      = c*c
            c3      = c2*c
            c4      = c3*c
            c5      = c4*c
            p1(0,0) = ang_pf(0,0)
            p1(1,0) = ang_pf(1,0)*c
            p1(2,0) = ang_pf(2,0)*0.5d0*(3.0d0*c2-1.0d0)
            p1(3,0) = ang_pf(3,0)*0.5d0*(5.0d0*c3-3.0d0*c)
            p1(4,0) = ang_pf(4,0)*0.125d0*(35.0d0*c4-30.0d0*c2+3.0d0)
            p1(1,1) = ang_pf(1,1)*s
            p1(2,1) = ang_pf(2,1)*3*s*c
            p1(3,1) = ang_pf(3,1)*1.5d0*s*(5.0d0*c2-1.0d0)
            p1(4,1) = ang_pf(4,1)*2.5d0*s*(7.0d0*c3-3.0d0*c)
            p1(2,2) = ang_pf(2,2)*3.0d0*s2
            p1(3,2) = ang_pf(3,2)*15.0d0*s2*c
            p1(4,2) = ang_pf(4,2)*7.5d0*s2*(7.0d0*c2-1.0d0)
            p1(3,3) = ang_pf(3,3)*15.0d0*s3
            p1(4,3) = ang_pf(4,3)*105.0d0*s3*c
            p1(4,4) = ang_pf(4,4)*105.0d0*s4
            s       = sin(th2)
            c       = cos(th2)
            s2      = s*s
            s3      = s2*s
            s4      = s3*s
            s5      = s4*s
            c2      = c*c
            c3      = c2*c
            c4      = c3*c
            c5      = c4*c
            p2(0,0) = ang_pf(0,0)
            p2(1,0) = ang_pf(1,0)*c
            p2(2,0) = ang_pf(2,0)*0.5d0*(3.0d0*c2-1.0d0)
            p2(3,0) = ang_pf(3,0)*0.5d0*(5.0d0*c3-3.0d0*c)
            p2(4,0) = ang_pf(4,0)*0.125d0*(35.0d0*c4-30.0d0*c2+3.0d0)
            p2(1,1) = ang_pf(1,1)*s
            p2(2,1) = ang_pf(2,1)*3.0d0*s*c
            p2(3,1) = ang_pf(3,1)*1.5d0*s*(5.0d0*c2-1.0d0)
            p2(4,1) = ang_pf(4,1)*2.5d0*s*(7.0d0*c3-3.0d0*c)
            p2(2,2) = ang_pf(2,2)*3.0d0*s2
            p2(3,2) = ang_pf(3,2)*15.0d0*s2*c
            p2(4,2) = ang_pf(4,2)*7.5d0*s2*(7.0d0*c2-1.0d0)
            p2(3,3) = ang_pf(3,3)*15.0d0*s3
            p2(4,3) = ang_pf(4,3)*105.0d0*s3*c
            p2(4,4) = ang_pf(4,4)*105.0d0*s4
            icu     = 31
            do j1=0, 4, 1
                do j2=j1, 4, 1
                    do j3=(j2-j1), (j1+j2), 2
                        potcm_ac(1,j1,j2,j3) = potcm_a(icu)
                        potcm_ac(2,j1,j2,j3) = potcm_a(icu+35)
                        potcm_ac(3,j1,j2,j3) = potcm_a(icu+70)
                        icu                  = icu+1
                    104 enddo
                103 enddo
            102 enddo
            !
            call gsrsqs
            !
            do j1=0, 4, 1
                do j2=0, 4, 1
                    do j3=abs(j2-j1), (j2+j1), 2
                        vp(j1,j2,j3) = ang_g(j1,j2,j3,0)*p1(j1,0)*p2(j2,0)
                        do m=1, min(j1, j2), 1
                            vp(j1,j2,j3) = vp(j1,j2,j3)+(-1)**m*2*ang_g(j1,j2,j3,m)* cs(m)*p1(j1,m)*p2(j2,m)
                        108 enddo
                        vp(j1,j2,j3) = vp(j1,j2,j3)*(-1)**(j1-j2)*(2*j3+1)*pif
                    107 enddo
                106 enddo
            105 enddo
            !
            e = 0.0d0
            !
            do j1=0, 4, 1
                do j2=0, 4, 1
                    do j3=abs(j2-j1), (j2+j1), 2
                        e = e+vp(j1,j2,j3)*potcm_ar(j1,j2,j3)
                    111 enddo
                110 enddo
            109 enddo
            !
            !     e=e+a(1)
            !
            if( rcm .lt. 7.2d0 )then
                e = e+vp(0,0,0)*exp(-((7.2d0-rcm)/rcm)**2)*potcm_a(2)/rcm**6
            else
                e = e+vp(0,0,0)*potcm_a(2)/rcm**6
            endif
            if( rcm .lt. 8.4d0 )then
                e = e+vp(0,0,0)*exp(-((8.4d0-rcm)/rcm)**2)*potcm_a(3)/rcm**8
            else
                e = e+vp(0,0,0)*potcm_a(3)/rcm**8
            endif
            if( rcm .lt. 9.6d0 )then
                e = e+vp(0,0,0)*exp(-((9.6d0-rcm)/rcm)**2)*potcm_a(4)/rcm**10
            else
                e = e+vp(0,0,0)*potcm_a(4)/rcm**10
            endif
            if( rcm .lt. 9.5d0 )then
                e = e+vp(1,1,2)*potcm_a(5)/rcm**3*exp(-((9.5d0-rcm)/rcm)**2)
            else
                e = e+vp(1,1,2)*potcm_a(5)/rcm**3
            endif
            if( rcm .lt. 10.5d0 )then
                e = e+(vp(2,1,3)-vp(1,2,3))*potcm_a(6)/rcm**4* exp(-((10.5d0-rcm)/rcm)**2)
            else
                e = e+(vp(2,1,3)-vp(1,2,3))*potcm_a(6)/rcm**4
            endif
            if( rcm .lt. 11.5d0 )then
                e = e+vp(2,2,4)*potcm_a(7)/rcm**5*exp(-((11.5d0-rcm)/rcm)**2)
                e = e+(vp(3,1,4)+vp(1,3,4))*potcm_a(8)/rcm**5*exp(-((11.5-rcm)/rcm)**2)
            else
                e = e+vp(2,2,4)*potcm_a(7)/rcm**5
                e = e+(vp(3,1,4)+vp(1,3,4))*potcm_a(8)/rcm**5
            endif
            e = e+potcm_a(9)*exp(-0.5d0*r12)
            e = e+potcm_a(10)*exp(-1.5d0*r12)
            e = e+potcm_a(13)*exp(-1.5d0*r34)
            e = e+potcm_a(17)*(exp(-10.0d0*r14)+exp(-10.0d0*r23))
            if( rcm .lt. 4.5d0 )then
                e = e+potcm_a(18)/rcm**6*((3.0d0*(cos(th1))**2+1.0d0)+ (3.0d0*(cos(th2))**2+1.0d0))* exp(-((4.5d0- &
                    & rcm)/rcm)**2)
            else
                e = e+potcm_a(18)/rcm**6*((3.0d0*(cos(th1))**2+1.0d0)+ (3.0d0*(cos(th2))**2+1.0d0))
            endif
            !
            ! monomer terms
            !
            e=e+47635d0*(1.0d0-potcm_a(19)*exp(-0.1d0*(r23-3.0d0)**2))* (1.0d0-exp(-1.1953d0* (1.0d0-potcm_a(26)* &
                & (1.0d0-tanh(potcm_a(27)*(r14- 3.0d0))) )* (r1-1.7374d0*(1.0d0 -potcm_a(20)* (1.0d0-tanh(potcm_a(21)* &
                & (r14-3.0d0))) ))))**2
            e=e+47635d0*(1.0d0-potcm_a(19)*exp(-0.1d0*(r14-3.0d0)**2))* (1.0d0-exp(-1.1953d0* (1.0d0-potcm_a(26)* &
                & (1.0d0-tanh(potcm_a(27)*(r23- 3.0d0))) )* (r2-1.7374d0*(1.0d0 -potcm_a(20)* (1.0d0-tanh(potcm_a(21)* &
                & (r23-3.0d0))) ))))**2
            !
            !     next line corrects for asymptotic diatomic potential which is
            !     added in separately.
            !
            e = e-(ptntl(r1-1.7374d0)+ptntl(r2-1.7374d0))*219474.63067d0
            e = e+potcm_a(28)*(exp(-4.0d0*(r14-r1)**2)+exp(-4.0d0*(r23-r2)**2))
            e = e+potcm_a(22)*exp(-4.0d0*r14)*exp(-4.0d0*r23)
            e = e+potcm_a(23)
            e = e+potcm_a(24)*(exp(-r14)+exp(-r23))* (exp(-15.0d0*r1)+exp(-15.0d0*r2))
            !
            !   Units conversion:  Conversion to atomic units
            !
            e = e/219474.63067d0
            !
            pt1cm_engygs = e
            !
            !
            call eunitzero
            !
            if( usricm_nder .ne. 0 )then
                call rtocart
                call dedcou
            endif
            !
        end subroutine pot
        !> @brief Evaluate PTNTL for the F2H2_GEN_SQSBDE_1995 potential.
        function ptntl( x )
            !
            !     HF DIATOMIC POTENTIAL CURVE FOR SQSBDE SURFACE, CALCULATED USING
            !     PARAMETERS GIVEN IN JCP, VOL. 95(5), P. 31.
            !     NOTICE THAT THE PARAMETER DM IS GIVEN IN ATOMIC UNITS HERE AND
            !     IN CM-1 IN THE ICP REFERENCE.
            !
            implicit none
            double precision :: ptntl, x
            !
            !
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            !
            double precision, parameter :: pi    = 3.141592653589793d0
            integer, parameter          :: natom = 25
            integer, parameter          :: ma    = 135
            double precision, parameter :: am    = 1.1953d0
            double precision, parameter :: dm    = 0.21704d0
            double precision, parameter :: one   = 1.0d0
            !
            !
            !
            !
            !
            !
            !      SAVE
            !
            ptntl = dm*(one-exp(-am*x))**2
            !
        end function ptntl
        !> @brief Evaluate POTINFO for the F2H2_GEN_SQSBDE_1995 potential.
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
                write(usricm_nflag(18), '(2x,''Cartesian coordinates are supplied by'',/, 2x,''the user in the array CART.'',//)')
                write(usricm_nflag(18), '(2x,''Provide cartesian coordinates in the'',/, 2x,''following order using the array CART'',' // &
                    '//, 2x,'' CART(1,1)...CART(1,3)   => ATOM 1'',/, 2x,'' CART(2,1)...CART(2,3)   => AT' // &
                    'OM 2'',/, 2x,'' CART(3,1)...CART(3,3)   => ATOM 3'',/, 2x,'' CART(N,1)...CART(N,3)  ' // &
                    ' => ATOM N'',/, 2x,''CART(25,1)...CART(25,3)  => ATOM 25'',/)')
                write(usricm_nflag(18), '(2x,''If the user wishes to relabel the atoms,'',/, 2x,''set the variable IREORDER equal to ' // &
                    '1'',/, 2x,''in the PARAMETER statement.  The user'',/, 2x,''must also specify the ne' // &
                    'w labeling'',/, 2x,''scheme.  This is done using the array'',/, 2x,''NULBL in the fo' // &
                    'llowing manner:'',//, 2x,''NULBL(i) = j'',/, 2x,''where:  i => old index'',/, 2x,'' ' // &
                    '       j => new index'',//)')
                write(usricm_nflag(18), '(2x,''Cartesian coordinates can be provided to'',/, 2x,''the potential routine in a variety ' // &
                    'of units.'',/, 2x,''The input units will be converted to Bohr'',/, 2x,''based on the' // &
                    ' following values of the NFLAG'',/, 2x,''variable:'',//, 2x,''NFLAG(1)  =  1  =>  CA' // &
                    'RTESIANS IN BOHR (no'',/, 2x,''                    conversion required)'',/, 2x,''NF' // &
                    'LAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS'',//)')
                write(usricm_nflag(18), '(2x,''The value of the energy and derivatives'',/, 2x,''(if computed) can be reported in a v' // &
                    'ariety'',/, 2x,''units.  A units conversion will take place'',/, 2x,''as specified b' // &
                    'y the following values of the'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(2) = 1 =>  ' // &
                    'ENERGIES REPORTED IN HARTEEE'',/, 2x,''NFLAG(2) = 2 =>  ENERGIES REPORTED IN mHARTRE' // &
                    'E'',/, 2x,''NFLAG(2) = 3 =>  ENERGIES REPORTED IN eV'',/, 2x,''NFLAG(2) = 4 =>  ENER' // &
                    'GIES REPORTED IN kcal/mol'',/, 2x,''NFLAG(2) = 5 =>  ENERGIES REPORTED IN cm**-1'',/' // &
                    '/)')
                write(usricm_nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
                    'he'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(1)=1 & NFLAG(2)=1 => DERIVATIVES REPOR' // &
                    'TED IN'',/, 2x,''                           HARTEEE/BOHR'',/, 2x,''NFLAG(1)=1 & NFLA' // &
                    'G(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/BOHR' // &
                    ''',/, 2x,''NFLAG(1)=1 & NFLAG(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''            ' // &
                    '               eV/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=4 => DERIVATIVES REPORTED IN' // &
                    ''',/, 2x,''                           kcal/mol/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=' // &
                    '5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm**-1/BOHR'',//)')
                write(usricm_nflag(18), '(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
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
        !> @brief Evaluate ANCVRT for the F2H2_GEN_SQSBDE_1995 potential.
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
        !> @brief Evaluate CARTOU for the F2H2_GEN_SQSBDE_1995 potential.
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
        !> @brief Evaluate CARTTOR for the F2H2_GEN_SQSBDE_1995 potential.
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
                write(usricm_nflag(18), '(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') infocm_icartr
                stop
            endif
        end subroutine carttor
        !> @brief Evaluate EUNITZERO for the F2H2_GEN_SQSBDE_1995 potential.
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
        !> @brief Evaluate RTOCART for the F2H2_GEN_SQSBDE_1995 potential.
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
                        20 continue
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
                write(usricm_nflag(18), '(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') infocm_icartr
                stop
            endif
        end subroutine rtocart
        !> @brief Evaluate DEDCOU for the F2H2_GEN_SQSBDE_1995 potential.
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
end module chempotpy_f2h2_gen_sqsbde_1995_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_f2h2_gen_sqsbde_1995_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 4
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates)
    double precision, intent(out) :: g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
