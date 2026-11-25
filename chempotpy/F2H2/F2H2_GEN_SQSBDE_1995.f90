!   System:                        F2H2
!   Common name:                   SQSBDE
!   Number of derivatives:         0
!   Number of bodies:              4
!   Number of electronic surfaces: 1
!   Interface:                     potlib2001
!   References:      M. Quack, M. Suhm Mol. Phys. 69 (1990) P791
!                    Kofranek et al Chem. Phys. 121 (1988) 137
!                    Rijks and Wormer J. Chem. Phys. 90(1989) 6507
!   Protocol:
!      prepare_local_context - initializes the potential's variables
!               must be called once before any calls to POT
!      POT    - driver for the evaluation of the energy
!   Units:
!      energies    - cm-1
!      coordinates - distances: bohr; angles: degrees
!   Surfaces:
!      ground electronic state
! M. A. SUHM  SEPTEMBER 1990
! INPUT : SIX INTERNAL COORDINATES
! --R1,R2,R IN ATOMIC UNITS (BOHR)
! --TH1,TH2,TAU IN DEGREES
! DEFINITION AS IN QUACK+SUHM MOL.PHYS. 69 (1990) P791
! OUTPUT : HF DIMER (SQSBDE) POTENTIAL ENERGY IN CM-1
! THE POTENTIAL PARAMETERS ARE EMBEDDED BELOW
! (FILENAME GA.IN)
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
! E(1.7445,1.7404,5.1440,9.0007,64.1361,180)=-1.307CM-1 IS
! THE ABSOLUTE MINIMUM OF THE POTENTIAL
! M. A. SUHM  SEPTEMBER 1989
! preparatory calculations for gssqs-routine

!> @brief Implement F2H2_GEN_SQSBDE_1995, a potential energy surface for F2H2.
module chempotpy_f2h2_gen_sqsbde_1995_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the F2H2_GEN_SQSBDE_1995 potential energy surface for F2H2.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: electronvolt_to_inverse_centimeter
        use chempotpy_units,  only: bohr_to_angstrom
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 4
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter    :: ma      = 135
        integer, parameter    :: n3atom  = 75
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: ang_g(0:4,0:4,0:2*4,0:4)
        real(wp)              :: ang_pf(0:4,0:4)
        real(wp)              :: cartesian_coordinates(natom,3)
        integer               :: coordinate_representation
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: atom_count
        character(len=75)     :: infocm_ref(5)
        real(wp)              :: potcm_a(ma)
        real(wp)              :: potcm_ac(3,0:4,0:4,0:2*4)
        real(wp)              :: potcm_ar(0:4,0:4,0:2*4)
        real(wp)              :: potcm_da(ma)
        real(wp)              :: potcm_efct(401)
        real(wp)              :: potcm_fac(0:2*4)
        real(wp)              :: potcm_fct(401)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: ground_state_energy
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: input_cartesian(natom,3)
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        ang_g                             = 0.e0_wp; ang_pf = 0.e0_wp; cartesian_coordinates = 0.e0_wp; coordinate_representation = 0
        maximum_derivative_order          = 0; maximum_excited_state = 0; atom_count = 0
        infocm_ref                        = ' '
        potcm_a                           = 0.e0_wp
        potcm_ac                          = 0.e0_wp
        potcm_ar                          = 0.e0_wp
        potcm_da                          = 0.e0_wp
        potcm_efct                        = 0.e0_wp
        potcm_fac                         = 0.e0_wp
        potcm_fct                         = 0.e0_wp
        ground_state_internal_gradient    = 0.e0_wp
        ground_state_energy               = 0.e0_wp
        pt1cm_r                           = 0.e0_wp
        pt3cm_ezero                       = 0.e0_wp
        excited_state_internal_gradients  = 0.e0_wp
        excited_state_energies            = 0.e0_wp
        diabatic_internal_gradients       = 0.e0_wp
        diabatic_couplings                = 0.e0_wp
        input_cartesian                   = 0.e0_wp
        surface_requests                  = 0
        requested_derivative_order        = 0
        output_excited_state_gradients    = 0.e0_wp
        output_ground_state_gradient      = 0.e0_wp
        output_diabatic_gradients         = 0.e0_wp
        output_excited_state_energies     = 0.e0_wp
        output_ground_state_energy        = 0.e0_wp
        output_diabatic_couplings         = 0.e0_wp
        excited_state_cartesian_gradients = 0.e0_wp
        ground_state_cartesian_gradient   = 0.e0_wp
        diabatic_cartesian_gradients      = 0.e0_wp
        last_requested_excited_state      = 0
        last_requested_coupling           = 0
        first_requested_excited_state     = 0
        first_requested_coupling          = 0
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 0
        coordinate_representation  = 4
        maximum_excited_state      = 0
        maximum_derivative_order   = 0
        atom_count                 = 4
        potcm_a(1)                 = 1558.0e0_wp
        potcm_da(1)                = 0.0e0_wp
        potcm_a(2)                 = -202900000.0e0_wp
        potcm_da(2)                = 0.0e0_wp
        potcm_a(3)                 = -3703000000.0e0_wp
        potcm_da(3)                = 0.0e0_wp
        potcm_a(4)                 = -11350000000.0e0_wp
        potcm_da(4)                = 0.0e0_wp
        potcm_a(5)                 = -1736693.09500000e0_wp
        potcm_da(5)                = 0.0e0_wp
        potcm_a(6)                 = -4416487.722000003e0_wp
        potcm_da(6)                = 0.0e0_wp
        potcm_a(7)                 = 16012389.5100000e0_wp
        potcm_da(7)                = 0.0e0_wp
        potcm_a(8)                 = -6389444.85500000e0_wp
        potcm_da(8)                = 0.0e0_wp
        potcm_a(9)                 = -866.867375551329e0_wp
        potcm_da(9)                = 117.589789927560e0_wp
        potcm_a(10)                = 129789.363470169e0_wp
        potcm_da(10)               = 1937.56316874928e0_wp
        potcm_a(11)                = 0.0e0_wp
        potcm_da(11)               = 0.0e0_wp
        potcm_a(12)                = 0.0e0_wp
        potcm_da(12)               = 0.0e0_wp
        potcm_a(13)                = 2622184.75665877e0_wp
        potcm_da(13)               = 29726.6301460444e0_wp
        potcm_a(14)                = 0.0e0_wp
        potcm_da(14)               = 0.0e0_wp
        potcm_a(15)                = 0.0e0_wp
        potcm_da(15)               = 0.0e0_wp
        potcm_a(16)                = 0.0e0_wp
        potcm_da(16)               = 0.0e0_wp
        potcm_a(17)                = 184254701843.311e0_wp
        potcm_da(17)               = 17553092934.4463e0_wp
        potcm_a(18)                = -320000.0e0_wp
        potcm_da(18)               = 0.0e0_wp
        potcm_a(19)                = 3.326438943860246e-2_wp
        potcm_da(19)               = 1.353645657747023e-3_wp
        potcm_a(20)                = -4.576887206887037e-3_wp
        potcm_da(20)               = 1.865017706495822e-4_wp
        potcm_a(21)                = 0.409730666039653e0_wp
        potcm_da(21)               = 2.777020547964313e-2_wp
        potcm_a(22)                = 32623551943283.8e0_wp
        potcm_da(22)               = 3630087056915.05e0_wp
        potcm_a(23)                = 0.0e0_wp
        potcm_da(23)               = 0.0e0_wp
        potcm_a(24)                = 71754497894201.8e0_wp
        potcm_da(24)               = 2525041881368.02e0_wp
        potcm_a(25)                = 0.0e0_wp
        potcm_da(25)               = 0.0e0_wp
        potcm_a(26)                = 5.032387872479119e-2_wp
        potcm_da(26)               = 1.860431838521344e-3_wp
        potcm_a(27)                = 1.46832057915281e0_wp
        potcm_da(27)               = 9.713248031677595e-2_wp
        potcm_a(28)                = 3384.63155928653e0_wp
        potcm_da(28)               = 85.5567822945890e0_wp
        potcm_a(29)                = 0.0e0_wp
        potcm_da(29)               = 0.0e0_wp
        potcm_a(30)                = 0.0e0_wp
        potcm_da(30)               = 0.0e0_wp
        potcm_a(31)                = -160267.101272529e0_wp
        potcm_da(31)               = 4055.33877671770e0_wp
        potcm_a(32)                = 0.0e0_wp
        potcm_da(32)               = 0.0e0_wp
        potcm_a(33)                = 0.0e0_wp
        potcm_da(33)               = 0.0e0_wp
        potcm_a(34)                = 0.0e0_wp
        potcm_da(34)               = 0.0e0_wp
        potcm_a(35)                = 0.0e0_wp
        potcm_da(35)               = 0.0e0_wp
        potcm_a(36)                = -16268.8946598104e0_wp
        potcm_da(36)               = 541.139677344219e0_wp
        potcm_a(37)                = 0.0e0_wp
        potcm_da(37)               = 0.0e0_wp
        potcm_a(38)                = 0.0e0_wp
        potcm_da(38)               = 0.0e0_wp
        potcm_a(39)                = 2802.57610197247e0_wp
        potcm_da(39)               = 56.8038283168248e0_wp
        potcm_a(40)                = 0.0e0_wp
        potcm_da(40)               = 0.0e0_wp
        potcm_a(41)                = 0.0e0_wp
        potcm_da(41)               = 0.0e0_wp
        potcm_a(42)                = 0.0e0_wp
        potcm_da(42)               = 0.0e0_wp
        potcm_a(43)                = 0.0e0_wp
        potcm_da(43)               = 0.0e0_wp
        potcm_a(44)                = 0.0e0_wp
        potcm_da(44)               = 0.0e0_wp
        potcm_a(45)                = 0.0e0_wp
        potcm_da(45)               = 0.0e0_wp
        potcm_a(46)                = 2860.06394339599e0_wp
        potcm_da(46)               = 44.3633384808498e0_wp
        potcm_a(47)                = 0.0e0_wp
        potcm_da(47)               = 0.0e0_wp
        potcm_a(48)                = 0.0e0_wp
        potcm_da(48)               = 0.0e0_wp
        potcm_a(49)                = 0.0e0_wp
        potcm_da(49)               = 0.0e0_wp
        potcm_a(50)                = 0.0e0_wp
        potcm_da(50)               = 0.0e0_wp
        potcm_a(51)                = 0.0e0_wp
        potcm_da(51)               = 0.0e0_wp
        potcm_a(52)                = 509.102022118840e0_wp
        potcm_da(52)               = 18.7375880046217e0_wp
        potcm_a(53)                = 0.0e0_wp
        potcm_da(53)               = 0.0e0_wp
        potcm_a(54)                = 0.0e0_wp
        potcm_da(54)               = 0.0e0_wp
        potcm_a(55)                = 0.0e0_wp
        potcm_da(55)               = 0.0e0_wp
        potcm_a(56)                = 0.0e0_wp
        potcm_da(56)               = 0.0e0_wp
        potcm_a(57)                = 0.0e0_wp
        potcm_da(57)               = 0.0e0_wp
        potcm_a(58)                = 0.0e0_wp
        potcm_da(58)               = 0.0e0_wp
        potcm_a(59)                = 0.0e0_wp
        potcm_da(59)               = 0.0e0_wp
        potcm_a(60)                = 0.0e0_wp
        potcm_da(60)               = 0.0e0_wp
        potcm_a(61)                = 0.0e0_wp
        potcm_da(61)               = 0.0e0_wp
        potcm_a(62)                = 0.0e0_wp
        potcm_da(62)               = 0.0e0_wp
        potcm_a(63)                = 0.0e0_wp
        potcm_da(63)               = 0.0e0_wp
        potcm_a(64)                = 0.0e0_wp
        potcm_da(64)               = 0.0e0_wp
        potcm_a(65)                = 0.0e0_wp
        potcm_da(65)               = 0.0e0_wp
        potcm_a(66)                = 0.0e0_wp
        potcm_da(66)               = 0.0e0_wp
        potcm_a(67)                = 0.0e0_wp
        potcm_da(67)               = 0.0e0_wp
        potcm_a(68)                = -3467.23417017878e0_wp
        potcm_da(68)               = 52.9105695619519e0_wp
        potcm_a(69)                = 0.0e0_wp
        potcm_da(69)               = 0.0e0_wp
        potcm_a(70)                = 0.0e0_wp
        potcm_da(70)               = 0.0e0_wp
        potcm_a(71)                = 0.0e0_wp
        potcm_da(71)               = 0.0e0_wp
        potcm_a(72)                = 0.0e0_wp
        potcm_da(72)               = 0.0e0_wp
        potcm_a(73)                = 2304.6328580383e0_wp
        potcm_da(73)               = 44.0232721532835e0_wp
        potcm_a(74)                = 0.0e0_wp
        potcm_da(74)               = 0.0e0_wp
        potcm_a(75)                = 0.0e0_wp
        potcm_da(75)               = 0.0e0_wp
        potcm_a(76)                = 0.0e0_wp
        potcm_da(76)               = 0.0e0_wp
        potcm_a(77)                = 0.0e0_wp
        potcm_da(77)               = 0.0e0_wp
        potcm_a(78)                = 0.0e0_wp
        potcm_da(78)               = 0.0e0_wp
        potcm_a(79)                = 0.0e0_wp
        potcm_da(79)               = 0.0e0_wp
        potcm_a(80)                = 0.0e0_wp
        potcm_da(80)               = 0.0e0_wp
        potcm_a(81)                = 0.0e0_wp
        potcm_da(81)               = 0.0e0_wp
        potcm_a(82)                = 0.0e0_wp
        potcm_da(82)               = 0.0e0_wp
        potcm_a(83)                = 0.0e0_wp
        potcm_da(83)               = 0.0e0_wp
        potcm_a(84)                = -644.815083021031e0_wp
        potcm_da(84)               = 18.1151965381963e0_wp
        potcm_a(85)                = 0.0e0_wp
        potcm_da(85)               = 0.0e0_wp
        potcm_a(86)                = 0.0e0_wp
        potcm_da(86)               = 0.0e0_wp
        potcm_a(87)                = 0.0e0_wp
        potcm_da(87)               = 0.0e0_wp
        potcm_a(88)                = 0.0e0_wp
        potcm_da(88)               = 0.0e0_wp
        potcm_a(89)                = 0.0e0_wp
        potcm_da(89)               = 0.0e0_wp
        potcm_a(90)                = 0.0e0_wp
        potcm_da(90)               = 0.0e0_wp
        potcm_a(91)                = 0.0e0_wp
        potcm_da(91)               = 0.0e0_wp
        potcm_a(92)                = 0.0e0_wp
        potcm_da(92)               = 0.0e0_wp
        potcm_a(93)                = 0.0e0_wp
        potcm_da(93)               = 0.0e0_wp
        potcm_a(94)                = 0.0e0_wp
        potcm_da(94)               = 0.0e0_wp
        potcm_a(95)                = 0.0e0_wp
        potcm_da(95)               = 0.0e0_wp
        potcm_a(96)                = 0.0e0_wp
        potcm_da(96)               = 0.0e0_wp
        potcm_a(97)                = 0.0e0_wp
        potcm_da(97)               = 0.0e0_wp
        potcm_a(98)                = 0.0e0_wp
        potcm_da(98)               = 0.0e0_wp
        potcm_a(99)                = 0.0e0_wp
        potcm_da(99)               = 0.0e0_wp
        potcm_a(100)               = 0.0e0_wp
        potcm_da(100)              = 0.0e0_wp
        potcm_a(101)               = 22530.9324112547e0_wp
        potcm_da(101)              = 989.768118754483e0_wp
        potcm_a(102)               = 1630.40722103571e0_wp
        potcm_da(102)              = 111.781936229583e0_wp
        potcm_a(103)               = 0.0e0_wp
        potcm_da(103)              = 0.0e0_wp
        potcm_a(104)               = 1042.90771467500e0_wp
        potcm_da(104)              = 18.4325305864792e0_wp
        potcm_a(105)               = -366.368726725276e0_wp
        potcm_da(105)              = 12.8434812694663e0_wp
        potcm_a(106)               = 2428.17244795167e0_wp
        potcm_da(106)              = 246.109466713404e0_wp
        potcm_a(107)               = -1882.37787794149e0_wp
        potcm_da(107)              = 70.8301781102542e0_wp
        potcm_a(108)               = 0.0e0_wp
        potcm_da(108)              = 0.0e0_wp
        potcm_a(109)               = 0.0e0_wp
        potcm_da(109)              = 0.0e0_wp
        potcm_a(110)               = -326.002341194829e0_wp
        potcm_da(110)              = 19.3580274207174e0_wp
        potcm_a(111)               = 0.0e0_wp
        potcm_da(111)              = 0.0e0_wp
        potcm_a(112)               = 0.0e0_wp
        potcm_da(112)              = 0.0e0_wp
        potcm_a(113)               = 116.704136523402e0_wp
        potcm_da(113)              = 11.3452311839242e0_wp
        potcm_a(114)               = 0.0e0_wp
        potcm_da(114)              = 0.0e0_wp
        potcm_a(115)               = 572.632667044664e0_wp
        potcm_da(115)              = 33.6908429103859e0_wp
        potcm_a(116)               = 0.0e0_wp
        potcm_da(116)              = 0.0e0_wp
        potcm_a(117)               = 0.0e0_wp
        potcm_da(117)              = 0.0e0_wp
        potcm_a(118)               = 0.0e0_wp
        potcm_da(118)              = 0.0e0_wp
        potcm_a(119)               = 0.0e0_wp
        potcm_da(119)              = 0.0e0_wp
        potcm_a(120)               = 0.0e0_wp
        potcm_da(120)              = 0.0e0_wp
        potcm_a(121)               = 0.0e0_wp
        potcm_da(121)              = 0.0e0_wp
        potcm_a(122)               = 0.0e0_wp
        potcm_da(122)              = 0.0e0_wp
        potcm_a(123)               = 0.0e0_wp
        potcm_da(123)              = 0.0e0_wp
        potcm_a(124)               = 0.0e0_wp
        potcm_da(124)              = 0.0e0_wp
        potcm_a(125)               = 0.0e0_wp
        potcm_da(125)              = 0.0e0_wp
        potcm_a(126)               = 0.0e0_wp
        potcm_da(126)              = 0.0e0_wp
        potcm_a(127)               = 0.0e0_wp
        potcm_da(127)              = 0.0e0_wp
        potcm_a(128)               = 0.0e0_wp
        potcm_da(128)              = 0.0e0_wp
        potcm_a(129)               = 0.0e0_wp
        potcm_da(129)              = 0.0e0_wp
        potcm_a(130)               = 0.0e0_wp
        potcm_da(130)              = 0.0e0_wp
        potcm_a(131)               = 0.0e0_wp
        potcm_da(131)              = 0.0e0_wp
        potcm_a(132)               = 0.0e0_wp
        potcm_da(132)              = 0.0e0_wp
        potcm_a(133)               = 0.0e0_wp
        potcm_da(133)              = 0.0e0_wp
        potcm_a(134)               = 0.0e0_wp
        potcm_da(134)              = 0.0e0_wp
        potcm_a(135)               = 0.0e0_wp
        potcm_da(135)              = 0.0e0_wp
        v                          = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        input_cartesian            = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                input_cartesian(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = output_ground_state_energy/electronvolt_to_inverse_centimeter
            enddo
        else
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives are unavailable')
        endif

        contains

        !> @brief Evaluate prepare_local_context for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: cg3j, e, esum, f, sum, x
            integer  :: i, i1, i2, ia, ib, ic, id, ie, if, j, j1, j2, j3, l, m, m1, m2, m3, nmax
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !   References:      M. Quack, M. Suhm Mol. Phys. 69 (1990) P791
            !                    Kofranek et al Chem. Phys. 121 (1988) 137
            !                    Rijks and Wormer J. Chem. Phys. 90(1989) 6507
            infocm_ref(1) = 'M. Quack, M. Suhm, Mol. Phys. 69, P791(1990)'
            infocm_ref(2) = 'Kofranek et al., Chem. Phys. 121, 137(1988)'
            infocm_ref(3) = 'Rijks and Wormer, J. Chem. Phys. 90, 6507(1989)'
            ! prefactors for polynomials
            potcm_fac(0) = 1.
            do i=1, 2*4
                potcm_fac(i) = potcm_fac(i-1)*i
            enddo
            do l=0, 4
                do m=0, l
                    ang_pf(l,m) = (-1)**m*sqrt((2*l+1)*potcm_fac(l-m)/(4*PI*potcm_fac(l+m)))
                enddo
            enddo
            ! 3j-symbols for potential energy surface
            ! (abc/d-d0)=g(a,b,c,d)
            potcm_fct(1)  = 1.0e0_wp
            potcm_efct(1) = 0.0e0_wp
            do i=1, 400
                x             = i
                j             = i+1
                potcm_efct(j) = potcm_efct(i)
                potcm_fct(j)  = potcm_fct(i)*x
                do while( potcm_fct(j) .ge. 10.0e0_wp )
                    potcm_fct(j)  = 0.1e0_wp*potcm_fct(j)
                    potcm_efct(j) = potcm_efct(j)+1.0e0_wp
                enddo
            enddo
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
                            f    = 1.0e0_wp
                            if( nmax .gt. 0 )then
                                do i=nmax, 1, -1
                                    x = (ib-i+1)*(ic-i+1)*(id-i+1)
                                    x = x/((ia+i)*(ie+i)*(if+i))
                                    f = 1.0e0_wp-x*f
                                enddo
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
                            cg3j                       = sum*sqrt(f)*10.0e0_wp**(esum+0.5e0_wp*e)* (-1)**((j2-j1-m3)/2)
                            ang_g(j1/2,j2/2,j3/2,m1/2) = cg3j
                        enddo
                    enddo
                enddo
            enddo
        end subroutine prepare_local_context

        !> @brief Evaluate gsrsqs for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine gsrsqs
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dab1, dab2, dab3
            integer  :: j1, j2, j3
            !      dab1=r-3.d0
            !      dab2=r-4.d0
            !      dab3=r-4.9d0
            dab1 = pt1cm_r(1)-3.0e0_wp
            dab2 = pt1cm_r(1)-4.0e0_wp
            dab3 = pt1cm_r(1)-4.9e0_wp
            do j1=0, 4, 1
                do j2=j1, 4, 1
                    do j3=(j2-j1), (j1+j2), 2
                        potcm_ar(j1,j2,j3)=potcm_ac(1,j1,j2,j3)*exp(-0.20e0_wp*dab1**2)+ potcm_ac(2,j1,j2,j3)*exp(- &
                            & 0.35e0_wp*dab2**2)+ potcm_ac(3,j1,j2,j3)* exp(-0.80e0_wp*dab3**2)
                        potcm_ar(j2,j1,j3) = (-1)**(j1+j2)*potcm_ar(j1,j2,j3)
                    enddo
                enddo
            enddo
        end subroutine gsrsqs

        !> @brief Evaluate pot for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine pot
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: amf, amh, c, c2, c3, c4, c5, cs(4), e, fac2, p1(0:4,0:4), p2(0:4,0:4), phicon, pif, r1, &
                & r12, r14, r2, r23, r34, rcm, s, s2, s3, s4, s5, tau, th1, th2, vp(0:4, 0:4, 0:2*4), x(4), y(4), z(4)
            integer :: icu, j, j1, j2, j3, m
            !     PUT COORDINATES IN PROPER ARRAYS
            amh                                   = 1837.152756e0_wp
            amf                                   = 34631.9678e0_wp
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !                                *
            !   Assign entries from R-array  *
            !   to variables used internally *
            !   by the POT subroutine        *
            !                                *
            rcm    = pt1cm_r(1)
            r1     = pt1cm_r(2)
            r2     = pt1cm_r(3)
            th1    = pt1cm_r(4)
            th2    = pt1cm_r(5)
            tau    = pt1cm_r(6)
            phicon = 0.01745329252e0_wp
            pif    = 1.0e0_wp/(2.0e0_wp*sqrt(PI))
            fac2   = amh/(amh+amf)
            x(3)   = 0.0e0_wp-r1*cos(th1)*fac2
            y(3)   = 0.0e0_wp-r1*sin(th1)*fac2
            z(3)   = 0.0e0_wp
            x(1)   = r1*cos(th1)*(1.0e0_wp-fac2)
            y(1)   = r1*sin(th1)*(1.0e0_wp-fac2)
            z(1)   = 0.0e0_wp
            x(4)   = rcm-r2*cos(th2)*fac2
            y(4)   = 0.0e0_wp-r2*sin(th2)*cos(tau)*fac2
            z(4)   = 0.0e0_wp-r2*sin(th2)*sin(tau)*fac2
            x(2)   = rcm+r2*cos(th2)*(1.0e0_wp-fac2)
            y(2)   = 0.0e0_wp+r2*sin(th2)*cos(tau)*(1.0e0_wp-fac2)
            z(2)   = 0.0e0_wp+r2*sin(th2)*sin(tau)*(1.0e0_wp-fac2)
            r12    = sqrt((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
            r23    = sqrt((x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2)
            r14    = sqrt((x(4)-x(1))**2+(y(4)-y(1))**2+(z(4)-z(1))**2)
            r34    = sqrt((x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2)
            do j=1, 4
                cs(j) = cos(j*tau)
            enddo
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
            p1(2,0) = ang_pf(2,0)*0.5e0_wp*(3.0e0_wp*c2-1.0e0_wp)
            p1(3,0) = ang_pf(3,0)*0.5e0_wp*(5.0e0_wp*c3-3.0e0_wp*c)
            p1(4,0) = ang_pf(4,0)*0.125e0_wp*(35.0e0_wp*c4-30.0e0_wp*c2+3.0e0_wp)
            p1(1,1) = ang_pf(1,1)*s
            p1(2,1) = ang_pf(2,1)*3*s*c
            p1(3,1) = ang_pf(3,1)*1.5e0_wp*s*(5.0e0_wp*c2-1.0e0_wp)
            p1(4,1) = ang_pf(4,1)*2.5e0_wp*s*(7.0e0_wp*c3-3.0e0_wp*c)
            p1(2,2) = ang_pf(2,2)*3.0e0_wp*s2
            p1(3,2) = ang_pf(3,2)*15.0e0_wp*s2*c
            p1(4,2) = ang_pf(4,2)*7.5e0_wp*s2*(7.0e0_wp*c2-1.0e0_wp)
            p1(3,3) = ang_pf(3,3)*15.0e0_wp*s3
            p1(4,3) = ang_pf(4,3)*105.0e0_wp*s3*c
            p1(4,4) = ang_pf(4,4)*105.0e0_wp*s4
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
            p2(2,0) = ang_pf(2,0)*0.5e0_wp*(3.0e0_wp*c2-1.0e0_wp)
            p2(3,0) = ang_pf(3,0)*0.5e0_wp*(5.0e0_wp*c3-3.0e0_wp*c)
            p2(4,0) = ang_pf(4,0)*0.125e0_wp*(35.0e0_wp*c4-30.0e0_wp*c2+3.0e0_wp)
            p2(1,1) = ang_pf(1,1)*s
            p2(2,1) = ang_pf(2,1)*3.0e0_wp*s*c
            p2(3,1) = ang_pf(3,1)*1.5e0_wp*s*(5.0e0_wp*c2-1.0e0_wp)
            p2(4,1) = ang_pf(4,1)*2.5e0_wp*s*(7.0e0_wp*c3-3.0e0_wp*c)
            p2(2,2) = ang_pf(2,2)*3.0e0_wp*s2
            p2(3,2) = ang_pf(3,2)*15.0e0_wp*s2*c
            p2(4,2) = ang_pf(4,2)*7.5e0_wp*s2*(7.0e0_wp*c2-1.0e0_wp)
            p2(3,3) = ang_pf(3,3)*15.0e0_wp*s3
            p2(4,3) = ang_pf(4,3)*105.0e0_wp*s3*c
            p2(4,4) = ang_pf(4,4)*105.0e0_wp*s4
            icu     = 31
            do j1=0, 4, 1
                do j2=j1, 4, 1
                    do j3=(j2-j1), (j1+j2), 2
                        potcm_ac(1,j1,j2,j3) = potcm_a(icu)
                        potcm_ac(2,j1,j2,j3) = potcm_a(icu+35)
                        potcm_ac(3,j1,j2,j3) = potcm_a(icu+70)
                        icu                  = icu+1
                    enddo
                enddo
            enddo
            call gsrsqs
            do j1=0, 4, 1
                do j2=0, 4, 1
                    do j3=abs(j2-j1), (j2+j1), 2
                        vp(j1, j2, j3) = ang_g(j1,j2,j3,0)*p1(j1,0)*p2(j2,0)
                        do m=1, min(j1, j2), 1
                            vp(j1, j2, j3) = vp(j1, j2, j3)+(-1)**m*2*ang_g(j1,j2,j3,m)* cs(m)*p1(j1,m)*p2(j2,m)
                        enddo
                        vp(j1, j2, j3) = vp(j1, j2, j3)*(-1)**(j1-j2)*(2*j3+1)*pif
                    enddo
                enddo
            enddo
            e = 0.0e0_wp
            do j1=0, 4, 1
                do j2=0, 4, 1
                    do j3=abs(j2-j1), (j2+j1), 2
                        e = e+vp(j1, j2, j3)*potcm_ar(j1,j2,j3)
                    enddo
                enddo
            enddo
            if( rcm .lt. 7.2e0_wp )then
                e = e+vp(0, 0, 0)*exp(-((7.2e0_wp-rcm)/rcm)**2)*potcm_a(2)/rcm**6
            else
                e = e+vp(0, 0, 0)*potcm_a(2)/rcm**6
            endif
            if( rcm .lt. 8.4e0_wp )then
                e = e+vp(0, 0, 0)*exp(-((8.4e0_wp-rcm)/rcm)**2)*potcm_a(3)/rcm**8
            else
                e = e+vp(0, 0, 0)*potcm_a(3)/rcm**8
            endif
            if( rcm .lt. 9.6e0_wp )then
                e = e+vp(0, 0, 0)*exp(-((9.6e0_wp-rcm)/rcm)**2)*potcm_a(4)/rcm**10
            else
                e = e+vp(0, 0, 0)*potcm_a(4)/rcm**10
            endif
            if( rcm .lt. 9.5e0_wp )then
                e = e+vp(1, 1, 2)*potcm_a(5)/rcm**3*exp(-((9.5e0_wp-rcm)/rcm)**2)
            else
                e = e+vp(1, 1, 2)*potcm_a(5)/rcm**3
            endif
            if( rcm .lt. 10.5e0_wp )then
                e = e+(vp(2, 1, 3)-vp(1, 2, 3))*potcm_a(6)/rcm**4* exp(-((10.5e0_wp-rcm)/rcm)**2)
            else
                e = e+(vp(2, 1, 3)-vp(1, 2, 3))*potcm_a(6)/rcm**4
            endif
            if( rcm .lt. 11.5e0_wp )then
                e = e+vp(2, 2, 4)*potcm_a(7)/rcm**5*exp(-((11.5e0_wp-rcm)/rcm)**2)
                e = e+(vp(3, 1, 4)+vp(1, 3, 4))*potcm_a(8)/rcm**5*exp(-((11.5-rcm)/rcm)**2)
            else
                e = e+vp(2, 2, 4)*potcm_a(7)/rcm**5
                e = e+(vp(3, 1, 4)+vp(1, 3, 4))*potcm_a(8)/rcm**5
            endif
            e = e+potcm_a(9)*exp(-0.5e0_wp*r12)
            e = e+potcm_a(10)*exp(-1.5e0_wp*r12)
            e = e+potcm_a(13)*exp(-1.5e0_wp*r34)
            e = e+potcm_a(17)*(exp(-10.0e0_wp*r14)+exp(-10.0e0_wp*r23))
            if( rcm .lt. 4.5e0_wp )then
                e = e+potcm_a(18)/rcm**6*((3.0e0_wp*(cos(th1))**2+1.0e0_wp)+ (3.0e0_wp*(cos(th2))**2+1.0e0_wp))* exp(- &
                    & ((4.5e0_wp- rcm)/rcm)**2)
            else
                e = e+potcm_a(18)/rcm**6*((3.0e0_wp*(cos(th1))**2+1.0e0_wp)+ (3.0e0_wp*(cos(th2))**2+1.0e0_wp))
            endif
            ! monomer terms
            e=e+47635e0_wp*(1.0e0_wp-potcm_a(19)*exp(-0.1e0_wp*(r23-3.0e0_wp)**2))* (1.0e0_wp-exp(-1.1953e0_wp* &
                & (1.0e0_wp-potcm_a(26)* (1.0e0_wp-tanh(potcm_a(27)*(r14- 3.0e0_wp))) )* (r1-1.7374e0_wp*(1.0e0_wp - &
                & potcm_a(20)* (1.0e0_wp-tanh(potcm_a(21)* (r14-3.0e0_wp))) ))))**2
            e=e+47635e0_wp*(1.0e0_wp-potcm_a(19)*exp(-0.1e0_wp*(r14-3.0e0_wp)**2))* (1.0e0_wp-exp(-1.1953e0_wp* &
                & (1.0e0_wp-potcm_a(26)* (1.0e0_wp-tanh(potcm_a(27)*(r23- 3.0e0_wp))) )* (r2-1.7374e0_wp*(1.0e0_wp - &
                & potcm_a(20)* (1.0e0_wp-tanh(potcm_a(21)* (r23-3.0e0_wp))) ))))**2
            !     next line corrects for asymptotic diatomic potential which is
            !     added in separately.
            e = e-(ptntl(r1-1.7374e0_wp)+ptntl(r2-1.7374e0_wp))*219474.63067e0_wp
            e = e+potcm_a(28)*(exp(-4.0e0_wp*(r14-r1)**2)+exp(-4.0e0_wp*(r23-r2)**2))
            e = e+potcm_a(22)*exp(-4.0e0_wp*r14)*exp(-4.0e0_wp*r23)
            e = e+potcm_a(23)
            e = e+potcm_a(24)*(exp(-r14)+exp(-r23))* (exp(-15.0e0_wp*r1)+exp(-15.0e0_wp*r2))
            !   Units conversion:  Conversion to atomic units
            e                   = e/219474.63067e0_wp
            ground_state_energy = e
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate PTNTL for the F2H2_GEN_SQSBDE_1995 potential.
        function ptntl( x )
            !     HF DIATOMIC POTENTIAL CURVE FOR SQSBDE SURFACE, CALCULATED USING
            !     PARAMETERS GIVEN IN JCP, VOL. 95(5), P. 31.
            !     NOTICE THAT THE PARAMETER DM IS GIVEN IN ATOMIC UNITS HERE AND
            !     IN CM-1 IN THE ICP REFERENCE.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: ptntl, x
            real(wp), parameter :: am  = 1.1953e0_wp
            real(wp), parameter :: dm  = 0.21704e0_wp
            real(wp), parameter :: one = 1.0e0_wp
            !      SAVE
            ptntl = dm*(one-exp(-am*x))**2
        end function ptntl

        !> @brief Evaluate POTINFO for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine potinfo
            implicit none
            integer            :: i, j, k, kmax
            integer, parameter :: isurf = 5
            do i =1, 5
            enddo
            kmax = 0
            do i = 1, isurf+1
                do j = 1, isurf+1
                    if( surface_requests(i,j) .ne. 0.and.kmax .lt. max(i, j) ) kmax = max(i, j)
                enddo
            enddo
            if( kmax-1 .gt. maximum_excited_state )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            first_requested_excited_state = 0
            last_requested_excited_state  = 0
            do i = 2, isurf+1
                if( surface_requests(i,i) .ne. 0 )then
                    last_requested_excited_state = i-1
                    if( first_requested_excited_state .eq. 0 ) first_requested_excited_state = i-1
                endif
            enddo
            first_requested_coupling = 0
            last_requested_coupling  = 0
            k                        = 0
            do i = 1, isurf
                do j = i+1, isurf+1
                    k = k+1
                    if( surface_requests(i,j)+surface_requests(j,i) .ne. 0 )then
                        last_requested_coupling = k
                        if( first_requested_coupling .eq. 0 ) first_requested_coupling = k
                    endif
                enddo
            enddo
            if( requested_derivative_order .gt. maximum_derivative_order )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine potinfo

        !> @brief Evaluate CARTTOR for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(cartesian_coordinates, atom_count, coordinate_representation, pt1cm_r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine eunitzero
            implicit none
            integer :: i, j
            output_ground_state_energy = ground_state_energy
            if( first_requested_excited_state .ne. 0 )then
                do i=first_requested_excited_state, last_requested_excited_state
                    output_excited_state_energies(i) = excited_state_energies(i)
                enddo
            endif
            if( first_requested_coupling .ne. 0 )then
                do j=first_requested_coupling, last_requested_coupling
                    output_diabatic_couplings(j) = diabatic_couplings(j)
                enddo
            endif
        end subroutine eunitzero

        !> @brief Evaluate RTOCART for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine rtocart
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            real(wp)           :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, &
                & termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes(n3atom, isurf), &
                & ygs(n3atom), yij(n3atom, jsurf), za, zb, zc, zd
            integer :: i, ind, j, j1, j2, k, m1
            if( coordinate_representation .eq. 1 )then
                do i = 1, atom_count
                    ind                                  = 3*i-2
                    ground_state_cartesian_gradient(i,1) = ground_state_internal_gradient(ind)
                    ground_state_cartesian_gradient(i,2) = ground_state_internal_gradient(ind+1)
                    ground_state_cartesian_gradient(i,3) = ground_state_internal_gradient(ind+2)
                    if( first_requested_excited_state .ne. 0 )then
                        do j = first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(i,1,j) = excited_state_internal_gradients(ind,j)
                            excited_state_cartesian_gradients(i,2,j) = excited_state_internal_gradients(ind+1,j)
                            excited_state_cartesian_gradients(i,3,j) = excited_state_internal_gradients(ind+2,j)
                        enddo
                    endif
                    if( last_requested_coupling .ne. 0 )then
                        do k = first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(i,1,k) = diabatic_internal_gradients(ind,k)
                            diabatic_cartesian_gradients(i,2,k) = diabatic_internal_gradients(ind+1,k)
                            diabatic_cartesian_gradients(i,3,k) = diabatic_internal_gradients(ind+2,k)
                        enddo
                    endif
                enddo
            elseif( coordinate_representation .eq. 2 )then
                do i = 1, atom_count
                    ground_state_cartesian_gradient(i,1) = 0.e0_wp
                    ground_state_cartesian_gradient(i,2) = 0.e0_wp
                    ground_state_cartesian_gradient(i,3) = 0.e0_wp
                    if( first_requested_excited_state .ne. 0 )then
                        do j1=first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(i,1,j1) = 0.e0_wp
                            excited_state_cartesian_gradients(i,2,j1) = 0.e0_wp
                            excited_state_cartesian_gradients(i,3,j1) = 0.e0_wp
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do j2=first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(i,1,j2) = 0.e0_wp
                            diabatic_cartesian_gradients(i,2,j2) = 0.e0_wp
                            diabatic_cartesian_gradients(i,3,j2) = 0.e0_wp
                        enddo
                    endif
                    do j = 1, atom_count
                        if( j .lt. i )then
                            m1 = atom_count*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
                            m1 = atom_count*(i-1) - (i*(i-1))/2 + j-i
                        else
                            cycle
                        endif
                        y = ground_state_internal_gradient(m1)
                        termx                                = (cartesian_coordinates(i,1)-cartesian_coordinates(j, &
                            & 1))/pt1cm_r(m1)
                        termy                                = (cartesian_coordinates(i,2)-cartesian_coordinates(j, &
                            & 2))/pt1cm_r(m1)
                        termz                                = (cartesian_coordinates(i,3)-cartesian_coordinates(j, &
                            & 3))/pt1cm_r(m1)
                        ground_state_cartesian_gradient(i,1) = ground_state_cartesian_gradient(i,1) + termx*y
                        ground_state_cartesian_gradient(i,2) = ground_state_cartesian_gradient(i,2) + termy*y
                        ground_state_cartesian_gradient(i,3) = ground_state_cartesian_gradient(i,3) + termz*y
                        if( first_requested_excited_state .gt. 0 )then
                            y = excited_state_internal_gradients(m1,j1)
                            do j1=first_requested_excited_state, last_requested_excited_state
                                excited_state_cartesian_gradients(i,1,j1) = excited_state_cartesian_gradients(i,1, &
                                    & j1) + termx*y
                                excited_state_cartesian_gradients(i,2,j1) = excited_state_cartesian_gradients(i,2, &
                                    & j1) + termy*y
                                excited_state_cartesian_gradients(i,3,j1) = excited_state_cartesian_gradients(i,3, &
                                    & j1) + termz*y
                            enddo
                        elseif( first_requested_coupling .gt. 0 )then
                            do j2=first_requested_coupling, last_requested_coupling
                                y                                    = diabatic_internal_gradients(m1,j2)
                                diabatic_cartesian_gradients(i,1,j2) = diabatic_cartesian_gradients(i,1,j2) + termx*y
                                diabatic_cartesian_gradients(i,2,j2) = diabatic_cartesian_gradients(i,2,j2) + termy*y
                                diabatic_cartesian_gradients(i,3,j2) = diabatic_cartesian_gradients(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            elseif( coordinate_representation .eq. 3 )then
                do i = 1, atom_count
                    ygs(i) = ground_state_internal_gradient(i)/pt1cm_r(i)
                    if( first_requested_excited_state .ne. 0 )then
                        do j=first_requested_excited_state, last_requested_excited_state
                            yes(i, j) = excited_state_internal_gradients(i,j)/pt1cm_r(i)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do k=first_requested_coupling, last_requested_coupling
                            yij(i, k) = diabatic_internal_gradients(i,k)/pt1cm_r(i)
                        enddo
                    endif
                enddo
                do k = 1, 3
                    term12                               = cartesian_coordinates(1,k)-cartesian_coordinates(2,k)
                    term23                               = cartesian_coordinates(2,k)-cartesian_coordinates(3,k)
                    term13                               = cartesian_coordinates(1,k)-cartesian_coordinates(3,k)
                    ground_state_cartesian_gradient(1,k) = term12*ygs(1) + term13*ygs(3)
                    ground_state_cartesian_gradient(2,k) = -term12*ygs(1) + term23*ygs(2)
                    ground_state_cartesian_gradient(3,k) = -term13*ygs(3) - term23*ygs(2)
                    if( first_requested_excited_state .ne. 0 )then
                        do j1=first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(1,k,j1) = term12*yes(1, j1) + term13*yes(3, j1)
                            excited_state_cartesian_gradients(2,k,j1) = -term12*yes(1, j1) + term23*yes(2, j1)
                            excited_state_cartesian_gradients(3,k,j1) = -term13*yes(3, j1) - term23*yes(2, j1)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do j2=first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(1,k,j2) = term12*yij(1, j2) + term13*yij(3, j2)
                            diabatic_cartesian_gradients(2,k,j2) = -term12*yij(1, j2) + term23*yij(2, j2)
                            diabatic_cartesian_gradients(3,k,j2) = -term13*yij(3, j2) - term23*yij(2, j2)
                        enddo
                    endif
                enddo
            elseif( coordinate_representation .eq. 4 )then
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
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the F2H2_GEN_SQSBDE_1995 potential.
        subroutine dedcou
            implicit none
            integer :: i, j, k
            do i = 1, atom_count
                output_ground_state_gradient(i,1) = ground_state_cartesian_gradient(i,1)
                output_ground_state_gradient(i,2) = ground_state_cartesian_gradient(i,2)
                output_ground_state_gradient(i,3) = ground_state_cartesian_gradient(i,3)
                if( first_requested_excited_state .ne. 0 )then
                    do j=first_requested_excited_state, last_requested_excited_state
                        output_excited_state_gradients(i,1,j) = excited_state_cartesian_gradients(i,1,j)
                        output_excited_state_gradients(i,2,j) = excited_state_cartesian_gradients(i,2,j)
                        output_excited_state_gradients(i,3,j) = excited_state_cartesian_gradients(i,3,j)
                    enddo
                endif
                if( first_requested_coupling .ne. 0 )then
                    do k=first_requested_coupling, last_requested_coupling
                        output_diabatic_gradients(i,1,k) = diabatic_cartesian_gradients(i,1,k)
                        output_diabatic_gradients(i,2,k) = diabatic_cartesian_gradients(i,2,k)
                        output_diabatic_gradients(i,3,k) = diabatic_cartesian_gradients(i,3,k)
                    enddo
                endif
            enddo
        end subroutine dedcou
    end subroutine evaluate
end module chempotpy_f2h2_gen_sqsbde_1995_pes

!> @brief Expose F2H2_GEN_SQSBDE_1995 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_f2h2_gen_sqsbde_1995_pes, only: evaluate
    use chempotpy_kinds,                    only: wp
    implicit none
    integer, parameter    :: natoms  = 4
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
