!   System: CH3OH + Cl --> CH2OH + HCl
!                         --> CH3O + HCl
!              Cipriano Rangel and Joaquin Espinosa-Garcia 2023
!   All the information passed to and from the potential energy surface
!   routine is in hartree atomic units.
!        This potential is written such that:
!                       X(1)  - X(3)  : X, Y, Z for C
!                       X(4)  - X(6)  : X, Y, Z for O
!                       X(7)  - X(9)  : X, Y, Z for H1
!                       X(10) - X(12) : X, Y, Z for H2
!                       X(13) - X(15) : X, Y, Z for H3
!                       X(16) - X(18) : X, Y, Z for H4
!                       X(19) - X(21) : X, Y, Z for Cl

!> @brief Implement CH4OCl_VBMM_2023, a potential energy surface for CH4OCl.
module chempotpy_ch4ocl_vbmm_2023_pes
    use chempotpy_geometry,  only: cross_product_3d
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the CH4OCl_VBMM_2023 potential energy surface for CH4OCl.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: hartree_to_electronvolt
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 7
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter    :: n3atom  = 75
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        integer               :: coord2_n1h(3,3)
        integer               :: coord2_n2h(1,3)
        integer               :: coord2_nc(3)
        integer               :: coord2_nc1(3)
        integer               :: coord2_nc2(3)
        integer               :: coord2_nhb(3)
        real(wp)              :: coord3_tc1b(3)
        real(wp)              :: coord3_tc2b(3)
        real(wp)              :: coord3_tcc(3)
        real(wp)              :: coord4_tb1h(3,3)
        real(wp)              :: coord4_tb2h(1,3)
        real(wp)              :: coord4_tc1h(3,3)
        real(wp)              :: coord4_tc2h(1,3)
        real(wp)              :: coord5_rb1h(3)
        real(wp)              :: coord5_rb2h(1)
        real(wp)              :: coord5_rc1b
        real(wp)              :: coord5_rc1h(3)
        real(wp)              :: coord5_rc2b
        real(wp)              :: coord5_rc2h(1)
        real(wp)              :: coord5_rcc
        real(wp)              :: coord6_r0c1h
        real(wp)              :: coord6_r0c2h
        real(wp)              :: coord6_r0cc
        real(wp)              :: cartesian_coordinates(natom,3)
        integer               :: coordinate_representation
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: atom_count
        character(len=75)     :: infocm_ref(5)
        real(wp)              :: ipforce1_f1(6)
        real(wp)              :: ipforce1_fk0(6,6)
        integer               :: potcm2_nc1h(3)
        integer               :: potcm2_nc2h(1)
        real(wp)              :: potcm_a1c
        real(wp)              :: potcm_a1c1h
        real(wp)              :: potcm_a1c2h
        real(wp)              :: potcm_a1s1
        real(wp)              :: potcm_a1s2
        real(wp)              :: potcm_a1scc
        real(wp)              :: potcm_a2c
        real(wp)              :: potcm_a2s1
        real(wp)              :: potcm_a2s2
        real(wp)              :: potcm_a2scc
        real(wp)              :: potcm_a3s
        real(wp)              :: potcm_a3scc
        real(wp)              :: potcm_aa11
        real(wp)              :: potcm_aa12
        real(wp)              :: potcm_aa21
        real(wp)              :: potcm_aa22
        real(wp)              :: potcm_aa31
        real(wp)              :: potcm_aa32
        real(wp)              :: potcm_aa41
        real(wp)              :: potcm_aa42
        real(wp)              :: potcm_ac1b
        real(wp)              :: potcm_ac2b
        real(wp)              :: potcm_acc
        real(wp)              :: potcm_ah1h
        real(wp)              :: potcm_ah2h
        real(wp)              :: potcm_ak1
        real(wp)              :: potcm_aphi
        real(wp)              :: potcm_aphi1
        real(wp)              :: potcm_aphi2
        real(wp)              :: potcm_atheta
        real(wp)              :: potcm_athetacc
        real(wp)              :: potcm_b1c1h
        real(wp)              :: potcm_b1c2h
        real(wp)              :: potcm_b1s1
        real(wp)              :: potcm_b1s2
        real(wp)              :: potcm_b1scc
        real(wp)              :: potcm_b2s1
        real(wp)              :: potcm_b2s2
        real(wp)              :: potcm_b2scc
        real(wp)              :: potcm_b3s
        real(wp)              :: potcm_b3scc
        real(wp)              :: potcm_bphi
        real(wp)              :: potcm_bphi1
        real(wp)              :: potcm_bphi2
        real(wp)              :: potcm_btheta
        real(wp)              :: potcm_bthetacc
        real(wp)              :: potcm_c1c1h
        real(wp)              :: potcm_c1c2h
        real(wp)              :: potcm_cphi
        real(wp)              :: potcm_cphi1
        real(wp)              :: potcm_cphi2
        real(wp)              :: potcm_ctheta
        real(wp)              :: potcm_cthetacc
        real(wp)              :: potcm_d1c1b
        real(wp)              :: potcm_d1c1h
        real(wp)              :: potcm_d1c2b
        real(wp)              :: potcm_d1c2h
        real(wp)              :: potcm_d1cc
        real(wp)              :: potcm_d1h1h
        real(wp)              :: potcm_d1h2h
        real(wp)              :: potcm_d3c1b
        real(wp)              :: potcm_d3c1h
        real(wp)              :: potcm_d3c2b
        real(wp)              :: potcm_d3c2h
        real(wp)              :: potcm_d3h1h
        real(wp)              :: potcm_d3h2h
        real(wp)              :: potcm_fch3
        real(wp)              :: potcm_fkinf1
        real(wp)              :: potcm_fkinf2
        real(wp)              :: potcm_hch3
        integer               :: potcm_nnb
        integer               :: potcm_nnc1
        integer               :: potcm_nnc2
        real(wp)              :: potcm_r0c1b
        real(wp)              :: potcm_r0c1hp
        real(wp)              :: potcm_r0c1hr
        real(wp)              :: potcm_r0c2b
        real(wp)              :: potcm_r0c2hp
        real(wp)              :: potcm_r0c2hr
        real(wp)              :: potcm_r0ccp
        real(wp)              :: potcm_r0ccr
        real(wp)              :: potcm_r0h1h
        real(wp)              :: potcm_r0h2h
        real(wp)              :: potcm_v3
        real(wp)              :: potcm_w1
        real(wp)              :: potcm_w10
        real(wp)              :: potcm_w2
        real(wp)              :: potcm_w3
        real(wp)              :: potcm_w4
        real(wp)              :: potcm_w5
        real(wp)              :: potcm_w6
        real(wp)              :: potcm_w7
        real(wp)              :: potcm_w8
        real(wp)              :: potcm_w9
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: ground_state_energy
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: qpdot_pdot(27)
        real(wp)              :: qpdot_q(27)
        real(wp)              :: refangles1_sphi(6)
        real(wp)              :: refangles1_stheta(6)
        real(wp)              :: refangles1_theta0(6,6)
        real(wp)              :: switch1_s1(6)
        real(wp)              :: switch1_s2(6)
        real(wp)              :: switch1_s3(4)
        real(wp)              :: torsion1_t2
        real(wp)              :: torsion1_thetahcch(3)
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
        real(wp)              :: vop6_fdelta(4)
        real(wp)              :: vop6_hdelta(4)
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        coord2_n1h                        = 0; coord2_n2h = 0; coord2_nc = 0; coord2_nc1 = 0; coord2_nc2 = 0; coord2_nhb = 0
        coord3_tc1b                       = 0.e0_wp; coord3_tc2b = 0.e0_wp; coord3_tcc = 0.e0_wp; coord4_tb1h = 0.e0_wp; coord4_tb2h = 0.e0_wp
        coord4_tc1h                       = 0.e0_wp; coord4_tc2h = 0.e0_wp; coord5_rb1h = 0.e0_wp; coord5_rb2h = 0.e0_wp; coord5_rc1b = 0.e0_wp
        coord5_rc1h                       = 0.e0_wp; coord5_rc2b = 0.e0_wp; coord5_rc2h = 0.e0_wp; coord5_rcc = 0.e0_wp; coord6_r0c1h = 0.e0_wp
        coord6_r0c2h                      = 0.e0_wp; coord6_r0cc = 0.e0_wp; cartesian_coordinates = 0.e0_wp; coordinate_representation = 0
        maximum_derivative_order          = 0; maximum_excited_state = 0; atom_count = 0
        infocm_ref                        = ' '
        ipforce1_f1                       = 0.e0_wp
        ipforce1_fk0                      = 0.e0_wp
        potcm2_nc1h                       = 0
        potcm2_nc2h                       = 0
        potcm_a1c                         = 0.e0_wp
        potcm_a1c1h                       = 0.e0_wp
        potcm_a1c2h                       = 0.e0_wp
        potcm_a1s1                        = 0.e0_wp
        potcm_a1s2                        = 0.e0_wp
        potcm_a1scc                       = 0.e0_wp
        potcm_a2c                         = 0.e0_wp
        potcm_a2s1                        = 0.e0_wp
        potcm_a2s2                        = 0.e0_wp
        potcm_a2scc                       = 0.e0_wp
        potcm_a3s                         = 0.e0_wp
        potcm_a3scc                       = 0.e0_wp
        potcm_aa11                        = 0.e0_wp
        potcm_aa12                        = 0.e0_wp
        potcm_aa21                        = 0.e0_wp
        potcm_aa22                        = 0.e0_wp
        potcm_aa31                        = 0.e0_wp
        potcm_aa32                        = 0.e0_wp
        potcm_aa41                        = 0.e0_wp
        potcm_aa42                        = 0.e0_wp
        potcm_ac1b                        = 0.e0_wp
        potcm_ac2b                        = 0.e0_wp
        potcm_acc                         = 0.e0_wp
        potcm_ah1h                        = 0.e0_wp
        potcm_ah2h                        = 0.e0_wp
        potcm_ak1                         = 0.e0_wp
        potcm_aphi                        = 0.e0_wp
        potcm_aphi1                       = 0.e0_wp
        potcm_aphi2                       = 0.e0_wp
        potcm_atheta                      = 0.e0_wp
        potcm_athetacc                    = 0.e0_wp
        potcm_b1c1h                       = 0.e0_wp
        potcm_b1c2h                       = 0.e0_wp
        potcm_b1s1                        = 0.e0_wp
        potcm_b1s2                        = 0.e0_wp
        potcm_b1scc                       = 0.e0_wp
        potcm_b2s1                        = 0.e0_wp
        potcm_b2s2                        = 0.e0_wp
        potcm_b2scc                       = 0.e0_wp
        potcm_b3s                         = 0.e0_wp
        potcm_b3scc                       = 0.e0_wp
        potcm_bphi                        = 0.e0_wp
        potcm_bphi1                       = 0.e0_wp
        potcm_bphi2                       = 0.e0_wp
        potcm_btheta                      = 0.e0_wp
        potcm_bthetacc                    = 0.e0_wp
        potcm_c1c1h                       = 0.e0_wp
        potcm_c1c2h                       = 0.e0_wp
        potcm_cphi                        = 0.e0_wp
        potcm_cphi1                       = 0.e0_wp
        potcm_cphi2                       = 0.e0_wp
        potcm_ctheta                      = 0.e0_wp
        potcm_cthetacc                    = 0.e0_wp
        potcm_d1c1b                       = 0.e0_wp
        potcm_d1c1h                       = 0.e0_wp
        potcm_d1c2b                       = 0.e0_wp
        potcm_d1c2h                       = 0.e0_wp
        potcm_d1cc                        = 0.e0_wp
        potcm_d1h1h                       = 0.e0_wp
        potcm_d1h2h                       = 0.e0_wp
        potcm_d3c1b                       = 0.e0_wp
        potcm_d3c1h                       = 0.e0_wp
        potcm_d3c2b                       = 0.e0_wp
        potcm_d3c2h                       = 0.e0_wp
        potcm_d3h1h                       = 0.e0_wp
        potcm_d3h2h                       = 0.e0_wp
        potcm_fch3                        = 0.e0_wp
        potcm_fkinf1                      = 0.e0_wp
        potcm_fkinf2                      = 0.e0_wp
        potcm_hch3                        = 0.e0_wp
        potcm_nnb                         = 0
        potcm_nnc1                        = 0
        potcm_nnc2                        = 0
        potcm_r0c1b                       = 0.e0_wp
        potcm_r0c1hp                      = 0.e0_wp
        potcm_r0c1hr                      = 0.e0_wp
        potcm_r0c2b                       = 0.e0_wp
        potcm_r0c2hp                      = 0.e0_wp
        potcm_r0c2hr                      = 0.e0_wp
        potcm_r0ccp                       = 0.e0_wp
        potcm_r0ccr                       = 0.e0_wp
        potcm_r0h1h                       = 0.e0_wp
        potcm_r0h2h                       = 0.e0_wp
        potcm_v3                          = 0.e0_wp
        potcm_w1                          = 0.e0_wp
        potcm_w10                         = 0.e0_wp
        potcm_w2                          = 0.e0_wp
        potcm_w3                          = 0.e0_wp
        potcm_w4                          = 0.e0_wp
        potcm_w5                          = 0.e0_wp
        potcm_w6                          = 0.e0_wp
        potcm_w7                          = 0.e0_wp
        potcm_w8                          = 0.e0_wp
        potcm_w9                          = 0.e0_wp
        ground_state_internal_gradient    = 0.e0_wp
        ground_state_energy               = 0.e0_wp
        pt1cm_r                           = 0.e0_wp
        pt3cm_ezero                       = 0.e0_wp
        excited_state_internal_gradients  = 0.e0_wp
        excited_state_energies            = 0.e0_wp
        diabatic_internal_gradients       = 0.e0_wp
        diabatic_couplings                = 0.e0_wp
        qpdot_pdot                        = 0.e0_wp
        qpdot_q                           = 0.e0_wp
        refangles1_sphi                   = 0.e0_wp
        refangles1_stheta                 = 0.e0_wp
        refangles1_theta0                 = 0.e0_wp
        switch1_s1                        = 0.e0_wp
        switch1_s2                        = 0.e0_wp
        switch1_s3                        = 0.e0_wp
        torsion1_t2                       = 0.e0_wp
        torsion1_thetahcch                = 0.e0_wp
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
        vop6_fdelta                       = 0.e0_wp
        vop6_hdelta                       = 0.e0_wp
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 1
        coordinate_representation  = 1
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 7
        potcm_nnc1                 = 1
        potcm_nnc2                 = 2
        potcm_nnb                  = 7
        potcm2_nc1h                = [3, 4, 5]
        potcm2_nc2h                = 6
        potcm_r0c1hr               = 1.09047e0_wp
        potcm_r0c1hp               = 1.08177e0_wp
        potcm_w3                   = 1.0000e0_wp
        potcm_w4                   = 1.09047e0_wp
        potcm_d1c1h                = 105.23000e0_wp
        potcm_d3c1h                = 27.30400e0_wp
        potcm_a1c1h                = 1.59300e0_wp
        potcm_b1c1h                = 0.52500e0_wp
        potcm_c1c1h                = 0.51404e0_wp
        potcm_r0c2hr               = 0.95897e0_wp
        potcm_r0c2hp               = 0.96077e0_wp
        potcm_w5                   = 1.000e0_wp
        potcm_w6                   = 0.96897e0_wp
        potcm_d1c2h                = 116.23000e0_wp
        potcm_d3c2h                = 37.00400e0_wp
        potcm_a1c2h                = 1.55300e0_wp
        potcm_b1c2h                = 1.42500e0_wp
        potcm_c1c2h                = 0.51404e0_wp
        potcm_r0ccr                = 1.42007e0_wp
        potcm_r0ccp                = 1.37500e0_wp
        potcm_w9                   = 1.5000e0_wp
        potcm_w10                  = 1.08507e0_wp
        potcm_d1cc                 = 109.94900e0_wp
        potcm_acc                  = 1.934e0_wp
        potcm_a1c                  = 1.20997e0_wp
        potcm_a2c                  = 1.949e0_wp
        potcm_r0h1h                = 1.29191e0_wp
        potcm_d1h1h                = 111.458e0_wp
        potcm_d3h1h                = 37.064e0_wp
        potcm_ah1h                 = 1.9157e0_wp
        potcm_r0h2h                = 1.29191e0_wp
        potcm_d1h2h                = 109.458e0_wp
        potcm_d3h2h                = 26.064e0_wp
        potcm_ah2h                 = 1.8057e0_wp
        potcm_r0c1b                = 2.05897e0_wp
        potcm_d1c1b                = 35.009e0_wp
        potcm_d3c1b                = 23.000e0_wp
        potcm_ac1b                 = 3.1090285e0_wp
        potcm_r0c2b                = 1.74897e0_wp
        potcm_d1c2b                = 68.609e0_wp
        potcm_d3c2b                = 17.400e0_wp
        potcm_ac2b                 = 1.0170285e0_wp
        potcm_aphi                 = 0.528790e0_wp
        potcm_bphi                 = 1.200664e0_wp
        potcm_cphi                 = 1.420074e0_wp
        potcm_aphi1                = 0.19700e0_wp
        potcm_bphi1                = 1.20600e0_wp
        potcm_cphi1                = 1.08200e0_wp
        potcm_aphi2                = 0.19700e0_wp
        potcm_bphi2                = 1.20600e0_wp
        potcm_cphi2                = 0.96800e0_wp
        potcm_a3s                  = 0.12770e0_wp
        potcm_b3s                  = 1.08200e0_wp
        potcm_a3scc                = 0.297000e0_wp
        potcm_b3scc                = 1.420664e0_wp
        potcm_atheta               = 1.65781e0_wp
        potcm_btheta               = 0.950600e0_wp
        potcm_ctheta               = 1.08201e0_wp
        potcm_athetacc             = 1.45781e0_wp
        potcm_bthetacc             = 0.950600e0_wp
        potcm_cthetacc             = 1.42001e0_wp
        potcm_fch3                 = 0.19770e0_wp
        potcm_hch3                 = 0.09600e0_wp
        potcm_fkinf1               = 0.44770e0_wp
        potcm_ak1                  = 0.12600e0_wp
        potcm_aa11                 = 0.15501e0_wp
        potcm_aa21                 = 1.59995e0_wp
        potcm_aa31                 = 2.16996e0_wp
        potcm_aa41                 = 11.56595e0_wp
        potcm_fkinf2               = 0.45770e0_wp
        potcm_aa12                 = 0.05001e0_wp
        potcm_aa22                 = 1.30095e0_wp
        potcm_aa32                 = 1.49996e0_wp
        potcm_aa42                 = 11.46595e0_wp
        potcm_v3                   = 0.03000e0_wp
        potcm_w1                   = 0.11000e0_wp
        potcm_w2                   = 1.08697e0_wp
        potcm_w7                   = 0.11000e0_wp
        potcm_w8                   = 1.08367e0_wp
        potcm_a1s1                 = 1.53137e0_wp
        potcm_b1s1                 = 1.09463e0_wp
        potcm_a2s1                 = 1.01474e0_wp
        potcm_b2s1                 = 0.96898e0_wp
        potcm_a1s2                 = 1.53137e0_wp
        potcm_b1s2                 = 1.09463e0_wp
        potcm_a2s2                 = 1.51474e0_wp
        potcm_b2s2                 = 0.96898e0_wp
        potcm_a1scc                = 1.53137e0_wp
        potcm_b1scc                = 1.42003e0_wp
        potcm_a2scc                = 1.51474e0_wp
        potcm_b2scc                = 1.42008e0_wp
        v                          = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        input_cartesian            = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                input_cartesian(iatom,idir) = x(iatom,idir)
            enddo
        enddo
        requested_derivative_order = igrad
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = output_ground_state_energy*hartree_to_electronvolt
            enddo
        else
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives are unavailable')
        endif

        contains

        !> @brief Evaluate POT for the CH4OCl_VBMM_2023 potential.
        subroutine pot
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: en, v, vip, vop, vstr, vtor
            integer  :: i, n
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            ! Start of changes by EMN.
            ! The lines below do not seem to be needed because the Cartesian
            ! coordinates are in common qpdot.
            ! Besides, BOHR_A in utility.f does not match the value below.
            ! As a consequence, q(i) changes and energy is not conserved
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !     PUT COORDINATES IN PROPER ARRAYS
            ! Inicialization de variables
            do i=1, 27
                qpdot_q(i)    = pt1cm_r(i)
                qpdot_pdot(i) = 0.e0_wp
            enddo
            ! Changing units from angstrom to bohr
            ! End of corrections by EMN
            !  calculate relative coordinates and bond lengths
            en = 0.0e0_wp
            call coorden
            call torsion(vtor)
            !  calculate switching functions
            call switchf
            !  calculate reference angles and their derivatives
            call refangles
            !  calculate stretching potential
            call stretch(vstr)
            !  calculate out of plane bending potential
            call opbend(vop)
            !  calculate in plane bending potential
            call ipbend(vip)
            !  total potential energy is vstr+vop+vip
            en = vstr+vop+vip+vtor
            !  convert from 10(5) j/mol to au
            en                  = en*0.03812e0_wp
            v                   = en
            ground_state_energy = en
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                !  Initialize derivatives
                do n=1, 21
                    ground_state_internal_gradient(n) = qpdot_pdot(n)*0.0201723e0_wp
                enddo
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate coorden for the CH4OCl_VBMM_2023 potential.
        subroutine coorden
            !  calculates relative coordinates and bond lengths
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: argmax, argp5, p5, rav4, t5tmp
            integer  :: i, ind
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            !  calculate relative coordinates
            do ind=1, 3
                coord3_tc1b(ind) = qpdot_q(coord2_nc1(ind))-qpdot_q(coord2_nhb(ind))
                coord3_tc2b(ind) = qpdot_q(coord2_nc2(ind))-qpdot_q(coord2_nhb(ind))
                coord3_tcc(ind)  = qpdot_q(coord2_nc1(ind))-qpdot_q(coord2_nc2(ind))
            enddo
            do ind=1, 3
                do i=1, 3
                    coord4_tc1h(i,ind) = qpdot_q(coord2_nc1(ind))-qpdot_q(coord2_n1h(i,ind))
                enddo
            enddo
            do ind=1, 3
                coord4_tc2h(1,ind) = qpdot_q(coord2_nc2(ind))-qpdot_q(coord2_n2h(1,ind))
            enddo
            do ind=1, 3
                do i=1, 3
                    coord4_tb1h(i,ind) = qpdot_q(coord2_nhb(ind))-qpdot_q(coord2_n1h(i,ind))
                enddo
                coord4_tb2h(1,ind) = qpdot_q(coord2_nhb(ind))-qpdot_q(coord2_n2h(1,ind))
            enddo
            coord5_rc1b = sqrt(coord3_tc1b(1)*coord3_tc1b(1)+coord3_tc1b(2)*coord3_tc1b(2)+coord3_tc1b(3)* &
                & coord3_tc1b(3))
            coord5_rc2b = sqrt(coord3_tc2b(1)*coord3_tc2b(1)+coord3_tc2b(2)*coord3_tc2b(2)+coord3_tc2b(3)* &
                & coord3_tc2b(3))
            coord5_rcc = sqrt(coord3_tcc(1)*coord3_tcc(1)+coord3_tcc(2)*coord3_tcc(2)+coord3_tcc(3)*coord3_tcc(3))
            do i=1, 3
                coord5_rc1h(i) = sqrt(coord4_tc1h(i,1)*coord4_tc1h(i,1)+coord4_tc1h(i,2)*coord4_tc1h(i,2)+ &
                    & coord4_tc1h(i,3)*coord4_tc1h(i,3))
            enddo
            coord5_rc2h(1) = sqrt(coord4_tc2h(1,1)*coord4_tc2h(1,1)+coord4_tc2h(1,2)*coord4_tc2h(1,2)+ coord4_tc2h(1, &
                & 3)*coord4_tc2h(1,3))
            do i=1, 3
                coord5_rb1h(i) = sqrt(coord4_tb1h(i,1)*coord4_tb1h(i,1)+coord4_tb1h(i,2)*coord4_tb1h(i,2)+ &
                    & coord4_tb1h(i,3)*coord4_tb1h(i,3))
            enddo
            coord5_rb2h(1) = sqrt(coord4_tb2h(1,1)*coord4_tb2h(1,1)+coord4_tb2h(1,2)*coord4_tb2h(1,2)+ coord4_tb2h(1, &
                & 3)*coord4_tb2h(1,3))
            ! crd 2013 modificaciones de las distancias
            ! CRD 2021 modificacion para hacer solo un parametro que
            ! modifique todos los r0c en funcion de rav
            argmax = 19.e0_wp
            rav4   = (coord5_rc1h(1)+coord5_rc1h(2)+coord5_rc1h(3)+coord5_rc2h(1))/4.0e0_wp
            p5     = 1.e0_wp
            argp5  = (potcm_w9*(rav4-potcm_w10))
            if( argp5 .lt. argmax )then
                t5tmp = 1-tanh(argp5)
            else
                t5tmp = 0.e0_wp
            endif
            p5           = p5*t5tmp
            coord6_r0cc  = p5*potcm_r0ccr+(1-p5)*potcm_r0ccp
            coord6_r0c1h = p5*potcm_r0c1hr+(1-p5)*potcm_r0c1hp
            coord6_r0c2h = p5*potcm_r0c2hr+(1-p5)*potcm_r0c2hp
            ! Fin CRD
            !c Dejamos las distancias fijas
            !c Dejamos las distancias fijas
        end subroutine coorden

        !> @brief Evaluate torsion for the CH4OCl_VBMM_2023 potential.
        subroutine torsion( vtor )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a, arg1, arg2, dihed, t1(3), vtor
            integer  :: i
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            !      pi=4.0d0*atan(1.0d0)
            ! CRD 2011 incluimos un apartado para calcular los angulos dihedro
            ! para poder calcular despues vtorsion
            call dihedro(3, 1, 2, 6, dihed)
            torsion1_thetahcch(1) = dihed*PI/180.0e0_wp
            call dihedro(4, 1, 2, 6, dihed)
            torsion1_thetahcch(2) = dihed*PI/180.0e0_wp
            call dihedro(5, 1, 2, 6, dihed)
            torsion1_thetahcch(3) = dihed*PI/180.0e0_wp
            ! CRD 2011 calculamos los terminos switching t1 y t2
            do i=1, 3
                t1(i) = (1-tanh(potcm_w1*(coord5_rc1h(i)-coord6_r0c1h)))
            enddo
            torsion1_t2 = (1-tanh(potcm_w7*(coord5_rc2h(1)-coord6_r0c2h)))
            vtor        = 0.0e0_wp
            ! CRD 2021 incluimos terminos para pasar de 3 a 2 minimos
            arg1 = (coord5_rb1h(1)-potcm_r0h1h)*(coord5_rb1h(2)-potcm_r0h1h)*(coord5_rb1h(3)-potcm_r0h1h)* &
                & (coord5_rb2h(1)-potcm_r0h2h)
            arg2 = tanh(arg1)
            if( arg2 .lt. 0.0 )then
                a = 2.0e0_wp
            else
                a = 2.0e0_wp+1.e0_wp*arg2
            endif
            a    = 3.00e0_wp
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(a*torsion1_thetahcch(1)))*t1(1)*torsion1_t2
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(a*torsion1_thetahcch(2)))*t1(2)*torsion1_t2
            vtor = vtor+(potcm_v3/3.0e0_wp)*(1+cos(a*torsion1_thetahcch(3)))*t1(3)*torsion1_t2
        end subroutine torsion

        !> @brief Evaluate dihedro for the CH4OCl_VBMM_2023 potential.
        subroutine dihedro( i, j, k, l, dihed )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dihed
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            integer  :: i, j, k, l
            integer  :: m
            real(wp) :: rij(3), rjk(3), rkl(3)
            real(wp) :: zijjk(3), zjkkl(3)
            real(wp) :: arg1, arg2
            !        pi=4.0d0*atan(1.0d0)
            do m = 1, 3
                rij(m) = qpdot_q((i-1)*3+m)-qpdot_q((j-1)*3+m)
                rjk(m) = qpdot_q((j-1)*3+m)-qpdot_q((k-1)*3+m)
                rkl(m) = qpdot_q((k-1)*3+m)-qpdot_q((l-1)*3+m)
            enddo
            zjkkl = cross_product_3d(rjk, rkl)
            zijjk = cross_product_3d(rij, rjk)
            arg2  = norm2(rjk)*dot_product(rij, zjkkl)
            arg1  = dot_product(zijjk, zjkkl)
            if( abs(arg1) .lt. 1.e-6_wp )then
                if( arg2 .gt. 0.e0_wp ) dihed = 90.e0_wp
                if( arg2 .lt. 0.e0_wp ) dihed = -90.e0_wp
                return
            endif
            if( arg1 .gt. 0.e0_wp )then
                dihed = atan(arg2/arg1) *180.e0_wp/PI
                return
            endif
            if( arg2 .ge. 0.e0_wp )then
                dihed = 180.e0_wp + atan(arg2/arg1) *180.e0_wp/PI
                return
            endif
            if( arg2 .lt. 0.e0_wp )then
                dihed = -180.e0_wp + atan(arg2/arg1) *180.e0_wp/PI
                return
            endif
            dihed = dihed*PI/180.e0_wp
        end subroutine dihedro

        !> @brief Evaluate refangles for the CH4OCl_VBMM_2023 potential.
        subroutine refangles
            !  subroutine calculates reference angles for the "in-plane" potential
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: halfpi, tau, tauco, tauoc, twopi
            integer  :: i, j
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            tau = acos(-1.0e0_wp/3.0e0_wp)
            ! CRD 2012 taucc para el angulo H-C-C
            tauco  = acos(-0.3616e0_wp)
            tauoc  = acos(-0.4558e0_wp)
            halfpi = 0.6555e0_wp*PI
            twopi  = 2.0e0_wp*PI
            ! jcc-2010
            ! valor para ajustar el angulo
            ! jcc-2010
            !  set diagonal elements to zero
            do i=1, 6
                refangles1_theta0(i,i) = 0.0e0_wp
            enddo
            !  calculate reference angles
            ! el atomo 4 es el c2
            refangles1_theta0(1,2) = tau+(tau-halfpi)*(refangles1_sphi(1)*refangles1_sphi(2)-1.0e0_wp) +(tau-twopi/ &
                & 3.0e0_wp)*(refangles1_stheta(3)*refangles1_stheta(4)-1.0e0_wp)
            refangles1_theta0(1,3) = tau+(tau-halfpi)*(refangles1_sphi(1)*refangles1_sphi(3)-1.0e0_wp) +(tau-twopi/ &
                & 3.0e0_wp)*(refangles1_stheta(2)*refangles1_stheta(4)-1.0e0_wp)
            refangles1_theta0(1,4) = tauco+(tau-halfpi)*(refangles1_sphi(1)*refangles1_sphi(4)-1.0e0_wp) +(tau-twopi/ &
                & 3.0e0_wp)*(refangles1_stheta(2)*refangles1_stheta(3)-1.0e0_wp)
            refangles1_theta0(2,3) = tau+(tau-halfpi)*(refangles1_sphi(2)*refangles1_sphi(3)-1.0e0_wp) +(tau-twopi/ &
                & 3.0e0_wp)*(refangles1_stheta(1)*refangles1_stheta(4)-1.0e0_wp)
            refangles1_theta0(2,4) = tauco+(tau-halfpi)*(refangles1_sphi(2)*refangles1_sphi(4)-1.0e0_wp) +(tau-twopi/ &
                & 3.0e0_wp)*(refangles1_stheta(1)*refangles1_stheta(3)-1.0e0_wp)
            refangles1_theta0(3,4) = tauco+(tau-halfpi)*(refangles1_sphi(3)*refangles1_sphi(4)-1.0e0_wp) +(tau-twopi/ &
                & 3.0e0_wp)*(refangles1_stheta(1)*refangles1_stheta(2)-1.0e0_wp)
            !  el atomo 6 es el c1
            refangles1_theta0(5,6) = tauoc+(tauoc-halfpi)*(refangles1_sphi(5)*refangles1_sphi(6)-1.0e0_wp)
            !  fill in the other half of the matrix
            do i=1, 3
                do j=i+1, 4
                    refangles1_theta0(j,i) = refangles1_theta0(i,j)
                enddo
            enddo
            refangles1_theta0(6,5) = refangles1_theta0(5,6)
        end subroutine refangles

        !> @brief Evaluate stretch for the CH4OCl_VBMM_2023 potential.
        subroutine stretch( vstr )
            !  subroutine to calculate leps-type stretching potential and its
            !  derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: ach1, ach2, arga1, arga2, dt, e1, e3, expterm, rav1, rav2, vcc, vj(4), vjb1h(3), vjb2h(1), vjc1b, &
                & vjc1h(3), vjc2b, vjc2h(1), vq(4), vqb1h(3), vqb2h(1), vqc1b, vqc1h(3), vqc2b, vqc2h(1), vstr
            integer :: i
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            !  calculate avergage bond length for the methane moiety
            rav1 = (coord5_rc1h(1)+coord5_rc1h(2)+coord5_rc1h(3))/3.0e0_wp
            rav2 = coord5_rc2h(1)
            !  initialise:
            vstr = 0.0e0_wp
            !  ach:
            !  nb: in double precision tanh(19.0d0)=1.0d0 and we put the if statemen
            !  in to avoid overflow/underflow errors
            arga1 = potcm_c1c1h*(rav1-coord6_r0c1h)
            arga2 = potcm_c1c2h*(rav2-coord6_r0c2h)
            if( arga1 .lt. 19.0e0_wp )then
                ach1 = potcm_a1c1h+potcm_b1c1h*(tanh(arga1)+1.0e0_wp)*0.5e0_wp
            else
                ach1 = potcm_a1c1h+potcm_b1c1h
            endif
            if( arga2 .lt. 19.0e0_wp )then
                ach2 = potcm_a1c2h+potcm_b1c2h*(tanh(arga2)+1.0e0_wp)*0.5e0_wp
            else
                ach2 = potcm_a1c2h+potcm_b1c2h
            endif
            !  calculate singlet: e1, triplet: e3 energies and vq and vj
            !  terms for each bond
            e1    = potcm_d1c1b*(exp(-2.0e0_wp*potcm_ac1b*(coord5_rc1b-potcm_r0c1b)) -2.0e0_wp*exp(-potcm_ac1b* &
                & (coord5_rc1b- potcm_r0c1b)))
            e3    = potcm_d3c1b*(exp(-2.0e0_wp*potcm_ac1b*(coord5_rc1b-potcm_r0c1b)) +2.0e0_wp*exp(-potcm_ac1b* &
                & (coord5_rc1b- potcm_r0c1b)))
            vqc1b = (e1+e3)*0.5e0_wp
            vjc1b = (e1-e3)*0.5e0_wp
            e1    = potcm_d1c2b*(exp(-2.0e0_wp*potcm_ac2b*(coord5_rc2b-potcm_r0c2b)) -2.0e0_wp*exp(-potcm_ac2b* &
                & (coord5_rc2b- potcm_r0c2b)))
            e3    = potcm_d3c2b*(exp(-2.0e0_wp*potcm_ac2b*(coord5_rc2b-potcm_r0c2b)) +2.0e0_wp*exp(-potcm_ac2b* &
                & (coord5_rc2b- potcm_r0c2b)))
            vqc2b = (e1+e3)*0.5e0_wp
            vjc2b = (e1-e3)*0.5e0_wp
            do i=1, 3
                e1       = potcm_d1c1h*(exp(-2.0e0_wp*ach1*(coord5_rc1h(i)-coord6_r0c1h)) -2.0e0_wp*exp(-ach1* &
                    & (coord5_rc1h(i)-coord6_r0c1h)))
                e3       = potcm_d3c1h*(exp(-2.0e0_wp*ach1*(coord5_rc1h(i)-coord6_r0c1h)) +2.0e0_wp*exp(-ach1* &
                    & (coord5_rc1h(i)-coord6_r0c1h)))
                vqc1h(i) = (e1+e3)*0.5e0_wp
                vjc1h(i) = (e1-e3)*0.5e0_wp
            enddo
            e1       = potcm_d1c2h*(exp(-2.0e0_wp*ach2*(coord5_rc2h(1)-coord6_r0c2h)) -2.0e0_wp*exp(-ach2* &
                & (coord5_rc2h(1)- coord6_r0c2h)))
            e3       = potcm_d3c2h*(exp(-2.0e0_wp*ach2*(coord5_rc2h(1)-coord6_r0c2h)) +2.0e0_wp*exp(-ach2* &
                & (coord5_rc2h(1)- coord6_r0c2h)))
            vqc2h(1) = (e1+e3)*0.5e0_wp
            vjc2h(1) = (e1-e3)*0.5e0_wp
            do i=1, 3
                e1       = potcm_d1h1h*(exp(-2.0e0_wp*potcm_ah1h*(coord5_rb1h(i)-potcm_r0h1h)) -2.0e0_wp*exp(- &
                    & potcm_ah1h* (coord5_rb1h(i)-potcm_r0h1h)))
                e3       = potcm_d3h1h*(exp(-2.0e0_wp*potcm_ah1h*(coord5_rb1h(i)-potcm_r0h1h)) +2.0e0_wp*exp(- &
                    & potcm_ah1h* (coord5_rb1h(i)-potcm_r0h1h)))
                vqb1h(i) = (e1+e3)*0.5e0_wp
                vjb1h(i) = (e1-e3)*0.5e0_wp
            enddo
            e1       = potcm_d1h2h*(exp(-2.0e0_wp*potcm_ah2h*(coord5_rb2h(1)-potcm_r0h2h)) -2.0e0_wp*exp(-potcm_ah2h* &
                & (coord5_rb2h(1)-potcm_r0h2h)))
            e3       = potcm_d3h2h*(exp(-2.0e0_wp*potcm_ah2h*(coord5_rb2h(1)-potcm_r0h2h)) +2.0e0_wp*exp(-potcm_ah2h* &
                & (coord5_rb2h(1)-potcm_r0h2h)))
            vqb2h(1) = (e1+e3)*0.5e0_wp
            vjb2h(1) = (e1-e3)*0.5e0_wp
            !  calculate 3 body potential
            do i=1, 3
                vq(i) = vqc1h(i)+vqc1b+vqb1h(i)
                vj(i) = -sqrt(((vjc1h(i)-vjc1b)**2+(vjc1b-vjb1h(i))**2 +(vjb1h(i)-vjc1h(i))**2)*0.5e0_wp)
                vstr  = vstr+vq(i)+vj(i)
            enddo
            vq(4) = vqc2h(1)+vqc2b+vqb2h(1)
            vj(4) = -sqrt(((vjc2h(1)-vjc2b)**2+(vjc2b-vjb2h(1))**2 +(vjb2h(1)-vjc2h(1))**2)*0.5e0_wp)
            vstr  = vstr+vq(4)+vj(4)
            ! CRD 2012 potencial morse para el enlace C-C
            dt      = (coord5_rcc-coord6_r0cc)
            expterm = exp(-potcm_acc*dt)
            vcc     = potcm_d1cc*(1.0e0_wp-expterm)**2.0e0_wp
            vstr    = vstr + vcc
        end subroutine stretch

        !> @brief Evaluate opbend for the CH4OCl_VBMM_2023 potential.
        subroutine opbend( vop )
            !  subroutine calculates symmetrized vop potential and derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: sum2, sum4, sumd2(4), sumd4(4), vop
            integer  :: i, j, k, l
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            vop = 0.0e0_wp
            !  calculate force constants
            call opforce
            !  calculate out-of-plane angle
            ! CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C2
            do i=1, 4
                j = i+1
                if( j .gt. 4 ) j=j-4
                k = j+1
                if( k .gt. 4 ) k=k-4
                l = k+1
                if( l .gt. 4 ) l=l-4
                call calcdelta1(i, j, k, l, sum2, sum4)
                sumd2(i) = sum2
                sumd4(i) = sum4
                vop      = vop+vop6_fdelta(i)*sumd2(i)+vop6_hdelta(i)*sumd4(i)
            enddo
        end subroutine opbend

        !> @brief Evaluate calcdelta1 for the CH4OCl_VBMM_2023 potential.
        subroutine calcdelta1( i, j, k, l, sum2, sum4 )
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3), argd(8), axb(3), b(3), c(8,3), delta(8), sum2, sum4
            integer  :: i, ii, in(3), ind, j, k, l
            real(wp) :: norma
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            !  initialise
            sum2 = 0.0e0_wp; sum4 = 0.0e0_wp
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
                    a(ind) = qpdot_q(coord2_n1h(k,ind))-qpdot_q(coord2_nc2(ind))
                    b(ind) = qpdot_q(coord2_n1h(l,ind))-qpdot_q(coord2_nc2(ind))
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
                    a(ind) = qpdot_q(coord2_nc2(ind))-qpdot_q(coord2_n1h(j,ind))
                    b(ind) = qpdot_q(coord2_n1h(l,ind))-qpdot_q(coord2_n1h(j,ind))
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
                    a(ind) = qpdot_q(coord2_n1h(k,ind))-qpdot_q(coord2_n1h(j,ind))
                    b(ind) = qpdot_q(coord2_nc2(ind))-qpdot_q(coord2_n1h(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                ! modificacion para que hiciese vxu
                norma = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma = sqrt(norma)
            else
            endif
            ! CRD 2012 eliminamos porque seria la salida del C2
            if( k .lt. 4.and.j .lt. 4.and.l .lt. 4 )then
                do ind=1, 3
                    a(ind) = qpdot_q(coord2_n1h(k,ind))-qpdot_q(coord2_n1h(j,ind))
                    b(ind) = qpdot_q(coord2_n1h(l,ind))-qpdot_q(coord2_n1h(j,ind))
                enddo
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
            else
            endif
            !  c is position vector of h(ii): calculate c(j),c(k),c(l)
            do ii=1, 3
                if( in(ii) .eq. 4 )then
                    do ind=1, 3
                        c(in(ii),ind) = -coord3_tcc(ind)/coord5_rcc
                    enddo
                else
                    do ind=1, 3
                        c(in(ii),ind) = -coord4_tc1h(in(ii),ind)/coord5_rc1h(in(ii))
                    enddo
                endif
            enddo
            !  argd is the dot product axb dot c
            do ii=1, 3
                argd(in(ii))  = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                argd(in(ii))  = argd(in(ii))/norma
                delta(in(ii)) = acos(argd(in(ii)))-refangles1_theta0(i,in(ii))
                sum2          = sum2+delta(in(ii))**2
                sum4          = sum4+delta(in(ii))**4
            enddo
        end subroutine calcdelta1

        !> @brief Evaluate opforce for the CH4OCl_VBMM_2023 potential.
        subroutine opforce
            !  calculates the out-of-plane bending force constants
            !  and their derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: switch(4)
            integer  :: i
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            !  calculate switching functions:
            switch(1) = (1.0e0_wp-switch1_s3(1))*switch1_s3(2)*switch1_s3(3)*switch1_s3(4)
            switch(2) = (1.0e0_wp-switch1_s3(2))*switch1_s3(3)*switch1_s3(4)*switch1_s3(1)
            switch(3) = (1.0e0_wp-switch1_s3(3))*switch1_s3(4)*switch1_s3(1)*switch1_s3(2)
            switch(4) = (1.0e0_wp-switch1_s3(4))*switch1_s3(1)*switch1_s3(2)*switch1_s3(3)
            !  calculate the force constants
            do i=1, 4
                vop6_fdelta(i) = switch(i)*potcm_fch3
                vop6_hdelta(i) = switch(i)*potcm_hch3
            enddo
        end subroutine opforce

        !> @brief Evaluate ipbend for the CH4OCl_VBMM_2023 potential.
        subroutine ipbend( vip )
            !  subroutine calculates symmetrised in plane bend term
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: costh(6,6), dth(6,6), theta(6,6), vip
            integer  :: i, j
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            !  initialise
            vip = 0.0e0_wp
            !  calculate force constants: fk0(i,j), f1(i)
            !  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
            call ipforce
            !  calculate theta(i,j) and in plane bend potential
            ! c1
            do i=1, 3
                do j=i+1, 3
                    costh(i,j) = coord4_tc1h(i,1)*coord4_tc1h(j,1)+coord4_tc1h(i,2)*coord4_tc1h(j,2) + coord4_tc1h(i, &
                        & 3)*coord4_tc1h(j,3)
                    costh(i,j) = costh(i,j)/coord5_rc1h(i)/coord5_rc1h(j)
                    theta(i,j) = acos(costh(i,j))
                    dth(i,j)   = theta(i,j)-refangles1_theta0(i,j)
                    vip        = vip+0.5e0_wp*ipforce1_fk0(i,j)*ipforce1_f1(i)*ipforce1_f1(j)*dth(i,j)**2
                enddo
                j = 4
                costh(i,j) = coord4_tc1h(i,1)*coord3_tcc(1)+coord4_tc1h(i,2)*coord3_tcc(2) +coord4_tc1h(i,3)* &
                    & coord3_tcc(3)
                costh(i,j) = costh(i,j)/coord5_rc1h(i)/coord5_rcc
                theta(i,j) = acos(costh(i,j))
                dth(i,j)   = theta(i,j)-refangles1_theta0(i,j)
                vip        = vip+0.5e0_wp*ipforce1_fk0(i,j)*ipforce1_f1(i)*ipforce1_f1(j)*dth(i,j)**2
            enddo
            ! c2
            costh(5,6) = coord4_tc2h(1,1)*(-coord3_tcc(1))+coord4_tc2h(1,2)*(-coord3_tcc(2)) +coord4_tc2h(1,3)*(- &
                & coord3_tcc(3))
            costh(5,6) = costh(5,6)/coord5_rc2h(1)/coord5_rcc
            theta(5,6) = acos(costh(5,6))
            dth(5,6)   = theta(5,6)-refangles1_theta0(5,6)
            vip        = vip+0.5e0_wp*ipforce1_fk0(5,6)*ipforce1_f1(5)*ipforce1_f1(6)*dth(5,6)**2
        end subroutine ipbend

        !> @brief Evaluate ipforce for the CH4OCl_VBMM_2023 potential.
        subroutine ipforce
            !  calculates the symmetrised in plane bend force constants and
            !  all partial derivatives involving them
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: a1, a2, arga1, arga2, f0, f2
            real(wp), parameter :: r0hh = 0.0e0_wp
            integer             :: i
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            !  set force constant at asymptotes
            f0 = potcm_fkinf1+potcm_ak1
            f2 = potcm_fkinf1
            ipforce1_fk0(1,2) = f0+f0*(switch1_s1(1)*switch1_s1(2)-1.0e0_wp)+(f0-f2)*(switch1_s2(3)*switch1_s2(4)- &
                & 1.0e0_wp)
            ipforce1_fk0(1,3) = f0+f0*(switch1_s1(1)*switch1_s1(3)-1.0e0_wp)+(f0-f2)*(switch1_s2(2)*switch1_s2(4)- &
                & 1.0e0_wp)
            ipforce1_fk0(1,4) = f0+f0*(switch1_s1(1)*switch1_s1(4)-1.0e0_wp)+(f0-f2)*(switch1_s2(2)*switch1_s2(3)- &
                & 1.0e0_wp)
            ipforce1_fk0(2,3) = f0+f0*(switch1_s1(2)*switch1_s1(3)-1.0e0_wp)+(f0-f2)*(switch1_s2(1)*switch1_s2(4)- &
                & 1.0e0_wp)
            ipforce1_fk0(2,4) = f0+f0*(switch1_s1(2)*switch1_s1(4)-1.0e0_wp)+(f0-f2)*(switch1_s2(1)*switch1_s2(3)- &
                & 1.0e0_wp)
            ipforce1_fk0(3,4) = f0+f0*(switch1_s1(3)*switch1_s1(4)-1.0e0_wp)+(f0-f2)*(switch1_s2(1)*switch1_s2(2)- &
                & 1.0e0_wp)
            ipforce1_fk0(5,6) = f0+f0*(switch1_s1(5)*switch1_s1(6)-1.0e0_wp)+(f0-f2)*(switch1_s2(5)*switch1_s2(6)- &
                & 1.0e0_wp)
            ! f1 para los h del c1
            do i=1, 3
                !  calc derivatives of fk0 wrt each of the rch(i) bonds
                !  calculate the terms f1(i)
                arga1          = potcm_aa11*coord5_rb1h(i)*coord5_rb1h(i)
                arga2          = potcm_aa41*(coord5_rb1h(i)-r0hh)*(coord5_rb1h(i)-r0hh)
                a1             = 1.0e0_wp-exp(-arga1)
                a2             = potcm_aa21+potcm_aa31*exp(-arga2)
                ipforce1_f1(i) = a1*exp(-a2*(coord5_rc1h(i)-coord6_r0c1h)**2)
            enddo
            ! f1 para el c2
            ! CRD 2012 eliminamos el calculo de a1c y a2c pq es constante
            ! lo incluimos en el CONST
            ipforce1_f1(4) = potcm_a1c*exp(-potcm_a2c*(coord5_rcc-coord6_r0cc)**2)
            ipforce1_f1(6) = ipforce1_f1(4)
            ! f1 para los h del c2
            !  calc derivatives of fk0 wrt each of the rch(i) bonds
            !  calculate the terms f1(i)
            arga1          = potcm_aa12*coord5_rb2h(1)*coord5_rb2h(1)
            arga2          = potcm_aa42*(coord5_rb2h(1)-r0hh)*(coord5_rb2h(1)-r0hh)
            a1             = 1.0e0_wp-exp(-arga1)
            a2             = potcm_aa22+potcm_aa32*exp(-arga2)
            ipforce1_f1(5) = a1*exp(-a2*(coord5_rc2h(1)-coord6_r0c2h)**2)
        end subroutine ipforce

        !> @brief Evaluate switchf for the CH4OCl_VBMM_2023 potential.
        subroutine switchf
            !  calculates switching functions: s3,sphi,stheta
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: argmax, args1, args2, args3, argsc1, argsc2, argsc3, argsphi, argsphic, argstheta, argsthetac
            integer  :: i
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
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
            argmax = 19.0e0_wp
            !  calculate s1
            do i=1, 3
                args1 = potcm_a1s1*(coord5_rc1h(i)-coord6_r0c1h)*(coord5_rc1h(i)-potcm_b1s1)**8
                if( args1 .lt. argmax )then
                    switch1_s1(i) = 1.0e0_wp-tanh(args1)
                else
                    switch1_s1(i) = 0.0e0_wp
                endif
            enddo
            args1 = potcm_a2s1*(coord5_rc2h(1)-coord6_r0c2h)*(coord5_rc2h(1)-potcm_b2s1)**8
            if( args1 .lt. argmax )then
                switch1_s1(5) = 1.0e0_wp-tanh(args1)
            else
                switch1_s1(5) = 0.0e0_wp
            endif
            ! CRD 2012 todos los S de los carbonos son 1
            argsc1 = potcm_a1scc*(coord5_rcc-coord6_r0cc)*(coord5_rcc-potcm_b1scc)**8
            if( argsc1 .lt. argmax )then
                switch1_s1(4) = 1.0e0_wp-tanh(argsc1)
                switch1_s1(6) = switch1_s1(4)
            else
                switch1_s1(4) = 0.0e0_wp
                switch1_s1(6) = switch1_s1(4)
            endif
            ! CRD correccion para dejar C-C fijo
            switch1_s1(4) = 1.0e0_wp
            switch1_s1(6) = 1.0e0_wp
            !  calculate s2 and ds2
            do i=1, 3
                args2 = potcm_a1s2*(coord5_rc1h(i)-coord6_r0c1h)*(coord5_rc1h(i)-potcm_b1s2)**6
                if( args2 .lt. argmax )then
                    switch1_s2(i) = 1.0e0_wp-tanh(args2)
                else
                    switch1_s2(i) = 0.0e0_wp
                endif
            enddo
            args2 = potcm_a2s2*(coord5_rc2h(1)-coord6_r0c2h)*(coord5_rc2h(1)-potcm_b2s2)**6
            if( args2 .lt. argmax )then
                switch1_s2(5) = 1.0e0_wp-tanh(args2)
            else
                switch1_s2(5) = 0.0e0_wp
            endif
            argsc2 = potcm_a2scc*(coord5_rcc-coord6_r0cc)*(coord5_rcc-potcm_b2scc)**6
            if( argsc2 .lt. argmax )then
                switch1_s2(4) = 1.0e0_wp-tanh(argsc2)
                switch1_s2(6) = switch1_s2(4)
            else
                switch1_s2(4) = 0.0e0_wp
                switch1_s2(6) = switch1_s2(4)
            endif
            ! CRD correccion para dejar C-C fijo
            switch1_s2(4) = 1.0e0_wp
            switch1_s2(6) = 1.0e0_wp
            !  calculate s3 and ds3
            do i=1, 3
                args3 = potcm_a3s*(coord5_rc1h(i)-coord6_r0c1h)*(coord5_rc1h(i)-potcm_b3s)**2
                if( args3 .lt. argmax )then
                    switch1_s3(i) = 1.0e0_wp-tanh(args3)
                else
                    switch1_s3(i) = 0.0e0_wp
                endif
            enddo
            argsc3 = potcm_a3scc*(coord5_rcc-coord6_r0cc)*(coord5_rcc-potcm_b3scc)**2
            if( argsc3 .lt. argmax )then
                switch1_s3(4) = 1.0e0_wp-tanh(argsc3)
            else
                switch1_s3(4) = 0.0e0_wp
            endif
            ! CRD correccion para dejar C-C fijo
            switch1_s3(4) = 1.0e0_wp
            !  calculate sphi
            !  condition here is on the bondlength rch(i)
            !  st argsphi is lt approx 19.0d0
            do i=1, 3
                if( coord5_rc1h(i) .lt. 3.8e0_wp )then
                    argsphi            = potcm_aphi1*(coord5_rc1h(i)-coord6_r0c1h)*exp(potcm_bphi1*(coord5_rc1h(i)- &
                        & potcm_cphi1)**3)
                    refangles1_sphi(i) = 1.0e0_wp-tanh(argsphi)
                else
                    refangles1_sphi(i) = 0.0e0_wp
                endif
            enddo
            if( coord5_rc2h(1) .lt. 3.8e0_wp )then
                argsphi            = potcm_aphi2*(coord5_rc2h(1)-coord6_r0c2h)*exp(potcm_bphi2*(coord5_rc2h(1)- &
                    & potcm_cphi2)**3)
                refangles1_sphi(5) = 1.0e0_wp-tanh(argsphi)
            else
                refangles1_sphi(5) = 0.0e0_wp
            endif
            if( coord5_rcc .lt. 3.8e0_wp )then
                argsphic           = potcm_aphi*(coord5_rcc-coord6_r0cc)*exp(potcm_bphi*(coord5_rcc-potcm_cphi)**3)
                refangles1_sphi(4) = 1.0e0_wp-tanh(argsphic)
                refangles1_sphi(6) = refangles1_sphi(4)
            else
                refangles1_sphi(4) = 0.0e0_wp
                refangles1_sphi(6) = refangles1_sphi(4)
            endif
            !  calculate stheta and dstheta
            do i=1, 3
                if( coord5_rc1h(i) .lt. 3.8e0_wp )then
                    argstheta = potcm_atheta*(coord5_rc1h(i)-coord6_r0c1h)*exp(potcm_btheta*(coord5_rc1h(i)- &
                        & potcm_ctheta)**3)
                    refangles1_stheta(i) = 1.0e0_wp-tanh(argstheta)
                else
                    refangles1_stheta(i) = 0.0e0_wp
                endif
            enddo
            ! CRD 2012 todos los S de los carbonos son 1
            if( coord5_rcc .lt. 3.8e0_wp )then
                argsthetac           = potcm_athetacc*(coord5_rcc-coord6_r0cc)*exp(potcm_bthetacc*(coord5_rcc- &
                    & potcm_cthetacc)**3)
                refangles1_stheta(4) = 1.0e0_wp-tanh(argsthetac)
                refangles1_stheta(6) = refangles1_stheta(4)
            else
                refangles1_stheta(4) = 0.0e0_wp
                refangles1_stheta(6) = refangles1_stheta(4)
            endif
        end subroutine switchf

        !> @brief Evaluate prepare_local_context for the CH4OCl_VBMM_2023 potential.
        subroutine prepare_local_context
            !   N3TMMN = 3 * NATOMS
            !   NATOMS = the number of atoms represented by this potential function
            !   The variable N3TMMN is the minimum value of N3TM allowed to be
            !   passed by the calling routine for the number of cartesian
            !   coordinates needed to represent the full system represented by this
            !   potential energy surface routine.
            !   N3TM must be greater than or equal to N3TMMN.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: fact1, fact2
            integer  :: i, icount, ind
            !  include file for CH3OH+H PES
            !  CONST
            ! Subroutine Coorden
            ! Subroutine torsion
            ! Subroutine refangles
            ! Subroutine vop
            ! Subroutine ipforce
            ! Subroutine Switch
            !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !      REF(1)='J. Espinosa-Garcia and J. C. Corchado'
            !  calculate indexes for coordinates
            do ind=1, 3
                icount          = ind-3
                coord2_nc1(ind) = 3*potcm_nnc1+icount
                coord2_nc2(ind) = 3*potcm_nnc2+icount
                coord2_nhb(ind) = 3*potcm_nnb+icount
                do i=1, 3
                    coord2_n1h(i,ind) = 3*potcm2_nc1h(i)+icount
                enddo
                coord2_n2h(1,ind) = 3*potcm2_nc2h(1)+icount
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
            fact1        = 0.041840e0_wp
            fact2        = 6.022045e0_wp
            potcm_d1c1h  = potcm_d1c1h*fact1
            potcm_d3c1h  = potcm_d3c1h*fact1
            potcm_d1c2h  = potcm_d1c2h*fact1
            potcm_d3c2h  = potcm_d3c2h*fact1
            potcm_d1cc   = potcm_d1cc*fact1
            potcm_d1c1b  = potcm_d1c1b*fact1
            potcm_d3c1b  = potcm_d3c1b*fact1
            potcm_d1c2b  = potcm_d1c2b*fact1
            potcm_d3c2b  = potcm_d3c2b*fact1
            potcm_d1h1h  = potcm_d1h1h*fact1
            potcm_d3h1h  = potcm_d3h1h*fact1
            potcm_d1h2h  = potcm_d1h2h*fact1
            potcm_d3h2h  = potcm_d3h2h*fact1
            potcm_fkinf1 = potcm_fkinf1*fact2
            potcm_ak1    = potcm_ak1*fact2
            potcm_fkinf2 = potcm_fkinf2*fact2
            potcm_a1s1   = potcm_a1s1*1.e-8_wp
            potcm_a2s1   = potcm_a2s1*1.e-8_wp
            potcm_a1s2   = potcm_a1s2*1.e-8_wp
            potcm_a2s2   = potcm_a2s2*1.e-8_wp
            potcm_a1scc  = potcm_a1scc*1.e-8_wp
            potcm_a2scc  = potcm_a2scc*1.e-8_wp
        end subroutine prepare_local_context

        !> @brief Evaluate POTINFO for the CH4OCl_VBMM_2023 potential.
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

        !> @brief Evaluate CARTTOR for the CH4OCl_VBMM_2023 potential.
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

        !> @brief Evaluate EUNITZERO for the CH4OCl_VBMM_2023 potential.
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

        !> @brief Evaluate RTOCART for the CH4OCl_VBMM_2023 potential.
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

        !> @brief Evaluate DEDCOU for the CH4OCl_VBMM_2023 potential.
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
end module chempotpy_ch4ocl_vbmm_2023_pes

!> @brief Expose CH4OCl_VBMM_2023 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_ch4ocl_vbmm_2023_pes, only: evaluate
    use chempotpy_kinds,                only: wp
    implicit none
    integer, parameter    :: natoms  = 7
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
