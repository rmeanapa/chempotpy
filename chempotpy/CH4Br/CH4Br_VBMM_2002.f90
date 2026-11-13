! System:          CH4Br
! Functional form:
! Common name:
! Reference  :     J. Espinosa-Garcia, J. Chem. Phys., Vol. 117, 2076 (2
! Number of bodies: 6
! Number of electronic states: 1
! Number of derivatives: 0
! Interface: Section-2
! Data file: immutable CONST.dat parameters embedded below
! Notes:  CH4+Br based on the surface by
!         Jordan and Gilbert, JCP, 102, 5669 (1995).

!> @brief Implement CH4Br_VBMM_2002, a potential energy surface for CH4Br.
module chempotpy_ch4br_vbmm_2002_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate
    ! Released CONST.dat data, SHA-256:
    ! a76be154bf20b50c18c1d0ccd572aa29f19769e16be54662c30cb2b3a36adc43
    integer, parameter  :: source_carbon_index        = 2
    integer, parameter  :: source_bromine_index       = 6
    integer, parameter  :: source_hydrogen_indices(4) = [3, 4, 5, 1]
    real(wp), parameter :: source_ch_equilibrium      = 1.09397_wp
    real(wp), parameter :: source_ch_dissociation_1   = 112.341_wp
    real(wp), parameter :: source_ch_dissociation_3   = 39.101_wp
    real(wp), parameter :: source_ch_attenuation_a    = 1.71000_wp
    real(wp), parameter :: source_ch_attenuation_b    = 0.12000_wp
    real(wp), parameter :: source_ch_attenuation_c    = 6.61404_wp
    real(wp), parameter :: source_brh_equilibrium     = 1.41443_wp
    real(wp), parameter :: source_brh_dissociation_1  = 90.226_wp
    real(wp), parameter :: source_brh_dissociation_3  = 43.309_wp
    real(wp), parameter :: source_brh_attenuation     = 1.8500_wp
    real(wp), parameter :: source_cb_equilibrium      = 1.89397_wp
    real(wp), parameter :: source_cb_dissociation_1   = 87.149_wp
    real(wp), parameter :: source_cb_dissociation_3   = 23.553_wp
    real(wp), parameter :: source_cb_attenuation      = 1.888555_wp
    real(wp), parameter :: source_out_switch_a        = 0.1419147_wp
    real(wp), parameter :: source_out_switch_b        = -0.306845_wp
    real(wp), parameter :: source_phi_a               = 0.5287903_wp
    real(wp), parameter :: source_phi_b               = 0.4006638_wp
    real(wp), parameter :: source_phi_c               = 1.9209937_wp
    real(wp), parameter :: source_theta_a             = 0.9078714_wp
    real(wp), parameter :: source_theta_b             = 0.3548859_wp
    real(wp), parameter :: source_theta_c             = 1.8915497_wp
    real(wp), parameter :: source_fch3                = 0.0957500_wp
    real(wp), parameter :: source_hch3                = 0.1915000_wp
    real(wp), parameter :: source_force_infinite      = 0.3215000_wp
    real(wp), parameter :: source_force_attenuation   = 0.1560000_wp
    real(wp), parameter :: source_force_b             = 10.7132_wp
    real(wp), parameter :: source_force_shape_1       = 0.313952_wp
    real(wp), parameter :: source_force_shape_2       = 2.099963_wp
    real(wp), parameter :: source_force_shape_3       = 3.265953_wp
    real(wp), parameter :: source_force_shape_4       = 11.569977_wp

    contains

    !> @brief Evaluate the CH4Br_VBMM_2002 potential energy surface for CH4Br.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units,    only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_geometry, only: flatten_cartesian_coordinates
        use chempotpy_units,    only: hartree_to_electronvolt
        use chempotpy_kinds,    only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 6
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: angles_dtheta0(4,4,4)
        real(wp)              :: angles_theta0(4,4)
        real(wp)              :: bonds_r0cb
        real(wp)              :: bonds_r0ch
        real(wp)              :: bonds_r0hh
        real(wp)              :: bonds_rbh(4)
        real(wp)              :: bonds_rcb
        real(wp)              :: bonds_rch(4)
        real(wp)              :: coords_tbh(4,3)
        real(wp)              :: coords_tcb(3)
        real(wp)              :: coords_tch(4,3)
        real(wp)              :: delta1_fdelta(4)
        real(wp)              :: delta1_hdelta(4)
        real(wp)              :: delta2_dfdelta(4,4)
        real(wp)              :: delta2_dhdelta(4,4)
        real(wp)              :: force1_dfdc(4,4,4)
        real(wp)              :: force1_dfdh(4,4,4)
        real(wp)              :: force1_f1(4)
        real(wp)              :: force1_fk0(4,4)
        real(wp)              :: force2_aa1
        real(wp)              :: force2_aa2
        real(wp)              :: force2_aa3
        real(wp)              :: force2_aa4
        real(wp)              :: force2_ak
        real(wp)              :: force2_bk
        real(wp)              :: force2_fkinf
        real(wp)              :: fsw1_a1s
        real(wp)              :: fsw1_a2s
        real(wp)              :: fsw1_a3s
        real(wp)              :: fsw1_b1s
        real(wp)              :: fsw1_b2s
        real(wp)              :: fsw1_b3s
        real(wp)              :: ip1_ds1(4)
        real(wp)              :: ip1_ds2(4)
        real(wp)              :: ip1_s1(4)
        real(wp)              :: ip1_s2(4)
        integer               :: ndx_nc(3)
        integer               :: ndx_nh(4,3)
        integer               :: ndx_nhb(3)
        integer               :: ndx_nnb
        integer               :: ndx_nnc
        integer               :: ndx_nnh(4)
        real(wp)              :: op1_ds3(4)
        real(wp)              :: op1_s3(4)
        real(wp)              :: params_fch3
        real(wp)              :: params_hch3
        real(wp)              :: qpdot_pdot(150)
        real(wp)              :: qpdot_q(150)
        real(wp)              :: sphi1_aphi
        real(wp)              :: sphi1_bphi
        real(wp)              :: sphi1_cphi
        real(wp)              :: stheta1_atheta
        real(wp)              :: stheta1_btheta
        real(wp)              :: stheta1_ctheta
        real(wp)              :: stret1_a1ch
        real(wp)              :: stret1_acb
        real(wp)              :: stret1_ahh
        real(wp)              :: stret1_b1ch
        real(wp)              :: stret1_c1ch
        real(wp)              :: stret1_d1cb
        real(wp)              :: stret1_d1ch
        real(wp)              :: stret1_d1hh
        real(wp)              :: stret1_d3cb
        real(wp)              :: stret1_d3ch
        real(wp)              :: stret1_d3hh
        real(wp)              :: switch1_dsphi(4)
        real(wp)              :: switch1_dstheta(4)
        real(wp)              :: switch1_sphi(4)
        real(wp)              :: switch1_stheta(4)
        integer               :: istate
        real(wp)              :: tx(18), dx(18), v
        integer               :: iatom, idir, i, j
        !initialize
        angles_dtheta0 = 0.e0_wp; angles_theta0 = 0.e0_wp; bonds_r0cb = 0.e0_wp; bonds_r0ch = 0.e0_wp; bonds_r0hh = 0.e0_wp
        bonds_rbh      = 0.e0_wp; bonds_rcb = 0.e0_wp; bonds_rch = 0.e0_wp; coords_tbh = 0.e0_wp; coords_tcb = 0.e0_wp; coords_tch = 0.e0_wp
        delta1_fdelta  = 0.e0_wp; delta1_hdelta = 0.e0_wp; delta2_dfdelta = 0.e0_wp; delta2_dhdelta = 0.e0_wp; force1_dfdc = 0.e0_wp
        force1_dfdh    = 0.e0_wp; force1_f1 = 0.e0_wp; force1_fk0 = 0.e0_wp; force2_aa1 = 0.e0_wp; force2_aa2 = 0.e0_wp; force2_aa3 = 0.e0_wp
        force2_aa4     = 0.e0_wp; force2_ak = 0.e0_wp; force2_bk = 0.e0_wp; force2_fkinf = 0.e0_wp; fsw1_a1s = 0.e0_wp; fsw1_a2s = 0.e0_wp
        fsw1_a3s       = 0.e0_wp; fsw1_b1s = 0.e0_wp; fsw1_b2s = 0.e0_wp; fsw1_b3s = 0.e0_wp; ip1_ds1 = 0.e0_wp; ip1_ds2 = 0.e0_wp
        ip1_s1         = 0.e0_wp; ip1_s2 = 0.e0_wp; ndx_nc = 0; ndx_nh = 0; ndx_nhb = 0; ndx_nnb = 0; ndx_nnc = 0; ndx_nnh = 0
        op1_ds3        = 0.e0_wp; op1_s3 = 0.e0_wp; params_fch3 = 0.e0_wp; params_hch3 = 0.e0_wp; qpdot_pdot = 0.e0_wp; qpdot_q = 0.e0_wp
        sphi1_aphi     = 0.e0_wp; sphi1_bphi = 0.e0_wp; sphi1_cphi = 0.e0_wp; stheta1_atheta = 0.e0_wp; stheta1_btheta = 0.e0_wp
        stheta1_ctheta = 0.e0_wp; stret1_a1ch = 0.e0_wp; stret1_acb = 0.e0_wp; stret1_ahh = 0.e0_wp; stret1_b1ch = 0.e0_wp
        stret1_c1ch    = 0.e0_wp; stret1_d1cb = 0.e0_wp; stret1_d1ch = 0.e0_wp; stret1_d1hh = 0.e0_wp; stret1_d3cb = 0.e0_wp
        stret1_d3ch    = 0.e0_wp; stret1_d3hh = 0.e0_wp; switch1_dsphi = 0.e0_wp; switch1_dstheta = 0.e0_wp; switch1_sphi = 0.e0_wp
        switch1_stheta = 0.e0_wp
        ndx_nnc        = source_carbon_index
        ndx_nnb        = source_bromine_index
        ndx_nnh        = source_hydrogen_indices
        bonds_r0ch     = source_ch_equilibrium
        stret1_d1ch    = source_ch_dissociation_1
        stret1_d3ch    = source_ch_dissociation_3
        stret1_a1ch    = source_ch_attenuation_a
        stret1_b1ch    = source_ch_attenuation_b
        stret1_c1ch    = source_ch_attenuation_c
        bonds_r0hh     = source_brh_equilibrium
        stret1_d1hh    = source_brh_dissociation_1
        stret1_d3hh    = source_brh_dissociation_3
        stret1_ahh     = source_brh_attenuation
        bonds_r0cb     = source_cb_equilibrium
        stret1_d1cb    = source_cb_dissociation_1
        stret1_d3cb    = source_cb_dissociation_3
        stret1_acb     = source_cb_attenuation
        fsw1_a3s       = source_out_switch_a
        fsw1_b3s       = source_out_switch_b
        sphi1_aphi     = source_phi_a
        sphi1_bphi     = source_phi_b
        sphi1_cphi     = source_phi_c
        stheta1_atheta = source_theta_a
        stheta1_btheta = source_theta_b
        stheta1_ctheta = source_theta_c
        params_fch3    = source_fch3
        params_hch3    = source_hch3
        force2_fkinf   = source_force_infinite
        force2_ak      = source_force_attenuation
        force2_bk      = source_force_b
        force2_aa1     = source_force_shape_1
        force2_aa2     = source_force_shape_2
        force2_aa3     = source_force_shape_3
        force2_aa4     = source_force_shape_4
        v              = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        call flatten_cartesian_coordinates(x, tx)
        call surf(v, tx, dx, 18)
        v  = v*hartree_to_electronvolt
        dx = dx*hartree_per_bohr_to_ev_per_angstrom
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = v
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = v
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    j               = (iatom-1)*3+idir
                    g(1,iatom,idir) = dx(j)
                enddo
            enddo
        elseif( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif

        contains

        !> @brief Evaluate SURF for the CH4Br_VBMM_2002 potential.
        subroutine surf( v, coord, dx, n3tm )
            !         CH4+Br based on the surface by
            !         Jordan and Gilbert, JCP, 102, 5669 (1995).
            !         JEG, 2002
            !        This potential is written such that:
            !                       X(1)  - X(3)  : X, Y, Z for H1
            !                       X(4)  - X(6)  : X, Y, Z for C
            !                       X(7)  - X(9)  : X, Y, Z for H3
            !                       X(10) - X(12) : X, Y, Z for H4
            !                       X(13) - X(15) : X, Y, Z for H2
            !                       X(16) - X(18) : X, Y, Z for Br
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: coord(n3tm), dx(n3tm), en, fact1, fact2, v, vip, vop, vstr
            integer  :: i, icount, ind, n, n3tm
            !  include file for CH4Br PES
            do i=1, 150
                qpdot_q(i)    = 0.e0_wp
                qpdot_pdot(i) = 0.e0_wp
            enddo
            do i = 1, 18
                qpdot_q(i) = coord(i)
            enddo
            !  calculate indexes for coordinates
            do ind=1, 3
                icount       = ind-3
                ndx_nc(ind)  = 3*ndx_nnc+icount
                ndx_nhb(ind) = 3*ndx_nnb+icount
                do i=1, 4
                    ndx_nh(i,ind) = 3*ndx_nnh(i)+icount
                enddo
            enddo
            !  convert to appropriate units:
            !  kcal/mol
            !  mdyn A-1        -> 1.0d+05 j/mol...
            !  mdyn A rad-1
            !  NB integration units are:
            !  energy   in 1.0d+05 j/mol
            !  time     in 1.0d-14 s
            !  distance in 1.0d-10 m
            !  angles   in radians
            !  mass     in amu
            fact1        = 0.041840e0_wp
            fact2        = 6.022045e0_wp
            stret1_d1ch  = stret1_d1ch*fact1
            stret1_d3ch  = stret1_d3ch*fact1
            stret1_d1cb  = stret1_d1cb*fact1
            stret1_d3cb  = stret1_d3cb*fact1
            stret1_d1hh  = stret1_d1hh*fact1
            stret1_d3hh  = stret1_d3hh*fact1
            params_fch3  = params_fch3*fact2
            params_hch3  = params_hch3*fact2
            force2_fkinf = force2_fkinf*fact2
            force2_ak    = force2_ak*fact2
            !  calculate relative coordinates and bond lengths
            en = 0.0e0_wp
            call coorden
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
            !  nb: total potential energy is vstr+vop+vip
            en = vstr+vop+vip
            ! 0.03812 factor from 10(5) j/mol to au
            en = en*0.03812e0_wp
            v  = en
            !  copy the pdots calculated using MARINER code to the
            !  appropriate elements in dv. dx en POLYRATE.
            do n=1, n3tm
                dx(n) = 0.0e0_wp
            enddo
            ! 0.0201723 factor from 10(5)j/mol/A to hartrees/bohr
            do i=1, 18
                dx(i) = qpdot_pdot(i)*0.0201723e0_wp
            enddo
        end subroutine surf

        !> @brief Evaluate coorden for the CH4Br_VBMM_2002 potential.
        subroutine coorden
            use chempotpy_kinds, only: wp
            implicit none
            integer :: i, ind
            !  calculates relative coordinates and bond lengths
            !  calculate relative coordinates
            do ind=1, 3
                coords_tcb(ind) = qpdot_q(ndx_nc(ind))-qpdot_q(ndx_nhb(ind))
                do i=1, 4
                    coords_tch(i,ind) = qpdot_q(ndx_nc(ind))-qpdot_q(ndx_nh(i,ind))
                    coords_tbh(i,ind) = qpdot_q(ndx_nhb(ind))-qpdot_q(ndx_nh(i,ind))
                enddo
            enddo
            !  calculate bond lengths
            bonds_rcb = sqrt(coords_tcb(1)*coords_tcb(1)+coords_tcb(2)*coords_tcb(2)+coords_tcb(3)*coords_tcb(3))
            do i=1, 4
                bonds_rch(i) = sqrt(coords_tch(i,1)*coords_tch(i,1)+coords_tch(i,2)*coords_tch(i,2)+ coords_tch(i,3)* &
                    & coords_tch(i,3))
                bonds_rbh(i) = sqrt(coords_tbh(i,1)*coords_tbh(i,1)+coords_tbh(i,2)*coords_tbh(i,2)+ coords_tbh(i,3)* &
                    & coords_tbh(i,3))
            enddo
        end subroutine coorden

        !> @brief Evaluate refangles for the CH4Br_VBMM_2002 potential.
        subroutine refangles
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: ddr(4,4), halfpi, pi, sumd2(4), sumd4(4), tau, twopi
            integer  :: i, j, k
            !  subroutine calculates reference angles for the "in-plane" potential
            tau    = acos(-1.0e0_wp/3.0e0_wp)
            pi     = 4.0e0_wp*atan(1.0e0_wp)
            halfpi = 0.5e0_wp*pi
            twopi  = 2.0e0_wp*pi
            !  set diagonal elements to zero
            do i=1, 4
                angles_theta0(i,i) = 0.0e0_wp
                do k=1, 4
                    angles_dtheta0(i,i,k) = 0.0e0_wp
                enddo
            enddo
            !  calculate reference angles
            angles_theta0(1,2) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(2)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(3)*switch1_stheta(4)-1.0e0_wp)
            angles_theta0(1,3) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(3)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(2)*switch1_stheta(4)-1.0e0_wp)
            angles_theta0(1,4) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(4)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(2)*switch1_stheta(3)-1.0e0_wp)
            angles_theta0(2,3) = tau+(tau-halfpi)*(switch1_sphi(2)*switch1_sphi(3)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(1)*switch1_stheta(4)-1.0e0_wp)
            angles_theta0(2,4) = tau+(tau-halfpi)*(switch1_sphi(2)*switch1_sphi(4)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(1)*switch1_stheta(3)-1.0e0_wp)
            angles_theta0(3,4) = tau+(tau-halfpi)*(switch1_sphi(3)*switch1_sphi(4)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(1)*switch1_stheta(2)-1.0e0_wp)
            !  calculate the derivatives of theta0(i,j) in terms of rch(k)
            !  quantity calulated is dtheta0(i,j,k)
            !  derivatives wrt rch(1)
            angles_dtheta0(1,2,1) = (tau-halfpi)*switch1_dsphi(1)*switch1_sphi(2)
            angles_dtheta0(1,3,1) = (tau-halfpi)*switch1_dsphi(1)*switch1_sphi(3)
            angles_dtheta0(1,4,1) = (tau-halfpi)*switch1_dsphi(1)*switch1_sphi(4)
            angles_dtheta0(2,3,1) = (tau-twopi/3.0e0_wp)*switch1_dstheta(1)*switch1_stheta(4)
            angles_dtheta0(2,4,1) = (tau-twopi/3.0e0_wp)*switch1_dstheta(1)*switch1_stheta(3)
            angles_dtheta0(3,4,1) = (tau-twopi/3.0e0_wp)*switch1_dstheta(1)*switch1_stheta(2)
            !  derivatives wrt rch(2)
            angles_dtheta0(1,2,2) = (tau-halfpi)*switch1_sphi(1)*switch1_dsphi(2)
            angles_dtheta0(1,3,2) = (tau-twopi/3.0e0_wp)*switch1_dstheta(2)*switch1_stheta(4)
            angles_dtheta0(1,4,2) = (tau-twopi/3.0e0_wp)*switch1_dstheta(2)*switch1_stheta(3)
            angles_dtheta0(2,3,2) = (tau-halfpi)*switch1_dsphi(2)*switch1_sphi(3)
            angles_dtheta0(2,4,2) = (tau-halfpi)*switch1_dsphi(2)*switch1_sphi(4)
            angles_dtheta0(3,4,2) = (tau-twopi/3.0e0_wp)*switch1_stheta(1)*switch1_dstheta(2)
            !  derivatives wrt rch(3)
            angles_dtheta0(1,2,3) = (tau-twopi/3.0e0_wp)*switch1_dstheta(3)*switch1_stheta(4)
            angles_dtheta0(1,3,3) = (tau-halfpi)*switch1_sphi(1)*switch1_dsphi(3)
            angles_dtheta0(1,4,3) = (tau-twopi/3.0e0_wp)*switch1_stheta(2)*switch1_dstheta(3)
            angles_dtheta0(2,3,3) = (tau-halfpi)*switch1_sphi(2)*switch1_dsphi(3)
            angles_dtheta0(2,4,3) = (tau-twopi/3.0e0_wp)*switch1_stheta(1)*switch1_dstheta(3)
            angles_dtheta0(3,4,3) = (tau-halfpi)*switch1_dsphi(3)*switch1_sphi(4)
            !  derivatives wrt rch(4)
            angles_dtheta0(1,2,4) = (tau-twopi/3.0e0_wp)*switch1_stheta(3)*switch1_dstheta(4)
            angles_dtheta0(1,3,4) = (tau-twopi/3.0e0_wp)*switch1_stheta(2)*switch1_dstheta(4)
            angles_dtheta0(1,4,4) = (tau-halfpi)*switch1_sphi(1)*switch1_dsphi(4)
            angles_dtheta0(2,3,4) = (tau-twopi/3.0e0_wp)*switch1_stheta(1)*switch1_dstheta(4)
            angles_dtheta0(2,4,4) = (tau-halfpi)*switch1_sphi(2)*switch1_dsphi(4)
            angles_dtheta0(3,4,4) = (tau-halfpi)*switch1_sphi(3)*switch1_dsphi(4)
            !  fill in the other half of the matrix
            do i=1, 3
                do j=i+1, 4
                    angles_theta0(j,i) = angles_theta0(i,j)
                    do k=1, 4
                        angles_dtheta0(j,i,k) = angles_dtheta0(i,j,k)
                    enddo
                enddo
            enddo
        end subroutine refangles

        !> @brief Evaluate stretch for the CH4Br_VBMM_2002 potential.
        subroutine stretch( vstr )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: ach, achdc(3), achdh(4,3), arga, dumach, dumjbh, dumjcb, dumjch, dumjhi, dumjhj, dumqbh, dumqcb, &
                & dumqch, dumqhi, dumqhj, e1, e3, factj, rav, vj(4), vjbh(4), vjcb, vjch(4), vq(4), vqbh(4), vqcb, vqch(4), vstr
            integer :: i, ind, j, k
            !  subroutine to calculate leps-type stretching potential and its
            !  derivatives
            !  calculate avergage bond length for the methane moiety
            rav = (bonds_rch(1)+bonds_rch(2)+bonds_rch(3)+bonds_rch(4))/4.0e0_wp
            !  initialise:
            vstr = 0.0e0_wp
            !  ach:
            !  nb: in double precision tanh(19.0d0)=1.0d0 and we put the if statemen
            !  in to avoid overflow/underflow errors
            arga = stret1_c1ch*(rav-bonds_r0ch)
            if( arga .lt. 19.0e0_wp )then
                ach    = stret1_a1ch+stret1_b1ch*(tanh(arga)+1.0e0_wp)*0.5e0_wp
                dumach = stret1_b1ch*stret1_c1ch/(2.0e0_wp*cosh(arga)**2)
            else
                ach    = stret1_a1ch+stret1_b1ch
                dumach = 0.0e0_wp
            endif
            !  calculate singlet: e1, triplet: e3 energies and vq and vj
            !  terms for each bond
            e1   = stret1_d1cb*(exp(-2.0e0_wp*stret1_acb*(bonds_rcb-bonds_r0cb))-2.0e0_wp*exp(-stret1_acb*(bonds_rcb- &
                & bonds_r0cb)))
            e3   = stret1_d3cb*(exp(-2.0e0_wp*stret1_acb*(bonds_rcb-bonds_r0cb))+2.0e0_wp*exp(-stret1_acb*(bonds_rcb- &
                & bonds_r0cb)))
            vqcb = (e1+e3)*0.5e0_wp
            vjcb = (e1-e3)*0.5e0_wp
            do i=1, 4
                e1      = stret1_d1ch*(exp(-2.0e0_wp*ach*(bonds_rch(i)-bonds_r0ch)) -2.0e0_wp*exp(-ach*(bonds_rch(i)- &
                    & bonds_r0ch)))
                e3      = stret1_d3ch*(exp(-2.0e0_wp*ach*(bonds_rch(i)-bonds_r0ch)) +2.0e0_wp*exp(-ach*(bonds_rch(i)- &
                    & bonds_r0ch)))
                vqch(i) = (e1+e3)*0.5e0_wp
                vjch(i) = (e1-e3)*0.5e0_wp
                e1      = stret1_d1hh*(exp(-2.0e0_wp*stret1_ahh*(bonds_rbh(i)-bonds_r0hh)) -2.0e0_wp*exp(-stret1_ahh* &
                    & (bonds_rbh(i)-bonds_r0hh)))
                e3      = stret1_d3hh*(exp(-2.0e0_wp*stret1_ahh*(bonds_rbh(i)-bonds_r0hh)) +2.0e0_wp*exp(-stret1_ahh* &
                    & (bonds_rbh(i)-bonds_r0hh)))
                vqbh(i) = (e1+e3)*0.5e0_wp
                vjbh(i) = (e1-e3)*0.5e0_wp
                !  calculate 3 body potential
                vq(i) = vqch(i)+vqcb+vqbh(i)
                vj(i) = -sqrt(((vjch(i)-vjcb)**2+(vjcb-vjbh(i))**2 +(vjbh(i)-vjch(i))**2)*0.5e0_wp)
                vstr  = vstr+vq(i)+vj(i)
            enddo
            !  partial derivatives
            !  first we need the derivative of ach:
            do ind=1, 3
                achdc(ind) = dumach*(coords_tch(1,ind)/bonds_rch(1)+coords_tch(2,ind)/bonds_rch(2) +coords_tch(3, &
                    & ind)/bonds_rch(3)+coords_tch(4,ind)/bonds_rch(4))/4.0e0_wp
                do i=1, 4
                    achdh(i,ind) = -dumach*coords_tch(i,ind)/bonds_rch(i)/4.0e0_wp
                enddo
            enddo
            dumqcb = -stret1_acb*((stret1_d1cb+stret1_d3cb)*exp(-2.0e0_wp*stret1_acb*(bonds_rcb-bonds_r0cb))- &
                & (stret1_d1cb-stret1_d3cb)*exp(-stret1_acb*(bonds_rcb-bonds_r0cb)))/bonds_rcb
            !  calculate cartesian derivatives:
            !  looping over ch(i) and bh(i)
            do i=1, 4
                dumqbh = -stret1_ahh*((stret1_d1hh+stret1_d3hh)*exp(-2.0e0_wp*stret1_ahh*(bonds_rbh(i)-bonds_r0hh))- &
                    & (stret1_d1hh-stret1_d3hh)*exp(-stret1_ahh*(bonds_rbh(i)-bonds_r0hh)))/bonds_rbh(i)
                factj = 0.5e0_wp/vj(i)
                dumjcb = -stret1_acb*((stret1_d1cb-stret1_d3cb)*exp(-2.0e0_wp*stret1_acb*(bonds_rcb-bonds_r0cb)) - &
                    & (stret1_d1cb+stret1_d3cb)*exp(-stret1_acb*(bonds_rcb-bonds_r0cb)))*factj/bonds_rcb
                dumjbh = -stret1_ahh*((stret1_d1hh-stret1_d3hh)*exp(-2.0e0_wp*stret1_ahh*(bonds_rbh(i)-bonds_r0hh)) - &
                    & (stret1_d1hh+stret1_d3hh)*exp(-stret1_ahh*(bonds_rbh(i)-bonds_r0hh)))*factj/bonds_rbh(i)
                do ind=1, 3
                    !  deriv wrt hb:
                    qpdot_pdot(ndx_nhb(ind))=qpdot_pdot(ndx_nhb(ind)) -coords_tcb(ind)*dumqcb+coords_tbh(i,ind)* &
                        & dumqbh +(vjch(i)-vjcb)*(dumjcb*coords_tcb(ind)) + (vjcb- vjbh(i))*(-dumjcb* coords_tcb(ind)- &
                        & dumjbh*coords_tbh(i,ind)) +(vjbh(i)-vjch(i))*dumjbh*coords_tbh(i,ind)
                    !  dvqch(i)/dc
                    dumqch=-(ach*coords_tch(i,ind)/bonds_rch(i)+achdc(ind)*(bonds_rch(i)-bonds_r0ch)) * ((stret1_d1ch+ &
                        & stret1_d3ch)*exp(-2.0e0_wp*ach*(bonds_rch(i)- bonds_r0ch)) - (stret1_d1ch- stret1_d3ch)* &
                        & exp(- ach*(bonds_rch(i)-bonds_r0ch)))
                    qpdot_pdot(ndx_nc(ind)) = qpdot_pdot(ndx_nc(ind))+dumqch+coords_tcb(ind)*dumqcb
                    !  dvqch(i)/dh(i)
                    dumqhi=(ach*coords_tch(i,ind)/bonds_rch(i)-achdh(i,ind)*(bonds_rch(i)-bonds_r0ch)) * &
                        & ((stret1_d1ch+stret1_d3ch)*exp(-2.0e0_wp*ach*(bonds_rch(i)- bonds_r0ch)) - (stret1_d1ch- &
                        & stret1_d3ch)*exp(-ach*(bonds_rch(i)-bonds_r0ch)))
                    qpdot_pdot(ndx_nh(i,ind)) = qpdot_pdot(ndx_nh(i,ind))+dumqhi-coords_tbh(i,ind)*dumqbh
                    !  dvjch(i)/dc
                    dumjch=-(ach*coords_tch(i,ind)/bonds_rch(i)+achdc(ind)*(bonds_rch(i)-bonds_r0ch)) * ((stret1_d1ch- &
                        & stret1_d3ch)*exp(-2.0e0_wp*ach*(bonds_rch(i)- bonds_r0ch)) - (stret1_d1ch+ stret1_d3ch)* &
                        & exp(- ach*(bonds_rch(i)-bonds_r0ch)))*factj
                    !  dvj(i)/dnc(ind)
                    qpdot_pdot(ndx_nc(ind))=qpdot_pdot(ndx_nc(ind)) +(vjch(i)-vjcb)*(dumjch-dumjcb* coords_tcb(ind)) + &
                        & (vjcb-vjbh(i))*dumjcb*coords_tcb(ind) - (vjbh(i)-vjch(i))*dumjch
                    !  dvjch(i)/dh(i)
                    dumjhi=(ach*coords_tch(i,ind)/bonds_rch(i)-achdh(i,ind)*(bonds_rch(i)-bonds_r0ch)) * &
                        & ((stret1_d1ch-stret1_d3ch)*exp(-2.0e0_wp*ach*(bonds_rch(i)- bonds_r0ch)) - (stret1_d1ch+ &
                        & stret1_d3ch)*exp(-ach*(bonds_rch(i)-bonds_r0ch)))*factj
                    !  dvj(i)/dnh(i,ind)
                    qpdot_pdot(ndx_nh(i,ind))=qpdot_pdot(ndx_nh(i,ind)) +(vjch(i)-vjcb)*dumjhi +(vjcb-vjbh(i))* &
                        & dumjbh*coords_tbh(i,ind) +(vjbh(i)- vjch(i))* (-dumjbh*coords_tbh(i,ind)-dumjhi)
                    !  dv(i)/dh(j)
                    do k=1, 3
                        j = i+k
                        if( j .gt. 4 ) j=j-4
                        dumqhj=-achdh(j,ind)*(bonds_rch(i)-bonds_r0ch) *((stret1_d1ch+stret1_d3ch)*exp(-2.0e0_wp*ach* &
                            & (bonds_rch(i)-bonds_r0ch)) -(stret1_d1ch-stret1_d3ch)*exp(- ach* (bonds_rch(i)- &
                            & bonds_r0ch)))
                        dumjhj=-achdh(j,ind)*(bonds_rch(i)-bonds_r0ch) *((stret1_d1ch-stret1_d3ch)*exp(-2.0e0_wp*ach* &
                            & (bonds_rch(i)-bonds_r0ch)) -(stret1_d1ch+stret1_d3ch)*exp(- ach* (bonds_rch(i)- &
                            & bonds_r0ch)))*factj
                        qpdot_pdot(ndx_nh(j,ind)) = qpdot_pdot(ndx_nh(j,ind))+dumqhj +(vjch(i)-vjcb)*dumjhj - &
                            & (vjbh(i)-vjch(i))*dumjhj
                    enddo
                enddo
            enddo
        end subroutine stretch

        !> @brief Evaluate opbend for the CH4Br_VBMM_2002 potential.
        subroutine opbend( vop )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3), argd(4), axb(3), b(3), c(4,3), ddr, sum2, sum4, sumd2(4), sumd4(4), vop
            integer  :: i, ii, in(3), ind, itemp, j, k, l
            !  subroutine calculates symmetrized vop potential and derivatives
            real(wp) :: norma
            vop = 0.0e0_wp
            !  calculate force constants and their derivatives
            call opforce
            !  calculate out-of-plane angle and derivatives
            do i=1, 4
                j = i+1
                if( j .gt. 4 ) j=j-4
                k = j+1
                if( k .gt. 4 ) k=k-4
                l = k+1
                if( l .gt. 4 ) l=l-4
                !  jcc 6/4/99
                !  modification to ensure that the set of methane CH bond vector
                !  (rj,rk,rl) is a right-handed set
                in(1) = j
                in(2) = k
                in(3) = l
                !  vector a is rk-rj, vector b is rl-rj
                do ind=1, 3
                    a(ind) = qpdot_q(ndx_nh(k,ind))-qpdot_q(ndx_nh(j,ind))
                    b(ind) = qpdot_q(ndx_nh(l,ind))-qpdot_q(ndx_nh(j,ind))
                enddo
                !  axb is vector a cross b
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
                !  c is position vector of h(ii): calculate c(j),c(k),c(l)
                do ii=1, 3
                    do ind=1, 3
                        c(in(ii),ind) = -coords_tch(in(ii),ind)/bonds_rch(in(ii))
                    enddo
                enddo
                !  argd is the dot product axb dot c
                do ii=1, 3
                    argd(in(ii)) = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                    argd(in(ii)) = argd(in(ii))/norma
                    !  if argd > 0 we need to switch vectors k and l around
                    if( argd(in(ii)) .gt. 0.e0_wp )then
                        itemp = k
                        k     = l
                        l     = itemp
                    endif
                enddo
                !  jcc 6/4/99
                !  if i is an even number then we have switched the vectors
                !  from a right handed set to a left handed set
                !  in this case we need to switch vectors k and l around
                !  subroutine performs sum over j, k, l
                !  sum2 = sum delta**2
                !  sum4 = sum delta**4
                call calcdelta(i, j, k, l, sum2, sum4)
                sumd2(i) = sum2
                sumd4(i) = sum4
                vop      = vop+delta1_fdelta(i)*sumd2(i)+delta1_hdelta(i)*sumd4(i)
            enddo
            do i=1, 4
                do j=1, 4
                    !  overall derivatives of force constants i wrt the bond-length rch(j)
                    ddr = delta2_dfdelta(i,j)*sumd2(i)+delta2_dhdelta(i,j)*sumd4(i)
                    !  calculate derivatives in terms of cartesian coordinates:
                    do ind=1, 3
                        qpdot_pdot(ndx_nh(j,ind)) = qpdot_pdot(ndx_nh(j,ind))-coords_tch(j,ind)*ddr/bonds_rch(j)
                        qpdot_pdot(ndx_nc(ind))   = qpdot_pdot(ndx_nc(ind))+coords_tch(j,ind)*ddr/bonds_rch(j)
                    enddo
                enddo
            enddo
        end subroutine opbend

        !> @brief Evaluate ipbend for the CH4Br_VBMM_2002 potential.
        subroutine ipbend( vip )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: costh(4,4), dth(4,4), dth0c, dth0k, dthc, dthi, dthj, termth, theta(4,4), vip
            integer  :: i, ind, j, k
            !  subroutine calculates symmetrised in plane bend term
            !  and its derivatives
            !  initialise
            vip = 0.0e0_wp
            !  calculate force constants: fk0(i,j), f1(i)
            !  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
            call ipforce
            !  calculate theta(i,j) and in plane bend potential
            do i=1, 3
                do j=i+1, 4
                    costh(i,j) = coords_tch(i,1)*coords_tch(j,1)+coords_tch(i,2)*coords_tch(j,2) +coords_tch(i,3)* &
                        & coords_tch(j,3)
                    costh(i,j) = costh(i,j)/bonds_rch(i)/bonds_rch(j)
                    theta(i,j) = acos(costh(i,j))
                    dth(i,j)   = theta(i,j)-angles_theta0(i,j)
                    vip        = vip+0.5e0_wp*force1_fk0(i,j)*force1_f1(i)*force1_f1(j)*dth(i,j)**2
                    !  calculate partial derivatives wrt cartesian coordinates
                    !  calculate pdots wrt theta:
                    termth = -1.0e0_wp/sqrt(1.0e0_wp-costh(i,j)*costh(i,j))
                    do ind=1, 3
                        dthi             = -coords_tch(j,ind)/bonds_rch(i)/bonds_rch(j) +costh(i,j)*coords_tch(i, &
                            & ind)/bonds_rch(i)/bonds_rch(i)
                        dthi = dthi*termth
                        dthj             = -coords_tch(i,ind)/bonds_rch(i)/bonds_rch(j) +costh(i,j)*coords_tch(j, &
                            & ind)/bonds_rch(j)/bonds_rch(j)
                        dthj = dthj*termth
                        dthc = -(dthi+dthj)
                        qpdot_pdot(ndx_nh(i,ind)) = qpdot_pdot(ndx_nh(i,ind)) +force1_fk0(i,j)*force1_f1(i)* &
                            & force1_f1(j)*dthi*dth(i,j)
                        qpdot_pdot(ndx_nh(j,ind)) = qpdot_pdot(ndx_nh(j,ind)) +force1_fk0(i,j)*force1_f1(i)* &
                            & force1_f1(j)*dthj*dth(i,j)
                        qpdot_pdot(ndx_nc(ind))    = qpdot_pdot(ndx_nc(ind)) +force1_fk0(i,j)*force1_f1(i)* &
                            & force1_f1(j)*dthc*dth(i,j)
                        do k=1, 4
                            !  calculate pdots wrt force constants and wrt theta0
                            dth0k = -angles_dtheta0(i,j,k)*coords_tch(k,ind)/bonds_rch(k)
                            dth0c = -dth0k
                            qpdot_pdot(ndx_nh(k,ind))=qpdot_pdot(ndx_nh(k,ind)) -0.5e0_wp*coords_tch(k,ind)* &
                                & force1_dfdc(i,j,k)*dth(i,j)**2/bonds_rch(k) -0.5e0_wp* coords_tbh(k,ind)* &
                                & force1_dfdh(i,j,k)*dth(i,j)**2/bonds_rbh(k) -force1_fk0(i,j)*force1_f1(i)* &
                                & force1_f1(j)*dth0k*dth(i,j)
                            qpdot_pdot(ndx_nc(ind))=qpdot_pdot(ndx_nc(ind)) +0.5e0_wp*coords_tch(k,ind)* &
                                & force1_dfdc(i,j,k)*dth(i,j)**2/bonds_rch(k) -force1_fk0(i,j)* force1_f1(i)* &
                                & force1_f1(j)* dth0c* dth(i,j)
                            qpdot_pdot(ndx_nhb(ind)) = qpdot_pdot(ndx_nhb(ind)) +0.5e0_wp*coords_tbh(k,ind)* &
                                & force1_dfdh(i,j,k)*dth(i,j)**2/bonds_rbh(k)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine ipbend

        !> @brief Evaluate calcdelta for the CH4Br_VBMM_2002 potential.
        subroutine calcdelta( i, j, k, l, sum2, sum4 )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3), argd(4), atemp1, atemp2(3), atemp3, atemp4, atemp5, axb(3), b(3), c(4,3), cdot(4,3,3), &
                & daxb(4, 3, 3), deldot, delta(4), sum2, sum4
            integer :: i, ii, in(3), ind, j, jj, k, l
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            real(wp) :: norma
            !  initialise
            sum2 = 0.0e0_wp; sum4 = 0.0e0_wp
            !  set j,k,l indices
            in(1) = j
            in(2) = k
            in(3) = l
            !  vector a is rk-rj, vector b is rl-rj
            do ind=1, 3
                a(ind) = qpdot_q(ndx_nh(k,ind))-qpdot_q(ndx_nh(j,ind))
                b(ind) = qpdot_q(ndx_nh(l,ind))-qpdot_q(ndx_nh(j,ind))
            enddo
            !  axb is vector a cross b
            axb(1) = a(2)*b(3)-a(3)*b(2)
            axb(2) = a(3)*b(1)-a(1)*b(3)
            axb(3) = a(1)*b(2)-a(2)*b(1)
            norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
            norma  = sqrt(norma)
            !  c is position vector of h(ii): calculate c(j),c(k),c(l)
            do ii=1, 3
                do ind=1, 3
                    c(in(ii),ind) = -coords_tch(in(ii),ind)/bonds_rch(in(ii))
                enddo
            enddo
            !  argd is the dot product axb dot c
            do ii=1, 3
                argd(in(ii))  = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                argd(in(ii))  = argd(in(ii))/norma
                delta(in(ii)) = acos(argd(in(ii)))-angles_theta0(i,in(ii))
                sum2          = sum2+delta(in(ii))**2
                sum4          = sum4+delta(in(ii))**4
            enddo
            !  derivatives of axb wrt hj:
            daxb(j, 1, 1) = 0.0e0_wp
            daxb(j, 1, 2) = b(3)-a(3)
            daxb(j, 1, 3) = -b(2)+a(2)
            daxb(j, 2, 1) = -b(3)+a(3)
            daxb(j, 2, 2) = 0.0e0_wp
            daxb(j, 2, 3) = b(1)-a(1)
            daxb(j, 3, 1) = b(2)-a(2)
            daxb(j, 3, 2) = -b(1)+a(1)
            daxb(j, 3, 3) = 0.0e0_wp
            !  derivatives of axb wrt hk:
            daxb(k, 1, 1) = 0.0e0_wp
            daxb(k, 1, 2) = -b(3)
            daxb(k, 1, 3) = b(2)
            daxb(k, 2, 1) = b(3)
            daxb(k, 2, 2) = 0.0e0_wp
            daxb(k, 2, 3) = -b(1)
            daxb(k, 3, 1) = -b(2)
            daxb(k, 3, 2) = b(1)
            daxb(k, 3, 3) = 0.0e0_wp
            !  derivatives of axb wrt hl:
            daxb(l, 1, 1) = 0.0e0_wp
            daxb(l, 1, 2) = a(3)
            daxb(l, 1, 3) = -a(2)
            daxb(l, 2, 1) = -a(3)
            daxb(l, 2, 2) = 0.0e0_wp
            daxb(l, 2, 3) = a(1)
            daxb(l, 3, 1) = a(2)
            daxb(l, 3, 2) = -a(1)
            daxb(l, 3, 3) = 0.0e0_wp
            !   loop over cdot(in(ii),ind,jind) where we consider deriv of c(in(ii))
            !   wrt h(in(ii),jind) with components jind
            do ii=1, 3
                !  deriv of cdot(in(ii),x) wrt x, y, z
                cdot(in(ii),1,1) = 1.0e0_wp/bonds_rch(in(ii)) +coords_tch(in(ii),1)*c(in(ii),1)/bonds_rch(in(ii))**2
                cdot(in(ii),1,2) = coords_tch(in(ii),2)*c(in(ii),1)/bonds_rch(in(ii))**2
                cdot(in(ii),1,3) = coords_tch(in(ii),3)*c(in(ii),1)/bonds_rch(in(ii))**2
                !  deriv of cdot(in(ii),y) wrt x, y, z
                cdot(in(ii),2,1) = coords_tch(in(ii),1)*c(in(ii),2)/bonds_rch(in(ii))**2
                cdot(in(ii),2,2) = 1.0e0_wp/bonds_rch(in(ii)) +coords_tch(in(ii),2)*c(in(ii),2)/bonds_rch(in(ii))**2
                cdot(in(ii),2,3) = coords_tch(in(ii),3)*c(in(ii),2)/bonds_rch(in(ii))**2
                !  deriv of cdot(in(ii),z) wrt x, y, z
                cdot(in(ii),3,1) = coords_tch(in(ii),1)*c(in(ii),3)/bonds_rch(in(ii))**2
                cdot(in(ii),3,2) = coords_tch(in(ii),2)*c(in(ii),3)/bonds_rch(in(ii))**2
                cdot(in(ii),3,3) = 1.0e0_wp/bonds_rch(in(ii)) +coords_tch(in(ii),3)*c(in(ii),3)/bonds_rch(in(ii))**2
            enddo
            do ii=1, 3
                do ind=1, 3
                    deldot = -angles_dtheta0(i,in(ii),i)
                    !  derivative wrt h(i,ind)
                    !  for  rch(i) only terms are from the derivatives of theta0
                    deldot = -deldot*coords_tch(i,ind)/bonds_rch(i)
                    qpdot_pdot(ndx_nh(i,ind))=qpdot_pdot(ndx_nh(i,ind)) +2.0e0_wp*delta1_fdelta(i)*delta(in(ii))* &
                        & deldot +4.0e0_wp*delta1_hdelta(i)* delta(in(ii))**3* deldot
                    !  derivative wrt c(ind)
                    deldot = -deldot
                    qpdot_pdot(ndx_nc(ind)) = qpdot_pdot(ndx_nc(ind)) +2.0e0_wp*delta1_fdelta(i)*delta(in(ii))* &
                        & deldot + 4.0e0_wp*delta1_hdelta(i)*delta(in(ii))**3* deldot
                    do jj=1, 3
                        !  partial derivatives wrt h(in(jj),ind), loop over delta(i,in(ii))
                        !   atemp1 is axb dot daxb wrt h(in(jj))
                        atemp1 = axb(1)*daxb(in(jj), ind, 1) +axb(2)*daxb(in(jj), ind, 2) +axb(3)*daxb(in(jj), ind, 3)
                        atemp1 = atemp1/(norma**3)
                        !  atemp2 is deriv of normalised axb
                        atemp2(1) = daxb(in(jj), ind, 1)/norma-atemp1*axb(1)
                        atemp2(2) = daxb(in(jj), ind, 2)/norma-atemp1*axb(2)
                        atemp2(3) = daxb(in(jj), ind, 3)/norma-atemp1*axb(3)
                        !  atemp3 is daxb dot c(in(ii))
                        atemp3 = atemp2(1)*c(in(ii),1)+atemp2(2)*c(in(ii),2) +atemp2(3)*c(in(ii),3)
                        !  atemp4 is axb dot cdot
                        atemp4 = 0.0e0_wp
                        if( ii .eq. jj )then
                            !  ie deriv of c(in(ii)) wrt h(in(jj)) is non zero only for ii = jj
                            atemp4 = axb(1)*cdot(in(ii),1,ind) +axb(2)*cdot(in(ii),2,ind) +axb(3)*cdot(in(ii),3,ind)
                            atemp4 = atemp4/norma
                        endif
                        !  atemp5 is deriv of theta0(i,in(ii)) wrt to nh(in(jj),ind)
                        atemp5 = -angles_dtheta0(i,in(ii),in(jj))
                        !  deriv wrt h(in(jj)),ind):
                        atemp5 = -atemp5*coords_tch(in(jj),ind)/bonds_rch(in(jj))
                        deldot = atemp3+atemp4
                        deldot = -1.0e0_wp/sqrt(1.0e0_wp-argd(in(ii))**2)*deldot
                        deldot = deldot+atemp5
                        qpdot_pdot(ndx_nh(in(jj),ind))=qpdot_pdot(ndx_nh(in(jj),ind)) +2.0e0_wp*delta1_fdelta(i)* &
                            & delta(in(ii))*deldot +4.0e0_wp* delta1_hdelta(i)* delta(in(ii))**3*deldot
                        !  for carbon the only contributions are from axb dot cdot term and
                        !  from theta0 and derivative cdot wrt carbon=-cdot wrt hydrogen
                        deldot = 1.0e0_wp/sqrt(1.0e0_wp-argd(in(ii))**2)*atemp4
                        deldot = deldot-atemp5
                        qpdot_pdot(ndx_nc(ind))=qpdot_pdot(ndx_nc(ind)) +2.0e0_wp*delta1_fdelta(i)*delta(in(ii))* &
                            & deldot +4.0e0_wp*delta1_hdelta(i)*delta(in(ii))** 3* deldot
                    enddo
                enddo
            enddo
        end subroutine calcdelta

        !> @brief Evaluate opforce for the CH4Br_VBMM_2002 potential.
        subroutine opforce
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dswitch(4,4), switch(4)
            integer  :: i, j
            !  calculates the out-of-plane bending force constants
            !  and their derivatives
            !  calculate switching functions:
            switch(1) = (1.0e0_wp-op1_s3(1))*op1_s3(2)*op1_s3(3)*op1_s3(4)
            switch(2) = (1.0e0_wp-op1_s3(2))*op1_s3(3)*op1_s3(4)*op1_s3(1)
            switch(3) = (1.0e0_wp-op1_s3(3))*op1_s3(4)*op1_s3(1)*op1_s3(2)
            switch(4) = (1.0e0_wp-op1_s3(4))*op1_s3(1)*op1_s3(2)*op1_s3(3)
            !  calculate derivatives:
            !  derivative of switch(1) wrt the 4 rch bond lengths
            dswitch(1,1) = -op1_ds3(1)*op1_s3(2)*op1_s3(3)*op1_s3(4)
            dswitch(1,2) = (1.0e0_wp-op1_s3(1))*op1_ds3(2)*op1_s3(3)*op1_s3(4)
            dswitch(1,3) = (1.0e0_wp-op1_s3(1))*op1_s3(2)*op1_ds3(3)*op1_s3(4)
            dswitch(1,4) = (1.0e0_wp-op1_s3(1))*op1_s3(2)*op1_s3(3)*op1_ds3(4)
            !  derivative of switch(2) wrt the 4 rch bond lengths
            dswitch(2,1) = (1.0e0_wp-op1_s3(2))*op1_s3(3)*op1_s3(4)*op1_ds3(1)
            dswitch(2,2) = -op1_ds3(2)*op1_s3(3)*op1_s3(4)*op1_s3(1)
            dswitch(2,3) = (1.0e0_wp-op1_s3(2))*op1_ds3(3)*op1_s3(4)*op1_s3(1)
            dswitch(2,4) = (1.0e0_wp-op1_s3(2))*op1_s3(3)*op1_ds3(4)*op1_s3(1)
            !  derivative of switch(3) wrt the 4 rch bond lengths
            dswitch(3,1) = (1.0e0_wp-op1_s3(3))*op1_s3(4)*op1_ds3(1)*op1_s3(2)
            dswitch(3,2) = (1.0e0_wp-op1_s3(3))*op1_s3(4)*op1_s3(1)*op1_ds3(2)
            dswitch(3,3) = -op1_ds3(3)*op1_s3(4)*op1_s3(1)*op1_s3(2)
            dswitch(3,4) = (1.0e0_wp-op1_s3(3))*op1_ds3(4)*op1_s3(1)*op1_s3(2)
            !  derivative of switch(3) wrt the 4 rch bond lengths
            dswitch(4,1) = (1.0e0_wp-op1_s3(4))*op1_ds3(1)*op1_s3(2)*op1_s3(3)
            dswitch(4,2) = (1.0e0_wp-op1_s3(4))*op1_s3(1)*op1_ds3(2)*op1_s3(3)
            dswitch(4,3) = (1.0e0_wp-op1_s3(4))*op1_s3(1)*op1_s3(2)*op1_ds3(3)
            dswitch(4,4) = -op1_ds3(4)*op1_s3(1)*op1_s3(2)*op1_s3(3)
            !  calculate the force constants and their derivatives
            do i=1, 4
                delta1_fdelta(i) = switch(i)*params_fch3
                delta1_hdelta(i) = switch(i)*params_hch3
                do j=1, 4
                    delta2_dfdelta(i,j) = dswitch(i,j)*params_fch3
                    delta2_dhdelta(i,j) = dswitch(i,j)*params_hch3
                enddo
            enddo
        end subroutine opforce

        !> @brief Evaluate ipforce for the CH4Br_VBMM_2002 potential.
        subroutine ipforce
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a1, a2, arga1, arga2, df1dc(4), df1dh(4), dfk0(4,4,4), duma1, duma2, f0, f2
            integer  :: i
            !  calculates the symmetrised in plane bend force constants and
            !  all partial derivatives involving them
            !  set force constant at asymptotes
            f0              = force2_fkinf+force2_ak
            f2              = force2_fkinf
            force1_fk0(1,2) = f0+f0*(ip1_s1(1)*ip1_s1(2)-1.0e0_wp)+(f0-f2)*(ip1_s2(3)*ip1_s2(4)-1.0e0_wp)
            force1_fk0(1,3) = f0+f0*(ip1_s1(1)*ip1_s1(3)-1.0e0_wp)+(f0-f2)*(ip1_s2(2)*ip1_s2(4)-1.0e0_wp)
            force1_fk0(1,4) = f0+f0*(ip1_s1(1)*ip1_s1(4)-1.0e0_wp)+(f0-f2)*(ip1_s2(2)*ip1_s2(3)-1.0e0_wp)
            force1_fk0(2,3) = f0+f0*(ip1_s1(2)*ip1_s1(3)-1.0e0_wp)+(f0-f2)*(ip1_s2(1)*ip1_s2(4)-1.0e0_wp)
            force1_fk0(2,4) = f0+f0*(ip1_s1(2)*ip1_s1(4)-1.0e0_wp)+(f0-f2)*(ip1_s2(1)*ip1_s2(3)-1.0e0_wp)
            force1_fk0(3,4) = f0+f0*(ip1_s1(3)*ip1_s1(4)-1.0e0_wp)+(f0-f2)*(ip1_s2(1)*ip1_s2(2)-1.0e0_wp)
            !  derivative of fk0
            dfk0(1,2,1) = f0*ip1_ds1(1)*ip1_s1(2)
            dfk0(1,2,2) = f0*ip1_s1(1)*ip1_ds1(2)
            dfk0(1,2,3) = (f0-f2)*ip1_ds2(3)*ip1_s2(4)
            dfk0(1,2,4) = (f0-f2)*ip1_s2(3)*ip1_ds2(4)
            dfk0(1,3,1) = f0*ip1_ds1(1)*ip1_s1(3)
            dfk0(1,3,2) = (f0-f2)*ip1_ds2(2)*ip1_s2(4)
            dfk0(1,3,3) = f0*ip1_s1(1)*ip1_ds1(3)
            dfk0(1,3,4) = (f0-f2)*ip1_s2(2)*ip1_ds2(4)
            dfk0(1,4,1) = f0*ip1_ds1(1)*ip1_s1(4)
            dfk0(1,4,2) = (f0-f2)*ip1_ds2(2)*ip1_s2(3)
            dfk0(1,4,3) = (f0-f2)*ip1_s2(2)*ip1_ds2(3)
            dfk0(1,4,4) = f0*ip1_s1(1)*ip1_ds1(4)
            dfk0(2,3,1) = (f0-f2)*ip1_ds2(1)*ip1_s2(4)
            dfk0(2,3,2) = f0*ip1_ds1(2)*ip1_s1(3)
            dfk0(2,3,3) = f0*ip1_s1(2)*ip1_ds1(3)
            dfk0(2,3,4) = (f0-f2)*ip1_s2(1)*ip1_ds2(4)
            dfk0(2,4,1) = (f0-f2)*ip1_ds2(1)*ip1_s2(3)
            dfk0(2,4,2) = f0*ip1_ds1(2)*ip1_s1(4)
            dfk0(2,4,3) = (f0-f2)*ip1_s2(1)*ip1_ds2(3)
            dfk0(2,4,4) = f0*ip1_s1(2)*ip1_ds1(4)
            dfk0(3,4,1) = (f0-f2)*ip1_ds2(1)*ip1_s2(2)
            dfk0(3,4,2) = (f0-f2)*ip1_s2(1)*ip1_ds2(2)
            dfk0(3,4,3) = f0*ip1_ds1(3)*ip1_s1(4)
            dfk0(3,4,4) = f0*ip1_s1(3)*ip1_ds1(4)
            !     *              +(rch(3)-r0ch)**2+(rch(4)-r0ch)**2)
            do i=1, 4
                !  calc derivatives of fk0 wrt each of the rch(i) bonds
                !  calculate the terms f1(i)
                arga1        = force2_aa1*bonds_rbh(i)*bonds_rbh(i)
                arga2        = force2_aa4*(bonds_rbh(i)-bonds_r0hh)*(bonds_rbh(i)-bonds_r0hh)
                a1           = 1.0e0_wp-exp(-arga1)
                a2           = force2_aa2+force2_aa3*exp(-arga2)
                force1_f1(i) = a1*exp(-a2*(bonds_rch(i)-bonds_r0ch)**2)
                !  and calculate the derivatives wrt rch(i) and rbh(i)
                duma1    = 2.0e0_wp*force2_aa1*bonds_rbh(i)*exp(-arga1)
                duma2    = -2.0e0_wp*force2_aa3*force2_aa4*(bonds_rbh(i)-bonds_r0hh)*exp(-arga2)
                df1dc(i) = -2.0e0_wp*(bonds_rch(i)-bonds_r0ch)*a1*a2*exp(-a2*(bonds_rch(i)-bonds_r0ch)**2)
                df1dh(i) = duma1*exp(-a2*(bonds_rch(i)-bonds_r0ch)**2) -duma2*(bonds_rch(i)-bonds_r0ch)**2*a1* exp(- &
                    & a2*(bonds_rch(i)-bonds_r0ch)**2)
            enddo
            !  derivative of total force constant f(i,j) wrt bond length rch(k)
            !  is given by dfdc(i,j,k)
            force1_dfdc(1,2,1) = dfk0(1,2,1)*force1_f1(1)*force1_f1(2)+force1_fk0(1,2)*df1dc(1)*force1_f1(2)
            force1_dfdc(1,2,2) = dfk0(1,2,2)*force1_f1(1)*force1_f1(2)+force1_fk0(1,2)*force1_f1(1)*df1dc(2)
            force1_dfdc(1,2,3) = dfk0(1,2,3)*force1_f1(1)*force1_f1(2)
            force1_dfdc(1,2,4) = dfk0(1,2,4)*force1_f1(1)*force1_f1(2)
            force1_dfdc(1,3,1) = dfk0(1,3,1)*force1_f1(1)*force1_f1(3)+force1_fk0(1,3)*df1dc(1)*force1_f1(3)
            force1_dfdc(1,3,2) = dfk0(1,3,2)*force1_f1(1)*force1_f1(3)
            force1_dfdc(1,3,3) = dfk0(1,3,3)*force1_f1(1)*force1_f1(3)+force1_fk0(1,3)*force1_f1(1)*df1dc(3)
            force1_dfdc(1,3,4) = dfk0(1,3,4)*force1_f1(1)*force1_f1(3)
            force1_dfdc(1,4,1) = dfk0(1,4,1)*force1_f1(1)*force1_f1(4)+force1_fk0(1,4)*df1dc(1)*force1_f1(4)
            force1_dfdc(1,4,2) = dfk0(1,4,2)*force1_f1(1)*force1_f1(4)
            force1_dfdc(1,4,3) = dfk0(1,4,3)*force1_f1(1)*force1_f1(4)
            force1_dfdc(1,4,4) = dfk0(1,4,4)*force1_f1(1)*force1_f1(4)+force1_fk0(1,4)*force1_f1(1)*df1dc(4)
            force1_dfdc(2,3,1) = dfk0(2,3,1)*force1_f1(2)*force1_f1(3)
            force1_dfdc(2,3,2) = dfk0(2,3,2)*force1_f1(2)*force1_f1(3)+force1_fk0(2,3)*df1dc(2)*force1_f1(3)
            force1_dfdc(2,3,3) = dfk0(2,3,3)*force1_f1(2)*force1_f1(3)+force1_fk0(2,3)*force1_f1(2)*df1dc(3)
            force1_dfdc(2,3,4) = dfk0(2,3,4)*force1_f1(2)*force1_f1(3)
            force1_dfdc(2,4,1) = dfk0(2,4,1)*force1_f1(2)*force1_f1(4)
            force1_dfdc(2,4,2) = dfk0(2,4,2)*force1_f1(2)*force1_f1(4)+force1_fk0(2,4)*df1dc(2)*force1_f1(4)
            force1_dfdc(2,4,3) = dfk0(2,4,3)*force1_f1(2)*force1_f1(4)
            force1_dfdc(2,4,4) = dfk0(2,4,4)*force1_f1(2)*force1_f1(4)+force1_fk0(2,4)*force1_f1(2)*df1dc(4)
            force1_dfdc(3,4,1) = dfk0(3,4,1)*force1_f1(3)*force1_f1(4)
            force1_dfdc(3,4,2) = dfk0(3,4,2)*force1_f1(3)*force1_f1(4)
            force1_dfdc(3,4,3) = dfk0(3,4,3)*force1_f1(3)*force1_f1(4)+force1_fk0(3,4)*df1dc(3)*force1_f1(4)
            force1_dfdc(3,4,4) = dfk0(3,4,4)*force1_f1(3)*force1_f1(4)+force1_fk0(3,4)*force1_f1(3)*df1dc(4)
            !  derivative of total force constant f(i,j) wrt bond length rbh(k)
            !  is given by dfdh(i,j,k)
            !  nb only non-zero derivatives are those from rbh(i) and rbh(j)
            force1_dfdh(1,2,1) = force1_fk0(1,2)*df1dh(1)*force1_f1(2)
            force1_dfdh(1,2,2) = force1_fk0(1,2)*force1_f1(1)*df1dh(2)
            force1_dfdh(1,2,3) = 0.0e0_wp
            force1_dfdh(1,2,4) = 0.0e0_wp
            force1_dfdh(1,3,1) = force1_fk0(1,3)*df1dh(1)*force1_f1(3)
            force1_dfdh(1,3,2) = 0.0e0_wp
            force1_dfdh(1,3,3) = force1_fk0(1,3)*force1_f1(1)*df1dh(3)
            force1_dfdh(1,3,4) = 0.0e0_wp
            force1_dfdh(1,4,1) = force1_fk0(1,4)*df1dh(1)*force1_f1(4)
            force1_dfdh(1,4,2) = 0.0e0_wp
            force1_dfdh(1,4,3) = 0.0e0_wp
            force1_dfdh(1,4,4) = force1_fk0(1,4)*force1_f1(1)*df1dh(4)
            force1_dfdh(2,3,1) = 0.0e0_wp
            force1_dfdh(2,3,2) = force1_fk0(2,3)*df1dh(2)*force1_f1(3)
            force1_dfdh(2,3,3) = force1_fk0(2,3)*force1_f1(2)*df1dh(3)
            force1_dfdh(2,3,4) = 0.0e0_wp
            force1_dfdh(2,4,1) = 0.0e0_wp
            force1_dfdh(2,4,2) = force1_fk0(2,4)*df1dh(2)*force1_f1(4)
            force1_dfdh(2,4,3) = 0.0e0_wp
            force1_dfdh(2,4,4) = force1_fk0(2,4)*force1_f1(2)*df1dh(4)
            force1_dfdh(3,4,1) = 0.0e0_wp
            force1_dfdh(3,4,2) = 0.0e0_wp
            force1_dfdh(3,4,3) = force1_fk0(3,4)*df1dh(3)*force1_f1(4)
            force1_dfdh(3,4,4) = force1_fk0(3,4)*force1_f1(3)*df1dh(4)
        end subroutine ipforce

        !> @brief Evaluate switchf for the CH4Br_VBMM_2002 potential.
        subroutine switchf
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: argmax, args1, args2, args3, argsphi, argstheta
            integer  :: i
            !  calculates switching functions: s3,sphi,stheta
            !  and their derivatives ds3,dsphi,dstheta
            !  nb remember that integration units are:
            !  energy in   1.0d+05 j/mol
            !  time in     1.0d-14 s
            fsw1_a1s = 1.5313681e-7_wp
            fsw1_b1s = -4.6696246e0_wp
            fsw1_a2s = 1.0147402e-7_wp
            fsw1_b2s = -12.363798e0_wp
            !  use double precision criterion:
            !  tanh(19.0d0)=1.0d0
            argmax = 19.0e0_wp
            do i=1, 4
                args1 = fsw1_a1s*(bonds_rch(i)-bonds_r0ch)*(bonds_rch(i)-fsw1_b1s)**8
                if( args1 .lt. argmax )then
                    ip1_s1(i) = 1.0e0_wp-tanh(args1)
                    ip1_ds1(i) = fsw1_a1s*((bonds_rch(i)-fsw1_b1s)**8 +8.0e0_wp*(bonds_rch(i)-bonds_r0ch)* &
                        & (bonds_rch(i)- fsw1_b1s)**7)
                    ip1_ds1(i) = -ip1_ds1(i)/cosh(args1)**2
                else
                    ip1_s1(i)  = 0.0e0_wp
                    ip1_ds1(i) = 0.0e0_wp
                endif
                args2 = fsw1_a2s*(bonds_rch(i)-bonds_r0ch)*(bonds_rch(i)-fsw1_b2s)**6
                if( args2 .lt. argmax )then
                    ip1_s2(i) = 1.0e0_wp-tanh(args2)
                    ip1_ds2(i) = fsw1_a2s*((bonds_rch(i)-fsw1_b2s)**6 +6.0e0_wp*(bonds_rch(i)-bonds_r0ch)* &
                        & (bonds_rch(i)- fsw1_b2s)**5)
                    ip1_ds2(i) = -ip1_ds2(i)/cosh(args2)**2
                else
                    ip1_s2(i)  = 0.0e0_wp
                    ip1_ds2(i) = 0.0e0_wp
                endif
                !  calculate s3 and ds3
                args3 = fsw1_a3s*(bonds_rch(i)-bonds_r0ch)*(bonds_rch(i)-fsw1_b3s)**2
                if( args3 .lt. argmax )then
                    op1_s3(i) = 1.0e0_wp-tanh(args3)
                    op1_ds3(i) = fsw1_a3s*(3.0e0_wp*bonds_rch(i)**2-2.0e0_wp*bonds_rch(i)*(bonds_r0ch+2.0e0_wp* &
                        & fsw1_b3s) + fsw1_b3s*(fsw1_b3s+2.0e0_wp*bonds_r0ch))
                    op1_ds3(i) = -op1_ds3(i)/cosh(args3)**2
                else
                    op1_s3(i)  = 0.0e0_wp
                    op1_ds3(i) = 0.0e0_wp
                endif
                !  calculate sphi and dsphi
                !  condition here is on the bondlength rch(i)
                !  st argsphi is lt approx 19.0d0
                if( bonds_rch(i) .lt. 3.8e0_wp )then
                    argsphi         = sphi1_aphi*(bonds_rch(i)-bonds_r0ch)*exp(sphi1_bphi*(bonds_rch(i)-sphi1_cphi)**3)
                    switch1_sphi(i) = 1.0e0_wp-tanh(argsphi)
                    switch1_dsphi(i) = sphi1_aphi*(1.0e0_wp+3.0e0_wp*sphi1_bphi*(bonds_rch(i)-bonds_r0ch) * &
                        & (bonds_rch(i)- sphi1_cphi)**2)
                    switch1_dsphi(i) = switch1_dsphi(i)*exp(sphi1_bphi*(bonds_rch(i)-sphi1_cphi)**3)
                    switch1_dsphi(i) = -switch1_dsphi(i)/cosh(argsphi)**2
                else
                    switch1_sphi(i)  = 0.0e0_wp
                    switch1_dsphi(i) = 0.0e0_wp
                endif
                !  calculate stheta and dstheta
                if( bonds_rch(i) .lt. 3.8e0_wp )then
                    argstheta  = stheta1_atheta*(bonds_rch(i)-bonds_r0ch)*exp(stheta1_btheta*(bonds_rch(i)- &
                        & stheta1_ctheta)**3)
                    switch1_stheta(i) = 1.0e0_wp-tanh(argstheta)
                    switch1_dstheta(i) = stheta1_atheta*(1.0e0_wp+3.0e0_wp*stheta1_btheta*(bonds_rch(i)-bonds_r0ch) * &
                        & (bonds_rch(i)-stheta1_ctheta)**2)
                    switch1_dstheta(i) = switch1_dstheta(i)*exp(stheta1_btheta*(bonds_rch(i)-stheta1_ctheta)**3)
                    switch1_dstheta(i) = -switch1_dstheta(i)/cosh(argstheta)**2
                else
                    switch1_stheta(i)  = 0.0e0_wp
                    switch1_dstheta(i) = 0.0e0_wp
                endif
            enddo
        end subroutine switchf

        !> @brief Evaluate SETUP for the CH4Br_VBMM_2002 potential.
        subroutine setup( n3tm )
            use chempotpy_kinds, only: wp
            implicit none
            integer :: n3tm
            !   N3TMMN = 3 * NATOMS
            !   NATOMS = the number of atoms represented by this potential function
            !   The variable N3TMMN is the minimum value of N3TM allowed to be
            !   passed by the calling routine for the number of cartesian
            !   coordinates needed to represent the full system represented by this
            !   potential energy surface routine.
            !   N3TM must be greater than or equal to N3TMMN.
            integer, parameter :: n3tmmn = 18
            !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
            if( n3tm .lt. n3tmmn )then
                write (6, '(/,2X,T5,''Warning: N3TM is set equal to '',I3, '' but this potential routine'', /,2X,T14,''requir' // &
                    'es N3TM be greater than or '', ''equal to '',I3,/)') n3tm, n3tmmn
                stop 'SETUP 1'
            endif
        end subroutine setup
    end subroutine evaluate
end module chempotpy_ch4br_vbmm_2002_pes

!> @brief Expose CH4Br_VBMM_2002 through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_ch4br_vbmm_2002_pes, only: evaluate
    use chempotpy_kinds,               only: wp
    implicit none
    integer, parameter              :: natoms  = 6
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
