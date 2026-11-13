!> @brief Implement C2H7_VBMM, a potential energy surface for C2H7.
module chempotpy_c2h7_vbmm_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the C2H7_VBMM potential energy surface for C2H7.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units, only: electronvolt_to_kilocalorie_per_mole
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 9
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = isurf*(isurf+1)/2
        integer, parameter    :: n3atom  = 75
        integer, parameter    :: natom   = 25
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: angles_dtheta0(4,4,4)
        real(wp)              :: angles_theta0(4,4)
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
        real(wp)              :: fsw1_a1s
        real(wp)              :: fsw1_a2s
        real(wp)              :: fsw1_b1s
        real(wp)              :: fsw1_b2s
        real(wp)              :: infocm_cartnu(natom,3)
        integer               :: infocm_icartr
        integer               :: infocm_indexes(natom)
        integer               :: infocm_irctnt
        integer               :: infocm_mder
        integer               :: infocm_msurf
        integer               :: infocm_natoms
        character             :: infocm_ref(5)
        real(wp)              :: ip1_ds1(4)
        real(wp)              :: ip1_ds2(4)
        real(wp)              :: ip1_s1(4)
        real(wp)              :: ip1_s2(4)
        integer               :: ndx_nc(3)
        integer               :: ndx_nh(4,3)
        integer               :: ndx_nhb(3)
        real(wp)              :: op1_ds3(4)
        real(wp)              :: op1_s3(4)
        real(wp)              :: potcm_a1ch
        real(wp)              :: potcm_a3s
        real(wp)              :: potcm_aa1
        real(wp)              :: potcm_aa2
        real(wp)              :: potcm_aa3
        real(wp)              :: potcm_aa4
        real(wp)              :: potcm_acb
        real(wp)              :: potcm_ahh
        real(wp)              :: potcm_ak
        real(wp)              :: potcm_aphi
        real(wp)              :: potcm_atheta
        real(wp)              :: potcm_b1ch
        real(wp)              :: potcm_b3s
        real(wp)              :: potcm_bk
        real(wp)              :: potcm_bphi
        real(wp)              :: potcm_btheta
        real(wp)              :: potcm_c1ch
        real(wp)              :: potcm_cphi
        real(wp)              :: potcm_ctheta
        real(wp)              :: potcm_d1cb
        real(wp)              :: potcm_d1ch
        real(wp)              :: potcm_d1hh
        real(wp)              :: potcm_d3cb
        real(wp)              :: potcm_d3ch
        real(wp)              :: potcm_d3hh
        real(wp)              :: potcm_fch3
        real(wp)              :: potcm_fkinf
        real(wp)              :: potcm_hch3
        integer               :: potcm_nnb
        integer               :: potcm_nnc
        integer               :: potcm_nnh(4)
        real(wp)              :: potcm_r0cb
        real(wp)              :: potcm_r0ch
        real(wp)              :: potcm_r0hh
        real(wp)              :: pt1cm_degsdr(n3atom)
        real(wp)              :: pt1cm_engygs
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: pt4cm_deesdr(n3atom,isurf)
        real(wp)              :: pt4cm_engyes(isurf)
        real(wp)              :: pt5cm_deijdr(n3atom,jsurf)
        real(wp)              :: pt5cm_engyij(jsurf)
        real(wp)              :: qpdot_pl_pdot(150)
        real(wp)              :: qpdot_pl_q(150)
        real(wp)              :: switch1_dsphi(4)
        real(wp)              :: switch1_dstheta(4)
        real(wp)              :: switch1_sphi(4)
        real(wp)              :: switch1_stheta(4)
        real(wp)              :: usricm_anuzero
        real(wp)              :: usricm_cart(natom,3)
        integer               :: usricm_nasurf(isurf+1,isurf+1)
        integer               :: usricm_nder
        integer               :: usricm_nflag(20)
        integer               :: usricm_nulbl(natom)
        real(wp)              :: usrocm_descart(natom,3,isurf)
        real(wp)              :: usrocm_dgscart(natom,3)
        real(wp)              :: usrocm_dijcart(natom,3,jsurf)
        real(wp)              :: usrocm_pengyes(isurf)
        real(wp)              :: usrocm_pengygs
        real(wp)              :: usrocm_pengyij(jsurf)
        real(wp)              :: xangles_dtheta0(4,4,4)
        real(wp)              :: xangles_theta0(4,4)
        real(wp)              :: v
        integer               :: iatom, idir, j, istate, jstate
        !initialize
        angles_dtheta0  = 0.e0_wp; angles_theta0 = 0.e0_wp; bonds_rbh = 0.e0_wp; bonds_rcb = 0.e0_wp; bonds_rch = 0.e0_wp
        coords_tbh      = 0.e0_wp; coords_tcb = 0.e0_wp; coords_tch = 0.e0_wp; delta1_fdelta = 0.e0_wp; delta1_hdelta = 0.e0_wp
        delta2_dfdelta  = 0.e0_wp; delta2_dhdelta = 0.e0_wp; force1_dfdc = 0.e0_wp; force1_dfdh = 0.e0_wp; force1_f1 = 0.e0_wp
        force1_fk0      = 0.e0_wp; fsw1_a1s = 0.e0_wp; fsw1_a2s = 0.e0_wp; fsw1_b1s = 0.e0_wp; fsw1_b2s = 0.e0_wp; infocm_cartnu = 0.e0_wp
        infocm_icartr   = 0; infocm_indexes = 0; infocm_irctnt = 0; infocm_mder = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref      = ' '
        ip1_ds1         = 0.e0_wp
        ip1_ds2         = 0.e0_wp
        ip1_s1          = 0.e0_wp
        ip1_s2          = 0.e0_wp
        ndx_nc          = 0
        ndx_nh          = 0
        ndx_nhb         = 0
        op1_ds3         = 0.e0_wp
        op1_s3          = 0.e0_wp
        potcm_a1ch      = 0.e0_wp
        potcm_a3s       = 0.e0_wp
        potcm_aa1       = 0.e0_wp
        potcm_aa2       = 0.e0_wp
        potcm_aa3       = 0.e0_wp
        potcm_aa4       = 0.e0_wp
        potcm_acb       = 0.e0_wp
        potcm_ahh       = 0.e0_wp
        potcm_ak        = 0.e0_wp
        potcm_aphi      = 0.e0_wp
        potcm_atheta    = 0.e0_wp
        potcm_b1ch      = 0.e0_wp
        potcm_b3s       = 0.e0_wp
        potcm_bk        = 0.e0_wp
        potcm_bphi      = 0.e0_wp
        potcm_btheta    = 0.e0_wp
        potcm_c1ch      = 0.e0_wp
        potcm_cphi      = 0.e0_wp
        potcm_ctheta    = 0.e0_wp
        potcm_d1cb      = 0.e0_wp
        potcm_d1ch      = 0.e0_wp
        potcm_d1hh      = 0.e0_wp
        potcm_d3cb      = 0.e0_wp
        potcm_d3ch      = 0.e0_wp
        potcm_d3hh      = 0.e0_wp
        potcm_fch3      = 0.e0_wp
        potcm_fkinf     = 0.e0_wp
        potcm_hch3      = 0.e0_wp
        potcm_nnb       = 0
        potcm_nnc       = 0
        potcm_nnh       = 0
        potcm_r0cb      = 0.e0_wp
        potcm_r0ch      = 0.e0_wp
        potcm_r0hh      = 0.e0_wp
        pt1cm_degsdr    = 0.e0_wp
        pt1cm_engygs    = 0.e0_wp
        pt1cm_r         = 0.e0_wp
        pt3cm_ezero     = 0.e0_wp
        pt4cm_deesdr    = 0.e0_wp
        pt4cm_engyes    = 0.e0_wp
        pt5cm_deijdr    = 0.e0_wp
        pt5cm_engyij    = 0.e0_wp
        qpdot_pl_pdot   = 0.e0_wp
        qpdot_pl_q      = 0.e0_wp
        switch1_dsphi   = 0.e0_wp
        switch1_dstheta = 0.e0_wp
        switch1_sphi    = 0.e0_wp
        switch1_stheta  = 0.e0_wp
        usricm_anuzero  = 0.e0_wp
        usricm_cart     = 0.e0_wp
        usricm_nasurf   = 0
        usricm_nder     = 0
        usricm_nflag    = 0
        usricm_nulbl    = 0
        usrocm_descart  = 0.e0_wp
        usrocm_dgscart  = 0.e0_wp
        usrocm_dijcart  = 0.e0_wp
        usrocm_pengyes  = 0.e0_wp
        usrocm_pengygs  = 0.e0_wp
        usrocm_pengyij  = 0.e0_wp
        xangles_dtheta0 = 0.e0_wp
        xangles_theta0  = 0.e0_wp
        usricm_nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0], shape(usricm_nasurf))
        usricm_nder    = 1
        usricm_nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero = 0.0e0_wp
        infocm_icartr  = 1
        infocm_msurf   = 0
        infocm_mder    = 1
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 6
        potcm_nnc      = 2
        potcm_nnb      = 6
        potcm_nnh      = [3, 4, 5, 1]
        potcm_r0ch     = 1.09397e0_wp
        potcm_d1ch     = 112.230e0_wp
        potcm_d3ch     = 38.834e0_wp
        potcm_a1ch     = 1.71300e0_wp
        potcm_b1ch     = 0.13500e0_wp
        potcm_c1ch     = 6.61404e0_wp
        potcm_r0hh     = 0.74191e0_wp
        potcm_d1hh     = 109.458e0_wp
        potcm_d3hh     = 39.664e0_wp
        potcm_ahh      = 1.9457e0_wp
        potcm_r0cb     = 1.09397e0_wp
        potcm_d1cb     = 26.409e0_wp
        potcm_d3cb     = 20.063e0_wp
        potcm_acb      = 1.8530285e0_wp
        potcm_a3s      = 0.1419147e0_wp
        potcm_b3s      = -0.3068450e0_wp
        potcm_aphi     = 0.5287903e0_wp
        potcm_bphi     = 0.4006638e0_wp
        potcm_cphi     = 1.9209937e0_wp
        potcm_atheta   = 0.9078714e0_wp
        potcm_btheta   = 0.3548859e0_wp
        potcm_ctheta   = 1.8915497e0_wp
        potcm_fch3     = 0.0957500e0_wp
        potcm_hch3     = 0.1915000e0_wp
        potcm_fkinf    = 0.4077000e0_wp
        potcm_ak       = 0.1260000e0_wp
        potcm_bk       = 10.7132e0_wp
        potcm_aa1      = 3.213952e0_wp
        potcm_aa2      = 1.599963e0_wp
        potcm_aa3      = 2.165953e0_wp
        potcm_aa4      = 11.569977e0_wp
        p              = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        if( igrad == 0 )then
            call c2h7pot(x, v)
            v = v/electronvolt_to_kilocalorie_per_mole
            do istate=1, nstates
                p(istate) = v
            enddo
        else
            write(*, *) 'Only energy is available'
        endif

        contains

        !> @brief Evaluate c2h7poth for the C2H7_VBMM potential.
        subroutine c2h7poth( x, eng_kcal, grad, xnorm, hess )
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter  :: npmm = 23
            integer, parameter  :: npjg = 34
            real(wp)            :: x(9,3)
            real(wp)            :: eng_kcal
            real(wp)            :: grad(9,3)
            real(wp)            :: xnorm(9)
            real(wp)            :: hess(27,27)
            integer             :: i, j
            integer             :: i1, j1, k1
            integer             :: i2, j2, k2
            real(wp)            :: z
            real(wp)            :: x1(9,3)
            real(wp)            :: x2(9,3)
            real(wp)            :: x3(9,3)
            real(wp)            :: x4(9,3)
            real(wp)            :: v1
            real(wp)            :: v2
            real(wp)            :: v3
            real(wp)            :: v4
            real(wp), parameter :: delta = 1.0e-4_wp
            ! get energy
            call c2h7pot(x, eng_kcal)
            ! get numerical gradients
            do i=1, 9
                do j=1, 3
                    x1(i,j) = x(i,j)
                    x2(i,j) = x(i,j)
                    x3(i,j) = x(i,j)
                    x4(i,j) = x(i,j)
                enddo
            enddo
            do i=1, 9
                xnorm(i) = 0.0e0_wp
                do j=1, 3
                    x1(i,j) = x(i,j)-delta
                    x2(i,j) = x(i,j)+delta
                    call c2h7pot(x1, v1)
                    call c2h7pot(x2, v2)
                    grad(i,j) = (v2-v1)*0.50e0_wp/delta
                    x1(i,j)   = x(i,j)
                    x2(i,j)   = x(i,j)
                    xnorm(i)  = xnorm(i)+(grad(i,j)**2)
                enddo
                xnorm(i) = sqrt(xnorm(i))
            enddo
            ! get numerical hessian
            k1 = 0
            do i2=1, 9
                do j2=1, 3
                    k2 = k2+1
                    k1 = 0
                    do i1=1, 9
                        do j1=1, 3
                            k1 = k1+1
                            if( k1 > k2 )then
                                x1(i1,j1) = x(i1,j1)+delta
                                x1(i2,j2) = x(i2,j2)+delta
                                x2(i1,j1) = x(i1,j1)-delta
                                x2(i2,j2) = x(i2,j2)+delta
                                x3(i1,j1) = x(i1,j1)+delta
                                x3(i2,j2) = x(i2,j2)-delta
                                x4(i1,j1) = x(i1,j1)-delta
                                x4(i2,j2) = x(i2,j2)-delta
                                call c2h7pot(x1, v1)
                                call c2h7pot(x2, v2)
                                call c2h7pot(x3, v3)
                                call c2h7pot(x4, v4)
                                z           = 0.250e0_wp*(v1-v2-v3+v4)
                                hess(k1,k2) = z/(delta*delta)
                                hess(k2,k1) = hess(k1,k2)
                                x1(i1,j1)   = x(i1,j1)
                                x1(i2,j2)   = x(i2,j2)
                                x2(i1,j1)   = x(i1,j1)
                                x2(i2,j2)   = x(i2,j2)
                                x3(i1,j1)   = x(i1,j1)
                                x3(i2,j2)   = x(i2,j2)
                                x4(i1,j1)   = x(i1,j1)
                                x4(i2,j2)   = x(i2,j2)
                            elseif( k1 == k2 )then
                                x1(i1,j1) = x(i1,j1)+delta
                                x2(i2,j2) = x(i2,j2)-delta
                                call c2h7pot(x1, v1)
                                call c2h7pot(x2, v2)
                                z           = (v1-eng_kcal-eng_kcal+v2)
                                hess(k1,k1) = z/(delta*delta)
                                x1(i1,j1)   = x(i1,j1)
                                x2(i2,j2)   = x(i2,j2)
                            endif
                        enddo
                    enddo
                enddo
            enddo
        end subroutine c2h7poth

        !> @brief Evaluate c2h7potg for the C2H7_VBMM potential.
        subroutine c2h7potg( x, eng_kcal, grad, xnorm )
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter  :: npmm = 23
            integer, parameter  :: npjg = 34
            real(wp)            :: x(9,3)
            real(wp)            :: eng_kcal
            real(wp)            :: grad(9,3)
            real(wp)            :: xnorm(9)
            integer             :: i, j
            real(wp)            :: x1(9,3)
            real(wp)            :: x2(9,3)
            real(wp)            :: v1
            real(wp)            :: v2
            real(wp), parameter :: delta = 1.0e-4_wp
            ! get energy
            call c2h7pot(x, eng_kcal)
            ! get numerical gradients
            do i=1, 9
                do j=1, 3
                    x1(i,j) = x(i,j)
                    x2(i,j) = x(i,j)
                enddo
            enddo
            do i=1, 9
                xnorm(i) = 0.0e0_wp
                do j=1, 3
                    x1(i,j) = x(i,j)-delta
                    x2(i,j) = x(i,j)+delta
                    call c2h7pot(x1, v1)
                    call c2h7pot(x2, v2)
                    grad(i,j) = (v2-v1)*0.50e0_wp/delta
                    x1(i,j)   = x(i,j)
                    x2(i,j)   = x(i,j)
                    xnorm(i)  = xnorm(i)+(grad(i,j)**2)
                enddo
                xnorm(i) = sqrt(xnorm(i))
            enddo
        end subroutine c2h7potg

        !> @brief Evaluate c2h7pot for the C2H7_VBMM potential.
        subroutine c2h7pot( x, eng_kcal )
            ! Coordinates in Ang
            ! Energy in kcal
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: npmm = 25
            integer, parameter :: npjg = 34
            real(wp)           :: x(9,3)
            real(wp)           :: c1(3)
            real(wp)           :: h2(3)
            real(wp)           :: h3(3)
            real(wp)           :: h4(3)
            real(wp)           :: c5(3)
            real(wp)           :: h6(3)
            real(wp)           :: h7(3)
            real(wp)           :: h8(3)
            real(wp)           :: h9(3)
            real(wp)           :: eng_kcal
            real(wp)           :: vmm3
            real(wp)           :: vjg
            real(wp)           :: etotalmm
            real(wp)           :: paramm(npmm)
            real(wp)           :: parajg(npjg)
            ! average bond length
            real(wp) :: avgch
            real(wp) :: avgbh
            real(wp) :: rxncrd
            real(wp) :: r1, r2, r3
            real(wp) :: b1, b2, b3
            ! scale=C-H/C-C
            !%%%  parameter(scale=0.7153679930d0)
            real(wp), parameter :: scale = 0.7171830610e0_wp
            real(wp)            :: h5prime(3)
            integer             :: i
            logical, parameter  :: lmm3 = .false.
            ! Convert geometry from:
            !     H2 H6
            !     |  |
            !  H4-C1-C5-H8.....H9
            !     |  |
            !     H3 H7
            ! to:
            !     H6 H3
            !     |  |
            !  H8-C5-C1-H2.....H9
            !     |  |
            !     H7 H4
            do i=1, 3
                c1(i)      = x(5,i)
                h2(i)      = x(8,i)
                h3(i)      = x(6,i)
                h4(i)      = x(7,i)
                c5(i)      = x(1,i)
                h6(i)      = x(2,i)
                h7(i)      = x(3,i)
                h8(i)      = x(4,i)
                h9(i)      = x(9,i)
                h5prime(i) = (scale*c5(i)) + ((1.0e0_wp-scale)*c1(i))
            enddo
            ! Calculate avg C-H bond length for 3 QM H
            r1 = 0.0e0_wp
            r2 = 0.0e0_wp
            r3 = 0.0e0_wp
            b1 = 0.0e0_wp
            b2 = 0.0e0_wp
            b3 = 0.0e0_wp
            do i=1, 3
                r1 = r1+( (c1(i)-h2(i))**2 )
                r2 = r2+( (c1(i)-h3(i))**2 )
                r3 = r3+( (c1(i)-h4(i))**2 )
                b1 = b1+( (h9(i)-h2(i))**2 )
                b2 = b2+( (h9(i)-h3(i))**2 )
                b3 = b3+( (h9(i)-h4(i))**2 )
            enddo
            r1     = sqrt(r1)
            r2     = sqrt(r2)
            r3     = sqrt(r3)
            b1     = sqrt(b1)
            b2     = sqrt(b2)
            b3     = sqrt(b3)
            avgch  = (r1+r2+r3)/3.0e0_wp
            avgbh  = (1.0e0_wp/b1)+(1.0e0_wp/b2)+(1.0e0_wp/b3)
            rxncrd = avgch*avgbh
            ! Get energy
            if( lmm3 )then
                call pmm3(npmm, paramm)
                call ethane(npmm, paramm, rxncrd, &
                    & c1, h2, h3, h4, c5, h6, h7, h8, vmm3, etotalmm)
                eng_kcal = etotalmm
            else
                call pmm3_mod(npmm, paramm, avgch, rxncrd)
                call pjg_mod(npjg, parajg, avgch)
                call ethane(npmm, paramm, rxncrd, &
                    & c1, h2, h3, h4, c5, h6, h7, h8, vmm3, etotalmm)
                call xch5pot(npjg, parajg, rxncrd, &
                    & c1, h2, h3, h4, h5prime, h9, vjg)
                eng_kcal = vjg+vmm3
            endif
        end subroutine c2h7pot

        !> @brief Evaluate pmm3_mod for the C2H7_VBMM potential.
        subroutine pmm3_mod( np, p, ravg, rxncrd )
            use chempotpy_kinds, only: wp
            implicit none
            ! Output Variables
            integer  :: np
            real(wp) :: p(np)
            real(wp) :: ravg
            real(wp) :: rxncrd
            ! Local Variables
            real(wp) :: r0ch, r0cc
            real(wp) :: t0hch, t0cch
            real(wp) :: t0hcrh, t0ccrh
            real(wp) :: r0crh, r0ccr
            !     parameter(r0ch=1.1120d0,r0cc=1.52470d0)
            !##   parameter(r0ch=1.09120d0,r0cc=1.51950d0)
            real(wp), parameter :: ksch = 4.740e0_wp
            real(wp), parameter :: kscc = 4.490e0_wp
            !%%%  parameter(t0hch=107.80d0,t0cch=110.70d0)
            !%%%  parameter(t0hcRh=107.60d0,t0ccRh=109.310d0)
            real(wp), parameter :: kthch   = 0.550e0_wp
            real(wp), parameter :: ktcch   = 0.590e0_wp
            real(wp), parameter :: ksbhch  = 0.0e0_wp
            real(wp), parameter :: ksbcch  = 0.080e0_wp
            real(wp), parameter :: kbbhch  = 0.0e0_wp
            real(wp), parameter :: kbbcch  = 0.30e0_wp**2
            real(wp), parameter :: v1hcch  = 0.0e0_wp
            real(wp), parameter :: v2hcch  = 0.0e0_wp
            real(wp), parameter :: v3hcch  = 0.2380e0_wp
            real(wp), parameter :: ksthcch = 0.0590e0_wp
            real(wp), parameter :: espc    = 0.0270e0_wp
            real(wp), parameter :: rvc     = 2.040e0_wp
            real(wp), parameter :: esph    = 0.020e0_wp
            real(wp), parameter :: rvh     = 1.620e0_wp
            real(wp), parameter :: scalech = 0.9230e0_wp
            !%%%  parameter(r0cRh=1.1120d0,r0ccR=1.52470d0)
            ! Switching
            real(wp) :: ravg0
            real(wp) :: cr, cp, cs
            ! Modified parameters
            integer  :: i, j
            integer  :: ipt
            real(wp) :: x
            real(wp) :: t1, t2, t3
            real(wp) :: t4, t5, t6
            real(wp) :: z1, z2, z3
            real(wp) :: z4, z5, z6
            real(wp) :: ar0ch, ar0cc
            real(wp) :: asch, ascc
            real(wp) :: av3hcch
            real(wp) :: aksthcch
            real(wp) :: akthch
            real(wp) :: aktcch
            real(wp) :: aesph
            real(wp) :: arvh
            real(wp) :: at0hch
            real(wp) :: ar0crh
            real(wp) :: ar0ccr
            real(wp) :: at0ccrh
            !%%%  parameter(r0ch=1.1120d0,r0cc=1.52470d0)
            !     parameter(r0ch=1.092510d0,r0cc=1.52720d0)
            !%%% for ethane
            ipt = 4
            if( ipt .eq. 1 )then
                !       ETHANE
                r0ch   = 1.0850e0_wp
                r0cc   = 1.510e0_wp
                t0hch  = 107.52790e0_wp
                t0cch  = 111.34990e0_wp
                r0crh  = r0ch
                r0ccr  = 1.47110e0_wp
                t0hcrh = 120.0e0_wp
                t0ccrh = 120.0e0_wp
            elseif( ipt .eq. 2 )then
                !       TS
                r0ch   = 1.0850e0_wp
                r0cc   = 1.510e0_wp
                t0hch  = 107.52790e0_wp
                t0cch  = 111.34990e0_wp
                r0crh  = r0ch
                r0ccr  = 1.47110e0_wp
                t0hcrh = 120.0e0_wp
                t0ccrh = 120.0e0_wp
            elseif( ipt .eq. 3 )then
                !       ETHYL
                r0ch   = 1.0850e0_wp
                r0cc   = 1.510e0_wp
                t0hch  = 107.52790e0_wp
                t0cch  = 111.34990e0_wp
                r0crh  = r0ch
                r0ccr  = 1.47110e0_wp
                t0hcrh = 120.0e0_wp
                t0ccrh = 120.0e0_wp
            elseif( ipt .eq. 4 )then
                ! geometry
                r0ch   = 1.0850e0_wp
                r0cc   = 1.50e0_wp
                t0hch  = 107.52790e0_wp
                t0cch  = 111.34990e0_wp
                r0crh  = r0ch
                r0ccr  = 1.47110e0_wp
                t0hcrh = 120.0e0_wp
                t0ccrh = 120.0e0_wp
                !%%%
                t0ccrh = 117.2989954203791430e0_wp
                t0ccrh = 121.0e0_wp
                t0cch  = t0ccrh-7.10e0_wp
                t1     = 6.0e0_wp
                t2     = 0.90e0_wp
                t3     = 0.50e0_wp*(1.0e0_wp+tanh(t1*(rxncrd-t2)))
                t4     = 1.0e0_wp-t3
                t0cch  = (t0cch*t3)+(111.34990e0_wp*t4)
                t0ccrh = (t0ccrh*t3)+(111.34990e0_wp*t4)
                !       t1=6.0d0
                !       t2=20.0d0
                !       t3=0.50d0*(1.0d0+dtanh(t1*(rxncrd-t2)))
                !       t4=1.0d0-t3
                !       t0cch=(t0cch*t4)+(120.0d0*t3)
                !       t0ccRh=(t0ccRh*t4)+(120.0d0*t3)
                ! freq and energy
                cs     = 10.0e0_wp
                ravg0  = 1.17320e0_wp
                cp     = 0.50e0_wp*(1.0e0_wp+(tanh(cs*(ravg-ravg0))))
                cr     = 1.0e0_wp-cp
                ascc   = 0.70e0_wp
                aktcch = 0.90e0_wp
                akthch = 0.90e0_wp
                !       z1=10.0d0
                !       z2=2.144106410736830d0
                !       z3=dexp(-z1*((rxncrd-z2)**2))
                !       Ascc=1.0d0-(0.30d0*z3)
                !       Aktcch=1.0d0-(0.10d0*z3)
                !       Akthch=1.0d0-(0.10d0*z3)
                t1     = 6.0e0_wp
                t2     = 0.950e0_wp
                t3     = 0.50e0_wp*(1.0e0_wp+tanh(t1*(rxncrd-t2)))
                t4     = (1.0e0_wp*(1.0e0_wp-t3))+(0.70e0_wp*t3)
                t5     = (1.0e0_wp*(1.0e0_wp-t3))+(0.90e0_wp*t3)
                t6     = (1.0e0_wp*(1.0e0_wp-t3))+(0.90e0_wp*t3)
                z1     = 6.0e0_wp
                z2     = 10.0e0_wp
                z3     = 0.50e0_wp*(1.0e0_wp+tanh(z1*(rxncrd-z2)))
                z4     = (t4*(1.0e0_wp-z3))+(1.0e0_wp*z3)
                z5     = (t5*(1.0e0_wp-z3))+(1.0e0_wp*z3)
                z6     = (t6*(1.0e0_wp-z3))+(1.0e0_wp*z3)
                ascc   = z4
                aktcch = z5
                akthch = z6
            endif
            !     Ar0ch=0.980d0
            !     Ar0cc=0.9950d0
            !     Av3hcch=1.00d0
            !     Aksthcch=1.0d0
            !     Aesph=1.0d0
            !     Arvh=1.0d0
            !     At0hch=1.00d0
            !     Akthch=1.00d0
            !     Aktcch=1.00d0
            !     Ar0cRh=1.0d0
            !     Ar0ccR=0.980d0
            !     At0ccRh=1.100d0
            !     Asch=1.10d0
            !%%%
            p(1)  = r0ch
            p(2)  = r0cc
            p(3)  = ksch
            p(4)  = kscc * ascc
            p(5)  = t0hch
            p(6)  = t0cch
            p(7)  = t0hcrh
            p(8)  = t0ccrh
            p(9)  = kthch * akthch
            p(10) = ktcch * aktcch
            p(11) = ksbhch
            p(12) = ksbcch
            p(13) = kbbhch
            p(14) = kbbcch
            p(15) = v1hcch
            p(16) = v2hcch
            p(17) = v3hcch
            p(18) = ksthcch
            p(19) = espc
            p(20) = rvc
            p(21) = esph
            p(22) = rvh
            p(23) = scalech
            p(24) = r0crh
            p(25) = r0ccr
        end subroutine pmm3_mod

        !> @brief Evaluate pjg_mod for the C2H7_VBMM potential.
        subroutine pjg_mod( np, para, ravg )
            use chempotpy_kinds, only: wp
            implicit none
            ! I/O variables
            integer  :: np
            real(wp) :: para(np)
            real(wp) :: ravg
            ! Parameter list
            real(wp)            :: r0cb
            real(wp)            :: r0hh
            real(wp)            :: alpch
            real(wp)            :: r0ch
            real(wp)            :: tau
            real(wp), parameter :: fact1 = 0.041840e0_wp
            real(wp), parameter :: fact2 = 6.022045e0_wp
            ! streching parameters
            real(wp), parameter :: a1ch  = 1.7130e0_wp
            real(wp), parameter :: b1ch  = 0.1350e0_wp
            real(wp), parameter :: c1ch  = 6.6140e0_wp
            real(wp), parameter :: alpcb = 1.8530285e0_wp
            real(wp), parameter :: d1cb  = 26.4090e0_wp*fact1
            real(wp), parameter :: d3cb  = 20.0630e0_wp*fact1
            !     parameter (r0cb=1.093970d0)
            !%%%  parameter (r0cb=1.092510d0)
            !%%%  parameter (r0ch=1.092510d0)
            !##   parameter (r0cb=1.09030d0)
            !##   parameter (r0ch=1.09070d0)
            real(wp), parameter :: alphh = 1.94570e0_wp
            real(wp), parameter :: d1hh  = 109.4580e0_wp*fact1
            real(wp), parameter :: d3hh  = 39.6640e0_wp*fact1
            !%%%  parameter (r0hh=0.741910d0)
            !     parameter (alphch=CALC_AT_RUNTIME)
            real(wp), parameter :: d1ch = 112.230e0_wp*fact1
            real(wp), parameter :: d3ch = 38.8340e0_wp*fact1
            !     parameter (r0ch=1.093970d0)
            !##   parameter (r0ch=1.092510d0)
            ! out-of-plane bending parameters
            !     parameter (tau=CALC_AT_RUNTIME)
            real(wp), parameter :: aphi   = 0.52879030e0_wp
            real(wp), parameter :: bphi   = 0.40066380e0_wp
            real(wp), parameter :: cphi   = 1.92099370e0_wp
            real(wp), parameter :: atheta = 0.90787140e0_wp
            real(wp), parameter :: btheta = 0.35488590e0_wp
            real(wp), parameter :: ctheta = 1.89154970e0_wp
            real(wp), parameter :: a3s    = 0.14191470e0_wp
            real(wp), parameter :: b3s    = -0.3068450e0_wp
            real(wp), parameter :: fd0    = 0.0957500e0_wp*fact2
            real(wp), parameter :: hd0    = 0.1915000e0_wp*fact2
            ! in-plane bending parameters
            real(wp), parameter :: a1s  = 1.53136810e-7_wp
            real(wp), parameter :: a2s  = 1.01474020e-7_wp
            real(wp), parameter :: b1s  = -4.6696246e0_wp
            real(wp), parameter :: b2s  = -12.363798e0_wp
            real(wp), parameter :: ak   = 0.1260e0_wp*fact2
            real(wp), parameter :: kch3 = 0.40770e0_wp*fact2
            real(wp), parameter :: aa1  = 3.2139520e0_wp
            real(wp), parameter :: aa2  = 1.5999630e0_wp
            real(wp), parameter :: aa3  = 2.1659530e0_wp
            real(wp), parameter :: aa4  = 11.5699770e0_wp
            ! Switching
            real(wp) :: ravg0
            real(wp) :: cr, cp, cs
            integer  :: i
            integer  :: ipt
            logical  :: labinitio
            real(wp) :: ad1hh
            real(wp) :: ad1ch
            real(wp) :: aalpcb
            !%%% QM parameters
            ipt    = 4
            ad1hh  = 1.0e0_wp
            ad1ch  = 1.0e0_wp
            aalpcb = 1.0e0_wp
            if( ipt .eq. 1 )then
                !       ETHANE
                r0ch = 1.0850e0_wp
                r0cb = 1.400e0_wp
                r0hh = 0.680e0_wp
            elseif( ipt .eq. 2 )then
                !       TS
                r0ch = 1.080e0_wp
                r0cb = 1.400e0_wp
                r0hh = 0.660e0_wp
            elseif( ipt .eq. 3 )then
                !       ETHYL
                r0ch = 1.0750e0_wp
                r0cb = 1.400e0_wp
                r0hh = 0.73590e0_wp
            elseif( ipt .eq. 4 )then
                ! for geometry
                cs    = 10.0e0_wp
                ravg0 = 1.17320e0_wp
                cp    = 0.50e0_wp*(1.0e0_wp+(tanh(cs*(ravg-ravg0))))
                cr    = 1.0e0_wp-cp
                r0ch  = (1.0850e0_wp*cr)+(1.0750e0_wp*cp)
                r0cb  = 1.400e0_wp
                !       r0hh=0.73590d0
                ! for rxn energy and freq
                cs    = 10.0e0_wp
                ravg0 = 1.7320e0_wp
                cp    = 0.50e0_wp*(1.0e0_wp+(tanh(cs*(ravg-ravg0))))
                cr    = 1.0e0_wp-cp
                ad1hh = (1.0e0_wp*cr)+(1.027420e0_wp*cp)
                ad1hh = (1.0e0_wp*cr)+(1.027500e0_wp*cp)
                r0hh  = (0.720e0_wp*cr)+(0.73590e0_wp*cp)
                r0hh  = (0.70e0_wp*cr)+(0.73590e0_wp*cp)
            endif
            ! Str
            para(1)  = a1ch
            para(2)  = b1ch
            para(3)  = c1ch
            para(4)  = alpcb
            para(5)  = d1cb
            para(6)  = d3cb
            para(7)  = r0cb
            para(8)  = alphh
            para(9)  = d1hh * ad1hh
            para(10) = d3hh
            para(11) = r0hh
            para(12) = d1ch
            para(13) = d3ch
            para(14) = r0ch
            ! out-of-plane
            para(15) = aphi
            para(16) = bphi
            para(17) = cphi
            para(18) = atheta
            para(19) = btheta
            para(20) = ctheta
            para(21) = a3s
            para(22) = b3s
            para(23) = fd0
            para(24) = hd0
            ! in-plane bending parameters
            para(25) = a1s
            para(26) = a2s
            para(27) = b1s
            para(28) = b2s
            para(29) = ak
            para(30) = kch3
            para(31) = aa1
            para(32) = aa2
            para(33) = aa3
            para(34) = aa4
        end subroutine pjg_mod

        !> @brief Evaluate pmm3 for the C2H7_VBMM potential.
        subroutine pmm3( np, p )
            use chempotpy_kinds, only: wp
            implicit none
            ! Output Variables
            integer  :: np
            real(wp) :: p(np)
            ! Local Variables
            real(wp), parameter :: r0ch    = 1.1120e0_wp
            real(wp), parameter :: r0cc    = 1.52470e0_wp
            real(wp), parameter :: ksch    = 4.740e0_wp
            real(wp), parameter :: kscc    = 4.490e0_wp
            real(wp), parameter :: t0hch   = 107.80e0_wp
            real(wp), parameter :: t0cch   = 110.70e0_wp
            real(wp), parameter :: t0hcrh  = 107.60e0_wp
            real(wp), parameter :: t0ccrh  = 109.310e0_wp
            real(wp), parameter :: kthch   = 0.550e0_wp
            real(wp), parameter :: ktcch   = 0.590e0_wp
            real(wp), parameter :: ksbhch  = 0.0e0_wp
            real(wp), parameter :: ksbcch  = 0.080e0_wp
            real(wp), parameter :: kbbhch  = 0.0e0_wp
            real(wp), parameter :: kbbcch  = 0.30e0_wp**2
            real(wp), parameter :: v1hcch  = 0.0e0_wp
            real(wp), parameter :: v2hcch  = 0.0e0_wp
            real(wp), parameter :: v3hcch  = 0.2380e0_wp
            real(wp), parameter :: ksthcch = 0.0590e0_wp
            real(wp), parameter :: espc    = 0.0270e0_wp
            real(wp), parameter :: rvc     = 2.040e0_wp
            real(wp), parameter :: esph    = 0.020e0_wp
            real(wp), parameter :: rvh     = 1.620e0_wp
            real(wp), parameter :: scalech = 0.9230e0_wp
            real(wp), parameter :: r0crh   = 1.1120e0_wp
            real(wp), parameter :: r0ccr   = 1.52470e0_wp
            p(1)  = r0ch
            p(2)  = r0cc
            p(3)  = ksch
            p(4)  = kscc
            p(5)  = t0hch
            p(6)  = t0cch
            p(7)  = t0hcrh
            p(8)  = t0ccrh
            p(9)  = kthch
            p(10) = ktcch
            p(11) = ksbhch
            p(12) = ksbcch
            p(13) = kbbhch
            p(14) = kbbcch
            p(15) = v1hcch
            p(16) = v2hcch
            p(17) = v3hcch
            p(18) = ksthcch
            p(19) = espc
            p(20) = rvc
            p(21) = esph
            p(22) = rvh
            p(23) = scalech
            p(24) = r0crh
            p(25) = r0ccr
        end subroutine pmm3

        !> @brief Evaluate pjg for the C2H7_VBMM potential.
        subroutine pjg( np, para )
            use chempotpy_kinds, only: wp
            implicit none
            ! I/O variables
            integer  :: np
            real(wp) :: para(np)
            ! Parameter list
            real(wp)            :: alpch
            real(wp)            :: tau
            real(wp), parameter :: fact1 = 0.041840e0_wp
            real(wp), parameter :: fact2 = 6.022045e0_wp
            ! streching parameters
            real(wp), parameter :: a1ch  = 1.7130e0_wp
            real(wp), parameter :: b1ch  = 0.1350e0_wp
            real(wp), parameter :: c1ch  = 6.6140e0_wp
            real(wp), parameter :: alpcb = 1.8530285e0_wp
            real(wp), parameter :: d1cb  = 26.4090e0_wp*fact1
            real(wp), parameter :: d3cb  = 20.0630e0_wp*fact1
            real(wp), parameter :: r0cb  = 1.093970e0_wp
            real(wp), parameter :: alphh = 1.94570e0_wp
            real(wp), parameter :: d1hh  = 109.4580e0_wp*fact1
            real(wp), parameter :: d3hh  = 39.6640e0_wp*fact1
            real(wp), parameter :: r0hh  = 0.741910e0_wp
            !     parameter (alphch=CALC_AT_RUNTIME)
            real(wp), parameter :: d1ch = 112.230e0_wp*fact1
            real(wp), parameter :: d3ch = 38.8340e0_wp*fact1
            real(wp), parameter :: r0ch = 1.093970e0_wp
            ! out-of-plane bending parameters
            !     parameter (tau=CALC_AT_RUNTIME)
            real(wp), parameter :: aphi   = 0.52879030e0_wp
            real(wp), parameter :: bphi   = 0.40066380e0_wp
            real(wp), parameter :: cphi   = 1.92099370e0_wp
            real(wp), parameter :: atheta = 0.90787140e0_wp
            real(wp), parameter :: btheta = 0.35488590e0_wp
            real(wp), parameter :: ctheta = 1.89154970e0_wp
            real(wp), parameter :: a3s    = 0.14191470e0_wp
            real(wp), parameter :: b3s    = -0.3068450e0_wp
            real(wp), parameter :: fd0    = 0.0957500e0_wp*fact2
            real(wp), parameter :: hd0    = 0.1915000e0_wp*fact2
            ! in-plane bending parameters
            real(wp), parameter :: a1s  = 1.53136810e-7_wp
            real(wp), parameter :: a2s  = 1.01474020e-7_wp
            real(wp), parameter :: b1s  = -4.6696246e0_wp
            real(wp), parameter :: b2s  = -12.363798e0_wp
            real(wp), parameter :: ak   = 0.1260e0_wp*fact2
            real(wp), parameter :: kch3 = 0.40770e0_wp*fact2
            real(wp), parameter :: aa1  = 3.2139520e0_wp
            real(wp), parameter :: aa2  = 1.5999630e0_wp
            real(wp), parameter :: aa3  = 2.1659530e0_wp
            real(wp), parameter :: aa4  = 11.5699770e0_wp
            ! Str
            para(1)  = a1ch
            para(2)  = b1ch
            para(3)  = c1ch
            para(4)  = alpcb
            para(5)  = d1cb
            para(6)  = d3cb
            para(7)  = r0cb
            para(8)  = alphh
            para(9)  = d1hh
            para(10) = d3hh
            para(11) = r0hh
            para(12) = d1ch
            para(13) = d3ch
            para(14) = r0ch
            ! out-of-plane
            para(15) = aphi
            para(16) = bphi
            para(17) = cphi
            para(18) = atheta
            para(19) = btheta
            para(20) = ctheta
            para(21) = a3s
            para(22) = b3s
            para(23) = fd0
            para(24) = hd0
            ! in-plane bending parameters
            para(25) = a1s
            para(26) = a2s
            para(27) = b1s
            para(28) = b2s
            para(29) = ak
            para(30) = kch3
            para(31) = aa1
            para(32) = aa2
            para(33) = aa3
            para(34) = aa4
        end subroutine pjg

        !> @brief Evaluate ethane for the C2H7_VBMM potential.
        subroutine ethane( np, para, rxncrd, &
                & c1, h2, h3, h4, c5, h6, h7, h8, vmm3, etotalmm )
            use chempotpy_kinds, only: wp
            implicit none
            ! Input variables............
            integer  :: np
            real(wp) :: para(np)
            real(wp) :: rxncrd
            real(wp) :: c1(3)
            real(wp) :: h2(3)
            real(wp) :: h3(3)
            real(wp) :: h4(3)
            real(wp) :: c5(3)
            real(wp) :: h6(3)
            real(wp) :: h7(3)
            real(wp) :: h8(3)
            ! Output variables............
            real(wp) :: vmm3
            real(wp) :: estr
            real(wp) :: ebend
            real(wp) :: estrbnd
            real(wp) :: ebndbnd
            real(wp) :: etor
            real(wp) :: estrtor
            real(wp) :: evdwaal
            real(wp) :: etotalmm
            ! parameters:
            ! Angles in Deg.
            ! Dist in Ang.
            real(wp) :: r0ch, r0cc
            real(wp) :: ksch, kscc
            real(wp) :: t0hch, t0cch
            real(wp) :: t0hcrh, t0ccrh
            real(wp) :: kthch, ktcch
            real(wp) :: ksbhch, ksbcch
            real(wp) :: kbbhch, kbbcch
            real(wp) :: v1hcch, v2hcch, v3hcch
            real(wp) :: ksthcch
            real(wp) :: esph, espc, rvh, rvc, scalech
            ! Func called:
            ! Local parameters:
            integer, parameter  :: iout     = 6
            logical, parameter  :: lreact   = .false.
            logical, parameter  :: lradical = .false.
            logical, parameter  :: lswitch  = .true.
            logical, parameter  :: ljg      = .true.
            logical, parameter  :: ldebugmm = .false.
            real(wp), parameter :: pi       = 3.1415926540e0_wp
            ! Stretch varaibles
            real(wp) :: rch(6), rcc
            real(wp) :: esch, escc
            ! Bend varaibles
            real(wp) :: thch(6), tcch(6)
            real(wp) :: ebhch, ebcch
            ! Stretch-Bend varaibles
            real(wp) :: esbcch
            ! Bend-Bend varaibles
            real(wp) :: ebbhch, ebbcch
            ! Torsional variables
            real(wp) :: phcch(9)
            real(wp) :: ethcch
            ! Stretch-xtorsion variables
            real(wp) :: esthcch
            ! van der Waal varaibles
            real(wp) :: sh2(3)
            real(wp) :: sh3(3)
            real(wp) :: sh4(3)
            real(wp) :: sh5(3)
            real(wp) :: sh6(3)
            real(wp) :: sh7(3)
            real(wp) :: sh8(3)
            real(wp) :: rvdwhh(9)
            real(wp) :: evdwhh
            ! energy compponents
            real(wp) :: astr(7)
            real(wp) :: abhch(6)
            real(wp) :: abcch(6)
            real(wp) :: asbcch(6)
            real(wp) :: abbhch(6)
            real(wp) :: abbcch(6)
            real(wp) :: athcch(9)
            real(wp) :: asthcch(9)
            real(wp) :: avdwhh(9)
            real(wp) :: tstr(7)
            real(wp) :: tbhch(6)
            real(wp) :: tbcch(6)
            real(wp) :: tsbcch(6)
            real(wp) :: tbbhch(6)
            real(wp) :: tbbcch(6)
            real(wp) :: tthcch(9)
            real(wp) :: tsthcch(9)
            real(wp) :: tvdwhh(9)
            ! Switching func
            real(wp)            :: scale
            real(wp)            :: sw1(3)
            real(wp)            :: sw2
            real(wp)            :: swt1(3)
            real(wp)            :: swt2
            real(wp)            :: r0crh
            real(wp)            :: r0ccr
            real(wp)            :: r0chsw
            real(wp)            :: r0ccsw
            real(wp)            :: t0hchsw
            real(wp)            :: t0cchsw
            real(wp), parameter :: t1 = 22.0e0_wp
            real(wp), parameter :: t2 = 1.354030e0_wp
            ! Misc varaibles
            integer  :: i, j, k
            real(wp) :: xrv, xesp
            real(wp) :: xsumrx
            real(wp) :: xsumjg
            real(wp) :: xstr
            real(wp) :: xbend
            real(wp) :: xstrbnd
            real(wp) :: xbndbnd
            real(wp) :: xtor
            real(wp) :: xstrtor
            real(wp) :: xvdw
            real(wp) :: z1
            real(wp) :: z2
            real(wp) :: z3
            real(wp) :: z4
            real(wp) :: y1
            real(wp) :: y2
            real(wp) :: y3
            ! Assign parameters
            r0ch    = para(1)
            r0cc    = para(2)
            ksch    = para(3)
            kscc    = para(4)
            t0hch   = para(5)
            t0cch   = para(6)
            t0hcrh  = para(7)
            t0ccrh  = para(8)
            kthch   = para(9)
            ktcch   = para(10)
            ksbhch  = para(11)
            ksbcch  = para(12)
            kbbhch  = para(13)
            kbbcch  = para(14)
            v1hcch  = para(15)
            v2hcch  = para(16)
            v3hcch  = para(17)
            ksthcch = para(18)
            espc    = para(19)
            rvc     = para(20)
            esph    = para(21)
            rvh     = para(22)
            scalech = para(23)
            r0crh   = para(24)
            r0ccr   = para(25)
            !     write(IOUT,*)
            !     write(IOUT,*)'LREACT   :',LREACT
            !     write(IOUT,*)'LRADICAL :',LRADICAL
            ! Evaluate bond xdist
            rch(1) = xdist(c1, h2)
            rch(2) = xdist(c1, h3)
            rch(3) = xdist(c1, h4)
            rch(4) = xdist(c5, h6)
            rch(5) = xdist(c5, h7)
            rch(6) = xdist(c5, h8)
            rcc    = xdist(c1, c5)
            ! calculate the swithing func
            swt2 = 1.0e0_wp
            do i=1, 3
                swt1(i) = (1.0e0_wp-tanh(t1*(rch(i)-t2)))*0.50e0_wp
                swt2    = swt2*swt1(i)
            enddo
            ! switch C-C xdist
            r0cc = (swt2*r0cc)+((1.0e0_wp-swt2)*r0ccr)
            !%%%
            ! switch force constants
            z1    = 10.0e0_wp
            z2    = 2.144106410736830e0_wp
            z3    = exp(-z1*((rxncrd-z2)**2))
            scale = 1.0e0_wp-(0.050e0_wp*z3)
            !     scale=0.950d0
            ksch  = ksch*scale
            kthch = kthch*scale
            ktcch = ktcch*scale
            ! calculate str energy
            call es(r0cc, kscc, rcc, escc)
            estr = escc
            do i=1, 6
                call es(r0ch, ksch, rch(i), esch)
                estr    = estr + esch
                astr(i) = esch
            enddo
            astr(7) = escc
            ! calculate bend energy
            ! there are 12 bending terms
            thch(1) = xangle(h2, c1, h3)*180.0e0_wp/pi
            thch(2) = xangle(h2, c1, h4)*180.0e0_wp/pi
            thch(3) = xangle(h3, c1, h4)*180.0e0_wp/pi
            thch(4) = xangle(h6, c5, h7)*180.0e0_wp/pi
            thch(5) = xangle(h6, c5, h8)*180.0e0_wp/pi
            thch(6) = xangle(h7, c5, h8)*180.0e0_wp/pi
            tcch(1) = xangle(h2, c1, c5)*180.0e0_wp/pi
            tcch(2) = xangle(h3, c1, c5)*180.0e0_wp/pi
            tcch(3) = xangle(h4, c1, c5)*180.0e0_wp/pi
            tcch(4) = xangle(h6, c5, c1)*180.0e0_wp/pi
            tcch(5) = xangle(h7, c5, c1)*180.0e0_wp/pi
            tcch(6) = xangle(h8, c5, c1)*180.0e0_wp/pi
            ebend   = 0.0e0_wp
            do i=1, 6
                call eb(t0hch, kthch, thch(i), ebhch)
                call eb(t0cch, ktcch, tcch(i), ebcch)
                ebend    = ebend + ebhch + ebcch
                abhch(i) = ebhch
                abcch(i) = ebcch
            enddo
            ! Calculate Stretch-bend terms
            ! there are 6 stretch-bend cross terms
            ! H2-C1-C5
            ! H3-C1-C5
            ! H4-C1-C5
            ! H6-C5-C1
            ! H7-C5-C1
            ! H8-C5-C1
            estrbnd = 0.0e0_wp
            do i=1, 6
                call esb(r0ch, r0cc, t0cch, ksbcch, rch(i), rcc, tcch(i), esbcch)
                estrbnd   = estrbnd + esbcch
                asbcch(i) = esbcch
            enddo
            ! Bend-Bend terms
            ! there are 12 bend-bend terms
            ! H2-C1-H3 & H2-C1-H4
            ! H2-C1-H3 & H3-C1-H4
            ! H2-C1-H4 & H3-C1-H4
            ! H6-C5-H7 & H6-C5-H8
            ! H6-C5-H7 & H7-C5-H8
            ! H6-C5-H8 & H7-C5-H8
            ! H2-C1-C5 & H3-C1-C5
            ! H2-C1-C5 & H4-C1-C5
            ! H3-C1-C5 & H4-C1-C5
            ! H6-C5-C1 & H6-C5-C1
            ! H6-C5-C1 & H7-C5-C1
            ! H7-C5-C1 & H8-C5-C1
            ebndbnd = 0.0e0_wp
            k       = 0
            do i=1, 3
                do j=i+1, 3
                    k = k + 1
                    call ebb(t0hch, t0hch, kbbhch, thch(i), thch(j), ebbhch)
                    call ebb(t0cch, t0cch, kbbcch, tcch(i), tcch(j), ebbcch)
                    ebndbnd   = ebndbnd + ebbhch + ebbcch
                    abbhch(k) = ebbhch
                    abbcch(k) = ebbcch
                enddo
            enddo
            do i=4, 6
                do j=i+1, 6
                    k = k + 1
                    call ebb(t0hch, t0hch, kbbhch, thch(i), thch(j), ebbhch)
                    call ebb(t0cch, t0cch, kbbcch, tcch(i), tcch(j), ebbcch)
                    ebndbnd   = ebndbnd + ebbhch + ebbcch
                    abbhch(k) = ebbhch
                    abbcch(k) = ebbcch
                enddo
            enddo
            ! there are 9 xtorsional terms
            ! H2-C1-C5-H6
            ! H2-C1-C5-H7
            ! H2-C1-C5-H8
            ! H3-C1-C5-H6
            ! H3-C1-C5-H7
            ! H3-C1-C5-H8
            ! H4-C1-C5-H6
            ! H4-C1-C5-H7
            ! H4-C1-C5-H8
            ! Scale the C-H bonds and evaluate the
            ! coordinates of the scaled H
            phcch(1) = xtorsion(h2, c1, c5, h6)
            phcch(2) = xtorsion(h2, c1, c5, h7)
            phcch(3) = xtorsion(h2, c1, c5, h8)
            phcch(4) = xtorsion(h3, c1, c5, h6)
            phcch(5) = xtorsion(h3, c1, c5, h7)
            phcch(6) = xtorsion(h3, c1, c5, h8)
            phcch(7) = xtorsion(h4, c1, c5, h6)
            phcch(8) = xtorsion(h4, c1, c5, h7)
            phcch(9) = xtorsion(h4, c1, c5, h8)
            etor     = 0.0e0_wp
            do i=1, 9
                call et(v1hcch, v2hcch, v3hcch, phcch(i), ethcch)
                etor      = etor + ethcch
                athcch(i) = ethcch
                if( ldebugmm )then
                    write(*, *) 'Torsion Angle:', phcch(i)
                endif
            enddo
            ! stretch-xtorsion terms
            ! interaction of the C-C str with xtorsion
            ! there are 9 such terms
            ! H2-C1-C5-H6
            ! H2-C1-C5-H7
            ! H2-C1-C5-H8
            ! H3-C1-C5-H6
            ! H3-C1-C5-H7
            ! H3-C1-C5-H8
            ! H4-C1-C5-H6
            ! H4-C1-C5-H7
            ! H4-C1-C5-H8
            estrtor = 0.0e0_wp
            do i=1, 9
                call est(r0cc, ksthcch, rcc, phcch(i), esthcch)
                estrtor    = estrtor + esthcch
                asthcch(i) = esthcch
            enddo
            ! Van der Waals terms
            ! all interaction which are greater than 1-3
            ! are van der Waal intereactions
            ! there are 9 vdw terms
            ! H2-H6
            ! H2-H7
            ! H2-H8
            ! H3-H6
            ! H3-H7
            ! H3-H8
            ! H4-H6
            ! H4-H7
            ! H4-H8
            do i=1, 3
                sh2(i) = (scalech*h2(i)) + ((1.0e0_wp-scalech)*c1(i))
                sh3(i) = (scalech*h3(i)) + ((1.0e0_wp-scalech)*c1(i))
                sh4(i) = (scalech*h4(i)) + ((1.0e0_wp-scalech)*c1(i))
                sh6(i) = (scalech*h6(i)) + ((1.0e0_wp-scalech)*c5(i))
                sh7(i) = (scalech*h7(i)) + ((1.0e0_wp-scalech)*c5(i))
                sh8(i) = (scalech*h8(i)) + ((1.0e0_wp-scalech)*c5(i))
            enddo
            rvdwhh(1) = xdist(sh2, sh6)
            rvdwhh(2) = xdist(sh2, sh7)
            rvdwhh(3) = xdist(sh2, sh8)
            rvdwhh(4) = xdist(sh3, sh6)
            rvdwhh(5) = xdist(sh3, sh7)
            rvdwhh(6) = xdist(sh3, sh8)
            rvdwhh(7) = xdist(sh4, sh6)
            rvdwhh(8) = xdist(sh4, sh7)
            rvdwhh(9) = xdist(sh4, sh8)
            evdwaal   = 0.0e0_wp
            do i=1, 9
                xrv  = rvh + rvh
                xesp = sqrt(esph*esph)
                call evdw(xrv, xesp, rvdwhh(i), evdwhh)
                evdwaal   = evdwaal + evdwhh
                avdwhh(i) = evdwhh
            enddo
            etotalmm = estr+ebend+estrbnd+ebndbnd+ etor+estrtor+evdwaal
            ! Switch between the parameters of ethane and ethyl radical
            ! swt2=1 for ethane and 0 for ethyl
            r0ccsw  = r0cc
            r0chsw  = (swt2*r0ch)+((1.0e0_wp-swt2)*r0crh)
            t0hchsw = (swt2*t0hch)+((1.0e0_wp-swt2)*t0hcrh)
            t0cchsw = (swt2*t0cch)+((1.0e0_wp-swt2)*t0ccrh)
            !%%%
            z2  = 2.141145010007599490e0_wp
            z2  = 1.50e0_wp
            z1  = 0.50e0_wp*(1.0e0_wp+tanh(6.0e0_wp*(rxncrd-z2)))
            z4  = ((1.0e0_wp-z1)*t0cch)+(z1*t0ccrh)
            z1  = 117.2989954203791430e0_wp
            sw2 = 1.0e0_wp
            do i=1, 3
                z1     = 10.0e0_wp
                sw1(i) = (1.0e0_wp-tanh(z1*(rch(i)-t2)))*0.50e0_wp
                sw2    = sw2*sw1(i)
            enddo
            t0cchsw = (sw2*t0cch)+((1.0e0_wp-sw2)*t0ccrh)
            if( lswitch )then
                do i=1, 3
                    call eb(t0hchsw, kthch, thch(i), ebhch)
                    call eb(t0cchsw, ktcch, tcch(i), ebcch)
                    call esb(r0chsw, r0ccsw, t0cchsw, ksbcch, &
                        & rch(i), rcc, tcch(i), esbcch)
                    abhch(i)  = ebhch
                    abcch(i)  = ebcch
                    asbcch(i) = esbcch
                enddo
                k = 0
                do i=1, 3
                    do j=i+1, 3
                        k = k + 1
                        call ebb(t0hchsw, t0hchsw, kbbhch, thch(i), thch(j), ebbhch)
                        call ebb(t0cchsw, t0cchsw, kbbcch, tcch(i), tcch(j), ebbcch)
                        abbhch(k) = ebbhch
                        abbcch(k) = ebbcch
                    enddo
                enddo
                ! compute engergies
                do i=1, 7
                    tstr(i) = astr(i)
                enddo
                do i=1, 6
                    tbhch(i)  = abhch(i)
                    tbcch(i)  = abcch(i)
                    tsbcch(i) = asbcch(i)
                    tbbhch(i) = abbhch(i)
                    tbbcch(i) = abbcch(i)
                enddo
                do i=1, 9
                    tthcch(i)  = athcch(i)
                    tsthcch(i) = asthcch(i)
                    tvdwhh(i)  = avdwhh(i)
                enddo
                !       switch C-H str term
                tstr(1) = tstr(1)*swt1(1)
                tstr(2) = tstr(2)*swt1(2)
                tstr(3) = tstr(3)*swt1(3)
                !       switch the following H-C-H bend terms
                !       H2-C1-H3
                !       H2-C1-H4
                !       H3-C1-H4
                tbhch(1) = tbhch(1)*swt1(1)*swt1(2)
                tbhch(2) = tbhch(2)*swt1(1)*swt1(3)
                tbhch(3) = tbhch(3)*swt1(2)*swt1(3)
                !       switch the following C-C-H bend terms
                !       C5-C1-H2
                !       C5-C1-H3
                !       C5-C1-H4
                tbcch(1) = tbcch(1)*swt1(1)
                tbcch(2) = tbcch(2)*swt1(2)
                tbcch(3) = tbcch(3)*swt1(3)
                !       switch the following C-C-H strbnd term
                !       C5-C1-H2
                !       C5-C1-H3
                !       C5-C1-H4
                tsbcch(1) = tsbcch(1)*swt1(1)
                tsbcch(2) = tsbcch(2)*swt1(2)
                tsbcch(3) = tsbcch(3)*swt1(3)
                !       remove three H-C-H2 bnd-bnd terms
                tbbhch(1) = tbbhch(1)*swt1(1)*swt1(2)
                tbbhch(2) = tbbhch(2)*swt1(1)*swt1(3)
                tbbhch(3) = tbbhch(3)*swt1(2)*swt1(3)
                !       remove two C-C-H2 bnd-bnd terms
                tbbcch(1) = tbbcch(1)*swt1(1)*swt1(2)
                tbbcch(2) = tbbcch(2)*swt1(1)*swt1(3)
                tbbcch(3) = tbbcch(3)*swt1(2)*swt1(3)
                !       switch the following  H-C-C-H xtorsional terms
                !       H2-C1-C5-H6
                !       H2-C1-C5-H7
                !       H2-C1-C5-H8
                !       H3-C1-C5-H6
                !       H3-C1-C5-H7
                !       H3-C1-C5-H8
                !       H4-C1-C5-H6
                !       H4-C1-C5-H7
                !       H4-C1-C5-H8
                tthcch(1) = tthcch(1)*swt1(1)
                tthcch(2) = tthcch(2)*swt1(1)
                tthcch(3) = tthcch(3)*swt1(1)
                tthcch(4) = tthcch(4)*swt1(2)
                tthcch(5) = tthcch(5)*swt1(2)
                tthcch(6) = tthcch(6)*swt1(2)
                tthcch(7) = tthcch(7)*swt1(3)
                tthcch(8) = tthcch(8)*swt1(3)
                tthcch(9) = tthcch(9)*swt1(3)
                !       switch the following H-C-C-H str-xtorsional terms
                !       H2-C1-C5-H6
                !       H2-C1-C5-H7
                !       H2-C1-C5-H8
                !       H3-C1-C5-H6
                !       H3-C1-C5-H7
                !       H3-C1-C5-H8
                !       H4-C1-C5-H6
                !       H4-C1-C5-H7
                !       H4-C1-C5-H8
                tsthcch(1) = tsthcch(1)*swt1(1)
                tsthcch(2) = tsthcch(2)*swt1(1)
                tsthcch(3) = tsthcch(3)*swt1(1)
                tsthcch(4) = tsthcch(4)*swt1(2)
                tsthcch(5) = tsthcch(5)*swt1(2)
                tsthcch(6) = tsthcch(6)*swt1(2)
                tsthcch(7) = tsthcch(7)*swt1(3)
                tsthcch(8) = tsthcch(8)*swt1(3)
                tsthcch(9) = tsthcch(9)*swt1(3)
                !       switch the following vdw terms
                !       H2-H6
                !       H2-H7
                !       H2-H8
                !       H3-H6
                !       H3-H7
                !       H3-H8
                !       H4-H6
                !       H4-H7
                !       H4-H8
                tvdwhh(1) = tvdwhh(1)*swt1(1)
                tvdwhh(2) = tvdwhh(2)*swt1(1)
                tvdwhh(3) = tvdwhh(3)*swt1(1)
                tvdwhh(4) = tvdwhh(4)*swt1(2)
                tvdwhh(5) = tvdwhh(5)*swt1(2)
                tvdwhh(6) = tvdwhh(6)*swt1(2)
                tvdwhh(7) = tvdwhh(7)*swt1(3)
                tvdwhh(8) = tvdwhh(8)*swt1(3)
                tvdwhh(9) = tvdwhh(9)*swt1(3)
                xsumrx    = 0.0e0_wp
                xstr      = 0.0e0_wp
                do i=1, 7
                    xstr = xstr+tstr(i)
                enddo
                xbend   = 0.0e0_wp
                xstrbnd = 0.0e0_wp
                xbndbnd = 0.0e0_wp
                do i=1, 6
                    xbend   = xbend+tbhch(i)+tbcch(i)
                    xstrbnd = xstrbnd+tsbcch(i)
                    xbndbnd = xbndbnd+tbbhch(i)+tbbcch(i)
                enddo
                xtor    = 0.0e0_wp
                xstrtor = 0.0e0_wp
                xvdw    = 0.0e0_wp
                do i=1, 9
                    xtor    = xtor+tthcch(i)
                    xstrtor = xstrtor+tsthcch(i)
                    xvdw    = xvdw+tvdwhh(i)
                enddo
                xsumrx = xstr+xbend+xstrbnd+xbndbnd+xtor+xstrtor+xvdw
                !       call xdisp(xstr,xbend,xstrbnd,xbndbnd,
                !    1             xtor,xstrtor,xvdw,xsumRX)
                etotalmm = xsumrx
            endif
            if( lradical )then
                !       modify the h-cR-h and h-cR-C bends and c-cR-h str-bnd
                !       for C-sp3 bonded to only two atoms
                do i=1, 3
                    call eb(t0hcrh, kthch, thch(i), ebhch)
                    call eb(t0ccrh, ktcch, tcch(i), ebcch)
                    call esb(r0ch, r0cc, t0ccrh, ksbcch, rch(i), rcc, tcch(i), esbcch)
                    abhch(i)  = ebhch
                    abcch(i)  = ebcch
                    asbcch(i) = esbcch
                enddo
                k = 0
                do i=1, 3
                    do j=i+1, 3
                        k = k + 1
                        call ebb(t0hcrh, t0hcrh, kbbhch, thch(i), thch(j), ebbhch)
                        call ebb(t0ccrh, t0ccrh, kbbcch, tcch(i), tcch(j), ebbcch)
                        abbhch(k) = ebbhch
                        abbcch(k) = ebbcch
                    enddo
                enddo
            endif
            if( lreact )then
                do i=1, 7
                    tstr(i) = astr(i)
                enddo
                do i=1, 6
                    tbhch(i)  = abhch(i)
                    tbcch(i)  = abcch(i)
                    tsbcch(i) = asbcch(i)
                    tbbhch(i) = abbhch(i)
                    tbbcch(i) = abbcch(i)
                enddo
                do i=1, 9
                    tthcch(i)  = athcch(i)
                    tsthcch(i) = asthcch(i)
                    tvdwhh(i)  = avdwhh(i)
                enddo
                !       remove one C-H2 str term
                tstr(1) = 0.0e0_wp
                !       remove two H-C-H2 bend terms
                tbhch(1) = 0.0e0_wp
                tbhch(2) = 0.0e0_wp
                !       remove one C-C-H2 bend term
                tbcch(1) = 0.0e0_wp
                !       remove one C-C-H2 strbnd term
                tsbcch(1) = 0.0e0_wp
                !       remove three H-C-H2 bnd-bnd terms
                tbbhch(1) = 0.0e0_wp
                tbbhch(2) = 0.0e0_wp
                tbbhch(3) = 0.0e0_wp
                !       remove two C-C-H2 bnd-bnd terms
                tbbcch(1) = 0.0e0_wp
                tbbcch(2) = 0.0e0_wp
                !       remove three H-C-C-H2 xtorsional terms
                tthcch(1) = 0.0e0_wp
                tthcch(2) = 0.0e0_wp
                tthcch(3) = 0.0e0_wp
                !       remove three H-C-C-H2 str-xtorsional terms
                tsthcch(1) = 0.0e0_wp
                tsthcch(2) = 0.0e0_wp
                tsthcch(3) = 0.0e0_wp
                !       remove three vdw terms
                tvdwhh(1) = 0.0e0_wp
                tvdwhh(2) = 0.0e0_wp
                tvdwhh(3) = 0.0e0_wp
                xsumrx    = 0.0e0_wp
                xstr      = 0.0e0_wp
                do i=1, 7
                    xstr = xstr+tstr(i)
                enddo
                xbend   = 0.0e0_wp
                xstrbnd = 0.0e0_wp
                xbndbnd = 0.0e0_wp
                do i=1, 6
                    xbend   = xbend+tbhch(i)+tbcch(i)
                    xstrbnd = xstrbnd+tsbcch(i)
                    xbndbnd = xbndbnd+tbbhch(i)+tbbcch(i)
                enddo
                xtor    = 0.0e0_wp
                xstrtor = 0.0e0_wp
                xvdw    = 0.0e0_wp
                do i=1, 9
                    xtor    = xtor+tthcch(i)
                    xstrtor = xstrtor+tsthcch(i)
                    xvdw    = xvdw+tvdwhh(i)
                enddo
                xsumrx = xstr+xbend+xstrbnd+xbndbnd+xtor+xstrtor+xvdw
                !       call xdisp(xstr,xbend,xstrbnd,xbndbnd,
                !    1             xtor,xstrtor,xvdw,xsumRX)
                etotalmm = xsumrx
            endif
            ! start of if-block label: jg_case
            if( ljg )then
                do i=1, 7
                    tstr(i) = astr(i)
                enddo
                do i=1, 6
                    tbhch(i)  = abhch(i)
                    tbcch(i)  = abcch(i)
                    tsbcch(i) = asbcch(i)
                    tbbhch(i) = abbhch(i)
                    tbbcch(i) = abbcch(i)
                enddo
                do i=1, 9
                    tthcch(i)  = athcch(i)
                    tsthcch(i) = asthcch(i)
                    tvdwhh(i)  = avdwhh(i)
                enddo
                !       remove the following C-H str terms
                !       C-H2
                !       C-H3
                !       C-H4
                tstr(1) = 0.0e0_wp
                tstr(2) = 0.0e0_wp
                tstr(3) = 0.0e0_wp
                !       remove the following H-C-H bend terms
                !       H2-C1-H3
                !       H2-C1-H4
                !       H3-C1-H4
                tbhch(1) = 0.0e0_wp
                tbhch(2) = 0.0e0_wp
                tbhch(3) = 0.0e0_wp
                !       switch the following C-C-H bend terms
                !       C5-C1-H2
                !       C5-C1-H3
                !       C5-C1-H4
                tbcch(1) = tbcch(1)*swt1(1)
                tbcch(2) = tbcch(2)*swt1(2)
                tbcch(3) = tbcch(3)*swt1(3)
                !       switch the following C-C-H strbnd term
                !       C5-C1-H2
                !       C5-C1-H3
                !       C5-C1-H4
                tsbcch(1) = tsbcch(1)*swt1(1)
                tsbcch(2) = tsbcch(2)*swt1(2)
                tsbcch(3) = tsbcch(3)*swt1(3)
                !       remove three H-C-H2 bnd-bnd terms
                tbbhch(1) = 0.0e0_wp
                tbbhch(2) = 0.0e0_wp
                tbbhch(3) = 0.0e0_wp
                !       switch C-C-H bnd-bnd terms involving H2, H3 and H4
                tbbcch(1) = tbbcch(1)*swt1(1)
                tbbcch(2) = tbbcch(2)*swt1(2)
                tbbcch(3) = tbbcch(3)*swt1(3)
                !       switch the following  H-C-C-H xtorsional terms
                !       H2-C1-C5-H6
                !       H2-C1-C5-H7
                !       H2-C1-C5-H8
                !       H3-C1-C5-H6
                !       H3-C1-C5-H7
                !       H3-C1-C5-H8
                !       H4-C1-C5-H6
                !       H4-C1-C5-H7
                !       H4-C1-C5-H8
                tthcch(1) = tthcch(1)*swt1(1)
                tthcch(2) = tthcch(2)*swt1(1)
                tthcch(3) = tthcch(3)*swt1(1)
                tthcch(4) = tthcch(4)*swt1(2)
                tthcch(5) = tthcch(5)*swt1(2)
                tthcch(6) = tthcch(6)*swt1(2)
                tthcch(7) = tthcch(7)*swt1(3)
                tthcch(8) = tthcch(8)*swt1(3)
                tthcch(9) = tthcch(9)*swt1(3)
                !       switch the following H-C-C-H str-xtorsional terms
                !       H2-C1-C5-H6
                !       H2-C1-C5-H7
                !       H2-C1-C5-H8
                !       H3-C1-C5-H6
                !       H3-C1-C5-H7
                !       H3-C1-C5-H8
                !       H4-C1-C5-H6
                !       H4-C1-C5-H7
                !       H4-C1-C5-H8
                tsthcch(1) = tsthcch(1)*swt1(1)
                tsthcch(2) = tsthcch(2)*swt1(1)
                tsthcch(3) = tsthcch(3)*swt1(1)
                tsthcch(4) = tsthcch(4)*swt1(2)
                tsthcch(5) = tsthcch(5)*swt1(2)
                tsthcch(6) = tsthcch(6)*swt1(2)
                tsthcch(7) = tsthcch(7)*swt1(3)
                tsthcch(8) = tsthcch(8)*swt1(3)
                tsthcch(9) = tsthcch(9)*swt1(3)
                !       switch the following vdw terms
                !       H2-H6
                !       H2-H7
                !       H2-H8
                !       H3-H6
                !       H3-H7
                !       H3-H8
                !       H4-H6
                !       H4-H7
                !       H4-H8
                tvdwhh(1) = tvdwhh(1)*swt1(1)
                tvdwhh(2) = tvdwhh(2)*swt1(1)
                tvdwhh(3) = tvdwhh(3)*swt1(1)
                tvdwhh(4) = tvdwhh(4)*swt1(2)
                tvdwhh(5) = tvdwhh(5)*swt1(2)
                tvdwhh(6) = tvdwhh(6)*swt1(2)
                tvdwhh(7) = tvdwhh(7)*swt1(3)
                tvdwhh(8) = tvdwhh(8)*swt1(3)
                tvdwhh(9) = tvdwhh(9)*swt1(3)
                xsumjg    = 0.0e0_wp
                xstr      = 0.0e0_wp
                do i=1, 7
                    xstr = xstr+tstr(i)
                enddo
                xbend   = 0.0e0_wp
                xstrbnd = 0.0e0_wp
                xbndbnd = 0.0e0_wp
                do i=1, 6
                    xbend   = xbend+tbhch(i)+tbcch(i)
                    xstrbnd = xstrbnd+tsbcch(i)
                    xbndbnd = xbndbnd+tbbhch(i)+tbbcch(i)
                enddo
                xtor    = 0.0e0_wp
                xstrtor = 0.0e0_wp
                xvdw    = 0.0e0_wp
                do i=1, 9
                    xtor    = xtor+tthcch(i)
                    xstrtor = xstrtor+tsthcch(i)
                    xvdw    = xvdw+tvdwhh(i)
                enddo
                xsumjg = xstr+xbend+xstrbnd+xbndbnd+xtor+xstrtor+xvdw
                if( ldebugmm )then
                    write(*, '(A,F14.8)')
                    write(*, '(A,F14.8)') '----------------------------'
                    write(*, '(A,F14.8)') 'Energy from MM3'
                    write(*, '(A,F14.8)') '----------------------------'
                    call xdisp(xstr, xbend, xstrbnd, xbndbnd, &
                        & xtor, xstrtor, xvdw, xsumjg)
                endif
            endif
            ! end of if-block label: jg_case
            vmm3 = xsumjg
        end subroutine ethane

        !> @brief Evaluate xdisp for the C2H7_VBMM potential.
        subroutine xdisp( xstr, xbend, xstrbnd, xbndbnd, &
                & xtor, xstrtor, xvdw, xsum )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: xstr
            real(wp) :: xbend
            real(wp) :: xstrbnd
            real(wp) :: xbndbnd
            real(wp) :: xtor
            real(wp) :: xstrtor
            real(wp) :: xvdw
            real(wp) :: xsum
            write(*, '(A,F14.4)')
            write(*, '(A,F14.4)') 'Energy in kcal/mol'
            write(*, '(A,F14.4)')'Stretch   :', xstr
            write(*, '(A,F14.4)')'Bend      :', xbend
            write(*, '(A,F14.4)')'Str-Bnd   :', xstrbnd
            write(*, '(A,F14.4)')'Bnd-Bnd   :', xbndbnd
            write(*, '(A,F14.4)')'Torsional :', xtor
            write(*, '(A,F14.4)')'Str-Tor   :', xstrtor
            write(*, '(A,F14.4)')'vdW       :', xvdw
            write(*, '(A,F14.4)')
            write(*, '(A,F14.4)')'Total MM eng.  :', xsum
        end subroutine xdisp

        !> @brief Evaluate xdist for the C2H7_VBMM potential.
        function xdist( x, y )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: x(3)
            real(wp) :: y(3)
            real(wp) :: xdist
            integer  :: i
            real(wp) :: z
            z = 0.0e0_wp
            do i=1, 3
                z = z + ( (y(i)-x(i))**2 )
            enddo
            xdist = sqrt(z)
        end function xdist

        !> @brief Evaluate xangle for the C2H7_VBMM potential.
        function xangle( a, b, c )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3)
            real(wp) :: b(3)
            real(wp) :: c(3)
            real(wp) :: xangle
            integer  :: i
            real(wp) :: sum, rab(3), rcb(3)
            do i=1, 3
                rab(i) = a(i)-b(i)
                rcb(i) = c(i)-b(i)
            enddo
            sum = 0.0e0_wp
            do i=1, 3
                sum = sum + (rab(i)*rcb(i))
            enddo
            sum    = sum/(xdist(a, b)*xdist(c, b))
            xangle = abs(acos(sum))
        end function xangle

        !> @brief Evaluate xtorsion for the C2H7_VBMM potential.
        function xtorsion( a, b, c, d )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3)
            real(wp) :: b(3)
            real(wp) :: c(3)
            real(wp) :: d(3)
            real(wp) :: xtorsion
            ! This func calculates the dihedral xangle
            ! between a-b and c-d
            ! a
            ! |
            ! b---c
            !     |
            !     d
            ! we will use the cross product definition of the
            ! dihedral xangle
            integer  :: i, j
            real(wp) :: x, y, z
            real(wp) :: s
            real(wp) :: phi
            real(wp) :: dpuv
            real(wp) :: rab(3)
            real(wp) :: rcb(3)
            real(wp) :: rcd(3)
            real(wp) :: u(3)
            real(wp) :: v(3)
            ! define bond-vectors
            do i=1, 3
                rab(i) = a(i)-b(i)
                rcb(i) = c(i)-b(i)
                rcd(i) = c(i)-d(i)
            enddo
            ! define the normal vectors
            ! u = rab X rcb
            ! v = rcb X rcd
            call crossprod(rab, rcb, u)
            call crossprod(rcb, rcd, v)
            ! evaluate the norm and dot product
            ! u and v
            ! dpuv = u . v
            ! s = rab . v and is used for sign of
            ! the dihedral xangle
            x    = 0.0e0_wp
            y    = 0.0e0_wp
            s    = 0.0e0_wp
            dpuv = 0.0e0_wp
            do i=1, 3
                x    = x+(u(i)*u(i))
                y    = y+(v(i)*v(i))
                dpuv = dpuv + (u(i)*v(i))
                s    = s + (rab(i)*v(i))
            enddo
            x = sqrt(x)
            y = sqrt(y)
            ! z is cos(phi) and should be in the
            ! range [-1,1]
            if( x*y .ne. 0.0e0_wp )then
                z = dpuv/(x*y)
            else
                z = dpuv
            endif
            if( z .gt. 1.0e0_wp ) z=1.0e0_wp
            if( z .lt. -1.0e0_wp ) z=-1.0e0_wp
            ! get dihedral xangle
            phi = acos(z)
            ! get the sign
            if( s /= 0.0e0_wp )then
                s = s/abs(s)
            else
                s = 1.0e0_wp
            endif
            xtorsion = s * phi
        end function xtorsion

        !> @brief Evaluate crossprod for the C2H7_VBMM potential.
        subroutine crossprod( a, b, c )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3)
            real(wp) :: b(3)
            real(wp) :: c(3)
            c(1) = (a(2)*b(3))-(a(3)*b(2))
            c(2) = (a(3)*b(1))-(a(1)*b(3))
            c(3) = (a(1)*b(2))-(a(2)*b(1))
        end subroutine crossprod

        !> @brief Evaluate Es for the C2H7_VBMM potential.
        subroutine es( r0, ks, r, eng )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: r0
            real(wp)            :: ks
            real(wp)            :: r
            real(wp)            :: eng
            real(wp)            :: dr
            real(wp), parameter :: c1 = 71.940e0_wp
            real(wp), parameter :: c2 = 2.550e0_wp
            real(wp), parameter :: c3 = 7.0e0_wp/12.0e0_wp
            eng = 0.0e0_wp
            dr  = (r-r0)
            eng = (c1*ks*dr*dr) * (1.0e0_wp-(c2*dr)+(c3*c2*c2*dr*dr))
            !    1 (1.0d0-(c2*dr)+(c3*c2*dr*dr))
        end subroutine es

        !> @brief Evaluate Eb for the C2H7_VBMM potential.
        subroutine eb( t0, kt, t, eng )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: t0
            real(wp) :: kt
            real(wp) :: t
            real(wp) :: eng
            real(wp) :: dt
            ! Alligers (JACS, vol 111, pg 8551) parameters
            !     parameter(c1=0.0219140d0)
            !     parameter(c2=0.0140d0)
            !     parameter(c3=5.60d-5)
            !     parameter(c4=7.0d-7)
            !     parameter(c5=9.0d-10)
            ! Tinker 4.1 MM3 parameter
            real(wp), parameter :: c1 = 0.021914180e0_wp
            real(wp), parameter :: c2 = 0.0140e0_wp
            real(wp), parameter :: c3 = 0.0000560e0_wp
            real(wp), parameter :: c4 = 0.00000070e0_wp
            real(wp), parameter :: c5 = 0.0000000220e0_wp
            dt  = t-t0
            eng = c1*kt*dt*dt* ( 1.0e0_wp-(c2*dt)+(c3*dt*dt)-(c4*dt*dt*dt) +(c5*dt*dt*dt*dt) )
        end subroutine eb

        !> @brief Evaluate Et for the C2H7_VBMM potential.
        subroutine et( v1, v2, v3, w, eng )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: v1
            real(wp) :: v2
            real(wp) :: v3
            real(wp) :: w
            real(wp) :: eng
            eng = 0.50e0_wp*( ( v1*(1.0e0_wp+cos(w)) ) + ( v2*(1.0e0_wp-cos(2.0e0_wp*w)) ) + ( v3*(1.0e0_wp+ &
                & cos(3.0e0_wp*w)) ) )
        end subroutine et

        !> @brief Evaluate Esb for the C2H7_VBMM potential.
        subroutine esb( ra0, rb0, t0, kst, ra, rb, t, eng )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: ra0
            real(wp)            :: rb0
            real(wp)            :: t0
            real(wp)            :: kst
            real(wp)            :: ra
            real(wp)            :: rb
            real(wp)            :: t
            real(wp)            :: eng
            real(wp), parameter :: c1 = 2.511180e0_wp
            eng = c1*kst*(t-t0)*((ra-ra0)+(rb-rb0))
        end subroutine esb

        !> @brief Evaluate Est for the C2H7_VBMM potential.
        subroutine est( r0, kws, r, w, eng )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: r0
            real(wp)            :: kws
            real(wp)            :: r
            real(wp)            :: w
            real(wp)            :: eng
            real(wp), parameter :: c1 = -11.9950e0_wp
            ! In Alligers paper(JACS, vol 111, pg 8551)
            ! c1 is +ve but TINKER (version 4.1) defines
            ! c1/2 as -ve
            eng = 0.50e0_wp*c1*kws*(r-r0)*(1.0e0_wp+cos(3.0e0_wp*w))
        end subroutine est

        !> @brief Evaluate Ebb for the C2H7_VBMM potential.
        subroutine ebb( ta0, tb0, kbb, ta, tb, eng )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: ta0
            real(wp)            :: tb0
            real(wp)            :: kbb
            real(wp)            :: ta
            real(wp)            :: tb
            real(wp)            :: eng
            real(wp), parameter :: c1 = -0.021914180e0_wp
            eng = c1*kbb*(ta-ta0)*(tb-tb0)
            !     write(*,1000) c1,kbb,(ta-ta0)*(tb-tb0),eng
        end subroutine ebb

        !> @brief Evaluate Evdw for the C2H7_VBMM potential.
        subroutine evdw( rv, esp, r, eng )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: rv
            real(wp)            :: esp
            real(wp)            :: r
            real(wp)            :: eng
            real(wp)            :: x
            real(wp), parameter :: c1 = 2.250e0_wp
            real(wp), parameter :: c2 = 1.840e+5_wp
            real(wp), parameter :: c3 = 12.00e0_wp
            x   = rv/r
            eng = esp*((-c1*x*x*x*x*x*x)+(c2*exp(-c3/x)))
        end subroutine evdw

        !> @brief Evaluate xch5pot for the C2H7_VBMM potential.
        subroutine xch5pot( np, para, rxncrd, c1, h2, h3, h4, h5, hb, xpot_kcal )
            ! Input : Coordintes are in Ang.
            ! Output: xpot_kcal is in kcal/mol
            ! All energy variables local to this
            ! routine are in 10(5)Joules/mol
            use chempotpy_kinds, only: wp
            implicit none
            ! Input variables
            integer  :: np
            real(wp) :: para(np)
            real(wp) :: rxncrd
            real(wp) :: c1(3)
            real(wp) :: h2(3)
            real(wp) :: h3(3)
            real(wp) :: h4(3)
            real(wp) :: h5(3)
            real(wp) :: hb(3)
            ! Output variables
            real(wp) :: xpot_kcal
            ! Parameter list
            real(wp)            :: a1ch
            real(wp)            :: b1ch
            real(wp)            :: c1ch
            real(wp)            :: alpcb
            real(wp)            :: d1cb
            real(wp)            :: d3cb
            real(wp)            :: r0cb
            real(wp)            :: alphh
            real(wp)            :: d1hh
            real(wp)            :: d3hh
            real(wp)            :: r0hh
            real(wp)            :: alpch
            real(wp)            :: d1ch
            real(wp)            :: d3ch
            real(wp)            :: r0ch
            real(wp)            :: tau
            real(wp)            :: aphi
            real(wp)            :: bphi
            real(wp)            :: cphi
            real(wp)            :: atheta
            real(wp)            :: btheta
            real(wp)            :: ctheta
            real(wp)            :: a3s
            real(wp)            :: b3s
            real(wp)            :: fd0
            real(wp)            :: hd0
            real(wp)            :: a1s
            real(wp)            :: a2s
            real(wp)            :: b1s
            real(wp)            :: b2s
            real(wp)            :: ak
            real(wp)            :: kch3
            real(wp)            :: aa1
            real(wp)            :: aa2
            real(wp)            :: aa3
            real(wp)            :: aa4
            real(wp), parameter :: fact1 = 0.041840e0_wp
            real(wp), parameter :: fact2 = 6.022045e0_wp
            ! Funct called
            ! Local variables
            !     Local parameters
            integer, parameter  :: iout       = 6
            logical, parameter  :: ldebug     = .false.
            real(wp), parameter :: kj_in_kcal = 0.2390057
            !     bond lenght and bond vectors
            real(wp) :: rch2
            real(wp) :: rch3
            real(wp) :: rch4
            real(wp) :: rch5
            real(wp) :: rbh2
            real(wp) :: rbh3
            real(wp) :: rbh4
            real(wp) :: rbh5
            real(wp) :: rcb
            real(wp) :: rchav
            real(wp) :: bvch2(3)
            real(wp) :: bvch3(3)
            real(wp) :: bvch4(3)
            real(wp) :: bvch5(3)
            real(wp) :: bvcb(3)
            !     Vstr varaibles
            real(wp) :: cbq
            real(wp) :: cbj
            real(wp) :: bh2q
            real(wp) :: bh2j
            real(wp) :: bh3q
            real(wp) :: bh3j
            real(wp) :: bh4q
            real(wp) :: bh4f
            real(wp) :: bh5q
            real(wp) :: bh5j
            real(wp) :: ch2q
            real(wp) :: ch2j
            real(wp) :: ch3q
            real(wp) :: ch3j
            real(wp) :: ch4q
            real(wp) :: ch4j
            real(wp) :: ch5q
            real(wp) :: ch5j
            real(wp) :: vqh2
            real(wp) :: vqh3
            real(wp) :: vqh4
            real(wp) :: vqh5
            real(wp) :: vjh2
            real(wp) :: vjh3
            real(wp) :: vjh4
            real(wp) :: vjh5
            real(wp) :: v3h2
            real(wp) :: v3h3
            real(wp) :: v3h4
            real(wp) :: v3h5
            real(wp) :: vstr
            !     Bond xangles and out-of-plane-bendings
            real(wp) :: pi
            real(wp) :: theta023
            real(wp) :: theta024
            real(wp) :: theta025
            real(wp) :: theta032
            real(wp) :: theta034
            real(wp) :: theta035
            real(wp) :: theta042
            real(wp) :: theta043
            real(wp) :: theta045
            real(wp) :: theta052
            real(wp) :: theta053
            real(wp) :: theta054
            real(wp) :: delta23
            real(wp) :: delta24
            real(wp) :: delta25
            real(wp) :: delta32
            real(wp) :: delta34
            real(wp) :: delta35
            real(wp) :: delta42
            real(wp) :: delta43
            real(wp) :: delta45
            real(wp) :: delta52
            real(wp) :: delta53
            real(wp) :: delta54
            real(wp) :: d23
            real(wp) :: d24
            real(wp) :: d25
            real(wp) :: d32
            real(wp) :: d34
            real(wp) :: d35
            real(wp) :: d42
            real(wp) :: d43
            real(wp) :: d45
            real(wp) :: d52
            real(wp) :: d53
            real(wp) :: d54
            real(wp) :: fd2
            real(wp) :: fd3
            real(wp) :: fd4
            real(wp) :: fd5
            real(wp) :: hd2
            real(wp) :: hd3
            real(wp) :: hd4
            real(wp) :: hd5
            real(wp) :: voph2
            real(wp) :: voph3
            real(wp) :: voph4
            real(wp) :: voph5
            real(wp) :: vop
            !     in-plane bending parameters
            real(wp) :: xk023
            real(wp) :: xk024
            real(wp) :: xk025
            real(wp) :: xk034
            real(wp) :: xk035
            real(wp) :: xk045
            real(wp) :: xk2
            real(wp) :: xk3
            real(wp) :: xk4
            real(wp) :: xk5
            real(wp) :: theta23
            real(wp) :: theta24
            real(wp) :: theta25
            real(wp) :: theta34
            real(wp) :: theta35
            real(wp) :: theta45
            real(wp) :: vip23
            real(wp) :: vip24
            real(wp) :: vip25
            real(wp) :: vip34
            real(wp) :: vip35
            real(wp) :: vip45
            real(wp) :: vip
            !     misc and temp variables
            integer  :: i, j, k
            real(wp) :: t1
            real(wp) :: t2
            real(wp) :: t3
            real(wp) :: t4
            real(wp) :: t5
            real(wp) :: t6
            real(wp) :: z1
            real(wp) :: z2
            real(wp) :: z3
            real(wp) :: z4
            real(wp) :: scat0
            real(wp) :: xpot
            !     CH5 surf
            real(wp) :: vsurf
            real(wp) :: r(18)
            real(wp) :: dv(18)
            ! Assign parameters
            a1ch  = para(1)
            b1ch  = para(2)
            c1ch  = para(3)
            alpcb = para(4)
            d1cb  = para(5)
            d3cb  = para(6)
            r0cb  = para(7)
            alphh = para(8)
            d1hh  = para(9)
            d3hh  = para(10)
            r0hh  = para(11)
            d1ch  = para(12)
            d3ch  = para(13)
            r0ch  = para(14)
            ! out-of-plane
            aphi   = para(15)
            bphi   = para(16)
            cphi   = para(17)
            atheta = para(18)
            btheta = para(19)
            ctheta = para(20)
            a3s    = para(21)
            b3s    = para(22)
            fd0    = para(23)
            hd0    = para(24)
            ! in-plane bending parameters
            a1s  = para(25)
            a2s  = para(26)
            b1s  = para(27)
            b2s  = para(28)
            ak   = para(29)
            kch3 = para(30)
            aa1  = para(31)
            aa2  = para(32)
            aa3  = para(33)
            aa4  = para(34)
            ! Define j,k,l indices wrt to each of the 4 H
            ! H2 : 3,4,5
            ! H3 : 5,4,2
            ! H4 : 5,2,3
            ! H5 : 2,4,3
            ! Calculate bond lenghts and bond vectors
            ! bv is the bond vectors and rch is the
            ! corresponding bond lenghts
            ! rchav is the avg of the four ch xdist
            rch2 = 0.0e0_wp
            rch3 = 0.0e0_wp
            rch4 = 0.0e0_wp
            rch5 = 0.0e0_wp
            rbh2 = 0.0e0_wp
            rbh3 = 0.0e0_wp
            rbh4 = 0.0e0_wp
            rbh5 = 0.0e0_wp
            rcb  = 0.0e0_wp
            do i=1, 3
                bvch2(i) = h2(i)-c1(i)
                bvch3(i) = h3(i)-c1(i)
                bvch4(i) = h4(i)-c1(i)
                bvch5(i) = h5(i)-c1(i)
                bvcb(i)  = hb(i)-c1(i)
                rch2     = rch2+(bvch2(i)*bvch2(i))
                rch3     = rch3+(bvch3(i)*bvch3(i))
                rch4     = rch4+(bvch4(i)*bvch4(i))
                rch5     = rch5+(bvch5(i)*bvch5(i))
                rcb      = rcb +( bvcb(i)*bvcb(i) )
                rbh2     = rbh2+((hb(i)-h2(i))*(hb(i)-h2(i)))
                rbh3     = rbh3+((hb(i)-h3(i))*(hb(i)-h3(i)))
                rbh4     = rbh4+((hb(i)-h4(i))*(hb(i)-h4(i)))
                rbh5     = rbh5+((hb(i)-h5(i))*(hb(i)-h5(i)))
            enddo
            rch2  = sqrt(rch2)
            rch3  = sqrt(rch3)
            rch4  = sqrt(rch4)
            rch5  = sqrt(rch5)
            rcb   = sqrt(rcb)
            rbh2  = sqrt(rbh2)
            rbh3  = sqrt(rbh3)
            rbh4  = sqrt(rbh4)
            rbh5  = sqrt(rbh5)
            rchav = 0.250e0_wp*(rch2+rch3+rch4+rch5)
            ! use switch func to change r0ch from
            ! ethane to ethyl value
            t1 = 1.09250e0_wp
            t2 = 1.08850e0_wp
            t3 = c1ch*10.0e0_wp
            t4 = 1.12750e0_wp
            t5 = 0.50e0_wp*(1.0e0_wp+tanh(t3*(rchav-t4)))
            !     r0ch=(t5*t2)+((1.0d0-t5)*t1)
            ! calculate the Morse parameter alphaCH (alpch)
            ! used in the Vstr using Equation 2 of
            ! Ref. [2]
            t1    = 0.50e0_wp*( 1.0e0_wp+tanh(c1ch*(rchav-r0ch)) )
            alpch = a1ch+(b1ch*t1)
            ! calculate Vstr
            ! first get the coulomb and resonance integrals
            ! involving Hb and all the other atoms
            ! cbq and cbj : coulomb and resonance integral between C1 and Hb
            call v3qj77(alpcb, d1cb, d3cb, r0cb, rcb, cbq, cbj)
            call v3qj77(alphh, d1hh, d3hh, r0hh, rbh2, bh2q, bh2j)
            call v3qj77(alphh, d1hh, d3hh, r0hh, rbh3, bh3q, bh3j)
            call v3qj77(alphh, d1hh, d3hh, r0hh, rbh4, bh4q, bh4f)
            call v3qj77(alphh, d1hh, d3hh, r0hh, rbh5, bh5q, bh5j)
            !%%%--------------------------------------------------
            !     t1 = 0.9950d0
            !     bh2q=bh2q*t1
            !     bh3q=bh3q*t1
            !     bh4q=bh4q*t1
            !     bh5q=bh5q*t1
            !     bh2j=bh2j*t1
            !     bh3j=bh3j*t1
            !     bh4f=bh4f*t1
            !     bh5j=bh5j*t1
            !%%%--------------------------------------------------
            ! C1 and (H2...H5) terms
            t1 = 0.50e0_wp*(1.0e0_wp+tanh(10.0e0_wp*(rbh2-1.50e0_wp)))
            t2 = 1.0e0_wp-t1
            t3 = (1.0e0_wp+(t2*0.080e0_wp))*r0ch
            !     call V3qj77(alpch,d1ch,d3ch,t3,rch2,ch2q,ch2j)
            call v3qj77(alpch, d1ch, d3ch, r0ch, rch2, ch2q, ch2j)
            call v3qj77(alpch, d1ch, d3ch, r0ch, rch3, ch3q, ch3j)
            call v3qj77(alpch, d1ch, d3ch, r0ch, rch4, ch4q, ch4j)
            call v3qj77(alpch, d1ch, d3ch, r0ch, rch5, ch5q, ch5j)
            ! coulomb and resonance integral for H2...H5
            ! Equation 2 of Ref. [1]
            vqh2 = cbq+ch2q+bh2q
            vqh3 = cbq+ch3q+bh3q
            vqh4 = cbq+ch4q+bh4q
            vqh5 = cbq+ch5q+bh5q
            t2   = ((ch2j-cbj)**2) + ((cbj-bh2j )**2) + ((bh2j-ch2j)**2)
            t3   = ((ch3j-cbj)**2) + ((cbj-bh3j )**2) + ((bh3j-ch3j)**2)
            t4   = ((ch4j-cbj)**2) + ((cbj-bh4f )**2) + ((bh4f-ch4j)**2)
            t5   = ((ch5j-cbj)**2) + ((cbj-bh5j )**2) + ((bh5j-ch5j)**2)
            vjh2 = -sqrt(0.50e0_wp*t2)
            vjh3 = -sqrt(0.50e0_wp*t3)
            vjh4 = -sqrt(0.50e0_wp*t4)
            vjh5 = -sqrt(0.50e0_wp*t5)
            v3h2 = vqh2+vjh2
            v3h3 = vqh3+vjh3
            v3h4 = vqh4+vjh4
            v3h5 = vqh5+vjh5
            vstr = v3h2+v3h3+v3h4+v3h5
            !%%%--------------------------------------------------
            !     t1=10.0d0
            !     t2=2.144106410736830d0
            !     t3=0.013850d0
            !     t4=1.0d0+(t3*dexp(-t1*((rxncrd-t2)**2)))
            t1 = 6.0e0_wp
            t2 = 0.950e0_wp
            t3 = 0.50e0_wp*(1.0e0_wp+tanh(t1*(rxncrd-t2)))
            t4 = 1.013850e0_wp
            call xgaussfit(rxncrd, t4)
            !     t4=t4*0.999930d0
            !     t4=t4*1.00205d0
            !     t4=t4*1.0159995120d0/1.016043066097427340d0
            t4   = t4*1.016043066097427340e0_wp/1.0159995120e0_wp
            t5   = (1.0e0_wp*(1.0e0_wp-t3))+(t3*t4)
            z1   = 6.0e0_wp
            z2   = 10.0e0_wp
            z3   = 0.50e0_wp*(1.0e0_wp+tanh(z1*(rxncrd-z2)))
            z4   = (t5*(1.0e0_wp-z3))+(1.0e0_wp*z3)
            vstr = vstr * z4
            !     vstr = vstr * 1.013850d0
            !%%%--------------------------------------------------
            if( ldebug )then
                write(*, *)
                write(*, *)'H2 Vstr :', v3h2*1.0e2_wp*kj_in_kcal
                write(*, *)'H3 Vstr :', v3h3*1.0e2_wp*kj_in_kcal
                write(*, *)'H4 Vstr :', v3h4*1.0e2_wp*kj_in_kcal
                write(*, *)'H5 Vstr :', v3h5*1.0e2_wp*kj_in_kcal
                write(*, *) 'VSTR FROM XPOT:', vstr*1.0e2_wp*kj_in_kcal
            endif
            ! calculate the out-of-plane bend terms
            ! evaluate tau and PI
            pi  = acos(-1.0e0_wp)
            tau = acos(-1.0e0_wp/3.0e0_wp)
            ! calculate theta0
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch2, rch3, rch4, rch5, theta023)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch2, rch4, rch5, rch3, theta024)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch2, rch5, rch3, rch4, theta025)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch3, rch2, rch5, rch4, theta032)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch3, rch4, rch2, rch5, theta034)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch3, rch5, rch4, rch2, theta035)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch4, rch2, rch3, rch5, theta042)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch4, rch3, rch5, rch2, theta043)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch4, rch5, rch2, rch3, theta045)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch5, rch2, rch4, rch3, theta052)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch5, rch3, rch2, rch4, theta053)
            call gjtheta0(pi, tau, aphi, bphi, cphi, &
                & atheta, btheta, ctheta, r0ch, &
                & rch5, rch4, rch3, rch2, theta054)
            !     scat0=1.0d0
            !     theta023=theta023*scat0
            !     theta024=theta024*scat0
            !     theta025=theta025*scat0
            !     theta032=theta032*scat0
            !     theta034=theta034*scat0
            !     theta035=theta035*scat0
            !     theta042=theta042*scat0
            !     theta043=theta043*scat0
            !     theta045=theta045*scat0
            !     theta052=theta052*scat0
            !     theta053=theta053*scat0
            !     theta054=theta054*scat0
            ! calculate delta
            call steckdel(bvch3, bvch4, bvch5, bvch3, d23)
            call steckdel(bvch3, bvch4, bvch5, bvch4, d24)
            call steckdel(bvch3, bvch4, bvch5, bvch5, d25)
            call steckdel(bvch5, bvch4, bvch2, bvch2, d32)
            call steckdel(bvch5, bvch4, bvch2, bvch4, d34)
            call steckdel(bvch5, bvch4, bvch2, bvch5, d35)
            call steckdel(bvch5, bvch2, bvch3, bvch2, d42)
            call steckdel(bvch5, bvch2, bvch3, bvch3, d43)
            call steckdel(bvch5, bvch2, bvch3, bvch5, d45)
            call steckdel(bvch2, bvch4, bvch3, bvch2, d52)
            call steckdel(bvch2, bvch4, bvch3, bvch3, d53)
            call steckdel(bvch2, bvch4, bvch3, bvch4, d54)
            delta23 = d23-theta023
            delta24 = d24-theta024
            delta25 = d25-theta025
            delta32 = d32-theta032
            delta34 = d34-theta034
            delta35 = d35-theta035
            delta42 = d42-theta042
            delta43 = d43-theta043
            delta45 = d45-theta045
            delta52 = d52-theta052
            delta53 = d53-theta053
            delta54 = d54-theta054
            if( ldebug )then
                write(iout, *)
                write(iout, *) 'theta0 and delta terms:'
                write(iout, *) 'H2-H3:', theta023, delta23
                write(iout, *) 'H2-H4:', theta024, delta24
                write(iout, *) 'H2-H5:', theta025, delta25
                write(iout, *) 'H3-H4:', theta034, delta34
                write(iout, *) 'H3-H5:', theta035, delta35
                write(iout, *) 'H4-H5:', theta045, delta45
            endif
            ! calculate fd and hd terms
            call calcfh(a3s, b3s, fd0, hd0, r0ch, rch2, rch3, rch4, rch5, fd2, hd2)
            call calcfh(a3s, b3s, fd0, hd0, r0ch, rch3, rch5, rch4, rch2, fd3, hd3)
            call calcfh(a3s, b3s, fd0, hd0, r0ch, rch4, rch5, rch2, rch3, fd4, hd4)
            call calcfh(a3s, b3s, fd0, hd0, r0ch, rch5, rch2, rch4, rch3, fd5, hd5)
            if( ldebug )then
                write(iout, *)
                write(iout, *)'FD and HD terms:'
                write(iout, *)'H2', fd2, hd2
                write(iout, *)'H3', fd3, hd3
                write(iout, *)'H4', fd4, hd4
                write(iout, *)'H5', fd5, hd5
            endif
            ! calculate vop contributions for each H
            t1    = delta23*delta23
            t2    = delta24*delta24
            t3    = delta25*delta25
            voph2 = (fd2*(t1+t2+t3))+ ( hd2*((t1*t1)+(t2*t2)+(t3*t3)) )
            t1    = delta32*delta32
            t2    = delta34*delta34
            t3    = delta35*delta35
            voph3 = (fd3*(t1+t2+t3))+ ( hd3*((t1*t1)+(t2*t2)+(t3*t3)) )
            t1    = delta42*delta42
            t2    = delta43*delta43
            t3    = delta45*delta45
            voph4 = (fd4*(t1+t2+t3))+ ( hd4*((t1*t1)+(t2*t2)+(t3*t3)) )
            t1    = delta52*delta52
            t2    = delta53*delta53
            t3    = delta54*delta54
            voph5 = (fd5*(t1+t2+t3))+ ( hd5*((t1*t1)+(t2*t2)+(t3*t3)) )
            ! total vop contribution
            vop = voph2+voph3+voph4+voph5
            if( ldebug )then
                write(*, *)
                write(*, *)'H2 Vop : ', voph2*1.0e2_wp*kj_in_kcal
                write(*, *)'H3 Vop : ', voph3*1.0e2_wp*kj_in_kcal
                write(*, *)'H4 Vop : ', voph4*1.0e2_wp*kj_in_kcal
                write(*, *)'H5 Vop : ', voph5*1.0e2_wp*kj_in_kcal
                write(iout, *) 'VOP FROM XPOT :', vop*1.0e2_wp*kj_in_kcal
            endif
            ! Calculate in-plane terms
            ! get the k0ij terms
            call calcvipk0(a1s, b1s, a2s, b2s, ak, kch3, r0ch, &
                & rch2, rch3, rch4, rch5, xk023)
            call calcvipk0(a1s, b1s, a2s, b2s, ak, kch3, r0ch, &
                & rch2, rch4, rch5, rch3, xk024)
            call calcvipk0(a1s, b1s, a2s, b2s, ak, kch3, r0ch, &
                & rch2, rch5, rch3, rch4, xk025)
            !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
            !    1               rch3,rch2,rch5,rch4,xk032)
            call calcvipk0(a1s, b1s, a2s, b2s, ak, kch3, r0ch, &
                & rch3, rch4, rch2, rch5, xk034)
            call calcvipk0(a1s, b1s, a2s, b2s, ak, kch3, r0ch, &
                & rch3, rch5, rch4, rch2, xk035)
            !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
            !    1               rch4,rch2,rch3,rch5,xk042)
            !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
            !    1               rch4,rch3,rch5,rch2,xk043)
            call calcvipk0(a1s, b1s, a2s, b2s, ak, kch3, r0ch, &
                & rch4, rch5, rch2, rch3, xk045)
            !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,rch5,rch2,rch4,rch3,xk
            !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,rch5,rch3,rch2,rch4,xk
            !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,rch5,rch4,rch3,rch2,xk
            ! get ki terms
            call calcvipk(aa1, aa2, aa3, aa4, r0ch, r0hh, rch2, rbh2, xk2)
            call calcvipk(aa1, aa2, aa3, aa4, r0ch, r0hh, rch3, rbh3, xk3)
            call calcvipk(aa1, aa2, aa3, aa4, r0ch, r0hh, rch4, rbh4, xk4)
            call calcvipk(aa1, aa2, aa3, aa4, r0ch, r0hh, rch5, rbh5, xk5)
            ! get bond xangles
            theta23 = xangle(h2, c1, h3)
            theta24 = xangle(h2, c1, h4)
            theta25 = xangle(h2, c1, h5)
            theta34 = xangle(h3, c1, h4)
            theta35 = xangle(h3, c1, h5)
            theta45 = xangle(h4, c1, h5)
            vip23   = xk023*xk2*xk3*(theta23-theta023)*(theta23-theta023)
            vip24   = xk024*xk2*xk4*(theta24-theta024)*(theta24-theta024)
            vip25   = xk025*xk2*xk5*(theta25-theta025)*(theta25-theta025)
            vip34   = xk034*xk3*xk4*(theta34-theta034)*(theta34-theta034)
            vip35   = xk035*xk3*xk5*(theta35-theta035)*(theta35-theta035)
            vip45   = xk045*xk4*xk5*(theta45-theta045)*(theta45-theta045)
            vip     = 0.50e0_wp*(vip23+vip24+vip25+vip34+vip35+vip45)
            ! remove H5 str term
            ! total potential in 10(5)Joules/mol and kcal/mol
            xpot      = vstr+vop+vip-v3h5-voph5-((vip25+vip35+vip45)*0.50e0_wp)
            xpot_kcal = xpot*100.0e0_wp*kj_in_kcal
            if( ldebug )then
                write(*, *)
                write(*, *)'H2-H3 Vip :', 0.50e0_wp*vip23*1.0e2_wp*kj_in_kcal
                write(*, *)'H2-H4 Vip :', 0.50e0_wp*vip24*1.0e2_wp*kj_in_kcal
                write(*, *)'H2-H5 Vip :', 0.50e0_wp*vip25*1.0e2_wp*kj_in_kcal
                write(*, *)'H3-H4 Vip :', 0.50e0_wp*vip34*1.0e2_wp*kj_in_kcal
                write(*, *)'H3-H5 Vip :', 0.50e0_wp*vip35*1.0e2_wp*kj_in_kcal
                write(*, *)'H4-H5 Vip :', 0.50e0_wp*vip45*1.0e2_wp*kj_in_kcal
                write(*, *) 'VIP FROM XPOT:', vip*1.0e2_wp*kj_in_kcal
                write(*, *) 'VSTR+VIP+VOP :', (vstr+vip+vop)*1.0e2_wp*kj_in_kcal
                write(*, *)
                write(*, *)'Energy after removing H5 term'
                write(*, *)'Vstr:', (vstr-v3h5)*1.0e2_wp*kj_in_kcal
                write(*, *)'Vop :', (vop-voph5)*1.0e2_wp*kj_in_kcal
                write(*, *)'Vip :', vip-((vip25+vip35+vip45)*0.50e0_wp)* &
                    & 1.0e2_wp*kj_in_kcal
                write(*, *)'Total xpot:', xpot_kcal
            endif
            if( .false. )then
                r(1)  = h2(1)/0.52918e0_wp
                r(2)  = h2(2)/0.52918e0_wp
                r(3)  = h2(3)/0.52918e0_wp
                r(4)  = c1(1)/0.52918e0_wp
                r(5)  = c1(2)/0.52918e0_wp
                r(6)  = c1(3)/0.52918e0_wp
                r(7)  = h3(1)/0.52918e0_wp
                r(8)  = h3(2)/0.52918e0_wp
                r(9)  = h3(3)/0.52918e0_wp
                r(10) = h4(1)/0.52918e0_wp
                r(11) = h4(2)/0.52918e0_wp
                r(12) = h4(3)/0.52918e0_wp
                r(13) = h5(1)/0.52918e0_wp
                r(14) = h5(2)/0.52918e0_wp
                r(15) = h5(3)/0.52918e0_wp
                r(16) = hb(1)/0.52918e0_wp
                r(17) = hb(2)/0.52918e0_wp
                r(18) = hb(3)/0.52918e0_wp
                call csurf(vsurf, r, dv, 18)
                !       vsurf=(vsurf-v3h5)*100.0d0*KJ_IN_KCAL
                write(*, *) 'VSURF:', vsurf*627.510e0_wp
                write(*, *) 'MYVSURF:', (vstr+vop+vip)*0.03812e0_wp*627.510e0_wp
            endif
        end subroutine xch5pot

        !> @brief Evaluate V3qj77 for the C2H7_VBMM potential.
        subroutine v3qj77( alpha, d1, d3, r0, r, vq, vj )
            ! this subroutine calculates the Coulomb (vq) and
            ! resonance (vj) integrals from singlet (e1) and
            ! triplet (E3) states
            ! Equations: 2, 3, 4, 5, and 8 of Ref. [1]
            use chempotpy_kinds, only: wp
            implicit none
            ! Input variables
            real(wp) :: alpha
            real(wp) :: d1
            real(wp) :: d3
            real(wp) :: r0
            real(wp) :: r
            ! Output variables
            real(wp) :: vq
            real(wp) :: vj
            ! Local variables
            real(wp) :: z1, z2
            real(wp) :: e1, e3
            z1 = exp(-alpha*(r-r0))
            z2 = exp(-2.0e0_wp*alpha*(r-r0))
            e1 = d1*(z2-(2.0e0_wp*z1))
            e3 = d3*(z2+(2.0e0_wp*z1))
            vq = 0.50e0_wp*(e1+e3)
            vj = 0.50e0_wp*(e1-e3)
        end subroutine v3qj77

        !> @brief Evaluate funcS3 for the C2H7_VBMM potential.
        function funcs3( alpha, beta, r0, r )
            ! This func. calculates the S3 switching func
            ! Equation 4 of Ref. [2]
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: alpha
            real(wp) :: beta
            real(wp) :: r0
            real(wp) :: r
            real(wp) :: funcs3
            funcs3 = 1.0e0_wp-tanh(alpha*(r-r0)*(r-beta)*(r-beta))
        end function funcs3

        !> @brief Evaluate calcfh for the C2H7_VBMM potential.
        subroutine calcfh( a3, b3, fd0, hd0, r0, ri, rj, rk, rl, fdi, hdi )
            ! This sub. calculates the fdelta_i and hdelta_i
            ! as defined by Equation 4 of Ref. [2]
            use chempotpy_kinds, only: wp
            implicit none
            ! Input variables
            real(wp) :: a3
            real(wp) :: b3
            real(wp) :: fd0
            real(wp) :: hd0
            real(wp) :: r0
            real(wp) :: ri
            real(wp) :: rj
            real(wp) :: rk
            real(wp) :: rl
            ! Output variables
            real(wp) :: fdi
            real(wp) :: hdi
            ! Local variables
            real(wp) :: s3i
            real(wp) :: s3j
            real(wp) :: s3k
            real(wp) :: s3l
            s3i = 1.0e0_wp-tanh( a3*(ri-r0)*(ri-b3)*(ri-b3) )
            s3j = 1.0e0_wp-tanh( a3*(rj-r0)*(rj-b3)*(rj-b3) )
            s3k = 1.0e0_wp-tanh( a3*(rk-r0)*(rk-b3)*(rk-b3) )
            s3l = 1.0e0_wp-tanh( a3*(rl-r0)*(rl-b3)*(rl-b3) )
            fdi = (1.0e0_wp-s3i)*s3j*s3k*s3l*fd0
            hdi = (1.0e0_wp-s3i)*s3j*s3k*s3l*hd0
        end subroutine calcfh

        !> @brief Evaluate gjtheta0 for the C2H7_VBMM potential.
        subroutine gjtheta0( pi, tau, aphi, bphi, cphi, atheta, btheta, ctheta, &
                & r0, ri, rj, rk, rl, theta0ij )
            ! This sub. computes the theta0 for bending terms
            ! Equation 5 in Ref. [2]
            use chempotpy_kinds, only: wp
            implicit none
            ! Input variables
            real(wp) :: pi
            real(wp) :: tau
            real(wp) :: aphi
            real(wp) :: bphi
            real(wp) :: cphi
            real(wp) :: atheta
            real(wp) :: btheta
            real(wp) :: ctheta
            real(wp) :: r0
            real(wp) :: ri
            real(wp) :: rj
            real(wp) :: rk
            real(wp) :: rl
            ! Output variables
            real(wp) :: theta0ij
            ! Local variables
            real(wp) :: si
            real(wp) :: sj
            real(wp) :: sk
            real(wp) :: sl
            real(wp) :: s1
            real(wp) :: s2
            real(wp) :: t1
            real(wp) :: t2
            real(wp) :: sphii
            real(wp) :: sphij
            real(wp) :: sthetak
            real(wp) :: sthetal
            si       = aphi*(ri-r0)*exp(bphi*(ri-cphi)*(ri-cphi)*(ri-cphi))
            sj       = aphi*(rj-r0)*exp(bphi*(rj-cphi)*(rj-cphi)*(rj-cphi))
            sk       = atheta*(rk-r0)*exp(btheta*(rk-ctheta)*(rk-ctheta)*(rk-ctheta))
            sl       = atheta*(rl-r0)*exp(btheta*(rl-ctheta)*(rl-ctheta)*(rl-ctheta))
            sphii    = 1.0e0_wp-tanh(si)
            sphij    = 1.0e0_wp-tanh(sj)
            sthetak  = 1.0e0_wp-tanh(sk)
            sthetal  = 1.0e0_wp-tanh(sl)
            s1       = (sphii*sphij)-1.0e0_wp
            s2       = (sthetak*sthetal)-1.0e0_wp
            t1       = tau-(0.50e0_wp*pi)
            t2       = tau-(2.0e0_wp*pi/3.0e0_wp)
            theta0ij = tau+(t1*s1)+(t2*s2)
        end subroutine gjtheta0

        !> @brief Evaluate steckdel for the C2H7_VBMM potential.
        subroutine steckdel( v1, v2, v3, vi, deltai )
            ! This sub. computes the out-of-plane delta variable
            ! as defined by Steckler et al. in Ref. [3] equation 18.
            ! Note that is *NOT* symmetric with respect to all
            ! the three methyl hydrogen.
            ! Symmetrization (as given in Ref. [2] equation 4)
            ! will be done in the calling suroutine.
            use chempotpy_kinds, only: wp
            implicit none
            ! Input variables
            real(wp) :: v1(3)
            real(wp) :: v2(3)
            real(wp) :: v3(3)
            real(wp) :: vi(3)
            ! Output variables
            real(wp) :: deltai
            ! Local variables
            integer  :: i
            real(wp) :: ri
            real(wp) :: rn
            real(wp) :: r21
            real(wp) :: r31
            real(wp) :: doti
            real(wp) :: v21(3)
            real(wp) :: v31(3)
            real(wp) :: pn(3)
            ! r21, r31 and ri are the lenghts of vectors
            ! v21, v31 and vi resp.
            ! v21=v2-v1
            ! v31=v3-v1
            do i=1, 3
                v21(i) = v2(i)-v1(i)
                v31(i) = v3(i)-v1(i)
            enddo
            ! the normal to the plane containing
            ! points 1, 2, and 3 is given by
            ! pn = r21 X r31
            ! rn is lenght of pn
            call crossprod(v21, v31, pn)
            ! The xangle between pn and vi is deltaI
            ! evalute the dot product:
            ! doti = pn . vi
            ri   = 0.0e0_wp
            rn   = 0.0e0_wp
            doti = 0.0e0_wp
            do i=1, 3
                ri   = ri+(vi(i)*vi(i))
                rn   = rn+(pn(i)*pn(i))
                doti = doti+(pn(i)*vi(i))
            enddo
            ! doti should be -ve by convension
            ! for right-handed set
            doti   = -abs(doti)
            ri     = sqrt(ri)
            rn     = sqrt(rn)
            deltai = acos(doti/(ri*rn))
        end subroutine steckdel

        !> @brief Evaluate calcvipk for the C2H7_VBMM potential.
        subroutine calcvipk( a1, a2, a3, a4, r0ch, r0bh, rch, rbh, vipk )
            ! The sub. calculates the force constant (vipk) for the in-plane
            ! bendinging as given by equation 6 of Ref. [2]
            use chempotpy_kinds, only: wp
            implicit none
            ! Input varaibles
            real(wp) :: a1
            real(wp) :: a2
            real(wp) :: a3
            real(wp) :: a4
            real(wp) :: r0ch
            real(wp) :: r0bh
            real(wp) :: rch
            real(wp) :: rbh
            ! Output variables
            real(wp) :: vipk
            ! Local variables
            real(wp) :: ua1
            real(wp) :: ua2
            ua1  = 1.0e0_wp-exp(-a1*rbh*rbh)
            ua2  = a2+(a3*exp(-a4*(rbh-r0bh)*(rbh-r0bh)))
            vipk = ua1*exp(-ua2*(rch-r0ch)*(rch-r0ch))
        end subroutine calcvipk

        !> @brief Evaluate calcvipk0 for the C2H7_VBMM potential.
        subroutine calcvipk0( a1, b1, a2, b2, ak, kch3, r0ch, ri, rj, rk, rl, vipk0ij )
            ! This sub. calculates the equilibrium force constant (k0ij) used
            ! for in-plane-bending terms
            ! Equation 9 of Ref. [2]
            use chempotpy_kinds, only: wp
            implicit none
            ! Input varaibles
            real(wp) :: a1
            real(wp) :: b1
            real(wp) :: a2
            real(wp) :: b2
            real(wp) :: kch3
            real(wp) :: ak
            real(wp) :: r0ch
            real(wp) :: ri
            real(wp) :: rj
            real(wp) :: rk
            real(wp) :: rl
            ! Output varaibles
            real(wp) :: vipk0ij
            ! Local varaibles
            real(wp) :: kch4
            real(wp) :: s1i
            real(wp) :: s1j
            real(wp) :: s2k
            real(wp) :: s2l
            kch4 = kch3+ak
            ! calculate switching func
            s1i = 1.0e0_wp-tanh( a1*(ri-r0ch)*((ri-b1)**8) )
            s1j = 1.0e0_wp-tanh( a1*(rj-r0ch)*((rj-b1)**8) )
            s2k = 1.0e0_wp-tanh( a2*(rk-r0ch)*((rk-b2)**6) )
            s2l = 1.0e0_wp-tanh( a2*(rl-r0ch)*((rl-b2)**6) )
            ! calculate k0 using Eq. 9 of Ref. [2]
            vipk0ij = kch4+( kch4*((s1i*s1j)-1.0e0_wp) ) +( (kch4-kch3)*((s2k*s2l)-1.0e0_wp) )
        end subroutine calcvipk0

        !> @brief Evaluate xlimit for the C2H7_VBMM potential.
        subroutine xlimit( a, b, s, x0, x, y )
            ! This subroutine uses a tanh function
            ! to smoothly connect two points (a,b)
            ! The tanh function is centered at x0
            ! with a scale factor of s.
            ! The output value is given by y
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a
            real(wp) :: b
            real(wp) :: s
            real(wp) :: x0
            real(wp) :: x
            real(wp) :: y
            real(wp) :: r
            real(wp) :: ca
            real(wp) :: cb
            ! When s is +ve:
            ! As x-->b
            ! cb--->1 & ca-->0
            ! As x-->a
            ! ca--->1 & cb-->0
            r  = x-x0
            cb = 0.50e0_wp*(1.0e0_wp+tanh(s*r))
            ca = 1.0e0_wp-cb
            y  = (ca*a)+(cb*b)
        end subroutine xlimit

        !> @brief Evaluate CSURF for the C2H7_VBMM potential.
        subroutine csurf( engygs, r, degsdr, n3tm )
            !        This potential is written such that:
            !                       X(1)  - X(3)  : X, Y, Z for H1
            !                       X(4)  - X(6)  : X, Y, Z for C
            !                       X(7)  - X(9)  : X, Y, Z for H3
            !                       X(10) - X(12) : X, Y, Z for H4
            !                       X(13) - X(15) : X, Y, Z for H2
            !                       X(16) - X(18) : X, Y, Z for H5
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: degsdr(n3tm), en, engygs, r(n3tm), vip, vop, vstr
            integer             :: i, n3tm
            real(wp), parameter :: pi = 3.141592653589793e0_wp
            call prepare_local_context
            ! Changing units from angstrom to bohr and initialize
            do i = 1, 18
                qpdot_pl_q(i)    = r(i)*0.52918e0_wp
                qpdot_pl_pdot(i) = 0.e0_wp
            enddo
            !  calculate relative coordinates and bond lengths
            en = 0.0e0_wp
            call coorden
            !  calculate switching functions
            call switchf
            !  calculate reference xangles and their derivatives
            call refxangles
            !  calculate stretching potential
            call stretch(vstr)
            !  calculate out of plane bending potential
            call opbend(vop)
            !  calculate in plane bending potential
            call ipbend(vip)
            !  total potential energy is vstr+vop+vip
            en = vstr+vop+vip
            ! 0.03812 conversion factor from 10(5) j/mol to hartrees
            !##--------------------------------------------------
            !       write(*,*)
            !       write(*,*) 'vstr 10(5)j/mol and au : ',vstr,vstr*0.03812D0
            !       write(*,*) 'vop  10(5)j/mol and au : ',vop,vop*0.03812D0
            !       write(*,*) 'vip  10(5)j/mol and au : ',vip,vip*0.03812D0
            !       write(*,*) 'pot  10(5)j/mol and au : ',en,en*0.03812D0
            !##--------------------------------------------------
            en     = en*0.03812e0_wp
            engygs = en
            ! 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
            do i=1, 18
                degsdr(i) = qpdot_pl_pdot(i)*0.0201723e0_wp
            enddo
        end subroutine csurf

        !> @brief Evaluate coorden for the C2H7_VBMM potential.
        subroutine coorden
            !  calculates relative coordinates and bond lengths
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i, ind
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !  calculate relative coordinates
            do ind=1, 3
                coords_tcb(ind) = qpdot_pl_q(ndx_nc(ind))-qpdot_pl_q(ndx_nhb(ind))
                do i=1, 4
                    coords_tch(i,ind) = qpdot_pl_q(ndx_nc(ind))-qpdot_pl_q(ndx_nh(i,ind))
                    coords_tbh(i,ind) = qpdot_pl_q(ndx_nhb(ind))-qpdot_pl_q(ndx_nh(i,ind))
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

        !> @brief Evaluate refxangles for the C2H7_VBMM potential.
        subroutine refxangles
            !  subroutine calculates reference xangles for the "in-plane" potential
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: halfpi, tau, twopi
            integer             :: i, j, k
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            tau = acos(-1.0e0_wp/3.0e0_wp)
            !       pi=4.0d0*atan(1.0d0)
            halfpi = 0.5e0_wp*pi
            twopi  = 2.0e0_wp*pi
            !  set diagonal elements to zero
            do i=1, 4
                xangles_theta0(i,i) = 0.0e0_wp
                do k=1, 4
                    xangles_dtheta0(i,i,k) = 0.0e0_wp
                enddo
            enddo
            !  calculate reference xangles
            !##
            xangles_theta0(1,2) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(2)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(3)*switch1_stheta(4)-1.0e0_wp)
            xangles_theta0(1,3) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(3)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(2)*switch1_stheta(4)-1.0e0_wp)
            xangles_theta0(1,4) = tau+(tau-halfpi)*(switch1_sphi(1)*switch1_sphi(4)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(2)*switch1_stheta(3)-1.0e0_wp)
            xangles_theta0(2,3) = tau+(tau-halfpi)*(switch1_sphi(2)*switch1_sphi(3)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(1)*switch1_stheta(4)-1.0e0_wp)
            xangles_theta0(2,4) = tau+(tau-halfpi)*(switch1_sphi(2)*switch1_sphi(4)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(1)*switch1_stheta(3)-1.0e0_wp)
            xangles_theta0(3,4) = tau+(tau-halfpi)*(switch1_sphi(3)*switch1_sphi(4)-1.0e0_wp) +(tau-twopi/3.0e0_wp)* &
                & (switch1_stheta(1)*switch1_stheta(2)-1.0e0_wp)
            !  calculate the derivatives of theta0(i,j) in terms of rch(k)
            !  quantity calulated is dtheta0(i,j,k)
            !  derivatives wrt rch(1)
            xangles_dtheta0(1,2,1) = (tau-halfpi)*switch1_dsphi(1)*switch1_sphi(2)
            xangles_dtheta0(1,3,1) = (tau-halfpi)*switch1_dsphi(1)*switch1_sphi(3)
            xangles_dtheta0(1,4,1) = (tau-halfpi)*switch1_dsphi(1)*switch1_sphi(4)
            xangles_dtheta0(2,3,1) = (tau-twopi/3.0e0_wp)*switch1_dstheta(1)*switch1_stheta(4)
            xangles_dtheta0(2,4,1) = (tau-twopi/3.0e0_wp)*switch1_dstheta(1)*switch1_stheta(3)
            xangles_dtheta0(3,4,1) = (tau-twopi/3.0e0_wp)*switch1_dstheta(1)*switch1_stheta(2)
            !  derivatives wrt rch(2)
            xangles_dtheta0(1,2,2) = (tau-halfpi)*switch1_sphi(1)*switch1_dsphi(2)
            xangles_dtheta0(1,3,2) = (tau-twopi/3.0e0_wp)*switch1_dstheta(2)*switch1_stheta(4)
            xangles_dtheta0(1,4,2) = (tau-twopi/3.0e0_wp)*switch1_dstheta(2)*switch1_stheta(3)
            xangles_dtheta0(2,3,2) = (tau-halfpi)*switch1_dsphi(2)*switch1_sphi(3)
            xangles_dtheta0(2,4,2) = (tau-halfpi)*switch1_dsphi(2)*switch1_sphi(4)
            xangles_dtheta0(3,4,2) = (tau-twopi/3.0e0_wp)*switch1_stheta(1)*switch1_dstheta(2)
            !  derivatives wrt rch(3)
            xangles_dtheta0(1,2,3) = (tau-twopi/3.0e0_wp)*switch1_dstheta(3)*switch1_stheta(4)
            xangles_dtheta0(1,3,3) = (tau-halfpi)*switch1_sphi(1)*switch1_dsphi(3)
            xangles_dtheta0(1,4,3) = (tau-twopi/3.0e0_wp)*switch1_stheta(2)*switch1_dstheta(3)
            xangles_dtheta0(2,3,3) = (tau-halfpi)*switch1_sphi(2)*switch1_dsphi(3)
            xangles_dtheta0(2,4,3) = (tau-twopi/3.0e0_wp)*switch1_stheta(1)*switch1_dstheta(3)
            xangles_dtheta0(3,4,3) = (tau-halfpi)*switch1_dsphi(3)*switch1_sphi(4)
            !  derivatives wrt rch(4)
            xangles_dtheta0(1,2,4) = (tau-twopi/3.0e0_wp)*switch1_stheta(3)*switch1_dstheta(4)
            xangles_dtheta0(1,3,4) = (tau-twopi/3.0e0_wp)*switch1_stheta(2)*switch1_dstheta(4)
            xangles_dtheta0(1,4,4) = (tau-halfpi)*switch1_sphi(1)*switch1_dsphi(4)
            xangles_dtheta0(2,3,4) = (tau-twopi/3.0e0_wp)*switch1_stheta(1)*switch1_dstheta(4)
            xangles_dtheta0(2,4,4) = (tau-halfpi)*switch1_sphi(2)*switch1_dsphi(4)
            xangles_dtheta0(3,4,4) = (tau-halfpi)*switch1_sphi(3)*switch1_dsphi(4)
            !  fill in the other half of the matrix
            do i=1, 3
                do j=i+1, 4
                    xangles_theta0(j,i) = xangles_theta0(i,j)
                    do k=1, 4
                        xangles_dtheta0(j,i,k) = xangles_dtheta0(i,j,k)
                    enddo
                enddo
            enddo
        end subroutine refxangles

        !> @brief Evaluate stretch for the C2H7_VBMM potential.
        subroutine stretch( vstr )
            !  subroutine to calculate leps-type stretching potential and its
            !  derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: ach, achdc(3), achdh(4,3), arga, dumach, dumjbh, dumjcb, dumjch, dumjhi, dumjhj, dumqbh, &
                & dumqcb, dumqch, dumqhi, dumqhj, e1, e3, factj, rav, vj(4), vjbh(4), vjcb, vjch(4), vq(4), vqbh(4), vqcb, vqch(4), &
                & vstr
            integer             :: i, ind, j, k
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !  calculate avergage bond length for the methane moiety
            rav = (bonds_rch(1)+bonds_rch(2)+bonds_rch(3)+bonds_rch(4))/4.0e0_wp
            !  initialise:
            vstr = 0.0e0_wp
            !  ach:
            !  in double precision tanh(19.0d0)=1.0d0 and we put the if statement
            !  in to avoid overflow/underflow errors
            arga = potcm_c1ch*(rav-potcm_r0ch)
            if( arga .lt. 19.0e0_wp )then
                ach    = potcm_a1ch+potcm_b1ch*(tanh(arga)+1.0e0_wp)*0.5e0_wp
                dumach = potcm_b1ch*potcm_c1ch/(2.0e0_wp*cosh(arga)**2)
            else
                ach    = potcm_a1ch+potcm_b1ch
                dumach = 0.0e0_wp
            endif
            !  calculate singlet: e1, triplet: e3 energies and vq and vj
            !  terms for each bond
            e1   = potcm_d1cb*(exp(-2.0e0_wp*potcm_acb*(bonds_rcb-potcm_r0cb))-2.0e0_wp*exp(-potcm_acb*(bonds_rcb- &
                & potcm_r0cb)))
            e3   = potcm_d3cb*(exp(-2.0e0_wp*potcm_acb*(bonds_rcb-potcm_r0cb))+2.0e0_wp*exp(-potcm_acb*(bonds_rcb- &
                & potcm_r0cb)))
            vqcb = (e1+e3)*0.5e0_wp
            vjcb = (e1-e3)*0.5e0_wp
            do i=1, 4
                e1      = potcm_d1ch*(exp(-2.0e0_wp*ach*(bonds_rch(i)-potcm_r0ch)) -2.0e0_wp*exp(-ach*(bonds_rch(i)- &
                    & potcm_r0ch)))
                e3      = potcm_d3ch*(exp(-2.0e0_wp*ach*(bonds_rch(i)-potcm_r0ch)) +2.0e0_wp*exp(-ach*(bonds_rch(i)- &
                    & potcm_r0ch)))
                vqch(i) = (e1+e3)*0.5e0_wp
                vjch(i) = (e1-e3)*0.5e0_wp
                e1      = potcm_d1hh*(exp(-2.0e0_wp*potcm_ahh*(bonds_rbh(i)-potcm_r0hh)) -2.0e0_wp*exp(-potcm_ahh* &
                    & (bonds_rbh(i)-potcm_r0hh)))
                e3      = potcm_d3hh*(exp(-2.0e0_wp*potcm_ahh*(bonds_rbh(i)-potcm_r0hh)) +2.0e0_wp*exp(-potcm_ahh* &
                    & (bonds_rbh(i)-potcm_r0hh)))
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
                achdc(ind) = dumach*(coords_tch(1,ind)/bonds_rch(1)+coords_tch(2,ind)/bonds_rch(2) +coords_tch(3,ind)/ &
                    & bonds_rch(3)+coords_tch(4,ind)/bonds_rch(4))/4.0e0_wp
                do i=1, 4
                    achdh(i,ind) = -dumach*coords_tch(i,ind)/bonds_rch(i)/4.0e0_wp
                enddo
            enddo
            dumqcb = -potcm_acb*((potcm_d1cb+potcm_d3cb)*exp(-2.0e0_wp*potcm_acb*(bonds_rcb-potcm_r0cb))- (potcm_d1cb- &
                & potcm_d3cb)*exp(-potcm_acb*(bonds_rcb-potcm_r0cb)))/bonds_rcb
            !  calculate cartesian derivatives:
            !  looping over ch(i) and bh(i)
            do i=1, 4
                dumqbh = -potcm_ahh*((potcm_d1hh+potcm_d3hh)*exp(-2.0e0_wp*potcm_ahh*(bonds_rbh(i)-potcm_r0hh))- &
                    & (potcm_d1hh-potcm_d3hh)*exp(-potcm_ahh*(bonds_rbh(i)-potcm_r0hh)))/bonds_rbh(i)
                factj = 0.5e0_wp/vj(i)
                dumjcb = -potcm_acb*((potcm_d1cb-potcm_d3cb)*exp(-2.0e0_wp*potcm_acb*(bonds_rcb-potcm_r0cb)) - &
                    & (potcm_d1cb+potcm_d3cb)*exp(-potcm_acb*(bonds_rcb-potcm_r0cb)))*factj/bonds_rcb
                dumjbh = -potcm_ahh*((potcm_d1hh-potcm_d3hh)*exp(-2.0e0_wp*potcm_ahh*(bonds_rbh(i)-potcm_r0hh)) - &
                    & (potcm_d1hh+potcm_d3hh)*exp(-potcm_ahh*(bonds_rbh(i)-potcm_r0hh)))*factj/bonds_rbh(i)
                do ind=1, 3
                    !  deriv wrt hb:
                    qpdot_pl_pdot(ndx_nhb(ind))=qpdot_pl_pdot(ndx_nhb(ind)) -coords_tcb(ind)*dumqcb+coords_tbh(i,ind)* &
                        & dumqbh +(vjch(i)-vjcb)*(dumjcb*coords_tcb(ind)) + (vjcb- vjbh(i))*(-dumjcb*coords_tcb(ind)- &
                        & dumjbh*coords_tbh(i,ind)) +(vjbh(i)-vjch(i))*dumjbh*coords_tbh(i,ind)
                    !  dvqch(i)/dc
                    dumqch=-(ach*coords_tch(i,ind)/bonds_rch(i)+achdc(ind)*(bonds_rch(i)-potcm_r0ch)) *((potcm_d1ch+ &
                        & potcm_d3ch)*exp(-2.0e0_wp*ach*(bonds_rch(i)-potcm_r0ch)) - (potcm_d1ch- potcm_d3ch)*exp(- &
                        & ach* (bonds_rch(i)-potcm_r0ch)))
                    qpdot_pl_pdot(ndx_nc(ind)) = qpdot_pl_pdot(ndx_nc(ind))+dumqch+coords_tcb(ind)*dumqcb
                    !  dvqch(i)/dh(i)
                    dumqhi=(ach*coords_tch(i,ind)/bonds_rch(i)-achdh(i,ind)*(bonds_rch(i)-potcm_r0ch)) *((potcm_d1ch+ &
                        & potcm_d3ch)*exp(-2.0e0_wp*ach*(bonds_rch(i)-potcm_r0ch)) - (potcm_d1ch-potcm_d3ch)*exp(-ach* &
                        & (bonds_rch(i)-potcm_r0ch)))
                    qpdot_pl_pdot(ndx_nh(i,ind)) = qpdot_pl_pdot(ndx_nh(i,ind))+dumqhi-coords_tbh(i,ind)*dumqbh
                    !  dvjch(i)/dc
                    dumjch=-(ach*coords_tch(i,ind)/bonds_rch(i)+achdc(ind)*(bonds_rch(i)-potcm_r0ch)) *((potcm_d1ch- &
                        & potcm_d3ch)*exp(-2.0e0_wp*ach*(bonds_rch(i)-potcm_r0ch)) - (potcm_d1ch+ potcm_d3ch)*exp(- &
                        & ach* (bonds_rch(i)-potcm_r0ch)))*factj
                    !  dvj(i)/dnc(ind)
                    qpdot_pl_pdot(ndx_nc(ind))=qpdot_pl_pdot(ndx_nc(ind)) +(vjch(i)-vjcb)*(dumjch-dumjcb* &
                        & coords_tcb(ind)) +(vjcb-vjbh(i))*dumjcb*coords_tcb(ind) - (vjbh(i)-vjch(i))*dumjch
                    !  dvjch(i)/dh(i)
                    dumjhi=(ach*coords_tch(i,ind)/bonds_rch(i)-achdh(i,ind)*(bonds_rch(i)-potcm_r0ch)) *((potcm_d1ch- &
                        & potcm_d3ch)*exp(-2.0e0_wp*ach*(bonds_rch(i)-potcm_r0ch)) - (potcm_d1ch+potcm_d3ch)*exp(-ach* &
                        & (bonds_rch(i)-potcm_r0ch)))*factj
                    !  dvj(i)/dnh(i,ind)
                    qpdot_pl_pdot(ndx_nh(i,ind))=qpdot_pl_pdot(ndx_nh(i,ind)) +(vjch(i)-vjcb)*dumjhi +(vjcb-vjbh(i))* &
                        & dumjbh*coords_tbh(i,ind) +(vjbh(i)- vjch(i))* (-dumjbh*coords_tbh(i,ind)-dumjhi)
                    !  dv(i)/dh(j)
                    do k=1, 3
                        j = i+k
                        if( j .gt. 4 ) j=j-4
                        dumqhj=-achdh(j,ind)*(bonds_rch(i)-potcm_r0ch) *((potcm_d1ch+potcm_d3ch)*exp(-2.0e0_wp*ach* &
                            & (bonds_rch(i)-potcm_r0ch)) -(potcm_d1ch-potcm_d3ch)*exp(-ach* (bonds_rch(i)-potcm_r0ch)))
                        dumjhj=-achdh(j,ind)*(bonds_rch(i)-potcm_r0ch) *((potcm_d1ch-potcm_d3ch)*exp(-2.0e0_wp*ach* &
                            & (bonds_rch(i)-potcm_r0ch)) -(potcm_d1ch+potcm_d3ch)*exp(-ach* (bonds_rch(i)- &
                            & potcm_r0ch)))*factj
                        qpdot_pl_pdot(ndx_nh(j,ind)) = qpdot_pl_pdot(ndx_nh(j,ind))+dumqhj +(vjch(i)-vjcb)*dumjhj - &
                            & (vjbh(i)-vjch(i))*dumjhj
                    enddo
                enddo
            enddo
        end subroutine stretch

        !> @brief Evaluate opbend for the C2H7_VBMM potential.
        subroutine opbend( vop )
            !  subroutine calculates symmetrized vop potential and derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: ddr, sum2, sum4, sumd2(4), sumd4(4), vop
            integer             :: i, ind, itemp, j, k, l
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            vop = 0.0e0_wp
            !  calculate force constants and their derivatives
            call opforce
            !  calculate out-of-plane xangle and derivatives
            do i=1, 4
                j = i+1
                if( j .gt. 4 ) j=j-4
                k = j+1
                if( k .gt. 4 ) k=k-4
                l = k+1
                if( l .gt. 4 ) l=l-4
                !  if i is an even number then we have switched the vectors
                !  from a right handed set to a left handed set
                !  in this case we need to switch vectors k and l around
                if( (i .eq. 2).or.(i .eq. 4) )then
                    itemp = k
                    k     = l
                    l     = itemp
                endif
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
                        qpdot_pl_pdot(ndx_nh(j,ind)) = qpdot_pl_pdot(ndx_nh(j,ind))-coords_tch(j,ind)*ddr/bonds_rch(j)
                        qpdot_pl_pdot(ndx_nc(ind))   = qpdot_pl_pdot(ndx_nc(ind))+coords_tch(j,ind)*ddr/bonds_rch(j)
                    enddo
                enddo
            enddo
        end subroutine opbend

        !> @brief Evaluate ipbend for the C2H7_VBMM potential.
        subroutine ipbend( vip )
            !  subroutine calculates symmetrised in plane bend term
            !  and its derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: costh(4,4), dth(4,4), dth0c, dth0k, dthc, dthi, dthj, termth, theta(4,4), vip
            integer             :: i, ind, j, k
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
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
                    dth(i,j)   = theta(i,j)-xangles_theta0(i,j)
                    vip        = vip+0.5e0_wp*force1_fk0(i,j)*force1_f1(i)*force1_f1(j)*dth(i,j)**2
                    !##
                    !                 write(*,*) i,j,f1(j)
                    !  calculate partial derivatives wrt cartesian coordinates
                    !  calculate pdots wrt theta:
                    termth = -1.0e0_wp/sqrt(1.0e0_wp-costh(i,j)*costh(i,j))
                    do ind=1, 3
                        dthi             = -coords_tch(j,ind)/bonds_rch(i)/bonds_rch(j) +costh(i,j)*coords_tch(i,ind)/ &
                            & bonds_rch(i)/bonds_rch(i)
                        dthi = dthi*termth
                        dthj             = -coords_tch(i,ind)/bonds_rch(i)/bonds_rch(j) +costh(i,j)*coords_tch(j,ind)/ &
                            & bonds_rch(j)/bonds_rch(j)
                        dthj = dthj*termth
                        dthc = -(dthi+dthj)
                        qpdot_pl_pdot(ndx_nh(i,ind)) = qpdot_pl_pdot(ndx_nh(i,ind)) +force1_fk0(i,j)*force1_f1(i)* &
                            & force1_f1(j)*dthi*dth(i,j)
                        qpdot_pl_pdot(ndx_nh(j,ind)) = qpdot_pl_pdot(ndx_nh(j,ind)) +force1_fk0(i,j)*force1_f1(i)* &
                            & force1_f1(j)*dthj*dth(i,j)
                        qpdot_pl_pdot(ndx_nc(ind))    = qpdot_pl_pdot(ndx_nc(ind)) +force1_fk0(i,j)*force1_f1(i)* &
                            & force1_f1(j)*dthc*dth(i,j)
                        do k=1, 4
                            !  calculate pdots wrt force constants and wrt theta0
                            dth0k = -xangles_dtheta0(i,j,k)*coords_tch(k,ind)/bonds_rch(k)
                            dth0c = -dth0k
                            qpdot_pl_pdot(ndx_nh(k,ind))=qpdot_pl_pdot(ndx_nh(k,ind)) -0.5e0_wp*coords_tch(k,ind)* &
                                & force1_dfdc(i,j,k)*dth(i,j)**2/bonds_rch(k) -0.5e0_wp* coords_tbh(k,ind)* &
                                & force1_dfdh(i,j,k)*dth(i,j)**2/bonds_rbh(k) -force1_fk0(i,j)*force1_f1(i)* &
                                & force1_f1(j)*dth0k*dth(i,j)
                            qpdot_pl_pdot(ndx_nc(ind))=qpdot_pl_pdot(ndx_nc(ind)) +0.5e0_wp*coords_tch(k,ind)* &
                                & force1_dfdc(i,j,k)*dth(i,j)**2/bonds_rch(k) -force1_fk0(i,j)*force1_f1(i)* &
                                & force1_f1(j)* dth0c*dth(i,j)
                            qpdot_pl_pdot(ndx_nhb(ind)) = qpdot_pl_pdot(ndx_nhb(ind)) +0.5e0_wp*coords_tbh(k,ind)* &
                                & force1_dfdh(i,j,k)*dth(i,j)**2/bonds_rbh(k)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine ipbend

        !> @brief Evaluate calcdelta for the C2H7_VBMM potential.
        subroutine calcdelta( i, j, k, l, sum2, sum4 )
            !  subroutine calculates out of plane xangle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3), argd(4), atemp1, atemp2(3), atemp3, atemp4, atemp5, axb(3), b(3), c(4,3), &
                & cdot(4, 3, 3), daxb(4, 3, 3), deldot, delta(4), sum2, sum4
            integer             :: i, ii, in(3), ind, j, jj, k, l
            real(wp)            :: norma
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !  initialise
            sum2 = 0.0e0_wp; sum4 = 0.0e0_wp
            !  set j,k,l indices
            in(1) = j
            in(2) = k
            in(3) = l
            !  vector a is rk-rj, vector b is rl-rj
            do ind=1, 3
                a(ind) = qpdot_pl_q(ndx_nh(k,ind))-qpdot_pl_q(ndx_nh(j,ind))
                b(ind) = qpdot_pl_q(ndx_nh(l,ind))-qpdot_pl_q(ndx_nh(j,ind))
            enddo
            !##
            !        write(*,*)
            !        write(*,*) ' i,j,k,l : ',i,j,k,l
            !        write(*,1000) ' i vec :',q(nh(i,1:3))
            !        write(*,1000) ' j vec :',q(nh(j,1:3))
            !        write(*,1000) ' k vec :',q(nh(k,1:3))
            !        write(*,1000) ' l vec :',q(nh(l,1:3))
            !        write(*,*)
            !1000    format(A,3F8.2)
            !  axb is vector a cross b
            axb(1) = a(2)*b(3)-a(3)*b(2)
            axb(2) = a(3)*b(1)-a(1)*b(3)
            axb(3) = a(1)*b(2)-a(2)*b(1)
            norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
            norma  = sqrt(norma)
            !##
            !        write(*,*) 'axb1 : ',i,axb(1)
            !        write(*,*) 'axb2 : ',i,axb(2)
            !        write(*,*) 'axb3 : ',i,axb(3)
            !  c is position vector of h(ii): calculate c(j),c(k),c(l)
            !##
            do ii=1, 3
                do ind=1, 3
                    c(in(ii),ind) = -coords_tch(in(ii),ind)/bonds_rch(in(ii))
                enddo
            enddo
            !##
            !  argd is the dot product axb dot c
            do ii=1, 3
                argd(in(ii))  = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                argd(in(ii))  = argd(in(ii))/norma
                delta(in(ii)) = acos(argd(in(ii)))-xangles_theta0(i,in(ii))
                !##        write(*,*) 'theta0 ',i,in(ii),theta0(i,in(ii))
                sum2 = sum2+delta(in(ii))**2
                sum4 = sum4+delta(in(ii))**4
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
                cdot(in(ii), 1, 1) = 1.0e0_wp/bonds_rch(in(ii)) +coords_tch(in(ii),1)*c(in(ii),1)/bonds_rch(in(ii))**2
                cdot(in(ii), 1, 2) = coords_tch(in(ii),2)*c(in(ii),1)/bonds_rch(in(ii))**2
                cdot(in(ii), 1, 3) = coords_tch(in(ii),3)*c(in(ii),1)/bonds_rch(in(ii))**2
                !  deriv of cdot(in(ii),y) wrt x, y, z
                cdot(in(ii), 2, 1) = coords_tch(in(ii),1)*c(in(ii),2)/bonds_rch(in(ii))**2
                cdot(in(ii), 2, 2) = 1.0e0_wp/bonds_rch(in(ii)) +coords_tch(in(ii),2)*c(in(ii),2)/bonds_rch(in(ii))**2
                cdot(in(ii), 2, 3) = coords_tch(in(ii),3)*c(in(ii),2)/bonds_rch(in(ii))**2
                !  deriv of cdot(in(ii),z) wrt x, y, z
                cdot(in(ii), 3, 1) = coords_tch(in(ii),1)*c(in(ii),3)/bonds_rch(in(ii))**2
                cdot(in(ii), 3, 2) = coords_tch(in(ii),2)*c(in(ii),3)/bonds_rch(in(ii))**2
                cdot(in(ii), 3, 3) = 1.0e0_wp/bonds_rch(in(ii)) +coords_tch(in(ii),3)*c(in(ii),3)/bonds_rch(in(ii))**2
            enddo
            do ii=1, 3
                do ind=1, 3
                    deldot = -xangles_dtheta0(i,in(ii),i)
                    !  derivative wrt h(i,ind)
                    !  for  rch(i) only terms are from the derivatives of theta0
                    deldot = -deldot*coords_tch(i,ind)/bonds_rch(i)
                    qpdot_pl_pdot(ndx_nh(i,ind))=qpdot_pl_pdot(ndx_nh(i,ind)) +2.0e0_wp*delta1_fdelta(i)* &
                        & delta(in(ii))* deldot +4.0e0_wp*delta1_hdelta(i)*delta(in(ii))**3* deldot
                    !  derivative wrt c(ind)
                    deldot = -deldot
                    qpdot_pl_pdot(ndx_nc(ind)) = qpdot_pl_pdot(ndx_nc(ind)) +2.0e0_wp*delta1_fdelta(i)*delta(in(ii))* &
                        & deldot +4.0e0_wp*delta1_hdelta(i)*delta(in(ii))**3*deldot
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
                            atemp4 = axb(1)*cdot(in(ii), 1, ind) +axb(2)*cdot(in(ii), 2, ind) +axb(3)*cdot(in(ii), 3, &
                                & ind)
                            atemp4 = atemp4/norma
                        endif
                        !  atemp5 is deriv of theta0(i,in(ii)) wrt to nh(in(jj),ind)
                        atemp5 = -xangles_dtheta0(i,in(ii),in(jj))
                        !  deriv wrt h(in(jj)),ind):
                        atemp5 = -atemp5*coords_tch(in(jj),ind)/bonds_rch(in(jj))
                        deldot = atemp3+atemp4
                        deldot = -1.0e0_wp/sqrt(1.0e0_wp-argd(in(ii))**2)*deldot
                        deldot = deldot+atemp5
                        qpdot_pl_pdot(ndx_nh(in(jj),ind))=qpdot_pl_pdot(ndx_nh(in(jj),ind)) +2.0e0_wp* &
                            & delta1_fdelta(i)* delta(in(ii))*deldot +4.0e0_wp*delta1_hdelta(i)* delta(in(ii))**3*deldot
                        !  for carbon the only contributions are from axb dot cdot term and
                        !  from theta0 and derivative cdot wrt carbon=-cdot wrt hydrogen
                        deldot = 1.0e0_wp/sqrt(1.0e0_wp-argd(in(ii))**2)*atemp4
                        deldot = deldot-atemp5
                        qpdot_pl_pdot(ndx_nc(ind))=qpdot_pl_pdot(ndx_nc(ind)) +2.0e0_wp*delta1_fdelta(i)* &
                            & delta(in(ii))* deldot +4.0e0_wp*delta1_hdelta(i)*delta(in(ii))**3* deldot
                    enddo
                enddo
            enddo
        end subroutine calcdelta

        !> @brief Evaluate opforce for the C2H7_VBMM potential.
        subroutine opforce
            !  calculates the out-of-plane bending force constants
            !  and their derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: dswitch(4,4), switch(4)
            integer             :: i, j
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
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
                delta1_fdelta(i) = switch(i)*potcm_fch3
                delta1_hdelta(i) = switch(i)*potcm_hch3
                do j=1, 4
                    delta2_dfdelta(i,j) = dswitch(i,j)*potcm_fch3
                    delta2_dhdelta(i,j) = dswitch(i,j)*potcm_hch3
                enddo
            enddo
        end subroutine opforce

        !> @brief Evaluate ipforce for the C2H7_VBMM potential.
        subroutine ipforce
            !  calculates the symmetrised in plane bend force constants and
            !  all partial derivatives involving them
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: a1, a2, arga1, arga2, df1dc(4), df1dh(4), dfk0(4,4,4), duma1, duma2, f0, f2
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !  set force constant at asymptotes
            f0 = potcm_fkinf+potcm_ak
            f2 = potcm_fkinf
            !##
            force1_fk0(1,2) = f0+f0*(ip1_s1(1)*ip1_s1(2)-1.0e0_wp)+(f0-f2)*(ip1_s2(3)*ip1_s2(4)-1.0e0_wp)
            force1_fk0(1,3) = f0+f0*(ip1_s1(1)*ip1_s1(3)-1.0e0_wp)+(f0-f2)*(ip1_s2(2)*ip1_s2(4)-1.0e0_wp)
            force1_fk0(1,4) = f0+f0*(ip1_s1(1)*ip1_s1(4)-1.0e0_wp)+(f0-f2)*(ip1_s2(2)*ip1_s2(3)-1.0e0_wp)
            force1_fk0(2,3) = f0+f0*(ip1_s1(2)*ip1_s1(3)-1.0e0_wp)+(f0-f2)*(ip1_s2(1)*ip1_s2(4)-1.0e0_wp)
            force1_fk0(2,4) = f0+f0*(ip1_s1(2)*ip1_s1(4)-1.0e0_wp)+(f0-f2)*(ip1_s2(1)*ip1_s2(3)-1.0e0_wp)
            force1_fk0(3,4) = f0+f0*(ip1_s1(3)*ip1_s1(4)-1.0e0_wp)+(f0-f2)*(ip1_s2(1)*ip1_s2(2)-1.0e0_wp)
            !##
            !           write(*,*)
            !           write(*,*) 'f0    : ',f0
            !           write(*,*) 'f2    : ',f2
            !           write(*,*) 's1(1) : ',s1(1)
            !           write(*,*) 's1(2) : ',s1(2)
            !           write(*,*) 's2(3) : ',s2(3)
            !           write(*,*) 's2(4) : ',s2(4)
            !           write(*,*)
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
            do i=1, 4
                !  calculate the terms f1(i)
                !##
                arga1        = potcm_aa1*bonds_rbh(i)*bonds_rbh(i)
                arga2        = potcm_aa4*(bonds_rbh(i)-potcm_r0hh)*(bonds_rbh(i)-potcm_r0hh)
                a1           = 1.0e0_wp-exp(-arga1)
                a2           = potcm_aa2+potcm_aa3*exp(-arga2)
                force1_f1(i) = a1*exp(-a2*(bonds_rch(i)-potcm_r0ch)**2)
                !           write(*,*) 'A1 :',a1
                !           write(*,*) 'A2 :',a2
                !  and calculate the derivatives wrt rch(i) and rbh(i)
                duma1    = 2.0e0_wp*potcm_aa1*bonds_rbh(i)*exp(-arga1)
                duma2    = -2.0e0_wp*potcm_aa3*potcm_aa4*(bonds_rbh(i)-potcm_r0hh)*exp(-arga2)
                df1dc(i) = -2.0e0_wp*(bonds_rch(i)-potcm_r0ch)*a1*a2*exp(-a2*(bonds_rch(i)-potcm_r0ch)**2)
                df1dh(i) = duma1*exp(-a2*(bonds_rch(i)-potcm_r0ch)**2) -duma2*(bonds_rch(i)-potcm_r0ch)**2*a1*exp(-a2* &
                    & (bonds_rch(i)-potcm_r0ch)**2)
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
            !  only non-zero derivatives are those from rbh(i) and rbh(j)
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

        !> @brief Evaluate switchf for the C2H7_VBMM potential.
        subroutine switchf
            !  calculates switching functions: s3,sphi,stheta
            !  and their derivatives ds3,dsphi,dstheta
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: argmax, args1, args2, args3, argsphi, argstheta
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !##
            fsw1_a1s = 1.5313681e-7_wp
            fsw1_b1s = -4.6696246e0_wp
            fsw1_a2s = 1.0147402e-7_wp
            fsw1_b2s = -12.363798e0_wp
            !  use double precision criterion:
            !  tanh(19.0d0)=1.0d0
            argmax = 19.0e0_wp
            do i=1, 4
                args1 = fsw1_a1s*(bonds_rch(i)-potcm_r0ch)*(bonds_rch(i)-fsw1_b1s)**8
                if( args1 .lt. argmax )then
                    ip1_s1(i) = 1.0e0_wp-tanh(args1)
                    ip1_ds1(i) = fsw1_a1s*((bonds_rch(i)-fsw1_b1s)**8 +8.0e0_wp*(bonds_rch(i)-potcm_r0ch)* &
                        & (bonds_rch(i)- fsw1_b1s)**7)
                    ip1_ds1(i) = -ip1_ds1(i)/cosh(args1)**2
                else
                    ip1_s1(i)  = 0.0e0_wp
                    ip1_ds1(i) = 0.0e0_wp
                endif
                args2 = fsw1_a2s*(bonds_rch(i)-potcm_r0ch)*(bonds_rch(i)-fsw1_b2s)**6
                if( args2 .lt. argmax )then
                    ip1_s2(i) = 1.0e0_wp-tanh(args2)
                    ip1_ds2(i) = fsw1_a2s*((bonds_rch(i)-fsw1_b2s)**6 +6.0e0_wp*(bonds_rch(i)-potcm_r0ch)* &
                        & (bonds_rch(i)- fsw1_b2s)**5)
                    ip1_ds2(i) = -ip1_ds2(i)/cosh(args2)**2
                else
                    ip1_s2(i)  = 0.0e0_wp
                    ip1_ds2(i) = 0.0e0_wp
                endif
                !##          write(*,*) i,s1(i),rch(i)
                !  calculate s3 and ds3
                args3 = potcm_a3s*(bonds_rch(i)-potcm_r0ch)*(bonds_rch(i)-potcm_b3s)**2
                if( args3 .lt. argmax )then
                    op1_s3(i) = 1.0e0_wp-tanh(args3)
                    op1_ds3(i) = potcm_a3s*(3.0e0_wp*bonds_rch(i)**2-2.0e0_wp*bonds_rch(i)*(potcm_r0ch+2.0e0_wp* &
                        & potcm_b3s) + potcm_b3s*(potcm_b3s+2.0e0_wp*potcm_r0ch))
                    op1_ds3(i) = -op1_ds3(i)/cosh(args3)**2
                else
                    op1_s3(i)  = 0.0e0_wp
                    op1_ds3(i) = 0.0e0_wp
                endif
                !  calculate sphi and dsphi
                !  condition here is on the bondlength rch(i)
                !  st argsphi is lt approx 19.0d0
                if( bonds_rch(i) .lt. 3.8e0_wp )then
                    argsphi         = potcm_aphi*(bonds_rch(i)-potcm_r0ch)*exp(potcm_bphi*(bonds_rch(i)-potcm_cphi)**3)
                    switch1_sphi(i) = 1.0e0_wp-tanh(argsphi)
                    switch1_dsphi(i) = potcm_aphi*(1.0e0_wp+3.0e0_wp*potcm_bphi*(bonds_rch(i)-potcm_r0ch) * &
                        & (bonds_rch(i)- potcm_cphi)**2)
                    switch1_dsphi(i) = switch1_dsphi(i)*exp(potcm_bphi*(bonds_rch(i)-potcm_cphi)**3)
                    switch1_dsphi(i) = -switch1_dsphi(i)/cosh(argsphi)**2
                else
                    switch1_sphi(i)  = 0.0e0_wp
                    switch1_dsphi(i) = 0.0e0_wp
                endif
                !  calculate stheta and dstheta
                if( bonds_rch(i) .lt. 3.8e0_wp )then
                    argstheta         = potcm_atheta*(bonds_rch(i)-potcm_r0ch)*exp(potcm_btheta*(bonds_rch(i)- &
                        & potcm_ctheta)**3)
                    switch1_stheta(i) = 1.0e0_wp-tanh(argstheta)
                    switch1_dstheta(i) = potcm_atheta*(1.0e0_wp+3.0e0_wp*potcm_btheta*(bonds_rch(i)-potcm_r0ch) * &
                        & (bonds_rch(i)-potcm_ctheta)**2)
                    switch1_dstheta(i) = switch1_dstheta(i)*exp(potcm_btheta*(bonds_rch(i)-potcm_ctheta)**3)
                    switch1_dstheta(i) = -switch1_dstheta(i)/cosh(argstheta)**2
                else
                    switch1_stheta(i)  = 0.0e0_wp
                    switch1_dstheta(i) = 0.0e0_wp
                endif
            enddo
        end subroutine switchf

        !> @brief Evaluate prepare_local_context for the C2H7_VBMM potential.
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
            real(wp)            :: fact1, fact2
            integer             :: i, icount, ind
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            integer, parameter  :: n3tmmn = 18
            !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            infocm_ref(1)     = 'M. J. T. Jordan and R. G. Gilbert'
            infocm_ref(2)     = 'J. Chem. Phys. Vol. 102, p. 5669, 1995'
            infocm_indexes(1) = 1
            infocm_indexes(2) = 6
            infocm_indexes(3) = 1
            infocm_indexes(4) = 1
            infocm_indexes(5) = 1
            infocm_indexes(6) = 1
            infocm_irctnt     = 6
            !  calculate indexes for coordinates
            do ind=1, 3
                icount       = ind-3
                ndx_nc(ind)  = 3*potcm_nnc+icount
                ndx_nhb(ind) = 3*potcm_nnb+icount
                do i=1, 4
                    ndx_nh(i,ind) = 3*potcm_nnh(i)+icount
                enddo
            enddo
            !  convert to appropriate units:
            !  energy   in 1.0d+05 j/mol
            !  time     in 1.0d-14 s
            !  xdistance in 1.0d-10 m
            !  xangles   in radians
            !  mass     in amu
            fact1 = 0.041840e0_wp
            fact2 = 6.022045e0_wp
            !##----------------------------------
            potcm_d1ch  = potcm_d1ch*fact1
            potcm_d3ch  = potcm_d3ch*fact1
            potcm_d1cb  = potcm_d1cb*fact1
            potcm_d3cb  = potcm_d3cb*fact1
            potcm_d1hh  = potcm_d1hh*fact1
            potcm_d3hh  = potcm_d3hh*fact1
            potcm_fch3  = potcm_fch3*fact2
            potcm_hch3  = potcm_hch3*fact2
            potcm_fkinf = potcm_fkinf*fact2
            potcm_ak    = potcm_ak*fact2
            !##----------------------------------
        end subroutine prepare_local_context

        !> @brief Evaluate CSETUP for the C2H7_VBMM potential.
        subroutine csetup( n3tm )
            !       dummy subroutine to make the code compatable with chk.f
            use chempotpy_kinds, only: wp
            implicit none
            integer :: n3tm
        end subroutine csetup

        !> @brief Evaluate xgaussfit for the C2H7_VBMM potential.
        subroutine xgaussfit( r, fit )
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)           :: r
            real(wp)           :: fit
            integer, parameter :: ng = 6
            integer            :: i, j, k
            real(wp)           :: d
            real(wp)           :: ans
            real(wp)           :: r0(ng)
            real(wp)           :: alpha(ng)
            real(wp)           :: cf(ng)
            real(wp)           :: fx
            r0(1) = 1.949230e0_wp
            r0(2) = 1.013860e0_wp
            r0(3) = 2.790270e0_wp
            r0(4) = (r0(1)+r0(2))*0.50e0_wp
            r0(5) = (r0(1)+r0(3))*0.50e0_wp
            r0(6) = (r0(2)+r0(3))*0.50e0_wp
            do i=1, ng
                alpha(i) = 0.250e0_wp
            enddo
            cf(1) = -1178.890067483491520e0_wp
            cf(2) = 32.80653355230102620e0_wp
            cf(3) = 10.31188564665787060e0_wp
            cf(4) = -159.8465369524321600e0_wp
            cf(5) = 17.33838608018073660e0_wp
            cf(6) = 1280.367802966489080e0_wp
            ans   = 0.0e0_wp
            do i=1, ng
                d   = r-r0(i)
                fx  = exp(-alpha(i)*d*d)
                ans = ans+(cf(i)*fx)
            enddo
            fit = ans
        end subroutine xgaussfit
    end subroutine evaluate
end module chempotpy_c2h7_vbmm_pes

!> @brief Expose C2H7_VBMM through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_c2h7_vbmm_pes, only: evaluate
    use chempotpy_kinds,         only: wp
    implicit none
    integer, parameter    :: natoms  = 9
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
