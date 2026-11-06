!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module chempotpy_ch5_gen_j1_1987_pes
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
        integer, parameter            :: natoms  = 6
        integer, parameter            :: natom   = 25
        integer, parameter            :: isurf   = 5
        integer, parameter            :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter            :: n3atom  = 75
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: alpdcm_dalp
        double precision              :: dervcp_deccr1
        double precision              :: dervcp_deccr5
        double precision              :: dervcp_deccr6
        double precision              :: icorcm_ar(15)
        double precision              :: icorcm_der(15)
        double precision              :: infocm_cartnu(natom,3)
        integer                       :: infocm_icartr
        integer                       :: infocm_indexes(natom)
        integer                       :: infocm_irctnt
        integer                       :: infocm_mder
        integer                       :: infocm_msurf
        integer                       :: infocm_natoms
        character                     :: infocm_ref(5)
        double precision              :: pdatcm_aa(3)
        double precision              :: pdatcm_alph(3)
        double precision              :: pdatcm_aparm(5)
        double precision              :: pdatcm_b1
        double precision              :: pdatcm_beta(3)
        double precision              :: pdatcm_bp1
        double precision              :: pdatcm_c1
        double precision              :: pdatcm_cc(3)
        double precision              :: pdatcm_cp
        double precision              :: pdatcm_cp1
        double precision              :: pdatcm_d1(3)
        double precision              :: pdatcm_d3(3)
        double precision              :: pdatcm_re(3)
        double precision              :: pdatcm_refv
        double precision              :: pdatcm_tau
        double precision              :: pdt2cm_a
        double precision              :: pdt2cm_a0
        double precision              :: pdt2cm_a3
        double precision              :: pdt2cm_aphi
        double precision              :: pdt2cm_b
        double precision              :: pdt2cm_b0
        double precision              :: pdt2cm_b3
        double precision              :: pdt2cm_bphi
        double precision              :: pdt2cm_c0
        double precision              :: pdt2cm_cphi
        double precision              :: pdt2cm_d
        double precision              :: pdt2cm_d0
        double precision              :: pdt2cm_delt
        double precision              :: pdt2cm_delta
        double precision              :: pdt2cm_dij
        double precision              :: pdt2cm_fch3
        double precision              :: pdt2cm_hch3
        double precision              :: pdt2cm_pch4
        double precision              :: pdt2cm_rex
        double precision              :: pdt2cm_rnot
        double precision              :: pdt3cm_ca
        double precision              :: pdt3cm_cb
        double precision              :: pdt3cm_fa
        double precision              :: pdt3cm_fk0
        double precision              :: potrcm_acs(6)
        double precision              :: potrcm_ag(6)
        double precision              :: potrcm_ass(6)
        double precision              :: potrcm_r2(15)
        double precision              :: potxcm_dedx(18)
        double precision              :: potxcm_drdx(15,18)
        double precision              :: potxcm_energy
        double precision              :: potxcm_x(6)
        double precision              :: potxcm_xx(6)
        double precision              :: potxcm_y(6)
        double precision              :: potxcm_z(6)
        double precision              :: pt1cm_degsdr(n3atom)
        double precision              :: pt1cm_engygs
        double precision              :: pt1cm_r(n3atom)
        double precision              :: pt3cm_ezero(isurf+1)
        double precision              :: pt4cm_deesdr(n3atom,isurf)
        double precision              :: pt4cm_engyes(isurf)
        double precision              :: pt5cm_deijdr(n3atom,jsurf)
        double precision              :: pt5cm_engyij(jsurf)
        double precision              :: temcom_dtemp(9)
        double precision              :: temcom_etem
        double precision              :: temcom_rtem(9)
        double precision              :: testcm_dgamma(3,18)
        double precision              :: testcm_gamma(3)
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
        alpdcm_dalp      = 0.d0; dervcp_deccr1 = 0.d0; dervcp_deccr5 = 0.d0; dervcp_deccr6 = 0.d0; icorcm_ar = 0.d0
        icorcm_der       = 0.d0; infocm_cartnu = 0.d0; infocm_icartr = 0; infocm_indexes = 0; infocm_irctnt = 0
        infocm_mder      = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref       = ' '
        pdatcm_aa        = 0.d0
        pdatcm_alph      = 0.d0
        pdatcm_aparm     = 0.d0
        pdatcm_b1        = 0.d0
        pdatcm_beta      = 0.d0
        pdatcm_bp1       = 0.d0
        pdatcm_c1        = 0.d0
        pdatcm_cc        = 0.d0
        pdatcm_cp        = 0.d0
        pdatcm_cp1       = 0.d0
        pdatcm_d1        = 0.d0
        pdatcm_d3        = 0.d0
        pdatcm_re        = 0.d0
        pdatcm_refv      = 0.d0
        pdatcm_tau       = 0.d0
        pdt2cm_a         = 0.d0
        pdt2cm_a0        = 0.d0
        pdt2cm_a3        = 0.d0
        pdt2cm_aphi      = 0.d0
        pdt2cm_b         = 0.d0
        pdt2cm_b0        = 0.d0
        pdt2cm_b3        = 0.d0
        pdt2cm_bphi      = 0.d0
        pdt2cm_c0        = 0.d0
        pdt2cm_cphi      = 0.d0
        pdt2cm_d         = 0.d0
        pdt2cm_d0        = 0.d0
        pdt2cm_delt      = 0.d0
        pdt2cm_delta     = 0.d0
        pdt2cm_dij       = 0.d0
        pdt2cm_fch3      = 0.d0
        pdt2cm_hch3      = 0.d0
        pdt2cm_pch4      = 0.d0
        pdt2cm_rex       = 0.d0
        pdt2cm_rnot      = 0.d0
        pdt3cm_ca        = 0.d0
        pdt3cm_cb        = 0.d0
        pdt3cm_fa        = 0.d0
        pdt3cm_fk0       = 0.d0
        potrcm_acs       = 0.d0
        potrcm_ag        = 0.d0
        potrcm_ass       = 0.d0
        potrcm_r2        = 0.d0
        potxcm_dedx      = 0.d0
        potxcm_drdx      = 0.d0
        potxcm_energy    = 0.d0
        potxcm_x         = 0.d0
        potxcm_xx        = 0.d0
        potxcm_y         = 0.d0
        potxcm_z         = 0.d0
        pt1cm_degsdr     = 0.d0
        pt1cm_engygs     = 0.d0
        pt1cm_r          = 0.d0
        pt3cm_ezero      = 0.d0
        pt4cm_deesdr     = 0.d0
        pt4cm_engyes     = 0.d0
        pt5cm_deijdr     = 0.d0
        pt5cm_engyij     = 0.d0
        temcom_dtemp     = 0.d0
        temcom_etem      = 0.d0
        temcom_rtem      = 0.d0
        testcm_dgamma    = 0.d0
        testcm_gamma     = 0.d0
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
        usricm_nder    = 1
        usricm_nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero = 0.0d0
        infocm_icartr  = 1
        infocm_msurf   = 0
        infocm_mder    = 1
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 6
        pdatcm_d1      = [4.8668d0, 4.74660d0, 1.1452d0]
        pdatcm_d3      = [1.6890d0, 1.75200d0, 0.8870d0]
        pdatcm_b1      = 3.00d0
        pdatcm_c1      = 11.00d0
        pdatcm_alph    = [0.939d0, 1.029955d0, 0.939d0]
        pdatcm_re      = [2.0673d0, 1.402d0, 2.0673d0]
        pdatcm_beta    = [0.939d0, 1.029955d0, 0.939d0]
        pdatcm_aa      = [0.0d0, 0.0d0, -2.3603d0]
        pdatcm_refv    = 0.71098329d0
        pdatcm_tau     = 1.9105588d0
        pdatcm_cp      = 0.080768d0
        pdatcm_aparm   = [4.0d0, 0.8999975d0, 0.448035d0, 0.606528d0, 3.239921d0]
        pdt2cm_aphi    = 2.7982366d-1
        pdt2cm_bphi    = 5.937217d-2
        pdt2cm_cphi    = 3.6301534d0
        pdt2cm_pch4    = 1.9106d0
        pdt2cm_fch3    = 2.19d-2
        pdt2cm_hch3    = 4.38d-02
        pdt2cm_rnot    = 2.052d0
        pdt2cm_a3      = 2.102956d-2
        pdt2cm_b3      = 5.798533d-1
        pdt2cm_delt    = 0.1d0
        pdt2cm_dij     = 3.1190d0
        pdt2cm_a       = -1.56944d+3
        pdt2cm_b       = -1.12008d+3
        pdt2cm_c0      = 1.5786151d+3
        pdt2cm_d       = 2.67146d0
        pdt3cm_ca      = 0.903d0
        pdt3cm_cb      = 0.072d0
        pdt2cm_rex     = 1.18584d0
        pdt3cm_fk0     = 2.55d0
        pdt3cm_fa      = 0.8d0
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
        !> @brief Evaluate POT for the CH5_GEN_J1_1987 potential.
        subroutine pot
            !
            !     SUBROUTINE TO CALCULATE THE CH4-H SURFACE
            !     ALSO COMPUTES DE/DR FOR ALL R AND STORES RESULTS IN ARRAY DER
            !     DEFINITION OF CARTESIAN COORDINATES X(I), Y(I), AND Z(I)
            !      I = 1, 2, 3, 4   REFERS TO H1, H2, H3, H4
            !      I = 5            REFERS TO C
            !      I = 6            REFERS TO H'
            !
            !     DEFINITION OF THE R-ARRAY  AR(I)
            !     H1 IS THE ATOM THAT WILL FORM H-H WITH H5
            !      I = 1, 2, 3, 4   REFERS TO C - H(I)
            !      I = 5            REFERS TO C - H'
            !      I = 6, 7, 8, 9   REFERS TO H' - H(I)
            !      I = 10, 11, 12   REFERS TO H1 - H(2,3,4)
            !      I = 13, 14       REFERS TO H2 - H(3,4)
            !      I = 15           REFERS TO H3 - H4
            !
            implicit none
            double precision :: a11, a12, a13, a14, a21, a22, a23, a24, asb, asf, coord, daa, dadr1, dasb, dasf, der2, der3, der4, df, dfr, dk, dr2, dtot, dvoop, dxx, e, ee, evcon, exa11, exa12, exa13, exa14, exa21, exa22, exa23, exa24, exf1, exf2, exf3, exf4, f, o, rhmax, stp, sum, v, voop, w, xa, xa2, xb, xb2, xc, xc2, xd, xd2, xxa, xxb, xxc, xxd, y1, y2
            integer          :: i, ib1, ib2, ib3, if1, if2, if3, ipt, klm, kmax
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
            dimension dfr(4), dvoop(18), y1(6), y2(6), f(4), df(4), dk(6)
            !
            dimension coord(18), dxx(18)
            !
            !
            !
            data evcon / 2.721161d+1 /
            !
            call cartou
            call carttor
            !
            !     PUT COORDINATES IN PROPER ARRAYS
            !
            do 1 i=1, 18
                coord(i) = pt1cm_r(i)
            1 continue
            !
            ipt = 1
            do 5 i = 1, 6
                potxcm_x(i) = coord(ipt)
                potxcm_y(i) = coord(ipt+1)
                potxcm_z(i) = coord(ipt+2)
                ipt         = ipt+3
            5 continue
            !
            ! COMPUTE THE INTERPARTICLE DISTANCES R
            !
            call dista
            !
            ! CHECK FOR THE FIFTH ATOM--I.E.--CH4---(H)
            ! BY DEFINITION THE FIFTH H ATOM IS THE ONE WHOSE C-H DISTANCE IS MAXIM
            !
            rhmax = icorcm_ar(1)
            kmax  = 1
            if( icorcm_ar(2) .lt. rhmax ) go to 10
            rhmax = icorcm_ar(2)
            kmax  = 2
            10 if (icorcm_ar(3).lt.rhmax) go to 20
            rhmax = icorcm_ar(3)
            kmax  = 3
            20 if (icorcm_ar(4).lt.rhmax) go to 30
            rhmax = icorcm_ar(4)
            kmax  = 4
            30 if (icorcm_ar(5).lt.rhmax) go to 40
            rhmax = icorcm_ar(5)
            kmax  = 5
            40 call switch (kmax, 1, v, coord, dxx)
            !
            ! COMPUTE THE HCH ANGLES AND STORE RESULTS IN ARRAY AG
            !
            call angle
            !
            do 50 i = 1, 15
                icorcm_der(i) = 0.0d0
            50 continue
            !
            !  COMPUTE THE MORSE RANGE PARAMETERS ALPH:
            !
            call calpha
            !
            ! COMPUTE THE FOUR TRIATOMIC TERMS
            !
            call tri (icorcm_ar(1), icorcm_ar(5), icorcm_ar(6), icorcm_der(1), icorcm_der(5), icorcm_der(6), 1, 3, 2, ee, dtot)
            icorcm_der(1) = icorcm_der(1)+dtot
            der2          = dtot
            der3          = dtot
            der4          = dtot
            call tri (icorcm_ar(2), icorcm_ar(5), icorcm_ar(7), icorcm_der(2), dr2, icorcm_der(7), 1, 3, 2, e, dtot)
            ee            = ee+e
            icorcm_der(2) = icorcm_der(2)+der2+dtot
            icorcm_der(1) = icorcm_der(1)+dtot
            der3          = der3+dtot
            der4          = der4+dtot
            icorcm_der(5) = icorcm_der(5)+dr2
            call tri (icorcm_ar(3), icorcm_ar(5), icorcm_ar(8), icorcm_der(3), dr2, icorcm_der(8), 1, 3, 2, e, dtot)
            ee            = ee+e
            icorcm_der(3) = icorcm_der(3)+der3+dtot
            icorcm_der(1) = icorcm_der(1)+dtot
            icorcm_der(2) = icorcm_der(2)+dtot
            der4          = der4+dtot
            icorcm_der(5) = icorcm_der(5)+dr2
            call tri (icorcm_ar(4), icorcm_ar(5), icorcm_ar(9), icorcm_der(4), dr2, icorcm_der(9), 1, 3, 2, e, dtot)
            ee            = ee+e
            icorcm_der(4) = icorcm_der(4)+der4+dtot
            icorcm_der(1) = icorcm_der(1)+dtot
            icorcm_der(2) = icorcm_der(2)+dtot
            icorcm_der(3) = icorcm_der(3)+dtot
            icorcm_der(5) = icorcm_der(5)+dr2
            !
            ! COMPUTE FORCE CONSTANTS
            ! COMPUTE ATTENUATION TERMS
            !
            xa    = icorcm_ar(1)-pdatcm_re(1)
            xb    = icorcm_ar(2)-pdatcm_re(1)
            xc    = icorcm_ar(3)-pdatcm_re(1)
            xd    = icorcm_ar(4)-pdatcm_re(1)
            xa2   = xa*xa
            xb2   = xb*xb
            xc2   = xc*xc
            xd2   = xd*xd
            xxa   = icorcm_ar(6)-pdatcm_re(2)
            xxb   = icorcm_ar(7)-pdatcm_re(2)
            xxc   = icorcm_ar(8)-pdatcm_re(2)
            xxd   = icorcm_ar(9)-pdatcm_re(2)
            exa11 = exp(-pdatcm_aparm(2)*potrcm_r2(6))
            exa12 = exp(-pdatcm_aparm(2)*potrcm_r2(7))
            exa13 = exp(-pdatcm_aparm(2)*potrcm_r2(8))
            exa14 = exp(-pdatcm_aparm(2)*potrcm_r2(9))
            exa21 = exp(-pdatcm_aparm(5)*xxa*xxa)
            exa22 = exp(-pdatcm_aparm(5)*xxb*xxb)
            exa23 = exp(-pdatcm_aparm(5)*xxc*xxc)
            exa24 = exp(-pdatcm_aparm(5)*xxd*xxd)
            a11   = 1.0d0-exa11
            a12   = 1.0d0-exa12
            a13   = 1.0d0-exa13
            a14   = 1.0d0-exa14
            a21   = pdatcm_aparm(3)+pdatcm_aparm(4)*exa21
            a22   = pdatcm_aparm(3)+pdatcm_aparm(4)*exa22
            a23   = pdatcm_aparm(3)+pdatcm_aparm(4)*exa23
            a24   = pdatcm_aparm(3)+pdatcm_aparm(4)*exa24
            exf1  = exp(-a21*xa2)
            exf2  = exp(-a22*xb2)
            exf3  = exp(-a23*xc2)
            exf4  = exp(-a24*xd2)
            f(1)  = a11*exf1
            f(2)  = a12*exf2
            f(3)  = a13*exf3
            f(4)  = a14*exf4
            !
            ! DK'S ARE THE FORCE CONSTANTS
            !     FIRST COMPUTE THE K0 SINCE IT IS NOW A FUNCTION OF AR(1)
            !
            call force (dadr1)
            dk(1) = pdatcm_aparm(1)*f(1)*f(2)
            dk(2) = pdatcm_aparm(1)*f(1)*f(3)
            dk(3) = pdatcm_aparm(1)*f(1)*f(4)
            dk(4) = pdatcm_aparm(1)*f(2)*f(3)
            dk(5) = pdatcm_aparm(1)*f(2)*f(4)
            dk(6) = pdatcm_aparm(1)*f(3)*f(4)
            !
            ! COMPUTE THE EQUILIBRIUM ANGLES--THETA AND D(THETA)/DR
            !
            call theta (asf, dasf, asb, dasb)
            daa = 0.0d0
            !
            ! COMPUTE CONTRIBUTION TO THE TOTAL ENERGY FROM THE ANGLE TERMS
            !
            sum = 0.0d0
            if1 = 1
            if2 = 2
            if3 = 3
            ib1 = 4
            ib2 = 5
            ib3 = 6
            do 60 i = 1, 3
                o     = potrcm_ag(i)-asf
                y1(i) = o
                y2(i) = o*o
                sum   = sum+dk(i)*y2(i)
            60 continue
            do 70 i = 4, 6
                y1(i) = potrcm_ag(i)-asb
                y2(i) = y1(i)*y1(i)
                sum   = sum+dk(i)*y2(i)
            70 continue
            ee = ee+0.5d0*sum
            !
            ! COMPUTE DERIVATIVES OF THE ANGLE ATTENUATION TERMS WITH RESPECT TO R6,
            !                       AND R9
            !
            dfr(1) = 2.0d0*pdatcm_aparm(2)*icorcm_ar(6)*exa11*exf1+a11*exf1*xa2*2.0d0* pdatcm_aparm(4)* &
                & pdatcm_aparm(5)*xxa*exa21
            dfr(2) = 2.0d0*pdatcm_aparm(2)*icorcm_ar(7)*exa12*exf2+a12*exf2*xb2*2.0d0* pdatcm_aparm(4)* &
                & pdatcm_aparm(5)*xxb*exa22
            dfr(3) = 2.0d0*pdatcm_aparm(2)*icorcm_ar(8)*exa13*exf3+a13*exf3*xc2*2.0d0* pdatcm_aparm(4)* &
                & pdatcm_aparm(5)*xxc*exa23
            dfr(4) = 2.0d0*pdatcm_aparm(2)*icorcm_ar(9)*exa14*exf4+a14*exf4*xd2*2.0d0* pdatcm_aparm(4)* &
                & pdatcm_aparm(5)*xxd*exa24
            !
            ! COMPUTE DERIVATIVES OF ATTENUATION TERMS WITH RESPECT TO R1,R2,R3, AND
            !
            df(1) = -2.0d0*a21*xa*f(1)
            df(2) = -2.0d0*a22*xb*f(2)
            df(3) = -2.0d0*a23*xc*f(3)
            df(4) = -2.0d0*a24*xd*f(4)
            !
            ! ADD IN CONTRIBUTION TO DER FROM DK/DR ANGLE TERMS AND D(THETA)/DR TERM
            ! FIRST THREE ANGLE TERMS
            !
            do 80 i = 1, 3
                stp = 0.5d0*dk(i)*y2(i)
                w   = dk(i)*y1(i)
                if( abs(f(1)) .gt. 1.0d-10 ) icorcm_der(1)=icorcm_der(1)+stp*(df(1)/f(1)+dadr1/ pdatcm_aparm(1))+w*(- &
                    & 1.0d0/(icorcm_ar(i+1)*potrcm_ass(i))+potrcm_acs(i)/ (icorcm_ar(1)*potrcm_ass(i)))
                if( abs(f(i+1)) .gt. 1.0d-10 ) icorcm_der(i+1)=icorcm_der(i+1)+stp*df(i+1)/f(i+ 1)+w*(-1.0d0/ &
                    & (icorcm_ar(1)*potrcm_ass(i))+potrcm_acs(i)/(potrcm_ass(i)* icorcm_ar(i+1)))
                if( abs(f(1)) .gt. 1.0d-10 ) icorcm_der(6)     = icorcm_der(6)+stp*dfr(1)/f(1)
                if( abs(f(i+1)) .gt. 1.0d-10 ) icorcm_der(i+6) = icorcm_der(i+6)+stp*dfr(i+1)/f(i +1)
                !
                icorcm_der(i+9) = icorcm_der(i+9)+w*icorcm_ar(i+9)/(icorcm_ar(1)*icorcm_ar(i+1)*potrcm_ass(i))
            80 continue
            !
            ! FOURTH AND FIFTH ANGLE TERMS
            !
            do 90 i = 4, 5
                stp           = 0.5d0*dk(i)*y2(i)
                w             = dk(i)*y1(i)
                icorcm_der(1) = icorcm_der(1)+stp*dadr1/pdatcm_aparm(1)
                icorcm_der(2) = icorcm_der(2)+stp*df(2)/f(2)+w*(-1.0d0/(potrcm_ass(i)*icorcm_ar(i-1))+ potrcm_acs(i)/ &
                    & (potrcm_ass(i)*icorcm_ar(2)))
                icorcm_der(7) = icorcm_der(7)+stp*dfr(2)/f(2)
                if( abs(f(i-1)) .gt. 1.0d-10 ) icorcm_der(i-1)=icorcm_der(i-1)+stp*df(i-1)/f(i- 1)+w*(-1.0d0/ &
                    & (potrcm_ass(i)*icorcm_ar(2))+potrcm_acs(i)/(potrcm_ass(i)* icorcm_ar(i-1)))
                !
                if( abs(f(i-1)) .gt. 1.0d-10 ) icorcm_der(i+4) = icorcm_der(i+4)+stp*dfr(i-1)/f(i -1)
                !
                icorcm_der(i+9) = icorcm_der(i+9)+w*icorcm_ar(i+9)/(icorcm_ar(2)*icorcm_ar(i-1)*potrcm_ass(i))
            90 continue
            !
            ! CONTRIBUTION FROM SIXTH ANGLE TERM
            !
            stp           = 0.5d0*dk(6)*y2(6)
            w             = dk(6)*y1(6)
            icorcm_der(1) = icorcm_der(1)+stp*dadr1/pdatcm_aparm(1)
            icorcm_der(3)  = icorcm_der(3)+stp*df(3)/f(3)+w*(-1.0d0/(potrcm_ass(6)*icorcm_ar(4))+potrcm_acs(6)/ &
                & (potrcm_ass(6)*icorcm_ar(3)))
            icorcm_der(8) = icorcm_der(8)+stp*dfr(3)/f(3)
            icorcm_der(4)  = icorcm_der(4)+stp*df(4)/f(4)+w*(-1.0d0/(potrcm_ass(6)*icorcm_ar(3))+potrcm_acs(6)/ &
                & (potrcm_ass(6)*icorcm_ar(4)))
            icorcm_der(9)  = icorcm_der(9)+stp*dfr(4)/f(4)
            icorcm_der(15) = icorcm_der(15)+w*icorcm_ar(15)/(potrcm_ass(6)*icorcm_ar(3)*icorcm_ar(4))
            !
            ! CONTRIBUTION FROM EQ. ANGLE VARIATION TO D(ANGLE)/DR
            !
            do 100 klm = 1, 4
                icorcm_der(klm) = icorcm_der(klm)-(dk(if1)*y1(if1)+dk(if2)*y1(if2)+dk(if3)*y1 (if3))*dasf-(dk(ib1)* &
                    & y1(ib1)+dk(ib2)*y1(ib2)+ dk(ib3)*y1(ib3) )*dasb
            100 continue
            icorcm_der(5) = icorcm_der(5)-(dk(if1)*y1(if1)+dk(if2)*y1(if2)+dk(if3)*y1(if3))* daa-(dk(ib1)*y1(ib1)+ &
                & dk(ib2)*y1(ib2)+dk(ib3)* y1(ib3))*3.0d0* sin(asf)*cos(asf)*daa/sin(asb)
            !
            call switch (kmax, 2, v, coord, dxx)
            !
            potxcm_energy = ee/evcon+pdatcm_refv
            do 110 i = 1, 15
                icorcm_der(i) = icorcm_der(i)/evcon
            110 continue
            !
            !      CALCULATE DRDX ARRAY
            !
            call drdxc
            !
            !      CALCULATE DEDX = DEDR(DER) * DRDX
            !
            call mult (icorcm_der, 1, 15, potxcm_drdx, 15, 18, potxcm_dedx)
            ipt = 1
            do 120 i = 1, 6
                dxx(ipt)   = potxcm_dedx(i)
                dxx(ipt+1) = potxcm_dedx(i+6)
                dxx(ipt+2) = potxcm_dedx(i+12)
                ipt        = ipt+3
            120 continue
            !
            !     NOW ADD IN THE OUT OF PLANE BEND TERM
            !
            call oopb (voop, dvoop, v, coord, dxx)
            potxcm_energy = potxcm_energy+voop
            do 130 i = 1, 18
                dxx(i) = dxx(i)+dvoop(i)
            130 continue
            !
            pt1cm_engygs = potxcm_energy
            call eunitzero
            if( usricm_nder .ne. 0 )then
                do i=1, 18
                    pt1cm_degsdr(i) = dxx(i)
                enddo
                call rtocart
                call dedcou
            endif
            !
        end subroutine pot
        !> @brief Evaluate ANGLE for the CH5_GEN_J1_1987 potential.
        subroutine angle
            !
            ! ROUTINE TO COMPUTE CH4 ANGLES FROM THE INTERPARTICLE DISTANCES
            !
            implicit none
            double precision :: x
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
            do 10 i = 1, 3
                x             = (potrcm_r2(1)+potrcm_r2(i+1)-potrcm_r2(i+9))/(2.0d0*icorcm_ar(1)*icorcm_ar(i+1))
                potrcm_acs(i) = x
                potrcm_ass(i) = sqrt(1.0d0-x*x)
                potrcm_ag(i)  = acos(x)
            10 continue
            do 20 i = 1, 2
                x               = (potrcm_r2(2)+potrcm_r2(i+2)-potrcm_r2(i+12))/(2.0d0*icorcm_ar(2)*icorcm_ar(i+2))
                potrcm_acs(i+3) = x
                potrcm_ass(i+3) = sqrt(1.0d0-x*x)
                potrcm_ag(i+3)  = acos(x)
            20 continue
            x             = (potrcm_r2(3)+potrcm_r2(4)-potrcm_r2(15))/(2.0d0*icorcm_ar(3)*icorcm_ar(4))
            potrcm_acs(6) = x
            potrcm_ass(6) = sqrt(1.0d0-x*x)
            potrcm_ag(6)  = acos(x)
        end subroutine angle
        !> @brief Evaluate CALPHA for the CH5_GEN_J1_1987 potential.
        subroutine calpha
            !
            !     COMPUTE THE C-H AND C-H' MORSE RANGE PARAMETERS FOR THE SINGLET
            !     AND TRIPLET ENERGY CURVES.  THESE WILL BE A FUNCTION OF THE
            !     AVERAGE C-H BOND LENGTH IN METHANE.
            !
            implicit none
            double precision :: htan, ravg
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
            !     FIRST COMPUTE THE AVERAGE METHANE BOND LENGTH
            !
            ravg = (icorcm_ar(1)+icorcm_ar(2)+icorcm_ar(3)+icorcm_ar(4))/4.0d0
            !
            htan           = tanh(pdatcm_c1*(ravg-pdatcm_re(1)))
            pdatcm_alph(1) = pdt3cm_ca+pdt3cm_cb*((htan+1.0d0)/2.0d0)
            pdatcm_alph(3) = pdatcm_alph(1)
            pdatcm_beta(1) = pdatcm_alph(1)
            pdatcm_beta(3) = pdatcm_alph(1)
            !
            !     NOW THE DERIVATIVES---DALPH WILL BE THE DERIVATIVE OF ALPHA WITH
            !     RESPECT TO R1,R2,R3 OR R4
            !
            alpdcm_dalp = (pdatcm_c1*pdt3cm_cb/8.0d0)*(1.0d0-htan**2)
            !
        end subroutine calpha
        !> @brief Evaluate DISTA for the CH5_GEN_J1_1987 potential.
        subroutine dista
            !
            ! ROUTINE TO COMPUTE THE INTERPARTICLE DISTANCES FROM THE CARTESIAN
            !     COORDINATES OF THE ATOMS
            !
            implicit none
            double precision :: df
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
            dimension df(15,3)
            !
            !
            !
            do 10 i = 1, 4
                df(i,1)      = potxcm_x(i)-potxcm_x(5)
                df(i,2)      = potxcm_y(i)-potxcm_y(5)
                df(i,3)      = potxcm_z(i)-potxcm_z(5)
                potrcm_r2(i) = df(i,1)*df(i,1)+df(i,2)*df(i,2)+df(i,3)*df(i,3)
                icorcm_ar(i) = sqrt(potrcm_r2(i))
            10 continue
            df(5,1)      = potxcm_x(5)-potxcm_x(6)
            df(5,2)      = potxcm_y(5)-potxcm_y(6)
            df(5,3)      = potxcm_z(5)-potxcm_z(6)
            potrcm_r2(5) = df(5,1)*df(5,1)+df(5,2)*df(5,2)+df(5,3)*df(5,3)
            icorcm_ar(5) = sqrt(potrcm_r2(5))
            do 20 i = 1, 4
                j            = i+5
                df(j,1)      = potxcm_x(i)-potxcm_x(6)
                df(j,2)      = potxcm_y(i)-potxcm_y(6)
                df(j,3)      = potxcm_z(i)-potxcm_z(6)
                potrcm_r2(j) = df(j,1)*df(j,1)+df(j,2)*df(j,2)+df(j,3)*df(j,3)
                icorcm_ar(j) = sqrt(potrcm_r2(j))
            20 continue
            do 30 i = 1, 3
                j            = i+1
                k            = i+9
                df(k,1)      = potxcm_x(1)-potxcm_x(j)
                df(k,2)      = potxcm_y(1)-potxcm_y(j)
                df(k,3)      = potxcm_z(1)-potxcm_z(j)
                potrcm_r2(k) = df(k,1)*df(k,1)+df(k,2)*df(k,2)+df(k,3)*df(k,3)
                icorcm_ar(k) = sqrt(potrcm_r2(k))
            30 continue
            do 40 i = 1, 2
                j            = i+2
                k            = i+12
                df(k,1)      = potxcm_x(2)-potxcm_x(j)
                df(k,2)      = potxcm_y(2)-potxcm_y(j)
                df(k,3)      = potxcm_z(2)-potxcm_z(j)
                potrcm_r2(k) = df(k,1)*df(k,1)+df(k,2)*df(k,2)+df(k,3)*df(k,3)
                icorcm_ar(k) = sqrt(potrcm_r2(k))
            40 continue
            df(15,1)      = potxcm_x(3)-potxcm_x(4)
            df(15,2)      = potxcm_y(3)-potxcm_y(4)
            df(15,3)      = potxcm_z(3)-potxcm_z(4)
            potrcm_r2(15) = df(15,1)*df(15,1)+df(15,2)*df(15,2)+df(15,3)*df(15,3)
            icorcm_ar(15) = sqrt(potrcm_r2(15))
        end subroutine dista
        !> @brief Evaluate DRDXC for the CH5_GEN_J1_1987 potential.
        subroutine drdxc
            !
            ! THIS SUBROUTINE CALCULATES DRDX
            !
            !     ADDED BY FBB AT THE U. OF MINN.  1/24/84
            !
            implicit none
            integer :: i, j
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
            !
            !
            ! ZERO OUT DRDX
            !
            do 10 i = 1, 18
                do 10 j = 1, 15
                    potxcm_drdx(j,i) = 0.0d0
            10 continue
            !
            ! CALCULATE DRDX FOR C - H(1,2,3,4) BONDS
            !
            do 20 i = 1, 4
                potxcm_drdx(i,i)    = (potxcm_x(i)-potxcm_x(5))/icorcm_ar(i)
                potxcm_drdx(i,5)    = -potxcm_drdx(i,i)
                potxcm_drdx(i,i+6)  = (potxcm_y(i)-potxcm_y(5))/icorcm_ar(i)
                potxcm_drdx(i,11)   = -potxcm_drdx(i,i+6)
                potxcm_drdx(i,i+12) = (potxcm_z(i)-potxcm_z(5))/icorcm_ar(i)
                potxcm_drdx(i,17)   = -potxcm_drdx(i,i+12)
            20 continue
            !
            ! CALCULATE DRDX FOR C - H' BOND
            !
            potxcm_drdx(5,5)  = (potxcm_x(5)-potxcm_x(6))/icorcm_ar(5)
            potxcm_drdx(5,6)  = -potxcm_drdx(5,5)
            potxcm_drdx(5,11) = (potxcm_y(5)-potxcm_y(6))/icorcm_ar(5)
            potxcm_drdx(5,12) = -potxcm_drdx(5,11)
            potxcm_drdx(5,17) = (potxcm_z(5)-potxcm_z(6))/icorcm_ar(5)
            potxcm_drdx(5,18) = -potxcm_drdx(5,17)
            !
            ! CALCULATE DRDX FOR H - H' BONDS
            !
            do 30 i = 1, 4
                potxcm_drdx(5+i,i)    = (potxcm_x(i)-potxcm_x(6))/icorcm_ar(5+i)
                potxcm_drdx(5+i,6)    = -potxcm_drdx(5+i,i)
                potxcm_drdx(5+i,i+6)  = (potxcm_y(i)-potxcm_y(6))/icorcm_ar(5+i)
                potxcm_drdx(5+i,12)   = -potxcm_drdx(5+i,i+6)
                potxcm_drdx(5+i,i+12) = (potxcm_z(i)-potxcm_z(6))/icorcm_ar(5+i)
                potxcm_drdx(5+i,18)   = -potxcm_drdx(5+i,i+12)
            30 continue
            !
            ! CALCULATE DRDX FOR H1 - H(2,3,4) BONDS
            !
            do 40 i = 1, 3
                potxcm_drdx(9+i,i+1)  = -(potxcm_x(1)-potxcm_x(i+1))/icorcm_ar(9+i)
                potxcm_drdx(9+i,1)    = -potxcm_drdx(9+i,i+1)
                potxcm_drdx(9+i,i+7)  = -(potxcm_y(1)-potxcm_y(i+1))/icorcm_ar(9+i)
                potxcm_drdx(9+i,7)    = -potxcm_drdx(9+i,i+7)
                potxcm_drdx(9+i,i+13) = -(potxcm_z(1)-potxcm_z(i+1))/icorcm_ar(9+i)
                potxcm_drdx(9+i,13)   = -potxcm_drdx(9+i,i+13)
            40 continue
            !
            ! CALCULATE DRDX FOR H2 - H(3,4)
            !
            do 50 i = 1, 2
                potxcm_drdx(12+i,i+2)  = -(potxcm_x(2)-potxcm_x(i+2))/icorcm_ar(12+i)
                potxcm_drdx(12+i,2)    = -potxcm_drdx(12+i,i+2)
                potxcm_drdx(12+i,i+8)  = -(potxcm_y(2)-potxcm_y(i+2))/icorcm_ar(12+i)
                potxcm_drdx(12+i,8)    = -potxcm_drdx(12+i,i+8)
                potxcm_drdx(12+i,i+14) = -(potxcm_z(2)-potxcm_z(i+2))/icorcm_ar(12+i)
                potxcm_drdx(12+i,14)   = -potxcm_drdx(12+i,i+14)
            50 continue
            !
            ! CALCULATE DRDX FOR H3 - H4
            !
            potxcm_drdx(15,4)  = -(potxcm_x(3)-potxcm_x(4))/icorcm_ar(15)
            potxcm_drdx(15,3)  = -potxcm_drdx(15,4)
            potxcm_drdx(15,10) = -(potxcm_y(3)-potxcm_y(4))/icorcm_ar(15)
            potxcm_drdx(15,9)  = -potxcm_drdx(15,10)
            potxcm_drdx(15,16) = -(potxcm_z(3)-potxcm_z(4))/icorcm_ar(15)
            potxcm_drdx(15,15) = -potxcm_drdx(15,16)
            !
        end subroutine drdxc
        !> @brief Evaluate FORCE for the CH5_GEN_J1_1987 potential.
        subroutine force ( dadr1 )
            !
            !     THIS SUBROUTINE COMPUTES THE BENDING FORCE CONSTANT AS A FUNCTION
            !     OF AR(1). (APARM(1))   IT ALSO COMPUTES THE DERIVATIVE OF THIS
            !     TERM WITH RESPECT TO AR(1).
            !
            implicit none
            double precision :: dadr1, fex
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
            !
            fex             = exp(-pdatcm_b1*(icorcm_ar(1)-pdatcm_re(1))**2)
            pdatcm_aparm(1) = pdt3cm_fk0+pdt3cm_fa*fex
            dadr1           = -2.0d0*pdt3cm_fa*pdatcm_b1*(icorcm_ar(1)-pdatcm_re(1))*fex
        end subroutine force
        !> @brief Evaluate MULT for the CH5_GEN_J1_1987 potential.
        subroutine mult ( a, na, ma, b, nb, mb, c )
            implicit none
            double precision :: a, b, c, sum
            integer          :: i, j, k, ma, mb, na, nb
            !
            ! THIS SUBROUTINE PERFORMS THE MATRIX MULTIPLICATION :
            !
            !         C = A X B
            !
            !
            dimension a(na,ma), b(nb,mb), c(na,mb)
            !
            do 30 i = 1, na
                !
                do 20 j = 1, mb
                    !
                    sum = 0.0d0
                    do 10 k = 1, ma
                        !
                        sum = sum+a(i,k)*b(k,j)
                    10 continue
                    !
                    c(i,j) = sum
                20 continue
                !
            30 continue
        end subroutine mult
        !> @brief Evaluate OOPB for the CH5_GEN_J1_1987 potential.
        subroutine oopb ( voop, dvoop, v, x, dx )
            !
            !     THIS SUBROUTINE CALCULTES THE OUT OF PLANE BEND POTENTIAL THAT
            !     RESULTS FROM THE METHYL.  THE FORM USED IS THE SAME AS THAT
            !     GIVEN BY DUCHOVIC ET. AL. (JPC 89, 1339, (1984)).
            !     THE DERIVATIVES WITH RESPECT TO THE CARTESIAN COORDINATES (X)
            !     ARE ALSO COMPUTED AND STORED IN DVOOP.
            !
            implicit none
            double precision :: crossx, crossy, crossz, cst, dcross, ddelta, delta, dfds, dfdx, dgcros, dhds, dhdx, dphi, drdx, drv, dsdr, dtop, dvoop, dx, dxmagr, erest, f, gcross, h, phinot, rest, rstar, rv, s3, sphi, term1, term2, term3, term4, texp, top, v, voop, x, xdif, xmagr, ydif, zdif
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
            dimension dvoop(18), xdif(2), ydif(2), zdif(2), rv(3,3), drv(3,3,18)
            dimension xmagr(3), dxmagr(3,18), delta(3), ddelta(3,18)
            dimension top(3), dtop(3,18), dcross(3,18)
            dimension dgcros(18), dphi(18), dfdx(18), dhdx(18), drdx(18)
            !
            !
            dimension x(18), dx(18)
            !
            !
            !     START BY COMPUTING THE NECESSARY DIFFERENCES IN THE R1,R2 AND
            !     R3 VECTORS.
            !
            xdif(1) = x(7)-x(4)
            ydif(1) = x(8)-x(5)
            zdif(1) = x(9)-x(6)
            xdif(2) = x(10)-x(4)
            ydif(2) = x(11)-x(5)
            zdif(2) = x(12)-x(6)
            !
            !     NOW COMPUTE THE INDIVIDUAL COMPONENTS OF THE R1,R2 AND R3 VECTORS
            !
            j = 1
            do 10 i = 1, 3
                rv(i,1) = x(j+3)-x(13)
                rv(i,2) = x(j+4)-x(14)
                rv(i,3) = x(j+5)-x(15)
                j       = j+3
            10 continue
            !
            !     RSTAR IS THE C-H4 BOND LENGTH
            !
            rstar = sqrt((x(1)-x(13))**2+(x(2)-x(14))**2+(x(3)-x(15))**2)
            !
            !     COMPUTE THE CROSS PRODUCT OF THESE DIFFERENCES
            !
            crossx = ydif(1)*zdif(2)-ydif(2)*zdif(1)
            crossy = zdif(1)*xdif(2)-xdif(1)*zdif(2)
            crossz = xdif(1)*ydif(2)-xdif(2)*ydif(1)
            gcross = sqrt(crossx**2+crossy**2+crossz**2)
            !
            !     COMPUTE THE THREE TRIPLE PRODUCTS--WITH R1,R2 AND R3 (TOP)
            !     AND THEN COMPUTE THE ANGLE GAMMA USED IN DELTA.
            !
            do 20 i = 1, 3
                xmagr(i)        = sqrt(rv(i,1)**2+rv(i,2)**2+rv(i,3)**2)
                top(i)          = crossx*rv(i,1)+crossy*rv(i,2)+crossz*rv(i,3)
                testcm_gamma(i) = +top(i)/(gcross*xmagr(i))
            20 continue
            !
            !     NOW THE SECOND TERM IN DELTA--PHINOT
            !     FIRST THE SWITCHING FUNCTION SPHI THEN PHINOT
            !
            sphi = 0.0d0
            rest = pdt2cm_bphi*(rstar-pdt2cm_cphi)**3
            if( rest .lt. 70.0d0 )then
                erest = exp(rest)
                sphi  = 1.0d0-tanh(pdt2cm_aphi*(rstar-pdt2cm_rnot)*erest)
            endif
            phinot = pdt2cm_pch4+(pdt2cm_pch4-pi/2.0d0)*(sphi-1.0d0)
            !
            !     NOW PUT THESE TWO TERMS TOGETHER TO GET EACH DELTA
            !
            do 30 i = 1, 3
                delta(i) = acos(testcm_gamma(i))-phinot
            30 continue
            !
            !     NEXT COMPUTE THE PARTIALS OF DELTA W/R TO EACH X(I)
            !      IN DOING THIS COMPUTATION FIRST FOR EACH TERM CALCULATED
            !      ABOVE THE PARTIAL WILL BE COMPUTED AND A STORED USING THE
            !      SAME VARIABLE WITH A 'D' IN FRONT.
            !
            !
            !     START WITH THE PARTIAL OF CROSS WITH RESPECT TO EACH X(I)
            !
            do 60 i = 1, 3
                do 40 j = 1, 3
                    dcross(i,j) = 0.0d0
                40 continue
                do 50 j = 13, 18
                    dcross(i,j) = 0.0d0
                50 continue
                if( i .eq. 1 )then
                    dcross(i,4)  = 0.0d0
                    dcross(i,5)  = zdif(1)-zdif(2)
                    dcross(i,6)  = ydif(2)-ydif(1)
                    dcross(i,7)  = 0.0d0
                    dcross(i,8)  = zdif(2)
                    dcross(i,9)  = -ydif(2)
                    dcross(i,10) = 0.0d0
                    dcross(i,11) = -zdif(1)
                    dcross(i,12) = ydif(1)
                else if( i .eq. 2 )then
                    dcross(i,4)  = zdif(2)-zdif(1)
                    dcross(i,5)  = 0.0d0
                    dcross(i,6)  = xdif(1)-xdif(2)
                    dcross(i,7)  = -zdif(2)
                    dcross(i,8)  = 0.0d0
                    dcross(i,9)  = xdif(2)
                    dcross(i,10) = zdif(1)
                    dcross(i,11) = 0.0d0
                    dcross(i,12) = -xdif(1)
                else if( i .eq. 3 )then
                    dcross(i,4)  = ydif(1)-ydif(2)
                    dcross(i,5)  = xdif(2)-xdif(1)
                    dcross(i,6)  = 0.0d0
                    dcross(i,7)  = ydif(2)
                    dcross(i,8)  = -xdif(2)
                    dcross(i,9)  = 0.0d0
                    dcross(i,10) = -ydif(1)
                    dcross(i,11) = xdif(1)
                    dcross(i,12) = 0.0d0
                endif
            60 continue
            !
            !     NOW COMPUTE THE PARTIALS FOR THE THREE R VECTORS
            !
            do 70 i = 1, 3
                do 70 j = 1, 3
                    !            DO 70 K = 1, 3
                    ! BUGFIX:  DRV was not fully initialized to zero causing trouble on some
                    !          platforms.
                    do 70 k = 1, 18
                        drv(i,j,k) = 0.0d0
            70 continue
            do 80 k = 4, 6
                drv(1,k-3,k) = 1.0d0
            80 continue
            do 90 k = 7, 9
                drv(2,k-6,k) = 1.0d0
            90 continue
            do 100 k = 10, 12
                drv(3,k-9,k) = 1.0d0
            100 continue
            do 110 i = 1, 3
                do 110 j = 1, 3
                    drv(i,j,j+12) = -1.0d0
            110 continue
            !
            !     NOW COMPUTE THE PARTIALS FOR THE TRIPLE PRODUCT  (TOP)
            !
            do 120 i = 1, 18
                do 120 j = 1, 3
                    dtop(j,i) = dcross(1,i)*rv(j,1)+crossx*drv(j,1,i)+dcross(2,i )*rv(j,2)+crossy*drv(j,2,i)+ &
                        & dcross(3,i)* rv(j,3)+crossz* drv(j,3,i)
            120 continue
            !
            !     PARTIALS FOR THE MAGNITUDE OF THE CROSS PRODUCT  (GCROSS)
            !
            do 130 i = 1, 18
                dgcros(i) = (1.0d0/gcross)*(dcross(1,i)*crossx+dcross(2,i)* crossy+dcross(3,i)*crossz)
            130 continue
            !
            !     PARTIALS FOR THE MAGNITUDE OF R1, R2 AND R3
            !
            do 150 i = 1, 3
                do 150 j = 1, 18
                    dxmagr(i,j) = 0.0d0
                    do 140 k = 1, 3
                        dxmagr(i,j) = dxmagr(i,j)+(1.0d0/xmagr(i))*(rv(i,k)*drv(i,k,j))
                    140 continue
            150 continue
            !
            !     NOW WE'RE READY FOR THE PARTIALS OF GAMMA
            !
            do 160 i = 1, 3
                do 160 j = 1, 18
                    testcm_dgamma(i,j) = +gcross*xmagr(i)*dtop(i,j)-top(i)*(xmagr(i)* dgcros(j)+gcross*dxmagr(i,j))
                    testcm_dgamma(i,j) = testcm_dgamma(i,j)/(gcross*xmagr(i))**2
                    testcm_dgamma(i,j) = testcm_dgamma(i,j)*(-1.0d0)/(sqrt(1.0d0-testcm_gamma(i)**2) )
            160 continue
            !
            !     NOW THE PARTIAL FOR THE SECOND TERM IN DELTA -- PHINOT
            !       FIRST COMPUTE THE PARTIAL OF S(SWITCHING FUNCTION) WITH
            !       RESPECT TO RSTAR, THE THE PARTIAL OF RSTAR W/R X(I) AND
            !       FINALLY THE PARTIAL OF PHINOT W/R TO S.
            !
            dsdr = 0.0d0
            rest = pdt2cm_bphi*(rstar-pdt2cm_cphi)**3
            if( rest .lt. 70.0d0 )then
                texp = exp(rest)
                dsdr = -((sech(pdt2cm_aphi*(rstar-pdt2cm_rnot)*texp))**2)
                dsdr = dsdr*pdt2cm_aphi*texp*(1.0d0+(rstar-pdt2cm_rnot)*pdt2cm_bphi*3.0d0*(rstar- pdt2cm_cphi)**2)
            endif
            do 170 i = 1, 18
                drdx(i) = 0.0d0
                dphi(i) = 0.0d0
            170 continue
            !
            cst = 1.0d0/(2.0d0*rstar)
            do 180 i = 1, 3
                drdx(i)    = cst*2.0d0*(x(i)-x(i+12))
                drdx(i+12) = -drdx(i)
                dphi(i)    = drdx(i)*dsdr*(pdt2cm_pch4-pi/2.0d0)
                dphi(i+12) = -dphi(i)
            180 continue
            !
            !     NOW PUT BOTH TERMS TOGETHER TO GET THE PARTIALS OF DELTA
            !
            do 190 i = 1, 3
                do 190 j = 1, 18
                    ddelta(i,j) = testcm_dgamma(i,j)-dphi(j)
            190 continue
            !
            !     WE HAVE NOW COMPUTED ALL THE DELTA TERMS.  NEXT COMPUTE THE
            !       FORCE CONSTANTS (F AND H) AND THEIR DERIVATIVES.  THE
            !       OUT-OF-PLANE BEND ENERGY WILL ALSO BE COMPUTED HERE (VOOP)
            !
            !     FIRST THE SWITCHING FUNCTION
            !
            s3   = 1.0d0-tanh(pdt2cm_a3*(rstar-pdt2cm_rnot)*(rstar-pdt2cm_b3)**2)
            f    = (1.0d0-s3)*pdt2cm_fch3
            h    = (1.0d0-s3)*pdt2cm_hch3
            voop = 0.0d0
            do 200 i = 1, 3
                voop = voop+f*(delta(i)**2)+h*(delta(i)**4)
            200 continue
            !
            !     PARTIAL OF THE SWITCHING FUNCTION W/R TO RSTAR
            !
            dsdr = -(sech(pdt2cm_a3*(rstar-pdt2cm_rnot)*(rstar-pdt2cm_b3)**2)**2)*pdt2cm_a3*((rstar-pdt2cm_b3)**2 + &
                & (rstar-pdt2cm_rnot)*2*(rstar-pdt2cm_b3))
            dfds = -pdt2cm_fch3
            dhds = -pdt2cm_hch3
            !
            !     PARTIAL OF EACH FORCE CONSTANT
            !
            do 210 i = 1, 18
                dfdx(i) = 0.0d0
                dhdx(i) = 0.0d0
            210 continue
            do 220 i = 1, 3
                dfdx(i)    = dfds*dsdr*drdx(i)
                dhdx(i)    = dhds*dsdr*drdx(i)
                dfdx(i+12) = -dfdx(i)
                dhdx(i+12) = -dhdx(i)
            220 continue
            !
            !    NOW WE ARE READY FOR THE FULL PARTIAL DERIVATIVES OF VOOP
            !
            do 250 i = 1, 18
                term1 = 0.0d0
                do 230 j = 1, 3
                    term1 = term1+dfdx(i)*delta(j)**2
                230 continue
                term2 = 2.0d0*delta(1)*ddelta(1,i)+2.0d0*delta(2)*ddelta(2,i) +2.0d0*delta(3)*ddelta(3,i)
                term2 = term2*f
                term3 = 0.0d0
                do 240 j = 1, 3
                    term3 = term3+dhdx(i)*delta(j)**4
                240 continue
                term4    = 4.0d0*delta(1)**3*ddelta(1,i)+4.0d0*delta(2)**3* ddelta(2,i)+4.0d0*delta(3)**3*ddelta(3,i)
                term4    = term4*h
                dvoop(i) = term1+term2+term3+term4
            250 continue
        end subroutine oopb
        !> @brief Evaluate SECH for the CH5_GEN_J1_1987 potential.
        function sech ( x )
            implicit none
            double precision :: sech, x
            if( x .gt. 70.0d0 )then
                sech = 0.0d0
            else
                sech = 1.0d0/cosh(x)
            endif
        end function sech
        !> @brief Evaluate SWITCH for the CH5_GEN_J1_1987 potential.
        subroutine switch ( kmax, lmg, v, coord, dxx )
            !
            ! THIS ROUTINE SWITCHES MATRICES R,R2,DER TO ALLOW FOR DIFFERENT H ATOMS
            !      BEING OFF THE CH4 MOLECULE.
            ! LMG IS A CONTROL PARAMETER
            ! KMAX IS THE LABEL OF THE H ATOM WITH MAXIMUM C-H DISTANCE.
            ! IF LMG=1, R AND R2 MATRICES ARE SWITCHED APPROPRIATELY.
            ! IF LMG=2, R AND R2 MATRICES ARE SWITCHED BACK AND MATRIX DER IS SWITCH
            !
            implicit none
            double precision :: coord, dxx, store, v
            integer          :: k, kk, kmax, l, ll, lmg, m, mm, n, nn
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
            dimension coord(18), dxx(18)
            !
            !
            !
            go to (10, 20, 30, 40, 50), kmax
            !
            ! H1 IS MAXIMUM
            !
            10 k = 1
            kk = 5
            l  = 7
            ll = 10
            m  = 8
            mm = 11
            n  = 9
            nn = 12
            go to 60
            !
            ! H2 IS MAXIMUM
            !
            20 k = 2
            kk = 5
            l  = 6
            ll = 10
            m  = 8
            mm = 13
            n  = 9
            nn = 14
            go to 60
            !
            ! H3 IS MAXIMUM
            !
            30 k = 3
            kk = 5
            l  = 6
            ll = 11
            m  = 7
            mm = 13
            n  = 9
            nn = 15
            go to 60
            !
            ! H4 IS MAXIMUM
            !
            40 k = 4
            kk = 5
            l  = 6
            ll = 12
            m  = 7
            mm = 14
            n  = 8
            nn = 15
            go to 60
            !
            ! H5 IS MAXIMUM
            ! NO SWITCHES ARE REQUIRED
            !
            50 go to 80
            !
            ! SWITCH R AND R2 MATRICES
            !
            60 store = icorcm_ar(k)
            icorcm_ar(k)  = icorcm_ar(kk)
            icorcm_ar(kk) = store
            store         = icorcm_ar(l)
            icorcm_ar(l)  = icorcm_ar(ll)
            icorcm_ar(ll) = store
            store         = icorcm_ar(m)
            icorcm_ar(m)  = icorcm_ar(mm)
            icorcm_ar(mm) = store
            store         = icorcm_ar(n)
            icorcm_ar(n)  = icorcm_ar(nn)
            icorcm_ar(nn) = store
            potrcm_r2(k)  = icorcm_ar(k)*icorcm_ar(k)
            potrcm_r2(kk) = icorcm_ar(kk)*icorcm_ar(kk)
            potrcm_r2(l)  = icorcm_ar(l)*icorcm_ar(l)
            potrcm_r2(ll) = icorcm_ar(ll)*icorcm_ar(ll)
            potrcm_r2(m)  = icorcm_ar(m)*icorcm_ar(m)
            potrcm_r2(mm) = icorcm_ar(mm)*icorcm_ar(mm)
            potrcm_r2(n)  = icorcm_ar(n)*icorcm_ar(n)
            potrcm_r2(nn) = icorcm_ar(nn)*icorcm_ar(nn)
            if( lmg-1 ) 80, 80, 70
            !
            ! SWITCH DER MATRIX *******
            !
            70 store = icorcm_der(k)
            icorcm_der(k)  = icorcm_der(kk)
            icorcm_der(kk) = store
            store          = icorcm_der(l)
            icorcm_der(l)  = icorcm_der(ll)
            icorcm_der(ll) = store
            store          = icorcm_der(m)
            icorcm_der(m)  = icorcm_der(mm)
            icorcm_der(mm) = store
            store          = icorcm_der(n)
            icorcm_der(n)  = icorcm_der(nn)
            icorcm_der(nn) = store
            80 return
        end subroutine switch
        !> @brief Evaluate THETA for the CH5_GEN_J1_1987 potential.
        subroutine theta ( asf, dasf, asb, dasb )
            !
            !     THIS SUBROUTINE COMPUTES THE EQUILIBRIUM THETA AS A FUNCTION OF
            !     R (THE AVERAGE C-H BOND LENGTH IN METHANE).
            !     THE DERIVATIVE OF THETA WITH RESPECT TO AR(I) IS ALSO COMPUTED.
            !     NOTE HERE THAT BECAUSE OF SYMMETRY EACH OF THESE PARTIALS IS
            !     THE SAME.       ADDED 11/22/85
            !
            implicit none
            double precision :: a, asb, asf, b, c, ctau, ctau2, dasb, dasf, delta2, f, fk1, fk2, ravg, rdif, rdif2, rdif3, rdif4, rdif5, xyz
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
            ravg = (icorcm_ar(1)+icorcm_ar(2)+icorcm_ar(3)+icorcm_ar(4))/4.0d0
            !
            !     COMPUTE THETA
            !
            if( ravg .gt. pdt2cm_dij )then
                asf  = pi/2.0d0
                dasf = 0.0d0
            else if( ravg .le. pdatcm_re(1) )then
                asf  = pdatcm_tau
                dasf = 0.0d0
            else if( ravg .gt. pdatcm_re(1).and.ravg .le. pdatcm_re(1)+pdt2cm_delta )then
                delta2 = pdt2cm_delta*pdt2cm_delta
                fk1    = -4.0d0*pdatcm_cp/(delta2)
                a      = 3.0d0*fk1/(delta2)
                b      = -(5.0d0*a*delta2+fk1)/(2.0d0*pdt2cm_delta)
                c      = -(20.0d0*a*delta2+12.0d0*b*pdt2cm_delta)/6.0d0
                f      = pdatcm_tau
                rdif   = ravg-pdatcm_re(1)
                rdif2  = rdif*rdif
                rdif3  = rdif2*rdif
                rdif4  = rdif3*rdif
                rdif5  = rdif4*rdif
                asf    = a*rdif5+b*rdif4+c*rdif3+f
                dasf   = 5.0d0*a*rdif4+4.0d0*b*rdif3+3.0d0*c*rdif2
            else if( ravg .le. pdt2cm_dij-pdt2cm_delta )then
                asf  = pdatcm_tau-4.0d0*pdatcm_cp*(ravg-pdatcm_re(1))
                dasf = -4.0d0*pdatcm_cp
            else if( ravg .le. pdt2cm_dij.and.ravg .gt. pdt2cm_dij-pdt2cm_delta )then
                fk1   = (pdatcm_tau-4.0d0*pdatcm_cp*(pdt2cm_dij-pdt2cm_delta-pdatcm_re(1))-f)/(pdt2cm_delta**3)
                fk2   = (-4.0d0*pdatcm_cp)/(pdt2cm_delta**2)
                a     = -3.0d0*(2.0d0*fk1+fk2)/(pdt2cm_delta**2)
                b     = (fk2+5.0d0*a*pdt2cm_delta**2)/(2.0d0*pdt2cm_delta)
                c     = (12.0d0*b*pdt2cm_delta-20.0d0*a*pdt2cm_delta**2)/6.0d0
                f     = pi/2.0d0
                rdif  = ravg-pdt2cm_dij
                rdif3 = rdif*rdif*rdif
                rdif4 = rdif3*rdif
                rdif5 = rdif4*rdif
                asf   = a*rdif5+b*rdif4+c*rdif3+f
                dasf  = 5.0d0*a*rdif4+4.0d0*b*rdif3+3.0d0*c*rdif*rdif
            endif
            !
            ! COMPUTE BACKSIDE EQUILIBRIUM ANGLE
            !
            xyz   = cos(asf)
            ctau  = cos(pdatcm_tau)
            ctau2 = ctau*ctau
            a     = (ctau+0.5d0)/ctau2
            asb   = acos(a*xyz*xyz-0.5d0)
            dasb  = 2.0d0*a*xyz*sin(asf)*dasf/sin(asb)
            !
            ! NOW COMPUTE THE DER W/R TO RI INSTEAD OF RAVG
            !
            dasf = dasf/4.0d0
            dasb = dasb/4.0d0
        end subroutine theta
        !> @brief Evaluate TRI for the CH5_GEN_J1_1987 potential.
        subroutine tri ( rr1, rr2, rr3, dr1, dr2, dr3, i, j, k, e, dtot )
            !
            ! ROUTINE TO CALCULATE THREE-BODY ENERGY FOR ABC SYSTEM
            ! RR1,RR2,RR3 ARE THE THREE INTERPARTICLE DISTANCES
            ! DR1,DR2,DR3 ARE THE THREE DERIVATIVES WITH RESPECT TO RR1,RR2, AND RR3
            ! I,J, AND K GIVE THE ARRAY NUMBER FOR THE APPROPRIATE POTENTIAL PARAMET
            ! E IS THE FINAL THREE-BODY ENERGY
            !
            implicit none
            double precision :: a, a1, a2, aa1, aa2, alab, alac, albc, b, b1, b2, bb, bb1, bb2, c, c1, c2, cc1, cc2, dab1, dab3, dac1, dac3, de1ab, de1ac, de1bc, de3ab, de3ac, de3bc, dqr, dr1, dr2, dr3, dtot, dur, dux, duy, duz, dxr, dyr, dzr, e, e1ab, e1ac, e1bc, e3ab, e3ac, e3bc, qab, qac, qbc, rr1, rr2, rr3, u, xx, yy, zz
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
            a     = rr1-pdatcm_re(i)
            b     = rr2-pdatcm_re(j)
            c     = rr3-pdatcm_re(k)
            a1    = exp(-pdatcm_alph(i)*a)
            b1    = exp(-pdatcm_alph(j)*b)
            c1    = exp(-pdatcm_alph(k)*c)
            a2    = a1*a1
            b2    = b1*b1
            c2    = c1*c1
            e1ab  = pdatcm_d1(i)*(a2-2.0d0*a1)
            e1ac  = pdatcm_d1(j)*(b2-2.0d0*b1)
            e1bc  = pdatcm_d1(k)*(c2-2.0d0*c1)
            de1ab = 2.0d0*pdatcm_d1(i)*(a1-a2)*pdatcm_alph(i)
            dab1  = 2.0d0*pdatcm_d1(i)*(a1-a2)*(a*alpdcm_dalp)
            de1ac = 2.0d0*pdatcm_d1(j)*(b1-b2)*pdatcm_alph(j)
            dac1  = 2.0d0*pdatcm_d1(j)*(b1-b2)*(b*alpdcm_dalp)
            de1bc = 2.0d0*pdatcm_alph(k)*pdatcm_d1(k)*(c1-c2)
            aa1   = exp(-pdatcm_beta(i)*a)
            aa2   = aa1*aa1
            e3ab  = pdatcm_d3(i)*(aa2+2.0d0*aa1)
            de3ab = -2.0d0*pdatcm_d3(i)*(aa1+aa2)*pdatcm_beta(i)
            dab3  = -2.0d0*pdatcm_d3(i)*(aa1+aa2)*a*alpdcm_dalp
            bb1   = exp(-pdatcm_beta(j)*b)
            bb2   = bb1*bb1
            e3ac  = pdatcm_d3(j)*(bb2+2.0d0*bb1)
            de3ac = -2.0d0*pdatcm_d3(j)*(bb1+bb2)*pdatcm_beta(j)
            dac3  = -2.0d0*pdatcm_d3(j)*(bb1+bb2)*b*alpdcm_dalp
            cc1   = exp(-pdatcm_beta(k)*c)
            cc2   = cc1*cc1
            e3bc  = pdatcm_d3(k)*(cc2+2.0d0*cc1)
            de3bc = -2.0d0*pdatcm_beta(k)*pdatcm_d3(k)*(cc1+cc2)
            qab   = (e1ab+e3ab)/2.0d0
            qac   = (e1ac+e3ac)/2.0d0
            qbc   = (e1bc+e3bc)/2.0d0
            alab  = (e1ab-e3ab)/2.0d0
            alac  = (e1ac-e3ac)/2.0d0
            albc  = (e1bc-e3bc)/2.0d0
            xx    = alab-albc
            yy    = albc-alac
            zz    = alac-alab
            u     = (xx*xx+yy*yy+zz*zz)/2.0d0
            u     = sqrt(u)
            e     = qab+qbc+qac-u
            bb    = 4.0d0*u
            dxr   = (dab1-dab3)/2.0d0
            dyr   = (dac3-dac1)/2.0d0
            dzr   = -(dxr+dyr)
            dux   = xx/(2.0d0*u)
            duy   = yy/(2.0d0*u)
            duz   = zz/(2.0d0*u)
            dur   = dux*dxr+duy*dyr+duz*dzr
            dqr   = (dab1+dac1+dab3+dac3)/2.0d0
            dtot  = dqr-dur
            dr1   = (de1ab+de3ab)/2.0d0-(2.0d0*alab-albc-alac)*(de1ab-de3ab)/bb
            dr2   = (de1ac+de3ac)/2.0d0-(2.0d0*alac-alab-albc)*(de1ac-de3ac)/bb
            dr3   = (de1bc+de3bc)/2.0d0-(2.0d0*albc-alac-alab)*(de1bc-de3bc)/bb
        end subroutine tri
        !> @brief Evaluate TRIPLE for the CH5_GEN_J1_1987 potential.
        subroutine triple
            !
            !     THIS SUBROUTINE COMPUTES THE C-H' SATO PARAMETER WHICH IS A
            !     FUNCTION OF AR(H'-H) AND THE C-H-H' BEND ANGLE.
            !
            implicit none
            double precision :: alpha, c1, check, dalr1, dalr5, dalr6, dc1r6, dcyx, deg20, dyxr, r6, tem, tem2, tem22, tembot, u, ui, yx
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
            !     COMPUTE C1--IT IS A FUNCTION OF R6--THE H-H' BOND DISTANCE
            !
            if( icorcm_ar(6) .gt. 15.0d0 )then
                r6 = icorcm_ar(6)/3.0d0
                yx = (exp(2.0d0*pi*(r6-pdt2cm_rex)/pdt2cm_d))**3
            else
                yx = exp(2.0d0*pi*(icorcm_ar(6)-pdt2cm_rex)/pdt2cm_d)
            endif
            c1    = -(pdt2cm_a*yx/(yx+1.0d0)+pdt2cm_b/yx)-pdt2cm_c0
            dyxr  = (2.0d0*pi/pdt2cm_d)*yx
            dcyx  = pdt2cm_b/yx-(pdt2cm_a*yx/(1+yx))/(1+yx)
            dc1r6 = dcyx*dyxr/yx
            !
            !     COMPUTE THE C-H4-H5 ANGLE
            !
            tem    = icorcm_ar(5)**2-icorcm_ar(1)**2-icorcm_ar(6)**2
            tembot = -2.0d0*icorcm_ar(1)*icorcm_ar(6)
            check  = tem/tembot
            if( check .le. -1.0d0 )then
                alpha         = 0.0d0
                pdatcm_cc(3)  = pdt2cm_c0
                dervcp_deccr1 = 0.0d0
                dervcp_deccr5 = 0.0d0
                dervcp_deccr6 = 0.0d0
                return
            else
                alpha = acos(tem/(-2.0d0*icorcm_ar(1)*icorcm_ar(6)))
                alpha = pi-alpha
            endif
            !
            !     NOW COMPUTE THE ANGLE DERIVATIVES (DERIV. ALPHA W.R. R1,R5,R6)
            !
            u     = tem/(-2.0d0*icorcm_ar(1)*icorcm_ar(6))
            ui    = -1.0d0/(sqrt(1.0d0-u**2))
            dalr5 = -ui*icorcm_ar(5)/(icorcm_ar(1)*icorcm_ar(6))
            tem2  = 1.0d0/(-2.0d0*icorcm_ar(1)*icorcm_ar(6))
            tem22 = tem2*tem2
            dalr1 = -2.0d0*ui*(icorcm_ar(1)*tem2-icorcm_ar(6)*tem*tem22)
            dalr6 = -2.0d0*ui*(icorcm_ar(6)*tem2-icorcm_ar(1)*tem*tem22)
            !
            !     COMBINE THE RADIAL AND ANGLE COMPONENT TO COMPUTE THE TRIPLET
            !     PARAMETER CC(3) AND THE DERIVATIVES.  HERE TEM IS JUST THE PARTIAL
            !     OF CC W.R.T. ALPHA.
            !
            deg20 = pi/9.0d0
            if( alpha .le. deg20 )then
                pdatcm_cc(3)  = pdt2cm_c0+c1*(sin(9.0d0*alpha/2.0d0))**4
                tem           = 18.0d0*(sin(9.0d0*alpha/2.0d0))**3*cos(9.0d0*alpha/ 2.0d0)
                dervcp_deccr1 = c1*tem*dalr1
                dervcp_deccr5 = c1*tem*dalr5
                dervcp_deccr6 = dc1r6*(sin(9.0d0*alpha/2.0d0))**4+c1*tem*dalr6
            else
                pdatcm_cc(3)  = pdt2cm_c0+c1
                dervcp_deccr1 = 0.0d0
                dervcp_deccr5 = 0.0d0
                dervcp_deccr6 = dc1r6
            endif
        end subroutine triple
        !> @brief Evaluate prepare_local_context for the CH5_GEN_J1_1987 potential.
        subroutine prepare_local_context
            !
            !   N3TMMN = 3 * NATOMS
            !   NATOMS = the number of atoms represented by this potential function
            !
            !   The variable N3TMMN is the minimum value of N3TM allowed to be
            !   passed by the calling routine for the number of cartesian
            !   coordinates needed to represent the full system represented by this
            !   potential energy surface routine.
            !   N3TM must be greater than or equal to N3TMMN.
            !
            implicit none
            double precision :: coord, df, dfr, dk, dvoop, dxx, evcon, f, y1, y2
            integer          :: i
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
            dimension dfr(4), dvoop(18), y1(6), y2(6), f(4), df(4), dk(6)
            !
            dimension coord(18), dxx(18)
            !
            !
            !
            data evcon / 2.721161d+1 /
            !
            !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
            !
            if( infocm_natoms .gt. 25 )then
                write(usricm_nflag(18), '(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
                stop
            endif
            !
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !
            infocm_ref(1) = 'T. Joseph, R. Steckler, and D. G. Truhlar'
            infocm_ref(2) = 'J. Chem. Phys. Vol. 87, p. 7036, 1987'
            !
            infocm_indexes(1) = 1
            infocm_indexes(2) = 1
            infocm_indexes(3) = 1
            infocm_indexes(4) = 1
            infocm_indexes(5) = 6
            infocm_indexes(6) = 1
            !
            !
            !
            infocm_irctnt = 6
            !
            !      CALL POTINFO
            !
            call ancvrt
            !
        end subroutine prepare_local_context
        !> @brief Evaluate POTINFO for the CH5_GEN_J1_1987 potential.
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
        !> @brief Evaluate ANCVRT for the CH5_GEN_J1_1987 potential.
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
        !> @brief Evaluate CARTOU for the CH5_GEN_J1_1987 potential.
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
        !> @brief Evaluate CARTTOR for the CH5_GEN_J1_1987 potential.
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
        !> @brief Evaluate EUNITZERO for the CH5_GEN_J1_1987 potential.
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
        !> @brief Evaluate RTOCART for the CH5_GEN_J1_1987 potential.
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
                            go to 20
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
        !> @brief Evaluate DEDCOU for the CH5_GEN_J1_1987 potential.
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
end module chempotpy_ch5_gen_j1_1987_pes

!> @brief Preserve the historical Python-facing pes entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_ch5_gen_j1_1987_pes, only: evaluate
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
