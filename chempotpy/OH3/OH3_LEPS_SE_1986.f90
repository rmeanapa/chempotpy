!   System:                        OH3
!   Functional form:
!   Common name:
!   Number of derivatives:         1
!   Number of bodies:              4
!   Number of electronic surfaces: 1
!   Interface:                     4-XS
!   Reference: G. C. Schatz and H. Elgersma
!              Chem. Phys. Lett. 73, 21 (1980).

!> @brief Implement OH3_LEPS_SE_1986, a potential energy surface for OH3.
module chempotpy_oh3_leps_se_1986_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the OH3_LEPS_SE_1986 potential energy surface for OH3.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_units, only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 4
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter    :: n3atom  = 75
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: concom_acon(2)
        real(wp)              :: concom_alp(4)
        real(wp)              :: concom_beta(3)
        real(wp)              :: concom_clam(4)
        real(wp)              :: concom_con(7)
        real(wp)              :: concom_de(3)
        real(wp)              :: concom_dum(22)
        real(wp)              :: concom_dum1(10)
        real(wp)              :: concom_dum2(10)
        real(wp)              :: concom_gam(3)
        real(wp)              :: concom_re(3)
        real(wp)              :: concom_rehh
        real(wp)              :: concom_reoh
        real(wp)              :: concom_sato
        real(wp)              :: infocm_cartnu(natom,3)
        integer               :: infocm_icartr
        integer               :: infocm_indexes(natom)
        integer               :: infocm_irctnt
        integer               :: infocm_mder
        integer               :: infocm_msurf
        integer               :: infocm_natoms
        character             :: infocm_ref(5)
        real(wp)              :: pot2cm_dedr(4)
        real(wp)              :: pot2cm_energy
        real(wp)              :: pot2cm_r(4)
        real(wp)              :: potcm_dvdr(6)
        real(wp)              :: potcm_r(6)
        real(wp)              :: potcm_vtot
        real(wp)              :: pt1cm_degsdr(n3atom)
        real(wp)              :: pt1cm_engygs
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: pt4cm_deesdr(n3atom,isurf)
        real(wp)              :: pt4cm_engyes(isurf)
        real(wp)              :: pt5cm_deijdr(n3atom,jsurf)
        real(wp)              :: pt5cm_engyij(jsurf)
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
        real(wp)              :: utilcm_cnvrtd
        real(wp)              :: utilcm_cnvrtde
        real(wp)              :: utilcm_cnvrte
        real(wp)              :: utilcm_descartnu(natom,3,isurf)
        real(wp)              :: utilcm_dgscartnu(natom,3)
        real(wp)              :: utilcm_dijcartnu(natom,3,jsurf)
        integer               :: utilcm_ireorder
        integer               :: utilcm_kediag
        integer               :: utilcm_keoffd
        integer               :: utilcm_ksdiag
        integer               :: utilcm_ksoffd
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        concom_acon      = 0.e0_wp; concom_alp = 0.e0_wp; concom_beta = 0.e0_wp; concom_clam = 0.e0_wp; concom_con = 0.e0_wp
        concom_de        = 0.e0_wp; concom_dum = 0.e0_wp; concom_dum1 = 0.e0_wp; concom_dum2 = 0.e0_wp; concom_gam = 0.e0_wp; concom_re = 0.e0_wp
        concom_rehh      = 0.e0_wp; concom_reoh = 0.e0_wp; concom_sato = 0.e0_wp; infocm_cartnu = 0.e0_wp; infocm_icartr = 0
        infocm_indexes   = 0; infocm_irctnt = 0; infocm_mder = 0; infocm_msurf = 0; infocm_natoms = 0
        infocm_ref       = ' '
        pot2cm_dedr      = 0.e0_wp
        pot2cm_energy    = 0.e0_wp
        pot2cm_r         = 0.e0_wp
        potcm_dvdr       = 0.e0_wp
        potcm_r          = 0.e0_wp
        potcm_vtot       = 0.e0_wp
        pt1cm_degsdr     = 0.e0_wp
        pt1cm_engygs     = 0.e0_wp
        pt1cm_r          = 0.e0_wp
        pt3cm_ezero      = 0.e0_wp
        pt4cm_deesdr     = 0.e0_wp
        pt4cm_engyes     = 0.e0_wp
        pt5cm_deijdr     = 0.e0_wp
        pt5cm_engyij     = 0.e0_wp
        usricm_anuzero   = 0.e0_wp
        usricm_cart      = 0.e0_wp
        usricm_nasurf    = 0
        usricm_nder      = 0
        usricm_nflag     = 0
        usricm_nulbl     = 0
        usrocm_descart   = 0.e0_wp
        usrocm_dgscart   = 0.e0_wp
        usrocm_dijcart   = 0.e0_wp
        usrocm_pengyes   = 0.e0_wp
        usrocm_pengygs   = 0.e0_wp
        usrocm_pengyij   = 0.e0_wp
        utilcm_cnvrtd    = 0.e0_wp
        utilcm_cnvrtde   = 0.e0_wp
        utilcm_cnvrte    = 0.e0_wp
        utilcm_descartnu = 0.e0_wp
        utilcm_dgscartnu = 0.e0_wp
        utilcm_dijcartnu = 0.e0_wp
        utilcm_ireorder  = 0
        utilcm_kediag    = 0
        utilcm_keoffd    = 0
        utilcm_ksdiag    = 0
        utilcm_ksoffd    = 0
        usricm_nasurf = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0], shape(usricm_nasurf))
        usricm_nder    = 1
        usricm_nflag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0]
        usricm_anuzero = 0.0e0_wp
        infocm_icartr  = 2
        infocm_msurf   = 0
        infocm_mder    = 1
        usricm_nulbl   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        infocm_natoms  = 4
        concom_de      = [0.148201e0_wp, 0.0275690e0_wp, 0.151548e0_wp]
        concom_beta    = [1.260580e0_wp, 0.924180e0_wp, 1.068620e0_wp]
        concom_re      = [1.863300e0_wp, 2.907700e0_wp, 1.428600e0_wp]
        concom_sato    = 0.10e0_wp
        concom_gam     = [2.399700e0_wp, 1.058350e0_wp, 2.399700e0_wp]
        concom_reoh    = 1.808090e0_wp
        concom_rehh    = 2.861590e0_wp
        concom_con     = [-.0015920e0_wp, 0.026963e0_wp, 0.0014689e0_wp, 0.080011e0_wp, 0.085816e0_wp, -0.063179e0_wp, &
            & 0.101380e0_wp]
        concom_alp    = [4.773e0_wp, 7.14e0_wp, 2.938e0_wp, 5.28e0_wp]
        concom_clam   = [0.10e0_wp, 0.10e0_wp, 0.20e0_wp, 0.03e0_wp]
        concom_acon   = [0.10e0_wp, 0.009e0_wp]
        v             = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        usricm_nder   = igrad
        usricm_nasurf = 1
        usricm_cart   = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                usricm_cart(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = usrocm_pengygs*hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = usrocm_dgscart(iatom,idir)*hartree_per_bohr_to_ev_per_angstrom
                enddo
            enddo
        elseif( igrad == 2 )then
            write (*, *) 'Only energy and gradient are available'
        endif

        contains

        !> @brief Evaluate prepare_local_context for the OH3_LEPS_SE_1986 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            infocm_ref(1)     = 'OH + H2 potential energy function'
            infocm_ref(2)     = 'G. C. Schatz and H. Elgersma'
            infocm_ref(3)     = 'Chem. Phys. Lett. 73, 21 (1980)'
            infocm_irctnt     = 3
            infocm_indexes(1) = 8
            infocm_indexes(2) = 1
            infocm_indexes(3) = 1
            infocm_indexes(4) = 1
            call ancvrt
        end subroutine prepare_local_context

        !> @brief Evaluate pot for the OH3_LEPS_SE_1986 potential.
        subroutine pot
            !   This subprogram evaluates the OH + H2 potential energy and derivativ
            !   of the potential with respect to the internal coordinates.
            !   The subprogram PREPEF must be called once before any calls to this
            !   subprogram.
            !   All calculations in this subprogram are in hartree atomic units.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: b, d, dvmor, rr, t, vmor
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 1
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !   Define the statement functions VMOR and DVMOR which evaluate the Mor
            !   diatomic potential energy  and the derivatives of the Morse potentia
            vmor(d, b, t, rr)  = d*(1.0e0_wp-exp(-b*(rr-t)))**2
            dvmor(d, b, t, rr) = 2.0e0_wp*b*d*(1.0e0_wp-exp(-b*(rr-t)))*exp(-b*(rr-t))
            call cartou
            call carttor
            do i = 1, 6
                potcm_r(i) = pt1cm_r(i)
            enddo
            !   Zero the array which will contain the derivatives of the energy
            !   with respect to the internal coordinates.
            if( usricm_nder .eq. 1 )then
                do i = 1, 6
                    potcm_dvdr(i) = 0.0e0_wp
                enddo
            endif
            !   Calculate the Morse part of the potential energy.
            potcm_vtot = vmor(concom_de(1), concom_beta(1), concom_re(1), potcm_r(1))+vmor(concom_de(2), &
                & concom_beta(2), concom_re(2), potcm_r(4 ))+vmor(concom_de(2), concom_beta(2), concom_re(2), &
                & potcm_r(5))
            !   Calculate the derivatives of the Morse part of the potential energy.
            if( usricm_nder .eq. 1 )then
                potcm_dvdr(1) = potcm_dvdr(1)+dvmor(concom_de(1), concom_beta(1), concom_re(1), potcm_r(1))
                potcm_dvdr(4) = potcm_dvdr(4)+dvmor(concom_de(2), concom_beta(2), concom_re(2), potcm_r(4))
                potcm_dvdr(5) = potcm_dvdr(5)+dvmor(concom_de(2), concom_beta(2), concom_re(2), potcm_r(5))
            endif
            !   Initialize the coordinates for the three-body LEPS part of the poten
            pot2cm_r(1) = potcm_r(2)
            pot2cm_r(2) = potcm_r(3)
            pot2cm_r(3) = potcm_r(6)
            !   Calculate the three-body LEPS portion of the potential and update th
            !   energy term.
            call v3pot
            potcm_vtot = potcm_vtot+pot2cm_energy
            !   Update the array containing the derivatives with the LEPS derivative
            if( usricm_nder .eq. 1 )then
                potcm_dvdr(2) = potcm_dvdr(2)+pot2cm_dedr(1)
                potcm_dvdr(3) = potcm_dvdr(3)+pot2cm_dedr(2)
                potcm_dvdr(6) = potcm_dvdr(6)+pot2cm_dedr(3)
            endif
            !   Initialize the coordinates for the H2O part of the potential for H1
            pot2cm_r(1) = potcm_r(1)
            pot2cm_r(2) = potcm_r(2)
            pot2cm_r(3) = potcm_r(4)
            !   Calculate the H2O part of the potential and update the energy term.
            call vh2o
            potcm_vtot = potcm_vtot+pot2cm_energy
            !   Update the array containing the derivatives with the H2O derivatives
            if( usricm_nder .eq. 1 )then
                potcm_dvdr(1) = potcm_dvdr(1)+pot2cm_dedr(1)
                potcm_dvdr(2) = potcm_dvdr(2)+pot2cm_dedr(2)
                potcm_dvdr(4) = potcm_dvdr(4)+pot2cm_dedr(3)
            endif
            !   Initialize the coordinates for the H2O part of the potential for H1
            pot2cm_r(1) = potcm_r(1)
            pot2cm_r(2) = potcm_r(3)
            pot2cm_r(3) = potcm_r(5)
            !   Calculate the H2O part of the potential and update the energy term.
            call vh2o
            potcm_vtot = potcm_vtot+pot2cm_energy
            !   Update the array containing the derivatives with the H2O derivatives
            if( usricm_nder .eq. 1 )then
                potcm_dvdr(1) = potcm_dvdr(1)+pot2cm_dedr(1)
                potcm_dvdr(3) = potcm_dvdr(3)+pot2cm_dedr(2)
                potcm_dvdr(5) = potcm_dvdr(5)+pot2cm_dedr(3)
            endif
            !   Initialize the coordinates for the four-body part of the potential.
            pot2cm_r(1) = potcm_r(2)
            pot2cm_r(2) = potcm_r(3)
            pot2cm_r(3) = potcm_r(4)
            pot2cm_r(4) = potcm_r(5)
            !   Calculate the four-body part of the potential and update the energy
            call v4pot
            potcm_vtot = potcm_vtot+pot2cm_energy
            !   Update the array containing the derivatives.
            if( usricm_nder .eq. 1 )then
                potcm_dvdr(2) = potcm_dvdr(2)+pot2cm_dedr(1)
                potcm_dvdr(3) = potcm_dvdr(3)+pot2cm_dedr(2)
                potcm_dvdr(4) = potcm_dvdr(4)+pot2cm_dedr(3)
                potcm_dvdr(5) = potcm_dvdr(5)+pot2cm_dedr(4)
            endif
            !   Adjust the potential to the correct zero, which corresponds to
            !   OH(R=RE) and H2(R=RE) at infinity.
            potcm_vtot   = potcm_vtot-2.0e0_wp*concom_de(2)+concom_de(3)
            pt1cm_engygs = potcm_vtot
            do i = 1, 6
                pt1cm_degsdr(i) = potcm_dvdr(i)
            enddo
            call eunitzero
            if( usricm_nder .ne. 0 )then
                call rtocart
                if( usricm_nflag(1)+usricm_nflag(2) .ne. 0 ) call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate V3POT for the OH3_LEPS_SE_1986 potential.
        subroutine v3pot
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: b(3), beta(3), coul(3), de(3), do4z(3), exch(3), op3z(3), rad, re(3), s, top3z(3), &
                & tzp3(3), x(3), z(3), zp3(3), zpo(3)
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !   Set up the values of DE, BETA, and RE for the three-body LEPS potent
            de(1)   = concom_de(1)
            beta(1) = concom_beta(1)
            re(1)   = concom_re(1)
            de(2)   = concom_de(1)
            beta(2) = concom_beta(1)
            re(2)   = concom_re(1)
            de(3)   = concom_de(3)
            beta(3) = concom_beta(3)
            re(3)   = concom_re(3)
            !   Compute useful constants.
            do i = 1, 3
                z(i)     = concom_sato
                zpo(i)   = 1.0e0_wp+z(i)
                op3z(i)  = 1.0e0_wp+3.0e0_wp*z(i)
                top3z(i) = 2.0e0_wp*op3z(i)
                zp3(i)   = z(i)+3.0e0_wp
                tzp3(i)  = 2.0e0_wp*zp3(i)
                do4z(i)  = de(i)/4.0e0_wp/zpo(i)
                b(i)     = beta(i)*do4z(i)*2.0e0_wp
            enddo
            !   Initialize the variable used for storing the energy.
            pot2cm_energy = 0.e0_wp
            do i = 1, 3
                x(i)          = exp(-beta(i)*(pot2cm_r(i)-re(i)))
                coul(i)       = do4z(i)*(zp3(i)*x(i)-top3z(i))*x(i)
                exch(i)       = do4z(i)*(op3z(i)*x(i)-tzp3(i))*x(i)
                pot2cm_energy = pot2cm_energy + coul(i)
            enddo
            rad           = sqrt((exch(1)-exch(2))**2+(exch(2)-exch(3))**2+(exch(3)-exch (1))**2)
            pot2cm_energy = pot2cm_energy - rad/sqrt(2.e0_wp)
            !   Compute the derivatives of the energy with respect to the internal
            !   coordinates.
            if( usricm_nder .eq. 1 )then
                s = exch(1) + exch(2) + exch(3)
                do i = 1, 3
                    pot2cm_dedr(i) = b(i)*x(i)*((3.0e0_wp*exch(i)-s)/sqrt(2.e0_wp)* (op3z(i)*x(i)-zp3(i))/rad- zp3(i)* &
                        & x(i)+ op3z(i))
                enddo
            endif
        end subroutine v3pot

        !> @brief Evaluate V4POT for the OH3_LEPS_SE_1986 potential.
        subroutine v4pot
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: t1, t2
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            t1 = exp(-concom_clam(1)*(pot2cm_r(1)-concom_alp(1))**2-concom_clam(1)*(pot2cm_r(2)-concom_alp(1))**2- &
                & concom_clam(3)*(pot2cm_r(3)-concom_alp(3)) **2-concom_clam(3)*(pot2cm_r(4)-concom_alp(3))**2)* &
                & concom_acon(1)
            t2 = exp(-concom_clam(2)*(pot2cm_r(1)-concom_alp(2))**2-concom_clam(2)*(pot2cm_r(2)-concom_alp(2))**2- &
                & concom_clam(4)*(pot2cm_r(3)-concom_alp(4)) **2-concom_clam(4)*(pot2cm_r(4)-concom_alp(4))**2)* &
                & concom_acon(2)
            pot2cm_energy = t1+t2
            if( usricm_nder .eq. 1 )then
                pot2cm_dedr(1) = -2.0e0_wp*(t1*concom_clam(1)*(pot2cm_r(1)-concom_alp(1))+t2*concom_clam(2)* &
                    & (pot2cm_r(1)- concom_alp(2)))
                pot2cm_dedr(2) = -2.0e0_wp*(t1*concom_clam(1)*(pot2cm_r(2)-concom_alp(1))+t2*concom_clam(2)* &
                    & (pot2cm_r(2)- concom_alp(2)))
                pot2cm_dedr(3) = -2.0e0_wp*(t1*concom_clam(3)*(pot2cm_r(3)-concom_alp(3))+t2*concom_clam(4)* &
                    & (pot2cm_r(3)- concom_alp(4)))
                pot2cm_dedr(4) = -2.0e0_wp*(t1*concom_clam(3)*(pot2cm_r(4)-concom_alp(3))+t2*concom_clam(4)* &
                    & (pot2cm_r(4)- concom_alp(4)))
            endif
        end subroutine v4pot

        !> @brief Evaluate VH2O for the OH3_LEPS_SE_1986 potential.
        subroutine vh2o
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: dp(3), dq(3), p, q(3), s(3), tanlg, temp, trm1, x(3), xmax1, xmax2, xx
            integer             :: i
            integer, parameter  :: n3atom = 75
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
            integer, parameter  :: natom  = 25
            !     XMAX1 FOR WHEN TANH SET=1.0 ON VAX
            !     XMAX2 FOR PREVENTING OVERFLOWS ON VAX
            !     STATEMENT FUNCTION
            tanlg(xx) = 2.0e0_wp/(1.0e0_wp+exp(2.0e0_wp*xx))
            s(1)      = pot2cm_r(1)-concom_reoh
            s(3)      = pot2cm_r(2)-concom_reoh
            s(2)      = pot2cm_r(3)-concom_rehh
            xmax1     = 15.0e0_wp; xmax2 = 43.0e0_wp
            do i = 1, 3
                x(i) = 0.5e0_wp*concom_gam(i)*s(i)
                q(i) = 1.0e0_wp-tanh(x(i))
                if( x(i) < xmax1 )then
                    dq(i) = -0.5e0_wp*concom_gam(i)/cosh(x(i))**2
                elseif( x(i) < xmax2 )then
                    q(i)  = tanlg(x(i))
                    dq(i) = -0.5e0_wp*concom_gam(i)/cosh(x(i))**2
                else
                    q(i)  = 0.0e0_wp
                    dq(i) = 0.0e0_wp
                endif
            enddo
            p = concom_con(1)+concom_con(2)*(s(1)+s(3))+concom_con(3)*s(2)+0.5e0_wp*concom_con(4)*(s(1)*s(1)+s(3)* &
                & s( 3))+ 0.5e0_wp*concom_con(5)*s(2)*s(2)+concom_con(6)*s(2)*(s(1)+s(3))+ concom_con(7)*s(1)*s(3)
            pot2cm_energy = q(1)*q(2)*q(3)*p
            if( usricm_nder .eq. 1 )then
                dp(1) = concom_con(2)+concom_con(4)*s(1)+concom_con(6)*s(2)+concom_con(7)*s(3)
                dp(2) = concom_con(3)+concom_con(5)*s(2)+concom_con(6)*(s(1)+s(3))
                dp(3) = concom_con(2)+concom_con(4)*s(3)+concom_con(6)*s(2)+concom_con(7)*s(1)
                do i = 1, 3
                    trm1 = 0.0e0_wp
                    if( q(i) /= 0.0e0_wp ) trm1 = dq(i)/q(i)
                    pot2cm_dedr(i) = pot2cm_energy*(trm1+dp(i)/p)
                enddo
                temp           = pot2cm_dedr(2)
                pot2cm_dedr(2) = pot2cm_dedr(3)
                pot2cm_dedr(3) = temp
            endif
        end subroutine vh2o

        !> @brief Evaluate POTINFO for the OH3_LEPS_SE_1986 potential.
        subroutine potinfo
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i, j, k, kmax
            integer, parameter  :: n3atom = 75
            integer, parameter  :: natom  = 25
            integer, parameter  :: isurf  = 5
            integer, parameter  :: jsurf  = isurf*(isurf+1)/2
            real(wp), parameter :: pi     = 3.141592653589793e0_wp
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
            if( kmax-1 .lt. infocm_msurf )then
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

        !> @brief Evaluate ANCVRT for the OH3_LEPS_SE_1986 potential.
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
            elseif( usricm_nflag(2) .eq. 3 )then
                units = 'EV                  '
            elseif( usricm_nflag(2) .eq. 4 )then
                units = 'KCAL PER MOLE       '
            elseif( usricm_nflag(2) .eq. 5 )then
                units = 'WAVENUMBERS         '
            elseif( usricm_nflag(2) .eq. 6 )then
                units = 'KILOJOULES PER MOLE '
            endif
            utilcm_cnvrtd  = 1.e0_wp
            utilcm_cnvrte  = 1.e0_wp
            utilcm_cnvrtde = 1.e0_wp
            if( usricm_nflag(1) .eq. 2 ) utilcm_cnvrtd = bohr_a
            if( usricm_nflag(2) .eq. 2 )then
                utilcm_cnvrte = utilcm_cnvrte*htomillh
            elseif( usricm_nflag(2) .eq. 3 )then
                utilcm_cnvrte = utilcm_cnvrte*htoev
            elseif( usricm_nflag(2) .eq. 4 )then
                utilcm_cnvrte = utilcm_cnvrte*htokcal
            elseif( usricm_nflag(2) .eq. 5 )then
                utilcm_cnvrte = utilcm_cnvrte*htowave
            elseif( usricm_nflag(2) .eq. 6 )then
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

        !> @brief Evaluate CARTOU for the OH3_LEPS_SE_1986 potential.
        subroutine cartou
            use chempotpy_legacy_coordinates, only: convert_cartesian_input
            use chempotpy_kinds,              only: wp
            implicit none
            integer            :: i, j
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            call convert_cartesian_input(usricm_cart, usricm_nulbl, infocm_natoms, utilcm_cnvrtd, utilcm_ireorder == 1, infocm_cartnu)
        end subroutine cartou

        !> @brief Evaluate CARTTOR for the OH3_LEPS_SE_1986 potential.
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
            call cartesian_to_potlib_coordinates(infocm_cartnu, infocm_natoms, infocm_icartr, pt1cm_r, coordinate_status)
            if( coordinate_status /= 0 )then
                write(usricm_nflag(18), '(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') infocm_icartr
                stop
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the OH3_LEPS_SE_1986 potential.
        subroutine eunitzero
            use chempotpy_kinds, only: wp
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

        !> @brief Evaluate RTOCART for the OH3_LEPS_SE_1986 potential.
        subroutine rtocart
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: natom  = 25
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            real(wp)           :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, &
                & termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes(n3atom, isurf), &
                & ygs(n3atom), yij(n3atom, jsurf), za, zb, zc, zd
            integer :: i, ind, j, j1, j2, k, m1
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
            elseif( infocm_icartr .eq. 2 )then
                do i = 1, infocm_natoms
                    utilcm_dgscartnu(i,1) = 0.e0_wp
                    utilcm_dgscartnu(i,2) = 0.e0_wp
                    utilcm_dgscartnu(i,3) = 0.e0_wp
                    if( utilcm_ksdiag .ne. 0 )then
                        do j1=utilcm_ksdiag, utilcm_kediag
                            utilcm_descartnu(i,1,j1) = 0.e0_wp
                            utilcm_descartnu(i,2,j1) = 0.e0_wp
                            utilcm_descartnu(i,3,j1) = 0.e0_wp
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do j2=utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(i,1,j2) = 0.e0_wp
                            utilcm_dijcartnu(i,2,j2) = 0.e0_wp
                            utilcm_dijcartnu(i,3,j2) = 0.e0_wp
                        enddo
                    endif
                    do j = 1, infocm_natoms
                        if( j .lt. i )then
                            m1 = infocm_natoms*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
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
                        elseif( utilcm_ksoffd .gt. 0 )then
                            do j2=utilcm_ksoffd, utilcm_keoffd
                                y                        = pt5cm_deijdr(m1,j2)
                                utilcm_dijcartnu(i,1,j2) = utilcm_dijcartnu(i,1,j2) + termx*y
                                utilcm_dijcartnu(i,2,j2) = utilcm_dijcartnu(i,2,j2) + termy*y
                                utilcm_dijcartnu(i,3,j2) = utilcm_dijcartnu(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            elseif( infocm_icartr .eq. 3 )then
                do i = 1, infocm_natoms
                    ygs(i) = pt1cm_degsdr(i)/pt1cm_r(i)
                    if( utilcm_ksdiag .ne. 0 )then
                        do j=utilcm_ksdiag, utilcm_kediag
                            yes(i, j) = pt4cm_deesdr(i,j)/pt1cm_r(i)
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do k=utilcm_ksoffd, utilcm_keoffd
                            yij(i, k) = pt5cm_deijdr(i,k)/pt1cm_r(i)
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
                            utilcm_descartnu(1,k,j1) = term12*yes(1, j1) + term13*yes(3, j1)
                            utilcm_descartnu(2,k,j1) = -term12*yes(1, j1) + term23*yes(2, j1)
                            utilcm_descartnu(3,k,j1) = -term13*yes(3, j1) - term23*yes(2, j1)
                        enddo
                    endif
                    if( utilcm_ksoffd .ne. 0 )then
                        do j2=utilcm_ksoffd, utilcm_keoffd
                            utilcm_dijcartnu(1,k,j2) = term12*yij(1, j2) + term13*yij(3, j2)
                            utilcm_dijcartnu(2,k,j2) = -term12*yij(1, j2) + term23*yij(2, j2)
                            utilcm_dijcartnu(3,k,j2) = -term13*yij(3, j2) - term23*yij(2, j2)
                        enddo
                    endif
                enddo
            elseif( infocm_icartr .eq. 4 )then
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
                write(usricm_nflag(18), '(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') infocm_icartr
                stop
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the OH3_LEPS_SE_1986 potential.
        subroutine dedcou
            use chempotpy_kinds, only: wp
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
end module chempotpy_oh3_leps_se_1986_pes

!> @brief Expose OH3_LEPS_SE_1986 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_oh3_leps_se_1986_pes, only: evaluate
    use chempotpy_kinds,                only: wp
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
