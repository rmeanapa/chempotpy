!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the C2H7_VBMM potential.
subroutine pes( x, igrad, p, g, d )

    implicit none

    ! number of electronic state
    integer, parameter            :: nstates = 1
    integer, parameter            :: natoms  = 9
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)

    double precision :: v
    integer          :: iatom, idir, j, istate, jstate
    !initialize
    p = 0.d0
    g = 0.d0
    d = 0.d0

    if (igrad==0) then
        call c2h7pot(x,v)
        v = v/23.0609
        do istate=1,nstates
            p(istate) = v
        enddo
    else
        write(*,*) 'Only energy is available'
    endif

end subroutine pes

!*******************************************************************
!> @brief Evaluate c2h7poth for the C2H7_VBMM potential.
subroutine c2h7poth( x, eng_kcal, grad, xnorm, hess )
    !*******************************************************************
    implicit none

    integer, parameter :: npmm = 23
    integer, parameter :: npjg = 34
    double precision   :: x(9,3)
    double precision   :: eng_kcal
    double precision   :: grad(9,3)
    double precision   :: xnorm(9)
    double precision   :: hess(27,27)

    integer                     :: i,j
    integer                     :: i1,j1,k1
    integer                     :: i2,j2,k2
    double precision            :: z
    double precision            :: x1(9,3)
    double precision            :: x2(9,3)
    double precision            :: x3(9,3)
    double precision            :: x4(9,3)
    double precision            :: v1
    double precision            :: v2
    double precision            :: v3
    double precision            :: v4
    double precision, parameter :: delta = 1.0d-4

    ! get energy
    call c2h7pot(x,eng_kcal)

    ! get numerical gradients
    do i=1,9
        do j=1,3
            x1(i,j) = x(i,j)
            x2(i,j) = x(i,j)
            x3(i,j) = x(i,j)
            x4(i,j) = x(i,j)
        end do
    end do

    do i=1,9
        xnorm(i) = 0.0d0
        do j=1,3
            x1(i,j) = x(i,j)-delta
            x2(i,j) = x(i,j)+delta
            call c2h7pot(x1,v1)
            call c2h7pot(x2,v2)
            grad(i,j) = (v2-v1)*0.50d0/delta
            x1(i,j)   = x(i,j)
            x2(i,j)   = x(i,j)
            xnorm(i)  = xnorm(i)+(grad(i,j)**2)
        end do
        xnorm(i) = dsqrt(xnorm(i))
    end do

    ! get numerical hessian
    k1 = 0
    do i2=1,9
        do j2=1,3
            k2 = k2+1
            k1 = 0
            do i1=1,9
                do j1=1,3
                    k1 = k1+1
                    if(k1>k2) then
                        x1(i1,j1) = x(i1,j1)+delta
                        x1(i2,j2) = x(i2,j2)+delta
                        x2(i1,j1) = x(i1,j1)-delta
                        x2(i2,j2) = x(i2,j2)+delta
                        x3(i1,j1) = x(i1,j1)+delta
                        x3(i2,j2) = x(i2,j2)-delta
                        x4(i1,j1) = x(i1,j1)-delta
                        x4(i2,j2) = x(i2,j2)-delta
                        call c2h7pot(x1,v1)
                        call c2h7pot(x2,v2)
                        call c2h7pot(x3,v3)
                        call c2h7pot(x4,v4)
                        z           = 0.250d0*(v1-v2-v3+v4)
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
                    elseif(k1==k2) then
                        x1(i1,j1) = x(i1,j1)+delta
                        x2(i2,j2) = x(i2,j2)-delta
                        call c2h7pot(x1,v1)
                        call c2h7pot(x2,v2)
                        z           = (v1-eng_kcal-eng_kcal+v2)
                        hess(k1,k1) = z/(delta*delta)
                        x1(i1,j1)   = x(i1,j1)
                        x2(i2,j2)   = x(i2,j2)
                    end if
                end do
            end do
        end do
    end do

end subroutine c2h7poth
! end of c2h7poth
!*******************************************************************
!> @brief Evaluate c2h7potg for the C2H7_VBMM potential.
subroutine c2h7potg( x, eng_kcal, grad, xnorm )
    !*******************************************************************
    implicit none

    integer, parameter :: npmm = 23
    integer, parameter :: npjg = 34
    double precision   :: x(9,3)
    double precision   :: eng_kcal
    double precision   :: grad(9,3)
    double precision   :: xnorm(9)

    integer                     :: i,j
    double precision            :: x1(9,3)
    double precision            :: x2(9,3)
    double precision            :: v1
    double precision            :: v2
    double precision, parameter :: delta = 1.0d-4

    ! get energy
    call c2h7pot(x,eng_kcal)

    ! get numerical gradients
    do i=1,9
        do j=1,3
            x1(i,j) = x(i,j)
            x2(i,j) = x(i,j)
        end do
    end do

    do i=1,9
        xnorm(i) = 0.0d0
        do j=1,3
            x1(i,j) = x(i,j)-delta
            x2(i,j) = x(i,j)+delta
            call c2h7pot(x1,v1)
            call c2h7pot(x2,v2)
            grad(i,j) = (v2-v1)*0.50d0/delta
            x1(i,j)   = x(i,j)
            x2(i,j)   = x(i,j)
            xnorm(i)  = xnorm(i)+(grad(i,j)**2)
        end do
        xnorm(i) = dsqrt(xnorm(i))
    end do

end subroutine c2h7potg
! end of routine c2h7potg
!**************************************************************
!> @brief Evaluate c2h7pot for the C2H7_VBMM potential.
subroutine c2h7pot( x, eng_kcal )
    !**************************************************************
    ! Coordinates in Ang
    ! Energy in kcal
    !
    implicit none

    integer, parameter :: npmm = 25
    integer, parameter :: npjg = 34
    double precision   :: x(9,3)
    double precision   :: c1(3)
    double precision   :: h2(3)
    double precision   :: h3(3)
    double precision   :: h4(3)
    double precision   :: c5(3)
    double precision   :: h6(3)
    double precision   :: h7(3)
    double precision   :: h8(3)
    double precision   :: h9(3)
    double precision   :: eng_kcal
    double precision   :: vmm3
    double precision   :: vjg
    double precision   :: etotalmm
    double precision   :: paramm(npmm)
    double precision   :: parajg(npjg)

    ! average bond length
    double precision :: avgch
    double precision :: avgbh
    double precision :: rxncrd
    double precision :: r1,r2,r3
    double precision :: b1,b2,b3

    ! scale=C-H/C-C
    !%%%  parameter(scale=0.7153679930d0)
    double precision, parameter :: scale = 0.7171830610d0
    double precision            :: h5prime(3)

    integer            :: i
    logical, parameter :: lmm3 = .false.
    !--------------------------------------------
    ! Convert geometry from:
    !
    !     H2 H6
    !     |  |
    !  H4-C1-C5-H8.....H9
    !     |  |
    !     H3 H7
    !
    ! to:
    !
    !     H6 H3
    !     |  |
    !  H8-C5-C1-H2.....H9
    !     |  |
    !     H7 H4
    !
    do i=1,3
        c1(i)      = x(5,i)
        h2(i)      = x(8,i)
        h3(i)      = x(6,i)
        h4(i)      = x(7,i)
        c5(i)      = x(1,i)
        h6(i)      = x(2,i)
        h7(i)      = x(3,i)
        h8(i)      = x(4,i)
        h9(i)      = x(9,i)
        h5prime(i) = (scale*c5(i)) + ((1.0d0-scale)*c1(i))
    end do

    ! Calculate avg C-H bond length for 3 QM H
    r1 = 0.0d0
    r2 = 0.0d0
    r3 = 0.0d0
    b1 = 0.0d0
    b2 = 0.0d0
    b3 = 0.0d0
    do i=1,3
        r1 = r1+( (c1(i)-h2(i))**2 )
        r2 = r2+( (c1(i)-h3(i))**2 )
        r3 = r3+( (c1(i)-h4(i))**2 )
        b1 = b1+( (h9(i)-h2(i))**2 )
        b2 = b2+( (h9(i)-h3(i))**2 )
        b3 = b3+( (h9(i)-h4(i))**2 )
    end do
    r1     = dsqrt(r1)
    r2     = dsqrt(r2)
    r3     = dsqrt(r3)
    b1     = dsqrt(b1)
    b2     = dsqrt(b2)
    b3     = dsqrt(b3)
    avgch  = (r1+r2+r3)/3.0d0
    avgbh  = (1.0d0/b1)+(1.0d0/b2)+(1.0d0/b3)
    rxncrd = avgch*avgbh
    !
    ! Get energy
    !
    if(lmm3) then
        call pmm3(npmm,paramm)
        call ethane(npmm,paramm,rxncrd, &
            & c1,h2,h3,h4,c5,h6,h7,h8,vmm3,etotalmm)
        eng_kcal = etotalmm
    else
        call pmm3_mod(npmm,paramm,avgch,rxncrd)
        call pjg_mod(npjg,parajg,avgch)
        call ethane(npmm,paramm,rxncrd, &
            & c1,h2,h3,h4,c5,h6,h7,h8,vmm3,etotalmm)
        call xch5pot(npjg,parajg,rxncrd, &
            & c1,h2,h3,h4,h5prime,h9,vjg)
        eng_kcal = vjg+vmm3
    end if

end subroutine c2h7pot
! end of routine c2h7pot

!**************************************************
!> @brief Evaluate pmm3_mod for the C2H7_VBMM potential.
subroutine pmm3_mod( np, p, ravg, rxncrd )
    !**************************************************
    implicit none

    ! Output Variables
    integer          :: np
    double precision :: p(np)
    double precision :: ravg
    double precision :: rxncrd
    ! Local Variables
    double precision :: r0ch,r0cc
    double precision :: t0hch,t0cch
    double precision :: t0hcrh,t0ccrh
    double precision :: r0crh,r0ccr
    !     parameter(r0ch=1.1120d0,r0cc=1.52470d0)
    !##   parameter(r0ch=1.09120d0,r0cc=1.51950d0)

    double precision, parameter :: ksch = 4.740d0
    double precision, parameter :: kscc = 4.490d0
    !%%%  parameter(t0hch=107.80d0,t0cch=110.70d0)
    !%%%  parameter(t0hcRh=107.60d0,t0ccRh=109.310d0)
    double precision, parameter :: kthch   = 0.550d0
    double precision, parameter :: ktcch   = 0.590d0
    double precision, parameter :: ksbhch  = 0.0d0
    double precision, parameter :: ksbcch  = 0.080d0
    double precision, parameter :: kbbhch  = 0.0d0
    double precision, parameter :: kbbcch  = 0.30d0**2
    double precision, parameter :: v1hcch  = 0.0d0
    double precision, parameter :: v2hcch  = 0.0d0
    double precision, parameter :: v3hcch  = 0.2380d0
    double precision, parameter :: ksthcch = 0.0590d0
    double precision, parameter :: espc    = 0.0270d0
    double precision, parameter :: rvc     = 2.040d0
    double precision, parameter :: esph    = 0.020d0
    double precision, parameter :: rvh     = 1.620d0
    double precision, parameter :: scalech = 0.9230d0
    !%%%  parameter(r0cRh=1.1120d0,r0ccR=1.52470d0)
    !----------------------------------------
    ! Switching
    double precision :: ravg0
    double precision :: cr,cp,cs
    !--------------------------------------
    ! Modified parameters
    !--------------------------------------
    integer          :: i,j
    integer          :: ipt
    double precision :: x
    double precision :: t1,t2,t3
    double precision :: t4,t5,t6
    double precision :: z1,z2,z3
    double precision :: z4,z5,z6
    double precision :: ar0ch,ar0cc
    double precision :: asch,ascc
    double precision :: av3hcch
    double precision :: aksthcch
    double precision :: akthch
    double precision :: aktcch
    double precision :: aesph
    double precision :: arvh
    double precision :: at0hch
    double precision :: ar0crh
    double precision :: ar0ccr
    double precision :: at0ccrh
    !%%%  parameter(r0ch=1.1120d0,r0cc=1.52470d0)
    !     parameter(r0ch=1.092510d0,r0cc=1.52720d0)

    !%%% for ethane

    ipt = 4

    if(ipt .eq. 1) then
        !       ETHANE
        r0ch   = 1.0850d0
        r0cc   = 1.510d0
        t0hch  = 107.52790d0
        t0cch  = 111.34990d0
        r0crh  = r0ch
        r0ccr  = 1.47110d0
        t0hcrh = 120.0d0
        t0ccrh = 120.0d0
    elseif(ipt .eq. 2) then
        !       TS
        r0ch   = 1.0850d0
        r0cc   = 1.510d0
        t0hch  = 107.52790d0
        t0cch  = 111.34990d0
        r0crh  = r0ch
        r0ccr  = 1.47110d0
        t0hcrh = 120.0d0
        t0ccrh = 120.0d0
    elseif(ipt .eq. 3) then
        !       ETHYL
        r0ch   = 1.0850d0
        r0cc   = 1.510d0
        t0hch  = 107.52790d0
        t0cch  = 111.34990d0
        r0crh  = r0ch
        r0ccr  = 1.47110d0
        t0hcrh = 120.0d0
        t0ccrh = 120.0d0
    elseif(ipt .eq. 4) then
        !--------------------------------------
        ! geometry
        r0ch   = 1.0850d0
        r0cc   = 1.50d0
        t0hch  = 107.52790d0
        t0cch  = 111.34990d0
        r0crh  = r0ch
        r0ccr  = 1.47110d0
        t0hcrh = 120.0d0
        t0ccrh = 120.0d0
        !%%%
        t0ccrh = 117.2989954203791430d0
        t0ccrh = 121.0d0
        t0cch  = t0ccrh-7.10d0

        t1     = 6.0d0
        t2     = 0.90d0
        t3     = 0.50d0*(1.0d0+dtanh(t1*(rxncrd-t2)))
        t4     = 1.0d0-t3
        t0cch  = (t0cch*t3)+(111.34990d0*t4)
        t0ccrh = (t0ccrh*t3)+(111.34990d0*t4)

        !       t1=6.0d0
        !       t2=20.0d0
        !       t3=0.50d0*(1.0d0+dtanh(t1*(rxncrd-t2)))
        !       t4=1.0d0-t3
        !       t0cch=(t0cch*t4)+(120.0d0*t3)
        !       t0ccRh=(t0ccRh*t4)+(120.0d0*t3)
        !--------------------------------------
        ! freq and energy
        cs    = 10.0d0
        ravg0 = 1.17320d0
        cp    = 0.50d0*(1.0d0+(dtanh(cs*(ravg-ravg0))))
        cr    = 1.0d0-cp

        ascc   = 0.70d0
        aktcch = 0.90d0
        akthch = 0.90d0
        !
        !       z1=10.0d0
        !       z2=2.144106410736830d0
        !       z3=dexp(-z1*((rxncrd-z2)**2))
        !       Ascc=1.0d0-(0.30d0*z3)
        !       Aktcch=1.0d0-(0.10d0*z3)
        !       Akthch=1.0d0-(0.10d0*z3)

        t1 = 6.0d0
        t2 = 0.950d0
        t3 = 0.50d0*(1.0d0+dtanh(t1*(rxncrd-t2)))
        t4 = (1.0d0*(1.0d0-t3))+(0.70d0*t3)
        t5 = (1.0d0*(1.0d0-t3))+(0.90d0*t3)
        t6 = (1.0d0*(1.0d0-t3))+(0.90d0*t3)

        z1 = 6.0d0
        z2 = 10.0d0
        z3 = 0.50d0*(1.0d0+dtanh(z1*(rxncrd-z2)))

        z4 = (t4*(1.0d0-z3))+(1.0d0*z3)
        z5 = (t5*(1.0d0-z3))+(1.0d0*z3)
        z6 = (t6*(1.0d0-z3))+(1.0d0*z3)

        ascc   = z4
        aktcch = z5
        akthch = z6

    end if
    !###
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
! end of routine pmm3_mod
!**************************************************
!> @brief Evaluate pjg_mod for the C2H7_VBMM potential.
subroutine pjg_mod( np, para, ravg )
    !**************************************************
    implicit none

    ! I/O variables
    integer          :: np
    double precision :: para(np)
    double precision :: ravg
    ! Parameter list
    double precision :: r0cb
    double precision :: r0hh
    double precision :: alpch
    double precision :: r0ch
    !
    double precision :: tau
    !
    !
    double precision, parameter :: fact1 = 0.041840d0
    double precision, parameter :: fact2 = 6.022045d0
    ! streching parameters
    double precision, parameter :: a1ch  = 1.7130d0
    double precision, parameter :: b1ch  = 0.1350d0
    double precision, parameter :: c1ch  = 6.6140d0
    double precision, parameter :: alpcb = 1.8530285d0
    double precision, parameter :: d1cb  = 26.4090d0*fact1
    double precision, parameter :: d3cb  = 20.0630d0*fact1
    !     parameter (r0cb=1.093970d0)
    !%%%  parameter (r0cb=1.092510d0)
    !%%%  parameter (r0ch=1.092510d0)

    !##   parameter (r0cb=1.09030d0)
    !##   parameter (r0ch=1.09070d0)

    double precision, parameter :: alphh = 1.94570d0
    double precision, parameter :: d1hh  = 109.4580d0*fact1
    double precision, parameter :: d3hh  = 39.6640d0*fact1
    !%%%  parameter (r0hh=0.741910d0)
    !     parameter (alphch=CALC_AT_RUNTIME)
    double precision, parameter :: d1ch = 112.230d0*fact1
    double precision, parameter :: d3ch = 38.8340d0*fact1
    !     parameter (r0ch=1.093970d0)
    !##   parameter (r0ch=1.092510d0)

    ! out-of-plane bending parameters
    !     parameter (tau=CALC_AT_RUNTIME)
    double precision, parameter :: aphi   = 0.52879030d0
    double precision, parameter :: bphi   = 0.40066380d0
    double precision, parameter :: cphi   = 1.92099370d0
    double precision, parameter :: atheta = 0.90787140d0
    double precision, parameter :: btheta = 0.35488590d0
    double precision, parameter :: ctheta = 1.89154970d0
    double precision, parameter :: a3s    = 0.14191470d0
    double precision, parameter :: b3s    = -0.3068450d0
    double precision, parameter :: fd0    = 0.0957500d0*fact2
    double precision, parameter :: hd0    = 0.1915000d0*fact2
    ! in-plane bending parameters
    double precision, parameter :: a1s  = 1.53136810d-7
    double precision, parameter :: a2s  = 1.01474020d-7
    double precision, parameter :: b1s  = -4.6696246d0
    double precision, parameter :: b2s  = -12.363798d0
    double precision, parameter :: ak   = 0.1260d0*fact2
    double precision, parameter :: kch3 = 0.40770d0*fact2
    double precision, parameter :: aa1  = 3.2139520d0
    double precision, parameter :: aa2  = 1.5999630d0
    double precision, parameter :: aa3  = 2.1659530d0
    double precision, parameter :: aa4  = 11.5699770d0
    !----------------------------------------
    ! Switching
    double precision :: ravg0
    double precision :: cr,cp,cs
    !----------------------------------------
    !###
    integer          :: i
    integer          :: ipt
    logical          :: labinitio
    double precision :: ad1hh
    double precision :: ad1ch
    double precision :: aalpcb
    !
    !%%% QM parameters

    ipt = 4

    ad1hh  = 1.0d0
    ad1ch  = 1.0d0
    aalpcb = 1.0d0

    if(ipt .eq. 1) then
        !       ETHANE
        r0ch = 1.0850d0
        r0cb = 1.400d0
        r0hh = 0.680d0
    elseif(ipt .eq. 2) then
        !       TS
        r0ch = 1.080d0
        r0cb = 1.400d0
        r0hh = 0.660d0
    elseif(ipt .eq. 3) then
        !       ETHYL
        r0ch = 1.0750d0
        r0cb = 1.400d0
        r0hh = 0.73590d0
    elseif(ipt .eq. 4) then
        !------------------------------------------------------
        ! for geometry
        cs    = 10.0d0
        ravg0 = 1.17320d0
        cp    = 0.50d0*(1.0d0+(dtanh(cs*(ravg-ravg0))))
        cr    = 1.0d0-cp
        r0ch  = (1.0850d0*cr)+(1.0750d0*cp)
        r0cb  = 1.400d0
        !       r0hh=0.73590d0
        !------------------------------------------------------
        ! for rxn energy and freq
        cs    = 10.0d0
        ravg0 = 1.7320d0
        cp    = 0.50d0*(1.0d0+(dtanh(cs*(ravg-ravg0))))
        cr    = 1.0d0-cp
        ad1hh = (1.0d0*cr)+(1.027420d0*cp)
        ad1hh = (1.0d0*cr)+(1.027500d0*cp)
        r0hh  = (0.720d0*cr)+(0.73590d0*cp)
        r0hh  = (0.70d0*cr)+(0.73590d0*cp)
    end if
    !======================================================
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
! end of routine pjp_mod
!**************************************************
!> @brief Evaluate pmm3 for the C2H7_VBMM potential.
subroutine pmm3( np, p )
    !**************************************************
    implicit none

    ! Output Variables
    integer          :: np
    double precision :: p(np)
    ! Local Variables

    double precision, parameter :: r0ch    = 1.1120d0
    double precision, parameter :: r0cc    = 1.52470d0
    double precision, parameter :: ksch    = 4.740d0
    double precision, parameter :: kscc    = 4.490d0
    double precision, parameter :: t0hch   = 107.80d0
    double precision, parameter :: t0cch   = 110.70d0
    double precision, parameter :: t0hcrh  = 107.60d0
    double precision, parameter :: t0ccrh  = 109.310d0
    double precision, parameter :: kthch   = 0.550d0
    double precision, parameter :: ktcch   = 0.590d0
    double precision, parameter :: ksbhch  = 0.0d0
    double precision, parameter :: ksbcch  = 0.080d0
    double precision, parameter :: kbbhch  = 0.0d0
    double precision, parameter :: kbbcch  = 0.30d0**2
    double precision, parameter :: v1hcch  = 0.0d0
    double precision, parameter :: v2hcch  = 0.0d0
    double precision, parameter :: v3hcch  = 0.2380d0
    double precision, parameter :: ksthcch = 0.0590d0
    double precision, parameter :: espc    = 0.0270d0
    double precision, parameter :: rvc     = 2.040d0
    double precision, parameter :: esph    = 0.020d0
    double precision, parameter :: rvh     = 1.620d0
    double precision, parameter :: scalech = 0.9230d0
    double precision, parameter :: r0crh   = 1.1120d0
    double precision, parameter :: r0ccr   = 1.52470d0
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
! end of routine pmm3
!**************************************************
!> @brief Evaluate pjg for the C2H7_VBMM potential.
subroutine pjg( np, para )
    !**************************************************
    implicit none

    ! I/O variables
    integer          :: np
    double precision :: para(np)
    ! Parameter list
    double precision :: alpch
    !
    double precision :: tau
    !
    !
    double precision, parameter :: fact1 = 0.041840d0
    double precision, parameter :: fact2 = 6.022045d0
    ! streching parameters
    double precision, parameter :: a1ch  = 1.7130d0
    double precision, parameter :: b1ch  = 0.1350d0
    double precision, parameter :: c1ch  = 6.6140d0
    double precision, parameter :: alpcb = 1.8530285d0
    double precision, parameter :: d1cb  = 26.4090d0*fact1
    double precision, parameter :: d3cb  = 20.0630d0*fact1
    double precision, parameter :: r0cb  = 1.093970d0
    double precision, parameter :: alphh = 1.94570d0
    double precision, parameter :: d1hh  = 109.4580d0*fact1
    double precision, parameter :: d3hh  = 39.6640d0*fact1
    double precision, parameter :: r0hh  = 0.741910d0
    !     parameter (alphch=CALC_AT_RUNTIME)
    double precision, parameter :: d1ch = 112.230d0*fact1
    double precision, parameter :: d3ch = 38.8340d0*fact1
    double precision, parameter :: r0ch = 1.093970d0
    ! out-of-plane bending parameters
    !     parameter (tau=CALC_AT_RUNTIME)
    double precision, parameter :: aphi   = 0.52879030d0
    double precision, parameter :: bphi   = 0.40066380d0
    double precision, parameter :: cphi   = 1.92099370d0
    double precision, parameter :: atheta = 0.90787140d0
    double precision, parameter :: btheta = 0.35488590d0
    double precision, parameter :: ctheta = 1.89154970d0
    double precision, parameter :: a3s    = 0.14191470d0
    double precision, parameter :: b3s    = -0.3068450d0
    double precision, parameter :: fd0    = 0.0957500d0*fact2
    double precision, parameter :: hd0    = 0.1915000d0*fact2
    ! in-plane bending parameters
    double precision, parameter :: a1s  = 1.53136810d-7
    double precision, parameter :: a2s  = 1.01474020d-7
    double precision, parameter :: b1s  = -4.6696246d0
    double precision, parameter :: b2s  = -12.363798d0
    double precision, parameter :: ak   = 0.1260d0*fact2
    double precision, parameter :: kch3 = 0.40770d0*fact2
    double precision, parameter :: aa1  = 3.2139520d0
    double precision, parameter :: aa2  = 1.5999630d0
    double precision, parameter :: aa3  = 2.1659530d0
    double precision, parameter :: aa4  = 11.5699770d0

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
! end of routine pjp

!**************************************************
!> @brief Evaluate ethane for the C2H7_VBMM potential.
subroutine ethane( np, para, rxncrd, &
        & c1, h2, h3, h4, c5, h6, h7, h8, vmm3, etotalmm )
    !**************************************************
    implicit none

    ! Input variables............
    integer          :: np
    double precision :: para(np)
    double precision :: rxncrd
    double precision :: c1(3)
    double precision :: h2(3)
    double precision :: h3(3)
    double precision :: h4(3)
    double precision :: c5(3)
    double precision :: h6(3)
    double precision :: h7(3)
    double precision :: h8(3)
    ! Output variables............
    double precision :: vmm3
    double precision :: estr
    double precision :: ebend
    double precision :: estrbnd
    double precision :: ebndbnd
    double precision :: etor
    double precision :: estrtor
    double precision :: evdwaal
    double precision :: etotalmm

    ! parameters:
    ! Angles in Deg.
    ! Dist in Ang.
    double precision :: r0ch,r0cc
    double precision :: ksch,kscc
    double precision :: t0hch,t0cch
    double precision :: t0hcrh,t0ccrh
    double precision :: kthch,ktcch
    double precision :: ksbhch,ksbcch
    double precision :: kbbhch,kbbcch
    double precision :: v1hcch,v2hcch,v3hcch
    double precision :: ksthcch
    double precision :: esph,espc,rvh,rvc,scalech

    ! Func called:
    double precision :: xdist
    double precision :: xangle
    double precision :: xtorsion

    ! Local parameters:
    integer, parameter          :: iout     = 6
    logical, parameter          :: lreact   = .false.
    logical, parameter          :: lradical = .false.
    logical, parameter          :: lswitch  = .true.
    logical, parameter          :: ljg      = .true.
    logical, parameter          :: ldebugmm = .false.
    double precision, parameter :: pi       = 3.1415926540d0

    ! Stretch varaibles
    double precision :: rch(6),rcc
    double precision :: esch,escc

    ! Bend varaibles
    double precision :: thch(6),tcch(6)
    double precision :: ebhch,ebcch

    ! Stretch-Bend varaibles
    double precision :: esbcch

    ! Bend-Bend varaibles
    double precision :: ebbhch,ebbcch

    ! Torsional variables
    double precision :: phcch(9)
    double precision :: ethcch

    ! Stretch-xtorsion variables
    double precision :: esthcch

    ! van der Waal varaibles
    double precision :: sh2(3)
    double precision :: sh3(3)
    double precision :: sh4(3)
    double precision :: sh5(3)
    double precision :: sh6(3)
    double precision :: sh7(3)
    double precision :: sh8(3)
    double precision :: rvdwhh(9)
    double precision :: evdwhh

    ! energy compponents
    double precision :: astr(7)
    double precision :: abhch(6)
    double precision :: abcch(6)
    double precision :: asbcch(6)
    double precision :: abbhch(6)
    double precision :: abbcch(6)
    double precision :: athcch(9)
    double precision :: asthcch(9)
    double precision :: avdwhh(9)
    double precision :: tstr(7)
    double precision :: tbhch(6)
    double precision :: tbcch(6)
    double precision :: tsbcch(6)
    double precision :: tbbhch(6)
    double precision :: tbbcch(6)
    double precision :: tthcch(9)
    double precision :: tsthcch(9)
    double precision :: tvdwhh(9)

    ! Switching func
    double precision            :: scale
    double precision            :: sw1(3)
    double precision            :: sw2
    double precision            :: swt1(3)
    double precision            :: swt2
    double precision            :: r0crh
    double precision            :: r0ccr
    double precision            :: r0chsw
    double precision            :: r0ccsw
    double precision            :: t0hchsw
    double precision            :: t0cchsw
    double precision, parameter :: t1 = 22.0d0
    double precision, parameter :: t2 = 1.354030d0

    ! Misc varaibles
    integer          :: i,j,k
    double precision :: xrv,xesp
    double precision :: xsumrx
    double precision :: xsumjg
    double precision :: xstr
    double precision :: xbend
    double precision :: xstrbnd
    double precision :: xbndbnd
    double precision :: xtor
    double precision :: xstrtor
    double precision :: xvdw
    double precision :: z1
    double precision :: z2
    double precision :: z3
    double precision :: z4
    double precision :: y1
    double precision :: y2
    double precision :: y3
    !-------------------------------------------
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
    !-------------------------------------------
    !     write(IOUT,*)
    !     write(IOUT,*)'LREACT   :',LREACT
    !     write(IOUT,*)'LRADICAL :',LRADICAL
    !-------------------------------------------
    ! Evaluate bond xdist
    rch(1) = xdist(c1,h2)
    rch(2) = xdist(c1,h3)
    rch(3) = xdist(c1,h4)
    rch(4) = xdist(c5,h6)
    rch(5) = xdist(c5,h7)
    rch(6) = xdist(c5,h8)
    rcc    = xdist(c1,c5)
    !-------------------------------------------
    !###
    ! calculate the swithing func
    swt2 = 1.0d0
    do i=1,3
        swt1(i) = (1.0d0-dtanh(t1*(rch(i)-t2)))*0.50d0
        swt2    = swt2*swt1(i)
    end do
    ! switch C-C xdist
    r0cc = (swt2*r0cc)+((1.0d0-swt2)*r0ccr)
    !%%%
    ! switch force constants
    z1    = 10.0d0
    z2    = 2.144106410736830d0
    z3    = dexp(-z1*((rxncrd-z2)**2))
    scale = 1.0d0-(0.050d0*z3)
    !     scale=0.950d0

    ksch  = ksch*scale
    kthch = kthch*scale
    ktcch = ktcch*scale
    !-------------------------------------------
    ! calculate str energy
    call es(r0cc,kscc,rcc,escc)
    estr = escc
    do i=1,6
        call es(r0ch,ksch,rch(i),esch)
        estr    = estr + esch
        astr(i) = esch
    end do
    astr(7) = escc
    !-------------------------------------------
    ! calculate bend energy
    ! there are 12 bending terms
    thch(1) = xangle(h2,c1,h3)*180.0d0/pi
    thch(2) = xangle(h2,c1,h4)*180.0d0/pi
    thch(3) = xangle(h3,c1,h4)*180.0d0/pi
    thch(4) = xangle(h6,c5,h7)*180.0d0/pi
    thch(5) = xangle(h6,c5,h8)*180.0d0/pi
    thch(6) = xangle(h7,c5,h8)*180.0d0/pi

    tcch(1) = xangle(h2,c1,c5)*180.0d0/pi
    tcch(2) = xangle(h3,c1,c5)*180.0d0/pi
    tcch(3) = xangle(h4,c1,c5)*180.0d0/pi
    tcch(4) = xangle(h6,c5,c1)*180.0d0/pi
    tcch(5) = xangle(h7,c5,c1)*180.0d0/pi
    tcch(6) = xangle(h8,c5,c1)*180.0d0/pi

    ebend = 0.0d0
    do i=1,6
        call eb(t0hch,kthch,thch(i),ebhch)
        call eb(t0cch,ktcch,tcch(i),ebcch)
        ebend    = ebend + ebhch + ebcch
        abhch(i) = ebhch
        abcch(i) = ebcch
    end do
    !-------------------------------------------
    ! Calculate Stretch-bend terms
    ! there are 6 stretch-bend cross terms
    ! H2-C1-C5
    ! H3-C1-C5
    ! H4-C1-C5
    ! H6-C5-C1
    ! H7-C5-C1
    ! H8-C5-C1
    !
    estrbnd = 0.0d0
    do i=1,6
        call esb(r0ch,r0cc,t0cch,ksbcch,rch(i),rcc,tcch(i),esbcch)
        estrbnd   = estrbnd + esbcch
        asbcch(i) = esbcch
    end do
    !-------------------------------------------
    ! Bend-Bend terms
    ! there are 12 bend-bend terms
    ! H2-C1-H3 & H2-C1-H4
    ! H2-C1-H3 & H3-C1-H4
    ! H2-C1-H4 & H3-C1-H4
    ! H6-C5-H7 & H6-C5-H8
    ! H6-C5-H7 & H7-C5-H8
    ! H6-C5-H8 & H7-C5-H8
    !
    ! H2-C1-C5 & H3-C1-C5
    ! H2-C1-C5 & H4-C1-C5
    ! H3-C1-C5 & H4-C1-C5
    ! H6-C5-C1 & H6-C5-C1
    ! H6-C5-C1 & H7-C5-C1
    ! H7-C5-C1 & H8-C5-C1

    ebndbnd = 0.0d0
    k       = 0
    do i=1,3
        do j=i+1,3
            k = k + 1
            call ebb(t0hch,t0hch,kbbhch,thch(i),thch(j),ebbhch)
            call ebb(t0cch,t0cch,kbbcch,tcch(i),tcch(j),ebbcch)
            ebndbnd   = ebndbnd + ebbhch + ebbcch
            abbhch(k) = ebbhch
            abbcch(k) = ebbcch
        end do
    end do
    do i=4,6
        do j=i+1,6
            k = k + 1
            call ebb(t0hch,t0hch,kbbhch,thch(i),thch(j),ebbhch)
            call ebb(t0cch,t0cch,kbbcch,tcch(i),tcch(j),ebbcch)
            ebndbnd   = ebndbnd + ebbhch + ebbcch
            abbhch(k) = ebbhch
            abbcch(k) = ebbcch
        end do
    end do
    !-------------------------------------------
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
    !
    ! Scale the C-H bonds and evaluate the
    ! coordinates of the scaled H
    phcch(1) = xtorsion(h2,c1,c5,h6)
    phcch(2) = xtorsion(h2,c1,c5,h7)
    phcch(3) = xtorsion(h2,c1,c5,h8)
    phcch(4) = xtorsion(h3,c1,c5,h6)
    phcch(5) = xtorsion(h3,c1,c5,h7)
    phcch(6) = xtorsion(h3,c1,c5,h8)
    phcch(7) = xtorsion(h4,c1,c5,h6)
    phcch(8) = xtorsion(h4,c1,c5,h7)
    phcch(9) = xtorsion(h4,c1,c5,h8)

    etor = 0.0d0
    do i=1,9
        call et(v1hcch,v2hcch,v3hcch,phcch(i),ethcch)
        etor      = etor + ethcch
        athcch(i) = ethcch
        if(ldebugmm) then
            write(*,*) 'Torsion Angle:',phcch(i)
        end if
    end do
    !-------------------------------------------
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
    !
    estrtor = 0.0d0
    do i=1,9
        call est(r0cc,ksthcch,rcc,phcch(i),esthcch)
        estrtor    = estrtor + esthcch
        asthcch(i) = esthcch
    end do
    !-------------------------------------------
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
    !
    do i=1,3
        sh2(i) = (scalech*h2(i)) + ((1.0d0-scalech)*c1(i))
        sh3(i) = (scalech*h3(i)) + ((1.0d0-scalech)*c1(i))
        sh4(i) = (scalech*h4(i)) + ((1.0d0-scalech)*c1(i))
        sh6(i) = (scalech*h6(i)) + ((1.0d0-scalech)*c5(i))
        sh7(i) = (scalech*h7(i)) + ((1.0d0-scalech)*c5(i))
        sh8(i) = (scalech*h8(i)) + ((1.0d0-scalech)*c5(i))
    end do

    rvdwhh(1) = xdist(sh2,sh6)
    rvdwhh(2) = xdist(sh2,sh7)
    rvdwhh(3) = xdist(sh2,sh8)
    rvdwhh(4) = xdist(sh3,sh6)
    rvdwhh(5) = xdist(sh3,sh7)
    rvdwhh(6) = xdist(sh3,sh8)
    rvdwhh(7) = xdist(sh4,sh6)
    rvdwhh(8) = xdist(sh4,sh7)
    rvdwhh(9) = xdist(sh4,sh8)

    evdwaal = 0.0d0
    do i=1,9
        xrv  = rvh + rvh
        xesp = dsqrt(esph*esph)
        call evdw(xrv,xesp,rvdwhh(i),evdwhh)
        evdwaal   = evdwaal + evdwhh
        avdwhh(i) = evdwhh
    end do

    etotalmm = estr+ebend+estrbnd+ebndbnd+ etor+estrtor+evdwaal
    !-------------------------------------------------------
    ! Switch between the parameters of ethane and ethyl radical
    ! swt2=1 for ethane and 0 for ethyl
    !###
    r0ccsw  = r0cc
    r0chsw  = (swt2*r0ch)+((1.0d0-swt2)*r0crh)
    t0hchsw = (swt2*t0hch)+((1.0d0-swt2)*t0hcrh)
    t0cchsw = (swt2*t0cch)+((1.0d0-swt2)*t0ccrh)
    !-----------------------------------------------------------
    !%%%
    z2 = 2.141145010007599490d0
    z2 = 1.50d0
    z1 = 0.50d0*(1.0d0+dtanh(6.0d0*(rxncrd-z2)))
    z4 = ((1.0d0-z1)*t0cch)+(z1*t0ccrh)
    z1 = 117.2989954203791430d0

    sw2 = 1.0d0
    do i=1,3
        z1     = 10.0d0
        sw1(i) = (1.0d0-dtanh(z1*(rch(i)-t2)))*0.50d0
        sw2    = sw2*sw1(i)
    end do
    t0cchsw = (sw2*t0cch)+((1.0d0-sw2)*t0ccrh)
    !-----------------------------------------------------------

    if(lswitch) then
        do i=1,3
            call eb(t0hchsw,kthch,thch(i),ebhch)
            call eb(t0cchsw,ktcch,tcch(i),ebcch)
            call esb(r0chsw,r0ccsw,t0cchsw,ksbcch, &
                & rch(i),rcc,tcch(i),esbcch)
            abhch(i)  = ebhch
            abcch(i)  = ebcch
            asbcch(i) = esbcch
        end do
        k = 0
        do i=1,3
            do j=i+1,3
                k = k + 1
                call ebb(t0hchsw,t0hchsw,kbbhch,thch(i),thch(j),ebbhch)
                call ebb(t0cchsw,t0cchsw,kbbcch,tcch(i),tcch(j),ebbcch)
                abbhch(k) = ebbhch
                abbcch(k) = ebbcch
            end do
        end do
        ! compute engergies
        do i=1,7
            tstr(i) = astr(i)
        end do
        do i=1,6
            tbhch(i)  = abhch(i)
            tbcch(i)  = abcch(i)
            tsbcch(i) = asbcch(i)
            tbbhch(i) = abbhch(i)
            tbbcch(i) = abbcch(i)
        end do
        do i=1,9
            tthcch(i)  = athcch(i)
            tsthcch(i) = asthcch(i)
            tvdwhh(i)  = avdwhh(i)
        end do

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
        !
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

        xsumrx = 0.0d0
        xstr   = 0.0d0
        do i=1,7
            xstr = xstr+tstr(i)
        end do
        xbend   = 0.0d0
        xstrbnd = 0.0d0
        xbndbnd = 0.0d0
        do i=1,6
            xbend   = xbend+tbhch(i)+tbcch(i)
            xstrbnd = xstrbnd+tsbcch(i)
            xbndbnd = xbndbnd+tbbhch(i)+tbbcch(i)
        end do
        xtor    = 0.0d0
        xstrtor = 0.0d0
        xvdw    = 0.0d0
        do i=1,9
            xtor    = xtor+tthcch(i)
            xstrtor = xstrtor+tsthcch(i)
            xvdw    = xvdw+tvdwhh(i)
        end do
        xsumrx = xstr+xbend+xstrbnd+xbndbnd+xtor+xstrtor+xvdw

        !       call xdisp(xstr,xbend,xstrbnd,xbndbnd,
        !    1             xtor,xstrtor,xvdw,xsumRX)

        etotalmm = xsumrx
    end if
    !-------------------------------------------------------
    if(lradical) then
        !       modify the h-cR-h and h-cR-C bends and c-cR-h str-bnd
        !       for C-sp3 bonded to only two atoms
        do i=1,3
            call eb(t0hcrh,kthch,thch(i),ebhch)
            call eb(t0ccrh,ktcch,tcch(i),ebcch)
            call esb(r0ch,r0cc,t0ccrh,ksbcch,rch(i),rcc,tcch(i),esbcch)
            abhch(i)  = ebhch
            abcch(i)  = ebcch
            asbcch(i) = esbcch
        end do
        k = 0
        do i=1,3
            do j=i+1,3
                k = k + 1
                call ebb(t0hcrh,t0hcrh,kbbhch,thch(i),thch(j),ebbhch)
                call ebb(t0ccrh,t0ccrh,kbbcch,tcch(i),tcch(j),ebbcch)
                abbhch(k) = ebbhch
                abbcch(k) = ebbcch
            end do
        end do
    end if

    if(lreact) then
        do i=1,7
            tstr(i) = astr(i)
        end do
        do i=1,6
            tbhch(i)  = abhch(i)
            tbcch(i)  = abcch(i)
            tsbcch(i) = asbcch(i)
            tbbhch(i) = abbhch(i)
            tbbcch(i) = abbcch(i)
        end do
        do i=1,9
            tthcch(i)  = athcch(i)
            tsthcch(i) = asthcch(i)
            tvdwhh(i)  = avdwhh(i)
        end do

        !       remove one C-H2 str term
        tstr(1) = 0.0d0

        !       remove two H-C-H2 bend terms
        tbhch(1) = 0.0d0
        tbhch(2) = 0.0d0

        !       remove one C-C-H2 bend term
        tbcch(1) = 0.0d0

        !       remove one C-C-H2 strbnd term
        tsbcch(1) = 0.0d0

        !       remove three H-C-H2 bnd-bnd terms
        tbbhch(1) = 0.0d0
        tbbhch(2) = 0.0d0
        tbbhch(3) = 0.0d0

        !       remove two C-C-H2 bnd-bnd terms
        tbbcch(1) = 0.0d0
        tbbcch(2) = 0.0d0

        !       remove three H-C-C-H2 xtorsional terms
        tthcch(1) = 0.0d0
        tthcch(2) = 0.0d0
        tthcch(3) = 0.0d0

        !       remove three H-C-C-H2 str-xtorsional terms
        tsthcch(1) = 0.0d0
        tsthcch(2) = 0.0d0
        tsthcch(3) = 0.0d0

        !       remove three vdw terms
        tvdwhh(1) = 0.0d0
        tvdwhh(2) = 0.0d0
        tvdwhh(3) = 0.0d0

        xsumrx = 0.0d0
        xstr   = 0.0d0
        do i=1,7
            xstr = xstr+tstr(i)
        end do
        xbend   = 0.0d0
        xstrbnd = 0.0d0
        xbndbnd = 0.0d0
        do i=1,6
            xbend   = xbend+tbhch(i)+tbcch(i)
            xstrbnd = xstrbnd+tsbcch(i)
            xbndbnd = xbndbnd+tbbhch(i)+tbbcch(i)
        end do
        xtor    = 0.0d0
        xstrtor = 0.0d0
        xvdw    = 0.0d0
        do i=1,9
            xtor    = xtor+tthcch(i)
            xstrtor = xstrtor+tsthcch(i)
            xvdw    = xvdw+tvdwhh(i)
        end do
        xsumrx = xstr+xbend+xstrbnd+xbndbnd+xtor+xstrtor+xvdw

        !       call xdisp(xstr,xbend,xstrbnd,xbndbnd,
        !    1             xtor,xstrtor,xvdw,xsumRX)

        etotalmm = xsumrx
    end if

    ! start of if-block label: jg_case
    !....................................
    if(ljg) then
        do i=1,7
            tstr(i) = astr(i)
        end do
        do i=1,6
            tbhch(i)  = abhch(i)
            tbcch(i)  = abcch(i)
            tsbcch(i) = asbcch(i)
            tbbhch(i) = abbhch(i)
            tbbcch(i) = abbcch(i)
        end do
        do i=1,9
            tthcch(i)  = athcch(i)
            tsthcch(i) = asthcch(i)
            tvdwhh(i)  = avdwhh(i)
        end do

        !       remove the following C-H str terms
        !       C-H2
        !       C-H3
        !       C-H4
        tstr(1) = 0.0d0
        tstr(2) = 0.0d0
        tstr(3) = 0.0d0

        !       remove the following H-C-H bend terms
        !       H2-C1-H3
        !       H2-C1-H4
        !       H3-C1-H4
        tbhch(1) = 0.0d0
        tbhch(2) = 0.0d0
        tbhch(3) = 0.0d0

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
        tbbhch(1) = 0.0d0
        tbbhch(2) = 0.0d0
        tbbhch(3) = 0.0d0

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
        !
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

        xsumjg = 0.0d0
        xstr   = 0.0d0
        do i=1,7
            xstr = xstr+tstr(i)
        end do
        xbend   = 0.0d0
        xstrbnd = 0.0d0
        xbndbnd = 0.0d0
        do i=1,6
            xbend   = xbend+tbhch(i)+tbcch(i)
            xstrbnd = xstrbnd+tsbcch(i)
            xbndbnd = xbndbnd+tbbhch(i)+tbbcch(i)
        end do
        xtor    = 0.0d0
        xstrtor = 0.0d0
        xvdw    = 0.0d0
        do i=1,9
            xtor    = xtor+tthcch(i)
            xstrtor = xstrtor+tsthcch(i)
            xvdw    = xvdw+tvdwhh(i)
        end do
        xsumjg = xstr+xbend+xstrbnd+xbndbnd+xtor+xstrtor+xvdw
        if(ldebugmm) then
            write(*,'(A,F14.8)')
            write(*,'(A,F14.8)') '----------------------------'
            write(*,'(A,F14.8)') 'Energy from MM3'
            write(*,'(A,F14.8)') '----------------------------'
            call xdisp(xstr,xbend,xstrbnd,xbndbnd, &
                & xtor,xstrtor,xvdw,xsumjg)

        end if

    end if
    ! end of if-block label: jg_case
    !....................................

    vmm3 = xsumjg

end subroutine ethane
! end sub ethane
!**************************************************
!> @brief Evaluate xdisp for the C2H7_VBMM potential.
subroutine xdisp( xstr, xbend, xstrbnd, xbndbnd, &
        & xtor, xstrtor, xvdw, xsum )
    !**************************************************
    implicit none

    double precision :: xstr
    double precision :: xbend
    double precision :: xstrbnd
    double precision :: xbndbnd
    double precision :: xtor
    double precision :: xstrtor
    double precision :: xvdw
    double precision :: xsum
    write(*,'(A,F14.4)')
    write(*,'(A,F14.4)') 'Energy in kcal/mol'
    write(*,'(A,F14.4)')'Stretch   :',xstr
    write(*,'(A,F14.4)')'Bend      :',xbend
    write(*,'(A,F14.4)')'Str-Bnd   :',xstrbnd
    write(*,'(A,F14.4)')'Bnd-Bnd   :',xbndbnd
    write(*,'(A,F14.4)')'Torsional :',xtor
    write(*,'(A,F14.4)')'Str-Tor   :',xstrtor
    write(*,'(A,F14.4)')'vdW       :',xvdw
    write(*,'(A,F14.4)')
    write(*,'(A,F14.4)')'Total MM eng.  :',xsum
end subroutine xdisp
! end sub xdisp
!**************************************************
!> @brief Evaluate xdist for the C2H7_VBMM potential.
function xdist( x, y )
    !**************************************************
    implicit none

    double precision :: x(3)
    double precision :: y(3)
    double precision :: xdist
    integer          :: i
    double precision :: z
    z = 0.0d0
    do i=1,3
        z = z + ( (y(i)-x(i))**2 )
    end do
    xdist = dsqrt(z)

end function xdist
! end of func xdist
!**************************************************
!> @brief Evaluate xangle for the C2H7_VBMM potential.
function xangle( a, b, c )
    !**************************************************
    implicit none

    double precision :: a(3)
    double precision :: b(3)
    double precision :: c(3)
    double precision :: xangle
    double precision :: xdist

    integer          :: i
    double precision :: sum,rab(3),rcb(3)

    do i=1,3
        rab(i) = a(i)-b(i)
        rcb(i) = c(i)-b(i)
    end do
    sum = 0.0d0
    do i=1,3
        sum = sum + (rab(i)*rcb(i))
    end do
    sum    = sum/(xdist(a,b)*xdist(c,b))
    xangle = dabs(dacos(sum))

end function xangle
! end of xangle

!***********************************************************
!> @brief Evaluate xtorsion for the C2H7_VBMM potential.
function xtorsion( a, b, c, d )
    !***********************************************************
    implicit none

    double precision :: a(3)
    double precision :: b(3)
    double precision :: c(3)
    double precision :: d(3)
    double precision :: xtorsion
    !
    ! This func calculates the dihedral xangle
    ! between a-b and c-d
    ! a
    ! |
    ! b---c
    !     |
    !     d
    ! we will use the cross product definition of the
    ! dihedral xangle
    !
    integer          :: i,j
    double precision :: x,y,z
    double precision :: s
    double precision :: phi
    double precision :: dpuv
    double precision :: rab(3)
    double precision :: rcb(3)
    double precision :: rcd(3)
    double precision :: u(3)
    double precision :: v(3)
    !
    ! define bond-vectors
    do i=1,3
        rab(i) = a(i)-b(i)
        rcb(i) = c(i)-b(i)
        rcd(i) = c(i)-d(i)
    end do
    !
    ! define the normal vectors
    ! u = rab X rcb
    ! v = rcb X rcd
    call crossprod(rab,rcb,u)
    call crossprod(rcb,rcd,v)
    !
    ! evaluate the norm and dot product
    ! u and v
    ! dpuv = u . v
    ! s = rab . v and is used for sign of
    ! the dihedral xangle
    !
    x    = 0.0d0
    y    = 0.0d0
    s    = 0.0d0
    dpuv = 0.0d0
    do i=1,3
        x    = x+(u(i)*u(i))
        y    = y+(v(i)*v(i))
        dpuv = dpuv + (u(i)*v(i))
        s    = s + (rab(i)*v(i))
    end do
    x = sqrt(x)
    y = sqrt(y)
    !
    ! z is cos(phi) and should be in the
    ! range [-1,1]
    !
    if(x*y .ne. 0.0d0) then
        z = dpuv/(x*y)
    else
        z = dpuv
    end if
    if( z .gt.  1.0d0) z=1.0d0
    if( z .lt. -1.0d0) z=-1.0d0
    ! get dihedral xangle
    phi = dacos(z)
    !
    ! get the sign
    if(s /= 0.0d0) then
        s = s/dabs(s)
    else
        s = 1.0d0
    end if

    xtorsion = s * phi

end function xtorsion
! end of sub xtorsion

!**************************************************
!> @brief Evaluate crossprod for the C2H7_VBMM potential.
subroutine crossprod( a, b, c )
    !**************************************************
    implicit none

    double precision :: a(3)
    double precision :: b(3)
    double precision :: c(3)
    c(1) = (a(2)*b(3))-(a(3)*b(2))
    c(2) = (a(3)*b(1))-(a(1)*b(3))
    c(3) = (a(1)*b(2))-(a(2)*b(1))
end subroutine crossprod

!23456789
!*****************************************
!> @brief Evaluate Es for the C2H7_VBMM potential.
subroutine es( r0, ks, r, eng )
    !*****************************************
    implicit none

    double precision :: r0
    double precision :: ks
    double precision :: r
    double precision :: eng

    double precision            :: dr
    double precision, parameter :: c1 = 71.940d0
    double precision, parameter :: c2 = 2.550d0
    double precision, parameter :: c3 = 7.0d0/12.0d0

    eng = 0.0d0
    dr  = (r-r0)
    eng = (c1*ks*dr*dr) * (1.0d0-(c2*dr)+(c3*c2*c2*dr*dr))
    !    1 (1.0d0-(c2*dr)+(c3*c2*dr*dr))

end subroutine es
! end of sub Es
!*****************************************
!> @brief Evaluate Eb for the C2H7_VBMM potential.
subroutine eb( t0, kt, t, eng )
    !*****************************************
    implicit none

    double precision :: t0
    double precision :: kt
    double precision :: t
    double precision :: eng

    double precision :: dt
    ! Alligers (JACS, vol 111, pg 8551) parameters
    !     parameter(c1=0.0219140d0)
    !     parameter(c2=0.0140d0)
    !     parameter(c3=5.60d-5)
    !     parameter(c4=7.0d-7)
    !     parameter(c5=9.0d-10)
    ! Tinker 4.1 MM3 parameter
    double precision, parameter :: c1 = 0.021914180d0
    double precision, parameter :: c2 = 0.0140d0
    double precision, parameter :: c3 = 0.0000560d0
    double precision, parameter :: c4 = 0.00000070d0
    double precision, parameter :: c5 = 0.0000000220d0

    dt  = t-t0
    eng = c1*kt*dt*dt* ( 1.0d0-(c2*dt)+(c3*dt*dt)-(c4*dt*dt*dt) +(c5*dt*dt*dt*dt) )

end subroutine eb
! end of sub Eb
!*****************************************
!> @brief Evaluate Et for the C2H7_VBMM potential.
subroutine et( v1, v2, v3, w, eng )
    !*****************************************
    implicit none

    double precision :: v1
    double precision :: v2
    double precision :: v3
    double precision :: w
    double precision :: eng

    eng = 0.50d0*( ( v1*(1.0d0+dcos(w)) ) + ( v2*(1.0d0-dcos(2.0d0*w)) ) + ( v3*(1.0d0+dcos(3.0d0*w)) ) )

end subroutine et
! end of Et
!*****************************************
!> @brief Evaluate Esb for the C2H7_VBMM potential.
subroutine esb( ra0, rb0, t0, kst, ra, rb, t, eng )
    !*****************************************
    implicit none

    double precision :: ra0
    double precision :: rb0
    double precision :: t0
    double precision :: kst
    double precision :: ra
    double precision :: rb
    double precision :: t
    double precision :: eng

    double precision, parameter :: c1 = 2.511180d0

    eng = c1*kst*(t-t0)*((ra-ra0)+(rb-rb0))

end subroutine esb
! end of Esb
!*****************************************
!> @brief Evaluate Est for the C2H7_VBMM potential.
subroutine est( r0, kws, r, w, eng )
    !*****************************************
    implicit none

    double precision :: r0
    double precision :: kws
    double precision :: r
    double precision :: w
    double precision :: eng

    double precision, parameter :: c1 = -11.9950d0
    ! In Alligers paper(JACS, vol 111, pg 8551)
    ! c1 is +ve but TINKER (version 4.1) defines
    ! c1/2 as -ve

    eng = 0.50d0*c1*kws*(r-r0)*(1.0d0+dcos(3.0d0*w))

end subroutine est
! end of Est
!*****************************************
!> @brief Evaluate Ebb for the C2H7_VBMM potential.
subroutine ebb( ta0, tb0, kbb, ta, tb, eng )
    !*****************************************
    implicit none

    double precision :: ta0
    double precision :: tb0
    double precision :: kbb
    double precision :: ta
    double precision :: tb
    double precision :: eng

    double precision, parameter :: c1 = -0.021914180d0

    eng = c1*kbb*(ta-ta0)*(tb-tb0)
    !     write(*,1000) c1,kbb,(ta-ta0)*(tb-tb0),eng
end subroutine ebb
! end of Ebb
!*****************************************
!> @brief Evaluate Evdw for the C2H7_VBMM potential.
subroutine evdw( rv, esp, r, eng )
    !*****************************************
    implicit none

    double precision :: rv
    double precision :: esp
    double precision :: r
    double precision :: eng

    double precision            :: x
    double precision, parameter :: c1 = 2.250d0
    double precision, parameter :: c2 = 1.840d+5
    double precision, parameter :: c3 = 12.00d0

    x   = rv/r
    eng = esp*((-c1*x*x*x*x*x*x)+(c2*dexp(-c3/x)))

end subroutine evdw
! end of vdw

!23456789
! References:
! [1]. Joseph et al.   : JCP, vol 87(12),  pg. 7036, 1987
! [2]. Jordan et al.   : JCP, vol 102(14), pg. 5669, 1995
! [3]. Steckler et al. : JCP, vol 87(15),  pg. 7024, 1987
!

!*****************************************************************
!> @brief Evaluate xch5pot for the C2H7_VBMM potential.
subroutine xch5pot( np, para, rxncrd, c1, h2, h3, h4, h5, hb, xpot_kcal )
    !*****************************************************************
    !
    ! Input : Coordintes are in Ang.
    ! Output: xpot_kcal is in kcal/mol
    !
    ! All energy variables local to this
    ! routine are in 10(5)Joules/mol
    !
    implicit none

    ! Input variables
    integer          :: np
    double precision :: para(np)
    double precision :: rxncrd
    double precision :: c1(3)
    double precision :: h2(3)
    double precision :: h3(3)
    double precision :: h4(3)
    double precision :: h5(3)
    double precision :: hb(3)
    ! Output variables
    double precision :: xpot_kcal
    ! Parameter list
    double precision :: a1ch
    double precision :: b1ch
    double precision :: c1ch
    double precision :: alpcb
    double precision :: d1cb
    double precision :: d3cb
    double precision :: r0cb
    double precision :: alphh
    double precision :: d1hh
    double precision :: d3hh
    double precision :: r0hh
    double precision :: alpch
    double precision :: d1ch
    double precision :: d3ch
    double precision :: r0ch
    !
    double precision :: tau
    double precision :: aphi
    double precision :: bphi
    double precision :: cphi
    double precision :: atheta
    double precision :: btheta
    double precision :: ctheta
    double precision :: a3s
    double precision :: b3s
    double precision :: fd0
    double precision :: hd0
    !
    double precision :: a1s
    double precision :: a2s
    double precision :: b1s
    double precision :: b2s
    double precision :: ak
    double precision :: kch3
    double precision :: aa1
    double precision :: aa2
    double precision :: aa3
    double precision :: aa4
    !
    double precision, parameter :: fact1 = 0.041840d0
    double precision, parameter :: fact2 = 6.022045d0
    ! Funct called
    double precision :: xangle
    ! Local variables
    !     Local parameters
    integer, parameter          :: iout       = 6
    logical, parameter          :: ldebug     = .false.
    double precision, parameter :: kj_in_kcal = 0.2390057
    !
    !     bond lenght and bond vectors
    !
    double precision :: rch2
    double precision :: rch3
    double precision :: rch4
    double precision :: rch5
    double precision :: rbh2
    double precision :: rbh3
    double precision :: rbh4
    double precision :: rbh5
    double precision :: rcb
    double precision :: rchav
    double precision :: bvch2(3)
    double precision :: bvch3(3)
    double precision :: bvch4(3)
    double precision :: bvch5(3)
    double precision :: bvcb(3)
    !
    !     Vstr varaibles
    !
    double precision :: cbq
    double precision :: cbj
    double precision :: bh2q
    double precision :: bh2j
    double precision :: bh3q
    double precision :: bh3j
    double precision :: bh4q
    double precision :: bh4f
    double precision :: bh5q
    double precision :: bh5j
    double precision :: ch2q
    double precision :: ch2j
    double precision :: ch3q
    double precision :: ch3j
    double precision :: ch4q
    double precision :: ch4j
    double precision :: ch5q
    double precision :: ch5j
    double precision :: vqh2
    double precision :: vqh3
    double precision :: vqh4
    double precision :: vqh5
    double precision :: vjh2
    double precision :: vjh3
    double precision :: vjh4
    double precision :: vjh5
    double precision :: v3h2
    double precision :: v3h3
    double precision :: v3h4
    double precision :: v3h5
    double precision :: vstr
    !
    !     Bond xangles and out-of-plane-bendings
    !
    double precision :: pi
    double precision :: theta023
    double precision :: theta024
    double precision :: theta025
    double precision :: theta032
    double precision :: theta034
    double precision :: theta035
    double precision :: theta042
    double precision :: theta043
    double precision :: theta045
    double precision :: theta052
    double precision :: theta053
    double precision :: theta054
    double precision :: delta23
    double precision :: delta24
    double precision :: delta25
    double precision :: delta32
    double precision :: delta34
    double precision :: delta35
    double precision :: delta42
    double precision :: delta43
    double precision :: delta45
    double precision :: delta52
    double precision :: delta53
    double precision :: delta54
    double precision :: d23
    double precision :: d24
    double precision :: d25
    double precision :: d32
    double precision :: d34
    double precision :: d35
    double precision :: d42
    double precision :: d43
    double precision :: d45
    double precision :: d52
    double precision :: d53
    double precision :: d54
    double precision :: fd2
    double precision :: fd3
    double precision :: fd4
    double precision :: fd5
    double precision :: hd2
    double precision :: hd3
    double precision :: hd4
    double precision :: hd5
    double precision :: voph2
    double precision :: voph3
    double precision :: voph4
    double precision :: voph5
    double precision :: vop
    !
    !     in-plane bending parameters
    !
    double precision :: xk023
    double precision :: xk024
    double precision :: xk025
    double precision :: xk034
    double precision :: xk035
    double precision :: xk045
    double precision :: xk2
    double precision :: xk3
    double precision :: xk4
    double precision :: xk5
    double precision :: theta23
    double precision :: theta24
    double precision :: theta25
    double precision :: theta34
    double precision :: theta35
    double precision :: theta45
    double precision :: vip23
    double precision :: vip24
    double precision :: vip25
    double precision :: vip34
    double precision :: vip35
    double precision :: vip45
    double precision :: vip
    !
    !     misc and temp variables
    !
    integer          :: i,j,k
    double precision :: t1
    double precision :: t2
    double precision :: t3
    double precision :: t4
    double precision :: t5
    double precision :: t6
    double precision :: z1
    double precision :: z2
    double precision :: z3
    double precision :: z4
    double precision :: scat0
    double precision :: xpot
    !
    !     CH5 surf
    !
    double precision :: vsurf
    double precision :: r(18)
    double precision :: dv(18)
    !
    ! Assign parameters
    !
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
    !
    ! Define j,k,l indices wrt to each of the 4 H
    ! H2 : 3,4,5
    ! H3 : 5,4,2
    ! H4 : 5,2,3
    ! H5 : 2,4,3
    !
    ! Calculate bond lenghts and bond vectors
    ! bv is the bond vectors and rch is the
    ! corresponding bond lenghts
    ! rchav is the avg of the four ch xdist
    !
    rch2 = 0.0d0
    rch3 = 0.0d0
    rch4 = 0.0d0
    rch5 = 0.0d0
    rbh2 = 0.0d0
    rbh3 = 0.0d0
    rbh4 = 0.0d0
    rbh5 = 0.0d0
    rcb  = 0.0d0
    do i=1,3
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
    end do
    rch2  = dsqrt(rch2)
    rch3  = dsqrt(rch3)
    rch4  = dsqrt(rch4)
    rch5  = dsqrt(rch5)
    rcb   = dsqrt(rcb)
    rbh2  = dsqrt(rbh2)
    rbh3  = dsqrt(rbh3)
    rbh4  = dsqrt(rbh4)
    rbh5  = dsqrt(rbh5)
    rchav = 0.250d0*(rch2+rch3+rch4+rch5)
    !
    ! use switch func to change r0ch from
    ! ethane to ethyl value
    !
    !####
    t1 = 1.09250d0
    t2 = 1.08850d0
    t3 = c1ch*10.0d0
    t4 = 1.12750d0
    t5 = 0.50d0*(1.0d0+dtanh(t3*(rchav-t4)))
    !     r0ch=(t5*t2)+((1.0d0-t5)*t1)

    !
    ! calculate the Morse parameter alphaCH (alpch)
    ! used in the Vstr using Equation 2 of
    ! Ref. [2]
    !
    t1    = 0.50d0*( 1.0d0+dtanh(c1ch*(rchav-r0ch)) )
    alpch = a1ch+(b1ch*t1)
    !
    ! calculate Vstr
    ! first get the coulomb and resonance integrals
    ! involving Hb and all the other atoms
    ! cbq and cbj : coulomb and resonance integral between C1 and Hb
    !
    call v3qj77(alpcb,d1cb,d3cb,r0cb,rcb,cbq,cbj)
    call v3qj77(alphh,d1hh,d3hh,r0hh,rbh2,bh2q,bh2j)
    call v3qj77(alphh,d1hh,d3hh,r0hh,rbh3,bh3q,bh3j)
    call v3qj77(alphh,d1hh,d3hh,r0hh,rbh4,bh4q,bh4f)
    call v3qj77(alphh,d1hh,d3hh,r0hh,rbh5,bh5q,bh5j)
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
    !
    ! C1 and (H2...H5) terms
    !
    !####
    t1 = 0.50d0*(1.0d0+dtanh(10.0d0*(rbh2-1.50d0)))
    t2 = 1.0d0-t1
    t3 = (1.0d0+(t2*0.080d0))*r0ch
    !     call V3qj77(alpch,d1ch,d3ch,t3,rch2,ch2q,ch2j)

    call v3qj77(alpch,d1ch,d3ch,r0ch,rch2,ch2q,ch2j)
    call v3qj77(alpch,d1ch,d3ch,r0ch,rch3,ch3q,ch3j)
    call v3qj77(alpch,d1ch,d3ch,r0ch,rch4,ch4q,ch4j)
    call v3qj77(alpch,d1ch,d3ch,r0ch,rch5,ch5q,ch5j)
    !
    ! coulomb and resonance integral for H2...H5
    ! Equation 2 of Ref. [1]
    !
    vqh2 = cbq+ch2q+bh2q
    vqh3 = cbq+ch3q+bh3q
    vqh4 = cbq+ch4q+bh4q
    vqh5 = cbq+ch5q+bh5q

    t2 = ((ch2j-cbj)**2) + ((cbj-bh2j )**2) + ((bh2j-ch2j)**2)

    t3 = ((ch3j-cbj)**2) + ((cbj-bh3j )**2) + ((bh3j-ch3j)**2)

    t4 = ((ch4j-cbj)**2) + ((cbj-bh4f )**2) + ((bh4f-ch4j)**2)

    t5 = ((ch5j-cbj)**2) + ((cbj-bh5j )**2) + ((bh5j-ch5j)**2)

    vjh2 = -dsqrt(0.50d0*t2)
    vjh3 = -dsqrt(0.50d0*t3)
    vjh4 = -dsqrt(0.50d0*t4)
    vjh5 = -dsqrt(0.50d0*t5)

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

    t1 = 6.0d0
    t2 = 0.950d0
    t3 = 0.50d0*(1.0d0+dtanh(t1*(rxncrd-t2)))
    t4 = 1.013850d0

    call xgaussfit(rxncrd,t4)
    !     t4=t4*0.999930d0
    !     t4=t4*1.00205d0
    !     t4=t4*1.0159995120d0/1.016043066097427340d0
    t4 = t4*1.016043066097427340d0/1.0159995120d0

    t5 = (1.0d0*(1.0d0-t3))+(t3*t4)

    z1 = 6.0d0
    z2 = 10.0d0
    z3 = 0.50d0*(1.0d0+dtanh(z1*(rxncrd-z2)))
    z4 = (t5*(1.0d0-z3))+(1.0d0*z3)

    vstr = vstr * z4

    !     vstr = vstr * 1.013850d0
    !%%%--------------------------------------------------

    if(ldebug) then
        write(*,*)
        write(*,*)'H2 Vstr :',v3h2*1.0d2*kj_in_kcal
        write(*,*)'H3 Vstr :',v3h3*1.0d2*kj_in_kcal
        write(*,*)'H4 Vstr :',v3h4*1.0d2*kj_in_kcal
        write(*,*)'H5 Vstr :',v3h5*1.0d2*kj_in_kcal
        write(*,*) 'VSTR FROM XPOT:',vstr*1.0d2*kj_in_kcal
    end if
    !
    ! calculate the out-of-plane bend terms
    !
    ! evaluate tau and PI
    pi  = dacos(-1.0d0)
    tau = dacos(-1.0d0/3.0d0)
    !
    ! calculate theta0
    !
    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch2,rch3,rch4,rch5,theta023)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch2,rch4,rch5,rch3,theta024)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch2,rch5,rch3,rch4,theta025)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch3,rch2,rch5,rch4,theta032)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch3,rch4,rch2,rch5,theta034)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch3,rch5,rch4,rch2,theta035)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch4,rch2,rch3,rch5,theta042)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch4,rch3,rch5,rch2,theta043)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch4,rch5,rch2,rch3,theta045)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch5,rch2,rch4,rch3,theta052)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch5,rch3,rch2,rch4,theta053)

    call gjtheta0(pi,tau,aphi,bphi,cphi, &
        & atheta,btheta,ctheta,r0ch, &
        & rch5,rch4,rch3,rch2,theta054)

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

    !
    ! calculate delta
    !
    call steckdel(bvch3,bvch4,bvch5,bvch3,d23)
    call steckdel(bvch3,bvch4,bvch5,bvch4,d24)
    call steckdel(bvch3,bvch4,bvch5,bvch5,d25)

    call steckdel(bvch5,bvch4,bvch2,bvch2,d32)
    call steckdel(bvch5,bvch4,bvch2,bvch4,d34)
    call steckdel(bvch5,bvch4,bvch2,bvch5,d35)

    call steckdel(bvch5,bvch2,bvch3,bvch2,d42)
    call steckdel(bvch5,bvch2,bvch3,bvch3,d43)
    call steckdel(bvch5,bvch2,bvch3,bvch5,d45)

    call steckdel(bvch2,bvch4,bvch3,bvch2,d52)
    call steckdel(bvch2,bvch4,bvch3,bvch3,d53)
    call steckdel(bvch2,bvch4,bvch3,bvch4,d54)

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

    if(ldebug) then
        write(iout,*)
        write(iout,*) 'theta0 and delta terms:'
        write(iout,*) 'H2-H3:',theta023,delta23
        write(iout,*) 'H2-H4:',theta024,delta24
        write(iout,*) 'H2-H5:',theta025,delta25
        write(iout,*) 'H3-H4:',theta034,delta34
        write(iout,*) 'H3-H5:',theta035,delta35
        write(iout,*) 'H4-H5:',theta045,delta45
    end if
    !
    ! calculate fd and hd terms
    call calcfh(a3s,b3s,fd0,hd0,r0ch,rch2,rch3,rch4,rch5,fd2,hd2)
    call calcfh(a3s,b3s,fd0,hd0,r0ch,rch3,rch5,rch4,rch2,fd3,hd3)
    call calcfh(a3s,b3s,fd0,hd0,r0ch,rch4,rch5,rch2,rch3,fd4,hd4)
    call calcfh(a3s,b3s,fd0,hd0,r0ch,rch5,rch2,rch4,rch3,fd5,hd5)

    if(ldebug) then
        write(iout,*)
        write(iout,*)'FD and HD terms:'
        write(iout,*)'H2',fd2,hd2
        write(iout,*)'H3',fd3,hd3
        write(iout,*)'H4',fd4,hd4
        write(iout,*)'H5',fd5,hd5
    end if
    !
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

    if(ldebug) then
        write(*,*)
        write(*,*)'H2 Vop : ',voph2*1.0d2*kj_in_kcal
        write(*,*)'H3 Vop : ',voph3*1.0d2*kj_in_kcal
        write(*,*)'H4 Vop : ',voph4*1.0d2*kj_in_kcal
        write(*,*)'H5 Vop : ',voph5*1.0d2*kj_in_kcal
        write(iout,*) 'VOP FROM XPOT :',vop*1.0d2*kj_in_kcal
    end if
    !
    ! Calculate in-plane terms

    ! get the k0ij terms
    call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch, &
        & rch2,rch3,rch4,rch5,xk023)
    call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch, &
        & rch2,rch4,rch5,rch3,xk024)
    call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch, &
        & rch2,rch5,rch3,rch4,xk025)

    !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
    !    1               rch3,rch2,rch5,rch4,xk032)
    call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch, &
        & rch3,rch4,rch2,rch5,xk034)
    call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch, &
        & rch3,rch5,rch4,rch2,xk035)

    !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
    !    1               rch4,rch2,rch3,rch5,xk042)
    !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
    !    1               rch4,rch3,rch5,rch2,xk043)
    call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch, &
        & rch4,rch5,rch2,rch3,xk045)

    !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,rch5,rch2,rch4,rch3,xk
    !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,rch5,rch3,rch2,rch4,xk
    !     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,rch5,rch4,rch3,rch2,xk

    !
    ! get ki terms
    call calcvipk(aa1,aa2,aa3,aa4,r0ch,r0hh,rch2,rbh2,xk2)
    call calcvipk(aa1,aa2,aa3,aa4,r0ch,r0hh,rch3,rbh3,xk3)
    call calcvipk(aa1,aa2,aa3,aa4,r0ch,r0hh,rch4,rbh4,xk4)
    call calcvipk(aa1,aa2,aa3,aa4,r0ch,r0hh,rch5,rbh5,xk5)
    !
    ! get bond xangles
    theta23 = xangle(h2,c1,h3)
    theta24 = xangle(h2,c1,h4)
    theta25 = xangle(h2,c1,h5)
    theta34 = xangle(h3,c1,h4)
    theta35 = xangle(h3,c1,h5)
    theta45 = xangle(h4,c1,h5)

    vip23 = xk023*xk2*xk3*(theta23-theta023)*(theta23-theta023)
    vip24 = xk024*xk2*xk4*(theta24-theta024)*(theta24-theta024)
    vip25 = xk025*xk2*xk5*(theta25-theta025)*(theta25-theta025)
    vip34 = xk034*xk3*xk4*(theta34-theta034)*(theta34-theta034)
    vip35 = xk035*xk3*xk5*(theta35-theta035)*(theta35-theta035)
    vip45 = xk045*xk4*xk5*(theta45-theta045)*(theta45-theta045)

    vip = 0.50d0*(vip23+vip24+vip25+vip34+vip35+vip45)

    !
    ! remove H5 str term
    ! total potential in 10(5)Joules/mol and kcal/mol
    xpot      = vstr+vop+vip-v3h5-voph5-((vip25+vip35+vip45)*0.50d0)
    xpot_kcal = xpot*100.0d0*kj_in_kcal

    if(ldebug) then
        write(*,*)
        write(*,*)'H2-H3 Vip :',0.50d0*vip23*1.0d2*kj_in_kcal
        write(*,*)'H2-H4 Vip :',0.50d0*vip24*1.0d2*kj_in_kcal
        write(*,*)'H2-H5 Vip :',0.50d0*vip25*1.0d2*kj_in_kcal
        write(*,*)'H3-H4 Vip :',0.50d0*vip34*1.0d2*kj_in_kcal
        write(*,*)'H3-H5 Vip :',0.50d0*vip35*1.0d2*kj_in_kcal
        write(*,*)'H4-H5 Vip :',0.50d0*vip45*1.0d2*kj_in_kcal
        write(*,*) 'VIP FROM XPOT:',vip*1.0d2*kj_in_kcal
        write(*,*) 'VSTR+VIP+VOP :',(vstr+vip+vop)*1.0d2*kj_in_kcal

        write(*,*)
        write(*,*)'Energy after removing H5 term'
        write(*,*)'Vstr:',(vstr-v3h5)*1.0d2*kj_in_kcal
        write(*,*)'Vop :',(vop-voph5)*1.0d2*kj_in_kcal
        write(*,*)'Vip :',vip-((vip25+vip35+vip45)*0.50d0)* &
            & 1.0d2*kj_in_kcal
        write(*,*)'Total xpot:',xpot_kcal
    end if

    if(.false.) then
        r(1)  = h2(1)/0.52918d0
        r(2)  = h2(2)/0.52918d0
        r(3)  = h2(3)/0.52918d0
        r(4)  = c1(1)/0.52918d0
        r(5)  = c1(2)/0.52918d0
        r(6)  = c1(3)/0.52918d0
        r(7)  = h3(1)/0.52918d0
        r(8)  = h3(2)/0.52918d0
        r(9)  = h3(3)/0.52918d0
        r(10) = h4(1)/0.52918d0
        r(11) = h4(2)/0.52918d0
        r(12) = h4(3)/0.52918d0
        r(13) = h5(1)/0.52918d0
        r(14) = h5(2)/0.52918d0
        r(15) = h5(3)/0.52918d0
        r(16) = hb(1)/0.52918d0
        r(17) = hb(2)/0.52918d0
        r(18) = hb(3)/0.52918d0
        call csurf(vsurf,r,dv,18)
        !       vsurf=(vsurf-v3h5)*100.0d0*KJ_IN_KCAL
        write(*,*) 'VSURF:',vsurf*627.510d0
        write(*,*) 'MYVSURF:',(vstr+vop+vip)*0.03812d0*627.510d0
    end if

end subroutine xch5pot
!23456789
! References:
! [1]. Joseph et al. : JCP, vol 87(12),  pg 7036, 1987
! [2]. Jordan et al. : JCP, vol 102(14), pg 5669, 1995
!
!************************************************
!> @brief Evaluate V3qj77 for the C2H7_VBMM potential.
subroutine v3qj77( alpha, d1, d3, r0, r, vq, vj )
    !************************************************
    ! this subroutine calculates the Coulomb (vq) and
    ! resonance (vj) integrals from singlet (e1) and
    ! triplet (E3) states
    ! Equations: 2, 3, 4, 5, and 8 of Ref. [1]
    !
    implicit none

    ! Input variables
    double precision :: alpha
    double precision :: d1
    double precision :: d3
    double precision :: r0
    double precision :: r
    ! Output variables
    double precision :: vq
    double precision :: vj
    ! Local variables
    double precision :: z1,z2
    double precision :: e1,e3

    z1 = dexp(-alpha*(r-r0))
    z2 = dexp(-2.0d0*alpha*(r-r0))

    e1 = d1*(z2-(2.0d0*z1))
    e3 = d3*(z2+(2.0d0*z1))
    vq = 0.50d0*(e1+e3)
    vj = 0.50d0*(e1-e3)
end subroutine v3qj77
! end of V3qj

!23456789
! References:
! [1]. Joseph et al. : JCP, vol 87(12),  pg. 7036, 1987
! [2]. Jordan et al. : JCP, vol 102(14), pg. 5669, 1995
! [3]. Steckler et al. : JCP, vol 87(15), pg. 7024, 1987
!
!************************************************
!> @brief Evaluate funcS3 for the C2H7_VBMM potential.
function funcs3( alpha, beta, r0, r )
    !************************************************
    ! This func. calculates the S3 switching func
    ! Equation 4 of Ref. [2]
    !
    implicit none

    double precision :: alpha
    double precision :: beta
    double precision :: r0
    double precision :: r
    double precision :: funcs3
    funcs3 = 1.0d0-dtanh(alpha*(r-r0)*(r-beta)*(r-beta))
end function funcs3
! end of func funcS3
!***********************************************************************
!> @brief Evaluate calcfh for the C2H7_VBMM potential.
subroutine calcfh( a3, b3, fd0, hd0, r0, ri, rj, rk, rl, fdi, hdi )
    !***********************************************************************
    ! This sub. calculates the fdelta_i and hdelta_i
    ! as defined by Equation 4 of Ref. [2]
    !
    implicit none

    ! Input variables
    double precision :: a3
    double precision :: b3
    double precision :: fd0
    double precision :: hd0
    double precision :: r0
    double precision :: ri
    double precision :: rj
    double precision :: rk
    double precision :: rl
    ! Output variables
    double precision :: fdi
    double precision :: hdi
    ! Local variables
    double precision :: s3i
    double precision :: s3j
    double precision :: s3k
    double precision :: s3l

    s3i = 1.0d0-tanh( a3*(ri-r0)*(ri-b3)*(ri-b3) )
    s3j = 1.0d0-tanh( a3*(rj-r0)*(rj-b3)*(rj-b3) )
    s3k = 1.0d0-tanh( a3*(rk-r0)*(rk-b3)*(rk-b3) )
    s3l = 1.0d0-tanh( a3*(rl-r0)*(rl-b3)*(rl-b3) )
    fdi = (1.0d0-s3i)*s3j*s3k*s3l*fd0
    hdi = (1.0d0-s3i)*s3j*s3k*s3l*hd0

end subroutine calcfh
! end of sub calcfh

!***********************************************************************
!23456789
!> @brief Evaluate gjtheta0 for the C2H7_VBMM potential.
subroutine gjtheta0( pi, tau, aphi, bphi, cphi, atheta, btheta, ctheta, &
        & r0, ri, rj, rk, rl, theta0ij )
    !***********************************************************************
    ! This sub. computes the theta0 for bending terms
    ! Equation 5 in Ref. [2]
    !
    implicit none

    ! Input variables
    double precision :: pi
    double precision :: tau
    double precision :: aphi
    double precision :: bphi
    double precision :: cphi
    double precision :: atheta
    double precision :: btheta
    double precision :: ctheta
    double precision :: r0
    double precision :: ri
    double precision :: rj
    double precision :: rk
    double precision :: rl
    ! Output variables
    double precision :: theta0ij
    ! Local variables
    double precision :: si
    double precision :: sj
    double precision :: sk
    double precision :: sl
    double precision :: s1
    double precision :: s2
    double precision :: t1
    double precision :: t2
    double precision :: sphii
    double precision :: sphij
    double precision :: sthetak
    double precision :: sthetal

    si       = aphi*(ri-r0)*dexp(bphi*(ri-cphi)*(ri-cphi)*(ri-cphi))
    sj       = aphi*(rj-r0)*dexp(bphi*(rj-cphi)*(rj-cphi)*(rj-cphi))
    sk       = atheta*(rk-r0)*dexp(btheta*(rk-ctheta)*(rk-ctheta)*(rk-ctheta))
    sl       = atheta*(rl-r0)*dexp(btheta*(rl-ctheta)*(rl-ctheta)*(rl-ctheta))
    sphii    = 1.0d0-dtanh(si)
    sphij    = 1.0d0-dtanh(sj)
    sthetak  = 1.0d0-dtanh(sk)
    sthetal  = 1.0d0-dtanh(sl)
    s1       = (sphii*sphij)-1.0d0
    s2       = (sthetak*sthetal)-1.0d0
    t1       = tau-(0.50d0*pi)
    t2       = tau-(2.0d0*pi/3.0d0)
    theta0ij = tau+(t1*s1)+(t2*s2)

end subroutine gjtheta0
! end of gjtheta0
!***********************************************************************
!> @brief Evaluate steckdel for the C2H7_VBMM potential.
subroutine steckdel( v1, v2, v3, vi, deltai )
    !***********************************************************************
    ! This sub. computes the out-of-plane delta variable
    ! as defined by Steckler et al. in Ref. [3] equation 18.
    ! Note that is *NOT* symmetric with respect to all
    ! the three methyl hydrogen.
    ! Symmetrization (as given in Ref. [2] equation 4)
    ! will be done in the calling suroutine.
    !
    implicit none

    ! Input variables
    double precision :: v1(3)
    double precision :: v2(3)
    double precision :: v3(3)
    double precision :: vi(3)
    ! Output variables
    double precision :: deltai
    ! Local variables
    integer          :: i
    double precision :: ri
    double precision :: rn
    double precision :: r21
    double precision :: r31
    double precision :: doti
    double precision :: v21(3)
    double precision :: v31(3)
    double precision :: pn(3)
    !
    ! r21, r31 and ri are the lenghts of vectors
    ! v21, v31 and vi resp.
    ! v21=v2-v1
    ! v31=v3-v1
    !
    do i=1,3
        v21(i) = v2(i)-v1(i)
        v31(i) = v3(i)-v1(i)
    end do
    !
    ! the normal to the plane containing
    ! points 1, 2, and 3 is given by
    ! pn = r21 X r31
    ! rn is lenght of pn
    !
    call crossprod(v21,v31,pn)
    !
    ! The xangle between pn and vi is deltaI
    ! evalute the dot product:
    ! doti = pn . vi
    !
    ri   = 0.0d0
    rn   = 0.0d0
    doti = 0.0d0
    do i=1,3
        ri   = ri+(vi(i)*vi(i))
        rn   = rn+(pn(i)*pn(i))
        doti = doti+(pn(i)*vi(i))
    end do
    ! doti should be -ve by convension
    ! for right-handed set
    doti = -dabs(doti)
    ri   = dsqrt(ri)
    rn   = dsqrt(rn)

    deltai = dacos(doti/(ri*rn))

end subroutine steckdel
! end of sub steckdel
!23456789
! References:
! [1]. Joseph et al. : JCP, vol 87(12),  pg. 7036, 1987
! [2]. Jordan et al. : JCP, vol 102(14), pg. 5669, 1995
! [3]. Steckler et al. : JCP, vol 87(15), pg. 7024, 1987
!
!***********************************************************************
!> @brief Evaluate calcvipk for the C2H7_VBMM potential.
subroutine calcvipk( a1, a2, a3, a4, r0ch, r0bh, rch, rbh, vipk )
    !***********************************************************************
    ! The sub. calculates the force constant (vipk) for the in-plane
    ! bendinging as given by equation 6 of Ref. [2]
    !
    implicit none

    ! Input varaibles
    double precision :: a1
    double precision :: a2
    double precision :: a3
    double precision :: a4
    double precision :: r0ch
    double precision :: r0bh
    double precision :: rch
    double precision :: rbh
    ! Output variables
    double precision :: vipk
    ! Local variables
    double precision :: ua1
    double precision :: ua2

    ua1  = 1.0d0-dexp(-a1*rbh*rbh)
    ua2  = a2+(a3*dexp(-a4*(rbh-r0bh)*(rbh-r0bh)))
    vipk = ua1*dexp(-ua2*(rch-r0ch)*(rch-r0ch))
end subroutine calcvipk
! end of sub calcvipk
!***********************************************************************
!> @brief Evaluate calcvipk0 for the C2H7_VBMM potential.
subroutine calcvipk0( a1, b1, a2, b2, ak, kch3, r0ch, ri, rj, rk, rl, vipk0ij )
    !***********************************************************************
    ! This sub. calculates the equilibrium force constant (k0ij) used
    ! for in-plane-bending terms
    ! Equation 9 of Ref. [2]
    implicit none

    ! Input varaibles
    double precision :: a1
    double precision :: b1
    double precision :: a2
    double precision :: b2
    double precision :: kch3
    double precision :: ak
    double precision :: r0ch
    double precision :: ri
    double precision :: rj
    double precision :: rk
    double precision :: rl
    ! Output varaibles
    double precision :: vipk0ij
    ! Local varaibles
    double precision :: kch4
    double precision :: s1i
    double precision :: s1j
    double precision :: s2k
    double precision :: s2l

    kch4 = kch3+ak
    ! calculate switching func
    s1i = 1.0d0-dtanh( a1*(ri-r0ch)*((ri-b1)**8) )
    s1j = 1.0d0-dtanh( a1*(rj-r0ch)*((rj-b1)**8) )
    s2k = 1.0d0-dtanh( a2*(rk-r0ch)*((rk-b2)**6) )
    s2l = 1.0d0-dtanh( a2*(rl-r0ch)*((rl-b2)**6) )

    ! calculate k0 using Eq. 9 of Ref. [2]
    vipk0ij = kch4+( kch4*((s1i*s1j)-1.0d0) ) +( (kch4-kch3)*((s2k*s2l)-1.0d0) )

end subroutine calcvipk0
! end of sub. calcvipk0
!23456789
!************************************************
!> @brief Evaluate xlimit for the C2H7_VBMM potential.
subroutine xlimit( a, b, s, x0, x, y )
    !************************************************
    ! This subroutine uses a tanh function
    ! to smoothly connect two points (a,b)
    ! The tanh function is centered at x0
    ! with a scale factor of s.
    !
    ! The output value is given by y
    !
    implicit none

    double precision :: a
    double precision :: b
    double precision :: s
    double precision :: x0
    double precision :: x
    double precision :: y

    double precision :: r
    double precision :: ca
    double precision :: cb
    !---------------------
    ! When s is +ve:
    !---------------------
    ! As x-->b
    ! cb--->1 & ca-->0
    !
    ! As x-->a
    ! ca--->1 & cb-->0
    !
    r  = x-x0
    cb = 0.50d0*(1.0d0+dtanh(s*r))
    ca = 1.0d0-cb
    y  = (ca*a)+(cb*b)

end subroutine xlimit
! end of xlimit
!==========================================================
!==========================================================
!==========================================================

!
! System:          CH5
! Functional form:
! Common name:     JG (2nd coding)
! Reference  :      M. J. T. Jordan, R. G. Gilbert, J. Chem. Phys.,
!                   Vol. 102, p. 5669, 1995.
!
! Number of bodies: 6
! Number of electronic states: 1
! Number of derivatives: 0
! Interface: Section-2
! Notes: This routine contains the same JG CH5 surface as the routine 'c
!        but has a different interface.

!> @brief Evaluate CSURF for the C2H7_VBMM potential.
subroutine csurf( engygs, r, degsdr, n3tm )
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

    double precision :: a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, aphi, atheta, b1ch, b1s, b2s, b3s, bk, &
        & bphi, btheta, c1ch, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, degsdr, dfdc, dfdelta, dfdh, dhdelta, &
        & ds1, ds2, ds3, dsphi, dstheta, dtheta0, en, engygs, f1, fch3, fdelta, fk0, fkinf, hch3, hdelta, pdot, q, r, &
        & r0cb, r0ch, r0hh, rbh, rcb, rch, s1, s2, s3, sphi, stheta, tbh, tcb, tch, theta0, vip, vop, vstr
    integer :: i, n3tm, nc, nh, nhb, nnb, nnc, nnh

    double precision, parameter :: pi = 3.141592653589793d0
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    dimension r(n3tm),degsdr(n3tm)
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)

    call prepot
    ! Changing units from angstrom to bohr and initialize
    !
    do 10 i = 1, 18
        q(i)    = r(i)*0.52918d0
        pdot(i) = 0.d0
    10 continue
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
    !  calculate reference xangles and their derivatives
    !
    call refxangles
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
    !##--------------------------------------------------
    !       write(*,*)
    !       write(*,*) 'vstr 10(5)j/mol and au : ',vstr,vstr*0.03812D0
    !       write(*,*) 'vop  10(5)j/mol and au : ',vop,vop*0.03812D0
    !       write(*,*) 'vip  10(5)j/mol and au : ',vip,vip*0.03812D0
    !       write(*,*) 'pot  10(5)j/mol and au : ',en,en*0.03812D0
    !##--------------------------------------------------
    en     = en*0.03812d0
    engygs = en
    !
    ! 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
    !
    do i=1,18
        degsdr(i) = pdot(i)*0.0201723d0
    enddo
    !
end subroutine csurf
!
!******************************************************
!
!> @brief Evaluate coorden for the C2H7_VBMM potential.
subroutine coorden
    !
    !  calculates relative coordinates and bond lengths
    !
    implicit none

    double precision :: a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, atheta, b1ch, b1s, b2s, &
        & b3s, bk, bphi, btheta, c1ch, cart, cartnu, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, deesdr, degsdr, &
        & deijdr, descart, dfdc, dfdelta, dfdh, dgscart, dhdelta, dijcart, ds1, ds2, ds3, dsphi, dstheta, dtheta0, &
        & engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, fkinf, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, &
        & r, r0cb, r0ch, r0hh, rbh, rcb, rch, s1, s2, s3, sphi, stheta, tbh, tcb, tch, theta0
    integer :: i, icartr, ind, indexes, irctnt, mder, msurf, nasurf, &
        & natoms, nc, nder, nflag, nh, nhb, nnb, nnc, nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    !  calculate relative coordinates
    !
    do ind=1,3
        tcb(ind) = q(nc(ind))-q(nhb(ind))
        do i=1,4
            tch(i,ind) = q(nc(ind))-q(nh(i,ind))
            tbh(i,ind) = q(nhb(ind))-q(nh(i,ind))
        enddo
    enddo
    !
    !  calculate bond lengths
    !
    rcb = sqrt(tcb(1)*tcb(1)+tcb(2)*tcb(2)+tcb(3)*tcb(3))
    do i=1,4
        rch(i) = sqrt(tch(i,1)*tch(i,1)+tch(i,2)*tch(i,2)+ tch(i,3)*tch(i,3))
        rbh(i) = sqrt(tbh(i,1)*tbh(i,1)+tbh(i,2)*tbh(i,2)+ tbh(i,3)*tbh(i,3))
    enddo
end subroutine coorden
!
!******************************************************
!
!
!> @brief Evaluate refxangles for the C2H7_VBMM potential.
subroutine refxangles
    !
    !  subroutine calculates reference xangles for the "in-plane" potential
    !
    implicit none

    double precision :: a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, atheta, b1ch, b1s, b2s, &
        & b3s, bk, bphi, btheta, c1ch, cart, cartnu, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, deesdr, degsdr, &
        & deijdr, descart, dfdc, dfdelta, dfdh, dgscart, dhdelta, dijcart, ds1, ds2, ds3, dsphi, dstheta, dtheta0, &
        & engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, fkinf, halfpi, hch3, hdelta, pdot, pengyes, pengygs, &
        & pengyij, q, r, r0cb, r0ch, r0hh, rbh, rcb, rch, s1, s2, s3, sphi, stheta, tau, tbh, tcb, tch, theta0, twopi
    integer :: i, icartr, indexes, irctnt, j, k, mder, msurf, nasurf, &
        & natoms, nc, nder, nflag, nh, nhb, nnb, nnc, nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    tau = acos(-1.0d0/3.0d0)
    !       pi=4.0d0*atan(1.0d0)
    halfpi = 0.5d0*pi
    twopi  = 2.0d0*pi
    !
    !  set diagonal elements to zero
    !
    do i=1,4
        theta0(i,i) = 0.0d0
        do k=1,4
            dtheta0(i,i,k) = 0.0d0
        enddo
    enddo
    !
    !  calculate reference xangles
    !
    !##
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
    do i=1,3
        do j=i+1,4
            theta0(j,i) = theta0(i,j)
            do k=1,4
                dtheta0(j,i,k) = dtheta0(i,j,k)
            enddo
        enddo
    enddo
end subroutine refxangles
!
!******************************************************
!
!
!> @brief Evaluate stretch for the C2H7_VBMM potential.
subroutine stretch( vstr )
    !
    !  subroutine to calculate leps-type stretching potential and its
    !  derivatives
    !
    implicit none

    double precision :: a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ach, achdc, achdh, ahh, ak, anuzero, aphi, arga, &
        & atheta, b1ch, b1s, b2s, b3s, bk, bphi, btheta, c1ch, cart, cartnu, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, &
        & d3hh, deesdr, degsdr, deijdr, descart, dfdc, dfdelta, dfdh, dgscart, dhdelta, dijcart, ds1, ds2, ds3, dsphi, &
        & dstheta, dtheta0, dumach, dumjbh, dumjcb, dumjch, dumjhi, dumjhj, dumqbh, dumqcb, dumqch, dumqhi, dumqhj, e1, &
        & e3, engyes, engygs, engyij, ezero, f1, factj, fch3, fdelta, fk0, fkinf, hch3, hdelta, pdot, pengyes, pengygs, &
        & pengyij, q, r, r0cb, r0ch, r0hh, rav, rbh, rcb, rch, s1, s2, s3, sphi, stheta, tbh, tcb, tch, theta0, vj, vjbh, &
        & vjcb, vjch, vq, vqbh, vqcb, vqch, vstr
    integer :: i, icartr, ind, indexes, irctnt, j, k, mder, msurf, nasurf, natoms, nc, nder, nflag, nh, nhb, nnb, nnc, &
        & nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    dimension vqch(4),vjch(4),vqbh(4),vjbh(4),vq(4),vj(4), &
        & achdc(3),achdh(4,3)
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
    if(arga.lt.19.0d0)then
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
    do i=1,4
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
    do ind=1,3
        achdc(ind) = dumach*(tch(1,ind)/rch(1)+tch(2,ind)/rch(2) +tch(3,ind)/rch(3)+tch(4,ind)/rch(4))/4.0d0
        do i=1,4
            achdh(i,ind) = -dumach*tch(i,ind)/rch(i)/4.0d0
        enddo
    enddo
    dumqcb = -acb*((d1cb+d3cb)*exp(-2.0d0*acb*(rcb-r0cb))- (d1cb-d3cb)*exp(-acb*(rcb-r0cb)))/rcb
    !
    !  calculate cartesian derivatives:
    !  looping over ch(i) and bh(i)
    !
    do i=1,4
        dumqbh = -ahh*((d1hh+d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh))- (d1hh-d3hh)*exp(-ahh*(rbh(i)-r0hh)))/rbh(i)
        factj  = 0.5d0/vj(i)
        dumjcb = -acb*((d1cb-d3cb)*exp(-2.0d0*acb*(rcb-r0cb)) -(d1cb+d3cb)*exp(-acb*(rcb-r0cb)))*factj/rcb
        dumjbh = -ahh*((d1hh-d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh)) -(d1hh+d3hh)*exp(-ahh*(rbh(i)-r0hh)))*factj/rbh(i)
        do ind=1,3
            !
            !  deriv wrt hb:
            !
            pdot(nhb(ind))=pdot(nhb(ind)) -tcb(ind)*dumqcb+tbh(i,ind)*dumqbh +(vjch(i)-vjcb)*(dumjcb*tcb(ind)) +(vjcb- &
                & vjbh(i))*(-dumjcb*tcb(ind)-dumjbh*tbh(i,ind)) +(vjbh(i)-vjch(i))*dumjbh*tbh(i,ind)
            !
            !  dvqch(i)/dc
            !
            dumqch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch)) *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch)) - &
                & (d1ch- d3ch)*exp(-ach*(rch(i)-r0ch)))
            pdot(nc(ind)) = pdot(nc(ind))+dumqch+tcb(ind)*dumqcb
            !
            !  dvqch(i)/dh(i)
            !
            dumqhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch)) *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch)) - &
                & (d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
            pdot(nh(i,ind)) = pdot(nh(i,ind))+dumqhi-tbh(i,ind)*dumqbh
            !
            !  dvjch(i)/dc
            !
            dumjch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch)) *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch)) - &
                & (d1ch+ d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
            !
            !  dvj(i)/dnc(ind)
            !
            pdot(nc(ind))=pdot(nc(ind)) +(vjch(i)-vjcb)*(dumjch-dumjcb*tcb(ind)) +(vjcb-vjbh(i))*dumjcb*tcb(ind) - &
                & (vjbh(i)-vjch(i))*dumjch
            !
            !  dvjch(i)/dh(i)
            !
            dumjhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch)) *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch)) - &
                & (d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
            !
            !  dvj(i)/dnh(i,ind)
            !
            pdot(nh(i,ind))=pdot(nh(i,ind)) +(vjch(i)-vjcb)*dumjhi +(vjcb-vjbh(i))*dumjbh*tbh(i,ind) +(vjbh(i)- &
                & vjch(i))* (-dumjbh*tbh(i,ind)-dumjhi)
            !
            !  dv(i)/dh(j)
            !
            do k=1,3
                j = i+k
                if(j.gt.4)j=j-4
                dumqhj=-achdh(j,ind)*(rch(i)-r0ch) *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch)) -(d1ch-d3ch)*exp(-ach* &
                    & (rch(i)-r0ch)))
                dumjhj=-achdh(j,ind)*(rch(i)-r0ch) *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch)) -(d1ch+d3ch)*exp(-ach* &
                    & (rch(i)-r0ch)))*factj
                pdot(nh(j,ind)) = pdot(nh(j,ind))+dumqhj +(vjch(i)-vjcb)*dumjhj -(vjbh(i)-vjch(i))*dumjhj
            enddo
        enddo
    enddo
end subroutine stretch
!
!******************************************************
!
!
!> @brief Evaluate opbend for the C2H7_VBMM potential.
subroutine opbend( vop )
    !
    !  subroutine calculates symmetrized vop potential and derivatives
    !
    implicit none

    double precision :: a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, atheta, b1ch, b1s, b2s, &
        & b3s, bk, bphi, btheta, c1ch, cart, cartnu, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, ddr, deesdr, &
        & degsdr, deijdr, descart, dfdc, dfdelta, dfdh, dgscart, dhdelta, dijcart, ds1, ds2, ds3, dsphi, dstheta, &
        & dtheta0, engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, fkinf, hch3, hdelta, pdot, pengyes, pengygs, &
        & pengyij, q, r, r0cb, r0ch, r0hh, rbh, rcb, rch, s1, s2, s3, sphi, stheta, sum2, sum4, sumd2, sumd4, tbh, tcb, &
        & tch, theta0, vop
    integer :: i, icartr, ind, indexes, irctnt, itemp, j, k, l, mder, msurf, nasurf, natoms, nc, nder, nflag, nh, nhb, &
        & nnb, nnc, nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    dimension sumd2(4),sumd4(4)
    !
    !
    vop = 0.0d0
    !
    !  calculate force constants and their derivatives
    !
    call opforce
    !
    !  calculate out-of-plane xangle and derivatives
    !
    do i=1,4
        j = i+1
        if(j.gt.4)j=j-4
        k = j+1
        if(k.gt.4)k=k-4
        l = k+1
        if(l.gt.4)l=l-4
        !
        !  if i is an even number then we have switched the vectors
        !  from a right handed set to a left handed set
        !
        !  in this case we need to switch vectors k and l around
        !
        if((i.eq.2).or.(i.eq.4))then
            itemp = k
            k     = l
            l     = itemp
        endif
        !
        !  subroutine performs sum over j, k, l
        !  sum2 = sum delta**2
        !  sum4 = sum delta**4
        !
        call calcdelta(i,j,k,l,sum2,sum4)
        sumd2(i) = sum2
        sumd4(i) = sum4
        vop      = vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
    enddo
    do i=1,4
        do j=1,4
            !
            !  overall derivatives of force constants i wrt the bond-length rch(j)
            !
            ddr = dfdelta(i,j)*sumd2(i)+dhdelta(i,j)*sumd4(i)
            !
            !  calculate derivatives in terms of cartesian coordinates:
            !
            do ind=1,3
                pdot(nh(j,ind)) = pdot(nh(j,ind))-tch(j,ind)*ddr/rch(j)
                pdot(nc(ind))   = pdot(nc(ind))+tch(j,ind)*ddr/rch(j)
            enddo
        enddo
    enddo
end subroutine opbend
!
!******************************************************
!
!
!> @brief Evaluate ipbend for the C2H7_VBMM potential.
subroutine ipbend( vip )
    !
    !  subroutine calculates symmetrised in plane bend term
    !  and its derivatives
    !
    implicit none

    double precision :: a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, atheta, b1ch, b1s, b2s, &
        & b3s, bk, bphi, btheta, c1ch, cart, cartnu, costh, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, deesdr, &
        & degsdr, deijdr, descart, dfdc, dfdelta, dfdh, dgscart, dhdelta, dijcart, ds1, ds2, ds3, dsphi, dstheta, dth, &
        & dth0c, dth0k, dthc, dtheta0, dthi, dthj, engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, fkinf, hch3, &
        & hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0cb, r0ch, r0hh, rbh, rcb, rch, s1, s2, s3, sphi, stheta, tbh, &
        & tcb, tch, termth, theta, theta0, vip
    integer :: i, icartr, ind, indexes, irctnt, j, k, mder, msurf, nasurf, natoms, nc, nder, nflag, nh, nhb, nnb, nnc, &
        & nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    dimension costh(4,4),theta(4,4),dth(4,4)
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
    do i=1,3
        do j=i+1,4
            costh(i,j) = tch(i,1)*tch(j,1)+tch(i,2)*tch(j,2) +tch(i,3)*tch(j,3)
            costh(i,j) = costh(i,j)/rch(i)/rch(j)
            theta(i,j) = acos(costh(i,j))
            dth(i,j)   = theta(i,j)-theta0(i,j)
            vip        = vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
            !##
            !                 write(*,*) i,j,f1(j)
            !
            !  calculate partial derivatives wrt cartesian coordinates
            !
            !  calculate pdots wrt theta:
            !
            termth = -1.0d0/sqrt(1.0d0-costh(i,j)*costh(i,j))
            do ind=1,3
                dthi            = -tch(j,ind)/rch(i)/rch(j) +costh(i,j)*tch(i,ind)/rch(i)/rch(i)
                dthi            = dthi*termth
                dthj            = -tch(i,ind)/rch(i)/rch(j) +costh(i,j)*tch(j,ind)/rch(j)/rch(j)
                dthj            = dthj*termth
                dthc            = -(dthi+dthj)
                pdot(nh(i,ind)) = pdot(nh(i,ind)) +fk0(i,j)*f1(i)*f1(j)*dthi*dth(i,j)
                pdot(nh(j,ind)) = pdot(nh(j,ind)) +fk0(i,j)*f1(i)*f1(j)*dthj*dth(i,j)
                pdot(nc(ind))   = pdot(nc(ind)) +fk0(i,j)*f1(i)*f1(j)*dthc*dth(i,j)
                do k=1,4
                    !
                    !  calculate pdots wrt force constants and wrt theta0
                    !
                    dth0k = -dtheta0(i,j,k)*tch(k,ind)/rch(k)
                    dth0c = -dth0k
                    pdot(nh(k,ind))=pdot(nh(k,ind)) -0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k) -0.5d0*tbh(k, &
                        & ind)* dfdh(i,j,k)*dth(i,j)**2/rbh(k) -fk0(i,j)*f1(i)*f1(j)*dth0k*dth(i,j)
                    pdot(nc(ind))=pdot(nc(ind)) +0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k) -fk0(i,j)*f1(i)* &
                        & f1(j)* dth0c*dth(i,j)
                    pdot(nhb(ind)) = pdot(nhb(ind)) +0.5d0*tbh(k,ind)*dfdh(i,j,k)*dth(i,j)**2/rbh(k)
                enddo
            enddo
        enddo
    enddo
end subroutine ipbend
!
!***********************************************************************
!
!> @brief Evaluate calcdelta for the C2H7_VBMM potential.
subroutine calcdelta( i, j, k, l, sum2, sum4 )
    !
    !  subroutine calculates out of plane xangle delta, loops
    !  through delta(i,j), delta(i,k), delta(i,l)
    !
    !   also calculates the derivatives wrt delta
    !
    implicit none

    double precision :: a, a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, argd, atemp1, atemp2, &
        & atemp3, atemp4, atemp5, atheta, axb, b, b1ch, b1s, b2s, b3s, bk, bphi, btheta, c, c1ch, cart, cartnu, cdot, &
        & cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, daxb, deesdr, degsdr, deijdr, deldot, delta, descart, dfdc, &
        & dfdelta, dfdh, dgscart, dhdelta, dijcart, ds1, ds2, ds3, dsphi, dstheta, dtheta0, engyes, engygs, engyij, &
        & ezero, f1, fch3, fdelta, fk0, fkinf, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0cb, r0ch, r0hh, &
        & rbh, rcb, rch, s1, s2, s3, sphi, stheta, sum2, sum4, tbh, tcb, tch, theta0
    integer :: i, icartr, ii, in, ind, indexes, irctnt, j, jj, k, l, mder, msurf, nasurf, natoms, nc, nder, nflag, nh, &
        & nhb, nnb, nnc, nnh, nulbl
    double precision :: norma
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    dimension  delta(4),in(3),a(3),b(3),axb(3),c(4,3),argd(4), &
        & daxb(4,3,3),cdot(4,3,3),atemp2(3)
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
    do ind=1,3
        a(ind) = q(nh(k,ind))-q(nh(j,ind))
        b(ind) = q(nh(l,ind))-q(nh(j,ind))
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
    !
    !  axb is vector a cross b
    !
    axb(1) = a(2)*b(3)-a(3)*b(2)
    axb(2) = a(3)*b(1)-a(1)*b(3)
    axb(3) = a(1)*b(2)-a(2)*b(1)
    norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
    norma  = sqrt(norma)
    !##
    !        write(*,*) 'axb1 : ',i,axb(1)
    !        write(*,*) 'axb2 : ',i,axb(2)
    !        write(*,*) 'axb3 : ',i,axb(3)
    !
    !  c is position vector of h(ii): calculate c(j),c(k),c(l)
    !
    !##
    do ii=1,3
        do ind=1,3
            c(in(ii),ind) = -tch(in(ii),ind)/rch(in(ii))
        enddo
    enddo
    !##
    !
    !  argd is the dot product axb dot c
    !
    do ii=1,3
        argd(in(ii)) = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
        argd(in(ii)) = argd(in(ii))/norma

        delta(in(ii)) = acos(argd(in(ii)))-theta0(i,in(ii))
        !##        write(*,*) 'theta0 ',i,in(ii),theta0(i,in(ii))
        sum2 = sum2+delta(in(ii))**2
        sum4 = sum4+delta(in(ii))**4
    enddo
    !
    !  derivatives of axb wrt hj:
    !
    daxb(j,1,1) = 0.0d0
    daxb(j,1,2) = b(3)-a(3)
    daxb(j,1,3) = -b(2)+a(2)
    daxb(j,2,1) = -b(3)+a(3)
    daxb(j,2,2) = 0.0d0
    daxb(j,2,3) = b(1)-a(1)
    daxb(j,3,1) = b(2)-a(2)
    daxb(j,3,2) = -b(1)+a(1)
    daxb(j,3,3) = 0.0d0
    !
    !  derivatives of axb wrt hk:
    !
    daxb(k,1,1) = 0.0d0
    daxb(k,1,2) = -b(3)
    daxb(k,1,3) = b(2)
    daxb(k,2,1) = b(3)
    daxb(k,2,2) = 0.0d0
    daxb(k,2,3) = -b(1)
    daxb(k,3,1) = -b(2)
    daxb(k,3,2) = b(1)
    daxb(k,3,3) = 0.0d0
    !
    !  derivatives of axb wrt hl:
    !
    daxb(l,1,1) = 0.0d0
    daxb(l,1,2) = a(3)
    daxb(l,1,3) = -a(2)
    daxb(l,2,1) = -a(3)
    daxb(l,2,2) = 0.0d0
    daxb(l,2,3) = a(1)
    daxb(l,3,1) = a(2)
    daxb(l,3,2) = -a(1)
    daxb(l,3,3) = 0.0d0
    !
    !   loop over cdot(in(ii),ind,jind) where we consider deriv of c(in(ii))
    !   wrt h(in(ii),jind) with components jind
    !
    do ii=1,3
        !
        !  deriv of cdot(in(ii),x) wrt x, y, z
        !
        cdot(in(ii),1,1) = 1.0d0/rch(in(ii)) +tch(in(ii),1)*c(in(ii),1)/rch(in(ii))**2
        cdot(in(ii),1,2) = tch(in(ii),2)*c(in(ii),1)/rch(in(ii))**2
        cdot(in(ii),1,3) = tch(in(ii),3)*c(in(ii),1)/rch(in(ii))**2
        !
        !  deriv of cdot(in(ii),y) wrt x, y, z
        !
        cdot(in(ii),2,1) = tch(in(ii),1)*c(in(ii),2)/rch(in(ii))**2
        cdot(in(ii),2,2) = 1.0d0/rch(in(ii)) +tch(in(ii),2)*c(in(ii),2)/rch(in(ii))**2
        cdot(in(ii),2,3) = tch(in(ii),3)*c(in(ii),2)/rch(in(ii))**2
        !
        !  deriv of cdot(in(ii),z) wrt x, y, z
        !
        cdot(in(ii),3,1) = tch(in(ii),1)*c(in(ii),3)/rch(in(ii))**2
        cdot(in(ii),3,2) = tch(in(ii),2)*c(in(ii),3)/rch(in(ii))**2
        cdot(in(ii),3,3) = 1.0d0/rch(in(ii)) +tch(in(ii),3)*c(in(ii),3)/rch(in(ii))**2
    enddo
    !
    do ii=1,3
        do ind=1,3
            deldot = -dtheta0(i,in(ii),i)
            !
            !  derivative wrt h(i,ind)
            !  for  rch(i) only terms are from the derivatives of theta0
            !
            deldot = -deldot*tch(i,ind)/rch(i)
            pdot(nh(i,ind))=pdot(nh(i,ind)) +2.0d0*fdelta(i)*delta(in(ii))*deldot +4.0d0*hdelta(i)*delta(in(ii))**3* &
                & deldot
            !
            !  derivative wrt c(ind)
            !
            deldot        = -deldot
            pdot(nc(ind)) = pdot(nc(ind)) +2.0d0*fdelta(i)*delta(in(ii))*deldot +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
            do jj=1,3
                !
                !  partial derivatives wrt h(in(jj),ind), loop over delta(i,in(ii))
                !
                !   atemp1 is axb dot daxb wrt h(in(jj))
                !
                atemp1 = axb(1)*daxb(in(jj),ind,1) +axb(2)*daxb(in(jj),ind,2) +axb(3)*daxb(in(jj),ind,3)
                atemp1 = atemp1/(norma**3)
                !
                !  atemp2 is deriv of normalised axb
                !
                atemp2(1) = daxb(in(jj),ind,1)/norma-atemp1*axb(1)
                atemp2(2) = daxb(in(jj),ind,2)/norma-atemp1*axb(2)
                atemp2(3) = daxb(in(jj),ind,3)/norma-atemp1*axb(3)
                !
                !  atemp3 is daxb dot c(in(ii))

                atemp3 = atemp2(1)*c(in(ii),1)+atemp2(2)*c(in(ii),2) +atemp2(3)*c(in(ii),3)
                !
                !  atemp4 is axb dot cdot
                !
                atemp4 = 0.0d0
                if(ii.eq.jj)then
                    !
                    !  ie deriv of c(in(ii)) wrt h(in(jj)) is non zero only for ii = jj
                    !
                    atemp4 = axb(1)*cdot(in(ii),1,ind) +axb(2)*cdot(in(ii),2,ind) +axb(3)*cdot(in(ii),3,ind)
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
                pdot(nh(in(jj),ind))=pdot(nh(in(jj),ind)) +2.0d0*fdelta(i)*delta(in(ii))*deldot +4.0d0*hdelta(i)* &
                    & delta(in(ii))**3*deldot
                !
                !  for carbon the only contributions are from axb dot cdot term and
                !  from theta0 and derivative cdot wrt carbon=-cdot wrt hydrogen
                !
                deldot = 1.0d0/sqrt(1.0d0-argd(in(ii))**2)*atemp4
                deldot = deldot-atemp5
                pdot(nc(ind))=pdot(nc(ind)) +2.0d0*fdelta(i)*delta(in(ii))*deldot +4.0d0*hdelta(i)*delta(in(ii))**3* &
                    & deldot
            enddo
        enddo
    enddo
end subroutine calcdelta

!******************************************************
!
!> @brief Evaluate opforce for the C2H7_VBMM potential.
subroutine opforce
    !
    !  calculates the out-of-plane bending force constants
    !  and their derivatives
    !
    implicit none

    double precision :: a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, atheta, b1ch, b1s, b2s, &
        & b3s, bk, bphi, btheta, c1ch, cart, cartnu, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, deesdr, degsdr, &
        & deijdr, descart, dfdc, dfdelta, dfdh, dgscart, dhdelta, dijcart, ds1, ds2, ds3, dsphi, dstheta, dswitch, &
        & dtheta0, engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, fkinf, hch3, hdelta, pdot, pengyes, pengygs, &
        & pengyij, q, r, r0cb, r0ch, r0hh, rbh, rcb, rch, s1, s2, s3, sphi, stheta, switch, tbh, tcb, tch, theta0
    integer :: i, icartr, indexes, irctnt, j, mder, msurf, nasurf, &
        & natoms, nc, nder, nflag, nh, nhb, nnb, nnc, nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    dimension switch(4),dswitch(4,4)
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
    do i=1,4
        fdelta(i) = switch(i)*fch3
        hdelta(i) = switch(i)*hch3
        do j=1,4
            dfdelta(i,j) = dswitch(i,j)*fch3
            dhdelta(i,j) = dswitch(i,j)*hch3
        enddo
    enddo
end subroutine opforce
!
!******************************************************
!
!
!> @brief Evaluate ipforce for the C2H7_VBMM potential.
subroutine ipforce
    !
    !  calculates the symmetrised in plane bend force constants and
    !  all partial derivatives involving them
    !
    implicit none

    double precision :: a1, a1ch, a1s, a2, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, arga1, arga2, &
        & atheta, b1ch, b1s, b2s, b3s, bk, bphi, btheta, c1ch, cart, cartnu, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, &
        & d3hh, deesdr, degsdr, deijdr, descart, df1dc, df1dh, dfdc, dfdelta, dfdh, dfk0, dgscart, dhdelta, dijcart, ds1, &
        & ds2, ds3, dsphi, dstheta, dtheta0, duma1, duma2, engyes, engygs, engyij, ezero, f0, f1, f2, fch3, fdelta, fk0, &
        & fkinf, hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0cb, r0ch, r0hh, rbh, rcb, rch, s1, s2, s3, sphi, &
        & stheta, tbh, tcb, tch, theta0
    integer :: i, icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nc, nder, nflag, nh, nhb, nnb, nnc, nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    dimension dfk0(4,4,4),df1dc(4),df1dh(4)
    !
    !  set force constant at asymptotes
    !
    f0 = fkinf+ak
    f2 = fkinf
    !
    !##
    fk0(1,2) = f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)*s2(4)-1.0d0)
    fk0(1,3) = f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)*s2(4)-1.0d0)
    fk0(1,4) = f0+f0*(s1(1)*s1(4)-1.0d0)+(f0-f2)*(s2(2)*s2(3)-1.0d0)
    fk0(2,3) = f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)*s2(4)-1.0d0)
    fk0(2,4) = f0+f0*(s1(2)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(3)-1.0d0)
    fk0(3,4) = f0+f0*(s1(3)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(2)-1.0d0)

    !##
    !           write(*,*)
    !           write(*,*) 'f0    : ',f0
    !           write(*,*) 'f2    : ',f2
    !           write(*,*) 's1(1) : ',s1(1)
    !           write(*,*) 's1(2) : ',s1(2)
    !           write(*,*) 's2(3) : ',s2(3)
    !           write(*,*) 's2(4) : ',s2(4)
    !           write(*,*)

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
    do i=1,4
        !
        !  calculate the terms f1(i)
        !
        !##
        arga1 = aa1*rbh(i)*rbh(i)
        arga2 = aa4*(rbh(i)-r0hh)*(rbh(i)-r0hh)
        a1    = 1.0d0-exp(-arga1)
        a2    = aa2+aa3*exp(-arga2)
        f1(i) = a1*exp(-a2*(rch(i)-r0ch)**2)

        !           write(*,*) 'A1 :',a1
        !           write(*,*) 'A2 :',a2
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
!
!******************************************************
!
!
!> @brief Evaluate switchf for the C2H7_VBMM potential.
subroutine switchf
    !
    !  calculates switching functions: s3,sphi,stheta
    !  and their derivatives ds3,dsphi,dstheta
    !
    implicit none

    double precision :: a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, argmax, args1, args2, &
        & args3, argsphi, argstheta, atheta, b1ch, b1s, b2s, b3s, bk, bphi, btheta, c1ch, cart, cartnu, cphi, ctheta, &
        & d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, deesdr, degsdr, deijdr, descart, dfdc, dfdelta, dfdh, dgscart, dhdelta, &
        & dijcart, ds1, ds2, ds3, dsphi, dstheta, dtheta0, engyes, engygs, engyij, ezero, f1, fch3, fdelta, fk0, fkinf, &
        & hch3, hdelta, pdot, pengyes, pengygs, pengyij, q, r, r0cb, r0ch, r0hh, rbh, rcb, rch, s1, s2, s3, sphi, stheta, &
        & tbh, tcb, tch, theta0
    integer :: i, icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nc, nder, nflag, nh, nhb, nnb, nnc, nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt1cm/ r(n3atom), engygs, degsdr(n3atom)
    common/pt3cm/ ezero(isurf+1)
    common/pt4cm/ engyes(isurf), deesdr(n3atom,isurf)
    common/pt5cm/ engyij(jsurf), deijdr(n3atom,jsurf)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /xangles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    !##
    a1s = 1.5313681d-7
    b1s = -4.6696246d0
    a2s = 1.0147402d-7
    b2s = -12.363798d0
    !
    !  use double precision criterion:
    !
    !  tanh(19.0d0)=1.0d0
    !
    argmax = 19.0d0
    do i=1,4
        args1 = a1s*(rch(i)-r0ch)*(rch(i)-b1s)**8
        if(args1.lt.argmax)then
            s1(i)  = 1.0d0-tanh(args1)
            ds1(i) = a1s*((rch(i)-b1s)**8 +8.0d0*(rch(i)-r0ch)*(rch(i)-b1s)**7)
            ds1(i) = -ds1(i)/cosh(args1)**2
        else
            s1(i)  = 0.0d0
            ds1(i) = 0.0d0
        endif
        !
        args2 = a2s*(rch(i)-r0ch)*(rch(i)-b2s)**6
        if(args2.lt.argmax)then
            s2(i)  = 1.0d0-tanh(args2)
            ds2(i) = a2s*((rch(i)-b2s)**6 +6.0d0*(rch(i)-r0ch)*(rch(i)-b2s)**5)
            ds2(i) = -ds2(i)/cosh(args2)**2
        else
            s2(i)  = 0.0d0
            ds2(i) = 0.0d0
        endif
        !##          write(*,*) i,s1(i),rch(i)

        !
        !  calculate s3 and ds3
        !
        args3 = a3s*(rch(i)-r0ch)*(rch(i)-b3s)**2
        if (args3.lt.argmax)then
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
        if(rch(i).lt.3.8d0)then
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
        if(rch(i).lt.3.8d0)then
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
!
!> @brief Evaluate PREPOT for the C2H7_VBMM potential.
subroutine prepot
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

    double precision :: a1ch, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, atheta, b1ch, b3s, bk, bphi, btheta, &
        & c1ch, cart, cartnu, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, ezero, fact1, fact2, fch3, fkinf, hch3, &
        & r0cb, r0ch, r0hh
    integer :: i, icartr, icount, ind, indexes, irctnt, mder, msurf, nasurf, natoms, nc, nder, nflag, nh, nhb, nnb, nnc, &
        & nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi     = 3.141592653589793d0
    integer, parameter          :: natom  = 25
    integer, parameter          :: n3tmmn = 18
    !
    common/pt3cm/ ezero(isurf+1)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    common /ndx/     nc(3),nhb(3),nh(4,3)
    !
    !
    !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
    !
    if(natoms.gt.25) then
        write(nflag(18),'(2X,''STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS'')')
        stop
    end if
    !
    do i=1,5
        ref(i) = ' '
    end do
    !
    ref(1) = 'M. J. T. Jordan and R. G. Gilbert'
    ref(2) = 'J. Chem. Phys. Vol. 102, p. 5669, 1995'
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
    !  calculate indexes for coordinates
    !
    do ind=1,3
        icount   = ind-3
        nc(ind)  = 3*nnc+icount
        nhb(ind) = 3*nnb+icount
        do i=1,4
            nh(i,ind) = 3*nnh(i)+icount
        enddo
    enddo
    !
    !  convert to appropriate units:
    !
    !  energy   in 1.0d+05 j/mol
    !  time     in 1.0d-14 s
    !  xdistance in 1.0d-10 m
    !  xangles   in radians
    !  mass     in amu
    !
    fact1 = 0.041840d0
    fact2 = 6.022045d0
    !
    !##----------------------------------
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
    !##----------------------------------
    !
end subroutine prepot
!> @brief Evaluate CSETUP for the C2H7_VBMM potential.
subroutine csetup( n3tm )
    !       dummy subroutine to make the code compatable with chk.f

    implicit none

    integer :: n3tm
end subroutine csetup
!
block data ptpacm
    !
    implicit none

    double precision :: a1ch, a1s, a2s, a3s, aa1, aa2, aa3, aa4, acb, ahh, ak, anuzero, aphi, atheta, b1ch, b1s, b2s, &
        & b3s, bk, bphi, btheta, c1ch, cart, cartnu, cphi, ctheta, d1cb, d1ch, d1hh, d3cb, d3ch, d3hh, dfdc, dfdelta, &
        & dfdh, dhdelta, ds1, ds2, ds3, dsphi, dstheta, dtheta0, ezero, f1, fch3, fdelta, fk0, fkinf, hch3, hdelta, pdot, &
        & q, r0cb, r0ch, r0hh, rbh, rcb, rch, s1, s2, s3, sphi, stheta, tbh, tcb, tch, theta0
    integer :: icartr, indexes, irctnt, mder, msurf, nasurf, natoms, &
        & nc, nder, nflag, nh, nhb, nnb, nnc, nnh, nulbl
    !
    character(len=75) :: ref(5)
    !
    integer, parameter :: n3atom = 75
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    !
    double precision, parameter :: pi    = 3.141592653589793d0
    integer, parameter          :: natom = 25
    !
    common/pt3cm/ ezero(isurf+1)
    !
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    !
    !
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    !
    common /potcm/ nnc,nnb,nnh(4), &
        & r0ch,d1ch,d3ch, &
        & a1ch,b1ch,c1ch, &
        & r0hh,d1hh,d3hh,ahh, &
        & r0cb,d1cb,d3cb,acb, &
        & a3s,b3s,aphi,bphi,cphi, &
        & atheta,btheta,ctheta, &
        & fch3,hch3, &
        & fkinf,ak,bk,aa1,aa2,aa3,aa4
    !
    common /angles/  theta0(4,4),dtheta0(4,4,4)
    common /bonds/   rcb,rch(4),rbh(4)
    common /coords/  tcb(3),tch(4,3),tbh(4,3)
    common /delta1/  fdelta(4),hdelta(4)
    common /delta2/  dfdelta(4,4),dhdelta(4,4)
    common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
    common /fsw1/    a1s,b1s,a2s,b2s
    common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
    common /ndx/     nc(3),nhb(3),nh(4,3)
    common /op1/     s3(4),ds3(4)
    common /qpdot_pl/   q(150),pdot(150)
    common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
    !
    data nasurf /1,35*0/
    data nder /1/
    data nflag /1,1,15*0,6,0,0/
    !
    data anuzero /0.0d0/
    data icartr,msurf,mder/1,0,1/
    data nulbl /25*0/
    data natoms /6/
    data nnc     /2/
    data nnb     /6/
    data nnh     /3,4,5,1/
    data r0ch    /   1.09397d0/
    data d1ch    / 112.230d0/
    data d3ch    /  38.834d0/
    data a1ch    /   1.71300d0/
    data b1ch    /   0.13500d0/
    data c1ch    /   6.61404d0/
    data r0hh    /   0.74191d0/
    data d1hh    / 109.458d0/
    data d3hh    /  39.664d0/
    data ahh     /   1.9457d0/
    data r0cb    /   1.09397d0/
    data d1cb    /  26.409d0/
    data d3cb    /  20.063d0/
    data acb     /   1.8530285d0/
    data a3s     /   0.1419147d0/
    data b3s     /  -0.3068450d0/
    data aphi    /   0.5287903d0/
    data bphi    /   0.4006638d0/
    data cphi    /   1.9209937d0/
    data atheta  /   0.9078714d0/
    data btheta  /   0.3548859d0/
    data ctheta  /   1.8915497d0/
    data fch3    /   0.0957500d0/
    data hch3    /   0.1915000d0/
    data fkinf   /   0.4077000d0/
    data ak      /   0.1260000d0/
    data bk      /  10.7132d0/
    data aa1     /   3.213952d0/
    data aa2     /   1.599963d0/
    data aa3     /   2.165953d0/
    data aa4     /  11.569977d0/
    !
end

!23456789
!> @brief Evaluate xgaussfit for the C2H7_VBMM potential.
subroutine xgaussfit( r, fit )
    implicit none

    double precision :: r
    double precision :: fit

    integer, parameter :: ng = 6

    integer          :: i,j,k
    double precision :: d
    double precision :: ans
    double precision :: r0(ng)
    double precision :: alpha(ng)
    double precision :: cf(ng)
    double precision :: fx

    r0(1) = 1.949230d0
    r0(2) = 1.013860d0
    r0(3) = 2.790270d0
    r0(4) = (r0(1)+r0(2))*0.50d0
    r0(5) = (r0(1)+r0(3))*0.50d0
    r0(6) = (r0(2)+r0(3))*0.50d0
    do i=1,ng
        alpha(i) = 0.250d0
    end do

    cf(1) = -1178.890067483491520d0
    cf(2) = 32.80653355230102620d0
    cf(3) = 10.31188564665787060d0
    cf(4) = -159.8465369524321600d0
    cf(5) = 17.33838608018073660d0
    cf(6) = 1280.367802966489080d0

    ans = 0.0d0
    do i=1,ng
        d   = r-r0(i)
        fx  = dexp(-alpha(i)*d*d)
        ans = ans+(cf(i)*fx)
    end do

    fit = ans

end subroutine xgaussfit
