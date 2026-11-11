!   System:                        NH3
!   Functional form:               polynomial 2x2 diabatic fit
!   Common name:                   NH3 (code 2, optimized)
!   Number of derivatives:         1
!   Number of bodies:              4
!   Number of electronic surfaces: 2
!   Interface:                     Section-2
!  Notes: This set of coupled potentials is from Ref. 2.  The earlier
!  potentials in  Ref. 1 are not included in POTLIB because they are
!  not recommended.  The potential routines give the ground (U11) and
!  first singlet excited (U22) diabats, the diabatic coupling (U12),
!  the ground (V1) and first singlet excited (V2) adiabats, and the
!  Cartesian gradients of U11, U22, U12, V1, and V2.
!  This potential has additional correction terms used for fitting
!  to the experimental dissociation limits.
!  This is code 2 for the ammonia potentials and diabatic coupling.
!  Code 2 differs from code 1 in the order in which the terms are
!  evaluated and so these two codes give slightly different answers
!  due to roundoff; both are correct within the limits of round off
!  error.  Code 2 is more efficient than code 1 (it runs faster).
!  For the convenience of users who may wish to modify the potential,
!  we include both versions in POTLIB.  We recommend using code 2 for
!  applications because it runs faster.
! References:
!      1) Nangia, S.; Truhlar, D.G. J. Chem. Phys. 2006, 124, 124309
!      2) Li, Z.H.; Valero, R., Truhlar, D.G. Theor. Chem. Acc. 2007,
!         118, 9-24.
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  All gradients are in hartree/bohr
!  Format for   Xcart
!  Xcart(1:3)   is N
!  Xcart(4:6)   is H1
!  Xcart(7:9)   is H2
!  Xcart(10:12) is H3

!> @brief Implement NH3_GEN_FFW2_2007, a diabatic potential matrix for NH3.
module chempotpy_nh3_gen_ffw2_2007_dpem_pes
    implicit none
    private
    public :: evaluate
    integer, parameter          :: coefficient_variable_count = 7
    double precision, parameter :: nhcoeu1(62)                = [-.1698267005d+01, -.3116975297d+01, .4989324584d+02, .6140536234d+02, &
        & -.9137903351d+02, -.4170308921d+02, -.2822911145d+03, -.1324724944d+03, &
        & .4554289660d+02, .6482512608d+02, .1622981341d+03, -.1418440334d+03, &
        & .3128761841d+03, .2291158383d+03, .4212981689d+03, .5904053030d+03, &
        & -.2011576383d+03, -.2847074833d+03, -.4207881721d+03, -.3401129381d+03, &
        & -.1017600665d+03, -.3135766858d+02, .7766878674d+02, -.1863622757d+03, &
        & -.2454112660d+03, -.2311744548d+03, -.7810834823d+02, .1476582353d+03, &
        & .2809941380d+03, .3033038990d+03, .6163060145d+02, .1175366540d+01, &
        & .3888738002d+01, -.3651170315d+02, -.3522094635d+02, .8722843200d+02, &
        & .7398300089d+02, .1872947617d+03, .1217028644d+03, -.3865307026d+02, &
        & -.9112045320d+02, -.8906173604d+02, .5130386478d+02, -.2935045013d+03, &
        & -.3045399445d+03, -.2456762591d+03, -.4827738566d+03, .1665284893d+03, &
        & .3478380563d+03, .1826085704d+03, .3027131358d+03, .1219769642d+03, &
        & .4171232808d+02, -.9647425161d+01, .2035053367d+03, .2297924112d+03, &
        & .2731855209d+03, .3746015525d+02, -.1574954942d+03, -.2708242232d+03, &
        & -.3486651786d+03, -.1927783168d+02]
    double precision, parameter :: nhcoeu2(71) = [-.5839297267d+01, .1948308393d+03, .1476399426d+03, -.1952770539d+03, &
        & -.4956220164d+03, .1455034861d+04, -.6852363440d+03, -.7980707978d+03, &
        & -.7150307241d+02, -.1898884973d+04, .6944385319d+03, .1278073813d+04, &
        & .2496161582d+03, -.6605419506d+03, -.2259061202d+03, .3208341955d+03, &
        & .1438880740d+04, -.2462117089d+04, .1264999815d+03, .1618808427d+04, &
        & -.4277636246d+03, .2889205530d+04, -.4525182746d+03, -.2216782595d+04, &
        & -.5474397436d+03, .9523771905d+03, .8001804807d+03, -.2413673183d+03, &
        & -.4689796908d+03, -.1226082702d+03, .2013096370d+03, .2497300043d+03, &
        & .2354985426d+03, -.1609355796d+03, -.5933048064d+03, -.2301959837d+03, &
        & -.1750594296d+03, -.1343884916d+03, .1997871612d+03, .4919735429d+03, &
        & -.1486208797d+04, .6482490789d+03, .7895412677d+03, .3391285547d+02, &
        & .1900289391d+04, -.6724668674d+03, -.1274908946d+04, -.2563702497d+03, &
        & .6673491569d+03, .2226064703d+03, -.3178753276d+03, -.1352904130d+04, &
        & .2450815486d+04, -.1269876680d+03, -.1626526640d+04, .4191476857d+03, &
        & -.2837798016d+04, .4598348283d+03, .2208871680d+04, .5496794539d+03, &
        & -.9765928983d+03, -.7638623676d+03, .2378197831d+03, .4359353643d+03, &
        & .9157000100d+02, -.2028118062d+03, -.2304177594d+03, -.2325038433d+03, &
        & .1797589106d+03, .5547260471d+03, .2316886253d+03]
    double precision, parameter :: nhcoeu12(45) = [-.1370792834d+01, .6836772012d+00, .6331956838d+00, .1048352694d+01, &
        & .1450863054d+01, -.6761575220d+00, -.1763525455d+01, -.2757812968d+00, &
        & -.1771035818d+00, .1521156999d+00, .3338335016d+00, -.7011767589d+00, &
        & .1057293370d+00, -.1499798375d+00, .8059440388d+00, .1100511326d+01, &
        & -.9054527734d+00, -.1502093390d+00, -.5755857753d+00, -.1859473317d+01, &
        & .1376357981d+01, .8781480771d+00, -.3398810498d+00, .1058482524d+01, &
        & -.8694281520d+00, .1168933914d+00, .9741155480d+00, -.5185373065d+00, &
        & .4463417737d+00, -.8328439796d+00, -.4614413328d+00, .3513046496d+00, &
        & -.6307182534d-01, .6508886978d+00, .3787965679d+00, -.6678165549d+00, &
        & .1524235094d+00, .1114227008d+00, -.3191010676d+00, .5041174980d+00, &
        & -.2955290687d+00, -.7789245841d+00, .5060705599d+00, -.2054781302d+00, &
        & .1836178250d+00]

    contains

    !> @brief Evaluate the NH3_GEN_FFW2_2007 diabatic potential matrix for NH3.
    subroutine evaluate( x, igrad, u, ug )
        use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 2
        integer, parameter            :: natoms  = 4
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: u(nstates,nstates)
        double precision, intent(out) :: ug(nstates,nstates,natoms,3)
        double precision              :: v1, v2, u11, u12, u22
        double precision              :: tx(12), gu11(12), gu12(12), gu22(12)
        double precision              :: gv1(12), gv2(12)
        double precision              :: dudx(nstates,nstates,12)
        double precision              :: t(nstates,nstates)
        double precision              :: tmpmat(nstates,nstates)
        double precision              :: hx(nstates,nstates,12), gx(nstates,12)
        integer                       :: iatom, idir, i, j, k, l
        !initialize
        u = 0.d0; ug = 0.d0
        do iatom=1, natoms
            do idir=1, 3
                j     = (iatom-1)*3+idir
                tx(j) = x(iatom,idir)/legacy_bohr_to_angstrom
            enddo
        enddo
        call pot(tx, u11, u22, u12, v1, v2, gu11, gu22, gu12, gv1, gv2)
        u(1,1)      = u11*legacy_hartree_to_electronvolt
        u(1,2)      = u12*legacy_hartree_to_electronvolt
        u(2,1)      = u12*legacy_hartree_to_electronvolt
        u(2,2)      = u22*legacy_hartree_to_electronvolt
        dudx(1,1,:) = gu11(:)*51.422067
        dudx(2,2,:) = gu22(:)*51.422067
        dudx(1,2,:) = gu12(:)*51.422067
        dudx(2,1,:) = gu12(:)*51.422067
        do iatom=1, natoms
            do idir=1, 3
                j                  = (iatom-1)*3+idir
                ug(:,:,iatom,idir) = dudx(:,:,j)
            enddo
        enddo
    end subroutine evaluate

    !> @brief Evaluate diagonalize for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine diagonalize( n, a_ss, ev_s, u_ss )
        implicit none
        integer, intent(in)           :: n
        double precision, intent(in)  :: a_ss(n,n)
        double precision, intent(out) :: u_ss(n,n)
        double precision, intent(out) :: ev_s(n)
        integer                       :: io, i
        double precision, allocatable :: work_d(:)
        integer                       :: lwork_d
        u_ss    = a_ss
        lwork_d = -1
        allocate( work_d(1) )
        call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
        lwork_d = int(work_d(1))
        deallocate( work_d )
        allocate( work_d(lwork_d) )
        call dsyev('V', 'L', n, u_ss, n, ev_s, work_d, lwork_d, io)
    end subroutine diagonalize

    !> @brief Evaluate pot for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine pot( xcart, u11, u22, u12, v1, v2, gu11, gu22, gu12, gv1, gv2 )
        implicit none
        ! I/O variables
        integer, parameter :: natom = 4
        double precision   :: xcart(12)
        double precision   :: u11
        double precision   :: u22
        double precision   :: u12
        double precision   :: v1
        double precision   :: v2
        double precision   :: gu11(12)
        double precision   :: gu22(12)
        double precision   :: gu12(12)
        double precision   :: gv1(12)
        double precision   :: gv2(12)
        ! Local variables
        ! cart(12) local Xcart in ANG, x1,y1,z1,...,x4,y4,z4
        double precision :: cart(12)
        integer          :: i, j
        ! Unit conversion parameters
        double precision, parameter :: autoang       = 0.5291772108d0
        double precision, parameter :: autoev        = 27.2113845d0
        double precision, parameter :: ev_to_hartree = 1.d0/autoev
        double precision, parameter :: ang_to_bohr   = 1.0d0/autoang
        double precision, parameter :: gconv         = ev_to_hartree/ang_to_bohr
        ! convert geometry from bohr to ang
        do i=1, 12
            cart(i) = xcart(i)*autoang
        enddo
        ! get energies and gradients and convert to atomic unit
        call packpot(cart, u11, u22, u12, v1, v2, gu11, gu22, gu12, gv1, gv2)
        u11 = u11 * ev_to_hartree
        u22 = u22 * ev_to_hartree
        u12 = u12 * ev_to_hartree
        v1  = v1  * ev_to_hartree
        v2  = v2  * ev_to_hartree
        do i=1, 12
            gu11(i) = gu11(i) * gconv
            gu22(i) = gu22(i) * gconv
            gu12(i) = gu12(i) * gconv
            gv1(i)  = gv1(i)  * gconv
            gv2(i)  = gv2(i)  * gconv
        enddo
    end subroutine pot

    !23456789
    !> @brief Evaluate packpot for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine packpot( coord, u11, u22, u12, v1, v2, &
            & gcu11, gcu22, gcu12, gcv1, gcv2 )
        implicit none
        ! I/O variables
        ! coord: x1,y1,z1...x4,y4,z4
        double precision :: coord(12)
        double precision :: u11, u22, u12, v1, v2
        double precision :: gcu11(12)
        double precision :: gcu22(12)
        double precision :: gcu12(12)
        double precision :: gcv1(12)
        double precision :: gcv2(12)
        ! Common blocks for coefficients and indices
        integer, parameter :: nv    = 7
        integer, parameter :: ncu1  = 62
        integer, parameter :: ncu2  = 71
        integer, parameter :: ncu12 = 45
        integer            :: nhindu1(ncu1,nv)
        integer            :: nhindu2(ncu2,nv)
        integer            :: nhindu12(ncu12,nv)
        ! Distances 1-3: NH, 4-6: HH
        double precision :: rr(6)
        ! gradients of distances to cartesian coordinates
        double precision :: grr(12,6)
        double precision :: derep(6)
        ! gradients of integernal coordinates
        double precision :: gu11(12)
        double precision :: gu22(12)
        double precision :: gu12(12)
        double precision :: gv1(12)
        double precision :: gv2(12)
        ! Parameter List
        double precision, parameter :: pi = 3.1415926535897930d0
        ! Internal coordinates used
        ! internal coodinates: three bond distances, three bond angles, one
        ! trisector angle
        double precision :: angint(7)
        ! gradients of angle internal coordinates to cartesian coordinates
        double precision :: gangint(12,7)
        ! Bowman internal coodinates: three bond distances, three bond angles
        ! (the three angles in the projected plane), one trisector angle
        double precision :: bowint(7)
        ! gradients of bowman internal coordinates to cartesian coordinates
        double precision :: gbowint(12,7)
        ! Misc
        integer          :: i, ii, j, k, n
        double precision :: erep
        ! to save the first three bond distances. If they are calculated in any
        ! subroutines, avoid to calculate them again
        ! This code computes the fitted pot using the polynomial funcs as
        ! described by Bowman for NH3
        ! INPUT FILES : TWO INPUT FILES ARE NEEDED
        !       [1]. Coefficient File
        !       [2]. Index File (contains the powers for the polynomial)
        ! Initialize energy derivatives
        call calcindu1(nhindu1, 2, 2, 2, 2, 2, 2, 2, ncu1)
        call calcindu2(nhindu2, 3, 3, 3, 2, 2, 2, 2, ncu2)
        call calcindu12(nhindu12, 2, 2, 2, 2, 2, 2, 5, ncu12)
        do i=1, 12
            gcu11(i) = 0.d0
            gcu22(i) = 0.d0
            gcu12(i) = 0.d0
            gcv1(i)  = 0.d0
            gcv2(i)  = 0.d0
        enddo
        ! Initialize gangint,and dbow
        do i=1, 7
            do j=1, 12
                gangint(j,i) = 0.d0
                gbowint(j,i) = 0.d0
            enddo
        enddo
        ! Calculate internal coodinates
        call tobowman(coord, bowint, gbowint)
        call carttoang(coord, angint, gangint)
        ! Calculate the distances
        call rarray(4, coord, 12, rr, grr, 6)
        ! Calculate bowman internal coodinates
        ! The 7th internal for angint and bowint are both beta
        angint(7) = bowint(7)
        ! Derivatives of 7 (beta) are the same for gangint and gbowint, which
        ! are calculated in toBowman
        do j=1, 12
            gangint(j,7) = gbowint(j,7)
        enddo
        !     Calculate U11,U22,wighout repulsive
        call sumpoly(ncu1, nhindu1, angint, nhcoeu1, u11, gu11, 'U11')
        call sumpoly(ncu2, nhindu2, angint, nhcoeu2, u22, gu22, 'U22')
        call sumpoly(ncu12, nhindu12, angint, nhcoeu12, u12, gu12, 'U12')
        !     Convert derivatives of internal coordinates to Catersian by chain
        !     rule U12 used bowint while other used angint, thus must convert to
        !     the same cartesian derivs before calling getAdiabat
        do j=1, 12
            do i=1, 7
                gcu11(j) = gcu11(j) + gu11(i) * gangint(j,i)
                gcu22(j) = gcu22(j) + gu22(i) * gangint(j,i)
                gcu12(j) = gcu12(j) + gu12(i) * gangint(j,i)
            enddo
        enddo
        call repulsive(rr, erep, derep, 'U11')
        u11 = u11 + erep
        do j=1, 12
            do i=1, 6
                gcu11(j) = gcu11(j) + derep(i) * grr(j,i)
            enddo
        enddo
        call repulsive(rr, erep, derep, 'U22')
        u22 = u22 + erep
        do j=1, 12
            do i=1, 6
                gcu22(j) = gcu22(j) + derep(i) * grr(j,i)
            enddo
        enddo
        !     Calculate V1,V2 from U11 and U22
        call getadiabat(12, u11, u22, u12, v1, v2, &
            & gcu11, gcu22, gcu12, gcv1, gcv2)
    end subroutine packpot

    !> @brief Evaluate sumpoly for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine sumpoly( nc, ind0, bowint, coeff, v, gv, cpot )
        !    This part are all the same for U,U1,U12 except the poly normial
        !    called
        implicit none
        integer          :: nc
        integer          :: ind0(nc,7)
        double precision :: bowint(7)
        double precision :: x(7), xin(7), dx(3), dxin(3)
        double precision :: g(3), gin(3)
        double precision :: g1(3), g1in(3)
        double precision :: x7(0:7)
        double precision :: dx7(0:7)
        double precision :: coeff(nc)
        double precision :: v
        double precision :: gv(7)
        character(len=3) :: cpot
        integer          :: i, j, k, p
        integer          :: i1, i2, i3, ip(3)
        double precision :: s(3,3)
        double precision :: fval
        double precision :: fsum
        double precision :: y(7)
        double precision :: geoint(7)
        !     intemediate storage of derivatives
        ! dy(i,j) = dy(i)/dx(j)
        double precision :: dgeoint(7,7)
        double precision :: dy(7,7)
        double precision :: dfval(7), dfsum(7)
        integer          :: indp(7)
        double precision :: cos1, cos2, cos3
        double precision :: dbeta
        ! Parameters
        double precision, parameter :: pi     = 3.1415926535897930d0
        double precision, parameter :: deg2au = pi/180.d0
        double precision, parameter :: hba    = 0.5d0*pi
        ! U11
        double precision :: frb, frc, alp, cosgta
        ! U22
        double precisionre, frd, dk
        ! U12
        double precision :: fru12b, fru12c
        if( cpot .ne. 'U12' )then
            cos1  = dcos(  bowint(4)*deg2au )
            cos2  = dcos(  bowint(5)*deg2au )
            cos3  = dcos(  bowint(6)*deg2au )
            dx(1) = dsin(  bowint(4)*deg2au )*deg2au
            dx(2) = dsin(  bowint(5)*deg2au )*deg2au
            dx(3) = dsin(  bowint(6)*deg2au )*deg2au
            x(7)  = bowint(7)
        endif
        dbeta = bowint(7)*deg2au-hba
        if( cpot .eq. 'U11' )then
            frb    = .1872506360d+01
            frc    = .4278236310d+00
            alp    = .1009053173d+00
            cosgta = -0.2870261592d0
            x(1)   = dexp(-frb*(bowint(1) - frc))
            x(2)   = dexp(-frb*(bowint(2) - frc))
            x(3)   = dexp(-frb*(bowint(3) - frc))
            x(4)   = cosgta - cos1
            x(5)   = cosgta - cos2
            x(6)   = cosgta - cos3
            g(1)   = dexp( -alp * bowint(1)*bowint(1) )
            g(2)   = dexp( -alp * bowint(2)*bowint(2) )
            g(3)   = dexp( -alp * bowint(3)*bowint(3) )
            do i=0, 2
                x7(i)  = dcos( dble(i) * dbeta )
                dx7(i) = -dble(i)*dsin( dble(i) * dbeta )*deg2au
            enddo
        elseif( cpot .eq. 'U22' )then
            frb    = .1633658498d+01
            frc    = .1987056751d+00
            alp    = .3037603850d+00
            re     = .7147714247d+00
            frd    = .4546981630d+01
            dk     = -.2617146615d+00
            cosgta = -0.5d0
            x(1)   = dexp(-frb*(bowint(1) - frc))
            x(2)   = dexp(-frb*(bowint(2) - frc))
            x(3)   = dexp(-frb*(bowint(3) - frc))
            x(4)   = cosgta - cos1
            x(5)   = cosgta - cos2
            x(6)   = cosgta - cos3
            g(1)   = dexp( -alp * bowint(1)*bowint(1) )
            g(2)   = dexp( -alp * bowint(2)*bowint(2) )
            g(3)   = dexp( -alp * bowint(3)*bowint(3) )
            g1(1)  = dexp(-frd*(bowint(1)-re)**2)
            g1(2)  = dexp(-frd*(bowint(2)-re)**2)
            g1(3)  = dexp(-frd*(bowint(3)-re)**2)
            do i=0, 2
                x7(i) = dcos( dble(i) * dbeta ) + 0.5d0*dk*dbeta*dbeta/(0.05d0+dbeta*dbeta)
                dx7(i) = -dble(i)*dsin( dble(i) * dbeta ) * deg2au + dk*dbeta*deg2au/(0.05d0+dbeta*dbeta) - dk*dbeta** &
                    & 3* deg2au/(0.05d0+dbeta*dbeta)**2
            enddo
        elseif( cpot .eq. 'U12' )then
            fru12c = .1501070926d+01
            fru12b = 1.00d0
            g(1)   = dtanh(fru12b*(bowint(1)-fru12c))
            g(2)   = dtanh(fru12b*(bowint(2)-fru12c))
            g(3)   = dtanh(fru12b*(bowint(3)-fru12c))
            do i=3, 5
                x7(i)  = -dsin( dble(i)*dbeta )
                dx7(i) = -dcos( dble(i)*dbeta )*dble(i)*deg2au
            enddo
        endif
        !     Initialize V,gV
        v = 0.0d0
        do i=1, 7
            gv(i) = 0.d0
        enddo
        do i1=1, 3
            do i2=1, 3
                do i3=1, 3
                    if( i2 == i1 ) cycle
                    if( i3 == i1 ) cycle
                    if( i3 == i2 ) cycle
                    !       symmetrize x to get y
                    ip(1) = i1
                    ip(2) = i2
                    ip(3) = i3
                    call  symopt(3, ip, s)
                    call  getsym7(s, x, xin)
                    call  getsym3(s, dx, dxin)
                    call  getsym3(s, g, gin)
                    call  getsym7g(s, bowint, geoint, dgeoint)
                    !       symmetrize x to get y
                    if( cpot .eq. 'U22' )then
                        call getsym3(s, g1, g1in)
                    endif
                    do k=1, nc
                        do i=1, 7
                            indp(i) = ind0(k,i)
                        enddo
                        p = indp(7)
                        if( cpot .eq. 'U11' )then
                            call u11poly(indp, geoint, xin, dxin, gin, x7(p), &
                                & dx7(p), fval, dfval)
                        elseif( cpot .eq. 'U22' )then
                            call u22poly(indp, geoint, xin, dxin, gin, g1in, &
                                & x7(p), dx7(p), fval, dfval)
                        elseif( cpot .eq. 'U12' )then
                            call u12poly(indp, geoint, gin, x7(p), &
                                & dx7(p), fval, dfval)
                        endif
                        v = v + fval * coeff(k)
                        do i=1, 7
                            !           apply chain rule
                            do j=1, 7
                                gv(i) = gv(i) + dfval(j)*dgeoint(j,i)*coeff(k)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine sumpoly

    !> @brief Evaluate symOpt for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine symopt( n, ip, s )
        implicit none
        integer          :: n, ip(n)
        double precision :: s(n,n)
        integer          :: i, j
        do i=1, n
            do j=1, n
                s(i,j) = 0.d0
            enddo
        enddo
        do i=1, n
            s(i,ip(i)) = 1.d0
        enddo
    end subroutine symopt

    !> @brief Evaluate getSym7g for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine getsym7g( s, x, y, dy )
        implicit none
        integer, parameter :: n = 7
        double precision   :: x(n)
        double precision   :: y(n)
        ! dy(i,j) = dy(i)/dx(j)
        double precision :: dy(n,n)
        integer          :: i, j, k
        double precision :: s(3,3)
        do i=1, n
            y(i) = 0.d0
            do j=1, n
                dy(i,j) = 0.d0
            enddo
        enddo
        do i=1, 3
            do j=1, 3
                y(i)        = y(i) + s(i,j)*x(j)
                y(i+3)      = y(i+3) + s(i,j)*x(j+3)
                dy(i,j)     = s(i,j)
                dy(i+3,j+3) = s(i,j)
            enddo
        enddo
        y(7)    = x(7)
        dy(7,7) = 1.d0
    end subroutine getsym7g

    !> @brief Evaluate getSym7 for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine getsym7( s, x, y )
        implicit none
        integer, parameter :: n = 7
        double precision   :: x(n)
        double precision   :: y(n)
        integer            :: i, j, k
        double precision   :: s(3,3)
        do i=1, n
            y(i) = 0.d0
        enddo
        do i=1, 3
            do j=1, 3
                y(i)   = y(i) + s(i,j)*x(j)
                y(i+3) = y(i+3) + s(i,j)*x(j+3)
            enddo
        enddo
        y(7) = x(7)
    end subroutine getsym7

    !> @brief Evaluate getSym3 for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine getsym3( s, x, y )
        implicit none
        double precision :: x(3)
        double precision :: y(3)
        integer          :: i, j, k
        double precision :: s(3,3)
        do i=1, 3
            y(i) = 0.d0
        enddo
        do i=1, 3
            do j=1, 3
                y(i) = y(i) + s(i,j)*x(j)
            enddo
        enddo
    end subroutine getsym3

    !> @brief Evaluate repulsive for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine repulsive( rr, erep, derep, cpot )
        ! Additional repuslive terms
        implicit none
        double precision :: rr(6), erep, derep(6)
        character(len=3) :: cpot
        double precision :: bnh1, bnh2, bhh
        !     exponental coefficient
        double precision :: alpnh1, alpnh2, alphh
        !     Center of Gaussian Function
        double precision :: r0nh, r0hh
        !      parameter (R00 = 0.0d0)
        double precision :: r0
        double precision :: erepi, erepnh, erephh
        ! second additional repulsive terms added to NH repulsion
        double precision :: erepnh1, erepnh2
        ! Used for asymtotic correction
        double precision :: ecorrlimit
        ! Derivatives: derepXH: 1-3 NH repul. derivs, 4-6, HH repul. derivs
        double precision :: derepxh(6), derepnh1, derepnh2, derepi
        double precision :: decorrlimit
        ! U11
        !     DE between the dissociation limit and exp
        double precision :: deu11
        double precision :: shiftu11
        ! U22
        !     DE between the dissociation limit and exp
        double precision :: deu22
        !     energy shift to make the U22 (==V2) planar minimum in agreement
        !     with exp. All the energy are relative to V1
        double precision :: shiftu22
        double precision :: deltae
        double precision :: eshift
        !     Re for the correction function
        !     1.024d0: experimental NH distance for NH2 (2B1, ground state)
        double precision :: ecor1
        double precision :: decor1(6)
        !     gamma
        !      parameter (gam  = 1.81d0 )
        double precision, parameter :: gam   = 3.00d0
        double precision, parameter :: recor = 2.0d0*1.024d0
        integer                     :: i, j
        !  misc
        double precision :: tmp
        double precision :: dtmp
        double precision :: dr, dr1, dr2
        double precision :: rr1
        !      Re0 = 1.004d0
        r0nh = -0.2d0
        r0hh = -0.2d0
        !              Exp        Theor
        shiftu11 = 0.00189d0
        deu11    = 0.21967d0
        shiftu22 = 0.04726d0
        deu22    = 0.19700d0
        if( cpot .eq. 'U11' )then
            ! Adjust AlpNH and BNH to make sure that  (U11)
            !     rNH = 0.4, erepNH= 10.00
            !     rNH = 0.5, erepNH=  1.91
            !     rNH = 0.6, erepNH=  0.31
            !     rNH = 0.7, erepNH=  0.04
            !     rNH = 0.8, erepNH=  0.00
            alpnh2 = 11.00d0
            bnh2   = 209.83d0
            eshift = shiftu11
            deltae = deu11
        elseif( cpot .eq. 'U22' )then
            ! Adjust AlpNH and BNH to make sure that  (U22)
            !     rNH = 0.4, erepNH= 20.00
            !     rNH = 0.5, erepNH=  3.83
            !     rNH = 0.6, erepNH=  0.61
            !     rNH = 0.7, erepNH=  0.08
            !     rNH = 0.8, erepNH=  0.01
            alpnh2 = 11.00d0
            bnh2   = 419.66d0
            eshift = shiftu22
            deltae = deu22
        endif
        ! Fitted
        bnh1       = .32363d+02
        bhh        = .12743d+02
        alpnh1     = .40668d+01
        alphh      = .41735d+01
        erep       = 0.d0
        erepnh1    = 0.d0
        erepnh2    = 0.d0
        erephh     = 0.d0
        ecorrlimit = 0.d0
        ecor1      = 0.d0
        if( cpot .ne. 'U12' )then
            do i=1, 3
                dr         = rr(i)-r0nh
                rr1        = 1.d0/rr(i)
                erepi      = dexp(-alpnh1*dr*dr) * rr1
                erepnh1    = erepnh1+ erepi
                derepnh1   = -erepi*( 2.d0*alpnh1*dr + rr1 )
                erepi      = dexp(-alpnh2*dr*dr) * rr1
                erepnh2    = erepnh2 + erepi
                derepnh2   = -erepi*( 2.d0*alpnh2*dr + rr1 )
                derepxh(i) = bnh1*derepnh1 + bnh2*derepnh2
                ! Correction to right limit
                !        erepi   = 0.5d0*(1.d0+dtanh(gam*sum[rr(i)-Recor]) )
                tmp = dtanh(gam*(rr(i)-recor))
                !       d(tanh) = (sech)**2 = 1-(tanh)**2
                dtmp      = gam*(1.d0-tmp*tmp)
                erepi     = 0.5d0*( 1.d0 + tmp )
                ecor1     = ecor1 + erepi
                decor1(i) = 0.5d0*dtmp
            enddo
            erepnh = bnh1*erepnh1 + bnh2*erepnh2
            do i=4, 6
                dr         = rr(i)-r0hh
                rr1        = 1.d0/rr(i)
                erepi      = dexp(-alphh*dr*dr) * rr1
                erephh     = erephh + erepi
                derepxh(i) = -bhh*erepi*( 2.d0*alphh*dr + rr1 )
                decor1(i)  = 0.d0
            enddo
            erephh = bhh * erephh
            do i=1, 6
                decorrlimit = deltae*decor1(i)
                derep(i)    = derepxh(i) + decorrlimit
            enddo
            ecorrlimit = eshift + deltae*ecor1
            erep       = erepnh + erephh + ecorrlimit
        endif
    end subroutine repulsive

    !> @brief Evaluate U11poly for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine u11poly( ind0, geoint, xin, dxin, gin, x7, dx7, ans, g_ans )
        implicit none
        integer          :: ind0(7)
        double precision :: geoint(7), gin(3)
        double precision :: xin(7), dxin(3)
        double precision :: ans
        double precision :: g_ans(7)
        ! Parameter List
        double precision            :: frb, frc, alp
        double precision, parameter :: pi     = 3.1415926535897930d0
        double precision, parameter :: deg2au = pi/180.d0
        ! Experimental value: Landolt-Bötein, Numerical Data and Function
        ! Relationships in Science and Technology II(7,2) (Springer, New York,
        ! 1976). Cited in J. Pesonen, A. Miani, and L. Halonen, J. Chem. Phys.
        ! 115, 1243 (2001)
        double precision, parameter :: gta = 106.68*pi/180.0d0
        ! CoS(GTa) = CoS(106.28)
        double precision, parameter :: cosgta = -0.2870261592d0
        double precision, parameter :: hba    = 0.5d0*pi
        ! Local Variables
        double precision :: g1
        double precision :: g2
        double precision :: g3
        double precision :: x1
        double precision :: x2
        double precision :: x3
        double precision :: x4
        double precision :: x5
        double precision :: x6
        double precision :: x7
        double precision :: y(7)
        double precision :: beta, r1, r2, r3, the1, the2, the3
        integer          :: i, j, k, l, m, n, p
        integer          :: ind
        ! Derivatives
        double precision :: y123, y4567
        double precision :: dg1, dg2, dg3
        double precision :: dx1, dx2, dx3, dx4, dx5, dx6, dx7
        double precision :: dy71, dy72, dy73, dy77
        double precision :: dy(3)
        double precision :: dy42, dy43, dy44
        double precision :: dy51, dy53, dy55
        double precision :: dy61, dy62, dy66
        double precision :: tmp
        !  parameters to be optimized:
        !                 GTa(theta0),FRb(gamma),FRc(r0),ALP(delta)
        frb = .1872506360d+01
        frc = .4278236310d+00
        alp = .1009053173d+00
        ! geoInt contains the internals in the following format
        ! geoInt(1) : r1
        ! geoInt(2) : r2
        ! geoInt(3) : r3
        ! geoInt(4) : theta1
        ! geoInt(5) : theta2
        ! geoInt(6) : theta3
        ! geoInt(7) : beta
        l = ind0(4)
        m = ind0(5)
        n = ind0(6)
        p = ind0(7)
        !        r1 = geoInt(1)
        !        r2 = geoInt(2)
        !        r3 = geoInt(3)
        !      the1 = geoInt(4)*deg2au
        !      the2 = geoInt(5)*deg2au
        !      the3 = geoInt(6)*deg2au
        !      beta = geoInt(7)*deg2au
        !      g1= dexp( -ALP * r1*r1 )
        !      g2= dexp( -ALP * r2*r2 )
        !      g3= dexp( -ALP * r3*r3 )
        g1 = gin(1)
        g2 = gin(2)
        g3 = gin(3)
        ! dgi = dgi/dri (i=1,2,3)
        dg1 = -2.d0*alp*geoint(1)*g1
        dg2 = -2.d0*alp*geoint(2)*g2
        dg3 = -2.d0*alp*geoint(3)*g3
        ! Energy and Derivatives
        ! f(r)  = 1 - exp[-B(r-C)]
        ! dxi = dxi/dri (i=1,2,3) A,B,C are constants
        ! df/dr = B exp[-B(r-C)]
        do i=1, 3
            ind = ind0(i)
            if( ind .eq. 0 )then
                dy(i) = 0.d0
                y(i)  = 1.d0
            else
                !          x1   = dexp(-FRb*(geoInt(i) - FRc))
                x1    = xin(i)
                dx1   = frb*x1
                x1    = 1.d0 - x1
                tmp   = x1**(ind-1)
                dy(i) = dble(ind)* tmp * dx1
                y(i)  = tmp*x1
            endif
        enddo
        ! De = 1.d0 in GT function, geoInt(4-6) are theta1,theta2,theta3
        !     x* = dcos(GTa) - dcos(  theta )
        ! dy4i = dy(4)/dri (i=1,2,3), dy44 = dy(4)/dtheta1
        if( l .eq. 0 )then
            dy42 = 0.d0
            dy43 = 0.d0
            dy44 = 0.d0
            y(4) = 1.d0
        else
            !         x4  = cosGTa - dcos(  the1 )
            !         dx4  = dsin(  the1 ) * deg2au
            x4   = xin(4)
            dx4  = dxin(1)
            tmp  = (x4*g2*g3)**(l-1)
            dy42 = dble(l)* tmp *  x4 * dg2 *  g3
            dy43 = dble(l)* tmp *  x4 *  g2 * dg3
            dy44 = dble(l)* tmp * dx4 *  g2 *  g3
            y(4) = tmp*x4*g2*g3
        endif
        ! dy5i = dy(5)/dri (i=1,2,3), dy55 = dy(5)/dtheta2
        if( m .eq. 0 )then
            dy51 = 0.d0
            dy53 = 0.d0
            dy55 = 0.d0
            y(5) = 1.d0
        else
            !         x5  = cosGTa - dcos(  the2 )
            !         dx5  = dsin(  the2 ) * deg2au
            x5   = xin(5)
            dx5  = dxin(2)
            tmp  = (x5*g1*g3)**(m-1)
            dy51 = dble(m)* tmp *  x5 * dg1 *  g3
            dy53 = dble(m)* tmp *  x5 *  g1 * dg3
            dy55 = dble(m)* tmp * dx5 *  g1 *  g3
            y(5) = tmp*x5*g1*g3
        endif
        ! dy6i = dy(6)/dri (i=1,2,3), dy66 = dy(6)/dtheta3
        if( n .eq. 0 )then
            dy61 = 0.d0
            dy62 = 0.d0
            dy66 = 0.d0
            y(6) = 1.d0
        else
            !         x6  = cosGTa - dcos(  the3 )
            !         dx6  = dsin(  the3 ) * deg2au
            x6   = xin(6)
            dx6  = dxin(3)
            tmp  = (x6*g1*g2)**(n-1)
            dy61 = dble(n)* tmp *  x6 * dg1 *  g2
            dy62 = dble(n)* tmp *  x6 *  g1 * dg2
            dy66 = dble(n)* tmp * dx6 *  g1 *  g2
            y(6) = tmp*x6*g1*g2
        endif
        !  geoInt(7) beta; HBa PI/2
        !      tmp = dble(p) * (beta-HBa)
        !      x7  = dcos( tmp )
        !   dx7/dbeta, not in degree
        !      dx7 = -dble(p)*dsin( tmp ) * deg2au
        ! dy7i = dy(7)/dri (i=1,2,3), dy77 = dy(7)/dbeta
        tmp   = 1.d0-x7
        dy71  = tmp* dg1 *  g2 *  g3
        dy72  = tmp*  g1 * dg2 *  g3
        dy73  = tmp*  g1 *  g2 * dg3
        dy77  = -dx7*  g1 *  g2 *  g3
        y(7)  = 1.0d0 + tmp*g1*g2*g3
        y123  = y(1)*y(2)*y(3)
        y4567 = y(4)*y(5)*y(6)*y(7)
        ans   = y123 * y4567
        ! Derivatives ans = y1 * y2 * y3 * y4 * y5 * y6 * y7
        !     dans / dr1
        g_ans(1) = dy(1)* y(2) * y(3) * y4567       + y123 * y(4) * dy51 * y(6) * y(7) + y123 * y(4) * y(5) * dy61 * &
            & y(7) + y123 * y(4) * y(5) * y(6) * dy71
        !     dans / dr2
        g_ans(2) = y(1) * dy(2)* y(3) * y4567       + y123 * dy42 * y(5) * y(6) * y(7) + y123 * y(4) * y(5) * dy62 * &
            & y(7) + y123 * y(4) * y(5) * y(6) * dy72
        !     dans / dr3
        g_ans(3) = y(1) * y(2) * dy(3)* y4567       + y123 * dy43 * y(5) * y(6) * y(7) + y123 * y(4) * dy53 * y(6) * &
            & y(7) + y123 * y(4) * y(5) * y(6) * dy73
        !     dans / dtheta1
        g_ans(4) = y123 * dy44 * y(5) * y(6) * y(7)
        !     dans / dtheta2
        g_ans(5) = y123 * y(4) * dy55 * y(6) * y(7)
        !     dans / dtheta3
        g_ans(6) = y123 * y(4) * y(5) * dy66 * y(7)
        !     dans / dbeta
        g_ans(7) = y123 * y(4) * y(5) * y(6) * dy77
    end subroutine u11poly

    !> @brief Evaluate U22poly for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine u22poly( ind0, geoint, xin, dxin, gin, g1in, x7, dx7, ans, g_ans )
        implicit none
        integer          :: ind0(7)
        double precision :: geoint(7), gin(3), g1in(3)
        double precision :: xin(7), dxin(3)
        double precision :: ans
        double precision :: g_ans(7)
        ! Parameter List
        double precision :: pi
        double precision :: frb, frc, gta, hba, alp, re, frd, cosgta
        !      double precision deg2au
        !      parameter (  PI     = 3.1415926535897930d0        )
        !      parameter (  deg2au = PI/180.d0                   )
        !      parameter (  GTa    = 2.0d0*PI/3.0d0              )
        ! CoS(GTa) = CoS(120.0) = -1/2
        !      parameter (  cosGTa = -0.5d0                      )
        !      parameter (  HBa    = 0.5d0*PI                    )
        ! Correction to umbralla bending vibration mode
        double precision :: dk
        ! FRd, and Re: used for Gaussian function
        ! Local varibales
        double precision :: g1
        double precision :: g2
        double precision :: g3
        double precision :: x1, x11
        double precision :: x2, x21
        double precision :: x3, x31
        double precision :: x4
        double precision :: x5
        double precision :: x6
        double precision :: x7
        double precision :: y(7)
        double precision :: beta, r1, r2, r3, the1, the2, the3
        integer          :: i, j, k, l, m, n, p
        integer          :: ind
        ! Derivatives
        double precision :: y123, y4567
        double precision :: dg1, dg2, dg3, g1g2, g2g3, g1g3
        double precision :: dx1, dx2, dx3, dx4, dx5, dx6, dx7
        double precision :: dy71, dy72, dy73, dy77
        double precision :: dy(3)
        double precision :: dy11, dy12, dy13
        double precision :: dy21, dy22, dy23
        double precision :: dy31, dy32, dy33
        double precision :: dy42, dy43, dy44
        double precision :: dy51, dy53, dy55
        double precision :: dy61, dy62, dy66
        double precision :: tmp, dtmp(3)
        double precision :: dr
        ! Wed Oct 18 14:57:07 CDT 2006, with additional NH2 (linear)
        frb = .1633658498d+01
        frc = .1987056751d+00
        alp = .3037603850d+00
        re  = .7147714247d+00
        frd = .4546981630d+01
        dk  = -.2617146615d+00
        ! geoInt contains the internals in the following format
        ! geoInt(1) : r1
        ! geoInt(2) : r2
        ! geoInt(3) : r3
        ! geoInt(4) : theta1
        ! geoInt(5) : theta2
        ! geoInt(6) : theta3
        ! geoInt(7) : beta
        l = ind0(4)
        m = ind0(5)
        n = ind0(6)
        p = ind0(7)
        !        r1 = geoInt(1)
        !        r2 = geoInt(2)
        !        r3 = geoInt(3)
        !      the1 = geoInt(4)*deg2au
        !      the2 = geoInt(5)*deg2au
        !      the3 = geoInt(6)*deg2au
        !      beta = geoInt(7)*deg2au
        !      g1 = dexp( -ALP * r1*r1 )
        !      g2 = dexp( -ALP * r2*r2 )
        !      g3 = dexp( -ALP * r3*r3 )
        g1   = gin(1)
        g2   = gin(2)
        g3   = gin(3)
        g1g2 = g1*g2
        g1g3 = g1*g3
        g2g3 = g2*g3
        ! dgi = dgi/dri (i=1,2,3)
        dg1 = -2.d0*alp*geoint(1)*g1
        dg2 = -2.d0*alp*geoint(2)*g2
        dg3 = -2.d0*alp*geoint(3)*g3
        ! Energy and Derivatives
        ! f(r)  = 1 - exp[-B(r-C)]
        ! dxi = dxi/dri (i=1,2,3) A,B,C are constants
        ! df/dr = B exp[-B(r-C)]
        do i=1, 3
            ind = ind0(i)
            if( ind .eq. 0 )then
                dy(i) = 0.d0
                y(i)  = 1.d0
            elseif( ind .le. 2 )then
                !          x1   = dexp(-FRb*(geoInt(i) - FRc))
                x1    = xin(i)
                dx1   = frb*x1
                x1    = 1.d0 - x1
                tmp   = x1**(ind-1)
                dy(i) = dble(ind)* tmp * dx1
                y(i)  = tmp*x1
            endif
        enddo
        ! Li on Suggested by Don. (R-Re)^2 exp(-a (R-Re)^2)
        if( ind0(1) .eq. 3 )then
            dr = geoint(1) - re
            !        x1   = dexp(-FRd*dr*dr)
            x1   = g1in(1)
            dx1  = -2.d0*frd*dr*x1
            y(1) = x1 * g2g3
            dy11 = dx1* g2g3
            dy12 = x1 * dg2*g3
            dy13 = x1 * g2*dg3
        else
            dy11 = dy(1)
            dy12 = 0.d0
            dy13 = 0.d0
        endif
        if( ind0(2) .eq. 3 )then
            dr = geoint(2) - re
            !        x1   = dexp(-FRd*dr*dr)
            x1   = g1in(2)
            dx1  = -2.d0*frd*dr*x1
            y(2) = x1 * g1g3
            dy21 = x1 * dg1*g3
            dy22 = dx1* g1g3
            dy23 = x1 * g1*dg3
        else
            dy21 = 0.d0
            dy22 = dy(2)
            dy23 = 0.d0
        endif
        if( ind0(3) .eq. 3 )then
            dr = geoint(3) - re
            !        x1   = dexp(-FRd*dr*dr)
            x1   = g1in(3)
            dx1  = -2.d0*frd*dr*x1
            y(3) = x1 * g1g2
            dy31 = x1 * dg1*g2
            dy32 = x1 * g1*dg2
            dy33 = dx1* g1g2
        else
            dy31 = 0.d0
            dy32 = 0.d0
            dy33 = dy(3)
        endif
        !     x* = dcos(GTa) - dcos(  theta )
        !     GTa = 120.0 degree, dcos(GTa) = -0.5d0
        ! dy5i = dy(5)/dri (i=1,2,3), dy55 = dy(5)/dtheta2
        if( l .eq. 0 )then
            dy42 = dg2*g3
            dy43 = g2*dg3
            dy44 = 0.d0
            y(4) = 1.d0 + g2g3
        else
            !          x4  = cosGTa - dcos(  the1 )
            x4  = xin(4)
            dx4 = dxin(1)
            !         dx4  = dsin(  the1 ) * deg2au
            tmp  = x4**(l-1)
            dy42 = tmp *  x4 * dg2 *  g3
            dy43 = tmp *  x4 *  g2 * dg3
            dy44 = dble(l)* tmp * dx4 *  g2g3
            y(4) = 1.d0+tmp*x4*g2g3
        endif
        if( m .eq. 0 )then
            dy51 = dg1*g3
            dy53 = g1*dg3
            dy55 = 0.d0
            y(5) = 1.d0 + g1g3
        else
            !          x5  = cosGTa - dcos(  the2 )
            !         dx5  = dsin(  the2 ) * deg2au
            x5   = xin(5)
            dx5  = dxin(2)
            tmp  = x5**(m-1)
            dy51 = tmp *  x5 * dg1 *  g3
            dy53 = tmp *  x5 *  g1 * dg3
            dy55 = dble(m)* tmp * dx5 *  g1g3
            y(5) = 1.d0+tmp*x5*g1g3
        endif
        ! dy6i = dy(6)/dri (i=1,2,3), dy66 = dy(6)/dtheta3
        if( n .eq. 0 )then
            dy61 = dg1*g2
            dy62 = g1*dg2
            dy66 = 0.d0
            y(6) = 1.d0 + g1g2
        else
            !         x6  = cosGTa - dcos(  the3 )
            !         dx6  = dsin(  the3 ) * deg2au
            x6   = xin(6)
            dx6  = dxin(3)
            tmp  = x6**(n-1)
            dy61 = tmp *  x6 * dg1 *  g2
            dy62 = tmp *  x6 *  g1 * dg2
            dy66 = dble(n)* tmp * dx6 *  g1g2
            y(6) = 1.d0+tmp*x6*g1g2
        endif
        !  geoInt(7) beta; HBa PI/2
        !      y(7) = 1 + (1-x7)*g1*g2*g3
        !      tmp = dble(p) * (beta-HBa)
        !      x7  = dcos( tmp )
        !   dx7/dbeta, not in degree
        !      dx7 = -dble(p)*dsin( tmp ) * deg2au
        !      tmp = beta-HBa
        !      x7 = dcos( dble(p) * tmp ) + 0.5d0*dk*tmp*tmp/(0.05d0+tmp*tmp)
        !      dx7 = -dble(p)*dsin( dble(p) * tmp ) * deg2au
        !     &      + dk*tmp*deg2au/(0.05d0+tmp*tmp)
        !     &      - dk*tmp**3*deg2au/(0.05d0+tmp*tmp)**2
        ! dy7i = dy(7)/dri (i=1,2,3), dy77 = dy(7)/dbeta
        tmp   = 1.d0-x7
        dy71  = tmp* dg1 *  g2g3
        dy72  = tmp*  g1g3 * dg2
        dy73  = tmp*  g1g2 * dg3
        dy77  = -dx7*  g1 *  g2g3
        y(7)  = 1.0d0 + tmp*g1*g2g3
        y123  = y(1)*y(2)*y(3)
        y4567 = y(4)*y(5)*y(6)*y(7)
        ans   = y123 * y4567
        ! Derivatives ans = y1 * y2 * y3 * y4 * y5 * y6 * y7
        !     dans / dr1
        g_ans(1) = dy11 * y(2) * y(3) * y4567       + y(1) * dy21 * y(3) * y4567       + y(1) * y(2) * dy31 * &
            & y4567       + y123 * y(4) * dy51 * y(6) * y(7) + y123 * y(4) * y(5) * dy61 * y(7) + y123 * y(4) * y(5) * &
            & y(6) * dy71
        !     dans / dr2
        g_ans(2) = dy12 * dy(2)* y(3) * y4567       + y(1) * dy22 * y(3) * y4567       + y(1) * y(2) * dy32 * &
            & y4567       + y123 * dy42 * y(5) * y(6) * y(7) + y123 * y(4) * y(5) * dy62 * y(7) + y123 * y(4) * y(5) * &
            & y(6) * dy72
        !     dans / dr3
        g_ans(3) = dy13 * y(2) * dy(3)* y4567       + y(1) * dy23 * y(3) * y4567       + y(1) * y(2) * dy33 * &
            & y4567       + y123 * dy43 * y(5) * y(6) * y(7) + y123 * y(4) * dy53 * y(6) * y(7) + y123 * y(4) * y(5) * &
            & y(6) * dy73
        !     dans / dtheta1
        g_ans(4) = y123 * dy44 * y(5) * y(6) * y(7)
        !     dans / dtheta2
        g_ans(5) = y123 * y(4) * dy55 * y(6) * y(7)
        !     dans / dtheta3
        g_ans(6) = y123 * y(4) * y(5) * dy66 * y(7)
        !     dans / dbeta
        g_ans(7) = y123 * y(4) * y(5) * y(6) * dy77
    end subroutine u22poly

    !> @brief Evaluate U12poly for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine u12poly( ind0, geoint, gin, hbu12, d7, ans, g_ans )
        implicit none
        integer          :: ind0(7)
        double precision :: geoint(7), gin(3)
        double precision :: ans
        ! Derivatives of 7 internal coordinates
        double precision :: g_ans(7)
        ! Parameter List
        double precision, parameter :: pi = 3.1415926535897930d0
        !! Parameter List (U12)
        ! FRU12a    parameters used to define angular for U22
        ! FRU12b    parameters used to define angular for U22
        ! FRU12c    parameters used to define angular for U22
        double precision            :: fru12a, fru12c
        double precision            :: alp
        double precision, parameter :: deg2au = pi/180.d0
        double precision, parameter :: fru12b = 1.00d0
        !      parameter ( FRU12c = 2.0150d0        )
        double precision, parameter :: gtu12a = 2.0d0*pi/3.0d0
        double precision, parameter :: hbu12a = 0.5d0*pi
        ! Misc
        integer          :: i, j, ind
        double precision :: z1(6), tz, tf, z123, z456
        double precision :: hbu12, fru121, fru122, fru123
        double precision :: theta, beta
        double precision :: y(7)
        double precision :: tmp
        ! Derivatives
        double precision :: dz1(6)
        double precision :: d7, dfr1, dfr2, dfr3
        !  parameters to be optimized:
        !                 FRU12c(r04)
        fru12c = .1501070926d+01
        ! geoInt contains the internals in the following format
        ! geoInt(1) : r1
        ! geoInt(2) : r2
        ! geoInt(3) : r3
        ! geoInt(4) : theta1
        ! geoInt(5) : theta2
        ! geoInt(6) : theta3
        ! geoInt(7) : beta
        do i=1, 3
            ind = ind0(i)
            if( ind .eq. 0 )then
                dz1(i) = 0.d0
                z1(i)  = 1.d0
            else
                z1(i)  = geoint(i)**(ind-1)
                dz1(i) = dble(ind)*z1(i)
                z1(i)  = z1(i)*geoint(i)
            endif
        enddo
        do i=4, 6
            ind = ind0(i)
            if( ind .eq. 0 )then
                dz1(i) = 0.d0
                z1(i)  = 1.d0
            else
                !  Unit degree/A
                theta  = geoint(i) * deg2au
                z1(i)  = theta**(ind-1)
                dz1(i) = dble(ind)*z1(i) * deg2au
                z1(i)  = z1(i)*theta
            endif
        enddo
        !       beta  = geoInt(7)*deg2au
        !       tmp   =  dble(ind0(7))*(beta-HBU12a)
        !       HBU12 = -dsin( tmp )
        !  d_HBU12 / d geoInt(7)
        !       d7  = -dcos( tmp ) * dble(ind0(7)) * deg2au
        ! f(r) = 1 - tanh[-B(r-C)]
        !  d tanh(x) = sech(x)^2 = 1 - tanh(x)^2
        ! df(r)/dr =  - {1 - tanh[-B(r-C)] }^2 * B
        !       tmp    = dtanh( FRU12b*(geoInt(1)-FRU12c) )
        tmp    = gin(1)
        fru121 = 1.d0 - tmp
        dfr1   = -(1.d0 - tmp*tmp)*fru12b
        !       tmp    = dtanh( FRU12b*(geoInt(2)-FRU12c) )
        tmp    = gin(2)
        fru122 = 1.d0 - tmp
        dfr2   = -(1.d0 - tmp*tmp)*fru12b
        !       tmp    = dtanh( FRU12b*(geoInt(3)-FRU12c) )
        tmp    = gin(3)
        fru123 = 1.d0 - tmp
        dfr3   = -(1.d0 - tmp*tmp)*fru12b
        z123   = z1(1)*z1(2)*z1(3)
        z456   = z1(4)*z1(5)*z1(6)
        tz     = z123*z456
        tf     = fru121 * fru122 * fru123
        ans    = tz * tf * hbu12
        !       ans =  z1* HBU12 * FRU121 * FRU122 * FRU123
        g_ans(1) = dz1(1)*z1(2)*z1(3)*z456*tf*hbu12 + tz * hbu12 * dfr1 * fru122 * fru123
        g_ans(2) = z1(1)*dz1(2)*z1(3)*z456*tf*hbu12 + tz * hbu12 * fru121 * dfr2 * fru123
        g_ans(3) = z1(1)*z1(2)*dz1(3)*z456*tf*hbu12 + tz * hbu12 * fru121 * fru122 * dfr3
        g_ans(4) = z123*dz1(4)*z1(5)*z1(6)*tf*hbu12
        g_ans(5) = z123*z1(4)*dz1(5)*z1(6)*tf*hbu12
        g_ans(6) = z123*z1(4)*z1(5)*dz1(6)*tf*hbu12
        g_ans(7) = tz * tf * d7
    end subroutine u12poly

    !> @brief Evaluate getAdiabat for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine getadiabat( nv, u11, u22, u12, v1, v2, &
            & gu11, gu22, gu12, gv1, gv2 )
        implicit none
        ! input
        integer          :: nv
        double precision :: u11, u22, u12
        double precision :: gu11(nv), gu22(nv), gu12(nv)
        ! Output
        double precision :: v1, v2
        double precision :: gv1(nv), gv2(nv)
        integer          :: i, j
        double precision :: x, y, dx, dy
        x  = u11+u22
        y  = dsqrt( ((u22-u11)**2) + (4.0d0*u12*u12) )
        v1 = 0.5d0*(x - y)
        v2 = 0.5d0*(x + y)
        ! Derivatives
        do i=1, nv
            dx     = gu11(i) + gu22(i)
            dy     = ( (u22-u11)*( gu22(i)-gu11(i) ) + 4.d0*u12*gu12(i) ) / y
            gv1(i) = 0.5d0 * (dx - dy)
            gv2(i) = 0.5d0 * (dx + dy)
        enddo
    end subroutine getadiabat

    !> @brief Evaluate calcindU1 for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine calcindu1( ind, m1, m2, m3, m4, m5, m6, m7, nc )
        implicit none
        integer :: nc
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: m6
        integer :: m7
        integer :: ind(nc,7)
        integer :: i1, i2, i3, i4, i5, i6, i7, inc
        inc = 0
        do i7=1, m7
            do i1=0, m1
                do i2=i1, m2
                    do i3=i2, m3
                        do i4=0, m4
                            do i5=i4, m5
                                do i6=i5, m6
                                    if( i1+i2+i3+i4+i5+i6 .eq. 0 ) cycle
                                    if( i1+i2+i3 .gt. 4 ) cycle
                                    if( i4+i5+i6 .gt. 2 ) cycle
                                    inc        = inc + 1
                                    ind(inc,1) = i1
                                    ind(inc,2) = i2
                                    ind(inc,3) = i3
                                    ind(inc,4) = i4
                                    ind(inc,5) = i5
                                    ind(inc,6) = i6
                                    ind(inc,7) = i7
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine calcindu1

    !> @brief Evaluate calcindU2 for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine calcindu2( ind, m1, m2, m3, m4, m5, m6, m7, nc )
        implicit none
        integer :: nc
        integer :: mnc
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: m6
        integer :: m7
        integer :: ind(nc,7)
        integer :: i1, i2, i3, i4, i5, i6, i7, inc
        inc = 0
        do i7=0, m7
            do i1=0, m1
                do i2=i1, m2
                    do i3=i2, m3
                        do i4=0, m4
                            do i5=i4, m5
                                do i6=i5, m6
                                    if( i1+i2+i3+i4+i5+i6 .gt. 0 .and. i7 .eq. 0 ) cycle
                                    if( i1+i2+i3+i4+i5+i6 .eq. 0 .and. i7 .gt. 0 ) cycle
                                    if( i1+i2+i3 .gt. 5 ) cycle
                                    if( i4+i5+i6 .gt. 2 ) cycle
                                    if( i1+i2+i3 .ge. 4 )then
                                        if( i1 .eq. 2 .or. i2 .eq. 2 .or.i3 .eq. 2 ) cycle
                                    endif
                                    inc        = inc + 1
                                    ind(inc,1) = i1
                                    ind(inc,2) = i2
                                    ind(inc,3) = i3
                                    ind(inc,4) = i4
                                    ind(inc,5) = i5
                                    ind(inc,6) = i6
                                    ind(inc,7) = i7
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !      NC=iNC
    end subroutine calcindu2

    !> @brief Evaluate calcindU12 for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine calcindu12( ind, m1, m2, m3, m4, m5, m6, m7, nc )
        implicit none
        integer :: nc
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: m6
        integer :: m7
        integer :: ind(nc,7)
        integer :: i1, i2, i3, i4, i5, i6, i7, inc
        inc = 0
        do i7=3, m7
            do i1=0, m1
                do i2=i1, m2
                    do i3=i2, m3
                        do i4=0, m4
                            do i5=i4, m5
                                do i6=i5, m6
                                    if( i1+i2+i3+i4+i5+i6 .eq. 0 ) cycle
                                    if( i1+i2+i3 .gt. 2 ) cycle
                                    if( i4+i5+i6 .gt. 2 ) cycle
                                    inc        = inc + 1
                                    ind(inc,1) = i1
                                    ind(inc,2) = i2
                                    ind(inc,3) = i3
                                    ind(inc,4) = i4
                                    ind(inc,5) = i5
                                    ind(inc,6) = i6
                                    ind(inc,7) = i7
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine calcindu12

    !> @brief Evaluate toBowman for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine tobowman( x, bow, g_bow )
        implicit none
        integer, parameter :: natom = 4
        !     X in one array with x1,y1,z1,...x4,y4,z4
        double precision            :: x(3 * natom)
        double precision            :: bow(7)
        double precision, parameter :: pi = 3.1415926535897930d0
        ! to save the first three bond distances. If they are calculated in any
        ! subroutines, avoid to calculate them again
        ! This function returns the Bowman
        ! coordinates from cartesian coordinates
        ! All angles are in DEGREES
        ! Output:
        ! bow(1) : r1
        ! bow(2) : r2
        ! bow(3) : r3
        ! bow(4) : theta1
        ! bow(5) : theta2
        ! bow(6) : theta3
        ! bow(7) : beta
        ! Input:
        ! coord(1,:) : N
        ! coord(2,:) : H1
        ! coord(3,:) : H2
        ! coord(4,:) : H3
        double precision :: r(3)
        double precision :: beta(3)
        double precision :: theta(3)
        ! Parameter list
        integer, parameter          :: iout1 = 0
        integer, parameter          :: inp1  = 5
        double precision, parameter :: tol   = 1.0d-4
        ! Geometry
        double precision :: cartn(natom,3)
        double precision :: cartp(natom,3)
        ! Unit Vectors
        double precision :: ua(3), ub(3), uc(3)
        double precision :: ut(3)
        double precision :: uba(3), uca(3), ucb(3)
        ! Eq. of Plane
        double precision :: a1
        double precision :: a2
        double precision :: a3
        double precision :: d
        double precision :: dist
        ! Misc
        integer          :: i, j, k, n, m
        double precision :: xt, yt, zt, w
        double precision :: b(3,3)
        double precision :: a(3,3)
        ! Get data from user
        ! Now we will evaluate the direction of the trisector vector
        ! the unit vectors along Ha, Hb and Hc bonds are labelled as
        ! uA, uB, and uC. The trisector direction is refered as uT
        ! cartN has the origin on the N atom
        integer, parameter :: g_pmax_ = 12
        integer            :: g_i_
        double precision   :: d5_b, d3_p, d7_b, d2_v, d3_v, d8_b, d2_b, d3_b, &
            & d1_w, d4_v
        double precision :: d5_v, d1_p, d7_v, d2_p, d4_b, g_cartn(g_pmax_,&
            & natom,3), g_x(g_pmax_,3 * natom), g_d1_w(g_pmax_), g_r(g_pmax_, &
            & 3), g_ua(g_pmax_, 3)
        double precision :: g_ub(g_pmax_,3), g_uc(g_pmax_,3), g_uba(g_pma&
            &x_,3), g_uca(g_pmax_, 3), g_ucb(g_pmax_, 3), g_a(g_pmax_, 3, 3), &
            & g_b(g_pmax_, 3, 3), g_ut(g_pmax_, 3), g_w(g_pmax_), g_beta(g_pmax_ &
            &, 3)
        double precision :: g_a1(g_pmax_), g_a2(g_pmax_), g_a3(g_pmax_), g_&
            &d(g_pmax_), g_xt(g_pmax_), g_yt(g_pmax_), g_zt(g_pmax_), g_cartp(g&
            &_pmax_, natom, 3), g_dist(g_pmax_), g_theta(g_pmax_, 3)
        double precision :: g_bow(g_pmax_,7)
        integer          :: g_ehfid
        !        call ehsfid(g_ehfid, 'toBowman','g_toBowman.f')
        g_ehfid = 0
        do i=1, 12
            do j=1, 12
                if( i .eq. j )then
                    g_x(i,j) = 1.d0
                else
                    g_x(i,j) = 0.d0
                endif
            enddo
        enddo
        do i = 1, natom
            do g_i_ = 1, g_pmax_
                g_cartn(g_i_,i,1) = -g_x(g_i_,1) + g_x(g_i_,3 * i - 2)
            enddo
            cartn(i,1) = x(3 * i - 2) - x(1)
            do g_i_ = 1, g_pmax_
                g_cartn(g_i_,i,2) = -g_x(g_i_,2) + g_x(g_i_,3 * i - 1)
            enddo
            cartn(i,2) = x(3 * i - 1) - x(2)
            do g_i_ = 1, g_pmax_
                g_cartn(g_i_,i,3) = -g_x(g_i_,3) + g_x(g_i_,3 * i)
            enddo
            cartn(i,3) = x(3 * i) - x(3)
        enddo
        ! Get the bond distances
        do i = 2, natom
            d2_v = cartn(i,1) * cartn(i,1)
            d3_p = 2.0d0 * cartn(i,1)
            d4_v = cartn(i,2) * cartn(i,2)
            d2_p = 2.0d0 * cartn(i,2)
            d7_v = cartn(i,3) * cartn(i,3)
            d1_p = 2.0d0 * cartn(i,3)
            do g_i_ = 1, g_pmax_
                g_d1_w(g_i_) = d1_p * g_cartn(g_i_,i,3) + d2_p * g_cartn(g_i_,i,2) + d3_p * g_cartn(g_i_,i,1)
            enddo
            d1_w = d2_v+d4_v+d7_v
            d2_v = sqrt(d1_w)
            ! RVM Corrected
            if( d1_w .gt. 0.0d0 )then
                !          if (d2_v .gt. 0.d0) then
                d1_p = 1.0d0 / (2.0d0 *  d2_v)
            endif
            do g_i_ = 1, g_pmax_
                g_r(g_i_, i - 1) = d1_p * g_d1_w(g_i_)
            enddo
            r(i - 1) = d2_v
        enddo
        ! get the units vectors along Ha, Hb, and Hc
        do i = 1, 3
            d3_v = cartn(2,i) / r(1)
            d2_b = 1.0d0 / r(1)
            d3_b = (-d3_v) / r(1)
            do g_i_ = 1, g_pmax_
                g_ua(g_i_, i) = d3_b * g_r(g_i_, 1) + d2_b * g_cartn(g_i_,2,i)
            enddo
            ua(i) = d3_v
            d3_v  = cartn(3,i) / r(2)
            d2_b  = 1.0d0 / r(2)
            d3_b  = (-d3_v) / r(2)
            do g_i_ = 1, g_pmax_
                g_ub(g_i_,i) = d3_b * g_r(g_i_, 2) + d2_b * g_cartn(g_i_,3,i)
            enddo
            ub(i) = d3_v
            d3_v  = cartn(4,i) / r(3)
            d2_b  = 1.0d0 / r(3)
            d3_b  = (-d3_v) / r(3)
            do g_i_ = 1, g_pmax_
                g_uc(g_i_,i) = d3_b * g_r(g_i_, 3) + d2_b * g_cartn(g_i_,4,i)
            enddo
            uc(i) = d3_v
        enddo
        ! The trisector direction satisfies the following eq
        ! (by definition of it being trisector)
        ! uT . uA = cos(beta)
        ! uT . uB = cos(beta)
        ! uT . uC = cos(beta)
        ! These three equation can be modified to:
        ! uT . (uB-uA) = 0
        ! uT . (uC-uA) = 0
        ! uT . (uC-uB) = 0
        do i = 1, 3
            do g_i_ = 1, g_pmax_
                g_uba(g_i_,i) = -g_ua(g_i_, i) + g_ub(g_i_,i)
            enddo
            uba(i) = ub(i) - ua(i)
            do g_i_ = 1, g_pmax_
                g_uca(g_i_, i) = -g_ua(g_i_, i) + g_uc(g_i_,i)
            enddo
            uca(i) = uc(i) - ua(i)
            do g_i_ = 1, g_pmax_
                g_ucb(g_i_, i) = -g_ub(g_i_,i) + g_uc(g_i_,i)
            enddo
            ucb(i) = uc(i) - ub(i)
        enddo
        ! Find the eq of the plane defined by
        ! points uBA, uCA, and uCB.
        ! The normal to that plane is the trisector
        ! directions
        ! The eq of a plane
        !      Ax + By + Cz + D = 0
        ! defined by points:
        ! (x1,y1,z1) (x2,y2,z2) and (x3,y3,z3) is
        !     | 1 y1 z1 |      | x1 1 z1 |      | x1 y1 1 |   | x1 y1 z1 |
        ! A = | 1 y2 z2 |  B = | x2 1 z2 |  C = | x2 y2 1 |  D = -| x2 y2 z2 |
        !     | 1 y3 z3 |      | x3 1 z3 |      | x3 y3 1 |       | x3 y3 z3 |
        ! The eq of the normal to the plane is given by the vector {A,B,C}
        do i = 1, 3
            do g_i_ = 1, g_pmax_
                g_a(g_i_, 1, i) = g_uba(g_i_,i)
            enddo
            a(1,i) = uba(i)
            do g_i_ = 1, g_pmax_
                g_a(g_i_, 2, i) = g_uca(g_i_, i)
            enddo
            a(2,i) = uca(i)
            do g_i_ = 1, g_pmax_
                g_a(g_i_, 3, i) = g_ucb(g_i_, i)
            enddo
            a(3,i) = ucb(i)
        enddo
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    do g_i_ = 1, g_pmax_
                        g_b(g_i_, j, k) = g_a(g_i_, j, k)
                    enddo
                    b(j,k) = a(j,k)
                enddo
            enddo
            do j = 1, 3
                do g_i_ = 1, g_pmax_
                    g_b(g_i_, j, i) = 0.0d0
                enddo
                b(j,i) = 1.0d0
            enddo
            call g_det3(b, g_b, ut(i), g_ut(1, i))
        enddo
        d4_b = ut(3) + ut(3)
        d7_b = ut(2) + ut(2)
        d8_b = ut(1) + ut(1)
        do g_i_ = 1, g_pmax_
            g_d1_w(g_i_) = d4_b * g_ut(g_i_, 3) + d7_b * g_ut(g_i_, 2) + d8_b * g_ut(g_i_, 1)
        enddo
        d1_w = ut(1) * ut(1) + ut(2) * ut(2) + ut(3) * ut(3)
        d2_v = sqrt(d1_w)
        if( d1_w .gt. 0.0d0 )then
            d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
            !           call ehufDO (9,d1_w, d2_v, d1_p,
            !     +g_ehfid,
            !     +302)
        endif
        do g_i_ = 1, g_pmax_
            g_w(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        w = d2_v
        do i = 1, 3
            d3_v = ut(i) / w
            d2_b = 1.0d0 / w
            d3_b = (-d3_v) / w
            do g_i_ = 1, g_pmax_
                g_ut(g_i_, i) = d3_b * g_w(g_i_) + d2_b * g_ut(g_i_, i)
            enddo
            ut(i) = d3_v
        enddo
        ! Get the value of  beta (angle between any H-N and {A,B,C})
        do g_i_ = 1, g_pmax_
            g_w(g_i_) = ut(3) * g_ua(g_i_, 3) + ua(3) * g_ut(g_i_, 3) + ut (2) * g_ua(g_i_, 2) + ua(2) * g_ut(g_i_, &
                & 2) + ut(1) * g_ua(g_i_, 1 ) + ua(1) * g_ut(g_i_, 1)
        enddo
        w    = ut(1) * ua(1) + ut(2) * ua(2) + ut(3) * ua(3)
        d2_v = acos(w)
        if( abs(w) .lt. 1.0d0 )then
            d1_p = -1.0d0 / sqrt ((1.0d0-w)*(1.0d0+w))
        else
            !           call ehufDO (14,w, d2_v, d1_p,
            !     +g_ehfid,
            !     +335)
        endif
        d4_b = 1.0d0 / pi * 180.0d0 * d1_p
        do g_i_ = 1, g_pmax_
            g_beta(g_i_, 1) = d4_b * g_w(g_i_)
        enddo
        beta(1) = d2_v * 180.0d0 / pi
        !     beta(2) = dacos( dot_product(uT,uB) ) * 180.0d0/PI
        !     beta(3) = dacos( dot_product(uT,uC) ) * 180.0d0/PI
        ! Now we define the eq. of the plane
        !     a1x + a2y + a3z + d = 0
        ! such that the normal to the plane is uT and
        ! the Ha hydrogen lies on that plane
        ! The values of {a1,a2,a3} are the same as
        ! the values of uT(1),uT(2),uT(3). We will use
        ! the coordinates of Ha to find the value of d
        do g_i_ = 1, g_pmax_
            g_a1(g_i_) = g_ut(g_i_, 1)
        enddo
        a1 = ut(1)
        do g_i_ = 1, g_pmax_
            g_a2(g_i_) = g_ut(g_i_, 2)
        enddo
        a2 = ut(2)
        do g_i_ = 1, g_pmax_
            g_a3(g_i_) = g_ut(g_i_, 3)
        enddo
        a3 = ut(3)
        do g_i_ = 1, g_pmax_
            g_w(g_i_) = ut(3) * g_cartn(g_i_,2,3) + cartn(2,3) * g_ut(g_i_, 3) + ut(2) * g_cartn(g_i_,2,2) + cartn(2, &
                & 2) * g_ut(g_i_, 2 ) + ut(1) * g_cartn(g_i_,2,1) + cartn(2,1) * g_ut(g_i_, 1)
        enddo
        w = ut(1) * cartn(2,1) + ut(2) * cartn(2,2) + ut(3) * cartn(2,3)
        do g_i_ = 1, g_pmax_
            g_d(g_i_) = -g_w(g_i_)
        enddo
        d = -w
        ! Project all the atoms on the plane
        ! The projected geometry is known as cartP
        ! The project of a point A on the plane P is denoted by point B
        ! The relation between rA and rB is given as :
        ! rB = rA - (dist*uT)
        ! where dist is the distance of point A from the plane.
        do i = 1, natom
            do g_i_ = 1, g_pmax_
                g_xt(g_i_) = g_cartn(g_i_,i,1)
            enddo
            xt = cartn(i,1)
            do g_i_ = 1, g_pmax_
                g_yt(g_i_) = g_cartn(g_i_,i,2)
            enddo
            yt = cartn(i,2)
            do g_i_ = 1, g_pmax_
                g_zt(g_i_) = g_cartn(g_i_,i,3)
            enddo
            zt = cartn(i,3)
            call g_distplane(a1, g_a1, a2, g_a2, a3, g_a3, d, g_d, xt, g_x&
                &t, yt, g_yt, zt, g_zt, dist, g_dist)
            do j = 1, 3
                do g_i_ = 1, g_pmax_
                    g_cartp(g_i_, i, j) = (-dist) * g_ut(g_i_, j) + (-ut(j)) * g_dist(g_i_) + g_cartn(g_i_,i,j)
                enddo
                cartp(i,j) = cartn(i,j) - dist * ut(j)
            enddo
        enddo
        do i = 2, natom
            do j = 1, 3
                do g_i_ = 1, g_pmax_
                    g_cartp(g_i_, i, j) = -g_cartp(g_i_, 1, j) + g_cartp(g_i_, i, j)
                enddo
                cartp(i,j) = cartp(i,j) - cartp(1,j)
            enddo
        enddo
        do i = 1, 3
            do g_i_ = 1, g_pmax_
                g_cartp(g_i_, 1, i) = 0.0d0
            enddo
            cartp(1,i) = 0.0d0
        enddo
        ! The theta are labelled as:
        ! H1-N-H2 : theta(1)
        ! H1-N-H3 : theta(2)
        ! H2-N-H3 : theta(3)
        ! Note that theta <= 180 because we are using acos
        ! and acos gives angle in the range [-PI,PI]
        ! Also all the angles are +ve
        k = 0
        do j = 2, natom
            do i = j + 1, natom
                k    = k + 1
                d2_v = cartp(i,1) * cartp(i,1)
                d3_p = 2.0d0 * cartp(i,1)
                d4_v = cartp(i,2) * cartp(i,2)
                d2_p = 2.0d0 * cartp(i,2)
                d7_v = cartp(i,3) * cartp(i,3)
                d1_p = 2.0d0 * cartp(i,3)
                do g_i_ = 1, g_pmax_
                    g_d1_w(g_i_) = d1_p * g_cartp(g_i_, i, 3) + d2_p * g_cartp (g_i_, i, 2) + d3_p * g_cartp(g_i_, i, 1)
                enddo
                d1_w = d2_v + d4_v + d7_v
                d2_v = sqrt(d1_w)
                if( d1_w .gt. 0.0d0 )then
                    d1_p = 1.0d0 / (2.0d0 *  d2_v)
                else
                    !               call ehufDO (9,d1_w, d2_v, d1_p,
                    !     +g_ehfid,
                    !     +473)
                endif
                do g_i_ = 1, g_pmax_
                    g_xt(g_i_) = d1_p * g_d1_w(g_i_)
                enddo
                xt   = d2_v
                d2_v = cartp(j,1) * cartp(j,1)
                d3_p = 2.0d0 * cartp(j,1)
                d4_v = cartp(j,2) * cartp(j,2)
                d2_p = 2.0d0 * cartp(j,2)
                d7_v = cartp(j,3) * cartp(j,3)
                d1_p = 2.0d0 * cartp(j,3)
                do g_i_ = 1, g_pmax_
                    g_d1_w(g_i_) = d1_p * g_cartp(g_i_, j, 3) + d2_p * g_cartp (g_i_, j, 2) + d3_p * g_cartp(g_i_, j, 1)
                enddo
                d1_w = d2_v + d4_v + d7_v
                d2_v = sqrt(d1_w)
                if( d1_w .gt. 0.0d0 )then
                    d1_p = 1.0d0 / (2.0d0 *  d2_v)
                else
                    !               call ehufDO (9,d1_w, d2_v, d1_p,
                    !     +g_ehfid,
                    !     +498)
                endif
                do g_i_ = 1, g_pmax_
                    g_yt(g_i_) = d1_p * g_d1_w(g_i_)
                enddo
                yt = d2_v
                do g_i_ = 1, g_pmax_
                    g_zt(g_i_) = cartp(i,3) * g_cartp(g_i_, j, 3) + cartp(j,3) * g_cartp(g_i_, i, 3) + cartp(i,2) * &
                        & g_cartp(g_i_, j, 2) + cartp(j,2) * g_cartp(g_i_, i, 2) + cartp(i,1) * g_cartp(g_i_, j, 1) + &
                        & cartp(j,1) * g_cartp(g_i_, i, 1)
                enddo
                zt   = cartp(i,1) * cartp(j,1) + cartp(i,2) * cartp(j,2) + cartp(i,3) * cartp(j,3)
                d4_v = xt * yt
                d5_v = zt / d4_v
                d2_b = 1.0d0 / d4_v
                d3_b = (-d5_v) / d4_v
                d4_b = d3_b * yt
                d5_b = d3_b * xt
                do g_i_ = 1, g_pmax_
                    g_w(g_i_) = d5_b * g_yt(g_i_) + d4_b * g_xt(g_i_) + d2_b * g_zt(g_i_)
                enddo
                w = d5_v
                if( w .gt. 1.0d0 )then
                    do g_i_ = 1, g_pmax_
                        g_w(g_i_) = 0.0d0
                    enddo
                    w = 1.0d0
                endif
                if( w .lt. (-1.0d0) )then
                    do g_i_ = 1, g_pmax_
                        g_w(g_i_) = 0.0d0
                    enddo
                    w = -1.0d0
                endif
                d2_v = acos(w)
                if( abs(w) .lt. 1.0d0 )then
                    d1_p = -1.0d0 / sqrt ((1.0d0-w)*(1.0d0+w))
                else
                    !               call ehufDO (14,w, d2_v, d1_p,
                    !     +g_ehfid,
                    !     +547)
                endif
                do g_i_ = 1, g_pmax_
                    g_d1_w(g_i_) = d1_p * g_w(g_i_)
                enddo
                d1_w = d2_v
                d2_v = abs(d1_w)
                if( d1_w .gt. 0.0d0 )then
                    d1_p = 1.0d0
                elseif( d1_w .lt. 0.0d0 )then
                    d1_p = -1.0d0
                else
                    !               call ehufDO (3,d1_w, d2_v, d1_p,
                    !     +g_ehfid,
                    !     +562)
                endif
                d4_b = 1.0d0 / pi * 180.0d0 * d1_p
                do g_i_ = 1, g_pmax_
                    g_theta(g_i_, k) = d4_b * g_d1_w(g_i_)
                enddo
                theta(k) = d2_v * 180.0d0 / pi
            enddo
        enddo
        ! Now we we will relabel theta such that
        ! H1-N-H2 : theta(3)
        ! H1-N-H3 : theta(2)
        ! H2-N-H3 : theta(1)
        do g_i_ = 1, g_pmax_
            g_xt(g_i_) = g_theta(g_i_, 1)
        enddo
        xt = theta(1)
        do g_i_ = 1, g_pmax_
            g_theta(g_i_, 1) = g_theta(g_i_, 3)
        enddo
        theta(1) = theta(3)
        do g_i_ = 1, g_pmax_
            g_theta(g_i_, 3) = g_xt(g_i_)
        enddo
        theta(3) = xt
        ! all theta should be +ve and <= 180
        !    do k=1,3
        !        if(theta(k) < 0.0d0) theta(k) = -theta(k)
        !        !if(theta(k) > 180.0d0) theta(k) = theta(k) - 180.0d0
        !    end do
        do g_i_ = 1, g_pmax_
            g_bow(g_i_,1) = g_r(g_i_, 1)
        enddo
        bow(1) = r(1)
        do g_i_ = 1, g_pmax_
            g_bow(g_i_,2) = g_r(g_i_, 2)
        enddo
        bow(2) = r(2)
        do g_i_ = 1, g_pmax_
            g_bow(g_i_,3) = g_r(g_i_, 3)
        enddo
        bow(3) = r(3)
        do g_i_ = 1, g_pmax_
            g_bow(g_i_,4) = g_theta(g_i_, 1)
        enddo
        bow(4) = theta(1)
        do g_i_ = 1, g_pmax_
            g_bow(g_i_,5) = g_theta(g_i_, 2)
        enddo
        bow(5) = theta(2)
        do g_i_ = 1, g_pmax_
            g_bow(g_i_,6) = g_theta(g_i_, 3)
        enddo
        bow(6) = theta(3)
        do g_i_ = 1, g_pmax_
            g_bow(g_i_,7) = g_beta(g_i_, 1)
        enddo
        bow(7) = beta(1)
    end subroutine tobowman

    !> @brief Evaluate g_distplane for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine g_distplane( a1, g_a1, a2, g_a2, a3, g_a3, d, g_d, x, g_&
            &x, y, g_y, z, g_z, dist, g_dist )
        implicit none
        double precision :: a1
        double precision :: a2
        double precision :: a3
        double precision :: d
        double precision :: x
        double precision :: y
        double precision :: z
        double precision :: dist
        ! This subroutine calculates the dist between the
        ! point (x,y,z) and the plane (a1x+a2y+a3z+d=0)
        ! dist = a1x + a2y + a3z + d
        !      sqrt(a**2+b**2+c**2)
        integer, parameter :: g_pmax_ = 12
        integer            :: g_i_
        double precision   :: d16_b, d15_b, d14_b, d13_b, d1_w, d3_p, d2_v, d&
            &1_p, d4_v, d10_b
        double precision :: d9_b, d7_v, d16_v, d15_v, d2_b, d2_p, d4_b, g_d&
            &1_w(g_pmax_), g_a1(g_pmax_), g_a2(g_pmax_)
        double precision :: g_a3(g_pmax_), g_dist(g_pmax_), g_x(g_pmax_), g&
            &_y(g_pmax_), g_z(g_pmax_), g_d(g_pmax_)
        integer :: g_ehfid
        !        call ehsfid(g_ehfid, 'distplane','g_toBowman.f')
        g_ehfid = 0
        d2_v    = a1 * a1
        d3_p    = 2.0d0 * a1
        d4_v    = a2 * a2
        d2_p    = 2.0d0 * a2
        d7_v    = a3 * a3
        d1_p    = 2.0d0 * a3
        do g_i_ = 1, g_pmax_
            g_d1_w(g_i_) = d1_p * g_a3(g_i_) + d2_p * g_a2(g_i_) + d3_p * g_a1(g_i_)
        enddo
        d1_w  = d2_v + d4_v + d7_v
        d15_v = sqrt(d1_w)
        if( d1_w .gt. 0.0d0 )then
            d1_p = 1.0d0 / (2.0d0 *  d15_v)
        else
            !           call ehufDO (9,d1_w, d15_v, d1_p,
            !     +g_ehfid,
            !     +695)
        endif
        d16_v = (a1 * x + a2 * y + a3 * z + d) / d15_v
        d2_b  = 1.0d0 / d15_v
        d4_b  = (-d16_v) / d15_v * d1_p
        d9_b  = d2_b * z
        d10_b = d2_b * a3
        d13_b = d2_b * y
        d14_b = d2_b * a2
        d15_b = d2_b * x
        d16_b = d2_b * a1
        do g_i_ = 1, g_pmax_
            g_dist(g_i_) = d4_b * g_d1_w(g_i_) + d2_b * g_d(g_i_) + d10_b * g_z(g_i_) + d9_b * g_a3(g_i_) + d14_b * &
                & g_y(g_i_) + d13_b * g_a2 (g_i_) + d16_b * g_x(g_i_) + d15_b * g_a1(g_i_)
        enddo
        dist = d16_v
    end subroutine g_distplane

    !> @brief Evaluate g_det3 for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine g_det3( a, g_a, ans, g_ans )
        implicit none
        double precision :: a(3,3)
        double precision :: ans
        ! Compute determinant of a 3x3 matrix
        double precision   :: x, y, z
        integer, parameter :: g_pmax_ = 12
        integer            :: g_i_
        double precision   :: d9_b, d8_b, d7_b, d6_b, d8_v, g_x(g_pmax_), g_a &
            & (g_pmax_, 3, 3), g_y(g_pmax_), g_z(g_pmax_), g_ans(g_pmax_)
        integer :: g_ehfid
        !        call ehsfid(g_ehfid, 'det3','g_toBowman.f')
        g_ehfid = 0
        d8_v    = a(2,2) * a(3,3) - a(2,3) * a(3,2)
        d6_b    = (-a(1,1)) * a(3,2)
        d7_b    = (-a(1,1)) * a(2,3)
        d8_b    = a(1,1) * a(3,3)
        d9_b    = a(1,1) * a(2,2)
        do g_i_ = 1, g_pmax_
            g_x(g_i_) = d7_b * g_a(g_i_, 3, 2) + d6_b * g_a(g_i_, 2, 3) + d9_b * g_a(g_i_, 3, 3) + d8_b * g_a(g_i_, 2, &
                & 2) + d8_v * g_a(g_i_, 1, 1)
        enddo
        x    = a(1,1) * d8_v
        d8_v = a(2,1) * a(3,3) - a(2,3) * a(3,1)
        d6_b = (-a(1,2)) * a(3,1)
        d7_b = (-a(1,2)) * a(2,3)
        d8_b = a(1,2) * a(3,3)
        d9_b = a(1,2) * a(2,1)
        do g_i_ = 1, g_pmax_
            g_y(g_i_) = d7_b * g_a(g_i_, 3, 1) + d6_b * g_a(g_i_, 2, 3) + d9_b * g_a(g_i_, 3, 3) + d8_b * g_a(g_i_, 2, &
                & 1) + d8_v * g_a(g_i_, 1, 2)
        enddo
        y    = a(1,2) * d8_v
        d8_v = a(2,1) * a(3,2) - a(2,2) * a(3,1)
        d6_b = (-a(1,3)) * a(3,1)
        d7_b = (-a(1,3)) * a(2,2)
        d8_b = a(1,3) * a(3,2)
        d9_b = a(1,3) * a(2,1)
        do g_i_ = 1, g_pmax_
            g_z(g_i_) = d7_b * g_a(g_i_, 3, 1) + d6_b * g_a(g_i_, 2, 2) + d9_b * g_a(g_i_, 3, 2) + d8_b * g_a(g_i_, 2, &
                & 1) + d8_v * g_a(g_i_, 1, 3)
        enddo
        z = a(1,3) * d8_v
        do g_i_ = 1, g_pmax_
            g_ans(g_i_) = g_z(g_i_) + (-g_y(g_i_)) + g_x(g_i_)
        enddo
        ans = x - y + z
    end subroutine g_det3

    !> @brief Evaluate carttoang for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine carttoang( cart, ans, gangint )
        ! This function converts cartesian coordinates
        ! to internal coordinates (dist and bond angles)
        ! Cartesian coordinate format: x1,y1,z1,x2,y2,z2....xn,yn,zn
        ! N  : cart(1-3)
        ! H1 : cart(4-6)
        ! H2 : cart(7-9)
        ! H3 : cart(10-11)
        ! Dist in Ang. Angles in degrees
        ! Internal coordinate format:
        ! H1-N-H2 : theta3
        ! H1-N-H3 : theta2
        ! H2-N-H3 : theta1
        ![r1,r2,r3,theta1,theta2,theta3]
        implicit none
        ! I/O variables
        integer, parameter          :: natom = 4
        double precision, parameter :: pi    = 3.1415926535897930d0
        double precision            :: cart(3*natom)
        double precision            :: gangint(12,7)
        ! Local variables
        integer          :: i, j, k
        double precision :: r(3), t(3)
        double precision :: dt(3)
        double precision :: theta(3)
        double precision :: carth(3,3)
        double precision :: ans(7)
        ! to save the first three bond distances. If they are calculated in any
        ! subroutines, avoid to calculate them again
        ! initialize gangint
        do i=1, 7
            do j=1, 12
                gangint(j,i) = 0.d0
            enddo
        enddo
        ! Put N at the origin and
        ! and get N-H bond vectors
        k = 0
        do i=2, natom
            k          = k+1
            carth(k,1) = cart(3*i-2)-cart(1)
            carth(k,2) = cart(3*i-1)-cart(2)
            carth(k,3) = cart(3*i  )-cart(3)
        enddo
        ! calculate bond distances
        do i=1, 3
            r(i) = dsqrt(carth(i,1)*carth(i,1) + carth(i,2)*carth(i,2) + carth(i,3)*carth(i,3) )
        enddo
        ! dr/dx,cartH(i) is a function of the cartesian coord. of the i+1 atom
        do i=1, 3
            do j=1, 3
                k = 3*i + j
                !         dr(i)/dx(2),dr(i)/dx(3),dr(i)/dx(4)
                gangint(k,i) = carth(i,j)/r(i)
                !         dr(i)/dx(1)
                gangint(j,i) = -gangint(k,i)
            enddo
        enddo
        ! calculate bond angles
        t(1) = (carth(2,1)*carth(3,1)+carth(2,2)*carth(3,2)+ carth(2,3)*carth(3,3))/(r(2)*r(3))
        t(2) = (carth(1,1)*carth(3,1)+carth(1,2)*carth(3,2)+ carth(1,3)*carth(3,3))/(r(1)*r(3))
        t(3) = (carth(1,1)*carth(2,1)+carth(1,2)*carth(2,2)+ carth(1,3)*carth(2,3))/(r(1)*r(2))
        do i=1, 3
            !    Due to numerical errors, t maybe larger than 1.0 by a tiny value
            !        by definition 0 <= dacos(t) <= PI
            if( t(i) .gt. 1.d0 )then
                t(i) = 1.d0
            elseif( t(i) .lt. -1.d0 )then
                t(i) = -1.d0
            endif
            theta(i) = dacos( t(i) )
            !       convert to degrees
            theta(i) = theta(i)*180.0d0/pi
        enddo
        ! return answer
        do i=1, 3
            ans(i)   = r(i)
            ans(3+i) = theta(i)
        enddo
        !   Derivatives of the three angles
        !    dtheta1/dt1  theta = acos(t)
        do i=1, 3
            if( dabs(t(i)) .ge. 1.d0 )then
                !         write(*,*) "Abnormal structure: N-H2-H3 is linear"
                !         stop
                dt(i) = 0.d0
            else
                dt(i) = -(1.d0/dsqrt( 1.d0-t(i)*t(i) ))*180.d0/pi
            endif
        enddo
        !    dtheta1/dx3,dtheta1/dx4,dtheta1/dx1. dtheta1/dx2=0
        !    theta1 is a function of atom 1,3,4
        do j=1, 3
            k            = 6+j
            gangint(k,4) = carth(3,j)/(r(2)*r(3)) - gangint(k,2)*t(1)/r(2)
            gangint(k,4) = dt(1) * gangint(k,4)
            k            = 9+j
            gangint(k,4) = carth(2,j)/(r(2)*r(3)) - gangint(k,3)*t(1)/r(3)
            gangint(k,4) = dt(1) * gangint(k,4)
            gangint(j,4) = -(carth(2,j)+carth(3,j))/(r(2)*r(3)) - ( r(3)*gangint(j,2)+r(2)*gangint(j,3) ) *  t(1)/ &
                & (r(2)* r(3))
            gangint(j,4) = dt(1) * gangint(j,4)
        enddo
        !    dtheta2/dt2  theta = acos(t)
        !    dtheta2/dx2,dtheta2/dx4,dtheta2/dx1. dtheta2/dx3=0
        !    theta2 is a function of atom 1,2,4
        do j=1, 3
            k            = 3+j
            gangint(k,5) = carth(3,j)/(r(1)*r(3)) - gangint(k,1)*t(2)/r(1)
            gangint(k,5) = dt(2) * gangint(k,5)
            k            = 9+j
            gangint(k,5) = carth(1,j)/(r(1)*r(3)) - gangint(k,3)*t(2)/r(3)
            gangint(k,5) = dt(2) * gangint(k,5)
            gangint(j,5) = -(carth(1,j)+carth(3,j))/(r(1)*r(3)) - ( r(1)*gangint(j,3)+r(3)*gangint(j,1) ) *  t(2)/ &
                & (r(1)* r(3))
            gangint(j,5) = dt(2) * gangint(j,5)
        enddo
        !    dtheta3/dt3  theta = acos(t)
        !    dtheta3/dx2,dtheta3/dx3,dtheta3/dx1. dtheta3/dx4=0
        !    theta3 is a function of atom 1,2,3
        do j=1, 3
            k            = 3+j
            gangint(k,6) = carth(2,j)/(r(1)*r(2)) - gangint(k,1)*t(3)/r(1)
            gangint(k,6) = dt(3) * gangint(k,6)
            k            = 6+j
            gangint(k,6) = carth(1,j)/(r(1)*r(2)) - gangint(k,2)*t(3)/r(2)
            gangint(k,6) = dt(3) * gangint(k,6)
            gangint(j,6) = -(carth(1,j)+carth(2,j))/(r(1)*r(2)) - ( r(1)*gangint(j,2)+r(2)*gangint(j,1) ) *  t(3)/ &
                & (r(1)* r(2))
            gangint(j,6) = dt(3) * gangint(j,6)
        enddo
    end subroutine carttoang

    !> @brief Evaluate cart_to_polar for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine cart_to_polar( x, ans )
        implicit none
        double precision            :: x(3)
        double precision            :: ans(3)
        double precision, parameter :: pi = 3.1415926535897930d0
        ! ans(1) : r
        ! ans(2) : phi (0-2PI)  !in degrees
        ! ans(3) : theta (0-PI) !in degrees
        double precision :: t
        ans(1) = dsqrt( x(1)*x(1) + x(2)*x(2) + x(3)*x(3) )
        if( ans(1) .ne. 0.d0 )then
            t      = dacos(x(3)/ans(1))      !angle is between [-PI,PI]
            ans(3) = dasin(dsin(t))          !angle is between [0,2PI]
            if( x(1) .ne. 0.0d0 )then
                ans(2) = datan(x(2)/x(1))
            else
                ans(2)                       = 0.50d0 * pi
                if( x(2) .eq. 0.0d0 ) ans(2) = 0.0d0
            endif
            ! Convert to degrees
            ans(2) = ans(2) * 180.0d0/pi
            ans(3) = ans(3) * 180.0d0/pi
        else
            ans(2) = 0.d0
            ans(3) = 0.d0
        endif
    end subroutine cart_to_polar

    !> @brief Evaluate polar_to_cart for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine polar_to_cart( p, ans )
        implicit none
        double precision            :: p(3)
        double precision            :: ans(3)
        double precision, parameter :: pi = 3.1415926535897930d0
        ! p(1) : r
        ! p(2) : phi (0-2PI)  !in degrees
        ! p(3) : theta (0-PI) !in degrees
        double precision :: phi, theta
        double precision :: sinth
        phi    = p(2) * pi/180.0d0
        theta  = p(3) * pi/180.0d0
        sinth  = dsin(theta)
        ans(1) = p(1) * dcos(phi) * sinth
        ans(2) = p(1) * dsin(phi) * sinth
        ans(3) = p(1) * dcos(theta)
    end subroutine polar_to_cart

    !> @brief Evaluate fromBowman for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine frombowman( bow, coord )
        implicit none
        double precision            :: bow(7)
        double precision            :: coord(12)
        double precision, parameter :: pi = 3.1415926535897930d0
        ! This function returns the cartesian
        ! coordinates from bowman coordinates for 4 atom systems ONLY
        ! All angles are in DEGREES
        ! Input:
        ! bow(1) : r1
        ! bow(2) : r2
        ! bow(3) : r3
        ! bow(4) : theta1
        ! bow(5) : theta2
        ! bow(5) : theta3
        ! bow(7) : beta
        ! Output:
        ! coord(1:3)   : N   x,y,z
        ! coord(4:6)   : H1  x,y,z
        ! coord(7:9)   : H2  x,y,z
        ! coord(10:12) : H3  x,y,z
        double precision :: r1, r2, r3
        double precision :: theta1, theta2, theta3
        double precision :: phi1, phi2, phi3
        double precision :: beta
        double precision :: ph1(3)
        double precision :: ph2(3)
        double precision :: ph3(3)
        double precision :: ans(3)
        integer          :: i, j
        double precision :: small
        small = 1.d-3
        ! define stuff
        r1     = bow(1)
        r2     = bow(2)
        r3     = bow(3)
        theta1 = bow(4)
        theta2 = bow(5)
        theta3 = bow(6)
        beta   = bow(7)
        ! z-axis is the trisector vector
        ! G1, G2, G3 are the projection of H1, H2 and H3
        ! in the xy plane
        !  G2
        !   \
        !    N--G1-------> x-axis
        !   /
        !  G3
        ! Angle definitions:
        ! G1-N-G2 : theta3
        ! G1-N-G3 : theta2
        ! G2-N-G3 : theta1
        ! define phi
        phi1 = 0.0d0
        phi2 = theta3
        ! the following line (Shikha) is correct only when the sum of
        ! the three theta is 360.0
        !      phi3 = theta1+theta3
        if( dabs(360.d0-(theta1+theta2+theta3)) .lt. small )then
            phi3 = theta1+theta3
        elseif( dabs(theta1-theta2-theta3) .lt. small )then
            phi3 = 360.d0-theta2
        elseif( dabs(theta2-theta1-theta3) .lt. small )then
            phi3 = theta2
        elseif( dabs(theta3-theta1-theta2) .lt. small )then
            phi3 = theta2
        else
            write(6, *) "Is this geometry possible?"
            write(6, *) bow
            stop
        endif
        ! vector for H1 in polar
        ph1(1) = r1
        ph1(2) = phi1     !in degrees
        ph1(3) = beta     !in degrees
        ! vector for H2 in polar
        ph2(1) = r2
        ph2(2) = phi2   !G2 makes theta2 with x-axis
        ph2(3) = beta
        ! vector for H3 in polar
        ph3(1) = r3
        ph3(2) = phi3   !G3 makes theta2+theta3 with x-axis
        ph3(3) = beta
        ! convert to cartesian vectors
        ! N is at the origin
        coord(1) = 0.0d0
        coord(2) = 0.0d0
        coord(3) = 0.0d0
        ! H1
        call polar_to_cart(ph1, ans)
        do i=1, 3
            coord(3+i) = ans(i)
        enddo
        ! H2
        call polar_to_cart(ph2, ans)
        do i=1, 3
            coord(6+i) = ans(i)
        enddo
        ! H3
        call polar_to_cart(ph3, ans)
        do i=1, 3
            coord(9+i) = ans(i)
        enddo
    end subroutine frombowman

    !> @brief Evaluate angle_to_cart for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine angle_to_cart( angint, coord )
        ! Convert Valence Internal Coordinates (VIC)
        ! to cartesian coordinates. for 4 atom systems only
        ! Cartesian coordinate format:
        ! N  : cart(1:3)
        ! H1 : cart(4:6)
        ! H2 : cart(7:9)
        ! H3 : cart(10:12)
        ! Dist in Ang. Angles in degrees
        ! Internal coordinate format:
        ! H1-N-H2 : theta3
        ! H1-N-H3 : theta2
        ! H2-N-H3 : theta1
        ![r1,r2,r3,theta1,theta2,theta3]
        implicit none
        ! I/O variables
        double precision :: angint(7)
        double precision :: coord(12)
        ! Local variables
        double precision :: r1
        double precision :: r2
        double precision :: r3
        double precision :: theta1
        double precision :: theta2
        double precision :: theta3
        double precision :: uh3(3)
        ! Misc
        double precision, parameter :: degtopi = 3.1415926535897930d0/180.d0
        double precision            :: x
        double precision            :: costh1, costh2, costh3, sinth3
        double precision, parameter :: small = 1.d-10
        integer                     :: i
        ! initialize
        r1     = angint(1)
        r2     = angint(2)
        r3     = angint(3)
        theta1 = angint(4)*degtopi
        theta2 = angint(5)*degtopi
        theta3 = angint(6)*degtopi
        do i = 1, 12
            coord(i) = 0.0d0
        enddo
        costh2 = dcos(theta2)
        costh3 = dcos(theta3)
        sinth3 = dsin(theta3)
        ! N is at origin
        ! H1 is along x
        coord(4) = r1
        ! H2 is in the xy plane
        coord(7) = r2*costh3
        coord(8) = r2*sinth3
        ! H3 is in the xyz p
        ! If H2 along x: sin(theta3) = 0), put H3 on xz plane
        ! theta3 = 0   theta1=theta2
        ! theta3 = PI  theta1=PI-theta2
        ! unit vector along the N-H3 bond (uH3)
        uh3(1) = costh2
        if( dabs(sinth3) .lt. small )then
            uh3(2) = 0.d0
            uh3(3) = dsin(theta2)
        else
            costh1 = dcos(theta1)
            uh3(2) = ( costh1 - costh2*costh3 )/sinth3
            x      = uh3(1)**2 + uh3(2)**2
            if( x .le. 1.d0 )then
                uh3(3) = dsqrt(1.0d0-x)
                !  Due to numerical errors, x may be just marginally larger than 1.d0
                !  (planar structure)
            elseif( (x-1.d0) .le. small )then
                uh3(3) = 0.d0
            else
                write(*, *)'ERROR: Set of three bond angles is unphysical'
                write(*, *) r1, r2, r3, theta1/degtopi, &
                    & theta2/degtopi, theta3/degtopi
                stop
            endif
        endif
        coord(10) = r3 * uh3(1)
        coord(11) = r3 * uh3(2)
        coord(12) = r3 * uh3(3)
    end subroutine angle_to_cart

    !> @brief Evaluate rarray for the NH3_GEN_FFW2_2007_DPEM potential.
    subroutine rarray( natoms, x, mx, r, dr, mr )
        ! calculate the rij array N(N-1)/2 elements
        ! x:  x1,y1,z1,...,xn,yn,zn
        ! dr(i,j): derivatives of jth r to i th cart coord
        implicit none
        integer :: natoms
        ! dimension of x and r, must be larger than 3*natoms and
        ! natoms*(natoms-1)/2 respectively
        integer          :: mx, mr
        double precision :: x(mx), r(mr), dr(mx,mr)
        integer          :: i, j, k, ii, ij, nr
        double precision :: dx, dy, dz
        ! to save the first three bond distances. If they are calculated in any
        ! subroutines, avoid to calculate them again
        nr = natoms*(natoms-1)/2
        do i=1, nr
            do j = 1, 3*natoms
                dr(j,i) = 0.d0
            enddo
        enddo
        k = 0
        do i = 1, natoms
            ii = 3*(i-1)
            do j = i+1, natoms
                k  = k+1
                ij = 3*(j-1)
                dx = x(ij+1) - x(ii+1)
                dy = x(ij+2) - x(ii+2)
                dz = x(ij+3) - x(ii+3)
                ! Since N is always the first atom, r(1) to r(3) are always the first
                ! three NH distances
                r(k)       = dsqrt(dx*dx+dy*dy+dz*dz)
                dr(ij+1,k) = dx/r(k)
                dr(ij+2,k) = dy/r(k)
                dr(ij+3,k) = dz/r(k)
                dr(ii+1,k) = -dr(ij+1,k)
                dr(ii+2,k) = -dr(ij+2,k)
                dr(ii+3,k) = -dr(ij+3,k)
            enddo
        enddo
    end subroutine rarray
end module chempotpy_nh3_gen_ffw2_2007_dpem_pes

!> @brief Expose NH3_GEN_FFW2_2007 through the Python-facing DPEM entry point.
subroutine dpem( x, igrad, u, ug )
    use chempotpy_nh3_gen_ffw2_2007_dpem_pes, only: evaluate
    implicit none
    integer, parameter            :: natoms  = 4
    integer, parameter            :: nstates = 2
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: u(nstates,nstates)
    double precision, intent(out) :: ug(nstates,nstates,natoms,3)
    call evaluate(x, igrad, u, ug)
end subroutine dpem
