!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the H2O2_GEN_1998 potential.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_legacy_units, only: legacy_hartree_to_electronvolt
    implicit none
    ! number of electronic state
    integer, parameter            :: nstates = 1
    integer, parameter            :: natoms  = 4
    integer, intent(in)           :: igrad
    double precision, intent(in)  :: x(natoms,3)
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    integer                       :: istate
    double precision              :: tx(12), dx(12), v
    integer                       :: iatom, idir, i, j
    logical, save                 :: first_time_data=.true.
    !initialize
    v = 0.d0; g = 0.d0; d = 0.d0

    do iatom=1, natoms
        do idir=1, 3
            j     = (iatom-1)*3+idir
            tx(j) = x(iatom,idir)
        end do
    end do

    call surf(v, tx, dx, 12)

    v  = v*legacy_hartree_to_electronvolt
    dx = dx*legacy_hartree_to_electronvolt
    !dx=dx*51.422067

    if (igrad==0) then
        do istate=1, nstates
            p(istate) = v
        end do
    else if (igrad==1) then
        do istate=1, nstates
            p(istate) = v
        end do
        do iatom=1, natoms
            do idir=1, 3
                j               = (iatom-1)*3+idir
                g(1,iatom,idir) = dx(j)
            end do
        end do
    else if (igrad==2) then
        write (*, *) 'Only energy and gradient are available'
    end if

end subroutine pes

!
!   System:                        H2O2
!   Functional Form:
!   Common name:                   Koput-Carter-Handy
!   Number of derivatives:         1
!   Number of bodies:              4
!   Number of electronic surfaces: 1
!   Interface:                     Section-2
!
!   References: Koput, J.; Carter, S.; Handy, N. C.;
!               J. Phys. Chem. A, 102 (1998) 6325.
!
!   Notes:
!
!
!> @brief Evaluate SURF for the H2O2_GEN_1998 potential.
subroutine surf( v, coord, dx, n3tm )
    implicit none
    double precision :: v, coord(n3tm), dx(n3tm)
    double precision :: dxx(n3tm), x(n3tm)
    integer          :: i, n3tm
    ! Changing units from bohr to angstroms
    !
    do 10 i = 1, 12
        x(i) = coord(i)*0.52918d0
    10 continue

    do i = 1, 12
        dx(i) = 0.0d0
    end do

    call potential(v, x, dx, dxx, 12)

    do i = 1, 12
        x(i)   = x(i)/0.52918d0
        dx(i)  = dx(i)*0.52918d0
        dxx(i) = dxx(i) * 0.52918d0 * 0.52918d0
    end do

end subroutine surf

!***********************************************************************

!> @brief Evaluate POTENTIAL for the H2O2_GEN_1998 potential.
subroutine potential( v, x, dx, dxx, n3tm )
    !
    !   System: H2O2, see Koput, J.; Carter, S.; Handy, N. C.;
    !              J. Phys. Chem. A, 102 (1998) 6325
    !              Prepared for Polyrate by VMA, October 2002.
    !
    !   Reference: Koput,Carter,&Handy, JPC A, 102, 6325 (1998).
    !
    !   All the information passed to and from the potential energy surface
    !   routine is in hartree atomic units.
    !
    !        This potential is written such that:
    !                       X(1)  - X(3)  : X, Y, Z for H1
    !                       X(4)  - X(6)  : X, Y, Z for O1
    !                       X(7)  - X(9)  : X, Y, Z for O2
    !                       X(10) - X(12) : X, Y, Z for H2
    !
    implicit none
    double precision :: coord(n3tm), dx(n3tm), dxx(n3tm), x(n3tm)
    double precision :: v, pass(6)
    double precision :: gx(12,12), gv(12,12)
    double precision :: gxx(12,12,12), gvv(12,12,12)
    double precision :: temp, temp1, temp2, temp3, temp4, temp5, &
        & temp6, temp7, temp8, temp9, temp10, temp11, temp12, &
        & temp13, temp14, temp15, temp16, temp17, temp18, &
        & temp19, temp20, temp21, temp22, temp23, temp24, &
        & temp25, temp26, temp27, temp28, temp29, temp30, &
        & temp31, temp32, temp33, temp34, temp35, temp36, &
        & temp37, temp38, temp39, temp40, temp41, temp42, &
        & temp43, temp44, temp45, temp46, temp47, temp48, &
        & temp49, temp50, temp51, temp52, temp53, temp54, &
        & temp55, temp56, temp57, temp58, temp59, temp60, &
        & temp61, temp62, temp63, temp64, temp65, temp66
    double precision :: tempx, temp1x, temp2x, temp3x, temp4x, temp5x, &
        & temp6x, temp7x, temp8x, temp9x, temp10x, temp11x, temp12x, &
        & temp13x, temp14x, temp15x, temp16x, temp17x, temp18x, &
        & temp19x, temp20x, temp21x, temp22x, temp23x, temp24x, &
        & temp25x, temp26x, temp27x, temp28x, temp29x, temp30x, &
        & temp31x, temp32x, temp33x, temp34x, temp35x, temp36x, &
        & temp37x, temp38x, temp39x, temp40x, temp41x, temp42x, &
        & temp43x, temp44x, temp45x, temp46x, temp47x, temp48x, &
        & temp49x, temp50x, temp51x, temp52x, temp53x, temp54x, &
        & temp55x, temp56x, temp57x, temp58x, temp59x, temp60x, &
        & temp61x, temp62x, temp63x, temp64x, temp65x, temp66x
    integer                     :: n3tm, i, j, ii
    double precision, parameter :: rohlow  = 1.0d0
    double precision, parameter :: rohhigh = 4.0d0
    double precision, parameter :: roolow  = 1.8d0
    double precision, parameter :: roohigh = 4.0d0
    double precision, parameter :: rhhlow  = 0.7d0
    !
    ! Call Koput, Carter, and Handy PES
    !
    call pot(v, x)

    !
    !  Initialize derivatives
    !
    do j=1, n3tm
        dx(j) = 0.0d0
    end do

    gx(1,1)   = 1.0d0
    gx(2,2)   = 1.0d0
    gx(3,3)   = 1.0d0
    gx(4,4)   = 1.0d0
    gx(5,5)   = 1.0d0
    gx(6,6)   = 1.0d0
    gx(7,7)   = 1.0d0
    gx(8,8)   = 1.0d0
    gx(9,9)   = 1.0d0
    gx(10,10) = 1.0d0
    gx(11,11) = 1.0d0
    gx(12,12) = 1.0d0

    gx(1,2)  = 0.0d0
    gx(1,3)  = 0.0d0
    gx(1,4)  = 0.0d0
    gx(1,5)  = 0.0d0
    gx(1,6)  = 0.0d0
    gx(1,7)  = 0.0d0
    gx(1,8)  = 0.0d0
    gx(1,9)  = 0.0d0
    gx(1,10) = 0.0d0
    gx(1,11) = 0.0d0
    gx(1,12) = 0.0d0

    gx(2,1)  = 0.0d0
    gx(2,3)  = 0.0d0
    gx(2,4)  = 0.0d0
    gx(2,5)  = 0.0d0
    gx(2,6)  = 0.0d0
    gx(2,7)  = 0.0d0
    gx(2,8)  = 0.0d0
    gx(2,9)  = 0.0d0
    gx(2,10) = 0.0d0
    gx(2,11) = 0.0d0
    gx(2,12) = 0.0d0

    gx(3,1)  = 0.0d0
    gx(3,2)  = 0.0d0
    gx(3,4)  = 0.0d0
    gx(3,5)  = 0.0d0
    gx(3,6)  = 0.0d0
    gx(3,7)  = 0.0d0
    gx(3,8)  = 0.0d0
    gx(3,9)  = 0.0d0
    gx(3,10) = 0.0d0
    gx(3,11) = 0.0d0
    gx(3,12) = 0.0d0

    gx(4,1)  = 0.0d0
    gx(4,2)  = 0.0d0
    gx(4,3)  = 0.0d0
    gx(4,5)  = 0.0d0
    gx(4,6)  = 0.0d0
    gx(4,7)  = 0.0d0
    gx(4,8)  = 0.0d0
    gx(4,9)  = 0.0d0
    gx(4,10) = 0.0d0
    gx(4,11) = 0.0d0
    gx(4,12) = 0.0d0

    gx(5,1)  = 0.0d0
    gx(5,2)  = 0.0d0
    gx(5,3)  = 0.0d0
    gx(5,4)  = 0.0d0
    gx(5,6)  = 0.0d0
    gx(5,7)  = 0.0d0
    gx(5,8)  = 0.0d0
    gx(5,9)  = 0.0d0
    gx(5,10) = 0.0d0
    gx(5,11) = 0.0d0
    gx(5,12) = 0.0d0

    gx(6,1)  = 0.0d0
    gx(6,2)  = 0.0d0
    gx(6,3)  = 0.0d0
    gx(6,4)  = 0.0d0
    gx(6,5)  = 0.0d0
    gx(6,7)  = 0.0d0
    gx(6,8)  = 0.0d0
    gx(6,9)  = 0.0d0
    gx(6,10) = 0.0d0
    gx(6,11) = 0.0d0
    gx(6,12) = 0.0d0

    gx(7,1)  = 0.0d0
    gx(7,2)  = 0.0d0
    gx(7,3)  = 0.0d0
    gx(7,4)  = 0.0d0
    gx(7,5)  = 0.0d0
    gx(7,6)  = 0.0d0
    gx(7,8)  = 0.0d0
    gx(7,9)  = 0.0d0
    gx(7,10) = 0.0d0
    gx(7,11) = 0.0d0
    gx(7,12) = 0.0d0

    gx(8,1)  = 0.0d0
    gx(8,2)  = 0.0d0
    gx(8,3)  = 0.0d0
    gx(8,4)  = 0.0d0
    gx(8,5)  = 0.0d0
    gx(8,6)  = 0.0d0
    gx(8,7)  = 0.0d0
    gx(8,9)  = 0.0d0
    gx(8,10) = 0.0d0
    gx(8,11) = 0.0d0
    gx(8,12) = 0.0d0

    gx(9,1)  = 0.0d0
    gx(9,2)  = 0.0d0
    gx(9,3)  = 0.0d0
    gx(9,4)  = 0.0d0
    gx(9,5)  = 0.0d0
    gx(9,6)  = 0.0d0
    gx(9,7)  = 0.0d0
    gx(9,8)  = 0.0d0
    gx(9,10) = 0.0d0
    gx(9,11) = 0.0d0
    gx(9,12) = 0.0d0

    gx(10,1)  = 0.0d0
    gx(10,2)  = 0.0d0
    gx(10,3)  = 0.0d0
    gx(10,4)  = 0.0d0
    gx(10,5)  = 0.0d0
    gx(10,6)  = 0.0d0
    gx(10,7)  = 0.0d0
    gx(10,8)  = 0.0d0
    gx(10,9)  = 0.0d0
    gx(10,11) = 0.0d0
    gx(10,12) = 0.0d0

    gx(11,1)  = 0.0d0
    gx(11,2)  = 0.0d0
    gx(11,3)  = 0.0d0
    gx(11,4)  = 0.0d0
    gx(11,5)  = 0.0d0
    gx(11,6)  = 0.0d0
    gx(11,7)  = 0.0d0
    gx(11,8)  = 0.0d0
    gx(11,9)  = 0.0d0
    gx(11,10) = 0.0d0
    gx(11,12) = 0.0d0

    gx(12,1)  = 0.0d0
    gx(12,2)  = 0.0d0
    gx(12,3)  = 0.0d0
    gx(12,4)  = 0.0d0
    gx(12,5)  = 0.0d0
    gx(12,6)  = 0.0d0
    gx(12,7)  = 0.0d0
    gx(12,8)  = 0.0d0
    gx(12,9)  = 0.0d0
    gx(12,10) = 0.0d0
    gx(12,11) = 0.0d0

    !
    !  Calculate first derivatives using a analytic gradients
    !

    call g_surf(12, v, gv, 1, x, gx, 12)

    !  TRANSPOSE GV

    tempx   = gv(2,1)
    gv(2,1) = gv(1,2)
    gv(1,2) = tempx

    temp2x  = gv(3,1)
    gv(3,1) = gv(1,3)
    gv(1,3) = temp2x

    temp3x  = gv(4,1)
    gv(4,1) = gv(1,4)
    gv(1,4) = temp3x

    temp4x  = gv(5,1)
    gv(5,1) = gv(1,5)
    gv(1,5) = temp4x

    temp5x  = gv(6,1)
    gv(6,1) = gv(1,6)
    gv(1,6) = temp5x

    temp6x  = gv(7,1)
    gv(7,1) = gv(1,7)
    gv(1,7) = temp6x

    temp7x  = gv(8,1)
    gv(8,1) = gv(1,8)
    gv(1,8) = temp7x

    temp8x  = gv(9,1)
    gv(9,1) = gv(1,9)
    gv(1,9) = temp8x

    temp9x   = gv(10,1)
    gv(10,1) = gv(1,10)
    gv(1,10) = temp9x

    temp10x  = gv(11,1)
    gv(11,1) = gv(1,11)
    gv(1,11) = temp10x

    temp11x  = gv(12,1)
    gv(12,1) = gv(1,12)
    gv(1,12) = temp11x

    temp12x = gv(2,3)
    gv(2,3) = gv(3,2)
    gv(3,2) = temp12x

    temp13x = gv(2,4)
    gv(2,4) = gv(4,2)
    gv(4,2) = temp13x

    temp14x = gv(2,5)
    gv(2,5) = gv(5,2)
    gv(5,2) = temp14x

    temp15x = gv(2,6)
    gv(2,6) = gv(6,2)
    gv(6,2) = temp15x

    temp16x = gv(2,7)
    gv(2,7) = gv(7,2)
    gv(7,2) = temp16x

    temp17x = gv(2,8)
    gv(2,8) = gv(8,2)
    gv(8,2) = temp17x

    temp18x = gv(2,9)
    gv(2,9) = gv(9,2)
    gv(9,2) = temp18x

    temp19x  = gv(2,10)
    gv(2,10) = gv(10,2)
    gv(10,2) = temp19x

    temp20x  = gv(2,11)
    gv(2,11) = gv(11,2)
    gv(11,2) = temp20x

    temp21x  = gv(2,12)
    gv(2,12) = gv(12,2)
    gv(12,2) = temp21x

    temp22x = gv(3,4)
    gv(3,4) = gv(4,3)
    gv(4,3) = temp22x

    temp23x = gv(3,5)
    gv(3,5) = gv(5,3)
    gv(5,3) = temp23x

    temp24x = gv(3,6)
    gv(3,6) = gv(6,3)
    gv(6,3) = temp24x

    temp25x = gv(3,7)
    gv(3,7) = gv(7,3)
    gv(7,3) = temp25x

    temp26x = gv(3,8)
    gv(3,8) = gv(8,3)
    gv(8,3) = temp26x

    temp27x = gv(3,9)
    gv(3,9) = gv(9,3)
    gv(9,3) = temp27x

    temp28x  = gv(3,10)
    gv(3,10) = gv(10,3)
    gv(10,3) = temp28x

    temp29x  = gv(3,11)
    gv(3,11) = gv(11,3)
    gv(11,3) = temp29x

    temp30x  = gv(3,12)
    gv(3,12) = gv(12,3)
    gv(12,3) = temp30x

    temp31x = gv(4,5)
    gv(4,5) = gv(5,4)
    gv(5,4) = temp31x

    temp32x = gv(4,6)
    gv(4,6) = gv(6,4)
    gv(6,4) = temp32x

    temp33x = gv(4,7)
    gv(4,7) = gv(7,4)
    gv(7,4) = temp33x

    temp34x = gv(4,8)
    gv(4,8) = gv(8,4)
    gv(8,4) = temp34x

    temp35x = gv(4,9)
    gv(4,9) = gv(9,4)
    gv(9,4) = temp35x

    temp36x  = gv(4,10)
    gv(4,10) = gv(10,4)
    gv(10,4) = temp36x

    temp37x  = gv(4,11)
    gv(4,11) = gv(11,4)
    gv(11,4) = temp37x

    temp38x  = gv(4,12)
    gv(4,12) = gv(12,4)
    gv(12,4) = temp38x

    temp39x = gv(5,6)
    gv(5,6) = gv(6,5)
    gv(6,5) = temp39x

    temp40x = gv(5,7)
    gv(5,7) = gv(7,5)
    gv(7,5) = temp40x

    temp41x = gv(5,8)
    gv(5,8) = gv(8,5)
    gv(8,5) = temp41x

    temp42x = gv(5,9)
    gv(5,9) = gv(9,5)
    gv(9,5) = temp42x

    temp43x  = gv(5,10)
    gv(5,10) = gv(10,5)
    gv(10,5) = temp43x

    temp44x  = gv(5,11)
    gv(5,11) = gv(11,5)
    gv(11,5) = temp44x

    temp45x  = gv(5,12)
    gv(5,12) = gv(12,5)
    gv(12,5) = temp45x

    temp46x = gv(6,7)
    gv(6,7) = gv(7,6)
    gv(7,6) = temp46x

    temp47x = gv(6,8)
    gv(6,8) = gv(8,6)
    gv(8,6) = temp47x

    temp48x = gv(6,9)
    gv(6,9) = gv(9,6)
    gv(9,6) = temp48x

    temp49x  = gv(6,10)
    gv(6,10) = gv(10,6)
    gv(10,6) = temp49x

    temp50x  = gv(6,11)
    gv(6,11) = gv(11,6)
    gv(11,6) = temp50x

    temp51x  = gv(6,12)
    gv(6,12) = gv(12,6)
    gv(12,6) = temp51x

    temp52x = gv(7,8)
    gv(7,8) = gv(8,7)
    gv(8,7) = temp52x

    temp53x = gv(7,9)
    gv(7,9) = gv(9,7)
    gv(9,7) = temp53x

    temp54x  = gv(7,10)
    gv(7,10) = gv(10,7)
    gv(10,7) = temp54x

    temp55x  = gv(7,11)
    gv(7,11) = gv(11,7)
    gv(11,7) = temp55x

    temp56x  = gv(7,12)
    gv(7,12) = gv(12,7)
    gv(12,7) = temp56x

    temp57x = gv(8,9)
    gv(8,9) = gv(9,8)
    gv(9,8) = temp57x

    temp58x  = gv(8,10)
    gv(8,10) = gv(10,8)
    gv(10,8) = temp58x

    temp59x  = gv(8,11)
    gv(8,11) = gv(11,8)
    gv(11,8) = temp59x

    temp60x  = gv(8,12)
    gv(8,12) = gv(12,8)
    gv(12,8) = temp60x

    temp61x  = gv(9,10)
    gv(9,10) = gv(10,9)
    gv(10,9) = temp61x

    temp62x  = gv(9,11)
    gv(9,11) = gv(11,9)
    gv(11,9) = temp62x

    temp63x  = gv(9,12)
    gv(9,12) = gv(12,9)
    gv(12,9) = temp63x

    temp64x   = gv(10,11)
    gv(10,11) = gv(11,10)
    gv(11,10) = temp64x

    temp65x   = gv(10,12)
    gv(10,12) = gv(12,10)
    gv(12,10) = temp65x

    temp66x   = gv(11,12)
    gv(11,12) = gv(12,11)
    gv(12,11) = temp66x

    do i = 1, 12
        dx(i) = gv(1,i)
    end do

    !
    !  Initialize second derivatives
    !
    do j=1, n3tm
        dxx(j) = 0.0d0
    end do

    do i = 1, 12

        gxx(1,2,i)  = 0.0d0
        gxx(1,3,i)  = 0.0d0
        gxx(1,4,i)  = 0.0d0
        gxx(1,5,i)  = 0.0d0
        gxx(1,6,i)  = 0.0d0
        gxx(1,7,i)  = 0.0d0
        gxx(1,8,i)  = 0.0d0
        gxx(1,9,i)  = 0.0d0
        gxx(1,10,i) = 0.0d0
        gxx(1,11,i) = 0.0d0
        gxx(1,12,i) = 0.0d0

        gxx(2,1,i)  = 0.0d0
        gxx(2,3,i)  = 0.0d0
        gxx(2,4,i)  = 0.0d0
        gxx(2,5,i)  = 0.0d0
        gxx(2,6,i)  = 0.0d0
        gxx(2,7,i)  = 0.0d0
        gxx(2,8,i)  = 0.0d0
        gxx(2,9,i)  = 0.0d0
        gxx(2,10,i) = 0.0d0
        gxx(2,11,i) = 0.0d0
        gxx(2,12,i) = 0.0d0

        gxx(3,1,i)  = 0.0d0
        gxx(3,2,i)  = 0.0d0
        gxx(3,4,i)  = 0.0d0
        gxx(3,5,i)  = 0.0d0
        gxx(3,6,i)  = 0.0d0
        gxx(3,7,i)  = 0.0d0
        gxx(3,8,i)  = 0.0d0
        gxx(3,9,i)  = 0.0d0
        gxx(3,10,i) = 0.0d0
        gxx(3,11,i) = 0.0d0
        gxx(3,12,i) = 0.0d0

        gxx(4,1,i)  = 0.0d0
        gxx(4,2,i)  = 0.0d0
        gxx(4,3,i)  = 0.0d0
        gxx(4,5,i)  = 0.0d0
        gxx(4,6,i)  = 0.0d0
        gxx(4,7,i)  = 0.0d0
        gxx(4,8,i)  = 0.0d0
        gxx(4,9,i)  = 0.0d0
        gxx(4,10,i) = 0.0d0
        gxx(4,11,i) = 0.0d0
        gxx(4,12,i) = 0.0d0

        gxx(5,1,i)  = 0.0d0
        gxx(5,2,i)  = 0.0d0
        gxx(5,3,i)  = 0.0d0
        gxx(5,4,i)  = 0.0d0
        gxx(5,6,i)  = 0.0d0
        gxx(5,7,i)  = 0.0d0
        gxx(5,8,i)  = 0.0d0
        gxx(5,9,i)  = 0.0d0
        gxx(5,10,i) = 0.0d0
        gxx(5,11,i) = 0.0d0
        gxx(5,12,i) = 0.0d0

        gxx(6,1,i)  = 0.0d0
        gxx(6,2,i)  = 0.0d0
        gxx(6,3,i)  = 0.0d0
        gxx(6,4,i)  = 0.0d0
        gxx(6,5,i)  = 0.0d0
        gxx(6,7,i)  = 0.0d0
        gxx(6,8,i)  = 0.0d0
        gxx(6,9,i)  = 0.0d0
        gxx(6,10,i) = 0.0d0
        gxx(6,11,i) = 0.0d0
        gxx(6,12,i) = 0.0d0

        gxx(7,1,i)  = 0.0d0
        gxx(7,2,i)  = 0.0d0
        gxx(7,3,i)  = 0.0d0
        gxx(7,4,i)  = 0.0d0
        gxx(7,5,i)  = 0.0d0
        gxx(7,6,i)  = 0.0d0
        gxx(7,8,i)  = 0.0d0
        gxx(7,9,i)  = 0.0d0
        gxx(7,10,i) = 0.0d0
        gxx(7,11,i) = 0.0d0
        gxx(7,12,i) = 0.0d0

        gxx(8,1,i)  = 0.0d0
        gxx(8,2,i)  = 0.0d0
        gxx(8,3,i)  = 0.0d0
        gxx(8,4,i)  = 0.0d0
        gxx(8,5,i)  = 0.0d0
        gxx(8,6,i)  = 0.0d0
        gxx(8,7,i)  = 0.0d0
        gxx(8,9,i)  = 0.0d0
        gxx(8,10,i) = 0.0d0
        gxx(8,11,i) = 0.0d0
        gxx(8,12,i) = 0.0d0

        gxx(9,1,i)  = 0.0d0
        gxx(9,2,i)  = 0.0d0
        gxx(9,3,i)  = 0.0d0
        gxx(9,4,i)  = 0.0d0
        gxx(9,5,i)  = 0.0d0
        gxx(9,6,i)  = 0.0d0
        gxx(9,7,i)  = 0.0d0
        gxx(9,8,i)  = 0.0d0
        gxx(9,10,i) = 0.0d0
        gxx(9,11,i) = 0.0d0
        gxx(9,12,i) = 0.0d0

        gxx(10,1,i)  = 0.0d0
        gxx(10,2,i)  = 0.0d0
        gxx(10,3,i)  = 0.0d0
        gxx(10,4,i)  = 0.0d0
        gxx(10,5,i)  = 0.0d0
        gxx(10,6,i)  = 0.0d0
        gxx(10,7,i)  = 0.0d0
        gxx(10,8,i)  = 0.0d0
        gxx(10,9,i)  = 0.0d0
        gxx(10,11,i) = 0.0d0
        gxx(10,12,i) = 0.0d0

        gxx(11,1,i)  = 0.0d0
        gxx(11,2,i)  = 0.0d0
        gxx(11,3,i)  = 0.0d0
        gxx(11,4,i)  = 0.0d0
        gxx(11,5,i)  = 0.0d0
        gxx(11,6,i)  = 0.0d0
        gxx(11,7,i)  = 0.0d0
        gxx(11,8,i)  = 0.0d0
        gxx(11,9,i)  = 0.0d0
        gxx(11,10,i) = 0.0d0
        gxx(11,12,i) = 0.0d0

        gxx(12,1,i)  = 0.0d0
        gxx(12,2,i)  = 0.0d0
        gxx(12,3,i)  = 0.0d0
        gxx(12,4,i)  = 0.0d0
        gxx(12,5,i)  = 0.0d0
        gxx(12,6,i)  = 0.0d0
        gxx(12,7,i)  = 0.0d0
        gxx(12,8,i)  = 0.0d0
        gxx(12,9,i)  = 0.0d0
        gxx(12,10,i) = 0.0d0
        gxx(12,11,i) = 0.0d0

    end do

    do i = 1, 12
        do j = 1, 12
            do ii = 1, 12
                gxx(i,j,ii) = 0.0d0
            end do
        end do
    end do

    gxx(1,1,1)    = 1.0d0
    gxx(2,2,2)    = 1.0d0
    gxx(3,3,3)    = 1.0d0
    gxx(4,4,4)    = 1.0d0
    gxx(5,5,5)    = 1.0d0
    gxx(6,6,6)    = 1.0d0
    gxx(7,7,7)    = 1.0d0
    gxx(8,8,8)    = 1.0d0
    gxx(9,9,9)    = 1.0d0
    gxx(10,10,10) = 1.0d0
    gxx(11,11,11) = 1.0d0
    gxx(12,12,12) = 1.0d0

    do i = 1, 12
        do j = 1, 12
            do ii = 1, 12
            end do
        end do
    end do

    do i = 1, 12
        do j = 1, 12
            gv(i,j) = 0.0d0
            gx(i,j) = 0.0d0
        end do
    end do

    !
    !  Calculate first derivatives using a analytic gradients
    !

    call g_analytic(12, 12, v, gv, gvv, 1, 1, x, gx, gxx, 12, 12)

    !  TRANSPOSE GVV

    do i = 1, 12
        dxx(i) = gvv(1,i,i)/0.529177d0
    end do

end subroutine potential

!***********************************************************************

!> @brief Evaluate POT for the H2O2_GEN_1998 potential.
subroutine pot( v, x )
    !
    !   System: H2O2, see Koput, J.; Carter, S.; Handy, N. C.;
    !              J. Phys. Chem. A, 102 (1998) 6325
    !              Prepared for Polyrate by VMA, October 2002.
    !
    !   Reference: Koput,Carter,&Handy, JPC A, 102, 6325 (1998).
    !
    !   All the information passed to and from the potential energy surface
    !   routine is in hartree atomic units.
    !
    !        This potential is written such that:
    !                       X(1)  - X(3)  : X, Y, Z for H1
    !                       X(4)  - X(6)  : X, Y, Z for O1
    !                       X(7)  - X(9)  : X, Y, Z for O2
    !                       X(10) - X(12) : X, Y, Z for H2
    !

    implicit none
    double precision            :: n, n1, n2
    double precision            :: v, pi
    double precision            :: coord(12), dx(12), x(12), pass(6)
    double precision            :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10
    double precision            :: c11, c12, c13, c14, c15, c16, c17, c18
    double precision            :: rab, rbc, rcd, rac, rbd, rad
    double precision            :: rh1o1, ro1o2, ro2h2
    double precision            :: d1, d2, d, ah1oo, aooh2, thooh
    double precision            :: xcomp1, ycomp1, zcomp1
    double precision            :: xcomp2, ycomp2, zcomp2
    double precision            :: xcomp3, ycomp3, zcomp3
    double precision            :: xcomp4, ycomp4, zcomp4
    double precision            :: x_normal1, y_normal1, z_normal1
    double precision            :: x_normal2, y_normal2, z_normal2
    double precision, parameter :: rohlow  = 1.0d0
    double precision, parameter :: rohhigh = 4.0d0
    double precision, parameter :: roolow  = 1.8d0
    double precision, parameter :: roohigh = 4.0d0
    double precision, parameter :: rhhlow  = 0.7d0
    !  Calculate energy at the configuration space point
    !  NOTE :: This is not necessarily a quadrature point
    ! Expecting input in cartesians
    pi = 4.d0*atan(1.d0)

    c1 = x(1) - x(4)
    c2 = x(2) - x(5)
    c3 = x(3) - x(6)
    c4 = x(4) - x(7)
    c5 = x(5) - x(8)
    c6 = x(6) - x(9)
    c7 = x(7) - x(10)
    c8 = x(8) - x(11)
    c9 = x(9) - x(12)

    c10 = x(1) - x(7)
    c11 = x(2) - x(8)
    c12 = x(3) - x(9)
    c13 = x(4) - x(10)
    c14 = x(5) - x(11)
    c15 = x(6) - x(12)
    c16 = x(1)-x(10)
    c17 = x(2)-x(11)
    c18 = x(3)-x(12)

    !
    ! H1 -- O1 -- O2 -- H2
    !   rH1O1 rO1O2 rO2H2
    !
    !  angle H1O1O2 = AH1OO
    !  angle O1O2H2 = AOOH2
    !
    !  TORSION HOOH = THOOH

    rab   = sqrt(c1**2+c2**2+c3**2)
    rbc   = sqrt(c4**2+c5**2+c6**2)
    rcd   = sqrt(c7**2+c8**2+c9**2)
    rac   = sqrt(c10**2+c11**2+c12**2)
    rbd   = sqrt(c13**2+c14**2+c15**2)
    rad   = sqrt(c16**2+c17**2+c18**2)
    rh1o1 = rab
    ro1o2 = rbc
    ro2h2 = rcd

    ! The potential can give weird values for unusual geometries
    ! replace with some high value in these unphysical regions and
    ! then exit
    !

    n1    = (-c4)*c1 + (-c5)*c2 + (-c6)*c3
    d1    = ro1o2*rh1o1
    ah1oo = acos(n1/d1)

    n2    = -(c7*c4 + c8*c5 + c9*c6)
    d2    = ro2h2*ro1o2
    aooh2 = acos(n2/d2)

    ! Calculate torsion angle HOOH

    xcomp1 = x(4) - x(1)
    ycomp1 = x(5) - x(2)
    zcomp1 = x(6) - x(3)

    xcomp2 = c10
    ycomp2 = c11
    zcomp2 = c12

    x_normal1 = ycomp1 * zcomp2 - zcomp1 * ycomp2
    y_normal1 = zcomp1 * xcomp2 - xcomp1 * zcomp2
    z_normal1 = xcomp1 * ycomp2 - ycomp1 * xcomp2

    xcomp3 = x(7) - x(4)
    ycomp3 = x(8) - x(5)
    zcomp3 = x(9) - x(6)

    xcomp4 = c13
    ycomp4 = c14
    zcomp4 = c15

    x_normal2 = ycomp3 * zcomp4 - zcomp3 * ycomp4
    y_normal2 = zcomp3 * xcomp4 - xcomp3 * zcomp4
    z_normal2 = xcomp3 * ycomp4 - ycomp3 * xcomp4

    n = x_normal1 * x_normal2 + y_normal1 * y_normal2 + z_normal1 * z_normal2

    d = (sqrt(x_normal1 * x_normal1 + y_normal1 * y_normal1 + z_normal1 * z_normal1)) * (sqrt(x_normal2 * x_normal2 + &
        & y_normal2 * y_normal2 + z_normal2 * z_normal2))

    thooh = acos(n/d)  ! in radians

    pass(1) = rh1o1
    pass(2) = ro1o2
    pass(3) = ro2h2
    pass(4) = ah1oo
    pass(5) = aooh2
    pass(6) = thooh

    !
    ! Call Koput, Carter, and Handy PES
    !

    call vibpot(pass, v)

end subroutine pot

!***********************************************************************

!> @brief Evaluate VIBPOT for the H2O2_GEN_1998 potential.
subroutine vibpot( coord, v )
    implicit none
    double precision :: v, pi
    double precision :: coord(6), q1p(0:4), q2p(0:4), q3p(0:4), q4p(0:4)
    double precision :: q5p(0:4), q6p(0:4)
    double precision :: ind1(152), ind2(152), ind3(152), ind4(152)
    double precision :: ind5(152), ind6(152), coeff(152)
    double precision :: rm_oo, rm_oh, angle_ooh, r1, r2, r3
    double precision :: a1, a2, q1, q2, q3, q4, q5, q6
    integer          :: i
    double precision :: zoe
    data zoe /3.995139982058d-3/
    data ind1 /0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, &
        & 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, &
        & 1, 2, 1, 0, 2, 0, 1, 0, 2, 0, 0, 0, 1, 1, 0, 1, &
        & 0, 0, 1, 1, 1, 0, 0, 4, 0, 0, 0, 2, 0, 0, 0, 2, &
        & 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 3, 1, 0, 3, 0, &
        & 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 2, 0, 2, 0, 0, 0, &
        & 1, 1, 2, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 2, 0, &
        & 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 3, 0, 0, 0, 2, &
        & 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, &
        & 1, 0, 0, 0, 1, 0, 0, 0/
    data ind2 /0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, &
        & 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, &
        & 2, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 0, 1, 0, 1, 0, &
        & 1, 0, 1, 1, 0, 1, 0, 0, 4, 0, 0, 0, 2, 0, 0, 0, &
        & 2, 0, 0, 1, 0, 0, 0, 3, 0, 0, 3, 1, 0, 1, 0, 3, &
        & 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 2, 0, 2, 0, 0, &
        & 1, 1, 0, 2, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 2, &
        & 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 3, 0, 0, 0, &
        & 2, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, &
        & 0, 1, 0, 0, 0, 1, 0, 0/
    data ind3 /0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, &
        & 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 2, 2, &
        & 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, &
        & 0, 0, 4, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, &
        & 0, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 1, 1, 1, &
        & 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, &
        & 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 1, 1, &
        & 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, &
        & 0, 0, 0, 0, 1, 0, 0, 0, 0/
    data ind4 /0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, &
        & 0, 1, 0, 0, 1, 1, 0, 0, 0, 3, 0, 0, 0, &
        & 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, 0, 0, &
        & 2, 0, 1, 1, 2, 0, 1, 0, 0, 1, 1, 1, 0, &
        & 1, 1, 0, 0, 0, 4, 0, 0, 0, 2, 0, 2, 0, &
        & 2, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 3, 0, &
        & 1, 0, 1, 3, 0, 1, 0, 0, 1, 1, 2, 0, 1, &
        & 0, 0, 1, 1, 2, 2, 0, 1, 1, 1, 2, 2, 1, &
        & 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, &
        & 0, 0, 1, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, &
        & 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, &
        & 0, 1, 0, 1, 0, 0, 0, 1, 0/
    data ind5 /0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, &
        & 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 3, 0, &
        & 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, &
        & 2, 0, 1, 0, 2, 1, 0, 0, 1, 1, 0, 1, 0, &
        & 1, 1, 1, 0, 0, 0, 0, 4, 0, 0, 0, 2, 0, &
        & 2, 2, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, &
        & 3, 0, 1, 3, 1, 0, 0, 1, 1, 0, 1, 0, 2, &
        & 0, 1, 1, 0, 2, 1, 0, 2, 1, 1, 2, 1, 1, &
        & 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, &
        & 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 3, 0, 0, &
        & 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, &
        & 1, 0, 0, 1, 1, 0, 0, 0, 0, 1/
    data ind6 /1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
        & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, &
        & 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
        & 2, 2, 2, 2, 3, 3, 3, 3, 3/
    data coeff /0.00482409d0, 0.00325127d0, 0.00026630d0, &
        & 0.00004135d0, 1.08154998d0, &
        & 0.86097723d0, 0.86097723d0, 0.11016112d0, &
        & 0.11016112d0, -0.03637740d0, &
        & -0.03637740d0, 0.17491854d0, 0.17491854d0, &
        & 0.00057054d0, -0.00137967d0, &
        & -0.00137967d0, 0.00062152d0, 0.00062152d0, &
        & 0.01375726d0, -1.15691421d0, &
        & -0.25495918d0, -0.25495918d0, -0.02272830d0, &
        & -0.02272830d0, -0.18381415d0, &
        & -0.18381415d0, -0.34627237d0, -0.34627237d0, &
        & 0.13974588d0, 0.13974588d0, &
        & -0.27290455d0, -0.27290455d0, -0.00674721d0, &
        & -0.00674721d0, -0.02179545d0, &
        & -0.02179545d0, -0.02125643d0, -0.02125643d0, &
        & -0.00491968d0, -0.00491968d0, &
        & 0.00233773d0, 0.00233773d0, -0.00050066d0, &
        & -0.00050066d0, 0.01817536d0, &
        & -0.04666009d0, -0.04666009d0, -0.02424748d0, &
        & -0.02424748d0, -0.01727148d0, &
        & -0.00420506d0, -0.00420506d0, -0.00647944d0, &
        & -0.00647944d0, -1.06749007d0, &
        & -0.35741007d0, -0.35741007d0, -0.00796836d0, &
        & -0.00796836d0, -0.42556742d0, &
        & -0.42556742d0, 0.06278896d0, 0.06278896d0, &
        & -0.04010419d0, -0.04010419d0, &
        & -0.00993912d0, 0.47562894d0, 0.47562894d0, &
        & -0.40830627d0, -0.40830627d0, &
        & 0.22073222d0, 0.22073222d0, 0.07828212d0, &
        & 0.07828212d0, -0.02954687d0, &
        & -0.02954687d0, 0.03057888d0, 0.03057888d0, &
        & -0.06363999d0, -0.06363999d0, &
        & -0.00373964d0, -0.00373964d0, -0.04114668d0, &
        & 0.11249614d0, 0.11249614d0, &
        & 0.02616679d0, 0.02616679d0, -0.07824425d0, &
        & 0.04266205d0, 0.04266205d0, &
        & -0.07420432d0, -0.07420432d0, -0.08251268d0, &
        & -0.08251268d0, 0.00270940d0, &
        & 0.00270940d0, 0.00199953d0, 0.00199953d0, &
        & -0.01292325d0, -0.01292325d0, &
        & -0.02074323d0, -0.02074323d0, -0.00789732d0, &
        & -0.00789732d0, -0.01435326d0, &
        & -0.00180710d0, -0.00180710d0, -0.01135671d0, &
        & -0.01135671d0, 0.00020655d0, &
        & -0.00492533d0, -0.00492533d0, 0.00270990d0, &
        & 0.00270990d0, 0.00376086d0, &
        & 0.00376086d0, 0.00044732d0, 0.00044732d0, &
        & 0.00569979d0, -0.00244774d0, &
        & -0.00244774d0, -0.02065564d0, 0.05249331d0, &
        & -0.02490299d0, -0.02490299d0, &
        & 0.00391460d0, 0.00391460d0, 0.08893744d0, &
        & 0.08893744d0, -0.01051618d0, &
        & 0.00120479d0, 0.00120479d0, -0.00111888d0, &
        & -0.00111888d0, 0.00884757d0, &
        & 0.00416289d0, 0.00416289d0, 0.00126763d0, &
        & 0.00126763d0, -0.00706563d0, &
        & -0.00706563d0, -0.00840146d0, -0.00840146d0, &
        & -0.00139219d0, 0.00801673d0, &
        & 0.00801673d0, 0.00463860d0, -0.00096051d0, &
        & 0.00019906d0, 0.00019906d0, &
        & -0.00057576d0, -0.00057576d0/
    !   convert offset parameters to bohr and radians
    rm_oo     = 1.456199d0
    rm_oh     = 0.962755d0
    angle_ooh = 100.9059d0

    pi        = 4.d0*atan(1.d0)
    angle_ooh = angle_ooh*pi/180.d0

    !
    ! Potential energy surface for H2O2
    ! Jacek Koput, Stuart Carter, and Nicholas Handy, J. Phys. Chem. A
    ! 1998, volume 102, pages 6325-6330
    !
    ! V is the potential energy in hartrees
    ! minimum configuration is at R(OH)=0.96265, R(OO)=1.45248,
    ! theta(OOH)=99.909, and dihed angle=112.456 degrees

    !      R1 = COORD(1)                ! radius from H1 to O1 (in bohr)
    !      R2 = COORD(2)                ! radius from O1 to O2 (in bohr)
    !      R3 = COORD(3)                ! radius from O2 to H2 (in bohr)
    !      A1 = COORD(4)                ! angle between H1-O1-O2 (in radians
    !      A2 = COORD(5)                ! angle between O1-O2-H2 (in radians
    !      Q6 = COORD(6)                ! HOOH torsion angle (in radians)

    r1 = coord(1)
    r2 = coord(2)
    r3 = coord(3)
    a1 = coord(4)
    a2 = coord(5)
    q6 = coord(6)

    ! Q1,Q2, and Q3 represent the stretching modes
    ! Q4 and Q5  are the bending modes
    ! Q6 is the torsional mode
    ! Q1, Q2, Q3 are dimensionless, Q4, Q5, and Q6 are in radians

    q3 = (r2 - rm_oo) / r2
    q1 = (r1 - rm_oh) / r1
    q2 = (r3 - rm_oh) / r3
    q4 = a1 - angle_ooh
    q5 = a2 - angle_ooh

    do i=0, 4
        q1p(i) = 1.d0
        q2p(i) = 1.d0
        q3p(i) = 1.d0
        q4p(i) = 1.d0
        q5p(i) = 1.d0
        q6p(i) = 1.d0
    end do

    do i=1, 4
        q1p(i) = q1p(i-1)*q1
        q2p(i) = q2p(i-1)*q2
        q3p(i) = q3p(i-1)*q3
        q4p(i) = q4p(i-1)*q4
        q5p(i) = q5p(i-1)*q5
        q6p(i) = cos(i*q6)
    end do

    v = zoe

    do i=1, 152
        v = v+coeff(i)*q1p(ind1(i))*q2p(ind2(i))*q3p(ind3(i))* q4p(ind4(i))*q5p(ind5(i))*q6p(ind6(i))
    end do

end subroutine vibpot

!                           DISCLAIMER
!
!   This file was generated on 12/09/02 by the version of
!   ADIFOR compiled on June, 1998.
!
!   ADIFOR was prepared as an account of work sponsored by an
!   agency of the United States Government, Rice University, and
!   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
!   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
!   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
!   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
!   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
!   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
!   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
!
!***********************************************************************
!
!> @brief Evaluate g_surf for the H2O2_GEN_1998 potential.
subroutine g_surf( g_p_, v, g_v, ldg_v, x, g_x, ldg_x )
    !
    !   System: H2O2, see Koput, J.; Carter, S.; Handy, N. C.;
    !              J. Phys. Chem. A, 102 (1998) 6325
    !              Prepared for Polyrate by VMA, October 2002.
    !
    !   Reference: Koput,Carter,&Handy, JPC A, 102, 6325 (1998).
    !
    !   All the information passed to and from the potential energy surface
    !   routine is in hartree atomic units.
    !
    !        This potential is written such that:
    !                       X(1)  - X(3)  : X, Y, Z for H1
    !                       X(4)  - X(6)  : X, Y, Z for O1
    !                       X(7)  - X(9)  : X, Y, Z for O2
    !                       X(10) - X(12) : X, Y, Z for H2
    !

    implicit none
    double precision            :: n, n1, n2
    double precision            :: v, pi
    double precision            :: coord(12), dx(12), x(12), pass(6)
    double precision            :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10
    double precision            :: c11, c12, c13, c14, c15, c16, c17, c18
    double precision            :: rab, rbc, rcd, rac, rbd, rad
    double precision            :: rh1o1, ro1o2, ro2h2
    double precision            :: d1, d2, d, ah1oo, aooh2, thooh
    double precision            :: xcomp1, ycomp1, zcomp1
    double precision            :: xcomp2, ycomp2, zcomp2
    double precision            :: xcomp3, ycomp3, zcomp3
    double precision            :: xcomp4, ycomp4, zcomp4
    double precision            :: x_normal1, y_normal1, z_normal1
    double precision            :: x_normal2, y_normal2, z_normal2
    double precision, parameter :: rohlow  = 1.0d0
    double precision, parameter :: rohhigh = 4.0d0
    double precision, parameter :: roolow  = 1.8d0
    double precision, parameter :: roohigh = 4.0d0
    double precision, parameter :: rhhlow  = 0.7d0
    integer, parameter          :: g_pmax_ = 12
    integer                     :: g_i_, g_p_, ldg_x, ldg_v
    double precision            :: d2_w, d1_p, d7_v, d2_p, d4_b, d5_b, d3_p, d2_v, &
        & d3_v, d7_b
    double precision :: d2_b, d3_b, d1_w, d4_v, d8_b, g_c1(g_pmax_), g_&
        &x(ldg_x,12), g_c2(g_pmax_), g_c3(g_pmax_), g_c4(g_pmax_)
    double precision :: g_c5(g_pmax_), g_c6(g_pmax_), g_c7(g_pmax_), g_&
        &c8(g_pmax_), g_c9(g_pmax_), g_c10(g_pmax_), g_c11(g_pmax_), g_c12( &
        & g_pmax_), g_c13(g_pmax_), g_c14(g_pmax_)
    double precision :: g_c15(g_pmax_), g_d1_w(g_pmax_), g_rab(g_pmax_) &
        &, g_rbc(g_pmax_), g_rcd(g_pmax_), g_rh1o1(g_pmax_), g_ro1o2(g_pmax&
        &_), g_ro2h2(g_pmax_), g_n1(g_pmax_), g_d1(g_pmax_)
    double precision :: g_ah1oo(g_pmax_), g_n2(g_pmax_), g_d2(g_pmax_), &
        & g_aooh2(g_pmax_), g_xcomp1(g_pmax_), g_ycomp1(g_pmax_), g_zcomp1( &
        & g_pmax_), g_xcomp2(g_pmax_), g_ycomp2(g_pmax_), g_zcomp2(g_pmax_)
    double precision :: g_x_normal1(g_pmax_), g_y_normal1(g_pmax_), g_z&
        &_normal1(g_pmax_), g_xcomp3(g_pmax_), g_ycomp3(g_pmax_), g_zcomp3( &
        & g_pmax_), g_xcomp4(g_pmax_), g_ycomp4(g_pmax_), g_zcomp4(g_pmax_), &
        & g_x_normal2(g_pmax_)
    double precision :: g_y_normal2(g_pmax_), g_z_normal2(g_pmax_), g_n &
        & (g_pmax_), g_d2_w(g_pmax_), g_d(g_pmax_), g_thooh(g_pmax_), g_pass &
        & (g_pmax_, 6), g_v(ldg_v)
    integer :: g_ehfid
    save g_pass
    save g_xcomp4, g_ycomp4, g_zcomp4, g_x_normal2, g_y_normal2, g_z&
        &_normal2, g_n, g_d2_w, g_d, g_thooh
    save g_zcomp1, g_xcomp2, g_ycomp2, g_zcomp2, g_x_normal1, g_y_no&
        &rmal1, g_z_normal1, g_xcomp3, g_ycomp3, g_zcomp3
    save g_ro1o2, g_ro2h2, g_n1, g_d1, g_ah1oo, g_n2, g_d2, g_aooh2, &
        & g_xcomp1, g_ycomp1
    save g_c11, g_c12, g_c13, g_c14, g_c15, g_d1_w, g_rab, g_rbc, g_&
        &rcd, g_rh1o1
    save g_c1, g_c2, g_c3, g_c4, g_c5, g_c6, g_c7, g_c8, g_c9, g_c10
    external g_vibpot
    !
    !  Calculate energy at the configuration space point
    !  NOTE :: This is not necessarily a quadrature point
    !
    ! Expecting input in cartesians
    !
    data g_ehfid /0/
    if (g_p_ .gt. g_pmax_) then
        print *, 'Parameter g_p_ is greater than g_pmax_'
        stop
    end if
    pi = 4.d0 * atan(1.d0)
    !
    do g_i_ = 1, g_p_
        g_c1(g_i_) = -g_x(g_i_, 4) + g_x(g_i_, 1)
    end do
    c1 = x(1) - x(4)
    !--------
    do g_i_ = 1, g_p_
        g_c2(g_i_) = -g_x(g_i_, 5) + g_x(g_i_, 2)
    end do
    c2 = x(2) - x(5)
    !--------
    do g_i_ = 1, g_p_
        g_c3(g_i_) = -g_x(g_i_, 6) + g_x(g_i_, 3)
    end do
    c3 = x(3) - x(6)
    !--------
    do g_i_ = 1, g_p_
        g_c4(g_i_) = -g_x(g_i_, 7) + g_x(g_i_, 4)
    end do
    c4 = x(4) - x(7)
    !--------
    do g_i_ = 1, g_p_
        g_c5(g_i_) = -g_x(g_i_, 8) + g_x(g_i_, 5)
    end do
    c5 = x(5) - x(8)
    !--------
    do g_i_ = 1, g_p_
        g_c6(g_i_) = -g_x(g_i_, 9) + g_x(g_i_, 6)
    end do
    c6 = x(6) - x(9)
    !--------
    do g_i_ = 1, g_p_
        g_c7(g_i_) = -g_x(g_i_, 10) + g_x(g_i_, 7)
    end do
    c7 = x(7) - x(10)
    !--------
    do g_i_ = 1, g_p_
        g_c8(g_i_) = -g_x(g_i_, 11) + g_x(g_i_, 8)
    end do
    c8 = x(8) - x(11)
    !--------
    do g_i_ = 1, g_p_
        g_c9(g_i_) = -g_x(g_i_, 12) + g_x(g_i_, 9)
    end do
    c9 = x(9) - x(12)
    !--------
    !
    do g_i_ = 1, g_p_
        g_c10(g_i_) = -g_x(g_i_, 7) + g_x(g_i_, 1)
    end do
    c10 = x(1) - x(7)
    !--------
    do g_i_ = 1, g_p_
        g_c11(g_i_) = -g_x(g_i_, 8) + g_x(g_i_, 2)
    end do
    c11 = x(2) - x(8)
    !--------
    do g_i_ = 1, g_p_
        g_c12(g_i_) = -g_x(g_i_, 9) + g_x(g_i_, 3)
    end do
    c12 = x(3) - x(9)
    !--------
    do g_i_ = 1, g_p_
        g_c13(g_i_) = -g_x(g_i_, 10) + g_x(g_i_, 4)
    end do
    c13 = x(4) - x(10)
    !--------
    do g_i_ = 1, g_p_
        g_c14(g_i_) = -g_x(g_i_, 11) + g_x(g_i_, 5)
    end do
    c14 = x(5) - x(11)
    !--------
    do g_i_ = 1, g_p_
        g_c15(g_i_) = -g_x(g_i_, 12) + g_x(g_i_, 6)
    end do
    c15 = x(6) - x(12)
    !--------
    c16 = x(1) - x(10)
    c17 = x(2) - x(11)
    c18 = x(3) - x(12)
    !
    !
    ! H1 -- O1 -- O2 -- H2
    !   rH1O1 rO1O2 rO2H2
    !
    !  angle H1O1O2 = AH1OO
    !  angle O1O2H2 = AOOH2
    !
    !  TORSION HOOH = THOOH
    !
    d2_v = c1 * c1
    d3_p = 2.0d0 * c1
    d4_v = c2 * c2
    d2_p = 2.0d0 * c2
    d7_v = c3 * c3
    d1_p = 2.0d0 * c3
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = d1_p * g_c3(g_i_) + d2_p * g_c2(g_i_) + d3_p * g_c1(g_i_)
    end do
    d1_w = d2_v + d4_v + d7_v
    d2_v = sqrt(d1_w)

    if ( d1_w .gt. 0.0d0 ) then
        d1_p = 1.0d0 / (2.0d0 *  d2_v)
    else
        !           call ehufDO (9,d1_w, d2_v, d1_p,
        !     +g_ehfid,
        !     +201)
    end if
    do g_i_ = 1, g_p_
        g_rab(g_i_) = d1_p * g_d1_w(g_i_)
    end do
    rab = d2_v
    !--------
    d2_v = c4 * c4
    d3_p = 2.0d0 * c4
    d4_v = c5 * c5
    d2_p = 2.0d0 * c5
    d7_v = c6 * c6
    d1_p = 2.0d0 * c6
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = d1_p * g_c6(g_i_) + d2_p * g_c5(g_i_) + d3_p * g_c4(g_i_)
    end do
    d1_w = d2_v + d4_v + d7_v
    d2_v = sqrt(d1_w)

    if ( d1_w .gt. 0.0d0 ) then
        d1_p = 1.0d0 / (2.0d0 *  d2_v)
    else
        !           call ehufDO (9,d1_w, d2_v, d1_p,
        !     +g_ehfid,
        !     +226)
    end if
    do g_i_ = 1, g_p_
        g_rbc(g_i_) = d1_p * g_d1_w(g_i_)
    end do
    rbc = d2_v
    !--------
    d2_v = c7 * c7
    d3_p = 2.0d0 * c7
    d4_v = c8 * c8
    d2_p = 2.0d0 * c8
    d7_v = c9 * c9
    d1_p = 2.0d0 * c9
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = d1_p * g_c9(g_i_) + d2_p * g_c8(g_i_) + d3_p * g_c7(g_i_)
    end do
    d1_w = d2_v + d4_v + d7_v
    d2_v = sqrt(d1_w)

    if ( d1_w .gt. 0.0d0 ) then
        d1_p = 1.0d0 / (2.0d0 *  d2_v)
    else
        !           call ehufDO (9,d1_w, d2_v, d1_p,
        !     +g_ehfid,
        !     +251)
    end if
    do g_i_ = 1, g_p_
        g_rcd(g_i_) = d1_p * g_d1_w(g_i_)
    end do
    rcd = d2_v
    !--------
    rac = sqrt(c10 ** 2 + c11 ** 2 + c12 ** 2)
    rbd = sqrt(c13 ** 2 + c14 ** 2 + c15 ** 2)
    rad = sqrt(c16 ** 2 + c17 ** 2 + c18 ** 2)
    do g_i_ = 1, g_p_
        g_rh1o1(g_i_) = g_rab(g_i_)
    end do
    rh1o1 = rab
    !--------
    do g_i_ = 1, g_p_
        g_ro1o2(g_i_) = g_rbc(g_i_)
    end do
    ro1o2 = rbc
    !--------
    do g_i_ = 1, g_p_
        g_ro2h2(g_i_) = g_rcd(g_i_)
    end do
    ro2h2 = rcd
    !--------
    !
    ! The potential can give weird values for unusual geometries
    ! replace with some high value in these unphysical regions and
    ! then exit
    !
    !
    !      if((min(rab,rcd,rac,rbd).gt.rohhigh) .or.
    !     !(min(rab,rcd,rac,rbd).lt.rohlow)  .or.
    !     !(rbc.gt.roohigh) .or. (rbc.lt.roolow) .or.
    !     !(rad.lt.rhhlow)) then
    !      vdc=1000.d0
    !      return
    !      endif
    !
    do g_i_ = 1, g_p_
        g_n1(g_i_) = (-c6) * g_c3(g_i_) + (-c3) * g_c6(g_i_) + (-c5) * g_c2(g_i_) + (-c2) * g_c5(g_i_) + (-c4) * &
            & g_c1(g_i_) + (-c1) * g_c4(g_i_)
    end do
    n1 = (-c4) * c1 + (-c5) * c2 + (-c6) * c3
    !--------
    do g_i_ = 1, g_p_
        g_d1(g_i_) = ro1o2 * g_rh1o1(g_i_) + rh1o1 * g_ro1o2(g_i_)
    end do
    d1 = ro1o2 * rh1o1
    !--------
    d3_v = n1 / d1
    d2_b = 1.0d0 / d1
    d3_b = (-d3_v) / d1
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = d3_b * g_d1(g_i_) + d2_b * g_n1(g_i_)
    end do
    d1_w = d3_v
    d2_v = acos(d1_w)

    if ( abs(d1_w) .lt. 1.0d0 ) then
        d1_p = -1.0d0 / sqrt ((1.0d0-d1_w)*(1.0d0+d1_w))
    else
        !           call ehufDO (14,d1_w, d2_v, d1_p,
        !     +g_ehfid,
        !     +316)
    end if
    do g_i_ = 1, g_p_
        g_ah1oo(g_i_) = d1_p * g_d1_w(g_i_)
    end do
    ah1oo = d2_v
    !--------
    !
    do g_i_ = 1, g_p_
        g_n2(g_i_) = (-c9) * g_c6(g_i_) + (-c6) * g_c9(g_i_) + (-c8) * g_c5(g_i_) + (-c5) * g_c8(g_i_) + (-c7) * &
            & g_c4(g_i_) + (-c4) * g_c7(g_i_)
    end do
    n2 = -(c7 * c4 + c8 * c5 + c9 * c6)
    !--------
    do g_i_ = 1, g_p_
        g_d2(g_i_) = ro2h2 * g_ro1o2(g_i_) + ro1o2 * g_ro2h2(g_i_)
    end do
    d2 = ro2h2 * ro1o2
    !--------
    d3_v = n2 / d2
    d2_b = 1.0d0 / d2
    d3_b = (-d3_v) / d2
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = d3_b * g_d2(g_i_) + d2_b * g_n2(g_i_)
    end do
    d1_w = d3_v
    d2_v = acos(d1_w)

    if ( abs(d1_w) .lt. 1.0d0 ) then
        d1_p = -1.0d0 / sqrt ((1.0d0-d1_w)*(1.0d0+d1_w))
    else
        !           call ehufDO (14,d1_w, d2_v, d1_p,
        !     +g_ehfid,
        !     +350)
    end if
    do g_i_ = 1, g_p_
        g_aooh2(g_i_) = d1_p * g_d1_w(g_i_)
    end do
    aooh2 = d2_v
    !--------
    !
    ! Calculate torsion angle HOOH
    !
    do g_i_ = 1, g_p_
        g_xcomp1(g_i_) = -g_x(g_i_, 1) + g_x(g_i_, 4)
    end do
    xcomp1 = x(4) - x(1)
    !--------
    do g_i_ = 1, g_p_
        g_ycomp1(g_i_) = -g_x(g_i_, 2) + g_x(g_i_, 5)
    end do
    ycomp1 = x(5) - x(2)
    !--------
    do g_i_ = 1, g_p_
        g_zcomp1(g_i_) = -g_x(g_i_, 3) + g_x(g_i_, 6)
    end do
    zcomp1 = x(6) - x(3)
    !--------
    !
    do g_i_ = 1, g_p_
        g_xcomp2(g_i_) = g_c10(g_i_)
    end do
    xcomp2 = c10
    !--------
    do g_i_ = 1, g_p_
        g_ycomp2(g_i_) = g_c11(g_i_)
    end do
    ycomp2 = c11
    !--------
    do g_i_ = 1, g_p_
        g_zcomp2(g_i_) = g_c12(g_i_)
    end do
    zcomp2 = c12
    !--------
    !
    do g_i_ = 1, g_p_
        g_x_normal1(g_i_) = (-zcomp1) * g_ycomp2(g_i_) + (-ycomp2) * g_zcomp1(g_i_) + ycomp1 * g_zcomp2(g_i_) + &
            & zcomp2 * g_ycomp1(g_i_)
    end do
    x_normal1 = ycomp1 * zcomp2 - zcomp1 * ycomp2
    !--------
    do g_i_ = 1, g_p_
        g_y_normal1(g_i_) = (-xcomp1) * g_zcomp2(g_i_) + (-zcomp2) * g_xcomp1(g_i_) + zcomp1 * g_xcomp2(g_i_) + &
            & xcomp2 * g_zcomp1(g_i_)
    end do
    y_normal1 = zcomp1 * xcomp2 - xcomp1 * zcomp2
    !--------
    do g_i_ = 1, g_p_
        g_z_normal1(g_i_) = (-ycomp1) * g_xcomp2(g_i_) + (-xcomp2) * g_ycomp1(g_i_) + xcomp1 * g_ycomp2(g_i_) + &
            & ycomp2 * g_xcomp1(g_i_)
    end do
    z_normal1 = xcomp1 * ycomp2 - ycomp1 * xcomp2
    !--------
    !
    do g_i_ = 1, g_p_
        g_xcomp3(g_i_) = -g_x(g_i_, 4) + g_x(g_i_, 7)
    end do
    xcomp3 = x(7) - x(4)
    !--------
    do g_i_ = 1, g_p_
        g_ycomp3(g_i_) = -g_x(g_i_, 5) + g_x(g_i_, 8)
    end do
    ycomp3 = x(8) - x(5)
    !--------
    do g_i_ = 1, g_p_
        g_zcomp3(g_i_) = -g_x(g_i_, 6) + g_x(g_i_, 9)
    end do
    zcomp3 = x(9) - x(6)
    !--------
    !
    do g_i_ = 1, g_p_
        g_xcomp4(g_i_) = g_c13(g_i_)
    end do
    xcomp4 = c13
    !--------
    do g_i_ = 1, g_p_
        g_ycomp4(g_i_) = g_c14(g_i_)
    end do
    ycomp4 = c14
    !--------
    do g_i_ = 1, g_p_
        g_zcomp4(g_i_) = g_c15(g_i_)
    end do
    zcomp4 = c15
    !--------
    !
    do g_i_ = 1, g_p_
        g_x_normal2(g_i_) = (-zcomp3) * g_ycomp4(g_i_) + (-ycomp4) * g_zcomp3(g_i_) + ycomp3 * g_zcomp4(g_i_) + &
            & zcomp4 * g_ycomp3(g_i_)
    end do
    x_normal2 = ycomp3 * zcomp4 - zcomp3 * ycomp4
    !--------
    do g_i_ = 1, g_p_
        g_y_normal2(g_i_) = (-xcomp3) * g_zcomp4(g_i_) + (-zcomp4) * g_xcomp3(g_i_) + zcomp3 * g_xcomp4(g_i_) + &
            & xcomp4 * g_zcomp3(g_i_)
    end do
    y_normal2 = zcomp3 * xcomp4 - xcomp3 * zcomp4
    !--------
    do g_i_ = 1, g_p_
        g_z_normal2(g_i_) = (-ycomp3) * g_xcomp4(g_i_) + (-xcomp4) * g_ycomp3(g_i_) + xcomp3 * g_ycomp4(g_i_) + &
            & ycomp4 * g_xcomp3(g_i_)
    end do
    z_normal2 = xcomp3 * ycomp4 - ycomp3 * xcomp4
    !--------
    !
    do g_i_ = 1, g_p_
        g_n(g_i_) = z_normal1 * g_z_normal2(g_i_) + z_normal2 * g_z_normal1(g_i_) + y_normal1 * g_y_normal2(g_i_) + &
            & y_normal2 * g_y_normal1(g_i_) + x_normal1 * g_x_normal2(g_i_) + x_normal2 * g_x_normal1(g_i_)
    end do
    n = x_normal1 * x_normal2 + y_normal1 * y_normal2 + z_normal1 * z_normal2
    !--------
    !
    d4_b = z_normal1 + z_normal1
    d7_b = y_normal1 + y_normal1
    d8_b = x_normal1 + x_normal1
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = d4_b * g_z_normal1(g_i_) + d7_b * g_y_normal1(g_i_) + d8_b * g_x_normal1(g_i_)
    end do
    d1_w = x_normal1 * x_normal1 + y_normal1 * y_normal1 + z_normal1 * z_normal1
    d4_b = z_normal2 + z_normal2
    d7_b = y_normal2 + y_normal2
    d8_b = x_normal2 + x_normal2
    do g_i_ = 1, g_p_
        g_d2_w(g_i_) = d4_b * g_z_normal2(g_i_) + d7_b * g_y_normal2(g_i_) + d8_b * g_x_normal2(g_i_)
    end do
    d2_w = x_normal2 * x_normal2 + y_normal2 * y_normal2 + z_normal2 * z_normal2
    d2_v = sqrt(d1_w)

    if ( d1_w .gt. 0.0d0 ) then
        d2_p = 1.0d0 / (2.0d0 *  d2_v)
    else
        !           call ehufDO (9,d1_w, d2_v, d2_p,
        !     +g_ehfid,
        !     +497)
    end if
    d4_v = sqrt(d2_w)

    if ( d2_w .gt. 0.0d0 ) then
        d1_p = 1.0d0 / (2.0d0 *  d4_v)
    else
        !           call ehufDO (9,d2_w, d4_v, d1_p,
        !     +g_ehfid,
        !     +506)
    end if
    d4_b = d2_v * d1_p
    d5_b = d4_v * d2_p
    do g_i_ = 1, g_p_
        g_d(g_i_) = d4_b * g_d2_w(g_i_) + d5_b * g_d1_w(g_i_)
    end do
    d = d2_v * d4_v
    !--------
    !
    !  We can occasionally get "NaNs" if we don't trap for special cases whi
    !lating
    !  the dihedral angle.  Acos wants an argument in [-1,1], but slight rou
    !n push us
    !  outside this range.  Also if certain triatoms are collinear then the
    ! is not defined
    !  and we set these cases to 0 degrees
    !
    !      if(d.eq.0.d0)then
    !       thooh=0.d0
    !      else
    !       argnd=n/d
    !      if(argnd.ge.1.d0)then
    !       thooh=0.d0
    !      elseif(argnd.le.-1.d0)then
    !       thooh=pi
    !      else
    d3_v = n / d
    d2_b = 1.0d0 / d
    d3_b = (-d3_v) / d
    do g_i_ = 1, g_p_
        g_d1_w(g_i_) = d3_b * g_d(g_i_) + d2_b * g_n(g_i_)
    end do
    d1_w = d3_v
    d2_v = acos(d1_w)

    if ( abs(d1_w) .lt. 1.0d0 ) then
        d1_p = -1.0d0 / sqrt ((1.0d0-d1_w)*(1.0d0+d1_w))
    else
        !           call ehufDO (14,d1_w, d2_v, d1_p,
        !     +g_ehfid,
        !     +547)
    end if
    do g_i_ = 1, g_p_
        g_thooh(g_i_) = d1_p * g_d1_w(g_i_)
    end do
    thooh = d2_v
    !--------
    !      endif
    !      endif
    !
    do g_i_ = 1, g_p_
        g_pass(g_i_, 1) = g_rh1o1(g_i_)
    end do
    pass(1) = rh1o1
    !--------
    do g_i_ = 1, g_p_
        g_pass(g_i_, 2) = g_ro1o2(g_i_)
    end do
    pass(2) = ro1o2
    !--------
    do g_i_ = 1, g_p_
        g_pass(g_i_, 3) = g_ro2h2(g_i_)
    end do
    pass(3) = ro2h2
    !--------
    do g_i_ = 1, g_p_
        g_pass(g_i_, 4) = g_ah1oo(g_i_)
    end do
    pass(4) = ah1oo
    !--------
    do g_i_ = 1, g_p_
        g_pass(g_i_, 5) = g_aooh2(g_i_)
    end do
    pass(5) = aooh2
    !--------
    do g_i_ = 1, g_p_
        g_pass(g_i_, 6) = g_thooh(g_i_)
    end do
    pass(6) = thooh
    !--------
    !
    !
    ! Call Koput, Carter, and Handy PES
    !
    !
    call g_vibpot(g_p_, pass, g_pass, g_pmax_, v, g_v, ldg_v)
    !
    !
end subroutine g_surf

!
!***********************************************************************
!
!> @brief Evaluate g_vibpot for the H2O2_GEN_1998 potential.
subroutine g_vibpot( g_p_, coord, g_coord, ldg_coord, v, g_v, ldg_v &
        & )
    !

    implicit none
    double precision   :: v, pi
    double precision   :: coord(6), q1p(0:4), q2p(0:4), q3p(0:4), q4p(0:4)
    double precision   :: q5p(0:4), q6p(0:4)
    double precision   :: ind1(152), ind2(152), ind3(152), ind4(152)
    double precision   :: ind5(152), ind6(152), coeff(152)
    double precision   :: rm_oo, rm_oh, angle_ooh, r1, r2, r3
    double precision   :: a1, a2, q1, q2, q3, q4, q5, q6
    integer            :: i
    double precision   :: zoe
    integer, parameter :: g_pmax_ = 12
    integer            :: g_i_, g_p_, ldg_coord, ldg_v
    double precision   :: d14_b, d13_b, d7_v, d11_b, d10_b, d9_b, d8_b, d&
        &7_b, d2_v, d3_v
    double precision :: d6_b, d2_b, d3_b, d1_w, d1_p, d11_v, d5_v, d9_v &
        &, g_r1(g_pmax_), g_coord(ldg_coord, 6)
    double precision :: g_r2(g_pmax_), g_r3(g_pmax_), g_a1(g_pmax_), g_&
        &a2(g_pmax_), g_q6(g_pmax_), g_q3(g_pmax_), g_q1(g_pmax_), g_q2(g_p&
        &max_), g_q4(g_pmax_), g_q5(g_pmax_)
    double precision :: g_q1p(g_pmax_,0:4), g_q2p(g_pmax_,0:4), g_q3p &
        & (g_pmax_, 0:4), g_q4p(g_pmax_, 0:4), g_q5p(g_pmax_, 0:4), g_q6p(g_&
        &pmax_, 0:4), g_d1_w(g_pmax_), g_v(ldg_v)
    integer :: g_ehfid
    save g_q5, g_q1p, g_q2p, g_q3p, g_q4p, g_q5p, g_q6p, g_d1_w
    save g_r1, g_r2, g_r3, g_a1, g_a2, g_q6, g_q3, g_q1, g_q2, g_q4
    intrinsic dble
    data zoe /3.995139982058d-3/
    !
    data ind1 /0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, &
        & 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 2, 1, 0, 2, 0, 1, 0, &
        & 2, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 4, 0, 0, 0, 2, 0, 0, &
        & 0, 2, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 3, 1, 0, 3, 0, 0, 0, 1, 1, &
        & 0, 1, 0, 0, 1, 0, 2, 0, 2, 0, 0, 0, 1, 1, 2, 0, 1, 0, 1, 0, 0, 1, &
        & 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 3, 0, 0, 0, 2, &
        & 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, &
        & 0, 0/
    !
    data ind2 /0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, &
        & 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 2, 1, 0, 1, 0, 2, 0, 1, &
        & 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 4, 0, 0, 0, 2, 0, &
        & 0, 0, 2, 0, 0, 1, 0, 0, 0, 3, 0, 0, 3, 1, 0, 1, 0, 3, 0, 0, 1, 0, &
        & 1, 0, 1, 0, 0, 1, 0, 2, 0, 2, 0, 0, 1, 1, 0, 2, 0, 1, 0, 1, 0, 0, &
        & 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 3, 0, 0, 0, &
        & 2, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, &
        & 0, 0/
    !
    data ind3 /0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, &
        & 0, 3, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 4, 0, 0, 0, 0, 2, 2, 2, &
        & 2, 0, 0, 0, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, &
        & 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, &
        & 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, &
        & 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, &
        & 0, 0/
    !
    data ind4 /0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, &
        & 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, 0, 0, 2, &
        & 0, 1, 1, 2, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 4, 0, 0, 0, 2, &
        & 0, 2, 0, 2, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 3, 0, 1, 0, 1, 3, 0, 1, &
        & 0, 0, 1, 1, 2, 0, 1, 0, 0, 1, 1, 2, 2, 0, 1, 1, 1, 2, 2, 1, 0, 0, &
        & 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 3, 0, 0, &
        & 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, &
        & 1, 0/
    !
    data ind5 /0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 1, 0, &
        & 1, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, 2, 0, &
        & 1, 0, 2, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 4, 0, 0, 0, &
        & 2, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 3, 0, 1, 3, 1, 0, 0, &
        & 1, 1, 0, 1, 0, 2, 0, 1, 1, 0, 2, 1, 0, 2, 1, 1, 2, 1, 1, 2, 0, 0, &
        & 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 3, 0, &
        & 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, &
        & 0, 1/
    !
    data ind6 /1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, &
        & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
        & 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, &
        & 3, 3/
    !
    data coeff /0.00482409d0, 0.00325127d0, 0.00026630d0, 0.00004135&
        &d0, 1.08154998d0, 0.86097723d0, 0.86097723d0, 0.11016112d0, 0.1101&
        &6112d0, -0.03637740d0, -0.03637740d0, 0.17491854d0, 0.17491854d0, &
        & 0.00057054d0, -0.00137967d0, -0.00137967d0, 0.00062152d0, 0.000621&
        &52d0, 0.01375726d0, -1.15691421d0, -0.25495918d0, -0.25495918d0, -&
        &0.02272830d0, -0.02272830d0, -0.18381415d0, -0.18381415d0, -0.3462&
        &7237d0, -0.34627237d0, 0.13974588d0, 0.13974588d0, -0.27290455d0, &
        & -0.27290455d0, -0.00674721d0, -0.00674721d0, -0.02179545d0, -0.021&
        &79545d0, -0.02125643d0, -0.02125643d0, -0.00491968d0, -0.00491968d&
        &0, 0.00233773d0, 0.00233773d0, -0.00050066d0, -0.00050066d0, 0.018&
        &17536d0, -0.04666009d0, -0.04666009d0, -0.02424748d0, -0.02424748d&
        &0, -0.01727148d0, -0.00420506d0, -0.00420506d0, -0.00647944d0, -0.&
        &00647944d0, -1.06749007d0, -0.35741007d0, -0.35741007d0, -0.007968&
        &36d0, -0.00796836d0, -0.42556742d0, -0.42556742d0, 0.06278896d0, 0&
        &.06278896d0, -0.04010419d0, -0.04010419d0, -0.00993912d0, 0.475628&
        &94d0, 0.47562894d0, -0.40830627d0, -0.40830627d0, 0.22073222d0, 0.&
        &22073222d0, 0.07828212d0, 0.07828212d0, -0.02954687d0, -0.02954687&
        &d0, 0.03057888d0, 0.03057888d0, -0.06363999d0, -0.06363999d0, -0.0&
        &0373964d0, -0.00373964d0, -0.04114668d0, 0.11249614d0, 0.11249614d&
        &0, 0.02616679d0, 0.02616679d0, -0.07824425d0, 0.04266205d0, 0.0426&
        &6205d0, -0.07420432d0, -0.07420432d0, -0.08251268d0, -0.08251268d0 &
        &, 0.00270940d0, 0.00270940d0, 0.00199953d0, 0.00199953d0, -0.01292&
        &325d0, -0.01292325d0, -0.02074323d0, -0.02074323d0, -0.00789732d0, &
        & -0.00789732d0, -0.01435326d0, -0.00180710d0, -0.00180710d0, -0.01&
        &135671d0, -0.01135671d0, 0.00020655d0, -0.00492533d0, -0.00492533d&
        &0, 0.00270990d0, 0.00270990d0, 0.00376086d0, 0.00376086d0, 0.00044&
        &732d0, 0.00044732d0, 0.00569979d0, -0.00244774d0, -0.00244774d0, -&
        &0.02065564d0, 0.05249331d0, -0.02490299d0, -0.02490299d0, 0.003914&
        &60d0, 0.00391460d0, 0.08893744d0, 0.08893744d0, -0.01051618d0, 0.0&
        &0120479d0, 0.00120479d0, -0.00111888d0, -0.00111888d0, 0.00884757d&
        &0, 0.00416289d0, 0.00416289d0, 0.00126763d0, 0.00126763d0, -0.0070&
        &6563d0, -0.00706563d0, -0.00840146d0, -0.00840146d0, -0.00139219d0 &
        &, 0.00801673d0, 0.00801673d0, 0.00463860d0, -0.00096051d0, 0.00019&
        &906d0, 0.00019906d0, -0.00057576d0, -0.00057576d0/
    !
    !      data rm_oo, rm_oh, ANGLE_OOH/1.456199d0, 0.962755d0, 100.9059d0/
    !
    !   convert offset parameters to bohr and radians
    !
    data g_ehfid /0/
    !
    !        call ehsfid(g_ehfid, 'vibpot','g_surf.f')
    !
    if (g_p_ .gt. g_pmax_) then
        print *, 'Parameter g_p_ is greater than g_pmax_'
        stop
    end if
    rm_oo     = 1.456199d0
    rm_oh     = 0.962755d0
    angle_ooh = 100.9059d0
    !
    pi = 4.d0 * atan(1.d0)
    !        rm_oo = rm_oo / 0.529177d0
    !        rm_oh = rm_oh / 0.529177d0
    angle_ooh = angle_ooh * pi / 180.d0
    !
    !
    ! Potential energy surface for H2O2
    ! Jacek Koput, Stuart Carter, and Nicholas Handy, J. Phys. Chem. A
    ! 1998, volume 102, pages 6325-6330
    !
    ! V is the potential energy in hartrees
    ! minimum configuration is at R(OH)=0.96265, R(OO)=1.45248,
    ! theta(OOH)=99.909, and dihed angle=112.456 degrees
    !
    !      R1 = COORD(1)                ! radius from H1 to O1 (in bohr)
    !      R2 = COORD(2)                ! radius from O1 to O2 (in bohr)
    !      R3 = COORD(3)                ! radius from O2 to H2 (in bohr)
    !      A1 = COORD(4)                ! angle between H1-O1-O2 (in radians
    !      A2 = COORD(5)                ! angle between O1-O2-H2 (in radians
    !      Q6 = COORD(6)                ! HOOH torsion angle (in radians)
    !
    do g_i_ = 1, g_p_
        g_r1(g_i_) = g_coord(g_i_, 1)
    end do
    r1 = coord(1)
    !--------
    do g_i_ = 1, g_p_
        g_r2(g_i_) = g_coord(g_i_, 2)
    end do
    r2 = coord(2)
    !--------
    do g_i_ = 1, g_p_
        g_r3(g_i_) = g_coord(g_i_, 3)
    end do
    r3 = coord(3)
    !--------
    do g_i_ = 1, g_p_
        g_a1(g_i_) = g_coord(g_i_, 4)
    end do
    a1 = coord(4)
    !--------
    do g_i_ = 1, g_p_
        g_a2(g_i_) = g_coord(g_i_, 5)
    end do
    a2 = coord(5)
    !--------
    do g_i_ = 1, g_p_
        g_q6(g_i_) = g_coord(g_i_, 6)
    end do
    q6 = coord(6)
    !--------
    !
    ! Q1,Q2, and Q3 represent the stretching modes
    ! Q4 and Q5  are the bending modes
    ! Q6 is the torsional mode
    ! Q1, Q2, Q3 are dimensionless, Q4, Q5, and Q6 are in radians
    !
    d3_v = (r2 - rm_oo) / r2
    d3_b = (-d3_v) / r2 + 1.0d0 / r2
    do g_i_ = 1, g_p_
        g_q3(g_i_) = d3_b * g_r2(g_i_)
    end do
    q3 = d3_v
    !--------
    d3_v = (r1 - rm_oh) / r1
    d3_b = (-d3_v) / r1 + 1.0d0 / r1
    do g_i_ = 1, g_p_
        g_q1(g_i_) = d3_b * g_r1(g_i_)
    end do
    q1 = d3_v
    !--------
    d3_v = (r3 - rm_oh) / r3
    d3_b = (-d3_v) / r3 + 1.0d0 / r3
    do g_i_ = 1, g_p_
        g_q2(g_i_) = d3_b * g_r3(g_i_)
    end do
    q2 = d3_v
    !--------
    do g_i_ = 1, g_p_
        g_q4(g_i_) = g_a1(g_i_)
    end do
    q4 = a1 - angle_ooh
    !--------
    do g_i_ = 1, g_p_
        g_q5(g_i_) = g_a2(g_i_)
    end do
    q5 = a2 - angle_ooh
    !--------
    !
    !
    !
    do i = 0, 4
        do g_i_ = 1, g_p_
            g_q1p(g_i_,i) = 0.0d0
        end do
        q1p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            g_q2p(g_i_,i) = 0.0d0
        end do
        q2p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            g_q3p(g_i_, i) = 0.0d0
        end do
        q3p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            g_q4p(g_i_, i) = 0.0d0
        end do
        q4p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            g_q5p(g_i_, i) = 0.0d0
        end do
        q5p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            g_q6p(g_i_, i) = 0.0d0
        end do
        q6p(i) = 1.d0
        !--------
    end do
    !
    !
    do i = 1, 4
        do g_i_ = 1, g_p_
            g_q1p(g_i_,i) = q1p(i - 1) * g_q1(g_i_) + q1 * g_q1p(g_i_,i - 1)
        end do
        q1p(i) = q1p(i - 1) * q1
        !--------
        do g_i_ = 1, g_p_
            g_q2p(g_i_,i) = q2p(i - 1) * g_q2(g_i_) + q2 * g_q2p(g_i_,i - 1)
        end do
        q2p(i) = q2p(i - 1) * q2
        !--------
        do g_i_ = 1, g_p_
            g_q3p(g_i_, i) = q3p(i - 1) * g_q3(g_i_) + q3 * g_q3p(g_i_, i - 1)
        end do
        q3p(i) = q3p(i - 1) * q3
        !--------
        do g_i_ = 1, g_p_
            g_q4p(g_i_, i) = q4p(i - 1) * g_q4(g_i_) + q4 * g_q4p(g_i_, i - 1)
        end do
        q4p(i) = q4p(i - 1) * q4
        !--------
        do g_i_ = 1, g_p_
            g_q5p(g_i_, i) = q5p(i - 1) * g_q5(g_i_) + q5 * g_q5p(g_i_, i - 1)
        end do
        q5p(i) = q5p(i - 1) * q5
        !--------
        d2_b = dble(i)
        do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d2_b * g_q6(g_i_)
        end do
        d1_w = dble(i) * q6
        d2_v = cos(d1_w)
        d1_p = -sin(d1_w)
        do g_i_ = 1, g_p_
            g_q6p(g_i_, i) = d1_p * g_d1_w(g_i_)
        end do
        q6p(i) = d2_v
        !--------
    end do
    !
    do g_i_ = 1, g_p_
        g_v(g_i_) = 0.0d0
    end do
    v = zoe
    !--------
    !
    do i = 1, 152
        d3_v  = coeff(i) * q1p(ind1(i))
        d5_v  = d3_v * q2p(ind2(i))
        d7_v  = d5_v * q3p(ind3(i))
        d9_v  = d7_v * q4p(ind4(i))
        d11_v = d9_v * q5p(ind5(i))
        d6_b  = q6p(ind6(i)) * q5p(ind5(i))
        d7_b  = q6p(ind6(i)) * d9_v
        d8_b  = d6_b * q4p(ind4(i))
        d9_b  = d6_b * d7_v
        d10_b = d8_b * q3p(ind3(i))
        d11_b = d8_b * d5_v
        d13_b = d10_b * d3_v
        d14_b = d10_b * q2p(ind2(i)) * coeff(i)
        do g_i_ = 1, g_p_
            g_v(g_i_) = d11_v * g_q6p(g_i_, ind6(i)) + d7_b * g_q5p(g_i_, ind5(i)) + d9_b * g_q4p(g_i_, ind4(i)) + &
                & d11_b * g_q3p(g_i_, ind3(i)) + d13_b * g_q2p(g_i_,ind2(i)) + d14_b * g_q1p(g_i_,ind1(i) ) + g_v(g_i_)
        end do
        v = v + d11_v * q6p(ind6(i))
        !--------
    end do
    !
    !
    !
end subroutine g_vibpot

!
!
!
!                           DISCLAIMER
!
!   This file was generated on 01/03/03 by the version of
!   ADIFOR compiled on June, 1998.
!
!   ADIFOR was prepared as an account of work sponsored by an
!   agency of the United States Government, Rice University, and
!   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
!   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
!   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
!   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
!   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
!   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
!   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
!
!                           DISCLAIMER
!
!   This file was generated on 01/03/03 by the version of
!   ADIFOR compiled on June, 1998.
!
!   ADIFOR was prepared as an account of work sponsored by an
!   agency of the United States Government, Rice University, and
!   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
!   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
!   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
!   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
!   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
!   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
!   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
!
!***********************************************************************
!

!      CALL G_ANALYTIC(GPP, GP, V, GV, GVV, LDGVV, LDGV, X,
!    * GX, GXX, LDGXX, LDGX)

!      CALL G_ANALYTIC(12, 12, V, GV, GVV, 1, 1, X, GX, GXX, 12, 12)

!> @brief Evaluate G_ANALYTIC for the H2O2_GEN_1998 potential.
subroutine g_analytic( g2p, g_p_, v, g_v, g_g_v, ldg_g_v, ldg_v, x, &
        & g_x, g_g_x, ldg_g_x, ldg_x )
    !
    !   System: H2O2, see Koput, J.; Carter, S.; Handy, N. C.;
    !              J. Phys. Chem. A, 102 (1998) 6325
    !              Prepared for Polyrate by VMA, October 2002.
    !
    !   Reference: Koput,Carter,&Handy, JPC A, 102, 6325 (1998).
    !
    !   All the information passed to and from the potential energy surface
    !   routine is in hartree atomic units.
    !
    !        This potential is written such that:
    !                       X(1)  - X(3)  : X, Y, Z for H1
    !                       X(4)  - X(6)  : X, Y, Z for O1
    !                       X(7)  - X(9)  : X, Y, Z for O2
    !                       X(10) - X(12) : X, Y, Z for H2
    !

    implicit none
    double precision            :: n, n1, n2
    double precision            :: v, pi
    double precision            :: coord(12), dx(12), x(12), pass(6)
    double precision            :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10
    double precision            :: c11, c12, c13, c14, c15, c16, c17, c18
    double precision            :: rab, rbc, rcd, rac, rbd, rad
    double precision            :: rh1o1, ro1o2, ro2h2
    double precision            :: d1, d2, d, ah1oo, aooh2, thooh
    double precision            :: xcomp1, ycomp1, zcomp1
    double precision            :: xcomp2, ycomp2, zcomp2
    double precision            :: xcomp3, ycomp3, zcomp3
    double precision            :: xcomp4, ycomp4, zcomp4
    double precision            :: x_normal1, y_normal1, z_normal1
    double precision            :: x_normal2, y_normal2, z_normal2
    double precision, parameter :: rohlow  = 1.0d0
    double precision, parameter :: rohhigh = 4.0d0
    double precision, parameter :: roolow  = 1.8d0
    double precision, parameter :: roohigh = 4.0d0
    double precision, parameter :: rhhlow  = 0.7d0
    integer, parameter          :: g_pmax_ = 12
    integer, parameter          :: g2pmax  = 12
    integer                     :: g_i_, g_p_, ldg_x, ldg_v
    double precision            :: d2_w, d1_p, d7_v, d2_p, d4_b, d5_b, d3_p, d2_v, &
        & d3_v, d7_b
    double precision :: d2_b, d3_b, d1_w, d4_v, d8_b, g_c1(g_pmax_), g_&
        &x(ldg_x,12), g_c2(g_pmax_), g_c3(g_pmax_), g_c4(g_pmax_)
    double precision :: g_c5(g_pmax_), g_c6(g_pmax_), g_c7(g_pmax_), g_&
        &c8(g_pmax_), g_c9(g_pmax_), g_c10(g_pmax_), g_c11(g_pmax_), g_c12( &
        & g_pmax_), g_c13(g_pmax_), g_c14(g_pmax_)
    double precision :: g_c15(g_pmax_), g_d1_w(g_pmax_), g_rab(g_pmax_) &
        &, g_rbc(g_pmax_), g_rcd(g_pmax_), g_rh1o1(g_pmax_), g_ro1o2(g_pmax&
        &_), g_ro2h2(g_pmax_), g_n1(g_pmax_), g_d1(g_pmax_)
    double precision :: g_ah1oo(g_pmax_), g_n2(g_pmax_), g_d2(g_pmax_), &
        & g_aooh2(g_pmax_), g_xcomp1(g_pmax_), g_ycomp1(g_pmax_), g_zcomp1( &
        & g_pmax_), g_xcomp2(g_pmax_), g_ycomp2(g_pmax_), g_zcomp2(g_pmax_)
    double precision :: g_x_normal1(g_pmax_), g_y_normal1(g_pmax_), g_z&
        &_normal1(g_pmax_), g_xcomp3(g_pmax_), g_ycomp3(g_pmax_), g_zcomp3( &
        & g_pmax_), g_xcomp4(g_pmax_), g_ycomp4(g_pmax_), g_zcomp4(g_pmax_), &
        & g_x_normal2(g_pmax_)
    double precision :: g_y_normal2(g_pmax_), g_z_normal2(g_pmax_), g_n &
        & (g_pmax_), g_d2_w(g_pmax_), g_d(g_pmax_), g_thooh(g_pmax_), g_pass &
        & (g_pmax_, 6), g_v(ldg_v)
    integer :: g_ehfid
    save g_pass
    save g_xcomp4, g_ycomp4, g_zcomp4, g_x_normal2, g_y_normal2, g_z&
        &_normal2, g_n, g_d2_w, g_d, g_thooh
    save g_zcomp1, g_xcomp2, g_ycomp2, g_zcomp2, g_x_normal1, g_y_no&
        &rmal1, g_z_normal1, g_xcomp3, g_ycomp3, g_zcomp3
    save g_ro1o2, g_ro2h2, g_n1, g_d1, g_ah1oo, g_n2, g_d2, g_aooh2, &
        & g_xcomp1, g_ycomp1
    save g_c11, g_c12, g_c13, g_c14, g_c15, g_d1_w, g_rab, g_rbc, g_&
        &rcd, g_rh1o1
    save g_c1, g_c2, g_c3, g_c4, g_c5, g_c6, g_c7, g_c8, g_c9, g_c10
    integer, parameter :: g__pmax_ = 12
    integer            :: g2i, g2p, ldg_g_x, ldg_g_v
    double precision   :: g_g_c1(g__pmax_,g_pmax_), g_g_x(ldg_g_x,ldg_x &
        &,12), g_g_c2(g__pmax_, g_pmax_), g_g_c3(g__pmax_, g_pmax_), g_g_c&
        &4(g__pmax_, g_pmax_), g_g_c5(g__pmax_, g_pmax_), g_g_c6(g__pmax_, &
        & g_pmax_), g_g_c7(g__pmax_, g_pmax_), g_g_c8(g__pmax_, g_pmax_), g_&
        &g_c9(g__pmax_, g_pmax_)
    double precision :: g_g_c10(g__pmax_,g_pmax_), g_g_c11(g__pmax_,g&
        &_pmax_), g_g_c12(g__pmax_, g_pmax_), g_g_c13(g__pmax_, g_pmax_), g&
        &_g_c14(g__pmax_, g_pmax_), g_g_c15(g__pmax_, g_pmax_), g_g_d1_w(g_&
        &_pmax_, g_pmax_), g_g_rab(g__pmax_, g_pmax_), g_g_rbc(g__pmax_, g_&
        &pmax_), g_g_rcd(g__pmax_, g_pmax_)
    double precision :: g_g_rh1o1(g__pmax_,g_pmax_), g_g_ro1o2(g__pmax&
        &_,g_pmax_), g_g_ro2h2(g__pmax_, g_pmax_), g_g_n1(g__pmax_, g_pmax&
        &_), g_g_d1(g__pmax_, g_pmax_), g_g_ah1oo(g__pmax_, g_pmax_), g_g_n&
        &2(g__pmax_, g_pmax_), g_g_d2(g__pmax_, g_pmax_), g_g_aooh2(g__pmax&
        &_, g_pmax_), g_g_xcomp1(g__pmax_, g_pmax_)
    double precision :: g_g_ycomp1(g__pmax_,g_pmax_), g_g_zcomp1(g__pm&
        &ax_,g_pmax_), g_g_xcomp2(g__pmax_, g_pmax_), g_g_ycomp2(g__pmax_, &
        & g_pmax_), g_g_zcomp2(g__pmax_, g_pmax_), g_g_x_normal1(g__pmax_, &
        & g_pmax_), g_g_y_normal1(g__pmax_, g_pmax_), g_g_z_normal1(g__pmax_ &
        &, g_pmax_), g_g_xcomp3(g__pmax_, g_pmax_), g_g_ycomp3(g__pmax_, g_&
        &pmax_)
    double precision :: g_g_zcomp3(g__pmax_,g_pmax_), g_g_xcomp4(g__pm&
        &ax_,g_pmax_), g_g_ycomp4(g__pmax_, g_pmax_), g_g_zcomp4(g__pmax_, &
        & g_pmax_), g_g_x_normal2(g__pmax_, g_pmax_), g_g_y_normal2(g__pmax&
        &_, g_pmax_), g_g_z_normal2(g__pmax_, g_pmax_), g_g_n(g__pmax_, g_p&
        &max_), g_g_d2_w(g__pmax_, g_pmax_), g_g_d(g__pmax_, g_pmax_)
    double precision :: g_g_thooh(g__pmax_,g_pmax_), g_g_pass(g__pmax_ &
        &,g_pmax_,6), g_g_v(ldg_g_v, ldg_v)
    save g_g_pass
    save g_g_xcomp4, g_g_ycomp4, g_g_zcomp4, g_g_x_normal2, g_g_y_no&
        &rmal2, g_g_z_normal2, g_g_n, g_g_d2_w, g_g_d, g_g_thooh
    save g_g_zcomp1, g_g_xcomp2, g_g_ycomp2, g_g_zcomp2, g_g_x_norma&
        &l1, g_g_y_normal1, g_g_z_normal1, g_g_xcomp3, g_g_ycomp3, g_g_zcom&
        &p3
    save g_g_ro1o2, g_g_ro2h2, g_g_n1, g_g_d1, g_g_ah1oo, g_g_n2, g_&
        &g_d2, g_g_aooh2, g_g_xcomp1, g_g_ycomp1
    save g_g_c11, g_g_c12, g_g_c13, g_g_c14, g_g_c15, g_g_d1_w, g_g_&
        &rab, g_g_rbc, g_g_rcd, g_g_rh1o1
    save g_g_c1, g_g_c2, g_g_c3, g_g_c4, g_g_c5, g_g_c6, g_g_c7, g_g&
        &_c8, g_g_c9, g_g_c10
    external g_g_vibpot
    !        data rohlow, rohhigh /1.d0, 4.d0/
    !        data roolow, roohigh, rhhlow /1.8d0, 4.d0, 0.7d0/
    !
    !  Calculate energy at the configuration space point
    !  NOTE :: This is not necessarily a quadrature point
    !
    ! Expecting input in cartesians
    !
    data g_ehfid /0/
    !
    if (g2p .gt. g2pmax) then
        print *, 'Parameter g2p is greater than g2pmax'
        stop
    end if
    if (g_p_ .gt. g_pmax_) then
        print *, 'Parameter g_p_ is greater than g_pmax_'
        stop
    end if
    pi = 4.d0 * atan(1.d0)
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c1(g2i,g_i_) = g_g_x(g2i,g_i_,1) + (-g_g_x(g2i,g_i_,4))
        end do
        g_c1(g_i_) = -g_x(g_i_, 4) + g_x(g_i_, 1)
        !--------
    end do
    c1 = x(1) - x(4)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c2(g2i, g_i_) = g_g_x(g2i,g_i_,2) + (-g_g_x(g2i,g_i_,5))
        end do
        g_c2(g_i_) = -g_x(g_i_, 5) + g_x(g_i_, 2)
        !--------
    end do
    c2 = x(2) - x(5)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c3(g2i, g_i_) = g_g_x(g2i,g_i_,3) + (-g_g_x(g2i,g_i_,6))
        end do
        g_c3(g_i_) = -g_x(g_i_, 6) + g_x(g_i_, 3)
        !--------
    end do
    c3 = x(3) - x(6)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c4(g2i, g_i_) = g_g_x(g2i,g_i_,4) + (-g_g_x(g2i,g_i_,7))
        end do
        g_c4(g_i_) = -g_x(g_i_, 7) + g_x(g_i_, 4)
        !--------
    end do
    c4 = x(4) - x(7)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c5(g2i, g_i_) = g_g_x(g2i,g_i_,5) + (-g_g_x(g2i,g_i_,8))
        end do
        g_c5(g_i_) = -g_x(g_i_, 8) + g_x(g_i_, 5)
        !--------
    end do
    c5 = x(5) - x(8)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c6(g2i, g_i_) = g_g_x(g2i,g_i_,6) + (-g_g_x(g2i,g_i_,9))
        end do
        g_c6(g_i_) = -g_x(g_i_, 9) + g_x(g_i_, 6)
        !--------
    end do
    c6 = x(6) - x(9)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c7(g2i, g_i_) = g_g_x(g2i,g_i_,7) + (-g_g_x(g2i,g_i_,10))
        end do
        g_c7(g_i_) = -g_x(g_i_, 10) + g_x(g_i_, 7)
        !--------
    end do
    c7 = x(7) - x(10)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c8(g2i, g_i_) = g_g_x(g2i,g_i_,8) + (-g_g_x(g2i,g_i_,11))
        end do
        g_c8(g_i_) = -g_x(g_i_, 11) + g_x(g_i_, 8)
        !--------
    end do
    c8 = x(8) - x(11)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c9(g2i, g_i_) = g_g_x(g2i,g_i_,9) + (-g_g_x(g2i,g_i_,12))
        end do
        g_c9(g_i_) = -g_x(g_i_, 12) + g_x(g_i_, 9)
        !--------
    end do
    c9 = x(9) - x(12)
    !--------
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c10(g2i,g_i_) = g_g_x(g2i,g_i_,1) + (-g_g_x(g2i,g_i_,7))
        end do
        g_c10(g_i_) = -g_x(g_i_, 7) + g_x(g_i_, 1)
        !--------
    end do
    c10 = x(1) - x(7)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c11(g2i,g_i_) = g_g_x(g2i,g_i_,2) + (-g_g_x(g2i,g_i_,8))
        end do
        g_c11(g_i_) = -g_x(g_i_, 8) + g_x(g_i_, 2)
        !--------
    end do
    c11 = x(2) - x(8)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c12(g2i, g_i_) = g_g_x(g2i,g_i_,3) + (-g_g_x(g2i,g_i_,9))
        end do
        g_c12(g_i_) = -g_x(g_i_, 9) + g_x(g_i_, 3)
        !--------
    end do
    c12 = x(3) - x(9)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c13(g2i, g_i_) = g_g_x(g2i,g_i_,4) + (-g_g_x(g2i,g_i_,10))
        end do
        g_c13(g_i_) = -g_x(g_i_, 10) + g_x(g_i_, 4)
        !--------
    end do
    c13 = x(4) - x(10)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c14(g2i, g_i_) = g_g_x(g2i,g_i_,5) + (-g_g_x(g2i,g_i_,11))
        end do
        g_c14(g_i_) = -g_x(g_i_, 11) + g_x(g_i_, 5)
        !--------
    end do
    c14 = x(5) - x(11)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_c15(g2i, g_i_) = g_g_x(g2i,g_i_,6) + (-g_g_x(g2i,g_i_,12))
        end do
        g_c15(g_i_) = -g_x(g_i_, 12) + g_x(g_i_, 6)
        !--------
    end do
    c15 = x(6) - x(12)
    !--------
    c16 = x(1) - x(10)
    c17 = x(2) - x(11)
    c18 = x(3) - x(12)
    !
    !
    ! H1 -- O1 -- O2 -- H2
    !   rH1O1 rO1O2 rO2H2
    !
    !  angle H1O1O2 = AH1OO
    !  angle O1O2H2 = AOOH2
    !
    !  TORSION HOOH = THOOH
    !
    d2_v = c1 * c1
    d3_p = 2.0d0 * c1
    d4_v = c2 * c2
    d2_p = 2.0d0 * c2
    d7_v = c3 * c3
    d1_p = 2.0d0 * c3
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d1_w(g2i, g_i_) = d3_p * g_g_c1(g2i,g_i_) + d2_p * g_g_c2(g2i, g_i_) + d1_p * g_g_c3(g2i, g_i_)
        end do
        g_d1_w(g_i_) = d1_p * g_c3(g_i_) + d2_p * g_c2(g_i_) + d3_p * g_c1(g_i_)
        !--------
    end do
    d1_w = d2_v + d4_v + d7_v
    d2_v = sqrt(d1_w)
    !
    d1_p = 1.0d0 / (2.0d0 * d2_v)
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_rab(g2i, g_i_) = d1_p * g_g_d1_w(g2i, g_i_)
        end do
        g_rab(g_i_) = d1_p * g_d1_w(g_i_)
        !--------
    end do
    rab = d2_v
    !--------
    d2_v = c4 * c4
    d3_p = 2.0d0 * c4
    d4_v = c5 * c5
    d2_p = 2.0d0 * c5
    d7_v = c6 * c6
    d1_p = 2.0d0 * c6
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d1_w(g2i, g_i_) = d3_p * g_g_c4(g2i, g_i_) + d2_p * g_g_c5(g2i, g_i_) + d1_p * g_g_c6(g2i, g_i_)
        end do
        g_d1_w(g_i_) = d1_p * g_c6(g_i_) + d2_p * g_c5(g_i_) + d3_p * g_c4(g_i_)
        !--------
    end do
    d1_w = d2_v + d4_v + d7_v
    d2_v = sqrt(d1_w)
    !
    d1_p = 1.0d0 / (2.0d0 * d2_v)
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_rbc(g2i, g_i_) = d1_p * g_g_d1_w(g2i, g_i_)
        end do
        g_rbc(g_i_) = d1_p * g_d1_w(g_i_)
        !--------
    end do
    rbc = d2_v
    !--------
    d2_v = c7 * c7
    d3_p = 2.0d0 * c7
    d4_v = c8 * c8
    d2_p = 2.0d0 * c8
    d7_v = c9 * c9
    d1_p = 2.0d0 * c9
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d1_w(g2i, g_i_) = d3_p * g_g_c7(g2i, g_i_) + d2_p * g_g_c8(g2i, g_i_) + d1_p * g_g_c9(g2i, g_i_)
        end do
        g_d1_w(g_i_) = d1_p * g_c9(g_i_) + d2_p * g_c8(g_i_) + d3_p * g_c7(g_i_)
        !--------
    end do
    d1_w = d2_v + d4_v + d7_v
    d2_v = sqrt(d1_w)
    !
    d1_p = 1.0d0 / (2.0d0 * d2_v)
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_rcd(g2i, g_i_) = d1_p * g_g_d1_w(g2i, g_i_)
        end do
        g_rcd(g_i_) = d1_p * g_d1_w(g_i_)
        !--------
    end do
    rcd = d2_v
    !--------
    rac = sqrt(c10 ** 2 + c11 ** 2 + c12 ** 2)
    rbd = sqrt(c13 ** 2 + c14 ** 2 + c15 ** 2)
    rad = sqrt(c16 ** 2 + c17 ** 2 + c18 ** 2)
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_rh1o1(g2i,g_i_) = g_g_rab(g2i, g_i_)
        end do
        g_rh1o1(g_i_) = g_rab(g_i_)
        !--------
    end do
    rh1o1 = rab
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_ro1o2(g2i,g_i_) = g_g_rbc(g2i, g_i_)
        end do
        g_ro1o2(g_i_) = g_rbc(g_i_)
        !--------
    end do
    ro1o2 = rbc
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_ro2h2(g2i, g_i_) = g_g_rcd(g2i, g_i_)
        end do
        g_ro2h2(g_i_) = g_rcd(g_i_)
        !--------
    end do
    ro2h2 = rcd
    !--------
    !
    ! The potential can give weird values for unusual geometries
    ! replace with some high value in these unphysical regions and
    ! then exit
    !
    !
    !      if((min(rab,rcd,rac,rbd).gt.rohhigh) .or.
    !     !(min(rab,rcd,rac,rbd).lt.rohlow)  .or.
    !     !(rbc.gt.roohigh) .or. (rbc.lt.roolow) .or.
    !     !(rad.lt.rhhlow)) then
    !      vdc=1000.d0
    !      return
    !      endif
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_n1(g2i, g_i_) = (-c1) * g_g_c4(g2i, g_i_) + (-c4) * g_g_c1(g2i,g_i_) + (-c2) * g_g_c5(g2i, g_i_) + (- &
                & c5) * g_g_c2 (g2i, g_i_) + (-c3) * g_g_c6(g2i, g_i_) + (-c6) * g_g_c3(g2i, g_i_)
        end do
        g_n1(g_i_) = (-c6) * g_c3(g_i_) + (-c3) * g_c6(g_i_) + (-c5) * g_c2(g_i_) + (-c2) * g_c5(g_i_) + (-c4) * &
            & g_c1(g_i_) + (-c1) * g_c4(g_i_)
        !--------
    end do
    n1 = (-c4) * c1 + (-c5) * c2 + (-c6) * c3
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d1(g2i, g_i_) = rh1o1 * g_g_ro1o2(g2i,g_i_) + ro1o2 * g_g_rh1o1(g2i,g_i_)
        end do
        g_d1(g_i_) = ro1o2 * g_rh1o1(g_i_) + rh1o1 * g_ro1o2(g_i_)
        !--------
    end do
    d1 = ro1o2 * rh1o1
    !--------
    d3_v = n1 / d1
    d2_b = 1.0d0 / d1
    d3_b = (-d3_v) / d1
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d1_w(g2i, g_i_) = d2_b * g_g_n1(g2i, g_i_) + d3_b * g_g_d1(g2i, g_i_)
        end do
        g_d1_w(g_i_) = d3_b * g_d1(g_i_) + d2_b * g_n1(g_i_)
        !--------
    end do
    d1_w = d3_v
    d2_v = acos(d1_w)
    !
    d1_p = (-1.0d0) / sqrt((1.0d0 - d1_w) * (1.0d0 + d1_w))
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_ah1oo(g2i, g_i_) = d1_p * g_g_d1_w(g2i, g_i_)
        end do
        g_ah1oo(g_i_) = d1_p * g_d1_w(g_i_)
        !--------
    end do
    ah1oo = d2_v
    !--------
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_n2(g2i, g_i_) = (-c4) * g_g_c7(g2i, g_i_) + (-c7) * g_g_c4(g2i, g_i_) + (-c5) * g_g_c8(g2i, g_i_) + (- &
                & c8) * g_g_c5 (g2i, g_i_) + (-c6) * g_g_c9(g2i, g_i_) + (-c9) * g_g_c6(g2i, g_i_)
        end do
        g_n2(g_i_) = (-c9) * g_c6(g_i_) + (-c6) * g_c9(g_i_) + (-c8) * g_c5(g_i_) + (-c5) * g_c8(g_i_) + (-c7) * &
            & g_c4(g_i_) + (-c4) * g_c7(g_i_)
        !--------
    end do
    n2 = -(c7 * c4 + c8 * c5 + c9 * c6)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d2(g2i, g_i_) = ro1o2 * g_g_ro2h2(g2i, g_i_) + ro2h2 * g_g_ro1o2(g2i,g_i_)
        end do
        g_d2(g_i_) = ro2h2 * g_ro1o2(g_i_) + ro1o2 * g_ro2h2(g_i_)
        !--------
    end do
    d2 = ro2h2 * ro1o2
    !--------
    d3_v = n2 / d2
    d2_b = 1.0d0 / d2
    d3_b = (-d3_v) / d2
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d1_w(g2i, g_i_) = d2_b * g_g_n2(g2i, g_i_) + d3_b * g_g_d2(g2i, g_i_)
        end do
        g_d1_w(g_i_) = d3_b * g_d2(g_i_) + d2_b * g_n2(g_i_)
        !--------
    end do
    d1_w = d3_v
    d2_v = acos(d1_w)
    !
    d1_p = (-1.0d0) / sqrt((1.0d0 - d1_w) * (1.0d0 + d1_w))
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_aooh2(g2i, g_i_) = d1_p * g_g_d1_w(g2i, g_i_)
        end do
        g_aooh2(g_i_) = d1_p * g_d1_w(g_i_)
        !--------
    end do
    aooh2 = d2_v
    !--------
    !
    ! Calculate torsion angle HOOH
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_xcomp1(g2i, g_i_) = g_g_x(g2i,g_i_,4) + (-g_g_x( g2i,g_i_,1))
        end do
        g_xcomp1(g_i_) = -g_x(g_i_, 1) + g_x(g_i_, 4)
        !--------
    end do
    xcomp1 = x(4) - x(1)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_ycomp1(g2i,g_i_) = g_g_x(g2i,g_i_,5) + (-g_g_x( g2i,g_i_,2))
        end do
        g_ycomp1(g_i_) = -g_x(g_i_, 2) + g_x(g_i_, 5)
        !--------
    end do
    ycomp1 = x(5) - x(2)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_zcomp1(g2i,g_i_) = g_g_x(g2i,g_i_,6) + (-g_g_x( g2i,g_i_,3))
        end do
        g_zcomp1(g_i_) = -g_x(g_i_, 3) + g_x(g_i_, 6)
        !--------
    end do
    zcomp1 = x(6) - x(3)
    !--------
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_xcomp2(g2i, g_i_) = g_g_c10(g2i,g_i_)
        end do
        g_xcomp2(g_i_) = g_c10(g_i_)
        !--------
    end do
    xcomp2 = c10
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_ycomp2(g2i, g_i_) = g_g_c11(g2i,g_i_)
        end do
        g_ycomp2(g_i_) = g_c11(g_i_)
        !--------
    end do
    ycomp2 = c11
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_zcomp2(g2i, g_i_) = g_g_c12(g2i, g_i_)
        end do
        g_zcomp2(g_i_) = g_c12(g_i_)
        !--------
    end do
    zcomp2 = c12
    !--------
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_x_normal1(g2i, g_i_) = zcomp2 * g_g_ycomp1(g2i,g_i_ ) + ycomp1 * g_g_zcomp2(g2i, g_i_) + (-ycomp2) * &
                & g_g_zcomp1(g2i,g_i_) + (-zcomp1) * g_g_ycomp2(g2i, g_i_)
        end do
        g_x_normal1(g_i_) = (-zcomp1) * g_ycomp2(g_i_) + (-ycomp2) * g_zcomp1(g_i_) + ycomp1 * g_zcomp2(g_i_) + &
            & zcomp2 * g_ycomp1(g_i_)
        !--------
    end do
    x_normal1 = ycomp1 * zcomp2 - zcomp1 * ycomp2
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_y_normal1(g2i, g_i_) = xcomp2 * g_g_zcomp1(g2i,g_i_ ) + zcomp1 * g_g_xcomp2(g2i, g_i_) + (-zcomp2) * &
                & g_g_xcomp1(g2i, g_i_) + (-xcomp1) * g_g_zcomp2(g2i, g_i_)
        end do
        g_y_normal1(g_i_) = (-xcomp1) * g_zcomp2(g_i_) + (-zcomp2) * g_xcomp1(g_i_) + zcomp1 * g_xcomp2(g_i_) + &
            & xcomp2 * g_zcomp1(g_i_)
        !--------
    end do
    y_normal1 = zcomp1 * xcomp2 - xcomp1 * zcomp2
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_z_normal1(g2i, g_i_) = ycomp2 * g_g_xcomp1(g2i, g_i_ ) + xcomp1 * g_g_ycomp2(g2i, g_i_) + (-xcomp2) * &
                & g_g_ycomp1(g2i,g_i_) + (-ycomp1) * g_g_xcomp2(g2i, g_i_)
        end do
        g_z_normal1(g_i_) = (-ycomp1) * g_xcomp2(g_i_) + (-xcomp2) * g_ycomp1(g_i_) + xcomp1 * g_ycomp2(g_i_) + &
            & ycomp2 * g_xcomp1(g_i_)
        !--------
    end do
    z_normal1 = xcomp1 * ycomp2 - ycomp1 * xcomp2
    !--------
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_xcomp3(g2i, g_i_) = g_g_x(g2i,g_i_,7) + (-g_g_x( g2i,g_i_,4))
        end do
        g_xcomp3(g_i_) = -g_x(g_i_, 4) + g_x(g_i_, 7)
        !--------
    end do
    xcomp3 = x(7) - x(4)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_ycomp3(g2i, g_i_) = g_g_x(g2i,g_i_,8) + (-g_g_x( g2i,g_i_,5))
        end do
        g_ycomp3(g_i_) = -g_x(g_i_, 5) + g_x(g_i_, 8)
        !--------
    end do
    ycomp3 = x(8) - x(5)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_zcomp3(g2i,g_i_) = g_g_x(g2i,g_i_,9) + (-g_g_x( g2i,g_i_,6))
        end do
        g_zcomp3(g_i_) = -g_x(g_i_, 6) + g_x(g_i_, 9)
        !--------
    end do
    zcomp3 = x(9) - x(6)
    !--------
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_xcomp4(g2i,g_i_) = g_g_c13(g2i, g_i_)
        end do
        g_xcomp4(g_i_) = g_c13(g_i_)
        !--------
    end do
    xcomp4 = c13
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_ycomp4(g2i, g_i_) = g_g_c14(g2i, g_i_)
        end do
        g_ycomp4(g_i_) = g_c14(g_i_)
        !--------
    end do
    ycomp4 = c14
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_zcomp4(g2i, g_i_) = g_g_c15(g2i, g_i_)
        end do
        g_zcomp4(g_i_) = g_c15(g_i_)
        !--------
    end do
    zcomp4 = c15
    !--------
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_x_normal2(g2i, g_i_) = zcomp4 * g_g_ycomp3(g2i, g_i_ ) + ycomp3 * g_g_zcomp4(g2i, g_i_) + (-ycomp4) * &
                & g_g_zcomp3(g2i,g_i_) + (-zcomp3) * g_g_ycomp4(g2i, g_i_)
        end do
        g_x_normal2(g_i_) = (-zcomp3) * g_ycomp4(g_i_) + (-ycomp4) * g_zcomp3(g_i_) + ycomp3 * g_zcomp4(g_i_) + &
            & zcomp4 * g_ycomp3(g_i_)
        !--------
    end do
    x_normal2 = ycomp3 * zcomp4 - zcomp3 * ycomp4
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_y_normal2(g2i, g_i_) = xcomp4 * g_g_zcomp3(g2i,g_i_ ) + zcomp3 * g_g_xcomp4(g2i,g_i_) + (-zcomp4) * &
                & g_g_xcomp3(g2i, g_i_) + (-xcomp3) * g_g_zcomp4(g2i, g_i_)
        end do
        g_y_normal2(g_i_) = (-xcomp3) * g_zcomp4(g_i_) + (-zcomp4) * g_xcomp3(g_i_) + zcomp3 * g_xcomp4(g_i_) + &
            & xcomp4 * g_zcomp3(g_i_)
        !--------
    end do
    y_normal2 = zcomp3 * xcomp4 - xcomp3 * zcomp4
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_z_normal2(g2i, g_i_) = ycomp4 * g_g_xcomp3(g2i, g_i_ ) + xcomp3 * g_g_ycomp4(g2i, g_i_) + (-xcomp4) * &
                & g_g_ycomp3(g2i, g_i_) + (-ycomp3) * g_g_xcomp4(g2i,g_i_)
        end do
        g_z_normal2(g_i_) = (-ycomp3) * g_xcomp4(g_i_) + (-xcomp4) * g_ycomp3(g_i_) + xcomp3 * g_ycomp4(g_i_) + &
            & ycomp4 * g_xcomp3(g_i_)
        !--------
    end do
    z_normal2 = xcomp3 * ycomp4 - ycomp3 * xcomp4
    !--------
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_n(g2i, g_i_) = x_normal2 * g_g_x_normal1(g2i, g_i_) + x_normal1 * g_g_x_normal2(g2i, g_i_) + &
                & y_normal2 * g_g_y_normal1(g2i, g_i_) + y_normal1 * g_g_y_normal2(g2i, g_i_) + z_normal2 * &
                & g_g_z_normal1(g2i, g_i_) + z_normal1 * g_g_z_normal2(g2i, g_i_)
        end do
        g_n(g_i_) = z_normal1 * g_z_normal2(g_i_) + z_normal2 * g_z_normal1(g_i_) + y_normal1 * g_y_normal2(g_i_) + &
            & y_normal2 * g_y_normal1(g_i_) + x_normal1 * g_x_normal2(g_i_) + x_normal2 * g_x_normal1(g_i_)
        !--------
    end do
    n = x_normal1 * x_normal2 + y_normal1 * y_normal2 + z_normal1 * z_normal2
    !--------
    !
    d4_b = z_normal1 + z_normal1
    d7_b = y_normal1 + y_normal1
    d8_b = x_normal1 + x_normal1
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d1_w(g2i, g_i_) = d8_b * g_g_x_normal1(g2i, g_i_) + d7_b * g_g_y_normal1(g2i, g_i_) + d4_b * &
                & g_g_z_normal1(g2i, g_i_)
        end do
        g_d1_w(g_i_) = d4_b * g_z_normal1(g_i_) + d7_b * g_y_normal1(g_i_) + d8_b * g_x_normal1(g_i_)
        !--------
    end do
    d1_w = x_normal1 * x_normal1 + y_normal1 * y_normal1 + z_normal1 * z_normal1
    d4_b = z_normal2 + z_normal2
    d7_b = y_normal2 + y_normal2
    d8_b = x_normal2 + x_normal2
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d2_w(g2i, g_i_) = d8_b * g_g_x_normal2(g2i, g_i_) + d7_b * g_g_y_normal2(g2i, g_i_) + d4_b * &
                & g_g_z_normal2(g2i, g_i_)
        end do
        g_d2_w(g_i_) = d4_b * g_z_normal2(g_i_) + d7_b * g_y_normal2(g_i_) + d8_b * g_x_normal2(g_i_)
        !--------
    end do
    d2_w = x_normal2 * x_normal2 + y_normal2 * y_normal2 + z_normal2 * z_normal2
    d2_v = sqrt(d1_w)
    !
    d2_p = 1.0d0 / (2.0d0 * d2_v)
    d4_v = sqrt(d2_w)
    !
    d1_p = 1.0d0 / (2.0d0 * d4_v)
    d4_b = d2_v * d1_p
    d5_b = d4_v * d2_p
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d(g2i, g_i_) = d5_b * g_g_d1_w(g2i, g_i_) + d4_b * g_g_d2_w(g2i, g_i_)
        end do
        g_d(g_i_) = d4_b * g_d2_w(g_i_) + d5_b * g_d1_w(g_i_)
        !--------
    end do
    d = d2_v * d4_v
    !--------
    !
    !  We can occasionally get "NaNs" if we don't trap for special cases whi
    !lating
    !  the dihedral angle.  Acos wants an argument in [-1,1], but slight rou
    !n push us
    !  outside this range.  Also if certain triatoms are collinear then the
    ! is not defined
    !  and we set these cases to 0 degrees
    !
    !      if(d.eq.0.d0)then
    !       thooh=0.d0
    !      else
    !       argnd=n/d
    !      if(argnd.ge.1.d0)then
    !       thooh=0.d0
    !      elseif(argnd.le.-1.d0)then
    !       thooh=pi
    !      else
    d3_v = n / d
    d2_b = 1.0d0 / d
    d3_b = (-d3_v) / d
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_d1_w(g2i, g_i_) = d2_b * g_g_n(g2i, g_i_) + d3_b * g_g_d(g2i, g_i_)
        end do
        g_d1_w(g_i_) = d3_b * g_d(g_i_) + d2_b * g_n(g_i_)
        !--------
    end do
    d1_w = d3_v
    d2_v = acos(d1_w)
    !
    d1_p = (-1.0d0) / sqrt((1.0d0 - d1_w) * (1.0d0 + d1_w))
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_thooh(g2i,g_i_) = d1_p * g_g_d1_w(g2i, g_i_)
        end do
        g_thooh(g_i_) = d1_p * g_d1_w(g_i_)
        !--------
    end do
    thooh = d2_v
    !--------
    !      endif
    !      endif
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_pass(g2i,g_i_,1) = g_g_rh1o1(g2i,g_i_)
        end do
        g_pass(g_i_, 1) = g_rh1o1(g_i_)
        !--------
    end do
    pass(1) = rh1o1
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_pass(g2i,g_i_,2) = g_g_ro1o2(g2i,g_i_)
        end do
        g_pass(g_i_, 2) = g_ro1o2(g_i_)
        !--------
    end do
    pass(2) = ro1o2
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_pass(g2i,g_i_,3) = g_g_ro2h2(g2i, g_i_)
        end do
        g_pass(g_i_, 3) = g_ro2h2(g_i_)
        !--------
    end do
    pass(3) = ro2h2
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_pass(g2i,g_i_,4) = g_g_ah1oo(g2i, g_i_)
        end do
        g_pass(g_i_, 4) = g_ah1oo(g_i_)
        !--------
    end do
    pass(4) = ah1oo
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_pass(g2i,g_i_,5) = g_g_aooh2(g2i, g_i_)
        end do
        g_pass(g_i_, 5) = g_aooh2(g_i_)
        !--------
    end do
    pass(5) = aooh2
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_pass(g2i,g_i_,6) = g_g_thooh(g2i,g_i_)
        end do
        g_pass(g_i_, 6) = g_thooh(g_i_)
        !--------
    end do
    pass(6) = thooh
    !--------
    !
    !
    ! Call Koput, Carter, and Handy PES
    !
    !
    call g_g_vibpot(g2p, g_p_, pass, g_pass, g_g_pass, g2pmax, g&
        &_pmax_, v, g_v, g_g_v, ldg_g_v, ldg_v)
    !
    !
end subroutine g_analytic

!
!***********************************************************************
!
!> @brief Evaluate g_g_vibpot for the H2O2_GEN_1998 potential.
subroutine g_g_vibpot( g2p, g_p_, coord, g_coord, g_g_coord, ldg_&
        &g_coord, ldg_coord, v, g_v, g_g_v, ldg_g_v, ldg_v )
    !
    implicit none
    double precision   :: v, pi
    double precision   :: coord(6), q1p(0:4), q2p(0:4), q3p(0:4), q4p(0:4)
    double precision   :: q5p(0:4), q6p(0:4)
    double precision   :: ind1(152), ind2(152), ind3(152), ind4(152)
    double precision   :: ind5(152), ind6(152), coeff(152)
    double precision   :: rm_oo, rm_oh, angle_ooh, r1, r2, r3
    double precision   :: a1, a2, q1, q2, q3, q4, q5, q6
    integer            :: i
    double precision   :: zoe
    integer, parameter :: g_pmax_ = 12
    integer            :: g_i_, g_p_, ldg_coord, ldg_v
    double precision   :: d14_b, d13_b, d7_v, d11_b, d10_b, d9_b, d8_b, d&
        &7_b, d2_v, d3_v
    double precision :: d6_b, d2_b, d3_b, d1_w, d1_p, d11_v, d5_v, d9_v &
        &, g_r1(g_pmax_), g_coord(ldg_coord, 6)
    double precision :: g_r2(g_pmax_), g_r3(g_pmax_), g_a1(g_pmax_), g_&
        &a2(g_pmax_), g_q6(g_pmax_), g_q3(g_pmax_), g_q1(g_pmax_), g_q2(g_p&
        &max_), g_q4(g_pmax_), g_q5(g_pmax_)
    double precision :: g_q1p(g_pmax_,0:4), g_q2p(g_pmax_,0:4), g_q3p &
        & (g_pmax_, 0:4), g_q4p(g_pmax_, 0:4), g_q5p(g_pmax_, 0:4), g_q6p(g_&
        &pmax_, 0:4), g_d1_w(g_pmax_), g_v(ldg_v)
    save g_q5, g_q1p, g_q2p, g_q3p, g_q4p, g_q5p, g_q6p, g_d1_w
    save g_r1, g_r2, g_r3, g_a1, g_a2, g_q6, g_q3, g_q1, g_q2, g_q4
    intrinsic dble
    integer            :: g_ehfid
    integer, parameter :: g__pmax_ = 12
    integer            :: g2i, g2p, ldg_g_coord, ldg_g_v
    double precision   :: g_g_r1(g__pmax_,g_pmax_), g_g_coord(ldg_g_coor&
        &d,ldg_coord,6), g_g_r2(g__pmax_, g_pmax_), g_g_r3(g__pmax_, g_pm&
        &ax_), g_g_a1(g__pmax_, g_pmax_), g_g_a2(g__pmax_, g_pmax_), g_g_q6 &
        & (g__pmax_, g_pmax_), g_g_q3(g__pmax_, g_pmax_), g_g_q1(g__pmax_, g&
        &_pmax_), g_g_q2(g__pmax_, g_pmax_)
    double precision :: g_g_q4(g__pmax_,g_pmax_), g_g_q5(g__pmax_,g_p&
        &max_), g_g_q1p(g__pmax_, g_pmax_, 0:4), g_g_q2p(g__pmax_, g_pmax_, &
        & 0:4), g_g_q3p(g__pmax_, g_pmax_, 0:4), g_g_q4p(g__pmax_, g_pmax_, &
        & 0:4), g_g_q5p(g__pmax_, g_pmax_, 0:4), g_g_q6p(g__pmax_, g_pmax_, &
        & 0:4), g_g_d1_w(g__pmax_, g_pmax_), g_g_v(ldg_g_v, ldg_v)
    save g_g_q5, g_g_q1p, g_g_q2p, g_g_q3p, g_g_q4p, g_g_q5p, g_g_q6&
        &p, g_g_d1_w
    save g_g_r1, g_g_r2, g_g_r3, g_g_a1, g_g_a2, g_g_q6, g_g_q3, g_g&
        &_q1, g_g_q2, g_g_q4
    data zoe /3.995139982058d-3/
    !
    data ind1 /0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, &
        & 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 2, 1, 0, 2, 0, 1, 0, &
        & 2, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 4, 0, 0, 0, 2, 0, 0, &
        & 0, 2, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 3, 1, 0, 3, 0, 0, 0, 1, 1, &
        & 0, 1, 0, 0, 1, 0, 2, 0, 2, 0, 0, 0, 1, 1, 2, 0, 1, 0, 1, 0, 0, 1, &
        & 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 3, 0, 0, 0, 2, &
        & 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, &
        & 0, 0/
    !
    data ind2 /0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, &
        & 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 2, 1, 0, 1, 0, 2, 0, 1, &
        & 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 4, 0, 0, 0, 2, 0, &
        & 0, 0, 2, 0, 0, 1, 0, 0, 0, 3, 0, 0, 3, 1, 0, 1, 0, 3, 0, 0, 1, 0, &
        & 1, 0, 1, 0, 0, 1, 0, 2, 0, 2, 0, 0, 1, 1, 0, 2, 0, 1, 0, 1, 0, 0, &
        & 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 3, 0, 0, 0, &
        & 2, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, &
        & 0, 0/
    !
    data ind3 /0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, &
        & 0, 3, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 4, 0, 0, 0, 0, 2, 2, 2, &
        & 2, 0, 0, 0, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, &
        & 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, &
        & 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, &
        & 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, &
        & 0, 0/
    !
    data ind4 /0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, &
        & 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, 0, 0, 2, &
        & 0, 1, 1, 2, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 4, 0, 0, 0, 2, &
        & 0, 2, 0, 2, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 3, 0, 1, 0, 1, 3, 0, 1, &
        & 0, 0, 1, 1, 2, 0, 1, 0, 0, 1, 1, 2, 2, 0, 1, 1, 1, 2, 2, 1, 0, 0, &
        & 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 3, 0, 0, &
        & 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, &
        & 1, 0/
    !
    data ind5 /0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 1, 0, &
        & 1, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, 2, 0, &
        & 1, 0, 2, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 4, 0, 0, 0, &
        & 2, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 3, 0, 1, 3, 1, 0, 0, &
        & 1, 1, 0, 1, 0, 2, 0, 1, 1, 0, 2, 1, 0, 2, 1, 1, 2, 1, 1, 2, 0, 0, &
        & 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 3, 0, &
        & 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, &
        & 0, 1/
    !
    data ind6 /1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, &
        & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
        & 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, &
        & 3, 3/
    !
    data coeff /0.00482409d0, 0.00325127d0, 0.00026630d0, 0.00004135&
        &d0, 1.08154998d0, 0.86097723d0, 0.86097723d0, 0.11016112d0, 0.1101&
        &6112d0, -0.03637740d0, -0.03637740d0, 0.17491854d0, 0.17491854d0, &
        & 0.00057054d0, -0.00137967d0, -0.00137967d0, 0.00062152d0, 0.000621&
        &52d0, 0.01375726d0, -1.15691421d0, -0.25495918d0, -0.25495918d0, -&
        &0.02272830d0, -0.02272830d0, -0.18381415d0, -0.18381415d0, -0.3462&
        &7237d0, -0.34627237d0, 0.13974588d0, 0.13974588d0, -0.27290455d0, &
        & -0.27290455d0, -0.00674721d0, -0.00674721d0, -0.02179545d0, -0.021&
        &79545d0, -0.02125643d0, -0.02125643d0, -0.00491968d0, -0.00491968d&
        &0, 0.00233773d0, 0.00233773d0, -0.00050066d0, -0.00050066d0, 0.018&
        &17536d0, -0.04666009d0, -0.04666009d0, -0.02424748d0, -0.02424748d&
        &0, -0.01727148d0, -0.00420506d0, -0.00420506d0, -0.00647944d0, -0.&
        &00647944d0, -1.06749007d0, -0.35741007d0, -0.35741007d0, -0.007968&
        &36d0, -0.00796836d0, -0.42556742d0, -0.42556742d0, 0.06278896d0, 0&
        &.06278896d0, -0.04010419d0, -0.04010419d0, -0.00993912d0, 0.475628&
        &94d0, 0.47562894d0, -0.40830627d0, -0.40830627d0, 0.22073222d0, 0.&
        &22073222d0, 0.07828212d0, 0.07828212d0, -0.02954687d0, -0.02954687&
        &d0, 0.03057888d0, 0.03057888d0, -0.06363999d0, -0.06363999d0, -0.0&
        &0373964d0, -0.00373964d0, -0.04114668d0, 0.11249614d0, 0.11249614d&
        &0, 0.02616679d0, 0.02616679d0, -0.07824425d0, 0.04266205d0, 0.0426&
        &6205d0, -0.07420432d0, -0.07420432d0, -0.08251268d0, -0.08251268d0 &
        &, 0.00270940d0, 0.00270940d0, 0.00199953d0, 0.00199953d0, -0.01292&
        &325d0, -0.01292325d0, -0.02074323d0, -0.02074323d0, -0.00789732d0, &
        & -0.00789732d0, -0.01435326d0, -0.00180710d0, -0.00180710d0, -0.01&
        &135671d0, -0.01135671d0, 0.00020655d0, -0.00492533d0, -0.00492533d&
        &0, 0.00270990d0, 0.00270990d0, 0.00376086d0, 0.00376086d0, 0.00044&
        &732d0, 0.00044732d0, 0.00569979d0, -0.00244774d0, -0.00244774d0, -&
        &0.02065564d0, 0.05249331d0, -0.02490299d0, -0.02490299d0, 0.003914&
        &60d0, 0.00391460d0, 0.08893744d0, 0.08893744d0, -0.01051618d0, 0.0&
        &0120479d0, 0.00120479d0, -0.00111888d0, -0.00111888d0, 0.00884757d&
        &0, 0.00416289d0, 0.00416289d0, 0.00126763d0, 0.00126763d0, -0.0070&
        &6563d0, -0.00706563d0, -0.00840146d0, -0.00840146d0, -0.00139219d0 &
        &, 0.00801673d0, 0.00801673d0, 0.00463860d0, -0.00096051d0, 0.00019&
        &906d0, 0.00019906d0, -0.00057576d0, -0.00057576d0/
    !
    !      data rm_oo, rm_oh, ANGLE_OOH/1.456199d0, 0.962755d0, 100.9059d0/
    !
    !   convert offset parameters to bohr and radians
    !
    data g_ehfid /0/
    !
    if (g2p .gt. g__pmax_) then
        print *, 'Parameter g2p is greater than g__pmax_'
        stop
    end if
    if (g_p_ .gt. g_pmax_) then
        print *, 'Parameter g_p_ is greater than g_pmax_'
        stop
    end if
    rm_oo     = 1.456199d0
    rm_oh     = 0.962755d0
    angle_ooh = 100.9059d0
    !
    pi = 4.d0 * atan(1.d0)
    !        rm_oo = rm_oo / 0.529177d0
    !        rm_oh = rm_oh / 0.529177d0
    angle_ooh = angle_ooh * pi / 180.d0
    !
    !
    ! Potential energy surface for H2O2
    ! Jacek Koput, Stuart Carter, and Nicholas Handy, J. Phys. Chem. A
    ! 1998, volume 102, pages 6325-6330
    !
    ! V is the potential energy in hartrees
    ! minimum configuration is at R(OH)=0.96265, R(OO)=1.45248,
    ! theta(OOH)=99.909, and dihed angle=112.456 degrees
    !
    !      R1 = COORD(1)                ! radius from H1 to O1 (in bohr)
    !      R2 = COORD(2)                ! radius from O1 to O2 (in bohr)
    !      R3 = COORD(3)                ! radius from O2 to H2 (in bohr)
    !      A1 = COORD(4)                ! angle between H1-O1-O2 (in radians
    !      A2 = COORD(5)                ! angle between O1-O2-H2 (in radians
    !      Q6 = COORD(6)                ! HOOH torsion angle (in radians)
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_r1(g2i,g_i_) = g_g_coord(g2i,g_i_,1)
        end do
        g_r1(g_i_) = g_coord(g_i_, 1)
        !--------
    end do
    r1 = coord(1)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_r2(g2i, g_i_) = g_g_coord(g2i,g_i_,2)
        end do
        g_r2(g_i_) = g_coord(g_i_, 2)
        !--------
    end do
    r2 = coord(2)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_r3(g2i, g_i_) = g_g_coord(g2i,g_i_,3)
        end do
        g_r3(g_i_) = g_coord(g_i_, 3)
        !--------
    end do
    r3 = coord(3)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_a1(g2i, g_i_) = g_g_coord(g2i,g_i_,4)
        end do
        g_a1(g_i_) = g_coord(g_i_, 4)
        !--------
    end do
    a1 = coord(4)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_a2(g2i, g_i_) = g_g_coord(g2i,g_i_,5)
        end do
        g_a2(g_i_) = g_coord(g_i_, 5)
        !--------
    end do
    a2 = coord(5)
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_q6(g2i, g_i_) = g_g_coord(g2i,g_i_,6)
        end do
        g_q6(g_i_) = g_coord(g_i_, 6)
        !--------
    end do
    q6 = coord(6)
    !--------
    !
    ! Q1,Q2, and Q3 represent the stretching modes
    ! Q4 and Q5  are the bending modes
    ! Q6 is the torsional mode
    ! Q1, Q2, Q3 are dimensionless, Q4, Q5, and Q6 are in radians
    !
    d3_v = (r2 - rm_oo) / r2
    d3_b = (-d3_v) / r2 + 1.0d0 / r2
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_q3(g2i, g_i_) = d3_b * g_g_r2(g2i, g_i_)
        end do
        g_q3(g_i_) = d3_b * g_r2(g_i_)
        !--------
    end do
    q3 = d3_v
    !--------
    d3_v = (r1 - rm_oh) / r1
    d3_b = (-d3_v) / r1 + 1.0d0 / r1
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_q1(g2i, g_i_) = d3_b * g_g_r1(g2i,g_i_)
        end do
        g_q1(g_i_) = d3_b * g_r1(g_i_)
        !--------
    end do
    q1 = d3_v
    !--------
    d3_v = (r3 - rm_oh) / r3
    d3_b = (-d3_v) / r3 + 1.0d0 / r3
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_q2(g2i, g_i_) = d3_b * g_g_r3(g2i, g_i_)
        end do
        g_q2(g_i_) = d3_b * g_r3(g_i_)
        !--------
    end do
    q2 = d3_v
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_q4(g2i,g_i_) = g_g_a1(g2i, g_i_)
        end do
        g_q4(g_i_) = g_a1(g_i_)
        !--------
    end do
    q4 = a1 - angle_ooh
    !--------
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_q5(g2i,g_i_) = g_g_a2(g2i, g_i_)
        end do
        g_q5(g_i_) = g_a2(g_i_)
        !--------
    end do
    q5 = a2 - angle_ooh
    !--------
    !
    !
    !
    do i = 0, 4
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q1p(g2i, g_i_, i) = 0.0d0
            end do
            g_q1p(g_i_,i) = 0.0d0
            !--------
        end do
        q1p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q2p(g2i, g_i_, i) = 0.0d0
            end do
            g_q2p(g_i_,i) = 0.0d0
            !--------
        end do
        q2p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q3p(g2i, g_i_, i) = 0.0d0
            end do
            g_q3p(g_i_, i) = 0.0d0
            !--------
        end do
        q3p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q4p(g2i, g_i_, i) = 0.0d0
            end do
            g_q4p(g_i_, i) = 0.0d0
            !--------
        end do
        q4p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q5p(g2i, g_i_, i) = 0.0d0
            end do
            g_q5p(g_i_, i) = 0.0d0
            !--------
        end do
        q5p(i) = 1.d0
        !--------
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q6p(g2i, g_i_, i) = 0.0d0
            end do
            g_q6p(g_i_, i) = 0.0d0
            !--------
        end do
        q6p(i) = 1.d0
        !--------
    end do
    !
    !
    do i = 1, 4
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q1p(g2i, g_i_, i) = q1 * g_g_q1p(g2i, g_i_, i - 1) + q1p(i - 1) * g_g_q1(g2i, g_i_)
            end do
            g_q1p(g_i_,i) = q1p(i - 1) * g_q1(g_i_) + q1 * g_q1p(g_i_,i - 1)
            !--------
        end do
        q1p(i) = q1p(i - 1) * q1
        !--------
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q2p(g2i, g_i_, i) = q2 * g_g_q2p(g2i, g_i_, i - 1) + q2p(i - 1) * g_g_q2(g2i, g_i_)
            end do
            g_q2p(g_i_,i) = q2p(i - 1) * g_q2(g_i_) + q2 * g_q2p(g_i_,i - 1)
            !--------
        end do
        q2p(i) = q2p(i - 1) * q2
        !--------
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q3p(g2i, g_i_, i) = q3 * g_g_q3p(g2i, g_i_, i - 1) + q3p(i - 1) * g_g_q3(g2i, g_i_)
            end do
            g_q3p(g_i_, i) = q3p(i - 1) * g_q3(g_i_) + q3 * g_q3p(g_i_, i - 1)
            !--------
        end do
        q3p(i) = q3p(i - 1) * q3
        !--------
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q4p(g2i, g_i_, i) = q4 * g_g_q4p(g2i, g_i_, i - 1) + q4p(i - 1) * g_g_q4(g2i,g_i_)
            end do
            g_q4p(g_i_, i) = q4p(i - 1) * g_q4(g_i_) + q4 * g_q4p(g_i_, i - 1)
            !--------
        end do
        q4p(i) = q4p(i - 1) * q4
        !--------
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q5p(g2i, g_i_, i) = q5 * g_g_q5p(g2i, g_i_, i - 1) + q5p(i - 1) * g_g_q5(g2i,g_i_)
            end do
            g_q5p(g_i_, i) = q5p(i - 1) * g_q5(g_i_) + q5 * g_q5p(g_i_, i - 1)
            !--------
        end do
        q5p(i) = q5p(i - 1) * q5
        !--------
        d2_b = dble(i)
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_d1_w(g2i, g_i_) = d2_b * g_g_q6(g2i, g_i_)
            end do
            g_d1_w(g_i_) = d2_b * g_q6(g_i_)
            !--------
        end do
        d1_w = dble(i) * q6
        d2_v = cos(d1_w)
        d1_p = -sin(d1_w)
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_q6p(g2i, g_i_, i) = d1_p * g_g_d1_w(g2i, g_i_)
            end do
            g_q6p(g_i_, i) = d1_p * g_d1_w(g_i_)
            !--------
        end do
        q6p(i) = d2_v
        !--------
    end do
    !
    do g_i_ = 1, g_p_
        do g2i = 1, g2p
            g_g_v(g2i, g_i_) = 0.0d0
        end do
        g_v(g_i_) = 0.0d0
        !--------
    end do
    v = zoe
    !--------
    !
    do i = 1, 152
        d3_v  = coeff(i) * q1p(ind1(i))
        d5_v  = d3_v * q2p(ind2(i))
        d7_v  = d5_v * q3p(ind3(i))
        d9_v  = d7_v * q4p(ind4(i))
        d11_v = d9_v * q5p(ind5(i))
        d6_b  = q6p(ind6(i)) * q5p(ind5(i))
        d7_b  = q6p(ind6(i)) * d9_v
        d8_b  = d6_b * q4p(ind4(i))
        d9_b  = d6_b * d7_v
        d10_b = d8_b * q3p(ind3(i))
        d11_b = d8_b * d5_v
        d13_b = d10_b * d3_v
        d14_b = d10_b * q2p(ind2(i)) * coeff(i)
        do g_i_ = 1, g_p_
            do g2i = 1, g2p
                g_g_v(g2i, g_i_) = g_g_v(g2i, g_i_) + d14_b * g_g_q1p( g2i, g_i_, ind1(i)) + d13_b * g_g_q2p(g2i, &
                    & g_i_, ind2(i)) + d11_b * g_g_q3p(g2i, g_i_, ind3(i)) + d9_b * g_g_q4p(g2i, g_i_, ind4(i)) + &
                    & d7_b * g_g_q5p(g2i, g_i_, ind5(i)) + d11_v * g_g_q6p( g2i, g_i_, ind6(i))
            end do
            g_v(g_i_) = d11_v * g_q6p(g_i_, ind6(i)) + d7_b * g_q5p(g_i_, ind5(i)) + d9_b * g_q4p(g_i_, ind4(i)) + &
                & d11_b * g_q3p(g_i_, ind3(i)) + d13_b * g_q2p(g_i_,ind2(i)) + d14_b * g_q1p(g_i_,ind1(i) ) + g_v(g_i_)
            !--------
        end do
        v = v + d11_v * q6p(ind6(i))
        !--------
    end do
    !
    !
    !
end subroutine g_g_vibpot

!
!
!

!> @brief Evaluate SETUP for the H2O2_GEN_1998 potential.
subroutine setup( n3tm )
    !
    implicit none
    integer :: n3tm
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
    integer, parameter :: n3tmmn = 12
    !
    !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
    !
    !WRITE (6, 1300)
    if (n3tm .lt. n3tmmn) then
        write (6, '(/,2X,T5,''Warning: N3TM is set equal to '',I3, '' but this potential routine'', /,2X,T14,''requir' // &
            'es N3TM be greater than or '', ''equal to '',I3,/)') n3tm, n3tmmn
        stop 'SETUP 1'
    end if
    !
    !
    !
end subroutine setup
