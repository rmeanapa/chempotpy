!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the ClH2_LEPS_TTGI2_1985 potential.
subroutine pes( x, igrad, path, p, g, d )

    implicit none

    ! number of electronic state
    integer, parameter :: nstates = 1
    integer, parameter :: natoms  = 3

    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    double precision, intent(in)    :: x(natoms,3)
    double precision, intent(out)   :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out)   :: d(nstates,nstates,natoms,3)

    double precision :: dq, q, v, xter
    integer          :: i, iatom, idebug, idir, idum, istate, j, nder, ndum
    common /pt1cm/ q(3),xter,dq(3)
    common /pt2cm/ nder, ndum, idebug, idum(19)
    double precision :: tx(9), r2(3), dedx(9), drdx(3,9)
    logical, save    :: first_time_data=.true.

    !initialize
    v = 0.d0
    g = 0.d0
    d = 0.d0

    do iatom=1, natoms
        do idir=1,3
            j     = (iatom-1)*3+idir
            tx(j) = x(iatom,idir)
        enddo
    enddo
    ! input cartesian is ClHH
    q(1) = sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2 +(x(1,3)-x(2,3))**2)/0.529177211
    q(2) = sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2 +(x(2,3)-x(3,3))**2)/0.529177211
    q(3) = sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2 +(x(1,3)-x(3,3))**2)/0.529177211

    nder = igrad
    dq   = 0.d0

    if(first_time_data) then
        call prepot(path)
        first_time_data = .false.
    endif

    call pot

    xter = xter*27.211386
    dq   = dq*51.422067

    r2 = q*0.529177211
    call evdrdx(tx, r2, drdx)
    dedx = 0.d0
    do i=1,9
        do j=1,3
            dedx(i) = dedx(i)+dq(j)*drdx(j,i)
        enddo
    enddo

    if (igrad==0) then
        do istate=1,nstates
            p(istate) = xter
        enddo
    else if (igrad==1) then
        do istate=1,nstates
            p(istate) = xter
        enddo
        do iatom=1,natoms
            do idir=1,3
                j               = (iatom-1)*3+idir
                g(1,iatom,idir) = dedx(j)
            enddo
        enddo
    else if (igrad==2) then
        write (*,*) 'Only energy and gradient are available'
    endif

end subroutine pes

!> @brief Evaluate EvdRdX for the ClH2_LEPS_TTGI2_1985 potential.
subroutine evdrdx( x, r, drdx )

    implicit none

    double precision, intent(in)  :: x(9), r(3)
    double precision, intent(out) :: drdx(3,9)

    integer :: i,j

    ! Initialize dRdX(3,9)
    do i=1,3
        do j=1,9
            drdx(i,j) = 0.0d0
        enddo
    enddo

    drdx(1,1) = (x(1)-x(4))/r(1)
    drdx(1,2) = (x(2)-x(5))/r(1)
    drdx(1,3) = (x(3)-x(6))/r(1)
    drdx(1,4) = -drdx(1,1)
    drdx(1,5) = -drdx(1,2)
    drdx(1,6) = -drdx(1,3)

    drdx(2,4) = (x(4)-x(7))/r(2)
    drdx(2,5) = (x(5)-x(8))/r(2)
    drdx(2,6) = (x(6)-x(9))/r(2)
    drdx(2,7) = -drdx(2,4)
    drdx(2,8) = -drdx(2,5)
    drdx(2,9) = -drdx(2,6)

    drdx(3,1) = (x(1)-x(7))/r(3)
    drdx(3,2) = (x(2)-x(8))/r(3)
    drdx(3,3) = (x(3)-x(9))/r(3)
    drdx(3,7) = -drdx(3,1)
    drdx(3,8) = -drdx(3,2)
    drdx(3,9) = -drdx(3,3)

end subroutine evdrdx

!
!> @brief Evaluate PREPOT for the ClH2_LEPS_TTGI2_1985 potential.
subroutine prepot( path )
    !
    !   System:          ClH2
    !   Functional form: Isaacson and Muckerman's diatomics-in-molecules ver
    !   Common name:     ClH2 DIMS
    !   Reference:       unpublished
    !   Cross Reference: S. C. Tucker, D. G. Truhlar,
    !                    B. C. Garrett, and A. D. Isaacson
    !                    J. Chem. Phys. 82, 4102 (1985)
    !
    !   Calling Sequence:
    !      PREPOT - initializes the potential's variables and
    !               must be called once before any calls to POT
    !      POT    - driver for the evaluation of the energy and the derivati
    !               of the energy with respect to the coordinates for a give
    !               geometric configuration
    !
    !   Units:
    !      energies    - hartrees
    !      coordinates - bohr
    !      derivatives - hartrees/bohr
    !
    !   Surfaces:
    !      ground electronic state
    !
    !   Zero of energy:
    !      The classical potential energy is set equal to zero for the Cl
    !      infinitely far from the H2 diatomic and R(H2) set equal to the
    !      H2 equilibrium diatomic value.
    !
    !   Parameters:
    !      Set in the BLOCK DATA subprogram PTPARM and read in from the file
    !      potclh2dims.dat that is linked to FORTRAN unit IPTRD1
    !
    !   Coordinates:
    !      Internal, Definition: R(1) = R(Cl-H)
    !                            R(2) = R(H-H)
    !                            R(3) = R(Cl-H)
    !
    !   Common Blocks (used between the calling program and this potential):
    !      /PT1CM/ R(3), ENERGY, DEDR(3)
    !        passes the coordinates, ground state electronic energy, and
    !        derivatives of the ground electronic state energy with respect
    !        to the coordinates.
    !      /PT2CM/ NDER, NDUM, IDEBUG, IDUM(19)
    !        passes the control flags where
    !        NDER   = 0 => no derivatives are computed
    !               = 1 => derivatives of the energy for the ground electron
    !                     state with respect to the coordinates are computed
    !        NDUM   - not used
    !        IDEBUG - control flag for priniting extra information
    !        IDEBUG = 0 => do not print extra debug information
    !        IDEBUG = 1 => print intermediate values to the file linked
    !                      to FORTRAN unit IPTPRT.
    !        IDUM   - not used
    !      /PT4CM/ IPTPRT, IPTRD1
    !        passes the FORTRAN unit number used for potential output
    !      /PT5CM/ EASYAB, EASYBC, EASYAC
    !        passes the energy in the three asymptotic valleys for an A + BC
    !        The energy in the AB valley, EASYAB, is equal to the energy of
    !        C atom "infinitely" far from the AB diatomic and R(AB) set equa
    !        Re(AB), the equilibrium bond length for the AB diatomic.
    !        In this potential the AB valley represents H infinitely far fro
    !        the ClH diatomic and R(ClH) equal to Re(ClH).  Similarly, the t
    !        EASYBC and EASYAC represent the energies in the H2 and the othe
    !        valleys, respectively.
    !
    !   Default Parameter Values:
    !      Variable      Default value
    !      NDER             1
    !      IDEBUG           0
    !      IPTPRT           6
    !      IPTRD1           4
    !
    !*****
    !
    implicit none

    character(len=1024), intent(in) :: path

    double precision :: arg, c, cs, ct, ct2, dd, ddgab, ddsxa, ddsxb, dduab, deldsx, delgab, dels11, &
        & dels12, dels22, delsab, delspx, deltab, deltpx, deltsx, deluab, dgab, dh, dh11, dh12, dh13, &
        & dh14, dh15, dh16, dh22, dh24, dh25, dh26, dh33, dh34, dh44, dh45, dh46, dh55, dh56, dh66, &
        & dnumb, dpux, dq, dsgh, dspxa, dspxb, dss11a, dss11b, dss12a, dss12b, dss22a, dss22b, dssab, &
        & dsxa, dsxb, dsxfit, dtdr1, dtdr2, dtdr3, dtpxa, dtpxb, dtsab, dtsxa, dtsxb, duab, easyab, &
        & easyac, easybc, en, eps, factr, forth, gabfit, gvec, half, one, q, r, s11fit, s12fit, &
        & s22fit, s2dsx, s2gab, s2s11, s2s12, s2s22, s2sab, s2spx, s2tab, s2tpx, s2tsx, s2uab, s3, &
        & s3dsx, s3gab, s3s11, s3s12, s3s22, s3sab, s3spx, s3tab, s3tpx, s3tsx, s3uab, sabfit, spxa, &
        & spxb, spxfit, sqrt3, ss11a, ss11b, ss12a, ss12b, ss22a, ss22b, ssab, ssgh, ssgx, st, st2, &
        & tabfit, test, th, tpxa, tpxb, tpxfit, tre, tsab, tsxa, tsxb, tsxfit, two, uabfit, vec, w, &
        & x, xter, xtra, zero
    integer :: i, idebug, idum, iflag, indx, iptprt, iptrd1, ix, j, k, nder, ndim, ndum, nroot, ns
    common /pt4cm/ iptprt, iptrd1
    common /pt2cm/ nder, ndum, idebug, idum(19)
    common /pt5cm/ easyab, easybc, easyac
    dimension th(6,6),vec(6,6),w(6),en(6),gvec(6),dh(6,6),r(3),dd(3)
    dimension dh11(3),dh12(3),dh13(3),dh14(3),dh15(3),dh16(3),dh22(3)
    dimension dh24(3),dh25(3),dh26(3),dh33(3),dh34(3),dh44(3),dh45(3)
    dimension dh46(3),dh55(3),dh56(3),dh66(3)
    dimension indx(3)
    character(len=80) :: title(2)
    common /pt1cm/ q(3),xter,dq(3)
    common /potcm2/ c(3)
    common/potcm3/zero,one,two,tre,half,forth,sqrt3, &
        & dpux,ssgx,dsgh,ssgh
    common/potcm4/x(121), &
        & sabfit(121),s2sab(121),s3sab(121),delsab(121), &
        & tabfit(121),s2tab(121),s3tab(121),deltab(121), &
        & gabfit(121),s2gab(121),s3gab(121),delgab(121), &
        & uabfit(121),s2uab(121),s3uab(121),deluab(121), &
        & s11fit(121),s2s11(121),s3s11(121),dels11(121), &
        & s12fit(121),s2s12(121),s3s12(121),dels12(121), &
        & s22fit(121),s2s22(121),s3s22(121),dels22(121), &
        & tsxfit(121),s2tsx(121),s3tsx(121),deltsx(121), &
        & spxfit(121),s2spx(121),s3spx(121),delspx(121), &
        & tpxfit(121),s2tpx(121),s3tpx(121),deltpx(121), &
        & dsxfit(121),s2dsx(121),s3dsx(121),deldsx(121),ns
    equivalence (th(1,1),dh(1,1))
    character(len=1024) :: file_path1

    data ndim/6/,nroot/1/
    data test/0.999999999d0/
    data indx /1, 2, 3/
    data title(1) /'    Final BNL DIM surface for Cl+H2 -- 5 and7 kca&
        &l/mole barriers'/
    data title(2) /'    Spline fits used in atomic units'/
    !
    !   Echo the name of the potential to the file linked to FORTRAN unit IP
    !WRITE (IPTPRT, 600)
    !WRITE (IPTPRT, 610) TITLE

    file_path1 = trim(path)//"/ClH2/potclh2dims.dat "
    !   Open the file with the potential input data
    open (unit=iptrd1, file=file_path1, form='FORMATTED', &
        & status='OLD', err=100)
    !     SET UP SPLINE FITS TO H2 CURVES
    call setspl
    !
    !   Close the input data file
    close (unit = iptrd1)
    !
    !   Initialize the energies in the three asymptotic valleys
    easyab = c(1)
    easybc = c(2)
    easyac = c(3)
    !
    return
    100 write (iptprt, '(/,2X,T5,''Error opening input data file'')')
    stop 'PREPOT 1'
    !
    entry pot
    !
    !***********************************************************************
    !
    !     ENTRY POT GIVES NROOT-TH POTENTIAL ENERGY FOR THE CLH2 SYSTEM
    !     ENTRY POT ALSO GIVES PARTIALS W.R.T. THREE DISTANCES
    !
    !     INTERNUCLEAR DISTANCES STORED IN Q(3) AND R(3)
    !     NROOT-TH STATE ENERGY STORED IN XTER (OR TERM FOR DIATOMIC)
    !     THREE PARTIALS STORED IN DD(3)
    !
    !***********************************************************************
    !
    !   Check the value of NDER
    if (nder .gt. 1) then
        write (iptprt, '(/,1X,T5,''Error: POT has been called with NDER = '', I5, /,1X,T12,''only the first derivativ' // &
            'es, NDER = 1, are '', ''coded in this potential'')') nder
        stop 'POT 1'
    endif
    !
    !     IAPP=1
    !     SET FLAG FOR EIGENVECTOR OUTPUT
    iflag = sign(1,nroot)
    nroot = abs(nroot)
    !
    !     FOLLOWING COORDINATE DEFINITIONS PROVIDE FOR PROPER SYSTEM
    do 15 ix=1,3
    15 r(ix)=q(indx(ix))
    !
    !     ZERO OUT UNUSED HAMILTONIAN ELEMENTS
    th(2,3) = zero
    th(3,5) = zero
    th(3,6) = zero
    !
    !     COMPUTE SINE AND COSINE OF ANGLE AXB, AND POWERS
    ct = (r(1)*r(1)+r(3)*r(3)-r(2)*r(2))/(two*r(1)*r(3))
    if(ct.gt.-test .and. ct.lt.test) go to 10
    ct = sign(test,ct)
    10 continue
    arg = one-ct*ct
    if(arg.lt.zero) arg=zero
    st  = sqrt(arg)
    ct2 = ct*ct
    st2 = st*st
    cs  = ct*st
    !     COMPUTE DIATOMIC ENERGIES FOR R(3) ARRAY
    call spline(ns,x,s11fit,r(1),ss11a,dss11a,s3,eps,s2s11,s3s11, &
        & dels11)
    call spline(ns,x,s11fit,r(3),ss11b,dss11b,s3,eps,s2s11,s3s11, &
        & dels11)
    call spline(ns,x,s12fit,r(1),ss12a,dss12a,s3,eps,s2s12,s3s12, &
        & dels12)
    call spline(ns,x,s12fit,r(3),ss12b,dss12b,s3,eps,s2s12,s3s12, &
        & dels12)
    call spline(ns,x,s22fit,r(1),ss22a,dss22a,s3,eps,s2s22,s3s22, &
        & dels22)
    call spline(ns,x,s22fit,r(3),ss22b,dss22b,s3,eps,s2s22,s3s22, &
        & dels22)
    call spline(ns,x,tsxfit,r(1),tsxa,dtsxa,s3,eps,s2tsx,s3tsx,deltsx)
    call spline(ns,x,tsxfit,r(3),tsxb,dtsxb,s3,eps,s2tsx,s3tsx,deltsx)
    call spline(ns,x,spxfit,r(1),spxa,dspxa,s3,eps,s2spx,s3spx,delspx)
    call spline(ns,x,spxfit,r(3),spxb,dspxb,s3,eps,s2spx,s3spx,delspx)
    call spline(ns,x,tpxfit,r(1),tpxa,dtpxa,s3,eps,s2tpx,s3tpx,deltpx)
    call spline(ns,x,tpxfit,r(3),tpxb,dtpxb,s3,eps,s2tpx,s3tpx,deltpx)
    call spline(ns,x,dsxfit,r(1),dsxa,ddsxa,s3,eps,s2dsx,s3dsx,deldsx)
    call spline(ns,x,dsxfit,r(3),dsxb,ddsxb,s3,eps,s2dsx,s3dsx,deldsx)
    call spline(ns,x,sabfit,r(2),ssab,dssab,s3,eps,s2sab,s3sab,delsab)
    call spline(ns,x,tabfit,r(2),tsab,dtsab,s3,eps,s2tab,s3tab,deltab)
    call spline(ns,x,gabfit,r(2),dgab,ddgab,s3,eps,s2gab,s3gab,delgab)
    call spline(ns,x,uabfit,r(2),duab,dduab,s3,eps,s2uab,s3uab,deluab)
    !
    !     CONSTRUCT DIM HAMILTONIAN IN ORTHOGONAL BASIS
    !
    th(1,1) = ss11a+forth*((ssab+tre*tsab)+ct2*(ss11b+tre*tsxb) +st2*(spxb+tre*tpxb))-dpux-two*dsgh
    th(1,2) = forth*sqrt3*(ssab-tsab-ct2*(ss11b-tsxb)-st2*(spxb-tpxb))
    th(1,3) = ss12a
    th(1,4) = -half*ct*ss12b
    th(1,5) = forth*cs*(ss11b+tre*tsxb-spxb-tre*tpxb)
    th(1,6) = -forth*sqrt3*cs*(ss11b-tsxb-spxb+tpxb)
    th(2,2) = tsxa+forth*((tre*ssab+tsab)+ct2*(tre*ss11b+tsxb) +st2*(tre*spxb+tpxb))-dpux-two*dsgh
    th(2,4) = half*sqrt3*ct*ss12b
    th(2,5) = th(1,6)
    th(2,6) = forth*cs*(tre*ss11b+tsxb-tre*spxb-tpxb)
    th(3,3) = ss22a+half*(dgab+duab)+dsxb-ssgx-ssgh-dsgh
    th(3,4) = half*(dgab-duab)
    th(4,4) = dsxa+half*(dgab+duab)+ss22b-ssgx-ssgh-dsgh
    th(4,5) = -half*st*ss12b
    th(4,6) = half*sqrt3*st*ss12b
    th(5,5) = spxa+forth*((ssab+tre*tsab)+st2*(ss11b+tre*tsxb) +ct2*(spxb+tre*tpxb))-dpux-two*dsgh
    th(5,6) = forth*sqrt3*(ssab-tsab-st2*(ss11b-tsxb)-ct2*(spxb-tpxb))
    th(6,6) = tpxa+forth*((tre*ssab+tsab)+st2*(tre*ss11b+tsxb) +ct2*(tre*spxb+tpxb))-dpux-two*dsgh
    !
    !     FILL OUT MATRIX AND DIAGONALIZE
    !
    do 30 i=1,ndim
        do 30 j=1,i
            th(i,j) = th(j,i)
    30 continue
    !
    call eign(ndim,th,vec,en,w,ndim)
    !     FOLLOWING STATEMENT ASSUMES ENERGY ZERO AT A + BC(R=RE)
    xter = en(nroot)+c(2)
    if(iflag.lt.zero) write(iptprt,'(/,2X,T5,''Vector:'',(T15,3E20.8))') (vec(i,nroot),i=1,ndim)
    !
    !************************* END OF FPOT3N *******************************
    !
    if (nder .ne. 1) return
    !     IF(IAPP.EQ.0) RETURN
    !
    !     SET  VECTOR FOR STATE OF INTEREST
    do 40 i=1,ndim
        gvec(i) = vec(i,nroot)
    40 continue
    !     COMPUTE PARTIALS OF U=COS(THETA)
    dtdr1 = (r(1)*r(1)+r(2)*r(2)-r(3)*r(3))/(two*r(1)*r(1)*r(3))
    dtdr2 = -r(2)/(r(1)*r(3))
    dtdr3 = (r(2)*r(2)+r(3)*r(3)-r(1)*r(1))/(two*r(3)*r(3)*r(1))
    factr = -one/st
    if(abs(st).lt.1.d-6) write(iptprt, &
        & '(/,2X,T5,''Warning: A serious error may result because '', ''sine(THETA) = '',' // '1PE20.10)') st
    !     COMPUTE PARTIALS OF HAMILTONIAN ELEMENTS
    xtra    = +half*ct*(ss11b+tre*tsxb-spxb-tre*tpxb)
    dh11(1) = dss11a+xtra*dtdr1
    dh11(2) = forth*(dssab+tre*dtsab)+xtra*dtdr2
    dh11(3) = forth*(ct2*(dss11b+tre*dtsxb)+st2*(dspxb+tre*dtpxb))+xtra* dtdr3
    xtra    = -half*sqrt3*ct*(ss11b-tsxb-spxb+tpxb)
    dh12(1) = xtra*dtdr1
    dh12(2) = forth*sqrt3*(dssab-dtsab)+xtra*dtdr2
    dh12(3) = -forth*sqrt3*(ct2*(dss11b-dtsxb)+st2*(dspxb-dtpxb))+xtra* dtdr3
    dh13(1) = dss12a
    dh13(2) = zero
    dh13(3) = zero
    xtra    = -half*ss12b
    dh14(1) = xtra*dtdr1
    dh14(2) = xtra*dtdr2
    dh14(3) = -half*ct*dss12b+xtra*dtdr3
    xtra    = forth*(ct2-st2)*(ss11b+tre*tsxb-spxb-tre*tpxb)*factr
    dh15(1) = xtra*dtdr1
    dh15(2) = xtra*dtdr2
    dh15(3) = forth*cs*(dss11b+tre*dtsxb-dspxb-tre*dtpxb)+xtra*dtdr3
    xtra    = -forth*sqrt3*(ct2-st2)*(ss11b-tsxb-spxb+tpxb)*factr
    dh16(1) = xtra*dtdr1
    dh16(2) = xtra*dtdr2
    dh16(3) = -forth*sqrt3*cs*(dss11b-dtsxb-dspxb+dtpxb)+xtra*dtdr3
    xtra    = +half*ct*(tre*ss11b+tsxb-tre*spxb-tpxb)
    dh22(1) = dtsxa+xtra*dtdr1
    dh22(2) = forth*(tre*dssab+dtsab)+xtra*dtdr2
    dh22(3) = forth*(ct2*(tre*dss11b+dtsxb)+st2*(tre*dspxb+dtpxb))+xtra* dtdr3
    xtra    = +half*sqrt3*ss12b
    dh24(1) = xtra*dtdr1
    dh24(2) = xtra*dtdr2
    dh24(3) = half*sqrt3*ct*dss12b+xtra*dtdr3
    dh25(1) = dh16(1)
    dh25(2) = dh16(2)
    dh25(3) = dh16(3)
    xtra    = forth*(ct2-st2)*(tre*ss11b+tsxb-tre*spxb-tpxb)*factr
    dh26(1) = xtra*dtdr1
    dh26(2) = xtra*dtdr2
    dh26(3) = forth*cs*(tre*dss11b+dtsxb-tre*dspxb-dtpxb)+xtra*dtdr3
    dh33(1) = dss22a
    dh33(2) = half*(ddgab+dduab)
    dh33(3) = ddsxb
    dh34(1) = zero
    dh34(2) = half*(ddgab-dduab)
    dh34(3) = zero
    dh44(1) = ddsxa
    dh44(2) = half*(ddgab+dduab)
    dh44(3) = dss22b
    xtra    = -half*ct*ss12b*factr
    dh45(1) = xtra*dtdr1
    dh45(2) = xtra*dtdr2
    dh45(3) = -half*st*dss12b+xtra*dtdr3
    xtra    = half*sqrt3*ct*ss12b*factr
    dh46(1) = xtra*dtdr1
    dh46(2) = xtra*dtdr2
    dh46(3) = half*sqrt3*st*dss12b+xtra*dtdr3
    xtra    = -half*ct*(ss11b+tre*tsxb-spxb-tre*tpxb)
    dh55(1) = dspxa+xtra*dtdr1
    dh55(2) = forth*(dssab+tre*dtsab)+xtra*dtdr2
    dh55(3) = forth*(st2*(dss11b+tre*dtsxb)+ct2*(dspxb+tre*dtpxb))+xtra* dtdr3
    xtra    = +half*sqrt3*ct*(ss11b-tsxb-spxb+tpxb)
    dh56(1) = xtra*dtdr1
    dh56(2) = forth*sqrt3*(dssab-dtsab)+xtra*dtdr2
    dh56(3) = -forth*sqrt3*(st2*(dss11b-dtsxb)+ct2*(dspxb-dtpxb))+xtra* dtdr3
    xtra    = -half*ct*(tre*ss11b+tsxb-tre*spxb-tpxb)
    dh66(1) = dtpxa+xtra*dtdr1
    dh66(2) = forth*(tre*dssab+dtsab)+xtra*dtdr2
    dh66(3) = forth*(st2*(tre*dss11b+dtsxb)+ct2*(tre*dspxb+dtpxb))+xtra* dtdr3
    !
    !     LOOP OVER R(I) DERIVATIVES
    !     FIRST ZERO OUT UNUSED ELEMENTS
    dh(2,3) = zero
    dh(3,5) = zero
    dh(3,6) = zero
    do 50 i=1,3
        !
        !     SET UP DERIVATIVE MATRIX
        !
        dh(1,1) = dh11(i)
        dh(1,2) = dh12(i)
        dh(1,3) = dh13(i)
        dh(1,4) = dh14(i)
        dh(1,5) = dh15(i)
        dh(1,6) = dh16(i)
        dh(2,2) = dh22(i)
        dh(2,4) = dh24(i)
        dh(2,5) = dh25(i)
        dh(2,6) = dh26(i)
        dh(3,3) = dh33(i)
        dh(3,4) = dh34(i)
        dh(4,4) = dh44(i)
        dh(4,5) = dh45(i)
        dh(4,6) = dh46(i)
        dh(5,5) = dh55(i)
        dh(5,6) = dh56(i)
        dh(6,6) = dh66(i)
        !
        !     FILL OUT DERIVATIVE MATRIX
        do 70 j=1,ndim
            do 70 k=1,j
                dh(j,k) = dh(k,j)
        70 continue
        !
        dnumb = zero
        do 80 j=1,ndim
            do 80 k=1,ndim
                dnumb = dnumb+gvec(j)*dh(j,k)*gvec(k)
        80 continue
        dd(i) = dnumb
        !     END OF DERIVATIVE LOOP
    50 continue
    !     PUT DERIVATIVES INTO CORRESPONDING SPOTS FOR CORRECT COORD CHOICE.
    do 90 ix=1,3
    90 dq(indx(ix))=dd(ix)
    !
    !************************ END OF DFPOT *********************************
    !
end subroutine prepot
!> @brief Evaluate SETSPL for the ClH2_LEPS_TTGI2_1985 potential.
subroutine setspl
    !     SETS UP SPLINE FITS FOR DIATOMIC CURVES
    implicit none

    double precision :: angta0, auterg, deldsx, delgab, dels11, dels12, dels22, delsab, delspx, &
        & deltab, deltpx, deltsx, deluab, dpux, dsgh, dsxfit, eps, fac, gabfit, half, one, quar, r, &
        & r1, r2, r3, rt3, s11fit, s12fit, s22fit, s2dsx, s2gab, s2s11, s2s12, s2s22, s2sab, s2spx, &
        & s2tab, s2tpx, s2tsx, s2uab, s3dsx, s3gab, s3s11, s3s12, s3s22, s3sab, s3spx, s3tab, s3tpx, &
        & s3tsx, s3uab, sabfit, spxfit, ssgh, ssgx, tabfit, tpxfit, tre, tsxfit, two, uabfit, x, &
        & zero
    integer :: i, iptprt, iptrd1, ns
    common /pt4cm/ iptprt, iptrd1
    common/potcm4/x(121), &
        & sabfit(121),s2sab(121),s3sab(121),delsab(121), &
        & tabfit(121),s2tab(121),s3tab(121),deltab(121), &
        & gabfit(121),s2gab(121),s3gab(121),delgab(121), &
        & uabfit(121),s2uab(121),s3uab(121),deluab(121), &
        & s11fit(121),s2s11(121),s3s11(121),dels11(121), &
        & s12fit(121),s2s12(121),s3s12(121),dels12(121), &
        & s22fit(121),s2s22(121),s3s22(121),dels22(121), &
        & tsxfit(121),s2tsx(121),s3tsx(121),deltsx(121), &
        & spxfit(121),s2spx(121),s3spx(121),delspx(121), &
        & tpxfit(121),s2tpx(121),s3tpx(121),deltpx(121), &
        & dsxfit(121),s2dsx(121),s3dsx(121),deldsx(121),ns
    common/potcm3/zero,one,two,tre,half,quar,rt3, &
        & dpux,ssgx,dsgh,ssgh
    data auterg/4.3598283d-11/,angta0/1.88972634d0/
    ns = 121
    !     SET UP X ARRAY IN BOHR
    fac = 0.1d0*angta0
    do 20 i=1,ns
    20 x(i)=dble(i-1)*fac
    !     READ IN FIT DATA
    read(iptrd1,'(2D20.13)') (sabfit(i),tabfit(i),i=1,ns)
    read(iptrd1,'(2D20.13)') (gabfit(i),uabfit(i),i=1,ns)
    read(iptrd1,'(4D20.13)') (s11fit(i),s12fit(i),s22fit(i),tsxfit(i),i=1,ns)
    read(iptrd1,'(3D20.13)') (spxfit(i),tpxfit(i),dsxfit(i),i=1,ns)
    !     CONVERT TO ATOMIC UNITS
    fac = one/auterg
    do 10 i=1,ns
        sabfit(i) = sabfit(i)*fac
        tabfit(i) = tabfit(i)*fac
        gabfit(i) = gabfit(i)*fac
        uabfit(i) = uabfit(i)*fac
        s11fit(i) = s11fit(i)*fac
        s12fit(i) = s12fit(i)*fac
        s22fit(i) = s22fit(i)*fac
        tsxfit(i) = tsxfit(i)*fac
        spxfit(i) = spxfit(i)*fac
        tpxfit(i) = tpxfit(i)*fac
        dsxfit(i) = dsxfit(i)*fac
    10 continue
    !     READ IN SECOND DERIVATIVES AT ENDS OF FIT CURVES
    !     AND CONVERT TO ATOMIC UNITS
    fac = fac*1.d-16/(angta0*angta0)
    read(iptrd1,'(2D20.13)') s2sab(1),s2tab(1),s2sab(ns),s2tab(ns)
    read(iptrd1,'(2D20.13)') s2gab(1),s2uab(1),s2gab(ns),s2uab(ns)
    read(iptrd1,'(4D20.13)') s2s11(1),s2s12(1),s2s22(1),s2tsx(1)
    read(iptrd1,'(4D20.13)') s2s11(ns),s2s12(ns),s2s22(ns),s2tsx(ns)
    read(iptrd1,'(3D20.13)') s2spx(1),s2tpx(1),s2dsx(1)
    read(iptrd1,'(3D20.13)') s2spx(ns),s2tpx(ns),s2dsx(ns)
    s2sab(1) = s2sab(1)*fac
    s2tab(1) = s2tab(1)*fac
    s2gab(1) = s2gab(1)*fac
    s2uab(1) = s2uab(1)*fac
    s2s11(1) = s2s11(1)*fac
    s2s12(1) = s2s12(1)*fac
    s2s22(1) = s2s22(1)*fac
    s2tsx(1) = s2tsx(1)*fac
    s2spx(1) = s2spx(1)*fac
    s2tpx(1) = s2tpx(1)*fac
    s2dsx(1) = s2dsx(1)*fac
    !     GENERATE SPLINE FITS
    call splnin(ns,x,sabfit,r,r1,r2,r3,eps,s2sab,s3sab,delsab)
    call splnin(ns,x,tabfit,r,r1,r2,r3,eps,s2tab,s3tab,deltab)
    call splnin(ns,x,gabfit,r,r1,r2,r3,eps,s2gab,s3gab,delgab)
    call splnin(ns,x,uabfit,r,r1,r2,r3,eps,s2uab,s3uab,deluab)
    call splnin(ns,x,s11fit,r,r1,r2,r3,eps,s2s11,s3s11,dels11)
    call splnin(ns,x,s12fit,r,r1,r2,r3,eps,s2s12,s3s12,dels12)
    call splnin(ns,x,s22fit,r,r1,r2,r3,eps,s2s22,s3s22,dels22)
    call splnin(ns,x,tsxfit,r,r1,r2,r3,eps,s2tsx,s3tsx,deltsx)
    call splnin(ns,x,spxfit,r,r1,r2,r3,eps,s2spx,s3spx,delspx)
    call splnin(ns,x,tpxfit,r,r1,r2,r3,eps,s2tpx,s3tpx,deltpx)
    call splnin(ns,x,dsxfit,r,r1,r2,r3,eps,s2dsx,s3dsx,deldsx)
end subroutine setspl
!> @brief Evaluate SPLINE for the ClH2_LEPS_TTGI2_1985 potential.
subroutine spline( n, x, y, t, ss, ss1, ss2, epsln, s2, s3, dely )
    !     JIM STINE'S SPLINE INTERPOLATOR
    implicit none

    double precision :: delsqs, dely, epsln, ht1, ht2, prod, s2, s3, ss, ss1, ss2, t, x, y
    integer          :: i, n
    dimension x(n),y(n),s2(n),s3(n),dely(n)
    i = 1
    if(t.gt.x(1)) go to 10
    ss  = y(1)
    ss1 = dely(1)
    return
    10 if(t.lt.x(n)) go to 20
    ss  = y(n)
    ss1 = dely(n)
    return
    20 continue
    56 if(t-x(i)) 60,17,57
    57 i=i+1
    go to 56
    59 i=n
    60 i=i-1
    17 ht1=t-x(i)
    ht2    = t-x(i+1)
    prod   = ht1*ht2
    ss2    = s2(i)+ht1*s3(i)
    delsqs = (s2(i)+s2(i+1)+ss2)/6.d-00
    ss     = y(i)+ht1*dely(i)+prod*delsqs
    ss1    = dely(i)+(ht1+ht2)*delsqs+prod*s3(i)/6.d-00
    61 return
end subroutine spline
!> @brief Evaluate SPLNIN for the ClH2_LEPS_TTGI2_1985 potential.
subroutine splnin( n, x, y, t, ss, ss1, ss2, epsln, s2, s3, dely )
    !     JIM STINE'S INITIALIZATION ROUTINE FOR SPLINE FITTING
    !     N=NO. POINTS (X,Y); T=ARGUMENT; SS=VALUE; SS1=1ST DERIV;
    !     SS2=2ND DERIV; S2,S3, AND DELY ARE COEFFICIENT ARRAYS
    !     SECOND DERIV'S AT ENDS STORED IN S2(1) AND S2(N) ON INPUT
    !     EPSLN IS OVERALL FIT CRITERION
    implicit none

    double precision :: b, c, delsqy, dely, epsln, eta, ex, h, h2, half, omega, s2, s3, ss, ss1, ss2, &
        & t, w, x, y
    integer :: i, idebug, idum, iptprt, iptrd1, knt, n, n1, nder, ndum
    common /pt4cm/ iptprt, iptrd1
    common /pt2cm/ nder, ndum, idebug, idum(19)
    dimension c(121),h(121),h2(121),b(121)
    dimension x(n),y(n),s2(n),s3(n),dely(n)
    data ex/1.1d0/,omega/1.0717968d0/,half/0.5d0/
    epsln = 6.4229d-13
    n1    = n-1
    knt   = 0
    do 51 i=1,n1
        h(i) = x(i+1)-x(i)
    51 dely(i)=(y(i+1)-y(i))/h(i)
    do 52 i=2,n1
        h2(i)  = h(i-1)+h(i)
        b(i)   = half*h(i-1)/h2(i)
        delsqy = (dely(i)-dely(i-1))/h2(i)
        s2(i)  = delsqy+delsqy
    52 c(i)=delsqy+s2(i)
    5 eta=0.0d0
    do 10 i=2,n1
        w = (c(i)-b(i)*s2(i-1)-(half-b(i))*s2(i+1)-s2(i))*omega
        if(abs(w)-eta) 10,10,9
        9 eta=abs(w)
    10 s2(i)=s2(i)+w
    knt = knt+1
    if (knt.gt.10) epsln=epsln*ex
    if(eta-epsln) 14,5,5
    14 do 53 i=1,n1
    53 s3(i)=(s2(i+1)-s2(i))/h(i)
    if (idebug .eq. 1) write(iptprt,'(2X,T5,''Final EPS = '', 1PE15.6)') epsln
end subroutine splnin
!> @brief Evaluate EIGN for the ClH2_LEPS_TTGI2_1985 potential.
subroutine eign( nn, a, vec, eig, w, nd )
    !     MATRIX DIAGONALIZATION ROUTINE FOR REAL SYMMETRIC CASE
    !     HOUSEHOLDER METHOD
    !     RHO=UPPERLIMIT FOR OFF-DIAGONAL ELEMENT
    !     NN=SIZE OF MATRIX
    !     A=MATRIX (ONLY LOWER TRIANGLE IS USED ,THIS IS DESTROYED)
    !     EIG=RETURNED EIGENVALUES IN ALGEBRAIC ASCENDING ORDER
    !     VEC=RETURNED EIGENVECTORS IN COLUMNS
    implicit none

    double precision :: a, a2, b, beta, betasq, cosa, cosap, d, dia, dif, eig, g, gamma, p, pp, ppbr, &
        & ppbs, q, qj, r1, r12, r2, rho, rhosq, s, sgn, shift, sina, sina2, sqrts, sum, temp, u, vec, &
        & w, wj, wtaw
    integer :: i, i1, iord, iposv, itemp, ivpos, j, k, l, lv, m, n, n1, n2, nd, nn, np, npas, nr, &
        & nrr, nt, nv
    dimension eig(nd),w(nd)
    dimension gamma(6),beta(6),betasq(6)
    dimension p(6),q(6)
    equivalence (p(1),beta(1)),(q(1),beta(1))
    dimension iposv(6),ivpos(6),iord(6)
    dimension a(nd,nd),vec(nd,nd)
    rho   = 1.0d-14
    rhosq = rho*rho
    n     = nn
    if(n.eq.0) go to 560
    1 n1=n-1
    n2       = n-2
    gamma(1) = a(1,1)
    if(n2) 280,270,120
    120 do 260 nr=1,n2
        b = a(nr+1,nr)
        s = 0.0d0
        do 130 i=nr,n2
            !     PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION
        130 s=s+a(i+2,nr)**2
        a(nr+1,nr) = 0.0d0
        if(s) 250,250,140
        140 s=s+b*b
        sgn = +1.0d0
        if(b) 150,160,160
        150 sgn=-1.0d0
        160 sqrts=sqrt(s)
        d          = sgn/(sqrts+sqrts)
        temp       = sqrt(.5d0+b*d)
        w(nr)      = temp
        a(nr+1,nr) = temp
        d          = d/temp
        b          = -sgn*sqrts
        !     D IS FACTOR OF PROPORTIONALITY  NOW COMPUTE AND SAVE W VECTOR
        !     EXTRA SINGLY SUBSCRIPTED W VECTOR USED FOR SPEED
        do 170 i=nr,n2
            temp   = d*a(i+2,nr)
            w(i+1) = temp
        170 a(i+2,nr)=temp
        !     PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR
        !     SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K)
        wtaw = 0.0d0
        do 220 i=nr,n1
            sum = 0.0d0
            do 180 j=nr,i
            180 sum=sum+a(i+1,j+1)*w(j)
            i1 = i+1
            if(n1-i1) 210,190,190
            190 do 200 j=i1,n1
            200 sum=sum+a(j+1,i+1)*w(j)
            210 p(i)=sum
        220 wtaw=wtaw+sum*w(i)
        !     P VECTOR AND SCALAR K NOW STORED, NEXT COMPUTE Q VECTOR.
        do 230 i=nr,n1
            !     NOW FORM PAP MATRIX, REQUIRED PART
        230 q(i)=p(i)-wtaw*w(i)
        do 240 j=nr,n1
            qj = q(j)
            wj = w(j)
            do 240 i=j,n1
        240 a(i+1,j+1)=a(i+1,j+1)-2.0d0*(w(i)*qj+wj*q(i))
        250 beta(nr)=b
        betasq(nr) = b*b
    260 gamma(nr+1)=a(nr+1,nr+1)
    270 b=a(n,n-1)
    beta(n-1)   = b
    betasq(n-1) = b*b
    gamma(n)    = a(n,n)
    280 betasq(n)=0.0d0
    !     ADJOIN AN IDENTITY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS
    do 290 j=1,n
        do 290 i=1,n
    290 vec(i,j)=0.0d0
    do 300 i=1,n
    300 vec(i,i)=1.0d0
    m    = n
    sum  = 0.0d0
    npas = 1
    go to 400
    310 sum=sum+shift
    cosa = 1.0d0
    g    = gamma(1)-shift
    pp   = g
    ppbs = pp*pp+betasq(1)
    ppbr = sqrt(ppbs)
    do 370 j=1,m
        cosap = cosa
        if(ppbs.ne.0.0d0) go to 320
        311 sina=0.0d0
        sina2 = 0.0d0
        cosa  = 1.0d0
        go to 350
        320 sina=beta(j)/ppbr
        sina2 = betasq(j)/ppbs
        cosa  = pp/ppbr
        !     POSTMULTIPLY IDAENTITY BY P-TRANSPOSE MATRIX
        nt = j+npas
        if(nt.lt.n) go to 330
        321 nt=n
        330 do 340 i=1,nt
            temp       = cosa*vec(i,j)+sina*vec(i,j+1)
            vec(i,j+1) = -sina*vec(i,j)+cosa*vec(i,j+1)
        340 vec(i,j)=temp
        350 dia=gamma(j+1)-shift
        u        = sina2*(g+dia)
        gamma(j) = g+u
        g        = dia-u
        pp       = dia*cosa-sina*cosap*beta(j)
        if(j.ne.m) go to 360
        351 beta(j)=sina*pp
        betasq(j) = sina2*pp*pp
        go to 380
        360 ppbs=pp*pp+betasq(j+1)
        ppbr    = sqrt(ppbs)
        beta(j) = sina*ppbr
    370 betasq(j)=sina2*ppbs
    380 gamma(m+1)=g
    !     TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT
    npas = npas+1
    if(betasq(m).gt.rhosq) go to 410
    390 eig(m+1)=gamma(m+1)+sum
    400 beta(m)=0.0d0
    betasq(m) = 0.0d0
    m         = m-1
    if(m.eq.0) go to 430
    401 if(betasq(m).le.rhosq) go to 390
    !     TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN VALUE
    !     AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT
    410 a2=gamma(m+1)
    r2   = .5d0*a2
    r1   = .5d0*gamma(m)
    r12  = r1+r2
    dif  = r1-r2
    temp = sqrt(dif*dif+betasq(m))
    r1   = r12+temp
    r2   = r12-temp
    dif  = abs(a2-r1)-abs(a2-r2)
    if(dif.lt.0.0d0) go to 420
    411 shift=r2
    go to 310
    420 shift=r1
    go to 310
    430 eig(1)=gamma(1)+sum
    !     INITIALIZE AUXILARY TABLES REQUIRED FOR REARRANGING THE VECTORS
    do 440 j=1,n
        iposv(j) = j
        ivpos(j) = j
    440 iord(j)=j
    !     USE A TRANSPOSITON SORT TO ORDER THE EIGENVALUES
    m = n
    go to 470
    450 do 460 j=1,m
        if(eig(j).le.eig(j+1)) go to 460
        451 temp=eig(j)
        eig(j)    = eig(j+1)
        eig(j+1)  = temp
        itemp     = iord(j)
        iord(j)   = iord(j+1)
        iord(j+1) = itemp
    460 continue
    470 m=m-1
    if(m.ne.0) go to 450
    471 if(n1.eq.0) go to 500
    472 do 490 l=1,n1
        nv = iord(l)
        np = iposv(nv)
        if(np.eq.l) go to 490
        473 lv=ivpos(l)
        ivpos(np) = lv
        iposv(lv) = np
        do 480 i=1,n
            temp     = vec(i,l)
            vec(i,l) = vec(i,np)
        480 vec(i,np)=temp
    490 continue
    !     BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX
    500 do 550 nrr=1,n
        k = n1
        510 k=k-1
        if(k.le.0) go to 550
        511 sum=0.0d0
        do 520 i=k,n1
        520 sum=sum+vec(i+1,nrr)*a(i+1,k)
        sum = sum+sum
        do 530 i=k,n1
        530 vec(i+1,nrr)=vec(i+1,nrr)-sum*a(i+1,k)
        go to 510
    550 continue
    560 return
    !  END OF EIGN
end subroutine eign
block data ptparm
    implicit none

    double precision :: c, dpux, dsgh, half, one, quar, rt3, ssgh, ssgx, tre, two, zero
    integer          :: idebug, idum, iptprt, iptrd1, nder, ndum
    common /pt4cm/iptprt, iptrd1
    common /pt2cm/ nder, ndum, idebug, idum(19)
    common /potcm2/ c(3)
    common/potcm3/zero,one,two,tre,half,quar,rt3, &
        & dpux,ssgx,dsgh,ssgh
    !   Initialize the control flags for the potenetial
    data iptprt, iptrd1 /6, 4/
    data ndum, nder /0, 1/
    data idebug /0/
    data idum /19*0/
    !   Initialize the parameters for the potential
    data zero,one,two,tre/0.d0,1.d0,2.d0,3.d0/
    data half,quar/0.5d0,0.25d0/
    data rt3/1.73205080756888d0/
    data dpux,ssgx,dsgh,ssgh/0.d0,-0.132774d0,0.d0,0.5d0/
    data c/0.169552d0,0.1744734d0,0.169552d0/
end

!> @brief Evaluate POTINFO for the ClH2_LEPS_TTGI2_1985 potential.
subroutine potinfo

    implicit none

    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, deesdr, degsdr, deijdr, descart, descartnu, &
        & dgscart, dgscartnu, dijcart, dijcartnu, engyes, engygs, engyij, ezero, pengyes, pengygs, pengyij, r
    integer :: i, icartr, indexes, irctnt, ireorder, j, k, kediag, keoffd, kmax, ksdiag, ksoffd, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    character(len=75)           :: ref(5)
    integer, parameter          :: n3atom = 75
    integer, parameter          :: natom  = 25
    integer, parameter          :: isurf  = 5
    integer, parameter          :: jsurf  = isurf*(isurf+1)/2
    double precision, parameter :: pi     = 3.141592653589793d0
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
    common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte, &
        & cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    write(nflag(18),'(/)')
    do i =1,5
        write(nflag(18),'(2x,a75)') ref(i)
    end do
    write(nflag(18),'(/)')
    kmax = 0
    do i = 1,isurf+1
        do j = 1,isurf+1
            if(nasurf(i,j).ne.0.and.kmax.lt.max(i,j)) kmax = max(i,j)
        enddo
    enddo
    write(nflag(18),'(2x,'' MAX. AND ACTUAL NO. OF EXCITED SURFACES: '',I3,5x,I3)') msurf,kmax-1
    if(kmax-1.gt.msurf) then
        write(6,*) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
        stop
    endif
    ksdiag = 0
    kediag = 0
    do i = 2,isurf+1
        if(nasurf(i,i).ne.0) then
            kediag = i-1
            if(ksdiag.eq.0) ksdiag = i-1
        endif
    enddo
    ksoffd = 0
    keoffd = 0
    k      = 0
    do i = 1,isurf
        do j = i+1,isurf+1
            k = k+1
            if(nasurf(i,j)+nasurf(j,i).ne.0) then
                keoffd = k
                if(ksoffd.eq.0) ksoffd = k
            endif
        enddo
    enddo
    write(nflag(18),'(2x,'' MAX. AND ACTUAL ORDER OF DERIVATIVES:    '',I3,5x,I3)') mder,nder
    if(nder.gt.mder) then
        write(6,*) ' WRONG INPUT ON ORDER OF DERIVATIVES'
        stop
    endif
    if(nflag(19).eq.1) then
        write(nflag(18),'(/)')
        write(nflag(18),'(2x,''Cartesian coordinates are supplied by'',/, 2x,''the user in the array CART.'',//)')
        write(nflag(18),'(2x,''Provide cartesian coordinates in the'',/, 2x,''following order using the array CART'',' // &
            '//, 2x,'' CART(1,1)...CART(1,3)   => ATOM 1'',/, 2x,'' CART(2,1)...CART(2,3)   => AT' // &
            'OM 2'',/, 2x,'' CART(3,1)...CART(3,3)   => ATOM 3'',/, 2x,'' CART(N,1)...CART(N,3)  ' // &
            ' => ATOM N'',/, 2x,''CART(25,1)...CART(25,3)  => ATOM 25'',/)')
        write(nflag(18),'(2x,''If the user wishes to relabel the atoms,'',/, 2x,''set the variable IREORDER equal to ' // &
            '1'',/, 2x,''in the PARAMETER statement.  The user'',/, 2x,''must also specify the ne' // &
            'w labeling'',/, 2x,''scheme.  This is done using the array'',/, 2x,''NULBL in the fo' // &
            'llowing manner:'',//, 2x,''NULBL(i) = j'',/, 2x,''where:  i => old index'',/, 2x,'' ' // &
            '       j => new index'',//)')
        write(nflag(18),'(2x,''Cartesian coordinates can be provided to'',/, 2x,''the potential routine in a variety ' // &
            'of units.'',/, 2x,''The input units will be converted to Bohr'',/, 2x,''based on the' // &
            ' following values of the NFLAG'',/, 2x,''variable:'',//, 2x,''NFLAG(1)  =  1  =>  CA' // &
            'RTESIANS IN BOHR (no'',/, 2x,''                    conversion required)'',/, 2x,''NF' // &
            'LAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS'',//)')
        write(nflag(18),'(2x,''The value of the energy and derivatives'',/, 2x,''(if computed) can be reported in a v' // &
            'ariety'',/, 2x,''units.  A units conversion will take place'',/, 2x,''as specified b' // &
            'y the following values of the'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(2) = 1 =>  ' // &
            'ENERGIES REPORTED IN HARTEEE'',/, 2x,''NFLAG(2) = 2 =>  ENERGIES REPORTED IN mHARTRE' // &
            'E'',/, 2x,''NFLAG(2) = 3 =>  ENERGIES REPORTED IN eV'',/, 2x,''NFLAG(2) = 4 =>  ENER' // &
            'GIES REPORTED IN kcal/mol'',/, 2x,''NFLAG(2) = 5 =>  ENERGIES REPORTED IN cm**-1'',/' // &
            '/)')
        write(nflag(18),'(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
            'he'',/, 2x,''NFLAG variable:'',//, 2x,''NFLAG(1)=1 & NFLAG(2)=1 => DERIVATIVES REPOR' // &
            'TED IN'',/, 2x,''                           HARTEEE/BOHR'',/, 2x,''NFLAG(1)=1 & NFLA' // &
            'G(2)=2 => DERIVATIVES REPORTED IN'',/, 2x,''                           mHARTREE/BOHR' // &
            ''',/, 2x,''NFLAG(1)=1 & NFLAG(2)=3 => DERIVATIVES REPORTED IN'',/, 2x,''            ' // &
            '               eV/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=4 => DERIVATIVES REPORTED IN' // &
            ''',/, 2x,''                           kcal/mol/BOHR'',/, 2x,''NFLAG(1)=1 & NFLAG(2)=' // &
            '5 => DERIVATIVES REPORTED IN'',/, 2x,''                           cm**-1/BOHR'',//)')
        write(nflag(18),'(2x,''A units conversion will take place'',/, 2x,''as specified by the following values of t' // &
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

!> @brief Evaluate ANCVRT for the ClH2_LEPS_TTGI2_1985 potential.
subroutine ancvrt
    implicit none

    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, deesdr, degsdr, deijdr, descart, descartnu, &
        & dgscart, dgscartnu, dijcart, dijcartnu, engyes, engygs, engyij, ezero, pengyes, pengygs, pengyij, r
    integer :: i, ianum, icartr, inc1, inc2, ind, ind2, indexes, inu, irctnt, ireorder, isave, isum, j, jsave, kediag, &
        & keoffd, ksdiag, ksoffd, mder, msurf, nasurf, natoms, nder, nflag, nulbl
    character(len=75)  :: ref(5)
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
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
    common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte, &
        & cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    dimension ianum(7,32)
    dimension isave(natom),jsave(natom)
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
    do i=1,7
        do j=1,32
            ianum(i,j)      = 0
            periodic_1(i,j) = ' '
        end do
    end do
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
    do i=1,natoms
        isave(i) = 0
        jsave(i) = 0
        name1(i) = '  '
        name2(i) = '  '
    end do
    iblank = ' '
    do ind=1,natoms
        do i=1,7
            do j=1,32
                if(indexes(ind).eq.ianum(i,j)) then
                    isave(ind) = i
                    jsave(ind) = j
                end if
            end do
        end do
    end do

    do ind=1,natoms
        ind2 = nulbl(ind)
        if(ind2.eq.0) ind2=ind
    end do
    inc1 = 0
    do ind=1,irctnt-1
        inc1        = inc1+1
        name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
    end do
    inc2 = 0
    do ind=irctnt,natoms
        inc2        = inc2+1
        name2(inc2) = periodic_1(isave(ind),jsave(ind))(:2)
    end do
    if(nflag(1).eq.2) distance = 'ANGSTROMS           '
    if(nflag(2).eq.2) then
        units = 'MILLIHARTREE        '
    elseif(nflag(2).eq.3) then
        units = 'EV                  '
    elseif(nflag(2).eq.4) then
        units = 'KCAL PER MOLE       '
    elseif(nflag(2).eq.5) then
        units = 'WAVENUMBERS         '
    elseif(nflag(2).eq.6) then
        units = 'KILOJOULES PER MOLE '
    endif
    cnvrtd  = 1.d0
    cnvrte  = 1.d0
    cnvrtde = 1.d0
    if(nflag(1).eq.2) cnvrtd = bohr_a
    if(nflag(2).eq.2) then
        cnvrte = cnvrte*htomillh
    elseif(nflag(2).eq.3) then
        cnvrte = cnvrte*htoev
    elseif(nflag(2).eq.4) then
        cnvrte = cnvrte*htokcal
    elseif(nflag(2).eq.5) then
        cnvrte = cnvrte*htowave
    elseif(nflag(2).eq.6) then
        cnvrte = cnvrte*htokj
    endif
    cnvrtde = cnvrte/cnvrtd
    isum    = 0
    do inu=1,25
        isum = isum + nulbl(inu)
    end do
    ireorder = 0
    if(isum.ne.0) ireorder = 1
end subroutine ancvrt
!> @brief Evaluate CARTOU for the ClH2_LEPS_TTGI2_1985 potential.
subroutine cartou
    use chempotpy_legacy_coordinates, only: convert_cartesian_input
    implicit none

    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, descartnu, dgscartnu, &
        & dijcartnu
    integer :: i, icartr, indexes, irctnt, ireorder, j, kediag, keoffd, ksdiag, ksoffd, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte, &
        & cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    call convert_cartesian_input( cart, nulbl, natoms, cnvrtd, ireorder == 1, cartnu )
end subroutine cartou

!> @brief Evaluate CARTTOR for the ClH2_LEPS_TTGI2_1985 potential.
subroutine carttor
    use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
    implicit none
    integer :: coordinate_status

    double precision :: anuzero, cart, cartnu, cmm, cphi, degsdr, engygs, flm, hhd, hym, phi, pi, q, &
        & q1, q2, q3, r, theta1, theta2, xcm1, xcm2, xcm3, xrm1, xrm2, ycm1, ycm2, ycm3, yrm1, yrm2, &
        & zcm1, zcm2, zcm3, zrm1, zrm2
    integer :: i, icartr, ind, indexes, irctnt, k, l, mder, msurf, nasurf, &
        & natoms, nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    call cartesian_to_potlib_coordinates( cartnu, natoms, icartr, r, coordinate_status )
    if (coordinate_status /= 0) then
        write(nflag(18),'(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') icartr
        stop
    end if
end subroutine carttor

!> @brief Evaluate EUNITZERO for the ClH2_LEPS_TTGI2_1985 potential.
subroutine eunitzero
    implicit none

    double precision :: anuzero, cart, cnvrtd, cnvrtde, cnvrte, deesdr, degsdr, deijdr, descart, &
        & descartnu, dgscart, dgscartnu, dijcart, dijcartnu, engyes, engygs, engyij, ezero, pengyes, &
        & pengygs, pengyij, r
    integer :: i, ireorder, j, kediag, keoffd, ksdiag, ksoffd, nasurf, &
        & nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2

    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
    common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
    common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte, &
        & cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    pengygs = engygs * cnvrte - anuzero
    if(ksdiag.ne.0) then
        do i=ksdiag,kediag
            pengyes(i) = engyes(i) * cnvrte - anuzero
        end do
    endif
    if(ksoffd.ne.0) then
        do j=ksoffd,keoffd
            pengyij(j) = engyij(j) * cnvrte
        end do
    endif
end subroutine eunitzero

!> @brief Evaluate RTOCART for the ClH2_LEPS_TTGI2_1985 potential.
subroutine rtocart
    implicit none

    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, deesdr, degsdr, deijdr, &
        & descartnu, dgscartnu, dijcartnu, engyes, engygs, engyij, eps, epsp, ezero, phi, r, r1, r2, &
        & rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, termy, termz, theta1, &
        & theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes, ygs, yij, za, zb, zc, &
        & zd
    integer :: i, icartr, ind, indexes, irctnt, ireorder, j, j1, j2, k, kediag, keoffd, ksdiag, ksoffd, m1, mder, msurf, &
        & nasurf, natoms, nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    common /pt3cm/  ezero(isurf+1)
    common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
    common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
    common /utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte,cnvrtde, &
        & ireorder,ksdiag,kediag,ksoffd,keoffd
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    dimension ygs(n3atom),yes(n3atom,isurf),yij(n3atom,jsurf)
    if(icartr.eq.1) then
        do i = 1, natoms
            ind            = 3*i-2
            dgscartnu(i,1) = degsdr(ind)
            dgscartnu(i,2) = degsdr(ind+1)
            dgscartnu(i,3) = degsdr(ind+2)
            if(ksdiag.ne.0) then
                do j = ksdiag,kediag
                    descartnu(i,1,j) = deesdr(ind,j)
                    descartnu(i,2,j) = deesdr(ind+1,j)
                    descartnu(i,3,j) = deesdr(ind+2,j)
                end do
            endif
            if(keoffd.ne.0) then
                do k = ksoffd,keoffd
                    dijcartnu(i,1,k) = deijdr(ind,k)
                    dijcartnu(i,2,k) = deijdr(ind+1,k)
                    dijcartnu(i,3,k) = deijdr(ind+2,k)
                end do
            endif
        end do
    elseif(icartr.eq.2) then
        do i = 1, natoms
            dgscartnu(i,1) = 0.d0
            dgscartnu(i,2) = 0.d0
            dgscartnu(i,3) = 0.d0
            if(ksdiag.ne.0) then
                do j1=ksdiag,kediag
                    descartnu(i,1,j1) = 0.d0
                    descartnu(i,2,j1) = 0.d0
                    descartnu(i,3,j1) = 0.d0
                enddo
            endif
            if(ksoffd.ne.0) then
                do j2=ksoffd,keoffd
                    dijcartnu(i,1,j2) = 0.d0
                    dijcartnu(i,2,j2) = 0.d0
                    dijcartnu(i,3,j2) = 0.d0
                enddo
            endif
            do j = 1,natoms
                if(j.lt.i) then
                    m1 = natoms*(j-1) - (j*(j-1))/2 + i-j
                elseif(j.gt.i) then
                    m1 = natoms*(i-1) - (i*(i-1))/2 + j-i
                else
                    go to 20
                endif
                y              = degsdr(m1)
                termx          = (cartnu(i,1)-cartnu(j,1))/r(m1)
                termy          = (cartnu(i,2)-cartnu(j,2))/r(m1)
                termz          = (cartnu(i,3)-cartnu(j,3))/r(m1)
                dgscartnu(i,1) = dgscartnu(i,1) + termx*y
                dgscartnu(i,2) = dgscartnu(i,2) + termy*y
                dgscartnu(i,3) = dgscartnu(i,3) + termz*y
                if(ksdiag.gt.0) then
                    y = deesdr(m1,j1)
                    do j1=ksdiag,kediag
                        descartnu(i,1,j1) = descartnu(i,1,j1) + termx*y
                        descartnu(i,2,j1) = descartnu(i,2,j1) + termy*y
                        descartnu(i,3,j1) = descartnu(i,3,j1) + termz*y
                    enddo
                elseif(ksoffd.gt.0) then
                    do j2=ksoffd,keoffd
                        y                 = deijdr(m1,j2)
                        dijcartnu(i,1,j2) = dijcartnu(i,1,j2) + termx*y
                        dijcartnu(i,2,j2) = dijcartnu(i,2,j2) + termy*y
                        dijcartnu(i,3,j2) = dijcartnu(i,3,j2) + termz*y
                    enddo
                endif
                20 continue
            enddo
        enddo
    elseif(icartr.eq.3) then
        do i = 1, natoms
            ygs(i) = degsdr(i)/r(i)
            if(ksdiag.ne.0) then
                do j=ksdiag,kediag
                    yes(i,j) = deesdr(i,j)/r(i)
                enddo
            endif
            if(ksoffd.ne.0) then
                do k=ksoffd,keoffd
                    yij(i,k) = deijdr(i,k)/r(i)
                enddo
            endif
        enddo
        do k = 1,3
            term12         = cartnu(1,k)-cartnu(2,k)
            term23         = cartnu(2,k)-cartnu(3,k)
            term13         = cartnu(1,k)-cartnu(3,k)
            dgscartnu(1,k) = term12*ygs(1) + term13*ygs(3)
            dgscartnu(2,k) = -term12*ygs(1) + term23*ygs(2)
            dgscartnu(3,k) = -term13*ygs(3) - term23*ygs(2)
            if(ksdiag.ne.0) then
                do j1=ksdiag,kediag
                    descartnu(1,k,j1) = term12*yes(1,j1) + term13*yes(3,j1)
                    descartnu(2,k,j1) = -term12*yes(1,j1) + term23*yes(2,j1)
                    descartnu(3,k,j1) = -term13*yes(3,j1) - term23*yes(2,j1)
                enddo
            endif
            if(ksoffd.ne.0) then
                do j2=ksoffd,keoffd
                    dijcartnu(1,k,j2) = term12*yij(1,j2) + term13*yij(3,j2)
                    dijcartnu(2,k,j2) = -term12*yij(1,j2) + term23*yij(2,j2)
                    dijcartnu(3,k,j2) = -term13*yij(3,j2) - term23*yij(2,j2)
                enddo
            endif
        enddo
    elseif(icartr.eq.4) then
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
        write(nflag(18),'(2X,'' WRONG ICARTR FOR DERIVATIVE; ICARTR ='',I5//)') icartr
        stop
    endif
end subroutine rtocart

!> @brief Evaluate DEDCOU for the ClH2_LEPS_TTGI2_1985 potential.
subroutine dedcou
    implicit none

    double precision :: anuzero, cart, cartnu, cnvrtd, cnvrtde, cnvrte, descart, descartnu, dgscart, &
        & dgscartnu, dijcart, dijcartnu, pengyes, pengygs, pengyij
    integer :: i, icartr, indexes, irctnt, ireorder, j, k, kediag, keoffd, ksdiag, ksoffd, mder, msurf, nasurf, natoms, &
        & nder, nflag, nulbl
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom = 75
    integer, parameter :: natom  = 25
    integer, parameter :: isurf  = 5
    integer, parameter :: jsurf  = isurf*(isurf+1)/2
    common /utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf), &
        & dijcartnu(natom,3,jsurf),cnvrtd,cnvrte,cnvrtde, &
        & ireorder,ksdiag,kediag,ksoffd,keoffd
    common/usrocm/ pengygs,pengyes(isurf), &
        & pengyij(jsurf), &
        & dgscart(natom,3),descart(natom,3,isurf), &
        & dijcart(natom,3,jsurf)
    common/infocm/ cartnu(natom,3),indexes(natom), &
        & irctnt,natoms,icartr,mder,msurf,ref
    common/usricm/ cart(natom,3),anuzero, &
        & nulbl(natom),nflag(20), &
        & nasurf(isurf+1,isurf+1),nder
    if (ireorder.eq.1) then
        do i = 1, natoms
            dgscart(i,1) = dgscartnu(nulbl(i),1) * cnvrtde
            dgscart(i,2) = dgscartnu(nulbl(i),2) * cnvrtde
            dgscart(i,3) = dgscartnu(nulbl(i),3) * cnvrtde
            if(ksdiag.ne.0) then
                do j=ksdiag,kediag
                    descart(i,1,j) = descartnu(nulbl(i),1,j) * cnvrtde
                    descart(i,2,j) = descartnu(nulbl(i),2,j) * cnvrtde
                    descart(i,3,j) = descartnu(nulbl(i),3,j) * cnvrtde
                end do
            endif
            if(ksoffd.ne.0) then
                do k=ksoffd,keoffd
                    dijcart(i,1,k) = dijcartnu(nulbl(i),1,k) * cnvrtde
                    dijcart(i,2,k) = dijcartnu(nulbl(i),2,k) * cnvrtde
                    dijcart(i,3,k) = dijcartnu(nulbl(i),3,k) * cnvrtde
                end do
            endif
        end do
    else
        do i = 1, natoms
            dgscart(i,1) = dgscartnu(i,1) * cnvrtde
            dgscart(i,2) = dgscartnu(i,2) * cnvrtde
            dgscart(i,3) = dgscartnu(i,3) * cnvrtde
            if(ksdiag.ne.0) then
                do j=ksdiag,kediag
                    descart(i,1,j) = descartnu(i,1,j) * cnvrtde
                    descart(i,2,j) = descartnu(i,2,j) * cnvrtde
                    descart(i,3,j) = descartnu(i,3,j) * cnvrtde
                end do
            endif
            if(ksoffd.ne.0) then
                do k=ksoffd,keoffd
                    dijcart(i,1,k) = dijcartnu(i,1,k) * cnvrtde
                    dijcart(i,2,k) = dijcartnu(i,2,k) * cnvrtde
                    dijcart(i,3,k) = dijcartnu(i,3,k) * cnvrtde
                end do
            endif
        end do
    endif
end subroutine dedcou
