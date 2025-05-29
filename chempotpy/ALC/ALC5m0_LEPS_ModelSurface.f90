module potential_energy_surface
implicit none
private
public :: pes
! constants and common variables go here
integer, parameter :: nstates = 1
integer, parameter :: natoms  = 3
real(8) :: rr(3), e, dedrr(3)
integer :: nsurf, nder, ndum(21)

! convert
!energies from hartree to ckcalmol
function hartree2kcalmol( e, eout )
function kcalmol2hartree( e, eout )
function ang2bohr
function borh2ang




common /pt1cm/ rr(3), e, dedrr(3)
common /pt2cm/ nsurf, nder, ndum(21)
common /pt3cm/  ezero(isurf+1)
common /pt4cm/ iptprt, idum(19)
common /pt5cm/ easyab,easybc, easyac
common /satocm/ d(3), re(3), beta(3), z(3) 
common /vacom/ an1, fns
parameter (n3atom=75)
parameter (natom=25)
parameter (isurf = 5)
parameter (jsurf = isurf*(isurf+1)/2)
common /utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf),dijcartnu(natom,3,jsurf),cnvrtd,cnvrte,cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
common/usrocm/ pengygs,pengyes(isurf),pengyij(jsurf),dgscart(natom,3),descart(natom,3,isurf),dijcart(natom,3,jsurf)
common/infocm/ cartnu(natom,3),indexes(natom),irctnt,natoms,icartr,mder,msurf,ref
common/usricm/ cart(natom,3),anuzero,nulbl(natom),nflag(20),nasurf(isurf+1,isurf+1),nder

contains
      ! core library
      ! xyz --> x
      ! xyz --> r
      ! workflow
      ! 1 check nstates and natoms
      ! 2 
    subroutine check_states_atoms( e, number_of_states, number_of_atoms )
        integer, intent(in) :: number_of_states, number_of_atoms 
        real(8), intent(in) :: e(:,:)
        if( size(x,dim=1) /= number_of_states ) THROW_HARD()
    end subroutine check_states_atoms


    subroutine pes( xyz, gradient_order, energy, gradient, hessian )
        implicit none
        integer, intent(in)  :: gradient_order 
        real(8), intent(in)  :: xyz(natoms,3)
        real(8), intent(out) :: energy(nstates), gradient(nstates,natoms,3)
        real(8), intent(out) :: hessian(nstates,nstates,natoms,3)
        real(8)       :: tx(9), r2(3), dedx(9), drdx(3,9)
        logical, save :: first_time_data=.true.
        !initialize 
        v        = 0.d0
        gradient = 0.d0
        hessian  = 0.d0
        ! interatomic distances
        do i_atom = 1, natoms
            do idir = 1,3
                j     = (i_atom-1)*3+idir
                tx(j) = x(iatom,idir)
            enddo
        enddo
        rr(1) = sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2+(x(1,3)-x(2,3))**2)/0.529177211
        rr(2) = sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2+(x(2,3)-x(3,3))**2)/0.529177211
        rr(3) = sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2+(x(1,3)-x(3,3))**2)/0.529177211
        ! input cartesian is clhh
        nsurf = 0
        nder  = igrad
        dedrr = 0.d0     
        if( first_time_data )then
            call prepot
            first_time_data = .false.
        endif
        call pot
        e     = e*27.211386
        dedrr = dedrr*51.422067
        r2    = rr*0.529177211
        call evdrdx(tx, r2, drdx)
        dedx  = 0.d0
        do i = 1, 9
            do j = 1, 3
                dedx(i) = dedx(i)+dedrr(j)*drdx(j,i)
            enddo     
        enddo
        if( igrad==0 )then
            do istate=1,nstates
                p(istate) = e
            enddo
        elseif( igrad==1 )then
            do istate=1,nstates
                p(istate)=e
            enddo
            do iatom=1,natoms
                do idir=1,3
                    j               = (iatom-1)*3+idir
                    g(1,iatom,idir) = dedx(j)
                enddo
            enddo
        elseif( igrad==2 )then
            write (*,*) 'only energy and gradient are available'
        endif
      end subroutine pes

  subroutine evdrdx( x, r, drdx )
      implicit none
      !real(8), intent(in)  :: x(9), r(3)
      !real(8), intent(out) :: drdx(3,9)
      real(8), intent(in)  :: x(:)
      real(8), intent(in)  :: r(:)
      real(8), intent(out) :: drdx(:)
      integer :: i,j
      ! initialize drdx(3,9)
      drdx      = 0.0d0
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

!   system:          alc
!   functional form: eleps + va
!                    Extended leps plus an auxiliary barrier term
!   common name:     chc5m0
!   reference:       unpublished
!   cross reference: M. M. Kreevoy, D. Ostovic, 
!                    D. G. Truhlar, and B. C. Garrett
!                    J. Phys. Chem. 90, 3766-3774 (1986)
!   note:            in a-l-c, a = acceptor group, l = h (d), 
!                    and c = methyl donor
!   calling sequence: 
!      prepot - initializes the potential's variables and
!               must be called once before any calls to pot
!      pot    - driver for the evaluation of the energy and the derivatives 
!               of the energy with respect to the coordinates for a given 
!               geometric configuration
!   units: 
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   surfaces: 
!      ground electronic state
!   zero of energy: 
!      the classical potential energy is set equal to zero for the al
!      infinitely far from the lc diatomic and r(lc) set equal to the
!      lc equilibrium diatomic value.
!   parameters:
!      set in the block data subprogram ptparm
!   coordinates:
!      internal, definition: r(1) = r(al-l)
!                            r(2) = r(l-c)
!                            r(3) = r(al-c)
!   common blocks (used between the calling program and this potential):
!      /pt1cm/ r(3), energy, dedr(3)
!        passes the coordinates, ground state electronic energy, and
!        derivatives of the ground electronic state energy with respect
!        to the coordinates.
!      /pt2cm/ nsurf, nder, ndum(21)
!        passes the control flags where
!        nder  = 0 => no derivatives are computed
!              = 1 => derivatives of the energy for the ground electronic
!                     state with respect to the coordinates are computed
!        ndum  - not used
!      /pt4cm/ iptprt, idum(19)
!        iptprt passes the fortran unit number used for potential output
!        idum   not used
!      /pt5cm/ easyab, easybc, easyac
!        passes the energy in the three asymptotic valleys for an a + bc system.
!        the energy in the ab valley, easyab, is equal to the energy of the
!        c atom "infinitely" far from the ab diatomic and r(ab) set equal to
!        re(ab), the equilibrium bond length for the ab diatomic.
!        similarly, the terms easybc and easyac represent the energies
!        in the bc and the ac valleys, respectively.
!   default parameter values:
!      variable      default value
!      nsurf            0
!      nder             1
!      iptprt           6
      subroutine prepot
      implicit none
      parameter (root2=1.41421356d0)
      parameter (ckcal = 627.5095d0)
      parameter (cangs =   0.529177106d0)
      dimension zpo(3),op3z(3),zp3(3),tzp3(3),top3z(3),do4z(3),b(3)
      dimension x(3),coul(3),exch(3)
      save
      ! echo the potential parameters to the file linked to fortran unit iptprt
      do i = 1,3
      ! convert to atomic units
      d(i)     = d(i)/ckcal 
      re(i)    = re(i)/cangs
      beta(i)  = beta(i)*cangs 
      ! compute useful constants
      zpo(i)   = 1.0d0 + z(i)
      op3z(i)  = 1.0d0 + 3.0d0*z(i)
      top3z(i) = 2.0d0*op3z(i)
      zp3(i)   = z(i) + 3.0d0
      tzp3(i)  = 2.0d0*zp3(i)
      do4z(i)  = d(i)/4.0d0/zpo(i)
      b(i)     = beta(i)*do4z(i)*2.0d0
      enddo
      ! initialize the energies in the asymptotic valleys
      easyab   = d(1)
      easybc   = d(2)
      easyac   = d(3)
      return
      entry pot
      ! check the value of nder
      ! compute the eleps part of the potential energy
      do i = 1,3
          x(i)    = exp(-beta(i)*(rr(i)-re(i)))
          coul(i) = do4z(i)*(zp3(i)*x(i)-top3z(i))*x(i)
          exch(i) = do4z(i)*(op3z(i)*x(i)-tzp3(i))*x(i)
      enddo
      rad = sqrt((exch(1)-exch(2))**2+(exch(2)-exch(3))**2+(exch(3)-exch(1))**2)
      e = coul(1) + coul(2) + coul(3) - (rad/root2) + d(2)
      ! compute the derivative of the eleps part of the potential, if nder=1
      if( nder .ne. 1) go to 700
      s = exch(1) + exch(2) + exch(3)
      do 22 i = 1,3
      dedrr(i) = 0.d0     
      if(x(i).lt.1.d-30) go to 22 
      t= (3.0d0*exch(i)-s)/root2*(op3z(i)*x(i)-zp3(i))
!   if rad is very small compared to t print a warning to the file linked 
!   to fortran iptprt and set t/rad equal to t.
      if(abs(rad).lt.1.d-32.and.abs(t).gt.1.d-12) then
        write(iptprt,6000) t,rad     
      else if(abs(rad).gt.1.d-32) then
        t = t/rad 
      end if
      dedrr(i) = b(i)*x(i)*(t-zp3(i)*x(i)+op3z(i)) 
22    continue
700   continue
      ! compute the auxiliary barrier term, va, where 
      ! va = a{[(r1-r2)*(r2-r3)*(r3-r1)]^2}*exp[-alpha*{(r1+r2+r3)^3}]
      r    = rr(1)+rr(2)+rr(3)
      r2   = r*r
      r3   = r2*r
      exns = exp(-fns*r3)
      ens  = 0.d0
      wnt  = (rr(1)-rr(2))*(rr(2)-rr(3))*(rr(3)-rr(1))
      wn   = abs(wnt)
      wn2  = wn*wn
      ens  = (an1*wn2)*exns
      ! add the energy terms together
      e    = e+ens
      ! compute the derivative of the va term with respect to the 
      ! coordinates, if nder = 1
      if( nder == 1 )then
          delta    =-1.d0
          if( wn .eq. wnt ) delta =1.d0
          ensp1    = 0.d0
          ensp2    = 0.d0
          ensp3    = 0.d0
          enspwn   = (an1*wn*2.d0)*exns
          enspr    = (-3.d0)*fns*r2*ens
          ! dens/dxi=(dens/dwn)(dwn/dxi)+(dens/dr)(dr/dxi)
          x21      = rr(1)*rr(1)
          x22      = rr(2)*rr(2)
          x23      = rr(3)*rr(3)
          wnp1     = (2.d0*rr(1)*(rr(3)-rr(2))+x22-x23)*delta
          wnp2     = (2.d0*rr(2)*(rr(1)-rr(3))+x23-x21)*delta
          wnp3     = (2.d0*rr(3)*(rr(2)-rr(1))+x21-x22)*delta
          ensp1    = enspwn*wnp1+enspr
          ensp2    = enspwn*wnp2+enspr
          ensp3    = enspwn*wnp3+enspr
          dedrr(1) = dedrr(1)+ensp1
          dedrr(2) = dedrr(2)+ensp2
          dedrr(3) = dedrr(3)+ensp3
      endif
600   format(/,1x,'*****','potential energy surface',1x,'*****', &
            //,1x,t5,'al + lc chc5m0 potential energy surface', &
            //,1x,t5,'parameters:', &
              /,2x,t5,'bond', t46, 'ai-l', t58, 'l-c', t69, 'c-ai', &
              /,2x,t5,'dissociation energies (kcal/mol):',  &
              t44, f10.5, t55, f10.5, t66, f10.5, &
              /,2x,t5,'equilibrium bond lengths (angstroms):',  &
              t44, f10.5, t55, f10.5, t66, f10.5, &
              /,2x,t5,'morse beta parameters (angstroms**-1):',  &
              t44, f10.5, t55, f10.5, t66, f10.5, &
              /,2x,t5,'sato parameters:',  &
              t44, f10.5, t55, f10.5, t66, f10.5)
610   format(/,2x,t5,'parameters for the auxiliary barrier term ', &
             '(a.u.)',/,2x,t20,'a = ',f15.8,t45,'alpha = ',t53,f15.8,
            //,1x,'*****')
900   format(/,1x,t5,'error: pot has been called with nder = ', i5, &
             /,1x,t12,'only the first derivatives, nder = 1, are ', &
                      'coded in this potential')
6000  format(/,2x,'in the potential routine: t, rad = ',1p,2e15.7, &
             /,2x,'t/rad has been set equal to t')  
      end subrotine 

block data ptparm
implicit double precision (a-h,o-z)
common /pt4cm/ iptprt, idum(19)
common /pt2cm/ nsurf, nder, ndum(21)
common /satocm/ de(3), re(3), beta(3), z(3) 
common /vacom/ an1, fns
! initialize the flags and the i/o unit numbers for the potential
data iptprt, nder, nsurf /6, 1, 0/
data ndum /21*0/ 
! initialize the eleps potential parameters; the energy parameters are in
! kcal/mol, and the lengths are in angstroms.
data de/ 73.00000d0, 73.00000d0, 80.00000d0/
data re/ 1.10000d0, 1.10000d0, 1.54000d0/
data beta/ 2.17500d0, 2.17500d0, 2.07800d0/
data z/ 0.00000d0, 0.00000d0, -0.18900d0/
! initialize the potential parameters for the auxiliary barrier term; 
! the parameters are in hartree atomic units.
data an1, fns / 0.01d0, 0.004d0/
end

subroutine potinfo
      implicit none
      character*75 ref(5)
      parameter (n3atom=75)
      parameter (natom=25)
      parameter (isurf = 5)
      parameter (jsurf = isurf*(isurf+1)/2)
      parameter (pi = 3.141592653589793d0)
      common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
      common /pt3cm/  ezero(isurf+1)
      common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
      common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
      common/infocm/ cartnu(natom,3),indexes(natom),
     +               irctnt,natoms,icartr,mder,msurf,ref
      common/usrocm/ pengygs,pengyes(isurf),
     +               pengyij(jsurf),
     +               dgscart(natom,3),descart(natom,3,isurf),
     +               dijcart(natom,3,jsurf)
      common/usricm/ cart(natom,3),anuzero,
     +               nulbl(natom),nflag(20),
     +               nasurf(isurf+1,isurf+1),nder
      common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf),
     +               dijcartnu(natom,3,jsurf),cnvrtd,cnvrte,
     +               cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
      write(nflag(18),'(/)')
      do i =1,5
         write(nflag(18),97) ref(i)
         97 format(2x,a75)
      enddo
      write(nflag(18),'(/)')
      kmax = 0
      do i = 1,isurf+1
         do j = 1,isurf+1
            if(nasurf(i,j).ne.0.and.kmax.lt.max(i,j)) kmax = max(i,j)
         enddo
      enddo
      write(nflag(18),'(2x,' max. and actual no. of excited surfaces: ',i3,5x,i3)') msurf,kmax-1
      if(kmax-1.gt.msurf) then
         write(6,*) ' wrong input on number of excited surfaces'
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
      k = 0
      do i = 1,isurf
         do j = i+1,isurf+1
            k = k+1
            if(nasurf(i,j)+nasurf(j,i).ne.0) then
               keoffd = k
               if(ksoffd.eq.0) ksoffd = k
            endif
         enddo
      enddo
      write(nflag(18),103) mder,nder
103   format(2x,' max. and actual order of derivatives:    ',i3,5x,i3)
      if(nder.gt.mder) then
         write(6,*) ' wrong input on order of derivatives'
         stop
      endif
      if(nflag(19).eq.1) then
         write(nflag(18),'(/)')
         write(nflag(18),120)
 120     format(2x,'cartesian coordinates are supplied by',/,               &
                2x,'the user in the array cart.',//)
         write(nflag(18),125)
 125     format(2x,'provide cartesian coordinates in the',/,               &
                2x,'following order using the array cart',//,               &
                2x,' cart(1,1)...cart(1,3)   => atom 1',/,               &
                2x,' cart(2,1)...cart(2,3)   => atom 2',/,               &
                2x,' cart(3,1)...cart(3,3)   => atom 3',/,               &
                2x,' cart(n,1)...cart(n,3)   => atom n',/,               &
                2x,'cart(25,1)...cart(25,3)  => atom 25',/)
         write(nflag(18),130)
 130     format(2x,'if the user wishes to relabel the atoms,',/,               &
                2x,'set the variable ireorder equal to 1',/,               &
                2x,'in the parameter statement.  the user',/,               &
                2x,'must also specify the new labeling',/,               &
                2x,'scheme.  this is done using the array',/,               &
                2x,'nulbl in the following manner:',//,               &
                2x,'nulbl(i) = j',/,               &
                2x,'where:  i => old index',/,               &
                2x,'        j => new index',//)
         write(nflag(18),150)
 150     format(2x,'cartesian coordinates can be provided to',/,               &
                2x,'the potential routine in a variety of units.',/,               &
                2x,'the input units will be converted to bohr',/,               &
                2x,'based on the following values of the nflag',/,               &
                2x,'variable:',//,               &
                2x,'nflag(1)  =  1  =>  cartesians in bohr (no',/,               &
                2x,'                    conversion required)',/,               &
                2x,'nflag(1)  =  2  =>  cartesians in angstroms',//)
         write(nflag(18),160)
 160     format(2x,'the value of the energy and derivatives',/,               &
                2x,'(if computed) can be reported in a variety',/,               &
                2x,'units.  a units conversion will take place',/,               &
                2x,'as specified by the following values of the',/,               &
                2x,'nflag variable:',//,               &
                2x,'nflag(2) = 1 =>  energies reported in harteee',/,               &
                2x,'nflag(2) = 2 =>  energies reported in mhartree',/,               &
                2x,'nflag(2) = 3 =>  energies reported in ev',/,               &
                2x,'nflag(2) = 4 =>  energies reported in kcal/mol',/,               &
                2x,'nflag(2) = 5 =>  energies reported in cm**-1',//)               &
         write(nflag(18),165)
 165     format(2x,'a units conversion will take place',/,               &
             2x,'as specified by the following values of the',/,               &
             2x,'nflag variable:',//,               &
             2x,'nflag(1)=1 & nflag(2)=1 => derivatives reported in',/,               &
             2x,'                           harteee/bohr',/,               &
             2x,'nflag(1)=1 & nflag(2)=2 => derivatives reported in',/,               &
             2x,'                           mhartree/bohr',/,               &
             2x,'nflag(1)=1 & nflag(2)=3 => derivatives reported in',/,               &
             2x,'                           ev/bohr',/,               &
             2x,'nflag(1)=1 & nflag(2)=4 => derivatives reported in',/,               &
             2x,'                           kcal/mol/bohr',/,               &
             2x,'nflag(1)=1 & nflag(2)=5 => derivatives reported in',/,               &
             2x,'                           cm**-1/bohr',//)
         write(nflag(18),170)
 170     format(2x,'a units conversion will take place',/,
             2x,'as specified by the following values of the',/,
             2x,'nflag variable:',//,
             2x,'nflag(1)=2 & nflag(2)=1 => derivatives reported in',/,               &
             2x,'                           harteee/angstrom',/,               &
             2x,'nflag(1)=2 & nflag(2)=2 => derivatives reported in',/,               &
             2x,'                           mhartree/angstrom',/,               &
             2x,'nflag(1)=2 & nflag(2)=3 => derivatives reported in',/,               &
             2x,'                           ev/angstrom',/,               &
             2x,'nflag(1)=2 & nflag(2)=4 => derivatives reported in',/,               &
             2x,'                           kcal/mol/angstrom',/,               &
             2x,'nflag(1)=2 & nflag(2)=5 => derivatives reported in',/,               &
             2x,'                           cm**-1/angstrom',//)
      endif
      return
      end


subroutine ancvrt
    implicit none
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom=75, natom=25, isurf=5
    integer, parameter :: jsurf=isurf*(isurf+1)/2
    integer            :: ind, i, k
    character(len=20)  :: distance, units
    character(len=1)   :: iblank 
    character(len=2)   :: name1(natom), name2(natom)
    character(len=3)   :: periodic_1(7,32)
    !   common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    !   common /pt3cm/  ezero(isurf+1)
    !   common /pt4cm/  engyes(isurf),deesdr(n3atom,isurf)
    !   common /pt5cm/  engyij(jsurf),deijdr(n3atom,jsurf)
    !   common/usrocm/ pengygs,pengyes(isurf),
    !  +               pengyij(jsurf),
    !  +               dgscart(natom,3),descart(natom,3,isurf),
    !  +               dijcart(natom,3,jsurf)
    !   common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf),
    !  +               dijcartnu(natom,3,jsurf),cnvrtd,cnvrte,
    !  +               cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    !   common/infocm/ cartnu(natom,3),indexes(natom),
    !  +               irctnt,natoms,icartr,mder,msurf,ref
    !   common/usricm/ cart(natom,3),anuzero,
    !  +               nulbl(natom),nflag(20),
    !  +               nasurf(isurf+1,isurf+1),nder
      dimension ianum(7,32)
      dimension isave(natom),jsave(natom)
      parameter(        pi = 3.141592653589793d0)
      parameter(    clight = 2.99792458d08)
      parameter(     cmu_0 = 4.0d0*pi*1.0d-07)
      parameter(cepsilon_0 = 1.0d0/(cmu_0*clight**2))
      parameter(        ce = 1.602176462d-19)
      parameter(   cplanck = 6.62606876d-34)
      parameter(      cm_e = 9.10938188d-31)
      parameter(      cang = 1.0d-10)
      parameter( cavogadro = 6.02214199d23)
      parameter(     ckcal = 4.184d10)
      parameter(  htomillh = 1000.d0)
      parameter(     htoev = 27.2113834d0)
      parameter(   htokcal = 627.509470d0)
      parameter(   htowave = 219474.631d0)
      parameter(     htokj = 2625.49962d0)
      parameter(    bohr_a = .5291772083d0)
      do i=1,7
          do j=1,32
             ianum(i,j)      = 0
             periodic_1(i,j) = ' '
          enddo 
      enddo 
      distance = 'bohr                '
      units    = 'hartree             '
      ianum(1,1)  =  1
      ianum(1,32) =  2
      ianum(2,1)  =  3
      ianum(2,2)  =  4
      ianum(2,27) =  5
      ianum(2,28) =  6
      ianum(2,29) =  7
      ianum(2,30) =  8
      ianum(2,31) =  9
      ianum(2,32) = 10
      ianum(3,1)  = 11
      ianum(3,2)  = 12
      ianum(3,27) = 13
      ianum(3,28) = 14
      ianum(3,29) = 15
      ianum(3,30) = 16
      ianum(3,31) = 17
      ianum(3,32) = 18
      ianum(4,1)  = 19
      ianum(4,2)  = 20
      ianum(4,17) = 21
      ianum(4,18) = 22
      ianum(4,19) = 23
      ianum(4,20) = 24
      ianum(4,21) = 25
      ianum(4,22) = 26
      ianum(4,23) = 27
      ianum(4,24) = 28
      ianum(4,25) = 29
      ianum(4,26) = 30
      ianum(4,27) = 31
      ianum(4,28) = 32
      ianum(4,29) = 33
      ianum(4,30) = 34
      ianum(4,31) = 35
      ianum(4,32) = 36
      ianum(5,1)  = 37
      ianum(5,2)  = 38
      ianum(5,17) = 39
      ianum(5,18) = 40
      ianum(5,19) = 41
      ianum(5,20) = 42
      ianum(5,21) = 43
      ianum(5,22) = 44
      ianum(5,23) = 45
      ianum(5,24) = 46
      ianum(5,25) = 47
      ianum(5,26) = 48
      ianum(5,27) = 49
      ianum(5,28) = 50
      ianum(5,29) = 51
      ianum(5,30) = 52
      ianum(5,31) = 53
      ianum(5,32) = 54
      ianum(6,1)  = 55
      ianum(6,2)  = 56
      ianum(6,3)  = 57
      ianum(6,4)  = 58
      ianum(6,5)  = 59
      ianum(6,6)  = 60
      ianum(6,7)  = 61
      ianum(6,8)  = 62
      ianum(6,9)  = 63
      ianum(6,10) = 64
      ianum(6,11) = 65
      ianum(6,12) = 66
      ianum(6,13) = 67
      ianum(6,14) = 68
      ianum(6,15) = 69
      ianum(6,16) = 70
      ianum(6,17) = 71
      ianum(6,18) = 72
      ianum(6,19) = 73
      ianum(6,20) = 74
      ianum(6,21) = 75
      ianum(6,22) = 76
      ianum(6,23) = 77
      ianum(6,24) = 78
      ianum(6,25) = 79
      ianum(6,26) = 80
      ianum(6,27) = 81
      ianum(6,28) = 82
      ianum(6,29) = 83
      ianum(6,30) = 84
      ianum(6,31) = 85
      ianum(6,32) = 86
      ianum(7,1)  = 87
      ianum(7,2)  = 88
      ianum(7,3)  = 89
      ianum(7,4)  = 90
      ianum(7,5)  = 91
      ianum(7,6)  = 92
      ianum(7,7)  = 93
      ianum(7,8)  = 94
      ianum(7,9)  = 95
      ianum(7,10) = 96
      ianum(7,11) = 97
      ianum(7,12) = 98
      ianum(7,13) = 99
      ianum(7,14) = 100
      ianum(7,15) = 101
      ianum(7,16) = 102
      ianum(7,17) = 103
      ianum(7,18) = 104
      ianum(7,19) = 105
      ianum(7,20) = 106
      ianum(7,21) = 107
      ianum(7,22) = 108
      ianum(7,23) = 109
      ianum(7,24) = 110
      ianum(7,25) = 111
      ianum(7,26) = 112
      ianum(7,27) = 113
      ianum(7,28) = 114
      ianum(7,29) = 115
      ianum(7,30) = 116
      ianum(7,31) = 117
      ianum(7,32) = 120
      periodic_1(1,1)   = 'H  '
      periodic_1(1,32)  = 'He '
      periodic_1(2,1)   = 'Li '
      periodic_1(2,2)   = 'Be '
      periodic_1(2,27)  = 'B  '
      periodic_1(2,28)  = 'C  '
      periodic_1(2,29)  = 'N  '
      periodic_1(2,30)  = 'O  '
      periodic_1(2,31)  = 'F  '
      periodic_1(2,32)  = 'Ne '
      periodic_1(3,1)   = 'Na '
      periodic_1(3,2)   = 'Mg '
      periodic_1(3,27)  = 'Al '
      periodic_1(3,28)  = 'Si '
      periodic_1(3,29)  = 'P  '
      periodic_1(3,30)  = 'S  '
      periodic_1(3,31)  = 'Cl '
      periodic_1(3,32)  = 'Ar '
      periodic_1(4,1)   = 'K  '
      periodic_1(4,2)   = 'Ca '
      periodic_1(4,17)  = 'Sc '
      periodic_1(4,18)  = 'Ti '
      periodic_1(4,19)  = 'V  '
      periodic_1(4,20)  = 'Cr '
      periodic_1(4,21)  = 'Mn '
      periodic_1(4,22)  = 'Fe '
      periodic_1(4,23)  = 'Co '
      periodic_1(4,24)  = 'Ni '
      periodic_1(4,25)  = 'Cu '
      periodic_1(4,26)  = 'Zn '
      periodic_1(4,27)  = 'Ga '
      periodic_1(4,28)  = 'Ge '
      periodic_1(4,29)  = 'As '
      periodic_1(4,30)  = 'Se '
      periodic_1(4,31)  = 'Br '
      periodic_1(4,32)  = 'Kr '
      periodic_1(5,1)   = 'Rb '
      periodic_1(5,2)   = 'Sr '
      periodic_1(5,17)  = 'Y  '
      periodic_1(5,18)  = 'Zr '
      periodic_1(5,19)  = 'Nb '
      periodic_1(5,20)  = 'Mo '
      periodic_1(5,21)  = 'Tc '
      periodic_1(5,22)  = 'Ru '
      periodic_1(5,23)  = 'Rh '
      periodic_1(5,24)  = 'Pd '
      periodic_1(5,25)  = 'Ag '
      periodic_1(5,26)  = 'Cd '
      periodic_1(5,27)  = 'In '
      periodic_1(5,28)  = 'Sn '
      periodic_1(5,29)  = 'Sb '
      periodic_1(5,30)  = 'Te '
      periodic_1(5,31)  = 'I  '
      periodic_1(5,32)  = 'Xe '
      do i=1,natoms
          isave(i) = 0
          jsave(i) = 0
          name1(i) = '  '
          name2(i) = '  '
      enddo
      iblank=' '
      do ind=1,natoms
          do i=1,7
              do j=1,32
                  if( indexes(ind) .eq. ianum(i,j) )then
                      isave(ind) = i
                      jsave(ind) = j
                  endif
              enddo
          enddo
      enddo
      do ind = 1, natoms
          ind2 = nulbl(ind)
          if( ind2 == 0 ) ind2 = ind
      enddo
      inc1 = 0
      do ind = 1, irctnt-1
          inc1        = inc1+1
          name1(inc1) = periodic_1(isave(ind),jsave(ind))(:2)
      enddo
      inc2 = 0
      do ind = irctnt, natoms
          inc2        = inc2+1
          name2(inc2) = periodic_1(isave(ind),jsave(ind))(:2)
      enddo
      if( nflag(1) == 2 ) distance = 'angstroms           '
      if( nflag(2) == 2 )then
          units = 'millihartree        '
      elseif( nflag(2) .eq. 3 )then
          units = 'ev                  '
      elseif( nflag(2) .eq. 4 )then
          units = 'kcal per mole       '
      elseif( nflag(2) .eq. 5 )then
          units = 'wavenumbers         '
      elseif( nflag(2) .eq. 6 )then
          units = 'kilojoules per mole '
      endif
      cnvrtd  = 1.d0
      cnvrte  = 1.d0
      cnvrtde = 1.d0
      if( nflag(1) == 2 ) cnvrtd = bohr_a
      if( nflag(2) == 2 )then
          cnvrte = cnvrte*htomillh
      elseif( nflag(2) == 3 )then
          cnvrte = cnvrte*htoev
      elseif( nflag(2) == 4 )then
          cnvrte = cnvrte*htokcal
      elseif( nflag(2) == 5 )then
          cnvrte = cnvrte*htowave
      elseif( nflag(2) == 6 )then
          cnvrte = cnvrte*htokj
      endif
      cnvrtde = cnvrte/cnvrtd
      isum    = 0
      do inu=1,25
          isum = isum + nulbl(inu)
      enddo
      ireorder = 0
      if( isum .ne. 0 ) ireorder = 1
end subroutine ancvrt

subroutine cartou
    implicit none
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom=75, natom=25, isurf=5
    integer, parameter :: jsurf=isurf*(isurf+1)/2
    integer            :: ind, i, k
    !   common/infocm/ cartnu(natom,3),indexes(natom),
    !  +               irctnt,natoms,icartr,mder,msurf,ref
    !   common/utilcm/ dgscartnu(natom,3),descartnu(natom,3,isurf),
    !  +               dijcartnu(natom,3,jsurf),cnvrtd,cnvrte,
    !  +               cnvrtde,ireorder,ksdiag,kediag,ksoffd,keoffd
    !   common/usricm/ cart(natom,3),anuzero,
    !  +               nulbl(natom),nflag(20),
    !  +               nasurf(isurf+1,isurf+1),nder
    if( ireorder .eq. 1 )then
        do i=1,natoms
            do j=1,3
                cartnu(nulbl(i),j)=cart(i,j)/cnvrtd
            enddo
        enddo
    else
        do i=1,natoms
            do j=1,3
                cartnu(i,j)=cart(i,j)/cnvrtd
            enddo
        enddo
      endif
end subroutine cartou
 
subroutine carttor
    implicit none
    character(len=75)  :: ref(5)
    integer, parameter :: n3atom=75, natom=25, isurf=5
    integer            :: ind, i, k
    !   common /pt1cm/  r(n3atom), engygs, degsdr(n3atom)
    !   common/infocm/ cartnu(natom,3),indexes(natom),
    !  +               irctnt,natoms,icartr,mder,msurf,ref
    !   common/usricm/ cart(natom,3),anuzero,
    !  +               nulbl(natom),nflag(20),
    !  +               nasurf(isurf+1,isurf+1),nder
    if( icartr == 1 )then
        do i = 1, natoms
            ind      = 3*i-2
            r(ind)   = cartnu(i,1)
            r(ind+1) = cartnu(i,2)
            r(ind+2) = cartnu(i,3)
        enddo
    elseif( icartr == 2 )then
        i = 1                                                       
        do k = 1, natoms-1
            do l = k+1,natoms                                  
               r(i) = sqrt( (cartnu(k,1)-cartnu(l,1))**2 + (cartnu(k,2)-cartnu(l,2))**2 + (cartnu(k,3)-cartnu(l,3))**2 )
               i = i + 1                  
            enddo
         enddo
    elseif( icartr == 3 )then
        r(1) = sqrt( (cartnu(1,1)-cartnu(2,1))**2 +(cartnu(1,2)-cartnu(2,2))**2 +(cartnu(1,3)-cartnu(2,3))**2 )
        r(2) = sqrt( (cartnu(2,1)-cartnu(3,1))**2 +(cartnu(2,2)-cartnu(3,2))**2 +(cartnu(2,3)-cartnu(3,3))**2 )
        r(3) = sqrt( (cartnu(1,1)-cartnu(3,1))**2 +(cartnu(1,2)-cartnu(3,2))**2 +(cartnu(1,3)-cartnu(3,3))**2 )
    elseif( icartr == 4 )then
        flm    = 18.99840d0
        hym    = 1.007825d0
        xcm1   = (hym*cartnu(1,1)+flm*cartnu(2,1))/(flm+hym)
        ycm1   = (hym*cartnu(1,2)+flm*cartnu(2,2))/(flm+hym)
        zcm1   = (hym*cartnu(1,3)+flm*cartnu(2,3))/(flm+hym)
        xcm2   = (hym*cartnu(3,1)+flm*cartnu(4,1))/(flm+hym)
        ycm2   = (hym*cartnu(3,2)+flm*cartnu(4,2))/(flm+hym)
        zcm2   = (hym*cartnu(3,3)+flm*cartnu(4,3))/(flm+hym)
        xcm3   = xcm2-xcm1
        ycm3   = ycm2-ycm1
        zcm3   = zcm2-zcm1
        xrm1   = cartnu(1,1)-xcm1
        yrm1   = cartnu(1,2)-ycm1
        zrm1   = cartnu(1,3)-zcm1
        theta1 = (xrm1*xcm3+yrm1*ycm3+zrm1*zcm3)
        theta1 = theta1/(sqrt(xrm1**2+yrm1**2+zrm1**2))
        theta1 = theta1/(sqrt(xcm3**2+ycm3**2+zcm3**2))
        if( theta1 .gt.  1.0d0 ) theta1=1.0d0
        if( theta1 .lt. -1.0d0 ) theta1=-1.0d0
        theta1 = acos(theta1)
        xrm2   = cartnu(3,1)-xcm2
        yrm2   = cartnu(3,2)-ycm2
        zrm2   = cartnu(3,3)-zcm2
        theta2 = (xrm2*(-xcm3)+yrm2*(-ycm3)+zrm2*(-zcm3))
        theta2 = theta2/(sqrt(xrm2**2+yrm2**2+zrm2**2))
        theta2 = theta2/(sqrt(xcm3**2+ycm3**2+zcm3**2))
        if( theta2 .gt.  1.0d0 ) theta2 = 1.0d0
        if( theta2 .lt. -1.0d0 ) theta2 =-1.0d0
        theta2 = acos(theta2)
        pi     = acos(-1.0d0)
        theta2 = pi-theta2
        q1     = sqrt(xrm1**2+yrm1**2+zrm1**2)
        q2     = sqrt(xrm2**2+yrm2**2+zrm2**2)
        cmm    = (xcm3**2+ycm3**2+zcm3**2)
        cmm    = sqrt(cmm)
        hhd    = (cartnu(1,1)-cartnu(3,1))**2+(cartnu(1,2)-cartnu(3,2))**2+(cartnu(1,3)-cartnu(3,3))**2
        hhd    = sqrt(hhd)
        q      = cmm-q1*cos(theta1)+q2*cos(theta2)
        q3     = sqrt(abs(hhd**2-q**2))
        q1     = q1*sin(theta1)
        q2     = q2*sin(theta2)
        cphi   = (q1**2+q2**2-q3**2)/(2.*q1*q2)
        if( cphi .lt. -1.0d0 ) cphi = -1.0d0
        if( cphi .gt.  1.0d0 ) cphi =  1.0d0
        phi    = acos(cphi)
        2001 format(6f12.8)
        r(1)=sqrt(xcm3**2+ycm3**2+zcm3**2)
        r(2)=(sqrt(xrm1**2+yrm1**2+zrm1**2))*(flm+hym)/flm
        r(3)=(sqrt(xrm2**2+yrm2**2+zrm2**2))*(flm+hym)/flm
        r(4)=theta1
        r(5)=theta2
        r(6)=phi
    elseif( icartr .ne. 0 )then
        write(nflag(18),1000) icartr
        1000    format(2x,'wrong icartr for cartnu; icartr =',i5//)
        stop
    endif
end subroutine carttor
 
subroutine eunitzero
    implicit none
    integer :: i, j
    pengygs = engygs * cnvrte - anuzero
    if( ksdiag /= 0 )then
        do i = ksdiag, kediag
            pengyes(i) = engyes(i) * cnvrte - anuzero
        enddo
    endif
    if( ksoffd /= 0 )then
        do j = ksoffd, keoffd
            pengyij(j) = engyij(j) * cnvrte
        enddo
    endif
end subroutine eunitzero
 
subroutine rtocart
    implicit none 
    real(8) :: ygs(n3atom), yes(n3atom,isurf), yij(n3atom,jsurf)
    integer :: i, j, j1, j2
    if( icartr == 1 )then
        do i = 1, natoms
            ind            = 3*i-2
            dgscartnu(i,1) = degsdr(ind)
            dgscartnu(i,2) = degsdr(ind+1)
            dgscartnu(i,3) = degsdr(ind+2)
            if( ksdiag /= 0 )then
                do j = ksdiag,kediag
                    descartnu(i,1,j) = deesdr(ind,j)
                    descartnu(i,2,j) = deesdr(ind+1,j)
                    descartnu(i,3,j) = deesdr(ind+2,j)
                enddo
            endif
            if( keoffd /= 0 ) then
               do k = ksoffd,keoffd
                  dijcartnu(i,1,k) = deijdr(ind,k)
                  dijcartnu(i,2,k) = deijdr(ind+1,k)
                  dijcartnu(i,3,k) = deijdr(ind+2,k)
               end do
            endif
        enddo
    elseif( icartr == 2 )then
        do i = 1, natoms         
            dgscartnu(i,1) = 0.d0
            dgscartnu(i,2) = 0.d0
            dgscartnu(i,3) = 0.d0
            if( ksdiag .ne. 0 )then
                do j1 = ksdiag, kediag
                  descartnu(i,1,j1) = 0.d0
                  descartnu(i,2,j1) = 0.d0
                  descartnu(i,3,j1) = 0.d0
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do j2 = ksoffd, keoffd
                    dijcartnu(i,1,j2) = 0.d0
                    dijcartnu(i,2,j2) = 0.d0
                    dijcartnu(i,3,j2) = 0.d0
                enddo
            endif
            do j = 1,natoms
                if( j .lt. i )then
                    m1 = natoms*(j-1) - (j*(j-1))/2 + i-j
                elseif( j .gt. i )then
                    m1 = natoms*(i-1) - (i*(i-1))/2 + j-i
                else
                    exit
                    !go to 20
                endif
                y              = degsdr(m1)
                termx          = (cartnu(i,1)-cartnu(j,1))/r(m1)
                termy          = (cartnu(i,2)-cartnu(j,2))/r(m1)
                termz          = (cartnu(i,3)-cartnu(j,3))/r(m1)
                dgscartnu(i,1) = dgscartnu(i,1) + termx*y
                dgscartnu(i,2) = dgscartnu(i,2) + termy*y
                dgscartnu(i,3) = dgscartnu(i,3) + termz*y
                if( ksdiag > 0 )then
                    y = deesdr(m1,j1)
                    do j1 = ksdiag, kediag
                        descartnu(i,1,j1) = descartnu(i,1,j1) + termx*y
                        descartnu(i,2,j1) = descartnu(i,2,j1) + termy*y
                        descartnu(i,3,j1) = descartnu(i,3,j1) + termz*y
                    enddo
                elseif( ksoffd > 0 )then
                    do j2 = ksoffd,keoffd
                        y                 = deijdr(m1,j2)
                        dijcartnu(i,1,j2) = dijcartnu(i,1,j2) + termx*y
                        dijcartnu(i,2,j2) = dijcartnu(i,2,j2) + termy*y
                        dijcartnu(i,3,j2) = dijcartnu(i,3,j2) + termz*y
                    enddo
                endif
                !20 continue
            enddo
        enddo
    elseif( icartr == 3 )then
        do i = 1, natoms
            ygs(i) = degsdr(i)/r(i)
            if( ksdiag .ne. 0 ) then
                do j = ksdiag, kediag
                   yes(i,j) = deesdr(i,j)/r(i)
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do k = ksoffd,keoffd
                    yij(i,k) = deijdr(i,k)/r(i)
                enddo
            endif
        enddo
        do k = 1,3
            term12         =  cartnu(1,k) - cartnu(2,k)
            term23         =  cartnu(2,k) - cartnu(3,k)
            term13         =  cartnu(1,k) - cartnu(3,k)
            dgscartnu(1,k) =  term12*ygs(1) + term13*ygs(3)
            dgscartnu(2,k) = -term12*ygs(1) + term23*ygs(2)
            dgscartnu(3,k) = -term13*ygs(3) - term23*ygs(2)
            if( ksdiag .ne. 0 )then
                do j1 = ksdiag, kediag
                    descartnu(1,k,j1) =  term12*yes(1,j1) + term13*yes(3,j1)
                    descartnu(2,k,j1) = -term12*yes(1,j1) + term23*yes(2,j1)
                    descartnu(3,k,j1) = -term13*yes(3,j1) - term23*yes(2,j1)
                enddo
            endif
            if( ksoffd .ne. 0 )then
                do j2 = ksoffd, keoffd
                    dijcartnu(1,k,j2) =  term12*yij(1,j2) + term13*yij(3,j2)
                    dijcartnu(2,k,j2) = -term12*yij(1,j2) + term23*yij(2,j2)
                    dijcartnu(3,k,j2) = -term13*yij(3,j2) - term23*yij(2,j2)
                enddo
            endif
        enddo
    elseif( icartr .eq. 4 )then
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
        write(nflag(18),'(2x,' wrong icartr for derivative; icartr =',i5//)') icartr
        stop
    endif
end subroutine rtocart
 
subroutine dedcou
    implicit none 
    integer :: i, cnt(2)
    if( ksdiag /= 0 )then
        cnt(1) = ksdiag
        cnt(2) = kediag 
    elseif( ksoffd /= 0 )then
        cnt(1) = ksoffd 
        cnt(2) = keoffd 
    endif 
    do i = 1, natoms
        if( ireorder == 1 )then
            dgscart(i,:) = dgscartnu(nulbl(i),:) * cnvrtde
        else
            dgscart(i,:) = dgscartnu(i,:) * cnvrtde
        endif
        do j = cnt(1), cnt(2)
            if( ireorder == 1 )then
                descart(i,:,j) = descartnu(nulbl(i),:,j) * cnvrtde
            else
                descart(i,:,j) = descartnu(i,:,j) * cnvrtde
            endif
        enddo
    enddo
end subroutine dedcou
