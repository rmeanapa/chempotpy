      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=2
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      integer :: nt
      double precision :: r(1,3), r2(3), v(1)
      double precision :: u11(1), u12(1), u22(1)
      double precision :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
      double precision :: u(nstates,nstates), dudr(nstates,nstates,3)
      double precision :: t(nstates,nstates)
      double precision :: tmpmat(nstates,nstates)
      double precision :: hr(nstates,nstates,3), gr(nstates,3)
      double precision :: tx(9), drdx(3,9)
      double precision :: hx(nstates,nstates,9), gx(nstates,9)
      logical, save :: first_time_data=.true.
      integer :: i, j, k, l

      !initialize 
      u=0.d0
      p=0.d0
      g=0.d0
      d=0.d0

      nt=1
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        tx(j)=x(iatom,idir)
      enddo
      enddo
      ! input cartesian is ClHH
      r(1,1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)/0.529177211
      r(1,2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)/0.529177211
      r(1,3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)/0.529177211

      if(first_time_data) then
      call prepot
      first_time_data=.false.
      endif

      call pot(r,u11,de_u11,1,1)
      call pot(r,u12,de_u12,1,2)
      call pot(r,u22,de_u22,1,3)
      u(1,1)=u11(1)*27.211386
      u(1,2)=u12(1)*27.211386
      u(2,1)=u12(1)*27.211386
      u(2,2)=u22(1)*27.211386
      dudr(1,1,:)=de_u11(:,1)*51.422067
      dudr(1,2,:)=de_u12(:,1)*51.422067
      dudr(2,1,:)=de_u12(:,1)*51.422067
      dudr(2,2,:)=de_u22(:,1)*51.422067
      call diagonalize(nstates,u,p,t)

      do i=1,3
        tmpmat(:,:)=dudr(:,:,i)
        tmpmat=matmul(transpose(t), matmul(tmpmat,t))
        do j=1,nstates
          do k=j+1,nstates
            if ((p(k)-p(j)) > 1.d-8) then
              hr(j,k,i)=tmpmat(j,k)/(p(k)-p(j))
            else
              hr(j,k,i)=tmpmat(j,k)/1.d-8
            endif
            hr(k,j,i)=-hr(j,k,i)
          enddo 
        enddo
        do j=1,nstates
          gr(j,i)=tmpmat(j,j)
        enddo
      enddo
      
      r2=r(1,:)*0.529177211
      call evdrdx(tx, r2, drdx)
      hx=0.d0
      gx=0.d0
      do i=1,9
      do j=1,3
        hx(:,:,i)=hx(:,:,i)+hr(:,:,j)*drdx(j,i)
        gx(:,i)=gx(:,i)+gr(:,j)*drdx(j,i)
      enddo
      enddo
  
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        g(:,iatom,idir)=gx(:,j)
        d(:,:,iatom,idir)=hx(:,:,j)
      enddo
      enddo


      endsubroutine

      subroutine EvdRdX(X,r,drdx)

      integer i,j
      double precision, intent(in) :: X(9), R(3)
      double precision, intent(out) :: dRdX(3,9)

! Initialize dRdX(3,9)
      do i=1,3
        do j=1,9
          dRdX(i,j)=0.0d0
        enddo
      enddo

      dRdX(1,1)=(x(1)-x(4))/r(1)
      dRdX(1,2)=(x(2)-x(5))/r(1)
      dRdX(1,3)=(x(3)-x(6))/r(1)
      dRdX(1,4)=-dRdX(1,1)
      dRdX(1,5)=-dRdX(1,2)
      dRdX(1,6)=-dRdX(1,3)

      dRdX(2,4)=(x(4)-x(7))/r(2)
      dRdX(2,5)=(x(5)-x(8))/r(2)
      dRdX(2,6)=(x(6)-x(9))/r(2)
      dRdX(2,7)=-dRdX(2,4)
      dRdX(2,8)=-dRdX(2,5)
      dRdX(2,9)=-dRdX(2,6)

      dRdX(3,1)=(x(1)-x(7))/r(3)
      dRdX(3,2)=(x(2)-x(8))/r(3)
      dRdX(3,3)=(x(3)-x(9))/r(3)
      dRdX(3,7)=-dRdX(3,1)
      dRdX(3,8)=-dRdX(3,2)
      dRdX(3,9)=-dRdX(3,3)

      endsubroutine

      subroutine diagonalize(n,A_ss,EV_s,U_ss)
      implicit none
      integer, intent(in) :: n
      real*8, intent(in)  :: A_ss(n,n)
      real*8, intent(out) :: U_ss(n,n)
      real*8, intent(out) :: EV_s(n)
      integer :: io,i
      real*8,allocatable :: work_d(:)
      integer, save :: lwork_d
      logical, save :: first_time_diag=.true.
      U_ss=A_ss
      if( first_time_diag )then
        lwork_d= -1
        allocate( work_d(2) )
        call dsyev('V','L',n,U_ss,n,EV_s,work_d,lwork_d,io)
        lwork_d=int(work_d(1))
        deallocate ( work_d )
        first_time_diag=.false.
      endif
      allocate( work_d(lwork_d) )
      call dsyev('V','L',n,U_ss,n,EV_s,work_d,lwork_d,io)
      return
      end subroutine diagonalize

!   System:          NaFH
!   Common name:     NaFH surface fit B
!   Functional form: Modified Extended LEPS with a LEPS correction function.
!   Number of derivatives: 1
!   Number of electronic surfaces: 2
!   Interface: 3-2V
!   References:  A. W. Jasper, M. D. Hack, A. Chakraborty, D. G. Truhlar, 
!                and P. Piecuch, J. Chem. Phys. 115, 7945 (2001).
!   Notes:  This is an improved version of NaFH surface fit A.  The diabatic 
!           coupling has been made to vanish asymptotically and the diatomic 
!           curves have been replaced with experimental curves.
!   Protocol:
!      PREPOT - initializes the potential's variables
!               must be called once before any calls to POT
!      POT    - driver for the evaluation of the energy
!   Units:
!      energies    - hartrees
!      coordinates - bohr
!      derivatives - hartrees/bohr
!   Surfaces:
!      ground and first excited electronic states, diabatic representation
!  This is version 2 of the NaFH potential energy matrix.
!  Details can be found in 
!  "Coupled Quasidiabatic Potential Energy Surfaces for LiFH", 
!  M. D. Hack, A. W. Jasper, D. G. Truhlar, P. Piecuch, prepared 
!  for publication in J. Chem. Phys.
!  This surface is distributed with the NAT code.
!  This surface uses version 1, but the diagonal diabats are corrected with
!  a correction function which is the difference of two LEPS functions.  The
!  surfaces are corrected such that they have the experimental bond energies.
!  The diabatic coupling is cut off in the HF bond direction so that it vanishes 
!  asymptotically.
subroutine prepot
      common /com_para/ g_myrank,g_nprocs
      integer g_myrank,g_nprocs
    call prepot12
    call prepot22
    call prepot11
  100 format(/71('*')/)
end subroutine prepot

!         R     is an array with dimensions R(NT,3)
!               this array contains the geometry(ies) to be computed  (input)
!               Geometries are specified in terms of interatomic distances so
!               that
!               R(i,1) should be the distance between atoms 1 and 2 (A-B), NaH
!               R(i,2) should be the distance between atoms 2 and 3 (B-C), HF
!               R(i,3) should be the distance between atoms 1 and 3 (A-C), NaF
!               where i is an index in the range 1 <= i <= NT.
!         E     is an array with dimensions E(NT)
!               this array contains the diabatic potential energy(ies)  (output)
!         DE    is an array of the derivatives of the diabatic potential
!               energy(ies) with dimensions DE(3,NT) where the second
!               index labels the derivative with respect to the specific
!               component of R (interatomic distance)  (output)
!         NT    is the number of geometries to be computed  (input)
!         NSURF indicates the surface for which the energy should be
!               computed  (input)
!               The parameter NSURF may have the values 1, 2, or 3.
!               These values select a particular diabatic surface:
!               NSURF = 1  diabatic surface which is the lower one
!                          in the reactant channel
!               NSURF = 2  diabatic coupling surface
!               NSURF = 3  diabatic surface which is the upper one
!                          in the reactant channel
subroutine pot( r, e, De, nt, nsurf )
    implicit none
    dimension r(nt,3),e(nt)
    dimension De(3,nt)
    if( nsurf .eq. 1 )then
        call pot11(r, e, De, nt)
    elseif( nsurf .eq. 2) then
        call pot12(r, e, De, nt)
    else
        call pot22(r, e, De, nt)
    endif
end subroutine pot

! U11 surface - lower diabatic surface
      subroutine prepot11
      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt),rdum(3)
      dimension De(3,nt),g_r(3,1,3),Dcfu11(3)
      common /com_para/ g_myrank,g_nprocs
      integer g_myrank,g_nprocs
      save /onetwo/,/uone/,/angular1/,rac2,af1,af2
      common/onetwo/CNF(6),CNF1(5),CHF(6),bc(0:5),
     & ads1,ars1,abs1,amh,amf,amna,amu,una2p,asmall,cevau,
     & DeNF,reNF,DeHF,reHF
      common/uone/y(4,3,4),c2a(4),c2b(4),c2c(4),dip1(4),
     &            adip1(4),rdip1(4),acut1(4),rcut1(4),
     &            dipn(4),al2(4),aln0(4),
     &            dipnn(4),al2n(4),aln0n1(4),aln0n2(4),rn0n(4)
      common/angular1/ac(4)
      rac2=1.d0/dsqrt(2.d0)
      af1 = amh/(amf+amh)
      af2 = amf/(amf+amh)
      if( g_myrank .eq. 0 )then
!        write(6,100)af1,af2
!     Excited NaF surface:
!        write(6,102)(cNF1(i),i=1,5)
!     Ground HF suface:
!        write(6,103)DeHF,reHF,(cHF(i),i=1,6)
!        write(6,104)ads1,ars1,abs1
!        write(6,105)una2p,asmall
!        write(6,106)(y(1,1,j),j=1,4)
!        write(6,107)(y(2,1,j),j=1,4)
!        write(6,108)(y(3,1,j),j=1,4)
!        write(6,109)(y(4,1,j),j=1,4)
!        write(6,110)(y(1,2,j),j=1,4)
!        write(6,111)(y(2,2,j),j=1,4)
!        write(6,112)(y(3,2,j),j=1,4)
!        write(6,113)(y(4,2,j),j=1,4)
!        write(6,114)(y(1,3,j),j=1,4)
!        write(6,115)(y(2,3,j),j=1,4)
!        write(6,116)(y(3,3,j),j=1,4)
!        write(6,117)(y(4,3,j),j=1,4)
!        write(6,118)(c2a(j),j=1,4)
!        write(6,119)(c2b(j),j=1,4)
!        write(6,120)(c2c(j),j=1,4)
!        write(6,121)(dip1(j),j=1,4)
!        write(6,122)(adip1(j),j=1,4)
!        write(6,123)(rdip1(j),j=1,4)
!        write(6,124)(acut1(j),j=1,4)
!        write(6,125)(rcut1(j),j=1,4)
!        write(6,126)(dipn(j),j=1,4)
!        write(6,127)(al2(j),j=1,4)
!        write(6,128)(aln0(j),j=1,4)
!        write(6,129)(dipnn(j),j=1,4)
!        write(6,130)(al2n(j),j=1,4)
!        write(6,131)(aln0n1(j),j=1,4)
!        write(6,132)(aln0n2(j),j=1,4)
!        write(6,133)(rn0n(j),j=1,4)
      endif
return
      entry pot11(r,e,De,nt)
      do 10 i=1,nt
      r1s      = r(i,1)*r(i,1)
      r2s      = r(i,2)*r(i,2)
      r3s      = r(i,3)*r(i,3)
      rnag2    = af1*r1s+af2*r3s-af1*af2*r2s
      Drnag2D1 = 2.d0*af1*r(i,1)
      Drnag2D2 = -2.d0*af1*af2*r(i,2)
      Drnag2D3 = 2.d0*af2*r(i,3)
      ! Cosine of Na-F-H bond angle:
      hlp      = 1.d0/(r(i,2)*r(i,3))
      csb      = 0.5d0*(r2s+r3s-r1s)*hlp
      DcsbD1   = -2.d0*r(i,1)*hlp
      DcsbD2   = 2.d0*(1.d0/r(i,3)-csb/r(i,2))
      DcsbD3   = 2.d0*(1.d0/r(i,2)-csb/r(i,3))
      csb      = 2.d0*(1.d0+csb)
      if( csb .ge. 4.d0 ) csb = 3.99999d0
      l1 = idint(csb)+1
      l2 = l1+1
      if( l1  .eq. 4       ) l2  = l1
      if( csb .le. bc(l1)  ) csb = bc(l1)+1.d-13
      if( csb .ge. bc(l1+1)) csb = bc(l1+1)-1.d-13
      axs1     = dexp(-abs1*(r(i,1)-ars1))
      s1       = ads1*axs1*(axs1-2.d0)
      Ds1D1    = 2.d0*ads1*abs1*axs1*(1.d0-axs1)
      y2       = r(i,2)-reHF
      hlp1     = cHF(1)*dexp(-cHF(2)*y2)
      hlp2     = dexp(-cHF(6)*y2)
      hlp22    = (-1.d0-cHF(1)+cHF(3)*y2+cHF(4)*y2**2+cHF(5)*y2**3)*hlp2
      s2       = DeHF*(hlp1+hlp22)
      Ds2D2    = DeHF*(-cHF(2)*hlp1 - cHF(6)*hlp22+ (cHF(3)+2.d0*cHF(4)*y2+3.d0*cHF(5)*y2**2)*hlp2)
      y3       = r(i,3)-reNF
      hlp1     = dexp(-cNF1(3)*y3)
      hlp11    = (cNF1(1)+cNF1(2)*y3)*hlp1
      hlp2     = cNF1(4)*dexp(-cNF1(5)*y3)
      s3       = hlp11+hlp2
      Ds3D3    = -cNF1(3)*hlp11+cNF1(2)*hlp1-cNF1(5)*hlp2
      e(i)     = 0.d0
      De(1,i)  = 0.d0
      De(2,i)  = 0.d0
      De(3,i)  = 0.d0
      anorm    = 0.d0
      DanormD1 = 0.d0
      DanormD2 = 0.d0
      DanormD3 = 0.d0
      do 20 l = l1, l2
      hlp1   = y(1,1,l)*dexp(-y(2,1,l)*r(i,1))
      hlp2   = y(3,1,l)*dexp(-y(4,1,l)*r(i,1))
      t1     = hlp1+hlp2+s1
      Dt1D1  = -y(2,1,l)*hlp1-y(4,1,l)*hlp2+Ds1D1
      hlp1   = y(1,2,l)*dexp(-y(2,2,l)*r(i,2))
      hlp2   = y(3,2,l)*dexp(-y(4,2,l)*r(i,2))
      t2     = hlp1+hlp2+s2
      Dt2D2  = -y(2,2,l)*hlp1-y(4,2,l)*hlp2+Ds2D2
      hlp1   = y(1,3,l)*dexp(-y(2,3,l)*r(i,3))
      hlp2   = y(3,3,l)*dexp(-y(4,3,l)*r(i,3))
      t3     = hlp1+hlp2+s3
      Dt3D3  = -y(2,3,l)*hlp1-y(4,3,l)*hlp2+Ds3D3
      coul1  = 0.5d0*(s1+t1)
      Dcl1D1 = 0.5d0*(Ds1D1+Dt1D1)
      coul2  = 0.5d0*(s2+t2)
      Dcl2D2 = 0.5d0*(Ds2D2+Dt2D2)
      coul3  = 0.5d0*(s3+t3)
      Dcl3D3 = 0.5d0*(Ds3D3+Dt3D3)
      exch1  = 0.5d0*(s1-t1)
      Dex1D1 = 0.5d0*(Ds1D1-Dt1D1)
      exch2  = 0.5d0*(s2-t2)
      Dex2D2 = 0.5d0*(Ds2D2-Dt2D2)
      exch3  = 0.5d0*(s3-t3)
      Dex3D3 = 0.5d0*(Ds3D3-Dt3D3)
      w      = (exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
      DwD1   = 2.d0*(2.d0*exch1-exch2-exch3)*Dex1D1
      DwD2   = 2.d0*(2.d0*exch2-exch1-exch3)*Dex2D2
      DwD3   = 2.d0*(2.d0*exch3-exch1-exch2)*Dex3D3
      cplg2  = c2a(l)*dexp(-c2b(l)*w-c2c(l)*(r(i,1)+r(i,2)+r(i,3)))
      DcD1   = (-c2b(l)*DwD1 - c2c(l))*cplg2
      DcD2   = (-c2b(l)*DwD2 - c2c(l))*cplg2
      DcD3   = (-c2b(l)*DwD3 - c2c(l))*cplg2
      rootw  = rac2*dsqrt(w+cplg2*cplg2)
      surf   = (coul1+coul2+coul3+DeHF-rootw)
! check for zero by zero division (a case when all three atoms are asymptotically
! separated)
      if( (DwD1+2.d0*cplg2*DcD1) .eq. 0.0d0 )then
         DsurfD1=Dcl1D1
      else
         DsurfD1=(Dcl1D1-0.25d0*(DwD1+2.d0*cplg2*DcD1)/rootw)
      endif
      if( (DwD2+2.d0*cplg2*DcD2) .eq. 0.0d0 )then
         DsurfD2=Dcl2D2
      else
         DsurfD2=(Dcl2D2-0.25d0*(DwD2+2.d0*cplg2*DcD2)/rootw)
      endif
      if( (DwD3+2.d0*cplg2*DcD3) .eq. 0.0d0 )then
         DsurfD3=Dcl3D3
      else
         DsurfD3=(Dcl3D3-0.25d0*(DwD3+2.d0*cplg2*DcD3)/rootw)
      endif
      hlp1    = adip1(l)*(rnag2-rdip1(l))
      hlp2    = acut1(l)*(r(i,2)-rcut1(l))
      hlp11   = 0.5d0/(dcosh(hlp1)**2)
      hlp22   = 0.5d0/(dcosh(hlp2)**2)
      hlp1    = 0.5d0*(dtanh(hlp1)-1.d0)
      hlp2    = 0.5d0*(1.d0-dtanh(hlp2))
      surf    = surf+dip1(l)*hlp1*hlp2
      DsurfD1 = DsurfD1+dip1(l)*adip1(l)*Drnag2D1*hlp11*hlp2
      DsurfD2 = DsurfD2+dip1(l)*adip1(l)*Drnag2D2*hlp11*hlp2
     &               -dip1(l)*hlp1*acut1(l)*hlp22
      DsurfD3 = DsurfD3+dip1(l)*adip1(l)*Drnag2D3*hlp11*hlp2
      hlp1    = al2(l)*(r(i,2)-reHF)
      hlp2    = aln0(l)*rnag2
      hlp11   = dtanh(hlp1)
      hlp22   = dtanh(hlp2)
      hlp1    = 1.d0/dcosh(hlp1)
      hlp2    = 1.d0/dcosh(hlp2)
      surf    = surf - dipn(l)*hlp1*hlp2
      DsurfD1 = DsurfD1 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*Drnag2D1
      DsurfD2 = DsurfD2 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*Drnag2D2
     &                + dipn(l)*hlp2*hlp1*hlp11*al2(l)
      DsurfD3 = DsurfD3 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*Drnag2D3
      hlp1    = al2n(l)*(r(i,2)-reHF)
      hlp2    = aln0n1(l)*(rnag2-rn0n(l)**2)
      ! try to prevent an overflow
      if( hlp2 .lt. 709.0d0 ) then
          hlp2 = dexp(hlp2)
      else
          hlp2 = 1.0D308
      endif
      hlp3     = dexp(-aln0n2(l)*(rnag2-rn0n(l)**2))
      hlp23    = 1.d0/(hlp2+hlp3)**2
      hlp11    = dtanh(hlp1)
      hlp1     = 1.d0/dcosh(hlp1)
      surf     = surf - 2.d0*dipnn(l)*hlp1/(hlp2+hlp3)
      DsurfD1  = DsurfD1 + 2.d0*dipnn(l)*hlp1*(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*Drnag2D1
      DsurfD2  = DsurfD2 + 2.d0*dipnn(l)*hlp1*(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*Drnag2D2+ 2.d0*dipnn(l)*hlp1*hlp11*al2n(l)/(hlp2+hlp3)
      DsurfD3  = DsurfD3 + 2.d0*dipnn(l)*hlp1*(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*Drnag2D3
      hlp      = 1.d0/((1.d0-(csb-bc(l))/(bc(l-1)-bc(l)))* (1.d0-(csb-bc(l))/(bc(l+1)-bc(l))))
      hlp1     = (csb-bc(l))**2*hlp
      cog      = dexp(-ac(l)*hlp1)
      Dcog     = -cog*ac(l)*(2.d0*(csb-bc(l))*hlp+hlp1*(1.d0/(bc(l-1)-csb)+1.d0/(bc(l+1)-csb)))
      DcogD1   = Dcog*DcsbD1
      DcogD2   = Dcog*DcsbD2
      DcogD3   = Dcog*DcsbD3
      e(i)     = e(i) + surf*cog
      De(1,i)  = De(1,i) + DsurfD1*cog+surf*DcogD1
      De(2,i)  = De(2,i) + DsurfD2*cog+surf*DcogD2
      De(3,i)  = De(3,i) + DsurfD3*cog+surf*DcogD3
      anorm    = anorm + cog
      DanormD1 = DanormD1 + DcogD1
      DanormD2 = DanormD2 + DcogD2
      DanormD3 = DanormD3 + DcogD3
   20 continue
!     if (csb.gt.3.99d0) write(6,77)l1,l2,csb,r(i,3),
!    &                         r(i,2),e(i),anorm
      hlp     = 1.d0/anorm
      e(i)    = e(i)*hlp
      De(1,i) = (De(1,i) - e(i)*DanormD1)*hlp
      De(2,i) = (De(2,i) - e(i)*DanormD2)*hlp
      De(3,i) = (De(3,i) - e(i)*DanormD3)*hlp
c  77 format(2(2x,i3),5(2x,e16.8))
      do j=1,3
          rdum(j) = r(i,j)
      enddo
      call lepscfu11(3,rdum,g_r,3,cfu11,Dcfu11,3)
      e(i)    = e(i) + cfu11
      De(1,i) = De(1,i) + Dcfu11(1)
      De(2,i) = De(2,i) + Dcfu11(2)
      De(3,i) = De(3,i) + Dcfu11(3)

      e(i) = e(i)*cevau
      De(1,i) = De(1,i)*cevau
      De(2,i) = De(2,i)*cevau
      De(3,i) = De(3,i)*cevau
10    continue
  100 format(/26('='),' NaFH - U11 diabat ',26('='),/
     &2x,'mH/(mF+mH) = ',f13.5,2x,'mF/(mF+mH) = ',f13.5)
  102 format(2x,'cNaF1 = ',9x,3e14.6/19x,2e14.6)
  103 format(2x,'DeHF = ',10x,e14.6/2x,'reHF = ',10x,e14.6
     & /2x,'cHF = ',11x,3e14.6/19x,3e14.6)
  104 format(2x,'ads1,ars1,abs1 = ',3e14.6)
  105 format(2x,'una2p,asmall = ',2x,2e14.6)
  106 format(2x,'y(1,1,i) = ',2x,4e13.5)
  107 format(2x,'y(2,1,i) = ',2x,4e13.5)
  108 format(2x,'y(3,1,i) = ',2x,4e13.5)
  109 format(2x,'y(4,1,i) = ',2x,4e13.5)
  110 format(2x,'y(1,2,i) = ',2x,4e13.5)
  111 format(2x,'y(2,2,i) = ',2x,4e13.5)
  112 format(2x,'y(3,2,i) = ',2x,4e13.5)
  113 format(2x,'y(4,2,i) = ',2x,4e13.5)
  114 format(2x,'y(1,3,i) = ',2x,4e13.5)
  115 format(2x,'y(2,3,i) = ',2x,4e13.5)
  116 format(2x,'y(3,3,i) = ',2x,4e13.5)
  117 format(2x,'y(4,3,i) = ',2x,4e13.5)
  118 format(2x,'c2a = ',7x,4e13.5)
  119 format(2x,'c2b = ',7x,4e13.5)
  120 format(2x,'c2c = ',7x,4e13.5)
  121 format(2x,'dip1 = ',6x,4e13.5)
  122 format(2x,'adip1 = ',5x,4e13.5)
  123 format(2x,'rdip1 = ',5x,4e13.5)
  124 format(2x,'acut1 = ',5x,4e13.5)
  125 format(2x,'rcut1 = ',5x,4e13.5)
  126 format(2x,'dipn = ',6x,4e13.5)
  127 format(2x,'al2 = ',7x,4e13.5)
  128 format(2x,'aln0 = ',6x,4e13.5)
  129 format(2x,'dipnn = ',5x,4e13.5)
  130 format(2x,'al2n = ',6x,4e13.5)
  131 format(2x,'aln0n1 = ',4x,4e13.5)
  132 format(2x,'aln0n2 = ',4x,4e13.5)
  133 format(2x,'rn0n = ',6x,4e13.5
     &    /71('=')/)
      return
      end

!  U12 surface - coupling surface
      subroutine prepot12
      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt)
      dimension De(3,nt)
      common /com_para/ g_myrank,g_nprocs
      integer g_myrank,g_nprocs

      save /onetwo/,/ucoupl/,/ucoupl1/,/u12nah/,af1,af2

      common/onetwo/CNF(6),CNF1(5),CHF(6),bc(0:5),
     & ads1,ars1,abs1,amh,amf,amna,amu,una2p,asmall,cevau,
     & DeNF,reNF,DeHF,reHF
      common/ucoupl/ahf(3),anaf(3),anah(3),r21(3),r23(3),r13(3)
      common/ucoupl1/ccNF(4),ccNHa(4),ccNHi(4),ccNHi1(4)
      common/u12nah/xpNH,xpNH1
      dimension ccNH(4)
      dimension DccNHD2(4),DccNHD3(4)

      af1 = amh/(amf+amh)
      af2 = amf/(amf+amh)
      if( g_myrank .eq. 0 )then
c         write(6,100)af1,af2

C     NaF coupling:
c         write(6,101)(ccNF(i),i=1,4)
c         write(6,102)(ccNHa(i),i=1,4)
c         write(6,103)(ccNHi(i),i=1,4)
c         write(6,104)(ccNHi1(i),i=1,4)
c         write(6,105)xpNH,xpNH1
c         write(6,106)(ahf(i),i=1,3)
c         write(6,107)(anaf(i),i=1,3)
c         write(6,108)(anah(i),i=1,3)
cj         write(6,109)(r21(i),i=1,3)
c         write(6,110)(r23(i),i=1,3)
c         write(6,111)(r13(i),i=1,3)
      endif

      return
c
      entry pot12(r,e,De,nt)

      do 10 i=1,nt

      hlp1 = dexp(-xpNH*r(i,2)**2)
      hlp2 = dexp(-xpNH1*r(i,3)**2)
      ccNH(1) = ccNHa(1)+ccNHi(1)*hlp1+ccNHi1(1)*hlp2
      ccNH(2) = ccNHa(2)+ccNHi(2)*hlp1+ccNHi1(2)*hlp2
      ccNH(3) = ccNHa(3)+ccNHi(3)*hlp1+ccNHi1(3)*hlp2
      ccNH(4) = ccNHa(4)+ccNHi(4)*hlp1+ccNHi1(4)*hlp2

      hlp1 = - 2.d0*xpNH*r(i,2)*hlp1
      hlp2 = - 2.d0*xpNH1*r(i,3)*hlp2
      DccNHD2(1) = ccNHi(1)*hlp1
      DccNHD3(1) = ccNHi1(1)*hlp2
      DccNHD2(2) = ccNHi(2)*hlp1
      DccNHD3(2) = ccNHi1(2)*hlp2
      DccNHD2(3) = ccNHi(3)*hlp1
      DccNHD3(3) = ccNHi1(3)*hlp2
      DccNHD2(4) = ccNHi(4)*hlp1
      DccNHD3(4) = ccNHi1(4)*hlp2

      hlp = ccNF(3)/(1.d0+ccNF(4)*r(i,3))
      unaf = dexp(-ccNF(1)-ccNF(2)*r(i,3)- hlp*r(i,3))
      unaf = unaf*27.2113961d0
      DunafD3 = -unaf*(ccNF(2)+hlp-ccNF(4)*hlp*hlp*r(i,3)/ccNF(3))


      hlp = 1.d0/(1.d0+ccNH(4)*r(i,1))
      unah = dexp(-ccNH(1)-ccNH(2)*r(i,1)- ccNH(3)*hlp*r(i,1))
      DunahD1 = -unah*(ccNH(2)+ccNH(3)*hlp
     &        -ccNH(4)*hlp*hlp*ccNH(3)*r(i,1))
      DunahD2 = -unah*(DccNHD2(1)+DccNHD2(2)*r(i,1)
     &+ DccNHD2(3)*r(i,1)*hlp-ccNH(3)*r(i,1)*r(i,1)*dccNHD2(4)*hlp*hlp)
      DunahD3 = -unah*(DccNHD3(1)+DccNHD3(2)*r(i,1)
     &+ DccNHD3(3)*r(i,1)*hlp-ccNH(3)*r(i,1)*r(i,1)*dccNHD3(4)*hlp*hlp)
C-----------------------------------------------------------------------

      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)
      rnag2=af1*r1s+af2*r3s-af1*af2*r2s
      Drnag2D1 = 2.d0*af1*r(i,1)
      Drnag2D2 = -2.d0*af1*af2*r(i,2)
      Drnag2D3 = 2.d0*af2*r(i,3)
      hlp2 = rnag2+asmall
      hlp = dsqrt(hlp2)
      cstsq=0.5d0*(r3s-r1s+(af2-af1)*r2s)
     &               /(r(i,2)*hlp)
      DcstsqD1 = -2.d0*r(i,1)/(r(i,2)*hlp)
     &         - cstsq*Drnag2D1/hlp2
      DcstsqD2 = 2.d0*((af2-af1)/hlp-cstsq/r(i,2))
     &         - cstsq*Drnag2D2/hlp2
      DcstsqD3 = 2.d0*r(i,3)/(r(i,2)*hlp)
     &         - cstsq*Drnag2D3/hlp2
      cstsq = 2.0d0*(1.d0+cstsq)

            hlp = r21(2)-r21(1)-r21(3)
            r21t = r21(1)+hlp*cstsq
     &                   +r21(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*r21(3)*cstsq
            Dr21tD1 = hlp1*DcstsqD1
            Dr21tD2 = hlp1*DcstsqD2
            Dr21tD3 = hlp1*DcstsqD3

            hlp = r23(2)-r23(1)-r23(3)
            r23t = r23(1)+hlp*cstsq
     &                   +r23(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*r23(3)*cstsq
            Dr23tD1 = hlp1*DcstsqD1
            Dr23tD2 = hlp1*DcstsqD2
            Dr23tD3 = hlp1*DcstsqD3

            hlp = r13(2)-r13(1)-r13(3)
            r13t = r13(1)+hlp*cstsq
     &                   +r13(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*r13(3)*cstsq
            Dr13tD1 = hlp1*DcstsqD1
            Dr13tD2 = hlp1*DcstsqD2
            Dr13tD3 = hlp1*DcstsqD3

            hlp = ahf(2)-ahf(1)-ahf(3)
            ahft = ahf(1)+hlp*cstsq
     &                   +ahf(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*ahf(3)*cstsq
            DahftD1 = hlp1*DcstsqD1
            DahftD2 = hlp1*DcstsqD2
            DahftD3 = hlp1*DcstsqD3

            hlp = anaf(2)-anaf(1)-anaf(3)
            anaft = anaf(1)+hlp*cstsq
     &                   +anaf(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*anaf(3)*cstsq
            DanaftD1 = hlp1*DcstsqD1
            DanaftD2 = hlp1*DcstsqD2
            DanaftD3 = hlp1*DcstsqD3

            hlp = anah(2)-anah(1)-anah(3)
            anaht = anah(1)+hlp*cstsq
     &                   +anah(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*anah(3)*cstsq
            DanahtD1 = hlp1*DcstsqD1
            DanahtD2 = hlp1*DcstsqD2
            DanahtD3 = hlp1*DcstsqD3
 
      hlp = anaht*r(i,1)-anaft*r(i,3)-r13t
      sw1 = 0.5d0*(1.d0+dtanh(hlp))
      hlp1 = 0.5d0/(dcosh(hlp)**2)
      Dsw1D1 = hlp1*(anaht+DanahtD1*r(i,1)-DanaftD1*r(i,3)-Dr13tD1)
      Dsw1D2 = hlp1*(DanahtD2*r(i,1)-DanaftD2*r(i,3)-Dr13tD2)
      Dsw1D3 = hlp1*(DanahtD3*r(i,1)-anaft-DanaftD3*r(i,3)-Dr13tD3)

      hlp = ahft*r(i,2)-anaht*r(i,1)-r21t
      sw21 = 0.5d0*(1.d0+dtanh(hlp))
      hlp1 = 0.5d0/(dcosh(hlp)**2)
      Dsw21D1 = hlp1*(DahftD1*r(i,2)-anaht-DanahtD1*r(i,1)-Dr21tD1)
      Dsw21D2 = hlp1*(ahft+DahftD2*r(i,2)-DanahtD2*r(i,1)-Dr21tD2)
      Dsw21D3 = hlp1*(DahftD3*r(i,2)-DanahtD3*r(i,1)-Dr21tD3)

      hlp = ahft*r(i,2)-anaft*r(i,3)-r23t
      sw23 = 0.5d0*(1.d0+dtanh(hlp))
      hlp1 = 0.5d0/(dcosh(hlp)**2)
      Dsw23D1 = hlp1*(DahftD1*r(i,2)-DanaftD1*r(i,3)-Dr23tD1)
      Dsw23D2 = hlp1*(ahft+DahftD2*r(i,2)-DanaftD2*r(i,3)-Dr23tD2)
      Dsw23D3 = hlp1*(DahftD3*r(i,2)-anaft-DanaftD3*r(i,3)-Dr23tD3)

      e(i) = unaf*sw1*sw23+unah*(1.d0-sw1)*sw21
      De(1,i) = unaf*Dsw1D1*sw23+unaf*sw1*Dsw23D1
     &        + DunahD1*(1.d0-sw1)*sw21-unah*Dsw1D1*sw21
     &        + unah*(1.d0-sw1)*Dsw21D1
      De(2,i) = unaf*Dsw1D2*sw23+unaf*sw1*Dsw23D2
     &        + DunahD2*(1.d0-sw1)*sw21-unah*Dsw1D2*sw21
     &        + unah*(1.d0-sw1)*Dsw21D2
      De(3,i) = DunafD3*sw1*sw23+unaf*Dsw1D3*sw23+unaf*sw1*Dsw23D3
     &        + DunahD3*(1.d0-sw1)*sw21-unah*Dsw1D3*sw21
     &        + unah*(1.d0-sw1)*Dsw21D3

c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     correction function
      cf_p = 2.d0
      cf_r0 = 3.5d0
      swit = 0.5d0*(1.d0+dtanh(cf_p*(cf_r0-r(i,2))))
      Dswit1 = 0.d0
      Dswit2 = -0.5d0*cf_p*1.d0/(dcosh(cf_p*(-r(i,2)+cf_r0))**2)
      Dswit3 = 0.d0

      De(1,i)=e(i)*Dswit1+swit*De(1,i)
      De(2,i)=e(i)*Dswit2+swit*De(2,i)
      De(3,i)=e(i)*Dswit3+swit*De(3,i)
      e(i) = e(i)*swit
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      e(i) = e(i)*cevau
      De(1,i) = De(1,i)*cevau
      De(2,i) = De(2,i)*cevau
      De(3,i) = De(3,i)*cevau
10    continue

      return

  100 format(/25('='),' NaFH - U12 coupling ',25('='),/
     &2x,'mF/(mF+mH) = ',f13.5,2x,'mH/(mF+mH) = ',f13.5)
  101 format(2x,'ccNaF = ',5x,4e14.6/10x,2e14.6)
  102 format(2x,'ccNaHa = ',4x,4e14.6/11x,2e14.6)
  103 format(2x,'ccNaHi = ',4x,4e14.6/11x,2e14.6)
  104 format(2x,'ccNaHi1 = ',3x,4e14.6/12x,2e14.6)
  105 format(2x,'xpNH,xpNH1 = ',2e14.6)
  106 format(2x,'ahf = ',7x,3e14.6)
  107 format(2x,'anaf = ',6x,3e14.6)
  108 format(2x,'anah = ',6x,3e14.6)
  109 format(2x,'r21 = ',7x,3e14.6)
  110 format(2x,'r23 = ',7x,3e14.6)
  111 format(2x,'r13 = ',7x,3e14.6
     &    /71('=')/)

      end
c **********************************************************************
c
c  U22 surface - upper diabatic surface
c
      subroutine prepot22
      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt),rdum(3)
      dimension De(3,nt),g_r(3,1,3),Dcfu22(3)
      common /com_para/ g_myrank,g_nprocs
      integer g_myrank,g_nprocs

      save /onetwo/,/utwo/,/angular2/,rac2,af1,af2

      common/onetwo/CNF(6),CNF1(5),CHF(6),bc(0:5),
     & ads1,ars1,abs1,amh,amf,amna,amu,una2p,asmall,cevau,
     & DeNF,reNF,DeHF,reHF
      common/utwo/y(4,3,4),c2a(4),c2b(4),c2c(4),r20i(4),
     &      r20a,al20i(4),al20i1(4),al20a,psi,rsi,al2a(4),
     &      rr0a(4),dip2(4),adip2(4),rdip2(4),acut2(4),
     &      rcut2(4),ry(4),aly(4),dy(4),dy1(4),dy2(4),
     &      dip2n(4),adip2n(4),rdip2n(4),acut2n(4),rcut2n(4)
      common/angular2/ac(4)
C-----------------------------------------------------------------------

      rac2=1.d0/dsqrt(2.d0)

      af1 = amh/(amf+amh)
      af2 = amf/(amf+amh)
      if( g_myrank .eq. 0 )then
c         write(6,100)af1,af2
C----------------------------------
C     Ground NaF surface
c         write(6,101)DeNF,reNF,(cNF(i),i=1,6)
C----------------------------------
C     Ground HF suface:
c         write(6,103)DeHF,reHF,(cHF(i),i=1,6)
C----------------------------------
c         write(6,105)una2p,asmall
c         write(6,106)(y(1,1,j),j=1,4)
c         write(6,107)(y(2,1,j),j=1,4)
c         write(6,108)(y(3,1,j),j=1,4)
c         write(6,109)(y(4,1,j),j=1,4)
c         write(6,110)(y(1,2,j),j=1,4)
c         write(6,111)(y(2,2,j),j=1,4)
c         write(6,112)(y(3,2,j),j=1,4)
c         write(6,113)(y(4,2,j),j=1,4)
c         write(6,114)(y(1,3,j),j=1,4)
c         write(6,115)(y(2,3,j),j=1,4)
c         write(6,116)(y(3,3,j),j=1,4)
c         write(6,117)(y(4,3,j),j=1,4)
c         write(6,118)(c2a(j),j=1,4)
c         write(6,119)(c2b(j),j=1,4)
c         write(6,120)(c2c(j),j=1,4)
c         write(6,121)(r20i(j),j=1,4)
c         write(6,122)r20a
c         write(6,123)(al20i(j),j=1,4)
c         write(6,124)(al20i1(j),j=1,4)
c         write(6,125)al20a
c         write(6,126)(al2a(j),j=1,4)
c         write(6,127)(rr0a(j),j=1,4)
c         write(6,128)psi
c         write(6,129)rsi
c         write(6,130)(dip2(j),j=1,4)
c         write(6,131)(adip2(j),j=1,4)
c         write(6,132)(rdip2(j),j=1,4)
c         write(6,133)(acut2(j),j=1,4)
c         write(6,134)(rcut2(j),j=1,4)
c         write(6,135)(ry(j),j=1,4)
c         write(6,136)(aly(j),j=1,4)
c         write(6,137)(dy(j),j=1,4)
c         write(6,138)(dip2n(j),j=1,4)
c         write(6,139)(adip2n(j),j=1,4)
c         write(6,140)(rdip2n(j),j=1,4)
c         write(6,141)(acut2n(j),j=1,4)
c         write(6,142)(rcut2n(j),j=1,4)
      endif
C
      return
c
      entry pot22(r,e,De,nt)
C
      do 10 i=1,nt

      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)
      rnag2=af1*r1s+af2*r3s-af1*af2*r2s
      Drnag2D1 = 2.d0*af1*r(i,1)
      Drnag2D2 = -2.d0*af1*af2*r(i,2)
      Drnag2D3 = 2.d0*af2*r(i,3)

C     Cosine of Na-F-H bond angle:
      hlp = 1.d0/(r(i,2)*r(i,3))
      csb = 0.5d0*(r2s+r3s-r1s)*hlp
      DcsbD1 = -2.d0*r(i,1)*hlp
      DcsbD2 = 2.d0*(1.d0/r(i,3)-csb/r(i,2))
      DcsbD3 = 2.d0*(1.d0/r(i,2)-csb/r(i,3))
      csb = 2.d0*(1.d0+csb)

      if (csb.ge.4.d0) csb = 3.99999d0
      l1 = idint(csb)+1
      l2 = l1+1
      if (l1.eq.4) l2=l1
      if (csb.le.bc(l1)) csb = bc(l1)+1.d-13
      if (csb.ge.bc(l1+1)) csb = bc(l1+1)-1.d-13

C-----------------------------------------------------------------------
      s1=0.d0
      Ds1D1 = 0.d0
C-----------------------------------------------------------------------
      y2 = r(i,2)-reHF
      hlp1 = cHF(1)*dexp(-cHF(2)*y2)
      hlp2 = dexp(-cHF(6)*y2)
      hlp22 = (-1.d0-cHF(1)+cHF(3)*y2+cHF(4)*y2**2+cHF(5)*y2**3)*hlp2 
      s2a = DeHF*(hlp1+hlp22) + una2p
      Ds2aD2 = DeHF*(-cHF(2)*hlp1-cHF(6)*hlp22
     &       + (cHF(3)+2.d0*cHF(4)*y2+3.d0*cHF(5)*y2**2)*hlp2)
C-----------------------------------------------------------------------
      y3 = r(i,3)-reNF
      hlp1 = cNF(1)*dexp(-cNF(2)*y3)
      hlp2 = dexp(-cNF(6)*y3)
      hlp22 = (-1.d0-cNF(1)+cNF(3)*y3+cNF(4)*y3*y3+cNF(5)*y3**3)*hlp2
      s3 = DeNF*(hlp1+hlp22)
      Ds3D3 = DeNF*(-cNF(2)*hlp1-cNF(6)*hlp22
     &       + (cNF(3)+2.d0*cNF(4)*y3+3.d0*cNF(5)*y3**2)*hlp2)
C-----------------------------------------------------------------------

      e(i) = 0.d0
      De(1,i) = 0.d0
      De(2,i) = 0.d0
      De(3,i) = 0.d0
      anorm = 0.d0
      DanormD1 = 0.d0
      DanormD2 = 0.d0
      DanormD3 = 0.d0

      do 20 l = l1, l2

      hlp1 = y(1,1,l)*dexp(-y(2,1,l)*r(i,1))
      hlp2 = y(3,1,l)*dexp(-y(4,1,l)*r(i,1))
      t1=hlp1+hlp2+s1
      Dt1D1 = -y(2,1,l)*hlp1-y(4,1,l)*hlp2+Ds1D1

      hlp = al2a(l)*(r(i,3)+r(i,1)-rr0a(l))
      hlp1 = al2a(l)/(dcosh(hlp)**2)
      hlp = 0.5d0*(1.d0+dtanh(hlp))
      hlp2 = r(i,3)+r(i,1)-rsi
      if (hlp2.gt.0.d0) then
c     write(6,*)'R1+R3=',(r(i,1)+r(i,3))
      ctanh = dexp(-psi/hlp2)
      dctanh = psi*ctanh/(r(i,3)+r(i,1)-rsi)**2
      else
      ctanh = 0.d0
      dctanh = 0.d0
      end if
      r20 = r20i(l)+(r20a-r20i(l))*hlp
      Dr20 = 0.5d0*(r20a-r20i(l))*hlp1
      al20 = al20i(l)+(al20i1(l)-al20i(l))*hlp
     &     + (al20a-al20i1(l))*ctanh
      Dal20 = 0.5d0*(al20i1(l)-al20i(l))*hlp1
     &      + (al20a-al20i1(l))*dctanh
   
      hlp1 = y(1,2,l)*dexp(-y(2,2,l)*r(i,2))
      hlp2 =  dexp(-(y(4,2,l)+dy2(l)*ctanh)*r(i,2))
      hlp3 = (y(3,2,l)+dy1(l)*ctanh)*hlp2
      t2=hlp1+hlp3
      Dt2D1=dy1(l)*dctanh*hlp2 - dy2(l)*dctanh*r(i,2)*hlp3
      Dt2D2=-y(2,2,l)*hlp1-(y(4,2,l)+dy2(l)*ctanh)*hlp3
      Dt2D3=Dt2D1

      hlp = 0.5d0*(1.d0+dtanh(al20*(r(i,2)-r20)))
      s2=s2a+(t2-s2a)*hlp
      hlp1 = 0.5d0/(dcosh(al20*(r(i,2)-r20))**2)
      Ds2D1 = (t2-s2a)*hlp1*(Dal20*(r(i,2)-r20)-al20*Dr20)
      Ds2D2 = Ds2aD2+(Dt2D2-Ds2aD2)*hlp
     &      + al20*hlp1*(t2-s2a)
      Ds2D3 = Ds2D1

      hlp = aly(l)*(r(i,3)-ry(l))
      cy=y(2,3,l)+0.5d0*dy(l)*(1.d0+dtanh(hlp))
      dcyD3 = 0.5d0*dy(l)*aly(l)/(dcosh(hlp)**2)
      hlp1 = y(1,3,l)*dexp(-cy*r(i,3))
      hlp2 = y(3,3,l)*dexp(-y(4,3,l)*r(i,3))
      t3 = hlp1+hlp2
      Dt3D3=-(cy+dcyD3*r(i,3))*hlp1-y(4,3,l)*hlp2

      coul1=0.5d0*(s1+t1)
      Dcl1D1=0.5d0*(Ds1D1+Dt1D1)
      coul2=0.5d0*(s2+t2)
      Dcl2D1=0.5d0*(Ds2D1+Dt2D1)
      Dcl2D2=0.5d0*(Ds2D2+Dt2D2)
      Dcl2D3=0.5d0*(Ds2D3+Dt2D3)
      coul3=0.5d0*(s3+t3)
      Dcl3D3=0.5d0*(Ds3D3+Dt3D3)

      exch1=0.5d0*(s1-t1)
      Dex1D1=0.5d0*(Ds1D1-Dt1D1)
      exch2=0.5d0*(s2-t2)
      Dex2D1=0.5d0*(Ds2D1-Dt2D1)
      Dex2D2=0.5d0*(Ds2D2-Dt2D2)
      Dex2D3=0.5d0*(Ds2D3-Dt2D3)
      exch3=0.5d0*(s3-t3)
      Dex3D3=0.5d0*(Ds3D3-Dt3D3)

      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
      DwD1 = 2.d0*(exch1-exch2)*(Dex1D1-Dex2D1)
     &     + 2.d0*(exch2-exch3)*Dex2D1
     &     - 2.d0*(exch3-exch1)*Dex1D1
      DwD2 = - 2.d0*(exch1-exch2)*Dex2D2
     &     + 2.d0*(exch2-exch3)*Dex2D2
      DwD3 = - 2.d0*(exch1-exch2)*Dex2D3
     &     + 2.d0*(exch2-exch3)*(Dex2D3-Dex3D3)
     &     + 2.d0*(exch3-exch1)*Dex3D3

      cplg2=c2a(l)*dexp(-c2b(l)*w-c2c(l)*(r(i,1)+r(i,2)+r(i,3)))
      DcD1 = (-c2b(l)*DwD1 - c2c(l))*cplg2
      DcD2 = (-c2b(l)*DwD2 - c2c(l))*cplg2
      DcD3 = (-c2b(l)*DwD3 - c2c(l))*cplg2

      rootw = rac2*dsqrt(w+cplg2*cplg2)
      surf=(coul1+coul2+coul3+DeHF-rootw)
      DsurfD1=(Dcl1D1+Dcl2D1-0.25d0*(DwD1+2.d0*cplg2*DcD1)/rootw)
      DsurfD2=(Dcl2D2-0.25d0*(DwD2+2.d0*cplg2*DcD2)/rootw)
      DsurfD3=(Dcl2D3+Dcl3D3-0.25d0*(DwD3+2.d0*cplg2*DcD3)/rootw)

      hlp1 = acut2(l)*(r(i,2)-rcut2(l))
      hlp2 = adip2(l)*(rnag2-rdip2(l)**2)
      hlp11 = dtanh(hlp1)
      hlp22 = dtanh(hlp2)
      hlp1 = 1.d0/dcosh(hlp1)
      hlp2 = 1.d0/dcosh(hlp2)
      surf = surf - dip2(l)*hlp1*hlp2
      DsurfD1=DsurfD1 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*Drnag2D1
      DsurfD2=DsurfD2 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*Drnag2D2
     &                + dip2(l)*hlp2*hlp1*hlp11*acut2(l)
      DsurfD3=DsurfD3 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*Drnag2D3

      hlp1 = acut2n(l)*(r(i,2)-rcut2n(l))
      hlp2 = adip2n(l)*(rnag2-rdip2n(l)**2)
      hlp11 = dtanh(hlp1)
      hlp22 = dtanh(hlp2)
      hlp1 = 1.d0/dcosh(hlp1)
      hlp2 = 1.d0/dcosh(hlp2)
      surf = surf - dip2n(l)*hlp1*hlp2
      DsurfD1=DsurfD1 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*Drnag2D1
      DsurfD2=DsurfD2 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*Drnag2D2
     &                + dip2n(l)*hlp2*hlp1*hlp11*acut2n(l)
      DsurfD3=DsurfD3 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*Drnag2D3

      hlp = 1.d0/(
     &              (1.d0-(csb-bc(l))/(bc(l-1)-bc(l)))
     &            * (1.d0-(csb-bc(l))/(bc(l+1)-bc(l)))
     &                                                 )
      hlp1 = (csb-bc(l))**2*hlp
      cog = dexp(-ac(l)*hlp1)
      Dcog = -cog*ac(l)*(2.d0*(csb-bc(l))*hlp
     & +hlp1*(1.d0/(bc(l-1)-csb)+1.d0/(bc(l+1)-csb)))
      DcogD1 = Dcog*DcsbD1
      DcogD2 = Dcog*DcsbD2
      DcogD3 = Dcog*DcsbD3

      e(i) = e(i) + surf*cog
      De(1,i) = De(1,i) + DsurfD1*cog+surf*DcogD1
      De(2,i) = De(2,i) + DsurfD2*cog+surf*DcogD2
      De(3,i) = De(3,i) + DsurfD3*cog+surf*DcogD3
      anorm = anorm + cog
      DanormD1 = DanormD1 + DcogD1
      DanormD2 = DanormD2 + DcogD2
      DanormD3 = DanormD3 + DcogD3
   20 continue
      hlp = 1.d0/anorm
      e(i) = e(i)*hlp
      De(1,i) = (De(1,i) - e(i)*DanormD1)*hlp
      De(2,i) = (De(2,i) - e(i)*DanormD2)*hlp
      De(3,i) = (De(3,i) - e(i)*DanormD3)*hlp

      do j=1,3
      rdum(j) = r(i,j)
      enddo

      call lepscfu22(3,rdum,g_r,3,cfu22,Dcfu22,3)

      e(i) = e(i) + cfu22
      De(1,i) = De(1,i) + Dcfu22(1)
      De(2,i) = De(2,i) + Dcfu22(2)
      De(3,i) = De(3,i) + Dcfu22(3)

      e(i) = e(i)*cevau
      De(1,i) = De(1,i)*cevau
      De(2,i) = De(2,i)*cevau
      De(3,i) = De(3,i)*cevau

10    continue

      return

  100 format(/26('='),' NaFH - U22 diabat ',26('='),
     &      /2x,'mF/(mF+mH) = ',f13.5,2x,'mH/(mF+mH) = ',f13.5)
  101 format(2x,'DeNaF = ',9x,e14.6/2x,'reNaF = ',9x,e14.6,
     & /2x,'cNaF = ',10x,3e14.6/19x,3e14.6)
  103 format(2x,'DeHF = ',10x,e14.6/2x,'reHF = ',10x,e14.6
     & /2x,'cHF = ',11x,3e14.6/19x,3e14.6)
  105 format(2x,'una2p,asmall = ',2x,2e14.6)
  106 format(2x,'y(1,1,i) = ',2x,4e13.5)
  107 format(2x,'y(2,1,i) = ',2x,4e13.5)
  108 format(2x,'y(3,1,i) = ',2x,4e13.5)
  109 format(2x,'y(4,1,i) = ',2x,4e13.5)
  110 format(2x,'y(1,2,i) = ',2x,4e13.5)
  111 format(2x,'y(2,2,i) = ',2x,4e13.5)
  112 format(2x,'y(3,2,i) = ',2x,4e13.5)
  113 format(2x,'y(4,2,i) = ',2x,4e13.5)
  114 format(2x,'y(1,3,i) = ',2x,4e13.5)
  115 format(2x,'y(2,3,i) = ',2x,4e13.5)
  116 format(2x,'y(3,3,i) = ',2x,4e13.5)
  117 format(2x,'y(4,3,i) = ',2x,4e13.5)
  118 format(2x,'c2a = ',7x,4e13.5)
  119 format(2x,'c2b = ',7x,4e13.5)
  120 format(2x,'c2c = ',7x,4e13.5)
  121 format(2x,'r20i = ',6x,4e13.5)
  122 format(2x,'r20a = ',6x,4e13.5)
  123 format(2x,'al20i = ',5x,4e13.5)
  124 format(2x,'al20i1 = ',4x,4e13.5)
  125 format(2x,'al20a = ',5x,4e13.5)
  126 format(2x,'al2a = ',6x,4e13.5)
  127 format(2x,'rr0a = ',6x,4e13.5)
  128 format(2x,'psi = ',7x,4e13.5)
  129 format(2x,'rsi = ',7x,4e13.5)
  130 format(2x,'dip2 = ',6x,4e13.5)
  131 format(2x,'adip2 = ',5x,4e13.5)
  132 format(2x,'rdip2 = ',5x,4e13.5)
  133 format(2x,'acut2 = ',5x,4e13.5)
  134 format(2x,'rcut2 = ',5x,4e13.5)
  135 format(2x,'ry = ',8x,4e13.5)
  136 format(2x,'aly = ',7x,4e13.5)
  137 format(2x,'dy = ',8x,4e13.5)
  138 format(2x,'dip2n = ',5x,4e13.5)
  139 format(2x,'adip2n = ',4x,4e13.5)
  140 format(2x,'rdip2n = ',4x,4e13.5)
  141 format(2x,'acut2n = ',4x,4e13.5)
  142 format(2x,'rcut2n = ',4x,4e13.5
     &    /71('=')/)

      end
c **********************************************************************

      BLOCK DATA GENERAL
      double precision bc
      double precision CNF,CNF1,CHF
      double precision amh,amf,amna,amu
      double precision DeNF,reNF,DeHF,reHF
      double precision ads1,ars1,abs1
      double precision una2p,asmall,cevau

      common/onetwo/CNF(6),CNF1(5),CHF(6),bc(0:5),
     & ads1,ars1,abs1,amh,amf,amna,amu,una2p,asmall,cevau,
     & DeNF,reNF,DeHF,reHF
      data amh/1.007825d0/,amf/18.998403d0/,amna/22.989767d0/,
     &     amu/1822.8885d0/
      data DeNF/4.3830525d0/,reNF/3.6395d0/,
     & CNF/0.168364d0,2.29333d0,-0.848620d-1,-0.109707d-1,
     &   0.246582d-2,0.377012d0/
      data cNF1/0.620115d0,0.694459d-1,0.757570d0,0.283055d0,2.69998d0/
      data DeHF/5.77096d0/,reHF/1.7328d0/,cHF/0.441258d0,3.88056d0,
     & -1.27110d0,-1.37766d0,0.168186d0,2.07230d0/
      data bc/-1.d0,0.d0,1.d0,2.d0,3.0,5.d0/
      data ads1/2.91087d-3/,ars1/7.37707d0/,abs1/7.13223d-1/
      data una2p/2.0973375d0/
      data asmall/0.09d0/
      data cevau/0.036749308867692d0/

      END

      BLOCK DATA UONEONE
      double precision ac
      double precision y,dip1,adip1,rdip1,acut1,rcut1,
     &            c2a,c2b,c2c,dipn,al2,aln0,
     &            dipnn,al2n,aln0n1,aln0n2,rn0n
      common/uone/y(4,3,4),c2a(4),c2b(4),c2c(4),dip1(4),
     &            adip1(4),rdip1(4),acut1(4),rcut1(4),
     &            dipn(4),al2(4),aln0(4),
     &            dipnn(4),al2n(4),aln0n1(4),aln0n2(4),rn0n(4)
      common/angular1/ac(4)
      data y(1,1,1)/0.00000d0/,y(2,1,1)/0.00000d0/,
     &     y(1,1,2)/0.00000d0/,y(2,1,2)/0.00000d0/,
     &     y(1,1,3)/0.00000d0/,y(2,1,3)/0.00000d0/,
     &     y(1,1,4)/0.00000d0/,y(2,1,4)/0.00000d0/,
     &     y(3,1,1)/6.6796d0/,y(4,1,1)/2.7839d0/,
     &     y(3,1,2)/6.6796d0/,y(4,1,2)/2.7839d0/,
     &     y(3,1,3)/7.1360d0/,y(4,1,3)/1.0376d0/,
     &     y(3,1,4)/8.8737d0/,y(4,1,4)/1.9627d0/
      data y(1,2,1)/6.0459d+1/,y(2,2,1)/3.7503d0/,
     &     y(1,2,2)/6.0459d+1/,y(2,2,2)/3.7503d0/,
     &     y(1,2,3)/6.0459d+1/,y(2,2,3)/3.7503d0/,
     &     y(1,2,4)/8.5430d0/,y(2,2,4)/9.00856d-1/,
     &     y(3,2,1)/8.4874d0/,y(4,2,1)/1.0916d+1/,
     &     y(3,2,2)/8.4874d0/,y(4,2,2)/9.9850d0/,
     &     y(3,2,3)/8.4874d0/,y(4,2,3)/5.5898d0/,
     &     y(3,2,4)/9.5899d0/,y(4,2,4)/1.1615d+1/
      data y(1,3,1)/7.3251d0/,y(2,3,1)/2.8774d0/,
     &     y(1,3,2)/7.3251d0/,y(2,3,2)/3.6094d0/,
     &     y(1,3,3)/7.3251d0/,y(2,3,3)/3.9288d0/,
     &     y(1,3,4)/7.5811d0/,y(2,3,4)/4.4272d0/,
     &     y(3,3,1)/1.9300d+3/,y(4,3,1)/2.7091d0/,
     &     y(3,3,2)/3.8969d+3/,y(4,3,2)/2.7591d0/,
     &     y(3,3,3)/5.9130d+3/,y(4,3,3)/2.6881d0/,
     &     y(3,3,4)/5.9007d+3/,y(4,3,4)/3.8929d0/
      data c2a(1)/8.5742d0/,c2b(1)/1.5953d0/,c2c(1)/3.3745d-1/
      data c2a(2)/9.8509d0/,c2b(2)/6.6359d0/,c2c(2)/3.0396d-1/
      data c2a(3)/9.5937d0/,c2b(3)/1.3441d+1/,c2c(3)/2.8812d-1/
      data c2a(4)/1.1386d0/,c2b(4)/8.1112d0/,c2c(4)/2.6061d0/
      data dip1(1)/1.9651d0/,rdip1(1)/8.3109d-1/,
     &     dip1(2)/1.9651d0/,rdip1(2)/8.3109d-1/,
     &     dip1(3)/2.1803d0/,rdip1(3)/7.1647d-1/,
     &     dip1(4)/2.2483d0/,rdip1(4)/9.7128d-1/,
     &     adip1(1)/4.8674d-2/,adip1(2)/5.0603d-2/,
     &     adip1(3)/4.8499d-2/,adip1(4)/8.5924d-2/
      data acut1(1)/6.3761d0/,rcut1(1)/5.2007d0/,
     &     acut1(2)/6.3761d0/,rcut1(2)/3.6585d0/,
     &     acut1(3)/6.0812d0/,rcut1(3)/3.4684d0/,
     &     acut1(4)/5.0516d0/,rcut1(4)/1.6434d0/
      data dipn(1)/1.5754d-1/,al2(1)/8.4809d-2/,
     &     dipn(2)/1.2274d-1/,al2(2)/8.4809d-2/,
     &     dipn(3)/4.1823d-2/,al2(3)/9.4867d-2/,
     &     dipn(4)/3.4350d-1/,al2(4)/6.5542d-2/,
     &      aln0(1)/3.1971d-2/,
     &      aln0(2)/3.5817d-2/,
     &      aln0(3)/1.6686d-1/,
     &      aln0(4)/5.9961d-2/
      data dipnn(1)/0.0000d0/,al2n(1)/0.0000d-2/,
     &     dipnn(2)/3.7501d-3/,al2n(2)/9.2376d0/,
     &     dipnn(3)/7.6323d-3/,al2n(3)/9.3515d-2/,
     &     dipnn(4)/0.0000d0/,al2n(4)/.0000d-2/,
     &      rn0n(1)/0.0d0/,
     &      rn0n(2)/6.6001d0/,
     &      rn0n(3)/5.9949d0/,
     &      rn0n(4)/0.0d0/,
     &      aln0n1(1)/0.0000d-2/,aln0n2(1)/0.0000d-2/,
     &      aln0n1(2)/3.2545d-2/,aln0n2(2)/4.9001d-1/,
     &      aln0n1(3)/8.5914d-2/,aln0n2(3)/2.1382d0/,
     &      aln0n1(4)/0.0000d-2/,aln0n2(4)/0.0000d-2/
      data ac/4.6913d0,5.4684d-1,7.5932d-1,8.0000d-1/

      END

      BLOCK DATA UTWOTWO
      double precision ac
      double precision y,r20i,r20a,al20i,al20a,al2a,
     &            rr0a,ry,aly,dy,c2a,c2b,c2c,al20i1,
     &            dip2,adip2,rdip2,acut2,rcut2,psi,rsi,
     &            dip2n,adip2n,rdip2n,acut2n,rcut2n,
     &            dy1,dy2
      common/utwo/y(4,3,4),c2a(4),c2b(4),c2c(4),r20i(4),
     &      r20a,al20i(4),al20i1(4),al20a,psi,rsi,al2a(4),
     &      rr0a(4),dip2(4),adip2(4),rdip2(4),acut2(4),
     &      rcut2(4),ry(4),aly(4),dy(4),dy1(4),dy2(4),
     &      dip2n(4),adip2n(4),rdip2n(4),acut2n(4),rcut2n(4)
      common/angular2/ac(4)

      data y(1,1,1)/1.0407d+3/,y(2,1,1)/4.5911d0/,
     &     y(1,1,2)/1.1608d+3/,y(2,1,2)/3.8886d0/,
     &     y(1,1,3)/4.9492d+3/,y(2,1,3)/2.6147d0/,
     &     y(1,1,4)/8.3093d+3/,y(2,1,4)/3.2431d0/,
     &     y(3,1,1)/1.4941d0/,y(4,1,1)/4.5525d0/,
     &     y(3,1,2)/1.4941d0/,y(4,1,2)/4.5525d0/,
     &     y(3,1,3)/1.4941d0/,y(4,1,3)/4.5525d0/,
     &     y(3,1,4)/2.8361d0/,y(4,1,4)/4.4686d0/
      data y(1,2,1)/0.0000d0/,y(2,2,1)/0.0000d0/,
     &     y(1,2,2)/0.0000d0/,y(2,2,2)/0.0000d0/,
     &     y(1,2,3)/0.0000d0/,y(2,2,3)/0.0000d0/,
     &     y(1,2,4)/0.0000d0/,y(2,2,4)/0.0000d0/,
     &     y(3,2,1)/1.8076d+3/,y(4,2,1)/2.8815d0/
     &     y(3,2,2)/1.8076d+3/,y(4,2,2)/2.8815d0/,
     &     y(3,2,3)/1.8076d+3/,y(4,2,3)/2.8815d0/,
     &     y(3,2,4)/2.3109d+3/,y(4,2,4)/3.1315d0/
      data y(1,3,1)/5.7330d0/,y(2,3,1)/3.2053d-1/,
     &     y(1,3,2)/5.7330d0/,y(2,3,2)/3.2053d-1/,
     &     y(1,3,3)/5.7330d0/,y(2,3,3)/3.2053d-1/,
     &     y(1,3,4)/5.0393d0/,y(2,3,4)/2.7567d-1/,
     &     y(3,3,1)/2.5774d+5/,y(4,3,1)/4.2267d0/,
     &     y(3,3,2)/2.5774d+5/,y(4,3,2)/4.2267d0/,
     &     y(3,3,3)/2.5774d+5/,y(4,3,3)/4.2267d0/,
     &     y(3,3,4)/1.6506d+5/,y(4,3,4)/3.9246d0/

      data c2a(1)/5.6406d0/,c2b(1)/7.7152d-2/
      data c2a(2)/8.2390d0/,c2b(2)/1.7133d-1/
      data c2a(3)/7.1504d0/,c2b(3)/1.1011d-1/
      data c2a(4)/8.0196d0/,c2b(4)/1.0017d-1/,
     &                         c2c(1)/8.4595d-2/,
     &                         c2c(2)/8.4595d-2/,
     &                         c2c(3)/8.4595d-2/,
     &                         c2c(4)/8.4411d-2/
      data r20i(1)/2.9013d0/,al20i(1)/1.7119d0/
      data r20i(2)/3.1196d0/,al20i(2)/9.4849d-1/
      data r20i(3)/3.0995d0/,al20i(3)/0.0d0/
      data r20i(4)/2.6572d0/,al20i(4)/2.3214d0/,
     &     r20a/3.1305d0/,al20i1(1)/3.5428d0/,
     &     al20a/3.024d+1/,al20i1(2)/6.7257d0/,
     &     psi/0.7d0/,al20i1(3)/3.0240d+1/,
     &     rsi/22.2d0/,al20i1(4)/9.6978d0/,
     &     al2a(1)/5.8088d-1/,rr0a(1)/9.5543d0/,
     &     al2a(2)/5.8088d-1/,rr0a(2)/9.5543d0/,
     &     al2a(3)/5.8088d-1/,rr0a(3)/9.5543d0/,
     &     al2a(4)/8.1525d-1/,rr0a(4)/8.2500d0/
      data dip2(1)/2.6475d-1/
      data dip2(2)/2.6475d-1/
      data dip2(3)/2.6475d-1/
      data dip2(4)/5.0336d-1/
      data adip2(1)/3.2317d-1/,rdip2(1)/3.7260d0/
      data adip2(2)/1.8513d-1/,rdip2(2)/4.5515d0/
      data adip2(3)/4.3999d-2/,rdip2(3)/2.1191d-1/
      data adip2(4)/2.3003d-1/,rdip2(4)/5.0160d0/
      data acut2(1)/2.0799d0/,rcut2(1)/1.2496d0/
      data acut2(2)/1.4204d0/,rcut2(2)/1.0817d0/
      data acut2(3)/7.9086d-1/,rcut2(3)/1.0575d0/
      data acut2(4)/2.8447d0/,rcut2(4)/3.0015d0/
      data dip2n(1)/2.2091d-1/,adip2n(1)/1.0288d-1/
      data dip2n(2)/2.2091d-1/,adip2n(2)/1.0288d-1/
      data dip2n(3)/2.2091d-1/,adip2n(3)/1.0288d-1/
      data dip2n(4)/4.5301d-1/,adip2n(4)/1.1009d-1/
      data rdip2n(1)/5.1554d0/
      data rdip2n(2)/5.1554d0/
      data rdip2n(3)/5.1554d0/
      data rdip2n(4)/6.7003d0/
      data acut2n(1)/3.5607d0/,rcut2n(1)/3.0661d0/
      data acut2n(2)/4.3956d0/,rcut2n(2)/3.0661d0/
      data acut2n(3)/4.1924d0/,rcut2n(3)/3.0661d0/
      data acut2n(4)/3.0765d0/,rcut2n(4)/2.8845d0/
      data ry(1)/1.4328d+1/,aly(1)/4.9044d-1/
      data ry(2)/1.2504d+1/,aly(2)/6.4718d-1/
      data ry(3)/1.2102d+1/,aly(3)/2.5833d-1/
      data ry(4)/7.9381d0/,aly(4)/6.4853d-1/
      data dy(1)/7.5126d-1/,dy1(1)/0.d0/,dy2(1)/0.d0/
      data dy(2)/2.0210d-1/,dy1(2)/0.d0/,dy2(2)/0.d0/
      data dy(3)/1.1749d-1/,dy1(3)/0.d0/,dy2(3)/0.d0/
      data dy(4)/8.1959d-2/,
     &               dy1(4)/-5.03d+2/,dy2(4)/-0.25d0/
      data ac/2.0018d0,4.1908d-1,5.4988d0,3.0178d-1/

      END

      BLOCK DATA UONETWO
      double precision ahf,anaf,anah,r21,r23,r13
      double precision ccNF,ccNHa,ccNHi,ccNHi1,xpNH,xpNH1
      common/ucoupl/ahf(3),anaf(3),anah(3),r21(3),r23(3),r13(3)
      common/ucoupl1/ccNF(4),ccNHa(4),ccNHi(4),ccNHi1(4)
      common/u12nah/xpNH,xpNH1
      data ccNF/3.9960D0,5.4730D-01,-1.1428D0,2.0440D-01/
      data ccNHa/1.00d0,0.80d0,-2.67d0,0.456d0/
      data ccNHi/1.0897d0,-6.5205d-1,7.7968d-1,-1.2933d-2/
      data ccNHi1/-6.1073d-1,-1.7425d-1,-8.3877d-1,
     & -3.5361d-1/
      data xpNH/5.6906d-2/,xpNH1/2.3185d-1/

      data ahf/1.7177d0,2.2707d0,1.2248d-1/,
     &    anaf/6.7786d-1,9.8809d-1,-1.0541d-1/,
     &    anah/6.5680d-1,1.0141d0,-5.1056d-2/
      data r21/3.5114d0,1.3540d0,3.2267d-1/,
     &     r23/-3.3006d-1,2.0804d-1,2.4499d-1/,
     &     r13/1.2431d0,9.9188d-1,1.1051d-1/
      END

c **********************************************************************
C                           DISCLAIMER
C
C   This file was generated on 07/12/00 by the version of
C   ADIFOR compiled on June, 1998.
C
C   ADIFOR was prepared as an account of work sponsored by an
C   agency of the United States Government, Rice University, and
C   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
C   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
C   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
C   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
C   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
C   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
      subroutine lepscfu11(g_p_, r, g_r, ldg_r, cfu11, g_cfu11, ldg_cf
     *u11)
C
C     correction function is the difference of two LEPS surfaces
C     each triplet is weighted by a-f, which was optimezed with a GA
C
        implicit none
C
        double precision a, b, c, d, e, f
C
        double precision dehf, j1_ap, rehf, shf_ac, yhf, j1_ac, snah, ch
     *i, tnah, q1_ap, q1_ac, j2_ap, j2_ac, w, leps_ac, leps_ap, q2_ap, c
     *7, c6, shf_ap, q2_ac, thf, r2, r1, r3, snaf, ronaf, j3_ac, q3_ap, 
     *j3_ap, denah, renah, q3_ac, t2, t1, t3, tnaf, t4, y, thf_ap, thf_a
     *c, tnaf_ap, tnaf_ac, tnah_ap, tnah_ac, cfu11, de, re, c5, beta, c1
     *, c2, c3, c4, r
        dimension r(1, 3)
C
        integer g_pmax_,ia,ib,ic
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_cfu11
        double precision g_r(ldg_r, 1, 3)
        double precision d4_p, d3_p, d17_b, d5_b, d6_b, d7_b, d8_b, d12_
     *v, d2_v, d11_b
        double precision d10_b, d1_w, d2_w, d3_v, d4_v, d5_v, d1_p, d7_v
     *, d8_v, d2_p
        double precision d10_v, d14_v, d4_b, g_r1(g_pmax_), g_r2(g_pmax_
     *), g_r3(g_pmax_), g_y(g_pmax_), g_d1_w(g_pmax_), g_d2_w(g_pmax_), 
     *g_snaf(g_pmax_)
        double precision g_tnaf_ac(g_pmax_), g_tnaf_ap(g_pmax_), g_q3_ac
     *(g_pmax_), g_q3_ap(g_pmax_), g_j3_ac(g_pmax_), g_j3_ap(g_pmax_), g
     *_beta(g_pmax_), g_chi(g_pmax_), g_snah(g_pmax_), g_tnah_ac(g_pmax_
     *)
        double precision g_tnah_ap(g_pmax_), g_q1_ac(g_pmax_), g_q1_ap(g
     *_pmax_), g_j1_ac(g_pmax_), g_j1_ap(g_pmax_), g_yhf(g_pmax_), g_shf
     *_ac(g_pmax_), g_shf_ap(g_pmax_), g_thf_ac(g_pmax_), g_thf_ap(g_pma
     *x_)
        double precision g_q2_ac(g_pmax_), g_q2_ap(g_pmax_), g_j2_ac(g_p
     *max_), g_j2_ap(g_pmax_), g_w(g_pmax_), g_leps_ac(g_pmax_), g_leps_
     *ap(g_pmax_), g_cfu11(ldg_cfu11)
        integer g_ehfid
        save g_j2_ap, g_w, g_leps_ac, g_leps_ap
        save g_j1_ac, g_j1_ap, g_yhf, g_shf_ac, g_shf_ap, g_thf_ac, g_th
     *f_ap, g_q2_ac, g_q2_ap, g_j2_ac
        save g_q3_ap, g_j3_ac, g_j3_ap, g_beta, g_chi, g_snah, g_tnah_ac
     *, g_tnah_ap, g_q1_ac, g_q1_ap
        save g_r1, g_r2, g_r3, g_y, g_d1_w, g_d2_w, g_snaf, g_tnaf_ac, g
     *_tnaf_ap, g_q3_ac
        intrinsic dble
C
C

        do ia = 1,3
        do ib = 1,1
        do ic = 1,3
        if(ia.eq.ic) g_r(ia,ib,ic)=1.d0
        if(ia.ne.ic) g_r(ia,ib,ic)=0.d0
        enddo
        enddo
        enddo

        do g_i_ = 1, g_p_
          g_r1(g_i_) = dble(g_r(g_i_, 1, 1))
        enddo
        r1 = dble(r(1, 1))
C--------
        do g_i_ = 1, g_p_
          g_r2(g_i_) = dble(g_r(g_i_, 1, 2))
        enddo
        r2 = dble(r(1, 2))
C--------
        do g_i_ = 1, g_p_
          g_r3(g_i_) = dble(g_r(g_i_, 1, 3))
        enddo
        r3 = dble(r(1, 3))
C--------
C
        a = 3.58749d0
        b = 2.91300d0
        c = 2.06256d0
        d = 3.58749d0
        e = 3.11828d0
        f = 2.05279d0
C
C     NaF stuff
C     singlets:  approx and acc the same (from paper)
C
        c1 = 0.620115d0
        c2 = 0.694459d-1
        c3 = 0.757570d0
        c4 = 0.283055d0
        c5 = 2.69998d0
        ronaf = 3.6395d0
C
        do g_i_ = 1, g_p_
          g_y(g_i_) = g_r3(g_i_)
        enddo
        y = r3 - ronaf
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-c3) * g_y(g_i_)
        enddo
        d1_w = (-c3) * y
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = (-c5) * g_y(g_i_)
        enddo
        d2_w = (-c5) * y
        d3_v = c1 + c2 * y
        d5_v = exp(d1_w)
        d2_p =  d5_v
        d8_v = exp(d2_w)
        d1_p =  d8_v
        d5_b = c4 * d1_p
        d8_b = d3_v * d2_p
        d10_b = d5_v * c2
        do g_i_ = 1, g_p_
          g_snaf(g_i_) = d5_b * g_d2_w(g_i_) + d8_b * g_d1_w(g_i_) + d10
     *_b * g_y(g_i_)
        enddo
        snaf = d3_v * d5_v + c4 * d8_v
C--------
C
C     triplets:  approx and acc the same
C     used 180 degree triplet from the paper
C
        t1 = 7.3251d0
        t2 = 2.8774d0
        t3 = 1930.0d0
        t4 = 2.7091d0
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-t2) * g_r3(g_i_)
        enddo
        d1_w = (-t2) * r3
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = (-t4) * g_r3(g_i_)
        enddo
        d2_w = (-t4) * r3
        d3_v = exp(d1_w)
        d2_p =  d3_v
        d7_v = exp(d2_w)
        d1_p =  d7_v
        d6_b = c * t3 * d1_p
        d10_b = c * t1 * d2_p
        do g_i_ = 1, g_p_
          g_tnaf_ac(g_i_) = d6_b * g_d2_w(g_i_) + d10_b * g_d1_w(g_i_) +
     * c * g_snaf(g_i_)
        enddo
        tnaf_ac = (snaf + t1 * d3_v + t3 * d7_v) * c
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-t2) * g_r3(g_i_)
        enddo
        d1_w = (-t2) * r3
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = (-t4) * g_r3(g_i_)
        enddo
        d2_w = (-t4) * r3
        d3_v = exp(d1_w)
        d2_p =  d3_v
        d7_v = exp(d2_w)
        d1_p =  d7_v
        d6_b = f * t3 * d1_p
        d10_b = f * t1 * d2_p
        do g_i_ = 1, g_p_
          g_tnaf_ap(g_i_) = d6_b * g_d2_w(g_i_) + d10_b * g_d1_w(g_i_) +
     * f * g_snaf(g_i_)
        enddo
        tnaf_ap = (snaf + t1 * d3_v + t3 * d7_v) * f
C--------
C
        do g_i_ = 1, g_p_
          g_q3_ac(g_i_) = 0.5d0 * g_tnaf_ac(g_i_) + 0.5d0 * g_snaf(g_i_)
        enddo
        q3_ac = 0.5d0 * (snaf + tnaf_ac)
C--------
        do g_i_ = 1, g_p_
          g_q3_ap(g_i_) = 0.5d0 * g_tnaf_ap(g_i_) + 0.5d0 * g_snaf(g_i_)
        enddo
        q3_ap = 0.5d0 * (snaf + tnaf_ap)
C--------
C
        do g_i_ = 1, g_p_
          g_j3_ac(g_i_) = (-0.5d0) * g_tnaf_ac(g_i_) + 0.5d0 * g_snaf(g_
     *i_)
        enddo
        j3_ac = 0.5d0 * (snaf - tnaf_ac)
C--------
        do g_i_ = 1, g_p_
          g_j3_ap(g_i_) = (-0.5d0) * g_tnaf_ap(g_i_) + 0.5d0 * g_snaf(g_
     *i_)
        enddo
        j3_ap = 0.5d0 * (snaf - tnaf_ap)
C--------
C
C     NaH asymptote
C     singlets:  approx and acc the same (from paper)
C
        renah = 7.37707d0
        denah = 0.00291087d0
        do g_i_ = 1, g_p_
          g_beta(g_i_) = 0.0d0
        enddo
        beta = 0.713223d0
C--------
        d4_v = r1 - renah
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-beta) * g_r1(g_i_) + (-d4_v) * g_beta(g_i_)
        enddo
        d1_w = (-beta) * d4_v
        d2_v = exp(d1_w)
        d1_p =  d2_v
        do g_i_ = 1, g_p_
          g_chi(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        chi = d2_v
C--------
        d2_v = denah * chi
        d3_v = chi - 2.d0
        d4_b = d2_v + d3_v * denah
        do g_i_ = 1, g_p_
          g_snah(g_i_) = d4_b * g_chi(g_i_)
        enddo
        snah = d2_v * d3_v
C--------
C
C     triplets:  same for approx and acc
C     used the 180 degree triplet
C
        t1 = 6.6796d0
        t2 = 2.7839d0
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-t2) * g_r1(g_i_)
        enddo
        d1_w = (-t2) * r1
        d3_v = exp(d1_w)
        d1_p =  d3_v
        d6_b = a * t1 * d1_p
        do g_i_ = 1, g_p_
          g_tnah_ac(g_i_) = d6_b * g_d1_w(g_i_) + a * g_snah(g_i_)
        enddo
        tnah_ac = (snah + t1 * d3_v) * a
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-t2) * g_r1(g_i_)
        enddo
        d1_w = (-t2) * r1
        d3_v = exp(d1_w)
        d1_p =  d3_v
        d6_b = d * t1 * d1_p
        do g_i_ = 1, g_p_
          g_tnah_ap(g_i_) = d6_b * g_d1_w(g_i_) + d * g_snah(g_i_)
        enddo
        tnah_ap = (snah + t1 * d3_v) * d
C--------
C
        do g_i_ = 1, g_p_
          g_q1_ac(g_i_) = 0.5d0 * g_tnah_ac(g_i_) + 0.5d0 * g_snah(g_i_)
        enddo
        q1_ac = 0.5d0 * (snah + tnah_ac)
C--------
        do g_i_ = 1, g_p_
          g_q1_ap(g_i_) = 0.5d0 * g_tnah_ap(g_i_) + 0.5d0 * g_snah(g_i_)
        enddo
        q1_ap = 0.5d0 * (snah + tnah_ap)
C--------
C
        do g_i_ = 1, g_p_
          g_j1_ac(g_i_) = (-0.5d0) * g_tnah_ac(g_i_) + 0.5d0 * g_snah(g_
     *i_)
        enddo
        j1_ac = 0.5d0 * (snah - tnah_ac)
C--------
        do g_i_ = 1, g_p_
          g_j1_ap(g_i_) = (-0.5d0) * g_tnah_ap(g_i_) + 0.5d0 * g_snah(g_
     *i_)
        enddo
        j1_ap = 0.5d0 * (snah - tnah_ap)
C--------
C
C     HF asymptote
C     accurate singlet
C
        dehf = 6.122d0
        rehf = 1.733d0
        c1 = 1.1622d0
        c2 = -0.025647d0
        c3 = 0.059062d0
        c4 = 2.1042d0
C
        do g_i_ = 1, g_p_
          g_yhf(g_i_) = g_r2(g_i_)
        enddo
        yhf = r2 - c4
C--------
        d4_v = yhf * yhf
        d1_p = 2.0d0 * yhf
        d5_b = c3 * d1_p + c2
        do g_i_ = 1, g_p_
          g_beta(g_i_) = d5_b * g_yhf(g_i_)
        enddo
        beta = c1 + c2 * yhf + c3 * d4_v
C--------
        d4_v = r2 - rehf
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
        enddo
        d1_w = (-beta) * d4_v
        d2_v = exp(d1_w)
        d2_p =  d2_v
        d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
        d1_p = 2.0d0 * (1.d0 - d2_v)
        d6_b = (-(dehf * d1_p)) * d2_p
        do g_i_ = 1, g_p_
          g_shf_ac(g_i_) = d6_b * g_d1_w(g_i_)
        enddo
        shf_ac = dehf * d4_v - dehf
C--------
C
C     approx singlet (from paper)
C
        dehf = 5.77096d0
        c1 = 0.441258d0
        c2 = 3.88056d0
        c3 = -1.441258d0
        c4 = -1.2711d0
        c5 = -1.37766d0
        c6 = 0.168186d0
        c7 = 2.07230d0
        rehf = 1.7328d0
C
        do g_i_ = 1, g_p_
          g_yhf(g_i_) = g_r2(g_i_)
        enddo
        yhf = r2 - rehf
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-c2) * g_yhf(g_i_)
        enddo
        d1_w = (-c2) * yhf
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = (-c7) * g_yhf(g_i_)
        enddo
        d2_w = (-c7) * yhf
        d2_v = exp(d1_w)
        d4_p =  d2_v
        d7_v = yhf * yhf
        d3_p = 2.0d0 * yhf
        d10_v = yhf ** ( 3 - 2)
        d10_v =  d10_v * yhf
        d2_p =  3 *  d10_v
        d10_v =  d10_v * yhf
        d12_v = c3 + c4 * yhf + c5 * d7_v + c6 * d10_v
        d14_v = exp(d2_w)
        d1_p =  d14_v
        d5_b = dehf * d14_v
        d7_b = dehf * d12_v * d1_p
        d11_b = d5_b * c6 * d2_p + d5_b * c5 * d3_p + d5_b * c4
        d17_b = dehf * c1 * d4_p
        do g_i_ = 1, g_p_
          g_shf_ap(g_i_) = d7_b * g_d2_w(g_i_) + d11_b * g_yhf(g_i_) + d
     *17_b * g_d1_w(g_i_)
        enddo
        shf_ap = dehf * (c1 * d2_v + d12_v * d14_v)
C--------
C
C     triplets:  same for acc and approx, taken from constant beta fit to MRDCI 
Cdata
C
        dehf = 5.77096d0
        rehf = 1.7328d0
        do g_i_ = 1, g_p_
          g_beta(g_i_) = 0.0d0
        enddo
        beta = 1.26686217d0
C--------
C
        d4_v = r2 - rehf
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
        enddo
        d1_w = (-beta) * d4_v
        d2_v = exp(d1_w)
        d2_p =  d2_v
        d4_v = (1.d0 + d2_v) * (1.d0 + d2_v)
        d1_p = 2.0d0 * (1.d0 + d2_v)
        d7_b = b * dehf * d1_p * d2_p
        do g_i_ = 1, g_p_
          g_thf_ac(g_i_) = d7_b * g_d1_w(g_i_)
        enddo
        thf_ac = (dehf * d4_v - dehf) * b
C--------
        d4_v = r2 - rehf
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
        enddo
        d1_w = (-beta) * d4_v
        d2_v = exp(d1_w)
        d2_p =  d2_v
        d4_v = (1.d0 + d2_v) * (1.d0 + d2_v)
        d1_p = 2.0d0 * (1.d0 + d2_v)
        d7_b = e * dehf * d1_p * d2_p
        do g_i_ = 1, g_p_
          g_thf_ap(g_i_) = d7_b * g_d1_w(g_i_)
        enddo
        thf_ap = (dehf * d4_v - dehf) * e
C--------
C
        do g_i_ = 1, g_p_
          g_q2_ac(g_i_) = 0.5d0 * g_thf_ac(g_i_) + 0.5d0 * g_shf_ac(g_i_
     *)
        enddo
        q2_ac = 0.5d0 * (shf_ac + thf_ac)
C--------
        do g_i_ = 1, g_p_
          g_q2_ap(g_i_) = 0.5d0 * g_thf_ap(g_i_) + 0.5d0 * g_shf_ap(g_i_
     *)
        enddo
        q2_ap = 0.5d0 * (shf_ap + thf_ap)
C--------
C
        do g_i_ = 1, g_p_
          g_j2_ac(g_i_) = (-0.5d0) * g_thf_ac(g_i_) + 0.5d0 * g_shf_ac(g
     *_i_)
        enddo
        j2_ac = 0.5d0 * (shf_ac - thf_ac)
C--------
        do g_i_ = 1, g_p_
          g_j2_ap(g_i_) = (-0.5d0) * g_thf_ap(g_i_) + 0.5d0 * g_shf_ap(g
     *_i_)
        enddo
        j2_ap = 0.5d0 * (shf_ap - thf_ap)
C--------
C
C     compute LEPS functions
C
        do g_i_ = 1, g_p_
          g_j1_ac(g_i_) = g_j1_ac(g_i_)
        enddo
        j1_ac = j1_ac
C--------
        do g_i_ = 1, g_p_
          g_j2_ac(g_i_) = g_j2_ac(g_i_)
        enddo
        j2_ac = j2_ac
C--------
        do g_i_ = 1, g_p_
          g_j3_ac(g_i_) = g_j3_ac(g_i_)
        enddo
        j3_ac = j3_ac
C--------
        d4_b = -j1_ac + (-j2_ac) + j3_ac + j3_ac
        d8_b = -j3_ac + (-j1_ac) + j2_ac + j2_ac
        d5_b = -j3_ac + (-j2_ac) + j1_ac + j1_ac
        do g_i_ = 1, g_p_
          g_w(g_i_) = d4_b * g_j3_ac(g_i_) + d8_b * g_j2_ac(g_i_) + d5_b
     * * g_j1_ac(g_i_)
        enddo
        w = j1_ac * j1_ac + j2_ac * j2_ac + j3_ac * j3_ac - j1_ac * j2_a
     *c - j2_ac * j3_ac - j3_ac * j1_ac
C--------
C
        d7_v = sqrt(w)

c        if ( w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d7_v)
c        else
c           call ehufDO (9,w, d7_v, d1_p,
c     +g_ehfid,
c     +482)
c        endif
        do g_i_ = 1, g_p_
          g_leps_ac(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ac(g_i_) + g_q2_a
     *c(g_i_) + g_q1_ac(g_i_)
        enddo
        leps_ac = q1_ac + q2_ac + q3_ac - d7_v
C--------
C
        do g_i_ = 1, g_p_
          g_j1_ap(g_i_) = g_j1_ap(g_i_)
        enddo
        j1_ap = j1_ap
C--------
        do g_i_ = 1, g_p_
          g_j2_ap(g_i_) = g_j2_ap(g_i_)
        enddo
        j2_ap = j2_ap
C--------
        do g_i_ = 1, g_p_
          g_j3_ap(g_i_) = g_j3_ap(g_i_)
        enddo
        j3_ap = j3_ap
C--------
        d4_b = -j1_ap + (-j2_ap) + j3_ap + j3_ap
        d8_b = -j3_ap + (-j1_ap) + j2_ap + j2_ap
        d5_b = -j3_ap + (-j2_ap) + j1_ap + j1_ap
        do g_i_ = 1, g_p_
          g_w(g_i_) = d4_b * g_j3_ap(g_i_) + d8_b * g_j2_ap(g_i_) + d5_b
     * * g_j1_ap(g_i_)
        enddo
        w = j1_ap * j1_ap + j2_ap * j2_ap + j3_ap * j3_ap - j1_ap * j2_a
     *p - j2_ap * j3_ap - j3_ap * j1_ap
C--------
C
        d7_v = sqrt(w)

c        if ( w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d7_v)
c        else
c           call ehufDO (9,w, d7_v, d1_p,
c     +g_ehfid,
c     +524)
c        endif
        do g_i_ = 1, g_p_
          g_leps_ap(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ap(g_i_) + g_q2_a
     *p(g_i_) + g_q1_ap(g_i_)
        enddo
        leps_ap = q1_ap + q2_ap + q3_ap - d7_v
C--------
C
        do g_i_ = 1, g_p_
          g_cfu11(g_i_) = -g_leps_ap(g_i_) + g_leps_ac(g_i_)
        enddo
        cfu11 = leps_ac - leps_ap + 6.122d0 - 5.77096d0
C--------
      end

c **********************************************************************
C                           DISCLAIMER
C
C   This file was generated on 08/31/00 by the version of
C   ADIFOR compiled on June, 1998.
C
C   ADIFOR was prepared as an account of work sponsored by an
C   agency of the United States Government, Rice University, and
C   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
C   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
C   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
C   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
C   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
C   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
      subroutine lepscfu22(g_p_, r, g_r, ldg_r, cfu22, g_cfu22, ldg_cf
     *u22)
C
        implicit none
C
        double precision dehf, j1_ap, rehf, shf_ac, yhf, j1_ac, snah, ch
     *i, tnah, q1_ap, q1_ac, j2_ap, j2_ac, w, leps_ac, leps_ap, q2_ap, c
     *7, c6, shf_ap, q2_ac, thf, r2, r1, r3, snaf, ronaf, a, b, c, d, e,
     * f, j3_ac, q3_ap, j3_ap, denah, renah, q3_ac, t2, t1, t3, tnaf, t4
     *, y, thf_ap, thf_ac, tnaf_ap, tnaf_ac, tnah_ap, tnah_ac, r
        double precision shf_u11, nap, t3nah, t4nah, r0, swit, shfa, alp
     *h0, ynaf, snaf_ap, denaf, renaf, t1nah, t2nah, c2naf, snaf_ac
        double precision c1nf, c2nf, c3nf, c4nf, c5nf, c6nf, c7nf, hlp1,
     * hlp2, hlp22, thf_ap_as, t1hf, dt1hf, t2hf, dt2hf, t, b1, b2, e1, 
     *e2, hlp, a3, b3, cc3, d3, e3, f3, cfu22, c5, c4, c3, c2, c1, cc1, 
     *cc2, f1, f2, beta
        dimension r(1, 3)
C
        integer g_pmax_,ia,ib,ic
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_cfu22
        double precision d4_p, d3_p, d17_b, d8_b, d12_v, d2_p, d2_w, d10
     *_v, d2_v, d11_b
        double precision d2_b, d1_w, d3_v, d3_b, d1_p, d4_v, d5_v, d6_v,
     * d7_v, d8_v
        double precision d14_v, d4_b, d5_b, d6_b, d7_b, g_r1(g_pmax_), g
     *_r(ldg_r, 1, 3), g_r2(g_pmax_), g_r3(g_pmax_), g_ynaf(g_pmax_)
        double precision g_d1_w(g_pmax_), g_hlp1(g_pmax_), g_hlp2(g_pmax
     *_), g_hlp22(g_pmax_), g_snaf_ap(g_pmax_), g_beta(g_pmax_), g_snaf_
     *ac(g_pmax_), g_hlp(g_pmax_), g_tnaf_ac(g_pmax_), g_tnaf_ap(g_pmax_
     *)
        double precision g_q3_ac(g_pmax_), g_q3_ap(g_pmax_), g_j3_ac(g_p
     *max_), g_j3_ap(g_pmax_), g_d2_w(g_pmax_), g_tnah_ac(g_pmax_), g_tn
     *ah_ap(g_pmax_), g_q1_ac(g_pmax_), g_q1_ap(g_pmax_), g_j1_ac(g_pmax
     *_)
        double precision g_j1_ap(g_pmax_), g_thf_ac(g_pmax_), g_thf_ap(g
     *_pmax_), g_thf(g_pmax_), g_yhf(g_pmax_), g_shf_u11(g_pmax_), g_shf
     *a(g_pmax_), g_swit(g_pmax_), g_shf_ac(g_pmax_), g_thf_ap_as(g_pmax
     *_)
        double precision g_shf_ap(g_pmax_), g_q2_ac(g_pmax_), g_q2_ap(g_
     *pmax_), g_j2_ac(g_pmax_), g_j2_ap(g_pmax_), g_w(g_pmax_), g_leps_a
     *c(g_pmax_), g_leps_ap(g_pmax_), g_cfu22(ldg_cfu22)
        integer g_ehfid
        save g_leps_ac, g_leps_ap
        save g_shfa, g_swit, g_shf_ac, g_thf_ap_as, g_shf_ap, g_q2_ac, g
     *_q2_ap, g_j2_ac, g_j2_ap, g_w
        save g_tnah_ap, g_q1_ac, g_q1_ap, g_j1_ac, g_j1_ap, g_thf_ac, g_
     *thf_ap, g_thf, g_yhf, g_shf_u11
        save g_snaf_ac, g_hlp, g_tnaf_ac, g_tnaf_ap, g_q3_ac, g_q3_ap, g
     *_j3_ac, g_j3_ap, g_d2_w, g_tnah_ac
        save g_r1, g_r2, g_r3, g_ynaf, g_d1_w, g_hlp1, g_hlp2, g_hlp22, 
     *g_snaf_ap, g_beta
C
C
        do ia = 1,3
        do ib = 1,1
        do ic = 1,3
        if(ia.eq.ic) g_r(ia,ib,ic)=1.d0
        if(ia.ne.ic) g_r(ia,ib,ic)=0.d0
        enddo
        enddo
        enddo
        do g_i_ = 1, g_p_
          g_r1(g_i_) = g_r(g_i_, 1, 1)
        enddo
        r1 = r(1, 1)
C--------
        do g_i_ = 1, g_p_
          g_r2(g_i_) = g_r(g_i_, 1, 2)
        enddo
        r2 = r(1, 2)
C--------
        do g_i_ = 1, g_p_
          g_r3(g_i_) = g_r(g_i_, 1, 3)
        enddo
        r3 = r(1, 3)
C--------
C
C
        b1 = 5.125d0
        b2 = 5.125d0
        cc1 = 1.92d0
        cc2 = 9.01d0
        e1 = 5.54d0
        e2 = 4.15d0
        f1 = 1.62d0
        f2 = 8.80d0
C
C     NaF stuff
C     singlets:  approx (from paper)
C
        c1nf = 0.168364d0
        c2nf = 2.29333d0
        c3nf = -1.168364d0
        c4nf = -0.084862d0
        c5nf = -0.0109707d0
        c6nf = 0.00246582d0
        c7nf = 0.377012d0
        denaf = 4.3830525d0
        renaf = 3.6359d0
C
        do g_i_ = 1, g_p_
          g_ynaf(g_i_) = g_r3(g_i_)
        enddo
        ynaf = r3 - renaf
C--------
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-c2nf) * g_ynaf(g_i_)
        enddo
        d1_w = (-c2nf) * ynaf
        d2_v = exp(d1_w)
        d1_p =  d2_v
        d3_b = c1nf * d1_p
        do g_i_ = 1, g_p_
          g_hlp1(g_i_) = d3_b * g_d1_w(g_i_)
        enddo
        hlp1 = c1nf * d2_v
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-c7nf) * g_ynaf(g_i_)
        enddo
        d1_w = (-c7nf) * ynaf
        d2_v = exp(d1_w)
        d1_p =  d2_v
        do g_i_ = 1, g_p_
          g_hlp2(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        hlp2 = d2_v
C--------
C
        d4_v = ynaf * ynaf
        d2_p = 2.0d0 * ynaf
        d7_v = ynaf ** ( 3 - 2)
        d7_v =  d7_v * ynaf
        d1_p =  3 *  d7_v
        d7_v =  d7_v * ynaf
        d5_b = c6nf * d1_p + c5nf * d2_p + c4nf
        do g_i_ = 1, g_p_
          g_hlp22(g_i_) = d5_b * g_ynaf(g_i_)
        enddo
        hlp22 = c3nf + c4nf * ynaf + c5nf * d4_v + c6nf * d7_v
C--------
C
        d5_b = denaf * hlp2
        d6_b = denaf * hlp22
        do g_i_ = 1, g_p_
          g_snaf_ap(g_i_) = d6_b * g_hlp2(g_i_) + d5_b * g_hlp22(g_i_) +
     * denaf * g_hlp1(g_i_)
        enddo
        snaf_ap = denaf * (hlp1 + hlp22 * hlp2)
C--------
C
C     singlet:  accurate from fit to experimental
C
        c1 = 0.32453d0
        c2 = 1.5102d0
        c3 = 3.0938d0
        c4 = 1.7107d0
C
        denaf = 4.94d0
        renaf = 3.6395d0
C

        if ( (r3 / c3) .ne. 0.0d0 ) then
           d3_v = (r3 / c3) ** ( c4 - 2.0d0)
           d3_v =  d3_v * (r3 / c3)
           d2_p =  c4 *  d3_v
           d3_v =  d3_v * (r3 / c3)
        else
C          ((r3 / c3) = 0)
           d3_v = (r3 / c3) **  c4

              d2_p = 0.0d0
        endif

        if ( (r3 / c3) .ne. 0.0d0 ) then
           d6_v = (r3 / c3) ** ( c4 - 2.0d0)
           d6_v =  d6_v * (r3 / c3)
           d1_p =  c4 *  d6_v
           d6_v =  d6_v * (r3 / c3)
        else
C          ((r3 / c3) = 0)
           d6_v = (r3 / c3) **  c4

              d1_p = 0.0d0
        endif
        d7_v = c1 + d6_v
        d8_v = (c2 + d3_v) / d7_v
        d7_b = c1 * ((-d8_v) / d7_v) * d1_p * (1.0d0 / c3) + c1 * (1.0d0
     * / d7_v) * d2_p * (1.0d0 / c3)
        do g_i_ = 1, g_p_
          g_beta(g_i_) = d7_b * g_r3(g_i_)
        enddo
        beta = c1 * d8_v
C--------
C
        d4_v = r3 - renaf
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-beta) * g_r3(g_i_) + (-d4_v) * g_beta(g_i_)
        enddo
        d1_w = (-beta) * d4_v
        d2_v = exp(d1_w)
        d2_p =  d2_v
        d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
        d1_p = 2.0d0 * (1.d0 - d2_v)
        d6_b = (-(denaf * d1_p)) * d2_p
        do g_i_ = 1, g_p_
          g_snaf_ac(g_i_) = d6_b * g_d1_w(g_i_)
        enddo
        snaf_ac = denaf * d4_v - denaf
C--------
C
C
C     triplets:  approx and acc the same
C     anti-Morse from a constant beta fit to the approx curve
C
        denaf = 4.49d0
        renaf = 3.6395d0
        do g_i_ = 1, g_p_
          g_beta(g_i_) = 0.0d0
        enddo
        beta = 0.69614d0
C--------
        d4_v = r3 - renaf
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-beta) * g_r3(g_i_) + (-d4_v) * g_beta(g_i_)
        enddo
        d1_w = (-beta) * d4_v
        d2_v = exp(d1_w)
        d1_p =  d2_v
        do g_i_ = 1, g_p_
          g_hlp(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        hlp = d2_v
C--------
C
        d2_v = hlp * hlp
        d1_p = 2.0d0 * hlp
        d2_b = 0.5d0 * denaf
        d5_b = d2_b * (cc2 * 2.d0) + d2_b * cc1 * d1_p
        do g_i_ = 1, g_p_
          g_tnaf_ac(g_i_) = d5_b * g_hlp(g_i_)
        enddo
        tnaf_ac = 0.5d0 * denaf * (cc1 * d2_v + cc2 * 2.d0 * hlp)
C--------
        d2_v = hlp * hlp
        d1_p = 2.0d0 * hlp
        d2_b = 0.5d0 * denaf
        d5_b = d2_b * (f2 * 2.d0) + d2_b * f1 * d1_p
        do g_i_ = 1, g_p_
          g_tnaf_ap(g_i_) = d5_b * g_hlp(g_i_)
        enddo
        tnaf_ap = 0.5d0 * denaf * (f1 * d2_v + f2 * 2.d0 * hlp)
C--------
C
        do g_i_ = 1, g_p_
          g_q3_ac(g_i_) = 0.5d0 * g_tnaf_ac(g_i_) + 0.5d0 * g_snaf_ac(g_
     *i_)
        enddo
        q3_ac = 0.5d0 * (snaf_ac + tnaf_ac)
C--------
        do g_i_ = 1, g_p_
          g_q3_ap(g_i_) = 0.5d0 * g_tnaf_ap(g_i_) + 0.5d0 * g_snaf_ap(g_
     *i_)
        enddo
        q3_ap = 0.5d0 * (snaf_ap + tnaf_ap)
C--------
C
        do g_i_ = 1, g_p_
          g_j3_ac(g_i_) = (-0.5d0) * g_tnaf_ac(g_i_) + 0.5d0 * g_snaf_ac
     *(g_i_)
        enddo
        j3_ac = 0.5d0 * (snaf_ac - tnaf_ac)
C--------
        do g_i_ = 1, g_p_
          g_j3_ap(g_i_) = (-0.5d0) * g_tnaf_ap(g_i_) + 0.5d0 * g_snaf_ap
     *(g_i_)
        enddo
        j3_ap = 0.5d0 * (snaf_ap - tnaf_ap)
C--------
C
C
C     NaH asymptote
C     singlets:  approx and acc the same (from paper)
C
        snah = 0.d0
C
C     triplets:  same for approx and acc
C     used the 180 degree triplet
C
        t1nah = 1040.7d0
        t2nah = 4.5911d0
        t3nah = 1.4941d0
        t4nah = 4.5525d0
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-t2nah) * g_r1(g_i_)
        enddo
        d1_w = (-t2nah) * r1
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = (-t4nah) * g_r1(g_i_)
        enddo
        d2_w = (-t4nah) * r1
        d2_v = exp(d1_w)
        d2_p =  d2_v
        d5_v = exp(d2_w)
        d1_p =  d5_v
        d5_b = t3nah * d1_p
        d7_b = t1nah * d2_p
        do g_i_ = 1, g_p_
          g_tnah_ac(g_i_) = d5_b * g_d2_w(g_i_) + d7_b * g_d1_w(g_i_)
        enddo
        tnah_ac = t1nah * d2_v + t3nah * d5_v
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-t2nah) * g_r1(g_i_)
        enddo
        d1_w = (-t2nah) * r1
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = (-t4nah) * g_r1(g_i_)
        enddo
        d2_w = (-t4nah) * r1
        d2_v = exp(d1_w)
        d2_p =  d2_v
        d5_v = exp(d2_w)
        d1_p =  d5_v
        d5_b = t3nah * d1_p
        d7_b = t1nah * d2_p
        do g_i_ = 1, g_p_
          g_tnah_ap(g_i_) = d5_b * g_d2_w(g_i_) + d7_b * g_d1_w(g_i_)
        enddo
        tnah_ap = t1nah * d2_v + t3nah * d5_v
C--------
C
        do g_i_ = 1, g_p_
          g_q1_ac(g_i_) = 0.5d0 * g_tnah_ac(g_i_)
        enddo
        q1_ac = 0.5d0 * (snah + tnah_ac)
C--------
        do g_i_ = 1, g_p_
          g_q1_ap(g_i_) = 0.5d0 * g_tnah_ap(g_i_)
        enddo
        q1_ap = 0.5d0 * (snah + tnah_ap)
C--------
C
        do g_i_ = 1, g_p_
          g_j1_ac(g_i_) = (-0.5d0) * g_tnah_ac(g_i_)
        enddo
        j1_ac = 0.5d0 * (snah - tnah_ac)
C--------
        do g_i_ = 1, g_p_
          g_j1_ap(g_i_) = (-0.5d0) * g_tnah_ap(g_i_)
        enddo
        j1_ap = 0.5d0 * (snah - tnah_ap)
C--------
C
C     HF asymptote
C
C     triplets:  
C     anti-morse from a 
C     constant beta fit to MRDCI data
C     same for approx and accurate
C
C
        dehf = 5.77096d0
        rehf = 1.7328d0
        do g_i_ = 1, g_p_
          g_beta(g_i_) = 0.0d0
        enddo
        beta = 1.2669d0
C--------
        d4_v = r2 - rehf
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
        enddo
        d1_w = (-beta) * d4_v
        d2_v = exp(d1_w)
        d1_p =  d2_v
        do g_i_ = 1, g_p_
          g_hlp(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        hlp = d2_v
C--------
        d3_v = hlp * hlp
        d1_p = 2.0d0 * hlp
        d2_b = 0.5d0 * dehf
        d6_b = d2_b * e2 * d1_p + d2_b * (e1 * 2.d0)
        do g_i_ = 1, g_p_
          g_thf_ac(g_i_) = d6_b * g_hlp(g_i_)
        enddo
        thf_ac = 0.5d0 * dehf * (e1 * 2.d0 * hlp + e2 * d3_v)
C--------
C
        dehf = 5.77096d0
        rehf = 1.7328d0
        do g_i_ = 1, g_p_
          g_beta(g_i_) = 0.0d0
        enddo
        beta = 1.2669d0
C--------
        d4_v = r2 - rehf
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
        enddo
        d1_w = (-beta) * d4_v
        d2_v = exp(d1_w)
        d1_p =  d2_v
        do g_i_ = 1, g_p_
          g_hlp(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        hlp = d2_v
C--------
        d3_v = hlp * hlp
        d1_p = 2.0d0 * hlp
        d2_b = 0.5d0 * dehf
        d6_b = d2_b * b2 * d1_p + d2_b * (b1 * 2.d0)
        do g_i_ = 1, g_p_
          g_thf_ap(g_i_) = d6_b * g_hlp(g_i_)
        enddo
        thf_ap = 0.5d0 * dehf * (b1 * 2.d0 * hlp + b2 * d3_v)
C--------
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-15.d0) * g_r2(g_i_)
        enddo
        d1_w = (-15.d0) * (r2 - 2.d0)
        d2_v = exp(d1_w)
        d1_p =  d2_v
        do g_i_ = 1, g_p_
          g_thf(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        thf = d2_v
C--------
C
C     accurate singlet
C     uses u11 singlet which is
C
        dehf = 6.122d0
        rehf = 1.733d0
        c1 = 1.1622d0
        c2 = -0.025647d0
        c3 = 0.059062d0
        c4 = 2.1042d0
C
        do g_i_ = 1, g_p_
          g_yhf(g_i_) = g_r2(g_i_)
        enddo
        yhf = r2 - c4
C--------
        d4_v = yhf * yhf
        d1_p = 2.0d0 * yhf
        d5_b = c3 * d1_p + c2
        do g_i_ = 1, g_p_
          g_beta(g_i_) = d5_b * g_yhf(g_i_)
        enddo
        beta = c1 + c2 * yhf + c3 * d4_v
C--------
        d4_v = r2 - rehf
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-beta) * g_r2(g_i_) + (-d4_v) * g_beta(g_i_)
        enddo
        d1_w = (-beta) * d4_v
        d2_v = exp(d1_w)
        d2_p =  d2_v
        d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
        d1_p = 2.0d0 * (1.d0 - d2_v)
        d6_b = (-(dehf * d1_p)) * d2_p
        do g_i_ = 1, g_p_
          g_shf_u11(g_i_) = d6_b * g_d1_w(g_i_)
        enddo
        shf_u11 = dehf * d4_v - dehf
C--------
C
        nap = 2.0973375d0
C
        do g_i_ = 1, g_p_
          g_shfa(g_i_) = g_shf_u11(g_i_)
        enddo
        shfa = shf_u11 + nap
C--------
C
        alph0 = 15.d0
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-alph0) * g_thf(g_i_) + alph0 * g_shfa(g_i_)
        enddo
        d1_w = alph0 * (shfa - thf)
        d2_v = tanh (d1_w)
        d1_p = 1.0d0 - ( d2_v *  d2_v)
        do g_i_ = 1, g_p_
          g_swit(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        swit = d2_v
C--------
        d6_v = 0.5d0 * (shfa - thf)
        d6_b = (-swit) * 0.5d0
        d7_b = d6_b + 0.5d0
        d8_b = -d6_b + 0.5d0
        do g_i_ = 1, g_p_
          g_shf_ac(g_i_) = (-d6_v) * g_swit(g_i_) + d8_b * g_thf(g_i_) +
     * d7_b * g_shfa(g_i_)
        enddo
        shf_ac = 0.5d0 * (shfa + thf) - d6_v * swit
C--------
C
C
C     approx singlet (from paper)
C     uses u11 singlet, which is
C
        dehf = 5.77096d0
        c1 = 0.441258d0
        c2 = 3.88056d0
        c3 = -1.441258d0
        c4 = -1.2711d0
        c5 = -1.37766d0
        c6 = 0.168186d0
        c7 = 2.07230d0
        rehf = 1.7328d0
C
        do g_i_ = 1, g_p_
          g_yhf(g_i_) = g_r2(g_i_)
        enddo
        yhf = r2 - rehf
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = (-c2) * g_yhf(g_i_)
        enddo
        d1_w = (-c2) * yhf
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = (-c7) * g_yhf(g_i_)
        enddo
        d2_w = (-c7) * yhf
        d2_v = exp(d1_w)
        d4_p =  d2_v
        d7_v = yhf * yhf
        d3_p = 2.0d0 * yhf
        d10_v = yhf ** ( 3 - 2)
        d10_v =  d10_v * yhf
        d2_p =  3 *  d10_v
        d10_v =  d10_v * yhf
        d12_v = c3 + c4 * yhf + c5 * d7_v + c6 * d10_v
        d14_v = exp(d2_w)
        d1_p =  d14_v
        d5_b = dehf * d14_v
        d7_b = dehf * d12_v * d1_p
        d11_b = d5_b * c6 * d2_p + d5_b * c5 * d3_p + d5_b * c4
        d17_b = dehf * c1 * d4_p
        do g_i_ = 1, g_p_
          g_shf_u11(g_i_) = d7_b * g_d2_w(g_i_) + d11_b * g_yhf(g_i_) + 
     *d17_b * g_d1_w(g_i_)
        enddo
        shf_u11 = dehf * (c1 * d2_v + d12_v * d14_v)
C--------
C
        nap = 2.0973375d0
C
        do g_i_ = 1, g_p_
          g_shfa(g_i_) = g_shf_u11(g_i_)
        enddo
        shfa = shf_u11 + nap
C--------
C
        alph0 = 15.d0
        r0 = 3.1305d0
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = alph0 * g_r2(g_i_)
        enddo
        d1_w = alph0 * (r2 - r0)
        d2_v = tanh (d1_w)
        d1_p = 1.0d0 - ( d2_v *  d2_v)
        d4_b = 0.5d0 * d1_p
        do g_i_ = 1, g_p_
          g_swit(g_i_) = d4_b * g_d1_w(g_i_)
        enddo
        swit = 0.5d0 * (1.d0 + d2_v)
C--------
C
        t1hf = 1807.6d0
        dt1hf = 0.d0
        t2hf = 2.8815d0
        dt2hf = 0.d0
        d2_b = -(t2hf + dt2hf)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d2_b * g_r2(g_i_)
        enddo
        d1_w = (-(t2hf + dt2hf)) * r2
        d2_v = exp(d1_w)
        d1_p =  d2_v
        d3_b = (t1hf + dt1hf) * d1_p
        do g_i_ = 1, g_p_
          g_thf_ap_as(g_i_) = d3_b * g_d1_w(g_i_)
        enddo
        thf_ap_as = (t1hf + dt1hf) * d2_v
C--------
C
        d3_v = thf_ap_as - shfa
        d2_b = 1.0d0 + (-swit)
        do g_i_ = 1, g_p_
          g_shf_ap(g_i_) = d3_v * g_swit(g_i_) + swit * g_thf_ap_as(g_i_
     *) + d2_b * g_shfa(g_i_)
        enddo
        shf_ap = shfa + d3_v * swit
C--------
C
C     compute coulomb and exchange integrals
C
        do g_i_ = 1, g_p_
          g_q2_ac(g_i_) = 0.5d0 * g_thf_ac(g_i_) + 0.5d0 * g_shf_ac(g_i_
     *)
        enddo
        q2_ac = 0.5d0 * (shf_ac + thf_ac)
C--------
        do g_i_ = 1, g_p_
          g_q2_ap(g_i_) = 0.5d0 * g_thf_ap(g_i_) + 0.5d0 * g_shf_ap(g_i_
     *)
        enddo
        q2_ap = 0.5d0 * (shf_ap + thf_ap)
C--------
C
        do g_i_ = 1, g_p_
          g_j2_ac(g_i_) = (-0.5d0) * g_thf_ac(g_i_) + 0.5d0 * g_shf_ac(g
     *_i_)
        enddo
        j2_ac = 0.5d0 * (shf_ac - thf_ac)
C--------
        do g_i_ = 1, g_p_
          g_j2_ap(g_i_) = (-0.5d0) * g_thf_ap(g_i_) + 0.5d0 * g_shf_ap(g
     *_i_)
        enddo
        j2_ap = 0.5d0 * (shf_ap - thf_ap)
C--------
C
C     compute LEPS functions
C
        do g_i_ = 1, g_p_
          g_j1_ac(g_i_) = g_j1_ac(g_i_)
        enddo
        j1_ac = j1_ac
C--------
        do g_i_ = 1, g_p_
          g_j2_ac(g_i_) = g_j2_ac(g_i_)
        enddo
        j2_ac = j2_ac
C--------
        do g_i_ = 1, g_p_
          g_j3_ac(g_i_) = g_j3_ac(g_i_)
        enddo
        j3_ac = j3_ac
C--------
        d4_b = -j1_ac + (-j2_ac) + j3_ac + j3_ac
        d8_b = -j3_ac + (-j1_ac) + j2_ac + j2_ac
        d5_b = -j3_ac + (-j2_ac) + j1_ac + j1_ac
        do g_i_ = 1, g_p_
          g_w(g_i_) = d4_b * g_j3_ac(g_i_) + d8_b * g_j2_ac(g_i_) + d5_b
     * * g_j1_ac(g_i_)
        enddo
        w = j1_ac * j1_ac + j2_ac * j2_ac + j3_ac * j3_ac - j1_ac * j2_a
     *c - j2_ac * j3_ac - j3_ac * j1_ac
C--------
C
        d7_v = sqrt(w)

c        if ( w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d7_v)
c        else
c           call ehufDO (9,w, d7_v, d1_p,
c     +g_ehfid,
c     +704)
c        endif
        do g_i_ = 1, g_p_
          g_leps_ac(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ac(g_i_) + g_q2_a
     *c(g_i_) + g_q1_ac(g_i_)
        enddo
        leps_ac = q1_ac + q2_ac + q3_ac - d7_v + 6.122d0
C--------
C
        do g_i_ = 1, g_p_
          g_j1_ap(g_i_) = g_j1_ap(g_i_)
        enddo
        j1_ap = j1_ap
C--------
        do g_i_ = 1, g_p_
          g_j2_ap(g_i_) = g_j2_ap(g_i_)
        enddo
        j2_ap = j2_ap
C--------
        do g_i_ = 1, g_p_
          g_j3_ap(g_i_) = g_j3_ap(g_i_)
        enddo
        j3_ap = j3_ap
C--------
        d4_b = -j1_ap + (-j2_ap) + j3_ap + j3_ap
        d8_b = -j3_ap + (-j1_ap) + j2_ap + j2_ap
        d5_b = -j3_ap + (-j2_ap) + j1_ap + j1_ap
        do g_i_ = 1, g_p_
          g_w(g_i_) = d4_b * g_j3_ap(g_i_) + d8_b * g_j2_ap(g_i_) + d5_b
     * * g_j1_ap(g_i_)
        enddo
        w = j1_ap * j1_ap + j2_ap * j2_ap + j3_ap * j3_ap - j1_ap * j2_a
     *p - j2_ap * j3_ap - j3_ap * j1_ap
C--------
C
        d7_v = sqrt(w)

c        if ( w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d7_v)
c        else
c           call ehufDO (9,w, d7_v, d1_p,
c     +g_ehfid,
c     +746)
c        endif
        do g_i_ = 1, g_p_
          g_leps_ap(g_i_) = (-d1_p) * g_w(g_i_) + g_q3_ap(g_i_) + g_q2_a
     *p(g_i_) + g_q1_ap(g_i_)
        enddo
        leps_ap = q1_ap + q2_ap + q3_ap - d7_v + 5.77096d0
C--------
        do g_i_ = 1, g_p_
          g_cfu22(g_i_) = -g_leps_ap(g_i_) + g_leps_ac(g_i_)
        enddo
        cfu22 = leps_ac - leps_ap
C--------
C
      end
c **********************************************************************
