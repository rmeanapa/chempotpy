!> @brief Implement H3ClO_PIPNN, a potential energy surface for H3ClO.
module chempotpy_h3clo_pipnn_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the H3ClO_PIPNN potential energy surface for H3ClO.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 5
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: v
        real(wp)              :: tx(3,natoms), tg(3,natoms)
        integer               :: iatom, idir, j, istate
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                tx(idir,iatom) = x(iatom,idir)
            enddo
        enddo
        if( igrad == 0 )then
            call evvdvdx(tx, v)
        else
            write(*, *) 'Only energy is available'
        endif
        do istate=1, nstates
            p(istate) = v
        enddo
    end subroutine evaluate

    !> @brief Evaluate evvdvdx.
    subroutine evvdvdx( xcart, v )
        !Subroutine to calculate the average potential energy v and analytical gradient dvdxa
        !Xcart     ==>Cartesian coordinates in angstrom
        !             H1: xcart(1,1), xcart(2,1), xcart(3,1)
        !             H2: xcart(1,2), xcart(2,2), xcart(3,2)
        !             H3: xcart(1,3), xcart(2,3), xcart(3,3)
        !             O4: xcart(1,4), xcart(2,4), xcart(3,4)
        !            Cl5: xcart(1,5), xcart(2,5), xcart(3,5)
        !v         ==>Average potential energy(in eV), v=sum(vpes)/npes
        !vpes      ==>Potential energy(in eV) for each PES
        !dvdxa     ==>Average gradient(in eV/ang), dvdxa=sum(dvdx)/npes
        !dvdx      ==>Gradient(in eV/ang) for each PES
        !ndriv     ==>ndriv=1 to compute the analytical gradient, otherwise ndriv=0
        !PES1: totrmse=1.463, errmax=78.909, trainrmse=1.026, validrmse=2.217, testrmse=4.350 meV
        !PES2: totrmse=1.551, errmax=70.610, trainrmse=1.310, validrmse=2.332, testrmse=3.440 meV
        !PES3: totrmse=1.724, errmax=89.420, trainrmse=1.438, validrmse=3.773, testrmse=2.827 meV
        !Average potential energy surface: totormse=1.266 meV
        use chempotpy_geometry,               only: all_pair_distances
        use chempotpy_h3clo_pipnn_parameters, only: natom, nbond, network_count, ninput, nmorse
        use chempotpy_kinds,                  only: wp
        implicit none
        integer  :: i, j, ndriv
        real(wp) :: v, vpes(network_count), c(0:ninput), dvdx(1:natom*3,network_count)
        real(wp) :: xcart(1:3,1:natom), m(0:nmorse-1), xmorse(1:nbond), r(1:nbond)
        real(wp) :: dvdr(1:nbond), x(1:natom*3), drdx(1:nbond,1:natom*3)
        real(wp) :: txinput(1:ninput), dvdg(1:ninput,network_count), p(0:ninput)
        real(wp) :: dpdr(1:nbond,0:ninput), dvdxa(1:natom*3)
        call all_pair_distances(transpose(xcart), r)
        xmorse = exp(-r)
        call bemsav(xmorse, m, p)
        do j=1, ninput
            txinput(j) = p(j)
        enddo
        call pot3a(txinput, vpes)
        v = 0e0_wp
        do i=1, network_count
            v = v+vpes(i)
        enddo
        v = v/network_count
    end subroutine evvdvdx

    !> @brief Evaluate pot3a.
    subroutine pot3a( x, vpot3 )
        use chempotpy_h3clo_pipnn_parameters, only: activation_function, hidden_1_biases, hidden_1_count, &
            & hidden_2_biases, hidden_2_count, hidden_hidden_weights, hidden_output_weights, &
            & input_hidden_weights, network_count, ninput, nscale, output_biases, scale_averages, &
            & scale_deltas
        use chempotpy_kinds, only: wp
        implicit none
        real(wp), intent(in)  :: x(ninput)
        real(wp), intent(out) :: vpot3(network_count)
        integer               :: i, inode1, inode2, m
        real(wp)              :: input_values(ninput), hidden_1_values(hidden_1_count)
        real(wp)              :: hidden_2_values(hidden_2_count), output_value
        do m=1, network_count
            do i=1, ninput
                input_values(i) = (x(i)-scale_averages(i, m))/scale_deltas(i, m)
            enddo
            do inode1=1, hidden_1_count
                hidden_1_values(inode1) = hidden_1_biases(inode1, m)
                do inode2=1, ninput
                    hidden_1_values(inode1) = hidden_1_values(inode1) +input_values(inode2)* &
                        & input_hidden_weights(inode2, inode1, m)
                enddo
                hidden_1_values(inode1) = tranfun(hidden_1_values(inode1), activation_function)
            enddo
            do inode1=1, hidden_2_count
                hidden_2_values(inode1) = hidden_2_biases(inode1, m)
                do inode2=1, hidden_1_count
                    hidden_2_values(inode1) = hidden_2_values(inode1) +hidden_1_values(inode2)* &
                        & hidden_hidden_weights(inode2, inode1, m)
                enddo
                hidden_2_values(inode1) = tranfun(hidden_2_values(inode1), activation_function)
            enddo
            output_value = output_biases(1, m)
            do inode2=1, hidden_2_count
                output_value = output_value+hidden_2_values(inode2)*hidden_output_weights(inode2, 1, m)
            enddo
            vpot3(m) = output_value*scale_deltas(nscale, m)+scale_averages(nscale, m)
        enddo
    end subroutine pot3a

    !> @brief Evaluate tranfun.
    function tranfun( x, ifunc )
        use chempotpy_kinds, only: wp
        implicit none
        integer  :: ifunc
        real(wp) :: tranfun, x
        if( ifunc .eq. 1 )then
            tranfun = tanh(x)
        elseif( ifunc .eq. 2 )then
            tranfun = 1e0_wp/(1e0_wp+exp(-x))
        elseif( ifunc .eq. 3 )then
            tranfun = x
        endif
    end function tranfun

    !> @brief Evaluate LOGSIG.
    function logsig( x )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x, logsig
        logsig = 1.e0_wp/(1.e0_wp+exp(-x))
    end function logsig

    !      real*8 dvdx(1:natom*3,npes),xcart(1:3,1:natom)
    !      real*8 drdx(1:nbond,1:natom*3),x(1:ninput)
    !      real*8 dvdr(1:nbond),dpdr(1:nbond,0:ninput),dvdp(0:ninput)
    !      real*8 dgdx(1:natom*3,1:ninput),dgdr(1:nbond,1:ninput)
    !      real*8 dvdxa(1:natom*3),dvdg(1:ninput,npes),rtmp
    !      dgdx=0d0
    !      dgdr(1:3,1:12)=dpdr(1:3,1:12)
    !         dgdx(i,j)=dgdx(i,j) + dgdr(k,j)*drdx(k,i)
    !      dvdx=0d0;dvdxa=0
    !         dvdx(i,k)=dvdx(i,k) + dvdg(j,k)*dgdx(i,j)
    !        dvdxa(i)=dvdxa(i)+dvdx(i,k)
    !      dvdxa=dvdxa/dble(npes)
    !      real*8 r(1:nbond),dmsdr(1:nbond,1:nbond),dmdr(1:nbond,0:4)
    !      real*8 dvdr(1:nbond),dpdr(1:nbond,0:ninput),c(0:ninput)
    !      real*8 m(0:nmorse-1),p(0:ninput)
    !      dvdr(:)=0d0
    !      real*8 r(1:nbond),xcart(1:3,1:natom),drdx(1:nbond,1:natom*3)
    !      real*8 x(1:natom*3)
    !      drdx(:,:)=0d0
    !!      do i=1,natom-1
    !!       do j=i+1,natom
    !!        k=k+1
    !!        drdx(k,(i-1)*3+1)=(x((i-1)*3+1)-x((j-1)*3+1))/r(k)
    !!        drdx(k,(i-1)*3+2)=(x((i-1)*3+2)-x((j-1)*3+2))/r(k)
    !!        drdx(k,(i-1)*3+3)=(x((i-1)*3+3)-x((j-1)*3+3))/r(k)
    !!        drdx(k,(j-1)*3+1)=-1d0*drdx(k,(i-1)*3+1)
    !!        drdx(k,(j-1)*3+2)=-1d0*drdx(k,(i-1)*3+2)
    !!        drdx(k,(j-1)*3+3)=-1d0*drdx(k,(i-1)*3+3)
    !!       enddo
    !!      enddo
    !!dr1dx
    !      drdx(1,1)=(x(1)-x(4))/r(1)
    !      drdx(1,2)=(x(2)-x(5))/r(1)
    !      drdx(1,3)=(x(3)-x(6))/r(1)
    !      drdx(1,4)=-drdx(1,1)
    !      drdx(1,5)=-drdx(1,2)
    !      drdx(1,6)=-drdx(1,3)
    !!dr2dx
    !      drdx(2,1)=(x(1)-x(7))/r(2)
    !      drdx(2,2)=(x(2)-x(8))/r(2)
    !      drdx(2,3)=(x(3)-x(9))/r(2)
    !      drdx(2,7)=-drdx(2,1)
    !      drdx(2,8)=-drdx(2,2)
    !      drdx(2,9)=-drdx(2,3)
    !!dr3dx
    !      drdx(3,4)=(x(4)-x(7))/r(3)
    !      drdx(3,5)=(x(5)-x(8))/r(3)
    !      drdx(3,6)=(x(6)-x(9))/r(3)
    !      drdx(3,7)=-drdx(3,4)
    !      drdx(3,8)=-drdx(3,5)
    !      drdx(3,9)=-drdx(3,6)
    !      real*8 dmsdr(1:nbond,1:nbond),r(1:nbond)
    !      dmsdr(:,:)=0d0
    !       dmsdr(i,i)=-(1d0/1d0)*dexp(-r(i)/1.0d0)
    !!      dmsdr(i,i)=-(1d0/alpha)*dexp(-r(i)/alpha)
    !      real*8 dmsdr(1:nbond,1:nbond),dmdr(1:nbond,0:4),m(0:nmorse-1)
    !       dmdr(i,0)=0.D0
    !       dmdr(i,1)=dmsdr(i,3)
    !       dmdr(i,2)=dmsdr(i,2)
    !       dmdr(i,3)=dmsdr(i,1)
    !       dmdr(i,4)=dmdr(i,1)*m(2)+m(1)*dmdr(i,2)
    !      real*8 dmdr(1:nbond,0:4),dpdr(1:nbond,0:ninput)
    !      real*8 p(0:ninput)
    !       dpdr(i,0)=dmdr(i,0)
    !       dpdr(i,1)=dmdr(i,1)+dmdr(i,2)
    !       dpdr(i,2)=dmdr(i,3)
    !       dpdr(i,3)=dmdr(i,4)
    !       dpdr(i,4)=dpdr(i,2)*P(1)+P(2)*dpdr(i,1)
    !       dpdr(i,5)=dpdr(i,1)*p(1)+p(1)*dpdr(i,1)-dpdr(i,3)-dpdr(i,3)
    !       dpdr(i,6)=dpdr(i,2)*P(2)+P(2)*dpdr(i,2)
    !       dpdr(i,7)=dpdr(i,2)*P(3)+P(2)*dpdr(i,3)
    !       dpdr(i,8)=dpdr(i,3)*P(1)+P(3)*dpdr(i,1)
    !       dpdr(i,9)=dpdr(i,2)*P(5)+P(2)*dpdr(i,5)
    !       dpdr(i,10)=dpdr(i,2)*P(4)+P(2)*dpdr(i,4)
    !       dpdr(i,11)=dpdr(i,1)*P(5)+P(1)*dpdr(i,5)-dpdr(i,8)
    !       dpdr(i,12)=dpdr(i,2)*P(6)+P(2)*dpdr(i,6)
    !> @brief Evaluate bemsav.
    subroutine bemsav( x, m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(1:10)
        real(wp) :: p(0:230)
        real(wp) :: m(0:109)
        call evmono(x, m)
        call evpoly(m, p)
    end subroutine bemsav

    !> @brief Evaluate evmono.
    subroutine evmono( x, m )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: x(1:10)
        real(wp) :: m(0:109)
        m(0)   = 1.0e0_wp
        m(1)   = x(10)
        m(2)   = x(9)
        m(3)   = x(7)
        m(4)   = x(4)
        m(5)   = x(8)
        m(6)   = x(6)
        m(7)   = x(3)
        m(8)   = x(5)
        m(9)   = x(2)
        m(10)  = x(1)
        m(11)  = m(2)*m(3)
        m(12)  = m(2)*m(4)
        m(13)  = m(3)*m(4)
        m(14)  = m(3)*m(5)
        m(15)  = m(2)*m(6)
        m(16)  = m(4)*m(5)
        m(17)  = m(4)*m(6)
        m(18)  = m(2)*m(7)
        m(19)  = m(3)*m(7)
        m(20)  = m(5)*m(6)
        m(21)  = m(5)*m(7)
        m(22)  = m(6)*m(7)
        m(23)  = m(4)*m(8)
        m(24)  = m(3)*m(9)
        m(25)  = m(2)*m(10)
        m(26)  = m(7)*m(8)
        m(27)  = m(6)*m(9)
        m(28)  = m(5)*m(10)
        m(29)  = m(8)*m(9)
        m(30)  = m(8)*m(10)
        m(31)  = m(9)*m(10)
        m(32)  = m(2)*m(13)
        m(33)  = m(3)*m(16)
        m(34)  = m(2)*m(17)
        m(35)  = m(2)*m(19)
        m(36)  = m(4)*m(20)
        m(37)  = m(3)*m(21)
        m(38)  = m(2)*m(22)
        m(39)  = m(5)*m(22)
        m(40)  = m(4)*m(26)
        m(41)  = m(3)*m(27)
        m(42)  = m(2)*m(28)
        m(43)  = m(2)*m(23)
        m(44)  = m(3)*m(23)
        m(45)  = m(2)*m(24)
        m(46)  = m(4)*m(24)
        m(47)  = m(3)*m(25)
        m(48)  = m(4)*m(25)
        m(49)  = m(5)*m(26)
        m(50)  = m(6)*m(26)
        m(51)  = m(5)*m(27)
        m(52)  = m(7)*m(27)
        m(53)  = m(6)*m(28)
        m(54)  = m(7)*m(28)
        m(55)  = m(3)*m(29)
        m(56)  = m(4)*m(29)
        m(57)  = m(2)*m(30)
        m(58)  = m(4)*m(30)
        m(59)  = m(2)*m(31)
        m(60)  = m(3)*m(31)
        m(61)  = m(6)*m(29)
        m(62)  = m(7)*m(29)
        m(63)  = m(5)*m(30)
        m(64)  = m(7)*m(30)
        m(65)  = m(5)*m(31)
        m(66)  = m(6)*m(31)
        m(67)  = m(8)*m(31)
        m(68)  = m(2)*m(36)
        m(69)  = m(3)*m(36)
        m(70)  = m(2)*m(37)
        m(71)  = m(3)*m(38)
        m(72)  = m(4)*m(37)
        m(73)  = m(4)*m(38)
        m(74)  = m(2)*m(40)
        m(75)  = m(3)*m(40)
        m(76)  = m(2)*m(41)
        m(77)  = m(4)*m(41)
        m(78)  = m(3)*m(42)
        m(79)  = m(4)*m(42)
        m(80)  = m(4)*m(49)
        m(81)  = m(4)*m(50)
        m(82)  = m(3)*m(51)
        m(83)  = m(3)*m(52)
        m(84)  = m(2)*m(53)
        m(85)  = m(2)*m(54)
        m(86)  = m(3)*m(49)
        m(87)  = m(2)*m(50)
        m(88)  = m(4)*m(51)
        m(89)  = m(2)*m(52)
        m(90)  = m(4)*m(53)
        m(91)  = m(3)*m(54)
        m(92)  = m(3)*m(56)
        m(93)  = m(2)*m(58)
        m(94)  = m(2)*m(60)
        m(95)  = m(4)*m(61)
        m(96)  = m(3)*m(62)
        m(97)  = m(4)*m(63)
        m(98)  = m(2)*m(64)
        m(99)  = m(3)*m(65)
        m(100) = m(2)*m(66)
        m(101) = m(6)*m(62)
        m(102) = m(5)*m(64)
        m(103) = m(5)*m(66)
        m(104) = m(3)*m(61)
        m(105) = m(4)*m(62)
        m(106) = m(2)*m(63)
        m(107) = m(4)*m(64)
        m(108) = m(2)*m(65)
        m(109) = m(3)*m(66)
    end subroutine evmono

    !> @brief Evaluate evpoly.
    subroutine evpoly( m, p )
        use chempotpy_kinds, only: wp
        implicit none
        real(wp) :: p(0:230)
        real(wp) :: m(0:109)
        p(0)   = m(0)
        p(1)   = m(1)
        p(2)   = m(2) + m(3) + m(4)
        p(3)   = m(5) + m(6) + m(7)
        p(4)   = m(8) + m(9) + m(10)
        p(5)   = p(1)*p(2)
        p(6)   = m(11) + m(12) + m(13)
        p(7)   = p(1)*p(3)
        p(8)   = m(14) + m(15) + m(16) + m(17) + m(18) + m(19)
        p(9)   = m(20) + m(21) + m(22)
        p(10)  = p(2)*p(3) - p(8)
        p(11)  = p(1)*p(4)
        p(12)  = m(23) + m(24) + m(25)
        p(13)  = m(26) + m(27) + m(28)
        p(14)  = p(2)*p(4) - p(12)
        p(15)  = p(3)*p(4) - p(13)
        p(16)  = m(29) + m(30) + m(31)
        p(17)  = p(1)*p(1)
        p(18)  = p(2)*p(2) - p(6) - p(6)
        p(19)  = p(3)*p(3) - p(9) - p(9)
        p(20)  = p(4)*p(4) - p(16) - p(16)
        p(21)  = p(1)*p(6)
        p(22)  = m(32)
        p(23)  = p(1)*p(8)
        p(24)  = m(33) + m(34) + m(35)
        p(25)  = p(1)*p(9)
        p(26)  = m(36) + m(37) + m(38)
        p(27)  = m(39)
        p(28)  = p(1)*p(10)
        p(29)  = p(3)*p(6) - p(24)
        p(30)  = p(2)*p(9) - p(26)
        p(31)  = p(1)*p(12)
        p(32)  = p(1)*p(13)
        p(33)  = m(40) + m(41) + m(42)
        p(34)  = p(1)*p(14)
        p(35)  = m(43) + m(44) + m(45) + m(46) + m(47) + m(48)
        p(36)  = p(2)*p(13) - p(33)
        p(37)  = p(4)*p(6) - p(35)
        p(38)  = p(1)*p(15)
        p(39)  = p(3)*p(12) - p(33)
        p(40)  = m(49) + m(50) + m(51) + m(52) + m(53) + m(54)
        p(41)  = p(4)*p(8) - p(39) - p(36)
        p(42)  = p(4)*p(9) - p(40)
        p(43)  = p(4)*p(10) - p(33)
        p(44)  = p(1)*p(16)
        p(45)  = m(55) + m(56) + m(57) + m(58) + m(59) + m(60)
        p(46)  = m(61) + m(62) + m(63) + m(64) + m(65) + m(66)
        p(47)  = m(67)
        p(48)  = p(2)*p(16) - p(45)
        p(49)  = p(3)*p(16) - p(46)
        p(50)  = p(1)*p(5)
        p(51)  = p(1)*p(18)
        p(52)  = p(2)*p(6) - p(22) - p(22) - p(22)
        p(53)  = p(1)*p(7)
        p(54)  = p(2)*p(8) - p(29) - p(24) - p(24)
        p(55)  = p(2)*p(10) - p(29)
        p(56)  = p(1)*p(19)
        p(57)  = p(3)*p(8) - p(30) - p(26) - p(26)
        p(58)  = p(3)*p(9) - p(27) - p(27) - p(27)
        p(59)  = p(2)*p(19) - p(57)
        p(60)  = p(1)*p(11)
        p(61)  = p(2)*p(12) - p(35)
        p(62)  = p(3)*p(13) - p(40)
        p(63)  = p(4)*p(18) - p(61)
        p(64)  = p(4)*p(19) - p(62)
        p(65)  = p(1)*p(20)
        p(66)  = p(4)*p(12) - p(45)
        p(67)  = p(4)*p(13) - p(46)
        p(68)  = p(2)*p(20) - p(66)
        p(69)  = p(3)*p(20) - p(67)
        p(70)  = p(4)*p(16) - p(47) - p(47) - p(47)
        p(71)  = p(1)*p(17)
        p(72)  = p(2)*p(18) - p(52)
        p(73)  = p(3)*p(19) - p(58)
        p(74)  = p(4)*p(20) - p(70)
        p(75)  = p(1)*p(22)
        p(76)  = p(1)*p(24)
        p(77)  = p(1)*p(26)
        p(78)  = p(1)*p(27)
        p(79)  = p(1)*p(29)
        p(80)  = p(22)*p(3)
        p(81)  = p(1)*p(30)
        p(82)  = m(68) + m(69) + m(70) + m(71) + m(72) + m(73)
        p(83)  = p(27)*p(2)
        p(84)  = p(6)*p(9) - p(82)
        p(85)  = p(1)*p(33)
        p(86)  = p(1)*p(35)
        p(87)  = p(1)*p(36)
        p(88)  = m(74) + m(75) + m(76) + m(77) + m(78) + m(79)
        p(89)  = p(1)*p(37)
        p(90)  = p(22)*p(4)
        p(91)  = p(6)*p(13) - p(88)
        p(92)  = p(1)*p(39)
        p(93)  = p(1)*p(40)
        p(94)  = m(80) + m(81) + m(82) + m(83) + m(84) + m(85)
        p(95)  = p(1)*p(41)
        p(96)  = p(4)*p(24) - p(91)
        p(97)  = m(86) + m(87) + m(88) + m(89) + m(90) + m(91)
        p(98)  = p(1)*p(42)
        p(99)  = p(4)*p(26) - p(97)
        p(100) = p(27)*p(4)
        p(101) = p(1)*p(43)
        p(102) = p(3)*p(35) - p(96) - p(88)
        p(103) = p(2)*p(40) - p(97) - p(94)
        p(104) = p(3)*p(37) - p(91)
        p(105) = p(2)*p(42) - p(99)
        p(106) = p(1)*p(45)
        p(107) = m(92) + m(93) + m(94)
        p(108) = p(1)*p(46)
        p(109) = m(95) + m(96) + m(97) + m(98) + m(99) + m(100)
        p(110) = m(101) + m(102) + m(103)
        p(111) = m(104) + m(105) + m(106) + m(107) + m(108) + m(109)
        p(112) = p(1)*p(47)
        p(113) = p(1)*p(48)
        p(114) = p(6)*p(16) - p(107)
        p(115) = p(2)*p(46) - p(111) - p(109)
        p(116) = p(47)*p(2)
        p(117) = p(1)*p(49)
        p(118) = p(3)*p(45) - p(111) - p(109)
        p(119) = p(9)*p(16) - p(110)
        p(120) = p(47)*p(3)
        p(121) = p(2)*p(49) - p(118)
        p(122) = p(1)*p(21)
        p(123) = p(1)*p(52)
        p(124) = p(22)*p(2)
        p(125) = p(1)*p(23)
        p(126) = p(1)*p(54)
        p(127) = p(2)*p(24) - p(80)
        p(128) = p(1)*p(25)
        p(129) = p(2)*p(26) - p(82)
        p(130) = p(1)*p(28)
        p(131) = p(6)*p(8) - p(80) - p(127) - p(80)
        p(132) = p(1)*p(55)
        p(133) = p(6)*p(10) - p(80)
        p(134) = p(9)*p(18) - p(129)
        p(135) = p(1)*p(57)
        p(136) = p(3)*p(24) - p(82)
        p(137) = p(1)*p(58)
        p(138) = p(3)*p(26) - p(83)
        p(139) = p(27)*p(3)
        p(140) = p(8)*p(9) - p(83) - p(138) - p(83)
        p(141) = p(1)*p(59)
        p(142) = p(6)*p(19) - p(136)
        p(143) = p(9)*p(10) - p(83)
        p(144) = p(1)*p(31)
        p(145) = p(1)*p(61)
        p(146) = p(1)*p(32)
        p(147) = p(2)*p(33) - p(88)
        p(148) = p(1)*p(62)
        p(149) = p(3)*p(33) - p(94)
        p(150) = p(1)*p(34)
        p(151) = p(6)*p(12) - p(90)
        p(152) = p(2)*p(62) - p(149)
        p(153) = p(1)*p(63)
        p(154) = p(2)*p(35) - p(90) - p(151) - p(90)
        p(155) = p(13)*p(18) - p(147)
        p(156) = p(2)*p(37) - p(90)
        p(157) = p(1)*p(38)
        p(158) = p(3)*p(61) - p(147)
        p(159) = p(9)*p(13) - p(100)
        p(160) = p(2)*p(41) - p(104) - p(96)
        p(161) = p(4)*p(55) - p(147)
        p(162) = p(1)*p(64)
        p(163) = p(12)*p(19) - p(149)
        p(164) = p(3)*p(40) - p(100) - p(159) - p(100)
        p(165) = p(3)*p(41) - p(105) - p(97)
        p(166) = p(3)*p(42) - p(100)
        p(167) = p(4)*p(59) - p(149)
        p(168) = p(1)*p(44)
        p(169) = p(2)*p(45) - p(114) - p(107) - p(107)
        p(170) = p(3)*p(46) - p(119) - p(110) - p(110)
        p(171) = p(2)*p(48) - p(114)
        p(172) = p(3)*p(49) - p(119)
        p(173) = p(1)*p(66)
        p(174) = p(1)*p(67)
        p(175) = p(4)*p(33) - p(111)
        p(176) = p(1)*p(68)
        p(177) = p(12)*p(14) - p(114) - p(169)
        p(178) = p(2)*p(67) - p(175)
        p(179) = p(4)*p(37) - p(114)
        p(180) = p(1)*p(69)
        p(181) = p(3)*p(66) - p(175)
        p(182) = p(13)*p(15) - p(119) - p(170)
        p(183) = p(4)*p(41) - p(118) - p(115)
        p(184) = p(4)*p(42) - p(119)
        p(185) = p(10)*p(20) - p(175)
        p(186) = p(1)*p(70)
        p(187) = p(12)*p(16) - p(116)
        p(188) = p(13)*p(16) - p(120)
        p(189) = p(4)*p(45) - p(116) - p(187) - p(116)
        p(190) = p(4)*p(46) - p(120) - p(188) - p(120)
        p(191) = p(47)*p(4)
        p(192) = p(4)*p(48) - p(116)
        p(193) = p(4)*p(49) - p(120)
        p(194) = p(1)*p(50)
        p(195) = p(1)*p(51)
        p(196) = p(6)*p(6) - p(124) - p(124)
        p(197) = p(1)*p(72)
        p(198) = p(6)*p(18) - p(124)
        p(199) = p(1)*p(53)
        p(200) = p(2)*p(54) - p(131) - p(127)
        p(201) = p(2)*p(55) - p(133)
        p(202) = p(1)*p(56)
        p(203) = p(2)*p(57) - p(142) - p(136) - p(136)
        p(204) = p(9)*p(9) - p(139) - p(139)
        p(205) = p(2)*p(59) - p(142)
        p(206) = p(1)*p(73)
        p(207) = p(3)*p(57) - p(140) - p(138)
        p(208) = p(9)*p(19) - p(139)
        p(209) = p(2)*p(73) - p(207)
        p(210) = p(1)*p(60)
        p(211) = p(2)*p(61) - p(151)
        p(212) = p(3)*p(62) - p(159)
        p(213) = p(4)*p(72) - p(211)
        p(214) = p(4)*p(73) - p(212)
        p(215) = p(1)*p(65)
        p(216) = p(2)*p(66) - p(177)
        p(217) = p(3)*p(67) - p(182)
        p(218) = p(18)*p(20) - p(216)
        p(219) = p(19)*p(20) - p(217)
        p(220) = p(16)*p(16) - p(191) - p(191)
        p(221) = p(1)*p(74)
        p(222) = p(4)*p(66) - p(187)
        p(223) = p(4)*p(67) - p(188)
        p(224) = p(2)*p(74) - p(222)
        p(225) = p(3)*p(74) - p(223)
        p(226) = p(16)*p(20) - p(191)
        p(227) = p(1)*p(71)
        p(228) = p(2)*p(72) - p(198)
        p(229) = p(3)*p(73) - p(208)
        p(230) = p(4)*p(74) - p(226)
    end subroutine evpoly
end module chempotpy_h3clo_pipnn_pes

!> @brief Expose H3ClO_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_h3clo_pipnn_pes, only: evaluate
    use chempotpy_kinds,           only: wp
    implicit none
    integer, parameter              :: natoms  = 5
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    real(wp), intent(in)            :: x(natoms,3)
    real(wp), intent(out)           :: p(nstates)
    real(wp), intent(out)           :: g(nstates,natoms,3)
    real(wp), intent(out)           :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
