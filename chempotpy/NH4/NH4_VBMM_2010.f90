!> @brief Implement NH4_VBMM_2010, a potential energy surface for NH4.
module chempotpy_nh4_vbmm_2010_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the NH4_VBMM_2010 potential energy surface for NH4.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units,  only: hartree_to_electronvolt
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 5
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter    :: n3atom  = 75
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: a1ch
        real(wp)              :: a1s
        real(wp)              :: a2s
        real(wp)              :: a3s
        real(wp)              :: aa1
        real(wp)              :: aa2
        real(wp)              :: aa3
        real(wp)              :: aa4
        real(wp)              :: acb
        real(wp)              :: ahh
        real(wp)              :: ak
        real(wp)              :: aphi
        real(wp)              :: atheta
        real(wp)              :: b1ch
        real(wp)              :: b1s
        real(wp)              :: b2s
        real(wp)              :: b3s
        real(wp)              :: bk
        real(wp)              :: bphi
        real(wp)              :: btheta
        real(wp)              :: c1ch
        real(wp)              :: input_cartesian(natom,3)
        real(wp)              :: cartesian_coordinates(natom,3)
        real(wp)              :: cphi
        real(wp)              :: ctheta
        real(wp)              :: d1cb
        real(wp)              :: d1ch
        real(wp)              :: d1hh
        real(wp)              :: d3cb
        real(wp)              :: d3ch
        real(wp)              :: d3hh
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: dfdc(4,4,4)
        real(wp)              :: dfdelta(4,4)
        real(wp)              :: dfdh(4,4,4)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: dhdelta(4,4)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        real(wp)              :: doptheta0(4,4,4)
        real(wp)              :: ds1(4)
        real(wp)              :: ds2(4)
        real(wp)              :: ds3(4)
        real(wp)              :: dsphi(4)
        real(wp)              :: dstheta(4)
        real(wp)              :: dtheta0(4,4,4)
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: ground_state_energy
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: ezero(isurf+1)
        real(wp)              :: f1(4)
        real(wp)              :: fch3
        real(wp)              :: fdelta(4)
        real(wp)              :: fk0(4,4)
        real(wp)              :: fkinf
        real(wp)              :: hch3
        real(wp)              :: hdelta(4)
        integer               :: coordinate_representation
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        integer               :: atom_count
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: nc(3)
        integer               :: requested_derivative_order
        integer               :: nh(4,3)
        integer               :: nhb(3)
        integer               :: nnb
        integer               :: nnc
        integer               :: nnh(4)
        real(wp)              :: optau
        real(wp)              :: optheta0(4,4)
        real(wp)              :: pdot(150)
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: ppito
        real(wp)              :: q(150)
        real(wp)              :: r(n3atom)
        real(wp)              :: r0cb
        real(wp)              :: r0ch
        real(wp)              :: r0hh
        real(wp)              :: rbh(4)
        real(wp)              :: rcb
        real(wp)              :: rch(4)
        character(len=75)     :: ref(5)
        real(wp)              :: s1(4)
        real(wp)              :: s2(4)
        real(wp)              :: s3(4)
        real(wp)              :: sphi(4)
        real(wp)              :: stheta(4)
        real(wp)              :: tau
        real(wp)              :: taunh2
        real(wp)              :: tbh(4,3)
        real(wp)              :: tcb(3)
        real(wp)              :: tch(4,3)
        real(wp)              :: theta0(4,4)
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        a1ch                         = 0.e0_wp; a1s = 0.e0_wp; a2s = 0.e0_wp; a3s = 0.e0_wp; aa1 = 0.e0_wp; aa2 = 0.e0_wp; aa3 = 0.e0_wp; aa4 = 0.e0_wp; acb = 0.e0_wp
        ahh                          = 0.e0_wp; ak = 0.e0_wp; aphi = 0.e0_wp; atheta = 0.e0_wp; b1ch = 0.e0_wp; b1s = 0.e0_wp; b2s = 0.e0_wp
        b3s                          = 0.e0_wp; bk = 0.e0_wp; bphi = 0.e0_wp; btheta = 0.e0_wp; c1ch = 0.e0_wp; input_cartesian = 0.e0_wp; cartesian_coordinates = 0.e0_wp
        cphi                         = 0.e0_wp; ctheta = 0.e0_wp; d1cb = 0.e0_wp; d1ch = 0.e0_wp; d1hh = 0.e0_wp; d3cb = 0.e0_wp
        d3ch                         = 0.e0_wp; d3hh = 0.e0_wp; excited_state_internal_gradients = 0.e0_wp; ground_state_internal_gradient = 0.e0_wp; diabatic_internal_gradients = 0.e0_wp; output_excited_state_gradients = 0.e0_wp; excited_state_cartesian_gradients = 0.e0_wp
        dfdc                         = 0.e0_wp; dfdelta = 0.e0_wp; dfdh = 0.e0_wp; output_ground_state_gradient = 0.e0_wp; ground_state_cartesian_gradient = 0.e0_wp; dhdelta = 0.e0_wp; output_diabatic_gradients = 0.e0_wp
        diabatic_cartesian_gradients = 0.e0_wp; doptheta0 = 0.e0_wp; ds1 = 0.e0_wp; ds2 = 0.e0_wp; ds3 = 0.e0_wp; dsphi = 0.e0_wp; dstheta = 0.e0_wp
        dtheta0                      = 0.e0_wp; excited_state_energies = 0.e0_wp; ground_state_energy = 0.e0_wp; diabatic_couplings = 0.e0_wp; ezero = 0.e0_wp; f1 = 0.e0_wp; fch3 = 0.e0_wp; fdelta = 0.e0_wp
        fk0                          = 0.e0_wp; fkinf = 0.e0_wp; hch3 = 0.e0_wp; hdelta = 0.e0_wp; coordinate_representation = 0
        last_requested_excited_state = 0; last_requested_coupling = 0; first_requested_excited_state = 0; first_requested_coupling = 0; atom_count = 0; maximum_derivative_order = 0; maximum_excited_state = 0; surface_requests = 0; nc = 0
        requested_derivative_order   = 0; nh = 0; nhb = 0; nnb = 0; nnc = 0; nnh = 0; optau = 0.e0_wp; optheta0 = 0.e0_wp
        pdot                         = 0.e0_wp; output_excited_state_energies = 0.e0_wp; output_ground_state_energy = 0.e0_wp; output_diabatic_couplings = 0.e0_wp; ppito = 0.e0_wp; q = 0.e0_wp; r = 0.e0_wp; r0cb = 0.e0_wp
        r0ch                         = 0.e0_wp; r0hh = 0.e0_wp; rbh = 0.e0_wp; rcb = 0.e0_wp; rch = 0.e0_wp
        ref                          = ' '
        s1                           = 0.e0_wp
        s2                           = 0.e0_wp
        s3                           = 0.e0_wp
        sphi                         = 0.e0_wp
        stheta                       = 0.e0_wp
        tau                          = 0.e0_wp
        taunh2                       = 0.e0_wp
        tbh                          = 0.e0_wp
        tcb                          = 0.e0_wp
        tch                          = 0.e0_wp
        theta0                       = 0.e0_wp
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 1
        coordinate_representation  = 1
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 5
        nnc                        = 2
        nnb                        = 5
        nnh                        = [3, 4, 1, 0]
        r0ch                       = 1.01162e0_wp
        d1ch                       = 116.458e0_wp
        d3ch                       = 49.002e0_wp
        a1ch                       = 2.160210e0_wp
        b1ch                       = -0.161998e0_wp
        c1ch                       = 22.99980e0_wp
        r0hh                       = 0.74191e0_wp
        d1hh                       = 111.458e0_wp
        d3hh                       = 42.596e0_wp
        ahh                        = 1.9607e0_wp
        r0cb                       = 1.24562e0_wp
        d1cb                       = 46.253e0_wp
        d3cb                       = 12.241e0_wp
        acb                        = 1.2959012e0_wp
        a3s                        = 0.1897001e0_wp
        b3s                        = +2.0088000e0_wp
        aphi                       = 0.1715781e0_wp
        bphi                       = 0.2836958e0_wp
        cphi                       = 1.7508321e0_wp
        atheta                     = 0.7526130e0_wp
        btheta                     = 0.6197119e0_wp
        ctheta                     = 1.2100210e0_wp
        fch3                       = -0.0200200e0_wp
        hch3                       = -0.0200900e0_wp
        fkinf                      = 0.6002176e0_wp
        ak                         = 0.0120105e0_wp
        bk                         = 0.1145600e0_wp
        aa1                        = 5.069420e0_wp
        aa2                        = 0.008463e0_wp
        aa3                        = 0.097621e0_wp
        aa4                        = 5.432432e0_wp
        tau                        = 1.9022630e0_wp
        taunh2                     = 1.8046704e0_wp
        optau                      = 1.1847670e0_wp
        v                          = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        input_cartesian            = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                input_cartesian(iatom,idir) = x(iatom,idir)
            enddo
        enddo
        call prepare_local_context
        call pot
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = output_ground_state_energy*hartree_to_electronvolt
            enddo
        elseif( igrad == 1 )then
            do istate=1, nstates
                p(istate) = output_ground_state_energy*hartree_to_electronvolt
            enddo
            do iatom=1, natoms
                do idir=1, 3
                    g(1,iatom,idir) = output_ground_state_gradient(iatom,idir)*hartree_per_bohr_to_ev_per_angstrom
                enddo
            enddo
        elseif( igrad == 2 )then
            if( present(status) ) call status%set(status_derivative_unavailable, &
                'requested derivatives above first order are unavailable')
        endif

        contains

        !> @brief Evaluate POT for the NH4_VBMM_2010 potential.
        subroutine pot
            !        This potential is written such that:
            !                       X(1)  - X(3)  : X, Y, Z for H1
            !                       X(4)  - X(6)  : X, Y, Z for N
            !                       X(7)  - X(9)  : X, Y, Z for H2
            !                       X(10) - X(12) : X, Y, Z for H3
            !                       X(13) - X(15) : X, Y, Z for H4
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: en, vip, vstr
            integer  :: i
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            ! Changing units from angstrom to bohr and initialize
            do i = 1, 15
                q(i)    = r(i)
                pdot(i) = 0.e0_wp
            enddo
            !  calculate relative coordinates and bond lengths
            en = 0.0e0_wp
            call coorden
            !  calculate switching functions
            call switchf
            !  calculate reference angles and their derivatives
            call refangles
            call oprefangles
            !  calculate stretching potential
            call stretch(vstr)
            !  calculate in plane bending potential
            call ipbend(vip)
            !  total potential energy is vstr+vop+vip
            en = vstr+vip
            ! 0.03812 conversion factor from 10(5) j/mol to hartrees
            en                  = en*0.03812e0_wp
            ground_state_energy = en
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                ! 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
                do i=1, 15
                    ground_state_internal_gradient(i) = pdot(i)*0.0201723e0_wp
                enddo
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate coorden for the NH4_VBMM_2010 potential.
        subroutine coorden
            !  calculates relative coordinates and bond lengths
            implicit none
            integer :: i, ind
            !  calculate relative coordinates
            do ind=1, 3
                tcb(ind) = q(nc(ind))-q(nhb(ind))
                do i=1, 3
                    tch(i,ind) = q(nc(ind))-q(nh(i,ind))
                    tbh(i,ind) = q(nhb(ind))-q(nh(i,ind))
                enddo
            enddo
            !  calculate bond lengths
            rcb = sqrt(tcb(1)*tcb(1)+tcb(2)*tcb(2)+tcb(3)*tcb(3))
            do i=1, 3
                rch(i) = sqrt(tch(i,1)*tch(i,1)+tch(i,2)*tch(i,2)+ tch(i,3)*tch(i,3))
                rbh(i) = sqrt(tbh(i,1)*tbh(i,1)+tbh(i,2)*tbh(i,2)+ tbh(i,3)*tbh(i,3))
            enddo
        end subroutine coorden

        !> @brief Evaluate refangles for the NH4_VBMM_2010 potential.
        subroutine refangles
            !  subroutine calculates reference angles for the "in-plane" potential
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: halfpi, twopi
            integer  :: i, j, k
            !      tau=acos(-1.0d0/3.0d0)
            !      pi=4.0d0*atan(1.0d0)
            halfpi = 0.5e0_wp*PI
            twopi  = 2.0e0_wp*PI
            ppito  = (twopi-taunh2)/2.0e0_wp
            !  set diagonal elements to zero
            do i=1, 3
                theta0(i,i) = 0.0e0_wp
                do k=1, 4
                    dtheta0(i,i,k) = 0.0e0_wp
                enddo
            enddo
            !  calculate reference angles
            theta0(1,2) = tau+(tau-ppito)*(sphi(1)*sphi(2)-1.0e0_wp) +(tau-taunh2)*(stheta(3)-1.0e0_wp)
            theta0(1,3) = tau+(tau-ppito)*(sphi(1)*sphi(3)-1.0e0_wp) +(tau-taunh2)*(stheta(2)-1.0e0_wp)
            theta0(2,3) = tau+(tau-ppito)*(sphi(2)*sphi(3)-1.0e0_wp) +(tau-taunh2)*(stheta(1)-1.0e0_wp)
            !  calculate the derivatives of theta0(i,j) in terms of rch(k)
            !  quantity calulated is dtheta0(i,j,k)
            !  derivatives wrt rch(1)
            dtheta0(1,2,1) = (tau-ppito)*dsphi(1)*sphi(2)
            dtheta0(1,3,1) = (tau-ppito)*dsphi(1)*sphi(3)
            dtheta0(2,3,1) = (tau-taunh2)*dstheta(1)
            !  derivatives wrt rch(2)
            dtheta0(1,2,2) = (tau-ppito)*sphi(1)*dsphi(2)
            dtheta0(1,3,2) = (tau-taunh2)*dstheta(2)
            dtheta0(2,3,2) = (tau-ppito)*dsphi(2)*sphi(3)
            !  derivatives wrt rch(3)
            dtheta0(1,2,3) = (tau-taunh2)*dstheta(3)
            dtheta0(1,3,3) = (tau-ppito)*sphi(1)*dsphi(3)
            dtheta0(2,3,3) = (tau-ppito)*sphi(2)*dsphi(3)
            !  fill in the other half of the matrix
            do i=1, 2
                do j=i+1, 3
                    theta0(j,i) = theta0(i,j)
                    do k=1, 3
                        dtheta0(j,i,k) = dtheta0(i,j,k)
                    enddo
                enddo
            enddo
        end subroutine refangles

        !> @brief Evaluate oprefangles for the NH4_VBMM_2010 potential.
        subroutine oprefangles
            !  subroutine calculates reference angles for the "out-plane" potential
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: halfpi, twopi
            integer  :: i, j, k
            !       pi=4.0d0*atan(1.0d0)
            halfpi = 0.5e0_wp*PI
            twopi  = 2.0e0_wp*PI
            ppito  = (twopi-taunh2)/2.0e0_wp
            !  set diagonal elements to zero
            do i=1, 2
                optheta0(i,i) = 0.0e0_wp
                do k=1, 3
                    doptheta0(i,i,k) = 0.0e0_wp
                enddo
            enddo
            !  calculate reference angles
            optheta0(1,1)=optau+(optau-halfpi)*(sphi(1) -1.0e0_wp)+(optau-halfpi)*(stheta(3)-1.0e0_wp) +(optau- &
                & halfpi)* (stheta(2)- 1.0e0_wp)
            optheta0(2,2)=optau+(optau-halfpi)*(sphi(2) -1.0e0_wp)+(optau-halfpi)*(stheta(3)-1.0e0_wp) +(optau- &
                & halfpi)* (stheta(1)- 1.0e0_wp)
            optheta0(3,3)=optau+(optau-halfpi)*(sphi(3) -1.0e0_wp)+(optau-halfpi)*(stheta(1)-1.0e0_wp) +(optau- &
                & halfpi)* (stheta(2)- 1.0e0_wp)
            optheta0(1,2) = optau+(optau-halfpi)*(sphi(1)*sphi(2) -1.0e0_wp)+(optau-halfpi)*(stheta(3)-1.0e0_wp)
            optheta0(1,3) = optau+(optau-halfpi)*(sphi(1)*sphi(3) -1.0e0_wp)+(optau-halfpi)*(stheta(2)-1.0e0_wp)
            optheta0(2,3) = optau+(optau-halfpi)*(sphi(2)*sphi(3) -1.0e0_wp)+(optau-halfpi)*(stheta(1)-1.0e0_wp)
            !  calculate the derivatives of optheta0(i,j) in terms of rch(k)
            !  quantity calulated is doptheta0(i,j,k)
            !  derivatives wrt rch(1)
            doptheta0(1,1,1) = (optau-halfpi)*dsphi(1)
            doptheta0(2,2,1) = (optau-halfpi)*dstheta(1)
            doptheta0(3,3,1) = (optau-halfpi)*dstheta(1)
            doptheta0(1,2,1) = (optau-halfpi)*dsphi(1)* sphi(2)
            doptheta0(1,3,1) = (optau-halfpi)*dsphi(1)* sphi(3)
            doptheta0(2,3,1) = (optau-halfpi)* dstheta(1)
            !  derivatives wrt rch(2)
            doptheta0(1,1,2) = (optau-halfpi)*dstheta(2)
            doptheta0(2,2,2) = (optau-halfpi)*dsphi(2)
            doptheta0(3,3,2) = (optau-halfpi)*dstheta(2)
            doptheta0(1,2,2) = (optau-halfpi)*sphi(1)* dsphi(2)
            doptheta0(1,3,2) = (optau-halfpi)* dstheta(2)
            doptheta0(2,3,2) = (optau-halfpi)*dsphi(2)* sphi(3)
            !  derivatives wrt rch(3)
            doptheta0(1,1,3) = (optau-halfpi)*dstheta(3)
            doptheta0(2,2,3) = (optau-halfpi)*dstheta(3)
            doptheta0(3,3,3) = (optau-halfpi)*dsphi(3)
            doptheta0(1,2,3) = (optau-halfpi)* dstheta(3)
            doptheta0(1,3,3) = (optau-halfpi)*sphi(1)* dsphi(3)
            doptheta0(2,3,3) = (optau-halfpi)*sphi(2)* dsphi(3)
            !  fill in the other half of the matrix
            do i=1, 2
                do j=i+1, 3
                    optheta0(j,i) = optheta0(i,j)
                    do k=1, 3
                        doptheta0(j,i,k) = doptheta0(i,j,k)
                    enddo
                enddo
            enddo
        end subroutine oprefangles

        !> @brief Evaluate stretch for the NH4_VBMM_2010 potential.
        subroutine stretch( vstr )
            !  subroutine to calculate leps-type stretching potential and its
            !  derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: ach, achdc(3), achdh(4,3), arga, dumach, dumjbh, dumjcb, dumjch, dumjhi, dumjhj, dumqbh, &
                & dumqcb, dumqch, dumqhi, dumqhj, e1, e3, factj, rav, vj(4), vjbh(4), vjcb, vjch(4), vq(4), vqbh(4), vqcb, vqch(4), &
                & vstr
            integer :: i, ind, j, k
            !  calculate avergage bond length for the nh3 moiety
            rav = (rch(1)+rch(2)+rch(3))/3.0e0_wp
            !  initialise:
            vstr = 0.0e0_wp
            !  ach:
            !  in double precision tanh(19.0d0)=1.0d0 and we put the if statement
            !  in to avoid overflow/underflow errors
            arga = c1ch*(rav-r0ch)
            if( arga .lt. 19.0e0_wp )then
                ach    = a1ch+b1ch*(tanh(arga)+1.0e0_wp)*0.5e0_wp
                dumach = b1ch*c1ch/(2.0e0_wp*cosh(arga)**2)
            else
                ach    = a1ch+b1ch
                dumach = 0.0e0_wp
            endif
            !  calculate singlet: e1, triplet: e3 energies and vq and vj
            !  terms for each bond
            e1   = d1cb*(exp(-2.0e0_wp*acb*(rcb-r0cb))-2.0e0_wp*exp(-acb*(rcb-r0cb)))
            e3   = d3cb*(exp(-2.0e0_wp*acb*(rcb-r0cb))+2.0e0_wp*exp(-acb*(rcb-r0cb)))
            vqcb = (e1+e3)*0.5e0_wp
            vjcb = (e1-e3)*0.5e0_wp
            do i=1, 3
                e1      = d1ch*(exp(-2.0e0_wp*ach*(rch(i)-r0ch)) -2.0e0_wp*exp(-ach*(rch(i)-r0ch)))
                e3      = d3ch*(exp(-2.0e0_wp*ach*(rch(i)-r0ch)) +2.0e0_wp*exp(-ach*(rch(i)-r0ch)))
                vqch(i) = (e1+e3)*0.5e0_wp
                vjch(i) = (e1-e3)*0.5e0_wp
                e1      = d1hh*(exp(-2.0e0_wp*ahh*(rbh(i)-r0hh)) -2.0e0_wp*exp(-ahh*(rbh(i)-r0hh)))
                e3      = d3hh*(exp(-2.0e0_wp*ahh*(rbh(i)-r0hh)) +2.0e0_wp*exp(-ahh*(rbh(i)-r0hh)))
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
                achdc(ind) = dumach*(tch(1,ind)/rch(1)+tch(2,ind)/rch(2) +tch(3,ind)/rch(3))/3.0e0_wp
                do i=1, 3
                    achdh(i,ind) = -dumach*tch(i,ind)/rch(i)/3.0e0_wp
                enddo
            enddo
            dumqcb = -acb*((d1cb+d3cb)*exp(-2.0e0_wp*acb*(rcb-r0cb))- (d1cb-d3cb)*exp(-acb*(rcb-r0cb)))/rcb
            !  calculate cartesian derivatives:
            !  looping over ch(i) and bh(i)
            do i=1, 3
                dumqbh = -ahh*((d1hh+d3hh)*exp(-2.0e0_wp*ahh*(rbh(i)-r0hh))- (d1hh-d3hh)*exp(-ahh*(rbh(i)-r0hh)))/rbh(i)
                factj  = 0.5e0_wp/vj(i)
                dumjcb = -acb*((d1cb-d3cb)*exp(-2.0e0_wp*acb*(rcb-r0cb)) -(d1cb+d3cb)*exp(-acb*(rcb-r0cb)))*factj/rcb
                dumjbh = -ahh*((d1hh-d3hh)*exp(-2.0e0_wp*ahh*(rbh(i)-r0hh)) -(d1hh+d3hh)*exp(-ahh*(rbh(i)-r0hh)))* &
                    & factj/ rbh(i)
                do ind=1, 3
                    !  deriv wrt hb:
                    pdot(nhb(ind))=pdot(nhb(ind)) -tcb(ind)*dumqcb+tbh(i,ind)*dumqbh +(vjch(i)-vjcb)*(dumjcb* &
                        & tcb(ind)) + (vjcb- vjbh(i))*(-dumjcb*tcb(ind)-dumjbh*tbh(i,ind)) +(vjbh(i)-vjch(i))*dumjbh* &
                        & tbh(i,ind)
                    !  dvqch(i)/dc
                    dumqch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch)) *((d1ch+d3ch)*exp(-2.0e0_wp*ach*(rch(i)- &
                        & r0ch)) - (d1ch- d3ch)*exp(-ach*(rch(i)-r0ch)))
                    pdot(nc(ind)) = pdot(nc(ind))+dumqch+tcb(ind)*dumqcb
                    !  dvqch(i)/dh(i)
                    dumqhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch)) *((d1ch+d3ch)*exp(-2.0e0_wp*ach*(rch(i)- &
                        & r0ch)) - (d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
                    pdot(nh(i,ind)) = pdot(nh(i,ind))+dumqhi-tbh(i,ind)*dumqbh
                    !  dvjch(i)/dc
                    dumjch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch)) *((d1ch-d3ch)*exp(-2.0e0_wp*ach*(rch(i)- &
                        & r0ch)) - (d1ch+ d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
                    !  dvj(i)/dnc(ind)
                    pdot(nc(ind))=pdot(nc(ind)) +(vjch(i)-vjcb)*(dumjch-dumjcb*tcb(ind)) +(vjcb-vjbh(i))*dumjcb* &
                        & tcb(ind) - (vjbh(i)-vjch(i))*dumjch
                    !  dvjch(i)/dh(i)
                    dumjhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch)) *((d1ch-d3ch)*exp(-2.0e0_wp*ach*(rch(i)- &
                        & r0ch)) - (d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
                    !  dvj(i)/dnh(i,ind)
                    pdot(nh(i,ind))=pdot(nh(i,ind)) +(vjch(i)-vjcb)*dumjhi +(vjcb-vjbh(i))*dumjbh*tbh(i,ind) + &
                        & (vjbh(i)- vjch(i))* (-dumjbh*tbh(i,ind)-dumjhi)
                    !  dv(i)/dh(j)
                    do k=1, 2
                        j = i+k
                        if( j .gt. 3 ) j=j-3
                        dumqhj=-achdh(j,ind)*(rch(i)-r0ch) *((d1ch+d3ch)*exp(-2.0e0_wp*ach*(rch(i)-r0ch)) -(d1ch- &
                            & d3ch)* exp(-ach* (rch(i)-r0ch)))
                        dumjhj=-achdh(j,ind)*(rch(i)-r0ch) *((d1ch-d3ch)*exp(-2.0e0_wp*ach*(rch(i)-r0ch)) -(d1ch+ &
                            & d3ch)* exp(-ach* (rch(i)-r0ch)))*factj
                        pdot(nh(j,ind)) = pdot(nh(j,ind))+dumqhj +(vjch(i)-vjcb)*dumjhj -(vjbh(i)-vjch(i))*dumjhj
                    enddo
                enddo
            enddo
        end subroutine stretch

        !> @brief Evaluate opbend for the NH4_VBMM_2010 potential.
        subroutine opbend( vop )
            !  subroutine calculates symmetrized vop potential and derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3), argd(4), axb(3), b(3), c(4,3), ddr, sum2, sum4, sumd2(4), sumd4(4), vop
            integer  :: i, ii, in(3), ind, itemp, j, k
            real(wp) :: norma
            vop = 0.0e0_wp
            !  calculate force constants and their derivatives
            call opforce
            !  calculate out-of-plane angle and derivatives
            do i=1, 3
                j = i+1
                if( j .gt. 3 ) j=j-3
                k = j+1
                if( k .gt. 3 ) k=k-3
                !  modification to ensure that the set of methane CH bond vector
                !  (rj,rk,rl) is a right-handed set
                in(1) = i
                in(2) = j
                in(3) = k
                !  vector a is rk-rj, vector b is rl-rj
                do ind=1, 3
                    a(ind) = q(nh(j,ind))-q(nh(i,ind))
                    b(ind) = q(nh(k,ind))-q(nh(i,ind))
                enddo
                !  axb is vector a cross b
                axb(1) = a(2)*b(3)-a(3)*b(2)
                axb(2) = a(3)*b(1)-a(1)*b(3)
                axb(3) = a(1)*b(2)-a(2)*b(1)
                norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
                norma  = sqrt(norma)
                !  c is position vector of h(ii): calculate c(j),c(k),c(l)
                do ii=1, 3
                    do ind=1, 3
                        c(in(ii),ind) = -tch(in(ii),ind)/rch(in(ii))
                    enddo
                enddo
                !  argd is the dot product axb dot c
                do ii=1, 3
                    argd(in(ii)) = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                    argd(in(ii)) = argd(in(ii))/norma
                    !  if argd < 0 we need to switch vectors k and l around
                    if( argd(in(ii)) .lt. 0.e0_wp )then
                        itemp = j
                        j     = k
                        k     = itemp
                    endif
                enddo
                !  sum2 = sum delta**2
                !  sum4 = sum delta**4
                call calcdelta(i, j, k, sum2, sum4)
                sumd2(i) = sum2
                sumd4(i) = sum4
                vop      = vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
            enddo
            do i=1, 3
                do j=1, 3
                    !  overall derivatives of force constants i wrt the bond-length rch(j)
                    ddr = dfdelta(i,j)*sumd2(i)+dhdelta(i,j)*sumd4(i)
                    !  calculate derivatives in terms of cartesian coordinates:
                    do ind=1, 3
                        pdot(nh(j,ind)) = pdot(nh(j,ind))-tch(j,ind)*ddr/rch(j)
                        pdot(nc(ind))   = pdot(nc(ind))+tch(j,ind)*ddr/rch(j)
                    enddo
                enddo
            enddo
        end subroutine opbend

        !> @brief Evaluate ipbend for the NH4_VBMM_2010 potential.
        subroutine ipbend( vip )
            !  subroutine calculates symmetrised in plane bend term
            !  and its derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: costh(4,4), dth(4,4), dth0c, dth0k, dthc, dthi, dthj, termth, theta(4,4), vip
            integer  :: i, ind, j, k
            !  initialise
            vip = 0.0e0_wp
            !  calculate force constants: fk0(i,j), f1(i)
            !  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
            call ipforce
            !  calculate theta(i,j) and in plane bend potential
            do i=1, 2
                do j=i+1, 3
                    costh(i,j) = tch(i,1)*tch(j,1)+tch(i,2)*tch(j,2) +tch(i,3)*tch(j,3)
                    costh(i,j) = costh(i,j)/rch(i)/rch(j)
                    theta(i,j) = acos(costh(i,j))
                    dth(i,j)   = theta(i,j)-theta0(i,j)
                    vip        = vip+0.5e0_wp*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
                    !  calculate partial derivatives wrt cartesian coordinates
                    !  calculate pdots wrt theta:
                    termth = -1.0e0_wp/sqrt(1.0e0_wp-costh(i,j)*costh(i,j))
                    do ind=1, 3
                        dthi            = -tch(j,ind)/rch(i)/rch(j) +costh(i,j)*tch(i,ind)/rch(i)/rch(i)
                        dthi            = dthi*termth
                        dthj            = -tch(i,ind)/rch(i)/rch(j) +costh(i,j)*tch(j,ind)/rch(j)/rch(j)
                        dthj            = dthj*termth
                        dthc            = -(dthi+dthj)
                        pdot(nh(i,ind)) = pdot(nh(i,ind)) +fk0(i,j)*f1(i)*f1(j)*dthi*dth(i,j)
                        pdot(nh(j,ind)) = pdot(nh(j,ind)) +fk0(i,j)*f1(i)*f1(j)*dthj*dth(i,j)
                        pdot(nc(ind))   = pdot(nc(ind)) +fk0(i,j)*f1(i)*f1(j)*dthc*dth(i,j)
                        do k=1, 3
                            !  calculate pdots wrt force constants and wrt theta0
                            dth0k = -dtheta0(i,j,k)*tch(k,ind)/rch(k)
                            dth0c = -dth0k
                            pdot(nh(k,ind))=pdot(nh(k,ind)) -0.5e0_wp*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k) - &
                                & 0.5e0_wp* tbh(k,ind)* dfdh(i,j,k)*dth(i,j)**2/rbh(k) -fk0(i,j)*f1(i)*f1(j)*dth0k* &
                                & dth(i,j)
                            pdot(nc(ind))=pdot(nc(ind)) +0.5e0_wp*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k) -fk0(i,j)* &
                                & f1(i)* f1(j)* dth0c*dth(i,j)
                            pdot(nhb(ind)) = pdot(nhb(ind)) +0.5e0_wp*tbh(k,ind)*dfdh(i,j,k)*dth(i,j)**2/rbh(k)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine ipbend

        !> @brief Evaluate calcdelta for the NH4_VBMM_2010 potential.
        subroutine calcdelta( i, j, k, sum2, sum4 )
            !  subroutine calculates out of plane angle delta, loops
            !  through delta(i,j), delta(i,k), delta(i,l)
            !   also calculates the derivatives wrt delta
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a(3), argd(4), atemp1, atemp2(3), atemp3, atemp4, atemp5, axb(3), b(3), c(4,3), &
                & cdot(4, 3, 3), daxb(4, 3, 3), deldot, delta(4), sum2, sum4
            integer  :: i, ii, in(3), ind, j, jj, k
            real(wp) :: norma
            !  set i,j,k indices
            in(1) = i
            in(2) = j
            in(3) = k
            !  vector a is rj-ri, vector b is rk-ri
            do ind=1, 3
                a(ind) = q(nh(j,ind))-q(nh(i,ind))
                b(ind) = q(nh(k,ind))-q(nh(i,ind))
            enddo
            !  axb is vector a cross b
            axb(1) = a(2)*b(3)-a(3)*b(2)
            axb(2) = a(3)*b(1)-a(1)*b(3)
            axb(3) = a(1)*b(2)-a(2)*b(1)
            norma  = axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
            norma  = sqrt(norma)
            !  c is position vector of h(ii): calculate c(j),c(k),c(l)
            do ii=1, 3
                do ind=1, 3
                    c(in(ii),ind) = -tch(in(ii),ind)/rch(in(ii))
                enddo
            enddo
            !  initialise
            sum2 = 0.0e0_wp; sum4 = 0.0e0_wp
            !  argd is the dot product axb dot c
            do ii=1, 3
                argd(in(ii))  = axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2) +axb(3)*c(in(ii),3)
                argd(in(ii))  = argd(in(ii))/norma
                delta(in(ii)) = acos(argd(in(ii)))-optheta0(i,in(ii))
                sum2          = sum2+delta(in(ii))**2
                sum4          = sum4+delta(in(ii))**4
            enddo
            !  derivatives of axb wrt hi:
            daxb(i, 1, 1) = 0.0e0_wp
            daxb(i, 1, 2) = b(3)-a(3)
            daxb(i, 1, 3) = -b(2)+a(2)
            daxb(i, 2, 1) = -b(3)+a(3)
            daxb(i, 2, 2) = 0.0e0_wp
            daxb(i, 2, 3) = b(1)-a(1)
            daxb(i, 3, 1) = b(2)-a(2)
            daxb(i, 3, 2) = -b(1)+a(1)
            daxb(i, 3, 3) = 0.0e0_wp
            !  derivatives of axb wrt hj:
            daxb(j, 1, 1) = 0.0e0_wp
            daxb(j, 1, 2) = -b(3)
            daxb(j, 1, 3) = b(2)
            daxb(j, 2, 1) = b(3)
            daxb(j, 2, 2) = 0.0e0_wp
            daxb(j, 2, 3) = -b(1)
            daxb(j, 3, 1) = -b(2)
            daxb(j, 3, 2) = b(1)
            daxb(j, 3, 3) = 0.0e0_wp
            !  derivatives of axb wrt hl:
            daxb(k, 1, 1) = 0.0e0_wp
            daxb(k, 1, 2) = a(3)
            daxb(k, 1, 3) = -a(2)
            daxb(k, 2, 1) = -a(3)
            daxb(k, 2, 2) = 0.0e0_wp
            daxb(k, 2, 3) = a(1)
            daxb(k, 3, 1) = a(2)
            daxb(k, 3, 2) = -a(1)
            daxb(k, 3, 3) = 0.0e0_wp
            !   loop over cdot(in(ii),ind,jind) where we consider deriv of c(in(ii))
            !   wrt h(in(ii),jind) with components jind
            do ii=1, 3
                !  deriv of cdot(in(ii),x) wrt x, y, z
                cdot(in(ii), 1, 1) = 1.0e0_wp/rch(in(ii)) +tch(in(ii),1)*c(in(ii),1)/rch(in(ii))**2
                cdot(in(ii), 1, 2) = tch(in(ii),2)*c(in(ii),1)/rch(in(ii))**2
                cdot(in(ii), 1, 3) = tch(in(ii),3)*c(in(ii),1)/rch(in(ii))**2
                !  deriv of cdot(in(ii),y) wrt x, y, z
                cdot(in(ii), 2, 1) = tch(in(ii),1)*c(in(ii),2)/rch(in(ii))**2
                cdot(in(ii), 2, 2) = 1.0e0_wp/rch(in(ii)) +tch(in(ii),2)*c(in(ii),2)/rch(in(ii))**2
                cdot(in(ii), 2, 3) = tch(in(ii),3)*c(in(ii),2)/rch(in(ii))**2
                !  deriv of cdot(in(ii),z) wrt x, y, z
                cdot(in(ii), 3, 1) = tch(in(ii),1)*c(in(ii),3)/rch(in(ii))**2
                cdot(in(ii), 3, 2) = tch(in(ii),2)*c(in(ii),3)/rch(in(ii))**2
                cdot(in(ii), 3, 3) = 1.0e0_wp/rch(in(ii)) +tch(in(ii),3)*c(in(ii),3)/rch(in(ii))**2
            enddo
            do ii=1, 3
                do ind=1, 3
                    deldot = -doptheta0(i,ii,i)
                    !  derivative wrt h(i,ind)
                    !  for  rch(i) only terms are from the derivatives of theta0
                    deldot = -deldot*tch(i,ind)/rch(i)
                    !  derivative wrt c(ind)
                    deldot = -deldot
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
                        atemp5 = -doptheta0(i,in(ii),in(jj))
                        !  deriv wrt h(in(jj)),ind):
                        atemp5 = -atemp5*tch(in(jj),ind)/rch(in(jj))
                        deldot = atemp3+atemp4
                        deldot = -1.0e0_wp/sqrt(1.0e0_wp-argd(in(ii))**2)*deldot
                        deldot = deldot+atemp5
                        pdot(nh(in(jj),ind))=pdot(nh(in(jj),ind)) +2.0e0_wp*fdelta(i)*delta(in(ii))*deldot +4.0e0_wp* &
                            & hdelta(i)* delta(in(ii))**3*deldot
                        !  for carbon the only contributions are from axb dot cdot term and
                        !  from theta0 and derivative cdot wrt carbon=-cdot wrt hydrogen
                        deldot = 1.0e0_wp/sqrt(1.0e0_wp-argd(in(ii))**2)*atemp4
                        deldot = deldot-atemp5
                        pdot(nc(ind))=pdot(nc(ind)) +2.0e0_wp*fdelta(i)*delta(in(ii))*deldot +4.0e0_wp*hdelta(i)* &
                            & delta(in(ii))**3* deldot
                    enddo
                enddo
            enddo
        end subroutine calcdelta

        !> @brief Evaluate opforce for the NH4_VBMM_2010 potential.
        subroutine opforce
            !  calculates the out-of-plane bending force constants
            !  and their derivatives
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: dswitch(4,4), switch(4)
            integer  :: i, j
            !  calculate switching functions:
            switch(1) = s3(1)*s3(2)*s3(3)
            switch(2) = s3(2)*s3(3)*s3(1)
            switch(3) = s3(3)*s3(1)*s3(2)
            !  calculate derivatives:
            !  derivative of switch(1) wrt the 3 rch bond lengths
            dswitch(1,1) = ds3(1)*s3(2)*s3(3)
            dswitch(1,2) = s3(1)*ds3(2)*s3(3)
            dswitch(1,3) = s3(1)*s3(2)*ds3(3)
            !  derivative of switch(2) wrt the 3 rch bond lengths
            dswitch(2,1) = s3(2)*s3(3)*ds3(1)
            dswitch(2,2) = ds3(2)*s3(3)*s3(1)
            dswitch(2,3) = s3(2)*ds3(3)*s3(1)
            !  derivative of switch(3) wrt the 3 rch bond lengths
            dswitch(3,1) = s3(3)*ds3(1)*s3(2)
            dswitch(3,2) = s3(3)*s3(1)*ds3(2)
            dswitch(3,3) = ds3(3)*s3(1)*s3(2)
            !  calculate the force constants and their derivatives
            do i=1, 3
                fdelta(i) = switch(i)*fch3
                hdelta(i) = switch(i)*hch3
                do j=1, 3
                    dfdelta(i,j) = dswitch(i,j)*fch3
                    dhdelta(i,j) = dswitch(i,j)*hch3
                enddo
            enddo
        end subroutine opforce

        !> @brief Evaluate ipforce for the NH4_VBMM_2010 potential.
        subroutine ipforce
            !  calculates the symmetrised in plane bend force constants and
            !  all partial derivatives involving them
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: a1, a2, arga1, arga2, argfk0, df1dc(4), df1dh(4), dfko(4), duma1, duma2, f0, f2, fko
            integer  :: i
            !  set force constant at asymptotes
            f0       = fkinf+ak
            f2       = fkinf
            argfk0   = bk*((rch(1)-r0ch)**2+(rch(2)-r0ch)**2 +(rch(3)-r0ch)**2)
            fko      = fkinf+ak*exp(-argfk0)
            fk0(1,2) = fko
            fk0(1,3) = fko
            fk0(1,4) = fko
            fk0(2,3) = fko
            fk0(2,4) = fko
            fk0(3,4) = fko
            do i=1, 3
                !  calc derivatives of fko wrt each of the rch(i) bonds
                dfko(i) = -2.0e0_wp*ak*bk*(rch(i)-r0ch)*exp(-argfk0)
                !  calculate the terms f1(i)
                arga1 = aa1*rbh(i)*rbh(i)
                arga2 = aa4*(rbh(i)-r0hh)*(rbh(i)-r0hh)
                a1    = 1.0e0_wp-exp(-arga1)
                a2    = aa2+aa3*exp(-arga2)
                f1(i) = a1*exp(-a2*(rch(i)-r0ch)**2)
                !  and calculate the derivatives wrt rch(i) and rbh(i)
                duma1    = 2.0e0_wp*aa1*rbh(i)*exp(-arga1)
                duma2    = -2.0e0_wp*aa3*aa4*(rbh(i)-r0hh)*exp(-arga2)
                df1dc(i) = -2.0e0_wp*(rch(i)-r0ch)*a1*a2*exp(-a2*(rch(i)-r0ch)**2)
                df1dh(i) = duma1*exp(-a2*(rch(i)-r0ch)**2) -duma2*(rch(i)-r0ch)**2*a1*exp(-a2*(rch(i)-r0ch)**2)
            enddo
            !  derivative of total force constant f(i,j) wrt bond length rch(k)
            !  is given by dfdc(i,j,k)
            dfdc(1,2,1) = dfko(1)*f1(1)*f1(2)+fko*df1dc(1)*f1(2)
            dfdc(1,2,2) = dfko(2)*f1(1)*f1(2)+fko*f1(1)*df1dc(2)
            dfdc(1,2,3) = dfko(3)*f1(1)*f1(2)
            dfdc(1,3,1) = dfko(1)*f1(1)*f1(3)+fko*df1dc(1)*f1(3)
            dfdc(1,3,2) = dfko(2)*f1(1)*f1(3)
            dfdc(1,3,3) = dfko(3)*f1(1)*f1(3)+fko*f1(1)*df1dc(3)
            dfdc(2,3,1) = dfko(1)*f1(2)*f1(3)
            dfdc(2,3,2) = dfko(2)*f1(2)*f1(3)+fko*df1dc(2)*f1(3)
            dfdc(2,3,3) = dfko(3)*f1(2)*f1(3)+fko*f1(2)*df1dc(3)
            !  derivative of total force constant f(i,j) wrt bond length rbh(k)
            !  is given by dfdh(i,j,k)
            !  only non-zero derivatives are those from rbh(i) and rbh(j)
            dfdh(1,2,1) = fko*df1dh(1)*f1(2)
            dfdh(1,2,2) = fko*f1(1)*df1dh(2)
            dfdh(1,2,3) = 0.0e0_wp
            dfdh(1,3,1) = fko*df1dh(1)*f1(3)
            dfdh(1,3,2) = 0.0e0_wp
            dfdh(1,3,3) = fko*f1(1)*df1dh(3)
            dfdh(2,3,1) = 0.0e0_wp
            dfdh(2,3,2) = fko*df1dh(2)*f1(3)
            dfdh(2,3,3) = fko*f1(2)*df1dh(3)
        end subroutine ipforce

        !> @brief Evaluate switchf for the NH4_VBMM_2010 potential.
        subroutine switchf
            !  calculates switching functions: s3,sphi,stheta
            !  and their derivatives ds3,dsphi,dstheta
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: argmax, args1, args2, args3, argsphi, argstheta
            integer  :: i
            a1s = 1.5313681e-7_wp
            b1s = -4.6696246e0_wp
            a2s = 1.0147402e-7_wp
            b2s = -12.363798e0_wp
            !  use double precision criterion:
            !  tanh(19.0d0)=1.0d0
            argmax = 19.0e0_wp
            do i=1, 3
                args1 = a1s*(rch(i)-r0ch)*(rch(i)-b1s)**8
                if( args1 .lt. argmax )then
                    s1(i)  = 1.0e0_wp-tanh(args1)
                    ds1(i) = a1s*((rch(i)-b1s)**8 +8.0e0_wp*(rch(i)-r0ch)*(rch(i)-b1s)**7)
                    ds1(i) = -ds1(i)/cosh(args1)**2
                else
                    s1(i)  = 0.0e0_wp
                    ds1(i) = 0.0e0_wp
                endif
                args2 = a2s*(rch(i)-r0ch)*(rch(i)-b2s)**6
                if( args2 .lt. argmax )then
                    s2(i)  = 1.0e0_wp-tanh(args2)
                    ds2(i) = a2s*((rch(i)-b2s)**6 +6.0e0_wp*(rch(i)-r0ch)*(rch(i)-b2s)**5)
                    ds2(i) = -ds2(i)/cosh(args2)**2
                else
                    s2(i)  = 0.0e0_wp
                    ds2(i) = 0.0e0_wp
                endif
                !  calculate s3 and ds3
                args3 = a3s*(rch(i)-r0ch)*(rch(i)-b3s)**2
                if( args3 .lt. argmax )then
                    s3(i)  = 1.0e0_wp-tanh(args3)
                    ds3(i) = a3s*(3.0e0_wp*rch(i)**2-2.0e0_wp*rch(i)*(r0ch+2.0e0_wp*b3s) +b3s*(b3s+2.0e0_wp*r0ch))
                    ds3(i) = -ds3(i)/cosh(args3)**2
                else
                    s3(i)  = 0.0e0_wp
                    ds3(i) = 0.0e0_wp
                endif
                !  calculate sphi and dsphi
                !  condition here is on the bondlength rch(i)
                !  st argsphi is lt approx 19.0d0
                if( rch(i) .lt. 3.8e0_wp )then
                    argsphi  = aphi*(rch(i)-r0ch)*exp(bphi*(rch(i)-cphi)**3)
                    sphi(i)  = 1.0e0_wp-tanh(argsphi)
                    dsphi(i) = aphi*(1.0e0_wp+3.0e0_wp*bphi*(rch(i)-r0ch) *(rch(i)-cphi)**2)
                    dsphi(i) = dsphi(i)*exp(bphi*(rch(i)-cphi)**3)
                    dsphi(i) = -dsphi(i)/cosh(argsphi)**2
                else
                    sphi(i)  = 0.0e0_wp
                    dsphi(i) = 0.0e0_wp
                endif
                !  calculate stheta and dstheta
                if( rch(i) .lt. 3.8e0_wp )then
                    argstheta  = atheta*(rch(i)-r0ch)*exp(btheta*(rch(i)-ctheta)**3)
                    stheta(i)  = 1.0e0_wp-tanh(argstheta)
                    dstheta(i) = atheta*(1.0e0_wp+3.0e0_wp*btheta*(rch(i)-r0ch) *(rch(i)-ctheta)**2)
                    dstheta(i) = dstheta(i)*exp(btheta*(rch(i)-ctheta)**3)
                    dstheta(i) = -dstheta(i)/cosh(argstheta)**2
                else
                    stheta(i)  = 0.0e0_wp
                    dstheta(i) = 0.0e0_wp
                endif
            enddo
        end subroutine switchf

        !> @brief Evaluate prepare_local_context for the NH4_VBMM_2010 potential.
        subroutine prepare_local_context
            !   N3TMMN = 3 * atom_count
            !   atom_count = the number of atoms represented by this potential function
            !   The variable N3TMMN is the minimum value of N3TM allowed to be
            !   passed by the calling routine for the number of cartesian
            !   coordinates needed to represent the full system represented by this
            !   potential energy surface routine.
            !   N3TM must be greater than or equal to N3TMMN.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: fact1, fact2
            integer  :: i, icount, ind
            !  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            do i=1, 5
                ref(i) = ' '
            enddo
            ref(1) = 'J. Espinosa-Garcia'
            ref(2) = 'unpublished results, 2008'
            !  calculate indexes for coordinates
            do ind=1, 3
                icount   = ind-3
                nc(ind)  = 3*nnc+icount
                nhb(ind) = 3*nnb+icount
                do i=1, 3
                    nh(i,ind) = 3*nnh(i)+icount
                enddo
            enddo
            !  convert to appropriate units:
            !  energy   in 1.0d+05 j/mol
            !  time     in 1.0d-14 s
            !  distance in 1.0d-10 m
            !  angles   in radians
            !  mass     in amu
            fact1 = 0.041840e0_wp
            fact2 = 6.022045e0_wp
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
        end subroutine prepare_local_context

        !> @brief Evaluate CARTTOR for the NH4_VBMM_2010 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(cartesian_coordinates, atom_count, coordinate_representation, r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the NH4_VBMM_2010 potential.
        subroutine eunitzero
            implicit none
            integer :: i, j
            output_ground_state_energy = ground_state_energy
            if( first_requested_excited_state .ne. 0 )then
                do i=first_requested_excited_state, last_requested_excited_state
                    output_excited_state_energies(i) = excited_state_energies(i)
                enddo
            endif
            if( first_requested_coupling .ne. 0 )then
                do j=first_requested_coupling, last_requested_coupling
                    output_diabatic_couplings(j) = diabatic_couplings(j)
                enddo
            endif
        end subroutine eunitzero

        !> @brief Evaluate RTOCART for the NH4_VBMM_2010 potential.
        subroutine rtocart
            use chempotpy_kinds, only: wp
            implicit none
            integer, parameter :: n3atom = 75
            integer, parameter :: isurf  = 5
            integer, parameter :: jsurf  = isurf*(isurf+1)/2
            real(wp)           :: eps, epsp, phi, r1, r2, rcm, rff, ss1, ss2, ss3, sum, t0, t1, term12, term13, term23, termx, &
                & termy, termz, theta1, theta2, u1, u2, u3, wf, wh, xa, xbb, xc, xd, y, ya, yb, yc, yd, yes(n3atom, isurf), &
                & ygs(n3atom), yij(n3atom, jsurf), za, zb, zc, zd
            integer :: i, ind, j, j1, j2, k, m1
            if( coordinate_representation .eq. 1 )then
                do i = 1, atom_count
                    ind                                  = 3*i-2
                    ground_state_cartesian_gradient(i,1) = ground_state_internal_gradient(ind)
                    ground_state_cartesian_gradient(i,2) = ground_state_internal_gradient(ind+1)
                    ground_state_cartesian_gradient(i,3) = ground_state_internal_gradient(ind+2)
                    if( first_requested_excited_state .ne. 0 )then
                        do j = first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(i,1,j) = excited_state_internal_gradients(ind,j)
                            excited_state_cartesian_gradients(i,2,j) = excited_state_internal_gradients(ind+1,j)
                            excited_state_cartesian_gradients(i,3,j) = excited_state_internal_gradients(ind+2,j)
                        enddo
                    endif
                    if( last_requested_coupling .ne. 0 )then
                        do k = first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(i,1,k) = diabatic_internal_gradients(ind,k)
                            diabatic_cartesian_gradients(i,2,k) = diabatic_internal_gradients(ind+1,k)
                            diabatic_cartesian_gradients(i,3,k) = diabatic_internal_gradients(ind+2,k)
                        enddo
                    endif
                enddo
            elseif( coordinate_representation .eq. 2 )then
                do i = 1, atom_count
                    ground_state_cartesian_gradient(i,1) = 0.e0_wp
                    ground_state_cartesian_gradient(i,2) = 0.e0_wp
                    ground_state_cartesian_gradient(i,3) = 0.e0_wp
                    if( first_requested_excited_state .ne. 0 )then
                        do j1=first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(i,1,j1) = 0.e0_wp
                            excited_state_cartesian_gradients(i,2,j1) = 0.e0_wp
                            excited_state_cartesian_gradients(i,3,j1) = 0.e0_wp
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do j2=first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(i,1,j2) = 0.e0_wp
                            diabatic_cartesian_gradients(i,2,j2) = 0.e0_wp
                            diabatic_cartesian_gradients(i,3,j2) = 0.e0_wp
                        enddo
                    endif
                    do j = 1, atom_count
                        if( j .lt. i )then
                            m1 = atom_count*(j-1) - (j*(j-1))/2 + i-j
                        elseif( j .gt. i )then
                            m1 = atom_count*(i-1) - (i*(i-1))/2 + j-i
                        else
                            cycle
                        endif
                        y = ground_state_internal_gradient(m1)
                        termx                                = (cartesian_coordinates(i,1)-cartesian_coordinates(j, &
                            & 1))/r(m1)
                        termy                                = (cartesian_coordinates(i,2)-cartesian_coordinates(j, &
                            & 2))/r(m1)
                        termz                                = (cartesian_coordinates(i,3)-cartesian_coordinates(j, &
                            & 3))/r(m1)
                        ground_state_cartesian_gradient(i,1) = ground_state_cartesian_gradient(i,1) + termx*y
                        ground_state_cartesian_gradient(i,2) = ground_state_cartesian_gradient(i,2) + termy*y
                        ground_state_cartesian_gradient(i,3) = ground_state_cartesian_gradient(i,3) + termz*y
                        if( first_requested_excited_state .gt. 0 )then
                            y = excited_state_internal_gradients(m1,j1)
                            do j1=first_requested_excited_state, last_requested_excited_state
                                excited_state_cartesian_gradients(i,1,j1) = excited_state_cartesian_gradients(i,1, &
                                    & j1) + termx*y
                                excited_state_cartesian_gradients(i,2,j1) = excited_state_cartesian_gradients(i,2, &
                                    & j1) + termy*y
                                excited_state_cartesian_gradients(i,3,j1) = excited_state_cartesian_gradients(i,3, &
                                    & j1) + termz*y
                            enddo
                        elseif( first_requested_coupling .gt. 0 )then
                            do j2=first_requested_coupling, last_requested_coupling
                                y                                    = diabatic_internal_gradients(m1,j2)
                                diabatic_cartesian_gradients(i,1,j2) = diabatic_cartesian_gradients(i,1,j2) + termx*y
                                diabatic_cartesian_gradients(i,2,j2) = diabatic_cartesian_gradients(i,2,j2) + termy*y
                                diabatic_cartesian_gradients(i,3,j2) = diabatic_cartesian_gradients(i,3,j2) + termz*y
                            enddo
                        endif
                    enddo
                enddo
            elseif( coordinate_representation .eq. 3 )then
                do i = 1, atom_count
                    ygs(i) = ground_state_internal_gradient(i)/r(i)
                    if( first_requested_excited_state .ne. 0 )then
                        do j=first_requested_excited_state, last_requested_excited_state
                            yes(i, j) = excited_state_internal_gradients(i,j)/r(i)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do k=first_requested_coupling, last_requested_coupling
                            yij(i, k) = diabatic_internal_gradients(i,k)/r(i)
                        enddo
                    endif
                enddo
                do k = 1, 3
                    term12                               = cartesian_coordinates(1,k)-cartesian_coordinates(2,k)
                    term23                               = cartesian_coordinates(2,k)-cartesian_coordinates(3,k)
                    term13                               = cartesian_coordinates(1,k)-cartesian_coordinates(3,k)
                    ground_state_cartesian_gradient(1,k) = term12*ygs(1) + term13*ygs(3)
                    ground_state_cartesian_gradient(2,k) = -term12*ygs(1) + term23*ygs(2)
                    ground_state_cartesian_gradient(3,k) = -term13*ygs(3) - term23*ygs(2)
                    if( first_requested_excited_state .ne. 0 )then
                        do j1=first_requested_excited_state, last_requested_excited_state
                            excited_state_cartesian_gradients(1,k,j1) = term12*yes(1, j1) + term13*yes(3, j1)
                            excited_state_cartesian_gradients(2,k,j1) = -term12*yes(1, j1) + term23*yes(2, j1)
                            excited_state_cartesian_gradients(3,k,j1) = -term13*yes(3, j1) - term23*yes(2, j1)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do j2=first_requested_coupling, last_requested_coupling
                            diabatic_cartesian_gradients(1,k,j2) = term12*yij(1, j2) + term13*yij(3, j2)
                            diabatic_cartesian_gradients(2,k,j2) = -term12*yij(1, j2) + term23*yij(2, j2)
                            diabatic_cartesian_gradients(3,k,j2) = -term13*yij(3, j2) - term23*yij(2, j2)
                        enddo
                    endif
                enddo
            elseif( coordinate_representation .eq. 4 )then
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
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine rtocart

        !> @brief Evaluate DEDCOU for the NH4_VBMM_2010 potential.
        subroutine dedcou
            implicit none
            integer :: i, j, k
            do i = 1, atom_count
                output_ground_state_gradient(i,1) = ground_state_cartesian_gradient(i,1)
                output_ground_state_gradient(i,2) = ground_state_cartesian_gradient(i,2)
                output_ground_state_gradient(i,3) = ground_state_cartesian_gradient(i,3)
                if( first_requested_excited_state .ne. 0 )then
                    do j=first_requested_excited_state, last_requested_excited_state
                        output_excited_state_gradients(i,1,j) = excited_state_cartesian_gradients(i,1,j)
                        output_excited_state_gradients(i,2,j) = excited_state_cartesian_gradients(i,2,j)
                        output_excited_state_gradients(i,3,j) = excited_state_cartesian_gradients(i,3,j)
                    enddo
                endif
                if( first_requested_coupling .ne. 0 )then
                    do k=first_requested_coupling, last_requested_coupling
                        output_diabatic_gradients(i,1,k) = diabatic_cartesian_gradients(i,1,k)
                        output_diabatic_gradients(i,2,k) = diabatic_cartesian_gradients(i,2,k)
                        output_diabatic_gradients(i,3,k) = diabatic_cartesian_gradients(i,3,k)
                    enddo
                endif
            enddo
        end subroutine dedcou
    end subroutine evaluate
end module chempotpy_nh4_vbmm_2010_pes

!> @brief Expose NH4_VBMM_2010 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_nh4_vbmm_2010_pes, only: evaluate
    use chempotpy_kinds,             only: wp
    implicit none
    integer, parameter    :: natoms  = 5
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
