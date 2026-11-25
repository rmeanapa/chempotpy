!   System:          FH2
!   Functional form: eLEPS (extended-London-Eyring-Polanyi-Sato)
!   Common name:     no. 2a
!   Reference:       unpublished
!   Cross reference: D. G. Truhlar, B. C. Garrett, and N. C. Blais
!                    J. Chem. Phys. 80, 232 (1984).
! Number of bodies: 3
! Interface: potlib2001
! Number of electronic states: 1
! Number of derivatives: 1
!   Notes: In this potential routine the FH2 surface no. 2 has been modi
!          so as to make the functional form continuous.
!          In this version of the FH2 potential the AB and AC
!          Sato parameters are functions of PHI, where PHI is
!          the angle of A-to-BC vector with the BC axis.
!   prepare_local_context must be called once before any calls to POT.
!   The potential parameters are included in DATA statements and in
!   the block data subprogram PTPACM.
!   Coordinates, potential energy, and derivatives are passed
!   The potential energy in the three asymptotic valleys are
!   stored in the common block ASYCM:
!   The potential energy in the AB valley, EASYAB, is equal to the poten
!   energy of the H "infinitely" far from the FH diatomic, with the
!   FH diatomic at its equilibrium configuration.  Similarly, the terms
!   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
!   respectively.
!   All the information passed through the common blocks PT1CM and ASYCM
!   is in Hartree atomic units.
!   This potential is written such that:
!                  R(1) = R(F-H)
!                  R(2) = R(H-H)
!                  R(3) = R(H-F)
!   The zero of energy is defined at F "infinitely" far from the H2 diat
!   The flags that indicate what calculations should be carried out in
!   the potential routine are passed through the common block PT2CM:
!   where:
!        surface_requests - which electronic states availalble
!                 (1,1) = 1 as only gs state available
!        requested_derivative_order  = 0 => no derivatives should be calculated
!        requested_derivative_order  = 1 => calculate first derivatives
!        NFLAG - these integer values can be used to flag options
!                within the potential.
!   Potential parameters' default settings
!                  Variable            Default value
!                  requested_derivative_order                1
!                  NFLAG(18)           6

!> @brief Implement FH2_LEPS_NO2a_1984, a potential energy surface for FH2.
module chempotpy_fh2_leps_no2a_1984_pes
    use chempotpy_units,     only: hartree_to_kilocalorie_per_mole
    use chempotpy_constants, only: PI
    use chempotpy_kinds,     only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the FH2_LEPS_NO2a_1984 potential energy surface for FH2.
    subroutine evaluate( x, igrad, p, g, d, status )
        use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input
        use chempotpy_units,  only: hartree_per_bohr_to_ev_per_angstrom
        use chempotpy_units,  only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds,  only: wp
        implicit none
        type(chempotpy_status_t), intent(out), optional :: status
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, parameter    :: natoms  = 3
        integer, parameter    :: natom   = 25
        integer, parameter    :: isurf   = 5
        integer, parameter    :: jsurf   = int(isurf*(isurf+1)/2)
        integer, parameter    :: n3atom  = 75
        integer, intent(in)   :: igrad
        real(wp), intent(in)  :: x(natoms,3)
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: asycm_easyab
        real(wp)              :: asycm_easyac
        real(wp)              :: asycm_easybc
        real(wp)              :: coscm_cosphi
        real(wp)              :: coscm_dphi(3)
        real(wp)              :: cartesian_coordinates(natom,3)
        integer               :: coordinate_representation
        integer               :: maximum_derivative_order
        integer               :: maximum_excited_state
        integer               :: atom_count
        character(len=75)     :: infocm_ref(5)
        real(wp)              :: phicm_deg
        real(wp)              :: phicm_deg1
        real(wp)              :: phicm_fder
        real(wp)              :: phicm_z0
        real(wp)              :: precm_b(3)
        real(wp)              :: precm_beta(3)
        real(wp)              :: precm_coul(3)
        real(wp)              :: precm_d(3)
        real(wp)              :: precm_dcoul(3,3)
        real(wp)              :: precm_dexch(3,3)
        real(wp)              :: precm_do4z(3)
        real(wp)              :: precm_exch(3)
        real(wp)              :: precm_op3z(3)
        real(wp)              :: precm_re(3)
        real(wp)              :: precm_top3z(3)
        real(wp)              :: precm_tzp3(3)
        real(wp)              :: precm_x(3)
        real(wp)              :: precm_zp3(3)
        real(wp)              :: precm_zpo(3)
        real(wp)              :: ground_state_internal_gradient(n3atom)
        real(wp)              :: ground_state_energy
        real(wp)              :: pt1cm_r(n3atom)
        real(wp)              :: pt3cm_ezero(isurf+1)
        real(wp)              :: excited_state_internal_gradients(n3atom,isurf)
        real(wp)              :: excited_state_energies(isurf)
        real(wp)              :: diabatic_internal_gradients(n3atom,jsurf)
        real(wp)              :: diabatic_couplings(jsurf)
        real(wp)              :: input_cartesian(natom,3)
        integer               :: surface_requests(isurf+1,isurf+1)
        integer               :: requested_derivative_order
        real(wp)              :: output_excited_state_gradients(natom,3,isurf)
        real(wp)              :: output_ground_state_gradient(natom,3)
        real(wp)              :: output_diabatic_gradients(natom,3,jsurf)
        real(wp)              :: output_excited_state_energies(isurf)
        real(wp)              :: output_ground_state_energy
        real(wp)              :: output_diabatic_couplings(jsurf)
        real(wp)              :: excited_state_cartesian_gradients(natom,3,isurf)
        real(wp)              :: ground_state_cartesian_gradient(natom,3)
        real(wp)              :: diabatic_cartesian_gradients(natom,3,jsurf)
        integer               :: last_requested_excited_state
        integer               :: last_requested_coupling
        integer               :: first_requested_excited_state
        integer               :: first_requested_coupling
        real(wp)              :: z3cm_z(3)
        real(wp)              :: v
        integer               :: iatom, idir, istate
        !initialize
        asycm_easyab                      = 0.e0_wp; asycm_easyac = 0.e0_wp; asycm_easybc = 0.e0_wp; coscm_cosphi = 0.e0_wp; coscm_dphi = 0.e0_wp
        cartesian_coordinates             = 0.e0_wp; coordinate_representation = 0; maximum_derivative_order = 0
        maximum_excited_state             = 0; atom_count = 0
        infocm_ref                        = ' '
        phicm_deg                         = 0.e0_wp
        phicm_deg1                        = 0.e0_wp
        phicm_fder                        = 0.e0_wp
        phicm_z0                          = 0.e0_wp
        precm_b                           = 0.e0_wp
        precm_beta                        = 0.e0_wp
        precm_coul                        = 0.e0_wp
        precm_d                           = 0.e0_wp
        precm_dcoul                       = 0.e0_wp
        precm_dexch                       = 0.e0_wp
        precm_do4z                        = 0.e0_wp
        precm_exch                        = 0.e0_wp
        precm_op3z                        = 0.e0_wp
        precm_re                          = 0.e0_wp
        precm_top3z                       = 0.e0_wp
        precm_tzp3                        = 0.e0_wp
        precm_x                           = 0.e0_wp
        precm_zp3                         = 0.e0_wp
        precm_zpo                         = 0.e0_wp
        ground_state_internal_gradient    = 0.e0_wp
        ground_state_energy               = 0.e0_wp
        pt1cm_r                           = 0.e0_wp
        pt3cm_ezero                       = 0.e0_wp
        excited_state_internal_gradients  = 0.e0_wp
        excited_state_energies            = 0.e0_wp
        diabatic_internal_gradients       = 0.e0_wp
        diabatic_couplings                = 0.e0_wp
        input_cartesian                   = 0.e0_wp
        surface_requests                  = 0
        requested_derivative_order        = 0
        output_excited_state_gradients    = 0.e0_wp
        output_ground_state_gradient      = 0.e0_wp
        output_diabatic_gradients         = 0.e0_wp
        output_excited_state_energies     = 0.e0_wp
        output_ground_state_energy        = 0.e0_wp
        output_diabatic_couplings         = 0.e0_wp
        excited_state_cartesian_gradients = 0.e0_wp
        ground_state_cartesian_gradient   = 0.e0_wp
        diabatic_cartesian_gradients      = 0.e0_wp
        last_requested_excited_state      = 0
        last_requested_coupling           = 0
        first_requested_excited_state     = 0
        first_requested_coupling          = 0
        z3cm_z                            = 0.e0_wp
        surface_requests = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            & 0, 0, 0, 0, 0, 0, 0, 0, 0], shape(surface_requests))
        requested_derivative_order = 0
        coordinate_representation  = 3
        maximum_excited_state      = 0
        maximum_derivative_order   = 1
        atom_count                 = 3
        phicm_deg1                 = 10.0e0_wp
        z3cm_z                     = [0.175e0_wp, 0.104e0_wp, 0.175e0_wp]
        precm_d                    = [141.196e0_wp, 109.449e0_wp, 141.196e0_wp]
        precm_re                   = [0.917e0_wp, 0.7419e0_wp, 0.917e0_wp]
        precm_beta                 = [2.21870e0_wp, 1.942e0_wp, 2.2187e0_wp]
        v                          = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        input_cartesian            = 0.e0_wp
        do iatom=1, natoms
            do idir=1, 3
                input_cartesian(iatom,idir) = x(iatom,idir)/bohr_to_angstrom
            enddo
        enddo
        requested_derivative_order = igrad
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

        !> @brief Evaluate prepare_local_context for the FH2_LEPS_NO2a_1984 potential.
        subroutine prepare_local_context
            use chempotpy_kinds, only: wp
            implicit none
            integer             :: i
            real(wp), parameter :: ckcau  = hartree_to_kilocalorie_per_mole
            real(wp), parameter :: cangau = 0.52917706e0_wp
            !   Echo potential parameters to UNIT NFLAG(18)
            if( atom_count .gt. 25 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            phicm_z0 = z3cm_z(1)
            do i = 1, 3
                !   CONVERT TO ATOMIC UNITS
                precm_d(i)    = precm_d(i)/ckcau
                precm_re(i)   = precm_re(i)/cangau
                precm_beta(i) = precm_beta(i)*cangau
                !   COMPUTE USEFUL CONSTANTS
                precm_zpo(i)   = 1.0e0_wp + z3cm_z(i)
                precm_op3z(i)  = 1.0e0_wp + 3.0e0_wp*z3cm_z(i)
                precm_top3z(i) = 2.0e0_wp*precm_op3z(i)
                precm_zp3(i)   = z3cm_z(i) + 3.0e0_wp
                precm_tzp3(i)  = 2.0e0_wp*precm_zp3(i)
                precm_do4z(i)  = precm_d(i)/4.0e0_wp/precm_zpo(i)
                precm_b(i)     = precm_beta(i)*precm_do4z(i)*2.0e0_wp
            enddo
            !   10 B(I) = BETA(I)*DO4Z(I)*2.0D0
            !    Set the values of the classical energy in the three asymptotic vall
            asycm_easyab   = precm_d(1)
            asycm_easybc   = precm_d(2)
            asycm_easyac   = precm_d(3)
            pt3cm_ezero(1) = precm_d(2)
            do i=1, 5
                infocm_ref(i) = ' '
            enddo
            !   Reference:
            infocm_ref(1) = 'D. G. Truhlar, B. C. Garrett, N. C. Blais'
            infocm_ref(2) = 'J. Chem. Phys. 80, 232(1984)'
        end subroutine prepare_local_context

        !> @brief Evaluate POT for the FH2_LEPS_NO2a_1984 potential.
        subroutine pot
            !   The potential energy in the AB valley, EASYAB, is equal to the poten
            !   energy of the H "infinitely" far from the FH diatomic, with the
            !   FH diatomic at its equilibrium configuration.  Similarly, the terms
            !   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,
            !   respectively.
            !   This potential is written such that:
            !                  R(1) = R(F-H)
            !                  R(2) = R(H-H)
            !                  R(3) = R(H-F)
            !   The zero of energy is defined at F "infinitely" far from the H2 diat
            use chempotpy_kinds, only: wp
            implicit none
            real(wp)            :: phi, rad, s, s1, t, t1, t2, t3, t4, t5, t6
            integer             :: i, j
            real(wp), parameter :: r2 = 1.41421356e0_wp
            cartesian_coordinates(1:atom_count,:) = input_cartesian(1:atom_count,:)
            call carttor
            !   Check the values of surface_requests and requested_derivative_order for validity.
            if( requested_derivative_order .eq. 1 )then
                do i=1, 3
                    ground_state_internal_gradient(i) = 0.e0_wp
                    do j=1, 3
                        precm_dcoul(i,j) = 0.e0_wp
                        precm_dexch(i,j) = 0.e0_wp
                    enddo
                enddo
            endif
            !   Compute the angle dependent Sato parameters
            call phical(phi)
            z3cm_z(1) = fphi(phi)
            z3cm_z(3) = z3cm_z(1)
            do i=1, 3, 2
                precm_zpo(i)   = 1.0e0_wp + z3cm_z(i)
                precm_op3z(i)  = 1.0e0_wp + 3.0e0_wp*z3cm_z(i)
                precm_top3z(i) = 2.0e0_wp*precm_op3z(i)
                precm_zp3(i)   = z3cm_z(i) + 3.0e0_wp
                precm_tzp3(i)  = 2.0e0_wp*precm_zp3(i)
                precm_do4z(i)  = precm_d(i)/4.0e0_wp/precm_zpo(i)
                precm_b(i)     = precm_do4z(i)*precm_beta(i)*2.e0_wp
            enddo
            do i = 1, 3
                precm_x(i)    = exp(-precm_beta(i)*(pt1cm_r(i)-precm_re(i)))
                precm_coul(i) = precm_do4z(i)*(precm_zp3(i)*precm_x(i)-precm_top3z(i))*precm_x(i)
                precm_exch(i) = precm_do4z(i)*(precm_op3z(i)*precm_x(i)-precm_tzp3(i))*precm_x(i)
            enddo
            rad    = sqrt((precm_exch(1)-precm_exch(2))**2+(precm_exch(2)-precm_exch(3))**2+ (precm_exch(3)- &
                & precm_exch(1))**2)
            ground_state_energy = -rad/r2
            ground_state_energy = ground_state_energy+precm_coul(1)+precm_coul(2)+precm_coul(3)+ pt3cm_ezero(1)
            !   If requested_derivative_order = 1, then compute the derivatives otherwise return
            if( requested_derivative_order .ne. 1 )then
                call eunitzero
                if( requested_derivative_order .ne. 0 )then
                    call rtocart
                    call dedcou
                endif
                return
            endif
            s = precm_exch(1)+precm_exch(2)+precm_exch(3)
            do i = 1, 3
                if( precm_x(i) .lt. 1.e-30_wp ) cycle
                t = (3.0e0_wp*precm_exch(i)-s)/r2*(precm_op3z(i)*precm_x(i)-precm_zp3(i))
                if( rad /= 0.0e0_wp )then
                    t = t/rad
                elseif( t /= 0.0e0_wp )then
                endif
                ground_state_internal_gradient(i) = precm_b(i)*precm_x(i)*(t-precm_zp3(i)*precm_x(i)+precm_op3z(i))
            enddo
            !   If the angle phi is less than or equal to 10 deg. then the
            !   Sato parameter is a constant and the derivative of the potential
            !   is equal to the standard leps derivative.  If phi is greater than
            !   10 deg. then the derivatives involve derivatives of the Sato
            !   parameter w.r.t. phi and the derivative of phi w.r.t. r1,r2,&r3.
            if( phicm_deg .le. phicm_deg1 )then
                call eunitzero
                if( requested_derivative_order .ne. 0 )then
                    call rtocart
                    call dedcou
                endif
                return
            endif
            do i=1, 3, 2
                t1 = precm_do4z(i)*precm_x(i)
                t2 = precm_coul(i)/precm_zpo(i)
                t3 = t1*(precm_x(i) - 6.e0_wp)
                t4 = (t3 - t2)*phicm_fder
                t5 = t1*(3.e0_wp*precm_x(i) - 2.0e0_wp)
                t6 = precm_exch(i) / precm_zpo(i)
                s1 = 3.e0_wp*precm_exch(i) - s
                do j=1, 3
                    precm_dcoul(i,j) = t4*coscm_dphi(j)
                    precm_dexch(i,j) = phicm_fder * s1* (t5 - t6) * coscm_dphi(j)
                enddo
            enddo
            t = -1.e0_wp / (r2 * rad)
            do i=1, 3
                do j=1, 3
                    ground_state_internal_gradient(i) = ground_state_internal_gradient(i) + precm_dcoul(j,i) + t* &
                        & precm_dexch(j,i)
                enddo
            enddo
            call eunitzero
            if( requested_derivative_order .ne. 0 )then
                call rtocart
                call dedcou
            endif
        end subroutine pot

        !> @brief Evaluate FPHI for the FH2_LEPS_NO2a_1984 potential.
        function fphi( phi )
            !   FH2 NEW SURFACE NO. 2   DON TRUHLAR
            !   AB,AC SATO PARAMETERS ARE Z(1)=Z(3) FOR DEG.LE.DEG1.
            !   DEG1 MUST BE SET IN prepare_local_context.
            !   Z0 IS SET EQUAL TO Z(1)=Z(3) IN prepare_local_context.
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: cosdeg, dd2, dd3, degdif, fphi, phi, s2, sindeg
            phicm_deg = 57.29577951e0_wp*phi
            if( phicm_deg .gt. 180.0e0_wp ) phicm_deg=phicm_deg-180.0e0_wp
            if( phicm_deg .lt. 0.0e0_wp ) phicm_deg=-phicm_deg
            if( phicm_deg .gt. 90.0e0_wp ) phicm_deg=180.0e0_wp-phicm_deg
            if( phicm_deg .lt. 0.0e0_wp )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            if( phicm_deg <= phicm_deg1 )then
                fphi       = phicm_z0
                phicm_fder = 0.0e0_wp
            elseif( phicm_deg <= 20.0e0_wp )then
                degdif     = phicm_deg-10.0e0_wp
                dd2        = degdif*degdif
                dd3        = degdif*dd2
                fphi       = 0.175e0_wp +1.15e-4_wp*dd2    -6.0e-6_wp*dd3
                phicm_fder = 2.30e-4_wp*degdif -1.8e-5_wp*dd2
            elseif( phicm_deg <= 60.0e0_wp )then
                degdif     = phicm_deg-41.0e0_wp
                fphi       = 0.191e0_wp +0.0005e0_wp*degdif
                phicm_fder = 0.0005e0_wp
            elseif( phicm_deg <= 90.0e0_wp )then
                sindeg     = sin(phi)
                s2         = sindeg*sindeg
                fphi       = 0.17569e0_wp +3.308e-2_wp*s2
                cosdeg     = cos(phi)
                phicm_fder = 6.616e-2_wp*sindeg*cosdeg
            else
                if( present(status) ) call status%set(status_invalid_input, &
                    'FPHI: angle exceeds 90 degrees')
                return
            endif
        end function fphi

        !> @brief Evaluate PHICAL for the FH2_LEPS_NO2a_1984 potential.
        subroutine phical( phi )
            !   FH2 NEW SURFACE NO. 2   DON TRUHLAR
            !   Derivatives of phi w.r.t. r(i) also calculated in this routine
            use chempotpy_kinds, only: wp
            implicit none
            real(wp) :: cos23, denom, fhh, fhh2, hr2, hr22, ovv2, phi, r12, r22, r32, sinphi, t, u, v
            integer  :: i
            hr2          = 0.5e0_wp*pt1cm_r(2)
            r12          = pt1cm_r(1)*pt1cm_r(1)
            r22          = pt1cm_r(2)*pt1cm_r(2)
            r32          = pt1cm_r(3)*pt1cm_r(3)
            hr22         = 0.25e0_wp*r22
            cos23        = (r12-r22-r32)/(-2.0e0_wp*pt1cm_r(2)*pt1cm_r(3))
            fhh2         = 0.5e0_wp*(r12 + r32 - 0.5e0_wp*r22)
            fhh          = sqrt(fhh2)
            u            = r32 - hr22 - fhh2
            v            = -2.0e0_wp * hr2 * fhh
            coscm_cosphi = u / v
            if( coscm_cosphi .gt. 1.0e0_wp ) coscm_cosphi=1.0e0_wp
            if( coscm_cosphi .lt. -1.0e0_wp ) coscm_cosphi=-1.0e0_wp
            phi = acos(coscm_cosphi)
            !   Calculate the derivatives of phi w.r.t. r(i)
            !   if requested_derivative_order = 1
            if( requested_derivative_order .ne. 1 ) return
            do i=1, 3
                coscm_dphi(i) = 0.e0_wp
            enddo
            if( abs(phi*57.29577951e0_wp) .le. 10.0e0_wp ) return
            sinphi        = sin(phi)
            ovv2          = 1.e0_wp/(v * v)
            denom         = -ovv2/sinphi
            t             = hr2 / fhh
            coscm_dphi(1) = pt1cm_r(1)*(u*t - v)*denom
            coscm_dphi(2) = u*(fhh - hr22/fhh)*denom
            coscm_dphi(3) = pt1cm_r(3)*(v + u*t)*denom
        end subroutine phical

        !> @brief Evaluate POTINFO for the FH2_LEPS_NO2a_1984 potential.
        subroutine potinfo
            implicit none
            integer            :: i, j, k, kmax
            integer, parameter :: isurf = 5
            do i =1, 5
            enddo
            kmax = 0
            do i = 1, isurf+1
                do j = 1, isurf+1
                    if( surface_requests(i,j) .ne. 0.and.kmax .lt. max(i, j) ) kmax = max(i, j)
                enddo
            enddo
            if( kmax-1 .gt. maximum_excited_state )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
            first_requested_excited_state = 0
            last_requested_excited_state  = 0
            do i = 2, isurf+1
                if( surface_requests(i,i) .ne. 0 )then
                    last_requested_excited_state = i-1
                    if( first_requested_excited_state .eq. 0 ) first_requested_excited_state = i-1
                endif
            enddo
            first_requested_coupling = 0
            last_requested_coupling  = 0
            k                        = 0
            do i = 1, isurf
                do j = i+1, isurf+1
                    k = k+1
                    if( surface_requests(i,j)+surface_requests(j,i) .ne. 0 )then
                        last_requested_coupling = k
                        if( first_requested_coupling .eq. 0 ) first_requested_coupling = k
                    endif
                enddo
            enddo
            if( requested_derivative_order .gt. maximum_derivative_order )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine potinfo

        !> @brief Evaluate CARTTOR for the FH2_LEPS_NO2a_1984 potential.
        subroutine carttor
            use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates
            implicit none
            integer :: coordinate_status
            call cartesian_to_potlib_coordinates(cartesian_coordinates, atom_count, coordinate_representation, pt1cm_r, coordinate_status)
            if( coordinate_status /= 0 )then
                if( present(status) ) call status%set(status_invalid_input, &
                    'legacy kernel rejected the evaluation request')
                return
            endif
        end subroutine carttor

        !> @brief Evaluate EUNITZERO for the FH2_LEPS_NO2a_1984 potential.
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

        !> @brief Evaluate RTOCART for the FH2_LEPS_NO2a_1984 potential.
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
                            & 1))/pt1cm_r(m1)
                        termy                                = (cartesian_coordinates(i,2)-cartesian_coordinates(j, &
                            & 2))/pt1cm_r(m1)
                        termz                                = (cartesian_coordinates(i,3)-cartesian_coordinates(j, &
                            & 3))/pt1cm_r(m1)
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
                    ygs(i) = ground_state_internal_gradient(i)/pt1cm_r(i)
                    if( first_requested_excited_state .ne. 0 )then
                        do j=first_requested_excited_state, last_requested_excited_state
                            yes(i, j) = excited_state_internal_gradients(i,j)/pt1cm_r(i)
                        enddo
                    endif
                    if( first_requested_coupling .ne. 0 )then
                        do k=first_requested_coupling, last_requested_coupling
                            yij(i, k) = diabatic_internal_gradients(i,k)/pt1cm_r(i)
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

        !> @brief Evaluate DEDCOU for the FH2_LEPS_NO2a_1984 potential.
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
end module chempotpy_fh2_leps_no2a_1984_pes

!> @brief Expose FH2_LEPS_NO2a_1984 through the Python-facing PES entry point.
subroutine pes( x, igrad, p, g, d )
    use chempotpy_fh2_leps_no2a_1984_pes, only: evaluate
    use chempotpy_kinds,                  only: wp
    implicit none
    integer, parameter    :: natoms  = 3
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: igrad
    real(wp), intent(in)  :: x(natoms,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
