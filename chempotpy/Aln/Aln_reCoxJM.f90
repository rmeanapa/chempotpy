!   System:                        Aln
!   Functional form:               Extended Rydberg plus generalized three-body
!   Common name:                   reCoxJM
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Many-body aluminum potential energy function.  The parameters
!           for this PEF have been parameterized against a data set of
!           aluminum cluster energies and bulk data.
!  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
!              and D. G. Truhlar, "Analytic Potentials Energy Functions
!              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
!              H. Cox, R. L. Johnston, and J. N. Murrell,
!              Surf. Sci., 373, 67 (1997).
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximal number of atoms (input)

!> @brief Own the complete implementation of this aluminum PES.
module chempotpy_aln_recoxjm_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate an aluminum-cluster potential through the Python contract.
    subroutine evaluate( coordinates_angstrom, atom_count, derivative_order, &
            energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
        use chempotpy_units, only: coordinates_to_bohr, &
            energy_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        integer, parameter    :: state_count = 1
        real(wp), intent(in)  :: coordinates_angstrom(10000,3)
        integer, intent(in)   :: atom_count, derivative_order
        real(wp), intent(out) :: energy_ev(state_count)
        real(wp), intent(out) :: gradient_ev_per_angstrom( &
            state_count,atom_count,3 )
        real(wp), intent(out) :: coupling_per_angstrom( &
            state_count,state_count,atom_count,3 )
        real(wp) :: coordinates_bohr(atom_count,3), energy_hartree
        energy_ev                = 0.0e0_wp
        gradient_ev_per_angstrom = 0.0e0_wp
        coupling_per_angstrom    = 0.0e0_wp
        call coordinates_to_bohr( &
            coordinates_angstrom(1:atom_count,:), coordinates_bohr)
        call pot(coordinates_bohr(:,1), coordinates_bohr(:,2), &
            coordinates_bohr(:,3), energy_hartree, atom_count, atom_count)
        if( derivative_order == 0 )then
            energy_ev(1) = energy_to_electronvolt(energy_hartree)
        else
            write(*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate pot for the Aln_reCoxJM potential.
    subroutine pot( x, y, z, v, natom, maxatom )
        use chempotpy_geometry, only: fill_distance_matrix
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp) :: a2, a3, c0, c1, c2, c3, c4, c5, c6, dd, dx, dy, dz, fsech, oneha, onesi, oneth, poly, q1, q2, q23, &
            & q3, q332, r(maxatom, maxatom), re, rij, rjk, rki, rr, twoth, v, v2, v3, x(maxatom), y(maxatom), z(maxatom)
        integer             :: i, j, k, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        !  calculate and tabulate the distances between all atoms
        ! original parameters
        !      a2=7.d0
        !      a3=8.d0
        !      dd=0.9073d0/autoev
        !      re=2.7568d0/autoang
        !      c0= 0.2525d0
        !      c1=-0.4671d0
        !      c2= 4.4903d0
        !      c3=-1.1717d0
        !      c4= 1.6498d0
        !      c5=-5.3579d0
        !      c6= 1.6327d0
        ! new parameters
        a2    = 4.75085490962383972
        a3    = 7.73864191499755716
        dd    = 1.25012212994626282/autoev
        re    = 2.61700048851978506/autoang
        c0    = 0.182247191011235959
        c1    = -0.319491939423546645
        c2    = 9.46849047386419151
        c3    = -0.502198339032730923
        c4    = 1.59037616023448947
        c5    = -4.89252564728871508
        c6    = 1.19980459208597945
        oneth = sqrt(1.e0_wp/3.e0_wp)
        oneha = sqrt(1.e0_wp/2.e0_wp)
        twoth = sqrt(2.e0_wp/3.e0_wp)
        onesi = sqrt(1.e0_wp/6.e0_wp)
        call fill_distance_matrix(x(1:natom), y(1:natom), z(1:natom), r(1:natom, 1:natom))
        !  calculate the potential energy
        v  = 0.e0_wp
        v2 = 0.e0_wp
        v3 = 0.e0_wp
        do i=1, natom
            do j=i+1, natom
                !  two body term
                rr  = r(i, j)
                rij = (rr-re)/re
                v2  = v2-dd*(1.e0_wp+a2*rij)*exp(-a2*rij)
                do k=j+1, natom
                    ! three body term
                    rjk   = (r(j, k)-re)/re
                    rki   = (r(k, i)-re)/re
                    q1    = oneth*(rij+rjk+rki)
                    q2    = oneha*(rjk-rki)
                    q3    = twoth*rij-onesi*(rjk+rki)
                    q23   = q2**2+q3**2
                    q332  = q3**3-3.e0_wp*q3*q2**2
                    fsech = 1.e0_wp/cosh(a3*q1)
                    poly  = c0+c1*q1+c2*q1**2+c3*q23+c4*q1**3+c5*q1*q23+c6*q332
                    v3    = v3+dd*poly*fsech
                enddo
            enddo
        enddo
        v = v2+v3
    end subroutine pot
end module chempotpy_aln_recoxjm_pes

!> @brief Expose Aln_reCoxJM through the Python-facing PES entry point.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_aln_recoxjm_pes, only: evaluate
    use chempotpy_kinds,           only: wp
    implicit none
    integer, parameter    :: state_count = 1
    real(wp), intent(in)  :: coordinates_angstrom(10000,3)
    integer, intent(in)   :: atom_count, derivative_order
    real(wp), intent(out) :: energy_ev(state_count)
    real(wp), intent(out) :: gradient_ev_per_angstrom(state_count,atom_count,3)
    real(wp), intent(out) :: coupling_per_angstrom(state_count,state_count,atom_count,3)
    call evaluate(coordinates_angstrom, atom_count, derivative_order, energy_ev, gradient_ev_per_angstrom, &
        & coupling_per_angstrom)
end subroutine pes
