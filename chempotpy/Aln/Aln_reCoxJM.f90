!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate an aluminum-cluster potential through the Python contract.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
        energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_legacy_units, only: legacy_coordinates_to_bohr, &
        legacy_energy_to_electronvolt
    use chempotpy_kinds, only: real64
    implicit none
    integer, parameter        :: state_count = 1
    real(real64), intent(in)  :: coordinates_angstrom(10000,3)
    integer, intent(in)       :: atom_count, derivative_order
    real(real64), intent(out) :: energy_ev(state_count)
    real(real64), intent(out) :: gradient_ev_per_angstrom( &
        state_count,atom_count,3 )
    real(real64), intent(out) :: coupling_per_angstrom( &
        state_count,state_count,atom_count,3 )
    real(real64) :: coordinates_bohr(atom_count,3), energy_hartree
    energy_ev                = 0.0d0
    gradient_ev_per_angstrom = 0.0d0
    coupling_per_angstrom    = 0.0d0
    call legacy_coordinates_to_bohr( &
        coordinates_angstrom(1:atom_count,:), coordinates_bohr )
    call pot( coordinates_bohr(:,1), coordinates_bohr(:,2), &
        coordinates_bohr(:,3), energy_hartree, atom_count, atom_count )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*, *) 'Only energy is available'
    end if
end subroutine pes

!> @brief Evaluate pot for the Aln_reCoxJM potential.
subroutine pot( x, y, z, v, natom, maxatom )
    !   System:                        Aln
    !   Functional form:               Extended Rydberg plus generalized three-body
    !   Common name:                   reCoxJM
    !   Number of derivatives:         0
    !   Number of bodies:              variable
    !   Number of electronic surfaces: 1
    !   Interface:                     HO-MM-0
    !
    !   Notes:  Many-body aluminum potential energy function.  The parameters
    !           for this PEF have been parameterized against a data set of
    !           aluminum cluster energies and bulk data.
    !
    !  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
    !              and D. G. Truhlar, "Analytic Potentials Energy Functions
    !              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
    !              H. Cox, R. L. Johnston, and J. N. Murrell,
    !              Surf. Sci., 373, 67 (1997).
    !
    !  Units:
    !       energies    - hartrees
    !       coordinates - bohrs
    !
    !  v       --- energy of the system (output)
    !  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
    !              of the system (input)
    !  natom   --- number of atoms in the system (input)
    !  maxatom --- maximal number of atoms (input)

    implicit none
    double precision :: a2, a3, c0, c1, c2, c3, c4, c5, c6, dd, dx, dy, dz, &
        & fsech, oneha, onesi, oneth, poly, q1, q2, q23, q3, q332, r, re, rij, rjk, rki, rr, &
        & twoth, v, v2, v3, x, y, z
    integer :: i, j, k, maxatom, natom
    dimension x(maxatom), y(maxatom), z(maxatom), r(maxatom,maxatom)
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
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
    a2 = 4.75085490962383972
    a3 = 7.73864191499755716
    dd = 1.25012212994626282/autoev
    re = 2.61700048851978506/autoang
    c0 = 0.182247191011235959
    c1 = -0.319491939423546645
    c2 = 9.46849047386419151
    c3 = -0.502198339032730923
    c4 = 1.59037616023448947
    c5 = -4.89252564728871508
    c6 = 1.19980459208597945

    oneth = dsqrt(1.d0/3.d0)
    oneha = dsqrt(1.d0/2.d0)
    twoth = dsqrt(2.d0/3.d0)
    onesi = dsqrt(1.d0/6.d0)

    do 1 i=1, natom
        r(i,i) = 0d0
        do 1 j=i+1, natom
            dx     = x(i)-x(j)
            dy     = y(i)-y(j)
            dz     = z(i)-z(j)
            r(i,j) = dsqrt(dx*dx+dy*dy+dz*dz)
    1 r(j,i)=r(i,j)

    !  calculate the potential energy
    v  = 0.d0
    v2 = 0.d0
    v3 = 0.d0
    do 2 i=1, natom
        do 2 j=i+1, natom
            !  two body term
            rr  = r(i,j)
            rij = (rr-re)/re
            v2  = v2-dd*(1.d0+a2*rij)*dexp(-a2*rij)
            do 2 k=j+1, natom
                ! three body term
                rjk   = (r(j,k)-re)/re
                rki   = (r(k,i)-re)/re
                q1    = oneth*(rij+rjk+rki)
                q2    = oneha*(rjk-rki)
                q3    = twoth*rij-onesi*(rjk+rki)
                q23   = q2**2+q3**2
                q332  = q3**3-3.d0*q3*q2**2
                fsech = 1.d0/dcosh(a3*q1)
                poly  = c0+c1*q1+c2*q1**2+c3*q23+c4*q1**3+c5*q1*q23+c6*q332
                v3    = v3+dd*poly*fsech
    2 continue

    v = v2+v3

end subroutine pot
