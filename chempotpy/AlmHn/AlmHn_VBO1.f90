!   System:                        Al, H, Al-H
!   Functional form:               Valence-Bond Order
!   Common name:                   VBO.beta1
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!   Notes:  Valence-Bond Order potential energy funciton. The functional
!           form is from Ref. 1.
!  References: (1)"valence-Bond Order(VBO): A New Approach to Modeling
!                 Reactive Potential Energy Surfaces for Materials and N
!                 Meiyu Zhao, Mark A. Iron, Przemys?aw Staszewski, Natha
!                Rosendo Valero, and Donald G. Truhlar, manuscript in pr
!  Units:
!       energies    - hartrees
!       coordinates - bohrs
!  v       --- energy of the system (output)
!  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!              of the system (input)
!  natom   --- number of atoms in the system (input)
!  maxatom --- maximal number of atoms (input)
!  symbol  --- character*2 array of dimension holds the atomic symbols.
! This version of the code corrects two typos around lines 64 and 65. Ja

!> @brief Implement AlmHn_VBO1, a potential energy surface for AlmHn.
module chempotpy_almhn_vbo1_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the AlmHn_VBO1 potential energy surface for AlmHn.
    subroutine evaluate( xyz, nal, nh, natoms, igrad, p, g, d )
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        use chempotpy_kinds, only: wp
        implicit none
        ! number of electronic state
        integer, parameter    :: nstates = 1
        integer, intent(in)   :: nal, nh, natoms
        real(wp), intent(in)  :: xyz(10000,3)
        integer, intent(in)   :: igrad
        real(wp), intent(out) :: p(nstates), g(nstates,natoms,3)
        real(wp), intent(out) :: d(nstates,nstates,natoms,3)
        real(wp)              :: cx(natoms), cy(natoms), cz(natoms)
        real(wp)              :: v
        character(len=2)      :: symbol(natoms)
        integer               :: istate, iatom
        !initialize
        v = 0.e0_wp; g = 0.e0_wp; d = 0.e0_wp
        do iatom=1, natoms
            cx(iatom) = xyz(iatom,1)/bohr_to_angstrom
            cy(iatom) = xyz(iatom,2)/bohr_to_angstrom
            cz(iatom) = xyz(iatom,3)/bohr_to_angstrom
        enddo
        do iatom=1, nal
            symbol(iatom) = 'Al'
        enddo
        do iatom=1, nh
            symbol(iatom+nal) = 'H'
        enddo
        call alhpot(cx, cy, cz, v, natoms, 10000, symbol)
        if( igrad == 0 )then
            do istate=1, nstates
                p(istate) = v*hartree_to_electronvolt
            enddo
        else
            write (*, *) 'Only energy is available'
        endif
    end subroutine evaluate

    !> @brief Evaluate alhpot for the AlmHn_VBO1 potential.
    subroutine alhpot( x, y, z, v, natom, maxatom, symbol )
        use chempotpy_geometry, only: fill_distance_matrix
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp) :: add, c1, c1_al, c1_h, c1_hal, c2, c2_al, c2_h, c2_hal, c3, c3_al, c3_h, c4, c4_al, c4_h, &
            & del1, den1, den2, den3, den4, dx, dy, dz, en, en1, en2, esn, esn_al, esn_h, fa, gam1_al, gam1_h, gam1_hal, &
            & gam2_al, gam2_h, gam2_hal, gam3_al, gam3_alh, gam3_h, gam3_hal, gam4_al, gam4_alh, gam4_h, gam4_hal, &
            & r(maxatom, maxatom), r1, rij, sig3_alh, sig3_hal, sig4_alh, sig4_hal, v, v1, vi3, vi4, x(maxatom), y(maxatom), &
            & z(maxatom)
        integer             :: i, id(maxatom), j, k, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        character(len=2)    :: symbol(maxatom)
        real(wp)            :: b(4,maxatom,maxatom), sigma(4,maxatom,maxatom)
        real(wp)            :: cn(4,maxatom,maxatom), gam(4,maxatom,maxatom)
        !   gam1_h to gam4_h are gammaij in Eq. (1) of Ref. 1
        !   c1_h to c4_h are ci,1 to ci,4 in Eqs. (6) and (7)
        !   esn_h is the exponent of Vi in Eq. (6)
        gam1_h = 2.7292
        gam2_h = 8.8100
        gam3_h = 16.1462
        gam4_h = 16.1330
        c1_h   = 2.0112/autoev
        c2_h   = 0.8084/autoev
        esn_h  = 0.1536
        c3_h   = 0.00001
        c4_h   = 0.00001
        !   The same definitions apply for Aluminum
        gam1_al = 0.3599
        gam2_al = 4.2538
        gam3_al = 0.5684
        gam4_al = 2.9562
        c1_al   = 0.2248/autoev
        c2_al   = 0.2172/autoev
        esn_al  = 0.7931
        c3_al   = 0.0133
        c4_al   = 0.0068
        !   gam1_hal to gam4_hal and gam3_alh, gam4_alh are gammaij in Eq. (1)
        !   c1_hal and c2_hal are ci,1 and ci,2 in Eqs. (9) and (10)
        !   esn_h is the exponent of Vi in Eq. (9)
        !   sig3_hal to sig4_alh are the sigma in Eq. (10)
        gam1_hal = 7.4891
        gam2_hal = 1.8993
        gam3_hal = 13.5669
        gam4_hal = 13.6902
        gam3_alh = 1.6911
        gam4_alh = 2.8674
        sig3_hal = 0.00001
        sig4_hal = 0.0191
        sig3_alh = 2.9306
        sig4_alh = 0.00008
        c1_hal   = 0.3509/autoev
        c2_hal   = 0.5174/autoev
        !   Determine which atoms in the data file are Aluminum and which Hydrog
        do i=1, natom
            if((symbol(i).eq.'H ').or.(symbol(i).eq.'h ').or. &
                & (symbol(i).eq.' H').or.(symbol(i).eq.' h'))then
            id(i) = 1
        elseif((symbol(i).eq.'AL').or.(symbol(i).eq.'Al').or. &
                & (symbol(i).eq.'al').or.(symbol(i).eq.'aL'))then
            id(i) = 2
        else
            write(13, *)'Atomic symbol ', symbol(i), ' undefined. Stopping.'
            stop "Undefined atomic symbol."
            endif
        enddo
        !    Calculate all internuclear distances
        call fill_distance_matrix(x(1:natom), y(1:natom), z(1:natom), r(1:natom, 1:natom))
        !    Calculate en1 which corresponds to the first term in parenthesis in
        en1 = 0e0_wp
        en2 = 0e0_wp
        do i=1, natom
            v1 = 0e0_wp
            do j=1, natom
                if( i .eq. j ) cycle
                if( (id(i) .eq. 1).and.(id(j) .eq. 1) )then
                    del1       = 5.29/autoang
                    gam(1,i,j) = gam1_h
                    gam(2,i,j) = gam2_h
                    gam(3,i,j) = gam3_h
                    gam(4,i,j) = gam4_h
                    en         = 0.5
                    c1         = c1_h
                    c2         = c2_h
                    r1         = 0.7414/autoang
                    !    Calculate Nij according to Eq. (2)
                    den1         = 1e0_wp-(r1/del1)**en
                    cn(1,i,j)    = exp(gam(1,i,j)/den1)
                    den2         = 1e0_wp-(r1/del1)**en
                    cn(2,i,j)    = exp(gam(2,i,j)/den2)
                    den3         = 1e0_wp-(r1/del1)**en
                    cn(3,i,j)    = exp(gam(3,i,j)/den3)
                    den4         = 1e0_wp-(r1/del1)**en
                    cn(4,i,j)    = exp(gam(4,i,j)/den4)
                    sigma(3,i,j) = 1.0
                    sigma(4,i,j) = 1.0
                elseif( (id(i) .eq. 2).and.(id(j) .eq. 2) )then
                    del1         = 6.88/autoang
                    gam(1,i,j)   = gam1_al
                    gam(2,i,j)   = gam2_al
                    gam(3,i,j)   = gam3_al
                    gam(4,i,j)   = gam4_al
                    en           = 0.5
                    c1           = c1_al
                    c2           = c2_al
                    r1           = 2.7306/autoang
                    den1         = 1e0_wp-(r1/del1)**en
                    cn(1,i,j)    = exp(gam(1,i,j)/den1)
                    den2         = 1e0_wp-(r1/del1)**en
                    cn(2,i,j)    = exp(gam(2,i,j)/den2)
                    den3         = 1e0_wp-(r1/del1)**en
                    cn(3,i,j)    = exp(gam(3,i,j)/den3)
                    den4         = 1e0_wp-(r1/del1)**en
                    cn(4,i,j)    = exp(gam(4,i,j)/den4)
                    sigma(3,i,j) = 1.0
                    sigma(4,i,j) = 1.0
                else
                    del1       = 6.88/autoang
                    gam(1,i,j) = gam1_hal
                    gam(2,i,j) = gam2_hal
                    if( (id(i) .eq. 1).and.(id(j) .eq. 2) )then
                        gam(3,i,j)   = gam3_hal
                        gam(4,i,j)   = gam4_hal
                        sigma(3,i,j) = sig3_hal
                        sigma(4,i,j) = sig4_hal
                    elseif( (id(i) .eq. 2).and.(id(j) .eq. 1) )then
                        gam(3,i,j)   = gam3_alh
                        gam(4,i,j)   = gam4_alh
                        sigma(3,i,j) = sig3_alh
                        sigma(4,i,j) = sig4_alh
                    endif
                    en        = 0.5
                    c1        = c1_hal
                    c2        = c2_hal
                    r1        = 1.6637/autoang
                    den1      = 1e0_wp-(r1/del1)**en
                    cn(1,i,j) = exp(gam(1,i,j)/den1)
                    den2      = 1e0_wp-(r1/del1)**en
                    cn(2,i,j) = exp(gam(2,i,j)/den2)
                    den3      = 1e0_wp-(r1/del1)**en
                    cn(3,i,j) = exp(gam(3,i,j)/den3)
                    den4      = 1e0_wp-(r1/del1)**en
                    cn(4,i,j) = exp(gam(4,i,j)/den4)
                endif
                !   Calculate bij according to Eq. (1)
                rij = r(i, j)
                if( rij .lt. del1 )then
                    b(1,i,j) = exp(-gam(1,i,j)/(1e0_wp-(rij/del1)**en))*cn(1,i,j)
                    b(2,i,j) = exp(-gam(2,i,j)/(1e0_wp-(rij/del1)**en))*cn(2,i,j)
                    b(3,i,j) = exp(-gam(3,i,j)/(1e0_wp-(rij/del1)**en))*cn(3,i,j)
                    b(4,i,j) = exp(-gam(4,i,j)/(1e0_wp-(rij/del1)**en))*cn(4,i,j)
                    fa       = c1*b(1,i,j)+c2*b(2,i,j)
                    v1       = v1+fa
                else
                    do k=1, 4
                        b(k,i,j) = 0.0e0_wp
                    enddo
                endif
            enddo
            en1 = en1+v1
        enddo
        !    Calculate en2 which corresponds to the last term in Eq. (6) (pure A
        do i=1, natom
            if( id(i) .eq. 1 )then
                c3  = c3_h
                c4  = c4_h
                esn = esn_h
            elseif( id(i) .eq. 2 )then
                c3  = c3_al
                c4  = c4_al
                esn = esn_al
            endif
            vi3 = 0.0e0_wp
            vi4 = 0.0e0_wp
            do j=1, natom
                if( i .eq. j ) cycle
                vi3 = vi3+c3*sigma(3,i,j)*b(3,i,j)
                vi4 = vi4+c4*sigma(4,i,j)*b(4,i,j)
            enddo
            add = (vi3+vi4)**esn
            en2 = en2-add
        enddo
        v = en1+en2
    end subroutine alhpot
end module chempotpy_almhn_vbo1_pes

!> @brief Expose AlmHn_VBO1 through the Python-facing PES entry point.
subroutine pes( xyz, nal, nh, natoms, igrad, p, g, d )
    use chempotpy_almhn_vbo1_pes, only: evaluate
    use chempotpy_kinds,          only: wp
    implicit none
    integer, parameter    :: nstates = 1
    integer, intent(in)   :: nal, nh, natoms, igrad
    real(wp), intent(in)  :: xyz(10000,3)
    real(wp), intent(out) :: p(nstates)
    real(wp), intent(out) :: g(nstates,natoms,3)
    real(wp), intent(out) :: d(nstates,nstates,natoms,3)
    call evaluate(xyz, nal, nh, natoms, igrad, p, g, d)
end subroutine pes
