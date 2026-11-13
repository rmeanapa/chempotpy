!   System:                        AlmHn
!   Functional form:               Valence--Bond Order
!   Common name:                   VBO2.beta1
!   Number of derivatives:         0
!   Number of bodies:              variable
!   Number of electronic surfaces: 1
!   Interface:                     HO-MM-0
!     Notes:  Valence--Bond Order potential energy function. The functio
!             form is from Ref. 1.
!    References: (1) Valence--Bond Order(VBO): A New Approach to Modelin
!                   Reactive Potential Energy Surfaces for Materials and
!                   Meiyu Zhao, Mark A. Iron, Przemyslaw Staszewski, Nat
!                  Rosendo Valero, and Donald G. Truhlar, manuscript in
!    Units:
!         energies    - hartrees
!         coordinates - bohrs
!    v       --- energy of the system (output)
!    x,y,z   --- one-dimensional arrays of the Cartesian coordinates
!                of the system (input)
!    natom   --- number of atoms in the system (input)
!    maxatom --- maximal number of atoms (input)
!    symbol  --- character*2 array of dimension holds the atomic symbols

!> @brief Implement AlmHn_VBO2, a potential energy surface for AlmHn.
module chempotpy_almhn_vbo2_pes
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the AlmHn_VBO2 potential energy surface for AlmHn.
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

    !> @brief Evaluate alhpot for the AlmHn_VBO2 potential.
    subroutine alhpot( x, y, z, v, natom, maxatom, symbol )
        use chempotpy_geometry, only: fill_distance_matrix
        use chempotpy_kinds,    only: wp
        implicit none
        real(wp) :: add, alpha, alpha_al, alpha_alh, alpha_h, beta, beta_al, beta_alh, beta_h, ca, ca_al, ca_bar, &
            & ca_bar_al, ca_bar_h, ca_h, cb, cb_al, cb_bar, cb_bar_al, cb_bar_h, cb_h, cgam, cgam_al, cgam_bar, cgam_bar_al, &
            & cgam_bar_h, cgam_h, del1, den1, den2, den3, den4, den5, dx, dy, dz, en, en1, en2, fa, gam1_al, gam1_alh, gam1_h, &
            & gam2_al, gam2_alh, gam2_h, gam3_al, gam3_h, gam4_al, gam4_alh, gam4_h, gam5_alh, r(maxatom, maxatom), r1, rij, &
            & sig5_alh, v, v1, vi3, vi4, vi5, x(maxatom), y(maxatom), z(maxatom)
        integer             :: i, id(maxatom), j, k, maxatom, natom
        real(wp), parameter :: autoev  = 27.2113961e0_wp
        real(wp), parameter :: autoang = 0.529177249e0_wp
        character(len=2)    :: symbol(maxatom)
        real(wp)            :: b(5,maxatom,maxatom), sigma(5,maxatom,maxatom)
        real(wp)            :: kap3(2,2,2), kap4(2,2,2), kap5(2,2,2)
        real(wp)            :: kap4_hhal, kap4_halal
        real(wp)            :: cn(5,maxatom,maxatom), gam(5,maxatom,maxatom)
        !     gam1_h to gam4_h are gammaij in Eq. (1) of Ref. 1
        !     alpha_h and beta_h correspond to alphaij and betaij in Eq. (11)
        !     ca_h is Ai, cb_h is Bi, ca_bar_h is Ai(bar), and cb_bar_h is Bi(ba
        gam1_h   = 2.2671
        gam2_h   = 8.5293
        gam3_h   = 0.9651
        gam4_h   = 17.1086
        alpha_h  = 2.2520/autoev
        beta_h   = 0.3981
        ca_h     = 0.0948/autoev
        ca_bar_h = 2.4844
        cb_h     = 57.1859
        cb_bar_h = 0.1358
        !      The same definitions apply to Aluminum
        gam1_al   = 0.2626
        gam2_al   = 4.6212
        gam3_al   = 0.6588
        gam4_al   = 1.4696
        alpha_al  = 0.0707/autoev
        beta_al   = 1.5626
        ca_al     = 0.6096/autoev
        ca_bar_al = 0.6848
        cb_al     = 0.7314
        cb_bar_al = 0.6504
        !     gam1_alh to gam4_alh are gamij in Eq. (1)
        !     alpha_alh and beta_alh correspond to alphaij and betaij in Eq. (13
        !     cgam_h(al) is Gammai, cgam_bar_h(al) is Gammai(bar)
        !     kap4_hhal and kap4_halal are kappa in Eq. (14), and sig5_alh are s
        gam1_alh = 2.1414
        gam2_alh = 7.1248
        !          gam3_alh=
        gam4_alh    = 7.0492
        gam5_alh    = 2.4011
        alpha_alh   = 3.6460/autoev
        beta_alh    = 0.1161
        cgam_h      = 215.0668
        cgam_bar_h  = 0.9037
        cgam_al     = 27.7186
        cgam_bar_al = 0.9035
        kap4_hhal   = 0.0016
        kap4_halal  = 0.0031
        sig5_alh    = 0.1905
        !     Determine which atoms in the data file are Aluminum and which Hydr
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
        ! parameters for kappa.
        ! kap3(4,5) is kappa for p=3,4,5 in Eq. (14)
        kap3(1,1,1) = 0.0e0_wp
        kap3(2,2,2) = 0.0e0_wp
        kap3(1,1,2) = 0.e0_wp
        kap3(1,2,1) = 0.e0_wp
        kap3(2,1,1) = 0.e0_wp
        kap3(1,2,2) = 0.e0_wp
        kap3(2,1,2) = 0.e0_wp
        kap3(2,2,1) = 0.e0_wp
        kap4(1,1,1) = 0.0e0_wp
        kap4(2,2,2) = 0.0e0_wp
        kap4(1,1,2) = kap4_hhal
        kap4(1,2,1) = kap4_hhal
        kap4(2,1,1) = kap4_hhal
        kap4(1,2,2) = kap4_halal
        kap4(2,1,2) = kap4_halal
        kap4(2,2,1) = kap4_halal
        kap5(1,1,1) = 0.0e0_wp
        kap5(2,2,2) = 0.0e0_wp
        kap5(1,1,2) = 0.e0_wp
        kap5(1,2,1) = 0.e0_wp
        kap5(2,1,1) = 0.e0_wp
        kap5(1,2,2) = 0.e0_wp
        kap5(2,1,2) = 0.e0_wp
        kap5(2,2,1) = 0.e0_wp
        !    Calculate all internuclear distances
        call fill_distance_matrix(x(1:natom), y(1:natom), z(1:natom), r(1:natom, 1:natom))
        !    Calculate en1 which corresponds to the term in brackets in Eq. (11)
        en1 = 0e0_wp
        en2 = 0e0_wp
        do i=1, natom
            v1 = 0e0_wp
            do j=1, natom
                if( i .eq. j ) cycle
                if( (id(i) .eq. 1).and.(id(j) .eq. 1) )then
                    del1         = 5.29/autoang
                    gam(1,i,j)   = gam1_h
                    gam(2,i,j)   = gam2_h
                    gam(3,i,j)   = gam3_h
                    gam(4,i,j)   = gam4_h
                    en           = 0.5
                    alpha        = alpha_h
                    beta         = beta_h
                    r1           = 0.7414/autoang
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
                    sigma(5,i,j) = 0.0
                elseif( (id(i) .eq. 2).and.(id(j) .eq. 2) )then
                    del1         = 6.88/autoang
                    gam(1,i,j)   = gam1_al
                    gam(2,i,j)   = gam2_al
                    gam(3,i,j)   = gam3_al
                    gam(4,i,j)   = gam4_al
                    en           = 0.5
                    alpha        = alpha_al
                    beta         = beta_al
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
                    sigma(5,i,j) = 0.0
                else
                    del1       = 6.88/autoang
                    gam(1,i,j) = gam1_alh
                    gam(2,i,j) = gam2_alh
                    !          gam(3,i,j)=gam3_alh
                    gam(4,i,j) = gam4_alh
                    gam(5,i,j) = gam5_alh
                    en         = 0.5
                    alpha      = alpha_alh
                    beta       = beta_alh
                    r1         = 1.6637/autoang
                    den1       = 1e0_wp-(r1/del1)**en
                    cn(1,i,j)  = exp(gam(1,i,j)/den1)
                    den2       = 1e0_wp-(r1/del1)**en
                    cn(2,i,j)  = exp(gam(2,i,j)/den2)
                    !           den3=1d0-(r1/del1)**en
                    !           cn(3,i,j)=dexp(gam(3,i,j)/den3)
                    den4         = 1e0_wp-(r1/del1)**en
                    cn(4,i,j)    = exp(gam(4,i,j)/den4)
                    den5         = 1e0_wp-(r1/del1)**en
                    cn(5,i,j)    = exp(gam(5,i,j)/den5)
                    sigma(3,i,j) = 0.0
                    sigma(4,i,j) = 0.0
                    sigma(5,i,j) = sig5_alh
                endif
                !   Calculate bij according to Eq. (1)
                rij = r(i, j)
                if( rij .lt. del1 )then
                    b(1,i,j) = exp(-gam(1,i,j)/(1e0_wp-(rij/del1)**en))*cn(1,i,j)
                    b(2,i,j) = exp(-gam(2,i,j)/(1e0_wp-(rij/del1)**en))*cn(2,i,j)
                    b(3,i,j) = exp(-gam(3,i,j)/(1e0_wp-(rij/del1)**en))*cn(3,i,j)
                    b(4,i,j) = exp(-gam(4,i,j)/(1e0_wp-(rij/del1)**en))*cn(4,i,j)
                    if( id(i) .ne. id(j) )then
                        b(3,i,j) = 0.0
                        b(5,i,j) = exp(-gam(5,i,j)/(1e0_wp-(rij/del1)**en))*cn(5,i,j)
                    endif
                    fa = alpha*(b(1,i,j)+beta*b(2,i,j))
                    v1 = v1+fa
                else
                    do k=1, 5
                        b(k,i,j) = 0.0e0_wp
                    enddo
                endif
            enddo
            en1 = en1+v1
        enddo
        !    Calculate en2 which corresponds to the second term in Eq. (11) (pur
        do i=1, natom
            if( id(i) .eq. 1 )then
                ca       = ca_h
                cb       = cb_h
                ca_bar   = ca_bar_h
                cb_bar   = cb_bar_h
                cgam     = cgam_h
                cgam_bar = cgam_bar_h
            elseif( id(i) .eq. 2 )then
                ca       = ca_al
                cb       = cb_al
                ca_bar   = ca_bar_al
                cb_bar   = cb_bar_al
                cgam     = cgam_al
                cgam_bar = cgam_bar_al
            endif
            vi3 = 0.0e0_wp
            vi4 = 0.0e0_wp
            vi5 = 0.0e0_wp
            do j=1, natom
                if( i .eq. j ) cycle
                vi3 = vi3+sigma(3,i,j)*b(3,i,j)
                vi4 = vi4+sigma(4,i,j)*b(4,i,j)
                vi5 = vi5+sigma(5,i,j)*b(5,i,j)
                do k=1, natom
                    if( i .eq. k ) cycle
                    if( kap3(id(i),id(j),id(k)) .ne. 0.0e0_wp )then
                        vi3 = vi3+kap3(id(i),id(j),id(k))*b(3,i,j)*b(3,i,k)
                    endif
                    if( kap4(id(i),id(j),id(k)) .ne. 0.0e0_wp )then
                        vi4 = vi4+kap4(id(i),id(j),id(k))*b(4,i,j)*b(4,i,k)
                    endif
                    if( kap5(id(i),id(j),id(k)) .ne. 0.0e0_wp )then
                        vi5 = vi5+kap5(id(i),id(j),id(k))*b(5,i,j)*b(5,i,k)
                    endif
                enddo
            enddo
            add = vi3**ca_bar+cb*vi4**cb_bar
            if( vi5 .ne. 0.0e0_wp ) add=add+cgam*vi5**cgam_bar
            en2 = en2-ca*add
        enddo
        v = en1+en2
    end subroutine alhpot
end module chempotpy_almhn_vbo2_pes

!> @brief Expose AlmHn_VBO2 through the Python-facing PES entry point.
subroutine pes( xyz, nal, nh, natoms, igrad, p, g, d )
    use chempotpy_almhn_vbo2_pes, only: evaluate
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
