!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate pes for the AlmHn_VBO2 potential.
subroutine pes( xyz, nal, nh, natoms, igrad, p, g, d )
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
    implicit none
    ! number of electronic state
    integer, parameter            :: nstates = 1
    integer, intent(in)           :: nal, nh, natoms
    double precision, intent(in)  :: xyz(10000,3)
    integer, intent(in)           :: igrad
    double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
    double precision, intent(out) :: d(nstates,nstates,natoms,3)
    double precision              :: cx(natoms), cy(natoms), cz(natoms)
    double precision              :: v
    character(len=2)              :: symbol(natoms)
    integer                       :: istate, iatom
    !initialize
    v = 0.d0; g = 0.d0; d = 0.d0

    do iatom=1, natoms
        cx(iatom) = xyz(iatom,1)/legacy_bohr_to_angstrom
        cy(iatom) = xyz(iatom,2)/legacy_bohr_to_angstrom
        cz(iatom) = xyz(iatom,3)/legacy_bohr_to_angstrom
    end do

    do iatom=1, nal
        symbol(iatom) = 'Al'
    end do
    do iatom=1, nh
        symbol(iatom+nal) = 'H'
    end do

    call alhpot(cx, cy, cz, v, natoms, 10000, symbol)

    if (igrad==0) then
        do istate=1, nstates
            p(istate) = v*legacy_hartree_to_electronvolt
        end do
    else
        write (*, *) 'Only energy is available'
    end if

end subroutine pes

!> @brief Evaluate alhpot for the AlmHn_VBO2 potential.
subroutine alhpot( x, y, z, v, natom, maxatom, symbol )
    !   System:                        AlmHn
    !   Functional form:               Valence--Bond Order
    !   Common name:                   VBO2.beta1
    !   Number of derivatives:         0
    !   Number of bodies:              variable
    !   Number of electronic surfaces: 1
    !   Interface:                     HO-MM-0
    !
    !     Notes:  Valence--Bond Order potential energy function. The functio
    !             form is from Ref. 1.
    !
    !    References: (1) Valence--Bond Order(VBO): A New Approach to Modelin
    !                   Reactive Potential Energy Surfaces for Materials and
    !                   Meiyu Zhao, Mark A. Iron, Przemyslaw Staszewski, Nat
    !                  Rosendo Valero, and Donald G. Truhlar, manuscript in
    !
    !    Units:
    !         energies    - hartrees
    !         coordinates - bohrs
    !
    !    v       --- energy of the system (output)
    !    x,y,z   --- one-dimensional arrays of the Cartesian coordinates
    !                of the system (input)
    !    natom   --- number of atoms in the system (input)
    !    maxatom --- maximal number of atoms (input)
    !    symbol  --- character*2 array of dimension holds the atomic symbols

    implicit none
    double precision :: add, alpha, alpha_al, alpha_alh, alpha_h, beta, beta_al, &
        & beta_alh, beta_h, ca, ca_al, ca_bar, ca_bar_al, ca_bar_h, ca_h, cb, cb_al, cb_bar, &
        & cb_bar_al, cb_bar_h, cb_h, cgam, cgam_al, cgam_bar, cgam_bar_al, cgam_bar_h, cgam_h, del1, &
        & den1, den2, den3, den4, den5, dx, dy, dz, en, en1, en2, fa, gam1_al, gam1_alh, gam1_h, &
        & gam2_al, gam2_alh, gam2_h, gam3_al, gam3_h, gam4_al, gam4_alh, gam4_h, gam5_alh, r, r1, &
        & rij, sig5_alh, v, v1, vi3, vi4, vi5, x, y, z
    integer                     :: i, id, j, k, maxatom, natom
    double precision, parameter :: autoev  = 27.2113961d0
    double precision, parameter :: autoang = 0.529177249d0
    character(len=2)            :: symbol(maxatom)
    dimension x(maxatom), y(maxatom), z(maxatom), r(maxatom,maxatom)
    dimension id(maxatom)
    double precision :: b(5,maxatom,maxatom), sigma(5,maxatom,maxatom)
    double precision :: kap3(2,2,2), kap4(2,2,2), kap5(2,2,2)
    double precision :: kap4_hhal, kap4_halal
    double precision :: cn(5,maxatom,maxatom), gam(5,maxatom,maxatom)
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
        end if
    end do

    ! parameters for kappa.
    ! kap3(4,5) is kappa for p=3,4,5 in Eq. (14)

    kap3(1,1,1) = 0.0d0
    kap3(2,2,2) = 0.0d0
    kap3(1,1,2) = 0.d0
    kap3(1,2,1) = 0.d0
    kap3(2,1,1) = 0.d0
    kap3(1,2,2) = 0.d0
    kap3(2,1,2) = 0.d0
    kap3(2,2,1) = 0.d0
    !
    kap4(1,1,1) = 0.0d0
    kap4(2,2,2) = 0.0d0
    kap4(1,1,2) = kap4_hhal
    kap4(1,2,1) = kap4_hhal
    kap4(2,1,1) = kap4_hhal
    kap4(1,2,2) = kap4_halal
    kap4(2,1,2) = kap4_halal
    kap4(2,2,1) = kap4_halal
    !
    kap5(1,1,1) = 0.0d0
    kap5(2,2,2) = 0.0d0
    kap5(1,1,2) = 0.d0
    kap5(1,2,1) = 0.d0
    kap5(2,1,1) = 0.d0
    kap5(1,2,2) = 0.d0
    kap5(2,1,2) = 0.d0
    kap5(2,2,1) = 0.d0

    !    Calculate all internuclear distances

    do 1 i=1, natom-1
        r(i,i) = 0.d0
        do 1 j=i+1, natom
            dx     = x(i)-x(j)
            dy     = y(i)-y(j)
            dz     = z(i)-z(j)
            r(i,j) = dsqrt(dx*dx+dy*dy+dz*dz)
    1 r(j,i)=r(i,j)
    r(natom,natom) = 0.d0

    !    Calculate en1 which corresponds to the term in brackets in Eq. (11)

    en1 = 0d0
    en2 = 0d0
    do 2 i=1, natom
        v1 = 0d0
        do 3 j=1, natom
            if(i.eq.j)go to 3
            if((id(i).eq.1).and.(id(j).eq.1))then
                del1 = 5.29/autoang

                gam(1,i,j) = gam1_h
                gam(2,i,j) = gam2_h
                gam(3,i,j) = gam3_h
                gam(4,i,j) = gam4_h
                en         = 0.5
                alpha      = alpha_h
                beta       = beta_h

                r1 = 0.7414/autoang

                den1      = 1d0-(r1/del1)**en
                cn(1,i,j) = dexp(gam(1,i,j)/den1)
                den2      = 1d0-(r1/del1)**en
                cn(2,i,j) = dexp(gam(2,i,j)/den2)
                den3      = 1d0-(r1/del1)**en
                cn(3,i,j) = dexp(gam(3,i,j)/den3)
                den4      = 1d0-(r1/del1)**en
                cn(4,i,j) = dexp(gam(4,i,j)/den4)

                sigma(3,i,j) = 1.0
                sigma(4,i,j) = 1.0
                sigma(5,i,j) = 0.0

            else if((id(i).eq.2).and.(id(j).eq.2))then
                del1 = 6.88/autoang

                gam(1,i,j) = gam1_al
                gam(2,i,j) = gam2_al
                gam(3,i,j) = gam3_al
                gam(4,i,j) = gam4_al
                en         = 0.5
                alpha      = alpha_al
                beta       = beta_al

                r1 = 2.7306/autoang

                den1      = 1d0-(r1/del1)**en
                cn(1,i,j) = dexp(gam(1,i,j)/den1)
                den2      = 1d0-(r1/del1)**en
                cn(2,i,j) = dexp(gam(2,i,j)/den2)
                den3      = 1d0-(r1/del1)**en
                cn(3,i,j) = dexp(gam(3,i,j)/den3)
                den4      = 1d0-(r1/del1)**en
                cn(4,i,j) = dexp(gam(4,i,j)/den4)

                sigma(3,i,j) = 1.0
                sigma(4,i,j) = 1.0
                sigma(5,i,j) = 0.0

            else
                del1 = 6.88/autoang

                gam(1,i,j) = gam1_alh
                gam(2,i,j) = gam2_alh
                !          gam(3,i,j)=gam3_alh
                gam(4,i,j) = gam4_alh
                gam(5,i,j) = gam5_alh
                en         = 0.5
                alpha      = alpha_alh
                beta       = beta_alh

                r1 = 1.6637/autoang

                den1      = 1d0-(r1/del1)**en
                cn(1,i,j) = dexp(gam(1,i,j)/den1)
                den2      = 1d0-(r1/del1)**en
                cn(2,i,j) = dexp(gam(2,i,j)/den2)
                !           den3=1d0-(r1/del1)**en
                !           cn(3,i,j)=dexp(gam(3,i,j)/den3)
                den4      = 1d0-(r1/del1)**en
                cn(4,i,j) = dexp(gam(4,i,j)/den4)
                den5      = 1d0-(r1/del1)**en
                cn(5,i,j) = dexp(gam(5,i,j)/den5)

                sigma(3,i,j) = 0.0
                sigma(4,i,j) = 0.0
                sigma(5,i,j) = sig5_alh

            end if

            !   Calculate bij according to Eq. (1)

            rij = r(i,j)
            if(rij.lt.del1)then
                b(1,i,j) = dexp(-gam(1,i,j)/(1d0-(rij/del1)**en))*cn(1,i,j)
                b(2,i,j) = dexp(-gam(2,i,j)/(1d0-(rij/del1)**en))*cn(2,i,j)
                b(3,i,j) = dexp(-gam(3,i,j)/(1d0-(rij/del1)**en))*cn(3,i,j)
                b(4,i,j) = dexp(-gam(4,i,j)/(1d0-(rij/del1)**en))*cn(4,i,j)
                if(id(i).ne.id(j)) then
                    b(3,i,j) = 0.0
                    b(5,i,j) = dexp(-gam(5,i,j)/(1d0-(rij/del1)**en))*cn(5,i,j)
                end if
                fa = alpha*(b(1,i,j)+beta*b(2,i,j))
                v1 = v1+fa
            else
                do k=1, 5
                    b(k,i,j) = 0.0d0
                end do
            end if
        3 continue
        en1 = en1+v1
    2 continue

    !    Calculate en2 which corresponds to the second term in Eq. (11) (pur

    do 4 i=1, natom
        if(id(i).eq.1)then
            ca       = ca_h
            cb       = cb_h
            ca_bar   = ca_bar_h
            cb_bar   = cb_bar_h
            cgam     = cgam_h
            cgam_bar = cgam_bar_h

        elseif(id(i).eq.2)then
            ca       = ca_al
            cb       = cb_al
            ca_bar   = ca_bar_al
            cb_bar   = cb_bar_al
            cgam     = cgam_al
            cgam_bar = cgam_bar_al

        end if

        vi3 = 0.0d0
        vi4 = 0.0d0
        vi5 = 0.0d0
        do 5 j=1, natom
            if(i.eq.j)goto 5
            vi3 = vi3+sigma(3,i,j)*b(3,i,j)
            vi4 = vi4+sigma(4,i,j)*b(4,i,j)
            vi5 = vi5+sigma(5,i,j)*b(5,i,j)

            do 6 k=1, natom
                if(i.eq.k)goto 6
                if(kap3(id(i),id(j),id(k)).ne.0.0d0)then
                    vi3 = vi3+kap3(id(i),id(j),id(k))*b(3,i,j)*b(3,i,k)
                end if
                if(kap4(id(i),id(j),id(k)).ne.0.0d0)then
                    vi4 = vi4+kap4(id(i),id(j),id(k))*b(4,i,j)*b(4,i,k)
                end if
                if(kap5(id(i),id(j),id(k)).ne.0.0d0)then
                    vi5 = vi5+kap5(id(i),id(j),id(k))*b(5,i,j)*b(5,i,k)
                end if
            6 continue
        5 continue
        add = vi3**ca_bar+cb*vi4**cb_bar
        if(vi5.ne.0.0d0)add=add+cgam*vi5**cgam_bar
        en2 = en2-ca*add
    4 continue

    v = en1+en2
end subroutine alhpot
