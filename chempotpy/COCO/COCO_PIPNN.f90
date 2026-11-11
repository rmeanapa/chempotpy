!! 2020/4/21
!! NN input: PIP only considering the permutation between molecules.
!! 7-45-45-1, fitting to 16000 energies
!! SCF/CVTZ-F12 + Correlation[AE-CCSD(T)-F12a]/CBS[CVDZ-F12,CVTZ-F12]
!! BSSE has not been considered
!! 2020/6/26
!! updated: fitting using exp(-r*ita), and analytical derivatives
!! input: r O1-C1-O2-C2 Bohr
!! output: v (interaction energy only) cm-1

!> @brief Implement COCO_PIPNN, a potential energy surface for COCO.
module chempotpy_coco_pipnn_pes
    implicit none
    private
    public :: evaluate

    contains

    !> @brief Evaluate the COCO_PIPNN potential energy surface for COCO.
    subroutine evaluate( x, igrad, p, g, d )
        use chempotpy_geometry,     only: all_pair_distances_and_jacobian, transform_internal_gradient_to_cartesian
        use chempotpy_legacy_units, only: legacy_bohr_to_angstrom
        implicit none
        ! number of electronic state
        integer, parameter            :: nstates = 1
        integer, parameter            :: natoms  = 4
        integer, intent(in)           :: igrad
        double precision, intent(in)  :: x(natoms,3)
        double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
        double precision, intent(out) :: d(nstates,nstates,natoms,3)
        double precision              :: v
        double precision              :: r_bond(6), dvdr(6), drdx(6,12)
        double precision              :: tx(natoms,3), dvdx(12)
        !initialize
        v  = 0.d0; p = 0.d0; g = 0.d0; d = 0.d0
        tx = x/legacy_bohr_to_angstrom
        call all_pair_distances_and_jacobian(tx, r_bond, drdx)
        call pes_ococ(r_bond, v, igrad, dvdr)
        dvdr = dvdr/8065.7112013/legacy_bohr_to_angstrom
        call transform_internal_gradient_to_cartesian(dvdr, drdx, dvdx)
        if( igrad == 0 )then
            p(1) = v/8065.7112013
        elseif( igrad == 1 )then
            p(1)     = v/8065.7112013
            g(1,:,:) = transpose(reshape(dvdx, [3, natoms]))
        endif
        d = 0.d0
    end subroutine evaluate

    !> @brief Evaluate pes_ococ.
    subroutine pes_ococ( r, v, idv, dvdr )
        implicit none
        double precision, parameter   :: pi = dacos(-1.d0)
        integer, parameter            :: np = 7
        double precision, intent(in)  :: r(6)
        integer, intent(in)           :: idv
        double precision, intent(out) :: v, dvdr(6)
        double precision              :: p(np,1), vp(1), dvdp(np,1), dpdr(np,6)
        if( idv .ne. 0 .and. idv .ne. 1 ) stop "Only energy and gradient are available"
        call pip_ococ(r, p(:,1), idv, dpdr)
        call nnfit_ococ(np, 1, p, vp, idv, dvdp)
        if( idv .eq. 1 ) call dgemv('t', 7, 6, 1.d0, dpdr, 7, dvdp, 1, 0.d0, dvdr, 1)
        v = vp(1)
        return

        contains

        !> @brief Evaluate pip_ococ.
        subroutine pip_ococ( r, p, idv, dpdr )
            ! molecular based fundamental invariants
            implicit none
            double precision, parameter   :: ita = 0.3d0
            double precision, intent(in)  :: r(6)   ! atom order: O1 C1 O2 C2
            integer, intent(in)           :: idv
            double precision, intent(out) :: p(7), dpdr(7,6)
            integer                       :: i
            double precision              :: g(6)
            do i=1, 6
                g(i) = exp(-r(i)*ita)
            enddo
            p(1) = g(1)+g(6)
            p(2) = g(3)+g(4)
            p(3) = g(1)**2+g(6)**2
            p(4) = g(3)**2+g(4)**2
            p(5) = g(1)*g(3)+g(4)*g(6)
            p(6) = g(2)
            p(7) = g(5)
            p(3) = dsqrt(p(3))
            p(4) = dsqrt(p(4))
            p(5) = dsqrt(p(5))
            if( idv .eq. 0 ) return
            dpdr      = 0.d0
            dpdr(1,1) = 1.d0
            dpdr(1,6) = 1.d0
            dpdr(2,3) = 1.d0
            dpdr(2,4) = 1.d0
            dpdr(3,1) = 0.5d0/p(3)*2.d0*g(1)
            dpdr(3,6) = 0.5d0/p(3)*2.d0*g(6)
            dpdr(4,3) = 0.5d0/p(4)*2.d0*g(3)
            dpdr(4,4) = 0.5d0/p(4)*2.d0*g(4)
            dpdr(5,1) = 0.5d0/p(5)*g(3)
            dpdr(5,3) = 0.5d0/p(5)*g(1)
            dpdr(5,4) = 0.5d0/p(5)*g(6)
            dpdr(5,6) = 0.5d0/p(5)*g(4)
            dpdr(6,2) = 1.d0
            dpdr(7,5) = 1.d0
            do i=1, 6
                dpdr(:,i) = dpdr(:,i)*g(i)*(-1.d0*ita)
            enddo
        end subroutine pip_ococ

        !> @brief Evaluate nnfit_ococ.
        subroutine nnfit_ococ( ndim, ntot, r0, vx, idv, dv )
            use chempotpy_coco_pipnn_parameters, only: b1 => hidden_1_biases, b2 => hidden_2_biases, b3 => output_bias, rg => input_ranges, vg => output_range, w1 => input_weights, w2 => hidden_weights, w3 => output_weights
            implicit none
            integer, parameter            :: s0 = 7, s1=45, s2=45
            integer, intent(in)           :: ndim, ntot
            double precision, intent(in)  :: r0(ndim,ntot)
            double precision, intent(out) :: vx(ntot)
            integer, intent(in)           :: idv
            double precision, intent(out) :: dv(ndim,ntot)
            integer                       :: i
            if( ndim .ne. s0 ) stop "ndim .ne. s0"
            if( ntot .ge. 24 )then
                call nsimx(ntot, r0, vx, idv, dv, s0, s1, s2, rg, w1, b1, w2, b2, w3, b3, vg)
            else
                do i=1, ntot
                    call nsim(r0(1,i), vx(i), idv, dv(:,i), s0, s1, s2, rg, w1, b1, w2, b2, w3, b3, vg)
                enddo
            endif
        end subroutine nnfit_ococ

        !> @brief Evaluate nsim.
        subroutine nsim( r0, v, idv, dv, n0, n1, n2, rg, w1, b1, w2, b2, w3, b3, vg )
            ! ---- simulate a neural network with two hidden layers
            ! ----      n0-n1-n2-1
            ! ---- blas routines used in this subroutine: dgemv, ddot
            implicit none
            integer, intent(in)           :: n0, n1, n2, idv
            double precision, intent(in)  :: r0(n0), rg(2,n0), vg(2)
            double precision, intent(in)  :: w1(n0,n1), b1(n1)
            double precision, intent(in)  :: w2(n1,n2), b2(n2)
            double precision, intent(in)  :: w3(n2), b3
            double precision, intent(out) :: v, dv(n0)
            integer                       :: i, j, k
            double precision              :: r(n0), rgg(n0), vgg, ax(n1), bx(n2), maxx(n1), w2maxx(n1,n2), w3mbxx(n2)
            double precision              :: dvtm(n0,n2), rt1(n1), rt2(n2)
            double precision, external    :: ddot
            v   = 0.d0
            r   = r0
            vgg = vg(2)-vg(1)
            ! mapminmax [-1,1]
            do i=1, n0
                rgg(i) = rg(2,i)-rg(1,i)
                r(i)   = 2.d0*(r(i)-rg(1,i))/rgg(i)-1.d0
            enddo
            ! 1st layer
            rt1 = b1
            call dgemv('t', n0, n1, 1.d0, w1, n0, r, 1, 1.d0, rt1, 1)
            ax = dtanh(rt1)
            ! 2nd layer
            rt2 = b2
            call dgemv('t', n1, n2, 1.d0, w2, n1, ax, 1, 1.d0, rt2, 1)
            bx = dtanh(rt2)
            ! output layer
            v = b3+ddot(n2, w3, 1, bx, 1)
            !reverse map
            v = vgg*(v+1.d0)/2+vg(1)
            if( idv .ne. 1 ) return
            ! calculate first derivatives, dv(i)=dv/dr(i)
            ! updated: 2020/6/26, speed up x5
            do j=1, n1
                maxx(j) = 1-ax(j)**2
            enddo
            do k=1, n2
                do j=1, n1
                    w2maxx(j,k) = maxx(j)*w2(j,k)
                enddo
            enddo
            do k=1, n2
                w3mbxx(k) = w3(k)*(1-bx(k)**2)
            enddo
            dv   = 0.d0
            dvtm = 0.d0
            call dgemm('n', 'n', n0, n2, n1, 1.d0, w1, n0, w2maxx, n1, 0.d0, dvtm, n0)
            call dgemv('n', n0, n2, 1.d0, dvtm, n0, w3mbxx, 1, 0.d0, dv, 1)
            do i=1, n0
                dv(i) = dv(i)*vgg/rgg(i)
            enddo
        end subroutine nsim

        !> @brief Evaluate nsimx.
        subroutine nsimx( nt, r, v, idv, dv, n0, n1, n2, rg, w1, b1, w2, b2, w3, b3, vg )
            ! ---- simulate a neural network with two hidden layers
            ! ----    n0-n1(tansig)-n2(tansig)-1(pureline)
            ! ---- blas routines used in this subroutine: dgemm, dgemv
            ! ---- the syntax keeps the same with nsim, except input nt at first
            ! ---- Chen Jun, 2013/11/26
            ! ---- description of parameters
            ! ---- nt  : total number of points in r,v and dv
            ! ---- r   : input of nn, r(n0,nt), r(1:n0,i) is the i-th input
            ! ---- v   : output of nn, v(nt), v(i) is the i-th input
            ! ---- idv : if idv==0, calc. v only, else calc. v and derivatives dv
            ! ---- dv  : derivatives of dv/dr, dv(n0,nt), dv(j,i)=d v(i) / d r(j,i)
            ! ---- n0  : neurons of the input layer, dimension of r
            ! ---- n1  : neurons of the first hidden layer
            ! ---- n2  : neurons of the second hidden layer
            ! ---- rg  : input ranges of the training set, rg(2,n0)
            ! ---- w1  : weights connect input layer with 1-st hidden layer
            ! ---- b1  : biases of 1-st hidden layer
            ! ---- w2  : weights connect 1-st hidden layer with 2-nd hidden layer
            ! ---- b2  : biases of 2-nd hidden layer
            ! ---- w3  : weights connect 2-nd hidden layer with output layer
            ! ---- b3  : biase of output layer
            ! ---- vg  : output range of the training set, vg(2)
            implicit none
            integer, intent(in)           :: nt, idv, n0, n1, n2
            double precision, intent(in)  :: r(n0,nt)
            double precision, intent(out) :: v(nt), dv(n0,nt)
            double precision, intent(in)  :: w1(n0,n1), b1(n1), w2(n1,n2), b2(n2), w3(n2), b3
            double precision, intent(in)  :: rg(2,n0), vg(2)
            double precision              :: x0(n0,nt), x1(nt,n1), x2(nt,n2), tmp, rgg(n0), vgg
            integer                       :: i, j, k, n
            v = 0.d0
            do i=1, n0
                rgg(i) = rg(2,i)-rg(1,i)
            enddo
            vgg = vg(2)-vg(1)
            x0  = r
            do i=1, n0
                x0(i,1:nt) = (x0(i,1:nt)-rg(1,i))/rgg(i)*2.d0-1.d0
            enddo
            do i=1, n1
                x1(1:nt,i) = b1(i)
            enddo
            ! ---- x1=trans(x0(n0,nt))*w1(n0,n1)+x1(nt,n1)
            !call gpu_dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)
            call dgemm('t', 'n', nt, n1, n0, 1.d0, x0, n0, w1, n0, 1.d0, x1, nt)
            x1 = dtanh(x1)
            do i=1, n2
                x2(1:nt,i) = b2(i)
            enddo
            ! ---- x2=x1(nt,n1)*w2(n1,n2)+x2(nt,n2)
            !call gpu_dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)
            call dgemm('n', 'n', nt, n2, n1, 1.d0, x1, nt, w2, n1, 1.d0, x2, nt)
            x2      = dtanh(x2)
            v(1:nt) = b3
            ! ---- v=x2(nt,n2)*w3(n2)+v(nt)
            call dgemv('n', nt, n2, 1.d0, x2, nt, w3, 1, 1.d0, v, 1)
            v = (v+1.d0)/2.d0*vgg+vg(1)
            if( idv .ne. 1 ) return
            ! ---- calculate first derivatives
            dv = 0.d0
            do n=1, nt
                do i=1, n0
                    do k=1, n2
                        tmp = 0.d0
                        do j=1, n1
                            tmp = tmp+w2(j,k)*w1(i,j)*(1-x1(n,j)**2)
                        enddo
                        dv(i,n) = dv(i,n)+w3(k)*tmp*(1-x2(n,k)**2)
                    enddo
                    dv(i,n) = dv(i,n)*vgg/(rgg(i))
                enddo
            enddo
        end subroutine nsimx
    end subroutine pes_ococ

    !> @brief Evaluate potco.
    subroutine potco( r, v, dv )
        ! 2020/4/21  C-O AE-CCSD(T)(short range) MRCI(long range) / ACVTZ Energy in cm-1
        ! updated: 2020/6/25, dv/dr
        implicit none
        integer, parameter          :: ns     = 8
        double precision, parameter :: w1(ns) = [  2.5731334302226245d+00, &
            -4.1635475296365723d+00, &
            6.2652688635816842d+00, &
            -8.1228824615517947d+00, &
            2.5634824563114837d+01, &
            -1.6643878666848999d+00, &
            1.2863481593265147d+01, &
            -5.1395051186435685d+00]
        double precision, parameter :: b1(ns) = [ -1.6421783363943476d+00, &
            1.5364774081764951d+00, &
            -1.1338455332512607d+00, &
            -1.4448051810696709d+00, &
            7.5217573991947644d+00, &
            -1.4005229119446290d+00, &
            1.1053854378210930d+01, &
            -5.9299626180269485d+00]
        double precision, parameter :: w2(ns) = [ -8.8664820626030757d-03, &
            7.8571245473773067d-03, &
            -3.6411047563733342d-03, &
            -4.0358215533209145d-03, &
            9.6640587889626451d-04, &
            -1.4325782866595651d+00, &
            1.2002568907875554d-02, &
            8.3983298757280007d+00]
        double precision, parameter   :: b2 = 6.8970075393140338d+00
        double precision, parameter   :: ra = 1.4000000000000000d+00, rb=7.0000000000000000d+00
        double precision, parameter   :: va = 9.8087308941410573d-02, vb=1.9558422718340193d+05
        double precision, intent(in)  :: r !! bohr
        double precision, intent(out) :: v, dv !! cm-1, cm-1/bohr
        double precision              :: x, y1j
        integer                       :: j
        v  = 0.d0
        dv = 0.d0
        x  = 2.d0*(r-ra)/(rb-ra)-1.d0
        do j=1, ns
            y1j = tanh(b1(j)+w1(j)*x)
            v   = v+w2(j)*y1j
            dv  = dv+w2(j)*(1.d0-y1j**2)*w1(j)
        enddo
        v  = v+b2
        v  = (v+1.d0)*(vb-va)/2.d0+va
        v  = v+0.560096850315234d0
        dv = dv*(vb-va)/(rb-ra)
    end subroutine potco
end module chempotpy_coco_pipnn_pes

!> @brief Expose COCO_PIPNN through the Python-facing PES entry point.
subroutine pes( x, igrad, path, p, g, d )
    use chempotpy_coco_pipnn_pes, only: evaluate
    implicit none
    integer, parameter              :: natoms  = 4
    integer, parameter              :: nstates = 1
    integer, intent(in)             :: igrad
    character(len=1024), intent(in) :: path
    double precision, intent(in)    :: x(natoms,3)
    double precision, intent(out)   :: p(nstates)
    double precision, intent(out)   :: g(nstates,natoms,3)
    double precision, intent(out)   :: d(nstates,nstates,natoms,3)
    call evaluate(x, igrad, p, g, d)
end subroutine pes
