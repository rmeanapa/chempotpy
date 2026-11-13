!> @brief Solve a symmetric eigenproblem at the selected working precision.
subroutine dsyev( jobz, uplo, n, a, lda, w, work, lwork, info )
    use chempotpy_kinds, only: wp
    implicit none
    character(len=1), intent(in) :: jobz, uplo
    integer, intent(in)          :: n, lda, lwork
    real(wp), intent(inout)      :: a(lda,*)
    real(wp), intent(out)        :: w(*)
    real(wp), intent(inout)      :: work(*)
    integer, intent(out)         :: info
    integer                      :: i, j, p, q, iteration, max_iterations
    integer                      :: minimum_work, eigen_index
    real(wp)                     :: app, aqq, apq, tau, tangent, cosine, sine
    real(wp)                     :: updated_ip, updated_iq, off_diagonal_max
    real(wp)                     :: scale, tolerance, temporary_value
    real(wp), allocatable        :: matrix(:,:), eigenvectors(:,:)
    info         = 0
    minimum_work = max(1, 3*n-1)
    if( lwork == -1 )then
        work(1) = real(minimum_work, wp)
        return
    endif
    if( jobz /= 'N' .and. jobz /= 'n' .and. &
        jobz /= 'V' .and. jobz /= 'v' ) then
    info = -1
    return
    endif
    if( uplo /= 'U' .and. uplo /= 'u' .and. &
        uplo /= 'L' .and. uplo /= 'l' ) then
    info = -2
    return
    endif
    if( n < 0 )then
        info = -3
        return
    endif
    if( lda < max(1, n) )then
        info = -5
        return
    endif
    if( lwork < minimum_work )then
        info = -8
        return
    endif
    work(1) = real(minimum_work, wp)
    if( n == 0 ) return
    allocate(matrix(n,n), eigenvectors(n,n))
    matrix       = 0.0e0_wp
    eigenvectors = 0.0e0_wp
    do i = 1, n
        eigenvectors(i,i) = 1.0e0_wp
        matrix(i,i)       = a(i,i)
    enddo
    if( uplo == 'U' .or. uplo == 'u' )then
        do j = 2, n
            do i = 1, j-1
                matrix(i,j) = a(i,j)
                matrix(j,i) = a(i,j)
            enddo
        enddo
    else
        do j = 1, n-1
            do i = j+1, n
                matrix(i,j) = a(i,j)
                matrix(j,i) = a(i,j)
            enddo
        enddo
    endif
    scale          = max(1.0e0_wp, maxval(abs(matrix)))
    tolerance      = epsilon(1.0e0_wp)*scale*real(max(1, n), wp)
    max_iterations = max(32, 100*n*n)
    do iteration = 1, max_iterations
        off_diagonal_max = 0.0e0_wp
        p                = 1
        q                = min(2, n)
        do j = 2, n
            do i = 1, j-1
                if( abs(matrix(i,j)) > off_diagonal_max )then
                    off_diagonal_max = abs(matrix(i,j))
                    p                = i
                    q                = j
                endif
            enddo
        enddo
        if( off_diagonal_max <= tolerance ) exit
        app = matrix(p,p)
        aqq = matrix(q,q)
        apq = matrix(p,q)
        tau = (aqq-app)/(2.0e0_wp*apq)
        if( tau >= 0.0e0_wp )then
            tangent = 1.0e0_wp/(tau+hypot(1.0e0_wp, tau))
        else
            tangent = -1.0e0_wp/(-tau+hypot(1.0e0_wp, tau))
        endif
        cosine = 1.0e0_wp/sqrt(1.0e0_wp+tangent*tangent)
        sine   = tangent*cosine
        do i = 1, n
            if( i /= p .and. i /= q )then
                updated_ip  = cosine*matrix(i,p)-sine*matrix(i,q)
                updated_iq  = sine*matrix(i,p)+cosine*matrix(i,q)
                matrix(i,p) = updated_ip
                matrix(p,i) = updated_ip
                matrix(i,q) = updated_iq
                matrix(q,i) = updated_iq
            endif
        enddo
        matrix(p,p) = cosine*cosine*app-2.0e0_wp*sine*cosine*apq +sine*sine*aqq
        matrix(q,q) = sine*sine*app+2.0e0_wp*sine*cosine*apq +cosine*cosine*aqq
        matrix(p,q) = 0.0e0_wp
        matrix(q,p) = 0.0e0_wp
        do i = 1, n
            updated_ip        = cosine*eigenvectors(i,p) -sine*eigenvectors(i,q)
            updated_iq        = sine*eigenvectors(i,p) +cosine*eigenvectors(i,q)
            eigenvectors(i,p) = updated_ip
            eigenvectors(i,q) = updated_iq
        enddo
    enddo
    if( off_diagonal_max > tolerance ) info = 1
    do i = 1, n
        w(i) = matrix(i,i)
    enddo
    do i = 1, n-1
        eigen_index = i
        do j = i+1, n
            if( w(j) < w(eigen_index) ) eigen_index = j
        enddo
        if( eigen_index /= i )then
            temporary_value = w(i)
            w(i)            = w(eigen_index)
            w(eigen_index)  = temporary_value
            do j = 1, n
                temporary_value             = eigenvectors(j,i)
                eigenvectors(j,i)           = eigenvectors(j,eigen_index)
                eigenvectors(j,eigen_index) = temporary_value
            enddo
        endif
    enddo
    if( jobz == 'V' .or. jobz == 'v' )then
        a(1:n,1:n) = eigenvectors
    else
        a(1:n,1:n) = 0.0e0_wp
        do i = 1, n
            a(i,i) = w(i)
        enddo
    endif
    deallocate(matrix, eigenvectors)
end subroutine dsyev
