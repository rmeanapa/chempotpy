!> @brief Provide working-precision linear-algebra operations used by PES kernels.
module chempotpy_linear_algebra
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: symmetric_eigendecomposition
    interface
        !> @brief Solve a symmetric eigenproblem with the bundled DSYEV-compatible kernel.
        subroutine dsyev( jobz, uplo, n, a, lda, w, work, lwork, info )
            import wp
            implicit none
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in)          :: n, lda, lwork
            real(wp), intent(inout)      :: a(lda,*)
            real(wp), intent(out)        :: w(*)
            real(wp), intent(inout)      :: work(*)
            integer, intent(out)         :: info
        end subroutine dsyev
    end interface

    contains

    !> @brief Return the ordered eigenvalues and optional eigenvectors of a real symmetric matrix.
    subroutine symmetric_eigendecomposition( matrix, eigenvalues, eigenvectors, info )
        real(wp), intent(in)            :: matrix(:,:)
        real(wp), intent(out)           :: eigenvalues(size(matrix, 1))
        real(wp), intent(out), optional :: eigenvectors(size(matrix, 1),size(matrix, 2))
        integer, intent(out), optional  :: info
        character(len=1)                :: job
        integer                         :: solver_info, work_size
        real(wp)                        :: result_matrix(size(matrix, 1),size(matrix, 2))
        real(wp), allocatable           :: work(:)
        solver_info = 0
        eigenvalues = 0.0_wp
        if( size(matrix, 1) /= size(matrix, 2) )then
            solver_info = -1
            if( present(eigenvectors) ) eigenvectors = 0.0_wp
            if( present(info) ) info = solver_info
            return
        endif
        if( size(matrix, 1) == 0 )then
            if( present(eigenvectors) ) eigenvectors = matrix
            if( present(info) ) info = solver_info
            return
        endif
        result_matrix = matrix
        job           = merge('V', 'N', present(eigenvectors))
        allocate(work(1))
        call dsyev(job, 'L', size(matrix, 1), result_matrix, size(matrix, 1), &
            & eigenvalues, work, -1, solver_info)
        if( solver_info == 0 )then
            work_size = max(1, int(work(1)))
            deallocate(work)
            allocate(work(work_size))
            call dsyev(job, 'L', size(matrix, 1), result_matrix, size(matrix, 1), &
                & eigenvalues, work, work_size, solver_info)
        endif
        if( present(eigenvectors) ) eigenvectors = result_matrix
        if( present(info) ) info = solver_info
    end subroutine symmetric_eigendecomposition
end module chempotpy_linear_algebra
