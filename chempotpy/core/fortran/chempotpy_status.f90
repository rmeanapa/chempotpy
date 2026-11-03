!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Define status values shared by modern PES interfaces.
module chempotpy_status
    implicit none

    private

    public :: status_derivative_unavailable
    public :: status_invalid_input
    public :: status_success

    integer, parameter :: status_success                = 0
    integer, parameter :: status_derivative_unavailable = 1
    integer, parameter :: status_invalid_input          = 2
end module chempotpy_status
