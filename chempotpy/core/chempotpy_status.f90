!> @brief Define structured, non-terminating status values for PES interfaces.
module chempotpy_status
    implicit none
    private
    public :: chempotpy_status_t
    public :: status_derivative_unavailable
    public :: status_internal_error
    public :: status_invalid_input
    public :: status_model_unavailable
    public :: status_numerical_error
    public :: status_success
    integer, parameter :: status_success                = 0
    integer, parameter :: status_derivative_unavailable = 1
    integer, parameter :: status_invalid_input          = 2
    integer, parameter :: status_internal_error         = 3
    integer, parameter :: status_model_unavailable      = 4
    integer, parameter :: status_numerical_error        = 5
    integer, parameter :: status_message_length         = 256
    !> @brief Carry an evaluator outcome without writing or terminating.
    type :: chempotpy_status_t
    integer                              :: code = status_success
    character(len=status_message_length) :: message = ''

    contains

    procedure :: failed => status_failed
    procedure :: set    => set_status
    end type chempotpy_status_t

    contains

    !> @brief Return true when an evaluation did not complete successfully.
    pure logical function status_failed( self )
        class(chempotpy_status_t), intent(in) :: self
        status_failed = self%code /= status_success
    end function status_failed

    !> @brief Record the first failure while preserving its diagnostic text.
    pure subroutine set_status( self, code, message )
        class(chempotpy_status_t), intent(inout) :: self
        integer, intent(in)          :: code
        character(len=*), intent(in) :: message
        if( self%code == status_success )then
            self%code    = code
            self%message = message
        endif
    end subroutine set_status
end module chempotpy_status
