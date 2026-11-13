!> @brief Define selectable numerical kinds for modern CHEMPOTPY Fortran code.
module chempotpy_kinds
    use iso_fortran_env, only: real32, real64
    implicit none
    private
    integer, parameter, public :: sp = selected_real_kind(6, 37)
    integer, parameter, public :: dp = selected_real_kind(15, 307)
    integer, parameter, public :: qp = selected_real_kind(33, 4931)
    integer, parameter, public :: wp = dp
    public :: real32
    public :: real64
end module chempotpy_kinds
