!> @brief Define the language-neutral request, result, capability, and metadata contract for PES evaluators.
module chempotpy_contract
    use chempotpy_kinds,  only: wp
    use chempotpy_status, only: chempotpy_status_t, status_derivative_unavailable, status_invalid_input, &
        & status_numerical_error
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: pes_capabilities_t
    public :: pes_metadata_t
    public :: pes_request_t
    public :: pes_result_t
    public :: representation_adiabatic
    public :: representation_diabatic
    integer, parameter :: representation_adiabatic = 1
    integer, parameter :: representation_diabatic  = 2
    integer, parameter :: metadata_name_length     = 128
    integer, parameter :: element_symbol_length    = 3
    !> @brief Describe the derivative and electronic-representation outputs supported by one PES.
    type :: pes_capabilities_t
    logical :: gradient                = .false.
    logical :: nonadiabatic_coupling   = .false.
    logical :: diabatic_matrix         = .false.
    logical :: diabatic_gradient       = .false.
    end type pes_capabilities_t
    !> @brief Store immutable descriptive and dimensional metadata for one PES.
    type :: pes_metadata_t
    character(len=metadata_name_length)               :: system_name  = ''
    character(len=metadata_name_length)               :: surface_name = ''
    character(len=element_symbol_length), allocatable :: atom_symbols(:)
    integer                                           :: state_count = 1
    type(pes_capabilities_t)                              :: capabilities

    contains

    procedure :: atom_count => metadata_atom_count
    procedure :: validate   => validate_metadata
    end type pes_metadata_t
    !> @brief Select the electronic representation and derivative outputs for one PES evaluation.
    type :: pes_request_t
    integer :: representation         = representation_adiabatic
    logical :: gradient               = .false.
    logical :: nonadiabatic_coupling  = .false.

    contains

    procedure :: validate => validate_request
    end type pes_request_t
    !> @brief Own adiabatic or diabatic PES outputs in atomic units.
    type :: pes_result_t
    real(wp), allocatable :: energy(:)
    real(wp), allocatable :: gradient(:,:,:)
    real(wp), allocatable :: nonadiabatic_coupling(:,:,:,:)
    real(wp), allocatable :: diabatic_matrix(:,:)
    real(wp), allocatable :: diabatic_gradient(:,:,:,:)

    contains

    procedure :: clear    => clear_result
    procedure :: validate => validate_result
    end type pes_result_t

    contains

    !> @brief Return the number of atoms described by the metadata.
    pure integer function metadata_atom_count( self )
        class(pes_metadata_t), intent(in) :: self
        metadata_atom_count = 0
        if( allocated(self%atom_symbols) ) metadata_atom_count = size(self%atom_symbols)
    end function metadata_atom_count

    !> @brief Validate PES metadata and report the first contract violation.
    pure subroutine validate_metadata( self, status )
        class(pes_metadata_t), intent(in)     :: self
        type(chempotpy_status_t), intent(inout) :: status
        if( status%failed() ) return
        if( len_trim(self%system_name) == 0 )then
            call status%set(status_invalid_input, 'PES system_name must not be empty')
        elseif( len_trim(self%surface_name) == 0 )then
            call status%set(status_invalid_input, 'PES surface_name must not be empty')
        elseif( .not. allocated(self%atom_symbols) )then
            call status%set(status_invalid_input, 'PES atom_symbols must be allocated')
        elseif( size(self%atom_symbols) == 0 )then
            call status%set(status_invalid_input, 'PES atom_symbols must contain at least one atom')
        elseif( any(len_trim(self%atom_symbols) == 0) )then
            call status%set(status_invalid_input, 'PES atom_symbols contains an empty symbol')
        elseif( self%state_count < 1 )then
            call status%set(status_invalid_input, 'PES state_count must be positive')
        endif
    end subroutine validate_metadata

    !> @brief Validate a PES request against the capabilities declared by its metadata.
    pure subroutine validate_request( self, metadata, status )
        class(pes_request_t), intent(in)      :: self
        type(pes_metadata_t), intent(in)      :: metadata
        type(chempotpy_status_t), intent(inout) :: status
        call metadata%validate(status)
        if( status%failed() ) return
        if( self%representation /= representation_adiabatic .and. &
            & self%representation /= representation_diabatic )then
        call status%set(status_invalid_input, 'unknown PES electronic representation')
    elseif( self%representation == representation_adiabatic )then
        if( self%gradient .and. .not. metadata%capabilities%gradient ) &
            call status%set(status_derivative_unavailable, 'adiabatic gradient is unavailable')
        if( self%nonadiabatic_coupling .and. .not. metadata%capabilities%nonadiabatic_coupling ) &
            call status%set(status_derivative_unavailable, 'nonadiabatic coupling is unavailable')
    else
        if( .not. metadata%capabilities%diabatic_matrix ) &
            call status%set(status_derivative_unavailable, 'diabatic matrix is unavailable')
        if( self%gradient .and. .not. metadata%capabilities%diabatic_gradient ) &
            call status%set(status_derivative_unavailable, 'diabatic gradient is unavailable')
        if( self%nonadiabatic_coupling ) &
            call status%set(status_invalid_input, 'nonadiabatic coupling requires the adiabatic representation')
        endif
    end subroutine validate_request

    !> @brief Release all arrays owned by a PES result.
    pure subroutine clear_result( self )
        class(pes_result_t), intent(inout) :: self
        if( allocated(self%energy) ) deallocate(self%energy)
        if( allocated(self%gradient) ) deallocate(self%gradient)
        if( allocated(self%nonadiabatic_coupling) ) deallocate(self%nonadiabatic_coupling)
        if( allocated(self%diabatic_matrix) ) deallocate(self%diabatic_matrix)
        if( allocated(self%diabatic_gradient) ) deallocate(self%diabatic_gradient)
    end subroutine clear_result

    !> @brief Validate PES result allocation, shape, and finiteness against metadata and a request.
    pure subroutine validate_result( self, metadata, request, status )
        class(pes_result_t), intent(in)       :: self
        type(pes_metadata_t), intent(in)      :: metadata
        type(pes_request_t), intent(in)       :: request
        type(chempotpy_status_t), intent(inout) :: status
        integer :: atom_count, state_count
        call request%validate(metadata, status)
        if( status%failed() ) return
        atom_count  = metadata%atom_count()
        state_count = metadata%state_count
        if( request%representation == representation_adiabatic )then
            if( .not. allocated(self%energy) )then
                call status%set(status_invalid_input, 'adiabatic energy was not returned')
            elseif( any(shape(self%energy) /= [state_count]) )then
                call status%set(status_invalid_input, 'adiabatic energy has the wrong shape')
            elseif( .not. all(ieee_is_finite(self%energy)) )then
                call status%set(status_numerical_error, 'adiabatic energy contains a non-finite value')
            endif
            if( status%failed() ) return
            if( request%gradient ) call validate_adiabatic_gradient(self, state_count, atom_count, status)
            if( status%failed() ) return
            if( request%nonadiabatic_coupling ) &
                call validate_nonadiabatic_coupling(self, state_count, atom_count, status)
        else
            call validate_diabatic_matrix(self, state_count, status)
            if( status%failed() ) return
            if( request%gradient ) call validate_diabatic_gradient(self, state_count, atom_count, status)
        endif
    end subroutine validate_result

    !> @brief Validate the adiabatic-gradient component of a PES result.
    pure subroutine validate_adiabatic_gradient( self, state_count, atom_count, status )
        class(pes_result_t), intent(in)       :: self
        integer, intent(in) :: state_count, atom_count
        type(chempotpy_status_t), intent(inout) :: status
        if( .not. allocated(self%gradient) )then
            call status%set(status_invalid_input, 'adiabatic gradient was not returned')
        elseif( any(shape(self%gradient) /= [state_count, atom_count, 3]) )then
            call status%set(status_invalid_input, 'adiabatic gradient has the wrong shape')
        elseif( .not. all(ieee_is_finite(self%gradient)) )then
            call status%set(status_numerical_error, 'adiabatic gradient contains a non-finite value')
        endif
    end subroutine validate_adiabatic_gradient

    !> @brief Validate the nonadiabatic-coupling component of a PES result.
    pure subroutine validate_nonadiabatic_coupling( self, state_count, atom_count, status )
        class(pes_result_t), intent(in)       :: self
        integer, intent(in) :: state_count, atom_count
        type(chempotpy_status_t), intent(inout) :: status
        if( .not. allocated(self%nonadiabatic_coupling) )then
            call status%set(status_invalid_input, 'nonadiabatic coupling was not returned')
        elseif( any(shape(self%nonadiabatic_coupling) /= [state_count, state_count, atom_count, 3]) )then
            call status%set(status_invalid_input, 'nonadiabatic coupling has the wrong shape')
        elseif( .not. all(ieee_is_finite(self%nonadiabatic_coupling)) )then
            call status%set(status_numerical_error, 'nonadiabatic coupling contains a non-finite value')
        endif
    end subroutine validate_nonadiabatic_coupling

    !> @brief Validate the diabatic-matrix component of a PES result.
    pure subroutine validate_diabatic_matrix( self, state_count, status )
        class(pes_result_t), intent(in)       :: self
        integer, intent(in) :: state_count
        type(chempotpy_status_t), intent(inout) :: status
        if( .not. allocated(self%diabatic_matrix) )then
            call status%set(status_invalid_input, 'diabatic matrix was not returned')
        elseif( any(shape(self%diabatic_matrix) /= [state_count, state_count]) )then
            call status%set(status_invalid_input, 'diabatic matrix has the wrong shape')
        elseif( .not. all(ieee_is_finite(self%diabatic_matrix)) )then
            call status%set(status_numerical_error, 'diabatic matrix contains a non-finite value')
        endif
    end subroutine validate_diabatic_matrix

    !> @brief Validate the diabatic-gradient component of a PES result.
    pure subroutine validate_diabatic_gradient( self, state_count, atom_count, status )
        class(pes_result_t), intent(in)       :: self
        integer, intent(in) :: state_count, atom_count
        type(chempotpy_status_t), intent(inout) :: status
        if( .not. allocated(self%diabatic_gradient) )then
            call status%set(status_invalid_input, 'diabatic gradient was not returned')
        elseif( any(shape(self%diabatic_gradient) /= [state_count, state_count, atom_count, 3]) )then
            call status%set(status_invalid_input, 'diabatic gradient has the wrong shape')
        elseif( .not. all(ieee_is_finite(self%diabatic_gradient)) )then
            call status%set(status_numerical_error, 'diabatic gradient contains a non-finite value')
        endif
    end subroutine validate_diabatic_gradient
end module chempotpy_contract
