!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Read versioned neural-network metadata used by bundled PES data.
module chempotpy_neural_network_io
    implicit none

    private

    public :: read_network_header

    contains

    !> @brief Read a network header with an optional leading input count.
    subroutine read_network_header( file_unit, expected_input_count, &
            hidden_layer_count, output_count )
        integer, intent(in)  :: file_unit, expected_input_count
        integer, intent(out) :: hidden_layer_count, output_count
        character(len=512)   :: header
        integer              :: input_count, status

        read(file_unit,'(a)') header
        read(header,*,iostat=status) input_count, hidden_layer_count, output_count
        if (status == 0) then
            if (input_count /= expected_input_count) then
                error stop "neural-network input count does not match the PES"
            endif
        else
            read(header,*) hidden_layer_count, output_count
        endif
    end subroutine read_network_header
end module chempotpy_neural_network_io
