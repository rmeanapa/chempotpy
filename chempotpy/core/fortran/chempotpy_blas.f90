!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Return a double-precision vector dot product through the BLAS ABI.
function ddot( element_count, first, first_stride, second, second_stride ) &
        result(product)
    implicit none
    integer, intent(in)          :: element_count, first_stride, second_stride
    double precision, intent(in) :: first(*), second(*)
    double precision             :: product
    integer                      :: element, first_index, second_index
    product = 0.0d0
    if( element_count <= 0 ) return
    first_index  = merge(1, 1-(element_count-1)*first_stride, first_stride > 0)
    second_index = merge(1, 1-(element_count-1)*second_stride, second_stride > 0)
    do element = 1, element_count
        product      = product+first(first_index)*second(second_index)
        first_index  = first_index+first_stride
        second_index = second_index+second_stride
    enddo
end function ddot

!> @brief Multiply a double-precision matrix and vector through the BLAS ABI.
subroutine dgemv( transpose, row_count, column_count, alpha, matrix, &
        leading_dimension, input, input_stride, beta, output, output_stride )
    implicit none
    character(len=1), intent(in)    :: transpose
    integer, intent(in)             :: row_count, column_count, leading_dimension
    integer, intent(in)             :: input_stride, output_stride
    double precision, intent(in)    :: alpha, beta
    double precision, intent(in)    :: matrix(leading_dimension,*)
    double precision, intent(in)    :: input(*)
    double precision, intent(inout) :: output(*)
    integer                         :: column, input_index, output_index, row
    integer                         :: input_length, output_length, start_input, start_output
    double precision                :: accumulated
    logical                         :: use_transpose
    use_transpose = transpose == 'T' .or. transpose == 't' .or. transpose == 'C' .or. transpose == 'c'
    input_length  = merge(row_count, column_count, use_transpose)
    output_length = merge(column_count, row_count, use_transpose)
    start_input   = merge(1, 1-(input_length-1)*input_stride, input_stride > 0)
    start_output  = merge(1, 1-(output_length-1)*output_stride, output_stride > 0)
    output_index  = start_output
    do row = 1, output_length
        output(output_index) = beta*output(output_index)
        output_index         = output_index+output_stride
    enddo
    if( alpha == 0.0d0 ) return
    if( use_transpose )then
        output_index = start_output
        do column = 1, column_count
            accumulated = 0.0d0
            input_index = start_input
            do row = 1, row_count
                accumulated = accumulated+matrix(row,column)*input(input_index)
                input_index = input_index+input_stride
            enddo
            output(output_index) = output(output_index)+alpha*accumulated
            output_index         = output_index+output_stride
        enddo
    else
        output_index = start_output
        do row = 1, row_count
            accumulated = 0.0d0
            input_index = start_input
            do column = 1, column_count
                accumulated = accumulated+matrix(row,column)*input(input_index)
                input_index = input_index+input_stride
            enddo
            output(output_index) = output(output_index)+alpha*accumulated
            output_index         = output_index+output_stride
        enddo
    endif
end subroutine dgemv

!> @brief Multiply two double-precision matrices through the BLAS ABI.
subroutine dgemm( transpose_left, transpose_right, row_count, column_count, &
        shared_count, alpha, left, leading_left, right, leading_right, beta, &
        output, leading_output )
    implicit none
    character(len=1), intent(in)    :: transpose_left, transpose_right
    integer, intent(in)             :: row_count, column_count, shared_count
    integer, intent(in)             :: leading_left, leading_right, leading_output
    double precision, intent(in)    :: alpha, beta
    double precision, intent(in)    :: left(leading_left,*), right(leading_right,*)
    double precision, intent(inout) :: output(leading_output,*)
    integer                         :: column, row, shared
    double precision                :: accumulated, left_value, right_value
    logical                         :: left_transposed, right_transposed
    left_transposed  = transpose_left == 'T' .or. transpose_left == 't' .or. transpose_left == 'C' .or. transpose_left == 'c'
    right_transposed = transpose_right == 'T' .or. transpose_right == 't' .or. transpose_right == 'C' .or. transpose_right == 'c'
    do column = 1, column_count
        do row = 1, row_count
            accumulated = 0.0d0
            do shared = 1, shared_count
                if( left_transposed )then
                    left_value = left(shared,row)
                else
                    left_value = left(row,shared)
                endif
                if( right_transposed )then
                    right_value = right(column,shared)
                else
                    right_value = right(shared,column)
                endif
                accumulated = accumulated+left_value*right_value
            enddo
            output(row,column) = alpha*accumulated+beta*output(row,column)
        enddo
    enddo
end subroutine dgemm
