module io
    use params

    implicit none

    private
    public output

contains
    subroutine output(output_array, filename, stride_in)
        real(dp), intent(in) :: output_array(:,:)
        character(len=*), intent(in) :: filename
        integer, optional :: stride_in
        integer :: i, stride

        if (present(stride_in)) then
            stride = stride_in
        else
            stride = 1
        end if

        open(1, file=filename)
        do i = 1, size(output_array, 1), stride
            write (1,*) output_array(i,:)
        end do
        close(1)
    end subroutine output
end module io
