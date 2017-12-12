module io
    use params

    implicit none

    private
    public output

contains
    subroutine output(time_axis, output_array, filename, stride_in)
        real(dp), intent(in) :: time_axis(:), output_array(:,:)
        character(len=*), intent(in) :: filename
        integer, optional :: stride_in
        integer :: i, stride

        if (present(stride_in)) then
            stride = stride_in
        else
            stride = 1
        end if

        open(1, file=filename)
        do i = 1, size(output_array, 1)
            write (1,*) time_axis(1+(i-1)*stride), output_array(i,:)
        end do
        close(1)
    end subroutine output
end module io
