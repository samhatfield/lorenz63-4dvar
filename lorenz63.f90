module lorenz63
    use params

    implicit none

    private
    public run_model, run_adjoint, a, b, r

    real(dp), parameter :: a = 10.0_dp
    real(dp), parameter :: b = 8.0_dp/3.0_dp
    real(dp), parameter :: r = 28.0_dp

contains
    function run_model(tstep, in) result(out)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in(3)
        real(dp), dimension(tstep, 3) :: out
        real(dp), dimension(3) :: k1, k2
        integer :: i

        out(1,:) = in

        do i = 1, tstep-1
            k1 = h*(/&
            & a*(out(i,2) - out(i,1)), &
            & r*out(i,1) - out(i,2) - out(i,1)*out(i,3), &
            & out(i,1)*out(i,2) - b*out(i,3) &
            & /)

            k2 = h*(/&
            & a*(out(i,2) + k1(2) - out(i,1) - k1(1)), &
            & r*(out(i,1) + k1(1)) - out(i,2) - out(i,2) - (out(i,1) + k1(1))*(out(i,3) + k1(3)), &
            & (out(i,1) + k1(1))*(out(i,2) + k1(2)) - b*(out(i,3) + k1(3)) &
            & /)

            out(i+1,:) = out(i,:) + 0.5_dp * (k1 + k2)
        end do
    end function run_model

    function run_adjoint(tstep, in, in_hat) result(hat)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in(tstep,3), in_hat(3)
        real(dp) :: hat(3)
        real(dp), dimension(3) :: k1, k1_hat, k2_hat
        integer :: i

        hat = in_hat

        do i = tstep, 1, -1
            k1 = h*(/&
            & a*(in(i,2) - in(i,1)), &
            & r*in(i,1) - in(i,2) - in(i,1)*in(i,3), &
            & in(i,1)*in(i,2) - b*in(i,3) &
            & /)

            k1_hat = 0.5_dp * hat
            k2_hat = 0.5_dp * hat

            hat(1)    = hat(1)    + h*k2_hat(3)*(in(i,2) + k1(2))
            k1_hat(1) = k1_hat(1) + h*k2_hat(3)*(in(i,2) + k1(2))
            hat(2)    = hat(2)    + h*k2_hat(3)*(in(i,1) + k1(1))
            k1_hat(2) = k1_hat(2) + h*k2_hat(3)*(in(i,1) + k1(1))
            hat(3)    = hat(3)    - h*k2_hat(3)*b
            k1_hat(3) = k1_hat(3) - h*k2_hat(3)*b

            hat(1)    = hat(1)    + h*k2_hat(2)*r
            k1_hat(1) = k1_hat(1) + h*k2_hat(2)*r
            hat(2)    = hat(2)    - h*k2_hat(2)
            k1_hat(2) = k1_hat(2) - h*k2_hat(2)
            hat(1)    = hat(1)    - h*k2_hat(2)*(in(i,3) + k1(3))
            k1_hat(1) = k1_hat(1) - h*k2_hat(2)*(in(i,3) + k1(3))
            hat(3)    = hat(3)    - h*(in(i,1) + k1(1)) * k2_hat(2)
            k1_hat(3) = k1_hat(3) - h*(in(i,1) + k1(1)) * k2_hat(2)

            hat(2)    = hat(2)    + h*a*k2_hat(1)
            k1_hat(2) = k1_hat(2) + h*a*k2_hat(1)
            hat(1)    = hat(1)    - h*a*k2_hat(1)
            k1_hat(1) = k1_hat(1) - h*a*k2_hat(1)

            hat(1) = hat(1) + h*k1_hat(3)*in(i,2) + h*r*k1_hat(2)
            hat(2) = hat(2) + h*in(i,1)*k1_hat(3) - h*k1_hat(2)
            hat(3) = hat(3) - h*b*k1_hat(3) - h*in(i,1)*k1_hat(2)
            hat(1) = hat(1) - h*k1_hat(2)*in(i,3)
            hat(2) = hat(2) + h*a*k1_hat(1)
            hat(1) = hat(1) - h*a*k1_hat(1)
        end do
    end function run_adjoint
end module lorenz63
