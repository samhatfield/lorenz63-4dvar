module lorenz63
    use params

    implicit none

    private
    public run_model, run_tangent_linear, run_adjoint, a, b, r

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

    function run_tangent_linear(tstep, in, in_p) result(out_p)
        integer, intent(in) :: tstep
        real(dp), dimension(3), intent(in) :: in, in_p
        real(dp), dimension(tstep, 3) :: out_p
        real(dp), dimension(3) :: k1, k1_p, k2_p
        real(dp) :: x, y, z, x_p, y_p, z_p
        integer :: i

        out_p(1,:) = in_p

        do i = 1, tstep
            x = in(i); x_p = out_p(i,1)
            y = in(i); y_p = out_p(i,2)
            z = in(i); z_p = out_p(i,3)

            k1(1) = h*a*(y - x)
            k1(2) = h*(r*x - y - x*z)
            k1(3) = h*(x*y - b*z)

            k1_p(1) = h*a*(y_p - x_p)
            k2_p(2) = h*(r*x_p - y_p - x_p*z - x*z_p)
            k2_p(3) = h*(x_p*y + x*y_p - b*z_p)

            k2_p(1) = h*a*(y_p + k1_p(2) - x_p - k1_p(1))
            k2_p(2) = h*(r*(x_p + k1_p(1)) - y_p - k1_p(2) &
                & - (x_p + k1_p(1))*(z + k1(3)) - (x + k1(1))*(z_p + k1_p(3)))
            k2_p(3) = h*((x_p + k1_p(1))*(y + k1(2)) + (x + k1(1))*(y_p + k1_p(2)) &
                & - b*(z_p + k1_p(3)))

            out_p(i+1,:) = out_p(i,:) + 0.5_dp * (k1_p + k2_p)
        end do
    end function run_tangent_linear

    function run_adjoint(tstep, in, in_hat) result(hat)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in(3), in_hat(3)
        real(dp) :: hat(3)
        real(dp), dimension(3) :: k1, k1_hat, k2_hat
        integer :: i

        hat = in_hat

        do i = tstep, 1, -1
            k1 = h*(/&
            & a*(in(2) - in(1)), &
            & r*in(1) - in(2) - in(1)*in(3), &
            & in(1)*in(2) - b*in(3) &
            & /)

            k1_hat = 0.5_dp * hat
            k2_hat = 0.5_dp * hat

            hat(1)    = hat(1)    + h*k2_hat(3)*(in(2) + k1(2))
            k1_hat(1) = k1_hat(1) + h*k2_hat(3)*(in(2) + k1(2))
            hat(2)    = hat(2)    + h*k2_hat(3)*(in(1) + k1(1))
            k1_hat(2) = k1_hat(2) + h*k2_hat(3)*(in(1) + k1(1))
            hat(3)    = hat(3)    - h*k2_hat(3)*b
            k1_hat(3) = k1_hat(3) - h*k2_hat(3)*b

            hat(1)    = hat(1)    + h*k2_hat(2)*r
            k1_hat(1) = k1_hat(1) + h*k2_hat(2)*r
            hat(2)    = hat(2)    - h*k2_hat(2)
            k1_hat(2) = k1_hat(2) - h*k2_hat(2)
            hat(1)    = hat(1)    - h*k2_hat(2)*(in(3) + k1(3))
            k1_hat(1) = k1_hat(1) - h*k2_hat(2)*(in(3) + k1(3))
            hat(3)    = hat(3)    - h*k2_hat(2)*(in(1) + k1(1))
            k1_hat(3) = k1_hat(3) - h*k2_hat(2)*(in(1) + k1(1))

            hat(2)    = hat(2)    + h*a*k2_hat(1)
            k1_hat(2) = k1_hat(2) + h*a*k2_hat(1)
            hat(1)    = hat(1)    - h*a*k2_hat(1)
            k1_hat(1) = k1_hat(1) - h*a*k2_hat(1)

            hat(1) = hat(1) + h*k1_hat(3)*in(2) + h*r*k1_hat(2)
            hat(2) = hat(2) + h*in(1)*k1_hat(3) - h*k1_hat(2)
            hat(3) = hat(3) - h*b*k1_hat(3) - h*in(1)*k1_hat(2)
            hat(1) = hat(1) - h*k1_hat(2)*in(3)
            hat(2) = hat(2) + h*a*k1_hat(1)
            hat(1) = hat(1) - h*a*k1_hat(1)
        end do
    end function run_adjoint
end module lorenz63
