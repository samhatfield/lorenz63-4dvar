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
        real(dp) :: x, y, z
        integer :: i

        out(1,:) = in

        do i = 1, tstep-1
            x = out(i,1); y = out(i,2); z = out(i,3)
            k1 = h*(/ a*(y - x), r*x - y - x*z, x*y - b*z /)

            k2 = h*(/&
            & a*(y + k1(2) - x - k1(1)), &
            & r*(x + k1(1)) - y - k1(2) - (x + k1(1))*(z + k1(3)), &
            & (x + k1(1))*(y + k1(2)) - b*(z + k1(3)) &
            & /)

            out(i+1,:) = out(i,:) + 0.5_dp * (k1 + k2)
        end do
    end function run_model

    function run_tangent_linear(tstep, in, in_p) result(out_p)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in(tstep, 3), in_p(3)
        real(dp) :: out_p(tstep,3)
        real(dp), dimension(3) :: k1, k1_p, k2_p
        real(dp) :: x, y, z, x_p, y_p, z_p
        integer :: i

        out_p(1,:) = in_p

        do i = 1, tstep-1
            x = in(i,1); x_p = out_p(i,1)
            y = in(i,2); y_p = out_p(i,2)
            z = in(i,3); z_p = out_p(i,3)

            k1 = h*(/ a*(y - x), r*x - y - x*z, x*y - b*z /)

            k1_p(1) = h*a*(y_p - x_p)
            k1_p(2) = h*(r*x_p - y_p - x_p*z - x*z_p)
            k1_p(3) = h*(x_p*y + x*y_p - b*z_p)

            k2_p(1) = h*a*(y_p + k1_p(2) - x_p - k1_p(1))
            k2_p(2) = h*(r*(x_p + k1_p(1)) - y_p - k1_p(2) &
                & - (x_p + k1_p(1))*(z + k1(3)) - (x + k1(1))*(z_p + k1_p(3)))
            k2_p(3) = h*((x_p + k1_p(1))*(y + k1(2)) + (x + k1(1))*(y_p + k1_p(2)) &
                & - b*(z_p + k1_p(3)))

            out_p(i+1,:) = out_p(i,:) + 0.5_dp * (k1_p + k2_p)
        end do
    end function run_tangent_linear

    function run_adjoint(in, in_hat) result(hat)
        real(dp), intent(in) :: in(:,:), in_hat(3)
        real(dp) :: hat(3), x, y, z
        real(dp), dimension(3) :: k1, k1_hat, k2_hat
        integer :: i

        hat = in_hat

        do i = size(in,1), 1, -1
            x = in(i,1); y = in(i,2); z = in(i,3)

            k1 = h*(/ a*(y - x), r*x - y - x*z, x*y - b*z /)

            k1_hat = 0.5_dp * hat
            k2_hat = 0.5_dp * hat

            hat(1)    = hat(1)    + h*k2_hat(3)*(y + k1(2))
            k1_hat(1) = k1_hat(1) + h*k2_hat(3)*(y + k1(2))
            hat(2)    = hat(2)    + h*k2_hat(3)*(x + k1(1))
            k1_hat(2) = k1_hat(2) + h*k2_hat(3)*(x + k1(1))
            hat(3)    = hat(3)    - h*k2_hat(3)*b
            k1_hat(3) = k1_hat(3) - h*k2_hat(3)*b

            hat(1)    = hat(1)    + h*k2_hat(2)*r
            k1_hat(1) = k1_hat(1) + h*k2_hat(2)*r
            hat(2)    = hat(2)    - h*k2_hat(2)
            k1_hat(2) = k1_hat(2) - h*k2_hat(2)
            hat(1)    = hat(1)    - h*k2_hat(2)*(z + k1(3))
            k1_hat(1) = k1_hat(1) - h*k2_hat(2)*(z + k1(3))
            hat(3)    = hat(3)    - h*k2_hat(2)*(x + k1(1))
            k1_hat(3) = k1_hat(3) - h*k2_hat(2)*(x + k1(1))

            hat(2)    = hat(2)    + h*a*k2_hat(1)
            k1_hat(2) = k1_hat(2) + h*a*k2_hat(1)
            hat(1)    = hat(1)    - h*a*k2_hat(1)
            k1_hat(1) = k1_hat(1) - h*a*k2_hat(1)

            hat(1) = hat(1) + h*y*k1_hat(3)
            hat(2) = hat(2) + h*x*k1_hat(3)
            hat(3) = hat(3) - h*b*k1_hat(3)
            hat(1) = hat(1) + h*r*k1_hat(2)
            hat(2) = hat(2) - h*k1_hat(2)
            hat(1) = hat(1) - h*z*k1_hat(2)
            hat(3) = hat(3) - h*x*k1_hat(2)
            hat(2) = hat(2) + h*a*k1_hat(1)
            hat(1) = hat(1) - h*a*k1_hat(1)
        end do
    end function run_adjoint
end module lorenz63
