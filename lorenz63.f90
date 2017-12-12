!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains functions for integrating the Lorenz '63 system, along with its
!> tangent linear and adjoint models.
module lorenz63
    use params

    implicit none

    private
    public run_model, run_tangent_linear, run_adjoint

    ! Model parameters
    real(dp), parameter :: a = 10.0_dp
    real(dp), parameter :: b = 8.0_dp/3.0_dp
    real(dp), parameter :: r = 28.0_dp

contains
    !> @brief
    !> The full nonlinear ODE.
    !> @param[in] state the state vector at which to evaluate the ODE
    !> @return dRdT the evaluated ODE
    function dRdT(state)
        real(dp), intent(in) :: state(3)
        real(dp) :: dRdT(3), x, y, z

        ! Extract vector elements
        x = state(1); y = state(2); z = state(3)

        dRdT = (/ a*(y - x), r*x - y - x*z, x*y - b*z /)
    end function dRdT

    !> @brief
    !> The Jacobian of the ODE multiplied by the given perturbation.
    !> @param[in] state the state vector at which to evaluate the Jacobian
    !> @param[in] dstate the perturbation vector
    !> @return jacob the evaluated Jacobian multiplied by the given
    !> perturbation
    function jacob(state, dstate)
        real(dp), intent(in) :: state(3), dstate(3)
        real(dp) :: jacob(3), x, y, z, dx, dy, dz

        ! Extract vector elements
        x = state(1); y = state(2); z = state(3)
        dx = dstate(1); dy = dstate(2); dz = dstate(3)

        jacob = (/ a*(dy - dx), r*dx - dy - dx*z - x*dz, dx*y + x*dy - b*dz /)
    end function jacob

    !> @brief
    !> The adjoint of the Jacobian of the ODE multiplied by the given
    !> perturbation.
    !> @param[in] state the state vector at which to evaluate the adjoint of
    !> the Jacobian
    !> @param[in] dstate_a the perturbation vector
    !> @return jacob the evaluated adjoint of the Jacobian multiplied by the
    !> given perturbation
    function jacob_a(state, dstate_a)
        real(dp), intent(in) :: state(3), dstate_a(3)
        real(dp) :: jacob_a(3), x, y, z, dx_a, dy_a, dz_a

        ! Extract vector elements
        x = state(1); y = state(2); z = state(3)
        dx_a = dstate_a(1); dy_a = dstate_a(2); dz_a = dstate_a(3)

        jacob_a = (/ -a*dx_a + (r-z)*dy_a + y*dz_a, a*dx_a - dy_a + x*dz_a, -x*dy_a - b*dz_a /)
    end function jacob_a

    !> @brief
    !> Run full nonlinear model using modified Euler scheme.
    !> @param[in] tstep the length of the integration
    !> @param[in] in the initial condition
    !> @return out the computed trajectory
    function run_model(tstep, in) result(out)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in(3)
        real(dp), dimension(tstep, 3) :: out
        real(dp), dimension(3) :: k1, k2
        integer :: i

        ! First step
        out(1,:) = in

        do i = 1, tstep-1
            k1 = h*dRdT(out(i,:))
            k2 = h*dRdT(out(i,:) + k1)

            out(i+1,:) = out(i,:) + 0.5_dp * (k1 + k2)
        end do
    end function run_model

    !> @brief
    !> Run tangent linear model.
    !> @param[in] tstep the length of the integration
    !> @param[in] in the trajectory around which to perform the integration
    !> @param[in] din the initial condition
    !> @return dout the computed trajectory
    function run_tangent_linear(tstep, in, din) result(dout)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in(tstep, 3), din(3)
        real(dp) :: dout(tstep,3)
        real(dp), dimension(3) :: k1, dk1, dk2
        integer :: i

        ! First step
        dout(1,:) = din

        do i = 1, tstep-1
            k1 = h*dRdT(in(i,:))
            dk1 = h*jacob(in(i,:), dout(i,:))
            dk2 = h*jacob(in(i,:) + k1, dout(i,:) + dk1)

            dout(i+1,:) = dout(i,:) + 0.5_dp * (dk1 + dk2)
        end do
    end function run_tangent_linear

    !> @brief
    !> Run adjoint model.
    !> @param[in] tstep the length of the integration
    !> @param[in] in the trajectory around which to perform the integration
    !> @param[in] din_a the initial condition
    !> @return dout_a the computed trajectory
    function run_adjoint(tstep, in, din_a) result(dout_a)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in(:,:), din_a(3)
        real(dp) :: dout_a(3)
        real(dp), dimension(3) :: k1, dk1_a, dk2_a
        integer :: i

        ! First step
        dout_a = din_a

        do i = tstep-1, 1, -1
            k1 = h*dRdT(in(i,:))
            dk1_a = h*jacob_a(in(i,:) + k1, dout_a)
            dk2_a = h*jacob_a(in(i,:), dout_a + dk1_a)

            dout_a = dout_a + 0.5_dp * (dk1_a + dk2_a)
        end do
    end function run_adjoint
end module lorenz63
