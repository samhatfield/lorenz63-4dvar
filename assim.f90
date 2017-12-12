!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains function for performing the 4DVar data assimilation technique.
module assim
    use params
    use lorenz63, only: run_tangent_linear, run_adjoint

    implicit none

    private
    public calc_cost, calc_cost_grad

contains
    !> @brief
    !> Calculates cost function for a given initial perturbation with respect
    !> to the given observations, using the tangent linear model to evolve
    !> this perturbation through the assimilation window.
    !> @param[in] tstep the length of the assimilation window in timesteps
    !> @param[in] traj the nonlinear trajectory
    !> @param[in] del the initial perturbation
    !> @param[in] innov the innovations
    !> @return J the cost function
    function calc_cost(tstep, traj, del, innov) result(J)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: traj(tstep,3), del(3), innov(tstep/freq,3)
        real(dp) :: linear_traj(tstep,3)
        real(dp) :: J
        integer :: i

        ! Compute linear evolution of initial perturbation
        linear_traj = run_tangent_linear(tstep, traj, del)

        ! Calculate cost function
        J = 0.0_dp
        do i = 1, tstep, freq
            J = J + 0.5 * sum((linear_traj(i,:) - innov(1+i/freq,:))**2)/obs_var
        end do
    end function calc_cost

    !> @brief
    !> Calculate gradient of cost function for a given initial perturbation
    !> with respect to the given observations, using the tangent linear model
    !> to evolve the initial perturbation through the assimilation window and
    !> the adjoint model to step the perturbation backwards through time.
    !> @param[in] tstep the length of the assimilation window in timesteps
    !> @param[in] traj the nonlinear trajectory
    !> @param[in] del the initial perturbation
    !> @param[in] innov the innovations
    !> @return hat the gradient of the cost function with respect to the
    !> initial perturbation
    function calc_cost_grad(tstep, traj, del, innov) result(hat)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: traj(tstep,3), del(3), innov(tstep/freq,3)
        real(dp) :: linear_traj(tstep,3)
        real(dp) :: hat(3)
        integer :: i

        ! Compute linear evolution of initial perturbation
        linear_traj = run_tangent_linear(tstep, traj, del)

        ! Calculate first normalised "innovation innovation"
        hat = (linear_traj(last,:) - innov(1+last/freq,:))/obs_var

        ! Step backwards through time, summing each normalised
        ! "innovation innovation" while using the adjoint model to evolve them
        ! backwards at each step
        do i = last-freq, 1, -freq
            hat = run_adjoint(freq+1, traj(i:i+freq,:), hat) &
                & + (linear_traj(i,:) - innov(1+i/freq,:))/obs_var
        end do
    end function calc_cost_grad
end module assim
