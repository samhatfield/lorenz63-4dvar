!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains function for performing the 4DVar data assimilation technique.
module assim
    use params
    use lorenz63, only: run_adjoint

    implicit none

    private
    public calc_cost, calc_cost_grad

contains
    !> @brief
    !> Calculates cost function for a given nonlinear trajectory with respect
    !> to the given observations
    !> @param[in] tstep the length of the assimilation window in timesteps
    !> @param[in] traj the nonlinear trajectory
    !> @param[in] obs the observations
    !> @return J the cost function
    function calc_cost(tstep, traj, obs) result(J)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: traj(tstep,3), obs(tstep/freq,3)
        real(dp) :: J
        integer :: i

        ! Calculate cost function
        J = 0.0_dp
        do i = 1, tstep, freq
            J = J + 0.5 * sum((traj(i,:) - obs(1+i/freq,:))**2)/obs_var
        end do
    end function calc_cost

    !> @brief
    !> Calculate gradient of cost function for a given nonlinear trajectory
    !> with respect to the given observations, using the adjoint model.
    !> @param[in] tstep the length of the assimilation window in timesteps
    !> @param[in] traj the nonlinear trajectory
    !> @param[in] obs the observations
    !> @return hat the gradient of the cost function with respect to the
    !> initial perturbation
    function calc_cost_grad(tstep, traj, obs) result(hat)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: traj(tstep,3), obs(tstep/freq,3)
        real(dp) :: hat(3)
        integer :: i

        ! Calculate first normalised innovation
        hat = (traj(last,:) - obs(1+last/freq,:))/obs_var

        ! Step backwards through time, summing each normalised innovation
        ! while using the adjoint model to evolve them backwards at each step
        do i = last-freq, 1, -freq
            hat = run_adjoint(freq+1, traj(i:i+freq,:), hat) &
                & + (traj(i,:) - obs(1+i/freq,:))/obs_var
        end do
    end function calc_cost_grad
end module assim
