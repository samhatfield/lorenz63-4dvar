module assim
    use params
    use lorenz63, only: a, b, r, run_adjoint

    implicit none

    private
    public calc_cost, calc_cost_grad

contains
    function calc_cost(tstep, traj, obs) result(J)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: traj(tstep,3), obs(tstep,3)
        real(dp) :: J
        integer :: i

        ! Calculate cost function
        J = 0.0_dp
        do i = 1, tstep
            J = J + 0.5 * sum((traj(i,:) - obs(i,:))**2)/obs_var
        end do
    end function calc_cost

    function calc_cost_grad(tstep, traj, obs) result(hat)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: traj(tstep,3), obs(tstep,3)
        real(dp) :: hat(3)
        integer :: i

        hat = (traj(tstep,:) - obs(tstep,:))/obs_var

        do i = tstep-1, 2, -1
            hat = run_adjoint(1, traj(i,:), hat) + (traj(i,:) - obs(i,:))/obs_var
        end do

        hat = hat + (traj(1,:) - obs(1,:))/obs_var
    end function calc_cost_grad
end module assim
