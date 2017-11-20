module assim
    use params
    use lorenz63, only: a, b, r, run_adjoint

    implicit none

    private
    public calc_cost, calc_cost_grad

contains
    function calc_cost(tstep, traj, obs, D) result(f)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: traj(tstep,3), obs(tstep,3), D(tstep)
        real(dp) :: f
        integer :: i

        ! Calculate cost function
        f = 0.0_dp
        do i = 1, tstep
            f = f + 0.5 * sum(D(i)*(traj(i,:) - obs(i,:))**2)
        end do
    end function calc_cost

    function calc_cost_grad(tstep, traj, obs, D) result(hat)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: traj(tstep,3), obs(tstep,3), D(tstep)
        real(dp) :: hat(3)
        integer :: i

        hat = (traj(tstep,:) - obs(tstep,:))*D(tstep)

        do i = tstep-1, 1, -1
            hat = run_adjoint(1, 0, traj(i,:), hat) + (traj(i,:) - obs(i,:))*D(i)
        end do
    end function calc_cost_grad
end module assim
