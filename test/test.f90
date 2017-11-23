program test
    use params
    use lorenz63, only: run_model, run_tangent_linear
    use assim, only: calc_cost, calc_cost_grad

    implicit none

    integer, parameter :: n_samples = 12
    integer, parameter :: nsteps = 100

    call test_tl
contains
    subroutine test_tl
        real(dp), dimension(3) :: pert_orig, orig_state, pert, pert_state, diff
        real(dp), dimension(nsteps,3) :: orig_traj, pert_traj, pert_traj_tl
        real(dp) :: rel_err
        real(dp) :: gamma = 1000_dp
        integer :: i

        pert_orig = (/ 1.0_dp, -1.0_dp, 0.5_dp /)

        print *, 'Tangent linear test'
        print *, '============================================================'
        print *, 'The relative error should tend towards zero as the magnitude'
        print *, 'of the perturbation decreases.'
        print *, '============================================================'

        write (*, '(A10, A16)') 'Magnitude', 'Relative error'
        do i = 1, 9
            gamma = gamma * 0.1_dp

            orig_state = (/ 1.0_dp, 2.0_dp, 1.5_dp /)
            pert = gamma * pert_orig
            pert_state = orig_state + pert

            orig_traj = run_model(nsteps, orig_state)
            pert_traj = run_model(nsteps, pert_state)
            pert_traj_tl = run_tangent_linear(nsteps, orig_traj, pert)

            diff = pert_traj(nsteps,:) - orig_traj(nsteps,:)
            rel_err = 100.0_dp * norm2(diff - pert_traj_tl(nsteps,:)) &
                & / norm2(pert_traj_tl(nsteps,:))

            write (*,'(E10.1, F16.9)') gamma, rel_err
        end do
        print *, '============================================================'
    end subroutine test_tl

    subroutine test_adj

    end subroutine test_adj

    subroutine test_grad
        real(dp), dimension(3) :: initial, grad1, grad2, rand_unit
        real(dp) :: grad1_norm(3), alpha = 1.0_dp, obs(tstep/freq,3), cost1, cost2
        real(dp) :: traj(tstep,3)
        real(dp) :: alpha_store(n_samples)
        integer :: i

        initial = (/ 1.0_dp, 2.0_dp, 5.0_dp /)

        traj = run_model(tstep, initial)

        do i = 1, tstep, freq
            obs(1+i/freq,:) = traj(i,:)
        end do

        ! Generate random unit vector
        call random_number(rand_unit)
        rand_unit = rand_unit / norm2(rand_unit)

        ! Trajectory from random vector
        traj = run_model(tstep, rand_unit)
        grad1 = calc_cost_grad(tstep, traj, obs)
        cost1 = calc_cost(tstep, traj, obs)

        grad1_norm = grad1 / norm2(grad1)

        do i = 1, n_samples
            initial = rand_unit + alpha*grad1_norm

            traj = run_model(tstep, initial)
            cost2 = calc_cost(tstep, traj, obs)

            print *,

            alpha_store(i) = alpha
            print *, alpha, (cost2 - cost1)/(alpha * dot_product(grad1_norm, grad1))

            alpha = alpha * 0.1_dp
        end do
    end subroutine test_grad

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
            write (1,*) time_axis(i*stride), output_array(i,:)
        end do
        close(1)
    end subroutine output
end program test
