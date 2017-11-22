program lorenz63_4dvar
    use params
    use lorenz63, only: run_model
    use utils, only: randn
    use io, only: output
    use assim, only: calc_cost, calc_cost_grad

    implicit none

    real(dp), dimension(tstep,3) :: truth = 0.0_dp
    real(dp), dimension(tstep,3) :: best_guess
    real(dp) :: obs(tstep/freq,3)
    real(dp) :: l(3), f, norm = 100.0_dp, initial(3)
    real(dp) :: time(tstep)
    logical :: converged = .false.
    integer :: i, j = 1

    ! Check whether observation frequency divides into total number of timsteps
    if (mod(tstep, freq) /= 0) then
        stop 'Number of timesteps must be a multiple of the observation frequency'
    end if
    
    ! Generate time axis
    time = (/ (real(i)*h, i = 0, tstep-1) /)

    ! Run truth
    truth = run_model(tstep, (/ 1.0_dp, 1.0_dp, 1.0_dp /))

    ! Calculate observations
    do i = 1, tstep, freq
        last = i
        obs(1+i/freq,:) = (/&
        & randn(truth(i,1), sqrt(obs_var)), &
        & randn(truth(i,2), sqrt(obs_var)), &
        & randn(truth(i,3), sqrt(obs_var)) &
        & /)
    end do

    ! Output truth and observations
    call output(time, truth, "truth.txt")
    call output(time, obs, "obs.txt", freq)

    ! Set initial best guess
    initial = (/ 8.0_dp, 8.0_dp, 8.0_dp /)

    ! Perform minimisation
    do while (.not. converged)
        ! Compute cost of current best guess
        best_guess = run_model(tstep, initial)
        f = calc_cost(tstep, best_guess, obs)
        print *, f

        ! Output first guess
        if (j == 1) then
            call output(time, best_guess, "first_guess.txt")
        end if

        ! Compute gradient of cost function
        l = calc_cost_grad(tstep, best_guess, obs)

        ! Compute norm of cost function gradient
        norm = sqrt(sum(l**2))

        ! Normalise gradient vector
        l = l/norm

        initial = initial - 0.5_dp*l

        ! Max number of iterations reached
        if (j > max_iterations) then
            converged = .true.
            print *, 'Failed to converge in maximum number of iterations'
        end if

        j = j + 1
    end do

    call output(time, best_guess, "final_guess.txt")
end program lorenz63_4dvar
