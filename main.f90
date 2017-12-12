!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> An observing system simulation experiment (OSSE) using the Lorenz '63 model
!> and 4DVar.
!> Based on code by Amos Lawless, University of Reading.
program lorenz63_4dvar
    use params
    use lorenz63, only: run_model
    use utils, only: time_seed, randn
    use io, only: output
    use assim, only: calc_cost, calc_cost_grad

    implicit none

    real(dp), dimension(tstep,3) :: truth = 0.0_dp
    real(dp), dimension(tstep,3) :: best_guess
    real(dp), dimension(tstep/freq,3) :: obs, innov
    real(dp) :: diagn(out_iter*in_iter,1)
    real(dp) :: l(3), f, norm, initial(3), del(3)
    real(dp) :: time(tstep)
    integer :: i, j = 1

    ! Check whether observation frequency divides into total number of timsteps
    if (mod(tstep, freq) /= 0) then
        stop 'Number of timesteps must be a multiple of the observation frequency'
    end if

    ! Seed random number generator
    call time_seed

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

    ! Perform outer loop
    do j = 1, out_iter
        ! Generate nonlinear forecast
        best_guess = run_model(tstep, initial)

        ! Output first guess
        if (j == 1) then
            call output(time, best_guess, "first_guess.txt")
        end if

        ! Calculate innovations
        innov = obs - best_guess(1:tstep:freq,:)

        ! Initial increment for inner loop
        del = 0.0_dp

        ! Perform inner loop minimisation
        do i = 1, in_iter
            ! Compute cost of current best guess
            diagn((j-1)*in_iter+i,1) = calc_cost(tstep, best_guess, del, innov)

            ! Compute gradient of cost function
            l = calc_cost_grad(tstep, best_guess, del, innov)

            ! Compute norm of cost function gradient
            norm = sqrt(sum(l**2))

            ! Normalise gradient vector
            l = l/norm

            ! Update initial perturbation estimate at beginning of window
            del = del - 0.5_dp*l
        end do

        initial = initial + del
    end do

    ! Output final best guess
    call output(time, best_guess, "final_guess.txt")

    ! Output diagnostics
    call output((/(real(i,dp), i = 1, out_iter*in_iter)/), diagn, "diagnostics.txt")
end program lorenz63_4dvar
