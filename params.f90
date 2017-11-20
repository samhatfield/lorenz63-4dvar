module params
    use, intrinsic :: iso_fortran_env

    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32

    real(dp), parameter :: h = 0.01_dp
    real(dp), parameter :: wind_len = 2.0_dp
    real(dp), parameter :: fore_len = 3.0_dp
    integer, parameter :: fcstep = (wind_len + fore_len)/h
    integer, parameter :: tstep = wind_len/h + 1
    integer, parameter :: freq = 2
    integer, parameter :: n_obs = (tstep - 1) / freq + 1
    real(dp), parameter :: obs_var = 2.5_dp
    integer, parameter :: max_iterations = 500
    real(dp), parameter :: tolerance = 0.5_dp
end module params
