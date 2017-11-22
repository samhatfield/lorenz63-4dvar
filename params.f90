module params
    use, intrinsic :: iso_fortran_env

    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32

    real(dp), parameter :: h = 0.01_dp
    real(dp), parameter :: wind_len = 3.0_dp
    integer, parameter :: tstep = wind_len/h
    integer, parameter :: freq = 6
    integer, parameter :: n_obs = (tstep - 1) / freq + 1
    real(dp), parameter :: obs_var = 2.5_dp
    integer, parameter :: max_iterations = 500
    real(dp), parameter :: tolerance = 0.5_dp
    integer :: last
end module params
