module utils
    use params

    implicit none

    private
    public randn

contains
        function randn(mean, stdev)
            real(dp), intent(in) :: mean, stdev
            real(dp) :: u, v, randn
            real(dp) :: rand(2)

            call random_number(rand)

            ! Box-Muller method
            u = (-2.0_dp * log(rand(1))) ** 0.5_dp
            v =   2.0_dp * 6.28318530718_dp * rand(2)
            randn = mean + stdev * u * sin(v)
        end function randn
end module
