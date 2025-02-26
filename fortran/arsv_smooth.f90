module arsv_smooth_mod
  use kind_mod, only: dp
  use constants_mod, only: pi
  implicit none
  private
  public :: sv_kalman_smoother

contains

  subroutine sv_kalman_smoother(y, mu, phi, sigma_v, h_smooth, ll)
    ! Given observed returns from a log SV model,
    ! estimate the smoothed log volatility process and compute the log likelihood.
    !
    ! The observation equation is derived from:
    !    s_t = log(y_t^2) = h_t + eta_t,
    ! where eta_t = log(epsilon_t^2) has mean m_eta and variance sigma_eta2
    !
    ! We define:
    !    z_t = s_t - m_eta = h_t + e_t,   e_t ~ N(0, sigma_eta2)
    !
    ! The state equation is:
    !    h_t = mu + phi*(h_{t-1} - mu) + sigma_v*v_t,    v_t ~ N(0,1)
    !
    ! Parameters:
    !   y       : real(dp) array, observed returns.
    !   mu      : real(dp), long-term mean of log volatility.
    !   phi     : real(dp), autoregressive coefficient (|phi|<1).
    !   sigma_v : real(dp), volatility of volatility.
    !
    ! Returns:
    !   h_smooth : real(dp) array, the smoothed estimates of the latent log volatility.
    !   ll       : real(dp), the (approximate) log likelihood.

    real(kind=dp), intent(in) :: y(:), mu, phi, sigma_v
    real(kind=dp), intent(out) :: h_smooth(:), ll
    integer :: T, i
    real(kind=dp), allocatable :: s(:), z(:), h_pred(:), P_pred(:), h_filt(:), P_filt(:)
    real(kind=dp) :: m_eta, sigma_eta2, innov_var, K, A

    ! Get length of input array
    T = size(y)
    allocate(s(T), z(T), h_pred(T), P_pred(T), h_filt(T), P_filt(T))

    ! Transform observations: s_t = log(y_t^2)
    s = log(y * y)

    ! Constants from the log-chi-square distribution of log(epsilon^2)
    m_eta = -1.963510026_dp + log(2.0_dp)  ! psi(0.5) ≈ -1.963510026 + log(2)
    sigma_eta2 = 4.934802201_dp             ! polygamma(1, 0.5) ≈ 4.934802201

    ! Define adjusted observations: z_t = s_t - m_eta
    z = s - m_eta

    ! Initialize arrays for Kalman filtering
    h_pred = 0.0_dp   ! prediction: E[h_t | y_{1:t-1}]
    P_pred = 0.0_dp   ! prediction variance
    h_filt = 0.0_dp   ! filtered estimate: E[h_t | y_{1:t}]
    P_filt = 0.0_dp   ! filtered variance
    ll = 0.0_dp       ! log likelihood accumulator

    ! Initialize state using the stationary distribution
    h_filt(1) = mu
    P_filt(1) = sigma_v**2 / (1.0_dp - phi**2)

    ! Kalman filtering recursion
    do i = 1, T
      if (i == 1) then
        h_pred(i) = h_filt(1)
        P_pred(i) = P_filt(1)
      else
        ! Predict step
        h_pred(i) = mu + phi * (h_filt(i-1) - mu)
        P_pred(i) = phi**2 * P_filt(i-1) + sigma_v**2
      end if

      ! Innovation variance
      innov_var = P_pred(i) + sigma_eta2
      ! Kalman gain
      K = P_pred(i) / innov_var
      ! Update step using observation z(i)
      h_filt(i) = h_pred(i) + K * (z(i) - h_pred(i))
      P_filt(i) = (1.0_dp - K) * P_pred(i)

      ! Contribution to log likelihood (from the prediction error)
      ll = ll - 0.5_dp * (log(2.0_dp * pi) + log(innov_var) + ((z(i) - h_pred(i))**2) / innov_var)
    end do

    ! Rauch-Tung-Striebel smoothing (backward pass)
    h_smooth(T) = h_filt(T)
    do i = T-1, 1, -1
      ! Smoothing gain
      A = P_filt(i) * phi / P_pred(i+1)
      h_smooth(i) = h_filt(i) + A * (h_smooth(i+1) - h_pred(i+1))
    end do

    deallocate(s, z, h_pred, P_pred, h_filt, P_filt)
  end subroutine sv_kalman_smoother

end module arsv_smooth_mod
