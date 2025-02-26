module arsv_module
  use kind_mod, only: dp
  use constants_mod, only: pi
  use random_mod, only: random_normal
  use statistics_mod, only: mean, variance
  implicit none
  private
  public :: arsv_returns_zero_corr, arsv_vol_and_returns_zero_corr, fit_arsv_zero_corr, &
            arsv_returns, arsv_vol_and_returns

contains

  ! Simulate daily returns with autoregressive stochastic log volatility with zero correlation between
  ! return and log volatility innovations.
  !
  ! Parameters:
  ! - mu: real(dp), the long-term mean of the log volatility process.
  ! - phi: real(dp), the autoregressive parameter (|phi| < 1).
  ! - sigma_v: real(dp), the standard deviation of the volatility shocks.
  ! - T: integer, the number of time periods to simulate.
  ! - burn_in: integer, the number of initial observations to discard to mitigate the effect of initial values (default: 1000).
  !
  ! Returns:
  ! - y: real(dp) array of length T, the simulated daily returns.
  subroutine arsv_returns_zero_corr(mu, phi, sigma_v, T, y, burn_in)
    real(kind=dp), intent(in) :: mu, phi, sigma_v
    integer, intent(in) :: T, burn_in
    real(kind=dp), intent(out) :: y(T)
    integer :: total_length, i
    real(kind=dp), allocatable :: h(:), epsilon(:), v(:)

    total_length = T + burn_in
    allocate(h(total_length), epsilon(total_length), v(total_length))

    call random_seed()
    do i = 1, total_length
      epsilon(i) = random_normal()
      v(i) = random_normal()
    end do

    h(1) = mu
    do i = 2, total_length
      h(i) = mu + phi * (h(i-1) - mu) + sigma_v * v(i)
    end do

    do i = 1, total_length
      h(i) = exp(h(i) / 2.0_dp) * epsilon(i)
    end do

    y = h(burn_in + 1:total_length)
    deallocate(h, epsilon, v)
  end subroutine arsv_returns_zero_corr

  ! Simulate daily returns and conditional standard deviations with autoregressive stochastic log volatility
  ! with zero correlation between return and log volatility innovations.
  !
  ! Parameters:
  ! - mu: real(dp), the long-term mean of the log volatility process.
  ! - phi: real(dp), the autoregressive parameter (|phi| < 1).
  ! - sigma_v: real(dp), the standard deviation of the volatility shocks.
  ! - T: integer, the number of time periods to simulate.
  ! - burn_in: integer, the number of initial observations to discard to mitigate the effect of initial values (default: 1000).
  !
  ! Returns:
  ! - result: real(dp) array of shape (T, 2), where the first column contains the conditional standard deviations
  !           and the second column contains the simulated daily returns.
  subroutine arsv_vol_and_returns_zero_corr(mu, phi, sigma_v, T, result, burn_in)
    real(kind=dp), intent(in) :: mu, phi, sigma_v
    integer, intent(in) :: T, burn_in
    real(kind=dp), intent(out) :: result(T, 2)
    integer :: total_length, i
    real(kind=dp), allocatable :: h(:), sigma_t(:), epsilon(:), v(:)

    total_length = T + burn_in
    allocate(h(total_length), sigma_t(total_length), epsilon(total_length), v(total_length))

    call random_seed()
    do i = 1, total_length
      epsilon(i) = random_normal()
      v(i) = random_normal()
    end do

    h(1) = mu
    do i = 2, total_length
      h(i) = mu + phi * (h(i-1) - mu) + sigma_v * v(i)
    end do

    do i = 1, total_length
      sigma_t(i) = exp(h(i) / 2.0_dp)
      h(i) = sigma_t(i) * epsilon(i)
    end do

    result(:, 1) = sigma_t(burn_in + 1:total_length)
    result(:, 2) = h(burn_in + 1:total_length)
    deallocate(h, sigma_t, epsilon, v)
  end subroutine arsv_vol_and_returns_zero_corr

  ! Estimate the parameters mu, phi, and sigma_v of an ARSV model given observed returns y_t
  ! with zero correlation between return and log volatility innovations.
  !
  ! Parameters:
  ! - y: real(dp) array, the observed returns.
  ! - acf_1_var_max: real(dp), optional maximum autocorrelation variance (default: 0.99).
  !
  ! Returns:
  ! - mu_hat: real(dp), estimated mu.
  ! - phi_hat: real(dp), estimated phi.
  ! - sigma_v_hat: real(dp), estimated sigma_v.
  subroutine fit_arsv_zero_corr(y, mu_hat, phi_hat, sigma_v_hat, acf_1_var_max)
    real(kind=dp), intent(in) :: y(:)
    real(kind=dp), intent(out) :: mu_hat, phi_hat, sigma_v_hat
    real(kind=dp), intent(in), optional :: acf_1_var_max
    real(kind=dp) :: acf_max, mean_eta, var_eta, Var_s, rho_1, sigma_h2_hat, sigma_v2_hat
    real(kind=dp), allocatable :: s_t(:), s_t_star(:), s_t_star_centered(:)
    integer :: n
    real(kind=dp) :: autocov, Var_s_empirical

    acf_max = 0.99_dp
    if (present(acf_1_var_max)) acf_max = acf_1_var_max

    n = size(y)
    allocate(s_t(n), s_t_star(n), s_t_star_centered(n))

    s_t = log(y * y)
    mean_eta = -1.963510026_dp + log(2.0_dp)  ! psi(0.5) + log(2)
    var_eta = 4.934802201_dp                   ! polygamma(1, 0.5)
    s_t_star = s_t - mean_eta

    s_t_star_centered = s_t_star - mean(s_t_star)
    Var_s = variance(s_t_star)

    autocov = sum(s_t_star_centered(2:n) * s_t_star_centered(1:n-1)) / real(n - 1, kind=dp)
    Var_s_empirical = sum(s_t_star_centered**2) / real(n - 1, kind=dp)
    rho_1 = autocov / Var_s_empirical

    phi_hat = rho_1 / (1.0_dp - var_eta / Var_s)
    phi_hat = max(-acf_max, min(acf_max, phi_hat))
    sigma_h2_hat = Var_s - var_eta
    sigma_v2_hat = (1.0_dp - phi_hat**2) * sigma_h2_hat
    mu_hat = mean(s_t_star)
    sigma_v_hat = sqrt(sigma_v2_hat)

    deallocate(s_t, s_t_star, s_t_star_centered)
  end subroutine fit_arsv_zero_corr

  ! Simulate daily returns with autoregressive stochastic log volatility allowing for nonzero correlation 
  ! between return and log volatility innovations.
  !
  ! Parameters:
  ! - mu: real(dp), the long-term mean of the log volatility process.
  ! - phi: real(dp), the autoregressive parameter (|phi| < 1).
  ! - sigma_v: real(dp), the standard deviation of the volatility shocks.
  ! - rho: real(dp), the correlation between return and log volatility innovations (must lie in [-1, 1]).
  ! - T: integer, the number of time periods to simulate.
  ! - burn_in: integer, the number of initial observations to discard (to reduce the effect of initial values, default: 1000).
  !
  ! Returns:
  ! - y: real(dp) array of length T, the simulated daily returns.
  subroutine arsv_returns(mu, phi, sigma_v, rho, T, y, burn_in)
    real(kind=dp), intent(in) :: mu, phi, sigma_v, rho
    integer, intent(in) :: T, burn_in
    real(kind=dp), intent(out) :: y(T)
    integer :: total_length, i
    real(kind=dp), allocatable :: h(:), z1(:), z2(:), v(:), epsilon(:)

    total_length = T + burn_in
    allocate(h(total_length), z1(total_length), z2(total_length), v(total_length), epsilon(total_length))

    call random_seed()
    do i = 1, total_length
      z1(i) = random_normal()
      z2(i) = random_normal()
    end do

    v = z1
    epsilon = rho * z1 + sqrt(1.0_dp - rho**2) * z2

    h(1) = mu
    do i = 2, total_length
      h(i) = mu + phi * (h(i-1) - mu) + sigma_v * v(i)
    end do

    do i = 1, total_length
      h(i) = exp(h(i) / 2.0_dp) * epsilon(i)
    end do

    y = h(burn_in + 1:total_length)
    deallocate(h, z1, z2, v, epsilon)
  end subroutine arsv_returns

  ! Simulate daily returns and the corresponding conditional standard deviations with autoregressive 
  ! stochastic log volatility allowing for nonzero correlation between return and log volatility innovations.
  !
  ! Parameters:
  ! - mu: real(dp), the long-term mean of the log volatility process.
  ! - phi: real(dp), the autoregressive parameter (|phi| < 1).
  ! - sigma_v: real(dp), the standard deviation of the volatility shocks.
  ! - rho: real(dp), the correlation between return and log volatility innovations (must lie in [-1, 1]).
  ! - T: integer, the number of time periods to simulate.
  ! - burn_in: integer, the number of initial observations to discard (default: 1000).
  !
  ! Returns:
  ! - result: real(dp) array of shape (T, 2), where the first column contains the conditional standard deviations 
  !           (sigma_t = exp(h_t/2)) and the second column contains the simulated daily returns.
  subroutine arsv_vol_and_returns(mu, phi, sigma_v, rho, T, result, burn_in)
    real(kind=dp), intent(in) :: mu, phi, sigma_v, rho
    integer, intent(in) :: T, burn_in
    real(kind=dp), intent(out) :: result(T, 2)
    integer :: total_length, i
    real(kind=dp), allocatable :: h(:), sigma_t(:), z1(:), z2(:), v(:), epsilon(:)

    total_length = T + burn_in
    allocate(h(total_length), sigma_t(total_length), z1(total_length), z2(total_length), v(total_length), epsilon(total_length))

    call random_seed()
    do i = 1, total_length
      z1(i) = random_normal()
      z2(i) = random_normal()
    end do

    v = z1
    epsilon = rho * z1 + sqrt(1.0_dp - rho**2) * z2

    h(1) = mu
    do i = 2, total_length
      h(i) = mu + phi * (h(i-1) - mu) + sigma_v * v(i)
    end do

    do i = 1, total_length
      sigma_t(i) = exp(h(i) / 2.0_dp)
      h(i) = sigma_t(i) * epsilon(i)
    end do

    result(:, 1) = sigma_t(burn_in + 1:total_length)
    result(:, 2) = h(burn_in + 1:total_length)
  end subroutine arsv_vol_and_returns

end module arsv_module
