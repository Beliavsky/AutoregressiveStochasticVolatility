program main
! Simulate returns from an ARSV model with correlation, fit the model ignoring it,
! estimate smoothed log volatility, and print statistics and correlation with true 
! log(vol).
  use kind_mod, only: dp
  use arsv_module, only: arsv_vol_and_returns, fit_arsv_zero_corr
  use arsv_smooth_mod, only: sv_kalman_smoother
  use statistics_mod, only: print_mean_var_min_max, correl
  implicit none

  ! Parameters
  real(kind=dp), parameter :: mu = -0.5_dp       ! Long-term mean of log volatility
  real(kind=dp), parameter :: phi = 0.9_dp       ! Autoregressive parameter
  real(kind=dp), parameter :: sigma_v = 0.4_dp   ! Volatility of volatility
  real(kind=dp), parameter :: rho = -0.5_dp      ! Nonzero correlation between return and log volatility innovations
  integer, parameter :: nobs = 10**6             ! Number of returns simulated
  integer, parameter :: niter = 2                ! Number of sets of returns to simulate
  integer, parameter :: burn_in = 1001           ! Burn-in period
  integer, parameter :: nparam = 4               ! Number of parameters

  ! Variables
  real(kind=dp), allocatable :: x(:,:), vol(:), log_vol(:), returns(:), returns_div_vol(:)
  real(kind=dp), allocatable :: h_smooth(:)      ! Smoothed log volatility
  real(kind=dp) :: mu_hat, phi_hat, sigma_v_hat  ! Fitted parameters
  real(kind=dp) :: ll                            ! Log likelihood from smoother
  integer :: iter, i
  character(len=10) :: param_names(nparam) = ["mu      ", "phi     ", "sigma   ", "rho     "]
  real(kind=dp) :: true_values(nparam)

  ! Simulate from a stochastic volatility model with nonzero correlation 
  ! and fit an SV model (ignoring the correlation) to the data.

  write(*,"(A,I0,A)") "#obs: ", nobs, ""
  write(*,*)

  allocate(x(nobs, 2))
  allocate(h_smooth(nobs))  ! Allocate for smoothed log volatility
  true_values = [mu, phi, sigma_v, rho]
  allocate(vol(nobs), log_vol(nobs), returns(nobs), returns_div_vol(nobs))

  do iter = 1, niter
    ! Simulate returns and volatility with nonzero correlation
    call arsv_vol_and_returns(mu, phi, sigma_v, rho, nobs, x, burn_in)

    ! Print parameter names and true values (including rho)
    print "(/,*(a12))", "", (trim(param_names(i)), i=1,nparam)
    print "(a12,*(f12.6))", "true", true_values

    ! Fit the ARSV model to the simulated returns (ignoring rho)
    call fit_arsv_zero_corr(x(:, 2), mu_hat, phi_hat, sigma_v_hat)
    print "(a12,*(f12.6))", "fit", mu_hat, phi_hat, sigma_v_hat

    ! Estimate smoothed log volatility and log likelihood
    call sv_kalman_smoother(x(:, 2), mu, phi, sigma_v, h_smooth, ll)
    print "(a12,f12.1)", "log lik", ll
    print*

    ! Compute and print descriptive statistics for the simulated series
    vol = x(:, 1)
    log_vol = log(x(:, 1))
    returns = x(:, 2)
    returns_div_vol = x(:, 2) / x(:, 1)
    call print_mean_var_min_max("vol", vol, print_labels=.true.)
    call print_mean_var_min_max("log(vol)", log_vol)
    call print_mean_var_min_max("returns", returns)
    call print_mean_var_min_max("returns/vol", returns_div_vol)
    call print_mean_var_min_max("h_smooth/2", h_smooth/2)
    print "(/,a,f10.4)", "correlation of log(vol) and h_smooth:", &
       correl(log_vol, h_smooth)
  end do

  deallocate(x, h_smooth, vol, log_vol, returns, returns_div_vol)

end program main
