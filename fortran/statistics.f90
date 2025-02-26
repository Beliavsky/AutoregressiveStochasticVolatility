module statistics_mod
  use kind_mod, only: dp
  implicit none
  private
  public :: corr_mat, mean, sd, variance, print_means_and_sds, print_mean_var_min_max, correl, covar
  real(kind=dp), parameter :: TOL = 1.0e-12_dp  ! Tolerance for checks

contains

  function mean(x) result(m)
    ! Compute mean of a vector
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: m
    m = sum(x) / size(x)
  end function mean

  function sd(x) result(s)
    ! Compute standard deviation of a vector
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: s, m
    integer :: n
    n = size(x)
    if (n < 2) then
      s = 0.0_dp
      return
    end if
    m = mean(x)
    s = sqrt(sum((x - m)**2) / (n - 1))
  end function sd

  function variance(x) result(s)
    ! Compute variance of a vector
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: s, m
    integer :: n
    n = size(x)
    if (n < 2) then
      s = 0.0_dp
      return
    end if
    m = mean(x)
    s = sum((x - m)**2) / (n - 1)
  end function variance

  function corr_mat(x) result(corr)
    ! Compute correlation matrix from data (cols = variables, rows = samples)
    real(kind=dp), intent(in) :: x(:,:)
    real(kind=dp), allocatable :: corr(:,:)
    integer :: n_vars, n_samples, i, j
    real(kind=dp) :: m_i, m_j, s_i, s_j, cov
    n_vars = size(x, 2)     ! Variables are now columns
    n_samples = size(x, 1)  ! Samples are now rows
    allocate(corr(n_vars, n_vars))
    do i = 1, n_vars
      m_i = mean(x(:,i))  ! Mean over column i
      s_i = sd(x(:,i))    ! Std dev over column i
      if (abs(s_i) < TOL) s_i = 1.0_dp  ! Avoid division by zero
      do j = 1, n_vars
        if (i == j) then
          corr(i,j) = 1.0_dp
        else
          m_j = mean(x(:,j))  ! Mean over column j
          s_j = sd(x(:,j))    ! Std dev over column j
          if (abs(s_j) < TOL) s_j = 1.0_dp
          cov = sum((x(:,i) - m_i) * (x(:,j) - m_j)) / (n_samples - 1)
          corr(i,j) = cov / (s_i * s_j)
        end if
      end do
    end do
  end function corr_mat

  subroutine print_means_and_sds(x)
    ! Print means and standard deviations of columns of x(:,:)
    real(kind=dp), intent(in) :: x(:,:)
    integer :: i
    print "(/,a)", "empirical means and standard deviations"
    do i = 1, size(x, 2)
      print "(A,I2,A,F12.8,A,F12.8)", "Variable ", i, &
        ": Mean = ", mean(x(:,i)), ", StdDev = ", sd(x(:,i))
    end do
  end subroutine print_means_and_sds

  subroutine print_mean_var_min_max(label, data, print_labels)
    ! Print mean, variance, min, and max of an array
    character(len=*), intent(in) :: label
    real(kind=dp), intent(in) :: data(:)
    logical, intent(in), optional :: print_labels
    real(kind=dp) :: mean_val, var_val
    integer :: n
    n = size(data)
    mean_val = sum(data) / real(n, kind=dp)
    var_val = sum((data - mean_val)**2) / real(n - 1, kind=dp)
    if (present(print_labels)) then
      if (print_labels) write (*,"(15x, *(a12))") "mean", "variance", "min", "max"
    end if
    write(*,"(a15, *(f12.6))") label, mean_val, var_val, minval(data), maxval(data)
  end subroutine print_mean_var_min_max

  function correl(x, y) result(corr)
    ! Compute the correlation coefficient between two vectors x and y
    !
    ! Parameters:
    !   x : real(dp) array, first input vector
    !   y : real(dp) array, second input vector
    !
    ! Returns:
    !   corr : real(dp), the Pearson correlation coefficient between x and y,
    !          or -2.0 if vectors have different lengths
    !
    ! Notes:
    !   - Returns 1.0 if either standard deviation is below tolerance to avoid division by zero

    real(kind=dp), intent(in) :: x(:), y(:)
    real(kind=dp) :: corr
    integer :: n_samples
    real(kind=dp) :: m_x, m_y, s_x, s_y, cov

    ! Check if vectors have the same length
    if (size(x) /= size(y)) then
      corr = -2.0_dp
      return
    end if

    ! Get sample size
    n_samples = size(x)

    ! Compute means
    m_x = mean(x)
    m_y = mean(y)

    ! Compute standard deviations
    s_x = sd(x)
    s_y = sd(y)

    ! Avoid division by zero
    if (abs(s_x) < TOL) s_x = 1.0_dp
    if (abs(s_y) < TOL) s_y = 1.0_dp

    ! Compute covariance
    cov = sum((x - m_x) * (y - m_y)) / real(n_samples - 1, kind=dp)

    ! Compute correlation
    corr = cov / (s_x * s_y)

  end function correl

  function covar(x, y) result(cov)
    ! Compute the covariance between two vectors x and y
    !
    ! Parameters:
    !   x : real(dp) array, first input vector
    !   y : real(dp) array, second input vector
    !
    ! Returns:
    !   cov : real(dp), the sample covariance between x and y,
    !         or -2.0 if vectors have different lengths

    real(kind=dp), intent(in) :: x(:), y(:)
    real(kind=dp) :: cov
    integer :: n_samples
    real(kind=dp) :: m_x, m_y

    ! Check if vectors have the same length
    if (size(x) /= size(y)) then
      cov = -2.0_dp
      return
    end if

    ! Get sample size
    n_samples = size(x)

    ! Compute means
    m_x = mean(x)
    m_y = mean(y)

    ! Compute covariance
    cov = sum((x - m_x) * (y - m_y)) / real(n_samples - 1, kind=dp)

  end function covar

end module statistics_mod
