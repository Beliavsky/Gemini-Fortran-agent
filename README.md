# Gemini-Fortran-agent
Python script that uses Gemini LLMs to create Fortran programs, iterating with compiler error messages until they compile. Sample output:

```
c:\python\code\gemini>python xgemini.py
prompt:
Do a Fortran simulation to find the optimal trimmed mean estimator of
the location of the Cauchy distribution, trying trimming proportions
of 0%, 10%, 20%, 30%, 40%, and 45%. Declare real variables as
real(kind=dp) with dp a module constant, and put procedures in a
module. Have the simulation use 100 samples of 1000 observations each.
Only output Fortran code. Do not give commentary.

model: gemini-2.5-pro-exp-03-25

Attempt 1 failed with error (generation time: 38.606 seconds, LOC=151):
cauchy.f90:128:30:

  128 |   real(dp) :: current_alpha, current_tmean
      |                              1~~~~~~~~~~~~
Error: Unused variable 'current_tmean' declared at (1) [-Werror=unused-variable]
f951.exe: some warnings being treated as errors

Code compiled successfully after 2 attempts (generation time: 18.096 seconds, LOC=158)!
Final version:

 module simulation_utilities
  implicit none
  private
  public :: dp, generate_cauchy, trimmed_mean

  integer, parameter :: dp = selected_real_kind(p=15, r=300)
  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  subroutine generate_cauchy(n, x)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n)
    real(dp) :: u(n)
    integer :: i

    ! Generate uniform random numbers (0, 1)
    call random_number(u)

    ! Apply inverse CDF method for standard Cauchy
    do i = 1, n
      ! Avoid u exactly 0 or 1, though probability is negligible
      if (u(i) == 0.0_dp) u(i) = epsilon(u(i))
      if (u(i) == 1.0_dp) u(i) = 1.0_dp - epsilon(u(i))
      x(i) = tan(pi * (u(i) - 0.5_dp))
    end do
  end subroutine generate_cauchy

  ! Simple selection sort (in-place)
  subroutine selection_sort(arr)
    real(dp), intent(inout) :: arr(:)
    integer :: n, i, j, min_idx
    real(dp) :: temp

    n = size(arr)
    if (n <= 1) return

    do i = 1, n - 1
      min_idx = i
      do j = i + 1, n
        if (arr(j) < arr(min_idx)) then
          min_idx = j
        end if
      end do
      ! Swap arr(i) and arr(min_idx)
      if (min_idx /= i) then
        temp = arr(i)
        arr(i) = arr(min_idx)
        arr(min_idx) = temp
      end if
    end do
  end subroutine selection_sort

  function trimmed_mean(x, alpha) result(t_mean)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: alpha
    real(dp) :: t_mean
    integer :: n, k, i
    real(dp), allocatable :: sorted_x(:)
    real(dp) :: sum_trimmed

    n = size(x)
    if (n == 0) then
        t_mean = 0.0_dp ! Or NaN, or handle error
        return
    end if

    ! Calculate number of elements to trim from each end
    k = floor(real(n, dp) * alpha)

    ! Ensure k is within valid bounds
    k = max(0, k)
    if (2 * k >= n) then
        ! Handle case where trimming leaves no elements or negative number
        ! For alpha < 0.5 this should not happen if n > 0
        if (n > 0) then
            ! If alpha >= 0.5, could return median or handle as error
            print *, "Warning: Trimming proportion alpha=", alpha, " is too high (>= 0.5 or leaves <= 0 elements). Returning median."
            allocate(sorted_x(n))
            sorted_x = x
            call selection_sort(sorted_x)
            if (mod(n, 2) == 1) then
                t_mean = sorted_x(n/2 + 1) ! Odd number of elements
            else
                t_mean = (sorted_x(n/2) + sorted_x(n/2 + 1)) / 2.0_dp ! Even number of elements
            end if
            deallocate(sorted_x)
        else
            t_mean = 0.0_dp ! Should not happen if n > 0 check passed
        end if
        return
    end if

    ! Allocate and copy data to sort
    allocate(sorted_x(n))
    sorted_x = x

    ! Sort the data
    call selection_sort(sorted_x)

    ! Sum the middle elements
    sum_trimmed = 0.0_dp
    do i = k + 1, n - k
      sum_trimmed = sum_trimmed + sorted_x(i)
    end do

    ! Calculate the mean of the trimmed sample
    if ( (n - 2 * k) > 0 ) then
      t_mean = sum_trimmed / real(n - 2 * k, dp)
    else
       ! This case should be caught by the check 'if (2 * k >= n)' above,
       ! but adding a safeguard here.
       print *, "Error: Zero elements remaining after trimming. Check alpha and n."
       t_mean = 0.0_dp ! Or NaN
    end if


    deallocate(sorted_x)

  end function trimmed_mean

end module simulation_utilities

program cauchy_simulation
  use simulation_utilities
  implicit none

  integer, parameter :: nsamp = 100    ! Number of samples (replications)
  integer, parameter :: nobs = 1000   ! Number of observations per sample
  real(dp), parameter :: alphas(*) = [0.0_dp, 0.10_dp, 0.20_dp, 0.30_dp, 0.40_dp, 0.45_dp]
  integer, parameter :: num_alphas = size(alphas)

  real(dp), allocatable :: sample(:)      ! Holds one sample of Cauchy deviates
  real(dp) :: estimates(nsamp)          ! Holds trimmed means for one alpha across samples
  real(dp) :: mse(num_alphas)           ! Holds Mean Squared Error for each alpha
  real(dp) :: current_alpha              ! Variable 'current_tmean' removed as it was unused
  real(dp) :: min_mse, optimal_alpha
  integer :: i, j, optimal_idx

  ! Allocate memory for one sample
  allocate(sample(nobs))

  ! Initialize random seed
  call random_seed()

  print *, "Starting Cauchy location estimation simulation..."
  print *, "Number of samples (replications): ", nsamp
  print *, "Observations per sample: ", nobs
  print *, "--------------------------------------------------"

  ! Loop over trimming proportions
  do i = 1, num_alphas
    current_alpha = alphas(i)
    print '(A, F4.2, A)', "Processing alpha = ", current_alpha * 100.0_dp, "% ..."

    ! Loop over samples (replications) for the current alpha
    do j = 1, nsamp
      ! Generate a new sample from the standard Cauchy distribution
      call generate_cauchy(nobs, sample)

      ! Calculate the trimmed mean for this sample
      estimates(j) = trimmed_mean(sample, current_alpha)
    end do

    ! Calculate the Mean Squared Error (MSE) for the current alpha
    ! Since the true location of the standard Cauchy is 0, MSE = mean(estimate^2)
    mse(i) = sum(estimates**2) / real(nsamp, dp)

    print '(A, F6.2, A, E12.5)', "  Alpha = ", current_alpha*100.0_dp, "% | MSE = ", mse(i)

  end do ! End loop over alphas

  print *, "--------------------------------------------------"
  print *, "Simulation complete. Finding optimal alpha..."

  ! Find the alpha with the minimum MSE
  min_mse = mse(1)
  optimal_alpha = alphas(1)
  optimal_idx = 1
  do i = 2, num_alphas
    if (mse(i) < min_mse) then
      min_mse = mse(i)
      optimal_alpha = alphas(i)
      optimal_idx = i
    end if
  end do

  print '(A, F6.2, A)', "Optimal trimming proportion: ", optimal_alpha * 100.0_dp, "%"
  print '(A, E12.5)', "Minimum Mean Squared Error: ", min_mse
  print *, "--------------------------------------------------"

  ! Deallocate memory
  deallocate(sample)

end program cauchy_simulation

Running executable: .\cauchy.exe

Output:
  Starting Cauchy location estimation simulation...
 Number of samples (replications):          100
 Observations per sample:         1000
 --------------------------------------------------
Processing alpha = 0.00% ...
  Alpha =   0.00% | MSE =  0.46403E+04
Processing alpha = ****% ...
  Alpha =  10.00% | MSE =  0.43341E-02
Processing alpha = ****% ...
  Alpha =  20.00% | MSE =  0.36191E-02
Processing alpha = ****% ...
  Alpha =  30.00% | MSE =  0.21073E-02
Processing alpha = ****% ...
  Alpha =  40.00% | MSE =  0.24019E-02
Processing alpha = ****% ...
  Alpha =  45.00% | MSE =  0.30734E-02
 --------------------------------------------------
 Simulation complete. Finding optimal alpha...
Optimal trimming proportion:  30.00%
Minimum Mean Squared Error:  0.21073E-02
 --------------------------------------------------


Total generation time: 56.702 seconds across 2 attempts

Compilation command: gfortran -O0 -fmax-errors=1 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -o cauchy cauchy.f90
```
