import numpy as np

n = 200
rate = 0.05
repeats = 2000
rng = np.random.default_rng(42)

for target in [0.0, 0.1, 0.5]:
    naive_estimates = []
    mle_estimates = []
    censoring_proportions = []

    # Independent exponential censoring: P(C < T) = rho / (rate + rho).
    rho = rate * target / (1 - target)

    for _ in range(repeats):
        # NumPy uses scale = 1 / rate, not the rate itself.
        event_times = rng.exponential(scale=1 / rate, size=n)
        if target == 0:
            censoring_times = np.full(n, np.inf)  # No censoring.
        else:
            censoring_times = rng.exponential(scale=1 / rho, size=n)

        observed_times = np.minimum(event_times, censoring_times)
        event_observed = event_times <= censoring_times
        events = event_observed.sum()
        total_time = observed_times.sum()

        naive_estimates.append(n / total_time)
        # The assignment asks us to record NA when no events are observed.
        mle_estimates.append(events / total_time if events > 0 else np.nan)
        censoring_proportions.append(1 - event_observed.mean())

    print(f"\nTarget censoring: {target:.0%}")
    print(f"Average observed censoring: {np.mean(censoring_proportions):.4%}")
    print(f"Replicates with no events: {np.mean(np.isnan(mle_estimates)):.2%}")
    print(f"{'Estimator':<10} {'Valid':>6} {'Bias':>14} {'Variance':>14} {'MSE':>14}")

    for name, estimates in [("Naive", naive_estimates), ("MLE", mle_estimates)]:
        estimates = np.asarray(estimates)
        estimates = estimates[~np.isnan(estimates)]

        if estimates.size > 0:
            bias = estimates.mean() - rate
            # ddof=0 makes empirical MSE = variance + bias**2 exactly.
            variance = estimates.var(ddof=0)
            mse = np.mean((estimates - rate) ** 2)
        else:
            bias = variance = mse = np.nan

        print(f"{name:<10} {estimates.size:>6} {bias:>14.6e} "
              f"{variance:>14.6e} {mse:>14.6e}")
