# Assignment 2: simulation design
#
# Time to breast cancer recurrence, with five-year follow-up and independent dropout.
#
# 1. How accurately can a parametric survival model estimate the probability of
#    breast cancer recurrence over time based on routinely recorded patient and
#    clinical characteristics?
# 2. How does increasing dropout affect estimation accuracy?
#
# Compare the full Weibull PH model with a model that omits tumour size.
#
# Run from this folder with `uv sync --locked` and
# `uv run --locked python 1_data_generation.py`.

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import WeibullAFTFitter
from lifelines.exceptions import ConvergenceError
from scipy.special import expit


SEED = 2026
N = 800
REPETITIONS = 2000
TAU = 5.0
TIMES = (1, 3, 5)
SHAPE = 0.7
RAW_RATE = 0.001

# Centre and scale covariates without changing the model in the report.
COVARIATES = ["age_c10", "size_c10", "chemo"]
BETA = np.array([0.5, 1.0, -0.7])
RATE = RAW_RATE * np.exp(0.05 * 65 + 0.1 * 20)


def generate_covariates(n, rng):
    age = rng.gamma(shape=20, scale=3.25, size=n)
    while np.any(age < 18):
        underage = age < 18
        age[underage] = rng.gamma(shape=20, scale=3.25, size=underage.sum())
    size = rng.lognormal(mean=2.8, sigma=0.5, size=n)
    chemo = rng.binomial(1, expit(0.25 * (size - 15)))
    return pd.DataFrame({
        "age": age,
        "tumor_size": size,
        "chemo": chemo,
        "age_c10": (age - 65) / 10,
        "size_c10": (size - 20) / 10,
    })


def generate_dataset(n=N, dropout=0.2, rng=None):
    if not 0 <= dropout < 1:
        raise ValueError("dropout must be between 0 (inclusive) and 1 (exclusive)")
    rng = np.random.default_rng() if rng is None else rng
    data = generate_covariates(n, rng)
    lp = data[COVARIATES].to_numpy() @ BETA
    event_time = (rng.exponential(size=n) / (RATE * np.exp(lp))) ** (1 / SHAPE)
    rho = -np.log1p(-dropout) / TAU
    censor_time = (np.full(n, np.inf) if dropout == 0
                   else rng.exponential(scale=1 / rho, size=n))
    data["time"] = np.minimum(event_time, np.minimum(censor_time, TAU))
    data["event"] = ((event_time <= censor_time) & (event_time <= TAU)).astype(int)
    data["dropout_observed"] = (censor_time < event_time) & (censor_time < TAU)
    return data


def fit_weibull_ph(data, covariates=None):
    covariates = list(COVARIATES if covariates is None else covariates)
    if data["event"].sum() == 0:
        return None
    model = WeibullAFTFitter(penalizer=0.0)
    try:
        model.fit(data[covariates + ["time", "event"]],
                  duration_col="time", event_col="event")
    except ConvergenceError:
        return None

    # Convert the library's AFT parameters to the PH parameters in the report.
    shape = np.exp(model.params_.loc["rho_", "Intercept"])
    aft = model.params_.loc["lambda_"]
    return {"rate": np.exp(-shape * aft["Intercept"]), "shape": shape,
            "beta": -shape * aft.loc[covariates].to_numpy(),
            "covariates": covariates}


# ---------------------------------------------------------------- Example data
# Here dropout means P(C<5), not the overall censored fraction.

data = generate_dataset(rng=np.random.default_rng(SEED))
outdir = Path("data")
outdir.mkdir(exist_ok=True)
observed = ["age", "tumor_size", "chemo", "age_c10", "size_c10", "time", "event"]
data[observed].to_csv(outdir / "recurrence_sim_n800.csv", index=False)

params = {
    "gamma": SHAPE, "lambda": RAW_RATE,
    "beta_age": 0.05, "beta_size": 0.1, "beta_chemo": -0.7,
    "age_gamma_shape": 20, "age_gamma_scale": 3.25, "age_min": 18,
    "size_log_mean": 2.8, "size_log_sd": 0.5,
    "chemo_logit_intercept": -3.75, "chemo_logit_size_slope": 0.25,
    "dropout_by_tau": 0.2, "dropout_rate": -np.log(0.8) / TAU,
    "tau": TAU, "time_unit": "years", "n": N, "seed": SEED,
    "covariates": COVARIATES, "age_center": 65, "size_center": 20,
    "covariate_unit": 10, "lambda_centered": RATE, "beta_scaled": BETA.tolist(),
    "evaluation_times": list(TIMES),
    "prediction_profile": {"age": 65, "tumor_size": 20, "chemo": 0},
}
(outdir / "dgp_params.json").write_text(json.dumps(params, indent=2) + "\n")
print(f"n={len(data)}, events={data['event'].sum()}, censored={1 - data['event'].mean():.1%}")


# ---------------------------------------------------------- Repeated simulations
# Each run generates a new cohort and fits both models to the same observed data.
# The comparison uses the chemotherapy log-hazard ratio and recurrence risks at
# 1, 3 and 5 years for a woman aged 65, with a 20 mm tumour and no chemotherapy.

MODELS = {"full": COVARIATES, "no_size": ["age_c10", "chemo"]}
TRUTH = {"beta_chemo": -0.7}
TRUTH.update({f"risk_{t}y": -np.expm1(-RATE * t ** SHAPE) for t in TIMES})


def run_simulation(repetitions=REPETITIONS, n=N, seed=SEED):
    rng = np.random.default_rng(seed)
    rows = []
    for dropout in [0.0, 0.2, 0.5]:
        for run in range(repetitions):
            data = generate_dataset(n, dropout, rng)
            for model, covariates in MODELS.items():
                fit = fit_weibull_ph(data, covariates)
                rows.append({
                    "dropout_target": dropout, "model": model, "run": run,
                    "censoring": 1 - data["event"].mean(),
                    "dropout_observed": data["dropout_observed"].mean(),
                    "converged": fit is not None,
                    "beta_chemo": fit["beta"][-1] if fit else np.nan,
                    **{f"risk_{t}y": -np.expm1(-fit["rate"] * t ** fit["shape"])
                       if fit else np.nan for t in TIMES},
                })
    return pd.DataFrame(rows)


def summarize(results):
    rows = []
    for (dropout, model), group in results.groupby(["dropout_target", "model"], sort=False):
        for parameter, truth in TRUTH.items():
            values = group[parameter].dropna()
            rows.append({
                "dropout_target": dropout, "model": model, "parameter": parameter,
                "valid_runs": len(values), "failed_fraction": 1 - group["converged"].mean(),
                "censoring": group["censoring"].mean(),
                "dropout_observed": group["dropout_observed"].mean(),
                "bias": values.mean() - truth,
                "variance": values.var(ddof=0),
                "mse": ((values - truth) ** 2).mean(),
            })
    return pd.DataFrame(rows)


results = run_simulation()
summary = summarize(results)
results.to_csv(outdir / "simulation_runs.csv", index=False)
summary.to_csv(outdir / "simulation_summary.csv", index=False)
print(summary.round(5).to_string(index=False))
