from __future__ import annotations

import math

import numpy as np
import pandas as pd


BASE_COVARIATES = [
    "INDFMPIR",
    "BMXBMI",
    "MEAN_SBP",
    "LOG_COT",
    "LOG_UCR",
    "LBDSCR",
]


def _invert_full_rank(matrix: np.ndarray) -> np.ndarray:
    """Invert a small full-rank square matrix with partial pivoting.

    The survey models use a low-dimensional design matrix. A compact
    Gauss-Jordan implementation avoids platform-specific linear-algebra
    backend stalls while keeping the numerical operation explicit.
    """
    array = np.asarray(matrix, dtype=float)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("matrix must be square")
    size = array.shape[0]
    augmented = np.concatenate([array.copy(), np.eye(size)], axis=1)
    tolerance = max(float(np.max(np.abs(array))), 1.0) * 1e-12

    for column in range(size):
        pivot = column + int(np.argmax(np.abs(augmented[column:, column])))
        if abs(augmented[pivot, column]) <= tolerance:
            raise np.linalg.LinAlgError("design matrix is not full rank")
        if pivot != column:
            augmented[[column, pivot]] = augmented[[pivot, column]]
        augmented[column] /= augmented[column, column]
        for row in range(size):
            if row == column:
                continue
            augmented[row] -= augmented[row, column] * augmented[column]
    return augmented[:, size:]


def _beta_continued_fraction(a: float, b: float, x: float) -> float:
    maximum_iterations, tolerance, minimum = 200, 3e-14, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = minimum if abs(d) < minimum else d
    d = 1.0 / d
    result = d
    for iteration in range(1, maximum_iterations + 1):
        doubled = 2 * iteration
        value = iteration * (b - iteration) * x / ((qam + doubled) * (a + doubled))
        d = 1.0 + value * d
        d = minimum if abs(d) < minimum else d
        c = 1.0 + value / c
        c = minimum if abs(c) < minimum else c
        d = 1.0 / d
        result *= d * c
        value = -(a + iteration) * (qab + iteration) * x / (
            (a + doubled) * (qap + doubled)
        )
        d = 1.0 + value * d
        d = minimum if abs(d) < minimum else d
        c = 1.0 + value / c
        c = minimum if abs(c) < minimum else c
        d = 1.0 / d
        change = d * c
        result *= change
        if abs(change - 1.0) < tolerance:
            break
    return result


def _regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    scale = math.exp(
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return scale * _beta_continued_fraction(a, b, x) / a
    return 1.0 - scale * _beta_continued_fraction(b, a, 1.0 - x) / b


def _two_sided_t_p(t_value: float, degrees_of_freedom: int) -> float:
    x = degrees_of_freedom / (degrees_of_freedom + t_value**2)
    return float(_regularized_incomplete_beta(degrees_of_freedom / 2.0, 0.5, x))


def _t_critical(degrees_of_freedom: int, alpha: float = 0.05) -> float:
    low, high = 0.0, 12.0
    for _ in range(80):
        midpoint = (low + high) / 2.0
        if _two_sided_t_p(midpoint, degrees_of_freedom) > alpha:
            low = midpoint
        else:
            high = midpoint
    return (low + high) / 2.0


def _design(frame: pd.DataFrame, term: str, covariates: list[str]) -> pd.DataFrame:
    age10 = (frame["RIDAGEYR"].to_numpy(float) - 40.0) / 10.0
    design = pd.DataFrame(
        {
            "intercept": 1.0,
            term: frame[term].to_numpy(float),
            "age10": age10,
            "age10_sq": age10**2,
            "female": frame["RIAGENDR"].eq(2).astype(float).to_numpy(),
        }
    )
    for name in covariates:
        design[name.lower()] = frame[name].to_numpy(float)
    race = pd.get_dummies(
        frame["RIDRETH1"].astype(int), prefix="race", drop_first=True, dtype=float
    ).reset_index(drop=True)
    return pd.concat([design, race], axis=1)


def survey_linear_regression(
    frame: pd.DataFrame,
    outcome: str,
    term: str = "MIX_PRIMARY_Z",
    covariates: list[str] | None = None,
    weight: str = "WTSPP2YR",
) -> dict:
    covariates = BASE_COVARIATES if covariates is None else covariates
    required = [
        outcome,
        term,
        "RIDAGEYR",
        "RIAGENDR",
        "RIDRETH1",
        weight,
        "SDMVSTRA",
        "SDMVPSU",
        *covariates,
    ]
    data = (
        frame[required]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .reset_index(drop=True)
    )
    design = _design(data, term, covariates)
    x = design.to_numpy(float)
    y = data[outcome].to_numpy(float)
    w = (data[weight] / data[weight].mean()).to_numpy(float)

    gram = np.einsum("ni,n,nj->ij", x, w, x, optimize=False)
    bread = _invert_full_rank(gram)
    right_hand_side = np.einsum("ni,n,n->i", x, w, y, optimize=False)
    beta = np.einsum("ij,j->i", bread, right_hand_side, optimize=False)
    residual = y - np.einsum("ij,j->i", x, beta, optimize=False)
    meat = np.zeros((x.shape[1], x.shape[1]))
    n_psu = 0
    n_strata = 0
    strata = data["SDMVSTRA"].astype(str).to_numpy()
    psu = data["SDMVPSU"].astype(str).to_numpy()
    for stratum in np.unique(strata):
        in_stratum = strata == stratum
        stratum_psus = np.unique(psu[in_stratum])
        if len(stratum_psus) < 2:
            continue
        n_strata += 1
        n_psu += len(stratum_psus)
        totals = []
        for cluster in stratum_psus:
            mask = in_stratum & (psu == cluster)
            totals.append(
                np.einsum(
                    "ni,n,n->i",
                    x[mask],
                    w[mask],
                    residual[mask],
                    optimize=False,
                )
            )
        totals = np.asarray(totals)
        centered = totals - totals.mean(axis=0, keepdims=True)
        meat += (
            len(stratum_psus)
            / (len(stratum_psus) - 1.0)
            * np.einsum("ni,nj->ij", centered, centered, optimize=False)
        )

    covariance = np.einsum(
        "ia,ab,jb->ij", bread, meat, bread, optimize=False
    )
    design_df = int(n_psu - n_strata)
    index = list(design.columns).index(term)
    estimate = float(beta[index])
    standard_error = float(math.sqrt(max(covariance[index, index], 0.0)))
    t_value = estimate / standard_error
    critical = _t_critical(design_df)
    return {
        "n": int(len(data)),
        "design_df": design_df,
        "beta": estimate,
        "se": standard_error,
        "t": t_value,
        "p": _two_sided_t_p(t_value, design_df),
        "ci_low": estimate - critical * standard_error,
        "ci_high": estimate + critical * standard_error,
    }


def fixed_effect_summary(rows: list[dict]) -> dict:
    estimates = np.asarray([row["beta"] for row in rows], dtype=float)
    errors = np.asarray([row["se"] for row in rows], dtype=float)
    weights = 1.0 / errors**2
    estimate = float(np.sum(weights * estimates) / np.sum(weights))
    standard_error = float(math.sqrt(1.0 / np.sum(weights)))
    z_value = estimate / standard_error
    return {
        "beta": estimate,
        "se": standard_error,
        "z": z_value,
        "p": float(math.erfc(abs(z_value) / math.sqrt(2.0))),
        "ci_low": estimate - 1.959963984540054 * standard_error,
        "ci_high": estimate + 1.959963984540054 * standard_error,
    }
