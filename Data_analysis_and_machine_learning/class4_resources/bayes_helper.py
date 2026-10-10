import numpy as np
from scipy import stats
from scipy.special import logsumexp

def _integrate_exponential(t, y, C_sd=None, resolution=160):
    nodes, weights = np.polynomial.legendre.leggauss(resolution)
    taus = np.exp((nodes + 1) * np.log(1000) / 2)
    noise_sd = np.exp((nodes + 1) * np.log(90) / 2)
    noise_var = noise_sd**2
    prior_var = np.array([60.**2] if C_sd is None else [60.**2, C_sd**2])
    logs, means, variances = [], [], []
    for tau in taus:
        h = np.exp(-t / tau)
        X = h[:, None] if C_sd is None else np.column_stack([h, np.ones(len(t))])
        W = X * np.sqrt(prior_var)
        B = np.eye(len(prior_var))[None] + (W.T @ W)[None] / noise_var[:, None, None]
        rhs = W.T @ y
        solution = np.linalg.solve(B, np.broadcast_to(rhs, (resolution, len(rhs)))[..., None])[..., 0]
        logdet = len(y) * np.log(noise_var) + np.linalg.slogdet(B)[1]
        quadratic = (y @ y) / noise_var - np.sum(rhs * solution, axis=1) / noise_var**2
        logs.append(-0.5 * (len(y) * np.log(2*np.pi) + logdet + quadratic))
        if C_sd is not None:
            means.append(C_sd * solution[:, 1] / noise_var)
            variances.append(C_sd**2 * np.linalg.inv(B)[:, 1, 1])
    weighted_logs = np.array(logs) + np.log(weights[:, None]/2) + np.log(weights[None, :]/2)
    logZ = logsumexp(weighted_logs)
    if C_sd is None:
        return float(logZ), None
    mixture = {"weights": np.exp(weighted_logs - logZ).ravel(),
               "mean": np.array(means).ravel(), "sd": np.sqrt(variances).ravel()}
    return float(logZ), mixture

def bayesian_comparison(t, y, C_sd=20., resolution=160):
    logZ0, _ = _integrate_exponential(t, y, resolution=resolution)
    logZ1, mixture = _integrate_exponential(t, y, C_sd=C_sd, resolution=resolution)
    C_grid = np.linspace(-20, 20, 801)
    density = np.zeros_like(C_grid)
    for start in range(0, len(mixture["weights"]), 512):
        sl = slice(start, start + 512)
        pdf = stats.norm.pdf(C_grid[None, :], loc=mixture["mean"][sl, None],
                             scale=mixture["sd"][sl, None])
        density += mixture["weights"][sl] @ pdf
    C_mean = mixture["weights"] @ mixture["mean"]
    return {"BF10": np.exp(logZ1-logZ0), "logZ0": logZ0, "logZ1": logZ1,
            "C_grid": C_grid, "C_density": density, "C_mean": C_mean,
            "C_prior_sd": C_sd}
