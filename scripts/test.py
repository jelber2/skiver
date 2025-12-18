import numpy as np
import math

# Your data
p_obs = np.array([
    0.9425859654582232, 0.9367778144602179, 0.9423788546255507,
    0.943717277486911,  0.9463047354864276, 0.9501675041876047,
    0.9539444689290436, 0.9544929544929545, 0.9571636011616651,
    0.9560050568900127, 0.9677334038614124
], dtype=float)

v = np.arange(21, 32, dtype=float)  # 21..31

def p_model(alpha: float, beta: float, v: np.ndarray) -> np.ndarray:
    return (beta + v) / (alpha + beta + v)

def objective(alpha: float, beta: float, lam: float = 2000.0, target_zero: float = 0.9) -> float:
    """
    Fit p_model to p_obs, with a soft constraint that beta/(alpha+beta) ~ target_zero.
    lam controls how hard you enforce that constraint.
    """
    if alpha <= 0.0 or beta <= 0.0:
        return 1e100

    pm = p_model(alpha, beta, v)
    sse = float(np.mean((pm - p_obs) ** 2))

    zero_rate = beta / (alpha + beta)  # = 1 - alpha/(alpha+beta)
    penalty = lam * (zero_rate - target_zero) ** 2

    return sse + penalty

def fit_with_scipy():
    from scipy.optimize import minimize

    # Optimize in log-space to enforce positivity
    def obj_log(x):
        a = math.exp(x[0])
        b = math.exp(x[1])
        return objective(a, b)

    # A reasonable starting point: beta ~ 9*alpha to make zero_rate ~ 0.9
    x0 = np.log([2.0, 18.0])

    res = minimize(obj_log, x0, method="Nelder-Mead", options={"maxiter": 20000})
    a_hat, b_hat = float(math.exp(res.x[0])), float(math.exp(res.x[1]))
    return a_hat, b_hat, res

def fit_grid_fallback():
    # Coarse grid search (works without SciPy)
    best = (float("inf"), None, None)
    for a in np.logspace(-2, 2, 800):        # alpha in [1e-2, 1e2]
        for f in np.logspace(-1, 1, 200):    # beta around 9*alpha times factor in [0.1,10]
            b = 9.0 * a * f
            val = objective(a, b)
            if val < best[0]:
                best = (val, a, b)
    return best[1], best[2], best[0]

if __name__ == "__main__":
    alpha, beta, res = fit_with_scipy()
    #alpha, beta, val = fit_grid_fallback()
    

    zero_rate = beta / (alpha + beta)
    pm = p_model(alpha, beta, v)

    print(f"alpha = {alpha:.6g}")
    print(f"beta  = {beta:.6g}")
    print(f"1 - alpha/(alpha+beta) = beta/(alpha+beta) = {zero_rate:.6g}")
    print(f"mean squared error on p_v = {float(np.mean((pm - p_obs)**2)):.6g}")
    print("p_model(v=21..31) =", pm)