
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

gridsizes = [40, 60, 80, 100]
chis = []
cvs = []
temperatures = []
Tc_chis = []
Tc_cvs = []

for L in gridsizes:
    filename = f"../task 8 data /results_zoomed_in{L}.csv"
    df = pd.read_csv(filename)

    chi = df["chi"].to_numpy()
    cv = df["cv"].to_numpy()
    temp = df["T"].to_numpy()

    chis.append(chi)
    cvs.append(cv)
    temperatures.append(temp)

for L, temp, chi, cv in zip(gridsizes, temperatures, chis, cvs):
    Tc_chi = temp[np.argmax(chi)]
    Tc_cv = temp[np.argmax(cv)]

    Tc_chis.append(Tc_chi)
    Tc_cvs.append(Tc_cv)

#building the inverse
inverse_L = 1 / np.array(gridsizes)

#scipy expects a 2d array 
X = inverse_L.reshape(-1, 1)

Tc_chis = np.array(Tc_chis)
Tc_cvs = np.array(Tc_cvs)

#do linear regression with scikit 
model_chi = LinearRegression()
model_chi.fit(X, Tc_chis)

model_cv = LinearRegression()
model_cv.fit(X, Tc_cvs)

a_chi = model_chi.coef_[0]
Tc_inf_chi = model_chi.intercept_

a_cv = model_cv.coef_[0]
Tc_inf_cv = model_cv.intercept_

Tc_exact = 2 / np.log(1 + np.sqrt(2))

error_chi = abs(Tc_inf_chi - Tc_exact)
error_cv = abs(Tc_inf_cv - Tc_exact)

relative_error_chi = error_chi / Tc_exact * 100
relative_error_cv = error_cv / Tc_exact * 100

r2_chi = model_chi.score(X, Tc_chis)
r2_cv = model_cv.score(X, Tc_cvs)

print("Susceptibility:")
print(f"  a              = {a_chi:.6f}")
print(f"  Tc(infinity)   = {Tc_inf_chi:.6f}")
print(f"  Absolute error = {error_chi:.6f}")
print(f"  Relative error = {relative_error_chi:.3f}%")
print(f"  R²             = {r2_chi:.4f}")

print("\nHeat capacity:")
print(f"  a              = {a_cv:.6f}")
print(f"  Tc(infinity)   = {Tc_inf_cv:.6f}")
print(f"  Absolute error = {error_cv:.6f}")
print(f"  Relative error = {relative_error_cv:.3f}%")
print(f"  R²             = {r2_cv:.4f}")

#fitted curves
x_fit = np.linspace(0, inverse_L.max() * 1.1, 200)
X_fit = x_fit.reshape(-1, 1)

y_fit_chi = model_chi.predict(X_fit)
y_fit_cv = model_cv.predict(X_fit)

#plot
fig, ax = plt.subplots(figsize=(10, 6.5), dpi=150)

#using same color set as before
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(gridsizes)))

chi_color = colors[0]
cv_color = colors[-1]

# Standard deviation of linear regression residuals
residuals_chi = Tc_chis - model_chi.predict(X)
residuals_cv = Tc_cvs - model_cv.predict(X)

sd_chi = np.sqrt(np.sum(residuals_chi**2) / (len(Tc_chis) - 2))
sd_cv = np.sqrt(np.sum(residuals_cv**2) / (len(Tc_cvs) - 2))

print(f"Regression SD (chi) = {sd_chi:.6f}")
print(f"Regression SD (cv)  = {sd_cv:.6f}")

# Constant-width error bands
ax.fill_between(
    x_fit,
    y_fit_chi - sd_chi,
    y_fit_chi + sd_chi,
    color=chi_color,
    alpha=0.14,
    label=r"$\chi/N$ $\pm 1\sigma$ regression band"
)

ax.fill_between(
    x_fit,
    y_fit_cv - sd_cv,
    y_fit_cv + sd_cv,
    color=cv_color,
    alpha=0.14,
    label=r"$C_V/N$ $\pm 1\sigma$ regression band"
)

ax.scatter(
    inverse_L,
    Tc_chis,
    color=chi_color,
    s=85,
    marker="o",
    zorder=4,
    label=r"$T_c(L)$ from $\chi/N$"
)

ax.plot(
    x_fit,
    y_fit_chi,
    color=chi_color,
    linewidth=2.5,
    label=fr"$\chi/N$: $R^2={r2_chi:.3f}$"
)

ax.scatter(
    inverse_L,
    Tc_cvs,
    color=cv_color,
    s=85,
    marker="s",
    zorder=4,
    label=r"$T_c(L)$ from $C_V/N$"
)

ax.plot(
    x_fit,
    y_fit_cv,
    color=cv_color,
    linewidth=2.5,
    label=fr"$C_V/N$: $R^2={r2_cv:.3f}$"
)

ax.scatter(
    [0, 0],
    [Tc_inf_chi, Tc_inf_cv],
    color=[chi_color, cv_color],
    marker="X",
    s=130,
    zorder=5
)

ax.axhline(
    Tc_exact,
    color="black",
    linestyle="--",
    linewidth=1.5,
    alpha=0.7,
    label=fr"Onsager: $T_c={Tc_exact:.4f}$"
)

ax.set_title(
    "Finite-Size Scaling of the 2D Ising Model",
    fontsize=16,
    fontweight="bold",
    pad=18
)

ax.set_xlabel(r"Inverse lattice size $1/L$", fontsize=12)
ax.set_ylabel(r"Critical temperature $T_c(L)$", fontsize=12)

ax.set_xlim(left=0, right=inverse_L.max() * 1.12)

ax.grid(alpha=0.2, linestyle="--")
ax.legend(frameon=False, fontsize=9)

ax.spines["left"].set_position(("data", 0))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.tick_params(axis="both", labelsize=10)

fig.tight_layout()

plt.savefig(
    "finite_size_scaling.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
