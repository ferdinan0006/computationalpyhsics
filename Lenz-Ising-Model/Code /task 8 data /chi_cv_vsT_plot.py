import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#load data 
gridsizes = [40, 60, 80, 100]

chis = []
cvs = []
temperatures = []

for L in gridsizes:
    filename = f"results_zoomed_in{L}.csv"
    df = pd.read_csv(filename)

    chi = df["chi"].to_numpy()
    cv = df["cv"].to_numpy()
    temp = df["T"].to_numpy()

    chis.append(chi)
    cvs.append(cv)
    temperatures.append(temp)

#create to subplots 
fig, axes = plt.subplots(
    2, 1,
    figsize=(10, 11),
    dpi=150,
    sharex=True
)

ax1, ax2 = axes

colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(gridsizes)))

# Exact critical temperature
Tc_exact = 2 / np.log(1 + np.sqrt(2))

#chi plot
for L, temp, chi,cv, color in zip(gridsizes, temperatures, chis, cvs, colors):
    Tc_chi = temp[np.argmax(chi)]
    Tc_cv = temp[np.argmax(cv)]

    print(
        f"L = {L:3d}: "
        f"Tc(chi) = {Tc_chi:.3f}, "
        f"Tc(Cv) = {Tc_cv:.3f}"
    )
    ax1.plot(
        temp,
        chi,
        marker="o",
        markersize=4,
        linewidth=2,
        color=color,
        label=f"L = {L}"
    )

ax1.axvline(
    Tc_exact,
    color="black",
    linestyle="--",
    linewidth=1.5,
    alpha=0.7,
    label=r"Exact $T_c$"
)

ax1.set_title(
    "Magnetic Susceptibility",
    fontsize=15,
    fontweight="bold",
    pad=15
)

ax1.set_ylabel(r"Susceptibility $\chi/N$", fontsize=18)
ax1.legend(frameon=False, fontsize=14)

#cv subplot
for L, temp, cv, color in zip(gridsizes, temperatures, cvs, colors):
    ax2.plot(
        temp,
        cv,
        marker="o",
        markersize=4,
        linewidth=2,
        color=color,
        label=f"L = {L}"
    )

ax2.axvline(
    Tc_exact,
    color="black",
    linestyle="--",
    linewidth=1.5,
    alpha=0.7,
    label=r"Exact $T_c$"
)

ax2.set_title(
    "Heat Capacity",
    fontsize=15,
    fontweight="bold",
    pad=15
)

ax2.set_xlabel(r"Temperature $T$", fontsize=18)
ax2.set_ylabel(r"Heat capacity $C_V/N$", fontsize=18)

#design
for ax in axes:
    ax.grid(alpha=0.2, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=14)

'''#title
fig.suptitle(
    "2D Ising Model: Critical Behavior for Different Lattice Sizes",
    fontsize=17,
    fontweight="bold"
)'''

fig.tight_layout(rect=[0, 0, 1, 0.97])

plt.savefig(
    "critical_behavior_comparison.png",
    dpi=300,
    bbox_inches="tight")    

plt.show()