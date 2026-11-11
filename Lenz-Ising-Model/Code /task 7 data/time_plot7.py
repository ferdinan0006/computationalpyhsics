import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

threads = [1, 2, 4, 6, 8, 10]

filenames = []
total_times = []

for i in threads:
    filename = f"results_{i}.csv"
    filenames.append(filename)

for number_of_threads, name in zip(threads, filenames):
    df = pd.read_csv(name)

    total_time = df["total_simulation_second"].iloc[0]

    total_times.append(total_time)


fig, ax = plt.subplots(figsize=(9, 5.5), dpi=150)

#plot curve
ax.plot(
    threads,
    total_times,
    color="#2563EB",
    linewidth=2.5,
    marker="o",
    markersize=8,
    markerfacecolor="white",
    markeredgewidth=2.5,
    label="Measured runtime",
    zorder=3
)

#runtime labels added
for n, t in zip(threads, total_times):
    ax.annotate(
        f"{t:.2f} s",
        (n, t),
        xytext=(0, 12),
        textcoords="offset points",
        ha="center",
        fontsize=10,
        fontweight="medium",
        color="#1E293B"
    )

'''ax.set_title(
    "OpenMP Parallel Performance",
    fontsize=17,
    fontweight="bold",
    pad=20
)'''

ax.set_xlabel("Number of threads", fontsize=12, labelpad=10)
ax.set_ylabel("Total simulation time (s)", fontsize=12, labelpad=10)

#axis formating
ax.set_xticks(threads)
ax.set_ylim(0, max(total_times) * 1.2)

#subtle horizontal grid
ax.grid(axis="y", linestyle="--", alpha=0.25)
ax.grid(axis="x", visible=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

#add experiment details
ax.text(
    0.98, 0.96,
    r"$L=20,\quad N_{\mathrm{MC}}=10^5$",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=18,
    color="#475569"
)

fig.tight_layout()
plt.savefig("time_comparison.png", dpi=300, bbox_inches="tight")
plt.show()