from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from plotting_funcs import Motion_Animation
from plotting_funcs import Line_plot
from plotting_funcs import Relative_Error_Plot
from format_data import data_multi

REPORT_COLUMN_FIGSIZE = (3.35, 2.5)
REPORT_FULL_WIDTH_FIGSIZE = (6.8, 4.8)

plt.rcParams.update(
    {
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "lines.linewidth": 1,
        "figure.figsize": REPORT_COLUMN_FIGSIZE,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
    }
)

# Animation generation is disabled by default because it is time intensive.
# Uncomment the Motion_Animation calls below to generate animations.

# Set up paths and simulation parameters.
project_dir = Path(__file__).resolve().parent.parent
data_dir = project_dir / "data"
figures_dir = project_dir / "figures"
animation_dir = figures_dir / "Animations"
animation_dir.mkdir(parents=True, exist_ok=True)

mass = 39.962590863
charge_number = 1.0
V0 = 2.41e6
d0 = 500.0


def add_total_energy(data):
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    velocity = data[:, 3:6]
    kinetic_energy = 0.5 * mass * np.sum(velocity**2, axis=1)
    electric_potential = V0 / (2 * d0**2) * (2 * z**2 - x**2 - y**2)
    potential_energy = charge_number * electric_potential
    return np.column_stack((data, kinetic_energy + potential_energy))

# Load data for a single particle.
dataRK4 = np.loadtxt(data_dir / "RK4_single_32000_0.txt")
dataFE = np.loadtxt(data_dir / "FE_single_32000_0.txt")
number_steps = len(dataRK4)

dataRK4 = add_total_energy(dataRK4)
dataFE = add_total_energy(dataFE)

# Uncomment these calls to generate single-particle animations.
# Motion_Animation(dataRK4, animation_dir / f"RungeKutta_{number_steps}.gif")
# Motion_Animation(dataFE, animation_dir / f"ForwardEuler_{number_steps}.gif")

# Plot the movement in the z-t plane.
Line_plot(
    dataRK4,
    x=6,
    y=2,
    filename=figures_dir / f"RungeKutta_{number_steps}.pdf",
    hline=[20, -20],
    xlabel=r"time in $\mu$s",
    ylabel=r"z coordinate in $\mu$m",
    ret=False,
)
Line_plot(
    dataFE,
    x=6,
    y=2,
    filename=figures_dir / f"ForwardEuler_{number_steps}.pdf",
    hline=[20, -20],
    xlabel=r"time in $\mu$s",
    ylabel=r"z coordinate in $\mu$m",
    ret=False,
)

# Compare RK4 and FE in z-t space for a single particle.
fig, ax = plt.subplots(
    2, 1, sharex="all", figsize=REPORT_FULL_WIDTH_FIGSIZE
)
plt.sca(ax[0])
ax[0] = Line_plot(
    dataFE,
    x=6,
    y=2,
    filename=figures_dir / f"ForwardEuler_{number_steps}.pdf",
    hline=[20, -20],
    xlabel=r"time in $\mu$s",
    ylabel=r"z coordinate in $\mu$m",
    ret=True,
)
plt.sca(ax[1])
ax[1] = Line_plot(
    dataRK4,
    x=6,
    y=2,
    filename=figures_dir / f"RungeKutta_{number_steps}.pdf",
    hline=[20, -20],
    xlabel=r"time in $\mu$s",
    ylabel=r"z coordinate in $\mu$m",
    ret=True,
)
fig.savefig(figures_dir / "RK4_vs_FE_z-t-plot.pdf")

# Compare total energy for RK4 and FE for a single particle.
fig, ax = plt.subplots(
    2, 1, sharex=True, sharey=True, figsize=REPORT_FULL_WIDTH_FIGSIZE
)
energy_unit = r"$\mathrm{u}\cdot\frac{(\mu\mathrm{m})^2}{(\mu\mathrm{s})^2}$"
ax[0].plot(dataFE[:, 6], dataFE[:, 7])
ax[0].set_ylabel("FE\nEnergy " + energy_unit)
ax[1].plot(dataRK4[:, 6], dataRK4[:, 7])
ax[1].set_ylabel("RK4\nEnergy " + energy_unit)
ax[1].set_xlabel(r"time in $\mu$s")
ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
fig.tight_layout()
fig.savefig(figures_dir / "RK4_vs_FE_total_energy.pdf")

# Compare the analytical and numerical solutions for RK4 and FE.
sizes = [4000, 8000, 16000, 32000]
datasets = [
    np.loadtxt(
        data_dir / f"analytical_single_{size}_0.txt",
        usecols=(3, 4, 5),
        ndmin=2,
    )
    for size in sizes]

max_rows = max(data.shape[0] for data in datasets)
rel_error = np.full((len(datasets), max_rows, 3), np.nan)

for index, data in enumerate(datasets):
    rel_error[index, :data.shape[0], :] = data

Relative_Error_Plot(
    rel_error,
    x=2, y=0,
    filename=figures_dir / "relative_error_RK4.pdf",
    xlabel=r"time in $\mu$s",
    ylabel="relative error",
    labels=[f"N={size}" for size in sizes],
    # ymin=10**(-12)
)

Relative_Error_Plot(
    rel_error,
    x=2, y=1,
    filename=figures_dir / "relative_error_FE.pdf",
    xlabel=r"time in $\mu$s",
    ylabel="relative error",
    labels=[f"N={size}" for size in sizes],
    # ymin=10**(-3)
)

# Create a table of final errors to confirm the global error order.
final_errors_RK4 = np.array([data[-1, 0] for data in datasets])
final_errors_FE = np.array([data[-1, 1] for data in datasets])
timesteps = np.array([data[1, 2] - data[0, 2] for data in datasets])
final_errors_RK4_dt = final_errors_RK4 / timesteps**4
final_errors_FE_dt = final_errors_FE / (timesteps)
Table = np.column_stack((timesteps, final_errors_FE, final_errors_RK4,
                         final_errors_FE_dt, final_errors_RK4_dt))
np.savetxt(data_dir / "Error_Table.txt", Table)

# Compare FE and RK in terms of relative error.
fig, ax = plt.subplots(
    2, 1, sharex="all", figsize=REPORT_FULL_WIDTH_FIGSIZE
)
plt.sca(ax[0])
ax[0] = Relative_Error_Plot(
    rel_error,
    x=2, y=1,
    filename=figures_dir / "relative_error_FE.pdf",
    xlabel=r"time in $\mu$s",
    ylabel="Relative error)",
    labels=[f"N={size}" for size in sizes],
    ret=True,
    ymin=10**(-3),
)
plt.sca(ax[1])
ax[1] = Relative_Error_Plot(
    rel_error,
    x=2, y=0,
    filename=figures_dir / "relative_error_RK.pdf",
    xlabel=r"time in $\mu$s",
    ylabel="Relative error)",
    labels=[f"N={size}" for size in sizes],
    ret=True,
    ymin=10**(-12),
)
fig.savefig(figures_dir / "RK4_vs_FE_rel_error-plot.pdf")

# Format data for two particles.
data_FE_2_Particles_x_y = data_multi(2, data_dir / "FE_2_w_inter_", (0,1))
data_RK_2_Particles_x_y = data_multi(2, data_dir / "RK4_2_w_inter_", (0,1))
data_FE_2_Particles_x_y_no_inter = data_multi(2, data_dir / "FE_2_no_inter_", (0,1))
data_RK_2_Particles_x_y_no_inter = data_multi(2, data_dir / "RK4_2_no_inter_", (0,1))

# Uncomment these calls to generate two-particle animations.
# Motion_Animation(data_RK_2_Particles_x_y, animation_dir / f"RungeKutta_{number_steps}_2Particles.gif")
# Motion_Animation(data_FE_2_Particles_x_y, animation_dir / f"ForwardEuler_{number_steps}_2Particles.gif")
# Motion_Animation(data_RK_2_Particles_x_y_no_inter, animation_dir / f"RungeKutta_{number_steps}_2Particles.gif")
# Motion_Animation(data_FE_2_Particles_x_y_no_inter, animation_dir / f"ForwardEuler_{number_steps}_2Particles.gif")

# Load RK4 trajectories for the two-particle phase-space plots.
rk4_x_vx_no_inter = data_multi(2, data_dir / "RK4_2_no_inter_", (0, 3))
rk4_x_vx_with_inter = data_multi(2, data_dir / "RK4_2_w_inter_", (0, 3))
rk4_z_vz_no_inter = data_multi(2, data_dir / "RK4_2_no_inter_", (2, 5))
rk4_z_vz_with_inter = data_multi(2, data_dir / "RK4_2_w_inter_", (2, 5))

phase_space_path = figures_dir / "RK4_two_particle_phase_space.pdf"
fig, axes = plt.subplots(2, 2, figsize=(6.8, 5.2))
phase_space_plots = [
    (
        axes[0, 0],
        rk4_x_vx_no_inter,
        "No interactions: $x$-$v_x$",
        r"x in $\mu$m",
        r"$v_x$ in $\mu$m / $\mu$s",
    ),
    (
        axes[0, 1],
        rk4_x_vx_with_inter,
        "With interactions: $x$-$v_x$",
        r"x in $\mu$m",
        r"$v_x$ in $\mu$m / $\mu$s",
    ),
    (
        axes[1, 0],
        rk4_z_vz_no_inter,
        "No interactions: $z$-$v_z$",
        r"z in $\mu$m",
        r"$v_z$ in $\mu$m / $\mu$s",
    ),
    (
        axes[1, 1],
        rk4_z_vz_with_inter,
        "With interactions: $z$-$v_z$",
        r"z in $\mu$m",
        r"$v_z$ in $\mu$m / $\mu$s",
    ),
]

for axis, data, title, xlabel, ylabel in phase_space_plots:
    plt.sca(axis)
    Line_plot(
        data,
        filename=phase_space_path,
        x=0,
        y=1,
        xlabel=xlabel,
        ylabel=ylabel,
        ret=True,
    )
    axis.set_title(title)

fig.tight_layout()
fig.savefig(phase_space_path)

# Create data in the (x, y) phase space.
data_FE_2_Particles_x_y = data_multi(2, data_dir / "FE_2_w_inter_", (0,1))
data_RK_2_Particles_x_y = data_multi(2, data_dir / "RK4_2_w_inter_", (0,1))
data_FE_2_Particles_x_y_no_inter = data_multi(2, data_dir / "FE_2_no_inter_", (0,1))
data_RK_2_Particles_x_y_no_inter = data_multi(2, data_dir / "RK4_2_no_inter_", (0,1))

# Plot data in the (x, y) phase space.
Line_plot(
    data_FE_2_Particles_x_y,
    x=0,
    y=1,
    filename=figures_dir / f"FE_x_y.pdf",
    hline=[],
    xlabel=r"x in $\mu$m",
    ylabel=r"y in $\mu$m",
    ret=False,
    equal_aspect=True,
)
Line_plot(
    data_RK_2_Particles_x_y,
    x=0,
    y=1,
    filename=figures_dir / f"RK_x_y.pdf",
    hline=[],
    xlabel=r"x in $\mu$m",
    ylabel=r"y in $\mu$m",
    ret=False,
    equal_aspect=True,
)
Line_plot(
    data_FE_2_Particles_x_y_no_inter,
    x=0,
    y=1,
    filename=figures_dir / f"FE_x_y_no_inter.pdf",
    hline=[],
    xlabel=r"x in $\mu$m",
    ylabel=r"y in $\mu$m",
    ret=False,
    equal_aspect=True,
)
Line_plot(
    data_RK_2_Particles_x_y_no_inter,
    x=0,
    y=1,
    filename=figures_dir / f"RK_x_y_no_inter.pdf",
    hline=[],
    xlabel=r"x in $\mu$m",
    ylabel=r"y in $\mu$m",
    ret=False,
    equal_aspect=True,
)

# Compare RK4 and FE in the (x, y) phase space for ttwo particles with and withour interactions.
fig, ax = plt.subplots(
    1, 2, sharex="all", figsize=(REPORT_FULL_WIDTH_FIGSIZE[0], 3.6)
)
plt.sca(ax[0])
ax[0] = Line_plot(
    data_RK_2_Particles_x_y_no_inter,
    x=0,
    y=1,
    filename=figures_dir / f"ForwardEuler_{number_steps}.pdf",
    hline=[],
    xlabel=r"x in $\mu$m",
    ylabel=r"y in $\mu$m",
    ret=True,
    equal_aspect=True,
)
plt.sca(ax[1])
ax[1] = Line_plot(
    data_RK_2_Particles_x_y,
    x=0,
    y=1,
    filename=figures_dir / f"RungeKutta_{number_steps}.pdf",
    hline=[],
    xlabel=r"x in $\mu$m",
    ylabel=r"y in $\mu$m",
    ret=True,
    equal_aspect=True,
)
fig.savefig(figures_dir / "RK4_vs_FE_x-y-plot.pdf")


# Task 9 - coarse resonance sweep
fs = [0.1, 0.4, 0.7]
w0 = 0.2
dw = 0.02
n_particles = 100

counts = np.loadtxt(data_dir / "coarse_resonance_sweep.txt")
w = w0 + np.arange(counts.shape[1]) * dw
fractions = counts / n_particles

plt.figure(figsize=(9, 5.5))
for f_val, frac in zip(fs, fractions):
    plt.plot(w, frac, lw=1.8, label=rf"$f = {f_val}$")

plt.xlabel(r" Perturbation frequency $\omega_V$ [MHz]", fontsize=12)
plt.ylabel(r"Fraction of trapped particles", fontsize=12)
plt.title("Fraction of Particles Remaining in Trap after $500\\,\\mu\\mathrm{s}$", fontsize=13)
plt.xlim(w[0], w[-1])
plt.ylim(-0.05, 1.05)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=11)
plt.tight_layout()

plt.savefig(figures_dir / "coarse_resonance_sweep.pdf")

# Task 9 - fine resonance sweep (no interactions)
w0 = 1.0
dw = 0.002

counts = np.loadtxt(data_dir / "fine_resonance_sweep_nointer.txt")
w = w0 + np.arange(counts.shape[1]) * dw
fractions = counts / n_particles

plt.figure(figsize=(9, 5.5))
for f_val, frac in zip(fs, fractions):
    plt.plot(w, frac, lw=1.8, label=rf"$f = {f_val}$")

plt.xlabel(r" Perturbation frequency $\omega_V$ [MHz]", fontsize=12)
plt.ylabel(r"Fraction of trapped particles", fontsize=12)
plt.title("Fraction of Particles Remaining in Trap after $500\\,\\mu\\mathrm{s}$", fontsize=13)
plt.xlim(w[0], w[-1])
plt.ylim(-0.05, 1.05)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=11)
plt.tight_layout()

plt.savefig(figures_dir / "fine_resonance_sweep_nointer.pdf")

# Task 9 - fine resonance sweep (with interactions)
w0 = 1.0
dw = 0.002

counts = np.loadtxt(data_dir / "fine_resonance_sweep_inter.txt")
w = w0 + np.arange(counts.shape[1]) * dw
fractions = counts / n_particles

plt.figure(figsize=(9, 5.5))
for f_val, frac in zip(fs, fractions):
    plt.plot(w, frac, lw=1.8, label=rf"$f = {f_val}$")

plt.xlabel(r" Perturbation frequency $\omega_V$ [MHz]", fontsize=12)
plt.ylabel(r"Fraction of trapped particles", fontsize=12)
plt.title("Fraction of Particles Remaining in Trap after $500\\,\\mu\\mathrm{s}$", fontsize=13)
plt.xlim(w[0], w[-1])
plt.ylim(-0.05, 1.05)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=11)
plt.tight_layout()

plt.savefig(figures_dir / "fine_resonance_sweep_inter.pdf")