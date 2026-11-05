import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

#only use oredered start, after reaching equilibrium the initial conditions do not matter anymore 
filename = ['task 5 data /burn_in_T1_ordered_per_spin.csv', 'task 5 data /burn_in_T24_ordered_per_spin.csv']
equilibrium = 10000

data = {}

for file in filename: 
    df = pd.read_csv(file)
    df.columns = df.columns.str.strip()

    equilibrium_epsilon = df[df['cycle'] > equilibrium]['current_energy']

    mean_e = equilibrium_epsilon.mean()
    var_e = equilibrium_epsilon.var()  # Variance sigma^2_epsilon

    data[file] = {
        'equilibrium_epsilon': equilibrium_epsilon,
        'mean_e': mean_e,
        'var_e': var_e
    }

fig, axes = plt.subplots(2, 1, figsize=(6, 8))

for ax, file in zip(axes, filename):
    equilibrium_epsilon = data[file]['equilibrium_epsilon']
    mean_e = data[file]['mean_e']
    var_e = data[file]['var_e']

    #histogram
    ax.hist(
        equilibrium_epsilon,
        bins=40,
        density=True,
        color='#5B8DB8',
        alpha=0.8,
        edgecolor='white',
        linewidth=0.7
    )

    #mean energy
    ax.axvline(
        mean_e,
        color='#D65A5A',
        linestyle='--',
        linewidth=2,
        label=fr'Mean $\langle \epsilon \rangle = {mean_e:.4f}$'
    )

    #formatting for readability
    ax.set_xlabel(r'Instantaneous energy per spin $\epsilon = E/N$', fontsize=13)
    ax.set_ylabel(r'Probability density $p_\epsilon(\epsilon; T)$', fontsize=13)

    T = 1.0 if file == filename[0] else 2.4
    ax.set_title(f'Equilibrium Energy Distribution (T={T})', fontsize=16, pad=15)

    #subtle grid
    ax.grid(axis='y', alpha=0.2, linestyle='--')
    ax.set_axisbelow(True)

    #remove unnecessary borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    #legend with variance
    ax.legend(
        title=fr'Variance: $\sigma^2 = {var_e:.6f}$',
        fontsize=11,
        title_fontsize=10,
        frameon=False
    )

plt.tight_layout()
plt.savefig('normalizedHistogram.png')
plt.show()