import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd

df = pd.read_csv('task 4 data/results_task4b.csv')

variables = ['e', 'abs_m', 'cv', 'chi']
x = np.arange(len(variables))
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(df)))

fig, ax = plt.subplots(figsize=(10, 6))

for color, (_, row) in zip(colors, df.iterrows()):
    values = row[variables].to_numpy(dtype=float)
    ax.bar(
        x,
        values,
        width=0.8,
        color=color,
        alpha=0.55,
        edgecolor='black',
        label=f'T = {row["T"]:g}'
    )

ax.set_xticks(x, variables)
ax.set_ylabel('Value')
ax.set_title('Lenz-Ising model observables for different temperatures')
ax.axhline(0, color='black', linewidth=0.8)
ax.legend(title='Temperature')
ax.grid(axis='y', alpha=0.3)
fig.tight_layout()
#'fig.savefig('plot4b.png', dpi=300)'
plt.show()