import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd

filenames = ['burn_in_T1_ordered_per_spin.csv', 
             'burn_in_T1_random_per_spin.csv', 
             'burn_in_T24_ordered_per_spin.csv', 
             'burn_in_T24_random_per_spin.csv']

# empty dic that can be filled 
data = {}

for file in filenames: 
    df = pd.read_csv(file)
    number_cycles = df.iloc[:, 0].to_numpy()
    cumulative_e = df.iloc[:, 1].to_numpy()
    cumulative_abs = df.iloc[:, 2].to_numpy()

    #calculating the current energy by subtracting (i)-i(-1), and for i = 0, just take the value 
    current_e = df.iloc[:, 3].to_numpy()

    #store everything in dic, using filename as key 
    data[file] = {
        'cycles': number_cycles,
        'cumulative_e': cumulative_e,
        'current_e': current_e
    }

#extract data 
n_MC_cycles = data[filenames[0]]['cycles']

# current energy
eps_T1_o = data[filenames[0]]['current_e']
eps_T1_u = data[filenames[1]]['current_e']  # u = unordered
eps_T24_o = data[filenames[2]]['current_e']
eps_T24_u = data[filenames[3]]['current_e']

#cumulative energy
mean_eps_T1_o = data[filenames[0]]['cumulative_e']
mean_eps_T1_u = data[filenames[1]]['cumulative_e']
mean_eps_T24_o = data[filenames[2]]['cumulative_e']
mean_eps_T24_u = data[filenames[3]]['cumulative_e']

#plots using structure from the pdf 
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

#current energy
ax1.plot(n_MC_cycles, eps_T1_u, '-', color='#377eb8', alpha=0.4, linewidth=2.0, label='Current $E$ (T=1.0, Random)')
ax1.plot(n_MC_cycles, eps_T1_o, '-', color='#4daf4a', alpha=0.4, linewidth=2.0, label='Current $E$ (T=1.0, Ordered)')
ax1.plot(n_MC_cycles, eps_T24_u, '-', color='#e41a1c', alpha=0.4, linewidth=2.0, label='Current $E$ (T=2.4, Random)')
ax1.plot(n_MC_cycles, eps_T24_o, '-', color='#984ea3', alpha=0.4, linewidth=2.0, label='Current $E$ (T=2.4, Ordered)')

#cumulative
ax2.plot(n_MC_cycles, mean_eps_T1_u, '-', linewidth=2.0, color='#377eb8', label='Mean $E$ (T=1.0, Random)')
ax2.plot(n_MC_cycles, mean_eps_T1_o, '-', linewidth=2.0, color='#4daf4a', label='Mean $E$ (T=1.0, Ordered)')
ax2.plot(n_MC_cycles, mean_eps_T24_u, '-', linewidth=2.0, color='#e41a1c', label='Mean $E$ (T=2.4, Random)')
ax2.plot(n_MC_cycles, mean_eps_T24_o, '-', linewidth=2.0, color='#984ea3', label='Mean $E$ (T=2.4, Ordered)')

#formatting for readablity 
ax1.set_ylabel('Current Energy per Spin $E/N$', fontsize = '18')
ax2.set_ylabel('Energy per Spin $\\langle E \\rangle / N$', fontsize = '18')
ax2.set_xlabel('Number of Monte Carlo Cycles', fontsize = '18')
#no title just caption
#fig.suptitle('Burn-in Time Analysis: Current vs. Cumulative Energy')
ax1.tick_params(axis='both', labelsize=14)
ax2.tick_params(axis='both', labelsize=14)

#logarithmic x axis 
ax2.set_xscale('log')

#legend outside the plot, nothing will be covered
ax1.legend(loc='upper right', bbox_to_anchor=(1.02, 1), fontsize='small', ncol=1)
ax2.legend(loc='upper right', bbox_to_anchor=(1.02, 1), fontsize='small', ncol=1)

ax1.grid(True, alpha=0.3)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('burnin_time_all_t.png')
plt.show()