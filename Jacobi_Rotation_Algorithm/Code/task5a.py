import numpy as np
import matplotlib.pyplot as plt
# Set global params for plots
plt.rcParams.update(
    {
        "font.family": "serif",  
        "mathtext.fontset": "cm",
        "font.size": 11,
        "axes.labelsize": 15,
        "axes.titlesize": 18,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.titlesize": 14,
    }
)

data = np.loadtxt('./Task5_iterations.txt', delimiter = ' ', skiprows = 1)
xdata = data[:, 0]
xdata2 = xdata**2
ydata = data[:, 1]
sqrty = np.sqrt(ydata)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(xdata, ydata, label='Iterations', color='blue', marker='o')
ax1.set_title('Iterations for different Matrix sizes')
ax1.set_xlabel('N')
ax1.set_ylabel('Iterations')
ax1.grid(True)
ax1.legend()
ax1.set_xlim([0, xdata.max()*1.1])
ax1.set_ylim([0, ydata.max()*1.1])

ax2.plot(xdata, sqrty, label=r'$\sqrt{\text{Iterations}}$', color='blue', marker='o')
ax2.set_title(r'$\sqrt{\text{Iterations}}$ for different Matrix sizes')
ax2.set_xlabel('Matrix size N')
ax2.set_ylabel(r'$\sqrt{\text{Iterations}}$')
ax2.grid(True)
ax2.legend()
ax2.set_xlim([0, xdata.max()*1.1])
ax2.set_ylim([0, sqrty.max()*1.1])

plt.tight_layout()
plt.savefig("./task5a.pdf")
plt.show()
