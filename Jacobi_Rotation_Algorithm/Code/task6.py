import numpy as np
import matplotlib.pyplot as plt

def parse_eigen_data(filename):
    eigenvalues = []
    eigenvectors = []
    
    parsing_vals = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            if not line:
                continue
                
            if line == "Eigenvalues:":
                parsing_vals = True
                continue
            elif line == "Eigenvectors:":
                parsing_vals = False
                continue
                
            if parsing_vals:
                eigenvalues.append(float(line))
            else:
                row = [float(val) for val in line.split()]
                eigenvectors.append(row)
    
    evals = np.array(eigenvalues)
    evecs = np.array(eigenvectors)
    
    return evals, evecs

def plot_comparison(evals_num, evecs_num, evals_an, evecs_an, title_suffix):
    N = evecs_num.shape[0]
    x = np.linspace(0, 1, N)
    
    plt.figure(figsize=(10, 6))
    
    for i in range(3):
        v_num = evecs_num[:, i]
        v_an = evecs_an[:, i]
        
        # Align signs
        if np.dot(v_num, v_an) < 0:
            v_num = -v_num
            
        plt.plot(x, v_an, linestyle='-', linewidth=2, 
                 label=f'Analytical $v_{i+1}$ ($\\lambda={evals_an[i]:.4f}$)')
        
        marker_style = ('o', 's', '^')[i] if title_suffix == "10" else None
        plt.plot(x, v_num, linestyle=('--', '-.', ':')[i], marker=marker_style, markersize=5, color="black",
                 label=f'Numerical $v_{i+1}$ ($\\lambda={evals_num[i]:.4f}$)')

    plt.title(f'Solution Comparison (n={title_suffix})')
    plt.xlabel('x')
    plt.ylabel('Amplitude')
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"task6_{title_suffix}.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    evals_10, evecs_10 = parse_eigen_data("task6_10.txt")
    evals_10_an, evecs_10_an = parse_eigen_data("task6_10_an.txt")

    evals_100, evecs_100 = parse_eigen_data("task6_100.txt")
    evals_100_an, evecs_100_an = parse_eigen_data("task6_100_an.txt")

    plot_comparison(evals_10, evecs_10, evals_10_an, evecs_10_an, "10")
    plot_comparison(evals_100, evecs_100, evals_100_an, evecs_100_an, "100")