import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path


def Motion_Animation(data: np.ndarray, filename: Path, x: int = 0, y: int = 1, xmin =-250, xmax = 250, ymin = -250, ymax = 250):
    """
    Creates a GIF animation of one or more particles.
    Input:
        data: Either (timesteps, columns) for one particle or
              (particles, timesteps, columns) for multiple particles
        filename: Path: Path of the to be saved image
        x: column of data where x values are stored
        y: column of data where y values are stored
        xmin,xmax,ymin,ymax are used to set the limits of the plot
    """
    particle_data = np.asarray(data)
    #Add a third dimension if there is only one particle
    if particle_data.ndim == 2:
        particle_data = particle_data[np.newaxis, :, :]
    elif particle_data.ndim != 3:
        raise ValueError("data must have shape (timesteps, columns) or (particles, timesteps, columns)")

    #Plot formatting
    fig, ax = plt.subplots(figsize=(5, 4))
    point = ax.scatter([], [])
    x_padding = 0.0125 *(xmin - xmax)
    y_padding = 0.0125 *(ymin - ymax)
    ax.set_xlim(xmin - x_padding, xmax + x_padding)
    ax.set_ylim(ymin - y_padding, ymax + y_padding)

    #Function to update the frame
    def update(frame):
        offsets = np.column_stack((particle_data[:, frame, x], particle_data[:, frame, y]))
        point.set_offsets(offsets)
        return (point,)

    n_timesteps = particle_data.shape[1]
    frame_indices = np.linspace(0, n_timesteps - 1, min(n_timesteps, 500), dtype=int)
    ani = FuncAnimation(fig, update, frames=frame_indices, blit=True)
    ani.save(filename, fps=15, dpi=300)
    plt.close(fig)

def Line_plot(
    data: np.ndarray,
    filename: Path,
    x: int = 0,
    y: int = 1,
    hline=[],
    xlabel="",
    ylabel="",
    ret=False,
    equal_aspect: bool = False,
):
    """
    Creates a Line plot to represent the motion of a particle in one dimension
    Input:
            data: np.array: Should contain the x and y values of the particles displayed in the animation
            filename: Path: Path of the to be saved image
            x: column of data where x values are stored
            y: column of data where y values are stored
            xmin,xmax,ymin,ymax are used to set the limits of the plot
    """
    particle_data = np.asarray(data)

    #Reformate data to make One particle and multiple particle simulations plottable
    if particle_data.ndim == 2:
        particle_data = particle_data[np.newaxis, :, :]
    elif particle_data.ndim != 3:
        raise ValueError(
            "data must have shape (timesteps, columns) or "
            "(series, timesteps, columns)"
        )

    if x >= particle_data.shape[2] or y >= particle_data.shape[2]:
        raise ValueError("x and y must be valid column indices")

    if ret:
        ax = plt.gca()
    else:
        fig, ax = plt.subplots()

    for i, series in enumerate(particle_data):
        ax.plot(series[:, x], series[:, y], label = "Particle " + str(i))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()

    for value in hline:
        ax.axhline(y=value, color="black")

    ax.margins(x=0.05, y=0.05)
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")

    if ret:
        return ax

    fig.savefig(filename)
    plt.close(fig)


def Relative_Error_Plot(data: np.ndarray, filename: Path, x: int = 0, y: int = 1,
                        xlabel: str = "", ylabel: str = "", labels=None, ret: bool = False, ymin = None, ymax = None):
    """Plot one or more relative-error series with a logarithmic y axis."""
    error_data = np.asarray(data)

    if error_data.ndim == 2:
        error_data = error_data[np.newaxis, :, :]
    elif error_data.ndim != 3:
        raise ValueError(
            "data must have shape (timesteps, columns) or "
            "(series, timesteps, columns)"
        )

    if x >= error_data.shape[2] or y >= error_data.shape[2]:
        raise ValueError("x and y must be valid column indices")
    if labels is not None and len(labels) != error_data.shape[0]:
        raise ValueError("labels must contain one entry per series")

    if ret:
        ax = plt.gca()
    else:
        fig, ax = plt.subplots()

    for index, series in enumerate(error_data):
        label = labels[index] if labels is not None else f"Series {index}"
        ax.plot(series[:, x], series[:, y], label=label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if(ymin != None or ymax != None): ax.set_ylim(ymin,ymax)
    ax.set_yscale("log")
    ax.legend()

    if ret:
        return ax

    fig.savefig(filename)
    plt.close(fig)
    

