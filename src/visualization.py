import numpy as np
import matplotlib.pyplot as plt

def illustrate_polymer(
    ax,
    polymer: np.ndarray,
    cmap: str = "binary",
    numbers: bool = False,
    title: str = ""
) -> None:
    """
    Uses matplotlib.pyplot.pcolormesh to illustrate a polymer.

    Args:
        ax: the Axes to which illustrate_polymer will plot the result.
        polymer: Nx2-dimensional array containing coordinates for the N monomers
        cmap: matplotlib colormap
        number: Defaults to False. If True the monomers will display their index
        title: Title
        
    Returns:
        None
    """
    # Make a NxN-grid
    N = len(polymer)
    x = y = np.arange(N + 1)
    Z = np.zeros((N, N))

    # Placing the polymer on the grid such that the middle monomer is in the center.
    middle_monomer = polymer[int(N / 2)]
    for i, monomer in enumerate(polymer):
        shifty = (1 if N % 2 == 0 and middle_monomer[1] < 0 else 0)
        shiftx = (1 if N % 2 == 0 and middle_monomer[0] < 0 else 0)
        x_coord = monomer[1] + int(N / 2) - middle_monomer[1] - shifty
        y_coord = monomer[0] + int(N / 2) - middle_monomer[0] - shiftx
        Z[x_coord, y_coord] = i + 2
        if numbers:
            ax.text(y_coord + 0.5, x_coord + 0.5, i + 1, size="x-large", color='red')

    ax.pcolormesh(x, y, Z, shading="flat", cmap=cmap)
    ax.set(title=title, xticks=x, xticklabels=[], yticks=y, yticklabels=[])
    ax.tick_params(axis="both", left=False, right=False, bottom=False, top=False)
    ax.grid(True)

if __name__ == "__main__":
    a = np.array([[0, 0], [0, 1], [1, 1]])
    illustrate_polymer(a, numbers=True, title=str(a))
