import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def illustrate_polymer(
    polymer: np.ndarray,
    cmap: str = "binary",
    numbers: bool = False,
    title: str = "",
    save_to_file: bool = False,
    filename: str | Path = "",
) -> None:
    """
    Uses matplotlib.pyplot.pcolormesh to illustrate a polymer.

    Args:
        polymer: Nx2-dimensional array containing coordinates for the N monomers
        cmap: matplotlib colormap
        number: Defaults to False. If True the monomers will display their index
        title: Title
        save_to_file: Defaults to False
        filename: if save_to_file is True filename has to be specified.

    Returns:
        None
    """
    fig, ax = plt.subplots()

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
        Z[y_coord, x_coord] = i + 2
        if numbers:
            ax.text(x_coord + 0.5, y_coord + 0.5, i + 1, size="x-large", color='red')

    ax.pcolormesh(x, y, Z, shading="flat", cmap=cmap)
    ax.set(title=title, xticks=x, xticklabels=[], yticks=y, yticklabels=[])
    plt.tick_params(axis="both", left=False, right=False, bottom=False, top=False)
    ax.grid()

    if not save_to_file:
        plt.show()
    else:
        print(f"function illustrate_polymer() in {__file__}:\n\tsaving a plot to {filename}")
        plt.savefig(filename)


if __name__ == "__main__":
    a = np.array([[0, 0], [0, 1], [0, 2]])
    illustrate_polymer(a, numbers=True)
