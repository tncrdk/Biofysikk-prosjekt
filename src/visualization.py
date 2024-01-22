import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def illustrate_polymer(
    polymer: np.ndarray,
    cmap: str = "binary",
    title: str = "",
    save_to_file: bool = False,
    filename: str | Path = "",
) -> None:
    """
    Uses matplotlib.pyplot.pcolormesh to illustrate a polymer.

    Args:
        polymer:
            Nx2-dimensional array containing coordinates for the N monomers.

        cmap:
            matplotlib colormap

        title:
            Title

        save_to_file:
            Defaults to False.

        filename:
            if save_to_file is True filename has to be specified. Example: "test.png"

    Returns:
        None
    """
    # Make a NxN-grid
    min_x, min_y = np.min(polymer, axis=0)
    max_x, max_y = np.max(polymer, axis=0)
    N = max(max_y - min_y, max_x - min_x) + 1
    x = np.arange(N + 1)
    y = np.arange(N + 1)
    Z = np.zeros((N, N))

    # Generate a color gradient by assigning an integer to every monomer
    i = 1
    for monomer in polymer:
        Z[monomer[1] - min_y, monomer[0] - min_x] = i
        i += 1

    fig, ax = plt.subplots()
    ax.pcolormesh(x, y, Z, shading="flat", cmap=cmap)
    ax.set(title=title, xticks=x, xticklabels=[], yticks=y, yticklabels=[])
    plt.tick_params(axis="both", left=False, right=False, bottom=False, top=False)
    ax.grid()

    if not save_to_file:
        plt.show()
    elif filename == "":
        raise Exception(
            f"{__file__}: if save_to_file==True you have to supply a filename"
        )
    else:
        print(f"{__file__}: saving to {filename}")
        plt.savefig(filename)
