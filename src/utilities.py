import numpy as np
from scipy.spatial.distance import cdist
from numba import njit


def gen_V_matrix(size: int, fill_value: float = -1.0) -> np.ndarray:
    """
    With fill_value = -1 ge_V_matrix generates a size*size matrix:
         0  0 -1 -1 ... -1 -1 -1
         0  0  0 -1 ... -1 -1 -1
        -1  0  0  0 ... -1 -1 -1
        .          .           .
        .            .         .
        .              .       .
        -1 -1 -1 -1 ...  0  0  0
        -1 -1 -1 -1 ... -1  0  0

    Args:
        size: size of array
        fill_value: float

    Returns:
        the matrix
    """
    V = np.full((size, size), fill_value)
    np.fill_diagonal(V, 0)
    np.fill_diagonal(V[:-1, 1:], 0)
    np.fill_diagonal(V[1:, :-1], 0)
    return V


@njit  # TODO: Diameter kan kanskje regnes ut samtidig som energien, siden de deler store deler av koden
def calculate_diameter(polymer: np.ndarray) -> float:
    """Finner diameteren til et polymer

    Args:
        polymer (np.ndarray): polymeren som diameter skal finnes

    Returns:
        float: diameteren tim polymeren
    """
    N = len(polymer)
    L = np.repeat(polymer, N).reshape(2 * N, N)
    return np.sqrt(
        np.max(
            (L[::2] - L[::2].transpose()) ** 2 + (L[1::2] - L[1::2].transpose()) ** 2
        )
    )
