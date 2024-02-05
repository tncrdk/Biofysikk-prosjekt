import numpy as np
from scipy.spatial.distance import cdist
from numba import njit


def gen_V_matrix(
    size: int, fill_value: float | tuple[float, float] = -1.0
) -> np.ndarray:
    """
    With fill_value = -1 gen_V_matrix generates a size*size matrix:
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
        fill_value: float | tuple[lower, upper]:
            if tuple fill_values are drawn from a uniform distribution on [lower, upper].

    Returns:
        the matrix
    """
    if type(fill_value) == float:
        V = np.full((size, size), fill_value)
        np.fill_diagonal(V, 0)
        np.fill_diagonal(V[:-1, 1:], 0)
        np.fill_diagonal(V[1:, :-1], 0)
        return V

    else:
        assert len(fill_value) == 2, "There is an error in the type of fill_value"
        V = np.zeros((size, size))
        for i in range(1, size):
            for j in range(i):
                value = np.random.default_rng().uniform(fill_value[0], fill_value[1])
                V[i, j] = value
                V[j, i] = value
        np.fill_diagonal(V, 0)
        np.fill_diagonal(V[:-1, 1:], 0)
        np.fill_diagonal(V[1:, :-1], 0)
        # this can be done in the for-loops, but this is easier to read and is more foolproof

        return V


@njit  # TODO: Diameter kan kanskje regnes ut samtidig som energien, siden de deler store deler av koden
def calculate_diameter(polymer: np.ndarray) -> float:
    """Finds the diameter of a polymer

    Args:
        polymer (np.ndarray): the polymer to find the diameter of

    Returns:
        float: diameter of the polymer
    """
    N = len(polymer)
    L = np.repeat(polymer, N).reshape(2 * N, N)
    return np.sqrt(
        np.max(
            (L[::2] - L[::2].transpose()) ** 2 + (L[1::2] - L[1::2].transpose()) ** 2
        )
    )
