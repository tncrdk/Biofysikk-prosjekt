import numpy as np

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