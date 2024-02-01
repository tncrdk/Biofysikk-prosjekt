import numpy as np
import time
from numba import njit
import utilities

"""
polymer: [ 
          x1, y1
          x2, y2
          x3, y3
          .
          .
          .
          xN, yN
          ]
"""


def check_if_intact(polymer: np.ndarray, polymer_length: int) -> bool:
    """Checks if polymer is intact

    Args:
        polymer (np.ndarray): The polymer to check
        polymer_length (int): Length of the polymer

    Returns:
        bool: True if the polymer is intact
    """
    # Checks that the polymer has N monomers, where each has a unique whole number representation
    if np.size(np.unique(polymer, axis=0), axis=0) != polymer_length:
        return False

    for i in range(1, polymer_length):
        # Don't have to take the square root (faster), as any value different from 1 indicates a broken polymer anyway.
        distance = (polymer[i - 1, 0] - polymer[i, 0]) ** 2 + (
            polymer[i - 1, 1] - polymer[i, 1]
        ) ** 2
        if distance != 1:
            return False
    return True


def check_if_intact_2(polymer: np.ndarray, polymer_length: int) -> bool:
    """Sjekker om en polymer er intakt

    Args:
        polymer (np.ndarray): polymeren som sjekkes
        polymer_length (int): lengden til polymeren

    Returns:
        bool: True hvis polymeren er intakt
    """

    if np.size(np.unique(polymer, axis=0), axis=0) != polymer_length:
        return False  # Sjekker både at den har N monomerer og at de alle har en unik heltallsrepresentasjon

    test = polymer[1:]
    test_mot = polymer[:-1]
    distance_array = (test[:, 0] - test_mot[:, 0]) ** 2 + (
        test[:, 1] - test_mot[:, 1]
    ) ** 2  ## Trenger ikke kvadratrot fordi alt annet enn 1 som verdi er ikke intakt(også raskere)
    if np.any(
        distance_array != 1
    ):  ## Hvis distance_array inneholder noe annet enn 1 er ikke polymer intakt
        return False
    return True

def check_if_intact_3(polymer: np.ndarray, polymer_length: int) -> bool:
    """Sjekker om en polymer er intakt

    Args:
        polymer (np.ndarray): polymeren som sjekkes
        polymer_length (int): lengden til polymeren

    Returns:
        bool: True hvis polymeren er intakt
    """
    if np.size(np.unique(polymer, axis=0), axis=0) != polymer_length:
        return False
        # Sjekker både at den har N monomerer og da har de alle en unik heltallsrepresentasjon

    test = polymer[1:]
    test_mot = polymer[:-1]
    distance_array = np.abs(test[:, 0] - test_mot[:, 0]) + np.abs(
        test[:, 1] - test_mot[:, 1]
    )  ## Trenger ikke kvadratrot fordi alt annet enn 1 som verdi er ikke intakt(også raskere)
    if np.any(
        distance_array != 1
    ):  ## Hvis distance_array inneholder noe annet enn 1 er ikke polymer intakt; Kan bli raskere???
        return False
    return True

@njit
def check_if_intact_4(polymer: np.ndarray, polymer_length: int) -> bool:
    """Checks if polymer is intact

    Args:
        polymer (np.ndarray): polymer that is checked
        polymer_length (int): length of the polymer

    Returns:
        bool: True if polymer is intact
    """

    if len(polymer) != polymer_length:
        return False

    unique_monomer = np.zeros_like(polymer)

    for i in range(polymer_length):  
        for j in range(i):
            if (                    #checks that the monomer's coordinates are not similar to a previous monomer
                polymer[i, 0] == unique_monomer[j, 0]
                and polymer[i, 1] == unique_monomer[j, 1]
            ):
                return False
            else:
                unique_monomer[i] = polymer[i]

    test = polymer[1:]
    test_against = polymer[:-1]
    distance_array = (test[:, 0] - test_against[:, 0]) ** 2 + (
        test[:, 1] - test_against[:, 1]
    ) ** 2  ## Don't need squareroot since all other values than 1 means that it is not intact 
    if np.any(
        distance_array != 1
    ):  ## If distance-array contains anything other than 1, the polymer is not intact
        return False
    return True

@njit
def rotate_polymer(
    polymer: np.ndarray, rotation_center: int, positive_direction: bool = True
) -> np.ndarray:
    """Rotates a polymer in the given direction around a monomer

    Args:
        polymer: A 2D numpy array with monomer coordinates

        rotation_center: Which monomer to rotate around
        `Note: It is not the index, but the monomer_number. [1, N]`

        positive_direction: Rotate in the positive direction if True, or negative direction if False

    Returns:
        a rotated copy of the polymer
    """
    # Make a slicing array to rotate the correct end of the polymer
    rotation_slice = np.full(len(polymer), False)

    # Choose to rotate the shortest tail of the polymer
    if rotation_center >= len(polymer) / 2:
        rotation_slice[rotation_center:] = True
    else:
        rotation_slice[:rotation_center] = True

    if positive_direction:
        direction = 1
    else:
        direction = -1

    # The coordinates in space of the rotation center
    rotation_position = polymer[rotation_center - 1]

    # Where _rel means the position relative to the rotation center
    # new_x = x_s + new_x_rel
    # new_y = y_s + new_y_rel
    # new_x_rel = - (y - y_s) * direction
    # new_y_rel = (x - x_s) * direction
    new_pos_rel = ((polymer[rotation_slice] - rotation_position) * direction)[:, ::-1]
    new_pos_rel[:, 0] *= -1  # Changes the sign of the x-values

    # Makes a copy of the polymer
    # TODO: om lettere å mutere for så å mutere tilbake igjen
    # new_polymer = polymer[:]
    new_polymer = polymer.copy()

    new_polymer[rotation_slice] = rotation_position + new_pos_rel
    return new_polymer


@njit
def rotate_polymer_mut(
    polymer: np.ndarray, rotation_center: int, positive_direction: bool = True
) -> np.ndarray:
    """Rotates a polymer in the given direction around a monomer

    Args:
        polymer: A 2D numpy array with monomer coordinates

        rotation_center: Which monomer to rotate around
        `Note: It is not the index, but the monomer_number. [1, N]`

        positive_direction: Rotate in the positive direction if True, or negative direction if False

    Returns:
        a rotated copy of the polymer
    """
    # Make a slicing array to rotate the correct end of the polymer
    rotation_slice = np.full(len(polymer), False)

    # Choose to rotate the shortest tail of the polymer
    if rotation_center >= len(polymer) / 2:
        rotation_slice[rotation_center:] = True
    else:
        rotation_slice[:rotation_center] = True

    if positive_direction:
        direction = 1
    else:
        direction = -1

    # The coordinates in space of the rotation center
    rotation_position = polymer[rotation_center - 1]

    # Where _rel means the position relative to the rotation center
    # new_x = x_s + new_x_rel
    # new_y = y_s + new_y_rel
    # new_x_rel = - (y - y_s) * direction
    # new_y_rel = (x - x_s) * direction
    new_pos_rel = ((polymer[rotation_slice] - rotation_position) * direction)[:, ::-1]
    new_pos_rel[:, 0] *= -1  # Changes the sign of the x-values

    polymer[rotation_slice] = rotation_position + new_pos_rel
    return polymer


@njit
def generate_flat_polymer(
    polymer_length: int, mid_of_polymer: np.ndarray = np.zeros(2)
) -> np.ndarray:
    """Generates a horizontal polymer with N (polymer_length) monomers

    Args:
        polymer_length (int): Number of monomers
        mid_of_polymer (np.ndarray, optional): The coordinates of the center monomer. Defaults to np.zeros(2).

    Returns:
        np.ndarray: the generated polymer
    """
    polymer_array = np.zeros((polymer_length, 2), dtype=np.int32)
    polymer_start = -int(polymer_length / 2) + mid_of_polymer[0]
    # + 1/2 to handle even numbers
    polymer_end = int((polymer_length + 1) / 2) + mid_of_polymer[0]
    polymer_array[:, 1] = mid_of_polymer[1]
    polymer_array[:, 0] = np.arange(polymer_start, polymer_end, 1, dtype=np.int32)

    return polymer_array


# The function can, (and should?), be JIT-compiled by numba.
@njit()
def calculate_energy(polymer: np.ndarray, V: np.ndarray) -> float:
    """Calculates the energy of the given polymer.

    Args:
        polymer: A 2D numpy array with monomer coordinates

        V: A matrix with the strength of the interaction between monomers of the given polymer.
        V[i, j] = V[j, i] = strength between monomer number (i+1) og (j+1)

    Returns:
        The energy of the polymer
    """
    N = len(polymer)
    # En matrise som angir om monomer (i+1) og (j+1) er naboer. b_matrix[i, j] = 1 dersom de er naboer.
    # A matrix which tells if monomer (i+1) and (j+1) are neighbours. b_matrix[i, j] = 1 if they are, and 0 otherwise
    b_matrix = np.zeros((N, N))
    for i in range(0, N):
        # We only need to look at monomer combinations which we have not checked yet.
        # The next monomer in the polymer does not interact with the given monomer, so we don't have to check it.
        for j in range(i + 2, N):
            # Kun nærmeste-nabo koordinater gir en euklidsk avstand på nøyaktig 1.
            # Only "closest neighbour"-coordinates will give a euclidic distance of exactly one, so we dont have to take the square root
            if np.sum((polymer[i] - polymer[j]) ** 2) == 1:
                # Trenger bare fylle den nedre trekanten av matrisen,
                # We only need to fill the lower triangle of the matrix, as it otherwise will just be symmetric.
                # Note: j>i
                b_matrix[j, i] = 1
    # Don't have to divide by two, because we only filled the lower triangle of the b_matrix.
    # We don't count the same interaction twice
    return float(np.sum(V * b_matrix))


@njit
def calculate_energy_2(polymer: np.ndarray, V: np.ndarray) -> float:
    """idk... think it works. maybe. see Oskar's notebook for details lol."""
    N = len(polymer)
    L = np.repeat(polymer, N).reshape(2 * N, N)
    b = np.where(
        ((L[::2] - L[::2].transpose()) ** 2 + (L[1::2] - L[1::2].transpose()) ** 2)
        == 1,
        1,
        0,
    )
    return 0.5 * (np.sum(V * b))


@njit
def calculate_energy_3(polymer: np.ndarray, V: np.ndarray) -> float:
    """idk... think it works. maybe. see Oskar's notebook for details lol."""
    N = len(polymer)
    L = np.repeat(polymer, N).reshape(2 * N, N)

    def is_neighbor(L):
        return (
            np.abs(L[::2] - L[::2].transpose()) + np.abs(L[1::2] - L[1::2].transpose())
            == 1
        )

    b = np.where(is_neighbor(L), 1, 0)
    return 0.5 * (np.sum(V * b))


if __name__ == "__main__":
    pol = np.array(
        [
            [-2, 0],
            [-1, 0],
            [0, 0],
            [1, 0],
            [1, -1],
            [0, -1],
            [-1, -1],
        ]
    )

    # dette er ikkje eit lovleg polymer.
    pol2 = np.array(
        [
            [0, 0],
            [1, 0],
            [1, 1],
            [0, 1],
            [-1, 1],
            [-2, 1],
            [-2, 0],
            [-1, 0],
            [-2, -1],
            [-2, -2],
            [-2, -3],
        ]
    )

    V = utilities.gen_V_matrix(7)
    V2 = utilities.gen_V_matrix(11, fill_value=1)
    print(calculate_energy_3(pol, V))
    print(calculate_energy_3(pol2, V2))
    print(calculate_energy_2(pol2, V2))
