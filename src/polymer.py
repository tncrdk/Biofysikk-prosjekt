import numpy as np
import time


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
    """Sjekker om en polymer er intakt

    Args:
        polymer (np.ndarray): polymeren som sjekkes
        polymer_length (int): lengden til polymeren

    Returns:
        bool: True hvis polymeren er intakt
    """

    if (
        np.size(polymer, axis=0) != polymer_length
    ):  # Sjekker både at den har N monomerer og da har de alle en unik heltallsrepresentasjon
        return False

    for i in range(1, polymer_length):
        distance = (polymer[i - 1, 0] - polymer[i, 0]) ** 2 + (
            polymer[i - 1, 1] - polymer[i, 1]
        ) ** 2  ## Trenger ikke kvadratrot fordi alt annet enn 1 som verdi er ikke intakt(også raskere)
        distance = (polymer[i - 1, 0] - polymer[i, 0]) ** 2 + (
            polymer[i - 1, 1] - polymer[i, 1]
        ) ** 2  ## Trenger ikke kvadratrot fordi alt annet enn 1 som verdi er ikke intakt(også raskere)
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

    if (
        np.size(polymer, axis=0) != polymer_length
    ):  # Sjekker både at den har N monomerer og da har de alle en unik heltallsrepresentasjon
        return False

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


def rotate_polymer(
    polymer: np.ndarray, rotation_center: int, positive_direction: bool = True
) -> np.ndarray:
    """Roterer et polymer rundt rotasjonssenteret i en gitt retning.

    Args:
        polymer: Et 2D numpy array med monomer-koordinatene

        rotation_center: Hvilket monomer man skal rotere rundt.
        Merk: Det er ikke indeks, men polymer nummeret [1, N]

        positive_direction: Rotere i positiv retning. Hvis False roteres det i negativ retning.

    Returns:
        en rotert kopi av polymeret
    """

    if (
        rotation_center >= len(polymer) / 2
    ):  # Velger å rotere den delen av polymeret som er kortest
        return rotate_polymer_tail(polymer, rotation_center, positive_direction)
    else:
        return rotate_polymer_head(polymer, rotation_center, positive_direction)


def rotate_polymer_tail(
    polymer: np.ndarray, rotation_center: int, positive_direction: bool
) -> np.ndarray:
    """Roterer halen til et polymer rundt rotasjonssenteret i en gitt retning

    Args:
        polymer: Et 2D numpy array med monomer-koordinatene

        rotation_center: Hvilket monomer man skal rotere rundt.
        Merk: Det er ikke indeks, men polymer nummeret [1, N]

        positive_direction: Rotere i positiv retning. Hvis False roteres det i negativ retning.

    Returns:
        en rotert kopi av polymeret
    """
    if positive_direction:
        direction = 1
    else:
        direction = -1

    # Rotasjonssentrum sine koordinater
    rotation_position = polymer[rotation_center - 1]

    # new_x = x_s + x_rel
    # x_rel = - (y - y_s) * direction
    # new_y = y_s + y_rel
    # y_rel = (x - x_s) * direction
    new_pos_rel = (polymer[rotation_center:] - rotation_position) * direction
    new_pos_rel[:, 1] *= -1
    # Lager kopi av polymeret
    new_polymer = polymer[:]
    new_polymer[rotation_center:] = rotation_position + new_pos_rel[:, ::-1]
    return new_polymer


def rotate_polymer_head(
    polymer: np.ndarray, rotation_center: int, positive_direction: bool
) -> np.ndarray:
    """Roterer halen til et polymer rundt rotasjonssenteret i en gitt retning

    Args:
        polymer: Et 2D numpy array med monomer-koordinatene

        rotation_center: Hvilket monomer man skal rotere rundt.
        Merk: Det er ikke indeks, men polymer nummeret [1, N]

        positive_direction: Rotere i positiv retning. Hvis False roteres det i negativ retning.

    Returns:
        en rotert kopi av polymeret
    """
    # For kommentarer se rotate_polymer_tail
    if positive_direction:
        direction = 1
    else:
        direction = -1

    rotation_position = polymer[rotation_center - 1]
    new_pos_rel = (polymer[:rotation_center] - rotation_position) * direction
    new_pos_rel[:, 1] *= -1
    new_polymer = polymer[:]
    new_polymer[:rotation_center] = rotation_position + new_pos_rel[:, ::-1]
    return new_polymer


def generate_flat_polymer(
    polymer_length: int, mid_of_polymer: np.ndarray = np.zeros(2)
) -> np.ndarray:
    """Genererer en horisontal polymer med N monomerer

    Args:
        polymer_length (int): N antall monomerer
        mid_of_polymer (np.ndarray, optional): Midtpunktet til polymeren, Defaults to np.zeros(2).

    Returns:
        np.ndarray: den genererte polymeren
    """ """"""
    polymer_array = np.zeros((polymer_length, 2), dtype = np.int32)
    polymer_start = -int(polymer_length / 2) + mid_of_polymer[0]
    # + 1/2 for å håndtere partall
    polymer_end = int((polymer_length + 1) / 2) + mid_of_polymer[0]
    polymer_array[:, 1] = mid_of_polymer[1]
    polymer_array[:, 0] = np.arange(polymer_start, polymer_end, 1, dtype = np.int32)

    return polymer_array


print(generate_flat_polymer(3))
