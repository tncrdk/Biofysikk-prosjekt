import numpy as np

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

    if rotation_center >= len(polymer) / 2:
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

    rotation_position = polymer[rotation_center - 1]
    new_pos_rel = (polymer[rotation_center:] - rotation_position) * direction
    new_pos_rel[:, 1] *= -1
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
