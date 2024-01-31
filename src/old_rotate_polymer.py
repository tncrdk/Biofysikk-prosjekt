import numpy as np

def rotate_polymer(
    polymer: np.ndarray, rotation_center: int, positive_direction: bool = True
) -> np.ndarray:
    """Roterer et polymer rundt rotasjonssenteret i en gitt retning.

    Args:
        polymer: Et 2D numpy array med monomer-koordinatene

        rotation_center: Hvilket monomer man skal rotere rundt.
        Merk at det er ikke indeks, men polymer nummeret [1, N]

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

    # Med _rel menes posisjonen relativt rotasjonssenteret
    # new_x = x_s + new_x_rel
    # new_y = y_s + new_y_rel
    # new_x_rel = - (y - y_s) * direction
    # new_y_rel = (x - x_s) * direction
    new_pos_rel = ((polymer[rotation_center:] - rotation_position) * direction)[:, ::-1]
    new_pos_rel[:, 0] *= -1  # Endrer fortegnet til x-ene

    # Lager kopi av polymeret
    # om lettere å mutere for så å mutere tilbake igjen
    # new_polymer = polymer[:]
    new_polymer = polymer.copy()

    new_polymer[rotation_center:] = rotation_position + new_pos_rel
    return new_polymer


def rotate_polymer_head(
    polymer: np.ndarray, rotation_center: int, positive_direction: bool
) -> np.ndarray:
    """Roterer halen til et polymer rundt rotasjonssenteret i en gitt retning

    Args:
        polymer: Et 2D numpy array med monomer-koordinatene

        rotation_center: Hvilket monomer man skal rotere rundt.
        Merk: Det er ikke indeks, men polymer nummeret. [1, N]

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
    new_pos_rel = ((polymer[:rotation_center] - rotation_position) * direction)[:, ::-1]
    new_pos_rel[:, 0] *= -1

    # om lettere å mutere for så å mutere tilbake igjen
    # new_polymer = polymer[:]
    new_polymer = polymer.copy()

    new_polymer[:rotation_center] = rotation_position + new_pos_rel
    return new_polymer
