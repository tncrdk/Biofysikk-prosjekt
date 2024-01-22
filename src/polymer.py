import numpy as np
import time

a = np.array([[-1, -4], [0, -1], [-3, 0]])
b = np.array([[-3, -5], [-2, -1], [-1, 3]])
c = np.array([[3, -4], [-1, 0], [-5, -4]])
d = np.array([[1, 0], [1, 1], [1, 2], [1, 1]])


def check_if_intact(polymer: np.ndarray, polymer_length: int) -> bool:
    """Sjekker om en polymer er intakt

    Args:
        polymer (np.ndarray): polymeren som sjekkes
        polymer_length (int): lengden til polymeren

    Returns:
        bool: True hvis polymeren er intakt
    """
    if np.size(np.unique(polymer, axis = 0), axis = 0) != polymer_length:
        return False 
    # Sjekker både at den har N monomerer og da har de alle en unik heltallsrepresentasjon
    for i in range(1, polymer_length):
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
    if np.size(np.unique(polymer, axis = 0), axis = 0) != polymer_length:
        return False 
        # Sjekker både at den har N monomerer og da har de alle en unik heltallsrepresentasjon

    test = polymer[1:]
    test_mot = polymer[:-1]
    distance_array = (test[:, 0] - test_mot[:, 0]) ** 2 + (
        test[:, 1] - test_mot[:, 1]
    ) ** 2  ## Trenger ikke kvadratrot fordi alt annet enn 1 som verdi er ikke intakt(også raskere)
    if np.any(
        distance_array != 1
    ):  ## Hvis distance_array inneholder noe annet enn 1 er ikke polymer intakt; Kan bli raskere???
        return False
    return True

print(check_if_intact(d, 4))