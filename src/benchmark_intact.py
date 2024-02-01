import numpy as np
from numba import njit
import benchmarks
import polymer

@njit
def unique(polymer:np.ndarray, polymer_length:int) -> bool:
    unique_monomer = np.zeros_like(polymer)

    for i in range(polymer_length):  ## den sjekker at alle monomerene er unike
        for j in range(i):

            if polymer[i, 0] == unique_monomer[j, 0] and polymer[i, 1] == unique_monomer[j, 1]:
                return False
            else:
                unique_monomer[i] = polymer[i]  
    return True


def check_if_intact_2(polymer: np.ndarray, polymer_length: int) -> bool:
    """Sjekker om en polymer er intakt

    Args:
        polymer (np.ndarray): polymeren som sjekkes
        polymer_length (int): lengden til polymeren

    Returns:
        bool: True hvis polymeren er intakt
    """
    unique_monomers = np.unique(polymer, axis=0)
    if np.size(unique_monomers, axis=0) != polymer_length:
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


def check_if_intact_setup():
    return (polymer.generate_flat_polymer(35),)


@benchmarks.benchmark(iterations=1000, setup_func=check_if_intact_setup, warmup=1)
def bench_check_if_intact(a: np.ndarray, b: int = 35):
    return unique(a, b)


time, results = bench_check_if_intact()
print(time)
