import numpy as np
import polymer
import benchmarks

a = np.array(([1,2], [3,3], [4,5], [3,6]))

def unique(polymer:np.ndarray, polymer_length:int) -> bool:
    unique_monomer = np.zeros_like(polymer)

    for i in range(polymer_length):  ## den sjekker at alle monomerene er unike
        for j in range(i):

            if polymer[i, 0] == unique_monomer[j, 0] and polymer[i, 1] == unique_monomer[j, 1]:
                return False
            else:
                unique_monomer[i] = polymer[i]  
    return True


def check_if_intact_4(polymer: np.ndarray, polymer_length: int) -> bool:
    """Sjekker om en polymer er intakt

    Args:
        polymer (np.ndarray): polymeren som sjekkes
        polymer_length (int): lengden til polymeren

    Returns:
        bool: True hvis polymeren er intakt
    """

    if len(polymer) != polymer_length:
        return False
    
    unique_monomer = np.zeros_like(polymer)

    for i in range(polymer_length):  ## den sjekker at alle monomerene er unike
        for j in range(i):

            if polymer[i, 0] == unique_monomer[j, 0] and polymer[i, 1] == unique_monomer[j, 1]:
                return False
            else:
                unique_monomer[i] = polymer[i]  
        
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


