import numpy as np
import polymer


def alg1(N: int, Ns: int) -> tuple[np.ndarray, int]:
    """
    Args:
        N: length of polymer
        Ns: number of twists (attempts) to be performed

    Returns:
        (polymer, counter)
            polymer: polymer.
            counter: number of succesful twists.

    """
    counter = 1
    pol = polymer.generate_flat_polymer(N)
    for i in range(Ns):
        rnd_monomer = np.random.randint(2, N)
        rnd_rotate = bool(int(np.random.uniform() + 0.5))
        twisted_pol = polymer.rotate_polymer(pol, rnd_monomer, rnd_rotate)
        if polymer.check_if_intact_2(twisted_pol, N):
            counter += 1
            pol = twisted_pol
    return pol, counter
