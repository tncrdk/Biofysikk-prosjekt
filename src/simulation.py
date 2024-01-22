import numpy as np
import polymer
import visualization


def alg1(N: int, Ns: int) -> tuple[np.ndarray, int]:
    """
    Args:
        N:
            length of polymer

        Ns:
            number of twists (attempts) to be performed

    Returns: (polymer, counter)
        polymer: polymer.
        counter: number of succesful twists.

    """
    counter = 1
    pol = polymer.generate_flat_polymer(N)
    for i in range(Ns):
        rnd_monomer = np.random.randint(0, N)
        rnd_rotate = True if np.random.randint(0, 2) == 1 else False
        twisted_pol = polymer.rotate_polymer(pol, rnd_monomer, rnd_rotate)
        if polymer.check_if_intact_2(twisted_pol, N):
            counter += 1
            pol = twisted_pol
    return pol, counter


"""
counter ← 1
polymer ← rett horisontal polymer
for
i = 1 to Ns do
Velg en tilfeldig monomer ˚a forsøke ˚a rotere om
Velg tilfeldig om det skal forsøkes ˚a rotere med eller mot klokken.
twisted polymer ← twist(polymer, random monomer, twist clockwise)
if check if legal twist(twisted polymer) then
counter = counter + 1
polymer ← twisted polymer
end if
end for
return polymer, counter
"""

if __name__ == "__main__":
    pol, _ = alg1(20, 10)
    visualization.illustrate_polymer(pol)
