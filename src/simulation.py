import numpy as np
import polymer
import visualization


def alg1(N: int, Ns: int) -> tuple[np.ndarray, int]:
    """Implementation of algorithm 1.
    ---
    Args:
        N: length of polymer.
        Ns: number of twists (attempts) to be performed.

    Returns:
        (polymer, counter)
            polymer: polymer.
            counter: number of successful twists.
    """
    counter = 1
    pol = polymer.generate_flat_polymer(N)
    for i in range(Ns):
        # random monomer and random twisting direction
        rnd_monomer = np.random.randint(2, N)
        rnd_rotate = bool(int(np.random.uniform() + 0.5))

        # TODO: possible to mutate the same array instead of copying?
        twisted_pol = polymer.rotate_polymer(pol, rnd_monomer, rnd_rotate)
        if polymer.check_if_intact_2(twisted_pol, N):
            counter += 1
            pol = twisted_pol

    return pol, counter


if __name__ == "__main__":
    N = 100
    Ns = 100
    pol, _ = alg1(N, Ns)
    visualization.illustrate_polymer(
        pol, cmap="Greens", title=f"length: {N}, # twists: {Ns}"
    )
    pol, _ = alg1(N, Ns)
    visualization.illustrate_polymer(
        pol, cmap="Blues", title=f"length: {N}, # twists: {Ns}"
    )
