import numpy as np
import polymer
import visualization
from scipy.constants import Boltzmann


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


def metropolis(
    pol: np.ndarray, N_s: int, V: np.ndarray, T: float
) -> tuple[np.ndarray, np.ndarray]:
    """Kjører metropolis-algoritmen som beskrevet i oppgaveheftet

    Args:
        pol: polymer starttilstand
        N_s: Antall forsøk på rotasjon
        V: Vekselvirkningene mellom to gitte monomerer
        T: temperaturen i kelvin

    Returns:
        (Siste polymer , array med alle energiene som ble simulert)
    """
    E_array = np.zeros(N_s)
    N = len(pol)
    E = polymer.calculate_energy(pol, V)
    i = 0
    while i < N_s - 1:
        # random monomer and random twisting direction
        rnd_monomer = np.random.randint(2, N)
        rnd_rotate = bool(int(np.random.uniform() + 0.5))

        # TODO: possible to mutate the same array instead of copying?
        twisted_pol = polymer.rotate_polymer(pol, rnd_monomer, rnd_rotate)
        if polymer.check_if_intact_2(twisted_pol, N):
            i += 1
            E_new = polymer.calculate_energy(twisted_pol, V)

            if E_new < E:
                pol = twisted_pol
                E = E_new
            # TODO: Bruke en annen distribusjon enn uniform?
            # TODO: Boltzmann-konstanten er liten. Sjekk at python håndterer det.
            elif np.random.uniform() < np.exp(-(E_new - E) / (T * Boltzmann)):
                pol = twisted_pol
                E = E_new
            E_array[i] = E

    return pol, E_array


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
