import polymer
import visualization
import numpy as np
import simulation
import utilities

"""
Tests for polymer.py
"""


def test_generate_flat_polymer():
    """simple test for generate_flat_polymer"""
    expected_result = np.array([[-1, 0], [0, 0], [1, 0]], dtype=np.int32)
    res = polymer.generate_flat_polymer(3)
    assert np.all(
        np.equal(res, expected_result)
    ), f"Expected\n\t{expected_result}\nGot\n\t{res}"


def test_check_if_intact2():
    """generates 50 random polymers of size 3.
    If check_if_intact2 returns True one can validate this"""
    test_cases = [np.random.randint(-5, 5, (3, 2)) for _ in range(50)]
    for test_case in test_cases:
        if polymer.check_if_intact_2(test_case, len(test_case)):
            print(f"The polymer\n\t{test_case}\nis intact.\nPlease check")
            visualization.illustrate_polymer(test_case)


def test_check_if_intact():
    """generates 50 random polymers of size 3.
    If check_if_intact returns True one can validate this"""
    test_cases = [np.random.randint(-5, 5, (3, 2)) for _ in range(50)]
    for test_case in test_cases:
        if polymer.check_if_intact(test_case, len(test_case)):
            print(f"The polymer\n\t{test_case}\nis intact.\nPlease check")
            visualization.illustrate_polymer(test_case)


def test_check_if_intact_explicit():
    """Checks three polymers that is defect."""
    a = np.array([[-1, -4], [0, -1], [-3, 0]])
    b = np.array([[-3, -5], [-2, -1], [-1, 3]])
    c = np.array([[3, -4], [-1, 0], [-5, -4]])
    c = np.array([[3, -4], [-1, 0], [-5, -4]])
    d = [a, b, c]
    for i in d:
        print(f"check_if_intact gir {polymer.check_if_intact(i, len(i))}")
        print(f"check_if_intact2 gir {polymer.check_if_intact_2(i, len(i))}")


def test_rotate_polymer():
    """does some rotations and prints the result"""
    a = np.array([[i, 0] for i in range(15)])
    visualization.illustrate_polymer(a)
    a = polymer.rotate_polymer(a, 9)
    visualization.illustrate_polymer(a)
    a = polymer.rotate_polymer(a, 7, False)
    visualization.illustrate_polymer(a)
    a = polymer.rotate_polymer(a, 10, False)
    visualization.illustrate_polymer(a)
    a = polymer.rotate_polymer(a, 3, False)
    visualization.illustrate_polymer(a)
    a = polymer.rotate_polymer(a, 4, False)
    visualization.illustrate_polymer(a)
    a = polymer.rotate_polymer(a, 5)
    visualization.illustrate_polymer(a)


def test_calculate_energy():
    """runs the calculate_energy function against
    pre-calculated energy values for a handful of polymers."""
    
    V4_weird = utilities.gen_V_matrix(4, fill_value=-1)
    V4_weird[-1, 0] = -2
    V4_weird[0, -1] = -2

    V6 = utilities.gen_V_matrix(6, fill_value=-1)
    V8_neg = utilities.gen_V_matrix(8, fill_value=-2)
    V8_pos = utilities.gen_V_matrix(8, fill_value=2)
    
    test_cases = [
        # (V matrix, polymer, expected energy)
        (V4_weird, np.array([[-1,0], [0,0], [0,1], [-1,1]]), -2),

        (V6, polymer.generate_flat_polymer(6), 0),
        (V6, np.array([[0, 0], [1, 0], [1, 1], [0, 1], [-1, 1], [-2, 1]]), -1),
        (V6, np.array([[0, 0], [1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0]]), -2),
        
        (V8_neg, polymer.generate_flat_polymer(8), 0),
        (V8_neg, np.array([[0, 0], [1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1]]), -6),
        (V8_neg, np.array([[0, 0], [1, 0], [1, 1], [0, 1], [-1, 1], [-1, 2], [0, 2], [1, 2]]), -6),

        (V8_pos, polymer.generate_flat_polymer(8), 0),
        (V8_pos, np.array([[0, 0], [1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1]]), 6),
        (V8_pos, np.array([[0, 0], [1, 0], [1, 1], [0, 1], [-1, 1], [-1, 2], [0, 2], [1, 2]]), 6)
    ]

    for test_case in test_cases:
        energy = polymer.calculate_energy(test_case[1], test_case[0])
        # energy = polymer.calculate_energy_2(test_case[1], test_case[0])
        # energy = polymer.calculate_energy_3(test_case[1], test_case[0])
        expected_energy = test_case[2]
        assert  energy == expected_energy, f"{test_case[1]} has energy {energy}, but expected {expected_energy}"


"""
Tests for visualization.py
"""


def test_visualization():
    """simple function for visualizing some polymers"""
    a = np.array([[0, 0], [1, 0], [2, 0], [3, 0]])
    b = np.array([[0, 0], [0, 1], [0, 2], [0, 3]])
    c = np.array([[0, 0], [-1, 0], [-2, 0], [-3, 0]])
    d = np.array([[0, 0], [0, -1], [0, -2], [0, -3]])
    e = np.array([[1, 0], [2, 0], [3, 0]])
    f = np.array([[0, 1], [0, 2], [0, 3]])
    g = np.array([[-1, 0], [-2, 0], [-3, 0]])
    h = np.array([[0, -1], [0, -2], [0, -3]])
    m = np.array([[12, 13], [15, 16]])
    polymers_to_test = [m, a, b, c, d, e, f, g, h]

    for p in polymers_to_test:
        print(p)
        visualization.illustrate_polymer(p, title=str(p))


"""
Tester for simulation.py
"""


def test_metropolis():
    V = utilities.gen_V_matrix(11, fill_value=-1)
    T = 1000
    pol = polymer.generate_flat_polymer(11)
    N_s = 1000
    pol, E_array = simulation.metropolis(pol, N_s, V, T)
    visualization.illustrate_polymer(pol, numbers=True)
    print(E_array)


if __name__ == "__main__":
    tests = [
        # test_generate_flat_polymer,
        # test_check_if_intact2,
        # test_check_if_intact,
        # test_check_if_intact_explicit,
        # test_visualization,
        test_rotate_polymer,
        # test_calculate_energy,
        # test_metropolis,
        # test_calculate_energy,
    ]

    for i, test in enumerate(tests):
        test()
