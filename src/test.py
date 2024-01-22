import polymer
import visualization
import numpy as np


def test_generate_flat_polymer():
    """simple test for generate_flat_polymer"""
    expected_result = np.array([[-1, 0], [0, 0], [1, 0]])
    res = polymer.generate_flat_polymer(3)
    assert not np.all(
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
    d = [a, b, c]
    for i in d:
        print(f"check_if_intact gir {polymer.check_if_intact(i, len(i))}")
        print(f"check_if_intact2 gir {polymer.check_if_intact_2(i, len(i))}")


if __name__ == "__main__":
    tests = [
        test_generate_flat_polymer,
        test_check_if_intact2,
        test_check_if_intact,
        test_check_if_intact_explicit,
    ]

    for i, test in enumerate(tests):
        test()
