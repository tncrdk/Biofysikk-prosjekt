import numpy as np
from time import perf_counter
from typing import Callable
import polymer


def benchmark(iterations: int = 10, setup_func=lambda x: None, warmup=1) -> Callable:
    def decorator(func):
        if setup_func is None:
            _ = [func() for _ in range(warmup)]

            def none_wrapper():
                stime = perf_counter()
                results = [0] * iterations
                for i in range(iterations):
                    results[i] = func()
                tot_time = perf_counter() - stime
                avg_time = tot_time / iterations
                return avg_time, results

            return none_wrapper
        else:
            _ = [func(*setup_func()) for _ in range(warmup)]

            def wrapper():
                init_setup = setup_func()
                stime = perf_counter()
                results = [0] * iterations
                for i in range(iterations):
                    results[i] = func(*init_setup)
                tot_time = perf_counter() - stime
                avg_time = tot_time / iterations
                return avg_time, results

            return wrapper

    return decorator


def calculate_energy_setup():
    V = np.array(
        [
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
        ]
    )

    a = np.array(
        [
            [0, 0],
            [1, 0],
            [1, 1],
            [0, 1],
            [-1, 1],
            [-2, 1],
            [-2, 0],
            [-1, 0],
            [-2, -1],
            [-2, -2],
            [-2, -3],
        ]
    )
    return a, V


def bench_rotate_polymer_setup():
    return (polymer.generate_flat_polymer(35), )


@benchmark(iterations=1000, setup_func=calculate_energy_setup, warmup=1)
def bench_calculate_energy(a: np.ndarray, V: np.ndarray):
    return polymer.calculate_energy(a, V)


@benchmark(iterations=1000, setup_func=calculate_energy_setup, warmup=1)
def bench_calculate_energy_2(a: np.ndarray, V: np.ndarray):
    return polymer.calculate_energy_2(a, V)

@benchmark(iterations=1000, setup_func=calculate_energy_setup, warmup=1)
def bench_calculate_energy_3(a: np.ndarray, V: np.ndarray):
    return polymer.calculate_energy_3(a, V)

@benchmark(iterations=1000, setup_func=bench_rotate_polymer_setup)
def bench_rotate_polymer(pol: np.ndarray):
    a = polymer.rotate_polymer(pol, 15)
    a = polymer.rotate_polymer(a, 10)
    a = polymer.rotate_polymer(a, 3, False)
    a = polymer.rotate_polymer(a, 20)
    a = polymer.rotate_polymer(a, 15, False)
    a = polymer.rotate_polymer(a, 25)

if __name__ == "__main__":
    # time1, results1 = bench_calculate_energy()
    # time2, results2 = bench_calculate_energy_2()
    # time3, results3 = bench_calculate_energy_3()
    #
    # print(f"func: bench_calculate_energy\ntime: {time1:>10.5e}\nresults[0]: {results1[0]}\n")
    # print(f"func: bench_calculate_energy2\ntime: {time2:>10.5e}\nresults[0]: {results2[0]}\n")
    # print(f"func: bench_calculate_energy3\ntime: {time3:>10.5e}\nresults[0]: {results3[0]}\n")

    time, _ = bench_rotate_polymer()
    print(f"Time: {time/6}")
    # benchmarks = [bench_calculate_energy]
    # for b in benchmarks:
    #     time, results = b()
    #     print(f"function:{b.__name__}\ntime: {time:>10.5e}")
