import numpy as np
from time import perf_counter
from typing import Callable
import polymer


def benchmark(iterations: int = 10, setup_func=lambda x: None, warmup = 1) -> Callable:
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


@benchmark(iterations=1000, setup_func=calculate_energy_setup, warmup=1)
def bench_calculate_energy(a: np.ndarray, V: np.ndarray):
    return polymer.calculate_energy(a, V)


if __name__ == "__main__":
    benchmarks = [bench_calculate_energy]
    for b in benchmarks:
        time, results = b()
        print(time)
