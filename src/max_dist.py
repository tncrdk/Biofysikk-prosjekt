from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist, pdist
import time
import numpy as np
import polymer 
import math
import benchmarks
from numba import njit
from fastdist import fastdist

def calculate_diameter(polymer:np.ndarray) -> float:
    hull = ConvexHull(polymer)
    hull_points = polymer[hull.vertices, :]

    distances = cdist(hull_points, hull_points, metric = 'euclidean')

    return np.max(distances)

@njit
def calculate_diameter_2(polymer:np.ndarray) -> float:
    answer = 0
    for i in range(len(polymer)):
        for j in range(i + 1, len(polymer)):
            #distance_between = np.sqrt(np.sum(np.square(polymer[i] - polymer[j])))
            distance_between = fastdist.euclidean(polymer[i], polymer[j])
            if  distance_between > answer:
                answer = distance_between

    return answer

def calculate_diameter_3(polymer:np.ndarray) -> float:
    return np.max(cdist(polymer, polymer))


def calculate_diameter_4(polymer:np.ndarray) -> float:
    return np.max(pdist(polymer))


def calculate_diameter_5(polymer:np.ndarray) -> float:
    return np.max(fastdist.matrix_pairwise_distance(polymer, fastdist.euclidean, "euclidean"))


#a = np.random.randint(0, 35, (35, 2))
# start = time.time()
# print(calculate_diameter(a))
# end = time.time()
# print(end - start)

    
# b = polymer.generate_flat_polymer(35)
# start_b = time.time()
# print(calculate_diameter_2(a))
# end_b = time.time()

# print(end_b - start_b)

# start_c = time.time()
# print(calculate_diameter_3(a))
# end_c = time.time()

# print(end_c - start_c)

def calculate_diameter_setup():
    out = np.array(np.random.randint(0, 35, (35,2)))
    return (out, )

@benchmarks.benchmark(iterations=1000, setup_func=calculate_diameter_setup, warmup=1)
def bench_dia(polymer:np.ndarray):
    return calculate_diameter_5(polymer)


print(bench_dia()[0])