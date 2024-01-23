from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
import time
import numpy as np

def calculate_diameter(polymer:np.ndarray) -> float:
    hull = ConvexHull(polymer)
    hull_points = polymer[hull.vertices, :]

    distances = cdist(hull_points, hull_points, metric = 'euclidean')

    return np.max(distances)
    
a = np.random.randint(0, 2000, (35, 2))
start = time.time()
print(calculate_diameter(a))
end = time.time()
print(end - start)