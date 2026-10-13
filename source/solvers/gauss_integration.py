import math
import numpy as np

def gauss(element, n):
    if element.n_node == 3 or element.n_node == 6:
        if n == 1:
            return ((1/3, 1/3),), (0.5,)
        elif n == 2:
            return ((2/3, 1/6), (1/6, 2/3), (1/6, 1/6)), (1/3, 1/3, 1/3)

    if element.n_node == 4 or element.n_node == 8:
        if n == 1:
            return ((0, 0),), (4,)
        elif n == 2:
            a = 1 / math.sqrt(3)
            return ((a, a),
                    (-a, a),
                    (-a, -a),
                    (a, -a)), (1, 1, 1, 1)
        elif n == 3:
            a = math.sqrt(0.6)
            pts = [-a, 0.0, a]
            w_1d = [5/9, 8/9, 5/9]
            
            coords = []
            weights = []
            for i, r in enumerate(pts):
                for j, s in enumerate(pts):
                    coords.append((r, s))
                    weights.append(w_1d[i] * w_1d[j])
                    
            return tuple(coords), tuple(weights)
        

         
    

