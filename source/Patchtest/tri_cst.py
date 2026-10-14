
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from source.elements.elements import tri_cst
from source.materials.linear_elastic import LinearElastic
from source.solvers.solvers import jacobiano, b_matrix, ke

# Nodos
nodos = np.array([
    [0, 0],
    [100, 0],
    [100, 100],
    [0, 100],
    [45, 55]
], dtype=float)

material = LinearElastic(21000, 0.25)
tipo = "Plane_stress"

# Elementos
elementos = [
    tri_cst(1, [0, 1, 4], material, tipo),
    tri_cst(2, [4, 1, 2], material, tipo),
    tri_cst(3, [2, 3, 4], material, tipo),
    tri_cst(4, [3, 0, 4], material, tipo)
]

# Matriz de rigidez global
K = np.zeros((10, 10))

matrices_B = []

# Ensamblaje
for elemento in elementos:

    coordenadas = nodos[elemento.node]

    # Jacobiano y matriz B
    J, detJ, invJ = jacobiano(elemento, 0, 0, coordenadas)
    B = b_matrix(elemento, 0, 0, invJ)

    matrices_B.append(B)

    # Matriz de rigidez elemental
    Ke = ke(elemento, 1, coordenadas)

    # Grados de libertad globales
    gdl = []

    for nodo in elemento.node:
        gdl.append(2*nodo)
        gdl.append(2*nodo + 1)

    # Ensamblaje global
    for i in range(6):
        for j in range(6):
            K[gdl[i], gdl[j]] += Ke[i, j]

# Condiciones de borde
D = np.array([0, 0, 0.1, 0, 0.1, 0, 0, 0, 0, 0], dtype=float)

restringidos = [0, 1, 2, 3, 4, 5, 6, 7]
libres = [8, 9]

# Sistema reducido
Kr = K[np.ix_(libres, libres)]

Fr = -K[np.ix_(libres, restringidos)] @ D[restringidos]

# Resolver desplazamientos
D[libres] = np.linalg.solve(Kr, Fr)

print("Desplazamientos del nodo interior:", D[8], D[9])

# Verificar deformaciones
for n, elemento in enumerate(elementos):

    gdl = []

    for nodo in elemento.node:
        gdl.append(2*nodo)
        gdl.append(2*nodo + 1)

    De = D[gdl]

    B = matrices_B[n]

    epsilon = B @ De

    print("Elemento:", elemento.id)
    print("Deformaciones:", epsilon)
