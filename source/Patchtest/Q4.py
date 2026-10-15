
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from source.elements.elements import quad_q4
from source.materials.linear_elastic import LinearElastic
from source.solvers.solvers import ke, node_stress, gauss_stress

# Coordenadas de los nodos
coordnodos = np.array([
    [0, 0],       # 0
    [100, 0],     # 1
    [100, 100],   # 2
    [0, 100],     # 3
    [45, 55],     # 4 interior
    [50, 0],      # 5
    [100, 50],    # 6
    [50, 100],    # 7
    [0, 50]       # 8
], dtype=float)

# Conectividades
nodos = [
    [0, 5, 4, 8],
    [5, 1, 6, 4],
    [4, 6, 2, 7],
    [8, 4, 7, 3]
]

# Material
material = LinearElastic(21000, 0.25)
tipo = "Plane_strain"
espesor = 1.0

# Elementos Q4
elementos = [
    quad_q4(1, [0, 5, 4, 8], material, tipo),
    quad_q4(2, [5, 1, 6, 4], material, tipo),
    quad_q4(3, [4, 6, 2, 7], material, tipo),
    quad_q4(4, [8, 4, 7, 3], material, tipo)
]

# Grados de libertad por nodo
n_nodos = len(coordnodos)
n_gdl = 2 * n_nodos

eledof = np.arange(n_gdl).reshape(n_nodos, 2)

# Matriz de rigidez global
K = np.zeros((n_gdl, n_gdl))


for elemento in elementos:

    node_id = elemento.node
    coordenadas = coordnodos[node_id]


    Ke = espesor * ke(elemento, 2, coordenadas)


    dofs = eledof[node_id].flatten()

 
    for i in range(elemento.ndofs):
        for j in range(elemento.ndofs):
            K[dofs[i], dofs[j]] += Ke[i, j]



F = np.zeros(n_gdl)

F[2*1] = 0.5   # Nodo 1
F[2*6] = 1.0   # Nodo 6
F[2*2] = 0.5   # Nodo 2
F[2*3] = -0.5
# Condiciones de borde
restringidos = [0, 1, 16]

libres = [
    i for i in range(n_gdl)
    if i not in restringidos
]

# Sistema reducido
Kr = K[np.ix_(libres, libres)]

D = np.zeros(n_gdl)

# Resolver desplazamientos
D[libres] = np.linalg.solve(Kr, F[libres])

print("Desplazamientos:")
print(D)


# Tensiones nodales promediadas
connexiones = []

for elemento in elementos:
    connexiones.append(elemento.node)

sigma = np.zeros((3, n_nodos))

# Estado teorico del Patch Test
sigma_teorica = np.array([0.02, 0, 0])

pasa_patch = True

for elemento in elementos:

    node_id = elemento.node
    coordenadas = coordnodos[node_id]

    dofs = eledof[node_id].flatten()
    De = D[dofs]

    
    sigma_gauss, q = gauss_stress(
        elemento, 2, coordenadas, De
    )

    print("\nElemento:", elemento.id)
    print("Tensiones en Gauss:")
    print(sigma_gauss)

    # Verificar cada elemento sin promediar
    if np.allclose(
        sigma_gauss,
        sigma_teorica[:, None],
        atol=1e-10,
        rtol=1e-8
    ):
        print("PASA")
    else:
        print("NO PASA")
        pasa_patch = False

    sigma_elemento = node_stress(
        elemento, coordenadas, De
    )


    for j in range(elemento.n_node):
        div = sum(
            sublist.count(node_id[j])
            for sublist in connexiones
        )

        sigma[:, node_id[j]] += sigma_elemento[:, j] / div


print("\nTensiones nodales promediadas:")
print(sigma)

# Comprobar equilibrio de los GDL libres
residuo = K @ D - F

equilibrio = np.allclose(
    residuo[libres], 0, atol=1e-9
)

print("\nEquilibrio en GDL libres:", equilibrio)

if pasa_patch and equilibrio:
    print("\nPATCH TEST Q4: APROBADO")
else:
    print("\nPATCH TEST Q4: NO APROBADO")
