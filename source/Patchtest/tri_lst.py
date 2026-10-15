
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from source.elements.elements import tri_lst
from source.materials.linear_elastic import LinearElastic
from source.solvers.solvers import ke, node_stress, gauss_stress

nodos = np.array([
    [0, 0],
    [100, 0],
    [100, 100],
    [0, 100],
    [45, 55],
    [50, 0],
    [72.5, 27.5],
    [22.5, 27.5],
    [100, 50],
    [72.5, 77.5],
    [50, 100],
    [22.5, 77.5],
    [0, 50]
], dtype=float)

conectividades = [
    [0, 1, 4, 5, 6, 7],
    [4, 1, 2, 6, 8, 9],
    [2, 3, 4, 10, 11, 9],
    [3, 0, 4, 12, 7, 11]
]

E = 21000
nu = 0.25
espesor = 1.0
sigma_x = 0.02

material = LinearElastic(E, nu)
tipo = "Plane_stress"

elementos = [
    tri_lst(i + 1, conectividades[i], material, tipo)
    for i in range(4)
]

n_nodos = len(nodos)
n_gdl = 2 * n_nodos

eledof = np.arange(n_gdl).reshape(n_nodos, 2)

K = np.zeros((n_gdl, n_gdl))

for elemento in elementos:
    node_id = elemento.node
    coordenadas = nodos[node_id]

    Ke = espesor * ke(elemento, 2, coordenadas)
    dofs = eledof[node_id].flatten()

    K[np.ix_(dofs, dofs)] += Ke

F = np.zeros(n_gdl)

bordes_derechos = [[1, 8, 2]]
bordes_izquierdos = [[0, 12, 3]]

for bordes, signo in [
    (bordes_derechos, 1),
    (bordes_izquierdos, -1)
]:
    for a, medio, b in bordes:
        longitud = np.linalg.norm(nodos[b] - nodos[a])
        carga = signo * sigma_x * espesor * longitud

        F[2*a] += carga / 6
        F[2*medio] += 4*carga / 6
        F[2*b] += carga / 6

restringidos = [0, 6, 7]

libres = [
    i for i in range(n_gdl)
    if i not in restringidos
]

Kr = K[np.ix_(libres, libres)]

D = np.zeros(n_gdl)
D[libres] = np.linalg.solve(Kr, F[libres])

print("\nDesplazamientos:")
print(D)

D_teorico = np.zeros(n_gdl)

D_teorico[0::2] = sigma_x / E * nodos[:, 0]
D_teorico[1::2] = -nu * sigma_x / E * (nodos[:, 1] - 100)

pasa_desplazamientos = np.allclose(
    D, D_teorico, atol=1e-10, rtol=1e-6
)

sigma = np.zeros((3, n_nodos))
contador = np.zeros(n_nodos)

sigma_teorica = np.array([sigma_x, 0, 0])

pasa_patch = True

for elemento in elementos:
    node_id = elemento.node
    coordenadas = nodos[node_id]
    dofs = eledof[node_id].flatten()
    De = D[dofs]

    sigma_gauss, q = gauss_stress(
        elemento, 2, coordenadas, De
    )

    print("\nElemento:", elemento.id)
    print("Tensiones en Gauss:")
    print(sigma_gauss)

    pasa_elemento = np.allclose(
        sigma_gauss,
        sigma_teorica[:, None],
        atol=1e-10,
        rtol=1e-8
    )

    if pasa_elemento:
        print("PASA")
    else:
        print("NO PASA")
        pasa_patch = False

    sigma_elemento = node_stress(
        elemento, coordenadas, De
    )

    for j, nodo in enumerate(node_id):
        sigma[:, nodo] += sigma_elemento[:, j]
        contador[nodo] += 1

sigma /= contador[np.newaxis, :]

print("\nTensiones nodales promediadas:")
print(sigma)

pasa_nodos = np.allclose(
    sigma,
    sigma_teorica[:, None],
    atol=1e-10,
    rtol=1e-8
)

residuo = K @ D - F

equilibrio = np.allclose(
    residuo[libres], 0, atol=1e-9
)

print("\nDesplazamientos correctos:", pasa_desplazamientos)
print("Tensiones nodales correctas:", pasa_nodos)
print("Equilibrio:", equilibrio)
print("Reacciones:", residuo[restringidos])

if (pasa_patch and pasa_nodos
        and pasa_desplazamientos and equilibrio):
    print("\nPATCH TEST T6: APROBADO")
else:
    print("\nPATCH TEST T6: DESAPROBADO")
