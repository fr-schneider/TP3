
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from source.elements.elements import quad_q9
from source.materials.linear_elastic import LinearElastic
from source.solvers.solvers import ke, node_stress, gauss_stress

coordnodos = np.array([
    [0, 0],
    [25, 0],
    [50, 0],
    [75, 0],
    [100, 0],
    [100, 25],
    [100, 50],
    [100, 75],
    [100, 100],
    [75, 100],
    [50, 100],
    [25, 100],
    [0, 100],
    [0, 75],
    [0, 50],
    [0, 25],
    [47.5, 27.5],
    [45, 55],
    [22.5, 52.5],
    [72.5, 52.5],
    [47.5, 77.5],
    [23.75, 26.25],
    [73.75, 26.25],
    [73.75, 76.25],
    [23.75, 76.25]
], dtype=float)

nodos_q9 = [
    [0, 2, 17, 14, 1, 16, 18, 15, 21],
    [2, 4, 6, 17, 3, 5, 19, 16, 22],
    [17, 6, 8, 10, 19, 7, 9, 20, 23],
    [14, 17, 10, 12, 18, 20, 11, 13, 24]
]

E = 21000
nu = 0.25
espesor = 1.0
sigma_x = 0.02

material = LinearElastic(E, nu)
tipo = "Plane_stress"

elementos = [
    quad_q9(i + 1, nodos_q9[i], material, tipo)
    for i in range(4)
]

n_nodos = len(coordnodos)
n_gdl = 2 * n_nodos

eledof = np.arange(n_gdl).reshape(n_nodos, 2)

K = np.zeros((n_gdl, n_gdl))

for elemento in elementos:
    node_id = elemento.node
    coordenadas = coordnodos[node_id]

    Ke = espesor * ke(elemento, 3, coordenadas)
    dofs = eledof[node_id].flatten()

    for i in range(elemento.ndofs):
        for j in range(elemento.ndofs):
            K[dofs[i], dofs[j]] += Ke[i, j]

F = np.zeros(n_gdl)

bordes_derechos = [
    [4, 5, 6],
    [6, 7, 8]
]

bordes_izquierdos = [
    [0, 15, 14],
    [14, 13, 12]
]

for bordes, signo in [
    (bordes_derechos, 1),
    (bordes_izquierdos, -1)
]:
    for a, medio, b in bordes:
        longitud = np.linalg.norm(
            coordnodos[b] - coordnodos[a]
        )

        carga = signo * sigma_x * espesor * longitud

        F[2*a] += carga / 6
        F[2*medio] += 4*carga / 6
        F[2*b] += carga / 6

print("Fuerzas nodales:")
print(F)

restringidos = [0, 1, 24]

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

D_teorico[0::2] = sigma_x / E * coordnodos[:, 0]
D_teorico[1::2] = -nu * sigma_x / E * coordnodos[:, 1]

pasa_desplazamientos = np.allclose(
    D,
    D_teorico,
    atol=1e-10,
    rtol=1e-6
)

print("\nDesplazamientos teoricos correctos:",
      pasa_desplazamientos)

sigma = np.zeros((3, n_nodos))
contador = np.zeros(n_nodos)

sigma_teorica = np.array([sigma_x, 0, 0])

pasa_patch = True

for elemento in elementos:
    node_id = elemento.node
    coordenadas = coordnodos[node_id]
    dofs = eledof[node_id].flatten()
    De = D[dofs]

    sigma_gauss, q = gauss_stress(
        elemento, 3, coordenadas, De
    )

    print("\nElemento:", elemento.id)
    print("Tensiones en Gauss:")
    print(sigma_gauss)

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
    residuo[libres],
    0,
    atol=1e-9
)

reacciones = residuo[restringidos]
reacciones_ok = np.allclose(reacciones, 0, atol=1e-9)

print("\nReacciones en GDL restringidos:")
print(reacciones)

print("\nDesplazamientos correctos:", pasa_desplazamientos)
print("Tensiones nodales correctas:", pasa_nodos)
print("Equilibrio:", equilibrio)
print("Reacciones nulas:", reacciones_ok)

if (pasa_patch and pasa_nodos and pasa_desplazamientos
        and equilibrio and reacciones_ok):
    print("\nPATCH TEST Q9: APROBADO")
else:
    print("\nPATCH TEST Q9: DESAPROBADO")
