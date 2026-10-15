
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from source.elements.elements import quad_q8
from source.materials.linear_elastic import LinearElastic
from source.solvers.solvers import ke, node_stress, gauss_stress

coordnodos = np.array([
    [0, 0],        # 0
    [25, 0],       # 1
    [50, 0],       # 2
    [75, 0],       # 3
    [100, 0],      # 4
    [100, 25],     # 5
    [100, 50],     # 6
    [100, 75],     # 7
    [100, 100],    # 8
    [75, 100],     # 9
    [50, 100],     # 10
    [25, 100],     # 11
    [0, 100],      # 12
    [0, 75],       # 13
    [0, 50],       # 14
    [0, 25],       # 15
    [47.5, 27.5],  # 16
    [45, 55],      # 17
    [22.5, 52.5],  # 18
    [72.5, 52.5],  # 19
    [47.5, 77.5]   # 20
], dtype=float)
#primero los nodos de las esquinas luego los del medio
nodos_q8 = [
    [0, 2, 17, 14, 1, 16, 18, 15],
    [2, 4, 6, 17, 3, 5, 19, 16],
    [17, 6, 8, 10, 19, 7, 9, 20],
    [14, 17, 10, 12, 18, 20, 11, 13]
]

nodos_medios = {}

def nodo_medio(a, b):
    clave = tuple(sorted((a, b)))

    if clave not in nodos_medios:
        x = (coordnodos[a][0] + coordnodos[b][0]) / 2
        y = (coordnodos[a][1] + coordnodos[b][1]) / 2

        nodos_medios[clave] = len(coordnodos)
        coordnodos.append([x, y])

    return nodos_medios[clave]

coordnodos = np.array(coordnodos, dtype=float)


material = LinearElastic(21000, 0.25)
tipo = "Plane_stress"
espesor = 1.0


elementos = [
    quad_q8(i + 1, nodos_q8[i], material, tipo)
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



sigma_x = 0.02
F = np.zeros(n_gdl)


bordes = [
    [4, 5, 6],
    [6, 7, 8]
]

for a, medio, b in bordes:

    longitud = np.linalg.norm(
        coordnodos[b] - coordnodos[a]
    )

    carga = sigma_x * espesor * longitud

    F[2*a] += carga / 6
    F[2*medio] += 4*carga / 6
    F[2*b] += carga / 6

# Condiciones de borde
# ux = 0 sobre todo el borde izquierdo
# uy = 0 en el nodo 0 para impedir traslacion vertical
restringidos = [
    2*i for i, (x, y) in enumerate(coordnodos)
    if np.isclose(x, 0)
]
restringidos.append(1)

restringidos = sorted(set(restringidos))

libres = [
    i for i in range(n_gdl)
    if i not in restringidos
]


Kr = K[np.ix_(libres, libres)]
D = np.zeros(n_gdl)

D[libres] = np.linalg.solve(Kr, F[libres])

print("Desplazamientos:")
print(D)


sigma = np.zeros((3, n_nodos))
contador = np.zeros(n_nodos)


sigma_teorica = np.array([0.02, 0, 0])

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

# Verificar equilibrio
residuo = K @ D - F

equilibrio = np.allclose(
    residuo[libres],
    0,
    atol=1e-9
)

print("\nEquilibrio:", equilibrio)

# Resultado
if pasa_patch and pasa_nodos and equilibrio:
    print("\nPATCH TEST Q8: APROBADO")
else:
    print("\nPATCH TEST Q8: DESAPROBADO")
