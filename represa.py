
from pathlib import Path

import numpy as np

from source.mesh.mesh import Mesh
from source.materials.linear_elastic import LinearElastic
from source.output.paraview_2d import write_results_2d
from source.systems.boundary_conditions import robin_soil_stiffness_2d
from source.systems.loads import surface_load_2d

# Script principal de analisis de la represa

malla = 'represa.txt'
output = 'represa.vtu'
orden = 2
g = 9.81
h_agua = 25
rho_agua = 1000
rho_hormigon = 2300
k_suelo = [200, 500]  # k_horizontal, k_vertical
tipo = 'Plane_strain'


## CONFIGURACIÓN DE RUTAS Y CARPETAS

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR

MESHES_DIR = REPO_ROOT / "meshes"
RESULTS_DIR = REPO_ROOT / "results"

input_file = MESHES_DIR / malla
output_file = RESULTS_DIR / output

input_path_str = str(input_file)
output_path_str = str(output_file)


## Lectura de malla y creacion de elementos

mesh = Mesh()
mesh.read(input_file)
mesh.dofs_per_node = 2

for mid in mesh.materials:
    material, _ = mesh.materials[mid]
    mesh.materials[mid] = (material, tipo)

elements = mesh.create_elements()

n_dofs = mesh.number_of_dofs()
coord_global = np.array([[node.x, node.y] for node in mesh.nodes])


# Creacion de Kglobal e incorporacion de Krobin

K_global = mesh.assemble_stiffness(elements, 2)
Krobin_Global = np.zeros_like(K_global)


# Condiciones de Robin sobre nodos y = 0

y_min = min(node.y for node in mesh.nodes)

for elemento in elements:

    coord_node = coord_global[elemento.node]

    base_local_nodes = [
        loc_idx for loc_idx, g_id in enumerate(elemento.node)
        if np.isclose(mesh.node(g_id).y, y_min)
    ]

    if len(base_local_nodes) >= 2:

        base_sorted = sorted(
            base_local_nodes,
            key=lambda i: mesh.node(elemento.node[i]).x
        )

        Ke_robin = robin_soil_stiffness_2d(
            elemento, coord_node, base_sorted, k_suelo, orden=orden
        )

        for i_loc, g_node_i in enumerate(elemento.node):

            dofs_i = mesh.node_global_dofs(g_node_i)

            for j_loc, g_node_j in enumerate(elemento.node):

                dofs_j = mesh.node_global_dofs(g_node_j)

                Krobin_Global[np.ix_(dofs_i, dofs_j)] += Ke_robin[
                    2*i_loc:2*i_loc+2,
                    2*j_loc:2*j_loc+2
                ]


# Sumar Robin una sola vez

K_global += Krobin_Global


## Cargas

F_global = np.zeros(n_dofs)


# Peso propio

for elemento in elements:

    nodes = elemento.node
    coords = coord_global[nodes]

    x, y = coords[:, 0], coords[:, 1]

    area = 0.5 * np.abs(
        np.dot(x, np.roll(y, 1)) -
        np.dot(y, np.roll(x, 1))
    )

    peso_nodo = (rho_hormigon * g * area) / len(nodes)

    for g_id in nodes:

        dof_y = mesh.node_global_dofs(g_id)[1]
        F_global[dof_y] -= peso_nodo


# Presion de agua

x_min = min(node.x for node in mesh.nodes)

for elem in elements:

    coord_node = coord_global[elem.node]

    water_local_nodes = [
        loc_idx for loc_idx, g_id in enumerate(elem.node)
        if np.isclose(mesh.node(g_id).x, x_min)
        and mesh.node(g_id).y <= h_agua
    ]

    if len(water_local_nodes) >= 2:

        water_sorted = sorted(
            water_local_nodes,
            key=lambda i: mesh.node(elem.node[i]).y
        )

        traction_nodes = []

        for loc_idx in water_sorted:

            y_node = mesh.node(elem.node[loc_idx]).y
            p_hydro = rho_agua * g * (h_agua - y_node)

            traction_nodes.append([p_hydro, 0.0])

        traction_nodes = np.array(traction_nodes)

        fe_water = surface_load_2d(
            elem, coord_node, water_sorted, traction_nodes
        )

        if fe_water is not None:

            fe_water = fe_water.flatten()

            for i_loc, g_id in enumerate(elem.node):

                dofs = mesh.node_global_dofs(g_id)

                F_global[dofs] += fe_water[
                    2*i_loc:2*i_loc+2
                ]


# Presion en los bordes parcialmente sumergidos
# BLOQUE MOVIDO: ahora esta antes de calcular D

for elem in elements:

    coord_node = coord_global[elem.node]

    for i, j in [(0, 1), (1, 2), (2, 0)]:

        if not (
            np.isclose(coord_node[i, 0], x_min)
            and np.isclose(coord_node[j, 0], x_min)
        ):
            continue

        yi = coord_node[i, 1]
        yj = coord_node[j, 1]

        if not (min(yi, yj) < h_agua < max(yi, yj)):
            continue

        if yi > yj:
            i, j = j, i
            yi, yj = yj, yi

        fraccion = (h_agua - yi) / (yj - yi)
        L = np.linalg.norm(coord_node[j] - coord_node[i])

        fe_water = np.zeros(elem.ndofs)

        puntos, pesos = np.polynomial.legendre.leggauss(2)

        for xi, wi in zip(puntos, pesos):

            s = fraccion * (xi + 1) / 2

            y = yi + s * (yj - yi)

            presion = rho_agua * g * (h_agua - y)

            dF = presion * L * fraccion * wi / 2

            fe_water[2*i] += (1-s) * dF
            fe_water[2*j] += s * dF

        for i_loc, g_id in enumerate(elem.node):

            dofs = mesh.node_global_dofs(g_id)

            F_global[dofs] += fe_water[
                2*i_loc:2*i_loc+2
            ]


## Resolucion K*D=F
# Ahora F_global incluye TODA la presion del agua

D = np.linalg.solve(K_global, F_global)


## Calculo tensiones

sigma_nodal = mesh.compute_nodal_stresses(elements, D, tipo)


## Exportacion a ParaView

write_results_2d(output_path_str, mesh, D, sigma_nodal)


# MOMENTOS TEORICOS DE LA REPRESA

H = 30.0
B = 18.0
C = 3.0

# Area de la seccion

A_rect = C * H
A_tri = (B - C) * H / 2
A_total = A_rect + A_tri

# Centro de gravedad

x_G = (
    A_rect * (C / 2)
    + A_tri * (C + (B - C) / 3)
) / A_total

# Peso de la represa por metro longitudinal

W = rho_hormigon * g * A_total

# Momento estabilizante respecto del pie

brazo_peso = B - x_G
M_estabilizante = W * brazo_peso

# Fuerza resultante del agua

F_agua = rho_agua * g * h_agua**2 / 2

# Altura de aplicacion de la fuerza del agua

y_agua = h_agua / 3

# Momento desestabilizante respecto del pie

M_desestabilizante = F_agua * y_agua

# Factor de seguridad al vuelco

FS = M_estabilizante / M_desestabilizante

print("\nRESULTADOS DE ESTABILIDAD")

print("Momento estabilizante:", M_estabilizante, "N.m/m")
print("Momento desestabilizante:", M_desestabilizante, "N.m/m")
print("Factor de seguridad:", FS)


# VERIFICACION CON REACCIONES ROBIN

R_robin = -Krobin_Global @ D

x_pie = max(
    node.x for node in mesh.nodes
    if np.isclose(node.y, y_min)
)

M_reacciones = 0.0

for nodo in mesh.nodes:

    if np.isclose(nodo.y, y_min):

        i = nodo.global_id

        brazo = nodo.x - x_pie

        Ry = R_robin[2*i + 1]

        M_reacciones += brazo * Ry

print("Momento de las reacciones:", M_reacciones)


# Verificar equilibrio de momentos

M_equilibrio = M_estabilizante - M_desestabilizante + M_reacciones

error_relativo = abs(M_equilibrio) / max(
    abs(M_estabilizante),
    abs(M_desestabilizante),
    1.0
)

print("Error relativo:", error_relativo * 100, "%")

if error_relativo < 0.01:
    print("Equilibrio de momentos verificado")
else:
    print("No se verifica el equilibrio de momentos")
