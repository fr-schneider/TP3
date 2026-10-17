
from pathlib import Path

import numpy as np

from source.mesh.mesh import Mesh
from source.materials.linear_elastic import LinearElastic
from source.output.paraview_2d import write_results_2d
from source.systems.boundary_conditions import robin_soil_stiffness_2d
from source.systems.loads import surface_load_2d,body_load_2d

# Script principal de analisis de la represa

malla = "represa_Q4_gruesa.txt"
output = "represa_Q4_gruesa.vtu" 
orden = 2
g = 9.81
h_agua = 27
rho_agua = 1000
rho_hormigon = 2300
k_suelo = [2e8, 2e9]  # k_horizontal, k_vertical
tipo = 'Plane_strain'
Robin_BC=False
orden=2
tol=10e-3


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


# Creacion de Kglobal

K_global = mesh.assemble_stiffness(elements, 2)

## Cargas

F_global = np.zeros(n_dofs)


# Peso propio
for elem in elements:
 
    coord_node = coord_global[elem.node]
 
    fe_peso = body_load_2d(
        elem, coord_node, [0.0, -rho_hormigon * g], 1, orden
    )
 
    for i_loc, g_id in enumerate(elem.node):
        F_global[mesh.node_global_dofs(g_id)] += fe_peso[2*i_loc:2*i_loc+2]


# Presion de agua
x_min = min(node.x for node in mesh.nodes)
y_min = min(node.y for node in mesh.nodes)
y_agua = y_min + h_agua
 
for elem in elements:
 
    coord_node = coord_global[elem.node]
 
    nv = 3 if elem.n_node in (3, 6) else 4
 
    for k in range(nv):
 
        # borde k: [vertice, vertice, medio] (el medio solo en elementos de orden 2)
        water_nodes = [k, (k + 1) % nv]
        if elem.n_node > nv:
            water_nodes.append(nv + k)
 
        y_borde = coord_node[water_nodes, 1]
 
        if all(abs(coord_node[water_nodes, 0] - x_min) < tol) and all(y_borde <= y_agua + tol):
 
            p_hydro = rho_agua * g * np.maximum(y_agua - y_borde, 0.0)
            traction_nodes = np.column_stack((p_hydro, np.zeros_like(p_hydro)))
 
            fe_water = surface_load_2d(
                elem, coord_node, water_nodes, traction_nodes, espesor=1
            )
 
            for i_loc, g_id in enumerate(elem.node):
                F_global[mesh.node_global_dofs(g_id)] += fe_water[2*i_loc:2*i_loc+2]


K_solido = K_global.copy()   # rigidez del sólido, sin Robin ni Dirichlet
F_cargas = F_global.copy()   # cargas externas (peso + agua), sin modificar


#Condiciones de borde

if Robin_BC:
    Krobin_Global = np.zeros_like(K_global)

    # Condiciones de Robin sobre nodos y = 0
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
                elemento, coord_node, base_sorted, k_suelo, orden
            )

            for i_loc, g_node_i in enumerate(elemento.node):

                dofs_i = mesh.node_global_dofs(g_node_i)

                for j_loc, g_node_j in enumerate(elemento.node):

                    dofs_j = mesh.node_global_dofs(g_node_j)

                    Krobin_Global[np.ix_(dofs_i, dofs_j)] += Ke_robin[
                        2*i_loc:2*i_loc+2,
                        2*j_loc:2*j_loc+2
                    ]


        # Sumar Robin

    K_global += Krobin_Global

else: #aplicar v=0 en todo x, y min, y u=0 en x_max, y min
    base_nodes = [node for node in mesh.nodes if np.isclose(node.y, y_min)]

        # Seleccionamos el nodo en x_máx sobre la base para restringir u
    node_xmax_ymin = max(base_nodes, key=lambda node: node.x)

    fixed_dofs = []

        #Fijar v = 0 para todos los nodos de la base
    for node in base_nodes:
        dofs = mesh.node_global_dofs(node.global_id)
        fixed_dofs.append(dofs[1])  # dof de v (dirección y)

        #  Fijar u = 0 en x_max, y_min para evitar traslación rígida horizontal
    dofs_pin = mesh.node_global_dofs(node_xmax_ymin.global_id)
    fixed_dofs.append(dofs_pin[0])  # dof de u (dirección x)

        # Aplicación de restricciones fijando 1.0 en la diagonal de K_global tipo apply_dirichlet_bc
    for dof in set(fixed_dofs):
        K_global[dof, :] = 0.0
        K_global[:, dof] = 0.0
        K_global[dof, dof] = 1.0
        F_global[dof] = 0.0



## Cargas




## Resolucion K*D=F


D = np.linalg.solve(K_global, F_global)


## Calculo tensiones

sigma_nodal = mesh.compute_nodal_stresses(elements, D, tipo)


## Exportacion a ParaView

write_results_2d(output_path_str, mesh, D, sigma_nodal)



## Estabilidad al vuelco: calculo teorico con la geometria de la malla
# Seccion: trapecio con paramento de aguas arriba vertical (x = x_min), corona de ancho
# "corona" y talud de aguas abajo recto hasta el pie (x_max, y_min)

x_pie = coord_global[:, 0].max()
y_top = coord_global[:, 1].max()

altura = y_top - y_min
base = x_pie - x_min
corona = coord_global[np.abs(coord_global[:, 1] - y_top) < tol, 0].max() - x_min
h_mojada = min(h_agua, altura)

# Area y centro de gravedad (x medido desde el paramento de aguas arriba)

A_rect = corona * altura
A_tri = (base - corona) * altura / 2
A_total = A_rect + A_tri

x_G = (A_rect * corona / 2 + A_tri * (corona + (base - corona) / 3)) / A_total

# Peso y empuje por metro longitudinal

W = rho_hormigon * g * A_total
F_agua = rho_agua * g * h_mojada**2 / 2

# Momentos respecto del pie

M_estabilizante = W * (base - x_G)
M_desestabilizante = F_agua * h_mojada / 3
FS = M_estabilizante / M_desestabilizante

print("\nRESULTADOS DE ESTABILIDAD")

print("Momento estabilizante:", M_estabilizante, "N.m/m")
print("Momento desestabilizante:", M_desestabilizante, "N.m/m")
print("Factor de seguridad:", FS)


# VERIFICACION CON REACCIONES DEL SUELO
# Reaccion = fuerza del suelo sobre la represa = K*D - F (sirve para Robin y para apoyo fijo)

R = K_solido @ D - F_cargas

base_idx = np.abs(coord_global[:, 1] - y_min) < tol

M_reacciones = np.sum((coord_global[base_idx, 0] - x_pie) * R[1::2][base_idx])

print("Suma de Ry:", R[1::2].sum(), "(peso:", W, ")")
print("Suma de Rx:", R[0::2].sum(), "(empuje:", -F_agua, ")")
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
