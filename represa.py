from pathlib import Path

import numpy as np

from source.mesh.mesh import Mesh
from source.materials.linear_elastic import LinearElastic
from source.output.paraview_2d import write_results_2d
from source.systems.boundary_conditions import robin_soil_stiffness_2d
from source.systems.loads import surface_load_2d 
#Script principal de analisisd de la represa

malla='represa.txt'
output='represa.vtu'
tipo='aca va el plane_strain'
orden= 2
g= 9.81
h_agua= 25
rho_agua= 1000
rho_hormigon= 2300
k_suelo= [200,500] #k_horizontal, k_vertical



## CONFIGURACIÓN DE RUTAS Y CARPETAS

# Ubicación actual del script (C:\Users\schne\Desktop\tp3\Represa.py)
SCRIPT_DIR = Path(__file__).resolve().parent

# La raíz del repositorio ES directamente SCRIPT_DIR
REPO_ROOT = SCRIPT_DIR

# Directorios de entrada y salida dentro de tp3
MESHES_DIR = REPO_ROOT / "meshes"
RESULTS_DIR = REPO_ROOT / "results"

# Definición de archivos
input_file = MESHES_DIR / malla       # Buscará en tp3/meshes/represa.txt
output_file = RESULTS_DIR / output     # Escribirá en tp3/results/solucion_2d_represa.vtu

input_path_str = str(input_file)
output_path_str = str(output_file)

# Convertimos a string por compatibilidad con métodos de lectura/escritura
input_path_str = str(input_file)
output_path_str = str(output_file)

## Lectura de malla y creacion de elementos
mesh = Mesh()
mesh.read(input_file)
mesh.dofs_per_node = 2  


elements = mesh.create_elements(tipo)

n_dofs = mesh.number_of_dofs()
coord_global = np.array([[node.x, node.y] for node in mesh.nodes])

#Creacion de Kglobal e incorporacion de Krobin
K_global = mesh.assemble_stiffness(elements,2)

#Condiciones de Robin Sobre nodos y=0
y_min = min(node.y for node in mesh.nodes)
for elemento in elements:
    coord_node= coord_global[elemento.node]
    base_local_nodes=[
        loc_idx for loc_idx, g_id in enumerate(elemento.node)
        if np.isclose(mesh.node(g_id).y, y_min)
    ]

    #Si tiene 2 nodos en la base
        
    if len(base_local_nodes) >= 2:
        base_sorted = sorted(base_local_nodes, key=lambda i: mesh.node(elem.node[i]).x)
        
        # Calculamos la matriz K_robin elemental para el suelo
        Ke_robin = robin_soil_stiffness_2d(
            elemento, coord_node, base_sorted, k_suelo, orden
        )
        
        # Ensamblado directo en K_global
        for i_loc, g_node_i in enumerate(elemento.node):
            dofs_i = mesh.node_global_dofs(g_node_i)
            for j_loc, g_node_j in enumerate(elemento.node):
                dofs_j = mesh.node_global_dofs(g_node_j)
                K_global[np.ix_(dofs_i, dofs_j)] += Ke_robin[
                    2 * i_loc : 2 * i_loc + 2, 2 * j_loc : 2 * j_loc + 2
                ]


##Cargas
F_global = np.zeros(n_dofs)

#peso propio
for elemento in elements:
    nodes = elemento.node
    coords = coord_global[nodes]
    
    
    x, y = coords[:, 0], coords[:, 1]
    area = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
    
    # Peso asignado a cada nodo del elemento (fuerza hacia abajo en Y)
    peso_nodo = (rho_hormigon * g * area) / len(nodes)
    
    for g_id in nodes:
        dof_y = mesh.node_global_dofs(g_id)[1]  # Grado de libertad en Y
        F_global[dof_y] -= peso_nodo


#Presion de agua
x_min = min(node.x for node in mesh.nodes)

for elem in elements:
    coord_node = coord_global[elem.node]
    
    # Nodos del elemento sobre el paramento de aguas arriba expuesto al agua
    water_local_nodes = [
        loc_idx for loc_idx, g_id in enumerate(elem.node)
        if np.isclose(mesh.node(g_id).x, x_min) and mesh.node(g_id).y <= h_agua
    ]
    
    if len(water_local_nodes) >= 2:
        water_sorted = sorted(water_local_nodes, key=lambda i: mesh.node(elem.node[i]).y)
        
        # Presión hidrostática p(y) = rho * g * (h_agua - y) empujando horizontalmente (+x)
        traction_nodes = []
        for loc_idx in water_sorted:
            y_node = mesh.node(elem.node[loc_idx]).y
            p_hydro = rho_agua * g * (h_agua - y_node)
            traction_nodes.append([p_hydro, 0.0])
            
        traction_nodes = np.array(traction_nodes)
        
        # Integración de la carga de superficie en el borde
        fe_water = surface_load_2d(elem, coord_node, water_sorted, traction_nodes, orden=orden)
        
        if fe_water is not None:
            fe_water = fe_water.flatten()
            for i_loc, g_id in enumerate(elem.node):
                dofs = mesh.node_global_dofs(g_id)
                F_global[dofs] += fe_water[2*i_loc : 2*i_loc + 2]


##Resolucion k*u=f
u_global = np.linalg.solve(K_global, F_global)

##Exportacion a paraview

write_results_2d(output_path_str, mesh, u_global)