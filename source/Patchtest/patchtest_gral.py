import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from source.mesh.mesh import Mesh
from source.materials.linear_elastic import LinearElastic
from source.systems.loads import surface_load_2d
from source.output.paraview_2d import write_results_2d


def main():
    # --------------------------------------------------------------------------
    # 1. CONFIGURACIÓN DE RUTAS Y CARPETAS
    # --------------------------------------------------------------------------
    # Ubicación actual del script: tp3/source/patchtest/main.py
    SCRIPT_DIR = Path(__file__).resolve().parent
    
    # Raíz del repositorio: subimos 2 niveles (patchtest -> source -> tp3)
    REPO_ROOT = SCRIPT_DIR.parents[1]

    # Directorios de entrada y salida
    MESHES_DIR = REPO_ROOT / "meshes"
    RESULTS_DIR = REPO_ROOT / "results"



    # Definición de archivos
    input_file = MESHES_DIR / "malla_q4_patch.txt"      # Lee desde tp3/meshes/
    output_file = RESULTS_DIR / "solucion_2d_q4_patch.vtu"    # Escribe en tp3/results/

    # Convertimos a string por compatibilidad con métodos de lectura/escritura
    input_path_str = str(input_file)
    output_path_str = str(output_file)

    E = 210000.0     # MPa
    nu = 0.3
    sigma_0 = 100.0  # MPa

    # 1. Lectura de malla y creacion de elementos
    mesh = Mesh()
    mesh.read(input_file)
    mesh.dofs_per_node = 2  

    material = LinearElastic(E=E, nu=nu)
    elements = mesh.create_elements(material, tipo="Plane_stress")

    n_dofs = mesh.number_of_dofs()
    coord_global = np.array([[node.x, node.y] for node in mesh.nodes])
    # 2. Ensamblado de K_global a traves de Mesh
    K_global = mesh.assemble_stiffness(elements,2)
    F_global = np.zeros(n_dofs)

        # 4. Carga de superficie en el borde derecho (X = X_max)
    x_max = max(node.x for node in mesh.nodes)
    
    for elem in elements:
        coord_node = coord_global[elem.node]
        load_local_nodes = [
            loc_idx for loc_idx, g_id in enumerate(elem.node)
            if np.isclose(mesh.node(g_id).x, x_max)
        ]
        
        if len(load_local_nodes) >= 2:
            load_sorted = sorted(load_local_nodes, key=lambda i: mesh.node(elem.node[i]).y)
            traction_nodes = np.tile([sigma_0, 0.0], (len(load_sorted), 1))

            # Calculamos la carga de superficie y la aplanamos a vector 1D
            fe_elem = surface_load_2d(elem, coord_node, load_sorted, traction_nodes).flatten()

            # Ensamblado seguro usando los dofs globales del nodo
            for i_loc, g_id in enumerate(elem.node):
                dofs = mesh.node_global_dofs(g_id)
                F_global[dofs] += fe_elem[2 * i_loc : 2 * i_loc + 2]

    # 4. Apoyos automaticos en X = 0 (fijacion u_x = 0, y apoyo vertical u_y = 0)
    fixed_dofs = []
    nodes_x0 = [node for node in mesh.nodes if np.isclose(node.x, 0.0)]
    for node in nodes_x0:
        fixed_dofs.append(2 * node.global_id)  # u_x = 0
    
    node_y0 = min(nodes_x0, key=lambda n: abs(n.y))
    fixed_dofs.append(2 * node_y0.global_id + 1)  # u_y = 0 en el nodo base

    fixed_dofs = sorted(list(set(fixed_dofs)))
    free_dofs = [dof for dof in range(n_dofs) if dof not in fixed_dofs]

    # 5. Resolucion K_libre * u_libre = F_libre
    K_free = K_global[np.ix_(free_dofs, free_dofs)]
    F_free = F_global[free_dofs]
    u_free = np.linalg.solve(K_free, F_free)

    u_global = np.zeros(n_dofs)
    u_global[free_dofs] = u_free

    # 6. Promediado nodal de tensiones en la malla
    sigma_nodal = mesh.compute_nodal_stresses(elements, u_global)

    # 7. Exportacion directa a ParaView
    write_results_2d(output_file, mesh, u_global, nodal_stresses=sigma_nodal)


if __name__ == "__main__":
    main()