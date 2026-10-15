import os
import numpy as np
import meshio


def build_cells_2d(mesh):
    """
    Agrupa la conectividad de la malla según el tipo de elemento 2D
    usando los nombres de celda VTK estándar de meshio.
    """
    cell_type_map = {
        2: "line",
        3: "triangle",
        4: "quad",
        6: "triangle6",
        8: "quad8",
        9: "quad9"
    }
    cells_dict = {}
    for conn in mesh.connectivity:
        n_nodes = len(conn)
        c_type = cell_type_map.get(n_nodes, "polygon")
        if c_type not in cells_dict:
            cells_dict[c_type] = []
        cells_dict[c_type].append(conn)

    return [(c_type, np.array(conns)) for c_type, conns in cells_dict.items()]





def write(filename, points, cells, displacement, cell_data_dict):
    directory = os.path.dirname(filename)
    if directory != "":
        os.makedirs(directory, exist_ok=True)

    output_mesh = meshio.Mesh(
        points=points,
        cells=cells,
        point_data={"displacement": displacement},
        cell_data=cell_data_dict,
    )

    output_mesh.write(filename)


def write_results_2d(filename, mesh, elements, u, state_type='Plane_stress', nu=0.3):
    """
    Exporta la solución 2D (desplazamientos y campos de tensión sx, sy, sxy, sz) a ParaView.
    """
    points = []
    displacement = []

    # 1. Nodos y desplazamientos en formato 3D para VTK
    for node in mesh.nodes:
        dofs = mesh.node_global_dofs(node.global_id)
        points.append([node.x, node.y, 0.0])
        displacement.append([u[dofs], u[dofs[4]], 0.0])

    # 2. Matriz de coordenadas globales
    coord_global = np.array([[n.x, n.y] for n in mesh.nodes])
    
    sx_list, sy_list, txy_list, sz_list = [], [], [], []

    for element in elements:
        coord_node = coord_global[element.node]
        sx, sy, txy, sz = compute_element_stresses(element, u, coord_node, state_type, nu)
        sx_list.append(sx)
        sy_list.append(sy)
        txy_list.append(txy)
        sz_list.append(sz)

    # 3. Formateo de celdas y diccionario de tensiones
    cells = build_cells_2d(mesh)
    
    cell_data = {
        "sigma_x": [np.array(sx_list)],
        "sigma_y": [np.array(sy_list)],
        "tau_xy":  [np.array(txy_list)],
        "sigma_z": [np.array(sz_list)],
    }

    write(filename, points, cells, displacement, cell_data)