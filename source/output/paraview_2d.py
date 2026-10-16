import os
import numpy as np
import meshio


def build_cells_2d(mesh):
    cell_type_map = {
        2: "line", 3: "triangle", 4: "quad",
        6: "triangle6", 8: "quad8", 9: "quad9"
    }
    cells_dict = {}
    for conn in mesh.connectivity:
        n_nodes = len(conn)
        c_type = cell_type_map.get(n_nodes, "polygon")
        if c_type not in cells_dict:
            cells_dict[c_type] = []
        cells_dict[c_type].append(conn)

    return [(c_type, np.array(conns)) for c_type, conns in cells_dict.items()]


def write_results_2d(filename, mesh, u, nodal_stresses=None):
    """Exporta la solución 2D a ParaView mediante meshio."""
    directory = os.path.dirname(filename)
    if directory != "":
        os.makedirs(directory, exist_ok=True)

    points = []
    displacement = []

    for node in mesh.nodes:
        dofs = mesh.node_global_dofs(node.global_id)
        points.append([node.x, node.y, 0.0])
        # Accedemos a los índices 0 y 1 de dofs
        displacement.append([u[dofs[0]], u[dofs[1]], 0.0])

    point_data = {"displacement": np.array(displacement)}

    # Integración con compute_nodal_stresses
    if nodal_stresses is not None:
        if isinstance(nodal_stresses, dict):
            # Agrega directamente: sigma_x, sigma_y, sigma_z, tau_xy, von_mises
            point_data.update(nodal_stresses)
        else:
            stresses = np.asarray(nodal_stresses)
            if stresses.ndim == 2 and stresses.shape[0] != 3:
                stresses = stresses.T

            names = ["sigma_x", "sigma_y", "tau_xy", "sigma_z", "von_mises"]
            for i in range(min(len(stresses), len(names))):
                point_data[names[i]] = stresses[i, :]

    cells = build_cells_2d(mesh)
    output_mesh = meshio.Mesh(
        points=np.array(points), cells=cells, point_data=point_data
    )

    output_mesh.write(filename)
    print(f"Resultados exportados para ParaView en: {filename}")