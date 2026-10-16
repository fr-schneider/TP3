import math
import numpy as np
from source.solvers.gauss_integration import gauss_1d

def apply_dirichlet_bc(K, f, dof, value):
    for i in range(K.number_of_dofs):
        f[i] = f[i] - K.values[i][dof] * value

    for i in range(K.number_of_dofs):
        K.values[dof][i] = 0.0
        K.values[i][dof] = 0.0

    K.values[dof][dof] = 1.0
    f[dof] = value




def robin_soil_stiffness_2d(element, node_coord, edge_local_nodes, k_soil, orden=None):
    """
    K_robin = ∫ N^T K_S N dΓ  (por unidad de espesor)
    edge_local_nodes: [v1, v2] o [v1, v2, medio]  (misma convención que surface_load_2d)
    k_soil: escalar (kx = ky) o (kx, ky) en fuerza/longitud^3, ejes globales
    """
    n_edge = len(edge_local_nodes)
    if n_edge not in (2, 3):
        raise ValueError("El borde debe tener 2 o 3 nodos")

    if np.ndim(k_soil) == 0:
        kx = ky = float(k_soil)
    else:
        kx, ky = float(k_soil[0]), float(k_soil[1])
    K_S = np.diag([kx, ky])

    node_coord = np.asarray(node_coord, dtype=float)
    ndofs = element.ndofs
    dofs_node = element.dofs_node
    if orden is None:
        orden = 2 if n_edge == 2 else 3     # N^T N: grado 2 (lineal) o 4 (cuadratico)

    # orden de las N 1D: [extremo xi=-1, extremo xi=+1, medio]
    def N_1d(xi):
        if n_edge == 2:
            return np.array([0.5 * (1 - xi), 0.5 * (1 + xi)])
        return np.array([0.5 * xi * (xi - 1), 0.5 * xi * (xi + 1), 1 - xi**2])

    def dN_1d(xi):
        if n_edge == 2:
            return np.array([-0.5, 0.5])
        return np.array([xi - 0.5, xi + 0.5, -2 * xi])

    edge_coords = node_coord[edge_local_nodes]
    gauss_pts, weights = gauss_1d(orden)
    K_robin = np.zeros((ndofs, ndofs))

    for xi, w in zip(gauss_pts, weights):
        N_val, dN_val = N_1d(xi), dN_1d(xi)
        J = math.hypot(dN_val @ edge_coords[:, 0], dN_val @ edge_coords[:, 1])

        N_mat = np.zeros((dofs_node, ndofs))
        for k, loc in enumerate(edge_local_nodes):
            for d in range(dofs_node):
                N_mat[d, loc * dofs_node + d] = N_val[k]

        K_robin += w * J * (N_mat.T @ K_S @ N_mat)
    return K_robin