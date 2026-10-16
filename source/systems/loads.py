
import numpy as np
from source.solvers.gauss_integration import gauss, gauss_1d
from source.solvers.solvers import jacobiano


def surface_load_2d(element, node_coord, load_node_id, traction_nodes, espesor=1.0):
    """
    Calcula las fuerzas nodales equivalentes de una carga superficial.

    Fe = integral_borde N.T @ traccion * espesor * dS
    """

    n_ndofs = element.ndofs
    n_nodos = element.n_node

    traction_nodes = np.asarray(traction_nodes, dtype=float)

    n_surface_node = len(load_node_id)

    if n_surface_node == 2:
        orden = 2
    elif n_surface_node == 3:
        orden = 3
    else:
        raise ValueError("El borde debe tener 2 o 3 nodos")

    r1, s1 = element.node_nat_coord[load_node_id[0]]
    r2, s2 = element.node_nat_coord[load_node_id[1]]

    dr_dt = (r2 - r1) / 2
    ds_dt = (s2 - s1) / 2

    q, w = gauss_1d(orden)

    fe = np.zeros(n_ndofs)

    for i, t in enumerate(q):

        r = r1 + (t + 1) / 2 * (r2 - r1)
        s = s1 + (t + 1) / 2 * (s2 - s1)

        j, det_j, inv_j = jacobiano(
            element, r, s, node_coord
        )

        tang = np.array([dr_dt, ds_dt]) @ j
        J = np.hypot(tang[0], tang[1])

        N = np.array([
            element.shape(r, s, k)
            for k in range(1, n_nodos + 1)
        ])

        traction_vec = N[load_node_id] @ traction_nodes

        for k in range(n_nodos):
            fe[2*k] += (
                N[k] * traction_vec[0] * J * w[i] * espesor
            )

            fe[2*k+1] += (
                N[k] * traction_vec[1] * J * w[i] * espesor
            )

    return fe


def body_load_2d(element, node_coord, b, espesor=1.0, orden=2):
    """
    Calcula las fuerzas nodales equivalentes de una carga de volumen.

    Fe = integral_area N.T @ b * espesor * dA

    b = [bx, by] representa la fuerza por unidad de volumen.
    Para el peso propio: b = [0, -rho*g].
    """

    q, w = gauss(element, orden)

    b = np.asarray(b, dtype=float)

    fe = np.zeros(element.ndofs)

    for (r, s), wi in zip(q, w):

        j, det_j, inv_j = jacobiano(
            element, r, s, node_coord
        )

        if det_j <= 0:
            raise ValueError("Jacobiano no positivo")

        N = np.array([
            element.shape(r, s, k)
            for k in range(1, element.n_node + 1)
        ])

        for k in range(element.n_node):

            fe[2*k] += (
                N[k] * b[0] * det_j * wi * espesor
            )

            fe[2*k+1] += (
                N[k] * b[1] * det_j * wi * espesor
            )

    return fe
