import numpy as np
from source.solvers.gauss_integration import gauss_1d
from source.solvers.solvers import jacobiano


def surface_load_2d(element, node_coord, load_node_id, traction_nodes):
    """
    Cargas nodales equivalentes de una traccion superficial sobre un borde
    del elemento:

        {f_e} = integral_borde  [N]^T {Phi}  * dS

    donde [N] son las funciones de forma del elemento 2D y {Phi} = (tx, ty)
    es la traccion (fuerza por unidad de superficie) en ejes globales.

    Parametros
    ----------
    element        : elemento
    node_coord     : coordenadas de los nodos del ELEMENTO (n_nodos, 2), igual que en ke
    load_node_id   : ids LOCALES (0..n_nodos-1) de los nodos del borde cargado.
                     Los dos primeros son los vertices del borde; si hay un
                     tercero, es el nodo medio.
    traction_nodes : (m, 2) traccion (tx, ty) evaluada en cada nodo del borde,
                     en el mismo orden que load_node_id


    Devuelve
    --------
    fe : vector de n_ndofs con las fuerzas nodales [u1, v1, u2, v2, ...]
    """

    n_ndofs = element.ndofs
    n_nodos = element.n_node

    traction_nodes = np.asarray(traction_nodes, dtype=float)    # (m, 2)


    # Tomo los dos vertices del borde; el borde es una recta en (r, s).

    r1, s1 = element.node_nat_coord[load_node_id[0]]
    r2, s2 = element.node_nat_coord[load_node_id[1]]

    # Parametrizo el borde con t en [-1, 1]:
    #   t = -1  ->  vertice 1,   t = +1  ->  vertice 2
    #   r(t) = r1 + (t+1)/2 * (r2 - r1),   s(t) = s1 + (t+1)/2 * (s2 - s1)
    # Entonces dr/dt y ds/dt son constantes sobre el borde.
    dr_dt = (r2 - r1) / 2
    ds_dt = (s2 - s1) / 2


    # Puntos de gauss
    # Borde de 2 nodos: traccion lineal -> 2 puntos.
    # Borde de 3 nodos: traccion hasta cuadratica -> 3 puntos.

    n_surface_node = len(load_node_id)
    if n_surface_node == 2:
        orden = 2
    elif n_surface_node == 3:
        orden = 3

    q, w = gauss_1d(orden)

    fe = np.zeros(n_ndofs)


    # integracion numerica  integral(...) dS = sum( ... * J * w )

    for i, t in enumerate(q):

        # punto de Gauss del borde expresado en (r, s) del elemento
        r = r1 + (t + 1) / 2 * (r2 - r1)
        s = s1 + (t + 1) / 2 * (s2 - s1)

        #  jacobiano del borde (dS = J dt).
        #         j = [[dx/dr, dy/dr],
        #              [dx/ds, dy/ds]]
        #     Regla de la cadena sobre el borde:
        #         [dx/dt, dy/dt] = [dr/dt, ds/dt] @ j
        #     J es el largo de ese vector tangente (para un borde recto = L/2).
        j, det_j, inv_j = jacobiano(element, r, s, node_coord)
        tang = np.array([dr_dt, ds_dt]) @ j
        J = np.hypot(tang[0], tang[1])

        # funciones de forma del elemento 2D evaluadas en (r, s).
        #     Sobre el borde, las N de los nodos que no estan en el borde valen 0.
        N = np.array([element.shape(r, s, k) for k in range(1, n_nodos + 1)])

        # traccion en este punto de Gauss: interpolo los valores nodales
        #     del borde con las mismas funciones de forma:
        #         Phi(t) = sum_k N_k * Phi_k     (k = nodos del borde)
        traction_vec = N[load_node_id] @ traction_nodes          #

        # acumulo  N_k * Phi * J * w  en los gdl del nodo k
        for k in range(n_nodos):
            fe[2 * k]     += N[k] * traction_vec[0] * J * w[i]
            fe[2 * k + 1] += N[k] * traction_vec[1] * J * w[i]


