import numpy as np

# Importaciones de tu librería
from mesh import Mesh
from source.materials.linear_elastic import LinearElastic
from source.systems.loads import surface_load_2d
from source.output.paraview_2d import write_results_2d


def main():
    # --------------------------------------------------------------------------
    # 1. PARÁMETROS DE ENTRADA Y ARCHIVOS
    # --------------------------------------------------------------------------
    input_file = "malla_input.txt"      # Archivo de malla (Q4, Q8, CST, LST, etc.)
    output_file = "solucion_2d.vtu"     # Archivo de salida para ParaView

    E = 210000.0     # Módulo de Elasticidad (MPa)
    nu = 0.3         # Coeficiente de Poisson
    thickness = 1.0  # Espesor (mm)
    sigma_0 = 100.0  # Tracción en X = X_max (MPa)

    # --------------------------------------------------------------------------
    # 2. LECTURA DE MALLA E INSTANCIACIÓN DE ELEMENTOS (GENÉRICO)
    # --------------------------------------------------------------------------
    print(f"--> Leyendo malla desde: {input_file}")
    mesh = Mesh(dofs_per_node=2)
    mesh.read(input_file)

    material = LinearElastic(E=E, nu=nu)
    elements = mesh.create_elements(material, tipo_elasticidad="Plane_stress")

    n_dofs = mesh.number_of_dofs()
    print(f"    Nodos: {mesh.number_of_nodes()} | Elementos: {mesh.number_of_elements()} | DOFs: {n_dofs}")

    # --------------------------------------------------------------------------
    # 3. ENSAMBLADO DE K_GLOBAL
    # --------------------------------------------------------------------------
    print("--> Ensamblando matriz de rigidez K_global...")
    K_global = mesh.assemble_stiffness(elements, thickness=thickness)
    F_global = np.zeros(n_dofs)

    # --------------------------------------------------------------------------
    # 4. CARGA DE SUPERFICIE AUTOMÁTICA EN X = X_max (e.g. X = 1.0 o X_max)
    # --------------------------------------------------------------------------
    # Identificamos la coordenada X máxima de la malla
    x_max = max(node.x for node in mesh.nodes)
    print(f"--> Aplicando tracción sigma_0 = {sigma_0} MPa en el borde X = {x_max:.2f}...")

    for elem in elements:
        # Buscamos los nodos locales del elemento que están sobre X = X_max
        load_local_nodes = []
        for local_idx, g_id in enumerate(elem.node):
            if np.isclose(mesh.node(g_id).x, x_max):
                load_local_nodes.append(local_idx)

        # Si el elemento tiene 2 o más nodos en el borde derecho, es una cara cargada
        if len(load_local_nodes) >= 2:
            # Ordenamos los nodos locales por su coordenada Y creciente para recorrer la cara en orden
            load_local_nodes_sorted = sorted(
                load_local_nodes,
                key=lambda idx: mesh.node(elem.node[idx]).y
            )

            # Calculamos la carga de superficie equivalente para este elemento
            fe_elem = surface_load_2d(
                elem,
                elem.coord_node,
                load_node_id=load_local_nodes_sorted,
                traction_vec=[sigma_0, 0.0],
                thickness=thickness
            )

            # Ensamblamos en el vector global F_global
            for i_loc, g_id in enumerate(elem.node):
                F_global[2 * g_id]     += fe_elem[2 * i_loc]
                F_global[2 * g_id + 1] += fe_elem[2 * i_loc + 1]

    # --------------------------------------------------------------------------
    # 5. RESTRICCIONES AUTOMÁTICAS EN X = 0 (FIJACIÓN)
    # --------------------------------------------------------------------------
    fixed_dofs = []

    # A) Bloqueamos u_x = 0 en TODOS los nodos sobre X = 0
    nodes_x0 = [node for node in mesh.nodes if np.isclose(node.x, 0.0)]
    for node in nodes_x0:
        fixed_dofs.append(2 * node.global_id)  # u_x = 0

    # B) Bloqueamos u_y = 0 en el nodo de X=0 más cercano a Y=0 (evita cuerpo rígido en Y)
    node_support = min(nodes_x0, key=lambda n: abs(n.y))
    fixed_dofs.append(2 * node_support.global_id + 1)  # u_y = 0

    fixed_dofs = sorted(list(set(fixed_dofs)))
    free_dofs = [dof for dof in range(n_dofs) if dof not in fixed_dofs]

    # --------------------------------------------------------------------------
    # 6. RESOLUCIÓN DEL SISTEMA REDUCIDO
    # --------------------------------------------------------------------------
    print("--> Resolviendo sistema de ecuaciones K * u = F...")
    K_free = K_global[np.ix_(free_dofs, free_dofs)]
    F_free = F_global[free_dofs]

    u_free = np.linalg.solve(K_free, F_free)

    u_global = np.zeros(n_dofs)
    u_global[free_dofs] = u_free

    # --------------------------------------------------------------------------
    # 7. PROMEDIADO NODAL DE TENSIONES Y EXPORTACIÓN A PARAVIEW
    # --------------------------------------------------------------------------
    print("--> Calculando tensiones promediadas y generando archivo para ParaView...")
    sigma_nodal = compute_nodal_stresses(mesh, elements, u_global)

    write_results_2d(output_file, mesh, u_global, nodal_stresses=sigma_nodal)
    print(f"¡Proceso finalizado! Resultado guardado en: {output_file}")


if __name__ == "__main__":
    main()