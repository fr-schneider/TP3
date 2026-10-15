from source.mesh.node import Node
from source.solvers.solvers import node_stress, ke
import numpy as np

class Mesh:
    def __init__(self):
        self.dim = 0
        self.dofs_per_node = 0
        self.nodes = []
        self.element_global_ids = []
        self.element_types = []
        self.connectivity = []

    def read(self, filename):
        file = open(filename, "r")

        lines = []

        for raw_line in file:
            line = raw_line.split("#")[0]
            line = line.strip()

            if line != "":
                lines.append(line)

        file.close()

        header = lines[0].split()

        n_nodes = int(header[0])
        n_elements = int(header[1])
        dim = int(header[2])

        self.dim = dim

        line_number = 1

        for i in range(n_nodes):
            values = lines[line_number].split()
            line_number = line_number + 1

            global_id = int(values[0]) - 1

            if global_id != i:
                raise ValueError("Unexpected node numbering")

            x = float(values[1])
            y = float(values[2])

            self.nodes.append(Node(global_id, x, y))

        for i in range(n_elements):
            values = lines[line_number].split()
            line_number = line_number + 1

            element_global_id = int(values[0]) - 1
            element_type = values[1]

            element_conn = []

            for j in range(2, len(values)):
                element_conn.append(int(values[j]) - 1)

            self.element_global_ids.append(element_global_id)
            self.element_types.append(element_type)
            self.connectivity.append(element_conn)

    def number_of_nodes(self):
        return len(self.nodes)

    def number_of_elements(self):
        return len(self.connectivity)

    def node(self, node_global_id):
        return self.nodes[node_global_id]

    def number_of_dofs(self):
        return self.number_of_nodes() * self.dofs_per_node

    def node_global_dofs(self, node_global_id):
        dofs = []

        for i in range(self.dofs_per_node):
            dofs.append(node_global_id * self.dofs_per_node + i)

        return dofs



    def compute_nodal_stresses(self, elements, u):
        """
        Calcula las tensiones continuas promediadas en los nodos de la malla
        extrapolando las tensiones de cada elemento (Cook, Sec. 3.12).
        
        Devuelve:
        - sigma_nodal: Matriz (3 x n_nodos) con [sigma_x, sigma_y, tau_xy] promediados.
        """
        n_nodes = self.number_of_nodes()
        sigma_global = np.zeros((3, n_nodes))
        node_count = np.zeros(n_nodes)

        for element in elements:
            # 1. Coordenadas y desplazamientos locales del elemento
            node_ids = element.node
            coord_node = np.array([[self.node(i).x, self.node(i).y] for i in node_ids])
            
            # Extraemos los DOFs del elemento desde u global
            dofs = np.hstack([[2 * g, 2 * g + 1] for g in node_ids])
            u_e = u[dofs]

            # 2. Tu funcion node_stress evalua las tensiones en los nodos del elemento (3 x n_nodos_elem)
            sigma_elemento = node_stress(element, coord_node, u_e)

            # 3. Acumulamos las tensiones y el conteo de apariciones de cada nodo
            for j, g_id in enumerate(node_ids):
                sigma_global[:, g_id] += sigma_elemento[:, j]
                node_count[g_id] += 1.0

        # 4. Promediamos dividiendo por la cantidad de elementos compartidos por cada nodo
        sigma_nodal = sigma_global / node_count

        return sigma_nodal
    
    def assemble_stiffness(self, elements, orden, ):
   
    #Ensambla la matriz de rigidez global K_global usando Ke existente.
   
        n_dofs = self.number_of_dofs()
        K_global = np.zeros((n_dofs, n_dofs))
        coord_global = np.array([[node.x, node.y] for node in self.nodes])

        for elem in elements:
            coord_node = coord_global[elem.node]
            
            
            Ke = ke(elem,orden,coord_node)  

            # Mapeo a K_global usando los DOFs globales de tu Mesh
            for i_loc, g_node_i in enumerate(elem.node):
                dofs_i = self.node_global_dofs(g_node_i)
                for j_loc, g_node_j in enumerate(elem.node):
                    dofs_j = self.node_global_dofs(g_node_j)

                    # Colocamos el bloque 2x2 del par de nodos en K_global
                    K_global[dofs_i[0]:dofs_i[-1]+1, dofs_j[0]:dofs_j[-1]+1] += Ke[2*i_loc:2*i_loc+2, 2*j_loc:2*j_loc+2]

        return K_global
    
    def create_elements(self, material, tipo="Plane_strain", elem_map=None):
       
        #Instancia y devuelve la lista de objetos elemento asociados a la malla.
        
        # Mapa por defecto si no se pasa uno personalizado
        if elem_map is None:
            from elements import tri_cst, tri_lst, quad_q4, quad_q8, quad_q9

            elem_map = {
                "cst": tri_cst,
                "tri3": tri_cst,
                "lst": tri_lst,
                "tri6": tri_lst,
                "q4": quad_q4,
                "quad4": quad_q4,
                "q8": quad_q8,
                "quad8": quad_q8,
                "q9": quad_q9,
                "quad9": quad_q9,
            }

        coord_global = np.array([[node.x, node.y] for node in self.nodes])
        elements = []

        for i in range(self.number_of_elements()):
            g_id = self.element_global_ids[i]
            elem_type = self.element_types[i].lower()
            conn = self.connectivity[i]

            cls = elem_map[elem_type]
            elem = cls(g_id, conn, material, tipo)
            elem.coord_node = coord_global[conn]
            elements.append(elem)

        return elements