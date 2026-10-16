from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from source.mesh.node import Node
from source.materials.linear_elastic import LinearElastic
from source.solvers.solvers import node_stress, ke

# Cantidad de nodos que exige cada tipo de elemento
NODOS_POR_TIPO = {
    "cst": 3, "tri3": 3,
    "lst": 6, "tri6": 6,
    "q4": 4, "quad4": 4,
    "q8": 8, "quad8": 8,
    "q9": 9, "quad9": 9,
}

# Formato del archivo:
#   n_nodos n_elementos dim n_materiales
#   (n_materiales lineas)  id E nu estado        ej: 1 21000 0.25 Plane_stress
#   (n_nodos lineas)       id x y
#   (n_elementos lineas)   id tipo id_material nodo1 nodo2 ...
# Los comentarios empiezan con '#'.


class Mesh:
    def __init__(self):
        self.dim = 0
        self.dofs_per_node = 2
        self.nodes = []
        self.materials = {}          # id -> (objeto material, estado)
        self.element_global_ids = []
        self.element_types = []
        self.element_materials = []  # id de material por elemento
        self.connectivity = []

    def read(self, filename):
        lines = []
        with open(filename, "r") as file:
            for raw_line in file:
                line = raw_line.split("#")[0].strip()
                if line != "":
                    lines.append(line)

        header = lines[0].split()
        n_nodes = int(header[0])
        n_elements = int(header[1])
        self.dim = int(header[2])
        n_materials = int(header[3])
        self.dofs_per_node = 2

        line_number = 1

        # Materiales
        for _ in range(n_materials):
            v = lines[line_number].split()
            line_number += 1
            mat_id = int(v[0])
            self.materials[mat_id] = (
                LinearElastic(float(v[1]), float(v[2])),
                v[3],
            )

        # Nodos
        for i in range(n_nodes):
            values = lines[line_number].split()
            line_number += 1

            global_id = int(values[0]) - 1
            if global_id != i:
                raise ValueError("Unexpected node numbering")

            self.nodes.append(Node(global_id, float(values[1]), float(values[2])))

        # Elementos
        for _ in range(n_elements):
            values = lines[line_number].split()
            line_number += 1

            element_type = values[1].lower()
            mat_id = int(values[2])
            conn = [int(x) - 1 for x in values[3:]]

            if element_type not in NODOS_POR_TIPO:
                raise ValueError(f"Elemento {values[0]}: tipo '{element_type}' no soportado")
            if len(conn) != NODOS_POR_TIPO[element_type]:
                raise ValueError(
                    f"Elemento {values[0]}: {element_type} requiere "
                    f"{NODOS_POR_TIPO[element_type]} nodos, hay {len(conn)}"
                )
            if mat_id not in self.materials:
                raise ValueError(f"Elemento {values[0]}: material {mat_id} no definido")

            self.element_global_ids.append(int(values[0]) - 1)
            self.element_types.append(element_type)
            self.element_materials.append(mat_id)
            self.connectivity.append(conn)

    def number_of_nodes(self):
        return len(self.nodes)

    def number_of_elements(self):
        return len(self.connectivity)

    def node(self, node_global_id):
        return self.nodes[node_global_id]

    def number_of_dofs(self):
        return self.number_of_nodes() * self.dofs_per_node

    def node_global_dofs(self, node_global_id):
        return [node_global_id * self.dofs_per_node + i
                for i in range(self.dofs_per_node)]

    def compute_nodal_stresses(self, elements, u, tipo="Plane_strain"):
        """Calcula el campo continuo de tensiones nodales promediado globalmente,

        recibiendo 'tipo' ('Plane_strain' / 'Plane_stress') y 'nu' como argumentos.
        """
        n_nodes = self.number_of_nodes()
        sigma_global = np.zeros((3, n_nodes))  # [sigma_x, sigma_y, tau_xy]
        node_count = np.zeros(n_nodes)

        for element in elements:
            node_ids = element.node
            coord_node = np.array(
                [[self.node(i).x, self.node(i).y] for i in node_ids]
            )

            dofs = np.hstack([self.node_global_dofs(g) for g in node_ids])
            u_e = u[dofs]

            # Llamada a node_stress original con sus 3 argumentos
            sigma_elemento = node_stress(element, coord_node, u_e)

            for j, g_id in enumerate(node_ids):
                sigma_global[:, g_id] += sigma_elemento[:, j]
                node_count[g_id] += 1.0

        # Promedio nodal (suavizado)
        sigma_nodal = np.zeros_like(sigma_global)
        mask = node_count > 0
        sigma_nodal[:, mask] = sigma_global[:, mask] / node_count[mask]

        sx = sigma_nodal[0, :]
        sy = sigma_nodal[1, :]
        txy = sigma_nodal[2, :]

        # Calculo de sigma_z según el tipo de análisis recibido
        nu=element.mat.nu
        es_deformacion_plana = str(tipo).lower() == "plane_strain"
        sz = nu * (sx + sy) if es_deformacion_plana else np.zeros_like(sx)

        # Tensión equivalente de Von Mises (3D completa)
        von_mises = np.sqrt(
            sx**2
            + sy**2
            + sz**2
            - sx * sy
            - sy * sz
            - sz * sx
            + 3.0 * (txy**2)
        )

        return {
            "sigma_x": sx,
            "sigma_y": sy,
            "sigma_z": sz,
            "tau_xy": txy,
            "von_mises": von_mises,
        }


    def assemble_stiffness(self, elements, orden=2, espesor=1.0):
        """Ensambla K global con el orden de integracion dado."""
        n_dofs = self.number_of_dofs()
        K_global = np.zeros((n_dofs, n_dofs))
        coord_global = np.array([[n.x, n.y] for n in self.nodes])

        for elem in elements:
            coord_node = coord_global[elem.node]
            Ke = espesor * ke(elem, orden, coord_node)

            for i_loc, g_i in enumerate(elem.node):
                dofs_i = self.node_global_dofs(g_i)
                for j_loc, g_j in enumerate(elem.node):
                    dofs_j = self.node_global_dofs(g_j)
                    K_global[np.ix_(dofs_i, dofs_j)] += Ke[
                        2 * i_loc: 2 * i_loc + 2,
                        2 * j_loc: 2 * j_loc + 2,
                    ]

        return K_global

    def create_elements(self, elem_map=None):
        """Instancia los elementos usando el material y tipo de cada uno."""
        if elem_map is None:
            from source.elements.elements import tri_cst, tri_lst, quad_q4, quad_q8, quad_q9

            elem_map = {
                "cst": tri_cst, "tri3": tri_cst,
                "lst": tri_lst, "tri6": tri_lst,
                "q4": quad_q4, "quad4": quad_q4,
                "q8": quad_q8, "quad8": quad_q8,
                "q9": quad_q9, "quad9": quad_q9,
            }

        coord_global = np.array([[node.x, node.y] for node in self.nodes])
        elements = []

        for i in range(self.number_of_elements()):
            material, tipo = self.materials[self.element_materials[i]]
            conn = self.connectivity[i]
            cls = elem_map[self.element_types[i]]

            elem = cls(self.element_global_ids[i], conn, material, tipo)
            elem.coord_node = coord_global[conn]
            elements.append(elem)

        return elements
    
if __name__ == "__main__":

    mesh = Mesh()

    archivo = Path(__file__).resolve().parents[2] / "meshes" / "represa.txt"

    mesh.read(archivo)

    print("Archivo leído correctamente")
    print("Cantidad de nodos:", mesh.number_of_nodes())
    print("Cantidad de elementos:", mesh.number_of_elements())
    print("Dimensión:", mesh.dim)
    print("Materiales:", mesh.materials)
    print("Tipos de elementos:", set(mesh.element_types))

    print("\nPrimeros 5 nodos:")
    for nodo in mesh.nodes[:5]:
        print(nodo.x, nodo.y)

    print("\nPrimeras 5 conectividades:")
    for conn in mesh.connectivity[:5]:
        print(conn)
