
import ast
from pathlib import Path
import pandas as pd

TIPOS = {
    "CTRIA3": "tri3",
    "CPLSTN3": "tri3",
    "CTRIA6": "tri6",
    "CPLSTN6": "tri6",
    "CQUAD4": "q4",
    "CQUAD8": "q8",
}

NODOS_POR_TIPO = {
    "tri3": 3,
    "tri6": 6,
    "q4": 4,
    "q8": 8,
}


def csv_to_mesh_txt(carpeta, prefijo, salida, estado="Plane_strain"):
    carpeta = Path(carpeta)
    salida = Path(salida)

    nodes = pd.read_csv(carpeta / f"{prefijo}_nodes.csv")
    elems = pd.read_csv(carpeta / f"{prefijo}_elements.csv")
    props = pd.read_csv(carpeta / f"{prefijo}_properties.csv")
    mats = pd.read_csv(carpeta / f"{prefijo}_materials.csv")

    nid_a_idx = {
        int(nid): i + 1
        for i, nid in enumerate(nodes["nid"])
    }

    pid_a_mid = {}

    for _, row in props.iterrows():
        if pd.isna(row["mid"]):
            raise ValueError(f"Propiedad {row['pid']} sin material")

        pid_a_mid[int(row["pid"])] = int(row["mid"])

    lineas_elem = []
    mids_usados = set()

    for _, row in elems.iterrows():
        tipo_nas = str(row["type"]).upper()

        if tipo_nas not in TIPOS:
            raise ValueError(f"Elemento no soportado: {tipo_nas}")

        tipo = TIPOS[tipo_nas]
        conn = ast.literal_eval(row["nodes"])

        if len(conn) != NODOS_POR_TIPO[tipo]:
            raise ValueError(f"Conectividad incorrecta para {tipo}")

        if any(n is None for n in conn):
            raise ValueError("Hay nodos faltantes")

        mid = pid_a_mid[int(row["pid"])]
        mids_usados.add(mid)

        nodos = " ".join(
            str(nid_a_idx[int(n)]) for n in conn
        )

        eid = len(lineas_elem) + 1
        lineas_elem.append(f"{eid} {tipo} {mid} {nodos}")

    lineas_mat = []
    materiales_encontrados = set()

    for _, row in mats.iterrows():
        mid = int(row["mid"])

        if mid not in mids_usados:
            continue

        E = float(row["E"])
        nu = row["nu"]

        if pd.isna(nu):
            nu = E / (2.0 * float(row["G"])) - 1.0

        lineas_mat.append(f"{mid} {E} {float(nu)} {estado}")
        materiales_encontrados.add(mid)

    if materiales_encontrados != mids_usados:
        raise ValueError("Faltan materiales en el CSV")

    lineas_nodo = [
        f"{nid_a_idx[int(row['nid'])]} {row['x']} {row['y']}"
        for _, row in nodes.iterrows()
    ]

    salida.parent.mkdir(parents=True, exist_ok=True)

    with open(salida, "w", encoding="utf-8") as f:
        f.write("# n_nodos n_elementos dim n_materiales\n")
        f.write(f"{len(lineas_nodo)} {len(lineas_elem)} 2 {len(lineas_mat)}\n")

        f.write("# materiales\n")
        f.write("\n".join(lineas_mat) + "\n")

        f.write("# nodos\n")
        f.write("\n".join(lineas_nodo) + "\n")

        f.write("# elementos\n")
        f.write("\n".join(lineas_elem) + "\n")

    print("Archivo generado:", salida)
    print("Nodos:", len(lineas_nodo))
    print("Elementos:", len(lineas_elem))
    print("Materiales:", len(lineas_mat))


if __name__ == "__main__":
    SCRIPT_DIR = Path(__file__).resolve().parent

    csv_to_mesh_txt(
        carpeta=SCRIPT_DIR / "mesh_dat_files",
        prefijo="represa",
        salida=SCRIPT_DIR.parent / "meshes" / "represa.txt",
        estado="Plane_strain"
    )
