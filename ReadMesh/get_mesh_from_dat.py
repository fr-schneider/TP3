
from pathlib import Path
from pyNastran.bdf.bdf import BDF
import pandas as pd
import numpy as np


E = 210000.0
nu = 0.3
espesor = 1.0
mid = 900000


def nastran_mesh_to_csv(bdf_file, output_dir, file_prefix):
    bdf_file = Path(bdf_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    lineas = bdf_file.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines(keepends=True)

    inicio = next(
        (i for i, linea in enumerate(lineas)
         if linea.strip().upper() == "BEGIN BULK"),
        None
    )

    if inicio is None:
        raise ValueError("No se encontro BEGIN BULK")

    conversion = {
        "CPLSTN3": "CTRIA3",
        "CPLSTN6": "CTRIA6",
        "CPLSTN4": "CQUAD4",
        "CPLSTN8": "CQUAD8",
    }

    tipos = set(conversion.values())

    bulk_lines = []
    propiedades = set()
    leyendo_grid = False

    for linea in lineas[inicio + 1:]:
        nombre = linea[:8].strip().upper()

        if nombre == "ENDDATA":
            break

        if nombre in ("GRID", "GRID*"):
            bulk_lines.append(linea)
            leyendo_grid = True
            continue

        if nombre in ("*", "+") and leyendo_grid:
            bulk_lines.append(linea)
            continue

        leyendo_grid = False

        if nombre in conversion:
            tipo = conversion[nombre]
        elif nombre in tipos:
            tipo = nombre
        else:
            continue

        pid = int(linea[16:24])
        propiedades.add(pid)

        bulk_lines.append(tipo.ljust(8) + linea[8:])

    for pid in sorted(propiedades):
        bulk_lines.append(
            f"PSHELL,{pid},{mid},{espesor}\n"
        )

    bulk_lines.append(
        f"MAT1,{mid},{E},,{nu}\n"
    )

    bulk_file = output_dir / f"{file_prefix}_bulk.dat"

    bulk_file.write_text(
        "".join(bulk_lines),
        encoding="utf-8"
    )

    model = BDF(debug=False)

    model.read_bdf(
        str(bulk_file),
        punch=True,
        xref=True,
        validate=True
    )

    nids = sorted(model.nodes.keys())

    xyz = np.array([
        model.nodes[nid].get_position()
        for nid in nids
    ])

    df_nodes = pd.DataFrame({
        "nid": nids,
        "x": xyz[:, 0],
        "y": xyz[:, 1],
        "z": xyz[:, 2]
    })

    df_nodes.to_csv(
        output_dir / f"{file_prefix}_nodes.csv",
        index=False
    )

    df_elements = pd.DataFrame([
        {
            "eid": eid,
            "type": elem.type,
            "pid": elem.Pid(),
            "nodes": list(elem.node_ids)
        }
        for eid, elem in sorted(model.elements.items())
    ])

    df_elements.to_csv(
        output_dir / f"{file_prefix}_elements.csv",
        index=False
    )

    df_properties = pd.DataFrame([
        {
            "pid": pid,
            "type": prop.type,
            "mid": prop.mid1,
            "t": prop.t
        }
        for pid, prop in sorted(model.properties.items())
    ])

    df_properties.to_csv(
        output_dir / f"{file_prefix}_properties.csv",
        index=False
    )

    df_materials = pd.DataFrame([
        {
            "mid": mid_id,
            "E": mat.E(),
            "G": mat.G(),
            "nu": mat.Nu()
        }
        for mid_id, mat in sorted(model.materials.items())
    ])

    df_materials.to_csv(
        output_dir / f"{file_prefix}_materials.csv",
        index=False
    )

    print("Nodos:", len(df_nodes))
    print("Elementos:", len(df_elements))
    print("Propiedades:", len(df_properties))
    print("Materiales:", len(df_materials))
    print("Conversion finalizada")


if __name__ == "__main__":
    SCRIPT_DIR = Path(__file__).resolve().parent

    nastran_mesh_to_csv(
        bdf_file=SCRIPT_DIR / "mesh_dat_files" / "represa.dat",
        output_dir=SCRIPT_DIR / "mesh_dat_files",
        file_prefix="represa"
    )
