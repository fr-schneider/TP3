from pathlib import Path
from pyNastran.bdf.bdf import BDF
import pandas as pd
import numpy as np


# Unidades del modelo convertido: m, N, Pa
E = 21000000000.0       # 21 GPa (hormigon)
nu = 0.2
espesor = 1.0          # 1 m para el analisis 2D
mid = 900000
ESCALA_COORDENADAS = 1.0  


def nastran_mesh_to_csv(bdf_file, output_dir, file_prefix):
    bdf_file = Path(bdf_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    lineas = bdf_file.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines(keepends=True)

    inicio = next(
        (
            i for i, linea in enumerate(lineas)
            if linea.strip().upper() == "BEGIN BULK"
        ),
        None,
    )
    if inicio is None:
        raise ValueError(f"No se encontro BEGIN BULK en {bdf_file}")

    conversion = {
        "CPLSTN3": "CTRIA3",
        "CPLSTN6": "CTRIA6",
        "CPLSTN4": "CQUAD4",
        "CPLSTN8": "CQUAD8",
    }
    tipos = set(conversion.values())

    bulk_lines = []
    propiedades = set()
    leyendo_tarjeta = False

    for linea in lineas[inicio + 1:]:
        nombre = linea[:8].strip().upper()

        if nombre == "ENDDATA":
            break

        if linea.lstrip().startswith("$"):
            continue

        # Conservar los nodos y sus lineas de continuacion.
        if nombre in ("GRID", "GRID*"):
            bulk_lines.append(linea)
            leyendo_tarjeta = True
            continue

        # Tambien conservar continuaciones de elementos Q8 y T6.
        # Pueden empezar por '+', '*' o tener el primer campo vacio.
        if leyendo_tarjeta and (
            nombre.startswith("+")
            or nombre.startswith("*")
            or nombre == ""
        ):
            bulk_lines.append(linea)
            continue

        leyendo_tarjeta = False

        if nombre in conversion:
            tipo = conversion[nombre]
        elif nombre in tipos:
            tipo = nombre
        else:
            continue

        pid = int(linea[16:24])
        propiedades.add(pid)
        bulk_lines.append(tipo.ljust(8) + linea[8:])
        leyendo_tarjeta = True

    for pid in sorted(propiedades):
        bulk_lines.append(f"PSHELL,{pid},{mid},{espesor}\n")

    bulk_lines.append(f"MAT1,{mid},{E},,{nu}\n")

    bulk_file = output_dir / f"{file_prefix}_bulk.dat"
    bulk_file.write_text("".join(bulk_lines), encoding="utf-8")

    model = BDF(debug=False)
    model.read_bdf(
        str(bulk_file),
        punch=True,
        xref=True,
        validate=True,
    )

    # Verificar que no se hayan perdido nodos de los Q8/T6.
    cantidad_nodos = {
        "CTRIA3": 3,
        "CTRIA6": 6,
        "CQUAD4": 4,
        "CQUAD8": 8,
    }
    for eid, elem in model.elements.items():
        conectividad = list(elem.node_ids)
        esperado = cantidad_nodos.get(elem.type)
        if esperado is not None and (
            len(conectividad) != esperado
            or any(n is None for n in conectividad)
        ):
            raise ValueError(
                f"Elemento {eid} ({elem.type}) con conectividad "
                f"incompleta: {conectividad}. Revisar las lineas "
                "de continuacion en el .dat original exportado por NX."
            )

    nids = sorted(model.nodes.keys())
    xyz = np.array([
        model.nodes[nid].get_position()
        for nid in nids
    ], dtype=float) * ESCALA_COORDENADAS

    df_nodes = pd.DataFrame({
        "nid": nids,
        "x": xyz[:, 0],
        "y": xyz[:, 1],
        "z": xyz[:, 2],
    })
    df_nodes.to_csv(
        output_dir / f"{file_prefix}_nodes.csv", index=False
    )

    df_elements = pd.DataFrame([
        {
            "eid": eid,
            "type": elem.type,
            "pid": elem.Pid(),
            "nodes": list(elem.node_ids),
        }
        for eid, elem in sorted(model.elements.items())
    ])
    df_elements.to_csv(
        output_dir / f"{file_prefix}_elements.csv", index=False
    )

    df_properties = pd.DataFrame([
        {
            "pid": pid,
            "type": prop.type,
            "mid": prop.mid1,
            "t": prop.t,
        }
        for pid, prop in sorted(model.properties.items())
    ])
    df_properties.to_csv(
        output_dir / f"{file_prefix}_properties.csv", index=False
    )

    df_materials = pd.DataFrame([
        {
            "mid": mid_id,
            "E": mat.E(),
            "G": mat.G(),
            "nu": mat.Nu(),
        }
        for mid_id, mat in sorted(model.materials.items())
    ])
    df_materials.to_csv(
        output_dir / f"{file_prefix}_materials.csv", index=False
    )

    print("Nodos:", len(df_nodes))
    print("Elementos:", len(df_elements))
    print("Propiedades:", len(df_properties))
    print("Materiales:", len(df_materials))
    print("Base:", np.ptp(xyz[:, 0]), "m")
    print("Altura:", np.ptp(xyz[:, 1]), "m")
    print("Conversion finalizada")


# CONVERTIR TODAS LAS MALLAS DE LA REPRESA
if __name__ == "__main__":
    SCRIPT_DIR = Path(__file__).resolve().parent

    archivos = [
        "Represa_T3.dat",
        "represa_Q4_gruesa.dat",
        "Represa_Q8.dat",
        "Represa_T6.dat",
        "Represa_Q8_ModificacionCombinada.dat",
        "Represa_Q8_Modificacion1.dat",
    ]

    for archivo in archivos:
        ruta = SCRIPT_DIR / "mesh_dat_files" / archivo
        prefijo = Path(archivo).stem

        print("\nConvirtiendo:", archivo)
        nastran_mesh_to_csv(
            bdf_file=ruta,
            output_dir=SCRIPT_DIR / "mesh_dat_files",
            file_prefix=prefijo,
        )
