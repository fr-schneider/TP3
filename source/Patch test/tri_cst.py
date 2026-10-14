import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from source.elements.elements import tri_cst
from source.materials.linear_elastic import LinearElastic

nodos = [
    [0, 0],
    [100, 0],
    [100, 100],
    [0, 100],
    [45, 55]
]

triangulo = [
    [0, 1, 4],
    [4, 1, 2],
    [2, 3, 4],
    [3, 0, 4]
]

material = LinearElastic(21000, 0.25)
tipo = "Plane_stress"

elemento1 = tri_cst(1, [0, 1, 4], material, tipo)
elemento2 = tri_cst(2, [4, 1, 2], material, tipo)
elemento3 = tri_cst(3, [2, 3, 4], material, tipo)
elemento4 = tri_cst(4, [3, 0, 4], material, tipo)

elementos = [elemento1, elemento2, elemento3, elemento4]

for elemento in elementos:
    dN1 = elemento.dshape(0, 0, 1)
    dN2 = elemento.dshape(0, 0, 2)
    dN3 = elemento.dshape(0, 0, 3)

    dN = [dN1, dN2, dN3]

    coordenadas = [nodos[i] for i in elemento.node]
    
    x1, y1 = coordenadas[0]
    x2, y2 = coordenadas[1]
    x3, y3 = coordenadas[2]
    dp=[[0,0],
        [0,0],
        [0,0]]
    J = [
        [
            dN1[0]*x1 + dN2[0]*x2 + dN3[0]*x3,
            dN1[0]*y1 + dN2[0]*y2 + dN3[0]*y3
        ],
        [
            dN1[1]*x1 + dN2[1]*x2 + dN3[1]*x3,
            dN1[1]*y1 + dN2[1]*y2 + dN3[1]*y3
        ]
    ]
    detJ=J[0][0]*J[1][1]-J[1][0]*J[0][1]
    invJ = [
    [J[1][1]/detJ, -J[0][1]/detJ],
    [-J[1][0]/detJ, J[0][0]/detJ]
     ]
    
    dp[0][0] = invJ[0][0]*dN1[0] + invJ[0][1]*dN1[1]
    dp[0][1] = invJ[1][0]*dN1[0] + invJ[1][1]*dN1[1]
    dp[1][0] = invJ[0][0]*dN2[0] + invJ[0][1]*dN2[1]
    dp[1][1] = invJ[1][0]*dN2[0] + invJ[1][1]*dN2[1]
    dp[2][0] = invJ[0][0]*dN3[0] + invJ[0][1]*dN3[1]
    dp[2][1] = invJ[1][0]*dN3[0] + invJ[1][1]*dN3[1]
    
    B = [
        [dp[0][0], 0, dp[1][0], 0, dp[2][0], 0],
        [0, dp[0][1], 0, dp[1][1], 0, dp[2][1]],
        [dp[0][1], dp[0][0], dp[1][1], dp[1][0], dp[2][1], dp[2][0]]
    ]
    

    print("Elemento:", elemento.id)
    print("Coordenadas:", coordenadas)
    print("Derivadas:", dN)
    print("Jacobiano",J)
    print("Matriz",B)

