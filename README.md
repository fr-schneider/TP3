Analisis por elementos finitos de una represa de gravedad
Para pensar
Condiciones de Borde en el suelo
fijo atras y friccion en el medio? hay deformaicon en el suelo?
    Si es friccion hay que ver como evaluar el peso que actua punto a punto

apoyo fijo a la derechay modvil en el resto de puntos

Modelado
Tension plana vs Deformacion plana?
si asumo el laterales infinitamente rigidos deformacion plana (creo que la mas logica)
Esperar F roce y futuras clases

Documentacion:
Elementos
    tri_cst: triangular de 3 nodos (cst)
    tri_lst: triangular de 6 noodos (lst) (1,2,3) vertices, (4,5,6) los del medio.
Funciones de forma:
    shape(self,r,s,N):
        Self: elemento
        r y s, coordenadas para q de integracion
        N nro natural N1 es 1, N2 es 2...
        Devuelve el resultado de evaluar la funcion en (r,s)

    dshape(self,r,s,N,):
        Self: elemento
        r y s, coordenadas para q de integracion
        N nro natural N1 es 1, N2 es 2...
        Devuelve un vector donde cada posicion son [dN1/dr(r,s) dN1/ds(r,s)] 