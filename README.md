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
    Material:
        LinearElastic(E,nu)
            E: Modulo de young
            nu: Coeficiente de Poisson
            self.(Tipo): devuelve el constitutivo tipo tension plano (tipo='Plane_stress') o deformacion plana (tipo='Plane_strain')

    Elementos:
        tri_cst: triangular de 3 nodos (cst)
        tri_lst: triangular de 6 noodos (lst) (1,2,3) vertices, (4,5,6) los del medio.
        quad_q4: Cuadrangular de 4 nodos
        quad_q8: Cuadrangular de 8 nodos
        quad_q9: Cuadrangular de 9 nodos


        Definicion:

            elemento(global_id,node_global_ids,material,tipo)
            global_id:nro de elemento (elemento.id)
            node_global_ids: id de los nodos (elemento.node)
            material: Clase de acuerdo a Linear_elastic
            tipo: analisis del tipo tension plano (tipo='Plane_stress') o deformacion plana (tipo='Plane_strain')

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

    Integracion

        q, w = gauss(element,orden)

            element es el elemento que se desea integrar
            orden es el orden de la aproximacion
            q son las coordenadas e los puntos de gauss es una tupla de tuplas, cada una son un par de coordenadas
            w peso de cada punto, es solo un numero

        ke(element,orden,coord_node)
            devuelve la matriz de rigidez del elemento usando gauss de orden n
            coord_node son las coordenadas de los nodos
    
    Tensiones
        gauss_stress(element,orden,coord_node,d)
            element es el elemento a calcular tensiones
            orden es el orden de integracion de gauss
            coord_node son las coordenadas globales de los nodos
            d son los desplazamientos nodales 

        node_stress(element,coord_node,d)

            element es el elemento a calcular tensiones
            coord_node son las coordenadas globales de los nodos
            d son los desplazamientos nodales
            devuelve una matriz con fila: sigmax,sigmay..., y cada columna un ndodo
            