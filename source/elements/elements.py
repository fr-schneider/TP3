class tri_cst:
#tri_cst: triangular de 3 nodos (cst).
    def __init__(self,global_id,node_global_ids,material,tipo):
        self.id=global_id
        self.node=node_global_ids
        self.n_node=3
        self.dofs_node=2
        self.ndofs=6
        self.mat=material
        self.C=material.constitutive(tipo)

    def shape(self,r,s,N):
        if N==1:
            return 1-r-s
        elif N==2:
            return r
        elif N==3:
            return s

    def dshape(self,r,s,N,):
        if N == 1:
            return [-1, -1]
        elif N == 2:
            return [1, 0]
        elif N == 3:
            return [0, 1]

class tri_lst:
#triangular de 6 nodos.
    def __init__(self, global_id, node_global_ids,material,tipo):
        self.id = global_id
        self.node = node_global_ids
        self.n_node=6
        self.dofs_node = 2
        self.ndofs = 12
        self.mat=material
        self.C=material.constitutive(tipo)

    def shape(self, r, s, N):
        L1 = 1 - r - s
        L2 = r
        L3 = s

        shapes = [
            L1 * (2 * L1 - 1),
            L2 * (2 * L2 - 1),  
            L3 * (2 * L3 - 1), 
            4 * L1 * L2,     
            4 * L2 * L3,      
            4 * L3 * L1,        
        ]
        return shapes[N - 1]

    def dshape(self, r, s, N):
        
        dshapes = [
            [4 * r + 4 * s - 3, 4 * r + 4 * s - 3],
            [4 * r - 1, 0],                           
            [0, 4 * s - 1],    
            [4 - 8 * r - 4 * s, -4 * r],       
            [4 * s, 4 * r],                           
            [-4 * s, 4 - 4 * r - 8 * s],              
        ]
        return dshapes[N - 1]
    
class quad_q4:
    # Cuadrilatero de 4 nodos (Q4)
    def __init__(self, global_id, node_global_ids,material,tipo):
        self.id = global_id
        self.node = node_global_ids
        self.n_node = 4
        self.dofs_node = 2
        self.ndofs = 8
        self.mat=material
        self.C=material.constitutive(tipo)


    def shape(self, r, s, N):
        shapes = [
            0.25 * (1-r) * (1-s),
            0.25 * (1+r) * (1-s),
            0.25 * (1+r) * (1+s),
            0.25 * (1-r) * (1+s)
        ]
        return shapes[N-1]

    def dshape(self, r, s, N):
        dshapes = [
            [-0.25*(1-s), -0.25*(1-r)],
            [ 0.25*(1-s), -0.25*(1+r)],
            [ 0.25*(1+s),  0.25*(1+r)],
            [-0.25*(1+s),  0.25*(1-r)]
        ]
        return dshapes[N-1]
    
    
class quad_q8:
    # Cuadrilatero de 8 nodos (Q8)

    def __init__(self, global_id, node_global_ids, material, tipo):
        self.id = global_id
        self.node = node_global_ids
        self.n_node = 8
        self.dofs_node = 2
        self.ndofs = 16
        self.mat = material
        self.C = material.constitutive(tipo)

    def shape(self, r, s, N):

        N5 = 0.5 * (1-r**2) * (1-s)
        N6 = 0.5 * (1+r) * (1-s**2)
        N7 = 0.5 * (1-r**2) * (1+s)
        N8 = 0.5 * (1-r) * (1-s**2)
        N1 = 0.25 * (1-r) * (1-s) - 0.5 * (N8+N5)
        N2 = 0.25 * (1+r) * (1-s) - 0.5 * (N5+N6)
        N3 = 0.25 * (1+r) * (1+s) - 0.5 * (N6+N7)
        N4 = 0.25 * (1-r) * (1+s) - 0.5 * (N7+N8)

        shapes = [N1, N2, N3, N4, N5, N6, N7, N8]
        return shapes[N-1]

    def dshape(self, r, s, N):

        dshapes = [
            [0.25*(1-s)*(2*r+s), 0.25*(1-r)*(r+2*s)],
            [0.25*(1-s)*(2*r-s), 0.25*(1+r)*(2*s-r)],
            [0.25*(1+s)*(2*r+s), 0.25*(1+r)*(r+2*s)],
            [0.25*(1+s)*(2*r-s), 0.25*(1-r)*(2*s-r)],
            [-r*(1-s), -0.5*(1-r**2)],
            [0.5*(1-s**2), -(1+r)*s],
            [-r*(1+s), 0.5*(1-r**2)],
            [-0.5*(1-s**2), -(1-r)*s]
        ]

        return dshapes[N-1]


class quad_q9:
    # Cuadrilatero de 9 nodos (Q9)

    def __init__(self, global_id, node_global_ids, material, tipo):
        self.id = global_id
        self.node = node_global_ids
        self.n_node = 9
        self.dofs_node = 2
        self.ndofs = 18
        self.mat = material
        self.C = material.constitutive(tipo)

    def shape(self, r, s, N):

        N9 = (1-r**2)*(1-s**2)
        N5 = 0.5*(1-r**2)*(1-s) - 0.5*N9
        N6 = 0.5*(1+r)*(1-s**2) - 0.5*N9
        N7 = 0.5*(1-r**2)*(1+s) - 0.5*N9
        N8 = 0.5*(1-r)*(1-s**2) - 0.5*N9
        N1 = 0.25*r*(r-1)*s*(s-1)
        N2 = 0.25*r*(r+1)*s*(s-1)
        N3 = 0.25*r*(r+1)*s*(s+1)
        N4 = 0.25*r*(r-1)*s*(s+1)

        shapes = [N1, N2, N3, N4, N5, N6, N7, N8, N9]

        return shapes[N-1]

    def dshape(self, r, s, N):

        dshapes = [
            [0.25*(2*r-1)*s*(s-1),
             0.25*r*(r-1)*(2*s-1)],
            [0.25*(2*r+1)*s*(s-1),
             0.25*r*(r+1)*(2*s-1)],
            [0.25*(2*r+1)*s*(s+1),
             0.25*r*(r+1)*(2*s+1)],
            [0.25*(2*r-1)*s*(s+1),
             0.25*r*(r-1)*(2*s+1)],
            [-r*s*(s-1),
             0.5*(1-r**2)*(2*s-1)],
            [0.5*(2*r+1)*(1-s**2),
             -r*(r+1)*s],
            [-r*s*(s+1),
             0.5*(1-r**2)*(2*s+1)],
            [0.5*(2*r-1)*(1-s**2),
             -r*(r-1)*s],
            [-2*r*(1-s**2),
             -2*s*(1-r**2)]
        ]

        return dshapes[N-1]