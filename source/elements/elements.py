class tri_cst:
#tri_cst: triangular de 3 nodos (cst).
    def __init__(self,global_id,node_global_ids):
        self.id=global_id
        self.node=node_global_ids
        self.n_node=3
        self.dofs_node=2
        self.ndofs=6


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
    def __init__(self, global_id, node_global_ids):
        self.id = global_id
        self.node = node_global_ids
        self.n_node=3
        self.dofs_node = 2
        self.ndofs = 12

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
    def __init__(self, global_id, node_global_ids):
        self.id = global_id
        self.node = node_global_ids
        self.n_node = 4
        self.dofs_node = 2
        self.ndofs = 8

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
    # Cuadrilatero de 8 nodos (Q4)
    def __init__(self, global_id, node_global_ids):
        self.id = global_id
        self.node = node_global_ids
        self.n_node = 8
        self.dofs_node = 2
        self.ndofs = 14

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

