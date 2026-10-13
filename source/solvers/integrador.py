from gauss_integration import gauss
import numpy as np

def ke(element,orden,coord_node):
    C=element.C
    ke=np.zeros((element.ndofs,element.ndofs))
    q,w=gauss(element,orden)
    for qi in q:
        r=qi[0]
        s=qi[1]
        j,det_j,inv_j=jacobiano(element,r,s,coord_node)
        B=b_matrix(element,r,s,inv_j)
        #aca esta la formula clave, hay que chequear que python entienda esto
        ke=ke+(B.T@C@B)*det_j*w
    return ke
        


def b_matrix(element,r,s,inv_j):
    n_nodos=element.n_node

    #Inicio la matriz B, 3 filas por las deformaciones y una columna por dof
    B=np.zeros((3,element.ndofs))

    #Para cada fuoncion de forma evaluo las derivadas y organizo en array
    for i in range(n_nodos):
        dN=np.array(element.dshape(r,s,i+1))

        #Multiplico por inv_j para tener en coordenadas globales
        dN=np.dot(inv_j,dN)

        #Introduzco en B
        col=2*i
        B[0,col]=dN[0]
        B[1,col+1]=dN[1]
        B[2,col]=dN[1]
        B[2,col+1]=dN[0]
    return B





def jacobiano(element,r,s,coord_node):
    n_nodos=element.n_node
    dN_dr=np.zeros(n_nodos)
    dN_ds=np.zeros(n_nodos)
    # Evaluo las derivadas
    for i in range(n_nodos):
        dN_dr[i],dN_ds[i]=element.dshape(r,s,i+1)

    #Preparo las coordenadas
    x=coord_node[:,0]
    y=coord_node[:,1] 

    #MUltiplico vectores (fila por columna)
    dx_dr = np.dot(dN_dr, x)
    dy_dr = np.dot(dN_dr, y)
    dx_ds = np.dot(dN_ds, x)
    dy_ds = np.dot(dN_ds, y)

    #Armo el J y detJ y J-1
    j=np.array([[dx_dr, dy_dr], [dx_ds, dy_ds]])
    det_j=np.linalg.det(j)
    inv_j=np.linalg.inv(j)
    return j,det_j,inv_j