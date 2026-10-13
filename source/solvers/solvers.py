from gauss_integration import gauss
import numpy as np

def ke(element,orden,coord_node):
    C=element.C
    ke=np.zeros((element.ndofs,element.ndofs))
    q,w=gauss(element,orden)
    for i, qi in enumerate(q):
        r=qi[0]
        s=qi[1]
        j,det_j,inv_j=jacobiano(element,r,s,coord_node)
        B=b_matrix(element,r,s,inv_j)
        #aca esta la formula clave, hay que chequear que python entienda esto
        ke=ke+(B.T@C@B)*det_j*w[i]
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

def gauss_stress(element,orden,coord_node,d):
    q, w=gauss(element,orden)
    
    #Calculo las tensiones en los puntos de superconvergencia(gauss)
    sigma=np.zeros((3,len(q)))
    for i, qi in enumerate(q):
        r=qi[0]
        s=qi[1]    
        #Calculo B
        j,det_j,inv_j=jacobiano(element,r,s,coord_node)
        B=b_matrix(element,r,s,inv_j)

        #Calculo la tension y la guardo en la posicion
        sigma[:,i]=element.C@B@d
    return sigma, q #Devuelve las tensiones de gaussen cada punto de gauss y las coordenadas




def node_stress(element,coord_node,d):
    n_nodos=element.n_node

    #Hago caso por caso para cada elemento ( esta medio hardcodedo esto pero bueno)
    if n_nodos==3:
        orden=1
        sigma_gauss, q=gauss_stress(element, orden, coord_node, d)
        return np.repeat(sigma_gauss, 3, axis=1) #repito para los 3 nodos la misma (constant)

    if n_nodos==4 or n_nodos==8 or n_nodos==9:

        # Para extrapolar uso segun el orden un factor a distinto
        if n_nodos==9:
            orden=3
            a=np.sqrt(5/3)
        else:
            orden=2
            a=np.sqrt(3)


        #Calculo las tensiones en los puntos de gauss
        sigma_gauss, q= gauss_stress(element, orden, coord_node, d)

        #preparo el bucle
        sigma=np.zeros((3,n_nodos))
        coord=element.node_nat_coord
        
        for i, (ksi,eta) in enumerate(coord):
            r=ksi*a #r es coord del nodo en el sistema*
            s=eta*a

            for j in range(len(q)):
                N_j=element.shape(r,s,j+1) #Evaluo r,s de sistema*
                sigma[:,i]+=sigma_gauss[:,j]*N_j #voy sumando
        return sigma


            

    