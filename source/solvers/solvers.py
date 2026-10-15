from .gauss_integration import gauss
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





def node_stress(element, coord_node, d):
    n_nodos = element.n_node

    # T3: tensiones constantes
    if n_nodos == 3:
        orden = 1
        sigma_gauss, q = gauss_stress(
            element, orden, coord_node, d
        )
        return np.repeat(sigma_gauss, 3, axis=1)

    # Q4, Q8 y Q9
    elif n_nodos in (4, 8, 9):
        orden = 3 if n_nodos == 9 else 2

        sigma_gauss, q = gauss_stress(
            element, orden, coord_node, d
        )

        q = np.asarray(q)
        coord = np.asarray(element.node_nat_coord)

        def base(x, y):
            if n_nodos == 9:
                return np.column_stack((
                    np.ones_like(x), x, y, x*y,
                    x**2, y**2, x**2*y,
                    x*y**2, x**2*y**2
                ))

            return np.column_stack((
                np.ones_like(x), x, y, x*y
            ))

        P = base(q[:, 0], q[:, 1])
        Pn = base(coord[:, 0], coord[:, 1])

        E = Pn @ np.linalg.inv(P)

        return sigma_gauss @ E.T

    # T6
    elif n_nodos == 6:
        orden = 2

        sigma_gauss, q = gauss_stress(
            element, orden, coord_node, d
        )

        q = np.asarray(q)
        coord = np.asarray(element.node_nat_coord)

        P = np.column_stack((
            np.ones(3), q[:, 0], q[:, 1]
        ))

        Pn = np.column_stack((
            np.ones(n_nodos),
            coord[:, 0],
            coord[:, 1]
        ))

        E = Pn @ np.linalg.inv(P)

        return sigma_gauss @ E.T

    else:
        raise ValueError("Tipo de elemento no soportado")
