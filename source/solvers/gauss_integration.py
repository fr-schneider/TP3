import math

def gauss(element, n):
    if element.n_node==3:
        return ((1/3, 1/3),), 1
        
    if element.n_node==6:
        if n==1:
            return ((1/3, 1/3),), 1
        elif n==3:
            return ((2/3, 1/6),(1/6, 2/3),(1/6, 1/6)), 1/3

    if element.n_node==4:
        return ((0,0),), 2
    
    if element.n_node==8:
        if n==1:
            return ((0,0),), 2
        elif n==2:
            return ((1/math.sqrt(3),1/math.sqrt(3)),
                    (-1/math.sqrt(3),1/math.sqrt(3))
                    (-1/math.sqrt(3),-1/math.sqrt(3))
                    (1/math.sqrt(3),-1/math.sqrt(3))), 1
         
    

