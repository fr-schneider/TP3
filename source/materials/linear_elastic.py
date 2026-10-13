import numpy as np
  
class LinearElastic:
    def __init__(self, E, nu):
        self.E = E
        self.nu = nu


    def constitutive(self, tipo):
        if tipo == 'Plane_stress':
            factor = self.E / (1.0 - self.nu**2)
            return factor * np.array([[1.0,       self.nu,    0.0],
                                      [self.nu,   1.0,        0.0],
                                      [0.0,       0.0,        (1.0 - self.nu) / 2.0]])
        
        elif tipo == 'Plane_strain':
            factor = self.E / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
            return factor * np.array([[1.0 - self.nu, self.nu,            0.0],
                                      [self.nu,       1.0 - self.nu,      0.0],
                                      [0.0,           0.0,                (1.0 - 2.0 * self.nu) / 2.0]])