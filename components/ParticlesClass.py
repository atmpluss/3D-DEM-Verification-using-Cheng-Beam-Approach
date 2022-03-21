import numpy as np
import math

class Particles():

    def __init__(self, ID, radius,position,rho):
        self.id = ID
        self.r = radius
        self.rho = rho
        self.pos = position
        self.phi = np.array([0,0,0])
        self.v = np.array([0,0,0])
        self.an_v = np.array([0,0,0])
        self.a = np.array([0,0,0])
        self.an_a = np.array([0,0,0])
        self.F = np.array([0.,0.,0.])
        self.M = np.array([0.,0.,0.])
    def mass(self):  # particle mass
        return 4. / 3. * math.pi * self.rho * self.r ** 3

    def inertia(self):  # particle moment of inertia
        return 2. / 5. * self.mass() * self.r ** 2

