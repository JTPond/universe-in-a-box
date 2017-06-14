import numpy as np

class particle (object):

    def __init__(self,mass,radius,charge,position,velocity,scale):
        self.mass = mass
        self.radius = radius
        self.charge = charge
        self.position = position
        self.velocity = velocity
        self.scale = scale

    def set_m(self,val):
        self.mass = val

    def set_rho(self,val):
        self.radius = val

    def set_q(self,val):
        self.charge = val

    def set_r(self,val):
        self.position = val

    def set_v(self,val):
        self.velocity = val
 

    def m(self):
        return self.mass

    def rho(self):
        return self.radius

    def q(self):
        return self.charge

    def r(self):
        return self.position

    def v(self):
        return self.velocity

    def force(self,particle2):
        G = 6.67408e-11 * 100.0**3.0
        k = 8.987e9 * 100.0**3.0
        rmrp = self.position - particle2.r()
        R2 = (rmrp*rmrp).sum()
        
        return (-(G*self.mass*particle2.m())/R2 + (k*self.charge*particle2.q())/R2)*(rmrp/np.sqrt(R2))
