import numpy as np
from numpy.typing import ArrayLike
from typing import Optional
import sys

sys.path.append("../")

import attitude_math as am

R_mars = 3396.19e3
grav_const = 42828.3

#@dataclass
class Circular_Orbit:
    def __init__(self,
                 h: float,
                 EA_i: ArrayLike,
                 simga_i: Optional[ArrayLike] = None,
                 omega_i: Optional[ArrayLike] = None,
                 MOI: Optional[ArrayLike] = None):
        self.h = h
        self.EA_i = EA_i
        self.sigma_i = simga_i  # Initial attitude
        self.omega_Bi = omega_i # Initial body rotation in body frame
        self.MOI = MOI
        self.r = R_mars + h
        self.theta_dot = np.sqrt(grav_const*1e9/(self.r**3))
        # omega of hill frame relative to inertial in H frame
        self.omega_H_h = [0, 0, self.theta_dot]

    orbit_euler_func = lambda self, t: self.EA_i + [0,0, self.theta_dot*t]
    R_HN = lambda self, t: am.eulerA([3,1,3], self.orbit_euler_func(t))
    R_RnN = lambda self, t: R_RnH@self.R_HN(t)
    omega_RnN_N = lambda self, t: self.R_HN(t).transpose() @ self.omega_H_h

    def inertial_state(self, t):
        # Returns inertial position in the N frame and
        # velocity in the Hill frame
        R = self.R_HN(t).transpose()
        #print(euler_a,R)
        body_pos = np.array([self.r, 0, 0])
        inertial_pos = R @ body_pos

        body_vel = np.array([0, self.theta_dot*self.r, 0])
        inertial_v = R @ body_vel
        return np.array([inertial_pos,inertial_v])

   
LMO_sat = Circular_Orbit(
    h = 400e3,
    EA_i = np.deg2rad([20, 30, 60]),
    simga_i = np.array([0.3, -0.4, 0.5]),
    omega_i = np.array([1, 1.75, -2.2]),
    MOI =  np.diag([10, 5, 7.5]),
)
GMO_sat = Circular_Orbit(h = 20424.2e3-R_mars,
                                       EA_i = np.deg2rad([0,0,250]),
)

# Nadir can be define in reference to Hill frame
R_RnH = am.eulerAd([2],[180])
R_RsN = am.eulerAd([1,2],[90,180])

omega_RsN = np.zeros(3)