import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from typing import Optional
import sys

sys.path.append("../")
sys.path.append("C:/Users/jaz43/Documents/2024/Coursera Spacecraft Dynamics/")

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
                 MOI: Optional[ArrayLike] = np.identity(3)):
        self.h = h
        self.EA_i = EA_i
        self.sigma_i = simga_i  # Initial attitude
        self.omega_Bi = omega_i # Initial body rotation in body frame
        self.sigma_BN = simga_i
        self.omega_BN_B = omega_i
        self.MOI = MOI
        self.MOI_inv = np.linalg.inv(MOI)
        self.r = R_mars + h
        self.theta_dot = np.sqrt(grav_const*1e9/(self.r**3))
        # omega of hill frame relative to inertial in H frame
        self.omega_H_h = [0, 0, self.theta_dot]

    X = lambda self: np.vstack([self.sigma_BN,self.omega_BN_B])
    orbit_euler_func = lambda self, t: self.EA_i + [0,0, self.theta_dot*t]
    R_HN = lambda self, t: am.eulerA([3,1,3], self.orbit_euler_func(t))
    R_RnN = lambda self, t: R_RnH@self.R_HN(t)
    omega_RnN_N = lambda self, t: self.R_HN(t).transpose() @ self.omega_H_h
    
    def reset_X(self):
        self.sigma_BN = self.sigma_i
        self.omega_BN_B = self.omega_Bi

    def set_X(self, X):
        self.sigma_BN = X[0]
        self.omega_BN_B = X[1]

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
    
    def reference_error_state(self, DCM_RrN, omega_RrN_N):
        DCM_BN = am.MRP_2_DCM(self.sigma_BN)
        DCM_BR = DCM_BN@DCM_RrN.transpose()
        MRP_BR = am.DCM_2_MRP(DCM_BR)
        omega_BR_B = self.omega_BN_B - DCM_BN@omega_RrN_N
        return np.vstack([MRP_BR,omega_BR_B])
    
    def X_dot(self, u, X=None):
        if X is None:
            X = self.X()
        sigma_BN = X[0]
        omega_BN_B = X[1]

        sigma_dot = am.MRP_omega_DKE(sigma_BN,omega_BN_B)
        omega_dot = self.MOI_inv @ (
            u - np.cross(omega_BN_B, self.MOI@omega_BN_B)
            ) 
        return np.vstack([sigma_dot,omega_dot])

def orbital_att_RK4int(sat_obj, u, L, X0, ts):
    columns = ["t",
               "sigma1","sigma2","sigma3",
               "omega1","omega2","omega3"]
    #df = pd.DataFrame(np.hstack([ts[0],X0[0],X0[1]]).reshape(1,7),columns=columns)
    df = pd.DataFrame(np.hstack([ts.reshape(-1,1), np.zeros((len(ts),6))]),
                      columns=columns)
    df.loc[0,columns[1:]]=X0.reshape(6)
    t_old = ts[0]
    for i, t in enumerate(ts[1:]):
        X = sat_obj.X()
        dt = (t-t_old)
        control = u(X,t)

        K1 = sat_obj.X_dot(control,X)
        K2 = sat_obj.X_dot(control,X+(dt*K1/2))
        K3 = sat_obj.X_dot(control,X+(dt*K2/2))
        K4 = sat_obj.X_dot(control,X+(dt*K3))
        X_dot = 1/8*(K1+3*K2+3*K3+K4)
        #X_dot = sat_obj.X_dot(u(X,t),X)
        X_new = X + X_dot*dt

        t_old=t
        X_new[0] = am.MRP_short(X_new[0])
        
        df.iloc[i+1]=(np.hstack([t,X_new.flatten()]))
        sat_obj.set_X(X_new)
    return df
   
LMO_sat = Circular_Orbit(
    h = 400e3,
    EA_i = np.deg2rad([20, 30, 60]),
    simga_i = np.array([0.3, -0.4, 0.5]),
    omega_i = np.deg2rad([1, 1.75, -2.2]),
    MOI =  np.diag([10, 5, 7.5]),
)
GMO_sat = Circular_Orbit(h = 20424.2e3-R_mars,
                                       EA_i = np.deg2rad([0,0,250]),
)

# Nadir can be define in reference to Hill frame
R_RnH = am.eulerAd([2],[180])
# Sun pointing reference frame
R_RsN = am.eulerAd([1,2],[90,180])

omega_RsN = np.zeros(3)

# Create a function to build a DCM from component vectors
def Get_Comm_Frame_DCM(t):
    delta_r_N = GMO_sat.inertial_state(t)[0] - LMO_sat.inertial_state(t)[0]
    r1 = -delta_r_N/np.linalg.norm(delta_r_N)
    r2 = np.cross(delta_r_N, [0,0,1])
    r2 = r2/np.linalg.norm(r2)
    r3 = np.cross(r1,r2)
    R_RcN = np.vstack([r1,r2,r3])
    return R_RcN

def Get_Comm_Frame_Inertial_Rate(t1,dt):
    t0 = t1-dt
    dcm1 = Get_Comm_Frame_DCM(t1)
    dcm0 = Get_Comm_Frame_DCM(t0)
    omega = am.dcmdot_2_oemga(dcm1,dcm0,t1-t0)
    omega_RcN_N = dcm1.transpose()@omega
    return omega_RcN_N


if __name__ == "__main__":

    ## Code to compare integrator against controls Module 3 Concept Check 2
    Test_sat = Circular_Orbit(h=400e3,
                              EA_i = np.deg2rad([0,0,0]),
                              simga_i=np.array([0.1,.2,-.1]),
                              omega_i=np.deg2rad([30,10,-20]),
                              MOI=np.diag([100,75,80]))

    K = 5
    P = 10 * np.eye(3)
    freq = 0.05
    t_step = .1
    # Reference Tracking State Functions
    sigma_r = lambda t, f :  [0.2*np.sin(f*t), 0.3*np.cos(f*t), -0.3*np.sin(f*t)]
    sigma_dot_r = lambda t, f : [0.2*f*np.cos(f*t), -0.3*f*np.sin(f*t), -0.3*f*np.cos(f*t)]
    def state_r_func(t, dt, f) : 
        sigma = sigma_r(t, f)
        sigma_dot = sigma_dot_r(t, f)
        sigma_old = sigma_r(t-dt, f)
        sigma_dot_old = sigma_dot_r(t-dt, f)

        # 
        omega_r = am.MRPdot_2_omega(sigma, sigma_dot)
        omega_r_old = am.MRPdot_2_omega(sigma_old, sigma_dot_old)

        omega_r_dot = (omega_r-omega_r_old)/dt
        return np.concatenate([sigma, omega_r, omega_r_dot])

    def u(X,t):
        state_r = state_r_func(t, 0.01, freq)

        sigma = X[0]
        omega = X[1]
        
        att_diff = am.DCM_2_MRP(am.MRP_2_DCM(sigma)@am.MRP_2_DCM(state_r[:3]).transpose())
        omega_r_br = am.MRP_2_DCM(att_diff)@state_r[3:6]
        omega_rdot_br = am.MRP_2_DCM(att_diff)@state_r[6:]
        control = - K*att_diff - P@(omega-omega_r_br)
        return control

    state = orbital_att_RK4int(Test_sat,
                            u,
                            0,
                            Test_sat.X(),
                            np.arange(0., 20+t_step, t_step))
    time_measure=20
    answer_i  = int(time_measure/t_step)
    state_val = state.loc[answer_i,["sigma1","sigma2","sigma3"]].to_numpy()
    del_sigma = am.DCM_2_MRP(am.MRP_2_DCM(state_val)@am.MRP_2_DCM(sigma_r(time_measure, freq)).transpose())
    print(np.linalg.norm(del_sigma), del_sigma)