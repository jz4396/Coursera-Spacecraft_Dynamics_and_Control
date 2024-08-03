import numpy as np

def tilde(x):
    """
    Gets the skew symmetric matrix of vector x
    """
    x = np.array(x).reshape(3)
    x_tilde = np.array([[0, -x[2], x[1]], [x[2], 0, -x[0]], [-x[1], x[0], 0]])
    return x_tilde

def rotated(dim, deg):
    """
    Gets the DCM matrix for a rotation of deg degrees about the dim axis
    """
    rad = np.deg2rad(deg)
    return rotate(dim, rad)

def rotate(dim, rad):
    """
    Gets the DCM matrix for a rotation of rad radians about the dim axis
    """
    orders = np.array([[1, 2, 3], [3, 1, 2], [2, 3, 1]])
    r = np.array([[1, 0, 0], [0, np.cos(rad), np.sin(rad)], [0, -np.sin(rad), np.cos(rad)]])
    R = r[orders[dim-1]-1][:, orders[dim-1]-1]
    return R

def eulerA(dims, rads):
    """
    Gets the dcm matrix equivalent of euler angle set in radians
    """
    R = np.eye(3)
    for i in range(len(dims)):
        R = rotate(dims[i], rads[i]).dot(R)
    return R

def eulerAd(dims, deg):
    """
    Gets the dcm matrix equivalent of euler angle set in degrees
    """
    R = eulerA(dims, np.deg2rad(deg))
    return R

def dcm_2_quat(r):
    """
    Gets the quaternion rotation equivalent of dcm matrix r
    """
    b0 = np.sqrt(1/4*(1+np.trace(r)))
    b1 = np.sqrt(1/4*(1-np.trace(r)+2*r[0][0]))
    b2 = np.sqrt(1/4*(1-np.trace(r)+2*r[1][1]))
    b3 = np.sqrt(1/4*(1-np.trace(r)+2*r[2][2]))
    s1 = r[1][2]-r[2][1]
    b1 = s1/abs(s1)*b1
    s2 = r[2][0]-r[0][2]
    b2 = s2/abs(s2)*b2
    s3 = r[0][1]-r[1][0]
    b3 = s3/abs(s3)*b3
    ep = [b0, b1, b2, b3]
    return ep

def quat_2_dcm(q):
    """
    Gets the dcm matrix equivalent of quaternion rotation q
    """
    dcm = [q[0]**2+q[1]**2-q[2]**2-q[3]**2, 2*(q[1]*q[2]+q[0]*q[3]), 2*(q[1]*q[3]-q[0]*q[2]),
           2*(q[1]*q[2]-q[0]*q[3]), q[0]**2-q[1]**2+q[2]**2-q[3]**2, 2*(q[2]*q[3]+q[0]*q[1]),
           2*(q[1]*q[3]+q[0]*q[2]), 2*(q[2]*q[3]-q[0]*q[1]), q[0]**2-q[1]**2-q[2]**2+q[3]**2]
    return np.array(dcm).reshape(3, 3)

def quat_add(qs,q2):
    """
    Rotates attitude qs by q2 to achieve a final atitude qf
    """
    qs = np.array(qs).reshape(4, 1)
    m = np.array([[q2[0], -q2[1], -q2[2], -q2[3]],
                  [q2[1], q2[0], q2[3], -q2[2]],
                  [q2[2], -q2[3], q2[0], q2[1]],
                  [q2[3], q2[2], -q2[1], q2[0]],
                  ])
    qf = m@qs
    return qf

def quat_diff(qs,qf):
    """
    Gets the quaternion rotation that results in  qs -> qf
    """
    m = np.array([[qs[0], -qs[1], -qs[2], -qs[3]],
                  [qs[1], qs[0], -qs[3], qs[2]],
                  [qs[2], qs[3], qs[0], -qs[1]],
                  [qs[3], -qs[2], qs[1], qs[0]],
                  ])
    q2 = m.transpose()@qf
    return q2

def quat_revert(qf,q2):
    """
    Undoes the rotation of q2 from qf to result in qs
    """
    return quat_add(qf,q2*np.array([1, -1, -1, -1]))

def quat_omega_DKE(q, w):
    """
    Get quaternion rates from body omega rates and 
    and current attitude quaternion
    """
    m = np.array([[-q[1], -q[2], -q[3]],
                  [q[0], -q[3], q[2]],
                  [q[3], q[0], -q[1]],
                  [-q[2], q[1], q[0]],
                  ])
    q_dot = 0.5*m@w
    return q_dot

def quatdot_2_omega(q, qdot):
    m = np.array([[q[0], -q[1], -q[2], -q[3]],
                  [q[1], q[0], -q[3], q[2]],
                  [q[2], q[3], q[0], -q[1]],
                  [q[3], -q[2], q[1], q[0]],
                  ])
    omega = 2*m.transpose()@qdot
    return omega

def normalized_integrator(qdot,q0,t0,tf,dt):
    """
    Integrate quaternion rates to get new attitude quaternion
    """
    length=int((tf-t0)/dt)+1
    q=np.zeros([length,4])
    q[0,:]=np.array([q0])
    for i,t in enumerate(np.arange(t0, tf, dt)):
        q_new = q[i] + qdot(q[i],t)*dt
        q[i+1] = q_new/np.linalg.norm(q_new)
    return q

def quat_2_CRP(q):
    """
    Gets the Classical Rodriguez Parameters (CRP) equivalent of quaternion q
    """
    q = np.array(q).reshape(4, 1)
    p = q[1:]/q[0]
    return p

def CRP_2_quat(p):
    """
    Gets the quaternion equivalent of classical rodriguez parameters (CRP) p
    """
    p = np.array(p).reshape(3, 1)
    q = np.array([1, p[0], p[1], p[2]])
    q = q/np.sqrt(1+np.inner(p,p))
    return q


def CRP_2_DCM(p):
    """
    Gets the dcm matrix equivalent of crp p
    """
    p = np.array(p).reshape(3)
    scale = 1/(1+np.linalg.norm(p)**2)
    p_tilde = np.array([[0, -p[2], p[1]], [p[2], 0, -p[0]], [-p[1], p[0], 0]])
    dcm = scale*((1-np.linalg.norm(p)**2) * np.eye(3)
                  + 2*np.outer(p, p) 
                  - 2*p_tilde)
    return dcm

def CRP_add(pi,p2):
    """
    Adds two classical rodriguez parameters (CRP) as if inital attitude is pi 
    and we rotate by p2
    """
    pi=np.array(pi)
    p2=np.array(p2)
    pi = np.array(pi)
    p2 = np.array(p2)
    p = (pi + p2 - np.cross(p2, pi))/(1-p2.dot(pi))
    return p

def CRP_diff(pi,pf):
    """
    Gets the Classical Rodriguez Parameters (CRP) would rotate between inital
    attitude pi and final attitude pf
    """
    pi=np.array(pi).reshape(3,1)
    pf=np.array(pf).reshape(3,1)
    p = (pf - pi + np.cross(pf, pi))/(1+pf.dot(pi))
    return p

def CRP_omega_DKE(p, w):
    """
    Get CRP rates from body omega rates and current attitude CRP
    """
    w = np.array(w)
    m=0.5*np.array([[1+p[0]**2, p[0]*p[1]-p[2], p[0]*p[2]+p[1]],
                [p[1]*p[0]+p[2], 1+p[1]**2, p[1]*p[2]-p[0]],
                [p[2]*p[0]-p[1], p[2]*p[1]+p[0], 1+p[2]**2]])
    p_dot = m@w
    return p_dot

def CRPdot_2_omega(p, pdot):
    p=np.array(p).reshape(3,1)
    pdot = np.array(pdot).reshape(3,1)
    p_tilde = np.array([[0, -p[2], p[1]], [p[2], 0, -p[0]], [-p[1], p[0], 0]])
    m=2./(1+np.inner(p,p))*(np.eye(3)-p_tilde)
    omega = m@pdot
    return omega

def quat_2_MRP(q):
    """
    Gets the Modified Rodriguez Parameters (MRP) equivalent of quaternion q
    """
    q = np.array(q)
    p = q[1:]/(1+q[0])
    return p

def MRP_2_quat(p):
    """
    Gets the quaternion equivalent of modified rodriguez parameters (MRP) p
    """
    p = np.array(p)
    q = np.array([1-np.inner(p,p), 2*p[0], 2*p[1], 2*p[2]])
    q = q/(1+np.inner(p,p))
    return q

def MRP_shadow_set(p):
    """
    Gets the shadow set of the modified rodriguez parameters (MRP) p
    """
    p = np.array(p)
    p = -p/(np.inner(p,p))
    return p

def MRP_2_DCM(p):
    """
    Gets the dcm matrix equivalent of mrp p
    """
    q = MRP_2_quat(p)
    dcm = quat_2_dcm(q)
    return dcm

def DCM_2_MRP(dcm):
    """
    Gets the Modified Rodriguez Parameters (MRP) equivalent of dcm matrix dcm
    """
    q = dcm_2_quat(dcm)
    p = quat_2_MRP(q)
    return p

def MRP_add(pi,p2):
    """
    Adds two modified rodriguez parameters (MRP) as if inital attitude is pi 
    and we rotate by p2
    """
    pi=np.array(pi)
    p2=np.array(p2)
    pi_squared = np.inner(pi,pi)
    p2_squared = np.inner(p2,p2)
    denom=1+pi_squared*p2_squared-2*np.dot(pi,p2)
    if denom<=1e-2:
        #print("reshaping")
        if pi_squared<p2_squared:
            p2 = MRP_shadow_set(p2)
            p2_squared = np.inner(p2,p2)
        else:
            pi = MRP_shadow_set(pi)
            pi_squared = np.inner(pi,pi)
        denom=1+pi_squared*p2_squared-2*np.dot(pi,p2)
    pf=((1-pi_squared)*p2+(1-p2_squared)*pi-2*np.cross(p2,pi))/denom
    return pf

def MRP_diff(pi,pf):
    """
    Finds the modified rodriguex parameters (MRP) that would rotate between
    inital attitude pi and final attitude pf
    """
    pi=np.array(pi)
    pf=np.array(pf)
    return MRP_add(-pi,pf)

def MRP_short(p):
    """
    Gets the shortest rotation equivalent of modified rodriguez parameters (MRP) p
    """
    p = np.array(p).reshape(3)
    if np.inner(p,p)>1:
        p = MRP_shadow_set(p)
    return p

def MRP_omega_DKE(p, w):
    """
    Get MRP rates from current attitude MRP p and body omega rates w
    """
    p = np.array(p)
    w = np.array(w).reshape(3)
    k = 1-np.inner(p,p)
    M=1/4* np.array([[k+2*p[0]**2, 2*(p[0]*p[1]-p[2]), 2*(p[0]*p[2]+p[1])],
                    [2*(p[1]*p[0]+p[2]), k+2*p[1]**2, 2*(p[1]*p[2]-p[0])],
                    [2*(p[2]*p[0]-p[1]), 2*(p[2]*p[1]+p[0]), k+2*p[2]**2]])
    p_dot = M@w
    return p_dot

def MRPdot_2_omega(p, p_dot):
    """
    Get MRP rates from current attitude MRP p and body omega rates w
    """
    p = np.array(p)
    p_dot = np.array(p_dot).reshape(3)
    k = 1-np.inner(p,p)
    B = np.array([[k+2*p[0]**2, 2*(p[0]*p[1]-p[2]), 2*(p[0]*p[2]+p[1])],
                    [2*(p[1]*p[0]+p[2]), k+2*p[1]**2, 2*(p[1]*p[2]-p[0])],
                    [2*(p[2]*p[0]-p[1]), 2*(p[2]*p[1]+p[0]), k+2*p[2]**2]])
    # w = 4*np.linalg.inv(B)@p_dot
    w = 4/(1+np.inner(p,p))**2*B.transpose()@p_dot
    return w

def MRP_short_integrator(pdot,p0,t0,tf,dt):
    """
    Integrate quaternion rates to get new attitude quaternion
    """
    length=int((tf-t0)/dt)+1
    p=np.zeros([length,3])
    p[0,:]=np.array([p0])
    for i,t in enumerate(np.arange(t0, tf, dt)):
        p_new = p[i] + pdot(p[i],t).transpose()*dt
        p[i+1,:] = MRP_short(p_new)
    return p

def triad_method_estimator(b1, b2, r1, r2):
    """
    Gets the dcm matrix equivalent of the triad method
    """
    b1 = np.array(b1).reshape(3)
    b2 = np.array(b2).reshape(3)
    r1 = np.array(r1).reshape(3)
    r2 = np.array(r2).reshape(3)
    b1 = b1/np.linalg.norm(b1)
    b2 = b2/np.linalg.norm(b2)
    r1 = r1/np.linalg.norm(r1)
    r2 = r2/np.linalg.norm(r2)
    b2 = np.cross(b1, b2)
    r2 = np.cross(r1, r2)
    b2 = b2/np.linalg.norm(b2)
    r2 = r2/np.linalg.norm(r2)
    b3 = np.cross(b1, b2)
    r3 = np.cross(r1, r2)
    BR = np.array([b1, b2, b3]).transpose()
    NR = np.array([r1, r2, r3]).transpose()
    R = BR@NR.transpose()
    return R

def get_q_matricies(b_vs, r_vs, ws):
    b_vs = np.array(b_vs).reshape(-1,3)
    r_vs = np.array(r_vs).reshape(-1,3)
    b_vs = b_vs/np.linalg.norm(b_vs, axis=1).reshape(-1,1)
    #r_vs = r_vs/np.linalg.norm(r_vs, axis=1).reshape(-1,1)
    ws = np.array(ws).reshape(-1,1)

    B = np.zeros([3,3])
    for i in range(len(ws)):
        B = B + ws[i]*np.outer(b_vs[i], r_vs[i])
    
    #print(B)
    S = B + B.transpose()
    #print(S)
    Z = np.array([[B[1][2]-B[2][1]], [B[2][0]-B[0][2]], [B[0][1]-B[1][0]]])
    #print(Z)
    sigma = np.array([np.trace(B)]).reshape(1,1)
    K = np.vstack((np.hstack((sigma, Z.transpose())), np.hstack((Z, S-sigma*np.eye(3)))))
    #print(K)
    return K, B, S, Z, sigma

def q_method_estimator(b_vs, r_vs, ws):
    """
    Gets the quaternion vector using the q method and weighted measurments
    """
    K,_,_,_,_=get_q_matricies(b_vs, r_vs, ws)
    w, v = np.linalg.eig(K)
    #print(w,"\n",v)
    #print(np.argmax(w))
    q = v[:, np.argmax(w)]
    return q

def quest_method_estimator(b_vs, r_vs, ws, iterations):
    """
    Gets the CRP vector using the quest method and weighted measurments
    """
    K, _, S, Z, sigma = get_q_matricies(b_vs, r_vs, ws)

    lambda_opt = np.sum(ws)
    for _ in range(iterations):
        A = K-lambda_opt*np.eye(4)
        delta = 1/np.trace(np.linalg.inv(A)*(-1*np.eye(4)))
        lambda_opt = lambda_opt - delta
    q = np.linalg.inv((lambda_opt+sigma)*np.eye(3)-S)@np.array(Z)
    return q

def OLAE_estimator(b_vs, r_vs, ws):
    """
    Gets the quaternion vector using the OLAE method and weighted measurments
    """
    b_vs = np.array(b_vs).reshape(-1,3)
    r_vs = np.array(r_vs).reshape(-1,3)
    ws = np.array(ws).flatten()
    size = len(ws)
    b_vs = b_vs/np.linalg.norm(b_vs, axis=1).reshape(-1,1)
    r_vs = r_vs/np.linalg.norm(r_vs, axis=1).reshape(-1,1)
    D = (b_vs - r_vs).flatten()
    s = b_vs + r_vs
    W = np.zeros([size*3,size*3])
    S = np.zeros([size*3,3])
    for i in range(len(s)):
        S[i*3:i*3+3] = tilde(s[i])
        W[i*3:i*3+3,i*3:i*3+3] = np.eye(3)*ws[i]


    #print(f"Big S:\n{S.shape}")
    #print(f"D:\n{D}")
    #print(f"W:\n{W.shape}")
    #print((np.linalg.inv(S.transpose()@W@S)@S.transpose()@W).shape)
    q = np.linalg.inv(S.transpose()@W@S)@S.transpose()@W@D
    return q
